//! Per-read alignment driver: seed → anchor → chain → (L1 / intron-chain /
//! spliced / EA-Myers / banded DP) → split/span rescue → recount → coverage gate.
//!
//! All production gates are frozen (no env switches): the splice gate (skip try_spliced on
//! <2 segments), the non-best-chain spliced Myers pre-gate (skip the spliced DP when the
//! pseudo-reference is farther than read_len/3 edits from the read), the tail-rescue gate
//! (small-tail rescue tightening), multitie marking, EA-Myers decisive path, window cap.
//! Profiling instrumentation is not ported.

use std::collections::HashSet;
use std::sync::Arc;

use crate::chain::{
    chain_anchors_into, second_score_with, tied_best, Chain, ChainScratch, TOP_CANDIDATES,
};
use crate::extend::{
    extend_hint, CigarOp, DiagHint, EditDir, ExtendBuffer, ExtendParams, Extension,
};
use crate::gtf::{Junction, JunctionLib, RefinedJunction, SpliceSignal};
use crate::index::{collect_anchors_into, Index};
use crate::mapq::ReadAlignment;
use crate::seed::{minimizers_into, Minimizer, SeedParams, SeedScratch, Strand};
use crate::splice::{
    align_spliced, intron_length_penalty, refine_junction_read, segment_chain,
    splice_signal, SpliceParams,
};

/// Aligner configuration (`default` = DNA preset; `rna_default` = RNA production).
#[derive(Debug, Clone, Copy)]
pub struct AlignConfig {
    pub seed: SeedParams,
    pub chain: crate::chain::ChainParams,
    pub extend: ExtendParams,
    /// RNA mode: dense seeding + intron-aware chaining + splice machinery.
    pub rna: bool,
}

impl Default for AlignConfig {
    /// DNA preset (k=15, w=10).
    fn default() -> Self {
        AlignConfig {
            seed: SeedParams { k: 15, w: 10 },
            chain: crate::chain::ChainParams::default(),
            extend: ExtendParams::default(),
            rna: false,
        }
    }
}

impl AlignConfig {
    /// RNA mode defaults (spec §parameter table).
    pub fn rna_default() -> Self {
        AlignConfig {
            seed: SeedParams { k: 15, w: 5 },
            chain: crate::chain::ChainParams::rna_default(),
            extend: ExtendParams::default(),
            rna: true,
        }
    }
}

/// Second-round (mid-occ) seed occurrence cap (minimap2 --mid-occ semantics).
pub const MID_OCC_CAP: u32 = 5000;
/// Max candidates entering base-level tie-breaking (deterministic truncation).
pub const MAX_TIED_CANDIDATES: usize = 8;
/// Extension window hard cap (garbage chains must not blow up the DP).
pub const MAX_EXTEND_WINDOW: u32 = 20_000;
/// Extension tie-break budget (top-N candidates extended by banded DP).
pub const EXTEND_BUDGET: usize = 2;

/// Outcome of one seeding+chaining round.
pub(crate) struct RoundOutcome {
    chains: Vec<Chain>,
    second: i32,
}

impl RoundOutcome {
    fn best(&self) -> Option<&Chain> {
        self.chains.first()
    }
}

/// Minimum fraction of a read's bases a placement must explain, below which
/// Tier-1 refuses it outright (STAR's `outFilterMatchNminOverLread`
/// community default). Refused reads are unmapped and feed the Tier-2
/// collapsed-alphabet rescue pool.
pub const REFUSE_MIN_ALIGNED_FRAC: f64 = 0.66;

/// Read-span coverage of the round-1 best chain below which the
/// mid-occurrence seeding round runs (see `genomic_candidates`).
pub const MID_ROUND_MIN_COVERAGE: f64 = 0.6;

/// Anchor budget for the mid-occurrence round; beyond it the round is
/// skipped rather than chained (chaining cost grows with the square of the
/// anchors sharing a contig and strand).
pub const MID_ROUND_MAX_ANCHORS: usize = 20_000;

/// Half-width of the mate-anchored window, on top of the insert estimate
/// (`est_insert + est_insert/2 + read_len` each side of the mate).
pub const MATE_WINDOW_INSERT: u32 = 500;

/// Strict-domain score of one placement candidate (see
/// [`Aligner::rescore_strict`]). All candidates of a read are scored by the
/// same function with the same parameters, so these numbers are comparable
/// across layers and decide both placement and MAPQ.
#[derive(Clone, Debug, Default)]
pub(crate) struct StrictScore {
    /// Total strict score (match/mismatch/indel/splice-signal sums).
    pub score: i32,
    /// Read bases in M/I columns.
    pub aligned: u32,
    /// Exact matches in M columns.
    pub matches: u32,
    /// Ordinary mismatches (N on either side excluded).
    pub mm: u32,
    /// A-to-G-direction substitutions (ref A / read G), charged as mismatches.
    pub ea_ag: u32,
    /// T-to-C-direction substitutions (ref T / read C), charged as mismatches.
    pub ea_tc: u32,
    /// Insertion events.
    pub ins: u32,
    /// Deletion events.
    pub del: u32,
    /// Soft-clipped bases.
    pub clipped: u32,
    /// Intron (RefSkip) count.
    pub n_splice: u32,
}

/// One placement candidate from either layer, pending strict arbitration.
#[derive(Clone)]
pub(crate) struct Placement {
    pub origin: crate::mapq::Origin,
    pub contig: u32,
    pub pos: u32,
    pub strand: Strand,
    pub cigar: Vec<CigarOp>,
    pub junctions: Vec<RefinedJunction>,
    /// Chain score of the source chain (genomic layer; 0 for L1).
    pub chain_score: i32,
    /// Score the producing extension reported (domain: extension scores).
    pub ext_score: i32,
    pub strict: StrictScore,
    /// Editing direction the candidate was last scored under (`Off` = strict).
    pub ea_dir: EditDir,
}

impl Placement {
    fn from_genomic(
        contig: u32,
        pos: u32,
        strand: Strand,
        cigar: Vec<CigarOp>,
        junctions: Vec<RefinedJunction>,
        chain_score: i32,
        ext_score: i32,
    ) -> Placement {
        Placement {
            origin: crate::mapq::Origin::Genomic,
            contig,
            pos,
            strand,
            cigar,
            junctions,
            chain_score,
            ext_score,
            strict: StrictScore::default(),
            ea_dir: EditDir::Off,
        }
    }
}

/// One fully-scored placement from the exhaustive path (verify stage): the
/// strict score plus the reference M blocks, everything the vote-withdrawal
/// checks need.
#[derive(Clone, Debug)]
pub struct ScoredPlacement {
    /// Contig index into the alignment index.
    pub contig: u32,
    /// Leftmost 0-based reference position.
    pub pos: u32,
    /// Read strand.
    pub strand: Strand,
    /// Strict-domain score.
    pub score: i32,
    /// Reference M blocks `(start, end)`, half-open.
    pub blocks: Vec<(u32, u32)>,
    /// Full CIGAR (reference orientation) for column mapping.
    pub cigar: Vec<CigarOp>,
}

/** Genomic layer output: per-chain placement candidates plus read evidence.
 */
pub(crate) struct GenomicBatch {
    pub placements: Vec<Placement>,
    pub n_seeds: usize,
}

/** Per-thread aligner state (reuses DP buffers).
 */
pub struct Aligner<'a> {
    pub index: &'a Index,
    pub config: AlignConfig,
    /// Junction library (sjdb / 2-pass discoveries; RNA mode).
    pub jlib: Option<Arc<JunctionLib>>,
    /// L1 transcriptome-first engine bundle (None = pure legacy path).
    pub l1: Option<Arc<esperanto_engine::L1Index>>,
    pub(crate) buf: ExtendBuffer,
    /// Reusable hot-path buffers (seeding / chaining / revcomp / windows).
    /// Allocation-elimination only: every buffer is cleared before refill,
    /// so results are byte-identical to the fresh-allocation path.
    pub(crate) scratch: MapScratch,
}

/** Reusable buffers for the seeding/chaining hot path. Owned by the
 * [`Aligner`] (one per rayon worker, via `map_init`) so per-call `Vec`
 * allocations become `clear()` + refill. No algorithm state lives here.
 */
#[derive(Default)]
pub(crate) struct MapScratch {
    seed: SeedScratch,
    chain: ChainScratch,
    second_keys: HashSet<(u32, u32)>,
    anchors: Vec<crate::seed::Anchor>,
    mins_out: Vec<Minimizer>,
    chains_out: Vec<Chain>,
    rc: Vec<u8>,
    win: Vec<u8>,
    myers_a: Vec<u8>,
    myers_b: Vec<u8>,
    rescore: Vec<u8>,
    kmer_scan: Vec<u8>,
}

/** The scored candidate list of one `align_read` call, kept for the
 * mate-anchored re-arbitration so the read is not aligned twice.
 */
#[derive(Clone)]
pub struct CachedPlacements {
    placements: Vec<Placement>,
    n_seeds: usize,
    genomic_any: bool,
}

impl<'a> Aligner<'a> {
    pub fn new(index: &'a Index, config: AlignConfig) -> Self {
        Aligner {
            index,
            config,
            jlib: None,
            l1: None,
            buf: ExtendBuffer::default(),
            scratch: MapScratch::default(),
        }
    }

    pub fn with_lib(mut self, lib: Arc<JunctionLib>) -> Self {
        self.jlib = Some(lib);
        self
    }

    /// Enable the L1 transcriptome-first fast path: hit → direct placement;
    /// miss (Fallback) → transparent drop to the legacy genomic path.
    pub fn with_l1(mut self, l1: Arc<esperanto_engine::L1Index>) -> Self {
        self.l1 = Some(l1);
        self
    }

    /// Take the reusable reverse-complement buffer, filled with
    /// `revcomp(seq)` (allocation-elimination only; contents identical).
    /// Give it back with `self.scratch.rc = v` when done.
    fn take_rc(&mut self, seq: &[u8]) -> Vec<u8> {
        let mut v = std::mem::take(&mut self.scratch.rc);
        revcomp_into(seq, &mut v);
        v
    }

    /// Take the reusable reference-window buffer, decoded over
    /// `[wstart, wend)` (identical bytes to `contig.slice_ascii`).
    fn take_win(&mut self, wstart: u32, wend: u32, contig: &crate::fasta::Contig) -> Vec<u8> {
        let mut v = std::mem::take(&mut self.scratch.win);
        v.clear();
        contig.decode_into(wstart, wend, &mut v);
        v
    }

    /// Convert one L1 engine candidate into a genomic-layer placement.
    ///
    /// Conservation guard: the projected CIGAR's query consumption (M+I+S)
    /// must equal read_len, otherwise the projection is discarded.
    fn l1_placement(&self, read_len: usize, c: esperanto_engine::L1Candidate) -> Option<Placement> {
        let genome_contig = self
            .l1
            .as_ref()
            .and_then(|l1| l1.contig_name(c.contig))
            .and_then(|name| self.index.reference.contig_index(name.as_bytes()))?;

        let strand = match c.strand {
            esperanto_engine::Strand::Plus => Strand::Plus,
            esperanto_engine::Strand::Minus => Strand::Minus,
        };
        let cigar: Vec<CigarOp> = c
            .cigar
            .iter()
            .map(|op| match op {
                esperanto_engine::CigarOp::Match(n) => CigarOp::Match(*n),
                esperanto_engine::CigarOp::Ins(n) => CigarOp::Ins(*n),
                esperanto_engine::CigarOp::Del(n) => CigarOp::Del(*n),
                esperanto_engine::CigarOp::RefSkip(n) => CigarOp::RefSkip(*n),
                esperanto_engine::CigarOp::SoftClip(n) => CigarOp::SoftClip(*n),
            })
            .collect();
        let qspan: usize = cigar
            .iter()
            .map(|op| match op {
                CigarOp::Match(n) | CigarOp::Ins(n) | CigarOp::SoftClip(n) => *n as usize,
                _ => 0,
            })
            .sum();
        if qspan != read_len {
            return None; // non-conserving projection
        }
        Some(Placement {
            origin: crate::mapq::Origin::L1,
            contig: genome_contig,
            pos: c.pos,
            strand,
            cigar,
            junctions: Vec::new(),
            chain_score: 0,
            ext_score: 0,
            strict: StrictScore::default(),
            ea_dir: EditDir::Off,
        })
    }

    /// Align one read (ASCII ACGTN).
    ///
    /// Tier-1 unified arbitration: both layers act as candidate generators —
    /// the L1 engine contributes every defensible transcriptome projection,
    /// the genomic layer one extension per chained locus — and every
    /// candidate is rescored by one strict four-letter function with the
    /// editing exemption off. The best strict score decides the placement
    /// (ties prefer the genomic layer, whose CIGAR was computed against the
    /// genome directly), the margin to the best non-overlapping runner-up
    /// decides MAPQ, and a placement explaining less than
    /// [`REFUSE_MIN_ALIGNED_FRAC`] of the read is refused outright so the
    /// read reaches the Tier-2 pool.
    pub fn align_read(&mut self, seq: &[u8]) -> Option<ReadAlignment> {
        let (aln, _) = self.align_read_full(seq);
        aln
    }

    /// `align_read` plus the scored candidate list, for the mate-anchored
    /// re-arbitration to reuse instead of re-aligning the read from scratch.
    pub fn align_read_full(
        &mut self,
        seq: &[u8],
    ) -> (Option<ReadAlignment>, Option<CachedPlacements>) {
        let Some((placements, n_seeds, genomic_any)) = self.scored_placements(seq, false)
        else {
            return (None, None);
        };
        let cache = CachedPlacements {
            placements: placements.clone(),
            n_seeds,
            genomic_any,
        };
        (self.arbitrate(seq, placements, n_seeds, genomic_any), Some(cache))
    }

    /// Collect, strictly rescore, dedup and order the candidate placements of
    /// a read from both layers. `exhaustive` widens the genomic extension to
    /// every chained locus (TOP_CANDIDATES) and disables the decisive fast
    /// lane — the verify stage's re-alignment mode; the map stage uses the
    /// budgeted mode. Returns None when no layer produced any candidate.
    fn scored_placements(
        &mut self,
        seq: &[u8],
        exhaustive: bool,
    ) -> Option<(Vec<Placement>, usize, bool)> {
        let mut placements: Vec<Placement> = Vec::new();

        // 1. Candidate generation — transcriptome layer.
        if let Some(l1) = self.l1.clone() {
            let mut stats = esperanto_engine::ReadStats::default();
            for c in esperanto_engine::align_read_candidates(
                seq,
                l1.as_ref(),
                l1.as_ref(),
                l1.as_ref(),
                &esperanto_engine::EngineConfig::default(),
                &mut stats,
            ) {
                if let Some(p) = self.l1_placement(seq.len(), c) {
                    placements.push(p);
                }
            }
        }

        // 2. Candidate generation — genomic layer (always consulted; the L1
        // full-branch early return is abolished).
        let batch = self.genomic_candidates(seq, exhaustive);
        let n_seeds = batch.n_seeds;
        let genomic_any = !batch.placements.is_empty();
        placements.extend(batch.placements);
        if placements.is_empty() {
            return None;
        }

        // 3. Unified strict rescoring: same function, same parameters, no
        // editing exemption.
        let rc_storage = self.take_rc(seq);
        for p in &mut placements {
            let q: &[u8] = if p.strand == Strand::Minus {
                &rc_storage
            } else {
                seq
            };
            self.rescore_strict(q, p);
        }

        // 4. Same-locus dedup, then the decision order.
        dedup_by_footprint(&mut placements);
        sort_placements(&mut placements);

        // 5. Hyperedited read: re-arbitrate in its own reduced alphabet (see
        // `rescore_with`).
        self.hyper_rearbitrate(seq, &rc_storage, &mut placements);
        self.scratch.rc = rc_storage;
        Some((placements, n_seeds, genomic_any))
    }

    /// Hyperedited-read re-arbitration (see `rescore_with`). The signature is
    /// looked for in every candidate, not only the strict winner: the strict
    /// extension soft-clips an edited tail (each edit costs 4), so the winner
    /// is often the clipped version whose remaining edit density no longer
    /// qualifies, while a full-length candidate of the same locus (L1
    /// projection, or the clipped one re-extended) shows the signature
    /// plainly. When found, every unspliced genomic candidate is re-extended
    /// and every candidate rescored under that direction, then dedup and
    /// order are recomputed. No-op otherwise.
    fn hyper_rearbitrate(&mut self, seq: &[u8], rc_storage: &[u8], placements: &mut Vec<Placement>) {
        let Some(dir) = placements.iter().find_map(hyper_direction) else {
            return;
        };
        for p in placements.iter_mut() {
            let q: &[u8] = if p.strand == Strand::Minus {
                rc_storage
            } else {
                seq
            };
            if p.origin == crate::mapq::Origin::Genomic {
                self.reextend_with(q, p, dir);
            }
            self.rescore_with(q, p, dir);
        }
        dedup_by_footprint(placements);
        sort_placements(placements);
    }

    /// Align a read whose independent placement is poor, adding the locus
    /// the mate points at as one more candidate. The mate's window is
    /// searched directly (bit-parallel infix search, then strict extension),
    /// so a read whose every k-mer sits above the occurrence caps still
    /// reaches the locus its fragment belongs to — the way STAR and BWA
    /// resolve a repeat copy by the pair. The candidate is arbitrated with
    /// the read's own candidates under the same strict score; when it wins,
    /// the placement is marked rescued (MAPQ capped at the pileup gate, `RE`
    /// tag), because it was won on the mate's evidence rather than the read's
    /// own seeds. Returns None when the mate-window candidate does not win.
    pub fn align_read_with_mate(
        &mut self,
        seq: &[u8],
        mate_contig: u32,
        mate_pos: u32,
        cached: Option<&CachedPlacements>,
    ) -> Option<ReadAlignment> {
        let mut win = self.mate_window_placement(seq, mate_contig, mate_pos)?;
        let rc_storage = self.take_rc(seq);
        let q: &[u8] = if win.strand == Strand::Minus {
            &rc_storage
        } else {
            seq
        };
        self.rescore_strict(q, &mut win);
        let win_probe = win.clone();
        let (mut placements, n_seeds, _) = match cached {
            // Reuse the caller's own scored candidates: the read was just
            // aligned once already.
            Some(c) => (c.placements.clone(), c.n_seeds, c.genomic_any),
            None => {
                crate::prof::bump(crate::prof::MATE_REALIGNED, 1);
                self.scored_placements(seq, false)
                    .unwrap_or((Vec::new(), 0, false))
            }
        };
        placements.push(win);
        dedup_by_footprint(&mut placements);
        self.hyper_rearbitrate(seq, &rc_storage, &mut placements);
        // The mate-window candidate must have won the locus (dedup keeps the
        // better-scored of two same-locus candidates, so compare footprints).
        if footprint_overlap(&placements[0], &win_probe) < 0.5 {
            self.scratch.rc = rc_storage;
            return None;
        }
        self.scratch.rc = rc_storage;
        // Pair consistency prunes the competition: a runner-up on another
        // contig, or outside the fragment's reach of the mate, is not a
        // placement of this fragment. MAPQ is then the margin over the
        // mate-consistent alternatives only, which is how a paired aligner
        // resolves a repeat copy by its pair (STAR, BWA-MEM).
        let reach = 2 * (MATE_WINDOW_INSERT + MATE_WINDOW_INSERT / 2 + seq.len() as u32);
        placements.retain(|p| p.contig == mate_contig && p.pos.abs_diff(mate_pos) <= reach);
        // The genomic layer was consulted (the window search is genomic), so
        // the L1-only last-resort policy does not apply here.
        let mut aln = self.arbitrate(seq, placements, n_seeds, true)?;
        aln.rescued = true;
        Some(aln)
    }

    /// Candidate at the mate's locus: the read (both orientations) is placed
    /// in a window around the mate by the EA-Myers infix search and extended
    /// strictly on the diagonal it finds. No seeds are consulted, so the
    /// occurrence caps cannot hide the locus. None when the window holds no
    /// placement within the distance gate.
    fn mate_window_placement(
        &mut self,
        seq: &[u8],
        mate_contig: u32,
        mate_pos: u32,
    ) -> Option<Placement> {
        let read_len = seq.len() as u32;
        if read_len == 0 || read_len > 256 {
            return None;
        }
        let m = read_len as usize;
        let contig = &self.index.reference.contigs[mate_contig as usize];
        let span = MATE_WINDOW_INSERT + MATE_WINDOW_INSERT / 2 + read_len;
        let wstart = mate_pos.saturating_sub(span);
        let wend = mate_pos.saturating_add(span + read_len).min(contig.len);
        if wend <= wstart || ((wend - wstart) as usize) < m {
            return None;
        }
        let window = self.take_win(wstart, wend, contig);
        let rc_storage = self.take_rc(seq);
        let (mut myers_a, mut myers_b) = (
            std::mem::take(&mut self.scratch.myers_a),
            std::mem::take(&mut self.scratch.myers_b),
        );
        let threshold = ((read_len as i32) / 3).max(10);
        let mut best: Option<(Strand, i32, usize)> = None;
        for (strand, q) in [(Strand::Plus, seq), (Strand::Minus, &rc_storage[..])] {
            let (d, start) = if m <= 128 {
                crate::myers_ea::infix_best_start_into(q, &window, &mut myers_a, &mut myers_b)
            } else {
                crate::myers_ea::long::infix_best_start_into(q, &window, &mut myers_a, &mut myers_b)
            };
            if d > threshold || start + m > window.len() {
                continue;
            }
            if best.as_ref().is_none_or(|b| d < b.1) {
                best = Some((strand, d, start));
            }
        }
        self.scratch.myers_a = myers_a;
        self.scratch.myers_b = myers_b;
        let Some((strand, _, start)) = best else {
            self.scratch.rc = rc_storage;
            self.scratch.win = window;
            return None;
        };
        let q: &[u8] = match strand {
            Strand::Plus => seq,
            Strand::Minus => &rc_storage,
        };
        let hint = DiagHint {
            offset: start as i64,
            num: 1,
            den: 1,
        };
        crate::prof::dp_call(crate::prof::DP_MATE);
        let ext = extend_hint(q, &window, &self.config.extend, &mut self.buf, hint);
        self.scratch.rc = rc_storage;
        self.scratch.win = window;
        if ext.score <= 0 {
            return None;
        }
        Some(Placement::from_genomic(
            mate_contig,
            wstart + ext.ref_start,
            strand,
            ext.cigar,
            Vec::new(),
            0,
            ext.score,
        ))
    }

    /// Exhaustive re-alignment for the verify stage: every chained locus
    /// extended, decisive shortcuts off, the full strict-scored candidate
    /// list returned without arbitration or refusal.
    pub fn align_read_exhaustive(&mut self, seq: &[u8]) -> Vec<ScoredPlacement> {
        let Some((placements, _, _)) = self.scored_placements(seq, true) else {
            return Vec::new();
        };
        placements
            .into_iter()
            .map(|p| ScoredPlacement {
                blocks: ref_m_blocks(p.pos, &p.cigar),
                cigar: p.cigar,
                contig: p.contig,
                pos: p.pos,
                strand: p.strand,
                score: p.strict.score,
            })
            .collect()
    }

    /// The Tier-1 decision over an ordered, deduplicated candidate list.
    fn arbitrate(
        &mut self,
        seq: &[u8],
        placements: Vec<Placement>,
        n_seeds: usize,
        genomic_any: bool,
    ) -> Option<ReadAlignment> {
        if placements.is_empty() {
            return None;
        }
        if std::env::var_os("ESP_PROBE").is_some() {
            for (i, p) in placements.iter().enumerate() {
                eprintln!(
                    "[probe] cand {i} rtag={} origin={} contig={} pos={} strand={:?} strict={} mm={} ea={}/{} clip={} cigar={}",
                    String::from_utf8_lossy(&seq[..16.min(seq.len())]),
                    p.origin.label(),
                    self.index.reference.contigs[p.contig as usize].name,
                    p.pos,
                    p.strand,
                    p.strict.score,
                    p.strict.mm,
                    p.strict.ea_ag,
                    p.strict.ea_tc,
                    p.strict.clipped,
                    p.cigar
                        .iter()
                        .map(|op| match op {
                            CigarOp::Match(n) => format!("{n}M"),
                            CigarOp::Ins(n) => format!("{n}I"),
                            CigarOp::Del(n) => format!("{n}D"),
                            CigarOp::RefSkip(n) => format!("{n}N"),
                            CigarOp::SoftClip(n) => format!("{n}S"),
                        })
                        .collect::<String>()
                );
            }
        }
        let best = &placements[0];
        // Runner-up must be a different locus (footprint overlap < 0.5).
        let second = placements[1..]
            .iter()
            .filter(|p| footprint_overlap(best, p) < 0.5)
            .map(|p| p.strict.score)
            .filter(|s| *s > 0)
            .max()
            .unwrap_or(0);

        // 5. Refusal: a placement that cannot explain two thirds of the read
        // is not a placement — the read belongs in the unmapped pool, where
        // the Tier-2 collapsed-alphabet rescue looks for it.
        if f64::from(best.strict.aligned) < REFUSE_MIN_ALIGNED_FRAC * seq.len() as f64 {
            return None;
        }

        // A best hit on a decoy contig is the decoy contract doing its job:
        // the read's source is absent from the primary reference (rDNA,
        // retrotransposed copies, contamination), and parking it on the decoy
        // keeps it out of the pileup while remaining visible and auditable in
        // the BAM. Decoy placements never enter the Tier-2 pool.
        let best_origin = if is_decoy_contig(self.index.reference, best.contig) {
            crate::mapq::Origin::Decoy
        } else {
            best.origin
        };

        let rc_storage = self.take_rc(seq);
        let q: &[u8] = if best.strand == Strand::Minus {
            &rc_storage
        } else {
            seq
        };
        let n_anchors = self.exact_kmer_support(q, best);

        let mut aln = ReadAlignment {
            contig: best.contig,
            pos: best.pos,
            strand: best.strand,
            score: best.strict.score,
            second_score: second,
            chain_score: if best.origin == crate::mapq::Origin::Genomic {
                best.chain_score
            } else {
                best.strict.score
            },
            second_chain_score: 0,
            cigar: best.cigar.clone(),
            n_anchors,
            junctions: best.junctions.clone(),
            ea_ag: best.strict.ea_ag,
            ea_tc: best.strict.ea_tc,
            mm_count: best.strict.mm,
            n_seeds,
            rescued: false,
            mapq_override: None,
            origin: best_origin,
            n_cands: placements.len(),
        };

        // Only the transcriptome layer could place the read. Its apparent
        // ambiguity has a cause when the read carries a hyperediting
        // signature — dense single-direction A>G (or T>C) substitution
        // destroys the exact-seed discrimination the genomic layer runs on,
        // which is the whole reason the Tier-2 rescue exists. Such a
        // placement takes the pileup-gate ceiling; anything else the
        // transcriptome layer alone placed is a last resort and scores 0.
        if !genomic_any && aln.origin == crate::mapq::Origin::L1 {
            if aln.hyperedited_signature() {
                aln.mapq_override = Some(crate::mapq::PILEUP_GATE);
            } else {
                crate::mapq::force_mapq0(&mut aln);
            }
        }
        self.scratch.rc = rc_storage;

        Some(aln)
    }

    /// Rescore one candidate in the strict domain: per-base match +2 /
    /// mismatch −4 (A>G and T>C are charged like any mismatch and counted
    /// separately for the density gates), affine gaps `−(open + ext·n)`,
    /// introns by splice signal plus known-junction bonus minus the intron
    /// length penalty, soft clips free. N on either side scores 0 and counts
    /// as neither match nor mismatch. Identical constants to the extender
    /// and splicer, so scores are comparable across every producing path.
    fn rescore_strict(&mut self, query: &[u8], p: &mut Placement) {
        self.rescore_with(query, p, EditDir::Off);
    }

    /// [`rescore_strict`] with one editing direction exempt: substitutions of
    /// that direction score 0 instead of `mismatch` (they are still counted
    /// in `ea_ag`/`ea_tc`). This is the reduced-alphabet scoring of the
    /// hyperediting literature (Porath, Carmi & Levanon, Nat Commun 5:4726,
    /// 2014) applied to placement arbitration: a read whose strict-domain
    /// winner carries a hyperediting signature is re-arbitrated in the
    /// alphabet its own edits define, so that its true locus is not out-
    /// scored by a paralog that merely lacks the edits, and its edited tail
    /// is not soft-clipped away. The four-letter `unexplained` density cap in
    /// `mapq_of` still applies, so a paralog exploiting both directions or
    /// ordinary mismatches is not rescued by this.
    fn rescore_with(&mut self, query: &[u8], p: &mut Placement, dir: EditDir) {
        let ep = self.config.extend;
        let sp = SpliceParams::default();
        let contig = &self.index.reference.contigs[p.contig as usize];
        let empty_lib = JunctionLib::default();
        let lib = self.jlib.as_deref().unwrap_or(&empty_lib);
        let mut s = StrictScore {
            ..Default::default()
        };
        let mut score: i32 = 0;
        let mut ro = 0u32;
        let mut rf = p.pos;
        let mut scratch = std::mem::take(&mut self.scratch.rescore);
        let base = |i: u32| -> u8 {
            if i < contig.len {
                contig.base(i).to_ascii()
            } else {
                b'N'
            }
        };
        for op in &p.cigar {
            match op {
                CigarOp::Match(n) => {
                    let avail = (*n).min(contig.len.saturating_sub(rf));
                    scratch.clear();
                    if avail > 0 {
                        contig.decode_append(rf, rf + avail, &mut scratch);
                    }
                    for k in 0..*n as usize {
                        let qb = query.get(ro as usize + k).copied().unwrap_or(b'N');
                        let rb = if k < scratch.len() {
                            scratch[k]
                        } else {
                            b'N'
                        };
                        if qb == b'N' || rb == b'N' {
                            continue; // 0 score, uncounted
                        }
                        s.aligned += 1;
                        if qb == rb {
                            score += ep.match_score;
                            s.matches += 1;
                        } else {
                            if !dir.exempt(rb, qb) {
                                score += ep.mismatch;
                            }
                            match crate::extend::exempt_dir(rb, qb) {
                                Some(true) => s.ea_ag += 1,
                                Some(false) => s.ea_tc += 1,
                                None => s.mm += 1,
                            }
                        }
                    }
                    ro += *n;
                    rf += *n;
                }
                CigarOp::Ins(n) => {
                    score -= ep.gap_open + ep.gap_ext * *n as i32;
                    s.aligned += *n;
                    s.ins += 1;
                    ro += *n;
                }
                CigarOp::Del(n) => {
                    score -= ep.gap_open + ep.gap_ext * *n as i32;
                    s.aligned += *n;
                    s.del += 1;
                    rf += *n;
                }
                CigarOp::SoftClip(n) => {
                    s.clipped += *n;
                    ro += *n;
                }
                CigarOp::RefSkip(n) => {
                    // An annotated intron earns the known bonus alone,
                    // regardless of the read's orientation: a read
                    // antisense to the transcript crosses the same
                    // annotated intron. Signal score and intron-length
                    // penalty exist to discriminate *novel* junctions from
                    // spurious ones — charging them on annotated junctions
                    // lets a one-base soft clip beat a real splice across
                    // every long or noncanonical intron (measured on the
                    // simulated corpus: junction-exact 92.8% vs the 96%
                    // gate; the loss concentrated exactly there).
                    if lib.contains_any_strand(p.contig, rf, rf + *n) {
                        score += sp.known_bonus;
                    } else {
                        // Donor/acceptor are the first and last two bases
                        // OF the intron (GT..AG), matching `splice::dinuc`'s
                        // convention. The transcript's strand is unknown
                        // here — an unstranded library sequences half of
                        // its reads antisense — so a canonical motif on
                        // either strand counts and the read's own strand
                        // is not consulted.
                        let donor = [base(rf), base(rf.saturating_add(1))];
                        let acc = [
                            base(rf + n.saturating_sub(2)),
                            base(rf + n.saturating_sub(1)),
                        ];
                        let signal = splice_signal(donor, acc, false)
                            .score(&sp)
                            .max(splice_signal(donor, acc, true).score(&sp));
                        score += signal - intron_length_penalty(*n);
                    }
                    s.n_splice += 1;
                    rf += *n;
                }
            }
        }
        s.score = score;
        p.strict = s;
        p.ea_dir = dir;
        self.scratch.rescore = scratch;
    }

    /// Re-extend a soft-clipped genomic candidate with one editing direction
    /// exempt, over the window its current footprint plus clips define.
    /// Replaces `pos`/`cigar` when the new extension aligns more of the read.
    /// Spliced candidates are left alone (the window would have to span the
    /// intron; their clips are refined elsewhere).
    fn reextend_with(&mut self, query: &[u8], p: &mut Placement, dir: EditDir) {
        if p.cigar.iter().any(|op| matches!(op, CigarOp::RefSkip(_))) {
            return;
        }
        let read_len = query.len() as u32;
        let (mut lead, mut tail, mut ref_span, mut aligned) = (0u32, 0u32, 0u32, 0u32);
        let mut seen_m = false;
        for op in &p.cigar {
            match op {
                CigarOp::SoftClip(n) => {
                    if seen_m {
                        tail += *n;
                    } else {
                        lead += *n;
                    }
                }
                CigarOp::Match(n) => {
                    seen_m = true;
                    ref_span += *n;
                    aligned += *n;
                }
                CigarOp::Ins(n) => aligned += *n,
                CigarOp::Del(n) => ref_span += *n,
                CigarOp::RefSkip(_) => {}
            }
        }
        let contig = &self.index.reference.contigs[p.contig as usize];
        let wstart = p.pos.saturating_sub(lead + 30);
        let wend = (p.pos + ref_span + tail + 30).min(contig.len);
        if wend <= wstart {
            return;
        }
        let window = self.take_win(wstart, wend, contig);
        let hint = DiagHint {
            offset: (p.pos - wstart) as i64 - i64::from(lead),
            num: 1,
            den: 1,
        };
        let mut ep = self.config.extend;
        ep.edit_dir = dir;
        crate::prof::dp_call(crate::prof::DP_HYPER);
        let ext = extend_hint(query, &window, &ep, &mut self.buf, hint);
        self.scratch.win = window;
        let new_aligned: u32 = ext
            .cigar
            .iter()
            .map(|op| match op {
                CigarOp::Match(n) | CigarOp::Ins(n) => *n,
                _ => 0,
            })
            .sum();
        if ext.score > 0 && new_aligned >= aligned && new_aligned <= read_len {
            p.pos = wstart + ext.ref_start;
            p.cigar = ext.cigar;
            p.ext_score = ext.score;
        }
    }

    /// Count non-overlapping exact k-mers along the placement's M blocks
    /// (k from the seed config). Feeds `pen_cm` in the MAPQ formula for L1
    /// and genomic candidates alike — an anchor count measured on the
    /// placement, not borrowed from a chain.
    fn exact_kmer_support(&mut self, query: &[u8], p: &Placement) -> usize {
        let k = self.config.seed.k as usize;
        if k == 0 || query.is_empty() {
            return 0;
        }
        let contig = &self.index.reference.contigs[p.contig as usize];
        let mut count = 0usize;
        let mut ro = 0usize;
        let mut rf = p.pos;
        let mut scratch = std::mem::take(&mut self.scratch.kmer_scan);
        for op in &p.cigar {
            match op {
                CigarOp::Match(n) => {
                    let n = *n as usize;
                    let avail = (n as u32).min(contig.len.saturating_sub(rf));
                    scratch.clear();
                    if avail > 0 {
                        contig.decode_append(rf, rf + avail, &mut scratch);
                    }
                    let mut i = 0usize;
                    while i + k <= n {
                        let mut ok = true;
                        for t in 0..k {
                            let qb = query.get(ro + i + t).copied().unwrap_or(b'N');
                            let rb = if i + t < scratch.len() {
                                scratch[i + t]
                            } else {
                                b'N'
                            };
                            // A k-mer is exact in the alphabet the candidate
                            // was scored in: for a hyperedited read that is
                            // its single-direction reduced alphabet, where its
                            // edits are identities. Counting four-letter exact
                            // k-mers on such a read gives structurally ~0
                            // anchors (30 edits in 150 bp break every 15-mer)
                            // and hence MAPQ 0 whatever the margin.
                            if qb != rb && !p.ea_dir.exempt(rb, qb) {
                                ok = false;
                                break;
                            }
                        }
                        if ok {
                            count += 1;
                            i += k;
                        } else {
                            i += 1;
                        }
                    }
                    ro += n;
                    rf += n as u32;
                }
                CigarOp::Ins(n) | CigarOp::SoftClip(n) => ro += *n as usize,
                CigarOp::Del(n) | CigarOp::RefSkip(n) => rf += *n,
            }
        }
        self.scratch.kmer_scan = scratch;
        count
    }


    /// Genomic layer driver, in candidate form: seed → anchor → chain →
    /// build a locus-deduplicated chain list → extend each chain
    /// (intron-chain / spliced / EA-Myers / fast lane / banded DP) into
    /// 0..2 placements → split/span rescue on the best placement. The
    /// single-winner logic of the old path is gone: every defensible locus
    /// is reported and the arbiter in [`Aligner::align_read`] decides.
    fn genomic_candidates(&mut self, seq: &[u8], exhaustive: bool) -> GenomicBatch {
        let cfg = self.config;
        let k = cfg.seed.k;
        let read_len = seq.len() as u32;

        let mut mins = std::mem::take(&mut self.scratch.mins_out);
        minimizers_into(seq, cfg.seed, &mut self.scratch.seed, &mut mins);
        if mins.is_empty() {
            self.scratch.mins_out = mins;
            return GenomicBatch {
                placements: Vec::new(),
                n_seeds: 0,
            };
        }

        // Round 1: default-density seeding, occurrence cap = freq_cutoff.
        let r1 = self.chain_round(&mins, read_len, k, self.index.freq_cutoff);
        if std::env::var_os("ESP_PROBE").is_some() {
            eprintln!(
                "[probe] mins={} r1_best={:?}",
                mins.len(),
                r1.as_ref().and_then(|r| r.best()).map(|b| b.score)
            );
        }

        // Round 2: mid-occurrence seeding when round 1 explains little of
        // the read. A read from a high-copy repeat has every k-mer above the
        // default occurrence cap, so round 1 sees only its rare (mutated)
        // k-mers and chains them at whichever paralog they happen to hit —
        // measured on the real sample, the true locus of such a read was
        // absent from the candidate list for nearly every read STAR placed
        // uniquely and esperanto placed elsewhere. minimap2 runs the same
        // second round (`--mid-occ`); anchors are capped so a read made of
        // nothing but Alu k-mers cannot turn chaining quadratic.
        let r1_cov = r1
            .as_ref()
            .and_then(|r| r.best())
            .map(|b| {
                let (qs, qe) = b.read_span(k);
                f64::from(qe.saturating_sub(qs)) / f64::from(read_len.max(1))
            })
            .unwrap_or(0.0);
        let r2 = if r1_cov < MID_ROUND_MIN_COVERAGE {
            let mut anchors = std::mem::take(&mut self.scratch.anchors);
            collect_anchors_into(self.index, &mins, read_len, k, MID_OCC_CAP, &mut anchors);
            if anchors.is_empty() || anchors.len() > MID_ROUND_MAX_ANCHORS {
                self.scratch.anchors = anchors;
                None
            } else {
                let mut chains = std::mem::take(&mut self.scratch.chains_out);
                chain_anchors_into(
                    &mut anchors,
                    &cfg.chain,
                    &mut chains,
                    &mut self.scratch.chain,
                );
                self.scratch.anchors = anchors;
                let second = second_score_with(&chains, &mut self.scratch.second_keys);
                Some(RoundOutcome { chains, second })
            }
        } else {
            None
        };

        // Editing-aware variant round: only when chaining fully failed (weak).
        let weak = r1
            .as_ref()
            .and_then(|r| r.best())
            .map(|b| b.score)
            .unwrap_or(0)
            == 0;
        let r3 = if cfg.extend.edit_dir.is_on() && weak {
            // Base anchors first (see the note in the old single-winner path:
            // edit-variant anchors must not pollute a clean chain).
            let mut base_anchors = std::mem::take(&mut self.scratch.anchors);
            collect_anchors_into(
                self.index,
                &mins,
                read_len,
                k,
                MID_OCC_CAP,
                &mut base_anchors,
            );
            if !base_anchors.is_empty() {
                let mut cp = cfg.chain;
                cp.min_chain_score = 25;
                let mut chains = std::mem::take(&mut self.scratch.chains_out);
                chain_anchors_into(&mut base_anchors, &cp, &mut chains, &mut self.scratch.chain);
                self.scratch.anchors = base_anchors;
                let second = second_score_with(&chains, &mut self.scratch.second_keys);
                Some(RoundOutcome { chains, second })
            } else {
                self.scratch.anchors = base_anchors;
                let (mut anchors, _var_hits) = crate::index::collect_anchors_edit_variants(
                    self.index,
                    &mins,
                    read_len,
                    k,
                    MID_OCC_CAP,
                );
                if anchors.is_empty() {
                    None
                } else {
                    let mut cp = cfg.chain;
                    cp.min_chain_score = 25;
                    let mut chains = std::mem::take(&mut self.scratch.chains_out);
                    chain_anchors_into(&mut anchors, &cp, &mut chains, &mut self.scratch.chain);
                    let second = second_score_with(&chains, &mut self.scratch.second_keys);
                    Some(RoundOutcome { chains, second })
                }
            }
        } else {
            None
        };

        // Merge rounds: winner = round with the highest best-chain score.
        let rounds: Vec<&RoundOutcome> = r1.iter().chain(r2.iter()).chain(r3.iter()).collect();
        let mut win_idx: Option<usize> = None;
        for (i, r) in rounds.iter().enumerate() {
            let better = match (win_idx, r.best()) {
                (None, Some(_)) => true,
                (Some(w), Some(b)) => b.score > rounds[w].best().map(|x| x.score).unwrap_or(0),
                _ => false,
            };
            if better {
                win_idx = Some(i);
            }
        }
        let Some(w) = win_idx else {
            let n_seeds = mins.len();
            self.scratch.mins_out = mins;
            return GenomicBatch {
                placements: Vec::new(),
                n_seeds,
            };
        };
        let cr = rounds[w];
        let second = cr.second;
        let extra_top: Vec<Chain> = rounds
            .iter()
            .flat_map(|r| r.chains.iter().take(TOP_CANDIDATES).cloned())
            .collect();
        let best = match cr.best() {
            Some(b) => b.clone(),
            None => {
                let n_seeds = mins.len();
                self.scratch.mins_out = mins;
                return GenomicBatch {
                    placements: Vec::new(),
                    n_seeds,
                }
            }
        };
        if std::env::var_os("ESP_PROBE").is_some() {
            eprintln!(
                "[probe] best chain: contig={} strand={:?} score={} anchors={} second={} tied={}",
                self.index.reference.contigs[best.contig as usize].name,
                best.strand,
                best.score,
                best.anchors.len(),
                second,
                tied_best(&cr.chains).len()
            );
        }

        // Candidate chains: deduplicated by locus (same contig+strand,
        // ref-span start within 50bp), in chain-score order.
        let mut cand_chains: Vec<Chain> = Vec::new();
        let push_dedup = |c: Chain, out: &mut Vec<Chain>| {
            let dup = out.iter().any(|x| {
                x.contig == c.contig
                    && x.strand == c.strand
                    && x.ref_span(k).0.abs_diff(c.ref_span(k).0) < 50
            });
            if !dup {
                out.push(c);
            }
        };
        for c in cr.chains.iter().take(TOP_CANDIDATES) {
            push_dedup(c.clone(), &mut cand_chains);
        }
        for c in extra_top {
            if cand_chains.len() >= MAX_TIED_CANDIDATES {
                break;
            }
            push_dedup(c, &mut cand_chains);
        }

        let tied = tied_best(&cr.chains);
        let decisive = tied.len() <= 1 && second * 5 < best.score * 4;
        let rc_storage = self.take_rc(seq);
        let mut placements: Vec<Placement> = Vec::new();
        let mut dp_budget = EXTEND_BUDGET;
        let width = if exhaustive {
            TOP_CANDIDATES
        } else {
            MAX_TIED_CANDIDATES
        };
        for (i, chain) in cand_chains.iter().take(width).enumerate() {
            let is_best = i == 0;
            self.extend_chain(
                seq,
                &mins,
                chain,
                is_best,
                decisive && is_best && !exhaustive,
                &rc_storage,
                &mut dp_budget,
                &mut placements,
            );
        }
        self.scratch.rc = rc_storage;
        if let Some(r1) = r1 {
            self.scratch.chains_out = r1.chains;
        }
        // Split/span rescue and micro-deletion conversion on the best
        // genomic placement (by extension score, chain order on ties) — the
        // same post-processing the single-winner path applied to its winner.
        let mut best_idx: Option<usize> = None;
        for (i, p) in placements.iter().enumerate() {
            let better = match best_idx {
                None => true,
                Some(j) => {
                    p.ext_score > placements[j].ext_score
                        || (p.ext_score == placements[j].ext_score
                            && p.chain_score > placements[j].chain_score)
                }
            };
            if better {
                best_idx = Some(i);
            }
        }
        if let Some(i) = best_idx {
            let p = placements[i].clone();
            let tmp = ReadAlignment {
                contig: p.contig,
                pos: p.pos,
                strand: p.strand,
                score: p.ext_score,
                chain_score: p.chain_score,
                cigar: p.cigar,
                junctions: p.junctions,
                ..Default::default()
            };
            let tmp = self.apply_split_rescue(seq, tmp);
            let tmp = self.apply_span_rescue(seq, tmp);
            let tmp = self.convert_annotated_micro_deletions(tmp);
            let q = &mut placements[i];
            q.pos = tmp.pos;
            q.cigar = tmp.cigar;
            q.junctions = tmp.junctions;
        }

        let n_seeds = mins.len();
        self.scratch.mins_out = mins;

        GenomicBatch {
            placements,
            n_seeds,
        }
    }

    /// Extend one chain into 0..2 placement candidates: the spliced
    /// interpretation first (RNA), then the unspliced one (EA-Myers window
    /// extend, with the banded DP as indel fallback while budget lasts).
    /// `fast_lane_ok` gates the Kadane fast lane (decisive best chain only).
    #[allow(clippy::too_many_arguments)]
    fn extend_chain(
        &mut self,
        seq: &[u8],
        mins: &[Minimizer],
        chain: &Chain,
        is_best: bool,
        fast_lane_ok: bool,
        rc_storage: &[u8],
        dp_budget: &mut usize,
        out: &mut Vec<Placement>,
    ) {
        let cfg = self.config;
        let k = cfg.seed.k;
        let read_len = seq.len() as u32;

        // Spliced interpretation (RNA; splice gate frozen on: <2 segments
        // has no spliced reading).
        if cfg.rna {
            let min_intron = SpliceParams::default().min_intron;
            let segs = segment_chain(chain, k, min_intron);
            if segs.len() == 2 && is_best {
                // Intron-chain fast path: verified intra-chain intron
                // algorithm, with read-aware boundary refinement.
                let empty_lib = JunctionLib::default();
                let lib = self.jlib.as_deref().unwrap_or(&empty_lib);
                let query: &[u8] = match chain.strand {
                    Strand::Plus => seq,
                    // Same bytes as a fresh revcomp(seq): the caller passes
                    // exactly that buffer.
                    Strand::Minus => rc_storage,
                };
                let mut ic_anchors = std::mem::take(&mut self.scratch.anchors);
                collect_anchors_into(
                    self.index,
                    mins,
                    read_len,
                    k,
                    MID_OCC_CAP,
                    &mut ic_anchors,
                );
                ic_anchors.retain(|a| a.contig == chain.contig && a.strand == chain.strand);
                let ic_hit = crate::intron_chain::try_intron_chain_placement(
                    self.index.reference,
                    lib,
                    query,
                    &ic_anchors,
                    k,
                    &crate::intron_chain::IntronParams::default(),
                    &cfg.extend,
                    &mut self.buf,
                );
                self.scratch.anchors = ic_anchors;
                if let Some(hit) = ic_hit {
                    let mut placement = Placement::from_genomic(
                        chain.contig,
                        hit.genomic_start(),
                        chain.strand,
                        hit.cigar,
                        hit.junctions,
                        chain.score,
                        hit.extension.score,
                    );
                    // Read-aware boundary refinement (±3bp snap).
                    if !placement.junctions.is_empty() {
                        let minus = matches!(chain.strand, Strand::Minus);
                        let sp = SpliceParams::default();
                        if let Some(refined) = refine_junction_read(
                            self.index.reference,
                            query,
                            placement.contig,
                            segs[0].r_end,
                            segs[1].r_start,
                            segs[0].q_end,
                            segs[1].q_start,
                            minus,
                            lib,
                            &sp,
                        ) {
                            let new_intron =
                                refined.junction.end.saturating_sub(refined.junction.start);
                            for op in &mut placement.cigar {
                                if let CigarOp::RefSkip(n) = op {
                                    *n = new_intron;
                                }
                            }
                            placement.junctions[0] = refined;
                        }
                    }
                    out.push(placement);
                    return;
                }
            }
            if segs.len() >= 2 {
                if let Some(hit) = self.try_spliced(seq, rc_storage, chain) {
                    out.push(Placement::from_genomic(
                        chain.contig,
                        hit.pos,
                        chain.strand,
                        hit.cigar,
                        hit.junctions,
                        chain.score,
                        hit.score,
                    ));
                    return;
                }
            }
        }

        // Unspliced interpretation — fast lane (decisive best chain only:
        // zero-mismatch Kadane, point-equivalent to the affine DP).
        if fast_lane_ok {
            if let Some((ext, wstart)) = self.try_fast_lane(seq, rc_storage, chain, read_len, k) {
                if ext.score > 0 {
                    crate::prof::bump(crate::prof::FAST_LANE_HITS, 1);
                    out.push(Placement::from_genomic(
                        chain.contig,
                        wstart + ext.ref_start,
                        chain.strand,
                        ext.cigar,
                        Vec::new(),
                        chain.score,
                        ext.score,
                    ));
                    return;
                }
            }
        }

        // Unspliced interpretation — EA-Myers window extend (≤256bp; the
        // bit-parallel infix search with Kadane-optimal soft clips).
        let mut ea_done_here = false;
        if read_len > 0 && read_len <= 256 {
            let dir = cfg.extend.edit_dir;
            let m = read_len as usize;
            let match_s = cfg.extend.match_score;
            let mm_s = cfg.extend.mismatch;
            let threshold = ((read_len as i32) / 3).max(10);
            let contig = &self.index.reference.contigs[chain.contig as usize];
            let flank = read_len / 2 + 30;
            let (rs, re) = chain.ref_span(k);
            let wstart = rs.saturating_sub(flank);
            let wcap = rs + read_len + flank;
            let wend = (re + flank).min(contig.len).min(wcap);
            let window = self.take_win(wstart, wend, contig);
            let q: &[u8] = match chain.strand {
                Strand::Plus => seq,
                Strand::Minus => rc_storage,
            };
            if window.len() >= m {
                let (mut myers_a, mut myers_b) = (
                    std::mem::take(&mut self.scratch.myers_a),
                    std::mem::take(&mut self.scratch.myers_b),
                );
                let (d, start) = if m <= 128 {
                    crate::myers_ea::infix_best_start_into(q, &window, &mut myers_a, &mut myers_b)
                } else {
                    crate::myers_ea::long::infix_best_start_into(
                        q,
                        &window,
                        &mut myers_a,
                        &mut myers_b,
                    )
                };
                self.scratch.myers_a = myers_a;
                self.scratch.myers_b = myers_b;
                let mut usable = d <= threshold && start + m <= window.len();
                if usable {
                    // Per-base verification: mm must exactly explain d,
                    // otherwise an indel is indicated.
                    let mut mm_full: i32 = 0;
                    let mut has_n = false;
                    for kk in 0..m {
                        let qb = q[kk];
                        let rb = window[start + kk];
                        if qb == b'N' || rb == b'N' {
                            has_n = true;
                            break;
                        }
                        if qb != rb && crate::extend::exempt_dir(rb, qb).is_none() {
                            mm_full += 1;
                        }
                    }
                    if has_n {
                        usable = false;
                    } else if mm_full != d {
                        usable = false; // indel implied → DP below
                    }
                }
                if usable {
                    // Optimal soft clip ≡ Kadane max subarray on the diagonal.
                    crate::prof::bump(crate::prof::MYERS_ONLY_HITS, 1);
                    let mut best_sum = 0i32;
                    let (mut bi, mut bj) = (0usize, 0usize);
                    let mut cur_sum = 0i32;
                    let mut cur_i = 0usize;
                    for kk in 0..m {
                        let qb = q[kk];
                        let rb = window[start + kk];
                        let s = if qb == rb {
                            match_s
                        } else if dir.exempt(rb, qb) {
                            0
                        } else {
                            mm_s
                        };
                        if cur_sum <= 0 {
                            cur_sum = s;
                            cur_i = kk;
                        } else {
                            cur_sum += s;
                        }
                        if cur_sum > best_sum {
                            best_sum = cur_sum;
                            bi = cur_i;
                            bj = kk + 1;
                        }
                    }
                    let seg_len = (bj - bi) as u32;
                    if best_sum > 0 && seg_len * 2 >= read_len {
                        let mut cigar = Vec::with_capacity(3);
                        if bi > 0 {
                            cigar.push(CigarOp::SoftClip(bi as u32));
                        }
                        cigar.push(CigarOp::Match(seg_len));
                        if m - bj > 0 {
                            cigar.push(CigarOp::SoftClip((m - bj) as u32));
                        }
                        let pos = wstart + (start + bi) as u32;
                        out.push(Placement::from_genomic(
                            chain.contig,
                            pos,
                            chain.strand,
                            cigar,
                            Vec::new(),
                            chain.score,
                            best_sum,
                        ));
                        ea_done_here = true;
                    }
                }
            }
            self.scratch.win = window;
        }
        if ea_done_here {
            return;
        }
        if ea_done_here {
            return;
        }

        // Unspliced interpretation — banded DP (indel fallback and long
        // reads), within the shared extension budget.
        if *dp_budget == 0 {
            return;
        }
        *dp_budget -= 1;
        let contig = &self.index.reference.contigs[chain.contig as usize];
        let flank = read_len / 2 + 30;
        let (rs, re) = chain.ref_span(k);
        let wstart = rs.saturating_sub(flank);
        let wcap = rs + read_len + flank;
        let wend = (re + flank).min(contig.len).min(wcap);
        let window = self.take_win(wstart, wend, contig);
        let q: &[u8] = match chain.strand {
            Strand::Plus => seq,
            Strand::Minus => rc_storage,
        };
        let (qs, _) = chain.read_span(k);
        let hint = DiagHint {
            offset: (rs - wstart) as i64 - qs as i64,
            num: 1,
            den: 1,
        };
        crate::prof::dp_call(crate::prof::DP_UNSPLICED);
        let ext = extend_hint(q, &window, &cfg.extend, &mut self.buf, hint);
        self.scratch.win = window;
        if ext.score > 0 {
            out.push(Placement::from_genomic(
                chain.contig,
                wstart + ext.ref_start,
                chain.strand,
                ext.cigar,
                Vec::new(),
                chain.score,
                ext.score,
            ));
        }
    }

    /// Splice path: segment → junction refine → pseudo-reference extend →
    /// N CIGAR. splice gate frozen on: <2 segments ⇒ None (align_spliced would
    /// return None anyway — skip the revcomp + entry cost).
    fn try_spliced(&mut self, seq: &[u8], rc_storage: &[u8], best: &Chain) -> Option<ReadAlignment> {
        let min_intron = SpliceParams::default().min_intron;
        let segs = segment_chain(best, self.config.seed.k, min_intron);
        if segs.len() < 2 {
            return None;
        }
        let empty_lib = JunctionLib::default();
        let lib = self.jlib.as_deref().unwrap_or(&empty_lib);
        // rc_storage is the caller's revcomp(seq) buffer (same bytes).
        let query: &[u8] = match best.strand {
            Strand::Plus => seq,
            Strand::Minus => rc_storage,
        };
        let hit = align_spliced(
            self.index.reference,
            lib,
            query,
            best,
            self.config.seed.k,
            &SpliceParams::default(),
            &self.config.extend,
            &mut self.buf,
        )?;
        if hit.extension.score <= 0 {
            return None;
        }
        let mut aln = ReadAlignment {
            contig: best.contig,
            pos: hit.genomic_start(),
            strand: best.strand,
            score: hit.extension.score,
            chain_score: best.score,
            second_chain_score: 0,
            cigar: hit.cigar,
            n_anchors: best.anchors.len(),
            junctions: hit.junctions,
            ..Default::default()
        };
        self.recount_mm_ea(query, &mut aln);
        Some(aln)
    }

    /// Rescue channel: mate anchored, A/G-masking reseeding in the anchored
    /// window. Success → rescued=true (RE tag).
    pub fn rescue_with_mate_anchor(
        &mut self,
        seq: &[u8],
        mate: &ReadAlignment,
        est_insert: u32,
    ) -> Option<ReadAlignment> {
        let read_len = seq.len() as u32;
        let mut mins = std::mem::take(&mut self.scratch.mins_out);
        minimizers_into(seq, self.config.seed, &mut self.scratch.seed, &mut mins);
        if mins.is_empty() {
            self.scratch.mins_out = mins;
            return None;
        }
        let k = self.config.seed.k;
        // Anchor window: mate pos ± (est_insert + 3×est_insert/2 + read_len).
        let span = est_insert + est_insert / 2 + read_len;
        let wlo = mate.pos.saturating_sub(span);
        let whi = mate.pos.saturating_add(span) + read_len;
        // Normal + single-edit-variant anchors (A/G masking), window-filtered.
        let (mut anchors, _) =
            crate::index::collect_anchors_edit_variants(self.index, &mins, read_len, k, 500);
        let mut base = std::mem::take(&mut self.scratch.anchors);
        collect_anchors_into(self.index, &mins, read_len, k, 500, &mut base);
        anchors.append(&mut base);
        self.scratch.anchors = base;
        let n_seeds = mins.len();
        anchors.retain(|a| a.contig == mate.contig && a.rpos >= wlo && a.rpos <= whi);
        if anchors.len() < 3 {
            return None;
        }
        let mut cp = self.config.chain;
        cp.min_chain_score = 20; // lower chaining threshold for rescue
        let mut chains = std::mem::take(&mut self.scratch.chains_out);
        chain_anchors_into(&mut anchors, &cp, &mut chains, &mut self.scratch.chain);
        self.scratch.mins_out = mins;
        let second = second_score_with(&chains, &mut self.scratch.second_keys);
        let best = chains.into_iter().next()?;
        // EA-scored extension (window = chain span ± read_len).
        let contig = &self.index.reference.contigs[best.contig as usize];
        let (rs, re) = best.ref_span(k);
        let flank = read_len + 50;
        let wstart = rs.saturating_sub(flank);
        let wend = (re + flank).min(contig.len).min(wstart + MAX_EXTEND_WINDOW);
        let window = self.take_win(wstart, wend, contig);
        let rc_storage = self.take_rc(seq);
        let query: &[u8] = match best.strand {
            Strand::Plus => seq,
            Strand::Minus => &rc_storage,
        };
        let hint = DiagHint {
            offset: (rs - wstart) as i64 - best.read_span(k).0 as i64,
            num: 1,
            den: 1,
        };
        // The mate window has already fixed the locus, so an editing exemption
        // here only has to explain the read's substitutions — but it still has
        // to pick a direction, since a read is one transcript and exempting
        // both at once would score a paralogous copy as an edited one. Extend
        // once per direction and keep the better explanation.
        let mut ep = self.config.extend;
        ep.edit_dir = EditDir::Ag;
        crate::prof::dp_call(crate::prof::DP_MATE);
        let ag = extend_hint(query, &window, &ep, &mut self.buf, hint);
        ep.edit_dir = EditDir::Tc;
        crate::prof::dp_call(crate::prof::DP_MATE);
        let tc = extend_hint(query, &window, &ep, &mut self.buf, hint);
        self.scratch.win = window;
        let ext = if tc.score > ag.score { tc } else { ag };
        if ext.score <= 0 {
            self.scratch.rc = rc_storage;
            return None;
        }
        let mut aln = ReadAlignment {
            contig: best.contig,
            pos: wstart + ext.ref_start,
            strand: best.strand,
            score: ext.score,
            chain_score: best.score,
            second_chain_score: second,
            cigar: ext.cigar,
            n_anchors: best.anchors.len(),
            n_seeds,
            rescued: true,
            ..Default::default()
        };
        self.recount_mm_ea(query, &mut aln);
        self.scratch.rc = rc_storage;
        Some(aln)
    }

    /// Split rescue (both tails; success rewrites CIGAR/junctions/pos).
    fn apply_split_rescue(&mut self, seq: &[u8], aln: ReadAlignment) -> ReadAlignment {
        let empty_lib = JunctionLib::default();
        // Query in extension orientation: (RC'd if minus) read; tails operate
        // in extension coordinates; result CIGAR is plus-reference isomorphic.
        let rc_storage = self.take_rc(seq);
        let lib = self.jlib.as_deref().unwrap_or(&empty_lib);
        let query: &[u8] = match aln.strand {
            Strand::Plus => seq,
            Strand::Minus => &rc_storage,
        };
        let mut rs = 0u32;
        let mut re = query.len() as u32;
        let mut seen_body = false;
        for op in &aln.cigar {
            match op {
                CigarOp::SoftClip(n) => {
                    if !seen_body {
                        rs = *n;
                    } else {
                        re -= *n;
                    }
                }
                _ => seen_body = true,
            }
        }
        let read_start = rs;
        let read_end = re;
        let mut out = aln;
        let primary_ref_end = {
            let mut e = out.pos;
            for op in &out.cigar {
                match op {
                    CigarOp::Match(n) | CigarOp::Del(n) | CigarOp::RefSkip(n) => e += *n,
                    _ => {}
                }
            }
            e
        };
        let params = SpliceParams::default();
        // Right tail = high-coordinate side (downstream intron), left tail =
        // low side; same on both strands in extension coordinates.
        if query.len() as u32 - read_end >= 1 {
            let r_tail = query.len() as u32 - read_end;
            let ctx = crate::split::SplitContext {
                reference: self.index.reference,
                index: Some(self.index),
                lib,
                read: query,
                contig: out.contig,
                strand: out.strand,
                pos: out.pos,
                ref_end: primary_ref_end,
                read_start,
                read_end,
                cigar: &out.cigar,
                seed_params: self.config.seed,
                splice_params: params,
                extend_params: self.config.extend,
            };
            if let Some(rescue) = crate::split::rescue_right_tail(&ctx, &mut self.buf) {
                if cigar_read_span(&rescue.cigar) == query.len() as u32 {
                    // tail-rescue gate (frozen on): tails <12bp require either
                    // known_support ≥ 2 OR (canonical signal + zero tail-vs-ref
                    // mismatch) — a de-novo short tail with a perfect anchor-side
                    // junction match is accepted.
                    let t3_reject = r_tail < 12 && {
                        let contig = &self.index.reference.contigs[out.contig as usize];
                        let mm = count_tail_mismatches(&rescue.cigar, out.pos, contig, query, true);
                        if rescue.junction.known_support < 2 {
                            mm > 0 || rescue.junction.signal == SpliceSignal::NonCanonical
                        } else {
                            mm > 0
                        }
                    };
                    if !t3_reject {
                        out.cigar = rescue.cigar;
                        out.junctions.push(rescue.junction);
                        out.score += rescue.score;
                        self.scratch.rc = rc_storage;
                        return out;
                    }
                }
            }
        }
        if read_start >= 1 {
            let ctx = crate::split::SplitContext {
                reference: self.index.reference,
                index: Some(self.index),
                lib,
                read: query,
                contig: out.contig,
                strand: out.strand,
                pos: out.pos,
                ref_end: primary_ref_end,
                read_start,
                read_end,
                cigar: &out.cigar,
                seed_params: self.config.seed,
                splice_params: params,
                extend_params: self.config.extend,
            };
            if let Some(rescue) = crate::split::rescue_left_tail(&ctx, &mut self.buf) {
                if cigar_read_span(&rescue.cigar) == query.len() as u32 {
                    // Mirror of the right-tail gate (frozen asymmetry removed):
                    // de-novo short tails accepted with canonical signal +
                    // zero tail-vs-ref mismatch.
                    let t3_reject = read_start < 12 && {
                        let contig = &self.index.reference.contigs[out.contig as usize];
                        let mm = count_tail_mismatches(
                            &rescue.cigar,
                            rescue.pos,
                            contig,
                            query,
                            false,
                        );
                        if rescue.junction.known_support < 2 {
                            mm > 0 || rescue.junction.signal == SpliceSignal::NonCanonical
                        } else {
                            mm > 0
                        }
                    };
                    if !t3_reject {
                        out.pos = rescue.pos;
                        out.cigar = rescue.cigar;
                        out.junctions.push(rescue.junction);
                        out.score += rescue.score;
                    }
                }
            }
        }
        self.scratch.rc = rc_storage;
        out
    }

    /// Library-driven run-through reinterpretation (span.rs; skipped without
    /// a library or when the alignment is already spliced).
    fn apply_span_rescue(&mut self, seq: &[u8], mut aln: ReadAlignment) -> ReadAlignment {
        let rc_storage = self.take_rc(seq);
        let Some(lib) = self.jlib.as_deref() else {
            self.scratch.rc = rc_storage;
            return aln;
        };
        if !aln.junctions.is_empty() {
            self.scratch.rc = rc_storage;
            return aln; // already has a spliced interpretation
        };
        let query: &[u8] = match aln.strand {
            Strand::Plus => seq,
            Strand::Minus => &rc_storage,
        };
        let ctx = crate::split::SplitContext {
            reference: self.index.reference,
            index: Some(self.index),
            lib,
            read: query,
            contig: aln.contig,
            strand: aln.strand,
            pos: aln.pos,
            ref_end: 0, // span recomputes ref_end from pos+cigar
            read_start: 0,
            read_end: 0,
            cigar: &aln.cigar,
            seed_params: self.config.seed,
            splice_params: SpliceParams::default(),
            extend_params: self.config.extend,
        };
        if let Some(rescue) = crate::span::rescue_span(&ctx, &mut self.buf) {
            if cigar_read_span(&rescue.cigar) == query.len() as u32 {
                aln.cigar = rescue.cigar;
                aln.junctions.push(rescue.junction);
                aln.score += rescue.score;
            }
        }
        self.scratch.rc = rc_storage;
        aln
    }

    /// Annotated micro-deletion → junction: a D of ≤50bp exactly matching a
    /// library junction is rewritten to N (1-40bp annotated introns and
    /// deletions are sequence-indistinguishable; only the library tells).
    fn convert_annotated_micro_deletions(&self, mut aln: ReadAlignment) -> ReadAlignment {
        let Some(lib) = self.jlib.as_deref() else {
            return aln;
        };
        let mut pos = aln.pos;
        let mut new_cigar = Vec::with_capacity(aln.cigar.len());
        for op in &aln.cigar {
            let adv = match op {
                CigarOp::Match(n) | CigarOp::Del(n) | CigarOp::RefSkip(n) => *n,
                _ => 0,
            };
            let converted = match op {
                CigarOp::Del(n) if *n <= 50 => {
                    let j = Junction {
                        contig: aln.contig,
                        start: pos,
                        end: pos + *n,
                        minus_strand: aln.strand == Strand::Minus,
                    };
                    if lib.contains(&j) {
                        aln.junctions.push(RefinedJunction {
                            junction: j,
                            signal: SpliceSignal::GtAg,
                            known_support: lib.support(&j),
                        });
                        CigarOp::RefSkip(*n)
                    } else {
                        *op
                    }
                }
                _ => *op,
            };
            new_cigar.push(converted);
            pos += adv;
        }
        aln.cigar = new_cigar;
        aln
    }

    /// No-DP fast lane for decisive, single-diagonal, fully-covered chains.
    ///
    /// Equivalence (zero-mismatch gate + Kadane max subarray): gap-free affine
    /// DP on one diagonal = Kadane; with zero mismatches the score is the
    /// theoretical maximum, so DP can improve neither via gaps nor off-diagonal.
    /// Returns (Extension, wstart) matching extend_hint's coordinates.
    fn try_fast_lane(
        &mut self,
        seq: &[u8],
        rc_storage: &[u8],
        best: &Chain,
        read_len: u32,
        k: u32,
    ) -> Option<(Extension, u32)> {
        let cfg = self.config;

        // Condition 2: all anchors on the same diagonal (ref_pos − read_pos).
        let first = best.anchors.first()?;
        let diag = first.rpos as i64 - first.qpos as i64;
        if !best
            .anchors
            .iter()
            .all(|a| a.rpos as i64 - a.qpos as i64 == diag)
        {
            return None;
        }

        // Condition 3: adjacent anchor read-gap ≤ k+w+5 (=25bp).
        let max_gap = k + cfg.seed.w + 5;
        for w in best.anchors.windows(2) {
            if w[1].qpos - w[0].qpos > max_gap {
                return None;
            }
        }

        // Condition 4: chain read span covers the full read (≤2bp overhangs).
        let (qstart, qend) = best.read_span(k);
        if qstart > 2 || read_len - qend > 2 {
            return None;
        }

        // rc_storage is the caller's revcomp(seq) buffer (same bytes).
        let query: &[u8] = match best.strand {
            Strand::Plus => seq,
            Strand::Minus => rc_storage,
        };

        // Window formula identical to the DP extend path.
        let contig = &self.index.reference.contigs[best.contig as usize];
        let (rs, _re) = best.ref_span(k);
        let flank = read_len / 2 + 30;
        let wstart = rs.saturating_sub(flank);
        let ep = cfg.extend;

        // Kadane scan along the diagonal. read_pos i ↔ contig_pos (i + diag).
        // Per-position scores match SubstLut: match=+2, mismatch=−4, EA=0, N=0.
        let mut best_score = 0i32;
        let mut best_start = 0usize;
        let mut best_end = 0usize;
        let mut cur_score = 0i32;
        let mut cur_start = 0usize;
        let mut cur_mm = 0u32;
        let mut best_mm = 0u32;

        for (i, &qb) in query.iter().enumerate() {
            let rpos = i as i64 + diag;
            let (s, is_mm) = if rpos < 0 || rpos >= contig.len as i64 {
                (0, false)
            } else {
                let rb = contig.base(rpos as u32).to_ascii();
                if qb == b'N' || rb == b'N' {
                    (0, false)
                } else if qb == rb {
                    (ep.match_score, false)
                } else if ep.edit_dir.exempt(rb, qb) {
                    (0, false)
                } else {
                    (ep.mismatch, true)
                }
            };
            // Kadane with ≤0 reset (matches DP H[i] = max(0, H[i−1] + s(i))).
            if cur_score <= 0 {
                cur_score = s;
                cur_start = i;
                cur_mm = if is_mm { 1 } else { 0 };
            } else {
                cur_score += s;
                if is_mm {
                    cur_mm += 1;
                }
            }
            if cur_score > best_score {
                best_score = cur_score;
                best_start = cur_start;
                best_end = i + 1;
                best_mm = cur_mm;
            }
        }

        if best_score <= 0 {
            return None;
        }

        // Condition 5 (zero-mismatch gate): any mismatch in the best subarray
        // opens a gap-improvement opportunity for DP → refuse the fast lane.
        if best_mm > 0 {
            return None;
        }

        let read_start = best_start as u32;
        let read_end = best_end as u32;
        let match_len = read_end - read_start;

        let mut cigar = Vec::with_capacity(3);
        if read_start > 0 {
            cigar.push(CigarOp::SoftClip(read_start));
        }
        cigar.push(CigarOp::Match(match_len));
        let trailing = read_len - read_end;
        if trailing > 0 {
            cigar.push(CigarOp::SoftClip(trailing));
        }

        let ref_start = (read_start as i64 + diag - wstart as i64) as u32;
        let ext = Extension {
            score: best_score,
            ref_start,
            read_start,
            read_end,
            cigar,
        };

        Some((ext, wstart))
    }

    /// mm/ea recount along the CIGAR, base-by-base (extension coordinates).
    /// ea = (ref A, read G)/(ref T, read C) editing-type mismatches; mm = rest.
    pub(crate) fn recount_mm_ea(&mut self, query: &[u8], aln: &mut ReadAlignment) {
        recount_mm_ea_in(self.index.reference, query, aln);
    }

    pub(crate) fn chain_round(
        &mut self,
        mins: &[Minimizer],
        read_len: u32,
        k: u32,
        occ_cap: u32,
    ) -> Option<RoundOutcome> {
        let mut anchors = std::mem::take(&mut self.scratch.anchors);
        collect_anchors_into(self.index, mins, read_len, k, occ_cap, &mut anchors);
        if anchors.is_empty() {
            self.scratch.anchors = anchors;
            return None;
        }
        let mut chains = std::mem::take(&mut self.scratch.chains_out);
        chain_anchors_into(
            &mut anchors,
            &self.config.chain,
            &mut chains,
            &mut self.scratch.chain,
        );
        self.scratch.anchors = anchors;
        let second = second_score_with(&chains, &mut self.scratch.second_keys);
        Some(RoundOutcome { chains, second })
    }
}


/// Set `aln.mm_count` / `aln.ea_count` by walking the CIGAR base-by-base
/// against `reference`. `query` must already be in the alignment's frame
/// (reverse-complemented for [`Strand::Minus`]). Callers that synthesise an
/// alignment without an extension pass through here so the counts — and the
/// NM tag derived from them — describe the placement that was actually
/// written.
pub(crate) fn recount_mm_ea_in(
    reference: &crate::fasta::Reference,
    query: &[u8],
    aln: &mut ReadAlignment,
) {
    let contig = &reference.contigs[aln.contig as usize];
    let mut scratch: Vec<u8> = Vec::new();
    let (mut ro, mut rf) = (0u32, aln.pos);
    let (mut mm, mut ag, mut tc) = (0u32, 0u32, 0u32);
    for op in &aln.cigar {
        match op {
            CigarOp::Match(n) => {
                // Defensive: L1 projection can place a read past a
                // contig end (bad annotation); clamp instead of panic.
                let avail = (*n).min(contig.len.saturating_sub(rf));
                let s0 = scratch.len();
                if avail > 0 {
                    contig.decode_append(rf, rf + avail, &mut scratch);
                }
                let got = scratch.len() - s0;
                for k in 0..(*n as usize) {
                    let qb = query[ro as usize + k];
                    let rb = if k < got { scratch[s0 + k] } else { b'N' };
                    if qb != b'N' && rb != b'N' && qb != rb {
                        match crate::extend::exempt_dir(rb, qb) {
                            Some(true) => ag += 1,
                            Some(false) => tc += 1,
                            None => mm += 1,
                        }
                    }
                }
                ro += *n;
                rf += *n;
            }
            CigarOp::Ins(n) | CigarOp::SoftClip(n) => ro += *n,
            CigarOp::Del(n) | CigarOp::RefSkip(n) => rf += *n,
        }
    }
    aln.mm_count = mm;
    aln.ea_ag = ag;
    aln.ea_tc = tc;
}

/** Reverse complement (N stays N).
 */
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    revcomp_into(seq, &mut out);
    out
}

/// [`revcomp`] into a reusable buffer (`out` is cleared, then filled).
pub fn revcomp_into(seq: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.extend(
        seq.iter()
            .rev()
            .map(|&b| crate::fasta::Base::from_ascii(b).complement().to_ascii()),
    );
}

/// Total read consumption (M/I/S) of a CIGAR.
fn cigar_read_span(cigar: &[CigarOp]) -> u32 {
    cigar
        .iter()
        .map(|op| match op {
            CigarOp::Match(n) | CigarOp::Ins(n) | CigarOp::SoftClip(n) => *n,
            _ => 0,
        })
        .sum()
}

/// Tail-rescue: count base mismatches (+ indels) in the tail portion of a rescue CIGAR.
/// `right_tail = true` → tail after RefSkip; `false` → tail before RefSkip.
/// Non-zero = imperfect tail (any mismatch, insertion, or deletion).
fn count_tail_mismatches(
    cigar: &[CigarOp],
    pos: u32,
    contig: &crate::fasta::Contig,
    query: &[u8],
    right_tail: bool,
) -> u32 {
    let mut ref_pos = pos;
    let mut read_pos = 0u32;
    let mut seen_skip = false;
    let mut mismatches = 0u32;
    for op in cigar {
        let in_tail = if right_tail { seen_skip } else { !seen_skip };
        match op {
            CigarOp::RefSkip(n) => {
                seen_skip = true;
                ref_pos += n;
            }
            CigarOp::Match(n) => {
                if in_tail {
                    for i in 0..*n {
                        let qb = query[(read_pos + i) as usize];
                        // Out-of-bounds ⇒ mismatch (old base() None semantics;
                        // new base() is total and would panic OOB).
                        if ref_pos + i >= contig.len {
                            mismatches += 1;
                            continue;
                        }
                        let ra = contig.base(ref_pos + i).to_ascii();
                        if ra != b'N' && !qb.eq_ignore_ascii_case(&ra) {
                            mismatches += 1;
                        }
                    }
                }
                ref_pos += n;
                read_pos += n;
            }
            CigarOp::Ins(n) => {
                if in_tail {
                    mismatches += n;
                }
                read_pos += n;
            }
            CigarOp::Del(n) => {
                if in_tail {
                    mismatches += n;
                }
                ref_pos += n;
            }
            CigarOp::SoftClip(n) => {
                read_pos += n;
            }
        }
    }
    mismatches
}

/// Decoy contract: contigs whose name carries the `decoy_` prefix are
/// decoys. The prefix travels in the FASTA itself, so the index, the BAM
/// header and every downstream stage see it without a sidecar.
pub fn is_decoy_contig(reference: &crate::fasta::Reference, contig: u32) -> bool {
    reference
        .contigs
        .get(contig as usize)
        .is_some_and(|c| c.name.starts_with("decoy_"))
}

/// Public wrapper of [`ref_m_blocks`] for sibling crates (verify's column
/// mapping tests).
pub fn ref_m_blocks_public(pos: u32, cigar: &[CigarOp]) -> Vec<(u32, u32)> {
    ref_m_blocks(pos, cigar)
}

/// Reference-side M blocks (spliced gaps excluded) of a placement, as
/// sorted disjoint `(start, end)` intervals.
fn ref_m_blocks(pos: u32, cigar: &[CigarOp]) -> Vec<(u32, u32)> {
    let mut out: Vec<(u32, u32)> = Vec::new();
    let mut rf = pos;
    for op in cigar {
        match op {
            CigarOp::Match(n) => {
                out.push((rf, rf + *n));
                rf += *n;
            }
            CigarOp::Del(n) | CigarOp::RefSkip(n) => rf += *n,
            _ => {}
        }
    }
    out
}

/// Footprint overlap of two placements: shared reference M bases divided by
/// the smaller footprint. Two candidates on the same contig and strand whose
/// overlap reaches 0.5 are the same locus; anything below is a genuine
/// alternative placement.
fn footprint_overlap(a: &Placement, b: &Placement) -> f64 {
    if a.contig != b.contig || a.strand != b.strand {
        return 0.0;
    }
    let ba = ref_m_blocks(a.pos, &a.cigar);
    let bb = ref_m_blocks(b.pos, &b.cigar);
    let la: u32 = ba.iter().map(|(s, e)| e - s).sum();
    let lb: u32 = bb.iter().map(|(s, e)| e - s).sum();
    if la == 0 || lb == 0 {
        return 0.0;
    }
    let mut inter: u32 = 0;
    let (mut i, mut j) = (0usize, 0usize);
    while i < ba.len() && j < bb.len() {
        let lo = ba[i].0.max(bb[j].0);
        let hi = ba[i].1.min(bb[j].1);
        if hi > lo {
            inter += hi - lo;
        }
        if ba[i].1 < bb[j].1 {
            i += 1;
        } else {
            j += 1;
        }
    }
    f64::from(inter) / f64::from(la.min(lb))
}

/// Decision order: score desc, then fewer ordinary mismatches, then the
/// genomic layer, then coordinates (deterministic).
fn sort_placements(ps: &mut [Placement]) {
    ps.sort_by(|a, b| {
        b.strict
            .score
            .cmp(&a.strict.score)
            .then(a.strict.mm.cmp(&b.strict.mm))
            .then(a.origin.rank().cmp(&b.origin.rank()))
            .then(a.contig.cmp(&b.contig))
            .then(a.pos.cmp(&b.pos))
    });
}

/// The editing direction a placement's mismatch spectrum implies when it
/// looks hyperedited: single-direction edits at or above
/// [`crate::mapq::HYPEREDIT_MIN_DENSITY`] of the aligned columns with at most
/// [`crate::mapq::HYPEREDIT_MAX_NON_EDIT_MM`] ordinary mismatches (the same
/// criterion as `ReadAlignment::hyperedited_signature`). None otherwise.
fn hyper_direction(p: &Placement) -> Option<EditDir> {
    let s = &p.strict;
    if s.mm > crate::mapq::HYPEREDIT_MAX_NON_EDIT_MM || (s.ea_ag > 0 && s.ea_tc > 0) {
        return None;
    }
    let ea = s.ea_ag + s.ea_tc;
    if s.aligned == 0 || f64::from(ea) / f64::from(s.aligned) < crate::mapq::HYPEREDIT_MIN_DENSITY {
        return None;
    }
    Some(EditDir::from_plus(s.ea_ag > 0))
}

fn dedup_by_footprint(ps: &mut Vec<Placement>) {
    sort_placements(ps);
    let mut kept: Vec<Placement> = Vec::with_capacity(ps.len());
    for p in ps.drain(..) {
        if kept.iter().all(|k| footprint_overlap(k, &p) < 0.5) {
            kept.push(p);
        }
    }
    *ps = kept;
}

impl ReadAlignment {
    #[cfg(test)]
    fn mapq_of_checked(&self) -> u8 {
        crate::mapq::mapq_of(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ref_from(seq: &str) -> crate::fasta::Reference {
        crate::fasta::parse_fasta_bytes(format!(">chrT\n{seq}\n").as_bytes()).expect("fasta")
    }

    fn aligner_on(seq: &str) -> Aligner<'static> {
        let reference = ref_from(seq);
        let index = Box::leak(Box::new(Index::build(
            reference,
            SeedParams { k: 15, w: 5 },
        )));
        Aligner::new(index, AlignConfig::rna_default())
    }

    fn placement(contig: u32, pos: u32, cigar: &[CigarOp]) -> Placement {
        Placement::from_genomic(
            contig,
            pos,
            Strand::Plus,
            cigar.to_vec(),
            Vec::new(),
            0,
            100,
        )
    }

    #[test]
    fn reused_buffers_two_reads_independent() {
        // One Aligner (one set of reusable buffers), consecutive reads with
        // different content: no state may leak between calls, and a replay
        // of the first read must reproduce the first result exactly
        // (allocation-elimination invariant).
        let contig = "GATTACAGGCATTACGACCAAGTTACAGGATCCGATTAAGGCTTACGGATCCTAGCATCGGATCCATGGCATTGCAAGTGGCATTGACCAAGTTTGCAAGTTACAGGATCCGATTAAGGCTTACGGATCCTAGCATCGGATCCATGGCATT";
        let mut aln = aligner_on(contig);
        let r1 = &contig.as_bytes()[0..60];
        let r2 = &contig.as_bytes()[91..151];

        let a1 = aln.align_read(r1).expect("first read placed");
        assert_eq!(a1.pos, 0);

        let mut r2n = r2.to_vec();
        r2n[30] = b'N'; // different content through the same scratch buffers
        let a2 = aln.align_read(&r2n).expect("second read placed");
        assert_eq!(a2.pos, 91);
        assert_ne!(a1.pos, a2.pos);

        // Replay of the first read: identical output to the first call.
        let a1b = aln.align_read(r1).expect("replay placed");
        let mq1b = a1b.mapq_of_checked();
        let mq1 = a1.mapq_of_checked();
        assert_eq!(
            (a1b.pos, a1b.strand, a1b.score, a1b.cigar),
            (a1.pos, a1.strand, a1.score, a1.cigar)
        );
        assert_eq!(mq1b, mq1);
    }

    #[test]
    fn rescore_strict_manual_cigars() {
        let mut aln = aligner_on(&"ACGT".repeat(50));
        let read = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
        // 50M perfect at pos 0: 50 matches.
        let mut p = placement(0, 0, &[CigarOp::Match(50)]);
        aln.rescore_strict(read, &mut p);
        assert_eq!(p.strict.score, 100);
        assert_eq!(p.strict.matches, 50);
        assert_eq!(p.strict.mm, 0);
        assert_eq!(p.strict.aligned, 50);
        // One ordinary mismatch: 49*2 - 4.
        let mut r2 = read.to_vec();
        r2[25] = if r2[25] == b'A' { b'C' } else { b'A' };
        let mut p = placement(0, 0, &[CigarOp::Match(50)]);
        aln.rescore_strict(&r2, &mut p);
        assert_eq!(p.strict.score, 49 * 2 - 4);
        assert_eq!(p.strict.mm, 1);
        // An A>G editing-type substitution is charged identically but
        // counted in the editing channel.
        let mut r3 = read.to_vec();
        let ref_b = aln.index.reference.contigs[0].base(25).to_ascii();
        if ref_b == b'A' {
            r3[25] = b'G'; // ref A / read G
        } else {
            // shift to a position with ref A: "ACGT"*n has A at multiples of 4
            let pos = 24; // A at 24 (0-based) in "ACGT" repeats
            let _ = pos;
            r3[24] = b'G';
        }
        let mut p = placement(0, 0, &[CigarOp::Match(50)]);
        aln.rescore_strict(&r3, &mut p);
        assert_eq!(p.strict.score, 49 * 2 - 4);
        assert_eq!(p.strict.mm, 0);
        assert!(p.strict.ea_ag + p.strict.ea_tc == 1);
        // Insertion of 2: 48 aligned matches, -(4 + 2*2).
        let mut ins_read: Vec<u8> = read[..24].to_vec();
        ins_read.extend_from_slice(b"TT");
        ins_read.extend_from_slice(&read[24..48]);
        let mut p = placement(0, 0, &[CigarOp::Match(24), CigarOp::Ins(2), CigarOp::Match(24)]);
        aln.rescore_strict(&ins_read, &mut p);
        assert_eq!(p.strict.score, 48 * 2 - (4 + 2 * 2));
        assert_eq!(p.strict.ins, 1);
    }

    #[test]
    fn rescore_strict_splice_signal() {
        // Exon1 | intron(GT + 56N + AG) | exon2; read spans 20+20 exon bases.
        let exon1 = "ACGTTGCAAGGCATTACGACCAAGTTACAGGATCCAAGTTGACGTTACAGGCATTAC"; // 57
        let intron = format!("GT{}AG", "N".repeat(56)); // 60 total, GT..AG inside
        let exon2 = "CATGGACTTGAACCTTGACGTACATGGACTTGAACCTTGACGTA"; // 44
        let seq = format!("{exon1}{intron}{exon2}");
        let mut aln = aligner_on(&seq);
        let start = (exon1.len() - 20) as u32; // last 20 exon1 bases
        let read: Vec<u8> = {
            let s = seq.bytes().collect::<Vec<u8>>();
            let mut r: Vec<u8> = s[(start as usize)..exon1.len()].to_vec();
            let e2 = exon1.len() + intron.len();
            r.extend_from_slice(&s[e2..e2 + 20]);
            r
        };
        let n = intron.len() as u32;
        let mut p = placement(0, start, &[CigarOp::Match(20), CigarOp::RefSkip(n), CigarOp::Match(20)]);
        aln.rescore_strict(&read, &mut p);
        let expected = 40 * 2 + SpliceParams::default().gt_ag_bonus
            - intron_length_penalty(n);
        assert_eq!(p.strict.score, expected);
        assert_eq!(p.strict.n_splice, 1);
        assert_eq!(p.strict.aligned, 40);
        // A read antisense to the transcript crosses the same GT..AG intron
        // (seen as CT..AC from its side); the strict score must not depend
        // on which strand the read was sequenced from.
        let read_rc = revcomp(&read);
        let mut pm = placement(0, start, &[CigarOp::Match(20), CigarOp::RefSkip(n), CigarOp::Match(20)]);
        pm.strand = Strand::Minus;
        aln.rescore_strict(&revcomp(&read_rc), &mut pm);
        assert_eq!(pm.strict.score, expected);
    }

    #[test]
    fn exact_kmer_support_counts() {
        let mut aln = aligner_on(&"ACGT".repeat(50));
        let read: Vec<u8> = b"ACGT".repeat(13)[..50].to_vec();
        let p = placement(0, 0, &[CigarOp::Match(50)]);
        assert_eq!(aln.exact_kmer_support(&read, &p), 3); // 0,15,30
        // One mismatch in the second 15-mer window pushes the scan past it.
        let mut r2 = read.clone();
        r2[20] = if r2[20] == b'A' { b'C' } else { b'A' };
        let n = aln.exact_kmer_support(&r2, &p);
        assert!(n < 3);
    }

    #[test]
    fn footprint_overlap_boundaries() {
        let a = placement(0, 0, &[CigarOp::Match(100)]);
        let half = placement(0, 50, &[CigarOp::Match(100)]);
        assert!(footprint_overlap(&a, &half) >= 0.5); // exactly 0.5: same locus
        let below = placement(0, 51, &[CigarOp::Match(100)]);
        assert!(footprint_overlap(&a, &below) < 0.5);
        let disjoint = placement(0, 200, &[CigarOp::Match(100)]);
        assert_eq!(footprint_overlap(&a, &disjoint), 0.0);
        // Different strand never overlaps.
        let mut other_strand = placement(0, 0, &[CigarOp::Match(100)]);
        other_strand.strand = Strand::Minus;
        assert_eq!(footprint_overlap(&a, &other_strand), 0.0);
        // Spliced placements overlap on M blocks only.
        let spliced = placement(0, 0, &[CigarOp::Match(50), CigarOp::RefSkip(1000), CigarOp::Match(50)]);
        let head = placement(0, 0, &[CigarOp::Match(50)]);
        assert!(footprint_overlap(&spliced, &head) > 0.9);
    }

    #[test]
    fn dedup_keeps_best_strict_per_locus() {
        let mut a = placement(0, 0, &[CigarOp::Match(100)]);
        a.strict.score = 150;
        let mut b = placement(0, 10, &[CigarOp::Match(100)]);
        b.strict.score = 200; // same locus, better
        let mut c = placement(0, 5000, &[CigarOp::Match(100)]);
        c.strict.score = 100; // different locus
        let mut ps = vec![a, b.clone(), c];
        dedup_by_footprint(&mut ps);
        assert_eq!(ps.len(), 2);
        assert_eq!(ps[0].strict.score, 200);
        assert_eq!(ps[0].pos, b.pos);
        assert_eq!(ps[1].strict.score, 100);
    }

    #[test]
    fn refuse_placement_below_two_thirds_coverage() {
        // 300bp contig; read = 75bp of the contig + 75bp absent sequence.
        let contig = "GATTACAGGCATTACGACCAAGTTACAGGATCCGATTAAGGCTTACGGATCCTAGCATCGGATCCATGGCATTGCAAGTGGCATTGACCAAGTT"
            .to_string()
            + "TGCAAGTTACAGGATCCGATTAAGGCTTACGGATCCTAGCATCGGATCCATGGCATT";
        let mut aln = aligner_on(&contig);
        let head = &contig.as_bytes()[..75];
        let tail = b"T".repeat(75).leak();
        let mut read = head.to_vec();
        read.extend_from_slice(&tail[..]);
        assert_eq!(read.len(), 150);
        assert!(aln.align_read(&read).is_none());
    }

    #[test]
    fn repeat_copy_scores_mapq_zero() {
        // Two identical 120bp copies: a 60bp read from the copy seeds both.
        let copy = "ACGTTGCAAGGCATTACGACCAAGTTACAGGATCCGATTAAGGCTTACGGATCCTAGCATCGGATCCATGGCATTGACCAAGTTACAGGATCCGATTAAGGC";
        let spacer = "A".repeat(120);
        let seq = format!("{copy}{spacer}{copy}");
        let mut aln = aligner_on(&seq);
        let read = &copy.as_bytes()[30..90];
        let aln_out = aln.align_read(read).expect("placed");
        assert_eq!(aln_out.mapq_of_checked(), 0);
    }

    /// Exon1 | GT..AG intron(60) | exon2 fixture with a 40bp read spanning
    /// 20+20 exon bases and one 15-mer anchor per exon (read-gap 5,
    /// ref-gap 65 ⇒ implied intron exactly 60).
    fn spliced_chain_fixture() -> (String, Vec<u8>, Chain) {
        let exon1 = "ACGTTGCAAGGCATTACGACCAAGTTACAGGATCCAAGTTGACGTTACAGGCATTAC"; // 57
        let intron = format!("GT{}AG", "N".repeat(56)); // 60, GT..AG inside
        let exon2 = "CATGGACTTGAACCTTGACGTACATGGACTTGAACCTTGACGTA"; // 44
        let seq = format!("{exon1}{intron}{exon2}");
        let s = seq.as_bytes();
        let start = (exon1.len() - 20) as u32; // last 20 exon1 bases
        let e2 = exon1.len() + intron.len();
        let mut read: Vec<u8> = s[start as usize..exon1.len()].to_vec();
        read.extend_from_slice(&s[e2..e2 + 20]);
        let chain = Chain {
            contig: 0,
            strand: Strand::Plus,
            anchors: vec![
                crate::seed::Anchor {
                    qpos: 0,
                    rpos: start,
                    contig: 0,
                    strand: Strand::Plus,
                },
                crate::seed::Anchor {
                    qpos: 20,
                    rpos: e2 as u32,
                    contig: 0,
                    strand: Strand::Plus,
                },
            ],
            score: 200,
        };
        (seq, read, chain)
    }

    #[test]
    fn try_spliced_stitches_two_segment_chain() {
        let (seq, read, chain) = spliced_chain_fixture();
        let mut aln = aligner_on(&seq);
        let rc = revcomp(&read);
        let hit = aln.try_spliced(&read, &rc, &chain).expect("spliced hit");
        assert_eq!(hit.pos, 37); // exon1 is 57bp: last 20 exon1 bases start at 37
        assert_eq!(
            hit.cigar,
            vec![CigarOp::Match(20), CigarOp::RefSkip(60), CigarOp::Match(20)]
        );
        assert_eq!(hit.junctions.len(), 1);
        assert_eq!(hit.junctions[0].signal, SpliceSignal::GtAg);
        assert_eq!(hit.junctions[0].junction.start, 57);
        assert_eq!(hit.junctions[0].junction.end, 117);
    }

}

