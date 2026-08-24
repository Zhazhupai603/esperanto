//! Per-read alignment driver: seed → anchor → chain → (L1 / intron-chain /
//! spliced / EA-Myers / banded DP) → split/span rescue → recount → coverage gate.
//!
//! All production gates are frozen (no env switches): the splice gate (skip try_spliced on
//! <2 segments), the tail-rescue gate (small-tail rescue tightening), multitie marking, EA-Myers
//! decisive path, window cap. Profiling instrumentation is not ported.

use std::sync::Arc;

use crate::chain::{chain_anchors, second_score, tied_best, Chain, TOP_CANDIDATES};
use crate::extend::{extend_hint, CigarOp, DiagHint, ExtendBuffer, ExtendParams, Extension};
use crate::gtf::{Junction, JunctionLib, RefinedJunction, SpliceSignal};
use crate::index::{collect_anchors, Index};
use crate::mapq::ReadAlignment;
use crate::seed::{minimizers, Minimizer, SeedParams, Strand};
use crate::splice::{align_spliced, refine_junction, segment_chain, SpliceParams};

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

/// EA tie-break winner candidate: (contig, strand, score, pos, n_anchors, ext, mm, ea).
type EaBest = Option<(u32, Strand, i32, u32, usize, Extension, u32, u32)>;

/// Per-thread aligner state (reuses DP buffers).
pub struct Aligner<'a> {
    pub index: &'a Index,
    pub config: AlignConfig,
    /// Junction library (sjdb / 2-pass discoveries; RNA mode).
    pub jlib: Option<Arc<JunctionLib>>,
    /// L1 transcriptome-first engine bundle (None = pure legacy path).
    pub l1: Option<Arc<esperanto_engine::L1Index>>,
    pub(crate) buf: ExtendBuffer,
}

impl<'a> Aligner<'a> {
    pub fn new(index: &'a Index, config: AlignConfig) -> Self {
        Aligner {
            index,
            config,
            jlib: None,
            l1: None,
            buf: ExtendBuffer::default(),
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

    /// Map an L1 engine outcome to a ReadAlignment.
    ///
    /// Conservation guard: the projected CIGAR's query consumption (M+I+S)
    /// must equal read_len, otherwise None (caller falls back to legacy).
    /// chain_score/second_chain_score/n_anchors are reverse-engineered so the
    /// legacy mapq() formula reproduces the engine's MAPQ (no formula change).
    #[allow(clippy::too_many_arguments)]
    fn l1_to_read_alignment(
        &self,
        contig: u32,
        pos: u32,
        strand: esperanto_engine::Strand,
        cigar: Vec<esperanto_engine::CigarOp>,
        score: i32,
        mapq: u8,
        read_len: usize,
    ) -> Option<ReadAlignment> {
        // txmap contig id → .paidx contig id via name lookup.
        let genome_contig = self
            .l1
            .as_ref()
            .and_then(|l1| l1.contig_name(contig))
            .and_then(|name| self.index.reference.contig_index(name.as_bytes()))
            .unwrap_or(contig);

        let map_strand = match strand {
            esperanto_engine::Strand::Plus => Strand::Plus,
            esperanto_engine::Strand::Minus => Strand::Minus,
        };

        let map_cigar: Vec<CigarOp> = cigar
            .iter()
            .map(|op| match op {
                esperanto_engine::CigarOp::Match(n) => CigarOp::Match(*n),
                esperanto_engine::CigarOp::Ins(n) => CigarOp::Ins(*n),
                esperanto_engine::CigarOp::Del(n) => CigarOp::Del(*n),
                esperanto_engine::CigarOp::RefSkip(n) => CigarOp::RefSkip(*n),
                esperanto_engine::CigarOp::SoftClip(n) => CigarOp::SoftClip(*n),
            })
            .collect();

        // SAM invariant guard: query consumption (M+I+S) must equal read_len.
        let qspan: usize = map_cigar
            .iter()
            .map(|op| match op {
                CigarOp::Match(n) | CigarOp::Ins(n) | CigarOp::SoftClip(n) => *n as usize,
                _ => 0,
            })
            .sum();
        if qspan != read_len {
            return None; // non-conserving → caller falls back to legacy
        }

        // mapq = round(60 × (1 − s2/s1) × min(1, n_anchors/10)).
        // With n_anchors=10 and s1=read_len: s2 = round(s1 × (1 − mapq/60)).
        let s1 = read_len.max(1) as i32;
        let n_anchors = 10usize;
        let s2 = if mapq >= 60 {
            0
        } else {
            let frac = 1.0 - (mapq as f64 / 60.0);
            (s1 as f64 * frac).round() as i32
        };
        let s2 = s2.clamp(0, s1);

        Some(ReadAlignment {
            contig: genome_contig,
            pos,
            strand: map_strand,
            score,
            chain_score: s1,
            second_chain_score: s2,
            cigar: map_cigar,
            n_anchors,
            junctions: Vec::new(),
            ea_count: 0,
            mm_count: 0,
            n_seeds: n_anchors,
            rescued: false,
        })
    }

    /// Align one read (ASCII ACGTN). No confident chain ⇒ None.
    pub fn align_read(&mut self, seq: &[u8]) -> Option<ReadAlignment> {
        // L1 transcriptome-first fast path.
        if let Some(l1) = self.l1.clone() {
            let outcome = {
                let mut stats = esperanto_engine::ReadStats::default();
                esperanto_engine::align_read(
                    seq,
                    l1.as_ref(),
                    l1.as_ref(),
                    l1.as_ref(),
                    &esperanto_engine::EngineConfig::default(),
                    &esperanto_engine::NoRepeats,
                    &mut stats,
                )
            };
            if let esperanto_engine::L1Outcome::Aligned {
                contig,
                pos,
                strand,
                cigar,
                score,
                mapq,
            } = outcome
            {
                if let Some(mut aln) =
                    self.l1_to_read_alignment(contig, pos, strand, cigar, score, mapq, seq.len())
                {
                    // L1 direct path needs recount (mm/ea otherwise 0 → EK tag
                    // would misread as a clean read).
                    let rc2;
                    let q: &[u8] = match aln.strand {
                        Strand::Plus => seq,
                        Strand::Minus => {
                            rc2 = revcomp(seq);
                            &rc2
                        }
                    };
                    self.recount_mm_ea(q, &mut aln);
                    return Some(aln);
                }
            }
            // Fallback → legacy path below.
        }

        let cfg = self.config;
        let k = cfg.seed.k;
        let read_len = seq.len() as u32;
        let mins = minimizers(seq, cfg.seed);
        if mins.is_empty() {
            return None;
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

        // Editing-aware variant round: only when chaining fully failed (weak).
        let weak = r1
            .as_ref()
            .and_then(|r| r.best())
            .map(|b| b.score)
            .unwrap_or(0)
            == 0;
        let r3 = if cfg.extend.editing_aware && weak {
            let (mut anchors, _var_hits) = crate::index::collect_anchors_edit_variants(
                self.index,
                &mins,
                read_len,
                k,
                MID_OCC_CAP,
            );
            if std::env::var_os("ESP_PROBE").is_some() {
                eprintln!("[probe] r3 variant_anchors={}", anchors.len());
            }
            if anchors.is_empty() {
                None
            } else {
                let (base_anchors, _) =
                    collect_anchors(self.index, &mins, read_len, k, MID_OCC_CAP);
                anchors.extend(base_anchors);
                let mut cp = cfg.chain;
                cp.min_chain_score = 25;
                let chains = chain_anchors(anchors, &cp);
                let second = second_score(&chains);
                if std::env::var_os("ESP_PROBE").is_some() {
                    eprintln!("[probe] r3 best={:?} n_chains={}", chains.first().map(|c| c.score), chains.len());
                }
                Some(RoundOutcome { chains, second })
            }
        } else {
            None
        };

        // Merge rounds: winner = round with the highest best-chain score
        // (strictly greater replaces; r1 first on ties).
        let rounds: Vec<&RoundOutcome> = r1.iter().chain(r3.iter()).collect();
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
            return None; // chain failure
        };
        let cr = rounds[w];
        let second = cr.second;
        let extra_top: Vec<Chain> = rounds
            .iter()
            .flat_map(|r| r.chains.iter().take(TOP_CANDIDATES).cloned())
            .collect();
        let unedit_hit = r3.as_ref().and_then(|r| r.best()).is_some()
            && cr.best().map(|b| b.score) == r3.as_ref().and_then(|r| r.best()).map(|b| b.score);
        let best = cr.best()?.clone();
        // TEMP parity probe: ESP_PROBE=1 dumps chain/segment/stage traces.
        if std::env::var_os("ESP_PROBE").is_some() {
            eprintln!(
                "[probe] best chain: contig={} strand={:?} score={} anchors={} ref_span={:?} read_span={:?} second={} tied={}",
                self.index.reference.contigs[best.contig as usize].name,
                best.strand,
                best.score,
                best.anchors.len(),
                best.ref_span(k),
                best.read_span(k),
                second,
                tied_best(&cr.chains).len()
            );
        }

        // Intron-chain fast path: RNA chain with exactly 2 segments and an
        // intron-sized gap → verified intra-chain intron algorithm first.
        if cfg.rna {
            let min_intron = SpliceParams::default().min_intron;
            let segs = segment_chain(&best, k, min_intron);
            if std::env::var_os("ESP_PROBE").is_some() {
                eprintln!("[probe] intron-chain segs={} (gate: ==2)", segs.len());
            }
            if segs.len() == 2 {
                let empty_lib = JunctionLib::default();
                let lib = self.jlib.as_deref().unwrap_or(&empty_lib);
                let rc_storage;
                let query: &[u8] = match best.strand {
                    Strand::Plus => seq,
                    Strand::Minus => {
                        rc_storage = revcomp(seq);
                        &rc_storage
                    }
                };
                // Old intron-chain parity: it re-collects base anchors
                // internally (occ_cap = MID_OCC_CAP) and filters to the
                // chain's contig+strand — NOT just the chain's own anchors.
                let (mut ic_anchors, _) =
                    collect_anchors(self.index, &mins, read_len, k, MID_OCC_CAP);
                ic_anchors.retain(|a| a.contig == best.contig && a.strand == best.strand);
                if let Some(hit) = crate::intron_chain::try_intron_chain_placement(
                    self.index.reference,
                    lib,
                    query,
                    &ic_anchors,
                    k,
                    &crate::intron_chain::IntronParams::default(),
                    &cfg.extend,
                    &mut self.buf,
                ) {
                    let mut aln = ReadAlignment {
                        contig: best.contig,
                        pos: hit.genomic_start(),
                        strand: best.strand,
                        score: hit.extension.score,
                        chain_score: read_len as i32,
                        second_chain_score: 0,
                        cigar: hit.cigar,
                        n_anchors: 2,
                        junctions: hit.junctions,
                        n_seeds: mins.len(),
                        ..Default::default()
                    };
                    // Post-process: snap junction to splice signal + library.
                    if !aln.junctions.is_empty() {
                        let raw = aln.junctions[0].junction;
                        let minus = matches!(best.strand, Strand::Minus);
                        let sp = SpliceParams::default();
                        if let Some(refined) = refine_junction(
                            self.index.reference,
                            aln.contig,
                            raw.start,
                            raw.end,
                            minus,
                            lib,
                            &sp,
                        ) {
                            let new_intron =
                                refined.junction.end.saturating_sub(refined.junction.start);
                            for op in &mut aln.cigar {
                                if let CigarOp::RefSkip(n) = op {
                                    *n = new_intron;
                                }
                            }
                            aln.junctions[0] = refined;
                        }
                    }
                    self.recount_mm_ea(query, &mut aln);
                    return Some(aln);
                }
            }
        }

        // RNA: splice path (chain segments ≥ 2; splice gate frozen on).
        if cfg.rna {
            let hit = self.try_spliced(seq, &best);
            if std::env::var_os("ESP_PROBE").is_some() {
                eprintln!("[probe] try_spliced -> {}", hit.as_ref().map(|h| h.score).map_or("None".to_string(), |s| format!("Some(score={s})")));
            }
            if let Some(hit) = hit {
                return Some(hit);
            }
        }

        // Fast path: best chain uncontested (no tie, second < 0.8×best) ⇒
        // extend only best; contested → top-N candidate tie-break.
        let tied = tied_best(&cr.chains);
        let decisive = tied.len() <= 1 && second * 5 < best.score * 4;
        let mut candidates: Vec<Chain> = if decisive {
            Vec::new()
        } else {
            cr.chains.iter().take(TOP_CANDIDATES).cloned().collect()
        };
        if candidates.is_empty() {
            candidates.push(best.clone());
        }
        for c in extra_top {
            if candidates.len() >= MAX_TIED_CANDIDATES {
                break;
            }
            let dup = candidates.iter().any(|x| {
                x.contig == c.contig
                    && x.strand == c.strand
                    && x.ref_span(k).0.abs_diff(c.ref_span(k).0) < 50
            });
            if !dup {
                candidates.push(c);
            }
        }
        // Multimap marking (STAR semantics): near-tied competing reads
        // (second ≥ 0.95×best) skip hard tie-breaking — extend best only.
        let near_tie = !decisive && second * 20 >= best.score * 19;
        if near_tie {
            candidates.truncate(1);
        }

        let rc_storage = revcomp(seq);
        let mut best_hit: Option<(u32, Strand, i32, u32, usize, Extension)> = None;
        let mut second_ext_score: i32 = 0;
        let mut ea_counts: Option<(u32, u32)> = None;
        let mut ea_cigar: Vec<CigarOp> = Vec::new();
        let mut ea_pos: u32 = 0;
        let mut ea_done = false;

        // Fast lane: decisive chains with a single diagonal + full coverage take the
        // Kadane no-DP fast lane (point-equivalent to affine-gap DP).
        let mut t4_hit = false;
        if decisive {
            if let Some((ext, wstart)) = self.try_fast_lane(seq, &best, read_len, k) {
                if ext.score > 0 {
                    let pos = wstart + ext.ref_start;
                    best_hit = Some((
                        best.contig,
                        best.strand,
                        ext.score,
                        pos,
                        best.anchors.len(),
                        ext,
                    ));
                    t4_hit = true;
                }
            }
        }

        // E1: EA-Myers candidate tie-break (≤128 single word, ≤256 two-block).
        if read_len > 0 && read_len <= 256 {
            let ea_w = cfg.extend.editing_aware;
            let m = read_len as usize;
            let match_s = cfg.extend.match_score;
            let mm_s = cfg.extend.mismatch;
            let threshold = ((read_len as i32) / 3).max(10);
            let mut ea_best: EaBest = None;
            let mut ea_second: i32 = 0;
            let mut ea_fallback = false;
            for (cand_idx, cand) in candidates.iter().take(MAX_TIED_CANDIDATES).enumerate() {
                if t4_hit && cand_idx == 0 {
                    continue;
                }
                let contig = &self.index.reference.contigs[cand.contig as usize];
                let flank = read_len / 2 + 30;
                let (rs, re) = cand.ref_span(k);
                let wstart = rs.saturating_sub(flank);
                // Window cap: chain head anchor + (read+flank) — repeat-region
                // chains inflate the span; garbage windows burn DP time.
                let wcap = rs + read_len + flank;
                let wend = (re + flank).min(contig.len).min(wcap);
                let window = contig.slice_ascii(wstart, wend);
                let q: &[u8] = match cand.strand {
                    Strand::Plus => seq,
                    Strand::Minus => &rc_storage,
                };
                if window.len() < m {
                    continue;
                }
                let (d, start) = if m <= 128 {
                    crate::myers_ea::infix_best_start(q, &window)
                } else {
                    crate::myers_ea::long::infix_best_start(q, &window)
                };
                if d > threshold {
                    continue;
                }
                if start + m > window.len() {
                    ea_fallback = true;
                    break;
                }
                // Per-base verification: mm (non-EA mismatches) must exactly
                // explain d, otherwise an indel is indicated. N on either side
                // abandons the candidate (recount's N semantics can't align).
                let mut mm_full: i32 = 0;
                let mut has_n = false;
                for kk in 0..m {
                    let qb = q[kk];
                    let rb = window[start + kk];
                    if qb == b'N' || rb == b'N' {
                        has_n = true;
                        break;
                    }
                    if qb != rb && !((rb == b'A' && qb == b'G') || (rb == b'T' && qb == b'C')) {
                        mm_full += 1;
                    }
                }
                if has_n {
                    continue;
                }
                if mm_full != d {
                    // Indel sign. Multimap reads (near_tie): CIGAR differences
                    // are invisible downstream (MAPQ≈0 filtered) — place all-M.
                    if near_tie {
                        let ext = Extension {
                            score: match_s * m as i32 + mm_s * d,
                            ref_start: start as u32,
                            read_start: 0,
                            read_end: m as u32,
                            cigar: vec![CigarOp::Match(read_len)],
                        };
                        let pos = wstart + ext.ref_start;
                        let better = match &ea_best {
                            None => true,
                            Some((_, _, bscore, _, _, _, _, _)) => ext.score > *bscore,
                        };
                        if better {
                            if let Some((_, _, old, _, _, _, _, _)) = &ea_best {
                                ea_second = ea_second.max(*old);
                            }
                            // mm/ea unknown (indel unresolved) — u32::MAX marks
                            // "force recount".
                            ea_best = Some((
                                cand.contig,
                                cand.strand,
                                ext.score,
                                pos,
                                cand.anchors.len(),
                                ext,
                                u32::MAX,
                                0,
                            ));
                        } else {
                            ea_second = ea_second.max(ext.score);
                        }
                        continue;
                    }
                    // Unique-position read: fall back to DP only if this
                    // candidate could still beat the current best.
                    let optimistic = match_s * m as i32 + mm_s * d;
                    let cur = ea_best.as_ref().map(|b| b.2).unwrap_or(i32::MIN);
                    if optimistic > cur {
                        ea_fallback = true;
                        break;
                    }
                    continue;
                }
                // No indel: optimal soft-clip ≡ Kadane max subarray on the
                // diagonal (per-base affine scores).
                let mut best_sum = 0i32;
                let (mut bi, mut bj) = (0usize, 0usize);
                let mut cur_sum = 0i32;
                let mut cur_i = 0usize;
                for kk in 0..m {
                    let qb = q[kk];
                    let rb = window[start + kk];
                    let s = if qb == rb {
                        match_s
                    } else if (rb == b'A' && qb == b'G') || (rb == b'T' && qb == b'C') {
                        if ea_w {
                            0
                        } else {
                            mm_s
                        }
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
                if best_sum <= 0 {
                    continue;
                }
                let seg_len = (bj - bi) as u32;
                if seg_len * 2 < read_len {
                    continue; // DP's acceptance gate wouldn't take it either
                }
                let mut mm2 = 0u32;
                let mut ea2 = 0u32;
                for kk in bi..bj {
                    let qb = q[kk];
                    let rb = window[start + kk];
                    if qb != rb {
                        if (rb == b'A' && qb == b'G') || (rb == b'T' && qb == b'C') {
                            ea2 += 1;
                        } else {
                            mm2 += 1;
                        }
                    }
                }
                let mut cigar = Vec::with_capacity(3);
                if bi > 0 {
                    cigar.push(CigarOp::SoftClip(bi as u32));
                }
                cigar.push(CigarOp::Match(seg_len));
                if m - bj > 0 {
                    cigar.push(CigarOp::SoftClip((m - bj) as u32));
                }
                let ext = Extension {
                    score: best_sum,
                    ref_start: (start + bi) as u32,
                    read_start: bi as u32,
                    read_end: bj as u32,
                    cigar,
                };
                let pos = wstart + ext.ref_start;
                let better = match &ea_best {
                    None => true,
                    Some((_, _, bscore, _, _, _, _, _)) => ext.score > *bscore,
                };
                if better {
                    if let Some((_, _, old, _, _, _, _, _)) = &ea_best {
                        ea_second = ea_second.max(*old);
                    }
                    ea_best = Some((
                        cand.contig,
                        cand.strand,
                        ext.score,
                        pos,
                        cand.anchors.len(),
                        ext,
                        mm2,
                        ea2,
                    ));
                } else {
                    ea_second = ea_second.max(ext.score);
                }
            }
            if !ea_fallback {
                if let Some((c, st, sc, pos, na, ext, mm2, ea2)) = ea_best {
                    ea_cigar = ext.cigar.clone();
                    ea_pos = pos;
                    // mm2==u32::MAX = near-tie multimap indel placement (not
                    // per-base counted) — mark None, force recount.
                    ea_counts = if mm2 == u32::MAX {
                        None
                    } else {
                        Some((mm2, ea2))
                    };
                    second_ext_score = second_ext_score.max(ea_second);
                    best_hit = Some((c, st, sc, pos, na, ext));
                    ea_done = true;
                }
            }
        }

        // Extension tie-break budget: top-2 by banded DP (multimap reads are
        // MAPQ≈0 anyway; exhaustive tie-break has no downstream yield).
        if !ea_done {
            for (cand_idx, cand) in candidates.iter().take(EXTEND_BUDGET).enumerate() {
                // Fast lane: skip best candidate if the fast lane handled it.
                if t4_hit && cand_idx == 0 {
                    continue;
                }
                let contig = &self.index.reference.contigs[cand.contig as usize];
                let flank = read_len / 2 + 30;
                let (rs, re) = cand.ref_span(k);
                let wstart = rs.saturating_sub(flank);
                let wcap = rs + read_len + flank;
                let wend = (re + flank).min(contig.len).min(wcap);
                let window = contig.slice_ascii(wstart, wend);
                let q: &[u8] = match cand.strand {
                    Strand::Plus => seq,
                    Strand::Minus => &rc_storage,
                };
                // Chain diagonal: read anchor q ↔ window j = (rs − wstart) + (i − qs).
                let (qs, _) = cand.read_span(k);
                let hint = DiagHint {
                    offset: (rs - wstart) as i64 - qs as i64,
                    num: 1,
                    den: 1,
                };
                let ext = extend_hint(q, &window, &cfg.extend, &mut self.buf, hint);
                if ext.score <= 0 {
                    continue;
                }
                let pos = wstart + ext.ref_start;
                // Strictly-better replaces; ties keep the earlier (higher
                // chain score — chain score is the repeat-region prior).
                let better = match &best_hit {
                    None => true,
                    Some((_, _, bscore, _, _, _)) => ext.score > *bscore,
                };
                if better {
                    if let Some((_, _, old_bscore, _, _, _)) = &best_hit {
                        second_ext_score = second_ext_score.max(*old_bscore);
                    }
                    best_hit = Some((
                        cand.contig,
                        cand.strand,
                        ext.score,
                        pos,
                        cand.anchors.len(),
                        ext,
                    ));
                } else {
                    second_ext_score = second_ext_score.max(ext.score);
                }
            }
        }

        let (contig_id, strand, score, pos, n_anchors, ext) = best_hit?;

        if cfg.rna {
            let aln = ReadAlignment {
                contig: contig_id,
                pos,
                strand,
                score,
                chain_score: best.score,
                second_chain_score: if second_ext_score > 0 {
                    second_ext_score
                } else {
                    second
                },
                cigar: ext.cigar.clone(),
                n_anchors,
                n_seeds: mins.len(),
                rescued: unedit_hit,
                ..Default::default()
            };
            // Split-DP draft layer skipped (frozen); tail split rescue kept.
            let aln = self.apply_split_rescue(seq, aln);
            let aln = self.apply_span_rescue(seq, aln);
            let mut aln = self.convert_annotated_micro_deletions(aln);
            let rc2;
            let q: &[u8] = match aln.strand {
                Strand::Plus => seq,
                Strand::Minus => {
                    rc2 = revcomp(seq);
                    &rc2
                }
            };
            // E1: EA path already counted mm/ea; equivalent skip when
            // split/rescue left CIGAR+pos untouched.
            match ea_counts {
                Some((mm, ea)) if aln.cigar == ea_cigar && aln.pos == ea_pos => {
                    aln.mm_count = mm;
                    aln.ea_count = ea;
                }
                _ => self.recount_mm_ea(q, &mut aln),
            }
            // Coverage gate (variant-round products must pass; no
            // rescue-of-rescue).
            if (aln.rescued || cfg.extend.editing_aware) && read_coverage(&aln.cigar) < 0.7 {
                return None;
            }
            return Some(aln);
        }

        let mut aln = ReadAlignment {
            contig: contig_id,
            pos,
            strand,
            score,
            chain_score: best.score,
            second_chain_score: if second_ext_score > 0 {
                second_ext_score
            } else {
                second
            },
            cigar: ext.cigar,
            n_anchors,
            n_seeds: mins.len(),
            rescued: unedit_hit,
            ..Default::default()
        };
        {
            let rc2;
            let q: &[u8] = match aln.strand {
                Strand::Plus => seq,
                Strand::Minus => {
                    rc2 = revcomp(seq);
                    &rc2
                }
            };
            match ea_counts {
                Some((mm, ea)) if aln.cigar == ea_cigar && aln.pos == ea_pos => {
                    aln.mm_count = mm;
                    aln.ea_count = ea;
                }
                _ => self.recount_mm_ea(q, &mut aln),
            }
        }
        if (aln.rescued || cfg.extend.editing_aware) && read_coverage(&aln.cigar) < 0.7 {
            return None;
        }
        Some(aln)
    }

    /// Splice path: segment → junction refine → pseudo-reference extend →
    /// N CIGAR. splice gate frozen on: <2 segments ⇒ None (align_spliced would
    /// return None anyway — skip the revcomp + entry cost).
    fn try_spliced(&mut self, seq: &[u8], best: &Chain) -> Option<ReadAlignment> {
        let min_intron = SpliceParams::default().min_intron;
        let segs = segment_chain(best, self.config.seed.k, min_intron);
        if segs.len() < 2 {
            return None;
        }
        let empty_lib = JunctionLib::default();
        let lib = self.jlib.as_deref().unwrap_or(&empty_lib);
        let rc_storage;
        let query: &[u8] = match best.strand {
            Strand::Plus => seq,
            Strand::Minus => {
                rc_storage = revcomp(seq);
                &rc_storage
            }
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
        let mins = minimizers(seq, self.config.seed);
        if mins.is_empty() {
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
        let (base, _) = collect_anchors(self.index, &mins, read_len, k, 500);
        anchors.extend(base);
        anchors.retain(|a| a.contig == mate.contig && a.rpos >= wlo && a.rpos <= whi);
        if anchors.len() < 3 {
            return None;
        }
        let mut cp = self.config.chain;
        cp.min_chain_score = 20; // lower chaining threshold for rescue
        let chains = chain_anchors(anchors, &cp);
        let second = second_score(&chains);
        let best = chains.into_iter().next()?;
        // EA-scored extension (window = chain span ± read_len).
        let contig = &self.index.reference.contigs[best.contig as usize];
        let (rs, re) = best.ref_span(k);
        let flank = read_len + 50;
        let wstart = rs.saturating_sub(flank);
        let wend = (re + flank).min(contig.len).min(wstart + MAX_EXTEND_WINDOW);
        let window = contig.slice_ascii(wstart, wend);
        let rc_storage;
        let query: &[u8] = match best.strand {
            Strand::Plus => seq,
            Strand::Minus => {
                rc_storage = revcomp(seq);
                &rc_storage
            }
        };
        let hint = DiagHint {
            offset: (rs - wstart) as i64 - best.read_span(k).0 as i64,
            num: 1,
            den: 1,
        };
        let mut ep = self.config.extend;
        ep.editing_aware = true;
        let ext = extend_hint(query, &window, &ep, &mut self.buf, hint);
        if ext.score <= 0 {
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
            n_seeds: mins.len(),
            rescued: true,
            ..Default::default()
        };
        self.recount_mm_ea(query, &mut aln);
        Some(aln)
    }

    /// Split rescue (both tails; success rewrites CIGAR/junctions/pos).
    fn apply_split_rescue(&mut self, seq: &[u8], aln: ReadAlignment) -> ReadAlignment {
        let empty_lib = JunctionLib::default();
        let lib = self.jlib.as_deref().unwrap_or(&empty_lib);
        // Query in extension orientation: (RC'd if minus) read; tails operate
        // in extension coordinates; result CIGAR is plus-reference isomorphic.
        let rc_storage;
        let query: &[u8] = match aln.strand {
            Strand::Plus => seq,
            Strand::Minus => {
                rc_storage = revcomp(seq);
                &rc_storage
            }
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
                    // tail-rescue gate (frozen on): tails <12bp require known_support ≥ 2
                    // AND zero tail-vs-ref mismatch.
                    let t3_reject = r_tail < 12 && {
                        if rescue.junction.known_support < 2 {
                            true
                        } else {
                            let contig = &self.index.reference.contigs[out.contig as usize];
                            count_tail_mismatches(&rescue.cigar, out.pos, contig, query, true) > 0
                        }
                    };
                    if !t3_reject {
                        out.cigar = rescue.cigar;
                        out.junctions.push(rescue.junction);
                        out.score += rescue.score;
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
                    let t3_reject = read_start < 12 && {
                        if rescue.junction.known_support < 2 {
                            true
                        } else {
                            let contig = &self.index.reference.contigs[out.contig as usize];
                            count_tail_mismatches(
                                &rescue.cigar,
                                rescue.pos,
                                contig,
                                query,
                                false,
                            ) > 0
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
        out
    }

    /// Library-driven run-through reinterpretation (span.rs; skipped without
    /// a library or when the alignment is already spliced).
    fn apply_span_rescue(&mut self, seq: &[u8], mut aln: ReadAlignment) -> ReadAlignment {
        let Some(lib) = self.jlib.as_deref() else {
            return aln;
        };
        if !aln.junctions.is_empty() {
            return aln; // already has a spliced interpretation
        }
        let rc_storage;
        let query: &[u8] = match aln.strand {
            Strand::Plus => seq,
            Strand::Minus => {
                rc_storage = revcomp(seq);
                &rc_storage
            }
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

        let rc_storage;
        let query: &[u8] = match best.strand {
            Strand::Plus => seq,
            Strand::Minus => {
                rc_storage = revcomp(seq);
                &rc_storage
            }
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
                } else if ep.editing_aware
                    && ((rb == b'A' && qb == b'G') || (rb == b'T' && qb == b'C'))
                {
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
        let contig = &self.index.reference.contigs[aln.contig as usize];
        let mut scratch: Vec<u8> = Vec::new();
        let (mut ro, mut rf) = (0u32, aln.pos);
        let (mut mm, mut ea) = (0u32, 0u32);
        for op in &aln.cigar {
            match op {
                CigarOp::Match(n) => {
                    let s0 = scratch.len();
                    contig.decode_append(rf, rf + n, &mut scratch);
                    let got = scratch.len() - s0;
                    for k in 0..(*n as usize) {
                        let qb = query[ro as usize + k];
                        let rb = if k < got { scratch[s0 + k] } else { b'N' };
                        if qb != b'N' && rb != b'N' && qb != rb {
                            if (rb == b'A' && qb == b'G') || (rb == b'T' && qb == b'C') {
                                ea += 1;
                            } else {
                                mm += 1;
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
        aln.ea_count = ea;
    }

    pub(crate) fn chain_round(
        &mut self,
        mins: &[Minimizer],
        read_len: u32,
        k: u32,
        occ_cap: u32,
    ) -> Option<RoundOutcome> {
        let (anchors, _seeds_hit) = collect_anchors(self.index, mins, read_len, k, occ_cap);
        if anchors.is_empty() {
            return None;
        }
        let chains = chain_anchors(anchors, &self.config.chain);
        let second = second_score(&chains);
        Some(RoundOutcome { chains, second })
    }
}

/// Read coverage (M+I over read length).
pub(crate) fn read_coverage(cigar: &[CigarOp]) -> f64 {
    let (mut body, mut total) = (0u32, 0u32);
    for op in cigar {
        match op {
            CigarOp::Match(n) | CigarOp::Ins(n) => body += *n,
            CigarOp::SoftClip(n) => total += *n,
            _ => {}
        }
    }
    total += body;
    if total == 0 {
        return 0.0;
    }
    body as f64 / total as f64
}

/// Reverse complement (N stays N).
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| crate::fasta::Base::from_ascii(b).complement().to_ascii())
        .collect()
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
                        if ra != b'N' && qb != ra {
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
