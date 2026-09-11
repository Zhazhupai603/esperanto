//! # esperanto-engine
//!
//! L1 transcriptome-first matching engine: seed a read in a transcriptome
//! k-mer index, greedily extend each seed under the editing-aware (EA)
//! predicate, verify interrupted candidates with a bit-parallel EA-Myers
//! infix search, and project the winner deterministically onto the genome.
//!
//! Public API contract
//! -------------------
//! * [`align_read`] — pure function over one read; no global state, reads
//!   are fully independent.
//! * [`Tidx`] / [`TxMap`] / [`TxSeqs`] / [`RepeatTrack`] — data-source
//!   traits; [`L1Index`] implements the first three over one shared
//!   `tx_id` space.
//! * [`L1Index::build`] / [`L1Index::open`] / [`L1Index::save`] — runtime
//!   build (GTF + FASTA + `.tidx`) and the deterministic `L1BNDL01`
//!   bundle format.
//! * [`myers`] — EA-redefined bit-parallel Myers verifier (`infix`,
//!   `global`, and the two-block `long` variants for 128 < m <= 256).
//!
//! Invariants
//! ----------
//! * Determinism: identical inputs produce identical outcomes; every sort
//!   uses an explicit total-order key and no unsorted hash iteration
//!   reaches any output.
//! * Frozen production parameters live in [`EngineConfig::default`]; no
//!   environment switches exist.
//! * `k` comes from the index (`Tidx::k`, production 31), never from the
//!   engine config.

#![deny(unsafe_code)]

pub mod extend;
pub mod index;
pub mod kmer;
pub mod myers;
pub mod placement;
pub mod repeat;

use std::collections::{BTreeMap, BTreeSet};
use thiserror::Error;

pub use crate::extend::{extend_ea, Extension};
pub use crate::index::L1Index;
pub use crate::kmer::KmerStream;
pub use crate::placement::to_genomic_strand;
pub use crate::repeat::{NoRepeats, RepeatBed};

/// Errors raised by this crate.
#[derive(Debug, Error)]
pub enum Error {
    /// Underlying I/O failure.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// Bundle file does not start with the `L1BNDL01` magic.
    #[error("bad bundle magic: expected {expected:?}, found {found:?}")]
    BadMagic {
        /// Expected magic constant (`"L1BNDL01"`).
        expected: [u8; 8],
        /// Magic bytes found in the file.
        found: [u8; 8],
    },

    /// Unsupported bundle format version.
    #[error("unsupported bundle version: {0}")]
    UnsupportedVersion(u32),

    /// Malformed bundle or BED content.
    #[error("format error: {0}")]
    Format(String),

    /// Inputs disagree with each other (e.g. GTF vs `.tidx` counts).
    #[error("inconsistent inputs: {0}")]
    Inconsistent(String),

    /// Failure raised by `esperanto-tidx`.
    #[error("tidx error: {0}")]
    Tidx(String),

    /// Failure raised by `esperanto-txmap`.
    #[error("txmap error: {0}")]
    TxMap(String),
}

// ---------------------------------------------------------------------------
// Public basic types
// ---------------------------------------------------------------------------

/// Alignment strand of the read.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    /// Forward orientation.
    Plus,
    /// Reverse-complement orientation.
    Minus,
}

impl Strand {
    /// SAM FLAG bit: 0 for `+`, 16 for `-`.
    pub fn sam_flag(self) -> u16 {
        match self {
            Strand::Plus => 0,
            Strand::Minus => 16,
        }
    }

    /// Single-letter representation: `"+"` or `"-"`.
    pub fn letter(self) -> &'static str {
        match self {
            Strand::Plus => "+",
            Strand::Minus => "-",
        }
    }

    /// Flip the strand.
    pub fn flip(self) -> Strand {
        match self {
            Strand::Plus => Strand::Minus,
            Strand::Minus => Strand::Plus,
        }
    }
}

/// One CIGAR operation with SAM semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CigarOp {
    /// Alignment match (matches or mismatches, SAM `M`).
    Match(u32),
    /// Insertion to the reference (SAM `I`).
    Ins(u32),
    /// Deletion from the reference (SAM `D`).
    Del(u32),
    /// Reference skip across an intron (SAM `N`).
    RefSkip(u32),
    /// Soft-clipped bases (SAM `S`).
    SoftClip(u32),
}

impl CigarOp {
    /// SAM character of this operation.
    pub fn code(self) -> char {
        match self {
            CigarOp::Match(_) => 'M',
            CigarOp::Ins(_) => 'I',
            CigarOp::Del(_) => 'D',
            CigarOp::RefSkip(_) => 'N',
            CigarOp::SoftClip(_) => 'S',
        }
    }

    /// True when the operation length is zero.
    pub fn is_empty(self) -> bool {
        self.len() == 0
    }

    /// Operation length.
    pub fn len(self) -> u32 {
        match self {
            CigarOp::Match(n)
            | CigarOp::Ins(n)
            | CigarOp::Del(n)
            | CigarOp::RefSkip(n)
            | CigarOp::SoftClip(n) => n,
        }
    }
}

/// Render a CIGAR op list in SAM text form, e.g. `25M100N25M`.
pub fn cigar_string(ops: &[CigarOp]) -> String {
    use std::fmt::Write as _;
    let mut out = String::with_capacity(ops.len() * 4);
    for op in ops {
        let _ = write!(out, "{}{}", op.len(), op.code());
    }
    out
}

/// Outcome of L1 transcriptome-first matching for one read.
#[derive(Debug, Clone, PartialEq)]
pub enum L1Outcome {
    /// Read placed on the genome through a transcript projection.
    Aligned {
        /// Contig id (index into the projection's contig name table).
        contig: u32,
        /// Leftmost 0-based genomic position of the alignment.
        pos: u32,
        /// Genomic strand of the read.
        strand: Strand,
        /// CIGAR in reference (left-to-right) orientation.
        cigar: Vec<CigarOp>,
        /// Alignment score: 0 on the full branch, EA edit distance on the
        /// interrupted branch (lower is better).
        score: i32,
        /// Mapping quality 0..=60.
        mapq: u8,
    },
    /// No confident placement; defer to L2.
    Fallback,
}

/// Pipeline branch taken for a read (diagnostic, written into
/// [`ReadStats`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Branch {
    /// A full-length extension existed; projected directly.
    Full,
    /// Only partial extensions existed; placed after EA-Myers verification.
    Interrupted,
    /// Verification or the coverage/distance gates rejected the candidate.
    GateFail,
    /// No seed hit produced an extension.
    NoHit,
    /// Read shorter than the index k.
    TooShort,
}

/// Per-read diagnostic counters.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadStats {
    /// K-mer queries issued (every valid k-mer is queried; no coverage skipping).
    pub queries: u32,
    /// Read bases covered by the winning candidate's extension.
    pub extension_bases: u32,
    /// Branch taken by the read.
    pub branch: Branch,
    /// Number of distinct transcripts that received an extension.
    pub extension_tx_count: u32,
    /// Full-extension candidates (dist-0 incl. EA-free edits).
    pub fulls_count: u32,
    /// Partial-extension candidates (broken by indels/non-EA mismatches).
    pub parts_count: u32,
    /// Best verified partial distance (u32::MAX when none verified).
    pub best_part_dist: u32,
    /// Best verified partial projects to a different locus than the winner.
    pub part_conflict: bool,
}

impl Default for ReadStats {
    fn default() -> Self {
        ReadStats {
            queries: 0,
            extension_bases: 0,
            branch: Branch::NoHit,
            extension_tx_count: 0,
            fulls_count: 0,
            parts_count: 0,
            best_part_dist: u32::MAX,
            part_conflict: false,
        }
    }
}

/// Engine configuration. Production defaults are frozen; there are no
/// environment switches.
#[derive(Debug, Clone)]
pub struct EngineConfig {
    /// Flank added on both sides of the candidate window for verification.
    pub flank: usize,
    /// Minimum fraction of the read covered after verification.
    pub coverage_gate: f64,
    /// Number of partial candidates verified before picking the minimum.
    pub max_tx_candidates: usize,
    /// Maximum distinct diagonals retained per transcript.
    pub max_diagonals_per_tx: usize,
    /// K-mers prefetched per batch (pure performance hint, no semantics).
    pub prefetch_batch: usize,
    /// Global cap on raw index hits consumed per read.
    pub max_raw_hits: usize,
    /// A k-mer with more hits than this is skipped entirely.
    pub max_hits_per_kmer: usize,
}

impl Default for EngineConfig {
    fn default() -> Self {
        EngineConfig {
            flank: 10,
            coverage_gate: 0.90,
            max_tx_candidates: 8,
            max_diagonals_per_tx: 4,
            prefetch_batch: 16,
            max_raw_hits: 16384,
            max_hits_per_kmer: 100,
        }
    }
}

// ---------------------------------------------------------------------------
// EA predicate (directional ADAR rules)
// ---------------------------------------------------------------------------

/// Editing-aware (EA) pairing predicate: a reference/read base pair is
/// free when it is an identity match of two ACGT bases (case-insensitive)
/// or one of the directional editing pairs `ref A / read G` and
/// `ref T / read C`. The reverse conversions (G->A, C->T), `N`, IUPAC
/// ambiguity codes and every other byte always mismatch.
#[inline]
pub fn ea_free(ref_b: u8, read_b: u8) -> bool {
    let r = ref_b.to_ascii_uppercase();
    let q = read_b.to_ascii_uppercase();
    matches!(
        (r, q),
        (b'A', b'G') | (b'T', b'C') | (b'A', b'A') | (b'C', b'C') | (b'G', b'G') | (b'T', b'T')
    )
}

// ---------------------------------------------------------------------------
// Data-source traits
// ---------------------------------------------------------------------------

/// Transcriptome k-mer index (L1 seeding source).
pub trait Tidx {
    /// K-mer length the index was built with.
    fn k(&self) -> u32;
    /// Number of transcripts in the dense id space.
    fn tx_count(&self) -> u32;
    /// Look up a canonical k-mer; `(tx_id, offset)` pairs in ascending
    /// order, empty on miss.
    fn lookup(&self, canonical_kmer: u64) -> &[(u32, u32)];
    /// Display name of transcript `tx_id` (empty when out of range).
    fn transcript_name(&self, tx_id: u32) -> &str;
    /// Performance hint only; default no-op, no semantic effect.
    fn prefetch(&self, _: u64) {}
}

/// Transcript-to-genome projection source.
pub trait TxMap {
    /// Project transcript interval `[tx_start, tx_start + len)` onto the
    /// forward genome strand: `(contig_id, genomic_start, cigar)` with
    /// `Match`/`RefSkip` over genomically sorted pieces.
    fn project(&self, tx_id: u32, tx_start: u32, len: u32) -> Option<(u32, u32, Vec<CigarOp>)>;
    /// Total spliced length of the transcript.
    fn tx_len(&self, tx_id: u32) -> Option<u32>;
    /// Annotated strand of the transcript (None treated as `Plus`).
    fn strand(&self, _tx_id: u32) -> Option<Strand> {
        None
    }
}

/// Transcript sequence source in transcription orientation (minus-strand
/// transcripts already reverse-complemented).
pub trait TxSeqs {
    /// Sequence of `tx_id`; empty slice when unknown.
    fn seq(&self, tx_id: u32) -> &[u8];
}

/// Repeat region oracle for the marginal-placement gate.
pub trait RepeatTrack {
    /// True when `[pos, pos + len)` on contig `contig` overlaps a repeat.
    fn overlaps(&self, contig: u32, pos: u32, len: u32) -> bool;
}

// ---------------------------------------------------------------------------
// Engine: align_read
// ---------------------------------------------------------------------------

/// One retained seed extension (engine-internal).
#[derive(Debug, Clone)]
pub(crate) struct ExtEntry {
    pub tx_id: u32,
    /// Orientation of the oriented (possibly reverse-complemented) read
    /// against the transcript sequence.
    pub strand: Strand,
    /// `tx_lo - read_lo` of the extension (equals the seed diagonal).
    pub diagonal: i64,
    pub ext: Extension,
}

/// Run the L1 transcriptome-first pipeline on one read.
///
/// Full-length extensions are projected directly (score 0); otherwise the
/// top partial candidates are verified with EA-Myers, gated on coverage
/// (`coverage_gate`) and distance (`max(read_len/33, 1)`), and the winner
/// is attributed deterministically. Everything else is
/// [`L1Outcome::Fallback`].
pub fn align_read(
    read: &[u8],
    tidx: &impl Tidx,
    txmap: &impl TxMap,
    txseqs: &impl TxSeqs,
    cfg: &EngineConfig,
    repeats: &impl RepeatTrack,
    stats: &mut ReadStats,
) -> L1Outcome {
    let k = tidx.k() as usize;
    let read_len = read.len();
    if read_len < k {
        stats.branch = Branch::TooShort;
        return L1Outcome::Fallback;
    }
    stats.branch = Branch::NoHit;

    let read_rc = kmer::revcomp(read);

    // ---- interleaved seeding: query + extend (full information) ---------
    let mut exts: Vec<ExtEntry> = Vec::new();
    let mut seen: BTreeSet<(u32, u16, i64)> = BTreeSet::new();
    let mut diag_per_tx: BTreeMap<u32, usize> = BTreeMap::new();
    let mut ext_per_tx: BTreeMap<u32, usize> = BTreeMap::new();
    let mut raw_hits: usize = 0;
    let per_tx_ext_cap = cfg.max_tx_candidates * cfg.max_diagonals_per_tx;
    let batch_size = cfg.prefetch_batch.max(1);

    let mut stream = KmerStream::new(read, k);
    let mut batch: Vec<(usize, u64)> = Vec::with_capacity(batch_size);
    'seeding: loop {
        batch.clear();
        while batch.len() < batch_size {
            match stream.next() {
                Some(item) => batch.push(item),
                None => break,
            }
        }
        if batch.is_empty() {
            break;
        }
        for &(_, code) in &batch {
            tidx.prefetch(code);
        }
        let exhausted = batch.len() < batch_size;
        for (pos, code) in batch.iter().copied() {
            let canon = kmer::canonical(code, k);
            let hits = tidx.lookup(canon);
            stats.queries += 1;
            if hits.len() > cfg.max_hits_per_kmer {
                continue;
            }
            for &(tx_id, tx_off) in hits {
                if raw_hits >= cfg.max_raw_hits {
                    break 'seeding;
                }
                raw_hits += 1;
                if ext_per_tx.get(&tx_id).copied().unwrap_or(0) >= per_tx_ext_cap {
                    continue;
                }
                let tx_seq = txseqs.seq(tx_id);
                let t_off = tx_off as usize;
                if t_off + k > tx_seq.len() {
                    continue;
                }
                let Some(tx_code) = kmer::window_code(tx_seq, t_off, k) else {
                    continue;
                };
                // Strand: forward code comparison at the hit offset.
                let (strand, oriented_pos) = if tx_code == code {
                    (Strand::Plus, pos)
                } else {
                    (Strand::Minus, read_len - k - pos)
                };
                let diagonal = tx_off as i64 - oriented_pos as i64;
                let key = (tx_id, strand.sam_flag(), diagonal);
                if seen.contains(&key) {
                    continue;
                }
                if diag_per_tx.get(&tx_id).copied().unwrap_or(0) >= cfg.max_diagonals_per_tx {
                    continue;
                }
                let oriented: &[u8] = if strand == Strand::Plus {
                    read
                } else {
                    &read_rc
                };
                let ext = extend_ea(oriented, oriented_pos, k, tx_seq, t_off);
                seen.insert(key);
                *diag_per_tx.entry(tx_id).or_insert(0) += 1;
                *ext_per_tx.entry(tx_id).or_insert(0) += 1;
                exts.push(ExtEntry {
                    tx_id,
                    strand,
                    diagonal,
                    ext,
                });
            }
        }
        if exhausted {
            break;
        }
    }
    stats.extension_tx_count = ext_per_tx.len() as u32;

    if exts.is_empty() {
        stats.branch = Branch::NoHit;
        return L1Outcome::Fallback;
    }

    let mut fulls: Vec<ExtEntry> = Vec::new();
    let mut parts: Vec<ExtEntry> = Vec::new();
    for e in exts {
        if e.ext.full {
            fulls.push(e);
        } else {
            parts.push(e);
        }
    }
    parts.sort_by(|a, b| {
        b.ext
            .read_cov()
            .cmp(&a.ext.read_cov())
            .then(a.tx_id.cmp(&b.tx_id))
            .then(a.diagonal.cmp(&b.diagonal))
    });

    // EA-Myers infix distance of one candidate against its anchored
    // transcript window (None when the window is empty).
    let verify = |e: &ExtEntry| -> Option<u32> {
        let tx_seq = txseqs.seq(e.tx_id);
        let anchor = e.ext.tx_lo as i64 - e.ext.read_lo as i64;
        let wlo = (anchor - cfg.flank as i64).max(0) as usize;
        let whi = ((anchor + read_len as i64 + cfg.flank as i64).max(0) as usize).min(tx_seq.len());
        if whi <= wlo {
            return None;
        }
        let oriented: &[u8] = if e.strand == Strand::Plus {
            read
        } else {
            &read_rc
        };
        Some(if read_len <= 128 {
            myers::infix(oriented, &tx_seq[wlo..whi])
        } else {
            myers::long::infix(oriented, &tx_seq[wlo..whi])
        })
    };
    stats.fulls_count = fulls.len() as u32;
    stats.parts_count = parts.len() as u32;

    // ---- branch 2: full extensions --------------------------------------
    if !fulls.is_empty() {
        stats.branch = Branch::Full;
        fulls.sort_by(|a, b| {
            (a.strand.sam_flag(), a.tx_id, a.diagonal).cmp(&(b.strand.sam_flag(), b.tx_id, b.diagonal))
        });
        let mut placed: Vec<placement::Placed> = fulls
            .iter()
            .filter_map(|e| placement::project_full(txmap, e, read_len))
            .collect();
        if placed.is_empty() {
            stats.branch = Branch::GateFail;
            return L1Outcome::Fallback;
        }
        placed.sort_by(|a, b| {
            (a.score, a.tx_id, a.diagonal).cmp(&(b.score, b.tx_id, b.diagonal))
        });
        let Some((best, mapq)) = placement::finalize_candidates(&placed, false) else {
            stats.branch = Branch::GateFail;
            return L1Outcome::Fallback;
        };
        // Cross-branch competition probe: a full (dist-0) winner can be a
        // repeat decoy when the read's true transcript sits in the PARTIAL
        // list (indel-broken extension). Verify the top partials; a verified
        // partial at a different locus within the dist gate contests the win.
        let max_dist = (read_len / 33).max(1) as u32;
        for entry in parts.iter().take(3) {
            let Some(dist) = verify(entry) else {
                continue;
            };
            if dist < stats.best_part_dist {
                stats.best_part_dist = dist;
            }
            if dist <= max_dist {
                if let Some((c, p)) = placement::partial_locus(txmap, entry, read_len) {
                    if c != best.contig || p.abs_diff(best.pos) > 1_000 {
                        stats.part_conflict = true;
                    }
                }
            }
        }
        stats.extension_bases = read_len as u32;
        return L1Outcome::Aligned {
            contig: best.contig,
            pos: best.pos,
            strand: best.strand,
            cigar: best.cigar.clone(),
            score: best.score,
            mapq,
        };
    }

    // ---- branch 3: partials, verification --------------------------------
    stats.branch = Branch::Interrupted;

    // Winner: the FIRST candidate in parts order (read_cov DESC, tx_id
    // ASC, diagonal ASC) whose dist is strictly smaller than every
    // previously verified one — the first minimum wins and the list is
    // never re-sorted afterwards.
    let take = cfg.max_tx_candidates.min(parts.len());
    let mut winner: Option<(&ExtEntry, u32)> = None;
    if read_len <= 256 {
        for entry in parts.iter().take(take) {
            if let Some(dist) = verify(entry) {
                let better = match winner {
                    Some((_, best)) => dist < best,
                    None => true,
                };
                if better {
                    winner = Some((entry, dist));
                }
            }
        }
    }
    let Some((winner_entry, winner_dist)) = winner else {
        stats.branch = Branch::GateFail;
        return L1Outcome::Fallback;
    };
    stats.best_part_dist = winner_dist;

    // ---- branch 4: gates ---------------------------------------------------
    let cov = (read_len - winner_dist as usize) as f64 / read_len as f64;
    if cov < cfg.coverage_gate {
        stats.branch = Branch::GateFail;
        return L1Outcome::Fallback;
    }
    let max_dist = (read_len / 33).max(1) as u32;
    if winner_dist > max_dist {
        stats.branch = Branch::GateFail;
        return L1Outcome::Fallback;
    }

    // ---- branch 5: R5 cluster competition ---------------------------------
    // Clusters come from the FULL sorted parts list, collapsed by first
    // occurrence of (tx_id, strand); cluster order = parts order (no
    // re-sort by distance). Competitor distances are computed lazily —
    // a cluster may sit far below the verified top-N.
    let mut force_mapq0 = false;
    let mut cluster_seen: BTreeSet<(u32, u16)> = BTreeSet::new();
    let mut clusters: Vec<&ExtEntry> = Vec::new();
    for e in &parts {
        if cluster_seen.insert((e.tx_id, e.strand.sam_flag())) {
            clusters.push(e);
        }
    }
    if clusters.len() >= 2 && winner_dist >= 2 {
        if let Some(best_locus) = placement::partial_locus(txmap, winner_entry, read_len) {
            for c in clusters.iter().skip(1) {
                let Some(d_c) = verify(c) else {
                    continue;
                };
                let cov_c = (read_len - d_c as usize) as f64 / read_len as f64;
                if cov_c < cfg.coverage_gate || d_c > max_dist {
                    continue;
                }
                let Some(c_locus) = placement::partial_locus(txmap, c, read_len) else {
                    continue;
                };
                if c_locus == best_locus {
                    continue;
                }
                if (winner_dist as i64 - d_c as i64).abs() <= 2 {
                    force_mapq0 = true;
                    break;
                }
            }
        }
    }

    // ---- branch 6: repeat-region marginal placement -----------------------
    if winner_dist >= max_dist {
        if let Some((contig, pos)) = placement::partial_locus(txmap, winner_entry, read_len) {
            if repeats.overlaps(contig, pos, read_len as u32) {
                force_mapq0 = true;
            }
        }
    }

    // ---- branch 7: finalize (the single winner) ----------------------------
    let Some(bestp) = placement::project_partial(txmap, winner_entry, read_len, winner_dist as i32)
    else {
        stats.branch = Branch::GateFail;
        return L1Outcome::Fallback;
    };
    let mapq = placement::finalize_candidates(std::slice::from_ref(&bestp), force_mapq0)
        .map_or(60, |(_, q)| q);
    stats.extension_bases = winner_entry.ext.read_cov() as u32;
    L1Outcome::Aligned {
        contig: bestp.contig,
        pos: bestp.pos,
        strand: bestp.strand,
        cigar: bestp.cigar.clone(),
        score: bestp.score,
        mapq,
    }
}

