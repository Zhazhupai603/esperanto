//! MEM intron-chain fast path (two-segment spliced reads).
//!
//! Reads whose anchors split into two collinear clusters with a reference gap
//! shaped like an intron are placed directly: anchors are extended into
//! maximal exact matches (MEMs), MEMs cluster into fragments by implied
//! read-start (±5bp), a dominant-fragment gate (≥55% of anchored coverage)
//! keeps only concentrated configurations, the best fragment pair is scored
//! `cov*4 + anchors*3 − bal`, the middle gap is scanned for splice signals
//! (with naive fallback), EA-Myers sharpens both breakpoints, the split is
//! verified under dual mismatch-rate thresholds, and the final CIGAR is
//! guarded for read-length conservation. Any failing stage returns `None`;
//! the caller falls through to later paths.

use crate::extend::{
    CigarOp, ExtendBuffer, ExtendParams, Extension,
};
use crate::fasta::Reference;
use crate::gtf::{Junction, JunctionLib, RefinedJunction};
use crate::myers_ea::{infix_best_end, infix_best_start};
use crate::seed::{Anchor, Strand};
use crate::splice::{
    junction_signal, refine_junction, SplicedAlignment, SpliceParams,
};

/// Intron-chain parameters (frozen; see parameter table).
#[derive(Clone, Copy, Debug)]
pub struct IntronParams {
    /// Minimum fragment read coverage to participate.
    pub min_frag_cov: u32,
    /// EA-Myers breakpoint refinement half-window.
    pub refine_window: u32,
    /// Breakpoint refinement pattern length (read bases per side).
    pub refine_pattern_len: u32,
    /// Maximum intron length.
    pub intron_max: u32,
}

impl Default for IntronParams {
    fn default() -> IntronParams {
        IntronParams {
            min_frag_cov: 15,
            refine_window: 30,
            refine_pattern_len: 50,
            intron_max: 105_000,
        }
    }
}

/// Minimum intron length (matches the chain RNA `min_intron`).
pub const MIN_INTRON: u32 = 20;

/// Fragment clustering tolerance on the implied read-start (frozen: ±5bp).
pub const FRAG_CLUSTER_TOL: i64 = 5;

/// Dominant-fragment gate fraction numerator/denominator (≥55%).
pub const DOMINANT_FRAC: (u64, u64) = (55, 100);

/// Per-segment mismatch-rate ceiling ×100 (≤0.17).
pub const SEG_MM_RATE_X100: u64 = 17;

/// Combined mismatch-rate ceiling ×100 (r1 + r2 ≤ 0.24).
pub const SUM_MM_RATE_X100: u64 = 24;

/// One maximal exact match between read and reference (alignment
/// orientation; `q`/`r` are half-open read/reference spans).
#[derive(Clone, Copy, Debug)]
pub struct MemHit {
    /// Contig index.
    pub contig: u32,
    /// Alignment strand.
    pub strand: Strand,
    /// Read span start (inclusive).
    pub q_lo: u32,
    /// Read span end (exclusive).
    pub q_hi: u32,
    /// Reference span start (inclusive).
    pub r_lo: u32,
    /// Reference span end (exclusive).
    pub r_hi: u32,
}

/// A cluster of MEMs sharing an implied read-start within ±5bp.
#[derive(Clone, Debug)]
pub struct Frag {
    /// Contig index.
    pub contig: u32,
    /// Alignment strand.
    pub strand: Strand,
    /// Member MEMs sorted by `q_lo`.
    pub mems: Vec<MemHit>,
    /// Leftmost read position covered.
    pub q_lo: u32,
    /// One past the rightmost read position covered.
    pub q_hi: u32,
    /// Leftmost reference position covered.
    pub r_lo: u32,
    /// One past the rightmost reference position covered.
    pub r_hi: u32,
}

impl Frag {
    /// Read bases covered (union of member MEM read spans).
    pub fn cov(&self) -> u32 {
        let mut total = 0u32;
        let mut end = 0u32;
        for m in &self.mems {
            if m.q_hi > end {
                total += m.q_hi - m.q_lo.max(end);
                end = m.q_hi;
            }
        }
        total
    }

    /// Number of anchors backing the fragment (MEM count).
    pub fn n_anchors(&self) -> usize {
        self.mems.len()
    }

    /// Reference offset of the fragment's leftmost MEM (`r_lo − q_lo`).
    fn off_lo(&self) -> i64 {
        self.r_lo as i64 - self.q_lo as i64
    }
}

fn base_eq(a: u8, b: u8) -> bool {
    let a = a & 0xDF;
    let b = b & 0xDF;
    (b'A'..=b'T').contains(&a) && a == b && matches!(a, b'A' | b'C' | b'G' | b'T')
}

/// Extend each anchor's seed match into a maximal exact match.
///
/// The read must be in alignment orientation (anchor `qpos` space). Extension
/// stops at read/reference bounds, at an `N`, or at the first mismatch. The
/// reference window is decoded once per anchor.
pub fn extend_mems(reference: &Reference, read: &[u8], anchors: &[Anchor], k: u32) -> Vec<MemHit> {
    let rl = read.len() as u32;
    let mut out = Vec::with_capacity(anchors.len());
    for a in anchors {
        let Some(ctg) = reference.contigs.get(a.contig as usize) else {
            continue;
        };
        let win_lo = a.rpos.saturating_sub(rl);
        let win_hi = (a.rpos as u64 + k as u64 + rl as u64).min(ctg.len as u64) as u32;
        if win_hi <= win_lo {
            continue;
        }
        let window = ctg.slice_ascii(win_lo, win_hi);
        let wf = |rpos: u32| -> usize { (rpos - win_lo) as usize };

        // left extension
        let mut q_lo = a.qpos;
        let mut r_lo = a.rpos;
        while q_lo > 0 && r_lo > win_lo && base_eq(read[(q_lo - 1) as usize], window[wf(r_lo - 1)])
        {
            q_lo -= 1;
            r_lo -= 1;
        }
        // right extension
        let seed_hi_q = a.qpos + k;
        let seed_hi_r = a.rpos + k;
        let mut q_hi = seed_hi_q;
        let mut r_hi = seed_hi_r;
        while q_hi < rl && r_hi < win_hi && base_eq(read[q_hi as usize], window[wf(r_hi)]) {
            q_hi += 1;
            r_hi += 1;
        }
        if q_hi > q_lo {
            out.push(MemHit {
                contig: a.contig,
                strand: a.strand,
                q_lo,
                q_hi,
                r_lo,
                r_hi,
            });
        }
    }
    out
}

/// Cluster MEMs into fragments by implied read-start (`q_lo − r_lo`) within
/// ±5bp, per (contig, strand). Fragments come out sorted by
/// (contig, strand, diag, q_lo).
pub fn cluster_frags(mems: Vec<MemHit>) -> Vec<Frag> {
    let mut mems = mems;
    mems.sort_by_key(|m| {
        (
            m.contig,
            m.strand as u8,
            m.q_lo as i64 - m.r_lo as i64,
            m.q_lo,
            m.r_lo,
        )
    });
    let mut frags: Vec<Frag> = Vec::new();
    for m in mems {
        let diag = m.q_lo as i64 - m.r_lo as i64;
        let start_new = match frags.last() {
            Some(f) if f.contig == m.contig && f.strand == m.strand => {
                diag - (f.q_lo as i64 - f.r_lo as i64) > FRAG_CLUSTER_TOL
            }
            _ => true,
        };
        if start_new {
            frags.push(Frag {
                contig: m.contig,
                strand: m.strand,
                mems: vec![m],
                q_lo: m.q_lo,
                q_hi: m.q_hi,
                r_lo: m.r_lo,
                r_hi: m.r_hi,
            });
        } else {
            let f = frags.last_mut().unwrap();
            f.q_lo = f.q_lo.min(m.q_lo);
            f.q_hi = f.q_hi.max(m.q_hi);
            f.r_lo = f.r_lo.min(m.r_lo);
            f.r_hi = f.r_hi.max(m.r_hi);
            f.mems.push(m);
        }
    }
    for f in &mut frags {
        f.mems.sort_by_key(|m| m.q_lo);
    }
    frags
}

/// Dominant-fragment gate: the largest fragment must hold ≥55% of the total
/// anchored coverage.
pub fn dominant_gate(frags: &[Frag]) -> bool {
    let mut total = 0u64;
    let mut max_cov = 0u64;
    for f in frags {
        let c = f.cov() as u64;
        total += c;
        max_cov = max_cov.max(c);
    }
    if total == 0 {
        return false;
    }
    max_cov * DOMINANT_FRAC.1 >= total * DOMINANT_FRAC.0
}

/// A selected fragment pair with its naive intron interval.
#[derive(Clone, Debug)]
pub struct FragPair {
    /// Read-upstream fragment (`A`).
    pub left: Frag,
    /// Read-downstream fragment (`B`).
    pub right: Frag,
    /// Pair score `cov*4 + anchors*3 − bal` (`bal = |covA − covB|`).
    pub score: i32,
    /// Read split point (`min(a.q_hi, b.q_lo)`).
    pub split: u32,
    /// Naive intron start (first intron base) from the diagonals.
    pub intron_start: u32,
    /// Naive intron end (first base of the right exon).
    pub intron_end: u32,
}

/// Pick the best fragment pair (same contig+strand, read-ordered, intron
/// within [MIN_INTRON, intron_max]) by `cov*4 + anchors*3 − bal`.
pub fn best_frag_pair(frags: &[Frag], intron_max: u32) -> Option<FragPair> {
    let mut best: Option<FragPair> = None;
    for i in 0..frags.len() {
        for j in 0..frags.len() {
            if i == j {
                continue;
            }
            let (a, b) = (&frags[i], &frags[j]);
            // MEM edges may overlap by a couple of chance-matching bases at
            // the seam; tolerate up to 3bp (frozen legacy value).
            if a.contig != b.contig
                || a.strand != b.strand
                || a.q_hi > b.q_lo + 3
            {
                continue;
            }
            // Split at the earlier fragment edge so both segments stay
            // inside their exons.
            let split = a.q_hi.min(b.q_lo);
            let off_a = a.r_hi as i64 - a.q_hi as i64;
            let off_b = b.off_lo();
            let len0 = off_b - off_a;
            if len0 < MIN_INTRON as i64 || len0 > intron_max as i64 {
                continue;
            }
            let (cov_a, cov_b) = (a.cov(), b.cov());
            let score = ((cov_a + cov_b) * 4) as i32 + ((a.n_anchors() + b.n_anchors()) * 3) as i32
                - (cov_a as i64 - cov_b as i64).abs() as i32;
            let start = (split as i64 + off_a) as u32;
            let end = (split as i64 + off_b) as u32;
            if best.as_ref().is_none_or(|p| score > p.score) {
                best = Some(FragPair {
                    left: a.clone(),
                    right: b.clone(),
                    score,
                    split,
                    intron_start: start,
                    intron_end: end,
                });
            }
        }
    }
    best
}

/// Scan the middle gap for splice signals via [`refine_junction`] with the
/// wider `refine_window` radius; fall back to the naive breakpoints when the
/// scan cannot produce an in-bounds candidate.
#[allow(clippy::too_many_arguments)]
pub fn scan_intron_gap(
    reference: &Reference,
    lib: &JunctionLib,
    contig: u32,
    minus: bool,
    naive_start: u32,
    naive_end: u32,
    splice_params: &SpliceParams,
    refine_window: u32,
) -> RefinedJunction {
    let mut scan_params = *splice_params;
    scan_params.refine_radius = refine_window;
    refine_junction(
        reference,
        contig,
        naive_start,
        naive_end,
        minus,
        lib,
        &scan_params,
    )
    .unwrap_or_else(|| {
        let ctg = &reference.contigs[contig as usize];
        let junction = Junction {
            contig,
            start: naive_start,
            end: naive_end,
            minus_strand: minus,
        };
        let signal = junction_signal(ctg, &junction);
        RefinedJunction {
            junction,
            signal,
            known_support: lib.support(&junction),
        }
    })
}

/// EA-Myers donor refinement (frozen): `infix_best_end` of the read pattern
/// ending at the junction over the window `[frag_ref_hi ± refine_window]`.
/// Pattern = `query[pat_start..junction_off]` with `pat_start = max(junction_off
/// − refine_pattern_len, text_lo_in_read)` where `text_lo_in_read = (lo −
/// frag_implied_start + 2).max(0)`. Empty window/pattern ⇒ `frag_ref_hi − 1`.
/// The returned coordinate is adopted directly (no signal comparison).
#[allow(clippy::too_many_arguments)]
pub fn refine_donor_end(
    ctg: &crate::fasta::Contig,
    read: &[u8],
    junction_off: u32,
    frag_ref_hi: u32,
    frag_implied_start: i64,
    params: &IntronParams,
) -> u32 {
    let w = params.refine_window as i64;
    let plen = params.refine_pattern_len as usize;
    let off = (junction_off as usize).min(read.len());
    let lo = (frag_ref_hi as i64 - w).max(0) as u32;
    let hi = (frag_ref_hi as i64 + w).min(ctg.len as i64) as u32;
    if hi <= lo || off == 0 {
        return frag_ref_hi.saturating_sub(1);
    }
    let text_lo_in_read = (lo as i64 - frag_implied_start + 2).max(0) as usize;
    let pat_start = off.saturating_sub(plen).max(text_lo_in_read);
    if pat_start >= off {
        return frag_ref_hi.saturating_sub(1);
    }
    let pattern = &read[pat_start..off];
    let text = ctg.slice_ascii(lo, hi);
    let (_dist, end) = infix_best_end(pattern, &text);
    lo + end as u32
}

/// Mirror of [`refine_donor_end`] for the acceptor: `infix_best_start` of the
/// read pattern starting at the junction over `[frag_ref_lo ± refine_window]`.
/// Pattern = `query[junction_off..pat_end]` with `pat_end = min(junction_off +
/// refine_pattern_len, read_len, text_hi_in_read)` where `text_hi_in_read =
/// (hi − frag_implied_start − 2).max(0)`. Empty window/pattern ⇒ `frag_ref_lo`.
#[allow(clippy::too_many_arguments)]
pub fn refine_acceptor_start(
    ctg: &crate::fasta::Contig,
    read: &[u8],
    junction_off: u32,
    frag_ref_lo: u32,
    frag_implied_start: i64,
    params: &IntronParams,
) -> u32 {
    let w = params.refine_window as i64;
    let plen = params.refine_pattern_len as usize;
    let off = (junction_off as usize).min(read.len());
    let lo = (frag_ref_lo as i64 - w).max(0) as u32;
    let hi = (frag_ref_lo as i64 + w).min(ctg.len as i64) as u32;
    if hi <= lo {
        return frag_ref_lo;
    }
    let text_hi_in_read = (hi as i64 - frag_implied_start - 2).max(0) as usize;
    let pat_end = (off + plen).min(read.len()).min(text_hi_in_read);
    if pat_end <= off {
        return frag_ref_lo;
    }
    let pattern = &read[off..pat_end];
    let text = ctg.slice_ascii(lo, hi);
    let (_dist, start) = infix_best_start(pattern, &text);
    lo + start as u32
}



/// Per-fragment EA-Myers split-read verification (frozen legacy): returns
/// (dist_donor, dist_acceptor) — real junctions score low, spurious high.
fn verify_split_read(
    read: &[u8],
    ctg: &crate::fasta::Contig,
    frag1: &Frag,
    frag2: &Frag,
    donor_end: i64,
    acceptor_start: i64,
) -> (i32, i32) {
    let mid = read.len() / 2;
    let d_start = (frag1.off_lo() - 30).max(0) as u32;
    let d_end = (donor_end + 1 + 30).max(1).min(ctg.len as i64) as u32;
    let a_start = (acceptor_start - 30).max(0).min(ctg.len as i64) as u32;
    let a_end = (frag2.r_hi as i64 + 30).max(a_start as i64 + 1).min(ctg.len as i64) as u32;
    if d_start >= d_end || a_start >= a_end || mid == 0 || mid >= read.len() {
        return (i32::MAX, i32::MAX);
    }
    let read_donor = &read[..mid];
    let read_acceptor = &read[mid..];
    let donor_text = ctg.slice_ascii(d_start, d_end);
    let acceptor_text = ctg.slice_ascii(a_start, a_end);
    // long::infix dispatches m<=128 to the single-block path internally.
    let dist_donor = crate::myers_ea::long::infix(read_donor, &donor_text);
    let dist_acceptor = crate::myers_ea::long::infix(read_acceptor, &acceptor_text);
    (dist_donor, dist_acceptor)
}

/// Conservation guard: CIGAR `M/I/S` ops consume exactly `read_len` bases.
pub fn cigar_conserving(cigar: &[CigarOp], read_len: usize) -> bool {
    let mut total = 0u64;
    for op in cigar {
        match op {
            CigarOp::Match(n) | CigarOp::Ins(n) | CigarOp::SoftClip(n) => total += *n as u64,
            CigarOp::Del(_) | CigarOp::RefSkip(_) => {}
        }
    }
    total == read_len as u64
}

/// Two-segment spliced placement over the anchor set; `None` at any failing
/// stage. The read must be in alignment orientation (anchor `qpos` space).
/// EA-refined breakpoints are adopted directly; `verify_split_read` is the
/// only acceptance stage.
#[allow(clippy::too_many_arguments)]
pub fn try_intron_chain_placement(
    reference: &Reference,
    _lib: &JunctionLib,
    read: &[u8],
    anchors: &[Anchor],
    k: u32,
    params: &IntronParams,
    _extend_params: &ExtendParams,
    _buf: &mut ExtendBuffer,
) -> Option<SplicedAlignment> {
    if anchors.is_empty() || read.is_empty() {
        return None;
    }
    let mems = extend_mems(reference, read, anchors, k);
    let probe = std::env::var_os("ESP_PROBE").is_some();
    if probe {
        eprintln!("[probe-ic] anchors={} mems={}", anchors.len(), mems.len());
    }
    let frags: Vec<Frag> = cluster_frags(mems)
        .into_iter()
        .filter(|f| f.cov() >= params.min_frag_cov)
        .collect();
    // Old parity: the dominant-fragment check is a CONTINUOUS fast path,
    // not a rejection gate — a fragment covering ≥55% of the READ means an
    // unspliced read (caller falls through to other paths); otherwise pair
    // finding proceeds regardless of dominance.
    if frags.is_empty() {
        return None;
    }
    if let Some(dominant) = frags.iter().max_by_key(|f| f.cov()) {
        if dominant.cov() as usize * 20 >= read.len() * 11 {
            if probe {
                eprintln!("[probe-ic] continuous fast path (dominant {} >= 55% of read)", dominant.cov());
            }
            return None;
        }
    }
    if probe {
        eprintln!("[probe-ic] frags_after_cov={} (no dominant; pair finding)", frags.len());
    }
    let pair = best_frag_pair(&frags, params.intron_max)?;
    if probe {
        eprintln!("[probe-ic] pair: split={} left_cov={} right_cov={}", pair.split, pair.left.cov(), pair.right.cov());
    }
    let left = &pair.left;
    let right = &pair.right;
    let minus = left.strand == Strand::Minus;
    let ctg = reference.contigs.get(left.contig as usize)?;

    // Junction offset for refinement: read midpoint clamped into the
    // inter-fragment read gap (frozen legacy).
    let j_est = (read.len() / 2) as u32;
    let junction_off = if left.q_hi >= right.q_lo {
        j_est
    } else {
        j_est.clamp(left.q_hi, right.q_lo)
    };
    let donor_end = refine_donor_end(ctg, read, junction_off, left.r_hi, left.off_lo(), params) as i64;
    let acceptor_start = refine_acceptor_start(ctg, read, junction_off, right.r_lo, right.off_lo(), params) as i64;
    let intron_len = acceptor_start - donor_end - 1;

    // Split-read verification (dual threshold, frozen legacy).
    let (vd, va) = verify_split_read(read, ctg, left, right, donor_end, acceptor_start);
    let total_anchors = left.n_anchors() + right.n_anchors();
    let max_v = vd.max(va);
    let sum_v = vd + va;
    let max_threshold = (read.len() as f64 * 0.17) as i32;
    let sum_threshold = (read.len() as f64 * 0.24) as i32;
    if total_anchors < 2 {
        return None;
    }
    if intron_len < MIN_INTRON as i64 || intron_len > params.intron_max as i64 {
        return None;
    }
    if max_v > max_threshold || sum_v > sum_threshold {
        return None;
    }
    if probe {
        eprintln!("[probe-ic] verify PASS vd={} va={} donor_end={} acceptor_start={}", vd, va, donor_end, acceptor_start);
    }

    // Assembly (frozen legacy): [Match(read_break)] N(intron) [Match(rest)]
    // — pure M blocks, no soft clips, no indels.
    let read_break = ((left.q_hi + right.q_lo) / 2) as usize;
    let frag1_len = read_break;
    let frag2_len = read.len() - read_break;
    let mut cigar = Vec::with_capacity(3);
    if frag1_len > 0 {
        cigar.push(CigarOp::Match(frag1_len as u32));
    }
    cigar.push(CigarOp::RefSkip(intron_len as u32));
    if frag2_len > 0 {
        cigar.push(CigarOp::Match(frag2_len as u32));
    }
    if !cigar_conserving(&cigar, read.len()) {
        return None;
    }

    let refined = RefinedJunction {
        junction: Junction {
            contig: left.contig,
            start: (donor_end + 1) as u32,
            end: acceptor_start as u32,
            minus_strand: minus,
        },
        signal: crate::gtf::SpliceSignal::NonCanonical,
        known_support: 0,
    };
    let pos = left.off_lo().max(0) as u32;
    let read_len32 = read.len() as u32;
    let exon_bounds = vec![
        (0u32, pos, pos + frag1_len as u32),
        (frag1_len as u32, acceptor_start as u32, acceptor_start as u32 + frag2_len as u32),
    ];
    Some(SplicedAlignment {
        extension: Extension {
            read_start: 0,
            read_end: read_len32,
            ref_start: 0,
            cigar: cigar.clone(),
            score: read_len32 as i32 * 2,
        },
        junctions: vec![refined],
        cigar,
        exon_bounds,
    })
}

