//! Tail split rescue (STAR-style split discovery).
//!
//! When the primary alignment leaves a soft-clipped end tail of at least
//! `TAIL_MIN` bases, the tail is relocated across a splice junction. Channel A
//! (library-driven, the main channel) enumerates library junctions in the
//! donor/acceptor ±20bp window and retries the tail start ±5bp, scoring both
//! segments by extension. Channels B (direct genomic scan) and C (tail
//! seeding) are de-novo discovery and are entirely disabled when the library
//! is non-empty (97% of false junctions came from support-1 de-novo calls).
//!
//! All CIGAR products conserve the read length exactly (`M/I/S` sums).

use crate::extend::{
    extend_hint, push_op, CigarOp, DiagHint, ExtendBuffer, ExtendParams, Extension,
};
use crate::fasta::Reference;
use crate::gtf::{Junction, JunctionLib, RefinedJunction, SpliceSignal};
use crate::index::Index;
use crate::intron_chain::cigar_conserving;
use crate::seed::{minimizers, SeedParams, Strand};
use crate::splice::{
    junction_signal, lib_window_candidates, refine_junction, s14_pass, splice_signal, tail_gate,
    SpliceParams, SPLIT_COST, SPLIT_LIB_WINDOW,
};

/// Library-driven tail lower bound (frozen).
pub const TAIL_MIN: usize = 7;

/// Tail-start retry radius (parity experiments, frozen: ±5bp).
pub const TAIL_RETRY_RADIUS: i64 = 5;

/// Direct-scan channel tail length range (5..19bp; tails ≥ 19 bp reach the
/// seed channel, which needs k+w-1 = 19 bp for one minimizer).
pub const DIRECT_TAIL_MIN: usize = 5;
pub const DIRECT_TAIL_MAX: usize = 19;

/// Direct-scan donor search radius (±12bp).
pub const DIRECT_DONOR_RADIUS: i64 = 12;

/// Direct-scan downstream search span.
pub const DIRECT_SCAN_SPAN: u32 = 50_000;

/// Exact prefix probe length cap.
pub const PROBE_MAX: usize = 10;

/// Seed channel minimum tail (≥ k = 15).
pub const SEED_TAIL_MIN: usize = 15;

/// Seed channel: first N tail minimizers.
pub const SEED_CHANNEL_MINS: usize = 8;

/// Seed channel occupancy cap (`0 < count ≤ 100`).
pub const SEED_CHANNEL_OCC_CAP: u32 = 100;

/// Seed channel positions cap (first 64).
pub const SEED_CHANNEL_POS_CAP: usize = 64;

/// Extra reference margin past the tail length in tail windows.
const TAIL_MARGIN: u32 = 30;

/// Successful split rescue: replacement placement pieces.
#[derive(Clone, Debug)]
pub struct SplitRescue {
    /// New genomic start of the whole stitched alignment.
    pub pos: u32,
    /// Conserving CIGAR over the whole read.
    pub cigar: Vec<CigarOp>,
    /// The junction the tail was placed across.
    pub junction: RefinedJunction,
    /// Rescue score (channel A retry total; tail extension score otherwise).
    pub score: i32,
}

/// Inputs describing the primary alignment of the read being rescued.
pub struct SplitContext<'a> {
    /// Reference genome.
    pub reference: &'a Reference,
    /// Minimizer index (needed by the seed channel; `None` disables it).
    pub index: Option<&'a Index>,
    /// Junction library (empty ⇒ de-novo channels enabled).
    pub lib: &'a JunctionLib,
    /// The read in alignment orientation (anchor `qpos` space).
    pub read: &'a [u8],
    /// Primary alignment contig.
    pub contig: u32,
    /// Primary alignment strand.
    pub strand: Strand,
    /// Primary genomic start.
    pub pos: u32,
    /// Primary genomic end (M/D/N advance), exclusive.
    pub ref_end: u32,
    /// First aligned read base of the primary alignment.
    pub read_start: u32,
    /// One past the last aligned read base of the primary alignment.
    pub read_end: u32,
    /// Primary CIGAR (whole read).
    pub cigar: &'a [CigarOp],
    /// Seed parameters (seed channel).
    pub seed_params: SeedParams,
    /// Splice parameters.
    pub splice_params: SpliceParams,
    /// Extension parameters.
    pub extend_params: ExtendParams,
}

impl SplitContext<'_> {
    fn minus(&self) -> bool {
        self.strand == Strand::Minus
    }

    /// Right-tail library candidates: junction starts near `ref_end` (plus)
    /// or junction ends near `pos` (minus), zipped with support.
    fn right_lib_candidates(&self) -> Vec<(&Junction, u32)> {
        let lo = self.ref_end.saturating_sub(SPLIT_LIB_WINDOW);
        let hi = self.ref_end + SPLIT_LIB_WINDOW;
        lib_window_candidates(self.lib, self.minus(), self.contig, lo, hi)
    }

    /// Left-tail library candidates: junction ends near `pos` (plus) or
    /// junction starts near `ref_end` (minus).
    fn left_lib_candidates(&self) -> Vec<(&Junction, u32)> {
        let lo = self.pos.saturating_sub(SPLIT_LIB_WINDOW);
        let hi = self.pos + SPLIT_LIB_WINDOW;
        lib_window_candidates(self.lib, !self.minus(), self.contig, lo, hi)
    }
}

/// One `retry_tail_start` outcome.
#[derive(Clone, Debug)]
pub struct RetryOutcome {
    /// Winning tail-start offset (split point in read coordinates).
    pub dro: u32,
    /// Left-window genomic start (for `new_pos` reconstruction).
    pub left_win_lo: u32,
    /// Left segment extension (`read[..dro]`).
    pub left: Extension,
    /// Right segment extension (`read[dro..]`).
    pub right: Extension,
    /// `le.score + re.score + signal_score + lib_bonus − SPLIT_COST`.
    pub total: i32,
}

/// Retry the tail start ±5bp around `tail_start` against one junction.
///
/// Per `dro`: the left segment extends in the exon window ending at `j.start`,
/// the right segment in the window starting at `j.end`; each
/// segment must align ≥ `len × 4/5`. The best total wins (ties keep the
/// lowest `dro`).
#[allow(clippy::too_many_arguments)]
pub fn retry_tail_start(
    reference: &Reference,
    read: &[u8],
    tail_start: u32,
    j: &Junction,
    signal: SpliceSignal,
    support: u32,
    params: &SpliceParams,
    extend_params: &ExtendParams,
    buf: &mut ExtendBuffer,
) -> Option<RetryOutcome> {
    let ctg = reference.contigs.get(j.contig as usize)?;
    let read_len = read.len() as u32;
    let lo_dro = tail_start.saturating_sub(TAIL_RETRY_RADIUS as u32).max(1);
    let hi_dro = (tail_start + TAIL_RETRY_RADIUS as u32).min(read_len.saturating_sub(1));
    let signal_score = signal.score(params);
    let lib_bonus = (support.min(10) * 5) as i32;

    let mut best: Option<RetryOutcome> = None;
    for dro in lo_dro..=hi_dro {
        let lseg = &read[..dro as usize];
        let rseg = &read[dro as usize..];
        let llen = lseg.len() as u32;
        let rlen = rseg.len() as u32;

        // Oriented-read geometry is strand-symmetric: the left segment
        // ends at the donor side (j.start), the right segment starts at the
        // acceptor side (j.end).
        let l_hi = j.start.min(ctg.len);
        let l_lo = l_hi.saturating_sub(llen + TAIL_MARGIN);
        let r_lo = j.end;
        let r_hi =
            ((r_lo as u64) + (rlen as u64) + (TAIL_MARGIN as u64)).min(ctg.len as u64) as u32;
        if l_hi <= l_lo || r_hi <= r_lo {
            continue;
        }
        let lwin = ctg.slice_ascii(l_lo, l_hi);
        let rwin = ctg.slice_ascii(r_lo, r_hi);
        let hint = DiagHint {
            offset: 0,
            num: 1,
            den: 1,
        };
        let le = extend_hint(lseg, &lwin, extend_params, buf, hint);
        let re = extend_hint(rseg, &rwin, extend_params, buf, hint);
        let l_aligned = le.read_end.saturating_sub(le.read_start) as usize;
        let r_aligned = re.read_end.saturating_sub(re.read_start) as usize;
        if l_aligned * 5 < lseg.len() * 4 || r_aligned * 5 < rseg.len() * 4 {
            continue;
        }
        let total = le.score + re.score + signal_score + lib_bonus - SPLIT_COST;
        if best.as_ref().is_none_or(|b| total > b.total) {
            best = Some(RetryOutcome {
                dro,
                left_win_lo: l_lo,
                left: le,
                right: re,
                total,
            });
        }
    }
    best
}

/// Channel-A stitch: `[left ext] N(intron) [right ext]`, fully independent of
/// the primary CIGAR; `new_pos` = left segment reference start.
pub fn build_rescue_ext(outcome: &RetryOutcome, refined: &RefinedJunction) -> SplitRescue {
    let intron = refined.junction.end - refined.junction.start;
    let mut cigar = outcome.left.cigar.clone();
    push_op(&mut cigar, CigarOp::RefSkip(intron));
    for op in &outcome.right.cigar {
        push_op(&mut cigar, *op);
    }
    SplitRescue {
        pos: outcome.left_win_lo + outcome.left.ref_start,
        cigar,
        junction: refined.clone(),
        score: outcome.total,
    }
}

/// Right-tail stitch against the primary body: body CIGAR truncated at
/// `body_read_end` read depth (dropping the trailing soft clip) + `N` + the
/// tail extension as a whole. Conserves the read length exactly.
pub fn build_rescue(
    body_cigar: &[CigarOp],
    body_read_end: u32,
    tail_ext: &Extension,
    intron_len: u32,
) -> Vec<CigarOp> {
    let mut cigar = cigar_before_read_depth(body_cigar, body_read_end);
    push_op(&mut cigar, CigarOp::RefSkip(intron_len));
    for op in &tail_ext.cigar {
        push_op(&mut cigar, *op);
    }
    cigar
}

/// Left-tail stitch: the tail extension as a whole + `N` + body CIGAR from
/// `body_read_start` read depth (dropping the leading soft clip). Conserves
/// the read length exactly.
pub fn build_rescue_left(
    tail_ext: &Extension,
    body_cigar: &[CigarOp],
    body_read_start: u32,
    intron_len: u32,
) -> Vec<CigarOp> {
    let mut cigar: Vec<CigarOp> = Vec::new();
    for op in &tail_ext.cigar {
        push_op(&mut cigar, *op);
    }
    push_op(&mut cigar, CigarOp::RefSkip(intron_len));
    for op in cigar_from_read_depth(body_cigar, body_read_start) {
        push_op(&mut cigar, op);
    }
    cigar
}

/// Prefix of `cigar` consuming at most `depth` read bases (mid-op split; ops
/// past the depth, including trailing reference-only ops, are dropped).
pub fn cigar_before_read_depth(cigar: &[CigarOp], depth: u32) -> Vec<CigarOp> {
    let mut out = Vec::new();
    let mut used = 0u32;
    'outer: for &op in cigar {
        match op {
            CigarOp::Match(n) | CigarOp::Ins(n) | CigarOp::SoftClip(n) => {
                if used >= depth {
                    break 'outer;
                }
                let room = depth - used;
                if n <= room {
                    push_op(&mut out, op);
                    used += n;
                } else {
                    push_op(&mut out, same_kind(op, room));
                    used = depth;
                }
            }
            CigarOp::Del(_) | CigarOp::RefSkip(_) => {
                if used < depth {
                    push_op(&mut out, op);
                }
            }
        }
    }
    out
}

/// Suffix of `cigar` from `depth` read bases on (leading op split if needed).
pub fn cigar_from_read_depth(cigar: &[CigarOp], depth: u32) -> Vec<CigarOp> {
    let mut out = Vec::new();
    let mut used = 0u32;
    for &op in cigar {
        match op {
            CigarOp::Match(n) | CigarOp::Ins(n) | CigarOp::SoftClip(n) => {
                let end = used + n;
                if end <= depth {
                    used = end;
                    continue;
                }
                if used < depth {
                    let rest = end - depth;
                    push_op(&mut out, same_kind(op, rest));
                } else {
                    push_op(&mut out, op);
                }
                used = end;
            }
            CigarOp::Del(_) | CigarOp::RefSkip(_) => {
                if used >= depth {
                    push_op(&mut out, op);
                }
            }
        }
    }
    out
}

fn same_kind(op: CigarOp, len: u32) -> CigarOp {
    match op {
        CigarOp::Match(_) => CigarOp::Match(len),
        CigarOp::Ins(_) => CigarOp::Ins(len),
        CigarOp::Del(_) => CigarOp::Del(len),
        CigarOp::RefSkip(_) => CigarOp::RefSkip(len),
        CigarOp::SoftClip(_) => CigarOp::SoftClip(len),
    }
}

pub(crate) fn extend_tail(
    ctg: &crate::fasta::Contig,
    tail: &[u8],
    win_lo: u32,
    win_hi: u32,
    params: &ExtendParams,
    buf: &mut ExtendBuffer,
) -> Option<(u32, Extension)> {
    if win_hi <= win_lo {
        return None;
    }
    let win = ctg.slice_ascii(win_lo, win_hi);
    let ext = extend_hint(
        tail,
        &win,
        params,
        buf,
        DiagHint {
            offset: 0,
            num: 1,
            den: 1,
        },
    );
    if !tail_gate(&ext, tail.len()) {
        return None;
    }
    Some((win_lo, ext))
}

/// Channel A (library-driven): best retry over the ±20bp window junctions.
fn channel_a(ctx: &SplitContext, right: bool, buf: &mut ExtendBuffer) -> Option<SplitRescue> {
    let cands = if right {
        ctx.right_lib_candidates()
    } else {
        ctx.left_lib_candidates()
    };
    let ctg = ctx.reference.contigs.get(ctx.contig as usize)?;
    let mut best: Option<(RetryOutcome, RefinedJunction)> = None;
    for (j, support) in cands {
        if j.minus_strand != ctx.minus() {
            continue;
        }
        let intron_len = j.end - j.start;
        if intron_len < ctx.splice_params.min_intron || intron_len > ctx.splice_params.max_intron {
            continue;
        }
        let signal = junction_signal(ctg, j);
        if signal == SpliceSignal::NonCanonical && support < 2 {
            continue;
        }
        let tail_start = if right { ctx.read_end } else { ctx.read_start };
        let outcome = match retry_tail_start(
            ctx.reference,
            ctx.read,
            tail_start,
            j,
            signal,
            support,
            &ctx.splice_params,
            &ctx.extend_params,
            buf,
        ) {
            Some(o) => o,
            None => continue,
        };
        // Keep the highest total; ties keep the first candidate (sorted).
        if best.as_ref().is_none_or(|(b, _)| outcome.total > b.total) {
            best = Some((
                outcome,
                RefinedJunction {
                    junction: *j,
                    signal,
                    known_support: support,
                },
            ));
        }
    }
    let (outcome, refined) = best?;
    let rescue = build_rescue_ext(&outcome, &refined);
    if !cigar_conserving(&rescue.cigar, ctx.read.len()) {
        return None;
    }
    Some(rescue)
}

/// Dinucleotide ending at `p` (i.e. `[p-2, p)`), read from a decoded window
/// starting at `lo`; `NN` when out of range.
const GT_SET: [[u8; 2]; 3] = [*b"GT", *b"GC", *b"AT"];
const AG_AC_SET: [[u8; 2]; 2] = [*b"AG", *b"AC"];
const CT_GT_SET: [[u8; 2]; 2] = [*b"CT", *b"GT"];
const AC_GC_AT_SET: [[u8; 2]; 3] = [*b"AC", *b"GC", *b"AT"];

fn dinuc_before(win: &[u8], lo: u32, p: u32) -> [u8; 2] {
    if p < lo + 2 {
        return *b"NN";
    }
    let i = (p - lo) as usize;
    if i > win.len() {
        return *b"NN";
    }
    [win[i - 2].to_ascii_uppercase(), win[i - 1].to_ascii_uppercase()]
}

fn dinuc_after(win: &[u8], lo: u32, p: u32) -> [u8; 2] {
    let i = (p - lo) as usize;
    if i + 2 > win.len() {
        return *b"NN";
    }
    [win[i].to_ascii_uppercase(), win[i + 1].to_ascii_uppercase()]
}

fn matches_at(win: &[u8], lo: u32, p: u32, pat: &[u8]) -> bool {
    let i = (p - lo) as usize;
    if p < lo || i + pat.len() > win.len() {
        return false;
    }
    win[i..i + pat.len()]
        .iter()
        .zip(pat.iter())
        .all(|(r, q)| r.eq_ignore_ascii_case(q))
}

/// Channel B (direct scan, empty library only).
///
/// Right tail: donors ±12 around `ref_end`, probe scan ≤50kb downstream for
/// an exact tail-prefix match ending at a plausible acceptor. Left tail:
/// acceptors ±12 around `pos`, probe scan ≤50kb upstream for an exact
/// tail-suffix match at a plausible donor. Candidate score =
/// `sig.score*100 - |body_boundary_deviation| - |intron - 5000|/1000`; the
/// best candidate's whole tail is then extension-validated (≥ 4/5, score > 0).
/// The discovered junction is recorded with the detected splice signal, support 0.
fn channel_b(ctx: &SplitContext, right: bool, buf: &mut ExtendBuffer) -> Option<SplitRescue> {
    let ctg = ctx.reference.contigs.get(ctx.contig as usize)?;
    let tl = if right {
        ctx.read.len() - ctx.read_end as usize
    } else {
        ctx.read_start as usize
    };
    if !(DIRECT_TAIL_MIN..DIRECT_TAIL_MAX).contains(&tl) {
        return None;
    }
    let probe_len = tl.min(PROBE_MAX);
    let min_intron = ctx.splice_params.min_intron;
    let span = DIRECT_SCAN_SPAN as u64;

    // Transcript-oriented dinucleotide sets: body-side boundary vs tail-side
    // boundary of the intron (mirrored on minus-strand junctions).
    let (body_set, tail_set): (&[[u8; 2]], &[[u8; 2]]) = match (right, ctx.minus()) {
        (true, false) => (&GT_SET, &AG_AC_SET),
        (true, true) => (&CT_GT_SET, &AC_GC_AT_SET),
        (false, false) => (&AG_AC_SET, &GT_SET),
        (false, true) => (&AC_GC_AT_SET, &CT_GT_SET),
    };

    // One decoded window covering the ±12 body scan and the full search span.
    let (scan_lo, scan_hi) = if right {
        let lo = ctx.ref_end.saturating_sub(DIRECT_DONOR_RADIUS as u32);
        let hi = (ctx.ref_end as u64 + DIRECT_DONOR_RADIUS as u64).min(ctg.len as u64) as u32;
        (lo, ((hi as u64) + span + 2).min(ctg.len as u64) as u32)
    } else {
        let hi = (ctx.pos as u64 + DIRECT_DONOR_RADIUS as u64 + 1).min(ctg.len as u64) as u32;
        let lo = ctx
            .pos
            .saturating_sub(DIRECT_DONOR_RADIUS as u32)
            .saturating_sub(span as u32 + 2);
        (lo, hi)
    };
    if scan_hi <= scan_lo || (scan_hi - scan_lo) as u64 <= min_intron as u64 {
        return None;
    }
    let win = ctg.slice_ascii(scan_lo, scan_hi);
    let body_anchor = if right { ctx.ref_end } else { ctx.pos };

    // Body-overshoot note: a continuous extension can run a few bp INTO the
    // intron (repeat similarity). The true junction boundary then sits `over`
    // bp before the body edge, and the exon fragment on the tail side starts
    // `over` bp earlier in the read. Each boundary candidate computes its own
    // overshoot and uses an overshoot-adjusted probe and tail.
    let mut best: Option<(i32, u32, u32, u32)> = None; // (score, intron_start, intron_end, over)
    if right {
        let d_lo = ctx.ref_end.saturating_sub(DIRECT_DONOR_RADIUS as u32);
        let d_hi = (ctx.ref_end + DIRECT_DONOR_RADIUS as u32).min(ctg.len.saturating_sub(2));
        for d in d_lo..=d_hi {
            if !body_set.contains(&dinuc_after(&win, scan_lo, d)) {
                continue;
            }
            let dev = (d as i64 - body_anchor as i64).abs() as i32;
            let over = (ctx.ref_end as i64 - d as i64).max(0) as u32;
            if over > ctx.read_end {
                continue;
            }
            let probe_start = (ctx.read_end - over) as usize;
            let probe_d = &ctx.read[probe_start..probe_start + probe_len];
            let q_lo = d as u64 + min_intron as u64;
            let q_hi = (d as u64 + span).min(scan_hi as u64);
            for q in q_lo..=q_hi {
                let qu = q as u32;
                if tail_set.contains(&dinuc_before(&win, scan_lo, qu))
                    && matches_at(&win, scan_lo, qu, probe_d)
                {
                    let intron = qu as i64 - d as i64;
                    let score = 4 * 100 - dev - ((intron - 5000).abs() / 1000) as i32;
                    if best.as_ref().is_none_or(|(s, _, _, _)| score > *s) {
                        best = Some((score, d, qu, over));
                    }
                }
            }
        }
    } else {
        let a_lo = ctx.pos.saturating_sub(DIRECT_DONOR_RADIUS as u32);
        let a_hi = (ctx.pos + DIRECT_DONOR_RADIUS as u32).min(ctg.len);
        for a in a_lo..=a_hi {
            if a < min_intron + 2 || !body_set.contains(&dinuc_before(&win, scan_lo, a)) {
                continue;
            }
            let dev = (a as i64 - body_anchor as i64).abs() as i32;
            let over = (a as i64 - ctx.pos as i64).max(0) as u32;
            if ctx.read_start + over > ctx.read.len() as u32 {
                continue;
            }
            let probe_end = (ctx.read_start + over) as usize;
            let probe_d = &ctx.read[probe_end - probe_len..probe_end];
            let s_lo = ((a as i64 - DIRECT_SCAN_SPAN as i64).max(scan_lo as i64).max(0)) as u64;
            let s_hi = a as u64 - min_intron as u64;
            for s in s_lo..=s_hi {
                if s < probe_len as u64 {
                    continue;
                }
                let su = s as u32;
                if tail_set.contains(&dinuc_after(&win, scan_lo, su))
                    && matches_at(&win, scan_lo, su - probe_len as u32, probe_d)
                {
                    let intron = a as i64 - s as i64;
                    let score = 4 * 100 - dev - ((intron - 5000).abs() / 1000) as i32;
                    if best.as_ref().is_none_or(|(s0, _, _, _)| score > *s0) {
                        best = Some((score, su, a, over));
                    }
                }
            }
        }
    }
    if std::env::var_os("ESP_PROBE").is_some() {
        eprintln!("[probe-cb] right={right} tail={tl} best={best:?}");
    }
    let (scan_score, start, end, over) = best?;
    let intron_len = end - start;

    // Overshoot-adjusted tail and body edge: the body ran `over` bp into the
    // intron, so the true exon fragment is `over` bp longer on the tail side.
    let (ext_tail, body_edge) = if right {
        let edge = ctx.read_end - over;
        (&ctx.read[edge as usize..], edge)
    } else {
        let edge = ctx.read_start + over;
        (&ctx.read[..edge as usize], edge)
    };
    let etl = ext_tail.len();

    // Whole-tail extension validation.
    let (win_lo, win_hi) = if right {
        let hi = ((end as u64) + (etl as u64) + (TAIL_MARGIN as u64)).min(ctg.len as u64) as u32;
        (end, hi)
    } else {
        let hi = start;
        let lo = hi.saturating_sub(etl as u32 + TAIL_MARGIN);
        (lo, hi)
    };
    let ext_res = extend_tail(ctg, ext_tail, win_lo, win_hi, &ctx.extend_params, buf);
    if std::env::var_os("ESP_PROBE").is_some() {
        eprintln!("[probe-cb] start={start} end={end} over={over} ext={}", ext_res.as_ref().map_or("None".to_string(), |(_, e)| format!("score={} ref_start={}", e.score, e.ref_start)));
    }
    let (_, ext) = ext_res?;
    let signal = splice_signal(
        dinuc_after(&win, scan_lo, start),
        dinuc_before(&win, scan_lo, end),
        ctx.minus(),
    );
    let junction = RefinedJunction {
        junction: Junction {
            contig: ctx.contig,
            start,
            end,
            minus_strand: ctx.minus(),
        },
        signal,
        known_support: 0,
    };
    let (cigar, pos) = if right {
        (
            build_rescue(ctx.cigar, body_edge, &ext, intron_len),
            ctx.pos,
        )
    } else {
        (
            build_rescue_left(&ext, ctx.cigar, body_edge, intron_len),
            win_lo + ext.ref_start,
        )
    };
    if !cigar_conserving(&cigar, ctx.read.len()) {
        return None;
    }
    Some(SplitRescue {
        pos,
        cigar,
        junction,
        score: ext.score + scan_score,
    })
}

/// Channel C (seed-driven, empty library only): seed the tail, map its first
/// minimizers through the index, keep hits on the same contig/strand falling
/// in the min/max intron window across the body edge, refine with
/// `refine_junction` + S14, extension-validate. First success (deterministic
/// candidate order) wins.
fn channel_c(ctx: &SplitContext, right: bool, buf: &mut ExtendBuffer) -> Option<SplitRescue> {
    let index = ctx.index?;
    let ctg = ctx.reference.contigs.get(ctx.contig as usize)?;
    let (tail, body_edge_read) = if right {
        (&ctx.read[ctx.read_end as usize..], ctx.read_end)
    } else {
        (&ctx.read[..ctx.read_start as usize], ctx.read_start)
    };
    let tl = tail.len();
    if tl < SEED_TAIL_MIN || tl < ctx.seed_params.k as usize {
        return None;
    }
    let k = ctx.seed_params.k;
    let mins: Vec<_> = minimizers(tail, ctx.seed_params)
        .into_iter()
        .take(SEED_CHANNEL_MINS)
        .collect();

    // Naive (start, end) intron candidates from seed placements.
    let mut cands: Vec<(u32, u32)> = Vec::new();
    for m in &mins {
        let hit = index.query(m.kmer);
        if hit.count == 0 || hit.count > SEED_CHANNEL_OCC_CAP {
            continue;
        }
        for &p in hit.positions.iter().take(SEED_CHANNEL_POS_CAP) {
            let (contig2, rpos, rstrand) = crate::index::unpack_pos(p);
            if contig2 != ctx.contig {
                continue;
            }
            let astrand = Strand::xor(m.strand, rstrand);
            if astrand != ctx.strand {
                continue;
            }
            let qpos = if astrand == Strand::Plus {
                m.pos
            } else {
                tl as u32 - m.pos - k
            };
            let g = rpos as i64 - qpos as i64; // genomic start of the tail
            let cand = if right {
                (ctx.ref_end as i64, g)
            } else {
                (g + tl as i64, ctx.pos as i64)
            };
            if cand.1 <= cand.0 {
                continue;
            }
            let len = cand.1 - cand.0;
            if len < ctx.splice_params.min_intron as i64
                || len > ctx.splice_params.max_intron as i64
            {
                continue;
            }
            cands.push((cand.0 as u32, cand.1 as u32));
        }
    }
    cands.sort_unstable();
    cands.dedup();

    for (s0, e0) in cands {
        let Some(refined) = refine_junction(
            ctx.reference,
            ctx.contig,
            s0,
            e0,
            ctx.minus(),
            ctx.lib,
            &ctx.splice_params,
        ) else {
            continue;
        };
        let intron_len = refined.junction.end - refined.junction.start;
        if !s14_pass(intron_len, refined.signal, refined.known_support) {
            continue;
        }
        let (win_lo, win_hi) = if right {
            let lo = refined.junction.end;
            let hi = ((lo as u64) + (tl as u64) + (TAIL_MARGIN as u64)).min(ctg.len as u64) as u32;
            (lo, hi)
        } else {
            let hi = refined.junction.start;
            let lo = hi.saturating_sub(tl as u32 + TAIL_MARGIN);
            (lo, hi)
        };
        let Some((win_lo, ext)) = extend_tail(ctg, tail, win_lo, win_hi, &ctx.extend_params, buf)
        else {
            continue;
        };
        let (cigar, pos) = if right {
            (
                build_rescue(ctx.cigar, body_edge_read, &ext, intron_len),
                ctx.pos,
            )
        } else {
            (
                build_rescue_left(&ext, ctx.cigar, body_edge_read, intron_len),
                win_lo + ext.ref_start,
            )
        };
        if !cigar_conserving(&cigar, ctx.read.len()) {
            continue;
        }
        return Some(SplitRescue {
            pos,
            cigar,
            junction: refined,
            score: ext.score,
        });
    }
    None
}

/// Rescue a 3' (right) soft-clipped tail across a splice junction.
///
/// Library non-empty ⇒ channel A only; otherwise direct scan (B) then tail
/// seeding (C). Tails shorter than [`TAIL_MIN`] are not rescued.
pub fn rescue_right_tail(ctx: &SplitContext, buf: &mut ExtendBuffer) -> Option<SplitRescue> {
    let tail_len = ctx.read.len() - ctx.read_end as usize;
    if tail_len < TAIL_MIN {
        return None;
    }
    if std::env::var_os("ESP_PROBE").is_some() && !ctx.lib.is_empty() {
        eprintln!("[probe-tail] rescue_right_tail tail_len={tail_len} read_end={} lib_empty={}", ctx.read_end, ctx.lib.is_empty());
    }
    if !ctx.lib.is_empty() {
        // Library first (high-confidence); fall back to de-novo scan/seed
        // (B/C) for junctions NOT in the library (novel splices).
        if let Some(a) = channel_a(ctx, true, buf) {
            return Some(a);
        }
    }
    channel_b(ctx, true, buf).or_else(|| channel_c(ctx, true, buf))
}

/// Rescue a 5' (left) soft-clipped tail across a splice junction (mirror of
/// [`rescue_right_tail`]; the alignment start moves left on success).
pub fn rescue_left_tail(ctx: &SplitContext, buf: &mut ExtendBuffer) -> Option<SplitRescue> {
    let tail_len = ctx.read_start as usize;
    if tail_len < TAIL_MIN {
        return None;
    }
    if !ctx.lib.is_empty() {
        if let Some(a) = channel_a(ctx, false, buf) {
            return Some(a);
        }
    }
    channel_b(ctx, false, buf).or_else(|| channel_c(ctx, false, buf))
}
