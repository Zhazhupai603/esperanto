//! Splice signals and pseudo-reference spliced alignment.
//!
//! Signal classification is orientation-aware: a plus-strand intron reads
//! GT..AG / GC..AG / AT..AC on the genome; a minus-strand intron is tested
//! through its plus-strand mirror (CT..AC etc.). `refine_junction` walks a
//! ±`refine_radius` box around naive breakpoints, `align_spliced` concatenates
//! the flanking exon windows into a pseudo-reference, extends with a direct
//! diagonal hint and stitches the CIGAR back to genomic coordinates with a
//! `RefSkip` at the seam. S14 guards (short/unmotivated introns without
//! library support) gate every de-novo junction.

use crate::extend::{
    extend_hint, push_op, CigarOp, DiagHint, ExtendBuffer, ExtendParams, Extension,
};
use crate::fasta::{Base, Reference};
use crate::gtf::{Junction, JunctionLib, RefinedJunction, SpliceSignal};
use crate::seed::Strand;

/// Splice scoring parameters (frozen; see parameter table).
#[derive(Clone, Copy, Debug)]
pub struct SpliceParams {
    /// Shortest intron accepted anywhere.
    pub min_intron: u32,
    /// Longest intron accepted anywhere.
    pub max_intron: u32,
    /// GT-AG bonus.
    pub gt_ag_bonus: i32,
    /// GC-AG bonus.
    pub gc_ag_bonus: i32,
    /// AT-AC bonus.
    pub at_ac_bonus: i32,
    /// NonCanonical penalty (stored positive; subtracted).
    pub noncanonical_penalty: i32,
    /// Bonus for a junction present in the library.
    pub known_bonus: i32,
    /// Breakpoint refinement search radius.
    pub refine_radius: u32,
}

impl Default for SpliceParams {
    fn default() -> SpliceParams {
        SpliceParams {
            min_intron: 20,
            max_intron: 1_000_000,
            gt_ag_bonus: 4,
            gc_ag_bonus: 2,
            at_ac_bonus: 1,
            noncanonical_penalty: 8,
            known_bonus: 6,
            refine_radius: 12,
        }
    }
}

/// Cost of rewriting an alignment into a split (frozen).
pub const SPLIT_COST: i32 = 4;

/// S14 guard: introns shorter than this need library support.
pub const S14_MIN_INTRON_NO_LIB: u32 = 50;

impl SpliceSignal {
    /// Signal score under `params` (bonuses positive, NonCanonical negative).
    pub fn score(self, p: &SpliceParams) -> i32 {
        match self {
            SpliceSignal::GtAg => p.gt_ag_bonus,
            SpliceSignal::GcAg => p.gc_ag_bonus,
            SpliceSignal::AtAc => p.at_ac_bonus,
            SpliceSignal::NonCanonical => -p.noncanonical_penalty,
        }
    }
}

/// Classify donor/acceptor dinucleotides (plus-strand genomic reading).
///
/// `minus` selects the mirror tests: a minus-strand GT..AG intron appears on
/// the plus strand as CT..AC (and analogously GC..AG → CT..GC, AT..AC →
/// GT..AT). Comparison is case-insensitive; anything else is NonCanonical.
pub fn splice_signal(donor: [u8; 2], acceptor: [u8; 2], minus: bool) -> SpliceSignal {
    let d = [donor[0].to_ascii_uppercase(), donor[1].to_ascii_uppercase()];
    let a = [
        acceptor[0].to_ascii_uppercase(),
        acceptor[1].to_ascii_uppercase(),
    ];
    if !minus {
        if d == *b"GT" && a == *b"AG" {
            SpliceSignal::GtAg
        } else if d == *b"GC" && a == *b"AG" {
            SpliceSignal::GcAg
        } else if d == *b"AT" && a == *b"AC" {
            SpliceSignal::AtAc
        } else {
            SpliceSignal::NonCanonical
        }
    } else if d == *b"CT" && a == *b"AC" {
        SpliceSignal::GtAg
    } else if d == *b"CT" && a == *b"GC" {
        SpliceSignal::GcAg
    } else if d == *b"GT" && a == *b"AT" {
        SpliceSignal::AtAc
    } else {
        SpliceSignal::NonCanonical
    }
}

/// S14 guard: reject library-unsupported introns that are too short or
/// unmotivated (NonCanonical with support < 2).
pub fn s14_pass(intron_len: u32, signal: SpliceSignal, support: u32) -> bool {
    if intron_len < S14_MIN_INTRON_NO_LIB && support == 0 {
        return false;
    }
    if signal == SpliceSignal::NonCanonical && support == 0 {
        return false;
    }
    true
}

/// Intron-length penalty, frozen formula: `8 × |log10(len) − 3|` as `i32` —
/// linear in log10 distance from the 1 kb mode (zero at 1 kb).
pub fn intron_length_penalty(len: u32) -> i32 {
    (8.0 * ((len.max(1) as f64).log10() - 3.0).abs()) as i32
}

/// Refine naive breakpoints against splice signals and the junction library.
///
/// Enumerates `(start, end)` in the ±`refine_radius` box around the naive
/// coordinates; score = `signal.score × 100 − dist + known_bonus` with
/// `dist = |Δstart| + |Δend|` and the bonus applied when the exact candidate
/// junction is in the library. Best (strictly greater) score wins; iteration
/// is in ascending offset order, so ties resolve deterministically to the
/// lowest offsets. Returns `None` when no in-bounds candidate exists.
pub fn refine_junction(
    reference: &Reference,
    contig: u32,
    naive_start: u32,
    naive_end: u32,
    minus: bool,
    lib: &JunctionLib,
    params: &SpliceParams,
) -> Option<RefinedJunction> {
    let ctg = reference.contigs.get(contig as usize)?;
    let len = ctg.len;
    let r = params.refine_radius as i64;

    let lo_start = naive_start as i64 - r;
    let hi_start = naive_start as i64 + r;
    let lo_end = naive_end as i64 - r;
    let hi_end = naive_end as i64 + r;

    // Pre-extract the dinucleotides once per distinct offset.
    let mut donors: Vec<[u8; 2]> = Vec::with_capacity((hi_start - lo_start + 1) as usize);
    for s in lo_start..=hi_start {
        donors.push(dinuc(ctg, s as u32, true));
    }
    let mut acceptors: Vec<[u8; 2]> = Vec::with_capacity((hi_end - lo_end + 1) as usize);
    for e in lo_end..=hi_end {
        acceptors.push(dinuc(ctg, e as u32, false));
    }

    let mut best: Option<(i32, RefinedJunction)> = None;
    for (si, s) in (lo_start..=hi_start).enumerate() {
        if s < 0 || s as u64 + 2 > len as u64 {
            continue;
        }
        let ds = (s - naive_start as i64).abs() as i32;
        for (ei, e) in (lo_end..=hi_end).enumerate() {
            if e < 2 || e as u64 > len as u64 || e <= s {
                continue;
            }
            let de = (e - naive_end as i64).abs() as i32;
            let signal = splice_signal(donors[si], acceptors[ei], minus);
            let junction = Junction {
                contig,
                start: s as u32,
                end: e as u32,
                minus_strand: minus,
            };
            let support = lib.support(&junction);
            let mut score = signal.score(params) * 100 - ds - de;
            if support > 0 {
                score += params.known_bonus;
            }
            if best.as_ref().is_none_or(|(bs, _)| score > *bs) {
                best = Some((
                    score,
                    RefinedJunction {
                        junction,
                        signal,
                        known_support: support,
                    },
                ));
            }
        }
    }
    best.map(|(_, rj)| rj)
}

/// Read-aware junction refinement (preferred over [`refine_junction`]).
///
/// The intron length is fixed by the seed anchors: `ref_gap − read_gap`,
/// where `ref_gap = naive_end − naive_start` and `read_gap = q_b − q_a`. The
/// boundary therefore has a single degree of freedom — where along the
/// un-anchored read gap the splice falls. Enumerate that shift and pick the
/// boundary that carries a canonical splice signal AND matches the read on
/// both flanks (editing-aware). This removes the proximity bias that drifted
/// off at competing GT/AG motifs.
#[allow(clippy::too_many_arguments)]
pub fn refine_junction_read(
    reference: &Reference,
    read: &[u8],
    contig: u32,
    naive_start: u32,
    naive_end: u32,
    q_a: u32,
    q_b: u32,
    minus: bool,
    lib: &JunctionLib,
    params: &SpliceParams,
) -> Option<RefinedJunction> {
    let ctg = reference.contigs.get(contig as usize)?;
    let len = ctg.len;
    let read_len = read.len() as u32;

    let ref_gap = naive_end as i64 - naive_start as i64;
    let read_gap = q_b as i64 - q_a as i64;
    if read_gap > 20 {
        return None; // not a tight junction
    }
    let intron_len = ref_gap - read_gap;
    if intron_len < params.min_intron as i64 || intron_len > params.max_intron as i64 {
        return None;
    }

    // Boundary shift delta = s − naive_start, constrained to the read gap
    // (plus a small margin for anchor-boundary imprecision). The read gap may
    // be negative when the two anchors overlap by a base at the junction (a
    // read base matching both the last exon base and the acceptor dinucleotide);
    // signed arithmetic keeps the intron length exact.
    let margin = 2i64;
    let lo = read_gap.min(0) - margin;
    let hi = read_gap.max(0) + margin;

    let mut best: Option<(i32, RefinedJunction)> = None;
    for delta in lo..=hi {
        let s = naive_start as i64 + delta;
        if s < 0 {
            continue;
        }
        let e = s + intron_len;
        if e < 2 || e as u64 > len as u64 {
            continue;
        }
        let qb = q_a as i64 + delta;
        if qb < 1 || qb as u32 >= read_len {
            continue;
        }
        let signal = splice_signal(
            dinuc(ctg, s as u32, true),
            dinuc(ctg, e as u32, false),
            minus,
        );
        let junction = Junction {
            contig,
            start: s as u32,
            end: e as u32,
            minus_strand: minus,
        };
        let support = lib.support(&junction);

        // Read-match on both flanks of the boundary (editing-aware).
        let mut m = 0i32;
        let w = 4i64;
        for i in 1..=w {
            let q1 = qb - i;
            let g1 = s - i;
            if q1 >= 0
                && g1 >= 0
                && edit_equiv(Base::from_ascii(read[q1 as usize]), ctg.base(g1 as u32))
            {
                m += 1;
            }
            let q2 = qb + i - 1;
            let g2 = e + i - 1;
            if (q2 as u32) < read_len
                && (g2 as u64) < len as u64
                && edit_equiv(Base::from_ascii(read[q2 as usize]), ctg.base(g2 as u32))
            {
                m += 1;
            }
        }

        let mut score = m * 50 + signal.score(params) * 5 - (delta.abs() as i32) * 5;
        if support > 0 {
            score += params.known_bonus;
        }
        if best.as_ref().is_none_or(|(bs, _)| score > *bs) {
            best = Some((
                score,
                RefinedJunction {
                    junction,
                    signal,
                    known_support: support,
                },
            ));
        }
    }
    best.map(|(_, rj)| rj)
}

/// Editing-equivalent bases: exact match, or A↔G / T↔C (A-to-I editing on
/// either strand).
fn edit_equiv(a: Base, b: Base) -> bool {
    use Base::{A, C, G, T};
    a == b || matches!((a, b), (A, G) | (G, A) | (T, C) | (C, T))
}

/// Read a donor (at `p`) or acceptor (at `p − 2`) dinucleotide, N-padded when
/// out of bounds so it can never match a signal.
fn dinuc(ctg: &crate::fasta::Contig, p: u32, donor: bool) -> [u8; 2] {
    if !donor && p < 2 {
        return *b"NN";
    }
    let range = if donor { (p, p + 2) } else { (p - 2, p) };
    let (a, b) = range;
    if a as u64 >= ctg.len as u64 {
        return *b"NN";
    }
    let b = b.min(ctg.len);
    let s = ctg.slice_ascii(a, b);
    let mut out = *b"NN";
    for (i, &c) in s.iter().enumerate() {
        out[i] = c;
    }
    out
}

/// Flank (in reference bases) added around the outermost segments of the
/// pseudo-reference (frozen: 30).
pub const ALIGN_SPLICED_FLANK: u32 = 30;

/// One collinear anchor run of a chain: a genomic span with the matching
/// read span (alignment orientation).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ChainSegment {
    /// First read base of the segment (inclusive).
    pub q_start: u32,
    /// One past the last read base of the segment.
    pub q_end: u32,
    /// First reference base of the segment (inclusive).
    pub r_start: u32,
    /// One past the last reference base of the segment.
    pub r_end: u32,
}

/// Split a chain's anchors (sorted by `(rpos, qpos)`) into segments: a
/// break lands where the reference hole exceeds `min_intron` AND the
/// read-side gap is small (≤ 10 — anchors read-adjacent across the hole);
/// otherwise the anchor extends the current segment.
pub fn segment_chain(chain: &crate::chain::Chain, k: u32, min_intron: u32) -> Vec<ChainSegment> {
    let mut segs: Vec<ChainSegment> = Vec::new();
    let anchors = &chain.anchors;
    if anchors.is_empty() {
        return segs;
    }
    let mut cur = ChainSegment {
        q_start: anchors[0].qpos,
        q_end: anchors[0].qpos + k,
        r_start: anchors[0].rpos,
        r_end: anchors[0].rpos + k,
    };
    for w in anchors.windows(2) {
        let (a1, a2) = (w[0], w[1]);
        let ref_gap = a2.rpos.saturating_sub(a1.rpos + k);
        let read_gap = a2.qpos.saturating_sub(a1.qpos + k);
        if ref_gap > min_intron && read_gap <= 10 {
            segs.push(cur);
            cur = ChainSegment {
                q_start: a2.qpos,
                q_end: a2.qpos + k,
                r_start: a2.rpos,
                r_end: a2.rpos + k,
            };
        } else {
            cur.q_end = a2.qpos + k;
            cur.r_end = a2.rpos + k;
        }
    }
    segs.push(cur);
    segs
}

/// A stitched multi-segment spliced alignment.
#[derive(Clone, Debug)]
pub struct SplicedAlignment {
    /// The raw extension against the pseudo-reference.
    pub extension: Extension,
    /// Refined junctions, one per adjacent segment pair.
    pub junctions: Vec<RefinedJunction>,
    /// CIGAR mapped back to genomic coordinates (`RefSkip` at each seam).
    pub cigar: Vec<CigarOp>,
    /// Per-segment `(pseudo_start, ref_start, ref_end)` bounds of the
    /// pseudo-reference slices.
    pub exon_bounds: Vec<(u32, u32, u32)>,
}

impl SplicedAlignment {
    /// Genomic start of the stitched alignment: the pseudo coordinate
    /// `extension.ref_start` mapped back through `exon_bounds`.
    pub fn genomic_start(&self) -> u32 {
        let p = self.extension.ref_start;
        for &(ps, rs, re) in &self.exon_bounds {
            if p >= ps && p < ps + (re - rs) {
                return rs + (p - ps);
            }
        }
        // Defensive: past the last bound — map off the last segment.
        let &(ps, rs, _re) = self.exon_bounds.last().expect("exon_bounds nonempty");
        rs + p.saturating_sub(ps)
    }
}

/// Map a pseudo-reference CIGAR back to genomic coordinates, inserting a
/// `RefSkip` of the real intron length at every exon boundary crossing.
///
/// `p_start` is the pseudo coordinate the CIGAR's aligned block starts at
/// (`extension.ref_start`). `Match`/`Del` ops split at exon boundaries
/// (`take = min(left, remain_in_exon)`); `Ins`/`SoftClip` pass through (the
/// leading soft clip is preserved); a walk past the last exon bound emits the
/// remainder as-is (defensive — shouldn't happen).
pub fn stitch_cigar(
    cigar: &[CigarOp],
    p_start: u32,
    exon_bounds: &[(u32, u32, u32)],
    junctions: &[RefinedJunction],
) -> Vec<CigarOp> {
    let mut out: Vec<CigarOp> = Vec::with_capacity(cigar.len() + exon_bounds.len());
    let mut p = p_start as u64;
    let exon_end = |ei: usize| -> u64 {
        let (ps, rs, re) = exon_bounds[ei];
        ps as u64 + (re - rs) as u64
    };
    if exon_bounds.is_empty() {
        return cigar.to_vec();
    }
    let mut ei = 0usize;
    while ei + 1 < exon_bounds.len() && p >= exon_end(ei) {
        ei += 1;
    }
    for &op in cigar {
        match op {
            CigarOp::Ins(_) | CigarOp::SoftClip(_) | CigarOp::RefSkip(_) => push_op(&mut out, op),
            CigarOp::Match(_) | CigarOp::Del(_) => {
                let kind = |n: u32| same_op(op, n);
                let mut remaining = op_len(op) as u64;
                while remaining > 0 {
                    if ei >= exon_bounds.len() {
                        // Past the last exon bound: emit the remainder.
                        push_op(&mut out, kind(remaining as u32));
                        p += remaining;
                        break;
                    }
                    if p < exon_bounds[ei].0 as u64 {
                        p = exon_bounds[ei].0 as u64; // snap into the exon slice
                    }
                    let remain_in_exon = exon_end(ei).saturating_sub(p);
                    let take = remaining.min(remain_in_exon);
                    // Zero-take ops ARE emitted: a Match/Del resuming exactly
                    // at an exon boundary yields a 0-length op before the
                    // RefSkip (old stitch parity).
                    push_op(&mut out, kind(take as u32));
                    p += take;
                    remaining -= take;
                    if remaining > 0 {
                        if ei + 1 >= exon_bounds.len() {
                            push_op(&mut out, kind(remaining as u32));
                            p += remaining;
                            break;
                        }
                        if ei < junctions.len() {
                            let ij = junctions[ei].junction;
                            push_op(&mut out, CigarOp::RefSkip(ij.end - ij.start));
                        }
                        ei += 1;
                    }
                }
            }
        }
    }
    out
}

fn op_len(op: CigarOp) -> u32 {
    match op {
        CigarOp::Match(n)
        | CigarOp::Ins(n)
        | CigarOp::Del(n)
        | CigarOp::RefSkip(n)
        | CigarOp::SoftClip(n) => n,
    }
}

fn same_op(op: CigarOp, len: u32) -> CigarOp {
    match op {
        CigarOp::Match(_) => CigarOp::Match(len),
        CigarOp::Ins(_) => CigarOp::Ins(len),
        CigarOp::Del(_) => CigarOp::Del(len),
        CigarOp::RefSkip(_) => CigarOp::RefSkip(len),
        CigarOp::SoftClip(_) => CigarOp::SoftClip(len),
    }
}

/// Multi-segment spliced alignment of `read` across a chain's segments
/// (pseudo-reference concatenation).
///
/// The read must be in alignment orientation (anchor `qpos` space).
/// `segment_chain(chain, k, min_intron)` must yield ≥2 segments. Every
/// adjacent segment pair is refined with [`refine_junction`] and gated by
/// min/max intron plus the per-pair S14 guards (short or NonCanonical intron
/// without library support ⇒ `None`). The pseudo-reference concatenates one
/// reference slice per segment — segment 0 starting `ALIGN_SPLICED_FLANK`
/// before its first anchor, the last ending `ALIGN_SPLICED_FLANK` after its
/// last, middle boundaries at the refined junction coordinates — and the
/// read is extended over it with a direct diagonal hint. Failure anywhere
/// returns `None`.
#[allow(clippy::too_many_arguments)]
pub fn align_spliced(
    reference: &Reference,
    lib: &JunctionLib,
    read: &[u8],
    chain: &crate::chain::Chain,
    k: u32,
    params: &SpliceParams,
    extend_params: &ExtendParams,
    buf: &mut ExtendBuffer,
) -> Option<SplicedAlignment> {
    let segments = segment_chain(chain, k, params.min_intron);
    if segments.len() < 2 {
        return None;
    }
    let contig = chain.contig;
    let minus = chain.strand == Strand::Minus;
    let ctg = reference.contigs.get(contig as usize)?;

    // Per-pair refinement + S14 guards.
    let mut junctions: Vec<RefinedJunction> = Vec::with_capacity(segments.len() - 1);
    for w in segments.windows(2) {
        let refined = refine_junction_read(
            reference,
            read,
            contig,
            w[0].r_end,
            w[1].r_start,
            w[0].q_end,
            w[1].q_start,
            minus,
            lib,
            params,
        );
        if std::env::var_os("ESP_PROBE").is_some() {
            eprintln!("[probe-sp] refine seg r=({},{}) q=({},{}) -> {}",
                w[0].r_end, w[1].r_start, w[0].q_end, w[1].q_start,
                refined.as_ref().map_or("None".to_string(), |r| format!("intron=({},{}) sig={:?} known={}", r.junction.start, r.junction.end, r.signal, r.known_support)));
        }
        let refined = refined?;
        let intron_len = refined.junction.end - refined.junction.start;
        if intron_len < params.min_intron || intron_len > params.max_intron {
            return None;
        }
        if !s14_pass(intron_len, refined.signal, refined.known_support) {
            return None;
        }
        junctions.push(refined);
    }

    // Pseudo-reference: one slice per segment, middle boundaries at the
    // refined junction coordinates.
    let mut pseudo: Vec<u8> = Vec::new();
    let mut exon_bounds: Vec<(u32, u32, u32)> = Vec::with_capacity(segments.len());
    let last = segments.len() - 1;
    for (i, seg) in segments.iter().enumerate() {
        let rs = if i == 0 {
            seg.r_start.saturating_sub(ALIGN_SPLICED_FLANK)
        } else {
            junctions[i - 1].junction.end
        };
        let re = if i == last {
            ((seg.r_end as u64) + (ALIGN_SPLICED_FLANK as u64)).min(ctg.len as u64) as u32
        } else {
            junctions[i].junction.start
        };
        if re <= rs {
            return None;
        }
        exon_bounds.push((pseudo.len() as u32, rs, re));
        ctg.decode_append(rs, re, &mut pseudo);
    }

    // TEMP parity probe: ESP_PROBE_PSEUDO=/path dumps pseudo-reference + bounds.
    if let Ok(p) = std::env::var("ESP_PROBE_PSEUDO") {
        for &(ps, rs, re) in &exon_bounds {
            eprintln!("[probe] bound {ps} {rs} {re}");
        }
        std::fs::write(p, &pseudo).ok();
    }
    let extension = extend_hint(
        read,
        &pseudo,
        extend_params,
        buf,
        DiagHint {
            offset: 0,
            num: 1,
            den: 1,
        },
    );
    let cigar = stitch_cigar(
        &extension.cigar,
        extension.ref_start,
        &exon_bounds,
        &junctions,
    );

    Some(SplicedAlignment {
        extension,
        junctions,
        cigar,
        exon_bounds,
    })
}

/// Score a tail against a reference window (direct-diagonal `extend_hint`).
pub fn score_tail_win(
    tail: &[u8],
    ref_window: &[u8],
    params: &ExtendParams,
    buf: &mut ExtendBuffer,
) -> Extension {
    extend_hint(
        tail,
        ref_window,
        params,
        buf,
        DiagHint {
            offset: 0,
            num: 1,
            den: 1,
        },
    )
}

/// One scored tail relocation candidate.
#[derive(Clone, Debug)]
pub struct SplitResolution {
    /// Library junction hosting the tail (with signal and support).
    pub junction: RefinedJunction,
    /// `score_tail_win` extension score of the tail in its window.
    pub tail_score: i32,
    /// `tail_score + signal.score + lib_bonus − SPLIT_COST − intron_length_penalty`.
    pub total: i32,
}

/// Tail-window half-width around the anchor point (frozen: ±20bp).
pub const SPLIT_LIB_WINDOW: u32 = 20;

/// Extra reference margin past the tail length in tail windows.
const TAIL_WIN_MARGIN: u32 = 30;

/// Whether an extension passes the tail acceptance gate
/// (`aligned ≥ len × 4/5` and `score > 0`).
pub fn tail_gate(ext: &Extension, tail_len: usize) -> bool {
    let aligned = ext.read_end.saturating_sub(ext.read_start) as usize;
    aligned * 5 >= tail_len * 4 && ext.score > 0
}

fn lib_bonus(support: u32) -> i32 {
    (support.min(10) * 5) as i32
}

/// Resolve the best library placement for a 3' (right) tail.
///
/// The read is in alignment orientation (for a minus-strand placement that is
/// the reverse complement of the original read), so the tail always extends
/// genomically right of the body: candidates are library junctions whose
/// start (donor) lies in `body_end ± 20`, junction strand matching the
/// placement. Each is filtered by min/max intron and S14, scored by
/// `score_tail_win` in the post-acceptor exon window under the
/// `aligned >= 4/5 && score > 0` gate. Total =
/// `tail + signal + lib_bonus - SPLIT_COST - intron_length_penalty`; the best
/// strictly-greater total wins (ties keep the first candidate in library
/// order, which is deterministic).
#[allow(clippy::too_many_arguments)]
pub fn resolve_split_right(
    reference: &Reference,
    lib: &JunctionLib,
    tail: &[u8],
    contig: u32,
    strand: Strand,
    _body_start: u32,
    body_end: u32,
    params: &SpliceParams,
    extend_params: &ExtendParams,
    buf: &mut ExtendBuffer,
) -> Option<SplitResolution> {
    let minus = strand == Strand::Minus;
    let ctg = reference.contigs.get(contig as usize)?;
    let tl = tail.len() as u32;
    let win = SPLIT_LIB_WINDOW;

    let cands = lib_window_candidates(
        lib,
        false,
        contig,
        body_end.saturating_sub(win),
        body_end + win,
    );

    let mut best: Option<SplitResolution> = None;
    for (j, support) in cands {
        if j.minus_strand != minus {
            continue;
        }
        let intron_len = j.end - j.start;
        if intron_len < params.min_intron || intron_len > params.max_intron {
            continue;
        }
        let signal = junction_signal(ctg, j);
        if !s14_pass(intron_len, signal, support) {
            continue;
        }
        let lo = j.end;
        let hi = ((lo as u64) + (tl as u64) + (TAIL_WIN_MARGIN as u64)).min(ctg.len as u64) as u32;
        if hi <= lo {
            continue;
        }
        let window = ctg.slice_ascii(lo, hi);
        let ext = score_tail_win(tail, &window, extend_params, buf);
        if !tail_gate(&ext, tail.len()) {
            continue;
        }
        let total = ext.score
            + signal.score(params)
            + lib_bonus(support)
            - SPLIT_COST
            - intron_length_penalty(intron_len);
        if best.as_ref().is_none_or(|b| total > b.total) {
            best = Some(SplitResolution {
                junction: RefinedJunction {
                    junction: *j,
                    signal,
                    known_support: support,
                },
                tail_score: ext.score,
                total,
            });
        }
    }
    best
}

/// Resolve the best library placement for a 5' (left) tail; mirror of
/// [`resolve_split_right`] (acceptor ends in `body_start ± 20`, tail window
/// upstream of the donor).
#[allow(clippy::too_many_arguments)]
pub fn resolve_split_left(
    reference: &Reference,
    lib: &JunctionLib,
    tail: &[u8],
    contig: u32,
    strand: Strand,
    body_start: u32,
    _body_end: u32,
    params: &SpliceParams,
    extend_params: &ExtendParams,
    buf: &mut ExtendBuffer,
) -> Option<SplitResolution> {
    let minus = strand == Strand::Minus;
    let ctg = reference.contigs.get(contig as usize)?;
    let tl = tail.len() as u32;
    let win = SPLIT_LIB_WINDOW;

    let cands = lib_window_candidates(
        lib,
        true,
        contig,
        body_start.saturating_sub(win),
        body_start + win,
    );

    let mut best: Option<SplitResolution> = None;
    for (j, support) in cands {
        if j.minus_strand != minus {
            continue;
        }
        let intron_len = j.end - j.start;
        if intron_len < params.min_intron || intron_len > params.max_intron {
            continue;
        }
        let signal = junction_signal(ctg, j);
        if !s14_pass(intron_len, signal, support) {
            continue;
        }
        let hi = j.start.min(ctg.len);
        let lo = hi.saturating_sub(tl + TAIL_WIN_MARGIN);
        if hi <= lo {
            continue;
        }
        let window = ctg.slice_ascii(lo, hi);
        let ext = score_tail_win(tail, &window, extend_params, buf);
        if !tail_gate(&ext, tail.len()) {
            continue;
        }
        let total = ext.score
            + signal.score(params)
            + lib_bonus(support)
            - SPLIT_COST
            - intron_length_penalty(intron_len);
        if best.as_ref().is_none_or(|b| total > b.total) {
            best = Some(SplitResolution {
                junction: RefinedJunction {
                    junction: *j,
                    signal,
                    known_support: support,
                },
                tail_score: ext.score,
                total,
            });
        }
    }
    best
}

/// Signal of a library junction read from the reference.
pub(crate) fn junction_signal(ctg: &crate::fasta::Contig, j: &Junction) -> SpliceSignal {
    let donor = dinuc(ctg, j.start, true);
    let acceptor = dinuc(ctg, j.end, false);
    splice_signal(donor, acceptor, j.minus_strand)
}

/// Junctions with `start` (or `end` when `by_end`) in `[lo, hi)`, zipped with
/// their support counts, in library order.
pub(crate) fn lib_window_candidates(
    lib: &JunctionLib,
    by_end: bool,
    contig: u32,
    lo: u32,
    hi: u32,
) -> Vec<(&Junction, u32)> {
    if by_end {
        lib.range_end(contig, lo, hi)
            .iter()
            .map(|&i| (&lib.junctions[i as usize], lib.counts[i as usize]))
            .collect()
    } else {
        let (js, cs) = lib.range_start(contig, lo, hi);
        js.iter().zip(cs.iter()).map(|(j, &c)| (j, c)).collect()
    }
}
