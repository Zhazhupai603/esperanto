//! Run-through reinterpretation (library-only).
//!
//! A 1–4bp overhang can be aligned straight through a junction by the
//! extension (no soft clip ⇒ split rescue never triggers). `rescue_span`
//! re-reads such alignments: every library junction whose donor falls inside
//! `(aln.pos + 2, ref_end − 1)` on the same strand is tried in start order;
//! the read offset of the donor is found by re-walking the CIGAR, and the
//! tail from there is either exact-matched at the acceptor (<4bp) or
//! extension-validated (≥4bp, `aligned ≥ 4/5 && score > 0`). The body is
//! truncated at the donor offset (head soft clip preserved; reference-only
//! ops falling inside the kept region are kept whole) and an `N` joins it to
//! the tail. First success returns; no de-novo run-through is attempted.

use crate::extend::{push_op, CigarOp, ExtendBuffer};
use crate::gtf::RefinedJunction;
use crate::intron_chain::cigar_conserving;
use crate::splice::junction_signal;
use crate::split::{cigar_before_read_depth, extend_tail, SplitContext, SplitRescue};

/// Micro-tail threshold: shorter tails are exact-matched, longer ones go
/// through extension.
pub const MICRO_TAIL: usize = 4;

/// Extra reference margin past the tail length in tail windows.
const TAIL_MARGIN: u32 = 30;

/// Reference end of a CIGAR (M/D/N advance) starting at `pos`.
pub fn cigar_ref_end(pos: u32, cigar: &[CigarOp]) -> u32 {
    let mut r = pos as u64;
    for op in cigar {
        match op {
            CigarOp::Match(n) | CigarOp::Del(n) | CigarOp::RefSkip(n) => r += *n as u64,
            CigarOp::Ins(_) | CigarOp::SoftClip(_) => {}
        }
    }
    r as u32
}

/// Read offset whose aligned reference position is `ref_pos` (absolute;
/// `aln_pos` is the alignment's reference start): within an M run the offset
/// interpolates; a donor falling in a D/N gap maps to the current read
/// cursor. `None` when `ref_pos` is outside the aligned reference span.
pub fn read_offset_at_ref(cigar: &[CigarOp], aln_pos: u32, ref_pos: u32) -> Option<u32> {
    let ref_pos = ref_pos.checked_sub(aln_pos)?;
    let mut r = 0u64;
    let mut q = 0u64;
    let target = ref_pos as u64;
    for &op in cigar {
        match op {
            CigarOp::Match(n) => {
                let n = n as u64;
                if target >= r && target <= r + n {
                    return Some((q + (target - r)) as u32);
                }
                r += n;
                q += n;
            }
            CigarOp::Ins(n) | CigarOp::SoftClip(n) => q += n as u64,
            CigarOp::Del(n) | CigarOp::RefSkip(n) => {
                let n = n as u64;
                if target >= r && target < r + n {
                    return Some(q as u32);
                }
                r += n;
            }
        }
    }
    None
}

/// Span rescue over the junction library. Requires a non-empty library; the
/// caller gates soft-clip split rescue separately.
pub fn rescue_span(ctx: &SplitContext, buf: &mut ExtendBuffer) -> Option<SplitRescue> {
    if ctx.lib.is_empty() {
        return None;
    }
    let ctg = ctx.reference.contigs.get(ctx.contig as usize)?;
    let ref_end = cigar_ref_end(ctx.pos, ctx.cigar);
    if ref_end <= ctx.pos + 3 {
        return None;
    }

    // Donors strictly inside (pos + 2, ref_end - 1), start ascending.
    let lo = ctx.pos + 3;
    let hi = ref_end - 1;
    if hi <= lo {
        return None;
    }
    let cands: Vec<(&crate::gtf::Junction, u32)> = ctx
        .lib
        .range_start(ctx.contig, lo, hi)
        .0
        .iter()
        .zip(ctx.lib.range_start(ctx.contig, lo, hi).1.iter())
        .map(|(j, &c)| (j, c))
        .collect();

    for (j, support) in cands {
        if j.minus_strand != (ctx.strand == crate::seed::Strand::Minus) {
            continue;
        }
        let intron_len = j.end - j.start;
        if intron_len < ctx.splice_params.min_intron || intron_len > ctx.splice_params.max_intron {
            continue;
        }
        let Some(dro) = read_offset_at_ref(ctx.cigar, ctx.pos, j.start) else {
            continue;
        };
        let tail = &ctx.read[dro as usize..];
        if tail.is_empty() {
            continue;
        }

        // Tail placement at the acceptor.
        let tail_ops: Vec<CigarOp> = if tail.len() < MICRO_TAIL {
            if j.end as u64 + tail.len() as u64 > ctg.len as u64 {
                continue;
            }
            let acceptor = ctg.slice_ascii(j.end, j.end + tail.len() as u32);
            let exact = tail
                .iter()
                .zip(acceptor.iter())
                .all(|(a, b)| a.eq_ignore_ascii_case(b));
            if !exact {
                continue;
            }
            vec![CigarOp::Match(tail.len() as u32)]
        } else {
            let hi = ((j.end as u64) + (tail.len() as u64) + (TAIL_MARGIN as u64))
                .min(ctg.len as u64) as u32;
            let Some((_, ext)) = extend_tail(ctg, tail, j.end, hi, &ctx.extend_params, buf) else {
                continue;
            };
            ext.cigar.clone()
        };

        // Body truncated at dro (head soft clip preserved), N, tail.
        let mut cigar = cigar_before_read_depth(ctx.cigar, dro);
        push_op(&mut cigar, CigarOp::RefSkip(intron_len));
        for op in &tail_ops {
            push_op(&mut cigar, *op);
        }
        if !cigar_conserving(&cigar, ctx.read.len()) {
            continue;
        }

        let signal = junction_signal(ctg, j);
        return Some(SplitRescue {
            pos: ctx.pos,
            cigar,
            junction: RefinedJunction {
                junction: *j,
                signal,
                known_support: support,
            },
            score: 0,
        });
    }
    None
}
