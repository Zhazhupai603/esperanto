//! Fixed-format evidence tag.
//!
//! Field order is frozen (consumers parse this string):
//! `seeds=..;chain=..;sub=..;dq=..;splice=..;mm=..;ea=..;mapq_src=..;origin=..;cand=..`.
//!
//! `chain`/`sub`/`dq` report the best/second pair that produced the final
//! MAPQ. `mapq_src` names it: `ext_margin` for the extension pair,
//! `chain_margin` for the chain pair, `unique` for no runner-up in either,
//! `policy` for a caller-fixed MAPQ ([`crate::mapq::ReadAlignment::mapq_override`]).
//! The pair comes from [`crate::mapq::margin`], so the tag cannot disagree
//! with the MAPQ it explains. `origin` names the layer whose candidate won
//! the unified strict arbitration and `cand` counts the distinct-locus
//! candidates that arbitration saw.

use crate::mapq::{margin, ReadAlignment};

/// Render the EK evidence tag for an alignment.
pub fn evidence_tag(aln: &ReadAlignment) -> String {
    let splice = aln
        .junctions
        .first()
        .map(|j| j.signal.label())
        .unwrap_or("-");
    let (best, second, mapq_src) = margin(aln);
    let dq = (best - second).max(0);
    format!(
        "seeds={};chain={};sub={};dq={};splice={};mm={};ea={};mapq_src={};origin={};cand={}",
        aln.n_seeds,
        best,
        second,
        dq,
        splice,
        aln.mm_count,
        aln.ea_count(),
        mapq_src,
        aln.origin.label(),
        aln.n_cands
    )
}
