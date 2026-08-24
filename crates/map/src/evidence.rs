//! Fixed-format evidence tag.
//!
//! Field order is frozen (consumers parse this string):
//! `seeds=..;chain=..;sub=..;dq=..;splice=..;mm=..;ea=..;mapq_src=..`.

use crate::mapq::ReadAlignment;

/// Render the EK evidence tag for an alignment.
pub fn evidence_tag(aln: &ReadAlignment) -> String {
    let splice = aln
        .junctions
        .first()
        .map(|j| j.signal.label())
        .unwrap_or("-");
    let mapq_src = if aln.second_chain_score > 0 {
        "chain_margin"
    } else {
        "unique"
    };
    let dq = (aln.chain_score - aln.second_chain_score).max(0);
    format!(
        "seeds={};chain={};sub={};dq={};splice={};mm={};ea={};mapq_src={}",
        aln.n_seeds,
        aln.chain_score,
        aln.second_chain_score,
        dq,
        splice,
        aln.mm_count,
        aln.ea_count,
        mapq_src
    )
}
