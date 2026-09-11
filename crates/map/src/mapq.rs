//! Mapping quality.
//!
//! `mapq = round(60 × (1 − second/best) × min(1, n_anchors/10))` with
//! `best <= 0 ⇒ 0` and `second = 0 ⇒ factor 1`. Track-2 synthetic records use
//! the same formula (see the pipeline quirk note in the spec).
//!
//! This module also hosts the crate-wide [`ReadAlignment`] record (spec:
//! shared core type), since MAPQ is its first consumer.

use crate::extend::CigarOp;
use crate::gtf::RefinedJunction;
use crate::seed::Strand;

/// A fully placed read alignment (all field semantics frozen by the spec).
#[derive(Clone, Debug, Default)]
pub struct ReadAlignment {
    /// Contig index.
    pub contig: u32,
    /// 0-based inclusive leftmost reference position.
    pub pos: u32,
    /// Alignment strand.
    pub strand: Strand,
    /// Extension score.
    pub score: i32,
    /// Best chain score (MAPQ input).
    pub chain_score: i32,
    /// Best non-overlapping second chain score (MAPQ input).
    pub second_chain_score: i32,
    /// CIGAR over the whole read.
    pub cigar: Vec<CigarOp>,
    /// Number of anchors in the best chain.
    pub n_anchors: usize,
    /// Splice junctions crossed (empty = unspliced).
    pub junctions: Vec<RefinedJunction>,
    /// Editing-aware tolerated A>G/T>C count (EA tag).
    pub ea_count: u32,
    /// Ordinary mismatch count (EK mm).
    pub mm_count: u32,
    /// Seed count (EK seeds).
    pub n_seeds: usize,
    /// Rescued / re-seeded / Track-2 record (RE tag).
    pub rescued: bool,
}

/// Compute MAPQ from chain scores and anchor count.
pub fn mapq(chain_score: i32, second_chain_score: i32, n_anchors: usize) -> u8 {
    if chain_score <= 0 {
        return 0;
    }
    let ratio = if second_chain_score <= 0 {
        1.0
    } else {
        1.0 - (second_chain_score as f64 / chain_score as f64)
    };
    let anchor_factor = (n_anchors as f64 / 10.0).min(1.0);
    let v = 60.0 * ratio * anchor_factor;
    let q = v.round();
    if q < 0.0 {
        0
    } else if q > 60.0 {
        60
    } else {
        q as u8
    }
}

/// MAPQ for a placed alignment.
pub fn mapq_of(aln: &ReadAlignment) -> u8 {
    mapq(aln.chain_score, aln.second_chain_score, aln.n_anchors)
}
