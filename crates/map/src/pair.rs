//! Paired-end relationship and insert-size statistics.

use crate::extend::CigarOp;
use crate::mapq::ReadAlignment;
use crate::seed::Strand;

/// Maximum reference span for a proper pair (frozen: 1000).
pub const PROPER_SPAN_MAX: u32 = 1000;

/// Relationship of two mates.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PairRelation {
    /// R1 left/Plus, R2 right/Minus on the same contig within the span cap.
    Proper,
    /// Anything else.
    Discordant,
}

/// Classify a mate pair: proper iff same contig, R1 on the left on the Plus
/// strand, R2 on the right on the Minus strand, and the span ≤ 1000.
pub fn relate(r1: &ReadAlignment, r2: &ReadAlignment) -> PairRelation {
    if r1.contig != r2.contig {
        return PairRelation::Discordant;
    }
    if r1.strand != Strand::Plus || r2.strand != Strand::Minus {
        return PairRelation::Discordant;
    }
    if r1.pos > r2.pos {
        return PairRelation::Discordant;
    }
    if ref_end(r2).saturating_sub(r1.pos) > PROPER_SPAN_MAX {
        return PairRelation::Discordant;
    }
    PairRelation::Proper
}

/// One-past-the-end reference position: `pos` advanced by M/D/N ops.
pub fn ref_end(aln: &ReadAlignment) -> u32 {
    let mut end = aln.pos;
    for op in &aln.cigar {
        match op {
            CigarOp::Match(n) | CigarOp::Del(n) | CigarOp::RefSkip(n) => end += n,
            CigarOp::Ins(_) | CigarOp::SoftClip(_) => {}
        }
    }
    end
}

/// Signed template length from R1's viewpoint: positive when R1 is the left
/// mate, negative (magnitude = right-pair view) when R1 is on the right.
pub fn template_length(r1: &ReadAlignment, r2: &ReadAlignment) -> i32 {
    if r1.contig != r2.contig {
        return 0;
    }
    if r1.pos <= r2.pos {
        (ref_end(r2).saturating_sub(r1.pos)) as i32
    } else {
        -((ref_end(r1).saturating_sub(r2.pos)) as i32)
    }
}

/// Integer-friendly insert statistics over |tlen| values (population
/// moments: `stdev = sqrt(max(0, E[x²] − m²))`).
#[derive(Clone, Copy, Debug, Default)]
pub struct InsertStats {
    /// Number of observations.
    pub count: u64,
    /// Σx.
    pub sum: f64,
    /// Σx².
    pub sum_sq: f64,
}

impl InsertStats {
    /// Record one |tlen| observation.
    pub fn push(&mut self, tlen: i32) {
        let x = tlen.unsigned_abs() as f64;
        self.count += 1;
        self.sum += x;
        self.sum_sq += x * x;
    }

    /// Mean of the observations (0 when empty).
    pub fn mean(&self) -> f64 {
        if self.count == 0 {
            0.0
        } else {
            self.sum / self.count as f64
        }
    }

    /// Population standard deviation (`sqrt(max(0, E[x²] − m²))`).
    pub fn stdev(&self) -> f64 {
        if self.count == 0 {
            return 0.0;
        }
        let n = self.count as f64;
        let ex2 = self.sum_sq / n;
        let m = self.sum / n;
        (ex2 - m * m).max(0.0).sqrt()
    }
}
