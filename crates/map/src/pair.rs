//! Paired-end relationship and insert-size statistics.

use crate::extend::CigarOp;
use crate::mapq::ReadAlignment;

/// Maximum reference span for a proper pair (frozen: 1000).
pub const PROPER_SPAN_MAX: u32 = 1000;

/// Relationship of two mates.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PairRelation {
    /// Mates on opposite strands (FR or RF) on the same contig within the
    /// transcript-aware span cap.
    Proper,
    /// Anything else.
    Discordant,
}

/// Classify a mate pair: proper iff same contig, mates on opposite strands
/// (FR or RF — library strand specificity flips which orientation dominates,
/// but either is a valid inward-facing pair), and the transcript-aware span
/// (genomic span minus intronic N-skips) is within the span cap.
pub fn relate(r1: &ReadAlignment, r2: &ReadAlignment) -> PairRelation {
    if r1.contig != r2.contig {
        return PairRelation::Discordant;
    }
    // Same-strand (FF/RR) is discordant; opposite-strand (FR/RF) is a
    // valid inward-facing pair.
    if r1.strand == r2.strand {
        return PairRelation::Discordant;
    }
    if transcript_span(r1, r2) > PROPER_SPAN_MAX {
        return PairRelation::Discordant;
    }
    PairRelation::Proper
}

/// Span of the pair on the transcript: reference span from R1 start to R2 end
/// with intronic RefSkip (N) intervals removed. RNA mates routinely bridge
/// kilobase introns; counting N toward the span would reject spliced proper
/// pairs that are adjacent on the transcript.
pub fn transcript_span(r1: &ReadAlignment, r2: &ReadAlignment) -> u32 {
    let left = r1.pos.min(r2.pos);
    let right = ref_end(r1).max(ref_end(r2));
    if right <= left {
        return 0;
    }
    let mut n_total = 0u64;
    for aln in [r1, r2] {
        let mut p = aln.pos;
        for op in &aln.cigar {
            match op {
                CigarOp::RefSkip(n) => {
                    let s = p.max(left);
                    let e = (p + *n).min(right);
                    if e > s {
                        n_total += u64::from(e - s);
                    }
                    p += *n;
                }
                CigarOp::Match(n) | CigarOp::Del(n) => p += *n,
                CigarOp::Ins(_) | CigarOp::SoftClip(_) => {}
            }
        }
    }
    ((right - left) as u64).saturating_sub(n_total) as u32
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
/// Magnitudes are transcript-aware (intronic N-skips excluded), consistent
/// with [`transcript_span`].
pub fn template_length(r1: &ReadAlignment, r2: &ReadAlignment) -> i32 {
    if r1.contig != r2.contig {
        return 0;
    }
    if r1.pos <= r2.pos {
        transcript_span(r1, r2) as i32
    } else {
        -(transcript_span(r1, r2) as i32)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seed::Strand;

    fn aln(pos: u32, strand: Strand, cigar: Vec<CigarOp>) -> ReadAlignment {
        ReadAlignment {
            contig: 0,
            pos,
            strand,
            cigar,
            ..Default::default()
        }
    }

    #[test]
    fn opposite_strand_pairs_are_proper_in_either_order() {
        let a = aln(1000, Strand::Plus, vec![CigarOp::Match(150)]);
        let b = aln(1200, Strand::Minus, vec![CigarOp::Match(150)]);
        assert_eq!(relate(&a, &b), PairRelation::Proper);
        assert_eq!(relate(&b, &a), PairRelation::Proper);
        assert_eq!(template_length(&a, &b), 350);
        assert_eq!(template_length(&b, &a), -350);
    }

    #[test]
    fn same_strand_pairs_are_discordant() {
        let a = aln(1000, Strand::Plus, vec![CigarOp::Match(150)]);
        let b = aln(1200, Strand::Plus, vec![CigarOp::Match(150)]);
        assert_eq!(relate(&a, &b), PairRelation::Discordant);
    }

    #[test]
    fn intron_skips_do_not_count_toward_the_span() {
        // R1 bridges a 5 kb intron; the transcript-level span is 350 bp.
        let a = aln(
            1000,
            Strand::Plus,
            vec![CigarOp::Match(100), CigarOp::RefSkip(5000), CigarOp::Match(50)],
        );
        let b = aln(6200, Strand::Minus, vec![CigarOp::Match(150)]);
        assert_eq!(transcript_span(&a, &b), 350);
        assert_eq!(relate(&a, &b), PairRelation::Proper);
    }

    #[test]
    fn distant_mates_are_discordant() {
        let a = aln(1000, Strand::Plus, vec![CigarOp::Match(150)]);
        let b = aln(5000, Strand::Minus, vec![CigarOp::Match(150)]);
        assert_eq!(relate(&a, &b), PairRelation::Discordant);
    }

    #[test]
    fn different_contigs_are_discordant() {
        let a = aln(1000, Strand::Plus, vec![CigarOp::Match(150)]);
        let mut b = aln(1200, Strand::Minus, vec![CigarOp::Match(150)]);
        b.contig = 1;
        assert_eq!(relate(&a, &b), PairRelation::Discordant);
    }
}
