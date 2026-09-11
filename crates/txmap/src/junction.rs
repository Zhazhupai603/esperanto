//! Known-splice-junction set used for attribution tie-breaking.

use crate::cigar::CigarOp;

/// A known intron on one contig, forward 0-based half-open `[start, end)`.
/// `contig_id` indexes into [`crate::TxMap::contigs`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Junction {
    pub contig_id: u32,
    pub start: u32,
    pub end: u32,
}

/// Set of known junctions over contig ids, sorted and deduplicated at build
/// time so iteration order is deterministic.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct JunctionSet {
    inner: Vec<Junction>,
}

impl JunctionSet {
    /// Build from an arbitrary iterator of junctions; the result is sorted
    /// by `(contig_id, start, end)` and deduplicated.
    pub fn from_junctions<I: IntoIterator<Item = Junction>>(junctions: I) -> Self {
        let mut inner: Vec<Junction> = junctions.into_iter().collect();
        inner.sort_unstable();
        inner.dedup();
        Self { inner }
    }

    /// True iff `(contig_id, start, end)` is a known junction.
    pub fn contains(&self, contig_id: u32, start: u32, end: u32) -> bool {
        self.inner
            .binary_search(&Junction {
                contig_id,
                start,
                end,
            })
            .is_ok()
    }

    /// All junctions in sorted order.
    pub fn junctions(&self) -> &[Junction] {
        &self.inner
    }

    /// Number of junctions.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// True iff the set is empty.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Junction-support score of a candidate placement: the number of
    /// `RefSkip` intervals of `cigar` placed at `pos` on `contig_id` that hit
    /// a known junction. Internal helper for attribution.
    pub(crate) fn support_score(&self, contig_id: u32, pos: u32, cigar: &[CigarOp]) -> u32 {
        let mut cursor: u64 = pos as u64;
        let mut score: u32 = 0;
        for op in cigar {
            match op {
                CigarOp::Match(len) => cursor += *len as u64,
                CigarOp::RefSkip(len) => {
                    let intron_start = cursor;
                    let intron_end = cursor + *len as u64;
                    if intron_start <= u32::MAX as u64
                        && intron_end <= u32::MAX as u64
                        && self.contains(contig_id, intron_start as u32, intron_end as u32)
                    {
                        score += 1;
                    }
                    cursor = intron_end;
                }
            }
        }
        score
    }
}
