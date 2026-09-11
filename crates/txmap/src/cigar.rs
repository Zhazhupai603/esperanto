//! CIGAR representation for projected placements.

use std::fmt;
use std::fmt::Write as _;

/// A single CIGAR operation. Only the two ops produced by projection exist
/// today: `Match` (aligned block, SAM `M`) and `RefSkip` (intron gap, SAM `N`).
/// The enum is extensible if later layers need more ops.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CigarOp {
    /// Aligned block of the given length.
    Match(u32),
    /// Reference skip (intron) of the given length.
    RefSkip(u32),
}

impl fmt::Display for CigarOp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CigarOp::Match(len) => write!(f, "{len}M"),
            CigarOp::RefSkip(len) => write!(f, "{len}N"),
        }
    }
}

/// Render a CIGAR as its SAM-style string, e.g. `25M100N25M`.
///
/// The rendering is injective: equal strings imply equal op sequences, so it
/// is safe to use as a comparison key (attribution tie-breaking).
pub fn cigar_string(ops: &[CigarOp]) -> String {
    let mut out = String::with_capacity(ops.len() * 4);
    for op in ops {
        let _ = write!(out, "{op}");
    }
    out
}
