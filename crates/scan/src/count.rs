//! Shared helper subset (the common pieces of the original count.rs; the legacy main scan engine is not ported).

use std::collections::BTreeSet;
use std::path::Path;

/// Whether the fasta contains a contig: read the .fai first column (htslib's behavior on a missing contig is panic-level,
/// so a name lookup must come first; the .fai is tiny and page-cached, so the cost is negligible).
pub(crate) fn fasta_has_contig(fa: &Path, chrom: &str) -> bool {
    let fai = format!("{}.fai", fa.display());
    let Ok(text) = std::fs::read_to_string(&fai) else {
        return false;
    };
    text.lines().any(|l| l.split('\t').next() == Some(chrom))
}

/// Shared contig reference-sequence cache: deduplicates whole-contig fasta loads across block tasks.
pub type SeqCache = std::sync::RwLock<std::collections::HashMap<String, std::sync::Arc<Vec<u8>>>>;

/// Maximum distance (bp) for junction evidence.
pub(crate) const JUNCTION_EVID_DIST: u32 = 8;

/// Minimum run length for the homopolymer evidence code.
pub(crate) const HP_MIN: u32 = 4;

fn base_idx(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Length of the same-base run containing pos (ACGT only, capped at 10).
pub(crate) fn hp_len(seq: &[u8], pos: usize) -> u32 {
    let b = seq[pos];
    if base_idx(b).is_none() {
        return 0;
    }
    let mut l = pos;
    while l > 0 && seq[l - 1] == b {
        l -= 1;
    }
    let mut r = pos;
    while r + 1 < seq.len() && seq[r + 1] == b {
        r += 1;
    }
    ((r - l + 1) as u32).min(10)
}

/// Distance from pos to the nearest boundary in the jbounds set.
pub(crate) fn nearest_dist(bounds: &BTreeSet<i64>, pos: i64) -> Option<u32> {
    let below = bounds.range(..=pos).next_back().map(|b| pos - b);
    let above = bounds.range(pos..).next().map(|b| b - pos);
    match (below, above) {
        (Some(a), Some(b)) => Some((a.min(b)) as u32),
        (Some(a), None) => Some(a as u32),
        (None, Some(b)) => Some(b as u32),
        (None, None) => None,
    }
}
