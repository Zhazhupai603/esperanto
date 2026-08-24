//! Genomic projection and deterministic attribution.

use crate::{CigarOp, ExtEntry, Strand, TxMap};

/// Genomic strand of a read mapping with `read_strand` on a transcript of
/// strand `tx_strand`: the read orientation flips on minus transcripts.
pub fn to_genomic_strand(read_strand: Strand, tx_strand: Strand) -> Strand {
    match tx_strand {
        Strand::Plus => read_strand,
        Strand::Minus => read_strand.flip(),
    }
}

/// A projected, scored candidate ready for [`finalize_candidates`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Placed {
    /// Contig id.
    pub contig: u32,
    /// Leftmost 0-based genomic position.
    pub pos: u32,
    /// Genomic strand of the read.
    pub strand: Strand,
    /// CIGAR in reference left-to-right orientation.
    pub cigar: Vec<CigarOp>,
    /// Alignment score (0 on the full branch, EA distance on partials).
    pub score: i32,
    /// Source transcript id (attribution tie-break key).
    pub tx_id: u32,
    /// Extension diagonal (attribution tie-break key).
    pub diagonal: i64,
    /// Read bases covered by the winning extension.
    pub read_cov: u32,
}

/// Annotated strand of `tx_id`; the trait's `None` default reads as `Plus`.
fn tx_strand_of(txmap: &impl TxMap, tx_id: u32) -> Strand {
    txmap.strand(tx_id).unwrap_or(Strand::Plus)
}

/// Project a full-length extension: `tx_start = ext.tx_lo` (full coverage
/// implies `read_lo == 0`), score 0.
pub(crate) fn project_full(txmap: &impl TxMap, e: &ExtEntry, read_len: usize) -> Option<Placed> {
    let (contig, pos, cigar) = txmap.project(e.tx_id, e.ext.tx_lo as u32, read_len as u32)?;
    Some(Placed {
        contig,
        pos,
        strand: to_genomic_strand(e.strand, tx_strand_of(txmap, e.tx_id)),
        cigar,
        score: 0,
        tx_id: e.tx_id,
        diagonal: e.diagonal,
        read_cov: read_len as u32,
    })
}

/// Partial-branch anchor start: `min(max(anchor, 0), tx_len - read_len)`
/// with `anchor = tx_lo - read_lo` (never negative).
pub(crate) fn partial_tx_start(txmap: &impl TxMap, e: &ExtEntry, read_len: usize) -> Option<u32> {
    let tx_len = txmap.tx_len(e.tx_id)? as i64;
    let hi = (tx_len - read_len as i64).max(0);
    Some(e.diagonal.max(0).min(hi) as u32)
}

/// Projected genomic locus `(contig, pos)` of a partial candidate.
pub(crate) fn partial_locus(
    txmap: &impl TxMap,
    e: &ExtEntry,
    read_len: usize,
) -> Option<(u32, u32)> {
    let ts = partial_tx_start(txmap, e, read_len)?;
    let (contig, pos, _) = txmap.project(e.tx_id, ts, read_len as u32)?;
    Some((contig, pos))
}

/// Project a verified partial candidate with the anchor-clamped start.
pub(crate) fn project_partial(
    txmap: &impl TxMap,
    e: &ExtEntry,
    read_len: usize,
    score: i32,
) -> Option<Placed> {
    let ts = partial_tx_start(txmap, e, read_len)?;
    let (contig, pos, cigar) = txmap.project(e.tx_id, ts, read_len as u32)?;
    Some(Placed {
        contig,
        pos,
        strand: to_genomic_strand(e.strand, tx_strand_of(txmap, e.tx_id)),
        cigar,
        score,
        tx_id: e.tx_id,
        diagonal: e.diagonal,
        read_cov: e.ext.read_cov() as u32,
    })
}

/// Deterministic attribution over candidates already sorted by
/// `(score asc, tx_id asc, diagonal asc)`.
///
/// MAPQ rules: a single candidate scores 60; a strictly better first
/// candidate scores 60; all candidates tied at the best score sharing the
/// winner's `(contig, pos, cigar)` score 60; anything else scores 0.
/// `force_mapq0` (cluster competition / repeat regions) pins the winner's
/// MAPQ to 0 while still placing it.
pub fn finalize_candidates(sorted: &[Placed], force_mapq0: bool) -> Option<(&Placed, u8)> {
    let best = sorted.first()?;
    let mapq = if sorted.len() == 1 || sorted[1].score != best.score {
        60
    } else {
        let tied_same = sorted
            .iter()
            .filter(|c| c.score == best.score)
            .all(|c| c.contig == best.contig && c.pos == best.pos && c.cigar == best.cigar);
        if tied_same {
            60
        } else {
            0
        }
    };
    let mapq = if force_mapq0 { 0 } else { mapq };
    Some((best, mapq))
}
