//! Core transcript model: strand, exons, transcript length, and the
//! transcript→genome projection that turns a transcript interval into a
//! forward-strand genomic placement with a Match/RefSkip CIGAR.

use crate::cigar::CigarOp;
use crate::Error;

/// Gene strand as annotated in the source GTF.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strand {
    /// Forward strand (`+` in GTF).
    Plus,
    /// Reverse strand (`-` in GTF).
    Minus,
}

/// One exon in forward-strand genomic coordinates,
/// 0-based half-open `[g_start, g_end)` with `g_start < g_end`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Exon {
    pub g_start: u32,
    pub g_end: u32,
}

/// A transcript with its exons stored in transcription order (5'→3'):
/// ascending genomic order on the plus strand, descending on the minus
/// strand, so `exons[0]` is always the transcription-start end.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TranscriptRecord {
    /// Transcript id (e.g. GENCODE `transcript_id`). Unique within a [`crate::TxMap`].
    pub name: String,
    /// Contig / chromosome name, forward reference strand.
    pub contig: String,
    pub strand: Strand,
    /// Exons in transcription order (see struct docs for the ordering invariant).
    pub exons: Vec<Exon>,
}

impl TranscriptRecord {
    /// Total transcript length in bases (sum of exon lengths).
    pub fn tx_len(&self) -> u64 {
        self.exons
            .iter()
            .map(|e| (e.g_end - e.g_start) as u64)
            .sum()
    }

    /// Validate the record invariants:
    /// non-empty `name`/`contig`, at least one exon, every exon a legal
    /// half-open interval, and adjacent exons (in transcription order)
    /// non-overlapping and correctly ordered. Returns `Error::Format` on
    /// the first violation.
    pub fn validate(&self) -> Result<(), Error> {
        if self.name.is_empty() {
            return Err(Error::Format("transcript with empty name".to_string()));
        }
        if self.contig.is_empty() {
            return Err(Error::Format(format!(
                "transcript '{}' with empty contig",
                self.name
            )));
        }
        if self.exons.is_empty() {
            return Err(Error::Format(format!(
                "transcript '{}' has no exons",
                self.name
            )));
        }
        for (i, exon) in self.exons.iter().enumerate() {
            if exon.g_start >= exon.g_end {
                return Err(Error::Format(format!(
                    "transcript '{}': exon {} is not a legal half-open interval [{}, {})",
                    self.name, i, exon.g_start, exon.g_end
                )));
            }
        }
        for i in 1..self.exons.len() {
            let (prev, cur) = (self.exons[i - 1], self.exons[i]);
            let ok = match self.strand {
                // Plus: exons ascend the genome without overlapping.
                Strand::Plus => cur.g_start >= prev.g_end,
                // Minus: exons descend the genome without overlapping.
                Strand::Minus => cur.g_end <= prev.g_start,
            };
            if !ok {
                return Err(Error::Format(format!(
                    "transcript '{}': exons {} and {} overlap or are out of transcription order",
                    self.name,
                    i - 1,
                    i
                )));
            }
        }
        Ok(())
    }

    /// Intron pairs between transcription-adjacent exons, stored RAW as
    /// `(exon[i].g_end, exon[i+1].g_start)` — the reference tool's exact
    /// representation. On the plus strand this is the legal forward intron
    /// `(start, end)`; on the minus strand (exons genomically descending) it
    /// is the inverted pair `(end, start)`. The inverted entries are kept
    /// bug-compatibly: they are part of the serialized junction table.
    /// Intron pairs between transcription-adjacent exons as FORWARD 0-based
    /// half-open intervals `(start, end)` with `start <= end`, independent of
    /// strand. The old tool stored raw transcription-order pairs, which for
    /// minus-strand transcripts (exons genomically descending) produced
    /// inverted `(end, start)` entries that could never match the forward
    /// intron queries used by junction-support scoring — minus-strand
    /// annotated junctions silently contributed zero support. ESPERANTO
    /// 1.0.0 normalizes (scientific-priority deviation, see
    /// docs/SCIENCE-DEVIATIONS.md); the reader still accepts legacy inverted
    /// entries and normalizes on load.
    pub fn introns(&self) -> Vec<(u32, u32)> {
        self.exons
            .windows(2)
            .map(|w| {
                let (a, b) = (w[0].g_end, w[1].g_start);
                if a <= b { (a, b) } else { (b, a) }
            })
            .collect()
    }

    /// Project the transcript interval `[tx_start, tx_start + len)` onto the
    /// forward genome strand.
    ///
    /// Returns `(genomic_start, cigar)` where `genomic_start` is the leftmost
    /// projected base and the CIGAR is one `Match` per exon piece separated by
    /// `RefSkip` intron gaps, or `None` when `len == 0` or
    /// `tx_start >= tx_len()`. Overhang past the transcript end is clamped to
    /// the transcript end (not an error): poly(A) tails and slight overhangs
    /// are real biology.
    pub fn project(&self, tx_start: u32, len: u32) -> Option<(u32, Vec<CigarOp>)> {
        if len == 0 {
            return None;
        }
        let tx_len = self.tx_len();
        let tx_start = tx_start as u64;
        if tx_start >= tx_len {
            return None;
        }
        // Overhang clamp: never read past the transcript end.
        let end = (tx_start + len as u64).min(tx_len);

        // Walk exons in transcription order, intersecting each exon's
        // transcript interval with [tx_start, end).
        let mut pieces: Vec<(u64, u64)> = Vec::new();
        let mut offset: u64 = 0; // transcript coordinate of the current exon start
        for exon in &self.exons {
            let exon_len = (exon.g_end - exon.g_start) as u64;
            let lo = offset.max(tx_start);
            let hi = (offset + exon_len).min(end);
            if lo < hi {
                let r_lo = lo - offset;
                let r_hi = hi - offset;
                // Reflect offsets inside the exon on the minus strand.
                let piece = match self.strand {
                    Strand::Plus => (exon.g_start as u64 + r_lo, exon.g_start as u64 + r_hi),
                    Strand::Minus => (exon.g_end as u64 - r_hi, exon.g_end as u64 - r_lo),
                };
                pieces.push(piece);
            }
            offset += exon_len;
            if offset >= end {
                break;
            }
        }

        // Minus-strand pieces come out in descending genomic order; restore
        // ascending forward-strand order (validated records never overlap).
        pieces.sort_unstable_by_key(|p| p.0);
        let genomic_start = pieces.first()?.0 as u32;

        let mut cigar: Vec<CigarOp> = Vec::with_capacity(pieces.len() * 2);
        let mut prev_end: Option<u64> = None;
        for &(start, stop) in &pieces {
            // Only emit RefSkip for a real intron gap (start > prev_end).
            // Contiguous exons (gap 0) yield adjacent Match ops; overlapping
            // pieces (invalid records) would underflow — skip the op instead.
            if let Some(pe) = prev_end {
                if start > pe {
                    cigar.push(CigarOp::RefSkip((start - pe) as u32));
                }
            }
            cigar.push(CigarOp::Match((stop - start) as u32));
            prev_end = Some(stop);
        }
        Some((genomic_start, cigar))
    }
}
