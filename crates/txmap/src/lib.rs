//! # esperanto-txmap
//!
//! Transcriptome→genome projection and deterministic multi-isoform
//! attribution (L1 coordinate-mapping crate).
//!
//! ## Coordinate conventions (strict)
//!
//! * GTF input: 1-based inclusive, converted on ingest.
//! * All genomic output: 0-based half-open `[start, end)` on the **forward**
//!   reference strand (SAM convention).
//! * Transcript coordinates (`tx_start`): 0-based along the transcript 5'→3';
//!   `tx = 0` is the transcription start on either strand.
//!
//! ## Public API contract
//!
//! * [`TranscriptRecord`] holds exons in transcription order (plus strand:
//!   ascending genomic; minus strand: descending), so `exons[0]` is always
//!   the transcription-start end. `project` maps a transcript interval to
//!   `(genomic_start, CIGAR)` with minus-strand reflection, overhang
//!   clamping at the transcript end, and `Match`/`RefSkip` CIGAR over
//!   genomically sorted pieces.
//! * [`TxMap`] is the transcript index: built from records or a GTF file
//!   (`from_gtf`, source hash = BLAKE3 of raw GTF bytes), persisted via the
//!   deterministic little-endian `.txmap` format (`save`/`open`), and
//!   queried through `project` / accessors.
//! * [`attribute`] (exposed as `TxMap::attribute`) collapses a read's
//!   candidate projections: merge identical `(contig, pos, cigar)` groups →
//!   junction-support score decides → ties resolved deterministically by
//!   lexicographic `(contig name, pos, cigar string)` with MAPQ 0.
//!   MAPQ otherwise = `min(60, floor(-10·log10(1 − w)))` with
//!   `w = |winner members| / Σ|all members|`.
//! * [`JunctionSet`] is the sorted, deduplicated intron set with
//!   `contains(contig_id, start, end)` lookups.
//!
//! ## Invariants
//!
//! * Determinism: identical inputs produce byte-identical `.txmap` files —
//!   transcripts sorted by name, contigs sorted lexicographically, junctions
//!   sorted; no timestamps, no padding, no unsorted map iteration in output.
//! * All stored exons are legal half-open intervals, non-overlapping and in
//!   transcription order (enforced by `TranscriptRecord::validate`).
//! * No `unsafe` code.

#![deny(unsafe_code)]

mod attribution;
mod cigar;
mod gtf;
mod junction;
mod serialize;
mod transcript;

pub use crate::attribution::{Attribution, AttributionCandidate, Placement};
pub use crate::cigar::{cigar_string, CigarOp};
pub use crate::junction::{Junction, JunctionSet};
pub use crate::serialize::TxMap;
pub use crate::transcript::{Exon, Strand, TranscriptRecord};

use thiserror::Error;

/// Errors raised by this crate.
#[derive(Debug, Error)]
pub enum Error {
    /// Underlying I/O failure.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    /// File does not start with the `.txmap` magic bytes.
    #[error(
        "bad magic: expected {}, found {}",
        String::from_utf8_lossy(expected),
        String::from_utf8_lossy(found)
    )]
    Magic {
        /// Expected magic constant (`"TXMAP001"`).
        expected: [u8; 8],
        /// Magic bytes found in the file.
        found: [u8; 8],
    },
    /// Unsupported `.txmap` format version.
    #[error("unsupported version: file has {file}, code has {code}")]
    Version {
        /// Version field read from the file.
        file: u32,
        /// Version this code reads/writes.
        code: u32,
    },
    /// Malformed input: bad GTF line, invalid transcript record, or corrupt
    /// binary section.
    #[error("format error: {0}")]
    Format(String),
}
