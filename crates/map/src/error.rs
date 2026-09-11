//! Error type for the `esperanto-map` crate.
//!
//! All user-input paths (FASTA/FASTQ/index files) report errors; nothing panics
//! on malformed input. Line numbers are 1-based where present.

use thiserror::Error;

/// Top-level error for reference parsing, index (de)serialization and read IO.
#[derive(Debug, Error)]
pub enum AlignError {
    /// Failed to open or read a FASTA file.
    #[error("fasta io error: {path}: {source}")]
    FastaIo {
        /// File path that failed.
        path: String,
        /// Underlying IO error.
        source: std::io::Error,
    },

    /// Malformed FASTA content (empty name, duplicate name, empty file).
    #[error("fasta format error at line {line}: {msg}")]
    FastaFormat {
        /// 1-based line number of the offending line.
        line: usize,
        /// Human-readable description.
        msg: String,
    },

    /// Failed to read or write a paidx index file.
    #[error("index io error")]
    IndexIo,

    /// Malformed index file (bad lengths, trailing garbage, bad alignment).
    #[error("index format error: {msg}")]
    IndexFormat {
        /// Human-readable description.
        msg: String,
    },

    /// Index file version mismatch.
    #[error("index version mismatch in {file}: supported version {supported}")]
    IndexVersion {
        /// File path that was rejected.
        file: String,
        /// Version this build supports.
        supported: u32,
    },

    /// Index was built against a different reference (sha256 mismatch).
    #[error("index reference mismatch: {msg}")]
    IndexReferenceMismatch {
        /// Human-readable description.
        msg: String,
    },

    /// Malformed FASTQ content (bad header, length mismatch, truncation).
    #[error("fastq format error at line {line}: {msg}")]
    FastqFormat {
        /// 1-based line number of the offending line.
        line: usize,
        /// Human-readable description.
        msg: String,
    },

    /// Failed to open or read a FASTQ/BFQ file.
    #[error("fastq io error: {path}: {source}")]
    FastqIo {
        /// File path that failed (or logical sink name).
        path: String,
        /// Underlying IO error.
        source: std::io::Error,
    },
}
