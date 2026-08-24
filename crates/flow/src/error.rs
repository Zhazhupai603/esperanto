//! Flow error type (spec §errors): downstream semantics are wrapped, never
//! swallowed.

use std::io;

/// Pipeline-level error.
#[derive(Debug, thiserror::Error)]
pub enum FlowError {
    /// Invalid input field combination or missing entry-required parameter.
    #[error("invalid entry: {0}")]
    Entry(String),
    /// Species/version guardrail rejection (chr1 length mismatch).
    #[error(
        "species/version mismatch: reference chr1 length = {len}, \
         expected 248956422 (hg38) or < 10000000 (synthetic/test); \
         provide an hg38 reference or a bundle matching your reference"
    )]
    SpeciesMismatch {
        /// Observed chr1 length from `<fasta>.fai`.
        len: u64,
    },
    /// Bam/BamSites entry BAM lacks a `.bai`/`.csi` index.
    #[error("BAM index missing for {path}: coordinate-sorted BAM with .bai/.csi is required")]
    MissingBamIndex {
        /// Input BAM path.
        path: String,
    },
    /// candidates.bed / sites parse failure.
    #[error("site parse error at line {line}: {msg}")]
    BedParse {
        /// 1-based line number.
        line: usize,
        /// What went wrong.
        msg: String,
    },
    /// Plain I/O failure.
    #[error("io: {0}")]
    Io(#[from] io::Error),
    /// A stage crate failed; the original error is preserved as source.
    #[error("stage {stage} failed: {source}")]
    Stage {
        /// Stage name (qc / map / sort / scan / score / vcf).
        stage: &'static str,
        /// Downstream error, unmodified.
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },
}
