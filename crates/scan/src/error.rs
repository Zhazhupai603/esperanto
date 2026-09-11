//! esperanto-scan error types.

use thiserror::Error;

#[derive(Debug, Error)]
pub enum CallError {
    #[error("htslib: {0}")]
    Hts(#[from] rust_htslib::errors::Error),

    #[error("io error at {path}: {source}")]
    Io {
        path: String,
        #[source]
        source: std::io::Error,
    },

    #[error("snp vcf {path}: {msg}")]
    Snp { path: String, msg: String },

    #[error("fasta {path}: {msg}")]
    Fasta { path: String, msg: String },

    #[error("contig {0} not found in BAM header")]
    ContigNotFound(String),

    #[error("thread pool build: {0}")]
    Pool(#[from] rayon::ThreadPoolBuildError),
}
