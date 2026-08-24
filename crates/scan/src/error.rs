//! esperanto-scan error types (thiserror enum; user-supplied paths never panic; style mirrors esperanto-map).

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

    #[error("call spec error at {path}: {msg}")]
    Spec { path: String, msg: String },

    #[error("annotation error at {path}: {msg}")]
    Annot { path: String, msg: String },

    #[error("contig {0} not found in BAM header")]
    ContigNotFound(String),

    #[error("thread pool build: {0}")]
    Pool(#[from] rayon::ThreadPoolBuildError),
}
