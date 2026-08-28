//! `.baln` reader — canonical implementation lives in `esperanto_bamio::baln`
//! (shared with score/pileup); this shim keeps the scan-side error type.

use crate::error::CallError;
use std::path::Path;

pub use esperanto_bamio::baln::{BalnIndex, BalnReader, BalnRecord};

fn conv(e: std::io::Error) -> CallError {
    CallError::Io {
        path: ".baln".into(),
        source: e,
    }
}

/// Read one record at a byte offset (see bamio::baln::read_record_at).
pub fn read_record_at(file: &std::fs::File, off: u64) -> Result<Option<BalnRecord>, CallError> {
    esperanto_bamio::baln::read_record_at(file, off).map_err(conv)
}

/// Build the coordinate index (see bamio::baln::BalnReader::build_index).
pub fn build_index(path: &Path) -> Result<BalnIndex, CallError> {
    BalnReader::build_index(path).map_err(conv)
}
