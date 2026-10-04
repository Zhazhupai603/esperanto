//! Human hg38 reference validation before pipeline execution.

use crate::FlowError;
use std::path::Path;

const HG38_CHR1: u64 = 248_956_422;
const SYNTHETIC_CAP: u64 = 10_000_000;

/// Require hg38 when chr1 is present; small synthetic references are accepted.
pub fn check_species(fasta: &Path) -> Result<(), FlowError> {
    let text = std::fs::read_to_string(format!("{}.fai", fasta.display()))?;
    let chr1_len = text.lines().find_map(|line| {
        line.strip_prefix("chr1\t")
            .and_then(|rest| rest.split('\t').next())
            .and_then(|v| v.parse::<u64>().ok())
    });
    match chr1_len {
        Some(len) if len == HG38_CHR1 || len < SYNTHETIC_CAP => Ok(()),
        Some(len) => Err(FlowError::SpeciesMismatch { len }),
        None => Ok(()),
    }
}
