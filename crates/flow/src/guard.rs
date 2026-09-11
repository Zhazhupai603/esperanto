//! Species guardrail (spec §species guardrail): runs before any stage starts, so a
//! mismatched reference fails before a single byte of compute is burned.

use std::path::Path;

use crate::FlowError;

/// hg38 chr1 length.
const HG38_CHR1: u64 = 248_956_422;
/// Synthetic/test references are exempt below this length.
const SYNTHETIC_CAP: u64 = 10_000_000;

/// Enforce the species guardrail against `<fasta>.fai`.
///
/// With a `chr1` line: its length must equal the hg38 value or fall below the
/// synthetic cap, otherwise [`FlowError::SpeciesMismatch`]. Without a `chr1`
/// line the check passes (mirror of the score-internal guardrail).
pub fn check_species(fasta: &Path) -> Result<(), FlowError> {
    let fai = format!("{}.fai", fasta.display());
    let text = std::fs::read_to_string(&fai)?;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("chr1\t") {
            let len: u64 = rest
                .split('\t')
                .next()
                .and_then(|v| v.parse().ok())
                .unwrap_or(0);
            if len == HG38_CHR1 || len < SYNTHETIC_CAP {
                return Ok(());
            }
            return Err(FlowError::SpeciesMismatch { len });
        }
    }
    Ok(())
}
