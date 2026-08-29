//! Species guardrail (spec §species guardrail): runs before any stage starts, so a
//! mismatched reference fails before a single byte of compute is burned.
//! With a `species.json` manifest the check follows the declared baseline;
//! without one it falls back to the hg38 chr1 heuristic.

use std::path::Path;

use crate::manifest::SpeciesManifest;
use crate::FlowError;

/// hg38 chr1 length.
const HG38_CHR1: u64 = 248_956_422;
/// GRCm39 chr1 length (mouse baseline for hybrid references).
const GRCM39_CHR1: u64 = 195_154_007;
/// mm10 chr1 length (user-placed legacy mouse references).
const MM10_CHR1: u64 = 195_471_971;
/// Synthetic/test references are exempt below this length.
const SYNTHETIC_CAP: u64 = 10_000_000;

/// Enforce the species guardrail against `<fasta>.fai`.
///
/// With a manifest: the baseline tag decides the expected chr1 length
/// (`grch38` → hg38, `grcm39`/`mm10` → the mouse value); `hybrid` checks the
/// mouse baseline chr1 (human loci are verified at build time).
/// Without a manifest (legacy): a `chr1` line must equal the hg38 value or
/// fall below the synthetic cap, otherwise [`FlowError::SpeciesMismatch`];
/// no `chr1` line passes directly.
pub fn check_species(fasta: &Path, manifest: Option<&SpeciesManifest>) -> Result<(), FlowError> {
    let fai = format!("{}.fai", fasta.display());
    let text = std::fs::read_to_string(&fai)?;
    let chr1_len = text.lines().find_map(|line| {
        line.strip_prefix("chr1	")
            .and_then(|rest| rest.split('\t').next())
            .and_then(|v| v.parse::<u64>().ok())
    });
    if let Some(m) = manifest {
        let expected = match m.baseline.as_str() {
            "grch38" => HG38_CHR1,
            "grcm39" => GRCM39_CHR1,
            "mm10" => MM10_CHR1,
            other => {
                return Err(FlowError::Entry(format!(
                    "species.json declares unknown baseline '{other}' (expected grch38 / grcm39 / mm10)"
                )))
            }
        };
        match chr1_len {
            Some(len) if len == expected || len < SYNTHETIC_CAP => Ok(()),
            Some(len) => Err(FlowError::SpeciesMismatch { len }),
            None => Ok(()),
        }
    } else {
        match chr1_len {
            Some(len) if len == HG38_CHR1 || len < SYNTHETIC_CAP => Ok(()),
            Some(len) => Err(FlowError::SpeciesMismatch { len }),
            None => Ok(()),
        }
    }
}
