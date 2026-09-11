//! Candidate hard-filter for the scan→score stage.
//!
//! Faithful port of the legacy call filter (see docs/DESIGN-DECISIONS.md): the
//! scoring model was trained and validated on a seven-step-filtered candidate
//! set, so the pipeline must apply the same hard filter before scoring.
//! Keeps a candidate when its evidence passes one of two signal arms, subject
//! to a minimum depth. dbSNP-common and homopolymer exclusions are deferred.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::FlowError;

/// Hard-filter rules over `candidates.bed` (10-column contract:
/// `chrom  pos0  pos0+1  strand  evid  call_score  depth  var_freq  fwd_freq  rev_freq`).
#[derive(Debug, Clone)]
pub struct CallFilter {
    /// Keep only sites with depth ≥ this value.
    pub min_depth: u64,
    /// High-score arm: `CALL_SCORE ≥ this` is kept directly (`None` = disabled).
    pub score_arm: Option<f64>,
    /// Recall arm: `var_freq ≥ this` is kept (low-frequency editing sites).
    pub min_vaf: f64,
}

impl Default for CallFilter {
    fn default() -> Self {
        Self {
            min_depth: 10,
            score_arm: Some(0.9),
            min_vaf: 0.05,
        }
    }
}

/// Per-reason drop counts for one filter pass.
#[derive(Debug, Default)]
pub struct FilterStats {
    /// Candidate rows read.
    pub input: usize,
    /// Rows kept.
    pub kept: usize,
    /// Rows dropped for depth < min_depth.
    pub low_depth: usize,
    /// Rows dropped for failing both signal arms.
    pub no_signal: usize,
}

impl CallFilter {
    fn keep(&self, cols: &[&str]) -> bool {
        let depth: u64 = cols.get(6).and_then(|s| s.parse().ok()).unwrap_or(0);
        if depth < self.min_depth {
            return false;
        }
        let score: f64 = cols.get(5).and_then(|s| s.parse().ok()).unwrap_or(0.0);
        let var_freq: f64 = cols.get(7).and_then(|s| s.parse().ok()).unwrap_or(0.0);
        let arm_score = self.score_arm.is_some_and(|t| score >= t);
        let arm_vaf = var_freq >= self.min_vaf;
        arm_score || arm_vaf
    }

    /// Filter a `candidates.bed` in place; returns drop stats.
    pub fn apply_to_bed(&self, bed: &Path) -> Result<FilterStats, FlowError> {
        let input = File::open(bed)?;
        let mut header: Vec<String> = Vec::new();
        let mut kept: Vec<String> = Vec::new();
        let mut st = FilterStats::default();
        for line in BufReader::new(input).lines() {
            let line = line?;
            if line.is_empty() || line.starts_with('#') {
                header.push(line);
                continue;
            }
            st.input += 1;
            let cols: Vec<&str> = line.split('\t').collect();
            let depth: u64 = cols.get(6).and_then(|s| s.parse().ok()).unwrap_or(0);
            if depth < self.min_depth {
                st.low_depth += 1;
                continue;
            }
            if !self.keep(&cols) {
                st.no_signal += 1;
                continue;
            }
            st.kept += 1;
            kept.push(line);
        }
        // Rewrite header + kept rows (scan output is already (chrom, pos)-sorted).
        let out = BufWriter::new(File::create(bed)?);
        let mut w = out;
        for h in header {
            writeln!(w, "{h}")?;
        }
        for k in kept {
            writeln!(w, "{k}")?;
        }
        w.flush()?;
        Ok(st)
    }
}
