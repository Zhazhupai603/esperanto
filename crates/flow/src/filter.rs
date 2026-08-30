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
    /// Recall arm: editing-consistent frequency >= this is kept (low-frequency
    /// editing sites; A>G forward / T>C reverse).
    pub min_vaf: f64,
    /// Recall arm: minimum mutation reads (any non-reference base) required.
    pub min_mutation_reads: u64,
}

impl Default for CallFilter {
    fn default() -> Self {
        Self {
            min_depth: 10,
            score_arm: Some(0.9),
            min_vaf: 0.05,
            min_mutation_reads: 2,
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
        let fwd_freq: f64 = cols.get(8).and_then(|s| s.parse().ok()).unwrap_or(0.0);
        let rev_freq: f64 = cols.get(9).and_then(|s| s.parse().ok()).unwrap_or(0.0);
        let arm_score = self.score_arm.is_some_and(|t| score >= t);
        // Recall arm: editing-consistent signal (A>G forward / T>C reverse)
        // with at least `min_mutation_reads` supporting reads. REF=C/G sites
        // have fwd_freq == rev_freq == 0 and are dropped here.
        let var_reads = var_freq * depth as f64;
        let arm_edit = (fwd_freq >= self.min_vaf || rev_freq >= self.min_vaf)
            && var_reads >= self.min_mutation_reads as f64;
        arm_score || arm_edit
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

#[cfg(test)]
mod tests {
    use super::*;

    fn keep(f: &CallFilter, score: f64, depth: u64, var: f64, fwd: f64, rev: f64) -> bool {
        let line = format!("chr1\t100\t101\tamb\tNONE\t{score}\t{depth}\t{var}\t{fwd}\t{rev}");
        let cols: Vec<&str> = line.split('\t').collect();
        f.keep(&cols)
    }

    #[test]
    fn keeps_editing_direction_site() {
        let f = CallFilter::default();
        // A>G on the forward strand: fwd_freq 0.3, 3 mutation reads at depth 10.
        assert!(keep(&f, 0.5, 10, 0.3, 0.3, 0.0));
    }

    #[test]
    fn keeps_reverse_strand_editing() {
        let f = CallFilter::default();
        // T>C on the reverse strand: rev_freq 0.4.
        assert!(keep(&f, 0.5, 10, 0.4, 0.0, 0.4));
    }

    #[test]
    fn rejects_non_editing_direction() {
        let f = CallFilter::default();
        // REF=C/G: fwd_freq == rev_freq == 0, any-mismatch 0.5.
        assert!(!keep(&f, 0.5, 10, 0.5, 0.0, 0.0));
    }

    #[test]
    fn rejects_single_mutation_read() {
        let f = CallFilter::default();
        // One mutation read at depth 20: var_freq 0.05 -> var_reads 1 < 2.
        assert!(!keep(&f, 0.5, 20, 0.05, 0.05, 0.0));
    }

    #[test]
    fn keeps_high_score_arm() {
        let f = CallFilter::default();
        // The score arm keeps a site regardless of direction or read count.
        assert!(keep(&f, 0.95, 10, 0.0, 0.0, 0.0));
    }

    #[test]
    fn rejects_low_depth() {
        let f = CallFilter::default();
        assert!(!keep(&f, 0.95, 5, 1.0, 1.0, 1.0));
    }
}
