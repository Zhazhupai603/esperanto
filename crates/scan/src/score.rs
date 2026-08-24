//! Soft score CALL_SCORE (0–1): logistic-style weighted sum; weights come from a versioned spec JSON.
//! To be clear: this is not the v05 model, only a candidate-ranking score (DESIGN §1.2c).

use crate::error::CallError;
use serde::Deserialize;
use std::path::Path;

/// Crate-bundled default spec: v2 = real-data recalibration (W1, bench/W1_CALIBRATION.md),
/// the default production profile for candidate ranking (real_v2); the v1 synthetic calibration is kept in
/// call_spec.v1.json (profile synthetic_v1, pinned by the synthetic tests).
/// --spec can override with any profile.
pub const DEFAULT_SPEC: &str = include_str!("../call_spec.v2.json");

/// v1 synthetic-data calibration (inflated on real BAMs).
pub const V1_SPEC: &str = include_str!("../call_spec.v1.json");

#[derive(Debug, Clone, Deserialize)]
pub struct CallSpec {
    pub version: u32,
    /// Calibration-environment profile: real_v2 (real BAM, default) / synthetic_v1 (synthetic calibration).
    /// v1 files lack this field → None (i.e. synthetic_v1).
    #[serde(default)]
    pub profile: Option<String>,
    pub intercept: f64,
    pub weights: Weights,
}

#[derive(Debug, Clone, Deserialize)]
pub struct Weights {
    pub depth: f64,
    pub var_freq: f64,
    pub mean_bq: f64,
    pub mean_mapq: f64,
    pub strand_bias: f64,
    pub gnomad_af: f64,
    pub homopolymer: f64,
    pub junction: f64,
}

/// Scoring inputs for a single candidate site (all also written to BED columns or evid; nothing is removed).
#[derive(Debug, Clone, Copy)]
pub struct ScoreFeatures {
    pub depth: u64,
    /// Strand-resolved edit-consistent variant frequency (plus strand A>G / minus strand T>C; amb takes the max of both; includes C>U when enable_cu).
    pub edit_frac: f64,
    pub mean_bq: f64,
    pub mean_mapq: f64,
    /// Strand concentration of variant reads |var_fwd - var_rev| / (var_fwd + var_rev), 0–1.
    pub strand_bias: f64,
    /// gnomAD AF (soft down-weighting, never removes).
    pub gnomad_af: Option<f64>,
    /// Homopolymer context length (same-base run containing the site itself; 0 without fasta).
    pub hp_len: u32,
    /// Distance (bp) to the nearest junction boundary; None when no junction evidence.
    pub junction_dist: Option<u32>,
}

impl CallSpec {
    pub fn load(path: Option<&Path>) -> Result<CallSpec, CallError> {
        let (loc, text) = match path {
            Some(p) => {
                let t = std::fs::read_to_string(p).map_err(|e| CallError::Io {
                    path: p.display().to_string(),
                    source: e,
                })?;
                (p.display().to_string(), t)
            }
            None => ("<builtin>".to_string(), DEFAULT_SPEC.to_string()),
        };
        serde_json::from_str(&text).map_err(|e| CallError::Spec {
            path: loc,
            msg: e.to_string(),
        })
    }

    /// Logistic-style weighted sum → 0–1.
    pub fn score(&self, f: &ScoreFeatures) -> f64 {
        let w = &self.weights;
        let depth_x = (f.depth.min(50) as f64) / 50.0;
        let bq_x = ((f.mean_bq - 10.0) / 30.0).clamp(0.0, 1.0);
        let mapq_x = (f.mean_mapq / 60.0).clamp(0.0, 1.0);
        // AF 5% → 1.0: common SNPs strongly down-weighted; missing AF = 0 term (neutral).
        let gnomad_x = f.gnomad_af.map(|af| (af * 20.0).min(1.0)).unwrap_or(0.0);
        let hp_x = (f.hp_len.min(10) as f64) / 10.0;
        let junc_x = match f.junction_dist {
            Some(d) if d <= 4 => (5 - d) as f64 / 5.0,
            _ => 0.0,
        };
        let z = self.intercept
            + w.depth * depth_x
            + w.var_freq * f.edit_frac
            + w.mean_bq * bq_x
            + w.mean_mapq * mapq_x
            + w.strand_bias * f.strand_bias
            + w.gnomad_af * gnomad_x
            + w.homopolymer * hp_x
            + w.junction * junc_x;
        sigmoid(z)
    }
}

fn sigmoid(z: f64) -> f64 {
    1.0 / (1.0 + (-z.clamp(-50.0, 50.0)).exp())
}
