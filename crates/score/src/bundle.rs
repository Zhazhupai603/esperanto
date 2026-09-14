//! Bundle loading for the CMC model bundles (v1.7 architecture):
//! `rust/net_main.safetensors`, `rust/net_hyper.safetensors` and the frozen
//! calibration `rust/threshold.json`.

use crate::cmcnet::CmcNet;
use ndarray::Array2;
use serde::Deserialize;
use std::path::Path;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ScoreError {
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
    #[error("json: {0}")]
    Json(#[from] serde_json::Error),
    #[error("safetensors: {0}")]
    Safe(#[from] safetensors::SafeTensorError),
    #[error("tensor {name} shape mismatch: got {got:?}, want {want:?}")]
    Shape {
        name: String,
        got: Vec<usize>,
        want: Vec<usize>,
    },
    #[error("missing tensor {0}")]
    Missing(String),
    #[error("htslib: {0}")]
    Hts(#[from] rust_htslib::errors::Error),
}

/// Frozen calibration thresholds shipped with the bundle (`rust/threshold.json`).
#[derive(Debug, Deserialize)]
struct ThresholdDoc {
    #[serde(flatten)]
    tasks: std::collections::HashMap<String, ThresholdTask>,
}

#[derive(Debug, Deserialize)]
struct ThresholdTask {
    threshold: f64,
}

/// A CMC model bundle: the main (DM vs RE) and hyper (rescue candidate)
/// networks plus their frozen decision thresholds.
pub struct CmcBundle {
    pub main: CmcNet,
    pub hyper: CmcNet,
    pub main_threshold: f64,
    pub hyper_threshold: f64,
}

/// True when `dir` is the `rust/` directory of a CMC bundle.
pub fn is_cmc_bundle(dir: &Path) -> bool {
    dir.join("net_main.safetensors").is_file()
}

/// Load a CMC bundle from its `rust/` directory.
pub fn load_cmc_bundle(dir: &Path) -> Result<CmcBundle, ScoreError> {
    let load_net = |name: &str| -> Result<CmcNet, ScoreError> {
        let bytes = std::fs::read(dir.join(name))?;
        let st = safetensors::SafeTensors::deserialize(&bytes)?;
        CmcNet::load(&st)
    };
    let doc: ThresholdDoc =
        serde_json::from_str(&std::fs::read_to_string(dir.join("threshold.json"))?)?;
    let task = |name: &str| -> Result<f64, ScoreError> {
        doc.tasks
            .get(name)
            .map(|t| t.threshold)
            .ok_or_else(|| ScoreError::Missing(format!("threshold task {name}")))
    };
    Ok(CmcBundle {
        main: load_net("net_main.safetensors")?,
        hyper: load_net("net_hyper.safetensors")?,
        main_threshold: task("main")?,
        hyper_threshold: task("hyper")?,
    })
}

/// Read only the frozen thresholds of a CMC bundle (no weights loaded).
pub fn load_cmc_thresholds(dir: &Path) -> Result<(f64, f64), ScoreError> {
    let doc: ThresholdDoc =
        serde_json::from_str(&std::fs::read_to_string(dir.join("threshold.json"))?)?;
    let task = |name: &str| -> Result<f64, ScoreError> {
        doc.tasks
            .get(name)
            .map(|t| t.threshold)
            .ok_or_else(|| ScoreError::Missing(format!("threshold task {name}")))
    };
    Ok((task("main")?, task("hyper")?))
}

pub(crate) fn take_matrix(
    st: &safetensors::SafeTensors,
    name: &str,
    rows: usize,
    cols: usize,
) -> Result<Array2<f32>, ScoreError> {
    let t = st
        .tensor(name)
        .map_err(|_| ScoreError::Missing(name.into()))?;
    let got: Vec<usize> = t.shape().to_vec();
    if got != vec![rows, cols] {
        return Err(ScoreError::Shape {
            name: name.into(),
            got,
            want: vec![rows, cols],
        });
    }
    let data: Vec<f32> = t
        .data()
        .chunks_exact(4)
        .map(|b| f32::from_le_bytes(b.try_into().unwrap()))
        .collect();
    Ok(Array2::from_shape_vec((rows, cols), data).unwrap())
}

pub(crate) fn take_vector(
    st: &safetensors::SafeTensors,
    name: &str,
    len: usize,
) -> Result<ndarray::Array1<f32>, ScoreError> {
    let t = st
        .tensor(name)
        .map_err(|_| ScoreError::Missing(name.into()))?;
    let got: Vec<usize> = t.shape().to_vec();
    if got != vec![len] {
        return Err(ScoreError::Shape {
            name: name.into(),
            got,
            want: vec![len],
        });
    }
    let data: Vec<f32> = t
        .data()
        .chunks_exact(4)
        .map(|b| f32::from_le_bytes(b.try_into().unwrap()))
        .collect();
    Ok(ndarray::Array1::from_vec(data))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cmc_bundle_detection_follows_the_network_file() {
        let tmp = std::env::temp_dir().join(format!("esp-cmc-detect-{}", std::process::id()));
        std::fs::create_dir_all(&tmp).unwrap();
        assert!(!is_cmc_bundle(&tmp));
        std::fs::write(tmp.join("net_main.safetensors"), b"x").unwrap();
        assert!(is_cmc_bundle(&tmp));
        std::fs::remove_dir_all(&tmp).unwrap();
    }

    #[test]
    fn thresholds_load_from_json() {
        let tmp = std::env::temp_dir().join(format!("esp-cmc-thr-{}", std::process::id()));
        std::fs::create_dir_all(&tmp).unwrap();
        std::fs::write(
            tmp.join("threshold.json"),
            r#"{"main": {"threshold": 0.25}, "hyper": {"threshold": 0.5}}"#,
        )
        .unwrap();
        let (main, hyper) = load_cmc_thresholds(&tmp).unwrap();
        assert_eq!(main, 0.25);
        assert_eq!(hyper, 0.5);
        std::fs::remove_dir_all(&tmp).unwrap();
    }
}
