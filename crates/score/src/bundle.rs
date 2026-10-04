//! Loading and validation for the unified reference-aware model bundle.

use std::path::Path;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ScoreError {
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
    #[error("json: {0}")]
    Json(#[from] serde_json::Error),
    #[error("missing model data: {0}")]
    Missing(String),
}

/// True when `dir` contains the native reference-aware unified model.
pub fn is_unified_bundle(dir: &Path) -> bool {
    dir.join("model.txt").is_file()
        && dir.join("freeze.json").is_file()
        && dir.join("alu.bed").is_file()
}

/// Load the frozen decision threshold from the unified bundle.
pub fn load_primary_threshold(dir: &Path) -> Result<f64, ScoreError> {
    if !is_unified_bundle(dir) {
        return Err(ScoreError::Missing(
            "unified model bundle (model.txt, freeze.json, alu.bed)".into(),
        ));
    }
    let text = std::fs::read_to_string(dir.join("freeze.json"))?;
    let doc: serde_json::Value = serde_json::from_str(&text)?;
    doc.get("threshold")
        .and_then(|value| value.as_f64())
        .ok_or_else(|| ScoreError::Missing("threshold".into()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unified_bundle_detection_requires_all_frozen_files() {
        let tmp = std::env::temp_dir().join(format!("esp-unified-detect-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&tmp);
        std::fs::create_dir_all(&tmp).unwrap();
        assert!(!is_unified_bundle(&tmp));
        std::fs::write(tmp.join("model.txt"), b"model").unwrap();
        std::fs::write(tmp.join("freeze.json"), b"{}").unwrap();
        assert!(!is_unified_bundle(&tmp));
        std::fs::write(tmp.join("alu.bed"), b"").unwrap();
        assert!(is_unified_bundle(&tmp));
        std::fs::remove_dir_all(&tmp).unwrap();
    }

    #[test]
    fn unified_threshold_loads_from_freeze_json() {
        let tmp = std::env::temp_dir().join(format!("esp-unified-thr-{}", std::process::id()));
        std::fs::create_dir_all(&tmp).unwrap();
        std::fs::write(tmp.join("model.txt"), b"model").unwrap();
        std::fs::write(tmp.join("alu.bed"), b"").unwrap();
        std::fs::write(tmp.join("freeze.json"), r#"{"threshold": 0.25}"#).unwrap();
        assert_eq!(load_primary_threshold(&tmp).unwrap(), 0.25);
        std::fs::remove_dir_all(&tmp).unwrap();
    }
}
