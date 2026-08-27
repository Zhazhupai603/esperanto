//! Bundle loading: safetensors heads + norm.json (artifact contract).

use crate::head::FoldHead;
use ndarray::Array2;
use serde::Deserialize;
use std::collections::HashMap;
use std::path::Path;
use thiserror::Error;

pub const N_FOLDS: usize = 5;
pub const PILEUP_DIM: usize = 8;

#[derive(Debug, Error)]
pub enum ScoreError {
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
    #[error("json: {0}")]
    Json(#[from] serde_json::Error),
    #[error("safetensors: {0}")]
    Safe(#[from] safetensors::SafeTensorError),
    #[error("tensor {name} shape mismatch: got {got:?}, want {want:?}")]
    Shape { name: String, got: Vec<usize>, want: Vec<usize> },
    #[error("missing tensor {0}")]
    Missing(String),
    #[error("htslib: {0}")]
    Hts(#[from] rust_htslib::errors::Error),
    #[error("token id {0} out of range [0,16)")]
    Token(i64),
    #[cfg(feature = "gpu")]
    #[error("gpu: {0}")]
    Gpu(#[from] candle_core::Error),
}

#[derive(Debug, Clone, Copy)]
pub struct NormStats {
    pub mean: [f64; PILEUP_DIM],
    pub std: [f64; PILEUP_DIM],
}

pub struct Bundle {
pub heads: [FoldHead; N_FOLDS],
pub norms: [NormStats; N_FOLDS],
/// Half-width of the encoding window (bp; feature_spec.sequence.half_window).
/// v1.0 = 500 (1001bp); v1.2 = 250 (501bp). Missing spec file = 500 (backwards compatible).
    pub half_window: i64,
    /// v1.3: pileup-only veto gate (Some only when all three components are present; partial presence = corrupted-bundle error).
    /// Required by the score pipeline (mandatory since v1.3, no switch).
    pub gate: Option<Gate>,
    /// v1.4.2: embedding cache identity (hash of feature_spec + window; cache files are checked for compatibility against it).
    pub cache_id: u64,
}

/// v1.3 veto gate: pileup-only 5-fold ensemble (= the ablate_seq EsperantoFusionHead with zero embedding).
pub struct Gate {
    pub heads: [FoldHead; N_FOLDS],
    pub norms: [NormStats; N_FOLDS],
    /// Veto threshold (feature_spec.gate.threshold): gate RE_PROB < threshold -> skip the encoder.
    pub threshold: f64,
}

#[derive(Deserialize)]
struct NormFold {
    mean: Vec<f64>,
    std: Vec<f64>,
}

/// `dir` = bundle/esperanto-model-v1.0.0/rust (contains heads/, norm.json)
pub fn load_bundle(dir: &Path) -> Result<Bundle, ScoreError> {
let norm_doc: HashMap<String, NormFold> =
        serde_json::from_str(&std::fs::read_to_string(dir.join("norm.json"))?)?;
    // v1.2: window half-width comes from feature_spec.sequence.half_window (v1.0 has no such file = 500).
    // The spec is looked up in two places: dir (= bundle/rust, embedded in v1.2 exports) and dir.parent
    // (v1.0 layout, spec at the bundle root). Observed failure: when only dir was checked, a v1.2 root-level
    // spec was missed and half_window silently fell back to 500 -> a 501 model scored with a 1001bp window
    // (no speedup + window mismatch).
    let read_half = |p: &Path| -> Option<i64> {
        std::fs::read_to_string(p)
            .ok()
            .and_then(|t| serde_json::from_str::<serde_json::Value>(&t).ok())
            .and_then(|v| {
                v.get("sequence")
                    .and_then(|sq| sq.get("half_window"))
                    .and_then(|h| h.as_i64())
            })
    };
    let spec_text = std::fs::read_to_string(dir.join("feature_spec.json")).unwrap_or_else(|_| {
        dir.parent()
            .and_then(|p| std::fs::read_to_string(p.join("feature_spec.json")).ok())
            .unwrap_or_default()
    });
    let half_window = read_half(&dir.join("feature_spec.json"))
        .or_else(|| dir.parent().and_then(|p| read_half(&p.join("feature_spec.json"))))
        .unwrap_or(500);
    // v1.4.2: cache identity = hash of full feature_spec text + window (any change to encoder/window/gate -> id changes).
    let cache_id = {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};
        let mut h = DefaultHasher::new();
        spec_text.hash(&mut h);
        half_window.hash(&mut h);
        h.finish()
    };

    // v1.3: gate three-piece set (gate_heads/ + gate_norm.json + feature_spec.gate.threshold).
    // All missing = old bundle (None); partially missing = corrupted (hard error, no guessing).
    let read_gate_threshold = |p: &Path| -> Option<f64> {
        std::fs::read_to_string(p)
            .ok()
            .and_then(|t| serde_json::from_str::<serde_json::Value>(&t).ok())
            .and_then(|v| v.get("gate")?.get("threshold")?.as_f64())
    };
    let gate_threshold = read_gate_threshold(&dir.join("feature_spec.json"))
        .or_else(|| dir.parent().and_then(|p| read_gate_threshold(&p.join("feature_spec.json"))));
    let gate_dir = dir.join("gate_heads");
    let gate_norm_path = dir.join("gate_norm.json");
    let n_present = [gate_dir.exists(), gate_norm_path.exists(), gate_threshold.is_some()]
        .iter()
        .filter(|b| **b)
        .count();
    let gate = if n_present == 0 {
        None
    } else {
        if n_present < 3 {
            return Err(ScoreError::Missing(format!(
                "incomplete gate triple (gate_heads={} gate_norm={} threshold={}) — corrupted bundle",
                gate_dir.exists(),
                gate_norm_path.exists(),
                gate_threshold.is_some()
            )));
        }
        let gate_norm_doc: HashMap<String, NormFold> =
            serde_json::from_str(&std::fs::read_to_string(&gate_norm_path)?)?;
        let mut g_heads: Vec<FoldHead> = Vec::with_capacity(N_FOLDS);
        let mut g_norms: Vec<NormStats> = Vec::with_capacity(N_FOLDS);
        for fold in 0..N_FOLDS {
            let bytes = std::fs::read(gate_dir.join(format!("fold_{fold}.safetensors")))?;
            let st = safetensors::SafeTensors::deserialize(&bytes)?;
            g_heads.push(FoldHead::from_tensors(&st)?);
            let nf = gate_norm_doc
                .get(&format!("fold_{fold}"))
                .ok_or_else(|| ScoreError::Missing(format!("gate norm fold_{fold}")))?;
            g_norms.push(NormStats {
                mean: nf.mean.clone().try_into().map_err(|_| ScoreError::Shape {
                    name: "gate_norm.mean".into(), got: vec![nf.mean.len()], want: vec![8],
                })?,
                std: nf.std.clone().try_into().map_err(|_| ScoreError::Shape {
                    name: "gate_norm.std".into(), got: vec![nf.std.len()], want: vec![8],
                })?,
            });
        }
        Some(Gate {
            heads: g_heads.try_into().map_err(|_| ScoreError::Missing("gate heads<5".into()))?,
            norms: g_norms.try_into().map_err(|_| ScoreError::Missing("gate norms<5".into()))?,
            threshold: gate_threshold.unwrap(),
        })
    };


    let mut heads: Vec<FoldHead> = Vec::with_capacity(N_FOLDS);
    let mut norms: Vec<NormStats> = Vec::with_capacity(N_FOLDS);
    for fold in 0..N_FOLDS {
        let bytes = std::fs::read(dir.join("heads").join(format!("fold_{fold}.safetensors")))?;
        let st = safetensors::SafeTensors::deserialize(&bytes)?;
        heads.push(FoldHead::from_tensors(&st)?);

        let nf = norm_doc
            .get(&format!("fold_{fold}"))
            .ok_or_else(|| ScoreError::Missing(format!("norm fold_{fold}")))?;
        norms.push(NormStats {
            mean: nf.mean.clone().try_into().map_err(|_| ScoreError::Shape {
                name: "norm.mean".into(), got: vec![nf.mean.len()], want: vec![8],
            })?,
            std: nf.std.clone().try_into().map_err(|_| ScoreError::Shape {
                name: "norm.std".into(), got: vec![nf.std.len()], want: vec![8],
            })?,
        });
    }
Ok(Bundle {
heads: heads.try_into().map_err(|_| ScoreError::Missing("heads<5".into()))?,
norms: norms.try_into().map_err(|_| ScoreError::Missing("norms<5".into()))?,
        half_window,
        gate,
        cache_id,
})
}

pub(crate) fn take_matrix(
    st: &safetensors::SafeTensors,
    name: &str,
    rows: usize,
    cols: usize,
) -> Result<Array2<f32>, ScoreError> {
    let t = st.tensor(name).map_err(|_| ScoreError::Missing(name.into()))?;
    let got: Vec<usize> = t.shape().to_vec();
    if got != vec![rows, cols] {
        return Err(ScoreError::Shape {
            name: name.into(), got, want: vec![rows, cols],
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
    let t = st.tensor(name).map_err(|_| ScoreError::Missing(name.into()))?;
    let got: Vec<usize> = t.shape().to_vec();
    if got != vec![len] {
        return Err(ScoreError::Shape {
            name: name.into(), got, want: vec![len],
        });
    }
    let data: Vec<f32> = t
        .data()
        .chunks_exact(4)
        .map(|b| f32::from_le_bytes(b.try_into().unwrap()))
        .collect();
    Ok(ndarray::Array1::from_vec(data))
}

// ===================== Reference embedding cache =====================
// Reference-genome embedding lookup table (the encoder only reads the reference sequence -> site
// embeddings are sample-independent -> permanently reusable across runs and samples). Industry
// analogues: Ensembl VEP SpliceAI precomputed lookup, Sei whole-genome annotations (Zenodo).
//
// Bit-equivalence design: production embed output is fp16 (feature_spec embedding_dtype); the cache
// stores the production-grade fp16 bit patterns -- a hit returns a value bit-identical to online
// computation (experiments show an embedding perturbation of 1e-2 can flip a probability of 0.69 at
// the decision boundary, so this cache applies no compression/quantization).
//
// File format (LE):
//   magic "ESPEMBC1" (8B)
//   count u64 -- deprecated; count is derived from file length; reserved field, always 0 (forward compatible)
//   records: { chrom_len u16 | chrom bytes | pos_1based u32 | emb 118 x f16 bits u16 }
// Records are appended (append-only); at startup the whole table is loaded into HashMap((String,u32) -> [u16;118]).
// 118x2B + ~12B key ~= 250B/site -> 100M sites ~= 25GB (upper bound; generated on demand in practice).

/// In-process embedding cache (table + target file path; explicit flush before drop).
pub struct EmbCache {
    map: std::collections::HashMap<(String, u32), [u16; crate::caduceus::D_MODEL]>,
    path: std::path::PathBuf,
    dirty: usize,
    cache_id: u64,
    half_window: i64,
}

impl EmbCache {
    /// Open (or create) the cache file; if the file does not exist -> empty table, new path.
    /// v1.4.2: the header carries the model identity (cache_id) + half_window -- incompatible (old format /
    /// different model / different window) -> warn and reopen as an empty table, overwriting on flush.
    /// The cache is a pure optimization: semantically, prefer recomputation over ever reading a stale embedding.
    pub fn open(path: &std::path::Path, cache_id: u64, half_window: i64) -> std::io::Result<Self> {
        let mut map = std::collections::HashMap::new();
        if path.exists() {
            let data = std::fs::read(path)?;
            let compatible = data.len() >= 24
                && &data[..8] == b"ESPEMBC2"
                && u64::from_le_bytes(data[8..16].try_into().unwrap()) == cache_id
                && i64::from_le_bytes(data[16..24].try_into().unwrap()) == half_window;
            if compatible {
                let mut off = 24usize;
                while off + 2 <= data.len() {
                    let cl = u16::from_le_bytes([data[off], data[off + 1]]) as usize;
                    off += 2;
                    if off + cl + 4 + crate::caduceus::D_MODEL * 2 > data.len() {
                        break; // trailing truncated record (interrupted append) -- ignore
                    }
                    let chrom = String::from_utf8_lossy(&data[off..off + cl]).into_owned();
                    off += cl;
                    let pos = u32::from_le_bytes([
                        data[off], data[off + 1], data[off + 2], data[off + 3],
                    ]);
                    off += 4;
                    let mut emb = [0u16; crate::caduceus::D_MODEL];
                    for (k, e) in emb.iter_mut().enumerate() {
                        *e = u16::from_le_bytes([data[off + k * 2], data[off + k * 2 + 1]]);
                    }
                    off += crate::caduceus::D_MODEL * 2;
                    map.insert((chrom, pos), emb);
                }
            } else if !data.is_empty() {
                eprintln!(
                    "[score] embedding cache incompatible with current model (old format / different model / different window), reopening empty: {}",
                    path.display()
                );
            }
        }
        Ok(EmbCache {
            map,
            path: path.to_path_buf(),
            dirty: 0,
            cache_id,
            half_window,
        })
    }

    pub fn len(&self) -> usize {
        self.map.len()
    }

    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// Lookup: a hit returns the fp16 bit pattern (bit-identical to online embed output).
    pub fn get(&self, chrom: &str, pos_1based: u32) -> Option<[u16; crate::caduceus::D_MODEL]> {
        self.map.get(&(chrom.to_string(), pos_1based)).copied()
    }

    /// Insert (into memory + mark dirty; flush persists). An existing key = redundant recomputation; overwritten with the same value.
    pub fn put(&mut self, chrom: &str, pos_1based: u32, emb_bits: &[u16; crate::caduceus::D_MODEL]) {
        self.map
            .insert((chrom.to_string(), pos_1based), *emb_bits);
        self.dirty += 1;
    }

    /// Persist (dirty records already entered memory via put; this implementation rewrites the whole table --
    /// simple and atomically safe: write a temp file then rename, avoiding corruption from an interrupted append).
    pub fn flush(&mut self) -> std::io::Result<()> {
        use std::io::Write;
        if self.dirty == 0 {
            return Ok(());
        }
        let tmp = self.path.with_extension("tmp");
        {
            let mut w = std::io::BufWriter::new(std::fs::File::create(&tmp)?);
            w.write_all(b"ESPEMBC2")?;
            w.write_all(&self.cache_id.to_le_bytes())?;
            w.write_all(&self.half_window.to_le_bytes())?;
            for ((chrom, pos), emb) in &self.map {
                let cb = chrom.as_bytes();
                w.write_all(&(cb.len() as u16).to_le_bytes())?;
                w.write_all(cb)?;
                w.write_all(&pos.to_le_bytes())?;
                for v in emb.iter() {
                    w.write_all(&v.to_le_bytes())?;
                }
            }
            w.flush()?;
        }
        std::fs::rename(&tmp, &self.path)?;
        self.dirty = 0;
        Ok(())
    }
}
