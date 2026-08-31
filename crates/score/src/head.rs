//! EsperantoFusionHead forward (feature_spec.json v1):
//! pileup MLP(8->32->32, ReLU); concat(emb118, 32)=150; MLP(150->128->2);
//! softmax[:,1] = RE_PROB. Ensemble = 5-fold mean.

use crate::bundle::{take_matrix, take_vector, NormStats, ScoreError};
use ndarray::{Array1, Array2, ArrayView1};

pub const CADUCEUS_DIM: usize = 118;
pub const PILEUP_HIDDEN: usize = 32;
pub const HEAD_HIDDEN: usize = 128;

pub struct FoldHead {
    pe0_w: Array2<f32>, pe0_b: Array1<f32>, // (32,9) (32,)
    pe2_w: Array2<f32>, pe2_b: Array1<f32>, // (32,32)
    h0_w: Array2<f32>, h0_b: Array1<f32>,   // (128,150)
    h2_w: Array2<f32>, h2_b: Array1<f32>,   // (2,128)
}

impl FoldHead {
    pub fn from_tensors(st: &safetensors::SafeTensors) -> Result<Self, ScoreError> {
        Ok(Self {
            pe0_w: take_matrix(st, "pileup_encoder.net.0.weight", 32, 9)?,
            pe0_b: take_vector(st, "pileup_encoder.net.0.bias", 32)?,
            pe2_w: take_matrix(st, "pileup_encoder.net.2.weight", 32, 32)?,
            pe2_b: take_vector(st, "pileup_encoder.net.2.bias", 32)?,
            h0_w: take_matrix(st, "head.0.weight", 128, 150)?,
            h0_b: take_vector(st, "head.0.bias", 128)?,
            h2_w: take_matrix(st, "head.2.weight", 2, 128)?,
            h2_b: take_vector(st, "head.2.bias", 2)?,
        })
    }
}

fn relu(v: &mut Array1<f32>) {
    v.mapv_inplace(|x| if x > 0.0 { x } else { 0.0 });
}

/// Single-fold RE_PROB. pileup is the raw 9-dim feature (z-score applied internally).
/// Returns ScoreError when emb length is not CADUCEUS_DIM (was a debug_assert, ineffective in
/// release; promoted to a hard error, effective in all build modes).
pub fn re_prob_fold(
    head: &FoldHead,
    norm: &NormStats,
    emb: &ArrayView1<f32>,
    pileup: &[f32; 9],
) -> Result<f64, ScoreError> {
    if emb.len() != CADUCEUS_DIM {
        return Err(ScoreError::Shape {
            name: "emb".into(),
            got: vec![emb.len()],
            want: vec![CADUCEUS_DIM],
        });
    }

    // per-fold z-score (feature_spec normalization)
    let z: Vec<f32> = (0..9)
        .map(|i| ((pileup[i] as f64 - norm.mean[i]) / norm.std[i]) as f32)
        .collect();
    let z = Array1::from_vec(z);

    // pileup MLP 8->32->32 (ReLU). PyTorch Linear: y = x @ W^T + b
    let mut h = head.pe0_w.dot(&z) + &head.pe0_b;
    relu(&mut h);
    let mut h = head.pe2_w.dot(&h) + &head.pe2_b;
    relu(&mut h);

    // concat(emb118, pileup32) = 150 -> MLP 150->128->2
    let mut x = Array1::<f32>::zeros(150);
    x.slice_mut(ndarray::s![..118]).assign(emb);
    x.slice_mut(ndarray::s![118..]).assign(&h);
    let mut h = head.h0_w.dot(&x) + &head.h0_b;
    relu(&mut h);
    let logits = head.h2_w.dot(&h) + &head.h2_b; // (2,)

    // softmax (fp64 accumulation to match torch behavior)
    let m = logits[0].max(logits[1]) as f64;
    let e0 = ((logits[0] as f64) - m).exp();
    let e1 = ((logits[1] as f64) - m).exp();
    Ok(e1 / (e0 + e1))
}

/// 5-fold ensemble mean RE_PROB.
pub fn re_prob_ensemble(
bundle: &crate::bundle::Bundle,
emb: &ArrayView1<f32>,
pileup: &[f32; 9],
) -> Result<f64, ScoreError> {
let sum = (0..5)
.map(|f| re_prob_fold(&bundle.heads[f], &bundle.norms[f], emb, pileup))
.sum::<Result<f64, ScoreError>>()?;
Ok(sum / 5.0)
}

/// v1.3 veto-gate probability: 5-fold ensemble, zero embedding (isomorphic to the v02 pileup_only model step by step).
pub fn gate_prob_ensemble(
    gate: &crate::bundle::Gate,
    pileup: &[f32; 9],
) -> Result<f64, ScoreError> {
    let zero = Array1::<f32>::zeros(CADUCEUS_DIM);
    let sum = (0..5)
        .map(|f| re_prob_fold(&gate.heads[f], &gate.norms[f], &zero.view(), pileup))
        .sum::<Result<f64, ScoreError>>()?;
    Ok(sum / 5.0)
}
