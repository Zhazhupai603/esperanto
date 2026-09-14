//! CMC scoring network forward pass (bundle v1.7 architecture).
//!
//! Frozen architecture (847,042 parameters, identical for the main and the
//! hyper net):
//!
//! ```text
//! (16,50) -> conv1d(64,k5,same) -> batchnorm -> gelu
//!         -> conv1d(128,k3,same) -> batchnorm -> gelu
//!         -> transpose to (50 tokens,128)
//!         -> token embedding (learned cls + position, d=128)
//!         -> transformer encoder x4 (8 heads, ff 512, pre-norm, final LN)
//!         -> cls pooling -> linear(128) -> gelu -> linear(2) -> softmax
//! ```
//!
//! Numerics: f32 storage and arithmetic; GELU uses the exact erf form
//! (evaluated in f64); softmax accumulates in f64. Deterministic — the
//! forward is a fixed sequence of dense ops with no data-dependent dispatch.

use crate::bundle::{take_matrix, take_vector, ScoreError};
use crate::cmc::{CMC_CHANNELS, CMC_COLUMNS};
use ndarray::{Array1, Array2};

const CONV1_OUT: usize = 64;
const CONV1_K: usize = 5;
const CONV2_OUT: usize = 128;
const CONV2_K: usize = 3;
const D_MODEL: usize = 128;
const N_TOKENS: usize = 50 + 1; // flank tokens + cls
const N_HEADS: usize = 8;
const HEAD_DIM: usize = D_MODEL / N_HEADS;
const N_LAYERS: usize = 4;
const FF_DIM: usize = 512;
const LN_EPS: f32 = 1e-5;
const BN_EPS: f32 = 1e-5;

struct Conv1d {
    // (out, in, k) row-major
    w: Vec<f32>,
    b: Vec<f32>,
    out_ch: usize,
    in_ch: usize,
    k: usize,
}

impl Conv1d {
    fn load(
        st: &safetensors::SafeTensors,
        prefix: &str,
        out_ch: usize,
        in_ch: usize,
        k: usize,
    ) -> Result<Self, ScoreError> {
        let w = take_3d(st, &format!("{prefix}.weight"), out_ch, in_ch, k)?;
        let b = take_vector(st, &format!("{prefix}.bias"), out_ch)?;
        Ok(Self {
            w,
            b: b.to_vec(),
            out_ch,
            in_ch,
            k,
        })
    }

    /// "same" padding (k odd, stride 1), zero-padded input of length `len`.
    fn forward_same(&self, x: &[f32], len: usize) -> Array2<f32> {
        let pad = self.k / 2;
        let mut out = Array2::<f32>::zeros((self.out_ch, len));
        for oc in 0..self.out_ch {
            let mut row = out.row_mut(oc);
            let wbase = oc * self.in_ch * self.k;
            let mut acc = vec![0.0f32; len];
            for ic in 0..self.in_ch {
                for (l, slot) in acc.iter_mut().enumerate() {
                    let mut s = 0.0f32;
                    for kk in 0..self.k {
                        let src = l + kk;
                        // zero padding on both sides
                        if src >= pad && src < len + pad {
                            s += self.w[wbase + ic * self.k + kk] * x[ic * len + src - pad];
                        }
                    }
                    *slot += s;
                }
            }
            for v in acc.iter_mut() {
                *v += self.b[oc];
            }
            row.assign(&Array1::from(acc));
        }
        out
    }
}

struct BatchNorm {
    scale: Vec<f32>,
    shift: Vec<f32>,
}

impl BatchNorm {
    fn load(st: &safetensors::SafeTensors, prefix: &str, ch: usize) -> Result<Self, ScoreError> {
        let w = take_vector(st, &format!("{prefix}.weight"), ch)?;
        let b = take_vector(st, &format!("{prefix}.bias"), ch)?;
        let mean = take_vector(st, &format!("{prefix}.running_mean"), ch)?;
        let var = take_vector(st, &format!("{prefix}.running_var"), ch)?;
        let mut scale = vec![0.0f32; ch];
        let mut shift = vec![0.0f32; ch];
        for c in 0..ch {
            let inv = 1.0f32 / (var[c] + BN_EPS).sqrt();
            scale[c] = w[c] * inv;
            shift[c] = b[c] - scale[c] * mean[c];
        }
        Ok(Self { scale, shift })
    }

    fn apply(&self, x: &mut Array2<f32>) {
        for c in 0..x.nrows() {
            x.row_mut(c).mapv_inplace(|v| v * self.scale[c] + self.shift[c]);
        }
    }
}

struct Linear {
    w: Array2<f32>, // (out, in)
    b: Array1<f32>,
}

impl Linear {
    fn load(
        st: &safetensors::SafeTensors,
        prefix: &str,
        out: usize,
        inp: usize,
    ) -> Result<Self, ScoreError> {
        Ok(Self {
            w: take_matrix(st, &format!("{prefix}.weight"), out, inp)?,
            b: take_vector(st, &format!("{prefix}.bias"), out)?,
        })
    }

    fn forward(&self, x: &Array1<f32>) -> Array1<f32> {
        self.w.dot(x) + &self.b
    }
}

struct LayerNorm {
    w: Array1<f32>,
    b: Array1<f32>,
}

impl LayerNorm {
    fn load(st: &safetensors::SafeTensors, prefix: &str, d: usize) -> Result<Self, ScoreError> {
        Ok(Self {
            w: take_vector(st, &format!("{prefix}.weight"), d)?,
            b: take_vector(st, &format!("{prefix}.bias"), d)?,
        })
    }

    /// Row-wise layer norm over the last axis (d = D_MODEL), f64 accumulation.
    fn norm(&self, x: &Array1<f32>) -> Array1<f32> {
        let n = x.len() as f64;
        let mean = x.iter().map(|v| f64::from(*v)).sum::<f64>() / n;
        let var = x
            .iter()
            .map(|v| (f64::from(*v) - mean).powi(2))
            .sum::<f64>()
            / n;
        let inv = 1.0 / (var + f64::from(LN_EPS)).sqrt();
        x.iter()
            .zip(&self.w)
            .zip(&self.b)
            .map(|((&v, &w), &b)| {
                (((f64::from(v) - mean) * inv) * f64::from(w) + f64::from(b)) as f32
            })
            .collect()
    }
}

struct EncoderLayer {
    qkv_w: Array2<f32>, // (3*D, D)
    qkv_b: Array1<f32>,
    out: Linear,
    ff1: Linear,
    ff2: Linear,
    ln1: LayerNorm,
    ln2: LayerNorm,
}

impl EncoderLayer {
    fn load(st: &safetensors::SafeTensors, prefix: &str) -> Result<Self, ScoreError> {
        Ok(Self {
            qkv_w: take_matrix(
                st,
                &format!("{prefix}.self_attn.in_proj_weight"),
                3 * D_MODEL,
                D_MODEL,
            )?,
            qkv_b: take_vector(st, &format!("{prefix}.self_attn.in_proj_bias"), 3 * D_MODEL)?,
            out: Linear::load(
                st,
                &format!("{prefix}.self_attn.out_proj"),
                D_MODEL,
                D_MODEL,
            )?,
            ff1: Linear::load(st, &format!("{prefix}.linear1"), FF_DIM, D_MODEL)?,
            ff2: Linear::load(st, &format!("{prefix}.linear2"), D_MODEL, FF_DIM)?,
            ln1: LayerNorm::load(st, &format!("{prefix}.norm1"), D_MODEL)?,
            ln2: LayerNorm::load(st, &format!("{prefix}.norm2"), D_MODEL)?,
        })
    }
}

/// The frozen CMC scoring network.
pub struct CmcNet {
    conv1: Conv1d,
    bn1: BatchNorm,
    conv2: Conv1d,
    bn2: BatchNorm,
    cls_token: Array1<f32>,
    pos_embed: Vec<f32>, // (N_TOKENS, D_MODEL) row-major
    layers: Vec<EncoderLayer>,
    final_ln: LayerNorm,
    dense: Linear,
    classifier: Linear,
}

/// Row-major f32 load of a 3-D tensor with a shape check.
fn take_3d(
    st: &safetensors::SafeTensors,
    name: &str,
    d0: usize,
    d1: usize,
    d2: usize,
) -> Result<Vec<f32>, ScoreError> {
    let view = st.tensor(name)?;
    let shape = view.shape().to_vec();
    if shape != vec![d0, d1, d2] {
        return Err(ScoreError::Shape {
            name: name.into(),
            got: shape,
            want: vec![d0, d1, d2],
        });
    }
    let data = view.data();
    let n = d0 * d1 * d2;
    if data.len() != n * std::mem::size_of::<f32>() {
        return Err(ScoreError::Shape {
            name: name.into(),
            got: vec![data.len()],
            want: vec![n * std::mem::size_of::<f32>()],
        });
    }
    // safetensors tensors are little-endian; the target is LE too, but parse
    // explicitly to stay defined on any host.
    Ok(data
        .chunks_exact(4)
        .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect())
}

fn gelu(v: f32) -> f32 {
    let x = f64::from(v);
    (0.5 * x * (1.0 + libm::erf(x * std::f64::consts::FRAC_1_SQRT_2))) as f32
}

/// Batched GELU over a whole matrix.
fn gelu_inplace(x: &mut Array2<f32>) {
    x.mapv_inplace(gelu);
}

impl CmcNet {
    /// Load from the safetensors weight file (names as exported by the
    /// training pipeline; every shape is validated).
    pub fn load(st: &safetensors::SafeTensors) -> Result<Self, ScoreError> {
        let cls = take_3d(st, "embed.cls_token", 1, 1, D_MODEL)?;
        Ok(Self {
            conv1: Conv1d::load(st, "conv1", CONV1_OUT, CMC_CHANNELS, CONV1_K)?,
            bn1: BatchNorm::load(st, "bn1", CONV1_OUT)?,
            conv2: Conv1d::load(st, "conv2", CONV2_OUT, CONV1_OUT, CONV2_K)?,
            bn2: BatchNorm::load(st, "bn2", CONV2_OUT)?,
            cls_token: Array1::from(cls),
            pos_embed: take_3d(st, "embed.position_embedding", 1, N_TOKENS, D_MODEL)?,
            layers: (0..N_LAYERS)
                .map(|i| EncoderLayer::load(st, &format!("transformer.layers.{i}")))
                .collect::<Result<_, _>>()?,
            final_ln: LayerNorm::load(st, "transformer.norm", D_MODEL)?,
            dense: Linear::load(st, "dense", D_MODEL, D_MODEL)?,
            classifier: Linear::load(st, "classifier", 2, D_MODEL)?,
        })
    }


    /// RE_PROB = softmax(logits)[1] for one CMC tensor.
    pub fn re_prob(&self, x: &crate::cmc::CmcFeatures) -> f64 {
        // conv stack over the (16,50) input laid out row-major
        let mut h = self.conv1.forward_same(x, CMC_COLUMNS);
        self.bn1.apply(&mut h);
        gelu_inplace(&mut h);
        let flat: Vec<f32> = h.iter().copied().collect();
        let mut h = self.conv2.forward_same(&flat, CMC_COLUMNS);
        self.bn2.apply(&mut h);
        gelu_inplace(&mut h);
        // transpose (128,50) -> tokens (51,128), cls token first, + position embedding
        let mut tokens = Array2::<f32>::zeros((N_TOKENS, D_MODEL));
        for t in 0..N_TOKENS - 1 {
            for c in 0..D_MODEL {
                tokens[[t + 1, c]] = h[[c, t]];
            }
        }
        for c in 0..D_MODEL {
            tokens[[0, c]] = self.cls_token[c];
        }
        for t in 0..N_TOKENS {
            let base = t * D_MODEL;
            for c in 0..D_MODEL {
                tokens[[t, c]] += self.pos_embed[base + c];
            }
        }
        // transformer encoder (pre-norm residual blocks) + final layer norm
        for layer in &self.layers {
            let mut normed = Array2::<f32>::zeros((N_TOKENS, D_MODEL));
            for t in 0..N_TOKENS {
                normed
                    .row_mut(t)
                    .assign(&layer.ln1.norm(&tokens.row(t).to_owned()));
            }
            let attn = attention(layer, &normed);
            tokens = &tokens + &attn;
            let mut normed = Array2::<f32>::zeros((N_TOKENS, D_MODEL));
            for t in 0..N_TOKENS {
                normed
                    .row_mut(t)
                    .assign(&layer.ln2.norm(&tokens.row(t).to_owned()));
            }
            let mut ff = Array2::<f32>::zeros((N_TOKENS, D_MODEL));
            for t in 0..N_TOKENS {
                let mut v = layer.ff1.forward(&normed.row(t).to_owned());
                v.mapv_inplace(gelu);
                ff.row_mut(t).assign(&layer.ff2.forward(&v));
            }
            tokens = &tokens + &ff;
        }
        let cls = self.final_ln.norm(&tokens.row(0).to_owned());
        let mut v = self.dense.forward(&cls);
        v.mapv_inplace(gelu);
        let logits = self.classifier.forward(&v);
        // softmax over 2 classes (f64 accumulation)
        let m = f64::from(logits[0].max(logits[1]));
        let e0 = (f64::from(logits[0]) - m).exp();
        let e1 = (f64::from(logits[1]) - m).exp();
        e1 / (e0 + e1)
    }
}

/// Multi-head self-attention over the token axis (no mask, eval mode).
fn attention(layer: &EncoderLayer, x: &Array2<f32>) -> Array2<f32> {
    // qkv: (N_TOKENS, 3*D)
    let mut qkv = Array2::<f32>::zeros((N_TOKENS, 3 * D_MODEL));
    for t in 0..N_TOKENS {
        let row = x.row(t).to_owned();
        let q = layer.qkv_w.dot(&row) + &layer.qkv_b;
        qkv.row_mut(t).assign(&q);
    }
    let mut out = Array2::<f32>::zeros((N_TOKENS, D_MODEL));
    for h in 0..N_HEADS {
        let off_q = h * HEAD_DIM;
        let off_k = D_MODEL + h * HEAD_DIM;
        let off_v = 2 * D_MODEL + h * HEAD_DIM;
        let mut scores = Array2::<f32>::zeros((N_TOKENS, N_TOKENS));
        for i in 0..N_TOKENS {
            for j in 0..N_TOKENS {
                let mut s = 0.0f32;
                for d in 0..HEAD_DIM {
                    s += qkv[[i, off_q + d]] * qkv[[j, off_k + d]];
                }
                scores[[i, j]] = s * (1.0 / (HEAD_DIM as f32).sqrt());
            }
        }
        let mut probs = Array2::<f32>::zeros((N_TOKENS, N_TOKENS));
        for i in 0..N_TOKENS {
            let row_max = f64::from(
                (0..N_TOKENS)
                    .map(|j| scores[[i, j]])
                    .fold(f32::NEG_INFINITY, f32::max),
            );
            let exps: Vec<f64> = (0..N_TOKENS)
                .map(|j| (f64::from(scores[[i, j]]) - row_max).exp())
                .collect();
            let sum: f64 = exps.iter().sum();
            for j in 0..N_TOKENS {
                probs[[i, j]] = (exps[j] / sum) as f32;
            }
        }
        for i in 0..N_TOKENS {
            for d in 0..HEAD_DIM {
                let mut s = 0.0f32;
                for j in 0..N_TOKENS {
                    s += probs[[i, j]] * qkv[[j, off_v + d]];
                }
                out[[i, h * HEAD_DIM + d]] = s;
            }
        }
    }
    // out projection
    let mut projected = Array2::<f32>::zeros((N_TOKENS, D_MODEL));
    for t in 0..N_TOKENS {
        let row = out.row(t).to_owned();
        projected.row_mut(t).assign(&layer.out.forward(&row));
    }
    projected
}

#[cfg(test)]
mod tests {
    use super::*;

    struct ZeroTensor {
        shape: Vec<usize>,
    }

    impl safetensors::View for &ZeroTensor {
        fn dtype(&self) -> safetensors::Dtype {
            safetensors::Dtype::F32
        }
        fn shape(&self) -> &[usize] {
            &self.shape
        }
        fn data(&self) -> std::borrow::Cow<'_, [u8]> {
            std::borrow::Cow::Owned(vec![0u8; self.data_len()])
        }
        fn data_len(&self) -> usize {
            self.shape.iter().product::<usize>() * 4
        }
    }

    /// A zero-weight network of the frozen architecture loads and yields
    /// softmax([0, 0]) = 0.5 for any input (all biases zero).
    #[test]
    fn zero_network_outputs_one_half() {
        let mut entries: Vec<(String, ZeroTensor)> = Vec::new();
        let zero = |name: String, shape: &[usize], entries: &mut Vec<(String, ZeroTensor)>| {
            entries.push((name, ZeroTensor { shape: shape.to_vec() }));
        };
        zero("conv1.weight".into(), &[CONV1_OUT, CMC_CHANNELS, CONV1_K], &mut entries);
        zero("conv1.bias".into(), &[CONV1_OUT], &mut entries);
        zero("bn1.weight".into(), &[CONV1_OUT], &mut entries);
        zero("bn1.bias".into(), &[CONV1_OUT], &mut entries);
        zero("bn1.running_mean".into(), &[CONV1_OUT], &mut entries);
        zero("bn1.running_var".into(), &[CONV1_OUT], &mut entries);
        zero("conv2.weight".into(), &[CONV2_OUT, CONV1_OUT, CONV2_K], &mut entries);
        zero("conv2.bias".into(), &[CONV2_OUT], &mut entries);
        zero("bn2.weight".into(), &[CONV2_OUT], &mut entries);
        zero("bn2.bias".into(), &[CONV2_OUT], &mut entries);
        zero("bn2.running_mean".into(), &[CONV2_OUT], &mut entries);
        zero("bn2.running_var".into(), &[CONV2_OUT], &mut entries);
        zero("embed.cls_token".into(), &[1, 1, D_MODEL], &mut entries);
        zero("embed.position_embedding".into(), &[1, N_TOKENS, D_MODEL], &mut entries);
        zero("dense.weight".into(), &[D_MODEL, D_MODEL], &mut entries);
        zero("dense.bias".into(), &[D_MODEL], &mut entries);
        zero("classifier.weight".into(), &[2, D_MODEL], &mut entries);
        zero("classifier.bias".into(), &[2], &mut entries);
        for layer in 0..N_LAYERS {
            let p = format!("transformer.layers.{layer}");
            zero(format!("{p}.self_attn.in_proj_weight"), &[3 * D_MODEL, D_MODEL], &mut entries);
            zero(format!("{p}.self_attn.in_proj_bias"), &[3 * D_MODEL], &mut entries);
            zero(format!("{p}.self_attn.out_proj.weight"), &[D_MODEL, D_MODEL], &mut entries);
            zero(format!("{p}.self_attn.out_proj.bias"), &[D_MODEL], &mut entries);
            zero(format!("{p}.linear1.weight"), &[FF_DIM, D_MODEL], &mut entries);
            zero(format!("{p}.linear1.bias"), &[FF_DIM], &mut entries);
            zero(format!("{p}.linear2.weight"), &[D_MODEL, FF_DIM], &mut entries);
            zero(format!("{p}.linear2.bias"), &[D_MODEL], &mut entries);
            zero(format!("{p}.norm1.weight"), &[D_MODEL], &mut entries);
            zero(format!("{p}.norm1.bias"), &[D_MODEL], &mut entries);
            zero(format!("{p}.norm2.weight"), &[D_MODEL], &mut entries);
            zero(format!("{p}.norm2.bias"), &[D_MODEL], &mut entries);
        }
        zero("transformer.norm.weight".into(), &[D_MODEL], &mut entries);
        zero("transformer.norm.bias".into(), &[D_MODEL], &mut entries);
        let doc = safetensors::serialize(entries.iter().map(|(n, v)| (n.as_str(), v)), &None)
            .expect("serialize");
        let st = safetensors::SafeTensors::deserialize(&doc).expect("parse");
        let net = CmcNet::load(&st).expect("load");
        let x = [0.25f32; CMC_CHANNELS * CMC_COLUMNS];
        let p = net.re_prob(&x);
        assert!((p - 0.5).abs() < 1e-12, "zero network must output 0.5, got {p}");
        let x2 = [0.0f32; CMC_CHANNELS * CMC_COLUMNS];
        assert!((net.re_prob(&x2) - 0.5).abs() < 1e-12);
    }

    /// A missing tensor is a hard load error (frozen shape contract).
    #[test]
    fn load_rejects_missing_tensors() {
        let doc = safetensors::serialize(
            [(
                "conv1.weight",
                &ZeroTensor { shape: vec![CONV1_OUT, CMC_CHANNELS, CONV1_K] },
            )],
            &None,
        )
        .expect("serialize");
        let st = safetensors::SafeTensors::deserialize(&doc).expect("parse");
        assert!(CmcNet::load(&st).is_err());
    }
}
