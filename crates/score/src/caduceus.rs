//! Caduceus-Ph 4L (d118) native forward: embedding + 4x BiMamba pre-norm block + final RMSNorm
//! + mean-pool (over L, fp32) -> fp16. Weight keys are prefixed with `caduceus.backbone.*` (all F32).
//!
//! Residual stream (mamba_ssm Block, fused_add_norm=False path).
//!
//! Per layer: residual = hidden + residual (first layer: residual = hidden); hidden = RMSNorm(residual);
//! hidden = mamba_fwd(x) + flip(mamba_rev(flip(x))) (flip along L only; rev's in/out_proj tied to fwd).
//! Finally: residual = hidden + residual; out = RMSNorm_f(residual).

use crate::bundle::{take_matrix, take_vector, ScoreError};
use crate::mamba::{MambaWeights, D_INNER, D_STATE, DT_RANK};
use half::f16;
use ndarray::{s, Array1, Array2, Axis};
use std::path::Path;

pub const D_MODEL: usize = 118;
pub const N_LAYER: usize = 4;
pub const VOCAB: usize = 16;
pub const NORM_EPS: f32 = 1e-5;

struct LayerWeights {
    norm: Array1<f32>, // (D_MODEL,)
    fwd: MambaWeights,
    rev: MambaWeights,
}

pub struct CaduceusEncoder {
    emb: Array2<f32>, // (VOCAB, D_MODEL)
    layers: Vec<LayerWeights>,
    norm_f: Array1<f32>, // (D_MODEL,)
}

/// Pure fp32 RMSNorm (matches `_TorchRMSNorm` in encoder.py): x*rsqrt(mean(x^2)+eps)*weight.
fn rms_norm(x: &Array2<f32>, w: &Array1<f32>) -> Array2<f32> {
    let mut out = x.clone();
    out.axis_iter_mut(Axis(0)).for_each(|mut row| {
        let ms = row.iter().map(|v| v * v).sum::<f32>() / row.len() as f32;
        let inv = (ms + NORM_EPS).sqrt().recip();
        for (v, &wv) in row.iter_mut().zip(w.iter()) {
            *v = *v * inv * wv;
        }
    });
    out
}

/// conv1d.weight [D_INNER, 1, D_CONV] -> (D_INNER, D_CONV).
fn take_conv(st: &safetensors::SafeTensors, name: &str) -> Result<Array2<f32>, ScoreError> {
    let t = st.tensor(name).map_err(|_| ScoreError::Missing(name.into()))?;
    let got: Vec<usize> = t.shape().to_vec();
    let want = vec![D_INNER, 1, crate::mamba::D_CONV];
    if got != want {
        return Err(ScoreError::Shape {
            name: name.into(),
            got,
            want,
        });
    }
    let data: Vec<f32> = t
        .data()
        .chunks_exact(4)
        .map(|b| f32::from_le_bytes(b.try_into().unwrap()))
        .collect();
    Ok(Array2::from_shape_vec((D_INNER, crate::mamba::D_CONV), data).unwrap())
}

fn load_mamba(
    st: &safetensors::SafeTensors,
    prefix: &str,
    tie: Option<(&Array2<f32>, &Array2<f32>)>,
) -> Result<MambaWeights, ScoreError> {
    let (in_proj, out_proj) = match tie {
        Some((i, o)) => (i.clone(), o.clone()),
        None => (
            take_matrix(st, &format!("{prefix}.in_proj.weight"), 2 * D_INNER, D_MODEL)?,
            take_matrix(st, &format!("{prefix}.out_proj.weight"), D_MODEL, D_INNER)?,
        ),
    };
    let conv_w = take_conv(st, &format!("{prefix}.conv1d.weight"))?;
    let conv_b = take_vector(st, &format!("{prefix}.conv1d.bias"), D_INNER)?;
    let x_proj = take_matrix(
        st,
        &format!("{prefix}.x_proj.weight"),
        DT_RANK + 2 * D_STATE,
        D_INNER,
    )?;
    let dt_w = take_matrix(st, &format!("{prefix}.dt_proj.weight"), D_INNER, DT_RANK)?;
    let dt_b = take_vector(st, &format!("{prefix}.dt_proj.bias"), D_INNER)?;
    let a_log = take_matrix(st, &format!("{prefix}.A_log"), D_INNER, D_STATE)?;
    let a = a_log.mapv(|v| -v.exp());
    let d = take_vector(st, &format!("{prefix}.D"), D_INNER)?;
    Ok(MambaWeights {
        in_proj,
        conv_w,
        conv_b,
        x_proj,
        dt_w,
        dt_b,
        a,
        d,
        out_proj,
    })
}

impl CaduceusEncoder {
    /// Load a checkpoint directory (containing model.safetensors). rev projections: load independently if standalone weights exist (v07 self-trained); otherwise (caduceus-ph weight_tie) tie to fwd. Compatible with both backbones.
    pub fn load(ckpt_dir: &Path) -> Result<Self, ScoreError> {
        let bytes = std::fs::read(ckpt_dir.join("model.safetensors"))?;
        let st = safetensors::SafeTensors::deserialize(&bytes)?;
        let emb = take_matrix(
            &st,
            "caduceus.backbone.embeddings.word_embeddings.weight",
            VOCAB,
            D_MODEL,
        )?;
        let norm_f = take_vector(&st, "caduceus.backbone.norm_f.weight", D_MODEL)?;
        let mut layers = Vec::with_capacity(N_LAYER);
        for i in 0..N_LAYER {
            let p = format!("caduceus.backbone.layers.{i}");
            let norm = take_vector(&st, &format!("{p}.norm.weight"), D_MODEL)?;
            let fwd = load_mamba(&st, &format!("{p}.mixer.mamba_fwd"), None)?;
            // rev: try standalone load first (v07 self-trained standalone rev); on failure (caduceus-ph weight_tie has no standalone keys) tie to fwd.
            let rev_prefix = format!("{p}.mixer.mamba_rev");
            let rev = match load_mamba(&st, &rev_prefix, None) {
                Ok(r) => r,
                Err(_) => load_mamba(&st, &rev_prefix, Some((&fwd.in_proj, &fwd.out_proj)))?,
            };
            layers.push(LayerWeights { norm, fwd, rev });
        }
        Ok(Self {
            emb,
            layers,
            norm_f,
        })
    }

    /// tokens (L=1002) -> mean-pooled fp16 embedding (118,). Read-only: &self, shareable via Arc.
    pub fn embed(&self, tokens: &[i64]) -> Result<[f16; D_MODEL], ScoreError> {
        let l = tokens.len();
        let mut hidden = Array2::<f32>::zeros((l, D_MODEL));
        for (i, &t) in tokens.iter().enumerate() {
            if !(0..VOCAB as i64).contains(&t) {
                return Err(ScoreError::Token(t));
            }
            hidden.row_mut(i).assign(&self.emb.row(t as usize));
        }

        let mut residual: Option<Array2<f32>> = None;
        for layer in &self.layers {
            let res = match residual {
                None => hidden.clone(),
                Some(r) => &hidden + &r,
            };
            let normed = rms_norm(&res, &layer.norm);
            let fwd_out = layer.fwd.forward(&normed);
            let flipped = normed.slice(s![..;-1, ..]).to_owned();
            let rev_out = layer.rev.forward(&flipped);
            let rev_out = rev_out.slice(s![..;-1, ..]).to_owned();
            hidden = &fwd_out + &rev_out;
            residual = Some(res);
        }
        let res = &hidden + &residual.expect("N_LAYER >= 1");
        let out = rms_norm(&res, &self.norm_f);

        // mean over L (fp32) -> fp16
        let mean = out.mean_axis(Axis(0)).unwrap();
        let mut emb = [f16::ZERO; D_MODEL];
        for (i, v) in mean.iter().enumerate() {
            emb[i] = f16::from_f32(*v);
        }
        Ok(emb)
    }
}

// =========================================================================
// Batched encoding (embed_batch) -- one batch of sites shares GEMM/scan; residual structure mirrors embed step by step
// =========================================================================

use crate::mamba::BatchBufs;

/// Working buffers for embed_batch (one per rayon worker).
#[derive(Default)]
pub struct EmbedBatchBufs {
    mamba: BatchBufs,
    x: Vec<f32>,       // (B*L, D_MODEL) embedding input
    hidden: Vec<f32>,  // (B*L, D_MODEL) fwd+rev mixed result
    residual: Vec<f32>,
    normed: Vec<f32>,
    flipped: Vec<f32>,
    rev: Vec<f32>,
    out: Vec<f32>,
}

impl EmbedBatchBufs {
    fn resize(&mut self, bl: usize) {
        let n = bl * D_MODEL;
        self.x.resize(n, 0.0);
        self.hidden.resize(n, 0.0);
        self.residual.resize(n, 0.0);
        self.normed.resize(n, 0.0);
        self.flipped.resize(n, 0.0);
        self.rev.resize(n, 0.0);
        self.out.resize(n, 0.0);
    }
}

/// Batched RMSNorm: per-row (site x timestep) ms summed in order, same formula as the single-sequence rms_norm.
fn rms_norm_batch(x: &[f32], w: &Array1<f32>, out: &mut [f32]) {
    for (row, orow) in x.chunks_exact(D_MODEL).zip(out.chunks_exact_mut(D_MODEL)) {
        let ms = row.iter().map(|v| v * v).sum::<f32>() / D_MODEL as f32;
        let inv = (ms + NORM_EPS).sqrt().recip();
        for (o, (&v, &wv)) in orow.iter_mut().zip(row.iter().zip(w.iter())) {
            *o = v * inv * wv;
        }
    }
}

/// Flip along L (per site): out[(b*L+i)*D+k] = x[(b*L+(L-1-i))*D+k].
fn flip_batch(x: &[f32], out: &mut [f32], b: usize, l: usize) {
    for bi in 0..b {
        let so = bi * l * D_MODEL;
        let src = &x[so..so + l * D_MODEL];
        let dst = &mut out[so..so + l * D_MODEL];
        for i in 0..l {
            dst[i * D_MODEL..(i + 1) * D_MODEL]
                .copy_from_slice(&src[(l - 1 - i) * D_MODEL..(l - i) * D_MODEL]);
        }
    }
}

impl CaduceusEncoder {
    /// Batched encoding: tokens_flat (B*L) -> per-site fp16 embedding (B, 118).
    /// Isomorphic to embed() step by step (residual stream / norm / flip / mean-pool order unchanged);
    /// GEMM is executed by faer (accumulation order differs from ndarray, fp32 difference ~1e-6,
    /// fp16 quantization happens after the difference).
    pub fn embed_batch(
        &self,
        tokens_flat: &[i64],
        b: usize,
        l: usize,
        bufs: &mut EmbedBatchBufs,
    ) -> Result<Vec<[f16; D_MODEL]>, ScoreError> {
        let bl = b * l;
        if tokens_flat.len() != bl {
            return Err(ScoreError::Token(-1));
        }
        bufs.resize(bl);
        // embedding lookup -> x (site-major rows)
        for (j, &t) in tokens_flat.iter().enumerate() {
            if !(0..VOCAB as i64).contains(&t) {
                return Err(ScoreError::Token(t));
            }
            let src = self.emb.row(t as usize);
            bufs.x[j * D_MODEL..(j + 1) * D_MODEL].copy_from_slice(src.as_slice().unwrap());
        }

        let mut have_residual = false;
        for layer in &self.layers {
            // res = hidden + residual (first layer: residual = x)
            if have_residual {
                for j in 0..bl * D_MODEL {
                    bufs.residual[j] += bufs.hidden[j];
                }
            } else {
                bufs.residual.copy_from_slice(&bufs.x);
                have_residual = true;
            }
            rms_norm_batch(&bufs.residual, &layer.norm, &mut bufs.normed);
            // fwd
            layer.fwd.forward_batch(&bufs.normed, b, l, &mut bufs.mamba);
            bufs
                .hidden
                .copy_from_slice(crate::mamba::MambaWeights::batch_output(&bufs.mamba));
            // rev: flip -> forward -> flip back
            flip_batch(&bufs.normed, &mut bufs.flipped, b, l);
            layer.rev.forward_batch(&bufs.flipped, b, l, &mut bufs.mamba);
            flip_batch(
                crate::mamba::MambaWeights::batch_output(&bufs.mamba),
                &mut bufs.rev,
                b,
                l,
            );
            for j in 0..bl * D_MODEL {
                bufs.hidden[j] += bufs.rev[j];
            }
        }
        // out = RMSNorm_f(hidden + residual)
        for j in 0..bl * D_MODEL {
            bufs.out[j] = bufs.hidden[j] + bufs.residual[j];
        }
        rms_norm_batch(&bufs.out, &self.norm_f, &mut bufs.normed);

        // mean over L (order 0..L, same as mean_axis) -> fp16
        let mut embs = Vec::with_capacity(b);
        for bi in 0..b {
            let base = bi * l * D_MODEL;
            let mut mean = [0f32; D_MODEL];
            for i in 0..l {
                let row = &bufs.normed[base + i * D_MODEL..base + (i + 1) * D_MODEL];
                for (m, &v) in mean.iter_mut().zip(row) {
                    *m += v;
                }
            }
            let mut e = [f16::ZERO; D_MODEL];
            for (k, m) in mean.iter().enumerate() {
                e[k] = f16::from_f32(*m / l as f32);
            }
            embs.push(e);
        }
        Ok(embs)
    }

    /// Per-position normed-hidden variant of embed_batch (for cluster-merge probes/experiments):
    /// returns the (b*l, D_MODEL) post-RMSNorm hidden (fp32, the direct input of mean-pool).
    /// Shares all forward code with embed_batch -- duplicated rather than refactored so the existing
    /// path stays untouched; pooling is performed by the caller over arbitrary sub-ranges
    /// (the production mean is the special case of sub-range [0,l)).
    pub fn embed_batch_hidden(
        &self,
        tokens_flat: &[i64],
        b: usize,
        l: usize,
        bufs: &mut EmbedBatchBufs,
    ) -> Result<Vec<f32>, ScoreError> {
        let bl = b * l;
        if tokens_flat.len() != bl {
            return Err(ScoreError::Token(-1));
        }
        bufs.resize(bl);
        for (j, &t) in tokens_flat.iter().enumerate() {
            if !(0..VOCAB as i64).contains(&t) {
                return Err(ScoreError::Token(t));
            }
            let src = self.emb.row(t as usize);
            bufs.x[j * D_MODEL..(j + 1) * D_MODEL].copy_from_slice(src.as_slice().unwrap());
        }
        let mut have_residual = false;
        for layer in &self.layers {
            if have_residual {
                for j in 0..bl * D_MODEL {
                    bufs.residual[j] += bufs.hidden[j];
                }
            } else {
                bufs.residual.copy_from_slice(&bufs.x);
                have_residual = true;
            }
            rms_norm_batch(&bufs.residual, &layer.norm, &mut bufs.normed);
            layer.fwd.forward_batch(&bufs.normed, b, l, &mut bufs.mamba);
            bufs
                .hidden
                .copy_from_slice(crate::mamba::MambaWeights::batch_output(&bufs.mamba));
            flip_batch(&bufs.normed, &mut bufs.flipped, b, l);
            layer.rev.forward_batch(&bufs.flipped, b, l, &mut bufs.mamba);
            flip_batch(
                crate::mamba::MambaWeights::batch_output(&bufs.mamba),
                &mut bufs.rev,
                b,
                l,
            );
            for j in 0..bl * D_MODEL {
                bufs.hidden[j] += bufs.rev[j];
            }
        }
        for j in 0..bl * D_MODEL {
            bufs.out[j] = bufs.hidden[j] + bufs.residual[j];
        }
        rms_norm_batch(&bufs.out, &self.norm_f, &mut bufs.normed);
        Ok(bufs.normed.clone())
    }
}
