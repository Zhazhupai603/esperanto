//! GPU (CUDA) forward for the Caduceus-Ph encoder via candle — the numerical twin of the CPU
//! `embed_batch` path: same weight mapping (reuses `CaduceusEncoder::load`, including the
//! rev weight-tie fallback and `A = -exp(A_log)`), same op semantics (causal conv1d k=4 with
//! left-pad 3, silu, softplus with torch's threshold=20, the scan recurrence, flip->forward->flip
//! for the rev direction, pre-norm residual stream, final RMSNorm, mean-pool over L, fp16 out).
//!
//! Scan strategy (decided by measurement on the target RTX 4090, B=64, L=502): the LITERAL
//! sequential 502-step recurrence, kept as `h_t = exp(delta_t·A) ⊙ h_{t-1} + (delta_t·u_t)·B_t`
//! with ~5 small (B, D_INNER, D_STATE) kernels per step. A blocked affine associative scan
//! (Hillis–Steele doubling, two levels) was implemented and measured ~2.3x SLOWER end-to-end
//! (527 vs 225 ms per batch): it streams ~40 GB per direction-layer through DRAM, while the
//! literal loop's small kernels stay launch-bound and pipeline behind the CUDA driver. Numerics:
//! `a = exp(delta·A)` is always in (0, 1) (delta > 0 from softplus, A < 0); the doubling variant
//! was overflow-free for the same reason (products stay in [0, 1]) — it lost on speed, not
//! correctness. Parity vs the CPU recurrence is measured by `examples/gpu_parity.rs`.

use crate::bundle::ScoreError;
use crate::caduceus::{CaduceusEncoder, D_MODEL, NORM_EPS, VOCAB};
use crate::mamba::{MambaWeights, DT_RANK, D_CONV, D_INNER, D_STATE};
use candle_core::{Device, Tensor};
use half::f16;
use ndarray::{Array1, Array2};
use std::path::Path;

type Cdl = Result<Tensor, candle_core::Error>;

fn matrix(dev: &Device, a: &Array2<f32>) -> Cdl {
    let (r, c) = a.dim();
    Tensor::from_vec(a.as_slice().unwrap().to_vec(), (r, c), dev)
}

fn row_vec(dev: &Device, a: &Array1<f32>) -> Cdl {
    let v = a.to_vec();
    let n = v.len();
    Tensor::from_vec(v, n, dev)
}

/// Pure fp32 RMSNorm, GPU twin of `caduceus::rms_norm_batch`: x * rsqrt(mean(x^2)+eps) * w.
fn rms_norm(x: &Tensor, w: &Tensor) -> Cdl {
    let ms = x.sqr()?.mean_keepdim(2)?;
    let inv = (ms + NORM_EPS as f64)?.sqrt()?.recip()?;
    x.broadcast_mul(&inv)?.broadcast_mul(w)
}

struct GpuMamba {
    in_proj: Tensor,     // (2*D_INNER, D_MODEL)
    conv_w: Vec<Tensor>, // D_CONV x (1, 1, D_INNER) split taps
    conv_b: Tensor,      // (1, 1, D_INNER)
    x_proj: Tensor,      // (DT_RANK + 2*D_STATE, D_INNER)
    dt_w: Tensor,        // (D_INNER, DT_RANK)
    dt_b: Tensor,        // (1, 1, D_INNER)
    a: Tensor,           // (1, 1, D_INNER, D_STATE) = -exp(A_log), precomputed at CPU load
    d: Tensor,           // (1, 1, D_INNER)
    out_proj: Tensor,    // (D_MODEL, D_INNER)
}

impl GpuMamba {
    fn load(w: &MambaWeights, dev: &Device) -> Result<Self, ScoreError> {
        let conv_tap = |k: usize| -> Cdl {
            let v: Vec<f32> = (0..D_INNER).map(|d| w.conv_w[(d, k)]).collect();
            Tensor::from_vec(v, (1, 1, D_INNER), dev)
        };
        let a_flat = w.a.as_slice().unwrap().to_vec();
        Ok(Self {
            in_proj: matrix(dev, &w.in_proj)?,
            conv_w: (0..D_CONV).map(conv_tap).collect::<Result<_, _>>()?,
            conv_b: row_vec(dev, &w.conv_b)?.reshape((1, 1, D_INNER))?,
            x_proj: matrix(dev, &w.x_proj)?,
            dt_w: matrix(dev, &w.dt_w)?,
            dt_b: row_vec(dev, &w.dt_b)?.reshape((1, 1, D_INNER))?,
            a: Tensor::from_vec(a_flat, (1, 1, D_INNER, D_STATE), dev)?,
            d: row_vec(dev, &w.d)?.reshape((1, 1, D_INNER))?,
            out_proj: matrix(dev, &w.out_proj)?,
        })
    }

    /// Batched single-direction forward: x (b, l, D_MODEL) contiguous -> (b, l, D_MODEL).
    /// Mirrors `MambaWeights::forward_batch` step by step (see mamba.rs header for the contract).
    fn forward(&self, x3: &Tensor, b: usize, l: usize) -> Cdl {
        let bl = b * l;

        // xz = x @ in_proj^T -> (bl, 472); split channels / gate (strided views)
        let xf = x3.reshape((bl, D_MODEL))?;
        let xz = xf.matmul(&self.in_proj.t()?)?;
        let xz3 = xz.reshape((b, l, 2 * D_INNER))?;
        let xv = xz3.narrow(2, 0, D_INNER)?;
        let zv = xz3.narrow(2, D_INNER, D_INNER)?;
        let z_gate = zv.silu()?; // silu(z), contiguous (b, l, D_INNER)

        // causal conv1d (k=4, left-pad 3, take first L) + bias, then silu
        let xp = xv.pad_with_zeros(1, D_CONV - 1, 0)?; // (b, l+3, D_INNER)
        let mut u = xp.narrow(1, 0, l)?.broadcast_mul(&self.conv_w[0])?;
        for (k, w) in self.conv_w.iter().enumerate().skip(1) {
            u = u.broadcast_add(&xp.narrow(1, k, l)?.broadcast_mul(w)?)?;
        }
        let u = u.broadcast_add(&self.conv_b)?.silu()?; // (b, l, D_INNER) contiguous

        // x_dbl = u @ x_proj^T -> (bl, 40); split dt/B/C
        let uf = u.reshape((bl, D_INNER))?;
        let xd = uf.matmul(&self.x_proj.t()?)?;
        let xd3 = xd.reshape((b, l, DT_RANK + 2 * D_STATE))?;
        let dt8 = xd3
            .narrow(2, 0, DT_RANK)?
            .contiguous()?
            .reshape((bl, DT_RANK))?;
        let bm = xd3.narrow(2, DT_RANK, D_STATE)?.contiguous()?; // (b, l, N)
        let cm = xd3.narrow(2, DT_RANK + D_STATE, D_STATE)?.contiguous()?;

        // delta = softplus(dt8 @ dt_w^T + b_dt), torch threshold=20 branch
        let delta = dt8.matmul(&self.dt_w.t()?)?.reshape((b, l, D_INNER))?;
        let delta = delta.broadcast_add(&self.dt_b)?;
        let e = delta.exp()?;
        let sp = (&e + 1.0)?.log()?; // ln(1 + exp(x))
        let delta = delta.gt(20f32)?.where_cond(&delta, &sp)?;

        // selective scan -> y (b, l, D_INNER)
        let y = self.scan(&delta, &u, &bm, &cm, b, l)?;

        // out = (y + D·u) ⊙ silu(z); out_proj -> (b, l, D_MODEL)
        let gated = y
            .broadcast_add(&u.broadcast_mul(&self.d)?)?
            .broadcast_mul(&z_gate)?;
        gated
            .reshape((bl, D_INNER))?
            .matmul(&self.out_proj.t()?)?
            .reshape((b, l, D_MODEL))
    }

    /// Selective scan: delta/u (b, l, D_INNER), bm/cm (b, l, D_STATE) -> y (b, l, D_INNER)
    /// with y_t = C_t · h_t, h_0 = 0 — the literal sequential recurrence, step by step.
    ///
    /// Structure decided by measurement (RTX 4090, B=64, L=502; see module docs): a blocked
    /// associative (Hillis–Steele) scan was ~2.4x SLOWER than this loop (527 vs 217 ms per
    /// batch end-to-end) — its ~40 GB/direction-layer of DRAM traffic dominates, while this
    /// loop's 5 small (b, d, n) kernels/step stay launch-bound and pipeline behind the CUDA
    /// driver. a_t = exp(delta_t·A) and du = delta·u are precomputed in two full-tensor passes
    /// so the loop body is 5 kernels (h update 3, y 2).
    fn scan(
        &self,
        delta: &Tensor,
        u: &Tensor,
        bm: &Tensor,
        cm: &Tensor,
        b: usize,
        l: usize,
    ) -> Cdl {
        let (d, n) = (D_INNER, D_STATE);
        let a4 = delta.unsqueeze(3)?.broadcast_mul(&self.a)?.exp()?; // (b, l, d, n)
        let du = delta.mul(u)?; // (b, l, d)
        let mut h = Tensor::zeros((b, d, n), delta.dtype(), delta.device())?;
        let mut ys: Vec<Tensor> = Vec::with_capacity(l);
        for t in 0..l {
            let at = a4.get_on_dim(1, t)?; // (b, d, n) views
            let du_t = du.get_on_dim(1, t)?;
            let bt = bm.get_on_dim(1, t)?;
            let ct = cm.get_on_dim(1, t)?;
            // h = a_t ⊙ h + du_t ⊗ B_t
            h = at
                .mul(&h)?
                .broadcast_add(&du_t.unsqueeze(2)?.broadcast_mul(&bt.unsqueeze(1)?)?)?;
            ys.push(h.broadcast_mul(&ct.unsqueeze(1)?)?.sum(2)?);
        }
        Tensor::stack(&ys, 1)
    }
}

struct GpuLayer {
    norm: Tensor, // (1, 1, D_MODEL)
    fwd: GpuMamba,
    rev: GpuMamba,
}

/// GPU twin of `CaduceusEncoder`: same weights uploaded to a CUDA device, batched forward only.
pub struct GpuCaduceusEncoder {
    dev: Device,
    emb: Tensor, // (VOCAB, D_MODEL)
    layers: Vec<GpuLayer>,
    norm_f: Tensor, // (1, 1, D_MODEL)
}

impl GpuCaduceusEncoder {
    /// Load a CPU checkpoint directory and upload it to `dev`.
    pub fn load(ckpt_dir: &Path, dev: &Device) -> Result<Self, ScoreError> {
        Self::from_cpu(&CaduceusEncoder::load(ckpt_dir)?, dev)
    }

    /// Upload an already-loaded CPU encoder (weight-tie resolution and A preprocessing included).
    pub fn from_cpu(cpu: &CaduceusEncoder, dev: &Device) -> Result<Self, ScoreError> {
        let norm3 = |v: &Array1<f32>| -> Cdl { row_vec(dev, v)?.reshape((1, 1, D_MODEL)) };
        let mut layers = Vec::with_capacity(cpu.layers.len());
        for ly in &cpu.layers {
            layers.push(GpuLayer {
                norm: norm3(&ly.norm)?,
                fwd: GpuMamba::load(&ly.fwd, dev)?,
                rev: GpuMamba::load(&ly.rev, dev)?,
            });
        }
        Ok(Self {
            dev: dev.clone(),
            emb: matrix(dev, &cpu.emb)?,
            layers,
            norm_f: norm3(&cpu.norm_f)?,
        })
    }

    /// Batched encoding, the GPU twin of `embed_batch`: tokens_flat (b*l, site-major) ->
    /// per-site fp16 embedding (b, 118). Same residual/norm/flip/mean-pool order as the CPU path.
    pub fn embed_batch(
        &self,
        tokens_flat: &[i64],
        b: usize,
        l: usize,
    ) -> Result<Vec<[f16; D_MODEL]>, ScoreError> {
        let bl = b * l;
        if tokens_flat.len() != bl {
            return Err(ScoreError::Token(-1));
        }
        let mut ids = Vec::with_capacity(bl);
        for &t in tokens_flat {
            if !(0..VOCAB as i64).contains(&t) {
                return Err(ScoreError::Token(t));
            }
            ids.push(t as u32);
        }
        let ids = Tensor::from_vec(ids, bl, &self.dev)?;

        let mut hidden = self.emb.embedding(&ids)?.reshape((b, l, D_MODEL))?;
        let mut residual: Option<Tensor> = None;
        for layer in &self.layers {
            // residual = hidden + residual (first layer: residual = x)
            let res = match &residual {
                None => hidden.clone(),
                Some(r) => r.add(&hidden)?,
            };
            let normed = rms_norm(&res, &layer.norm)?;
            // fwd
            let mut h = layer.fwd.forward(&normed, b, l)?;
            // rev: flip -> forward -> flip back
            let rev = layer.rev.forward(&normed.flip(&[1])?, b, l)?;
            h = h.add(&rev.flip(&[1])?)?;
            hidden = h;
            residual = Some(res);
        }

        // out = RMSNorm_f(hidden + residual); mean over L (fp32) -> fp16
        let res = hidden.add(&residual.expect("N_LAYER >= 1"))?;
        let out = rms_norm(&res, &self.norm_f)?;
        let rows = out.mean(1)?.to_vec2::<f32>()?; // (b, D_MODEL), synced download
        Ok(rows
            .into_iter()
            .map(|row| {
                let mut e = [f16::ZERO; D_MODEL];
                for (k, v) in row.iter().enumerate() {
                    e[k] = f16::from_f32(*v);
                }
                e
            })
            .collect())
    }
}
