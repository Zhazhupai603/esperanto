//! Single-direction Mamba block forward (mirrors mamba_simple with use_fast_path=False
//! and selective_scan_ref, step by step).
//!
//! xz = in_proj(x); x, z = split(xz); u = silu(causal_conv1d(x));
//! x_dbl = x_proj(u) -> dt, B, C; delta = softplus(W_dt*dt + b_dt) (bias added exactly once, inside the scan);
//! scan: h = exp(delta*A)*h + delta*B*u; y = C*h; out = (y + D*u) . silu(z); out = out_proj(out).
//! All in fp32. Channel-major (D, L) layout; the scan is parallelized across channels with rayon.

use ndarray::{s, Array1, Array2};
use rayon::prelude::*;

pub const D_INNER: usize = 236;
pub const D_STATE: usize = 16;
pub const D_CONV: usize = 4;
pub const DT_RANK: usize = 8;

pub struct MambaWeights {
    pub in_proj: Array2<f32>,  // (2*D_INNER, D_MODEL) = (472, 118)
    pub conv_w: Array2<f32>,   // (D_INNER, D_CONV)
    pub conv_b: Array1<f32>,   // (D_INNER,)
    pub x_proj: Array2<f32>,   // (DT_RANK + 2*D_STATE, D_INNER) = (40, 236)
    pub dt_w: Array2<f32>,     // (D_INNER, DT_RANK)
    pub dt_b: Array1<f32>,     // (D_INNER,)
    pub a: Array2<f32>,        // (D_INNER, D_STATE) = -exp(A_log), precomputed at load
    pub d: Array1<f32>,        // (D_INNER,)
    pub out_proj: Array2<f32>, // (D_MODEL, D_INNER)
}

#[inline]
fn silu(x: f32) -> f32 {
    x / (1.0 + (-x).exp())
}

#[inline]
fn softplus(x: f32) -> f32 {
    // torch F.softplus(beta=1, threshold=20)
    if x > 20.0 {
        x
    } else {
        x.exp().ln_1p()
    }
}

/// x (L,K) @ w(N,K)^T -> (L,N); rayon over row blocks (for large L).
fn dot_t(x: &Array2<f32>, w: &Array2<f32>) -> Array2<f32> {
    let (l, k) = x.dim();
    let n = w.nrows();
    debug_assert_eq!(w.ncols(), k);
    let mut out = vec![0f32; l * n];
    const CHUNK: usize = 32;
    out.par_chunks_mut(CHUNK * n)
        .enumerate()
        .for_each(|(ci, chunk)| {
            let r0 = ci * CHUNK;
            let r1 = (r0 + CHUNK).min(l);
            let prod = x.slice(s![r0..r1, ..]).dot(&w.t());
            chunk.copy_from_slice(prod.as_slice().unwrap());
        });
    Array2::from_shape_vec((l, n), out).unwrap()
}

/// out (N,L) = w (N,K) @ x (K,L); rayon over rows (for small N).
fn dot_l(w: &Array2<f32>, x: &Array2<f32>) -> Array2<f32> {
    let (n, k) = w.dim();
    let l = x.ncols();
    debug_assert_eq!(x.nrows(), k);
    let xs = x.as_slice().unwrap();
    let rows: Vec<Vec<f32>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let wr = w.row(i);
            let mut r = vec![0f32; l];
            for (j, &wv) in wr.iter().enumerate() {
                let xr = &xs[j * l..(j + 1) * l];
                for c in 0..l {
                    r[c] += wv * xr[c];
                }
            }
            r
        })
        .collect();
    Array2::from_shape_fn((n, l), |(i, j)| rows[i][j])
}

impl MambaWeights {
    /// x: (L, D_MODEL) -> out: (L, D_MODEL). fp32 throughout.
    pub fn forward(&self, x: &Array2<f32>) -> Array2<f32> {
        let l = x.nrows();

        // xz = in_proj(x) -> (L, 472); split into channel-major xc / zc (236, L)
        let xz = dot_t(x, &self.in_proj);
        let xc = xz.slice(s![.., ..D_INNER]).t().to_owned();
        let zc = xz.slice(s![.., D_INNER..]).t().to_owned();
        let xc_s = xc.as_slice().unwrap();

        // causal conv1d (groups=236, k=4, left-pad 3, take first L) + bias -> silu
        let mut u = Array2::<f32>::zeros((D_INNER, l));
        u.as_slice_mut()
            .unwrap()
            .par_chunks_mut(l)
            .enumerate()
            .for_each(|(d, ud)| {
                let w = self.conv_w.row(d);
                let bias = self.conv_b[d];
                let xrow = &xc_s[d * l..(d + 1) * l];
                for (i, v) in ud.iter_mut().enumerate() {
                    let mut acc = bias;
                    for k in 0..D_CONV {
                        if i + k >= D_CONV - 1 {
                            acc += w[k] * xrow[i + k - (D_CONV - 1)];
                        }
                    }
                    *v = silu(acc);
                }
            });

        // x_dbl^T = x_proj @ u -> (40, L); split into dt8 (8,L) / B (16,L) / C (16,L)
        let x_dbl = dot_l(&self.x_proj, &u);
        let dt8 = x_dbl.slice(s![..DT_RANK, ..]).to_owned();
        let bmat = x_dbl
            .slice(s![DT_RANK..DT_RANK + D_STATE, ..])
            .to_owned();
        let cmat = x_dbl.slice(s![DT_RANK + D_STATE.., ..]).to_owned();
        let b_s = bmat.as_slice().unwrap();
        let c_s = cmat.as_slice().unwrap();

        // delta = softplus(W_dt @ dt8 + b_dt) (dt_proj bias added exactly once, here)
        let mut delta = dot_l(&self.dt_w, &dt8); // (236, L)
        delta
            .as_slice_mut()
            .unwrap()
            .par_chunks_mut(l)
            .enumerate()
            .for_each(|(d, ch)| {
                let b = self.dt_b[d];
                ch.iter_mut().for_each(|v| *v = softplus(*v + b));
            });
        let delta_s = delta.as_slice().unwrap();
        let u_s = u.as_slice().unwrap();
        let a_s = self.a.as_slice().unwrap();
        let z_s = zc.as_slice().unwrap();

        // selective scan (parallel across channels): h = deltaA*h + deltaB_u; y = C*h
        let mut out = Array2::<f32>::zeros((D_INNER, l));
        out.as_slice_mut()
            .unwrap()
            .par_chunks_mut(l)
            .enumerate()
            .for_each(|(d, outd)| {
                let a_row = &a_s[d * D_STATE..(d + 1) * D_STATE];
                let delta_c = &delta_s[d * l..(d + 1) * l];
                let u_c = &u_s[d * l..(d + 1) * l];
                let z_c = &z_s[d * l..(d + 1) * l];
                let dd = self.d[d];
                let mut h = [0f32; D_STATE];
                for i in 0..l {
                    let dl = delta_c[i];
                    let uu = u_c[i];
                    for n in 0..D_STATE {
                        h[n] = (dl * a_row[n]).exp() * h[n] + dl * b_s[n * l + i] * uu;
                    }
                    let mut acc = 0f32;
                    for n in 0..D_STATE {
                        acc += c_s[n * l + i] * h[n];
                    }
                    // out = (y + D·u) ∘ silu(z)
                    outd[i] = (acc + dd * uu) * silu(z_c[i]);
                }
            });

        // out_proj: out (L,118) = (W_out @ out)^T
        dot_l(&self.out_proj, &out).reversed_axes()
    }
}


// =========================================================================
// Batched forward (one batch per forward; faer GEMM + AVX-512 vectorized scan)
// =========================================================================
//
// Layout contract (the key to zero-copy):
// - External input/output: site-major (B*L, D_MODEL) row-major.
// - in_proj writes directly into channel-major (2*D_INNER, B*L) row-major:
//   dst = W(472,118) @ x^T(118, B*L); x^T is a transpose view, no manual transpose.
// - conv / scan read a contiguous L-length row per (d, b) from the channel-major buffers.
// - out_proj: dst(B*L,118) = out^T(B*L,236) @ W_out^T(236,118), two transpose views.
// - faer Par::Seq: single-threaded GEMM (parallelism lives in batch-level rayon); fixed shapes -> deterministic.

use crate::caduceus::D_MODEL;
use faer::linalg::matmul::matmul;
use faer::mat::{MatMut, MatRef};
use faer::{Accum, Par};

/// Working buffers for the batched forward (one per rayon worker, reused across batches to avoid repeated large allocations).
#[derive(Default)]
pub struct BatchBufs {
    xz: Vec<f32>,    // (2*D_INNER, B*L) channel-major
    u: Vec<f32>,     // (D_INNER, B*L)
    x_dbl: Vec<f32>, // (DT_RANK+2*D_STATE, B*L)
    delta: Vec<f32>, // (D_INNER, B*L)
    out: Vec<f32>,   // (D_INNER, B*L)
    bt: Vec<f32>,    // (B*L, D_STATE) contiguous within a site, for SIMD loads
    ct: Vec<f32>,    // (B*L, D_STATE)
    y: Vec<f32>,     // (B*L, D_MODEL) site-major output
}

impl BatchBufs {
    fn resize(&mut self, bl: usize) {
        self.xz.resize(2 * D_INNER * bl, 0.0);
        self.u.resize(D_INNER * bl, 0.0);
        self.x_dbl.resize((DT_RANK + 2 * D_STATE) * bl, 0.0);
        self.delta.resize(D_INNER * bl, 0.0);
        self.out.resize(D_INNER * bl, 0.0);
        self.bt.resize(bl * D_STATE, 0.0);
        self.ct.resize(bl * D_STATE, 0.0);
        self.y.resize(bl * D_MODEL, 0.0);
    }
}

/// Row views for a single scan channel (packed to reduce the parameter count).
struct ScanRows<'a> {
    delta: &'a [f32],
    u: &'a [f32],
    /// Gate values with silu(z) precomputed
    z: &'a [f32],
    bt: &'a [f32],
    ct: &'a [f32],
    a: &'a [f32],
    dd: f32,
}

/// Scan a single channel (scalar fallback; mirrors the single-sequence inner loop step by step).
#[allow(clippy::needless_range_loop)] // multiple slices indexed in lockstep; an iterator rewrite would be less readable
fn scan_channel_scalar(r: &ScanRows, out: &mut [f32], l: usize) {
    let mut h = [0f32; D_STATE];
    for i in 0..l {
        let dl = r.delta[i];
        let uu = r.u[i];
        for n in 0..D_STATE {
            h[n] = (dl * r.a[n]).exp() * h[n] + dl * r.bt[i * D_STATE + n] * uu;
        }
        let mut acc = 0f32;
        for n in 0..D_STATE {
            acc += r.ct[i * D_STATE + n] * h[n];
        }
        out[i] = (acc + r.dd * uu) * r.z[i]; // z already has silu applied
    }
}

/// Scan a single channel (AVX-512: 16 states = one vector; exp uses the vector polynomial).
/// Same formula and order as the scalar version element by element (only the exp implementation differs, ~1 ulp).
/// # Safety
/// The caller must guarantee the CPU supports AVX-512F (gated by is_x86_feature_detected!).
#[cfg(target_arch = "x86_64")]
#[allow(clippy::needless_range_loop)] // SIMD loads need raw indices
#[target_feature(enable = "avx512f")]
unsafe fn scan_channel_avx512(r: &ScanRows, out: &mut [f32], l: usize) {
    use std::arch::x86_64::*;
    unsafe {
        let av = _mm512_loadu_ps(r.a.as_ptr());
        let mut h = _mm512_setzero_ps();
        for i in 0..l {
            let dl = r.delta[i];
            let uu = r.u[i];
            let e = crate::simdexp::exp16_avx512(_mm512_mul_ps(_mm512_set1_ps(dl), av));
            let coeff = _mm512_set1_ps(dl * uu);
            let bv = _mm512_loadu_ps(r.bt.as_ptr().add(i * D_STATE));
            h = _mm512_fmadd_ps(e, h, _mm512_mul_ps(coeff, bv));
            let cv = _mm512_loadu_ps(r.ct.as_ptr().add(i * D_STATE));
            let acc = _mm512_reduce_add_ps(_mm512_mul_ps(cv, h));
            out[i] = (acc + r.dd * uu) * r.z[i]; // z already has silu applied
        }
    }
}

fn scan_channel(r: &ScanRows, out: &mut [f32], l: usize) {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx512f") {
        unsafe {
            scan_channel_avx512(r, out, l);
        }
        return;
    }
    scan_channel_scalar(r, out, l);
}

/// Four-channel interleaved scan (ILP): 4 consecutive bi channels of the same d advance in one L loop.
/// Each channel's floating-point operation sequence is identical to scan_channel (same formula, same
/// order, same exp implementation); only the scheduling interleaves to hide the exp/fmadd latency
/// chain -- output is bit-identical.
fn scan4(r: [&ScanRows; 4], outs: [&mut [f32]; 4], l: usize) {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx512f") {
        unsafe {
            scan4_avx512(&r, outs, l);
        }
        return;
    }
    for (ri, o) in r.iter().zip(outs) {
        scan_channel_scalar(ri, o, l);
    }
}

/// AVX-512 implementation of scan4: 4 independent h chains interleaved. Same d -> shared av vector.
/// # Safety
/// The caller must guarantee the CPU supports AVX-512F.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
#[allow(clippy::needless_range_loop)] // 4-channel ILP kernel; index-based port kept frozen
unsafe fn scan4_avx512(r: &[&ScanRows; 4], outs: [&mut [f32]; 4], l: usize) {
    use std::arch::x86_64::*;
    unsafe {
        let av = _mm512_loadu_ps(r[0].a.as_ptr()); // same d, shared by all four channels
        let mut h0 = _mm512_setzero_ps();
        let mut h1 = _mm512_setzero_ps();
        let mut h2 = _mm512_setzero_ps();
        let mut h3 = _mm512_setzero_ps();
        for i in 0..l {
            let dl0 = r[0].delta[i];
            let uu0 = r[0].u[i];
            let dl1 = r[1].delta[i];
            let uu1 = r[1].u[i];
            let dl2 = r[2].delta[i];
            let uu2 = r[2].u[i];
            let dl3 = r[3].delta[i];
            let uu3 = r[3].u[i];
            let e0 = crate::simdexp::exp16_avx512(_mm512_mul_ps(_mm512_set1_ps(dl0), av));
            let e1 = crate::simdexp::exp16_avx512(_mm512_mul_ps(_mm512_set1_ps(dl1), av));
            let e2 = crate::simdexp::exp16_avx512(_mm512_mul_ps(_mm512_set1_ps(dl2), av));
            let e3 = crate::simdexp::exp16_avx512(_mm512_mul_ps(_mm512_set1_ps(dl3), av));
            let bv0 = _mm512_loadu_ps(r[0].bt.as_ptr().add(i * D_STATE));
            let bv1 = _mm512_loadu_ps(r[1].bt.as_ptr().add(i * D_STATE));
            let bv2 = _mm512_loadu_ps(r[2].bt.as_ptr().add(i * D_STATE));
            let bv3 = _mm512_loadu_ps(r[3].bt.as_ptr().add(i * D_STATE));
            h0 = _mm512_fmadd_ps(e0, h0, _mm512_mul_ps(_mm512_set1_ps(dl0 * uu0), bv0));
            h1 = _mm512_fmadd_ps(e1, h1, _mm512_mul_ps(_mm512_set1_ps(dl1 * uu1), bv1));
            h2 = _mm512_fmadd_ps(e2, h2, _mm512_mul_ps(_mm512_set1_ps(dl2 * uu2), bv2));
            h3 = _mm512_fmadd_ps(e3, h3, _mm512_mul_ps(_mm512_set1_ps(dl3 * uu3), bv3));
            let cv0 = _mm512_loadu_ps(r[0].ct.as_ptr().add(i * D_STATE));
            let cv1 = _mm512_loadu_ps(r[1].ct.as_ptr().add(i * D_STATE));
            let cv2 = _mm512_loadu_ps(r[2].ct.as_ptr().add(i * D_STATE));
            let cv3 = _mm512_loadu_ps(r[3].ct.as_ptr().add(i * D_STATE));
            let acc0 = _mm512_reduce_add_ps(_mm512_mul_ps(cv0, h0));
            let acc1 = _mm512_reduce_add_ps(_mm512_mul_ps(cv1, h1));
            let acc2 = _mm512_reduce_add_ps(_mm512_mul_ps(cv2, h2));
            let acc3 = _mm512_reduce_add_ps(_mm512_mul_ps(cv3, h3));
            outs[0][i] = (acc0 + r[0].dd * uu0) * r[0].z[i];
            outs[1][i] = (acc1 + r[1].dd * uu1) * r[1].z[i];
            outs[2][i] = (acc2 + r[2].dd * uu2) * r[2].z[i];
            outs[3][i] = (acc3 + r[3].dd * uu3) * r[3].z[i];
        }
    }
}

impl MambaWeights {
    /// Batched forward: x (B*L, D_MODEL) site-major -> y (B*L, D_MODEL) site-major.
    /// Single-threaded (faer Par::Seq); parallelism is provided by batch-level rayon above. fp32 throughout.
    pub fn forward_batch(&self, x: &[f32], b: usize, l: usize, bufs: &mut BatchBufs) {
        let bl = b * l;
        debug_assert_eq!(x.len(), bl * D_MODEL);
        bufs.resize(bl);
        let BatchBufs { xz, u, x_dbl, delta, out, bt, ct, y } = bufs;

        // xz = in_proj @ x^T → (472, BL) channel-major
        // With ESPERANTO_BF16=1 use the bf16 microkernel (Zen4 VDPBF16PS); otherwise faer f32.
        if !crate::bf16gemm::in_proj_bf16(x, self.in_proj.as_slice().unwrap(), bl, xz) {
            let dst = MatMut::from_row_major_slice_mut(xz.as_mut_slice(), 2 * D_INNER, bl);
            let lhs = MatRef::from_row_major_slice(self.in_proj.as_slice().unwrap(), 2 * D_INNER, D_MODEL);
            let rhs = MatRef::from_row_major_slice(x, bl, D_MODEL).transpose();
            matmul(dst, Accum::Replace, lhs, rhs, 1.0f32, Par::Seq);
        }
        let (xc, zc) = xz.split_at_mut(D_INNER * bl);
        // silu of z pre-vectorized (used directly by the scan gate; differs from elementwise libm by ~1 ulp)
        crate::simdexp::silu_slice_inplace(zc);

        // causal conv1d (k=4, left-pad 3) + bias (pure multiply-add; silu is vectorized uniformly afterwards).
        // The 3 boundary points use the k-loop (term-by-term identical to the single-sequence version); the main loop is branch-free and auto-vectorizable.
        for d in 0..D_INNER {
            let w = [
                self.conv_w[(d, 0)],
                self.conv_w[(d, 1)],
                self.conv_w[(d, 2)],
                self.conv_w[(d, 3)],
            ];
            let bias = self.conv_b[d];
            for bi in 0..b {
                let base = d * bl + bi * l;
                let xrow = &xc[base..base + l];
                let urow = &mut u[base..base + l];
                let edge = (D_CONV - 1).min(l);
                for i in 0..edge {
                    let mut acc = bias;
                    for (k, wk) in w.iter().enumerate() {
                        if i + k >= D_CONV - 1 {
                            acc += wk * xrow[i + k - (D_CONV - 1)];
                        }
                    }
                    urow[i] = acc;
                }
                for i in edge..l {
                    urow[i] = bias
                        + w[0] * xrow[i - 3]
                        + w[1] * xrow[i - 2]
                        + w[2] * xrow[i - 1]
                        + w[3] * xrow[i];
                }
            }
        }
        crate::simdexp::silu_slice_inplace(u);

        // x_dbl = x_proj @ u → (40, BL)
        if !crate::bf16gemm::x_proj_bf16(u, self.x_proj.as_slice().unwrap(), bl, x_dbl) {
            let dst = MatMut::from_row_major_slice_mut(x_dbl.as_mut_slice(), DT_RANK + 2 * D_STATE, bl);
            let lhs = MatRef::from_row_major_slice(self.x_proj.as_slice().unwrap(), DT_RANK + 2 * D_STATE, D_INNER);
            let rhs = MatRef::from_row_major_slice(u.as_slice(), D_INNER, bl);
            matmul(dst, Accum::Replace, lhs, rhs, 1.0f32, Par::Seq);
        }
        let (dt8, bc) = x_dbl.split_at(DT_RANK * bl);
        let (bmat, cmat) = bc.split_at(D_STATE * bl);

        // delta = softplus(dt_w @ dt8 + b_dt)
        // GEMM picks one of two paths (bf16 microkernel / faer f32); bias+softplus is shared by both.
        if !crate::bf16gemm::dt_bf16(dt8, self.dt_w.as_slice().unwrap(), bl, delta) {
            let dst = MatMut::from_row_major_slice_mut(delta.as_mut_slice(), D_INNER, bl);
            let lhs = MatRef::from_row_major_slice(self.dt_w.as_slice().unwrap(), D_INNER, DT_RANK);
            let rhs = MatRef::from_row_major_slice(dt8, DT_RANK, bl);
            matmul(dst, Accum::Replace, lhs, rhs, 1.0f32, Par::Seq);
        }
        for d in 0..D_INNER {
            let bias = self.dt_b[d];
            let row = &mut delta[d * bl..(d + 1) * bl];
            for v in row.iter_mut() {
                *v += bias;
            }
        }
        crate::simdexp::softplus_slice_inplace(delta);

        // Transpose B/C to (L,16) contiguous within a site (for SIMD)
        for j in 0..bl {
            for n in 0..D_STATE {
                bt[j * D_STATE + n] = bmat[n * bl + j];
                ct[j * D_STATE + n] = cmat[n * bl + j];
            }
        }

        // selective scan: each (d, b) is independent; same-d groups of 4 consecutive bi are interleaved (ILP, bit-identical output)
        let a_s = self.a.as_slice().unwrap();
        for d in 0..D_INNER {
            let a_row = &a_s[d * D_STATE..(d + 1) * D_STATE];
            let dd = self.d[d];
            let dbase = d * bl;
            let mk = |bi: usize| ScanRows {
                delta: &delta[dbase + bi * l..dbase + (bi + 1) * l],
                u: &u[dbase + bi * l..dbase + (bi + 1) * l],
                z: &zc[dbase + bi * l..dbase + (bi + 1) * l],
                bt: &bt[bi * l * D_STATE..(bi + 1) * l * D_STATE],
                ct: &ct[bi * l * D_STATE..(bi + 1) * l * D_STATE],
                a: a_row,
                dd,
            };
            let out_chan = &mut out[dbase..dbase + bl];
            let mut quads = out_chan.chunks_exact_mut(4 * l);
            let mut bi = 0;
            for quad in &mut quads {
                let (o0, r) = quad.split_at_mut(l);
                let (o1, r) = r.split_at_mut(l);
                let (o2, o3) = r.split_at_mut(l);
                let (r0, r1, r2, r3) = (mk(bi), mk(bi + 1), mk(bi + 2), mk(bi + 3));
                scan4([&r0, &r1, &r2, &r3], [o0, o1, o2, o3], l);
                bi += 4;
            }
            for rem in quads.into_remainder().chunks_exact_mut(l) {
                scan_channel(&mk(bi), rem, l);
                bi += 1;
            }
        }

        // y = out^T @ out_proj^T → (BL, 118) site-major
        {
            let dst = MatMut::from_row_major_slice_mut(y.as_mut_slice(), bl, D_MODEL);
            let lhs = MatRef::from_row_major_slice(out.as_slice(), D_INNER, bl).transpose();
            let rhs = MatRef::from_row_major_slice(self.out_proj.as_slice().unwrap(), D_MODEL, D_INNER).transpose();
            matmul(dst, Accum::Replace, lhs, rhs, 1.0f32, Par::Seq);
        }
    }

    /// Batched output (valid after the last forward_batch call).
    pub fn batch_output(bufs: &BatchBufs) -> &[f32] {
        &bufs.y
    }
}
