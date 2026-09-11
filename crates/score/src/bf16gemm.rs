//! bf16 GEMM microkernels (Zen4 AVX512-BF16, VDPBF16PS).
//!
//! Two kernels, each taking the natural form for the actual memory layouts of the
//! four GEMMs in forward_batch:
//! - dot kernel: C[i,j] = Sum_k A[i,k]*B[j,k] (both A and B are K-contiguous) -- used for in_proj
//!   (W(472,118) row-major x x(BL,118) row-major -> xz(472,BL) row-major).
//! - axpy kernel: C[i,n] = Sum_k A[i,k]*B[k,n] (B row-contiguous) -- used for x_proj/dt
//!   (W(M,K) x activations(K,BL) row-major). B is pre-packed in a k-pair interleaved layout (see below).
//!
//! VDPBF16PS semantics: acc[f32 lane l] += a.bf16[2l]*b.bf16[2l] + a.bf16[2l+1]*b.bf16[2l+1].
//! - dot: lane = partial sums of output (i,j) (k pairs spread across 16 lanes, horizontal reduction at the end);
//!   K is zero-padded to a multiple of 32 (weights/activations are zero-padded at packing time; zero elements contribute nothing).
//! - axpy: lane = output n (16 per vector); a = broadcast of the pair (A[i][2p], A[i][2p+1]),
//!   b = (B[2p][n], B[2p+1][n]) pre-packed interleaved -- B_pack[p][n][2] layout.
//!
//! out_proj stays on faer f32 (the transpose cost of its (BL,118) output layout eats the gains; see bench review).
//! Numerics: bf16 inputs (8-bit mantissa, ~0.4% element-level relative error) + f32 accumulation; K<=236 ->
//! typical dot-product relative error <1%; after mean-pool (502 sites) the partially independent errors cancel,
//! so embedding perturbation is on the same order as fp16 quantization noise -- guarded by the embed
//! equivalence test and the probability-drift probe.

// Port frozen: keep this kernel verbatim; newer clippy style lints are exempted module-wide instead of rewriting the SIMD code.
#![allow(clippy::manual_is_multiple_of, clippy::missing_transmute_annotations, clippy::needless_range_loop, clippy::manual_div_ceil)]

use crate::caduceus::D_MODEL;
use crate::mamba::{D_INNER, DT_RANK, D_STATE};

/// f32 -> bf16 bit pattern (RNE rounding). Returns u16 (bits).
#[inline]
pub fn f32_to_bf16_bits(f: f32) -> u16 {
    let u = f.to_bits();
    let bias = 0x7FFFu32 + ((u >> 16) & 1);
    ((u + bias) >> 16) as u16
}

/// Row-major f32 -> bf16, with K zero-padded to k_pad (a multiple of 32). Output is (rows, k_pad) row-major.
pub fn pack_rows_padded(src: &[f32], rows: usize, k: usize, k_pad: usize) -> Vec<u16> {
    debug_assert!(k_pad % 32 == 0 && k <= k_pad);
    let mut out = vec![0u16; rows * k_pad];
    for r in 0..rows {
        let s = &src[r * k..r * k + k];
        let d = &mut out[r * k_pad..r * k_pad + k];
        for (d, &v) in d.iter_mut().zip(s.iter()) {
            *d = f32_to_bf16_bits(v);
        }
    }
    out
}

/// Row-major f32 (k_rows, n) -> k-pair interleaved bf16 layout (k_rows/2, n, 2):
/// out[p][n][2] = {bf16(B[2p][n]), bf16(B[2p+1][n])}. k_rows must be even.
pub fn pack_pair_interleaved(src: &[f32], k_rows: usize, n: usize) -> Vec<u16> {
    debug_assert!(k_rows % 2 == 0);
    let pairs = k_rows / 2;
    let mut out = vec![0u16; pairs * n * 2];
    for p in 0..pairs {
        let r0 = &src[(2 * p) * n..(2 * p + 1) * n];
        let r1 = &src[(2 * p + 1) * n..(2 * p + 2) * n];
        let d = &mut out[p * n * 2..(p + 1) * n * 2];
        for j in 0..n {
            d[2 * j] = f32_to_bf16_bits(r0[j]);
            d[2 * j + 1] = f32_to_bf16_bits(r1[j]);
        }
    }
    out
}

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

/// dot kernel: out[i*n+j] = Sum_k a[i,k]*b[j,k].
/// a: (m, k_pad) row-major bf16; b: (n, k_pad) row-major bf16; out: (m, n) f32.
/// 4-way output parallelism (4 independent accumulation chains to fill the ports).
/// # Safety
/// The caller must guarantee the CPU supports AVX512-BF16 and that a/b lengths match the declared dimensions.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bf16")]
pub unsafe fn gemm_dot_bf16(
    a: &[u16],
    b: &[u16],
    out: &mut [f32],
    m: usize,
    n: usize,
    k_pad: usize,
) {
    debug_assert!(a.len() == m * k_pad && b.len() == n * k_pad && out.len() == m * n);
    let kchunks = k_pad / 32;
    for i in 0..m {
        let arow = &a[i * k_pad..];
        let orow = &mut out[i * n..];
        let mut j = 0usize;
        while j + 4 <= n {
            // independent accumulation chains for 4 outputs
            let mut acc0 = _mm512_setzero_ps();
            let mut acc1 = _mm512_setzero_ps();
            let mut acc2 = _mm512_setzero_ps();
            let mut acc3 = _mm512_setzero_ps();
            for kc in 0..kchunks {
                let av = _mm512_loadu_si512(arow[kc * 32..].as_ptr() as *const _);
                let b0 = _mm512_loadu_si512(b[(j) * k_pad + kc * 32..].as_ptr() as *const _);
                let b1 = _mm512_loadu_si512(b[(j + 1) * k_pad + kc * 32..].as_ptr() as *const _);
                let b2 = _mm512_loadu_si512(b[(j + 2) * k_pad + kc * 32..].as_ptr() as *const _);
                let b3 = _mm512_loadu_si512(b[(j + 3) * k_pad + kc * 32..].as_ptr() as *const _);
                let av: __m512bh = core::mem::transmute(av);
                acc0 = _mm512_dpbf16_ps(acc0, av, core::mem::transmute(b0));
                acc1 = _mm512_dpbf16_ps(acc1, av, core::mem::transmute(b1));
                acc2 = _mm512_dpbf16_ps(acc2, av, core::mem::transmute(b2));
                acc3 = _mm512_dpbf16_ps(acc3, av, core::mem::transmute(b3));
            }
            orow[j] = _mm512_reduce_add_ps(acc0);
            orow[j + 1] = _mm512_reduce_add_ps(acc1);
            orow[j + 2] = _mm512_reduce_add_ps(acc2);
            orow[j + 3] = _mm512_reduce_add_ps(acc3);
            j += 4;
        }
        for jj in j..n {
            let mut acc = _mm512_setzero_ps();
            for kc in 0..kchunks {
                let av = _mm512_loadu_si512(arow[kc * 32..].as_ptr() as *const _);
                let bv = _mm512_loadu_si512(b[jj * k_pad + kc * 32..].as_ptr() as *const _);
                acc = _mm512_dpbf16_ps(acc, core::mem::transmute(av), core::mem::transmute(bv));
            }
            orow[jj] = _mm512_reduce_add_ps(acc);
        }
    }
}

/// axpy kernel: out[i,n] = Sum_k a[i,k]*b[k,n].
/// a: (m, k) row-major **f32** (weights, packed into scalar pairs inline; small);
/// b: output of pack_pair_interleaved, (k/2, n, 2); out: (m, n) f32 row-major.
/// # Safety
/// The caller must guarantee the CPU supports AVX512-BF16 and that dimensions match.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bf16")]
pub unsafe fn gemm_axpy_bf16(a: &[f32], b: &[u16], out: &mut [f32], m: usize, n: usize, k: usize) {
    debug_assert!(k % 2 == 0 && a.len() == m * k && b.len() == k / 2 * n * 2 && out.len() == m * n);
    let pairs = k / 2;
    let nvec = n / 16; // n is BL (even, multiple of 64); scalar fallback for the tail
    for i in 0..m {
        let arow = &a[i * k..];
        let orow = &mut out[i * n..];
        // one accumulator per group of 16 outputs; pair count k is usually > 64 -> the accumulation chain is long enough that rotating groups fills the ports
        for nv in 0..nvec {
            let mut acc = _mm512_setzero_ps();
            for p in 0..pairs {
                // VDPBF16PS per 32-bit lane: low16(a)xlow16(b) + high16(a)xhigh16(b).
                // B interleaved layout b[p][n][2]: output group nv (16 columns) of pair p starts at p*(n*2) + nv*32.
                // A low16 = A[i][2p], high16 = A[i][2p+1] (paired with B's (2p, 2p+1)).
                let alo = f32_to_bf16_bits(arow[2 * p]) as u32;
                let ahi = f32_to_bf16_bits(arow[2 * p + 1]) as u32;
                let av = _mm512_set1_epi32(((ahi << 16) | alo) as i32);
                let base = p * n * 2 + nv * 32;
                let bv = _mm512_loadu_si512(b[base..].as_ptr() as *const _);
                acc = _mm512_dpbf16_ps(acc, core::mem::transmute(av), core::mem::transmute(bv));
            }
            let mut tmp = [0f32; 16];
            _mm512_storeu_ps(tmp.as_mut_ptr(), acc);
            orow[nv * 16..nv * 16 + 16].copy_from_slice(&tmp);
        }
        // tail n%16: scalar fallback (last few columns when bl is not a multiple of 16)
        for j in (nvec * 16)..n {
            let mut acc = 0f32;
            for t in 0..k {
                acc += arow[t] * bf16_at(b, t, j, n);
            }
            orow[j] = acc;
        }
    }
}

/// Read from the interleaved layout: the f32 value (bf16-quantized) of logical element (k=t, column j) in b[p][n][2].
#[inline]
fn bf16_at(b: &[u16], t: usize, j: usize, n: usize) -> f32 {
    let bits = b[(t / 2) * n * 2 + j * 2 + (t % 2)];
    f32::from_bits((bits as u32) << 16)
}

/// Whether the CPU supports AVX512-BF16 (runtime detection).
pub fn bf16_supported() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::arch::is_x86_feature_detected!("avx512bf16")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// Switch for the bf16 path in forward_batch (off by default; enable with ESPERANTO_BF16).
pub fn bf16_enabled() -> bool {
    bf16_supported() && std::env::var_os("ESPERANTO_BF16").is_some()
}

/// in_proj hook: xz(472, BL) = W_in(472,118) @ x(BL,118)^T, dot kernel.
/// false = not enabled (caller falls back to the original faer f32 path).
pub fn in_proj_bf16(x: &[f32], w_in: &[f32], bl: usize, xz: &mut [f32]) -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        if !bf16_enabled() {
            return false;
        }
        let k_pad = (D_MODEL + 31) / 32 * 32; // 118 → 128
        let a = pack_rows_padded(w_in, 2 * D_INNER, D_MODEL, k_pad);
        let bx = pack_rows_padded(x, bl, D_MODEL, k_pad);
        unsafe { gemm_dot_bf16(&a, &bx, xz, 2 * D_INNER, bl, k_pad) };
        true
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let _ = (x, w_in, bl, xz);
        false
    }
}

/// x_proj hook: x_dbl(40, BL) = W_xp(40,236) @ u(236,BL), axpy kernel.
pub fn x_proj_bf16(u: &[f32], w_xp: &[f32], bl: usize, x_dbl: &mut [f32]) -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        if !bf16_enabled() {
            return false;
        }
        let m = DT_RANK + 2 * D_STATE; // 40
        let bp = pack_pair_interleaved(u, D_INNER, bl);
        unsafe { gemm_axpy_bf16(w_xp, &bp, x_dbl, m, bl, D_INNER) };
        true
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let _ = (u, w_xp, bl, x_dbl);
        false
    }
}

/// dt hook: delta(236, BL) = W_dt(236,8) @ dt8(8,BL), axpy kernel (DT_RANK=8, even).
pub fn dt_bf16(dt8: &[f32], w_dt: &[f32], bl: usize, delta: &mut [f32]) -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        if !bf16_enabled() {
            return false;
        }
        let bp = pack_pair_interleaved(dt8, DT_RANK, bl);
        unsafe { gemm_axpy_bf16(w_dt, &bp, delta, D_INNER, bl, DT_RANK) };
        true
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let _ = (dt8, w_dt, bl, delta);
        false
    }
}

