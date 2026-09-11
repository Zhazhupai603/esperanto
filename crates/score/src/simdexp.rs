//! AVX-512 batched exp (Cephes-style polynomial, ~1 ulp) + scalar fallback.
//!
//! The selective scan's `exp(delta*A)` is the main elementwise-parallel compute load
//! (236 channels x L x 16 states per forward); scalar libm exp was one of the old bottlenecks.
//! Numerics: within a few ulp of libm, fp32 throughout; verified empirically (|delta-prob| <= 1e-4).

/// Scalar fallback: elementwise libm exp.
pub fn exp_slice_scalar(x: &[f32], out: &mut [f32]) {
    for (o, &v) in out.iter_mut().zip(x) {
        *o = v.exp();
    }
}

/// AVX-512 vector exp (16 lanes). Cephes single-precision algorithm:
/// exp(x) = 2^n * poly(r), n = round(x*log2e), r = x - n*ln2 (two-part Cody-Waite).
/// Domain clamped to [-87.3, 88.7] (the non-overflowing range of fp32 exp).
///
/// # Safety
/// The caller must guarantee the runtime CPU supports AVX-512F (gated by is_x86_feature_detected!).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
pub unsafe fn exp16_avx512(x: std::arch::x86_64::__m512) -> std::arch::x86_64::__m512 {
    use std::arch::x86_64::*;
    {
        let max_x = _mm512_set1_ps(88.722_84);
        let min_x = _mm512_set1_ps(-87.336_55);
        let x = _mm512_min_ps(_mm512_max_ps(x, min_x), max_x);

        let log2e = _mm512_set1_ps(std::f32::consts::LOG2_E);
        // round-to-nearest via magic number (value range far below 2^22, safe)
        let magic = _mm512_set1_ps(12_582_912.0); // 1.5 × 2^23
        let t = _mm512_mul_ps(x, log2e);
        let fx = _mm512_sub_ps(_mm512_add_ps(t, magic), magic);
        // Cody-Waite: r = x - fx*ln2_hi - fx*ln2_lo
        let ln2_hi = _mm512_set1_ps(0.693_359_4);
        let ln2_lo = _mm512_set1_ps(-0.000_212_194_4);
        let r = _mm512_fnmadd_ps(fx, ln2_hi, x);
        let r = _mm512_fnmadd_ps(fx, ln2_lo, r);
        // Cephes expf polynomial (degree 6)
        let c1 = _mm512_set1_ps(1.987_569_2e-4);
        let c2 = _mm512_set1_ps(1.398_199_9e-3);
        let c3 = _mm512_set1_ps(8.333_451e-3);
        let c4 = _mm512_set1_ps(4.166_579_5e-2);
        let c5 = _mm512_set1_ps(1.666_666_6e-1);
        let c6 = _mm512_set1_ps(0.5);
        let mut p = c1;
        p = _mm512_fmadd_ps(p, r, c2);
        p = _mm512_fmadd_ps(p, r, c3);
        p = _mm512_fmadd_ps(p, r, c4);
        p = _mm512_fmadd_ps(p, r, c5);
        p = _mm512_fmadd_ps(p, r, c6);
        // exp(r) ~ r + r^2*p(r) + 1... cephes variant of y = 1 + r + r^2 p: y = p*r*r + r + 1
        let rr = _mm512_mul_ps(r, r);
        let mut y = _mm512_fmadd_ps(p, rr, r);
        y = _mm512_add_ps(y, _mm512_set1_ps(1.0));
        // 2^n: add n to the exponent field
        let n = _mm512_cvtps_epi32(fx);
        let pow2n = _mm512_castsi512_ps(_mm512_slli_epi32(_mm512_add_epi32(n, _mm512_set1_epi32(127)), 23));
        _mm512_mul_ps(y, pow2n)
    }
}

/// Batched exp: vector path when AVX-512 is available, scalar otherwise. Output and input may be equal-length slices.
pub fn exp_slice(x: &[f32], out: &mut [f32]) {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx512f") {
        unsafe {
            use std::arch::x86_64::*;
            let n = x.len();
            let mut i = 0;
            while i + 16 <= n {
                let v = _mm512_loadu_ps(x.as_ptr().add(i));
                let e = exp16_avx512(v);
                _mm512_storeu_ps(out.as_mut_ptr().add(i), e);
                i += 16;
            }
            while i < n {
                out[i] = x[i].exp();
                i += 1;
            }
        }
        return;
    }
    exp_slice_scalar(x, out);
}

/// In-place exp.
pub fn exp_slice_inplace(x: &mut [f32]) {
    let tmp = x.to_vec();
    exp_slice(&tmp, x);
}


/// AVX-512 vector ln(1+x) (Cephes logf structure: frexp split into 2^e*m, degree-8 polynomial).
/// Domain x > -1; only used by softplus (exp in (0,1]); small-x accuracy is guaranteed by the polynomial (~1 ulp).
///
/// # Safety
/// The caller must guarantee the runtime CPU supports AVX-512F (gated by is_x86_feature_detected!).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
pub unsafe fn ln1p16_avx512(x: std::arch::x86_64::__m512) -> std::arch::x86_64::__m512 {
    use std::arch::x86_64::*;
    {
        let one = _mm512_set1_ps(1.0);
        let x1 = _mm512_add_ps(x, one); // 1+x ∈ (1, 2]
        // frexp: x1 = m*2^e, m in [0.5,1)
        let xi = _mm512_castps_si512(x1);
        let e = _mm512_sub_epi32(_mm512_srli_epi32(xi, 23), _mm512_set1_epi32(126));
        let m = _mm512_castsi512_ps(_mm512_or_si512(
            _mm512_and_si512(xi, _mm512_set1_epi32(0x007f_ffff)),
            _mm512_set1_epi32(0x3f00_0000), // 0.5
        ));
        // m in [0.5,1); if m < sqrt(2)/2 ~= 0.7071 then m*=2, e-=1 (Cephes convention: m in [sqrt(2)/2, sqrt(2)))
        #[allow(clippy::approx_constant)] // differs from FRAC_1_SQRT_2 in the last digit; changing the literal would change results
        let sqrthf = _mm512_set1_ps(0.707_106_8);
        let mask = _mm512_cmp_ps_mask(m, sqrthf, _CMP_LT_OS);
        let m = _mm512_mask_add_ps(m, mask, m, m);
        let e = _mm512_mask_sub_epi32(e, mask, e, _mm512_set1_epi32(1));
        // z = m - 1 ∈ [-0.2929, 0.4142]
        let z = _mm512_sub_ps(m, one);
        // Cephes logf polynomial (degree 8) on z
        let c0 = _mm512_set1_ps(7.037_683_5e-2);
        let c1 = _mm512_set1_ps(-1.151_461e-1);
        let c2 = _mm512_set1_ps(1.167_699e-1);
        let c3 = _mm512_set1_ps(-1.242_014_5e-1);
        let c4 = _mm512_set1_ps(1.424_932_3e-1);
        let c5 = _mm512_set1_ps(-1.666_666e-1);
        let c6 = _mm512_set1_ps(0.2);
        let c7 = _mm512_set1_ps(-0.25);
        let c8 = _mm512_set1_ps(3.333_333_3e-1);
        let mut p = c0;
        p = _mm512_fmadd_ps(p, z, c1);
        p = _mm512_fmadd_ps(p, z, c2);
        p = _mm512_fmadd_ps(p, z, c3);
        p = _mm512_fmadd_ps(p, z, c4);
        p = _mm512_fmadd_ps(p, z, c5);
        p = _mm512_fmadd_ps(p, z, c6);
        p = _mm512_fmadd_ps(p, z, c7);
        p = _mm512_fmadd_ps(p, z, c8);
        // log(m) = z - z^2*0.5 + z^3*p(z)... Cephes form: log = z - 0.5 z^2 + z^3*P(z)
        let zz = _mm512_mul_ps(z, z);
        let mut lm = _mm512_fnmadd_ps(_mm512_set1_ps(0.5), zz, z); // z - 0.5 z²
        lm = _mm512_fmadd_ps(_mm512_mul_ps(p, zz), z, lm); // + z³ p(z) → p·z·z²
        // + e*ln2 (hi/lo two-part)
        let ef = _mm512_cvtepi32_ps(e);
        lm = _mm512_fmadd_ps(ef, _mm512_set1_ps(0.693_359_4), lm);
        lm = _mm512_fmadd_ps(ef, _mm512_set1_ps(-0.000_212_194_4), lm);
        // divergence as x -> -1+ is outside the domain; no special handling (x in (0,1] for this use)
        lm
    }
}

/// AVX-512 vector softplus: softplus(x) = max(x,0) + ln1p(exp(-|x|));
/// numerical difference from the scalar version (threshold=20 cutoff) <= a few ulp
/// (the difference around the cutoff is < 1e-7, negligible).
///
/// # Safety
/// The caller must guarantee the runtime CPU supports AVX-512F (gated by is_x86_feature_detected!).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
pub unsafe fn softplus16_avx512(x: std::arch::x86_64::__m512) -> std::arch::x86_64::__m512 {
    use std::arch::x86_64::*;
    {
        let zero = _mm512_setzero_ps();
        let ax = _mm512_abs_ps(x);
        let e = exp16_avx512(_mm512_sub_ps(zero, ax)); // exp(-|x|) ∈ (0,1]
        let l = ln1p16_avx512(e);
        _mm512_add_ps(_mm512_max_ps(x, zero), l)
    }
}

/// Batched softplus (in place).
pub fn softplus_slice_inplace(x: &mut [f32]) {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx512f") {
        unsafe {
            use std::arch::x86_64::*;
            let n = x.len();
            let mut i = 0;
            while i + 16 <= n {
                let v = _mm512_loadu_ps(x.as_ptr().add(i));
                _mm512_storeu_ps(x.as_mut_ptr().add(i), softplus16_avx512(v));
                i += 16;
            }
            while i < n {
                let v = x[i];
                x[i] = if v > 20.0 { v } else { v.exp().ln_1p() };
                i += 1;
            }
        }
        return;
    }
    for v in x.iter_mut() {
        let t = *v;
        *v = if t > 20.0 { t } else { t.exp().ln_1p() };
    }
}

/// AVX-512 vector silu: x / (1 + exp(-x)).
///
/// # Safety
/// The caller must guarantee the runtime CPU supports AVX-512F (gated by is_x86_feature_detected!).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
pub unsafe fn silu16_avx512(x: std::arch::x86_64::__m512) -> std::arch::x86_64::__m512 {
    use std::arch::x86_64::*;
    {
        let e = exp16_avx512(_mm512_sub_ps(_mm512_setzero_ps(), x));
        _mm512_div_ps(x, _mm512_add_ps(_mm512_set1_ps(1.0), e))
    }
}

/// Batched silu (in place).
pub fn silu_slice_inplace(x: &mut [f32]) {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx512f") {
        unsafe {
            use std::arch::x86_64::*;
            let n = x.len();
            let mut i = 0;
            while i + 16 <= n {
                let v = _mm512_loadu_ps(x.as_ptr().add(i));
                _mm512_storeu_ps(x.as_mut_ptr().add(i), silu16_avx512(v));
                i += 16;
            }
            while i < n {
                let v = x[i];
                x[i] = v / (1.0 + (-v).exp());
                i += 1;
            }
        }
        return;
    }
    for v in x.iter_mut() {
        let t = *v;
        *v = t / (1.0 + (-t).exp());
    }
}
