//! EA-Myers bit-parallel verification (decisive fast lane).
//!
//! Word width is u128: a single block covers patterns up to 128 bases (no
//! truncation), and patterns of 129..=256 bases run through the [`long`]
//! module's 2-block chain. Patterns longer than 256 hit the guard and report
//! their own length as the distance — no verification threshold can pass, so
//! callers fall back to the banded DP.
//!
//! `build_peq` implements the editing-aware merge: read A matches {A}; read C
//! matches {C, T}; read G matches {G, A}; read T matches {T} (A-to-I editing
//! tolerance — read G vs reference A and read C vs reference T are free).

/// Maximum pattern length handled by bit-parallel blocks (2 × 128).
const LONG_MAX: usize = 256;

/// Low `m` bits set; `m >= 128` special-cased to all-ones (avoids `1 << 128`).
#[inline]
pub fn mask_for(m: usize) -> u128 {
    if m >= 128 {
        u128::MAX
    } else {
        (1u128 << m) - 1
    }
}

/// EA equivalence: does this read base match this reference base for free?
#[inline]
fn ea_matches(read: u8, reference: u8) -> bool {
    let r = read.to_ascii_uppercase();
    let f = reference.to_ascii_uppercase();
    match r {
        b'A' => f == b'A',
        b'C' => f == b'C' || f == b'T',
        b'G' => f == b'G' || f == b'A',
        b'T' => f == b'T',
        _ => false,
    }
}

/// Pattern equality table indexed by raw reference byte: bit `i` of
/// `peq[ref_byte]` is set iff `ea_matches(read[i], ref_byte)`.
pub fn build_peq(read: &[u8]) -> [u128; 256] {
    let mut peq = [0u128; 256];
    for (i, &r) in read.iter().enumerate().take(128) {
        for f in b"ACGT" {
            if ea_matches(r, *f) {
                peq[*f as usize] |= 1u128 << i;
            }
        }
    }
    // lowercase reference bytes share their uppercase masks
    for (lo, up) in [(b'a', b'A'), (b'c', b'C'), (b'g', b'G'), (b't', b'T')] {
        peq[lo as usize] = peq[up as usize];
    }
    peq
}

/// One block-advance step of the Myers algorithm, edlib-form carry.
///
/// `hin` is the horizontal delta carried into row 0 (+1/0/−1; block 0 always
/// takes 0, a chained block takes the previous block's `hout`). `mask` trims
/// bits beyond the pattern; `high_bit` selects the block's top pattern
/// position for the outgoing carry. Returns `(pv_out, mv_out, hout)`.
pub fn calculate_block(
    pv: u128,
    mv: u128,
    eq_in: u128,
    hin: i32,
    mask: u128,
    high_bit: u128,
) -> (u128, u128, i32) {
    let hin_is_neg = (hin < 0) as u128;
    let ph_bit0 = (hin > 0) as u128;
    let xv = eq_in | mv;
    let eq = eq_in | hin_is_neg;
    let xh = (((eq & pv).wrapping_add(pv)) ^ pv) | eq;
    let ph = mv | !(xh | pv);
    let mh = pv & xh;
    let hout = if ph & high_bit != 0 {
        1
    } else if mh & high_bit != 0 {
        -1
    } else {
        0
    };
    let ph_sh = (ph << 1) | ph_bit0;
    let mh_sh = (mh << 1) | hin_is_neg;
    let pv_out = (mh_sh | !(xv | ph_sh)) & mask;
    let mv_out = (ph_sh & xv) & mask;
    (pv_out, mv_out, hout)
}

fn scan(
    peq: &[u128; 256],
    m: usize,
    text: &[u8],
    mut f: impl FnMut(&mut i32, &mut usize, i32, usize),
) -> (i32, usize) {
    let mask = mask_for(m);
    let high = 1u128 << (m - 1);
    let mut vp = mask;
    let mut vn = 0u128;
    let mut score = m as i32;
    let mut best = score;
    let mut best_at = 0usize;
    for (j, &c) in text.iter().enumerate() {
        let (nvp, nvn, hout) = calculate_block(vp, vn, peq[c as usize], 0, mask, high);
        vp = nvp;
        vn = nvn;
        score += hout;
        f(&mut best, &mut best_at, score, j + 1);
    }
    (best, best_at)
}

/// Minimum EA edit distance between the read and any substring of `text`
/// (text prefix and suffix free). Patterns longer than 256 return `m`.
pub fn infix(read: &[u8], text: &[u8]) -> i32 {
    let m = read.len();
    if m == 0 {
        return 0;
    }
    if m > LONG_MAX {
        return m as i32;
    }
    let peq = build_peq(read);
    scan(&peq, m, text, |best, _at, score, _j| {
        if score < *best {
            *best = score;
        }
    })
    .0
}

/// Like [`infix`], also reporting the best end point: the number of text
/// bases consumed (exclusive end index; earliest optimum wins ties by the
/// strict-`<` update).
pub fn infix_best_end(read: &[u8], text: &[u8]) -> (i32, usize) {
    let m = read.len();
    if m == 0 {
        return (0, 0);
    }
    if m > LONG_MAX {
        return (m as i32, 0);
    }
    let peq = build_peq(read);
    scan(&peq, m, text, |best, at, score, j| {
        if score < *best {
            *best = score;
            *at = j;
        }
    })
}

/// Best start point: [`infix_best_end`] on the doubly reversed pair; the
/// returned `start` is the inclusive 0-based start of the best-scoring
/// alignment in the original text (`text.len() − rev_end`).
pub fn infix_best_start(read: &[u8], text: &[u8]) -> (i32, usize) {
    let m = read.len();
    if m == 0 {
        return (0, 0);
    }
    if m > LONG_MAX {
        return (m as i32, 0);
    }
    let rread: Vec<u8> = read.iter().rev().copied().collect();
    let rtext: Vec<u8> = text.iter().rev().copied().collect();
    let (score, rev_end) = infix_best_end(&rread, &rtext);
    (score, text.len() - rev_end)
}

/// Long-pattern (128 < m ≤ 256) 2-block chain; m ≤ 128 dispatches to the
/// single-block versions, m > 256 hits the guard.
pub mod long {
    use super::{calculate_block, ea_matches, mask_for};

    /// 2-block pattern table: positions i < 128 live in `peq0` bit i,
    /// positions i ≥ 128 in `peq1` bit i − 128.
    pub fn build_peq_2block(read: &[u8]) -> ([u128; 256], [u128; 256]) {
        let mut peq0 = [0u128; 256];
        let mut peq1 = [0u128; 256];
        for (i, &r) in read.iter().enumerate().take(256) {
            for f in b"ACGT" {
                if ea_matches(r, *f) {
                    if i < 128 {
                        peq0[*f as usize] |= 1u128 << i;
                    } else {
                        peq1[*f as usize] |= 1u128 << (i - 128);
                    }
                }
            }
        }
        for (lo, up) in [(b'a', b'A'), (b'c', b'C'), (b'g', b'G'), (b't', b'T')] {
            peq0[lo as usize] = peq0[up as usize];
            peq1[lo as usize] = peq1[up as usize];
        }
        (peq0, peq1)
    }

    fn scan2(
        peq0: &[u128; 256],
        peq1: &[u128; 256],
        m: usize,
        text: &[u8],
        mut f: impl FnMut(&mut i32, &mut usize, i32, usize),
    ) -> (i32, usize) {
        let m1 = m - 128;
        let mask1 = mask_for(m1);
        let high0 = 1u128 << 127;
        let high1 = 1u128 << (m1 - 1);
        let mut vp0 = u128::MAX;
        let mut vn0 = 0u128;
        let mut vp1 = mask1;
        let mut vn1 = 0u128;
        let mut score = m as i32;
        let mut best = score;
        let mut best_at = 0usize;
        for (j, &c) in text.iter().enumerate() {
            let (nvp0, nvn0, hout0) =
                calculate_block(vp0, vn0, peq0[c as usize], 0, u128::MAX, high0);
            let (nvp1, nvn1, hout1) =
                calculate_block(vp1, vn1, peq1[c as usize], hout0, mask1, high1);
            vp0 = nvp0;
            vn0 = nvn0;
            vp1 = nvp1;
            vn1 = nvn1;
            score += hout1;
            f(&mut best, &mut best_at, score, j + 1);
        }
        (best, best_at)
    }

    /// 2-block [`super::infix`](crate::myers_ea::infix).
    pub fn infix(read: &[u8], text: &[u8]) -> i32 {
        let m = read.len();
        if m <= 128 {
            return super::infix(read, text);
        }
        if m > 256 {
            return m as i32;
        }
        let (peq0, peq1) = build_peq_2block(read);
        scan2(&peq0, &peq1, m, text, |best, _at, score, _j| {
            if score < *best {
                *best = score;
            }
        })
        .0
    }

    /// 2-block [`super::infix_best_end`](crate::myers_ea::infix_best_end).
    pub fn infix_best_end(read: &[u8], text: &[u8]) -> (i32, usize) {
        let m = read.len();
        if m <= 128 {
            return super::infix_best_end(read, text);
        }
        if m > 256 {
            return (m as i32, 0);
        }
        let (peq0, peq1) = build_peq_2block(read);
        scan2(&peq0, &peq1, m, text, |best, at, score, j| {
            if score < *best {
                *best = score;
                *at = j;
            }
        })
    }

    /// 2-block [`super::infix_best_start`](crate::myers_ea::infix_best_start).
    pub fn infix_best_start(read: &[u8], text: &[u8]) -> (i32, usize) {
        let m = read.len();
        if m <= 128 {
            return super::infix_best_start(read, text);
        }
        if m > 256 {
            return (m as i32, 0);
        }
        let rread: Vec<u8> = read.iter().rev().copied().collect();
        let rtext: Vec<u8> = text.iter().rev().copied().collect();
        let (score, rev_end) = infix_best_end(&rread, &rtext);
        (score, text.len() - rev_end)
    }
}
