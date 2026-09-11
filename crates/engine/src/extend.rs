//! EA seed extension: greedy bidirectional walk under [`crate::ea_free`].

/// Result of one bidirectional EA extension from a seed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Extension {
    /// Query interval `[read_lo, read_hi)` in oriented-read coordinates.
    pub read_lo: usize,
    pub read_hi: usize,
    /// Transcript interval `[tx_lo, tx_hi)`.
    pub tx_lo: usize,
    pub tx_hi: usize,
    /// True when the extension spans the whole read.
    pub full: bool,
}

impl Extension {
    /// Number of read bases covered by the extension.
    pub fn read_cov(&self) -> usize {
        self.read_hi - self.read_lo
    }
}

/// Greedy bidirectional extension of the seed
/// `read[a..a + k] == tx_seq[t..t + k]` (EA-equal by index construction).
///
/// Both sides advance one base at a time while each new pair satisfies
/// the EA predicate; the first true mismatch stops that direction. The
/// seed window itself is always accepted.
pub fn extend_ea(read: &[u8], a: usize, k: usize, tx_seq: &[u8], t: usize) -> Extension {
    let read_len = read.len();
    let mut lo_r = a;
    let mut lo_t = t;
    while lo_r > 0 && lo_t > 0 && crate::ea_free(tx_seq[lo_t - 1], read[lo_r - 1]) {
        lo_r -= 1;
        lo_t -= 1;
    }
    let mut hi_r = a + k;
    let mut hi_t = t + k;
    while hi_r < read_len && hi_t < tx_seq.len() && crate::ea_free(tx_seq[hi_t], read[hi_r]) {
        hi_r += 1;
        hi_t += 1;
    }
    Extension {
        read_lo: lo_r,
        read_hi: hi_r,
        tx_lo: lo_t,
        tx_hi: hi_t,
        full: hi_r - lo_r == read_len,
    }
}
