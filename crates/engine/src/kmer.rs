//! K-mer encoding and the N-skip reseed stream.
//!
//! 2-bit codes: `A=0 C=1 G=2 T=3`, case-insensitive, with the 5'
//! base in the highest bits of the window. The canonical form of a k-mer
//! is `min(forward, revcomp)`. Any non-ACGT byte (N / IUPAC / other)
//! invalidates every window overlapping it; after an invalid byte the
//! stream reseeds from the byte after it (no phantom codes).

/// 2-bit code of one base byte, case-insensitive; `None` for any non-ACGT.
#[inline]
pub fn base_code(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Reverse complement of a packed 2-bit k-mer code of length `k`
/// (complement every code and reverse their order).
pub fn revcomp_code(code: u64, k: usize) -> u64 {
    let mut rc = 0u64;
    for i in 0..k {
        let c = (code >> (2 * i)) & 3;
        rc = (rc << 2) | (3 - c);
    }
    rc
}

/// Canonical form of a forward k-mer code: `min(forward, revcomp)`.
pub fn canonical(code: u64, k: usize) -> u64 {
    code.min(revcomp_code(code, k))
}

/// Forward 2-bit code of `seq[off..off + k]`; `None` when the window is
/// out of bounds or contains a non-ACGT byte.
pub fn window_code(seq: &[u8], off: usize, k: usize) -> Option<u64> {
    if off + k > seq.len() {
        return None;
    }
    let mut code = 0u64;
    for &b in &seq[off..off + k] {
        code = (code << 2) | base_code(b)?;
    }
    Some(code)
}

/// Case-preserving reverse complement: `a<->t`, `c<->g`; non-ACGT bytes
/// pass through unchanged. Used for minus-strand oriented reads — k-mer
/// encoding is case-insensitive, this is consistency only.
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'a' => b't',
            b'C' => b'G',
            b'c' => b'g',
            b'G' => b'C',
            b'g' => b'c',
            b'T' => b'A',
            b't' => b'a',
            other => other,
        })
        .collect()
}

/// Streaming forward k-mer codes over a read with N-skip reseed.
///
/// Yields `(offset, forward_code)` for every ACGT-only window in ascending
/// offset order. When a window contains a non-ACGT byte, the stream
/// reseeds at the byte after the first invalid byte in that window.
#[derive(Debug, Clone)]
pub struct KmerStream<'a> {
    seq: &'a [u8],
    k: usize,
    pos: usize,
}

impl<'a> KmerStream<'a> {
    /// Create the stream; `k` must be in `1..=32`.
    pub fn new(seq: &'a [u8], k: usize) -> Self {
        assert!((1..=32).contains(&k), "k-mer size must be in 1..=32");
        KmerStream { seq, k, pos: 0 }
    }
}

impl<'a> Iterator for KmerStream<'a> {
    type Item = (usize, u64);

    fn next(&mut self) -> Option<(usize, u64)> {
        loop {
            if self.pos + self.k > self.seq.len() {
                return None;
            }
            match window_code(self.seq, self.pos, self.k) {
                Some(code) => {
                    let p = self.pos;
                    self.pos += 1;
                    return Some((p, code));
                }
                None => {
                    // reseed: skip to the byte after the first invalid one
                    let stop = self.pos + self.k;
                    let mut next = stop;
                    for (i, &b) in self.seq[self.pos..stop].iter().enumerate() {
                        if base_code(b).is_none() {
                            next = self.pos + i + 1;
                            break;
                        }
                    }
                    self.pos = next;
                }
            }
        }
    }
}
