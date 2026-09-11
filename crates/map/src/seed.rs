//! Minimizer seeding (RNA high-density parameters).
//!
//! K-mers are encoded in 2-bit canonical form (min of forward and reverse
//! complement); window sort keys are the splitmix64 finalizer of the canonical
//! k-mer (not the raw k-mer value). K-mers overlapping an N (code 0xFF) are
//! invalid; a window with invalid members still selects among its valid ones.
//! Ties pick the leftmost k-mer; a minimizer shared by adjacent windows is
//! emitted only once.

use crate::fasta::Base;

/// Minimizer extraction parameters.
#[derive(Clone, Copy, Debug)]
pub struct SeedParams {
    /// K-mer length in bases.
    pub k: u32,
    /// Window size in k-mers.
    pub w: u32,
}

impl SeedParams {
    /// RNA defaults: k=15, w=5.
    pub fn rna_default() -> SeedParams {
        SeedParams { k: 15, w: 5 }
    }
}

/// Strand of a seed relative to its source sequence.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum Strand {
    /// Forward orientation.
    #[default]
    Plus,
    /// Reverse-complement orientation.
    Minus,
}

impl Strand {
    /// Opposite strand.
    pub fn flip(self) -> Strand {
        match self {
            Strand::Plus => Strand::Minus,
            Strand::Minus => Strand::Plus,
        }
    }
}

/// One minimizer: canonical 2-bit k-mer at a position.
#[derive(Clone, Copy, Debug)]
pub struct Minimizer {
    /// Canonical (min of k-mer and revcomp) 2-bit encoding.
    pub kmer: u64,
    /// 0-based position of the k-mer start.
    pub pos: u32,
    /// Orientation of the original k-mer (Plus if it equals the canonical).
    pub strand: Strand,
}

/// A seed hit pairing a read minimizer with a reference position.
#[derive(Clone, Copy, Debug)]
pub struct Anchor {
    /// Read position (minimizer start).
    pub qpos: u32,
    /// Reference position (minimizer start).
    pub rpos: u32,
    /// Contig index.
    pub contig: u32,
    /// Seed strand (read-relative XOR reference-relative).
    pub strand: Strand,
}

/// splitmix64 finalizer — window sort key for a canonical k-mer.
#[inline]
pub fn mix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E3779B97F4A7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
    x ^ (x >> 31)
}

const N_CODE: u64 = 0xFF;

/// Extract minimizers from an ASCII read sequence.
pub fn minimizers(seq: &[u8], params: SeedParams) -> Vec<Minimizer> {
    let codes: Vec<u64> = seq
        .iter()
        .map(|&b| match Base::from_ascii(b).code() {
            Some(c) => c as u64,
            None => N_CODE,
        })
        .collect();
    minimizers_from_codes(&codes, params.k, params.w)
}

/// Extract minimizers from pre-encoded 2-bit codes (`N` must be 0xFF).
///
/// Shared by read seeding and index construction. Uses a monotonic queue so
/// the whole pass is O(n).
pub fn minimizers_from_codes(codes: &[u64], k: u32, w: u32) -> Vec<Minimizer> {
    let k = k as usize;
    let w = w as usize;
    if codes.len() < k || w == 0 || k == 0 || k > 32 {
        return Vec::new();
    }
    let n_kmers = codes.len() - k + 1;
    // Per-k-mer (sort key, kmer, strand); k-mers overlapping an N are invalid.
    let mut keys: Vec<Option<(u64, u64, Strand)>> = Vec::with_capacity(n_kmers);
    let mut fwd: u64 = 0;
    let mut rev: u64 = 0;
    let mut valid_run = 0usize; // current N-free run length
    let mask: u64 = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    for (i, &c) in codes.iter().enumerate() {
        if c == N_CODE {
            valid_run = 0;
            fwd = 0;
            rev = 0;
        } else {
            fwd = ((fwd << 2) | c) & mask;
            // revcomp: new base complemented into the high slot
            rev = (rev >> 2) | ((3 - c) << (2 * (k - 1)));
            valid_run += 1;
        }
        if i + 1 >= k {
            if valid_run >= k {
                let (kmer, strand) = if fwd <= rev {
                    (fwd, Strand::Plus)
                } else {
                    (rev, Strand::Minus)
                };
                keys.push(Some((mix64(kmer), kmer, strand)));
            } else {
                keys.push(None);
            }
        }
    }
    debug_assert_eq!(keys.len(), n_kmers);

    // Sliding-window minimum (monotonic queue); a window containing invalid
    // k-mers still selects among its VALID members (only fully-invalid
    // windows emit nothing).
    let mut out: Vec<Minimizer> = Vec::new();
    let mut dq: std::collections::VecDeque<usize> = std::collections::VecDeque::new();
    for i in 0..n_kmers {
        if keys[i].is_some() {
            while let Some(&b) = dq.back() {
                let kb = keys[b].map(|x| x.0).unwrap_or(u64::MAX);
                let ki = keys[i].map(|x| x.0).unwrap_or(u64::MAX);
                // ties keep the leftmost (strict > pops)
                if kb > ki {
                    dq.pop_back();
                } else {
                    break;
                }
            }
            dq.push_back(i);
        }
        if i + 1 >= w {
            while let Some(&f) = dq.front() {
                if f + w <= i {
                    dq.pop_front();
                } else {
                    break;
                }
            }
            if let Some(&f) = dq.front() {
                if let Some((_, kmer, strand)) = keys[f] {
                    // Adjacent windows selecting the same minimizer emit once.
                    if out.last().map(|m: &Minimizer| (m.kmer, m.pos)) != Some((kmer, f as u32)) {
                        out.push(Minimizer {
                            kmer,
                            pos: f as u32,
                            strand,
                        });
                    }
                }
            }
        }
    }
    out
}

