//! Sequence entry: tokenize (A7 C8 G9 T10 N11, others -> 6 [UNK], trailing SEP=1, no BOS) +
//! 1001bp (+/-500) reference window extraction (center anchor = 1-based pos, N-pad at contig edges).
//! Mirrors `_LUT`/`SEP_ID` in src/esperanto/model/encoder.py and `fetch_window` in
//! src/esperanto/features/sequence.py point by point.

use crate::bundle::ScoreError;
use rust_htslib::faidx;

pub use crate::caduceus::CaduceusEncoder;

pub const WINDOW: i64 = 500;
pub const TOTAL_BP: usize = (2 * WINDOW + 1) as usize; // 1001

/// Parameterized window (v1.2: half_window decided by the bundle feature_spec; v1.0 = 500).
pub fn total_bp_for(half: i64) -> usize {
    (2 * half + 1) as usize
}
pub const SEP_ID: i64 = 1;

#[inline]
fn lut(b: u8) -> i64 {
    match b {
        b'A' | b'a' => 7,
        b'C' | b'c' => 8,
        b'G' | b'g' => 9,
        b'T' | b't' => 10,
        b'N' | b'n' => 11,
        _ => 6, // [UNK]
    }
}

/// 1001bp window -> 1002 token ids (trailing SEP=1).
pub fn tokenize(window: &[u8]) -> Vec<i64> {
    let mut out = Vec::with_capacity(window.len() + 1);
    out.extend(window.iter().map(|&b| lut(b)));
    out.push(SEP_ID);
    out
}

/// 1001bp (2*WINDOW+1) reference window, center = pos (1-based); N-pad at contig edges and uppercase.
pub fn fetch_window(
    fasta: &mut faidx::Reader,
    chrom: &str,
    pos_1based: i64,
) -> Result<Vec<u8>, ScoreError> {
    fetch_window_hw(fasta, chrom, pos_1based, WINDOW)
}

/// Parameterized half-width version (v1.2: bundle half_window; 501bp = 250).
pub fn fetch_window_hw(
    fasta: &mut faidx::Reader,
    chrom: &str,
    pos_1based: i64,
    half: i64,
) -> Result<Vec<u8>, ScoreError> {
    let start = pos_1based - 1 - half; // 0-based inclusive
    let end_excl = pos_1based + half; // 0-based exclusive
    let chrom_len = fasta.fetch_seq_len(chrom) as i64;
    let c_start = start.max(0);
    let c_end = end_excl.min(chrom_len);

    let mut out = vec![b'N'; total_bp_for(half)];
    if c_end > c_start {
        // rust-htslib fetch_seq's end is inclusive
        let seq = fasta.fetch_seq(chrom, c_start as usize, (c_end - 1) as usize)?;
        let off = (c_start - start) as usize; // >0 when left-clipped
        for (i, &b) in seq.iter().enumerate() {
            out[off + i] = b.to_ascii_uppercase();
        }
    }
    Ok(out)
}

/// In-memory-reference version of fetch_window (semantics point-by-point identical to the faidx version: N-pad, uppercase, +/-500).
/// seq is the whole contig's uppercase sequence; chrom_len = seq.len().
pub fn fetch_window_mem(seq: &[u8], pos_1based: i64) -> Vec<u8> {
    fetch_window_mem_hw(seq, pos_1based, WINDOW)
}

/// Parameterized half-width version (v1.2).
pub fn fetch_window_mem_hw(seq: &[u8], pos_1based: i64, half: i64) -> Vec<u8> {
    let start = pos_1based - 1 - half; // 0-based inclusive (may be negative)
    let end_excl = pos_1based + half; // 0-based exclusive
    let chrom_len = seq.len() as i64;
    let c_start = start.max(0);
    let c_end = end_excl.min(chrom_len);
    let mut out = vec![b'N'; total_bp_for(half)];
    if c_end > c_start {
        let off = (c_start - start) as usize;
        let n = (c_end - c_start) as usize;
        out[off..off + n].copy_from_slice(&seq[c_start as usize..c_end as usize]);
    }
    out.to_vec()
}
