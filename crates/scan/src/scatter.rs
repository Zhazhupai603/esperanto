//! Block pileup: scatter each read's counted bases into block-local arrays.
//!
//! Per position: total depth (all counted bases) and, for positions where at
//! least one counted base differs from the reference, the per-base mismatch
//! counts. The reference count is `depth - mismatches`. Reads on the reverse
//! strand are stored reference-forward in the BAM, so bases are compared to
//! the reference as stored, with no orientation split.

use crate::error::CallError;
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::HashMap;
use std::path::Path;

pub fn base_idx(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Block-local pileup result over `[cs, ce)`.
pub struct BlockAcc {
    pub cs: u64,
    pub depth: Vec<u32>,
    /// genomic pos0 -> mismatch counts by base (A, C, G, T)
    pub mism: HashMap<u32, [u32; 4]>,
}

impl BlockAcc {
    fn new(cs: u64, ce: u64) -> Self {
        BlockAcc {
            cs,
            depth: vec![0u32; (ce - cs) as usize],
            mism: HashMap::new(),
        }
    }
}

/// Read-end trim: bases within this distance of either end of the read are
/// not counted (Ramaswami et al., Nat Methods 9:490, 2012). Alignment
/// artifacts concentrate at read ends; the trim is symmetric, so the BAM's
/// reference-forward storage needs no orientation handling.
pub const READ_END_TRIM: usize = 6;

/// Add one read to the block accumulator.
#[allow(clippy::too_many_arguments)]
pub fn scatter_read(
    pos: i64,
    cigar: &[(char, i64)],
    seq: &[u8],
    qual: &[u8],
    min_bq: u8,
    cs: u64,
    ce: u64,
    refseq: &[u8],
    acc: &mut BlockAcc,
) {
    let mut rp = pos;
    let mut qp = 0usize;
    for &(op, l) in cigar {
        match op {
            'M' | '=' | 'X' => {
                let os = rp.max(cs as i64);
                let oe = (rp + l).min(ce as i64);
                if os < oe {
                    let q0 = (qp as i64 + (os - rp)) as usize;
                    for k in 0..(oe - os) as usize {
                        let qi = q0 + k;
                        if qi >= seq.len() {
                            break;
                        }
                        // Read-end trim: ends of the read carry the bulk of
                        // the alignment-driven artifacts.
                        if qi < READ_END_TRIM || qi + READ_END_TRIM >= seq.len() {
                            continue;
                        }
                        if qual.get(qi).copied().unwrap_or(0) < min_bq {
                            continue;
                        }
                        let Some(bi) = base_idx(seq[qi]) else { continue };
                        let gpos = os as usize + k;
                        let Some(&rb) = refseq.get(gpos) else { continue };
                        if base_idx(rb).is_none() {
                            continue;
                        }
                        acc.depth[gpos - cs as usize] += 1;
                        if seq[qi] != rb {
                            acc.mism.entry(gpos as u32).or_default()[bi] += 1;
                        }
                    }
                }
                rp += l;
                qp += l as usize;
            }
            'I' | 'S' => qp += l as usize,
            'D' | 'N' => rp += l,
            _ => {}
        }
    }
}

/// Pileup of one contig block `[cs, ce)` from an indexed BAM.
pub fn scatter_block(
    bam_path: &Path,
    chrom: &str,
    cs: u64,
    ce: u64,
    refseq: &[u8],
    min_mapq: u8,
    min_bq: u8,
) -> Result<BlockAcc, CallError> {
    let ioe = |e: rust_htslib::errors::Error| CallError::Io {
        path: bam_path.display().to_string(),
        source: std::io::Error::other(e.to_string()),
    };
    let mut reader = bam::IndexedReader::from_path(bam_path).map_err(ioe)?;
    let tid = reader
        .header()
        .tid(chrom.as_bytes())
        .ok_or_else(|| CallError::ContigNotFound(chrom.to_string()))?;
    reader.fetch((tid, cs, ce)).map_err(ioe)?;

    let mut acc = BlockAcc::new(cs, ce);
    for r in reader.records() {
        let rec = r.map_err(ioe)?;
        if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() {
            continue;
        }
        if rec.mapq() < min_mapq {
            continue;
        }
        let pos = rec.pos();
        if pos < 0 || pos as u64 >= ce {
            continue;
        }
        let cigar = rec.cigar();
        let ops: Vec<(char, i64)> = cigar.iter().map(|c| (c.char(), c.len() as i64)).collect();
        let seq: Vec<u8> = rec.seq().as_bytes();
        scatter_read(pos, &ops, &seq, rec.qual(), min_bq, cs, ce, refseq, &mut acc);
    }
    Ok(acc)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_depth_and_mismatches_with_quality_gate() {
        let refseq = b"ACGTACGTACGTACGTACGT".to_vec();
        let mut acc = BlockAcc::new(0, 20);
        // 24-base read placed at ref 2..26 would overflow; use ref 0..20
        // padding plus a read long enough that the end trim leaves a middle.
        // Read covers ref 0..20; query base 8 (ref pos 8, A-slot) reads G;
        // query base 9 is below the quality gate.
        let mut seq = b"GTGCGT".to_vec();
        seq.extend_from_slice(b"TTTTTTTTTTTTTTTTTT");
        let mut qual = [40u8; 24];
        qual[9] = 10;
        scatter_read(0, &[('M', 20)], &seq, &qual, 30, 0, 20, &refseq, &mut acc);
        // ends trimmed: only query 6..18 counted (ref 6..18)
        assert_eq!(acc.depth[6], 1);
        assert_eq!(acc.depth[18], 0, "beyond trim");
        assert_eq!(acc.depth[8], 1);
        assert_eq!(acc.depth[9], 0, "low-quality base not counted");
        assert_eq!(acc.mism.get(&8), Some(&[0, 0, 0, 1]));
    }

    #[test]
    fn read_end_trim_excludes_six_bases_each_side() {
        let refseq = vec![b'A'; 40];
        let mut acc = BlockAcc::new(0, 40);
        let seq = vec![b'C'; 30]; // all mismatches vs A reference
        let qual = [40u8; 30];
        scatter_read(0, &[('M', 30)], &seq, &qual, 30, 0, 40, &refseq, &mut acc);
        // counted: query 6..24 -> ref 6..24
        for i in 0..6 {
            assert_eq!(acc.depth[i], 0, "head trimmed at {i}");
            assert_eq!(acc.depth[29 - i], 0, "tail trimmed at {}", 29 - i);
        }
        for i in 6..24 {
            assert_eq!(acc.depth[i], 1, "middle counted at {i}");
        }
    }

    #[test]
    fn skips_refskip_and_softclip() {
        let refseq = b"AAAAAAAAAAAAAAAAAAAA".to_vec();
        let mut acc = BlockAcc::new(0, 20);
        // 2S 3M 5N 3M: query bases 2..5 -> ref 0..3, query 5..8 -> ref 8..11
        let seq = b"GGGGGGGG";
        let qual = [40u8; 8];
        scatter_read(
            0,
            &[('S', 2), ('M', 3), ('N', 5), ('M', 3)],
            seq,
            &qual,
            30,
            0,
            20,
            &refseq,
            &mut acc,
        );
        // query 2..8 after a 2-base clip: query 2..8 are read positions
        // 2..8, all inside the 6-base end trim of an 8-base read — nothing
        // counts. (Ends of short reads are uninformative by design.)
        assert_eq!(&acc.depth[0..3], &[0, 0, 0]);
        assert_eq!(&acc.depth[8..11], &[0, 0, 0]);
        assert_eq!(acc.mism.len(), 0);
    }
}
