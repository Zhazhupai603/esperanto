//! Scatter-count pileup engine (v14 call_v2) — replaces htslib bam_plp.
//!
//! Math: instead of a per-column state machine (~100ns/entry), scatter each
//! read's contribution directly into block-local arrays (~2-3ns/base):
//!   - depth/bq_sum/mapq_sum scattered per base (no state machine);
//!   - junction strand votes scattered per base (per-read intron vote counts);
//!   - mismatch localization only for dirty reads (EK mm+ea>0, ~26%).
//!
//! Semantics mirror the legacy htslib path (CALL_V2_DESIGN.md S1-S12):
//! per-site counting includes an entry iff bq>=13 AND base in ACGT;
//! D/N CIGAR blocks occupy no column; junction votes gated the same way.

use crate::error::CallError;
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::HashMap;
use std::path::Path;

/// Same threshold as the legacy column processor (count.rs MIN_BASE_QUALITY).
pub const MIN_BQ: u8 = 13;

fn base_idx(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Parse EK aux tag for (mm, ea). Absent → None (caller treats read as dirty).
fn ek_counts(rec: &bam::record::Record) -> Option<(u32, u32)> {
    let aux = rec.aux(b"EK").ok()?;
    let s = match aux {
        bam::record::Aux::String(v) => v,
        _ => return None,
    };
    ek_str_counts(s)
}

/// EK:Z string → (mm, ea) (same format for BAM aux and .baln aux).
fn ek_str_counts(s: &str) -> Option<(u32, u32)> {
    let mut mm = None;
    let mut ea = None;
    for part in s.split(';') {
        if let Some(v) = part.strip_prefix("mm=") {
            mm = v.parse::<u32>().ok();
        } else if let Some(v) = part.strip_prefix("ea=") {
            ea = v.parse::<u32>().ok();
        }
    }
    Some((mm?, ea?))
}

/// .baln CIGAR u32 ((len<<4)|op) → (char, len).
fn baln_cigar_op(c: u32) -> Option<(char, i64)> {
    let ch = match c & 0xF {
        0 => 'M',
        1 => 'I',
        2 => 'D',
        3 => 'N',
        4 => 'S',
        5 => 'H',
        6 => 'P',
        7 => '=',
        8 => 'X',
        _ => return None,
    };
    Some((ch, (c >> 4) as i64))
}

/// Per-site tallies for sites with ≥1 mismatch: 4 bases × 2 strands.
#[derive(Default, Clone, Copy)]
pub struct SiteAcc {
    pub cnt: [u32; 8], // A,C,G,T on fwd then rev
}

/// Per-block scatter result. All arrays are block-local ([0, ce-cs)).
pub struct BlockAcc {
    pub depth_fwd: Vec<u32>,
    pub depth_rev: Vec<u32>,
    pub bq_sum: Vec<u64>,
    pub mapq_sum: Vec<u64>,
    pub junc_plus: Vec<u32>,
    pub junc_minus: Vec<u32>,
    pub sites: HashMap<u32, SiteAcc>, // key = genomic pos (0-based)
    pub jbounds: Vec<i32>,
}

impl BlockAcc {
    fn new(bn: usize) -> Self {
        BlockAcc {
            depth_fwd: vec![0u32; bn],
            depth_rev: vec![0u32; bn],
            bq_sum: vec![0u64; bn],
            mapq_sum: vec![0u64; bn],
            junc_plus: vec![0u32; bn],
            junc_minus: vec![0u32; bn],
            sites: HashMap::new(),
            jbounds: Vec::new(),
        }
    }
}

/// Scatter kernel for a single record (BAM/baln dual sources share the same code path — structural parity).
/// Input fields use as-sequenced convention; RC(seq)+rev(qual) for minus-strand reads is done here
/// (htslib pileup forward-strand view). cigar is an (op, len) sequence.
#[allow(clippy::too_many_arguments)]
fn scatter_one_record(
    pos: i64,
    is_rev: bool,
    mapq: u64,
    dirty: bool,
    cigar: &[(char, i64)],
    seq_raw: &[u8],
    qual_raw: &[u8],
    cs: u64,
    ce: u64,
    refseq: &[u8],
    has_fasta: bool,
    acc: &mut BlockAcc,
) {
    // v14.2 orientation final review (decided by per-site control measurements): per the SAM spec, BAM SEQ for
    // minus-strand reads is already reference-forward (SAM spec "segments are represented on the forward genomic
    // strand"; STAR output verified by direct comparison against the reference via samtools view), so the htslib
    // pileup reads it forward directly. v14's "RC fix" rested on the wrong premise that "SEQ is as-sequenced":
    // after RC, minus-strand read bases were mirror-shifted against reference positions and complement-flipped →
    // measured minus-strand rev_freq stuck at ~1.0 (every minus-strand read miscounted as a mismatch) and var_freq
    // inflated (0.75 vs true 0.5). The correct semantics = pass-through (consistent with htslib pileup / esperanto-pile).
    let seqb: &[u8] = seq_raw;
    let qual_v: &[u8] = qual_raw;
    let seq_len = seqb.len();

    // Walk 1: intron strand votes for this read (plus/minus counts) + jbounds.
    let (mut n_plus, mut n_minus) = (0u32, 0u32);
    {
        let mut rp = pos;
        for &(op, l) in cigar {
            match op {
                'M' | '=' | 'X' | 'D' => rp += l,
                'N' => {
                    let is = rp;
                    let ie = rp + l - 1;
                    let vote = if has_fasta {
                        let d0 = refseq.get(is as usize).copied().unwrap_or(b'N');
                        let d1 = refseq.get(is as usize + 1).copied().unwrap_or(b'N');
                        let a0 = refseq
                            .get((ie - 1).max(0) as usize)
                            .copied()
                            .unwrap_or(b'N');
                        let a1 = refseq.get(ie as usize).copied().unwrap_or(b'N');
                        match ([d0, d1], [a0, a1]) {
                            // Byte-for-byte aligned with legacy intron_strand: GT-AG/GC-AG/AT-AC → plus strand.
                            ([b'G', b'T'], [b'A', b'G'])
                            | ([b'G', b'C'], [b'A', b'G'])
                            | ([b'A', b'T'], [b'A', b'C']) => Some(true),
                            ([b'C', b'T'], [b'A', b'C']) => Some(false),
                            _ => None,
                        }
                    } else {
                        None
                    };
                    acc.jbounds.push(is as i32);
                    acc.jbounds.push(ie as i32);
                    match vote {
                        Some(true) => n_plus += 1,
                        Some(false) => n_minus += 1,
                        None => {}
                    }
                    rp += l;
                }
                _ => {}
            }
        }
    }

    // Walk 2: per-base scatter over M/=/X positions within [cs, ce).
    let mut rp = pos;
    let mut qp = 0usize;
    for &(op, l) in cigar {
        match op {
            'M' | '=' | 'X' => {
                // overlap with block
                let os = rp.max(cs as i64);
                let oe = (rp + l).min(ce as i64);
                if os < oe {
                    let qoff = (qp as i64 + (os - rp)) as usize;
                    let qend = (qp as i64 + (oe - rp)) as usize;
                    let base0 = (os - cs as i64) as usize;
                    let depth_vec = if is_rev {
                        &mut acc.depth_rev
                    } else {
                        &mut acc.depth_fwd
                    };
                    for (k, qi) in (qoff..qend).enumerate() {
                        let p = base0 + k;
                        // bq: qual[qpos] if within seq_len else 0 (legacy parity);
                        // missing qual array → 0xFF (kept), mirroring legacy.
                        let bq = if qi < seq_len {
                            qual_v.get(qi).copied().unwrap_or(0xFF)
                        } else {
                            0
                        };
                        if bq < MIN_BQ {
                            continue;
                        }
                        let Some(bi) = base_idx(seqb[qi]) else {
                            continue;
                        };
                        depth_vec[p] += 1;
                        acc.bq_sum[p] += u64::from(bq);
                        acc.mapq_sum[p] += mapq;
                        if n_plus > 0 {
                            acc.junc_plus[p] += n_plus;
                        }
                        if n_minus > 0 {
                            acc.junc_minus[p] += n_minus;
                        }
                        if dirty {
                            let rb = if has_fasta {
                                refseq.get(os as usize + k).copied().unwrap_or(b'N')
                            } else {
                                b'N'
                            };
                            let e = acc.sites.entry((os + k as i64) as u32).or_default();
                            let slot = if is_rev { bi + 4 } else { bi };
                            if rb == b'A' || rb == b'C' || rb == b'G' || rb == b'T' {
                                // With reference: count mismatches only
                                if seqb[qi] != rb {
                                    e.cnt[slot] += 1;
                                }
                            } else {
                                // N site / no fasta: count everything (majority pseudo-ref path)
                                e.cnt[slot] += 1;
                            }
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

/// Scatter one contig block [cs, ce) — BAM source (htslib IndexedReader + fetch).
/// refseq = whole contig sequence. `has_fasta` gates junction voting exactly
/// like the legacy path (votes only when contig_seq is Some).
pub fn scatter_block(
    bam_path: &Path,
    chrom: &str,
    cs: u64,
    ce: u64,
    refseq: &[u8],
    has_fasta: bool,
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

    let mut acc = BlockAcc::new((ce - cs) as usize);
    for r in reader.records() {
        let rec = r.map_err(|e| CallError::Io {
            path: bam_path.display().to_string(),
            source: std::io::Error::other(e.to_string()),
        })?;
        if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() {
            continue;
        }
        let pos = rec.pos();
        if pos < 0 || pos as u64 >= ce {
            continue;
        }
        let is_rev = rec.is_reverse();
        let mapq = u64::from(rec.mapq());
        let dirty = !has_fasta || ek_counts(&rec).is_none_or(|(mm, ea)| mm + ea > 0);
        let cigar = rec.cigar();
        let ops: Vec<(char, i64)> = cigar.iter().map(|c| (c.char(), c.len() as i64)).collect();
        // rec.seq() returns a borrowed view; as_bytes() yields a temporary Vec — bind first, then pass the reference.
        let seq_bytes: Vec<u8> = rec.seq().as_bytes();
        scatter_one_record(
            pos,
            is_rev,
            mapq,
            dirty,
            &ops,
            &seq_bytes,
            rec.qual(),
            cs,
            ce,
            refseq,
            has_fasta,
            &mut acc,
        );
    }
    Ok(acc)
}

/// Scatter one contig block [cs, ce) — .baln source (v14 binary input).
/// idx/max_span come from `BalnReader::build_index` (shared read-only across block tasks); each task has
/// its own file handle + seek. The block window uses overlap semantics (pos ∈ [cs−max_span, ce) and
/// pos+span > cs), aligned with htslib fetch overlap detection; scattering goes through the same
/// scatter_one_record, structurally equivalent to the BAM path.
#[allow(clippy::too_many_arguments)]
pub fn scatter_block_baln(
    baln_path: &Path,
    idx: &[(i32, i64, u64, i64)],
    max_span: i64,
    tid: i32,
    cs: u64,
    ce: u64,
    refseq: &[u8],
    has_fasta: bool,
) -> Result<BlockAcc, CallError> {
    let mut acc = BlockAcc::new((ce - cs) as usize);
    let lo_pos = cs as i64 - max_span;
    let lo = idx.partition_point(|e| (e.0, e.1) < (tid, lo_pos));
    let hi = idx.partition_point(|e| (e.0, e.1) < (tid, ce as i64));
    let file = std::fs::File::open(baln_path).map_err(|e| CallError::Io {
        path: baln_path.display().to_string(),
        source: e,
    })?;
    for &(_t, pos, off, span) in &idx[lo..hi] {
        if pos + span <= cs as i64 {
            continue; // inside the index window but does not overlap the block (short read entirely outside)
        }
        let Some(rec) = crate::baln::read_record_at(&file, off)? else {
            continue;
        };
        // Same filters as the BAM path (unmapped already excluded at index time; belt-and-braces + secondary/supplementary)
        if rec.flag & 0x4 != 0 || rec.flag & 0x100 != 0 || rec.flag & 0x800 != 0 {
            continue;
        }
        if rec.pos < 0 || rec.pos as u64 >= ce {
            continue;
        }
        let dirty = !has_fasta
            || rec
                .ek
                .as_deref()
                .and_then(ek_str_counts)
                .is_none_or(|(mm, ea)| mm + ea > 0);
        let ops: Vec<(char, i64)> = rec.cigar.iter().filter_map(|&c| baln_cigar_op(c)).collect();
        scatter_one_record(
            rec.pos,
            rec.flag & 0x10 != 0,
            u64::from(rec.mapq),
            dirty,
            &ops,
            &rec.seq_ascii,
            &rec.qual,
            cs,
            ce,
            refseq,
            has_fasta,
            &mut acc,
        );
    }
    Ok(acc)
}
