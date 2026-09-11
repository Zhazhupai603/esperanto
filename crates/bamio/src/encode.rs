//! Direct BAM record binary encoder — zero per-record allocation.
//!
//! Serializes [`OutRecord`] straight into a reusable `Vec<u8>`, bypassing the
//! noodles `RecordBuf` builder. The output is byte-identical to
//! `noodles::bam::io::Writer::write_alignment_record` for the record shapes
//! produced here. Unsupported shapes fall back to the caller's RecordBuf path.

use std::io;

use crate::{AlnView, CigarOp, OutRecord, TagValue};

// § 4.2.3 "SEQ and QUAL encoding" — case-insensitive base → 4-bit code table.
const SEQ_CODES: [u8; 256] = build_seq_codes();

const fn build_seq_codes() -> [u8; 256] {
    const BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
    const N: u8 = 0x0f;
    let mut codes = [N; 256];
    let mut i = 0;
    while i < BASES.len() {
        codes[BASES[i] as usize] = i as u8;
        codes[BASES[i].to_ascii_lowercase() as usize] = i as u8;
        i += 1;
    }
    codes
}

// § 5.3 reg2bin — start/end are 0-based, inclusive.
const UNMAPPED_BIN: u16 = 4680;

#[allow(clippy::eq_op)] // § 5.3 reg2bin — constant formula per BAM spec
fn reg2bin(start: usize, end: usize) -> u16 {
    let bin = if start >> 14 == end >> 14 {
        ((1 << 15) - 1) / 7 + (start >> 14)
    } else if start >> 17 == end >> 17 {
        ((1 << 12) - 1) / 7 + (start >> 17)
    } else if start >> 20 == end >> 20 {
        ((1 << 9) - 1) / 7 + (start >> 20)
    } else if start >> 23 == end >> 23 {
        ((1 << 6) - 1) / 7 + (start >> 23)
    } else if start >> 26 == end >> 26 {
        ((1 << 3) - 1) / 7 + (start >> 26)
    } else {
        0
    };
    bin as u16
}

/// Attempts to encode `rec` directly into `buf`.
///
/// Returns `Some(Ok(()))` on success, `Some(Err(_))` on a recoverable encoding
/// error (e.g. seq/qual mismatch), or `None` if the record shape is not
/// supported by this fast path (caller falls back to noodles RecordBuf).
pub fn try_encode(buf: &mut Vec<u8>, rec: &OutRecord) -> Option<io::Result<()>> {
    // Name validity (mirrors noodles name::is_valid).
    let name = rec.name.as_bytes();
    if !(1..=254).contains(&name.len()) || name == b"*" {
        return None;
    }
    if !name.iter().all(|&b| b.is_ascii_graphic() && b != b'@') {
        return None;
    }
    // CIGAR oversized-cg handling (>65535 ops) not implemented in fast path.
    let n_cigar = rec.aln.as_ref().map_or(0, |a| a.cigar.len());
    if n_cigar > 0xffff {
        return None;
    }
    // Quality scores must be raw Phred [0, 93] (§ 4.2.3).
    if rec.qual.iter().any(|&q| q > 93) {
        return None;
    }
    // Contig ID must fit in i32.
    if rec.aln.as_ref().is_some_and(|a| a.contig > i32::MAX as u32) {
        return None;
    }
    Some(encode(buf, rec, n_cigar))
}

fn encode(buf: &mut Vec<u8>, rec: &OutRecord, n_cigar: usize) -> io::Result<()> {
    let (ref_id, pos, bin): (i32, i32, u16) = match &rec.aln {
        Some(aln) => {
            let p = aln.pos as usize;
            let span: usize = aln
                .cigar
                .iter()
                .map(|op| match op {
                    CigarOp::Match(n) | CigarOp::Del(n) | CigarOp::RefSkip(n) => *n as usize,
                    CigarOp::Ins(_) | CigarOp::SoftClip(_) => 0,
                })
                .sum();
            let bin = if span > 0 {
                reg2bin(p, p + span - 1)
            } else {
                UNMAPPED_BIN
            };
            (aln.contig as i32, p as i32, bin)
        }
        None => (-1, -1, UNMAPPED_BIN),
    };

    let (next_ref, next_pos, tlen): (i32, i32, i32) = match rec.mate {
        Some((mc, mp, t)) if mc >= 0 => (mc, mp, t),
        _ => (-1, -1, 0),
    };

    // --- fixed 32-byte alignment header ---
    buf.extend_from_slice(&ref_id.to_le_bytes());
    buf.extend_from_slice(&pos.to_le_bytes());
    buf.push(rec.name.len() as u8 + 1); // l_read_name = name + NUL (guarded by try_encode)
    buf.push(if rec.aln.is_some() { rec.mapq } else { 255 });
    buf.extend_from_slice(&bin.to_le_bytes());
    buf.extend_from_slice(&(n_cigar as u16).to_le_bytes());
    buf.extend_from_slice(&rec.flag.to_le_bytes());
    buf.extend_from_slice(&(rec.seq.len() as u32).to_le_bytes());
    buf.extend_from_slice(&next_ref.to_le_bytes());
    buf.extend_from_slice(&next_pos.to_le_bytes());
    buf.extend_from_slice(&tlen.to_le_bytes());

    // --- read_name (NUL-terminated) ---
    buf.extend_from_slice(rec.name.as_bytes());
    buf.push(0);

    // --- cigar ops (u32 LE, kind | len<<4) ---
    if let Some(aln) = &rec.aln {
        for op in &aln.cigar {
            let (kind, n) = cigar_code(op);
            buf.extend_from_slice(&((n << 4) | kind).to_le_bytes());
        }
    }

    // --- seq (4-bit packed, high nibble first; odd last low nibble = 0) ---
    let seq = &rec.seq;
    let mut i = 0;
    while i < seq.len() {
        let hi = SEQ_CODES[seq[i] as usize];
        let lo = if i + 1 < seq.len() {
            SEQ_CODES[seq[i + 1] as usize]
        } else {
            0
        };
        buf.push((hi << 4) | lo);
        i += 2;
    }

    // --- qual ---
    let l_seq = seq.len();
    if rec.qual.len() == l_seq {
        buf.extend_from_slice(&rec.qual);
    } else if rec.qual.is_empty() {
        buf.resize(buf.len() + l_seq, 0xff);
    } else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "seq/qual length mismatch",
        ));
    }

    // --- tags (written verbatim in the order given by RawTag) ---
    if let Some(aln) = &rec.aln {
        write_tags(buf, aln);
    }

    Ok(())
}

#[inline]
fn cigar_code(op: &CigarOp) -> (u32, u32) {
    match op {
        CigarOp::Match(n) => (0, *n),
        CigarOp::Ins(n) => (1, *n),
        CigarOp::Del(n) => (2, *n),
        CigarOp::RefSkip(n) => (3, *n),
        CigarOp::SoftClip(n) => (4, *n),
    }
}

fn write_tags(buf: &mut Vec<u8>, aln: &AlnView) {
    for crate::RawTag(tag, value) in &aln.tags {
        match value {
            TagValue::Char(c) => {
                buf.extend_from_slice(tag);
                buf.push(b'A');
                buf.push(*c);
            }
            TagValue::Int(i) => {
                buf.extend_from_slice(tag);
                buf.push(b'i');
                buf.extend_from_slice(&i.to_le_bytes());
            }
            TagValue::Str(s) => {
                buf.extend_from_slice(tag);
                buf.push(b'Z');
                buf.extend_from_slice(s.as_bytes());
                buf.push(0);
            }
        }
    }
}
