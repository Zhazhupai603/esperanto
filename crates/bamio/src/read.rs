//! Read side: sequential BAM record view + original-orientation restore.
//!
//! Manual parsing (SAM §4.2 record layout over the BGZF decompression stream); strictly inverse to encode.rs.
//! No index/region fetch — region consumers (pile/scan) are pinned to rust-htslib semantics and do not use this module.

use std::fs::File;
use std::io::{self, BufReader, Read};
use std::path::Path;

use crate::{CigarOp, RawTag, TagValue};

/// A record read back (SEQ is in stored orientation = reference-forward).
#[derive(Clone, Debug)]
pub struct InRecord {
    pub name: String,
    pub flag: u16,
    pub mapq: u8,
    /// Header SQ order; -1 = no reference.
    pub ref_id: i32,
    /// 0-based; -1 = no reference.
    pub pos: i64,
    pub cigar: Vec<CigarOp>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    /// File order (only A/i/Z types are parsed; other types are skipped, not entered into the table).
    pub tags: Vec<RawTag>,
}

/// Original-orientation restore: 0x10 set → revcomp + QUAL reversed back; otherwise unchanged.
pub fn restore_original(flag: u16, seq: &[u8], qual: &[u8]) -> (Vec<u8>, Vec<u8>) {
    crate::apply_t13(flag & crate::flag::REVERSE != 0, seq, qual)
}

/// BAM read-side error.
#[derive(Debug, thiserror::Error)]
pub enum ReadError {
    #[error("io: {0}")]
    Io(#[from] io::Error),
    #[error("bad magic: not BAM")]
    BadMagic,
    #[error("truncated record stream")]
    Truncated,
    #[error("invalid CIGAR op kind {0}")]
    BadCigarKind(u32),
}

const SEQ_DECODE: [u8; 16] = *b"=ACMGRSVTWYHKDBN";

struct Parser {
    buf: Vec<u8>,
    off: usize,
}

impl Parser {
    fn take(&mut self, n: usize) -> Result<&[u8], ReadError> {
        if self.off + n > self.buf.len() {
            return Err(ReadError::Truncated);
        }
        let s = &self.buf[self.off..self.off + n];
        self.off += n;
        Ok(s)
    }
    fn u8(&mut self) -> Result<u8, ReadError> {
        Ok(self.take(1)?[0])
    }
    fn u16(&mut self) -> Result<u16, ReadError> {
        let b = self.take(2)?;
        Ok(u16::from_le_bytes([b[0], b[1]]))
    }
    fn i32(&mut self) -> Result<i32, ReadError> {
        let b = self.take(4)?;
        Ok(i32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }
    fn u32(&mut self) -> Result<u32, ReadError> {
        Ok(self.i32()? as u32)
    }
}

fn read_exact_stream(r: &mut impl Read, n: usize) -> Result<Vec<u8>, ReadError> {
    let mut buf = vec![0u8; n];
    let mut filled = 0;
    while filled < n {
        let got = r.read(&mut buf[filled..])?;
        if got == 0 {
            return Err(ReadError::Truncated);
        }
        filled += got;
    }
    Ok(buf)
}

/// Open a BAM for whole-file sequential iteration.
///
/// Returns (contig name table, record iterator). The iterator yields `Result<InRecord>` items.
pub fn open_sequential(
    path: &Path,
) -> Result<
    (
        Vec<String>,
        impl Iterator<Item = Result<InRecord, ReadError>>,
    ),
    ReadError,
> {
    let file = File::open(path)?;
    let mut bgzf = noodles_bgzf::io::Reader::new(BufReader::new(file));

    // --- header ---
    if read_exact_stream(&mut bgzf, 4)? != b"BAM\x01" {
        return Err(ReadError::BadMagic);
    }
    let l_text = i32::from_le_bytes(read_exact_stream(&mut bgzf, 4)?.try_into().unwrap());
    let _text = read_exact_stream(&mut bgzf, l_text as usize)?;
    let n_ref = i32::from_le_bytes(read_exact_stream(&mut bgzf, 4)?.try_into().unwrap());
    let mut contigs = Vec::with_capacity(n_ref as usize);
    for _ in 0..n_ref {
        let l_name = i32::from_le_bytes(read_exact_stream(&mut bgzf, 4)?.try_into().unwrap());
        let name = read_exact_stream(&mut bgzf, l_name as usize)?;
        let _l_ref = read_exact_stream(&mut bgzf, 4)?;
        contigs.push(String::from_utf8_lossy(&name[..name.len().saturating_sub(1)]).into_owned());
    }

    let iter = std::iter::from_fn(move || {
        // block_size
        let mut sz = [0u8; 4];
        match bgzf.read(&mut sz) {
            Ok(0) => return None,
            Ok(4) => {}
            Ok(_) => return Some(Err(ReadError::Truncated)),
            Err(e) => return Some(Err(ReadError::Io(e))),
        }
        let block_size = u32::from_le_bytes(sz) as usize;
        match read_block(&mut bgzf, block_size) {
            Ok(rec) => Some(Ok(rec)),
            Err(e) => Some(Err(e)),
        }
    });
    Ok((contigs, iter))
}

fn read_block(r: &mut impl Read, block_size: usize) -> Result<InRecord, ReadError> {
    let raw = read_exact_stream(r, block_size)?;
    let mut p = Parser { buf: raw, off: 0 };

    let ref_id = p.i32()?;
    let pos = p.i32()? as i64;
    let l_read_name = p.u8()? as usize;
    let mapq = p.u8()?;
    let _bin = p.u16()?;
    let n_cigar = p.u16()? as usize;
    let flag = p.u16()?;
    let l_seq = p.u32()? as usize;
    let _next_ref = p.i32()?;
    let _next_pos = p.i32()?;
    let _tlen = p.i32()?;

    let name_bytes = p.take(l_read_name)?;
    let name = String::from_utf8_lossy(&name_bytes[..l_read_name.saturating_sub(1)]).into_owned();

    let mut cigar = Vec::with_capacity(n_cigar);
    for _ in 0..n_cigar {
        let v = p.u32()?;
        let (kind, len) = (v & 0xf, v >> 4);
        let op = match kind {
            0 => CigarOp::Match(len),
            1 => CigarOp::Ins(len),
            2 => CigarOp::Del(len),
            3 => CigarOp::RefSkip(len),
            4 => CigarOp::SoftClip(len),
            k => return Err(ReadError::BadCigarKind(k)),
        };
        cigar.push(op);
    }

    let packed = p.take(l_seq.div_ceil(2))?;
    let mut seq = Vec::with_capacity(l_seq);
    for i in 0..l_seq {
        let byte = packed[i / 2];
        let nib = if i % 2 == 0 { byte >> 4 } else { byte & 0xf };
        seq.push(SEQ_DECODE[nib as usize]);
    }

    let qual = p.take(l_seq)?.to_vec();

    // --- aux tags (parse A/i/Z; skip other types by width) ---
    let mut tags = Vec::new();
    while p.off < p.buf.len() {
        let t0 = p.u8()?;
        let t1 = p.u8()?;
        let ty = p.u8()?;
        match ty {
            b'A' => tags.push(RawTag([t0, t1], TagValue::Char(p.u8()?))),
            b'i' | b'I' => tags.push(RawTag([t0, t1], TagValue::Int(p.i32()?))),
            b'Z' | b'H' => {
                let start = p.off;
                let mut end = start;
                while end < p.buf.len() && p.buf[end] != 0 {
                    end += 1;
                }
                if end >= p.buf.len() {
                    return Err(ReadError::Truncated);
                }
                let s = String::from_utf8_lossy(&p.buf[start..end]).into_owned();
                p.off = end + 1;
                tags.push(RawTag([t0, t1], TagValue::Str(s)));
            }
            b'c' | b'C' => {
                p.u8()?;
            }
            b's' | b'S' => {
                p.u16()?;
            }
            b'f' => {
                p.u32()?;
            }
            b'd' => {
                p.take(8)?;
            }
            b'B' => {
                let sub = p.u8()?;
                let n = p.u32()? as usize;
                let width = match sub {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => return Err(ReadError::Truncated),
                };
                p.take(n * width)?;
            }
            _ => return Err(ReadError::Truncated),
        }
    }

    Ok(InRecord {
        name,
        flag,
        mapq,
        ref_id,
        pos,
        cigar,
        seq,
        qual,
        tags,
    })
}
