//! `.baln` reader for the call side. Reads the flat binary stream written by
//! align (magic + contig table + `block_size u32 + BAM-disk-record bytes`),
//! parses each record into fields, ready to fill `bam1_t` via `bam_set1`.
//!
//! BAM-disk record layout (little-endian, 32-byte core + data):
//!   core: refID i32 | pos i32 | l_read_name u8 | mapq u8 | bin u16 |
//!         n_cigar u16 | flag u16 | l_seq i32 | next_refID i32 | next_pos i32 | tlen i32
//!   data: qname(l_read_name, NUL-term) | cigar[n_cigar] u32((len<<4)|code) |
//!         seq(4-bit packed) | qual(l_seq bytes)

use crate::error::CallError;
use std::io::{BufReader, Read};
use std::path::Path;

pub struct BalnReader {
    r: BufReader<std::fs::File>,
    pub contigs: Vec<String>,
    /// Consumed byte offset (including header); a record's start offset = the value before next() returns it.
    pub pos_bytes: u64,
}

/// One parsed .baln record, fields ready for `bam_set1`.
pub struct BalnRecord {
    pub flag: u16,
    pub tid: i32,
    pub pos: i64,
    pub mapq: u8,
    pub mtid: i32,
    pub mpos: i64,
    pub isize: i64,
    pub l_seq: usize,
    pub name: Vec<u8>,      // without trailing NUL
    pub cigar: Vec<u32>,    // (len<<4)|code, as bam_set1 expects
    pub seq_ascii: Vec<u8>, // decoded to ASCII ACGTN
    pub qual: Vec<u8>,      // raw phred [0,93] (BAM-disk convention, same as htslib qual())
    /// EK:Z aux (edit-evidence string; None when absent). The caller's dirty-read gating depends on it.
    pub ek: Option<String>,
    /// RE:Z aux (rescue provenance; `collapsed` marks an alphabet-ambiguous placement).
    pub re: Option<String>,
}

/// Parse a `Z`-type string aux field from the aux region (BAM aux TLV:
/// tag[2]+type[1]+value). Returns the string value of the first `want` tag.
fn aux_string(aux: &[u8], want: &[u8; 2]) -> Option<String> {
    let mut i = 0usize;
    while i + 3 <= aux.len() {
        let tag = &aux[i..i + 2];
        let ty = aux[i + 2];
        let vstart = i + 3;
        let vlen = match ty {
            b'A' | b'c' | b'C' => 1,
            b's' | b'S' => 2,
            b'i' | b'I' | b'f' => 4,
            b'd' => 8,
            b'Z' | b'H' => {
                let mut e = vstart;
                while e < aux.len() && aux[e] != 0 {
                    e += 1;
                }
                if tag == want {
                    return Some(
                        String::from_utf8_lossy(&aux[vstart..e.min(aux.len())]).into_owned(),
                    );
                }
                e - vstart + 1
            }
            b'B' => {
                if vstart + 5 > aux.len() {
                    return None;
                }
                let sub = aux[vstart];
                let n =
                    u32::from_le_bytes(aux[vstart + 1..vstart + 5].try_into().unwrap()) as usize;
                let esz = match sub {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => return None,
                };
                1 + 4 + n * esz
            }
            _ => return None,
        };
        i = vstart + vlen;
    }
    None
}

const SEQ4: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

impl BalnReader {
    pub fn open(path: &Path) -> Result<Self, CallError> {
        let f = std::fs::File::open(path).map_err(|e| CallError::Io {
            path: path.display().to_string(),
            source: e,
        })?;
        let mut r = BufReader::with_capacity(8 << 20, f);
        let mut magic = [0u8; 8];
        r.read_exact(&mut magic).map_err(|e| CallError::Io {
            path: path.display().to_string(),
            source: e,
        })?;
        if &magic != b"ESPBALN\x01" {
            return Err(CallError::Io {
                path: path.display().to_string(),
                source: std::io::Error::new(io::ErrorKind::InvalidData, "bad .baln magic"),
            });
        }
        let n = read_u32(&mut r, path)?;
        let mut contigs = Vec::with_capacity(n as usize);
        for _ in 0..n {
            let mut lb = [0u8; 1];
            r.read_exact(&mut lb).map_err(ioe(path))?;
            let mut name = vec![0u8; lb[0] as usize];
            r.read_exact(&mut name).map_err(ioe(path))?;
            contigs.push(String::from_utf8_lossy(&name).into_owned());
        }
        let pos_bytes = 8 + 4 + contigs.iter().map(|c| 1 + c.len()).sum::<usize>() as u64;
        Ok(Self {
            r,
            contigs,
            pos_bytes,
        })
    }

    /// Read the next record; Ok(None) at EOF.
    pub fn next_record(&mut self) -> Result<Option<BalnRecord>, CallError> {
        let mut sz = [0u8; 4];
        match self.r.read_exact(&mut sz) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
            Err(e) => {
                return Err(CallError::Io {
                    path: ".baln".into(),
                    source: e,
                })
            }
        }
        let block = u32::from_le_bytes(sz) as usize;
        let mut buf = vec![0u8; block];
        self.r.read_exact(&mut buf).map_err(ioe_path())?;
        self.pos_bytes += 4 + block as u64;

        Ok(Some(parse_block(&buf)))
    }
}

/// Parse one .baln record body (block_size already stripped; 32-byte core + data).
fn parse_block(buf: &[u8]) -> BalnRecord {
    {
        // --- 32-byte core ---
        let ref_id = i32::from_le_bytes(buf[0..4].try_into().unwrap());
        let pos = i32::from_le_bytes(buf[4..8].try_into().unwrap());
        let l_read_name = buf[8] as usize;
        let mapq = buf[9];
        let n_cigar = u16::from_le_bytes(buf[12..14].try_into().unwrap()) as usize;
        let flag = u16::from_le_bytes(buf[14..16].try_into().unwrap());
        let l_seq = i32::from_le_bytes(buf[16..20].try_into().unwrap()) as usize;
        let next_ref = i32::from_le_bytes(buf[20..24].try_into().unwrap());
        let next_pos = i32::from_le_bytes(buf[24..28].try_into().unwrap());
        let tlen = i32::from_le_bytes(buf[28..32].try_into().unwrap());

        // --- data: qname | cigar | seq(4bit) | qual ---
        let d = &buf[32..];
        let name = d[..l_read_name.saturating_sub(1)].to_vec();
        let mut off = l_read_name;
        let mut cigar = Vec::with_capacity(n_cigar);
        for i in 0..n_cigar {
            let b = off + i * 4;
            cigar.push(u32::from_le_bytes(d[b..b + 4].try_into().unwrap()));
        }
        off += n_cigar * 4;
        let seq_packed_len = l_seq.div_ceil(2);
        let packed = &d[off..off + seq_packed_len];
        off += seq_packed_len;
        let qual = d[off..off + l_seq].to_vec();
        let ek = aux_string(&d[off + l_seq..], b"EK");
        let re = aux_string(&d[off + l_seq..], b"RE");
        let mut seq_ascii = Vec::with_capacity(l_seq);
        for i in 0..l_seq {
            let byte = packed[i / 2];
            let nib = if i % 2 == 0 { byte >> 4 } else { byte & 0xF };
            seq_ascii.push(SEQ4[nib as usize]);
        }

        BalnRecord {
            flag,
            tid: ref_id,
            pos: pos as i64,
            mapq,
            mtid: next_ref,
            mpos: next_pos as i64,
            isize: tlen as i64,
            l_seq,
            name,
            cigar,
            seq_ascii,
            qual,
            ek,
            re,
        }
    }
}

/// Random-access read of one record (pread-style, disturbs no cursor; off comes from build_index).
/// EOF (off past end of file) → Ok(None).
pub fn read_record_at(file: &std::fs::File, off: u64) -> Result<Option<BalnRecord>, CallError> {
    use std::os::unix::fs::FileExt;
    let mut sz = [0u8; 4];
    match file.read_exact_at(&mut sz, off) {
        Ok(()) => {}
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => {
            return Err(CallError::Io {
                path: ".baln".into(),
                source: e,
            });
        }
    }
    let block = u32::from_le_bytes(sz) as usize;
    let mut buf = vec![0u8; block];
    file.read_exact_at(&mut buf, off + 4)
        .map_err(|e| CallError::Io {
            path: ".baln".into(),
            source: e,
        })?;
    Ok(Some(parse_block(&buf)))
}

/// .baln coordinate index (single-pass scan product, shared read-only across block tasks).
pub struct BalnIndex {
    /// (tid, pos, block start byte offset, reference span) sorted by (tid,pos).
    /// A block task [cs,ce) selects `pos < ce && pos + span > cs` (aligned with htslib fetch
    /// overlap semantics; span = M/D/N/=/X reference consumption, same as bam_endpos).
    pub idx: Vec<(i32, i64, u64, i64)>,
    /// Global maximum span (block-window backward extension).
    pub max_span: i64,
    /// contig names (tid order).
    pub contigs: Vec<String>,
    /// Per-contig derived length = max(pos+span) (coverage upper bound; read-free regions have no sites,
    /// so the tail-block read set is equivalent to the BAM header length convention — interior block boundaries stay 32M-aligned).
    pub tid_len: Vec<u32>,
}

impl BalnReader {
    /// Build the coordinate index in a single pass. The pos_bytes value before next() is the record start offset;
    /// the read-ordered stream does not guarantee coordinate order → collect then sort.
    pub fn build_index(path: &Path) -> Result<BalnIndex, CallError> {
        let mut r = BalnReader::open(path)?;
        let contigs = std::mem::take(&mut r.contigs);
        let mut idx: Vec<(i32, i64, u64, i64)> = Vec::new();
        let mut max_span: i64 = 0;
        loop {
            let start = r.pos_bytes;
            let Some(rec) = r.next_record()? else { break };
            if rec.tid >= 0 && rec.flag & 0x4 == 0 {
                let mut span = 0i64;
                for &c in &rec.cigar {
                    // BAM op codes: 0=M 1=I 2=D 3=N 4=S 5=H 6=P 7== 8=X
                    if matches!(c & 0xF, 0 | 2 | 3 | 7 | 8) {
                        span += (c >> 4) as i64;
                    }
                }
                max_span = max_span.max(span);
                idx.push((rec.tid, rec.pos, start, span));
            }
        }
        let mut tid_len = vec![0u32; contigs.len()];
        for e in &idx {
            let end = (e.1 + e.3).max(0) as u32;
            let t = e.0 as usize;
            if t < tid_len.len() && end > tid_len[t] {
                tid_len[t] = end;
            }
        }
        idx.sort_unstable();
        Ok(BalnIndex {
            idx,
            max_span,
            contigs,
            tid_len,
        })
    }
}

fn read_u32(r: &mut impl Read, path: &Path) -> Result<u32, CallError> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b).map_err(|e| CallError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    Ok(u32::from_le_bytes(b))
}

fn ioe(path: &Path) -> impl FnOnce(io::Error) -> CallError + '_ {
    move |e| CallError::Io {
        path: path.display().to_string(),
        source: e,
    }
}
fn ioe_path() -> impl FnOnce(io::Error) -> CallError {
    |e| CallError::Io {
        path: ".baln".into(),
        source: e,
    }
}

use std::io;
