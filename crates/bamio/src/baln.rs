//! `.baln` reader — the read side of the flat binary alignment channel
//! written by align (magic + contig table + `block_size u32` + BAM-disk-record
//! bytes). BAM stays the user-facing artifact; `.baln` is the fast internal
//! interface shared by scan and score/pileup.
//!
//! BAM-disk record layout (little-endian, 32-byte core + data):
//!   core: refID i32 | pos i32 | l_read_name u8 | mapq u8 | bin u16 |
//!         n_cigar u16 | flag u16 | l_seq i32 | next_refID i32 | next_pos i32 | tlen i32
//!   data: qname(l_read_name, NUL-term) | cigar[n_cigar] u32((len<<4)|code) |
//!         seq(4-bit packed) | qual(l_seq bytes) | aux

use std::io::{self, BufReader, Read};
use std::path::Path;

/// Format magic + version byte (writer contract in the map crate).
pub const MAGIC: &[u8; 8] = b"ESPBALN\x01";

const SEQ4: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

fn invalid(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.to_string())
}

/// Sequential `.baln` reader.
pub struct BalnReader {
    r: BufReader<std::fs::File>,
    /// Contig name table (tid order).
    pub contigs: Vec<String>,
    /// Consumed byte offset (including header); a record's start offset = the
    /// value before `next_record` returns it.
    pub pos_bytes: u64,
}

/// One parsed `.baln` record.
#[derive(Clone)]
pub struct BalnRecord {
    /// SAM flag.
    pub flag: u16,
    /// Reference id (header order; -1 = unmapped).
    pub tid: i32,
    /// 0-based leftmost position.
    pub pos: i64,
    /// Mapping quality.
    pub mapq: u8,
    /// Mate reference id.
    pub mtid: i32,
    /// Mate position.
    pub mpos: i64,
    /// Template length.
    pub isize: i64,
    /// Sequence length.
    pub l_seq: usize,
    /// Read name (without trailing NUL).
    pub name: Vec<u8>,
    /// CIGAR as raw `(len<<4)|code` words.
    pub cigar: Vec<u32>,
    /// Sequence decoded to ASCII (IUPAC).
    pub seq_ascii: Vec<u8>,
    /// Raw phred qualities (BAM-disk convention, same as htslib `qual()`).
    pub qual: Vec<u8>,
    /// EK:Z aux (edit-evidence string; None when absent).
    pub ek: Option<String>,
    /// RE:Z aux (rescue provenance; `collapsed` marks an alphabet-ambiguous placement).
    pub re: Option<String>,
}

/// Parse a `Z`-type string aux field from the aux region (BAM aux TLV:
/// tag[2]+type[1]+value). Returns the value of the first `want` tag.
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

/// Parse one record body (block_size already stripped; 32-byte core + data).
fn parse_block(buf: &[u8]) -> BalnRecord {
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

    // --- data: qname | cigar | seq(4bit) | qual | aux ---
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
    let aux = &d[off + l_seq..];
    let ek = aux_string(aux, b"EK");
    let re = aux_string(aux, b"RE");
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

impl BalnReader {
    /// Open a `.baln` stream and read the header.
    pub fn open(path: &Path) -> Result<Self, io::Error> {
        let f = std::fs::File::open(path)?;
        let mut r = BufReader::with_capacity(8 << 20, f);
        let mut magic = [0u8; 8];
        r.read_exact(&mut magic)?;
        if &magic != MAGIC {
            return Err(invalid("bad .baln magic"));
        }
        let mut n = [0u8; 4];
        r.read_exact(&mut n)?;
        let n = u32::from_le_bytes(n);
        let mut contigs = Vec::with_capacity(n as usize);
        for _ in 0..n {
            let mut lb = [0u8; 1];
            r.read_exact(&mut lb).map_err(io::Error::other)?;
            let mut name = vec![0u8; lb[0] as usize];
            r.read_exact(&mut name).map_err(io::Error::other)?;
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
    pub fn next_record(&mut self) -> Result<Option<BalnRecord>, io::Error> {
        let mut sz = [0u8; 4];
        match self.r.read_exact(&mut sz) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
            Err(e) => return Err(e),
        }
        let block = u32::from_le_bytes(sz) as usize;
        let mut buf = vec![0u8; block];
        self.r.read_exact(&mut buf)?;
        self.pos_bytes += 4 + block as u64;

        Ok(Some(parse_block(&buf)))
    }
}

/// Random-access read of one record (pread-style, disturbs no cursor; `off`
/// comes from [`BalnReader::build_index`]). EOF → Ok(None).
pub fn read_record_at(file: &std::fs::File, off: u64) -> Result<Option<BalnRecord>, io::Error> {
    use std::os::unix::fs::FileExt;
    let mut sz = [0u8; 4];
    match file.read_exact_at(&mut sz, off) {
        Ok(()) => {}
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e),
    }
    let block = u32::from_le_bytes(sz) as usize;
    let mut buf = vec![0u8; block];
    file.read_exact_at(&mut buf, off + 4)?;
    Ok(Some(parse_block(&buf)))
}

/// `.baln` coordinate index (single-pass build, shared read-only across
/// region tasks).
pub struct BalnIndex {
    /// (tid, pos, block start byte offset, reference span) sorted by (tid,pos).
    /// A region [cs,ce) selects `pos < ce && pos + span > cs` (htslib fetch
    /// overlap semantics; span = M/D/N/=/X reference consumption).
    pub idx: Vec<(i32, i64, u64, i64)>,
    /// Global maximum span (backward window extension).
    pub max_span: i64,
    /// Contig names (tid order).
    pub contigs: Vec<String>,
    /// Per-contig derived length = max(pos+span) (coverage upper bound).
    pub tid_len: Vec<u32>,
}

impl BalnReader {
    /// Build the coordinate index in a single pass. The read-ordered stream
    /// does not guarantee coordinate order → collect then sort.
    pub fn build_index(path: &Path) -> Result<BalnIndex, io::Error> {
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

/// Convert a parsed `.baln` record into an htslib record (for consumers
/// pinned to rust-htslib semantics, e.g. the pileup engine). The EK/RE aux
/// tags are re-attached; mate fields follow the source record.
impl BalnRecord {
    /// Build the htslib view. Cigar codes map 1:1 from BAM op codes.
    pub fn to_htslib_record(&self) -> Result<rust_htslib::bam::Record, io::Error> {
        use rust_htslib::bam::record::{Cigar, CigarString};
        let mut ops: Vec<Cigar> = Vec::with_capacity(self.cigar.len());
        for &c in &self.cigar {
            let len = c >> 4;
            let op = match c & 0xF {
                0 => Cigar::Match(len),
                1 => Cigar::Ins(len),
                2 => Cigar::Del(len),
                3 => Cigar::RefSkip(len),
                4 => Cigar::SoftClip(len),
                5 => Cigar::HardClip(len),
                6 => Cigar::Pad(len),
                7 => Cigar::Equal(len),
                8 => Cigar::Diff(len),
                other => {
                    return Err(invalid(&format!(
                        "baln record: unknown cigar op code {other}"
                    )));
                }
            };
            ops.push(op);
        }
        let mut rec = rust_htslib::bam::Record::new();
        rec.set(
            &self.name,
            Some(&CigarString(ops)),
            &self.seq_ascii,
            &self.qual,
        );
        rec.set_tid(self.tid);
        rec.set_pos(self.pos);
        rec.set_mapq(self.mapq);
        rec.set_flags(self.flag);
        rec.set_mtid(self.mtid);
        rec.set_mpos(self.mpos);
        if let Some(re) = &self.re {
            rec.push_aux(b"RE", rust_htslib::bam::record::Aux::String(re))
                .map_err(|e| invalid(&format!("baln record RE aux: {e}")))?;
        }
        if let Some(ek) = &self.ek {
            rec.push_aux(b"EK", rust_htslib::bam::record::Aux::String(ek))
                .map_err(|e| invalid(&format!("baln record EK aux: {e}")))?;
        }
        Ok(rec)
    }
}
