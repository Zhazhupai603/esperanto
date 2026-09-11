//! paidx v1 index serialization (byte-level deterministic layout).
//!
//! ```text
//! magic[8]="PAIDXfmt" | version u32 | k u32 | w u32 | freq_cutoff u32
//! contig_count u32 | kmer_count u64 | positions_count u64 | fasta_sha256[32]
//! per contig: name_len u32 | name | seq_len u32 | packed_len u32 | packed
//!           | n_interval_count u32 | n×(start u32,end u32)
//! pad1 = (8 − (60+contig_bytes)%8)%8
//! kmers ×u64 | offsets ×u64 | counts ×u32
//! pad2 = (8 − (60+contig_bytes+pad1+kmer_count*16+kmer_count*4)%8)%8
//! positions ×u64
//! ```
//!
//! All little-endian, no timestamps or randomness: two builds from the same
//! input are byte-identical. Loading mmaps and leaks the file; tables point
//! straight into mapped pages via `bytemuck::cast_slice` (alignment is
//! guaranteed by the pad computations).

use crate::error::AlignError;
use crate::fasta::{Contig, Interval, Reference};
use crate::index::Index;
use crate::seed::SeedParams;
use memmap2::Mmap;
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// Format magic (frozen).
pub const PAIDX_MAGIC: &[u8; 8] = b"PAIDXfmt";
/// Supported format version.
pub const PAIDX_VERSION: u32 = 1;

const HEADER_FIXED: u64 = 60; // 8+4+4+4+4+4+32 (two u64 counts ≡ 0 mod 8)
const COUNTS_SIZE: u64 = 16; // kmer_count + positions_count

/// Serialize an index to `path` in paidx v1 layout.
pub fn save(index: &Index, path: &Path) -> Result<(), AlignError> {
    let bytes = encode(index)?;
    let mut file = File::create(path).map_err(|_| AlignError::IndexIo)?;
    file.write_all(&bytes).map_err(|_| AlignError::IndexIo)?;
    file.flush().map_err(|_| AlignError::IndexIo)?;
    Ok(())
}

/// Encode an index into in-memory paidx v1 bytes.
pub fn encode(index: &Index) -> Result<Vec<u8>, AlignError> {
    let contig_count = index.reference.contigs.len() as u32;
    let kmer_count = index.kmers.len() as u64;
    let positions_count = index.positions.len() as u64;

    // contig_bytes computed analytically
    let mut contig_bytes: u64 = 0;
    for c in &index.reference.contigs {
        contig_bytes += 4
            + c.name.len() as u64
            + 4
            + 4
            + c.packed.len() as u64
            + 4
            + (c.n_intervals.len() as u64) * 8;
    }
    let pad1 = (8 - (HEADER_FIXED + contig_bytes) % 8) % 8;
    let pad2 =
        (8 - (HEADER_FIXED + contig_bytes + pad1 + kmer_count * 16 + kmer_count * 4) % 8) % 8;

    let total = HEADER_FIXED
        + COUNTS_SIZE
        + contig_bytes
        + pad1
        + kmer_count * 16
        + kmer_count * 4
        + pad2
        + positions_count * 8;

    let mut out: Vec<u8> = Vec::with_capacity(total as usize);
    out.extend_from_slice(PAIDX_MAGIC);
    out.extend_from_slice(&index.version.to_le_bytes());
    out.extend_from_slice(&index.params.k.to_le_bytes());
    out.extend_from_slice(&index.params.w.to_le_bytes());
    out.extend_from_slice(&index.freq_cutoff.to_le_bytes());
    out.extend_from_slice(&contig_count.to_le_bytes());
    out.extend_from_slice(&kmer_count.to_le_bytes());
    out.extend_from_slice(&positions_count.to_le_bytes());
    out.extend_from_slice(&index.reference.fasta_sha256);
    for c in &index.reference.contigs {
        out.extend_from_slice(&(c.name.len() as u32).to_le_bytes());
        out.extend_from_slice(c.name.as_bytes());
        out.extend_from_slice(&c.len.to_le_bytes());
        out.extend_from_slice(&(c.packed.len() as u32).to_le_bytes());
        out.extend_from_slice(&c.packed);
        out.extend_from_slice(&(c.n_intervals.len() as u32).to_le_bytes());
        for iv in &c.n_intervals {
            out.extend_from_slice(&iv.start.to_le_bytes());
            out.extend_from_slice(&iv.end.to_le_bytes());
        }
    }
    out.resize(out.len() + pad1 as usize, 0);
    for kmer in index.kmers {
        out.extend_from_slice(&kmer.to_le_bytes());
    }
    for off in index.offsets {
        out.extend_from_slice(&off.to_le_bytes());
    }
    for cnt in index.counts {
        out.extend_from_slice(&cnt.to_le_bytes());
    }
    out.resize(out.len() + pad2 as usize, 0);
    for p in index.positions {
        out.extend_from_slice(&p.to_le_bytes());
    }
    debug_assert_eq!(out.len() as u64, total);
    Ok(out)
}

/// Load a paidx v1 index: mmap + leak; strict validation throughout.
pub fn load(path: &Path) -> Result<Index, AlignError> {
    let file = File::open(path).map_err(|_| AlignError::IndexIo)?;
    let mmap = unsafe { Mmap::map(&file) }.map_err(|_| AlignError::IndexIo)?;
    let mmap: &'static Mmap = Box::leak(Box::new(mmap));
    let data: &'static [u8] = &mmap[..];
    decode(data, &path.display().to_string())
}

/// Decode paidx v1 bytes (already leaked/mapped, `'static`) into an `Index`.
pub fn decode(data: &'static [u8], file: &str) -> Result<Index, AlignError> {
    let bad = |msg: &str| AlignError::IndexFormat {
        msg: msg.to_string(),
    };
    let mut cur = Cursor::new(data);

    let magic = cur.take(8, &bad)?;
    if magic != PAIDX_MAGIC {
        return Err(bad("bad magic"));
    }
    let version = u32::from_le_bytes(cur.take_u32(&bad)?);
    if version != PAIDX_VERSION {
        return Err(AlignError::IndexVersion {
            file: file.to_string(),
            supported: PAIDX_VERSION,
        });
    }
    let k = u32::from_le_bytes(cur.take_u32(&bad)?);
    let w = u32::from_le_bytes(cur.take_u32(&bad)?);
    let freq_cutoff = u32::from_le_bytes(cur.take_u32(&bad)?);
    let contig_count = u32::from_le_bytes(cur.take_u32(&bad)?) as usize;
    let kmer_count = u64::from_le_bytes(cur.take_u64(&bad)?) as usize;
    let positions_count = u64::from_le_bytes(cur.take_u64(&bad)?) as usize;
    let sha_slice = cur.take(32, &bad)?;
    let mut fasta_sha256 = [0u8; 32];
    fasta_sha256.copy_from_slice(sha_slice);

    // Contigs.
    let mut contigs: Vec<Contig> = Vec::with_capacity(contig_count);
    let mut contig_bytes: u64 = 0;
    for _ in 0..contig_count {
        let name_len = u32::from_le_bytes(cur.take_u32(&bad)?) as usize;
        let name_bytes = cur.take(name_len, &bad)?;
        let seq_len = u32::from_le_bytes(cur.take_u32(&bad)?);
        let packed_len = u32::from_le_bytes(cur.take_u32(&bad)?) as usize;
        if packed_len != seq_len.div_ceil(4) as usize {
            return Err(bad("packed_len != ceil(seq_len/4)"));
        }
        let packed = cur.take(packed_len, &bad)?.to_vec();
        // trailing bits of the final byte must be zero
        if seq_len % 4 != 0 {
            let used_bits = 2 * (seq_len % 4);
            if let Some(last) = packed.last() {
                if last & (0xFFu8 >> used_bits) != 0 {
                    return Err(bad("non-zero tail bits in packed contig"));
                }
            }
        }
        let n_count = u32::from_le_bytes(cur.take_u32(&bad)?) as usize;
        let mut n_intervals = Vec::with_capacity(n_count);
        for _ in 0..n_count {
            let start = u32::from_le_bytes(cur.take_u32(&bad)?);
            let end = u32::from_le_bytes(cur.take_u32(&bad)?);
            n_intervals.push(Interval { start, end });
        }
        contig_bytes += 4 + name_len as u64 + 4 + 4 + packed_len as u64 + 4 + (n_count as u64) * 8;
        contigs.push(Contig {
            name: String::from_utf8_lossy(name_bytes).into_owned(),
            len: seq_len,
            packed,
            n_intervals,
        });
    }

    // pad1 must exist and be zero.
    let pad1 = (8 - (HEADER_FIXED + contig_bytes) % 8) % 8;
    let pad1_bytes = cur.take(pad1 as usize, &bad)?;
    if pad1_bytes.iter().any(|&b| b != 0) {
        return Err(bad("pad1 not zero"));
    }

    // Tables: kmers/offsets/counts (u64/u64/u32).
    let kmers: &'static [u64] = cast_table(&mut cur, kmer_count * 8, &bad, 8)?;
    let offsets: &'static [u64] = cast_table(&mut cur, kmer_count * 8, &bad, 8)?;
    let counts: &'static [u32] = cast_table(&mut cur, kmer_count * 4, &bad, 4)?;

    // pad2 must exist and be zero.
    let pad2 = (8
        - (HEADER_FIXED + contig_bytes + pad1 + kmer_count as u64 * 16 + kmer_count as u64 * 4)
            % 8)
        % 8;
    let pad2_bytes = cur.take(pad2 as usize, &bad)?;
    if pad2_bytes.iter().any(|&b| b != 0) {
        return Err(bad("pad2 not zero"));
    }

    let positions: &'static [u64] = cast_table(&mut cur, positions_count * 8, &bad, 8)?;

    if cur.remaining() != 0 {
        return Err(bad("trailing bytes after positions table"));
    }

    // Table consistency.
    if kmers.windows(2).any(|w| w[0] >= w[1]) {
        return Err(bad("kmers not strictly ascending"));
    }
    if offsets.len() != counts.len() || offsets.len() != kmers.len() {
        return Err(bad("table length mismatch"));
    }
    for i in 0..offsets.len() {
        let o = offsets[i];
        let c = counts[i] as u64;
        if o + c > positions.len() as u64 {
            return Err(bad("offset+count exceeds positions table"));
        }
    }

    let reference: &'static Reference = Box::leak(Box::new(Reference {
        contigs,
        fasta_sha256,
    }));

    Ok(Index {
        params: SeedParams { k, w },
        version,
        freq_cutoff,
        reference,
        kmers,
        offsets,
        counts,
        positions,
    })
}

struct Cursor {
    data: &'static [u8],
    pos: usize,
}

impl Cursor {
    fn new(data: &'static [u8]) -> Cursor {
        Cursor { data, pos: 0 }
    }

    fn remaining(&self) -> usize {
        self.data.len() - self.pos
    }

    fn take(
        &mut self,
        n: usize,
        bad: &dyn Fn(&str) -> AlignError,
    ) -> Result<&'static [u8], AlignError> {
        if self.remaining() < n {
            return Err(bad("unexpected end of file"));
        }
        let out = &self.data[self.pos..self.pos + n];
        self.pos += n;
        Ok(out)
    }

    fn take_u32(&mut self, bad: &dyn Fn(&str) -> AlignError) -> Result<[u8; 4], AlignError> {
        let bytes = self.take(4, bad)?;
        Ok([bytes[0], bytes[1], bytes[2], bytes[3]])
    }

    fn take_u64(&mut self, bad: &dyn Fn(&str) -> AlignError) -> Result<[u8; 8], AlignError> {
        let bytes = self.take(8, bad)?;
        Ok([
            bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5], bytes[6], bytes[7],
        ])
    }
}

/// Take `n` bytes aligned to `align` and cast them in place to the element
/// type `T` (caller chooses the output slice type).
fn cast_table<T: bytemuck::AnyBitPattern>(
    cur: &mut Cursor,
    n: usize,
    bad: &dyn Fn(&str) -> AlignError,
    align: usize,
) -> Result<&'static [T], AlignError> {
    if !n.is_multiple_of(std::mem::size_of::<T>()) {
        return Err(bad("table size not a multiple of element size"));
    }
    if !(cur.data.as_ptr() as usize + cur.pos).is_multiple_of(align) {
        return Err(bad("table not aligned"));
    }
    let bytes = cur.take(n, bad)?;
    Ok(bytemuck::cast_slice(bytes))
}
