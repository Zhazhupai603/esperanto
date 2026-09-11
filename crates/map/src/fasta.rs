//! Reference genome storage in 2-bit packed form.
//!
//! Each byte packs 4 bases, most-significant base first: base at `pos` lives
//! in byte `pos / 4` at shift `6 - 2 * (pos % 4)`. Tail bits of the final byte
//! (when `len % 4 != 0`) are zero. `N` cannot be encoded in 2 bits, so N runs
//! are kept as half-open intervals per contig and consulted on every decode.

use crate::error::AlignError;
use sha2::{Digest, Sha256};
use std::fs;
use std::path::Path;

/// A single nucleotide.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Base {
    /// Adenine.
    A,
    /// Cytosine.
    C,
    /// Guanine.
    G,
    /// Thymine.
    T,
    /// Any ambiguity code.
    N,
}

impl Base {
    /// Case-insensitive ASCII decoder; anything unrecognized becomes `N`.
    pub fn from_ascii(b: u8) -> Base {
        match b {
            b'A' | b'a' => Base::A,
            b'C' | b'c' => Base::C,
            b'G' | b'g' => Base::G,
            b'T' | b't' => Base::T,
            _ => Base::N,
        }
    }

    /// 2-bit code; `None` for `N`.
    pub fn code(self) -> Option<u8> {
        match self {
            Base::A => Some(0),
            Base::C => Some(1),
            Base::G => Some(2),
            Base::T => Some(3),
            Base::N => None,
        }
    }
}

/// Half-open `[start, end)` interval of N bases.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Interval {
    /// First N position (0-based).
    pub start: u32,
    /// One past the last N position.
    pub end: u32,
}

/// One contig: name, length, 2-bit packed bases and N intervals.
pub struct Contig {
    /// Contig name (first whitespace-delimited token of the FASTA header).
    pub name: String,
    /// Sequence length in bases.
    pub len: u32,
    /// Packed bases; `ceil(len / 4)` bytes, high bases first within a byte.
    pub packed: Vec<u8>,
    /// Merged half-open N runs, sorted ascending.
    pub n_intervals: Vec<Interval>,
}

impl Contig {
    /// Decode `[start, end)` into ASCII bases, writing `N` over N intervals.
    ///
    /// Panics if the range is out of bounds; callers use validated ranges.
    pub fn decode_into(&self, start: u32, end: u32, out: &mut Vec<u8>) {
        for pos in start..end {
            out.push(self.base(pos).to_ascii());
        }
    }

    /// Append `[start, end)` ASCII bases to `out` (same as [`decode_into`]).
    pub fn decode_append(&self, start: u32, end: u32, out: &mut Vec<u8>) {
        self.decode_into(start, end, out);
    }

    /// ASCII slice of `[start, end)`; `N` wherever an N interval covers.
    pub fn slice_ascii(&self, start: u32, end: u32) -> Vec<u8> {
        let mut out = Vec::with_capacity((end - start) as usize);
        self.decode_into(start, end, &mut out);
        out
    }

    /// Bounds-checked base access: None when `pos` is past the contig end.
    pub fn base_checked(&self, pos: u32) -> Option<Base> {
        if pos >= self.len {
            None
        } else {
            Some(self.base(pos))
        }
    }

    /// Base at `pos`: N-interval table first, then 2-bit unpack.
    pub fn base(&self, pos: u32) -> Base {
        if self.is_n(pos) {
            return Base::N;
        }
        let byte = self.packed[(pos >> 2) as usize];
        let shift = 6 - 2 * (pos & 3);
        match (byte >> shift) & 3 {
            0 => Base::A,
            1 => Base::C,
            2 => Base::G,
            _ => Base::T,
        }
    }

    fn is_n(&self, pos: u32) -> bool {
        let mut lo = 0usize;
        let mut hi = self.n_intervals.len();
        while lo < hi {
            let mid = (lo + hi) / 2;
            let iv = &self.n_intervals[mid];
            if pos < iv.start {
                hi = mid;
            } else if pos >= iv.end {
                lo = mid + 1;
            } else {
                return true;
            }
        }
        false
    }
}

impl Base {
    pub fn to_ascii(self) -> u8 {
        match self {
            Base::A => b'A',
            Base::C => b'C',
            Base::G => b'G',
            Base::T => b'T',
            Base::N => b'N',
        }
    }
    /// Complement (A↔T, C↔G, N→N).
    pub fn complement(self) -> Base {
        match self {
            Base::A => Base::T,
            Base::C => Base::G,
            Base::G => Base::C,
            Base::T => Base::A,
            Base::N => Base::N,
        }
    }
}

/// A parsed reference: contigs plus the FASTA file's sha256.
pub struct Reference {
    /// Contigs in file order.
    pub contigs: Vec<Contig>,
    /// sha256 of the raw FASTA file bytes.
    pub fasta_sha256: [u8; 32],
}

impl Reference {
    /// Look up a contig by name; `None` if absent.
    pub fn contig_index(&self, name: &[u8]) -> Option<u32> {
        self.contigs
            .iter()
            .position(|c| c.name.as_bytes() == name)
            .map(|i| i as u32)
    }
}

/// Parse a FASTA file from disk.
///
/// Errors (with line numbers) on: empty file, empty or duplicate contig names.
/// Whitespace inside sequence lines is stripped; illegal characters become `N`.
pub fn parse_fasta(path: &Path) -> Result<Reference, AlignError> {
    let bytes = fs::read(path).map_err(|source| AlignError::FastaIo {
        path: path.display().to_string(),
        source,
    })?;
    parse_fasta_bytes(&bytes)
}

/// Parse a FASTA file from in-memory bytes; see [`parse_fasta`].
pub fn parse_fasta_bytes(bytes: &[u8]) -> Result<Reference, AlignError> {
    let mut contigs: Vec<Contig> = Vec::new();
    let mut cur_name: Option<String> = None;
    let mut cur_len: u32 = 0;
    let mut cur_packed: Vec<u8> = Vec::new();
    let mut cur_n: Vec<Interval> = Vec::new();
    let mut saw_any = false;

    let finish = |contigs: &mut Vec<Contig>,
                  name: Option<String>,
                  len: u32,
                  mut packed: Vec<u8>,
                  n_intervals: Vec<Interval>| {
        if let Some(name) = name {
            // pad the tail byte even when the final bases are all N
            packed.resize(len.div_ceil(4) as usize, 0);
            contigs.push(Contig {
                name,
                len,
                packed,
                n_intervals,
            });
        }
    };

    for (idx, raw_line) in split_lines(bytes).into_iter().enumerate() {
        let lineno = idx + 1;
        let line: Vec<u8> = raw_line
            .iter()
            .copied()
            .filter(|b| !b.is_ascii_whitespace())
            .collect();
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            saw_any = true;
            finish(
                &mut contigs,
                cur_name.take(),
                cur_len,
                std::mem::take(&mut cur_packed),
                std::mem::take(&mut cur_n),
            );
            cur_len = 0;
            let header = &raw_line[1..];
            let token: Vec<u8> = header
                .iter()
                .copied()
                .take_while(|b| !b.is_ascii_whitespace())
                .collect();
            if token.is_empty() {
                return Err(AlignError::FastaFormat {
                    line: lineno,
                    msg: "empty contig name".to_string(),
                });
            }
            cur_name = Some(String::from_utf8_lossy(&token).into_owned());
        } else {
            if cur_name.is_none() {
                return Err(AlignError::FastaFormat {
                    line: lineno,
                    msg: "sequence before any '>' header".to_string(),
                });
            }
            for &b in &line {
                push_base(b, cur_len, &mut cur_packed, &mut cur_n);
                cur_len += 1;
            }
        }
    }

    if !saw_any {
        return Err(AlignError::FastaFormat {
            line: 0,
            msg: "empty fasta file".to_string(),
        });
    }
    finish(&mut contigs, cur_name, cur_len, cur_packed, cur_n);

    for (i, c) in contigs.iter().enumerate() {
        for d in contigs.iter().skip(i + 1) {
            if c.name == d.name {
                return Err(AlignError::FastaFormat {
                    line: 0,
                    msg: format!("duplicate contig name '{}'", c.name),
                });
            }
        }
    }

    let mut hasher = Sha256::new();
    hasher.update(bytes);
    let fasta_sha256: [u8; 32] = hasher.finalize().into();

    Ok(Reference {
        contigs,
        fasta_sha256,
    })
}

fn push_base(b: u8, pos: u32, packed: &mut Vec<u8>, n_intervals: &mut Vec<Interval>) {
    let base = Base::from_ascii(b);
    match base.code() {
        Some(code) => {
            let byte_idx = (pos >> 2) as usize;
            while byte_idx >= packed.len() {
                packed.push(0);
            }
            let shift = 6 - 2 * (pos & 3);
            packed[byte_idx] |= code << shift;
        }
        None => match n_intervals.last_mut() {
            Some(iv) if iv.end == pos => iv.end += 1,
            _ => n_intervals.push(Interval {
                start: pos,
                end: pos + 1,
            }),
        },
    }
}

fn split_lines(bytes: &[u8]) -> Vec<&[u8]> {
    let mut lines = Vec::new();
    let mut start = 0usize;
    for (i, b) in bytes.iter().enumerate() {
        if *b == b'\n' {
            lines.push(&bytes[start..i]);
            start = i + 1;
        }
    }
    if start < bytes.len() {
        lines.push(&bytes[start..]);
    }
    lines
}
