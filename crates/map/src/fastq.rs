//! FASTQ / BFQ record sources.
//!
//! `RecordSource` unifies plain and gzipped FASTQ (streaming) with the
//! memory-mapped EBFQ format. All FASTQ format errors carry 1-based line
//! numbers.

use crate::error::AlignError;
use flate2::read::MultiGzDecoder;
use memmap2::Mmap;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

/// One sequencing read: name (first header token), sequence, qualities.
#[derive(Debug)]
pub struct FastqRecord {
    /// Read name, first whitespace-delimited token of the `@` header.
    pub name: Vec<u8>,
    /// Sequence bases (ASCII).
    pub seq: Vec<u8>,
    /// Quality values (ASCII, phred+offset as stored).
    pub qual: Vec<u8>,
}

/// Source of reads, consumable one record at a time.
pub trait RecordSource {
    /// Next record, or `None` at end of stream. Malformed input errors.
    fn next_record(&mut self) -> Result<Option<FastqRecord>, AlignError>;
}

/// Streaming reader for plain and gzipped FASTQ.
pub struct FastqReader {
    inner: Box<dyn Read + Send>,
    buf: Vec<u8>,
    eof: bool,
    line_no: usize,
}

impl FastqReader {
    /// Open a FASTQ file; a `.gz` extension routes through `MultiGzDecoder`
    /// (1 MiB buffer), anything else is read as plain text (64 KiB buffer).
    pub fn open(path: &Path) -> Result<Self, AlignError> {
        let file = File::open(path).map_err(|e| AlignError::FastqIo {
            path: path.display().to_string(),
            source: e,
        })?;
        let inner: Box<dyn Read + Send> = if path.extension().is_some_and(|e| e == "gz") {
            Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::with_capacity(1 << 16, file))
        };
        Ok(FastqReader {
            inner,
            buf: Vec::new(),
            eof: false,
            line_no: 0,
        })
    }

    /// Fill `buf` with the next `\n`-terminated chunk; strips `\r`/`\n`.
    /// Returns `false` at clean EOF with no bytes read.
    fn read_line(&mut self) -> Result<bool, AlignError> {
        self.buf.clear();
        let mut byte = [0u8; 1];
        loop {
            match self.inner.read(&mut byte) {
                Ok(0) => {
                    self.eof = true;
                    break;
                }
                Ok(_) => {
                    if byte[0] == b'\n' {
                        break;
                    }
                    self.buf.push(byte[0]);
                }
                Err(e) => return Err(AlignError::FastqIo {
                    path: "<fastq stream>".to_string(),
                    source: e,
                }),
            }
        }
        if self.buf.last() == Some(&b'\r') {
            self.buf.pop();
        }
        self.line_no += 1;
        Ok(!self.buf.is_empty() || !self.eof)
    }
}

impl RecordSource for FastqReader {
    fn next_record(&mut self) -> Result<Option<FastqRecord>, AlignError> {
        // Header.
        let header_line = self.line_no + 1;
        if !self.read_line()? {
            return Ok(None);
        }
        if self.buf.is_empty() || self.buf[0] != b'@' {
            return Err(AlignError::FastqFormat {
                line: header_line,
                msg: "expected '@' header line".to_string(),
            });
        }
        let name: Vec<u8> = self.buf[1..]
            .iter()
            .copied()
            .take_while(|b| !b.is_ascii_whitespace())
            .collect();
        if name.is_empty() {
            return Err(AlignError::FastqFormat {
                line: header_line,
                msg: "empty read name".to_string(),
            });
        }

        // Sequence.
        let seq_line = self.line_no + 1;
        if !self.read_line()? {
            return Err(AlignError::FastqFormat {
                line: seq_line,
                msg: "truncated record: missing sequence".to_string(),
            });
        }
        let seq = self.buf.clone();

        // Separator.
        let plus_line = self.line_no + 1;
        if !self.read_line()? {
            return Err(AlignError::FastqFormat {
                line: plus_line,
                msg: "truncated record: missing '+' separator".to_string(),
            });
        }
        if self.buf.first() != Some(&b'+') {
            return Err(AlignError::FastqFormat {
                line: plus_line,
                msg: "expected '+' separator line".to_string(),
            });
        }

        // Quality.
        let qual_line = self.line_no + 1;
        if !self.read_line()? {
            return Err(AlignError::FastqFormat {
                line: qual_line,
                msg: "truncated record: missing quality".to_string(),
            });
        }
        let qual = self.buf.clone();
        if seq.len() != qual.len() {
            return Err(AlignError::FastqFormat {
                line: qual_line,
                msg: format!(
                    "sequence length {} != quality length {}",
                    seq.len(),
                    qual.len()
                ),
            });
        }

        Ok(Some(FastqRecord { name, seq, qual }))
    }
}

/// Memory-mapped EBFQ reader.
///
/// Layout: magic `EBFQ` at [0..4], [4..8] reserved, [8..16] read_count as
/// u64 LE; each record is u16 name_len + name + u16 seq_len + seq + qual.
pub struct BfqReader {
    mmap: Mmap,
    read_count: u64,
    seen: u64,
    pos: usize,
}

impl BfqReader {
    /// Map `path` and validate the EBFQ header.
    pub fn open(path: &Path) -> Result<Self, AlignError> {
        let file = File::open(path).map_err(|e| AlignError::FastqIo {
            path: path.display().to_string(),
            source: e,
        })?;
        let mmap = unsafe { Mmap::map(&file) }.map_err(|e| AlignError::FastqIo {
            path: path.display().to_string(),
            source: e,
        })?;
        if mmap.len() < 16 || &mmap[0..4] != b"EBFQ" {
            return Err(AlignError::FastqFormat {
                line: 0,
                msg: "bad EBFQ magic".to_string(),
            });
        }
        let read_count = u64::from_le_bytes(mmap[8..16].try_into().unwrap());
        Ok(BfqReader {
            mmap,
            read_count,
            seen: 0,
            pos: 16,
        })
    }
}

impl RecordSource for BfqReader {
    fn next_record(&mut self) -> Result<Option<FastqRecord>, AlignError> {
        if self.seen >= self.read_count {
            return Ok(None);
        }
        let data = &self.mmap[self.pos..];
        if data.len() < 2 {
            return Err(AlignError::FastqFormat {
                line: 0,
                msg: "truncated EBFQ record".to_string(),
            });
        }
        let name_len = u16::from_le_bytes([data[0], data[1]]) as usize;
        if data.len() < 2 + name_len + 2 {
            return Err(AlignError::FastqFormat {
                line: 0,
                msg: "truncated EBFQ record".to_string(),
            });
        }
        let raw_name = &data[2..2 + name_len];
        let name: Vec<u8> = raw_name
            .iter()
            .copied()
            .take_while(|b| !b.is_ascii_whitespace())
            .collect();
        let off = 2 + name_len;
        let seq_len = u16::from_le_bytes([data[off], data[off + 1]]) as usize;
        if data.len() < off + 2 + seq_len * 2 {
            return Err(AlignError::FastqFormat {
                line: 0,
                msg: "truncated EBFQ record".to_string(),
            });
        }
        let seq = data[off + 2..off + 2 + seq_len].to_vec();
        let qual = data[off + 2 + seq_len..off + 2 + seq_len * 2].to_vec();
        self.pos += off + 2 + seq_len * 2;
        self.seen += 1;
        Ok(Some(FastqRecord { name, seq, qual }))
    }
}
