//! Coordinate sort + BAI index (contract covers map-produced BAMs only:
//! CIGAR subset M/I/D/N/S, tag types Char/Int/Str).
//!
//! Algorithm: sequential read → chunk by record count → stable in-chunk sort
//! by `(ref_id, pos)` (ties keep input order) → temp shards → stable k-way
//! merge. Total order + stable merge ⇒ byte-identical output for the same
//! input regardless of chunk size. The header is copied without semantic
//! edits (the SO line is left untouched). After write-out, a BAI index is
//! built via rust-htslib at `<output>.bai`.

use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::fs::{self, File};
use std::io::{self, BufReader, Write};
use std::path::{Path, PathBuf};

use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::Header;

/// Coordinate-sort options.
#[derive(Clone, Debug)]
pub struct SortOptions {
    /// Maximum records held in memory per chunk (default 2_000_000).
    pub max_in_memory_records: usize,
    /// Temp shard directory; defaults to `<output>.sorttmp/`.
    pub temp_dir: Option<PathBuf>,
    /// BGZF compression threads for the output (clamped to 1..=4, same as
    /// [`crate::create_writer`]); also used for index building.
    pub threads: usize,
}

impl Default for SortOptions {
    fn default() -> Self {
        Self {
            max_in_memory_records: 2_000_000,
            temp_dir: None,
            threads: 1,
        }
    }
}

/// Coordinate-sort statistics.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SortStats {
    /// Total records moved.
    pub records: u64,
    /// Shards produced (1 = whole input fit in memory, no temp files).
    pub chunks: usize,
}

type InputReader = bam::io::Reader<noodles_bgzf::io::Reader<BufReader<File>>>;

/// Sort key: `(ref_id, pos)`; unmapped records (no reference) sort after all
/// mapped records and keep their relative input order via stable merge.
fn sort_key(rec: &RecordBuf) -> (u64, u64) {
    let r = rec.reference_sequence_id().map_or(u64::MAX, |i| i as u64);
    let p = rec
        .alignment_start()
        .map_or(u64::MAX, |pos| usize::from(pos) as u64);
    (r, p)
}

fn open_reader(path: &Path) -> io::Result<InputReader> {
    let file = File::open(path)?;
    Ok(bam::io::Reader::new(BufReader::new(file)))
}

fn write_record<W: Write>(
    writer: &mut bam::io::Writer<W>,
    header: &Header,
    rec: &RecordBuf,
) -> io::Result<()> {
    writer.write_alignment_record(header, rec)
}

/// Sort a chunk in memory and write it as one shard BAM.
fn flush_chunk(
    dir: &Path,
    idx: usize,
    header: &Header,
    chunk: &mut [RecordBuf],
) -> io::Result<PathBuf> {
    chunk.sort_by_key(sort_key);
    let path = dir.join(format!("chunk_{idx:06}"));
    let mut writer = bam::io::Writer::new(File::create(&path)?);
    writer.write_header(header)?;
    for rec in chunk.iter() {
        write_record(&mut writer, header, rec)?;
    }
    writer.try_finish()?;
    Ok(path)
}

struct ChunkReader {
    reader: InputReader,
    header: Header,
    next: Option<RecordBuf>,
}

impl ChunkReader {
    fn open(path: &Path) -> io::Result<Self> {
        let mut reader = open_reader(path)?;
        let header = reader.read_header()?;
        let mut this = Self {
            reader,
            header,
            next: None,
        };
        this.advance()?;
        Ok(this)
    }

    fn advance(&mut self) -> io::Result<()> {
        self.next = self.reader.record_bufs(&self.header).next().transpose()?;
        Ok(())
    }
}

struct Node {
    key: (u64, u64),
    chunk: usize,
    record: RecordBuf,
}

impl PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key && self.chunk == other.chunk
    }
}

impl Eq for Node {}

impl PartialOrd for Node {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

// Inverted: BinaryHeap is a max-heap; the smallest (key, chunk) pops first.
// Equal keys resolve to the lower chunk index, preserving global input order.
impl Ord for Node {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .key
            .cmp(&self.key)
            .then_with(|| other.chunk.cmp(&self.chunk))
    }
}

fn merge_shards(
    paths: &[PathBuf],
    header: &Header,
    output: &Path,
    threads: usize,
) -> io::Result<u64> {
    let mut chunks = Vec::with_capacity(paths.len());
    let mut heap = BinaryHeap::new();
    for (i, path) in paths.iter().enumerate() {
        let mut cr = ChunkReader::open(path)?;
        if let Some(rec) = cr.next.take() {
            heap.push(Node {
                key: sort_key(&rec),
                chunk: i,
                record: rec,
            });
        }
        chunks.push(cr);
    }

    let mut writer = crate::create_writer(File::create(output)?, header, threads)?;
    let mut n = 0u64;
    while let Some(Node { chunk, record, .. }) = heap.pop() {
        write_record(&mut writer, header, &record)?;
        n += 1;
        let cr = &mut chunks[chunk];
        cr.advance()?;
        if let Some(rec) = cr.next.take() {
            heap.push(Node {
                key: sort_key(&rec),
                chunk,
                record: rec,
            });
        }
    }
    writer.into_inner().finish()?;
    Ok(n)
}

fn build_bai(output: &Path, threads: usize) -> io::Result<()> {
    rust_htslib::bam::index::build(
        output,
        None::<&Path>,
        rust_htslib::bam::index::Type::Bai,
        threads.clamp(1, 255) as u32,
    )
    .map_err(io::Error::other)
}

/// Sort a BAM by coordinate and build a BAI index at `<output>.bai`.
///
/// Byte-identical output for identical input, independent of chunk size and
/// of `threads`. Temp shards live in `opts.temp_dir` (default
/// `<output>.sorttmp/`) and are removed on success; residue from a failed run
/// is overwritten by the next run.
pub fn coordinate_sort(input: &Path, output: &Path, opts: &SortOptions) -> io::Result<SortStats> {
    let mut reader = open_reader(input)?;
    let header = reader.read_header()?;
    let cap = opts.max_in_memory_records.max(1);

    // First chunk stays in memory; one look-ahead record decides whether any
    // temp shard is needed at all.
    let mut chunk: Vec<RecordBuf> = Vec::new();
    let mut pending: Option<RecordBuf> = None;
    {
        let mut it = reader.record_bufs(&header);
        for rec in it.by_ref() {
            chunk.push(rec?);
            if chunk.len() >= cap {
                pending = it.next().transpose()?;
                break;
            }
        }
    }

    if pending.is_none() {
        chunk.sort_by_key(sort_key);
        let mut writer = crate::create_writer(File::create(output)?, &header, opts.threads)?;
        let n = chunk.len() as u64;
        for rec in &chunk {
            write_record(&mut writer, &header, rec)?;
        }
        writer.into_inner().finish()?;
        build_bai(output, opts.threads)?;
        return Ok(SortStats {
            records: n,
            chunks: 1,
        });
    }

    let tmp = opts.temp_dir.clone().unwrap_or_else(|| {
        let mut s = output.as_os_str().to_os_string();
        s.push(".sorttmp");
        PathBuf::from(s)
    });
    fs::create_dir_all(&tmp)?;

    let mut paths: Vec<PathBuf> = Vec::new();
    let mut total = 0u64;
    loop {
        if let Some(rec) = pending.take() {
            chunk.push(rec);
        }
        total += chunk.len() as u64;
        paths.push(flush_chunk(&tmp, paths.len(), &header, &mut chunk)?);
        chunk.clear();
        {
            let mut it = reader.record_bufs(&header);
            for rec in it.by_ref() {
                chunk.push(rec?);
                if chunk.len() >= cap {
                    pending = it.next().transpose()?;
                    break;
                }
            }
        }
        if chunk.is_empty() && pending.is_none() {
            break;
        }
    }

    merge_shards(&paths, &header, output, opts.threads)?;
    fs::remove_dir_all(&tmp)?;
    build_bai(output, opts.threads)?;

    Ok(SortStats {
        records: total,
        chunks: paths.len(),
    })
}

