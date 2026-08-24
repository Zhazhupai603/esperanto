//! esperanto-tidx — transcriptome k-mer index (`.tidx`).
//!
//! Crate responsibility
//! --------------------
//! Build and query a k-mer → (transcript, offset) index over a GENCODE-style
//! GTF annotation and a reference FASTA. The builder extracts every ACGT-only
//! k-mer window of every transcript (exons concatenated in genomic ascending
//! order, reverse-complemented whole on the minus strand), maps each window
//! to its canonical form `min(fwd, revcomp)`, and writes a bucketed, sorted
//! CSR file. The reader serves random lookups of a canonical k-mer as a slice
//! of `(tx_id, offset)` pairs in ascending `(tx_id, offset)` order.
//!
//! Public API
//! ----------
//! * [`build`] — build an index from `GTF + FASTA(.fai)`; deterministic:
//!   identical inputs and an identical [`BuildOptions::timestamp`] produce
//!   byte-identical `.tidx` files.
//! * [`BuildOptions`] — `k = 31`, `prefix_bits = 20`, `timestamp = now()`,
//!   `threads = 16`; every field overridable.
//! * [`BuildStats`] — build outcome summary (transcript count, indexed bp,
//!   entry count, distinct canonical k-mers, file size, wall time, reference
//!   BLAKE3, ...).
//! * [`Tidx`] — in-memory reader: `open`, `k`, `tx_count`, `lookup`,
//!   `transcript_name`, `verify_reference`.
//! * [`Error`] — `Io`, `BadMagic`, `UnsupportedVersion`, `Truncated`,
//!   `Inconsistent`, `BadGtf`, `UnknownContig`, `RefMismatch`.
//!
//! k-mer encoding
//! --------------
//! `A=00 C=01 G=10 T=11` (case-insensitive); a k-mer occupies the low `2k`
//! bits of a `u64` with the 5'-most base in the highest bits of the window.
//! `revcomp` complements each code (`3 - code`) and reverses their order;
//! `canonical = min(forward, revcomp)`. Any non-ACGT byte (N / IUPAC)
//! invalidates every window overlapping it. The builder maintains `fwd` and
//! `rc` incrementally: on an entering base `c`,
//! `fwd' = ((fwd << 2) | c) & mask` and `rc' = (rc >> 2) | ((3 - c) << (2k - 2))`
//! where `mask = (1 << 2k) - 1` (`u64::MAX` when `k = 32`).
//!
//! File format (all little-endian)
//! -------------------------------
//! ```text
//! offset  size        content
//! 0       128         header
//! 128     B           name blob (transcript names, UTF-8, concatenated)
//!         (N+1)*4     name offsets (u32; name i = blob[off[i]..off[i+1]])
//!         (2^P+1)*8   bucket offsets (u64; exclusive CSR suffix)
//!         E*w         keys (w = 4 for version 2, 8 for version 1)
//!         E*8         payloads (u32 tx_id || u32 offset)
//! ```
//! Header fields: magic `"TIDX"`, version (2 when `2k - P <= 32`, else 1),
//! `k` (15..=32), `m = 16`, `prefix_bits P`, reserved 0, fixed seed
//! `0x7469647873656564`, 32-byte BLAKE3 of the raw reference FASTA file,
//! `tx_count`, reserved 0, unix-seconds timestamp, `total_entries E`,
//! `bucket_count = 1 << P`, `name_blob_len`, `name_offsets_len = tx_count + 1`,
//! reserved 0.
//!
//! Records are globally sorted by `(canonical, tx_id, offset)`; the bucket of
//! a canonical k-mer is `canonical >> (2k - P)` and the stored key is its low
//! `2k - P` bits. Lookup = bucket locate → binary search on the low bits →
//! contiguous payload slice. The reader supports both versions.
//!
//! GTF semantics
//! -------------
//! Only `transcript` and `exon` rows are consumed. Coordinates are 1-based
//! inclusive in the file and stored 0-based half-open `[start-1, end)`. Rows
//! aggregate by the `transcript_id` attribute; a `transcript` row refreshes
//! name (`transcript_name`, falling back to `transcript_id`), contig and
//! strand; an `exon` row only adds a span (exon rows appearing before the
//! `transcript` row are handled). Transcripts without exons are dropped;
//! exons are sorted by `(start, end)`. Bad rows (fewer than 9 columns,
//! missing `transcript_id`, strand not `+`/`-`, `start == 0`, `end < start`)
//! fail with the 1-based line number.
//!
//! Invariants
//! ----------
//! * Determinism: identical inputs + pinned timestamp → byte-identical file.
//!   All output-affecting iteration is over sorted vectors; rayon only fills
//!   per-transcript buffers that are flattened in transcript order, and the
//!   final record sort is a total order.
//! * `tx_id`s are dense `0..N-1` assigned in `transcript_id` lexicographic
//!   order over transcripts that survived GTF filtering. Transcripts dropped
//!   later (shorter than k, exon beyond contig) keep their id slot.
//! * Lookup results are sorted by `(tx_id, offset)` as a consequence of the
//!   global record order; misses return an empty slice.
//! * `Tidx::open` returns `io::Result`, wrapping crate errors as
//!   `io::ErrorKind::InvalidData`.

#![deny(unsafe_code)]

use std::collections::{BTreeMap, HashMap};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use rayon::prelude::*;
use thiserror::Error;

/// File magic.
const MAGIC: [u8; 4] = *b"TIDX";
/// Header size in bytes.
const HEADER_LEN: usize = 128;
/// Reserved `m` header field (fixed value 16).
const M_FIELD: u32 = 16;
/// Fixed hash seed constant (ASCII `"tidxseed"`).
const SEED: u64 = 0x7469647873656564;

/// Smallest legal k-mer length.
pub const MIN_K: u32 = 15;
/// Largest legal k-mer length (2-bit codes fill a u64 exactly).
pub const MAX_K: u32 = 32;
/// Default k-mer length.
pub const DEFAULT_K: u32 = 31;
/// Default bucket prefix bits.
pub const DEFAULT_PREFIX_BITS: u32 = 20;
/// Default worker thread count for the build.
pub const DEFAULT_THREADS: usize = 16;

/// Crate result alias.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors produced by this crate.
#[derive(Debug, Error)]
pub enum Error {
    /// Underlying I/O failure.
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),

    /// File does not start with the `"TIDX"` magic.
    #[error("bad magic {magic:?} (version field {version})")]
    BadMagic {
        /// The four magic bytes found in the file.
        magic: [u8; 4],
        /// The version field as read (may be garbage).
        version: u32,
    },

    /// Header version is neither 1 nor 2.
    #[error("unsupported .tidx version: {0}")]
    UnsupportedVersion(u32),

    /// A declared segment does not fit the file (wrong segment lengths).
    #[error("truncated or malformed .tidx: {0}")]
    Truncated(String),

    /// Header/segment metadata contradicts itself (k out of range, mismatched
    /// counts, non-monotone offsets, ...).
    #[error("inconsistent .tidx metadata: {0}")]
    Inconsistent(String),

    /// Malformed GTF row; `line` is 1-based.
    #[error("bad GTF line {line}: {reason}")]
    BadGtf {
        /// 1-based line number in the GTF file.
        line: usize,
        /// What is wrong with the row.
        reason: String,
    },

    /// A GTF contig is absent from the reference `.fai`.
    #[error("contig not present in reference .fai: {0}")]
    UnknownContig(String),

    /// `verify_reference` found a BLAKE3 mismatch.
    #[error("reference FASTA does not match the index (expected blake3 {expected:02x?}, got {actual:02x?})")]
    RefMismatch {
        /// BLAKE3 recorded in the index.
        expected: [u8; 32],
        /// BLAKE3 of the FASTA passed to `verify_reference`.
        actual: [u8; 32],
    },
}

/// Options for [`build`]. All fields overridable; see the `Default` values.
#[derive(Debug, Clone)]
pub struct BuildOptions {
    /// k-mer length, legal range 15..=32.
    pub k: u32,
    /// Bucket prefix bits P, legal range 1..=min(2k, 63).
    pub prefix_bits: u32,
    /// Unix-seconds timestamp written into the header. Pin it for
    /// byte-reproducible builds; the default is the current wall clock.
    pub timestamp: u64,
    /// Worker threads for per-transcript k-mer extraction (passed to a local
    /// rayon pool; 0 lets rayon choose).
    pub threads: usize,
}

impl Default for BuildOptions {
    fn default() -> Self {
        Self {
            k: DEFAULT_K,
            prefix_bits: DEFAULT_PREFIX_BITS,
            timestamp: unix_now(),
            threads: DEFAULT_THREADS,
        }
    }
}

/// Summary of a completed [`build`].
#[derive(Debug, Clone)]
pub struct BuildStats {
    /// Number of transcripts that survived GTF filtering (dense id space).
    pub tx_count: u32,
    /// Total transcript bases indexed (sum of retained transcript lengths).
    pub total_bp: u64,
    /// Total (canonical, tx_id, offset) records written.
    pub total_entries: u64,
    /// Number of distinct canonical k-mers.
    pub distinct_canonical: u64,
    /// Bucket count (`1 << prefix_bits`).
    pub bucket_count: u64,
    /// Size of the written `.tidx` file in bytes.
    pub file_size: u64,
    /// Wall-clock build time in seconds (statistic only; never written).
    pub build_seconds: f64,
    /// Prefix bits used.
    pub prefix_bits: u32,
    /// BLAKE3 of the raw reference FASTA file bytes.
    pub ref_blake3: [u8; 32],
}

/// In-memory reader for a `.tidx` index. The whole file is parsed at
/// [`Tidx::open`] time; lookups are pure in-memory binary searches.
#[derive(Debug, Clone)]
pub struct Tidx {
    k: u32,
    tx_count: u32,
    ref_blake3: [u8; 32],
    names: Vec<String>,
    bucket_offsets: Vec<u64>,
    /// Stored key low bits, widened to u64.
    keys: Vec<u64>,
    payloads: Vec<(u32, u32)>,
    key_mask: u64,
    bucket_shift: u32,
}

impl Tidx {
    /// Open and fully parse a `.tidx` file.
    ///
    /// Crate-level failures (bad magic, unsupported version, truncated or
    /// inconsistent segments) are reported as `io::ErrorKind::InvalidData`.
    pub fn open(path: &Path) -> io::Result<Self> {
        let data = fs::read(path)?;
        parse_index(&data).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    /// k-mer length the index was built with.
    pub fn k(&self) -> u32 {
        self.k
    }

    /// Number of transcripts in the dense id space.
    pub fn tx_count(&self) -> u32 {
        self.tx_count
    }

    /// Look up a canonical k-mer; returns its `(tx_id, offset)` pairs in
    /// ascending `(tx_id, offset)` order, or an empty slice on miss.
    /// `offset` is the 0-based transcript coordinate of the window start.
    pub fn lookup(&self, canonical_kmer: u64) -> &[(u32, u32)] {
        let bucket = (canonical_kmer >> self.bucket_shift) as usize;
        if bucket + 1 >= self.bucket_offsets.len() {
            return &[];
        }
        let s = self.bucket_offsets[bucket] as usize;
        let e = self.bucket_offsets[bucket + 1] as usize;
        let low = canonical_kmer & self.key_mask;
        let range = &self.keys[s..e];
        let start = range.partition_point(|&key| key < low);
        let end = range.partition_point(|&key| key <= low);
        &self.payloads[s + start..s + end]
    }

    /// Transcript display name for `tx_id`; empty string when out of range.
    pub fn transcript_name(&self, tx_id: u32) -> &str {
        self.names
            .get(tx_id as usize)
            .map(String::as_str)
            .unwrap_or("")
    }

    /// Verify that `fasta` hashes (BLAKE3) to the reference digest recorded
    /// in the index; mismatches fail with [`Error::RefMismatch`].
    pub fn verify_reference(&self, fasta: &Path) -> Result<()> {
        let actual = blake3_file(fasta)?;
        if actual == self.ref_blake3 {
            Ok(())
        } else {
            Err(Error::RefMismatch {
                expected: self.ref_blake3,
                actual,
            })
        }
    }
}

/// Build a `.tidx` index from a GTF annotation plus an indexed reference
/// FASTA (a samtools-style `.fai` must exist next to `ref_path`).
///
/// Deterministic: identical inputs with an identical `opts.timestamp`
/// produce byte-identical output files.
pub fn build(
    gtf_path: &Path,
    ref_path: &Path,
    out_path: &Path,
    opts: &BuildOptions,
) -> Result<BuildStats> {
    let t0 = Instant::now();

    // ---- option validation -------------------------------------------------
    let k = opts.k;
    if !(MIN_K..=MAX_K).contains(&k) {
        return Err(Error::Inconsistent(format!(
            "k = {k} outside legal range {MIN_K}..={MAX_K}"
        )));
    }
    let prefix_bits = opts.prefix_bits;
    let max_p = (2 * k).min(63);
    if prefix_bits == 0 || prefix_bits > max_p {
        return Err(Error::Inconsistent(format!(
            "prefix_bits = {prefix_bits} outside legal range 1..={max_p}"
        )));
    }
    let low_bits = 2 * k - prefix_bits; // in 1..=63
    let version = if low_bits <= 32 { 2u32 } else { 1u32 };
    let key_mask: u64 = (1u64 << low_bits) - 1;
    let bucket_count: u64 = 1u64 << prefix_bits;
    let bucket_table_bytes = bucket_count
        .checked_add(1)
        .and_then(|n| n.checked_mul(8))
        .filter(|n| *n <= usize::MAX as u64)
        .ok_or_else(|| {
            Error::Inconsistent(format!(
                "bucket table for prefix_bits = {prefix_bits} exceeds addressable memory"
            ))
        })?;

    // ---- 1. GTF: aggregate, filter, dense ids in transcript_id order -------
    let recs = parse_gtf(gtf_path)?;
    if recs.len() > u32::MAX as usize {
        return Err(Error::Inconsistent(format!(
            "{} transcripts exceed the u32 tx_id space",
            recs.len()
        )));
    }
    let tx_count = recs.len() as u32;

    // ---- 2. .fai: contig order, UnknownContig check, grouping --------------
    let fai = read_fai(ref_path)?;
    let mut fai_index: HashMap<&str, usize> = HashMap::with_capacity(fai.len());
    for (i, rec) in fai.iter().enumerate() {
        fai_index.insert(rec.name.as_str(), i);
    }
    let mut by_fai: Vec<Vec<u32>> = vec![Vec::new(); fai.len()];
    for (tx_id, rec) in recs.iter().enumerate() {
        let &fi = fai_index
            .get(rec.contig.as_str())
            .ok_or_else(|| Error::UnknownContig(rec.contig.clone()))?;
        by_fai[fi].push(tx_id as u32);
    }

    // ---- 3. per-contig fetch + parallel per-transcript extraction ----------
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(opts.threads)
        .build()
        .map_err(|e| io::Error::other(format!("failed to create thread pool: {e}")))?;
    let mut entries: Vec<(u64, u32, u32)> = Vec::new();
    // Reference semantics: total_bp counts the spliced length (sum of exon
    // spans) of EVERY exon-bearing transcript, including ones later dropped
    // for being shorter than k or running off the contig.
    let total_bp: u64 = recs
        .iter()
        .map(|r| r.exons.iter().map(|&(s, e)| e - s).sum::<u64>())
        .sum();
    let mut fasta = File::open(ref_path)?;
    for (fi, frec) in fai.iter().enumerate() {
        if by_fai[fi].is_empty() {
            continue;
        }
        let contig = fetch_contig(&mut fasta, frec)?;
        let results: Vec<(KmerEntries, u64)> = pool.install(|| {
            by_fai[fi]
                .par_iter()
                .map(|&tx_id| {
                    extract_for_transcript(&recs[tx_id as usize], &contig, frec.len, k, tx_id)
                })
                .collect()
        });
        for (mut transcript_entries, _bp) in results {
            entries.append(&mut transcript_entries);
        }
    }

    // ---- 4. global sort by (canonical, tx_id, offset) ----------------------
    entries.sort_unstable();
    let total_entries = entries.len() as u64;

    // ---- 5. reference BLAKE3 (streamed over the raw FASTA file) ------------
    let ref_blake3 = blake3_file(ref_path)?;

    // ---- 6. CSR bucket offsets ---------------------------------------------
    let mut bucket_offsets: Vec<u64> = Vec::with_capacity(bucket_table_bytes as usize / 8);
    bucket_offsets.resize(bucket_count as usize + 1, 0);
    for &(canonical, _, _) in &entries {
        bucket_offsets[(canonical >> low_bits) as usize + 1] += 1;
    }
    for i in 1..bucket_offsets.len() {
        bucket_offsets[i] += bucket_offsets[i - 1];
    }

    // ---- distinct canonical k-mers ------------------------------------------
    let mut distinct_canonical: u64 = 0;
    let mut prev: Option<u64> = None;
    for &(canonical, _, _) in &entries {
        if prev != Some(canonical) {
            distinct_canonical += 1;
            prev = Some(canonical);
        }
    }

    // ---- 7. name blob -------------------------------------------------------
    let mut blob: Vec<u8> = Vec::new();
    let mut name_offsets: Vec<u32> = Vec::with_capacity(tx_count as usize + 1);
    name_offsets.push(0);
    for rec in &recs {
        blob.extend_from_slice(rec.name.as_bytes());
        let len32 = u32::try_from(blob.len())
            .map_err(|_| Error::Inconsistent("name blob exceeds u32 length space".to_string()))?;
        name_offsets.push(len32);
    }

    // ---- 8. stream write ----------------------------------------------------
    let out = File::create(out_path)?;
    let mut w = BufWriter::with_capacity(1 << 20, out);
    let mut header = [0u8; HEADER_LEN];
    header[0..4].copy_from_slice(&MAGIC);
    header[4..8].copy_from_slice(&version.to_le_bytes());
    header[8..12].copy_from_slice(&k.to_le_bytes());
    header[12..16].copy_from_slice(&M_FIELD.to_le_bytes());
    header[16..20].copy_from_slice(&prefix_bits.to_le_bytes());
    // 20..24 reserved = 0
    header[24..32].copy_from_slice(&SEED.to_le_bytes());
    header[32..64].copy_from_slice(&ref_blake3);
    header[64..68].copy_from_slice(&tx_count.to_le_bytes());
    // 68..72 reserved = 0
    header[72..80].copy_from_slice(&opts.timestamp.to_le_bytes());
    header[80..88].copy_from_slice(&total_entries.to_le_bytes());
    header[88..96].copy_from_slice(&bucket_count.to_le_bytes());
    header[96..104].copy_from_slice(&(blob.len() as u64).to_le_bytes());
    header[104..112].copy_from_slice(&(name_offsets.len() as u64).to_le_bytes());
    // 112..128 reserved = 0
    w.write_all(&header)?;
    w.write_all(&blob)?;
    for off in &name_offsets {
        w.write_all(&off.to_le_bytes())?;
    }
    for off in &bucket_offsets {
        w.write_all(&off.to_le_bytes())?;
    }
    if version == 2 {
        for &(canonical, _, _) in &entries {
            w.write_all(&((canonical & key_mask) as u32).to_le_bytes())?;
        }
    } else {
        for &(canonical, _, _) in &entries {
            w.write_all(&(canonical & key_mask).to_le_bytes())?;
        }
    }
    for &(_, tx_id, offset) in &entries {
        w.write_all(&tx_id.to_le_bytes())?;
        w.write_all(&offset.to_le_bytes())?;
    }
    w.flush()?;

    let file_size = fs::metadata(out_path)?.len();

    Ok(BuildStats {
        tx_count,
        total_bp,
        total_entries,
        distinct_canonical,
        bucket_count,
        file_size,
        build_seconds: t0.elapsed().as_secs_f64(),
        prefix_bits,
        ref_blake3,
    })
}

// ---------------------------------------------------------------------------
// k-mer encoding
// ---------------------------------------------------------------------------

/// One transcript's extracted records: `(canonical, tx_id, offset)`.
type KmerEntries = Vec<(u64, u32, u32)>;

/// 2-bit code of one base byte, case-insensitive; `None` for any non-ACGT.
#[inline]
fn base_code(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Reverse complement of a byte sequence (non-ACGT bytes pass through; they
/// invalidate windows either way).
/// Case-preserving reverse complement : a<->t, c<->g;
/// non-ACGT bytes pass through unchanged. Zero effect on index bytes — k-mer
/// encoding is case-insensitive — this is internal consistency only.
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'a' => b't',
            b'C' => b'G',
            b'c' => b'g',
            b'G' => b'C',
            b'g' => b'c',
            b'T' => b'A',
            b't' => b'a',
            other => other,
        })
        .collect()
}

/// Slide a k-window over `seq`, emitting `(canonical, tx_id, offset)` for
/// every ACGT-only window. Rolling updates keep `fwd` and the reverse
/// complement `rc` in lockstep; a non-ACGT byte invalidates the window.
fn extract_kmers(seq: &[u8], k: u32, tx_id: u32, out: &mut KmerEntries) {
    let kn = k as usize;
    if seq.len() < kn {
        return;
    }
    let mask: u64 = if kn == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * kn)) - 1
    };
    let rc_top = 2 * kn - 2;
    let mut fwd: u64 = 0;
    let mut rc: u64 = 0;
    let mut valid: usize = 0;
    for (i, &b) in seq.iter().enumerate() {
        match base_code(b) {
            Some(c) => {
                fwd = ((fwd << 2) | c) & mask;
                rc = (rc >> 2) | ((3 - c) << rc_top);
                valid += 1;
            }
            None => {
                fwd = 0;
                rc = 0;
                valid = 0;
            }
        }
        if valid >= kn {
            out.push((fwd.min(rc), tx_id, (i + 1 - kn) as u32));
        }
    }
}

/// Build the 5'→3' transcript sequence and extract its k-mers. Returns the
/// entries plus the retained sequence length; dropped transcripts (any exon
/// beyond the contig, sequence shorter than k, sequence too long for u32
/// offsets) return no entries.
fn extract_for_transcript(
    rec: &TxRecord,
    contig: &[u8],
    contig_len: u64,
    k: u32,
    tx_id: u32,
) -> (KmerEntries, u64) {
    if rec.exons.iter().any(|&(_, end)| end > contig_len) {
        return (Vec::new(), 0);
    }
    let mut seq: Vec<u8> = Vec::new();
    for &(start, end) in &rec.exons {
        seq.extend_from_slice(&contig[start as usize..end as usize]);
    }
    if seq.len() < k as usize || seq.len() > u32::MAX as usize {
        return (Vec::new(), 0);
    }
    let seq = if rec.strand == b'-' {
        revcomp(&seq)
    } else {
        seq
    };
    let mut out = Vec::with_capacity(seq.len() - k as usize + 1);
    extract_kmers(&seq, k, tx_id, &mut out);
    (out, seq.len() as u64)
}

// ---------------------------------------------------------------------------
// GTF parsing
// ---------------------------------------------------------------------------

/// Aggregation state while scanning the GTF.
struct TxBuilder {
    name: Option<String>,
    contig: Option<String>,
    strand: Option<u8>,
    exons: Vec<(u64, u64)>,
}

/// A finalized transcript record, exons 0-based half-open and sorted.
/// A finalized transcript record: id, display name, contig, strand byte, exons.
#[derive(Debug, Clone)]
struct TxRecord {
    id: String,
    name: String,
    contig: String,
    strand: u8,
    exons: Vec<(u64, u64)>,
}

/// Parse a GTF: consume `transcript`/`exon` rows, aggregate by
/// `transcript_id`, drop transcripts without exons (or without a
/// `transcript` row supplying contig/strand), sort exons, and return records
/// in `transcript_id` lexicographic order.
fn parse_gtf(path: &Path) -> Result<Vec<TxRecord>> {
    let file = File::open(path)?;
    let mut map: BTreeMap<String, TxBuilder> = BTreeMap::new();
    for (idx, line) in BufReader::new(file).lines().enumerate() {
        let line = line?;
        let lineno = idx + 1;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            return Err(Error::BadGtf {
                line: lineno,
                reason: format!(
                    "expected at least 9 tab-separated columns, found {}",
                    cols.len()
                ),
            });
        }
        match cols[2] {
            "transcript" | "exon" => {}
            _ => continue,
        }
        let attrs = parse_attrs(cols[8]);
        let tx_id = match attrs.get("transcript_id") {
            Some(v) if !v.is_empty() => *v,
            _ => {
                return Err(Error::BadGtf {
                    line: lineno,
                    reason: "missing transcript_id attribute".to_string(),
                })
            }
        };
        let strand = match cols[6] {
            "+" => b'+',
            "-" => b'-',
            other => {
                return Err(Error::BadGtf {
                    line: lineno,
                    reason: format!("strand must be + or -, found {other:?}"),
                })
            }
        };
        let start: u64 = cols[3].parse().map_err(|_| Error::BadGtf {
            line: lineno,
            reason: format!("invalid start coordinate {:?}", cols[3]),
        })?;
        let end: u64 = cols[4].parse().map_err(|_| Error::BadGtf {
            line: lineno,
            reason: format!("invalid end coordinate {:?}", cols[4]),
        })?;
        if start == 0 {
            return Err(Error::BadGtf {
                line: lineno,
                reason: "start coordinate is 0 (GTF is 1-based)".to_string(),
            });
        }
        if end < start {
            return Err(Error::BadGtf {
                line: lineno,
                reason: format!("end {end} < start {start}"),
            });
        }
        let entry = map.entry(tx_id.to_string()).or_insert_with(|| TxBuilder {
            name: None,
            contig: None,
            strand: None,
            exons: Vec::new(),
        });
        if cols[2] == "exon" && entry.contig.is_none() {
            // Reference semantics: an exon row supplies contig/strand when no
            // transcript row has been seen yet (exon-only transcripts).
            entry.contig = Some(cols[0].to_string());
            entry.strand = Some(strand);
        }
        if cols[2] == "transcript" {
            // A transcript row refreshes name/contig/strand.
            entry.name = Some(
                attrs
                    .get("transcript_name")
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| tx_id.to_string()),
            );
            entry.contig = Some(cols[0].to_string());
            entry.strand = Some(strand);
        } else {
            // An exon row only adds a span (converted to 0-based half-open).
            entry.exons.push((start - 1, end));
        }
    }
    Ok(map
        .into_iter()
        .filter_map(|(transcript_id, b)| {
            if b.exons.is_empty() {
                return None;
            }
            // Reference semantics: an exon row may supply contig/strand when
            // the transcript row never appears (exon-only transcripts are
            // indexed, with name falling back to transcript_id).
            let contig = b.contig.clone()?;
            let strand = b.strand?;
            let mut exons = b.exons;
            exons.sort_unstable();
            let name = b.name.unwrap_or_else(|| transcript_id.clone());
            Some(TxRecord {
                id: transcript_id.clone(),
                name,
                contig,
                strand,
                exons,
            })
        })
        .collect())
}

/// Parse a GTF attribute column into `key → value` pairs, tolerating both
/// `key "value";` and `key=value;` styles; quotes are stripped.
/// Reference attribute parser (ported byte-for-byte in behavior): pairs are
/// `key "value";` — the value MUST be double-quoted; unquoted pairs are
/// skipped entirely (a file whose transcript_id is unquoted therefore fails
/// with missing-transcript_id, exactly like the oracle). Empty quoted values
/// are legal. `key=value` (GFF style) is not recognized.
fn parse_attrs(field: &str) -> HashMap<&str, &str> {
    let b = field.as_bytes();
    let mut out = HashMap::new();
    let mut i = 0usize;
    while i < b.len() {
        // skip leading spaces / semicolons
        while i < b.len() && (b[i] == b' ' || b[i] == b';') {
            i += 1;
        }
        // read key up to the first space
        let kstart = i;
        while i < b.len() && b[i] != b' ' {
            i += 1;
        }
        let key = &field[kstart..i];
        // skip one space
        if i < b.len() && b[i] == b' ' {
            i += 1;
        }
        // expect an opening quote; otherwise skip to the next semicolon
        if i >= b.len() || b[i] != b'"' {
            while i < b.len() && b[i] != b';' {
                i += 1;
            }
            continue;
        }
        i += 1; // consume the opening quote
        let vstart = i;
        while i < b.len() && b[i] != b'"' {
            i += 1;
        }
        let value = &field[vstart..i];
        if i < b.len() {
            i += 1; // consume the closing quote
        }
        if !key.is_empty() {
            out.insert(key, value);
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Reference FASTA (.fai) access
// ---------------------------------------------------------------------------

/// One `.fai` row.
struct FaiRec {
    name: String,
    len: u64,
    offset: u64,
    linebases: u64,
    linewidth: u64,
}

/// Read the samtools-style `.fai` beside `ref_path` (i.e. `ref_path` with
/// `.fai` appended), preserving file order.
fn read_fai(ref_path: &Path) -> Result<Vec<FaiRec>> {
    let fai_path = PathBuf::from(format!("{}.fai", ref_path.display()));
    let file = File::open(&fai_path)?;
    let mut recs = Vec::new();
    for (idx, line) in BufReader::new(file).lines().enumerate() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let lineno = idx + 1;
        let cols: Vec<&str> = line.split('\t').collect();
        let bad = |what: String| {
            Error::Io(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(".fai line {lineno}: {what}"),
            ))
        };
        if cols.len() < 5 {
            return Err(bad(format!(
                "expected at least 5 tab-separated columns, found {}",
                cols.len()
            )));
        }
        let num = |s: &str| -> Result<u64> {
            s.parse::<u64>()
                .map_err(|_| bad(format!("malformed integer field {s:?}")))
        };
        recs.push(FaiRec {
            name: cols[0].to_string(),
            len: num(cols[1])?,
            offset: num(cols[2])?,
            linebases: num(cols[3])?,
            linewidth: num(cols[4])?,
        });
    }
    Ok(recs)
}

/// Fetch one contig's bases by seeking into the FASTA per the `.fai` record
/// and stripping line terminators. Returns exactly `rec.len` bytes or an
/// I/O error (truncated FASTA / `.fai` mismatch).
fn fetch_contig(file: &mut File, rec: &FaiRec) -> io::Result<Vec<u8>> {
    if rec.len == 0 {
        return Ok(Vec::new());
    }
    // Guard against a malformed .fai (linebases 0).
    let linebases = if rec.linebases == 0 {
        rec.len
    } else {
        rec.linebases
    };
    let breaks = (rec.len - 1) / linebases;
    let newline_bytes = rec.linewidth.saturating_sub(rec.linebases);
    let span = rec.len + breaks * newline_bytes;
    file.seek(SeekFrom::Start(rec.offset))?;
    let span = usize::try_from(span).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            ".fai span exceeds addressable size",
        )
    })?;
    let mut raw = vec![0u8; span];
    file.read_exact(&mut raw)?;
    let want = rec.len as usize;
    let mut seq: Vec<u8> = Vec::with_capacity(want);
    for &b in &raw {
        if b != b'\n' && b != b'\r' {
            seq.push(b);
            if seq.len() == want {
                break;
            }
        }
    }
    if seq.len() != want {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "FASTA contig {} is shorter than the .fai length {}",
                rec.name, rec.len
            ),
        ));
    }
    Ok(seq)
}

// ---------------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------------

/// BLAKE3 of a whole file, streamed in 64 KiB chunks.
fn blake3_file(path: &Path) -> Result<[u8; 32]> {
    let mut f = File::open(path)?;
    let mut hasher = blake3::Hasher::new();
    let mut buf = vec![0u8; 1 << 16];
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(*hasher.finalize().as_bytes())
}

fn unix_now() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0)
}

/// Read a little-endian u32 from an exact-size chunk (chunks_exact output).
fn le32(c: &[u8]) -> u32 {
    u32::from_le_bytes([c[0], c[1], c[2], c[3]])
}

/// Read a little-endian u64 from an exact-size chunk (chunks_exact output).
fn le64(c: &[u8]) -> u64 {
    u64::from_le_bytes([c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]])
}

/// Read a little-endian u32 header field.
fn field_u32(data: &[u8], off: usize) -> Result<u32> {
    data.get(off..off + 4)
        .map(le32)
        .ok_or_else(|| Error::Truncated("header field out of bounds".to_string()))
}

/// Read a little-endian u64 header field.
fn field_u64(data: &[u8], off: usize) -> Result<u64> {
    data.get(off..off + 8)
        .map(le64)
        .ok_or_else(|| Error::Truncated("header field out of bounds".to_string()))
}

/// Take the next `len`-byte segment starting at `*off`, advancing `*off`.
fn take<'a>(data: &'a [u8], off: &mut usize, len: u64, what: &str) -> Result<&'a [u8]> {
    let len = usize::try_from(len)
        .map_err(|_| Error::Inconsistent(format!("{what} length exceeds the addressable size")))?;
    let end = off
        .checked_add(len)
        .ok_or_else(|| Error::Truncated(format!("{what} segment overruns the file")))?;
    let seg = data
        .get(*off..end)
        .ok_or_else(|| Error::Truncated(format!("{what} segment overruns the file")))?;
    *off = end;
    Ok(seg)
}

// ---------------------------------------------------------------------------
// Reader-side parsing
// ---------------------------------------------------------------------------

/// Validate and parse a whole `.tidx` byte buffer.
fn parse_index(data: &[u8]) -> Result<Tidx> {
    if data.len() < HEADER_LEN {
        return Err(Error::Truncated(format!(
            "file is {} bytes, shorter than the {HEADER_LEN}-byte header",
            data.len()
        )));
    }
    let magic: [u8; 4] = [data[0], data[1], data[2], data[3]];
    let version = field_u32(data, 4)?;
    if magic != MAGIC {
        return Err(Error::BadMagic { magic, version });
    }
    if version != 1 && version != 2 {
        return Err(Error::UnsupportedVersion(version));
    }
    let k = field_u32(data, 8)?;
    if !(MIN_K..=MAX_K).contains(&k) {
        return Err(Error::Inconsistent(format!(
            "k = {k} outside legal range {MIN_K}..={MAX_K}"
        )));
    }
    let m = field_u32(data, 12)?;
    if m != M_FIELD {
        return Err(Error::Inconsistent(format!(
            "header m field is {m}, expected {M_FIELD}"
        )));
    }
    let prefix_bits = field_u32(data, 16)?;
    let max_p = (2 * k).min(63);
    if prefix_bits == 0 || prefix_bits > max_p {
        return Err(Error::Inconsistent(format!(
            "prefix_bits = {prefix_bits} outside legal range 1..={max_p}"
        )));
    }
    let seed = field_u64(data, 24)?;
    if seed != SEED {
        return Err(Error::Inconsistent(format!(
            "header seed is {seed:#x}, expected {SEED:#x}"
        )));
    }
    let tx_count = field_u32(data, 64)?;
    let total_entries = field_u64(data, 80)?;
    let bucket_count = field_u64(data, 88)?;
    let name_blob_len = field_u64(data, 96)?;
    let name_offsets_len = field_u64(data, 104)?;
    if bucket_count != 1u64 << prefix_bits {
        return Err(Error::Inconsistent(format!(
            "bucket_count {bucket_count} != 1 << {prefix_bits}"
        )));
    }
    if name_offsets_len != u64::from(tx_count) + 1 {
        return Err(Error::Inconsistent(format!(
            "name_offsets_len {name_offsets_len} != tx_count + 1 ({})",
            u64::from(tx_count) + 1
        )));
    }
    let low_bits = 2 * k - prefix_bits; // in 1..=63
    let (want_version, key_width) = if low_bits <= 32 {
        (2u32, 4u64)
    } else {
        (1u32, 8u64)
    };
    if version != want_version {
        return Err(Error::Inconsistent(format!(
            "version {version} disagrees with the key width for 2k - P = {low_bits}"
        )));
    }

    // ---- segments -----------------------------------------------------------
    let mut off = HEADER_LEN;
    let blob = take(data, &mut off, name_blob_len, "name blob")?;
    let noff_len = name_offsets_len
        .checked_mul(4)
        .ok_or_else(|| Error::Inconsistent("name offsets length overflows".to_string()))?;
    let noff_seg = take(data, &mut off, noff_len, "name offsets")?;
    let boff_len = bucket_count
        .checked_add(1)
        .and_then(|n| n.checked_mul(8))
        .ok_or_else(|| Error::Inconsistent("bucket offsets length overflows".to_string()))?;
    let boff_seg = take(data, &mut off, boff_len, "bucket offsets")?;
    let keys_len = total_entries
        .checked_mul(key_width)
        .ok_or_else(|| Error::Inconsistent("keys length overflows".to_string()))?;
    let keys_seg = take(data, &mut off, keys_len, "keys")?;
    let pay_len = total_entries
        .checked_mul(8)
        .ok_or_else(|| Error::Inconsistent("payloads length overflows".to_string()))?;
    let pay_seg = take(data, &mut off, pay_len, "payloads")?;
    if off != data.len() {
        return Err(Error::Truncated(format!(
            "file has {} bytes beyond the declared segments",
            data.len() - off
        )));
    }

    // ---- names ----------------------------------------------------------------
    if name_blob_len > u64::from(u32::MAX) {
        return Err(Error::Inconsistent(
            "name blob too large for u32 offsets".to_string(),
        ));
    }
    let name_offsets: Vec<u32> = noff_seg.chunks_exact(4).map(le32).collect();
    if name_offsets.first() != Some(&0) {
        return Err(Error::Inconsistent(
            "first name offset is not 0".to_string(),
        ));
    }
    if name_offsets.last() != Some(&(name_blob_len as u32)) {
        return Err(Error::Inconsistent(
            "last name offset does not equal name_blob_len".to_string(),
        ));
    }
    for pair in name_offsets.windows(2) {
        if pair[0] > pair[1] {
            return Err(Error::Inconsistent(
                "name offsets are not monotone".to_string(),
            ));
        }
    }
    let mut names = Vec::with_capacity(tx_count as usize);
    for i in 0..tx_count as usize {
        let s = name_offsets[i] as usize;
        let e = name_offsets[i + 1] as usize;
        let sl = blob.get(s..e).ok_or_else(|| {
            Error::Inconsistent(format!("transcript name {i} escapes the name blob"))
        })?;
        names.push(
            String::from_utf8(sl.to_vec()).map_err(|_| {
                Error::Inconsistent(format!("transcript name {i} is not valid UTF-8"))
            })?,
        );
    }

    // ---- buckets ----------------------------------------------------------------
    let bucket_offsets: Vec<u64> = boff_seg.chunks_exact(8).map(le64).collect();
    if bucket_offsets.first() != Some(&0) {
        return Err(Error::Inconsistent(
            "first bucket offset is not 0".to_string(),
        ));
    }
    if bucket_offsets.last() != Some(&total_entries) {
        return Err(Error::Inconsistent(
            "last bucket offset does not equal total_entries".to_string(),
        ));
    }
    for pair in bucket_offsets.windows(2) {
        if pair[0] > pair[1] {
            return Err(Error::Inconsistent(
                "bucket offsets are not monotone".to_string(),
            ));
        }
    }

    // ---- keys & payloads -----------------------------------------------------
    let keys: Vec<u64> = if version == 2 {
        keys_seg
            .chunks_exact(4)
            .map(|c| u64::from(le32(c)))
            .collect()
    } else {
        keys_seg.chunks_exact(8).map(le64).collect()
    };
    let payloads: Vec<(u32, u32)> = pay_seg
        .chunks_exact(8)
        .map(|c| {
            let v = le64(c);
            (v as u32, (v >> 32) as u32)
        })
        .collect();

    Ok(Tidx {
        k,
        tx_count,
        ref_blake3: data
            .get(32..64)
            .and_then(|s| <[u8; 32]>::try_from(s).ok())
            .ok_or_else(|| Error::Truncated("reference digest out of bounds".to_string()))?,
        names,
        bucket_offsets,
        keys,
        payloads,
        key_mask: (1u64 << low_bits) - 1,
        bucket_shift: low_bits,
    })
}

// ---------------------------------------------------------------------------
// Public GTF parse (runtime L1Index build path)
// ---------------------------------------------------------------------------

/// One GTF-parsed transcript exposed by [`TranscriptSet`].
#[derive(Debug, Clone, Copy)]
pub struct GtfTranscript<'a> {
    /// GENCODE `transcript_id` (defines the dense tx_id sort order).
    pub transcript_id: &'a str,
    /// Display name (`transcript_name`, falling back to `transcript_id`).
    pub name: &'a str,
    /// Contig name.
    pub contig: &'a str,
    /// Strand byte: `b'+'` or `b'-'`.
    pub strand: u8,
    /// Exons, 0-based half-open `[start, end)`, genomic ascending.
    pub exons: &'a [(u64, u64)],
}

/// Full GTF parse for downstream runtime builds: every exon-bearing
/// transcript in `transcript_id` lexicographic order — exactly the dense
/// `tx_id` space of a `.tidx` built from the same GTF.
#[derive(Debug, Clone)]
pub struct TranscriptSet {
    records: Vec<TxRecord>,
}

impl TranscriptSet {
    /// Parse a GTF file (same semantics as the index builder).
    pub fn parse(path: &Path) -> Result<Self> {
        Ok(TranscriptSet {
            records: parse_gtf(path)?,
        })
    }

    /// Number of transcripts.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// True when the set is empty.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Transcript `i` (index = dense `tx_id`).
    pub fn transcript(&self, i: usize) -> Option<GtfTranscript<'_>> {
        self.records.get(i).map(|r| GtfTranscript {
            transcript_id: &r.id,
            name: &r.name,
            contig: &r.contig,
            strand: r.strand,
            exons: &r.exons,
        })
    }
}
