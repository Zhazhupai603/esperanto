//! esperanto-qc — FASTQ quality control with fastp-compatible semantics.
//!
//! Crate responsibility
//! --------------------
//! Clean RNA-seq FASTQ: fixed-position trimming, optional BWA-style quality
//! trimming, mode-dependent polyG tail trimming, fastp-style adapter
//! trimming (gapless + single-indel tolerant, with paired-end overlap
//! analysis), whole-read filtering (length → N → low quality), paired-end
//! synchronisation, and before/after statistics with a JSON + HTML report.
//!
//! The pipeline is streaming: memory usage is independent of input size.
//! Output is byte-level deterministic.
//!
//! Implementation note
//! -------------------
//! This build runs a *sequential single-threaded* pipeline (no reader thread,
//! no channels). Records are processed in fixed 2048-record chunks; every
//! chunk becomes one independent gzip member in `.fq.gz` outputs (a legal
//! multi-member gzip stream), so output bytes are identical to the chunked
//! threaded design regardless of `QcParams::threads`.
//!
//! Public API
//! ----------
//! * [`QcParams`] — all tunables with fastp defaults.
//! * [`run`] — execute the whole pipeline, writing outputs + `qc.json` +
//!   `qc.html` into `out_dir`.

#![deny(unsafe_code)]

mod trim;
mod detect;

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use serde::Serialize;

/// Number of input records per chunk; one gzip member per chunk.
const CHUNK_RECORDS: usize = 2048;

/// Builtin Illumina adapter table for R1 (TruSeq R1 + Nextera).
pub const BUILTIN_ADAPTERS_R1: [&str; 2] =
    ["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "CTGTCTCTTATACACATCT"];

/// Builtin Illumina adapter table for R2 (TruSeq R2 + Nextera).
pub const BUILTIN_ADAPTERS_R2: [&str; 2] =
    ["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", "CTGTCTCTTATACACATCT"];

/// Minimum overlap for paired-end read-through detection (fastp default).
const OVERLAP_MIN_LEN: usize = 30;

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// All ways the QC pipeline can fail.
#[derive(Debug, thiserror::Error)]
pub enum QcError {
    /// Underlying filesystem / stream error.
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    /// JSON serialization failure while writing the report.
    #[error("json error: {0}")]
    Json(#[from] serde_json::Error),
    /// Parameter combination is not usable.
    #[error("invalid params: {0}")]
    Params(String),
    /// Input looks Phred+64 encoded; only Phred+33 is supported.
    #[error(
        "input appears to be Phred+64 encoded (lowest quality byte {byte} >= 64 in \
         first 10000 reads); only Phred+33 is supported"
    )]
    Phred64 {
        /// Offending ASCII quality byte.
        byte: u8,
    },
    /// FASTQ record is malformed (missing line, bad header, seq/qual
    /// mismatch), or the input contains zero records.
    #[error("malformed FASTQ: {0}")]
    Fastq(String),
    /// R1/R2 record names do not describe the same fragment.
    #[error("paired name mismatch: {r1:?} vs {r2:?}")]
    NameMismatch {
        /// R1 header (sans '@').
        r1: String,
        /// R2 header (sans '@').
        r2: String,
    },
    /// R1 and R2 streams contain different numbers of records.
    #[error("R1/R2 input streams have different record counts")]
    PairLengthMismatch,
}

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

/// polyG trimming mode.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum PolygMode {
    /// Enable iff the instrument detected from the first R1 name is a
    /// two-colour machine (NextSeq / NovaSeq family).
    #[default]
    Auto,
    /// Always trim polyG tails.
    On,
    /// Never trim polyG tails.
    Off,
}

/// Output format for cleaned reads.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum OutFormat {
    /// Gzip-compressed FASTQ (`.clean.fq.gz`), one gzip member per chunk.
    #[default]
    Fqgz,
    /// Uncompressed little-endian binary FASTQ (`.clean.bfq`).
    Bfq,
}

/// All QC pipeline tunables. Defaults follow fastp.
#[derive(Debug, Clone)]
pub struct QcParams {
    /// R1 inputs (plain or `.gz`); multiple lanes are concatenated in order.
    pub r1: Vec<PathBuf>,
    /// R2 inputs; empty = single-end mode, otherwise `r2.len() == r1.len()`.
    pub r2: Vec<PathBuf>,
    /// Directory receiving all outputs.
    pub out_dir: PathBuf,
    /// Enable adapter trimming (default true).
    pub adapter_trim: bool,
    /// User-supplied R1 adapter table; empty = builtin Illumina table.
    pub adapters_r1: Vec<String>,
/// User-supplied R2 adapter table; empty = builtin Illumina table.
    pub adapters_r2: Vec<String>,
    /// SE adapter auto-detection when no user table is given (SE only,
    /// default false). When false every output is byte-identical to the
    /// feature being absent.
    pub detect_adapter_se: bool,
    /// PE overlap analysis before table matching (default true, PE only).
    pub pe_overlap: bool,
    /// Enable BWA-style 3' quality trimming (default false).
    pub qtrim: bool,
    /// Quality cutoff for qtrim (default 20).
    pub qtrim_cutoff: u8,
    /// Fixed bases to drop from the 5' end of R1.
    pub trim_front1: usize,
    /// Fixed bases to drop from the 3' end of R1.
    pub trim_tail1: usize,
    /// Fixed bases to drop from the 5' end of R2.
    pub trim_front2: usize,
    /// Fixed bases to drop from the 3' end of R2.
    pub trim_tail2: usize,
    /// polyG tail trimming mode (default Auto).
    pub polyg: PolygMode,
    /// Reads shorter than this fail the length filter (default 15).
    pub min_len: usize,
    /// Reads with more N bases than this fail (default 5).
    pub n_max: usize,
    /// Max fraction of bases with Q < 15 (default 0.4; strictly greater fails).
    pub q15_frac_max: f64,
    /// Keep the passing mate of a failed pair in `unpaired1/2` outputs.
    pub keep_unpaired: bool,
    /// Worker thread hint (0 = all cores). Sequential build: accepted,
    /// recorded in the report, not used for scheduling.
    pub threads: usize,
    /// Output container for cleaned reads.
    pub out_format: OutFormat,
}

impl Default for QcParams {
    fn default() -> Self {
        Self {
            r1: Vec::new(),
            r2: Vec::new(),
            out_dir: PathBuf::new(),
            adapter_trim: true,
            adapters_r1: Vec::new(),
            adapters_r2: Vec::new(),
            detect_adapter_se: false,
            pe_overlap: true,
            qtrim: false,
            qtrim_cutoff: 20,
            trim_front1: 0,
            trim_tail1: 0,
            trim_front2: 0,
            trim_tail2: 0,
            polyg: PolygMode::Auto,
            min_len: 15,
            n_max: 5,
            q15_frac_max: 0.4,
            keep_unpaired: false,
            threads: 0,
            out_format: OutFormat::Fqgz,
        }
    }
}

// ---------------------------------------------------------------------------
// FASTQ reading
// ---------------------------------------------------------------------------

/// One parsed FASTQ record. `name` keeps the full header after '@'
/// including any comment; sequence and quality are equal-length raw bytes.
struct Record {
    name: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

/// Reader over a single FASTQ file (plain or gzip, multi-member safe).
struct FastqReader {
    inner: BufReader<Box<dyn Read>>,
    file: String,
    line: u64,
}

impl FastqReader {
    fn open(path: &Path) -> Result<Self, QcError> {
        let f = File::open(path)?;
        let reader: Box<dyn Read> = if path.extension().and_then(|e| e.to_str()) == Some("gz") {
            Box::new(MultiGzDecoder::new(f))
        } else {
            Box::new(f)
        };
        Ok(Self {
            inner: BufReader::with_capacity(1 << 16, reader),
            file: path.display().to_string(),
            line: 0,
        })
    }

    /// Read one raw line (without trailing \n / \r). Ok(false) = clean EOF.
    fn read_line(&mut self, out: &mut Vec<u8>) -> Result<bool, QcError> {
        out.clear();
        let n = self.inner.read_until(b'\n', out)?;
        if n == 0 {
            return Ok(false);
        }
        self.line += 1;
        while matches!(out.last(), Some(b'\n') | Some(b'\r')) {
            out.pop();
        }
        Ok(true)
    }

    /// Next record, or None at end of this file.
    fn next_record(&mut self) -> Result<Option<Record>, QcError> {
        let rec_start = self.line + 1;
        let mut line = Vec::with_capacity(300);

        if !self.read_line(&mut line)? {
            return Ok(None);
        }
        if line.first() != Some(&b'@') {
            return Err(QcError::Fastq(format!(
                "{}:{}: bad record",
                self.file.clone(),
                rec_start
            )));
        }
        let name = line[1..].to_vec();

        if !self.read_line(&mut line)? {
            return Err(self.truncated(rec_start));
        }
        let seq = line.clone();

        if !self.read_line(&mut line)? {
            return Err(self.truncated(rec_start));
        }
        if line.first() != Some(&b'+') {
            return Err(QcError::Fastq(format!(
                "{}:{}: '+' line missing",
                self.file.clone(),
                rec_start + 2
            )));
        }

        if !self.read_line(&mut line)? {
            return Err(self.truncated(rec_start));
        }
        let qual = line;

        if seq.len() != qual.len() {
            return Err(QcError::Fastq(format!(
                "{}:{}: bad record",
                self.file.clone(),
                rec_start
            )));
        }
        Ok(Some(Record { name, seq, qual }))
    }

    fn truncated(&self, rec_start: u64) -> QcError {
        QcError::Fastq(format!(
            "{}:{}: truncated record",
            self.file.clone(),
            rec_start
        ))
    }
}

/// Logical stream over several FASTQ files concatenated in order.
struct FastqSet {
    files: Vec<PathBuf>,
    idx: usize,
    cur: Option<FastqReader>,
}

impl FastqSet {
    fn new(files: &[PathBuf]) -> Self {
        Self {
            files: files.to_vec(),
            idx: 0,
            cur: None,
        }
    }

    fn next(&mut self) -> Result<Option<Record>, QcError> {
        loop {
            if self.cur.is_none() {
                if self.idx >= self.files.len() {
                    return Ok(None);
                }
                let path = self.files[self.idx].clone();
                self.idx += 1;
                self.cur = Some(FastqReader::open(&path)?);
            }
            if let Some(r) = self.cur.as_mut() {
                if let Some(rec) = r.next_record()? {
                    return Ok(Some(rec));
                }
            }
            self.cur = None;
        }
    }
}

// ---------------------------------------------------------------------------
// Adapter matching (fastp trimBySequence semantics)
// ---------------------------------------------------------------------------

/// Try to align `adapter` inside `read`; return the cut point (number of
/// leading read bases to keep), `Some(0)` = the whole read is adapter.
///
/// Scan order: gapless, then read-has-one-inserted-base, then
/// read-is-missing-one-adapter-base. First hit wins.
pub(crate) fn match_adapter(read: &[u8], adapter: &[u8]) -> Option<usize> {
    let alen = adapter.len();
    let rlen = read.len();
    if alen < 4 || rlen < 4 {
        return None;
    }
    // Negative start offsets tolerate reads that begin inside the adapter
    // (A-tailing / adapter dimers): leading positions that fall before the
    // read start are free matches.
    let start: isize = if alen >= 16 {
        -4
    } else if alen >= 12 {
        -3
    } else if alen >= 8 {
        -2
    } else {
        0
    };

    // --- gapless: compare adapter prefix of cmplen = min(rlen-pos, alen) ---
    let mut pos = start;
    while rlen as isize - pos > 4 {
        let cmplen = ((rlen as isize - pos) as usize).min(alen);
        let budget = cmplen / 8;
        let mut mm = 0usize;
        for (i, &ab) in adapter.iter().enumerate().take(cmplen) {
            let ridx = pos + i as isize;
            if ridx >= 0 && !read[ridx as usize].eq_ignore_ascii_case(&ab) {
                mm += 1;
                if mm > budget {
                    break;
                }
            }
        }
        if mm <= budget {
            return Some(pos.max(0) as usize);
        }
        pos += 1;
    }

    // --- single-base indel variants, budget cmplen/8 - 1 (saturating) ---
    if let Some(hit) = match_indel(read, adapter, start, true) {
        return Some(hit);
    }
    match_indel(read, adapter, start, false)
}

/// Single-indel scan. `read_has_extra = true`: the read carries one inserted
/// base relative to the adapter (skip one read base); `false`: one adapter
/// base is missing from the read (skip one adapter base).
fn match_indel(read: &[u8], adapter: &[u8], start: isize, read_has_extra: bool) -> Option<usize> {
    let rlen = read.len() as isize;
    let alen = adapter.len() as isize;
    // Old-tool bounds: indel variants scan non-negative offsets only; the last
    // position must leave at least 5 compared bases (insertion consumes one
    // extra read base, so its bound is one tighter than the deletion's).
    let mut pos = start.max(0);
    loop {
        // Old-tool position bounds: insertion pos < rlen-5, deletion pos < rlen-4
        // (checked before the window is tried, so the minimum compared span in
        // range is 5 bases).
        let pos_limit = if read_has_extra { 5 } else { 4 };
        if rlen - pos <= pos_limit {
            break;
        }
        let read_span = if read_has_extra {
            rlen - pos - 1
        } else {
            rlen - pos
        };
        let ad_span = if read_has_extra { alen } else { alen - 1 };
        let cmplen = read_span.min(ad_span);
        if cmplen < 4 {
            break;
        }
        let budget = (cmplen as usize / 8).saturating_sub(1);
        // Old-tool fast-path gate (documented deviation): the indel variant is
        // only attempted where the UNGAPPED first min(8, cmplen) bases already
        // match within the variant budget. Ported as a hard filter for parity.
        let pre = 8.min(cmplen as usize);
        let mut ungapped = 0usize;
        let mut gated = false;
        for (i, &ab) in adapter.iter().enumerate().take(pre) {
            let ridx = pos + i as isize;
            if ridx < 0 {
                continue;
            }
            if !read[ridx as usize].eq_ignore_ascii_case(&ab) {
                ungapped += 1;
                if ungapped > budget {
                    gated = true;
                    break;
                }
            }
        }
        if gated {
            pos += 1;
            continue;
        }
        // Old-tool k range: the single indel must be strictly interior to the
        // compared window (k = 1..cmplen); skipping the first/last base would
        // shift the window instead.
        for k in 1..cmplen as usize {
            let mut mm = 0usize;
            for i in 0..cmplen as usize {
                let (ridx, aidx) = if read_has_extra {
                    (pos + i as isize + isize::from(i >= k), i)
                } else {
                    (pos + i as isize, i + usize::from(i >= k))
                };
                if ridx < 0 {
                    continue; // read starts inside the adapter: free match
                }
                if !read[ridx as usize].eq_ignore_ascii_case(&adapter[aidx]) {
                    mm += 1;
                    if mm > budget {
                        break;
                    }
                }
            }
            if mm <= budget {
                return Some(pos.max(0) as usize);
            }
        }
        pos += 1;
    }
    None
}

/// Apply the adapter table sequentially: a hit truncates the read, then the
/// next adapter is tried against the truncated sequence. Returns whether any
/// adapter hit (and counts reads/bases removed through `cnt`).
fn adapter_stage(rec: &mut Record, table: &[Vec<u8>], cnt: &mut Counters) -> bool {
    if table.is_empty() {
        return false;
    }
    let mut hit_any = false;
    let mut removed = 0usize;
    for ad in table {
        if let Some(cut) = match_adapter(&rec.seq, ad) {
            hit_any = true;
            if cut < rec.seq.len() {
                removed += rec.seq.len() - cut;
                rec.seq.truncate(cut);
                rec.qual.truncate(cut);
            }
        }
    }
    if hit_any {
        cnt.adapter_reads += 1;
        cnt.adapter_bases += removed as u64;
    }
    hit_any
}

// ---------------------------------------------------------------------------
// Paired-end overlap analysis
// ---------------------------------------------------------------------------

/// Reverse complement (case-preserving; unknown symbols pass through).
/// Scientific semantics: a base is complemented regardless of case,
/// preserving case (a<->t, c<->g); non-ACGT bytes pass through unchanged.
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b {
        b'A' => b'T',
        b'a' => b't',
        b'C' => b'G',
        b'c' => b'g',
        b'G' => b'C',
        b'g' => b'c',
        b'T' => b'A',
        b't' => b'a',
        other => other,
    }).collect()
}

/// Align rc(R2) against R1: offsets 0→+ first, then 0→− (read-through).
/// An offset wins when `ol = min(l1,l2,l1-o,l2+o) >= 30` and mismatches
/// `<= min(5, ol/10)`. Returns `(offset, overlap_len)`.
fn find_overlap(r1: &[u8], r2: &[u8]) -> Option<(isize, usize)> {
    let l1 = r1.len() as isize;
    let l2 = r2.len() as isize;
    if l1 <= OVERLAP_MIN_LEN as isize || l2 <= OVERLAP_MIN_LEN as isize {
        return None;
    }
    let rc2 = revcomp(r2);

    let check = |o: isize| -> Option<usize> {
        let ol = l1.min(l2).min(l1 - o).min(l2 + o);
        if ol < OVERLAP_MIN_LEN as isize {
            return None;
        }
        let max_mm = 5usize.min(ol as usize / 10);
        let (s1, s2): (&[u8], &[u8]) = if o >= 0 {
            (&r1[o as usize..], &rc2[..])
        } else {
            (r1, &rc2[(-o) as usize..])
        };
        let mut mm = 0usize;
        for i in 0..ol as usize {
            if !s1[i].eq_ignore_ascii_case(&s2[i]) {
                mm += 1;
                if mm > max_mm {
                    return None;
                }
            }
        }
        Some(ol as usize)
    };

    // Old-tool scan bounds (exclusive upper limits): positive offsets run
    // 0..l1-30, negative offsets run 1..l2-30. Boundary positions where the
    // overlap would shrink to exactly 30 are NOT tried — spurious 30-mer
    // overlaps in periodic adapter tails land exactly there.
    let mut o = 0isize;
    while o < l1 - OVERLAP_MIN_LEN as isize {
        if let Some(ol) = check(o) {
            return Some((o, ol));
        }
        o += 1;
    }
    let mut o = -1isize;
    while o > -(l2 - OVERLAP_MIN_LEN as isize) {
        if let Some(ol) = check(o) {
            return Some((o, ol));
        }
        o -= 1;
    }
    None
}

// ---------------------------------------------------------------------------
// Per-read trimming pipeline
// ---------------------------------------------------------------------------

/// Per-end trim options (fixed positions + qtrim + polyG). Adapter handling
/// lives with the caller because PE overlap analysis decides the table.
struct EndOpts<'a> {
    front: usize,
    tail: usize,
    qtrim: bool,
    qtrim_cutoff: u8,
    polyg: bool,
    adapters: &'a [Vec<u8>],
}

/// Fixed trim → qtrim → polyG → (if `opts.adapters` non-empty) table trim.
/// Returns whether adapter trimming removed anything (for `AdapterOnly`).
fn process_end(rec: &mut Record, opts: &EndOpts, cnt: &mut Counters) -> bool {
    trim::fixed_trim(&mut rec.seq, &mut rec.qual, opts.front, opts.tail);

    if opts.qtrim {
        let keep = trim::qtrim_tail(&rec.qual, opts.qtrim_cutoff);
        cnt.qtrim_bases += (rec.qual.len() - keep) as u64;
        rec.seq.truncate(keep);
        rec.qual.truncate(keep);
    }

    if opts.polyg {
        let keep = trim::poly_g_trim(&rec.seq);
        if keep < rec.seq.len() {
            cnt.polyg_reads += 1;
            cnt.polyg_bases += (rec.seq.len() - keep) as u64;
            rec.seq.truncate(keep);
            rec.qual.truncate(keep);
        }
    }

    adapter_stage(rec, opts.adapters, cnt)
}

// ---------------------------------------------------------------------------
// Whole-read filtering
// ---------------------------------------------------------------------------

/// First failing check wins; order: length → N → low quality.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FailReason {
    LowQuality,
    TooManyN,
    TooShort,
    AdapterOnly,
}

fn filter_end(
    seq: &[u8],
    qual: &[u8],
    adapter_trimmed: bool,
    min_len: usize,
    n_max: usize,
    q15_frac_max: f64,
) -> Option<FailReason> {
    if seq.len() < min_len {
        return Some(if seq.is_empty() && adapter_trimmed {
            FailReason::AdapterOnly
        } else {
            FailReason::TooShort
        });
    }
    let n = seq
        .iter()
        .filter(|&&c| c == b'N' || c == b'n')
        .count();
    if n > n_max {
        return Some(FailReason::TooManyN);
    }
    // Q < 15 ⇔ Phred+33 byte < 48; strictly greater than the allowed
    // fraction fails, exactly-equal is kept.
    let q15 = qual.iter().filter(|&&q| q < 48).count();
    if (q15 as f64) > seq.len() as f64 * q15_frac_max {
        return Some(FailReason::LowQuality);
    }
    None
}

// ---------------------------------------------------------------------------
// Statistics
// ---------------------------------------------------------------------------

/// Integer accumulators for one "side" (before or after). Per-cycle arrays
/// grow to the longest read seen.
#[derive(Default)]
struct Stats {
    reads: u64,
    bases: u64,
    q20: u64,
    q30: u64,
    gc: u64,
    n_bases: u64,
    cycle_qual_sum: Vec<u64>,
    cycle_count: Vec<u64>,
    /// Per-cycle base counts in order [A, C, G, T, N].
    cycle_bases: Vec<[u64; 5]>,
}

impl Stats {
    fn add(&mut self, seq: &[u8], qual: &[u8]) {
        self.reads += 1;
        let n = seq.len();
        if self.cycle_count.len() < n {
            self.cycle_qual_sum.resize(n, 0);
            self.cycle_count.resize(n, 0);
            self.cycle_bases.resize_with(n, || [0; 5]);
        }
        for (i, (&s, &q)) in seq.iter().zip(qual.iter()).enumerate() {
            self.bases += 1;
            let ph = q.saturating_sub(33) as u64;
            if ph >= 20 {
                self.q20 += 1;
            }
            if ph >= 30 {
                self.q30 += 1;
            }
            // Scientific semantics: sequence case carries no
            // information — lowercase acgt are the same bases. They land in
            // their letter buckets; only N/n (and other codes) land in the N
            // bucket. gc counts C/c/G/g.
            let bucket = match s {
                b'A' | b'a' => 0,
                b'C' | b'c' => {
                    self.gc += 1;
                    1
                }
                b'G' | b'g' => {
                    self.gc += 1;
                    2
                }
                b'T' | b't' => 3,
                _ => {
                    self.n_bases += 1;
                    4
                }
            };
            self.cycle_qual_sum[i] += ph;
            self.cycle_count[i] += 1;
            self.cycle_bases[i][bucket] += 1;
        }
    }
}

/// Pipeline-wide counters.
#[derive(Default)]
struct Counters {
    adapter_reads: u64,
    adapter_bases: u64,
    polyg_reads: u64,
    polyg_bases: u64,
    qtrim_bases: u64,
    low_quality: u64,
    too_many_n: u64,
    too_short: u64,
    adapter_only: u64,
    unpaired_written: u64,
}

impl Counters {
    fn fail(&mut self, reason: FailReason) {
        match reason {
            FailReason::LowQuality => self.low_quality += 1,
            FailReason::TooManyN => self.too_many_n += 1,
            FailReason::TooShort => self.too_short += 1,
            FailReason::AdapterOnly => self.adapter_only += 1,
        }
    }
}

// ---------------------------------------------------------------------------
// Output sinks
// ---------------------------------------------------------------------------

/// Buffered fq.gz writer: one independent gzip member per input chunk.
struct GzChunkWriter {
    file: BufWriter<File>,
    buf: Vec<u8>,
    wrote_any: bool,
}

impl GzChunkWriter {
    fn create(path: &Path) -> Result<Self, QcError> {
        Ok(Self {
            file: BufWriter::new(File::create(path)?),
            buf: Vec::with_capacity(1 << 20),
            wrote_any: false,
        })
    }

    fn push(&mut self, rec: &Record) {
        self.buf.extend_from_slice(b"@");
        self.buf.extend_from_slice(&rec.name);
        self.buf.extend_from_slice(b"\n");
        self.buf.extend_from_slice(&rec.seq);
        self.buf.extend_from_slice(b"\n+\n");
        self.buf.extend_from_slice(&rec.qual);
        self.buf.extend_from_slice(b"\n");
    }

    /// Close the current chunk as one gzip member. Only chunks with records
    /// emit a member; a stream that produced zero records leaves the file at
    /// 0 bytes (reference behavior: old tool writes an empty file, not an
    /// empty gzip member).
    fn end_chunk(&mut self) -> Result<(), QcError> {
        if !self.buf.is_empty() {
            let mut enc = GzEncoder::new(Vec::new(), Compression::default());
            enc.write_all(&self.buf)?;
            let member = enc.finish()?;
            self.file.write_all(&member)?;
            self.wrote_any = true;
        }
        self.buf.clear();
        Ok(())
    }

    fn finish(mut self) -> Result<(), QcError> {
        self.end_chunk()?;
        self.file.flush()?;
        Ok(())
    }
}

/// Binary FASTQ writer. 16-byte header: `EBFQ` magic, version u8 = 1,
/// 3 reserved zero bytes, read_count u64 LE (patched via seek on close).
/// Each record: name_len u16 LE + name + seq_len u16 LE + seq + qual.
/// All fields little-endian, no compression.
struct BfqWriter {
    w: BufWriter<File>,
    count: u64,
}

impl BfqWriter {
    fn create(path: &Path) -> Result<Self, QcError> {
        let mut w = BufWriter::new(File::create(path)?);
        w.write_all(b"EBFQ")?;
        w.write_all(&[1, 0, 0, 0])?;
        w.write_all(&0u64.to_le_bytes())?;
        Ok(Self { w, count: 0 })
    }

    fn push(&mut self, rec: &Record) -> Result<(), QcError> {
        let name_len = fit_u16(rec.name.len())?;
        let seq_len = fit_u16(rec.seq.len())?;
        self.w.write_all(&name_len.to_le_bytes())?;
        self.w.write_all(&rec.name)?;
        self.w.write_all(&seq_len.to_le_bytes())?;
        self.w.write_all(&rec.seq)?;
        self.w.write_all(&rec.qual)?;
        self.count += 1;
        Ok(())
    }

    /// Close WITHOUT patching the count (placeholder 0 remains) — mirrors the
    /// reference tool's R2/unpaired sinks.
    fn finish(mut self) -> Result<(), QcError> {
        self.w.flush()?;
        self.w.get_mut().flush()?;
        Ok(())
    }

    /// Close patching `read_count = count` (R1/SE sink; see call site for the
    /// bug-compat total-count note).
    fn finish_with_count(mut self, count: u64) -> Result<(), QcError> {
        self.w.flush()?;
        let file = self.w.get_mut();
        file.seek(SeekFrom::Start(8))?;
        file.write_all(&count.to_le_bytes())?;
        file.flush()?;
        Ok(())
    }
}

fn fit_u16(n: usize) -> Result<u16, QcError> {
    n.try_into().map_err(|_| {
        QcError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "bfq field longer than 65535 bytes",
        ))
    })
}

enum OutSink {
    Fq(GzChunkWriter),
    Bfq(BfqWriter),
}

impl OutSink {
    fn create(fmt: OutFormat, path: &Path) -> Result<Self, QcError> {
        Ok(match fmt {
            OutFormat::Fqgz => OutSink::Fq(GzChunkWriter::create(path)?),
            OutFormat::Bfq => OutSink::Bfq(BfqWriter::create(path)?),
        })
    }

    fn write(&mut self, rec: &Record) -> Result<(), QcError> {
        match self {
            OutSink::Fq(w) => {
                w.push(rec);
                Ok(())
            }
            OutSink::Bfq(w) => w.push(rec),
        }
    }

    fn end_chunk(&mut self) -> Result<(), QcError> {
        match self {
            OutSink::Fq(w) => w.end_chunk(),
            OutSink::Bfq(_) => Ok(()),
        }
    }

    fn finish(self) -> Result<(), QcError> {
        match self {
            OutSink::Fq(w) => w.finish(),
            OutSink::Bfq(w) => w.finish(),
        }
    }

    /// Bfq sinks only: close patching read_count (reference bug-compat).
    fn finish_with_count(self, count: u64) -> Result<(), QcError> {
        match self {
            OutSink::Fq(w) => w.finish(),
            OutSink::Bfq(w) => w.finish_with_count(count),
        }
    }
}

// ---------------------------------------------------------------------------
// Report (qc.json / qc.html)
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct BaseFrac {
    a: f64,
    c: f64,
    g: f64,
    t: f64,
    n: f64,
}

#[derive(Serialize)]
struct CycleEntry {
    position: usize,
    mean_quality: f64,
    base_frac: BaseFrac,
}

#[derive(Serialize)]
struct ParamsReport {
    adapter_trim: bool,
    pe_overlap: bool,
    qtrim: bool,
    polyg: String,
    min_len: usize,
    n_max: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    detect_adapter_se: Option<bool>,
    q15_frac_max: f64,
}

#[derive(Serialize)]
struct InputReport {
    r1: Vec<String>,
    r2: Vec<String>,
    instrument_polyg: bool,
}

/// After-side summary: counts on both sides, rates computed on "after" data.
#[derive(Serialize)]
struct Summary {
    reads_before: u64,
    reads_after: u64,
    bases_before: u64,
    bases_after: u64,
    q20_rate: f64,
    q30_rate: f64,
    gc_rate: f64,
    n_rate: f64,
}

/// Before-side summary (single-sided: everything is "before").
#[derive(Serialize)]
struct SummaryBefore {
    q20_rate: f64,
    q30_rate: f64,
    gc_rate: f64,
    n_rate: f64,
}

#[derive(Serialize)]
struct FilterReasons {
    low_quality: u64,
    too_many_n: u64,
    too_short: u64,
    adapter_only: u64,
}

#[derive(Serialize)]
struct Trimming {
    adapter_reads: u64,
    adapter_bases: u64,
    polyg_reads: u64,
    polyg_bases: u64,
    qtrim_bases: u64,
}

#[derive(Serialize)]
struct QcReport {
    esperanto_qc_version: String,
    params: ParamsReport,
    input: InputReport,
    summary: Summary,
    filter_reasons: FilterReasons,
    trimming: Trimming,
    per_cycle: Vec<CycleEntry>,
    per_cycle_before: Vec<CycleEntry>,
    summary_before: SummaryBefore,
    elapsed_seconds: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    adapter_source: Option<&'static str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    detected_adapter_se: Option<String>,
}

fn div(n: u64, d: u64) -> f64 {
    if d == 0 {
        0.0
    } else {
        n as f64 / d as f64
    }
}

fn cycle_entries(stats: &Stats) -> Vec<CycleEntry> {
    (0..stats.cycle_count.len())
        .map(|i| {
            let n = stats.cycle_count[i] as f64;
            let b = &stats.cycle_bases[i];
            CycleEntry {
                position: i + 1,
                mean_quality: stats.cycle_qual_sum[i] as f64 / n,
                base_frac: BaseFrac {
                    a: b[0] as f64 / n,
                    c: b[1] as f64 / n,
                    g: b[2] as f64 / n,
                    t: b[3] as f64 / n,
                    n: b[4] as f64 / n,
                },
            }
        })
        .collect()
}

fn summary_of(before: &Stats, after: &Stats) -> Summary {
    Summary {
        reads_before: before.reads,
        reads_after: after.reads,
        bases_before: before.bases,
        bases_after: after.bases,
        q20_rate: div(after.q20, after.bases),
        q30_rate: div(after.q30, after.bases),
        gc_rate: div(after.gc, after.bases),
        n_rate: div(after.n_bases, after.bases),
    }
}

fn summary_before_of(before: &Stats) -> SummaryBefore {
    SummaryBefore {
        q20_rate: div(before.q20, before.bases),
        q30_rate: div(before.q30, before.bases),
        gc_rate: div(before.gc, before.bases),
        n_rate: div(before.n_bases, before.bases),
    }
}

/// Minimal standalone HTML report embedding the full qc.json payload.
fn render_html(json: &str, before: &Stats, after: &Stats) -> String {
    let row = |label: &str, b: String, a: String| {
        format!("<tr><td>{label}</td><td>{b}</td><td>{a}</td></tr>\n")
    };
    let mut table = String::from("<tr><th>metric</th><th>before</th><th>after</th></tr>\n");
    table.push_str(&row(
        "reads",
        before.reads.to_string(),
        after.reads.to_string(),
    ));
    table.push_str(&row(
        "bases",
        before.bases.to_string(),
        after.bases.to_string(),
    ));
    table.push_str(&row(
        "q20_rate",
        (div(before.q20, before.bases)).to_string(),
        (div(after.q20, after.bases)).to_string(),
    ));
    table.push_str(&row(
        "q30_rate",
        (div(before.q30, before.bases)).to_string(),
        (div(after.q30, after.bases)).to_string(),
    ));
    table.push_str(&row(
        "gc_rate",
        (div(before.gc, before.bases)).to_string(),
        (div(after.gc, after.bases)).to_string(),
    ));
    table.push_str(&row(
        "n_rate",
        (div(before.n_bases, before.bases)).to_string(),
        (div(after.n_bases, after.bases)).to_string(),
    ));
    // Escape '<' so "</script>" cannot terminate the payload early.
    let embedded = json.replace('<', "\\u003c");
    format!(
        "<!doctype html>\n<html>\n<head>\n<meta charset=\"utf-8\">\n\
         <title>esperanto-qc report</title>\n\
         <style>body{{font-family:monospace;margin:2em}}\
         table{{border-collapse:collapse}}\
         td,th{{border:1px solid #999;padding:4px 12px;text-align:right}}</style>\n\
         </head>\n<body>\n<h1>esperanto-qc report</h1>\n<table>\n{table}</table>\n\
         <script type=\"application/json\" id=\"qc-json\">{embedded}</script>\n\
         </body>\n</html>\n"
    )
}

// ---------------------------------------------------------------------------
// Pre-scan: phred encoding check + instrument detection
// ---------------------------------------------------------------------------

/// Reject Phred+64-looking input (any quality byte >= 75 in the first
/// 10000 reads) and detect the sequencing instrument from the first R1 name.
fn prescan(params: &QcParams, paired: bool) -> Result<Option<String>, QcError> {
    let mut s1 = FastqSet::new(&params.r1);
    let mut s2 = if paired {
        Some(FastqSet::new(&params.r2))
    } else {
        None
    };
    let mut first_name: Option<Vec<u8>> = None;
    let mut min_char = u8::MAX;
    let mut n_scanned = 0usize;
    for _ in 0..10000 {
        let a = s1.next()?;
        let b = if paired {
            match s2.as_mut() {
                Some(s) => s.next()?,
                None => None,
            }
        } else {
            None
        };
        if !paired {
            match a {
                None => break,
                Some(a) => {
                    if first_name.is_none() {
                        first_name = Some(a.name.clone());
                    }
                    check_phred33(&a.qual, &mut min_char)?;
                    n_scanned += 1;
                }
            }
        } else {
            match (a, b) {
                (None, None) => break,
                (Some(a), Some(b)) => {
                    if first_name.is_none() {
                        first_name = Some(a.name.clone());
                    }
                    check_phred33(&a.qual, &mut min_char)?;
                    check_phred33(&b.qual, &mut min_char)?;
                    n_scanned += 1;
                }
                (Some(_), None) | (None, Some(_)) => return Err(QcError::PairLengthMismatch),
            }
        }
    }
    if n_scanned == 0 {
        // Reference behavior: an input with zero records is a parse error.
        return Err(QcError::Fastq("empty input: no FASTQ records".to_string()));
    }
    if min_char >= 64 {
        return Err(QcError::Phred64 { byte: min_char });
    }
    Ok(first_name.as_deref().and_then(detect_instrument))
}

/// Reference rule: every quality byte must lie in [33,126]; the caller also
/// rejects the file when the MINIMUM byte over the prescan window is >= 64
/// (all-'@'-or-higher = suspected Phred+64). A single byte >= 75 alone is
/// NOT an error (legitimate Q42+ exists in +33 files).
fn check_phred33(qual: &[u8], min_char: &mut u8) -> Result<(), QcError> {
    for &q in qual {
        if !(33..=126).contains(&q) {
            return Err(QcError::Phred64 { byte: q });
        }
        *min_char = (*min_char).min(q);
    }
    Ok(())
}

/// Two-colour instruments (NextSeq / NovaSeq families) suffer dark cycles
/// that read as G, so polyG Auto mode enables trimming for them. The
/// instrument token is the first ':'-separated field of the header.
fn detect_instrument(name: &[u8]) -> Option<String> {
    // Ported from the reference tool: token = first ':' field of the header;
    // A0<digit>* => NovaSeq family, NS5*/NB5* => NextSeq family, or the FULL
    // header containing the case-sensitive literal "NextSeq"/"NovaSeq".
    let head = before_space(name);
    let token = match head.iter().position(|&c| c == b':') {
        Some(i) => head[..i].to_vec(),
        None => head.to_vec(),
    };
    if token.is_empty() {
        return None;
    }
    let tok = token.as_slice();
    let two_colour = (tok.starts_with(b"A0")
        && tok.len() >= 3
        && tok[2].is_ascii_digit())
        || tok.starts_with(b"NS5")
        || tok.starts_with(b"NB5")
        || windows_contains(name, b"NextSeq")
        || windows_contains(name, b"NovaSeq");
    if two_colour {
        Some(String::from_utf8_lossy(&token).into_owned())
    } else {
        None
    }
}

/// Case-sensitive substring test on byte slices.
fn windows_contains(hay: &[u8], needle: &[u8]) -> bool {
    hay.windows(needle.len()).any(|w| w == needle)
}

fn before_space(x: &[u8]) -> &[u8] {
    match x.iter().position(|&c| c == b' ') {
        Some(i) => &x[..i],
        None => x,
    }
}

/// Reference rule (core_name): compare the first whitespace-delimited token
/// with a single trailing "/1" or "/2" stripped from EACH side independently.
fn paired_names_match(a: &[u8], b: &[u8]) -> bool {
    core_name(before_space(a)) == core_name(before_space(b))
}

fn core_name(tok: &[u8]) -> &[u8] {
    if tok.ends_with(b"/1") || tok.ends_with(b"/2") {
        &tok[..tok.len() - 2]
    } else {
        tok
    }
}

/// Output stem: R1 file name with trailing `.gz` / `.fastq` / `.fq`
/// extensions stripped (repeatedly, so `x.fastq.gz` → `x`).
fn stem_of(path: &Path) -> String {
    let mut name = path
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    loop {
        if name.ends_with(".gz") {
            name.truncate(name.len() - 3);
        } else if name.ends_with(".fastq") {
            name.truncate(name.len() - 6);
        } else if name.ends_with(".fq") {
            name.truncate(name.len() - 3);
        } else {
            break;
        }
    }
    if name.is_empty() {
        "reads".to_string()
    } else {
        name
    }
}

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------

/// Run the full QC pipeline described by `params`.
///
/// Writes `<stem>.clean.fq.gz` (or `.clean.bfq`) for R1 — plus
/// `<stem>.clean.2.fq.gz` for R2 in PE mode and `<stem>.unpaired1/2.*`
/// when `keep_unpaired` is set — and `qc.json` / `qc.html` into `out_dir`.
/// `stem` derives from the first R1 file name.
pub fn run(params: &QcParams) -> Result<(), QcError> {
    let t0 = Instant::now();

    // --- validation ---
    if params.r1.is_empty() {
        return Err(QcError::Params("r1 must not be empty".into()));
    }
    let paired = !params.r2.is_empty();
    if paired && params.r2.len() != params.r1.len() {
        return Err(QcError::Params(
            "r2 must be empty or have the same length as r1".into(),
        ));
    }
    if !params.q15_frac_max.is_finite() || params.q15_frac_max < 0.0 {
        return Err(QcError::Params(
            "q15_frac_max must be finite and >= 0".into(),
        ));
    }
    std::fs::create_dir_all(&params.out_dir)?;
    let _ = params.threads; // sequential build: not used for scheduling

    // --- prescan: phred encoding + instrument for polyG Auto ---
    let instrument = prescan(params, paired)?;
    let polyg_on = match params.polyg {
        PolygMode::On => true,
        PolygMode::Off => false,
        PolygMode::Auto => instrument.is_some(),
    };

    // --- resolve adapter tables ---
    let user1: Vec<Vec<u8>> = params
        .adapters_r1
        .iter()
        .map(|s| s.as_bytes().to_vec())
        .collect();
    let user2: Vec<Vec<u8>> = params
        .adapters_r2
        .iter()
        .map(|s| s.as_bytes().to_vec())
        .collect();
    let builtin1: Vec<Vec<u8>> = BUILTIN_ADAPTERS_R1
        .iter()
        .map(|s| s.as_bytes().to_vec())
        .collect();
    let builtin2: Vec<Vec<u8>> = BUILTIN_ADAPTERS_R2
        .iter()
        .map(|s| s.as_bytes().to_vec())
        .collect();
    // Direct table (SE, or PE with pe_overlap disabled): user table falling
    // back to the builtin Illumina table. The PE-overlap fallback uses user
    // tables only (fastp does not auto-scan known sequences in PE mode).
    let direct1: &[Vec<u8>] = if user1.is_empty() { &builtin1 } else { &user1 };
    let direct2: &[Vec<u8>] = if user2.is_empty() { &builtin2 } else { &user2 };
    let adapter_enabled = params.adapter_trim;

    // --- outputs ---
    let out_dir = &params.out_dir;
    let stem = stem_of(&params.r1[0]);
    let stem2 = if paired {
        stem_of(&params.r2[0])
    } else {
        stem.clone()
    };
    let ext = match params.out_format {
        OutFormat::Fqgz => "fq.gz",
        OutFormat::Bfq => "bfq",
    };
    // Reference naming: SE = <stem1>.clean; PE = <stem1>.clean_R1 +
    // <stem2>.clean_R2 (R2 uses ITS OWN file's stem).
    let r1_name = if paired {
        format!("{stem}.clean_R1.{ext}")
    } else {
        format!("{stem}.clean.{ext}")
    };
    let mut sink1 = OutSink::create(params.out_format, &out_dir.join(r1_name))?;
    let mut sink2 = if paired {
        Some(OutSink::create(
            params.out_format,
            &out_dir.join(format!("{stem2}.clean_R2.{ext}")),
        )?)
    } else {
        None
    };
    // Reference bug-compat: unpaired files are ALWAYS named `.fq.gz`, but in
    // Bfq mode their CONTENT is still EBFQ (header + raw records, count never
    // patched) — the old writer hardcoded the extension while the encoder
    // followed the format flag.
    let mut sink_u1 = if paired && params.keep_unpaired {
        Some(OutSink::create(
            params.out_format,
            &out_dir.join(format!("{stem}.unpaired_R1.fq.gz")),
        )?)
    } else {
        None
    };
    let mut sink_u2 = if paired && params.keep_unpaired {
        Some(OutSink::create(
            params.out_format,
            &out_dir.join(format!("{stem2}.unpaired_R2.fq.gz")),
        )?)
    } else {
        None
    };

    // --- main pass ---
    let mut s1 = FastqSet::new(&params.r1);
    let mut s2 = if paired {
        Some(FastqSet::new(&params.r2))
    } else {
        None
    };
    // SE adapter auto-detection (spec: "SE adapter auto-detection"). Default off; when
    // off, `backlog` stays empty and `se_table` is the resolved table, so
    // output is byte-identical to the feature being absent.
    let detection_ran =
        params.detect_adapter_se && adapter_enabled && !paired && user1.is_empty();
    let mut backlog: std::collections::VecDeque<Record> = std::collections::VecDeque::new();
    let mut adapter_source: Option<&'static str> = None;
    let mut detected_adapter: Option<String> = None;
    let mut detected_table: Vec<Vec<u8>> = Vec::new();
    if params.detect_adapter_se {
        adapter_source = if !adapter_enabled || paired {
            Some("none")
        } else if !user1.is_empty() {
            Some("table")
        } else {
            // Buffer up to 5000 records from the LIVE stream (no re-read);
            // they re-enter the main loop through `backlog`.
            for _ in 0..5000 {
                match s1.next()? {
                    Some(r) => backlog.push_back(r),
                    None => break,
                }
            }
            let seqs: Vec<&[u8]> = backlog.iter().map(|r| r.seq.as_slice()).collect();
            let table_hits = seqs
                .iter()
                .filter(|s| direct1.iter().any(|ad| match_adapter(s, ad).is_some()))
                .count();
            if table_hits * 100 >= seqs.len() {
                Some("table")
            } else if let Some(cand) = detect::detect_adapter(&seqs) {
                detected_adapter = Some(String::from_utf8_lossy(&cand).into_owned());
                detected_table = vec![cand];
                Some("detected")
            } else {
                Some("none")
            }
        };
    }
    let se_table: &[Vec<u8>] = match adapter_source {
        Some("detected") => &detected_table,
        Some("none") if detection_ran => &[],
        _ => direct1,
    };
let mut before = Stats::default();
    let mut after = Stats::default();
    let mut cnt = Counters::default();

    // Adapter tables are attached per-end below; the PE overlap path
    // overrides them, so `adapters` starts empty here.
    let opts1 = EndOpts {
        front: params.trim_front1,
        tail: params.trim_tail1,
        qtrim: params.qtrim,
        qtrim_cutoff: params.qtrim_cutoff,
        polyg: polyg_on,
        adapters: &[],
    };
    let opts2 = EndOpts {
        front: params.trim_front2,
        tail: params.trim_tail2,
        qtrim: params.qtrim,
        qtrim_cutoff: params.qtrim_cutoff,
        polyg: polyg_on,
        adapters: &[],
    };

    let mut chunk_pos = 0usize;
    loop {
        if !paired {
            let next = match backlog.pop_front() {
                Some(r) => Some(r),
                None => s1.next()?,
            };
            let mut rec = match next {
                Some(r) => r,
                None => break,
            };
            before.add(&rec.seq, &rec.qual);
            let mut flag = process_end(&mut rec, &opts1, &mut cnt);
            if adapter_enabled {
                flag |= adapter_stage(&mut rec, se_table, &mut cnt);
            }
            match filter_end(
                &rec.seq,
                &rec.qual,
                flag,
                params.min_len,
                params.n_max,
                params.q15_frac_max,
            ) {
                Some(reason) => cnt.fail(reason),
                None => {
                    sink1.write(&rec)?;
                    after.add(&rec.seq, &rec.qual);
                }
            }
        } else {
            let a = s1.next()?;
            let b = match s2.as_mut() {
                Some(s) => s.next()?,
                None => None,
            };
            let (mut r1, mut r2) = match (a, b) {
                (None, None) => break,
                (Some(a), Some(b)) => (a, b),
                (Some(_), None) | (None, Some(_)) => return Err(QcError::PairLengthMismatch),
            };
            if !paired_names_match(&r1.name, &r2.name) {
                return Err(QcError::NameMismatch {
                    r1: String::from_utf8_lossy(&r1.name).into_owned(),
                    r2: String::from_utf8_lossy(&r2.name).into_owned(),
                });
            }
            before.add(&r1.seq, &r1.qual);
            before.add(&r2.seq, &r2.qual);

            let mut flag1 = process_end(&mut r1, &opts1, &mut cnt);
            let mut flag2 = process_end(&mut r2, &opts2, &mut cnt);

            if adapter_enabled {
                if params.pe_overlap {
                    match find_overlap(&r1.seq, &r2.seq) {
                        Some((o, ol)) => {
                            // Only read-through (negative offset) trims: both
                            // ends cut to the overlap length. A positive
                            // offset is a reliable overlap, but nothing to
                            // trim — and no table fallback either.
                            if o < 0 {
                                // Per-end counting (reference semantics): each
                                // end counts +1 only if IT lost bases; ends
                                // already shorter than ol lose nothing.
                                let cut1 = r1.seq.len() - ol;
                                let cut2 = r2.seq.len() - ol;
                                r1.seq.truncate(ol);
                                r1.qual.truncate(ol);
                                r2.seq.truncate(ol);
                                r2.qual.truncate(ol);
                                cnt.adapter_reads += u64::from(cut1 > 0) + u64::from(cut2 > 0);
                                cnt.adapter_bases += (cut1 + cut2) as u64;
                                flag1 |= cut1 > 0;
                                flag2 |= cut2 > 0;
                            }
                        }
                        // No reliable overlap: fall back to the known tables
                        // only when the user supplied them explicitly.
                        None => {
                            if !user1.is_empty() {
                                flag1 |= adapter_stage(&mut r1, &user1, &mut cnt);
                            }
                            if !user2.is_empty() {
                                flag2 |= adapter_stage(&mut r2, &user2, &mut cnt);
                            }
                        }
                    }
                } else {
                    flag1 |= adapter_stage(&mut r1, direct1, &mut cnt);
                    flag2 |= adapter_stage(&mut r2, direct2, &mut cnt);
                }
            }

            let f1 = filter_end(
                &r1.seq,
                &r1.qual,
                flag1,
                params.min_len,
                params.n_max,
                params.q15_frac_max,
            );
            let f2 = filter_end(
                &r2.seq,
                &r2.qual,
                flag2,
                params.min_len,
                params.n_max,
                params.q15_frac_max,
            );
            match (f1, f2) {
                (None, None) => {
                    sink1.write(&r1)?;
                    if let Some(s) = sink2.as_mut() {
                        s.write(&r2)?;
                    }
                    after.add(&r1.seq, &r1.qual);
                    after.add(&r2.seq, &r2.qual);
                }
                (Some(ra), Some(rb)) => {
                    cnt.fail(ra);
                    cnt.fail(rb);
                }
                (Some(ra), None) => {
                    cnt.fail(ra);
                    if let Some(s) = sink_u2.as_mut() {
                        s.write(&r2)?;
                        cnt.unpaired_written += 1;
                        after.add(&r2.seq, &r2.qual);
                    }
                }
                (None, Some(rb)) => {
                    cnt.fail(rb);
                    if let Some(s) = sink_u1.as_mut() {
                        s.write(&r1)?;
                        cnt.unpaired_written += 1;
                        after.add(&r1.seq, &r1.qual);
                    }
                }
            }
        }

        chunk_pos += 1;
        if chunk_pos == CHUNK_RECORDS {
            sink1.end_chunk()?;
            if let Some(s) = sink2.as_mut() {
                s.end_chunk()?;
            }
            if let Some(s) = sink_u1.as_mut() {
                s.end_chunk()?;
            }
            if let Some(s) = sink_u2.as_mut() {
                s.end_chunk()?;
            }
            chunk_pos = 0;
        }
    }

    // --- close outputs (bfq read counts patched here) ---
    // Reference-tool bug-compat: only the R1/SE sink is patched, with the
    // TOTAL kept-read count (both mates in PE); the R2/unpaired sinks keep
    // their placeholder 0 (the old writer moved them away before patching).
    let total_kept = after.reads;
    sink1.finish_with_count(total_kept)?;
    if let Some(s) = sink2 {
        s.finish()?;
    }
    if let Some(s) = sink_u1 {
        s.finish()?;
    }
    if let Some(s) = sink_u2 {
        s.finish()?;
    }

    // --- report ---
    let report = QcReport {
        esperanto_qc_version: env!("CARGO_PKG_VERSION").to_string(),
        params: ParamsReport {
            adapter_trim: params.adapter_trim,
            pe_overlap: params.pe_overlap,
            qtrim: params.qtrim,
            polyg: match params.polyg {
                PolygMode::Auto => "auto".to_string(),
                PolygMode::On => "on".to_string(),
                PolygMode::Off => "off".to_string(),
            },
            min_len: params.min_len,
            n_max: params.n_max,
            detect_adapter_se: params.detect_adapter_se.then_some(true),
            q15_frac_max: params.q15_frac_max,
        },
        input: InputReport {
            r1: params.r1.iter().map(|p| p.display().to_string()).collect(),
            r2: params.r2.iter().map(|p| p.display().to_string()).collect(),
            instrument_polyg: instrument.is_some(),
        },
        summary: summary_of(&before, &after),
        filter_reasons: FilterReasons {
            low_quality: cnt.low_quality,
            too_many_n: cnt.too_many_n,
            too_short: cnt.too_short,
            adapter_only: cnt.adapter_only,
        },
        trimming: Trimming {
            adapter_reads: cnt.adapter_reads,
            adapter_bases: cnt.adapter_bases,
            polyg_reads: cnt.polyg_reads,
            polyg_bases: cnt.polyg_bases,
            qtrim_bases: cnt.qtrim_bases,
        },
        per_cycle: cycle_entries(&after),
        per_cycle_before: cycle_entries(&before),
        summary_before: summary_before_of(&before),
        elapsed_seconds: t0.elapsed().as_secs_f64(),
        adapter_source,
        detected_adapter_se: detected_adapter,
    };

    let json = serde_json::to_string_pretty(&report)?;
    std::fs::write(out_dir.join("qc.json"), &json)?;
    std::fs::write(out_dir.join("qc.html"), render_html(&json, &before, &after))?;
    Ok(())
}
