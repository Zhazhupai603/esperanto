//! esperanto-pile — 8-dimensional pileup feature extraction from BAM files.
//!
//! Crate responsibility
//! --------------------
//! Extract per-site pileup features from an indexed BAM, replicating the exact
//! semantics of `pysam.AlignmentFile.pileup(stepper="nofilter")` with default
//! settings (feature_spec.json v1). This includes two behaviors that live
//! inside htslib's `bam_plp` engine rather than in pysam itself:
//!
//! * the PE overlap quality tweak (`overlap_push` / `tweak_overlap_quality`),
//!   and
//! * the `maxcnt = 8000` silent read drop.
//!
//! The extracted features feed a frozen ML model whose scoring stage compares
//! against a golden reference at `rtol = 0`; every constant, ordering rule and
//! timing decision below is contractual and must not be changed.
//!
//! Public API
//! ----------
//! * [`N_FEATURES`] / [`FEATURE_NAMES`] — feature vector length and names:
//!   `depth`, `A/C/G/T counts`, `mean_base_quality`, `strand_bias`,
//!   `mean_mapq`.
//! * [`extract_pileup_features`] — features for one genomic site.
//! * [`extract_pileup_features_batch`] — features for many sites; sites within
//!   `MERGE_GAP` are served from a single region fetch via a sweep-line over
//!   streamed records. Results are bit-identical to calling the single-site
//!   API on the same sites.
//!
//! Invariants
//! ----------
//! * A "column" for site `pos0` contains exactly the reads pushed in file
//!   order that are live and whose `beg <= pos0` at the moment the engine
//!   emits column `pos0` (push order preserved).
//! * The push/emit interleaving of htslib's `bam_plp64_auto` is reproduced:
//!   each auto step calls `next()` first and only enters the push loop when no
//!   column can be emitted (`max_pos <= pos`); EOF is signalled by an explicit
//!   end-of-stream push; after EOF only `next()` drains the buffer.
//! * A read is silently dropped at push time when the engine position equals
//!   the read start and `alive + 1 > MAXCNT` (htslib compares its pool count,
//!   sentinel included, against `maxcnt`; this is the same predicate).
//! * A read whose reference end is `<= pos` at push time never enters the
//!   buffer (but still updates `max_tid`/`max_pos`).
//! * Overlap tweaking only ever modifies the two records of a paired hit, via
//!   quality overlays; untouched records keep their original qualities.
//! * Determinism: identical inputs produce byte-identical outputs. All
//!   iteration that influences output is over vectors in file order; the
//!   waiting table is only used for point lookups.
//!
//! Numerical conventions
//! ---------------------
//! * Integer sums (qualities, MAPQs) are exact; means are computed as
//!   `(sum as f64 / depth as f64) as f32`.
//! * The mismatch/equal-qual quality penalty `0.8 * q` is evaluated in f64 and
//!   truncated toward zero on conversion to `u8`, matching the C semantics of
//!   `uint8_t = 0.8 * qual`.

#![deny(unsafe_code)]

use std::collections::{HashMap, VecDeque};

use rust_htslib::bam::record::{Cigar, Record};
use rust_htslib::bam::{IndexedReader, Read as _};
use thiserror::Error;

/// Number of features per site.
pub const N_FEATURES: usize = 8;

/// Feature names, in output order.
pub const FEATURE_NAMES: [&str; N_FEATURES] = [
    "depth",
    "A_count",
    "C_count",
    "G_count",
    "T_count",
    "mean_base_quality",
    "strand_bias",
    "mean_mapq",
];

/// Minimum (possibly overlap-tweaked) base quality for a read to be counted.
pub const MIN_BASE_QUALITY: u8 = 13;

/// htslib pileup engine default maximum depth; excess reads starting at the
/// engine position are silently dropped at push time.
pub const MAXCNT: usize = 8000;

/// Batch mode: consecutive sites (same reference, ascending position) whose
/// distance is at most this are served from one region fetch.
pub const MERGE_GAP: i64 = 2000;

/// Errors produced by this crate.
#[derive(Debug, Error)]
pub enum PileError {
    /// Reference name not present in the BAM header.
    #[error("contig not found in BAM header: {0}")]
    ContigNotFound(String),

    /// Site position is not a valid 1-based coordinate (< 1).
    #[error("invalid 1-based position: {0}")]
    InvalidPosition(i64),

    /// Region fetch failed (missing/corrupt index, I/O error).
    #[error("failed to fetch region on tid {tid} at {start}-{end}: {source}")]
    Fetch {
        tid: i32,
        start: i64,
        end: i64,
        source: rust_htslib::errors::Error,
    },

    /// Record iteration failed (corrupt BAM stream).
    #[error("failed to read BAM record: {0}")]
    Read(#[from] rust_htslib::errors::Error),

    /// I/O failure (includes the `.baln` channel).
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
}

/// True when the record carries the collapsed-rescue provenance tag
/// (`RE:Z:collapsed`): a repeat-family placement whose bases are
/// alphabet-ambiguous (A==G, T==C). These records are excluded from pileup
/// feature extraction — consistent with the scan evidence rule.
fn is_collapsed(rec: &Record) -> bool {
    matches!(rec.aux(b"RE"), Ok(rust_htslib::bam::record::Aux::String(v)) if v == "collapsed")
}

/// Extract pileup features at a single site.
///
/// * `bam` — indexed BAM reader (a `.bai`/`.csi` index must exist).
/// * `chrom` — reference name as in the BAM header.
/// * `pos_1based` — 1-based site position.
///
/// Returns `[depth, A, C, G, T, mean_base_quality, strand_bias, mean_mapq]`,
/// the all-zero vector when no qualifying read covers the site. Errors with
/// [`PileError::ContigNotFound`] when `chrom` is absent from the header.
pub fn extract_pileup_features(
    bam: &mut IndexedReader,
    chrom: &str,
    pos_1based: i64,
) -> Result<[f32; N_FEATURES], PileError> {
    let (tid, pos0) = resolve_site(bam, chrom, pos_1based)?;
    let records = collect_site_records(bam, tid, pos0)?;
    Ok(run_site(pos0, records))
}

/// Extract pileup features for many sites in one pass.
///
/// Sites are sorted by `(reference, position)`; runs of consecutive sites on
/// the same reference with gap `<= MERGE_GAP` are served from a single region
/// fetch using a sweep-line over the streamed records. For every site the
/// record set handed to the pileup engine is `{pos <= pos0 && end > pos0}` in
/// file order — exactly the set the per-site fetch delivers — so the results
/// are bit-identical to [`extract_pileup_features`] called per site. Output
/// order follows the input `sites` order.
pub fn extract_pileup_features_batch(
    bam: &mut IndexedReader,
    sites: &[(&str, i64)],
) -> Result<Vec<[f32; N_FEATURES]>, PileError> {
    if sites.is_empty() {
        return Ok(Vec::new());
    }

    // Resolve all sites up front so failures happen before any region I/O.
    let mut order: Vec<(i32, i64, usize)> = Vec::with_capacity(sites.len());
    for (i, (chrom, pos_1based)) in sites.iter().enumerate() {
        let (tid, pos0) = resolve_site(bam, chrom, *pos_1based)?;
        order.push((tid, pos0, i));
    }
    // Stable sort: duplicates keep input order; results are identical anyway.
    order.sort_by_key(|&(tid, pos0, _)| (tid, pos0));

    let mut out: Vec<[f32; N_FEATURES]> = vec![[0.0; N_FEATURES]; sites.len()];

    let mut group_start = 0usize;
    while group_start < order.len() {
        let (g_tid, g_first_pos0, _) = order[group_start];
        let mut group_end = group_start + 1;
        while group_end < order.len() {
            let &(tid, pos0, _) = &order[group_end];
            let &(_, prev_pos0, _) = &order[group_end - 1];
            if tid != g_tid || pos0 - prev_pos0 > MERGE_GAP {
                break;
            }
            group_end += 1;
        }

        let last_pos0 = order[group_end - 1].1;
        sweep_group_bam(
            bam,
            g_tid,
            g_first_pos0,
            last_pos0,
            &order[group_start..group_end],
            &mut out,
        )?;

        group_start = group_end;
    }

    Ok(out)
}

/// Streaming pileup over one site group: one engine instance walks the
/// record stream once, emitting columns in position order; every matching
/// site takes its features from the emitted column. Collapsed-rescue
/// records never enter the engine. Returns true when a MAXCNT drop
/// occurred — streaming and per-site modes can diverge at saturated
/// columns, so the caller must fall back to the per-site sweep.
fn stream_core<R: PileRecord>(
    records: impl Iterator<Item = Result<R, PileError>>,
    group: &[(i32, i64, usize)],
    out: &mut [[f32; N_FEATURES]],
) -> Result<bool, PileError> {
    let mut plp: EventPlp<R> = EventPlp::new();
    let mut records = records.peekable();
    for &(tid, pos0, orig_idx) in group {
        // Feed every record starting at or before the site (stream is
        // (tid,pos)-sorted, same order the per-column engine consumes).
        loop {
            let take = match records.peek() {
                Some(Ok(r)) => r.tid() < tid || (r.tid() == tid && r.pos() <= pos0),
                Some(Err(_)) => true,
                None => false,
            };
            if !take {
                break;
            }
            let r = records.next().expect("peeked")?;
            plp.retire_until(r.tid(), r.pos());
            if !r.is_collapsed() {
                plp.push(r);
            }
        }
        plp.retire_until(tid, pos0);
        out[orig_idx] = plp.features_at(pos0);
        if plp.dropped_maxcnt {
            return Ok(true);
        }
    }
    Ok(plp.dropped_maxcnt)
}

/// BAM sweep for one site group: single streaming pass; on a MAXCNT drop
/// (rare repeat-region saturation) the region is re-fetched and the group
/// is redone with the exact per-site semantics.
fn sweep_group_bam(
    bam: &mut IndexedReader,
    tid: i32,
    first_pos0: i64,
    last_pos0: i64,
    group: &[(i32, i64, usize)],
    out: &mut [[f32; N_FEATURES]],
) -> Result<(), PileError> {
    let fetch = |bam: &mut IndexedReader| {
        bam.fetch((tid, first_pos0, last_pos0 + 1))
            .map_err(|source| PileError::Fetch {
                tid,
                start: first_pos0,
                end: last_pos0 + 1,
                source,
            })
    };
    fetch(bam)?;
    let dropped = stream_core(
        bam.records().map(|r| r.map_err(PileError::Read)),
        group,
        out,
    )?;
    if dropped {
        fetch(bam)?;
        let records = bam
            .records()
            .map(|r| r.map_err(PileError::Read))
            .collect::<Result<Vec<_>, _>>()?;
        legacy_sweep_group(records, group, out)?;
    }
    Ok(())
}

/// `.baln` sweep for one site group: same structure as
/// [`sweep_group_bam`]; the fallback re-reads the window.
fn sweep_group_baln(
    index: &esperanto_bamio::baln::BalnIndex,
    file: &std::fs::File,
    tid: i32,
    first_pos0: i64,
    last_pos0: i64,
    group: &[(i32, i64, usize)],
    out: &mut [[f32; N_FEATURES]],
) -> Result<(), PileError> {
    let raw = baln_window_records(index, file, tid, first_pos0, last_pos0)?;
    let dropped = stream_core(raw.into_iter().map(Ok), group, out)?;
    if dropped {
        let raw = baln_window_records(index, file, tid, first_pos0, last_pos0)?;
        legacy_sweep_group(raw, group, out)?;
    }
    Ok(())
}

/// Records overlapping `[first_pos0, last_pos0]` from the `.baln` channel
/// (index window `pos < last+1 && pos+span > first`, htslib fetch overlap
/// semantics).
fn baln_window_records(
    index: &esperanto_bamio::baln::BalnIndex,
    file: &std::fs::File,
    tid: i32,
    first_pos0: i64,
    last_pos0: i64,
) -> Result<Vec<BalnPileRecord>, PileError> {
    let lo_pos = first_pos0 - index.max_span;
    let lo = index
        .idx
        .partition_point(|e| (e.0, e.1) < (tid, lo_pos));
    let hi = index
        .idx
        .partition_point(|e| (e.0, e.1) < (tid, last_pos0 + 1));
    let mut raw: Vec<BalnPileRecord> = Vec::new();
    for &(t, pos, off, span) in &index.idx[lo..hi] {
        if t != tid || pos + span <= first_pos0 {
            continue;
        }
        let Some(rec) = esperanto_bamio::baln::read_record_at(file, off)? else {
            continue;
        };
        raw.push(BalnPileRecord::new(rec)?);
    }
    Ok(raw)
}

fn legacy_sweep_group<R: PileRecord>(
    records: impl IntoIterator<Item = R>,
    group: &[(i32, i64, usize)],
    out: &mut [[f32; N_FEATURES]],
) -> Result<(), PileError> {
    let mut active: Vec<R> = Vec::new();
    let mut pending: Option<R> = None;
    let mut records = records.into_iter();
    for &(_tid, pos0, orig_idx) in group {
        loop {
            let rec = match pending.take() {
                Some(r) => r,
                None => match records.next() {
                    Some(r) => r,
                    None => break,
                },
            };
            if rec.pos() <= pos0 {
                if !rec.is_collapsed() {
                    active.push(rec);
                }
            } else {
                pending = Some(rec);
                break;
            }
        }
        active.retain(|r| r.pos() + r.ref_len() > pos0);
        out[orig_idx] = run_site(pos0, active.iter().cloned());
    }
    Ok(())
}

/// Extract pileup features for many sites from the `.baln` fast channel.
///
/// Same grouping, sweep and engine as the BAM path; the coordinate index
/// selects exactly the records htslib fetch would deliver, so features are
/// bit-identical to [`extract_pileup_features_batch`]. Output order follows
/// the input `sites` order.
pub fn extract_pileup_features_batch_baln(
    index: &esperanto_bamio::baln::BalnIndex,
    file: &std::fs::File,
    sites: &[(&str, i64)],
) -> Result<Vec<[f32; N_FEATURES]>, PileError> {
    if sites.is_empty() {
        return Ok(Vec::new());
    }
    let tid_of = |chrom: &str| -> Result<i32, PileError> {
        index
            .contigs
            .iter()
            .position(|c| c == chrom)
            .map(|i| i as i32)
            .ok_or_else(|| PileError::ContigNotFound(chrom.to_string()))
    };
    let mut order: Vec<(i32, i64, usize)> = Vec::with_capacity(sites.len());
    for (i, (chrom, pos_1based)) in sites.iter().enumerate() {
        if *pos_1based < 1 {
            return Err(PileError::InvalidPosition(*pos_1based));
        }
        order.push((tid_of(chrom)?, pos_1based - 1, i));
    }
    order.sort_by_key(|&(tid, pos0, _)| (tid, pos0));

    let mut out: Vec<[f32; N_FEATURES]> = vec![[0.0; N_FEATURES]; sites.len()];

    let mut group_start = 0usize;
    while group_start < order.len() {
        let (g_tid, g_first_pos0, _) = order[group_start];
        let mut group_end = group_start + 1;
        while group_end < order.len() {
            let &(tid, pos0, _) = &order[group_end];
            let &(_, prev_pos0, _) = &order[group_end - 1];
            if tid != g_tid || pos0 - prev_pos0 > MERGE_GAP {
                break;
            }
            group_end += 1;
        }
        let last_pos0 = order[group_end - 1].1;
        sweep_group_baln(
            index,
            file,
            g_tid,
            g_first_pos0,
            last_pos0,
            &order[group_start..group_end],
            &mut out,
        )?;

        group_start = group_end;
    }
    Ok(out)
}

/// Resolve a (chrom, 1-based pos) site into (tid, 0-based pos).
fn resolve_site(
    bam: &mut IndexedReader,
    chrom: &str,
    pos_1based: i64,
) -> Result<(i32, i64), PileError> {
    let tid = bam
        .header()
        .tid(chrom.as_bytes())
        .map(|t| t as i32)
        .ok_or_else(|| PileError::ContigNotFound(chrom.to_string()))?;
    if pos_1based < 1 {
        return Err(PileError::InvalidPosition(pos_1based));
    }
    Ok((tid, pos_1based - 1))
}

/// Fetch all records overlapping the single base at `pos0`, in file order.
fn collect_site_records(
    bam: &mut IndexedReader,
    tid: i32,
    pos0: i64,
) -> Result<Vec<Record>, PileError> {
    bam.fetch((tid, pos0, pos0 + 1))
        .map_err(|source| PileError::Fetch {
            tid,
            start: pos0,
            end: pos0 + 1,
            source,
        })?;
    let mut records = Vec::new();
    for rec in bam.records() {
        let rec = rec?;
        if !is_collapsed(&rec) {
            records.push(rec);
        }
    }
    Ok(records)
}

/// Record abstraction for the pileup engine: implemented by htslib
/// `Record` (BAM source) and by [`BalnPileRecord`] (`.baln` source), so the
/// engine runs without re-encoding either representation.
pub trait PileRecord: Clone {
    /// Reference id (-1 = unmapped).
    fn tid(&self) -> i32;
    /// Mate reference id.
    fn mtid(&self) -> i32;
    /// 0-based leftmost position.
    fn pos(&self) -> i64;
    /// Mate position.
    fn mpos(&self) -> i64;
    /// Template length (signed).
    fn insert_size(&self) -> i64;
    /// SAM flags.
    fn flags(&self) -> u16;
    /// Mapping quality.
    fn mapq(&self) -> u8;
    /// Read name (no trailing NUL).
    fn qname(&self) -> &[u8];
    /// Sequence length.
    fn seq_len(&self) -> usize;
    /// Raw phred qualities (BAM-disk convention).
    fn qual(&self) -> &[u8];
    /// ASCII base at query position `i`.
    fn base_at(&self, i: usize) -> u8;
    /// CIGAR op count.
    fn cigar_len(&self) -> usize;
    /// CIGAR op at `i`.
    fn cigar_at(&self, i: usize) -> Cigar;
    /// Reference span (M/D/N/=/X consumption).
    fn ref_len(&self) -> i64;
    /// True for collapsed-rescue records (alphabet-ambiguous bases).
    fn is_collapsed(&self) -> bool;
}

impl PileRecord for Record {
    fn tid(&self) -> i32 { self.tid() }
    fn mtid(&self) -> i32 { self.mtid() }
    fn pos(&self) -> i64 { self.pos() }
    fn mpos(&self) -> i64 { self.mpos() }
    fn insert_size(&self) -> i64 { self.insert_size() }
    fn flags(&self) -> u16 { self.flags() }
    fn mapq(&self) -> u8 { self.mapq() }
    fn qname(&self) -> &[u8] { self.qname() }
    fn seq_len(&self) -> usize { self.seq_len() }
    fn qual(&self) -> &[u8] { self.qual() }
    fn base_at(&self, i: usize) -> u8 { self.seq()[i] }
    fn cigar_len(&self) -> usize { self.cigar_len() }
    fn cigar_at(&self, i: usize) -> Cigar { self.cigar()[i] }
    fn ref_len(&self) -> i64 { cigar_rlen(self) }
    fn is_collapsed(&self) -> bool { is_collapsed(self) }
}

/// `.baln` record prepared for the pileup engine: fields plus the decoded
/// CIGAR (computed once at conversion; no sequence re-encoding).
#[derive(Clone)]
pub struct BalnPileRecord {
    rec: esperanto_bamio::baln::BalnRecord,
    cigar: Vec<Cigar>,
}

impl BalnPileRecord {
    /// Convert a parsed `.baln` record (cheap: CIGAR words to enums only).
    pub fn new(rec: esperanto_bamio::baln::BalnRecord) -> Result<BalnPileRecord, PileError> {
        let mut cigar = Vec::with_capacity(rec.cigar.len());
        for &c in &rec.cigar {
            let len = c >> 4;
            cigar.push(match c & 0xF {
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
                    return Err(PileError::Io(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("baln record: unknown cigar op code {other}"),
                    )))
                }
            });
        }
        Ok(BalnPileRecord { rec, cigar })
    }
}

impl PileRecord for BalnPileRecord {
    fn tid(&self) -> i32 { self.rec.tid }
    fn mtid(&self) -> i32 { self.rec.mtid }
    fn pos(&self) -> i64 { self.rec.pos }
    fn mpos(&self) -> i64 { self.rec.mpos }
    fn insert_size(&self) -> i64 { self.rec.isize }
    fn flags(&self) -> u16 { self.rec.flag }
    fn mapq(&self) -> u8 { self.rec.mapq }
    fn qname(&self) -> &[u8] { &self.rec.name }
    fn seq_len(&self) -> usize { self.rec.l_seq }
    fn qual(&self) -> &[u8] { &self.rec.qual }
    fn base_at(&self, i: usize) -> u8 { self.rec.seq_ascii[i] }
    fn cigar_len(&self) -> usize { self.cigar.len() }
    fn cigar_at(&self, i: usize) -> Cigar { self.cigar[i] }
    fn ref_len(&self) -> i64 {
        self.cigar
            .iter()
            .map(|c| match c {
                Cigar::Match(n) | Cigar::Del(n) | Cigar::RefSkip(n) | Cigar::Equal(n)
                | Cigar::Diff(n) => *n as i64,
                _ => 0,
            })
            .sum()
    }
    fn is_collapsed(&self) -> bool { self.rec.re.as_deref() == Some("collapsed") }
}

/// Run the pileup engine over one site's record sequence and extract features.
fn run_site<R: PileRecord>(pos0: i64, records: impl IntoIterator<Item = R>) -> [f32; N_FEATURES] {
    let mut plp: Plp<R> = Plp::new();
    let mut it = records.into_iter();
    loop {
        match plp.next() {
            Some(col) => {
                if col.pos == pos0 {
                    return features_from_column(&plp, &col);
                }
                // Columns are emitted at strictly increasing positions;
                // passing the target means it had no coverage.
                if col.pos > pos0 {
                    return [0.0; N_FEATURES];
                }
            }
            None => {
                if plp.is_eof && plp.order.is_empty() {
                    return [0.0; N_FEATURES];
                }
                match it.next() {
                    Some(r) => plp.push(r),
                    None => plp.is_eof = true,
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Event-driven engine (fast path)
// ---------------------------------------------------------------------------

/// Event-driven pileup engine: records stream in (tid,pos) order; nodes
/// retire through an end-keyed heap (O(log n) per event) and columns are
/// built ONLY at site positions (O(depth) per site, not per position).
/// Semantics equal the per-column engine outside MAXCNT saturation:
/// pairing/tweaks happen at push time identically, and the live set at a
/// site column is `start <= pos && end > pos` in both engines. Feature
/// arithmetic is integer accumulation, so entry order is irrelevant.
/// Any buffer-full push sets `dropped_maxcnt` (a superset of the exact
/// drop rule); the caller must then redo the group with the exact
/// per-site engine.
struct EventPlp<R: PileRecord> {
    nodes: Vec<Option<Node<R>>>,
    free: Vec<usize>,
    /// (Reverse(tid), Reverse(end), idx) min-heap for retirement.
    by_end: std::collections::BinaryHeap<std::cmp::Reverse<(i32, i64, usize)>>,
    /// qname -> node waiting for its mate (same table as Plp).
    waiting: HashMap<Vec<u8>, usize>,
    live: usize,
    tid: i32,
    max_tid: i32,
    max_pos: i64,
    dropped_maxcnt: bool,
}

impl<R: PileRecord> EventPlp<R> {
    fn new() -> Self {
        EventPlp {
            nodes: Vec::new(),
            free: Vec::new(),
            by_end: std::collections::BinaryHeap::new(),
            waiting: HashMap::new(),
            live: 0,
            tid: -1,
            max_tid: -1,
            max_pos: -1,
            dropped_maxcnt: false,
        }
    }

    /// Retire every node whose (tid, end) passed `pos` on `tid`.
    fn retire_until(&mut self, tid: i32, pos: i64) {
        loop {
            let should_retire = match self.by_end.peek() {
                Some(std::cmp::Reverse((t, e, _))) => {
                    *t < tid || (*t == tid && *e <= pos)
                }
                None => false,
            };
            if !should_retire {
                break;
            }
            let std::cmp::Reverse((_, _, idx)) = self.by_end.pop().expect("peeked");
            if let Some(node) = self.nodes[idx].take() {
                overlap_remove(&mut self.waiting, node.rec.qname());
                self.free.push(idx);
                self.live -= 1;
            }
        }
    }

    /// `bam_plp_push` minus the per-column bookkeeping.
    fn push(&mut self, rec: R) {
        if rec.tid() < 0 || rec.flags() & FLAG_UNMAP != 0 {
            overlap_remove(&mut self.waiting, rec.qname());
            return;
        }
        // Superset of the exact MAXCNT drop rule (never misses a saturated
        // push): any buffer-full push flags the group for the exact
        // per-site redo.
        if self.live + 1 > MAXCNT {
            self.dropped_maxcnt = true;
        }
        let beg = rec.pos();
        let end = beg + rec.ref_len();
        self.max_tid = rec.tid();
        self.max_pos = beg;
        let idx = match self.free.pop() {
            Some(i) => i,
            None => {
                self.nodes.push(None);
                self.nodes.len() - 1
            }
        };
        let tid = rec.tid();
        self.tid = tid;
        self.by_end
            .push(std::cmp::Reverse((tid, end, idx)));
        self.nodes[idx] = Some(Node {
            tid,
            beg,
            end,
            rec,
            qual: None,
            cstate: CigarState { k: -1, x: 0, y: 0 },
        });
        self.live += 1;
        self.overlap_push(idx);
    }

    /// Pairing identical to Plp::overlap_push.
    fn overlap_push(&mut self, idx: usize) {
        let (rec_flags, rec_mtid, rec_tid, rec_isize, rec_mpos, rec_seq_len, rec_qname, node_end) = {
            let node = match self.nodes[idx].as_ref() {
                Some(n) => n,
                None => return,
            };
            (
                node.rec.flags(),
                node.rec.mtid(),
                node.rec.tid(),
                node.rec.insert_size(),
                node.rec.mpos(),
                node.rec.seq_len(),
                node.rec.qname().to_vec(),
                node.end,
            )
        };
        if rec_flags & FLAG_MATE_UNMAP != 0 || rec_flags & FLAG_PROPER_PAIR == 0 {
            return;
        }
        if (rec_mtid >= 0 && rec_tid != rec_mtid)
            || (rec_isize.abs() >= 2 * rec_seq_len as i64 && rec_mpos >= node_end)
        {
            return;
        }
        if let Some(&a_idx) = self.waiting.get(&rec_qname) {
            self.waiting.remove(&rec_qname);
            let (a, b) = if a_idx < idx {
                let (x, y) = self.nodes.split_at_mut(idx);
                (&mut x[a_idx], &mut y[0])
            } else {
                let (x, y) = self.nodes.split_at_mut(a_idx);
                (&mut y[0], &mut x[idx])
            };
            if let (Some(an), Some(bn)) = (a.as_mut(), b.as_mut()) {
                tweak_overlap_quality(an, bn);
            }
        } else {
            self.waiting.insert(rec_qname, idx);
        }
    }

    /// Build the feature vector at one site column (live set already
    /// retired to `end > pos0`).
    fn features_at(&mut self, pos0: i64) -> [f32; N_FEATURES] {
        let mut depth: u64 = 0;
        let mut base_counts: [u64; 4] = [0; 4];
        let mut qual_sum: u64 = 0;
        let mut mapq_sum: u64 = 0;
        let mut forward: u64 = 0;
        for slot in self.nodes.iter_mut() {
            let node = match slot.as_mut() {
                Some(n) => n,
                None => continue,
            };
            if node.tid != self.tid || node.beg > pos0 {
                continue;
            }
            let qpos = match resolve_at_site(&node.rec, pos0, &mut node.cstate) {
                Some(q) => q as usize,
                None => continue,
            };
            let rec = &node.rec;
            let seq_len = rec.seq_len();
            if qpos >= seq_len {
                continue;
            }
            let qual_slice = node.qual.as_deref().unwrap_or_else(|| rec.qual());
            let qual = if qpos < qual_slice.len() {
                qual_slice[qpos]
            } else if qpos < seq_len {
                0xFF
            } else {
                0
            };
            if qual < MIN_BASE_QUALITY {
                continue;
            }
            depth += 1;
            match rec.base_at(qpos) {
                b'A' | b'a' => base_counts[0] += 1,
                b'C' | b'c' => base_counts[1] += 1,
                b'G' | b'g' => base_counts[2] += 1,
                b'T' | b't' => base_counts[3] += 1,
                _ => {}
            }
            if qpos < qual_slice.len() {
                qual_sum += u64::from(qual_slice[qpos]);
            }
            mapq_sum += u64::from(rec.mapq());
            if rec.flags() & FLAG_REVERSE == 0 {
                forward += 1;
            }
        }
        if depth == 0 {
            return [0.0; N_FEATURES];
        }
        let n = depth as f32;
        [
            depth as f32,
            base_counts[0] as f32,
            base_counts[1] as f32,
            base_counts[2] as f32,
            base_counts[3] as f32,
            qual_sum as f32 / n,
            forward as f32 / n,
            mapq_sum as f32 / n,
        ]
    }
}

/// Compute the 8 features from one emitted column (entries in push order).
fn features_from_column<R: PileRecord>(plp: &Plp<R>, col: &Column) -> [f32; N_FEATURES] {
    let mut depth: u64 = 0;
    let mut base_counts: [u64; 4] = [0; 4];
    let mut qual_sum: u64 = 0;
    let mut mapq_sum: u64 = 0;
    let mut forward: u64 = 0;

    for entry in &col.entries {
        // Reads on D/N (or otherwise without an aligned query base at this
        // column) never contribute.
        let qpos = match entry.qpos {
            Some(q) => q as usize,
            None => continue,
        };
        let node = match plp.nodes[entry.node].as_ref() {
            Some(n) => n,
            None => continue, // unreachable: entries reference live nodes
        };
        let rec = &node.rec;
        let seq_len = rec.seq_len();
        // Spec quality rule: qpos >= l_qseq -> quality treated as 0 (filtered
        // below); this also guards sequence/quality slice access.
        if qpos >= seq_len {
            continue;
        }
        // Reference quality gate: beyond-l_qseq reads as 0 (filtered); a
        // quality slice shorter than expected reads as 0xFF (kept) for the
        // gate, and contributes 0 to the mean — never panics on odd BAMs.
        let qual_slice = node.qual.as_deref().unwrap_or_else(|| rec.qual());
        let qual = if qpos < qual_slice.len() {
            qual_slice[qpos]
        } else if qpos < seq_len {
            0xFF
        } else {
            0
        };
        if qual < MIN_BASE_QUALITY {
            continue;
        }
        // Reference depth semantics: EVERY read passing the quality gate
        // counts toward depth / mean quality / strand / MAPQ — including N
        // and other non-ACGT codes. Only the four base buckets are ACGT-only.
        depth += 1;
        let base = rec.base_at(qpos);
        match base {
            b'A' | b'a' => base_counts[0] += 1,
            b'C' | b'c' => base_counts[1] += 1,
            b'G' | b'g' => base_counts[2] += 1,
            b'T' | b't' => base_counts[3] += 1,
            _ => {}
        }
        if qpos < qual_slice.len() {
            qual_sum += u64::from(qual_slice[qpos]);
        }
        mapq_sum += u64::from(rec.mapq());
        if rec.flags() & FLAG_REVERSE == 0 {
            forward += 1;
        }
    }

    if depth == 0 {
        return [0.0; N_FEATURES];
    }

    // Reference arithmetic: means are f32 divisions of the integer sums
    // (kept bit-identical to the oracle; a f64 divide + cast can double-round
    // one ulp differently on tie-boundary quotients).
    let n = depth as f32;
    [
        depth as f32,
        base_counts[0] as f32,
        base_counts[1] as f32,
        base_counts[2] as f32,
        base_counts[3] as f32,
        qual_sum as f32 / n,
        forward as f32 / n,
        mapq_sum as f32 / n,
    ]
}

const FLAG_PAIRED: u16 = 0x1;
const FLAG_PROPER_PAIR: u16 = 0x2;
const FLAG_UNMAP: u16 = 0x4;
const FLAG_MATE_UNMAP: u16 = 0x8;
const FLAG_REVERSE: u16 = 0x10;

/// Raw reference length consumed by a CIGAR: sum of M/D/N/=/X lengths
/// (htslib `bam_cigar2rlen`; insertions and clips consume no reference).
fn cigar_rlen(rec: &Record) -> i64 {
    let mut rlen: i64 = 0;
    for c in rec.cigar().iter() {
        match c {
            Cigar::Match(_)
            | Cigar::Del(_)
            | Cigar::RefSkip(_)
            | Cigar::Equal(_)
            | Cigar::Diff(_) => rlen += c.len() as i64,
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    rlen
}

/// Per-node incremental CIGAR cursor (htslib `cstate_t`).
/// k: index of the CIGAR op last processed (-1 = never processed);
/// x/y: reference/query coordinate of the start of op k.
#[derive(Clone, Copy)]
struct CigarState {
    k: i64,
    x: i64,
    y: i64,
}

/// Port of htslib `resolve_cigar2`: resolve the aligned query position of a
/// buffered node at column `pos`, advancing the persistent per-node cursor
/// by at most one ref-consuming op per column. Returns `Some(qpos)` on M/=/X
/// ops, `None` on D/N (deletion / reference skip).
///
/// The incremental stepping is contractual: a zero-length D/N op makes the
/// cursor rest on it for one column (reporting a deletion), a quirk that a
/// random-access "which op contains pos" walk does not reproduce.
fn resolve_cigar2<R: PileRecord>(rec: &R, pos: i64, s: &mut CigarState) -> Option<u32> {
    let n = rec.cigar_len() as i64;
    if s.k == -1 {
        // Find the first M/D/N/=/X op, accumulating query offsets over I/S.
        let mut k: i64 = 0;
        let x = rec.pos();
        let mut y: i64 = 0;
        let mut found = false;
        while k < n {
            match rec.cigar_at(k as usize) {
                Cigar::Match(_) | Cigar::Del(_) | Cigar::RefSkip(_)
                | Cigar::Equal(_) | Cigar::Diff(_) => {
                    found = true;
                    break;
                }
                Cigar::Ins(_) | Cigar::SoftClip(_) => y += rec.cigar_at(k as usize).len() as i64,
                Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
            k += 1;
        }
        if !found {
            // Buffered reads have rlen > 0, so a ref-consuming op exists;
            // tolerate pathological CIGARs by reporting no alignment.
            return None;
        }
        s.k = k;
        s.x = x;
        s.y = y;
    } else {
        let cur = rec.cigar_at(s.k as usize);
        let l = cur.len() as i64;
        if pos - s.x >= l {
            // Advance exactly one ref-consuming op, scanning past I/S/H/P.
            // (Live nodes satisfy pos < end, so the cursor never legitimately
            // runs past the last ref-consuming op.)
            if matches!(cur, Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_)) {
                s.y += l;
            }
            s.x += l;
            let mut k = s.k + 1;
            while k < n {
                match rec.cigar_at(k as usize) {
                    Cigar::Match(_) | Cigar::Del(_) | Cigar::RefSkip(_)
                    | Cigar::Equal(_) | Cigar::Diff(_) => break,
                    Cigar::Ins(_) | Cigar::SoftClip(_) => s.y += rec.cigar_at(k as usize).len() as i64,
                    Cigar::HardClip(_) | Cigar::Pad(_) => {}
                }
                k += 1;
            }
            if k >= n {
                // C asserts here; unreachable for live nodes. Tolerate.
                return None;
            }
            s.k = k;
        }
    }
    match rec.cigar_at(s.k as usize) {
        Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => Some((s.y + (pos - s.x)) as u32),
        // D / N
        _ => None,
    }
}

/// Site-granular resolve: `resolve_cigar2` advances at most one
/// ref-consuming op per call (its per-column contract), so a jump to a far
/// site may need several calls. Returns the final answer once the cursor's
/// current op covers `pos` (or cannot advance further).
fn resolve_at_site<R: PileRecord>(rec: &R, pos: i64, s: &mut CigarState) -> Option<u32> {
    loop {
        let prev = (s.k, s.x);
        let out = resolve_cigar2(rec, pos, s);
        let cur_len = rec.cigar_at(s.k as usize).len() as i64;
        if pos - s.x < cur_len || (s.k, s.x) == prev {
            return out;
        }
    }
}

/// One buffered read. `qual` is the (lazily copied) quality overlay written
/// by the PE overlap tweak; `None` means original qualities.
struct Node<R: PileRecord> {
    tid: i32,
    beg: i64,
    end: i64,
    rec: R,
    qual: Option<Vec<u8>>,
    cstate: CigarState,
}

/// One column member: index into the engine's node slab plus the resolved
/// query position (None for del/refskip).
struct ColumnEntry {
    node: usize,
    qpos: Option<u32>,
}

/// An emitted pileup column.
struct Column {
    pos: i64,
    entries: Vec<ColumnEntry>,
}

/// The `bam_plp` engine: push/emit state machine with the PE overlap waiting
/// table, interleaved exactly per `bam_plp64_auto`.
struct Plp<R: PileRecord> {
    nodes: Vec<Option<Node<R>>>,
    free: Vec<usize>,
    /// Buffered node indices in push order (mirrors the C linked list).
    order: VecDeque<usize>,
    /// qname -> node waiting for its mate (khash olap_hash).
    waiting: HashMap<Vec<u8>, usize>,
    tid: i32,
    pos: i64,
    max_tid: i32,
    max_pos: i64,
    is_eof: bool,
    /// Set when the MAXCNT rule discarded a record (streaming vs per-site
    /// modes can diverge at saturated columns; callers must fall back).
    dropped_maxcnt: bool,
}

impl<R: PileRecord> Plp<R> {
    fn new() -> Self {
        Plp {
            nodes: Vec::new(),
            free: Vec::new(),
            order: VecDeque::new(),
            waiting: HashMap::new(),
            tid: 0,
            pos: 0,
            max_tid: -1,
            max_pos: -1,
            is_eof: false,
            dropped_maxcnt: false,
        }
    }

    /// `bam_plp64_auto`: try next() first; only when no column can come out
    /// (`max_pos <= pos`, pre-EOF) push records one at a time, retrying
    /// next() after each push; at end of stream push EOF once, then drain.
    /// `bam_plp64_next`: emit the next non-empty column, retiring finished
    /// nodes; returns None while the current position may still receive
    /// reads (`max_pos <= pos`, or the drained EOF buffer).
    fn next(&mut self) -> Option<Column> {
        if self.is_eof && self.order.is_empty() {
            return None;
        }
        while self.is_eof
            || self.max_tid > self.tid
            || (self.max_tid == self.tid && self.max_pos > self.pos)
        {
            let col_pos = self.pos;
            let mut entries: Vec<ColumnEntry> = Vec::new();

            // Build the column at self.pos, retiring dead nodes (this also
            // removes them from the overlap waiting table).
            let mut i = 0;
            while i < self.order.len() {
                let idx = self.order[i];
                // Decide retirement, resolve the CIGAR cursor, or drop stale
                // slots — in one mutable pass, in push order.
                let mut retired_qname: Option<Vec<u8>> = None;
                let mut resolved: Option<Option<u32>> = None;
                match self.nodes[idx].as_mut() {
                    Some(node) => {
                        if node.tid < self.tid || (node.tid == self.tid && node.end <= self.pos) {
                            retired_qname = Some(node.rec.qname().to_vec());
                        } else if node.tid == self.tid && node.beg <= self.pos {
                            resolved = Some(resolve_cigar2(&node.rec, self.pos, &mut node.cstate));
                        }
                    }
                    None => {
                        // Unreachable: the order list only holds live slots.
                        self.order.remove(i);
                        continue;
                    }
                }
                if let Some(qname) = retired_qname {
                    overlap_remove(&mut self.waiting, &qname);
                    self.nodes[idx] = None;
                    self.free.push(idx);
                    self.order.remove(i);
                } else {
                    if let Some(qpos) = resolved {
                        entries.push(ColumnEntry { node: idx, qpos });
                    }
                    i += 1;
                }
            }

            // Advance tid/pos. With an empty buffer the C code reads its
            // zeroed sentinel node (tid 0, beg 0), which degenerates to
            // `pos += 1` for any engine state reachable here.
            match self
                .order
                .front()
                .copied()
                .and_then(|i| self.nodes[i].as_ref())
            {
                Some(head) => {
                    if self.tid < head.tid {
                        self.tid = head.tid;
                        self.pos = head.beg;
                    } else if self.pos < head.beg {
                        self.pos = head.beg;
                    } else {
                        self.pos += 1;
                    }
                }
                None => {
                    self.pos += 1;
                }
            }

            if !entries.is_empty() {
                return Some(Column {
                    pos: col_pos,
                    entries,
                });
            }
            if self.is_eof && self.order.is_empty() {
                break;
            }
        }
        None
    }

    /// `bam_plp_push` with a real record.
    fn push(&mut self, rec: R) {
        if rec.tid() < 0 {
            overlap_remove(&mut self.waiting, rec.qname());
            return;
        }
        // Skip only unmapped reads; any further filtering belongs to the feed
        // (stepper="nofilter" applies none).
        if rec.flags() & FLAG_UNMAP != 0 {
            overlap_remove(&mut self.waiting, rec.qname());
            return;
        }
        // maxcnt silent drop: the engine position pins at the read start
        // while same-start reads stream in; htslib compares its pool count
        // (sentinel included) against maxcnt, i.e. alive + 1 > MAXCNT.
        if self.tid == rec.tid() && self.pos == rec.pos() && self.order.len() + 1 > MAXCNT {
            self.dropped_maxcnt = true;
            overlap_remove(&mut self.waiting, rec.qname());
            return;
        }

        let beg = rec.pos();
        let end = beg + rec.ref_len();
        // Input comes from a region iterator, so the C unsorted-input checks
        // cannot trigger and are not replicated.
        self.max_tid = rec.tid();
        self.max_pos = beg;

        // Reads ending at/before the engine position never enter the buffer
        // (they still updated max_tid/max_pos above, as in C).
        if end > self.pos || rec.tid() > self.tid {
            let idx = match self.free.pop() {
                Some(i) => i,
                None => {
                    self.nodes.push(None);
                    self.nodes.len() - 1
                }
            };
            self.nodes[idx] = Some(Node {
                tid: rec.tid(),
                beg,
                end,
                rec,
                qual: None,
                cstate: CigarState { k: -1, x: 0, y: 0 },
});
            self.order.push_back(idx);
            self.overlap_push(idx);
        }
    }

    /// `overlap_push`: gate on pairing flags and geometry, then either pair
    /// with the waiting first end (applying the quality tweak to both) or
    /// enter the waiting table.
    fn overlap_push(&mut self, idx: usize) {
        let node = match self.nodes[idx].as_ref() {
            Some(n) => n,
            None => return,
        };
        let rec = &node.rec;
        let flags = rec.flags();
        // mapped mates in proper pairs only
        if flags & FLAG_MATE_UNMAP != 0 || flags & FLAG_PROPER_PAIR == 0 {
            return;
        }
        // no overlap possible, unless some wild cigar
        if (rec.mtid() >= 0 && rec.tid() != rec.mtid())
            || (rec.insert_size().abs() >= 2 * rec.seq_len() as i64 && rec.mpos() >= node.end)
        {
            return;
        }

        let qname = rec.qname().to_vec();
        if let Some(&a_idx) = self.waiting.get(&qname) {
            if a_idx != idx {
                // Pair hit: tweak both ends, then the pair leaves the table.
                let mut a_node = self.nodes[a_idx].take();
                let mut b_node = self.nodes[idx].take();
                if let (Some(a), Some(b)) = (&mut a_node, &mut b_node) {
                    tweak_overlap_quality(a, b);
                }
                if a_node.is_some() {
                    self.nodes[a_idx] = a_node;
                }
                if b_node.is_some() {
                    self.nodes[idx] = b_node;
                }
            }
            self.waiting.remove(&qname);
        } else if rec.mpos() >= rec.pos() || (flags & FLAG_PAIRED != 0 && rec.mpos() == -1) {
            // Only reads whose mate is still to arrive enter the table.
            self.waiting.insert(qname, idx);
        }
    }
}

/// `overlap_remove`: drop any waiting-table entry under this qname
/// (miss-tolerant, as in C).
fn overlap_remove(waiting: &mut HashMap<Vec<u8>, usize>, qname: &[u8]) {
    waiting.remove(qname);
}

/// CIGAR cursor for the overlap tweak (`cigar_iref2iseq_set`/`_next`).
/// `op` indexes the CIGAR; `icig` is the offset inside the current op;
/// `iseq`/`iref` are offsets relative to the read start.
#[derive(Clone, Copy)]
struct TweakCursor {
    op: usize,
    icig: i64,
    iseq: i64,
    iref: i64,
}

/// `cigar_iref2iseq_set`: find the first aligned (M/=/X) base at the given
/// reference offset. Returns None when the position is not covered by an
/// aligned op (htslib returns -1: "no overlap").
fn cig_set<R: PileRecord>(rec: &R, pos_offset: i64) -> Option<TweakCursor> {
    if pos_offset < 0 {
        return None;
    }
    let mut pos = pos_offset;
    let mut iseq: i64 = 0;
    let mut iref: i64 = 0;
    for k in 0..rec.cigar_len() {
        let c = rec.cigar_at(k);
        match c {
            Cigar::SoftClip(_) => {
                iseq += c.len() as i64;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
            Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => {
                let len = c.len() as i64;
                pos -= len;
                if pos < 0 {
                    let icig = len + pos;
                    iseq += icig;
                    iref += icig;
                    return Some(TweakCursor {
                        op: k,
                        icig,
                        iseq,
                        iref,
                    });
                }
                iseq += len;
                iref += len;
            }
            Cigar::Ins(_) => {
                iseq += c.len() as i64;
            }
            Cigar::Del(_) | Cigar::RefSkip(_) => {
                pos -= c.len() as i64;
                if pos < 0 {
                    pos = 0;
                }
                iref += c.len() as i64;
            }
        }
    }
    None
}

/// `cigar_iref2iseq_next`: advance the cursor to the next aligned base.
/// Returns false when the CIGAR is exhausted (htslib returns -1).
fn cig_next<R: PileRecord>(rec: &R, cur: &mut TweakCursor) -> bool {
    while cur.op < rec.cigar_len() {
        let op = rec.cigar_at(cur.op);
        match op {
            Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => {
                let len = op.len() as i64;
                if cur.icig >= len - 1 {
                    cur.icig = -1;
                    cur.op += 1;
                    continue;
                }
                cur.iseq += 1;
                cur.icig += 1;
                cur.iref += 1;
                return true;
            }
            Cigar::Del(_) | Cigar::RefSkip(_) => {
                cur.iref += op.len() as i64;
                cur.icig = -1;
                cur.op += 1;
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) => {
                cur.iseq += op.len() as i64;
                cur.icig = -1;
                cur.op += 1;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {
                cur.icig = -1;
                cur.op += 1;
            }
        }
    }
    cur.iseq = -1;
    cur.iref = -1;
    false
}

/// `tweak_overlap_quality`: adjust qualities of overlapping mates.
/// `a` is the first-arriving (left) end, `b` the second (right). Both nodes'
/// quality overlays are written — including partial writes on early exit,
/// matching the in-place C mutation. Bad CIGARs are tolerated by returning
/// early (never panicking).
fn tweak_overlap_quality<R: PileRecord>(a: &mut Node<R>, b: &mut Node<R>) {
    let a_rec = &a.rec;
    let b_rec = &b.rec;

    // Start at the right read's start; both cursors walk to that position.
    let iref_start = b_rec.pos();
    let mut a_cur = match cig_set(a_rec, iref_start - a_rec.pos()) {
        Some(c) => c,
        None => return, // no overlap
    };
    let mut b_cur = match cig_set(b_rec, 0) {
        Some(c) => c,
        None => return,
    };

    // Semi-random end selection: wang(x31(qname)) & 1.
    let pick = wang_hash(x31_hash(a_rec.qname())) & 1 != 0;
    let (amul, bmul) = if pick { (1u8, 0u8) } else { (0u8, 1u8) };

    // Lazy quality copies: only paired hits pay for the clone.
    let mut a_qual: Vec<u8> = a_rec.qual().to_vec();
    let mut b_qual: Vec<u8> = b_rec.qual().to_vec();

    let mut a_ok = true;
    let mut b_ok = true;
    let mut iref = iref_start;
    loop {
        // Del-chasing: advance a lagging cursor (D/N ops advance the
        // reference without consuming query) until it aligns to `iref`.
        // Reference order (pysam htslib): step BOTH cursors to the OLD
        // iref first, THEN raise iref to the max of the two — so a D-skip in
        // one end leaves the other end lagging (drives the chase below).
        // (Interleaving the max-raise between the steps would let the
        // lagging end catch up silently and skip the chase writes.)
        while a_ok && a_cur.iref >= 0 && a_cur.iref < iref - a_rec.pos() {
            a_ok = cig_next(a_rec, &mut a_cur);
        }
        if !a_ok {
            break;
        }
        while b_ok && b_cur.iref >= 0 && b_cur.iref < iref - b_rec.pos() {
            b_ok = cig_next(b_rec, &mut b_cur);
        }
        if !b_ok {
            break;
        }
        if iref < a_cur.iref + a_rec.pos() {
            iref = a_cur.iref + a_rec.pos();
        }
        if iref < b_cur.iref + b_rec.pos() {
            iref = b_cur.iref + b_rec.pos();
        }

        iref += 1;

        if a_rec.pos() + a_cur.iref != b_rec.pos() + b_cur.iref {
            // Deletion catch-up (pysam's bundled htslib keeps this chase;
            // upstream htslib replaced it with a plain continue — pysam is
            // our contract). When one end carries a D and has moved ahead,
            // the lagging end is stepped forward base by base, its quality
            // scaled by 0.8 (selected end) or zeroed at every chased base.
            // RefSkip (N) does NOT trigger the chase: it falls through to
            // `continue` like upstream.
            let a_abs = a_rec.pos() + a_cur.iref;
            let b_abs = b_rec.pos() + b_cur.iref;
            let prev_is_del = |c: &TweakCursor, rec: &R| -> bool {
                c.op > 0 && matches!(rec.cigar_at(c.op - 1), Cigar::Del(_))
            };
            if a_abs < b_abs && prev_is_del(&b_cur, b_rec) {
                // Del in B moved it ahead of A: chase A forward.
                loop {
                    let i = a_cur.iseq as usize;
                    if a_cur.iseq >= 0 && i < a_qual.len() {
                        a_qual[i] = if amul == 1 {
                            (a_qual[i] as f64 * 0.8) as u8
                        } else {
                            0
                        };
                    }
                    if !cig_next(a_rec, &mut a_cur) {
                        return;
                    }
                    if a_cur.iseq >= a_rec.seq_len() as i64 {
                        return; // bad-CIGAR tolerance (never panic)
                    }
                    if a_rec.pos() + a_cur.iref >= b_rec.pos() + b_cur.iref {
                        break;
                    }
                }
            } else if prev_is_del(&a_cur, a_rec) {
                // Del in A moved it ahead of B: chase B forward.
                loop {
                    let i = b_cur.iseq as usize;
                    if b_cur.iseq >= 0 && i < b_qual.len() {
                        b_qual[i] = if bmul == 1 {
                            (b_qual[i] as f64 * 0.8) as u8
                        } else {
                            0
                        };
                    }
                    if !cig_next(b_rec, &mut b_cur) {
                        return;
                    }
                    if b_cur.iseq >= b_rec.seq_len() as i64 {
                        return; // bad-CIGAR tolerance (never panic)
                    }
                    if b_rec.pos() + b_cur.iref >= a_rec.pos() + a_cur.iref {
                        break;
                    }
                }
            } else {
                // Ref-skip and other desyncs: unsupported, skip the position.
                continue;
            }
        }

        // Bad-CIGAR tolerance: htslib checks iseq > l_qseq; we also reject
        // iseq == l_qseq, which htslib reads (and writes!) out of bounds.
        // Either way: tolerate and stop, never panic.
        if a_cur.iseq < 0
            || b_cur.iseq < 0
            || a_cur.iseq >= a_rec.seq_len() as i64
            || b_cur.iseq >= b_rec.seq_len() as i64
        {
            break;
        }

        let ai = a_cur.iseq as usize;
        let bi = b_cur.iseq as usize;
        if a_rec.base_at(ai) == b_rec.base_at(bi) {
            // Confident match: sum of qualities capped at 200; keep the
            // selected end's sum, zero the other.
            let sum = a_qual[ai] as i32 + b_qual[bi] as i32;
            let cap: u8 = if sum > 200 { 200 } else { sum as u8 };
            a_qual[ai] = amul * cap;
            b_qual[bi] = bmul * cap;
        } else if a_qual[ai] > b_qual[bi] {
            // Mismatch: scale down the high-quality end, zero the other.
            a_qual[ai] = (0.8f64 * a_qual[ai] as f64) as u8;
            b_qual[bi] = 0;
        } else if a_qual[ai] < b_qual[bi] {
            b_qual[bi] = (0.8f64 * b_qual[bi] as f64) as u8;
            a_qual[ai] = 0;
        } else {
            // Equal mismatched qualities: keep the selected end at 0.8x.
            a_qual[ai] = (amul as f64 * 0.8f64 * a_qual[ai] as f64) as u8;
            b_qual[bi] = (bmul as f64 * 0.8f64 * b_qual[bi] as f64) as u8;
        }
    }

    a.qual = Some(a_qual);
    b.qual = Some(b_qual);
}

/// khash `__ac_X31_hash_string` over the qname bytes (u32 wrapping).
fn x31_hash(s: &[u8]) -> u32 {
    let first = match s.first() {
        Some(&c) => c as u32,
        None => 0,
    };
    let mut h: u32 = first;
    if first != 0 {
        for &c in &s[1..] {
            h = (h << 5).wrapping_sub(h).wrapping_add(c as u32);
        }
    }
    h
}

/// khash `__ac_Wang_hash` (u32 wrapping).
fn wang_hash(mut key: u32) -> u32 {
    key = key.wrapping_add(!(key << 15));
    key ^= key >> 10;
    key = key.wrapping_add(key << 3);
    key ^= key >> 6;
    key = key.wrapping_add(!(key << 11));
    key ^= key >> 16;
    key
}
