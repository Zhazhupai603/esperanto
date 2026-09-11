//! Junction library built from annotation.
//!
//! `Junction` is a 0-based half-open intron interval on a contig with strand
//! orientation. `JunctionLib` stores junctions sorted with parallel support
//! counts plus an end-point index for acceptor-side range queries.
//!
//! `SpliceSignal` and `RefinedJunction` live here because `RefinedJunction`
//! references the signal enum; the splice module (later wave) adds its
//! scoring machinery on top.

use crate::error::AlignError;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

/// A 0-based half-open intron `[start, end)` on `contig`.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct Junction {
    /// Contig index.
    pub contig: u32,
    /// Intron start (donor side, plus-strand coordinates).
    pub start: u32,
    /// Intron end (acceptor side, exclusive).
    pub end: u32,
    /// Transcribed on the minus strand.
    pub minus_strand: bool,
}

/// Splice signal classification at donor/acceptor dinucleotides.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SpliceSignal {
    /// GT..AG.
    GtAg,
    /// GC..AG.
    GcAg,
    /// AT..AC.
    AtAc,
    /// Anything else.
    NonCanonical,
}

impl SpliceSignal {
    /// Short label used in tags and evidence strings.
    pub fn label(self) -> &'static str {
        match self {
            SpliceSignal::GtAg => "GT-AG",
            SpliceSignal::GcAg => "GC-AG",
            SpliceSignal::AtAc => "AT-AC",
            SpliceSignal::NonCanonical => "-",
        }
    }
}

/// A junction with refined breakpoints and provenance.
#[derive(Clone, Debug)]
pub struct RefinedJunction {
    /// The underlying intron interval.
    pub junction: Junction,
    /// Signal classification at the refined breakpoints.
    pub signal: SpliceSignal,
    /// Support count from the library (0 = de-novo).
    pub known_support: u32,
}

/// Sorted junction table with support counts and an end-point index.
#[derive(Clone, Debug, Default)]
pub struct JunctionLib {
    /// Junctions sorted by (contig, start, end, minus_strand).
    pub junctions: Vec<Junction>,
    /// Support counts, parallel to `junctions`.
    pub counts: Vec<u32>,
    /// Indices into `junctions` sorted by (contig, end, start).
    pub by_end: Vec<u32>,
}

impl JunctionLib {
    /// Build from an unsorted stream of junctions: sort, merge identical
    /// junctions summing their counts.
    pub fn build<I: IntoIterator<Item = Junction>>(junctions: I) -> JunctionLib {
        let mut items: Vec<(Junction, u32)> = Vec::new();
        for j in junctions {
            items.push((j, 1));
        }
        items.sort_by_key(|(j, _)| (j.contig, j.start, j.end, j.minus_strand));
        let mut merged: Vec<(Junction, u32)> = Vec::new();
        for (j, c) in items {
            match merged.last_mut() {
                Some((lj, lc)) if *lj == j => *lc += c,
                _ => merged.push((j, c)),
            }
        }
        let junctions: Vec<Junction> = merged.iter().map(|(j, _)| *j).collect();
        let counts: Vec<u32> = merged.iter().map(|(_, c)| *c).collect();
        let mut by_end: Vec<u32> = (0..junctions.len() as u32).collect();
        by_end.sort_by_key(|&i| {
            let j = &junctions[i as usize];
            (j.contig, j.end, j.start)
        });
        JunctionLib {
            junctions,
            counts,
            by_end,
        }
    }

    /// Build from (junction, count) pairs, preserving counts (2-pass merge path).
    pub fn build_with_counts(items: Vec<(Junction, u32)>) -> JunctionLib {
        let mut items = items;
        items.sort_by_key(|(j, _)| (j.contig, j.start, j.end, j.minus_strand));
        let mut merged: Vec<(Junction, u32)> = Vec::new();
        for (j, c) in items {
            match merged.last_mut() {
                Some((lj, lc)) if *lj == j => *lc += c,
                _ => merged.push((j, c)),
            }
        }
        let junctions: Vec<Junction> = merged.iter().map(|(j, _)| *j).collect();
        let counts: Vec<u32> = merged.iter().map(|(_, c)| *c).collect();
        let mut by_end: Vec<u32> = (0..junctions.len() as u32).collect();
        by_end.sort_by_key(|&i| {
            let j = &junctions[i as usize];
            (j.contig, j.end, j.start)
        });
        JunctionLib {
            junctions,
            counts,
            by_end,
        }
    }

    /// Whether the library holds no junctions.
    pub fn is_empty(&self) -> bool {
        self.junctions.is_empty()
    }

    /// Exact junction support (0 if absent).
    pub fn support(&self, j: &Junction) -> u32 {
        self.junctions
            .binary_search_by_key(&(j.contig, j.start, j.end, j.minus_strand), |x| {
                (x.contig, x.start, x.end, x.minus_strand)
            })
            .map(|i| self.counts[i])
            .unwrap_or(0)
    }

    /// Exact membership check.
    pub fn contains(&self, j: &Junction) -> bool {
        self.junctions
            .binary_search_by_key(&(j.contig, j.start, j.end, j.minus_strand), |x| {
                (x.contig, x.start, x.end, x.minus_strand)
            })
            .is_ok()
    }

    /// Junctions on `contig` with `start` in `[lo, hi)`, as a contiguous
    /// subslice of `(junctions, counts)`.
    pub fn range_start(&self, contig: u32, lo: u32, hi: u32) -> (&[Junction], &[u32]) {
        let (a, b) = self.range_contig_bounds(contig);
        let run = &self.junctions[a..b];
        let lo_i = run.partition_point(|j| j.start < lo);
        let hi_i = run.partition_point(|j| j.start < hi);
        (
            &self.junctions[a + lo_i..a + hi_i],
            &self.counts[a + lo_i..a + hi_i],
        )
    }

    /// Junctions on `contig` with `end` in `[lo, hi)`, as a contiguous slice
    /// of the `by_end` index (indices into `junctions`).
    pub fn range_end(&self, contig: u32, lo: u32, hi: u32) -> &[u32] {
        let run = self.by_end_contig_run(contig);
        let lo_i = run.partition_point(|&i| self.junctions[i as usize].end < lo);
        let hi_i = run.partition_point(|&i| self.junctions[i as usize].end < hi);
        &run[lo_i..hi_i]
    }

    fn range_contig_bounds(&self, contig: u32) -> (usize, usize) {
        let lo = self.junctions.partition_point(|j| j.contig < contig);
        let hi = self.junctions.partition_point(|j| j.contig <= contig);
        (lo, hi)
    }

    fn by_end_contig_run(&self, contig: u32) -> &[u32] {
        let lo = self
            .by_end
            .partition_point(|&i| self.junctions[i as usize].contig < contig);
        let hi = self
            .by_end
            .partition_point(|&i| self.junctions[i as usize].contig <= contig);
        &self.by_end[lo..hi]
    }
}


/// Extract `key "value"` from a GTF attributes column.
pub fn extract_attr(attrs: &str, key: &str) -> Option<String> {
    for part in attrs.split(';') {
        let part = part.trim();
        let mut fields = part.split_whitespace();
        if fields.next() == Some(key) {
            if let Some(v) = fields.next() {
                return Some(v.trim_matches('"').to_string());
            }
        }
    }
    None
}

/// Build a junction library from a GTF file.
///
/// Exon lines are grouped by `transcript_id`; each pair of adjacent exons
/// (sorted by start) contributes one 0-based half-open intron between
/// `exon_i.end` and `exon_(i+1).start − 1`. The feature strand determines
/// `minus_strand`. Non-exon lines are skipped; malformed exon lines error
/// with the 1-based line number.
pub fn from_gtf<F: Fn(&str) -> u32>(
    path: &Path,
    contig_id_fn: F,
) -> Result<JunctionLib, AlignError> {
    let bytes = fs::read(path).map_err(|source| AlignError::FastaIo {
        path: path.display().to_string(),
        source,
    })?;
    let text = String::from_utf8_lossy(&bytes);
    let mut junctions: Vec<Junction> = Vec::new();
    // transcript_id -> exons (contig, minus_strand, 0-based half-open list)
    let mut tx: HashMap<String, TranscriptExons> = HashMap::new();

    for (idx, line) in text.lines().enumerate() {
        let lineno = idx + 1;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }
        if cols[2] != "exon" {
            continue;
        }
        let contig = contig_id_fn(cols[0]);
        let start: u32 = match cols[3].parse() {
            Ok(v) => v,
            Err(_) => {
                return Err(AlignError::FastaFormat {
                    line: lineno,
                    msg: format!("bad exon start '{}'", cols[3]),
                })
            }
        };
        let end: u32 = match cols[4].parse() {
            Ok(v) => v,
            Err(_) => {
                return Err(AlignError::FastaFormat {
                    line: lineno,
                    msg: format!("bad exon end '{}'", cols[4]),
                })
            }
        };
        let minus_strand = cols[6] == "-";
        let transcript = extract_attr(cols[8], "transcript_id")
            .unwrap_or_else(|| format!("__line_{}", lineno));
        let entry = tx.entry(transcript).or_insert_with(|| TranscriptExons {
            contig,
            minus_strand,
            exons: Vec::new(),
        });
        entry.exons.push((start.saturating_sub(1), end)); // 0-based half-open exon
    }

    for (_tid, t) in tx {
        if t.exons.len() < 2 {
            continue;
        }
        let mut exons = t.exons;
        exons.sort_unstable();
        for w in exons.windows(2) {
            let (.., e1_end) = w[0];
            let (e2_start, ..) = w[1];
            if e2_start > e1_end {
                junctions.push(Junction {
                    contig: t.contig,
                    start: e1_end,
                    end: e2_start,
                    minus_strand: t.minus_strand,
                });
            }
        }
    }

    Ok(JunctionLib::build(junctions))
}

struct TranscriptExons {
    contig: u32,
    minus_strand: bool,
    exons: Vec<(u32, u32)>,
}
