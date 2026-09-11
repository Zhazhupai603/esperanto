//! Repeat-region tracking: the [`RepeatTrack`] oracle, the empty track
//! and the BED loader.

use std::collections::BTreeMap;
use std::io::Read;
use std::path::Path;

use crate::{Error, RepeatTrack};

/// Empty repeat track: nothing is a repeat (default / tests / rejection
/// disabled).
#[derive(Debug, Clone, Copy, Default)]
pub struct NoRepeats;

impl RepeatTrack for NoRepeats {
    fn overlaps(&self, _contig: u32, _pos: u32, _len: u32) -> bool {
        false
    }
}

/// Repeat regions loaded from a BED file (first three columns; plain text
/// or `.gz`). Rows whose contig is absent from the id mapping are
/// silently skipped.
///
/// Intervals are stored per contig sorted by start with a parallel start
/// array; `overlaps` binary-searches the rightmost interval with
/// `start < pos + len` and checks the 3 preceding entries for
/// `end > pos` (RepeatMasker intervals may overlap).
#[derive(Debug, Clone, Default)]
pub struct RepeatBed {
    regions: BTreeMap<u32, (Vec<u32>, Vec<u32>)>,
}

impl RepeatBed {
    /// Load a BED file (plain or gzip by `.gz` extension), mapping contig
    /// names to ids via `contig_name_to_id`.
    pub fn load(path: &Path, contig_name_to_id: &BTreeMap<String, u32>) -> Result<Self, Error> {
        let file = std::fs::File::open(path)?;
        let mut text = Vec::new();
        if path
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("gz"))
        {
            let mut reader = flate2::read::GzDecoder::new(file);
            reader.read_to_end(&mut text)?;
        } else {
            let mut reader = file;
            reader.read_to_end(&mut text)?;
        }
        RepeatBed::parse(&text, contig_name_to_id)
    }

    /// Parse BED bytes (same semantics as [`RepeatBed::load`]).
    pub fn parse(text: &[u8], contig_name_to_id: &BTreeMap<String, u32>) -> Result<Self, Error> {
        let mut raw: BTreeMap<u32, Vec<(u32, u32)>> = BTreeMap::new();
        for line in text.split(|&b| b == b'\n') {
            let line = trim_ascii(line);
            if line.is_empty() || line.first() == Some(&b'#') {
                continue;
            }
            let mut cols = line.split(|&b| b == b'\t');
            let (Some(contig), Some(start), Some(end)) = (cols.next(), cols.next(), cols.next())
            else {
                continue;
            };
            let contig = std::str::from_utf8(contig)
                .map_err(|_| Error::Format("BED contig is not valid UTF-8".to_string()))?;
            let Some(&contig_id) = contig_name_to_id.get(contig) else {
                continue; // contig not in the mapping: silently skipped
            };
            let (Ok(start), Ok(end)) = (
                std::str::from_utf8(start)
                    .map_err(|_| Error::Format("BED start is not valid UTF-8".to_string()))?
                    .parse::<u32>(),
                std::str::from_utf8(end)
                    .map_err(|_| Error::Format("BED end is not valid UTF-8".to_string()))?
                    .parse::<u32>(),
            ) else {
                continue;
            };
            if end > start {
                raw.entry(contig_id).or_default().push((start, end));
            }
        }
        let mut regions = BTreeMap::new();
        for (contig_id, mut ivs) in raw {
            ivs.sort_unstable();
            let mut starts = Vec::with_capacity(ivs.len());
            let mut ends = Vec::with_capacity(ivs.len());
            for (s, e) in ivs {
                starts.push(s);
                ends.push(e);
            }
            regions.insert(contig_id, (starts, ends));
        }
        Ok(RepeatBed { regions })
    }

    /// True when nothing was loaded.
    pub fn is_empty(&self) -> bool {
        self.regions.is_empty()
    }
}

impl RepeatTrack for RepeatBed {
    fn overlaps(&self, contig: u32, pos: u32, len: u32) -> bool {
        let Some((starts, ends)) = self.regions.get(&contig) else {
            return false;
        };
        let q_end = pos.saturating_add(len);
        let idx = starts.partition_point(|&s| s < q_end);
        if idx == 0 {
            return false;
        }
        // check the 3 intervals preceding the rightmost start < pos + len
        let lo = idx.saturating_sub(3);
        (lo..idx).any(|i| ends[i] > pos)
    }
}

/// Trim ASCII whitespace from both ends.
fn trim_ascii(line: &[u8]) -> &[u8] {
    let mut s = 0;
    let mut e = line.len();
    while s < e && line[s].is_ascii_whitespace() {
        s += 1;
    }
    while e > s && line[e - 1].is_ascii_whitespace() {
        e -= 1;
    }
    &line[s..e]
}
