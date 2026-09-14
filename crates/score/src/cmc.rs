//! CMC feature extraction: 16x50 flank-center joint base frequencies from BAM
//! (schema `cmc-cigar-window-1`).
//!
//! Frozen contract, ported bit-for-bit from the training extractor:
//!
//! * window = reference interval `[center-25, center+26)` (0-based, half-open),
//!   `center = pos-1`; the center column itself is excluded from the tensor.
//! * reads: `flag & 3844 == 0` and `MAPQ >= 20`; a read's per-position base
//!   contributes when the base quality is `>= 13` and the base is one of
//!   `A C T G` (uppercase). A missing quality byte (0xFF) counts as quality 0.
//! * CIGAR walk: `M`/`=`/`X` consume both tracks, `I`/`S` consume the query,
//!   `D`/`N` consume the reference, `H`/`P` consume neither; the walk stops
//!   once the reference cursor reaches the window end.
//! * template deduplication: read names are stripped of a trailing `/1` or
//!   `/2`; within a template the highest-quality base at each reference
//!   position wins and an equal-quality disagreement stores `N`.
//! * tensor row = `4 * flank_index + center_index` over `ACTG` order; column
//!   `c` holds the joint frequency for offset `-25..-1, +1..+25`; a covered
//!   column sums to one, an uncovered column is all zeros (the normalization
//!   divides integer counts in f32, matching the reference implementation
//!   bit-for-bit).
//! * a read tagged `RE:Z:collapsed` that yields the center base marks its
//!   template as collapsed (quality control; not part of the tensor).
//!
//! Determinism: counts are integers accumulated in template order, and the
//! per-column normalization is a single f32 division — the output is
//! thread-count independent and identical across runs.

use std::collections::HashMap;

use rust_htslib::bam::record::{Aux, Cigar, Record};
use rust_htslib::bam::{IndexedReader, Read as _};
use thiserror::Error;

/// Tensor rows: 4 (flank) x 4 (center) joint channels.
pub const CMC_CHANNELS: usize = 16;
/// Tensor columns: offsets -25..-1 and +1..+25.
pub const CMC_COLUMNS: usize = 50;
/// Half-width of the flank window in reference bases.
pub const CMC_HALF_WINDOW: i64 = 25;
/// Records with any of these flags set never contribute (unmapped, secondary,
/// failed QC, duplicate, supplementary).
pub const CMC_EXCLUDE_FLAGS: u16 = 3844;
/// Minimum mapping quality (same gate as the candidate scan).
pub const CMC_MIN_MAPQ: u8 = 20;
/// Minimum base quality for a base to contribute.
pub const CMC_MIN_BASE_QUALITY: u8 = 13;
/// Batch mode: consecutive sites on one reference within this gap are served
/// from a single region fetch (matches the pileup feature pass).
pub const CMC_MERGE_GAP: i64 = 2000;

/// Flat tensor, row-major `[channel * CMC_COLUMNS + column]`.
pub type CmcFeatures = [f32; CMC_CHANNELS * CMC_COLUMNS];

/// Errors produced by CMC extraction.
#[derive(Debug, Error)]
pub enum CmcError {
    /// Reference name not present in the BAM header.
    #[error("contig not found in BAM header: {0}")]
    ContigNotFound(String),
    /// Site position is not a valid 1-based coordinate.
    #[error("invalid 1-based position: {0}")]
    InvalidPosition(i64),
    /// Region fetch failed.
    #[error("failed to fetch region on tid {tid} at {start}-{end}: {source}")]
    Fetch {
        tid: i32,
        start: i64,
        end: i64,
        source: rust_htslib::errors::Error,
    },
    /// Record iteration failed.
    #[error("failed to read BAM record: {0}")]
    Read(#[from] rust_htslib::errors::Error),
}

/// True when the record carries the collapsed-rescue provenance tag.
fn is_collapsed(rec: &Record) -> bool {
    matches!(rec.aux(b"RE"), Ok(Aux::String(v)) if v == "collapsed")
}

/// Template name: the read name with a trailing `/1` or `/2` removed.
fn template_name(qname: &[u8]) -> &[u8] {
    match qname {
        [rest @ .., b'/', b'1' | b'2'] => rest,
        _ => qname,
    }
}

/// Channel index of a base byte, or None when the base is not one of ACTG.
fn base_index(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'T' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// Per-site accumulator: template -> best base per reference position.
struct SiteAcc {
    templates: Vec<HashMap<i64, (u8, u8)>>,
    collapsed_center: Vec<bool>,
    index: HashMap<Vec<u8>, usize>,
}

impl SiteAcc {
    fn new() -> Self {
        Self {
            templates: Vec::new(),
            collapsed_center: Vec::new(),
            index: HashMap::new(),
        }
    }

    fn intern(&mut self, name: &[u8]) -> usize {
        if let Some(&i) = self.index.get(name) {
            return i;
        }
        let i = self.templates.len();
        self.templates.push(HashMap::new());
        self.collapsed_center.push(false);
        self.index.insert(name.to_vec(), i);
        i
    }

    /// Feed one record's CIGAR walk over `[start, end)` (0-based, half-open).
    /// `center` marks the candidate position (0-based).
    fn feed(&mut self, rec: &Record, start: i64, end: i64, center: i64) {
        let name = template_name(rec.qname()).to_vec();
        let collapsed = is_collapsed(rec);
        let ti = self.intern(&name);
        let seq = rec.seq();
        let qual = rec.qual();
        let mut q: i64 = 0;
        let mut r: i64 = rec.pos();
        for op in rec.cigar().iter() {
            if r >= end {
                break;
            }
            let len = i64::from(op.len());
            match op {
                Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => {
                    let lo = r.max(start);
                    let hi = (r + len).min(end);
                    for rp in lo..hi {
                        let qp = (q + rp - r) as usize;
                        let qbyte = if qp < qual.len() { qual[qp] } else { 0 };
                        // 0xFF marks a missing quality; treat it as quality 0.
                        let qb = if qbyte == 0xFF { 0 } else { qbyte };
                        if qb < CMC_MIN_BASE_QUALITY {
                            continue;
                        }
                        if qp >= seq.len() {
                            continue;
                        }
                        let base = seq[qp];
                        if base_index(base).is_none() {
                            continue;
                        }
                        let tpl = &mut self.templates[ti];
                        match tpl.get(&rp) {
                            None => {
                                tpl.insert(rp, (base, qb));
                            }
                            Some(&(pb, pq)) => {
                                if qb > pq {
                                    tpl.insert(rp, (base, qb));
                                } else if qb == pq && base != pb {
                                    tpl.insert(rp, (b'N', qb));
                                }
                            }
                        }
                        if rp == center && collapsed {
                            self.collapsed_center[ti] = true;
                        }
                    }
                    q += len;
                    r += len;
                }
                Cigar::Ins(_) | Cigar::SoftClip(_) => q += len,
                Cigar::Del(_) | Cigar::RefSkip(_) => r += len,
                Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }
    }

    /// Collapse the accumulator into the normalized 16x50 tensor.
    fn finish(&self, center: i64) -> CmcFeatures {
        let mut counts = [[0u32; CMC_COLUMNS]; CMC_CHANNELS];
        let offsets = column_offsets();
        for tpl in &self.templates {
            let Some(&(base, _)) = tpl.get(&center) else {
                continue;
            };
            let Some(bi) = base_index(base) else {
                continue;
            };
            for (col, &offset) in offsets.iter().enumerate() {
                let Some(&(flank, _)) = tpl.get(&(center + offset)) else {
                    continue;
                };
                if let Some(fi) = base_index(flank) {
                    counts[fi * 4 + bi][col] += 1;
                }
            }
        }
        // Per-column pair depth (sum over the 16 channels), then divide.
        let mut col_depth = [0u32; CMC_COLUMNS];
        for channel in &counts {
            for (col, &c) in channel.iter().enumerate() {
                col_depth[col] += c;
            }
        }
        let mut out = [0.0f32; CMC_CHANNELS * CMC_COLUMNS];
        for (ch, channel) in counts.iter().enumerate() {
            for (col, &c) in channel.iter().enumerate() {
                if col_depth[col] > 0 {
                    out[ch * CMC_COLUMNS + col] = c as f32 / col_depth[col] as f32;
                }
            }
        }
        out
    }
}

/// Flank offsets for the 50 tensor columns: -25..-1 then +1..+25.
fn column_offsets() -> [i64; CMC_COLUMNS] {
    let mut a = [0i64; CMC_COLUMNS];
    for (i, slot) in a.iter_mut().enumerate() {
        *slot = if i < 25 {
            i as i64 - 25
        } else {
            (i - 24) as i64
        };
    }
    a
}

fn resolve_site(bam: &mut IndexedReader, chrom: &str, pos_1based: i64) -> Result<(i32, i64), CmcError> {
    if pos_1based < 1 {
        return Err(CmcError::InvalidPosition(pos_1based));
    }
    let tid = bam
        .header()
        .tid(chrom.as_bytes())
        .map(|t| t as i32)
        .ok_or_else(|| CmcError::ContigNotFound(chrom.to_string()))?;
    Ok((tid, pos_1based - 1))
}

/// Extract CMC tensors for many sites in one pass per site group.
///
/// Sites are sorted by `(reference, position)`; runs on one reference within
/// `CMC_MERGE_GAP` share a single region fetch and a streaming record pass
/// (each record is fed to every site window it overlaps — the same record set
/// a per-site fetch would deliver). Output order follows the input `sites`.
pub fn extract_cmc_batch(
    bam: &mut IndexedReader,
    sites: &[(&str, i64)],
) -> Result<Vec<CmcFeatures>, CmcError> {
    if sites.is_empty() {
        return Ok(Vec::new());
    }
    let mut order: Vec<(i32, i64, usize)> = Vec::with_capacity(sites.len());
    for (i, (chrom, pos)) in sites.iter().enumerate() {
        let (tid, pos0) = resolve_site(bam, chrom, *pos)?;
        order.push((tid, pos0, i));
    }
    order.sort_by_key(|&(tid, pos0, _)| (tid, pos0));

    let mut out: Vec<CmcFeatures> = vec![[0.0; CMC_CHANNELS * CMC_COLUMNS]; sites.len()];
    let mut group_start = 0usize;
    while group_start < order.len() {
        let (g_tid, g_first, _) = order[group_start];
        let mut group_end = group_start + 1;
        while group_end < order.len() {
            let &(tid, pos0, _) = &order[group_end];
            let &(_, prev, _) = &order[group_end - 1];
            if tid != g_tid || pos0 - prev > CMC_MERGE_GAP {
                break;
            }
            group_end += 1;
        }
        let group = &order[group_start..group_end];
        let last = group[group.len() - 1].1;
        let fetch_start = (g_first - CMC_HALF_WINDOW).max(0);
        let fetch_end = last + CMC_HALF_WINDOW + 1;
        bam.fetch((g_tid, fetch_start, fetch_end))
            .map_err(|source| CmcError::Fetch {
                tid: g_tid,
                start: fetch_start,
                end: fetch_end,
                source,
            })?;

        // Streaming pass: active sites are those whose window intersects the
        // current record. Sites leave when the record stream passed their
        // window end; sites enter while their window starts before the
        // record's reference end.
        let mut accs: Vec<SiteAcc> = group.iter().map(|_| SiteAcc::new()).collect();
        let mut next_site = 0usize; // index into `group` (sorted by pos0)
        let mut active: Vec<usize> = Vec::new(); // group-relative site indexes
        for rec in bam.records() {
            let rec = rec?;
            if rec.pos() > last + CMC_HALF_WINDOW {
                break;
            }
            let rec_start = rec.pos();
            let rec_end = rec_end(&rec);
            while next_site < group.len() && group[next_site].1 - CMC_HALF_WINDOW < rec_end {
                active.push(next_site);
                next_site += 1;
            }
            if rec.flags() & CMC_EXCLUDE_FLAGS != 0 || rec.mapq() < CMC_MIN_MAPQ {
                // Window bookkeeping still advances; the record never feeds.
                active.retain(|&si| group[si].1 + CMC_HALF_WINDOW + 1 > rec_start);
                continue;
            }
            for &si in &active {
                let (pos0, _) = (group[si].1, group[si].2);
                let ws = (pos0 - CMC_HALF_WINDOW).max(0);
                let we = pos0 + CMC_HALF_WINDOW + 1;
                if rec_start < we && rec_end > ws {
                    accs[si].feed(&rec, ws, we, pos0);
                }
            }
            active.retain(|&si| group[si].1 + CMC_HALF_WINDOW + 1 > rec_start);
        }
        for (si, &(.., orig)) in group.iter().enumerate() {
            out[orig] = accs[si].finish(group[si].1);
        }
        group_start = group_end;
    }
    Ok(out)
}

/// Reference end (0-based, exclusive) of a record: start + reference-consuming
/// CIGAR length.
fn rec_end(rec: &Record) -> i64 {
    let mut rlen: i64 = 0;
    for c in rec.cigar().iter() {
        match c {
            Cigar::Match(_)
            | Cigar::Del(_)
            | Cigar::RefSkip(_)
            | Cigar::Equal(_)
            | Cigar::Diff(_) => rlen += i64::from(c.len()),
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    rec.pos() + rlen
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn template_names_strip_paired_suffix_only() {
        assert_eq!(template_name(b"read/1"), b"read");
        assert_eq!(template_name(b"read/2"), b"read");
        assert_eq!(template_name(b"read/3"), b"read/3");
        assert_eq!(template_name(b"read"), b"read");
        assert_eq!(template_name(b"r/1/1"), b"r/1");
    }

    #[test]
    fn channel_order_is_actg() {
        assert_eq!(base_index(b'A'), Some(0));
        assert_eq!(base_index(b'C'), Some(1));
        assert_eq!(base_index(b'T'), Some(2));
        assert_eq!(base_index(b'G'), Some(3));
        assert_eq!(base_index(b'N'), None);
        assert_eq!(base_index(b'a'), None);
    }

    #[test]
    fn offsets_skip_the_center_column() {
        let o = column_offsets();
        assert_eq!(o.len(), CMC_COLUMNS);
        assert_eq!(o[0], -25);
        assert_eq!(o[24], -1);
        assert_eq!(o[25], 1);
        assert_eq!(o[49], 25);
    }

    #[test]
    fn finish_normalizes_each_covered_column() {
        let mut acc = SiteAcc::new();
        let center = 1000i64;
        // two templates: (flank A, center G) and (flank A, center G) at
        // offset +1 -> channel 0*4+3; one template (flank C, center A) at
        // offset +2 -> channel 1*4+0. Offset +3 stays uncovered.
        for (name, entries) in [
            (b"t0", vec![(center + 1, b'A'), (center, b'G')]),
            (b"t1", vec![(center + 1, b'A'), (center, b'G')]),
            (b"t2", vec![(center + 2, b'C'), (center, b'A')]),
        ] {
            let ti = acc.intern(name);
            for (pos, base) in entries {
                acc.templates[ti].insert(pos, (base, 40));
            }
        }
        let out = acc.finish(center);
        // column for +1 (index 25): channel 3 (A flank x G center) = 1.0
        assert_eq!(out[3 * CMC_COLUMNS + 25], 1.0);
        // column for +2 (index 26): channels 0 (C,A) and 3 (C,G)... t2 is
        // (C flank, A center) = channel 4; t0/t1 have no base at +2.
        assert_eq!(out[4 * CMC_COLUMNS + 26], 1.0);
        // column +3 uncovered -> zeros
        for ch in 0..CMC_CHANNELS {
            assert_eq!(out[ch * CMC_COLUMNS + 27], 0.0);
        }
    }

    #[test]
    fn finish_ignores_templates_without_a_valid_center() {
        let mut acc = SiteAcc::new();
        let center = 500i64;
        acc.intern(b"t0");
        acc.templates[0].insert(center + 5, (b'A', 40)); // center itself missing
        let out = acc.finish(center);
        assert!(out.iter().all(|&v| v == 0.0));
    }
}
