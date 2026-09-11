//! candidates.bed output (10 columns):
//! `chrom  pos0  pos0+1  ref  depth  A  C  G  T  edit_freq`
//!
//! `ref` is the reference base (A or T); the editing-consistent alternate is
//! G for A and C for T. Base counts and depth cover the reads that passed the
//! mapping- and base-quality gates. `edit_freq` = alternate reads / depth.
//! Sorted by (chrom, pos) before writing.

use crate::error::CallError;
use crate::Channel;
use std::io::Write as _;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Candidate {
    pub chrom: String,
    pub pos0: u32,
    pub ref_base: u8,
    pub depth: u32,
    /// Counts of A, C, G, T.
    pub counts: [u32; 4],
    pub edit_freq: f64,
    /// Standard evidence-filter tags (comma-separated; empty = clean).
    pub filters: String,
    /// Mismatch channel this candidate belongs to.
    pub channel: Channel,
}

pub fn write_bed(path: &Path, cands: &mut [Candidate]) -> Result<(), CallError> {
    cands.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.pos0.cmp(&b.pos0)));
    let io = |e: std::io::Error| CallError::Io {
        path: path.display().to_string(),
        source: e,
    };
    let mut w = std::io::BufWriter::new(std::fs::File::create(path).map_err(io)?);
    for c in cands.iter() {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}",
            c.chrom,
            c.pos0,
            c.pos0 + 1,
            c.ref_base as char,
            c.depth,
            c.counts[0],
            c.counts[1],
            c.counts[2],
            c.counts[3],
            c.edit_freq,
            c.filters
        )
        .map_err(io)?;
    }
    w.flush().map_err(io)
}

/// control_candidates.bed output: the candidates.bed columns plus the
/// channel label (`E` never appears here; rows are `C1` or `C2`).
pub fn write_control_bed(path: &Path, cands: &mut [Candidate]) -> Result<(), CallError> {
    cands.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.pos0.cmp(&b.pos0)));
    let io = |e: std::io::Error| CallError::Io {
        path: path.display().to_string(),
        source: e,
    };
    let mut w = std::io::BufWriter::new(std::fs::File::create(path).map_err(io)?);
    for c in cands.iter() {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}",
            c.chrom,
            c.pos0,
            c.pos0 + 1,
            c.ref_base as char,
            c.depth,
            c.counts[0],
            c.counts[1],
            c.counts[2],
            c.counts[3],
            c.edit_freq,
            c.filters,
            c.channel.label()
        )
        .map_err(io)?;
    }
    w.flush().map_err(io)
}
