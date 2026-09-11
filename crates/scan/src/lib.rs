//! esperanto-scan — A-to-I candidate site discovery.
//!
//! Genome-wide pileup over a coordinate-sorted BAM: every read that passes the
//! mapping-quality gate contributes its bases (above the base-quality gate) to
//! per-position depth and base counts. A position is a candidate when the
//! reference base is A or T and the editing-consistent alternate base (G for
//! A, C for T) is supported by enough reads. Known SNP positions are removed
//! when a SNP VCF is given.
//!
//! The counting convention matches the standard REDItools workflow (MAPQ >= 20,
//! base quality >= 30, coverage >= 10, no read-orientation split), so
//! candidate tables are directly comparable across pipelines.
//!
//! One filter is applied on top: a candidate whose read-length neighbourhood
//! holds both editing directions is not one transcript's editing and is
//! dropped (see [`strand`]).
//!
//! Parallel over (contig, 32 Mbp block); output sorted by (chrom, pos) and
//! independent of the thread count.

pub mod call;
pub mod filters;
pub mod error;
pub mod out;
pub mod scatter;
pub mod snp;

pub use error::CallError;

use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct CallParams {
    /// Input BAM (coordinate-sorted, with `.bai`/`.csi`).
    pub bam: PathBuf,
    /// Output candidates.bed.
    pub out: PathBuf,
    /// Reference FASTA (with `.fai`).
    pub fasta: PathBuf,
    /// Known-SNP VCF (plain or gzip); listed positions are excluded.
    pub snp: Option<PathBuf>,
    /// Reads below this MAPQ do not count.
    pub min_mapq: u8,
    /// Bases below this Phred quality do not count.
    pub min_bq: u8,
    /// Minimum counted depth at a candidate.
    pub min_depth: u32,
    /// Minimum reads carrying the editing-consistent alternate base.
    pub min_edit_reads: u32,
    /// Minimum editing-consistent alternate frequency.
    pub min_edit_freq: f64,
    /// Thread count; 0 = all cores.
    pub threads: usize,
}

impl CallParams {
    /// Standard thresholds (REDItools convention): MAPQ 20, BQ 30, depth 10,
    /// 2 alternate reads, 1% frequency.
    pub fn standard(bam: PathBuf, out: PathBuf, fasta: PathBuf) -> Self {
        Self {
            bam,
            out,
            fasta,
            snp: None,
            min_mapq: 20,
            min_bq: 30,
            min_depth: 10,
            min_edit_reads: 2,
            min_edit_freq: 0.01,
            threads: 0,
        }
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct CallStats {
    /// Candidates written.
    pub candidates: usize,
    /// Positions with editing-consistent evidence that failed a threshold.
    pub below_threshold: usize,
    /// Positions removed as known SNPs.
    pub snp_removed: usize,
    /// Per-channel site / alt-read totals over filter-clean positions,
    /// indexed [Edit, Ctrl1, Ctrl2].
    pub channels: [ChannelStat; 3],
}

/// Mismatch channel over an A/T reference column: the editing class plus
/// the two control classes counted with identical thresholds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Channel {
    /// Editing: A->G, T->C.
    Edit,
    /// Control 1: A->C, T->G.
    Ctrl1,
    /// Control 2: A->T, T->A.
    Ctrl2,
}

impl Channel {
    /// Alternate-base index (A, C, G, T order) for this channel over
    /// `ref_base`; None when the reference base is not A/T.
    pub fn alt_index(self, ref_base: u8) -> Option<usize> {
        Some(match (self, ref_base) {
            (Channel::Edit, b'A') => 2,
            (Channel::Edit, b'T') => 1,
            (Channel::Ctrl1, b'A') => 1,
            (Channel::Ctrl1, b'T') => 2,
            (Channel::Ctrl2, b'A') => 3,
            (Channel::Ctrl2, b'T') => 0,
            _ => return None,
        })
    }
    /// Bed column label.
    pub fn label(self) -> &'static str {
        match self {
            Channel::Edit => "E",
            Channel::Ctrl1 => "C1",
            Channel::Ctrl2 => "C2",
        }
    }
}

/// Per-channel totals over filter-clean positions.
#[derive(Debug, Clone, Copy, Default)]
pub struct ChannelStat {
    /// Candidate positions passing every gate.
    pub sites: usize,
    /// Alternate-base reads summed over those positions.
    pub alt_reads: u64,
}

pub fn run_call(params: &CallParams) -> Result<CallStats, CallError> {
    call::run_call(params)
}
