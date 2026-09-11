//! Pipeline parameters and entry derivation (spec §entry derivation / §parameters).

use std::path::PathBuf;

use crate::FlowError;

/// Entry type, derived purely from which input fields are populated.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Entry {
    /// `--r1` only: qc → map → sort → scan → score → vcf.
    FastqSe,
    /// `--r1` + `--r2`: same stage sequence as FastqSe.
    FastqPe,
    /// `--bam` without `--sites`: scan → score → vcf.
    Bam,
    /// `--bam` + `--sites`: score → vcf.
    BamSites,
}

/// Full parameter set for one pipeline run.
#[derive(Clone, Debug)]
pub struct RunParams {
    /// R1 FASTQ inputs (qc merges multiple lanes in order).
    pub r1: Vec<PathBuf>,
    /// R2 FASTQ inputs; empty = single-end.
    pub r2: Vec<PathBuf>,
    /// Input BAM (Bam / BamSites entries; must be coordinate-sorted + indexed).
    pub bam: Option<PathBuf>,
    /// User sites file `chrom\tpos` 1-based (BamSites entry).
    pub sites: Option<PathBuf>,
    /// paidx index path (required for FASTQ entries).
    pub index: Option<PathBuf>,
    /// Reference FASTA (requires `<fasta>.fai`).
    pub fasta: PathBuf,
    /// Optional GTF (junction library for map; annotation for scan).
    pub gtf: Option<PathBuf>,
    /// Optional gnomAD VCF for scan.
    pub gnomad: Option<PathBuf>,
    /// score bundle root.
    pub bundle: PathBuf,
    /// Optional caduceus encoder path; defaults to bundle resolution.
    pub caduceus: Option<PathBuf>,
    /// Optional L1 engine bundle (`engine::L1Index::open`).
    pub l1_bundle: Option<PathBuf>,
    /// Library strandedness (passed through to scan).
    pub lib: esperanto_scan::LibType,
    /// Output root; stages write under `<out_dir>/<stage>/`.
    pub out_dir: PathBuf,
    /// Worker threads.
    pub threads: usize,
    /// score batch size (default 64).
    pub batch: usize,
}

impl RunParams {
    /// Derive the entry type from populated fields (pure field branching).
    pub fn entry(&self) -> Result<Entry, FlowError> {
        let has_r1 = !self.r1.is_empty();
        let has_r2 = !self.r2.is_empty();
        let has_bam = self.bam.is_some();
        let has_sites = self.sites.is_some();
        match (has_r1, has_r2, has_bam, has_sites) {
            (true, false, false, false) => Ok(Entry::FastqSe),
            (true, true, false, false) => Ok(Entry::FastqPe),
            (false, false, true, false) => Ok(Entry::Bam),
            (false, false, true, true) => Ok(Entry::BamSites),
            _ => Err(FlowError::Entry(format!(
                "invalid input combination (r1={}, r2={}, bam={}, sites={}); \
                 expected r1[/r2], bam, or bam+sites",
                self.r1.len(),
                self.r2.len(),
                has_bam,
                has_sites
            ))),
        }
    }
}
