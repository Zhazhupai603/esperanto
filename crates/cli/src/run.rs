//! `esperanto run` — full pipeline (→ `esperanto_flow::run_pipeline`).

use std::path::PathBuf;

use clap::Args;
use esperanto_scan::LibType;

#[derive(Args)]
pub struct RunArgs {
    /// R1 FASTQ inputs (repeatable / comma-separated; FASTQ entry).
    #[arg(long, value_delimiter = ',')]
    r1: Vec<PathBuf>,
    /// R2 FASTQ inputs (empty = single-end).
    #[arg(long, value_delimiter = ',')]
    r2: Vec<PathBuf>,
    /// Input BAM (Bam entry; must be coordinate-sorted + indexed).
    #[arg(long)]
    bam: Option<PathBuf>,
    /// User sites file `chrom<TAB>pos` (BamSites entry).
    #[arg(long)]
    sites: Option<PathBuf>,
    /// paidx index (required for FASTQ entry).
    #[arg(long)]
    index: Option<PathBuf>,
    /// Reference FASTA (default: refs discovery).
    #[arg(long)]
    fasta: Option<PathBuf>,
    /// Optional GTF (default: refs discovery).
    #[arg(long)]
    gtf: Option<PathBuf>,
    /// gnomAD VCF (default: refs discovery).
    #[arg(long)]
    gnomad: Option<PathBuf>,
    /// Model bundle root (default: zero-config 5-level fallback).
    #[arg(long)]
    bundle: Option<PathBuf>,
    /// Caduceus encoder dir (default: resolved inside the bundle).
    #[arg(long)]
    caduceus: Option<PathBuf>,
    /// Optional L1 engine bundle.
    #[arg(long)]
    l1_bundle: Option<PathBuf>,
    /// Library strandedness.
    #[arg(long, default_value = "unstranded", value_parser = ["unstranded", "stranded"])]
    lib: String,
    /// Output root directory.
    #[arg(long)]
    out: PathBuf,
    /// Worker threads (0 = all cores).
    #[arg(long, default_value_t = 0)]
    threads: usize,
    /// score batch size.
    #[arg(long, default_value_t = 64)]
    batch: usize,
}

pub fn run(a: RunArgs) -> anyhow::Result<()> {
    let bundle = crate::resolve::bundle(&a.bundle)?;
    let l1_bundle = a
        .index
        .as_ref()
        .and_then(|idx| crate::resolve::l1_bundle(&a.l1_bundle, idx));
    let mut fasta = a.fasta;
    let mut gtf = a.gtf;
    let mut gnomad = a.gnomad;
    if fasta.is_none() || gtf.is_none() || gnomad.is_none() {
        if let Some(refs) = crate::resolve::refs() {
            refs.fill_fasta(&mut fasta);
            refs.fill_gtf(&mut gtf);
            refs.fill_gnomad(&mut gnomad);
        }
    }
    let fasta = crate::resolve::require_fasta(fasta)?;
    let params = esperanto_flow::RunParams {
        r1: a.r1,
        r2: a.r2,
        bam: a.bam,
        sites: a.sites,
        index: a.index,
        fasta,
        gtf,
        gnomad,
        bundle,
        caduceus: a.caduceus,
        l1_bundle,
        lib: if a.lib == "stranded" {
            LibType::Stranded
        } else {
            LibType::Unstranded
        },
        out_dir: a.out,
        threads: crate::resolve::threads(a.threads),
        batch: a.batch,
    };
    esperanto_flow::run_pipeline(&params)?;
    Ok(())
}
