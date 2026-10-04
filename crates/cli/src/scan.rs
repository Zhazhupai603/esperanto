//! Candidate discovery using the unified molecular evidence reader.

use clap::Args;
use std::path::PathBuf;

#[derive(Args)]
pub struct ScanArgs {
    /// Input coordinate-sorted BAM with an index.
    #[arg(long)]
    bam: PathBuf,
    /// Output candidate BED.
    #[arg(long)]
    out: PathBuf,
    /// Reference FASTA (default: refs discovery).
    #[arg(long)]
    fasta: Option<PathBuf>,
    /// Worker threads (0 = all cores).
    #[arg(long, default_value_t = 0)]
    threads: usize,
}

pub fn run(a: ScanArgs) -> anyhow::Result<()> {
    let mut fasta = a.fasta;
    if let Some(refs) = crate::resolve::refs() {
        refs.fill_fasta(&mut fasta);
    }
    let fasta = crate::resolve::require_fasta(fasta)?;
    esperanto_flow::guard::check_species(&fasta)?;
    let count = esperanto_score::unified::discover_sites(
        &a.bam,
        &fasta,
        &a.out,
        crate::resolve::threads(a.threads),
    )?;
    eprintln!("[scan] {count} candidate sites");
    Ok(())
}
