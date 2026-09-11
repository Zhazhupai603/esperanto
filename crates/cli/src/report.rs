//! `esperanto report` - regenerate the standalone HTML report for a
//! finished run directory (-> `esperanto_report::generate`).

use std::path::PathBuf;

use anyhow::anyhow;
use clap::Args;

#[derive(Args)]
pub struct ReportArgs {
    /// Finished run directory (sites.vcf, qc/, map/).
    #[arg(long)]
    out: PathBuf,
    /// Reference FASTA (default: refs discovery).
    #[arg(long)]
    fasta: Option<PathBuf>,
    /// GTF annotation (default: refs discovery).
    #[arg(long)]
    gtf: Option<PathBuf>,
}

pub fn run(a: ReportArgs) -> anyhow::Result<()> {
    let mut fasta = a.fasta;
    let mut gtf = a.gtf;
    if fasta.is_none() || gtf.is_none() {
        if let Some(refs) = crate::resolve::refs() {
            refs.fill_fasta(&mut fasta);
            refs.fill_gtf(&mut gtf);
        }
    }
    let fasta = crate::resolve::require_fasta(fasta)?;
    let gtf = gtf.ok_or_else(|| {
        anyhow!("--gtf not given and no refs directory found (set ESPERANTO_REFS)")
    })?;
    let path = esperanto_report::generate(&a.out, &fasta, &gtf)?;
    println!("report: {}", path.display());
    Ok(())
}
