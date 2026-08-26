//! `esperanto index` — build a paidx alignment index from a reference FASTA.

use std::path::PathBuf;

use clap::Args;
use esperanto_engine::Tidx;
use esperanto_map::seed::SeedParams;
use esperanto_map::{fasta, index, index_io};

#[derive(Args)]
pub struct IndexArgs {
    /// Reference FASTA.
    #[arg(long)]
    fasta: PathBuf,
    /// Optional transcript annotation GTF. When given, the L1 engine bundle
    /// (`<out stem>.bndl`) is built alongside the paidx index.
    #[arg(long)]
    gtf: Option<PathBuf>,
    /// Output paidx index path.
    #[arg(long)]
    out: PathBuf,
    /// K-mer length (must match the value used at alignment time, default 15).
    #[arg(long, default_value_t = 15)]
    k: u32,
    /// Window size in k-mers.
    #[arg(long, default_value_t = 5)]
    w: u32,
}

pub fn run(a: IndexArgs) -> anyhow::Result<()> {
    build_all(&a.fasta, a.gtf.as_deref(), &a.out, a.k, a.w)
}

/// Build the paidx index and — when a GTF is given — the L1 engine bundle
/// alongside it. Shared by `index` and `setup`.
pub(crate) fn build_all(
    fasta: &std::path::Path,
    gtf: Option<&std::path::Path>,
    out: &std::path::Path,
    k: u32,
    w: u32,
) -> anyhow::Result<()> {
    let reference =
        fasta::parse_fasta(fasta).map_err(|e| anyhow::anyhow!("parse fasta: {e}"))?;
    let idx = index::Index::build(reference, SeedParams { k, w });
    index_io::save(&idx, out).map_err(|e| anyhow::anyhow!("save paidx: {e}"))?;
    eprintln!("[index] wrote {}", out.display());
    if let Some(gtf) = gtf {
        // The L1 runtime needs BOTH the .bndl (projection + sequences) and its
        // .tidx sidecar (k-mer index) — keep all three artifacts together.
        let tidx_path = out.with_extension("tidx");
        esperanto_tidx::build(
            gtf,
            fasta,
            &tidx_path,
            &esperanto_tidx::BuildOptions::default(),
        )
        .map_err(|e| anyhow::anyhow!("build tidx: {e}"))?;
        let bndl = out.with_extension("bndl");
        let l1 = esperanto_engine::L1Index::build(&tidx_path, gtf, fasta)
            .map_err(|e| anyhow::anyhow!("build L1 bundle: {e}"))?;
        l1.save(&bndl)
            .map_err(|e| anyhow::anyhow!("save bndl: {e}"))?;
        eprintln!(
            "[index] wrote {} + {} ({} transcripts)",
            bndl.display(),
            tidx_path.display(),
            l1.tx_count()
        );
    }
    Ok(())
}
