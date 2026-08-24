//! `esperanto index` — build a paidx alignment index from a reference FASTA.

use std::path::PathBuf;

use clap::Args;
use esperanto_map::seed::SeedParams;
use esperanto_map::{fasta, index, index_io};

#[derive(Args)]
pub struct IndexArgs {
    /// Reference FASTA.
    #[arg(long)]
    fasta: PathBuf,
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
    let reference =
        fasta::parse_fasta(&a.fasta).map_err(|e| anyhow::anyhow!("parse fasta: {e}"))?;
    let idx = index::Index::build(reference, SeedParams { k: a.k, w: a.w });
    index_io::save(&idx, &a.out).map_err(|e| anyhow::anyhow!("save paidx: {e}"))?;
    eprintln!("[index] wrote {}", a.out.display());
    Ok(())
}
