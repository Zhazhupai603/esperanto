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

/// Collapse the A/G and T/C pairs (both cases) in sequence lines; headers
/// pass through. In collapsed space an A-to-I edit is an identity, so
/// hyperedited reads seed normally.
fn collapse_fasta_bytes(bytes: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(bytes.len());
    for line in bytes.split(|&b| b == b'\n') {
        if line.starts_with(b">") {
            out.extend_from_slice(line);
        } else {
            out.extend(line.iter().map(|&b| match b {
                b'A' | b'a' => b'G',
                b'G' | b'g' => b'G',
                b'C' | b'c' => b'C',
                b'T' | b't' => b'C',
                _ => b,
            }));
        }
        out.push(b'\n');
    }
    out
}

/// Build the paidx index and — when a GTF is given — the L1 engine bundle
/// alongside it. Always also builds the collapsed-alphabet rescue index
/// (`<out stem>.cpaidx`, k=31/w=10) used by the map-stage rescue pass.
/// Shared by `index` and `setup`.
pub(crate) fn build_all(
    fasta: &std::path::Path,
    gtf: Option<&std::path::Path>,
    out: &std::path::Path,
    k: u32,
    w: u32,
) -> anyhow::Result<()> {
    let bytes = std::fs::read(fasta).map_err(|e| anyhow::anyhow!("read fasta: {e}"))?;
    let reference =
        fasta::parse_fasta_bytes(&bytes).map_err(|e| anyhow::anyhow!("parse fasta: {e}"))?;
    let idx = index::Index::build(reference, SeedParams { k, w });
    index_io::save(&idx, out).map_err(|e| anyhow::anyhow!("save paidx: {e}"))?;
    eprintln!("[index] wrote {}", out.display());
    let creference = fasta::parse_fasta_bytes(&collapse_fasta_bytes(&bytes))
        .map_err(|e| anyhow::anyhow!("parse collapsed fasta: {e}"))?;
    let cidx = index::Index::build(creference, SeedParams { k: 31, w: 10 });
    let cpath = out.with_extension("cpaidx");
    index_io::save(&cidx, &cpath).map_err(|e| anyhow::anyhow!("save cpaidx: {e}"))?;
    eprintln!("[index] wrote {}", cpath.display());
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
