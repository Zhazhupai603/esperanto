//! `esperanto score` — RE_PROB scoring (→
//! `esperanto_score::pipeline::score_sites_batched`).

use std::fmt::Write as _;
use std::path::PathBuf;

use clap::Args;
use esperanto_score::pipeline as score_pipeline;

use crate::confirm::DeviceArg;

#[derive(Args)]
pub struct ScoreArgs {
    /// Input BAM (requires index).
    #[arg(long)]
    bam: PathBuf,
    /// Sites file `chrom<TAB>pos` (1-based, pos >= 1).
    #[arg(long)]
    sites: PathBuf,
    /// Reference FASTA (default: refs discovery).
    #[arg(long)]
    fasta: Option<PathBuf>,
    /// Optional .baln fast channel for the pileup pass (default: BAM).
    #[arg(long)]
    baln: Option<PathBuf>,
    /// Model bundle root (default: zero-config 5-level fallback).
    #[arg(long)]
    bundle: Option<PathBuf>,
    /// Caduceus encoder dir (default: resolved inside the bundle).
    #[arg(long)]
    caduceus: Option<PathBuf>,
    /// Output TSV `chrom<TAB>pos<TAB>prob`.
    #[arg(long)]
    out: PathBuf,
    /// Worker threads (0 = all cores).
    #[arg(long, default_value_t = 0)]
    threads: usize,
    /// Batch size.
    #[arg(long, default_value_t = 256)]
    batch: usize,
    /// Encoder device: auto (ask when a CUDA GPU is detected), cpu, or gpu.
    #[arg(long, value_enum, default_value_t = DeviceArg::Auto)]
    device: DeviceArg,
}

pub fn run(a: ScoreArgs) -> anyhow::Result<()> {
    let bundle = crate::resolve::bundle(&a.bundle)?;
    let mut fasta = a.fasta;
    if fasta.is_none() {
        if let Some(refs) = crate::resolve::refs() {
            refs.fill_fasta(&mut fasta);
        }
    }
    let fasta = crate::resolve::require_fasta(fasta)?;
    let caduceus = match a.caduceus {
        Some(c) => c,
        None => score_pipeline::resolve_encoder_from_bundle(&bundle)?,
    };
    let text = std::fs::read_to_string(&a.sites)?;
    let sites = score_pipeline::parse_sites(&text)?;
    let ask: fn() -> bool = crate::confirm::ask_use_gpu;
    let probs = score_pipeline::score_sites_batched(
        &a.bam,
        &fasta,
        &caduceus,
        &bundle,
        &sites,
        crate::resolve::threads(a.threads),
        a.batch.max(1),
        None,
        a.device.resolve(),
        Some(&ask),
        a.baln.as_deref(),
        score_pipeline::ReferenceCheck::Guardrail,
    )?;
    let mut buf = String::new();
    for ((chrom, pos), prob) in sites.iter().zip(&probs) {
        let _ = writeln!(buf, "{chrom}\t{pos}\t{prob}");
    }
    std::fs::write(&a.out, buf)?;
    eprintln!("[score] {} sites", probs.len());
    Ok(())
}
