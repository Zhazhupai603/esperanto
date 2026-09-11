//! `esperanto map` — read alignment (→ `esperanto_flow::map_stage`; RNA
//! 2-pass wiring lives in the flow crate, cli is a thin front).

use std::path::PathBuf;

use clap::Args;

#[derive(Args)]
pub struct MapArgs {
    /// R1 FASTQ input.
    #[arg(long)]
    r1: PathBuf,
    /// R2 FASTQ input (omit for single-end).
    #[arg(long)]
    r2: Option<PathBuf>,
    /// paidx index path.
    #[arg(long)]
    index: PathBuf,
    /// Optional GTF (junction library).
    #[arg(long)]
    gtf: Option<PathBuf>,
    /// Optional L1 engine bundle.
    #[arg(long)]
    l1_bundle: Option<PathBuf>,
    /// Output directory (raw.bam / unmapped.fq.gz / align_qc.json / align.baln).
    #[arg(long)]
    out: PathBuf,
    /// Worker threads (0 = all cores).
    #[arg(long, default_value_t = 0)]
    threads: usize,
    /// Also run the negative-control (A>C / T>G) rescue passes, count-only.
    #[arg(long)]
    rescue_control: bool,
}

pub fn run(a: MapArgs) -> anyhow::Result<()> {
    let l1 = crate::resolve::l1_bundle(&a.l1_bundle, &a.index);
    let (raw, _index, _l1) = esperanto_flow::map_stage(
        &a.index,
        a.gtf.as_deref(),
        l1.as_deref(),
        &a.r1,
        a.r2.as_deref(),
        &a.out,
        crate::resolve::threads(a.threads),
        a.rescue_control,
    )?;
    eprintln!("[map] wrote {}", raw.display());
    Ok(())
}
