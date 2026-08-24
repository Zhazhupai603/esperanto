//! `esperanto scan` — candidate editing-site discovery (→
//! `esperanto_scan::run_call`).

use std::path::PathBuf;

use clap::Args;
use esperanto_scan::{CallParams, LibType};

#[derive(Args)]
pub struct ScanArgs {
    /// Input BAM.
    #[arg(long)]
    bam: PathBuf,
    /// .baln dual-source channel (when given, takes priority over --bam).
    #[arg(long)]
    baln: Option<PathBuf>,
    /// Output candidates.bed.
    #[arg(long)]
    out: PathBuf,
    /// Reference FASTA (default: refs discovery; without any, majority
    /// pseudo-reference is used per scan contract).
    #[arg(long)]
    fasta: Option<PathBuf>,
    /// Optional GTF (strand evidence).
    #[arg(long)]
    gtf: Option<PathBuf>,
    /// gnomAD VCF (soft down-weight; default: refs discovery).
    #[arg(long)]
    gnomad: Option<PathBuf>,
    /// Library strandedness.
    #[arg(long, default_value = "unstranded", value_parser = ["unstranded", "stranded"])]
    lib: String,
    /// Enable C>U symmetric mode.
    #[arg(long)]
    enable_cu: bool,
    /// Minimum call score (marking only, per scan contract).
    #[arg(long)]
    min_call_score: Option<f64>,
    /// Override scoring spec JSON (builtin v2 default).
    #[arg(long)]
    spec: Option<PathBuf>,
    /// Worker threads (0 = all cores).
    #[arg(long, default_value_t = 0)]
    threads: usize,
}

pub fn run(a: ScanArgs) -> anyhow::Result<()> {
    let mut fasta = a.fasta;
    let mut gnomad = a.gnomad;
    if fasta.is_none() || gnomad.is_none() {
        if let Some(refs) = crate::resolve::refs() {
            refs.fill_fasta(&mut fasta);
            refs.fill_gnomad(&mut gnomad);
        }
    }
    let cp = CallParams {
        bam: a.bam,
        out: a.out,
        fasta,
        gtf: a.gtf,
        gnomad,
        lib: if a.lib == "stranded" {
            LibType::Stranded
        } else {
            LibType::Unstranded
        },
        enable_cu: a.enable_cu,
        min_call_score: a.min_call_score,
        spec: a.spec,
        threads: crate::resolve::threads(a.threads),
        baln: a.baln,
    };
    let stats = esperanto_scan::run_call(&cp)?;
    eprintln!("[scan] {} candidates", stats.candidates);
    Ok(())
}
