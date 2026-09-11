//! `esperanto scan` — A-to-I candidate site discovery (→
//! `esperanto_scan::run_call`).

use std::path::PathBuf;

use clap::Args;
use esperanto_scan::CallParams;

#[derive(Args)]
pub struct ScanArgs {
    /// Input BAM (coordinate-sorted, indexed).
    #[arg(long)]
    bam: PathBuf,
    /// Output candidates.bed.
    #[arg(long)]
    out: PathBuf,
    /// Reference FASTA (default: refs discovery).
    #[arg(long)]
    fasta: Option<PathBuf>,
    /// Known-SNP VCF; listed positions are removed (default: refs discovery;
    /// `--snp none` disables).
    #[arg(long)]
    snp: Option<String>,
    /// Reads below this MAPQ do not count.
    #[arg(long, default_value_t = 20)]
    min_mapq: u8,
    /// Bases below this Phred quality do not count.
    #[arg(long, default_value_t = 30)]
    min_bq: u8,
    /// Minimum counted depth.
    #[arg(long, default_value_t = 10)]
    min_depth: u32,
    /// Minimum reads carrying the editing-consistent alternate base.
    #[arg(long, default_value_t = 2)]
    min_edit_reads: u32,
    /// Minimum editing-consistent alternate frequency.
    #[arg(long, default_value_t = 0.01)]
    min_edit_freq: f64,
    /// Worker threads (0 = all cores).
    #[arg(long, default_value_t = 0)]
    threads: usize,
}

pub fn run(a: ScanArgs) -> anyhow::Result<()> {
    let mut fasta = a.fasta;
    let mut snp: Option<PathBuf> = match a.snp.as_deref() {
        Some("none") => None,
        Some(p) => Some(PathBuf::from(p)),
        None => None,
    };
    let resolve_snp = a.snp.is_none();
    if fasta.is_none() || resolve_snp {
        if let Some(refs) = crate::resolve::refs() {
            refs.fill_fasta(&mut fasta);
            if resolve_snp {
                refs.fill_snp(&mut snp);
            }
        }
    }
    let fasta = crate::resolve::require_fasta(fasta)?;
    let cp = CallParams {
        bam: a.bam,
        out: a.out,
        fasta,
        snp,
        min_mapq: a.min_mapq,
        min_bq: a.min_bq,
        min_depth: a.min_depth,
        min_edit_reads: a.min_edit_reads,
        min_edit_freq: a.min_edit_freq,
        threads: crate::resolve::threads(a.threads),
    };
    let stats = esperanto_scan::run_call(&cp)?;
    eprintln!(
        "[scan] {} candidates ({} below threshold, {} known SNPs removed)",
        stats.candidates, stats.below_threshold, stats.snp_removed
    );
    let [e, c1, c2] = stats.channels;
    let ctrl_mean = (c1.sites + c2.sites) as f64 / 2.0;
    let (snr, fdr) = channel_snr_fdr(e.sites, ctrl_mean);
    eprintln!("[scan] channel S/N {snr:.2}x (implied FDR {:.1}%)", fdr * 100.0);
    Ok(())
}

/// S/N and implied FDR from the three mismatch channels.
pub fn channel_snr_fdr(editing: usize, control_mean: f64) -> (f64, f64) {
    let snr = if control_mean > 0.0 {
        editing as f64 / control_mean
    } else {
        f64::INFINITY
    };
    let fdr = if editing > 0 {
        control_mean / editing as f64
    } else {
        0.0
    };
    (snr, fdr)
}
