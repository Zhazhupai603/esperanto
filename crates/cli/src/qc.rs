//! `esperanto qc` — FASTQ quality control (→ `esperanto_qc::run`).

use std::path::PathBuf;

use clap::Args;
use esperanto_qc::{OutFormat, PolygMode, QcParams};

#[derive(Args)]
pub struct QcArgs {
    /// R1 FASTQ inputs (repeatable / comma-separated; lanes merged in order).
    #[arg(long, required = true, value_delimiter = ',')]
    r1: Vec<PathBuf>,
    /// R2 FASTQ inputs (empty = single-end; must match r1 count).
    #[arg(long, value_delimiter = ',')]
    r2: Vec<PathBuf>,
    /// Output directory.
    #[arg(long)]
    out: PathBuf,
    /// Override R1 adapter table.
    #[arg(long)]
    adapter_r1: Vec<String>,
    /// Override R2 adapter table.
    #[arg(long)]
    adapter_r2: Vec<String>,
    /// Disable adapter trimming.
    #[arg(long)]
    disable_adapter_trim: bool,
    /// Disable PE overlap trimming.
    #[arg(long)]
    disable_pe_overlap: bool,
    /// Enable BWA-style 3' quality trimming.
    #[arg(long)]
    qtrim: bool,
    /// qtrim quality cutoff.
    #[arg(long, default_value_t = 20)]
    qtrim_cutoff: u8,
    /// Fixed bases to cut from R1 front.
    #[arg(long, default_value_t = 0)]
    trim_front1: usize,
    /// Fixed bases to cut from R1 tail.
    #[arg(long, default_value_t = 0)]
    trim_tail1: usize,
    /// Fixed bases to cut from R2 front.
    #[arg(long, default_value_t = 0)]
    trim_front2: usize,
    /// Fixed bases to cut from R2 tail.
    #[arg(long, default_value_t = 0)]
    trim_tail2: usize,
    /// polyG mode.
    #[arg(long, default_value = "auto", value_parser = ["auto", "on", "off"])]
    polyg: String,
    /// Minimum read length.
    #[arg(long, default_value_t = 15)]
    min_len: usize,
    /// Maximum N count.
    #[arg(long, default_value_t = 5)]
    n_max: usize,
    /// Maximum fraction of Q<15 bases.
    #[arg(long, default_value_t = 0.4)]
    q15_frac_max: f64,
    /// Keep the passing end of a failed pair.
    #[arg(long)]
    keep_unpaired: bool,
    /// Auto-detect SE adapter when the table misses.
    #[arg(long)]
    detect_adapter_se: bool,
    /// Worker threads (0 = all cores).
    #[arg(long, default_value_t = 0)]
    threads: usize,
    /// Output format.
    #[arg(long, default_value = "fqgz", value_parser = ["fqgz", "bfq"])]
    format: String,
}

pub fn run(a: QcArgs) -> anyhow::Result<()> {
    let mut p = QcParams {
        r1: a.r1,
        r2: a.r2,
        out_dir: a.out,
        ..Default::default()
    };
    if !a.adapter_r1.is_empty() {
        p.adapters_r1 = a.adapter_r1;
    }
    if !a.adapter_r2.is_empty() {
        p.adapters_r2 = a.adapter_r2;
    }
    p.adapter_trim = !a.disable_adapter_trim;
    p.pe_overlap = !a.disable_pe_overlap;
    p.qtrim = a.qtrim;
    p.qtrim_cutoff = a.qtrim_cutoff;
    p.trim_front1 = a.trim_front1;
    p.trim_tail1 = a.trim_tail1;
    p.trim_front2 = a.trim_front2;
    p.trim_tail2 = a.trim_tail2;
    p.polyg = match a.polyg.as_str() {
        "on" => PolygMode::On,
        "off" => PolygMode::Off,
        _ => PolygMode::Auto,
    };
    p.min_len = a.min_len;
    p.n_max = a.n_max;
    p.q15_frac_max = a.q15_frac_max;
    p.keep_unpaired = a.keep_unpaired;
    p.detect_adapter_se = a.detect_adapter_se;
    p.threads = crate::resolve::threads(a.threads);
    p.out_format = if a.format == "bfq" {
        OutFormat::Bfq
    } else {
        OutFormat::Fqgz
    };
    esperanto_qc::run(&p)?;
    Ok(())
}
