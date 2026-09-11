//! `esperanto run` — full pipeline (→ `esperanto_flow::run_pipeline`).

use std::path::PathBuf;

use clap::Args;
use esperanto_flow::DeviceAsk;
use esperanto_scan::LibType;

use crate::confirm::DeviceArg;

#[derive(Args)]
pub struct RunArgs {
    /// R1 FASTQ inputs (repeatable / comma-separated; FASTQ entry).
    #[arg(long, value_delimiter = ',')]
    r1: Vec<PathBuf>,
    /// R2 FASTQ inputs (empty = single-end).
    #[arg(long, value_delimiter = ',')]
    r2: Vec<PathBuf>,
    /// Input BAM (Bam entry; must be coordinate-sorted + indexed).
    #[arg(long)]
    bam: Option<PathBuf>,
    /// User sites file `chrom<TAB>pos` (BamSites entry).
    #[arg(long)]
    sites: Option<PathBuf>,
    /// paidx index (required for FASTQ entry).
    #[arg(long)]
    index: Option<PathBuf>,
    /// Reference FASTA (default: refs discovery).
    #[arg(long)]
    fasta: Option<PathBuf>,
    /// Optional GTF (default: refs discovery).
    #[arg(long)]
    gtf: Option<PathBuf>,
    /// gnomAD VCF (default: refs discovery).
    #[arg(long)]
    gnomad: Option<PathBuf>,
    /// Model bundle root (default: zero-config 5-level fallback).
    #[arg(long)]
    bundle: Option<PathBuf>,
    /// Caduceus encoder dir (default: resolved inside the bundle).
    #[arg(long)]
    caduceus: Option<PathBuf>,
    /// Optional L1 engine bundle.
    #[arg(long)]
    l1_bundle: Option<PathBuf>,
    /// Library strandedness.
    #[arg(long, default_value = "unstranded", value_parser = ["unstranded", "stranded"])]
    lib: String,
    /// Output root directory.
#[arg(long)]
    out: PathBuf,
    /// Sample name for the run directory (default: derived from the R1/BAM file name).
    #[arg(long)]
    sample: Option<String>,
    /// Knock-in mouse run: splice these human genes onto the mouse reference
    /// (comma-separated symbols; bare --hybrid opens the interactive picker).
    #[arg(long, num_args = 0..=1, default_missing_value = "")]
    hybrid: Option<String>,
    /// Worker threads (0 = all cores).
    #[arg(long, default_value_t = 0)]
    threads: usize,
    /// score batch size.
    #[arg(long, default_value_t = 256)]
    batch: usize,
    /// Encoder device: auto (ask when a CUDA GPU is detected), cpu, or gpu.
    #[arg(long, value_enum, default_value_t = DeviceArg::Auto)]
    device: DeviceArg,
}

pub fn run(a: RunArgs) -> anyhow::Result<()> {
    let sample = a.sample.clone().unwrap_or_else(|| derive_sample(&a));
let bundle = crate::resolve::bundle(&a.bundle)?;
    let mut index = a.index;
    let mut fasta = a.fasta;
    let mut gtf = a.gtf;
    let mut gnomad = a.gnomad;
    let mut run_manifest: Option<esperanto_flow::manifest::SpeciesManifest> = None;
    if let Some(h) = &a.hybrid {
        let genes = if h.is_empty() {
            None
        } else {
            Some(h.split(',').map(|s| s.trim().to_string()).filter(|s| !s.is_empty()).collect())
        };
        let dir = crate::hybrid::resolve(genes)?;
        index = Some(dir.join("hybrid.paidx"));
        fasta = Some(dir.join("hybrid.fa"));
        let g = dir.join("hybrid.gtf");
        gtf = g.is_file().then_some(g);
        run_manifest = esperanto_flow::manifest::SpeciesManifest::read(&dir);
    }
    if index.is_none() || fasta.is_none() || gtf.is_none() || gnomad.is_none() {
        if let Some(refs) = crate::resolve::refs() {
            refs.fill_index(&mut index);
            refs.fill_fasta(&mut fasta);
            refs.fill_gtf(&mut gtf);
            refs.fill_gnomad(&mut gnomad);
        }
    }
    let l1_bundle = index
        .as_ref()
        .and_then(|idx| crate::resolve::l1_bundle(&a.l1_bundle, idx));
    let fasta = crate::resolve::require_fasta(fasta)?;
    let out_dir = a.out.join(&sample);
    if out_dir.join("run.json").exists() {
        anyhow::bail!(
            "run directory {} already exists; use `esperanto resume {}` to continue it",
            out_dir.display(),
            out_dir.display()
        );
    }
    std::fs::create_dir_all(&out_dir)?;
    if let Some(m) = &run_manifest {
        m.write(&out_dir)?; // preflight guardrail + report/score contig split read it here
    }
    let params = esperanto_flow::RunParams {
        r1: a.r1,
        r2: a.r2,
        bam: a.bam,
        sites: a.sites,
        index,
        fasta,
        gtf,
        gnomad,
        bundle,
        caduceus: a.caduceus,
        l1_bundle,
        lib: if a.lib == "stranded" {
            LibType::Stranded
        } else {
            LibType::Unstranded
        },
        out_dir: out_dir.clone(),
        threads: crate::resolve::threads(a.threads),
        batch: a.batch,
        device: a.device.resolve(),
        device_ask: Some(DeviceAsk::new(crate::confirm::ask_use_gpu)),
    };
    let _lock = esperanto_flow::resume::RunLock::acquire(&out_dir)?;
    esperanto_flow::resume::write_run_json(&out_dir, &sample, &params)?;
    esperanto_flow::run_pipeline(&params)?;
    eprintln!("[run] sample directory: {}", out_dir.display());
    Ok(())
}

/// Sample name from the first input file name: strip compression/format
/// suffixes, a trailing `.bam`, and a trailing mate marker (`_R1/_R2/_1/_2`).
fn derive_sample(a: &RunArgs) -> String {
    let src = a.r1.first().or(a.bam.as_ref());
    let Some(src) = src else { return "sample".to_string() };
    let mut name = src
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    while let Some(s) = name
        .strip_suffix(".gz")
        .or_else(|| name.strip_suffix(".fastq"))
        .or_else(|| name.strip_suffix(".fq"))
        .or_else(|| name.strip_suffix(".bam"))
    {
        name = s.to_string();
    }
    for suf in ["_R1", "_R2", "_1", "_2"] {
        if let Some(s) = name.strip_suffix(suf) {
            name = s.to_string();
            break;
        }
    }
    if name.is_empty() {
        "sample".to_string()
    } else {
        name
    }
}

#[cfg(test)]
mod tests {
    use super::derive_sample;
    use super::RunArgs;

    fn args(r1: &str) -> RunArgs {
        RunArgs {
            r1: vec![r1.into()],
            r2: Vec::new(),
            bam: None,
            sites: None,
            index: None,
            fasta: None,
            gtf: None,
            gnomad: None,
            bundle: None,
            caduceus: None,
            l1_bundle: None,
            lib: "unstranded".to_string(),
            out: "out".into(),
            sample: None,
            hybrid: None,
            threads: 0,
            batch: 64,
            device: crate::confirm::DeviceArg::Auto,
        }
    }

    #[test]
    fn derive_sample_strips_formats_and_mate_markers() {
        assert_eq!(derive_sample(&args("APOE4_070938_1_1.fq.gz")), "APOE4_070938_1");
        assert_eq!(derive_sample(&args("sample_R1.fastq.gz")), "sample");
        assert_eq!(derive_sample(&args("lane.fq")), "lane");
        let mut b = args("x");
        b.r1 = Vec::new();
        b.bam = Some("tumor.sorted.bam".into());
        assert_eq!(derive_sample(&b), "tumor.sorted");
    }
}
