//! `esperanto resume` — continue an interrupted run from the first broken
//! stage (params and input fingerprints come from run.json, no flags needed).

use std::path::PathBuf;

use clap::Args;
use esperanto_flow::resume;
use esperanto_flow::{DeviceAsk, RunParams};
use esperanto_scan::LibType;
use esperanto_score::pipeline::DeviceChoice;

#[derive(Args)]
pub struct ResumeArgs {
    /// Run directory (the sample directory holding run.json).
    dir: PathBuf,
}

pub fn run(a: ResumeArgs) -> anyhow::Result<()> {
    let doc = resume::read_run_json(&a.dir)?;
    resume::verify_inputs(&doc)?;
    resume::verify_param_paths(&doc)?;
    let entry = resume::parse_entry(&doc.entry)?;
    if doc.esperanto != env!("CARGO_PKG_VERSION") {
        eprintln!(
            "[resume] note: run started with esperanto {}, resuming with {}",
            doc.esperanto,
            env!("CARGO_PKG_VERSION")
        );
    }
    let paths = |role: &str| -> Vec<PathBuf> {
        doc.inputs
            .iter()
            .filter(|i| i.role == role)
            .map(|i| i.path.clone())
            .collect()
    };
    let params = RunParams {
        r1: paths("r1"),
        r2: paths("r2"),
        bam: paths("bam").into_iter().next(),
        sites: paths("sites").into_iter().next(),
        index: doc.params.index.clone(),
        fasta: doc.params.fasta.clone(),
        gtf: doc.params.gtf.clone(),
        gnomad: doc.params.gnomad.clone(),
        bundle: doc.params.bundle.clone(),
        caduceus: doc.params.caduceus.clone(),
        l1_bundle: doc.params.l1_bundle.clone(),
        lib: if doc.params.lib == "stranded" {
            LibType::Stranded
        } else {
            LibType::Unstranded
        },
        out_dir: a.dir.clone(),
        threads: doc.params.threads,
        batch: doc.params.batch,
        device: match doc.params.device.as_str() {
            "cpu" => DeviceChoice::Cpu,
            "gpu" => DeviceChoice::Gpu,
            _ => DeviceChoice::Auto,
        },
        device_ask: Some(DeviceAsk::new(crate::confirm::ask_use_gpu)),
    };
    esperanto_flow::stages::preflight(&params, entry)?;
    // Lock before any artifact reads: a concurrent run must not race the walk.
    let _lock = resume::RunLock::acquire(&a.dir)?;
let has_cpaidx = params
        .index
        .as_ref()
        .is_some_and(|i| i.with_extension("cpaidx").is_file());
    let start = resume::walk(&a.dir, entry, has_cpaidx, params.gtf.is_some());
    let report_done = params.gtf.is_none() || {
        let sample = a
            .dir
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        resume::valid_nonempty(&a.dir.join(format!("{sample}.report.html")))
    };
    if start == resume::Stage::Report && report_done {
        eprintln!("[resume] nothing to do: all stage artifacts are valid");
        return Ok(());
    }
    let skipped: Vec<&str> = [
        resume::Stage::Qc,
        resume::Stage::Map,
        resume::Stage::Rescue,
        resume::Stage::Sort,
        resume::Stage::Scan,
        resume::Stage::Score,
        resume::Stage::Vcf,
    ]
    .into_iter()
    .filter(|s| *s < start && *s >= resume::Stage::first_for(entry))
    .map(|s| s.name())
    .collect();
    eprintln!(
        "[resume] {}: reusing {} stage(s); continuing at {}",
        doc.sample,
        skipped.join("/"),
        start.name()
    );
    resume::clean_stage(&a.dir, start);
    esperanto_flow::run_from(&params, entry, start)?;
    eprintln!("[resume] done: {}", a.dir.display());
    Ok(())
}
