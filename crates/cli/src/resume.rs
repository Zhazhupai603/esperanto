//! `esperanto resume` — continue an interrupted run from the first broken
//! stage (params and input fingerprints come from run.json, no flags needed).

use std::path::PathBuf;

use clap::Args;
use esperanto_flow::resume;
use esperanto_flow::RunParams;

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
    let (params, ignored_snp) = params_of(&doc, &a.dir);
    if let Some(snp) = &ignored_snp {
        eprintln!(
            "[resume] warning: run.json freezes a known-SNP VCF ({}); known-SNP exclusion is an internal channel that is always off, so the setting is ignored",
            snp.display()
        );
    }
    esperanto_flow::stages::preflight(&params, entry)?;
    // Lock before any artifact reads: a concurrent run must not race the walk.
    let _lock = resume::RunLock::acquire(&a.dir)?;
let has_rescue_index = params.index.as_ref().is_some_and(|i| {
        i.with_extension("ag.paidx").is_file() || i.with_extension("tc.paidx").is_file()
    });
    let start = resume::walk(&a.dir, entry, has_rescue_index, params.gtf.is_some());
    let start = start_stage_after_legacy_snp(start, entry, ignored_snp.is_some());
    if ignored_snp.is_some() && start == resume::Stage::Scan {
        eprintln!(
            "[resume] scan-stage artifacts may have been frozen with the SNP filter on; rebuilding from scan"
        );
    }
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

/// Rebuild the frozen RunParams from run.json. Returns the params plus the
/// legacy known-SNP VCF path the run had frozen, when present: known-SNP
/// exclusion is an internal channel that is always off, so a setting frozen
/// by an older build is reported and dropped rather than re-enabled.
fn params_of(doc: &resume::RunJson, dir: &std::path::Path) -> (RunParams, Option<PathBuf>) {
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
        snp: None,
        bundle: doc.params.bundle.clone(),
        l1_bundle: doc.params.l1_bundle.clone(),
        out_dir: dir.to_path_buf(),
        threads: doc.params.threads,
        batch: doc.params.batch,
        // Frozen runs predate the control pass; a resumed run completes the
        // pass set it started with.
        rescue_control: false,
        // run.json does not freeze --stop-after: it is a one-off inspection
        // switch; `resume` always means "finish the run".
        stop_after: None,
    };
    (params, doc.params.snp.clone())
}

/// A run frozen while the known-SNP filter was on may hold scan-stage
/// artifacts produced under it. Resume rebuilds from scan onward in that
/// case (the filter is now always off) instead of reusing the filtered
/// cache; invalidation beats refusal here because the cascade cleanup
/// already exists and a scan rerun is cheap next to the mapping stages a
/// refusal would force the user to redo by hand. The BamSites entry has no
/// scan stage, so nothing can be stale there.
fn start_stage_after_legacy_snp(
    start: resume::Stage,
    entry: esperanto_flow::Entry,
    legacy_snp: bool,
) -> resume::Stage {
    if legacy_snp && entry != esperanto_flow::Entry::BamSites && start > resume::Stage::Scan {
        resume::Stage::Scan
    } else {
        start
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn doc_with_snp(snp: Option<PathBuf>) -> resume::RunJson {
        resume::RunJson {
            version: resume::RUN_JSON_VERSION,
            esperanto: env!("CARGO_PKG_VERSION").to_string(),
            sample: "s1".to_string(),
            entry: "bam-sites".to_string(),
            inputs: Vec::new(),
            params: resume::ParamsJson {
                index: None,
                fasta: PathBuf::from("ref.fa"),
                gtf: None,
                snp,
                bundle: PathBuf::from("bundle"),
                l1_bundle: None,
                threads: 4,
                batch: 64,
            },
        }
    }

    #[test]
    fn legacy_snp_setting_is_ignored_and_reported() {
        let doc = doc_with_snp(Some(PathBuf::from("legacy.vcf")));
        let (params, ignored) = params_of(&doc, std::path::Path::new("run-dir"));
        assert!(
            params.snp.is_none(),
            "resumed runs never re-enable SNP exclusion"
        );
        assert_eq!(ignored, Some(PathBuf::from("legacy.vcf")));
    }

    #[test]
    fn snp_free_run_json_resumes_without_a_warning_payload() {
        let doc = doc_with_snp(None);
        let (params, ignored) = params_of(&doc, std::path::Path::new("run-dir"));
        assert!(params.snp.is_none());
        assert!(ignored.is_none());
        assert_eq!(params.threads, 4);
        assert_eq!(params.batch, 64);
    }

    #[test]
    fn legacy_snp_invalidates_cached_scan_artifacts() {
        use esperanto_flow::resume::Stage;
        // A run interrupted after scan, verify, or completion restarts at scan.
        for later in [Stage::Sort, Stage::Scan, Stage::Verify, Stage::Report] {
            let start = start_stage_after_legacy_snp(later, esperanto_flow::Entry::Bam, true);
            if later > Stage::Scan {
                assert_eq!(start, Stage::Scan, "{later:?} must rebuild from scan");
            } else {
                assert_eq!(start, later, "{later:?} already rebuilds scan");
            }
        }
        // Fastq entries behave the same way.
        assert_eq!(
            start_stage_after_legacy_snp(Stage::Report, esperanto_flow::Entry::FastqSe, true),
            Stage::Scan
        );
        // No legacy setting: the walk's stage stands.
        assert_eq!(
            start_stage_after_legacy_snp(Stage::Verify, esperanto_flow::Entry::Bam, false),
            Stage::Verify
        );
        // BamSites has no scan stage: nothing to invalidate.
        assert_eq!(
            start_stage_after_legacy_snp(Stage::Report, esperanto_flow::Entry::BamSites, true),
            Stage::Report
        );
    }
}
