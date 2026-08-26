//! Stage wiring (spec §stage wiring): contract paths between crates, artifacts
//! under `<out>/<stage>/`.

use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use esperanto_bamio::sort::{coordinate_sort, SortOptions};
use esperanto_engine::L1Index;
use esperanto_map::align::AlignConfig;
use esperanto_map::index::Index;
use esperanto_map::pipeline::{run_pe_2pass, run_se_2pass, PipelineOut};
use esperanto_map::{gtf, index_io};
use esperanto_qc::{OutFormat, QcParams};
use esperanto_scan::CallParams;
use esperanto_score::pipeline as score_pipeline;

use crate::params::{Entry, RunParams};
use crate::{guard, vcf, FlowError};

/// Wrap a downstream error without swallowing its semantics.
fn stage_err<E>(stage: &'static str) -> impl FnOnce(E) -> FlowError
where
    E: std::error::Error + Send + Sync + 'static,
{
    move |source| FlowError::Stage {
        stage,
        source: Box::new(source),
    }
}

/// anyhow-flavored variant for the score crate.
fn anyhow_stage_err(stage: &'static str) -> impl FnOnce(anyhow::Error) -> FlowError {
    move |source| FlowError::Stage {
        stage,
        source: source.into(),
    }
}

/// Top-level orchestration: entry derivation and every early failure happens
/// before `out_dir` is created.
pub fn run(params: &RunParams) -> Result<(), FlowError> {
    let entry = params.entry()?;
    guard::check_species(&params.fasta)?;
    match entry {
        Entry::FastqSe | Entry::FastqPe => {
            if params.index.is_none() {
                return Err(FlowError::Entry(
                    "FASTQ entry requires --index (paidx)".into(),
                ));
            }
        }
        Entry::Bam | Entry::BamSites => {
            if let Some(bam) = &params.bam {
                check_bam_index(bam)?;
            }
        }
    }
    fs::create_dir_all(&params.out_dir)?;

    let mut current_bam: PathBuf;
    let mut baln: Option<PathBuf> = None;
    match entry {
        Entry::FastqSe => {
            let (clean1, _) = stage_qc(params, entry)?;
            current_bam = stage_map(params, entry, &clean1, None)?;
            current_bam = stage_sort(&current_bam, params)?;
            // The SE mapper writes .baln; use it as the scan fast channel.
            baln = Some(params.out_dir.join("map").join("align.baln"));
        }
        Entry::FastqPe => {
            let (clean1, clean2) = stage_qc(params, entry)?;
            current_bam = stage_map(params, entry, &clean1, clean2.as_deref())?;
            current_bam = stage_sort(&current_bam, params)?;
            // The PE mapper does not write .baln (legacy parity); scan reads the BAM.
            baln = None;
        }
        Entry::Bam | Entry::BamSites => {
            if let Some(bam) = &params.bam {
                current_bam = bam.clone();
            } else {
                return Err(FlowError::Entry("Bam entry requires --bam".into()));
            }
        }
    }

    let sites = match entry {
        Entry::BamSites => match &params.sites {
            Some(p) => read_sites_file(p)?,
            None => return Err(FlowError::Entry("BamSites entry requires --sites".into())),
        },
        _ => {
            let bed = stage_scan(params, entry, &current_bam, baln)?;
            let fstats = crate::filter::CallFilter::default().apply_to_bed(&bed)?;
            eprintln!(
                "[filter] {} candidates -> {} kept (low_depth {} no_signal {})",
                fstats.input, fstats.kept, fstats.low_depth, fstats.no_signal
            );
            bed_to_sites(&bed)?
        }
    };

    let probs = stage_score(params, &current_bam, &sites)?;
    vcf::write_vcf(params, entry, &sites, &probs)?;
    Ok(())
}

/// Bam/BamSites entry contract: input BAM must already be coordinate-sorted
/// with a `.bai`/`.csi` alongside; flow never re-sorts user input.
fn check_bam_index(bam: &Path) -> Result<(), FlowError> {
    let candidates = [
        PathBuf::from(format!("{}.bai", bam.display())),
        bam.with_extension("bai"),
        PathBuf::from(format!("{}.csi", bam.display())),
        bam.with_extension("csi"),
    ];
    if candidates.iter().any(|p| p.exists()) {
        Ok(())
    } else {
        Err(FlowError::MissingBamIndex {
            path: bam.display().to_string(),
        })
    }
}

/// qc naming contract (spec: `<stem>.clean[_R1/_R2].fq.gz`), stem rules
/// mirrored from the qc crate.
fn stem_of(path: &Path) -> String {
    let mut name = path
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    loop {
        if name.ends_with(".gz") {
            name.truncate(name.len() - 3);
        } else if name.ends_with(".fastq") {
            name.truncate(name.len() - 6);
        } else if name.ends_with(".fq") {
            name.truncate(name.len() - 3);
        } else {
            break;
        }
    }
    name
}

fn stage_qc(params: &RunParams, entry: Entry) -> Result<(PathBuf, Option<PathBuf>), FlowError> {
    eprintln!("[qc] running");
    let dir = params.out_dir.join("qc");
    let qp = QcParams {
        r1: params.r1.clone(),
        r2: params.r2.clone(),
        out_dir: dir.clone(),
        out_format: OutFormat::Fqgz,
        threads: params.threads,
        ..Default::default()
    };
    esperanto_qc::run(&qp).map_err(stage_err("qc"))?;
    let stem1 = stem_of(&params.r1[0]);
    match entry {
        Entry::FastqSe => Ok((dir.join(format!("{stem1}.clean.fq.gz")), None)),
        _ => {
            let stem2 = stem_of(&params.r2[0]);
            Ok((
                dir.join(format!("{stem1}.clean_R1.fq.gz")),
                Some(dir.join(format!("{stem2}.clean_R2.fq.gz"))),
            ))
        }
    }
}

/// contig name → SQ id for the junction library; unknown contigs map to
/// `u32::MAX` (their junctions can never match a real alignment).
fn contig_id(index: &Index, name: &str) -> u32 {
    index
        .reference
        .contigs
        .iter()
        .position(|c| c.name == name)
        .map(|i| i as u32)
        .unwrap_or(u32::MAX)
}

fn stage_map(
    params: &RunParams,
    entry: Entry,
    clean1: &Path,
    clean2: Option<&Path>,
) -> Result<PathBuf, FlowError> {
    eprintln!("[map] running");
    let dir = params.out_dir.join("map");
    let paidx = match &params.index {
        Some(p) => p.clone(),
        None => return Err(FlowError::Entry("FASTQ entry requires --index".into())),
    };
    let r2 = match (entry, clean2) {
        (Entry::FastqSe, _) => None,
        (Entry::FastqPe, Some(c2)) => Some(c2),
        _ => {
            return Err(FlowError::Entry(
                "PE entry requires clean R2 from qc".into(),
            ))
        }
    };
    map_stage(
        &paidx,
        params.gtf.as_deref(),
        params.l1_bundle.as_deref(),
        clean1,
        r2,
        &dir,
        params.threads,
    )
}

/// Map stage wiring shared with the cli `map` subcommand: RNA 2-pass
/// (editing-aware), optional sjdb + L1, contract artifact names under
/// `out_dir`. Returns the raw (unsorted) BAM path.
pub fn map_stage(
    paidx: &Path,
    gtf_path: Option<&Path>,
    l1_bundle: Option<&Path>,
    r1: &Path,
    r2: Option<&Path>,
    out_dir: &Path,
    threads: usize,
) -> Result<PathBuf, FlowError> {
    let dir = out_dir;
    fs::create_dir_all(dir)?;
    let index = index_io::load(paidx).map_err(stage_err("map"))?;
    let mut config = AlignConfig::rna_default();
    config.extend.editing_aware = true;
    let jlib = match gtf_path {
        Some(g) => Some(Arc::new(
            gtf::from_gtf(g, |name| contig_id(&index, name)).map_err(stage_err("map"))?,
        )),
        None => None,
    };
    let l1 = match l1_bundle {
        Some(p) => {
            let l1 = L1Index::open(p).map_err(stage_err("map"))?;
            // Reference agreement guard: every L1 projection contig must
            // exist in the paidx reference — a mismatched bundle/index pair
            // would project reads onto wrong contigs.
            for name in l1.txmap().contigs() {
                if index.reference.contig_index(name.as_bytes()).is_none() {
                    return Err(FlowError::Entry(format!(
                        "L1 bundle contig '{name}' is absent from the index reference; the .bndl/.tidx pair must be built from the same reference as the paidx"
                    )));
                }
            }
            Some(Arc::new(l1))
        }
        None => None,
    };

    let raw_bam = dir.join("raw.bam");
    let mut out = PipelineOut {
        bam: Some(Box::new(File::create(&raw_bam)?)),
        unmapped_fq: Box::new(File::create(dir.join("unmapped.fq.gz"))?),
        index: &index,
        config,
        jlib,
        jkmer: None,
        l1,
        baln: Some(Box::new(File::create(dir.join("align.baln"))?)),
    };
    let stats = match r2 {
        None => run_se_2pass(&mut out, r1, threads),
        Some(c2) => run_pe_2pass(&mut out, r1, c2, threads),
    }
    .map_err(stage_err("map"))?
    .0;

    let json = serde_json::to_string_pretty(&stats).map_err(stage_err("map"))?;
    fs::write(dir.join("align_qc.json"), format!("{json}\n"))?;
    eprintln!(
        "[map] {} / {} mapped",
        stats.mapped_reads, stats.total_reads
    );
    Ok(raw_bam)
}

fn stage_sort(raw: &Path, params: &RunParams) -> Result<PathBuf, FlowError> {
    let sorted = raw.with_file_name("sorted.bam");
    let opts = SortOptions {
        threads: params.threads,
        ..Default::default()
    };
    let stats = coordinate_sort(raw, &sorted, &opts).map_err(stage_err("sort"))?;
    eprintln!(
        "[sort] {} records, {} chunk(s)",
        stats.records, stats.chunks
    );
    Ok(sorted)
}

fn stage_scan(
    params: &RunParams,
    entry: Entry,
    bam: &Path,
    baln: Option<PathBuf>,
) -> Result<PathBuf, FlowError> {
    eprintln!("[scan] running");
    let dir = params.out_dir.join("scan");
    fs::create_dir_all(&dir)?;
    let out = dir.join("candidates.bed");
    let cp = CallParams {
        bam: bam.to_path_buf(),
        out: out.clone(),
        fasta: Some(params.fasta.clone()),
        gtf: params.gtf.clone(),
        gnomad: params.gnomad.clone(),
        lib: params.lib,
        enable_cu: false,
        min_call_score: None,
        spec: None,
        threads: params.threads,
        baln: match entry {
            Entry::Bam | Entry::BamSites => None,
            _ => baln,
        },
    };
    let stats = esperanto_scan::run_call(&cp).map_err(stage_err("scan"))?;
    eprintln!("[scan] {} candidates", stats.candidates);
    Ok(out)
}

/// candidates.bed → sites bridge (spec: col 1 = chrom, col 3 = 1-based pos).
pub fn bed_to_sites(path: &Path) -> Result<Vec<(String, i64)>, FlowError> {
    let text = fs::read_to_string(path)?;
    let mut sites = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        let mut cols = line.split('\t');
        let chrom = cols.next();
        let pos = cols.nth(1);
        match (chrom, pos) {
            (Some(c), Some(p)) => {
                let pos: i64 = p.parse().map_err(|_| FlowError::BedParse {
                    line: i + 1,
                    msg: format!("bad pos '{p}'"),
                })?;
                sites.push((c.to_string(), pos));
            }
            _ => {
                return Err(FlowError::BedParse {
                    line: i + 1,
                    msg: "expected >= 3 tab-separated columns".into(),
                })
            }
        }
    }
    Ok(sites)
}

/// `--sites` file for BamSites entry: `chrom\tpos` (1-based).
fn read_sites_file(path: &Path) -> Result<Vec<(String, i64)>, FlowError> {
    let text = fs::read_to_string(path)?;
    let mut sites = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        let (chrom, pos) = line.split_once('\t').ok_or_else(|| FlowError::BedParse {
            line: i + 1,
            msg: "expected 'chrom<TAB>pos'".into(),
        })?;
        let pos: i64 = pos.parse().map_err(|_| FlowError::BedParse {
            line: i + 1,
            msg: format!("bad pos '{pos}'"),
        })?;
        sites.push((chrom.to_string(), pos));
    }
    Ok(sites)
}

fn stage_score(
    params: &RunParams,
    bam: &Path,
    sites: &[(String, i64)],
) -> Result<Vec<f64>, FlowError> {
    let dir = params.out_dir.join("score");
    fs::create_dir_all(&dir)?;
    if sites.is_empty() {
        // Empty candidates: empty scores.tsv + header-only VCF, exit 0.
        fs::write(dir.join("scores.tsv"), "")?;
        return Ok(Vec::new());
    }
    eprintln!("[score] running");
    let caduceus = match &params.caduceus {
        Some(c) => c.clone(),
        None => score_pipeline::resolve_encoder_from_bundle(&params.bundle)
            .map_err(anyhow_stage_err("score"))?,
    };
    let probs = score_pipeline::score_sites_batched(
        bam,
        &params.fasta,
        &caduceus,
        &params.bundle,
        sites,
        params.threads,
        params.batch.max(1),
        None,
    )
    .map_err(anyhow_stage_err("score"))?;

    let mut text = String::new();
    for ((chrom, pos), prob) in sites.iter().zip(&probs) {
        use std::fmt::Write as _;
        let _ = writeln!(text, "{chrom}\t{pos}\t{prob}");
    }
    fs::write(dir.join("scores.tsv"), text)?;
    eprintln!("[score] {} sites", probs.len());
    Ok(probs)
}
