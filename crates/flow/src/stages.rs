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
use rust_htslib::bam::Read as _;
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
    preflight(params, entry)?;
    fs::create_dir_all(&params.out_dir)?;
    run_from(params, entry, crate::resume::Stage::first_for(entry))
}

/// Entry + guardrail checks shared by fresh runs and resumes.
pub fn preflight(params: &RunParams, entry: Entry) -> Result<(), FlowError> {
    // Reference manifest copied into the run dir at run start (hybrid refs);
    // absent → legacy hg38 heuristic.
    let manifest = crate::manifest::SpeciesManifest::read(&params.out_dir);
    guard::check_species(&params.fasta, manifest.as_ref())?;
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
    Ok(())
}

/// Stage machine: execute from `start` onward, reloading upstream artifacts
/// instead of rebuilding them. `run` covers the full sequence; `resume`
/// enters at the first invalid stage (see `crate::resume::walk`).
pub fn run_from(
    params: &RunParams,
    entry: Entry,
    start: crate::resume::Stage,
) -> Result<(), FlowError> {
    use crate::resume::Stage;
    let map_dir = params.out_dir.join("map");
    let mut current_bam: PathBuf;
    let mut baln: Option<PathBuf> = None;
    match entry {
        Entry::FastqSe | Entry::FastqPe => {
            let (clean1, clean2) = if start <= Stage::Qc {
                stage_qc(params, entry)?
            } else {
                qc_outputs(params, entry)
            };
            if start <= Stage::Map {
                // map_stage runs the collapsed rescue before returning.
                current_bam = stage_map(params, entry, &clean1, clean2.as_deref())?;
            } else {
                current_bam = map_dir.join("raw.bam");
                if start == Stage::Rescue {
                    stage_rescue_collapsed(params, &current_bam)?;
                }
            }
            if start <= Stage::Sort {
                current_bam = stage_sort(&current_bam, params)?;
            } else {
                current_bam = map_dir.join("sorted.bam");
            }
            // The mapper writes .baln for SE and PE alike; use it as the
            // scan fast channel. A missing or stub (0-byte, pre-channel
            // run directories) .baln falls back to the BAM.
            let baln_path = map_dir.join("align.baln");
            baln = baln_path
                .metadata()
                .map(|m| m.len() > 12)
                .unwrap_or(false)
                .then_some(baln_path);
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
            let bed = params.out_dir.join("scan").join("candidates.bed");
            if start <= Stage::Scan {
                let bed = stage_scan(params, entry, &current_bam, baln)?;
                let fstats = crate::filter::CallFilter::default().apply_to_bed(&bed)?;
                eprintln!(
                    "[filter] {} candidates -> {} kept (low_depth {} no_signal {})",
                    fstats.input, fstats.kept, fstats.low_depth, fstats.no_signal
                );
                bed_to_sites(&bed)?
            } else {
                bed_to_sites(&bed)?
            }
        }
    };

    if start <= Stage::Vcf {
        let probs = if start <= Stage::Score {
            stage_score(params, &current_bam, &sites)?
        } else {
            load_scores(&params.out_dir.join("score").join("scores.tsv"))?
        };
        vcf::write_vcf(params, entry, &sites, &probs)?;
    }
    stage_report(params);
    Ok(())
}

/// scores.tsv → probs bridge for resumes that skip the score stage.
/// `NA` marks unscored rows (hybrid runs: mouse-contig sites).
fn load_scores(path: &Path) -> Result<Vec<Option<f64>>, FlowError> {
    let text = fs::read_to_string(path)?;
    let mut probs = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        let p = line
            .split('\t')
            .nth(2)
            .map(|v| {
                if v == "NA" {
                    Ok(None)
                } else {
                    v.parse::<f64>().map(Some).map_err(|_| ())
                }
            })
            .and_then(|r| r.ok())
            .ok_or_else(|| FlowError::BedParse {
                line: i + 1,
                msg: "scores.tsv expects 'chrom<TAB>pos<TAB>prob' (or NA)".into(),
            })?;
        probs.push(p);
    }
    Ok(probs)
}

/// Report stage (best-effort): pack the finished run directory into a
/// standalone `<out>/report.html`. A report failure must never lose the
/// science outputs — errors are printed and the pipeline still succeeds.
fn stage_report(params: &RunParams) {
    match params.gtf.as_deref() {
        Some(gtf) => match esperanto_report::generate(&params.out_dir, &params.fasta, gtf) {
            Ok(p) => eprintln!("[report] written: {}", p.display()),
            Err(e) => eprintln!("[report] warning: report generation failed: {e}"),
        },
        None => eprintln!("[report] warning: no GTF available; skipped the HTML report"),
    }
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
    Ok(qc_outputs(params, entry))
}

/// qc output paths without running qc (resume reload path).
fn qc_outputs(params: &RunParams, entry: Entry) -> (PathBuf, Option<PathBuf>) {
    let dir = params.out_dir.join("qc");
    let stem1 = stem_of(&params.r1[0]);
    match entry {
        Entry::FastqSe => (dir.join(format!("{stem1}.clean.fq.gz")), None),
        _ => {
            let stem2 = stem_of(&params.r2[0]);
            (
                dir.join(format!("{stem1}.clean_R1.fq.gz")),
                Some(dir.join(format!("{stem2}.clean_R2.fq.gz"))),
            )
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
    rescue_collapsed(paidx, Some(index.reference), dir, threads, &raw_bam)?;
    Ok(raw_bam)
}

// Collapsed-alphabet (A==G, T==C) rescue of unmapped reads. When a
// `<index stem>.cpaidx` sits next to the alignment index, the unmapped set
// is re-aligned against it. Every placement is then verified against the
// true (four-letter) reference: a hyperedited read must carry a dense,
// high-quality, editing-dominated mismatch cluster with sane cluster
// geometry and, for paired data, a concordant mapped mate (Porath,
// Carmi & Levanon, Nat Commun 5:4726, 2014). Verified survivors are
// written back into raw.bam with MAPQ 0 and an `RE:Z:collapsed` tag
// (repeat-family placement, never confident); everything else stays
// unmapped, and unmapped.fq.gz is rewritten without the verified reads
// only. No-op when no collapsed index exists.

/// Maximum non-editing-class mismatches tolerated in a verified
/// placement (collapsed-space edit distance).
const RESCUE_MAX_NON_EDIT_MM: usize = 2;
/// Reads at or below this length need a stricter editing purity.
const RESCUE_SHORT_READ_MAX_LEN: usize = 60;
/// Editing-class share of all mismatches required (long reads).
const RESCUE_PURITY_LONG: f64 = 0.60;
/// Editing-class share of all mismatches required (reads <= 60 bp).
const RESCUE_PURITY_SHORT: f64 = 0.80;
/// Minimum high-quality editing-class mismatch density over the read.
const RESCUE_MIN_EDIT_DENSITY: f64 = 0.05;
/// Phred floor for an editing-class mismatch to count toward the cluster.
const RESCUE_MIN_EDIT_PHRED: u8 = 30;
/// Mean Phred floor (after dropping the lowest decile) for a read to
/// enter the realignment at all.
const RESCUE_MIN_MEAN_PHRED: f64 = 25.0;
/// Per-nucleotide fraction bounds for a read to enter the realignment.
const RESCUE_MAX_NT_FRAC: f64 = 0.60;
const RESCUE_MIN_NT_FRAC: f64 = 0.10;
/// Ambiguous-base fraction bound for a read to enter the realignment.
const RESCUE_MAX_N_FRAC: f64 = 0.10;
/// Longest tolerated homopolymer run, in bases.
const RESCUE_MAX_HOMOPOLYMER: usize = 20;
/// Longest tolerated tandem dinucleotide repeat, in units.
const RESCUE_MAX_DINUC_UNITS: usize = 10;
/// Minimum first-to-last editing-mismatch span, as a fraction of read
/// length.
const RESCUE_MIN_CLUSTER_SPAN_FRAC: f64 = 0.10;
/// Clusters fully inside the first or last fraction of the read are
/// rejected (splice-junction artifacts).
const RESCUE_END_EXCLUSION_FRAC: f64 = 0.20;
/// Maximum single-nucleotide share inside the editing cluster.
const RESCUE_CLUSTER_MAX_NT_FRAC: f64 = 0.60;
/// Paired-end concordance window: the mate must already be mapped within
/// this distance on the same contig, in the opposite orientation.
const RESCUE_PE_MATE_WINDOW: i64 = 500_000;

fn stage_rescue_collapsed(params: &RunParams, raw_bam: &Path) -> Result<(), FlowError> {
    let Some(paidx) = params.index.clone() else {
        return Ok(());
    };
    rescue_collapsed(
        &paidx,
        None,
        &params.out_dir.join("map"),
        params.threads,
        raw_bam,
    )
}

fn rescue_collapsed(
    paidx: &Path,
    live_reference: Option<&'static esperanto_map::fasta::Reference>,
    map_dir: &Path,
    threads: usize,
    raw_bam: &Path,
) -> Result<(), FlowError> {
    let cpath = Some(paidx.with_extension("cpaidx")).filter(|p| p.is_file());
    let Some(cpath) = cpath else { return Ok(()) };
    let unm_path = map_dir.join("unmapped.fq.gz");
    if !unm_path.is_file() {
        return Ok(());
    }
    let cidx = index_io::load(&cpath).map_err(stage_err("map"))?;
    // RNA preset (dense seeding; the de-novo splice paths find no canonical
    // signal in collapsed space and fall through to contiguous placement);
    // seed k/w and the chain anchor reward must match the collapsed index.
    let mut cfg = esperanto_map::align::AlignConfig {
        seed: cidx.params,
        ..esperanto_map::align::AlignConfig::rna_default()
    };
    cfg.chain.k = cidx.params.k as i32;
    cfg.extend.editing_aware = true;

    // Read the unmapped set.
    let fq: Vec<(String, Vec<u8>, Vec<u8>)> = {
        let f = File::open(&unm_path)?;
        let mut dec = flate2::read::GzDecoder::new(f);
        let mut text = String::new();
        use std::io::Read as _;
        dec.read_to_string(&mut text)?;
        let mut out = Vec::new();
        let mut lines = text.lines();
        while let (Some(n), Some(q), Some(_), Some(qv)) =
            (lines.next(), lines.next(), lines.next(), lines.next())
        {
            out.push((
                n.trim_start_matches('@').to_string(),
                q.as_bytes().to_vec(),
                qv.as_bytes().to_vec(),
            ));
        }
        out
    };
    if fq.is_empty() {
        return Ok(());
    }

    // Artifact screen before any realignment work: composition and
    // quality gates discard reads that cannot carry a trustworthy editing
    // signal.
    let cand_idx: Vec<usize> = fq
        .iter()
        .enumerate()
        .filter(|(_, (_, seq, qual))| rescue_prefilter(seq, qual))
        .map(|(i, _)| i)
        .collect();
    let prefiltered = (fq.len() - cand_idx.len()) as u64;

    // Honest accounting into align_qc.json; the rescue runs after the
    // stats document is written, so the keys are patched in here (the
    // report stage treats absence as 0).
    let patch_stats = |verified: u64, rejected: u64| -> Result<(), FlowError> {
        let p = map_dir.join("align_qc.json");
        let mut stats: esperanto_map::stats::AlignStats =
            serde_json::from_str(&fs::read_to_string(&p)?).map_err(stage_err("map"))?;
        stats.rescued_collapsed = Some(verified);
        stats.rescue_rejected_collapsed = Some(rejected);
        stats.rescue_prefiltered_collapsed = Some(prefiltered);
        let json = serde_json::to_string_pretty(&stats).map_err(stage_err("map"))?;
        fs::write(&p, format!("{json}\n"))?;
        Ok(())
    };

    let collapse = |seq: &[u8]| -> Vec<u8> {
        seq.iter()
            .map(|&b| match b.to_ascii_uppercase() {
                b'A' | b'G' => b'G',
                b'T' | b'C' => b'C',
                _ => b'N',
            })
            .collect()
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| FlowError::Entry(format!("rescue pool: {e}")))?;
    // Two passes (novel-junction discovery): pass 1 collects spliced
    // placements, discoveries with support >= 2 form a junction library,
    // pass 2 re-aligns against it.
    let align_all = |lib: &Option<Arc<esperanto_map::gtf::JunctionLib>>| {
        pool.install(|| {
            use rayon::prelude::*;
            cand_idx
                .par_iter()
                .map_init(
                    || {
                        let mut a = esperanto_map::align::Aligner::new(&cidx, cfg);
                        a.jlib = lib.clone();
                        a
                    },
                    |al, &i| {
                        let c = collapse(&fq[i].1);
                        al.align_read(&c)
                    },
                )
                .collect::<Vec<Option<esperanto_map::mapq::ReadAlignment>>>()
        })
    };
    let pass1 = align_all(&None);
    let mut counts: std::collections::BTreeMap<
        (u32, u32, u32, bool),
        (esperanto_map::gtf::Junction, u32),
    > = std::collections::BTreeMap::new();
    for aln in pass1.iter().flatten() {
        for j in &aln.junctions {
            let key = (
                j.junction.contig,
                j.junction.start,
                j.junction.end,
                j.junction.minus_strand,
            );
            match counts.get_mut(&key) {
                Some((_, c)) => *c += 1,
                None => {
                    counts.insert(key, (j.junction, 1));
                }
            }
        }
    }
    let merged: Vec<(esperanto_map::gtf::Junction, u32)> =
        counts.into_values().filter(|(_, c)| *c >= 2).collect();
    let results = if merged.is_empty() {
        pass1
    } else {
        let lib = esperanto_map::gtf::JunctionLib::build_with_counts(merged);
        align_all(&Some(Arc::new(lib)))
    };

    let placed: Vec<(usize, esperanto_map::mapq::ReadAlignment)> = cand_idx
        .into_iter()
        .zip(results)
        .filter_map(|(i, aln)| aln.map(|a| (i, a)))
        .collect();
    let placed_count = placed.len();
    if placed.is_empty() {
        patch_stats(0, 0)?;
        eprintln!(
            "[rescue] collapsed: 0 of {} unmapped placed ({} artifact reads screened)",
            fq.len(),
            prefiltered
        );
        return Ok(());
    }

    // True (four-letter) reference for post-placement verification. The
    // map stage hands its live index reference over; standalone callers
    // load the paidx once here.
    let true_ref: &'static esperanto_map::fasta::Reference = match live_reference {
        Some(r) => r,
        None => index_io::load(paidx).map_err(stage_err("map"))?.reference,
    };

    // Paired-end concordance table: mate placements from the primary
    // alignment (before any merge), both template slots per name.
    let pe = fq
        .iter()
        .any(|(n, _, _)| n.ends_with("/1") || n.ends_with("/2"));
    let mates: Option<std::collections::HashMap<Vec<u8>, [Option<MateLoc>; 2]>> = if pe {
        Some(collect_mate_placements(raw_bam)?)
    } else {
        None
    };

    let mut rescued_names: std::collections::HashSet<&str> = std::collections::HashSet::new();
    let mut rescued: Vec<(String, Vec<u8>, Vec<u8>, esperanto_map::mapq::ReadAlignment)> =
        Vec::new();
    let mut rejected = 0u64;
    for (i, mut a) in placed {
        let (name, seq, qual) = &fq[i];
        if verify_rescue_placement(name, &a, seq, qual, true_ref, mates.as_ref()) {
            a.second_chain_score = a.chain_score; // MAPQ -> 0 (repeat-family placement)
            a.rescued = true;
            rescued.push((name.clone(), seq.clone(), qual.clone(), a));
            rescued_names.insert(name.as_str());
        } else {
            rejected += 1;
        }
    }
    patch_stats(rescued.len() as u64, rejected)?;
    if rescued.is_empty() {
        eprintln!(
            "[rescue] collapsed: 0 verified ({} placed, {} rejected, {} screened) of {} unmapped",
            placed_count, rejected, prefiltered, fq.len()
        );
        return Ok(());
    }

    // Placement sidecar for the report's hyperedited-region track: one
    // `chrom<TAB>pos` row (0-based) per rescued read, written before any
    // merge step so an interruption still leaves the raw data on disk.
    {
        let mut text = String::new();
        use std::fmt::Write as _;
        for (_, _, _, a) in &rescued {
            let _ = writeln!(
                text,
                "{}\t{}",
                cidx.reference.contigs[a.contig as usize].name, a.pos
            );
        }
        fs::write(map_dir.join("rescued.bed"), text)?;
    }

    // Rewrite unmapped.fq.gz without the rescued reads.
    {
        let mut out = Vec::new();
        for (name, seq, qual) in &fq {
            if !rescued_names.contains(name.as_str()) {
                out.extend_from_slice(format!("@{name}\n").as_bytes());
                out.extend_from_slice(seq);
                out.extend_from_slice(b"\n+\n");
                out.extend_from_slice(qual);
                out.push(b'\n');
            }
        }
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(1));
        use std::io::Write as _;
        enc.write_all(&out)?;
        fs::write(&unm_path, enc.finish()?)?;
    }

    // Merge rescued records into raw.bam (same contig order: both indices
    // derive from the same FASTA) -> replace raw.bam before sort.
    let header = {
        let rdr = rust_htslib::bam::Reader::from_path(raw_bam).map_err(stage_err("map"))?;
        rust_htslib::bam::Header::from_template(rdr.header())
    };
    let merged = raw_bam.with_extension("merged.bam");
    {
        let mut w = rust_htslib::bam::Writer::from_path(
            &merged,
            &header,
            rust_htslib::bam::Format::Bam,
        )
        .map_err(stage_err("map"))?;
        {
            let mut rdr = rust_htslib::bam::Reader::from_path(raw_bam).map_err(stage_err("map"))?;
            for r in rdr.records() {
                let rec = r.map_err(stage_err("map"))?;
                // Drop the rescued reads' original unmapped records (one
                // record per read name survives).
                if rec.is_unmapped()
                    && rescued_names.contains(String::from_utf8_lossy(rec.qname()).as_ref())
                {
                    continue;
                }
                w.write(&rec).map_err(stage_err("map"))?;
            }
        }
        let header_view = rust_htslib::bam::HeaderView::from_header(&header);
        for (name, seq, qual, a) in &rescued {
            let mut rec = rust_htslib::bam::Record::new();
            let qname = name.as_bytes();
            let rev = a.strand == esperanto_map::seed::Strand::Minus;
            let (s2, q2) = esperanto_bamio::apply_t13(rev, seq, qual);
            let tid = header_view
                .tid(cidx.reference.contigs[a.contig as usize].name.as_bytes())
                .ok_or_else(|| FlowError::Entry("rescued contig missing from BAM header".into()))?;
            rec.set(qname, None, &s2, &q2);
            rec.set_tid(tid as i32);
            rec.set_pos(a.pos as i64);
            rec.set_mapq(0);
            let mut flag = 0u16;
            if rev {
                flag |= 0x10;
            }
            rec.set_flags(flag);
            let cigar: Vec<rust_htslib::bam::record::Cigar> = a
                .cigar
                .iter()
                .map(|op| match op {
                    esperanto_map::extend::CigarOp::Match(n) => {
                        rust_htslib::bam::record::Cigar::Match(*n)
                    }
                    esperanto_map::extend::CigarOp::Ins(n) => {
                        rust_htslib::bam::record::Cigar::Ins(*n)
                    }
                    esperanto_map::extend::CigarOp::Del(n) => {
                        rust_htslib::bam::record::Cigar::Del(*n)
                    }
                    esperanto_map::extend::CigarOp::RefSkip(n) => {
                        rust_htslib::bam::record::Cigar::RefSkip(*n)
                    }
                    esperanto_map::extend::CigarOp::SoftClip(n) => {
                        rust_htslib::bam::record::Cigar::SoftClip(*n)
                    }
                })
                .collect();
            rec.set(qname, Some(&rust_htslib::bam::record::CigarString(cigar)), &s2, &q2);
            rec.push_aux(b"RE", rust_htslib::bam::record::Aux::String("collapsed"))
                .map_err(stage_err("map"))?;
            w.write(&rec).map_err(stage_err("map"))?;
        }
    }
    // The SE fast channel (.baln) must carry the rescued records too —
    // baln has no trailer, records append cleanly.
    let baln_path = map_dir.join("align.baln");
    if baln_path.is_file() {
        let mut w = std::io::BufWriter::new(
            std::fs::OpenOptions::new()
                .append(true)
                .open(&baln_path)?,
        );
        for (name, seq, qual, a) in &rescued {
            let mut a2 = a.clone();
            a2.rescued = false; // the RE tag below carries the provenance
            let mut br = esperanto_map::bam::record_se(name, seq, qual, Some(a2));
            br.mapq = 0;
            if let Some(view) = br.aln.as_mut() {
                view.tags.push(esperanto_bamio::RawTag(
                    *b"RE",
                    esperanto_bamio::TagValue::Str("collapsed".into()),
                ));
            }
            esperanto_map::baln::write_record(&mut w, &br).map_err(stage_err("map"))?;
        }
        use std::io::Write as _;
        w.flush()?;
    }
    fs::rename(&merged, raw_bam)?;
    eprintln!(
        "[rescue] collapsed: {} verified of {} placed ({} rejected, {} screened) / {} unmapped (MAPQ 0)",
        rescued.len(),
        placed_count,
        rejected,
        prefiltered,
        fq.len()
    );
    Ok(())
}

/// Mapped-mate placement slot for the paired-end concordance check.
struct MateLoc {
    tid: i32,
    pos: i64,
    reverse: bool,
}

/// Build the name -> both-template-slots table of mapped placements from
/// the primary alignment output.
fn collect_mate_placements(
    raw_bam: &Path,
) -> Result<std::collections::HashMap<Vec<u8>, [Option<MateLoc>; 2]>, FlowError> {
    let mut rdr = rust_htslib::bam::Reader::from_path(raw_bam).map_err(stage_err("map"))?;
    let mut map: std::collections::HashMap<Vec<u8>, [Option<MateLoc>; 2]> =
        std::collections::HashMap::new();
    for r in rdr.records() {
        let rec = r.map_err(stage_err("map"))?;
        if rec.is_unmapped() {
            continue;
        }
        let slot = if rec.is_first_in_template() { 0 } else { 1 };
        let loc = MateLoc {
            tid: rec.tid(),
            pos: rec.pos(),
            reverse: rec.is_reverse(),
        };
        map.entry(rec.qname().to_vec()).or_insert([None, None])[slot] = Some(loc);
    }
    Ok(map)
}

/// Strip the `/1`/`/2` mate marker and report the template slot (0 = first
/// in template). Names without a marker take slot 0.
fn split_mate_marker(name: &str) -> (&str, usize) {
    if let Some(base) = name.strip_suffix("/1") {
        (base, 0)
    } else if let Some(base) = name.strip_suffix("/2") {
        (base, 1)
    } else {
        (name, 0)
    }
}

/// Artifact screen for the collapsed realignment: composition and quality
/// gates applied before any alignment work. Returns `true` when the read
/// may enter the realignment.
fn rescue_prefilter(seq: &[u8], qual: &[u8]) -> bool {
    let len = seq.len();
    if len == 0 || qual.len() != len {
        return false;
    }
    let mut counts = [0usize; 5]; // A C G T other
    let mut run = 1usize;
    let mut max_homo = 1usize;
    let mut prev: Option<u8> = None;
    for &b in seq {
        let u = b.to_ascii_uppercase();
        match u {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' => counts[3] += 1,
            _ => counts[4] += 1,
        }
        run = if prev == Some(u) { run + 1 } else { 1 };
        max_homo = max_homo.max(run);
        prev = Some(u);
    }
    let acgt = counts[0] + counts[1] + counts[2] + counts[3];
    if (len - acgt) as f64 / len as f64 > RESCUE_MAX_N_FRAC {
        return false;
    }
    let max_nt = counts[0..4].iter().copied().max().unwrap_or(0) as f64 / len as f64;
    if !(RESCUE_MIN_NT_FRAC..=RESCUE_MAX_NT_FRAC).contains(&max_nt) {
        return false;
    }
    if max_homo > RESCUE_MAX_HOMOPOLYMER {
        return false;
    }
    // Tandem dinucleotide runs: the longest stretch with period 2 counts
    // as (stretch bases) / 2 units.
    let mut best = 0usize;
    let mut cur = 0usize;
    for i in 0..len.saturating_sub(2) {
        if seq[i].eq_ignore_ascii_case(&seq[i + 2]) {
            cur += 1;
            best = best.max(cur);
        } else {
            cur = 0;
        }
    }
    if (best + 1).div_ceil(2) > RESCUE_MAX_DINUC_UNITS {
        return false;
    }
    // Mean quality after dropping the lowest decile.
    let mut qs = qual.to_vec();
    qs.sort_unstable();
    let drop = len / 10;
    let kept = &qs[drop..];
    let mean = if kept.is_empty() {
        qs.iter().map(|&q| (q.saturating_sub(33)) as f64).sum::<f64>() / len as f64
    } else {
        kept.iter().map(|&q| (q.saturating_sub(33)) as f64).sum::<f64>() / kept.len() as f64
    };
    mean >= RESCUE_MIN_MEAN_PHRED
}

/// Post-placement verification against the true reference. A placement
/// survives only when its mismatch pattern is that of a hyperedited read:
/// a dense cluster of high-quality A-to-G / T-to-C differences dominating
/// every other mismatch class, sane cluster geometry, and — for paired
/// data — a concordant already-mapped mate.
fn verify_rescue_placement(
    name: &str,
    aln: &esperanto_map::mapq::ReadAlignment,
    seq: &[u8],
    qual: &[u8],
    true_ref: &esperanto_map::fasta::Reference,
    mates: Option<&std::collections::HashMap<Vec<u8>, [Option<MateLoc>; 2]>>,
) -> bool {
    use esperanto_map::extend::CigarOp;
    use esperanto_map::fasta::Base;

    let read_len = seq.len();
    if read_len == 0 || qual.len() != read_len {
        return false;
    }
    let Some(contig) = true_ref.contigs.get(aln.contig as usize) else {
        return false;
    };
    // The collapsed realignment is an ungapped placement channel: any
    // insertion or deletion disqualifies the placement outright.
    let mut ref_span = 0u32;
    for op in &aln.cigar {
        match op {
            CigarOp::Ins(_) | CigarOp::Del(_) => return false,
            CigarOp::Match(n) | CigarOp::RefSkip(n) => ref_span = ref_span.saturating_add(*n),
            CigarOp::SoftClip(_) => {}
        }
    }
    if aln.pos as u64 + ref_span as u64 > contig.len as u64 {
        return false;
    }

    // Read bases and qualities in reference orientation.
    let minus = aln.strand == esperanto_map::seed::Strand::Minus;
    let (rseq, rqual) = esperanto_bamio::apply_t13(minus, seq, qual);

    let mut non_edit = 0usize;
    let mut total_mm = 0usize;
    let mut edit_pos: Vec<usize> = Vec::new();
    let mut rpos = aln.pos;
    let mut qpos = 0usize;
    for op in &aln.cigar {
        match op {
            CigarOp::Match(n) => {
                for _ in 0..*n {
                    if qpos >= read_len {
                        return false; // malformed CIGAR: consumes more read than exists
                    }
                    let rb = Base::from_ascii(rseq[qpos]);
                    let fb = contig.base(rpos);
                    if rb == Base::N || fb == Base::N {
                        total_mm += 1;
                        non_edit += 1;
                    } else if rb != fb {
                        total_mm += 1;
                        let editing = (fb == Base::A && rb == Base::G)
                            || (fb == Base::T && rb == Base::C);
                        if editing
                            && rqual[qpos].saturating_sub(33) >= RESCUE_MIN_EDIT_PHRED
                        {
                            edit_pos.push(qpos);
                        }
                    }
                    rpos += 1;
                    qpos += 1;
                }
            }
            CigarOp::RefSkip(n) => rpos += *n,
            CigarOp::SoftClip(n) => qpos += *n as usize,
            CigarOp::Ins(_) | CigarOp::Del(_) => {
                unreachable!("gapped placements are rejected before the walk")
            }
        }
    }

    if non_edit > RESCUE_MAX_NON_EDIT_MM {
        return false;
    }
    let ne = edit_pos.len();
    if ne == 0 {
        return false;
    }
    let nlen = read_len as f64;
    if (ne as f64) < RESCUE_MIN_EDIT_DENSITY * nlen {
        return false;
    }
    let purity = if read_len <= RESCUE_SHORT_READ_MAX_LEN {
        RESCUE_PURITY_SHORT
    } else {
        RESCUE_PURITY_LONG
    };
    if ne as f64 <= purity * total_mm as f64 {
        return false;
    }
    let first = edit_pos[0];
    let last = edit_pos[ne - 1];
    if ((last - first) as f64) < RESCUE_MIN_CLUSTER_SPAN_FRAC * nlen {
        return false;
    }
    if (last as f64) < RESCUE_END_EXCLUSION_FRAC * nlen {
        return false;
    }
    if (first as f64) > (1.0 - RESCUE_END_EXCLUSION_FRAC) * nlen {
        return false;
    }
    // Cluster composition: a single nucleotide dominating the cluster
    // region flags alignment slippage rather than editing.
    let mut cc = [0usize; 4];
    let window = last - first + 1;
    for &b in &rseq[first..=last] {
        match b.to_ascii_uppercase() {
            b'A' => cc[0] += 1,
            b'C' => cc[1] += 1,
            b'G' => cc[2] += 1,
            b'T' => cc[3] += 1,
            _ => {}
        }
    }
    if cc
        .iter()
        .any(|&c| c as f64 > RESCUE_CLUSTER_MAX_NT_FRAC * window as f64)
    {
        return false;
    }

    // Paired-end concordance: the mate must be mapped on the same contig
    // within the insert window, in the opposite orientation.
    if let Some(mates) = mates {
        let (base, slot) = split_mate_marker(name);
        let Some(mate) = mates
            .get(base.as_bytes())
            .and_then(|slots| slots[1 - slot].as_ref())
        else {
            return false;
        };
        if mate.tid != aln.contig as i32 {
            return false;
        }
        if (mate.pos - aln.pos as i64).abs() > RESCUE_PE_MATE_WINDOW {
            return false;
        }
        if mate.reverse == minus {
            return false;
        }
    }
    true
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
) -> Result<Vec<Option<f64>>, FlowError> {
    let dir = params.out_dir.join("score");
    fs::create_dir_all(&dir)?;
    if sites.is_empty() {
        // Empty candidates: empty scores.tsv + header-only VCF, exit 0.
        fs::write(dir.join("scores.tsv"), "")?;
        return Ok(Vec::new());
    }
    eprintln!("[score] running");
    // Hybrid runs: human-locus sites score with the bundle; mouse-contig
    // sites are never scored with the human model (scientific contract) —
    // they reach the VCF marked UNSCORED with VAF/DEPTH kept.
    let manifest = crate::manifest::SpeciesManifest::read(&params.out_dir);
    let hybrid = manifest.as_ref().is_some_and(|m| m.kind == "hybrid");
    let caduceus = match &params.caduceus {
        Some(c) => c.clone(),
        None => score_pipeline::resolve_encoder_from_bundle(&params.bundle)
            .map_err(anyhow_stage_err("score"))?,
    };
    let ask: Option<&(dyn Fn() -> bool + Send + Sync)> = params
        .device_ask
        .as_ref()
        .map(|a| a.fn_ref());
    let scored: Vec<f64> = if hybrid {
        let m = manifest.as_ref().expect("hybrid checked");
        let human_sites: Vec<(String, i64)> = sites
            .iter()
            .filter(|(c, _)| m.owner_of(c) == Some("human"))
            .cloned()
            .collect();
        let mouse_sites: Vec<(String, i64)> = sites
            .iter()
            .filter(|(c, _)| m.owner_of(c) == Some("mouse"))
            .cloned()
            .collect();
        let mut by_site: std::collections::HashMap<(String, i64), f64> =
            std::collections::HashMap::with_capacity(sites.len());
        if !human_sites.is_empty() {
            let probs = score_pipeline::score_sites_batched(
                bam,
                &params.fasta,
                &caduceus,
                &params.bundle,
                &human_sites,
                params.threads,
                params.batch.max(1),
                None,
                params.device,
                ask,
                None,
                score_pipeline::ReferenceCheck::TrustedHybrid,
            )
            .map_err(anyhow_stage_err("score"))?;
            for (s, p) in human_sites.iter().zip(probs) {
                by_site.insert(s.clone(), p);
            }
        }
        if !mouse_sites.is_empty() {
            // Mouse-contig sites score with the mouse bundle (manifest
            // frozen at build time); without one they stay UNSCORED.
            let mouse_bundle: Option<std::path::PathBuf> = m.bundles.get("mouse").cloned();
            match mouse_bundle {
                Some(mb) => {
                    let m_caduceus = match &params.caduceus {
                        // the encoder is shared (same frozen backbone)
                        Some(c) => c.clone(),
                        None => score_pipeline::resolve_encoder_from_bundle(&mb)
                            .map_err(anyhow_stage_err("score"))?,
                    };
                    let probs = score_pipeline::score_sites_batched(
                        bam,
                        &params.fasta,
                        &m_caduceus,
                        &mb,
                        &mouse_sites,
                        params.threads,
                        params.batch.max(1),
                        None,
                        params.device,
                        ask,
                        None,
                        score_pipeline::ReferenceCheck::TrustedHybrid,
                    )
                    .map_err(anyhow_stage_err("score"))?;
                    for (s, p) in mouse_sites.iter().zip(probs) {
                        by_site.insert(s.clone(), p);
                    }
                    eprintln!(
                        "[score] hybrid: {} human + {} mouse sites scored (mouse model)",
                        human_sites.len(),
                        mouse_sites.len()
                    );
                }
                None => {
                    eprintln!(
                        "[score] hybrid: {} mouse-contig sites marked UNSCORED (no mouse model installed)",
                        mouse_sites.len()
                    );
                }
            }
        }
        sites
            .iter()
            .map(|s| by_site.get(s).copied().unwrap_or(f64::NAN))
            .collect::<Vec<f64>>()
    } else {
        score_pipeline::score_sites_batched(
            bam,
            &params.fasta,
            &caduceus,
            &params.bundle,
            sites,
            params.threads,
            params.batch.max(1),
            None,
            params.device,
            ask,
            None,
            score_pipeline::ReferenceCheck::Guardrail,
        )
        .map_err(anyhow_stage_err("score"))?
    };

    // Reassemble in site order (None = unscored rows: mouse-contig sites
    // without a mouse model surface as NAN in the hybrid branch).
    let mut probs: Vec<Option<f64>> = Vec::with_capacity(sites.len());
    if hybrid {
        for v in scored {
            probs.push(if v.is_nan() { None } else { Some(v) });
        }
    } else {
        for v in scored {
            probs.push(Some(v));
        }
    }
    let mut text = String::new();
    for ((chrom, pos), prob) in sites.iter().zip(&probs) {
        use std::fmt::Write as _;
        match prob {
            Some(p) => {
                let _ = writeln!(text, "{chrom}\t{pos}\t{p}");
            }
            None => {
                let _ = writeln!(text, "{chrom}\t{pos}\tNA");
            }
        }
    }
    fs::write(dir.join("scores.tsv"), text)?;
    eprintln!("[score] {} sites", probs.len());
    Ok(probs)
}


#[cfg(test)]
mod rescue_tests {
    use super::*;

    /// Deterministic reference: an `A` every `period` bases inside a mixed
    /// background (all four bases present, no long runs).
    fn reference_with_period(period: usize) -> &'static esperanto_map::fasta::Reference {
        let bg = b"CGTATGCCTAGG";
        let mut seq = Vec::with_capacity(400);
        for i in 0..400 {
            seq.push(if i % period == 0 {
                b'A'
            } else {
                bg[i % bg.len()]
            });
        }
        let text = format!(">c1\n{}\n", String::from_utf8(seq).unwrap());
        let r = esperanto_map::fasta::parse_fasta_bytes(text.as_bytes()).unwrap();
        Box::leak(Box::new(r))
    }

    fn alignment_at(pos: u32, n: u32) -> esperanto_map::mapq::ReadAlignment {
        esperanto_map::mapq::ReadAlignment {
            contig: 0,
            pos,
            strand: esperanto_map::seed::Strand::Plus,
            cigar: vec![esperanto_map::extend::CigarOp::Match(n)],
            ..esperanto_map::mapq::ReadAlignment::default()
        }
    }

    /// 150 bp read over `start` with every reference `A` inside read
    /// window [lo, hi] converted to `G`.
    fn hyper_read(r: &esperanto_map::fasta::Reference, start: usize, lo: usize, hi: usize) -> Vec<u8> {
        let mut seq = r.contigs[0].slice_ascii(start as u32, (start + 150) as u32);
        for (qi, b) in seq.iter_mut().enumerate() {
            if (lo..=hi).contains(&qi) && (start + qi).is_multiple_of(7) && *b == b'A' {
                *b = b'G';
            }
        }
        seq
    }

    #[test]
    fn split_mate_marker_slots() {
        assert_eq!(split_mate_marker("read/1"), ("read", 0));
        assert_eq!(split_mate_marker("read/2"), ("read", 1));
        assert_eq!(split_mate_marker("plain"), ("plain", 0));
    }

    #[test]
    fn prefilter_accepts_clean_read() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGT".repeat(6);
        let seq = &seq[..150];
        let qual = vec![b'I'; 150];
        assert!(rescue_prefilter(seq, &qual));
    }

    #[test]
    fn prefilter_rejects_artifacts() {
        let qual = vec![b'I'; 150];
        // Homopolymer run (25 A's, overall composition still balanced).
        let mut homo = vec![b'A'; 25];
        homo.extend_from_slice(&b"CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".repeat(3)[..125]);
        assert!(!rescue_prefilter(&homo, &qual));
        // Ambiguous bases: 20 Ns of 150.
        let mut nheavy = vec![b'N'; 20];
        nheavy.extend_from_slice(&b"ACGTACGT".repeat(17)[..130]);
        assert!(!rescue_prefilter(&nheavy, &qual));
        // Tandem dinucleotide repeat (AT x75).
        let dinuc = b"AT".repeat(75);
        assert!(!rescue_prefilter(&dinuc, &qual));
        // Biased composition: 62% one nucleotide.
        let mut biased = vec![b'A'; 93];
        biased.extend_from_slice(&b"CGTGCACTGGTCACGTGACTGGTCA".repeat(3)[..57]);
        assert!(!rescue_prefilter(&biased, &qual));
        // Uniformly low quality: trimmed mean stays below 25.
        let seq = b"ACGTACGTACGTACGTACGTACGTACGT".repeat(6);
        let lowq = vec![b'5'; 150]; // Phred 20
        assert!(!rescue_prefilter(&seq[..150], &lowq));
    }

    #[test]
    fn verify_accepts_hyperedited_read() {
        let r = reference_with_period(7);
        let seq = hyper_read(r, 20, 20, 110);
        assert!(seq.iter().filter(|&b| *b == b'G').count() >= 8);
        let aln = alignment_at(20, 150);
        let qual = vec![b'I'; 150];
        assert!(verify_rescue_placement("r", &aln, &seq, &qual, r, None));
    }

    #[test]
    fn verify_accepts_hyperedited_minus_strand_read() {
        let r = reference_with_period(7);
        let oriented = hyper_read(r, 20, 20, 110);
        // FASTQ stores the sequencer orientation: reverse-complement the
        // reference-forward view.
        let (seq, _) = esperanto_bamio::apply_t13(true, &oriented, &[b'I'; 150]);
        let mut aln = alignment_at(20, 150);
        aln.strand = esperanto_map::seed::Strand::Minus;
        let qual = vec![b'I'; 150];
        assert!(verify_rescue_placement("r", &aln, &seq, &qual, r, None));
    }

    #[test]
    fn verify_rejects_random_mismatches() {
        let r = reference_with_period(7);
        let mut seq = r.contigs[0].slice_ascii(20, 170);
        for qi in (5..145).step_by(4) {
            seq[qi] = match seq[qi] {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                _ => b'C',
            };
        }
        let aln = alignment_at(20, 150);
        let qual = vec![b'I'; 150];
        assert!(!verify_rescue_placement("r", &aln, &seq, &qual, r, None));
    }

    #[test]
    fn verify_rejects_gapped_placement() {
        let r = reference_with_period(7);
        let seq = hyper_read(r, 20, 20, 110);
        let mut aln = alignment_at(20, 0);
        aln.cigar = vec![
            esperanto_map::extend::CigarOp::Match(70),
            esperanto_map::extend::CigarOp::Del(1),
            esperanto_map::extend::CigarOp::Match(79),
        ];
        let qual = vec![b'I'; 150];
        assert!(!verify_rescue_placement("r", &aln, &seq, &qual, r, None));
    }

    #[test]
    fn verify_rejects_end_cluster() {
        // Dense period-3 editing packed into the first 20% of the read:
        // only the end-exclusion rule fails.
        let r = reference_with_period(3);
        let seq = hyper_read(r, 20, 0, 29);
        assert!(seq.iter().filter(|&b| *b == b'G').count() >= 8);
        let aln = alignment_at(20, 150);
        let qual = vec![b'I'; 150];
        assert!(!verify_rescue_placement("r", &aln, &seq, &qual, r, None));
    }

    #[test]
    fn verify_short_read_purity_threshold() {
        let r = reference_with_period(7);
        // 60 bp read: 4 editing mismatches, 0 non-editing -> accepted.
        let mut good = r.contigs[0].slice_ascii(20, 80);
        for qi in [8usize, 15, 22, 29] {
            if good[qi] == b'A' {
                good[qi] = b'G';
            }
        }
        let aln60 = alignment_at(20, 60);
        let qual60 = vec![b'I'; 60];
        assert!(verify_rescue_placement("r", &aln60, &good, &qual60, r, None));
        // 3 editing + 2 non-editing mismatches: purity 3/5 fails the 80%
        // short-read bar.
        let mut poor = good.clone();
        poor[5] = if poor[5] == b'A' { b'G' } else { b'A' }; // non-editing class
        poor[35] = if poor[35] == b'A' { b'G' } else { b'A' };
        // keep exactly three edits
        poor[29] = if poor[29] == b'G' { b'A' } else { poor[29] };
        assert!(!verify_rescue_placement("r", &aln60, &poor, &qual60, r, None));
    }

    #[test]
    fn verify_pe_mate_concordance() {
        let r = reference_with_period(7);
        let seq = hyper_read(r, 20, 20, 110);
        let aln = alignment_at(20, 150);
        let qual = vec![b'I'; 150];
        let mut mates = std::collections::HashMap::new();
        mates.insert(
            b"t".to_vec(),
            [
                Some(MateLoc {
                    tid: 0,
                    pos: 300,
                    reverse: true,
                }),
                None,
            ],
        );
        assert!(verify_rescue_placement("t/2", &aln, &seq, &qual, r, Some(&mates)));
        // Same orientation: rejected.
        let mut same = std::collections::HashMap::new();
        same.insert(
            b"t".to_vec(),
            [
                Some(MateLoc {
                    tid: 0, pos: 300, reverse: false }),
                None,
            ],
        );
        assert!(!verify_rescue_placement("t/2", &aln, &seq, &qual, r, Some(&same)));
        // Mate beyond the concordance window: rejected.
        let mut far = std::collections::HashMap::new();
        far.insert(
            b"t".to_vec(),
            [
                Some(MateLoc {
                    tid: 0, pos: 600_000, reverse: true }),
                None,
            ],
        );
        assert!(!verify_rescue_placement("t/2", &aln, &seq, &qual, r, Some(&far)));
        // Mate never mapped: rejected.
        let empty = std::collections::HashMap::new();
        assert!(!verify_rescue_placement("t/2", &aln, &seq, &qual, r, Some(&empty)));
    }
}