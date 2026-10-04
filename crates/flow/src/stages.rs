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
use esperanto_qc::QcParams;
use esperanto_score::load_primary_threshold;
use esperanto_score::pipeline as score_pipeline;
use rust_htslib::bam::Read as _;

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
    let stop_at: Option<Stage> = params
        .stop_after
        .as_deref()
        .map(|s| {
            Stage::parse(s).ok_or_else(|| {
                FlowError::Entry(format!(
                    "bad --stop-after '{s}': expected qc|map|sort|scan|score|vcf"
                ))
            })
        })
        .transpose()?;
    if let Some(st) = stop_at {
        if start > st {
            return Ok(());
        }
    }
    let map_dir = params.out_dir.join("map");
    let mut current_bam: PathBuf;
    match entry {
        Entry::FastqSe | Entry::FastqPe => {
            let (clean1, clean2) = if start <= Stage::Qc {
                stage_qc(params, entry)?
            } else {
                qc_outputs(params, entry)
            };
            if stop_at == Some(Stage::Qc) {
                eprintln!("[run] stop-after qc");
                return Ok(());
            }
            if start <= Stage::Map {
                // map_stage runs the collapsed rescue before returning.
                let (bam, _index, _l1) = stage_map(params, entry, &clean1, clean2.as_deref())?;
                current_bam = bam;
            } else {
                current_bam = map_dir.join("raw.bam");
                if start == Stage::Rescue {
                    stage_rescue_collapsed(params, &current_bam)?;
                }
            }
            if stop_at == Some(Stage::Map) {
                eprintln!("[run] stop-after map");
                return Ok(());
            }
            if start <= Stage::Sort {
                current_bam = stage_sort(&current_bam, params)?;
            } else {
                current_bam = map_dir.join("sorted.bam");
            }
            if stop_at == Some(Stage::Sort) {
                eprintln!("[run] stop-after sort");
                return Ok(());
            }
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
            let sites = if start <= Stage::Scan {
                let bed = stage_scan(params, &current_bam)?;
                bed_to_sites(&bed)?
            } else {
                bed_to_sites(&bed)?
            };
            if stop_at == Some(Stage::Scan) {
                eprintln!("[run] stop-after scan");
                return Ok(());
            }
            sites
        }
    };

    if start <= Stage::Vcf {
        let (scores, threshold) = if start <= Stage::Score {
            let scores = stage_score(params, &current_bam, &sites)?;
            if stop_at == Some(Stage::Score) {
                eprintln!("[run] stop-after score");
                return Ok(());
            }
            scores
        } else {
            load_score_state(params)?
        };
        vcf::write_vcf(params, entry, &sites, &scores, threshold)?;
    }
    if stop_at == Some(Stage::Vcf) {
        eprintln!("[run] stop-after vcf");
        return Ok(());
    }
    stage_report(params);
    Ok(())
}

/// scores.tsv → scores bridge for resumes that skip the score stage.
fn load_scores(path: &Path) -> Result<Vec<f64>, FlowError> {
    let text = fs::read_to_string(path)?;
    let mut scores = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        let p = line
            .split('\t')
            .nth(2)
            .map(|v| v.parse::<f64>())
            .and_then(|r| r.ok())
            .ok_or_else(|| FlowError::BedParse {
                line: i + 1,
                msg: "scores.tsv expects 'chrom<TAB>pos<TAB>score'".into(),
            })?;
        scores.push(p);
    }
    Ok(scores)
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

/// (raw BAM, leaked shared alignment index, L1 bundle when present). The
/// index handles remain process-lifetime for the alignment stage.
fn stage_map(
    params: &RunParams,
    entry: Entry,
    clean1: &Path,
    clean2: Option<&Path>,
) -> Result<
    (
        PathBuf,
        &'static esperanto_map::index::Index,
        Option<Arc<esperanto_engine::L1Index>>,
    ),
    FlowError,
> {
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
/// `out_dir`. Returns the raw (unsorted) BAM path, the process-lifetime
/// alignment index (leaked once so downstream stages share the mapping),
/// and the L1 bundle when present.
#[allow(clippy::too_many_arguments)]
pub fn map_stage(
    paidx: &Path,
    gtf_path: Option<&Path>,
    l1_bundle: Option<&Path>,
    r1: &Path,
    r2: Option<&Path>,
    out_dir: &Path,
    threads: usize,
) -> Result<
    (
        PathBuf,
        &'static esperanto_map::index::Index,
        Option<Arc<L1Index>>,
    ),
    FlowError,
> {
    let dir = out_dir;
    fs::create_dir_all(dir)?;
    let index: &'static esperanto_map::index::Index =
        Box::leak(Box::new(index_io::load(paidx).map_err(stage_err("map"))?));
    // Primary alignment uses all four bases. Hyperediting rescue uses
    // separate three-letter transforms (Porath, Carmi & Levanon,
    // Nat Commun 5:4726, 2014).
    let config = AlignConfig::rna_default();
    let jlib = match gtf_path {
        Some(g) => Some(Arc::new(
            gtf::from_gtf(g, |name| contig_id(index, name)).map_err(stage_err("map"))?,
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
        index,
        config,
        jlib,
        l1: l1.clone(),
    };
    let n_annotated = out.jlib.as_ref().map_or(0, |l| l.junctions.len());
    let (stats, discoveries, _merged, n_annotated_replaced) = match r2 {
        None => run_se_2pass(&mut out, r1, threads),
        Some(c2) => run_pe_2pass(&mut out, r1, c2, threads),
    }
    .map_err(stage_err("map"))?;

    let json = serde_json::to_string_pretty(&stats).map_err(stage_err("map"))?;
    fs::write(dir.join("align_qc.json"), format!("{json}\n"))?;
    eprintln!(
        "[map] {} / {} mapped",
        stats.mapped_reads, stats.total_reads
    );
    let n_merged_disc = discoveries
        .iter()
        .filter(|(_, _, f)| *f != esperanto_map::pipeline::DiscoveryFate::OutvotedByAnnotated)
        .count();
    let n_outvoted = discoveries
        .iter()
        .filter(|(_, _, f)| *f == esperanto_map::pipeline::DiscoveryFate::OutvotedByAnnotated)
        .count();
    let n_tie = discoveries
        .iter()
        .filter(|(_, _, f)| *f == esperanto_map::pipeline::DiscoveryFate::TieKept)
        .count();
    eprintln!(
        "[map] junction library: {n_annotated} annotated + {n_merged_disc} discovered merged ({n_outvoted} outvoted by annotated support, {n_annotated_replaced} annotated replaced, {n_tie} tie-kept)"
    );
    {
        let mut bed = String::new();
        for (j, support, fate) in &discoveries {
            let name = &index.reference.contigs[j.contig as usize].name;
            let strand = if j.minus_strand { '-' } else { '+' };
            let status = match fate {
                esperanto_map::pipeline::DiscoveryFate::Kept => "kept",
                esperanto_map::pipeline::DiscoveryFate::TieKept => "tie_kept",
                esperanto_map::pipeline::DiscoveryFate::OutvotesAnnotated => "outvotes_annotated",
                esperanto_map::pipeline::DiscoveryFate::OutvotedByAnnotated => "outvoted",
            };
            bed.push_str(&format!(
                "{name}\t{}\t{}\t{strand}\t{support}\t{status}\n",
                j.start, j.end
            ));
        }
        fs::write(dir.join("discovered_junctions.bed"), bed)?;
    }
    rescue_collapsed(paidx, Some(index.reference), dir, threads, &raw_bam)?;
    Ok((raw_bam, index, l1))
}

// Three-letter hyperediting rescue follows Porath, Carmi & Levanon,
// Nat Commun 5:4726, 2014. Each pass collapses one editing direction.
// Placements are checked against the four-letter reference; reads that
// satisfy both editing directions are rejected as ambiguous. Accepted reads
// receive RE:Z:collapsed and the mapping-quality cap in mapq_of.
// Rescue is skipped when neither transformed index exists.
/// Rescue passes: (editing direction, index extension, read from, read to).
/// Both indexes are required for each direction because canonical seeding
/// obeys revcomp(t_AG(revcomp(F))) == t_TC(F).
fn rescue_dirs() -> [((u8, u8), &'static str, u8, u8); 4] {
    [
        ((b'A', b'G'), "ag.paidx", b'A', b'G'),
        ((b'A', b'G'), "tc.paidx", b'A', b'G'),
        ((b'T', b'C'), "tc.paidx", b'T', b'C'),
        ((b'T', b'C'), "ag.paidx", b'T', b'C'),
    ]
}

/// Maximum mismatches outside the pass's own editing direction tolerated in a
/// verified placement — the transformed-space edit distance, since under the
/// matching transform that direction's edits cost nothing.
const RESCUE_MAX_NON_EDIT_MM: usize = esperanto_map::mapq::HYPEREDIT_MAX_NON_EDIT_MM as usize;
/// Gap penalty for ungapped collapsed realignment, as in Porath, Carmi &
/// Levanon (Nat Commun 5:4726, 2014). Local extension can soft-clip instead
/// of paying this penalty. The value stays within the i16 accumulator's
/// saturating arithmetic; accepted placements must also pass the gap check.
const RESCUE_GAP_PENALTY: i32 = 400;
/// Reads at or below this length need a stricter editing purity.
const RESCUE_SHORT_READ_MAX_LEN: usize = 60;
/// Editing-class share of all mismatches required (long reads).
const RESCUE_PURITY_LONG: f64 = 0.60;
/// Editing-class share of all mismatches required (reads <= 60 bp).
const RESCUE_PURITY_SHORT: f64 = 0.80;
/// Minimum high-quality editing-class mismatch density over the read.
const RESCUE_MIN_EDIT_DENSITY: f64 = esperanto_map::mapq::HYPEREDIT_MIN_DENSITY;
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
    if let Some(paidx) = &params.index {
        rescue_collapsed(
            paidx,
            None,
            &params.out_dir.join("map"),
            params.threads,
            raw_bam,
        )?;
    }
    Ok(())
}

fn rescue_collapsed(
    paidx: &Path,
    live_reference: Option<&'static esperanto_map::fasta::Reference>,
    map_dir: &Path,
    threads: usize,
    raw_bam: &Path,
) -> Result<u64, FlowError> {
    let mut dirs: Vec<((u8, u8), PathBuf, u8, u8)> = Vec::new();
    for &(dir, ext, from, to) in rescue_dirs().iter() {
        let path = paidx.with_extension(ext);
        if !path.is_file() {
            return Err(FlowError::Entry(format!(
                "rescue index {} not found; rebuild with `esperanto index` (it always builds the rescue indexes; the collapsed-rescue channel refuses to degrade silently)",
                path.display()
            )));
        }
        dirs.push((dir, path, from, to));
    }
    let unm_path = map_dir.join("unmapped.fq.gz");
    if !unm_path.is_file() {
        return Ok(0);
    }

    // Read the unmapped set with the strict FASTQ reader: a truncated or
    // malformed record is a hard error (file, 1-based line, reason), never a
    // silently shortened pool. Read names keep the writer's first-token
    // contract, so `/1`/`/2` mate markers survive intact.
    let mut fq: Vec<RescueRead> = {
        use esperanto_map::fastq::RecordSource as _;
        let mut rdr =
            esperanto_map::fastq::FastqReader::open(&unm_path).map_err(stage_err("map"))?;
        let mut out = Vec::new();
        while let Some(rec) = rdr.next_record().map_err(|e| {
            FlowError::Entry(format!(
                "malformed unmapped FASTQ {}: {e}",
                unm_path.display()
            ))
        })? {
            out.push((
                String::from_utf8_lossy(&rec.name).into_owned(),
                rec.seq,
                rec.qual,
            ));
        }
        out
    };
    // Include primary placements with insufficient aligned coverage or
    // excessive unexplained mismatches. Replace a deferred placement only
    // after a rescue direction passes all evidence checks.
    let n_unmapped = fq.len();
    fq.extend(collect_deferred(raw_bam)?);
    let dropped = dedup_pool(&mut fq);
    if dropped > 0 {
        eprintln!(
            "[rescue] pool: dropped {dropped} duplicate entrant(s) left by an interrupted earlier pass"
        );
    }
    let n_deferred = fq.len() - n_unmapped;
    if fq.is_empty() {
        return Ok(0);
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

    // Persist rescue counts after the mapping statistics have been written.
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

    // Same aligner as the map stage runs here: size worker stacks alike
    // (see WORKER_STACK in esperanto_map::pipeline).
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .stack_size(64 << 20)
        .build()
        .map_err(|e| FlowError::Entry(format!("rescue pool: {e}")))?;

    // One pass per editing direction. A read verified in both directions is
    // direction-ambiguous and is dropped rather than placed twice.
    // BTree, not Hash: the rescued records' order reaches raw.bam and
    // rescued.bed, and both artifacts are byte-compared between runs.
    #[allow(clippy::type_complexity)]
    let mut winner: std::collections::BTreeMap<
        usize,
        (esperanto_map::mapq::ReadAlignment, Vec<u32>, (u8, u8)),
    > = std::collections::BTreeMap::new();
    let mut ambiguous: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut placed_count = 0usize;
    let mut rejected = 0u64;
    // Count rejection reasons separately by strand.
    let mut why: std::collections::BTreeMap<(&'static str, bool), u64> =
        std::collections::BTreeMap::new();
    for (dir, ipath, from, to) in &dirs {
        let idx = index_io::load(ipath).map_err(stage_err("map"))?;
        // Seed k/w and the chain anchor reward must match the transformed
        // index. `edit_dir` stays off: under the matching transform this
        // direction's edits are already identities, and exempting anything
        // further would re-open the alphabet collapse this pass exists to
        // avoid.
        let mut cfg = esperanto_map::align::AlignConfig {
            seed: idx.params,
            ..esperanto_map::align::AlignConfig::rna_default()
        };
        cfg.chain.k = idx.params.k as i32;
        cfg.extend.edit_dir = esperanto_map::extend::EditDir::Off;
        // Ungapped by construction rather than by refusing the gapped answer
        // afterwards -- see [`RESCUE_GAP_PENALTY`].
        cfg.extend.gap_open = RESCUE_GAP_PENALTY;
        cfg.extend.gap_ext = RESCUE_GAP_PENALTY;
        let (lo_from, lo_to) = (from.to_ascii_lowercase(), to.to_ascii_lowercase());
        let results: Vec<Option<esperanto_map::mapq::ReadAlignment>> = pool.install(|| {
            use rayon::prelude::*;
            cand_idx
                .par_iter()
                .map_init(
                    || esperanto_map::align::Aligner::new(&idx, cfg),
                    |al, &i| {
                        let t: Vec<u8> = fq[i]
                            .1
                            .iter()
                            .map(|&b| {
                                if b == *from {
                                    *to
                                } else if b == lo_from {
                                    lo_to
                                } else {
                                    b
                                }
                            })
                            .collect();
                        al.align_read(&t)
                    },
                )
                .collect()
        });
        for (&i, aln) in cand_idx.iter().zip(results) {
            let Some(mut a) = aln else { continue };
            placed_count += 1;
            let (_, seq, qual) = &fq[i];
            let proof = match verify_rescue_placement(&a, seq, qual, true_ref, *dir) {
                Ok(edit_pos0) => edit_pos0,
                Err(gate) => {
                    let minus = a.strand == esperanto_map::seed::Strand::Minus;
                    *why.entry((gate, minus)).or_insert(0) += 1;
                    rejected += 1;
                    continue;
                }
            };
            // `rescued` caps MAPQ at the pileup gate (see `mapq::mapq_of`).
            // The margin the transformed index reported is kept rather than
            // zeroed: it is the only uniqueness evidence this placement has,
            // and zeroing it put the read below the caller's own gate, which
            // made the whole rescue unable to produce a call.
            a.rescued = true;
            if winner.insert(i, (a, proof, *dir)).is_some() {
                ambiguous.insert(i);
            }
        }
    }
    for i in &ambiguous {
        winner.remove(i);
        rejected += 2; // both passes' placements discarded
        *why.entry(("direction_ambiguous", false)).or_insert(0) += 2;
    }

    // Phase 2: paired-end concordance, once, over the four-letter placements
    // *plus* the rescue's own (see `rescue_mate_gate`).
    //
    // The table outlives the gate because the BAM writer needs the same view:
    // a fragment whose mates were *both* rescued has no four-letter mate to
    // point at, and there is nowhere else to learn where its mate landed.
    let mut placed = mates.as_ref().map(|m| {
        let mut placed = m.clone();
        for (i, (a, _, _)) in &winner {
            let (base, slot) = split_mate_marker(&fq[*i].0);
            placed
                .entry(base.as_bytes().to_vec())
                .or_insert([None, None])[slot] = Some(MateLoc {
                tid: a.contig as i32,
                pos: a.pos as i64,
                reverse: a.strand == esperanto_map::seed::Strand::Minus,
            });
        }
        placed
    });
    if let Some(placed) = placed.as_mut() {
        let discordant: Vec<usize> = winner
            .iter()
            .filter_map(|(i, (a, _, _))| {
                rescue_mate_gate(a, &fq[*i].0, placed).err().map(|gate| {
                    let minus = a.strand == esperanto_map::seed::Strand::Minus;
                    *why.entry((gate, minus)).or_insert(0) += 1;
                    *i
                })
            })
            .collect();
        rejected += discordant.len() as u64;
        for i in discordant {
            // Undo this read's contribution to the table, restoring whatever
            // the four-letter pass had in the slot, or the surviving mate gets
            // written pointing at a record that was dropped.
            let (base, slot) = split_mate_marker(&fq[i].0);
            let key = base.as_bytes();
            let four_letter = mates
                .as_ref()
                .and_then(|m| m.get(key))
                .and_then(|s| s[slot]);
            if let Some(slots) = placed.get_mut(key) {
                slots[slot] = four_letter;
            }
            winner.remove(&i);
        }
    }

    // Two keys, two jobs: `rescued_names` decides which records leave
    // unmapped.fq.gz (names as that file writes them), while
    // `rescued_slots` decides which raw.bam records the rescued placement
    // replaces — mate-resolved, since the two mates share a qname and only
    // one of them is being replaced.
    let mut rescued_names: std::collections::HashSet<&str> = std::collections::HashSet::new();
    let mut rescued_slots: std::collections::HashSet<(Vec<u8>, usize)> =
        std::collections::HashSet::new();
    #[allow(clippy::type_complexity)]
    let mut rescued: Vec<(
        String,
        Vec<u8>,
        Vec<u8>,
        esperanto_map::mapq::ReadAlignment,
        Vec<u32>,
        (u8, u8),
    )> = Vec::new();
    for (i, (a, proof, dir)) in winner {
        let (name, seq, qual) = &fq[i];
        let (base, slot) = split_mate_marker(name);
        rescued_slots.insert((base.as_bytes().to_vec(), slot));
        rescued.push((name.clone(), seq.clone(), qual.clone(), a, proof, dir));
        rescued_names.insert(name.as_str());
    }
    let rescued_n = rescued.len() as u64;
    eprintln!(
        "[rescue] collapsed: {} verified ({} placed, {} rejected, {} screened) \
         of {} offered ({} unmapped + {} deferred)",
        rescued.len(),
        placed_count,
        rejected,
        prefiltered,
        fq.len(),
        n_unmapped,
        n_deferred
    );
    if !why.is_empty() {
        let mut rows: Vec<_> = why.iter().collect();
        rows.sort_by_key(|(_, n)| std::cmp::Reverse(**n));
        let parts: Vec<String> = rows
            .iter()
            .map(|((gate, minus), n)| format!("{gate}{}={n}", if *minus { "(-)" } else { "(+)" }))
            .collect();
        eprintln!("[rescue] rejected by: {}", parts.join(" "));
    }
    if rescued.is_empty() {
        // Nothing merged and nothing evicted: the seal key is the whole
        // product of this pass.
        patch_stats(0, rejected)?;
        return Ok(0);
    }
    // Whether the verified placements clear the pileup gate is a separate
    // question from whether they are verified, and the two failure modes look
    // identical downstream: a placement held at 0 by the formula (too few
    // anchors, near-tied runner-up) contributes no call, exactly like one held
    // at 0 by policy. Report the split rather than inferring it from the
    // caller's output.
    {
        let gate = esperanto_map::mapq::PILEUP_GATE;
        let (mut at_gate, mut mid, mut zero) = (0u64, 0u64, 0u64);
        for (_, _, _, a, _, _) in &rescued {
            match esperanto_map::mapq::mapq_of(a) {
                q if q >= gate => at_gate += 1,
                0 => zero += 1,
                _ => mid += 1,
            }
        }
        eprintln!(
            "[rescue] mapq: {at_gate} at gate ({gate}), {mid} in 1..{}, {zero} at 0",
            gate - 1
        );
    }

    // Persist reference-based edit positions and support for each rescued
    // read. The resume validator checks this record against raw.bam.
    {
        let mut text = String::new();
        use std::fmt::Write as _;
        text.push_str("qname\tchrom\tpos0\tstrand\tdir\tn_edits\tedit_pos0\tmapq\n");
        for (name, _, _, a, proof, dir) in &rescued {
            let positions: Vec<String> = proof.iter().map(|p| p.to_string()).collect();
            let _ = writeln!(
                text,
                "{}\t{}\t{}\t{}\t{}{}\t{}\t{}\t{}",
                name,
                true_ref.contigs[a.contig as usize].name,
                a.pos,
                if a.strand == esperanto_map::seed::Strand::Minus {
                    '-'
                } else {
                    '+'
                },
                dir.0 as char,
                dir.1 as char,
                proof.len(),
                positions.join(","),
                esperanto_map::mapq::mapq_of(a),
            );
        }
        fs::write(map_dir.join("rescued.tsv"), text)?;
    }

    // Placement sidecar for the report's hyperedited-region track: one
    // `chrom<TAB>pos` row (0-based) per rescued read, written before any
    // merge step so an interruption still leaves the raw data on disk.
    {
        let mut text = String::new();
        use std::fmt::Write as _;
        for (_, _, _, a, _, _) in &rescued {
            let _ = writeln!(
                text,
                "{}\t{}",
                true_ref.contigs[a.contig as usize].name, a.pos
            );
        }
        fs::write(map_dir.join("rescued.bed"), text)?;
    }

    // Merge rescued records into raw.bam (same contig order: both indices
    // derive from the same FASTA) -> replace raw.bam before sort.
    let header = {
        let rdr = rust_htslib::bam::Reader::from_path(raw_bam).map_err(stage_err("map"))?;
        rust_htslib::bam::Header::from_template(rdr.header())
    };
    let merged = raw_bam.with_extension("merged.bam");
    {
        let mut w =
            rust_htslib::bam::Writer::from_path(&merged, &header, rust_htslib::bam::Format::Bam)
                .map_err(stage_err("map"))?;
        {
            let mut rdr = rust_htslib::bam::Reader::from_path(raw_bam).map_err(stage_err("map"))?;
            for r in rdr.records() {
                let rec = r.map_err(stage_err("map"))?;
                // Drop the rescued reads' original records — unmapped ones and
                // clipped placements alike — so exactly one record per read
                // survives. Keyed by mate slot, not by name: the mate keeps
                // its own record. A stale RE:Z:collapsed record from an
                // interrupted earlier pass whose read did not verify this
                // time is dropped too, so raw.bam, rescued.tsv and the seal
                // key always converge to the same rescued set.
                if merge_drops_record(&rec, &rescued_slots) {
                    continue;
                }
                w.write(&rec).map_err(stage_err("map"))?;
            }
        }
        let header_view = rust_htslib::bam::HeaderView::from_header(&header);
        for (name, seq, qual, a, _, _) in &rescued {
            let mut rec = rust_htslib::bam::Record::new();
            // The `/1`,`/2` marker belongs to unmapped.fq.gz, not to a BAM
            // qname: a record carrying it would never associate with its mate.
            // The pairing goes into the flags instead.
            let (base, slot) = split_mate_marker(name);
            let qname = base.as_bytes();
            let rev = a.strand == esperanto_map::seed::Strand::Minus;
            let (s2, q2) = esperanto_bamio::apply_t13(rev, seq, qual);
            let tid = header_view
                .tid(true_ref.contigs[a.contig as usize].name.as_bytes())
                .ok_or_else(|| FlowError::Entry("rescued contig missing from BAM header".into()))?;
            rec.set(qname, None, &s2, &esperanto_map::bam::phred_from_ascii(&q2));
            rec.set_tid(tid as i32);
            rec.set_pos(a.pos as i64);
            // Preserve the transformed-index margin, capped by mapq_of at
            // PILEUP_GATE for rescued placements.
            rec.set_mapq(esperanto_map::mapq::mapq_of(a));
            let mut flag = 0u16;
            if rev {
                flag |= 0x10;
            }
            // Pair flags follow the input type and template slot even when
            // the mate has no accepted placement.
            if pe {
                flag |= 0x1 | if slot == 0 { 0x40 } else { 0x80 };
                match placed
                    .as_ref()
                    .and_then(|m| m.get(qname))
                    .and_then(|slots| slots[1 - slot].as_ref())
                {
                    Some(m) => {
                        if m.reverse {
                            flag |= 0x20;
                        }
                        rec.set_mtid(m.tid);
                        rec.set_mpos(m.pos);
                    }
                    None => {
                        // Convention for an unplaced mate: flag it, and park
                        // it on this record's own coordinate.
                        flag |= 0x8;
                        rec.set_mtid(tid as i32);
                        rec.set_mpos(a.pos as i64);
                    }
                }
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
            // Record::set requires Phred values without the ASCII offset.
            rec.set(
                qname,
                Some(&rust_htslib::bam::record::CigarString(cigar)),
                &s2,
                &esperanto_map::bam::phred_from_ascii(&q2),
            );
            rec.push_aux(b"RE", rust_htslib::bam::record::Aux::String("collapsed"))
                .map_err(stage_err("map"))?;
            w.write(&rec).map_err(stage_err("map"))?;
        }
    }
    fs::rename(&merged, raw_bam)?;

    // Rewrite unmapped.fq.gz without the rescued reads. This eviction happens
    // only after the merge is durable: until the rename above lands, the
    // unmapped pool is the only place an unmapped-origin rescued read still
    // exists, and an interrupted pass must leave the pool complete for the
    // re-run that resume schedules.
    {
        let mut out = Vec::new();
        for (name, seq, qual) in &fq[..n_unmapped] {
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

    // The seal key is written last. A key in align_qc.json must mean the
    // merge, the sidecars and the pool eviction have all completed; writing
    // it earlier let a killed run resume past a raw.bam that never received
    // the rescued records, and every downstream stage then read the short
    // BAM as if it were whole.
    patch_stats(rescued_n, rejected)?;
    Ok(rescued_n)
}

/// Whether the rescue merge replaces a raw.bam record: the original record of
/// a read (re-)verified this pass — keyed by mate slot, not by name, since the
/// two mates share a qname and only one of them is being replaced — or a stale
/// RE:Z:collapsed record left by an interrupted earlier pass whose read did
/// not verify this time. raw.bam, rescued.tsv and the seal key always converge
/// to the same rescued set.
fn merge_drops_record(
    rec: &rust_htslib::bam::Record,
    rescued_slots: &std::collections::HashSet<(Vec<u8>, usize)>,
) -> bool {
    let slot = usize::from(rec.is_paired() && !rec.is_first_in_template());
    if rescued_slots.contains(&(rec.qname().to_vec(), slot)) {
        return true;
    }
    matches!(
        rec.aux(b"RE"),
        Ok(rust_htslib::bam::record::Aux::String("collapsed"))
    )
}

/// Mapped-mate placement slot for the paired-end concordance check.
#[derive(Clone, Copy)]
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

/// Aligned fraction below which a Tier-1 placement counts as only partly
/// explaining its read, and is therefore offered to the rescue passes.
const RESCUE_CLIPPED_COVERAGE: f64 = 0.7;

/// A read offered to the rescue: name, sequence and quality as sequenced. The
/// two sources — `unmapped.fq.gz` and [`collect_deferred`] — produce the same
/// shape so the passes below cannot tell them apart.
type RescueRead = (String, Vec<u8>, Vec<u8>);

/// Return primary placements for rescue in their sequenced orientation.
/// Eligible reads have aligned coverage below RESCUE_CLIPPED_COVERAGE or
/// ordinary-mismatch density at least HYPEREDIT_MIN_DENSITY. Rescued records
/// from an interrupted merge are included so the pass can be repeated.
/// The EK ordinary-mismatch count excludes editing-class substitutions and
/// is a lower bound on the unexplained density used by the MAPQ cap.
fn collect_deferred(raw_bam: &Path) -> Result<Vec<RescueRead>, FlowError> {
    use rust_htslib::bam::record::Cigar;
    let mut rdr = rust_htslib::bam::Reader::from_path(raw_bam).map_err(stage_err("map"))?;
    let header = rdr.header().clone();
    let mut out = Vec::new();
    for r in rdr.records() {
        let rec = r.map_err(stage_err("map"))?;
        if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() {
            continue;
        }
        // Decoy placements are explained reads: their source is absent from
        // the primary reference, so the collapsed-alphabet rescue has nothing
        // to find for them.
        if rec.tid() >= 0 {
            let name = header.tid2name(rec.tid() as u32);
            if name.starts_with(b"decoy_") {
                continue;
            }
        }
        let len = rec.seq_len();
        if len == 0 {
            continue;
        }
        // A record sealed by an earlier rescue pass carries RE:Z:collapsed.
        // A first pass has none; on a rescue re-run after an interrupted
        // merge these records are the only place the rescued reads still
        // exist, so they re-enter the pool and are verified again like any
        // other offered read.
        if matches!(
            rec.aux(b"RE"),
            Ok(rust_htslib::bam::record::Aux::String("collapsed"))
        ) {
            out.push(pool_record(&rec));
            continue;
        }
        // Aligned columns, matching `ReadAlignment::aligned_and_matching` so
        // the density below is the same quantity the MAPQ cap used.
        let aligned: u32 = rec
            .cigar()
            .iter()
            .filter_map(|op| match op {
                Cigar::Match(n)
                | Cigar::Equal(n)
                | Cigar::Diff(n)
                | Cigar::Ins(n)
                | Cigar::Del(n) => Some(*n),
                _ => None,
            })
            .sum();
        let covered: u32 = rec
            .cigar()
            .iter()
            .filter_map(|op| match op {
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Ins(n) => Some(*n),
                _ => None,
            })
            .sum();
        if f64::from(covered) >= RESCUE_CLIPPED_COVERAGE * len as f64
            && !unexplained_by_editing(&rec, aligned)
        {
            continue;
        }
        out.push(pool_record(&rec));
    }
    Ok(out)
}

/// Keep the first pool entry per (template name, mate slot). Unmapped reads
/// precede deferred reads, preserving the prefix used during FASTQ eviction.
/// This also prevents duplicate records after an interrupted BAM merge.
fn dedup_pool(fq: &mut Vec<RescueRead>) -> usize {
    let mut seen = std::collections::HashSet::new();
    let before = fq.len();
    fq.retain(|(name, _, _)| {
        let (base, slot) = split_mate_marker(name);
        seen.insert((base.as_bytes().to_vec(), slot))
    });
    before - fq.len()
}

/// A raw.bam record as a pool entry: sequence and quality as sequenced
/// (BAM stores reference orientation; `apply_t13` is an involution, so
/// applying it again recovers the read as sequenced — which is what the
/// rescue transform and `verify_rescue_placement` both expect), named the
/// way `unmapped.fq.gz` names its records so both sources share one code
/// path downstream.
fn pool_record(rec: &rust_htslib::bam::Record) -> RescueRead {
    let qual: Vec<u8> = rec.qual().iter().map(|q| q.saturating_add(33)).collect();
    let (seq, qual) = esperanto_bamio::apply_t13(rec.is_reverse(), &rec.seq().as_bytes(), &qual);
    let name = String::from_utf8_lossy(rec.qname()).into_owned();
    let name = if rec.is_paired() {
        format!("{name}/{}", if rec.is_first_in_template() { 1 } else { 2 })
    } else {
        name
    };
    (name, seq, qual)
}

/// Whether a placement's ordinary-mismatch density reaches
/// [`esperanto_map::mapq::HYPEREDIT_MIN_DENSITY`], read from the `EK` tag's
/// `mm` field. A record without `EK` — nothing in the map stage writes one,
/// but a re-run over a foreign BAM could — counts as explained, so an absent
/// tag never enlarges the pool.
fn unexplained_by_editing(rec: &rust_htslib::bam::Record, aligned: u32) -> bool {
    if aligned == 0 {
        return false;
    }
    let Ok(rust_htslib::bam::record::Aux::String(ek)) = rec.aux(b"EK") else {
        return false;
    };
    let Some(mm) = ek
        .split(';')
        .find_map(|f| f.strip_prefix("mm="))
        .and_then(|v| v.parse::<u32>().ok())
    else {
        return false;
    };
    f64::from(mm) / f64::from(aligned) >= esperanto_map::mapq::HYPEREDIT_MIN_DENSITY
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
        qs.iter()
            .map(|&q| (q.saturating_sub(33)) as f64)
            .sum::<f64>()
            / len as f64
    } else {
        kept.iter()
            .map(|&q| (q.saturating_sub(33)) as f64)
            .sum::<f64>()
            / kept.len() as f64
    };
    mean >= RESCUE_MIN_MEAN_PHRED
}

/// Verify a placement against the four-letter reference. Only the pass's
/// editing direction is exempt from ordinary-mismatch counting. On the
/// minus strand that direction is complemented to reference orientation.
/// Accepted reads must satisfy the quality, purity, density, geometry,
/// and ungapped-alignment conditions below. Errors name the rejection gate.
fn verify_rescue_placement(
    aln: &esperanto_map::mapq::ReadAlignment,
    seq: &[u8],
    qual: &[u8],
    true_ref: &esperanto_map::fasta::Reference,
    dir: (u8, u8),
) -> Result<Vec<u32>, &'static str> {
    use esperanto_map::extend::CigarOp;
    use esperanto_map::fasta::Base;

    let read_len = seq.len();
    if read_len == 0 || qual.len() != read_len {
        return Err("malformed_read");
    }
    let Some(contig) = true_ref.contigs.get(aln.contig as usize) else {
        return Err("bad_contig");
    };
    // Reject gaps before walking reference coordinates.
    let mut ref_span = 0u32;
    for op in &aln.cigar {
        match op {
            CigarOp::Ins(_) | CigarOp::Del(_) => return Err("gapped"),
            CigarOp::Match(n) | CigarOp::RefSkip(n) => ref_span = ref_span.saturating_add(*n),
            CigarOp::SoftClip(_) => {}
        }
    }
    if aln.pos as u64 + ref_span as u64 > contig.len as u64 {
        return Err("off_contig");
    }

    // Read bases and qualities in reference orientation, and the editing class
    // that orientation puts the pass's own direction into.
    let minus = aln.strand == esperanto_map::seed::Strand::Minus;
    let (rseq, rqual) = esperanto_bamio::apply_t13(minus, seq, qual);
    // The expectation pair is stated for a plus-frame placement; a minus
    // placement flips it to the complement pair.
    let (expect_from, expect_to) = if minus {
        (complement_base(dir.0), complement_base(dir.1))
    } else {
        dir
    };

    let mut non_edit = 0usize;
    let mut total_mm = 0usize;
    let mut edit_pos: Vec<usize> = Vec::new();
    let mut edit_rpos: Vec<u32> = Vec::new();
    let mut rpos = aln.pos;
    let mut qpos = 0usize;
    for op in &aln.cigar {
        match op {
            CigarOp::Match(n) => {
                for _ in 0..*n {
                    if qpos >= read_len {
                        // malformed CIGAR: consumes more read than exists
                        return Err("cigar_overrun");
                    }
                    let rb = Base::from_ascii(rseq[qpos]);
                    let fb = contig.base(rpos);
                    if rb == Base::N || fb == Base::N {
                        total_mm += 1;
                        non_edit += 1;
                    } else if rb != fb {
                        total_mm += 1;
                        let editing = fb == Base::from_ascii(expect_from)
                            && rb == Base::from_ascii(expect_to);
                        if !editing {
                            // Anything that is not this direction's edit is an
                            // ordinary mismatch, and the transformed-space edit
                            // distance is exactly the count of these.
                            non_edit += 1;
                        } else if rqual[qpos].saturating_sub(33) >= RESCUE_MIN_EDIT_PHRED {
                            edit_pos.push(qpos);
                            edit_rpos.push(rpos);
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
        return Err("non_edit_mm");
    }
    let ne = edit_pos.len();
    if ne == 0 {
        return Err("no_edits");
    }
    let nlen = read_len as f64;
    if (ne as f64) < RESCUE_MIN_EDIT_DENSITY * nlen {
        return Err("edit_density");
    }
    let purity = if read_len <= RESCUE_SHORT_READ_MAX_LEN {
        RESCUE_PURITY_SHORT
    } else {
        RESCUE_PURITY_LONG
    };
    if ne as f64 <= purity * total_mm as f64 {
        return Err("purity");
    }
    let first = edit_pos[0];
    let last = edit_pos[ne - 1];
    if ((last - first) as f64) < RESCUE_MIN_CLUSTER_SPAN_FRAC * nlen {
        return Err("cluster_span");
    }
    if (last as f64) < RESCUE_END_EXCLUSION_FRAC * nlen {
        return Err("cluster_at_start");
    }
    if (first as f64) > (1.0 - RESCUE_END_EXCLUSION_FRAC) * nlen {
        return Err("cluster_at_end");
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
        return Err("cluster_composition");
    }
    Ok(edit_rpos)
}

fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// Require the mate on the same contig within RESCUE_PE_MATE_WINDOW and
/// in the opposite orientation. Check after all passes using both primary
/// and rescued placements so two rescued mates can support each other.
fn rescue_mate_gate(
    aln: &esperanto_map::mapq::ReadAlignment,
    name: &str,
    placed: &std::collections::HashMap<Vec<u8>, [Option<MateLoc>; 2]>,
) -> Result<(), &'static str> {
    let (base, slot) = split_mate_marker(name);
    let minus = aln.strand == esperanto_map::seed::Strand::Minus;
    let Some(mate) = placed
        .get(base.as_bytes())
        .and_then(|slots| slots[1 - slot].as_ref())
    else {
        return Err("mate_unmapped");
    };
    if mate.tid != aln.contig as i32 {
        return Err("mate_contig");
    }
    if (mate.pos - aln.pos as i64).abs() > RESCUE_PE_MATE_WINDOW {
        return Err("mate_distance");
    }
    if mate.reverse == minus {
        return Err("mate_orientation");
    }
    Ok(())
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

fn stage_scan(params: &RunParams, bam: &Path) -> Result<PathBuf, FlowError> {
    eprintln!("[scan] running");
    let dir = params.out_dir.join("scan");
    fs::create_dir_all(&dir)?;
    let out = dir.join("candidates.bed");
    let count = esperanto_score::unified::discover_sites(bam, &params.fasta, &out, params.threads)
        .map_err(anyhow_stage_err("scan"))?;
    eprintln!("[scan] {count} observation sites");
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

/// Score sites in input order and persist the unified classifier output.
fn stage_score(
    params: &RunParams,
    bam: &Path,
    sites: &[(String, i64)],
) -> Result<(Vec<f64>, f64), FlowError> {
    let dir = params.out_dir.join("score");
    fs::create_dir_all(&dir)?;
    let threshold = load_primary_threshold(&params.bundle)
        .map_err(|e| anyhow_stage_err("score")(anyhow::anyhow!("{e}")))?;
    let scores = if sites.is_empty() {
        Vec::new()
    } else {
        eprintln!("[score] running");
        score_pipeline::score_sites_batched(
            bam,
            &params.fasta,
            &params.bundle,
            sites,
            params.threads,
        )
        .map_err(anyhow_stage_err("score"))?
    };
    let mut text = String::new();
    for ((chrom, pos), score) in sites.iter().zip(&scores) {
        use std::fmt::Write as _;
        let _ = writeln!(text, "{chrom}\t{pos}\t{score}");
    }
    fs::write(dir.join("scores.tsv"), text)?;
    eprintln!("[score] {} sites", scores.len());
    Ok((scores, threshold))
}

/// Reload unified scores when resuming after the score stage.
fn load_score_state(params: &RunParams) -> Result<(Vec<f64>, f64), FlowError> {
    let scores = load_scores(&params.out_dir.join("score").join("scores.tsv"))?;
    let threshold = load_primary_threshold(&params.bundle)
        .map_err(|e| anyhow_stage_err("score")(anyhow::anyhow!("{e}")))?;
    Ok((scores, threshold))
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
    fn hyper_read(
        r: &esperanto_map::fasta::Reference,
        start: usize,
        lo: usize,
        hi: usize,
    ) -> Vec<u8> {
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
        homo.extend_from_slice(
            &b"CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".repeat(3)[..125],
        );
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
        assert!(verify_rescue_placement(&aln, &seq, &qual, r, (b'A', b'G')).is_ok());
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
        // Reference-forward the edits are A>G, so the reverse-complemented
        // sequencer view shows T>C: `ag = false`, and `expect_ag = ag XOR
        // minus` puts the comparison back into the A>G class. Taking the class
        // from the pass's own direction instead — `ag = true` here — must
        // fail, or the XOR is doing nothing.
        assert!(verify_rescue_placement(&aln, &seq, &qual, r, (b'T', b'C')).is_ok());
        assert!(verify_rescue_placement(&aln, &seq, &qual, r, (b'A', b'G')).is_err());
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
        assert!(verify_rescue_placement(&aln, &seq, &qual, r, (b'A', b'G')).is_err());
    }

    #[test]
    fn verify_rejects_gapped_placement() {
        // Indel placements fail the ungapped-alignment requirement.
        let r = reference_with_period(7);
        let seq = hyper_read(r, 20, 20, 110);
        let mut aln = alignment_at(20, 0);
        aln.cigar = vec![
            esperanto_map::extend::CigarOp::Match(70),
            esperanto_map::extend::CigarOp::Del(1),
            esperanto_map::extend::CigarOp::Match(79),
        ];
        let qual = vec![b'I'; 150];
        assert_eq!(
            verify_rescue_placement(&aln, &seq, &qual, r, (b'A', b'G')),
            Err("gapped")
        );
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
        assert!(verify_rescue_placement(&aln, &seq, &qual, r, (b'A', b'G')).is_err());
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
        assert!(verify_rescue_placement(&aln60, &good, &qual60, r, (b'A', b'G')).is_ok());
        // 3 editing + 2 non-editing mismatches: purity 3/5 fails the 80%
        // short-read bar.
        let mut poor = good.clone();
        poor[5] = if poor[5] == b'A' { b'G' } else { b'A' }; // non-editing class
        poor[35] = if poor[35] == b'A' { b'G' } else { b'A' };
        // keep exactly three edits
        poor[29] = if poor[29] == b'G' { b'A' } else { poor[29] };
        assert!(verify_rescue_placement(&aln60, &poor, &qual60, r, (b'A', b'G')).is_err());
    }

    #[test]
    fn verify_pe_mate_concordance() {
        let aln = alignment_at(20, 150);
        let one = |loc: MateLoc| {
            let mut m = std::collections::HashMap::new();
            m.insert(b"t".to_vec(), [Some(loc), None]);
            m
        };
        let concordant = one(MateLoc {
            tid: 0,
            pos: 300,
            reverse: true,
        });
        assert!(rescue_mate_gate(&aln, "t/2", &concordant).is_ok());
        // Same orientation, beyond the window, other contig, absent.
        let same = one(MateLoc {
            tid: 0,
            pos: 300,
            reverse: false,
        });
        assert_eq!(
            rescue_mate_gate(&aln, "t/2", &same),
            Err("mate_orientation")
        );
        let far = one(MateLoc {
            tid: 0,
            pos: 600_000,
            reverse: true,
        });
        assert_eq!(rescue_mate_gate(&aln, "t/2", &far), Err("mate_distance"));
        let other = one(MateLoc {
            tid: 3,
            pos: 300,
            reverse: true,
        });
        assert_eq!(rescue_mate_gate(&aln, "t/2", &other), Err("mate_contig"));
        let empty = std::collections::HashMap::new();
        assert_eq!(rescue_mate_gate(&aln, "t/2", &empty), Err("mate_unmapped"));
    }

    /// A minimal BAM writer for pool/merge tests: one contig, records as
    /// specified. `(qname, flags, cigar, re_sealed)`.
    fn write_pool_bam(
        path: &std::path::Path,
        contig: &str,
        contig_len: u32,
        specs: &[(&str, u16, Vec<rust_htslib::bam::record::Cigar>, bool)],
    ) {
        let mut hdr = rust_htslib::bam::Header::new();
        let mut sq = rust_htslib::bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", contig);
        sq.push_tag(b"LN", contig_len);
        hdr.push_record(&sq);
        let mut w =
            rust_htslib::bam::Writer::from_path(path, &hdr, rust_htslib::bam::Format::Bam).unwrap();
        for (i, (name, flags, cigar, re)) in specs.iter().enumerate() {
            let mut rec = rust_htslib::bam::Record::new();
            let seq = b"AAACCCGG";
            let qual = [30u8; 8];
            let cs = rust_htslib::bam::record::CigarString(cigar.clone());
            rec.set(name.as_bytes(), Some(&cs), seq, &qual);
            rec.set_tid(0);
            rec.set_pos(100 + i as i64);
            rec.set_flags(*flags);
            if *re {
                rec.push_aux(b"RE", rust_htslib::bam::record::Aux::String("collapsed"))
                    .unwrap();
            }
            w.write(&rec).unwrap();
        }
    }

    #[test]
    fn merge_drop_rule_replaces_winners_and_stale_seals() {
        let d = std::env::temp_dir().join(format!("esp-merge-drop-{}", std::process::id()));
        fs::create_dir_all(&d).unwrap();
        let bam = d.join("raw.bam");
        use rust_htslib::bam::record::Cigar;
        write_pool_bam(
            &bam,
            "chr1",
            10000,
            &[
                ("keep", 0, vec![Cigar::Match(8)], false),
                ("orig", 0, vec![Cigar::Match(8)], false),
                ("stale", 0, vec![Cigar::Match(8)], true),
                ("rewin", 0, vec![Cigar::Match(8)], true),
            ],
        );
        let mut slots: std::collections::HashSet<(Vec<u8>, usize)> =
            std::collections::HashSet::new();
        slots.insert((b"orig".to_vec(), 0));
        slots.insert((b"rewin".to_vec(), 0));
        use rust_htslib::bam::Read as _;
        let mut rdr = rust_htslib::bam::Reader::from_path(&bam).unwrap();
        let mut kept = Vec::new();
        for r in rdr.records() {
            let rec = r.unwrap();
            let name = String::from_utf8_lossy(rec.qname()).into_owned();
            assert_eq!(merge_drops_record(&rec, &slots), name != "keep", "{name}");
            if !merge_drops_record(&rec, &slots) {
                kept.push(name);
            }
        }
        assert_eq!(kept, vec!["keep".to_string()]);
        let _ = fs::remove_dir_all(&d);
    }

    #[test]
    fn collect_deferred_recovers_sealed_records() {
        let d = std::env::temp_dir().join(format!("esp-deferred-re-{}", std::process::id()));
        fs::create_dir_all(&d).unwrap();
        let bam = d.join("raw.bam");
        use rust_htslib::bam::record::Cigar;
        // flags: 0x1 paired, 0x40 first in template, 0x10 reverse.
        write_pool_bam(
            &bam,
            "chr1",
            10000,
            &[
                ("sealed", 0x1 | 0x40 | 0x10, vec![Cigar::Match(8)], true),
                ("clean", 0, vec![Cigar::Match(8)], false),
                (
                    "clipped",
                    0,
                    vec![Cigar::Match(2), Cigar::SoftClip(6)],
                    false,
                ),
            ],
        );
        let pool = collect_deferred(&bam).unwrap();
        let names: Vec<&str> = pool.iter().map(|(n, _, _)| n.as_str()).collect();
        // The clean placement stays out; the sealed record and the clipped
        // placement are offered.
        assert!(names.contains(&"sealed/1"));
        assert!(names.contains(&"clipped"));
        assert!(!names.iter().any(|n| n.starts_with("clean")));
        // The sealed record was stored reverse-complemented; the pool entry
        // must carry the read as sequenced: revcomp(AAACCCGG) = CCGGGTTT,
        // qualities back in ASCII offset.
        let (_, seq, qual) = pool.iter().find(|(n, _, _)| n == "sealed/1").unwrap();
        assert_eq!(seq.as_slice(), b"CCGGGTTT");
        assert_eq!(qual.as_slice(), &[63u8; 8]);
        let _ = fs::remove_dir_all(&d);
    }

    #[test]
    fn dedup_pool_drops_cross_source_duplicates() {
        let read = |n: &str| (n.to_string(), b"ACGT".to_vec(), b"IIII".to_vec());
        // The unmapped half comes first; the deferred half re-offers x/1
        // (merged by an interrupted pass) and both mates of y.
        let mut fq = vec![
            read("x/1"),
            read("y/1"),
            read("y/2"),
            read("x/1"),
            read("z"),
        ];
        let dropped = dedup_pool(&mut fq);
        assert_eq!(dropped, 1);
        let names: Vec<&str> = fq.iter().map(|(n, _, _)| n.as_str()).collect();
        assert_eq!(names, ["x/1", "y/1", "y/2", "z"]);
    }

    /// RE:Z:collapsed record count in a BAM.
    fn re_count(path: &std::path::Path) -> u64 {
        use rust_htslib::bam::Read as _;
        let mut rdr = rust_htslib::bam::Reader::from_path(path).unwrap();
        let mut n = 0u64;
        for r in rdr.records() {
            let rec = r.unwrap();
            if matches!(
                rec.aux(b"RE"),
                Ok(rust_htslib::bam::record::Aux::String("collapsed"))
            ) {
                n += 1;
            }
        }
        n
    }

    fn rescued_tsv_rows(map_dir: &std::path::Path) -> u64 {
        fs::read_to_string(map_dir.join("rescued.tsv"))
            .unwrap()
            .lines()
            .filter(|l| !l.is_empty())
            .count() as u64
            - 1
    }

    fn sealed_value(map_dir: &std::path::Path) -> Option<u64> {
        let t = fs::read_to_string(map_dir.join("align_qc.json")).unwrap();
        serde_json::from_str::<serde_json::Value>(&t).unwrap()["rescued_collapsed"].as_u64()
    }

    /// End-to-end over the rescue stage: a first pass rescues a hyperedited
    /// read and seals; then the hand-built window state "raw.bam merged,
    /// unmapped.fq.gz eviction pending, seal absent" is rerun, offering the
    /// same read from both pool halves. Every counter must stay singular.
    #[test]
    fn rescue_rerun_after_merge_without_eviction_does_not_double() {
        // Same period-7 fixture as the verify tests.
        let bg = b"CGTATGCCTAGG";
        let mut seq = Vec::with_capacity(400);
        for i in 0..400 {
            seq.push(if i % 7 == 0 { b'A' } else { bg[i % bg.len()] });
        }
        let fasta = format!(">c1\n{}\n", String::from_utf8(seq.clone()).unwrap());
        let true_ref: &'static esperanto_map::fasta::Reference = Box::leak(Box::new(
            esperanto_map::fasta::parse_fasta_bytes(fasta.as_bytes()).unwrap(),
        ));
        let read = hyper_read(true_ref, 20, 20, 110);
        let qual = vec![b'I'; 150];

        let d = std::env::temp_dir().join(format!("esp-rescue-e2e-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        let map_dir = d.join("map");
        fs::create_dir_all(&map_dir).unwrap();
        // Collapsed-alphabet rescue indexes next to a nominal main paidx.
        let paidx = d.join("ref.paidx");
        for (from, to, ext) in [(b'A', b'G', "ag.paidx"), (b'T', b'C', "tc.paidx")] {
            let tseq: Vec<u8> = seq
                .iter()
                .map(|&b| if b == from { to } else { b })
                .collect();
            let tfasta = format!(">c1\n{}\n", String::from_utf8(tseq).unwrap());
            let tref = esperanto_map::fasta::parse_fasta_bytes(tfasta.as_bytes()).unwrap();
            let tidx = esperanto_map::index::Index::build(
                tref,
                esperanto_map::seed::SeedParams { k: 31, w: 10 },
            );
            esperanto_map::index_io::save(&tidx, &paidx.with_extension(ext)).unwrap();
        }
        // Map-stage artifacts: the unmapped pool holds the read; raw.bam
        // starts empty (header only).
        let raw_bam = map_dir.join("raw.bam");
        write_pool_bam(&raw_bam, "c1", 400, &[]);
        let unm = map_dir.join("unmapped.fq.gz");
        let write_unmapped = |with_read: bool| {
            let text = if with_read {
                format!(
                    "@hyper1\n{}\n+\n{}\n",
                    String::from_utf8(read.clone()).unwrap(),
                    String::from_utf8(qual.clone()).unwrap()
                )
            } else {
                String::new()
            };
            let f = fs::File::create(&unm).unwrap();
            let mut e = flate2::write::GzEncoder::new(f, flate2::Compression::new(1));
            use std::io::Write as _;
            e.write_all(text.as_bytes()).unwrap();
            e.finish().unwrap();
        };
        let write_stats = |sealed: Option<u64>| {
            let stats = esperanto_map::stats::AlignStats {
                esperanto_map_version: "test".into(),
                mode: "rna-se".into(),
                total_reads: 1,
                mapped_reads: 0,
                unmapped_reads: 1,
                mapping_rate: 0.0,
                proper_pairs: None,
                insert_mean: None,
                insert_stdev: None,
                mapq_hist: vec![0; 61],
                junctions_total: 0,
                rescued_total: 0,
                rescue_fail_total: 0,
                rescued_collapsed: sealed,
                rescue_rejected_collapsed: None,
                rescue_prefiltered_collapsed: None,
                decoy_reads: 0,
                elapsed_seconds: 0.0,
            };
            fs::write(
                map_dir.join("align_qc.json"),
                serde_json::to_string_pretty(&stats).unwrap(),
            )
            .unwrap();
        };
        write_unmapped(true);
        write_stats(None);

        // First pass: the read is rescued and the seal is written.
        let n1 = rescue_collapsed(&paidx, Some(true_ref), &map_dir, 1, &raw_bam).unwrap();
        assert_eq!(n1, 1);
        assert_eq!(re_count(&raw_bam), 1);
        assert_eq!(rescued_tsv_rows(&map_dir), 1);
        assert_eq!(sealed_value(&map_dir), Some(1));

        // Hand-built window: merge done, eviction pending, seal absent. The
        // rerun sees the read in unmapped.fq.gz AND as an RE record.
        write_unmapped(true);
        write_stats(None);
        let n2 = rescue_collapsed(&paidx, Some(true_ref), &map_dir, 1, &raw_bam).unwrap();
        assert_eq!(n2, 1);
        assert_eq!(re_count(&raw_bam), 1);
        assert_eq!(rescued_tsv_rows(&map_dir), 1);
        assert_eq!(sealed_value(&map_dir), Some(1));
        let _ = fs::remove_dir_all(&d);
    }
}
