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
    let stop_at: Option<Stage> = params
        .stop_after
        .as_deref()
        .map(|s| {
            Stage::parse(s).ok_or_else(|| {
                FlowError::Entry(format!(
                    "bad --stop-after '{s}': expected qc|map|sort|scan|verify|hyper|score|vcf"
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
    // When the map stage ran in this process, its (leaked) index and L1
    // bundle are reused by the verify stage instead of a second 12 GB read.
    let mut map_handles: Option<(&'static esperanto_map::index::Index, Option<Arc<esperanto_engine::L1Index>>)> = None;
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
                let (bam, index, l1) = stage_map(params, entry, &clean1, clean2.as_deref())?;
                current_bam = bam;
                map_handles = Some((index, l1));
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
            if start <= Stage::Verify {
                stage_verify(params, &current_bam, &bed, map_handles)?;
            }
            if stop_at == Some(Stage::Verify) {
                eprintln!("[run] stop-after verify");
                return Ok(());
            }
            if start <= Stage::Hyper {
                stage_hyper(params)?;
            }
            if stop_at == Some(Stage::Hyper) {
                eprintln!("[run] stop-after hyper");
                return Ok(());
            }
            sites
        }
    };

    if start <= Stage::Vcf {
        let probs = if start <= Stage::Score {
            let probs = stage_score(params, &current_bam, &sites)?;
            if stop_at == Some(Stage::Score) {
                eprintln!("[run] stop-after score");
                return Ok(());
            }
            probs
        } else {
            load_scores(&params.out_dir.join("score").join("scores.tsv"))?
        };
        vcf::write_vcf(params, entry, &sites, &probs)?;
    }
    if stop_at == Some(Stage::Vcf) {
        eprintln!("[run] stop-after vcf");
        return Ok(());
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

/// (raw BAM, leaked shared alignment index, L1 bundle when present). The
/// index handle is process-lifetime so `stage_verify` reuses it instead of
/// re-reading 12 GB from disk.
fn stage_map(
    params: &RunParams,
    entry: Entry,
    clean1: &Path,
    clean2: Option<&Path>,
) -> Result<(PathBuf, &'static esperanto_map::index::Index, Option<Arc<esperanto_engine::L1Index>>), FlowError> {
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
        params.rescue_control,
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
    rescue_control: bool,
) -> Result<(PathBuf, &'static esperanto_map::index::Index, Option<Arc<L1Index>>), FlowError> {
    let dir = out_dir;
    fs::create_dir_all(dir)?;
    let index: &'static esperanto_map::index::Index =
        Box::leak(Box::new(index_io::load(paidx).map_err(stage_err("map"))?));
    // Tier 1 is strict four-letter: every substitution costs. Exempting an
    // editing direction here would apply the cost to all 38M reads to benefit
    // the hyperedited fraction, and exempting both would collapse the alphabet
    // outright. Hyperediting is handled by the dedicated three-letter rescue
    // passes below, on the reads Tier 1 refuses — the division of labour in
    // Porath, Carmi & Levanon, Nat Commun 5:4726, 2014.
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
        jkmer: None,
        l1: l1.clone(),
        baln: Some(Box::new(File::create(dir.join("align.baln"))?)),
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
    if rescue_control {
        let n = rescue_collapsed(paidx, Some(index.reference), dir, threads, &raw_bam, true)?;
        let rescued = rescue_collapsed(paidx, Some(index.reference), dir, threads, &raw_bam, false)?;
        let ratio = if rescued > 0 { n as f64 / rescued as f64 } else { 0.0 };
        eprintln!("[rescue] control alphabet: {n} verified (ratio {ratio:.2})");
    } else {
        rescue_collapsed(paidx, Some(index.reference), dir, threads, &raw_bam, false)?;
    }
    Ok((raw_bam, index, l1))
}

// Three-letter hyperediting rescue of unmapped reads (Porath, Carmi &
// Levanon, Nat Commun 5:4726, 2014). The unmapped set is re-aligned once per
// editing direction, each against its own transformed index: reads and
// reference both carry A->G for the plus-strand pass (`<stem>.ag.paidx`) and
// both carry T->C for the minus-strand pass (`<stem>.tc.paidx`). Under the
// matching transform an edit of that direction is an identity, so a
// hyperedited read seeds and aligns normally, while the other three bases
// stay distinct and keep separating paralogs.
//
// One direction per pass, never both: a single index merging A/G *and* T/C
// leaves a two-letter alphabet, in which seeds stop being locators and the two
// transitions — which is what A>G and T>C are — become free, making a
// paralogous repeat copy indistinguishable from an edited sibling.
//
// Every placement is verified against the true (four-letter) reference: a
// hyperedited read must carry at most `RESCUE_MAX_NON_EDIT_MM` mismatches
// outside its own editing direction, plus a dense, high-quality,
// direction-pure mismatch cluster with sane geometry and, for paired data, a
// concordant mapped mate. A read that verifies under *both* directions is
// direction-ambiguous — the signature of a paralog, not a transcript — and is
// discarded, which is the same negative-control logic the reference method
// applies to its five non-A>G transforms.
//
// Verified survivors are written back into raw.bam with MAPQ 0 and an
// `RE:Z:collapsed` tag (repeat-family placement, never confident); everything
// else stays unmapped, and unmapped.fq.gz is rewritten without the verified
// reads only. No-op when neither transformed index exists.

/// The rescue passes: `(read shows A>G?, index extension, transform from,
/// transform to)`. Exactly one transform per pass.
///
/// Four passes, not two, because the transform and the reverse complement do
/// not commute: `t_AG(revcomp(F))` sends a genomic A to T while
/// `revcomp(t_AG(F))` sends it to C. Seeding is canonical — a minimizer is
/// `min(kmer, revcomp(kmer))` (`esperanto_map::seed`) — so the orientation a
/// read is found in is decided by which collapsed alphabet its k-mers live in,
/// and the governing identity is
///
/// ```text
/// revcomp( t_AG(revcomp(F)) ) == t_TC(F)
/// ```
///
/// An A>G-edited read that is *antisense* to the genome therefore only seeds
/// against the T>C index, and vice versa. Each transform consequently needs
/// both indexes: the same-name pass finds the sense-oriented reads, the
/// opposite-name pass the antisense ones. With only the two same-name passes
/// every placement is plus-strand — measured on the simulated corpus, 7987 of
/// 7987 — so one orientation of every transcript was unreachable, and since
/// the two mates of a pair are always in opposite orientations, no hyperedited
/// fragment could ever have both mates placed.
///
/// The cross passes are not a weaker test. A read is still collapsed under one
/// transform only, so the alphabet stays three-letter and the negative-control
/// logic is intact: a read that verifies under both the A>G and the T>C
/// hypothesis is direction-ambiguous and is dropped.
/// Rescue passes as (plus-frame expectation pair, index, read from, read to).
/// `control = true` swaps the editing transforms for the negative-control
/// alphabets (A>C, T>G): transforms no editing chemistry can produce.
fn rescue_dirs(control: bool) -> [((u8, u8), &'static str, u8, u8); 4] {
    if control {
        [
            ((b'A', b'C'), "ac.paidx", b'A', b'C'),
            ((b'A', b'C'), "tg.paidx", b'A', b'C'),
            ((b'T', b'G'), "tg.paidx", b'T', b'G'),
            ((b'T', b'G'), "ac.paidx", b'T', b'G'),
        ]
    } else {
        [
            ((b'A', b'G'), "ag.paidx", b'A', b'G'),
            ((b'A', b'G'), "tc.paidx", b'A', b'G'),
            ((b'T', b'C'), "tc.paidx", b'T', b'C'),
            ((b'T', b'C'), "ag.paidx", b'T', b'C'),
        ]
    }
}

/// Maximum mismatches outside the pass's own editing direction tolerated in a
/// verified placement — the transformed-space edit distance, since under the
/// matching transform that direction's edits cost nothing.
const RESCUE_MAX_NON_EDIT_MM: usize =
    esperanto_map::mapq::HYPEREDIT_MAX_NON_EDIT_MM as usize;
/// Gap penalty that makes the collapsed realignment gapless by construction.
///
/// Porath, Carmi & Levanon (Nat Commun 5:4726, 2014) realign hyperedited reads
/// ungapped, and the reason is structural rather than stylistic: with A and G
/// collapsed to one letter the dynamic program has many equal-scoring paths,
/// so a gap that buys a match it would have had anyway opens for free. Leaving
/// the gap penalties at their four-letter values and then refusing any output
/// that contains a gap answers the wrong question -- it discards the candidate
/// locus instead of asking what the best ungapped placement there is. On the
/// real sample that refusal was 464,703 of 840,066 rejections, 55%, the
/// channel's largest single loss, and the evidence offered for it was measured
/// on the simulated corpus only. Charging the gap up front instead lets the
/// same loci be scored and then policed by [`RESCUE_MAX_NON_EDIT_MM`]: a read
/// that really does carry an indel goes out of register, accumulates ordinary
/// mismatches and is refused on that evidence, while a read whose gap was an
/// artifact of the reduced alphabet keeps its ungapped placement.
///
/// The extension is local, so soft-clipping always scores at least 0 and no
/// gap can ever be on the optimal path once one costs more than the whole read
/// can earn -- a few hundred bases at `match_score` 2. 400 is prohibitive by
/// that measure and still small enough that `Acc for i16` never oversteps its
/// `-30000` floor: every accumulator step moves by at most 400 against 2768
/// bases of headroom, so the clamp in `add` is reached before the type is.
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

/// Collapsed rescue with the optional negative-control pass. The control
/// runs FIRST, while unmapped.fq.gz and raw.bam are pristine; it counts only
/// and never merges. The ratio line lands after the normal pass so both
/// numbers are final.
#[allow(clippy::too_many_arguments)]
fn rescue_with_optional_control(
    params: &RunParams,
    paidx: &Path,
    live_reference: Option<&'static esperanto_map::fasta::Reference>,
    map_dir: &Path,
    raw_bam: &Path,
) -> Result<(), FlowError> {
    let control_verified = if params.rescue_control {
        Some(rescue_collapsed(
            paidx,
            live_reference,
            map_dir,
            params.threads,
            raw_bam,
            true,
        )?)
    } else {
        None
    };
    let rescued = rescue_collapsed(paidx, live_reference, map_dir, params.threads, raw_bam, false)?;
    if let Some(n) = control_verified {
        let ratio = if rescued > 0 {
            n as f64 / rescued as f64
        } else {
            0.0
        };
        eprintln!("[rescue] control alphabet: {n} verified (ratio {ratio:.2})");
    }
    Ok(())
}

fn stage_rescue_collapsed(params: &RunParams, raw_bam: &Path) -> Result<(), FlowError> {
    let Some(paidx) = params.index.clone() else {
        return Ok(());
    };
    rescue_with_optional_control(
        params,
        &paidx,
        None,
        &params.out_dir.join("map"),
        raw_bam,
    )
}

fn rescue_collapsed(
    paidx: &Path,
    live_reference: Option<&'static esperanto_map::fasta::Reference>,
    map_dir: &Path,
    threads: usize,
    raw_bam: &Path,
    control: bool,
) -> Result<u64, FlowError> {
    let dirs: Vec<((u8, u8), PathBuf, u8, u8)> = rescue_dirs(control)
        .iter()
        .map(|&(dir, ext, from, to)| (dir, paidx.with_extension(ext), from, to))
        .filter(|(_, p, _, _)| p.is_file())
        .collect();
    if dirs.is_empty() {
        return Ok(0);
    }
    let unm_path = map_dir.join("unmapped.fq.gz");
    if !unm_path.is_file() {
        return Ok(0);
    }

    // Read the unmapped set.
    let mut fq: Vec<RescueRead> = {
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
    // Tier 1 almost never refuses a read outright — it soft-clips it, or
    // places it whole with a mismatch spectrum it does not stand behind — so
    // `unmapped.fq.gz` is the residue rather than the population this channel
    // is looking for. [`collect_deferred`] adds the placements Tier 1 declines
    // to explain; a deferred read's original record is replaced only if a
    // direction verifies, so nothing is lost by trying.
    let n_unmapped = fq.len();
    fq.extend(collect_deferred(raw_bam)?);
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

    // Honest accounting into align_qc.json; the rescue runs after the
    // stats document is written, so the keys are patched in here (the
    // report stage treats absence as 0).
    let patch_stats = |verified: u64, rejected: u64| -> Result<(), FlowError> {
        let p = map_dir.join("align_qc.json");
        let mut stats: esperanto_map::stats::AlignStats =
            serde_json::from_str(&fs::read_to_string(&p)?).map_err(stage_err("map"))?;
        if control {
            stats.rescue_control_verified = Some(verified);
        } else {
            stats.rescued_collapsed = Some(verified);
            stats.rescue_rejected_collapsed = Some(rejected);
            stats.rescue_prefiltered_collapsed = Some(prefiltered);
        }
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
    // Rejection attribution, keyed by (gate, placement was minus-strand). The
    // strand dimension is not decoration: `apply_t13` presents a minus-strand
    // placement in reference orientation, which relabels the pass's own
    // editing direction, so a gate that closes only on one strand is a
    // different defect from one that closes on both.
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
            let four_letter = mates.as_ref().and_then(|m| m.get(key)).and_then(|s| s[slot]);
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
    if control {
        // Negative-control pass: count only. Nothing is merged, moved, or
        // rewritten; the verified count is the whole product.
        patch_stats(winner.len() as u64, rejected)?;
        eprintln!(
            "[rescue] control alphabet: {} verified ({} rejected)",
            winner.len(),
            rejected
        );
        return Ok(winner.len() as u64);
    }
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
    patch_stats(rescued_n, rejected)?;
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
            .map(|((gate, minus), n)| {
                format!("{gate}{}={n}", if *minus { "(-)" } else { "(+)" })
            })
            .collect();
        eprintln!("[rescue] rejected by: {}", parts.join(" "));
    }
    if rescued.is_empty() {
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

    // Per-read rescue record for the hyperedited-region caller: edit
    // positions were computed by the verification walk and are reference
    // 0-based.
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
                if a.strand == esperanto_map::seed::Strand::Minus { '-' } else { '+' },
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

    // Rewrite unmapped.fq.gz without the rescued reads.
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
                // Drop the rescued reads' original records — unmapped ones and
                // clipped placements alike — so exactly one record per read
                // survives. Keyed by mate slot, not by name: the mate keeps
                // its own record.
                let slot = usize::from(rec.is_paired() && !rec.is_first_in_template());
                if rescued_slots.contains(&(rec.qname().to_vec(), slot)) {
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
            // Not 0. `mapq_of` already caps a rescued placement at
            // `PILEUP_GATE` because it was won in a reduced alphabet, and that
            // cap is the whole policy; overwriting it with 0 put every rescued
            // read below the caller's own `min_mapq`, so the channel could not
            // contribute a call however strictly it had been verified. The
            // margin term still applies underneath the cap, so a rescued
            // placement with a near-tied runner-up stays low on its own.
            rec.set_mapq(esperanto_map::mapq::mapq_of(a));
            let mut flag = 0u16;
            if rev {
                flag |= 0x10;
            }
            // Paired input makes this record one segment of a template
            // whatever became of its mate, so the pair flags follow `pe` and
            // the read's own slot. Conditioning them on a *found* mate instead
            // wrote both mates of an all-rescued fragment as unpaired
            // single-end records: same qname, neither flagged read1 or read2,
            // hence one `(qname, mate)` key for two primaries. That was
            // unreachable while the rescue needed an already-mapped mate and
            // appeared the moment phase 2 let both mates through — 418
            // duplicate primaries on the simulated corpus, every one of them
            // rescue-produced.
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
        let mut w = std::io::BufWriter::with_capacity(
            8 << 20,
            std::fs::OpenOptions::new()
                .append(true)
                .open(&baln_path)?,
        );
        for (name, seq, qual, a, _, _) in &rescued {
            let mapq = esperanto_map::mapq::mapq_of(a);
            let mut a2 = a.clone();
            a2.rescued = false; // the RE tag below carries the provenance
            let mut br = esperanto_map::bam::record_se(name.as_bytes(), seq, qual, Some(&a2));
            // Same value the BAM record carries: `a2.rescued` was cleared to
            // move the provenance into the RE tag, which also removes the cap
            // `mapq_of` would apply, so take the number before clearing it.
            br.mapq = mapq;
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
    Ok(rescued_n)
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

/// Primary placements Tier 1 does not stand behind, offered to the rescue
/// passes in the read's own orientation and named the way `unmapped.fq.gz`
/// names its records, so both sources share one code path downstream.
///
/// `unmapped.fq.gz` alone is the wrong input pool for this channel, and not by
/// a little. A local Smith-Waterman always finds a positive-scoring window in
/// 3 Gb, so esperanto maps 98.9% of reads and its unmapped set is 1.1% of the
/// library — the residue too degenerate to place at all. Porath, Carmi &
/// Levanon (Nat Commun 5:4726, 2014) run the same realignment on the unmapped
/// output of a *refusing* aligner, where the pool is 10-20% of the library and
/// enriched for hyperediting precisely because the mismatch ceiling threw
/// those reads out. Downstream of a permissive aligner the pool is instead
/// anti-selected for what the channel looks for: the hyperedited reads were
/// already placed, badly, and never arrive. That is measurable in the
/// rejection attribution — 28.5% of real rejections are `no_edits`, reads with
/// no editing signal at all at the locus offered.
///
/// So the pool is the two ways Tier 1 declines to explain a read:
///
///   - it aligned only part of the read, leaving more than
///     `1 - RESCUE_CLIPPED_COVERAGE` soft-clipped; the clipped remainder is
///     where a run of edits would sit.
///   - it aligned the whole read but with a mismatch spectrum at or above
///     [`esperanto_map::mapq::HYPEREDIT_MIN_DENSITY`] that A-to-I editing does
///     not account for. That is the same threshold, read off the same `EK`
///     field, that already capped the placement's MAPQ below the pileup gate,
///     so these records contribute no call as they stand and a verified
///     editing placement can only improve on them.
///
/// `EK` reports `mm` (ordinary mismatches) but sums the two editing
/// directions, so the density here is `mm` over aligned columns — a lower
/// bound on `ReadAlignment::unexplained_density`, which also charges the minor
/// editing direction. The pool is therefore slightly under-inclusive relative
/// to what the MAPQ cap distrusted, which is the safe direction.
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
        // Aligned columns, matching `ReadAlignment::aligned_and_matching` so
        // the density below is the same quantity the MAPQ cap used.
        let aligned: u32 = rec
            .cigar()
            .iter()
            .filter_map(|op| match op {
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Ins(n)
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
        // BAM stores reference orientation; `apply_t13` is an involution, so
        // applying it again recovers the read as sequenced — which is what the
        // rescue transform and `verify_rescue_placement` both expect.
        let qual: Vec<u8> = rec.qual().iter().map(|q| q.saturating_add(33)).collect();
        let (seq, qual) = esperanto_bamio::apply_t13(rec.is_reverse(), &rec.seq().as_bytes(), &qual);
        let name = String::from_utf8_lossy(rec.qname()).into_owned();
        let name = if rec.is_paired() {
            format!("{name}/{}", if rec.is_first_in_template() { 1 } else { 2 })
        } else {
            name
        };
        out.push((name, seq, qual));
    }
    Ok(out)
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
        qs.iter().map(|&q| (q.saturating_sub(33)) as f64).sum::<f64>() / len as f64
    } else {
        kept.iter().map(|&q| (q.saturating_sub(33)) as f64).sum::<f64>() / kept.len() as f64
    };
    mean >= RESCUE_MIN_MEAN_PHRED
}

/// Post-placement verification against the true reference. A placement
/// survives only when its mismatch pattern is that of a hyperedited read of
/// the pass's own direction: a dense cluster of high-quality edits of *that
/// one* direction dominating every other mismatch class, with sane cluster
/// geometry.
///
/// `ag` is the direction the read showed *as sequenced*. The comparison here
/// runs in reference orientation, and `apply_t13` reverse-complements a
/// minus-strand placement to get there, which relabels the direction: an
/// A>G-edited read placed on the minus strand presents as ref T/read C. The
/// expected class is therefore `ag XOR minus`, not `ag`. Taking it from `ag`
/// alone charged every edit of a minus-strand placement as an ordinary
/// mismatch, so no such placement could ever clear
/// [`RESCUE_MAX_NON_EDIT_MM`] — latent while the pass table produced only
/// plus-strand placements, live now that it produces both (see
/// [`rescue_dirs`]).
///
/// The opposite direction is not a second editing class here, it is an
/// ordinary mismatch. Counting both would restore the two-letter collapse at
/// the verification step and let a paralog that exploits transitions in both
/// directions pass as an edited transcript.
///
/// `Err` carries the gate that refused the placement, so the stage can report
/// *why* it rejected rather than only how many. A channel that rejects
/// everything is indistinguishable from a channel that is switched off unless
/// it says which gate closed.
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
    // The collapsed realignment is gapless by construction, not by refusal:
    // [`RESCUE_GAP_PENALTY`] prices a gap above what the whole read can earn,
    // so a local extension soft-clips instead of opening one. This arm is the
    // assertion that it did, kept as a returned gate rather than a panic so
    // that a gap arriving from somewhere else — a future splice or chain path
    // that does not consult the extension's penalties — shows up in the
    // rejection attribution instead of aborting the run. It is expected to be
    // zero; the walk below relies on it.
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
                        let editing =
                            fb == Base::from_ascii(expect_from) && rb == Base::from_ascii(expect_to);
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

/// Paired-end concordance for a rescued placement: the mate must sit on the
/// same contig within the insert window, in the opposite orientation.
///
/// `placed` must already include the placements the rescue itself produced,
/// not only the four-letter mapped ones. Both mates of a hyperedited fragment
/// are typically refused by Tier 1 together, so a gate that demands an
/// already-mapped mate rejects exactly the fragments this channel exists to
/// recover — measured on the simulated corpus, 2761 of 7987 placements passed
/// every editing-signature gate and died here alone. Each mate's placement is
/// evidence for the other, so the check runs once after all passes.
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
    let mut cp = CallParams::standard(bam.to_path_buf(), out.clone(), params.fasta.clone());
    cp.snp = params.snp.clone();
    cp.threads = params.threads;
    let stats = esperanto_scan::run_call(&cp).map_err(stage_err("scan"))?;
    eprintln!(
        "[scan] {} candidates ({} below threshold, {} known SNPs removed)",
        stats.candidates, stats.below_threshold, stats.snp_removed
    );
    let [e, c1, c2] = stats.channels;
    let ctrl_mean = (c1.sites + c2.sites) as f64 / 2.0;
    let (snr, fdr) = crate::stages::channel_snr_fdr(e.sites, ctrl_mean);
    eprintln!("[scan] channel S/N {snr:.2}x (implied FDR {:.1}%)", fdr * 100.0);
    Ok(out)
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

/// Supporting-read re-alignment (verify): withdraw votes whose exhaustive
/// re-alignment prefers another placement. Writes `candidates.verified.bed`
/// (input columns + alt_verified + withdrawn); sites are never dropped
/// here, only annotated — the VCF stage turns low-verification sites into
/// `FILTER=REALIGN` rows.
fn stage_verify(
    params: &RunParams,
    bam: &Path,
    bed: &Path,
    map_handles: Option<(
        &'static esperanto_map::index::Index,
        Option<Arc<esperanto_engine::L1Index>>,
    )>,
) -> Result<PathBuf, FlowError> {
    eprintln!("[verify] running");
    let out = bed.with_extension("verified.bed");
    let Some(index_path) = &params.index else {
        // Honest refusal: without the alignment index there is nothing to
        // re-align against, and a pass-through verified bed would be read as
        // real evidence downstream.
        return Err(FlowError::Entry(
            "verify requires an alignment index; pass --index or run setup".into(),
        ));
    };
    let vp = esperanto_verify::VerifyParams {
        bam: bam.to_path_buf(),
        candidates: bed.to_path_buf(),
        index: index_path.clone(),
        l1_bundle: params.l1_bundle.clone(),
        out: out.clone(),
        threads: params.threads,
        ..Default::default()
    };
    let stats = match &map_handles {
        // FASTQ entry this process: reuse the map stage's handles.
        Some((index, l1)) => esperanto_verify::run_with_index(&vp, index, l1.clone()),
        // BAM entry (no map stage ran here): load our own.
        None => esperanto_verify::run(&vp),
    }
    .map_err(stage_err("verify"))?;
    eprintln!(
        "[verify] {} candidates, {} sites kept, {} votes withdrawn of {} reads",
        stats.candidates, stats.sites_kept, stats.votes_withdrawn, stats.reads_examined
    );
    verify_control_channels(params, bam, bed, &out, map_handles)?;
    Ok(out)
}

/// Verify-output rows passing the keep rule (alt_verified >= 2 and >= 1% of
/// depth), counting filter-clean candidates only.
fn verify_kept_count(path: &Path) -> Result<usize, FlowError> {
    let text = fs::read_to_string(path)?;
    let mut n = 0;
    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        let raw: Vec<&str> = line.split('\t').collect();
        if raw.len() < 12 || !raw[10].is_empty() {
            continue;
        }
        let f: Vec<&str> = raw.into_iter().filter(|s| !s.is_empty()).collect();
        let depth: f64 = f[4].parse().unwrap_or(0.0);
        let altv: f64 = f[f.len() - 2].parse().unwrap_or(0.0);
        if altv >= 2.0 && (depth == 0.0 || altv / depth >= 0.01) {
            n += 1;
        }
    }
    Ok(n)
}

/// Control channels through the same withdrawal machinery: C1/C2 candidates
/// are verified with their own alternates, and only the implied-FDR
/// bookkeeping consumes the results (plus the kept beds for audit).
fn verify_control_channels(
    params: &RunParams,
    bam: &Path,
    bed: &Path,
    verified: &Path,
    map_handles: Option<(
        &'static esperanto_map::index::Index,
        Option<Arc<esperanto_engine::L1Index>>,
    )>,
) -> Result<(), FlowError> {
    let dir = bed.parent().unwrap_or_else(|| Path::new(".")).to_path_buf();
    let control_bed = dir.join("control_candidates.bed");
    let channels_path = dir.join("channels.json");
    if !control_bed.exists() || !channels_path.exists() {
        return Ok(());
    }
    let Some(index_path) = &params.index else {
        return Ok(());
    };
    let text = fs::read_to_string(&control_bed)?;
    let mut inputs = Vec::new();
    for (name, tag, ch) in [
        ("control_c1", "C1", esperanto_verify::Channel::Ctrl1),
        ("control_c2", "C2", esperanto_verify::Channel::Ctrl2),
    ] {
        let mut body = String::new();
        for line in text.lines() {
            if line.ends_with(tag) {
                body.push_str(line);
                body.push('\n');
            }
        }
        let inp = dir.join(format!("{name}.bed"));
        fs::write(&inp, body)?;
        inputs.push((inp, ch));
    }
    let mut kept = [verify_kept_count(verified)?, 0, 0];
    for (i, (inp, ch)) in inputs.iter().enumerate() {
        let vout = inp.with_extension("verified.bed");
        let vp = esperanto_verify::VerifyParams {
            bam: bam.to_path_buf(),
            candidates: inp.clone(),
            index: index_path.clone(),
            l1_bundle: params.l1_bundle.clone(),
            out: vout.clone(),
            threads: params.threads,
            channel: *ch,
            ..Default::default()
        };
        match &map_handles {
            Some((index, l1)) => esperanto_verify::run_with_index(&vp, index, l1.clone()),
            None => esperanto_verify::run(&vp),
        }
        .map_err(stage_err("verify"))?;
        kept[i + 1] = verify_kept_count(&vout)?;
    }
    let mut doc: serde_json::Value = serde_json::from_str(&fs::read_to_string(&channels_path)?)
        .map_err(|e| FlowError::Json(e.to_string()))?;
    let pre_fdr = doc["implied_fdr"].as_f64().unwrap_or(0.0);
    let post_mean = (kept[1] + kept[2]) as f64 / 2.0;
    let post_fdr = if kept[0] > 0 {
        post_mean / kept[0] as f64
    } else {
        0.0
    };
    doc["editing_post"] = serde_json::json!(kept[0]);
    doc["control_1_post"] = serde_json::json!(kept[1]);
    doc["control_2_post"] = serde_json::json!(kept[2]);
    doc["implied_fdr_post"] = serde_json::json!(post_fdr);
    fs::write(&channels_path, serde_json::to_string_pretty(&doc).map_err(|e| FlowError::Json(e.to_string()))?)?;
    let c1_pre = doc["control_1"]["sites"].as_u64().unwrap_or(0);
    let c2_pre = doc["control_2"]["sites"].as_u64().unwrap_or(0);
    eprintln!(
        "[verify] control channels: C1 {} -> {} sites, C2 {} -> {} sites; implied FDR pre {:.1}% -> post {:.1}%",
        c1_pre,
        kept[1],
        c2_pre,
        kept[2],
        pre_fdr * 100.0,
        post_fdr * 100.0
    );
    Ok(())
}

/// Hyperedited-region caller over the collapsed-rescue record. FASTQ entries
/// only: without map/rescued.tsv the stage is a no-op.
fn stage_hyper(params: &RunParams) -> Result<(), FlowError> {
    let map_dir = params.out_dir.join("map");
    if !map_dir.join("rescued.tsv").is_file() {
        return Ok(());
    }
    eprintln!("[hyper] running");
    let (nr, ns) = crate::hyper::run(&map_dir, &params.out_dir, &params.fasta)?;
    eprintln!("[hyper] {nr} regions, {ns} sites");
    Ok(())
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

    /// Same fixture as `hyper_read` but A>C edits (the control alphabet).
    fn control_read_ac(
        r: &esperanto_map::fasta::Reference,
        start: usize,
        lo: usize,
        hi: usize,
    ) -> Vec<u8> {
        let mut seq = r.contigs[0].slice_ascii(start as u32, (start + 150) as u32);
        for (qi, b) in seq.iter_mut().enumerate() {
            if (lo..=hi).contains(&qi) && (start + qi).is_multiple_of(7) && *b == b'A' {
                *b = b'C';
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
    fn verify_control_alphabet_accepts_ac_not_ag() {
        // Same hyperedited-read fixture, but the edits are A>C: the control
        // pair (A,C) must accept, the editing pair (A,G) must not.
        let r = reference_with_period(7);
        let seq = control_read_ac(r, 20, 20, 110);
        assert!(seq.iter().filter(|&b| *b == b'C').count() >= 8);
        let aln = alignment_at(20, 150);
        let qual = vec![b'I'; 150];
        assert!(verify_rescue_placement(&aln, &seq, &qual, r, (b'A', b'C')).is_ok());
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
        // A gap disqualifies the placement, and re-testing it ungapped is not
        // a way round that: measured on the simulated corpus, all 982
        // pure-indel placements that could be re-tested that way went on to
        // fail `non_edit_mm` (see the gate's own comment). The gap is real.
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
}