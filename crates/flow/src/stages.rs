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
fn load_scores(path: &Path) -> Result<Vec<f64>, FlowError> {
    let text = fs::read_to_string(path)?;
    let mut probs = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        let p = line
            .split('\t')
            .nth(2)
            .and_then(|v| v.parse::<f64>().ok())
            .ok_or_else(|| FlowError::BedParse {
                line: i + 1,
                msg: "scores.tsv expects 'chrom<TAB>pos<TAB>prob'".into(),
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
    rescue_collapsed(paidx, dir, threads, &raw_bam)?;
    Ok(raw_bam)
}

/// Collapsed-alphabet (A==G, T==C) rescue of unmapped reads. When a
/// `<index stem>.cpaidx` sits next to the alignment index, the unmapped set
/// is re-aligned against it; survivors are written back into raw.bam with
/// MAPQ 0 and an `RE:Z:collapsed` tag (repeat-family placement, never
/// confident), and unmapped.fq.gz is rewritten to the truly unplaced.
/// No-op when no collapsed index exists.
fn stage_rescue_collapsed(params: &RunParams, raw_bam: &Path) -> Result<(), FlowError> {
    let Some(paidx) = params.index.clone() else {
        return Ok(());
    };
    rescue_collapsed(&paidx, &params.out_dir.join("map"), params.threads, raw_bam)
}

fn rescue_collapsed(
    paidx: &Path,
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
            fq.par_iter()
                .map_init(
                    || {
                        let mut a = esperanto_map::align::Aligner::new(&cidx, cfg);
                        a.jlib = lib.clone();
                        a
                    },
                    |al, (_, seq, _)| {
                        let c = collapse(seq);
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

    let mut rescued_names: std::collections::HashSet<&str> = std::collections::HashSet::new();
    let mut rescued: Vec<(String, Vec<u8>, Vec<u8>, esperanto_map::mapq::ReadAlignment)> =
        Vec::new();
    for ((name, seq, qual), aln) in fq.iter().zip(results) {
        if let Some(mut a) = aln {
            a.second_chain_score = a.chain_score; // MAPQ -> 0 (repeat-family placement)
            a.rescued = true;
            rescued.push((name.clone(), seq.clone(), qual.clone(), a));
            rescued_names.insert(name);
        }
    }
if rescued.is_empty() {
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
    // Record the placement count into align_qc.json (the rescue runs after
    // the stats document is written, so the key is patched in here; the
    // report stage reads it, treating absence as 0).
    {
        let p = map_dir.join("align_qc.json");
        let mut stats: esperanto_map::stats::AlignStats =
            serde_json::from_str(&fs::read_to_string(&p)?).map_err(stage_err("map"))?;
        stats.rescued_collapsed = Some(rescued.len() as u64);
        let json = serde_json::to_string_pretty(&stats).map_err(stage_err("map"))?;
        fs::write(&p, format!("{json}\n"))?;
    }
    eprintln!(
        "[rescue] collapsed: {} of {} unmapped placed (MAPQ 0)",
        rescued.len(),
        fq.len()
    );
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
    let ask: Option<&(dyn Fn() -> bool + Send + Sync)> = params
        .device_ask
        .as_ref()
        .map(|a| a.fn_ref());
    let probs = score_pipeline::score_sites_batched(
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
