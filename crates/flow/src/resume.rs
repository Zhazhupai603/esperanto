//! Resume support (spec §resume): `run.json` params freeze, input
//! fingerprints, per-stage artifact validation, and the stage walk that
//! finds the first incomplete stage of an interrupted run.

use std::fs;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::params::{Entry, RunParams};
use crate::FlowError;

/// run.json schema version (bump on layout change).
pub const RUN_JSON_VERSION: u32 = 2;

/// BGZF end-of-file marker: presence proves a BAM/fq.gz writer finished.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// One input file frozen at run start (path + size + mtime seconds).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct InputFp {
    /// Input role: r1 / r2 / bam / sites.
    pub role: String,
    /// Path as given at run start.
    pub path: PathBuf,
    /// File size in bytes at run start.
    pub size: u64,
    /// Modification time (unix seconds) at run start.
    pub mtime: u64,
}

/// The parameter block frozen into run.json.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ParamsJson {
    /// paidx index path (FASTQ entries).
    pub index: Option<PathBuf>,
    /// Reference FASTA.
    pub fasta: PathBuf,
    /// Optional GTF.
    pub gtf: Option<PathBuf>,
    /// score bundle root.
    pub bundle: PathBuf,
    /// Optional L1 engine bundle.
    pub l1_bundle: Option<PathBuf>,
    /// Worker threads.
    pub threads: usize,
}

/// run.json: everything `esperanto resume` needs to reconstruct a run.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RunJson {
    /// Schema version (RUN_JSON_VERSION).
    pub version: u32,
    /// esperanto version that started the run.
    pub esperanto: String,
    /// Sample name (the run directory name).
    pub sample: String,
    /// Entry type: fastq-se | fastq-pe | bam | bam-sites.
    pub entry: String,
    /// Frozen input fingerprints.
    pub inputs: Vec<InputFp>,
    /// Frozen resolved parameters.
    pub params: ParamsJson,
}

fn fp(role: &str, path: &Path) -> Result<InputFp, FlowError> {
    let meta = fs::metadata(path)?;
    let mtime = meta
        .modified()?
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    Ok(InputFp {
        role: role.to_string(),
        path: path.to_path_buf(),
        size: meta.len(),
        mtime,
    })
}

fn entry_name(entry: Entry) -> &'static str {
    match entry {
        Entry::FastqSe => "fastq-se",
        Entry::FastqPe => "fastq-pe",
        Entry::Bam => "bam",
        Entry::BamSites => "bam-sites",
    }
}

/// Parse the entry tag written by [`entry_name`].
pub fn parse_entry(tag: &str) -> Result<Entry, FlowError> {
    match tag {
        "fastq-se" => Ok(Entry::FastqSe),
        "fastq-pe" => Ok(Entry::FastqPe),
        "bam" => Ok(Entry::Bam),
        "bam-sites" => Ok(Entry::BamSites),
        other => Err(FlowError::Entry(format!(
            "run.json holds unknown entry '{other}'"
        ))),
    }
}

/// Freeze the resolved run into `<dir>/run.json` (atomic tmp + rename).
pub fn write_run_json(dir: &Path, sample: &str, params: &RunParams) -> Result<(), FlowError> {
    let entry = params.entry()?;
    let mut inputs = Vec::new();
    for p in &params.r1 {
        inputs.push(fp("r1", p)?);
    }
    for p in &params.r2 {
        inputs.push(fp("r2", p)?);
    }
    if let Some(b) = &params.bam {
        inputs.push(fp("bam", b)?);
    }
    if let Some(s) = &params.sites {
        inputs.push(fp("sites", s)?);
    }
    let doc = RunJson {
        version: RUN_JSON_VERSION,
        esperanto: env!("CARGO_PKG_VERSION").to_string(),
        sample: sample.to_string(),
        entry: entry_name(entry).to_string(),
        inputs,
        params: ParamsJson {
            index: params.index.clone(),
            fasta: params.fasta.clone(),
            gtf: params.gtf.clone(),
            bundle: params.bundle.clone(),
            l1_bundle: params.l1_bundle.clone(),
            threads: params.threads,
        },
    };
    let text = serde_json::to_string_pretty(&doc)
        .map_err(|e| FlowError::Entry(format!("run.json serialize: {e}")))?;
    let tmp = dir.join(".run.json.tmp");
    fs::write(&tmp, format!("{text}\n"))?;
    fs::rename(&tmp, dir.join("run.json"))?;
    Ok(())
}

/// Read and version-check `<dir>/run.json`.
pub fn read_run_json(dir: &Path) -> Result<RunJson, FlowError> {
    let path = dir.join("run.json");
    let text = fs::read_to_string(&path).map_err(|_| {
        FlowError::Entry(format!(
            "not an ESPERANTO run directory: {} (no run.json); start with `esperanto run`",
            dir.display()
        ))
    })?;
    let doc: RunJson = serde_json::from_str(&text).map_err(|e| {
        FlowError::Entry(format!("run.json in {} does not parse: {e}", dir.display()))
    })?;
    if doc.version != RUN_JSON_VERSION {
        return Err(FlowError::Entry(format!(
            "run.json version {} in {} is not supported (requires {RUN_JSON_VERSION}); start a new run to regenerate unified candidates and scores",
            doc.version,
            dir.display()
        )));
    }
    Ok(doc)
}

/// Refuse to resume when an input changed since the run started.
pub fn verify_inputs(doc: &RunJson) -> Result<(), FlowError> {
    for i in &doc.inputs {
        let meta = fs::metadata(&i.path).map_err(|_| {
            FlowError::Entry(format!(
                "input {} is missing (role {}); resume refused",
                i.path.display(),
                i.role
            ))
        })?;
        let mtime = meta
            .modified()?
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_secs())
            .unwrap_or(0);
        if meta.len() != i.size || mtime != i.mtime {
            return Err(FlowError::Entry(format!(
                "input {} changed since the run started (size {} -> {}, mtime {} -> {}); \
                 resume refused — delete the run directory and start over",
                i.path.display(),
                i.size,
                meta.len(),
                i.mtime,
                mtime
            )));
        }
    }
    Ok(())
}

/// Refuse to resume when a frozen parameter path no longer exists.
pub fn verify_param_paths(doc: &RunJson) -> Result<(), FlowError> {
    let mut required: Vec<&PathBuf> = vec![&doc.params.fasta, &doc.params.bundle];
    if let Some(p) = &doc.params.index {
        required.push(p);
    }
    for p in required {
        if !p.exists() {
            return Err(FlowError::Entry(format!(
                "parameter path {} from run.json no longer exists; resume refused",
                p.display()
            )));
        }
    }
    Ok(())
}

/// Pipeline stages in execution order (Ord matches the cascade).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum Stage {
    /// qc (FASTQ entries).
    Qc,
    /// map (alignment; includes the collapsed rescue).
    Map,
    /// Collapsed rescue only (map seal: alignment valid, rescue pending).
    Rescue,
    /// coordinate sort + BAI.
    Sort,
    /// scan + candidate filter (one artifact, one stage).
    Scan,
    /// score.
    Score,
    /// vcf.
    Vcf,
    /// HTML report.
    Report,
}

impl Stage {
    /// Parse a `--stop-after` stage name (qc|map|sort|scan|score|vcf).
    pub fn parse(s: &str) -> Option<Stage> {
        Some(match s {
            "qc" => Stage::Qc,
            "map" => Stage::Map,
            "sort" => Stage::Sort,
            "scan" => Stage::Scan,
            "score" => Stage::Score,
            "vcf" => Stage::Vcf,
            _ => return None,
        })
    }

    /// First stage executed for an entry type.
    pub fn first_for(entry: Entry) -> Stage {
        match entry {
            Entry::FastqSe | Entry::FastqPe => Stage::Qc,
            Entry::Bam => Stage::Scan,
            Entry::BamSites => Stage::Score,
        }
    }

    /// Stage tag for stderr notes.
    pub fn name(self) -> &'static str {
        match self {
            Stage::Qc => "qc",
            Stage::Map => "map",
            Stage::Rescue => "rescue",
            Stage::Sort => "sort",
            Stage::Scan => "scan",
            Stage::Score => "score",
            Stage::Vcf => "vcf",
            Stage::Report => "report",
        }
    }
}

/// File ends with the BGZF EOF marker (writer finished cleanly).
pub fn valid_bgzf(path: &Path) -> bool {
    use std::io::{Read, Seek, SeekFrom};
    let Ok(mut f) = fs::File::open(path) else {
        return false;
    };
    let Ok(len) = f.metadata().map(|m| m.len()) else {
        return false;
    };
    if len < BGZF_EOF.len() as u64 {
        return false;
    }
    if f.seek(SeekFrom::End(-(BGZF_EOF.len() as i64))).is_err() {
        return false;
    }
    let mut tail = [0u8; 28];
    f.read_exact(&mut tail).is_ok() && tail == BGZF_EOF
}

/// Plain-gzip streams (qc clean reads, unmapped.fq.gz): integrity is the
/// gzip trailer check — the stream must decompress end-to-end. Multi-member
/// aware (qc writes one member per chunk).
pub fn valid_gzip(path: &Path) -> bool {
    let Ok(f) = fs::File::open(path) else {
        return false;
    };
    let mut dec = flate2::read::MultiGzDecoder::new(f);
    let mut sink = std::io::sink();
    std::io::copy(&mut dec, &mut sink).is_ok()
}

/// File parses as JSON.
pub fn valid_json(path: &Path) -> bool {
    fs::read_to_string(path)
        .map(|t| serde_json::from_str::<serde_json::Value>(&t).is_ok())
        .unwrap_or(false)
}

/// File exists and is non-empty.
pub fn valid_nonempty(path: &Path) -> bool {
    fs::metadata(path).map(|m| m.len() > 0).unwrap_or(false)
}

/// candidates.bed parses end-to-end through the sites bridge.
fn valid_bed(path: &Path) -> bool {
    crate::stages::bed_to_sites(path).is_ok()
}

/// scores.tsv: empty is legal (zero candidates); otherwise every data line
/// is `chrom<TAB>pos<TAB>score` with a parseable classifier score.
fn valid_scores(path: &Path) -> bool {
    let Ok(text) = fs::read_to_string(path) else {
        return false;
    };
    text.lines().filter(|l| !l.is_empty()).all(|l| {
        match l.split('\t').collect::<Vec<_>>().as_slice() {
            [_, pos, score] => {
                pos.parse::<i64>().is_ok() && score.parse::<f64>().is_ok_and(|v| v.is_finite())
            }
            _ => false,
        }
    })
}

/// sites.vcf: starts with the VCF format header.
fn valid_vcf(path: &Path) -> bool {
    fs::read_to_string(path)
        .map(|t| t.starts_with("##fileformat=VCF"))
        .unwrap_or(false)
}

/// Require the persisted rescued count, rescued.tsv rows, and
/// RE:Z:collapsed records in raw.bam to agree. An incomplete merge or sidecar
/// write restarts rescue before later stages consume the BAM.
fn rescue_recorded(map_dir: &Path) -> bool {
    let sealed = fs::read_to_string(map_dir.join("align_qc.json"))
        .ok()
        .and_then(|t| serde_json::from_str::<serde_json::Value>(&t).ok())
        .and_then(|v| v["rescued_collapsed"].as_u64());
    match sealed {
        Some(n) => rescue_seal_consistent(map_dir, n),
        None => {
            // rescued.tsv is written before the seal key, so the sidecar
            // without the key means the rescue was interrupted mid-write and
            // must rerun — even when unmapped.fq.gz already lost the reads.
            if map_dir.join("rescued.tsv").is_file() {
                return false;
            }
            // An empty unmapped set means the rescue was a no-op by
            // definition (the key is only written when there is something to
            // re-align).
            let Ok(f) = fs::File::open(map_dir.join("unmapped.fq.gz")) else {
                return false;
            };
            use std::io::Read as _;
            let mut dec = flate2::read::GzDecoder::new(f);
            let mut byte = [0u8; 1];
            matches!(dec.read(&mut byte), Ok(0))
        }
    }
}

/// The sealed rescued count must equal the rescued.tsv data rows and the
/// RE:Z:collapsed records present in raw.bam. All three counters are checked
/// even when the seal is zero: a pass killed between a partial merge and the
/// seal write can leave RE records behind with no key at all, and only a
/// clean triple zero — key, sidecar and BAM — counts as recorded.
fn rescue_seal_consistent(map_dir: &Path, sealed: u64) -> bool {
    let tsv_rows = match fs::read_to_string(map_dir.join("rescued.tsv")) {
        Ok(t) => t
            .lines()
            .filter(|l| !l.is_empty())
            .count()
            .saturating_sub(1) as u64,
        Err(_) => 0,
    };
    if tsv_rows != sealed {
        return false;
    }
    let Ok(mut rdr) = rust_htslib::bam::Reader::from_path(map_dir.join("raw.bam")) else {
        return false;
    };
    use rust_htslib::bam::Read as _;
    let mut re = 0u64;
    for r in rdr.records() {
        let Ok(rec) = r else {
            return false;
        };
        if matches!(
            rec.aux(b"RE"),
            Ok(rust_htslib::bam::record::Aux::String("collapsed"))
        ) {
            re += 1;
        }
    }
    re == sealed
}

/// Map-stage three-state check: done / only the rescue is pending / rerun.
/// The seal is raw.bam + align_qc.json + unmapped.fq.gz intact;
/// `has_rescue_index` decides whether the rescue is expected at all.
fn map_state(out_dir: &Path, has_rescue_index: bool) -> Stage {
    let map_dir = out_dir.join("map");
    let seal = valid_bgzf(&map_dir.join("raw.bam"))
        && valid_json(&map_dir.join("align_qc.json"))
        && valid_gzip(&map_dir.join("unmapped.fq.gz"));
    if !seal {
        return Stage::Map;
    }
    if has_rescue_index && !rescue_recorded(&map_dir) {
        return Stage::Rescue;
    }
    Stage::Sort
}

/// Exactly one clean file per expected suffix, each a finished gzip stream.
fn qc_clean_ok(qc_dir: &Path, suffixes: &[&str]) -> bool {
    suffixes.iter().all(|s| {
        qc_dir
            .read_dir()
            .map(|d| {
                let hits: Vec<PathBuf> = d
                    .filter_map(|e| e.ok())
                    .map(|e| e.path())
                    .filter(|p| {
                        p.file_name()
                            .is_some_and(|n| n.to_string_lossy().ends_with(s))
                    })
                    .collect();
                hits.len() == 1 && valid_gzip(&hits[0])
            })
            .unwrap_or(false)
    })
}

/// The report artifact: `<sample>.report.html`, non-empty. The sample name
/// is the directory name (report contract).
fn valid_report(out_dir: &Path) -> bool {
    let Some(sample) = out_dir.file_name() else {
        return false;
    };
    valid_nonempty(&out_dir.join(format!("{}.report.html", sample.to_string_lossy())))
}

/// Walk qc → report and return the first stage whose artifacts fail
/// validation; `Report` when everything is intact.
pub fn walk(out_dir: &Path, entry: Entry, has_rescue_index: bool, has_gtf: bool) -> Stage {
    let first = Stage::first_for(entry);
    if first <= Stage::Qc {
        let qc_dir = out_dir.join("qc");
        let suffixes: &[&str] = match entry {
            Entry::FastqSe => &[".clean.fq.gz"],
            Entry::FastqPe => &[".clean_R1.fq.gz", ".clean_R2.fq.gz"],
            _ => &[],
        };
        if !valid_json(&qc_dir.join("qc.json"))
            || !valid_nonempty(&qc_dir.join("qc.html"))
            || !qc_clean_ok(&qc_dir, suffixes)
        {
            return Stage::Qc;
        }
    }
    if first <= Stage::Map {
        let st = map_state(out_dir, has_rescue_index);
        if st != Stage::Sort {
            return st;
        }
    }
    if first <= Stage::Sort {
        let sorted = out_dir.join("map").join("sorted.bam");
        if !valid_bgzf(&sorted)
            || !valid_nonempty(&PathBuf::from(format!("{}.bai", sorted.display())))
        {
            return Stage::Sort;
        }
    }
    if first <= Stage::Scan && !valid_bed(&out_dir.join("scan").join("candidates.bed")) {
        return Stage::Scan;
    }
    if first <= Stage::Score && !valid_scores(&out_dir.join("score").join("scores.tsv")) {
        return Stage::Score;
    }
    if first <= Stage::Vcf && !valid_vcf(&out_dir.join("sites.vcf")) {
        return Stage::Vcf;
    }
    if has_gtf && !valid_report(out_dir) {
        return Stage::Report;
    }
    Stage::Report
}

/// Delete a stage's artifacts so its re-run starts from a clean slate;
/// later stages cascade (their inputs change). Earlier stages stay.
pub fn clean_stage(out_dir: &Path, stage: Stage) {
    let map_dir = out_dir.join("map");
    let sorted = map_dir.join("sorted.bam");
    if stage <= Stage::Qc {
        let _ = fs::remove_dir_all(out_dir.join("qc"));
    }
    if stage <= Stage::Map {
        let _ = fs::remove_dir_all(&map_dir);
    }
    if stage == Stage::Rescue {
        let _ = fs::remove_file(map_dir.join("raw.merged.bam"));
    }
    if stage <= Stage::Sort && stage != Stage::Map {
        let _ = fs::remove_file(&sorted);
        let _ = fs::remove_file(format!("{}.bai", sorted.display()));
    }
    if stage <= Stage::Scan {
        let _ = fs::remove_dir_all(out_dir.join("scan"));
    }
    if stage <= Stage::Score {
        let _ = fs::remove_dir_all(out_dir.join("score"));
    }
    if stage <= Stage::Vcf {
        let _ = fs::remove_file(out_dir.join("sites.vcf"));
    }
    if stage <= Stage::Report {
        if let Ok(d) = fs::read_dir(out_dir) {
            for e in d.filter_map(|e| e.ok()) {
                if e.file_name().to_string_lossy().ends_with(".report.html") {
                    let _ = fs::remove_file(e.path());
                }
            }
        }
    }
}

/// Advisory lock against two processes writing the same run directory.
pub struct RunLock {
    path: PathBuf,
}

impl RunLock {
    /// Create `<dir>/.lock` exclusively; fails while another holder exists.
    pub fn acquire(dir: &Path) -> Result<RunLock, FlowError> {
        let path = dir.join(".lock");
        fs::OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(&path)
            .map_err(|_| {
                FlowError::Entry(format!(
                    "run directory {} is locked (another process, or a stale .lock from a killed run); \
                     remove {} to proceed",
                    dir.display(),
                    path.display()
                ))
            })?;
        Ok(RunLock { path })
    }
}

impl Drop for RunLock {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    fn tmp(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!(
            "esperanto-resume-test-{tag}-{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        d
    }

    fn bgzf(path: &Path) {
        let mut f = fs::File::create(path).unwrap();
        f.write_all(b"payload").unwrap();
        f.write_all(&BGZF_EOF).unwrap();
    }

    fn gzip(path: &Path, text: &str) {
        let f = fs::File::create(path).unwrap();
        let mut e = flate2::write::GzEncoder::new(f, flate2::Compression::new(1));
        e.write_all(text.as_bytes()).unwrap();
        e.finish().unwrap();
    }

    /// A fully valid PE run directory (minus inputs; walk does not read them).
    fn full_pe_dir(root: &Path) {
        let qc = root.join("qc");
        fs::create_dir_all(&qc).unwrap();
        fs::write(qc.join("qc.json"), "{}").unwrap();
        fs::write(qc.join("qc.html"), "<html>").unwrap();
        gzip(&qc.join("r1.clean_R1.fq.gz"), "@a/1\nACGT\n+\nIIII\n");
        gzip(&qc.join("r2.clean_R2.fq.gz"), "@a/2\nTGCA\n+\nIIII\n");
        let map = root.join("map");
        fs::create_dir_all(&map).unwrap();
        // raw.bam must parse as BAM: the rescue seal check scans it for
        // RE:Z:collapsed records even when the sealed count is zero.
        bam_with_records(&map.join("raw.bam"), 0, 0);
        fs::write(map.join("align_qc.json"), "{\"rescued_collapsed\": 0}").unwrap();
        gzip(&map.join("unmapped.fq.gz"), "");
        bgzf(&map.join("sorted.bam"));
        fs::write(map.join("sorted.bam.bai"), "x").unwrap();
        let scan = root.join("scan");
        fs::create_dir_all(&scan).unwrap();
        fs::write(
            scan.join("candidates.bed"),
            "chr1\t0\t100\t+\tx\ty\t10\tz\n",
        )
        .unwrap();
        fs::write(
            scan.join("candidates.verified.bed"),
            "chr1\t0\t100\tA\t10\t0\t0\t2\t0\t0.2\t\t2\t0\n",
        )
        .unwrap();
        let score = root.join("score");
        fs::create_dir_all(&score).unwrap();
        fs::write(score.join("scores.tsv"), "chr1\t100\t0.9\n").unwrap();
        fs::write(root.join("sites.vcf"), "##fileformat=VCFv4.2\n").unwrap();
        fs::write(
            root.join(format!(
                "{}.report.html",
                root.file_name().unwrap().to_string_lossy()
            )),
            "<html>",
        )
        .unwrap();
    }

    #[test]
    fn validators_accept_and_reject() {
        let d = tmp("validators");
        let good = d.join("good.bam");
        bgzf(&good);
        assert!(valid_bgzf(&good));
        let mut bad = good.clone();
        bad.set_file_name("bad.bam");
        fs::write(&bad, b"payload-without-marker").unwrap();
        assert!(!valid_bgzf(&bad));
        let g = d.join("g.gz");
        gzip(&g, "hello");
        assert!(valid_gzip(&g));
        let raw = fs::read(&g).unwrap();
        fs::write(d.join("t.gz"), &raw[..raw.len() - 4]).unwrap();
        assert!(!valid_gzip(&d.join("t.gz")));
        let _ = fs::remove_dir_all(&d);
    }

    #[test]
    fn walk_intact_dir_reports_done() {
        let d = tmp("walk-intact");
        full_pe_dir(&d);
        assert_eq!(walk(&d, Entry::FastqPe, true, true), Stage::Report);
        let _ = fs::remove_dir_all(&d);
    }

    #[test]
    fn walk_finds_first_broken_stage() {
        let d = tmp("walk-stages");
        full_pe_dir(&d);
        // Broken BAI -> Sort.
        fs::write(d.join("map").join("sorted.bam.bai"), "").unwrap();
        assert_eq!(walk(&d, Entry::FastqPe, true, true), Stage::Sort);
        // Broken unmapped -> Map (seal fails).
        full_pe_dir(&d);
        fs::write(d.join("map").join("unmapped.fq.gz"), b"not-gzip").unwrap();
        assert_eq!(walk(&d, Entry::FastqPe, true, true), Stage::Map);
        // Seal intact, rescue key missing, unmapped non-empty -> Rescue.
        full_pe_dir(&d);
        fs::write(d.join("map").join("align_qc.json"), "{}").unwrap();
        gzip(
            &d.join("map").join("unmapped.fq.gz"),
            "@a/1\nACGT\n+\nIIII\n",
        );
        assert_eq!(walk(&d, Entry::FastqPe, true, true), Stage::Rescue);
        // Same but empty unmapped: rescue counts as done and the rest of the
        // directory is intact, so the walk reaches the end.
        gzip(&d.join("map").join("unmapped.fq.gz"), "");
        assert_eq!(walk(&d, Entry::FastqPe, true, true), Stage::Report);
        // Missing vcf -> Vcf.
        full_pe_dir(&d);
        fs::remove_file(d.join("sites.vcf")).unwrap();
        assert_eq!(walk(&d, Entry::FastqPe, true, true), Stage::Vcf);
        // Missing report with a GTF -> Report.
        full_pe_dir(&d);
        fs::remove_file(d.join(format!(
            "{}.report.html",
            d.file_name().unwrap().to_string_lossy()
        )))
        .unwrap();
        assert_eq!(walk(&d, Entry::FastqPe, true, true), Stage::Report);
        // No GTF -> report optional -> still Report (nothing to do).
        assert_eq!(walk(&d, Entry::FastqPe, true, false), Stage::Report);
        let _ = fs::remove_dir_all(&d);
    }

    #[test]
    fn clean_stage_cascades_but_keeps_earlier() {
        let d = tmp("clean");
        full_pe_dir(&d);
        clean_stage(&d, Stage::Scan);
        assert!(d.join("qc").exists());
        assert!(d.join("map").join("raw.bam").exists());
        assert!(!d.join("scan").exists());
        assert!(!d.join("score").exists());
        assert!(!d.join("sites.vcf").exists());
        assert!(!d
            .join(format!(
                "{}.report.html",
                d.file_name().unwrap().to_string_lossy()
            ))
            .exists());
        // Rescue keeps raw.bam but drops the merge leftover and cascades.
        full_pe_dir(&d);
        fs::write(d.join("map").join("raw.merged.bam"), b"x").unwrap();
        clean_stage(&d, Stage::Rescue);
        assert!(d.join("map").join("raw.bam").exists());
        assert!(!d.join("map").join("raw.merged.bam").exists());
        assert!(!d.join("map").join("sorted.bam").exists());
        let _ = fs::remove_dir_all(&d);
    }

    #[test]
    fn run_json_roundtrip_and_fingerprint() {
        let d = tmp("runjson");
        let input = d.join("reads.fq");
        fs::write(&input, b"@a\nACGT\n+\nIIII\n").unwrap();
        let params = RunParams {
            r1: vec![input.clone()],
            r2: Vec::new(),
            bam: None,
            sites: None,
            index: Some(d.join("ref.paidx")),
            fasta: d.join("ref.fa"),
            gtf: None,
            bundle: d.join("bundle"),
            l1_bundle: None,
            out_dir: d.clone(),
            threads: 4,
            stop_after: None,
        };
        write_run_json(&d, "s1", &params).unwrap();
        let doc = read_run_json(&d).unwrap();
        assert_eq!(doc.sample, "s1");
        assert_eq!(doc.entry, "fastq-se");
        assert_eq!(doc.inputs.len(), 1);
        verify_inputs(&doc).unwrap();
        assert_eq!(doc.version, RUN_JSON_VERSION);
        let mut old = serde_json::to_value(&doc).unwrap();
        old["version"] = serde_json::json!(1);
        fs::write(d.join("run.json"), serde_json::to_string(&old).unwrap()).unwrap();
        let error = read_run_json(&d).unwrap_err().to_string();
        assert!(error.contains("start a new run"));
        // Size change -> refusal.
        fs::write(&input, b"@a\nACGTACGT\n+\nIIIIIIII\n").unwrap();
        assert!(verify_inputs(&doc).is_err());
        let _ = fs::remove_dir_all(&d);
    }

    #[test]
    fn lock_is_exclusive_and_released() {
        let d = tmp("lock");
        let l1 = RunLock::acquire(&d).unwrap();
        assert!(RunLock::acquire(&d).is_err());
        drop(l1);
        assert!(RunLock::acquire(&d).is_ok());
        let _ = fs::remove_dir_all(&d);
    }

    /// A minimal valid BAM (BGZF EOF marker included) with `re` records
    /// carrying RE:Z:collapsed and `plain` records without the tag.
    fn bam_with_records(path: &Path, re: u64, plain: u64) {
        let mut hdr = rust_htslib::bam::Header::new();
        let mut sq = rust_htslib::bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1000);
        hdr.push_record(&sq);
        let mut w =
            rust_htslib::bam::Writer::from_path(path, &hdr, rust_htslib::bam::Format::Bam).unwrap();
        for i in 0..(re + plain) {
            let mut rec = rust_htslib::bam::Record::new();
            let qname = format!("r{i}");
            rec.set(qname.as_bytes(), None, b"ACGT", &[30u8; 4]);
            rec.set_tid(0);
            rec.set_pos(i as i64);
            if i < re {
                rec.push_aux(b"RE", rust_htslib::bam::record::Aux::String("collapsed"))
                    .unwrap();
            }
            w.write(&rec).unwrap();
        }
    }

    fn rescued_tsv(path: &Path, rows: u64) {
        let mut text = String::from("qname\tchrom\tpos0\tstrand\tdir\tn_edits\tedit_pos0\tmapq\n");
        for i in 0..rows {
            text.push_str(&format!("r{i}\tchr1\t{i}\t+\tAG\t2\t1,2\t20\n"));
        }
        fs::write(path, text).unwrap();
    }

    /// Seal key, rescued.tsv rows and raw.bam RE:Z:collapsed records agree:
    /// the rescue output counts as recorded.
    #[test]
    fn rescue_recorded_accepts_consistent_seal() {
        let d = tmp("rescue-consistent");
        let map = d.join("map");
        fs::create_dir_all(&map).unwrap();
        fs::write(map.join("align_qc.json"), "{\"rescued_collapsed\": 2}").unwrap();
        rescued_tsv(&map.join("rescued.tsv"), 2);
        bam_with_records(&map.join("raw.bam"), 2, 3);
        assert!(rescue_recorded(&map));
        // A zero seal is accepted only as a clean triple zero: no sidecar
        // rows and no RE records in raw.bam.
        fs::write(map.join("align_qc.json"), "{\"rescued_collapsed\": 0}").unwrap();
        fs::remove_file(map.join("rescued.tsv")).unwrap();
        bam_with_records(&map.join("raw.bam"), 0, 3);
        assert!(rescue_recorded(&map));
        let _ = fs::remove_dir_all(&d);
    }

    /// A zero seal with RE residue in raw.bam (partial merge before the seal
    /// write, then killed) is inconsistent and must rerun the rescue.
    #[test]
    fn rescue_recorded_rejects_zero_seal_with_re_residue() {
        let d = tmp("rescue-zero-residue");
        let map = d.join("map");
        fs::create_dir_all(&map).unwrap();
        fs::write(map.join("align_qc.json"), "{\"rescued_collapsed\": 0}").unwrap();
        bam_with_records(&map.join("raw.bam"), 1, 3);
        assert!(!rescue_recorded(&map));
        // Same with the sidecar left behind: still rejected (row count
        // disagrees with the seal even before the BAM check).
        rescued_tsv(&map.join("rescued.tsv"), 1);
        assert!(!rescue_recorded(&map));
        let _ = fs::remove_dir_all(&d);
    }

    /// The defect window: seal key written, raw.bam never received the rescued
    /// records. The key alone must not count as recorded.
    #[test]
    fn rescue_recorded_rejects_keyed_but_unmerged() {
        let d = tmp("rescue-unmerged");
        let map = d.join("map");
        fs::create_dir_all(&map).unwrap();
        fs::write(map.join("align_qc.json"), "{\"rescued_collapsed\": 2}").unwrap();
        rescued_tsv(&map.join("rescued.tsv"), 2);
        bam_with_records(&map.join("raw.bam"), 0, 5);
        assert!(!rescue_recorded(&map));
        // Key value disagreeing with the sidecar alone is also a mismatch.
        fs::write(map.join("align_qc.json"), "{\"rescued_collapsed\": 5}").unwrap();
        bam_with_records(&map.join("raw.bam"), 2, 5);
        assert!(!rescue_recorded(&map));
        // A zero seal with leftover sidecar rows is inconsistent.
        fs::write(map.join("align_qc.json"), "{\"rescued_collapsed\": 0}").unwrap();
        assert!(!rescue_recorded(&map));
        let _ = fs::remove_dir_all(&d);
    }

    /// rescued.tsv without the seal key: the rescue was interrupted mid-write
    /// and must rerun even when unmapped.fq.gz is already empty.
    #[test]
    fn rescue_recorded_rejects_interrupted_before_seal() {
        let d = tmp("rescue-interrupted");
        let map = d.join("map");
        fs::create_dir_all(&map).unwrap();
        fs::write(map.join("align_qc.json"), "{}").unwrap();
        rescued_tsv(&map.join("rescued.tsv"), 2);
        gzip(&map.join("unmapped.fq.gz"), "");
        assert!(!rescue_recorded(&map));
        let _ = fs::remove_dir_all(&d);
    }

    /// End-to-end at the resume decision: a run directory killed inside the
    /// defect window (key written, raw.bam unmerged) walks back to the rescue
    /// stage; once the rescue output is made consistent, the same directory
    /// walks to done. A consistently sealed directory never revisits rescue.
    #[test]
    fn walk_resumes_rescue_on_window_state() {
        let d = tmp("walk-window");
        full_pe_dir(&d);
        let map = d.join("map");
        // Hand-built kill-in-window state: sealed 2, sidecar 2 rows, raw.bam
        // holds no rescued record.
        fs::write(map.join("align_qc.json"), "{\"rescued_collapsed\": 2}").unwrap();
        rescued_tsv(&map.join("rescued.tsv"), 2);
        bam_with_records(&map.join("raw.bam"), 0, 5);
        assert_eq!(walk(&d, Entry::FastqPe, true, true), Stage::Rescue);
        // After the rescue re-run completes, the three counters agree and the
        // walk proceeds past the map stage again.
        bam_with_records(&map.join("raw.bam"), 2, 5);
        assert_eq!(walk(&d, Entry::FastqPe, true, true), Stage::Report);
        let _ = fs::remove_dir_all(&d);
    }

    /// A directory sealed by a run whose rescue produced records and finished
    /// is not sent back to the rescue stage.
    #[test]
    fn walk_leaves_consistent_rescued_dir_alone() {
        let d = tmp("walk-rescued-intact");
        full_pe_dir(&d);
        let map = d.join("map");
        fs::write(map.join("align_qc.json"), "{\"rescued_collapsed\": 2}").unwrap();
        rescued_tsv(&map.join("rescued.tsv"), 2);
        bam_with_records(&map.join("raw.bam"), 2, 5);
        assert_eq!(walk(&d, Entry::FastqPe, true, true), Stage::Report);
        let _ = fs::remove_dir_all(&d);
    }
}
