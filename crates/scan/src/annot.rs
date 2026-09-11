//! Annotation indexes: --gtf gene annotation (strand-call evidence priority 3) and
//! --gnomad frequencies (soft down-weighting). Both are optional; missing = empty
//! index/None, no error. plain / .gz both accepted.

use crate::error::CallError;
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Read as _, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, RwLock};
fn read_text(path: &Path) -> Result<String, CallError> {
    let p = path.display().to_string();
    let f = std::fs::File::open(path).map_err(|e| CallError::Io {
        path: p.clone(),
        source: e,
    })?;
    let mut s = String::new();
    if path.extension().is_some_and(|e| e == "gz") {
        MultiGzDecoder::new(f)
            .read_to_string(&mut s)
            .map_err(|e| CallError::Io { path: p, source: e })?;
    } else {
        std::io::Read::read_to_string(&mut std::io::BufReader::new(f), &mut s)
            .map_err(|e| CallError::Io { path: p, source: e })?;
    }
    Ok(s)
}

// ---------------- GTF ----------------

#[derive(Debug, Clone, Copy)]
pub struct Gene {
    pub start: i64, // 0-based half-open
    pub end: i64,
    pub plus: bool,
}

/// Per-chromosome gene intervals (sorted by start).
#[derive(Debug, Default)]
pub struct GtfIndex {
    by_chrom: HashMap<String, Vec<Gene>>,
}

impl GtfIndex {
    pub fn load(path: &Path) -> Result<GtfIndex, CallError> {
        let text = read_text(path)?;
        let mut by_chrom: HashMap<String, Vec<Gene>> = HashMap::new();
        for (ln, line) in text.lines().enumerate() {
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 9 || f[2] != "gene" {
                continue;
            }
            let parse = |s: &str, what: &str| -> Result<i64, CallError> {
                s.parse().map_err(|_| CallError::Annot {
                    path: path.display().to_string(),
                    msg: format!("line {}: bad {what}: {s}", ln + 1),
                })
            };
            // GTF 1-based closed interval → 0-based half-open
            let start = parse(f[3], "start")? - 1;
            let end = parse(f[4], "end")?;
            let plus = match f[6] {
                "+" => true,
                "-" => false,
                _ => continue, // genes without strand info are not indexed
            };
            by_chrom
                .entry(f[0].to_string())
                .or_default()
                .push(Gene { start, end, plus });
        }
        for v in by_chrom.values_mut() {
            v.sort_by_key(|g| g.start);
        }
        Ok(GtfIndex { by_chrom })
    }

    /// Site (0-based) in gene → (has plus-strand gene, has minus-strand gene).
    pub fn strands_at(&self, chrom: &str, pos0: i64) -> (bool, bool) {
        let mut hit = (false, false);
        if let Some(genes) = self.by_chrom.get(chrom) {
            for g in genes {
                if g.start > pos0 {
                    break; // sorted by start, early exit
                }
                if pos0 >= g.start && pos0 < g.end {
                    if g.plus {
                        hit.0 = true;
                    } else {
                        hit.1 = true;
                    }
                }
            }
        }
        hit
    }
}

// ---------------- gnomAD ----------------

/// filename template for per-chrom directory mode; `{chrom}` is replaced by the contig name.
const DIR_PATTERN: &str = "gnomad.joint.v4.1.sites.{chrom}.vcf.bgz";
/// `.afidx` binary cache magic (includes version).
const CACHE_MAGIC: u64 = 0x4541_4649_4458_3031; // "EAFIDX01"

/// Single-chromosome AF index: pos0 ascending + parallel AF array, binary search.
#[derive(Debug, Default)]
struct ChromAf {
    pos: Vec<u32>,
    af: Vec<f32>,
}

impl ChromAf {
    fn af_at(&self, pos0: i64) -> Option<f64> {
        if pos0 < 0 || pos0 > u32::MAX as i64 {
            return None;
        }
        let p = pos0 as u32;
        match self.pos.binary_search(&p) {
            Ok(i) => Some(self.af[i] as f64),
            Err(_) => None,
        }
    }
}

/// Directory-mode state: lazily loaded per contig; builds/reads the `.afidx` cache on first use.
#[derive(Debug)]
struct GnomadDir {
    dir: PathBuf,
    cache_dir: PathBuf,
    chroms: RwLock<HashMap<String, Arc<ChromAf>>>,
}

/// Soft down-weighting index; miss = None (neutral).
///
/// Accepts two input forms:
/// - Single-file VCF (plain/.gz): fully loaded into memory (existing behavior);
/// - Directory: `gnomad.joint.v4.1.sites.{chrom}.vcf.bgz` per-chrom lazy loading,
///   contig names chr1↔1 auto-mapped; on first use of a contig, builds an `.afidx`
///   in the cache directory (pos0:u32 + AF:f32 ascending binary), then direct reads + binary search.
///
/// Missing path → hard error (no longer silently skipped).
#[derive(Debug, Default)]
pub struct GnomadIndex {
    by_chrom: HashMap<String, HashMap<i64, f64>>,
    dir: Option<GnomadDir>,
}

/// Primary chromosome name after stripping the chr prefix: 1-22 / X / Y (chrM/MT have no gnomAD file; treated as neutral).
fn is_primary(contig: &str) -> bool {
    let c = contig.strip_prefix("chr").unwrap_or(contig);
    matches!(c, "X" | "Y")
        || (!c.is_empty() && c.len() <= 2 && c.bytes().all(|b| b.is_ascii_digit()))
}

/// Resolve the per-chrom VCF for a contig within the directory: try the original name first, then the chr-prefix mapping.
fn resolve_vcf(dir: &Path, contig: &str) -> Option<PathBuf> {
    let alt = match contig.strip_prefix("chr") {
        Some(s) => s.to_string(),
        None => format!("chr{contig}"),
    };
    [contig.to_string(), alt].iter().find_map(|name| {
        let p = dir.join(DIR_PATTERN.replace("{chrom}", name));
        p.is_file().then_some(p)
    })
}

/// Cache directory: $ESPERANTO_CACHE_DIR/gnomad, else ~/.cache/esperanto/gnomad.
fn gnomad_cache_dir() -> Result<PathBuf, CallError> {
    if let Some(d) = std::env::var_os("ESPERANTO_CACHE_DIR") {
        return Ok(PathBuf::from(d).join("gnomad"));
    }
    let home = std::env::var_os("HOME").ok_or_else(|| CallError::Annot {
        path: "<env>".into(),
        msg: "HOME is not set; cannot locate the gnomad index cache directory".into(),
    })?;
    Ok(PathBuf::from(home).join(".cache/esperanto/gnomad"))
}

fn annot_err(path: &Path, msg: impl Into<String>) -> CallError {
    CallError::Annot {
        path: path.display().to_string(),
        msg: msg.into(),
    }
}

/// per-chrom VCF → ChromAf (streaming parse; keeps only lines with a parseable AF; 0 matches = wrong contig name, hard error).
fn build_chrom(vcf: &Path, file_chrom: &str) -> Result<ChromAf, CallError> {
    let f = std::fs::File::open(vcf).map_err(|e| CallError::Io {
        path: vcf.display().to_string(),
        source: e,
    })?;
    let mut rdr = BufReader::with_capacity(1 << 20, MultiGzDecoder::new(f));
    let mut idx = ChromAf::default();
    let mut line = String::new();
    loop {
        line.clear();
        let n = rdr.read_line(&mut line).map_err(|e| CallError::Io {
            path: vcf.display().to_string(),
            source: e,
        })?;
        if n == 0 {
            break;
        }
        let l = line.trim_end();
        if l.starts_with('#') || l.is_empty() {
            continue;
        }
        // Only keep lines for the target contig (col0 == file_chrom).
        if !l.starts_with(file_chrom) || l.as_bytes().get(file_chrom.len()) != Some(&b'\t') {
            continue;
        }
        let rest = &l[file_chrom.len() + 1..];
        let t = rest
            .find('\t')
            .ok_or_else(|| annot_err(vcf, format!("bad line (no POS): {l:.80}")))?;
        let pos1: i64 = rest[..t]
            .parse()
            .map_err(|_| annot_err(vcf, format!("bad POS: {l:.80}")))?;
        // First AF value in INFO (same semantics as file mode: take the first comma-separated value).
        // gnomAD joint v4.1 uses the key AF_joint (plus AF_genomes/AF_exomes/subset AF_joint_XX etc.);
        // bare AF from older formats is also accepted. "\t" anchors the INFO-start boundary.
        let af = [";AF_joint=", "\tAF_joint=", ";AF=", "\tAF="]
            .iter()
            .find_map(|pat| l.find(pat).map(|i| &l[i + pat.len()..]))
            .and_then(|s| {
                let end = s.find([',', ';', '\t']).unwrap_or(s.len());
                s[..end].parse::<f64>().ok()
            });
        if let (Some(af), Ok(p0)) = (af, u32::try_from(pos1 - 1)) {
            idx.pos.push(p0);
            idx.af.push(af as f32);
        }
    }
    if idx.pos.is_empty() {
        return Err(annot_err(
            vcf,
            format!("0 lines match contig {file_chrom}: contig name inconsistent with file contents"),
        ));
    }
    // tabix files should already be coordinate-sorted; safety: sort if out of order.
    if !idx.pos.windows(2).all(|w| w[0] <= w[1]) {
        let mut order: Vec<usize> = (0..idx.pos.len()).collect();
        order.sort_by_key(|&i| idx.pos[i]);
        idx.pos = order.iter().map(|&i| idx.pos[i]).collect();
        idx.af = order.iter().map(|&i| idx.af[i]).collect();
    }
    Ok(idx)
}

/// Write the `.afidx` cache: magic + vcf size + vcf mtime + n + n×(pos0:u32, af:f32). tmp+rename is atomic.
fn write_cache(cache: &Path, vcf: &Path, idx: &ChromAf) -> Result<(), CallError> {
    let meta = std::fs::metadata(vcf).map_err(|e| CallError::Io {
        path: vcf.display().to_string(),
        source: e,
    })?;
    let mtime = meta
        .modified()
        .ok()
        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let tmp = cache.with_extension("afidx.tmp");
    let mut w = BufWriter::new(std::fs::File::create(&tmp).map_err(|e| CallError::Io {
        path: tmp.display().to_string(),
        source: e,
    })?);
    let mut wr = |b: &[u8]| {
        w.write_all(b).map_err(|e| CallError::Io {
            path: tmp.display().to_string(),
            source: e,
        })
    };
    wr(&CACHE_MAGIC.to_le_bytes())?;
    wr(&meta.len().to_le_bytes())?;
    wr(&mtime.to_le_bytes())?;
    wr(&(idx.pos.len() as u64).to_le_bytes())?;
    for (&p, &a) in idx.pos.iter().zip(&idx.af) {
        wr(&p.to_le_bytes())?;
        wr(&a.to_le_bytes())?;
    }
    w.flush().map_err(|e| CallError::Io {
        path: tmp.display().to_string(),
        source: e,
    })?;
    drop(w);
    std::fs::rename(&tmp, cache).map_err(|e| CallError::Io {
        path: cache.display().to_string(),
        source: e,
    })?;
    Ok(())
}

/// Read the `.afidx` cache; magic/size/mtime mismatch → None (triggers rebuild).
fn read_cache(cache: &Path, vcf: &Path) -> Option<ChromAf> {
    let buf = std::fs::read(cache).ok()?;
    let meta = std::fs::metadata(vcf).ok()?;
    let mtime = meta
        .modified()
        .ok()?
        .duration_since(std::time::UNIX_EPOCH)
        .ok()?
        .as_secs();
    let hdr = buf.get(..32)?;
    let magic = u64::from_le_bytes(hdr[0..8].try_into().ok()?);
    let vlen = u64::from_le_bytes(hdr[8..16].try_into().ok()?);
    let vmt = u64::from_le_bytes(hdr[16..24].try_into().ok()?);
    let n = u64::from_le_bytes(hdr[24..32].try_into().ok()?) as usize;
    if magic != CACHE_MAGIC || vlen != meta.len() || vmt != mtime {
        return None;
    }
    let body = buf.get(32..32 + n * 8)?;
    let mut idx = ChromAf {
        pos: Vec::with_capacity(n),
        af: Vec::with_capacity(n),
    };
    for rec in body.chunks_exact(8) {
        idx.pos.push(u32::from_le_bytes(rec[0..4].try_into().ok()?));
        idx.af.push(f32::from_le_bytes(rec[4..8].try_into().ok()?));
    }
    Some(idx)
}

/// Cache filename = VCF filename + .afidx; the contig name in the vcf filename is used for line filtering.
fn load_or_build(vcf: &Path, cache_dir: &Path) -> Result<ChromAf, CallError> {
    let fname = vcf.file_name().and_then(|s| s.to_str()).unwrap_or_default();
    let cache = cache_dir.join(format!("{fname}.afidx"));
    if let Some(idx) = read_cache(&cache, vcf) {
        eprintln!("[gnomad] index cache hit {}", cache.display());
        return Ok(idx);
    }
    let file_chrom = fname
        .strip_prefix("gnomad.joint.v4.1.sites.")
        .and_then(|s| s.strip_suffix(".vcf.bgz"))
        .ok_or_else(|| annot_err(vcf, format!("filename does not match the {DIR_PATTERN} pattern")))?;
    eprintln!(
        "[gnomad] first use of {file_chrom}: streaming parse {} to build the index cache (one-time, then read cache directly)",
        vcf.display()
    );
    let idx = build_chrom(vcf, file_chrom)?;
    write_cache(&cache, vcf, &idx)?;
    eprintln!(
        "[gnomad] {file_chrom} index done: {} sites → {}",
        idx.pos.len(),
        cache.display()
    );
    Ok(idx)
}

fn load_chrom(d: &GnomadDir, contig: &str) -> Result<(), CallError> {
    if d.chroms
        .read()
        .map(|g| g.contains_key(contig))
        .unwrap_or(false)
    {
        return Ok(());
    }
    // Build outside the lock: index building may take minutes (15GB+ VCF parse), and holding the write lock would freeze queries for all contigs.
    // Concurrent duplicate builds of the same contig only waste CPU; no on-disk conflict (tmp+rename is atomic), acceptable.
    let vcf = resolve_vcf(&d.dir, contig).ok_or_else(|| {
        annot_err(
            &d.dir,
            format!("directory missing the per-chrom VCF for contig {contig} (pattern {DIR_PATTERN}; chr-prefix mapping already tried)"),
        )
    })?;
    let idx = Arc::new(load_or_build(&vcf, &d.cache_dir)?);
    let mut w = d.chroms.write().map_err(|_| CallError::Annot {
        path: d.dir.display().to_string(),
        msg: "gnomad index lock poisoned".into(),
    })?;
    w.entry(contig.to_string()).or_insert(idx);
    Ok(())
}

impl GnomadIndex {
    pub fn load(path: &Path) -> Result<GnomadIndex, CallError> {
        if !path.exists() {
            return Err(annot_err(
                path,
                "gnomad path does not exist (configured, so this is a hard error — not silently skipped)",
            ));
        }
        if path.is_dir() {
            let cache_dir = gnomad_cache_dir()?;
            std::fs::create_dir_all(&cache_dir).map_err(|e| CallError::Io {
                path: cache_dir.display().to_string(),
                source: e,
            })?;
            return Ok(GnomadIndex {
                by_chrom: HashMap::new(),
                dir: Some(GnomadDir {
                    dir: path.to_path_buf(),
                    cache_dir,
                    chroms: RwLock::new(HashMap::new()),
                }),
            });
        }
        // Single file: full load (existing behavior)
        let text = read_text(path)?;
        let mut by_chrom: HashMap<String, HashMap<i64, f64>> = HashMap::new();
        for (ln, line) in text.lines().enumerate() {
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 8 {
                return Err(CallError::Annot {
                    path: path.display().to_string(),
                    msg: format!("line {}: expected 8 VCF columns", ln + 1),
                });
            }
            let pos1: i64 = f[1].parse().map_err(|_| CallError::Annot {
                path: path.display().to_string(),
                msg: format!("line {}: bad POS: {}", ln + 1, f[1]),
            })?;
            let af = f[7]
                .split(';')
                .find_map(|kv| kv.strip_prefix("AF="))
                .and_then(|v| v.split(',').next())
                .and_then(|v| v.parse::<f64>().ok());
            if let Some(af) = af {
                by_chrom
                    .entry(f[0].to_string())
                    .or_default()
                    .insert(pos1 - 1, af); // VCF 1-based → 0-based
            }
        }
        Ok(GnomadIndex {
            by_chrom,
            dir: None,
        })
    }

    /// Directory mode: verify per-chrom VCFs exist for the BAM's primary chromosomes; missing file = hard error (fail fast).
    /// Validates only, does not load — actual loading happens lazily on the first af_at query (BAM headers often list all contigs;
    /// only chromosomes that actually have data are worth indexing). No-op in file mode.
    pub fn prepare(&self, contigs: &[String]) -> Result<(), CallError> {
        let Some(d) = &self.dir else { return Ok(()) };
        for c in contigs {
            if is_primary(c) && resolve_vcf(&d.dir, c).is_none() {
                return Err(annot_err(
                    &d.dir,
                    format!(
                        "directory missing the per-chrom VCF for contig {c} (pattern {DIR_PATTERN}; chr-prefix mapping already tried)"
                    ),
                ));
            }
        }
        Ok(())
    }

    /// Query the AF at pos0. In directory mode the index is lazily built/read on first query of a contig (double-check inside the write lock).
    /// Missing/unparseable file for a primary chromosome = hard error; no file for a non-primary chromosome = Ok(None) (neutral).
    pub fn af_at(&self, chrom: &str, pos0: i64) -> Result<Option<f64>, CallError> {
        if let Some(d) = &self.dir {
            let hit = d
                .chroms
                .read()
                .ok()
                .and_then(|g| g.get(chrom).and_then(|c| c.af_at(pos0)));
            if hit.is_some() {
                return Ok(hit);
            }
            if d.chroms
                .read()
                .map(|g| g.contains_key(chrom))
                .unwrap_or(false)
            {
                return Ok(None); // loaded; site absent from gnomAD
            }
            // Not loaded: non-primary chromosomes usually have no file — neutral; primary chromosomes go through load_chrom (missing = Err).
            if !is_primary(chrom) && resolve_vcf(&d.dir, chrom).is_none() {
                return Ok(None);
            }
            load_chrom(d, chrom)?;
            let g = d.chroms.read().map_err(|_| CallError::Annot {
                path: d.dir.display().to_string(),
                msg: "gnomad index lock poisoned".into(),
            })?;
            return Ok(g.get(chrom).and_then(|c| c.af_at(pos0)));
        }
        Ok(self.by_chrom.get(chrom).and_then(|m| m.get(&pos0).copied()))
    }
}
