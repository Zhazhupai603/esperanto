//! Zero-config resolution (spec §zero-config resolution): bundle 5-level fallback, refs
//! 4-level discovery (fill-empty only), user data dir. Every candidate path
//! is recorded so failures can list what was tried.

use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail};

/// Bundle location relative to a layout root (v1.4.1 contract name).
const BUNDLE_REL: &str = "bundle/human/esperanto-model-v1.4.1-501_40ep/rust";
/// Preferred FASTA name inside a refs directory.
const PREFERRED_FASTA: &str = "hg38.fa";

/// `0` means all cores; resolved once at dispatch (spec convention).
pub fn threads(t: usize) -> usize {
    if t == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    } else {
        t
    }
}

/// The user-data refs directory (`~/.local/share/esperanto/refs`,
/// `ESPERANTO_HOME`/`XDG_DATA_HOME` aware) — the fixed `setup` folder.
pub fn home_refs_dir() -> PathBuf {
    if let Ok(v) = std::env::var("ESPERANTO_REFS") {
        return PathBuf::from(v);
    }
    home_data_dir()
        .unwrap_or_else(|| PathBuf::from("refs"))
        .join("refs")
}

/// The user-data model bundle directory
/// (`~/.local/share/esperanto/bundle`) — the `setup` install target.
pub fn home_bundle_dir() -> PathBuf {
    home_data_dir()
        .unwrap_or_else(|| PathBuf::from("bundle"))
        .join("bundle")
}

fn home_data_dir() -> Option<PathBuf> {
    if let Ok(v) = std::env::var("ESPERANTO_HOME") {
        return Some(PathBuf::from(v));
    }
    if let Ok(v) = std::env::var("XDG_DATA_HOME") {
        return Some(PathBuf::from(v).join("esperanto"));
    }
    std::env::var("HOME")
        .ok()
        .map(|h| PathBuf::from(h).join(".local/share/esperanto"))
}

fn valid_bundle(p: &Path) -> bool {
    p.join("norm.json").is_file()
}

/// Resolve the model bundle: explicit flag wins, else 5-level fallback.
pub fn bundle(explicit: &Option<PathBuf>) -> anyhow::Result<PathBuf> {
    if let Some(b) = explicit {
        return Ok(b.clone());
    }
    let mut tried: Vec<PathBuf> = Vec::new();
    if let Ok(env) = std::env::var("ESPERANTO_BUNDLE") {
        tried.push(PathBuf::from(env));
    }
    tried.push(
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../")
            .join(BUNDLE_REL),
    );
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            tried.push(dir.join(BUNDLE_REL));
            if let Some(up) = dir.parent() {
                tried.push(up.join(BUNDLE_REL));
            }
        }
    }
    if let Some(home) = home_data_dir() {
        tried.push(home.join(BUNDLE_REL));
    }
    for cand in &tried {
        if valid_bundle(cand) {
            return Ok(cand.clone());
        }
    }
    bail!(
        "model bundle not found (need norm.json); tried:\n  {}",
        tried
            .iter()
            .map(|p| p.display().to_string())
            .collect::<Vec<_>>()
            .join("\n  ")
    )
}

/// Resolved files inside a refs directory.
pub struct Refs {
    /// `hg38.fa` or the first `*.fa` with a `.fai` alongside.
    pub fasta: Option<PathBuf>,
    /// First `*.gtf`.
    pub gtf: Option<PathBuf>,
    /// First `dbsnp*`/`gnomad*` `*.vcf.gz`.
    pub gnomad: Option<PathBuf>,
    /// First `*.paidx` alignment index.
    pub index: Option<PathBuf>,
}

fn refs_in(dir: &Path) -> Option<Refs> {
    let index = {
        let mut v: Vec<PathBuf> = std::fs::read_dir(dir)
            .ok()?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| p.extension().is_some_and(|x| x == "paidx"))
            .collect();
        v.sort();
        v.into_iter().next()
    };
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .ok()?
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    let has_fai = |fa: &str| dir.join(format!("{fa}.fai")).is_file();
    let fasta = names
        .iter()
        .find(|n| n.as_str() == PREFERRED_FASTA && has_fai(n))
        .or_else(|| names.iter().find(|n| n.ends_with(".fa") && has_fai(n)))
        .map(|n| dir.join(n))?; // a refs dir without a usable fasta is not a refs dir
    let gtf = names
        .iter()
        .find(|n| n.ends_with(".gtf"))
        .map(|n| dir.join(n));
    let gnomad = names
        .iter()
        .find(|n| n.ends_with(".vcf.gz") && (n.starts_with("dbsnp") || n.starts_with("gnomad")))
        .map(|n| dir.join(n));
    Some(Refs {
        fasta: Some(fasta),
        gtf,
        gnomad,
        index,
    })
}

/// Discover a refs directory (4 levels; first valid wins).
pub fn refs() -> Option<Refs> {
    let mut candidates: Vec<PathBuf> = Vec::new();
    if let Ok(env) = std::env::var("ESPERANTO_REFS") {
        candidates.push(PathBuf::from(env));
    }
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            candidates.push(dir.join("refs"));
            if let Some(up) = dir.parent() {
                candidates.push(up.join("refs"));
            }
        }
    }
    if let Some(home) = home_data_dir() {
        candidates.push(home.join("refs"));
    }
    candidates.push(PathBuf::from("./refs"));
    candidates.iter().find_map(|d| refs_in(d))
}

impl Refs {
    /// Fill `fasta` when unset (explicit values are never overwritten).
    pub fn fill_fasta(&self, fasta: &mut Option<PathBuf>) {
        if fasta.is_none() {
            *fasta = self.fasta.clone();
        }
    }

    /// Fill `gtf` when unset.
    pub fn fill_gtf(&self, gtf: &mut Option<PathBuf>) {
        if gtf.is_none() {
            *gtf = self.gtf.clone();
        }
    }

    /// Fill `gnomad` when unset.
    pub fn fill_gnomad(&self, gnomad: &mut Option<PathBuf>) {
        if gnomad.is_none() {
            *gnomad = self.gnomad.clone();
        }
    }

    /// Fill `index` when unset.
    pub fn fill_index(&self, index: &mut Option<PathBuf>) {
        if index.is_none() {
            *index = self.index.clone();
        }
    }
}

/// A required fasta after resolution, or an actionable error.
pub fn require_fasta(fasta: Option<PathBuf>) -> anyhow::Result<PathBuf> {
    fasta.ok_or_else(|| {
        anyhow!("--fasta not given and no refs directory found (set ESPERANTO_REFS)")
    })
}

/// Resolve the L1 engine bundle (default-on): explicit `--l1-bundle` wins,
/// else `ESPERANTO_L1_BUNDLE`, else `<index stem>.bndl` next to the index.
/// Returns None (with a stderr note) when nothing is found — the pipeline
/// then runs the pure genomic (G) layer.
pub fn l1_bundle(explicit: &Option<PathBuf>, index: &Path) -> Option<PathBuf> {
    if let Some(b) = explicit {
        return Some(b.clone());
    }
    if let Ok(env) = std::env::var("ESPERANTO_L1_BUNDLE") {
        let p = PathBuf::from(env);
        if p.is_file() {
            return Some(p);
        }
    }
    let sibling = index.with_extension("bndl");
    if sibling.is_file() {
        // The runtime needs the .tidx sidecar next to the .bndl.
        if sibling.with_extension("tidx").is_file() {
            return Some(sibling);
        }
        eprintln!(
            "[resolve] {} has no .tidx sidecar; running pure G layer",
            sibling.display()
        );
        return None;
    }
    eprintln!(
        "[resolve] L1 bundle not found (tried ESPERANTO_L1_BUNDLE, {}); running pure G layer",
        sibling.display()
    );
    None
}
