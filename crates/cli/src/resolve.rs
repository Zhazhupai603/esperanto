//! Zero-config resolution (spec §zero-config resolution): bundle 5-level fallback, refs
//! 4-level discovery (fill-empty only), user data dir. Every candidate path
//! is recorded so failures can list what was tried.

use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail};

/// Unified model bundle location relative to a layout root.
const BUNDLE_REL: &str = "bundle/human/esperanto-model-v1.1.0/rust";
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
    p.join("model.txt").is_file() && p.join("freeze.json").is_file() && p.join("alu.bed").is_file()
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
        "unified model bundle not found (need model.txt, freeze.json, alu.bed); tried:\n  {}",
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
    /// First `*.paidx` alignment index.
    pub index: Option<PathBuf>,
}

fn refs_in(dir: &Path) -> Option<Refs> {
    let index = {
        let mut v: Vec<PathBuf> = std::fs::read_dir(dir)
            .ok()?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| {
                p.extension().is_some_and(|x| x == "paidx")
                    && !p.file_name().is_some_and(|n| {
                        let n = n.to_string_lossy();
                        [".ag.paidx", ".tc.paidx"]
                            .iter()
                            .any(|suffix| n.ends_with(suffix))
                    })
            })
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
        .or_else(|| {
            names
                .iter()
                .find(|n| (n.ends_with(".fa") || n.ends_with(".fasta")) && has_fai(n))
        })
        .map(|n| dir.join(n))?; // a refs dir without a usable fasta is not a refs dir
    let gtf = names
        .iter()
        .find(|n| n.ends_with(".gtf"))
        .map(|n| dir.join(n));
    Some(Refs {
        fasta: Some(fasta),
        gtf,
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
/// The L1 engine is required: a missing `.bndl` or a `.bndl` without its
/// `.tidx` sidecar is a hard error rather than a silent downgrade to the
/// pure genomic (G) layer.
pub fn l1_bundle(explicit: &Option<PathBuf>, index: &Path) -> anyhow::Result<PathBuf> {
    let candidate = if let Some(b) = explicit {
        b.clone()
    } else if let Ok(env) = std::env::var("ESPERANTO_L1_BUNDLE") {
        let p = PathBuf::from(env);
        if p.is_file() {
            p
        } else {
            let sibling = index.with_extension("bndl");
            if sibling.is_file() {
                sibling
            } else {
                anyhow::bail!(
                    "L1 bundle not found (tried ESPERANTO_L1_BUNDLE, {}); build it with `esperanto index --gtf <gtf>`",
                    sibling.display()
                );
            }
        }
    } else {
        let sibling = index.with_extension("bndl");
        if sibling.is_file() {
            sibling
        } else {
            anyhow::bail!(
                "L1 bundle not found (tried ESPERANTO_L1_BUNDLE, {}); build it with `esperanto index --gtf <gtf>`",
                sibling.display()
            );
        }
    };
    let tidx = candidate.with_extension("tidx");
    if !tidx.is_file() {
        anyhow::bail!(
            "{} has no .tidx sidecar; the L1 engine refuses to degrade to the pure G layer — rebuild with `esperanto index --gtf <gtf>`",
            candidate.display()
        );
    }
    Ok(candidate)
}

/// All files in `dir` matching one of `exts` (sorted).
pub fn find_all(dir: &Path, exts: &[&str]) -> anyhow::Result<Vec<PathBuf>> {
    let mut hits: Vec<PathBuf> = std::fs::read_dir(dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.is_file()
                && exts.iter().any(|x| {
                    p.file_name()
                        .is_some_and(|n| n.to_string_lossy().ends_with(x))
                })
        })
        .collect();
    hits.sort();
    Ok(hits)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn discovery_selects_genomic_index_and_accepts_fasta_extension() {
        let dir = std::env::temp_dir().join(format!("esperanto-refs-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        for name in [
            "reference.fasta",
            "reference.fasta.fai",
            "reference.paidx",
            "reference.ag.paidx",
            "reference.tc.paidx",
        ] {
            std::fs::write(dir.join(name), "").unwrap();
        }
        let refs = refs_in(&dir).unwrap();
        assert_eq!(refs.index, Some(dir.join("reference.paidx")));
        assert_eq!(refs.fasta, Some(dir.join("reference.fasta")));
        std::fs::remove_file(dir.join("reference.paidx")).unwrap();
        assert!(refs_in(&dir).unwrap().index.is_none());
        std::fs::remove_dir_all(dir).unwrap();
    }
}
