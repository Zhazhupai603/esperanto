//! `species.json` — reference manifest (spec §reference manifest). Written by
//! the reference builder into the refs dir and copied into each run dir; the
//! single source of truth for reference kind and per-contig ownership.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

/// Reference kind: single species or a hybrid (baseline + human loci).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SpeciesManifest {
    /// "human" | "mouse" | "hybrid".
    pub kind: String,
    /// Baseline assembly tag: "grch38" | "grcm39".
    pub baseline: String,
    /// Human locus contigs (hybrid only), e.g. ["hAPOE4"].
    #[serde(default)]
    pub human_loci: Vec<String>,
    /// contig → owner ("human" | "mouse").
    #[serde(default)]
    pub contig_owner: BTreeMap<String, String>,
    /// owner → model bundle path (rust dir); written by the reference
    /// builder so score stages resolve the right model per contig.
    #[serde(default)]
    pub bundles: BTreeMap<String, PathBuf>,
}

impl SpeciesManifest {
    /// Plain single-species manifest.
    pub fn single(kind: &str, baseline: &str) -> Self {
        SpeciesManifest {
            kind: kind.to_string(),
            baseline: baseline.to_string(),
            human_loci: Vec::new(),
            contig_owner: BTreeMap::new(),
            bundles: BTreeMap::new(),
        }
    }

    /// Owner tag of a contig ("human" | "mouse"); None when unlisted.
    pub fn owner_of(&self, contig: &str) -> Option<&str> {
        self.contig_owner.get(contig).map(String::as_str)
    }

    /// Write `<dir>/species.json` atomically (tmp + rename).
    pub fn write(&self, dir: &Path) -> std::io::Result<()> {
        let text = serde_json::to_string_pretty(self)
            .map_err(|e| std::io::Error::other(format!("species.json serialize: {e}")))?;
        let tmp = dir.join(".species.json.tmp");
        std::fs::write(&tmp, format!("{text}\n"))?;
        std::fs::rename(&tmp, dir.join("species.json"))
    }

    /// Read `<dir>/species.json`; None when absent or unparsable.
    pub fn read(dir: &Path) -> Option<Self> {
        let text = std::fs::read_to_string(dir.join("species.json")).ok()?;
        serde_json::from_str(&text).ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn manifest_roundtrip() {
        let dir = std::env::temp_dir().join(format!("esperanto-manifest-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut m = SpeciesManifest::single("hybrid", "grcm39");
        m.human_loci = vec!["hAPOE4".to_string()];
        m.contig_owner.insert("hAPOE4".into(), "human".into());
        m.contig_owner.insert("chr1".into(), "mouse".into());
        m.write(&dir).unwrap();
        let got = SpeciesManifest::read(&dir).unwrap();
        assert_eq!(got.kind, "hybrid");
        assert_eq!(got.owner_of("hAPOE4"), Some("human"));
        assert_eq!(got.owner_of("chr1"), Some("mouse"));
        assert!(SpeciesManifest::read(&dir.join("nonexistent")).is_none());
        let _ = std::fs::remove_dir_all(&dir);
    }
}
