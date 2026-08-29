//! `--hybrid` reference builder (spec §hybrid reference): splices selected
//! human gene loci onto the mouse baseline, cached per gene set under
//! `refs/hybrid/<key>/`. A cache hit (species.json + paidx present) is reused
//! with zero rebuild.

use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context};

/// One gene locus on the human reference (1-based inclusive span).
struct Locus {
    symbol: String,
    contig: String,
    start: u64,
    end: u64,
}

/// Resolve the hybrid reference directory for `genes` (None → interactive
/// picker), building it on first use. Returns the directory holding the
/// index (`hybrid.paidx`), fasta, gtf, and species.json.
pub fn resolve(genes: Option<Vec<String>>) -> anyhow::Result<PathBuf> {
    let refs = crate::resolve::home_refs_dir();
    let human_fa = tagged(&refs, "human", &[".fa", ".fasta"])?;
    let human_gtf = tagged(&refs, "human", &[".gtf"])?;
    let mouse_fa = tagged(&refs, "mouse", &[".fa", ".fasta"])?;
    let mouse_gtf = tagged(&refs, "mouse", &[".gtf"]).ok(); // optional
    let genes = match genes {
        Some(g) => g,
        None => {
            let names = gene_names(&human_gtf)?;
            crate::confirm::multi_pick("Select human gene(s) knocked into the mouse:", &names)?
                .ok_or_else(|| anyhow!("no gene selected (aborted)"))?
        }
    };
    if genes.is_empty() {
        bail!("--hybrid needs at least one gene symbol");
    }
    let mut sorted = genes.clone();
    sorted.sort();
    sorted.dedup();
    // Baseline tag follows the staged mouse assembly (name-carried):
    // mm10-named files keep mm10 coordinates, grcm39-named keep grcm39.
    let baseline = if mouse_fa
        .file_name()
        .is_some_and(|n| n.to_string_lossy().to_lowercase().contains("mm10"))
    {
        "mm10"
    } else {
        "grcm39"
    };
    let key = format!("{baseline}+{}", sorted.join("+"));
    let dir = refs.join("hybrid").join(&key);
    if dir.join("species.json").is_file() && dir.join("hybrid.paidx").is_file() {
        eprintln!("[hybrid] {} found; reusing the built index", dir.display());
        return Ok(dir);
    }
    build(&dir, baseline, &human_fa, &human_gtf, &mouse_fa, mouse_gtf, &sorted)?;
    eprintln!("[hybrid] built {} -> {}", key, dir.display());
    Ok(dir)
}

/// The single fasta/gtf tagged `species` in `dir`; actionable error when
/// missing (setup must have run first) or ambiguous.
fn tagged(dir: &Path, species: &str, exts: &[&str]) -> anyhow::Result<PathBuf> {
    let hits: Vec<PathBuf> = crate::resolve::find_all(dir, exts)?
        .into_iter()
        .filter(|p| {
            p.file_name()
                .is_some_and(|n| crate::resolve::species_of(&n.to_string_lossy()) == Some(species))
        })
        .collect();
    if hits.is_empty() && species == "human" {
        // Default installs carry an unmarked human GTF (gencode.vXX);
        // accept a single unmarked file when nothing species-tagged exists.
        let unmarked: Vec<PathBuf> = crate::resolve::find_all(dir, exts)?
            .into_iter()
            .filter(|p| {
                p.file_name()
                    .is_some_and(|n| crate::resolve::species_of(&n.to_string_lossy()).is_none())
            })
            .collect();
        if unmarked.len() == 1 {
            return Ok(unmarked.into_iter().next().expect("len 1"));
        }
    }
    match hits.len() {
        1 => Ok(hits.into_iter().next().expect("len 1")),
        0 => bail!(
            "no {species} {} found in {}; run `esperanto setup` first",
            exts.join("/"),
            dir.display()
        ),
        _ => bail!(
            "multiple {species} {} files in {}; keep exactly one",
            exts.join("/"),
            dir.display()
        ),
    }
}

/// All gene symbols in the human GTF (for the picker), sorted.
fn gene_names(gtf: &Path) -> anyhow::Result<Vec<String>> {
    let text = std::fs::read_to_string(gtf).context("read human gtf")?;
    let mut out: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let mut cols = line.split('\t');
        if cols.nth(2) != Some("gene") {
            continue;
        }
        if let Some(attrs) = line.split('\t').nth(8) {
            if let Some(name) = attr_value(attrs, "gene_name") {
                out.insert(name);
            }
        }
    }
    Ok(out.into_iter().collect())
}

/// `attr_value("gene_id \"APOE4\"; gene_name \"APOE4\";", "gene_name")` → `APOE4`.
fn attr_value(attrs: &str, key: &str) -> Option<String> {
    for part in attrs.split(';') {
        let part = part.trim();
        let mut kv = part.splitn(2, ' ');
        if kv.next() != Some(key) {
            continue;
        }
        return kv.next().map(|v| v.trim().trim_matches('"').to_string());
    }
    None
}

/// Gene spans from the human GTF (gene rows), keyed by symbol.
fn gene_spans(gtf: &Path, wanted: &[String]) -> anyhow::Result<Vec<Locus>> {
    let text = std::fs::read_to_string(gtf).context("read human gtf")?;
    let mut out = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 || cols[2] != "gene" {
            continue;
        }
        let Some(name) = attr_value(cols[8], "gene_name") else {
            continue;
        };
        if !wanted.iter().any(|w| w == &name) {
            continue;
        }
        out.push(Locus {
            symbol: name,
            contig: cols[0].to_string(),
            start: cols[3].parse().context("gtf start")?,
            end: cols[4].parse().context("gtf end")?,
        });
    }
    let mut missing: Vec<&String> = wanted
        .iter()
        .filter(|w| !out.iter().any(|l| &l.symbol == *w))
        .collect();
    if !missing.is_empty() {
        missing.sort();
        bail!(
            "gene symbol(s) not found in the human GTF: {}",
            missing
                .iter()
                .map(|s| s.as_str())
                .collect::<Vec<_>>()
                .join(", ")
        );
    }
    Ok(out)
}

/// Build the hybrid reference into `dir`: mouse baseline + one contig per
/// human locus (`h<SYMBOL>`), annotation rebased onto the new contigs,
/// species.json, and the full index set.
fn build(
    dir: &Path,
    baseline: &str,
    human_fa: &Path,
    human_gtf: &Path,
    mouse_fa: &Path,
    mouse_gtf: Option<PathBuf>,
    genes: &[String],
) -> anyhow::Result<()> {
    std::fs::create_dir_all(dir)?;
    let loci = gene_spans(human_gtf, genes)?;

    // --- hybrid fasta: mouse baseline + h* loci ---
    let fa = rust_htslib::faidx::Reader::from_path(human_fa)
        .map_err(|e| anyhow!("open human fasta: {e}"))?;
    let hybrid_fa = dir.join("hybrid.fa");
    let mouse_bytes = std::fs::read(mouse_fa).context("read mouse fasta")?;
    let mut out: Vec<u8> = mouse_bytes;
    let mut owner = Vec::new();
    for l in &loci {
        let seq = fa
            .fetch_seq(&l.contig, (l.start - 1) as usize, (l.end - 1) as usize)
            .map_err(|e| anyhow!("fetch {}:{}-{}: {e}", l.contig, l.start, l.end))?;
        out.extend_from_slice(format!(">h{}\n", l.symbol).as_bytes());
        for chunk in seq.chunks(60) {
            out.extend_from_slice(chunk);
            out.push(b'\n');
        }
        owner.push(format!("h{}", l.symbol));
    }
    std::fs::write(&hybrid_fa, &out).context("write hybrid fasta")?;
    crate::setup::write_fai(&hybrid_fa)?;

    // --- hybrid gtf: mouse baseline + rebased human-locus rows ---
    let hybrid_gtf = dir.join("hybrid.gtf");
    let mut gtf_out = Vec::new();
    if let Some(mg) = &mouse_gtf {
        gtf_out.extend_from_slice(&std::fs::read(mg).context("read mouse gtf")?);
    }
    let hg = std::fs::read_to_string(human_gtf).context("read human gtf")?;
    for line in hg.lines() {
        if line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }
        let Some(name) = attr_value(cols[8], "gene_name") else {
            continue;
        };
        let Some(l) = loci.iter().find(|l| l.symbol == name) else {
            continue;
        };
        if cols[0] != l.contig {
            continue; // a gene must live on one contig
        }
        let s: u64 = cols[3].parse().context("gtf start")?;
        let e: u64 = cols[4].parse().context("gtf end")?;
        if s < l.start || e > l.end {
            continue; // feature hangs outside the extracted locus
        }
        let rebased = format!(
            "h{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            l.symbol,
            cols[1],
            cols[2],
            s - l.start + 1,
            e - l.start + 1,
            cols[5],
            cols[6],
            cols[7],
            cols[8]
        );
        gtf_out.extend_from_slice(rebased.as_bytes());
    }
    std::fs::write(&hybrid_gtf, &gtf_out).context("write hybrid gtf")?;

    // --- species.json ---
    let mut manifest = esperanto_flow::manifest::SpeciesManifest::single("hybrid", baseline);
    manifest.human_loci = owner.clone();
    let fai_text = std::fs::read_to_string(dir.join("hybrid.fa.fai"))?;
    for line in fai_text.lines() {
        let Some(name) = line.split('\t').next() else {
            continue;
        };
        let owner_tag = if owner.iter().any(|o| o == name) {
            "human"
        } else {
            "mouse"
        };
        manifest
            .contig_owner
            .insert(name.to_string(), owner_tag.to_string());
    }
    let human_bundle = crate::resolve::bundle(&None)?;
    let human_bundle = human_bundle
        .canonicalize()
        .unwrap_or(human_bundle);
    manifest.bundles.insert("human".to_string(), human_bundle);
    if let Some(mb) = crate::resolve::mouse_bundle() {
        manifest.bundles.insert("mouse".to_string(), mb);
    }
    manifest.write(dir)?;

    // --- guardrail + full index set ---
    esperanto_flow::guard::check_species(&hybrid_fa, Some(&manifest))
        .map_err(|e| anyhow!("{}: {e}", hybrid_fa.display()))?;
    crate::index::build_all(
        &hybrid_fa,
        Some(&hybrid_gtf),
        &dir.join("hybrid.paidx"),
        15,
        5,
    )?;
    Ok(())
}
