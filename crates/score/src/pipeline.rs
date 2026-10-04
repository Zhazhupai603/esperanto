//! Unified reference-aware scoring pipeline: BAM + sites -> classifier score.

use anyhow::{Context, Result};
use std::path::Path;

/// sites TSV parsing: chrom<TAB>pos (1-based, >=1); empty lines skipped.
pub fn parse_sites(text: &str) -> Result<Vec<(String, i64)>> {
    let mut out = Vec::new();
    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        let (chrom, pos) = line
            .split_once('\t')
            .with_context(|| format!("bad site line: {line}"))?;
        let pos: i64 = pos
            .parse()
            .with_context(|| format!("parse pos in: {line}"))?;
        anyhow::ensure!(
            pos >= 1,
            "bad site line: {line} — pos must be >= 1 (1-based)"
        );
        out.push((chrom.to_string(), pos));
    }
    Ok(out)
}

/// Species guardrail on the FASTA index text: the model is a human hg38
/// model; a real-size reference whose chr1 does not match is refused, small
/// (<10 Mb) references are treated as synthetic/test and pass.
fn species_guardrail(fai_text: &str) -> Result<()> {
    if let Some(l) = fai_text.lines().find(|l| l.starts_with("chr1\t")) {
        let len: u64 = l
            .split('\t')
            .nth(1)
            .and_then(|v| v.parse().ok())
            .unwrap_or(0);
        anyhow::ensure!(
            len < 10_000_000 || len == 248_956_422,
            "species/version mismatch: bundle is a human hg38 model, but reference chr1 length = {len} (hg38 should be 248956422). Use an hg38 reference, or a model bundle matching this reference."
        );
    }
    Ok(())
}

fn score_unified_batched(
    bam: &Path,
    fasta: &Path,
    bundle_dir: &Path,
    sites: &[(String, i64)],
    threads: usize,
) -> Result<Vec<f64>> {
    let fai_text = std::fs::read_to_string(format!("{}.fai", fasta.display()))
        .with_context(|| format!("read {}.fai", fasta.display()))?;
    species_guardrail(&fai_text)?;
    let alu = crate::unified::AluIndex::from_bed(&bundle_dir.join("alu.bed"))
        .context("load unified Alu annotation")?;
    crate::unified::score_sites(bam, fasta, bundle_dir, &alu, sites, threads)
}

/// Score candidate sites with the unified classifier.
pub fn score_sites_batched(
    bam: &Path,
    fasta: &Path,
    bundle: &Path,
    sites: &[(String, i64)],
    threads: usize,
) -> Result<Vec<f64>> {
    if !crate::bundle::is_unified_bundle(bundle) {
        anyhow::bail!(
            "unified model bundle not found or incomplete (need model.txt, freeze.json, alu.bed): {}",
            bundle.display()
        );
    }
    score_unified_batched(bam, fasta, bundle, sites, threads)
}
