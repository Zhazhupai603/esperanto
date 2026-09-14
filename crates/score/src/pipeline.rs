//! score pipeline: BAM + sites -> RE_PROB (row order = input order).
//!
//! CMC feature extraction straight from the BAM (region-grouped sweep, one
//! IndexedReader per worker) followed by the frozen network forward, parallel
//! across sites. Deterministic: thread count does not change the numerics —
//! every site is scored independently and written back by its original index.
//!
//! The FASTA is used for the species guardrail only; the model reads no
//! reference sequence.

use anyhow::{anyhow, Context, Result};
use rayon::prelude::*;
use std::path::Path;

/// Species-guard mode for [`score_sites_batched`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ReferenceCheck {
    /// Default: the fasta's chr1 must match hg38 (or be a small synthetic reference).
    #[default]
    Guardrail,
    /// Hybrid reference verified upstream against its species.json manifest.
    TrustedHybrid,
}

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
fn species_guardrail(fai_text: &str, reference_check: ReferenceCheck) -> Result<()> {
    if reference_check == ReferenceCheck::TrustedHybrid {
        return Ok(());
    }
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

/// Which network of the bundle scores the sites.
#[derive(Clone, Copy)]
enum CmcTask {
    Main,
    Hyper,
}

/// Score sites with the CMC bundle: region-sweep extraction (BAM only)
/// followed by the network forward, parallel across sites. Output order
/// follows the input `sites`.
fn score_cmc_batched(
    bam: &Path,
    fasta: &Path,
    bundle_dir: &Path,
    sites: &[(String, i64)],
    threads: usize,
    reference_check: ReferenceCheck,
    task: CmcTask,
) -> Result<Vec<f64>> {
    let bundle = crate::bundle::load_cmc_bundle(bundle_dir).context("load bundle")?;
    let net = match task {
        CmcTask::Main => &bundle.main,
        CmcTask::Hyper => &bundle.hyper,
    };
    let fai_text = std::fs::read_to_string(format!("{}.fai", fasta.display()))
        .with_context(|| format!("read {}.fai", fasta.display()))?;
    species_guardrail(&fai_text, reference_check)?;
    if sites.is_empty() {
        return Ok(Vec::new());
    }
    let mut order: Vec<usize> = (0..sites.len()).collect();
    order.sort_by(|&a, &b| {
        sites[a]
            .0
            .cmp(&sites[b].0)
            .then(sites[a].1.cmp(&sites[b].1))
    });
    let sorted: Vec<(String, i64)> = order.iter().map(|&i| sites[i].clone()).collect();

    let t = std::time::Instant::now();
    let feats = cmc_extract_all(bam, &sorted, threads)?;
    eprintln!(
        "[score] cmc features: {} sites, {:.1}s",
        sorted.len(),
        t.elapsed().as_secs_f64()
    );

    let t0 = std::time::Instant::now();
    let pool = match threads {
        0 => rayon::ThreadPoolBuilder::new().build(),
        n => rayon::ThreadPoolBuilder::new().num_threads(n).build(),
    }
    .context("build thread pool")?;
    let probs_sorted = pool.install(|| {
        feats
            .par_iter()
            .map(|x| net.re_prob(x))
            .collect::<Vec<f64>>()
    });
    eprintln!(
        "[score] cmc forward: {} sites, device=cpu, {:.1}s",
        sorted.len(),
        t0.elapsed().as_secs_f64()
    );

    let mut probs = vec![0.0f64; sites.len()];
    for (i, p) in probs_sorted.into_iter().enumerate() {
        probs[order[i]] = p;
    }
    Ok(probs)
}

/// Region-grouped CMC extraction (one fetch per run of sites within
/// `CMC_MERGE_GAP`, capped group size, one IndexedReader per worker).
fn cmc_extract_all(
    bam: &Path,
    sorted: &[(String, i64)],
    threads: usize,
) -> Result<Vec<crate::cmc::CmcFeatures>> {
    if sorted.is_empty() {
        return Ok(Vec::new());
    }
    const GROUP_CAP: usize = 16_384;
    let mut groups: Vec<(usize, usize)> = Vec::new();
    let mut start = 0usize;
    for i in 1..=sorted.len() {
        let boundary = i == sorted.len()
            || i - start >= GROUP_CAP
            || sorted[i].0 != sorted[i - 1].0
            || sorted[i].1 - sorted[i - 1].1 > crate::cmc::CMC_MERGE_GAP;
        if boundary {
            groups.push((start, i));
            start = i;
        }
    }
    let pool = match threads {
        0 => rayon::ThreadPoolBuilder::new().build(),
        n => rayon::ThreadPoolBuilder::new().num_threads(n).build(),
    }
    .context("build cmc thread pool")?;
    let per_group: Vec<Result<Vec<crate::cmc::CmcFeatures>, String>> = pool.install(|| {
        groups
            .par_iter()
            .map_init(
                || -> Result<rust_htslib::bam::IndexedReader, String> {
                    rust_htslib::bam::IndexedReader::from_path(bam)
                        .map_err(|e| format!("open bam: {e}"))
                },
                |r, &(lo, hi)| {
                    let r = match r.as_mut() {
                        Ok(r) => r,
                        Err(e) => return Err(e.clone()),
                    };
                    let refs: Vec<(&str, i64)> = sorted[lo..hi]
                        .iter()
                        .map(|(c, p)| (c.as_str(), *p))
                        .collect();
                    crate::cmc::extract_cmc_batch(r, &refs).map_err(|e| format!("cmc group: {e}"))
                },
            )
            .collect()
    });
    let mut out = Vec::with_capacity(sorted.len());
    for (gi, g) in per_group.into_iter().enumerate() {
        let g = g.map_err(|e| anyhow!("cmc group {gi}: {e}"))?;
        out.extend_from_slice(&g);
    }
    Ok(out)
}

/// Score candidate sites with the main (DM vs RE) network.
#[allow(clippy::too_many_arguments)]
pub fn score_sites_batched(
    bam: &Path,
    fasta: &Path,
    bundle: &Path,
    sites: &[(String, i64)],
    threads: usize,
    reference_check: ReferenceCheck,
) -> Result<Vec<f64>> {
    if !crate::bundle::is_cmc_bundle(bundle) {
        anyhow::bail!(
            "model bundle not found or unsupported (need rust/net_main.safetensors): {}",
            bundle.display()
        );
    }
    score_cmc_batched(
        bam,
        fasta,
        bundle,
        sites,
        threads,
        reference_check,
        CmcTask::Main,
    )
}

/// Score Tier-2-only (rescue) sites with the hyper network of the bundle.
#[allow(clippy::too_many_arguments)]
pub fn score_hyper_sites_batched(
    bam: &Path,
    fasta: &Path,
    bundle: &Path,
    sites: &[(String, i64)],
    threads: usize,
    reference_check: ReferenceCheck,
) -> Result<Vec<f64>> {
    if !crate::bundle::is_cmc_bundle(bundle) {
        anyhow::bail!(
            "model bundle not found or unsupported (need rust/net_main.safetensors): {}",
            bundle.display()
        );
    }
    score_cmc_batched(
        bam,
        fasta,
        bundle,
        sites,
        threads,
        reference_check,
        CmcTask::Hyper,
    )
}
