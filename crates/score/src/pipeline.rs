//! score pipeline: BAM + sites -> RE_PROB (row order = input order).
//!
//! Reference preload + species guard + (chrom,pos) sorted batching + batch-level rayon;
//! within a batch: pileup -> veto gate -> embed the kept-site sub-batch -> fusion head.
//! Deterministic: thread count / batch size do not change the numerics.
//! Semantics from flow::score_sites_batched (the corrected feature source is out of scope for 1.0.0).

use anyhow::{anyhow, Context, Result};
use rayon::prelude::*;
use std::path::Path;
use std::sync::{Arc, Mutex};

use crate::bundle::{load_bundle, EmbCache};
use crate::caduceus::{CaduceusEncoder, EmbedBatchBufs};
use crate::encoder::{fetch_window_mem_hw, tokenize};
use crate::head::{gate_prob_ensemble, re_prob_ensemble};

/// Batch-level worker: one per rayon thread (BAM reader + reused encoding buffers).
struct ScoreBatchWorker {
    bam: rust_htslib::bam::IndexedReader,
    bufs: EmbedBatchBufs,
}

/// Encoder resolution: `bundle/encoder` or the bundle parent's `encoder` (must contain model.safetensors).
pub fn resolve_encoder_from_bundle(bundle: &Path) -> Result<std::path::PathBuf> {
    let parent = bundle
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| bundle.to_path_buf());
    let cands = [bundle.join("encoder"), parent.join("encoder")];
    cands
        .iter()
        .find(|c| c.join("model.safetensors").exists())
        .cloned()
        .ok_or_else(|| anyhow!("encoder not found in bundle (need encoder/model.safetensors)"))
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

#[allow(clippy::too_many_arguments)]
pub fn score_sites_batched(
    bam: &Path,
    fasta: &Path,
    caduceus: &Path,
    bundle: &Path,
    sites: &[(String, i64)],
    threads: usize,
    batch: usize,
    emb_cache: Option<&Arc<Mutex<EmbCache>>>,
) -> Result<Vec<f64>> {
    let bundle = load_bundle(bundle).context("load bundle")?;
    // The veto gate is mandatory for score (bundle contract; no switch, no fallback).
    let gate = bundle
        .gate
        .as_ref()
        .context("bundle missing the gate triple (required for score since v1.3)")?;
    let enc = Arc::new(CaduceusEncoder::load(caduceus).context("load caduceus")?);
    let batch = batch.max(1);
    // Window token length is decided by the bundle (v1.4.1: 2x250+2 = 502).
    let tok_len = (2 * bundle.half_window + 2) as usize;

    // Reference preload (only the contigs involved in sites).
    let mut need: std::collections::BTreeSet<&str> = std::collections::BTreeSet::new();
    for (c, _) in sites {
        need.insert(c.as_str());
    }
    let mut refmap: std::collections::HashMap<String, Arc<Vec<u8>>> =
        std::collections::HashMap::new();
    {
        // Validate names against .fai first (htslib aborts on a missing contig), then load into memory in one pass.
        let fai_text = std::fs::read_to_string(format!("{}.fai", fasta.display()))
            .with_context(|| format!("read {}.fai", fasta.display()))?;
        // Species guard: the bundle is a human hg38 model; if a real-size reference's chr1 length
        // does not match, refuse to run (hg38 chr1 = 248,956,422); small (<10Mb) is treated as a
        // synthetic/test reference and the check is skipped.
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
        let have: std::collections::HashSet<&str> = fai_text
            .lines()
            .filter_map(|l| l.split('\t').next())
            .collect();
        let fa = rust_htslib::faidx::Reader::from_path(fasta).context("open fasta")?;
        let mut missing = 0usize;
        for c in need {
            if !have.contains(c) {
                // fasta missing this contig -> empty seq; fetch_window_mem degrades to an all-N window; not fatal.
                missing += 1;
                refmap.insert(c.to_string(), Arc::new(Vec::new()));
                continue;
            }
            let len = fa.fetch_seq_len(c) as usize;
            let seq = fa
                .fetch_seq(c, 0, len - 1)
                .with_context(|| format!("fetch {c}"))?;
            refmap.insert(
                c.to_string(),
                Arc::new(seq.iter().map(|b| b.to_ascii_uppercase()).collect()),
            );
        }
        eprintln!(
            "[score] reference preload: {} contigs ({} missing from fasta, degrading to N-pad windows)",
            refmap.len(),
            missing
        );
    }

    // Sorting permutation (contract: output follows input order; internally batched after (chrom,pos) sorting).
    let mut order: Vec<usize> = (0..sites.len()).collect();
    order.sort_by(|&a, &b| {
        sites[a]
            .0
            .cmp(&sites[b].0)
            .then(sites[a].1.cmp(&sites[b].1))
    });
    let sorted: Vec<(String, i64)> = order.iter().map(|&i| sites[i].clone()).collect();

    let pool = match threads {
        0 => rayon::ThreadPoolBuilder::new().build(),
        n => rayon::ThreadPoolBuilder::new().num_threads(n).build(),
    }
    .context("build thread pool")?;
    let nb = sorted.len().div_ceil(batch);
    let t0 = std::time::Instant::now();
    let results: Vec<Result<Vec<f64>, String>> = pool.install(|| {
        (0..nb)
            .into_par_iter()
            .map_init(
                || -> Result<ScoreBatchWorker, String> {
                    Ok(ScoreBatchWorker {
                        bam: rust_htslib::bam::IndexedReader::from_path(bam)
                            .map_err(|e| format!("open bam: {e}"))?,
                        bufs: Default::default(),
                    })
                },
                |w, bi| {
                    let w = match w.as_mut() {
                        Ok(w) => w,
                        Err(e) => return Err(e.clone()),
                    };
                    let lo = bi * batch;
                    let hi = (lo + batch).min(sorted.len());
                    let chunk = &sorted[lo..hi];
                    // Sites in a chunk are fetched in one batch grouped by coordinate (same record set as per-site);
                    // pileup is moved ahead of caching/encoding -- it is the veto gate's input.
                    let refs: Vec<(&str, i64)> =
                        chunk.iter().map(|(c, p)| (c.as_str(), *p)).collect();
                    let piles = esperanto_pile::extract_pileup_features_batch(&mut w.bam, &refs)
                        .map_err(|e| format!("pileup batch: {e}"))?;
                    // Veto gate: gate RE_PROB < threshold -> skip the encoder and emit the gate probability directly;
                    // kept sites take the original path (embed + fusion head, numerically bit-identical to no-gate).
                    let mut final_probs: Vec<Option<f64>> = Vec::with_capacity(chunk.len());
                    let mut keep_idx: Vec<usize> = Vec::new();
                    for (bi2, pl) in piles.iter().enumerate() {
                        let gp = gate_prob_ensemble(gate, pl).map_err(|e| format!("gate: {e}"))?;
                        if gp < gate.threshold {
                            final_probs.push(Some(gp));
                        } else {
                            final_probs.push(None);
                            keep_idx.push(bi2);
                        }
                    }
                    // Cache split (kept sites only) -- a hit yields the fp16 bit pattern (bit-identical to online
                    // embed); misses are collected into the toks batch (order preserved).
                    let mut kept_hit: Vec<Option<[u16; 118]>> =
                        keep_idx.iter().map(|_| None).collect();
                    let mut toks: Vec<i64> = Vec::with_capacity(keep_idx.len() * tok_len);
                    {
                        let mut miss_ki: Vec<usize> = Vec::new();
                        if let Some(c) = emb_cache {
                            let guard = c.lock().unwrap();
                            for (ki, &bi2) in keep_idx.iter().enumerate() {
                                let (chrom, pos) = &chunk[bi2];
                                if let Some(bits) = guard.get(chrom, *pos as u32) {
                                    kept_hit[ki] = Some(bits);
                                } else {
                                    miss_ki.push(ki);
                                }
                            }
                        } else {
                            miss_ki.extend(0..keep_idx.len());
                        }
                        for ki in miss_ki {
                            let (chrom, pos) = &chunk[keep_idx[ki]];
                            let seq = refmap
                                .get(chrom)
                                .ok_or_else(|| format!("fasta missing contig {chrom}"))?;
                            let window = fetch_window_mem_hw(seq, *pos, bundle.half_window);
                            toks.extend(tokenize(&window));
                        }
                    }
                    // Batched encoding: forward in EMBED_SLICE sub-batch slices (numerics are per-site
                    // independent; slicing only changes GEMM shapes).
                    const EMBED_SLICE: usize = 128;
                    let n_miss = toks.len() / tok_len;
                    let mut miss_embs: Vec<[half::f16; 118]> = Vec::with_capacity(n_miss);
                    for sub_toks in toks.chunks(EMBED_SLICE * tok_len) {
                        let n_sub = sub_toks.len() / tok_len;
                        let sub = enc
                            .embed_batch(sub_toks, n_sub, tok_len, &mut w.bufs)
                            .map_err(|e| format!("embed_batch slice: {e}"))?;
                        miss_embs.extend(sub);
                    }
                    // Kept sites: cache hit / freshly computed embedding -> fusion head; write back to cache (new computations only).
                    let mut miss_iter = miss_embs.into_iter();
                    for (ki, &bi2) in keep_idx.iter().enumerate() {
                        let e: [half::f16; 118] = match kept_hit[ki] {
                            Some(bits) => {
                                let mut e = [half::f16::ZERO; 118];
                                for (k, b) in bits.iter().enumerate() {
                                    e[k] = half::f16::from_bits(*b);
                                }
                                e
                            }
                            None => {
                                let e = miss_iter.next().expect("miss count mismatch");
                                if let Some(c) = emb_cache {
                                    let bits: [u16; 118] = std::array::from_fn(|k| e[k].to_bits());
                                    let (chrom, pos) = &chunk[bi2];
                                    c.lock().unwrap().put(chrom, *pos as u32, &bits);
                                }
                                e
                            }
                        };
                        let emb32: Vec<f32> = e.iter().map(|v| v.to_f32()).collect();
                        let emb = ndarray::Array1::from_vec(emb32);
                        let p = re_prob_ensemble(&bundle, &emb.view(), &piles[bi2])
                            .map_err(|e| format!("head: {e}"))?;
                        final_probs[bi2] = Some(p);
                    }
                    let out: Vec<f64> = final_probs
                        .into_iter()
                        .map(|p| p.expect("all sites scored"))
                        .collect();
                    Ok(out)
                },
            )
            .collect()
    });
    eprintln!(
        "[score] batched: {} sites, batch={}, {} batches, {:.1}s",
        sorted.len(),
        batch,
        nb,
        t0.elapsed().as_secs_f64()
    );

    let mut probs = vec![0.0f64; sites.len()];
    for (bi, r) in results.into_iter().enumerate() {
        let ps = r.map_err(|e| anyhow!("batch {bi}: {e}"))?;
        for (j, p) in ps.into_iter().enumerate() {
            probs[order[bi * batch + j]] = p;
        }
    }
    Ok(probs)
}
