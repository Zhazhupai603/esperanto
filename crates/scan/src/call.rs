//! Call driver — wires the scatter engine to scoring + candidates.bed output.
//!
//! Per (contig, 32Mbp chunk): load refseq (cached) → scatter_block → per-site
//! feature assembly → spec.score → Candidate. Semantics mirror process_column
//! (legacy htslib path) point for point.

use crate::count::{self, SeqCache};
use crate::error::CallError;
use crate::out::{self, Candidate};
use crate::scatter::{self, BlockAcc};
use crate::score::{CallSpec, ScoreFeatures};
use crate::strand::{infer_strand, StrandCall};
use crate::{annot, CallParams, LibType};
use rayon::prelude::*;
static T_REFSEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static T_SCATTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static T_EMIT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
use rust_htslib::bam::Read as BamRead;
use std::collections::BTreeSet;
use std::path::Path;
use std::sync::Arc;

/// Load one contig's reference sequence (shared cache; None when no fasta or
/// contig missing — majority pseudo-ref fallback, mirroring legacy).
fn load_refseq(
    chrom: &str,
    len: u32,
    fasta: Option<&Path>,
    cache: &SeqCache,
) -> Option<Arc<Vec<u8>>> {
    let fa = fasta?;
    if let Some(s) = cache.read().ok().and_then(|c| c.get(chrom).cloned()) {
        return Some(s);
    }
    // P0 (observed on STAR BAMs): when the fasta lacks this contig, htslib faidx_fetch_seq64 sets a negative length,
    // and rust-htslib 0.49 fetch_seq calls Vec::from_raw_parts(len as usize) internally, panicking on
    // capacity overflow — the Result semantics break down (.ok() cannot catch it). Pre-check against the .fai list
    // (missing → None; this contig degrades to
    // majority pseudo-ref counting, never panics.
    if !count::fasta_has_contig(fa, chrom) {
        return None;
    }
    let faidx = rust_htslib::faidx::Reader::from_path(fa).ok()?;
    let s = faidx
        .fetch_seq(chrom, 0usize, (len as usize).saturating_sub(1))
        .ok()?;
    let seq: Vec<u8> = s.iter().map(|b| b.to_ascii_uppercase()).collect();
    let arc = Arc::new(seq);
    if let Ok(mut w) = cache.write() {
        w.entry(chrom.to_string()).or_insert_with(|| arc.clone());
    }
    Some(arc)
}

/// First two .fai columns → true contig length table (correction source for baln-derived lengths).
fn fai_lengths(fasta: Option<&Path>) -> std::collections::HashMap<String, u32> {
    let mut m = std::collections::HashMap::new();
    if let Some(fa) = fasta {
        let fai = format!("{}.fai", fa.display());
        if let Ok(text) = std::fs::read_to_string(&fai) {
            for line in text.lines() {
                let mut it = line.split('\t');
                if let (Some(name), Some(len)) = (it.next(), it.next()) {
                    if let Ok(l) = len.parse::<u32>() {
                        m.insert(name.to_string(), l);
                    }
                }
            }
        }
    }
    m
}

fn base_idx(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Assemble candidates for one chunk from its scatter result.
#[allow(clippy::too_many_arguments)]
fn emit_chunk(
    chrom: &str,
    cs: u64,
    acc: &BlockAcc,
    refseq: Option<&[u8]>,
    lib: LibType,
    enable_cu: bool,
    gtf: Option<&annot::GtfIndex>,
    gnomad: Option<&annot::GnomadIndex>,
    spec: &CallSpec,
    min_call_score: Option<f64>,
    out: &mut Vec<Candidate>,
) -> Result<(), CallError> {
    let jbounds: BTreeSet<i64> = acc.jbounds.iter().map(|&x| x as i64).collect();
    let mut positions: Vec<u32> = acc.sites.keys().copied().collect();
    positions.sort_unstable();

    for gpos in positions {
        let li = (gpos as u64 - cs) as usize;
        let site = &acc.sites[&gpos];
        let fwd_d = acc.depth_fwd[li] as u64;
        let rev_d = acc.depth_rev[li] as u64;
        let depth = fwd_d + rev_d;
        if depth == 0 {
            continue;
        }

        // Tally per strand. full_tally mode (N-ref / no fasta): site.cnt holds
        // every base; normal mode: mismatches only, matched = depth - mm.
        let ref_base = refseq.and_then(|s| base_idx(s[gpos as usize]));
        let full_tally = ref_base.is_none();
        let mm_total: u64 = site.cnt.iter().map(|&x| x as u64).sum::<u64>();

        let majority = {
            let mut tc = [0u64; 4];
            for (b, t) in tc.iter_mut().enumerate() {
                *t = site.cnt[b] as u64 + site.cnt[b + 4] as u64;
            }
            if !full_tally {
                let rb = ref_base.unwrap();
                tc[rb] += depth - mm_total;
            }
            tc.iter()
                .enumerate()
                .max_by_key(|(_, c)| *c)
                .map(|(i, _)| i)
                .unwrap_or(0)
        };
        let ref_idx = ref_base.unwrap_or(majority);

        let mut fwd_cnt = [0u64; 4];
        let mut rev_cnt = [0u64; 4];
        for b in 0..4 {
            fwd_cnt[b] = site.cnt[b] as u64;
            rev_cnt[b] = site.cnt[b + 4] as u64;
        }
        if !full_tally {
            fwd_cnt[ref_idx] += fwd_d - fwd_cnt.iter().sum::<u64>().min(fwd_d);
            rev_cnt[ref_idx] += rev_d - rev_cnt.iter().sum::<u64>().min(rev_d);
        }

        let total_cnt = [
            fwd_cnt[0] + rev_cnt[0],
            fwd_cnt[1] + rev_cnt[1],
            fwd_cnt[2] + rev_cnt[2],
            fwd_cnt[3] + rev_cnt[3],
        ];
        let var_reads = depth - total_cnt[ref_idx];
        if var_reads == 0 {
            continue; // pure-reference site is not a candidate
        }

        // A>I-consistent alt frequency per strand (legacy `ai` closure).
        let ai =
            |ref_i: usize, var_i: usize, cu_ref: usize, cu_var: usize, cnt: &[u64; 4], d: u64| {
                if d == 0 {
                    return 0.0;
                }
                let mut v = if ref_idx == ref_i {
                    cnt[var_i] as f64 / d as f64
                } else {
                    0.0
                };
                if enable_cu && ref_idx == cu_ref {
                    v = v.max(cnt[cu_var] as f64 / d as f64);
                }
                v
            };
        let fwd_ai = ai(0, 2, 1, 3, &fwd_cnt, fwd_d);
        let rev_ai = ai(3, 1, 2, 0, &rev_cnt, rev_d);

        let gtf_hit = gtf.map(|g| g.strands_at(chrom, gpos as i64));
        let (strand, evid_primary) = infer_strand(
            lib,
            fwd_d,
            rev_d,
            acc.junc_plus[li],
            acc.junc_minus[li],
            gtf_hit,
        );
        let edit_frac = match strand {
            StrandCall::Plus => fwd_ai,
            StrandCall::Minus => rev_ai,
            StrandCall::Amb => fwd_ai.max(rev_ai),
        };

        let var_fwd = fwd_d - fwd_cnt[ref_idx];
        let var_rev = rev_d - rev_cnt[ref_idx];
        let strand_bias =
            (var_fwd as f64 - var_rev as f64).abs() / (var_fwd + var_rev).max(1) as f64;

        let hp = refseq.map(|s| count::hp_len(s, gpos as usize)).unwrap_or(0);
        let jdist = count::nearest_dist(&jbounds, gpos as i64);
        let gnomad_af = match gnomad {
            Some(g) => g.af_at(chrom, gpos as i64)?,
            None => None,
        };

        let feats = ScoreFeatures {
            depth,
            edit_frac,
            mean_bq: acc.bq_sum[li] as f64 / depth as f64,
            mean_mapq: acc.mapq_sum[li] as f64 / depth as f64,
            strand_bias,
            gnomad_af,
            hp_len: hp,
            junction_dist: jdist.filter(|d| *d <= count::JUNCTION_EVID_DIST),
        };
        let score = spec.score(&feats);

        let mut evid = evid_primary.to_string();
        if let Some(d) = feats.junction_dist {
            evid.push_str(&format!(",JD{d}"));
        }
        if hp >= count::HP_MIN {
            evid.push_str(&format!(",HP{hp}"));
        }
        if gnomad_af.is_some() {
            evid.push_str(",GNOMAD");
        }
        if min_call_score.is_some_and(|t| score >= t) {
            evid.push_str(",MS");
        }

        out.push(Candidate {
            chrom: chrom.to_string(),
            pos0: gpos as i64,
            strand,
            evid,
            score,
            depth,
            var_freq: var_reads as f64 / depth as f64,
            fwd_freq: var_fwd as f64 / fwd_d.max(1) as f64,
            rev_freq: var_rev as f64 / rev_d.max(1) as f64,
        });
    }
    Ok(())
}

/// Engine entry: scatter + scoring/output contract.
pub fn run_call(params: &CallParams) -> Result<crate::CallStats, CallError> {
    let spec = CallSpec::load(params.spec.as_deref())?;
    let gtf = params
        .gtf
        .as_deref()
        .map(annot::GtfIndex::load)
        .transpose()?;
    let gnomad = params
        .gnomad
        .as_deref()
        .map(annot::GnomadIndex::load)
        .transpose()?;

    // .baln fast path — build the coordinate index in a single pass (Arc shared across block tasks); the contig list
    // and derived lengths come from the index; the BAM path is unchanged (both sources share the emit/score/output contract).
    let baln_idx: Option<std::sync::Arc<crate::baln::BalnIndex>> = params
        .baln
        .as_deref()
        .map(crate::baln::BalnReader::build_index)
        .transpose()?
        .map(std::sync::Arc::new);
    let contigs: Vec<(String, u32)> = match &baln_idx {
        Some(bi) => {
            // .baln derived length = coverage upper bound, which can truncate a terminal homopolymer run
            // (4/5.3M rows differ in HP annotation). Prefer true lengths from the fasta .fai when available;
            // with no fasta (hp always 0) or a missing contig, the derived value suffices.
            let fai = fai_lengths(params.fasta.as_deref());
            bi.contigs
                .iter()
                .zip(bi.tid_len.iter())
                .map(|(n, dl)| (n.clone(), fai.get(n).copied().unwrap_or(*dl)))
                .collect()
        }
        None => {
            let bam = rust_htslib::bam::IndexedReader::from_path(&params.bam)?;
            let h = bam.header();
            (0..h.target_count())
                .filter_map(|tid| {
                    let name = String::from_utf8_lossy(h.tid2name(tid)).into_owned();
                    h.target_len(tid).map(|l| (name, l as u32))
                })
                .collect()
        }
    };
    if let Some(g) = gnomad.as_ref() {
        let names: Vec<String> = contigs.iter().map(|(n, _)| n.clone()).collect();
        g.prepare(&names)?;
    }

    let pool = match params.threads {
        0 => rayon::ThreadPoolBuilder::new().build(),
        n => rayon::ThreadPoolBuilder::new().num_threads(n).build(),
    }
    .map_err(CallError::Pool)?;

    let seq_cache = SeqCache::default();
    const CHUNK: u64 = 32_000_000;
    let mut tasks: Vec<(usize, u64, u64)> = Vec::new();
    for (ci, (_, len)) in contigs.iter().enumerate() {
        let n_chunks = (*len as u64).div_ceil(CHUNK);
        for k in 0..n_chunks {
            tasks.push((ci, k * CHUNK, ((k + 1) * CHUNK).min(*len as u64)));
        }
    }

    let results: Vec<Result<Vec<Candidate>, CallError>> = pool.install(|| {
        tasks
            .par_iter()
            .map(|&(ci, cs, ce)| {
                let (chrom, len) = &contigs[ci];
                let t0 = std::time::Instant::now();
                let refseq = load_refseq(chrom, *len, params.fasta.as_deref(), &seq_cache);
                let t_ref = t0.elapsed();
                let t1 = std::time::Instant::now();
                let acc = if let Some(bi) = &baln_idx {
                    scatter::scatter_block_baln(
                        params.baln.as_ref().expect("baln set iff index built"),
                        &bi.idx,
                        bi.max_span,
                        ci as i32,
                        cs,
                        ce,
                        match &refseq {
                            Some(v) => v.as_slice(),
                            None => &[],
                        },
                        refseq.is_some(),
                    )?
                } else {
                    scatter::scatter_block(
                        &params.bam,
                        chrom,
                        cs,
                        ce,
                        match &refseq {
                            Some(v) => v.as_slice(),
                            None => &[],
                        },
                        refseq.is_some(),
                    )?
                };
                let t_scat = t1.elapsed();
                let t2 = std::time::Instant::now();
                let mut out = Vec::new();
                emit_chunk(
                    chrom,
                    cs,
                    &acc,
                    refseq.as_ref().map(|v| v.as_slice()),
                    params.lib,
                    params.enable_cu,
                    gtf.as_ref(),
                    gnomad.as_ref(),
                    &spec,
                    params.min_call_score,
                    &mut out,
                )?;
                let t_emit = t2.elapsed();
                T_REFSEQ.fetch_add(
                    t_ref.as_nanos() as u64,
                    std::sync::atomic::Ordering::Relaxed,
                );
                T_SCATTER.fetch_add(
                    t_scat.as_nanos() as u64,
                    std::sync::atomic::Ordering::Relaxed,
                );
                T_EMIT.fetch_add(
                    t_emit.as_nanos() as u64,
                    std::sync::atomic::Ordering::Relaxed,
                );
                Ok(out)
            })
            .collect()
    });

    eprintln!(
        "[scan timing] refseq={:.1}s scatter={:.1}s emit={:.1}s",
        T_REFSEQ.load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9,
        T_SCATTER.load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9,
        T_EMIT.load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9
    );
    let mut all: Vec<Candidate> = Vec::new();
    for r in results {
        all.extend(r?);
    }
    let n = all.len();
    let gnomad_hits = if gnomad.is_some() {
        all.iter().filter(|c| c.evid.contains("GNOMAD")).count()
    } else {
        0
    };
    out::write_bed(&params.out, &mut all)?;
    Ok(crate::CallStats {
        candidates: n,
        contigs: contigs.len(),
        gnomad_hits,
    })
}
