//! Scan driver: per (contig, block) pileup -> candidate assembly -> bed.

use crate::error::CallError;
use crate::filters;
use crate::out::{self, Candidate};
use crate::scatter::{self, BlockAcc};
use crate::snp::SnpSet;
use crate::{CallParams, CallStats, Channel};
use rayon::prelude::*;
use rust_htslib::bam::Read as BamRead;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

const BLOCK: u64 = 32_000_000;

/// Contig list from the BAM header (name, length).
fn bam_contigs(bam: &Path) -> Result<Vec<(String, u64)>, CallError> {
    let reader = rust_htslib::bam::IndexedReader::from_path(bam)?;
    let h = reader.header();
    Ok((0..h.target_count())
        .filter_map(|tid| {
            let name = String::from_utf8_lossy(h.tid2name(tid)).into_owned();
            h.target_len(tid).map(|l| (name, l))
        })
        .collect())
}

/// Contig names present in the fasta index.
fn fai_contigs(fasta: &Path) -> Result<std::collections::HashSet<String>, CallError> {
    let fai = format!("{}.fai", fasta.display());
    let text = std::fs::read_to_string(&fai).map_err(|e| CallError::Io {
        path: fai.clone(),
        source: e,
    })?;
    Ok(text
        .lines()
        .filter_map(|l| l.split('\t').next().map(str::to_string))
        .collect())
}

fn load_contig(fasta: &Path, chrom: &str, len: u64) -> Result<Arc<Vec<u8>>, CallError> {
    let fa = rust_htslib::faidx::Reader::from_path(fasta)?;
    let s = fa
        .fetch_seq(chrom, 0, (len as usize).saturating_sub(1))
        .map_err(|e| CallError::Fasta {
            path: fasta.display().to_string(),
            msg: format!("fetch {chrom}: {e}"),
        })?;
    Ok(Arc::new(s.iter().map(|b| b.to_ascii_uppercase()).collect()))
}

/// Assemble candidates for one block in one mismatch channel.
///
/// Decoy contigs never emit: their reads are already explained (the source
/// is absent from the primary reference), and any column called there would
/// be a placement artifact, not editing.
#[allow(clippy::too_many_arguments)]
fn emit_block(
    chrom: &str,
    acc: &BlockAcc,
    refseq: &[u8],
    snp: Option<&SnpSet>,
    p: &CallParams,
    stats: &BlockStats,
    channel: Channel,
    out: &mut Vec<Candidate>,
) {
    if chrom.starts_with("decoy_") {
        return;
    }
    let mut positions: Vec<u32> = acc.mism.keys().copied().collect();
    positions.sort_unstable();
    for gpos in positions {
        let rb = refseq[gpos as usize];
        let Some(alt_i) = channel.alt_index(rb) else {
            continue;
        };
        let ref_i = match rb {
            b'A' => 0,
            b'T' => 3,
            _ => continue,
        };
        let m = acc.mism[&gpos];
        let depth = acc.depth[(gpos as u64 - acc.cs) as usize];
        let n_alt = m[alt_i];
        if n_alt == 0 {
            continue;
        }
        let freq = f64::from(n_alt) / f64::from(depth);
        if depth < p.min_depth || n_alt < p.min_edit_reads || freq < p.min_edit_freq {
            if channel == Channel::Edit {
                stats.below_threshold.fetch_add(1, Ordering::Relaxed);
            }
            continue;
        }
        if snp.is_some_and(|s| s.contains(chrom, gpos)) {
            if channel == Channel::Edit {
                stats.snp_removed.fetch_add(1, Ordering::Relaxed);
            }
            continue;
        }
        let mut counts = m;
        counts[ref_i] = depth - m.iter().sum::<u32>();
        out.push(Candidate {
            chrom: chrom.to_string(),
            pos0: gpos,
            ref_base: rb,
            depth,
            counts,
            edit_freq: freq,
            filters: String::new(),
            channel,
        });
    }
}

#[derive(Default)]
struct BlockStats {
    below_threshold: AtomicUsize,
    snp_removed: AtomicUsize,
}

pub fn run_call(params: &CallParams) -> Result<CallStats, CallError> {
    let snp = params.snp.as_deref().map(SnpSet::load).transpose()?;
    let contigs = bam_contigs(&params.bam)?;
    let in_fasta = fai_contigs(&params.fasta)?;
    let pool = match params.threads {
        0 => rayon::ThreadPoolBuilder::new().build(),
        n => rayon::ThreadPoolBuilder::new().num_threads(n).build(),
    }
    .map_err(CallError::Pool)?;
    let stats = BlockStats::default();
    let filter_set = filters::FilterSet::load(
        &params
            .fasta
            .parent()
            .map(|p| p.to_path_buf())
            .unwrap_or_else(|| std::path::Path::new(".").to_path_buf()),
    );

    // One task per (contig, block); the contig sequence is loaded once per
    // contig and shared by its blocks.
    let per_contig: Vec<Result<Vec<Candidate>, CallError>> = pool.install(|| {
        contigs
            .par_iter()
            .filter(|(name, _)| in_fasta.contains(name))
            .map(|(chrom, len)| {
                let refseq = load_contig(&params.fasta, chrom, *len)?;
                let blocks: Vec<(u64, u64)> = (0..len.div_ceil(BLOCK))
                    .map(|k| (k * BLOCK, ((k + 1) * BLOCK).min(*len)))
                    .collect();
                let per_block: Vec<Result<Vec<Candidate>, CallError>> = blocks
                    .par_iter()
                    .map(|&(cs, ce)| {
                        let acc = scatter::scatter_block(
                            &params.bam,
                            chrom,
                            cs,
                            ce,
                            &refseq,
                            params.min_mapq,
                            params.min_bq,
                        )?;
                        let mut out = Vec::new();
                        for ch in [Channel::Edit, Channel::Ctrl1, Channel::Ctrl2] {
                            emit_block(
                                chrom, &acc, &refseq, snp.as_ref(), params, &stats, ch, &mut out,
                            );
                        }
                        Ok(out)
                    })
                    .collect();
                let mut all = Vec::new();
                for b in per_block {
                    all.extend(b?);
                }
                // Standard evidence filters (visible tags, never silent
                // drops) — annotated here where the contig sequence lives.
                for c in &mut all {
                    // The site sits at index 8 of the 17 bp context window;
                    // near the contig start the window clamps and the center
                    // follows it (was hardcoded 8, pointing at the wrong base).
                    let center = (c.pos0 as usize).min(8);
                    let ctx_start = c.pos0.saturating_sub(8) as usize;
                    let ctx_end = ((c.pos0 + 9) as usize).min(refseq.len());
                    let ctx = refseq.get(ctx_start..ctx_end).unwrap_or(&[]);
                    let f = filter_set.site(chrom, c.pos0, ctx, center);
                    c.filters = f.tags().join(",");
                }
                Ok(all)
            })
            .collect()
    });

    let mut edit_cands: Vec<Candidate> = Vec::new();
    let mut control_cands: Vec<Candidate> = Vec::new();
    let mut channel_stats = [crate::ChannelStat::default(); 3];
    for c in per_contig {
        for cand in c? {
            // Channel accounting rides the filter-clean subset, matching the
            // editing channel's own candidate contract.
            if cand.filters.is_empty() {
                let ci = match cand.channel {
                    Channel::Edit => 0,
                    Channel::Ctrl1 => 1,
                    Channel::Ctrl2 => 2,
                };
                channel_stats[ci].sites += 1;
                if let Some(ai) = cand.channel.alt_index(cand.ref_base) {
                    channel_stats[ci].alt_reads += u64::from(cand.counts[ai]);
                }
            }
            match cand.channel {
                Channel::Edit => edit_cands.push(cand),
                _ => control_cands.push(cand),
            }
        }
    }
    let mut cands = edit_cands;
    cands.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.pos0.cmp(&b.pos0)));
    let parent = params.out.parent().map(|p| p.to_path_buf()).unwrap_or_default();
    out::write_control_bed(&parent.join("control_candidates.bed"), &mut control_cands)?;
    let (e, c1, c2) = (channel_stats[0], channel_stats[1], channel_stats[2]);
    let ctrl_mean = (c1.sites + c2.sites) as f64 / 2.0;
    let snr = if ctrl_mean > 0.0 {
        e.sites as f64 / ctrl_mean
    } else {
        f64::INFINITY
    };
    let implied_fdr = if e.sites > 0 {
        ctrl_mean / e.sites as f64
    } else {
        0.0
    };
    // f64::INFINITY formats as `inf`, which is not valid JSON and crashed
    // the downstream verify/score parse whenever a control channel was
    // empty. Emit null for the unbounded case instead.
    let snr_json = if snr.is_finite() {
        format!("{snr:.4}")
    } else {
        "null".to_string()
    };
    let channels_doc = format!(
        "{{\n  \"editing\": {{\"sites\": {}, \"alt_reads\": {}}},\n  \"control_1\": {{\"sites\": {}, \"alt_reads\": {}}},\n  \"control_2\": {{\"sites\": {}, \"alt_reads\": {}}},\n  \"snr\": {},\n  \"implied_fdr\": {:.6}\n}}\n",
        e.sites, e.alt_reads, c1.sites, c1.alt_reads, c2.sites, c2.alt_reads, snr_json, implied_fdr
    );
    let channels_path = parent.join("channels.json");
    if let Err(e2) = std::fs::write(&channels_path, channels_doc) {
        return Err(CallError::Io {
            path: channels_path.display().to_string(),
            source: std::io::Error::other(e2.to_string()),
        });
    }

    // Filter bookkeeping: per-tag counts land next to the bed for the run
    // report, and the tags themselves ride the bed's trailing column (the
    // verify stage appends after them; the VCF stage merges everything).
    let mut counts: std::collections::BTreeMap<&str, u64> = std::collections::BTreeMap::new();
    for c in &cands {
        for tag in c.filters.split(',').filter(|t| !t.is_empty()) {
            *counts.entry(tag).or_insert(0) += 1;
        }
    }
    let json_path = params.out.with_extension("filters.json");
    let mut doc = String::from("{\n");
    for (i, (k, v)) in counts.iter().enumerate() {
        let comma = if i + 1 < counts.len() { "," } else { "" };
        doc.push_str(&format!("  \"{k}\": {v}{comma}\n"));
    }
    doc.push_str("}\n");
    if let Err(e) = std::fs::write(&json_path, doc) {
        return Err(CallError::Io {
            path: json_path.display().to_string(),
            source: std::io::Error::other(e.to_string()),
        });
    }
    out::write_bed(&params.out, &mut cands)?;
    Ok(CallStats {
        candidates: cands.len(),
        below_threshold: stats.below_threshold.load(Ordering::Relaxed),
        snp_removed: stats.snp_removed.load(Ordering::Relaxed),
        channels: channel_stats,
    })
}

/// Reference base at `pos0`, used by tests and by the VCF writer through the
/// bed contract; `N` when off the contig.
pub fn ref_base(seq: &[u8], pos0: usize) -> u8 {
    seq.get(pos0).copied().map(|b| b.to_ascii_uppercase()).unwrap_or(b'N')
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn params() -> CallParams {
        CallParams::standard("x.bam".into(), "o.bed".into(), "r.fa".into())
    }

    fn acc_with(cs: u64, depth: &[u32], mism: &[(u32, [u32; 4])]) -> BlockAcc {
        BlockAcc {
            cs,
            depth: depth.to_vec(),
            mism: mism.iter().copied().collect::<HashMap<_, _>>(),
        }
    }

    #[test]
    fn emits_a_to_g_and_t_to_c_only() {
        let refseq = b"ATCG".to_vec();
        let acc = acc_with(
            0,
            &[20, 20, 20, 20],
            &[
                (0, [0, 0, 5, 0]), // A>G: candidate
                (1, [0, 5, 0, 0]), // T>C: candidate
                (2, [0, 0, 0, 5]), // C>T: not A-to-I
                (3, [5, 0, 0, 0]), // G>A: not A-to-I
            ],
        );
        let stats = BlockStats::default();
        let mut out = Vec::new();
        emit_block("chr1", &acc, &refseq, None, &params(), &stats, Channel::Edit, &mut out);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].counts, [15, 0, 5, 0]);
        assert_eq!(out[1].counts, [0, 5, 0, 15]);
        assert!((out[0].edit_freq - 0.25).abs() < 1e-12);
    }

    #[test]
    fn applies_thresholds_and_counts_them() {
        let refseq = b"AAA".to_vec();
        let acc = acc_with(
            0,
            &[9, 200, 20],
            &[
                (0, [0, 0, 5, 0]), // depth 9 < 10
                (1, [0, 0, 1, 0]), // 1 read < 2
                (2, [0, 0, 2, 0]), // 2/20 = 10%: kept
            ],
        );
        let stats = BlockStats::default();
        let mut out = Vec::new();
        emit_block("chr1", &acc, &refseq, None, &params(), &stats, Channel::Edit, &mut out);
        assert_eq!(out.len(), 1);
        assert_eq!(stats.below_threshold.load(Ordering::Relaxed), 2);
    }

    #[test]
    fn a_site_with_only_non_editing_mismatch_is_not_a_candidate() {
        let refseq = b"A".to_vec();
        let acc = acc_with(0, &[20], &[(0, [0, 5, 0, 5])]); // A>C and A>T only
        let stats = BlockStats::default();
        let mut out = Vec::new();
        emit_block("chr1", &acc, &refseq, None, &params(), &stats, Channel::Edit, &mut out);
        assert!(out.is_empty());
        assert_eq!(stats.below_threshold.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn control_channels_count_the_same_thresholds() {
        let refseq = b"AT".to_vec();
        let acc = acc_with(
            0,
            &[20, 20],
            &[
                (0, [0, 3, 5, 4]), // A: C=3, G=5, T=4 -> E/C1/C2 all pass
                (1, [4, 5, 3, 0]), // T: A=4, C=5, G=3 -> all pass
            ],
        );
        let stats = BlockStats::default();
        for (ch, alt_idx) in [(Channel::Edit, [2usize, 1usize]), (Channel::Ctrl1, [1, 2]), (Channel::Ctrl2, [3, 0])] {
            let mut out = Vec::new();
            emit_block("chr1", &acc, &refseq, None, &params(), &stats, ch, &mut out);
            assert_eq!(out.len(), 2, "{ch:?}");
            assert_eq!(out[0].counts[alt_idx[0]], acc.mism[&0][alt_idx[0]]);
            assert_eq!(out[1].counts[alt_idx[1]], acc.mism[&1][alt_idx[1]]);
            assert!(out.iter().all(|c| c.channel == ch));
        }
        // below-threshold accounting stays on the editing channel only.
        assert_eq!(stats.below_threshold.load(Ordering::Relaxed), 0);
    }
}
