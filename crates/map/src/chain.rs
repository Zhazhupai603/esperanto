//! Anchor chaining (minimap2-style DP).
//!
//! Anchors are grouped by (contig, strand) and chained with
//! `f(i) = k + max_j[ f(j) + min(dq, dr) − gap_pen ]`, where j ranges over
//! earlier anchors with STRICTLY increasing rpos/qpos, `dr ≤ max_gap` (scan
//! breaks), `dq ≤ max_gap` (skip). DNA gap penalty is `dd / 2`; in RNA mode an
//! intron-shaped hole pays the flat `intron_penalty` — only when
//! `dr > dq && dd > min_intron && dq <= k+10` (big reference hole with
//! read-adjacent anchors; read-side holes pay dd/2 so the DP cannot earn
//! min(dq,dr) by skipping whole read segments).

use crate::seed::{Anchor, Strand};
use std::collections::HashSet;

/// Chaining parameters (see frozen-parameter table).
#[derive(Clone, Copy, Debug)]
pub struct ChainParams {
    /// Maximum gap in query or reference between chainable anchors.
    pub max_gap: i64,
    /// Minimum score for a chain to be reported.
    pub min_chain_score: i32,
    /// Anchor match reward (equals seed k).
    pub k: i32,
    /// RNA mode: enable intron-shaped gap penalty.
    pub rna: bool,
    /// Minimum intron length that qualifies as intron-shaped.
    pub min_intron: i64,
    /// Flat penalty for intron-shaped gaps (RNA only).
    pub intron_penalty: i64,
}

impl Default for ChainParams {
    fn default() -> ChainParams {
        ChainParams {
            max_gap: 1_000_000,
            min_chain_score: 40,
            k: 15,
            rna: false,
            min_intron: 20,
            intron_penalty: 8,
        }
    }
}

impl ChainParams {
    /// RNA defaults (frozen): max_gap 1M, min score 40, k 15, intron 20/8.
    pub fn rna_default() -> ChainParams {
        ChainParams {
            rna: true,
            ..ChainParams::default()
        }
    }
}

/// Maximum number of chain candidates advanced to the tie-break stage.
pub const TOP_CANDIDATES: usize = 24;

/// One chained candidate: anchors sorted by `(rpos, qpos)` ascending.
#[derive(Clone, Debug)]
pub struct Chain {
    /// Contig index.
    pub contig: u32,
    /// Chain strand (read-relative).
    pub strand: Strand,
    /// Anchors in `(rpos, qpos)` ascending order.
    pub anchors: Vec<Anchor>,
    /// DP score of the chain.
    pub score: i32,
}

impl Chain {
    /// Reference span covered by the chain: (first anchor rpos, last anchor rpos + k).
    pub fn ref_span(&self, k: u32) -> (u32, u32) {
        let first = self.anchors.first().map(|a| a.rpos).unwrap_or(0);
        let last = self.anchors.last().map(|a| a.rpos + k).unwrap_or(first);
        (first, last)
    }

    /// Read span covered by the chain: (first anchor qpos, last anchor qpos + k).
    pub fn read_span(&self, k: u32) -> (u32, u32) {
        let first = self.anchors.first().map(|a| a.qpos).unwrap_or(0);
        let last = self.anchors.last().map(|a| a.qpos + k).unwrap_or(first);
        (first, last)
    }
}

/// Chain anchors into candidates, highest score first.
///
/// Groups are processed in (contig, strand) order; results are merged and
/// sorted by score descending (ties: contig, strand, first rpos, first qpos —
/// byte-level reproducible). Chains scoring below `params.min_chain_score`
/// are dropped.
pub fn chain_anchors(mut anchors: Vec<Anchor>, params: &ChainParams) -> Vec<Chain> {
    anchors.sort_by_key(|a| (a.contig, a.strand as u8, a.rpos, a.qpos));

    let mut chains: Vec<Chain> = Vec::new();
    let mut group_start = 0usize;
    while group_start < anchors.len() {
        let mut group_end = group_start + 1;
        let key = (
            anchors[group_start].contig,
            anchors[group_start].strand as u8,
        );
        while group_end < anchors.len()
            && (anchors[group_end].contig, anchors[group_end].strand as u8) == key
        {
            group_end += 1;
        }
        chain_group(&anchors[group_start..group_end], params, &mut chains);
        group_start = group_end;
    }

    chains.retain(|c| c.score >= params.min_chain_score);
    chains.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.contig.cmp(&b.contig))
            .then((a.strand as u8).cmp(&(b.strand as u8)))
            .then(
                a.anchors
                    .first()
                    .map_or(0, |x| x.rpos)
                    .cmp(&b.anchors.first().map_or(0, |x| x.rpos)),
            )
            .then(
                a.anchors
                    .first()
                    .map_or(0, |x| x.qpos)
                    .cmp(&b.anchors.first().map_or(0, |x| x.qpos)),
            )
    });
    chains
}

/// Single-group DP chaining (emits all above-threshold, non-overlapping
/// chains, highest score first).
fn chain_group(group: &[Anchor], params: &ChainParams, out: &mut Vec<Chain>) {
    let n = group.len();
    if n == 0 {
        return;
    }
    let k = params.k;
    let mut f = vec![k; n];
    let mut parent = vec![usize::MAX; n];
    for i in 0..n {
        let a = group[i];
        // Look-back window: stop when the reference gap exceeds max_gap
        // (anchors are rpos-sorted within the group).
        let mut j = i;
        while j > 0 {
            j -= 1;
            let b = group[j];
            if a.rpos <= b.rpos {
                continue; // duplicate-position anchors never link
            }
            let dr = (a.rpos - b.rpos) as i64;
            if dr > params.max_gap {
                break;
            }
            if a.qpos <= b.qpos {
                continue; // read side must strictly increase too
            }
            let dq = (a.qpos - b.qpos) as i64;
            if dq > params.max_gap {
                continue;
            }
            let gain = dq.min(dr) as i32;
            let dd = (dr - dq).abs();
            // RNA: the flat intron penalty applies only to true junction
            // geometry (big ref hole, read-side anchors adjacent); read-side
            // holes pay dd/2 — otherwise the DP earns min(dq,dr) by skipping
            // whole read segments.
            let pen =
                if params.rna && dr > dq && dd > params.min_intron && dq <= params.k as i64 + 10 {
                    params.intron_penalty as i32
                } else {
                    (dd / 2) as i32
                };
            let cand = f[j] + gain - pen;
            if cand > f[i] {
                f[i] = cand;
                parent[i] = j;
            }
        }
    }

    // Backtrack: from the highest-scoring unused end, walk parents; a walk
    // hitting a used anchor DISCARDS the whole chain (marks nothing).
    let mut used = vec![false; n];
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| std::cmp::Reverse(f[i]));
    for &end in &order {
        if used[end] || f[end] < params.min_chain_score {
            continue;
        }
        let mut idx = Vec::new();
        let mut cur = end;
        let mut overlapped = false;
        loop {
            if used[cur] {
                overlapped = true;
                break;
            }
            idx.push(cur);
            if parent[cur] == usize::MAX {
                break;
            }
            cur = parent[cur];
        }
        if overlapped {
            continue;
        }
        for &i in &idx {
            used[i] = true;
        }
        idx.reverse();
        out.push(Chain {
            contig: group[end].contig,
            strand: group[end].strand,
            anchors: idx.iter().map(|&i| group[i]).collect(),
            score: f[end],
        });
    }
}

/// Best score among chains sharing no (qpos, rpos) anchor with `chains[0]`;
/// 0 if none. Chains are score-descending, so the first disjoint chain is
/// the maximum.
pub fn second_score(chains: &[Chain]) -> i32 {
    let Some(best) = chains.first() else {
        return 0;
    };
    let best_keys: HashSet<(u32, u32)> = best.anchors.iter().map(|a| (a.qpos, a.rpos)).collect();
    for c in chains.iter().skip(1) {
        if c.anchors
            .iter()
            .all(|a| !best_keys.contains(&(a.qpos, a.rpos)))
        {
            return c.score;
        }
    }
    0
}

/// All chains scoring exactly what `chains[0]` scores (multi-mapping check).
pub fn tied_best(chains: &[Chain]) -> Vec<&Chain> {
    if chains.is_empty() {
        return Vec::new();
    }
    let best = chains[0].score;
    chains.iter().filter(|c| c.score == best).collect()
}
