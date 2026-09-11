//! Deterministic multi-isoform attribution: collapse a read's candidate
//! projections into a single placement with MAPQ.

use crate::cigar::{cigar_string, CigarOp};
use crate::junction::JunctionSet;

/// One candidate placement of a read: a transcript id plus its projection.
/// Field types match [`crate::TxMap::project`] output (`contig_id` indexes
/// `TxMap::contigs()`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AttributionCandidate {
    pub tx_id: u32,
    pub contig_id: u32,
    pub pos: u32,
    pub cigar: Vec<CigarOp>,
}

/// The unique winning placement for a read.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Placement {
    pub contig_id: u32,
    pub pos: u32,
    pub cigar: Vec<CigarOp>,
    /// Mapping quality, 0..=60.
    pub mapq: u8,
    /// Sorted, deduplicated transcript ids that support this placement.
    pub tx_ids: Vec<u32>,
}

/// Attribution result: `placement` is `None` iff the candidate input was empty.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Attribution {
    pub placement: Option<Placement>,
}

/// One merged group of candidates sharing an identical projection
/// `(contig_id, pos, cigar)`; `members` holds the candidate tx_ids with
/// multiplicity.
struct Group<'a> {
    contig_id: u32,
    pos: u32,
    cigar: &'a [CigarOp],
    members: Vec<u32>,
    score: u32,
}

/// Attribute a read to a single placement. Rules, strictly in order:
///
/// 1. Merge candidates with identical projections `(contig, pos, cigar)`
///    into groups; group members are the candidate tx_id lists.
/// 2. Rank groups by junction support: every `RefSkip` interval of the
///    group's CIGAR (placed at `pos`) that hits `junctions` scores 1.
/// 3. Highest score wins. On a tie the winner is the lexicographically first
///    `(contig name, pos, cigar string)` and `mapq = 0`, guaranteeing
///    reproducibility.
///
/// MAPQ (non-tied): `w = |winner members| / Σ|all members|`,
/// `MAPQ = min(60, floor(-10·log10(1 - w)))`, with `w >= 1 → 60` and
/// `w <= 0 → 0`. `tx_ids` is sorted and deduplicated before output.
pub fn attribute(
    candidates: &[AttributionCandidate],
    junctions: &JunctionSet,
    contig_names: &[String],
) -> Attribution {
    if candidates.is_empty() {
        return Attribution { placement: None };
    }

    // Rule 1: merge identical projections (first-appearance scan; the final
    // winner selection below re-orders deterministically, so scan order here
    // does not affect output).
    let mut groups: Vec<Group> = Vec::new();
    for cand in candidates {
        if let Some(group) = groups.iter_mut().find(|g| {
            g.contig_id == cand.contig_id && g.pos == cand.pos && g.cigar == cand.cigar.as_slice()
        }) {
            group.members.push(cand.tx_id);
        } else {
            groups.push(Group {
                contig_id: cand.contig_id,
                pos: cand.pos,
                cigar: &cand.cigar,
                members: vec![cand.tx_id],
                score: 0,
            });
        }
    }

    // Rule 2: junction-support score per group.
    for group in &mut groups {
        group.score = junctions.support_score(group.contig_id, group.pos, group.cigar);
    }

    // Rule 3: highest score wins; ties broken lexicographically.
    let best_score = groups.iter().map(|g| g.score).max().unwrap_or(0);
    let winners: Vec<&Group> = groups.iter().filter(|g| g.score == best_score).collect();
    let tied = winners.len() > 1;
    let winner = winners
        .into_iter()
        .min_by(|a, b| sort_key(a, contig_names).cmp(&sort_key(b, contig_names)))
        .expect("winners is non-empty by construction");

    let w = winner.members.len() as f64 / candidates.len() as f64;
    let mapq = if tied { 0 } else { mapq_from_weight(w) };

    let mut tx_ids = winner.members.clone();
    tx_ids.sort_unstable();
    tx_ids.dedup();

    Attribution {
        placement: Some(Placement {
            contig_id: winner.contig_id,
            pos: winner.pos,
            cigar: winner.cigar.to_vec(),
            mapq,
            tx_ids,
        }),
    }
}

/// Deterministic total-order key for tie-breaking:
/// `(contig name, pos, cigar string)` compared lexicographically.
fn sort_key(group: &Group<'_>, contig_names: &[String]) -> (String, u32, String) {
    (
        contig_names
            .get(group.contig_id as usize)
            .map_or(String::new(), |s| s.clone()),
        group.pos,
        cigar_string(group.cigar),
    )
}

/// MAPQ from winner weight `w = |winner members| / |all members|`:
/// `min(60, floor(-10·log10(1 - w)))`, with `w >= 1 → 60`, `w <= 0 → 0`.
fn mapq_from_weight(w: f64) -> u8 {
    if w >= 1.0 {
        return 60;
    }
    if w <= 0.0 {
        return 0;
    }
    let raw = -10.0 * (1.0 - w).log10();
    if !raw.is_finite() {
        return 0;
    }
    let floored = raw.floor();
    if floored <= 0.0 {
        return 0;
    }
    floored.min(60.0) as u8
}
