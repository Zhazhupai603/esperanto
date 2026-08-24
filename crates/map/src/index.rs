//! Minimizer index over a reference genome.
//!
//! Tables are `Box::leak`-ed at build time (or leaked from an mmap at load
//! time) so an `Index` is a read-only, shareable, zero-copy structure.

use crate::fasta::{Base, Reference};
use crate::seed::{minimizers_from_codes, Minimizer, SeedParams, Strand};
use rayon::prelude::*;

/// Share fraction `HIGH_FREQ_FRAC` (frozen parameter): a k-mer whose hit
/// count exceeds [`Index::freq_cutoff`] is treated as high-frequency.
pub const HIGH_FREQ_FRAC: f64 = 0.0002;

/// Position packing: `contig(31b) << 33 | strand(1b) << 32 | pos(32b)`.
#[inline]
pub fn pack_pos(contig: u32, strand: Strand, pos: u32) -> u64 {
    ((contig as u64) << 33) | ((strand as u64) << 32) | (pos as u64)
}

/// Inverse of [`pack_pos`]: `(contig, pos, strand)`.
#[inline]
pub fn unpack_pos(p: u64) -> (u32, u32, Strand) {
    let contig = (p >> 33) as u32;
    let pos = p as u32;
    let strand = if (p >> 32) & 1 == 1 {
        Strand::Minus
    } else {
        Strand::Plus
    };
    (contig, pos, strand)
}

impl Strand {
    /// XOR combination used to derive anchor strand from read/ref strands.
    #[inline]
    pub fn xor(a: Strand, b: Strand) -> Strand {
        if a == b {
            Strand::Plus
        } else {
            Strand::Minus
        }
    }
}

/// Read-only minimizer index.
pub struct Index {
    /// Seed parameters used at build time.
    pub params: SeedParams,
    /// Format version (paidx v1).
    pub version: u32,
    /// High-frequency cutoff: k-mers with `count > freq_cutoff` are skipped
    /// by anchor collection.
    pub freq_cutoff: u32,
    /// The reference these tables describe.
    pub reference: &'static Reference,
    /// Distinct canonical k-mers, ascending.
    pub kmers: &'static [u64],
    /// Start offset of each k-mer's position run in `positions`.
    pub offsets: &'static [u64],
    /// Length of each k-mer's position run.
    pub counts: &'static [u32],
    /// Packed positions, grouped by k-mer, ascending within each run.
    pub positions: &'static [u64],
}

/// Zero-allocation query result: a slice into the position table.
pub struct IndexHit {
    /// Number of reference occurrences of the queried k-mer.
    pub count: u32,
    /// Packed positions (ascending).
    pub positions: &'static [u64],
}

impl Index {
    /// Build an index from a parsed reference (per-contig rayon parallel).
    ///
    /// All tables and the reference itself are leaked and never freed.
    pub fn build(reference: Reference, params: SeedParams) -> Index {
        let k = params.k;
        let w = params.w;

        let per_contig: Vec<Vec<(u64, u64)>> = reference
            .contigs
            .par_iter()
            .enumerate()
            .map(|(ci, contig)| {
                let codes = contig_codes(contig);
                let mins = minimizers_from_codes(&codes, k, w);
                mins.into_iter()
                    .map(|m| (m.kmer, pack_pos(ci as u32, m.strand, m.pos)))
                    .collect()
            })
            .collect();

        let mut all: Vec<(u64, u64)> = per_contig.into_iter().flatten().collect();
        all.sort_unstable();
        all.dedup();


        let mut kmers: Vec<u64> = Vec::new();
        let mut offsets: Vec<u64> = Vec::new();
        let mut counts: Vec<u32> = Vec::new();
        let mut positions: Vec<u64> = Vec::new();
        let mut i = 0usize;
        while i < all.len() {
            let kmer = all[i].0;
            let start = i;
            while i < all.len() && all[i].0 == kmer {
                positions.push(all[i].1);
                i += 1;
            }
            kmers.push(kmer);
            offsets.push(start as u64);
            counts.push((i - start) as u32);
        }

        let freq_cutoff = compute_cutoff(&counts);

        let reference: &'static Reference = Box::leak(Box::new(reference));
        let kmers: &'static [u64] = Box::leak(kmers.into_boxed_slice());
        let offsets: &'static [u64] = Box::leak(offsets.into_boxed_slice());
        let counts: &'static [u32] = Box::leak(counts.into_boxed_slice());
        let positions: &'static [u64] = Box::leak(positions.into_boxed_slice());

        Index {
            params,
            version: 1,
            freq_cutoff,
            reference,
            kmers,
            offsets,
            counts,
            positions,
        }
    }

    /// Query a canonical k-mer; binary search, zero allocation.
    pub fn query(&self, kmer: u64) -> IndexHit {
        let pos = self.kmers.binary_search(&kmer);
        match pos {
            Ok(i) => {
                let start = self.offsets[i] as usize;
                let count = self.counts[i] as usize;
                IndexHit {
                    count: self.counts[i],
                    positions: &self.positions[start..start + count],
                }
            }
            Err(_) => IndexHit {
                count: 0,
                positions: &[],
            },
        }
    }

    /// Whether `count` marks a k-mer as high-frequency (skip in anchoring).
    pub fn is_high_freq(&self, count: u32) -> bool {
        count > self.freq_cutoff
    }
}

/// Frequency cutoff from per-k-mer counts: sort descending, take the value at
/// rank `len * HIGH_FREQ_FRAC` (capped at `len - 1`), minimum 1. Empty input
/// yields `u32::MAX` (nothing is high-frequency).
pub fn compute_cutoff(counts: &[u32]) -> u32 {
    if counts.is_empty() {
        return u32::MAX;
    }
    let mut sorted: Vec<u32> = counts.to_vec();
    sorted.sort_unstable_by(|a, b| b.cmp(a));
    let rank = ((sorted.len() as f64) * HIGH_FREQ_FRAC) as usize;
    let rank = rank.min(sorted.len() - 1);
    sorted[rank].max(1)
}

fn contig_codes(contig: &crate::fasta::Contig) -> Vec<u64> {
    const N_CODE: u64 = 0xFF;
    contig
        .slice_ascii(0, contig.len)
        .into_iter()
        .map(|b| match Base::from_ascii(b).code() {
            Some(c) => c as u64,
            None => N_CODE,
        })
        .collect()
}

/// qpos for an anchor, keyed on the ANCHOR strand (not the minimizer's):
/// Plus → `m.pos`, Minus → `read_len − m.pos − k`.
#[inline]
fn anchor_qpos(astrand: Strand, m: &Minimizer, read_len: u32, k: u32) -> u32 {
    if astrand == Strand::Plus {
        m.pos
    } else {
        read_len - m.pos - k
    }
}

/// Collect anchors for plain read minimizers.
///
/// Skips k-mers with `count == 0` or `count > occ_cap` (all positions of a
/// passing k-mer are emitted). Anchor strand is the XOR of the read and
/// reference minimizer strands (same → Plus, different → Minus), and qpos is
/// derived from the anchor strand: Plus anchor → `m.pos`, Minus anchor →
/// `read_len − m.pos − k`. Returns the anchors and the number of read
/// minimizers that passed the cap.
pub fn collect_anchors(
    index: &Index,
    read_mins: &[Minimizer],
    read_len: u32,
    k: u32,
    occ_cap: u32,
) -> (Vec<crate::seed::Anchor>, usize) {
    let mut anchors = Vec::new();
    let mut seeds_hit = 0usize;
    for m in read_mins {
        let hit = index.query(m.kmer);
        if hit.count == 0 || hit.count > occ_cap {
            continue;
        }
        seeds_hit += 1;
        for &p in hit.positions {
            let (contig, rpos, rstrand) = unpack_pos(p);
            let astrand = Strand::xor(m.strand, rstrand);
            anchors.push(crate::seed::Anchor {
                qpos: anchor_qpos(astrand, m, read_len, k),
                rpos,
                contig,
                strand: astrand,
            });
        }
    }
    (anchors, seeds_hit)
}

/// Reverse-complement a 2-bit k-mer of length `k`.
fn revcomp_kmer(kmer: u64, k: u32) -> u64 {
    let mut out: u64 = 0;
    for i in 0..k {
        let base = (kmer >> (2 * i)) & 3;
        out |= (3 - base) << (2 * (k - 1 - i));
    }
    out
}

/// Canonical form (min of k-mer and revcomp) with orientation flag.
fn canonical(kmer: u64, k: u32) -> (u64, Strand) {
    let rc = revcomp_kmer(kmer, k);
    if kmer <= rc {
        (kmer, Strand::Plus)
    } else {
        (rc, Strand::Minus)
    }
}

/// Emit anchors for one hitting variant round.
///
/// `take` caps the positions consumed (64 for single-site, 32 for double);
/// `rev` selects the round: fwd round anchors carry the reference strand,
/// rev round anchors carry its flip. qpos always follows the anchor strand.
fn emit_variant_round(
    anchors: &mut Vec<crate::seed::Anchor>,
    positions: &[u64],
    m: &Minimizer,
    read_len: u32,
    k: u32,
    take: usize,
    rev: bool,
) {
    for &p in positions.iter().take(take) {
        let (contig, rpos, rstrand) = unpack_pos(p);
        let astrand = if rev { rstrand.flip() } else { rstrand };
        anchors.push(crate::seed::Anchor {
            qpos: anchor_qpos(astrand, m, read_len, k),
            rpos,
            contig,
            strand: astrand,
        });
    }
}

/// Single-site G→A / C→T edit variants (editing-aware anchor collection).
///
/// For each read minimizer, the k-mer is packed in fwd orientation (site i
/// lives at 2-bit slot `k − 1 − i`); every editable site runs TWO rounds —
/// fwd-orientation variant (want = alt) and rev-orientation variant (mirrored
/// site, want = 3 − alt). Both rounds canonicalize to the same key and BOTH
/// query and BOTH emit: a hitting site contributes its first 64 positions as
/// anchors with the reference strand (fwd round), then the same first 64 as
/// anchors with the flipped strand (rev round), so every hit position yields
/// exactly one `(m.pos, Plus)` and one `(read_len − m.pos − k, Minus)` anchor
/// Single-site G→A / C→T edit variants (editing-aware anchor collection).
///
/// Per editable site i (read-orientation k-mer base G/C → alt A/T), TWO
/// rounds run, each producing a DIFFERENT variant key:
/// - fwd round: edit slot (k−1−i) of the read-forward packing ← alt (edits
///   read position i); anchors carry the reference strand;
/// - rev round: edit slot (k−1−i) of the rev packing ← 3−alt (edits the
///   MIRROR read position k−1−i); skipped when the rev slot already holds
///   3−alt (no-op: mirror read base already equals alt); anchors carry the
///   flipped reference strand.
///
/// qpos follows the anchor strand (Plus → m.pos, Minus → read_len−m.pos−k).
/// Each hitting round takes the first 64 positions and increments var_hits.
pub fn collect_anchors_edit_variants(
    index: &Index,
    read_mins: &[Minimizer],
    read_len: u32,
    k: u32,
    occ_cap: u32,
) -> (Vec<crate::seed::Anchor>, usize) {
    let mut anchors = Vec::new();
    let mut var_hits = 0usize;
    for m in read_mins {
        // k-mer as it reads on the read's forward strand
        let fwd = if m.strand == Strand::Minus {
            revcomp_kmer(m.kmer, k)
        } else {
            m.kmer
        };
        let rev = revcomp_kmer(fwd, k);
        for i in 0..k as usize {
            let slot = (k as usize - 1) - i;
            let base = (fwd >> (2 * slot)) & 3;
            let alt = match base {
                2 => 0u64, // G -> A
                1 => 3u64, // C -> T
                _ => continue,
            };
            // fwd round
            let v = (fwd & !(3 << (2 * slot))) | (alt << (2 * slot));
            let key = canonical(v, k).0;
            let hit = index.query(key);
            if hit.count > 0 && hit.count <= occ_cap {
                var_hits += 1;
                emit_variant_round(&mut anchors, hit.positions, m, read_len, k, 64, false);
            }
            // rev round (mirror-position edit; skip no-ops)
            let want = 3 - alt;
            let cur = (rev >> (2 * slot)) & 3;
            if cur == want {
                continue;
            }
            let v = (rev & !(3 << (2 * slot))) | (want << (2 * slot));
            let key = canonical(v, k).0;
            let hit = index.query(key);
            if hit.count > 0 && hit.count <= occ_cap {
                var_hits += 1;
                emit_variant_round(&mut anchors, hit.positions, m, read_len, k, 64, true);
            }
        }
    }
    (anchors, var_hits)
}

/// Double-site edit variants (deep fallback): pairs of editable sites
/// (i < j, each G→A / C→T), two orientation rounds per pair, first 32
/// positions per round, no `var_hits` accounting, no deduplication (the same
/// canonical key produced by a different site pair emits again).
/// Double-site edit variants (deep fallback): pairs of editable sites
/// (i < j, each G→A / C→T). Per orientation round both slots of that
/// packing are replaced (rev round carries 3−alt at both, editing the
/// mirror read positions); no no-op skip; first 32 positions per hitting
/// round; no var_hits accounting; no deduplication.
pub fn collect_anchors_edit_variants_double(
    index: &Index,
    read_mins: &[Minimizer],
    read_len: u32,
    k: u32,
    occ_cap: u32,
) -> Vec<crate::seed::Anchor> {
    let mut anchors = Vec::new();
    for m in read_mins {
        let fwd = if m.strand == Strand::Minus {
            revcomp_kmer(m.kmer, k)
        } else {
            m.kmer
        };
        let rev = revcomp_kmer(fwd, k);
        // editable sites in left-to-right (i) order: (slot, alt)
        let sites: Vec<(usize, u64)> = (0..k as usize)
            .filter_map(|i| {
                let slot = (k as usize - 1) - i;
                match (fwd >> (2 * slot)) & 3 {
                    2 => Some((slot, 0u64)), // G -> A
                    1 => Some((slot, 3u64)), // C -> T
                    _ => None,
                }
            })
            .collect();
        for a in 0..sites.len() {
            for b in (a + 1)..sites.len() {
                let (pa, wa) = sites[a];
                let (pb, wb) = sites[b];
                for (packing, rev_round) in [(fwd, false), (rev, true)] {
                    let (wa2, wb2) = if rev_round {
                        (3 - wa, 3 - wb)
                    } else {
                        (wa, wb)
                    };
                    let mut v = packing;
                    v = (v & !(3 << (2 * pa))) | (wa2 << (2 * pa));
                    v = (v & !(3 << (2 * pb))) | (wb2 << (2 * pb));
                    let key = canonical(v, k).0;
                    let hit = index.query(key);
                    if hit.count == 0 || hit.count > occ_cap {
                        continue;
                    }
                    emit_variant_round(&mut anchors, hit.positions, m, read_len, k, 32, rev_round);
                }
            }
        }
    }
    anchors
}
