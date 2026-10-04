//! Mapping quality.
//!
//! [`mapq_of`] uses minimap2's alignment-stage MAPQ expression for short
//! spliced reads (`mm_set_mapq2`, Li, Bioinformatics 34:3094, 2018):
//!
//! ```text
//! mapq = identity × pen_cm × 40 × (1 − x²) × ln(score / match_score)
//! pen_cm = min(pen_s1, if cnt > 10 { 1 } else { 0.1 × cnt })
//! pen_s1 = if score > 100 { 1 } else { 0.01 × score }
//! x      = second_score / score
//! ```
//!
//! Two terms of the original are deliberately absent, since esperanto tracks
//! neither input: `uniq_ratio = sum_sc / (sum_sc + rep_len)`, which needs the
//! read's repeat-covered length, and `− 4.343 ln(n_sub + 1)`, which needs a
//! count of sub-optimal placements. Both would only lower MAPQ; the `(1 − x²)`
//! margin is what guards repeats here, since a repeat copy raises
//! `second_score` towards `score` and drives MAPQ to 0 on its own.
//!
//! `best` and `second` must come from the *same* scoring domain. Chain DP
//! scores and extension DP scores are not comparable: a chain score is the
//! anchor-covered span, whereas extension scores include base-level
//! substitutions and gaps. Mixing domains would distort the score ratio.
//! [`margin`] therefore picks a matched pair, and prefers the extension one:
//! minimap2 takes the same branch for short spliced reads (`is_sr &&
//! is_splice ⇒ x = dp_max2 / dp_max`) because an unspliced chain outscores the
//! correct spliced one often enough that the chain ratio misleads.
//!
//! This module also hosts the crate-wide [`ReadAlignment`] record (spec:
//! shared core type), since MAPQ is its first consumer.

use crate::extend::CigarOp;
use crate::gtf::RefinedJunction;
use crate::seed::Strand;

/// Which layer produced a placement.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Origin {
    /// The genomic layer (seeding/chaining/extension against the genome).
    #[default]
    Genomic,
    /// The transcriptome-first engine's projection (candidate generator).
    L1,
    /// A decoy contig in the reference contract.
    Decoy,
}

impl Origin {
    /// Lower sorts first: on equal strict scores the genomic layer's
    /// candidate wins (its CIGAR was computed against the genome directly,
    /// with no projection error).
    pub fn rank(self) -> u8 {
        match self {
            Origin::Genomic => 0,
            Origin::L1 => 1,
            Origin::Decoy => 2,
        }
    }

    /// Short name for the `EK` evidence tag.
    pub fn label(self) -> &'static str {
        match self {
            Origin::L1 => "l1",
            Origin::Genomic => "genomic",
            Origin::Decoy => "decoy",
        }
    }
}

/// A fully placed read alignment (field semantics defined by the spec).
#[derive(Clone, Debug, Default)]
pub struct ReadAlignment {
    /// Contig index.
    pub contig: u32,
    /// 0-based inclusive leftmost reference position.
    pub pos: u32,
    /// Alignment strand.
    pub strand: Strand,
    /// Extension score.
    pub score: i32,
    /// Runner-up extension score, same domain as `score`; 0 = none found.
    /// Preferred MAPQ input, because it compares the two placements the
    /// extender actually scored rather than their anchor coverage.
    pub second_score: i32,
    /// Best chain score (MAPQ input).
    pub chain_score: i32,
    /// Best non-overlapping second chain score, same domain as `chain_score`;
    /// 0 = none found. Used when no runner-up survived to extension.
    pub second_chain_score: i32,
    /// CIGAR over the whole read.
    pub cigar: Vec<CigarOp>,
    /// Number of anchors in the best chain.
    pub n_anchors: usize,
    /// Splice junctions crossed (empty = unspliced).
    pub junctions: Vec<RefinedJunction>,
    /// Editing-aware tolerated substitutions, split by the transcript strand
    /// each implies: ref A/read G and ref T/read C (see
    /// [`crate::extend::exempt_dir`]). Their sum is the EA tag; whether one of
    /// them is zero is what [`mapq_of`] tests.
    pub ea_ag: u32,
    /// Editing-aware ref T/read C count (minus-strand transcript).
    pub ea_tc: u32,
    /// Ordinary mismatch count (EK mm).
    pub mm_count: u32,
    /// Seed count (EK seeds).
    pub n_seeds: usize,
    /// Rescued / re-seeded record (RE tag).
    pub rescued: bool,
    /// MAPQ already decided by a layer that does not produce the evidence
    /// [`mapq_of`] reads, or by a caller-side policy. `None` = compute it.
    ///
    /// Carries an explicit quality estimate from a placement layer without
    /// encoding it into chain scores or creating a synthetic runner-up.
    pub mapq_override: Option<u8>,
    /// Layer that produced the winning placement.
    pub origin: Origin,
    /// Number of distinct-locus candidates the arbiter saw.
    pub n_cands: usize,
}

impl ReadAlignment {
    /// Total editing-aware tolerated substitutions (EA tag, NM contribution).
    #[inline]
    pub fn ea_count(&self) -> u32 {
        self.ea_ag + self.ea_tc
    }

    /// Are the tolerated substitutions all from one transcript strand?
    ///
    /// A-to-I editing on one transcript can only produce one direction, so a
    /// placement carrying both is not an edited read but a paralogous locus
    /// whose base differences happen to land in both exempt classes.
    ///
    /// Note that this stays reachable under a strict (`EditDir::Off`) round:
    /// `align::recount_mm_ea_in` classifies every substitution by
    /// [`crate::extend::exempt_dir`] whatever the extender charged for it, so
    /// the two counters describe the read, not the scoring.
    #[inline]
    pub fn ea_single_direction(&self) -> bool {
        self.ea_ag == 0 || self.ea_tc == 0
    }

    /// Does this placement's mismatch spectrum look hyperedited?
    ///
    /// Three conditions, all of them properties of the read against the
    /// reference rather than of how it was placed: the tolerated substitutions
    /// run in a single direction, they reach [`HYPEREDIT_MIN_DENSITY`] of the
    /// aligned length, and no more than [`HYPEREDIT_MAX_NON_EDIT_MM`] ordinary
    /// mismatches accompany them. It is the criterion the Tier-2 rescue
    /// verifies raw reads against (`esperanto_flow`'s
    /// `RESCUE_MIN_EDIT_DENSITY` / `RESCUE_MAX_NON_EDIT_MM`, which alias the
    /// constants below), and the mismatch-spectrum test of the hyperediting
    /// literature (Porath, Carmi & Levanon, Nat Commun 5:4726, 2014): one
    /// direction must dominate, densely, without a second explanation.
    ///
    /// Density is taken over aligned columns, not read length, so a clipped
    /// placement is judged on the part that was actually compared.
    pub fn hyperedited_signature(&self) -> bool {
        if !self.ea_single_direction() || self.mm_count > HYPEREDIT_MAX_NON_EDIT_MM {
            return false;
        }
        let blen: u32 = self
            .cigar
            .iter()
            .map(|op| match *op {
                CigarOp::Match(n) | CigarOp::Ins(n) | CigarOp::Del(n) => n,
                CigarOp::SoftClip(_) | CigarOp::RefSkip(_) => 0,
            })
            .sum();
        blen > 0 && f64::from(self.ea_count()) / f64::from(blen) >= HYPEREDIT_MIN_DENSITY
    }

    /// Aligned columns (`M`/`I`/`D`) and, of those, matching ones.
    ///
    /// The minor editing direction counts against the match total, the major
    /// one does not. One transcript can only be edited in one direction, so
    /// `min(ea_ag, ea_tc)` is the part of the tolerated substitutions that
    /// editing cannot explain — the same discriminator the hyperediting
    /// literature applies to a read's mismatch spectrum (Porath, Carmi &
    /// Levanon, Nat Commun 5:4726, 2014, require one direction to dominate),
    /// here charged continuously instead of as a threshold. A genuinely
    /// edited read pays nothing; a paralogous copy, whose differences scatter
    /// across both classes, pays for the smaller class.
    fn identity(&self) -> f32 {
        let (blen, m) = self.aligned_and_matching();
        if blen == 0 {
            return 0.0;
        }
        m.saturating_sub(self.unexplained()) as f32 / blen as f32
    }

    /// Aligned columns (`M`/`I`/`D`) and, of those, matching ones (`M`).
    fn aligned_and_matching(&self) -> (u32, u32) {
        let mut blen = 0u32;
        let mut m = 0u32;
        for op in &self.cigar {
            match *op {
                CigarOp::Match(n) => {
                    blen += n;
                    m += n;
                }
                CigarOp::Ins(n) | CigarOp::Del(n) => blen += n,
                CigarOp::SoftClip(_) | CigarOp::RefSkip(_) => {}
            }
        }
        (blen, m)
    }

    /// Differences from the reference that A-to-I editing cannot account for:
    /// ordinary mismatches plus the minor editing direction (see
    /// [`Self::identity`], which charges this continuously).
    ///
    /// The exemption of the major editing direction is conditional on the
    /// read otherwise looking like an edited transcript: with more than
    /// [`HYPEREDIT_MAX_NON_EDIT_MM`] ordinary mismatches the substitution
    /// spectrum is that of a diverged paralog, in which the A>G-type
    /// differences are just the transitions that happen to fall in the exempt
    /// class, and every difference is charged.
    #[inline]
    fn unexplained(&self) -> u32 {
        if self.mm_count > HYPEREDIT_MAX_NON_EDIT_MM {
            self.mm_count + self.ea_ag + self.ea_tc
        } else {
            self.mm_count + self.ea_ag.min(self.ea_tc)
        }
    }

    /// Their share of the aligned columns. Read against the reference, so it
    /// says nothing about how the read was placed — which is what makes it
    /// usable as independent evidence in [`mapq_of`].
    fn unexplained_density(&self) -> f64 {
        let (blen, _) = self.aligned_and_matching();
        if blen == 0 {
            return 0.0;
        }
        f64::from(self.unexplained()) / f64::from(blen)
    }
}

/// The MAPQ required for downstream pileup evidence: 20 corresponds to
/// an error probability of 0.01. It also caps rescued placements.
pub const PILEUP_GATE: u8 = 20;

/// Substitution density at which a placement's mismatch spectrum stops being
/// noise and becomes a signal. On a 150 bp read this is 8 substitutions.
///
/// [`mapq_of`] and [`ReadAlignment::hyperedited_signature`] read it in the two
/// directions the same threshold has: reached by the *editing* classes in one
/// direction, the signal is A-to-I editing and the placement is credible
/// despite its divergence; reached by the differences editing cannot explain
/// ([`ReadAlignment::unexplained_density`]), the signal is that the read
/// belongs somewhere else.
pub const HYPEREDIT_MIN_DENSITY: f64 = 0.05;

/// Maximum ordinary mismatches compatible with the editing signature.
/// Exempted editing substitutions alone do not resolve competing paralogs;
/// ordinary substitutions remain part of placement scoring.
pub const HYPEREDIT_MAX_NON_EDIT_MM: u32 = 2;

/// minimap2's `q_coef`.
const Q_COEF: f32 = 40.0;

/// The extender's match score — the unit `score` is denominated in, needed to
/// turn it back into an alignment length for the logarithm.
const MATCH_SCORE: f32 = 2.0;

/// The best/second pair MAPQ is computed from, and the name the evidence tag
/// reports for it.
///
/// The extension pair is preferred whenever it exists (see module docs on the
/// short-read spliced case); the chain pair answers only when no runner-up
/// survived to extension. `unique` means neither domain found a competitor;
/// `policy` identifies an explicit [`ReadAlignment::mapq_override`].
pub fn margin(aln: &ReadAlignment) -> (i32, i32, &'static str) {
    if aln.mapq_override.is_some() {
        (aln.score, 0, "policy")
    } else if aln.second_score > 0 {
        (aln.score, aln.second_score, "ext_margin")
    } else if aln.second_chain_score > 0 {
        (aln.chain_score, aln.second_chain_score, "chain_margin")
    } else {
        (aln.chain_score, 0, "unique")
    }
}

/// MAPQ for a placed alignment (formula and provenance in the module docs).
///
/// A rescued placement is capped at [`PILEUP_GATE`]: it was won on evidence
/// outside its own seeds — a mate's window, or a reduced alphabet in which
/// uniqueness is weaker than the four-letter reference can attest — so it may
/// enter the pileup but must not outrank a directly seeded placement.
pub fn mapq_of(aln: &ReadAlignment) -> u8 {
    let q = match aln.mapq_override {
        Some(q) => q,
        None => scored_mapq(aln),
    };
    if aln.rescued {
        return q.min(PILEUP_GATE);
    }
    // Excess unexplained mismatch density places the read below the pileup
    // gate, regardless of which layer assigned its MAPQ. Single-direction
    // editing does not add to this density; rescued placements return above.
    if aln.unexplained_density() >= HYPEREDIT_MIN_DENSITY {
        return q.min(PILEUP_GATE - 1);
    }
    q
}

/// The scoring formula's own estimate (formula and provenance in the module
/// docs), before the ceilings [`mapq_of`] applies to every placement whatever
/// produced it.
fn scored_mapq(aln: &ReadAlignment) -> u8 {
    if aln.score <= 0 {
        return 0;
    }
    let (best, second, _) = margin(aln);
    // (1 − x²) on the extension pair, (1 − x) on the chain pair, matching the
    // two branches of `mm_set_mapq2`.
    let x = if second > 0 && best > 0 {
        f64::from(second) / f64::from(best)
    } else {
        0.0
    };
    let squared = aln.second_score > 0;
    let margin_term = (1.0 - if squared { x * x } else { x }).max(0.0) as f32;

    let pen_s1 = if aln.score > 100 {
        1.0
    } else {
        0.01 * aln.score as f32
    };
    let cnt = aln.n_anchors as f32;
    let pen_cm = pen_s1.min(if cnt > 10.0 { 1.0 } else { 0.1 * cnt });
    // `score / match_score` is the alignment length the score was earned over.
    let len_term = (aln.score as f32 / MATCH_SCORE).max(1.0).ln();

    let base = aln.identity() * pen_cm * Q_COEF * len_term;
    let v = base * margin_term;
    let q = if v.is_finite() { v as i32 } else { 0 };
    // minimap2's floor: a strictly better best than second is worth at least 1.
    let mut q = if q <= 0 && second < best { 1 } else { q.clamp(0, 60) };
    // A sufficiently separated runner-up permits the pileup floor when
    // identity, anchor support and length also clear that floor.
    if second > 0
        && best - second >= UNIQUE_MARGIN
        && q < i32::from(PILEUP_GATE)
        && base >= f32::from(PILEUP_GATE)
    {
        q = i32::from(PILEUP_GATE);
    }
    q as u8
}

/// Strict-score margin at which a runner-up stops counting as a competitor
/// for the pileup gate.
pub const UNIQUE_MARGIN: i32 = 6;

/// Pin [`mapq_of`] to 0 — a caller-side policy for last-resort placements
/// (L1-only, repeat-family rescue) that must not be trusted as unique.
///
/// Sets an explicit override while preserving the placement evidence.
pub fn force_mapq0(aln: &mut ReadAlignment) {
    aln.mapq_override = Some(0);
}
