# Design Decisions

This document records scientific and engineering decisions made during the
rewrite — especially where the implementation intentionally deviates from
legacy behavior, or where behavior is frozen to preserve consistency with the
scoring model.

## Principle

Where scientific correctness and model consistency conflict, correctness wins
for components that feed our own aligner. But behaviors the frozen scoring
model was trained on are preserved *exactly*: changing them would shift the
model's input distribution and break the evaluation contract.

## Intentional deviations from legacy behavior (fixed)

These are places where the legacy implementation was unscientific and the
rewrite deliberately fixes it:

- **QC phred clamping** — quality bytes below 33 produced negative phred scores
  and amplified trimming artifacts. Now clamped to 0 (saturating), matching the
  intended arithmetic.
- **tidx reverse-complement casing** — legacy lowercased the complement (a→T,
  asymmetric). Now case-preserving (a→t); the index bytes are unaffected since
  k-mer encoding is case-insensitive.
- **engine base casing** — the legacy engine treated lowercase bases as N in
  k-mer / extension / validation. Now case-insensitive (lowercase recognized),
  consistent with the QC casing decision.
- **engine covered-skip removal** — the legacy default skipped k-mer queries
  whose start fell in a covered region, which missed true sites and produced
  over-confident MAPQ for RNA editing. Now a full-information query; this is
  strictly better on both legacy failure modes (verified byte-identical on a
  chromosome-scale corpus). Cost: roughly 2× L1 queries (the engine is not the
  pipeline bottleneck).
- **txmap junction strand normalization** — legacy stored negative-strand
  introns inverted (end, start), so forward strand queries never matched and
  negative-strand isoform splice evidence was lost. Now normalized to forward
  intervals on load (legacy files are normalized on read).
- **tidx k=32 mask** — legacy `1u64 << 64` overflowed, producing all-zero keys
  for k=32. Now `u64::MAX` with a build-time parameter check.

## Frozen model contract (unchanged by design)

These behaviors are what the frozen scoring model was trained on; they are kept
exact so the model's input distribution does not shift:

- **pile: non-ACGT counted in depth** — non-ACGT bases go to an N bucket but
  still contribute to depth / mean / strand / MAPQ.
- **pile: maxcnt = 8000 depth cap** — reads beyond 8000 at one start position
  are dropped silently; quantified negligible (0.001%) on the evaluation
  corpus, and the model was trained with this cap.
- **pile: f32 means** — per-base quality / MAPQ means use f32 division (not f64
  then cast), required for bit-exact feature parity.
- **map: synthetic-record MAPQ** — `round(60 × min(1, hits/10))`, not a
  junction-evidence-aware formula, because MAPQ feeds the pile features into
  the frozen model.
- **score: pathological CIGAR edge case** — the legacy pile mishandled a
  pathological 0-length `M` op; the new pile follows pysam exactly (the model's
  training semantics), so the legacy behavior is intentionally not reproduced.

## Collapsed-rescue verification (frozen 2026-08-30)

Hyperedited-read rescue follows Porath, Carmi & Levanon, Nat Commun
5:4726 (2014): unmapped reads are realigned in a collapsed alphabet
(A==G, T==C), and **every placement is verified against the true
four-letter reference before it is accepted**. A placement survives only
when it looks like a hyperedited read:

- at most 2 non-editing-class mismatches (the collapsed-space analog of
  the paper's ungapped edit distance 2);
- editing-class mismatches (A-to-G / T-to-C, strand-agnostic) dominate:
  >60% of all mismatches (>80% for reads <=60 bp), each at Phred >=30;
- editing density >=5% of read length;
- cluster geometry: first-to-last editing span >=10% of read length,
  not contained in the outer 20% of the read, single-nucleotide share
  inside the cluster <=60%;
- ungapped CIGAR; paired runs additionally require an already-mapped
  mate within 500 kbp on the same contig in the opposite orientation.

Read-level artifact screens (composition bounds, ambiguous-base
fraction, homopolymer/dinucleotide repeats, trimmed mean Phred >=25) run
before the realignment. Rejected placements stay unmapped; the mapping
rate in `align_qc.json` always reflects the primary alignment.

Known deviations from the paper: multi-mapping tie-breaking by
A-to-G-share margin is not implemented (the best chain is kept at MAPQ 0;
verified reads never vote on variants, only depth and the
hyperedited-region track). Motivation: without verification, the
collapsed realignment accepts near-random placements in two-letter space
and inflates both the mapping rate and downstream candidate counts
(measured on a mismatched-species sample: 70% of reads falsely
"rescued").

## Candidate direction gate (frozen 2026-08-30)

The scan->score candidate filter gates the recall arm on the
editing-consistent direction instead of the any-mismatch frequency.
`candidates.bed` now carries the per-strand editing frequency (forward
A>G, reverse T>C) in place of the per-strand any-mismatch frequency, and
the recall arm keeps a site only when that editing signal is >=5% *and*
at least 2 mutation reads support it. Motivation: the scoring model is a
frozen A-to-I (RE vs germline-mutation) classifier trained on A>G/T>C
candidates; the previous any-mismatch recall arm fed REF=C/G sites (which
cannot be A-to-I edited) and single-read noise into it, and the model
scored those out-of-distribution sites as if they were edits (measured on
a human sample: 51% of passing candidates were REF=C/G). The direction
gate restores the A-to-I candidate distribution the model was trained
on. After the VAF-balanced retrain the model itself rejects REF=C/G
sites (2.9% false-positive rate, matching the 3.4% on REF=A/T), so the
gate is now a performance optimization — it halves the candidate set —
rather than a correctness requirement.

## Hyperedited signal (frozen 2026-08-31)

The scoring model takes a ninth pileup feature marking whether a site is
covered by a collapsed-rescue (hyperedited) read. Motivation: hyperedited
reads are excluded from the pileup vote (their bases are
alphabet-ambiguous), so their editing sites show a pure-reference pileup
and were rejected by both scan (var_reads == 0) and score. The ninth
feature lets the model recognize these sites directly; retraining on
19,463 hyperedited sites as RE positives raised their recall from 2.2% to
100% without degrading AUROC (0.9959 protocol A). The veto gate is not
retrained — its threshold is low enough that hyperedited sites pass
through, so its 8-dim weights are zero-padded to 9 dims.

## Pending decisions

Deferred until wet-lab / further data is available:

- **scan amb-strand sites** use `edit_frac = max(fwd, rev)`, slightly inflating
  scores for overlapping dual-strand sites (ordering only, not gating).
- **refine r3 variant reseeding** narrows the fast legacy path; awaiting an
  editing-recall review.
- **align rescued reads** without strand orientation are counted as forward;
  documented, to be revisited.

## Historical lessons

- Direction contracts (reverse-complement on read side vs. stored sequence on
  write side) must be verified on *both* sides.
- When two authoritative references diverge (upstream htslib vs. pysam
  del-chase), the contract follows the one that produced the training data
  (pysam), and the divergence point is documented.
- A performance fast-path that changes a decision is *semantics*, not
  optimization, and must be documented.

## Integration notes

- **tx_id dual-track** — the legacy engine used a full-transcript projection
  for read placement (id space = the transcript index ENST order), separate
  from the filtered transcript map used for attribution. The rewrite reproduces
  both tracks; mixing the filtered map with the transcript index ids would
  mis-map roughly 22% of transcripts.
- **transcript naming** — the transcript index exposes a display name while the
  transcript map keys on the ENST id; cross-crate joins must use the id.
