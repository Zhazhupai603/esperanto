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
