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

## Scoring model v1.7.0 (frozen 2026-09-14)

Scoring is replaced by two donor-disjoint trained networks (identical
architecture, 847k parameters) over 16x50 flank-center joint base-frequency
features extracted straight from the BAM (`cmc-cigar-window-1`; MAPQ >= 20,
base quality >= 13, flags 3844 excluded, per-template best-quality base
consensus, collapsed-rescue reads included with their original query
bases). The model reads no reference sequence. Labeling: DNA-mutation
negatives require the exact reference>alternate allele in the donor's
high-coverage WGS genotype; RNA-editing labels are weak candidate labels
(WGS noncalls do not establish callable reference genotypes). The main
network (DNA mutation vs ordinary editing) reaches AUROC 0.9996 on the
held-out 20-donor test set; the hyper network (ordinary editing vs
rescue-supported editing candidate) scores Tier-2-only sites. Decision
thresholds are validation-donor calibrated (main: DM FPR <= 1%, 0.2854;
hyper: ordinary-editing FPR <= 5%, 0.3624) and ship inside the bundle.
The previous encoder + fusion-head path (and its GPU channel, embedding
cache and mouse bundle) is removed.

## Candidate convention (2026-09-04)

Candidate discovery follows the REDItools counting convention: reads with
MAPQ >= 20, bases with Phred >= 30, no read-orientation split, coverage
>= 10, >= 2 reads carrying the editing-consistent alternate (G at a
reference A, C at a reference T), frequency >= 1%. On the same BAM the
resulting table is identical to REDItools `analyze -q 20 -bq 30` under the
same rule (89,837 sites on the validation sample, every depth and base
count equal). Motivation: the previous caller counted MAPQ 0 reads, split
evidence by read orientation (which discards half of it for unstranded
libraries and zeroes low-frequency sites whose variant reads all fall in
one orientation), ranked with a logistic score whose threshold arm let
non-editing substitutions through, and lost exactly-two-read sites to a
decimal-rounding artifact in the downstream filter. The transcript strand
of an A-to-I site is implied by its reference base (A → plus, T → minus)
and is no longer inferred from junctions or annotation. The scoring
pileup uses the same MAPQ gate, so the model sees the reads the scan
counted.

Known SNPs (dbSNP common) are removed from the candidates before scoring
when a SNP VCF is present; the retrained model rejects most homozygous
sites on its own but heterozygous common SNPs in the editing direction
remained the main source of high-VAF false positives.

The pileup-only veto gate was removed from scoring: it skipped the encoder
for 3.6% of real candidates at the cost of a second head set in the
bundle and a zero-padding shim, and every site now takes the full model
path.

## Strand consistency, neither level (2026-09-05)

A-to-I editing on a single transcript is single-stranded. Seen against the
forward genome, a plus-strand transcript's edits are only ref A / read G and
a minus-strand transcript's only ref T / read C. Two consequences were
implemented and measured; neither survived, and the site-level one taught
the more useful lesson.

**Site level (implemented, then dropped).** A candidate whose
one-read-length neighbourhood held both directions was dropped together
with the opposite-direction neighbours that made it mixed. On the simulated
corpus this cut false candidates 905 -> 172, FDR 3.84% -> 0.76%, for 292
true sites (1.9% of hyper-edited ones), and the same filter on a STAR BAM
moved 0.43% -> 0.27% FDR at no cost to hyper-edited sites — which was read
as evidence that it was not compensating for this aligner.

Real data falsified that reading. It removed 313,136 esperanto candidates
(80.3%) against 13,609 of STAR's (17.5%), and kept 41.78% of esperanto's
REDIportal calls against 83.31% of STAR's. A filter that costs one aligner
half its real calls and the other a sixth is compensating for the first
one. The filter is gone (`esperanto_scan::strand`,
`BlockAcc::max_read_len` and `CallStats::strand_conflict` with it).

The cause was first diagnosed as alignment exempting ref-A/read-G and
ref-T/read-C *simultaneously*, letting a mis-placed read carry tolerated
substitutions in both directions at once, so that much of the mixing the
filter detected was manufactured. **That diagnosis is retracted.** With
strict four-letter Tier-1 alignment, one direction per rescue index, and
`min(ea_ag, ea_tc)` charged against `identity`, the same filter re-measured
post-hoc on the new architecture's real-data candidates removes 80.6% of
them against STAR's unchanged 17.5% (`mapbench/report-esp-final-20260907.md`,
task 3) — the asymmetry is intact to within 0.3 pp of the old aligner's
80.3%. Removing the simultaneous exemption changed nothing, so it was not
the source. The mixing the filter detects is real and its origin is
elsewhere; the mechanism is open.

A second candidate cause was raised by the same run, measured, rejected —
and then the rejection was itself retracted on re-measurement. On chr1,
403,659 gate-passing reads exceed a 1%-of-aligned-length mismatch ceiling;
94.8% sit at exactly MAPQ 60 and 53.5% carry no editing exemption, and the
divergence ceiling in `mapq::mapq_of` exempted them by construction
through a `q < 60` guard. That looked like a hole. The first measurement
said it was not: the `(q = 60, over ceiling, ea = 0)` population is
275,049 reads of which 732 are mis-placed — 0.27%, the same rate as the
whole gate-passing population — so capping it would discard 274,317
correct reads to remove 732 wrong ones.

That measurement was taken at the wrong threshold, and threshold choice
was its whole content. A ceiling of 1% of read length is 1.5
substitutions on a 150 bp read, inside sequencing-error range, so the
population it selected was about a hundredfold larger than the one the
gate acts on and correspondingly diluted. Re-measured at the gate's own
threshold (`mapbench/ceil60.py --ceiling 0.05`): gate-passing reads at
exactly MAPQ 60 whose unexplained density reaches
`HYPEREDIT_MIN_DENSITY` are 676 mis-placed of 1,592 — **42.5%**. MAPQ 60
asserts an error rate of 1e-6; 42.5% is MAPQ 4. The guard is gone, and
the gate now runs in `mapq_of` after both MAPQ paths converge rather than
inside the genomic formula's tail — 54.5% of the real sample's
gate-passing over-dense reads carry `mapq_src=l1` and every one sat at 60,
so a gate on the formula's own output reached under half the population it
is for. Hyperediting is unaffected by construction: `unexplained` charges
`min(ea_ag, ea_tc)`, zero for a single-direction read however densely
edited. Real-sample effect, chr1: candidates 413,854 -> 116,978, and the
esperanto-only partition's REDIportal enrichment goes from 3.99% — below
its own 5.46% displaced control, i.e. no evidence of editing — to 10.21%
against the same control, 1.87x.

The measurement carries a second, larger consequence. The simulated corpus
cannot discriminate here, because it is built from annotated transcripts
and therefore contains no read whose true locus is absent from the
reference — no rRNA, no contamination, no unassembled or unannotated
sequence. Real data does. A read from an absent locus places on its
nearest paralogue with no competitor, which is exactly the input for which
MAPQ 60 is wrong and exactly the input the corpus omits. That matches the
read-level finding below (these reads fail because the true locus is not
in the candidate list) and it explains why a site-level win of 3.8x on
simulated hyper-edited sites coexists with real-data candidates carrying
no enrichment over their own size-matched control. Specificity on real
data is therefore not a scoring problem and is not measurable on this
corpus; both the mechanism and the corpus are open.

**Read level (rejected).** Ranking competing placements one direction at a
time was implemented and reverted. The premise holds — single-direction
placements are 0.27% mis-placed against 15.98% for those with two
minority-direction exemptions and 74.15% for three or more — but the
affected population is 5,903 of 215,418 ambiguous reads and re-ranking left
it unchanged (`mapbench/eadir.log` vs `eadir7.log`). These reads fail
because the true locus is absent from the candidate list, not because it
was out-ranked; no scoring change can fix that. Only the single definition
of the exemption (`extend::EditDir` / `exempt_dir`, replacing four
hardcoded copies of the base test) was kept.

An earlier claim that the read-level split addressed the root cause of the
215k ambiguous population was wrong and is retracted: direction mixing
accounts for 8.4% of that population's errors, and 58% of them come from
reads carrying no exemptions at all.

## MAPQ floor for uncontested placements (2026-09-05)

`mapq = 60 x (1 - second/best) x min(1, n_anchors/10)` holds a placement
below the downstream pileup gate whenever the read yielded few anchors,
even when nothing competes with it in either scoring domain. Such
placements are now floored at MAPQ 20 provided their editing-aware
substitutions are single-directional.

Both halves of the condition are derived rather than fitted. 20 is
`-10*log10(0.01)`, the error rate the pileup gate already assumes, so the
floor claims no more accuracy than the gate requires; the admitted
population measures 0.844% wrong, i.e. honest MAPQ 20.7. Direction-mixing
is the exclusion because a placement carrying both exempt directions is a
paralogous locus rather than an edited transcript, and those reads measure
18.9-98.5% wrong. An `ea <= 5` magnitude cutoff was considered and rejected:
it would admit 922 more reads at the price of a tuned number.

The floor is not a waiver of the anchor divisor. Waiving it would yield
MAPQ 60 — a claimed error rate of 1e-6 for a population measured at
0.844% — which the honest-MAPQ criterion forbids.

Measured: MAPQ >= 20 read coverage 90.56% -> 91.87% at 0.216% -> 0.226%
error; +223 true edited sites (151 hyper-edited) for +4 false ones. With
the site filter in place the hyper-edited yield (5,753) exceeds the
unfiltered baseline (5,709), so the two changes together cost no
hyper-editing sensitivity and cut FDR five-fold.

## Known-SNP exclusion: internal channel, never exposed (2026-09-22)

The `--snp` option is removed from `run` and `scan`, and `resume` ignores a
known-SNP setting frozen by an older run.json (with a warning). The channel
itself (`CallParams.snp`, the VCF set-intersection in `run_call`) stays in the
code but is `None` on every entry point.

Motivation, two converging reasons:

- **Benchmark protocol contamination.** With refs discovery enabled, a bare
  `esperanto run` silently picked up whatever dbSNP VCF happened to sit in the
  refs directory and excluded those positions from the candidates. Two
  evaluations of the same build on the same data could therefore produce
  different candidate sets depending on local refs state — an undocumented
  protocol change no benchmark can tolerate.
- **Claim consistency.** The product's precision story is the in-caller
  RE-vs-variant scoring, not a blacklist. Shipping a silent exclusion path
  alongside that claim invites the reading that the numbers depend on the
  blacklist; closing the channel keeps every reported number attributable to
  the caller itself.

Future work: variant-awareness may return as a designed, three-layer scheme
(variant-aware alignment seeds, per-site variant scoring features, and an
explicit opt-in exclusion list with provenance reporting). Until that design
exists, the channel stays off rather than half-exposed.

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
