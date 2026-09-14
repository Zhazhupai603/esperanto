# Changelog

All notable changes to ESPERANTO are documented in this file.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), versioning follows [SemVer](https://semver.org/).

## [1.0.6] - 2026-09-14

### Changed

- Site scoring is replaced by a new frozen model bundle (v1.7.0): a
donor-disjoint trained CNN + Transformer over 16x50 flank-center joint base
frequency features computed directly from the alignment (schema
`cmc-cigar-window-1`; reads with MAPQ >= 20, base quality >= 13, no
secondary/supplementary/duplicate records; per-template best-quality base
consensus). The model reads no reference sequence — the FASTA is only used
for the species guardrail. On the held-out 20-donor test set the classifier
reaches AUROC 0.9996 (input-matched DEMINING controls: 0.917); the frozen
decision threshold (DM false-positive rate <= 1% on validation donors)
ships inside the bundle and the VCF FILTER criteria follow it
(PASS = RE_PROB >= threshold, printed in the header) instead of the fixed
0.5.
- Tier-2-only sites (hyperedited-region calls from rescued reads) are now
model-scored instead of reported unscored: the bundle carries a second
rescue-trained network whose frozen threshold (5% ordinary-editing
false-positive rate on validation donors) decides PASS; sites below it (or
without a score) keep FILTER=HYPER. HYPER_READS/HYPER_REGION annotations
are unchanged.
- The score stage runs on the CPU — the previous Caduceus-Mamba encoder and
  its fusion heads, the embedding cache, the CUDA/GPU channel and the
  `--device`/`--caduceus` options are removed entirely. The new network is
  small enough that no accelerator is needed (16k sites score in ~10 s on 16
  threads).
- Hybrid knock-in runs: the previous mouse bundle was an encoder-based model
  and is gone with the encoder path; mouse-contig sites are now always
  reported UNSCORED with VAF/depth kept. Human-locus sites score with the
  new model as before.
- `score/hyper_scores.tsv` persists the rescue-site probabilities so
  resumed runs reproduce the VCF verbatim.


## [1.0.5] - 2026-09-11

### Changed

- Candidate discovery now follows the standard REDItools counting
  convention (MAPQ >= 20, base quality >= 30, coverage >= 10, >= 2
  editing-consistent reads, frequency >= 1%, no read-orientation split);
  on the same BAM the candidate table is identical to REDItools, site for
  site and count for count. `candidates.bed` carries the reference base,
  depth, A/C/G/T counts and the editing frequency; the transcript strand
  follows from the reference base. The former call score, evidence codes,
  junction/annotation strand inference, C>U mode, `--lib`, the `.baln`
  scan channel and the separate candidate filter stage are gone.
- Known SNPs are removed from the candidates: `run`/`scan` take `--snp`
  (default: a `*snp*.vcf.gz` in the refs directory, e.g. dbSNP common).
  `--gnomad` is replaced by `--snp`.
- `sites.vcf` reports the edited base pair in REF/ALT (A>G or T>C) instead
  of a symbolic `<RE>` allele, and carries per-base counts in `AC`.
- Scoring no longer runs the pileup veto gate; every site takes the full
  encoder + head path. The pileup features are computed over reads with
  MAPQ >= 20, the same read set the candidate scan counted. `score --baln`
  is removed.
- Paired-end mates on opposite strands are proper in either order (FR or
  RF), and the pair span is measured on the transcript (intronic skips
  excluded), so spliced mates adjacent on the transcript are proper pairs.
- Map-stage throughput is ~3.9x on the production sample (151 bp
  paired-end, 38.6M pairs, 32 threads: 34.8 min measured, against a ~133
  min extrapolation of the previous build) with byte-identical output. The
  Gotoh DP kernels — banded and full-matrix — are AVX2 anti-diagonal
  vectorized with a runtime scalar fallback; the scalar kernel was
  rewritten (pre-encoded bases, flat traceback storage, no per-row
  clearing); paired-end chunks deduplicate identical sequences before
  primary alignment (~50% of primary calls eliminated on the production
  library); seeding/chaining/reverse-complement buffers are reused per
  worker; BAM and .baln share one binary record encode; read names stay
  bytes end to end; the unmapped FASTQ writer is buffered; the embedding
  cache writes back once per chunk instead of once per miss. Two-pass
  discovery, band escalation, rescue, and MAPQ semantics are unchanged —
  every gate compared raw.bam, discovered junctions, rescued records and
  unmapped reads byte-for-byte against the previous build. `ESP_PROF=1`
  prints per-run DP call/cell counters, `ESP_SHADOW=1` measures the
  escalation ladder read-only, `ESP_NO_SIMD=1` forces the scalar kernel;
  all three are inert by default.
- Hyperediting is now a two-tier architecture. Tier 1 aligns strictly, with
  no editing exemption, so an ordinarily-edited read is placed on its own
  four-letter evidence. Tier 2 takes the reads Tier 1 does not stand behind
  and realigns them under the four reduced alphabets the method requires
  (A>G and T>C, each in both read orientations), following Porath, Carmi &
  Levanon (Nat Commun 5:4726, 2014). Every collapsed placement is
  re-verified against the raw four-letter reference before it is accepted,
  so a reduced-alphabet alignment is only ever a candidate locus; a read
  that verifies under both directions is ambiguous and dropped. Rescued
  placements are capped at the pileup gate (MAPQ 20) rather than forced to
  0, so a verified rescue can contribute a call.
- Tier 2's input pool is the placements Tier 1 declines to explain, not the
  unmapped file. Porath's pipeline runs on the unmapped output of a
  *refusing* aligner (BWA, ~16% of reads), where the pool is enriched for
  hyperediting because the mismatch ceiling threw those reads out. A local
  Smith-Waterman always finds a positive-scoring window in 3 Gb, so
  esperanto's unmapped set is 1.5% of the library and is anti-selected for
  what the channel looks for — the hyperedited reads were already placed,
  badly, and never arrived. The pool now also carries every placement whose
  ordinary-mismatch density reaches `HYPEREDIT_MIN_DENSITY`, which is the
  same population the density gate had already capped below the pileup
  gate, so those records contribute no call as they stand. Real sample:
  pool 575,479 -> 1,681,786 reads, rescued 220 -> 866, mapping rate
  unchanged at 0.98511 (deferred reads are re-placed, not refused).
- The collapsed realignment is gapless by construction rather than by
  refusing any output that contains a gap. With A and G collapsed to one
  letter the dynamic program has many equal-scoring paths, so a gap that
  buys a match it would have had anyway opens for free; refusing the result
  discarded the candidate locus instead of asking what the best ungapped
  placement there is — 55% of all real-sample rejections. The gap penalty
  is now priced above what a whole read can earn, so the extension
  soft-clips instead, and the re-scored loci are policed by the
  non-editing-mismatch ceiling on their own evidence.
- MAPQ follows minimap2's alignment-stage function (`mm_set_mapq2`) for its
  short-read spliced case, replacing a chaining-stage estimate that had lost
  two of its three factors. Three sites that had been choosing a
  placement's score and anchor fields to shape the old formula's output
  — the intron-chain fast path, a chain-domain near-tie truncation, and the
  all-M CIGAR shortcut that depended on it — now report what was measured;
  every candidate is extended, and a tie that survives extension is
  evidence rather than an artifact of deciding the competition in the one
  domain that cannot see the answer.
- Editing enters MAPQ through the identity term instead of exempting a read
  from scrutiny. `min(ea_ag, ea_tc)` is the part of the tolerated
  substitutions that A-to-I editing cannot explain, since one transcript is
  edited in one direction only, and it is charged continuously; at
  `HYPEREDIT_MIN_DENSITY` it also caps MAPQ below the pileup gate. The cap
  applies to every placement, whatever produced it: it is read against the
  reference and so says nothing about how the read was placed, and it no
  longer requires the formula to have doubted the read first. Both
  exclusions were measured, not assumed. Placements at exactly MAPQ 60
  whose divergence reaches that density are 42.5% misplaced (676 of 1,592),
  against the 1e-6 error rate MAPQ 60 asserts; and 54.5% of the real
  sample's over-dense reads come from the transcriptome layer, which
  returns before the genomic formula runs and was therefore never gated.
  Real sample, chr1: candidates 413,854 -> 116,978, and the
  esperanto-only partition's REDIportal enrichment goes from 3.99% — below
  its own 5.46% size-matched displaced control, i.e. no evidence of
  editing — to 10.21% against the same control.
- The strand-consistency filter on candidates is removed. It was licensed
  on the simulated corpus, where it cut FDR 3.84% -> 0.76%, but on real
  data it removed 80.3% of esperanto's candidates against 17.5% of STAR's
  and kept 41.78% of esperanto's REDIportal calls against 83.31% of STAR's
  — an asymmetry that means it was compensating for this aligner rather
  than measuring editing. `BlockAcc::max_read_len` and
  `CallStats::strand_conflict` go with it. The asymmetry's cause is not
  the simultaneous two-direction exemption: re-measured on this
  architecture the same filter still removes 80.6% against STAR's 17.5%,
  so removing that exemption changed nothing and the mechanism is open.

Measured on `mapbench/sim1` (3,999,972 paired reads) against the
STAR + REDItools pipeline on the same sample: reads mapped 99.96% vs
99.04%, correctly placed 98.54% vs 97.88%, mis-placed among reads passing
the pileup gate 0.257% vs 0.289%, junction-exact 96.05% vs 94.32%. On
hyper-edited reads (6+ edits, 97,733 of them) usable and correctly placed
87.95% vs 60.63%, with 1.69% left unplaced against 29.82%. STAR leads on
one metric, reads passing the gate (97.07% vs 96.35%), because it assigns
255 to every uniquely-mapped read and has no intermediate MAPQ to report;
esperanto's own 3-19 bucket is honestly below the gate.

Measured on a real sample (230628-11-1, 38,649,558 paired reads) against
STAR + REDItools on the same FASTQ, scoring candidates by REDIportal
membership: 117,291 candidates carrying 30,474 REDIportal sites against
STAR's 77,759 carrying 29,469. A raw hit count rewards casting a bigger
net, so the esperanto-only partition is scored against its own
size-matched, position-displaced control: 10.57% REDIportal against
5.50%, 1.92x. Of the 346 candidates the two-tier rescue adds over the same
binary with the deferred pool disabled, 89.6% are in REDIportal (STAR's
entire output is 37.9%), 82.7% fall in clusters of three or more within
100 bp across 43 regions, and 73.1% are absent from STAR's output at any
threshold — the cluster signature the method predicts, in the regime
STAR's mismatch ceiling cannot enter. Strand-symmetric database-free
signal-to-noise on chr1 (A>G/T>C against the A>C/T>G and A>T/T>A channels
over the same reference A/T denominator) is unchanged by the rescue at
2.55x, so the added reads land in clusters rather than in the noise floor.

At the site level on the same corpus, scanned with the same caller over
both BAMs so the aligner is the only variable: 8,279 of 8,665 callable
hyper-edited sites recovered against STAR's 2,198 of 2,492 — 3.8x the
hyper-edited sites found, and a callable pool 3.5x larger, meaning most of
the gap is coverage STAR does not produce at all rather than sites it
fails to call. Recall over all injected edits 92.40% vs 90.59%, VAF error
-0.0280 +/- 0.1481 vs -0.0459 +/- 0.1581, zero SNP leakage on both. This
costs false-discovery rate: 1,132 false candidates of 27,310 (4.14%)
against STAR's 81 of 19,119 (0.42%), and only 25 of the 1,132 lie near a
hyper-edited window, so they come from Tier-1 placement rather than from
the rescue.

That partition-level result does not carry to the candidate set as a
whole, and the strand-symmetric panel says so once it is run against
STAR rather than against a previous esperanto build. Genome-wide on the
same sample (23 chromosomes, `mapbench/wgsum.py` over `eindex.py`):
esperanto's editing channel holds 94,537 sites against a mean control
channel of 38,594, a ratio of 2.45x, while STAR on the same reference and
the same caller settings holds 64,867 against 11,969, a ratio of 5.42x.
Every one of the 23 chromosomes agrees in the same direction. The
decisive form of it is the increment: esperanto adds 29,670 editing-channel
sites over STAR and 26,625 control-channel ones, so the marginal sites it
reports over STAR carry a signal-to-noise of 1.11x — indistinguishable
from what the controls say the pipeline manufactures. Implied false share
of the editing channel is 40.8% against STAR's 18.5%, on 8% *less*
covered A/T sequence. The excess is a concentration, not a larger error
budget: on chr1 esperanto's control channel carries 72,568 alternate
reads against STAR's 59,737, 21% more, but crosses the site threshold at
3,974 columns against 1,236, 3.2x as often. The same substitution
burden stacked onto fewer columns is what clears a two-read floor, and
independent sequencing error does not stack — it is one substitution per
column, which is why STAR's larger share of it never becomes a site.

This reproduces, on an independent and stronger control, the earlier
finding in `mapbench/report-esp-final-20260907.md` that the
esperanto-only candidates carry no editing signal (4.05% REDIportal
against 4.37% and 4.47% for their displaced controls). That finding is
therefore not superseded by the partition numbers above; both hold, of
different populations. The 346 sites the rescue adds are real editing and
the rest of the excess is not, and 346 against 29,670 is the proportion
that matters for a precision claim. Tier 1 placement, not the rescue, is
the open defect: a local Smith-Waterman that always emits an alignment
places the reads STAR refuses, and those placements contribute
mismatches to all three channels equally.

The rescue channel and the quantitative output are unaffected by this.
Editing level at the 30,166 REDIportal sites shared with the pre-rescue
build is identical at 98.0% of them (mean difference +0.0004, mean depth
difference +0.043 reads), so re-placing the deferred pool does not
dilute the reported editing level even though a rescued read counts
toward depth without voting on the variant. Against STAR at the 25,023
shared REDIportal sites, median depth and median editing level are equal
(24 reads, 0.1818) with esperanto 0.0064 higher on the mean — the
direction the simulated corpus, where truth is known, attributes to
STAR under-calling (VAF error -0.0280 +/- 0.1481 against -0.0459 +/-
0.1581). The 1%-of-length mismatch ceiling under which alt-read depth
appeared to halve is an analysis-side filter that is not in the pipeline
and is not proposed for it; the pipeline's ceiling is the
`HYPEREDIT_MIN_DENSITY` gate at 0.05, numerically STAR's own
`--outFilterMismatchNoverLmax 0.05`, differing only in discounting the
substitutions A-to-I can explain.

Scope: one sample; the read-level measurements are chr1. This branch is
not ready for a production or publication precision claim against STAR
on real data.

### Fixed

- The VCF FILTER column no longer mixes PASS with blocking tags
  (RETROCOPY/SIMPLE_REPEAT/REALIGN): a row is PASS/LOW_SCORE/UNSCORED only
  when nothing blocks it, and TAGS is written once.
- `run --stop-after score` stopped before the score stage instead of after
  it, so the stage never executed.
- `verify` without an alignment index fails with a clear error instead of
  passing candidates through.
- The pileup stage applies the same minimum base quality (30) as the
  candidate scan (was 13).

- Sites supported by exactly two editing-consistent reads at many depths
  were dropped by a decimal-rounding artifact in the candidate filter.
- `pile` printed eight of the nine pileup features.

## [1.0.4] - 2026-08-31

### Changed

- Collapsed-alphabet rescue now verifies every placement against the true
  reference before accepting it, rejecting near-random matches in
  two-letter space that inflated the mapping rate and candidate counts on
  contaminated samples. `align_qc.json` now reports the verified, rejected,
  and prefiltered rescue counts separately.
- Candidate filtering now gates the recall arm on the editing-consistent
  direction: `candidates.bed` carries the per-strand editing frequency
  (forward A>G, reverse T>C), and the recall arm requires that signal plus
  at least 2 mutation reads. This drops REF=C/G sites (which cannot be
  A-to-I edited) and single-read noise before scoring, restoring the
  A-to-I candidate distribution the model was trained on.
- Site scoring now uses a ninth pileup feature marking hyperedited
  (collapsed-rescue) coverage, recovering heavily-edited sites that
  were previously missed because their pileup showed only the reference
  base.
- The report's hyperedited-region track is now interactive: clicking a
  region opens the per-base view of that window, where sites near a
  rescued hyperedited read are highlighted. Genes and recoded proteins
  whose span overlaps a hyperedited region are flagged in the gene
  search and the recoded-protein table.

### Fixed

- Bundle auto-discovery still pointed at the v1.4.1 model (8-dim pileup)
  after scoring moved to a 9-dim pileup input, so `run` and `score` failed
  at the score stage with "load bundle". The default bundle path now
  resolves to the v1.6.0 model (9-dim).

## [1.0.3] - 2026-08-29

### Added

- Hybrid references for knock-in mouse models: `run --hybrid
  GENE[,GENE...]` stages the mouse reference on demand and splices the
  selected human gene loci (full locus incl. introns) onto the mouse
  baseline, building the hybrid index on first use and reusing it on later
  runs (bare `--hybrid` opens an interactive gene picker). Scoring routes
  per contig via the species.json manifest: human-locus sites use the
  human model bundle, mouse-contig sites use the mouse model bundle
  (donor-disjoint held-out AUROC 0.985 on 10 ENCODE B6CASTF1 brain
  samples); without a mouse bundle installed, mouse-contig sites are
  written UNSCORED with VAF/DEPTH. The manifest also drives the
  guardrail and the report's contig list.

### Changed
- Large reference downloads run as parallel resumable chunks (6 streams,
  per-chunk retry, mirror fallback); a dropped connection resumes where
  it stopped instead of restarting. The mouse reference is staged on
  demand at the first `--hybrid` run instead of during setup, and setup
  no longer downloads the (optional) mouse annotation.

- The HTML report is named `<sample>.report.html` and gains a red
  region-level hyperedited-density track (rescued-read counts per window)
  under the signal tracks in the genome explorer.

### Changed

- Collapsed-rescue reads (`RE:Z:collapsed`) no longer contribute variant
  evidence to candidate calling — their bases are alphabet-ambiguous and
  produced millions of spurious low-quality candidates on rescue-heavy
  samples. They still count toward depth. Candidate volume drops to the
  real scale and the score stage shrinks from hours to minutes on such
  samples.

- The pileup engine is now event-driven: records stream once per region
  (nodes retire through an end-keyed heap, columns are built only at site
  positions, and a MAXCNT-saturated group falls back to the exact per-site
  sweep). Pileup time on candidate-dense real data drops by ~20x with
  bit-identical features elsewhere.

### Added

- Paired-end runs now write the `.baln` fast channel too (previously a
  zero-byte stub), and the scan stage uses it for PE entries with a BAM
  fallback for stub/absent files; `score --baln` is available as an
  experimental pileup source.

### Changed

- Score runs its pileup pass as a single region sweep per worker instead
  of per-64-site fetches, and the default score batch is 256 (numerics are
  unchanged by construction: identical sites, identical record sets,
  identical model). Pileup features now skip `RE:Z:collapsed` records,
  matching the scan evidence rule; on rescue-heavy samples the depth
  features (and downstream gate/score inputs) change accordingly.
 dev/pile-stream
- The installer detects a working NVIDIA driver and downloads the
  GPU-enabled build automatically, falling back to the CPU build when the
  GPU asset is unavailable. On a GPU-enabled build, `--device auto` asks
  once whether to use a detected GPU; machines without a GPU are never
  asked and stay on the CPU path.

## [1.0.2] - 2026-08-27

### Added

- `--device auto|cpu|gpu` on `score` and `run`: a CUDA channel for the score
  encoder (cargo feature `gpu`). On `auto`, a detected GPU triggers a one-time
  interactive ask; `cpu` forces the CPU path; `gpu` errors clearly when the
  build lacks GPU support or no CUDA device initializes. On a 4-core box the
  GPU channel measured ~1.9× faster than CPU scoring; numerics match the CPU
  path within 3e-4 per site.
- `delete` and `update` subcommands (interactive confirm for both).
- `resume` subcommand: continues an interrupted `run` from the first broken
  stage. `run` now writes into a sample-scoped directory `<out>/<sample>/`
  (`--sample`, default derived from the input file name) and freezes inputs
  and resolved parameters into `run.json`; `resume` re-validates every stage
  artifact (BGZF EOF for BAM, full-stream integrity for gzip, parse checks
  for JSON/BED/TSV/VCF), refuses on changed inputs, and re-executes only
  from the first invalid stage. An intact alignment is never re-run: with a
  valid `raw.bam` + `align_qc.json`, at most the collapsed rescue repeats.

### Changed

- `run --out` now takes an output root: artifacts land in `<out>/<sample>/`
  instead of directly in `<out>/` (breaking layout change).

### Fixed

- PE `unmapped.fq.gz` was malformed (name line glued to the sequence),
  crashing the rescue stage with a seq/qual length mismatch.
- The pipeline no longer runs the coordinate sort and the collapsed rescue
  twice per run (both were executed redundantly after the map stage).
- The Phred+64 rejection message now names the actual threshold (lowest
  quality byte >= 64 in the first 10000 reads).

## [1.0.1] - 2026-08-27

### Added

- **setup** — one-step reference environment: detects a reference FASTA + GTF in the
  refs directory (decompressing `.gz` in place), downloads the GENCODE GRCh38
  primary assembly + v44 annotation when empty, writes the `.fai` when missing,
  validates both files (parse, contig-name agreement, species guardrail), builds the
  index set in place, and installs the scoring model bundle when none is present.
- **report** — every run ends with a self-contained HTML report (out/report.html):
  sample metrics, a drill-down genome explorer (chromosome → 1 Mb window → per-base
  site track with zoom/pan), gene search, and the recoded-protein table. Also
  available as `esperanto report --out <dir>` for existing run directories.
- **collapsed-alphabet rescue** — unmapped reads are re-aligned against a
  collapsed (A==G, T==C) index when `<index stem>.cpaidx` sits next to the paidx,
  recovering hyperedited reads; survivors are written back with MAPQ 0 and an
  `RE:Z:collapsed` tag.

### Changed

- `index --gtf` builds the L1 engine bundle (.bndl + .tidx) and the collapsed rescue
  index (.cpaidx) in the same pass as the paidx.
- The L1 engine is default-on (auto-detected next to the index; genomic-layer
  fallback with a note when absent). Single-pass alignment when L1 is present.
- The refs directory now also supplies the paidx index to `run` — after `setup`,
  `run` needs no reference flags.

### Fixed

- Channel-B tail rescue in repeats: GC/AT donors on the plus strand,
  body-overshoot probe placement, 15–19 bp tail coverage, case-insensitive tail
  verification, symmetric left-tail gate, detected splice signal recorded.
- L1↔G arbitration: variation-bearing L1 placements (indel or non-EA mismatch)
  defer to the genomic layer; L1 placement kept as a MAPQ-0 last resort.
- Guard against mismatched L1 bundle/index pairs (clear error instead of a crash).
- Paired-end entries no longer pass an empty .baln to scan; scan→score candidate
  hard-filter restored.

## [1.0.0] - 2026-08-24

### Added

Initial release: a full Rust rewrite of ESPERANTO as a single self-contained binary with no runtime dependencies.

- **qc** — FASTQ quality control: adapter trimming (including single-end auto-detection), quality filtering, per-base statistics, JSON/HTML reports.
- **tidx / txmap** — transcript index and transcript-to-genome projection for the transcriptome-first engine.
- **pile** — 8-dimensional pileup feature extraction, bit-exact against pysam semantics.
- **engine** — splice-aware alignment core with an optional transcriptome-first (L1) index.
- **bamio** — BAM record I/O, plus a coordinate sort module (chunked stable sort with BAI indexing).
- **map** — editing-aware RNA alignment with 2-pass junction discovery.
- **scan** — strand-resolved candidate editing-site discovery from BAM or `.baln`.
- **score** — RE_PROB scoring from a frozen Caduceus-Mamba model with a pileup veto gate and a 5-fold ensemble.
- **flow** — pipeline orchestration (qc → map → sort → scan → score → vcf) with FASTQ/BAM entry modes and a species guardrail.
- **cli** — the `esperanto` binary with `index`/`qc`/`map`/`pile`/`scan`/`score`/`run` subcommands and zero-config resolution (bundle auto-detection, reference discovery).

### Guarantees

- Deterministic output: identical input produces byte-identical results across runs and thread counts.
- Frozen evaluation corpus anchors hold at AUROC **A 0.9982 / B 0.9961**, re-verified at release freeze.
- Species guardrail rejects mismatched references (e.g. hg19) before computation.
