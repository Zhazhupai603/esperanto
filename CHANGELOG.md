# Changelog

All notable changes to ESPERANTO are documented in this file.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), versioning follows [SemVer](https://semver.org/).

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
