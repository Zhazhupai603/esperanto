# Changelog

All notable changes to ESPERANTO are documented in this file.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), versioning follows [SemVer](https://semver.org/).

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
