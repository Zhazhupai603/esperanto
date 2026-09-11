# Changelog

All notable changes to ESPERANTO are documented in this file.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), versioning follows [SemVer](https://semver.org/).

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
