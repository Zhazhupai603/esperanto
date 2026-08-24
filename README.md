# ESPERANTO

<p align="center">
  <b>RNA editing analysis, end to end.</b><br/>
  <em>RNA-seq in → editing-vs-germline answers out.</em>
</p>

<p align="center">
  <img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-blue.svg">
  <img alt="Language: Rust" src="https://img.shields.io/badge/language-Rust-orange.svg">
  <img alt="Version" src="https://img.shields.io/badge/version-1.0.0-green.svg">
</p>

ESPERANTO is a purpose-built toolkit for **RNA editing (A-to-I) analysis**. It takes
RNA-seq reads (or an existing BAM) through a fully deterministic pipeline — quality
control, splice-aware alignment, candidate discovery, and scoring by a frozen deep
model — and reports editing sites with calibrated probabilities. It is a **single
self-contained binary** with **no runtime dependencies**.

---

## What it does

- **QC** — adapter trimming (incl. single-end auto-detection), quality filtering, per-base statistics.
- **Alignment** — splice-aware, editing-aware RNA alignment with 2-pass junction discovery and an optional transcriptome-first (L1) engine.
- **Candidate discovery** — strand-resolved editing-site calling with dual BAM / `.baln` input.
- **Scoring** — RE_PROB (RNA-editing probability) from a frozen Caduceus-Mamba model with a pileup veto gate and 5-fold ensemble.
- **Reporting** — a minimal VCF of editing calls with `RE_PROB`, `VAF`, `DEPTH`, strand and evidence annotations.

## Highlights

- **Deterministic**: identical input → byte-identical output, across runs *and* thread counts.
- **Golden-validated**: model anchors AUROC **A 0.9982 / B 0.9961** on the frozen evaluation corpus, re-verified at release freeze.
- **Species guardrail**: refuses to run against a mismatched reference (e.g. hg19) *before* burning compute.
- **Zero-config**: the bundled model is auto-detected from the package layout — no environment variables required.
- **Contract discipline**: every crate has a frozen specification; scientific deviations are registered and reasoned about.

---

## Installation

### One-line install (recommended)

```sh
curl -fsSL https://github.com/Zhazhupai603/esperanto/releases/latest/download/install.sh | sh
```

Installs the `esperanto` binary to `~/.local/bin` and the model bundle to
`~/.local/share/esperanto/bundle` (no sudo). The binary auto-detects the bundle.
Re-run to upgrade in place.

### Pre-built binary (Linux x86_64)

Alternatively, download the release tarball and unpack it:

```sh
tar -xzf esperanto-1.0.0-linux-x86_64.tar.gz
cd esperanto-1.0.0-linux-x86_64
./bin/esperanto --help
```

The model bundle ships inside the tarball and is located automatically.

### Build from source

Requires a Rust toolchain (stable):

```sh
git clone https://github.com/Zhazhupai603/esperanto.git
cd esperanto
cargo build --release
./target/release/esperanto --help
```

---

## Quick start

### 1. Prepare a reference index (once)

```sh
# A reference FASTA + .fai are required. Create the .fai with:
samtools faidx ref.fa

# Build the alignment index:
esperanto index --fasta ref.fa --out ref.paidx
```

> The bundled model targets **human hg38**. The pipeline refuses references whose
> `chr1` length is neither hg38 (`248956422`) nor a small synthetic/test genome
> (`< 10 Mb`).

### 2a. Run the full pipeline from FASTQ

```sh
# paired-end
esperanto run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
    --index ref.paidx --fasta ref.fa --out out/

# single-end (drop --r2); multiple lanes: --r1 lane1.fq.gz,lane2.fq.gz
```

### 2b. Run from an existing BAM

The BAM must be **coordinate-sorted and indexed** (`.bai`/`.csi`):

```sh
samtools sort -o sample.sorted.bam sample.bam
samtools index sample.sorted.bam

esperanto run --bam sample.sorted.bam --fasta ref.fa --out out/
```

### 2c. Score your own site list (skip discovery)

```sh
esperanto run --bam sample.sorted.bam --sites my_sites.tsv \
    --fasta ref.fa --out out/
```

`my_sites.tsv` is one site per line: `chr<TAB>pos` (1-based).

### 3. Read the results

| File | Meaning |
|---|---|
| `out/scan/candidates.bed` | candidate sites (coordinate, strand, evidence, soft score, depth, allele frequency) |
| `out/score/scores.tsv` | `chrom  pos  RE_PROB` for every site |
| **`out/sites.vcf`** | final calls — `FILTER=PASS` when `RE_PROB ≥ 0.5` (editing), `LOW_SCORE` otherwise |

---

## Usage

```
esperanto <COMMAND> [OPTIONS]
```

| Command | Purpose |
|---|---|
| `index` | Build a paidx alignment index from a reference FASTA |
| `qc` | FASTQ quality control (trimming, filtering, `qc.json`/`qc.html`) |
| `map` | Splice-aware alignment (RNA 2-pass); writes `raw.bam` + `align.baln` |
| `pile` | 8-dim pileup features for a single site or a site list |
| `scan` | Strand-resolved candidate editing-site discovery (BAM or `.baln`) |
| `score` | RE_PROB scoring: encoder + pileup veto gate + 5-fold ensemble |
| `run` | Full pipeline: qc → map → sort → scan → score → vcf |

Every subcommand has `--help` with a full flag reference, e.g. `esperanto run --help`.

### Common options

| Option | Applies to | Meaning |
|---|---|---|
| `--threads N` | most | worker threads (0 = all cores) |
| `--lib stranded` | `scan`, `run` | stranded (dUTP-style) libraries; default `unstranded` |
| `--bundle PATH` | `score`, `run` | override the auto-detected model bundle |
| `--batch N` | `score`, `run` | score batch size (default 64) |
| `--l1-bundle PATH` | `map`, `run` | optional transcriptome-first (L1) engine |

### Zero-config resolution

When flags are omitted, ESPERANTO resolves inputs in this order (first hit wins):

- **Model bundle** — `ESPERANTO_BUNDLE` → package layout (`bundle/...` next to the binary) → user data dir (`~/.local/share/esperanto/bundle`).
- **Reference files** — `ESPERANTO_REFS` → package `refs/` → user data dir (`~/.local/share/esperanto/refs`) → `./refs`. A refs directory provides `hg38.fa`, a GTF, and a gnomAD VCF when present.

---

## Input requirements

| Entry | Required | Optional |
|---|---|---|
| FASTQ (PE/SE) | `--r1` [`--r2`], `--index` (paidx), `--fasta` (+`.fai`) | `--gtf`, `--l1-bundle`, `--gnomad` |
| BAM | `--bam` (sorted + indexed), `--fasta` | `--gtf`, `--gnomad` |
| BAM + sites | `--bam`, `--sites`, `--fasta` | `--gtf`, `--gnomad` |

---

## Documentation

- `docs/specs/` — engineering conventions and per-crate contracts (single source of truth).
- `docs/specs/crates/cli.md` — full CLI flag reference.
- `docs/DESIGN-DECISIONS.md` — scientific decisions, deviations from legacy, and frozen model contracts.
- `docs/BACKLOG.md` — planned, non-blocking enhancements.
- `CHANGELOG.md` — release history.

## References

The scoring model uses a bidirectional state-space (Mamba) architecture with
reverse-complement equivariance, following:

- Schiff et al., *Caduceus: Bi-Directional Equivariant Long-Range DNA Sequence Modeling*, 2024.
- Gu & Dao, *Mamba: Linear-Time Sequence Modeling with Selective State Spaces*, 2023.

The encoder is pretrained from scratch on a transcript corpus; model weights, training
labels and all pipeline code are developed for ESPERANTO.

## License

MIT — see [LICENSE](LICENSE). Model weights are released under MIT (see the bundled
`model_card.json`).
