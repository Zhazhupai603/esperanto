# ESPERANTO

![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)
![Language: Rust](https://img.shields.io/badge/language-Rust-orange.svg)
![Version](https://img.shields.io/badge/version-1.0.5-green.svg)

ESPERANTO is a command-line toolkit for detecting RNA editing (A-to-I) sites from RNA-seq data. It runs quality control, splice-aware alignment, candidate calling, and deep-model scoring as a single deterministic pipeline, and reports editing sites in VCF format with calibrated probabilities. The whole pipeline ships as one self-contained binary with no runtime dependencies beyond a reference FASTA and its index.

**Pipeline:** `qc` → `map` → `sort` → `scan` → `verify` → `hyper` → `score` → `vcf` (+ HTML report)

## Features

- **End-to-end processing:** Adapter trimming and quality filtering for single-end and paired-end FASTQ.
- **Splice-aware alignment:** 2-pass junction discovery; a transcriptome-first (L1) engine is enabled by default when an L1 bundle is found next to the alignment index.
- **Two-tier hyperediting rescue:**
  - *Tier 1* aligns strictly, with no editing exemption, so an ordinarily-edited read is placed on its own four-letter evidence.
  - *Tier 2* takes the reads Tier 1 refused (unmapped, or placed with most of the read soft-clipped) and realigns them against collapsed-alphabet indexes (A==G, T==C) in the four passes the method requires. Each placement is re-verified against the true four-letter reference, and a read that verifies under both directions is ambiguous and dropped. Verified reads are written back with an `RE:Z:collapsed` tag, capped at MAPQ 20.
- **Candidate-vote re-alignment (`verify` stage):** Every candidate's supporting reads are re-aligned across the whole reference in exhaustive mode; a read's vote is kept only when the site-supporting placement is unambiguously best.
- **Hyperedited-region calling (`hyper` stage):** Collapsed-rescue reads are clustered into regions (dir-stranded), reported as a region track in the HTML report, and sites supported only by Tier-2 (collapsed) evidence are marked `FILTER=HYPER` in the VCF.
- **Standardized candidate calling:** Uses the standard REDItools counting convention (MAPQ >= 20, base quality >= 30, coverage >= 10, >= 2 editing-consistent reads, frequency >= 1%). On the same BAM, the candidate table is identical to REDItools; known SNPs (dbSNP common) are removed when a SNP VCF is present.
- **Deep-model scoring:** Site scoring with a frozen Caduceus-Mamba model (5-fold ensemble) over the same read set the candidate scan counted.
- **Interactive HTML report:** Every run ends with a self-contained HTML report (`<out>/<sample>/<sample>.report.html`, opens in any browser) featuring sample metrics, a drill-down genome explorer, gene search, and a recoded-protein table.
- **Performance:** The `map` stage processes a full 38.6M-pair 151 bp PE RNA-seq sample in about 35 minutes on 32 threads (deterministic output).
- **Deterministic:** Identical input produces byte-identical output at any thread count.
- **Reference guardrail:** Refuses mismatched references (e.g., hg19) before running.

## Installation

### Pre-built binary (Linux x86_64)

```sh
curl -fsSL https://github.com/Zhazhupai603/esperanto/releases/latest/download/install.sh | sh
```

The installer places:

- the `esperanto` binary — `~/.local/bin/`
- the model bundle — `~/.local/share/esperanto/bundle/`

Two CPU builds are published for every release: an AVX2-optimized one
(x86-64-v3, essentially all x86 servers since 2013) and a portable baseline
one. The installer checks the local CPU and picks the faster build when AVX2
is available, so the download never lands on a machine it cannot run on.

Re-run the same command to upgrade.

### Pre-built binary (Linux x86_64)

```sh
curl -fsSL https://github.com/Zhazhupai603/esperanto/releases/latest/download/install.sh | sh
```

The installer places:
- the `esperanto` binary — `~/.local/bin/`
- the model bundle — `~/.local/share/esperanto/bundle/`

Re-run the same command to upgrade.

### GPU acceleration

The installer detects a working NVIDIA driver and picks the GPU-enabled build automatically (CPU-only machines get the CPU build). At runtime, the score stage detects a CUDA GPU and asks once whether to use it; without a GPU nothing is asked and the CPU path runs (measured ~1.9× over a 4-core CPU at 20k sites; numerics match the CPU path within 3e-4 per site):

- `--device auto` (default): ask once when a GPU is detected, CPU otherwise
- `--device cpu`: force the CPU path (never asks)
- `--device gpu`: require a CUDA GPU (clear error when none initializes)

### Build from source

Release builds target x86-64-v3 (AVX2, essentially all x86 servers since 2013) by default; a runtime fallback keeps the packed DP kernels scalar when AVX2 is unavailable (`ESP_NO_SIMD=1` forces scalar), but the default build profile assumes an AVX2-capable CPU. Include the CUDA channel with `cargo build --release --features gpu` (needs the CUDA toolkit at build time and an NVIDIA driver at run time).

```sh
git clone https://github.com/Zhazhupai603/esperanto.git
cd esperanto
cargo build --release
./target/release/esperanto --help
```

## Quick start

### First-time reference setup

```sh
esperanto setup
```

Setup works in the refs directory — default `~/.local/share/esperanto/refs`, override with `ESPERANTO_REFS`. It detects a reference FASTA + GTF placed there, or downloads the GENCODE GRCh38 primary-assembly reference and the GENCODE v48 primary annotation when the directory is empty, fetches the standard decoy set (45S/5S rDNA + hs38d1) by default — repeat-derived reads park on decoy contigs instead of force-mapping onto paralogs (skip with `--no-decoys`) — validates the files before building, and builds the index in place. It also installs the scoring model bundle when none is present.

After setup, no reference flags are needed — the index, FASTA, GTF, and L1 engine bundle are all auto-resolved from the refs directory:

```sh
esperanto run --r1 reads.fq.gz --out out/
```

To use your own reference instead, place `<name>.fa` and `<name>.gtf` in the refs directory (compressed `.gz` files are decompressed automatically) and run `esperanto setup` again.

Alternatively, build the alignment index by hand (once per reference):

```sh
samtools faidx ref.fa
esperanto index --fasta ref.fa --out ref.paidx

# with a transcript annotation, the L1 engine bundle is built in the same pass
esperanto index --fasta ref.fa --gtf genes.gtf --out ref.paidx
# -> ref.paidx + ref.bndl + ref.tidx + ref.cpaidx (collapsed rescue index)
```

### Run the full pipeline

```sh
# paired-end
esperanto run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
    --index ref.paidx --fasta ref.fa --out out/

# single-end: drop --r2
# multiple lanes: --r1 lane1.fq.gz,lane2.fq.gz
```

Or start from a coordinate-sorted, indexed BAM:

```sh
esperanto run --bam sample.sorted.bam --fasta ref.fa --out out/
```

Each run lands in a sample-scoped directory `<out>/<sample>/` (override the sample name with `--sample`). The main output is `<out>/<sample>/sites.vcf`. Open `<out>/<sample>/<sample>.report.html` in a browser for the interactive report.

### Knock-in mouse models (e.g. APOE knock-in)

Setup stages the human reference (the mouse reference is fetched automatically on the first hybrid run, with parallel resumable downloads), and `run --hybrid GENE[,GENE...]` splices the selected human gene loci onto the mouse baseline — building the hybrid index on first use and reusing it on later runs. Scoring routes per contig: human-locus sites use the human model, mouse-contig sites use the mouse model (installed bundles only; without a mouse bundle, mouse-contig sites are reported UNSCORED with VAF/depth):

```sh
esperanto run --r1 reads_R1.fq.gz --r2 reads_R2.fq.gz --out runs/ --hybrid APOE
```

Symbols follow the annotation's gene names: the human APOE locus is `APOE` (the ε2/ε3/ε4 alleles are haplotypes of the same gene, not separate symbols).

An interrupted run resumes with one command — no flags, it re-reads the frozen parameters and continues from the first broken stage (an intact alignment is never re-run):

```sh
esperanto resume <out>/<sample>/
```

The bundled model targets human hg38. References whose `chr1` length matches neither hg38 (248956422) nor a small test genome (< 10 Mb) are rejected.

## Usage

```
esperanto <COMMAND> [OPTIONS]
```

| Command | Purpose |
|---|---|
| `index` | Build a paidx alignment index from a reference FASTA |
| `setup` | One-step reference environment: fetch/detect references + build the index |
| `qc` | FASTQ quality control (trimming, filtering, `qc.json`/`qc.html`) |
| `map` | Splice-aware alignment; single-pass with the L1 engine, 2-pass junction discovery otherwise |
| `pile` | 9-dim pileup features for a single site or a site list |
| `scan` | A-to-I candidate site discovery (REDItools counting convention, SNP exclusion) |
| `verify` | Re-align candidate-supporting reads; withdraw ambiguous votes |
| `score` | RE_PROB scoring: encoder + 5-fold ensemble |
| `run` | Full pipeline: qc → map → sort → scan → verify → hyper → score → vcf + HTML report |
| `report` | Regenerate the HTML report for an existing run directory |
| `resume` | Continue an interrupted run from the first broken stage |
| `delete` | Remove local reference/index data or the model bundle (interactive confirm) |
| `update` | Check for a newer release and install it (interactive confirm) |

Each subcommand supports `--help`.

### Common options

| Option | Applies to | Meaning |
|---|---|---|
| `--threads N` | most | worker threads (0 = all cores) |
| `--snp PATH` | `scan`, `run` | known-SNP VCF whose positions are removed from the candidates (default: `*snp*.vcf.gz` in the refs directory; `scan --snp none` disables) |
| `--min-mapq`, `--min-bq`, `--min-depth`, `--min-edit-reads`, `--min-edit-freq` | `scan` | candidate thresholds (defaults 20 / 30 / 10 / 2 / 0.01) |
| `--bundle PATH` | `score`, `run` | override the auto-detected model bundle |
| `--l1-bundle PATH` | `map`, `run` | override the auto-detected L1 engine bundle |
| `--batch N` | `score`, `run` | score batch size (default 256) |
| `--device auto\|cpu\|gpu` | `score`, `run` | score device; `auto` uses a detected CUDA GPU, CPU otherwise |
| `--sample NAME` | `run` | sample directory name under `--out` (default: derived from the R1/BAM file name) |
| `--hybrid [GENES]` | `run` | knock-in mouse: splice these human gene loci onto the mouse reference (bare flag opens the interactive gene picker) |
| `--stop-after <STAGE>` | `run` | stop after the named stage (`qc\|map\|sort\|scan\|verify\|hyper\|score\|vcf`) (inspection switch; not recorded for resume) |
| `--rescue-control` | `index`, `map`, `run` | additionally build/run the A>C / T>G negative-control collapsed alphabets; control verifications are counted only |
| `--decoy <FASTA>` | `index`, `setup` | extra decoy contigs (renamed with a `decoy_` prefix) included in the alignment index but excluded from the rescue indexes and the L1 bundle (repeatable) |

When an option is omitted, inputs are resolved automatically (first hit wins):

- Model bundle: `ESPERANTO_BUNDLE`, package layout (`bundle/` next to the binary), user data dir (`~/.local/share/esperanto/bundle`).
- L1 bundle: `--l1-bundle`, `ESPERANTO_L1_BUNDLE`, `<index stem>.bndl` (with its `<index stem>.tidx` sidecar) next to the alignment index. If none is found, alignment runs the genomic layer only (a note is printed).
- Rescue index: `<index stem>.cpaidx` next to the alignment index. When present, the map stage re-aligns the unmapped set against it; absent, the pass is skipped silently.
- Reference files: `ESPERANTO_REFS`, package `refs/`, user data dir, `./refs`. A refs directory provides the FASTA, a GTF, the alignment index, and a known-SNP VCF (`*snp*.vcf.gz`, e.g. dbSNP common) when present.

## Input requirements

| Entry | Required | Optional |
|---|---|---|
| FASTQ (PE/SE) | `--r1` (plus `--r2` for paired-end); `--index` and `--fasta` are auto-resolved after `setup` | `--gtf`, `--l1-bundle`, `--snp` |
| BAM | `--bam` (sorted + indexed), `--fasta` | `--gtf`, `--snp` |
| BAM + sites | `--bam`, `--sites`, `--fasta` | `--gtf` |

`--sites` is one site per line, `chr<TAB>pos` (1-based), and skips candidate discovery.

## Output

| File | Content |
|---|---|
| `<out>/<sample>/scan/candidates.bed` | candidate sites: `chrom pos0 pos1 ref depth A C G T edit_freq` |
| `<out>/<sample>/score/scores.tsv` | `chrom`, `pos`, `RE_PROB` for every scored site |
| `<out>/<sample>/sites.vcf` | final calls. `FILTER` is `PASS` / `LOW_SCORE` / `UNSCORED` only when no analysis tag blocks the site; otherwise `FILTER` lists the blocking tags (`SPLICE_ADJ`, `HOMOPOLYMER`, `SIMPLE_REPEAT`, `RETROCOPY`, `REALIGN`). `HYPER` marks Tier-2-only sites. `INFO` carries `RE_PROB`, `VAF`, `DEPTH`, `AC` per-base counts, strand annotations, `TAGS` (filter annotations that never remove rows), and `HYPER_READS`/`HYPER_REGION` where applicable. |
| `<out>/<sample>/<sample>.report.html` | self-contained interactive report (metrics, genome explorer with a red hyperedited-region track, gene search, recoded proteins) |

## Guarantees

- **Deterministic:** identical input produces byte-identical output at any thread count.
- **Reference guardrail:** mismatched references are rejected before compute.
- **Frozen specifications:** every crate has a frozen specification under `docs/specs/`; scientific decisions are registered in `docs/DESIGN-DECISIONS.md`.

## Documentation

- `docs/specs/` — per-crate contracts and engineering conventions
- `docs/DESIGN-DECISIONS.md` — scientific decisions and frozen model contracts
- `CHANGELOG.md` — release history

## References

- Schiff et al. Caduceus: bi-directional equivariant long-range DNA sequence modeling. 2024.
- Gu and Dao. Mamba: linear-time sequence modeling with selective state spaces. 2023.

## License

MIT (see [LICENSE](LICENSE)). Model weights are released under the same terms (see `model_card.json` in the bundle).
