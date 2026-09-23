# ESPERANTO

[![Version](https://img.shields.io/badge/version-1.0.8-green.svg)](https://github.com/Zhazhupai603/esperanto/releases)
![Language: Rust](https://img.shields.io/badge/language-Rust-orange.svg)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

ESPERANTO detects candidate A-to-I RNA editing sites from human RNA-seq data.
It accepts single-end or paired-end FASTQ, or an existing BAM, and writes VCF results.

## Features

- **Hyperediting rescue:** Recovers heavily edited reads missed during initial alignment and reports the editing sites and regions they support.
- **Alignment verification:** Re-aligns candidate-supporting reads to check for alternative placements before accepting their support.
- **RNA editing scores:** Scores ordinary candidates for RNA editing versus DNA variation, directly from the read alignments.
- **Interactive results:** Browse editing sites and hyperedited regions, search by gene, and inspect predicted amino acid changes in an HTML report.

## Installation

Linux x86-64 packages include the executable and scoring models. No GPU is required.

```sh
curl -fsSL https://github.com/Zhazhupai603/esperanto/releases/latest/download/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"
esperanto --version
```

The installer selects a compatible CPU build and installs under `~/.local/`.
Run the same command to update, or download a package from [Releases](https://github.com/Zhazhupai603/esperanto/releases).

<details>
<summary>Build from source</summary>

Requires Rust, a C/C++ build toolchain, CMake, pkg-config, libclang, and the
OpenSSL, zlib, bzip2 and xz development libraries. Install the release package
above first to obtain the scoring models; compiling the source builds only the executable.

```sh
git clone https://github.com/Zhazhupai603/esperanto.git
cd esperanto
cargo build --release --locked
export PATH="$PWD/target/release:$PATH"
esperanto --version
```

The default build targets x86-64-v3. For older x86-64 CPUs, build with
`RUSTFLAGS='-C target-cpu=x86-64' cargo build --release --locked`.

</details>

## Quick start

### Prepare the reference

After installation, choose an empty directory for the reference and run setup once:

```sh
export ESPERANTO_REFS="$PWD/refs"
esperanto setup
export ESP_INDEX="$ESPERANTO_REFS/GRCh38.primary_assembly.genome.paidx"
```

This downloads the human GRCh38 reference and GTF annotation and builds the
alignment indexes. Keep these environment variables set when running the examples below.
To use your own GRCh38 files, place a `.fa` file (for example, `hg38.fa`) and its
matching GTF in the reference directory before setup. Set `ESP_INDEX` to the
generated main `.paidx` file.

### Run from FASTQ

Replace the example filenames with your reads:

```sh
# Paired-end
esperanto run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
    --index "$ESP_INDEX" --sample sample --out results

# Single-end
esperanto run --r1 single.fq.gz \
    --index "$ESP_INDEX" --sample single --out results
```

Use `--threads 16` to limit CPU use. Multiple lanes can be supplied as
comma-separated filenames to `--r1` and `--r2`, in matching order.

### Run from BAM

Use a coordinate-sorted BAM with a `.bai` or `.csi` index and the same reference
used for alignment:

```sh
esperanto run --bam sample.sorted.bam \
    --index "$ESP_INDEX" --sample aligned --out results
```

Add `--sites sites.tsv` to score an existing site list instead of discovering
candidates. The file has two tab-separated columns, chromosome and 1-based
position, with no header. Use `--fasta` and `--gtf` to override reference discovery.

### Resume or regenerate a report

```sh
esperanto resume results/sample
esperanto report --out results/sample
```

HTML reports require a matching GTF annotation, supplied by the default setup or
with `--gtf`. Without it, the pipeline still writes the VCF but skips the report.

### Human knock-in mouse samples

```sh
esperanto run --r1 knockin_R1.fq.gz --r2 knockin_R2.fq.gz \
    --hybrid APOE --sample knockin --out results
```

The first hybrid run downloads the mouse reference and builds the combined index.
Human gene loci are scored; mouse-contig sites have no model score.

## Output

For `--out results --sample sample`, files are written under `results/sample/`:

| File | Contents |
|---|---|
| `sites.vcf` | Candidate sites, scores, filter status and annotations |
| `sample.report.html` | HTML report; open in a browser |
| `map/sorted.bam` | Sorted alignments from FASTQ input |
| `scan/candidates.bed` | Discovered candidates and base counts |
| `score/scores.tsv` | Chromosome, position and editing score |
| `hyper/regions.bed` | Hyperedited regions from rescued reads |
| `score/hyper_scores.tsv` | Scores for rescued sites |

Files depend on the input and completed stages. Hyperediting output requires
ESPERANTO rescue records; it is not generated from an external BAM alone.
Scores are model predictions. Check the VCF `FILTER` field before using calls;
`UNSCORED` means no model score is available.

## Commands

| Command | Purpose |
|---|---|
| `setup` | Prepare reference files and indexes |
| `run` | Run the analysis from FASTQ or BAM |
| `resume` | Continue an interrupted run |
| `report` | Regenerate an HTML report |
| `index` | Build alignment indexes |
| `qc` | Trim and filter FASTQ reads |
| `map` | Align reads |
| `scan` | Find candidate editing sites |
| `verify` | Check candidate-supporting alignments |
| `score` | Score a supplied site list |
| `pile` | Extract per-site pileup features |
| `update` | Update the installation |
| `delete` | Remove installed reference or model files |

Run `esperanto --help` or `esperanto run --help` for options.
See the [changelog](CHANGELOG.md) for release updates.

## License

[MIT](LICENSE).
