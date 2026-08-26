# ESPERANTO

![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)
![Language: Rust](https://img.shields.io/badge/language-Rust-orange.svg)
![Version](https://img.shields.io/badge/version-1.0.0-green.svg)

ESPERANTO is a command-line toolkit for detecting RNA editing (A-to-I) sites from
RNA-seq data. It runs quality control, splice-aware alignment, candidate calling,
and deep-model scoring as a single deterministic pipeline, and reports editing
sites in VCF format with calibrated probabilities.

The whole pipeline ships as one self-contained binary with no runtime
dependencies beyond a reference FASTA and its index.

## Features

- Adapter trimming and quality filtering for single-end and paired-end FASTQ
- Splice-aware, editing-aware read alignment with 2-pass junction discovery;
  a transcriptome-first (L1) engine is enabled by default when an L1 bundle is
  found next to the alignment index
- Unmapped reads get a second pass against a collapsed-alphabet index (A==G,
  T==C), which recovers hyperedited reads; survivors are written back with
  MAPQ 0 and an `RE:Z:collapsed` tag
- Strand-resolved candidate calling from BAM or the binary `.baln` channel
- Site scoring with a frozen Caduceus-Mamba model (pileup veto gate, 5-fold
  ensemble); AUROC 0.998 on the frozen evaluation corpus
- Identical input produces byte-identical output, independent of thread count
- Reference guardrail: refuses mismatched references (e.g. hg19) before running

## Installation

Pre-built binary (Linux x86_64):

```sh
curl -fsSL https://github.com/Zhazhupai603/esperanto/releases/latest/download/install.sh | sh
```

This installs `esperanto` into `~/.local/bin` and the model bundle into
`~/.local/share/esperanto/bundle`. Re-run the same command to upgrade.

Build from source (stable Rust toolchain):

```sh
git clone https://github.com/Zhazhupai603/esperanto.git
cd esperanto
cargo build --release
./target/release/esperanto --help
```

## Quick start

First-time reference setup — one command:

```sh
esperanto setup
```

Setup works in the refs directory (`~/.local/share/esperanto/refs`; override
with `ESPERANTO_REFS`). It detects a reference FASTA + GTF placed there, or
downloads the GENCODE GRCh38 reference + annotation when the directory is
empty, validates the files before building (both must parse, GTF contig names
must match FASTA contig names, and the species guardrail must hold), and
builds the index in place. It also installs the scoring model bundle when
none is present.

After setup, the analysis needs no reference flags at all:

```sh
esperanto run --r1 reads.fq.gz --out out/
```

To use your own reference instead, place `<name>.fa` and `<name>.gtf` in the
refs directory (compressed `.gz` files are decompressed automatically) and run
`esperanto setup` again.

Alternatively, build the alignment index by hand (once per reference):

```sh
samtools faidx ref.fa
esperanto index --fasta ref.fa --out ref.paidx

# with a transcript annotation, the L1 engine bundle is built in the same pass
esperanto index --fasta ref.fa --gtf genes.gtf --out ref.paidx
# -> ref.paidx + ref.bndl + ref.tidx + ref.cpaidx (collapsed rescue index)
```

Run the full pipeline:

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

The main output is `out/sites.vcf`: one row per called site, `FILTER=PASS`
when `RE_PROB >= 0.5`, with `RE_PROB`, `VAF`, `DEPTH`, strand, and evidence
annotations.

The bundled model targets human hg38. References whose `chr1` length matches
neither hg38 (248956422) nor a small test genome (< 10 Mb) are rejected.

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
| `pile` | 8-dim pileup features for a single site or a site list |
| `scan` | Strand-resolved candidate editing-site discovery (BAM or `.baln`) |
| `score` | RE_PROB scoring: encoder + pileup veto gate + 5-fold ensemble |
| `run` | Full pipeline: qc → map → sort → scan → filter → score → vcf |

Each subcommand supports `--help`.

Common options:

| Option | Applies to | Meaning |
|---|---|---|
| `--threads N` | most | worker threads (0 = all cores) |
| `--lib stranded` | `scan`, `run` | stranded (dUTP-style) libraries; default `unstranded` |
| `--bundle PATH` | `score`, `run` | override the auto-detected model bundle |
| `--l1-bundle PATH` | `map`, `run` | override the auto-detected L1 engine bundle |
| `--batch N` | `score`, `run` | score batch size (default 64) |

When an option is omitted, inputs are resolved automatically (first hit wins):

- Model bundle: `ESPERANTO_BUNDLE`, package layout (`bundle/` next to the
  binary), user data dir (`~/.local/share/esperanto/bundle`).
- L1 bundle: `--l1-bundle`, `ESPERANTO_L1_BUNDLE`, `<index stem>.bndl` (with its
  `<index stem>.tidx` sidecar) next to the alignment index. If none is found,
  alignment runs the genomic layer only (a note is printed).
- Rescue index: `<index stem>.cpaidx` next to the alignment index. When
  present, the map stage re-aligns the unmapped set against it; absent, the
  pass is skipped silently.
- Reference files: `ESPERANTO_REFS`, package `refs/`, user data dir, `./refs`.
  A refs directory provides `hg38.fa`, a GTF, and a gnomAD VCF when present.

## Input requirements

| Entry | Required | Optional |
|---|---|---|
| FASTQ (PE/SE) | `--r1` (plus `--r2` for paired-end); `--index` and `--fasta` are auto-resolved after `setup` | `--gtf`, `--l1-bundle`, `--gnomad` |
| BAM | `--bam` (sorted + indexed), `--fasta` | `--gtf`, `--gnomad` |
| BAM + sites | `--bam`, `--sites`, `--fasta` | `--gtf`, `--gnomad` |

`--sites` is one site per line, `chr<TAB>pos` (1-based), and skips candidate
discovery.

## Output

| File | Content |
|---|---|
| `out/scan/candidates.bed` | candidate sites (coordinate, strand, evidence, score, depth, allele frequency) |
| `out/score/scores.tsv` | `chrom`, `pos`, `RE_PROB` for every scored site |
| `out/sites.vcf` | final calls |

## Guarantees

- Deterministic: identical input produces byte-identical output at any thread
  count.
- Reference guardrail: mismatched references are rejected before compute.
- Every crate has a frozen specification under `docs/specs/`; scientific
  decisions are registered in `docs/DESIGN-DECISIONS.md`.

## Documentation

- `docs/specs/` — per-crate contracts and engineering conventions
- `docs/DESIGN-DECISIONS.md` — scientific decisions and frozen model contracts
- `CHANGELOG.md` — release history

## References

- Schiff et al. Caduceus: bi-directional equivariant long-range DNA sequence
  modeling. 2024.
- Gu and Dao. Mamba: linear-time sequence modeling with selective state
  spaces. 2023.

## License

MIT (see [LICENSE](LICENSE)). Model weights are released under the same terms
(see `model_card.json` in the bundle).
