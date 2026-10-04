# cli — command-line interface

The `esperanto` executable parses arguments, resolves reference and model paths, and delegates to library crates. `--threads 0` selects the available worker count.

## Commands

| Command | Required arguments | Optional arguments and behavior |
|---|---|---|
| `setup` | None | `--force`, `--no-decoys`; prepare reference files and indexes |
| `index` | `--fasta`, `--out` | `--gtf`, `--k` (15), `--w` (5), repeatable `--decoy`; build alignment indexes |
| `run` | `--out` and FASTQ or BAM input | `--r1`, `--r2`, `--bam`, `--sites`, `--index`, `--fasta`, `--gtf`, `--bundle`, `--l1-bundle`, `--sample`, `--threads`, `--stop-after` |
| `resume` | Run directory | Validate and continue a compatible interrupted run |
| `map` | `--r1`, `--index`, `--out` | `--r2`, `--gtf`, `--l1-bundle`, `--threads`; write unsorted BAM |
| `scan` | `--bam`, `--out` | `--fasta`, `--threads`; unified candidate discovery |
| `score` | `--bam`, `--sites`, `--out` | `--fasta`, `--bundle`, `--threads`; unified classifier scores |
| `report` | `--out` | `--fasta`, `--gtf`; regenerate HTML from a finished run |
| `update` | None | `--check` reports availability; `--yes` confirms installation |
| `delete` | `refs`, `bundle`, or `all` | `--yes` confirms removal |

`index --gtf` builds the L1 `.bndl` and `.tidx` files alongside the `.paidx`. The A-to-G and T-to-C transformed rescue indexes use `.ag.paidx` and `.tc.paidx`. `map` writes `raw.bam`, `unmapped.fq.gz`, `align_qc.json`, and discovered junctions. `run` performs the subsequent coordinate sort and indexing.

## QC command

`qc --r1 <FASTQ> --out <directory>` accepts optional `--r2`; multiple lanes are supplied as repeated or comma-separated input arguments. Reads are written as gzip-compressed FASTQ.

| Option | Default |
|---|---|
| `--adapter-r1`, `--adapter-r2` | Built-in adapter handling |
| `--disable-adapter-trim`, `--disable-pe-overlap` | false |
| `--qtrim` | false |
| `--qtrim-cutoff` | 20 |
| `--trim-front1`, `--trim-front2`, `--trim-tail1`, `--trim-tail2` | 0 |
| `--polyg` | auto |
| `--min-len` | 15 |
| `--n-max` | 5 |
| `--q15-frac-max` | 0.4 |
| `--keep-unpaired`, `--detect-adapter-se` | false |
| `--threads` | 0 |

## Discovery and scoring

`scan` uses the same observation filters as the scorer: MAPQ at least 20, BQ at least 13, and exclusion of unmapped, secondary, and supplementary alignments. Reference A/T positions with one or more editing-consistent observations are emitted. Thresholds are part of the model contract.

`score` reads a headerless `chrom<TAB>pos` file with 1-based positive positions and writes `chrom<TAB>pos<TAB>score` in input order. The model bundle must contain `model.txt`, `freeze.json`, and `alu.bed`.

## Path resolution

Explicit paths take precedence. Discovery only fills missing values.

Model resolution checks `ESPERANTO_BUNDLE`, the source-tree bundle, bundle directories beside the executable or its parent, and the user-data bundle. The bundled model path is `bundle/human/esperanto-model-v1.1.0/rust`. An unresolved bundle is an error listing attempted paths.

Reference resolution checks `ESPERANTO_REFS`, `refs` beside the executable or its parent, the user-data `refs` directory, and `./refs`. FASTA selection prefers `hg38.fa`; other eligible files are selected in sorted filename order. FASTA files require matching `.fai` indexes.

User data resides under `ESPERANTO_HOME` when set, otherwise under `$XDG_DATA_HOME/esperanto`, or `$HOME/.local/share/esperanto`. `setup` honors `ESPERANTO_REFS` for reference installation.

An explicit L1 bundle or `ESPERANTO_L1_BUNDLE` takes precedence over the bundle next to the alignment index. The bundle needs its matching `.tidx` sidecar. Missing automatic discovery leaves genomic alignment available.

## Exit behavior

Exit status is 0 on success, 1 on a runtime error, and 2 on an argument parsing error. Error messages identify the relevant path or value. `update` and `delete` require confirmation unless `--yes` is supplied; noninteractive invocation without confirmation does not mutate the installation.
