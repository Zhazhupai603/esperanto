# cli — command-line shell (bin `esperanto`, zero-config resolution)

Thin orchestration layer: parse arguments → zero-config fill-in → dispatch to each crate entry. **Zero scientific semantics**; all numeric/byte contracts belong to the stage crates and flow. 1.0.0 ships 7 subcommands: index / qc / map / pile / scan / score / run.

## Subcommands and parameters

Conventions: `--threads` 0 = all cores (converted to the actual core count before dispatch); repeatable parameters may be comma-separated or given multiple times.

### index — build the paidx alignment index (→ `map::index_io::save`)

| flag | type | default | meaning |
|---|---|---|---|
| `--fasta` | PathBuf | **required** | reference FASTA |
| `--out` | PathBuf | **required** | output paidx index path |
| `--gtf` | PathBuf | optional | transcript annotation; also builds the L1 bundle (`<out stem>.bndl` + `<out stem>.tidx`) |
| `--k` | u32 | 15 | k-mer length (must match alignment time) |
| `--w` | u32 | 5 | window size (number of k-mers) |

The produced paidx connects directly to `run`/`map`'s `--index`; `rna_default`'s seed parameters are exactly k=15/w=5. Every build also writes `<out stem>.cpaidx` (collapsed-alphabet index, k=31/w=10) — the map-stage rescue index.

### qc — FASTQ quality control (→ `qc::run`)

| flag | type | default | meaning |
|---|---|---|---|
| `--r1` | Vec\<PathBuf\> | **required** | R1 input (multiple lanes merged in order) |
| `--r2` | Vec\<PathBuf\> | empty=SE | R2 (same length as r1) |
| `--out` | PathBuf | **required** | output directory |
| `--adapter-r1` / `--adapter-r2` | Vec\<String\> | built-in table | override adapter table |
| `--disable-adapter-trim` | bool | false | disable adapter trimming |
| `--disable-pe-overlap` | bool | false | disable PE overlap trimming |
| `--qtrim` / `--qtrim-cutoff` | bool / u8 | false / 20 | 3' quality trimming |
| `--trim-front1/2` `--trim-tail1/2` | usize | 0 | fixed-length front/back trimming |
| `--polyg` | auto\|on\|off | auto | polyG mode |
| `--min-len` | usize | 15 | minimum length |
| `--n-max` | usize | 5 | max N count |
| `--q15-frac-max` | f64 | 0.4 | max fraction of Q<15 |
| `--keep-unpaired` | bool | false | keep the passing mate when one PE mate fails |
| `--detect-adapter-se` | bool | false | SE adapter auto-detection |
| `--threads` | usize | 0 | threads |
| `--format` | fqgz\|bfq | fqgz | output format |

### map — alignment (→ `map::pipeline::run_se_2pass` / `run_pe_2pass`)

| flag | type | default | meaning |
|---|---|---|---|
| `--r1` / `--r2` | PathBuf | r1 **required**, r2 empty=SE | input FASTQ |
| `--index` | PathBuf | **required** | paidx index |
| `--gtf` | PathBuf | None | sjdb junction library |
| `--l1-bundle` | PathBuf | None | L1 engine bundle (`engine::L1Index::open`) |
| `--out` | PathBuf | **required** | output directory |
| `--threads` | usize | 0 | threads |

Behavior contract: `AlignConfig::rna_default()` + `extend.editing_aware=true`; RNA always 2-pass; `jkmer=None`. Fixed artifact names: `raw.bam` (streamed, unsorted; sorting happens in run's sort step) + `unmapped.fq.gz` + `align_qc.json` + `align.baln`.

### pile — pileup features (→ `pile::extract_pileup_features[_batch]`)

| flag | type | default | meaning |
|---|---|---|---|
| `--bam` | PathBuf | **required** | BAM (requires .bai/.csi) |
| `--chrom` + `--pos` | String / i64 | single-site mode | 1-based site |
| `--sites` | PathBuf | batch mode | `chrom\tpos` (1-based) file |
| `--out` | PathBuf | stdout | output TSV |

Single-site and batch are mutually exclusive; giving both is an error. Output columns: `chrom pos depth A_count C_count G_count T_count mean_base_quality strand_bias mean_mapq` (FEATURE_NAMES order, values Display); batch output row order = input order.

### scan — candidate-site discovery (→ `scan::run_call`)

| flag | type | default | meaning |
|---|---|---|---|
| `--bam` | PathBuf | **required** | input BAM |
| `--baln` | PathBuf | None | .baln dual source (takes precedence when given, ignores bam) |
| `--out` | PathBuf | **required** | candidates.bed |
| `--fasta` | PathBuf | None (refs probing by default; if still absent → majority pseudo-ref) | reference |
| `--gtf` | PathBuf | None | strand-decision evidence |
| `--gnomad` | PathBuf | None (refs probing by default) | soft down-weighting |
| `--lib` | unstranded\|stranded | unstranded | library strand specificity |
| `--enable-cu` | bool | false | C>U symmetric mode |
| `--min-call-score` | f64 | None | threshold (marks only, never deletes; contract owned by scan) |
| `--spec` | PathBuf | None (built-in v2) | scoring spec JSON |
| `--threads` | usize | 0 | threads |

### score — RE_PROB scoring (→ `score::pipeline::score_sites_batched`)

| flag | type | default | meaning |
|---|---|---|---|
| `--bam` | PathBuf | **required** | BAM (requires index) |
| `--sites` | PathBuf | **required** | `chrom\tpos` 1-based (pos ≥ 1, validated at library level) |
| `--fasta` | PathBuf | refs probing by default | reference (needs .fai) |
| `--bundle` | PathBuf | zero-config 5-level | model bundle root |
| `--caduceus` | PathBuf | resolved within bundle | encoder directory |
| `--out` | PathBuf | **required** | TSV `chrom\tpos\tprob` |
| `--threads` | usize | 0 | threads |
| `--batch` | usize | 64 | batch size |

### run — full-pipeline orchestration (→ `flow::run_pipeline`)

| flag | type | default | meaning |
|---|---|---|---|
| `--r1` / `--r2` | Vec\<PathBuf\> | entry derivation | FASTQ (multi-lane) |
| `--bam` | PathBuf | entry derivation | input BAM (must be sorted + indexed) |
| `--sites` | PathBuf | entry derivation | user sites (BamSites) |
| `--index` | PathBuf | required for FASTQ entries | paidx |
| `--fasta` | PathBuf | refs probing by default | reference (needs .fai) |
| `--gtf` | PathBuf | None (refs probing by default) | annotation |
| `--gnomad` | PathBuf | None (refs probing by default) | gnomAD |
| `--bundle` | PathBuf | zero-config 5-level | model bundle |
| `--caduceus` | PathBuf | resolved within bundle | encoder |
| `--l1-bundle` | PathBuf | None | L1 engine |
| `--lib` | unstranded\|stranded | unstranded | strand specificity |
| `--out` | PathBuf | **required** | output root |
| `--threads` | usize | 0 | threads |
| `--batch` | usize | 64 | score batch |

Entry derivation / species guardrail / stage paths / VCF contract all belong to the flow spec.

## Zero-config resolution

Everything follows "explicit flag first; probing only fills blanks". The first level that hits stops the search; if all levels miss → English error listing the attempted paths, exit code 1.

### bundle (5 levels, used by run / score)

1. `ESPERANTO_BUNDLE` environment variable;
2. `CARGO_MANIFEST_DIR/../../../bundle/human/esperanto-model-v1.4.1-501_40ep/rust` (source tree);
3. `<exe>/bundle/human/esperanto-model-v1.4.1-501_40ep/rust` and `<exe>/../bundle/human/esperanto-model-v1.4.1-501_40ep/rust` (release package);
4. `~/.local/share/esperanto/bundle/bundle/human/esperanto-model-v1.4.1-501_40ep/rust`;
5. No hit → error.

Validity at each level: the directory contains `norm.json`.

### refs (4 levels, used by run's fasta/gtf/gnomad and score's fasta)

1. `ESPERANTO_REFS`; 2. `<exe>/refs` and `<exe>/../refs`; 3. `~/.local/share/esperanto/refs`; 4. `./refs`.

Directory validity: contains a `*.fa` with a same-named `.fai`. File location: fasta = `hg38.fa` preferred, otherwise the first `*.fa` in lexicographic order; gtf = first `*.gtf`; gnomad = a `*.vcf.gz` starting with `dbsnp*` or `gnomad*`. **Only fills blanks; never overrides explicit values.**

### caduceus

When `--caduceus` is omitted → `score::pipeline::resolve_encoder_from_bundle(bundle)` (`bundle/encoder` → `bundle/../encoder`, requires `model.safetensors`).

## Exit codes and errors

- 0 success; 1 runtime error (downstream crate errors propagate verbatim through anyhow without swallowing semantics; zero-config failures list candidate paths); 2 argument error (clap).
- All user-facing strings are in English; `--help` on every subcommand exits 0.

## Determinism

cli itself writes no wall-clock/random quantities into any artifact; artifact byte contracts belong to each crate / flow.

## Out of scope (1.0.0)

- realign / report (not ported; defined in the flow spec); ui / web / fetch-refs / build-ref / bfq2fq (BACKLOG).
- config/preset YAML four-layer merge (BACKLOG); resume / --only / --stop-after (BACKLOG); `--keep-*` cleanup switches (BACKLOG); run's --name.
- The map subcommand does not emit sorted.bam (sorting belongs to run's sort step; use run when a sorted BAM is needed).

## Self-checks

- Every subcommand `--help` exits 0; illegal entry combinations (e.g., r1 and bam both given) exit 1 with a message identifying the fields.
- `ESPERANTO_BUNDLE` pointing to a valid bundle → score/run can run without `--bundle`; all 5 levels missing → error listing paths, exit 1.
- `ESPERANTO_REFS` pointing to a valid directory → run can run without `--fasta`; an explicit `--fasta` is never overridden.
- run FastqSe / BamSites and direct library calls with the same input → final artifacts byte-identical.
- `cargo clippy -p esperanto-cli -- -D warnings` zero warnings.
