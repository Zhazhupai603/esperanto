# flow — pipeline orchestration (FASTQ PE/SE/BAM entries, species guardrail)

Single-library entry `run_pipeline`: orchestrates qc → map → sort → scan → score by entry type; artifacts land in `<out_dir>/<stage>/`; final product `sites.vcf`. Orchestration carries zero scientific semantics: all numeric contracts belong to the stage crates; flow only does wiring, fail-fast, and deterministic persistence. 1.0.0 excludes realign (not ported, see BACKLOG) and the HTML report (deferred to 1.1).

## Entry derivation (Entry)

Pure field branching, `RunParams::entry()`:

| Field combination | Entry | Stage sequence |
|---|---|---|
| `r1` non-empty, `r2` empty | FastqSe | qc → map → sort → scan → score → vcf |
| `r1` + `r2` both non-empty | FastqPe | same as above |
| `bam` present, `sites` absent | Bam | scan → score → vcf |
| both `bam` + `sites` present | BamSites | score → vcf |
| anything else (r1 and bam both given, sites without bam, all empty) | — | `FlowError::Entry` |

Bam entry contract: the input BAM must be coordinate-sorted and have `.bai`/`.csi`; missing index = `FlowError::MissingBamIndex` (flow does not silently copy-and-sort for the user).

## Species guardrail (hard contract, fail fast at entry)

Before `run_pipeline` starts any stage: read `<fasta>.fai`; **if a `chr1` line exists**, its length must be `== 248956422` (hg38) or `< 10_000_000` (synthetic/test reference), otherwise refuse with `FlowError::SpeciesMismatch{len}`; no `chr1` line passes directly. Same rule as score's internal guardrail (score spec) as a double safeguard — flow guards before burning compute, score guards before scoring. The error message is in English and gives the actual length, the hg38 expected value, and a bundle-matching hint.

## Parameters (RunParams)

`r1: Vec<PathBuf>` / `r2: Vec<PathBuf>` (qc natively merges multiple lanes), `bam` / `sites`, `index` (paidx, required for FASTQ entries), `fasta` (required, needs `.fai`), `gtf` / `gnomad` (optional), `bundle` (score bundle root), `caduceus` (optional; defaults to in-bundle resolution), `l1_bundle` (optional, `engine::L1Index::open`), `lib: LibType` (passed through to scan), `out_dir`, `threads`, `batch` (score batch, default 64).

## Stage wiring (contract paths)

| Stage | Call | Artifacts |
|---|---|---|
| qc | `qc::run(QcParams{r1,r2,out_dir:<out>/qc,…})`, OutFormat::Fqgz | `<stem>.clean[_R1/_R2].fq.gz` (qc naming contract) + `qc.json` |
| map | `index_io::load(paidx)`; `AlignConfig::rna_default()` + `extend.editing_aware=true`; `jlib = gtf.map(gtf::from_gtf…)`; `jkmer = None` (not wired in 1.0.0, see BACKLOG); `l1 = l1_bundle.map(L1Index::open)`; single-pass when L1 present, else 2-pass; then collapsed-alphabet rescue of unmapped reads when `<paidx stem>.cpaidx` exists (survivors appended to raw.bam with MAPQ 0 + `RE:Z:collapsed`, unmapped.fq.gz rewritten; placements also land in `rescued.bed`, one `chrom\tpos` 0-based row each, for the report's hyperedited-region track) | `<out>/map/raw.bam` + `unmapped.fq.gz` + `align_qc.json` + `align.baln` + `rescued.bed` |
| sort | `bamio::sort::coordinate_sort(raw.bam → sorted.bam)` + builds `.bai` (see bamio spec addendum) | `<out>/map/sorted.bam(.bai)` |
| scan | `scan::run_call(CallParams{…})`; FASTQ entries (SE and PE) use `baln = align.baln` when it is a non-stub file (byte-identical dual-source contract with BAM), falling back to the sorted BAM otherwise; Bam entries always read the input BAM | `<out>/scan/candidates.bed` |
| score | bed→sites bridge (next section) → `score::pipeline::score_sites_batched(bam=sorted.bam,…)`; BamSites entries parse `--sites` directly (same bridge) | `<out>/score/scores.tsv` (`chrom\tpos\tprob`, row order = input order) |
| vcf | merge candidates.bed + scores.tsv (next section) | `<out>/sites.vcf` |

bed→sites bridge: column 1 of candidates.bed = chrom, column 3 (pos0+1) = 1-based pos; the `--sites` file is parsed as `chrom\tpos` (1-based). Empty candidates = empty scores.tsv + header-only vcf; no error.

## Resume (run.json + stage walk)

Every `run` lands in a sample-scoped directory `<out>/<sample>/` (`--sample` or derived from the R1/BAM file name: strip `.gz/.fastq/.fq/.bam`, then a trailing `_R1/_R2/_1/_2`). Before qc starts, flow atomically writes `<out>/<sample>/run.json` (tmp + rename):

- `version` (schema, currently 1), `esperanto` (binary version), `sample`, `entry` (`fastq-se|fastq-pe|bam|bam-sites`);
- `inputs`: per input file `{role, path, size, mtime}` (mtime in unix seconds);
- `params`: fully resolved `index/fasta/gtf/gnomad/bundle/caduceus/l1_bundle/lib/threads/batch/device`.

`run` refuses a directory that already holds `run.json` (points at `resume`). `esperanto resume <sample-dir>` needs no flags: it reads run.json, refuses on a changed input (size or mtime mismatch, naming the file) or a missing parameter path, notes an esperanto version mismatch on stderr and proceeds (determinism guarantee), then validates stage artifacts in pipeline order and re-executes from the first invalid stage, wiping that stage's artifacts first (later stages cascade; earlier stages are never touched).

Artifact validators (integrity, not presence): BAM = trailing BGZF EOF marker; plain-gzip streams (qc clean reads, `unmapped.fq.gz`) = full multi-member decompression succeeds; JSON = parses; `candidates.bed` = parses through the bed→sites bridge; `scores.tsv` = every data line is `chrom\tpos\tprob` with parseable numbers (empty is legal); `sites.vcf` = starts with `##fileformat=VCF`; `<sample>.report.html` = non-empty (skipped entirely when no GTF). Stage artifacts validated per stage: qc = `qc.json` + `qc.html` + exactly one clean file per mate; map = `raw.bam` + `align_qc.json` + `unmapped.fq.gz` + `align.baln`; sort = `sorted.bam` + non-empty `sorted.bam.bai`; scan/score/vcf/report as above.

Map seal: intact `raw.bam` + `align_qc.json` means the alignment completed; when a `.cpaidx` exists and `align_qc.json` lacks the `rescued_collapsed` key while `unmapped.fq.gz` is non-empty, resume re-runs only the collapsed rescue (never the alignment). An empty `unmapped.fq.gz` counts as rescue-done (the rescue is a no-op there by definition).

Concurrency: `<sample>/.lock` is created exclusively before any artifact is read or written (run and resume both); a held lock fails fast with the lock path; the lock is removed on exit. A stale lock from a killed process is removed by hand (the error says so).
## Output VCF (1.0.0 minimal contract)

VCF v4.2; `##reference=<fasta file name>`; CHROM=chrom, POS=1-based, ID=`.`, REF=reference base at that position from fasta (uppercase, `N` if not found), ALT=symbolic allele `<RE>` (precise ALT allele inference deferred to 1.1, BACKLOG), QUAL=`.`; FILTER: `RE_PROB ≥ 0.5` → `PASS`, otherwise `LOW_SCORE`; INFO: `RE_PROB` (score probability, Display), `VAF`/`DEPTH`/`STRAND`/`EVID` (passed through from candidates.bed). Row order = candidates.bed order (already sorted by (chrom,pos)). Contigs without candidates emit no rows.

## Determinism

- flow itself writes no wall-clock/random quantities into any artifact; no run_state/resume (out of scope for 1.0.0, BACKLOG).
- Two runs with the same input: `candidates.bed` / `scores.tsv` / `sites.vcf` byte-identical (sort-stability contract, see bamio addendum).
- Thread count (1 vs 8) does not change any artifact bytes.

## Errors

`FlowError` (thiserror): `Entry(String)` / `SpeciesMismatch{len}` / `MissingBamIndex{path}` / `BedParse{line,msg}` / `Io` / `Stage{stage:&'static str, source:Box<dyn Error+Send+Sync>}` (wraps downstream crate errors without swallowing semantics).

## Dependencies

esperanto-qc / map / scan / score / bamio / engine (L1Index); rust-htslib (bai validation probe), thiserror. No new workspace dependencies.

## Out of scope (1.0.0)

- run_state.json / resume / --only / --stop-after (BACKLOG).
- jkmer (Track 2) index construction and wiring (BACKLOG; L1 + 2-pass is the main path).
- realign / corrected_sites (not ported), HTML report (deferred to 1.1), VCF precise ALT alleles (BACKLOG).
- Automatic sorting of user BAMs (missing index errors out; flow does not sort on the user's behalf).

## Self-checks

- Species guardrail: chr1=249250621 (hg19) → run refused and no stage directory created; chr1 < 10 Mb synthetic reference → allowed.
- Bam entry missing .bai → `MissingBamIndex`, no intermediates produced.
- FastqSe/FastqPe full chain: sorted.bam + .bai exist and are consumed normally by scan/score.
- Empty candidates.bed → empty scores.tsv, header-only sites.vcf, exit code 0.
- Same input twice + threads 1 vs 8: the three final artifacts byte-identical; `cargo clippy -p esperanto-flow -- -D warnings` zero warnings.
