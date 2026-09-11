# flow — pipeline orchestration (FASTQ PE/SE/BAM entries, species guardrail)

Single-library entry `run_pipeline`: orchestrates qc → map → sort → scan → score → vcf → report by entry type; artifacts land in `<out_dir>/<stage>/`; final products `sites.vcf` and `<sample>.report.html`. Orchestration carries zero scientific semantics: all numeric contracts belong to the stage crates; flow only does wiring, fail-fast, and deterministic persistence.

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

Before `run_pipeline` starts any stage: read `<fasta>.fai`. **With a `species.json` manifest** in the run directory (hybrid runs): the declared baseline decides the expected chr1 length (`grch38` → 248956422, `grcm39` → 195154007, `mm10` → 195471971); human loci are verified at reference build time. **Without a manifest (legacy)**: if a `chr1` line exists, its length must be `== 248956422` (hg38) or `< 10_000_000` (synthetic/test reference), otherwise refuse with `FlowError::SpeciesMismatch{len}`; no `chr1` line passes directly. Same rule as score's internal guardrail (score spec) as a double safeguard — flow guards before burning compute, score guards before scoring. The error message is in English and gives the actual length, the hg38 expected value, and a bundle-matching hint.

## Parameters (RunParams)

`r1: Vec<PathBuf>` / `r2: Vec<PathBuf>` (qc natively merges multiple lanes), `bam` / `sites`, `index` (paidx, required for FASTQ entries), `fasta` (required, needs `.fai`), `gtf` (optional; junction library for map, gene annotation for the report), `snp` (optional known-SNP VCF for scan), `bundle` (score bundle root), `caduceus` (optional; defaults to in-bundle resolution), `l1_bundle` (optional, `engine::L1Index::open`), `out_dir`, `threads`, `batch` (score batch, default 256).

## Stage wiring (contract paths)

| Stage | Call | Artifacts |
|---|---|---|
| qc | `qc::run(QcParams{r1,r2,out_dir:<out>/qc,…})`, OutFormat::Fqgz | `<stem>.clean[_R1/_R2].fq.gz` (qc naming contract) + `qc.json` |
| map | **Tier 1**: `index_io::load(paidx)`; `AlignConfig::rna_default()` taken unchanged, so `extend.edit_dir = EditDir::Off` — strict four-letter alignment, no editing exemption; `jlib = gtf.map(gtf::from_gtf…)`; `jkmer = None` (not wired in 1.0.0, see BACKLOG); `l1 = l1_bundle.map(L1Index::open)`; single-pass when L1 present, else 2-pass. **Tier 2** (collapsed-alphabet rescue) runs when the `ag.paidx`/`tc.paidx` transformed indexes exist beside `paidx`. Its input is every read Tier 1 refused: `unmapped.fq.gz` plus the primary placements leaving more than `1 - RESCUE_CLIPPED_COVERAGE` (30%) of the read soft-clipped, which Tier 1 writes as *mapped* and which therefore never reach `unmapped.fq.gz`. Artifact reads are screened out first (composition/quality gates). Survivors are realigned in **four passes** — `RESCUE_DIRS`: A>G against each index, T>C against each index. Both indexes are needed per transform because `revcomp(t_AG(revcomp(F))) == t_TC(F)`: an A>G-edited read antisense to the genome only seeds against the T>C index, so with same-name passes alone every placement is plus-strand (measured 7987/7987) and no hyperedited fragment could place both mates. Each pass keeps `edit_dir = Off` — under the matching transform that direction's edits are already identities, and exempting further would re-open the collapse the pass exists to avoid. Every placement is then verified against the true four-letter reference (`verify_rescue_placement`, gates in order): well-formed read → known contig → **ungapped CIGAR** (measured: all 982 re-testable pure-indel placements went on to fail the mismatch gate, so the gaps are real) → on-contig span → at most `RESCUE_MAX_NON_EDIT_MM` (2) non-editing-class mismatches → at least one edit → editing-class density ≥ `RESCUE_MIN_EDIT_DENSITY` (5%) of aligned length → editing-class purity >60% (>80% for reads ≤60 bp) at Phred ≥30, cluster span ≥10% of read length not contained in the outer 20%, cluster single-nucleotide share ≤60%. A read that verifies under **two** directions is direction-ambiguous and dropped rather than arbitrated — it matches a paralogue, not an edited transcript. Paired runs then apply mate concordance once, over the four-letter placements plus the rescue's own. Verified reads are appended to raw.bam with `rescued = true`, which caps MAPQ at `PILEUP_GATE` (20) via `mapq::mapq_of` rather than forcing 0 — the transformed index's own margin is the only uniqueness evidence the placement has, and zeroing it put the read below the caller's gate, making the rescue unable to produce a call — plus `RE:Z:collapsed`; align_qc.json gains `rescued_collapsed` (verified), `rescue_rejected_collapsed`, `rescue_prefiltered_collapsed`; unmapped.fq.gz is rewritten without the verified reads; placements also land in `rescued.bed`, one `chrom\tpos` 0-based row each, for the report's hyperedited-region track | `<out>/map/raw.bam` + `unmapped.fq.gz` + `align_qc.json` + `align.baln` + `rescued.bed` |
| sort | `bamio::sort::coordinate_sort(raw.bam → sorted.bam)` + builds `.bai` (see bamio spec addendum) | `<out>/map/sorted.bam(.bai)` |
| scan | `scan::run_call(CallParams::standard(sorted.bam, candidates.bed, fasta))` with `snp` = the run's SNP VCF: MAPQ >= 20, BQ >= 30, depth >= 10, >= 2 editing-consistent reads, frequency >= 1%, known SNPs removed (scan spec). Bam entries read the input BAM. No further candidate filtering happens in flow | `<out>/scan/candidates.bed` |
| score | bed→sites bridge (next section) → `score::pipeline::score_sites_batched(bam=sorted.bam,…)` (pileup features over reads with MAPQ >= 20, the same read set scan counted); BamSites entries parse `--sites` directly (same bridge). Hybrid runs (species.json): sites route per contig owner — human-locus sites score with the human bundle, mouse-contig sites with the mouse bundle declared in the manifest (`ReferenceCheck::TrustedHybrid` skips the chr1 heuristic for both); with no mouse bundle in the manifest, mouse-contig sites reach the VCF as `FILTER=UNSCORED` with VAF/DEPTH kept and `NA` in scores.tsv | `<out>/score/scores.tsv` (`chrom\tpos\tprob`, row order = input order) |
| vcf | merge candidates.bed + scores.tsv (next section) | `<out>/sites.vcf` |

bed→sites bridge: column 1 of candidates.bed = chrom, column 3 (pos0+1) = 1-based pos; the `--sites` file is parsed as `chrom\tpos` (1-based). Empty candidates = empty scores.tsv + header-only vcf; no error.

## Resume (run.json + stage walk)

Every `run` lands in a sample-scoped directory `<out>/<sample>/` (`--sample` or derived from the R1/BAM file name: strip `.gz/.fastq/.fq/.bam`, then a trailing `_R1/_R2/_1/_2`). Before qc starts, flow atomically writes `<out>/<sample>/run.json` (tmp + rename):

- `version` (schema, currently 1), `esperanto` (binary version), `sample`, `entry` (`fastq-se|fastq-pe|bam|bam-sites`);
- `inputs`: per input file `{role, path, size, mtime}` (mtime in unix seconds);
- `params`: fully resolved `index/fasta/gtf/snp/bundle/caduceus/l1_bundle/threads/batch/device` (`snp` also accepts the older `gnomad` key when reading).

`run` refuses a directory that already holds `run.json` (points at `resume`). `esperanto resume <sample-dir>` needs no flags: it reads run.json, refuses on a changed input (size or mtime mismatch, naming the file) or a missing parameter path, notes an esperanto version mismatch on stderr and proceeds (determinism guarantee), then validates stage artifacts in pipeline order and re-executes from the first invalid stage, wiping that stage's artifacts first (later stages cascade; earlier stages are never touched).

Artifact validators (integrity, not presence): BAM = trailing BGZF EOF marker; plain-gzip streams (qc clean reads, `unmapped.fq.gz`) = full multi-member decompression succeeds; JSON = parses; `candidates.bed` = parses through the bed→sites bridge; `scores.tsv` = every data line is `chrom\tpos\tprob` with parseable numbers (empty is legal); `sites.vcf` = starts with `##fileformat=VCF`; `<sample>.report.html` = non-empty (skipped entirely when no GTF). Stage artifacts validated per stage: qc = `qc.json` + `qc.html` + exactly one clean file per mate; map = `raw.bam` + `align_qc.json` + `unmapped.fq.gz` + `align.baln`; sort = `sorted.bam` + non-empty `sorted.bam.bai`; scan/score/vcf/report as above.

Map seal: intact `raw.bam` + `align_qc.json` means the alignment completed; when a `.cpaidx` exists and `align_qc.json` lacks the `rescued_collapsed` key while `unmapped.fq.gz` is non-empty, resume re-runs only the collapsed rescue (never the alignment). An empty `unmapped.fq.gz` counts as rescue-done (the rescue is a no-op there by definition).

Concurrency: `<sample>/.lock` is created exclusively before any artifact is read or written (run and resume both); a held lock fails fast with the lock path; the lock is removed on exit. A stale lock from a killed process is removed by hand (the error says so).
## Output VCF

VCF v4.2; `##reference=<fasta file name>`; CHROM=chrom, POS=1-based, ID=`.`, REF=reference base (from candidates.bed; for BamSites entries from the fasta, uppercase, `N` if not found), ALT=the editing-consistent alternate (G for REF A, C for REF T, `N` otherwise), QUAL=`.`; FILTER: `RE_PROB ≥ 0.5` → `PASS`, otherwise `LOW_SCORE`, `UNSCORED` for hybrid mouse-contig sites without a model; INFO: `RE_PROB` (score probability, Display), `VAF` (editing-consistent alternate frequency), `DEPTH`, `AC` (counted A,C,G,T reads) passed through from candidates.bed, `STRAND` (`+` for REF A, `-` for REF T). Row order = candidates.bed order (already sorted by (chrom,pos)). Contigs without candidates emit no rows.

## Determinism

- flow itself writes no wall-clock/random quantities into any artifact.
- Two runs with the same input: `candidates.bed` / `scores.tsv` / `sites.vcf` byte-identical (sort-stability contract, see bamio addendum).
- Thread count (1 vs 8) does not change any artifact bytes.

## Errors

`FlowError` (thiserror): `Entry(String)` / `SpeciesMismatch{len}` / `MissingBamIndex{path}` / `BedParse{line,msg}` / `Io` / `Stage{stage:&'static str, source:Box<dyn Error+Send+Sync>}` (wraps downstream crate errors without swallowing semantics).

## Dependencies

esperanto-qc / map / scan / score / bamio / engine (L1Index) / report; rust-htslib (bai validation probe), thiserror, serde/serde_json (run.json). No new workspace dependencies.

## Out of scope

- `--only` / `--stop-after` stage selection (BACKLOG).
- jkmer (Track 2) index construction and wiring (BACKLOG; L1 + 2-pass is the main path).
- Automatic sorting of user BAMs (missing index errors out; flow does not sort on the user's behalf).

## Self-checks

- Species guardrail: chr1=249250621 (hg19) → run refused and no stage directory created; chr1 < 10 Mb synthetic reference → allowed.
- Bam entry missing .bai → `MissingBamIndex`, no intermediates produced.
- FastqSe/FastqPe full chain: sorted.bam + .bai exist and are consumed normally by scan/score.
- Empty candidates.bed → empty scores.tsv, header-only sites.vcf, exit code 0.
- Same input twice + threads 1 vs 8: the three final artifacts byte-identical; `cargo clippy -p esperanto-flow -- -D warnings` zero warnings.
