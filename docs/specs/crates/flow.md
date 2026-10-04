# flow — pipeline orchestration

`run_pipeline` orchestrates a single library and writes its artifacts under the sample directory.

## Entry types

| Input | Stage sequence |
|---|---|
| `--r1`, optionally `--r2` | QC, mapping, sorting, discovery, scoring, VCF, report |
| `--bam` | Discovery, scoring, VCF, report |
| `--bam --sites` | Scoring, VCF, report |

FASTQ inputs may contain multiple lanes, merged in the supplied order. BAM inputs must be coordinate-sorted and indexed. A report is generated when a matching GTF is available. Invalid combinations of FASTQ, BAM, and supplied sites are rejected.

## Inputs and parameters

`RunParams` contains input paths, reference FASTA, optional GTF, alignment index, optional L1 bundle, model bundle, output directory, worker count, and optional stopping stage. The FASTA requires a matching `.fai`. For a reference containing `chr1`, its length must equal 248956422, or be below 10000000 for a small reference; an incompatible length is an error.

The CLI writes each run under `<out>/<sample>`. The sample defaults to the input filename with compression, format, and mate suffixes removed. An existing `run.json` causes `run` to refuse the directory and direct the user to `resume`.

## Stage artifacts

| Stage | Operation | Artifacts |
|---|---|---|
| QC | Trim and filter reads | `qc/*.clean*.fq.gz`, `qc/qc.json`, `qc/qc.html` |
| Mapping | Splice-aware alignment and supported read rescue | `map/raw.bam`, `map/unmapped.fq.gz`, `map/align_qc.json`, `map/discovered_junctions.bed` |
| Sorting | Stable coordinate sort and BAM indexing | `map/sorted.bam`, `map/sorted.bam.bai` |
| Discovery | `score::unified::discover_sites` | `scan/candidates.bed` |
| Scoring | `score::pipeline::score_sites_batched` | `score/scores.tsv` |
| VCF | Combine candidates and scores | `sites.vcf` |
| Report | Reference and GTF annotation | `<sample>.report.html` |

Mapping uses `AlignConfig::rna_default()` and the GTF junction library when provided. The primary aligner uses strict four-letter matching. A complete optional L1 bundle supplies transcriptome placements. When the A-to-G and T-to-C transformed indexes are present, collapsed-read rescue validates placements against the untransformed reference before accepting them. Rescue annotations describe alignment evidence; they do not add a separate site-scoring model.

Discovery reports reference A/T sites with at least one editing-consistent observation after MAPQ 20 and BQ 13 filtering. BED columns 1 and 3 provide chromosome and 1-based position for scoring. Supplied sites are read directly as `chrom<TAB>pos`. Empty input produces an empty score table and a header-only VCF.

## VCF

VCF version 4.2. REF is the reference base; ALT is G for reference A and C for reference T. `PASS` means the score meets the threshold in `freeze.json`; `LOW_SCORE` means it does not. INFO fields are `ESP_SCORE`, `VAF`, `DEPTH`, `AC`, and `STRAND`. `ESP_SCORE` is a classifier score. `AC` contains A, C, G, T counts in that order. `STRAND` is `+` for reference A and `-` for reference T.

Discovered candidates retain their chromosome/position order. Supplied-site entries retain their input order. For supplied sites, the reference supplies REF/ALT; candidate count fields are unavailable.

## Resume

`run.json` uses schema version 2 and records the executable version, sample, entry type, input fingerprints, and resolved parameters. Each input fingerprint contains path, size, and modification time. Incompatible schema versions are rejected; cached outputs from such directories are not reused.

`resume` verifies the inputs and required reference/model paths, then walks the pipeline in stage order. It reuses completed stages and restarts at the first incomplete or invalid stage. A requested stopping stage leaves a resumable directory.

Artifact validation includes BAM completion and index availability, complete gzip decompression, parseable JSON, candidate and score table structure, VCF headers, and nonempty HTML reports when required. Collapsed rescue additionally validates its completion records against the merged BAM and rescued-read table.

An exclusive `<sample>/.lock` prevents concurrent access by `run` and `resume`. Failure to acquire the lock reports the path. The lock is removed on normal exit.

## Errors and determinism

Missing inputs, malformed sites, incompatible reference lengths, missing BAM indexes, and stage failures propagate as errors. Candidate tables, scores, and VCF records are byte-stable for identical inputs and model files across worker counts. QC and run metadata may contain elapsed times and file timestamps.
