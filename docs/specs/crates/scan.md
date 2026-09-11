# scan — A-to-I candidate site discovery

Genome-wide pileup from a coordinate-sorted BAM, emitting `candidates.bed`. The counting convention is the standard REDItools one (MAPQ >= 20, base quality >= 30, no read-orientation split), so candidate tables are directly comparable across pipelines. The only engine is the block scatter.

## Inputs

- `--bam`: aligned BAM (requires `.bai`/`.csi`). Skips unmapped (0x4) / secondary (0x100) / supplementary (0x800) records and records with MAPQ < `min_mapq`. Reads tagged `RE:Z:collapsed` are MAPQ 0 and therefore never counted.
- `--fasta` (required): reference; contigs absent from the `.fai` are skipped.
- `--snp` (optional): known-SNP VCF (plain or gzip). Only CHROM/POS are read; every data row is an excluded position. Contig names match with or without a `chr` prefix.
- Thresholds: `min_mapq` 20, `min_bq` 30, `min_depth` 10, `min_edit_reads` 2, `min_edit_freq` 0.01.

## Counting (hard contract)

Per the SAM specification, SEQ of a reverse-strand record is stored reference-forward, so bases are compared to the reference as stored, with no reverse complement and no orientation split. Per position: `depth` = counted bases (BQ >= `min_bq`, base in ACGT, reference base in ACGT) over CIGAR M/=/X columns; D/N consume reference only, I/S consume query only. Mismatching bases are tallied per base; the reference count is `depth − mismatches`.

## Candidate rule

A position is a candidate when its reference base is A or T, the editing-consistent alternate (G for A, C for T) is carried by >= `min_edit_reads` counted reads, `depth >= min_depth`, `alt / depth >= min_edit_freq`, and the position is not in the SNP set. Sites whose reference base is C or G cannot be A-to-I edited and never appear. The transcript strand is implied by the reference base (A → `+`, T → `-`) and is not computed by scan.

## No strand-consistency filter

A-to-I editing on one transcript is single-stranded, and a filter dropping candidates whose read-length neighbourhood mixed both directions was implemented on that basis. It is removed: on real data it deleted 80.3% of esperanto's candidates against 17.5% of STAR's, keeping 41.78% of esperanto's REDIportal calls against 83.31% of STAR's. That asymmetry means it was compensating for the aligner, not for the biology — mixing was being manufactured upstream by exempting both editing directions at once, which `map`'s identity term now charges at its source. Scan applies no positional filtering; the candidate rule above is the whole rule.

## Output (10-column contract)

`chrom  pos0  pos0+1  ref  depth  A  C  G  T  edit_freq(%.4f)` — `ref` is the reference base (A or T); base counts are integers over the counted reads; `edit_freq` = alternate reads / depth. Sorted by (chrom, pos0) before writing → thread-count independent, byte-identical. An empty contig yields zero rows.

`CallStats`: `candidates`, `below_threshold` (positions with >= 1 editing-consistent read that failed a threshold), `snp_removed`.

## Parallelism

Tasks = (contig, 32 Mbp block), rayon over both levels; the contig sequence is loaded once per contig (uppercase) and shared by its blocks. Every candidate is decided from its own column, so the result is independent of thread count.

## Dependencies

rust-htslib (BAM + faidx), thiserror, rayon, flate2 (gzip SNP VCF).

## Self-checks

- On a STAR BAM the candidate set and every per-site `ref/depth/A/C/G/T` equal the REDItools `analyze -q 20 -bq 30` table filtered by the same rule (verified: 89,837 sites, all identical).
- Same input twice, threads 1 vs 8: `candidates.bed` byte-identical.
- A reverse-strand read carrying a stored G at a reference A votes as A>G (no complement flip).
