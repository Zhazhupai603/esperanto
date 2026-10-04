# scan — unified candidate discovery

The `scan` command and the pipeline discovery stage both call `esperanto_score::unified::discover_sites`. Discovery consumes a coordinate-sorted, indexed BAM and a matching reference FASTA with `.fai`.

## Observation filters

- Skip unmapped, secondary, and supplementary alignments.
- Require MAPQ at least 20 and base quality at least 13.
- Exclude missing base qualities and non-ACGT observations.
- Compare stored BAM bases directly to the forward reference over aligned CIGAR columns.
- Deletions and reference skips do not contribute observations; inserted and soft-clipped bases do not consume reference coordinates.

A reference A is emitted when at least one G is observed; a reference T is emitted when at least one C is observed. Candidate discovery imposes no additional depth or frequency threshold. It does not exclude known variant positions.

## Output

The BED contains chromosome, 0-based start, 1-based end, reference base, depth, A count, C count, G count, T count, and alternate frequency. An empty trailing field is permitted. Alternate frequency is written to eight decimal places.

Depth and base counts are read-level observations. The scorer separately merges observations by template. Rows are sorted by chromosome and position. Empty input yields an empty file.

Discovery processes reference blocks in parallel and restores deterministic output ordering. BAM contigs absent from the FASTA index are skipped.
