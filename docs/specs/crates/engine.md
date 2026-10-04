# engine — transcriptome matching and projection

The optional L1 engine matches a read to indexed transcript sequences and projects the placement onto genomic coordinates. Reads are independent, and tied placements use deterministic ordering.

## Index contract

`L1Index::build(tidx_path, gtf_path, ref_path)` constructs the transcript sequence store and genomic projection map. `save` writes the bundle; `open` loads it together with the matching `.tidx` sidecar. The k-mer index, sequence store, and projection map share the same transcript identifier space, sorted by `transcript_id`.

L1 projection includes every indexed transcript. The filtered transcript map used for attribution must not replace that projection map. Transcript joins use identifiers, not display names.

## Matching

K-mers use case-insensitive two-bit A/C/G/T encoding and canonical forward/reverse-complement keys. Windows containing non-ACGT bases are skipped. The index supplies k, normally 31.

The engine visits read k-mers in position order and extends indexed hits. Repeated hits on the same transcript, strand, and diagonal are deduplicated. `EngineConfig` defaults are a 10-base flank, 0.90 coverage threshold, 8 transcript candidates, 4 diagonals per transcript, 16384 raw hits, and 100 hits per k-mer.

Extension and Myers verification recognize reference A/read G and reference T/read C as editing-consistent matches, in addition to identical A/C/G/T bases. Reverse conversions remain mismatches. Input case does not alter matching.

Full-length extensions are projected directly. Interrupted extensions are ranked by covered length, transcript identifier, and diagonal; up to eight candidates undergo Myers infix verification. Accepted interrupted placements require coverage at least 0.90 and distance no greater than `max(read_length / 33, 1)`. Reads shorter than k or without an accepted candidate fall back to genomic alignment.

## Placement and ambiguity

Transcript placements are converted to forward genomic coordinates and CIGAR operations. Read strand is combined with transcript strand. Equivalent projections are treated as the same genomic placement; equally supported distinct loci receive MAPQ 0. Marginal placements competing with another accepted locus also receive MAPQ 0. The mapper incorporates L1 evidence into its final alignment decision.

## Bundle format

The bundle begins with `L1BNDL01`, followed by little-endian version 1, a length-prefixed transcript map, and transcript sequences. Sequence entries contain a length and sequence bytes. Contig names are stored in the transcript-map component. Identical input produces identical bundle content.

## Errors and determinism

Invalid bundle magic or version, malformed records, inconsistent transcript counts, and I/O failures are reported as errors. All ordering that affects placements is explicit; worker scheduling does not choose among tied candidates.
