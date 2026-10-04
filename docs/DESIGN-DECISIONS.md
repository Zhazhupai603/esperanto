# Scientific contracts

## Reference and coordinates

The released model targets human GRCh38/hg38. Reference FASTA, alignments, indexes, and annotations must use matching contig names and coordinates. BAM and BED positions use the standard 0-based coordinate convention; supplied site lists, score tables, and VCF positions are 1-based.

## Candidate discovery

Discovery and scoring use the same BAM observation filters: MAPQ at least 20, base quality at least 13, and exclusion of unmapped, secondary, and supplementary alignments. A reference A with an observed G, or a reference T with an observed C, is a candidate. Candidate discovery does not impose a minimum depth, alternate-read count beyond one, or alternate frequency.

## Site scoring

The model reads the BAM, reference sequence, and GRCh38 Alu annotation. Its 27 features describe template-level support, quality, alignment geometry, neighborhood evidence, sequence context, and Alu overlap. Overlapping mates contribute a single template observation during scoring. Discovery base counts are read-level counts and need not equal template-level model depth.

The released bundle fixes feature semantics, model trees, and the decision threshold. `ESP_SCORE` is a classifier score. `PASS` and `LOW_SCORE` are determined only by the bundle threshold.

## Alignment orientation

Mapped minus-strand BAM records store reverse-complemented sequence and reversed qualities. Downstream discovery and scoring compare these stored bases directly to the forward reference. Reference A implies the plus editing direction; reference T implies the minus editing direction.

## Reproducibility

Scientific result tables have deterministic ordering and values across worker counts. Run metadata records input fingerprints and parameters so that `resume` can validate an interrupted run before reusing its artifacts.
