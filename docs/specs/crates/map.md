# map — splice-aware read alignment

The mapper aligns single-end or paired-end FASTQ reads to a `.paidx` genomic index and writes BAM alignments. RNA alignment uses minimizer seeds, chaining, splice-aware extension, and paired-end evidence.

## Inputs

- Plain or gzip-compressed FASTQ with synchronized paired-end records.
- A `.paidx` built from the matching reference; default minimizer parameters are k=15 and w=5.
- An optional GTF junction library.
- An optional complete transcriptome L1 bundle and `.tidx` sidecar built from the same reference.

`AlignConfig::rna_default()` supplies the production parameters. Primary alignment uses strict four-letter matching. Junction discovery and a second alignment pass support unannotated splice junctions. Annotation and discovered junctions are combined deterministically.

## Coordinates and orientation

Alignment positions and reference intervals use 0-based coordinates. CIGAR operations are `M`, `I`, `D`, `N`, and `S`; skipped introns use `N`. BAM sequence for a mapped minus-strand record is reverse-complemented and its qualities are reversed. Unmapped reads retain their original orientation.

## Placement evidence

Candidate placement ranking considers alignment score, alternative locations, splice structure, and paired-end agreement. MAPQ is bounded by 60. Ambiguous placements remain representable with low MAPQ; downstream observations require MAPQ at least 20.

The aligner can recover a read using an anchored mate. Alignment rescue and its `RE` tags are read-level evidence and do not bypass the downstream model's filters. Decoy reference contigs use the `decoy_` prefix and compete during alignment.

## Collapsed-read rescue

When `.ag.paidx` and `.tc.paidx` are present beside the primary index, the pipeline can realign unmapped and insufficiently covered reads using separate A-to-G and T-to-C transformations. These indexes use k=31 and w=10 and exclude decoy contigs.

Placements must pass verification against the true four-letter reference. Required evidence includes an ungapped alignment, at most two non-editing mismatches, editing density at least 5%, editing-class purity above 60% (above 80% for reads of at most 60 bases), high-quality edited bases, and internal cluster geometry. Paired-end placements additionally require mate concordance. Accepted records carry `RE:Z:collapsed`; rescue outputs are sealed only after the merged BAM is complete.

## Outputs

`raw.bam` is unsorted alignment output. `unmapped.fq.gz` contains unresolved reads. `align_qc.json` records total and mapped reads, proper pairs, insert-size summaries, MAPQ histogram, junction counts, read-rescue counts, and elapsed time. `discovered_junctions.bed` reports the junction-discovery result. The flow crate sorts and indexes BAM output for discovery and scoring.

## Index format

The `.paidx` magic is `PAIDXfmt` and the format version is 1. The file stores seed parameters, reference sequences and contig metadata, minimizer keys, and positions. Loading validates the file structure and integrity information before alignment.

## Errors and determinism

Malformed FASTQ, inconsistent mate records, invalid indexes, and I/O errors are reported with input context. Placement decisions and record order are deterministic. Worker scheduling does not choose between tied placements or alter scientific result ordering.
