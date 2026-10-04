# score — unified molecular evidence classifier

The scorer evaluates sites from a coordinate-sorted, indexed BAM using a reference-aware LightGBM model. It returns one classifier score per input site in the original order.

## Model bundle

The bundle contains `model.txt`, `freeze.json`, and `alu.bed`. The released model has 400 trees and 27 features. `freeze.json` provides the decision threshold. All three files must be present.

## Features and observations

Reads with unmapped, secondary, or supplementary flags, MAPQ below 20, or base quality below 13 are excluded. Missing base qualities are excluded. Reads are grouped by read-group and template name. Overlapping mates contribute the higher-quality base; equal-quality conflicting bases are ambiguous.

Features describe template depth and alternate support, base and mapping qualities, read-end position, soft clipping, alignment gaps, strand balance, likelihood evidence, exponentially weighted neighboring observations within 200 bases, template reuse, adjacent reference bases, and GRCh38 Alu overlap. Alu intervals use 0-based half-open coordinates.

The model preserves native missing-value branches. Summed tree outputs are converted through the logistic function to `ESP_SCORE`. This is a classifier score, not a calibrated probability.

## Input and output

Site lists contain headerless `chrom<TAB>pos` rows with positive 1-based positions. Empty lines are ignored. The FASTA requires `.fai`; the reference guard checks hg38 `chr1` length 248956422 when present, allowing small references below 10000000 bases.

Nearby sites share BAM and FASTA access. Results are restored to input order. Output rows are `chrom<TAB>pos<TAB>score`. Identical inputs and model files produce identical scores across worker counts.

## VCF decision

`PASS` means `ESP_SCORE` meets the bundle threshold; otherwise the filter is `LOW_SCORE`. `VAF`, `DEPTH`, `AC`, and `STRAND` describe candidate evidence. No second model or additional site filter modifies the classifier decision.
