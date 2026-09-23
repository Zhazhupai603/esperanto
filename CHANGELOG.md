# Changelog

## [1.0.8] - 2026-09-23

### Fixed

- Report generation is now deterministic for multi-copy gene annotation: a gene name present on more than one contig (e.g. PAR genes on chrX/chrY) is attributed to the last GTF row in file order, as documented, instead of whichever copy happened to win a hash-map iteration. Repeated report generation over the same run directory is byte-identical; site tables and scores are unaffected.
- Rescue resume is crash-safe: the rescue completion key is written only after the merged BAM is durable, and `resume` verifies that the sealed rescued count, the rescued sidecar and the `RE:Z:collapsed` records in the BAM all agree before accepting a rescue stage as done. A run interrupted around the merge re-runs the rescue to a consistent state instead of continuing over a half-merged BAM. Resume semantics change: a previously "complete" rescue whose artifacts disagree is now rebuilt (rescued reads are recovered from the unmapped pool or from their own earlier `RE:Z:collapsed` records, with duplicate pool entries dropped), and a stale `RE:Z:collapsed` record whose read no longer verifies is removed at the next merge.

## [1.0.7] - 2026-09-22

### Fixed

- Hyperedited site aggregation now normalizes each rescued read to its genomic polarity before grouping: mate pairs rescued by opposite passes over the same event merge into one cluster and one site row (previously split into AG/TC groups with halved support and half the VCF rows written with the wrong REF/ALT polarity). A tripwire cross-checks every kept site against the reference, reports the outcome in `hyper/polarity_report.json`, and prints a warning listing any inconsistent coordinates (the run still completes); a reference that cannot be read is a hard error.
- Duplicate hyper sites are rejected instead of silently merged or overwritten: the hyper stage refuses to write a sites table with a duplicated (chrom, pos) key, and readers of `hyper/sites.tsv` fail on duplicates, unknown directions, or malformed fields.
- Rescued records written to the merged BAM now carry offset-free Phred quality values; the second `Record::set` overwrote the converted values with raw ASCII bytes, inflating every rescued read's quality by +33.
- `scan` writes `"snr": null` instead of invalid `inf` JSON when a control channel is empty (crashed the downstream verify/score parse).
- Strict input handling: a truncated or malformed `map/unmapped.fq.gz` record and a malformed `map/rescued.tsv` row are hard errors naming the file, the 1-based line and the reason (previously such content was silently dropped); `hyper/sites.tsv` rows with a wrong column count are rejected; a present-but-unreadable `sites.tsv`, `candidates.verified.bed`, `hyper_scores.tsv` or `species.json` is an explicit error instead of being treated as absent.

### Changed

- Removed the `--snp` option from `run` and `scan` (known-SNP exclusion remains an internal channel; the refs-discovery dbSNP default silently altered calling protocols). `resume` now also ignores a known-SNP setting frozen in an older run.json and prints a warning, so no entry point can re-enable the channel; scan-stage artifacts that may have been produced under the filter are rebuilt from scan onward.

## [1.0.6] - 2026-09-14

### Changed

- Updated site scoring and added scores for rescued editing sites.
- Scoring now runs on the CPU; `--device` and `--caduceus` are no longer supported.

## [1.0.5] - 2026-09-11

### Changed

- Improved alignment, candidate calling and hyperediting detection.

## [1.0.4] - 2026-08-31

### Changed

- Improved read rescue, site scoring and report navigation.

## [1.0.3] - 2026-08-29

### Added

- Added hybrid-reference support for human knock-in mouse samples.

## [1.0.2] - 2026-08-27

### Added

- Added resumable runs, per-sample output directories, and update and delete commands.

## [1.0.1] - 2026-08-27

### Added

- Added reference setup, HTML reports and hyperedited-read recovery.

## [1.0.0] - 2026-08-24

### Added

- Initial release with read processing, alignment, editing-site detection and scoring.
