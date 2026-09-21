# Changelog

## [Unreleased] - fix/c2-rescue-reporting

### Fixed

- Hyperedited site aggregation now normalizes each rescued read to its genomic polarity before grouping: mate pairs rescued by opposite passes over the same event merge into one cluster and one site row (previously split into AG/TC groups with halved support and half the VCF rows written with the wrong REF/ALT polarity). A tripwire counts reference-base mismatches at kept sites.
- Rescued records written to the merged BAM now carry offset-free Phred quality values; the second `Record::set` overwrote the converted values with raw ASCII bytes, inflating every rescued read's quality by +33.
- `scan` writes `"snr": null` instead of invalid `inf` JSON when a control channel is empty (crashed the downstream verify/score parse).

### Changed

- Removed the `--snp` option from `run` and `scan` (known-SNP exclusion remains an internal channel; the refs-discovery dbSNP default silently altered calling protocols).

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
