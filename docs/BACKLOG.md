# Enhancement Backlog

Non-blocking enhancements tracked for future releases. Mainline correctness and evaluation gates take priority; every item below must be validated against reference data before landing.

## Stage selection

Stage selection flags (`--only`, `--stop-after`) on top of the existing `resume` subcommand.

## Track-2 (jkmer) wiring in the pipeline

The mapper's Track-2 junction k-mer index is currently not wired into `flow`; the transcriptome-first engine and 2-pass alignment are the primary paths. Wiring the jkmer index build and pipeline hook is deferred.

## Automatic sorting of user-provided BAMs

BAM entry points require a coordinate-sorted, indexed input and refuse to run otherwise. Sorting user BAMs automatically (via an explicit `--sort-input` flag) is a candidate enhancement — `bamio` already provides the sort primitive.

## Layered config / presets

A YAML config layer with preset profiles (standard / strict / exploratory) and a Default → config → preset → CLI precedence chain. Requires adding a YAML dependency to the workspace.

## `build-ref` / `fetch-refs` subcommands

A `build-ref` subcommand to produce the full reference bundle (paidx + transcript index + transcript map), and a `fetch-refs` subcommand to download and install a reference set. The `index` subcommand (paidx build) is already available.

## `ui` / `web` / `bfq2fq` subcommands

Optional TUI dashboard, web service, and a `.bfq`→FASTQ conversion utility.

## `--keep-*` intermediate cleanup flags

Per-stage flags to prune intermediate artifacts after a run. The current release keeps all stage outputs for reproducibility.
