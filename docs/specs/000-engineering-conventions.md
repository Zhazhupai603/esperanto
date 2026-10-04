# Engineering conventions

ESPERANTO is a Rust command-line application for RNA editing analysis. Each crate has a focused responsibility, and the CLI delegates scientific computations to library crates.

## Public contracts

- Crate specifications document supported inputs, outputs, coordinate conventions, parameters, and errors.
- Public APIs, documentation, command help, and error messages use English.
- Errors identify the affected file or input value and propagate to the caller. Panics are reserved for violated internal invariants.
- Dependencies are declared in the workspace. Crates communicate through their public APIs.
- Unsafe operations require a documented memory-safety justification.

## Determinism

- Candidate tables, classifier scores, and VCF records have stable ordering and deterministic values for identical inputs.
- Hash-map iteration order must not determine serialized record order or floating-point accumulation order.
- Timing and input modification times belong to run metadata and QC reports, not scientific result fields.
- Index formats document timestamps explicitly where present.

## Validation

Changes must compile without warnings and preserve the applicable input/output contracts. Verification should cover the affected scientific behavior, error handling, and reproducibility across worker counts. Model feature extraction must remain consistent with the released model bundle.

## Workspace layout

Each library resides in `crates/<name>/` and is named `esperanto-<name>`. The executable is `esperanto`, built by `esperanto-cli`. The workspace uses Rust edition 2021 and Cargo resolver 2.
