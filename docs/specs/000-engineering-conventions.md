# ESPERANTO Engineering Conventions (mandatory for all crates)

## Purpose

ESPERANTO is a command-line tool designed for RNA editing analysis (pure Rust, single CLI). This document defines the unified engineering conventions for all crates.

## Hard rules

1. **Specs are the source of truth**: do not record porting provenance in code comments or repo docs. For details not covered by the specs, list them explicitly rather than guessing.
2. **No tests in the repo**: no `#[cfg(test)]` modules, fixtures, or CI configuration. Verification is handled by external differential tests.
3. **Determinism**: identical input must produce byte-identical output. Forbidden: emitting unsorted HashMap iteration order directly, wall-clock entering artifacts (unless explicitly passed in as a parameter), floating-point accumulation order depending on parallel scheduling.
4. **Error handling**: library crates define error types with `thiserror`; do not use `unwrap`/`expect` for anticipated runtime failures; `panic!` is limited to internal invariants.
5. **unsafe**: `#![deny(unsafe_code)]` by default; any required unsafe must be justified in code.
6. **Dependencies**: use only dependencies declared in the workspace root `Cargo.toml` `[workspace.dependencies]`; a new dependency requires a deliberate, documented addition.
7. **Public API comments in English**; internal implementation comments may use either language; user-facing strings (CLI help, error messages) in English.

## Acceptance gates

1. Two runs on the same input produce byte-identical output, and `cargo clippy -p <crate> -- -D warnings` yields zero warnings;
2. Where model/scoring semantics are involved: gold-standard metrics (AUROC/AUPRC) must not fall below the anchors (A 0.9982 / B 0.9961);
3. New behavior paths: positive/negative constructed controls (external, not in the repo).

## Workspace layout

```
crates/<name>/          short directory names; package name esperanto-<name>
  Cargo.toml
  src/lib.rs
```

- edition = "2021", resolver = "2"
- release profile: lto = "thin", codegen-units = 1 (already configured in the root Cargo.toml)
- Each crate has one clear responsibility; crates communicate only through public type contracts; copying implementations across crates is forbidden

## Deliverable definition (per crate)

1. Compiles: `cargo build -p esperanto-<name>` succeeds, and `cargo clippy -p esperanto-<name> -- -D warnings` produces no warnings
2. Doc comment at the top of `src/lib.rs`: crate responsibility, public API contract, invariants
