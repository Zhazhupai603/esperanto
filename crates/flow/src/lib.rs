//! esperanto-flow — pipeline orchestration (entries FASTQ PE/SE / BAM /
//! BAM+sites, species guardrail).
//!
//! Spec: docs/specs/crates/flow.md (single source of truth).
//! Single entry [`run_pipeline`]: wires qc → map → sort → scan → score → vcf
//! by entry type, artifacts landing in `<out_dir>/<stage>/`. Orchestration
//! carries zero scientific semantics: every numeric contract belongs to the
//! stage crates; flow only wires, fails early, and writes deterministically.

pub mod error;
pub mod filter;
pub mod guard;
pub mod params;
pub mod stages;
pub mod vcf;

pub use error::FlowError;
pub use params::{DeviceAsk, Entry, RunParams};
pub use stages::map_stage;

/// Run the pipeline selected by `params.entry()`.
pub fn run_pipeline(params: &RunParams) -> Result<(), FlowError> {
    stages::run(params)
}
