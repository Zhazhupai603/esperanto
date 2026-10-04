//! Native reference-aware molecular evidence scoring for ESPERANTO 1.1.
pub mod bundle;
pub mod pipeline;
pub mod unified;

pub use bundle::{is_unified_bundle, load_primary_threshold, ScoreError};
