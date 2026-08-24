//! esperanto-score — native Caduceus encoder + 5-fold ensemble head (feature_spec.json v1).
pub mod bf16gemm;
pub mod bundle;
pub mod caduceus;
pub mod encoder;
pub mod head;
pub mod mamba;
pub mod pipeline;
pub mod simdexp;

pub use bundle::{load_bundle, Bundle, EmbCache, NormStats, ScoreError};
pub use head::{re_prob_ensemble, re_prob_fold};
