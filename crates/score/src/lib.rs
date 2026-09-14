//! esperanto-score — CMC feature extraction + frozen CNN/Transformer scoring
//! (model bundle v1.7 contract).
pub mod bundle;
pub mod cmc;
pub mod cmcnet;
pub mod pipeline;

pub use bundle::{is_cmc_bundle, load_cmc_bundle, load_cmc_thresholds, CmcBundle, ScoreError};
