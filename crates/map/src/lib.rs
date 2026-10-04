//! esperanto-map — splice-aware genomic aligner.
//!
//! Spec: docs/specs/crates/map.md.

pub mod align;
pub mod bam;
pub mod intron_chain;
pub mod chain;
pub mod error;
pub mod evidence;
pub mod extend;
pub mod fasta;
pub mod fastq;
pub mod gtf;
pub mod index;
pub mod index_io;
pub mod mapq;
pub mod myers_ea;
pub mod pair;
pub mod pipeline;
pub mod seed;
pub mod sinks;
pub mod stats;
pub mod span;
pub mod splice;
pub mod split;

pub use mapq::ReadAlignment;
