//! esperanto-map — splice-aware genomic aligner (Track 1 legacy path + Track 2 junction-kmer).
//!
//! Spec: docs/specs/crates/map.md (sole source of truth). Modules are added
//! by the implementation wave; keep this file's `pub mod` list in sync.

pub mod align;
pub mod baln;
pub mod bam;
pub mod intron_chain;
pub mod jkmer;
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
pub mod stats;
pub mod span;
pub mod splice;
pub mod split;

pub use mapq::ReadAlignment;
