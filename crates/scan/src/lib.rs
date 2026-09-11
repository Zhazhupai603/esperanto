//! esperanto-scan — strand-resolved candidate editing-site caller (scatter engine).
//!
//! Role (spec: docs/specs/crates/scan.md): genome-wide scan from BAM/.baln, emitting
//! candidates.bed (10-column contract). No hard filters: strand-resolved statistics, the soft score
//! CALL_SCORE (0–1), and the evidence codes EVID all pass through, consumed downstream with thresholds.
//! The threshold parameter only marks (evid gains `,MS`); by default all candidates are emitted.
//!
//! Chromosome-level parallelism (rayon); output sorted by (chrom, pos), thread-count independent, byte-identical.

pub mod annot;
pub mod baln;
pub mod call;
pub mod count;
pub mod error;
pub mod out;
pub mod scatter;
pub mod score;
pub mod strand;

pub use error::CallError;
pub use score::CallSpec;
pub use strand::StrandCall;

use std::path::PathBuf;

/// Library strandedness: Stranded (dUTP-type; read direction maps directly to the transcribed strand) / Unstranded.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LibType {
    Unstranded,
    Stranded,
}

#[derive(Debug, Clone)]
pub struct CallParams {
    /// Input BAM (requires .csi/.bai index).
    pub bam: PathBuf,
    /// Output candidates.bed.
    pub out: PathBuf,
    /// Reference FASTA (requires .fai): ref base / homopolymer / junction orientation; without it, falls back to majority pseudo-ref.
    pub fasta: Option<PathBuf>,
    /// Gene annotation GTF (strand-call evidence priority 3).
    pub gtf: Option<PathBuf>,
    /// gnomAD VCF single file or per-chrom directory (soft down-weighting, never removes); configured but unusable = hard error.
    pub gnomad: Option<PathBuf>,
    pub lib: LibType,
    /// C>U symmetric-mode switch (off by default).
    pub enable_cu: bool,
    /// Threshold marks only, never removes; None = emit all candidates.
    pub min_call_score: Option<f64>,
    /// Override scoring spec JSON; None = crate-bundled call_spec.v2.json.
    pub spec: Option<PathBuf>,
    /// Thread count; 0 = all cores.
    pub threads: usize,
    /// Binary input: map-output .baln (skips BAM decompression/indexing). When Some, the bam field is not opened;
    /// the contig list and read set both come from the .baln coordinate index (single-pass scan); output is
    /// byte-identical to the BAM path (both sources share the same scatter kernel scatter_one_record).
    pub baln: Option<PathBuf>,
}

#[derive(Debug, Clone, Copy)]
pub struct CallStats {
    pub candidates: usize,
    pub contigs: usize,
    /// Number of sites hit by gnomAD soft down-weighting (0 when gnomad is not enabled).
    pub gnomad_hits: usize,
}

/// Sole engine entry point: the scatter caller.
pub fn run_call(params: &CallParams) -> Result<CallStats, CallError> {
    call::run_call(params)
}
