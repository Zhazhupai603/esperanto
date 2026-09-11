//! Transcript-strand calling (DESIGN §1.2b): each candidate site outputs STRAND = +/-/amb + evidence code.
//! Evidence priority: 1) LIB=stranded user parameter → 2) junction orientation (GT-AG) →
//! 3) --gtf gene annotation (genes on both strands → amb) → 4) none of the above → amb.

use crate::LibType;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrandCall {
    Plus,
    Minus,
    Amb,
}

impl StrandCall {
    pub fn as_str(self) -> &'static str {
        match self {
            StrandCall::Plus => "+",
            StrandCall::Minus => "-",
            StrandCall::Amb => "amb",
        }
    }
}

/// Infer the transcript strand and primary evidence code for one candidate site.
///
/// - `fwd_depth`/`rev_depth`: plus/minus-strand read counts at the site (direction maps directly when LIB=stranded)
/// - `junc_plus`/`junc_minus`: intron-orientation votes from split reads covering the site
/// - `gtf`: GTF hit (has plus-strand gene, has minus-strand gene); None when --gtf is not provided
pub fn infer_strand(
    lib: LibType,
    fwd_depth: u64,
    rev_depth: u64,
    junc_plus: u32,
    junc_minus: u32,
    gtf: Option<(bool, bool)>,
) -> (StrandCall, &'static str) {
    // 1) User-declared strand-specific library (dUTP-type): read direction maps directly
    if lib == LibType::Stranded {
        return match fwd_depth.cmp(&rev_depth) {
            std::cmp::Ordering::Greater => (StrandCall::Plus, "LIB"),
            std::cmp::Ordering::Less => (StrandCall::Minus, "LIB"),
            std::cmp::Ordering::Equal => (StrandCall::Amb, "LIB"),
        };
    }
    // 2) Junction orientation (GT-AG dinucleotides): only unanimous votes decide; conflict → amb
    if junc_plus > 0 && junc_minus == 0 {
        return (StrandCall::Plus, "JUNC");
    }
    if junc_minus > 0 && junc_plus == 0 {
        return (StrandCall::Minus, "JUNC");
    }
    if junc_plus > 0 && junc_minus > 0 {
        return (StrandCall::Amb, "JUNC");
    }
    // 3) Gene annotation: site in a known gene → that gene's strand; genes on both strands → amb
    match gtf {
        Some((true, false)) => (StrandCall::Plus, "GTF"),
        Some((false, true)) => (StrandCall::Minus, "GTF"),
        Some((true, true)) => (StrandCall::Amb, "GTF"),
        // 4) none of the above → amb
        _ => (StrandCall::Amb, "NONE"),
    }
}
