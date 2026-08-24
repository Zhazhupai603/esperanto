//! candidates.bed output (10-column contract, DESIGN §1.3):
//! chrom  pos  pos+1  strand(+/-/amb)  evid  call_score  depth  var_freq  fwd_freq  rev_freq
//!
//! pos is 0-based (BED half-open interval [pos, pos+1)). Sorted by (chrom, pos) before writing:
//! deterministic output, thread-count independent.

use crate::error::CallError;
use crate::strand::StrandCall;
use std::io::Write as _;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Candidate {
    pub chrom: String,
    pub pos0: i64,
    pub strand: StrandCall,
    pub evid: String,
    pub score: f64,
    pub depth: u64,
    pub var_freq: f64,
    pub fwd_freq: f64,
    pub rev_freq: f64,
}

pub fn write_bed(path: &Path, cands: &mut [Candidate]) -> Result<(), CallError> {
    cands.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.pos0.cmp(&b.pos0)));
    let mut w =
        std::io::BufWriter::new(std::fs::File::create(path).map_err(|e| CallError::Io {
            path: path.display().to_string(),
            source: e,
        })?);
    for c in cands.iter() {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{:.4}\t{:.4}\t{:.4}",
            c.chrom,
            c.pos0,
            c.pos0 + 1,
            c.strand.as_str(),
            c.evid,
            c.score,
            c.depth,
            c.var_freq,
            c.fwd_freq,
            c.rev_freq
        )
        .map_err(|e| CallError::Io {
            path: path.display().to_string(),
            source: e,
        })?;
    }
    w.flush().map_err(|e| CallError::Io {
        path: path.display().to_string(),
        source: e,
    })
}
