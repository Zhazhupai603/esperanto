//! `esperanto pile` — pileup feature extraction (→
//! `pile::extract_pileup_features_batch`).

use std::io::Write as _;
use std::path::PathBuf;

use anyhow::bail;
use clap::Args;

#[derive(Args)]
pub struct PileArgs {
    /// Input BAM (requires .bai/.csi).
    #[arg(long)]
    bam: PathBuf,
    /// Single-site mode: contig name (with --pos).
    #[arg(long)]
    chrom: Option<String>,
    /// Single-site mode: 1-based position (with --chrom).
    #[arg(long)]
    pos: Option<i64>,
    /// Batch mode: sites file `chrom<TAB>pos` (1-based).
    #[arg(long)]
    sites: Option<PathBuf>,
    /// Output TSV (default stdout).
    #[arg(long)]
    out: Option<PathBuf>,
}

pub fn run(a: PileArgs) -> anyhow::Result<()> {
    let sites: Vec<(String, i64)> = match (&a.chrom, a.pos, &a.sites) {
        (Some(c), Some(p), None) => vec![(c.clone(), p)],
        (None, None, Some(path)) => {
            let text = std::fs::read_to_string(path)?;
            let mut v = Vec::new();
            for (i, line) in text.lines().enumerate() {
                if line.is_empty() {
                    continue;
                }
                let (chrom, pos) = line
                    .split_once('\t')
                    .ok_or_else(|| anyhow::anyhow!("bad site line {}: {}", i + 1, line))?;
                let pos: i64 = pos
                    .parse()
                    .map_err(|_| anyhow::anyhow!("bad pos at line {}: {}", i + 1, line))?;
                v.push((chrom.to_string(), pos));
            }
            v
        }
        _ => bail!("give exactly one of --chrom+--pos or --sites"),
    };

    let mut bam = rust_htslib::bam::IndexedReader::from_path(&a.bam)?;
    let refs: Vec<(&str, i64)> = sites.iter().map(|(c, p)| (c.as_str(), *p)).collect();
    let feats = esperanto_pile::extract_pileup_features_batch(&mut bam, &refs)?;

    let mut buf = String::new();
    use std::fmt::Write as _;
    for ((chrom, pos), f) in sites.iter().zip(&feats) {
        let _ = writeln!(
            buf,
            "{chrom}\t{pos}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7]
        );
    }
    match &a.out {
        Some(path) => std::fs::write(path, buf)?,
        None => {
            let stdout = std::io::stdout();
            let mut lock = stdout.lock();
            lock.write_all(buf.as_bytes())?;
        }
    }
    Ok(())
}
