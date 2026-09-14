//! CMC extraction parity helper: reads `chrom<TAB>pos` sites and an indexed
//! BAM, writes the 16x50 f32 tensors (little-endian, row-major) to stdout file.
//!
//! Usage: `cargo run --release -p esperanto-score --example cmc_extract_parity
//! -- <bam> <sites.tsv> <out.bin>`

use std::io::Write as _;

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    anyhow::ensure!(
        args.len() == 4,
        "usage: cmc_extract_parity <bam> <sites.tsv> <out.bin>"
    );
    let text = std::fs::read_to_string(&args[2])?;
    let mut sites: Vec<(String, i64)> = Vec::new();
    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        let (c, p) = line
            .split_once('\t')
            .ok_or_else(|| anyhow::anyhow!("bad line: {line}"))?;
        sites.push((c.to_string(), p.parse::<i64>()?));
    }
    let refs: Vec<(&str, i64)> = sites.iter().map(|(c, p)| (c.as_str(), *p)).collect();
    let mut bam = rust_htslib::bam::IndexedReader::from_path(&args[1])?;
    let feats = esperanto_score::cmc::extract_cmc_batch(&mut bam, &refs)
        .map_err(|e| anyhow::anyhow!("{e}"))?;
    let mut out = std::fs::File::create(&args[3])?;
    for f in &feats {
        for v in f.iter() {
            out.write_all(&v.to_le_bytes())?;
        }
    }
    eprintln!("{} sites -> {}", feats.len(), args[3]);
    Ok(())
}
