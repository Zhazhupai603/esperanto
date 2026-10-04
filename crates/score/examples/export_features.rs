//! Export the frozen 27-feature matrix for a set of sites.
//!
//! Usage:
//!   export_features --bam in.bam --fasta ref.fa --alu alu.bed \
//!       --sites sites.tsv --out features.tsv --threads 8
//!
//! `--sites` is a headerless TSV with columns `chrom pos` (1-based), extra
//! columns are ignored. Output rows follow the input site order:
//! `chrom pos f1..f27`, with per-site errors reported as `ERROR <message>`.

use anyhow::{bail, Context, Result};
use esperanto_score::unified::{feature_matrix, AluIndex, FEATURE_COUNT};
use std::io::Write as _;
use std::path::PathBuf;

fn main() -> Result<()> {
    let mut bam = None;
    let mut fasta = None;
    let mut alu = None;
    let mut sites_path = None;
    let mut out = None;
    let mut threads = 8usize;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        let mut take = |name: &str| -> Result<String> {
            args.next().with_context(|| format!("missing value for {name}"))
        };
        match arg.as_str() {
            "--bam" => bam = Some(PathBuf::from(take("--bam")?)),
            "--fasta" => fasta = Some(PathBuf::from(take("--fasta")?)),
            "--alu" => alu = Some(PathBuf::from(take("--alu")?)),
            "--sites" => sites_path = Some(PathBuf::from(take("--sites")?)),
            "--out" => out = Some(PathBuf::from(take("--out")?)),
            "--threads" => threads = take("--threads")?.parse()?,
            other => bail!("unknown argument {other}"),
        }
    }
    let bam = bam.context("--bam required")?;
    let fasta = fasta.context("--fasta required")?;
    let alu = alu.context("--alu required")?;
    let sites_path = sites_path.context("--sites required")?;
    let out = out.context("--out required")?;

    let text = std::fs::read_to_string(&sites_path)?;
    let mut sites = Vec::new();
    for (line_no, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let chrom = fields.next().context("missing chrom")?.to_owned();
        let pos: i64 = fields
            .next()
            .context("missing pos")?
            .parse()
            .with_context(|| format!("invalid pos at line {}", line_no + 1))?;
        sites.push((chrom, pos));
    }
    eprintln!("loaded {} sites", sites.len());

    let alu_index = AluIndex::from_bed(&alu)?;
    let matrix = feature_matrix(&bam, &fasta, &alu_index, &sites, threads)?;

    let mut writer = std::io::BufWriter::new(std::fs::File::create(&out)?);
    let mut errors = 0usize;
    for ((chrom, pos), row) in sites.iter().zip(matrix.iter()) {
        match row {
            Ok(values) => {
                write!(writer, "{chrom}\t{pos}")?;
                debug_assert_eq!(values.len(), FEATURE_COUNT);
                for value in values {
                    write!(writer, "\t{value}")?;
                }
                writeln!(writer)?;
            }
            Err(message) => {
                errors += 1;
                writeln!(writer, "{chrom}\t{pos}\tERROR\t{message}")?;
            }
        }
    }
    eprintln!("wrote {} rows ({} errors) to {}", sites.len(), errors, out.display());
    Ok(())
}
