//! VCF output: merges candidates.bed with unified classifier scores into
//! `<out>/sites.vcf`.

use std::fmt::Write as _;
use std::fs;
use std::path::Path;

use crate::params::{Entry, RunParams};
use crate::FlowError;

/// Passthrough columns from one candidates.bed row (original text preserved
/// byte-for-byte; no re-formatting).
struct BedInfo {
    ref_base: char,
    depth: String,
    counts: [String; 4],
    edit_freq: String,
}

/// Editing-consistent alternate for a reference base; `None` for a base
/// that cannot be A-to-I edited.
fn edit_alt(ref_base: char) -> Option<char> {
    match ref_base {
        'A' => Some('G'),
        'T' => Some('C'),
        _ => None,
    }
}

/// Transcript strand implied by the reference base of an A-to-I site: an A
/// on the forward view is edited on the plus strand, a T on the minus strand.
fn strand_of(ref_base: char) -> char {
    match ref_base {
        'A' => '+',
        'T' => '-',
        _ => '.',
    }
}

fn parse_bed_info(path: &Path) -> Result<Vec<BedInfo>, FlowError> {
    let text = fs::read_to_string(path)?;
    let mut rows = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 10 {
            return Err(FlowError::BedParse {
                line: i + 1,
                msg: format!("expected 10 tab-separated columns, got {}", cols.len()),
            });
        }
        let ref_base = cols[3].chars().next().unwrap_or('N').to_ascii_uppercase();
        rows.push(BedInfo {
            ref_base,
            depth: cols[4].to_string(),
            counts: [
                cols[5].to_string(),
                cols[6].to_string(),
                cols[7].to_string(),
                cols[8].to_string(),
            ],
            edit_freq: cols[9].to_string(),
        });
    }
    Ok(rows)
}

/// Reference base at 1-based `pos`, uppercased; `N` when unavailable.
fn ref_base(fai: &rust_htslib::faidx::Reader, chrom: &str, pos: i64) -> char {
    if pos < 1 {
        return 'N';
    }
    match fai.fetch_seq(chrom, (pos - 1) as usize, (pos - 1) as usize) {
        Ok(seq) if !seq.is_empty() => (seq[0] as char).to_ascii_uppercase(),
        _ => 'N',
    }
}

/// Write unified scores and observed read evidence to VCF v4.2.
pub fn write_vcf(
    params: &RunParams,
    entry: Entry,
    sites: &[(String, i64)],
    scores: &[f64],
    threshold: f64,
) -> Result<(), FlowError> {
    if scores.len() != sites.len() {
        return Err(FlowError::BedParse {
            line: 0,
            msg: format!("score count {} != site count {}", scores.len(), sites.len()),
        });
    }
    let bed_rows = match entry {
        Entry::BamSites => None,
        _ => Some(parse_bed_info(
            &params.out_dir.join("scan").join("candidates.bed"),
        )?),
    };
    if let Some(rows) = &bed_rows {
        if rows.len() != sites.len() {
            return Err(FlowError::BedParse {
                line: 0,
                msg: format!(
                    "candidates.bed row count {} != site count {}",
                    rows.len(),
                    sites.len()
                ),
            });
        }
    }
    let fai = rust_htslib::faidx::Reader::from_path(&params.fasta).map_err(|source| {
        FlowError::Stage {
            stage: "vcf",
            source: Box::new(source),
        }
    })?;
    let refname = params
        .fasta
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| params.fasta.display().to_string());
    let mut out = String::from("##fileformat=VCFv4.2\n");
    let _ = writeln!(out, "##reference={refname}");
    out.push_str("##INFO=<ID=ESP_SCORE,Number=1,Type=Float,Description=\"Unified molecular evidence classifier score\">\n");
    out.push_str("##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Editing-consistent alternate allele frequency\">\n");
    out.push_str("##INFO=<ID=DEPTH,Number=1,Type=Integer,Description=\"Counted depth (MAPQ >= 20, base quality >= 13)\">\n");
    out.push_str(
        "##INFO=<ID=AC,Number=4,Type=Integer,Description=\"Counted reads by base: A,C,G,T\">\n",
    );
    out.push_str("##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Transcript strand implied by the reference base\">\n");
    let _ = writeln!(
        out,
        "##FILTER=<ID=PASS,Description=\"ESP_SCORE >= {threshold}\">"
    );
    let _ = writeln!(
        out,
        "##FILTER=<ID=LOW_SCORE,Description=\"ESP_SCORE < {threshold}\">"
    );
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for (i, ((chrom, pos), score)) in sites.iter().zip(scores).enumerate() {
        let filter = if *score >= threshold {
            "PASS"
        } else {
            "LOW_SCORE"
        };
        let reference = match &bed_rows {
            Some(rows) => rows[i].ref_base,
            None => ref_base(&fai, chrom, *pos),
        };
        let alternate = edit_alt(reference).unwrap_or('N');
        let mut info = format!("ESP_SCORE={score}");
        if let Some(rows) = &bed_rows {
            let row = &rows[i];
            let _ = write!(
                info,
                ";VAF={};DEPTH={};AC={},{},{},{}",
                row.edit_freq,
                row.depth,
                row.counts[0],
                row.counts[1],
                row.counts[2],
                row.counts[3]
            );
        }
        let _ = write!(info, ";STRAND={}", strand_of(reference));
        let _ = writeln!(
            out,
            "{chrom}\t{pos}\t.\t{reference}\t{alternate}\t.\t{filter}\t{info}"
        );
    }
    fs::write(params.out_dir.join("sites.vcf"), out)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn unified_vcf_uses_only_current_scores_and_observations() {
        let dir =
            std::env::temp_dir().join(format!("esperanto-unified-vcf-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(dir.join("scan")).unwrap();
        fs::create_dir_all(dir.join("hyper")).unwrap();
        let fasta = dir.join("ref.fa");
        fs::write(&fasta, ">chr1\nAT\n").unwrap();
        fs::write(dir.join("ref.fa.fai"), "chr1\t2\t6\t2\t3\n").unwrap();
        fs::write(dir.join("scan/candidates.bed"),
            "chr1\t0\t1\tA\t10\t8\t0\t2\t0\t0.20000000\t\nchr1\t1\t2\tT\t10\t0\t2\t0\t8\t0.20000000\t\n").unwrap();
        fs::write(
            dir.join("scan/candidates.verified.bed"),
            "invalid retired data",
        )
        .unwrap();
        fs::write(dir.join("hyper/sites.tsv"), "invalid retired data").unwrap();
        fs::write(dir.join("species.json"), "invalid retired data").unwrap();
        let params = RunParams {
            r1: Vec::new(),
            r2: Vec::new(),
            bam: Some(PathBuf::from("input.bam")),
            sites: None,
            index: None,
            fasta,
            gtf: None,
            bundle: PathBuf::from("bundle"),
            l1_bundle: None,
            out_dir: dir.clone(),
            threads: 1,
            stop_after: None,
        };
        let sites = vec![("chr1".into(), 1), ("chr1".into(), 2)];
        write_vcf(&params, Entry::Bam, &sites, &[0.95, 0.8], 0.9).unwrap();
        let text = fs::read_to_string(dir.join("sites.vcf")).unwrap();
        assert!(text.contains("base quality >= 13"));
        assert!(text.contains(
            "chr1\t1\t.\tA\tG\t.\tPASS\tESP_SCORE=0.95;VAF=0.20000000;DEPTH=10;AC=8,0,2,0;STRAND=+"
        ));
        assert!(text.contains("chr1\t2\t.\tT\tC\t.\tLOW_SCORE\tESP_SCORE=0.8"));
        assert!(!text.contains("REALIGN"));
        assert!(!text.contains("HYPER"));
        assert!(!text.contains("UNSCORED"));
        assert!(write_vcf(&params, Entry::Bam, &sites, &[0.95], 0.9).is_err());
        fs::remove_dir_all(dir).unwrap();
    }
}
