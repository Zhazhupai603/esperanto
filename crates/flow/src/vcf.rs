//! Minimal VCF output (spec §output VCF): merges candidates.bed columns with
//! score probabilities into `<out>/sites.vcf`.

use std::fmt::Write as _;
use std::fs;
use std::path::Path;

use crate::params::{Entry, RunParams};
use crate::FlowError;

/// Passthrough columns from one candidates.bed row (original text preserved
/// byte-for-byte; no re-formatting).
struct BedInfo {
    strand: String,
    evid: String,
    depth: String,
    vaf: String,
}

fn parse_bed_info(path: &Path) -> Result<Vec<BedInfo>, FlowError> {
    let text = fs::read_to_string(path)?;
    let mut rows = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 8 {
            return Err(FlowError::BedParse {
                line: i + 1,
                msg: format!("expected >= 8 tab-separated columns, got {}", cols.len()),
            });
        }
        rows.push(BedInfo {
            strand: cols[3].to_string(),
            evid: cols[4].to_string(),
            depth: cols[6].to_string(),
            vaf: cols[7].to_string(),
        });
    }
    Ok(rows)
}

/// Reference base at 1-based `pos`, uppercased; `N` when unavailable (spec:
/// not found = N).
fn ref_base(fai: &rust_htslib::faidx::Reader, chrom: &str, pos: i64) -> char {
    if pos < 1 {
        return 'N';
    }
    match fai.fetch_seq(chrom, (pos - 1) as usize, (pos - 1) as usize) {
        Ok(seq) if !seq.is_empty() => (seq[0] as char).to_ascii_uppercase(),
        _ => 'N',
    }
}

/// Write `<out>/sites.vcf` (VCF v4.2 minimal contract).
pub fn write_vcf(
    params: &RunParams,
    entry: Entry,
    sites: &[(String, i64)],
    probs: &[Option<f64>],
) -> Result<(), FlowError> {
    let bed_rows: Option<Vec<BedInfo>> = match entry {
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
                    "candidates.bed row count {} != sites count {}",
                    rows.len(),
                    sites.len()
                ),
            });
        }
    }

    let fai =
        rust_htslib::faidx::Reader::from_path(&params.fasta).map_err(|e| FlowError::Stage {
            stage: "vcf",
            source: Box::new(e),
        })?;
    let refname = params
        .fasta
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| params.fasta.display().to_string());

    let mut out = String::new();
    out.push_str("##fileformat=VCFv4.2\n");
    let _ = writeln!(out, "##reference={refname}");
    out.push_str(
        "##INFO=<ID=RE_PROB,Number=1,Type=Float,Description=\"RNA editing probability (score)\">\n",
    );
    out.push_str(
        "##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Variant allele frequency (scan)\">\n",
    );
    out.push_str("##INFO=<ID=DEPTH,Number=1,Type=Integer,Description=\"Site depth (scan)\">\n");
    out.push_str(
        "##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Strand assignment (scan)\">\n",
    );
    out.push_str("##INFO=<ID=EVID,Number=1,Type=String,Description=\"Evidence type (scan)\">\n");
    out.push_str("##FILTER=<ID=PASS,Description=\"RE_PROB >= 0.5\">\n");
    out.push_str("##FILTER=<ID=LOW_SCORE,Description=\"RE_PROB < 0.5\">\n");
    out.push_str("##FILTER=<ID=UNSCORED,Description=\"site on a contig with no matching model (hybrid reference)\">\n");
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

    for (i, ((chrom, pos), prob)) in sites.iter().zip(probs).enumerate() {
        let (filter, mut info) = match prob {
            Some(p) => (
                if *p >= 0.5 { "PASS" } else { "LOW_SCORE" },
                format!("RE_PROB={p}"),
            ),
            None => ("UNSCORED", String::new()),
        };
        if let Some(rows) = &bed_rows {
            let r = &rows[i];
            let sep = if info.is_empty() { "" } else { ";" };
            let _ = write!(
                info,
                "{sep}VAF={};DEPTH={};STRAND={};EVID={}",
                r.vaf, r.depth, r.strand, r.evid
            );
        }
        let refb = ref_base(&fai, chrom, *pos);
        let _ = writeln!(out, "{chrom}\t{pos}\t.\t{refb}\t<RE>\t.\t{filter}\t{info}");
    }

    fs::write(params.out_dir.join("sites.vcf"), out)?;
    Ok(())
}
