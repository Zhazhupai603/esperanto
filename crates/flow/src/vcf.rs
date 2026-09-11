//! VCF output: merges candidates.bed with score probabilities into
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
    /// Standard evidence-filter tags (column 11; empty when the scan stage
    /// produced the bed before this column existed).
    site_filters: String,
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
        let site_filters = cols.get(10).unwrap_or(&"").trim().to_string();
        rows.push(BedInfo {
            site_filters,
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

/// Parse the verify stage's alt_verified column. Verify appends
/// (alt_verified, withdrawn) after scan's trailing filters column, whose
/// empty value breaks fixed positional indexing; the read is therefore from
/// the row's end: second-to-last field.
fn parse_alt_verified(text: &str) -> Vec<u32> {
    text.lines()
        .filter(|l| !l.is_empty())
        .filter_map(|l| l.rsplit('\t').nth(1).and_then(|v| v.parse().ok()))
        .collect()
}

/// Evidence filters that block PASS (FILTER column). Everything else the
/// scan annotates rides INFO/TAGS and never gates.
const PASS_BLOCKING: [&str; 2] = ["RETROCOPY", "SIMPLE_REPEAT"];

/// VCF 4.2 FILTER column: PASS/LOW_SCORE/UNSCORED appear alone; when any
/// blocking tag fired, the FILTER is exactly the blocking list.
fn filter_column(base: &str, blocking: &[&str]) -> String {
    if blocking.is_empty() {
        base.to_string()
    } else {
        blocking.join(";")
    }
}

/// Write `<out>/sites.vcf` (VCF v4.2).
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
    // Verified-vote counts (verify stage): sites whose re-aligned support
    // fell below the floor carry FILTER=REALIGN alongside the score filter.
    let verified: Option<Vec<u32>> = match entry {
        Entry::BamSites => None,
        _ => std::fs::read_to_string(params.out_dir.join("scan").join("candidates.verified.bed"))
            .ok()
            .map(|text| parse_alt_verified(&text)),
    };
    // Hyperedited-region sites (Stage::Hyper): rescued-read edit positions.
    // Sites already in candidates annotate INFO; the rest append as
    // FILTER=HYPER rows (never scored).
    #[allow(clippy::type_complexity)]
    let (hyper_info, hyper_extra): (
        std::collections::BTreeMap<(String, i64), (u32, String)>,
        Vec<((String, i64), String, String, u32, String)>,
    ) = match std::fs::read_to_string(params.out_dir.join("hyper").join("sites.tsv")) {
        Ok(text) => {
            let known: std::collections::HashSet<(String, i64)> =
                sites.iter().map(|(c, p)| (c.clone(), *p)).collect();
            let mut info = std::collections::BTreeMap::new();
            let mut extra = Vec::new();
            for line in text.lines() {
                if line.is_empty() {
                    continue;
                }
                let f: Vec<&str> = line.split('\t').collect();
                if f.len() < 5 {
                    continue;
                }
                let Ok(pos0) = f[1].parse::<i64>() else { continue };
                let (Ok(n),) = (f[2].parse::<u32>(),) else { continue };
                let key = (f[0].to_string(), pos0 + 1);
                info.insert(key.clone(), (n, f[3].to_string()));
                if !known.contains(&key) {
                    extra.push((key, f[4].to_string(), f[3].to_string(), n, String::new()));
                }
            }
            extra.sort_by(|a, b| a.0.cmp(&b.0));
            (info, extra)
        }
        Err(_) => (std::collections::BTreeMap::new(), Vec::new()),
    };

    if let (Some(v), Some(rows)) = (&verified, &bed_rows) {
        if v.len() != rows.len() {
            return Err(FlowError::BedParse {
                line: 0,
                msg: format!(
                    "candidates.verified.bed row count {} != candidates.bed {}",
                    v.len(),
                    rows.len()
                ),
            });
        }
    }
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
        "##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Editing-consistent alternate allele frequency\">\n",
    );
    out.push_str("##INFO=<ID=DEPTH,Number=1,Type=Integer,Description=\"Counted depth (MAPQ >= 20, base quality >= 30)\">\n");
    out.push_str("##INFO=<ID=AC,Number=4,Type=Integer,Description=\"Counted reads by base: A,C,G,T\">\n");
    out.push_str(
        "##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Transcript strand implied by the reference base\">\n",
    );
    out.push_str("##FILTER=<ID=PASS,Description=\"RE_PROB >= 0.5\">\n");
    out.push_str("##FILTER=<ID=LOW_SCORE,Description=\"RE_PROB < 0.5\">\n");
    out.push_str("##FILTER=<ID=UNSCORED,Description=\"site on a contig with no matching model (hybrid reference)\">\n");
    out.push_str("##FILTER=<ID=REALIGN,Description=\"fewer than 2 verified supporting reads after exhaustive re-alignment\">\n");
    out.push_str("##FILTER=<ID=HYPER,Description=\"site in a hyperedited region called from rescued reads (not model-scored)\">\n");
    out.push_str("##INFO=<ID=HYPER_READS,Number=1,Type=Integer,Description=\"rescued reads carrying this edit in the cluster\">\n");
    out.push_str("##INFO=<ID=HYPER_REGION,Number=1,Type=String,Description=\"hyperedited region chrom:start-end\">\n");
    out.push_str("##FILTER=<ID=SIMPLE_REPEAT,Description=\"inside a RepeatMasker simple-repeat or low-complexity region\">\n");
    out.push_str("##FILTER=<ID=RETROCOPY,Description=\"inside a mapped intronless retrocopy interval\">\n");
    out.push_str("##INFO=<ID=TAGS,Number=.,Type=String,Description=\"evidence annotations that never affect PASS (SPLICE_ADJ, HOMOPOLYMER, MULTICOPY)\">\n");
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

    let mut extras_i = 0usize;
    for (i, ((chrom, pos), prob)) in sites.iter().zip(probs).enumerate() {
        let (filter, mut info) = match prob {
            Some(p) => (
                if *p >= 0.5 { "PASS" } else { "LOW_SCORE" },
                format!("RE_PROB={p}"),
            ),
            None => ("UNSCORED", String::new()),
        };
        let refb = match &bed_rows {
            Some(rows) => rows[i].ref_base,
            None => ref_base(&fai, chrom, *pos),
        };
        let alt = edit_alt(refb).unwrap_or('N');
        if let Some(rows) = &bed_rows {
            let r = &rows[i];
            let sep = if info.is_empty() { "" } else { ";" };
            let _ = write!(
                info,
                "{sep}VAF={};DEPTH={};AC={},{},{},{};STRAND={}",
                r.edit_freq,
                r.depth,
                r.counts[0],
                r.counts[1],
                r.counts[2],
                r.counts[3],
                strand_of(refb)
            );
        } else {
            let sep = if info.is_empty() { "" } else { ";" };
            let _ = write!(info, "{sep}STRAND={}", strand_of(refb));
        }
        let realign = match (&verified, prob) {
            (Some(v), Some(_)) if v.get(i).is_some_and(|kept| *kept < 2) => {
                Some("REALIGN")
            }
            _ => None,
        };
        let site_tags = bed_rows
            .as_ref()
            .and_then(|rows| rows.get(i))
            .map(|r| r.site_filters.as_str())
            .unwrap_or("");
        // Only the two calibrated blacklists block PASS; the other evidence
        // tags annotate via INFO/TAGS (a filter may
        // block only when it removes >= 5x more false than true sites and
        // its true-site loss is not hyperedited-dominated).
        let mut blocking: Vec<&str> = Vec::new();
        let mut tag_info: Vec<&str> = Vec::new();
        for t in site_tags.split(',').filter(|t| !t.is_empty()) {
            if PASS_BLOCKING.contains(&t) {
                blocking.push(t);
            } else {
                tag_info.push(t);
            }
        }
        if let Some(r) = realign {
            blocking.push(r);
        }
        if !tag_info.is_empty() {
            let sep = if info.is_empty() { "" } else { ";" };
            let _ = write!(info, "{sep}TAGS={}", tag_info.join(","));
        }
        if let Some((n, region)) = hyper_info.get(&(chrom.clone(), *pos)) {
            let sep = if info.is_empty() { "" } else { ";" };
            let _ = write!(info, "{sep}HYPER_READS={n};HYPER_REGION={region}");
        }
        // VCF 4.2: PASS/LOW_SCORE/UNSCORED stand alone; any blocking tag
        // replaces the base verdict wholesale.
        let all = filter_column(filter, &blocking);
        // Flush sorted hyper-only rows that belong before this row.
        while extras_i < hyper_extra.len() && hyper_extra[extras_i].0 < (chrom.clone(), *pos) {
            let (key, dir, region, n, _) = &hyper_extra[extras_i];
            let (refb, alt) = (dir.chars().next().unwrap_or('N'), dir.chars().nth(1).unwrap_or('N'));
            let _ = writeln!(
                out,
                "{}\t{}\t.\t{}\t{}\t.\tHYPER\tHYPER_READS={};HYPER_REGION={}",
                key.0, key.1, refb, alt, n, region
            );
            extras_i += 1;
        }
        let _ = writeln!(out, "{chrom}\t{pos}\t.\t{refb}\t{alt}\t.\t{all}\t{info}");
    }
    for (key, dir, region, n, _) in &hyper_extra[extras_i.min(hyper_extra.len())..] {
        let (refb, alt) = (dir.chars().next().unwrap_or('N'), dir.chars().nth(1).unwrap_or('N'));
        let _ = writeln!(
            out,
            "{}\t{}\t.\t{}\t{}\t.\tHYPER\tHYPER_READS={};HYPER_REGION={}",
            key.0, key.1, refb, alt, n, region
        );
    }

    fs::write(params.out_dir.join("sites.vcf"), out)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn filter_column_never_mixes_pass_with_blockers() {
        assert_eq!(filter_column("PASS", &[]), "PASS");
        assert_eq!(filter_column("LOW_SCORE", &[]), "LOW_SCORE");
        assert_eq!(filter_column("UNSCORED", &[]), "UNSCORED");
        // A PASS-worthy site inside a retrocopy blacklist: FILTER is exactly
        // RETROCOPY — no PASS;RETROCOPY (VCF 4.2: PASS appears alone).
        assert_eq!(filter_column("PASS", &["RETROCOPY"]), "RETROCOPY");
        assert_eq!(
            filter_column("LOW_SCORE", &["RETROCOPY", "REALIGN"]),
            "RETROCOPY;REALIGN"
        );
    }

    #[test]
    fn alt_and_strand_follow_reference_base() {
        assert_eq!(edit_alt('A'), Some('G'));
        assert_eq!(edit_alt('T'), Some('C'));
        assert_eq!(edit_alt('C'), None);
        assert_eq!(strand_of('A'), '+');
        assert_eq!(strand_of('T'), '-');
    }

    #[test]
    fn parse_alt_verified_reads_production_layout() {
        // Production rows: scan's 11 columns (incl. a possibly empty filters
        // column) + verify's appended alt_verified + withdrawn.
        let text = "chr1\t100\t101\tA\t30\t20\t0\t10\t0\t0.3333\t\t7\t1\n\
                    chr1\t200\t201\tT\t25\t0\t15\t0\t10\t0.6000\tSPLICE_ADJ,HOMOPOLYMER\t3\t2\n";
        assert_eq!(parse_alt_verified(text), vec![7, 3]);
        // Legacy 10-column candidates + 2 appended (early fixtures).
        let text = "chr1\t100\t101\tA\t30\t20\t0\t10\t0\t0.3333\t5\t0\n";
        assert_eq!(parse_alt_verified(text), vec![5]);
    }
}
