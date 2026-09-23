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

/// One Tier-2-only (rescued) site row from `hyper/sites.tsv`.
#[derive(Clone)]
pub(crate) struct HyperRow {
    pub key: (String, i64),
    /// Editing direction, `AG` or `TC`.
    pub dir: String,
    pub region: String,
    pub n: u32,
}

/// Split `hyper/sites.tsv` into (a) INFO annotations for sites already in the
/// candidate table and (b) Tier-2-only rows, sorted by coordinate.
///
/// Strict on the two failure modes that would otherwise pass silently: a
/// duplicated (chrom, pos) key (the key excludes the direction; both the
/// candidate-annotation path and the tier-2 row path consume this table, and
/// neither may overwrite or double-emit) and a direction other than AG/TC
/// (the REF/ALT letters are taken from the direction, so a malformed value
/// would masquerade as a base call). Malformed numeric fields are errors for
/// the same reason. An absent file (no rescue output) is an empty table.
#[allow(clippy::type_complexity)]
pub(crate) fn hyper_rows(
    params: &RunParams,
    sites: &[(String, i64)],
) -> Result<
    (
        std::collections::BTreeMap<(String, i64), (u32, String)>,
        Vec<HyperRow>,
    ),
    FlowError,
> {
    let path = params.out_dir.join("hyper").join("sites.tsv");
    let text = match std::fs::read_to_string(&path) {
        Ok(text) => text,
        // Absent file = no rescue output (legal empty table); any other read
        // failure (permissions, I/O) is a hard error, never a fake empty set.
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            return Ok((std::collections::BTreeMap::new(), Vec::new()))
        }
        Err(e) => {
            return Err(FlowError::BedParse {
                line: 0,
                msg: format!("cannot read {}: {e}", path.display()),
            })
        }
    };
    let bad = |line: usize, msg: String| FlowError::BedParse {
        line,
        msg: format!("{}: {msg}", path.display()),
    };
    let known: std::collections::HashSet<(String, i64)> =
        sites.iter().map(|(c, p)| (c.clone(), *p)).collect();
    let mut info = std::collections::BTreeMap::new();
    let mut extra = Vec::new();
    for (i, line) in text.lines().enumerate() {
        let ln = i + 1;
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != 5 {
            return Err(bad(
                ln,
                format!("expected 5 tab-separated columns, got {}", f.len()),
            ));
        }
        let pos0 = f[1]
            .parse::<i64>()
            .map_err(|_| bad(ln, format!("pos '{}' is not an integer", f[1])))?;
        let n = f[2]
            .parse::<u32>()
            .map_err(|_| bad(ln, format!("support '{}' is not a non-negative integer", f[2])))?;
        if f[4] != "AG" && f[4] != "TC" {
            return Err(bad(
                ln,
                format!("direction '{}' (expected AG|TC)", f[4]),
            ));
        }
        let key = (f[0].to_string(), pos0 + 1);
        if let Some((prev_n, prev_region)) =
            info.insert(key.clone(), (n, f[3].to_string()))
        {
            return Err(bad(
                ln,
                format!(
                    "duplicate hyper site {}:{} ({} reads in {} vs {} reads in {}); refusing to overwrite, merge or drop rows",
                    key.0, key.1, prev_n, prev_region, n, f[3]
                ),
            ));
        }
        if !known.contains(&key) {
            extra.push(HyperRow {
                key,
                dir: f[4].to_string(),
                region: f[3].to_string(),
                n,
            });
        }
    }
    extra.sort_by(|a, b| a.key.cmp(&b.key));
    Ok((info, extra))
}

/// Decision thresholds the score stage hands to the VCF writer.
pub struct ScoreThresholds {
    /// PASS threshold for human-model RE_PROB (0.5 on legacy fusion bundles).
    pub pass: f64,
    /// Frozen rescue-model threshold for Tier-2-only rows (CMC bundles only).
    pub hyper: Option<f64>,
}


/// One Tier-2-only row: PASS above the frozen rescue threshold, FILTER=HYPER
/// below it or when unscored; HYPER_READS/HYPER_REGION always carried.
fn write_hyper_row(
    out: &mut String,
    row: &HyperRow,
    prob: Option<f64>,
    thresholds: &ScoreThresholds,
) {
    let (refb, alt) = (
        row.dir.chars().next().unwrap_or('N'),
        row.dir.chars().nth(1).unwrap_or('N'),
    );
    let pass = prob
        .zip(thresholds.hyper)
        .is_some_and(|(p, t)| p >= t);
    let (filter, mut info) = match prob {
        Some(p) => (if pass { "PASS" } else { "HYPER" }, format!("RE_PROB={p}")),
        None => ("HYPER", String::new()),
    };
    let sep = if info.is_empty() { "" } else { ";" };
    let _ = write!(
        info,
        "{sep}HYPER_READS={};HYPER_REGION={}",
        row.n, row.region
    );
    let _ = writeln!(
        out,
        "{}\t{}\t.\t{}\t{}\t.\t{filter}\t{info}",
        row.key.0, row.key.1, refb, alt
    );
}

/// Write `<out>/sites.vcf` (VCF v4.2).
pub fn write_vcf(
    params: &RunParams,
    entry: Entry,
    sites: &[(String, i64)],
    probs: &[Option<f64>],
    hyper_probs: &[Option<f64>],
    thresholds: &ScoreThresholds,
) -> Result<(), FlowError> {
    let bed_rows: Option<Vec<BedInfo>> = match entry {
        Entry::BamSites => None,
        _ => Some(parse_bed_info(
            &params.out_dir.join("scan").join("candidates.bed"),
        )?),
    };
    // Verified-vote counts (verify stage): sites whose re-aligned support
    // fell below the floor carry FILTER=REALIGN alongside the score filter.
    // Absent file = legacy candidates.bed without the verified columns
    // (legal); a present-but-unreadable file is a hard error.
    let verified: Option<Vec<u32>> = match entry {
        Entry::BamSites => None,
        _ => {
            let vpath = params.out_dir.join("scan").join("candidates.verified.bed");
            match std::fs::read_to_string(&vpath) {
                Ok(text) => Some(parse_alt_verified(&text)),
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => None,
                Err(e) => {
                    return Err(FlowError::BedParse {
                        line: 0,
                        msg: format!("cannot read {}: {e}", vpath.display()),
                    })
                }
            }
        }
    };
    // Hyperedited-region sites (Stage::Hyper): rescued-read edit positions.
    // Sites already in candidates annotate INFO; the rest append as Tier-2-only
    // rows, scored by the rescue model on CMC bundles.
    let (hyper_info, hyper_extra) = hyper_rows(params, sites)?;
    // Hybrid runs: mouse-contig sites keep the legacy mouse-model threshold.
    let manifest = crate::manifest::SpeciesManifest::read(&params.out_dir)?;
    let row_threshold = |chrom: &str| -> f64 {
        if manifest
            .as_ref()
            .is_some_and(|m| m.kind == "hybrid" && m.owner_of(chrom) == Some("mouse"))
        {
            0.5
        } else {
            thresholds.pass
        }
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
    let pass_desc = if thresholds.pass == 0.5 {
        "0.5".to_string()
    } else {
        format!("{:.4}", thresholds.pass)
    };
    let _ = writeln!(
        out,
        "##FILTER=<ID=PASS,Description=\"RE_PROB >= {pass_desc}\">"
    );
    let _ = writeln!(
        out,
        "##FILTER=<ID=LOW_SCORE,Description=\"RE_PROB < {pass_desc}\">"
    );
    out.push_str("##FILTER=<ID=UNSCORED,Description=\"site on a contig with no matching model (hybrid reference)\">\n");
    out.push_str("##FILTER=<ID=REALIGN,Description=\"fewer than 2 verified supporting reads after exhaustive re-alignment\">\n");
    out.push_str("##FILTER=<ID=HYPER,Description=\"tier-2-only site in a hyperedited region called from rescued reads, below the rescue model threshold (or unscored)\">\n");
    out.push_str("##INFO=<ID=HYPER_READS,Number=1,Type=Integer,Description=\"rescued reads carrying this edit in the cluster\">\n");
    out.push_str("##INFO=<ID=HYPER_REGION,Number=1,Type=String,Description=\"hyperedited region chrom:start-end\">\n");
    out.push_str("##FILTER=<ID=SIMPLE_REPEAT,Description=\"inside a RepeatMasker simple-repeat or low-complexity region\">\n");
    out.push_str("##FILTER=<ID=RETROCOPY,Description=\"inside a mapped intronless retrocopy interval\">\n");
    out.push_str("##INFO=<ID=TAGS,Number=.,Type=String,Description=\"evidence annotations that never affect PASS (SPLICE_ADJ, HOMOPOLYMER, MULTICOPY)\">\n");
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

    let mut extras_i = 0usize;
    for (i, ((chrom, pos), prob)) in sites.iter().zip(probs).enumerate() {
        let row_t = row_threshold(chrom);
        let (filter, mut info) = match prob {
            Some(p) => (
                if *p >= row_t { "PASS" } else { "LOW_SCORE" },
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
        while extras_i < hyper_extra.len()
            && hyper_extra[extras_i].key < (chrom.clone(), *pos)
        {
            let row = &hyper_extra[extras_i];
            write_hyper_row(
                &mut out,
                row,
                hyper_probs.get(extras_i).copied().flatten(),
                thresholds,
            );
            extras_i += 1;
        }
        let _ = writeln!(out, "{chrom}\t{pos}\t.\t{refb}\t{alt}\t.\t{all}\t{info}");
    }
    let leftover = extras_i.min(hyper_extra.len());
    for (i, row) in hyper_extra[leftover..].iter().enumerate() {
        write_hyper_row(
            &mut out,
            row,
            hyper_probs.get(leftover + i).copied().flatten(),
            thresholds,
        );
    }

    fs::write(params.out_dir.join("sites.vcf"), out)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hyper_rows_render_against_the_frozen_threshold() {
        let row = HyperRow {
            key: ("chr1".into(), 1000),
            dir: "AG".into(),
            region: "chr1:900-1100".into(),
            n: 4,
        };
        let t = ScoreThresholds {
            pass: 0.2854,
            hyper: Some(0.3624),
        };
        let mut out = String::new();
        write_hyper_row(&mut out, &row, Some(0.9), &t);
        assert_eq!(
            out,
            "chr1\t1000\t.\tA\tG\t.\tPASS\tRE_PROB=0.9;HYPER_READS=4;HYPER_REGION=chr1:900-1100\n"
        );
        let mut out = String::new();
        write_hyper_row(&mut out, &row, Some(0.1), &t);
        assert!(out.contains("\tHYPER\t"));
        assert!(out.contains("RE_PROB=0.1"));
        // unscored: legacy rendering, no RE_PROB
        let mut out = String::new();
        write_hyper_row(&mut out, &row, None, &t);
        assert_eq!(
            out,
            "chr1\t1000\t.\tA\tG\t.\tHYPER\tHYPER_READS=4;HYPER_REGION=chr1:900-1100\n"
        );
    }

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

#[cfg(test)]
mod hyper_rows_strict_tests {
    use super::*;
    use std::path::PathBuf;

    fn params_with_sites_tsv(tag: &str, body: &str) -> (PathBuf, RunParams) {
        let d = std::env::temp_dir().join(format!("vcf_hr_{tag}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(d.join("hyper")).unwrap();
        std::fs::write(d.join("hyper").join("sites.tsv"), body).unwrap();
        let params = RunParams {
            r1: Vec::new(),
            r2: Vec::new(),
            bam: None,
            sites: None,
            index: None,
            fasta: PathBuf::from("ref.fa"),
            gtf: None,
            snp: None,
            bundle: PathBuf::from("bundle"),
            l1_bundle: None,
            out_dir: d.clone(),
            threads: 1,
            batch: 64,
            rescue_control: false,
            stop_after: None,
        };
        (d, params)
    }

    #[test]
    fn duplicate_key_is_a_hard_error_on_the_read_side() {
        let (d, p) = params_with_sites_tsv(
            "dup",
            "chr1\t100\t3\tchr1:100-200\tTC\nchr1\t100\t7\tchr1:90-210\tTC\n",
        );
        let err = hyper_rows(&p, &[]).err().unwrap();
        assert!(err.to_string().contains("duplicate hyper site chr1:101"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn unknown_direction_is_rejected_not_truncated_into_bases() {
        let (d, p) = params_with_sites_tsv("dir", "chr1\t100\t3\tchr1:100-200\tAGx\n");
        let err = hyper_rows(&p, &[]).err().unwrap();
        assert!(err.to_string().contains("direction 'AGx'"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn malformed_numeric_fields_are_errors() {
        let (d, p) = params_with_sites_tsv("num", "chr1\tx\t3\tchr1:100-200\tAG\n");
        assert!(hyper_rows(&p, &[]).is_err());
        let _ = std::fs::remove_dir_all(&d);
        let (d, p) = params_with_sites_tsv("num2", "chr1\t100\tx\tchr1:100-200\tAG\n");
        assert!(hyper_rows(&p, &[]).is_err());
        let _ = std::fs::remove_dir_all(&d);
        let (d, p) = params_with_sites_tsv("cols", "chr1\t100\t3\n");
        assert!(hyper_rows(&p, &[]).is_err());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn extra_columns_are_rejected() {
        let (d, p) = params_with_sites_tsv("extra", "chr1\t100\t3\tchr1:100-200\tTC\tsurplus\n");
        let err = hyper_rows(&p, &[]).err().unwrap();
        assert!(err.to_string().contains("expected 5 tab-separated columns, got 6"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn unreadable_sites_tsv_is_an_error_not_an_empty_table() {
        let (d, p) = params_with_sites_tsv("unread", "chr1\t100\t3\tchr1:100-200\tTC\n");
        let path = d.join("hyper").join("sites.tsv");
        // chmod 000: present but unreadable (no-op when running as root).
        std::process::Command::new("chmod").arg("000").arg(&path).status().unwrap();
        match hyper_rows(&p, &[]) {
            Err(e) => assert!(e.to_string().contains("cannot read"), "{e}"),
            Ok(_) if nix_free_root() => {} // running as root: chmod 000 still readable
            Ok(_) => panic!("unreadable sites.tsv must not parse as an empty table"),
        }
        std::process::Command::new("chmod").arg("644").arg(&path).status().unwrap();
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn absent_sites_tsv_is_a_legal_empty_table() {
        let d = std::env::temp_dir().join(format!("vcf_hr_absent_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        let mut params = RunParams {
            r1: Vec::new(),
            r2: Vec::new(),
            bam: None,
            sites: None,
            index: None,
            fasta: PathBuf::from("ref.fa"),
            gtf: None,
            snp: None,
            bundle: PathBuf::from("bundle"),
            l1_bundle: None,
            out_dir: d.clone(),
            threads: 1,
            batch: 64,
            rescue_control: false,
            stop_after: None,
        };
        let (info, extra) = hyper_rows(&params, &[]).unwrap();
        assert!(info.is_empty() && extra.is_empty());
        params.out_dir = d.clone();
        let _ = std::fs::remove_dir_all(&d);
    }

    fn nix_free_root() -> bool {
        // true when euid != 0 (chmod 000 actually blocks reads)
        std::process::Command::new("id")
            .arg("-u")
            .output()
            .map(|o| String::from_utf8_lossy(&o.stdout).trim() != "0")
            .unwrap_or(true)
    }

    #[test]
    fn clean_table_splits_annotation_and_tier2() {
        let (d, p) = params_with_sites_tsv(
            "clean",
            "chr1\t100\t3\tchr1:100-200\tTC\nchr1\t300\t4\tchr1:250-350\tAG\n",
        );
        let known = vec![("chr1".to_string(), 101i64)];
        let (info, extra) = hyper_rows(&p, &known).unwrap();
        assert_eq!(info.len(), 2);
        assert_eq!(extra.len(), 1);
        assert_eq!(extra[0].key, ("chr1".to_string(), 301));
        assert_eq!(extra[0].dir, "AG");
        let _ = std::fs::remove_dir_all(&d);
    }
}
