//! Hyperedited-region caller over the collapsed-rescue record
//! (`map/rescued.tsv`, written by the rescue stage with per-read edit
//! positions). Reads sharing (chrom, edit direction) are clustered when
//! their reference spans sit within 150 bp; clusters with >= 3 reads become
//! regions. Products: `hyper/regions.bed` (region track for the report) and
//! `hyper/sites.tsv` (per-position edit support for the VCF's HYPER rows).

use std::io::Error as IoError;
use std::path::Path;

/// Maximum reference-span gap between reads in one cluster.
pub const CLUSTER_GAP_BP: u32 = 150;
/// Minimum reads per cluster.
pub const CLUSTER_MIN_READS: usize = 3;
/// Minimum reads carrying an edit position for a HYPER site row.
pub const SITE_MIN_READS: usize = 2;

struct RescuedRow {
    chrom: String,
    pos0: u32,
    /// Placement strand; load-bearing for canonical-polarity grouping.
    strand: char,
    dir: String,
    edits: Vec<u32>,
}

/// Header of `map/rescued.tsv` (writer: the rescue stage).
const RESCUED_HEADER: &str = "qname\tchrom\tpos0\tstrand\tdir\tn_edits\tedit_pos0\tmapq";

/// Strict reader for the rescue record: every field is validated and any
/// malformed row is a hard error carrying the file, the 1-based line number
/// and the reason. A header-only file is a legal empty record set; a
/// truncated or corrupted one must never be read as fewer records than it
/// holds (the previous lenient reader skipped bad rows and dropped
/// unparseable edit coordinates).
fn read_rescued_tsv(path: &Path) -> Result<Vec<RescuedRow>, IoError> {
    let text = std::fs::read_to_string(path)?;
    let bad = |line: usize, msg: String| {
        IoError::other(format!(
            "malformed rescued.tsv {} at line {line}: {msg}",
            path.display()
        ))
    };
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| bad(1, "missing header".to_string()))?;
    if header != RESCUED_HEADER {
        return Err(bad(1, format!("bad header (expected '{RESCUED_HEADER}')")));
    }
    let mut out = Vec::new();
    for (i, line) in lines.enumerate() {
        let ln = i + 2;
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != 8 {
            return Err(bad(
                ln,
                format!("expected 8 tab-separated columns, got {}", f.len()),
            ));
        }
        let pos0 = f[2]
            .parse::<u32>()
            .map_err(|_| bad(ln, format!("pos0 '{}' is not a non-negative integer", f[2])))?;
        let strand = match f[3] {
            "+" => '+',
            "-" => '-',
            other => return Err(bad(ln, format!("strand '{other}' (expected +|-)"))),
        };
        let dir = match f[4] {
            "AG" | "TC" => f[4].to_string(),
            other => return Err(bad(ln, format!("dir '{other}' (expected AG|TC)"))),
        };
        let n_edits = f[5]
            .parse::<usize>()
            .map_err(|_| bad(ln, format!("n_edits '{}' is not a non-negative integer", f[5])))?;
        let edits: Vec<u32> = if f[6].is_empty() {
            Vec::new()
        } else {
            f[6]
                .split(',')
                .map(|x| {
                    x.parse::<u32>().map_err(|_| {
                        bad(ln, format!("edit coordinate '{x}' is not a non-negative integer"))
                    })
                })
                .collect::<Result<_, _>>()?
        };
        if n_edits != edits.len() {
            return Err(bad(
                ln,
                format!("n_edits {n_edits} != {} edit coordinates", edits.len()),
            ));
        }
        f[7]
            .parse::<u8>()
            .map_err(|_| bad(ln, format!("mapq '{}' is not a 0-255 integer", f[7])))?;
        out.push(RescuedRow {
            chrom: f[1].to_string(),
            pos0,
            strand,
            dir,
            edits,
        });
    }
    Ok(out)
}

struct Cluster {
    start: u32,
    end: u32,
    reads: Vec<usize>, // indices into rows
}

/// Run the caller; returns (n_regions, n_sites). `fasta` (with .fai) feeds
/// the edit-index denominator and the polarity tripwire; without it the
/// edit-index column is written as NaN and the tripwire report marks the
/// reference as unavailable. A present-but-unreadable reference region is a
/// hard error. Products: `hyper/regions.bed` (region track for the report),
/// `hyper/sites.tsv` (per-position edit support for the VCF's HYPER rows)
/// and `hyper/polarity_report.json` (tripwire outcome; mismatches also raise
/// a stderr warning and never fail the run).
pub fn run(map_dir: &Path, out_dir: &Path, fasta: &Path) -> Result<(usize, usize), IoError> {
    let fai = rust_htslib::faidx::Reader::from_path(fasta).ok();
    // faidx_fetch_seq aborts on an unknown contig instead of returning an
    // error (htslib contract), so guard every fetch with an explicit
    // presence/range check: a reference that cannot be read is a hard error,
    // never a silent 0.
    let fai_contigs: Option<std::collections::HashSet<String>> = fai.as_ref().map(|f| {
        f.seq_names().unwrap_or_default().into_iter().collect()
    });
    let fetch_ref = |chrom: &str, begin: usize, end: usize, what: &str| -> Result<Vec<u8>, IoError> {
        let f = fai.as_ref().expect("caller checked");
        let names = fai_contigs.as_ref().expect("caller checked");
        let len = match names.contains(chrom) {
            true => f.fetch_seq_len(chrom) as usize,
            false => {
                return Err(IoError::other(format!(
                    "reference read failed for {chrom}:{begin}-{end} ({what}): contig absent from the fasta index"
                )))
            }
        };
        if begin >= len || end >= len {
            return Err(IoError::other(format!(
                "reference read failed for {chrom}:{begin}-{end} ({what}): outside contig length {len}"
            )));
        }
        f.fetch_seq(chrom, begin, end).map_err(|e| {
            IoError::other(format!(
                "reference read failed for {chrom}:{begin}-{end} ({what}): {e}"
            ))
        })
    };
    let rows = read_rescued_tsv(&map_dir.join("rescued.tsv"))?;
    // Genomic polarity of a rescued read is fixed by (dir, strand): the
    // verifier complements the expectation pair on minus placements, so an
    // AG-pass minus read and a TC-pass plus read over the same locus carry
    // the SAME genomic event. Grouping by the raw pass label split mate
    // pairs of one fragment into two groups (double-counted read-end
    // support, half of the rows written with the wrong REF/ALT polarity).
    // Normalize each read to its genomic polarity before grouping; the
    // sidecar keeps the raw pass label as the audit trail.
    let mut sites_checked = 0usize;
    let mut mismatch_sites: Vec<(String, u32)> = Vec::new();
    let canonical_dir = |dir: &str, strand: char| -> Result<String, IoError> {
        match (dir, strand) {
            ("AG", '+') | ("TC", '-') => Ok("AG".to_string()),
            ("TC", '+') | ("AG", '-') => Ok("TC".to_string()),
            (other, _) => Err(IoError::other(format!(
                "unknown rescue direction '{other}' (expected AG|TC)"
            ))),
        }
    };
    let mut by_group: std::collections::BTreeMap<(String, String), Vec<usize>> =
        std::collections::BTreeMap::new();
    for (i, r) in rows.iter().enumerate() {
        let canonical = canonical_dir(&r.dir, r.strand)?;
        by_group
            .entry((r.chrom.clone(), canonical))
            .or_default()
            .push(i);
    }

    let mut regions: Vec<String> = Vec::new();
    let mut sites: Vec<String> = Vec::new();
    let mut n_regions = 0usize;
    let mut n_sites = 0usize;
    // Duplicate guard state: one row per (chrom, pos) across all groups —
    // the key deliberately excludes the direction. Value: (dir, region) of
    // the first occurrence, for the error message.
    let mut seen_sites: std::collections::BTreeMap<(String, u32), (String, String)> =
        std::collections::BTreeMap::new();

    for ((chrom, dir), mut idxs) in by_group {
        idxs.sort_by_key(|&i| rows[i].pos0);
        let mut clusters: Vec<Cluster> = Vec::new();
        for i in idxs {
            let r = &rows[i];
            let r_end = r
                .edits
                .iter()
                .max()
                .map(|&e| e + 1)
                .unwrap_or(r.pos0 + 1)
                .max(r.pos0 + 1);
            match clusters.last_mut() {
                Some(c) if r.pos0 <= c.end + CLUSTER_GAP_BP => {
                    c.end = c.end.max(r_end);
                    c.reads.push(i);
                }
                _ => clusters.push(Cluster {
                    start: r.pos0,
                    end: r_end,
                    reads: vec![i],
                }),
            }
        }
        for c in clusters
            .iter()
            .filter(|c| c.reads.len() >= CLUSTER_MIN_READS)
        {
            // Per-position edit support across the cluster's reads.
            let mut pos_count: std::collections::BTreeMap<u32, usize> =
                std::collections::BTreeMap::new();
            for &i in &c.reads {
                for &e in &rows[i].edits {
                    *pos_count.entry(e).or_insert(0) += 1;
                }
            }
            let kept: Vec<u32> = pos_count
                .iter()
                .filter(|(_, &n)| n >= SITE_MIN_READS)
                .map(|(&p, _)| p)
                .collect();
            // Edit index: edited source bases / covered source bases of the
            // editing source letter (A for AG, T for TC), from the reference.
            let src = match dir.as_str() {
                "AG" => b'A',
                "TC" => b'T',
                other => {
                    return Err(IoError::other(format!(
                        "unknown rescue direction '{other}' (expected AG|TC)"
                    )))
                }
            };
            let edit_index = match &fai {
                Some(_) => {
                    // A reference that cannot be read here is an explicit
                    // error, never a silent 0 denominator (which would print
                    // as a 0.0000 edit index indistinguishable from "no
                    // source bases covered").
                    let span = fetch_ref(
                        &chrom,
                        c.start as usize,
                        (c.end - 1) as usize,
                        "edit-index denominator",
                    )?
                    .iter()
                    .filter(|&&b| b.eq_ignore_ascii_case(&src))
                    .count();
                    if span > 0 {
                        kept.len() as f64 / span as f64
                    } else {
                        0.0
                    }
                }
                None => f64::NAN,
            };
            n_regions += 1;
            let region_id = format!("{}:{}-{}", chrom, c.start, c.end);
            regions.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\n",
                chrom,
                c.start,
                c.end,
                dir,
                c.reads.len(),
                kept.len(),
                edit_index
            ));
            for &p in &kept {
                n_sites += 1;
                // Tripwire: with canonical grouping the reference base at a
                // kept position must be the group's source base. The verifier
                // contract makes a mismatch impossible unless upstream
                // semantics drift; count and report rather than fail the run.
                // A reference read failure is a hard error — it must never
                // masquerade as "0 mismatches".
                if fai.is_some() {
                    let base = fetch_ref(&chrom, p as usize, p as usize, "polarity tripwire")?;
                    sites_checked += 1;
                    if !base.first().is_some_and(|b| b.eq_ignore_ascii_case(&src)) {
                        mismatch_sites.push((chrom.clone(), p));
                    }
                }
                match seen_sites.entry((chrom.clone(), p)) {
                    std::collections::btree_map::Entry::Vacant(v) => {
                        v.insert((dir.clone(), region_id.clone()));
                    }
                    std::collections::btree_map::Entry::Occupied(e) => {
                        let (pd, pr) = e.get();
                        let kind = if *pd == dir && pr == &region_id {
                            "identical"
                        } else {
                            "conflicting"
                        };
                        return Err(IoError::other(format!(
                            "duplicate hyper site {chrom}:{p} ({kind}): {pd} {pr} vs {dir} {region_id}; refusing to overwrite, merge or drop rows"
                        )));
                    }
                }
                sites.push(format!(
                    "{}\t{}\t{}\t{}\t{}\n",
                    chrom, p, pos_count[&p], region_id, dir
                ));
            }
        }
    }

    let dir = out_dir.join("hyper");
    std::fs::create_dir_all(&dir)?;
    // Machine-readable polarity tripwire report. When the reference index is
    // unavailable the report says so explicitly (`sites_checked: 0` is then
    // "not checked", never "0 mismatches"); the NaN edit-index column is the
    // documented contract for that case.
    let report = serde_json::json!({
        "reference": if fai.is_some() { "ok" } else { "unavailable" },
        "sites_checked": sites_checked,
        "mismatches": mismatch_sites.len(),
        "mismatch_sites": mismatch_sites
            .iter()
            .map(|(c, p)| format!("{c}:{p}"))
            .collect::<Vec<_>>(),
    });
    std::fs::write(
        dir.join("polarity_report.json"),
        format!("{}\n", serde_json::to_string_pretty(&report)?),
    )?;
    if !mismatch_sites.is_empty() {
        const SHOW: usize = 20;
        let listed = mismatch_sites
            .iter()
            .take(SHOW)
            .map(|(c, p)| format!("{c}:{p}"))
            .collect::<Vec<_>>()
            .join(", ");
        let more = mismatch_sites.len().saturating_sub(SHOW);
        eprintln!(
            "[hyper] warning: {} of {} kept sites carry a reference base inconsistent with their group's editing polarity: {}{} (full list in hyper/polarity_report.json)",
            mismatch_sites.len(),
            sites_checked,
            listed,
            if more > 0 {
                format!(" and {more} more")
            } else {
                String::new()
            }
        );
    }
    std::fs::write(dir.join("regions.bed"), regions.concat())?;
    std::fs::write(dir.join("sites.tsv"), sites.concat())?;
    Ok((n_regions, n_sites))
}

#[cfg(test)]
mod canonical_tests {
    use super::*;
    use std::io::Write;

    /// Mates of one fragment rescued by opposite passes (TC/+ and AG/-)
    /// over the same genomic T->C event must group into ONE canonical
    /// polarity group with merged per-position support (C2-② regression).
    #[test]
    fn opposite_pass_mates_merge_into_one_canonical_group() {
        let tmp = std::env::temp_dir().join(format!("hyper_canon_{}", std::process::id()));
        std::fs::create_dir_all(&tmp).unwrap();
        let map = tmp.join("map");
        std::fs::create_dir_all(&map).unwrap();
        // 4 fragments x 2 mates over the same 3 edit positions; /1 mates
        // rescued by the TC pass on plus, /2 mates by the AG pass on minus.
        let mut tsv = String::from("qname\tchrom\tpos0\tstrand\tdir\tn_edits\tedit_pos0\tmapq\n");
        for i in 0..4 {
            tsv.push_str(&format!(
                "frag{i}/1\tctg1\t100\t+\tTC\t3\t200,210,220\t15\n"
            ));
            tsv.push_str(&format!(
                "frag{i}/2\tctg1\t100\t-\tAG\t3\t200,210,220\t15\n"
            ));
        }
        std::fs::write(map.join("rescued.tsv"), tsv).unwrap();
        // fasta without .fai: grouping still runs (fai only feeds edit-index)
        let fa = tmp.join("ref.fa");
        let mut f = std::fs::File::create(&fa).unwrap();
        writeln!(f, ">ctg1").unwrap();
        writeln!(f, "{}", "A".repeat(300)).unwrap();

        let (n_regions, n_sites) = run(&map, &tmp, &fa).unwrap();
        let sites = std::fs::read_to_string(tmp.join("hyper").join("sites.tsv")).unwrap();
        let lines = sites.lines();
        // one canonical group (TC: genome T->C polarity), 3 positions,
        // 8 read-ends of support each, single dir per position
        let mut pos_dir = std::collections::BTreeMap::new();
        for l in lines {
            let f: Vec<&str> = l.split('\t').collect();
            assert_eq!(f.len(), 5);
            assert_eq!(f[4], "TC", "canonical polarity must be TC for genomic T");
            pos_dir.insert(f[1].to_string(), f[2].parse::<u32>().unwrap());
        }
        assert_eq!(pos_dir.len(), 3, "three edit positions, one row each");
        assert_eq!(n_sites, 3);
        assert_eq!(n_regions, 1);
        assert!(pos_dir.values().all(|&n| n == 8), "8 read-ends per position");
        let _ = std::fs::remove_dir_all(&tmp);
    }

    /// A plus AG-pass read normalizes to AG (genome A->G polarity).
    #[test]
    fn plus_ag_read_keeps_ag_polarity() {
        let tmp = std::env::temp_dir().join(format!("hyper_plus_{}", std::process::id()));
        std::fs::create_dir_all(tmp.join("map")).unwrap();
        let mut tsv = String::from("qname\tchrom\tpos0\tstrand\tdir\tn_edits\tedit_pos0\tmapq\n");
        for i in 0..4 {
            tsv.push_str(&format!("r{i}/1\tctg1\t100\t+\tAG\t3\t200,210,220\t15\n"));
        }
        std::fs::write(tmp.join("map").join("rescued.tsv"), tsv).unwrap();
        let fa = tmp.join("ref.fa");
        let mut f = std::fs::File::create(&fa).unwrap();
        writeln!(f, ">ctg1").unwrap();
        writeln!(f, "{}", "A".repeat(300)).unwrap();
        run(&tmp.join("map"), &tmp, &fa).unwrap();
        let sites = std::fs::read_to_string(tmp.join("hyper").join("sites.tsv")).unwrap();
        assert!(sites.lines().all(|l| l.split('\t').nth(4) == Some("TC") || l.split('\t').nth(4) == Some("AG")));
        // plus AG stays AG
        assert!(sites.lines().all(|l| l.split('\t').nth(4) == Some("AG")));
        let _ = std::fs::remove_dir_all(&tmp);
    }
}

#[cfg(test)]
mod strict_tsv_tests {
    use super::*;
    use std::io::Write;

    fn tsv_dir(tag: &str, body: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("hyper_tsv_{tag}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        let mut f = std::fs::File::create(d.join("rescued.tsv")).unwrap();
        f.write_all(body.as_bytes()).unwrap();
        d
    }

    const HEADER: &str = "qname\tchrom\tpos0\tstrand\tdir\tn_edits\tedit_pos0\tmapq\n";
    const ROW: &str = "r1\tctg1\t100\t+\tAG\t2\t200,210\t15\n";

    #[test]
    fn header_only_is_a_legal_empty_set() {
        let d = tsv_dir("empty", HEADER);
        assert_eq!(read_rescued_tsv(&d.join("rescued.tsv")).unwrap().len(), 0);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn well_formed_rows_parse() {
        let d = tsv_dir("good", &format!("{HEADER}{ROW}"));
        let rows = read_rescued_tsv(&d.join("rescued.tsv")).unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].edits, vec![200, 210]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn malformed_rows_error_with_line_and_reason() {
        let cases: [(&str, &str, &str); 9] = [
            ("empty_file", "", "missing header"),
            ("bad_header", "qname\tchrom\n", "bad header"),
            ("few_cols", "r1\tctg1\t100\n", "expected 8 tab-separated columns, got 3"),
            ("bad_pos0", "r1\tctg1\tx\t+\tAG\t2\t200,210\t15\n", "pos0 'x'"),
            ("bad_strand", "r1\tctg1\t100\t.\tAG\t2\t200,210\t15\n", "strand '.'"),
            ("bad_dir", "r1\tctg1\t100\t+\tAC\t2\t200,210\t15\n", "dir 'AC'"),
            ("bad_nedits", "r1\tctg1\t100\t+\tAG\tx\t200,210\t15\n", "n_edits 'x'"),
            ("bad_coord", "r1\tctg1\t100\t+\tAG\t2\t200,x\t15\n", "edit coordinate 'x'"),
            ("bad_mapq", "r1\tctg1\t100\t+\tAG\t2\t200,210\tx\n", "mapq 'x'"),
        ];
        for (tag, row, want) in cases {
            let body = if tag == "empty_file" || tag == "bad_header" {
                row.to_string()
            } else {
                format!("{HEADER}{row}")
            };
            let d = tsv_dir(tag, &body);
            let err = read_rescued_tsv(&d.join("rescued.tsv")).err().unwrap();
            let msg = err.to_string();
            assert!(msg.contains("rescued.tsv"), "{tag}: path missing: {msg}");
            assert!(msg.contains("line"), "{tag}: line missing: {msg}");
            assert!(msg.contains(want), "{tag}: reason '{want}' missing: {msg}");
            let _ = std::fs::remove_dir_all(&d);
        }
    }

    #[test]
    fn n_edits_must_match_coordinate_count() {
        let d = tsv_dir("count", &format!("{HEADER}r1\tctg1\t100\t+\tAG\t3\t200,210\t15\n"));
        let err = read_rescued_tsv(&d.join("rescued.tsv")).err().unwrap();
        assert!(err.to_string().contains("n_edits 3 != 2 edit coordinates"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn empty_edit_list_pairs_with_zero_n_edits() {
        let d = tsv_dir("zeroedits", &format!("{HEADER}r1\tctg1\t100\t+\tAG\t0\t\t15\n"));
        let rows = read_rescued_tsv(&d.join("rescued.tsv")).unwrap();
        assert_eq!(rows.len(), 1);
        assert!(rows[0].edits.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }
}

#[cfg(test)]
mod tripwire_tests {
    use super::*;
    use std::io::Write;

    fn dir(tag: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("hyper_trip_{tag}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(d.join("map")).unwrap();
        d
    }

    fn fasta_with_fai(d: &std::path::Path) -> std::path::PathBuf {
        let fa = d.join("ref.fa");
        let mut f = std::fs::File::create(&fa).unwrap();
        // T at 200/210/220 (tripwire bait for an AG group), A elsewhere.
        let mut seq = vec![b'A'; 300];
        for p in [200usize, 210, 220] {
            seq[p] = b'T';
        }
        writeln!(f, ">ctg1").unwrap();
        f.write_all(&seq).unwrap();
        writeln!(f).unwrap();
        drop(f);
        rust_htslib::faidx::build(&fa).unwrap();
        fa
    }

    const HEADER: &str = "qname\tchrom\tpos0\tstrand\tdir\tn_edits\tedit_pos0\tmapq\n";

    #[test]
    fn tripwire_reports_mismatches_and_continues() {
        let d = dir("fires");
        // Canonical AG group (AG/+) whose kept positions hold reference T.
        let mut tsv = String::from(HEADER);
        for i in 0..3 {
            tsv.push_str(&format!("r{i}\tctg1\t100\t+\tAG\t3\t200,210,220\t15\n"));
        }
        std::fs::write(d.join("map").join("rescued.tsv"), tsv).unwrap();
        let fa = fasta_with_fai(&d);
        let (nr, ns) = run(&d.join("map"), &d, &fa).unwrap();
        assert_eq!((nr, ns), (1, 3));
        let report: serde_json::Value = serde_json::from_str(
            &std::fs::read_to_string(d.join("hyper").join("polarity_report.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(report["reference"], "ok");
        assert_eq!(report["sites_checked"], 3);
        assert_eq!(report["mismatches"], 3);
        assert_eq!(
            report["mismatch_sites"],
            serde_json::json!(["ctg1:200", "ctg1:210", "ctg1:220"])
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn tripwire_clean_run_reports_zero_checked_mismatches() {
        let d = dir("clean");
        let mut tsv = String::from(HEADER);
        for i in 0..3 {
            tsv.push_str(&format!("r{i}\tctg1\t100\t+\tAG\t3\t240,250,260\t15\n"));
        }
        std::fs::write(d.join("map").join("rescued.tsv"), tsv).unwrap();
        let fa = fasta_with_fai(&d);
        run(&d.join("map"), &d, &fa).unwrap();
        let report: serde_json::Value = serde_json::from_str(
            &std::fs::read_to_string(d.join("hyper").join("polarity_report.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(report["sites_checked"], 3);
        assert_eq!(report["mismatches"], 0);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn unreadable_reference_region_is_a_hard_error() {
        let d = dir("badref");
        let mut tsv = String::from(HEADER);
        for i in 0..3 {
            tsv.push_str(&format!("r{i}\tctgNOPE\t100\t+\tAG\t3\t200,210,220\t15\n"));
        }
        std::fs::write(d.join("map").join("rescued.tsv"), tsv).unwrap();
        let fa = fasta_with_fai(&d);
        let err = run(&d.join("map"), &d, &fa).err().unwrap();
        assert!(
            err.to_string().contains("reference read failed"),
            "unexpected error: {err}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn out_of_range_region_is_a_hard_error() {
        let d = dir("oor");
        let mut tsv = String::from(HEADER);
        for i in 0..3 {
            tsv.push_str(&format!("r{i}\tctg1\t280\t+\tAG\t3\t295,310,320\t15\n"));
        }
        std::fs::write(d.join("map").join("rescued.tsv"), tsv).unwrap();
        let fa = fasta_with_fai(&d);
        let err = run(&d.join("map"), &d, &fa).err().unwrap();
        assert!(
            err.to_string().contains("reference read failed"),
            "unexpected error: {err}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn missing_fai_marks_the_reference_unavailable() {
        let d = dir("nofai");
        let mut tsv = String::from(HEADER);
        for i in 0..3 {
            tsv.push_str(&format!("r{i}\tctg1\t100\t+\tAG\t3\t200,210,220\t15\n"));
        }
        std::fs::write(d.join("map").join("rescued.tsv"), tsv).unwrap();
        // htslib auto-builds a missing .fai when the directory is writable,
        // so the unavailable path is exercised with an unreadable reference.
        let fa = d.join("no-such-reference.fa");
        run(&d.join("map"), &d, &fa).unwrap();
        let report: serde_json::Value = serde_json::from_str(
            &std::fs::read_to_string(d.join("hyper").join("polarity_report.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(report["reference"], "unavailable");
        assert_eq!(report["sites_checked"], 0);
        assert_eq!(report["mismatches"], 0);
        let _ = std::fs::remove_dir_all(&d);
    }
}

#[cfg(test)]
mod duplicate_guard_tests {
    use super::*;
    use std::io::Write;

    const HEADER: &str = "qname\tchrom\tpos0\tstrand\tdir\tn_edits\tedit_pos0\tmapq\n";

    fn dir(tag: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("hyper_dup_{tag}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(d.join("map")).unwrap();
        d
    }

    fn fasta(d: &std::path::Path) -> std::path::PathBuf {
        let fa = d.join("ref.fa");
        let mut f = std::fs::File::create(&fa).unwrap();
        writeln!(f, ">ctg1").unwrap();
        writeln!(f, "{}", "A".repeat(300)).unwrap();
        drop(f);
        fa
    }

    /// Two canonical groups keeping the same position -> hard error, no write.
    #[test]
    fn duplicate_position_across_groups_is_a_hard_error() {
        let d = dir("conflict");
        let mut tsv = String::from(HEADER);
        // AG/+ reads edit at 200,210,220; TC/- reads (canonical AG) would
        // merge — so force a second group with TC/+ (canonical TC) over the
        // same positions: both groups keep them -> duplicate.
        for i in 0..3 {
            tsv.push_str(&format!("a{i}\tctg1\t100\t+\tAG\t3\t200,210,220\t15\n"));
        }
        for i in 0..3 {
            tsv.push_str(&format!("t{i}\tctg1\t100\t+\tTC\t3\t200,210,220\t15\n"));
        }
        std::fs::write(d.join("map").join("rescued.tsv"), tsv).unwrap();
        let fa = fasta(&d);
        let err = run(&d.join("map"), &d, &fa).err().unwrap();
        let msg = err.to_string();
        assert!(msg.contains("duplicate hyper site ctg1:200"), "{msg}");
        assert!(msg.contains("conflicting"), "{msg}");
        assert!(!d.join("hyper").join("sites.tsv").exists());
        let _ = std::fs::remove_dir_all(&d);
    }
}
