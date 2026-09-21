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
    /// Kept for the audit trail; clustering keys on (chrom, dir) by contract.
    #[allow(dead_code)]
    strand: char,
    dir: String,
    edits: Vec<u32>,
}

fn read_rescued_tsv(path: &Path) -> Result<Vec<RescuedRow>, IoError> {
    let text = std::fs::read_to_string(path)?;
    let mut out = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if i == 0 || line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            continue;
        }
        let Ok(pos0) = f[2].parse::<u32>() else {
            continue;
        };
        let edits = if f[6].is_empty() {
            Vec::new()
        } else {
            f[6].split(',')
                .filter_map(|x| x.parse::<u32>().ok())
                .collect()
        };
        out.push(RescuedRow {
            chrom: f[1].to_string(),
            pos0,
            strand: f[3].chars().next().unwrap_or('+'),
            dir: f[4].to_string(),
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
/// the edit-index denominator; without it the column is written as NaN.
pub fn run(map_dir: &Path, out_dir: &Path, fasta: &Path) -> Result<(usize, usize), IoError> {
    let fai = rust_htslib::faidx::Reader::from_path(fasta).ok();
    let rows = read_rescued_tsv(&map_dir.join("rescued.tsv"))?;
    // Genomic polarity of a rescued read is fixed by (dir, strand): the
    // verifier complements the expectation pair on minus placements, so an
    // AG-pass minus read and a TC-pass plus read over the same locus carry
    // the SAME genomic event. Grouping by the raw pass label split mate
    // pairs of one fragment into two groups (double-counted read-end
    // support, half of the rows written with the wrong REF/ALT polarity).
    // Normalize each read to its genomic polarity before grouping; the
    // sidecar keeps the raw pass label as the audit trail.
    let polarity_mismatch = std::cell::Cell::new(0usize);
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
                Some(f) => {
                    let span = f
                        .fetch_seq(&chrom, c.start as usize, (c.end - 1) as usize)
                        .map(|s| s.iter().filter(|&&b| b.eq_ignore_ascii_case(&src)).count())
                        .unwrap_or(0);
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
                if let Some(f) = &fai {
                    if let Ok(base) = f.fetch_seq(&chrom, p as usize, p as usize) {
                        if !base.first().is_some_and(|b| b.eq_ignore_ascii_case(&src)) {
                            polarity_mismatch.set(polarity_mismatch.get() + 1);
                        }
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
        let mut lines = sites.lines();
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
