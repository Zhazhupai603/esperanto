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
    let mut by_group: std::collections::BTreeMap<(String, String), Vec<usize>> =
        std::collections::BTreeMap::new();
    for (i, r) in rows.iter().enumerate() {
        by_group
            .entry((r.chrom.clone(), r.dir.clone()))
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
