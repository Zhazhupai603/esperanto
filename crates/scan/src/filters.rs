//! Standard evidence filters for candidate sites.
//!
//! Each rule is a visible FILTER tag in the final VCF, never a silent
//! deletion. Sources: Ramaswami et al. (Nat Methods 9:490, 2012; splice
//! proximity, read ends), REDItools/Picardi (homopolymers), Ramaswami et
//! al. (Nat Methods 10:998, 2013; simple repeats), Kleinman & Majewski
//! (RNA 18:87, 2012; intronless retrocopy blacklists).

use std::collections::BTreeMap;
use std::path::Path;

/// FILTER tag: site within this distance of an annotated exon boundary
/// (Ramaswami 2012 used 4 bp; Toung 2014, 6 bp).
pub const SPLICE_ADJ_BP: u32 = 4;

/// FILTER tag: reference homopolymer run at least this long containing the
/// site (REDItools default).
pub const HOMOPOLYMER_RUN: u32 = 5;

/// A site's filter tags (empty = none).
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct SiteFilters {
    /// Site sits in an intron within [`SPLICE_ADJ_BP`] of an exon boundary.
    pub splice_adj: bool,
    /// Reference homopolymer run >= [`HOMOPOLYMER_RUN`] at the site.
    pub homopolymer: bool,
    /// Site inside a RepeatMasker Simple_repeat / Low_complexity region.
    pub simple_repeat: bool,
    /// Site inside a generated intronless-retrocopy interval.
    pub retrocopy: bool,
    /// Site inside a generated multicopy-exon window.
    pub multicopy: bool,
}

impl SiteFilters {
    /// The FILTER tags, in fixed order.
    pub fn tags(&self) -> Vec<&'static str> {
        let mut out = Vec::new();
        if self.splice_adj {
            out.push("SPLICE_ADJ");
        }
        if self.homopolymer {
            out.push("HOMOPOLYMER");
        }
        if self.simple_repeat {
            out.push("SIMPLE_REPEAT");
        }
        if self.retrocopy {
            out.push("RETROCOPY");
        }
        if self.multicopy {
            out.push("MULTICOPY");
        }
        out
    }

    pub fn is_clean(&self) -> bool {
        self.tags().is_empty()
    }
}

/// Sorted interval set over one contig, binary-searched.
struct Intervals {
    starts: Vec<u32>,
    ends: Vec<u32>,
}

impl Intervals {
    fn from_bed(path: &Path, label: &str) -> Option<(BTreeMap<String, Intervals>, usize)> {
        let text = std::fs::read_to_string(path).ok()?;
        let mut by_contig: BTreeMap<String, Vec<(u32, u32)>> = BTreeMap::new();
        let mut n = 0usize;
        for line in text.lines() {
            if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 3 {
                continue;
            }
            let (Ok(s), Ok(e)) = (f[1].parse::<u32>(), f[2].parse::<u32>()) else {
                continue;
            };
            // The label filter keeps the bed generic (RepeatMasker style has
            // the family in column 4; generated blacklists do not).
            if let Some(fam) = f.get(3) {
                if label == "Simple_repeat"
                    && !fam.contains("Simple_repeat")
                    && !fam.contains("Low_complexity")
                {
                    continue;
                }
            }
            by_contig.entry(f[0].to_string()).or_default().push((s, e));
            n += 1;
        }
        let out = by_contig
            .into_iter()
            .map(|(c, mut ivs)| {
                ivs.sort_unstable();
                let mut merged: Vec<(u32, u32)> = Vec::with_capacity(ivs.len());
                for (s, e) in ivs {
                    match merged.last_mut() {
                        Some((_, le)) if s <= *le => *le = (*le).max(e),
                        _ => merged.push((s, e)),
                    }
                }
                let (starts, ends): (Vec<_>, Vec<_>) = merged.into_iter().unzip();
                (c, Intervals { starts, ends })
            })
            .collect();
        Some((out, n))
    }

    fn contains(&self, pos: u32) -> bool {
        let i = self.starts.partition_point(|&s| s <= pos);
        if i == 0 {
            return false;
        }
        self.ends[i - 1] > pos
    }
}

/// Exon boundary positions (0-based, both ends) per contig, from a GTF.
struct ExonBoundaries {
    by_contig: BTreeMap<String, Vec<u32>>,
}

impl ExonBoundaries {
    fn from_gtf(path: &Path) -> Option<(Self, usize)> {
        let text = std::fs::read_to_string(path).ok()?;
        let mut by_contig: BTreeMap<String, Vec<u32>> = BTreeMap::new();
        let mut n = 0usize;
        for line in text.lines() {
            if line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 8 || f[2] != "exon" {
                continue;
            }
            let (Ok(s), Ok(e)) = (f[3].parse::<u32>(), f[4].parse::<u32>()) else {
                continue;
            };
            let v = by_contig.entry(f[0].to_string()).or_default();
            v.push(s - 1); // 1-based inclusive -> 0-based boundary
            v.push(e); // exon end boundary (1-based inclusive == 0-based exclusive end)
            n += 1;
        }
        for v in by_contig.values_mut() {
            v.sort_unstable();
            v.dedup();
        }
        Some((ExonBoundaries { by_contig }, n))
    }

    fn near(&self, contig: &str, pos: u32, dist: u32) -> bool {
        let Some(v) = self.by_contig.get(contig) else {
            return false;
        };
        // Closest boundary at or after pos - dist.
        let lo = pos.saturating_sub(dist);
        let i = v.partition_point(|&b| b < lo);
        i < v.len() && v[i] <= pos + dist
    }
}

/// Site-filter oracle: loads the optional inputs once, answers per site.
pub struct FilterSet {
    exon: Option<ExonBoundaries>,
    rmsk: Option<BTreeMap<String, Intervals>>,
    retrocopy: Option<BTreeMap<String, Intervals>>,
    multicopy: Option<BTreeMap<String, Intervals>>,
    /// `ESP_SCAN_FILTERS` allow-list (None = all enabled; empty = all off).
    only: Option<Vec<String>>,
}

impl FilterSet {
    /// Load from the refs directory conventions: `genes.gtf` (any `*.gtf`),
    /// `*.rmsk.bed` / `*.rmsk.bed.gz` (RepeatMasker), `*.retrocopy.bed`,
    /// `*.multicopy.bed`. Every input is optional; a missing source simply
    /// disables its tag.
    pub fn load(refs_dir: &Path) -> FilterSet {
        let pick = |exts: &[&str]| -> Option<std::path::PathBuf> {
            let mut hits: Vec<std::path::PathBuf> = std::fs::read_dir(refs_dir)
                .ok()?
                .flatten()
                .map(|e| e.path())
                .filter(|p| {
                    let name = match p.file_name() {
                        Some(n) => n.to_string_lossy().into_owned(),
                        None => return false,
                    };
                    exts.iter().any(|x| name.ends_with(x))
                        && !name.ends_with(".tidx")
                        && !name.ends_with(".paidx")
                })
                .collect();
            hits.sort();
            hits.into_iter().next()
        };
        let exon = pick(&[".gtf"]).and_then(|p| ExonBoundaries::from_gtf(&p).map(|(e, _)| e));
        let rmsk = pick(&[".rmsk.bed", "rmsk.bed"])
            .and_then(|p| Intervals::from_bed(&p, "Simple_repeat").map(|(m, _)| m));
        let retrocopy =
            pick(&[".retrocopy.bed"]).and_then(|p| Intervals::from_bed(&p, "").map(|(m, _)| m));
        let multicopy =
            pick(&[".multicopy.bed"]).and_then(|p| Intervals::from_bed(&p, "").map(|(m, _)| m));
        // A missing input must be visible: silent filter loss is how
        // "disabled" and "enabled but empty" become indistinguishable.
        if exon.is_none() {
            eprintln!("[scan] filter SPLICE_ADJ disabled: no *.gtf in {}", refs_dir.display());
        }
        if rmsk.is_none() {
            eprintln!("[scan] filter SIMPLE_REPEAT disabled: no *.rmsk.bed in {}", refs_dir.display());
        }
        if retrocopy.is_none() {
            eprintln!("[scan] filter RETROCOPY disabled: no *.retrocopy.bed in {}", refs_dir.display());
        }
        if multicopy.is_none() {
            eprintln!("[scan] filter MULTICOPY disabled: no *.multicopy.bed in {}", refs_dir.display());
        }
        // Per-filter cost measurement (WP4 review): ESP_SCAN_FILTERS
        // restrict which tags fire ("none" = all off, comma list otherwise).
        let only = std::env::var("ESP_SCAN_FILTERS").ok().map(|v| {
            if v == "none" {
                Vec::new()
            } else {
                v.split(',').map(|x| x.trim().to_string()).collect()
            }
        });
        FilterSet {
            exon,
            rmsk,
            retrocopy,
            multicopy,
            only,
        }
    }

    /// The empty set (no inputs, no tags).
    pub fn empty() -> FilterSet {
        FilterSet {
            exon: None,
            rmsk: None,
            retrocopy: None,
            multicopy: None,
            only: None,
        }
    }

    /// Filter tags for one site. `ref_context` is the reference slice
    /// centered on the site (site at index `center`), used for the
    /// homopolymer run test.
    pub fn site(&self, contig: &str, pos0: u32, ref_context: &[u8], center: usize) -> SiteFilters {
        let mut f = SiteFilters::default();
        if let Some(exon) = &self.exon {
            f.splice_adj = exon.near(contig, pos0, SPLICE_ADJ_BP);
        }
        f.homopolymer = homopolymer_at(ref_context, center);
        if let Some(r) = &self.rmsk {
            f.simple_repeat = r.get(contig).is_some_and(|iv| iv.contains(pos0));
        }
        if let Some(r) = &self.retrocopy {
            f.retrocopy = r.get(contig).is_some_and(|iv| iv.contains(pos0));
        }
        if let Some(m) = &self.multicopy {
            f.multicopy = m.get(contig).is_some_and(|iv| iv.contains(pos0));
        }
        if let Some(only) = &self.only {
            let on = |t: &str| only.iter().any(|x| x == t);
            f.splice_adj &= on("SPLICE_ADJ");
            f.homopolymer &= on("HOMOPOLYMER");
            f.simple_repeat &= on("SIMPLE_REPEAT");
            f.retrocopy &= on("RETROCOPY");
            f.multicopy &= on("MULTICOPY");
        }
        f
    }
}

/// Reference homopolymer run >= [`HOMOPOLYMER_RUN`] through the site base.
fn homopolymer_at(ctx: &[u8], center: usize) -> bool {
    let Some(&b) = ctx.get(center) else {
        return false;
    };
    if !b.is_ascii_alphabetic() {
        return false;
    }
    let same = |i: usize| ctx.get(i).is_some_and(|&x| x.eq_ignore_ascii_case(&b));
    let mut run = 1u32;
    let mut i = center;
    while i > 0 && same(i - 1) {
        run += 1;
        i -= 1;
    }
    let mut i = center + 1;
    while same(i) {
        run += 1;
        i += 1;
    }
    run >= HOMOPOLYMER_RUN
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn homopolymer_runs() {
        let ctx = b"AAAAAGC".to_vec();
        assert!(homopolymer_at(&ctx, 1));
        assert!(homopolymer_at(&ctx, 4), "the fifth base of the run");
        assert!(!homopolymer_at(&ctx, 5));
        assert!(!homopolymer_at(&ctx, 6));
        let four = b"AAAAGC".to_vec();
        assert!(!homopolymer_at(&four, 1));
        let long = b"GGGGG".to_vec();
        assert!(homopolymer_at(&long, 2));
        let four = b"AGGGGC".to_vec();
        assert!(!homopolymer_at(&four, 2));
    }

    #[test]
    fn boundaries_near() {
        let mut by = BTreeMap::new();
        by.insert("chr1".to_string(), vec![100, 200]);
        let e = ExonBoundaries { by_contig: by };
        assert!(e.near("chr1", 103, 4));
        assert!(e.near("chr1", 96, 4));
        assert!(!e.near("chr1", 105, 4));
        assert!(!e.near("chr1", 94, 4));
        assert!(!e.near("chr2", 100, 4));
    }

    #[test]
    fn tag_order() {
        let f = SiteFilters {
            multicopy: true,
            splice_adj: true,
            ..Default::default()
        };
        assert_eq!(f.tags(), vec!["SPLICE_ADJ", "MULTICOPY"]);
        assert!(!f.is_clean());
        assert!(SiteFilters::default().is_clean());
    }
}
