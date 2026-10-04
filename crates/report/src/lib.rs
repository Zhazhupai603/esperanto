//! esperanto-report - self-contained HTML report for a finished run
//! directory.
//!
//! Reads `<out>/sites.vcf`, `<out>/qc/qc.json` and `<out>/map/align_qc.json`
//! plus the reference FASTA (`.fai` required) and GTF, then packs per-gene
//! aggregation, genome heat maps (5 Mb / 1 Mb) and recoding (amino-acid
//! change) annotations into the embedded report template. The data is
//! inlined as a JS object literal so `report.html` renders standalone from
//! `file://`.
//!
//! Aggregation uses deterministic ordering and decimal rounding.

mod annotate;

use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};

use annotate::{Annotation, Recoder};
use anyhow::{bail, Context};
use fasta::FastaIndex;
use serde::Serialize;

mod fasta;

/// Embedded UI template; the `__ESPERANTO_DATA__` marker line is replaced
/// with the serialized data pack.
const TEMPLATE: &str = include_str!("report_template.html");
/// Marker replaced by the data pack (must occur exactly once).
const MARKER: &str = "null; // __ESPERANTO_DATA__";
/// Per-gene top-sites cap.
const GENE_SITES_CAP: usize = 60;
/// Minimum score for a site to enter a gene's site list.
const GENE_SITE_MIN_SCORE: f64 = 0.3;
/// Minimum score for the recoding check.
const RECODE_MIN_SCORE: f64 = 0.2;
/// A site is flagged "hyperedited" when a rescued read lands within this
/// many base pairs of it.
const HYPER_PAD: i64 = 500;

/// Round the exact binary value to decimal places, ties to even, and return
/// the nearest f64.
fn round_dp(x: f64, dp: usize) -> f64 {
    format!("{x:.dp$}").parse().unwrap_or(x)
}

/// One VCF site row.
struct Site {
    chrom: String,
    pos: i64,
    score: f64,
    vaf: f64,
    depth: i64,
    pass_: bool,
    gene: String,
}

/// Embedded report document.
#[derive(Serialize)]
struct DataPack {
    sample: String,
    metrics: Metrics,
    chroms: Vec<ChromLen>,
    heat5: BTreeMap<String, Vec<f64>>,
    heat1: BTreeMap<String, Vec<f64>>,
    /// Rescued (hyperedited) read counts per window, same binning as heat5/heat1.
    resc5: BTreeMap<String, Vec<u64>>,
    resc1: BTreeMap<String, Vec<u64>>,
    sites: BTreeMap<String, Vec<SiteRow>>,
    recodings: Vec<Recoding>,
    genes: Vec<GeneOut>,
}

#[derive(Serialize)]
struct Metrics {
    events: u64,
    pass: u64,
    reads_in: u64,
    map_rate: f64,
    rescued: u64,
    recodings: u64,
}

#[derive(Serialize)]
struct ChromLen {
    chrom: String,
    len: i64,
}

/// Site row serialized as `[pos, score, vaf, depth, pass, gene, hyper]`.
#[derive(Serialize)]
struct SiteRow(i64, f64, f64, i64, i64, String, bool);

#[derive(Serialize)]
struct Recoding {
    gene: String,
    chrom: String,
    pos: i64,
    change: String,
    score: f64,
    vaf: f64,
    depth: i64,
    pass: i64,
    has_hyper: bool,
}

#[derive(Serialize)]
struct GeneOut {
    gene: String,
    chrom: String,
    start: i64,
    end: i64,
    n: u64,
    pass: u64,
    maxscore: f64,
    n_rec: u64,
    sites: Vec<GeneSiteOut>,
    has_hyper: bool,
}

#[derive(Clone, Serialize)]
struct GeneSiteOut {
    pos: i64,
    score: f64,
    vaf: f64,
    depth: i64,
    pass: i64,
}

/// Per-gene accumulator (insertion order preserved for stable sorting).
struct GeneAcc {
    gene: String,
    n: u64,
    pass: u64,
    maxscore: f64,
    n_rec: u64,
    sites: Vec<GeneSiteOut>,
}

/// Parse `<out>/sites.vcf`: INFO keys ESP_SCORE / VAF / DEPTH (defaults 0),
/// FILTER column `PASS`.
fn read_sites(vcf: &Path) -> anyhow::Result<Vec<Site>> {
    let text =
        std::fs::read_to_string(vcf).with_context(|| format!("reading {}", vcf.display()))?;
    let mut sites = Vec::new();
    for (ln, line) in text.lines().enumerate() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            bail!("{} line {}: expected >= 8 columns", vcf.display(), ln + 1);
        }
        let pos: i64 = f[1]
            .parse()
            .with_context(|| format!("{} line {}: bad POS '{}'", vcf.display(), ln + 1, f[1]))?;
        let (mut score, mut vaf, mut depth) = (0.0f64, 0.0f64, 0i64);
        for part in f[7].split(';') {
            if let Some((k, v)) = part.split_once('=') {
                match k {
                    "ESP_SCORE" => {
                        score = v.parse().with_context(|| {
                            format!("{} line {}: bad ESP_SCORE '{v}'", vcf.display(), ln + 1)
                        })?
                    }
                    "VAF" => {
                        vaf = v.parse().with_context(|| {
                            format!("{} line {}: bad VAF '{v}'", vcf.display(), ln + 1)
                        })?
                    }
                    "DEPTH" => {
                        depth = v.parse().with_context(|| {
                            format!("{} line {}: bad DEPTH '{v}'", vcf.display(), ln + 1)
                        })?
                    }
                    _ => {}
                }
            }
        }
        sites.push(Site {
            chrom: f[0].to_string(),
            pos,
            score,
            vaf,
            depth,
            pass_: f[6] == "PASS",
            gene: String::new(),
        });
    }
    Ok(sites)
}

/// metrics inputs: `qc/qc.json` `summary.reads_before`, `map/align_qc.json`
/// `mapping_rate` and `rescued_collapsed`. Missing artifacts (Bam-entry
/// runs have no qc/ or map/) count as zero; an `align_qc.json` predating
/// the `rescued_collapsed` key reads as 0 (backward compatible).
fn read_metrics(out_dir: &Path) -> anyhow::Result<(u64, f64, u64)> {
    let mut reads_in = 0u64;
    if let Ok(text) = std::fs::read_to_string(out_dir.join("qc").join("qc.json")) {
        let v: serde_json::Value = serde_json::from_str(&text).context("parsing qc/qc.json")?;
        if let Some(n) = v["summary"]["reads_before"].as_u64() {
            reads_in = n;
        }
    }
    let mut map_rate = 0.0f64;
    let mut rescued = 0u64;
    if let Ok(text) = std::fs::read_to_string(out_dir.join("map").join("align_qc.json")) {
        let v: serde_json::Value =
            serde_json::from_str(&text).context("parsing map/align_qc.json")?;
        if let Some(r) = v["mapping_rate"].as_f64() {
            map_rate = r;
        }
        if let Some(n) = v["rescued_collapsed"].as_u64() {
            rescued = n;
        }
    }
    Ok((reads_in, map_rate, rescued))
}

/// Reportable primary chromosomes in FASTA index order.
fn report_chroms(fa: &FastaIndex) -> Vec<ChromLen> {
    fa.contigs()
        .iter()
        .filter(|(n, _)| n.starts_with("chr") && !n.contains('_') && n != "chrM")
        .map(|(n, l)| ChromLen {
            chrom: n.clone(),
            len: *l,
        })
        .collect()
}

/// PASS-site score mass per `win`-sized window, rounded to 2 dp.
fn mkheat(sites: &[Site], chroms: &[ChromLen], win: i64) -> BTreeMap<String, Vec<f64>> {
    let mut vals: HashMap<&str, Vec<f64>> = chroms
        .iter()
        .map(|c| (c.chrom.as_str(), vec![0.0; (c.len / win + 1) as usize]))
        .collect();
    for s in sites {
        if !s.pass_ {
            continue;
        }
        if let Some(v) = vals.get_mut(s.chrom.as_str()) {
            v[(s.pos / win) as usize] += s.score;
        }
    }
    vals.into_iter()
        .map(|(k, v)| {
            (
                k.to_string(),
                v.into_iter().map(|x| round_dp(x, 2)).collect(),
            )
        })
        .collect()
}

/// Rescued-read placements from the mapping stage, grouped by chromosome.
fn read_rescued(out_dir: &Path) -> BTreeMap<String, Vec<i64>> {
    let Ok(text) = std::fs::read_to_string(out_dir.join("map").join("rescued.bed")) else {
        return BTreeMap::new();
    };
    let mut pos: BTreeMap<String, Vec<i64>> = BTreeMap::new();
    for line in text.lines() {
        let Some((chrom, p)) = line.split_once('\t') else {
            continue;
        };
        let Ok(p) = p.parse::<i64>() else {
            continue;
        };
        pos.entry(chrom.to_owned()).or_default().push(p);
    }
    for v in pos.values_mut() {
        v.sort_unstable();
    }
    pos
}

/// Rescued-read counts per `win`-sized window from `read_rescued` positions.
fn rescued_heat(
    rescued: &BTreeMap<String, Vec<i64>>,
    chroms: &[ChromLen],
    win: i64,
) -> BTreeMap<String, Vec<u64>> {
    let mut vals: BTreeMap<String, Vec<u64>> = chroms
        .iter()
        .map(|c| (c.chrom.clone(), vec![0u64; (c.len / win + 1) as usize]))
        .collect();
    for (chrom, positions) in rescued {
        if let Some(v) = vals.get_mut(chrom) {
            for &pos in positions {
                let i = (pos / win) as usize;
                if i < v.len() {
                    v[i] += 1;
                }
            }
        }
    }
    vals
}

/// True when any rescued placement lands in `[start, end]`.
fn overlaps_rescued(
    rescued: &BTreeMap<String, Vec<i64>>,
    chrom: &str,
    start: i64,
    end: i64,
) -> bool {
    let Some(positions) = rescued.get(chrom) else {
        return false;
    };
    let i = positions.partition_point(|&p| p < start);
    i < positions.len() && positions[i] <= end
}

/// True when a rescued placement lands within `pad` bp of `pos`.
fn near_rescued(rescued: &BTreeMap<String, Vec<i64>>, chrom: &str, pos: i64, pad: i64) -> bool {
    overlaps_rescued(rescued, chrom, pos - pad, pos + pad)
}

/// Generate `<out>/<sample>.report.html`; returns its path.
pub fn generate(out_dir: &Path, fasta: &Path, gtf: &Path) -> anyhow::Result<PathBuf> {
    let mut sites = read_sites(&out_dir.join("sites.vcf"))?;
    let ann = Annotation::load(gtf)?;
    let mut fa = FastaIndex::open(fasta)?;
    let chroms = report_chroms(&fa);
    let (reads_in, map_rate, rescued) = read_metrics(out_dir)?;
    let rescued_pos = read_rescued(out_dir);

    // Gene annotation + aggregation + recoding in one VCF-order pass.
    let mut per_gene: Vec<GeneAcc> = Vec::new();
    let mut gene_idx: HashMap<String, usize> = HashMap::new();
    let mut recodings: Vec<Recoding> = Vec::new();
    let mut recoder = Recoder::new(&mut fa, &ann);
    for s in &mut sites {
        let g = ann.gene_at(&s.chrom, s.pos).map(str::to_owned);
        if let Some(ref gname) = g {
            s.gene = gname.clone();
            let gi = *gene_idx.entry(gname.clone()).or_insert_with(|| {
                per_gene.push(GeneAcc {
                    gene: gname.clone(),
                    n: 0,
                    pass: 0,
                    maxscore: 0.0,
                    n_rec: 0,
                    sites: Vec::new(),
                });
                per_gene.len() - 1
            });
            let pg = &mut per_gene[gi];
            pg.n += 1;
            pg.pass += u64::from(s.pass_);
            pg.maxscore = pg.maxscore.max(s.score);
            if s.score >= GENE_SITE_MIN_SCORE {
                pg.sites.push(GeneSiteOut {
                    pos: s.pos,
                    score: round_dp(s.score, 3),
                    vaf: round_dp(s.vaf, 3),
                    depth: s.depth,
                    pass: i64::from(s.pass_),
                });
            }
        }
        if s.score < RECODE_MIN_SCORE {
            continue;
        }
        // First valid transcript hit per site (transcripts sorted by CDS
        // start per chrom); candidate span check uses the raw VCF pos.
        let mut hit: Option<(char, char, i64, String)> = None;
        for (c0, c1, tid) in ann.transcripts(&s.chrom) {
            if *c0 > s.pos {
                break;
            }
            if *c0 <= s.pos && s.pos < *c1 {
                if let Some((ra, aa, ap)) = recoder.recode_at(tid, s.pos)? {
                    hit = Some((ra, aa, ap, tid.clone()));
                    break;
                }
            }
        }
        if let Some((ra, aa, ap, tid)) = hit {
            let tgene = recoder.gene_of(&tid);
            recodings.push(Recoding {
                gene: tgene,
                chrom: s.chrom.clone(),
                pos: s.pos,
                change: format!("{ra}{ap}{aa}"),
                score: round_dp(s.score, 3),
                vaf: round_dp(s.vaf, 3),
                depth: s.depth,
                pass: i64::from(s.pass_),
                has_hyper: near_rescued(&rescued_pos, &s.chrom, s.pos, HYPER_PAD),
            });
            if let Some(&gi) = gene_idx.get(s.gene.as_str()) {
                per_gene[gi].n_rec += 1;
            }
        }
    }

    // Sites arrays per chrom, sorted lexicographically by the full row.
    let mut per_chrom: BTreeMap<String, Vec<SiteRow>> = BTreeMap::new();
    for s in &sites {
        per_chrom.entry(s.chrom.clone()).or_default().push(SiteRow(
            s.pos,
            round_dp(s.score, 3),
            round_dp(s.vaf, 3),
            s.depth,
            i64::from(s.pass_),
            s.gene.clone(),
            near_rescued(&rescued_pos, &s.chrom, s.pos, HYPER_PAD),
        ));
    }
    for v in per_chrom.values_mut() {
        v.sort_by(cmp_site_row);
    }

    // Recodings by score descending (stable).
    recodings.sort_by(|a, b| b.score.total_cmp(&a.score));

    // Genes by unrounded maxscore descending (stable, insertion order on
    // ties), dropping names without gene coordinates.
    let mut gene_pairs: Vec<(f64, GeneOut)> = per_gene
        .iter()
        .filter_map(|pg| {
            let (gchrom, gstart, gend) = ann.gene_coord(&pg.gene)?;
            let mut top = pg.sites.clone();
            top.sort_by(|a, b| b.score.total_cmp(&a.score));
            Some((
                pg.maxscore,
                GeneOut {
                    gene: pg.gene.clone(),
                    chrom: gchrom.to_string(),
                    start: gstart,
                    end: gend,
                    n: pg.n,
                    pass: pg.pass,
                    maxscore: round_dp(pg.maxscore, 3),
                    n_rec: pg.n_rec,
                    sites: top.into_iter().take(GENE_SITES_CAP).collect(),
                    has_hyper: overlaps_rescued(&rescued_pos, gchrom, gstart, gend),
                },
            ))
        })
        .collect();
    gene_pairs.sort_by(|a, b| b.0.total_cmp(&a.0));
    let gene_out: Vec<GeneOut> = gene_pairs.into_iter().map(|(_, g)| g).collect();

    let heat5 = mkheat(&sites, &chroms, 5_000_000);
    let heat1 = mkheat(&sites, &chroms, 1_000_000);
    let resc5 = rescued_heat(&rescued_pos, &chroms, 5_000_000);
    let resc1 = rescued_heat(&rescued_pos, &chroms, 1_000_000);
    let sample = out_dir
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| "run".to_string());
    let pack = DataPack {
        sample: sample.clone(),
        metrics: Metrics {
            events: sites.len() as u64,
            pass: sites.iter().filter(|s| s.pass_).count() as u64,
            reads_in,
            map_rate,
            rescued,
            recodings: recodings.len() as u64,
        },
        chroms,
        heat5,
        heat1,
        resc5,
        resc1,
        sites: per_chrom,
        recodings,
        genes: gene_out,
    };
    let html = render(&pack)?;
    let out_path = out_dir.join(format!("{sample}.report.html"));
    std::fs::write(&out_path, &html).with_context(|| format!("writing {}", out_path.display()))?;
    Ok(out_path)
}

/// Lexicographic row comparison:
/// pos, score, vaf, depth, pass, gene. The trailing `hyper` flag is display
/// metadata and never participates in ordering.
fn cmp_site_row(a: &SiteRow, b: &SiteRow) -> std::cmp::Ordering {
    a.0.cmp(&b.0)
        .then(a.1.total_cmp(&b.1))
        .then(a.2.total_cmp(&b.2))
        .then(a.3.cmp(&b.3))
        .then(a.4.cmp(&b.4))
        .then(a.5.cmp(&b.5))
}

/// Serialize the pack and splice it into the template at the marker line.
fn render(pack: &DataPack) -> anyhow::Result<String> {
    if TEMPLATE.matches(MARKER).count() != 1 {
        bail!("report template marker missing or duplicated");
    }
    let json = serde_json::to_string(pack)?;
    // `</` guard: keep any string content from terminating the script tag.
    let json = json.replace("</", "<\\/");
    Ok(TEMPLATE.replacen(MARKER, &format!("{json};"), 1))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_dp_matches_reference_conventions() {
        // Ties to even on exactly representable halves (0.0625 and 0.1875
        // are exact in binary; 0.0615 is not - its exact value sits below
        // the tie, so it rounds down, matching the reference).
        assert_eq!(round_dp(0.0625, 3), 0.062);
        assert_eq!(round_dp(0.1875, 3), 0.188);
        assert_eq!(round_dp(0.0615, 3), 0.061);
        assert_eq!(round_dp(0.9662381, 3), 0.966);
        assert_eq!(round_dp(85.529999, 2), 85.53);
    }

    #[test]
    fn render_embeds_data_with_script_guard() {
        let pack = DataPack {
            sample: "s</script>x".into(),
            metrics: Metrics {
                events: 1,
                pass: 1,
                reads_in: 2,
                map_rate: 0.5,
                rescued: 0,
                recodings: 0,
            },
            chroms: vec![],
            heat5: BTreeMap::new(),
            heat1: BTreeMap::new(),
            resc5: BTreeMap::new(),
            resc1: BTreeMap::new(),
            sites: BTreeMap::new(),
            recodings: vec![],
            genes: vec![],
        };
        let html = render(&pack).unwrap();
        assert!(html.contains("const EMBEDDED_DATA = {"));
        assert!(!html.contains("null; // __ESPERANTO_DATA__"));
        // The literal `</script>` inside the data must be escaped.
        assert!(!html[..html.find("let reportData").unwrap()].contains("</script>x"));
        assert!(html.contains("s<\\/script>x"));
    }

    /// Same inputs, same binary: the report is byte-identical across repeated
    /// generation, and a gene present on two contigs (PAR-style multi-copy)
    /// is attributed to the last GTF row in file order — including its
    /// hyperoverlap flag.
    #[test]
    fn report_is_deterministic_with_multicopy_genes() {
        let dir = std::env::temp_dir().join(format!("esperanto-reportdet-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let out = dir.join("sample");
        std::fs::create_dir_all(out.join("qc")).unwrap();
        std::fs::create_dir_all(out.join("map")).unwrap();
        std::fs::write(
            out.join("sites.vcf"),
            "##fileformat=VCFv4.2\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chrX\t500\t.\tA\tG\t.\tPASS\tESP_SCORE=0.95;VAF=0.4;DEPTH=40\n\
             chrY\t700\t.\tA\tG\t.\tPASS\tESP_SCORE=0.8;VAF=0.3;DEPTH=30\n",
        )
        .unwrap();
        std::fs::write(
            out.join("qc").join("qc.json"),
            "{\"summary\":{\"reads_before\":4}}",
        )
        .unwrap();
        std::fs::write(
            out.join("map").join("align_qc.json"),
            "{\"mapping_rate\":1.0,\"rescued_collapsed\":2}",
        )
        .unwrap();
        // Rescued region overlaps the chrY copy only: has_hyper follows the
        // attributed copy's coordinates.
        std::fs::write(out.join("map").join("rescued.bed"), "chrY\t500\n").unwrap();

        // Two-contig reference with a hand-computed .fai.
        let seq = "ACGT".repeat(500); // 2000 bp per contig
        let mut fa_text = String::new();
        let mut fai = String::new();
        for name in ["chrX", "chrY"] {
            let offset = fa_text.len();
            fa_text.push_str(&format!(">{name}\n"));
            for i in (0..2000).step_by(60) {
                fa_text.push_str(&seq[i..(i + 60).min(2000)]);
                fa_text.push('\n');
            }
            fai.push_str(&format!("{name}\t2000\t{offset}\t60\t61\n"));
        }
        let fa = dir.join("ref.fa");
        std::fs::write(&fa, &fa_text).unwrap();
        std::fs::write(dir.join("ref.fa.fai"), &fai).unwrap();

        // PARFX on chrX first, chrY second -> the chrY copy wins. One CDS
        // transcript on chrX exercises the recoding walk.
        let gtf = dir.join("ann.gtf");
        std::fs::write(
            &gtf,
            "chrX\tsrc\tgene\t101\t900\t.\t+\t.\tgene_id \"gx\"; gene_name \"PARFX\";\n\
             chrX\tsrc\tCDS\t401\t600\t.\t+\t0\tgene_id \"gx\"; gene_name \"PARFX\"; transcript_id \"tx1\";\n\
             chrY\tsrc\tgene\t101\t900\t.\t+\t.\tgene_id \"gy\"; gene_name \"PARFX\";\n",
        )
        .unwrap();

        let mut digests = Vec::new();
        for _ in 0..3 {
            let p = generate(&out, &fa, &gtf).unwrap();
            digests.push(std::fs::read(&p).unwrap());
        }
        assert!(digests[0] == digests[1] && digests[1] == digests[2]);
        let html = String::from_utf8(digests[0].clone()).unwrap();
        assert!(
            html.contains("\"gene\":\"PARFX\",\"chrom\":\"chrY\""),
            "attribution must be the chrY copy"
        );
        assert!(
            html.contains("\"has_hyper\":true"),
            "chrY copy overlaps the rescued region"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }
}
