//! Native reference-aware molecular evidence scoring.
//!
//! This module implements the frozen 27-feature LightGBM contract used by the
//! unified evidence model. It accepts a model directory containing `model.txt`,
//! `freeze.json`, and `alu.bed`, and computes scores directly from
//! a coordinate-sorted indexed BAM.

use anyhow::{bail, ensure, Context, Result};
use rayon::prelude::*;
use rust_htslib::bam::{record::Cigar, IndexedReader, Read as _};
use std::collections::HashMap;
use std::path::Path;

pub const FEATURE_COUNT: usize = 27;
pub const MIN_MAPQ: u8 = 20;
pub const MIN_BQ: u8 = 13;
pub const RADIUS: i64 = 200;

const BASES: [u8; 4] = *b"ACGT";
const EXCLUDE_FLAGS: u16 = 0x4 | 0x100 | 0x800;

#[derive(Clone, Debug)]
struct Observation {
    base: u8,
    bq: u8,
    mapq: f64,
    end_fraction: f64,
    orientation: i8,
    softclip_fraction: f64,
    gap_nearness: f64,
}

type Window = HashMap<i64, HashMap<Vec<u8>, Observation>>;

#[derive(Clone, Debug)]
struct IntervalIndex {
    starts: Vec<i64>,
    ends: Vec<i64>,
}

#[derive(Clone, Debug, Default)]
pub struct AluIndex {
    intervals: HashMap<String, IntervalIndex>,
}

impl AluIndex {
    pub fn from_bed(path: &Path) -> Result<Self> {
        let mut raw: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
        for (line_no, line) in std::fs::read_to_string(path)?.lines().enumerate() {
            if line.trim().is_empty() || line.starts_with('#') || line.starts_with("track") {
                continue;
            }
            let fields: Vec<_> = line.split_whitespace().collect();
            ensure!(fields.len() >= 3, "invalid Alu BED row {}", line_no + 1);
            raw.entry(fields[0].to_owned()).or_default().push((
                fields[1].parse().with_context(|| format!("Alu start at line {}", line_no + 1))?,
                fields[2].parse().with_context(|| format!("Alu end at line {}", line_no + 1))?,
            ));
        }
        let intervals = raw
            .into_iter()
            .map(|(chrom, mut values)| {
                values.sort_unstable();
                let starts = values.iter().map(|x| x.0).collect();
                let mut ends = Vec::with_capacity(values.len());
                for (_, end) in values {
                    ends.push(ends.last().copied().unwrap_or(0).max(end));
                }
                (chrom, IntervalIndex { starts, ends })
            })
            .collect();
        Ok(Self { intervals })
    }

    pub fn contains(&self, chrom: &str, pos0: i64) -> bool {
        let Some(index) = self.intervals.get(chrom) else { return false };
        let n = index.starts.partition_point(|&start| start <= pos0);
        n > 0 && index.ends[n - 1] > pos0
    }
}

#[derive(Clone, Debug)]
struct Tree {
    feature: Vec<usize>,
    threshold: Vec<f64>,
    decision: Vec<u8>,
    left: Vec<i32>,
    right: Vec<i32>,
    leaf: Vec<f64>,
}

impl Tree {
    fn score(&self, features: &[f32; FEATURE_COUNT]) -> f64 {
        let mut node = 0_i32;
        while node >= 0 {
            let i = node as usize;
            let value = features[self.feature[i]];
            // LightGBM encodes missing-value handling in decision_type. In
            // addition to NaN, type 1 treats an exact zero as missing (the
            // serialized values are 8/10 rather than plain 2).
            let missing_type = self.decision[i] >> 2;
            let missing = value.is_nan() || (missing_type == 1 && value == 0.0);
            let goes_left = if missing {
                self.decision[i] & 2 != 0
            } else {
                f64::from(value) <= self.threshold[i]
            };
            node = if goes_left { self.left[i] } else { self.right[i] };
        }
        self.leaf[(-node - 1) as usize]
    }
}

#[derive(Clone, Debug)]
pub struct Model {
    trees: Vec<Tree>,
}

#[derive(Clone, Debug)]
struct ObservationSite {
    chrom: String,
    pos: i64,
    ref_base: u8,
    counts: [u32; 4],
    depth: u32,
}

fn parse_values<T: std::str::FromStr>(block: &str, key: &str) -> Result<Vec<T>>
where
    T::Err: std::fmt::Debug,
{
    let prefix = format!("{key}=");
    block
        .lines()
        .find_map(|line| line.strip_prefix(&prefix))
        .with_context(|| format!("missing LightGBM field {key}"))?
        .split_whitespace()
        .map(|value| value.parse().map_err(|e| anyhow::anyhow!("invalid {key}: {e:?}")))
        .collect()
}

impl Model {
    pub fn load(model_dir: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(model_dir.join("model.txt"))?;
        let mut trees = Vec::new();
        for block in text.split("\nTree=").skip(1) {
            if !block.contains("\nsplit_feature=") {
                break;
            }
            let tree = Tree {
                feature: parse_values(block, "split_feature")?,
                threshold: parse_values(block, "threshold")?,
                decision: parse_values(block, "decision_type")?,
                left: parse_values(block, "left_child")?,
                right: parse_values(block, "right_child")?,
                leaf: parse_values(block, "leaf_value")?,
            };
            ensure!(tree.feature.len() == tree.threshold.len());
            ensure!(tree.feature.len() == tree.decision.len());
            ensure!(tree.feature.len() == tree.left.len());
            ensure!(tree.feature.len() == tree.right.len());
            ensure!(tree.leaf.len() == tree.feature.len() + 1);
            ensure!(tree.feature.iter().all(|&x| x < FEATURE_COUNT));
            trees.push(tree);
        }
        ensure!(trees.len() == 400, "expected 400 LightGBM trees, found {}", trees.len());
        Ok(Self { trees })
    }

    pub fn score(&self, features: &[f32; FEATURE_COUNT]) -> f64 {
        let raw: f64 = self.trees.iter().map(|tree| tree.score(features)).sum();
        1.0 / (1.0 + (-raw).exp())
    }
}

fn template_name(qname: &[u8]) -> Vec<u8> {
    if qname.len() >= 2 && qname[qname.len() - 2] == b'/' && matches!(qname[qname.len() - 1], b'1' | b'2') {
        qname[..qname.len() - 2].to_vec()
    } else {
        qname.to_vec()
    }
}

fn base_index(base: u8) -> Option<usize> {
    BASES.iter().position(|&value| value == base)
}

fn merge_observation(old: &mut Observation, new: Observation) {
    if new.bq > old.bq {
        *old = new;
    } else if new.bq == old.bq && old.base != new.base {
        old.base = b'N';
        old.mapq = if old.mapq.is_finite() { old.mapq.max(new.mapq) } else { new.mapq };
        old.end_fraction = old.end_fraction.max(new.end_fraction);
        old.orientation = if old.orientation == new.orientation { old.orientation } else { -1 };
        old.softclip_fraction = old.softclip_fraction.min(new.softclip_fraction);
        old.gap_nearness = old.gap_nearness.min(new.gap_nearness);
    }
}

fn read_observations(rec: &rust_htslib::bam::Record, left: i64, right: i64) -> Vec<(i64, Observation)> {
    if rec.flags() & EXCLUDE_FLAGS != 0 || rec.mapq() < MIN_MAPQ || rec.pos() < 0 {
        return Vec::new();
    }
    let seq = rec.seq().as_bytes();
    let qual = rec.qual();
    if seq.is_empty() { return Vec::new(); }
    let length = seq.len();
    let soft = rec.cigar().iter().filter_map(|c| match c { Cigar::SoftClip(n) => Some(*n as f64), _ => None }).sum::<f64>() / length.max(1) as f64;
    let mut gaps = Vec::new();
    let mut q = 0_i64;
    for op in rec.cigar().iter() {
        let n = op.len() as i64;
        match op {
            Cigar::Ins(_) => {
                gaps.push(q as f64 - 0.5);
                gaps.push((q + n) as f64 - 0.5);
                q += n;
            }
            Cigar::Del(_) | Cigar::RefSkip(_) => {
                gaps.push(q as f64 - 0.5);
            }
            Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) | Cigar::SoftClip(_) => q += n,
            _ => {}
        }
    }
    let orientation = if rec.is_reverse() ^ rec.is_last_in_template() { 1 } else { 0 };
    let mapq = if rec.mapq() == 255 { f64::NAN } else { f64::from(rec.mapq()) };
    let mut out = Vec::new();
    let mut qpos = 0_i64;
    let mut rpos = rec.pos();
    for op in rec.cigar().iter() {
        let n = op.len() as i64;
        match op {
            Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => {
                let lo = rpos.max(left);
                let hi = (rpos + n).min(right);
                for rp in lo..hi {
                    let qp = (qpos + rp - rpos) as usize;
                    let bq = qual.get(qp).copied().unwrap_or(255);
                    if bq == 255 || bq < MIN_BQ { continue; }
                    let Some(&raw_base) = seq.get(qp) else { continue; };
                    let base = raw_base.to_ascii_uppercase();
                    if base_index(base).is_none() { continue; }
                    let end_fraction = (qp.min(length - 1 - qp) as f64) / ((length - 1).max(1) as f64 / 2.0);
                    let gap = gaps.iter().map(|g| ((qp as f64) - g).abs()).fold(f64::INFINITY, f64::min);
                    out.push((rp, Observation { base, bq, mapq, end_fraction, orientation, softclip_fraction: soft, gap_nearness: if gap.is_finite() { 1.0 / (1.0 + gap) } else { 0.0 } }));
                }
                qpos += n; rpos += n;
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) => qpos += n,
            Cigar::Del(_) | Cigar::RefSkip(_) => rpos += n,
            _ => {}
        }
    }
    out
}

fn read_window(bam: &mut IndexedReader, chrom: &str, left: i64, right: i64) -> Result<Window> {
    let tid = bam.header().tid(chrom.as_bytes()).context("contig absent from BAM")?;
    bam.fetch((tid, left, right))?;
    let mut templates: HashMap<Vec<u8>, HashMap<i64, Observation>> = HashMap::new();
    for record in bam.records() {
        let record = record?;
        let mut name = record.aux(b"RG").ok().and_then(|aux| match aux {
            rust_htslib::bam::record::Aux::String(value) => Some(value.as_bytes().to_vec()),
            _ => None,
        }).unwrap_or_default();
        name.push(0);
        name.extend_from_slice(&template_name(record.qname()));
        let observations = templates.entry(name).or_default();
        for (pos, observation) in read_observations(&record, left, right) {
            if let Some(old) = observations.get_mut(&pos) { merge_observation(old, observation); } else { observations.insert(pos, observation); }
        }
    }
    let mut positions: Window = HashMap::new();
    for (name, observations) in templates {
        for (pos, observation) in observations {
            if base_index(observation.base).is_some() {
                positions.entry(pos).or_default().insert(name.clone(), observation);
            }
        }
    }
    Ok(positions)
}

fn mean(values: impl Iterator<Item = f64>) -> f64 {
    let values: Vec<_> = values.filter(|value| value.is_finite()).collect();
    if values.is_empty() { f64::NAN } else { values.iter().sum::<f64>() / values.len() as f64 }
}

fn log_likelihood(observations: &[&Observation], ref_base: u8, alt: u8, p: f64) -> f64 {
    observations.iter().map(|obs| {
        let error = 10_f64.powf(-(obs.bq as f64) / 10.0);
        let (a, b) = if obs.base == ref_base { (1.0 - error, -(1.0 - 4.0 * error / 3.0)) }
            else if obs.base == alt { (error / 3.0, 1.0 - 4.0 * error / 3.0) }
            else { (error / 3.0, 0.0) };
        (a + b * p).max(1e-300).ln()
    }).sum()
}

fn likelihood_derivative(observations: &[&Observation], ref_base: u8, alt: u8, p: f64) -> f64 {
    observations.iter().map(|obs| {
        let error = 10_f64.powf(-(obs.bq as f64) / 10.0);
        let (a, b) = if obs.base == ref_base {
            (1.0 - error, -(1.0 - 4.0 * error / 3.0))
        } else if obs.base == alt {
            (error / 3.0, 1.0 - 4.0 * error / 3.0)
        } else {
            (error / 3.0, 0.0)
        };
        b / (a + b * p).max(1e-300)
    }).sum()
}

fn likelihood_curvature(observations: &[&Observation], ref_base: u8, alt: u8, p: f64) -> f64 {
    observations.iter().map(|obs| {
        let error = 10_f64.powf(-(obs.bq as f64) / 10.0);
        let (a, b) = if obs.base == ref_base {
            (1.0 - error, -(1.0 - 4.0 * error / 3.0))
        } else if obs.base == alt {
            (error / 3.0, 1.0 - 4.0 * error / 3.0)
        } else {
            (error / 3.0, 0.0)
        };
        let denominator = (a + b * p).max(1e-300);
        (b / denominator).powi(2)
    }).sum()
}

/// 32-point Gauss-Legendre nodes and weights on [0, 1]. Used for every
/// quadrature in this module so that the cost of one site scales linearly
/// with its observation count, regardless of coverage depth.
#[allow(clippy::excessive_precision)]
const GL32_X: [f64; 32] = [0.0013680690752592151,0.007194244227365809,0.017618872206246805,0.03254696203113017,0.05183942211697395,0.07531619313371502,0.10275810201602881,0.13390894062985514,0.16847786653489238,0.20614212137961885,0.24655004553388532,0.28932436193468236,0.33406569885893617,0.38035631887393145,0.42776401920860174,0.4758461671561308,0.5241538328438692,0.5722359807913983,0.6196436811260685,0.6659343011410639,0.7106756380653176,0.7534499544661146,0.7938578786203812,0.8315221334651076,0.8660910593701449,0.8972418979839711,0.924683806866285,0.948160577883026,0.9674530379688698,0.9823811277937532,0.9928057557726342,0.9986319309247408];
#[allow(clippy::excessive_precision)]
const GL32_W: [f64; 32] = [0.003509305004735253,0.008137197365452872,0.01269603265463101,0.01713693145651088,0.02141794901111342,0.02549902963118805,0.02934204673926778,0.03291111138818084,0.03617289705442417,0.03909694789353511,0.04165596211347335,0.04382604650220189,0.04558693934788189,0.04692219954040226,0.04781936003963735,0.04827004425736383,0.04827004425736383,0.04781936003963735,0.04692219954040226,0.04558693934788189,0.04382604650220189,0.04165596211347335,0.03909694789353511,0.03617289705442417,0.03291111138818084,0.02934204673926778,0.02549902986561342,0.02141794901111342,0.01713693145651088,0.01269603265463101,0.008137197365452872,0.003509305004735253];

/// Fixed-node quadrature of `f` over [a, b].
fn integrate<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> f64 {
    if b <= a { return 0.0; }
    let width = b - a;
    width * GL32_X.iter().zip(GL32_W.iter()).map(|(&x, &w)| w * f(a + width * x)).sum::<f64>()
}

fn likelihood_features(observations: &[&Observation], ref_base: u8, alt: u8) -> (f64, f64) {
    if observations.is_empty() { return (0.0, 0.0); }
    let log_zero = log_likelihood(observations, ref_base, alt, 0.0);
    let log_one = log_likelihood(observations, ref_base, alt, 1.0);
    if observations.len() <= 63 {
        let log_uniform = GL32_X.iter().zip(GL32_W.iter()).map(|(&x, w)| (log_likelihood(observations, ref_base, alt, x) + w.ln()).exp()).sum::<f64>().ln();
        return (log_uniform - log_zero, log_uniform - log_one);
    }
    let derivative_zero = likelihood_derivative(observations, ref_base, alt, 0.0);
    let derivative_one = likelihood_derivative(observations, ref_base, alt, 1.0);
    let mode = if derivative_zero <= 0.0 {
        0.0
    } else if derivative_one >= 0.0 {
        1.0
    } else {
        let mut lo = 0.0;
        let mut hi = 1.0;
        for _ in 0..100 {
            let mid = (lo + hi) / 2.0;
            if likelihood_derivative(observations, ref_base, alt, mid) > 0.0 { lo = mid; } else { hi = mid; }
        }
        (lo + hi) / 2.0
    };
    let peak = log_likelihood(observations, ref_base, alt, mode);
    let scale = if mode == 0.0 || mode == 1.0 {
        1.0 / derivative_zero.abs().max(derivative_one.abs()).max(1.0)
    } else {
        1.0 / likelihood_curvature(observations, ref_base, alt, mode).max(1.0).sqrt()
    };
    let mut points = vec![0.0, 1.0];
    for multiple in [1.0, 4.0, 16.0, 64.0] {
        for point in [mode - multiple * scale, mode + multiple * scale] {
            if point > 0.0 && point < 1.0 { points.push(point); }
        }
    }
    points.sort_by(f64::total_cmp);
    points.dedup_by(|a, b| a == b);
    let integral = points.windows(2).map(|window| {
        integrate(|p| (log_likelihood(observations, ref_base, alt, p) - peak).exp(), window[0], window[1])
    }).sum::<f64>();
    (peak + integral.ln() - log_zero, peak + integral.ln() - log_one)
}

fn feature_vector(site: &(String, i64), ref_seq: &[u8], reference_start: i64, window: &Window, alu: &AluIndex) -> Result<[f32; FEATURE_COUNT]> {
    let (chrom, pos1) = site;
    let pos = *pos1 - 1;
    let local = usize::try_from(pos - reference_start)?;
    let ref_base = ref_seq[local].to_ascii_uppercase();
    let alt = match ref_base { b'A' => b'G', b'T' => b'C', _ => bail!("unsupported reference base at {chrom}:{pos1}") };
    let obs: Vec<&Observation> = window.get(&pos).into_iter().flat_map(|x| x.values()).collect();
    let alt_obs: Vec<_> = obs.iter().copied().filter(|o| o.base == alt).collect();
    let ref_obs: Vec<_> = obs.iter().copied().filter(|o| o.base == ref_base).collect();
    let depth = obs.len() as f64;
    let alt_count = alt_obs.len() as f64;
    let other = depth - alt_count - ref_obs.len() as f64;
    let (bf_error, bf_homo) = likelihood_features(&obs, ref_base, alt);
    let strand = if !alt_obs.is_empty() && !ref_obs.is_empty() {
        let count = |items: &[&Observation], direction: i8| items.iter().filter(|o| o.orientation == direction).count() as f64;
        Some(((count(&alt_obs, 0) + 0.5).ln() - (count(&alt_obs, 1) + 0.5).ln() - (count(&ref_obs, 0) + 0.5).ln() + (count(&ref_obs, 1) + 0.5).ln()).abs())
    } else { None };
    let mut weighted_opportunities = 0.0; let mut weighted_depth = 0.0; let mut weighted_alt = 0.0; let mut weighted_alt_sites = 0.0; let mut all_weighted_depth = 0.0; let mut all_weighted_other = 0.0; let mut template_contributions: HashMap<Vec<u8>, u32> = HashMap::new();
    for p in (pos - RADIUS).max(reference_start)..=(pos + RADIUS).min(reference_start + ref_seq.len() as i64 - 1) {
        if p == pos { continue; }
        let Some(local_obs) = window.get(&p) else { continue; };
        let lr = ref_seq[(p-reference_start) as usize].to_ascii_uppercase(); if base_index(lr).is_none() { continue; }
        let n = local_obs.len() as f64; let w = (-(p-pos).abs() as f64 / RADIUS as f64).exp();
        let local_ref = local_obs.values().filter(|o| o.base == lr).count() as f64;
        let local_alt = if lr == ref_base { local_obs.values().filter(|o| o.base == alt).count() as f64 } else { 0.0 };
        let local_other = n - local_ref - local_alt; all_weighted_depth += w*n; all_weighted_other += w*local_other;
        if lr != ref_base { continue; }
        weighted_opportunities += w; weighted_depth += w*n; weighted_alt += w*local_alt;
        if local_alt > 0.0 { weighted_alt_sites += w; for (name, o) in local_obs { if o.base == alt { *template_contributions.entry(name.clone()).or_default() += 1; } } }
    }
    let alt_observations: f64 = template_contributions.values().map(|&n| f64::from(n)).sum();
    let reuse = if alt_observations > 1.0 { template_contributions.values().map(|&n| { let n = f64::from(n); n * (n - 1.0) }).sum::<f64>() / (alt_observations * (alt_observations - 1.0)) } else { 0.0 };
    let rate = |n: f64, d: f64| if d > 0.0 { n/d } else { f64::NAN };
    let previous = if pos > reference_start { ref_seq[local-1] } else { b'N' }; let following = if local + 1 < ref_seq.len() { ref_seq[local+1] } else { b'N' };
    let (previous, following) = if ref_base == b'T' || ref_base == b'G' { (complement(following), complement(previous)) } else { (previous, following) };
    let mut values = [f64::NAN; FEATURE_COUNT];
    values[..20].copy_from_slice(&[depth.max(0.0).ln_1p(), alt_count.ln_1p(), rate(alt_count, depth), rate(other, depth), bf_error, bf_homo, mean(alt_obs.iter().map(|o| f64::from(o.bq))), mean(ref_obs.iter().map(|o| f64::from(o.bq))), mean(alt_obs.iter().map(|o| o.mapq)), mean(alt_obs.iter().map(|o| o.end_fraction)), strand.unwrap_or(f64::NAN), mean(alt_obs.iter().map(|o| o.softclip_fraction)), mean(alt_obs.iter().map(|o| o.gap_nearness)), rate(weighted_alt, weighted_depth), rate(weighted_alt_sites, weighted_opportunities), rate(all_weighted_other, all_weighted_depth), (template_contributions.len() as f64).ln_1p(), reuse, weighted_opportunities.ln_1p(), weighted_depth.ln_1p()]);
    for (i, base) in [previous, following].into_iter().enumerate() { for (j, target) in (*b"ACG").into_iter().enumerate() { values[20 + i*3 + j] = if base == target { 1.0 } else if base_index(base).is_some() { 0.0 } else { f64::NAN }; } }
    values[26] = if alu.contains(chrom, pos) { 1.0 } else { 0.0 };
    Ok(values.map(|value| value as f32))
}

fn complement(base: u8) -> u8 { match base { b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', _ => b'N' } }

pub fn score_sites(bam_path: &Path, fasta_path: &Path, model_dir: &Path, alu: &AluIndex, sites: &[(String, i64)], threads: usize) -> Result<Vec<f64>> {
    let model = Model::load(model_dir)?;
    let mut order: Vec<usize> = (0..sites.len()).collect(); order.sort_by(|&a,&b| sites[a].cmp(&sites[b]));
    let mut groups = Vec::new();
    let mut start = 0;
    while start < order.len() {
        let first = &sites[order[start]];
        let mut end = start + 1;
        let mut right = first.1 - 1 + RADIUS + 1;
        while end < order.len() {
            let next = &sites[order[end]];
            let left = next.1 - 1 - RADIUS;
            if next.0 != first.0 || left > right { break; }
            right = right.max(next.1 - 1 + RADIUS + 1);
            end += 1;
        }
        groups.push((start, end, first.0.clone(), first.1 - 1 - RADIUS));
        start = end;
    }
    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads.max(1)).build()?;
    let scores = pool.install(|| groups.par_iter().map(|&(lo, hi, ref chrom, left)| -> Result<Vec<(usize,f64)>> {
        let right = sites[order[hi - 1]].1 - 1 + RADIUS + 1;
        let mut bam = IndexedReader::from_path(bam_path)?; let fasta = rust_htslib::faidx::Reader::from_path(fasta_path)?;
        let ref_start = left.max(0); let reference = fasta.fetch_seq(chrom, ref_start as usize, (right - 1) as usize)?;
        let window = read_window(&mut bam, chrom, ref_start, right)?;
        let mut result = Vec::with_capacity(hi - lo);
        for &i in &order[lo..hi] { let features = feature_vector(&sites[i], &reference, ref_start, &window, alu)?; result.push((i, model.score(&features))); }
        Ok(result)
    }).collect::<Vec<_>>());
    let mut output = vec![0.0; sites.len()]; for result in scores { for (i, score) in result? { output[i] = score; } } Ok(output)
}

/// Feature vector or per-site failure, as produced by [`feature_matrix`].
pub type FeatureRow = std::result::Result<[f32; FEATURE_COUNT], String>;

/// Export raw 27-feature vectors for arbitrary sites, mirroring the grouping
/// and window logic of `score_sites`. Per-site failures (e.g. unsupported
/// reference base) are captured as `Err(message)` instead of aborting.
///
/// Groups are additionally capped in size/span for rayon load balance; this
/// does not change any feature value (a site's window depends only on its own
/// +-RADIUS neighbourhood).
pub fn feature_matrix(bam_path: &Path, fasta_path: &Path, alu: &AluIndex, sites: &[(String, i64)], threads: usize) -> Result<Vec<FeatureRow>> {
    let mut order: Vec<usize> = (0..sites.len()).collect(); order.sort_by(|&a,&b| sites[a].cmp(&sites[b]));
    const MAX_GROUP_SITES: usize = 512;
    const MAX_GROUP_SPAN: i64 = 400_000;
    let mut groups = Vec::new();
    let mut start = 0;
    while start < order.len() {
        let first = &sites[order[start]];
        let mut end = start + 1;
        let mut right = first.1 - 1 + RADIUS + 1;
        while end < order.len() {
            let next = &sites[order[end]];
            let left = next.1 - 1 - RADIUS;
            if next.0 != first.0 || left > right { break; }
            if end - start >= MAX_GROUP_SITES || next.1 - 1 + RADIUS + 1 - (first.1 - 1 - RADIUS) > MAX_GROUP_SPAN { break; }
            right = right.max(next.1 - 1 + RADIUS + 1);
            end += 1;
        }
        groups.push((start, end, first.0.clone(), first.1 - 1 - RADIUS));
        start = end;
    }
    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads.max(1)).build()?;
    let results = pool.install(|| groups.par_iter().map(|&(lo, hi, ref chrom, left)| -> Result<Vec<(usize, FeatureRow)>> {
        let right = sites[order[hi - 1]].1 - 1 + RADIUS + 1;
        let mut bam = IndexedReader::from_path(bam_path)?; let fasta = rust_htslib::faidx::Reader::from_path(fasta_path)?;
        let ref_start = left.max(0); let reference = fasta.fetch_seq(chrom, ref_start as usize, (right - 1) as usize)?;
        let window = read_window(&mut bam, chrom, ref_start, right)?;
        let mut result = Vec::with_capacity(hi - lo);
        for &i in &order[lo..hi] { result.push((i, feature_vector(&sites[i], &reference, ref_start, &window, alu).map_err(|e| e.to_string()))); }
        Ok(result)
    }).collect::<Vec<_>>());
    let mut output = vec![Err("not computed".to_string()); sites.len()];
    for result in results { for (i, value) in result? { output[i] = value; } }
    Ok(output)
}

fn observe_block(bam_path: &Path, fasta_path: &Path, chrom: &str, start: i64, end: i64) -> Result<Vec<ObservationSite>> {
    let fasta = rust_htslib::faidx::Reader::from_path(fasta_path)?;
    let reference = fasta.fetch_seq(chrom, start as usize, (end - 1) as usize)?;
    let mut bam = IndexedReader::from_path(bam_path)?;
    let tid = bam.header().tid(chrom.as_bytes()).context("contig absent from BAM")?;
    bam.fetch((tid, start, end))?;
    let mut depth = vec![0u32; reference.len()]; let mut counts = vec![[0u32; 4]; reference.len()];
    for record in bam.records() {
        let record = record?;
        for (pos, observation) in read_observations(&record, start, end) {
            let local = usize::try_from(pos - start)?; let Some(base) = base_index(observation.base) else { continue; }; if base_index(reference[local]).is_none() { continue; }
            depth[local] += 1; counts[local][base] += 1;
        }
    }
    let mut sites = Vec::new();
    for (local, &ref_base) in reference.iter().enumerate() {
        let (alt, alt_i) = match ref_base.to_ascii_uppercase() { b'A' => (b'G', 2), b'T' => (b'C', 1), _ => continue };
        if counts[local][alt_i] == 0 { continue; }
        sites.push(ObservationSite { chrom: chrom.to_owned(), pos: start + local as i64 + 1, ref_base, counts: counts[local], depth: depth[local] });
        let _ = alt;
    }
    Ok(sites)
}

pub fn discover_sites(bam_path: &Path, fasta_path: &Path, out: &Path, threads: usize) -> Result<usize> {
    let fai = std::fs::read_to_string(format!("{}.fai", fasta_path.display()))?;
    let lengths: HashMap<String, i64> = fai.lines().filter_map(|line| { let mut f = line.split('\t'); Some((f.next()?.to_owned(), f.next()?.parse().ok()?)) }).collect();
    let reader = IndexedReader::from_path(bam_path)?; let mut contigs = Vec::new();
    for tid in 0..reader.header().target_count() { let chrom = String::from_utf8(reader.header().tid2name(tid).to_vec())?; if let Some(&len) = lengths.get(&chrom) { contigs.push((chrom, len)); } }
    let mut blocks = Vec::new(); for (chrom, len) in contigs { let mut start = 0; while start < len { blocks.push((chrom.clone(), start, (start + 8_000_000).min(len))); start += 8_000_000; } }
    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads.max(1)).build()?;
    let mut all = pool.install(|| blocks.par_iter().map(|(chrom, start, end)| observe_block(bam_path, fasta_path, chrom, *start, *end)).collect::<Result<Vec<_>>>())?.into_iter().flatten().collect::<Vec<_>>();
    all.sort_by(|a,b| a.chrom.cmp(&b.chrom).then(a.pos.cmp(&b.pos)));
    let mut text = String::new(); for site in &all { let alt = if site.ref_base == b'A' { b'G' } else { b'C' }; let freq = site.counts[if alt == b'G' { 2 } else { 1 }] as f64 / site.depth.max(1) as f64; use std::fmt::Write as _; writeln!(text, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{freq:.8}\t", site.chrom, site.pos-1, site.pos, site.ref_base as char, site.depth, site.counts[0], site.counts[1], site.counts[2], site.counts[3])?; }
    std::fs::write(out, text)?; Ok(all.len())
}
