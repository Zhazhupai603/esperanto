//! # esperanto-verify
//!
//! Candidate-supporting read re-alignment: the BLAT re-alignment step of the
//! RNA-editing calling literature (Ramaswami et al., Nat Methods 9:490,
//! 2012; REDItools ships the same idea), performed by the pipeline's own
//! aligner instead of an external one.
//!
//! Every read that supported a candidate site is re-aligned against the
//! whole reference in exhaustive mode — every chained locus extended, no
//! decisive shortcuts — and its vote is kept only when the site-supporting
//! placement is *obviously best* (RASER's criterion, Ahn & Xiao,
//! Bioinformatics 31:866, 2015) and no competing placement explains the
//! very column that carried the edit as a plain reference match. A site
//! stays on the PASS track only when at least two verified votes remain.
//!
//! Invariants
//! ----------
//!  * Deterministic: reads are examined in (MAPQ desc, qname asc) order and
//!    every placement decision is a pure function of the sequence.
//!  * The input `candidates.bed` is preserved verbatim; the output appends
//!    two columns (`alt_verified`, `withdrawn`).

use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use rayon::prelude::*;
use rust_htslib::bam::{IndexedReader, Read as BamRead};

use esperanto_map::align::{AlignConfig, Aligner, ScoredPlacement};
use esperanto_map::extend::CigarOp;
use esperanto_map::index::Index;
use esperanto_map::index_io;
use esperanto_map::seed::Strand;

/// Errors raised by this crate.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Underlying I/O failure.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    /// BAM access failure.
    #[error("bam error: {0}")]
    Bam(String),
    /// Malformed candidates line.
    #[error("bad candidates line {line}: {reason}")]
    BadLine {
        /// 1-based line number.
        line: usize,
        /// What was wrong.
        reason: String,
    },
    /// Alignment index failure.
    #[error("index error: {0}")]
    Index(String),
}

/// Verify parameters (frozen defaults per the spec).
#[derive(Debug, Clone)]
pub struct VerifyParams {
    /// Coordinate-sorted, indexed BAM.
    pub bam: PathBuf,
    /// `scan`'s candidates.bed (10 columns).
    pub candidates: PathBuf,
    /// Alignment index (paidx).
    pub index: PathBuf,
    /// L1 engine bundle, when the map stage had one.
    pub l1_bundle: Option<PathBuf>,
    /// Output bed (input columns + alt_verified + withdrawn).
    pub out: PathBuf,
    /// Worker threads (0 = all).
    pub threads: usize,
    /// Maximum supporting reads examined per site (highest MAPQ first).
    pub max_reads: usize,
    /// "Obviously best" margin: best − second ≥ this (two mismatch
    /// equivalents) keeps a vote.
    pub ob_margin: i32,
    /// A competing placement within this score of the best withdraws the
    /// vote when it explains the edited column as a reference match (the
    /// BLAT rule). Kept equal to `ob_margin`: the BLAT check runs first, so
    /// the whole band up to the margin boundary is judged on the column.
    pub alt_margin: i32,
    /// Minimum verified votes for a site to stay on the PASS track.
    pub min_verified: u32,
    /// Minimum verified-vote fraction of depth.
    pub min_verified_frac: f64,
    /// Mismatch channel the candidates belong to (controls the verified
    /// alternate base).
    pub channel: Channel,
    /// Diagnostic probe sink; populated from `ESP_VERIFY_PROBE` by the runner.
    pub probe: Option<Arc<std::sync::Mutex<ProbeSink>>>,
}

impl Default for VerifyParams {
    fn default() -> Self {
        VerifyParams {
            bam: PathBuf::new(),
            candidates: PathBuf::new(),
            index: PathBuf::new(),
            l1_bundle: None,
            out: PathBuf::new(),
            threads: 0,
            max_reads: 50,
            ob_margin: 12,
            alt_margin: 12,
            min_verified: 2,
            min_verified_frac: 0.01,
            channel: Channel::Edit,
            probe: None,
        }
    }
}

/// Mismatch channel: which alternate base a candidate's votes must carry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Channel {
    /// Editing: A->G, T->C.
    #[default]
    Edit,
    /// Control 1: A->C, T->G.
    Ctrl1,
    /// Control 2: A->T, T->A.
    Ctrl2,
}

impl Channel {
    /// The verified alternate for a reference base (identity when the
    /// reference is not A/T, which then never verifies).
    pub fn alt_of(self, ref_base: u8) -> u8 {
        match (self, ref_base) {
            (Channel::Edit, b'A') => b'G',
            (Channel::Edit, b'T') => b'C',
            (Channel::Ctrl1, b'A') => b'C',
            (Channel::Ctrl1, b'T') => b'G',
            (Channel::Ctrl2, b'A') => b'T',
            (Channel::Ctrl2, b'T') => b'A',
            (_, other) => other,
        }
    }
}

/// Diagnostic probe sink (`ESP_VERIFY_PROBE`): sampled withdrawal detail
/// lines for the attribution reports.
/// lines for the attribution reports.
#[derive(Debug, Default)]
pub struct ProbeSink {
    /// Output path for the probe lines.
    pub path: PathBuf,
    /// Margin-rule samples (cap 500).
    pub margin_lines: Vec<String>,
    /// no_site samples (cap 100).
    pub nosite_lines: Vec<String>,
}

/// Run summary.
#[derive(Debug, Default, Clone)]
pub struct VerifyStats {
    /// Candidate sites examined.
    pub candidates: u64,
    /// Sites meeting the verified-vote floor.
    pub sites_kept: u64,
    /// Supporting reads re-aligned.
    pub reads_examined: u64,
    /// Votes withdrawn by the obviously-best / BLAT rules.
    pub votes_withdrawn: u64,
    /// Withdrawals where the best re-placement lost the site itself.
    pub wd_no_site: u64,
    /// Withdrawals on the obviously-best margin rule.
    pub wd_margin: u64,
    /// Withdrawals on the BLAT column-explanation rule.
    pub wd_blat: u64,
    /// Pair-anchored (`RE:Z:unmapped`) votes kept without solo re-alignment.
    pub pair_anchored_kept: u64,
    /// no_site subclass: exhaustive re-alignment returned no candidates.
    pub ns_empty: u64,
    /// no_site subclass: best re-placement sits on another locus.
    pub ns_other_locus: u64,
    /// no_site subclass: the site column falls inside a best-placement skip.
    pub ns_refskip: u64,
    /// no_site subclass: the site column falls inside a best-placement deletion.
    pub ns_del: u64,
    /// no_site subclass: the site column falls outside the best placement's span.
    pub ns_outside_span: u64,
    /// no_site subclass: best covers the column but not with the alt base.
    pub ns_base_mismatch: u64,
}

struct Site {
    chrom: String,
    pos0: u32,
    alt_base: u8,
    line: String,
}

/// Run the verify stage.
pub fn run(p: &VerifyParams) -> Result<VerifyStats, Error> {
    // The alignment index is process-lifetime (the map stage treats it the
    // same way): leak the handle so every worker shares one mapping.
    let index: &'static Index = Box::leak(Box::new(
        index_io::load(&p.index).map_err(|e| Error::Index(e.to_string()))?,
    ));
    let l1 = match &p.l1_bundle {
        Some(path) => Some(Arc::new(
            esperanto_engine::L1Index::open(path)
                .map_err(|e| Error::Index(format!("L1 bundle: {e}")))?,
        )),
        None => None,
    };
    run_with_index(p, index, l1)
}

/// Run the verify stage on an already-loaded index and L1 bundle (the flow
/// runtime hands the map stage's handles over instead of re-reading them).
pub fn run_with_index(
    p: &VerifyParams,
    index: &'static Index,
    l1: Option<Arc<esperanto_engine::L1Index>>,
) -> Result<VerifyStats, Error> {
    let p_owned;
    let p = match std::env::var_os("ESP_VERIFY_PROBE") {
        Some(path) => {
            let mut q = p.clone();
            q.probe = Some(Arc::new(std::sync::Mutex::new(ProbeSink {
                path: PathBuf::from(path),
                ..Default::default()
            })));
            p_owned = q;
            &p_owned
        }
        None => p,
    };
    let sites = parse_candidates(&p.candidates, p.channel)?;
    let bam_header = {
        let r = IndexedReader::from_path(&p.bam).map_err(|e| Error::Bam(e.to_string()))?;
        r.header().clone()
    };
    let tid_of = |chrom: &str| -> Option<u32> {
        (0..bam_header.target_count()).find(|&t| bam_header.tid2name(t) == chrom.as_bytes())
    };

    // Group sites by contig (scan emits per-contig blocks, so a plain
    // split keeps order); one BAM reader and one aligner per worker.
    let mut groups: Vec<(u32, Vec<usize>)> = Vec::new();
    for (i, s) in sites.iter().enumerate() {
        if let Some(tid) = tid_of(&s.chrom) {
            match groups.last_mut() {
                Some((t, v)) if *t == tid => v.push(i),
                _ => groups.push((tid, vec![i])),
            }
        }
    }

    let pool = match p.threads {
        0 => rayon::ThreadPoolBuilder::new().build(),
        n => rayon::ThreadPoolBuilder::new().num_threads(n).build(),
    }
    .map_err(|e| Error::Bam(format!("thread pool: {e}")))?;

    let results: Vec<(usize, u32, u32, [u32; 9], u32)> = pool.install(|| {
        groups
            .par_iter()
            .flat_map_iter(|(tid, idxs)| {
                let mut reader = IndexedReader::from_path(&p.bam)
                    .map_err(|e| Error::Bam(e.to_string()))
                    .ok()?;
                let mut aligner = Aligner::new(index, AlignConfig::rna_default());
                if let Some(l1) = l1.clone() {
                    aligner = aligner.with_l1(l1);
                }
                let mut out: Vec<(usize, u32, u32, [u32; 9], u32)> = Vec::with_capacity(idxs.len());
                for &i in idxs {
                    let s = &sites[i];
                    if reader
                        .fetch((*tid, i64::from(s.pos0), i64::from(s.pos0) + 1))
                        .is_err()
                    {
                        continue;
                    }
                    let (kept, withdrawn, reasons, pair) = verify_reads(s, &mut reader, &mut aligner, p);
                    out.push((i, kept, withdrawn, reasons, pair));
                }
                Some(out)
            })
            .flatten()
            .collect()
    });

    let mut kept_by_site = vec![(0u32, 0u32); sites.len()];
    let mut wd_reasons = [0u64; 9];
    let mut pair_total = 0u64;
    for (i, kept, withdrawn, reasons, pair) in results {
        kept_by_site[i] = (kept, withdrawn);
        pair_total += u64::from(pair);
        for r in 0..9 {
            wd_reasons[r] += u64::from(reasons[r]);
        }
    }

    let mut stats = VerifyStats {
        candidates: sites.len() as u64,
        wd_no_site: wd_reasons[0],
        wd_margin: wd_reasons[1],
        wd_blat: wd_reasons[2],
        pair_anchored_kept: pair_total,
        ns_empty: wd_reasons[3],
        ns_other_locus: wd_reasons[4],
        ns_refskip: wd_reasons[5],
        ns_del: wd_reasons[6],
        ns_outside_span: wd_reasons[7],
        ns_base_mismatch: wd_reasons[8],
        ..Default::default()
    };
    let mut out = std::io::BufWriter::new(std::fs::File::create(&p.out)?);
    for (s, (kept, withdrawn)) in sites.iter().zip(kept_by_site) {
        stats.reads_examined += (kept + withdrawn) as u64;
        stats.votes_withdrawn += withdrawn as u64;
        let depth: u32 = s
            .line
            .split('\t')
            .nth(4)
            .and_then(|d| d.parse().ok())
            .unwrap_or(0);
        let passes = f64::from(kept) / f64::from(depth.max(1)) >= p.min_verified_frac
            && kept >= p.min_verified;
        if passes {
            stats.sites_kept += 1;
        }
        writeln!(out, "{}\t{kept}\t{withdrawn}", s.line)
            .map_err(|e| Error::Io(std::io::Error::other(e.to_string())))?;
    }
    out.flush()?;
    if let Some(sink) = &p.probe {
        let g = sink
            .lock()
            .map_err(|e| Error::Io(std::io::Error::other(e.to_string())))?;
        let mut text = String::new();
        for l in g.margin_lines.iter().chain(g.nosite_lines.iter()) {
            text.push_str(l);
            text.push('\n');
        }
        std::fs::write(&g.path, text)?;
    }
    Ok(stats)
}

fn parse_candidates(path: &Path, channel: Channel) -> Result<Vec<Site>, Error> {
    let text = std::fs::read_to_string(path)?;
    let mut out = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let bad = |reason: &str| Error::BadLine {
            line: i + 1,
            reason: reason.to_string(),
        };
        if f.len() < 11 {
            return Err(bad("expected 11 columns (scan output)"));
        }
        let chrom = f[0].to_string();
        let pos0: u32 = f[1].parse().map_err(|_| bad("pos0 not a number"))?;
        let rb = f[3].as_bytes();
        if rb.is_empty() {
            return Err(bad("empty reference base"));
        }
        let ref_base = rb[0].to_ascii_uppercase();
        // The alternate is the channel-consistent base, not whichever
        // non-reference base happens to count highest. Anything else is
        // outside the channel's contract and never verifies.
        let alt_base = channel.alt_of(ref_base);
        out.push(Site {
            chrom,
            pos0,
            alt_base,
            line: line.to_string(),
        });
    }
    Ok(out)
}

/// Re-align the site's supporting reads and split kept / withdrawn votes.
fn verify_reads(
    s: &Site,
    reader: &mut IndexedReader,
    aligner: &mut Aligner<'static>,
    p: &VerifyParams,
) -> (u32, u32, [u32; 9], u32) {
    let mut supporting: Vec<(i32, String, Vec<u8>, u32)> = Vec::new();
    let mut kept_pair = 0u32;
    for r in reader.records() {
        let Ok(rec) = r else { continue };
        if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() {
            continue;
        }
        if rec.mapq() < esperanto_map::mapq::PILEUP_GATE {
            continue;
        }
        let Some(qoff) = query_offset_at(&rec, s.pos0) else {
            continue;
        };
        let seq_ref_frame = rec.seq().as_bytes();
        if qoff >= seq_ref_frame.len() {
            continue;
        }
        if seq_ref_frame[qoff].to_ascii_uppercase() != s.alt_base {
            continue;
        }
        if rec.qual().get(qoff).map_or(0, |q| *q) < 30 {
            continue;
        }
        // A mate-rescued placement already passed pair arbitration at map
        // time (mate confidence, in-window re-alignment, same strict scoring
        // domain, mate-consistent competitors only); its MAPQ is capped at
        // the pileup gate. Re-judging it solo uses weaker information to
        // overturn a stronger verdict, so the vote stands as kept.
        if matches!(
            rec.aux(b"RE"),
            Ok(rust_htslib::bam::record::Aux::String("unmapped"))
        ) {
            kept_pair += 1;
            continue;
        }
        // SEQ is reference-forward; recover the read as sequenced
        // (reverse-complementing is an involution).
        let original = if rec.is_reverse() {
            esperanto_map::align::revcomp(&seq_ref_frame)
        } else {
            seq_ref_frame
        };
        supporting.push((
            rec.mapq() as i32,
            String::from_utf8_lossy(rec.qname()).into_owned(),
            original,
            rec.pos() as u32,
        ));
    }
    supporting.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    supporting.truncate(p.max_reads);

    let mut kept = 0u32;
    let mut withdrawn = 0u32;
    let mut reasons = [0u32; 9];
    for (_, _, seq, orig_pos) in supporting {
        let cands = aligner.align_read_exhaustive(&seq);
        let index: &'static Index = aligner.index;
        let chrom = s.chrom.clone();
        let is_site_contig = move |c: u32| -> bool {
            index.reference.contigs[c as usize].name == chrom
        };
        let ref_base_of =
            |c: u32, pos: u32| -> u8 { ref_base(index, c, pos) };
        let outcome = vote_kept(
            &cands,
            s,
            orig_pos,
            &seq,
            p.ob_margin,
            p.alt_margin,
            &is_site_contig,
            &ref_base_of,
        );
        match outcome {
            VoteOutcome::Kept => kept += 1,
            VoteOutcome::NoBestAtSite => {
                withdrawn += 1;
                reasons[0] += 1;
                let sub = classify_no_site(&cands, s, orig_pos, &is_site_contig);
                reasons[3 + sub] += 1;
                if let Some(sink) = &p.probe {
                    let mut g = sink.lock().unwrap();
                    if g.nosite_lines.len() < 100 {
                        g.nosite_lines.push(no_site_probe_line(sub, &cands, s, index));
                    }
                }
            }
            VoteOutcome::Margin => {
                withdrawn += 1;
                reasons[1] += 1;
                if let Some(sink) = &p.probe {
                    let mut g = sink.lock().unwrap();
                    if g.margin_lines.len() < 500 {
                        g.margin_lines.push(margin_probe_line(&cands, s, index));
                    }
                }
            }
            VoteOutcome::Blat => {
                withdrawn += 1;
                reasons[2] += 1;
            }
        }
    }
    (kept + kept_pair, withdrawn, reasons, kept_pair)
}

/// no_site withdrawal subclasses, in counter order.
const NS_LABELS: [&str; 6] = [
    "empty",
    "other_locus",
    "refskip",
    "del",
    "outside_span",
    "base_mismatch",
];

/// Which way the best re-placement lost the site (counter index into
/// the tail of the reasons array).
fn classify_no_site(
    cands: &[ScoredPlacement],
    s: &Site,
    orig_pos: u32,
    is_site_contig: &dyn Fn(u32) -> bool,
) -> usize {
    let Some(best) = cands.first() else {
        return 0;
    };
    let same_locus = is_site_contig(best.contig) && best.pos.abs_diff(orig_pos) <= 5;
    if !same_locus {
        return 1;
    }
    if placement_offset_at(best, s.pos0).is_some() {
        // Same locus and the column is covered: the vote can only have
        // failed on the alt base itself.
        return 5;
    }
    let mut r = u64::from(best.pos);
    for op in &best.cigar {
        let n = u64::from(op_len(op));
        match op {
            CigarOp::Match(_) => r += n,
            CigarOp::RefSkip(_) => {
                if r <= u64::from(s.pos0) && u64::from(s.pos0) < r + n {
                    return 2;
                }
                r += n;
            }
            CigarOp::Del(_) => {
                if r <= u64::from(s.pos0) && u64::from(s.pos0) < r + n {
                    return 3;
                }
                r += n;
            }
            _ => {}
        }
    }
    4
}

fn span_of(p: &ScoredPlacement) -> (u32, u32) {
    let lo = p.blocks.iter().map(|b| b.0).min().unwrap_or(p.pos);
    let hi = p.blocks.iter().map(|b| b.1).max().unwrap_or(p.pos);
    (lo, hi)
}

/// Reference distance between two spans (0 when they overlap).
fn span_gap(a: (u32, u32), b: (u32, u32)) -> u32 {
    a.0.max(b.0).saturating_sub(a.1.min(b.1))
}

fn has_refskip(p: &ScoredPlacement) -> bool {
    p.cigar.iter().any(|op| matches!(op, CigarOp::RefSkip(_)))
}

fn cigar_string(cigar: &[CigarOp]) -> String {
    let mut s = String::new();
    for op in cigar {
        let c = match op {
            CigarOp::Match(_) => 'M',
            CigarOp::Ins(_) => 'I',
            CigarOp::Del(_) => 'D',
            CigarOp::RefSkip(_) => 'N',
            CigarOp::SoftClip(_) => 'S',
        };
        s.push_str(&format!("{}{}", op_len(op), c));
    }
    s
}

/// Query-side intervals covered by M columns.
fn query_m_intervals(cigar: &[CigarOp]) -> Vec<(u32, u32)> {
    let mut q = 0u32;
    let mut v = Vec::new();
    for op in cigar {
        match op {
            CigarOp::Match(n) => {
                v.push((q, q + *n));
                q += *n;
            }
            CigarOp::Ins(n) | CigarOp::SoftClip(n) => q += *n,
            _ => {}
        }
    }
    v
}

/// Overlap of two query-side M-cover sets, normalized by the smaller.
fn q_overlap_frac(a: &[(u32, u32)], b: &[(u32, u32)]) -> f64 {
    let mut ov = 0u64;
    for &(s1, e1) in a {
        for &(s2, e2) in b {
            ov += u64::from(e1.min(e2).saturating_sub(s1.max(s2)));
        }
    }
    let la: u64 = a.iter().map(|&(s, e)| u64::from(e - s)).sum();
    let lb: u64 = b.iter().map(|&(s, e)| u64::from(e - s)).sum();
    let d = la.min(lb);
    if d == 0 {
        0.0
    } else {
        ov as f64 / d as f64
    }
}

fn describe(p: &ScoredPlacement, index: &Index) -> String {
    let sp = span_of(p);
    format!(
        "{}:{}-{} score={} cigar={} refskip={}",
        index.reference.contigs[p.contig as usize].name,
        sp.0,
        sp.1,
        p.score,
        cigar_string(&p.cigar),
        has_refskip(p) as u8
    )
}

/// One sampled margin withdrawal: best and the margin-breaking runner-up.
fn margin_probe_line(cands: &[ScoredPlacement], s: &Site, index: &Index) -> String {
    let best = &cands[0];
    let second = cands
        .iter()
        .skip(1)
        .filter(|c| !same_locus_as(best, c))
        .max_by_key(|c| c.score);
    match second {
        None => format!(
            "M\t{}:{}\tbest={}\tsecond=none",
            s.chrom,
            s.pos0,
            describe(best, index)
        ),
        Some(sec) => {
            let same_c = best.contig == sec.contig;
            let gap = if same_c {
                span_gap(span_of(best), span_of(sec)).to_string()
            } else {
                "NA".to_string()
            };
            let qov = q_overlap_frac(
                &query_m_intervals(&best.cigar),
                &query_m_intervals(&sec.cigar),
            );
            format!(
                "M\t{}:{}\tbest={}\tsecond={}\tsame_contig={} gap={} q_overlap={:.3}",
                s.chrom,
                s.pos0,
                describe(best, index),
                describe(sec, index),
                same_c as u8,
                gap,
                qov
            )
        }
    }
}

/// One sampled no_site withdrawal: the subclass and the best re-placement.
fn no_site_probe_line(sub: usize, cands: &[ScoredPlacement], s: &Site, index: &Index) -> String {
    let best = cands
        .first()
        .map(|p| describe(p, index))
        .unwrap_or_else(|| "none".to_string());
    format!(
        "NS\t{}:{}\tsub={}\tbest={}",
        s.chrom,
        s.pos0,
        NS_LABELS[sub],
        best
    )
}

/// Query offset of reference position `pos0` inside a BAM record (pos0 must
/// fall in an M/=/X column).
fn query_offset_at(rec: &rust_htslib::bam::Record, pos0: u32) -> Option<usize> {
    use rust_htslib::bam::record::Cigar::*;
    let mut qpos = 0usize;
    let mut rpos = rec.pos() as u64;
    for op in rec.cigar().iter() {
        match op {
            Match(n) | Equal(n) | Diff(n) => {
                if rpos <= u64::from(pos0) && u64::from(pos0) < rpos + u64::from(*n) {
                    return Some(qpos + (u64::from(pos0) - rpos) as usize);
                }
                qpos += *n as usize;
                rpos += u64::from(*n);
            }
            Ins(n) | SoftClip(n) => qpos += *n as usize,
            Del(n) | RefSkip(n) => rpos += u64::from(*n),
            _ => {}
        }
    }
    None
}

/// Why a single vote was kept or withdrawn (per-rule counters feed the
/// run summary's attribution line).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VoteOutcome {
    /// Vote stands.
    Kept,
    /// The best re-placement does not cover the site with the alt base at
    /// the original locus.
    NoBestAtSite,
    /// A non-overlapping runner-up is within `ob_margin` of the best.
    Margin,
    /// A competitor within `alt_margin` explains the edited column as a
    /// reference match.
    Blat,
}

/// The withdrawal rules over the exhaustive candidate list:
/// 1. the best placement must still support the site (cover pos0 with the
///    alt base at the original locus ±5bp);
/// 2. the best must beat the best non-overlapping runner-up by `ob_margin`;
/// 3. no competitor within `alt_margin` may explain the edited column as a
///    reference match (the BLAT rule).
///
/// Frame convention: a placement's CIGAR consumes the read in the
/// placement's own orientation — the original read for plus placements,
/// its reverse complement for minus ones. Column offsets therefore
/// translate between placements through the original-read index:
/// plus qoff = i, minus qoff = len − 1 − i.
#[allow(clippy::too_many_arguments)]
fn vote_kept(
    cands: &[ScoredPlacement],
    s: &Site,
    orig_pos: u32,
    seq: &[u8],
    ob_margin: i32,
    alt_margin: i32,
    is_site_contig: &dyn Fn(u32) -> bool,
    ref_base_of: &dyn Fn(u32, u32) -> u8,
) -> VoteOutcome {
    let Some(best) = cands.first() else {
        return VoteOutcome::NoBestAtSite;
    };
    let Some(best_qoff) = placement_offset_at(best, s.pos0) else {
        return VoteOutcome::NoBestAtSite;
    };
    let same_locus = is_site_contig(best.contig) && best.pos.abs_diff(orig_pos) <= 5;
    let alt_supported = base_at(seq, best.strand, best_qoff) == s.alt_base;
    if !same_locus || !alt_supported {
        return VoteOutcome::NoBestAtSite;
    }
    // Original-read index of the site column under the winning placement.
    let i_orig = match best.strand {
        Strand::Plus => best_qoff,
        Strand::Minus => seq.len().saturating_sub(best_qoff + 1),
    };
    // BLAT first: a competitor within alt_margin that explains the edited
    // column as a reference match withdraws the vote. alt_margin == ob_margin,
    // so the whole band up to the margin boundary is judged on the column
    // before the margin rule can swallow it.
    for c in cands.iter().skip(1) {
        if same_locus_as(best, c) || c.score < best.score - alt_margin {
            continue;
        }
        // The same original-read base, re-expressed in the competitor's frame.
        let qoff_c = match c.strand {
            Strand::Plus => i_orig,
            Strand::Minus => seq.len().saturating_sub(i_orig + 1),
        };
        if let Some(rpos) = placement_ref_at_offset(c, qoff_c) {
            let rb = ref_base_of(c.contig, rpos);
            if rb != 0 && rb == base_at(seq, c.strand, qoff_c) {
                return VoteOutcome::Blat;
            }
        }
    }
    let second = cands
        .iter()
        .skip(1)
        .filter(|c| !same_locus_as(best, c))
        .map(|c| c.score)
        .max();
    // No runner-up at all: the placement is obviously best by default.
    // (Sentinel arithmetic would overflow: score - i32::MIN wraps negative
    // and withdraws every single-candidate vote.)
    if let Some(second) = second {
        if best.score - second < ob_margin {
            return VoteOutcome::Margin;
        }
    }
    VoteOutcome::Kept
}

/// Two placements are the same locus when their M blocks overlap at all.
fn same_locus_as(a: &ScoredPlacement, b: &ScoredPlacement) -> bool {
    a.contig == b.contig && blocks_overlap(&a.blocks, &b.blocks)
}

fn blocks_overlap(a: &[(u32, u32)], b: &[(u32, u32)]) -> bool {
    a.iter()
        .any(|(s1, e1)| b.iter().any(|(s2, e2)| s1.max(s2) < e1.min(e2)))
}

/// Query offset of reference position `pos0` in a placement (M columns
/// only); None when pos0 is outside its aligned blocks.
fn placement_offset_at(p: &ScoredPlacement, pos0: u32) -> Option<usize> {
    let mut q = 0usize;
    let mut r = u64::from(p.pos);
    for op in &p.cigar {
        let n = u64::from(op_len(op));
        match op {
            CigarOp::Match(_) => {
                if r <= u64::from(pos0) && u64::from(pos0) < r + n {
                    return Some(q + (u64::from(pos0) - r) as usize);
                }
                q += n as usize;
                r += n;
            }
            CigarOp::Ins(_) | CigarOp::SoftClip(_) => q += n as usize,
            CigarOp::Del(_) | CigarOp::RefSkip(_) => r += n,
        }
    }
    None
}

/// Reference position a query offset maps to inside a placement (M columns
/// only); None when the offset is not in an aligned block.
fn placement_ref_at_offset(p: &ScoredPlacement, qoff: usize) -> Option<u32> {
    let mut q = 0usize;
    let mut r = u64::from(p.pos);
    for op in &p.cigar {
        let n = u64::from(op_len(op));
        match op {
            CigarOp::Match(_) => {
                if q <= qoff && qoff < q + n as usize {
                    return Some((r + (qoff - q) as u64) as u32);
                }
                q += n as usize;
                r += n;
            }
            CigarOp::Ins(_) | CigarOp::SoftClip(_) => q += n as usize,
            CigarOp::Del(_) | CigarOp::RefSkip(_) => r += n,
        }
    }
    None
}

fn op_len(op: &CigarOp) -> u32 {
    match op {
        CigarOp::Match(n)
        | CigarOp::Ins(n)
        | CigarOp::Del(n)
        | CigarOp::RefSkip(n)
        | CigarOp::SoftClip(n) => *n,
    }
}

/// The read base at a placement-frame query offset, read in the reference
/// frame. Plus placements consume the original read directly; minus
/// placements consume its reverse complement, so offset `qoff` addresses
/// original index `len − 1 − qoff` and the base complements.
fn base_at(seq: &[u8], strand: Strand, qoff: usize) -> u8 {
    match strand {
        Strand::Plus => seq.get(qoff).copied().unwrap_or(b'N').to_ascii_uppercase(),
        Strand::Minus => seq
            .len()
            .checked_sub(qoff + 1)
            .and_then(|i| seq.get(i).copied())
            .map(|b| complement(b.to_ascii_uppercase()))
            .unwrap_or(b'N'),
    }
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

fn ref_base(index: &'static Index, contig: u32, pos: u32) -> u8 {
    let c = &index.reference.contigs[contig as usize];
    if pos < c.len {
        c.base(pos).to_ascii()
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn placement(pos: u32, cigar: &[CigarOp]) -> ScoredPlacement {
        ScoredPlacement {
            contig: 0,
            pos,
            strand: Strand::Plus,
            score: 100,
            blocks: esperanto_map::align::ref_m_blocks_public(pos, cigar),
            cigar: cigar.to_vec(),
        }
    }

    #[test]
    fn column_mapping_round_trip() {
        let p = placement(
            100,
            &[CigarOp::Match(10), CigarOp::RefSkip(50), CigarOp::Match(10)],
        );
        assert_eq!(placement_offset_at(&p, 105), Some(5));
        assert_eq!(placement_offset_at(&p, 75), None); // before the placement
        assert_eq!(placement_offset_at(&p, 115), None); // inside the skip
        assert_eq!(placement_offset_at(&p, 165), Some(15)); // second block
        assert_eq!(placement_ref_at_offset(&p, 15), Some(165));
        assert_eq!(placement_ref_at_offset(&p, 5), Some(105));
    }

    #[test]
    fn overlap_detection() {
        let a = placement(0, &[CigarOp::Match(100)]);
        let b = placement(150, &[CigarOp::Match(100)]);
        assert!(same_locus_as(&a, &placement(50, &[CigarOp::Match(100)])));
        assert!(!same_locus_as(&a, &b));
    }

    fn site(pos0: u32, alt: u8) -> Site {
        Site {
            chrom: "chrT".to_string(),
            pos0,
            alt_base: alt,
            line: String::new(),
        }
    }

    /// Minus-strand read: base_at must address the original read mirrored.
    /// Original read of length 20, site column at reference pos 100 under a
    /// minus 20M placement starting at 90 → placement-frame offset 10 →
    /// original index 9. A 'G' at original index 9 reads 'C' on the
    /// reference frame (a T>C edit on a minus transcript).
    #[test]
    fn minus_strand_vote_kept() {
        let mut seq = vec![b'A'; 20];
        seq[9] = b'G'; // original index 9: complements to C at the site column
        seq[10] = b'T'; // the buggy index would land here (T -> A != C)
        let mut best = placement(90, &[CigarOp::Match(20)]);
        best.strand = Strand::Minus;
        best.score = 200;
        let mut rival = placement(5000, &[CigarOp::Match(20)]);
        rival.score = 180; // outside the alt margin, inside nothing
        let cands = vec![best, rival];
        let s = site(100, b'C');
        assert_eq!(
            vote_kept(
                &cands,
                &s,
                90,
                &seq,
                12,
                6,
                &|c| c == 0,
                &|_, _| 0,
            ),
            VoteOutcome::Kept
        );
    }

    /// BLAT withdrawal across frames: a plus-strand competitor within the
    /// margin maps the same original base to a matching reference base —
    /// the vote is withdrawn even though the competitor sits on the other
    /// strand.
    #[test]
    fn blat_withdrawal_across_strands() {
        let mut seq = vec![b'A'; 20];
        seq[9] = b'G';
        let mut best = placement(90, &[CigarOp::Match(20)]);
        best.strand = Strand::Minus;
        best.score = 200;
        let mut rival = placement(5000, &[CigarOp::Match(20)]);
        rival.strand = Strand::Plus;
        rival.contig = 1; // another contig, so the BLAT closure's G lands
        rival.score = 197; // within alt_margin (200 - 6 <= 197)
        // best is minus: i_orig = 20 - 1 - 10 = 9; rival is plus: qoff_c = 9.
        // Reference at the rival's mapped position returns 'G' == seq[9]
        // read plus-frame → the column matches elsewhere → withdraw.
        let cands = vec![best, rival];
        let s = site(100, b'C');
        // BLAT runs before the margin rule, so production margins (12, 6)
        // reach this branch: the rival is inside alt_margin (200 - 6 <= 197).
        assert_eq!(
            vote_kept(
                &cands,
                &s,
                90,
                &seq,
                12,
                6,
                &|c| c == 0,
                &|c, _| if c == 1 { b'G' } else { 0 },
            ),
            VoteOutcome::Blat
        );
    }
}
