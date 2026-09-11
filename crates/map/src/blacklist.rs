//! Retrocopy / multicopy blacklist generation (Kleinman & Majewski 2012,
//! RNA 18:87 — the intronless-retrocopy screen; multicopy-exon windows per
//! the production design).
//!
//! Simulated 100 bp queries are scanned against the packed index in
//! gap-free, no-splice mode with at most `MAX_MM` mismatches. Seeding uses
//! every 15-mer of the query (the index's own k): with at most 2 mismatches
//! the three clean segments between them include one of >= 33 bp, so a clean
//! 15-mer survives in every accepted hit and recall is exact for the
//! mm <= 2 contract. High-frequency 15-mers are capped (repeat loci are the
//! point of the blacklist, but a single 15-mer may not contribute more than
//! `KMER_POS_CAP` candidates).

use std::collections::{BTreeMap, HashSet};
use std::path::Path;

use rayon::prelude::*;

use crate::fasta::Base;
use crate::gtf::extract_attr;
use crate::index::{unpack_pos, Index};
use crate::seed::Strand;

/// Simulated read length.
pub const READ_LEN: usize = 100;
/// Maximum mismatches for a hit (gap-free).
pub const MAX_MM: u32 = 2;
/// Per-15-mer position cap.
const KMER_POS_CAP: u32 = 2048;
/// Per-query candidate cap (diagonal space).
const CAND_CAP: usize = 8192;
/// Retrocopy interval padding around a verified ectopic span.
pub const RETRO_PAD: u32 = 200;
/// Multicopy window padding.
pub const MULTI_PAD: u32 = 10;
/// Multicopy: a second locus within this strict-score band (two mismatch
/// equivalents) of the best marks the window.
pub const MULTI_BAND: i32 = 12;

/// One gap-free placement of a query.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Hit {
    /// Contig index.
    pub contig: u32,
    /// 0-based start of the aligned span.
    pub pos: u32,
    /// Strand the query aligned in.
    pub strand: Strand,
    /// Mismatch count over the full span.
    pub mm: u32,
}

fn encode_kmer(seq: &[u8], i: usize, k: usize) -> Option<(u64, Strand)> {
    // 2-bit canonical encoding of seq[i..i+k], matching the index's kmer
    // space (canonical = min(forward, revcomp); strand = whether forward won).
    let mut fwd: u64 = 0;
    let mut rev: u64 = 0;
    for j in 0..k {
        let c = u64::from(Base::from_ascii(seq[i + j]).code()?);
        fwd = (fwd << 2) | c;
        rev = (rev >> 2) | ((3 - c) << (2 * (k - 1)));
    }
    Some(if fwd <= rev {
        (fwd, Strand::Plus)
    } else {
        (rev, Strand::Minus)
    })
}

fn hamming(query: &[u8], refw: &[u8]) -> Option<u32> {
    let mut mm = 0u32;
    for (&q, &r) in query.iter().zip(refw) {
        if !q.eq_ignore_ascii_case(&r) {
            mm += 1;
            if mm > MAX_MM {
                return None;
            }
        }
    }
    Some(mm)
}

/// All gap-free placements of `query` with at most `MAX_MM` mismatches.
pub fn scan_query(index: &Index, query: &[u8]) -> Vec<Hit> {
    let k = index.params.k as usize;
    let l = query.len();
    if l < k {
        return Vec::new();
    }
    let mut cands: HashSet<(u32, u32, bool)> = HashSet::new();
    'kmers: for i in 0..=(l - k) {
        let Some((kmer, qstrand)) = encode_kmer(query, i, k) else {
            continue;
        };
        let hit = index.query(kmer);
        if hit.count > KMER_POS_CAP {
            continue;
        }
        for &p in hit.positions {
            let (contig, rpos, rstrand) = unpack_pos(p);
            let strand = Strand::xor(qstrand, rstrand);
            // Reference start of the full query span, from the diagonal.
            let start = match strand {
                Strand::Plus => rpos as i64 - i as i64,
                Strand::Minus => rpos as i64 - (l - i - k) as i64,
            };
            if start < 0 {
                continue;
            }
            let contig_len = index.reference.contigs[contig as usize].len as i64;
            if start + l as i64 > contig_len {
                continue;
            }
            cands.insert((contig, start as u32, strand == Strand::Minus));
            if cands.len() > CAND_CAP {
                break 'kmers;
            }
        }
    }
    let mut out = Vec::new();
    for (contig, pos, minus) in cands {
        let strand = if minus { Strand::Minus } else { Strand::Plus };
        let window = index.reference.contigs[contig as usize].slice_ascii(pos, pos + l as u32);
        let ok = match strand {
            Strand::Plus => hamming(query, &window),
            Strand::Minus => hamming(&revcomp_of(query), &window),
        };
        if let Some(mm) = ok {
            out.push(Hit {
                contig,
                pos,
                strand,
                mm,
            });
        }
    }
    out.sort_by_key(|h| h.mm);
    out
}

fn revcomp_of(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

/// A transcript's exon list (0-based half-open, ascending), with strand.
struct TxExons {
    contig: String,
    minus: bool,
    exons: Vec<(u32, u32)>,
}

fn read_gtf_transcripts(gtf: &Path) -> Vec<TxExons> {
    let text = std::fs::read_to_string(gtf).expect("read gtf");
    let mut map: BTreeMap<String, TxExons> = BTreeMap::new();
    for line in text.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 9 || c[2] != "exon" {
            continue;
        }
        let (Ok(s), Ok(e)) = (c[3].parse::<u32>(), c[4].parse::<u32>()) else {
            continue;
        };
        let Some(tid) = extract_attr(c[8], "transcript_id") else {
            continue;
        };
        let tx = map.entry(tid).or_insert_with(|| TxExons {
            contig: c[0].to_string(),
            minus: c[6] == "-",
            exons: Vec::new(),
        });
        tx.exons.push((s - 1, e));
    }
    for tx in map.values_mut() {
        tx.exons.sort_unstable();
        tx.exons.dedup();
    }
    map.into_values().collect()
}

fn merge_write(path: &Path, mut ivs: Vec<(String, u32, u32)>) -> (usize, u64) {
    ivs.sort();
    let mut merged: Vec<(String, u32, u32)> = Vec::new();
    for (c, s, e) in ivs {
        if let Some(last) = merged.last_mut() {
            if last.0 == c && s <= last.2 {
                last.2 = last.2.max(e);
                continue;
            }
        }
        merged.push((c, s, e));
    }
    let mut bp = 0u64;
    let mut text = String::new();
    for (c, s, e) in &merged {
        bp += u64::from(e - s);
        text.push_str(&format!("{c}\t{s}\t{e}\n"));
    }
    std::fs::write(path, text).expect("write blacklist bed");
    (merged.len(), bp)
}


/// Repeat family at a locus (RepeatMasker bed column 5), for the multicopy
/// same-family exemption: two AluY copies do not make a multicopy window.
struct FamilyMap {
    by_contig: BTreeMap<String, Vec<(u32, u32, String)>>,
}

impl FamilyMap {
    fn load(path: &Path) -> FamilyMap {
        let text = std::fs::read_to_string(path).expect("read rmsk bed");
        let mut by_contig: BTreeMap<String, Vec<(u32, u32, String)>> = BTreeMap::new();
        for line in text.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 5 {
                continue;
            }
            let (Ok(s0), Ok(e0)) = (f[1].parse::<u32>(), f[2].parse::<u32>()) else {
                continue;
            };
            by_contig
                .entry(f[0].to_string())
                .or_default()
                .push((s0, e0, f[4].to_string()));
        }
        for v in by_contig.values_mut() {
            v.sort();
        }
        FamilyMap { by_contig }
    }

    fn family_at(&self, contig: &str, pos: u32) -> Option<&str> {
        let v = self.by_contig.get(contig)?;
        // intervals are start-sorted; walk back from the last start <= pos
        let i = v.partition_point(|iv| iv.0 <= pos);
        (0..i).rev().take(4).find_map(|j| {
            let iv = &v[j];
            (iv.0 <= pos && pos < iv.1).then_some(iv.2.as_str())
        })
    }
}

/// Build the retrocopy blacklist: junction-crossing simulated reads whose
/// unspliced ectopic placement is unique-best and misses the source exons.
pub fn build_retrocopy(index: &Index, gtf: &Path, out: &Path) -> (usize, u64) {
    let txs = read_gtf_transcripts(gtf);
    let contig_of: BTreeMap<String, u32> = index
        .reference
        .contigs
        .iter()
        .enumerate()
        .map(|(i, c)| (c.name.clone(), i as u32))
        .collect();
    let intervals: Vec<(String, u32, u32)> = txs
        .par_iter()
        .flat_map_iter(|tx| {
            let mut out = Vec::new();
            let Some(&ci) = contig_of.get(&tx.contig) else {
                return out;
            };
            let contig = &index.reference.contigs[ci as usize];
            for w in tx.exons.windows(2) {
                let (e1, e2) = (w[0], w[1]);
                // Junction: e1 end -> e2 start. Simulated reads crossing it.
                for offset in (15..=85usize).step_by(5) {
                    // The read covers [junction - offset, junction - offset +
                    // 100) in transcript space: the tail of e1 + head of e2.
                    let take1 = offset as u32;
                    let take2 = READ_LEN as u32 - take1;
                    if e1.1 < take1 || e2.0 + take2 > contig.len {
                        continue;
                    }
                    let mut parts = contig.slice_ascii(e1.1 - take1, e1.1);
                    parts.extend(contig.slice_ascii(e2.0, e2.0 + take2));
                    if parts.len() != READ_LEN {
                        continue;
                    }
                    let query = if tx.minus { revcomp_of(&parts) } else { parts };
                    let hits = scan_query(index, &query);
                    if hits.is_empty() {
                        continue;
                    }
                    let best_mm = hits[0].mm;
                    let n_best = hits.iter().filter(|h| h.mm == best_mm).count();
                    if n_best != 1 {
                        continue;
                    }
                    let h = hits[0];
                    // The ectopic copy must not overlap the source exons.
                    let (hs, he) = (h.pos, h.pos + READ_LEN as u32);
                    let spans_e1 = h.contig == ci && hs < e1.1 && he > e1.0;
                    let spans_e2 = h.contig == ci && hs < e2.1 && he > e2.0;
                    if spans_e1 || spans_e2 {
                        continue;
                    }
                    out.push((
                        index.reference.contigs[h.contig as usize].name.clone(),
                        h.pos.saturating_sub(RETRO_PAD),
                        (h.pos + READ_LEN as u32 + RETRO_PAD)
                            .min(index.reference.contigs[h.contig as usize].len),
                    ));
                }
            }
            out
        })
        .collect();
    merge_write(out, intervals)
}

/// Build the multicopy blacklist: exonic 100 bp windows (step 50) with a
/// second near-equal gap-free locus. A second locus inside the SAME repeat
/// family (RepeatMasker family column, e.g. AluY<->AluY) does not count —
/// that is one repeat family's ordinary self-similarity, not a multicopy
/// exon. `rmsk` may be None (no family information: every distinct locus
/// counts).
pub fn build_multicopy(index: &Index, gtf: &Path, rmsk: Option<&Path>, out: &Path) -> (usize, u64) {
    let txs = read_gtf_transcripts(gtf);
    let fams = rmsk.map(FamilyMap::load);
    let contig_of: BTreeMap<String, u32> = index
        .reference
        .contigs
        .iter()
        .enumerate()
        .map(|(i, c)| (c.name.clone(), i as u32))
        .collect();
    let intervals: Vec<(String, u32, u32)> = txs
        .par_iter()
        .flat_map_iter(|tx| {
            let mut out = Vec::new();
            let Some(&ci) = contig_of.get(&tx.contig) else {
                return out;
            };
            let contig = &index.reference.contigs[ci as usize];
            let mut seen = HashSet::new();
            for &(es, ee) in &tx.exons {
                if (ee - es) < READ_LEN as u32 {
                    continue;
                }
                let mut ws = es;
                while ws + READ_LEN as u32 <= ee {
                    if seen.insert(ws) {
                        let query = contig.slice_ascii(ws, ws + READ_LEN as u32);
                        let query = if tx.minus { revcomp_of(&query) } else { query };
                        let hits = scan_query(index, &query);
                        if hits.len() >= 2 {
                            let best = hits[0];
                            let win_fam = fams
                                .as_ref()
                                .and_then(|f| f.family_at(contig.name.as_str(), ws));
                            let second_ok = hits[1..].iter().any(|h| {
                                let distinct =
                                    h.contig != best.contig || h.pos.abs_diff(best.pos) > 100;
                                // score band of 12 = two mismatch equivalents
                                let score_diff = 6 * (h.mm as i32 - best.mm as i32);
                                if !distinct || score_diff > MULTI_BAND {
                                    return false;
                                }
                                // Same-family copies are the repeat's own
                                // self-similarity, not a multicopy exon.
                                match (&fams, win_fam) {
                                    (Some(f), Some(wf)) => {
                                        f.family_at(
                                            index.reference.contigs[h.contig as usize]
                                                .name
                                                .as_str(),
                                            h.pos,
                                        ) != Some(wf)
                                    }
                                    _ => true,
                                }
                            });
                            if second_ok {
                                out.push((
                                    contig.name.clone(),
                                    ws.saturating_sub(MULTI_PAD),
                                    (ws + READ_LEN as u32 + MULTI_PAD).min(contig.len),
                                ));
                            }
                        }
                    }
                    ws += 50;
                }
            }
            out
        })
        .collect();
    merge_write(out, intervals)
}


#[cfg(test)]
mod family_tests {
    use super::*;
    use std::io::Write as _;

    #[test]
    fn family_lookup_point_containment() {
        let dir = std::env::temp_dir();
        let path = dir.join("esperanto_family_test.rmsk.bed");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "chr1\t100\t200\tSINE\tAluY").unwrap();
        writeln!(f, "chr1\t500\t600\tLINE\tL1HS").unwrap();
        drop(f);
        let fm = FamilyMap::load(&path);
        assert_eq!(fm.family_at("chr1", 150), Some("AluY"));
        assert_eq!(fm.family_at("chr1", 550), Some("L1HS"));
        assert_eq!(fm.family_at("chr1", 300), None);
        assert_eq!(fm.family_at("chr9", 150), None);
        std::fs::remove_file(&path).ok();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::parse_fasta_bytes;
    use crate::seed::SeedParams;

    /// xorshift64 pseudo-random DNA (no periodic k-mers).
    fn prand_seq(n: usize, mut x: u64) -> Vec<u8> {
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            out.push(b"ACGT"[(x % 4) as usize]);
        }
        out
    }

    /// chrA 4000bp random; chrB carries a 1-mismatch copy of chrA[500..600]
    /// at 1000..1100.
    fn toy_index() -> Index {
        let a = prand_seq(4000, 0x12345);
        let mut b = prand_seq(4000, 0x99999);
        let mut copy = a[500..600].to_vec();
        copy[50] = match copy[50] {
            b'A' => b'C',
            _ => b'A',
        };
        b[1000..1100].copy_from_slice(&copy);
        let fa = format!(
            ">chrA\n{}\n>chrB\n{}\n",
            String::from_utf8(a).unwrap(),
            String::from_utf8(b).unwrap()
        );
        let reference = parse_fasta_bytes(fa.as_bytes()).unwrap();
        Index::build(reference, SeedParams::rna_default())
    }

    #[test]
    fn scan_finds_both_copies_with_mm_counts() {
        let index = toy_index();
        let a = index.reference.contigs[0].slice_ascii(500, 600);
        let hits = scan_query(&index, &a);
        assert!(hits.len() >= 2, "hits: {hits:?}");
        assert_eq!(hits[0].mm, 0);
        assert_eq!(hits[0].contig, 0);
        assert_eq!(hits[0].pos, 500);
        assert!(hits.iter().any(|h| h.contig == 1 && h.pos == 1000 && h.mm == 1));
    }

    #[test]
    fn scan_rejects_beyond_two_mismatches() {
        let index = toy_index();
        let mut q = index.reference.contigs[0].slice_ascii(500, 600);
        // Avoid index 50, where the chrB copy carries its own edit.
        for i in [10, 30, 90] {
            q[i] = match q[i] {
                b'A' => b'C',
                _ => b'A',
            };
        }
        let hits = scan_query(&index, &q);
        assert!(hits.is_empty(), "3 mismatches must not scan: {hits:?}");
    }
}
