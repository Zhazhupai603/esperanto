//! Track 2: junction-kmer direct relocation index (independent system).
//!
//! Reads that Track 1 leaves unmapped or heavily soft-clipped are relocated
//! through annotated junctions via exact k-mer hits: each junction contributes
//! the crossing 15-mers of its donor-tail + acceptor-head pseudo-sequence
//! (transcript orientation; minus-strand junctions reverse-complemented),
//! stored with the split offset of the breakpoint inside the k-mer and an
//! A-position mask whose subsets enumerate A→G editing variants. Queries
//! aggregate hits per junction (≥2 required), infer the read breakpoint as a
//! ±1bp mode, and confirm locally with an affine Smith-Waterman variant
//! before assembling a Track-2 record. Nothing here reuses the extend DP.
//!
//! File format (magic `JKMER01\0`, version 1, all little-endian):
//! `magic[8] version u32 gtf_sha256[32] fasta_sha256[32] n_junc u32
//! (junctions: contig,id,intron_start,intron_end u32, strand u8)
//! n_entries u32 (entry: packed u64, n_hits u32, hits: junction_id u32,
//! split_offset u8, a_mask u16) trailing sha256(body)[32]` where `body` is
//! every byte before the trailing hash. Same input ⇒ byte-identical file.

use crate::error::AlignError;
use crate::gtf;
use crate::seed::Strand;
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;

/// Junction k-mer length.
pub const K: usize = 15;
/// Flank length per junction side in the pseudo-sequence.
pub const FLANK: usize = 16;
/// Maximum A count for full A→G variant enumeration.
pub const MAX_A_FOR_VARIANTS: u32 = 6;

/// File magic.
pub const JKMER_MAGIC: [u8; 8] = *b"JKMER01\0";
/// File version.
pub const JKMER_VERSION: u32 = 1;

/// Minimum candidate hits on one junction (frozen).
pub const MIN_HITS: usize = 2;

/// A junction in the Track-2 system (NOT `gtf::Junction`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Junction {
    /// Contig index.
    pub contig: u32,
    /// Dense junction id (index into `JkmerIndex::junctions`).
    pub id: u32,
    /// Intron start (first intron base, plus-strand coordinate).
    pub intron_start: u32,
    /// Intron end (first base after the intron).
    pub intron_end: u32,
    /// Transcript strand.
    pub strand: Strand,
}

impl Junction {
    /// Intron length.
    pub fn intron_len(&self) -> u32 {
        self.intron_end - self.intron_start
    }
}

/// One k-mer hit registered under a packed k-mer.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct JkmerHit {
    /// Owning junction id.
    pub junction_id: u32,
    /// Breakpoint position inside the k-mer (donor-side base count).
    pub split_offset: u8,
    /// Bit `i` set iff k-mer base `i` (leftmost = 0) is `A`.
    pub a_mask: u16,
}

/// The junction-kmer index.
#[derive(Clone, Debug)]
pub struct JkmerIndex {
    /// File magic.
    pub magic: [u8; 8],
    /// File version.
    pub version: u32,
    /// sha256 of the source GTF file bytes.
    pub gtf_sha256: [u8; 32],
    /// sha256 of the source FASTA file bytes.
    pub fasta_sha256: [u8; 32],
    /// Junctions, dense by `id`.
    pub junctions: Vec<Junction>,
    /// Packed k-mer → hits.
    pub entries: BTreeMap<u64, Vec<JkmerHit>>,
}

/// 2-bit code of an ASCII base; `None` for anything else.
fn base_code(b: u8) -> Option<u64> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn code_base(c: u64) -> u8 {
    match c {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

/// Pack a K-length base string into 2 bits per base, high-order first.
pub fn pack_kmer(seq: &[u8]) -> u64 {
    let mut packed = 0u64;
    for &b in seq {
        packed = (packed << 2) | base_code(b).unwrap_or(0);
    }
    packed
}

/// Inverse of [`pack_kmer`] (always K bases).
pub fn unpack_kmer(packed: u64) -> [u8; K] {
    let mut out = [0u8; K];
    for (i, slot) in out.iter_mut().enumerate() {
        let shift = 2 * (K - 1 - i);
        *slot = code_base((packed >> shift) & 3);
    }
    out
}

/// Every K-window of `seq` without N, as `(read_pos, packed)`.
pub fn extract_kmers(seq: &[u8]) -> Vec<(u32, u64)> {
    if seq.len() < K {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(seq.len() - K + 1);
    for i in 0..=seq.len() - K {
        let win = &seq[i..i + K];
        if win.iter().all(|&b| base_code(b).is_some()) {
            out.push((i as u32, pack_kmer(win)));
        }
    }
    out
}

/// A-position mask of a k-mer (bit `i` = leftmost base position `i` is `A`).
pub fn compute_a_mask(seq: &[u8]) -> u16 {
    let mut mask = 0u16;
    for (i, &b) in seq.iter().enumerate().take(16) {
        if b.eq_ignore_ascii_case(&b'A') {
            mask |= 1 << i;
        }
    }
    mask
}

/// All A→G variants of a packed k-mer (subset enumeration, `2^n` including
/// the original; k-mers with more than `MAX_A_FOR_VARIANTS` A's yield only
/// the original).
pub fn enumerate_a_to_g_variants(packed: u64) -> Vec<u64> {
    let bases = unpack_kmer(packed);
    let a_pos: Vec<usize> = bases
        .iter()
        .enumerate()
        .filter(|(_, &b)| b.eq_ignore_ascii_case(&b'A'))
        .map(|(i, _)| i)
        .collect();
    if a_pos.len() as u32 > MAX_A_FOR_VARIANTS {
        return vec![packed];
    }
    let n = a_pos.len();
    let mut out = Vec::with_capacity(1usize << n);
    for subset in 0..(1u32 << n) {
        let mut v = packed;
        for (bit, &pos) in a_pos.iter().enumerate() {
            if subset >> bit & 1 == 1 {
                let shift = 2 * (K - 1 - pos);
                v = (v & !(3 << shift)) | (2 << shift); // A(0) -> G(2)
            }
        }
        out.push(v);
    }
    out
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match base_code(b) {
            Some(c) => code_base(3 - c),
            None => b,
        })
        .collect()
}

/// Fetch the donor tail and acceptor head flanks of a junction in
/// read-forward (transcript) orientation; minus-strand junctions
/// reverse-complement both (and swap the genomic sides).
///
/// `fetch_base(contig, start, end)` returns the requested reference slice and
/// may return fewer bases near contig boundaries.
pub fn fetch_flanks<F: Fn(u32, u32, u32) -> Vec<u8>>(
    junc: &Junction,
    fetch_base: F,
    donor: &mut Vec<u8>,
    acceptor: &mut Vec<u8>,
) {
    donor.clear();
    acceptor.clear();
    if junc.strand == Strand::Plus {
        donor.extend(fetch_base(
            junc.contig,
            junc.intron_start.saturating_sub(FLANK as u32),
            junc.intron_start,
        ));
        acceptor.extend(fetch_base(
            junc.contig,
            junc.intron_end,
            junc.intron_end + FLANK as u32,
        ));
    } else {
        donor.extend(revcomp(&fetch_base(
            junc.contig,
            junc.intron_end,
            junc.intron_end + FLANK as u32,
        )));
        acceptor.extend(revcomp(&fetch_base(
            junc.contig,
            junc.intron_start.saturating_sub(FLANK as u32),
            junc.intron_start,
        )));
    }
}

/// Build the k-mer entry table over junctions (crossing k-mers only, with
/// A→G variant expansion). Hits per k-mer are sorted by
/// (junction_id, split_offset, a_mask).
pub fn build_junction_kmers<F: Fn(u32, u32, u32) -> Vec<u8>>(
    junctions: &[Junction],
    fetch_base: F,
) -> BTreeMap<u64, Vec<JkmerHit>> {
    let mut entries: BTreeMap<u64, Vec<JkmerHit>> = BTreeMap::new();
    let mut donor: Vec<u8> = Vec::with_capacity(FLANK);
    let mut acceptor: Vec<u8> = Vec::with_capacity(FLANK);
    for junc in junctions {
        fetch_flanks(junc, &fetch_base, &mut donor, &mut acceptor);
        let dlen = donor.len();
        let mut pseudo = donor.clone();
        pseudo.extend_from_slice(&acceptor);
        if pseudo.len() < K {
            continue;
        }
        for p in 0..=pseudo.len() - K {
            // crossing k-mers only: the breakpoint (donor length) falls
            // strictly inside the window
            if p >= dlen || p + K <= dlen {
                continue;
            }
            let win = &pseudo[p..p + K];
            let packed = pack_kmer(win);
            let split_offset = (dlen - p) as u8;
            let a_mask = compute_a_mask(win);
            for var in enumerate_a_to_g_variants(packed) {
                entries.entry(var).or_default().push(JkmerHit {
                    junction_id: junc.id,
                    split_offset,
                    a_mask,
                });
            }
        }
    }
    for hits in entries.values_mut() {
        hits.sort_unstable_by_key(|h| (h.junction_id, h.split_offset, h.a_mask));
    }
    entries
}

/// Extract Track-2 junctions from a GTF: one per annotated intron, ids
/// assigned in library (sorted) order.
pub fn extract_junctions_from_gtf<F: Fn(&str) -> u32>(
    path: &Path,
    contig_id_fn: F,
) -> Result<Vec<Junction>, AlignError> {
    let lib = gtf::from_gtf(path, contig_id_fn)?;
    Ok(lib
        .junctions
        .iter()
        .enumerate()
        .map(|(id, j)| Junction {
            contig: j.contig,
            id: id as u32,
            intron_start: j.start,
            intron_end: j.end,
            strand: if j.minus_strand {
                Strand::Minus
            } else {
                Strand::Plus
            },
        })
        .collect())
}

impl JkmerIndex {
    /// Build from junctions and a reference fetch closure.
    pub fn build<F: Fn(u32, u32, u32) -> Vec<u8>>(
        junctions: Vec<Junction>,
        fetch_base: F,
        gtf_sha256: [u8; 32],
        fasta_sha256: [u8; 32],
    ) -> JkmerIndex {
        let entries = build_junction_kmers(&junctions, fetch_base);
        JkmerIndex {
            magic: JKMER_MAGIC,
            version: JKMER_VERSION,
            gtf_sha256,
            fasta_sha256,
            junctions,
            entries,
        }
    }

    /// Serialize (byte-deterministic; `BTreeMap` iteration is sorted).
    pub fn save<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        let mut body: Vec<u8> = Vec::new();
        body.extend_from_slice(&self.magic);
        body.extend_from_slice(&self.version.to_le_bytes());
        body.extend_from_slice(&self.gtf_sha256);
        body.extend_from_slice(&self.fasta_sha256);
        body.extend_from_slice(&(self.junctions.len() as u32).to_le_bytes());
        for j in &self.junctions {
            body.extend_from_slice(&j.contig.to_le_bytes());
            body.extend_from_slice(&j.id.to_le_bytes());
            body.extend_from_slice(&j.intron_start.to_le_bytes());
            body.extend_from_slice(&j.intron_end.to_le_bytes());
            body.extend_from_slice(&[if j.strand == Strand::Minus { 1u8 } else { 0u8 }]);
        }
        body.extend_from_slice(&(self.entries.len() as u32).to_le_bytes());
        for (packed, hits) in &self.entries {
            body.extend_from_slice(&packed.to_le_bytes());
            body.extend_from_slice(&(hits.len() as u32).to_le_bytes());
            for h in hits {
                body.extend_from_slice(&h.junction_id.to_le_bytes());
                body.push(h.split_offset);
                body.extend_from_slice(&h.a_mask.to_le_bytes());
            }
        }
        let sha: [u8; 32] = Sha256::digest(&body).into();
        w.write_all(&body)?;
        w.write_all(&sha)
    }

    /// Load and strictly validate (trailing sha256 verified before parsing;
    /// bad files error, never panic).
    pub fn load(path: &Path) -> Result<JkmerIndex, AlignError> {
        let bytes = std::fs::read(path).map_err(|_| AlignError::IndexIo)?;
        let min_len = 8 + 4 + 32 + 32 + 4 + 4 + 32;
        if bytes.len() < min_len {
            return Err(AlignError::IndexFormat {
                msg: format!("jkmer file too short: {} bytes", bytes.len()),
            });
        }
        let body = &bytes[..bytes.len() - 32];
        let trailer = &bytes[bytes.len() - 32..];
        let sha: [u8; 32] = Sha256::digest(body).into();
        if sha != trailer {
            return Err(AlignError::IndexFormat {
                msg: "jkmer trailing sha256 mismatch".to_string(),
            });
        }
        if body[0..8] != JKMER_MAGIC {
            return Err(AlignError::IndexFormat {
                msg: "jkmer bad magic".to_string(),
            });
        }
        let mut cur = 8usize;
        macro_rules! take {
            ($n:expr) => {{
                let n = $n;
                if cur + n > body.len() {
                    return Err(AlignError::IndexFormat {
                        msg: "jkmer truncated body".to_string(),
                    });
                }
                let s = &body[cur..cur + n];
                cur += n;
                s
            }};
        }
        let version = u32::from_le_bytes(take!(4).try_into().unwrap());
        if version != JKMER_VERSION {
            return Err(AlignError::IndexVersion {
                file: path.display().to_string(),
                supported: JKMER_VERSION,
            });
        }
        let gtf_sha256: [u8; 32] = take!(32).try_into().unwrap();
        let fasta_sha256: [u8; 32] = take!(32).try_into().unwrap();

        let n_junc = u32::from_le_bytes(take!(4).try_into().unwrap()) as usize;

        let n_junc_usize = n_junc;
        if cur
            + n_junc_usize
                .checked_mul(21)
                .ok_or_else(|| AlignError::IndexFormat {
                    msg: "jkmer junction count overflow".to_string(),
                })?
            > body.len()
        {
            return Err(AlignError::IndexFormat {
                msg: "jkmer truncated junction table".to_string(),
            });
        }
        let mut junctions = Vec::with_capacity(n_junc_usize);
        for _ in 0..n_junc_usize {
            let contig = u32::from_le_bytes(take!(4).try_into().unwrap());
            let id = u32::from_le_bytes(take!(4).try_into().unwrap());
            let intron_start = u32::from_le_bytes(take!(4).try_into().unwrap());
            let intron_end = u32::from_le_bytes(take!(4).try_into().unwrap());
            let strand_b = take!(1)[0];
            let strand = match strand_b {
                0 => Strand::Plus,
                1 => Strand::Minus,
                _ => {
                    return Err(AlignError::IndexFormat {
                        msg: format!("jkmer bad strand byte {strand_b}"),
                    })
                }
            };
            junctions.push(Junction {
                contig,
                id,
                intron_start,
                intron_end,
                strand,
            });
        }

        let n_entries = u32::from_le_bytes(take!(4).try_into().unwrap()) as usize;
        let mut entries = BTreeMap::new();
        for _ in 0..n_entries {
            let packed_s = take!(8);
            let packed = u64::from_le_bytes(packed_s.try_into().unwrap());
            let n_hits = u32::from_le_bytes(take!(4).try_into().unwrap()) as usize;
            let mut hits = Vec::with_capacity(n_hits.min(1024));
            for _ in 0..n_hits {
                let junction_id = u32::from_le_bytes(take!(4).try_into().unwrap());
                let split_offset = take!(1)[0];
                let a_mask_s = take!(2);
                let a_mask = u16::from_le_bytes(a_mask_s.try_into().unwrap());
                if junction_id as usize >= junctions.len() {
                    return Err(AlignError::IndexFormat {
                        msg: "jkmer hit references unknown junction".to_string(),
                    });
                }
                hits.push(JkmerHit {
                    junction_id,
                    split_offset,
                    a_mask,
                });
            }
            entries.insert(packed, hits);
        }
        if cur != body.len() {
            return Err(AlignError::IndexFormat {
                msg: "jkmer trailing bytes in body".to_string(),
            });
        }

        Ok(JkmerIndex {
            magic: JKMER_MAGIC,
            version: JKMER_VERSION,
            gtf_sha256,
            fasta_sha256,
            junctions,
            entries,
        })
    }

    /// Query a read: aggregate exact k-mer hits per junction, keep junctions
    /// with ≥2 hits, hits sorted by read position; candidates ordered by hit
    /// count descending then junction id ascending (deterministic).
    pub fn query_read(&self, read: &[u8]) -> Vec<JunctionCandidate<'_>> {
        let mut by_junction: BTreeMap<u32, Vec<JkmerHitInfo>> = BTreeMap::new();
        for (pos, packed) in extract_kmers(read) {
            if let Some(hits) = self.entries.get(&packed) {
                for h in hits {
                    by_junction
                        .entry(h.junction_id)
                        .or_default()
                        .push(JkmerHitInfo {
                            read_pos: pos,
                            split_offset: h.split_offset,
                            a_mask: h.a_mask,
                        });
                }
            }
        }
        let mut cands: Vec<JunctionCandidate<'_>> = Vec::new();
        for (jid, mut hits) in by_junction {
            if hits.len() < MIN_HITS {
                continue;
            }
            let Some(junction) = self.junctions.get(jid as usize) else {
                continue;
            };
            hits.sort_unstable_by_key(|h| h.read_pos);
            cands.push(JunctionCandidate { junction, hits });
        }
        cands.sort_by(|a, b| {
            b.hits
                .len()
                .cmp(&a.hits.len())
                .then(a.junction.id.cmp(&b.junction.id))
        });
        cands
    }
}

/// One aggregated hit on a junction (read coordinates).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct JkmerHitInfo {
    /// Read position of the hitting k-mer.
    pub read_pos: u32,
    /// Breakpoint offset inside that k-mer.
    pub split_offset: u8,
    /// A-mask of the registered k-mer.
    pub a_mask: u16,
}

/// One junction candidate for a queried read.
#[derive(Clone, Debug)]
pub struct JunctionCandidate<'a> {
    /// The junction (borrowed from the index).
    pub junction: &'a Junction,
    /// Aggregated hits, sorted by `read_pos`.
    pub hits: Vec<JkmerHitInfo>,
}

/// Breakpoint confidence.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Breakpoint {
    /// ≥2 votes within ±1bp.
    HighConf(u32),
    /// Fewer than 2 votes.
    LowConf(u32),
}

impl JunctionCandidate<'_> {
    /// Breakpoint estimate = `read_pos + split_offset` per hit; the mode over
    /// ±1bp windows wins (smallest position on ties). ≥2 votes ⇒ HighConf,
    /// otherwise LowConf; a negative mode yields `None`.
    pub fn infer_breakpoint(&self) -> Option<Breakpoint> {
        let mut est: Vec<i64> = self
            .hits
            .iter()
            .map(|h| h.read_pos as i64 + h.split_offset as i64)
            .collect();
        if est.is_empty() {
            return None;
        }
        est.sort_unstable();
        let mut best_pos = est[0];
        let mut best_votes = 0usize;
        for &v in &est {
            let votes = est.iter().filter(|&&e| (e - v).abs() <= 1).count();
            if votes > best_votes {
                best_votes = votes;
                best_pos = v;
            }
        }
        if best_pos < 0 {
            return None;
        }
        let pos = best_pos as u32;
        if best_votes >= 2 {
            Some(Breakpoint::HighConf(pos))
        } else {
            Some(Breakpoint::LowConf(pos))
        }
    }
}

/// Local confirmation scores (affine SW variant).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct LocalConfirm {
    /// Left segment vs donor tail.
    pub left: i32,
    /// Right segment vs acceptor head.
    pub right: i32,
    /// `left + right`.
    pub total: i32,
}

/// Confirm a read split locally: `left` = `local_align_score(read[..split],
/// donor_tail)` (0 when `split == 0`); `right` = `local_align_score(
/// read[split..], acceptor_head)` (0 when `split >= read.len()`).
pub fn local_confirm(
    read: &[u8],
    split: usize,
    donor_tail: &[u8],
    acceptor_head: &[u8],
) -> LocalConfirm {
    let left = if split == 0 {
        0
    } else {
        local_align_score(&read[..split.min(read.len())], donor_tail)
    };
    let right = if split >= read.len() {
        0
    } else {
        local_align_score(&read[split..], acceptor_head)
    };
    LocalConfirm {
        left,
        right,
        total: left + right,
    }
}

/// Affine Smith-Waterman variant: three matrices, `gap_open = −4`,
/// `gap_ext = −1`, H floored at 0; global maximum. Substitution:
/// match 5, `(read G, ref A)` and `(read C, ref T)` = 0 (editing-aware),
/// transitions (A↔G, C↔T) −1, everything else (incl. any N) −4.
pub fn local_align_score(read_seg: &[u8], ref_seg: &[u8]) -> i32 {
    let m = read_seg.len();
    let n = ref_seg.len();
    if m == 0 || n == 0 {
        return 0;
    }
    const GAP_OPEN: i32 = -4;
    const GAP_EXT: i32 = -1;
    let mut h_prev = vec![0i32; n + 1];
    let mut h_cur = vec![0i32; n + 1];
    let mut e_cur = vec![i32::MIN / 2; n + 1]; // gap over ref (horizontal)
    let mut best = 0i32;
    for i in 1..=m {
        h_cur[0] = 0;
        e_cur[0] = i32::MIN / 2;
        let mut f = i32::MIN / 2; // gap over read (vertical), per row
        for j in 1..=n {
            let diag = h_prev[j - 1] + sub_score(read_seg[i - 1], ref_seg[j - 1]);
            // E: gap consuming reference (from left)
            let e = (h_cur[j - 1] + GAP_OPEN).max(e_cur[j - 1] + GAP_EXT);
            // F: gap consuming read (from above)
            f = (h_prev[j] + GAP_OPEN).max(f + GAP_EXT);
            let h = diag.max(e).max(f).max(0);
            h_cur[j] = h;
            e_cur[j] = e;
            if h > best {
                best = h;
            }
        }
        std::mem::swap(&mut h_prev, &mut h_cur);
    }
    best
}

/// Substitution score of the Track-2 local alignment.
pub fn sub_score(read: u8, reference: u8) -> i32 {
    let r = read.to_ascii_uppercase();
    let f = reference.to_ascii_uppercase();
    let is_base = |b: u8| matches!(b, b'A' | b'C' | b'G' | b'T');
    if !is_base(r) || !is_base(f) {
        return -4;
    }
    if r == f {
        return 5;
    }
    if (r == b'G' && f == b'A') || (r == b'C' && f == b'T') {
        return 0;
    }
    let purine = |b: u8| b == b'A' || b == b'G';
    if purine(r) == purine(f) {
        -1 // transition
    } else {
        -4 // transversion
    }
}

/// Track-2 CIGAR op (M/N only).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CigarOpT2 {
    /// Aligned bases.
    Match(u32),
    /// Intron skip.
    RefSkip(u32),
}

/// `[M?] N [M?]` — zero-length segments emit nothing.
pub fn build_track2_cigar(left: u32, intron: u32, right: u32) -> Vec<CigarOpT2> {
    let mut out = Vec::with_capacity(3);
    if left > 0 {
        out.push(CigarOpT2::Match(left));
    }
    out.push(CigarOpT2::RefSkip(intron));
    if right > 0 {
        out.push(CigarOpT2::Match(right));
    }
    out
}

/// SAM-style string ("15M100N30M").
pub fn cigar_to_string(cigar: &[CigarOpT2]) -> String {
    let mut s = String::new();
    for op in cigar {
        match op {
            CigarOpT2::Match(n) => {
                s.push_str(&format!("{n}M"));
            }
            CigarOpT2::RefSkip(n) => {
                s.push_str(&format!("{n}N"));
            }
        }
    }
    s
}

/// Track-2 MAPQ: ties (`runner_up == best > 0`) ⇒ 0; otherwise
/// `log2(best)*10 + clamp(margin/10*5, 0..20) − min(ratio*20, 20)`, rounded
/// and clamped to 0..60 (`ratio = runner_up / best`). The pipeline does not
/// call this (see the spec quirk note); it is kept for parity tooling.
pub fn compute_track2_mapq(best_hits: u32, runner_up_hits: u32, score_margin: i32) -> u8 {
    if best_hits == 0 {
        return 0;
    }
    if runner_up_hits == best_hits {
        return 0;
    }
    let ratio = runner_up_hits as f64 / best_hits as f64;
    let margin_term = ((score_margin as f64) / 10.0 * 5.0).clamp(0.0, 20.0);
    let q = (best_hits as f64).log2() * 10.0 + margin_term - (ratio * 20.0).min(20.0);
    q.round().clamp(0.0, 60.0) as u8
}

/// An assembled Track-2 placement.
#[derive(Clone, Debug)]
pub struct Track2Record {
    /// 1-based leftmost reference position (saturating).
    pub pos_1based: u32,
    /// `[M?] N [M?]` CIGAR.
    pub cigar: Vec<CigarOpT2>,
    /// `compute_track2_mapq` of the same inputs (reference only).
    pub mapq: u8,
}

/// Assemble a Track-2 record: `left = read_split`, `right = read_len −
/// read_split`; plus strand pos = `intron_start − read_split + 1`, minus
/// strand pos = `intron_start − right_match + 1` (both saturating).
pub fn assemble_track2_record(
    junc: &Junction,
    read_len: u32,
    read_split: u32,
    best_hits: u32,
    runner_up: u32,
    margin: i32,
) -> Track2Record {
    let right_match = read_len.saturating_sub(read_split);
    let cigar = build_track2_cigar(read_split, junc.intron_len(), right_match);
    let pos_1based = match junc.strand {
        Strand::Plus => junc.intron_start.saturating_sub(read_split) + 1,
        Strand::Minus => junc.intron_start.saturating_sub(right_match) + 1,
    };
    Track2Record {
        pos_1based,
        cigar,
        mapq: compute_track2_mapq(best_hits, runner_up, margin),
    }
}
