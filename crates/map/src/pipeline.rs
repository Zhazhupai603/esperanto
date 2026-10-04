//! Alignment pipeline: chunked streaming FASTQ → rayon align → ordered write of
//! BAM / unmapped.fq.gz / align_qc.json.
//!
//! Engineering invariants:
//! - Determinism: parallel within a chunk, merged by chunk order; output record
//!   order == input order (byte-identical at any thread count).
//! - Constant memory: CHUNK pairs per block; stats and output advance per chunk.
//! - unmapped.fq.gz hard contract: every unplaced read (incl. single end of a
//!   PE pair) is written, name keeps the /1 /2 mate marker.
//! - BAM unmapped records written with FLAG 0x4 (+ mate flags).

use std::io::Write;
use std::path::Path;
use std::sync::Arc;

use rayon::prelude::*;

use crate::align::{AlignConfig, Aligner, CachedPlacements};
use crate::bam::{self, flag, BamRecord};
use crate::error::AlignError;
use crate::extend::CigarOp;
use crate::fastq::{FastqReader, FastqRecord, RecordSource};
use crate::gtf::{Junction, JunctionLib};
use crate::index::Index;
use crate::mapq::{mapq_of, ReadAlignment, PILEUP_GATE};
use crate::pair::{relate, template_length, PairRelation};
use crate::seed::Strand;
use crate::stats::{AlignStats, StatsAcc};

/// Streaming block size (pairs / reads).
pub const CHUNK: usize = 100_000;

/// Pipeline output sinks (bam sink ownership is taken at run time —
/// multithreaded BGZF requires 'static).
pub struct PipelineOut<'a> {
    pub bam: Option<Box<dyn Write + Send>>,
    pub unmapped_fq: Box<dyn Write + Send>,
    pub index: &'a Index,
    pub config: AlignConfig,
    /// Junction library (sjdb; RNA mode).
    pub jlib: Option<Arc<JunctionLib>>,
    /// L1 transcriptome-first engine bundle (None = genomic alignment only).
    pub l1: Option<Arc<esperanto_engine::L1Index>>,
}

/// pass1 junction discoveries: (junction, support, summed MAPQ of supporters).
pub type Discoveries = Vec<(Junction, u32, u64)>;

/// What pass 1 decided for one filtered discovery.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum DiscoveryFate {
    /// No annotated junction within 3 bp: genuinely novel, merged.
    Kept,
    /// Within 3 bp of an annotated junction; pass-1 support within a 2:1 band
    /// of the annotated boundary's own observation count. Both stay merged.
    TieKept,
    /// Within 3 bp of an annotated junction and observed >2x as often:
    /// merged, and the annotated copy leaves the merged library (the sample
    /// does not express it as annotated).
    OutvotesAnnotated,
    /// Within 3 bp of an annotated junction observed >2x as often: the
    /// curated boundary wins; reported but not merged.
    OutvotedByAnnotated,
}

/// A pass-1 discovery that passed the quality gates: (junction, support, fate).
pub type FilteredDiscovery = (Junction, u32, DiscoveryFate);

/// 2-pass outputs: (stats, pass1 filtered discoveries, final library list,
/// count of annotated junctions removed by the support vote).
pub type TwoPassOut = (
    AlignStats,
    Vec<FilteredDiscovery>,
    Vec<(Junction, u32)>,
    usize,
);

fn make_aligner<'a>(
    index: &'a Index,
    cfg: AlignConfig,
    lib: &Option<Arc<JunctionLib>>,
    l1: &Option<Arc<esperanto_engine::L1Index>>,
) -> Aligner<'a> {
    let a = Aligner::new(index, cfg);
    let a = match lib {
        Some(l) => a.with_lib(l.clone()),
        None => a,
    };
    match l1 {
        Some(l1) => a.with_l1(l1.clone()),
        None => a,
    }
}

/// Alignment worker stack size. On repeat-dense reads the aligner's
/// extension/arbitration chain comes close to rayon's 2 MiB default worker
/// stack, and the crossbeam epoch collector riding on the same stack then
/// tips it over (observed as an abort at 128 threads on the full real
/// sample). 64 MiB is reserved address space only — committed lazily.
const WORKER_STACK: usize = 64 << 20;

/// Single-end pipeline (returns stats + pass1 discoveries, empty in DNA mode).
pub fn run_se(
    out: &mut PipelineOut,
    r1: &Path,
    threads: usize,
) -> Result<(AlignStats, Discoveries), AlignError> {
    let t0 = std::time::Instant::now();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .stack_size(WORKER_STACK)
        .build()
        .map_err(|e| AlignError::FastqFormat {
            line: 0,
            msg: format!("thread pool: {e}"),
        })?;
    let mut sinks =
        crate::sinks::RecordSinks::new(out.bam.take(), out.index, threads)?;
    let mut unm = crate::sinks::UnmappedFq::new(&mut out.unmapped_fq);
    let mut stats = StatsAcc::new();
    let mut discoveries: Discoveries = Vec::new();

    let mut reader = FastqReader::open(r1)?;

    loop {
        let mut chunk = Vec::with_capacity(CHUNK);
        for _ in 0..CHUNK {
            match reader.next_record()? {
                Some(r) => chunk.push(r),
                None => break,
            }
        }
        if chunk.is_empty() {
            break;
        }
        let index = out.index;
        let cfg = out.config;
        let lib = out.jlib.clone();
        let l1_clone = out.l1.clone();
        // align_read(seq) depends only on seq, not qual/name.
        // Within-chunk dedup: align each unique seq once, fan out to all reads
        // sharing it (each keeps its own name/qual). Deterministic algorithm +
        // ordered write ⇒ byte-identical.
        let (unique_seqs, uid_of_read) = dedup_se(&chunk);
        let cached_alns: Vec<Option<ReadAlignment>> = {
            let miss_seqs: Vec<&[u8]> = unique_seqs.to_vec();
            pool.install(|| {
                let lib = lib.clone();
                miss_seqs
                    .par_iter()
                    .map_init(
                        move || make_aligner(index, cfg, &lib, &l1_clone),
                        |aligner, &seq| aligner.align_read(seq),
                    )
                    .collect()
            })
        };
        // Build BamRecord per read in input order, cloning the cached
        // alignment and pairing it with the read's own name/qual.
        let mut results: Vec<BamRecord> = Vec::with_capacity(chunk.len());
        for (i, rec) in chunk.iter().enumerate() {
            results.push(bam::record_se(
                &rec.name,
                &rec.seq,
                &rec.qual,
                cached_alns[uid_of_read[i]].as_ref(),
            ));
        }
        for (i, (rec, br)) in chunk.iter().zip(&results).enumerate() {
            let aln = &cached_alns[uid_of_read[i]];
            stats.push_read(aln.is_some(), br.mapq);
            if let Some(x) = aln {
                if x.origin == crate::mapq::Origin::Decoy {
                    stats.push_decoy();
                }
                if x.rescued {
                    stats.push_rescued();
                }
                stats.push_junctions(x.junctions.len() as u64);
                for j in &x.junctions {
                    discoveries.push((j.junction, 1, u64::from(br.mapq)));
                }
            }
            if aln.is_none() {
                unm.write(rec, None)?;
            }
            sinks.write(br)?;
        }
    }
    sinks.finish()?;
    unm.finish()?;
    let mode = if out.config.rna { "rna-se" } else { "dna-se" };
    Ok((
        stats.finalize(mode, false, t0.elapsed().as_secs_f64()),
        discoveries,
    ))
}

/// Single-end deduplication: group reads in a chunk by identical seq → unique seq
/// slices (by uid) + per-read → uid map.
fn dedup_se(chunk: &[FastqRecord]) -> (Vec<&[u8]>, Vec<usize>) {
    let mut seq_to_uid: std::collections::HashMap<&[u8], usize> =
        std::collections::HashMap::with_capacity(chunk.len());
    let mut unique_seqs: Vec<&[u8]> = Vec::new();
    let mut uid_of_read: Vec<usize> = Vec::with_capacity(chunk.len());
    for r in chunk {
        let seq: &[u8] = r.seq.as_slice();
        let uid = *seq_to_uid.entry(seq).or_insert_with(|| {
            let u = unique_seqs.len();
            unique_seqs.push(seq);
            u
        });
        uid_of_read.push(uid);
    }
    (unique_seqs, uid_of_read)
}

/// Paired-end deduplication: group both mates of a chunk by identical seq → unique
/// seq slices (by uid, first-seen order) + per-pair → uid maps for R1/R2.
/// One unified table: a sequence appearing as both an R1 and an R2 aligns once.
fn dedup_pe(chunk: &[(FastqRecord, FastqRecord)]) -> (Vec<&[u8]>, Vec<usize>, Vec<usize>) {
    let mut seq_to_uid: std::collections::HashMap<&[u8], usize> =
        std::collections::HashMap::with_capacity(2 * chunk.len());
    let mut unique_seqs: Vec<&[u8]> = Vec::new();
    let mut uid1_of_pair: Vec<usize> = Vec::with_capacity(chunk.len());
    let mut uid2_of_pair: Vec<usize> = Vec::with_capacity(chunk.len());
    for (a, b) in chunk {
        uid1_of_pair.push(intern_seq(&mut seq_to_uid, &mut unique_seqs, &a.seq));
        uid2_of_pair.push(intern_seq(&mut seq_to_uid, &mut unique_seqs, &b.seq));
    }
    (unique_seqs, uid1_of_pair, uid2_of_pair)
}

/// Intern one sequence into the uid table; uid = first-seen insertion order.
fn intern_seq<'a>(
    seq_to_uid: &mut std::collections::HashMap<&'a [u8], usize>,
    unique_seqs: &mut Vec<&'a [u8]>,
    seq: &'a [u8],
) -> usize {
    let next_uid = unique_seqs.len();
    *seq_to_uid.entry(seq).or_insert_with(|| {
        unique_seqs.push(seq);
        next_uid
    })
}

/// Paired-end pipeline (returns stats + pass1 discoveries, empty in DNA mode).
pub fn run_pe(
    out: &mut PipelineOut,
    r1: &Path,
    r2: &Path,
    threads: usize,
) -> Result<(AlignStats, Discoveries), AlignError> {
    let t0 = std::time::Instant::now();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .stack_size(WORKER_STACK)
        .build()
        .map_err(|e| AlignError::FastqFormat {
            line: 0,
            msg: format!("thread pool: {e}"),
        })?;
    let mut sinks =
        crate::sinks::RecordSinks::new(out.bam.take(), out.index, threads)?;
    let mut unm = crate::sinks::UnmappedFq::new(&mut out.unmapped_fq);
    let mut stats = StatsAcc::new();
    let mut discoveries: Discoveries = Vec::new();
    let mut rd1 = FastqReader::open(r1)?;
    let mut rd2 = FastqReader::open(r2)?;

    loop {
        let mut chunk: Vec<(FastqRecord, FastqRecord)> = Vec::with_capacity(CHUNK);
        for _ in 0..CHUNK {
            match (rd1.next_record()?, rd2.next_record()?) {
                (Some(a), Some(b)) => chunk.push((a, b)),
                (None, None) => break,
                // Pair sync: either end hitting EOF first is an error.
                (a, b) => {
                    return Err(AlignError::FastqFormat {
                        line: 0,
                        msg: format!(
                            "R1/R2 record count mismatch (r1 has more: {}, r2 has more: {})",
                            a.is_some(),
                            b.is_some()
                        ),
                    });
                }
            }
        }
        if chunk.is_empty() {
            break;
        }
        let index = out.index;
        let cfg = out.config;
        let lib = out.jlib.clone();
        let l1_clone = out.l1.clone();
        // Paired-end deduplication — two layers, chunk-local only (the table dies
        // with the chunk to keep memory bounded; no cross-chunk cache):
        //   1. primary alignment per UNIQUE sequence. Full seq bytes are the
        //      key: name/quality never participate in primary alignment
        //      today — if alignment ever consumes quality, this dedup must
        //      be revisited (the key must then include quality).
        //   2. per-pair mate-rescue arbitration is NEVER cached (rescue
        //      depends on the concrete pair): each pair clones the cached
        //      primary of its two sequences and runs the existing rescue
        //      logic unchanged.
        let (unique_seqs, uid1_of_pair, uid2_of_pair) = dedup_pe(&chunk);
        let unique_alns: Vec<(Option<ReadAlignment>, Option<CachedPlacements>)> = pool.install(|| {
            let lib = lib.clone();
            let l1_clone = l1_clone.clone();
            unique_seqs
                .par_iter()
                .map_init(
                    move || make_aligner(index, cfg, &lib, &l1_clone),
                    |aligner, &seq| aligner.align_read_full(seq),
                )
                .collect()
        });
        // Parallel per-pair rescue arbitration reads the uid-indexed primary
        // results (order preserved by collect); BAM records are built in the
        // serial write section.
        let results: Vec<(Option<ReadAlignment>, Option<ReadAlignment>)> = pool.install(|| {
            let lib = lib.clone();
            let l1_clone = l1_clone.clone();
            chunk
                .par_iter()
                .enumerate()
                .map_init(
                    move || make_aligner(index, cfg, &lib, &l1_clone),
                    |aligner, (i, (a, b))| {
                        let (base1, cache1) = &unique_alns[uid1_of_pair[i]];
                        let (base2, cache2) = &unique_alns[uid2_of_pair[i]];
                        let mut r1 = base1.clone();
                        let mut r2 = base2.clone();
                        // Pair-anchored arbitration: a mate whose own
                        // placement is poor (unplaced, mostly clipped, or
                        // below the pileup gate) is re-arbitrated with the
                        // locus its confidently placed mate points at as an
                        // extra candidate (see `Aligner::align_read_with_mate`).
                        let confident = |r: &Option<ReadAlignment>| {
                            r.as_ref()
                                .is_some_and(|m| !m.rescued && mapq_of(m) >= PILEUP_GATE)
                        };
                        if placement_poor(&r1, a.seq.len()) && confident(&r2) {
                            let (mc, mp) = r2.as_ref().map(|m| (m.contig, m.pos)).unwrap();
                            if let Some(n) =
                                aligner.align_read_with_mate(&a.seq, mc, mp, cache1.as_ref())
                            {
                                r1 = Some(n);
                            }
                        }
                        if placement_poor(&r2, b.seq.len()) && confident(&r1) {
                            let (mc, mp) = r1.as_ref().map(|m| (m.contig, m.pos)).unwrap();
                            if let Some(n) =
                                aligner.align_read_with_mate(&b.seq, mc, mp, cache2.as_ref())
                            {
                                r2 = Some(n);
                            }
                        }
                        (r1, r2)
                    },
                )
                .collect()
        });
        for ((a, b), (aln1, aln2)) in chunk.iter().zip(&results) {
            // Records carry every number the stats/discovery section needs
            // (MAPQ, proper-pair flag, tlen) — one derivation per field,
            // no re-computation from the alignment.
            let (ra, rb) = pair_records(
                &a.name, &a.seq, &a.qual, aln1,
                &b.name, &b.seq, &b.qual, aln2,
            );
            if let Some(x) = aln1 {
                stats.push_junctions(x.junctions.len() as u64);
                for j in &x.junctions {
                    discoveries.push((j.junction, 1, u64::from(ra.mapq)));
                }
            }
            if let Some(x) = aln2 {
                stats.push_junctions(x.junctions.len() as u64);
                for j in &x.junctions {
                    discoveries.push((j.junction, 1, u64::from(rb.mapq)));
                }
            }
            stats.push_read(aln1.is_some(), ra.mapq);
            stats.push_read(aln2.is_some(), rb.mapq);
            if let Some(x) = aln1 {
                if x.origin == crate::mapq::Origin::Decoy {
                    stats.push_decoy();
                }
                if x.rescued {
                    stats.push_rescued();
                }
            }
            if let Some(x) = aln2 {
                if x.origin == crate::mapq::Origin::Decoy {
                    stats.push_decoy();
                }
                if x.rescued {
                    stats.push_rescued();
                }
            }
            if aln1.is_none() {
                if aln2.is_some() {
                    stats.push_rescue_fail();
                }
                unm.write(a, Some(1))?;
            }
            if aln2.is_none() {
                if aln1.is_some() {
                    stats.push_rescue_fail();
                }
                unm.write(b, Some(2))?;
            }
            // proper ⇔ both placed and relate()==Proper (mk set the flag from
            // exactly that condition); tlen is the value mk stored for R1.
            if ra.flag & flag::PROPER_PAIR != 0 {
                stats.push_proper_pair();
                if let Some((_, _, tlen)) = ra.mate {
                    stats.push_insert(tlen);
                }
            }
            sinks.write(&ra)?;
            sinks.write(&rb)?;
        }
    }
    sinks.finish()?;
    unm.finish()?;
    let mode = if out.config.rna { "rna-pe" } else { "dna-pe" };
    Ok((
        stats.finalize(mode, true, t0.elapsed().as_secs_f64()),
        discoveries,
    ))
}

/// Build PE paired records (FLAG/mate fields per spec).
///
/// Alignments are borrowed: no per-pair `ReadAlignment` clone — the records
/// copy out every field they need.
#[allow(clippy::too_many_arguments)]
fn pair_records(
    n1: &[u8],
    s1: &[u8],
    q1: &[u8],
    a1: &Option<ReadAlignment>,
    n2: &[u8],
    s2: &[u8],
    q2: &[u8],
    a2: &Option<ReadAlignment>,
) -> (BamRecord, BamRecord) {
    let proper = match (a1, a2) {
        (Some(x), Some(y)) => relate(x, y) == PairRelation::Proper,
        _ => false,
    };
    let tlen = match (a1, a2) {
        (Some(x), Some(y)) if proper => template_length(x, y),
        _ => 0,
    };
    let mk = |name: &[u8],
              seq: &[u8],
              qual: &[u8],
              aln: &Option<ReadAlignment>,
              mate_aln: &Option<ReadAlignment>,
              self_bit: u16,
              tlen: i32| {
        let mut f = flag::PAIRED | self_bit;
        if proper {
            f |= flag::PROPER_PAIR;
        }
        let (mate_c, mate_p) = match mate_aln {
            Some(m) => {
                if m.strand == Strand::Minus {
                    f |= flag::MATE_REVERSE;
                }
                (m.contig as i32, m.pos as i32)
            }
            None => {
                f |= flag::MATE_UNMAPPED;
                (-1, -1)
            }
        };
        let q = match aln {
            Some(a) => {
                if a.strand == Strand::Minus {
                    f |= flag::REVERSE;
                }
                mapq_of(a)
            }
            None => {
                f |= flag::UNMAPPED;
                0
            }
        };
        // Reference-forward SEQ direction (same as record_se): minus strand stores the
        // reference-forward SEQ (revcomp) with reversed QUAL; unmapped keeps
        // the original orientation.
        let (s_out, q_out) = crate::bam::t13_phred(
            aln.as_ref().is_some_and(|a| a.strand == Strand::Minus),
            seq,
            qual,
        );
        BamRecord {
            name: name.to_vec(),
            flag: f,
            aln: aln.as_ref().map(crate::bam::aln_view),
            mapq: q,
            seq: s_out,
            qual: q_out,
            mate: Some((mate_c, mate_p, tlen)),
        }
    };
    (
        mk(n1, s1, q1, a1, a2, flag::READ1, tlen),
        mk(n2, s2, q2, a2, a1, flag::READ2, -tlen),
    )
}


/// SE 2-pass (same semantics as the PE version). Runs pass-1 discovery even
/// with the L1 engine on: the annotation covers only ~75% of the junctions
/// expressed in a typical sample, and the L1 path cannot place a read across
/// a junction neither the bundle nor the library knows.
pub fn run_se_2pass(
    out: &mut PipelineOut,
    r1: &Path,
    threads: usize,
) -> Result<TwoPassOut, AlignError> {
    let base_lib = out.jlib.take();
    out.jlib = base_lib.clone();
    let real_bam = out.bam.take();
    out.bam = Some(Box::new(std::io::sink()));
    let real_unm = std::mem::replace(&mut out.unmapped_fq, Box::new(std::io::sink()));
    let pass1 = run_se(out, r1, threads);
    let sink_unm = std::mem::replace(&mut out.unmapped_fq, real_unm);
    drop(sink_unm);
    out.bam = real_bam;
    let (_, discovered) = pass1?;
    // Merge: full sjdb + discoveries with ≥2 support.
    let mut merged: Vec<(Junction, u32)> = match &base_lib {
        Some(lib) => lib
            .junctions
            .iter()
            .copied()
            .zip(lib.counts.iter().copied())
            .collect(),
        None => Vec::new(),
    };
    let (kept, removed_annotated) = filter_discoveries(discovered, base_lib.as_deref());
    if !removed_annotated.is_empty() {
        let removed: std::collections::BTreeSet<&Junction> = removed_annotated.iter().collect();
        merged.retain(|(j, _)| !removed.contains(j));
    }
    merged.extend(
        kept.iter()
            .filter(|(_, _, fate)| *fate != DiscoveryFate::OutvotedByAnnotated)
            .map(|(j, c, _)| (*j, *c)),
    );
    let final_lib = JunctionLib::build_with_counts(merged);
    let final_list: Vec<(Junction, u32)> = final_lib
        .junctions
        .iter()
        .copied()
        .zip(final_lib.counts.iter().copied())
        .collect();
    out.jlib = Some(Arc::new(final_lib));
    let (stats, _) = run_se(out, r1, threads)?;
    Ok((stats, kept, final_list, removed_annotated.len()))
}

/// RNA 2-pass: pass1 discovers junctions (>=2 support into the library) →
/// merge with sjdb → pass2 aligns with the merged library. Runs even with
/// the L1 engine on: the annotation covers only ~75% of the junctions
/// expressed in a typical sample, and the L1 path cannot place a read across
/// a junction neither the bundle nor the library knows.
pub fn run_pe_2pass(
    out: &mut PipelineOut,
    r1: &Path,
    r2: &Path,
    threads: usize,
) -> Result<TwoPassOut, AlignError> {
    // pass1: discovery round (BAM/unmapped discarded; real artifacts only
    // come from pass2).
    let base_lib = out.jlib.take();
    out.jlib = base_lib.clone();
    let real_bam = out.bam.take();
    out.bam = Some(Box::new(std::io::sink()));
    let real_unm = std::mem::replace(&mut out.unmapped_fq, Box::new(std::io::sink()));
    let pass1 = run_pe(out, r1, r2, threads);
    let sink_unm = std::mem::replace(&mut out.unmapped_fq, real_unm);
    drop(sink_unm);
    out.bam = real_bam;
    let (_, discovered) = pass1?;
    let mut merged: Vec<(Junction, u32)> = match &base_lib {
        Some(lib) => lib
            .junctions
            .iter()
            .copied()
            .zip(lib.counts.iter().copied())
            .collect(),
        None => Vec::new(),
    };
    let (kept, removed_annotated) = filter_discoveries(discovered, base_lib.as_deref());
    if !removed_annotated.is_empty() {
        let removed: std::collections::BTreeSet<&Junction> = removed_annotated.iter().collect();
        merged.retain(|(j, _)| !removed.contains(j));
    }
    merged.extend(
        kept.iter()
            .filter(|(_, _, fate)| *fate != DiscoveryFate::OutvotedByAnnotated)
            .map(|(j, c, _)| (*j, *c)),
    );
    let final_lib = JunctionLib::build_with_counts(merged);
    let final_list: Vec<(Junction, u32)> = final_lib
        .junctions
        .iter()
        .copied()
        .zip(final_lib.counts.iter().copied())
        .collect();
    out.jlib = Some(Arc::new(final_lib));
    // pass2: production alignment with the merged library.
    let (stats, _) = run_pe(out, r1, r2, threads)?;
    Ok((stats, kept, final_list, removed_annotated.len()))
}

/// Junction-level filter for 2-pass: support >= 2 AND mean supporter MAPQ >= 20
/// AND span <= 500kb (low-quality cross-repeat-family junctions stay out).
///
/// A surviving discovery within 3 bp of an annotated junction (same contig,
/// either strand) is a shifted re-observation of a curated boundary. Keeping
/// both creates a same-score tie that flips boundaries at scale; always
/// dropping the discovery loses real signal where the annotation is the side
/// the sample does not express. So the pass-1 observation counts vote: the
/// side observed >2x as often stays and the other leaves the merged library;
/// inside a 1:2..2:1 band both stay and arbitration decides. The annotated
/// side's pass-1 support is its own exact re-observation count in the
/// discovery stream (0 when never observed). The strand check is deliberately
/// blind: a shifted copy can surface under the opposite strand label
/// (observed on sim1).
pub fn filter_discoveries(
    discovered: Discoveries,
    annotated: Option<&JunctionLib>,
) -> (Vec<FilteredDiscovery>, Vec<Junction>) {
    let mut acc: std::collections::BTreeMap<Junction, (u32, u64)> =
        std::collections::BTreeMap::new();
    for (j, c, mq) in discovered {
        let e = acc.entry(j).or_default();
        e.0 += c;
        e.1 += mq;
    }
    let mut out = Vec::new();
    let mut removed: std::collections::BTreeSet<Junction> = std::collections::BTreeSet::new();
    for (j, (c, mqsum)) in &acc {
        if !(*c >= 2 && *mqsum >= u64::from(*c) * 20 && j.end.saturating_sub(j.start) <= 500_000) {
            continue;
        }
        let fate = match annotated {
            None => DiscoveryFate::Kept,
            Some(lib) => {
                let (js, _) = lib.range_start(
                    j.contig,
                    j.start.saturating_sub(3),
                    j.start.saturating_add(4),
                );
                let nbrs: Vec<&Junction> = js
                    .iter()
                    .filter(|a| a.start.abs_diff(j.start) <= 3 && a.end.abs_diff(j.end) <= 3)
                    .collect();
                if nbrs.is_empty() {
                    DiscoveryFate::Kept
                } else {
                    let sup_a = nbrs
                        .iter()
                        .map(|a| acc.get(a).map_or(0, |(s, _)| *s))
                        .max()
                        .unwrap_or(0);
                    if *c > 2 * sup_a {
                        for a in nbrs {
                            if *c > 2 * acc.get(a).map_or(0, |(s, _)| *s) {
                                removed.insert(*a);
                            }
                        }
                        DiscoveryFate::OutvotesAnnotated
                    } else if sup_a > 2 * *c {
                        DiscoveryFate::OutvotedByAnnotated
                    } else {
                        DiscoveryFate::TieKept
                    }
                }
            }
        };
        out.push((*j, *c, fate));
    }
    (out, removed.into_iter().collect())
}

/// A placement the pair-anchored arbitration should second-guess: absent,
/// explaining less than 90% of the read, or below the pileup gate.
fn placement_poor(r: &Option<ReadAlignment>, read_len: usize) -> bool {
    match r {
        None => true,
        Some(a) => {
            let aligned: u32 = a
                .cigar
                .iter()
                .map(|op| match op {
                    CigarOp::Match(n) | CigarOp::Ins(n) => *n,
                    _ => 0,
                })
                .sum();
            (aligned as f64) < 0.9 * read_len as f64 || mapq_of(a) < PILEUP_GATE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// pair_records field equality after the borrow-only refactor: flags,
    /// mate coordinates, tlen sign and MAPQ must match the relation
    /// functions exactly (historical construction), for a mapped pair and a
    /// mate-unmapped pair.
    #[test]
    fn pair_records_flags_mates_tlen_match_relation() {
        use crate::extend::CigarOp;
        let aln = |strand: Strand, pos: u32| ReadAlignment {
            contig: 0,
            pos,
            strand,
            score: 90,
            second_score: 0,
            chain_score: 90,
            second_chain_score: 0,
            cigar: vec![CigarOp::Match(50)],
            n_anchors: 10,
            junctions: Vec::new(),
            ea_ag: 0,
            ea_tc: 0,
            mm_count: 0,
            n_seeds: 10,
            rescued: false,
            mapq_override: None,
            origin: crate::mapq::Origin::Genomic,
            n_cands: 1,
        };
        let seq = b"ACGT".repeat(13); // 52
        let qual = vec![b'I'; 52];
        let a1 = aln(Strand::Plus, 100);
        let a2 = aln(Strand::Minus, 400);
        let (ra, rb) = pair_records(
            b"t/1", &seq, &qual, &Some(a1.clone()),
            b"t/2", &seq, &qual, &Some(a2.clone()),
        );
        let proper = relate(&a1, &a2) == PairRelation::Proper;
        let expected_tlen = if proper { template_length(&a1, &a2) } else { 0 };
        assert_eq!(ra.name, b"t/1");
        assert_eq!(rb.name, b"t/2");
        assert!(ra.flag & flag::PAIRED != 0);
        assert!(ra.flag & flag::READ1 != 0);
        assert!(rb.flag & flag::READ2 != 0);
        assert_eq!(ra.flag & flag::PROPER_PAIR != 0, proper);
        assert_eq!(rb.flag & flag::PROPER_PAIR != 0, proper);
        assert_eq!(ra.mate, Some((0, 400, expected_tlen)));
        assert_eq!(rb.mate, Some((0, 100, -expected_tlen)));
        // a2 minus: MATE_REVERSE on R1, REVERSE on R2; R1 plus keeps orientation
        assert!(ra.flag & flag::MATE_REVERSE != 0);
        assert!(rb.flag & flag::REVERSE != 0);
        assert_eq!(ra.flag & flag::REVERSE, 0);
        assert_eq!(rb.seq, esperanto_bamio::revcomp(&seq));
        assert_eq!(ra.seq, seq);
        assert_eq!(ra.mapq, mapq_of(&a1));
        assert_eq!(rb.mapq, mapq_of(&a2));

        // mate unmapped: MATE_UNMAPPED on the partner, no mate coords, tlen 0
        let (rc, rd) = pair_records(
            b"t/1", &seq, &qual, &Some(a1),
            b"t/2", &seq, &qual, &None,
        );
        assert!(rc.flag & flag::MATE_UNMAPPED != 0);
        assert_eq!(rc.mate, Some((-1, -1, 0)));
        assert!(rd.flag & flag::UNMAPPED != 0);
        assert_eq!(rd.mapq, 0);
        assert!(rd.aln.is_none());
    }

    /// L1-on must not skip pass-1 junction discovery: the annotation covers
    /// only ~75% of the junctions a sample expresses, so a novel junction is
    /// found by pass 1 and merged into the library for pass 2. Builds a tiny
    /// world (one 20 kb contig, a one-transcript GTF for the L1 bundle, three
    /// pairs whose R1 spans an UNANNOTATED GT..AG intron) and asserts the
    /// discovery is returned in `kept`.
    #[test]
    fn l1_on_still_discovers_novel_junction_in_pass1() {
        // deterministic pseudo-random reference, no repeats
        let n = 20_000usize;
        let mut seq = Vec::with_capacity(n);
        let mut state: u64 = 0x9E3779B97F4A7C15;
        for _ in 0..n {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            seq.push(b"ACGT"[(state >> 33) as usize & 3]);
        }
        // novel intron [8000, 9000): GT donor, AG acceptor (plus strand)
        seq[8000] = b'G'; seq[8001] = b'T';
        seq[8998] = b'A'; seq[8999] = b'G';
        let dir = std::env::temp_dir().join(format!("esp-2pass-l1-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("ref.fa");
        let mut fasta_txt = b">chrT\n".to_vec();
        for chunk in seq.chunks(60) {
            fasta_txt.extend_from_slice(chunk);
            fasta_txt.push(b'\n');
        }
        std::fs::write(&fa, &fasta_txt).unwrap();
        std::fs::write(dir.join("ref.fa.fai"), b"chrT\t20000\t6\t60\t61\n").unwrap();
        let gtf = dir.join("anno.gtf");
        std::fs::write(
            &gtf,
            b"chrT\tt\ttranscript\t12001\t16000\t.\t+\t.\ttranscript_id \"txA\";\n\
              chrT\tt\texon\t12001\t13000\t.\t+\t.\ttranscript_id \"txA\";\n\
              chrT\tt\texon\t15001\t16000\t.\t+\t.\ttranscript_id \"txA\";\n",
        )
        .unwrap();
        let tidx = dir.join("ref.tidx");
        esperanto_tidx::build(&gtf, &fa, &tidx, &esperanto_tidx::BuildOptions::default()).unwrap();
        let l1 = esperanto_engine::L1Index::build(&tidx, &gtf, &fa).unwrap();

        let reference = crate::fasta::parse_fasta_bytes(&fasta_txt).unwrap();
        let index = crate::index::Index::build(reference, crate::seed::SeedParams { k: 15, w: 5 });

        // three PE pairs; R1 spans the novel junction 7925..8000 + 9000..9075
        let mut r1_txt = Vec::new();
        let mut r2_txt = Vec::new();
        for i in 0..3 {
            let mut span = seq[7925..8000].to_vec();
            span.extend_from_slice(&seq[9000..9075]);
            let mate = &seq[9100 + i * 200..9250 + i * 200];
            let q = b"I".repeat(150);
            r1_txt.extend_from_slice(&[&format!("@novel{i}/1\n").into_bytes()[..], &span, b"\n+\n", &q, b"\n"].concat());
            r2_txt.extend_from_slice(&[&format!("@novel{i}/2\n").into_bytes()[..], mate, b"\n+\n", &q, b"\n"].concat());
        }
        let r1 = dir.join("r1.fq");
        let r2 = dir.join("r2.fq");
        std::fs::write(&r1, &r1_txt).unwrap();
        std::fs::write(&r2, &r2_txt).unwrap();

        let mut out = PipelineOut {
            bam: None,
            unmapped_fq: Box::new(std::io::sink()),
            index: &index,
            config: crate::align::AlignConfig::rna_default(),
            jlib: None,
            l1: Some(std::sync::Arc::new(l1)),
        };
        let (_stats, kept, _merged, _removed) = run_pe_2pass(&mut out, &r1, &r2, 2).unwrap();
        assert!(
            kept.iter().any(|(j, support, fate)| {
                *support >= 2 && *fate != DiscoveryFate::OutvotedByAnnotated && !j.minus_strand &&
                (j.start as i64 - 8000).abs() <= 2 && (j.end as i64 - 9000).abs() <= 2
            }),
            "novel junction [8000, 9000) must be discovered under L1; kept = {kept:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The pass-1 support vote at an annotated boundary: a shifted
    /// re-observation that out-votes the annotated copy (>2x pass-1 support)
    /// replaces it in the merged library; one that loses is reported but not
    /// merged; within a 2:1 band both stay. Novel junctions (>3 bp from any
    /// annotated one) always stay.
    #[test]
    fn filter_discoveries_support_vote() {
        let annotated = JunctionLib::build([Junction {
            contig: 0,
            start: 100,
            end: 200,
            minus_strand: false,
        }]);
        let j = |start, end| Junction {
            contig: 0,
            start,
            end,
            minus_strand: false,
        };
        // annotated (100,200) observed 40 times in pass 1 (exact re-observation)
        let mut discovered: Discoveries = vec![(j(100, 200), 40, 40 * 60)];
        // shifted +1 with support 7: loses the vote (40 > 2*7)
        discovered.push((j(101, 200), 7, 7 * 60));
        // novel junction: kept
        discovered.push((j(100, 900), 3, 3 * 60));
        // support 1: quality-gated out entirely
        discovered.push((j(5000, 6000), 1, 60));
        let (got, removed) = filter_discoveries(discovered, Some(&annotated));
        assert_eq!(got.len(), 3, "support-1 junction must be gated out");
        assert!(got
            .iter()
            .any(|(j, _, f)| j.start == 100 && j.end == 200 && *f == DiscoveryFate::TieKept));
        assert!(got
            .iter()
            .any(|(j, _, f)| j.start == 101 && *f == DiscoveryFate::OutvotedByAnnotated));
        assert!(got
            .iter()
            .any(|(j, _, f)| j.end == 900 && *f == DiscoveryFate::Kept));
        assert!(removed.is_empty());

        // the shifted copy out-votes an unobserved annotated boundary
        let discovered2: Discoveries = vec![(j(101, 200), 40, 40 * 60)];
        let (got2, removed2) = filter_discoveries(discovered2, Some(&annotated));
        assert!(got2
            .iter()
            .any(|(j, _, f)| j.start == 101 && *f == DiscoveryFate::OutvotesAnnotated));
        assert_eq!(removed2.len(), 1);
        assert_eq!(removed2[0], j(100, 200));

        // without an annotation every survivor is novel-kept
        let (got3, removed3) = filter_discoveries(vec![(j(101, 200), 7, 7 * 60)], None);
        assert!(got3.iter().all(|(_, _, f)| *f == DiscoveryFate::Kept));
        assert!(removed3.is_empty());
    }

    /// Shared unmapped.fq sink: keeps the bytes run_pe wrote so the test can
    /// inspect them after the pipeline finished.
    #[derive(Clone)]
    struct SharedSink(Arc<std::sync::Mutex<Vec<u8>>>);

    impl Write for SharedSink {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            self.0.lock().expect("sink lock").extend_from_slice(buf);
            Ok(buf.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    /// dedup_pe contract: unique table in first-seen order, per-pair uid maps
    /// in chunk order, one uid per sequence across pairs AND across mates.
    #[test]
    fn dedup_pe_uid_mapping_stable_and_order_preserved() {
        let rec = |name: &str, seq: &str| FastqRecord {
            name: name.as_bytes().to_vec(),
            seq: seq.as_bytes().to_vec(),
            qual: b"I".repeat(seq.len()),
        };
        // pair0/pair2 share the R1 seq; pair0's R1 seq also appears as
        // pair1's R2 seq — the unified table must give it ONE uid.
        let chunk = vec![
            (rec("a/1", "AAAA"), rec("a/2", "CCCC")),
            (rec("b/1", "GGGG"), rec("b/2", "AAAA")),
            (rec("c/1", "AAAA"), rec("c/2", "TTTT")),
        ];
        let (unique_seqs, uid1, uid2) = dedup_pe(&chunk);
        // unique table: first-seen order, each seq exactly once
        assert_eq!(
            unique_seqs,
            vec![&b"AAAA"[..], &b"CCCC"[..], &b"GGGG"[..], &b"TTTT"[..]]
        );
        // one uid per pair position, in chunk order
        assert_eq!(uid1, vec![0, 2, 0]);
        assert_eq!(uid2, vec![1, 0, 3]);
        // same seq → same uid, across pairs (R1 dup) and across mates
        assert_eq!(uid1[0], uid1[2]);
        assert_eq!(uid1[0], uid2[1]);
        // empty chunk: empty tables
        let (u, m1, m2) = dedup_pe(&[]);
        assert!(u.is_empty());
        assert!(m1.is_empty());
        assert!(m2.is_empty());
    }

    /// Paired-end deduplication must not change per-pair outputs: a chunk whose
    /// sequences are duplicated (a second block of pairs with identical
    /// R1/R2 seqs but different names and quals) must yield exactly the same
    /// per-pair results as the unduplicated chunk — every stats counter
    /// doubles, the discovery stream is the exact block concatenation, and
    /// the unmapped.fq seq multiset doubles (quality must not move the
    /// primary alignment).
    #[test]
    fn pe_dedup_duplicated_chunk_matches_unduplicated_outputs() {
        // deterministic pseudo-random reference, no repeats (same recipe as
        // the L1 discovery test)
        let n = 20_000usize;
        let mut seq = Vec::with_capacity(n);
        let mut state: u64 = 0x9E3779B97F4A7C15;
        for _ in 0..n {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            seq.push(b"ACGT"[(state >> 33) as usize & 3]);
        }
        // GT..AG intron [8000, 9000) so one pair exercises the junction path
        seq[8000] = b'G';
        seq[8001] = b'T';
        seq[8998] = b'A';
        seq[8999] = b'G';
        // an R2 with no k-mer in the reference (LCG continuation — fixed
        // bytes, so non-matching is a fixed fact, not a flake)
        let mut absent = Vec::with_capacity(150);
        for _ in 0..150 {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            absent.push(b"ACGT"[(state >> 33) as usize & 3]);
        }
        // pairs: (R1, R2). p1's R1 duplicates p0's R1; p2's short R1 makes
        // its own placement poor (rescue attempt consuming the cached
        // placements); p4's absent R2 gives the unplaced-mate rescue path
        // (cache-free).
        let spliced = {
            let mut s = seq[7925..8000].to_vec();
            s.extend_from_slice(&seq[9000..9075]);
            s
        };
        let pairs: Vec<(Vec<u8>, Vec<u8>)> = vec![
            (seq[1000..1150].to_vec(), seq[1200..1350].to_vec()),
            (seq[1000..1150].to_vec(), seq[3000..3150].to_vec()),
            (seq[5000..5080].to_vec(), seq[4000..4150].to_vec()),
            (spliced, seq[9100..9250].to_vec()),
            (seq[7000..7150].to_vec(), absent),
        ];
        let dir = std::env::temp_dir().join(format!("esp-pe-dedup-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let rec_bytes = |name: &str, s: &[u8], qual: u8| -> Vec<u8> {
            [
                &format!("@{name}\n").into_bytes()[..],
                s,
                b"\n+\n",
                &vec![qual; s.len()],
                b"\n",
            ]
            .concat()
        };
        let mut base1 = Vec::new();
        let mut base2 = Vec::new();
        let mut dup1 = Vec::new();
        let mut dup2 = Vec::new();
        for (i, (s1, s2)) in pairs.iter().enumerate() {
            base1.extend(rec_bytes(&format!("p{i}"), s1, b'I'));
            base2.extend(rec_bytes(&format!("p{i}"), s2, b'I'));
            dup1.extend(rec_bytes(&format!("p{i}"), s1, b'I'));
            dup2.extend(rec_bytes(&format!("p{i}"), s2, b'I'));
        }
        // duplicated block: same seqs, different names and quals
        for (i, (s1, s2)) in pairs.iter().enumerate() {
            dup1.extend(rec_bytes(&format!("q{i}"), s1, b'J'));
            dup2.extend(rec_bytes(&format!("q{i}"), s2, b'J'));
        }
        let base_r1 = dir.join("base_1.fq");
        let base_r2 = dir.join("base_2.fq");
        let dup_r1 = dir.join("dup_1.fq");
        let dup_r2 = dir.join("dup_2.fq");
        std::fs::write(&base_r1, &base1).unwrap();
        std::fs::write(&base_r2, &base2).unwrap();
        std::fs::write(&dup_r1, &dup1).unwrap();
        std::fs::write(&dup_r2, &dup2).unwrap();

        let reference = crate::fasta::parse_fasta_bytes(&{
            let mut t = b">chrT\n".to_vec();
            for chunk in seq.chunks(60) {
                t.extend_from_slice(chunk);
                t.push(b'\n');
            }
            t
        })
        .unwrap();
        let index = crate::index::Index::build(reference, crate::seed::SeedParams { k: 15, w: 5 });

        let run = |r1: &Path, r2: &Path| -> (AlignStats, Discoveries, Vec<u8>) {
            let cell = Arc::new(std::sync::Mutex::new(Vec::new()));
            let mut out = PipelineOut {
                bam: None,
                unmapped_fq: Box::new(SharedSink(cell.clone())),
                index: &index,
                config: crate::align::AlignConfig::rna_default(),
                jlib: None,
                l1: None,
            };
            let (stats, disc) = run_pe(&mut out, r1, r2, 2).unwrap();
            let raw = std::mem::take(&mut *cell.lock().unwrap());
            (stats, disc, raw)
        };
        let (s1, d1, raw1) = run(&base_r1, &base_r2);
        let (s2, d2, raw2) = run(&dup_r1, &dup_r2);

        // the base run is non-trivial: something maps, the spliced pair
        // discovers the junction, the absent-R2 pair leaves unmapped output
        assert!(s1.mapped_reads > 0, "base run must map something");
        assert!(!d1.is_empty(), "spliced pair must discover the junction");
        let unm_seq_lines = |raw: &[u8]| -> Vec<Vec<u8>> {
            use std::io::Read as _;
            let mut txt = Vec::new();
            flate2::read::GzDecoder::new(raw)
                .read_to_end(&mut txt)
                .unwrap();
            txt.split(|&b| b == b'\n')
                .enumerate()
                .filter(|(i, _)| i % 4 == 1)
                .map(|(_, l)| l.to_vec())
                .collect()
        };
        let mut u1 = unm_seq_lines(&raw1);
        assert!(!u1.is_empty(), "absent R2 must stay unmapped");

        // every observable output doubles exactly
        assert_eq!(s2.total_reads, 2 * s1.total_reads);
        assert_eq!(s2.mapped_reads, 2 * s1.mapped_reads);
        assert_eq!(s2.unmapped_reads, 2 * s1.unmapped_reads);
        assert_eq!(s2.proper_pairs, s1.proper_pairs.map(|x| 2 * x));
        assert_eq!(s2.junctions_total, 2 * s1.junctions_total);
        assert_eq!(s2.rescued_total, 2 * s1.rescued_total);
        assert_eq!(s2.rescue_fail_total, 2 * s1.rescue_fail_total);
        assert_eq!(s2.decoy_reads, 2 * s1.decoy_reads);
        for (a, b) in s1.mapq_hist.iter().zip(s2.mapq_hist.iter()) {
            assert_eq!(2 * a, *b, "mapq histogram must double");
        }
        // discovery stream: exact block duplication (pair order preserved)
        assert_eq!(d2.len(), 2 * d1.len());
        assert_eq!(&d2[..d1.len()], &d1[..]);
        assert_eq!(&d2[d1.len()..], &d1[..]);
        // unmapped.fq: same seq multiset, doubled
        let mut u2 = unm_seq_lines(&raw2);
        u1.sort();
        u2.sort();
        let mut doubled = u1.clone();
        doubled.extend_from_slice(&u1);
        assert_eq!(u2, doubled);

        let _ = std::fs::remove_dir_all(&dir);
    }

}
