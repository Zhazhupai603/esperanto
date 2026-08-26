//! Alignment pipeline: chunked streaming FASTQ → rayon align → ordered write of
//! BAM / unmapped.fq.gz / align_qc.json / align.baln.
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

use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;

use crate::align::{AlignConfig, Aligner};
use crate::bam::{self, flag, BamRecord};
use crate::error::AlignError;
use crate::extend::CigarOp;
use crate::fastq::{FastqReader, FastqRecord, RecordSource};
use crate::gtf::{Junction, JunctionLib};
use crate::index::Index;
use crate::jkmer::{self, JkmerIndex};
use crate::mapq::{mapq_of, ReadAlignment};
use crate::pair::{relate, template_length, PairRelation};
use crate::seed::Strand;
use crate::stats::{AlignStats, StatsAcc};

/// Streaming block size (pairs / reads).
pub const CHUNK: usize = 100_000;

/// Track 2 routing: total soft-clip above this triggers the jkmer path.
pub const SOFTCLIP_THRESHOLD: u32 = 20;
/// Track 2 acceptance: local_confirm total below this rejects the candidate.
pub const TRACK2_LOCAL_THRESHOLD: i32 = 30;

/// Pipeline output sinks (bam sink ownership is taken at run time —
/// multithreaded BGZF requires 'static).
pub struct PipelineOut<'a> {
    pub bam: Option<Box<dyn Write + Send>>,
    pub unmapped_fq: Box<dyn Write + Send>,
    pub index: &'a Index,
    pub config: AlignConfig,
    /// Junction library (sjdb; RNA mode).
    pub jlib: Option<Arc<JunctionLib>>,
    /// Junction k-mer index (Track 2; None = DNA mode or no --jkmer).
    pub jkmer: Option<Arc<JkmerIndex>>,
    /// L1 transcriptome-first engine bundle (None = pure legacy path).
    pub l1: Option<Arc<esperanto_engine::L1Index>>,
    /// .baln fast internal channel (for call; None = BAM only).
    pub baln: Option<Box<dyn Write + Send>>,
}

/// pass1 junction discoveries: (junction, support, summed MAPQ of supporters).
pub type Discoveries = Vec<(Junction, u32, u64)>;

/// 2-pass outputs: (stats, pass1 kept discoveries, final library list).
pub type TwoPassOut = (AlignStats, Vec<(Junction, u32)>, Vec<(Junction, u32)>);

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

/// Single-end pipeline (returns stats + pass1 discoveries, empty in DNA mode).
pub fn run_se(
    out: &mut PipelineOut,
    r1: &Path,
    threads: usize,
) -> Result<(AlignStats, Discoveries), AlignError> {
    let t0 = std::time::Instant::now();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| AlignError::FastqFormat {
            line: 0,
            msg: format!("thread pool: {e}"),
        })?;
    let header = bam::build_header(out.index);
    // --no-bam (profiling mode): bam=None → skip per-record BAM encoding.
    let mut bam_writer = match out.bam.take() {
        Some(sink) => Some(bam::create_writer(sink, &header, threads).map_err(io_err)?),
        None => None,
    };
    // .baln fast channel (for call) — header first.
    let mut baln_writer = match out.baln.take() {
        Some(w) => {
            let mut bw = std::io::BufWriter::with_capacity(8 << 20, w);
            let names: Vec<String> = out
                .index
                .reference
                .contigs
                .iter()
                .map(|c| c.name.clone())
                .collect();
            crate::baln::write_header(&mut bw, &names).map_err(io_err)?;
            std::io::Write::flush(&mut bw).map_err(io_err)?;
            Some(bw)
        }
        None => None,
    };
    let mut unm = GzEncoder::new(&mut out.unmapped_fq, Compression::new(1));
    let mut stats = StatsAcc::new();
    let mut discoveries: Discoveries = Vec::new();
    let mut reader: Box<dyn RecordSource> = if r1.extension().is_some_and(|e| e == "bfq") {
        Box::new(crate::fastq::BfqReader::open(r1)?)
    } else {
        Box::new(FastqReader::open(r1)?)
    };

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
        // READCACHE: align_read(seq) depends only on seq, not qual/name.
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
        // Track 2 dispatch per unique seq (deterministic: same seq + same
        // jkmer + same reference → same mutation).
        let mut cached_alns = cached_alns;
        if let Some(jk) = &out.jkmer {
            let reference = out.index.reference;
            for (uid, seq) in unique_seqs.iter().enumerate() {
                try_track2(&mut cached_alns[uid], seq, Some(jk), reference);
            }
        }
        // Build BamRecord per read in input order, cloning the cached
        // alignment and pairing it with the read's own name/qual.
        let mut results: Vec<BamRecord> = Vec::with_capacity(chunk.len());
        for (i, rec) in chunk.iter().enumerate() {
            let aln = cached_alns[uid_of_read[i]].clone();
            let name = String::from_utf8_lossy(&rec.name);
            results.push(bam::record_se(&name, &rec.seq, &rec.qual, aln));
        }
        for (i, (rec, br)) in chunk.iter().zip(&results).enumerate() {
            let aln = &cached_alns[uid_of_read[i]];
            stats.push_read(aln.is_some(), br.mapq);
            if let Some(x) = aln {
                if x.rescued {
                    stats.push_rescued();
                }
                stats.push_junctions(x.junctions.len() as u64);
                for j in &x.junctions {
                    discoveries.push((j.junction, 1, u64::from(br.mapq)));
                }
            }
            if aln.is_none() {
                write_unmapped(&mut unm, rec)?;
            }
            if let Some(w) = &mut bam_writer {
                bam::write_record(w, &header, br).map_err(io_err)?;
            }
            if let Some(bw) = &mut baln_writer {
                crate::baln::write_record(&mut *bw, br).map_err(io_err)?;
            }
        }
    }
    if let Some(w) = bam_writer {
        w.into_inner().finish().map_err(io_err)?;
    }
    if let Some(mut bw) = baln_writer {
        std::io::Write::flush(&mut bw).map_err(io_err)?;
    }
    unm.finish().map_err(io_err2)?;
    let mode = if out.config.rna { "rna-se" } else { "dna-se" };
    Ok((
        stats.finalize(mode, false, t0.elapsed().as_secs_f64()),
        discoveries,
    ))
}

/// READCACHE SE dedup: group reads in a chunk by identical seq → unique seq
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
        .build()
        .map_err(|e| AlignError::FastqFormat {
            line: 0,
            msg: format!("thread pool: {e}"),
        })?;
    let header = bam::build_header(out.index);
    let mut bam_writer = match out.bam.take() {
        Some(sink) => Some(bam::create_writer(sink, &header, threads).map_err(io_err)?),
        None => None,
    };
    // Legacy parity: run_pe never writes .baln — the driver still creates
    // the sink file (0 bytes), matching the old pipeline.
    drop(out.baln.take());
    let mut unm = GzEncoder::new(&mut out.unmapped_fq, Compression::new(1));
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
        // Parallel section aligns only (zero string/seq clones); BAM records
        // are built in the serial write section.
        let mut results: Vec<(Option<ReadAlignment>, Option<ReadAlignment>)> = pool.install(|| {
            let lib = lib.clone();
            chunk
                .par_iter()
                .map_init(
                    move || make_aligner(index, cfg, &lib, &l1_clone),
                    |aligner, (a, b)| {
                        let mut r1 = aligner.align_read(&a.seq);
                        let mut r2 = aligner.align_read(&b.seq);
                        // Rescue channel: one end unmapped + mate placed →
                        // anchored-window A/G-masking relocation.
                        let mate_for_r2 = r1.as_ref().map(|m| (m.contig, m.pos));
                        let mate_for_r1 = r2.as_ref().map(|m| (m.contig, m.pos));
                        if r1.is_none() {
                            if let Some((mc, mp)) = mate_for_r1 {
                                let stub = ReadAlignment {
                                    contig: mc,
                                    pos: mp,
                                    ..Default::default()
                                };
                                r1 = aligner.rescue_with_mate_anchor(&a.seq, &stub, 500);
                            }
                        }
                        if r2.is_none() {
                            if let Some((mc, mp)) = mate_for_r2 {
                                let stub = ReadAlignment {
                                    contig: mc,
                                    pos: mp,
                                    ..Default::default()
                                };
                                r2 = aligner.rescue_with_mate_anchor(&b.seq, &stub, 500);
                            }
                        }
                        (r1, r2)
                    },
                )
                .collect()
        });
        // Track 2 dispatch for PE (serial, post-align; same pattern as run_se).
        if let Some(jk) = &out.jkmer {
            let reference = out.index.reference;
            for ((a, b), (aln1, aln2)) in chunk.iter().zip(results.iter_mut()) {
                try_track2(aln1, &a.seq, Some(jk), reference);
                try_track2(aln2, &b.seq, Some(jk), reference);
            }
        }
        for ((a, b), (aln1, aln2)) in chunk.iter().zip(&results) {
            if let Some(x) = aln1 {
                stats.push_junctions(x.junctions.len() as u64);
                for j in &x.junctions {
                    discoveries.push((j.junction, 1, u64::from(mapq_of(x))));
                }
            }
            if let Some(x) = aln2 {
                stats.push_junctions(x.junctions.len() as u64);
                for j in &x.junctions {
                    discoveries.push((j.junction, 1, u64::from(mapq_of(x))));
                }
            }
            let (ra, rb) = pair_records(
                &String::from_utf8_lossy(&a.name),
                &a.seq,
                &a.qual,
                aln1.clone(),
                &String::from_utf8_lossy(&b.name),
                &b.seq,
                &b.qual,
                aln2.clone(),
            );
            stats.push_read(aln1.is_some(), ra.mapq);
            stats.push_read(aln2.is_some(), rb.mapq);
            if let Some(x) = aln1 {
                if x.rescued {
                    stats.push_rescued();
                }
            }
            if let Some(x) = aln2 {
                if x.rescued {
                    stats.push_rescued();
                }
            }
            if aln1.is_none() {
                if aln2.is_some() {
                    stats.push_rescue_fail();
                }
                write_unmapped_mate(&mut unm, a, 1)?;
            }
            if aln2.is_none() {
                if aln1.is_some() {
                    stats.push_rescue_fail();
                }
                write_unmapped_mate(&mut unm, b, 2)?;
            }
            if let (Some(r1a), Some(r2a)) = (aln1, aln2) {
                if relate(r1a, r2a) == PairRelation::Proper {
                    stats.push_proper_pair();
                    stats.push_insert(template_length(r1a, r2a));
                }
            }
            if let Some(w) = &mut bam_writer {
                bam::write_record(w, &header, &ra).map_err(io_err)?;
                bam::write_record(w, &header, &rb).map_err(io_err)?;
            }
        }
    }
    if let Some(w) = bam_writer {
        w.into_inner().finish().map_err(io_err)?;
    }
    unm.finish().map_err(io_err2)?;
    let mode = if out.config.rna { "rna-pe" } else { "dna-pe" };
    Ok((
        stats.finalize(mode, true, t0.elapsed().as_secs_f64()),
        discoveries,
    ))
}

/// Build PE paired records (FLAG/mate fields per spec).
#[allow(clippy::too_many_arguments)]
fn pair_records(
    n1: &str,
    s1: &[u8],
    q1: &[u8],
    a1: Option<ReadAlignment>,
    n2: &str,
    s2: &[u8],
    q2: &[u8],
    a2: Option<ReadAlignment>,
) -> (BamRecord, BamRecord) {
    let proper = match (&a1, &a2) {
        (Some(x), Some(y)) => relate(x, y) == PairRelation::Proper,
        _ => false,
    };
    let tlen = match (&a1, &a2) {
        (Some(x), Some(y)) if proper => template_length(x, y),
        _ => 0,
    };
    let mk = |name: &str,
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
        let (s_out, q_out) = crate::bam::apply_t13(
            aln.as_ref().map(|a| a.strand) == Some(Strand::Minus),
            seq,
            qual,
        );
        BamRecord {
            name: name.to_string(),
            flag: f,
            aln: aln.as_ref().map(crate::bam::aln_view),
            mapq: q,
            seq: s_out,
            qual: q_out,
            mate: Some((mate_c, mate_p, tlen)),
        }
    };
    (
        mk(n1, s1, q1, &a1, &a2, flag::READ1, tlen),
        mk(n2, s2, q2, &a2, &a1, flag::READ2, -tlen),
    )
}

/// unmapped.fq write (SE).
fn write_unmapped<W: Write>(w: &mut W, rec: &FastqRecord) -> Result<(), AlignError> {
    w.write_all(b"@")
        .and_then(|_| w.write_all(&rec.name))
        .and_then(|_| w.write_all(b"\n"))
        .and_then(|_| w.write_all(&rec.seq))
        .and_then(|_| w.write_all(b"\n+\n"))
        .and_then(|_| w.write_all(&rec.qual))
        .and_then(|_| w.write_all(b"\n"))
        .map_err(io_err2)
}

/// unmapped.fq write (PE; name carries the /1 /2 mate marker — rescue contract).
fn write_unmapped_mate<W: Write>(
    w: &mut W,
    rec: &FastqRecord,
    which: u8,
) -> Result<(), AlignError> {
    w.write_all(b"@")
        .and_then(|_| w.write_all(&rec.name))
        .and_then(|_| write!(w, "/{which}"))
        .and_then(|_| w.write_all(&rec.seq))
        .and_then(|_| w.write_all(b"\n+\n"))
        .and_then(|_| w.write_all(&rec.qual))
        .and_then(|_| w.write_all(b"\n"))
        .map_err(io_err2)
}

fn io_err(e: std::io::Error) -> AlignError {
    AlignError::FastqIo {
        path: "<bam>".into(),
        source: e,
    }
}

fn io_err2(e: std::io::Error) -> AlignError {
    AlignError::FastqIo {
        path: "<unmapped.fq.gz>".into(),
        source: e,
    }
}

/// SE 2-pass (same semantics as the PE version).
pub fn run_se_2pass(
    out: &mut PipelineOut,
    r1: &Path,
    threads: usize,
) -> Result<TwoPassOut, AlignError> {
    // With L1 (transcriptome-first), annotated junctions are already
    // available and the genome fallback detects novel introns from the seed
    // chain — pass1 discovery adds nothing. Single-pass keeps the same map
    // rate at half the cost.
    if false { // TEMP: force 2-pass to test novel recovery
        let (stats, _) = run_se(out, r1, threads)?;
        return Ok((stats, Vec::new(), Vec::new()));
    }
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
    let kept = filter_discoveries(discovered);
    merged.extend(kept.iter().copied());
    let final_lib = JunctionLib::build_with_counts(merged);
    let final_list: Vec<(Junction, u32)> = final_lib
        .junctions
        .iter()
        .copied()
        .zip(final_lib.counts.iter().copied())
        .collect();
    out.jlib = Some(Arc::new(final_lib));
    let (stats, _) = run_se(out, r1, threads)?;
    Ok((stats, kept, final_list))
}

/// RNA 2-pass: pass1 discovers junctions (≥2 support into the library) →
/// merge with sjdb → pass2 aligns with the merged library.
pub fn run_pe_2pass(
    out: &mut PipelineOut,
    r1: &Path,
    r2: &Path,
    threads: usize,
) -> Result<TwoPassOut, AlignError> {
    // With L1 (transcriptome-first), annotated junctions are already
    // available and the genome fallback detects novel introns from the seed
    // chain — pass1 discovery adds nothing. Single-pass keeps the same map
    // rate at half the cost.
    if false { // TEMP: force 2-pass to test novel recovery
        let (stats, _) = run_pe(out, r1, r2, threads)?;
        return Ok((stats, Vec::new(), Vec::new()));
    }
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
    let kept = filter_discoveries(discovered);
    merged.extend(kept.iter().copied());
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
    Ok((stats, kept, final_list))
}

/// Junction-level filter for 2-pass: support ≥ 2 AND mean supporter MAPQ ≥ 20
/// AND span ≤ 500kb (low-quality cross-repeat-family junctions stay out).
pub fn filter_discoveries(discovered: Discoveries) -> Vec<(Junction, u32)> {
    let mut acc: std::collections::BTreeMap<Junction, (u32, u64)> =
        std::collections::BTreeMap::new();
    for (j, c, mq) in discovered {
        let e = acc.entry(j).or_default();
        e.0 += c;
        e.1 += mq;
    }
    acc.into_iter()
        .filter(|(j, (c, mqsum))| {
            *c >= 2 && *mqsum >= u64::from(*c) * 20 && j.end.saturating_sub(j.start) <= 500_000
        })
        .map(|(j, (c, _))| (j, c))
        .collect()
}

// =========================================================================
// Track 2 routing hook (dual-track)
// =========================================================================

/// Total soft-clip of an alignment's CIGAR.
fn total_softclip(aln: &ReadAlignment) -> u32 {
    aln.cigar
        .iter()
        .map(|op| match op {
            CigarOp::SoftClip(n) => *n,
            _ => 0,
        })
        .sum()
}

/// Route one read through Track 2 (junction k-mer direct placement).
///
/// Returns true when `*aln` was replaced (Track 2 hit passing local_confirm).
/// False = keep the Track 1 result.
///
/// Note: `*aln` is in-out; a Track 2 hit replaces it with a synthetic
/// ReadAlignment (rescued=true, junctions empty — Track 2 must not pollute
/// Track 1 junction statistics).
///
/// Quirk (parity-frozen): the Track 2 MAPQ model (`compute_track2_mapq`) is
/// NOT used for the written record; the written MAPQ comes from the standard
/// mapq formula applied to the synthetic alignment (chain_score=total,
/// second=0, n_anchors=hits ⇒ round(60 × min(1, hits/10))).
pub fn try_track2(
    aln: &mut Option<ReadAlignment>,
    read_seq: &[u8],
    jkmer: Option<&Arc<JkmerIndex>>,
    reference: &crate::fasta::Reference,
) -> bool {
    let Some(jk) = jkmer else {
        return false; // DNA mode or no --jkmer
    };
    if read_seq.len() < jkmer::K {
        return false;
    }
    // Routing condition: Track 1 unmapped or soft-clip above threshold.
    let should_try = match aln {
        None => true,
        Some(a) => total_softclip(a) > SOFTCLIP_THRESHOLD,
    };
    if !should_try {
        return false;
    }
    let candidates = jk.query_read(read_seq);
    if candidates.is_empty() {
        return false;
    }
    // Evaluate each candidate: local_confirm must pass the threshold; the
    // highest total wins (strictly greater replaces — deterministic).
    let mut best: Option<(usize, usize, i32, usize, jkmer::Breakpoint)> = None;
    for (cand_idx, cand) in candidates.iter().enumerate() {
        let junc = cand.junction;
        let Some(verdict) = cand.infer_breakpoint() else {
            continue;
        };
        let split = match verdict {
            jkmer::Breakpoint::HighConf(p) => p as usize,
            jkmer::Breakpoint::LowConf(p) => p as usize,
        };
        // Fetch flanks in read-forward orientation.
        let mut donor: Vec<u8> = Vec::new();
        let mut acceptor: Vec<u8> = Vec::new();
        let fetch = |contig: u32, start: u32, end: u32| -> Vec<u8> {
            match reference.contigs.get(contig as usize) {
                Some(c) => {
                    let s = start.min(c.len);
                    c.slice_ascii(s, end.min(c.len))
                }
                None => Vec::new(),
            }
        };
        jkmer::fetch_flanks(junc, fetch, &mut donor, &mut acceptor);
        let lc = jkmer::local_confirm(read_seq, split, &donor, &acceptor);
        if lc.total < TRACK2_LOCAL_THRESHOLD {
            continue;
        }
        let runner_up = candidates.get(1).map(|c| c.hits.len()).unwrap_or(0);
        let score_key = lc.total;
        if best.as_ref().is_none_or(|b| b.2 < score_key) {
            best = Some((cand_idx, cand.hits.len(), lc.total, runner_up, verdict));
        }
    }
    let Some((cand_idx, hit_count, total, _runner_up, verdict)) = best else {
        return false;
    };
    let cand = &candidates[cand_idx];
    let junc = cand.junction;
    let split = match verdict {
        jkmer::Breakpoint::HighConf(p) => p as usize,
        jkmer::Breakpoint::LowConf(p) => p as usize,
    };
    let left_match = split as u32;
    let right_match = (read_seq.len() as u32).saturating_sub(split as u32);
    let intron = junc.intron_len();
    // Track 1-compatible output: Match + RefSkip + Match, rescued=true.
    let cigar: Vec<CigarOp> = {
        let mut v = Vec::with_capacity(3);
        if left_match > 0 {
            v.push(CigarOp::Match(left_match));
        }
        v.push(CigarOp::RefSkip(intron));
        if right_match > 0 {
            v.push(CigarOp::Match(right_match));
        }
        v
    };
    let pos_0based = match junc.strand {
        // plus: read[0] at intron_start − left_match (0-based)
        Strand::Plus => junc.intron_start.saturating_sub(left_match),
        // minus: left_match consumes the high-coordinate exon; the genomic
        // left end of the transcript-forward alignment is intron_start − right_match
        Strand::Minus => junc.intron_start.saturating_sub(right_match),
    };
    let strand_seed = junc.strand;
    let new_aln = ReadAlignment {
        contig: junc.contig,
        pos: pos_0based,
        strand: strand_seed,
        score: total,
        chain_score: total,
        second_chain_score: 0,
        cigar,
        n_anchors: hit_count,
        junctions: Vec::new(), // Track 2 must not pollute Track 1 junction stats
        ea_count: 0,
        mm_count: 0,
        n_seeds: hit_count,
        rescued: true,
    };
    *aln = Some(new_aln);
    true
}
