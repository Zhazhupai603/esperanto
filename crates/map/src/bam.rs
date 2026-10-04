//! BAM write-out adapter layer: record construction (reference-forward/flag/mate/tag order) lives in
//! the map alignment domain; encoding, writer, and read side live in esperanto-bamio
//! (spec: docs/specs/crates/bamio.md).
//!
//! - Header SQ lines come from the index contig table;
//! - Aligned records: FLAG/RNAME/POS/MAPQ/CIGAR/SEQ/QUAL + mate fields (PE);
//! - Unmapped records: FLAG 0x4 (+ mate info if PE), written to BAM and to unmapped.fq.gz (rescue-channel hard contract);
//! - Determinism: record order = input order (sorting left to downstream); fixed tag order (NH → NM → XS → EA → EK → RE).

use esperanto_bamio::{AlnView, OutRecord, RawTag, TagValue};
use noodles::sam::Header;

use crate::evidence::evidence_tag;
use crate::extend::CigarOp;
use crate::index::Index;
use crate::mapq::{mapq_of, ReadAlignment};
use crate::seed::Strand;

pub use esperanto_bamio::{apply_t13, create_writer, flag, write_record};

/// Map-side record type = bamio output record.
pub type BamRecord = OutRecord;

/// Build the BAM header from the index (stable comment string).
pub fn build_header(index: &Index) -> Header {
    let contigs: Vec<(String, u64)> = index
        .reference
        .contigs
        .iter()
        .map(|c| (c.name.clone(), c.len as u64))
        .collect();
    esperanto_bamio::build_header(&contigs, "esperanto-map M2 (v0.1)")
}

/// CIGAR conversion (both = and X fold into Match).
pub(crate) fn convert_cigar(cigar: &[CigarOp]) -> Vec<esperanto_bamio::CigarOp> {
    cigar
        .iter()
        .map(|op| match op {
            CigarOp::Match(n) => esperanto_bamio::CigarOp::Match(*n),
            CigarOp::Ins(n) => esperanto_bamio::CigarOp::Ins(*n),
            CigarOp::Del(n) => esperanto_bamio::CigarOp::Del(*n),
            CigarOp::RefSkip(n) => esperanto_bamio::CigarOp::RefSkip(*n),
            CigarOp::SoftClip(n) => esperanto_bamio::CigarOp::SoftClip(*n),
        })
        .collect()
}

/// Edit distance to the reference (SAM `NM`): every base that differs plus
/// every inserted and deleted base. The editing-aware matches counted in
/// `ea_count` are real reference mismatches, so they belong in NM — an A-to-I
/// read is not a zero-distance alignment.
fn nm_of(a: &ReadAlignment) -> i32 {
    let indel: u32 = a
        .cigar
        .iter()
        .map(|op| match op {
            CigarOp::Ins(n) | CigarOp::Del(n) => *n,
            _ => 0,
        })
        .sum();
    (a.mm_count + a.ea_count() + indel) as i32
}

/// ReadAlignment → bamio reference-frame view; fixed tag order
/// NH → NM → XS → EA → EK → RE.
///
/// NH/NM lead because samtools, GATK and REDItools all expect them; the
/// esperanto-specific tags follow. NH is always 1: only the primary placement
/// is written (secondary and supplementary records are never emitted), so
/// there is no multi-hit count to report. MD is deliberately absent — every
/// consumer in this pipeline reads the reference FASTA directly, and a
/// spliced-CIGAR MD is easy to get subtly wrong.
pub(crate) fn aln_view(a: &ReadAlignment) -> AlnView {
    let mut tags = Vec::with_capacity(6);
    tags.push(RawTag(*b"NH", TagValue::Int(1)));
    tags.push(RawTag(*b"NM", TagValue::Int(nm_of(a))));
    if !a.junctions.is_empty() {
        let xs = if a.junctions[0].junction.minus_strand {
            b'-'
        } else {
            b'+'
        };
        tags.push(RawTag(*b"XS", TagValue::Char(xs)));
    }
    if a.ea_count() > 0 {
        tags.push(RawTag(*b"EA", TagValue::Int(a.ea_count() as i32)));
    }
    tags.push(RawTag(*b"EK", TagValue::Str(evidence_tag(a))));
    if a.rescued {
        tags.push(RawTag(*b"RE", TagValue::Str("unmapped".into())));
    }
    AlnView {
        contig: a.contig,
        pos: a.pos,
        cigar: convert_cigar(&a.cigar),
        tags,
    }
}

/// FASTQ quality lines are ASCII (Phred+33); the BAM QUAL field stores the
/// raw Phred value. Reads are carried through the aligner in the FASTQ
/// encoding, so the offset is removed once, at the record boundary.
pub fn phred_from_ascii(qual: &[u8]) -> Vec<u8> {
    qual.iter().map(|&q| q.saturating_sub(33)).collect()
}

/// Reference-forward SEQ + Phred-decoded QUAL in one pass — byte-identical to
/// `apply_t13(..)` followed by `phred_from_ascii(..)` (minus strand: revcomp
/// SEQ + reversed QUAL with the +33 offset removed in the same pass).
pub(crate) fn t13_phred(is_reverse: bool, seq: &[u8], qual: &[u8]) -> (Vec<u8>, Vec<u8>) {
    if is_reverse {
        (
            esperanto_bamio::revcomp(seq),
            qual.iter().rev().map(|&q| q.saturating_sub(33)).collect(),
        )
    } else {
        (
            seq.to_vec(),
            qual.iter().map(|&q| q.saturating_sub(33)).collect(),
        )
    }
}

/// Build a record (SE) from an alignment result.
///
/// Reference-forward SEQ direction (hard contract): minus-strand aligned records store SEQ = revcomp(original
/// read) with QUAL reversed in sync (STAR/BWA/htslib convention; required by direct samtools
/// pileup reads); unmapped keeps the original orientation.
///
/// `aln` is borrowed: the record copies out every field it needs (cigar/tags
/// via [`aln_view`], scores via `mapq_of`), so the caller's cached alignment
/// is never cloned or consumed.
pub fn record_se(name: &[u8], seq: &[u8], qual: &[u8], aln: Option<&ReadAlignment>) -> BamRecord {
    let rev = aln.is_some_and(|a| a.strand == Strand::Minus);
    let flag_bits = match aln {
        None => flag::UNMAPPED,
        Some(_) if rev => flag::REVERSE,
        Some(_) => 0,
    };
    let q = aln.map_or(0, mapq_of);
    let (seq, qual) = t13_phred(rev, seq, qual);
    BamRecord {
        name: name.to_vec(),
        flag: flag_bits,
        mapq: q,
        aln: aln.map(aln_view),
        seq,
        qual,
        mate: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::extend::CigarOp;
    use crate::gtf::RefinedJunction;

    fn synthetic_aln(strand: Strand, rescued: bool) -> ReadAlignment {
        ReadAlignment {
            contig: 0,
            pos: 100,
            strand,
            score: 140,
            second_score: 0,
            chain_score: 140,
            second_chain_score: 0,
            cigar: vec![
                CigarOp::SoftClip(2),
                CigarOp::Match(30),
                CigarOp::Ins(4),
                CigarOp::Del(2),
                CigarOp::RefSkip(60),
                CigarOp::Match(52),
            ],
            n_anchors: 12,
            junctions: vec![RefinedJunction {
                junction: crate::gtf::Junction {
                    contig: 0,
                    start: 130,
                    end: 190,
                    minus_strand: true,
                },
                signal: crate::gtf::SpliceSignal::GtAg,
                known_support: 0,
            }],
            ea_ag: 3,
            ea_tc: 0,
            mm_count: 2,
            n_seeds: 14,
            rescued,
            mapq_override: None,
            origin: crate::mapq::Origin::Genomic,
            n_cands: 1,
        }
    }

    const SEQ: &[u8] = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    const QUAL: &[u8] = b"IJIJIJIJIJIJIJIJIJIJIJIJIJIJIJIJIJIJIJIJIJIJI";

    /// t13_phred must equal the historical two-step path
    /// (apply_t13 then phred_from_ascii) byte for byte, both strands.
    #[test]
    fn t13_phred_matches_two_step_path() {
        for rev in [false, true] {
            let (fused_s, fused_q) = t13_phred(rev, SEQ, QUAL);
            let (s, q) = esperanto_bamio::apply_t13(rev, SEQ, QUAL);
            assert_eq!(fused_s, s, "seq mismatch rev={rev}");
            assert_eq!(fused_q, phred_from_ascii(&q), "qual mismatch rev={rev}");
        }
    }

    /// record_se field equality against the historical construction path:
    /// every field (flag/mapq/seq/qual/name/cigar/tags) must come out
    /// identical for plus, minus and unmapped records.
    #[test]
    fn record_se_fields_match_historical_path() {
        for (strand, rescued) in [
            (Strand::Plus, false),
            (Strand::Minus, false),
            (Strand::Minus, true),
        ] {
            let aln = synthetic_aln(strand, rescued);
            let rec = record_se(b"read\xC3\xA9", SEQ, QUAL, Some(&aln));
            let rev = strand == Strand::Minus;
            // old path: apply_t13 + phred_from_ascii
            let (old_s, old_q) = esperanto_bamio::apply_t13(rev, SEQ, QUAL);
            assert_eq!(rec.seq, old_s);
            assert_eq!(rec.qual, phred_from_ascii(&old_q));
            assert_eq!(rec.name, b"read\xC3\xA9");
            assert_eq!(rec.flag, if rev { flag::REVERSE } else { 0 });
            assert_eq!(rec.mapq, mapq_of(&aln));
            assert_eq!(rec.mate, None);
            let view = rec.aln.as_ref().expect("aligned record has AlnView");
            assert_eq!(view.contig, aln.contig);
            assert_eq!(view.pos, aln.pos);
            assert_eq!(view.cigar, convert_cigar(&aln.cigar));
            // tags in fixed order NH → NM → XS → EA → EK → RE
            let nm = nm_of(&aln);
            let mut expected: Vec<esperanto_bamio::RawTag> = vec![
                esperanto_bamio::RawTag(*b"NH", TagValue::Int(1)),
                esperanto_bamio::RawTag(*b"NM", TagValue::Int(nm)),
                esperanto_bamio::RawTag(*b"XS", TagValue::Char(b'-')),
                esperanto_bamio::RawTag(*b"EA", TagValue::Int(aln.ea_count() as i32)),
                esperanto_bamio::RawTag(*b"EK", TagValue::Str(evidence_tag(&aln))),
            ];
            if rescued {
                expected.push(esperanto_bamio::RawTag(
                    *b"RE",
                    TagValue::Str("unmapped".into()),
                ));
            }
            assert_eq!(view.tags, expected, "tags for {strand:?} rescued={rescued}");
        }
        // unmapped: original orientation, FLAG 0x4, MAPQ 0, no AlnView
        let rec = record_se(b"u", SEQ, QUAL, None);
        assert_eq!(rec.flag, flag::UNMAPPED);
        assert_eq!(rec.mapq, 0);
        assert_eq!(rec.seq, SEQ);
        assert_eq!(rec.qual, phred_from_ascii(QUAL));
        assert!(rec.aln.is_none());
    }
}
