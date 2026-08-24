//! BAM write-out adapter layer: record construction (reference-forward/flag/mate/tag order) lives in
//! the map alignment domain; encoding, writer, and read side live in esperanto-bamio
//! (spec: docs/specs/crates/bamio.md).
//!
//! - Header SQ lines come from the index contig table;
//! - Aligned records: FLAG/RNAME/POS/MAPQ/CIGAR/SEQ/QUAL + mate fields (PE);
//! - Unmapped records: FLAG 0x4 (+ mate info if PE), written to BAM and to unmapped.fq.gz (rescue-channel hard contract);
//! - Determinism: record order = input order (sorting left to downstream); fixed tag order (XS → EA → EK → RE).

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

/// Build the BAM header from the index (comment string frozen; any change shifts artifact bytes).
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

/// ReadAlignment → bamio reference-frame view; fixed tag order XS → EA → EK → RE.
pub(crate) fn aln_view(a: &ReadAlignment) -> AlnView {
    let mut tags = Vec::with_capacity(4);
    if !a.junctions.is_empty() {
        let xs = if a.junctions[0].junction.minus_strand {
            b'-'
        } else {
            b'+'
        };
        tags.push(RawTag(*b"XS", TagValue::Char(xs)));
    }
    if a.ea_count > 0 {
        tags.push(RawTag(*b"EA", TagValue::Int(a.ea_count as i32)));
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

/// Build a record (SE) from an alignment result.
///
/// Reference-forward SEQ direction (hard contract): minus-strand aligned records store SEQ = revcomp(original
/// read) with QUAL reversed in sync (STAR/BWA/htslib convention; required by direct samtools
/// pileup reads); unmapped keeps the original orientation.
pub fn record_se(name: &str, seq: &[u8], qual: &[u8], aln: Option<ReadAlignment>) -> BamRecord {
    let (flag_bits, q, out_seq, out_qual) = match &aln {
        Some(a) => {
            let rev = a.strand == Strand::Minus;
            let (s, q2) = apply_t13(rev, seq, qual);
            (if rev { flag::REVERSE } else { 0 }, mapq_of(a), s, q2)
        }
        None => (flag::UNMAPPED, 0, seq.to_vec(), qual.to_vec()),
    };
    BamRecord {
        name: name.to_string(),
        flag: flag_bits,
        mapq: q,
        aln: aln.as_ref().map(aln_view),
        seq: out_seq,
        qual: out_qual,
        mate: None,
    }
}
