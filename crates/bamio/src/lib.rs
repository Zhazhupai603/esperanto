//! esperanto-bamio — BAM record-level read/write (reference-forward orientation contract).
//!
//! Spec: docs/specs/crates/bamio.md (single source of truth).
//! Write side: record model + binary encoder + multithreaded BGZF writer + header construction;
//! Read side: sequential record view (read module) + original-orientation restore helpers.

pub mod baln;
pub mod encode;
pub mod read;
pub mod sort;

use std::cell::Cell;
use std::io::{self, Write as IoWrite};

use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::cigar::op::{Kind, Op};
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::header::record::value::map::Map;
use noodles::sam::header::record::value::map::ReferenceSequence;
use noodles::sam::Header;

pub use read::{restore_original, InRecord};

/// SAM FLAG bits.
pub mod flag {
    pub const PAIRED: u16 = 0x1;
    pub const PROPER_PAIR: u16 = 0x2;
    pub const UNMAPPED: u16 = 0x4;
    pub const MATE_UNMAPPED: u16 = 0x8;
    pub const REVERSE: u16 = 0x10;
    pub const MATE_REVERSE: u16 = 0x20;
    pub const READ1: u16 = 0x40;
    pub const READ2: u16 = 0x80;
    pub const SECONDARY: u16 = 0x100;
    pub const SUPPLEMENTARY: u16 = 0x800;
}

/// CIGAR operation (SAM; `=`/`X` folded into `Match`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CigarOp {
    /// Aligned bases (match or mismatch).
    Match(u32),
    /// Insertion to the reference (consumes read only).
    Ins(u32),
    /// Deletion from the reference (consumes ref only).
    Del(u32),
    /// Reference skip (intron; consumes ref only).
    RefSkip(u32),
    /// Soft-clipped read bases (consume read only).
    SoftClip(u32),
}

/// Aux tag value (the only three types bamio produces).
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TagValue {
    /// A: printable character.
    Char(u8),
    /// i: 32-bit integer.
    Int(i32),
    /// Z/H: NUL-terminated string (stored without the NUL).
    Str(String),
}

/// Ordered aux tag: encoded verbatim in Vec order.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RawTag(pub [u8; 2], pub TagValue);

/// Reference-frame view of an aligned record.
#[derive(Clone, Debug)]
pub struct AlnView {
    /// Contig index in header SQ order.
    pub contig: u32,
    /// 0-based leftmost reference position.
    pub pos: u32,
    pub cigar: Vec<CigarOp>,
    /// Ordered tags (map artifacts use fixed order XS → EA → EK → RE, built by the caller).
    pub tags: Vec<RawTag>,
}

/// A record to write out (single SE record or one end of a PE pair).
#[derive(Clone, Debug)]
pub struct OutRecord {
    pub name: String,
    pub flag: u16,
    /// ≤60; unmapped records encode as 255 (BAM convention).
    pub mapq: u8,
    /// None = unmapped.
    pub aln: Option<AlnView>,
    /// Stored orientation (reference-forward already applied).
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    /// (mate contig id, mate pos0, tlen); contig < 0 means absent.
    pub mate: Option<(i32, i32, i32)>,
}

/// revcomp: shared by the reference-forward write side and restore. Case-insensitive decode, A↔T, C↔G;
/// all other characters (including IUPAC degenerate codes) emit `N`; output is all uppercase.
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

/// Reference-forward write side: minus-strand aligned records store reference-forward SEQ + reversed QUAL; plus-strand/unmapped unchanged.
pub fn apply_t13(is_reverse: bool, seq: &[u8], qual: &[u8]) -> (Vec<u8>, Vec<u8>) {
    if is_reverse {
        (revcomp(seq), qual.iter().rev().copied().collect())
    } else {
        (seq.to_vec(), qual.to_vec())
    }
}

/// Build a BAM header from a contig table. Comment supplied by the caller (preserves per-artifact historical strings).
pub fn build_header(contigs: &[(String, u64)], comment: &str) -> Header {
    let mut builder = Header::builder();
    for (name, len) in contigs {
        let rs = Map::<ReferenceSequence>::new(
            std::num::NonZeroUsize::new(*len as usize).expect("contig len > 0"),
        );
        builder = builder.add_reference_sequence(name.as_str(), rs);
    }
    builder.add_comment(comment).build()
}

/// Create a BAM writer (multithreaded BGZF compression — block-order merge preserves byte-level
/// determinism; workers = clamp(threads, 1, 4), scaling with alignment threads without overtaking them).
pub fn create_writer<W: IoWrite + Send + 'static>(
    w: W,
    header: &Header,
    threads: usize,
) -> io::Result<bam::io::Writer<noodles::bgzf::io::MultithreadedWriter<W>>> {
    let workers = std::num::NonZeroUsize::new(threads.clamp(1, 4)).expect("nonzero");
    let bgzf = noodles::bgzf::io::MultithreadedWriter::with_worker_count(workers, w);
    let mut writer = bam::io::Writer::from(bgzf);
    writer.write_header(header)?;
    Ok(writer)
}

/// Write a single record to a BAM writer.
///
/// Fast path: direct binary encode into a thread-local reusable buffer, bypassing
/// the noodles `RecordBuf` builder. Falls back to `write_record_slow` for record
/// shapes the fast encoder does not support.
pub fn write_record<W: IoWrite>(
    writer: &mut bam::io::Writer<W>,
    header: &Header,
    rec: &OutRecord,
) -> io::Result<()> {
    let mut buf = take_enc_buf();
    buf.clear();

    let result = match encode::try_encode(&mut buf, rec) {
        Some(Ok(())) => {
            let block_size = u32::try_from(buf.len())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
            let inner = writer.get_mut();
            inner.write_all(&block_size.to_le_bytes())?;
            inner.write_all(&buf)
        }
        Some(Err(e)) => Err(e),
        None => write_record_slow(writer, header, rec),
    };

    return_enc_buf(buf);
    result
}

thread_local! {
    static ENC_BUF: Cell<Vec<u8>> = Cell::new(Vec::with_capacity(8192));
}

fn take_enc_buf() -> Vec<u8> {
    ENC_BUF.with(Cell::take)
}

fn return_enc_buf(buf: Vec<u8>) {
    ENC_BUF.with(|cell| cell.set(buf));
}

/// Fallback: noodles RecordBuf path for unsupported record shapes.
#[doc(hidden)] // public only for differential testing of the fast/slow paths; normal calls go through write_record.
pub fn write_record_slow<W: IoWrite>(
    writer: &mut bam::io::Writer<W>,
    header: &Header,
    rec: &OutRecord,
) -> io::Result<()> {
    let mut rb = RecordBuf::builder()
        .set_name(rec.name.as_str())
        .set_flags(noodles::sam::alignment::record::Flags::from(rec.flag))
        .set_sequence(rec.seq.clone().into())
        .set_quality_scores(rec.qual.clone().into());
    if let Some(aln) = &rec.aln {
        rb = rb
            .set_reference_sequence_id(aln.contig as usize)
            .set_alignment_start(((aln.pos as usize + 1).try_into()).expect("pos>0"))
            .set_mapping_quality(
                noodles::sam::alignment::record::MappingQuality::try_from(rec.mapq)
                    .expect("mapq <= 60"),
            )
            .set_cigar(aln.cigar.iter().map(noodles_op).collect::<Vec<_>>().into_iter().collect());
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;
        let mut data = noodles::sam::alignment::record_buf::Data::default();
        for RawTag(tag, value) in &aln.tags {
            let t = Tag::new(tag[0], tag[1]);
            match value {
                TagValue::Char(c) => data.insert(t, Value::Character(*c)),
                TagValue::Int(i) => data.insert(t, Value::Int32(*i)),
                TagValue::Str(s) => data.insert(t, Value::String(s.clone().into())),
            };
        }
        rb = rb.set_data(data);
    }
    if let Some((mc, mpos, tlen)) = rec.mate {
        if mc >= 0 {
            rb = rb
                .set_mate_reference_sequence_id(mc as usize)
                .set_mate_alignment_start(((mpos as usize + 1).try_into()).expect("pos>0"))
                .set_template_length(tlen);
        }
    }
    let record = rb.build();
    AlignmentWrite::write_alignment_record(writer, header, &record)?;
    Ok(())
}

/// Convert CIGAR to a noodles Op sequence.
fn noodles_op(op: &CigarOp) -> Op {
    let (kind, n) = match op {
        CigarOp::Match(n) => (Kind::Match, *n),
        CigarOp::Ins(n) => (Kind::Insertion, *n),
        CigarOp::Del(n) => (Kind::Deletion, *n),
        CigarOp::SoftClip(n) => (Kind::SoftClip, *n),
        CigarOp::RefSkip(n) => (Kind::Skip, *n),
    };
    Op::new(kind, n as usize)
}
