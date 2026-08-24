//! [`L1Index`]: the adapter fusing a `.tidx` k-mer index, a full-projection
//! `TxMap`, and a transcript sequence store into one shared `tx_id` space
//! (resolves the runtime-GTF dual-track problem: production never parses
//! GTF).

use std::collections::BTreeMap;
use std::io::Read;
use std::path::Path;

use esperanto_tidx::Tidx as TidxFile;
use esperanto_txmap::Strand as TxStrand;
use esperanto_txmap::TranscriptRecord;
use esperanto_txmap::TxMap as TxMapFile;

use crate::kmer;
use crate::{CigarOp, Error, Strand, Tidx, TxMap, TxSeqs};

/// Bundle magic.
const MAGIC: &[u8; 8] = b"L1BNDL01";
/// Bundle format version.
const VERSION: u32 = 1;

/// The L1 engine data bundle: a `.tidx` k-mer index plus the
/// full-projection `TxMap` and the transcript sequence store over the same
/// dense `tx_id` space (`transcript_id` lexicographic).
///
/// Projection and sequences cover ALL transcripts — no biotype filtering.
/// Implements [`Tidx`], [`TxMap`] and [`TxSeqs`].
#[derive(Debug, Clone)]
pub struct L1Index {
    tidx: TidxFile,
    txmap: TxMapFile,
    txseqs: Vec<Vec<u8>>,
}

impl L1Index {
    /// Build at runtime from a `.tidx` file, the source GTF, and the
    /// reference FASTA (plain or `.gz`).
    ///
    /// The GTF must parse to exactly the transcript set of the `.tidx`
    /// (same dense id space); otherwise [`Error::Inconsistent`] is
    /// returned. A transcript with any exon beyond its contig gets an
    /// empty sequence (it simply never seeds).
    pub fn build(tidx_path: &Path, gtf_path: &Path, ref_path: &Path) -> Result<Self, Error> {
        let tidx = TidxFile::open(tidx_path)?;
        let set = esperanto_tidx::TranscriptSet::parse(gtf_path)
            .map_err(|e| Error::Tidx(e.to_string()))?;
        if set.len() as u64 != u64::from(tidx.tx_count()) {
            return Err(Error::Inconsistent(format!(
                "GTF has {} transcripts, .tidx has {}",
                set.len(),
                tidx.tx_count()
            )));
        }

        // Full-projection records: exons in transcription order (minus
        // strand reversed), name = transcript_id.
        let mut records = Vec::with_capacity(set.len());
        for i in 0..set.len() {
            let t = set
                .transcript(i)
                .ok_or_else(|| Error::Inconsistent("transcript index out of range".to_string()))?;
            let strand = match t.strand {
                b'+' => TxStrand::Plus,
                b'-' => TxStrand::Minus,
                other => {
                    return Err(Error::Inconsistent(format!(
                        "transcript {} has strand byte {other:?}",
                        t.transcript_id
                    )))
                }
            };
            let mut exons = Vec::with_capacity(t.exons.len());
            for &(start, end) in t.exons {
                let (Ok(start), Ok(end)) = (u32::try_from(start), u32::try_from(end)) else {
                    return Err(Error::Inconsistent(format!(
                        "transcript {} has an exon outside u32 coordinates",
                        t.transcript_id
                    )));
                };
                exons.push(esperanto_txmap::Exon {
                    g_start: start,
                    g_end: end,
                });
            }
            if strand == TxStrand::Minus {
                exons.reverse(); // transcription order runs genomic descending
            }
            records.push(TranscriptRecord {
                name: t.transcript_id.to_string(),
                contig: t.contig.to_string(),
                strand,
                exons,
            });
        }
        // Source hash is opaque here (engine does not depend on a hasher);
        // zeros keep the bundle deterministic.
        let txmap =
            TxMapFile::from_records(records, [0u8; 32]).map_err(|e| Error::TxMap(e.to_string()))?;

        let genome = read_fasta(ref_path)?;
        let mut txseqs = Vec::with_capacity(set.len());
        for i in 0..set.len() {
            let t = set
                .transcript(i)
                .ok_or_else(|| Error::Inconsistent("transcript index out of range".to_string()))?;
            txseqs.push(transcript_sequence(t, &genome));
        }
        Ok(L1Index {
            tidx,
            txmap,
            txseqs,
        })
    }

    /// Open a bundle previously written by [`L1Index::save`].
    ///
    /// The `.tidx` k-mer index is opened from the sibling file
    /// `<bundle-stem>.tidx` (the bundle itself stores only the projection
    /// and the sequences).
    pub fn open(bundle_path: &Path) -> Result<Self, Error> {
        let tidx_path = bundle_path.with_extension("tidx");
        let tidx = TidxFile::open(&tidx_path)?;
        let data = std::fs::read(bundle_path)?;
        if data.len() < 20 {
            return Err(Error::Format("bundle shorter than its header".to_string()));
        }
        let mut found = [0u8; 8];
        found.copy_from_slice(&data[0..8]);
        if found != *MAGIC {
            return Err(Error::BadMagic {
                expected: *MAGIC,
                found,
            });
        }
        let version = u32::from_le_bytes([data[8], data[9], data[10], data[11]]);
        if version != VERSION {
            return Err(Error::UnsupportedVersion(version));
        }
        let blob_len = u64::from_le_bytes([
            data[12], data[13], data[14], data[15], data[16], data[17], data[18], data[19],
        ]) as usize;
        let mut off = 20usize;
        let blob_end = off
            .checked_add(blob_len)
            .ok_or_else(|| Error::Format("txmap section length overflows".to_string()))?;
        let blob = data
            .get(off..blob_end)
            .ok_or_else(|| Error::Format("truncated txmap section".to_string()))?;
        off = blob_end;
        let txmap = TxMapFile::from_bytes(blob).map_err(|e| Error::TxMap(e.to_string()))?;

        let read_u32 = |off: &mut usize| -> Result<u32, Error> {
            let b = data
                .get(*off..*off + 4)
                .ok_or_else(|| Error::Format("truncated txseqs section".to_string()))?;
            *off += 4;
            Ok(u32::from_le_bytes([b[0], b[1], b[2], b[3]]))
        };
        let tx_count = read_u32(&mut off)? as usize;
        if tx_count != txmap.tx_count() {
            return Err(Error::Inconsistent(format!(
                "bundle txseqs count {tx_count} != txmap count {}",
                txmap.tx_count()
            )));
        }
        let mut txseqs = Vec::with_capacity(tx_count);
        for _ in 0..tx_count {
            let len = read_u32(&mut off)? as usize;
            let end = off
                .checked_add(len)
                .ok_or_else(|| Error::Format("txseqs entry length overflows".to_string()))?;
            let seq = data
                .get(off..end)
                .ok_or_else(|| Error::Format("truncated txseqs entry".to_string()))?
                .to_vec();
            off = end;
            txseqs.push(seq);
        }
        if off != data.len() {
            return Err(Error::Format(
                "trailing bytes after the txseqs section".to_string(),
            ));
        }
        Ok(L1Index {
            tidx,
            txmap,
            txseqs,
        })
    }

    /// Write the deterministic bundle: magic `L1BNDL01`, version u32 = 1,
    /// length-prefixed `.txmap` blob, then the txseqs section
    /// (`tx_count` u32, per transcript `len` u32 + bytes). Identical
    /// inputs produce byte-identical files.
    pub fn save(&self, bundle_path: &Path) -> Result<(), Error> {
        let mut out = Vec::new();
        out.extend_from_slice(MAGIC);
        out.extend_from_slice(&VERSION.to_le_bytes());
        let blob = self.txmap.to_bytes();
        out.extend_from_slice(&(blob.len() as u64).to_le_bytes());
        out.extend_from_slice(&blob);
        out.extend_from_slice(&(self.txseqs.len() as u32).to_le_bytes());
        for seq in &self.txseqs {
            out.extend_from_slice(&(seq.len() as u32).to_le_bytes());
            out.extend_from_slice(seq);
        }
        std::fs::write(bundle_path, out)?;
        Ok(())
    }

    /// Contig name for a contig id (projection name table).
    pub fn contig_name(&self, contig_id: u32) -> Option<&str> {
        self.txmap
            .contigs()
            .get(contig_id as usize)
            .map(String::as_str)
    }

    /// The underlying projection map (for tooling).
    pub fn txmap(&self) -> &TxMapFile {
        &self.txmap
    }
}

impl Tidx for L1Index {
    fn k(&self) -> u32 {
        self.tidx.k()
    }

    fn tx_count(&self) -> u32 {
        self.tidx.tx_count()
    }

    fn lookup(&self, canonical_kmer: u64) -> &[(u32, u32)] {
        self.tidx.lookup(canonical_kmer)
    }

    fn transcript_name(&self, tx_id: u32) -> &str {
        self.tidx.transcript_name(tx_id)
    }
}

/// Convert a projection CIGAR op into the engine's op set.
fn convert_op(op: esperanto_txmap::CigarOp) -> CigarOp {
    match op {
        esperanto_txmap::CigarOp::Match(n) => CigarOp::Match(n),
        esperanto_txmap::CigarOp::RefSkip(n) => CigarOp::RefSkip(n),
    }
}

impl TxMap for L1Index {
    fn project(&self, tx_id: u32, tx_start: u32, len: u32) -> Option<(u32, u32, Vec<CigarOp>)> {
        self.txmap
            .project(tx_id, tx_start, len)
            .map(|(contig, pos, ops)| (contig, pos, ops.into_iter().map(convert_op).collect()))
    }

    fn tx_len(&self, tx_id: u32) -> Option<u32> {
        u32::try_from(self.txmap.tx_len(tx_id)?).ok()
    }

    fn strand(&self, tx_id: u32) -> Option<Strand> {
        self.txmap.strand(tx_id).map(|s| match s {
            TxStrand::Plus => Strand::Plus,
            TxStrand::Minus => Strand::Minus,
        })
    }
}

impl TxSeqs for L1Index {
    fn seq(&self, tx_id: u32) -> &[u8] {
        self.txseqs
            .get(tx_id as usize)
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }
}

/// Read a FASTA (plain or `.gz`) into memory, preserving sequence case.
fn read_fasta(path: &Path) -> Result<BTreeMap<String, Vec<u8>>, Error> {
    let file = std::fs::File::open(path)?;
    let mut text = Vec::new();
    if path
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("gz"))
    {
        flate2::read::GzDecoder::new(file).read_to_end(&mut text)?;
    } else {
        let mut file = file;
        file.read_to_end(&mut text)?;
    }
    let mut genome: BTreeMap<String, Vec<u8>> = BTreeMap::new();
    let mut name: Option<String> = None;
    let mut cur: Vec<u8> = Vec::new();
    for line in text.split(|&b| b == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.first() == Some(&b'>') {
            if let Some(n) = name.take() {
                genome.insert(n, std::mem::take(&mut cur));
            }
            let header = String::from_utf8_lossy(&line[1..]).into_owned();
            name = Some(header.split_whitespace().next().unwrap_or("").to_string());
        } else if name.is_some() {
            cur.extend_from_slice(line);
        }
    }
    if let Some(n) = name {
        genome.insert(n, cur);
    }
    Ok(genome)
}

/// Splice one transcript's sequence: exons concatenated in genomic
/// ascending order, whole reverse-complement on the minus strand, empty
/// when any exon leaves the contig or the contig is unknown.
fn transcript_sequence(
    t: esperanto_tidx::GtfTranscript<'_>,
    genome: &BTreeMap<String, Vec<u8>>,
) -> Vec<u8> {
    let Some(contig) = genome.get(t.contig) else {
        return Vec::new();
    };
    if t.exons.iter().any(|&(_, end)| end > contig.len() as u64) {
        return Vec::new();
    }
    let mut seq = Vec::new();
    for &(start, end) in t.exons {
        seq.extend_from_slice(&contig[start as usize..end as usize]);
    }
    if t.strand == b'-' {
        kmer::revcomp(&seq)
    } else {
        seq
    }
}
