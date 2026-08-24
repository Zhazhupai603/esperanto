//! [`TxMap`]: the transcript index plus the deterministic little-endian
//! `.txmap` binary format.
//!
//! Format (all integers little-endian, no timestamps, no padding; every
//! collection is written in sorted order so identical inputs produce
//! byte-identical files):
//!
//! ```text
//! magic        8 bytes "TXMAP001"
//! version      u32 = 1
//! source_hash [u8;32]   BLAKE3 of the source GTF bytes
//! tx_count     u32
//! contig_count u32
//!   repeat: name_len u32, name UTF-8          (contigs lexicographic)
//! junction_count u32
//!   repeat: contig_id u32, intron_start u32, intron_end u32   (sorted)
//! transcripts (lexicographic by name; index = tx_id):
//!   repeat: name_len u32, name; contig_id u32; strand u8 (0=+, 1=-);
//!           exon_count u32; repeat exon: g_start u32, g_end u32
//! ```

use std::collections::HashMap;
use std::path::Path;

use crate::attribution::{attribute, Attribution, AttributionCandidate};
use crate::cigar::CigarOp;
use crate::junction::{Junction, JunctionSet};
use crate::transcript::{Exon, Strand, TranscriptRecord};
use crate::Error;

const MAGIC: &[u8; 8] = b"TXMAP001";
const VERSION: u32 = 1;

/// Transcript index over a reference transcriptome.
///
/// Invariants: `transcripts` sorted by name (index = tx_id, names unique),
/// `contigs` sorted lexicographically (index = contig_id), `junctions` sorted
/// and deduplicated, and serialization is byte-deterministic.
#[derive(Clone, Debug)]
pub struct TxMap {
    contigs: Vec<String>,
    transcripts: Vec<TranscriptRecord>,
    junctions: JunctionSet,
    source_hash: [u8; 32],
    tx_index: HashMap<String, u32>,
    contig_index: HashMap<String, u32>,
}

impl TxMap {
    /// Build from transcript records (any order; records are validated,
    /// then sorted by name). `source_hash` is opaque here; use
    /// [`TxMap::from_gtf`] to get the BLAKE3 of the source GTF bytes.
    pub fn from_records(
        records: Vec<TranscriptRecord>,
        source_hash: [u8; 32],
    ) -> Result<Self, Error> {
        for record in &records {
            record.validate()?;
        }
        let mut records = records;
        records.sort_unstable_by(|a, b| a.name.cmp(&b.name));
        for pair in records.windows(2) {
            if pair[0].name == pair[1].name {
                return Err(Error::Format(format!(
                    "duplicate transcript name '{}'",
                    pair[0].name
                )));
            }
        }
        // Counts and name lengths must fit the u32 fields of the format.
        u32::try_from(records.len())
            .map_err(|_| Error::Format("transcript count exceeds u32".to_string()))?;
        for record in &records {
            u32::try_from(record.name.len()).map_err(|_| {
                Error::Format(format!("transcript name too long: '{}'", record.name))
            })?;
            u32::try_from(record.exons.len()).map_err(|_| {
                Error::Format(format!("exon count exceeds u32 for '{}'", record.name))
            })?;
        }

        let mut contigs: Vec<String> = records.iter().map(|r| r.contig.clone()).collect();
        contigs.sort();
        contigs.dedup();
        u32::try_from(contigs.len())
            .map_err(|_| Error::Format("contig count exceeds u32".to_string()))?;
        for name in &contigs {
            u32::try_from(name.len())
                .map_err(|_| Error::Format(format!("contig name too long: '{name}'")))?;
        }

        let contig_index: HashMap<String, u32> = contigs
            .iter()
            .enumerate()
            .map(|(i, name)| (name.clone(), i as u32))
            .collect();
        let mut junction_list: Vec<Junction> = Vec::new();
        for record in &records {
            let contig_id = contig_index[&record.contig]; // invariant: from records
            for (start, end) in record.introns() {
                junction_list.push(Junction {
                    contig_id,
                    start,
                    end,
                });
            }
        }
        let junctions = JunctionSet::from_junctions(junction_list);
        u32::try_from(junctions.len())
            .map_err(|_| Error::Format("junction count exceeds u32".to_string()))?;

        let tx_index: HashMap<String, u32> = records
            .iter()
            .enumerate()
            .map(|(i, record)| (record.name.clone(), i as u32))
            .collect();

        Ok(Self {
            contigs,
            transcripts: records,
            junctions,
            source_hash,
            tx_index,
            contig_index,
        })
    }

    /// Parse a GENCODE-flavored GTF file and build the index; the source
    /// hash is the BLAKE3 of the raw GTF bytes.
    pub fn from_gtf(path: &Path) -> Result<Self, Error> {
        let (records, hash) = crate::gtf::parse_gtf(path)?;
        Self::from_records(records, hash)
    }

    /// Serialize to the deterministic `.txmap` byte stream.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(MAGIC);
        out.extend_from_slice(&VERSION.to_le_bytes());
        out.extend_from_slice(&self.source_hash);
        put_u32(&mut out, self.transcripts.len() as u32);
        put_u32(&mut out, self.contigs.len() as u32);
        for name in &self.contigs {
            put_str(&mut out, name);
        }
        put_u32(&mut out, self.junctions.len() as u32);
        for junction in self.junctions.junctions() {
            put_u32(&mut out, junction.contig_id);
            put_u32(&mut out, junction.start);
            put_u32(&mut out, junction.end);
        }
        for record in &self.transcripts {
            put_str(&mut out, &record.name);
            put_u32(&mut out, self.contig_index[&record.contig]); // invariant
            out.push(match record.strand {
                Strand::Plus => 0u8,
                Strand::Minus => 1u8,
            });
            put_u32(&mut out, record.exons.len() as u32);
            for exon in &record.exons {
                put_u32(&mut out, exon.g_start);
                put_u32(&mut out, exon.g_end);
            }
        }
        out
    }

    /// Parse the `.txmap` byte stream (strict: rejects bad magic/version,
    /// truncated or trailing bytes, out-of-range ids, and records failing
    /// [`TranscriptRecord::validate`]).
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
        let mut reader = Reader { buf: bytes, pos: 0 };
        let magic = reader.bytes(8)?;
        if magic != MAGIC.as_slice() {
            let mut found = [0u8; 8];
            found.copy_from_slice(magic); // len 8 guaranteed by bytes()
            return Err(Error::Magic {
                expected: *MAGIC,
                found,
            });
        }
        let version = reader.u32()?;
        if version != VERSION {
            return Err(Error::Version {
                file: version,
                code: VERSION,
            });
        }
        let mut source_hash = [0u8; 32];
        source_hash.copy_from_slice(reader.bytes(32)?); // len 32 guaranteed

        let tx_count = reader.u32()? as usize;
        let contig_count = reader.u32()? as usize;
        if contig_count.saturating_mul(4) > reader.remaining() {
            return Err(Error::Format("contig count exceeds file size".to_string()));
        }
        let mut contigs: Vec<String> = Vec::new();
        for _ in 0..contig_count {
            contigs.push(reader.string()?);
        }
        for pair in contigs.windows(2) {
            if pair[0] >= pair[1] {
                return Err(Error::Format(
                    "contig table not strictly sorted / unique".to_string(),
                ));
            }
        }

        let junction_count = reader.u32()? as usize;
        if junction_count.saturating_mul(12) > reader.remaining() {
            return Err(Error::Format(
                "junction count exceeds file size".to_string(),
            ));
        }
        let mut junction_list: Vec<Junction> = Vec::new();
        for _ in 0..junction_count {
            let contig_id = reader.u32()?;
            let start = reader.u32()?;
            let end = reader.u32()?;
            if contig_id as usize >= contigs.len() {
                return Err(Error::Format(format!(
                    "junction contig id {contig_id} out of range"
                )));
            }
            // Legacy files written by the old tool may contain inverted
            // `(end, start)` entries for minus-strand transcripts; normalize
            // on load (new files always store forward intervals).
            let (start, end) = if start <= end { (start, end) } else { (end, start) };
            junction_list.push(Junction {
                contig_id,
                start,
                end,
            });
        }

        if tx_count.saturating_mul(13) > reader.remaining() {
            return Err(Error::Format("tx count exceeds file size".to_string()));
        }
        let mut transcripts: Vec<TranscriptRecord> = Vec::new();
        for _ in 0..tx_count {
            let name = reader.string()?;
            let contig_id = reader.u32()?;
            if contig_id as usize >= contigs.len() {
                return Err(Error::Format(format!(
                    "transcript '{name}' contig id {contig_id} out of range"
                )));
            }
            let strand = match reader.u8()? {
                0 => Strand::Plus,
                1 => Strand::Minus,
                other => {
                    return Err(Error::Format(format!(
                        "transcript '{name}': invalid strand byte {other}"
                    )))
                }
            };
            let exon_count = reader.u32()? as usize;
            if exon_count.saturating_mul(8) > reader.remaining() {
                return Err(Error::Format(format!(
                    "transcript '{name}': exon count exceeds file size"
                )));
            }
            let mut exons: Vec<Exon> = Vec::with_capacity(exon_count);
            for _ in 0..exon_count {
                let g_start = reader.u32()?;
                let g_end = reader.u32()?;
                exons.push(Exon { g_start, g_end });
            }
            let record = TranscriptRecord {
                name,
                contig: contigs[contig_id as usize].clone(),
                strand,
                exons,
            };
            record.validate()?;
            transcripts.push(record);
        }
        for pair in transcripts.windows(2) {
            if pair[0].name >= pair[1].name {
                return Err(Error::Format(
                    "transcript table not strictly sorted by name".to_string(),
                ));
            }
        }
        if reader.remaining() != 0 {
            return Err(Error::Format(
                "trailing bytes after transcript table".to_string(),
            ));
        }

        let contig_index: HashMap<String, u32> = contigs
            .iter()
            .enumerate()
            .map(|(i, name)| (name.clone(), i as u32))
            .collect();
        let tx_index: HashMap<String, u32> = transcripts
            .iter()
            .enumerate()
            .map(|(i, record)| (record.name.clone(), i as u32))
            .collect();
        Ok(Self {
            contigs,
            transcripts,
            junctions: JunctionSet::from_junctions(junction_list),
            source_hash,
            tx_index,
            contig_index,
        })
    }

    /// Write to `path` (`.txmap` binary format).
    pub fn save(&self, path: &Path) -> Result<(), Error> {
        std::fs::write(path, self.to_bytes())?;
        Ok(())
    }

    /// Read a `.txmap` file.
    pub fn open(path: &Path) -> Result<Self, Error> {
        let bytes = std::fs::read(path)?;
        Self::from_bytes(&bytes)
    }

    /// Number of transcripts.
    pub fn tx_count(&self) -> usize {
        self.transcripts.len()
    }

    /// Look up a transcript id by name.
    pub fn tx_id(&self, name: &str) -> Option<u32> {
        self.tx_index.get(name).copied()
    }

    /// Name of transcript `tx_id`.
    pub fn tx_name(&self, tx_id: u32) -> Option<&str> {
        self.transcripts
            .get(tx_id as usize)
            .map(|t| t.name.as_str())
    }

    /// Transcript length in bases.
    pub fn tx_len(&self, tx_id: u32) -> Option<u64> {
        self.transcripts.get(tx_id as usize).map(|t| t.tx_len())
    }

    /// Annotated strand of transcript `tx_id`.
    pub fn strand(&self, tx_id: u32) -> Option<Strand> {
        self.transcripts.get(tx_id as usize).map(|t| t.strand)
    }

    /// Contig name of transcript `tx_id`.
    pub fn contig_of(&self, tx_id: u32) -> Option<&str> {
        self.transcripts
            .get(tx_id as usize)
            .map(|t| t.contig.as_str())
    }

    /// All contig names, lexicographically sorted (index = contig_id).
    pub fn contigs(&self) -> &[String] {
        &self.contigs
    }

    /// All transcript records, sorted by name (index = tx_id).
    pub fn transcripts(&self) -> &[TranscriptRecord] {
        &self.transcripts
    }

    /// The known-junction set.
    pub fn junctions(&self) -> &JunctionSet {
        &self.junctions
    }

    /// BLAKE3 hash of the source GTF bytes.
    pub fn source_hash(&self) -> &[u8; 32] {
        &self.source_hash
    }

    /// Project transcript interval `[tx_start, tx_start + len)` of `tx_id`
    /// onto the forward genome strand. Returns
    /// `(contig_id, genomic_start, cigar)`; `None` when `tx_id` is unknown,
    /// `len == 0`, or `tx_start >= tx_len`. See
    /// [`TranscriptRecord::project`] for clamping semantics.
    pub fn project(&self, tx_id: u32, tx_start: u32, len: u32) -> Option<(u32, u32, Vec<CigarOp>)> {
        let record = self.transcripts.get(tx_id as usize)?;
        let contig_id = self.contig_index[&record.contig]; // invariant: from records
        let (genomic_start, cigar) = record.project(tx_start, len)?;
        Some((contig_id, genomic_start, cigar))
    }

    /// Attribute a read's candidate projections to a single placement using
    /// this map's junction set and contig table. See [`attribute`].
    pub fn attribute(&self, candidates: &[AttributionCandidate]) -> Attribution {
        attribute(candidates, &self.junctions, &self.contigs)
    }
}

fn put_u32(out: &mut Vec<u8>, value: u32) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn put_str(out: &mut Vec<u8>, value: &str) {
    put_u32(out, value.len() as u32);
    out.extend_from_slice(value.as_bytes());
}

/// Bounds-checked cursor over the input bytes; every read either returns
/// exactly the requested bytes or `Error::Format`.
struct Reader<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    fn remaining(&self) -> usize {
        self.buf.len() - self.pos
    }

    fn u8(&mut self) -> Result<u8, Error> {
        let byte = self
            .buf
            .get(self.pos)
            .copied()
            .ok_or_else(|| Error::Format("truncated file".to_string()))?;
        self.pos += 1;
        Ok(byte)
    }

    fn u32(&mut self) -> Result<u32, Error> {
        let bytes = self.bytes(4)?;
        let mut array = [0u8; 4];
        array.copy_from_slice(bytes);
        Ok(u32::from_le_bytes(array))
    }

    fn bytes(&mut self, n: usize) -> Result<&'a [u8], Error> {
        let end = self
            .pos
            .checked_add(n)
            .ok_or_else(|| Error::Format("file length overflow".to_string()))?;
        let slice = self
            .buf
            .get(self.pos..end)
            .ok_or_else(|| Error::Format("truncated file".to_string()))?;
        self.pos = end;
        Ok(slice)
    }

    fn string(&mut self) -> Result<String, Error> {
        let len = self.u32()? as usize;
        let bytes = self.bytes(len)?;
        String::from_utf8(bytes.to_vec())
            .map_err(|_| Error::Format("name field is not valid UTF-8".to_string()))
    }
}
