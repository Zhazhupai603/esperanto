//! GENCODE GTF ingestion, ported to match the reference tool's parser:
//!
//! * Both `transcript` and `exon` rows are consumed. A transcript only exists
//!   if it has a `transcript` row (metadata) plus at least one exon.
//! * `transcript_type` allowlist: exactly `protein_coding` and `lncRNA`
//!   (falling back to `gene_type` when `transcript_type` is absent) — a
//!   documented biotype choice, authoritative per-transcript.
//! * Record name = `transcript_id` (NOT `transcript_name`).
//! * Exons are ordered by `exon_number` (transcription order; ascending on
//!   plus, descending genomically on minus), deduplicated by identical span,
//!   with a genomic-order fallback when `exon_number` is unreliable.
//! * Strand `.` is accepted and treated as plus.
//! * Coordinates: 1-based inclusive input converted to 0-based half-open.

use std::collections::BTreeMap;
use std::path::Path;

use crate::transcript::{Exon, Strand, TranscriptRecord};
use crate::Error;

/// Transcript biotypes retained by the reference parser.
pub const ALLOWED_TYPES: &[&str] = &["protein_coding", "lncRNA"];

struct TxMeta {
    contig: String,
    strand: Strand,
    tx_type: String,
}

#[derive(Clone, Copy)]
struct RawExon {
    g_start: u32,
    g_end: u32,
    exon_number: u32,
}

/// Parse a GTF file into records (sorted by name = transcript_id) plus the
/// BLAKE3 hash of the raw file bytes.
pub(crate) fn parse_gtf(path: &Path) -> Result<(Vec<TranscriptRecord>, [u8; 32]), Error> {
    let raw = std::fs::read(path)?;
    let source_hash: [u8; 32] = blake3::hash(&raw).into();
    let text = std::str::from_utf8(&raw)
        .map_err(|e| Error::Format(format!("GTF file is not valid UTF-8: {e}")))?;

    let mut meta: BTreeMap<String, TxMeta> = BTreeMap::new();
    let mut exons: BTreeMap<String, Vec<RawExon>> = BTreeMap::new();

    for (idx, raw_line) in text.lines().enumerate() {
        let lineno = idx + 1;
        let line = raw_line.trim_end_matches('\r');
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            return Err(Error::Format(format!(
                "line {lineno}: expected at least 9 tab-separated fields"
            )));
        }
        let contig = fields[0];
        let feature = fields[2];
        if feature != "transcript" && feature != "exon" {
            continue;
        }
        let start: u64 = fields[3]
            .parse()
            .map_err(|_| Error::Format(format!("line {lineno}: invalid start '{}'", fields[3])))?;
        let end: u64 = fields[4]
            .parse()
            .map_err(|_| Error::Format(format!("line {lineno}: invalid end '{}'", fields[4])))?;
        if start < 1 {
            return Err(Error::Format(format!(
                "line {lineno}: start must be >= 1 (1-based inclusive)"
            )));
        }
        if end < start {
            return Err(Error::Format(format!(
                "line {lineno}: end {end} < start {start}"
            )));
        }
        if end > u32::MAX as u64 {
            return Err(Error::Format(format!(
                "line {lineno}: end {end} exceeds u32 coordinate range"
            )));
        }
        let strand = match fields[6] {
            "+" | "." => Strand::Plus,
            "-" => Strand::Minus,
            other => {
                return Err(Error::Format(format!(
                    "line {lineno}: unsupported strand '{other}'"
                )))
            }
        };
        let attrs = fields[8];
        let tx_id = parse_attr(attrs, "transcript_id")
            .ok_or_else(|| {
                Error::Format(format!("line {lineno}: missing transcript_id attribute"))
            })?
            .to_string();
        match feature {
            "transcript" => {
                let tx_type = parse_attr(attrs, "transcript_type")
                    .or_else(|| parse_attr(attrs, "gene_type"))
                    .unwrap_or("")
                    .to_string();
                // The transcript row is authoritative; refresh on repeats.
                let entry = meta.entry(tx_id).or_insert(TxMeta {
                    contig: String::new(),
                    strand,
                    tx_type: String::new(),
                });
                entry.contig = contig.to_string();
                entry.strand = strand;
                entry.tx_type = tx_type;
            }
            _ => {
                let exon_number = parse_attr(attrs, "exon_number")
                    .and_then(|s| s.parse::<u32>().ok())
                    .unwrap_or(0);
                exons.entry(tx_id).or_default().push(RawExon {
                    g_start: (start - 1) as u32,
                    g_end: end as u32,
                    exon_number,
                });
            }
        }
    }

    let mut records: Vec<TranscriptRecord> = Vec::new();
    for (id, m) in &meta {
        if !ALLOWED_TYPES.contains(&m.tx_type.as_str()) {
            continue;
        }
        let Some(raws) = exons.get(id) else {
            continue;
        };
        if raws.is_empty() {
            continue;
        }
        records.push(TranscriptRecord {
            name: id.clone(),
            contig: m.contig.clone(),
            strand: m.strand,
            exons: assemble_exons(raws, m.strand)?,
        });
    }
    records.sort_unstable_by(|a, b| a.name.cmp(&b.name));
    Ok((records, source_hash))
}

/// Order raw exons into transcription order: `exon_number` ascending (which is
/// genomic-ascending on plus, genomic-descending on minus), dedup identical
/// spans, with a genomic-order fallback when the numbers are not monotone.
fn assemble_exons(raws: &[RawExon], strand: Strand) -> Result<Vec<Exon>, Error> {
    let mut sorted: Vec<&RawExon> = raws.iter().collect();
    sorted.sort_by_key(|r| r.exon_number);
    sorted.dedup_by(|a, b| a.g_start == b.g_start && a.g_end == b.g_end);
    if sorted.is_empty() {
        return Err(Error::Format("transcript has no exons after dedup".into()));
    }
    let monotonic = match strand {
        Strand::Plus => sorted.windows(2).all(|w| w[0].g_start <= w[1].g_start),
        Strand::Minus => sorted.windows(2).all(|w| w[0].g_start >= w[1].g_start),
    };
    if !monotonic {
        sorted.sort_by_key(|r| r.g_start);
        if strand == Strand::Minus {
            sorted.reverse();
        }
    }
    Ok(sorted
        .iter()
        .map(|r| Exon {
            g_start: r.g_start,
            g_end: r.g_end,
        })
        .collect())
}

/// Extract `key "value";` (value optionally unquoted) from a GTF attribute
/// column. The key must be followed by whitespace or `=` so `transcript_id`
/// never matches `transcript_idfoo`.
fn parse_attr<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    for part in attrs.split(';') {
        let part = part.trim();
        if let Some(rest) = part.strip_prefix(key) {
            if rest.starts_with(' ') || rest.starts_with('\t') || rest.starts_with('=') {
                let value = rest.trim().trim_matches('"');
                if !value.is_empty() {
                    return Some(value);
                }
            }
        }
    }
    None
}
