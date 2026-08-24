//! `.baln` — lightweight binary alignment format for the internal align→call
//! channel. BAM stays the user-facing artifact; .baln is the fast internal
//! interface.
//!
//! Layout (all little-endian):
//!   header: magic[8]="ESPBALN\x01" | contig_count u32 | per contig: name_len u8, name
//!   record: block_size u32 | bamio encode::try_encode bytes (32-byte core + data)
//!
//! Each record is the exact BAM on-disk record layout, just without BGZF
//! framing — the call side memcpys straight into `bam1_t`.

use std::io::{self, Write};

use crate::bam::BamRecord;

pub const MAGIC: &[u8; 8] = b"ESPBALN\x01";

/// Write the .baln header (contig name table).
pub fn write_header(mut w: impl Write, contig_names: &[String]) -> io::Result<()> {
    w.write_all(MAGIC)?;
    w.write_all(&(contig_names.len() as u32).to_le_bytes())?;
    for name in contig_names {
        w.write_all(&(name.len() as u8).to_le_bytes())?;
        w.write_all(name.as_bytes())?;
    }
    Ok(())
}

/// Write one record: block_size + BAM record bytes (core+data, no BGZF).
///
/// Fast-path decline (rare: exotic name / qual>93 / cigar overflow) surfaces as
/// an error — never silently drop a record.
pub fn write_record(mut w: impl Write, rec: &BamRecord) -> io::Result<()> {
    let mut buf = Vec::with_capacity(256);
    match esperanto_bamio::encode::try_encode(&mut buf, rec) {
        Some(Ok(())) => {
            w.write_all(&(buf.len() as u32).to_le_bytes())?;
            w.write_all(&buf)
        }
        Some(Err(e)) => Err(e),
        None => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("baln: fast-path declined for record {}", rec.name),
        )),
    }
}
