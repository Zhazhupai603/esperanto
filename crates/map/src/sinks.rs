//! Output sinks for the map serial write section.
//!
//! - [`RecordSinks`]: BAM (BGZF) + `.baln` sharing ONE binary encode per
//!   record. The `.baln` record is exactly the BAM record bytes (`block_size`
//!   prefix included) without BGZF framing, so the whole encode — CIGAR
//!   serialization, reg2bin, 4-bit SEQ packing, tag rendering — runs once per
//!   record and the byte string is written to whichever sinks are present.
//! - [`UnmappedFq`]: gzipped unmapped FASTQ (level 1, rescue-channel contract)
//!   with an 8 MiB staging buffer in front of the file and one encoder write
//!   per record.

use std::io::{self, Write};

use flate2::write::GzEncoder;
use flate2::Compression;
use noodles::sam::Header;

use crate::bam::{self, BamRecord};
use crate::error::AlignError;
use crate::fastq::FastqRecord;
use crate::index::Index;

fn io_err(e: io::Error) -> AlignError {
    AlignError::FastqIo {
        path: "<bam>".into(),
        source: e,
    }
}

fn io_err_unm(e: io::Error) -> AlignError {
    AlignError::FastqIo {
        path: "<unmapped.fq.gz>".into(),
        source: e,
    }
}

/// BAM + `.baln` writers with a shared per-record encode buffer.
pub(crate) struct RecordSinks {
    bam: Option<
        noodles::bam::io::Writer<noodles::bgzf::io::MultithreadedWriter<Box<dyn Write + Send>>>,
    >,
    baln: Option<std::io::BufWriter<Box<dyn Write + Send>>>,
    /// Retained for the rare fast-path decline that must fall back to the
    /// noodles RecordBuf builder.
    header: Header,
    buf: Vec<u8>,
}

impl RecordSinks {
    /// Take the pipeline's sinks (`out.bam.take()` / `out.baln.take()`) and
    /// write both headers, in the same order the pipeline historically did:
    /// BAM header into BGZF first, then the `.baln` contig-name table.
    pub(crate) fn new(
        bam_sink: Option<Box<dyn Write + Send>>,
        baln_sink: Option<Box<dyn Write + Send>>,
        index: &Index,
        threads: usize,
    ) -> Result<Self, AlignError> {
        let header = bam::build_header(index);
        let bam = match bam_sink {
            Some(sink) => Some(bam::create_writer(sink, &header, threads).map_err(io_err)?),
            None => None,
        };
        let baln = match baln_sink {
            Some(w) => {
                let mut bw = std::io::BufWriter::with_capacity(8 << 20, w);
                let names: Vec<String> = index
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
        Ok(RecordSinks {
            bam,
            baln,
            header,
            buf: Vec::with_capacity(8192),
        })
    }

    /// Encode `rec` once and write it to every present sink.
    ///
    /// Record order and bytes are identical to writing each sink separately:
    /// both formats store `block_size u32 | record bytes`, the BAM writer just
    /// wraps the same byte string in BGZF.
    pub(crate) fn write(&mut self, rec: &BamRecord) -> Result<(), AlignError> {
        self.buf.clear();
        self.buf.extend_from_slice(&0u32.to_le_bytes()); // block_size placeholder
        match esperanto_bamio::encode::try_encode(&mut self.buf, rec) {
            Some(Ok(())) => {
                let block_size = u32::try_from(self.buf.len() - 4)
                    .map_err(|e| io_err(io::Error::new(io::ErrorKind::InvalidInput, e)))?;
                self.buf[..4].copy_from_slice(&block_size.to_le_bytes());
                if let Some(w) = &mut self.bam {
                    w.get_mut().write_all(&self.buf).map_err(io_err)?;
                }
                if let Some(w) = &mut self.baln {
                    w.write_all(&self.buf).map_err(io_err)?;
                }
                Ok(())
            }
            Some(Err(e)) => Err(io_err(e)),
            None => {
                // Fast-path decline (exotic name / qual>93 / cigar overflow):
                // BAM falls back to the noodles RecordBuf builder; .baln has
                // no slow path and surfaces the decline as an error — never
                // silently drop a record. BAM write first, matching the
                // historical per-sink write order.
                if let Some(w) = &mut self.bam {
                    esperanto_bamio::write_record_slow(w, &self.header, rec).map_err(io_err)?;
                }
                if self.baln.is_some() {
                    return Err(io_err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "baln: fast-path declined for record {}",
                            String::from_utf8_lossy(&rec.name)
                        ),
                    )));
                }
                Ok(())
            }
        }
    }

    /// Close both sinks (BGZF EOF block / final `.baln` flush).
    pub(crate) fn finish(self) -> Result<(), AlignError> {
        if let Some(w) = self.bam {
            w.into_inner().finish().map_err(io_err)?;
        }
        if let Some(mut bw) = self.baln {
            std::io::Write::flush(&mut bw).map_err(io_err)?;
        }
        Ok(())
    }
}

/// Gzipped unmapped-FASTQ sink: 8 MiB staging in front of the encoder's
/// output, one encoder write per record.
pub(crate) struct UnmappedFq<W: Write> {
    gz: GzEncoder<std::io::BufWriter<W>>,
    scratch: Vec<u8>,
}

impl<W: Write> UnmappedFq<W> {
    pub(crate) fn new(inner: W) -> Self {
        UnmappedFq {
            gz: GzEncoder::new(
                std::io::BufWriter::with_capacity(8 << 20, inner),
                Compression::new(1),
            ),
            scratch: Vec::with_capacity(1024),
        }
    }

    /// Write one unmapped record; `mate` carries the `/1` / `/2` marker
    /// (PE rescue-channel contract: the rescue stage re-parses these
    /// line-wise and mates records by that marker).
    pub(crate) fn write(&mut self, rec: &FastqRecord, mate: Option<u8>) -> Result<(), AlignError> {
        render_unmapped(&mut self.scratch, rec, mate);
        self.gz.write_all(&self.scratch).map_err(io_err_unm)
    }

    /// Flush the deflate stream and the staging buffer. Always called, even
    /// for zero records (an empty gzip stream is part of the artifact
    /// contract).
    pub(crate) fn finish(self) -> Result<(), AlignError> {
        let mut buf = self.gz.finish().map_err(io_err_unm)?;
        buf.flush().map_err(io_err_unm)
    }
}

/// Render one unmapped FASTQ record into `out` (four-line well-formed form).
pub(crate) fn render_unmapped(out: &mut Vec<u8>, rec: &FastqRecord, mate: Option<u8>) {
    out.clear();
    out.reserve(rec.name.len() + rec.seq.len() + rec.qual.len() + 8);
    out.push(b'@');
    out.extend_from_slice(&rec.name);
    if let Some(m) = mate {
        out.push(b'/');
        out.push(b'0' + m);
    }
    out.push(b'\n');
    out.extend_from_slice(&rec.seq);
    out.extend_from_slice(b"\n+\n");
    out.extend_from_slice(&rec.qual);
    out.push(b'\n');
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fastq::FastqRecord;
    use std::io::Read as _;
    use std::sync::{Arc, Mutex};

    /// Collecting sink shared with the writer (sinks take ownership).
    #[derive(Clone)]
    struct SharedSink(Arc<Mutex<Vec<u8>>>);

    impl Write for SharedSink {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            self.0.lock().expect("sink lock").extend_from_slice(buf);
            Ok(buf.len())
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    fn rec(name: &str, seq: &str, qual: &str) -> FastqRecord {
        FastqRecord {
            name: name.as_bytes().to_vec(),
            seq: seq.as_bytes().to_vec(),
            qual: qual.as_bytes().to_vec(),
        }
    }

    /// Four-line well-formed output, with and without the mate marker — the
    /// rescue stage re-parses unmapped.fq.gz line-wise.
    #[test]
    fn render_unmapped_is_four_line_wellformed() {
        let r = rec("READ", "ACGT", "IIII");
        let mut buf = Vec::new();
        render_unmapped(&mut buf, &r, None);
        assert_eq!(buf, b"@READ\nACGT\n+\nIIII\n");
        render_unmapped(&mut buf, &r, Some(1));
        assert_eq!(buf, b"@READ/1\nACGT\n+\nIIII\n");
        render_unmapped(&mut buf, &r, Some(2));
        assert_eq!(buf, b"@READ/2\nACGT\n+\nIIII\n");
    }

    /// UnmappedFq content roundtrip through gzip: exact bytes, /1 /2 markers
    /// intact, stream readable as a plain gzip file.
    #[test]
    fn unmapped_fq_roundtrip_exact_content() {
        let cell = Arc::new(Mutex::new(Vec::new()));
        let mut w = UnmappedFq::new(SharedSink(cell.clone()));
        w.write(&rec("a", "ACGT", "IIII"), None).expect("write");
        w.write(&rec("b", "TTTT", "JJJJ"), Some(2)).expect("write");
        w.finish().expect("finish");
        let raw = Arc::try_unwrap(cell)
            .map_err(|_| "still shared")
            .expect("unique")
            .into_inner()
            .expect("sink lock");
        let mut txt = Vec::new();
        flate2::read::GzDecoder::new(&raw[..])
            .read_to_end(&mut txt)
            .expect("gunzip");
        assert_eq!(txt, b"@a\nACGT\n+\nIIII\n@b/2\nTTTT\n+\nJJJJ\n");
    }

    fn tiny_index() -> Index {
        let fasta = b">chrT\nACGTACGTACGTACGTACGT\n".to_vec();
        let reference = crate::fasta::parse_fasta_bytes(&fasta).expect("fasta");
        Index::build(reference, crate::seed::SeedParams { k: 5, w: 2 })
    }

    fn synthetic_records() -> Vec<BamRecord> {
        use crate::extend::CigarOp;
        use crate::mapq::ReadAlignment;
        use crate::seed::Strand;
        let aln = |strand: Strand, pos: u32| ReadAlignment {
            contig: 0,
            pos,
            strand,
            score: 140,
            second_score: 0,
            chain_score: 140,
            second_chain_score: 0,
            cigar: vec![
                CigarOp::SoftClip(3),
                CigarOp::Match(20),
                CigarOp::RefSkip(79),
                CigarOp::Match(27),
            ],
            n_anchors: 12,
            junctions: Vec::new(),
            ea_ag: 2,
            ea_tc: 0,
            mm_count: 1,
            n_seeds: 14,
            rescued: false,
            mapq_override: None,
            origin: crate::mapq::Origin::Genomic,
            n_cands: 1,
        };
        let seq: Vec<u8> = b"ACGT".repeat(12);
        let qual: Vec<u8> = vec![b'I'; seq.len()];
        let mut out = vec![
            crate::bam::record_se(b"plus", &seq, &qual, Some(&aln(Strand::Plus, 10))),
            crate::bam::record_se(b"minus", &seq, &qual, Some(&aln(Strand::Minus, 4))),
            crate::bam::record_se(b"unmapped", &seq, &qual, None),
        ];
        // one PE-shaped record (mate fields + tlen) to cover the mate encode
        out[0].mate = Some((0, 55, 120));
        out
    }

    /// The shared encode must be byte-identical to the historical
    /// per-sink encodes: `.baln` output vs `baln::write_record` per record,
    /// and the BAM byte stream vs `bamio::create_writer` +
    /// `bamio::write_record` per record (gunzipped, @PG-free by construction).
    #[test]
    fn record_sinks_bytes_match_per_sink_encodes() {
        let index = tiny_index();
        let records = synthetic_records();

        // --- .baln side: old per-record encode (fresh Vec per call) ---
        let mut old_baln = Vec::new();
        let names: Vec<String> = index
            .reference
            .contigs
            .iter()
            .map(|c| c.name.clone())
            .collect();
        crate::baln::write_header(&mut old_baln, &names).expect("old baln header");
        for r in &records {
            crate::baln::write_record(&mut old_baln, r).expect("old baln write");
        }

        let baln_cell = Arc::new(Mutex::new(Vec::new()));
        let mut sinks = RecordSinks::new(
            None,
            Some(Box::new(SharedSink(baln_cell.clone()))),
            &index,
            1,
        )
        .expect("sinks");
        for r in &records {
            sinks.write(r).expect("sink write");
        }
        sinks.finish().expect("finish");
        let new_baln = baln_cell.lock().expect("sink lock").clone();
        assert_eq!(new_baln, old_baln, ".baln byte stream must be unchanged");

        // --- BAM side: old thread-local-buffer path via bamio::write_record ---
        let header = crate::bam::build_header(&index);
        let old_bam_cell = Arc::new(Mutex::new(Vec::new()));
        {
            let mut w = esperanto_bamio::create_writer(
                Box::new(SharedSink(old_bam_cell.clone())),
                &header,
                1,
            )
            .expect("old writer");
            for r in &records {
                esperanto_bamio::write_record(&mut w, &header, r).expect("old bam write");
            }
            w.into_inner().finish().expect("old finish");
        }
        let new_bam_cell = Arc::new(Mutex::new(Vec::new()));
        {
            let mut sinks = RecordSinks::new(
                Some(Box::new(SharedSink(new_bam_cell.clone()))),
                None,
                &index,
                1,
            )
            .expect("sinks");
            for r in &records {
                sinks.write(r).expect("sink write");
            }
            sinks.finish().expect("finish");
        }
        let old_raw = old_bam_cell.lock().expect("sink lock").clone();
        let new_raw = new_bam_cell.lock().expect("sink lock").clone();
        let gunzip = |raw: &[u8]| -> Vec<u8> {
            let mut out = Vec::new();
            flate2::read::MultiGzDecoder::new(raw)
                .read_to_end(&mut out)
                .expect("gunzip bgzf");
            out
        };
        assert_eq!(
            gunzip(&new_raw),
            gunzip(&old_raw),
            "BAM record byte stream must be unchanged"
        );
    }
}
