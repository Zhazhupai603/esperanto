//! Random-access FASTA reader driven by a samtools-style `.fai` index
//! (name, length, offset, line-bases, line-width per contig).
//!
//! Sequences are fetched as uppercased bytes; positions are 0-based,
//! half-open `[s, e)`, clamped to the contig length exactly like the
//! reference implementation.

use std::collections::HashMap;
use std::fs::File;
use std::io::{Read as _, Seek as _, SeekFrom};
use std::path::Path;

use anyhow::{anyhow, bail, Context};

/// One `.fai` row.
struct FaiEntry {
    /// Contig length in bases.
    len: i64,
    /// Byte offset of the first base.
    offset: u64,
    /// Bases per line.
    line_bases: i64,
    /// Bytes per line (bases plus the terminator).
    line_width: i64,
}

/// Open FASTA via its `.fai` sidecar (`<fasta>.fai`).
pub struct FastaIndex {
    file: File,
    entries: HashMap<String, FaiEntry>,
    order: Vec<(String, i64)>,
}

impl FastaIndex {
    /// Open `fasta` and parse `<fasta>.fai`.
    pub fn open(fasta: &Path) -> anyhow::Result<Self> {
        if !fasta.is_file() {
            bail!("fasta {} not found", fasta.display());
        }
        let fai_path = format!("{}.fai", fasta.display());
        let text = std::fs::read_to_string(&fai_path)
            .with_context(|| format!("reading fasta index {}", fai_path))?;
        let mut entries = HashMap::new();
        let mut order = Vec::new();
        for line in text.lines() {
            if line.is_empty() {
                continue;
            }
            let c: Vec<&str> = line.split('\t').collect();
            if c.len() < 5 {
                bail!("bad .fai line: {line}");
            }
            let entry = FaiEntry {
                len: c[1]
                    .parse()
                    .with_context(|| format!("bad .fai length '{}'", c[1]))?,
                offset: c[2]
                    .parse()
                    .with_context(|| format!("bad .fai offset '{}'", c[2]))?,
                line_bases: c[3]
                    .parse()
                    .with_context(|| format!("bad .fai line-bases '{}'", c[3]))?,
                line_width: c[4]
                    .parse()
                    .with_context(|| format!("bad .fai line-width '{}'", c[4]))?,
            };
            order.push((c[0].to_string(), entry.len));
            entries.insert(c[0].to_string(), entry);
        }
        Ok(Self {
            file: File::open(fasta).with_context(|| format!("opening {}", fasta.display()))?,
            entries,
            order,
        })
    }

    /// Contigs in `.fai` file order with their lengths.
    pub fn contigs(&self) -> &[(String, i64)] {
        &self.order
    }

    /// Fetch `[s, e)` (0-based, clamped to the contig) as uppercased bases.
    pub fn fetch(&mut self, chrom: &str, s: i64, e: i64) -> anyhow::Result<Vec<u8>> {
        let ent = self
            .entries
            .get(chrom)
            .ok_or_else(|| anyhow!("contig '{chrom}' missing from the fasta index"))?;
        let s = s.max(0);
        let e = e.min(ent.len);
        if e <= s {
            return Ok(Vec::new());
        }
        let want = (e - s) as usize;
        let lb = ent.line_bases.max(1) as usize;
        let lw = ent.line_width.max(lb as i64 + 1) as u64;
        let start = ent.offset + (s as u64 / lb as u64) * lw + (s as u64 % lb as u64);
        let mut raw = vec![0u8; want + want / lb + 2];
        self.file
            .seek(SeekFrom::Start(start))
            .with_context(|| format!("seeking {chrom}:{s}"))?;
        // Short reads (end of file) are fine - same as the reference,
        // which pads nothing and truncates to `want` after filtering.
        let n = self
            .file
            .read(&mut raw)
            .with_context(|| format!("reading {chrom}:{s}-{e}"))?;
        raw.truncate(n);
        Ok(raw
            .into_iter()
            .filter(|&b| b != b'\n' && b != b'\r')
            .map(|b| b.to_ascii_uppercase())
            .take(want)
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp_pair() -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "esperanto-report-test-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("t.fa");
        // 12 bases per line, 13 bytes per line (LF).
        std::fs::write(&fa, b"ACGTACGTACGT\nTTTTAAAACCCC\nGGGGTTTTAAAA\n").unwrap();
        std::fs::write(dir.join("t.fa.fai"), b"t\t36\t0\t12\t13\n").unwrap();
        dir
    }

    #[test]
    fn fetch_spans_line_breaks_and_clamps() {
        let dir = tmp_pair();
        let fa = dir.join("t.fa");
        let mut idx = FastaIndex::open(&fa).unwrap();
        assert_eq!(idx.fetch("t", 0, 4).unwrap(), b"ACGT");
        assert_eq!(idx.fetch("t", 10, 16).unwrap(), b"GTTTTT");
        assert_eq!(idx.fetch("t", 10, 17).unwrap(), b"GTTTTTA");
        assert_eq!(idx.fetch("t", -5, 3).unwrap(), b"ACG");
        assert_eq!(idx.fetch("t", 30, 999).unwrap(), b"TTAAAA");
        assert_eq!(idx.fetch("t", 40, 50).unwrap(), Vec::<u8>::new());
        assert_eq!(idx.contigs().len(), 1);
        assert_eq!(idx.contigs()[0].1, 36);
        assert!(idx.fetch("nope", 0, 4).is_err());
        std::fs::remove_dir_all(&dir).ok();
    }
}
