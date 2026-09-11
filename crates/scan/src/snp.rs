//! Known-SNP position set from a VCF (plain or gzip).
//!
//! Only the contig and position columns are used; every data row is a
//! position to exclude. Contig names are matched with and without a `chr`
//! prefix so a `1`-style VCF works against a `chr1`-style BAM.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use crate::error::CallError;

pub struct SnpSet {
    /// contig -> sorted 0-based positions
    by_chrom: HashMap<String, Vec<u32>>,
}

impl SnpSet {
    pub fn load(path: &Path) -> Result<SnpSet, CallError> {
        let file = std::fs::File::open(path).map_err(|e| CallError::Io {
            path: path.display().to_string(),
            source: e,
        })?;
        let is_gz = path
            .extension()
            .is_some_and(|x| x == "gz" || x == "bgz");
        let reader: Box<dyn Read> = if is_gz {
            Box::new(flate2::read::MultiGzDecoder::new(file))
        } else {
            Box::new(file)
        };
        let mut by_chrom: HashMap<String, Vec<u32>> = HashMap::new();
        let mut rows = 0usize;
        for (ln, line) in BufReader::new(reader).lines().enumerate() {
            let line = line.map_err(|e| CallError::Io {
                path: path.display().to_string(),
                source: e,
            })?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut it = line.split('\t');
            let (Some(chrom), Some(pos)) = (it.next(), it.next()) else {
                return Err(CallError::Snp {
                    path: path.display().to_string(),
                    msg: format!("line {}: expected CHROM<TAB>POS", ln + 1),
                });
            };
            let pos1: u32 = pos.parse().map_err(|_| CallError::Snp {
                path: path.display().to_string(),
                msg: format!("line {}: bad POS '{pos}'", ln + 1),
            })?;
            if pos1 == 0 {
                continue;
            }
            by_chrom
                .entry(chrom.to_string())
                .or_default()
                .push(pos1 - 1);
            rows += 1;
        }
        if rows == 0 {
            return Err(CallError::Snp {
                path: path.display().to_string(),
                msg: "no data rows".into(),
            });
        }
        for v in by_chrom.values_mut() {
            v.sort_unstable();
            v.dedup();
        }
        Ok(SnpSet { by_chrom })
    }

    /// Positions stored for `chrom` (tries the name as given, then with the
    /// `chr` prefix added or removed).
    pub fn positions(&self, chrom: &str) -> Option<&[u32]> {
        if let Some(v) = self.by_chrom.get(chrom) {
            return Some(v);
        }
        let alt = match chrom.strip_prefix("chr") {
            Some(bare) => bare.to_string(),
            None => format!("chr{chrom}"),
        };
        self.by_chrom.get(&alt).map(|v| v.as_slice())
    }

    pub fn contains(&self, chrom: &str, pos0: u32) -> bool {
        self.positions(chrom)
            .is_some_and(|v| v.binary_search(&pos0).is_ok())
    }

    pub fn len(&self) -> usize {
        self.by_chrom.values().map(|v| v.len()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write_tmp(text: &str) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(format!(
            "esperanto-snp-{}-{}.vcf",
            std::process::id(),
            text.len()
        ));
        std::fs::write(&p, text).unwrap();
        p
    }

    #[test]
    fn loads_positions_and_matches_chr_prefix() {
        let p = write_tmp("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\n1\t100\trs1\nchr2\t5\trs2\n1\t100\trs3\n");
        let s = SnpSet::load(&p).unwrap();
        assert_eq!(s.len(), 2);
        assert!(s.contains("chr1", 99));
        assert!(s.contains("1", 99));
        assert!(s.contains("2", 4));
        assert!(!s.contains("chr1", 100));
        std::fs::remove_file(p).unwrap();
    }

    #[test]
    fn rejects_empty_vcf() {
        let p = write_tmp("#CHROM\tPOS\n");
        assert!(SnpSet::load(&p).is_err());
        std::fs::remove_file(p).unwrap();
    }
}
