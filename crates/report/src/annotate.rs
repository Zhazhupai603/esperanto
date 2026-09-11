//! GTF annotation (gene rows + transcript CDS rows) and the A-to-I
//! recoding (amino-acid change) check.
//!
//! Semantics are ported verbatim from the validated reference
//! implementation:
//!
//! * CDS rows per transcript are ALWAYS sorted ascending by coordinate and
//!   concatenated in that order; minus-strand transcripts are then built by
//!   reverse-complementing the whole concatenation (never by concatenating
//!   in descending order).
//! * The reading frame is the phase of the first CDS row in transcription
//!   order: `cds[0]` for plus strand, `cds[last]` (largest coordinate) for
//!   minus strand.
//! * A site is recoding iff it lies inside a CDS, the transcript-frame base
//!   is `A`, the A->G substitution changes the amino acid, and the site
//!   probability is >= 0.2 (checked by the caller).

use std::collections::HashMap;
use std::path::Path;

use anyhow::Context;

use crate::fasta::FastaIndex;

/// Standard genetic code over TCAG base order (`*` = stop).
const AA1: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRLFIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

fn base_idx(b: u8) -> Option<usize> {
    match b {
        b'T' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// Amino acid of a 3-base codon; `None` mirrors the reference `"?"` (any
/// non-ACGT base makes the codon untranslatable).
fn codon_aa(codon: [u8; 3]) -> Option<char> {
    let i = base_idx(codon[0])? * 16 + base_idx(codon[1])? * 4 + base_idx(codon[2])?;
    Some(AA1[i] as char)
}

/// Gene row: 0-based start, end-exclusive end, gene name.
struct GeneRow {
    start: i64,
    end: i64,
    name: String,
}

/// Transcript CDS rows: (0-based start, end-exclusive end, GTF phase).
struct Tx {
    chrom: String,
    strand: String,
    gene: String,
    cds: Vec<(i64, i64, i32)>,
}

/// GTF genes + transcripts with lookup indexes.
pub struct Annotation {
    /// chrom -> gene rows sorted by (start, end, name).
    genes: HashMap<String, Vec<GeneRow>>,
    /// transcript_id -> transcript.
    txs: HashMap<String, Tx>,
    /// gene name -> (chrom, start, end); the last GTF row wins.
    gcoord: HashMap<String, (String, i64, i64)>,
    /// chrom -> (CDS start, CDS end, transcript_id) sorted.
    tx_by_chrom: HashMap<String, Vec<(i64, i64, String)>>,
}

/// Value of the first double-quoted string in `key ...` attribute parts.
fn quoted_attr(attrs: &str, key: &str) -> Option<String> {
    for part in attrs.split(';') {
        let part = part.trim();
        if let Some(rest) = part.strip_prefix(key) {
            let q1 = rest.find('"')?;
            let after = &rest[q1 + 1..];
            let q2 = after.find('"')?;
            return Some(after[..q2].to_string());
        }
    }
    None
}

impl Annotation {
    /// Parse gene and CDS rows from a GTF file.
    pub fn load(gtf: &Path) -> anyhow::Result<Self> {
        let text = std::fs::read_to_string(gtf)
            .with_context(|| format!("reading GTF {}", gtf.display()))?;
        let mut genes: HashMap<String, Vec<GeneRow>> = HashMap::new();
        let mut txs: HashMap<String, Tx> = HashMap::new();
        for (ln, line) in text.lines().enumerate() {
            if line.starts_with('#') {
                continue;
            }
            let c: Vec<&str> = line.split('\t').collect();
            if c.len() < 9 {
                continue;
            }
            match c[2] {
                "gene" => {
                    if let Some(name) = quoted_attr(c[8], "gene_name") {
                        genes.entry(c[0].to_string()).or_default().push(GeneRow {
                            start: c[3].parse::<i64>().with_context(|| {
                                format!("GTF line {}: bad gene start '{}'", ln + 1, c[3])
                            })? - 1,
                            end: c[4].parse::<i64>().with_context(|| {
                                format!("GTF line {}: bad gene end '{}'", ln + 1, c[4])
                            })?,
                            name,
                        });
                    }
                }
                "CDS" => {
                    let tid = match quoted_attr(c[8], "transcript_id") {
                        Some(t) => t,
                        None => continue,
                    };
                    let gene = quoted_attr(c[8], "gene_name").unwrap_or_default();
                    let start = c[3].parse::<i64>().with_context(|| {
                        format!("GTF line {}: bad CDS start '{}'", ln + 1, c[3])
                    })? - 1;
                    let end = c[4]
                        .parse::<i64>()
                        .with_context(|| format!("GTF line {}: bad CDS end '{}'", ln + 1, c[4]))?;
                    let phase = c[7].parse::<i32>().with_context(|| {
                        format!("GTF line {}: bad CDS phase '{}'", ln + 1, c[7])
                    })?;
                    let t = txs.entry(tid).or_insert_with(|| Tx {
                        chrom: c[0].to_string(),
                        strand: c[6].to_string(),
                        gene,
                        cds: Vec::new(),
                    });
                    t.cds.push((start, end, phase));
                }
                _ => {}
            }
        }
        for rows in genes.values_mut() {
            rows.sort_by(|a, b| (&a.start, &a.end, &a.name).cmp(&(&b.start, &b.end, &b.name)));
        }
        for t in txs.values_mut() {
            t.cds.sort();
        }
        let mut gcoord = HashMap::new();
        for (chrom, rows) in &genes {
            for r in rows {
                gcoord.insert(r.name.clone(), (chrom.clone(), r.start, r.end));
            }
        }
        let mut tx_by_chrom: HashMap<String, Vec<(i64, i64, String)>> = HashMap::new();
        for (tid, t) in &txs {
            tx_by_chrom.entry(t.chrom.clone()).or_default().push((
                t.cds[0].0,
                t.cds[t.cds.len() - 1].1,
                tid.clone(),
            ));
        }
        for v in tx_by_chrom.values_mut() {
            v.sort();
        }
        Ok(Self {
            genes,
            txs,
            gcoord,
            tx_by_chrom,
        })
    }

    /// Gene name covering `pos` (0-based expected by the caller contract).
    /// Candidates are the gene rows whose start is <= pos, walking back at
    /// most 20 starts from the insertion point; the smallest-start
    /// overlapping row wins (reference semantics).
    pub fn gene_at(&self, chrom: &str, pos: i64) -> Option<&str> {
        let lst = self.genes.get(chrom)?;
        // bisect_right over the sorted starts.
        let i = lst.partition_point(|g| g.start <= pos) as i64 - 1;
        if i < 0 {
            return None;
        }
        for j in (i - 20).max(0)..=i {
            let g = &lst[j as usize];
            if g.start <= pos && pos < g.end {
                return Some(&g.name);
            }
        }
        None
    }

    /// Genome coordinates for a gene name (the last GTF row wins).
    pub fn gene_coord(&self, name: &str) -> Option<(&str, i64, i64)> {
        self.gcoord.get(name).map(|(c, s, e)| (c.as_str(), *s, *e))
    }

    /// Transcripts of a chromosome sorted by (CDS start, CDS end, id) as
    /// `(start, end, transcript_id)` triples.
    pub fn transcripts(&self, chrom: &str) -> &[(i64, i64, String)] {
        self.tx_by_chrom
            .get(chrom)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    }

/// Lazily-built transcript CDS sequences + the recoding check.
pub struct Recoder<'a> {
    fa: &'a mut FastaIndex,
    txs: &'a HashMap<String, Tx>,
    /// transcript_id -> (CDS sequence in transcription order, frame).
    cache: HashMap<String, (Vec<u8>, i32)>,
}

impl<'a> Recoder<'a> {
    pub fn new(fa: &'a mut FastaIndex, ann: &'a Annotation) -> Self {
        Self {
            fa,
            txs: &ann.txs,
            cache: HashMap::new(),
        }
    }

    /// Build the transcript CDS sequence (ascending concatenation, whole
    /// reverse-complement for minus strand) and the transcription-order
    /// frame; cached per transcript.
    fn cds_seq(&mut self, tid: &str) -> anyhow::Result<(&[u8], i32)> {
        if !self.cache.contains_key(tid) {
            let t = &self.txs[tid];
            let mut seq: Vec<u8> = Vec::new();
            for &(s, e, _) in &t.cds {
                seq.extend_from_slice(&self.fa.fetch(&t.chrom, s, e)?);
            }
            let (frame, minus) = if t.strand == "-" {
                (t.cds[t.cds.len() - 1].2, true)
            } else {
                (t.cds[0].2, false)
            };
            if minus {
                revcomp(&mut seq);
            }
            self.cache.insert(tid.to_string(), (seq, frame));
        }
        let (seq, frame) = self.cache.get(tid).expect("just inserted");
        Ok((seq, *frame))
    }

    /// Gene name of a transcript.
    pub fn gene_of(&self, tid: &str) -> String {
        self.txs[tid].gene.clone()
    }

    /// `(ref_aa, alt_aa, aa_pos)` when `pos0` (as passed by the caller,
    /// the raw 1-based VCF position) is a non-synonymous A-to-G site of
    /// transcript `tid`.
    pub fn recode_at(&mut self, tid: &str, pos0: i64) -> anyhow::Result<Option<(char, char, i64)>> {
        let t = &self.txs[tid];
        let plus = t.strand != "-";
        let n = t.cds.len();
        let mut acc: i64 = 0;
        // CDS rows in transcription order (plus: ascending, minus: descending).
        for i in 0..n {
            let (s, e, _) = t.cds[if plus { i } else { n - 1 - i }];
            if s <= pos0 && pos0 < e {
                let off_in_row = if plus { pos0 - s } else { e - 1 - pos0 };
                let cds_off = acc + off_in_row;
                let (seq, frame) = self.cds_seq(tid)?;
                let ci = cds_off - frame as i64;
                if ci < 0 {
                    return Ok(None);
                }
                let codon_i = ci / 3;
                let ph = (ci % 3) as usize;
                let cstart = frame as i64 + codon_i * 3;
                if cstart < 0 || cstart as usize + 3 > seq.len() {
                    return Ok(None);
                }
                let codon = [
                    seq[cstart as usize],
                    seq[cstart as usize + 1],
                    seq[cstart as usize + 2],
                ];
                if codon[ph] != b'A' {
                    return Ok(None);
                }
                let mut alt = codon;
                alt[ph] = b'G';
                return match (codon_aa(codon), codon_aa(alt)) {
                    (Some(ra), Some(aa)) if ra != aa => Ok(Some((ra, aa, codon_i + 1))),
                    _ => Ok(None),
                };
            } else {
                acc += e - s;
            }
        }
        Ok(None)
    }
}

/// Reverse-complement ACGT in place (other bases pass through).
fn revcomp(seq: &mut [u8]) {
    seq.reverse();
    for b in seq.iter_mut() {
        *b = match *b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codon_table_basics() {
        assert_eq!(codon_aa(*b"ATG"), Some('M'));
        assert_eq!(codon_aa(*b"TGG"), Some('W'));
        assert_eq!(codon_aa(*b"TAA"), Some('*'));
        assert_eq!(codon_aa(*b"GGC"), Some('G'));
        assert_eq!(codon_aa(*b"NTG"), None);
        // A->G in codon position 2 of CAT (H) -> CGT (R): recoding.
        assert_ne!(codon_aa(*b"CAT"), codon_aa(*b"CGT"));
        // Synonymous: AGC -> GGC is not an A-to-G single-base edit; a true
        // synonymous A->G case: TTA (L) -> TTG (L).
        assert_eq!(codon_aa(*b"TTA"), codon_aa(*b"TTG"));
    }

    #[test]
    fn attr_parsing() {
        let a = "gene_id \"ENSG1\"; gene_name \"RPS24\"; transcript_id \"ENST1\";";
        assert_eq!(quoted_attr(a, "gene_name").as_deref(), Some("RPS24"));
        assert_eq!(quoted_attr(a, "transcript_id").as_deref(), Some("ENST1"));
        assert_eq!(quoted_attr(a, "gene_type"), None);
    }

    #[test]
    fn gene_at_walks_back() {
        // Hand-built annotation via a small GTF.
        let dir = std::env::temp_dir().join(format!("esperanto-geneat-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let gtf = dir.join("t.gtf");
        std::fs::write(
            &gtf,
            "chr1\tsrc\tgene\t101\t200\t.\t+\t.\tgene_id \"g1\"; gene_name \"FIRST\";\n\
             chr1\tsrc\tgene\t101\t400\t.\t+\t.\tgene_id \"g2\"; gene_name \"NESTED\";\n\
             chr1\tsrc\tgene\t1001\t1100\t.\t+\t.\tgene_id \"g3\"; gene_name \"AFTER\";\n",
        )
        .unwrap();
        let ann = Annotation::load(&gtf).unwrap();
        // pos 150 is inside both FIRST (100..200) and NESTED (100..400);
        // walking up from the newest start, but returning the smallest
        // start that overlaps => FIRST.
        assert_eq!(ann.gene_at("chr1", 150), Some("FIRST"));
        assert_eq!(ann.gene_at("chr1", 250), Some("NESTED"));
        assert_eq!(ann.gene_at("chr1", 90), None);
        assert_eq!(ann.gene_at("chr1", 1050), Some("AFTER"));
        assert_eq!(ann.gene_coord("AFTER"), Some(("chr1", 1000, 1100)));
        std::fs::remove_dir_all(&dir).ok();
    }
}
