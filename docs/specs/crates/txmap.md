# txmap — transcript→genome projection + multi-isoform attribution

Coordinate-mapping crate of the L1 layer: projects transcript-coordinate intervals onto forward genomic coordinates (with CIGAR), and performs deterministic attribution among multiple transcript hits.

## Coordinate systems (strict)

- GTF input: 1-based inclusive
- All output genomic coordinates: **0-based half-open `[start, end)`, forward reference strand** (SAM convention)
- Transcript coordinate `tx_start`: 0-based, along the transcript direction 5'→3'; for minus-strand transcripts, tx=0 corresponds to transcription start

## Core types

```rust
pub enum Strand { Plus, Minus }
pub struct Exon { pub g_start: u32, pub g_end: u32 }        // 0-based half-open, forward genomic
pub struct TranscriptRecord { pub name: String, pub contig: String, pub strand: Strand, pub exons: Vec<Exon> }
pub enum CigarOp { Match(u32), RefSkip(u32) }               // extensible; two kinds for now
```

## project: transcript interval → genome

`TxMap::project(tx_id, tx_start, len) -> Option<(contig_id, genomic_start, Vec<CigarOp>)>`

- `len==0` or `tx_start >= tx_len` → `None`
- **Overhang clamping**: when `tx_start+len > tx_len`, truncate at the transcript end (does not return None); design rationale: reads with poly(A) tails / slight overextension are real
- Walk the exons and compute the intersection of `[tx_start, end)` with each exon's transcript interval:
  - plus strand: `(exon.g_start + r_lo, exon.g_start + r_hi)`
  - minus strand: reflect within the exon `(exon.g_end - r_hi, exon.g_end - r_lo)`
- All pieces are **sorted ascending by genomic start** (minus-strand pieces are produced in descending order; sorting restores forward order)
- CIGAR construction: insert `RefSkip(gap)` for the gap between adjacent pieces, one `Match(len)` per piece; `genomic_start` = the leftmost piece's start

## attribute: multi-isoform attribution

Input: a set of candidate projections for the same read, `AttributionCandidate { tx_id, contig, pos, cigar }`; output: a unique `Placement { contig, pos, cigar, mapq, tx_ids }` (empty input → `placement: None`).

Rules (strict order):
1. **Merge identical projections**: candidates with exactly identical (contig, pos, cigar) are merged into one group; members = the group's tx_id list
2. **Junction support count decides**: within a group's CIGAR, each `RefSkip`'s forward intron interval scores 1 point if it hits the known junction set (the full intron set of all GTF transcripts); the highest score wins
3. **Tie**: multiple groups tied at the highest score → `MAPQ = 0`; take the first in lexicographic order of (contig, pos, cigar_string) to guarantee reproducibility

MAPQ (when not tied): `w = |winner_members| / Σ|all_members|`; `MAPQ = min(60, floor(-10·log10(1-w)))`; `w>=1 → 60`, `w<=0 → 0`. `tx_ids` are sorted and deduplicated before output.

## JunctionSet

- The intron set of all transcripts: adjacent exons `(exon[i].g_end, exon[i+1].g_start)` (forward, 0-based half-open)
- Supports `contains(contig_id, start, end)` queries; sorted and deduplicated at build time

## .txmap file format (all little-endian; deterministic: all sets in sorted order)

```
magic        8 bytes "TXMAP001"
version      u32 = 1
source_hash  [u8;32]  BLAKE3 of the source GTF bytes
tx_count     u32
contig_count u32
  repeated: name_len u32, name UTF-8
junction_count u32
  repeated: contig_id u32, intron_start u32, intron_end u32
transcripts (in name lexicographic order = tx_id order):
  repeated: name_len u32, name; contig_id u32; strand u8 (0=+,1=-); exon_count u32;
        repeated exon: g_start u32, g_end u32
```

No timestamps, no padding, no unstable iteration order: two builds from the same input GTF are byte-identical.

## Other API

`from_records / from_gtf / open / save / tx_count / tx_id(name) / tx_name(id) / tx_len(id) / strand(id) / contig_of(id) / contigs() / transcripts() / junctions() / source_hash()`.
`TranscriptRecord::validate()`: adjacent exons (in transcript order) overlapping or out of order → Format error.

## Error types (thiserror)

Io, Magic{expected,found}, Version{file,code}, Format(String).

## GTF parsing

Same GENCODE semantics as the tidx crate (1-based→0-based half-open, aggregate by transcript_id, drop transcripts without exons, assign ids in transcript-name sort order). **Note: txmap's TranscriptRecord.exons are stored in transcript direction (5'→3')** — after reading, exons of minus-strand transcripts must be sorted in descending genomic coordinate order (ascending for the plus strand), so that `exons[0]` is always the transcription-start end; `introns()` takes adjacent-exon windows to obtain (exon[i].g_end, exon[i+1].g_start); the implementation must ensure the produced intervals are valid forward intervals (start<=end).

## Self-check points

- Single-exon plus strand: project(0, len) → single Match
- Two exons: Match + RefSkip + Match
- Minus-strand reflection: tx=0 maps to the end of the rightmost exon
- Attribution: identical projections merged; junction support count decides; ties → MAPQ=0

## Corrections from differential testing against the legacy binary

1. transcript_type whitelist: keep only `protein_coding` and `lncRNA` (fall back to gene_type when transcript_type is missing).
2. TranscriptRecord.name = transcript_id (not transcript_name).
3. Exons are stored in ascending exon_number order (= transcript order), deduplicated across rows with the same (start,end); when exon_number is non-monotonic, fall back to genomic sort (reversed for the minus strand).
4. The junction table stores raw transcript-order pairs (exon[i].g_end, exon[i+1].g_start) — the minus strand produces inverted (end,start) intervals, **stored as-is without normalization** (bug compatibility, part of the serialization format).
5. strand "." is accepted and treated as the plus strand.
6. Both transcript lines and exon lines are consumed; exon-only transcripts without a transcript line are dropped.

## Correction to fix #4 (science-first ruling)

Inverted junctions were originally ruled "keep for bug compatibility" — after scientific review this was changed to **normalized forward intervals**: under the legacy semantics, minus-strand annotated splice sites could never match in junction support, which is a genuine scientific defect; model features do not pass through txmap, so fixing this does not affect the gold standard. The read side stays compatible with old files (normalize on load). chr22 verification: the old and new sets are element-wise equivalent after normalization (6422/6422, 3424 inverted intervals normalized, 0 lost).
