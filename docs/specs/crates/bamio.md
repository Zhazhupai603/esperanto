# bamio — BAM Reading and Writing (SEQ Orientation Contract)

The single owning crate for record-level BAM reading and writing. Write side: record model + binary encoder + BGZF multithreaded writer + header construction; read side: sequential record view + orientation-restoration helper. Downstream consumers see only BAM and do not know which aligner produced it.

**Consumption boundary**: the feature semantics of pile / scan are pinned to rust-htslib / pysam (see their respective specs) and do not go through bamio; bamio's read side serves flow entry-point checks and record-level inspection only — no region fetch / pileup.

## SEQ Orientation Contract (hard contract)

SAM specification rule (both write and read sides are reference-forward):

- **Write side**: for an aligned minus-strand record, SEQ = `revcomp(original read)`, QUAL reversed in step, FLAG has 0x10 set; plus-strand and **unmapped records keep their original orientation** (unmapped FLAG 0x4, MAPQ 0).
- **Read side**: SEQ in a BAM is always stored in reference-forward orientation and consumers use it directly; when the original read orientation is needed (e.g., falling back to FASTQ), use `restore_original(flag, seq, qual)`: 0x10 set → revcomp + QUAL reversed again, otherwise unchanged.
- SE/PE consistent; PE is constructed by the caller per end under this contract (see the map spec).

`revcomp` semantics (bit-identical to the aligner's internal implementation): case-insensitive decoding, A↔T, C↔G, **all other characters (including IUPAC degenerate codes and lowercase) output as `N`**; output is all uppercase. Applying revcomp twice to the same input does not restore sequences containing non-ACGT characters (lossy, consistent with current behavior, frozen).

## Record Model

```rust
pub struct OutRecord {
    pub name: String,          // 1..=254 visible ASCII (not '@'), otherwise the slow fallback path reports an error
    pub flag: u16,
    pub mapq: u8,              // ≤60; unmapped records are encoded as 255 (BAM convention)
    pub aln: Option<AlnView>,  // None = unmapped
    pub seq: Vec<u8>,          // storage orientation (orientation contract already applied)
    pub qual: Vec<u8>,         // raw Phred [0,93]; empty = missing quality (encoded as 0xFF fill)
    pub mate: Option<(i32, i32, i32)>, // (mate contig id, mate pos0, tlen); contig<0 treated as absent
}

pub struct AlnView {
    pub contig: u32,           // header SQ order
    pub pos: u32,              // 0-based leftmost reference position
    pub cigar: Vec<CigarOp>,
    pub tags: Vec<RawTag>,     // ordered; encoded verbatim in the given order
}

pub enum CigarOp { Match(u32), Ins(u32), Del(u32), RefSkip(u32), SoftClip(u32) } // =/X folded into Match
pub struct RawTag(pub [u8; 2], pub TagValue);
pub enum TagValue { Char(u8), Int(i32), Str(String) }
```

- FLAG bit constants: `PAIRED 0x1 / PROPER_PAIR 0x2 / UNMAPPED 0x4 / MATE_UNMAPPED 0x8 / REVERSE 0x10 / MATE_REVERSE 0x20 / READ1 0x40 / READ2 0x80 / SECONDARY 0x100 / SUPPLEMENTARY 0x800`.
- **Fixed tag order**: map output is `XS → EA → EK → RE` (XS only when spliced junctions present, EA only when >0, RE only for rescued records); bamio does not sort or validate semantics — it encodes in the order given by `tags`; the fixed order is guaranteed by the caller's construction.

## Header and Writer

- `build_header(contigs: &[(name, len)], comment: &str) -> Header`: SQ lines in the given order; comment supplied by the caller (preserves each artifact's historical string, byte-identical).
- `create_writer(w, header, threads)`: BGZF multithreaded compression, worker count `clamp(threads, 1, 4)`; block-ordered merge, byte-deterministic output (same input + same worker count ⇒ same bytes).
- Record order = input order (sorting is left to downstream).

## Encoder

- Fast path `try_encode(buf, rec) -> Option<io::Result<()>>`: writes the BAM record binary directly (32-byte fixed header + read_name + cigar + 4-bit SEQ + QUAL + tags), byte-identical to noodles `write_alignment_record`; `None` = shape not supported (illegal name / cigar op count >65535 / qual >93 / contig > i32::MAX), and the caller falls back to the RecordBuf slow path.
- SEQ 4-bit packing: high nibble first, odd-length tail padded with 0 in the low nibble; the base-code table is case-insensitive, characters not in `=ACMGRSVTWYHKDBN` are treated as N (0xF).
- reg2bin per SAM §5.3 (start/end 0-based closed interval); span 0 and unmapped use bin 4680.
- qual same length as seq → written as-is; qual empty → 0xFF fill; any other length mismatch → `InvalidInput` error.

## Read Side

```rust
pub struct InRecord {
    pub name: String, pub flag: u16, pub mapq: u8,
    pub ref_id: i32, pub pos: i64,          // pos 0-based; -1 = no reference
    pub cigar: Vec<CigarOp>, pub seq: Vec<u8>, pub qual: Vec<u8>,
    pub tags: Vec<RawTag>,                   // file order
}
```

- `open_sequential(path) -> iterator of io::Result<InRecord>`: sequential read of the whole file; no index/region queries.
- 4-bit SEQ decoded back to ASCII (uppercase code table); QUAL 0xFF segments stay 0xFF (the caller decides on missing quality).
- `restore_original(flag, seq, qual) -> (Vec<u8>, Vec<u8>)` restores per the orientation contract.

## Dependencies

`noodles` (bam/sam/bgzf), `noodles-bgzf`, rust-htslib (only for the sort module's .bai index building), thiserror. No dependency on any other esperanto crate.

map's raw.bam is streaming (read-order) output, whereas pile/score use IndexedReader for region random access, which requires coordinate sorting + `.bai`. The sort module is that bridge. **The contract covers only map-produced BAM** (CIGAR subset M/I/D/N/S, tag types Char/Int/Str); generic sorting of arbitrary third-party BAMs is outside this contract.

- `coordinate_sort(input, output, opts) -> io::Result<SortStats>`; `SortOptions{max_in_memory_records: usize (default 2_000_000), temp_dir: Option<PathBuf>, threads}`; `SortStats{records, chunks}`.
- Algorithm: sequential read → split into chunks by record count → stable sort within each chunk by `(ref_id, pos)` (equal keys keep input order) → write temporary shards → k-way stable merge to output. Total-order sort key + stable merge ⇒ byte-identical output for the same input, independent of chunk size.
- Records are moved via noodles `RecordBuf` (original fields preserved, no re-encoding of semantics); the header is copied byte-for-byte from the input header (the SO line is not modified — map's historical header string stays byte-identical).
- Temporary shards live in `temp_dir` (default `<output>.sorttmp/`), named `chunk_%06d` deterministically; the whole directory is removed on success; residue after failure is allowed and overwritten on rerun.
- After writing, build the index: `rust-htslib bam::index::build` (BAI format), producing `<output>.bai`. Index contents are determined by htslib; byte-deterministic for the same input and same htslib version.
- Unmapped records (ref_id < 0) sort after all aligned records, preserving input order among themselves.

Self-check: after sorting an archived map artifact, two runs on the same input + different chunk caps (e.g., 1000 vs all-in-memory) produce byte-identical BAM; pysam read-back verifies monotonically non-decreasing coordinates + multiset equality of the record set; `.bai` is consumable by IndexedReader fetch.

## sort — Coordinate Sorting + Indexing (in-pipeline step of flow)

## Self-Check Points

- Minus-strand aligned record: after writing, pysam reads back SEQ == revcomp(original), QUAL reversed, FLAG 0x10; plus strand unchanged; unmapped unchanged + FLAG 0x4.
- revcomp: mixed-case + IUPAC input → all-uppercase output, non-ACGT becomes N.
- Fast path vs RecordBuf slow path: the same record produces identical bytes on both paths.
- Read-back loop: write → `open_sequential` read → field-by-field equality (tag order preserved).
- BGZF worker 1 vs 4: same input → identical output bytes.
- sort: same input, two runs + two chunk-cap settings → byte-identical BAM; monotonic coordinates; multiset equality of records; `.bai` fetchable.
