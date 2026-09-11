# tidx — transcriptome k-mer index (.tidx)

Index crate of the L1 layer (transcript-first engine): builds a k-mer → (transcript, offset) index from a GENCODE GTF + reference FASTA, supporting fast queries.

## Public API contract

```rust
pub struct Tidx { /* private */ }
impl Tidx {
    pub fn open(path: &Path) -> io::Result<Self>;
    pub fn k(&self) -> u32;
    pub fn tx_count(&self) -> u32;
    /// Look up a canonical k-mer; returns &[(tx_id, offset)] sorted ascending by (tx_id, offset); empty slice on miss
    pub fn lookup(&self, canonical_kmer: u64) -> &[(u32, u32)];
    pub fn transcript_name(&self, tx_id: u32) -> &str;
    /// Verify that the reference BLAKE3 embedded in the index matches the given FASTA file; error on mismatch
    pub fn verify_reference(&self, fasta: &Path) -> Result<()>;
}
```

Build entry: `build(gtf_path, ref_path, out_path, &BuildOptions) -> Result<BuildStats>`.
`BuildOptions { k=31, prefix_bits=20, timestamp=now(), threads=16 }` (all overridable).
`BuildStats { tx_count, total_bp, total_entries, distinct_canonical, bucket_count, file_size, build_seconds, prefix_bits, ref_blake3 }`.

## k-mer encoding

- k=31 (default; A=00 C=01 G=10 T=11, either case accepted); the forward k-mer occupies the low 62 bits of a u64, top 2 bits always 0
- revcomp: take `3 - code` per 2-bit and reverse the order; canonical = min(forward, revcomp)
- Any non-ACGT byte (N/IUPAC) invalidates the window
- Sliding-window increment: `fwd' = ((fwd << 2) | entering) & ((1<<62)-1)`

## File format (all little-endian)

```
offset  size        content
0       128         Header
128     B           name_blob (transcript names, UTF-8, concatenated)
        (N+1)*4     name_offsets (u32; N=tx_count; name i = blob[off[i]..off[i+1]])
        (2^P+1)*8   bucket_offsets (u64; CSR exclusive suffix)
        E*w         keys (see version)
        E*8         payloads (u32 tx_id || u32 offset)
```

Header (128 bytes):

```
0   4   magic "TIDX"
4   4   version: 2 when 2k-P <= 32 (keys stored as u32, w=4), otherwise 1 (keys stored as u64, w=8)
8   4   k (valid range 15..=32)
12  4   m = 16 (reserved field)
16  4   prefix_bits P
20  4   reserved = 0
24  8   seed = 0x7469647873656564 (fixed constant)
32  32  ref_blake3 (BLAKE3 of the reference FASTA file bytes)
64  4   tx_count
68  4   reserved = 0
72  8   timestamp (unix seconds; build wall clock, can be pinned via --timestamp → byte-level reproducible for identical input)
80  8   total_entries E
88  8   bucket_count = 1<<P
96  8   name_blob_len
104 8   name_offsets_len = tx_count+1
112 16  reserved = 0
```

## Bucket structure

- The top P data bits of a canonical k-mer (bits [2k-P, 2k)) determine the bucket number; only the low 2k-P bits are stored within a bucket
- After the global records are sorted by (canonical, tx_id, offset), each bucket's records are naturally contiguous and ordered by the low bits within the bucket
- lookup = locate the bucket number → binary search the low bits within the bucket → return the payload slice
- version 2 (2k-P<=32) stores keys as u32; version 1 stores u64; the read side must support both versions

## Build pipeline (deterministic)

1. Parse the GTF (see below) → TranscriptSet; assign `tx_id` in lexicographic order of `transcript_id`
2. Group by contig; take the full sequence contig by contig in `.fai` order; within a contig, extract each transcript in parallel (rayon):
   - Concatenate exons (genomic ascending) → revcomp as a whole for the minus strand → obtain the 5'→3' transcript sequence
   - Slide a k=31 window; every N-free window yields (canonical, tx_id, offset), where offset = transcript coordinate
   - Transcript length < k, or any exon exceeding the contig → drop that transcript
3. Sort the global Vec by (canonical, tx_id, offset)
4. Write CSR directly; stream BLAKE3 over the reference file
5. Write the file (large BufWriter)

## GTF parsing semantics

- Only consume `transcript` and `exon` lines; coordinates 1-based inclusive → 0-based half-open `[start-1, end)`
- Aggregate by the `transcript_id` attribute; fall back to `transcript_id` when `transcript_name` is missing
- A transcript line refreshes name/contig/strand; an exon line only adds a span (must still be correct when exon lines appear before the transcript line)
- Transcripts without exons are dropped; each transcript's exons are sorted ascending by (start, end)
- Bad lines (fewer than 9 columns, missing transcript_id, strand not +-, start==0 or end<start) → error with line number
- Contig not present in the .fai → build error

## Error types (thiserror)

Io, BadMagic{magic,version}, UnsupportedVersion, Truncated (segment length mismatch), Inconsistent (k out of range, etc.), BadGtf{line,reason}, UnknownContig, RefMismatch (verify_reference).

## Self-checks (not written as tests; manual verification points after implementation)

- encode/revcomp/canonical are mutually consistent on 31-base "ACGT…" windows
- bucket + key_low can reconstruct the canonical k-mer
- Two builds with identical input and identical timestamp produce byte-identical artifacts

## Corrections from differential testing against the legacy binary

- Attribute parser aligned with the legacy version: values must be wrapped in double quotes (unquoted → the pair is skipped → missing transcript_id becomes an error); empty quoted values are legal; `=` separators (GFF style) are unsupported; semicolons inside quotes do not terminate the value.
- Intentional deviations (recorded): the read side is stricter than the legacy version (m==16 / seed constant / trailing bytes / monotonicity checks — the legacy version checks none of these); the k=32 mask is correct (the legacy version has a shift overflow); build-parameter validation (k∈15..=32, P bounds); unquoted attributes were rejected by the legacy version and, after alignment, are likewise rejected now.
