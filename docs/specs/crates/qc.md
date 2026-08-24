# qc — RNA-seq FASTQ quality control

FASTQ cleaning with fastp semantics: adapter trimming, optional quality trimming, polyG trimming, whole-read filtering, paired-end synchronization, before/after statistics. Streaming pipeline; memory independent of input size; byte-level deterministic output.

## Parameters and defaults

```rust
pub struct QcParams {
    pub r1: Vec<PathBuf>,          // R1 inputs (plain/.gz; multiple lanes concatenated in order)
    pub r2: Vec<PathBuf>,          // empty = single-end; when non-empty, count must equal r1
    pub out_dir: PathBuf,
    pub adapter_trim: bool,        // default true
    pub adapters_r1: Vec<String>,  // empty = built-in Illumina table
    pub adapters_r2: Vec<String>,
    pub detect_adapter_se: bool,   // default false, SE only; see "SE adapter auto-detection"
    pub pe_overlap: bool,          // default true (PE only)
    pub qtrim: bool,               // default false
    pub qtrim_cutoff: u8,          // default 20
    pub trim_front1/trim_tail1/trim_front2/trim_tail2: usize, // default 0
    pub polyg: PolygMode,          // Auto|On|Off, default Auto
    pub min_len: usize,            // default 15
    pub n_max: usize,              // default 5
    pub q15_frac_max: f64,         // default 0.4
    pub keep_unpaired: bool,       // default false
    pub threads: usize,            // 0 = all cores
    pub out_format: OutFormat,     // Fqgz (default) | Bfq
}
```

Built-in adapter table: R1 = TruSeq R1 `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` + Nextera `CTGTCTCTTATACACATCT`; R2 = TruSeq R2 `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` + Nextera.

## Processing order (per read, strict)

1. **Fixed trimming**: `[front, len-tail)`
2. **qtrim (optional)**: BWA-style 3' quality trimming — accumulate a deficit from the 3' end, `sum += cutoff - (q-33)`, stopping when the deficit goes negative; cut at the position of the maximum accumulated value; an all-low-quality read may be cut to 0
3. **polyG (per mode)**: scan the G tail backward from 3'; non-G count budget = min(5, scanned/8), stop when the budget is exceeded; only cut when the tail length is ≥10
4. **Adapter trimming** (see below)
5. **Whole-read filtering** (see below)

## Adapter trimming semantics (aligned with fastp trimBySequence)

Single-adapter match `match_adapter(seq, adapter) -> Option<usize>` (returns the cut point; 0 = the whole read is adapter; None = no match):
- alen<4 or rlen<4 → None
- **Gapless scan**: start offset from -4 when alen≥16 (from -3 when alen≥12, from -2 when ≥8, otherwise 0), tolerating A-tailing adapter dimers; shift right one base at a time; comparison length `cmplen = min(rlen-pos, alen)`, mismatch budget `cmplen/8`; return the **first** hit position
- **Single-insertion-in-read variant**: budget `cmplen/8 - 1` (saturating); a skip k exists such that `a[..k]+a[k+1..]` matches the adapter's first cmplen bases
- **Single-deletion-in-read variant**: symmetric to the above (the adapter has one more base than the read)
- Multiple adapters = applied one by one in table order (truncate on hit, then proceed to the next)

Paired-end:
- When `pe_overlap=true`, **overlap analysis first**: rc(R2) against R1, scanning offsets 0→+ then 0→- (read-through); the first offset with `mm <= min(5, ol/10)` wins; minimum overlap 30nt; **only when offset<0**, trim both ends to overlap_len
- When no reliable overlap is found: fall back to known-table matching **only if the user explicitly provided an adapter table** (fastp PE does not auto-scan known sequences by default, to prevent random mis-trimming); for SE or pe_overlap=false, use the known table directly (empty = built-in table)

## SE adapter auto-detection (optional, off by default)

Participates only for SE with adapter_trim=true and no explicit user-provided table. `detect_adapter_se` defaults to false; when false, all outputs (including qc.json) are byte-identical to the disabled case; the `adapter_source`/`detected_adapter_se` fields are written to qc.json only when enabled (when enabled but the table path is taken, the two fields are `"table"`/null).

Prescan: buffer the first N=5000 reads before streaming starts (bounded memory); after the decision, the buffered reads go directly into the main stream without re-reading the input.

Decision order:

1. **Table hit rate**: if the built-in table's hit fraction on the 5000 reads is ≥1% → use the table (`adapter_source="table"`, do not output `detected_adapter_se`); no detection.
2. **<1% triggers detection** (k-mer seed + anchored extension algorithm; scientific correction: the cross-read shared invariant under read-through is the adapter **5' prefix** — internal substrings drift with insert length — while the 3' suffix drifts with the insert and yields no stable shared suffix, so no suffix tree is used):
   a. Count all 12-mers within each read's 3' tail window (the last min(rlen,36) bases) (2-bit encoding; skip 12-mers containing non-ACGT);
   b. The most frequent 12-mer becomes the seed (ties broken by the smallest 2-bit encoding, for determinism); count < max(20, N/200) → no candidate, `adapter_source="none"`, no trimming;
   c. Anchored extension: align each read at the seed's leftmost occurrence within its tail window; extend base by base toward 5' and 3' — continue at a position while the majority base count there is ≥ max(20, N/200) and covers ≥60% of the reads participating at that position (ties broken A<C<G<T<N), otherwise stop (the 60% relative threshold prevents over-extension into the insert's genomic region under high support); sequences diverge across reads at the insert boundary, so 5'-ward extension terminates naturally;
   d. Extended sequence total length ≥10 → candidate adapter; the candidate is re-validated with `match_adapter` over the prescan buffer: it holds only with ≥ max(20, N/200) hits, giving `adapter_source="detected"`, and the candidate enters main-stream matching (replacing the table); otherwise → `adapter_source="none"`, no trimming.
3. No randomness and no floating-point ordering dependencies anywhere; identical input always yields identical results. Case rules match the rest of the library (eq_ignore_ascii_case; candidates are output in uppercase).

## Whole-read filtering (decision order: length → N → low quality; only the first failure reason is counted)

1. `len < min_len`: empty and adapter/overlap trimming occurred → `AdapterOnly`, otherwise `TooShort`
2. N count > n_max → `TooManyN`
3. Bases with Q<15 (byte < 33+15=48) > len × q15_frac_max → `LowQuality` (strictly greater; exactly equal is not dropped)

## Paired-end rules

- Mate-name check: R1/R2 name prefixes (before the first space) must match; 1/2 suffixes are allowed
- Either end failing → drop the whole pair (each end counts its own failure reason); with `keep_unpaired=true`, the passing end is written to the unpaired output
- Phred prescan: inspect the first 10000 reads; non-+33 encoding (a +64 signature observed) → error out and refuse
- polyG Auto: read the instrument name of R1's first read; if it indicates the NextSeq/NovaSeq family (two-color instruments) → enable

## Outputs

- SE: `<stem>.clean.fq.gz` (or `.bfq`); PE: plus R2 and optional unpaired1/2; stem = R1 file name stripped of `.gz`/`.fastq`/`.fq`
- **Fqgz**: chunked processing; each chunk is an independent gzip member, concatenated in order (multi-member gzip is legal) → thread count does not affect the output bytes
- **Bfq binary format**: 16-byte header (`EBFQ` + version u8=1 + 3 reserved + read_count u64); per read: name_len u16 + name + seq_len u16 + seq + qual; all little-endian, uncompressed; byte-identical to decompressed fq after conversion to text
- `qc.json` + `qc.html` reports

## Statistics (all-integer accumulation; converted to ratios once at report time)

One set each for before/after: reads, bases, Q20, Q30, GC, N, per-cycle quality sums, per-cycle coverage counts, per-cycle base counts [A,C,G,T,N].
Counts: adapter_reads/bases, polyg_reads/bases, qtrim_bases, each of the four failure reasons, unpaired_written.

## qc.json fields

`esperanto_qc_version, params{adapter_trim,pe_overlap,qtrim,polyg,min_len,n_max, detect_adapter_se?}, input{…,instrument_polyg}, summary{reads_before/after,bases_before/after,q20_rate,q30_rate,gc_rate,n_rate,duplication_estimate?}, filter_reasons{low_quality,too_many_n,too_short,adapter_only}, trimming{adapter_reads,adapter_bases,polyg_reads,polyg_bases,qtrim_bases}, per_cycle{position,mean_quality,base_frac{a,c,g,t,n}}, per_cycle_before, summary_before, elapsed_seconds, adapter_source?, detected_adapter_se?` (? = output only when detect_adapter_se=true)

## Pipeline structure

The reader thread chunks input into fixed blocks (2048 reads/block) → rayon workers parse + process + compress → the writer thread writes out in order. Chunked parsing happens on the worker side (raw byte blocks are transferred); PE's two blocks are zipped in order. Errors are passed through a shared slot; the writer skips on the sentinel.

## Self-check points

- polyG: 16nt non-G + 25G → cut to ≤16; a 5G tail is not cut; G10+A+G14 → cut to 0
- qtrim: 80×Q38+20×Q2, cutoff 20 → cut to 80
- adapter: TruSeq R1 full length appended to a read tail → cut at the adapter start; whole read = adapter → Some(0)
- filtering: exactly 40% low quality is not dropped, >40% is dropped
- overlap: 60nt insert + 93 read length read-through → offset<0, both ends cut to 60

## Corrections from differential testing against the legacy binary

1. Gapless scan upper bound: `pos < rlen-4` (shortest alignment 5 bases), not `<= rlen-4` (4-mer random-tail false positives).
2. Gate before the indel variants: only attempt k-skip when "the first min(8,cmplen) bases gapless-align within the variant budget" (the legacy v12 fast path is intentional semantics, not a pure optimization).
3. Indel k range: `1..cmplen` (the insertion/deletion must fall strictly inside the window; including either end equals shifting the window).
4. Indel scan start: pos=0 (no negative starts); insertion-variant position bound `pos < rlen-5`, deletion `pos < rlen-4`.
5. qc.json numbers are output at full precision (the legacy does no rounding); instrument_polyg is a bool; summary_before contains ratios only.
6. Two-color instrument detection (aligned verbatim with the legacy): token = the header's first field before the colon; enabled when any of: `A0<digit>` prefix (NovaSeq), `NS5`/`NB5` prefix (NextSeq), or the **entire header** case-sensitively containing the literal `NextSeq`/`NovaSeq`. Fuzzy substring matching (NS*/NB*) is not used.
7. PE overlap scan bounds (differential-testing correction): positive offsets `0..l1-30` (open interval), negative offsets `1..l2-30` (open interval); boundary positions where ol is exactly 30 are **not tried** (30-mer false overlaps from periodic adapter tails land exactly there). If either read is ≤30, directly no overlap.
8. bfq legacy bug compatibility (differential-testing correction): (a) read_count is backfilled only for R1/SE files, with value = total kept reads of both ends (for PE = 2× the pair count); R2/unpaired files always have count 0. (b) unpaired file names are hardcoded `.fq.gz`, but in bfq mode their content is still EBFQ (header + raw records).
9. Zero-record output files = 0 bytes (legacy semantics); no empty gzip member is written.
10. PE overlap trimming is counted per end independently: each end increments adapter_reads += 1 only when it was itself truly trimmed (bases>0); when one end is already shorter than ol, that end is not counted.

## Corrections from differential testing against the legacy binary (round 3)

11. PE output naming aligned with the legacy: `<stem1>.clean_R1.<ext>` + `<stem2>.clean_R2.<ext>` (R2 uses its own file's stem); unpaired = `<stem1>.unpaired_R1.fq.gz` / `<stem2>.unpaired_R2.fq.gz`; SE remains `<stem>.clean.<ext>`.
12. params gains the `q15_frac_max` field.
13. revcomp complements only uppercase ACGT; all other bytes pass through unchanged (including lowercase).
14. stem stripping: strip one .gz, then one .fastq/.fq (stop at the first hit); no repeated stripping.
15. Intentional deviations (recorded): the gzip compression backend and level differ (zlib-6 vs libdeflate-1), with decompressed content byte-identical; bfq fields >65535 error out (the legacy silently truncates, producing corrupt records); ≥64nt custom adapters panic in the legacy version but work in the new one; q15=NaN is rejected (the legacy silently disables low-quality filtering); run returns () (reporting goes through qc.json); html is a minimal rendering; the legacy bfq.rs conversion tool was not ported.
16. **Case = scientific semantics (ruling)**: sequence letter case carries no information — lowercase acgt are the same bases as uppercase. Statistic buckets are case-insensitive (gc counts C/c/G/g; the N bucket takes only N/n plus other codes); polyG continues the tail on both g/G; adapter matching, PE overlap, and the indel gate all use `eq_ignore_ascii_case`; revcomp complements case-insensitively and preserves case; the N filter counts both N and n. **Intentional deviation from the legacy** (the legacy is byte-sensitive throughout: lowercase goes into the N bucket, lowercase g is not tail-trimmed) — rationale: qc output feeds our own aligner rather than a frozen model, so no model contract constrains it; standard Illumina data is all-uppercase, where behavior is byte-identical to the legacy (regression guarantee). Verification = case-equivariance invariant: after randomly changing the case of any data, the output clean sequences (ignoring case), qualities, kept set, and qc.json numbers must be exactly identical (20 datasets × 40% lowercased, 0 failures).
