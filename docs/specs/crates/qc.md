# qc — RNA-seq FASTQ quality control

FASTQ cleaning with fastp semantics: adapter trimming, optional quality trimming, polyG trimming, whole-read filtering, paired-end synchronization, before/after statistics. Streaming pipeline; memory independent of input size; deterministic read output. Timing fields in reports vary by run.

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
2. **<1% triggers detection** using a k-mer seed and anchored extension:
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

- SE: `<stem>.clean.fq.gz`; PE: `<stem1>.clean_R1.fq.gz` and `<stem2>.clean_R2.fq.gz`; optional unpaired outputs use `<stem1>.unpaired_R1.fq.gz` and `<stem2>.unpaired_R2.fq.gz`. Each stem is the input filename with one `.gz` and one `.fastq` or `.fq` suffix removed.
- **Fqgz**: chunked processing; each chunk is an independent gzip member, concatenated in order (multi-member gzip is legal) → thread count does not affect the output bytes
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


## Boundary and case rules

- Adapter scans require at least five aligned bases. Insertion/deletion alternatives operate strictly inside the matching window.
- Paired overlap search uses positive offsets `0..len1-30` and negative offsets `1..len2-30`; reads of at most 30 bases do not enter this search.
- Empty read outputs contain no gzip members.
- Sequence comparisons and base statistics are case-insensitive. Reverse complementation preserves base case.
- Non-finite quality-fraction parameters are rejected.
- Each trimmed mate contributes independently to trimming statistics.
