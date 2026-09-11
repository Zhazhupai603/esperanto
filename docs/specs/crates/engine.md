# engine — L1 Transcript-First Matching Engine

The L1 engine of ESPERANTO v13: read a read → locate it in the transcriptome k-mer index → EA exact extension → (when necessary) EA-Myers verification → project to the genome → deterministic assignment. Purely functional, per-read independent, no global state of any kind; the parallel driver and output belong to flow/cli (Wave 4).

## Hard-Won Lessons (mandatory constraints at the head of this crate's spec)

1. **Validate direction-related conventions on both the write and read sides** (legacy write-side defects).
2. **When two references diverge, follow the one that generates the training data** (the htslib vs pysam fork precedent).
3. **A performance fast path is semantics**: any optimization with behavioral consequences (covered-skip, prefix gating) must be written into the spec.
4. **Spec errors propagate** (the N-base depth lesson): the spec is the source of truth; inconsistencies discovered during implementation must be written back into it.

## Public API Contract

```rust
pub enum CigarOp { Match(u32), Ins(u32), Del(u32), RefSkip(u32), SoftClip(u32) }  // SAM semantics
pub enum Strand { Plus, Minus }                    // sam_flag(): 0/16; letter(): "+"/"-"
pub enum L1Outcome {
    Aligned { contig: u32, pos: u32, strand: Strand, cigar: Vec<CigarOp>, score: i32, mapq: u8 },
    Fallback,
}
pub enum Branch { Full, Interrupted, GateFail, NoHit, TooShort }
pub struct ReadStats { pub queries: u32, pub extension_bases: u32, pub branch: Branch, pub extension_tx_count: u32 }

pub struct EngineConfig {
    pub flank: usize,            // 10 — Myers window flanks
    pub coverage_gate: f64,      // 0.90
    pub max_tx_candidates: usize,// 8
    pub max_diagonals_per_tx: usize, // 4
    pub prefetch_batch: usize,   // 16 (pure performance hint, no semantic effect)
    pub max_raw_hits: usize,     // 16384
    pub max_hits_per_kmer: usize,// 100
}

pub fn align_read(read, tidx: &impl Tidx, txmap: &impl TxMap, txseqs: &impl TxSeqs,
                  cfg: &EngineConfig, repeats: &impl RepeatTrack, stats: &mut ReadStats) -> L1Outcome;
```

Trait contracts:

```rust
pub trait Tidx { fn k(&self)->u32; fn tx_count(&self)->u32;
    fn lookup(&self, canonical_kmer:u64)->&[(u32,u32)];   // (tx_id, offset) ascending
    fn transcript_name(&self, tx_id:u32)->&str;
    fn prefetch(&self, _:u64) {} }                        // pure performance hint, default no-op
pub trait TxMap { fn project(&self, tx_id:u32, tx_start:u32, len:u32)->Option<(u32,u32,Vec<CigarOp>)>;
    fn tx_len(&self, tx_id:u32)->Option<u32>; fn strand(&self, _tx_id:u32)->Option<Strand> { None } }
pub trait TxSeqs { fn seq(&self, tx_id:u32)->&[u8]; }    // transcript forward sequence (minus strand already RC'd)
pub trait RepeatTrack { fn overlaps(&self, contig:u32, pos:u32, len:u32)->bool; }
pub struct NoRepeats;  // always false (default / testing / repeat rejection off)
```

`RepeatBed::load(path, contig_name_to_id)`: loads a BED (plain text or .gz) by its first 3 columns,
stores ascending intervals per contig plus a parallel start array; `overlaps(c,pos,len)`: binary-search
for the rightmost position with `start < pos+len`, check up to 3 entries backward (RepeatMasker
intervals may overlap) for any `end > pos`; BED rows whose contig is not in the mapping table are
silently skipped.

## k-mer Stream and Encoding

- 2-bit: A=0 C=1 G=2 T=3; first base in the high bits; canonical = min(fwd, revcomp); revcomp complements each 2-bit digit as 3-c and reverses the order
- Windows containing any non-ACGT byte (including N) are not emitted; re-seed after an N (the first window after an N is computed from N+1 onward, leaving no phantom codes)
- **Case-insensitive** (scientific decision): lowercase acgt are treated as valid bases and participate in encoding/extension/verification; revcomp preserves case-sign complementation (a↔t). The legacy version accepted uppercase only — behavior is bit-identical on standard all-uppercase input
- k is determined by `tidx.k()` (production = 31)

## Interleaved Query + Extension (Full-Information Query — No Covered-Skip)

Process k-mers in read-position order; extend immediately on each hit. **1.0.0 decision (differential-testing evidence, see SCIENCE-DEVIATIONS): covered-skip removed** — the legacy default "skip any k-mer whose start falls in a covered region" skipped boundary k-mers carrying true-locus/paralog evidence, producing false placements (salvaged only by MAPQ=0, losing recall) or overconfidence (missed paralogs, MAPQ wrongly set to 60). Full-information query strictly dominates on both failure modes; behavior is byte-identical to the legacy acceptance run with `L1_SKIP_COVERED=0` (chr22, 4002 reads, differential test passed).

For each k-mer: `tidx.lookup(canonical)`; hit count > `max_hits_per_kmer` → skip this k-mer.
For each hit `(tx_id, tx_off)`:
- Global raw_hits ≥ `max_raw_hits` → stop
- This tx already has ≥ `max_tx_candidates * max_diagonals_per_tx` extensions → skip
- **Strand determination**: compare the read's forward code at pos with the transcript's forward code at tx_off; equal → Plus (oriented_pos = pos), otherwise Minus (oriented_pos = read_len - k - pos)
- diagonal = tx_off - oriented_pos; deduplicate on (tx_id, strand, diagonal); per-tx diagonal count ≥ max_diagonals_per_tx → skip
- `extend_ea(oriented_read, oriented_pos, k, tx_seq, tx_off)` (Minus uses the RC read)
- The extension's covered interval (mapped back to original read coordinates) enters covered; the extension enters the list

## EA Exact Extension

`extend_ea(read, a, k, tx_seq, t) -> Extension {read_lo,read_hi,tx_lo,tx_hi,full}`:
the seed window is treated as EA-equal (guaranteed by the index hit); extend base by base on both sides of the seed, stopping a direction at the first true mismatch (EA predicate false). `full = (read_hi - read_lo == read_len)`. `read_cov = hi-lo`.
diagonal = tx_lo - read_lo.

## EA Predicate (Directional ADAR Rule)

Zero cost if and only if: (ref=A and read=G) or (ref=T and read=C) or (ref==read and both are ACGT).
Reverse conversions (G→A, C→T) are not free; N/lowercase/IUPAC/other bytes always mismatch.
(Note: lowercase is treated as the same uppercase base per the case-insensitivity decision — convert to uppercase before testing.)

## EA-Myers Bit-Parallel Verifier (bit-identical to the legacy implementation)

`myers::infix(read, text)` (m≤128, single u128) and `myers::long::infix` (128<m≤256, two blocks,
carry passed across blocks via hout). infix semantics: the full length of the pattern (read) must match, while text overhangs on both sides are free — block 0 takes hin=0 at every column (text left edge free) and tracks the running column-score minimum (text right edge free).
The core `calculate_block` recurrence and the EA-redefined `Peq` table must be bit-equivalent to the legacy implementation (this implementation was validated with 125K fuzz pairs at zero error). A global variant is also provided (m≤128 single word, else two blocks).

## Branches and Assignment (strict order)

1. **read length < k** → `Branch::TooShort` + Fallback.
2. **A full extension exists** → sort by (sam_flag, tx_id, diagonal) ascending; project each:
   `tx_start = ext.tx_lo` (full-length coverage implies read_lo==0), `txmap.project(tx_id,
   tx_start, read_len)`; score=0; strand = `to_genomic_strand` (Plus↔Minus flipped when the
   transcript is on the minus strand); finalize.
3. **Otherwise partials**: sort by (read_cov descending, tx_id ascending, diagonal ascending). Verify the top
   **min(8, len)**: for each, compute the EA-Myers infix distance (oriented read vs transcript window
   `[anchor - flank, anchor + read_len + flank)` clipped to sequence bounds, anchor =
   `ext.tx_lo - ext.read_lo`; skip if the window is empty); take the one with strictly minimal dist.
   No valid candidate → GateFail + Fallback.
4. cov = (read_len - dist)/read_len; **cov < 0.90** → GateFail.
   dist > **max(read_len/33, 1)** → GateFail.
5. **Cluster competition** (always on in production): cluster_bests folds partials by (tx_id, strand)
   into best-first clusters; when there are ≥2 clusters and dist ≥ 2: project the best cluster's locus
   (`tx_start = min(max(anchor,0), tx_len - read_len)`); compute dist/cov for each of the other clusters;
   if one passes the gates, maps to a different locus, and |dist difference| ≤ 2 → **place as usual but MAPQ=0** (option A / STAR semantics,
   authorized by users of the legacy project: the uncertainty of multi-mapping reads is expressed via MAPQ, with downstream filtering as the backstop).
6. **Repeat-region rejection** (always on in production): dist ≥ max_dist and the projection hits a RepeatTrack →
   likewise place + MAPQ=0.
7. finalize: project the best candidate (same tx_start formula as 5), score=dist; deterministic assignment below;
   if force_mapq0 → set MAPQ to 0.

## Deterministic Assignment (finalize_candidates)

Candidates sorted by (score ascending, tx_id ascending, diagonal ascending); best = first. MAPQ: unique candidate → 60;
runner-up score ≠ best → 60; all tied candidates have (contig, pos, cigar) exactly identical to the best → 60;
otherwise 0. `stats.extension_bases`: prefer the best candidate's cov (full branch) or ext.read_cov
(interrupted branch, written earlier by the engine).

Output TSV (for engine evaluation/differential testing): `name\tAligned\tcontig\tpos\tstrand\tcigar\tscore\tmapq`
or `name\tFallback\t*\t0\t*\t*\t255\t0` (column order: name, outcome, contig, pos,
strand, cigar, score, mapq; Fallback prints 255 in the score column and 0 in the mapq column).

## Adapter Layer: L1Index (resolves the I-1/I-2 dual-track minefield)

```rust
pub struct L1Index { /* tidx + full-coverage projection TxMap + transcriptome sequence store (same tx_id space) */ }
impl L1Index {
    pub fn build(tidx_path, gtf_path, ref_path) -> Result<Self>;   // build at runtime (differential testing/tests)
    pub fn open(bundle_path) -> Result<Self>;                      // production: open the build-ref artifact
    pub fn save(&self, bundle_path) -> Result<()>;
}
impl Tidx for L1Index { /* forwards to .tidx */ }
impl TxMap for L1Index { /* full-coverage projection */ }
impl TxSeqs for L1Index { /* rebuilt sequence store */ }
```

**Key invariant**: the projection and the sequence store cover **all** transcripts (no biotype filtering), and the id space =
tidx's (transcript_id in lexicographic order). `build()` implementation: parse the GTF with tidx's TranscriptSet
→ assert `tx_count == tidx.tx_count()`; per transcript → TranscriptRecord (exonic transcript sequence:
reversed for minus strand); `TxMap::from_records`; rebuild all transcript sequences from the reference FASTA (exons concatenated in genomic ascending order, whole-sequence RC for minus strand; any exon out of contig bounds → that transcript's sequence is empty). The filtered `.txmap` produced at build-ref time is **only** for aligner attribution — it does not participate in L1 projection.

bundle format (all little-endian, deterministic): magic "L1BNDL01" + version u32=1 +
[txmap binary (self-contained format, length-prefixed)] + [txseqs: tx_count u32, then per entry len u32 + bytes] +
[contig name table (contained inside the txmap blob)]. Writing the same input twice → byte-identical.

## Frozen Production Parameters (all environment switches removed)

| Legacy environment switch | Frozen production value |
|---|---|
| L1_SKIP_COVERED | Removed (full-information query) |
| L1_COVERAGE_GATE / L1_PLACE_INTERRUPTED | 0.90 |
| ESP_L1_VERIFY_TOP1 | Take minimum dist among top 8 before verification |
| L1_MAX_DIST | max(read_len/33, 1) |
| L1_CLUSTER_COMPETE(_BP/_MIN_DIST) | On; 2; 2 |
| L1_REPEAT_REJECT | On (RepeatTrack-driven) |
| L1_PLACE_MARGINAL | Place + MAPQ=0 (no legacy fallback) |
| PROF_PHASE / L1_DEBUG_* | Removed (no instrumentation in 1.0) |

## Explicitly Excluded

- PROF_PHASE phase timing (zero scientific logic, removed)
- synth/mock (test infrastructure, removed; verification runs in /tmp programs)
- engine-run binary and the io.rs output driver (belong to the cli crate, Wave 4)
- runtime.rs GTF runtime parsing (replaced by the L1Index bundle; zero parsing in production)

## Determinism

Same input → same output: no HashMap iteration reaches the output (diag_count is only queried, never iterated); all sorting uses explicit total-order keys;
reads are fully independent of each other.

## Self-Check Points (verified via /tmp programs, not in the repo)

- EA predicate direction table: 16 combinations + N/lowercase
- infix cross-validated against naive EA-DP (random pairs ≥10K, including N, indels, >128bp via long)
- Full-extension read → score 0, MAPQ 60
- Junction-spanning read (extension interrupted) → infix verification → correct projection + RefSkip CIGAR
- Two tied candidates with the same projection → MAPQ 60; different projections → MAPQ 0
- Cluster competition: two different loci at the same score → MAPQ 0
- Repeat region + borderline dist → MAPQ 0
- read < k → Fallback
