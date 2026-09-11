# map — Genome-Level Splice-Aware Aligner (Track 1 Main Track + Track 2 Junction-kmer Direct Track)

The alignment crate of ESPERANTO 1.0.0. Responsibilities: FASTQ → (L1 transcriptome engine fast path, engine crate) →
L2/L3 genomic fallback path (minimap2-style seed/chain/extend + RNA splice-aware extension) →
BAM / unmapped.fq.gz / align_qc.json / align.baln. Track 2 (jkmer) re-locates reads directly
through annotated junctions when Track 1 leaves them unaligned or heavily soft-clipped.

This spec is the single source of truth. The legacy research tree is a read-only reference/porting source (engineering convention v1.1); porting provenance is recorded outside this document. When spec and implementation disagree, the spec is written back.

## Hard-Won Lessons (mandatory constraints at the head of the spec)

1. **Validate direction-related conventions on both the write and read sides** (legacy: minus-strand SEQ must be revcomp'd + QUAL reversed — a write-side root cause).
2. **When two references diverge, follow the one that generates the training data** (the htslib vs pysam fork precedent; BAM encoding follows the htslib/noodles layout).
3. **A performance fast path is semantics**: optimizations with behavioral consequences (dedup cache, Kadane zero-mismatch fast lane, EA-Myers decisively skipping DP) must be written into the spec.
4. **Spec errors propagate**: inconsistencies discovered during implementation must be written back into this file — no silent "fixes".

## Frozen Production Parameter Table (every value is a hard acceptance gate; changes require re-running the differential tests)

| Module | Parameter | Value |
|---|---|---|
| seed | RNA SeedParams | k=15, w=5 |
| index | HIGH_FREQ_FRAC (high-frequency seed truncation rank fraction) | 0.0002 |
| index | position packing | `contig(31b)<<33 \| strand(1b)<<32 \| pos(32b)` |
| index | paidx version | 1 (magic `PAIDXfmt`) |
| chain | ChainParams::rna_default | max_gap=1_000_000, min_chain_score=40, k=15, rna=true, min_intron=20, intron_penalty=8 |
| chain | TOP_CANDIDATES (decisive-round candidate cap) | 24 |
| extend | ExtendParams | match=2, mismatch=−4, gap_open=4, gap_ext=2 |
| extend | production editing_aware | true (explicitly enabled by align-prod) |
| extend | acceptance threshold | aligned segment ≥85% perfect (perfect = full-score bases within M) |
| extend | banded initial width | 30, doubled on failure up to the full band |
| align | coverage gate | 0.7 (aligned fraction of the read below this ⇒ None) |
| splice | SpliceParams | min_intron=20, max_intron=1_000_000, gt_ag_bonus=4, gc_ag_bonus=2, at_ac_bonus=1, noncanonical_penalty=8, known_bonus=6, refine_radius=12 |
| splice | SPLIT_COST (split rewrite cost) | 4 |
| intron_chain | IntronParams | min_frag_cov=15, refine_window=30, refine_pattern_len=50, intron_max=105_000 |
| intron_chain | frag clustering tolerance | implied read-start ±5bp |
| intron_chain | dominant-frag gate | ≥55% |
| intron_chain | best-pair score | `cov*4 + anchors*3 − bal` |
| intron_chain | verify_split_read | per-segment mismatch rate ≤0.17 and sum of both segments ≤0.24 |
| split | TAIL_MIN (library-driven tail floor) | 7 |
| split | library window | donor/acceptor ±20bp |
| split | tail-start retry | ±5bp |
| split | lib_bonus | `min(support,10)*5` |
| split | seed channel | tail ≥15, take the first 8 minimizers, `0<count≤100`, take the first 64 positions |
| split | direct-scan channel | tail 5..15bp, donor ±12, downstream ≤50_000, exact prefix probe ≤10bp |
| split | de-novo channel disabled when | the library is non-empty (both seed and direct-scan channels disabled; 97% of false junctions come from support=1 de-novo) |
| split/span | tail extension acceptance | `aligned ≥ tail_len*4/5 and score>0` |
| span | donor range | `(aln.pos+2, ref_end−1)`, candidates in ascending start order, first success returns |
| span | micro tail | <4bp full match accepted directly; ≥4bp goes through extension |
| mapq | formula | `round(60 × (1 − s2/s1) × min(1, n_anchors/10))`; chain_score≤0 ⇒ 0 |
| pair | proper-pair span cap | 1000 |
| pair | InsertStats | count/sum/sum_sq (f64), push takes \|tlen\|, stdev = population standard deviation `sqrt(E[x²]−m²)` (max 0) |
| pipeline | CHUNK (streaming block) | 100_000 reads |
| pipeline | 2-pass filter | support≥2 and mean MAPQ≥20 (mqsum ≥ c×20) and span ≤500_000 |
| pipeline | SOFTCLIP_THRESHOLD (Track 2 routing) | 20 (total soft-clip > 20 triggers) |
| pipeline | TRACK2_LOCAL_THRESHOLD | 30 (local_confirm.total floor) |
| jkmer | K / FLANK / MAX_A_FOR_VARIANTS | 15 / 16 / 6 |
| jkmer | candidate gate | ≥2 k-mer hits on the same junction |
| jkmer | breakpoint inference | mode of `read_pos+split_offset`, ≥2 votes within ±1bp ⇒ HighConf |
| jkmer | local confirmation | affine SW: gap_open=−4, gap_ext=−1; sub: (G,A)=0, (C,T)=0, match=5, transition=−1, transversion=−4; SW clamped at 0 |
| jkmer | file | magic `JKMER01\0`, version=1, trailing sha256(body) |
| bam | UNMAPPED_BIN | 4680 |
| bam | fixed tag order | XS (only when junctions present) → EA (only >0) → EK (always present when aligned) → RE (only rescued, value `unmapped`) |
| bam | BGZF compression workers | `clamp(threads,1,4)` |
| fastq | .gz detection | by extension; MultiGzDecoder |
| bfq | header | magic `EBFQ` + [4..8] reserved + [8..16] read_count u64 LE; record = u16 name_len + name + u16 seq_len + seq + qual |
| baln | header | magic `ESPBALN\x01` + u32 contig_count + per contig (u8 name_len + name); record = u32 block_size + BAM record bytes (no BGZF) |

## Module Graph and Dependencies

```
fasta ── seed ── index(+index_io) ── chain ── extend ── myers_ea
  │                                        │
gtf(JunctionLib) ── splice ── split / span / intron_chain
  │                                        │
jkmer (Track 2 standalone system, does not reuse extend)   mapq / pair / evidence / stats
  └────────── align (main flow) ── pipeline (run_se/run_pe) ── bam(+bam_encode) / baln / fastq
```

- The `map` crate depends on the `engine` crate (L1 fast path) and workspace-shared types; no reverse dependency.
- `prof`/`phase_legacy` (legacy instrumentation, ESPERANTO_PROF/PROF_PHASE) is **not ported**: 1.0.0 has no global atomic counters; performance observation uses external tools.
- The legacy `bin/esperanto-sim*` (data simulators) is not ported.
- Repository rule: **no tests, no fixtures inside the crate**; verification runs externally via differential testing against the legacy binary.

## Public Core Types

```rust
pub enum Strand { Plus, Minus }        // flip(); bam flag 0/16
pub enum CigarOp { Match(u32), Ins(u32), Del(u32), RefSkip(u32), SoftClip(u32) }  // SAM; =/X folded into M

pub struct ReadAlignment {             // all field semantics frozen
    pub contig: u32, pub pos: u32,     // 0-based inclusive
    pub strand: Strand,
    pub score: i32,                    // extension score
    pub chain_score: i32, pub second_chain_score: i32,  // MAPQ inputs
    pub cigar: Vec<CigarOp>,
    pub n_anchors: usize,              // number of chain anchors
    pub junctions: Vec<RefinedJunction>, // empty = non-spliced
    pub ea_count: u32,                 // number of editing-aware tolerated A>G/T>C (EA tag)
    pub mm_count: u32,                 // ordinary mismatch count (EK mm)
    pub n_seeds: usize,                // seed count (EK seeds)
    pub rescued: bool,                 // rescue/re-seed/Track2 (RE tag)
}
```

## fasta — Reference Genome 2bit

- `Base::{A,C,G,T,N}`; `from_ascii` case-insensitive; `code()` A=0 C=1 G=2 T=3, N→None.
- Packing: 4 bases per byte, **high bits first** (byte pos/4, shift `6−2×(pos%4)`); a partial tail byte is zero-padded in the low bits.
- N intervals: adjacent N runs are collected and **merged** at parse time; the 2bit layer cannot encode N, so the query side `base(pos)` checks the interval table first (hit → N) before unpacking.
- `Contig{name, len:u32, packed:&[u8], n_intervals:Vec<Interval{start,end}>}` (half-open).
- `slice_ascii(start,end)` / `decode_into` / `decode_append`: decoding goes through a 4-base LUT; positions covered by N intervals output `N`.
- `Reference{contigs, fasta_sha256:[u8;32]}`; `contig_index(name)`.
- `parse_fasta`: the `>` header's **first token before whitespace** is the name; empty name/empty file/duplicate name are errors with line numbers; sequence lines have whitespace stripped, illegal characters → N (not an error).

## seed — Minimizers (RNA High Density)

- `SeedParams{k,w}`; `Minimizer{kmer:u64, pos:u32, strand:Strand}`; `Anchor{qpos,rpos,contig,strand}`.
- Encoding: 2bit canonical (the smaller of kmer and revcomp); strand records the original orientation.
- Sort key: **mix64 = splitmix64 finalizer mix** (not the raw kmer value).
- Windows: the minimum mix64 among w adjacent k-mers; **monotonic queue** O(n); k-mers containing N (code 0xFF) are disqualified as candidates, **but the window still picks its minimum among valid members** (only fully invalid windows emit nothing); ties take the **leftmost**; the same minimizer across adjacent windows is emitted only once.
- `minimizers(seq, params)` for reads; `minimizers_from_codes(codes, k, w)` for index building (codes pre-decoded 2bit, N=0xFF).

## index — Minimizer Index

- `Index{params, version, freq_cutoff, reference, kmers:&[u64], offsets:&[u64], counts:&[u32], positions:&[u64]}`; `Box::leak` at build time, leaked mmap at load time; read-only shared.
- Build: per contig, **rayon-parallel** decode + N masking (intervals stamped back to 0xFF) → minimizers_from_codes → (kmer, packed_pos); merged with `sort_unstable` (tuple order ⇒ positions within the same kmer ascending by packed value); deduplicated into tables.
- `freq_cutoff = compute_cutoff(counts)`: deduplicated counts in descending order, the value at rank `len×0.0002` (capped at len−1), at least 1; empty table ⇒ u32::MAX.
- `query(kmer)`: binary search, returns an `IndexHit{count, positions}` slice, zero allocation.
- `unpack_pos(p) → (contig=u32(p>>33), pos=u32(p), strand=bit32)`.
- `collect_anchors(index, read_mins, read_len, k, occ_cap)`: skips `count==0 || count>occ_cap`; anchor orientation = XOR of read and ref strands (same→Plus, different→Minus); **qpos is determined by the anchor orientation (not m.strand)**: Plus anchor qpos = m.pos, Minus anchor qpos = `read_len − m.pos − k`. Returns (anchors, seeds_hit); seeds_hit = number of read minimizers that passed the cap. Positions are **not count-limited** (all emitted).
- `collect_anchors_edit_variants` (single-site): for each read minimizer's k-mer (read forward-packed fwd; if canonical is Minus, first revcomp to restore), for each G(→A)/C(→T) position i (slot = k−1−i) run **two rounds, each producing a different key**: the fwd round writes alt into the slot of the fwd packing (editing read position i), queries canonical, anchor astrand = ref_strand; the rev round writes 3−alt into the **same slot of the rev packing** (i.e., editing the **mirror read position k−1−i**, a different key), **skipped as no-op** (rev slot already equals 3−alt, i.e., mirror-position read base == alt), anchor astrand = flip(ref_strand). When both rounds each satisfy `0<count≤occ_cap`, var_hits += 1; take the **first 64** positions; qpos follows the anchor-orientation rule (Plus=m.pos, Minus=read_len−m.pos−k); no deduplication. Returns (anchors, var_hits).
- `collect_anchors_edit_variants_double` (two-site deep water): pairwise combinations of editable positions ((slot,alt) pairs taken in ascending read-position order), **each orientation's packing has both positions rewritten** (rev round: both 3−alt = editing the mirror position pair); canonical query; **no no-op skip, no deduplication**; take the **first 32** positions; var_hits not counted; returns Vec<Anchor>.

## index_io — paidx v1 Serialization (layout frozen, byte-level deterministic)

```
magic[8]="PAIDXfmt" | version u32 | k u32 | w u32 | freq_cutoff u32
contig_count u32 | kmer_count u64 | positions_count u64 | fasta_sha256[32]
per contig: name_len u32 | name | seq_len u32 | packed_len u32 | packed
          | n_interval_count u32 | n×(start u32,end u32)
pad1 = (8 − (60+contig_bytes)%8)%8          // contig_bytes summed per the parse formula
kmers ×u64 | offsets ×u64 | counts ×u32
pad2 = (8 − (60+contig_bytes+pad1+kmer_count*16+kmer_count*4)%8)%8
positions ×u64
```

- All little-endian; no timestamps/random ordering ⇒ two builds of the same input are **byte-identical** (acceptance gate).
- load: mmap + leak; strict validation of magic/version (mismatch ⇒ IndexVersion{file,supported})/per-segment lengths (Cursor-centralized bounds checks)/`packed_len == ceil(seq_len/4)`/trailing bytes must be 0; cross-table consistency (kmers ascending, offsets/counts equal length, o+c ≤ positions.len).
- Zero-copy: kmers/offsets/counts/positions point into the mapped pages via `bytemuck::cast_slice` (alignment guaranteed by the pads).

## chain — Anchor Chaining (minimap2 DP)

- `ChainParams` (see parameter table); `Chain{contig, strand, anchors (ascending by (rpos,qpos)), score}`.
- Input anchors are first sorted and grouped by (contig, strand, rpos, qpos); DP runs within each (contig,strand).
- DP transition: `f(i) = k + max_j[ f(j) + min(dq, dr) − gap_pen ]`; j ranges over anchors before i with **both rpos and qpos strictly increasing** (`a.rpos ≤ b.rpos` or `a.qpos ≤ b.qpos` → not chainable); `dr > max_gap` → scan terminates (rpos ascending), `dq > max_gap` → skip; no legal j ⇒ f(i)=k.
  - `dq = qpos_i − qpos_j`, `dr = rpos_i − rpos_j`, `dd = |dr − dq|`.
  - gap_pen: DNA = dd/2; **the RNA flat intron penalty applies only to true junction geometry**: `rna && dr > dq && dd > min_intron && dq ≤ k+10` ⇒ intron_penalty=8; otherwise dd/2 (read-side gaps are heavily penalized, preventing the DP from earning min(dq,dr) for free by skipping whole read segments).
  - `dq = qpos_i − qpos_j`, `dr = rpos_i − rpos_j`, `dd = |dr − dq|`.
  - gap_pen: DNA = dd/2; **RNA with dr−dq > min_intron (intron-type gap) ⇒ fixed intron_penalty=8** (replacing dd/2, which would be devastating for a 1M gap).
- Backtracking: in descending f order, backtrack from unused endpoints; **hitting an already-used anchor during backtracking ⇒ the whole chain is discarded (no anchor is marked)**; otherwise mark and emit; `f[end] < min_chain_score` skips. Global ordering: score descending → contig → strand → first-anchor rpos → first-anchor qpos.
- `second_score` = the highest-scoring chain with **no (qpos,rpos) anchor overlap** with the best chain (MAPQ input, 0 if none); the chain list is descending, so the first non-overlapping one is the maximum.
- `tied_best`: the set of chains tied with the best score (for multi-mapping determination).
- The caller takes only the first TOP_CANDIDATES=24 candidates into the decisive round.

## extend — Gotoh Banded DP (editing-aware)

- `ExtendParams` (see parameter table); when `editing_aware=true`, substitution goes through **SubstLut 5×5**: read relative to reference (plus-strand coordinates) A>G / T>C penalized 0 (A>I editing not penalized), all other mismatches = −4.
- `Extension{read_start, read_end, ref_start, cigar, score}`; unaligned parts at both read ends are written as SoftClip (end-trimmed, conserved: the M/I/S consumption of the CIGAR equals exactly the read length).
- `extend(read, ref_window, params, buf)` full window; `extend_hint(read, ref_window, params, buf, DiagHint{offset,num,den})` banded along a diagonal hint.
- The band schedule has **two tiers** `[30, usize::MAX]` (not doubling): the first tier has half-bandwidth 30; **tier-upgrade gate**: `ext.score*20 < perfect*17` (perfect = read_len × match_score; i.e., score < 85% of perfect) ⇒ second tier recomputes the full matrix; the second-tier result is returned unconditionally. band ≥ ref_window.len() counts as already full-matrix and is accepted directly.
- **Path dispatch**: if every byte of read and ref_window ∈ {A,C,G,T,N} (uppercase) ⇒ `run_banded_packed` (i16 compact layout + substitution-score LUT); any other byte (lowercase/degenerate) ⇒ `run_banded_legacy` byte-wise path. Both paths produce the same substitution score for every input byte pair, and the output Extension is identical score-by-score and CIGAR-by-CIGAR (acceptance gate; packed is a pure performance rewrite).
- DP equations: H/E/F three-state affine (open=4, ext=2); the i16 packed implementation (run_banded_packed, LUT-accelerated) and the legacy implementation **must agree score-by-score** (packed is a pure performance rewrite); backtracking prefers H, soft clipping allowed at the ends.
- Acceptance: **the extension result has no independent acceptance gate** — the caller consumes Extension{score, cigar, …} directly; the only gate is the in-band tier-upgrade gate above. (The `aligned ≥ tail×4/5` of split/span is an alignment-length gate of those modules, not part of extend.)
- `ExtendBuffer` reuses DP row buffers (thread-local, no semantics).
- `push_op` merges adjacent ops of the same kind.

## myers_ea — EA-Myers Bit-Parallel Verification (decisive fast lane)

- Word size: **u128** (single block, pattern ≤128bp). `mask_for(m)`: sets the low m bits; special-case m≥128 as u128::MAX (avoids `1<<128` overflow).
- `build_peq(read) -> [u128;256]`: bit i of peq[ref_byte] is set if and only if ea_match(ref_byte, read[i]); EA merge rules frozen: read=A→{A}; read=C→{C,**T**}; read=G→{G,**A**}; read=T→{T}; others→empty.
- `calculate_block(pv, mv, eq_in, hin, mask, high_bit)`: edlib-style carried arithmetic: hin_is_neg=(hin<0)?1:0; ph_bit0=(hin>0)?1:0; xv=eq_in|mv; eq=eq_in|hin_is_neg; xh=(((eq&pv)+pv)^pv)|eq (wrapping); ph=mv|!(xh|pv); mh=pv&xh; hout: ph&high→+1, mh&high→−1; ph_sh=(ph<<1)|ph_bit0; mh_sh=(mh<<1)|hin_is_neg; pv_out=(mh_sh|!(xv|ph_sh))&mask; mv_out=(ph_sh&xv)&mask.
- `infix(read,text)`: init vp=mask, vn=0, score=m, best=m; advance per text character (hin=0), score+=hout, best takes the strictly smaller value; returns best. `infix_best_end` traces the same way and records the best end (strict < updates). `infix_best_start`: reverse both read and text and reuse infix_best_end, start = len−1−rev_end.
- **Long patterns (128<m≤256) go through the long module, 2-block**: peq0/peq1 split by i<128 / i−128; block0 mask=u128::MAX, high0=1<<127; block1 mask1=mask_for(m−128), high1=1<<(m1−1); init vp0=MAX, vn0=0, vp1=mask1, vn1=0; per character, block0 (hin=0) → **hout0 chained as block1's hin**; score+=hout1. `long::infix_best_end` is isomorphic; `long::infix_best_start` reuses via reversal. For m≤128, long dispatches back to the single-block version.
- m>256: the guard returns m (distance = pattern length, which no verification threshold will pass ⇒ the caller falls back to banded DP); the legacy version would shift-overflow here, the new version guards explicitly (in practice reads are ≤250bp, unreachable).
- Purpose: distance verification of decisive candidates at the decisive stage (accepted if ≤ max_dist, skipping banded DP); on failure, fall back to extend.

## mapq — Mapping Quality

```
mapq = round( 60 × (1 − second/best) × min(1, n_anchors/10) )
best ≤ 0 ⇒ 0; second=0 ⇒ factor 1.
```
- Inputs come from ReadAlignment{chain_score, second_chain_score, n_anchors}; **Track 2 synthetic records also go through this same formula** (see the pipeline::try_track2 note).

## pair — PE Relation and Insert Statistics

- `PairRelation`; `relate(r1, r2)`: R1 on the left and Plus, R2 on the right and Minus, span ≤1000 ⇒ proper pair.
- `template_length`: signed, from R1's perspective (left end positive, right end negative).
- `ref_end(aln)`: advanced by CIGAR M/D/N.
- `InsertStats{count, sum, sum_sq}`: push takes |tlen|; mean/stdev per the parameter table (population standard deviation).

## gtf — Junction Library

- `Junction{contig:u32, start:u32, end:u32, minus_strand:bool}` (0-based half-open intron interval); `RefinedJunction{junction, signal:SpliceSignal, known_support:u32}`.
- `JunctionLib{junctions (sorted), counts (parallel), by_end (endpoint index)}`;
  - `build`: sort + merge counts of identical junctions;
  - `range_start(contig,lo,hi)` / `range_end(contig,lo,hi)`: range slices (same order, same set; replaces full-library linear scans);
  - `support(j)` / `contains(j)`; `is_empty()`.
- `from_gtf(path, contig_id_fn)`: parses GTF exon lines, aggregates per transcript into introns between adjacent exons; strand determines minus_strand; attribute parsing takes transcript_id.

## splice — Splice Signals and Pseudo-Reference Alignment

- `SpliceSignal::{GtAg,GcAg,AtAc,NonCanonical}`; `score(p)` = gt_ag_bonus/gc_ag_bonus/at_ac_bonus/−noncanonical_penalty; `label()` = "GT-AG"/"GC-AG"/"AT-AC"/"-".
- `splice_signal(donor[2], acceptor[2], minus)`: plus strand GT-AG/GC-AG/AT-AC; minus strand tests the mirror (CT-AC etc.); none matching ⇒ NonCanonical.
- `refine_junction(reference, contig, naive_start, naive_end, minus, lib, params) -> Option<RefinedJunction>`: enumerates breakpoints within ±refine_radius(12), score = `signal.score × 100 − dist + known_bonus (library hit)`, takes the best; output carries signal and known_support.
- `align_spliced` (**multi-segment**, pseudo-reference stitching): `segment_chain(chain, k=15, min_intron)` cuts segments; <2 segments ⇒ None (caller falls back to the DNA path). For each adjacent segment pair, `refine_junction` (**supports ≥2 junctions**); each pair passes the **S14 guard**: intron <50 with no library support ⇒ None; NonCanonical with no library support ⇒ None. Pseudo-reference = concatenation of each segment's ref slice (first segment extended left by flank=30, last segment extended right by flank=30, middle-segment boundaries use the refined junction coordinates; re ≤ rs ⇒ None); main-diagonal hint extend_hint; `stitch_cigar` maps back to genomic coordinates per exon_bounds, inserting RefSkip at segment crossings (true intron length); SoftClip passed through verbatim; returns SplicedAlignment{extension, junctions, cigar, exon_bounds}.
- `resolve_split_right/left`: tail-relocation decisive round; the total formula includes SPLIT_COST=4 and `intron_length_penalty`; `score_tail_win` validates the tail by extend_hint within the candidate acceptor/donor window. **Exact intron_length_penalty formula (frozen)**: `(8.0 × ((len.max(1) as f64).log10() − 3.0).abs()) as i32` — linear |log10| distance, 0 at 1kb (not quadratic).

## intron_chain — MEM Intron-Chain Fast Path (two-segment splicing)

- `IntronParams` (see parameter table).
- Flow: MEM extension on the read → Frag (clustered by shared implied read-start ±5bp) → dominant-frag ≥55% gate → best segment-pair score `cov*4 + anchors*3 − bal` (bal = |cov_a − cov_b|) → middle intron-gap scan (with fallback) → EA-Myers breakpoint refinement (refine_donor_end uses `infix_best_end` / refine_acceptor_start uses `infix_best_start`, window frag end ±refine_window=30, pattern length refine_pattern_len=50; **refinement results are adopted directly, with no extra acceptance gate such as signal-score comparison**) → `verify_split_read` dual thresholds (per-segment mismatch rate ≤0.17, sum of both ≤0.24; the only acceptance stage) → conservativeness guard → `try_intron_chain_placement` output.
- Any stage failing returns None; the caller proceeds to later paths (no error).

## split — Tail Split Rescue (STAR-style split discovery)

- Trigger: main alignment has terminal soft-clip ≥ TAIL_MIN (library-driven; the seed channel additionally requires tail ≥k=15). Three channels; **when the library is non-empty, the de-novo channels (seed/direct-scan) are all disabled**.
- Channel A (library-driven, main channel): `lib_window` (right tail, start±20)/`lib_window_end` (left tail, end±20) enumerates library junctions; filters min/max_intron and strand; NonCanonical with support<2 ⇒ skip; `retry_tail_start` retries around the tail start ±5bp:
  - For each dro: left segment read[..dro] extends over window `[j.start−len−30, j.start)`, right segment read[dro..] extends over window `[j.end, j.end+len+30]`;
  - Each segment counts only if aligned ≥ len×4/5; total score = le.score + re.score + signal_score + lib_bonus − SPLIT_COST; take the maximum.
  - Success ⇒ `build_rescue_ext`: CIGAR = [left extension] N(intron) [right extension] (completely independent of the original CIGAR), new_pos = left segment's ref start.
- Channel B (direct-scan, only when the library is empty): tail 5..15bp; donor ±12 scan, searching downstream ≤50kb for an exact prefix match of the tail probe (≤10bp) + acceptor AG/AC ending; score = `sig.score×100 − |d−primary_end| − |intron−5000|/1000`; the best one gets whole-tail extension validation (≥4/5 and score>0); the emitted signal is recorded as GtAg, known_support=0.
- Channel C (seed-driven, only when the library is empty): tail ≥15; query the first 8 minimizers (count ≤100, positions ≤64); same contig, same strand, landing in the downstream min/max_intron window; refine_junction refinement + S14 guard (NonCanonical without library rejected; <50bp intron without library rejected); extension validation as before.
- Left tail = symmetric implementation (acceptor = body start, donor upstream; `build_rescue_left` = [tail extension] N [body with head soft-clip removed], pos shifted forward).
- `build_rescue` (right-tail stitching): the body has its **tail-end** soft clip removed (judged by read_end), everything else kept + N + the tail extension as a whole (conserved).

## span — Run-Through Reinterpretation (library-driven only)

- Scenario: a read with a 1–4bp overhang is aligned straight across a junction by extension (no soft clip, so split cannot trigger).
- Only when the library is non-empty; de-novo run-through is **not done**.
- Candidates: donors falling in (aln.pos+2, ref_end−1) (ref_end advanced by CIGAR M/D/N), same strand, intron ≤ max_intron, ascending by start, **first success returns**.
- The read offset dro corresponding to the donor is found by re-walking the CIGAR; tail = read[dro..]:
  - <4bp: base-by-base full match at the acceptor ⇒ accepted directly, CIGAR tail segment Match(tail_len);
  - ≥4bp: extension window [j.end, j.end+len+30], aligned ≥ 4/5 and score>0.
- `rebuild`: the body is truncated at dro (including truncation across an op; Del/RefSkip do not consume read, so landing inside one keeps the whole op), the head soft clip is kept and placed first, + N + tail CIGAR.

## align — Single-Read Main Flow (order frozen)

- `AlignConfig{seed:SeedParams, chain:ChainParams, extend:ExtendParams, rna:bool}`; `rna_default()` = k15/w5 + chain rna + extend defaults + rna=true; production additionally sets `extend.editing_aware=true` explicitly.
- `Aligner::new(&Index, AlignConfig)`; `align_read(read: &[u8]) -> Option<ReadAlignment>`. Flow (any stage succeeding returns immediately, failures fall through in order):

1. **L1 transcriptome fast path** (only when `l1: Option<Arc<L1Runtime>>` is provided): call engine::align_read; on successful placement, wrap into a ReadAlignment via the **conservativeness guard** (CIGAR conservative validation) and chain_score reverse-engineering (MAPQ compatibility); on `Fallback`/validation failure, fall to the genomic path.
2. **Normal seed round**: `minimizers(read, seed)` → `collect_anchors` (occ_cap = index.freq_cutoff) → `chain_anchors` (chain_round): yields best / second_score / tied_best.
3. **Edit-variant round** (when the normal round is weak): `collect_anchors_edit_variants` single-site → double-site when necessary; hit anchors are merged in and re-chained.
4. **Decisive test**: `tied_best.len() ≤ 1 and second_score×5 < best.score×4` ⇒ decisive; otherwise a competition set (top 24).
5. **intron-chain fast path** (two-segment splicing; returns directly on success).
6. **try_spliced** (pseudo-reference stitching; RNA mode).
7. **EA-Myers decisive verification**: for a decisive candidate, myers_ea::infix distance ≤ max_dist ⇒ accept directly, skipping banded DP; failure (including indel signs/out-of-window/total rejection) falls to 8.
8. **Banded DP decisive round**: loop candidates with extend_hint; first `try_fast_lane` (Kadane zero-mismatch scan, full-match reads skip DP); take the best extension.
9. **split/span rescue**: terminal soft-clip ≥7 ⇒ rescue_right/left_tail (library-driven first); if the library is non-empty, also try rescue_span; success ⇒ rewrite CIGAR, junction recorded.
10. **recount + coverage gate**: `recount_mm_ea` recounts mm/ea from the final CIGAR; read coverage < 0.7 ⇒ None.

- In the PE context there is also `rescue_with_mate_anchor`: when one end is unaligned, locally re-align using the mate anchor.
- Chain failure (no legal chain) ⇒ None (unmapped).

## pipeline — Streaming Parallelism and Output Contract

- `run_se(pipe, r1_path, threads)` / `run_pe(pipe, r1, r2, threads)`; CHUNK=100_000 reads/block, rayon-parallel alignment, **output order = input order** (collected within a block, then written in order).
- `PipelineOut{bam: Option<Box<Write+Send>>, unmapped_fq: Box<Write+Send>, index, config, jlib: Option<Arc<JunctionLib>>, jkmer: Option<Arc<JkmerIndex>>, l1, baln: Option<Box<Write+Send>>}`.
- **SE dedup**: identical seq within a block is aligned only once (readcache), sharing the result (semantics: same sequence, same output).
- **PE**: processed as pairs; mate-anchor rescue; `pair_records` assembles FLAG (PAIRED/PROPER_PAIR/READ1/READ2/REVERSE/MATE_REVERSE/UNMAPPED/MATE_UNMAPPED) and mate fields (mc, mpos, tlen).
- **SEQ orientation (hard contract)**: minus-strand aligned records have SEQ = `revcomp(original read)`, QUAL reversed in step; **unmapped keeps original orientation**; SE/PE consistent. Direct samtools pileup reading depends on this convention.
- **unmapped double-write**: into BAM (FLAG 0x4) and into unmapped.fq.gz (hard contract of the rescue channel).
- **2-pass junction discovery** (without --gtf): pass 1 runs with an empty library (de-novo channels open) collecting (junction, support, mapq) → `filter_discoveries` (support≥2 ∧ mean MAPQ≥20 ∧ span≤500_000) → pass 2 reruns with the library; with --gtf, a single pass.
- **Track 2 routing (try_track2)**: queried when jkmer is provided and (aln=None or total soft-clip>20); for each candidate, infer_breakpoint + fetch_flanks + local_confirm, accepted only at `total ≥ 30`, taking the maximum total; a hit ⇒ replace with a synthetic ReadAlignment{cigar = [M(left)] N(intron) [M(right)], pos per the jkmer::assemble rules, chain_score=total, second=0, n_anchors=hits, junctions=empty, rescued=true}.
  - **Note (quirk, frozen by differential testing)**: the legacy implementation computed compute_track2_mapq and then **discarded it**, writing MAPQ = the standard mapq formula applied to the synthetic record (= round(60 × min(1, hits/10))). The new implementation preserves this behavior; compute_track2_mapq remains a public jkmer function (test/debug semantics complete) but the pipeline does not call it.
- **Statistics**: StatsAcc accumulates integers, finalized once at the end; mapq histogram has 61 buckets; elapsed_seconds is the sole determinism exemption.
- **baln**: written alongside the bam output (header + per record block_size+BAM bytes); fast-path encoding refusal ⇒ hard error (records are never silently dropped).

## jkmer — Track 2 Junction-kmer Index (standalone system, does not reuse extend)

Purpose: reads left unaligned or heavily soft-clipped by Track 1 are re-located directly through **annotated junctions**. Self-contained: index building, query, breakpoint inference, local confirmation, CIGAR/MAPQ, record assembly.

### Types and Constants

- K=15, FLANK=16, MAX_A_FOR_VARIANTS=6.
- `Junction{contig:u32, id:u32, intron_start:u32, intron_end:u32, strand:Strand}` (jkmer-local type, **not** gtf::Junction); `intron_len() = intron_end − intron_start`.
- `JkmerHit{junction_id:u32, split_offset:u8, a_mask:u16}`; `JkmerIndex{magic, version, gtf_sha256, fasta_sha256, junctions:Vec<Junction>, entries:BTreeMap<u64, Vec<JkmerHit>>}`.

### Building

- `pack_kmer`/`unpack_kmer`: 2bit (A=0 C=1 G=2 T=3), high bits first; `extract_kmers(read) -> Vec<(read_pos, packed)>`: sliding windows of K=15, windows containing N skipped.
- `compute_a_mask`: bitmask of A positions within the flank (u16); `enumerate_a_to_g_variants`: when A count ≤6, enumerate all A→G combinations (2^n) for editing-aware hits.
- `extract_junctions_from_gtf` → `build_junction_kmers`: for each junction, concatenate donor tail FLANK + acceptor head FLANK into a pseudo-sequence, extract all **breakpoint-crossing** k-mers, recording split_offset (the breakpoint's position within the k-mer) and a_mask; `fetch_flanks(junc, fetch_base, donor, acceptor)` reads the reference for both flanks (**read orientation**: minus strand revcomp'd).

### File Format (magic `JKMER01\0`, version=1, all LE)

```
magic[8] | version u32 | gtf_sha256[32] | fasta_sha256[32]
n_junc u32 | per junction: contig u32, id u32, intron_start u32, intron_end u32, strand u8 (0=+,1=-)
n_entries u32 | per entry: packed u64, n_hits u32, per hit: junction_id u32, split_offset u8, a_mask u16
trailing sha256(body)[32]
```
- load verifies the trailing sha256 first, then parses; magic/version/strand bytes/lengths strictly validated; corrupt files raise errors, never panic.

### Query and Decision

- `query_read(read) -> Vec<JunctionCandidate>`: extract_kmers queries entries one by one, aggregating (read_pos, split_offset, a_mask) by junction_id; **filter ≥2 hits**; hits within a candidate sorted by read_pos ascending; candidates ordered by **hit count descending, junction_id ascending** (deterministic).
- `JunctionCandidate::infer_breakpoint()`: estimate = read_pos+split_offset; after sorting, find the mode with the most votes within ±1bp; **≥2 votes ⇒ HighConf(pos)**, otherwise LowConf(pos); negative mode ⇒ None.
- `local_confirm(read, split, donor_tail, acceptor_head) -> LocalConfirm{left,right,total}`: left = local_align_score(read[..split], donor_tail) (split=0 ⇒ 0); right = read[split..] vs acceptor_head (split≥len ⇒ 0).
- `local_align_score`: three-matrix affine SW (M/X/Y), gap_open=−4, gap_ext=−1, clamped at 0 from below; `sub_score`: (read G, ref A)=0, (read C, ref T)=0, same base=5, purine/pyrimidine transitions (A↔G, C↔T)=−1, all others=−4.
- `build_track2_cigar(left, intron, right)`: [M?] N [M?] (zero-length segments not emitted); `cigar_to_string`.
- `compute_track2_mapq(best_hits, runner_up_hits, score_margin) -> u8`: tie (runner_up==best>0) ⇒ 0; `log2(hits)×10 + clamp(margin/10×5, 0..20) − min(ratio×20, 20)`, rounded, then clamped to 0..60. **(Not called by the pipeline; see the pipeline note.)**
- `assemble_track2_record(junc, read_len, read_split, best_hits, runner_up, margin) -> Track2Record`: cigar = build_track2_cigar; pos_1based: Plus = `intron_start − read_split + 1`; Minus = `intron_start − right_match + 1` (right_match = read_len − read_split); saturation guards against negatives.

## bam — Write-Out Contract (noodles)

- `build_header(index)`: SQ lines from the index contig table (name + len); `@CO` comment `esperanto-map M2 (v0.1)`.
- `BamRecord{name, flag, aln:Option<ReadAlignment>, mapq:u8, seq, qual, mate:Option<(i32,i32,i32)>}`.
- `record_se(name, seq, qual, aln)`: SEQ orientation rule (minus strand revcomp + QUAL reversed; unmapped original orientation); Plus flag=0, Minus flag=0x10, unmapped flag=0x4/mapq=0.
- `write_record`: **fast path first** (bam_encode::try_encode into a thread-local reuse buffer; `Some(Ok)` ⇒ hand-written block_size+bytes; `None` ⇒ fall to the noodles RecordBuf slow path `write_record_slow`; `Some(Err)` ⇒ error).
- Slow-path tags are **same order, same values** as the fast path: XS:A (only when junctions non-empty, takes junctions[0]'s strand) → EA:i32 (only >0) → EK:Z (always present when aligned) → RE:Z=`unmapped` (only rescued). Mate fields written only when mc≥0.
- `create_writer(w, header, threads)`: multithreaded BGZF, workers = clamp(threads,1,4); writes the header.

## bam_encode — Zero-Allocation Direct Encoding (byte-identical to noodles)

- Goal: output **byte-identical** to `noodles::bam::io::Writer::write_alignment_record`; unsupported shapes return None (caller falls to the slow path).
- Refusal conditions (⇒None): name length ∉ 1..=254 or name is `*`; name contains non-ASCII-graphic or `@`; cigar op count >65535; qual has bytes >93; contig > i32::MAX.
- Layout: 32-byte fixed header (ref_id, pos, l_read_name=name length+1, mapq (unmapped=255), bin, n_cigar u16, flag u16, l_seq u32, next_ref, next_pos, tlen) → name+NUL → cigar (`len<<4|kind`, kind: M=0 I=1 D=2 N=3 S=4) → SEQ 4-bit high-nibble-first (code table `=ACMGRSVTWYHKDBN`, same code for either case, unknown=N=0xF, odd tail padded 0 in the low nibble) → QUAL (empty ⇒ all 0xFF; length mismatch ⇒ Err) → tags (same order as the bam section).
- `reg2bin(start,end)` (0-based inclusive): standard §5.3 hierarchy (14/17/20/23/26 shifts); span=0 ⇒ UNMAPPED_BIN=4680.

## baln — Internal align→call Fast Channel

- magic `ESPBALN\x01` + u32 contig_count + per contig (u8 name_len+name); record = u32 block_size + bam_encode bytes (32 core + data, no BGZF).
- fast-path encoding refusal ⇒ **error** (never silently dropped). The call side can memcpy straight into `bam1_t`.

## fastq — Reading

- `RecordSource` trait unifies `FastqReader` / `BfqReader`.
- `FastqReader::open`: .gz by extension ⇒ MultiGzDecoder (1MB Buf); strict four lines: `@` header / seq / `+` / qual; `seq.len≠qual.len` ⇒ error; **name = first token before whitespace in the header, empty name ⇒ error**; all errors carry line numbers; EOF (empty) ⇒ None; truncation ⇒ error. Line endings strip `\n`/`\r`.
- `BfqReader`: mmap; magic `EBFQ` + [8..16] read_count u64; record = u16 name_len + name (first token) + u16 seq_len + seq + qual; reaching read_count ⇒ None.

## evidence / stats / error

- `evidence_tag(aln) -> String`, fixed field order:
  `seeds={n_seeds};chain={chain_score};sub={second_chain_score};dq={max(chain−second,0)};splice={junctions[0].signal.label() or -};mm={mm_count};ea={ea_count};mapq_src={second>0 ? chain_margin : unique}`
- `AlignStats` (align_qc.json contract; consumers ignore unknown keys): esperanto_map_version, mode, total_reads, mapped_reads, unmapped_reads, mapping_rate, proper_pairs (PE only, otherwise null), insert_mean/insert_stdev (same), mapq_hist[61], junctions_total, rescued_total, rescue_fail_total, elapsed_seconds.
- `StatsAcc`: push_read(mapped, mapq) integer accumulation; finalize(mode, pe, elapsed) emits the table once; `CARGO_PKG_VERSION` goes into the version field.
- `AlignError` (thiserror): FastaIo{path,source} / FastaFormat{line,msg} / IndexIo / IndexFormat{msg} / IndexVersion{file,supported} / IndexReferenceMismatch{msg} / FastqFormat{line,msg} / FastqIo. User-supplied paths **never panic**.

## CLI Contract (align-prod form; wired into esperanto-cli at Wave 4)

Arguments: `--index` (paidx) `--r1` [`--r2`] [`--gtf`] `--out` [`--threads`=0 all cores] [`--no-bam`] [`--l1-tidx --l1-txmap --l1-gtf --l1-ref` **all four given together or all absent**; any missing ⇒ error].
Output directory: `raw.bam` (omitted with --no-bam), `unmapped.fq.gz`, `align_qc.json` (pretty + trailing newline), `align.baln`. Exit codes: 0 on success, 2 on error (message to stderr).

## Exclusions (explicitly not ported)

- `prof.rs` / `phase_legacy.rs`: development-time instrumentation (ESPERANTO_PROF/PROF_PHASE global atomics), not ported in 1.0.0; performance observation uses external profilers.
- `bin/esperanto-sim*`: legacy data simulators; the external differential-testing setup has its own generators.
- In-crate tests/fixtures: forbidden by repository rule; verification is entirely external differential testing.

## Deviation Register and Acceptance Gates

- **Rescued-read strand (maintained)**: rescued reads without strand evidence are counted as fwd — a documented legacy decision, maintained, not reopened.
- **Track-2 MAPQ quirk (frozen by differential testing)**: see the pipeline note; the new implementation replicates the behavior exactly; registered as a model contract in docs/SCIENCE-DEVIATIONS.md.
- **I-1/I-2**: the L1 fast path goes through the engine crate's L1Runtime; the tidx/txmap tx_id dual-track is guaranteed by the engine-side contract; map only consumes `Aligned/Fallback`.

Acceptance (standing gates):
1. `cargo build --workspace` + `cargo clippy -p <map-crate> -- -D warnings` with zero warnings.
2. **Determinism**: same input, same threads, two runs: `align.baln` byte-identical, `unmapped.fq.gz` identical after decompression, `align_qc.json` equal key-by-key except `elapsed_seconds`; same checks across thread counts (1 vs 8).
3. SE/PE, with/without GTF (2-pass), with/without jkmer, with/without L1 flag combinations each pass once under the external differential-testing pipeline.
