# pile — 8-Dimensional Pileup Feature Extraction

Extracts 8-dimensional pileup features at given sites from a BAM. Semantics = an exact replica of the **pysam pileup(stepper="nofilter") default behavior** (feature_spec.json v1), including two hidden htslib stages: the PE overlap quality tweak and the maxcnt=8000 cap. This is model input; values must match the reference implementation bit-for-bit (rtol=0).

## Feature Definition

`[depth, A_count, C_count, G_count, T_count, mean_base_quality, strand_bias, mean_mapq]`

- depth = number of reads at the site column whose query_position is valid (CIGAR falls on M/=/X) and whose (possibly tweaked) quality is ≥13
- A/C/G/T_count = base counts of these reads at qpos (both upper- and lowercase counted); **N and other non-ACGT characters do not enter the four buckets, but do count toward depth/mean qualities/strand/MAPQ** (audit fix: the initial spec incorrectly stated N does not count toward depth; both the legacy version and pysam count it)
- mean_base_quality = mean phred at qpos (**f32 single-precision division**, bit-identical to the legacy version; missing qualities handled per the rules below)
- strand_bias = fraction of forward-strand reads (forward / depth)
- mean_mapq = mean MAPQ
- depth==0 → all-zero vector

## Constants

`MIN_BASE_QUALITY = 13`, `MAXCNT = 8000`, `MERGE_GAP = 2000` (batch mode).

## CIGAR Column State

For each read, compute the state at pos0 (0-based): M/=/X → Match(qpos) (consumes query+ref); D/N → DelOrRefskip (consumes ref only); I/S consumes query only; H/P consumes nothing. Reference length rlen = sum of M+D+N+=+X.

## htslib bam_plp State Machine Replica (Core)

The push/emit interleaved timing must be replicated, because it determines two silent-drop paths:

1. **maxcnt drop**: at push time, if `iter->pos == b.pos` and `alive_count + 1 > 8000` → the read is silently dropped (and simultaneously removed from the overlap waiting table)
2. **reads with end <= pos never enter the buffer at all**

State machine structure (exact structure of bam_plp64_auto):
- Each auto call FIRST calls next(); only when next() has no column to emit (`max_pos <= pos`) does it enter the push loop, pushing while retrying next(), until a column is emitted or the stream is exhausted (push(0) sets EOF); after EOF no more pushes — next() alone drains the remainder
- next() = a single emit + advance: retire all nodes with `end <= pos` (decrement alive count, synchronously remove from the overlap waiting table); non-empty test = there exists a live node with `beg <= pos`; advance: `pos < hbeg (beg of first live node) → pos = hbeg`, otherwise `pos += 1`
- When the target column pos0 is emitted, capture its contents (live nodes with `beg <= pos0`, in push order)

The implementation may use heap optimizations (retirement heap on min-end, lazy cleanup of a min-beg heap, head pointers), but the semantics must not change.

## PE Overlap Quality Tweak (htslib overlap_push/tweak_overlap_quality)

- **Candidate condition**: proper-pair (0x2) and not mate-unmapped (0x8 clear), mtid==tid (or mtid<0), `abs(isize) < 2*l_qseq or mpos < rec_end`
- **Waiting table**: when a candidate read enters the buffer, if `mpos >= pos or mpos == -1` → enter the waiting table keyed by qname; when the second end with the same qname arrives, pair them, remove from the table, and tweak the qualities of both ends
- **Tweak rules** (a = left/arrives first, b = right/arrives later; directly rewrites the qpos quality of both ends):
  - Which end gets multiplied (amul/bmul) is decided by `wang_hash(x31_hash(qname)) & 1`; x31 = khash string hash (h = h*31 + c, u32 wrapping), wang = standard wang_hash(u32)
  - Both ends' CIGAR cursors advance synchronously into the overlap region; **del catch-up**: when one end's ref position lags behind and its preceding cigar op is D, the other end's quality is multiplied by 0.8 (on the amul end) or zeroed, until caught up
  - Bases equal: qualities summed and capped at 200; the amul end gets `amul*cap`, the bmul end gets `bmul*cap` (i.e., only the selected end keeps the sum, the other end is set to 0)
  - Bases differ: the higher-quality end ×0.8, the lower-quality end zeroed; when equal, the amul end ×0.8×amul and the bmul end ×0.8×bmul (equivalent to selected end ×0.8, other end 0)
  - CIGAR cursor semantics (cig_set/cig_next): S/I advance query only, H/P skipped, M/=/X advance query+ref synchronously, D/N advance ref only; malformed CIGAR (iseq out of bounds) exits tolerantly, must not panic
  - Unsupported cases such as refskip skip this round

## Quality Filtering Details

- When qpos < l_qseq take `qual[qpos]`, otherwise treat as 0; **BAM stores missing quality as 0xFF** (rust-htslib `qual()` returns a 0xFF-filled slice) → value 255, kept (≥13)
- The tweak only affects the two ends of a paired hit; unmatched reads use their original qualities
- Lazy copying is permitted for performance (clone qualities only for hit pairs); semantics unchanged

## Batch Interface

`extract_pileup_features_batch(bam, sites) -> Vec<[f32;8]>`: sort and group by (chrom, pos), merge groups with adjacent spacing ≤2000 into one fetch; scanline implementation (records stream in by coordinate with `pos <= pos0`, expired records with `end <= pos0` are evicted); the record set for each site = `{pos <= pos0 && end > pos0}` (file order), bit-identical to per-site fetch results. Output order = input sites order.

## Single-Site Interface

`extract_pileup_features(bam: IndexedReader, chrom, pos_1based) -> [f32;8]`; chrom not in BAM header → `ContigNotFound` error.

## Dependencies

`rust-htslib` (IndexedReader + .bai), thiserror.

## Self-Check Points

- Single read, single M covering the site: depth=1, counts per the qpos base, strand_bias 0 or 1
- qpos quality 12 → the read is excluded; 0xFF → included and contributes 255 to the mean
- PE fully overlapping, same base: selected-end quality = min(200, sum), other end 0 (usually filtered out by 13)
- 8001 reads with the same start: the 8001st is dropped by maxcnt
- D/N covering the site: the read is excluded

## Differential-Testing Corrections (post code-level audit + targeted differential testing against the legacy binary)

1. **del catch-up quality rewriting must be implemented**: the contract is the tweak_overlap_quality of the htslib bundled with pysam 0.24 (a 1.23.1 fork) — it preserves the legacy do-while catch-up: when one end advances due to D and the other end lags, the lagging end steps base by base to catch up, each step multiplying quality by 0.8 (selected end) or zeroing it; N (RefSkip) does not trigger catch-up and goes through continue. Note: upstream htslib (≥ some version) removed the catch-up (plain continue); the pysam-bundled version is authoritative.
2. **Main-loop stepping order**: both ends first step to the old iref, then max is raised and iref++ — if the a end's max raise is inserted before b's step, b silently catches up, the catch-up never triggers, and the comparison is misaligned.
3. Edge alignment: reaching the end of cig_next inside catch-up → return (not break); iseq out of bounds → tolerant return without panic (the new version is safer than C's out-of-bounds write, behaviorally equivalent).
4. N/non-ACGT bases count toward depth and all means (only the four buckets exclude them) — already reflected above.
5. Means use f32 division (legacy `sum as f32 / n as f32`); f64 division then conversion to f32 can differ by 1 ulp at rounding boundaries.
6. Quality slice bounds safety: gate on qual.len() (short slice → 0xFF kept, bq contributes 0), do not index directly by seq_len.
7. Intentional deviations from the legacy version (recorded, not replicated): reads that are FUNMAP but have coordinates are skipped at push time per htslib/pysam (the legacy version counted them); site pos<1 raises a clean error (the legacy version relied on fetch failure).
