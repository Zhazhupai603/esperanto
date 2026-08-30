# scan — strand-resolved candidate editing-site discovery (scatter engine)

Genome-wide scan from BAM or .baln, emitting candidates.bed (10-column contract). **No hard filtering**: strand-resolved statistics + soft score call_score (0–1) + evidence code EVID all pass through (var_freq/fwd_freq/rev_freq pass through), to be consumed downstream with thresholds. Semantic source: the scatter caller is the only engine.

## Inputs

- `--bam`: aligned BAM (requires .bai/.csi). Skips unmapped (0x4) / secondary (0x100) / supplementary (0x800).
- Reads tagged `RE:Z:collapsed` (collapsed-rescue placements) count toward depth/bq/mapq/junction
  accumulators only; their bases are alphabet-ambiguous (A==G, T==C) and never contribute variant
  evidence (`sites` counts). Same rule on the .baln path.
- `--baln`: binary fast path from the map product; **byte-identical output** with the BAM path (both sources share the same scatter kernel; block-window overlap semantics aligned with htslib fetch: pos ∈ [cs−max_span, ce) and pos+span > cs). contig-derived lengths are corrected against .fai (derived values used when fasta is absent or the contig is missing).
- `--fasta` (optional): reference base / homopolymer / junction direction voting; absent or contig missing → majority pseudo-ref (all sites on that contig fully counted, hp=0, junction does not vote). A fasta missing a contig must be pre-checked against the .fai list; it must not panic.
- `--gtf` (optional): strand-decision evidence, priority 3.
- `--gnomad` (optional): AF soft down-weighting (never deletes); single file (plain/.gz) or per-chrom directory (`gnomad.joint.v4.1.sites.{chrom}.vcf.bgz`, chr-prefix mapping, `.afidx` cache; a missing file for a primary chromosome is a hard error, non-primary chromosomes are neutral). Configured but path does not exist = hard error.
- `--lib`: unstranded (default) / stranded (dUTP-type; read direction maps directly).
- `--enable-cu`: C>U symmetric mode (off by default).
- `--min-call-score`: threshold only marks (appends `,MS` to evid); never deletes.

## Direction handling (hard contract)

Per the SAM specification, the BAM SEQ of an aligned record on the negative strand is already stored in reference-forward orientation (guaranteed by the writer side). The scan read side **passes through**: the scatter takes bases from the stored SEQ position by position, with no revcomp and no complement flipping; strand assignment splits into fwd/rev buckets solely by FLAG 0x10. Semantics match htslib pileup / the pile crate.

## Scatter-engine semantics (O(N+G) replacement for the bam_plp state machine)

Parallel over (contig, 32 Mbp blocks); each read's contribution is scattered directly into per-block arrays:

- **depth**: counted only when bq ≥ 13 (MIN_BQ) and the qpos base ∈ ACGT; split into depth_fwd/depth_rev buckets by 0x10. CIGAR M/=/X occupy a column; I/S consume query only; D/N consume reference but no column (not counted in depth).
- **bq_sum / mapq_sum**: accumulated in sync for entries counted into depth (bq missing quality counts as 0xFF; qpos out of bounds counts as 0).
- **Mismatch localization (dirty gating)**: only reads with EK tag `mm+ea > 0` (reads missing EK are treated as dirty; without fasta all reads are treated as dirty) are compared base-by-base against the reference: positions with a reference base (ACGT) record only mismatches; positions where the reference is N are fully counted. Counts go into sites[pos].cnt[4 bases × 2 strands].
- **Intron direction voting**: each N (refskip) segment of a read votes by donor/acceptor dinucleotides: GT-AG / GC-AG / AT-AC → plus; CT-AC → minus; anything else does not vote. Voting only occurs when a fasta is present. The read's n_plus/n_minus votes are scattered to every site covered by that read (junc_plus/junc_minus); the two end positions of N segments go into jbounds (for junction_dist).
- **maxcnt**: the scatter has no 8000 cap (unrelated to the legacy differential-testing acceptance corpus; by design, under scatter semantics site depth does not exceed the cap — if it ever triggers, record it faithfully in SCIENCE-DEVIATIONS).

## Candidate assembly (per site)

- depth==0 is skipped; full_tally (no reference) = full-count mode of sites, with matched filled by majority; in normal mode matched = depth − mm_total, added into the ref bucket per strand.
- var_reads = depth − total_cnt[ref]; **var_reads == 0 is skipped**.
- edit_frac: plus strand looks at A>G, minus strand looks at T>C (amb takes the max of both directions; with enable_cu, union with C>U).
- Strand-decision priority: 1) LIB=stranded → fwd/rev depth ratio; 2) junction votes (consistent votes rule; conflicts → amb); 3) GTF (genes on both strands → amb); 4) otherwise amb. Primary evidence code LIB/JUNC/GTF/NONE.
- Features: depth, edit_frac, mean_bq=bq_sum/depth, mean_mapq, strand_bias=|var_fwd−var_rev|/max(1,var_fwd+var_rev), gnomad_af, hp_len (same-base run, capped at 10, only ACGT counted), junction_dist (nearest distance in jbounds; used as evidence only when ≤8).
- call_score: logistic weighted sum; weights come from a versioned spec JSON (built-in call_spec.v2.json real_v2 profile; overridable with --spec). depth_x=min(d,50)/50; bq_x=(bq−10)/30 clip01; mapq_x=mapq/60 clip01; gnomad_x=min(af×20,1) (0 when absent); hp_x=min(hp,10)/10; junc_x=(5−d)/5 only for d≤4; z clipped to ±50 then sigmoid.
- EVID assembly: primary code[,JDd][,HPhp][,GNOMAD][,MS].

## Output (10-column contract)

`chrom  pos0  pos0+1  strand(+/-/amb)  evid  call_score(%.4f)  depth  var_freq(%.4f)  fwd_freq(%.4f)  rev_freq(%.4f)`. `var_freq` = any-mismatch frequency (`var_reads/depth`). `fwd_freq` = editing-consistent frequency on the forward strand (A>G), `rev_freq` = editing-consistent frequency on the reverse strand (T>C); both are zero for a site whose reference base cannot be A-to-I edited (REF=C or REF=G). Sorted by (chrom, pos0) before writing → thread-count independent, byte-identical. An empty contig yields zero candidates.

## Dependencies

rust-htslib (BAM path + faidx), thiserror, serde/serde_json (spec), rayon, flate2 (annotation gz), esperanto-bamio (canonical `.baln` reader, shared with score/pileup).

## Self-checks

- A negative-strand read carrying a stored variant → site goes into the rev bucket, var_rev>0, var_freq correct (pass-through; bases are not complemented).
- Reads missing EK are treated as dirty; clean reads with mm+ea=0 are not compared base-by-base.
- Without fasta: majority pseudo-ref counting does not crash; hp=0; junction does not vote.
- Same input twice, threads 1 vs 8: candidates.bed byte-identical.
- BAM and .baln dual sources: same data yields byte-identical output.
