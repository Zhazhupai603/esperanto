# score — RE_PROB model forward pass (v1.4.1 bundle, mamba ILP + pileup veto gate)

BAM + sites → RE_PROB TSV (chrom pos prob; row order = input sites row order). Semantic source: a wholesale port of the legacy esperanto-score scoring pipeline; model and scoring semantics frozen; the gold-standard metrics (AUROC/AUPRC) must not regress.

## bundle (v1.4.1-501_40ep, zero-config resolution)

- Layout: bundle root (encoder/model.safetensors + config.json) + `rust/` (heads/fold_{0..4}.safetensors, norm.json, gate_heads/, gate_norm.json, feature_spec.json).
- `load_bundle(rust_dir)`: 5-fold FoldHead + per-fold NormStats (8-dim mean/std); `half_window` read from feature_spec.sequence.half_window (v1.4.1 = 250 → 501 bp window; looked up in both dir and dir.parent, missing = 500 for compatibility); `cache_id` = hash of the full feature_spec text + half_window.
- Encoder resolution: when `--caduceus` is omitted, probe bundle/encoder then the bundle parent's encoder (model.safetensors required), in that order.
- Veto-gate trio (gate_heads/ + gate_norm.json + feature_spec.gate.threshold=0.004817): all missing = None; partially missing = hard error for a corrupt bundle. The score pipeline requires them to be present.

## Feature flow (per site)

1. **pileup 8-dim** (esperanto-pile batch interface; rtol=0 contract, see pile spec): `[depth, A,C,G,T count, mean_bq, strand_bias, mean_mapq]`.
2. **Veto gate**: gate RE_PROB = 5-fold zero-embedding ensemble (same head form, emb all zeros); `< threshold` → final probability = gate probability, **encoder skipped**; otherwise the site is kept.
3. **Sequence window**: center = pos (1-based), ±half_window; N-pad at contig edges; uppercase; if the fasta lacks the contig → all-N window (no error).
4. **tokenizer**: A7 C8 G9 T10 N11, others → 6 [UNK], trailing SEP=1, no BOS (total length 2·half+2 = 502).
5. **encoder**: CaduceusEncoder (Mamba 4L d118, BiMamba RC-equiv, faer GEMM single-threaded Par::Seq; numerical determinism guarded by sub-batch slice comparison), mean-pool → f16[118].
6. **head**: per-fold z-score (pileup) → MLP (8→32→32, ReLU) → concat(emb118, 32)=150 → MLP (150→128→2) → softmax[1] (fp64 accumulation); ensemble = mean over 5 folds.

## Pipeline (determinism contract)

- Reference preload: contigs involved in sites are loaded into memory once (uppercase); **species guardrail**: chr1 length in .fai must be == 248956422 (hg38) or < 10 Mb (synthetic/test reference), otherwise the run is refused.
- Sites are sorted by (chrom, pos); pileup features for all sites come from one region-sweep pass (MERGE_GAP groups, one IndexedReader per worker), then sites are cut into batches (default 256), batch-level rayon (0 = all cores); within a batch: gate → kept sites embedded in EMBED_SLICE=128 sub-batches → fused head; results written back by original index — **thread count and batch size do not change values**.
- Output: TSV `chrom\tpos\tprob` (f64 Display); row order = input order.

## Embedding cache (optional optimization)

ESPEMBC2 format (magic + cache_id + half_window + append records); incompatible → warn and reopen as an empty table; a hit returns the fp16 bit pattern (bit-identical to online embed); flush rewrites atomically (tmp+rename). Pure optimization, no semantics.

## Out of scope (1.0.0)

- corrected_sites.tsv feature source (the EM closed loop depends on the realign crate): not ported; registered in BACKLOG.
- report stage: deferred to 1.1 (design finalized).

## Dependencies

New workspace-level: ndarray 0.16, half 2.4, safetensors 0.4, faer 0.24.4; existing: rust-htslib, rayon, thiserror, serde/serde_json, anyhow, esperanto-pile (path).

## Self-checks

- Sites with gate probability < threshold: output = gate probability, encoder untouched (provable via timing/cache counters).
- z-score/softmax/ensemble match the torch reference with pointwise tolerance at the 1e-6 level (gold-standard probes).
- Same input twice, threads 1 vs 8, batch 64 vs 1024: output TSV byte-identical.
- Cold vs warm cache runs: output TSV byte-identical.
