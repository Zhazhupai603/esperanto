# score — RE_PROB model forward pass (CMC CNN + Transformer, bundle v1.7)

BAM + sites -> RE_PROB (chrom pos prob; row order = input sites row order). Model and scoring semantics are frozen; the gold-standard metrics (AUROC/AUPRC) must not regress.

## bundle (zero-config resolution)

- Layout: bundle root (`feature_spec.json`, `model_card.json`, `manifest.sha256`) + `rust/` (`net_main.safetensors`, `net_hyper.safetensors`, `threshold.json`). A directory with `rust/net_main.safetensors` is a CMC bundle (`is_cmc_bundle`).
- `load_cmc_bundle(rust_dir)`: both networks (shape-validated safetensors loads; every tensor is checked against the frozen architecture) + the frozen thresholds. `load_cmc_thresholds` reads the calibration without loading weights.
- `threshold.json`: `{"main": {"threshold": T_main}, "hyper": {"threshold": T_hyper}}` — validation-donor calibrated (main: DM false-positive rate <= 1%; hyper: ordinary-editing false-positive rate <= 5%); the VCF FILTER criteria follow these numbers.

## Feature flow (per site)

1. **CMC 16x50 tensor** (`cmc.rs`, schema `cmc-cigar-window-1`): window = `[center-25, center+26)` (0-based, half-open), `center = pos-1`; reads with `flag & 3844 == 0` and `MAPQ >= 20`; per-position base quality `>= 13`, base in `ACTG` (a 0xFF quality byte counts as 0); CIGAR walk (M/=/X consume both tracks, I/S the query, D/N the reference, H/P neither, stop once the reference cursor passes the window end); template names strip a trailing `/1`//`/2`; within a template the highest-quality base at each reference position wins, equal-quality disagreements store `N`; a read tagged `RE:Z:collapsed` that yields the center base marks its template (QC only).
2. Tensor row = `4*flank_index + center_index` over `ACTG` order; column `c` = offset `-25..-1,+1..+25` (center column excluded); a covered column sums to one (integer counts divided in f32 — bit-identical to the reference extractor), an uncovered column is all zeros.
3. **network** (`cmcnet.rs`): conv1d(16->64,k5,same) -> batchnorm -> gelu -> conv1d(64->128,k3,same) -> batchnorm -> gelu -> transpose to (50,128) -> token embedding (learned cls + position, d=128) -> transformer encoder x4 (8 heads, ff 512, pre-norm, final LayerNorm) -> cls pooling -> linear(128) -> gelu -> linear(2) -> softmax[1] (f64 accumulation). GELU uses the exact erf form evaluated in f64; all storage is f32. 847,042 parameters, identical architecture for both networks.
4. Output: `RE_PROB` = p(positive class). Main network classes: 0 = DNA mutation (WGS-confirmed exact allele carrier), 1 = RNA editing candidate. Hyper network classes: 0 = ordinary editing candidate, 1 = rescue-supported editing candidate.

## Pipeline (determinism contract)

- **Species guardrail**: the FASTA index's chr1 length must be 248956422 (hg38) or < 10 Mb (synthetic/test), otherwise refused; `ReferenceCheck::TrustedHybrid` (flow only, after the species.json manifest verified the hybrid reference) skips the chr1 heuristic. The model reads no reference sequence — the FASTA is used for this check only.
- Sites are sorted by (chrom, pos); CMC extraction runs as region sweeps (one fetch per run of sites within `CMC_MERGE_GAP` = 2000, group cap 16k, one IndexedReader per rayon worker); the network forward runs across sites on the same pool. Results are written back by original index — **thread count does not change values** (every site is scored independently).
- `score_sites_batched` (main network, ordinary candidates) and `score_hyper_sites_batched` (hyper network, Tier-2-only sites) share the sweep + forward path.
- Output: TSV `chrom\tpos\tprob` (f64 Display); row order = input order.

## Parity contracts

- Rust CMC extraction is bit-identical to the reference Python extractor on the same BAM (integer counts + single f32 division).
- The Rust forward matches the exported PyTorch reference on golden vectors (fixed inputs -> reference probabilities) within 5e-4 absolute; measured worst-case on release builds: 4.7e-7. The `cmc_parity` example checks the installed bundle against the golden file; `cmc_extract_parity` dumps raw tensors for extractor comparisons.

## Dependencies

Workspace-level: ndarray 0.16, libm 0.2, safetensors 0.4; existing: rust-htslib, rayon, thiserror, serde/serde_json, anyhow.
