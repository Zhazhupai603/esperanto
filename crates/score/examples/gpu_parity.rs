//! GPU-vs-CPU parity check for the Caduceus score encoder (feature `gpu`).
//!
//! Usage: gpu_parity <sorted.bam> <sites.vcf> <ref.fa>
//!
//! Loads the installed model bundle + encoder, builds 64 site windows from the
//! given BAM + VCF through the production window/tokenize path, then:
//!   1. compares CPU `embed_batch` vs GPU `embed_batch` embeddings (fp16 -> fp32),
//!   2. compares RE_PROB from the shared 5-fold head on both embeddings,
//!   3. times both paths on the same 64-site batch (median of 10).
//!
//! Targets: embedding max abs diff < 1e-2, RE_PROB max abs diff < 1e-3, GPU faster per batch.

use anyhow::{anyhow, Context, Result};
use esperanto_score::bundle::load_bundle;
use esperanto_score::caduceus::gpu_encoder::GpuCaduceusEncoder;
use esperanto_score::caduceus::{CaduceusEncoder, EmbedBatchBufs, D_MODEL};
use esperanto_score::encoder::{fetch_window_mem_hw, tokenize};
use esperanto_score::head::re_prob_ensemble;
use esperanto_score::pipeline::resolve_encoder_from_bundle;
use std::path::PathBuf;
use std::time::Instant;

const N_SITES: usize = 64;

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn main() -> Result<()> {
    let home = std::env::var("HOME").context("HOME")?;
    let bundle_dir = PathBuf::from(home)
        .join(".local/share/esperanto/bundle/human/esperanto-model-v1.6.0/rust");
    let mut args = std::env::args().skip(1);
    let (bam_path, vcf_path, fasta_path) = match (args.next(), args.next(), args.next()) {
        (Some(b), Some(v), Some(f)) => (PathBuf::from(b), PathBuf::from(v), PathBuf::from(f)),
        _ => return Err(anyhow!("usage: gpu_parity <sorted.bam> <sites.vcf> <ref.fa>")),
    };

    // ---- bundle + encoders (existing loaders) ----
    let bundle = load_bundle(&bundle_dir).context("load bundle")?;
    let enc_dir = resolve_encoder_from_bundle(&bundle_dir)?;
    let cpu = CaduceusEncoder::load(&enc_dir).context("load cpu encoder")?;
    let dev = candle_core::Device::new_cuda(0).context("cuda device 0")?;
    let gpu = GpuCaduceusEncoder::from_cpu(&cpu, &dev).context("upload encoder to gpu")?;

    // ---- 64 real sites from the vcf (first N_SITES data rows) ----
    let vcf = std::fs::read_to_string(&vcf_path).context("read sites.vcf")?;
    let mut sites: Vec<(String, i64)> = Vec::with_capacity(N_SITES);
    for line in vcf.lines() {
        if line.starts_with('#') {
            continue;
        }
        let mut cols = line.split('\t');
        let chrom = cols.next().context("chrom column")?;
        let pos: i64 = cols
            .next()
            .context("pos column")?
            .parse()
            .with_context(|| format!("parse pos in: {line}"))?;
        sites.push((chrom.to_string(), pos));
        if sites.len() == N_SITES {
            break;
        }
    }
    if sites.len() != N_SITES {
        return Err(anyhow!("only {} sites parsed from vcf", sites.len()));
    }

    // ---- reference preload (score-pipeline style) ----
    let fa = rust_htslib::faidx::Reader::from_path(&fasta_path).context("open fasta")?;
    let mut refmap: std::collections::HashMap<String, Vec<u8>> = std::collections::HashMap::new();
    for (chrom, _) in &sites {
        if refmap.contains_key(chrom) {
            continue;
        }
        let len = fa.fetch_seq_len(chrom) as usize;
        let seq = fa
            .fetch_seq(chrom, 0, len - 1)
            .with_context(|| format!("fetch {chrom}"))?;
        refmap.insert(
            chrom.clone(),
            seq.iter().map(|b| b.to_ascii_uppercase()).collect(),
        );
    }

    // ---- tokens via the production window/tokenize path ----
    let hw = bundle.half_window;
    let tok_len = (2 * hw + 2) as usize;
    let mut tokens_flat: Vec<i64> = Vec::with_capacity(N_SITES * tok_len);
    for (chrom, pos) in &sites {
        let seq = refmap.get(chrom).context("contig missing")?;
        let window = fetch_window_mem_hw(seq, *pos, hw);
        tokens_flat.extend(tokenize(&window));
    }

    // ---- pileup features for the shared head (production pile path) ----
    let mut bam = rust_htslib::bam::IndexedReader::from_path(&bam_path).context("open bam")?;
    let refs: Vec<(&str, i64)> = sites.iter().map(|(c, p)| (c.as_str(), *p)).collect();
    let piles = esperanto_pile::extract_pileup_features_batch(&mut bam, &refs)
        .map_err(|e| anyhow!("pileup: {e}"))
        .context("pileup batch")?;

    // ---- embeddings: CPU vs GPU ----
    let mut bufs = EmbedBatchBufs::default();
    let cpu_emb = cpu
        .embed_batch(&tokens_flat, N_SITES, tok_len, &mut bufs)
        .context("cpu embed_batch")?;
    let gpu_emb = gpu
        .embed_batch(&tokens_flat, N_SITES, tok_len)
        .context("gpu embed_batch")?;

    let (mut emb_max, mut emb_sum) = (0f32, 0f32);
    for (ce, ge) in cpu_emb.iter().zip(gpu_emb.iter()) {
        for (cv, gv) in ce.iter().zip(ge.iter()) {
            let (cv, gv) = (cv.to_f32(), gv.to_f32());
            let d = (cv - gv).abs();
            emb_max = emb_max.max(d);
            emb_sum += d;
        }
    }
    let emb_mean = emb_sum / (N_SITES * D_MODEL) as f32;

    // ---- RE_PROB from the shared 5-fold head, on both embeddings ----
    let mut re_max = 0f64;
    let mut re_cpu_sum = 0f64;
    for (i, (ce, ge)) in cpu_emb.iter().zip(gpu_emb.iter()).enumerate() {
        let cpu_v: Vec<f32> = ce.iter().map(|v| v.to_f32()).collect();
        let gpu_v: Vec<f32> = ge.iter().map(|v| v.to_f32()).collect();
        let cpu_v = ndarray::Array1::from_vec(cpu_v);
        let gpu_v = ndarray::Array1::from_vec(gpu_v);
        let pc = re_prob_ensemble(&bundle, &cpu_v.view(), &piles[i]).context("head cpu")?;
        let pg = re_prob_ensemble(&bundle, &gpu_v.view(), &piles[i]).context("head gpu")?;
        re_max = re_max.max((pc - pg).abs());
        re_cpu_sum += pc;
    }

    // ---- timing: median per-batch ms over 10 iterations (3 warmups each) ----
    const ITERS: usize = 10;
    const WARM: usize = 3;
    let mut cpu_ms = Vec::new();
    for it in 0..(WARM + ITERS) {
        let t0 = Instant::now();
        cpu.embed_batch(&tokens_flat, N_SITES, tok_len, &mut bufs)
            .context("cpu embed_batch timing")?;
        let ms = t0.elapsed().as_secs_f64() * 1e3;
        if it >= WARM {
            cpu_ms.push(ms);
        }
    }
    let mut gpu_ms = Vec::new();
    for it in 0..(WARM + ITERS) {
        let t0 = Instant::now();
        gpu.embed_batch(&tokens_flat, N_SITES, tok_len)
            .context("gpu embed_batch timing")?;
        dev.synchronize().context("sync")?;
        let ms = t0.elapsed().as_secs_f64() * 1e3;
        if it >= WARM {
            gpu_ms.push(ms);
        }
    }
    let (cpu_med, gpu_med) = (median(&mut cpu_ms), median(&mut gpu_ms));

    // ---- report ----
    println!("sites={N_SITES} tok_len={tok_len} half_window={hw}");
    println!("embedding abs diff: max={emb_max:.6} mean={emb_mean:.8} (target max < 1e-2)");
    println!("RE_PROB abs diff:   max={re_max:.8} (target < 1e-3)");
    println!(
        "RE_PROB mean (cpu head): {:.6}",
        re_cpu_sum / N_SITES as f64
    );
    println!(
        "per-batch median: cpu={cpu_med:.2} ms  gpu={gpu_med:.2} ms  speedup={:.2}x (target gpu < cpu)",
        cpu_med / gpu_med
    );
    let pass = emb_max < 1e-2 && re_max < 1e-3 && gpu_med < cpu_med;
    println!("RESULT: {}", if pass { "PASS" } else { "FAIL" });
    if !pass {
        std::process::exit(1);
    }
    Ok(())
}
