//! Parity check for the CMC scoring network against the exported golden
//! vectors (PyTorch reference, eval mode, f32, CPU).
//!
//! Usage: `cargo run --release -p esperanto-score --example cmc_parity`
//! Reads the installed v1.7.0 bundle and the golden JSON produced by the
//! training-side exporter; fails on any probability outside tolerance.

use std::path::PathBuf;

use esperanto_score::cmcnet::CmcNet;

fn main() -> anyhow::Result<()> {
    let bundle: PathBuf = std::env::var("ESPERANTO_BUNDLE")
        .map(PathBuf::from)
        .unwrap_or_else(|_| {
            std::env::var("HOME")
                .map(|h| {
                    PathBuf::from(h)
                        .join(".local/share/esperanto/bundle/human/esperanto-model-v1.7.0/rust")
                })
                .expect("HOME")
        });
    let golden_path = std::env::var("ESP_GOLDEN")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("run/export_v17/golden/golden.json"));

    let golden: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&golden_path)?)?;
    let inputs = golden["inputs"]
        .as_array()
        .expect("inputs array")
        .iter()
        .map(|site| {
            let mut x = [0.0f32; 800];
            for ch in 0..16usize {
                for col in 0..50usize {
                    x[ch * 50 + col] = site[ch][col].as_f64().expect("f32") as f32;
                }
            }
            x
        })
        .collect::<Vec<_>>();

    let atol: f64 = std::env::var("ESP_ATOL")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(5e-4);

    let mut worst_overall = 0.0f64;
    for task in ["main", "hyper"] {
        let bytes = std::fs::read(bundle.join(format!("net_{task}.safetensors")))?;
        let st = safetensors::SafeTensors::deserialize(&bytes)?;
        let net = CmcNet::load(&st)?;
        let want = golden["nets"][task]
            .as_array()
            .expect("reference probabilities")
            .iter()
            .map(|p| p[1].as_f64().expect("p_re"))
            .collect::<Vec<_>>();
        let mut worst = 0.0f64;
        let mut worst_i = 0usize;
        for (i, x) in inputs.iter().enumerate() {
            let got = net.re_prob(x);
            let d = (got - want[i]).abs();
            if d > worst {
                worst = d;
                worst_i = i;
            }
        }
        println!(
            "[{task}] sites={} worst_abs_diff={worst:.3e} (row {worst_i}) atol={atol:.1e}",
            inputs.len()
        );
        worst_overall = worst_overall.max(worst);
        if worst > atol {
            anyhow::bail!("{task} parity failed: worst {worst:.3e} > atol {atol:.1e}");
        }
    }
    println!("cmc_parity PASS (worst {worst_overall:.3e})");
    Ok(())
}
