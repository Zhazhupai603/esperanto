//! `esperanto verify` — candidate-supporting read re-alignment (→
//! `esperanto_verify::run`). Keeps the vote of a supporting read only when
//! its site placement survives exhaustive re-alignment.

use std::path::PathBuf;

use clap::Args;
use esperanto_verify::VerifyParams;

#[derive(Args)]
pub struct VerifyArgs {
    /// Input BAM (coordinate-sorted, indexed).
    #[arg(long)]
    bam: PathBuf,
    /// Candidate sites (scan's candidates.bed, 11 columns; verify's own
    /// verified beds, 13, also accepted).
    #[arg(long)]
    candidates: PathBuf,
    /// Alignment index (default: refs discovery).
    #[arg(long)]
    index: Option<PathBuf>,
    /// L1 engine bundle (default: discovered next to the index).
    #[arg(long)]
    l1_bundle: Option<PathBuf>,
    /// Output verified bed (input columns + alt_verified + withdrawn).
    #[arg(long)]
    out: PathBuf,
    /// Worker threads (0 = all cores).
    #[arg(long, default_value_t = 0)]
    threads: usize,
    /// "Obviously best" margin (best - second strict score) that keeps a
    /// vote. Default 12; the WP3 calibration plan sanctions 6 as the first
    /// relaxation step.
    #[arg(long)]
    ob_margin: Option<i32>,
    /// Mismatch channel of the candidates (edit | ctrl1 | ctrl2).
    #[arg(long, default_value = "edit")]
    channel: String,
}

pub fn run(a: VerifyArgs) -> anyhow::Result<()> {
    let mut index = a.index;
    if index.is_none() {
        if let Some(refs) = crate::resolve::refs() {
            refs.fill_index(&mut index);
        }
    }
    let index = index.ok_or_else(|| {
        anyhow::anyhow!("--index not given and no refs directory found (set ESPERANTO_REFS)")
    })?;
    let l1 = crate::resolve::l1_bundle(&a.l1_bundle, &index);
    let channel = match a.channel.as_str() {
        "edit" => esperanto_verify::Channel::Edit,
        "ctrl1" => esperanto_verify::Channel::Ctrl1,
        "ctrl2" => esperanto_verify::Channel::Ctrl2,
        other => anyhow::bail!("--channel must be edit|ctrl1|ctrl2, got '{other}'"),
    };
    let vp = VerifyParams {
        bam: a.bam,
        candidates: a.candidates,
        index,
        l1_bundle: l1,
        out: a.out,
        threads: crate::resolve::threads(a.threads),
        ob_margin: a.ob_margin.unwrap_or(12),
        channel,
        ..Default::default()
    };
    let stats = esperanto_verify::run(&vp)?;
    eprintln!(
        "[verify] {} candidates: {} kept (>= {} verified votes), {} reads re-aligned, {} votes withdrawn (no_site={} [empty={} other_locus={} refskip={} del={} outside={} base_mm={}], margin={}, blat={}, pair_anchored_kept={})",
        stats.candidates,
        stats.sites_kept,
        esperanto_verify::VerifyParams::default().min_verified,
        stats.reads_examined,
        stats.votes_withdrawn,
        stats.wd_no_site,
        stats.ns_empty,
        stats.ns_other_locus,
        stats.ns_refskip,
        stats.ns_del,
        stats.ns_outside_span,
        stats.ns_base_mismatch,
        stats.wd_margin,
        stats.wd_blat,
        stats.pair_anchored_kept
    );
    Ok(())
}
