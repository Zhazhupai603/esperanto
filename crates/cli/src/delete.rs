//! `esperanto delete` — remove local reference/index data (refs directory)
//! and/or the model bundle. Destructive: requires typing `yes` on a TTY
//! (or `--yes` for scripts); refuses otherwise.

use std::path::PathBuf;

use clap::{Args, ValueEnum};

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum DeleteTarget {
    /// The refs directory (reference FASTA/GTF + all index artifacts).
    Refs,
    /// The scoring model bundle.
    Bundle,
    /// Both of the above.
    All,
}

#[derive(Args)]
pub struct DeleteArgs {
    /// What to delete: refs | bundle | all.
    #[arg(value_enum)]
    target: DeleteTarget,
    /// Skip the confirmation prompt (scripts only).
    #[arg(long)]
    yes: bool,
}

fn dir_summary(p: &std::path::Path) -> (u64, u64) {
    let mut files = 0u64;
    let mut bytes = 0u64;
    let mut stack = vec![p.to_path_buf()];
    while let Some(d) = stack.pop() {
        if let Ok(rd) = std::fs::read_dir(&d) {
            for e in rd.flatten() {
                let p = e.path();
                if p.is_dir() {
                    stack.push(p);
                } else if p.is_file() {
                    files += 1;
                    bytes += e.metadata().map(|m| m.len()).unwrap_or(0);
                }
            }
        }
    }
    (files, bytes)
}

fn human_bytes(n: u64) -> String {
    if n >= 1 << 30 {
        format!("{:.1} GB", n as f64 / (1u64 << 30) as f64)
    } else if n >= 1 << 20 {
        format!("{:.1} MB", n as f64 / (1u64 << 20) as f64)
    } else {
        format!("{n} B")
    }
}

fn confirm(prompt_paths: &[(String, u64, u64)], yes: bool) -> anyhow::Result<bool> {
    for (path, files, bytes) in prompt_paths {
        eprintln!("  {path}  ({files} files, {})", human_bytes(*bytes));
    }
    if !yes && !std::io::IsTerminal::is_terminal(&std::io::stdin()) {
        eprintln!("refusing to delete without a TTY; pass --yes to confirm");
        return Ok(false);
    }
    crate::confirm::confirmed(
        "Delete the above? This cannot be undone.",
        ("No, keep everything", "Yes, delete permanently"),
        yes,
    )
}

pub fn run(a: DeleteArgs) -> anyhow::Result<()> {
    let mut targets: Vec<PathBuf> = Vec::new();
    match a.target {
        DeleteTarget::Refs => targets.push(crate::resolve::home_refs_dir()),
        DeleteTarget::Bundle => targets.push(crate::resolve::home_bundle_dir()),
        DeleteTarget::All => {
            targets.push(crate::resolve::home_refs_dir());
            targets.push(crate::resolve::home_bundle_dir());
        }
    }
    let existing: Vec<PathBuf> = targets.into_iter().filter(|p| p.exists()).collect();
    if existing.is_empty() {
        eprintln!("nothing to delete");
        return Ok(());
    }
    let summary: Vec<(String, u64, u64)> = existing
        .iter()
        .map(|p| {
            let (f, b) = dir_summary(p);
            (p.display().to_string(), f, b)
        })
        .collect();
    eprintln!("The following will be deleted permanently:");
    if !confirm(&summary, a.yes)? {
        eprintln!("aborted; nothing deleted");
        return Ok(());
    }
    for p in &existing {
        std::fs::remove_dir_all(p)?;
        eprintln!("deleted {}", p.display());
    }
    Ok(())
}
