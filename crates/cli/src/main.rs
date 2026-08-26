//! esperanto — CLI for the ESPERANTO RNA-editing pipeline.
//!
//! Spec: docs/specs/crates/cli.md (single source of truth). Thin dispatch:
//! parse args → zero-config resolution → crate entry points. No scientific
//! semantics live here.

mod index;
mod map;
mod pile;
mod qc;
mod report;
mod resolve;
mod run;
mod scan;
mod score;
mod setup;

use clap::{Parser, Subcommand};

/// ESPERANTO RNA editing analysis (1.0.0).
#[derive(Parser)]
#[command(name = "esperanto", version, about, propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// FASTQ quality control (trim/filter/report).
    Qc(qc::QcArgs),
/// Build a paidx alignment index from a reference FASTA.
Index(index::IndexArgs),
    /// Read alignment (RNA 2-pass; unsorted raw.bam — sort happens in `run`).
    Map(map::MapArgs),
    /// Pileup feature extraction for single sites or a site list.
    Pile(pile::PileArgs),
    /// Candidate editing-site discovery from BAM/.baln.
    Scan(scan::ScanArgs),
    /// RE_PROB scoring for a site list.
    Score(score::ScoreArgs),
    /// Full pipeline: qc → map → sort → scan → score → vcf.
    Run(run::RunArgs),
    /// Regenerate the standalone HTML report for a finished run directory.
    Report(report::ReportArgs),
    /// One-step reference environment: detect or download reference files
    /// in the refs directory, then build the index in place.
    Setup(setup::SetupArgs),
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Qc(a) => qc::run(a),
        Cmd::Index(a) => index::run(a),
        Cmd::Map(a) => map::run(a),
        Cmd::Pile(a) => pile::run(a),
        Cmd::Scan(a) => scan::run(a),
        Cmd::Score(a) => score::run(a),
        Cmd::Run(a) => run::run(a),
        Cmd::Report(a) => report::run(a),
        Cmd::Setup(a) => setup::run(a),
    }
}
