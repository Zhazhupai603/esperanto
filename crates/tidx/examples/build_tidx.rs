//! One-off: build a `.tidx` transcript k-mer index from GTF + FASTA.
//! Usage: cargo run -p esperanto-tidx --example build_tidx -- <gtf> <fasta> <out.tidx>
use std::path::Path;
use esperanto_tidx::{build, BuildOptions};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() != 4 {
        eprintln!("usage: build_tidx <gtf> <fasta> <out.tidx>");
        std::process::exit(2);
    }
    let stats = build(Path::new(&a[1]), Path::new(&a[2]), Path::new(&a[3]), &BuildOptions::default())
        .expect("tidx build");
    eprintln!(
        "[tidx] {} transcripts, {} bp, {} entries, {} bytes, {:.1}s -> {}",
        stats.tx_count, stats.total_bp, stats.total_entries, stats.file_size, stats.build_seconds, a[3]
    );
}
