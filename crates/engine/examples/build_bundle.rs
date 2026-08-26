//! One-off: fuse an existing `.tidx` + GTF + FASTA into an L1 bundle
//! (`L1Index::build` + `save`). Usage:
//!   cargo run -p esperanto-engine --example build_bundle -- <tidx> <gtf> <fasta> <out.bndl>
use std::path::Path;
use esperanto_engine::{L1Index, Tidx};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() != 5 {
        eprintln!("usage: build_bundle <tidx> <gtf> <fasta> <out.bndl>");
        std::process::exit(2);
    }
    let idx = L1Index::build(Path::new(&a[1]), Path::new(&a[2]), Path::new(&a[3]))
        .expect("L1Index::build");
    idx.save(Path::new(&a[4])).expect("save");
    eprintln!("[l1] built {} transcripts -> {}", idx.tx_count(), a[4]);
}
