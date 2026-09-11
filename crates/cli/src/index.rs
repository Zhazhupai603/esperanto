//! `esperanto index` — build a paidx alignment index from a reference FASTA.

use std::path::PathBuf;

use clap::Args;
use esperanto_engine::Tidx;
use esperanto_map::seed::SeedParams;
use esperanto_map::{fasta, index, index_io};

#[derive(Args)]
pub struct IndexArgs {
    /// Reference FASTA.
    #[arg(long)]
    fasta: PathBuf,
    /// Optional transcript annotation GTF. When given, the L1 engine bundle
    /// (`<out stem>.bndl`) is built alongside the paidx index.
    #[arg(long)]
    gtf: Option<PathBuf>,
    /// Output paidx index path.
    #[arg(long)]
    out: PathBuf,
    /// K-mer length (must match the value used at alignment time, default 15).
    #[arg(long, default_value_t = 15)]
    k: u32,
    /// Window size in k-mers.
    #[arg(long, default_value_t = 5)]
    w: u32,
    /// Decoy contig source FASTA (repeatable). Contigs are renamed with the
    /// `decoy_` prefix and indexed alongside the main reference; the
    /// three-letter rescue indexes and the L1 bundle exclude them.
    #[arg(long = "decoy")]
    decoys: Vec<PathBuf>,
    /// Also build the negative-control rescue alphabet indexes
    /// (`<out stem>.ac.paidx` A->C, `<out stem>.tg.paidx` T->G, k=31/w=10).
    #[arg(long)]
    rescue_control: bool,
    /// Build only the rescue alphabet indexes (all four), leaving an
    /// existing paidx untouched.
    #[arg(long)]
    rescue_only: bool,
}

pub fn run(a: IndexArgs) -> anyhow::Result<()> {
    if a.rescue_only {
        // Rescue indexes only (the paidx already exists; decoys never enter
        // the rescue alphabets, so the plain FASTA is the right source).
        let bytes = std::fs::read(&a.fasta).map_err(|e| anyhow::anyhow!("read fasta: {e}"))?;
        return build_rescue_indexes(&bytes, &a.out, true);
    }
    build_all(&a.fasta, a.gtf.as_deref(), &a.decoys, &a.out, a.k, a.w, a.rescue_control)
}

/// Concatenate decoy FASTAs onto the main one, forcing the `decoy_` contig
/// name prefix (the contract every downstream stage keys off).
fn with_decoys(main: &[u8], decoys: &[PathBuf]) -> anyhow::Result<Vec<u8>> {
    let mut out = main.to_vec();
    for d in decoys {
        let db = std::fs::read(d)
            .map_err(|e| anyhow::anyhow!("read decoy fasta {}: {e}", d.display()))?;
        for line in db.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            if line.starts_with(b">") {
                let header = &line[1..];
                let name_end = header
                    .iter()
                    .position(|b| b.is_ascii_whitespace())
                    .unwrap_or(header.len());
                let name = &header[..name_end];
                out.push(b'>');
                if !name.starts_with(b"decoy_") {
                    out.extend_from_slice(b"decoy_");
                }
                out.extend_from_slice(header);
            } else {
                out.extend_from_slice(line);
            }
            out.push(b'\n');
        }
    }
    Ok(out)
}

/// Rewrite `from` to `to` (both cases) in sequence lines; headers pass
/// through. Under the same transform applied to a read, an edit of that
/// direction becomes an identity, so a hyperedited read seeds normally.
///
/// One direction per index, never both: merging A/G *and* T/C at once leaves a
/// two-letter alphabet, in which a 31-mer carries 31 bits (~2×10⁹ distinct
/// keys) against a 3.1×10⁹-base genome — seeds stop being locators. A
/// three-letter 31-mer carries ~49 bits, which is ample. It is also the
/// specificity argument: A>G and T>C are the two transitions, so exempting
/// both makes a paralog indistinguishable from an edited sibling (Porath,
/// Carmi & Levanon, Nat Commun 5:4726, 2014, use exactly one transform per
/// index and treat the reverse transforms as negative controls).
fn transform_fasta_bytes(bytes: &[u8], from: u8, to: u8) -> Vec<u8> {
    let lo_from = from.to_ascii_lowercase();
    let lo_to = to.to_ascii_lowercase();
    let mut out = Vec::with_capacity(bytes.len());
    for line in bytes.split(|&b| b == b'\n') {
        if line.starts_with(b">") {
            out.extend_from_slice(line);
        } else {
            out.extend(line.iter().map(|&b| {
                if b == from {
                    to
                } else if b == lo_from {
                    lo_to
                } else {
                    b
                }
            }));
        }
        out.push(b'\n');
    }
    out
}

/// Build the paidx index and — when a GTF is given — the L1 engine bundle
/// alongside it. Always also builds the two three-letter hyperediting rescue
/// indexes (`<out stem>.ag.paidx` with A→G, `<out stem>.tc.paidx` with T→C,
/// both k=31/w=10) used by the map-stage rescue pass. Shared by `index` and
/// `setup`.
pub(crate) fn build_all(
    fasta: &std::path::Path,
    gtf: Option<&std::path::Path>,
    decoys: &[PathBuf],
    out: &std::path::Path,
    k: u32,
    w: u32,
    rescue_control: bool,
) -> anyhow::Result<()> {
    let bytes = std::fs::read(fasta).map_err(|e| anyhow::anyhow!("read fasta: {e}"))?;
    // Decoys join the paidx only: the rescue indexes below transform the
    // main reference bytes and the L1 bundle is built from the GTF alone.
    let combined = if decoys.is_empty() {
        bytes.clone()
    } else {
        with_decoys(&bytes, decoys)?
    };
    let reference = fasta::parse_fasta_bytes(&combined)
        .map_err(|e| anyhow::anyhow!("parse fasta: {e}"))?;
    let idx = index::Index::build(reference, SeedParams { k, w });
    index_io::save(&idx, out).map_err(|e| anyhow::anyhow!("save paidx: {e}"))?;
    if !decoys.is_empty() {
        let n: usize = idx
            .reference
            .contigs
            .iter()
            .filter(|c| c.name.starts_with("decoy_"))
            .count();
        eprintln!("[index] wrote {} (+{n} decoy contigs)", out.display());
    } else {
        eprintln!("[index] wrote {}", out.display());
    }
    build_rescue_indexes(&bytes, out, rescue_control)?;
    if let Some(gtf) = gtf {
        // The L1 runtime needs BOTH the .bndl (projection + sequences) and its
        // .tidx sidecar (k-mer index) — keep all three artifacts together.
        let tidx_path = out.with_extension("tidx");
        esperanto_tidx::build(
            gtf,
            fasta,
            &tidx_path,
            &esperanto_tidx::BuildOptions::default(),
        )
        .map_err(|e| anyhow::anyhow!("build tidx: {e}"))?;
        let bndl = out.with_extension("bndl");
        let l1 = esperanto_engine::L1Index::build(&tidx_path, gtf, fasta)
            .map_err(|e| anyhow::anyhow!("build L1 bundle: {e}"))?;
        l1.save(&bndl)
            .map_err(|e| anyhow::anyhow!("save bndl: {e}"))?;
        eprintln!(
            "[index] wrote {} + {} ({} transcripts)",
            bndl.display(),
            tidx_path.display(),
            l1.tx_count()
        );
    }
    Ok(())
}

/// The three-letter rescue indexes (`ag`/`tc`, plus `ac`/`tg` when
/// `rescue_control`): the main reference bytes transformed per alphabet,
/// k=31/w=10, decoys never included.
fn build_rescue_indexes(
    bytes: &[u8],
    out: &std::path::Path,
    rescue_control: bool,
) -> anyhow::Result<()> {
    let mut rescue: Vec<(u8, u8, &str)> = vec![(b'A', b'G', "ag.paidx"), (b'T', b'C', "tc.paidx")];
    if rescue_control {
        // Negative-control alphabets: transforms editing cannot produce.
        rescue.extend([(b'A', b'C', "ac.paidx"), (b'T', b'G', "tg.paidx")]);
    }
    for (from, to, ext) in rescue {
        let tref = fasta::parse_fasta_bytes(&transform_fasta_bytes(bytes, from, to))
            .map_err(|e| anyhow::anyhow!("parse {ext} fasta: {e}"))?;
        let tidx = index::Index::build(tref, SeedParams { k: 31, w: 10 });
        let tpath = out.with_extension(ext);
        index_io::save(&tidx, &tpath).map_err(|e| anyhow::anyhow!("save {ext}: {e}"))?;
        eprintln!("[index] wrote {}", tpath.display());
    }
    Ok(())
}
