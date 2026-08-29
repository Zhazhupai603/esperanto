//! `esperanto setup` — one-step reference environment: detect (or download)
//! reference files in the refs directory, generate the `.fai` when missing,
//! and build the alignment index (paidx + L1 bundle) in place.

use std::path::{Path, PathBuf};

use clap::Args;

/// Default human reference FASTA (GENCODE GRCh38 primary assembly),
/// downloaded when the refs directory holds neither a FASTA nor a GTF.
const HG38_FA_URL: &str = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz";
/// Default transcript annotation (GENCODE v44 basic), same trigger.
const GENCODE_GTF_URL: &str =
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz";
/// Release tarball providing the model bundle when none is installed
/// (same pin as install.sh).
const RELEASE_TARBALL_URL: &str = "https://github.com/Zhazhupai603/esperanto/releases/download/v1.0.0/esperanto-1.0.0-linux-x86_64.tar.gz";

#[derive(Args)]
pub struct SetupArgs {
    /// Rebuild even when the index artifacts already exist.
    #[arg(long)]
    force: bool,
}

/// The fixed setup folder (`ESPERANTO_REFS` honored): shared resolver.
fn refs_dir() -> PathBuf {
    crate::resolve::home_refs_dir()
}

/// Decompress a `.gz` file next to itself; the archive is removed.
pub(crate) fn gunzip_in_place(gz: &Path) -> anyhow::Result<PathBuf> {
    let plain = gz.with_extension("");
    let input = std::fs::File::open(gz)?;
    let mut dec = flate2::read::GzDecoder::new(input);
    let mut out = std::io::BufWriter::new(std::fs::File::create(&plain)?);
    std::io::copy(&mut dec, &mut out)?;
    std::fs::remove_file(gz)?;
    Ok(plain)
}

/// Download `url` to `dest` (streaming; honors http_proxy/https_proxy).
fn download(url: &str, dest: &Path) -> anyhow::Result<()> {
    crate::fetch::file(url, dest)
}

/// Write the standard `<fasta>.fai` (name, length, offset, line bases,
/// line width) without depending on external tools.
pub(crate) fn write_fai(fasta: &Path) -> anyhow::Result<()> {
    let data = std::fs::read(fasta)?;
    let mut fai = String::new();
    let mut i = 0usize;
    while i < data.len() {
        if data[i] != b'>' {
            i += 1;
            continue;
        }
        let header_end = data[i..]
            .iter()
            .position(|&b| b == b'\n')
            .map(|p| i + p)
            .ok_or_else(|| anyhow::anyhow!("unterminated fasta header"))?;
        let name_end = data[i + 1..]
            .iter()
            .position(|b| b.is_ascii_whitespace())
            .map(|p| i + 1 + p)
            .unwrap_or(header_end);
        let name = String::from_utf8_lossy(&data[i + 1..name_end]).into_owned();
        let mut offset = header_end + 1;
        let start = offset;
        let mut linebases = 0usize;
        let mut linewidth = 0usize;
        let mut length = 0usize;
        while offset < data.len() && data[offset] != b'>' {
            let line_end = data[offset..]
                .iter()
                .position(|&b| b == b'\n')
                .map(|p| offset + p)
                .unwrap_or(data.len());
            let mut lb = line_end - offset;
            if lb > 0 && data[line_end - 1] == b'\r' {
                lb -= 1;
            }
            let lw = line_end - offset + 1;
            if linebases == 0 && lb > 0 {
                linebases = lb;
                linewidth = lw;
            }
            length += lb;
            offset = line_end + 1;
        }
        fai.push_str(&format!(
            "{name}\t{length}\t{start}\t{linebases}\t{linewidth}\n"
        ));
        i = offset;
    }
    std::fs::write(fasta.with_extension("fa.fai"), fai)?;
    Ok(())
}

pub fn run(a: SetupArgs) -> anyhow::Result<()> {
    let dir = refs_dir();
    std::fs::create_dir_all(&dir)?;
    eprintln!("[setup] refs directory: {}", dir.display());

    // --- decompress any gzipped reference files in place (archives removed) ---
    for gz in crate::resolve::find_all(&dir, &[".fa.gz", ".fasta.gz", ".gtf.gz"])? {
        eprintln!("[setup] decompressing {}", gz.display());
        gunzip_in_place(&gz)?;
    }

    // --- split FASTA files by species (names carry the species tag; an
    // unmarked file is accepted only as the single fasta in the directory) ---
    let mut human_fa: Option<PathBuf> = None;
    let mut mouse_fa: Option<PathBuf> = None;
    let mut untagged: Vec<PathBuf> = Vec::new();
    for p in crate::resolve::find_all(&dir, &[".fa", ".fasta"])? {
        let name = p
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        match crate::resolve::species_of(&name) {
            Some("human") => {
                if human_fa.replace(p).is_some() {
                    anyhow::bail!(
                        "multiple human FASTA files in {}; keep exactly one",
                        dir.display()
                    );
                }
            }
            Some("mouse") => {
                if mouse_fa.replace(p).is_some() {
                    anyhow::bail!(
                        "multiple mouse FASTA files in {}; keep exactly one",
                        dir.display()
                    );
                }
            }
            _ => untagged.push(p),
        }
    }
    if !untagged.is_empty() {
        if untagged.len() == 1 && human_fa.is_none() && mouse_fa.is_none() {
            human_fa = untagged.pop(); // legacy single-species directory
        } else {
            anyhow::bail!(
                "unrecognizable FASTA name(s) in {}: {}; use a species-tagged name (e.g. hg38.fa / grcm39.fa)",
                dir.display(),
                untagged
                    .iter()
                    .map(|p| p.display().to_string())
                    .collect::<Vec<_>>()
                    .join(", ")
            );
        }
    }
    let mut human_gtf: Option<PathBuf> = None;
    let mut mouse_gtf: Option<PathBuf> = None;
    for p in crate::resolve::find_all(&dir, &[".gtf"])? {
        let name = p
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        match crate::resolve::species_of(&name) {
            Some("human") => {
                if human_gtf.replace(p).is_some() {
                    anyhow::bail!(
                        "multiple human GTF files in {}; keep exactly one",
                        dir.display()
                    );
                }
            }
            Some("mouse") => {
                if mouse_gtf.replace(p).is_some() {
                    anyhow::bail!(
                        "multiple mouse GTF files in {}; keep exactly one",
                        dir.display()
                    );
                }
            }
            _ => {
                // Unmarked GTF attaches to the only species present.
                if mouse_fa.is_some() && human_fa.is_none() {
                    mouse_gtf = Some(p);
                } else if human_fa.is_some() && mouse_gtf.is_none() {
                    human_gtf = Some(p);
                } else {
                    anyhow::bail!(
                        "unrecognizable GTF name {} in {}; use a species-tagged name",
                        p.display(),
                        dir.display()
                    );
                }
            }
        }
    }

    // --- fetch whatever is missing (human default as today; the mouse
    // reference is always staged so `run --hybrid` can splice onto it) ---
    let mut fetched_human = false;
    if human_fa.is_none() {
        let fa_gz = dir.join("GRCh38.primary_assembly.genome.fa.gz");
        download(HG38_FA_URL, &fa_gz)?;
        human_fa = Some(gunzip_in_place(&fa_gz)?);
        fetched_human = true;
    }
    if fetched_human && human_gtf.is_none() {
        let gtf_gz = dir.join("gencode.v44.basic.annotation.gtf.gz");
        download(GENCODE_GTF_URL, &gtf_gz)?;
        human_gtf = Some(gunzip_in_place(&gtf_gz)?);
    }
    if human_gtf.is_none() {
        eprintln!(
            "[setup] no GTF in {}; building the genomic index only (add a GTF and re-run to enable the L1 engine)",
            dir.display()
        );
    }
    // The mouse reference is staged on demand: the first
    // `run --hybrid` fetches it into this directory automatically
    // (genomic-only without a mouse GTF; place one for annotations).

    // --- human reference: fai, validation, index (unchanged behavior) ---
    let fasta = human_fa.expect("human fasta set above");
    let fai = fasta.with_extension("fa.fai");
    if !fai.is_file() {
        eprintln!("[setup] writing {}", fai.display());
        write_fai(&fasta)?;
    }
    validate_pair(&fasta, human_gtf.as_deref(), None)?;

    // --- index build ---
    let stem = fasta
        .file_stem()
        .and_then(|s| s.to_str())
        .ok_or_else(|| anyhow::anyhow!("bad fasta name {}", fasta.display()))?;
    let paidx = dir.join(format!("{stem}.paidx"));
    if paidx.is_file() && !a.force {
        eprintln!(
            "[setup] {} exists; skipping (use --force to rebuild)",
            paidx.display()
        );
    } else {
        crate::index::build_all(&fasta, human_gtf.as_deref(), &paidx, 15, 5)?;
    }

    // --- user-placed mouse baseline: fai + validation only (hybrid
    // indexes are built on demand at run time) ---
    if let Some(mouse_fasta) = mouse_fa {
        let mfai = mouse_fasta.with_extension("fa.fai");
        if !mfai.is_file() {
            eprintln!("[setup] writing {}", mfai.display());
            write_fai(&mouse_fasta)?;
        }
        let baseline = if mouse_fasta
            .file_name()
            .is_some_and(|n| n.to_string_lossy().to_lowercase().contains("mm10"))
        {
            "mm10"
        } else {
            "grcm39"
        };
        validate_pair(
            &mouse_fasta,
            mouse_gtf.as_deref(),
            Some(&esperanto_flow::manifest::SpeciesManifest::single("mouse", baseline)),
        )?;
    }

    // --- model bundle (scoring) ---
    ensure_model_bundle()?;

    eprintln!("[setup] done. Run with:");
    eprintln!("  esperanto run --r1 <reads.fq.gz> --out out/");
    eprintln!("  esperanto run --r1 <reads.fq.gz> --out out/ --hybrid APOE4   (knock-in mouse)");
    Ok(())
}

/// Validate a fasta(+gtf) pair: both parse, contig names agree, the species
/// guardrail holds (manifest = Some for known non-human baselines).
fn validate_pair(
    fasta: &Path,
    gtf: Option<&Path>,
    manifest: Option<&esperanto_flow::manifest::SpeciesManifest>,
) -> anyhow::Result<()> {
    let reference = esperanto_map::fasta::parse_fasta(fasta)
        .map_err(|e| anyhow::anyhow!("FASTA {} does not parse: {e}", fasta.display()))?;
    if reference.contigs.is_empty() {
        anyhow::bail!("FASTA {} contains no sequences", fasta.display());
    }
    esperanto_flow::guard::check_species(fasta, manifest)
        .map_err(|e| anyhow::anyhow!("{}: {e}", fasta.display()))?;
    if let Some(g) = gtf {
        let set = esperanto_tidx::TranscriptSet::parse(g)
            .map_err(|e| anyhow::anyhow!("GTF {} does not parse: {e}", g.display()))?;
        if set.is_empty() {
            anyhow::bail!("GTF {} has no exon-bearing transcripts", g.display());
        }
        // Contig-name agreement: a systematic miss (e.g. FASTA uses `1` but
        // the GTF uses `chr1`) means the two files come from different
        // conventions and must not be built together.
        let mut missing: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
        let mut on_missing = 0usize;
        for i in 0..set.len() {
            let t = set.transcript(i).expect("i < len");
            if reference.contig_index(t.contig.as_bytes()).is_none() {
                missing.insert(t.contig.to_string());
                on_missing += 1;
            }
        }
        if on_missing * 2 > set.len() {
            let fa_names: Vec<&str> = reference
                .contigs
                .iter()
                .take(3)
                .map(|c| c.name.as_str())
                .collect();
            anyhow::bail!(
                "GTF contigs are absent from the FASTA (e.g. GTF has '{}', FASTA has '{}'): the two files use different naming conventions",
                missing.iter().next().expect("nonempty"),
                fa_names.join("', '")
            );
        }
        if !missing.is_empty() {
            eprintln!(
                "[setup] note: {on_missing}/{} transcripts sit on contigs not in the FASTA (e.g. '{}'); they are skipped at build",
                set.len(),
                missing.iter().next().expect("nonempty")
            );
        }
        eprintln!("[setup] {} transcripts, contigs agree", set.len());
    }
    Ok(())
}

/// The model bundle is required by the scoring stage. When no bundle
/// resolves (fresh source build without the installer), fetch the release
/// tarball, unpack only its `bundle/` tree into the user data dir, and
/// remove the archive.
fn ensure_model_bundle() -> anyhow::Result<()> {
    if crate::resolve::bundle(&None).is_ok() {
        eprintln!("[setup] model bundle: present");
        return Ok(());
    }
    let bundle_dir = crate::resolve::home_bundle_dir();
    let tmp = std::env::temp_dir().join("esperanto-release.tar.gz");
    download(RELEASE_TARBALL_URL, &tmp)?;
    let gz = flate2::read::GzDecoder::new(std::fs::File::open(&tmp)?);
    let mut ar = tar::Archive::new(gz);
    let mut unpacked = 0usize;
    for entry in ar.entries()? {
        let mut entry = entry?;
        let path = entry.path()?.into_owned();
        // Release layout: <pkg>/bundle/<...>. Keep the tree below bundle/.
        if path
            .components()
            .any(|c| matches!(c, std::path::Component::ParentDir | std::path::Component::RootDir))
        {
            anyhow::bail!("unsafe path in release tarball: {}", path.display());
        }
        let mut comps = path.components().map(|c| c.as_os_str().to_owned());
        if comps.next().is_none() || comps.next().as_deref() != Some("bundle".as_ref()) {
            continue;
        }
        let rel: PathBuf = comps.collect();
        if rel.as_os_str().is_empty() {
            continue;
        }
        let dest = bundle_dir.join(&rel);
        // Guard against path traversal in archives.
        if !dest.starts_with(&bundle_dir) {
            anyhow::bail!("unsafe path in release tarball: {}", path.display());
        }
        if let Some(parent) = dest.parent() {
            std::fs::create_dir_all(parent)?;
        }
        entry.unpack(&dest)?;
        unpacked += 1;
    }
    std::fs::remove_file(&tmp)?;
    if unpacked == 0 {
        anyhow::bail!("release tarball held no bundle/ tree");
    }
    eprintln!(
        "[setup] model bundle installed -> {} ({unpacked} files)",
        bundle_dir.display()
    );
    Ok(())
}
