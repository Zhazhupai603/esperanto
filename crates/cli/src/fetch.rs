//! Parallel chunked downloader with per-chunk resume and mirror fallback.
//!
//! Large reference downloads dominate first-run setup time on throttled
//! links (per-connection limits). `file()` splits the transfer into
//! fixed-size chunks fetched over concurrent connections, tracks completed
//! chunks in a sidecar file, and resumes where a previous attempt stopped.
//! Mirrors are tried in order; the first reachable one wins.

use std::io::{Read, Seek, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use anyhow::{anyhow, Context};

/// Chunk size for parallel fetches.
const CHUNK: u64 = 8 << 20;
/// Concurrent connections.
const THREADS: usize = 6;
/// Per-chunk retry attempts.
const RETRIES: u32 = 6;

fn agent() -> ureq::Agent {
    ureq::Agent::config_builder()
        .timeout_per_call(Some(std::time::Duration::from_secs(600)))
        .proxy(ureq::Proxy::try_from_env())
        .build()
        .new_agent()
}

/// Content length + range support probe. Ok(None) length = unknown
/// (falls back to a single sequential stream).
fn probe(url: &str) -> anyhow::Result<Option<u64>> {
    let a = agent();
    let len = a.head(url).call().ok().and_then(|r| {
        r.headers()
            .get("content-length")
            .and_then(|v| v.to_str().ok())
            .and_then(|v| v.parse::<u64>().ok())
    });
    // Range check: a 1-byte ranged GET must answer 206.
    if let Some(total) = len {
        if total > CHUNK {
            let ok = a
                .get(url)
                .header("Range", "bytes=0-0")
                .call()
                .map(|r| r.status().as_u16() == 206)
                .unwrap_or(false);
            if !ok {
                return Ok(None);
            }
        }
    }
    Ok(len)
}

/// Sidecar progress: one byte per chunk ('1' = done).
struct Progress {
    path: PathBuf,
    marks: Vec<u8>,
}

impl Progress {
    fn load(dest: &Path, total_chunks: usize) -> Self {
        let path = dest.with_extension("part.progress");
        let marks = std::fs::read(&path)
            .ok()
            .filter(|v| v.len() == total_chunks)
            .unwrap_or_else(|| vec![b'0'; total_chunks]);
        Progress { path, marks }
    }
    fn done_count(&self) -> usize {
        self.marks.iter().filter(|&&m| m == b'1').count()
    }
    fn mark(&mut self, i: usize) {
        self.marks[i] = b'1';
        let tmp = self.path.with_extension("part.progress.tmp");
        if std::fs::write(&tmp, &self.marks).is_ok() {
            let _ = std::fs::rename(&tmp, &self.path);
        }
    }
    fn clear(&self) {
        let _ = std::fs::remove_file(&self.path);
    }
}

/// Download `url` to `dest` (parallel chunks + resume). Prints progress
/// lines to stderr.
pub fn file(url: &str, dest: &Path) -> anyhow::Result<()> {
    let name = url.rsplit('/').next().unwrap_or(url);
    let total = probe(url).context("probe download")?;
    let Some(total) = total else {
        // No length or no Range support: plain sequential stream.
        eprintln!("[fetch] {name} (sequential)");
        let resp = agent().get(url).call()?;
        let mut reader = resp.into_body().into_reader();
        let mut out = std::io::BufWriter::new(std::fs::File::create(dest)?);
        std::io::copy(&mut reader, &mut out)?;
        out.flush()?;
        return Ok(());
    };
    let total_chunks = total.div_ceil(CHUNK) as usize;
    let progress = Progress::load(dest, total_chunks);
    if progress.done_count() == total_chunks
        && dest.metadata().map(|m| m.len() == total).unwrap_or(false)
    {
        progress.clear();
        return Ok(());
    }
    eprintln!(
        "[fetch] {name} ({} MB, {THREADS} streams, {}/{} chunks done)",
        total / 1048576,
        progress.done_count(),
        total_chunks
    );
    std::fs::write(dest, b"")?; // fresh handle; chunks seek to their offsets
    let file = std::fs::OpenOptions::new().write(true).open(dest)?;
    file.set_len(total)?;

    let next = AtomicUsize::new(0);
    let progress = Mutex::new(progress);
    let failed: Mutex<Vec<(usize, String)>> = Mutex::new(Vec::new());
    let url = &url.to_string();
    std::thread::scope(|s| {
        for _ in 0..THREADS {
            s.spawn(|| {
                let a = agent();
                loop {
                    let i = next.fetch_add(1, Ordering::SeqCst);
                    if i >= total_chunks {
                        break;
                    }
                    if progress.lock().unwrap().marks[i] == b'1' {
                        continue;
                    }
                    let start = i as u64 * CHUNK;
                    let end = (start + CHUNK - 1).min(total - 1);
                    let mut last_err = String::new();
                    let mut ok = false;
                    for _ in 0..RETRIES {
                        match try_chunk(&a, url.as_str(), dest, start, end) {
                            Ok(()) => {
                                ok = true;
                                break;
                            }
                            Err(e) => last_err = e.to_string(),
                        }
                        std::thread::sleep(std::time::Duration::from_secs(2));
                    }
                    if ok {
                        let mut marks = progress.lock().unwrap();
                        marks.mark(i);
                        let done = marks.done_count();
                        if done.is_multiple_of(8) || done == total_chunks {
                            eprintln!(
                                "[fetch]   {done}/{total_chunks} chunks (~{} MB)",
                                done as u64 * CHUNK / 1048576
                            );
                        }
                    } else {
                        failed.lock().unwrap().push((i, last_err));
                        break;
                    }
                }
            });
        }
    });
    let failed = failed.into_inner().unwrap();
    if !failed.is_empty() {
        let (i, err) = &failed[0];
        anyhow::bail!("chunk {i} failed after {RETRIES} attempts: {err}");
    }
    // Definitive completeness check: every chunk must be marked done (a
    // worker exiting early would otherwise leave silent sparse holes).
    let marks = progress.lock().unwrap();
    let done = marks.done_count();
    if done != total_chunks {
        anyhow::bail!("incomplete download: {done}/{total_chunks} chunks (rerun to resume)");
    }
    marks.clear();
    Ok(())
}

/// Fetch one chunk into its file offset. Read errors (connection drops
/// mid-body) and short bodies are errors -> the caller retries the chunk.
fn try_chunk(
    a: &ureq::Agent,
    url: &str,
    dest: &Path,
    start: u64,
    end: u64,
) -> anyhow::Result<()> {
    let resp = a
        .get(url)
        .header("Range", format!("bytes={start}-{end}"))
        .call()
        .map_err(|e| anyhow!("range get: {e}"))?;
    if resp.status().as_u16() != 206 {
        anyhow::bail!("server ignored the range request ({})", resp.status());
    }
    let mut reader = resp.into_body().into_reader();
    let mut f = std::fs::OpenOptions::new()
        .write(true)
        .open(dest)
        .with_context(|| format!("open {}", dest.display()))?;
    f.seek(std::io::SeekFrom::Start(start))?;
    let mut buf = [0u8; 1 << 16];
    let mut copied = 0u64;
    loop {
        let n = reader
            .read(&mut buf)
            .map_err(|e| anyhow!("read at +{copied}: {e}"))?;
        if n == 0 {
            break;
        }
        f.write_all(&buf[..n])?;
        copied += n as u64;
    }
    if copied != end - start + 1 {
        anyhow::bail!("short body: {copied}/{} bytes", end - start + 1);
    }
    Ok(())
}

/// Try each mirror in order; the first successful download wins.
pub fn from_mirrors(mirrors: &[&str], dest: &Path) -> anyhow::Result<()> {
    let mut last = String::new();
    for url in mirrors {
        match file(url, dest) {
            Ok(()) => return Ok(()),
            Err(e) => {
                last = format!("{url}: {e}");
                eprintln!("[fetch] mirror failed, trying next: {last}");
            }
        }
    }
    Err(anyhow!("all mirrors failed; last: {last}"))
}

