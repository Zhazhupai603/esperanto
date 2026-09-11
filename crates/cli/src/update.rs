//! `esperanto update` — check the latest GitHub release against the running
//! version and, with confirmation, replace the binary (+ model bundle) in
//! place. `--check` reports without touching anything.

use std::path::PathBuf;

use clap::Args;

const LATEST_API: &str =
    "https://api.github.com/repos/Zhazhupai603/esperanto/releases/latest";

#[derive(Args)]
pub struct UpdateArgs {
    /// Only check; do not download or replace anything.
    #[arg(long)]
    check: bool,
    /// Skip the confirmation prompt (scripts only).
    #[arg(long)]
    yes: bool,
}

/// Parse "v1.0.1" / "1.0.1" into a comparable tuple.
fn version_tuple(v: &str) -> (u64, u64, u64) {
    let v = v.trim_start_matches('v');
    let mut it = v.split('.');
    let n = |s: Option<&str>| s.and_then(|x| x.parse().ok()).unwrap_or(0);
    (n(it.next()), n(it.next()), n(it.next()))
}

fn latest_release() -> anyhow::Result<(String, String)> {
    let agent = ureq::Agent::config_builder()
        .proxy(ureq::Proxy::try_from_env())
        .build()
        .new_agent();
    let mut resp = agent
        .get(LATEST_API)
        .header("User-Agent", "esperanto-cli")
        .header("Accept", "application/vnd.github+json")
        .call()?;
    let text = resp.body_mut().read_to_string()?;
    // minimal parse: first "tag_name": "<value>" pair
    let tag = text
        .split("\"tag_name\":")
        .nth(1)
        .and_then(|rest| rest.trim_start().strip_prefix('"'))
        .and_then(|rest| rest.split('"').next())
        .ok_or_else(|| anyhow::anyhow!("release payload has no tag_name"))?
        .to_string();
    Ok((tag.clone(), tag))
}

fn confirm_update(current: &str, latest: &str, yes: bool) -> anyhow::Result<bool> {
    if !yes && !std::io::IsTerminal::is_terminal(&std::io::stdin()) {
        eprintln!("refusing to update without a TTY; pass --yes to confirm");
        return Ok(false);
    }
    crate::confirm::confirmed(
        &format!("Update esperanto {current} -> {latest}?"),
        ("No, stay on current version", "Yes, update now"),
        yes,
    )
}

/// Replace the running binary with the tarball's `bin/esperanto` and merge
/// the tarball's `bundle/` tree into the user data dir.
fn install_tarball(bytes: &[u8], tag: &str) -> anyhow::Result<()> {
    let gz = flate2::read::GzDecoder::new(bytes);
    let mut ar = tar::Archive::new(gz);
    let exe = std::env::current_exe()?;
    let bundle_dir = crate::resolve::home_bundle_dir();
    let tmp = std::env::temp_dir().join(format!("esperanto-update-{tag}"));
    std::fs::create_dir_all(&tmp)?;
    ar.unpack(&tmp)?;
    // locate <tmp>/<pkg>/bin/esperanto and <tmp>/<pkg>/bundle
    let pkg = std::fs::read_dir(&tmp)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .find(|p| p.is_dir())
        .ok_or_else(|| anyhow::anyhow!("update tarball had no package directory"))?;
    let new_bin = pkg.join("bin").join("esperanto");
    if !new_bin.is_file() {
        anyhow::bail!("update tarball held no bin/esperanto");
    }
    // binary swap: copy + rename (atomic on the same filesystem)
    let staged = exe.with_extension("new");
    std::fs::copy(&new_bin, &staged)?;
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&staged, std::fs::Permissions::from_mode(0o755))?;
    }
    std::fs::rename(&staged, &exe)?;
    let new_bundle = pkg.join("bundle");
    if new_bundle.is_dir() {
        let dest_parent = bundle_dir.parent().map(|p| p.to_path_buf()).unwrap_or(bundle_dir.clone());
        std::fs::create_dir_all(&dest_parent)?;
        copy_tree(&new_bundle, &bundle_dir)?;
    }
    std::fs::remove_dir_all(&tmp).ok();
    Ok(())
}

fn copy_tree(src: &PathBuf, dst: &PathBuf) -> anyhow::Result<()> {
    std::fs::create_dir_all(dst)?;
    for e in std::fs::read_dir(src)? {
        let e = e?;
        let p = e.path();
        let d = dst.join(e.file_name());
        if p.is_dir() {
            copy_tree(&p, &d)?;
        } else {
            std::fs::copy(&p, &d)?;
        }
    }
    Ok(())
}

pub fn run(a: UpdateArgs) -> anyhow::Result<()> {
    let current = env!("CARGO_PKG_VERSION").to_string();
    let (tag, _name) = latest_release()?;
    let latest_v = tag.trim_start_matches('v');
    eprintln!("current: {current}  latest: {latest_v}");
    if version_tuple(&current) >= version_tuple(latest_v) {
        eprintln!("already up to date");
        return Ok(());
    }
    if a.check {
        eprintln!("a newer release is available: {latest_v}");
        return Ok(());
    }
    if !confirm_update(&current, latest_v, a.yes)? {
        eprintln!("aborted; nothing changed");
        return Ok(());
    }
    let url = format!(
        "https://github.com/Zhazhupai603/esperanto/releases/download/{tag}/esperanto-{latest_v}-linux-x86_64.tar.gz"
    );
    eprintln!("[update] downloading {url}");
    let agent = ureq::Agent::config_builder()
        .proxy(ureq::Proxy::try_from_env())
        .build()
        .new_agent();
    let mut resp = agent.get(&url).call()?;
    let bytes = resp.body_mut().read_to_vec()?;
    install_tarball(&bytes, &tag)?;
    eprintln!("[update] esperanto is now {latest_v}");
    Ok(())
}
