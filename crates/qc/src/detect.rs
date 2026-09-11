//! SE adapter auto-detection: 12-mer seed counting over 3' tail windows,
//! then anchored consensus extension.
//!
//! Invariant exploited: read-through makes the adapter's 5' prefix the
//! high-frequency internal substring of read tails. The plain 3' suffix
//! drifts with insert length, so a suffix tree is unsuitable; a seed plus
//! anchored extension reconstructs the adapter from its 5' end.

/// Tail window length scanned per read.
const TAIL: usize = 36;
/// Seed k-mer length.
const K: usize = 12;
/// Minimum accepted candidate length.
const MIN_ADAPTER: usize = 10;

/// Bucket order A, C, G, T, N — also the deterministic tie-break order.
const IDX_BASE: &[u8; 5] = b"ACGTN";

fn base_idx(b: u8) -> usize {
    match b.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

/// 2-bit encode a k-mer (A=0,C=1,G=2,T=3); None when any base is non-ACGT.
fn encode(w: &[u8]) -> Option<u64> {
    let mut key = 0u64;
    for &b in w {
        let c = match b.to_ascii_uppercase() {
            b'A' => 0u64,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        key = (key << 2) | c;
    }
    Some(key)
}

fn decode(mut key: u64, k: usize) -> Vec<u8> {
    let mut out = vec![0u8; k];
    for slot in out.iter_mut().rev() {
        *slot = IDX_BASE[(key & 3) as usize];
        key >>= 2;
    }
    out
}

/// Absolute support floor: max(20, n/200).
fn support_of(n: usize) -> u32 {
    (n / 200).max(20) as u32
}

/// Index of the max bucket; ties resolve to the lowest index (A<C<G<T<N).
fn majority(cnt: &[u32; 5]) -> usize {
    let mut bi = 0usize;
    for (i, &c) in cnt.iter().enumerate().skip(1) {
        if c > cnt[bi] {
            bi = i;
        }
    }
    bi
}

/// Detect a consensus adapter from prescan read sequences.
///
/// Returns the adapter (uppercase) when a dominant seed extends into a
/// candidate of length >= 10 that `crate::match_adapter` verifies on
/// >= max(20, n/200) buffered reads; otherwise None.
pub(crate) fn detect_adapter(seqs: &[&[u8]]) -> Option<Vec<u8>> {
    let n = seqs.len();
    if n == 0 {
        return None;
    }
    let support = support_of(n);

    // Seed: most frequent 12-mer across tail windows (ties -> smallest code).
    let mut counts: std::collections::HashMap<u64, u32> = std::collections::HashMap::new();
    for &s in seqs {
        let tail = &s[s.len().saturating_sub(TAIL)..];
        if tail.len() < K {
            continue;
        }
        for w in 0..=(tail.len() - K) {
            if let Some(key) = encode(&tail[w..w + K]) {
                *counts.entry(key).or_insert(0) += 1;
            }
        }
    }
    let mut seed: Option<(u64, u32)> = None;
    for (&k, &c) in &counts {
        if c < support {
            continue;
        }
        match seed {
            Some((sk, sc)) if sc > c || (sc == c && sk < k) => {}
            _ => seed = Some((k, c)),
        }
    }
    let (seed_key, _) = seed?;
    let seed_seq = decode(seed_key, K);

    // Anchor each read at the leftmost seed occurrence inside its tail window.
    let mut anchors: Vec<(&[u8], usize)> = Vec::new();
    for &s in seqs {
        let t0 = s.len().saturating_sub(TAIL);
        if s.len() < K {
            continue;
        }
        for p in t0..=(s.len() - K) {
            if s[p..p + K].eq_ignore_ascii_case(&seed_seq) {
                anchors.push((s, p));
                break;
            }
        }
    }

    // Anchored consensus extension. A layer survives while its majority base
    // keeps >= support counts AND >= 60% of the reads contributing to that
    // layer (the relative share stops the walk from running into divergent
    // insert sequence past the adapter boundary; the absolute floor stops it
    // at read ends where contributors decay).
    let mut consensus = seed_seq;

    let mut left: Vec<u8> = Vec::new();
    let mut d = 1usize;
    loop {
        let mut cnt = [0u32; 5];
        for &(s, p) in &anchors {
            if p >= d {
                cnt[base_idx(s[p - d])] += 1;
            }
        }
        let bi = majority(&cnt);
        let total: u32 = cnt.iter().sum();
        if cnt[bi] < support || cnt[bi] * 5 < total * 3 {
            break;
        }
        left.push(IDX_BASE[bi]);
        d += 1;
    }
    left.reverse();
    left.append(&mut consensus);
    let mut full = left;

    let mut e = K;
    loop {
        let mut cnt = [0u32; 5];
        for &(s, p) in &anchors {
            if p + e < s.len() {
                cnt[base_idx(s[p + e])] += 1;
            }
        }
        let bi = majority(&cnt);
        let total: u32 = cnt.iter().sum();
        if cnt[bi] < support || cnt[bi] * 5 < total * 3 {
            break;
        }
        full.push(IDX_BASE[bi]);
        e += 1;
    }

    if full.len() < MIN_ADAPTER {
        return None;
    }
    // Verify with the production matcher over the prescan buffer.
    let hits = seqs
        .iter()
        .filter(|s| crate::match_adapter(s, &full).is_some())
        .count() as u32;
    if hits >= support {
        Some(full)
    } else {
        None
    }
}
