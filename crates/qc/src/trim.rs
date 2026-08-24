//! Trimming stages applied per read, in fastp order:
//! fixed-position trim → BWA-style quality trim → polyG tail trim.
//!
//! All quality bytes are Phred+33 encoded ASCII.

/// Minimum polyG tail length for trimming to apply (fastp default).
pub const POLY_G_MIN_LEN: usize = 10;
/// Hard cap on tolerated non-G bases inside a polyG tail.
const POLY_G_BUDGET_CAP: usize = 5;
/// Budget grows by one for every this-many bases scanned.
const POLY_G_BUDGET_DIV: usize = 8;

/// Remove `front` bases from the 5' end and `tail` bases from the 3' end.
/// Sequence and quality are truncated identically; clamped to read length.
pub fn fixed_trim(seq: &mut Vec<u8>, qual: &mut Vec<u8>, front: usize, tail: usize) {
    let len = seq.len();
    let start = front.min(len);
    let end = len.saturating_sub(tail).max(start);
    let keep = end - start;
    seq.drain(..start);
    seq.truncate(keep);
    qual.drain(..start);
    qual.truncate(keep);
}

/// Phred value of a Phred+33 ASCII quality byte.
#[inline]
/// Phred value of a raw quality byte; bytes below 33 clamp to 0 (matches the
/// oracle's saturating arithmetic — stray control chars past the prescan
/// window must not produce negative phreds that skew the trim deficit).
fn phred(q: u8) -> i64 {
    i64::from(q.saturating_sub(33))
}

/// BWA-style 3' quality trim: number of bases to KEEP.
///
/// Scans from the tail accumulating `cutoff - phred(q)` (add-then-check);
/// stops as soon as the running sum goes negative. `cut` tracks the index at
/// which the running sum reached its maximum. A read that never dips below
/// `cutoff` (e.g. an all-Q20 read with cutoff 20) stays untrimmed because the
/// sum never exceeds the initial maximum of zero.
pub fn qtrim_tail(qual: &[u8], cutoff: u8) -> usize {
    let mut sum: i64 = 0;
    let mut max_sum: i64 = 0;
    let mut cut = qual.len();
    for pos in (0..qual.len()).rev() {
        sum += i64::from(cutoff) - phred(qual[pos]);
        if sum < 0 {
            break;
        }
        if sum > max_sum {
            max_sum = sum;
            cut = pos;
        }
    }
    cut
}

/// BWA-style 5' quality trim: number of bases to SKIP (mirror of [`qtrim_tail`]).
/// Unused by the pipeline (spec trims the 3' end only); kept as the
/// documented mirror of `qtrim_tail`.
#[allow(dead_code, clippy::needless_range_loop)]
pub fn qtrim_front(qual: &[u8], cutoff: u8) -> usize {
    let mut sum: i64 = 0;
    let mut max_sum: i64 = 0;
    let mut skip = 0usize;
    for pos in 0..qual.len() {
        sum += i64::from(cutoff) - phred(qual[pos]);
        if sum < 0 {
            break;
        }
        if sum > max_sum {
            max_sum = sum;
            skip = pos + 1;
        }
    }
    skip
}

/// Detect a polyG tail and return the number of bases to KEEP.
///
/// Scans backward from the tail; after scanning `k` bases the tolerated
/// number of non-G bases is `min(5, k / 8)`. Scanning stops at the first base
/// that pushes the mismatch count over budget. The detected run is the scanned
/// region excluding that breaking base; it is only trimmed when its length
/// reaches [`POLY_G_MIN_LEN`].
pub fn poly_g_trim(seq: &[u8]) -> usize {
    let len = seq.len();
    if len < POLY_G_MIN_LEN {
        return len;
    }
    let mut mismatches: usize = 0;
    let mut scanned: usize = 0;
    let mut break_pos: Option<usize> = None;
    for pos in (0..len).rev() {
        scanned += 1;
        // Scientific semantics: G in either case continues the tail.
        if !seq[pos].eq_ignore_ascii_case(&b'G') {
            mismatches += 1;
        }
        let budget = (scanned / POLY_G_BUDGET_DIV).min(POLY_G_BUDGET_CAP);
        if mismatches > budget {
            break_pos = Some(pos);
            break;
        }
    }
    let run = match break_pos {
        Some(pos) => len - pos - 1,
        None => len,
    };
    if run >= POLY_G_MIN_LEN {
        len - run
    } else {
        len
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quals(phreds: &[u8]) -> Vec<u8> {
        phreds.iter().map(|q| q + 33).collect()
    }

    #[test]
    fn q20_read_stays_untrimmed() {
        let q = quals(&[20; 50]);
        assert_eq!(qtrim_tail(&q, 20), 50);
        assert_eq!(qtrim_front(&q, 20), 0);
    }

    #[test]
    fn tail_garbage_is_cut() {
        let mut phreds = vec![30; 45];
        phreds.extend(vec![5; 5]);
        let q = quals(&phreds);
        assert_eq!(qtrim_tail(&q, 20), 45);
    }

    #[test]
    fn poly_g_pure_tail_is_trimmed() {
        // The mismatch budget (min(5, scanned/8) = 2 after the 20 G's) lets
        // the run absorb two trailing non-G bases (C, A) before the third
        // (T) breaks it: run = 22, so 8 leading bases are kept.
        let mut seq = b"ACGTACGTAC".to_vec();
        seq.extend(vec![b'G'; 20]);
        assert_eq!(poly_g_trim(&seq), 8);
    }

    #[test]
    fn poly_g_short_run_untouched() {
        let mut seq = b"ACGTACGTAC".to_vec();
        seq.extend(vec![b'G'; 6]);
        assert_eq!(poly_g_trim(&seq), seq.len());
    }

    #[test]
    fn poly_g_all_g_is_fully_trimmed() {
        let seq = vec![b'G'; 40];
        assert_eq!(poly_g_trim(&seq), 0);
    }
}
