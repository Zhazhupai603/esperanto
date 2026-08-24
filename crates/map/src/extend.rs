//! Affine (Gotoh) banded extension with editing-aware substitution.
//!
//! The DP is Smith-Waterman-shaped: an alignment may start (score 0) and end
//! (global max) anywhere, with unaligned read prefixes/suffixes written as
//! soft clips — CIGAR `M/I/S` ops consume exactly the read length.
//!
//! Two execution paths produce **identical scores and CIGARs** on any input:
//! `run_banded_packed` (i16 accumulation + substitution LUT) runs when every
//! byte of both sequences is one of {A,C,G,T,N} (uppercase); any other byte
//! (lowercase, ambiguity codes) routes to `run_banded_legacy`, the byte-wise
//! i32 path. Band scheduling is exactly two attempts — half-width 30, then
//! the full matrix — with the second attempt gated on the first attempt's
//! score falling below 85% of the perfect read score (see `extend_hint`).

/// Extension scoring parameters (frozen: match 2, mismatch −4, gap open 4,
/// extension 2).
#[derive(Clone, Copy, Debug)]
pub struct ExtendParams {
    /// Match reward.
    pub match_score: i32,
    /// Mismatch penalty (stored negative).
    pub mismatch: i32,
    /// Gap opening penalty (stored positive; subtracted).
    pub gap_open: i32,
    /// Gap extension penalty per base (stored positive; subtracted).
    pub gap_ext: i32,
    /// Editing-aware substitution: read A vs ref G and read T vs ref C score
    /// 0 (A-to-I RNA editing tolerance).
    pub editing_aware: bool,
}

impl Default for ExtendParams {
    fn default() -> ExtendParams {
        ExtendParams {
            match_score: 2,
            mismatch: -4,
            gap_open: 4,
            gap_ext: 2,
            editing_aware: false,
        }
    }
}

/// First (and only) banded attempt's half-width; the second attempt is the
/// full matrix.
pub const INITIAL_BAND: i64 = 30;

/// CIGAR operation (SAM; `=`/`X` folded into `Match`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CigarOp {
    /// Aligned bases (match or mismatch).
    Match(u32),
    /// Insertion to the reference (consumes read only).
    Ins(u32),
    /// Deletion from the reference (consumes ref only).
    Del(u32),
    /// Reference skip (intron; consumes ref only).
    RefSkip(u32),
    /// Soft-clipped read bases (consume read only).
    SoftClip(u32),
}

/// 5×5 substitution LUT over base codes {A,C,G,T,N} × {A,C,G,T,N}, indexed
/// `[read_code * 5 + ref_code]`.
#[derive(Clone, Copy)]
pub struct SubstLut {
    scores: [i16; 25],
}

/// Base code: A=0 C=1 G=2 T=3, anything else (incl. N) = 4.
#[inline]
pub fn base_code(b: u8) -> usize {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

impl SubstLut {
    /// Build the LUT from parameters; editing-aware mode zeroes (A,G) and
    /// (T,C) read-vs-ref scores.
    pub fn new(params: &ExtendParams) -> SubstLut {
        let m = params.match_score as i16;
        let mm = params.mismatch as i16;
        let mut scores = [mm; 25];
        for c in 0..4 {
            scores[c * 5 + c] = m;
        }
        if params.editing_aware {
            // EA-free pairs (frozen): ref A/read G and ref T/read C
            // (A>I editing appears as G in the read). LUT layout is
            // scores[read_code*5 + ref_code].
            scores[10] = 0; // read G(2), ref A(0): 2*5+0
            scores[8] = 0; // read C(1), ref T(3): 1*5+3
        }
        SubstLut { scores }
    }

    /// Score for read base code vs ref base code.
    #[inline]
    pub fn score(&self, read_code: usize, ref_code: usize) -> i16 {
        self.scores[read_code * 5 + ref_code]
    }
}

/// Diagonal hint from chaining: the anchor diagonal `ref − read` plus the
/// supporting-anchor fraction (num/den).
#[derive(Clone, Copy, Debug)]
pub struct DiagHint {
    /// Band center: `j − i` (ref index minus read index).
    pub offset: i64,
    /// Anchors supporting this diagonal.
    pub num: u32,
    /// Total anchors considered.
    pub den: u32,
}

/// Result of one extension.
#[derive(Clone, Debug)]
pub struct Extension {
    /// First aligned read base (leading soft clip before this).
    pub read_start: u32,
    /// One past the last aligned read base.
    pub read_end: u32,
    /// First aligned reference base, relative to the window start.
    pub ref_start: u32,
    /// CIGAR over the whole read (soft clips included); `M/I/S` consume
    /// exactly `read_len` bases.
    pub cigar: Vec<CigarOp>,
    /// DP score (global maximum).
    pub score: i32,
}

/// Reusable per-thread DP buffers (no semantic effect).
pub struct ExtendBuffer {
    packed: BandBuf<i16>,
    legacy: BandBuf<i32>,
}

impl Default for ExtendBuffer {
    fn default() -> ExtendBuffer {
        ExtendBuffer {
            packed: BandBuf::new(),
            legacy: BandBuf::new(),
        }
    }
}

impl ExtendBuffer {
    /// Fresh buffer.
    pub fn new() -> ExtendBuffer {
        ExtendBuffer::default()
    }
}

/// Append `op` to `ops`, merging with the previous op when kinds match.
pub fn push_op(ops: &mut Vec<CigarOp>, op: CigarOp) {
    let same = matches!(
        (ops.last(), op),
        (Some(CigarOp::Match(_)), CigarOp::Match(_))
            | (Some(CigarOp::Ins(_)), CigarOp::Ins(_))
            | (Some(CigarOp::Del(_)), CigarOp::Del(_))
            | (Some(CigarOp::RefSkip(_)), CigarOp::RefSkip(_))
            | (Some(CigarOp::SoftClip(_)), CigarOp::SoftClip(_))
    );
    if same {
        let len = match (ops.last_mut().unwrap(), op) {
            (CigarOp::Match(a), CigarOp::Match(b))
            | (CigarOp::Ins(a), CigarOp::Ins(b))
            | (CigarOp::Del(a), CigarOp::Del(b))
            | (CigarOp::RefSkip(a), CigarOp::RefSkip(b))
            | (CigarOp::SoftClip(a), CigarOp::SoftClip(b)) => *a + b,
            _ => unreachable!(),
        };
        match ops.last_mut().unwrap() {
            CigarOp::Match(a)
            | CigarOp::Ins(a)
            | CigarOp::Del(a)
            | CigarOp::RefSkip(a)
            | CigarOp::SoftClip(a) => *a = len,
        }
    } else {
        ops.push(op);
    }
}

/// Whether a byte qualifies for the packed path: uppercase {A,C,G,T,N} only.
fn pure_byte(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'N')
}

/// Path dispatch: pure {A,C,G,T,N} inputs run the i16+LUT packed path; any
/// other byte (lowercase, ambiguity codes) routes to the byte-wise legacy
/// path. Both paths score every byte pair identically, so the produced
/// `Extension` is the same either way.
fn run_dispatch(
    read: &[u8],
    ref_window: &[u8],
    params: &ExtendParams,
    buf: &mut ExtendBuffer,
    band_spec: Option<(i64, i64, i64, i64)>,
) -> (Extension, bool) {
    let core = CoreParams::new(params);
    if read.iter().all(|&b| pure_byte(b)) && ref_window.iter().all(|&b| pure_byte(b)) {
        run_band(read, ref_window, &core, &mut buf.packed, band_spec)
    } else {
        run_band(read, ref_window, &core, &mut buf.legacy, band_spec)
    }
}

/// Full-window extension.
pub fn extend(
    read: &[u8],
    ref_window: &[u8],
    params: &ExtendParams,
    buf: &mut ExtendBuffer,
) -> Extension {
    run_dispatch(read, ref_window, params, buf, None).0
}

/// Banded extension along a diagonal hint. The schedule is exactly two
/// attempts, `[30, usize::MAX]`: attempt 1 runs half-band 30 around the hint
/// diagonal; if its score passes the escalation gate
/// `score × 20 ≥ perfect × 17` (with `perfect = read_len × match_score`, i.e.
/// at least 85% of the perfect read score) it is returned immediately —
/// otherwise attempt 2 recomputes the full matrix and is returned
/// unconditionally. A band that already covers the whole window
/// (`30 ≥ ref_window.len()`) is accepted directly as the full matrix.
pub fn extend_hint(
    read: &[u8],
    ref_window: &[u8],
    params: &ExtendParams,
    buf: &mut ExtendBuffer,
    hint: DiagHint,
) -> Extension {
    let m = read.len() as i64;
    let n = ref_window.len() as i64;
    if m == 0 || n == 0 {
        return empty_extension(m as u32);
    }
    // Band already covers the whole window: full matrix, accept directly.
    if INITIAL_BAND >= n {
        return run_dispatch(read, ref_window, params, buf, None).0;
    }
    // Attempt 1: banded.
    let (ext, _) = run_dispatch(
        read,
        ref_window,
        params,
        buf,
        Some((hint.offset, hint.num as i64, hint.den.max(1) as i64, INITIAL_BAND)),
    );
    // Escalation gate: score < 85% of the perfect read score ⇒ full matrix.
    let perfect = read.len() as i32 * params.match_score;
    if ext.score * 20 < perfect * 17 {
        return run_dispatch(read, ref_window, params, buf, None).0;
    }
    ext
}

/// i16+LUT banded DP. `band_spec = None` runs the full window; `Some((offset,
/// half_width))` restricts cells to `|j − i − offset| <= half_width`. Also
/// reports whether any row's best H sat on a strict band edge.
pub fn run_banded_packed(
    read: &[u8],
    ref_window: &[u8],
    params: &ExtendParams,
    buf: &mut ExtendBuffer,
    band_spec: Option<(i64, i64, i64, i64)>,
) -> (Extension, bool) {
    let core = CoreParams::new(params);
    run_band(read, ref_window, &core, &mut buf.packed, band_spec)
}

/// Legacy byte-wise i32 path — taken automatically by [`extend`] /
/// [`extend_hint`] on non-{A,C,G,T,N} input, and directly callable; must
/// produce identical scores and CIGARs to [`run_banded_packed`] for the same
/// band specification.
pub fn run_banded_legacy(
    read: &[u8],
    ref_window: &[u8],
    params: &ExtendParams,
    buf: &mut ExtendBuffer,
    band_spec: Option<(i64, i64, i64, i64)>,
) -> (Extension, bool) {
    let core = CoreParams::new(params);
    run_band(read, ref_window, &core, &mut buf.legacy, band_spec)
}
fn empty_extension(read_len: u32) -> Extension {
    Extension {
        read_start: 0,
        read_end: 0,
        ref_start: 0,
        cigar: if read_len > 0 {
            vec![CigarOp::SoftClip(read_len)]
        } else {
            Vec::new()
        },
        score: 0,
    }
}

struct CoreParams {
    lut: SubstLut,
    open_ext: i32,
    ext: i32,
}

impl CoreParams {
    fn new(params: &ExtendParams) -> CoreParams {
        CoreParams {
            lut: SubstLut::new(params),
            open_ext: -params.gap_open, // first gap base pays gap_open only (frozen legacy model)
            ext: -params.gap_ext,
        }
    }
}

trait Acc: Copy {
    const NEG: Self;
    fn is_neg(x: Self) -> bool;
    fn same(a: Self, b: Self) -> bool;
    fn better(a: Self, b: Self) -> bool;
    fn add(self, o: Self) -> Self;
    fn from_i32(v: i32) -> Self;
}

impl Acc for i16 {
    const NEG: i16 = -30000;
    #[inline]
    fn is_neg(x: i16) -> bool {
        x == i16::NEG
    }
    #[inline]
    fn same(a: i16, b: i16) -> bool {
        a == b
    }
    #[inline]
    fn better(a: i16, b: i16) -> bool {
        a > b
    }
    #[inline]
    fn add(self, o: i16) -> i16 {
        if Acc::is_neg(self) || Acc::is_neg(o) {
            i16::NEG
        } else {
            (self + o).max(i16::NEG)
        }
    }
    #[inline]
    fn from_i32(v: i32) -> i16 {
        v as i16
    }
}

impl Acc for i32 {
    const NEG: i32 = -1_000_000_000;
    #[inline]
    fn is_neg(x: i32) -> bool {
        x == i32::NEG
    }
    #[inline]
    fn same(a: i32, b: i32) -> bool {
        a == b
    }
    #[inline]
    fn better(a: i32, b: i32) -> bool {
        a > b
    }
    #[inline]
    fn add(self, o: i32) -> Self {
        if Acc::is_neg(self) || Acc::is_neg(o) {
            i32::NEG
        } else {
            self + o
        }
    }
    #[inline]
    fn from_i32(v: i32) -> i32 {
        v
    }
}

struct RowTrace {
    lo: i64,
    d: Vec<u8>, // per cell: hdir | edir<<2 | fdir<<4
}

struct BandBuf<T: Acc> {
    h_prev: Vec<T>,
    h_cur: Vec<T>,
    e_prev: Vec<T>,
    e_cur: Vec<T>,
    f_cur: Vec<T>,
    rows: Vec<RowTrace>,
}

impl<T: Acc> BandBuf<T> {
    fn new() -> BandBuf<T> {
        BandBuf {
            h_prev: Vec::new(),
            h_cur: Vec::new(),
            e_prev: Vec::new(),
            e_cur: Vec::new(),
            f_cur: Vec::new(),
            rows: Vec::new(),
        }
    }
}

const H_DIAG: u8 = 0;
const H_FROM_E: u8 = 1;
const H_FROM_F: u8 = 2;
const H_START: u8 = 3;

fn run_band<T: Acc>(
    read: &[u8],
    rf: &[u8],
    core: &CoreParams,
    buf: &mut BandBuf<T>,
    band_spec: Option<(i64, i64, i64, i64)>,
) -> (Extension, bool) {
    let m = read.len() as i64;
    let n = rf.len() as i64;
    if m == 0 || n == 0 {
        return (empty_extension(m as u32), false);
    }

    let (offset, num, den, half) = match band_spec {
        Some((o, nm, dn, h)) => (o, nm, dn.max(1), h.max(0)),
        None => (0, 1, 1, m + n),
    };
    let width = (2 * half + 1) as usize;
    if buf.h_prev.len() != width {
        buf.h_prev = vec![T::NEG; width];
        buf.h_cur = vec![T::NEG; width];
        buf.e_prev = vec![T::NEG; width];
        buf.e_cur = vec![T::NEG; width];
        buf.f_cur = vec![T::NEG; width];
    } else {
        for v in buf.h_prev.iter_mut() {
            *v = T::NEG;
        }
        for v in buf.h_cur.iter_mut() {
            *v = T::NEG;
        }
        for v in buf.e_prev.iter_mut() {
            *v = T::NEG;
        }
        for v in buf.e_cur.iter_mut() {
            *v = T::NEG;
        }
        for v in buf.f_cur.iter_mut() {
            *v = T::NEG;
        }
    }
    buf.rows.clear();

    // Band center: offset + i*num/den (truncating division, legacy formula).
    // Compute rows (i>=1) clamp lo at 1 — column 0 is the pinned H=0
    // boundary and never enters the DP; row 0 itself is zeroed from col 0.
    let row_window = |i: i64| -> (i64, i64) {
        let center = offset + i * num / den;
        if i == 0 {
            ((center - half).max(0), (center + half).min(n))
        } else {
            let c = center.clamp(1, n);
            ((c - half).max(1), (c + half).min(n))
        }
    };

    // Row 0: H = 0 within the window (free reference prefix).
    {
        let (lo, hi) = row_window(0);
        if lo <= hi {
            let d = vec![H_START; (hi - lo + 1) as usize];
            for c in 0..d.len() {
                buf.h_prev[c] = T::from_i32(0);
            }
            buf.rows.push(RowTrace { lo, d });
        } else {
            buf.rows.push(RowTrace {
                lo,
                d: Vec::new(),
            });
        }
    }

    let open_ext = T::from_i32(core.open_ext);
    let ext_p = T::from_i32(core.ext);
    let zero = T::from_i32(0);

    let mut best_score = zero;
    let mut best_i = 0i64;
    let mut best_j = 0i64;
    let mut touched = false;

    for i in 1..=m {
        let (lo, hi) = row_window(i);
        let (plo, phi) = row_window(i - 1);
        if lo > hi {
            buf.rows.push(RowTrace { lo, d: Vec::new() });
            for v in buf.h_prev.iter_mut() {
                *v = T::NEG;
            }
            for v in buf.e_prev.iter_mut() {
                *v = T::NEG;
            }
            continue;
        }
        let mut d = vec![0u8; (hi - lo + 1) as usize];
        let mut row_max = T::NEG;
        let mut row_max_cell = 0i64;

        for jj in lo..=hi {
            let c = (jj - lo) as usize;

            // E: vertical (insertion): from H[i-1][j] or E[i-1][j].
            // Out-of-band H reads see the margin/boundary value 0 (legacy
            // writes margins and col 0 as 0); E reads stay NEG.
            let in_prev = jj >= plo && jj <= phi;
            let (h_up, e_up) = if in_prev {
                (
                    buf.h_prev[(jj - plo) as usize],
                    buf.e_prev[(jj - plo) as usize],
                )
            } else {
                (zero, T::NEG)
            };
            let e_open = h_up.add(open_ext);
            let e_ext = e_up.add(ext_p);
            let (e_val, edir) = if T::better(e_ext, e_open) {
                (e_ext, 1u8)
            } else {
                (e_open, 0u8)
            };

            // F: horizontal (deletion): from H[i][j-1] or F[i][j-1].
            // F: horizontal (deletion): from H[i][j-1] or F[i][j-1]; the
            // left margin/boundary cell reads H=0 (legacy), F stays NEG.
            let (h_left, f_left) = if jj > lo {
                (buf.h_cur[c - 1], buf.f_cur[c - 1])
            } else {
                (zero, T::NEG)
            };
            let f_open = h_left.add(open_ext);
            let f_ext = f_left.add(ext_p);
            let (f_val, fdir) = if T::better(f_ext, f_open) {
                (f_ext, 1u8)
            } else {
                (f_open, 0u8)
            };

            // Diagonal: the H read never fails — in-band cells are stored,
            // col 0 is the pinned 0 boundary, margin cells read 0 (legacy).
            let h_prev_idx = jj - 1;
            let h_diag = if h_prev_idx >= plo && h_prev_idx <= phi {
                buf.h_prev[(h_prev_idx - plo) as usize]
            } else {
                zero
            };
            let s = core.lut.score(
                base_code(read[(i - 1) as usize]),
                base_code(rf[(jj - 1) as usize]),
            );
            let diag = h_diag.add(T::from_i32(s as i32));

            // H: preference order on ties: diag > E > F > start.
            let mut h_val = zero;
            let mut hdir = H_START;
            if !T::is_neg(f_val) && !T::better(h_val, f_val) {
                h_val = f_val;
                hdir = H_FROM_F;
            }
            if !T::is_neg(e_val) && !T::better(h_val, e_val) {
                h_val = e_val;
                hdir = H_FROM_E;
            }
            if !T::is_neg(diag) && !T::better(h_val, diag) {
                h_val = diag;
                hdir = H_DIAG;
            }
            // Local-alignment zero floor: a zero-valued cell stops the walk
            // regardless of source (checked after the tie cascade).
            if T::same(h_val, zero) {
                hdir = H_START;
            }

            buf.h_cur[c] = h_val;
            buf.e_cur[c] = e_val;
            buf.f_cur[c] = f_val;
            d[c] = hdir | (edir << 2) | (fdir << 4);

            if T::better(h_val, row_max) {
                row_max = h_val;
                row_max_cell = c as i64;
            }
        }

        // Strict band edge touch (edge exists and holds the row maximum).
        let left_strict = lo > 0;
        let right_strict = hi < n;
        if (left_strict && T::same(buf.h_cur[0], row_max))
            || (right_strict && T::same(buf.h_cur[(hi - lo) as usize], row_max))
        {
            touched = true;
        }

        if T::better(row_max, best_score) {
            best_score = row_max;
            best_i = i;
            best_j = lo + row_max_cell;
        }

        buf.rows.push(RowTrace { lo, d });
        std::mem::swap(&mut buf.h_prev, &mut buf.h_cur);
        std::mem::swap(&mut buf.e_prev, &mut buf.e_cur);
        for v in buf.h_cur.iter_mut() {
            *v = T::NEG;
        }
        for v in buf.e_cur.iter_mut() {
            *v = T::NEG;
        }
        for v in buf.f_cur.iter_mut() {
            *v = T::NEG;
        }
    }

    // State-machine traceback (H/E/F) recording per-step deltas for end
    // trimming. Cell bits: low2 = hdir; bit2 = E-from-extension; bit4 =
    // F-from-extension. Gap-open cost lands on the first step of a run.
    #[derive(PartialEq)]
    enum Tb {
        H,
        E,
        F,
    }
    let mut steps: Vec<(CigarOp, i32)> = Vec::new(); // (op, delta), reversed
    let mut state = Tb::H;
    let mut i = best_i;
    let mut j = best_j;
    loop {
        if i <= 0 || j <= 0 {
            break;
        }
        let row = &buf.rows[i as usize];
        if j < row.lo || j >= row.lo + row.d.len() as i64 {
            break;
        }
        let cell = row.d[(j - row.lo) as usize];
        match state {
            Tb::H => match cell & 3 {
                H_DIAG => {
                    let s = core.lut.score(
                        base_code(read[(i - 1) as usize]),
                        base_code(rf[(j - 1) as usize]),
                    );
                    steps.push((CigarOp::Match(1), s as i32));
                    i -= 1;
                    j -= 1;
                }
                H_FROM_E => {
                    steps.push((CigarOp::Ins(1), core.open_ext));
                    state = if (cell >> 2) & 1 == 1 { Tb::E } else { Tb::H };
                    i -= 1;
                }
                H_FROM_F => {
                    steps.push((CigarOp::Del(1), core.open_ext));
                    state = if (cell >> 4) & 1 == 1 { Tb::F } else { Tb::H };
                    j -= 1;
                }
                _ => break, // H_START
            },
            Tb::E => {
                steps.push((CigarOp::Ins(1), core.ext));
                if (cell >> 2) & 1 == 0 {
                    state = Tb::H;
                }
                i -= 1;
            }
            Tb::F => {
                steps.push((CigarOp::Del(1), core.ext));
                if (cell >> 4) & 1 == 0 {
                    state = Tb::H;
                }
                j -= 1;
            }
        }
    }
    steps.reverse();

    // End trimming (mm2/bwa convention): terminal runs whose cumulative
    // delta stays <= 0 become soft-clip.
    let n_steps = steps.len();
    let mut lead = 0usize;
    {
        let mut acc = 0i32;
        for (k, (_, d)) in steps.iter().enumerate() {
            acc += d;
            if acc <= 0 {
                lead = k + 1;
            }
        }
    }
    let mut trail = 0usize;
    {
        let mut acc = 0i32;
        for (k, (_, d)) in steps.iter().enumerate().rev() {
            acc += d;
            if acc <= 0 {
                trail = n_steps - k;
            }
        }
    }
    if n_steps > 0 && lead + trail >= n_steps {
        // Degenerate all-trimmed case: keep the middle step.
        lead = n_steps / 2;
        trail = n_steps - lead - 1;
    }
    let kept = &steps[lead..n_steps - trail];
    let (mut read_start, mut ref_start) = (i.max(0) as u32, j.max(0) as u32);
    let lead_read = steps[..lead]
        .iter()
        .filter(|(op, _)| matches!(op, CigarOp::Match(_) | CigarOp::Ins(_)))
        .count() as u32;
    let lead_ref = steps[..lead]
        .iter()
        .filter(|(op, _)| matches!(op, CigarOp::Match(_) | CigarOp::Del(_)))
        .count() as u32;
    read_start += lead_read;
    ref_start += lead_ref;

    // Rebuild coordinates/CIGAR/score from the kept steps.
    let mut score = 0i32;
    let mut ops: Vec<CigarOp> = Vec::new();
    let mut rm = 0u32;
    for (op, d) in kept {
        score += d;
        push_op(&mut ops, *op);
        match op {
            CigarOp::Match(_) => {
                rm += 1;
            }
            CigarOp::Ins(_) => rm += 1,
            CigarOp::Del(_) => {}
            CigarOp::SoftClip(_) | CigarOp::RefSkip(_) => {}
        }
    }
    let read_end = read_start + rm;
    let mut cigar = Vec::new();
    if read_start > 0 {
        cigar.push(CigarOp::SoftClip(read_start));
    }
    for op in ops {
        push_op(&mut cigar, op);
    }
    if read_end < m as u32 {
        push_op(&mut cigar, CigarOp::SoftClip(m as u32 - read_end));
    }

    (
        Extension {
            read_start,
            read_end,
            ref_start,
            cigar,
            score,
        },
        touched,
    )
}
