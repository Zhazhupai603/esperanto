//! Affine (Gotoh) banded extension with editing-aware substitution.
//!
//! The DP is Smith-Waterman-shaped: an alignment may start (score 0) and end
//! (global max) anywhere, with unaligned read prefixes/suffixes written as
//! soft clips — CIGAR `M/I/S` ops consume exactly the read length.
//!
//! Two execution paths produce **identical scores and CIGARs** on any input:
//! Packed i16 accumulation with a substitution LUT runs when every
//! byte of both sequences is one of {A,C,G,T,N} (uppercase); any other byte
//! (lowercase, ambiguity codes) routes to the scalar byte-wise
//! i32 path. Band scheduling is exactly two attempts — half-width 30, then
//! the full matrix — with the second attempt gated on the first attempt's
//! score falling below 85% of the perfect read score (see `extend_hint`).
//! The full-matrix packed case and the unit-slope banded attempt both
//! dispatch to AVX2 anti-diagonal kernels (`avx2_full` / `avx2_band`) that
//! must stay cell-for-cell identical to the scalar i16 kernel, with a
//! runtime feature gate and scalar fallback.

/// Which A-to-I substitutions an alignment scores as free.
///
/// A-to-I editing on one transcript is single-stranded: read against the
/// forward genome a plus-strand transcript's edits appear only as ref A/read
/// G, a minus-strand transcript's only as ref T/read C, for both mates. **At
/// most one of the two can therefore be free in any one alignment.**
///
/// Exempting both at once collapses the four-letter alphabet into two letters:
/// free A↔G merges the purines, free T↔C merges the pyrimidines. A>G and T>C
/// are exactly the *transitions*, which outnumber transversions ~2:1 in genomic
/// divergence, so under a both-directions exemption a paralogous Alu copy is
/// indistinguishable from an edited version of its sibling — the aligner has
/// deleted the very cost that separates them. That is a specificity loss paid
/// on every read, not only on edited ones.
///
/// So the exemption is per-direction, and a read that needs it is aligned once
/// under `Ag` and once under `Tc`; the better-scoring pass wins. This mirrors
/// Porath et al. 2014 (Nat Commun 5:4726), where each hyper-editing index
/// applies exactly one three-letter transform and the reverse transforms serve
/// as negative controls.
///
/// This type exists so the exemption is defined once. Every scoring path
/// routes through [`EditDir::exempt`] rather than repeating the base test.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EditDir {
    /// No exemption: every substitution costs `mismatch`.
    Off,
    /// Plus-strand transcript: ref A/read G scores 0, ref T/read C does not.
    Ag,
    /// Minus-strand transcript: ref T/read C scores 0, ref A/read G does not.
    Tc,
}

impl EditDir {
    /// Is editing-aware scoring active?
    #[inline]
    pub fn is_on(self) -> bool {
        self != EditDir::Off
    }

    /// The exemption implying transcript strand `plus`.
    #[inline]
    pub fn from_plus(plus: bool) -> Self {
        if plus {
            EditDir::Ag
        } else {
            EditDir::Tc
        }
    }

    /// Does this (ref, read) base pair score 0?
    #[inline]
    pub fn exempt(self, ref_b: u8, read_b: u8) -> bool {
        matches!(
            (self, exempt_dir(ref_b, read_b)),
            (EditDir::Ag, Some(true)) | (EditDir::Tc, Some(false))
        )
    }
}

/// Which transcript strand an exempted substitution implies, or `None` if the
/// pair is an ordinary mismatch. The single definition of the exemption: every
/// scoring and counting path routes here, so the rule cannot drift between
/// them.
///
/// `true` = plus-strand transcript (ref A / read G), `false` = minus-strand
/// (ref T / read C, the reverse-strand transcript's A-to-I seen forward).
#[inline]
pub fn exempt_dir(ref_b: u8, read_b: u8) -> Option<bool> {
    match (ref_b, read_b) {
        (b'A', b'G') => Some(true),
        (b'T', b'C') => Some(false),
        _ => None,
    }
}

/// Extension scoring parameters (match 2, mismatch −4, gap open 4,
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
    /// Which A-to-I substitutions score 0 (see [`EditDir`]).
    pub edit_dir: EditDir,
}

impl Default for ExtendParams {
    fn default() -> ExtendParams {
        ExtendParams {
            match_score: 2,
            mismatch: -4,
            gap_open: 4,
            gap_ext: 2,
            edit_dir: EditDir::Off,
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
    /// Build the LUT from parameters; the exempted pairs are whichever
    /// [`EditDir`] declares free (A>I editing appears as G in the read).
    pub fn new(params: &ExtendParams) -> SubstLut {
        let m = params.match_score as i16;
        let mm = params.mismatch as i16;
        let mut scores = [mm; 25];
        for c in 0..4 {
            scores[c * 5 + c] = m;
        }
        // LUT layout is scores[read_code*5 + ref_code].
        if params.edit_dir.exempt(b'A', b'G') {
            scores[10] = 0; // read G(2), ref A(0): 2*5+0
        }
        if params.edit_dir.exempt(b'T', b'C') {
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

/// First strict band-edge touch. The `touched` flag controls escalation;
/// the coordinates describe the boundary condition and support kernel
/// equivalence checks. A clamped edge lies on the reference window boundary.
#[derive(Clone, Copy, Debug, Default)]
pub struct TouchRec {
    /// Any row maximum sat on a strict band edge.
    pub touched: bool,
    /// First touched row (1-based read index); 0 when untouched.
    pub row: i64,
    /// Left edge (true) or right edge (false).
    pub left: bool,
    /// Touched edge was the window boundary, not the band.
    pub clamped: bool,
}

/// Reusable per-thread DP buffers (no semantic effect).
pub struct ExtendBuffer {
    packed: BandBuf<i16>,
    scalar: BandBuf<i32>,
}

impl Default for ExtendBuffer {
    fn default() -> ExtendBuffer {
        ExtendBuffer {
            packed: BandBuf::new(),
            scalar: BandBuf::new(),
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
/// other byte (lowercase, ambiguity codes) routes to the byte-wise scalar
/// path. Both paths score every byte pair identically, so the produced
/// `Extension` is the same either way.
fn run_dispatch(
    read: &[u8],
    ref_window: &[u8],
    params: &ExtendParams,
    buf: &mut ExtendBuffer,
    band_spec: Option<(i64, i64, i64, i64)>,
) -> (Extension, TouchRec) {
    let core = CoreParams::new(params);
    let pure = read.iter().all(|&b| pure_byte(b)) && ref_window.iter().all(|&b| pure_byte(b));
    // Full-matrix packed calls run the AVX2 anti-diagonal kernel when the
    // CPU (and the ESP_NO_SIMD escape hatch) allow: it reproduces
    // `run_band::<i16>` cell-for-cell (see `avx2_full`). Banded attempts
    // run the banded AVX2 kernel (`avx2_band`) under the same gate when the
    // band is the unit-slope form (num == den == 1) every call site uses;
    // other slopes, the scalar i32 path, and tiny windows stay scalar.
    #[cfg(target_arch = "x86_64")]
    if pure
        && band_spec.is_none()
        && read.len() >= 8
        && ref_window.len() >= 8
        && avx2_full::available()
    {
        return avx2_full::run(read, ref_window, &core, &mut buf.packed);
    }
    #[cfg(target_arch = "x86_64")]
    if pure
        && band_spec.is_some()
        && read.len() >= 8
        && ref_window.len() >= 8
        && avx2_full::available()
    {
        if let Some((offset, num, den, half)) = band_spec {
            if num == 1 && den.max(1) == 1 {
                return avx2_band::run(
                    read, ref_window, &core, &mut buf.packed, offset, num, den, half,
                );
            }
        }
    }
    if pure {
        run_band(read, ref_window, &core, &mut buf.packed, band_spec)
    } else {
        run_band(read, ref_window, &core, &mut buf.scalar, band_spec)
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
/// diagonal, and attempt 2 recomputes the full matrix and is returned
/// unconditionally. A band that already covers the whole window
/// (`30 ≥ ref_window.len()`) is accepted directly as the full matrix.
///
/// Attempt 1 is kept only when it is safe to keep, and that takes both of
/// the DP's own signals:
///
///   score      `score × 20 ≥ perfect × 17` (with `perfect = read_len ×
///              match_score`, i.e. at least 85% of the perfect read score).
///              A badly scoring in-band path may have a much better
///              out-of-band alternative.
///   band edge  no row's best `H` sat on a strict band edge. This one is
///              what makes the score gate sound rather than a guess: while
///              the maximum stays strictly inside the band, widening it
///              cannot change the answer, so attempt 1 *is* the full
///              matrix's answer. Once the maximum sits on the edge, the
///              optimal path wanted a diagonal the band does not contain and
///              attempt 1 is suboptimal by construction.
///
/// Score alone is not enough, and the gap between the two is the whole
/// defect: at 85% of perfect a 150 bp read may carry seven mismatches, so a
/// read whose true alignment needs an indel shifting the diagonal by more
/// than 30 -- or a junction the hint missed -- is returned in-band, with a
/// clean-looking score and the wrong CIGAR. Measured on the real sample, the
/// reads esperanto and STAR align *identically* score 10.46x on the
/// strand-symmetric editing control (`mapbench/eindex.py`), matching STAR's
/// own 10.22x, while the 5.3% they shape differently score 3.96x and carry
/// 60% of chr1's control-channel sites. That partition is not reachable by
/// any per-read filter tried against a size-matched random control -- soft
/// clip or indel presence, `ext_margin`, multi-mapping, read coverage, and
/// cross-contig placement all came back null -- because the reads are
/// individually clean and only their shape is wrong. `touched` is the signal
/// that says so from inside the DP, and it was already being computed and
/// discarded here.
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
    let (ext, touch) = run_dispatch(
        read,
        ref_window,
        params,
        buf,
        Some((hint.offset, hint.num as i64, hint.den.max(1) as i64, INITIAL_BAND)),
    );
    // Escalation: the band held the maximum on its edge, so a wider band can
    // still improve the path; or the score is under 85% of perfect.
    let perfect = read.len() as i32 * params.match_score;
    if touch.touched || ext.score * 20 < perfect * 17 {
        return run_dispatch(read, ref_window, params, buf, None).0;
    }
    ext
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
            open_ext: -params.gap_open, // first gap base pays gap_open only
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
    // Long-window E/F chains may sink below -30000. Exact sentinel
    // equality and the `.max(NEG)` clamp merge these unreachable states
    // while preserving deterministic ties in the gap recurrences.
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

/// One DP row's direction bytes, as a span into `BandBuf::rows_flat`
/// (flat trace storage: one allocation per DP call, not one per row).
struct RowTrace {
    lo: i64,
    off: usize,
    len: usize,
}

struct BandBuf<T: Acc> {
    h_prev: Vec<T>,
    h_cur: Vec<T>,
    e_prev: Vec<T>,
    e_cur: Vec<T>,
    f_cur: Vec<T>,
    rows: Vec<RowTrace>,
    rows_flat: Vec<u8>,
    read_codes: Vec<usize>,
    ref_codes: Vec<usize>,
    // AVX2 full-matrix scratch (see `avx2_full`); unused by the scalar
    // paths, concrete-typed so the generic buffer carries them for free.
    avx_h: [Vec<i16>; 3],
    avx_e: [Vec<i16>; 2],
    avx_f: [Vec<i16>; 2],
    avx_row_max: Vec<i16>,
    avx_row_arg: Vec<i32>,
    avx_row_h1: Vec<i16>,
    avx_read: Vec<i16>,
    avx_refrev: Vec<i16>,
    // AVX2 banded scratch (see `avx2_band`): per-row window bounds,
    // anti-diagonal extent, and flat-scatter base; row_hl/hr are the row's
    // first/last in-band H for the touch record.
    avx_band_lo: Vec<i64>,
    avx_band_hi: Vec<i64>,
    avx_band_alo: Vec<i64>,
    avx_band_ahi: Vec<i64>,
    avx_band_base: Vec<i64>,
    avx_row_hl: Vec<i16>,
    avx_row_hr: Vec<i16>,
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
            rows_flat: Vec::new(),
            read_codes: Vec::new(),
            ref_codes: Vec::new(),
            avx_h: std::array::from_fn(|_| Vec::new()),
            avx_e: std::array::from_fn(|_| Vec::new()),
            avx_f: std::array::from_fn(|_| Vec::new()),
            avx_row_max: Vec::new(),
            avx_row_arg: Vec::new(),
            avx_row_h1: Vec::new(),
            avx_read: Vec::new(),
            avx_refrev: Vec::new(),
            avx_band_lo: Vec::new(),
            avx_band_hi: Vec::new(),
            avx_band_alo: Vec::new(),
            avx_band_ahi: Vec::new(),
            avx_band_base: Vec::new(),
            avx_row_hl: Vec::new(),
            avx_row_hr: Vec::new(),
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
) -> (Extension, TouchRec) {
    // Identical recurrence, boundary, tie preference (diag > E > F > start),
    // zero floor, traceback, and end trimming as the reference
    // implementation this replaced. Implementation changes only:
    // pre-encoded bases (no per-cell base_code call), one flat traceback
    // allocation per call (no per-row Vec), no per-row clearing (every read
    // of the previous row is guarded by that row's written window), and the
    // band center without per-row division (num == den == 1 at every call
    // site, asserted).
    let m = read.len() as i64;
    let n = rf.len() as i64;
    if m == 0 || n == 0 {
        return (empty_extension(m as u32), TouchRec::default());
    }

    let (offset, num, den, half) = match band_spec {
        Some((o, nm, dn, h)) => (o, nm, dn.max(1), h.max(0)),
        // Full window: half = n covers every diagonal whenever m <= n.
        None => {
            debug_assert!(m <= n);
            (0, 1, 1, n)
        }
    };
    debug_assert!(num == 1 && den == 1);
    let width = (2 * half + 1) as usize;
    if buf.h_prev.len() != width {
        buf.h_prev = vec![T::NEG; width];
        buf.h_cur = vec![T::NEG; width];
        buf.e_prev = vec![T::NEG; width];
        buf.e_cur = vec![T::NEG; width];
        buf.f_cur = vec![T::NEG; width];
    }
    // Once per CALL, not per row: every per-row read is guarded by the
    // written window, but row 1 reads E/F of the never-written row 0, and a
    // reused buffer otherwise carries the previous call's last row.
    buf.h_prev.fill(T::NEG);
    buf.h_cur.fill(T::NEG);
    buf.e_prev.fill(T::NEG);
    buf.e_cur.fill(T::NEG);
    buf.f_cur.fill(T::NEG);
    buf.rows.clear();
    buf.rows_flat.clear();

    // Band center at row i: offset + i (num/den is 1/1 at every call site).
    // Rows (i>=1) clamp lo at 1 — column 0 is the pinned H=0 boundary and
    // never enters the DP; row 0 itself is zeroed from col 0.
    let center0 = offset;
    let (lo0, hi0) = ((center0 - half).max(0), (center0 + half).min(n));

    buf.read_codes.clear();
    buf.read_codes.extend(read.iter().map(|&b| base_code(b)));
    buf.ref_codes.clear();
    buf.ref_codes.extend(rf.iter().map(|&b| base_code(b)));

    // Row 0: H = 0 within the window (free reference prefix).
    {
        if lo0 <= hi0 {
            let off = buf.rows_flat.len();
            let len = (hi0 - lo0 + 1) as usize;
            buf.rows_flat.resize(off + len, H_START);
            for c in 0..len {
                buf.h_prev[c] = T::from_i32(0);
            }
            buf.rows.push(RowTrace { lo: lo0, off, len });
        } else {
            buf.rows.push(RowTrace {
                lo: lo0,
                off: buf.rows_flat.len(),
                len: 0,
            });
        }
    }

    let open_ext = T::from_i32(core.open_ext);
    let ext_p = T::from_i32(core.ext);
    let zero = T::from_i32(0);

    let mut best_score = zero;
    let mut best_i = 0i64;
    let mut best_j = 0i64;
    let mut touch = TouchRec::default();
    // Previous row's window, tracked across iterations (no recompute).
    let (mut plo, mut phi) = (lo0, hi0);
    for i in 1..=m {
        let center = offset + i;
        let c = center.clamp(1, n);
        let (lo, hi) = ((c - half).max(1), (c + half).min(n));
        if lo > hi {
            buf.rows.push(RowTrace {
                lo,
                off: buf.rows_flat.len(),
                len: 0,
            });
            plo = lo;
            phi = hi;
            continue;
        }
        let off = buf.rows_flat.len();
        let len = (hi - lo + 1) as usize;
        buf.rows_flat.resize(off + len, 0);

        let rc = buf.read_codes[(i - 1) as usize];
        let mut row_max = T::NEG;
        let mut row_max_cell = 0i64;


        for jj in lo..=hi {
            let cidx = (jj - lo) as usize;

            // E: vertical (insertion): from H[i-1][j] or E[i-1][j].
            // Out-of-band H reads see the zero boundary; E reads stay NEG.
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

            // F: horizontal (deletion): from H[i][j-1] or F[i][j-1]; the
            // left margin/boundary cell reads H=0; F stays NEG.
            let (h_left, f_left) = if jj > lo {
                (buf.h_cur[cidx - 1], buf.f_cur[cidx - 1])
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
            // col 0 is the pinned 0 boundary, margin cells read 0.
            let h_prev_idx = jj - 1;
            let h_diag = if h_prev_idx >= plo && h_prev_idx <= phi {
                buf.h_prev[(h_prev_idx - plo) as usize]
            } else {
                zero
            };
            let s = core.lut.score(rc, buf.ref_codes[(jj - 1) as usize]);
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

            buf.h_cur[cidx] = h_val;
            buf.e_cur[cidx] = e_val;
            buf.f_cur[cidx] = f_val;
            buf.rows_flat[off + cidx] = hdir | (edir << 2) | (fdir << 4);

            if T::better(h_val, row_max) {
                row_max = h_val;
                row_max_cell = cidx as i64;
            }
        }

        // Strict band edge touch (edge exists and holds the row maximum).
        // Attribution showed 88% of first touches are first-decile tie
        // artifacts and 88% of escalations reproduce the band score -- but
        // suppressing tie touches (and/or requiring `lo > 1`) is NOT
        // byte-identical on real data: the remaining escalations are
        // load-bearing, the full matrix finds genuinely better out-of-band
        // paths for those reads. Escalation semantics stay exactly as they
        // were; only the touch is recorded for profiling.
        let left_strict = lo > 0;
        let right_strict = hi < n;
        if !touch.touched {
            if left_strict && T::same(buf.h_cur[0], row_max) {
                touch = TouchRec { touched: true, row: i, left: true, clamped: lo <= 1 };
            } else if right_strict && T::same(buf.h_cur[len - 1], row_max) {
                touch = TouchRec { touched: true, row: i, left: false, clamped: hi >= n };
            }
        }

        if T::better(row_max, best_score) {
            best_score = row_max;
            best_i = i;
            best_j = lo + row_max_cell;
        }

        buf.rows.push(RowTrace { lo, off, len });
        std::mem::swap(&mut buf.h_prev, &mut buf.h_cur);
        std::mem::swap(&mut buf.e_prev, &mut buf.e_cur);
        // No per-row clearing: every read of the previous row is guarded by
        // that row's written window (in_prev / diagonal), and the current
        // row is written before any read within it (f_left reads cidx-1).
        // The one necessary clear is once per call, at the top.
        plo = lo;
        phi = hi;
    }

    let ext = finish_extension(
        m,
        core,
        &buf.rows,
        &buf.rows_flat,
        &buf.read_codes,
        &buf.ref_codes,
        best_i,
        best_j,
    );
    (ext, touch)
}

/// Shared post-DP tail of both kernels: state-machine traceback (H/E/F) from
/// `(best_i, best_j)` recording per-step deltas for end trimming, mm2/bwa
/// end trimming, and the Extension rebuild. Scalar and AVX2 full-matrix
/// kernels call this on identically-filled trace storage, so the produced
/// Extension is byte-identical by construction.
///
/// Cell bits: low2 = hdir; bit2 = E-from-extension; bit4 = F-from-extension.
/// Gap-open cost lands on the first step of a run.
// Traceback walks rows/columns one cell at a time; the argument list mirrors
// the DP state it consumes and splitting it would only obscure the flow.
#[allow(clippy::too_many_arguments)]
fn finish_extension(
    m: i64,
    core: &CoreParams,
    rows: &[RowTrace],
    rows_flat: &[u8],
    read_codes: &[usize],
    ref_codes: &[usize],
    best_i: i64,
    best_j: i64,
) -> Extension {
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
        let row = &rows[i as usize];
        if j < row.lo || j >= row.lo + row.len as i64 {
            break;
        }
        let cell = rows_flat[row.off + (j - row.lo) as usize];
        match state {
            Tb::H => match cell & 3 {
                H_DIAG => {
                    let s = core.lut.score(
                        read_codes[(i - 1) as usize],
                        ref_codes[(j - 1) as usize],
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

    Extension {
        read_start,
        read_end,
        ref_start,
        cigar,
        score,
    }
}

/// AVX2 anti-diagonal kernel for the full-matrix i16 path (`band_spec == None`
/// at [`run_dispatch`]).
///
/// Per 16-lane block it evaluates the exact scalar recurrence — sentinel
/// `add` (exact-equality NEG + `max(NEG)` clamp), strict `>` selects, tie
/// cascade F → E → diag with the zero-floor rewrite to `H_START` — over the
/// cells of one anti-diagonal d = i + j, whose only dependencies sit on
/// d−1 (E reads (i−1, j)'s H and E; F reads (i, j−1)'s H and F) and d−2
/// (the diagonal H of (i−1, j−1)). Boundary cells (row 0, column 0) carry
/// H = 0 and E = F = NEG, exactly the scalar margins, so no lane needs
/// special casing inside the vector body.
///
/// Substitution scores come from the same 5×5 LUT: match = `scores[0]`,
/// mismatch = `scores[1]`, and the two editing-exempt cells `scores[10]`
/// (G/A) and `scores[8]` (C/T) blended on their (read, ref) code equality —
/// a pure function of the LUT, so any parameter set scores identically.
///
/// The global traceback start is re-derived per row after the sweep
/// (`avx_row_max`, first strict-improving column kept), which is the first
/// cell in row-major order attaining the strict maximum — the same cell
/// the scalar loop keeps. The touch record likewise: first row whose
/// H(i, 1) equals its row max (left edge, clamped — the full matrix's lo=1
/// is always strict, its hi=n never is). `dp_cells` is bumped by
/// `m × (2n + 1)`, the scalar full-matrix width.
#[cfg(target_arch = "x86_64")]
mod avx2_full {
    use super::{
        base_code, empty_extension, finish_extension, run_band, Acc, BandBuf, CoreParams,
        Extension, RowTrace, TouchRec, H_START,
    };
    use std::arch::x86_64::*;
    use std::sync::OnceLock;

    /// i16 sentinel (see `impl Acc for i16`; exact equality + clamp).
    const NEG: i16 = <i16 as Acc>::NEG;

    /// Runtime gate: AVX2 detected, unless `ESP_NO_SIMD=1` masks it (debug
    /// escape hatch forcing the scalar kernel).
    pub fn available() -> bool {
        static OK: OnceLock<bool> = OnceLock::new();
        *OK.get_or_init(|| {
            if std::env::var_os("ESP_NO_SIMD").is_some_and(|v| v == "1") {
                return false;
            }
            is_x86_feature_detected!("avx2")
        })
    }

    /// Full-matrix packed DP on AVX2. Fills `buf.rows` / `buf.rows_flat`
    /// exactly as `run_band::<i16>(…, None)` does and shares its traceback
    /// tail (`finish_extension`); the returned `Extension` and `TouchRec`
    /// are field-identical to the scalar kernel's.
    pub fn run(
        read: &[u8],
        rf: &[u8],
        core: &CoreParams,
        buf: &mut BandBuf<i16>,
    ) -> (Extension, TouchRec) {
        let m = read.len() as i64;
        let n = rf.len() as i64;
        if m == 0 || n == 0 {
            return (empty_extension(m as u32), TouchRec::default());
        }
        debug_assert!(m <= n);

        // LUT-derived substitution constants. The vector `add` assumes the
        // constant operand is never the sentinel; default scoring never is,
        // and anything exotic falls back to the scalar kernel (which bumps
        // dp_cells itself, so the census below stays exact).
        let s_diag = core.lut.score(0, 0);
        let s_mm = core.lut.score(0, 1);
        let s_ga = core.lut.score(2, 0); // read G, ref A (LUT index 10)
        let s_tc = core.lut.score(1, 3); // read C, ref T (LUT index 8)
        let open_i = core.open_ext as i16;
        let ext_i = core.ext as i16;
        if [s_diag, s_mm, s_ga, s_tc, open_i, ext_i].contains(&NEG) {
            return run_band(read, rf, core, buf, None);
        }

        buf.read_codes.clear();
        buf.read_codes.extend(read.iter().map(|&b| base_code(b)));
        buf.ref_codes.clear();
        buf.ref_codes.extend(rf.iter().map(|&b| base_code(b)));

        // Trace storage laid out exactly as the scalar full matrix: row 0
        // spans columns 0..=n (H_START), rows 1..=m span 1..=n.
        buf.rows.clear();
        buf.rows_flat.clear();
        buf.rows_flat.resize((n + 1 + m * n) as usize, 0);
        buf.rows_flat[0..=(n as usize)].fill(H_START);
        buf.rows.push(RowTrace {
            lo: 0,
            off: 0,
            len: (n + 1) as usize,
        });
        for i in 1..=m {
            buf.rows.push(RowTrace {
                lo: 1,
                off: ((n + 1) + (i - 1) * n) as usize,
                len: n as usize,
            });
        }

        let cap = m as usize + 33; // vector stores run one block past ihi
        for v in buf.avx_h.iter_mut() {
            v.resize(cap, 0);
        }
        for v in buf.avx_e.iter_mut() {
            v.resize(cap, 0);
        }
        for v in buf.avx_f.iter_mut() {
            v.resize(cap, 0);
        }
        buf.avx_row_max.clear();
        buf.avx_row_max.resize(m as usize + 1, NEG);
        buf.avx_row_arg.clear();
        buf.avx_row_arg.resize(m as usize + 1, 1);
        buf.avx_row_h1.clear();
        buf.avx_row_h1.resize(m as usize + 1, NEG);
        buf.avx_read.clear();
        buf.avx_read.extend(buf.read_codes.iter().map(|&c| c as i16));
        buf.avx_read.resize(m as usize + 16, 0);
        // Ref codes reversed: lane i (ascending) reads ref column d − i, so
        // `refrev[q] = ref_codes[n − 1 − q]` turns the descending gather
        // into the contiguous load at `i + n − d`.
        buf.avx_refrev.clear();
        buf.avx_refrev.extend(
            (0..n as usize).map(|q| buf.ref_codes[(n as usize) - 1 - q] as i16),
        );
        buf.avx_refrev.resize(n as usize + 16, 0);

        let hp = [
            buf.avx_h[0].as_mut_ptr(),
            buf.avx_h[1].as_mut_ptr(),
            buf.avx_h[2].as_mut_ptr(),
        ];
        let ep = [buf.avx_e[0].as_mut_ptr(), buf.avx_e[1].as_mut_ptr()];
        let fp = [buf.avx_f[0].as_mut_ptr(), buf.avx_f[1].as_mut_ptr()];
        // SAFETY: AVX2 was checked by `available()` at dispatch; all
        // pointers stay inside the resized scratch (cap ≥ m + 33, code
        // arrays padded 16, rows_flat length n + 1 + mn covers every
        // scatter index i·(n−1) + d), and the eleven buffers never alias.
        unsafe {
            dp_diagonals(
                m,
                n,
                buf.avx_read.as_ptr(),
                buf.avx_refrev.as_ptr(),
                hp,
                ep,
                fp,
                buf.rows_flat.as_mut_ptr(),
                buf.avx_row_max.as_mut_ptr(),
                buf.avx_row_arg.as_mut_ptr(),
                buf.avx_row_h1.as_mut_ptr(),
                open_i,
                ext_i,
                s_diag,
                s_mm,
                s_ga,
                s_tc,
            );
        }

        // Row-major first strict maximum + first left-edge tie: rows
        // ascending, strict `>` — the scalar loop's exact selection order
        // (equal-value cells never displace an earlier row/column winner).
        let mut best_score = 0i16;
        let mut best_i = 0i64;
        let mut best_j = 0i64;
        let mut touch = TouchRec::default();
        for i in 1..=m {
            let rm = buf.avx_row_max[i as usize];
            if !touch.touched && buf.avx_row_h1[i as usize] == rm {
                // Full matrix: lo = 1 is a strict left edge and clamped,
                // hi = n is never a strict right edge.
                touch = TouchRec {
                    touched: true,
                    row: i,
                    left: true,
                    clamped: true,
                };
            }
            if rm > best_score {
                best_score = rm;
                best_i = i;
                best_j = buf.avx_row_arg[i as usize] as i64;
            }
        }
        let ext = finish_extension(
            m,
            core,
            &buf.rows,
            &buf.rows_flat,
            &buf.read_codes,
            &buf.ref_codes,
            best_i,
            best_j,
        );
        (ext, touch)
    }

    /// The anti-diagonal sweep. `h`/`e`/`f` are per-diagonal buffers indexed
    /// by row i; roles rotate one step per diagonal (H needs d−2, d−1, d;
    /// E and F need d−1, d). Interior lanes are computed in 16-wide blocks;
    /// lanes past `ihi` compute garbage that is stored only where a later
    /// diagonal never reads (every read of diagonal d touches valid cells of
    /// d−1/d−2 only).
    #[target_feature(enable = "avx2")]
    #[allow(clippy::too_many_arguments)]
    unsafe fn dp_diagonals(
        m: i64,
        n: i64,
        read_i16: *const i16,
        refrev: *const i16,
        h: [*mut i16; 3],
        e: [*mut i16; 2],
        f: [*mut i16; 2],
        rows_flat: *mut u8,
        row_max: *mut i16,
        row_arg: *mut i32,
        row_h1: *mut i16,
        open_ext: i16,
        ext: i16,
        s_diag: i16,
        s_mm: i16,
        s_ga: i16,
        s_tc: i16,
    ) {
        let neg_v = _mm256_set1_epi16(NEG);
        let zero_v = _mm256_set1_epi16(0);
        let ones_v = _mm256_set1_epi16(-1);
        let one_v = _mm256_set1_epi16(1);
        let two_v = _mm256_set1_epi16(2);
        let three_v = _mm256_set1_epi16(3);
        let four_v = _mm256_set1_epi16(4);
        let open_v = _mm256_set1_epi16(open_ext);
        let ext_v = _mm256_set1_epi16(ext);
        let diag_v = _mm256_set1_epi16(s_diag);
        let mm_v = _mm256_set1_epi16(s_mm);
        let ga_v = _mm256_set1_epi16(s_ga);
        let tc_v = _mm256_set1_epi16(s_tc);
        let code_g = _mm256_set1_epi16(2);
        let code_a = _mm256_set1_epi16(0);
        let code_c = _mm256_set1_epi16(1);
        let code_t = _mm256_set1_epi16(3);

        // `Acc::add` with a constant non-NEG operand: NEG propagates,
        // otherwise max(x + c, NEG) — wrapping add, exactly the release
        // scalar (real inputs never wrap: H is zero-floored and small).
        macro_rules! addc {
            ($x:expr, $c:expr) => {{
                let isn = _mm256_cmpeq_epi16($x, neg_v);
                let r = _mm256_max_epi16(_mm256_add_epi16($x, $c), neg_v);
                _mm256_blendv_epi8(r, neg_v, isn)
            }};
        }

        let stride = (n - 1) as isize;
        let mut h2 = 0usize; // roles: h2 = d−2, h1 = d−1, h0 = d
        let mut h1 = 1usize;
        let mut h0 = 2usize;
        let mut e1 = 0usize; // e1 = d−1, e0 = d
        let mut e0 = 1usize;
        let mut f1 = 0usize;
        let mut f0 = 1usize;

        for d in 0..=(m + n) {
            let h2p = h[h2];
            let h1p = h[h1];
            let h0p = h[h0];
            let e1p = e[e1];
            let e0p = e[e0];
            let f1p = f[f1];
            let f0p = f[f0];
            let ilo = 1i64.max(d - n);
            let ihi = m.min(d - 1);
            let mut v = ilo;
            while v <= ihi {
                // Lane i covers cell (i, d − i): read code [i − 1], ref code
                // [d − i − 1] via the reversed array at [i + n − d].
                let rc = _mm256_loadu_si256(read_i16.add((v - 1) as usize) as *const __m256i);
                let rfc =
                    _mm256_loadu_si256(refrev.add((v + (n - d)) as usize) as *const __m256i);
                // Substitution LUT as blends: match (equal codes < 4) beats
                // the exempt cells (G/A, C/T), which beat mismatch.
                let eq = _mm256_cmpeq_epi16(rc, rfc);
                let is_match = _mm256_and_si256(eq, _mm256_cmpgt_epi16(four_v, rc));
                let ga = _mm256_and_si256(
                    _mm256_cmpeq_epi16(rc, code_g),
                    _mm256_cmpeq_epi16(rfc, code_a),
                );
                let tc = _mm256_and_si256(
                    _mm256_cmpeq_epi16(rc, code_c),
                    _mm256_cmpeq_epi16(rfc, code_t),
                );
                let mut s = _mm256_blendv_epi8(mm_v, ga_v, ga);
                s = _mm256_blendv_epi8(s, tc_v, tc);
                s = _mm256_blendv_epi8(s, diag_v, is_match);

                // E (vertical): from H[i−1][j] / E[i−1][j], one diagonal up.
                let h_up = _mm256_loadu_si256(h1p.add((v - 1) as usize) as *const __m256i);
                let e_up = _mm256_loadu_si256(e1p.add((v - 1) as usize) as *const __m256i);
                let e_open = addc!(h_up, open_v);
                let e_ext = addc!(e_up, ext_v);
                let egt = _mm256_cmpgt_epi16(e_ext, e_open);
                let e_val = _mm256_blendv_epi8(e_open, e_ext, egt);
                let edir = _mm256_and_si256(egt, one_v);

                // F (horizontal): from H[i][j−1] / F[i][j−1], same diagonal.
                let h_left = _mm256_loadu_si256(h1p.add(v as usize) as *const __m256i);
                let f_left = _mm256_loadu_si256(f1p.add(v as usize) as *const __m256i);
                let f_open = addc!(h_left, open_v);
                let f_ext = addc!(f_left, ext_v);
                let fgt = _mm256_cmpgt_epi16(f_ext, f_open);
                let f_val = _mm256_blendv_epi8(f_open, f_ext, fgt);
                let fdir = _mm256_and_si256(fgt, one_v);

                // Diagonal: H[i−1][j−1] + s.
                let h_diag = _mm256_loadu_si256(h2p.add((v - 1) as usize) as *const __m256i);
                let diag = addc!(h_diag, s);

                // H: zero floor, then F → E → diag, each winning only on
                // ≥ against the running value (later source wins ties),
                // then the ==0 rewrite to H_START. Lane-for-lane the scalar
                // cascade.
                let mut hv = zero_v;
                let mut hd = three_v;
                let fmask = _mm256_andnot_si256(
                    _mm256_or_si256(_mm256_cmpeq_epi16(f_val, neg_v), _mm256_cmpgt_epi16(hv, f_val)),
                    ones_v,
                );
                hv = _mm256_blendv_epi8(hv, f_val, fmask);
                hd = _mm256_blendv_epi8(hd, two_v, fmask);
                let emask = _mm256_andnot_si256(
                    _mm256_or_si256(_mm256_cmpeq_epi16(e_val, neg_v), _mm256_cmpgt_epi16(hv, e_val)),
                    ones_v,
                );
                hv = _mm256_blendv_epi8(hv, e_val, emask);
                hd = _mm256_blendv_epi8(hd, one_v, emask);
                let dmask = _mm256_andnot_si256(
                    _mm256_or_si256(_mm256_cmpeq_epi16(diag, neg_v), _mm256_cmpgt_epi16(hv, diag)),
                    ones_v,
                );
                hv = _mm256_blendv_epi8(hv, diag, dmask);
                hd = _mm256_blendv_epi8(hd, zero_v, dmask);
                hd = _mm256_blendv_epi8(hd, three_v, _mm256_cmpeq_epi16(hv, zero_v));

                let dir16 = _mm256_or_si256(
                    hd,
                    _mm256_or_si256(_mm256_slli_epi16(edir, 2), _mm256_slli_epi16(fdir, 4)),
                );

                _mm256_storeu_si256(h0p.add(v as usize) as *mut __m256i, hv);
                _mm256_storeu_si256(e0p.add(v as usize) as *mut __m256i, e_val);
                _mm256_storeu_si256(f0p.add(v as usize) as *mut __m256i, f_val);

                // Per-lane tail: direction byte scatter into row-major
                // trace storage (flat index i·(n−1) + d), row max with
                // first-strict-improving column, and H(i, 1) for the touch
                // record. Cells arrive in ascending j per row across
                // diagonals, so strict `>` keeps the scalar's first-j rule.
                let mut hb = [0i16; 16];
                let mut db = [0i16; 16];
                _mm256_storeu_si256(hb.as_mut_ptr() as *mut __m256i, hv);
                _mm256_storeu_si256(db.as_mut_ptr() as *mut __m256i, dir16);
                let mut addr = rows_flat.add((v * (n - 1) + d) as usize);
                for k in 0..16usize {
                    let i = v + k as i64;
                    if i > ihi {
                        break;
                    }
                    *addr = db[k] as u8;
                    addr = addr.offset(stride);
                    let hv16 = hb[k];
                    let rm = row_max.add(i as usize);
                    if hv16 > *rm {
                        *rm = hv16;
                        *row_arg.add(i as usize) = (d - i) as i32;
                    }
                    if d - i == 1 {
                        *row_h1.add(i as usize) = hv16;
                    }
                }
                v += 16;
            }

            // Boundary cells of this diagonal (row 0 / column 0): H = 0,
            // E = F = NEG — the scalar margins. Written after the vector
            // stores so the garbage lanes they leave past ihi cannot leak.
            if d <= n {
                *h0p = 0;
                *e0p = NEG;
                *f0p = NEG;
            }
            if d <= m {
                *h0p.add(d as usize) = 0;
                *e0p.add(d as usize) = NEG;
                *f0p.add(d as usize) = NEG;
            }

            // Rotate roles.
            let t = h2;
            h2 = h1;
            h1 = h0;
            h0 = t;
            std::mem::swap(&mut e1, &mut e0);
            std::mem::swap(&mut f1, &mut f0);
        }
    }
}

/// AVX2 anti-diagonal kernel for the BANDED i16 path (`band_spec == Some`
/// with the unit-slope band `num == den == 1` at [`run_dispatch`]).
///
/// Same sweep structure as [`avx2_full`]: per 16-lane block it evaluates the
/// exact scalar recurrence over the cells of one anti-diagonal d = i + j,
/// with H/E/F rotating per-diagonal buffers indexed by row i. The band
/// changes only which cells are valid: on diagonal d the valid rows form the
/// contiguous range [blo, bhi] = [min{i : hi_i + i ≥ d}, max{i : lo_i + i ≤ d}]
/// (row windows shift right by 0/1 per row, so both extents are monotone and
/// the range is an interval), and blocks run exactly over it.
///
/// Boundary/margin semantics. The scalar kernel reads three kinds of
/// out-of-window values: E from (i−1, j) when j sits past row i−1's right
/// edge (windows only ever shift right: a fresh-right-edge cell, one per row
/// at most), F from (i, j−1) at every row's left edge, and the diagonal H
/// from (i−1, j−1) when row i−1's window did not shift (lo_i == lo_{i−1}).
/// All three read as (H = 0, E/F = NEG). The kernel materialises them the
/// same way the scalar margins read: the per-diagonal buffers are
/// initialised H = 0 / E = F = NEG once per call (row 0 and never-touched
/// positions then carry exactly the margin values — row 0's window is
/// irrelevant since position 0 is never stored to), and a per-diagonal
/// margin schedule overwrites the three possible stale/garbage positions
/// AFTER the vector body (stores run one block past bhi; a margin position
/// is by construction never a valid cell of its own diagonal, so this
/// cannot clobber real data). Rows 1..=m are never empty (lo = max(c−half, 1)
/// ≤ n and hi = min(c+half, n) ≥ 1 for c ∈ [1, n]), so only row 0 can have
/// an empty window, and it is pure boundary.
///
/// Traceback start and touch record are re-derived per row after the sweep
/// (`avx_row_max` first strict-improving column, `avx_row_hl`/`avx_row_hr`
/// the row's first/last in-band H) in a scalar pass that applies the exact
/// scalar selection order: left edge before right, value-equality tie test,
/// `clamped = lo ≤ 1 / hi ≥ n`, strict `>` against the running best.
/// `dp_cells` is bumped by `m × (2·half + 1)`, the scalar band width.
#[cfg(target_arch = "x86_64")]
mod avx2_band {
    use super::{
        base_code, empty_extension, finish_extension, run_band, Acc, BandBuf, CoreParams,
        Extension, RowTrace, TouchRec, H_START,
    };
    use std::arch::x86_64::*;

    /// i16 sentinel (see `impl Acc for i16`; exact equality + clamp).
    const NEG: i16 = <i16 as Acc>::NEG;

    /// Banded packed DP on AVX2, `band_spec = Some((offset, 1, 1, half))`.
    /// Fills `buf.rows` / `buf.rows_flat` exactly as `run_band::<i16>(…, Some(..))`
    /// does and shares its traceback tail (`finish_extension`); the returned
    /// `Extension` and `TouchRec` are field-identical to the scalar kernel's.
    #[allow(clippy::too_many_arguments)]
    pub fn run(
        read: &[u8],
        rf: &[u8],
        core: &CoreParams,
        buf: &mut BandBuf<i16>,
        offset: i64,
        num: i64,
        den: i64,
        half: i64,
    ) -> (Extension, TouchRec) {
        let m = read.len() as i64;
        let n = rf.len() as i64;
        if m == 0 || n == 0 {
            return (empty_extension(m as u32), TouchRec::default());
        }
        let den = den.max(1);
        let half = half.max(0);
        debug_assert!(num == 1 && den == 1);

        // LUT-derived substitution constants, with the same scalar fallback
        // as the full kernel (which bumps dp_cells itself).
        let s_diag = core.lut.score(0, 0);
        let s_mm = core.lut.score(0, 1);
        let s_ga = core.lut.score(2, 0); // read G, ref A (LUT index 10)
        let s_tc = core.lut.score(1, 3); // read C, ref T (LUT index 8)
        let open_i = core.open_ext as i16;
        let ext_i = core.ext as i16;
        if [s_diag, s_mm, s_ga, s_tc, open_i, ext_i].contains(&NEG) {
            return run_band(read, rf, core, buf, Some((offset, num, den, half)));
        }

        buf.read_codes.clear();
        buf.read_codes.extend(read.iter().map(|&b| base_code(b)));
        buf.ref_codes.clear();
        buf.ref_codes.extend(rf.iter().map(|&b| base_code(b)));

        // Row windows, laid out exactly as the scalar band: row 0 spans
        // [lo0, hi0] with H_START bytes (possibly empty), rows 1..=m span
        // [max(c−half,1), min(c+half,n)] with c = clamp(offset+i,1,n) — never
        // empty. Also the anti-diagonal extents a_lo = lo+i / a_hi = hi+i
        // (strictly increasing: windows shift right by 0/1 per row) and the
        // flat-scatter base = off − lo so lane i's dir byte lands at
        // base_i + j.
        buf.rows.clear();
        buf.rows_flat.clear();
        let lo0 = (offset - half).max(0);
        let hi0 = (offset + half).min(n);
        let len0 = if lo0 <= hi0 { (hi0 - lo0 + 1) as usize } else { 0 };
        buf.rows_flat.resize(len0, H_START);
        buf.rows.push(RowTrace {
            lo: lo0,
            off: 0,
            len: len0,
        });

        let m_us = m as usize;
        buf.avx_band_lo.clear();
        buf.avx_band_lo.resize(m_us + 2, 0);
        buf.avx_band_hi.clear();
        buf.avx_band_hi.resize(m_us + 2, 0);
        buf.avx_band_alo.clear();
        buf.avx_band_alo.resize(m_us + 2, 0);
        buf.avx_band_ahi.clear();
        buf.avx_band_ahi.resize(m_us + 2, 0);
        buf.avx_band_base.clear();
        buf.avx_band_base.resize(m_us + 2, 0);
        let mut off_acc = len0;
        for i in 1..=m_us {
            let c = (offset + i as i64).clamp(1, n);
            let lo = (c - half).max(1);
            let hi = (c + half).min(n);
            debug_assert!(lo <= hi, "rows 1..=m are never empty");
            let len = (hi - lo + 1) as usize;
            buf.avx_band_lo[i] = lo;
            buf.avx_band_hi[i] = hi;
            buf.avx_band_alo[i] = lo + i as i64;
            buf.avx_band_ahi[i] = hi + i as i64;
            buf.avx_band_base[i] = off_acc as i64 - lo;
            buf.rows.push(RowTrace {
                lo,
                off: off_acc,
                len,
            });
            off_acc += len;
        }
        buf.rows_flat.resize(off_acc, 0);

        let dmin = buf.avx_band_alo[1]; // = 1 + lo_1 ≥ 2: first cell (1, lo_1)
        let dmax = buf.avx_band_ahi[m_us]; // = m + hi_m

        // Scratch. Unlike the full kernel every buffer is re-filled each
        // call: never-written positions must read as the scalar margins
        // (H = 0, E/F = NEG), and a reused buffer carries the previous
        // call's tail.
        let cap = m_us + 33; // vector stores run one block past bhi
        for v in buf.avx_h.iter_mut() {
            v.resize(cap, 0);
            v.fill(0);
        }
        for v in buf.avx_e.iter_mut() {
            v.resize(cap, NEG);
            v.fill(NEG);
        }
        for v in buf.avx_f.iter_mut() {
            v.resize(cap, NEG);
            v.fill(NEG);
        }
        buf.avx_row_max.clear();
        buf.avx_row_max.resize(m_us + 1, NEG);
        buf.avx_row_arg.clear();
        buf.avx_row_arg.resize(m_us + 1, 1);
        buf.avx_row_hl.clear();
        buf.avx_row_hl.resize(m_us + 1, NEG);
        buf.avx_row_hr.clear();
        buf.avx_row_hr.resize(m_us + 1, NEG);
        buf.avx_read.clear();
        buf.avx_read.extend(buf.read_codes.iter().map(|&c| c as i16));
        buf.avx_read.resize(m_us + 16, 0);
        // Ref codes reversed (lane i of diagonal d reads ref column d−i):
        // refrev[q] = ref_codes[n−1−q] turns the gather into a load.
        buf.avx_refrev.clear();
        buf.avx_refrev.extend(
            (0..n as usize).map(|q| buf.ref_codes[(n as usize) - 1 - q] as i16),
        );
        buf.avx_refrev.resize(n as usize + 16, 0);

        let hp = [
            buf.avx_h[0].as_mut_ptr(),
            buf.avx_h[1].as_mut_ptr(),
            buf.avx_h[2].as_mut_ptr(),
        ];
        let ep = [buf.avx_e[0].as_mut_ptr(), buf.avx_e[1].as_mut_ptr()];
        let fp = [buf.avx_f[0].as_mut_ptr(), buf.avx_f[1].as_mut_ptr()];
        // SAFETY: AVX2 was checked by `available()` at dispatch; all
        // pointers stay inside the resized scratch (cap ≥ m + 33, code
        // arrays padded 16, base_i + j lands inside rows_flat because j ∈
        // [lo_i, hi_i] is in-band), and the buffers never alias.
        unsafe {
            dp_diagonals(
                m,
                n,
                dmin,
                dmax,
                buf.avx_read.as_ptr(),
                buf.avx_refrev.as_ptr(),
                hp,
                ep,
                fp,
                buf.rows_flat.as_mut_ptr(),
                buf.avx_band_lo.as_ptr(),
                buf.avx_band_hi.as_ptr(),
                buf.avx_band_alo.as_ptr(),
                buf.avx_band_ahi.as_ptr(),
                buf.avx_band_base.as_ptr(),
                buf.avx_row_max.as_mut_ptr(),
                buf.avx_row_arg.as_mut_ptr(),
                buf.avx_row_hl.as_mut_ptr(),
                buf.avx_row_hr.as_mut_ptr(),
                open_i,
                ext_i,
                s_diag,
                s_mm,
                s_ga,
                s_tc,
            );
        }

        // Row-major first strict maximum + first edge tie: rows ascending,
        // strict `>`, left edge before right — the scalar loop's exact
        // selection order (equal-value cells never displace an earlier
        // winner; `left_strict` = lo > 0 always holds since lo ≥ 1).
        let mut best_score = 0i16;
        let mut best_i = 0i64;
        let mut best_j = 0i64;
        let mut touch = TouchRec::default();
        for i in 1..=m_us {
            let rm = buf.avx_row_max[i];
            if !touch.touched {
                if buf.avx_band_lo[i] > 0 && buf.avx_row_hl[i] == rm {
                    touch = TouchRec {
                        touched: true,
                        row: i as i64,
                        left: true,
                        clamped: buf.avx_band_lo[i] <= 1,
                    };
                } else if buf.avx_band_hi[i] < n && buf.avx_row_hr[i] == rm {
                    touch = TouchRec {
                        touched: true,
                        row: i as i64,
                        left: false,
                        clamped: buf.avx_band_hi[i] >= n,
                    };
                }
            }
            if rm > best_score {
                best_score = rm;
                best_i = i as i64;
                best_j = buf.avx_row_arg[i] as i64;
            }
        }
        let ext = finish_extension(
            m,
            core,
            &buf.rows,
            &buf.rows_flat,
            &buf.read_codes,
            &buf.ref_codes,
            best_i,
            best_j,
        );
        (ext, touch)
    }

    /// The anti-diagonal sweep over the band. `h`/`e`/`f` are per-diagonal
    /// buffers indexed by row i with the same rotation as the full kernel
    /// (written at d, read at d+1/d+2, rewritten at d+3 / d+2). Per diagonal
    /// the valid rows are the contiguous range [blo, bhi] maintained by two
    /// monotone row cursors; blocks run exactly over it, so every lane the
    /// tail loop scatters is in-band. After the body, the margin schedule
    /// writes (H=0, E=F=NEG) at the three positions later diagonals read
    /// out-of-window: fresh right edge (row me whose hi stepped: its E read
    /// of (me−1, hi_me)), left edge (row mf: its F read of (mf, lo_mf−1)),
    /// and unshifted left edge (row md with lo_md == lo_{md−1}: its diagonal
    /// read of (md−1, lo_md−1)). Rows whose margin diagonal falls before the
    /// sweep never fire — those positions still hold the call-start fill,
    /// which is the margin value.
    #[target_feature(enable = "avx2")]
    #[allow(clippy::too_many_arguments)]
    unsafe fn dp_diagonals(
        m: i64,
        n: i64,
        dmin: i64,
        dmax: i64,
        read_i16: *const i16,
        refrev: *const i16,
        h: [*mut i16; 3],
        e: [*mut i16; 2],
        f: [*mut i16; 2],
        rows_flat: *mut u8,
        lo_arr: *const i64,
        hi_arr: *const i64,
        a_lo: *const i64,
        a_hi: *const i64,
        base_flat: *const i64,
        row_max: *mut i16,
        row_arg: *mut i32,
        row_hl: *mut i16,
        row_hr: *mut i16,
        open_ext: i16,
        ext: i16,
        s_diag: i16,
        s_mm: i16,
        s_ga: i16,
        s_tc: i16,
    ) {
        let neg_v = _mm256_set1_epi16(NEG);
        let zero_v = _mm256_set1_epi16(0);
        let ones_v = _mm256_set1_epi16(-1);
        let one_v = _mm256_set1_epi16(1);
        let two_v = _mm256_set1_epi16(2);
        let three_v = _mm256_set1_epi16(3);
        let four_v = _mm256_set1_epi16(4);
        let open_v = _mm256_set1_epi16(open_ext);
        let ext_v = _mm256_set1_epi16(ext);
        let diag_v = _mm256_set1_epi16(s_diag);
        let mm_v = _mm256_set1_epi16(s_mm);
        let ga_v = _mm256_set1_epi16(s_ga);
        let tc_v = _mm256_set1_epi16(s_tc);
        let code_g = _mm256_set1_epi16(2);
        let code_a = _mm256_set1_epi16(0);
        let code_c = _mm256_set1_epi16(1);
        let code_t = _mm256_set1_epi16(3);

        // `Acc::add` with a constant non-NEG operand (as in `avx2_full`).
        macro_rules! addc {
            ($x:expr, $c:expr) => {{
                let isn = _mm256_cmpeq_epi16($x, neg_v);
                let r = _mm256_max_epi16(_mm256_add_epi16($x, $c), neg_v);
                _mm256_blendv_epi8(r, neg_v, isn)
            }};
        }

        let mut h2 = 0usize; // roles: h2 = d−2, h1 = d−1, h0 = d
        let mut h1 = 1usize;
        let mut h0 = 2usize;
        let mut e1 = 0usize; // e1 = d−1, e0 = d
        let mut e0 = 1usize;
        let mut f1 = 0usize;
        let mut f0 = 1usize;
        // Monotone row cursors: blo_c = min{i : a_hi[i] ≥ d} and
        // bhi_c = max{i : a_lo[i] ≤ d} (both move only forward in d); and
        // the three margin-schedule cursors (events keyed on a_lo/a_hi,
        // strictly increasing, so each fires on exactly one diagonal).
        let mut blo_c = 1i64;
        let mut bhi_c = 1i64;
        let mut me_r = 1i64; // E margin: row with a_hi == d+1, hi step
        let mut mf_r = 1i64; // F margin: row with a_lo == d+1 (every row)
        let mut md_r = 1i64; // diag margin: row with a_lo == d+2, lo flat

        for d in dmin..=dmax {
            let h2p = h[h2];
            let h1p = h[h1];
            let h0p = h[h0];
            let e1p = e[e1];
            let e0p = e[e0];
            let f1p = f[f1];
            let f0p = f[f0];
            while blo_c <= m && *a_hi.add(blo_c as usize) < d {
                blo_c += 1;
            }
            while bhi_c < m && *a_lo.add((bhi_c + 1) as usize) <= d {
                bhi_c += 1;
            }
            let blo = blo_c;
            let bhi = bhi_c;
            let mut v = blo;
            while v <= bhi {
                // Lane i covers cell (i, d − i); all lanes in [blo, bhi]
                // are in-band, so no per-lane validity check is needed.
                let rc = _mm256_loadu_si256(read_i16.add((v - 1) as usize) as *const __m256i);
                let rfc =
                    _mm256_loadu_si256(refrev.add((v + (n - d)) as usize) as *const __m256i);
                let eq = _mm256_cmpeq_epi16(rc, rfc);
                let is_match = _mm256_and_si256(eq, _mm256_cmpgt_epi16(four_v, rc));
                let ga = _mm256_and_si256(
                    _mm256_cmpeq_epi16(rc, code_g),
                    _mm256_cmpeq_epi16(rfc, code_a),
                );
                let tc = _mm256_and_si256(
                    _mm256_cmpeq_epi16(rc, code_c),
                    _mm256_cmpeq_epi16(rfc, code_t),
                );
                let mut s = _mm256_blendv_epi8(mm_v, ga_v, ga);
                s = _mm256_blendv_epi8(s, tc_v, tc);
                s = _mm256_blendv_epi8(s, diag_v, is_match);

                // E (vertical): from H[i−1][j] / E[i−1][j], one diagonal up.
                let h_up = _mm256_loadu_si256(h1p.add((v - 1) as usize) as *const __m256i);
                let e_up = _mm256_loadu_si256(e1p.add((v - 1) as usize) as *const __m256i);
                let e_open = addc!(h_up, open_v);
                let e_ext = addc!(e_up, ext_v);
                let egt = _mm256_cmpgt_epi16(e_ext, e_open);
                let e_val = _mm256_blendv_epi8(e_open, e_ext, egt);
                let edir = _mm256_and_si256(egt, one_v);

                // F (horizontal): from H[i][j−1] / F[i][j−1], same diagonal.
                let h_left = _mm256_loadu_si256(h1p.add(v as usize) as *const __m256i);
                let f_left = _mm256_loadu_si256(f1p.add(v as usize) as *const __m256i);
                let f_open = addc!(h_left, open_v);
                let f_ext = addc!(f_left, ext_v);
                let fgt = _mm256_cmpgt_epi16(f_ext, f_open);
                let f_val = _mm256_blendv_epi8(f_open, f_ext, fgt);
                let fdir = _mm256_and_si256(fgt, one_v);

                // Diagonal: H[i−1][j−1] + s.
                let h_diag = _mm256_loadu_si256(h2p.add((v - 1) as usize) as *const __m256i);
                let diag = addc!(h_diag, s);

                // H: zero floor, then F → E → diag (later source wins ties),
                // then the ==0 rewrite to H_START. Lane-for-lane the scalar
                // cascade (identical to `avx2_full`).
                let mut hv = zero_v;
                let mut hd = three_v;
                let fmask = _mm256_andnot_si256(
                    _mm256_or_si256(_mm256_cmpeq_epi16(f_val, neg_v), _mm256_cmpgt_epi16(hv, f_val)),
                    ones_v,
                );
                hv = _mm256_blendv_epi8(hv, f_val, fmask);
                hd = _mm256_blendv_epi8(hd, two_v, fmask);
                let emask = _mm256_andnot_si256(
                    _mm256_or_si256(_mm256_cmpeq_epi16(e_val, neg_v), _mm256_cmpgt_epi16(hv, e_val)),
                    ones_v,
                );
                hv = _mm256_blendv_epi8(hv, e_val, emask);
                hd = _mm256_blendv_epi8(hd, one_v, emask);
                let dmask = _mm256_andnot_si256(
                    _mm256_or_si256(_mm256_cmpeq_epi16(diag, neg_v), _mm256_cmpgt_epi16(hv, diag)),
                    ones_v,
                );
                hv = _mm256_blendv_epi8(hv, diag, dmask);
                hd = _mm256_blendv_epi8(hd, zero_v, dmask);
                hd = _mm256_blendv_epi8(hd, three_v, _mm256_cmpeq_epi16(hv, zero_v));

                let dir16 = _mm256_or_si256(
                    hd,
                    _mm256_or_si256(_mm256_slli_epi16(edir, 2), _mm256_slli_epi16(fdir, 4)),
                );

                _mm256_storeu_si256(h0p.add(v as usize) as *mut __m256i, hv);
                _mm256_storeu_si256(e0p.add(v as usize) as *mut __m256i, e_val);
                _mm256_storeu_si256(f0p.add(v as usize) as *mut __m256i, f_val);

                // Per-lane tail: dir byte scatter into the per-row flat
                // trace window (address base_i + j), row max with
                // first-strict-improving column, and the row's first/last
                // in-band H for the touch record.
                let mut hb = [0i16; 16];
                let mut db = [0i16; 16];
                _mm256_storeu_si256(hb.as_mut_ptr() as *mut __m256i, hv);
                _mm256_storeu_si256(db.as_mut_ptr() as *mut __m256i, dir16);
                for k in 0..16usize {
                    let i = v + k as i64;
                    if i > bhi {
                        break;
                    }
                    let j = d - i;
                    let iu = i as usize;
                    *rows_flat.add((*base_flat.add(iu) + j) as usize) = db[k] as u8;
                    let hv16 = hb[k];
                    let rm = row_max.add(iu);
                    if hv16 > *rm {
                        *rm = hv16;
                        *row_arg.add(iu) = j as i32;
                    }
                    if j == *lo_arr.add(iu) {
                        *row_hl.add(iu) = hv16;
                    }
                    if j == *hi_arr.add(iu) {
                        *row_hr.add(iu) = hv16;
                    }
                }
                v += 16;
            }

            // Margin writes for this diagonal, after the body (garbage
            // lanes past bhi may have stored here; a margin position is
            // never a valid cell of its own diagonal). All three values
            // are written uniformly — only one field is read per case.
            while me_r <= m && *a_hi.add(me_r as usize) < d + 1 {
                me_r += 1;
            }
            if me_r <= m
                && *a_hi.add(me_r as usize) == d + 1
                && me_r >= 2
                && *hi_arr.add(me_r as usize) == *hi_arr.add((me_r - 1) as usize) + 1
            {
                let p = (me_r - 1) as usize;
                *h0p.add(p) = 0;
                *e0p.add(p) = NEG;
                *f0p.add(p) = NEG;
            }
            while mf_r <= m && *a_lo.add(mf_r as usize) < d + 1 {
                mf_r += 1;
            }
            if mf_r <= m && *a_lo.add(mf_r as usize) == d + 1 {
                let p = mf_r as usize;
                *h0p.add(p) = 0;
                *e0p.add(p) = NEG;
                *f0p.add(p) = NEG;
            }
            while md_r <= m && *a_lo.add(md_r as usize) < d + 2 {
                md_r += 1;
            }
            if md_r <= m
                && *a_lo.add(md_r as usize) == d + 2
                && md_r >= 2
                && *lo_arr.add(md_r as usize) == *lo_arr.add((md_r - 1) as usize)
            {
                let p = (md_r - 1) as usize;
                *h0p.add(p) = 0;
                *e0p.add(p) = NEG;
                *f0p.add(p) = NEG;
            }

            // Rotate roles.
            let t = h2;
            h2 = h1;
            h1 = h0;
            h0 = t;
            std::mem::swap(&mut e1, &mut e0);
            std::mem::swap(&mut f1, &mut f0);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The i16 sentinel contract: exact-equality NEG with a max() clamp —
    /// deep E/F chains on long windows all collapse to NEG (load-bearing for
    /// tie behavior; see the impl's NOTE).
    #[test]
    fn i16_sentinel_clamp_semantics() {
        type A = i16;
        assert_eq!(<A as Acc>::NEG.add(<A as Acc>::from_i32(6)), <A as Acc>::NEG);
        assert_eq!(<A as Acc>::NEG.add(<A as Acc>::NEG), <A as Acc>::NEG);
        assert_eq!(<A as Acc>::from_i32(0).add(<A as Acc>::NEG), <A as Acc>::NEG);
        assert_eq!(<A as Acc>::from_i32(300).add(<A as Acc>::from_i32(2)), 302);
        // Deep chains clamp to the sentinel, not saturate below it.
        assert_eq!(<A as Acc>::from_i32(-31000).add(<A as Acc>::from_i32(-6)), <A as Acc>::NEG);
    }

    // ---- AVX2 full-matrix kernel vs scalar kernel (bit-identity) ----

    #[cfg(target_arch = "x86_64")]
    fn assert_ext_same(a: &Extension, b: &Extension, ctx: &str) {
        assert_eq!(a.score, b.score, "score {ctx}");
        assert_eq!(a.read_start, b.read_start, "read_start {ctx}");
        assert_eq!(a.read_end, b.read_end, "read_end {ctx}");
        assert_eq!(a.ref_start, b.ref_start, "ref_start {ctx}");
        assert_eq!(a.cigar, b.cigar, "cigar {ctx}");
    }

    #[cfg(target_arch = "x86_64")]
    fn assert_touch_same(a: &TouchRec, b: &TouchRec, ctx: &str) {
        assert_eq!(a.touched, b.touched, "touched {ctx}");
        assert_eq!(a.row, b.row, "row {ctx}");
        assert_eq!(a.left, b.left, "left {ctx}");
        assert_eq!(a.clamped, b.clamped, "clamped {ctx}");
    }

    /// Deterministic xorshift over ACGTN.
    #[cfg(target_arch = "x86_64")]
    struct Rng(u64);
    #[cfg(target_arch = "x86_64")]
    impl Rng {
        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
        fn base(&mut self) -> u8 {
            b"ACGTN"[self.below(5) as usize]
        }
    }

    #[cfg(target_arch = "x86_64")]
    fn rand_seq(rng: &mut Rng, len: usize) -> Vec<u8> {
        (0..len).map(|_| rng.base()).collect()
    }

    /// Window with the read embedded at a random offset (with ~1/8
    /// substitutions and 1/16 single-base drops) inside random flanks —
    /// real alignment structure, not just noise.
    #[cfg(target_arch = "x86_64")]
    fn rand_window(rng: &mut Rng, read: &[u8], n: usize) -> Vec<u8> {
        let mut w = rand_seq(rng, n);
        if n > read.len() {
            let off = rng.below((n - read.len()) as u64 + 1) as usize;
            let mut k = 0usize;
            while k < read.len() && off + k < n {
                let b = read[k];
                w[off + k] = if rng.below(8) == 0 { rng.base() } else { b };
                if rng.below(16) == 0 {
                    k += 2; // dropped base ⇒ indel
                } else {
                    k += 1;
                }
            }
        }
        w
    }

    #[cfg(target_arch = "x86_64")]
    fn cmp_case(read: &[u8], win: &[u8], dir: EditDir, ctx: &str) {
        let params = ExtendParams {
            edit_dir: dir,
            ..ExtendParams::default()
        };
        let core = CoreParams::new(&params);
        let mut buf = ExtendBuffer::new();
        let (sim_ext, sim_touch) = avx2_full::run(read, win, &core, &mut buf.packed);
        let (sca_ext, sca_touch) = run_band(read, win, &core, &mut buf.packed, None);
        assert_ext_same(&sim_ext, &sca_ext, ctx);
        assert_touch_same(&sim_touch, &sca_touch, ctx);
    }

    /// Randomized ACGTN read/window pairs across lengths (including edge
    /// shapes m = 1, n < 16, n not a multiple of 16), all three editing
    /// directions: the AVX2 kernel must return the exact scalar Extension.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx2_full_matrix_matches_scalar() {
        if !avx2_full::available() {
            eprintln!("skipping: AVX2 not available");
            return;
        }
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        let lens_m = [1usize, 2, 3, 7, 8, 15, 16, 17, 31, 32, 33, 63, 64, 100, 150, 151, 200];
        let lens_n = [1, 2, 3, 7, 8, 15, 16, 17, 31, 48, 63, 64, 100, 127, 255, 257, 1000, 3000];
        let dirs = [EditDir::Off, EditDir::Ag, EditDir::Tc];
        let mut cases = 0usize;
        for &m in lens_m.iter() {
            for &n in lens_n.iter() {
                if n < m {
                    continue;
                }
                for rep in 0..3 {
                    let read = rand_seq(&mut rng, m);
                    let win = if rep == 2 {
                        rand_seq(&mut rng, n)
                    } else {
                        rand_window(&mut rng, &read, n)
                    };
                    let dir = dirs[rng.below(3) as usize];
                    cmp_case(&read, &win, dir, &format!("m={m} n={n} rep={rep} dir={dir:?}"));
                    cases += 1;
                }
            }
        }
        // Long mate-anchor-scale window (n = 20,000, the documented max).
        for rep in 0..3 {
            let read = rand_seq(&mut rng, 150);
            let win = if rep == 2 {
                rand_seq(&mut rng, 20000)
            } else {
                rand_window(&mut rng, &read, 20000)
            };
            cmp_case(&read, &win, EditDir::Ag, &format!("long rep={rep}"));
            cases += 1;
        }
        // Degenerate: empty read / empty window (both kernels shortcut).
        cmp_case(b"", b"ACGT", EditDir::Off, "m=0");
        cmp_case(b"ACGT", b"", EditDir::Off, "n=0");
        assert!(cases > 300, "expected a broad sweep, got {cases}");
    }

    /// Exhaustive small shapes around the 16-lane block boundaries.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx2_full_matrix_block_boundaries() {
        if !avx2_full::available() {
            eprintln!("skipping: AVX2 not available");
            return;
        }
        let mut rng = Rng(0x0F1E_2D3C_4B5A_6978);
        for m in 14usize..=20 {
            for n in m..=20 { // production contract: read fits the window (m <= n)
                for rep in 0..4 {
                    let read = rand_seq(&mut rng, m);
                    let win = if rep & 1 == 0 {
                        rand_window(&mut rng, &read, n)
                    } else {
                        rand_seq(&mut rng, n)
                    };
                    let dir = [EditDir::Off, EditDir::Ag, EditDir::Tc][(rep % 3) as usize];
                    cmp_case(&read, &win, dir, &format!("m={m} n={n} rep={rep}"));
                }
            }
        }
        for &n in [30usize, 31, 32, 33, 47, 63, 64, 65].iter() {
            let read = rand_seq(&mut rng, 17);
            let win = rand_window(&mut rng, &read, n);
            cmp_case(&read, &win, EditDir::Tc, &format!("tail n={n}"));
        }
    }

    /// The dispatch path itself (extend → run_dispatch → AVX2 on this
    /// machine) must equal the scalar kernel, including below-threshold
    /// windows that stay scalar.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx2_dispatch_path_matches_scalar() {
        if !avx2_full::available() {
            eprintln!("skipping: AVX2 not available");
            return;
        }
        let mut rng = Rng(0xDEAD_BEEF_CAFE_1234);
        for (m, n) in [(6usize, 6usize), (8, 8), (9, 40), (50, 300), (120, 500)] {
            for rep in 0..4 {
                let read = rand_seq(&mut rng, m);
                let win = if rep == 3 {
                    rand_seq(&mut rng, n)
                } else {
                    rand_window(&mut rng, &read, n)
                };
                let dir = [EditDir::Off, EditDir::Ag, EditDir::Tc][(rep % 3) as usize];
                let params = ExtendParams {
                    edit_dir: dir,
                    ..ExtendParams::default()
                };
                let mut buf = ExtendBuffer::new();
                let e = extend(&read, &win, &params, &mut buf);
                let core = CoreParams::new(&params);
                let (sca, _) = run_band(&read, &win, &core, &mut buf.packed, None);
                assert_ext_same(&e, &sca, &format!("dispatch m={m} n={n} rep={rep}"));
            }
        }
    }


    // ---- AVX2 banded kernel vs scalar kernel (bit-identity) ----

    #[cfg(target_arch = "x86_64")]
    fn cmp_band_case(read: &[u8], win: &[u8], dir: EditDir, offset: i64, half: i64, ctx: &str) {
        let params = ExtendParams {
            edit_dir: dir,
            ..ExtendParams::default()
        };
        let core = CoreParams::new(&params);
        let mut buf = ExtendBuffer::new();
        let (sim_ext, sim_touch) =
            avx2_band::run(read, win, &core, &mut buf.packed, offset, 1, 1, half);
        let (sca_ext, sca_touch) =
            run_band(read, win, &core, &mut buf.packed, Some((offset, 1, 1, half)));
        assert_ext_same(&sim_ext, &sca_ext, ctx);
        assert_touch_same(&sim_touch, &sca_touch, ctx);
    }

    /// Randomized banded equivalence: offsets negative/zero/mid/right-edge/
    /// beyond-right/huge (including empty row-0 windows on both sides),
    /// half-widths 0/1/2/30/huge (huge ⇒ the band degenerates to the full
    /// window through the banded code path), all three editing directions,
    /// m ≤ n shapes plus a couple m > n, structured and noise windows.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx2_band_matches_scalar() {
        if !avx2_full::available() {
            eprintln!("skipping: AVX2 not available");
            return;
        }
        let mut rng = Rng(0xA5A5_5A5A_C3C3_3C3C);
        let shapes = [
            (8usize, 9usize),
            (8, 40),
            (15, 16),
            (16, 16),
            (17, 33),
            (20, 200),
            (31, 31),
            (32, 64),
            (60, 61),
            (60, 300),
            (100, 300),
            (150, 151),
            (150, 1000),
            (40, 20),  // m > n: scalar supports it; kernel must too
            (17, 9),
        ];
        let dirs = [EditDir::Off, EditDir::Ag, EditDir::Tc];
        let mut cases = 0usize;
        for &(m, n) in shapes.iter() {
            let mn = (m as i64).max(n as i64);
            let halfs: [i64; 5] = if n <= 64 {
                [0, 1, 2, 30, mn + 5] // huge covers the whole window
            } else {
                [0, 1, 2, 30, 300]
            };
            for rep in 0..2 {
                let read = rand_seq(&mut rng, m);
                let win = if rep == 1 {
                    rand_seq(&mut rng, n)
                } else {
                    rand_window(&mut rng, &read, n)
                };
                for &half in halfs.iter() {
                    let hh = half.min(mn);
                    for &offset in [
                        -hh - 5,             // row 0 empty (left)
                        -1i64,
                        0,
                        5,
                        (n as i64 - m as i64), // band along the read diagonal
                        n as i64,           // right window edge
                        n as i64 + hh + 5,  // row 0 empty (right)
                        4 * n as i64,       // huge
                    ]
                    .iter()
                    {
                        let dir = dirs[rng.below(3) as usize];
                        cmp_band_case(
                            &read,
                            &win,
                            dir,
                            offset,
                            half,
                            &format!(
                                "m={m} n={n} rep={rep} half={half} offset={offset} dir={dir:?}"
                            ),
                        );
                        cases += 1;
                    }
                }
            }
        }
        assert!(cases > 500, "expected a broad sweep, got {cases}");
        // Long mate-anchor-scale windows (n = 20,000, the documented max),
        // band offsets spanning the window.
        for rep in 0..3 {
            let read = rand_seq(&mut rng, 150);
            let win = if rep == 2 {
                rand_seq(&mut rng, 20000)
            } else {
                rand_window(&mut rng, &read, 20000)
            };
            for &(offset, half) in [
                (-40i64, 30i64),
                (0, 30),
                (5000, 30),
                (19850, 30),
                (-20000, 20000),
                (20000, 20000),
            ]
            .iter()
            {
                cmp_band_case(
                    &read,
                    &win,
                    EditDir::Ag,
                    offset,
                    half,
                    &format!("long rep={rep} off={offset} half={half}"),
                );
            }
        }
    }

    /// Exhaustive small shapes around the 16-lane block boundaries: window
    /// widths crossing 16/32/48 (half 0..=48), reads 8..=24, n just past m.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx2_band_block_boundaries() {
        if !avx2_full::available() {
            eprintln!("skipping: AVX2 not available");
            return;
        }
        let mut rng = Rng(0x1357_9BDF_2468_ACE0);
        let halfs = [0i64, 1, 2, 7, 8, 15, 16, 17, 30, 31, 32, 33, 48];
        for m in 8usize..=24 {
            for &dn in [0usize, 1, 15, 16, 17, 31, 32, 33].iter() {
                let n = m + dn;
                for rep in 0..2 {
                    let read = rand_seq(&mut rng, m);
                    let win = if rep == 1 {
                        rand_seq(&mut rng, n)
                    } else {
                        rand_window(&mut rng, &read, n)
                    };
                    for &half in halfs.iter() {
                        // offsets around the empty-row-0 boundaries and the
                        // clamp regions on both sides.
                        let offset = (rng.below(120) as i64) - 60 + (n as i64 - m as i64);
                        let dir = [EditDir::Off, EditDir::Ag, EditDir::Tc][(rep % 3) as usize];
                        cmp_band_case(
                            &read,
                            &win,
                            dir,
                            offset,
                            half,
                            &format!("blk m={m} n={n} rep={rep} half={half} off={offset}"),
                        );
                    }
                }
            }
        }
    }

    /// The banded dispatch path itself (run_dispatch with Some(..) → avx2_band
    /// on this machine) must equal the scalar kernel, including the windows
    /// that stay scalar (num ≠ den).
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx2_band_dispatch_path_matches_scalar() {
        if !avx2_full::available() {
            eprintln!("skipping: AVX2 not available");
            return;
        }
        let mut rng = Rng(0x0BAD_C0DE_1234_5678);
        for (m, n) in [(8usize, 9usize), (10, 300), (60, 400), (120, 500)] {
            for rep in 0..3 {
                let read = rand_seq(&mut rng, m);
                let win = if rep == 2 {
                    rand_seq(&mut rng, n)
                } else {
                    rand_window(&mut rng, &read, n)
                };
                let dir = [EditDir::Off, EditDir::Ag, EditDir::Tc][(rep % 3) as usize];
                let params = ExtendParams {
                    edit_dir: dir,
                    ..ExtendParams::default()
                };
                let offset = (rng.below(n as u64) as i64) - 30;
                let half = [0i64, 1, 2, 30][(rep % 4) as usize];
                let mut buf = ExtendBuffer::new();
                let (dsp_ext, dsp_touch) =
                    run_dispatch(&read, &win, &params, &mut buf, Some((offset, 1, 1, half)));
                let core = CoreParams::new(&params);
                let (sca_ext, sca_touch) =
                    run_band(&read, &win, &core, &mut buf.packed, Some((offset, 1, 1, half)));
                assert_ext_same(
                    &dsp_ext,
                    &sca_ext,
                    &format!("band dispatch m={m} n={n} rep={rep} off={offset} half={half}"),
                );
                assert_touch_same(
                    &dsp_touch,
                    &sca_touch,
                    &format!("band dispatch m={m} n={n} rep={rep} off={offset} half={half}"),
                );
            }
        }
    }

    /// ESP_NO_SIMD=1 masks the SIMD gate (escape hatch / debug aid): checked
    /// in a fresh child process so the parent's OnceLock state (and the
    /// other tests' dispatching) cannot interfere.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn esp_no_simd_escape_hatch() {
        if std::env::var_os("ESP_CHILD_NO_SIMD").is_some() {
            // Child: ESP_NO_SIMD=1 was set before anything initialized the
            // gate; dispatch must have fallen back to the scalar kernel.
            assert!(!avx2_full::available());
            let read = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
            let win = b"TTTACGTACGTACGTACGTAAGGTACGTACGTACGTACGTAA";
            let mut buf = ExtendBuffer::new();
            let e = extend(read, win, &ExtendParams::default(), &mut buf);
            assert!(e.score > 0, "scalar fallback must still align: {e:?}");
            // The banded dispatch (extend_hint attempt 1) shares the same
            // gate and must also have run scalar without incident.
            let h = extend_hint(
                read,
                win,
                &ExtendParams::default(),
                &mut buf,
                DiagHint { offset: 3, num: 1, den: 1 },
            );
            assert!(h.score > 0, "scalar banded fallback must still align: {h:?}");
            return;
        }
        let exe = std::env::current_exe().unwrap();
        let st = std::process::Command::new(exe)
            .arg("extend::tests::esp_no_simd_escape_hatch")
            .arg("--exact")
            .arg("--nocapture")
            .env("ESP_NO_SIMD", "1")
            .env("ESP_CHILD_NO_SIMD", "1")
            .status()
            .unwrap();
        assert!(st.success());
    }

}
