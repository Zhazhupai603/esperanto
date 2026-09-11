//! Optional DP profiling counters. Always-on atomic increments (negligible
//! cost); one JSON line on stderr at run end when ESP_PROF=1 is set.
//!
//! The shadow-ladder block below (ESP_SHADOW=1) is a parallel stats
//! sink of the same shape; the work order pins it beside the counters it
//! reports next to.
//! allow: SIZE_OK — one cohesive counter-sink module; the shadow block is
//! deliberately parallel to (and disjoint from) the table above.
//! cost); one JSON line on stderr at run end when ESP_PROF=1 is set.

use std::sync::atomic::{AtomicU64, Ordering};

/// Counter indices (stable order; report prints names).
pub const READS: usize = 0;
pub const DP_UNSPLICED: usize = 1;
pub const DP_SPLICED: usize = 2;
pub const DP_TAIL: usize = 3;
pub const DP_SPAN: usize = 4;
pub const DP_HYPER: usize = 5;
pub const DP_MATE: usize = 6;
pub const DP_FULL_MATRIX: usize = 7;
pub const DP_CELLS: usize = 8;
pub const FAST_LANE_HITS: usize = 9;
pub const MYERS_ONLY_HITS: usize = 10;
pub const MATE_REALIGNED: usize = 11;
pub const DP_MID_BAND: usize = 12;
pub const DP_FULL_EDGE: usize = 13;
pub const DP_FULL_SCORE: usize = 14;
pub const ESC_LEFT_CLAMPED: usize = 15;
pub const ESC_LEFT_TRUE: usize = 16;
pub const ESC_RIGHT: usize = 17;
const NAMES: [&str; 18] = [
    "reads",
    "dp_calls_unspliced",
    "dp_calls_spliced",
    "dp_calls_tail",
    "dp_calls_span",
    "dp_calls_hyper",
    "dp_calls_mate",
    "dp_full_matrix",
    "dp_cells",
    "fast_lane_hits",
    "myers_only_hits",
    "mate_realigned",
    "dp_mid_band",
    "dp_full_edge",
    "dp_full_score",
    "esc_left_clamped",
    "esc_left_true",
    "esc_right",
];

static COUNTERS: [AtomicU64; 18] = [const { AtomicU64::new(0) }; 18];

pub fn bump(i: usize, n: u64) {
    COUNTERS[i].fetch_add(n, Ordering::Relaxed);
}

/// Bump a DP call-site counter and mark it as this thread's current site.
pub fn dp_call(i: usize) {
    bump(i, 1);
    set_site(i);
}

pub fn enabled() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("ESP_PROF").is_some())
}

/// One JSON object of the counter table (empty string when disabled).
pub fn report_line() -> String {
    if !enabled() {
        return String::new();
    }
    let mut s = String::from("[prof] {");
    for (i, name) in NAMES.iter().enumerate() {
        if i > 0 {
            s.push_str(", ");
        }
        s.push_str(&format!(
            "\"{}\": {}",
            name,
            COUNTERS[i].load(Ordering::Relaxed)
        ));
    }
    s.push('}');
    s
}

// ---- Escalation touch attribution (instrumentation only) ----

use std::cell::Cell;
use std::sync::Mutex;

thread_local! {
    static SITE: Cell<usize> = const { Cell::new(usize::MAX) };
}

/// Mark the current thread's DP call site (use the DP_* counter id).
pub fn set_site(s: usize) {
    SITE.with(|c| c.set(s));
}

/// Current thread's DP call site (`usize::MAX` when unset).
pub fn site() -> usize {
    SITE.with(|c| c.get())
}

/// One escalated DP, sampled for attribution. `row` is the first touched row
/// (0 when the escalation was score-gated, not edge-gated); `clamped` marks a
/// "touch" sitting on the window boundary (column 1 or n) rather than on a
/// true band edge; `start_diag` is the full matrix's `ref_start - read_start`.
pub struct TouchEvent {
    pub site: usize,
    pub row: i64,
    pub m: i64,
    pub n: i64,
    pub left: bool,
    pub clamped: bool,
    pub by_score: bool,
    pub offset: i64,
    pub start_diag: i64,
    pub band_score: i32,
    pub full_score: i32,
}

const TOUCH_CAP: usize = 10_000;
static TOUCH_EVENTS: Mutex<Vec<TouchEvent>> = Mutex::new(Vec::new());
static TOUCH_SEEN: AtomicU64 = AtomicU64::new(0);
static TOUCH_RNG: AtomicU64 = AtomicU64::new(0x9E37_79B9_7F4A_7C15);

/// Reservoir-sample an escalation event (cap 10k; deterministic xorshift).
pub fn record_touch(ev: TouchEvent) {
    if !enabled() {
        return;
    }
    let seen = TOUCH_SEEN.fetch_add(1, Ordering::Relaxed);
    let mut g = match TOUCH_EVENTS.lock() {
        Ok(g) => g,
        Err(p) => p.into_inner(),
    };
    if (g.len()) < TOUCH_CAP {
        g.push(ev);
        return;
    }
    let mut x = TOUCH_RNG.load(Ordering::Relaxed);
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    TOUCH_RNG.store(x, Ordering::Relaxed);
    let j = (x % (seen + 1)) as usize;
    if j < TOUCH_CAP {
        g[j] = ev;
    }
}

/// One `[touch] {...}` JSON line per sampled event (empty when disabled).
pub fn touch_lines() -> String {
    if !enabled() {
        return String::new();
    }
    let mut g = match TOUCH_EVENTS.lock() {
        Ok(g) => g,
        Err(p) => p.into_inner(),
    };
    let g = std::mem::take(&mut *g);
    let mut s = String::new();
    for ev in g.iter() {
        s.push_str(&format!(
            "\n[touch] {{\"site\": {}, \"row\": {}, \"m\": {}, \"n\": {}, \"left\": {}, \"clamped\": {}, \"by_score\": {}, \"offset\": {}, \"start_diag\": {}, \"band_score\": {}, \"full_score\": {}}}",
            ev.site,
            ev.row,
            ev.m,
            ev.n,
            ev.left,
            ev.clamped,
            ev.by_score,
            ev.offset,
            ev.start_diag,
            ev.band_score,
            ev.full_score
        ));
    }
    s
}

// ---- Shadow ladder stats (ESP_SHADOW=1; measurement only) ----
//
// Parallel stats block for the would-be band-64 middle tier of the
// escalation ladder: `extend_hint` re-runs escalated DPs at half=64
// behind ESP_SHADOW and records the class here. Nothing in this block
// feeds any production decision; these counters are separate from
// `COUNTERS` above and nothing above them changes.

/// band64 == full, field-for-field: a `band30 → band64 → full` ladder
/// could accept the band-64 answer here and skip the full matrix.
pub const SHADOW_CLASS_A: usize = 0;
/// band64 == band30 != full: stability is NOT proof — a ladder that
/// accepted "band-64 agrees with band-30" would return the wrong answer.
pub const SHADOW_CLASS_B: usize = 1;
/// All three answers differ.
pub const SHADOW_CLASS_C: usize = 2;
const SHADOW_CLASSES: usize = 3;

/// Trigger split mirrors DP_FULL_EDGE / DP_FULL_SCORE.
const SHADOW_TRIGGERS: usize = 2;
/// One row per DP_* site id (NAMES.len() of them) plus one overflow row
/// for an unset site thread-local.
const SHADOW_SITE_UNSET: usize = 18;
const SHADOW_SITE_ROWS: usize = 19;

/// Read-length bin edges (upper-exclusive; last bin open-ended).
const SHADOW_READ_EDGES: [i64; 9] = [50, 75, 100, 125, 150, 200, 250, 300, 400];
/// Window-length bin edges: the full-matrix cost a ladder would avoid is
/// ~ read_len × window_len, so the payoff input is the window length.
const SHADOW_WIN_EDGES: [i64; 9] = [100, 200, 300, 500, 1000, 2000, 5000, 10_000, 20_000];
const SHADOW_LEN_BINS: usize = 10;

static SHADOW_BY_SITE: [AtomicU64; SHADOW_SITE_ROWS * SHADOW_TRIGGERS * SHADOW_CLASSES] =
    [const { AtomicU64::new(0) }; SHADOW_SITE_ROWS * SHADOW_TRIGGERS * SHADOW_CLASSES];
static SHADOW_BY_READ: [AtomicU64; SHADOW_LEN_BINS * SHADOW_CLASSES] =
    [const { AtomicU64::new(0) }; SHADOW_LEN_BINS * SHADOW_CLASSES];
static SHADOW_BY_WIN: [AtomicU64; SHADOW_LEN_BINS * SHADOW_CLASSES] =
    [const { AtomicU64::new(0) }; SHADOW_LEN_BINS * SHADOW_CLASSES];
static SHADOW_TOTAL: AtomicU64 = AtomicU64::new(0);

/// §五: min(left, right) unused window margin of the full-matrix optimum,
/// 8-bp buckets (0-7 .. 120-127, 128+). Per site.
const SHADOW_SLACK_BINS: usize = 17;
static SHADOW_SLACK: [AtomicU64; SHADOW_SITE_ROWS * SHADOW_SLACK_BINS] =
    [const { AtomicU64::new(0) }; SHADOW_SITE_ROWS * SHADOW_SLACK_BINS];
static SHADOW_CELLS: AtomicU64 = AtomicU64::new(0);

/// ESP_SHADOW=1 gate — same shape as [`enabled`]; caches on first call.
pub fn shadow_enabled() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("ESP_SHADOW").is_some())
}

/// One escalated DP measured by the shadow ladder.
pub struct ShadowEvent {
    /// DP call-site id (`DP_*`; `usize::MAX` when the thread-local is unset).
    pub site: usize,
    /// Read length.
    pub read_len: i64,
    /// Reference window length.
    pub win_len: i64,
    /// Escalation triggered by a band-edge touch (vs the score gate).
    pub by_edge: bool,
    /// SHADOW_CLASS_{A,B,C}.
    pub class: usize,
    /// Window margin unused by the full-matrix optimum on the tighter side:
    /// min(ref_start, win_len - ref_end) in bp (§五 window-shrink ceiling).
    pub min_slack: i64,
    /// Full-matrix cells this escalation paid (m * n).
    pub full_cells: u64,
}

fn shadow_site_row(site: usize) -> usize {
    if site < NAMES.len() {
        site
    } else {
        SHADOW_SITE_UNSET
    }
}

fn shadow_bin(len: i64, edges: &[i64; 9]) -> usize {
    edges
        .iter()
        .position(|&e| len < e)
        .unwrap_or(SHADOW_LEN_BINS - 1)
}

/// Record one shadow-ladder event (callers gate on [`shadow_enabled`]).
pub fn shadow_record(ev: ShadowEvent) {
    let row = shadow_site_row(ev.site);
    let trig = usize::from(!ev.by_edge);
    SHADOW_BY_SITE[(row * SHADOW_TRIGGERS + trig) * SHADOW_CLASSES + ev.class]
        .fetch_add(1, Ordering::Relaxed);
    SHADOW_BY_READ[shadow_bin(ev.read_len, &SHADOW_READ_EDGES) * SHADOW_CLASSES + ev.class]
        .fetch_add(1, Ordering::Relaxed);
    SHADOW_BY_WIN[shadow_bin(ev.win_len, &SHADOW_WIN_EDGES) * SHADOW_CLASSES + ev.class]
        .fetch_add(1, Ordering::Relaxed);
    let slack_bin = (ev.min_slack / 8).clamp(0, SHADOW_SLACK_BINS as i64 - 1) as usize;
    SHADOW_SLACK[row * SHADOW_SLACK_BINS + slack_bin].fetch_add(1, Ordering::Relaxed);
    SHADOW_CELLS.fetch_add(ev.full_cells, Ordering::Relaxed);
    SHADOW_TOTAL.fetch_add(1, Ordering::Relaxed);
}

/// Test-only raw-cell counter read-back: `trigger` 0 = edge, 1 = score.
#[cfg(test)]
pub fn shadow_count(site: usize, trigger: usize, class: usize) -> u64 {
    SHADOW_BY_SITE[(shadow_site_row(site) * SHADOW_TRIGGERS + trigger) * SHADOW_CLASSES + class]
        .load(Ordering::Relaxed)
}

/// One `[shadow] {\u{2026}}` JSON line (empty string when ESP_SHADOW is
/// unset), printed at the same points as [`report_line`]. Class arrays
/// are `[A, B, C]`.
pub fn shadow_lines() -> String {
    if !shadow_enabled() {
        return String::new();
    }
    let cls = |slice: &[AtomicU64]| {
        format!(
            "[{}, {}, {}]",
            slice[0].load(Ordering::Relaxed),
            slice[1].load(Ordering::Relaxed),
            slice[2].load(Ordering::Relaxed)
        )
    };
    let mut s = format!(
        "[shadow] {{\"total\": {}",
        SHADOW_TOTAL.load(Ordering::Relaxed)
    );
    s.push_str(", \"site\": {");
    let mut first = true;
    for (row, chunk) in SHADOW_BY_SITE
        .chunks(SHADOW_TRIGGERS * SHADOW_CLASSES)
        .enumerate()
    {
        if chunk.iter().map(|c| c.load(Ordering::Relaxed)).sum::<u64>() == 0 {
            continue;
        }
        if !first {
            s.push_str(", ");
        }
        first = false;
        let name = if row == SHADOW_SITE_UNSET {
            "site_unset"
        } else {
            NAMES[row]
        };
        s.push_str(&format!(
            "\"{}\": {{\"e\": {}, \"s\": {}}}",
            name,
            cls(&chunk[..SHADOW_CLASSES]),
            cls(&chunk[SHADOW_CLASSES..])
        ));
    }
    s.push('}');
    let bin_key = |i: usize, edges: &[i64; 9]| {
        let lo = if i == 0 { 0 } else { edges[i - 1] };
        match edges.get(i) {
            Some(&hi) => format!("{lo}-{}", hi - 1),
            None => format!("{lo}+"),
        }
    };
    s.push_str(&format!(
        ", \"cells\": {}",
        SHADOW_CELLS.load(Ordering::Relaxed)
    ));
    s.push_str(", \"slack\": {");
    let mut first = true;
    for (row, chunk) in SHADOW_SLACK.chunks(SHADOW_SLACK_BINS).enumerate() {
        if chunk.iter().map(|c| c.load(Ordering::Relaxed)).sum::<u64>() == 0 {
            continue;
        }
        if !first {
            s.push_str(", ");
        }
        first = false;
        let name = if row == SHADOW_SITE_UNSET { "site_unset" } else { NAMES[row] };
        let buckets: Vec<String> = chunk
            .iter()
            .enumerate()
            .map(|(i, c)| {
                let lo = i * 8;
                let key = if i == SHADOW_SLACK_BINS - 1 {
                    format!("{lo}+")
                } else {
                    format!("{lo}-{}", lo + 7)
                };
                format!("\"{key}\": {}", c.load(Ordering::Relaxed))
            })
            .collect();
        s.push_str(&format!("\"{}\": {{{}}}", name, buckets.join(", ")));
    }
    s.push('}');
    for (label, edges, table) in [
        ("readlen", &SHADOW_READ_EDGES, &SHADOW_BY_READ),
        ("winlen", &SHADOW_WIN_EDGES, &SHADOW_BY_WIN),
    ] {
        s.push_str(&format!(", \"{label}\": {{"));
        for i in 0..SHADOW_LEN_BINS {
            if i > 0 {
                s.push_str(", ");
            }
            s.push_str(&format!(
                "\"{}\": {}",
                bin_key(i, edges),
                cls(&table[i * SHADOW_CLASSES..(i + 1) * SHADOW_CLASSES])
            ));
        }
        s.push('}');
    }
    s.push('}');
    s
}
