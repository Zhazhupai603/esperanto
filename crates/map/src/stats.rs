//! Alignment statistics (align_qc.json contract).
//!
//! Consumers ignore unknown keys; `elapsed_seconds` is the only value exempt
//! from byte-level parity.

use crate::pair::InsertStats;
use serde::Serialize;

/// Final statistics document (serialized pretty with a trailing newline by
/// the pipeline).
#[derive(Debug, Serialize)]
pub struct AlignStats {
    /// Crate version (`CARGO_PKG_VERSION`).
    pub esperanto_map_version: String,
    /// Run mode label.
    pub mode: String,
    /// Total reads processed.
    pub total_reads: u64,
    /// Mapped reads.
    pub mapped_reads: u64,
    /// Unmapped reads.
    pub unmapped_reads: u64,
    /// mapped / total.
    pub mapping_rate: f64,
    /// Properly paired count (PE only; `null` for SE).
    pub proper_pairs: Option<u64>,
    /// Mean |tlen| (PE only; `null` for SE).
    pub insert_mean: Option<f64>,
    /// Population stdev of |tlen| (PE only; `null` for SE).
    pub insert_stdev: Option<f64>,
    /// MAPQ histogram, 61 buckets (0..=60).
    pub mapq_hist: Vec<u64>,
    /// Junctions reported (pass-2 library hits plus discoveries in use).
    pub junctions_total: u64,
    /// Reads rescued (split/span/Track-2/anchor-rescue).
    pub rescued_total: u64,
    /// Failed rescue attempts.
    pub rescue_fail_total: u64,
    /// Wall-clock seconds (parity-exempt).
    pub elapsed_seconds: f64,
}

/// Integer accumulator; `finalize` renders the JSON contract in one shot.
#[derive(Debug)]
pub struct StatsAcc {
    total: u64,
    mapped: u64,
    mapq_hist: [u64; 61],
    junctions: u64,
    rescued: u64,
    rescue_fail: u64,
    proper_pairs: u64,
    inserts: InsertStats,
}

impl Default for StatsAcc {
    fn default() -> StatsAcc {
        StatsAcc {
            total: 0,
            mapped: 0,
            mapq_hist: [0; 61],
            junctions: 0,
            rescued: 0,
            rescue_fail: 0,
            proper_pairs: 0,
            inserts: InsertStats::default(),
        }
    }
}

impl StatsAcc {
    /// Fresh accumulator.
    pub fn new() -> StatsAcc {
        StatsAcc::default()
    }

    /// Record one read outcome.
    pub fn push_read(&mut self, mapped: bool, mapq: u8) {
        self.total += 1;
        if mapped {
            self.mapped += 1;
            let bin = (mapq as usize).min(60);
            self.mapq_hist[bin] += 1;
        }
    }

    /// Record junctions reported for one alignment.
    pub fn push_junctions(&mut self, n: u64) {
        self.junctions += n;
    }

    /// Record one rescued read.
    pub fn push_rescued(&mut self) {
        self.rescued += 1;
    }

    /// Record one failed rescue attempt.
    pub fn push_rescue_fail(&mut self) {
        self.rescue_fail += 1;
    }

    /// Record one proper pair (PE runs only).
    pub fn push_proper_pair(&mut self) {
        self.proper_pairs += 1;
    }

    /// Record one template-length observation (PE runs only).
    pub fn push_insert(&mut self, tlen: i32) {
        self.inserts.push(tlen);
    }

    /// Render the statistics document. `pe = false` emits `null` for the
    /// pair-only fields.
    pub fn finalize(&self, mode: &str, pe: bool, elapsed_seconds: f64) -> AlignStats {
        let mapping_rate = if self.total == 0 {
            0.0
        } else {
            self.mapped as f64 / self.total as f64
        };
        AlignStats {
            esperanto_map_version: env!("CARGO_PKG_VERSION").to_string(),
            mode: mode.to_string(),
            total_reads: self.total,
            mapped_reads: self.mapped,
            unmapped_reads: self.total - self.mapped,
            mapping_rate,
            proper_pairs: if pe { Some(self.proper_pairs) } else { None },
            insert_mean: if pe { Some(self.inserts.mean()) } else { None },
            insert_stdev: if pe { Some(self.inserts.stdev()) } else { None },
            mapq_hist: self.mapq_hist.to_vec(),
            junctions_total: self.junctions,
            rescued_total: self.rescued,
            rescue_fail_total: self.rescue_fail,
            elapsed_seconds,
        }
    }
}
