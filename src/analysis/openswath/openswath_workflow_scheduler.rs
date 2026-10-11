// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Memory-bounded work planning for the OpenSWATH workflow.
//!
//! Ports `ANALYSIS/OPENSWATH/OpenSwathWorkflowScheduler.h` with its `.cpp`.
//! The scheduler keeps the policy for memory-aware SWATH concurrency outside
//! the extraction and scoring code: it estimates how many SWATH maps may be
//! resident at once ([`OpenSwathWorkflowScheduler::estimate_concurrency`]),
//! groups the non-MS1 maps into waves that are loaded as a unit
//! ([`OpenSwathWorkflowScheduler::plan_waves`]) and sizes the scoring batches
//! that keep the workers busy
//! ([`OpenSwathWorkflowScheduler::choose_inner_batch_size`]). The class is a
//! namespace in disguise in the source, with two nested classes: the slot
//! limiter [`ConcurrencyLimiter`] and its guard [`ScopedSlot`].
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_WORKFLOW_SUPPORT.md).
//!
//! # Example
//!
//! ```
//! use openms::analysis::openswath::openswath_workflow_scheduler::{
//!     OpenSwathWorkflowScheduler, SchedulerOptions,
//! };
//!
//! let options = SchedulerOptions::default();
//! // total 60000 compounds on 4 threads: 12 target jobs of 5000 each
//! assert_eq!(
//!     OpenSwathWorkflowScheduler::choose_inner_batch_size(60000, 1, 4, -1, &options),
//!     5000
//! );
//! ```

use super::swath_map::SwathMap;
use std::sync::{Condvar, Mutex};

/// The free system memory in KiB, or `None` when the platform does not report
/// it here.
///
/// Linux reads `MemAvailable` of `/proc/meminfo`, falling back to `MemFree`
/// (the source's own primary path); other platforms report nothing, as
/// `system::sys_info::free_system_memory` does. That function is the same query:
/// this module cannot import it, since a dependency of `analysis` on `system`
/// would close the cycle `analysis` - `system` - `format` - `analysis`
/// (`tools/check_module_cycles.py`).
pub fn free_memory_kib() -> Option<u64> {
    #[cfg(target_os = "linux")]
    {
        let text = std::fs::read_to_string("/proc/meminfo").ok()?;
        let field = |key: &str| {
            text.lines().find_map(|line| {
                let rest = line.strip_prefix(key)?;
                rest.split_whitespace().next()?.parse::<u64>().ok()
            })
        };
        field("MemAvailable:").or_else(|| field("MemFree:"))
    }
    #[cfg(not(target_os = "linux"))]
    {
        None
    }
}

const GIB: u64 = 1024 * 1024 * 1024;
const MIB: u64 = 1024 * 1024;

/// The memory and concurrency knobs of wave planning and batch sizing
/// (source `OpenSwathWorkflowScheduler::Options`). All quantities are byte
/// counts unless noted.
#[derive(Clone, Debug, PartialEq)]
pub struct SchedulerOptions {
    /// OSW writer memory reserved outside the scheduler's budget. Default 2 GiB.
    pub osw_buffer_bytes: u64,
    /// Fraction of the remaining free memory the scheduler may reserve;
    /// clamped to [0.05, 0.95] by
    /// [`estimate_available_memory_for_scoring`](OpenSwathWorkflowScheduler::estimate_available_memory_for_scoring).
    /// Default 0.90.
    pub memory_usage_fraction: f64,
    /// Estimated memory of one cached spectrum. Default 600 KiB.
    pub bytes_per_spectrum: u64,
    /// Fixed estimated overhead per loaded SWATH. Default 100 MiB.
    pub per_swath_overhead_bytes: u64,
    /// Average number of chromatograms extracted per non-MS1 SWATH; 0 drops
    /// the transition-density term from the per-SWATH estimate.
    pub avg_transitions_per_swath: usize,
    /// Estimated bytes retained per extracted chromatogram point. Default 64.
    pub bytes_per_chromatogram_point: u64,
    /// User cap on concurrent SWATHs: a positive value caps the estimate, a
    /// value of zero or below plans by memory, bounded only by
    /// [`scoring_threads`](Self::scoring_threads). Default -1.
    pub max_concurrent_swaths: i32,
    /// The number of scoring worker threads. Default 1.
    pub scoring_threads: usize,
    /// Smallest scoring batch size picked automatically. Default 2000.
    pub min_inner_batch_size: usize,
    /// Largest scoring batch size picked automatically. Default 10000.
    pub max_inner_batch_size: usize,
    /// Target number of queued scoring jobs per worker thread. Default 3.
    pub target_jobs_per_thread: usize,
}

impl Default for SchedulerOptions {
    fn default() -> Self {
        Self {
            osw_buffer_bytes: 2 * GIB,
            memory_usage_fraction: 0.90,
            bytes_per_spectrum: 600 * 1024,
            per_swath_overhead_bytes: 100 * MIB,
            avg_transitions_per_swath: 0,
            bytes_per_chromatogram_point: 64,
            max_concurrent_swaths: -1,
            scoring_threads: 1,
            min_inner_batch_size: 2000,
            max_inner_batch_size: 10000,
            target_jobs_per_thread: 3,
        }
    }
}

/// The derived limits of one extraction and scoring run (source
/// `OpenSwathWorkflowScheduler::ConcurrencyEstimate`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConcurrencyEstimate {
    /// Non-MS1 SWATH maps with spectra attached.
    pub non_ms1_swath_count: usize,
    /// Mean spectra per non-MS1 SWATH (integer division).
    pub avg_spectra_per_swath: usize,
    /// The cap on simultaneously resident non-MS1 SWATH maps; at least 1.
    pub max_concurrent_swaths: usize,
    /// The memory budget; 0 when it could not be determined, which
    /// [`plan_waves`](OpenSwathWorkflowScheduler::plan_waves) treats as unbounded.
    pub memory_budget_bytes: u64,
    /// The predicted memory of one non-MS1 SWATH.
    pub estimated_bytes_per_swath: u64,
}

impl Default for ConcurrencyEstimate {
    fn default() -> Self {
        Self {
            non_ms1_swath_count: 0,
            avg_spectra_per_swath: 0,
            max_concurrent_swaths: 1,
            memory_budget_bytes: 0,
            estimated_bytes_per_swath: 0,
        }
    }
}

/// A group of SWATH indices that may be resident at the same time (source
/// `OpenSwathWorkflowScheduler::Wave`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Wave {
    /// Indices into the SWATH map vector; never an MS1 map.
    pub swath_indices: Vec<usize>,
    /// The sum of the per-SWATH byte estimates of the indices.
    pub estimated_bytes: u64,
}

fn saturating_mul(lhs: u64, rhs: u64) -> u64 {
    lhs.saturating_mul(rhs)
}

fn transition_chromatogram_bytes(
    spectra_count: u64,
    avg_transitions_per_swath: usize,
    bytes_per_chromatogram_point: u64,
) -> u64 {
    if avg_transitions_per_swath == 0 || bytes_per_chromatogram_point == 0 {
        return 0;
    }
    let points = saturating_mul(spectra_count, avg_transitions_per_swath as u64);
    saturating_mul(points, bytes_per_chromatogram_point)
}

/// The planner (source `OpenSwathWorkflowScheduler`, all static members).
#[derive(Clone, Copy, Debug, Default)]
pub struct OpenSwathWorkflowScheduler;

impl OpenSwathWorkflowScheduler {
    /// The memory available for scoring after the OS, I/O and OSW writer
    /// reserves (source `estimateAvailableMemoryForScoring`).
    ///
    /// Queries the free system memory, subtracts 2 GiB for the OS, 500 MiB for
    /// I/O buffers and [`SchedulerOptions::osw_buffer_bytes`], and multiplies
    /// the rest by the usage fraction clamped to [0.05, 0.95]. Returns 0 when
    /// the memory cannot be queried or the reserves exceed it.
    pub fn estimate_available_memory_for_scoring(options: &SchedulerOptions) -> u64 {
        Self::estimate_available_memory_from(free_memory_kib(), options)
    }

    /// As [`estimate_available_memory_for_scoring`](Self::estimate_available_memory_for_scoring)
    /// for a given amount of free memory in KiB (`None`: not reported), so
    /// the arithmetic can be checked without a particular machine.
    pub fn estimate_available_memory_from(
        free_kib: Option<u64>,
        options: &SchedulerOptions,
    ) -> u64 {
        let Some(available_kb) = free_kib else {
            return 0;
        };
        let available_bytes = available_kb.saturating_mul(1024);
        let os_reserved = 2 * GIB;
        let io_buffer = 500 * MIB;
        let reserved = os_reserved
            .saturating_add(io_buffer)
            .saturating_add(options.osw_buffer_bytes);
        if available_bytes <= reserved {
            return 0;
        }
        let fraction = options.memory_usage_fraction.clamp(0.05, 0.95);
        ((available_bytes - reserved) as f64 * fraction) as u64
    }

    /// The per-run concurrency limits of a SWATH map vector (source
    /// `estimateConcurrency`).
    ///
    /// Counts the non-MS1 maps with spectra attached, takes the mean spectra
    /// count and estimates one SWATH as `avg_spectra * bytes_per_spectrum +
    /// per_swath_overhead + transition term`, the latter only when both
    /// `avg_transitions_per_swath` and `bytes_per_chromatogram_point` are
    /// non-zero. The cap is the minimum of the SWATH count, the memory budget
    /// over the estimate (skipped when either is zero or the budget does not
    /// exceed the estimate) and the user cap or `max(1, scoring_threads)`,
    /// clamped to at least 1. Without non-MS1 maps the default estimate
    /// (everything zero, cap 1) is returned.
    pub fn estimate_concurrency(
        swath_maps: &[SwathMap],
        options: &SchedulerOptions,
    ) -> ConcurrencyEstimate {
        Self::estimate_concurrency_with_budget(swath_maps, options, || {
            Self::estimate_available_memory_for_scoring(options)
        })
    }

    /// As [`estimate_concurrency`](Self::estimate_concurrency) with the memory
    /// budget supplied by the caller (queried only when there is a non-MS1
    /// map, as in the source).
    pub fn estimate_concurrency_with_budget(
        swath_maps: &[SwathMap],
        options: &SchedulerOptions,
        budget: impl FnOnce() -> u64,
    ) -> ConcurrencyEstimate {
        let mut estimate = ConcurrencyEstimate::default();
        let mut total_spectra: u64 = 0;
        for swath_map in swath_maps {
            if !swath_map.ms1 {
                if let Some(sptr) = &swath_map.sptr {
                    total_spectra += sptr.nr_spectra() as u64;
                    estimate.non_ms1_swath_count += 1;
                }
            }
        }
        if estimate.non_ms1_swath_count == 0 {
            return estimate;
        }
        estimate.memory_budget_bytes = budget();
        estimate.avg_spectra_per_swath =
            (total_spectra / estimate.non_ms1_swath_count as u64) as usize;
        estimate.estimated_bytes_per_swath = (estimate.avg_spectra_per_swath as u64)
            .wrapping_mul(options.bytes_per_spectrum)
            .wrapping_add(options.per_swath_overhead_bytes)
            .wrapping_add(transition_chromatogram_bytes(
                estimate.avg_spectra_per_swath as u64,
                options.avg_transitions_per_swath,
                options.bytes_per_chromatogram_point,
            ));

        let mut max_concurrent = estimate.non_ms1_swath_count;
        if estimate.estimated_bytes_per_swath > 0
            && estimate.memory_budget_bytes > estimate.estimated_bytes_per_swath
        {
            max_concurrent = (estimate.memory_budget_bytes / estimate.estimated_bytes_per_swath)
                .min(usize::MAX as u64) as usize;
        }
        if options.max_concurrent_swaths > 0 {
            max_concurrent = max_concurrent.min(options.max_concurrent_swaths as usize);
        } else {
            max_concurrent = max_concurrent.min(options.scoring_threads.max(1));
        }
        estimate.max_concurrent_swaths = max_concurrent.min(estimate.non_ms1_swath_count).max(1);
        estimate
    }

    /// Group the non-MS1 SWATHs into waves that honour the count and memory
    /// caps (source `planWaves`).
    ///
    /// Walks the maps in input order, skips MS1 maps and appends each other
    /// map to the current wave. The wave is closed first when it already holds
    /// `estimate.max_concurrent_swaths` maps, or when it is non-empty and
    /// adding the map would push its bytes past `estimate.memory_budget_bytes`
    /// (0 is unbounded). A map costs `spectra * bytes_per_spectrum +
    /// per_swath_overhead + transition term`, from the actual spectra count
    /// of its data or, with none attached, `estimate.avg_spectra_per_swath`.
    pub fn plan_waves(
        swath_maps: &[SwathMap],
        options: &SchedulerOptions,
        estimate: &ConcurrencyEstimate,
    ) -> Vec<Wave> {
        let mut waves = Vec::new();
        let mut current = Wave::default();
        let max_concurrent = estimate.max_concurrent_swaths.max(1);
        let budget = if estimate.memory_budget_bytes > 0 {
            estimate.memory_budget_bytes
        } else {
            u64::MAX
        };
        for (index, swath_map) in swath_maps.iter().enumerate() {
            if swath_map.ms1 {
                continue;
            }
            let spectra = swath_map
                .sptr
                .as_ref()
                .map_or(estimate.avg_spectra_per_swath as u64, |s| {
                    s.nr_spectra() as u64
                });
            let swath_bytes = spectra
                .wrapping_mul(options.bytes_per_spectrum)
                .wrapping_add(options.per_swath_overhead_bytes)
                .wrapping_add(transition_chromatogram_bytes(
                    spectra,
                    options.avg_transitions_per_swath,
                    options.bytes_per_chromatogram_point,
                ));
            let full_by_count = current.swath_indices.len() >= max_concurrent;
            let full_by_memory = !current.swath_indices.is_empty()
                && current.estimated_bytes.wrapping_add(swath_bytes) > budget;
            if full_by_count || full_by_memory {
                waves.push(std::mem::take(&mut current));
            }
            current.swath_indices.push(index);
            current.estimated_bytes = current.estimated_bytes.wrapping_add(swath_bytes);
        }
        if !current.swath_indices.is_empty() {
            waves.push(current);
        }
        waves
    }

    /// A scoring batch size that keeps the worker pool occupied (source
    /// `chooseInnerBatchSize`).
    ///
    /// No compounds give 0. A positive `user_inner_batch_size` gives
    /// `min(user, total_compounds)`. Otherwise `max(active_swaths,
    /// scoring_threads * target_jobs_per_thread)` jobs are targeted (threads
    /// and active SWATHs count as at least 1), the batch is
    /// `ceil(total_compounds / jobs)` clamped to
    /// `[min_inner_batch_size, max_inner_batch_size]` and finally capped at
    /// `total_compounds`.
    pub fn choose_inner_batch_size(
        total_compounds: usize,
        active_swaths: usize,
        scoring_threads: usize,
        user_inner_batch_size: i32,
        options: &SchedulerOptions,
    ) -> usize {
        if total_compounds == 0 {
            return 0;
        }
        if user_inner_batch_size > 0 {
            return (user_inner_batch_size as usize).min(total_compounds);
        }
        let workers = scoring_threads.max(1);
        let active_jobs = active_swaths.max(1);
        let jobs_per_thread = options.target_jobs_per_thread.max(1);
        let target_jobs = active_jobs.max(workers.saturating_mul(jobs_per_thread));
        let mut batch = total_compounds.div_ceil(target_jobs);
        let min_batch = options.min_inner_batch_size.max(1);
        let max_batch = options.max_inner_batch_size.max(min_batch);
        batch = batch.clamp(min_batch, max_batch);
        batch.min(total_compounds)
    }
}

/// A counting semaphore for paths that stream one SWATH per worker (source
/// `OpenSwathWorkflowScheduler::ConcurrencyLimiter`).
///
/// The requested cap is clamped to at least 1. Pair it with [`ScopedSlot`]
/// instead of calling [`acquire_slot`](Self::acquire_slot) and
/// [`release_slot`](Self::release_slot) directly.
#[derive(Debug)]
pub struct ConcurrencyLimiter {
    max_concurrent_swaths: usize,
    active: Mutex<usize>,
    available: Condvar,
}

impl ConcurrencyLimiter {
    /// A limiter of at most `max_concurrent_swaths` slots; zero counts as one.
    pub fn new(max_concurrent_swaths: usize) -> Self {
        Self {
            max_concurrent_swaths: max_concurrent_swaths.max(1),
            active: Mutex::new(0),
            available: Condvar::new(),
        }
    }

    /// Block until a slot is free and take it.
    ///
    /// A poisoned lock (a panic while a slot was taken or released) is
    /// recovered, since the counter stays consistent.
    pub fn acquire_slot(&self) {
        let mut active = self.active.lock().unwrap_or_else(|e| e.into_inner());
        while *active >= self.max_concurrent_swaths {
            active = self
                .available
                .wait(active)
                .unwrap_or_else(|e| e.into_inner());
        }
        *active += 1;
    }

    /// Release a slot and wake one waiter. Safe (a no-op) when no slot is
    /// taken.
    pub fn release_slot(&self) {
        {
            let mut active = self.active.lock().unwrap_or_else(|e| e.into_inner());
            if *active > 0 {
                *active -= 1;
            }
        }
        self.available.notify_one();
    }

    /// The configured cap (at least 1).
    pub fn capacity(&self) -> usize {
        self.max_concurrent_swaths
    }
}

/// An RAII guard of a [`ConcurrencyLimiter`] slot (source
/// `OpenSwathWorkflowScheduler::ScopedSlot`): taken on construction, released
/// on drop. Without a limiter both are no-ops, so guard sites that may or may
/// not be limited need no branch.
#[derive(Debug)]
pub struct ScopedSlot<'a> {
    limiter: Option<&'a ConcurrencyLimiter>,
}

impl<'a> ScopedSlot<'a> {
    /// Acquire a slot of `limiter` (blocking); a `None` limiter does nothing.
    pub fn new(limiter: Option<&'a ConcurrencyLimiter>) -> Self {
        if let Some(limiter) = limiter {
            limiter.acquire_slot();
        }
        Self { limiter }
    }
}

impl Drop for ScopedSlot<'_> {
    fn drop(&mut self) {
        if let Some(limiter) = self.limiter {
            limiter.release_slot();
        }
    }
}
