// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! OpenSwathWorkflowScheduler.
//!
//! The wave planning and batch sizing cases are transcribed from the upstream
//! class test `OpenSwathWorkflowScheduler_test.cpp` at the core pin `bc9cc12`
//! (evidence tier 3). The memory arithmetic, which the upstream test does not
//! reach because it reads the machine's free memory, is checked against values
//! derived from the source formula by hand (tier 4); the concurrency limiter
//! is checked with real threads.

use openms::analysis::openswath::openswath_workflow_scheduler::{
    ConcurrencyEstimate, ConcurrencyLimiter, OpenSwathWorkflowScheduler as Scheduler,
    SchedulerOptions, ScopedSlot, Wave,
};
use openms::analysis::openswath::swath_map::SwathMap;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

const MIB: u64 = 1024 * 1024;

/// One leading MS1 map, then `n_ms2` non-MS1 maps, with a second MS1 map
/// spliced into the middle (index `n_ms2 / 2 + 1`).
fn make_maps(n_ms2: usize) -> Vec<SwathMap> {
    let mut maps = vec![SwathMap::new(0.0, 0.0, 0.0, true)];
    let mut lo = 400.0;
    for i in 0..n_ms2 {
        if i == n_ms2 / 2 {
            maps.push(SwathMap::new(0.0, 0.0, 0.0, true));
        }
        maps.push(SwathMap::new(lo, lo + 25.0, lo + 12.5, false));
        lo += 25.0;
    }
    maps
}

fn non_ms1_indices(maps: &[SwathMap]) -> Vec<usize> {
    maps.iter()
        .enumerate()
        .filter(|(_, m)| !m.ms1)
        .map(|(i, _)| i)
        .collect()
}

fn estimate(max_concurrent: usize, budget: u64) -> ConcurrencyEstimate {
    ConcurrencyEstimate {
        max_concurrent_swaths: max_concurrent,
        memory_budget_bytes: budget,
        avg_spectra_per_swath: 0,
        ..ConcurrencyEstimate::default()
    }
}

/// A valid partition: every non-MS1 index once, in input order, no MS1 index,
/// no wave above the cap and no empty wave.
fn check_partition(maps: &[SwathMap], est: &ConcurrencyEstimate) -> Vec<Wave> {
    let waves = Scheduler::plan_waves(maps, &SchedulerOptions::default(), est);
    let cap = est.max_concurrent_swaths.max(1);
    let mut collected = Vec::new();
    for wave in &waves {
        assert!(wave.swath_indices.len() <= cap);
        assert!(!wave.swath_indices.is_empty());
        for &idx in &wave.swath_indices {
            assert!(!maps[idx].ms1, "MS1 map {idx} in a wave");
            collected.push(idx);
        }
    }
    assert_eq!(collected, non_ms1_indices(maps));
    waves
}

#[test]
fn plan_waves_count_cap() {
    let maps = make_maps(7);
    assert_eq!(non_ms1_indices(&maps).len(), 7);
    // memory budget 0 is unbounded: only the count cap closes waves
    assert_eq!(check_partition(&maps, &estimate(1, 0)).len(), 7);
    assert_eq!(check_partition(&maps, &estimate(3, 0)).len(), 3);
    let waves = check_partition(&maps, &estimate(100, 0));
    assert_eq!(waves.len(), 1);
    assert_eq!(waves[0].swath_indices.len(), 7);
}

#[test]
fn plan_waves_memory_cap() {
    // with no spectra each SWATH costs the 100 MiB overhead; a 250 MiB budget
    // fits two (three would overflow), so 7 maps give 4 waves
    let maps = make_maps(7);
    let waves = check_partition(&maps, &estimate(100, 250 * MIB));
    assert_eq!(waves.len(), 4);
    assert_eq!(waves[0].swath_indices.len(), 2);
    assert_eq!(waves[3].swath_indices.len(), 1);
    assert_eq!(waves[0].estimated_bytes, 200 * MIB);
    assert_eq!(waves[3].estimated_bytes, 100 * MIB);
}

#[test]
fn plan_waves_without_non_ms1_maps() {
    let mut ms1_only = vec![SwathMap::new(0.0, 0.0, 0.0, true); 2];
    ms1_only[1].ms1 = true;
    let est = ConcurrencyEstimate {
        max_concurrent_swaths: 4,
        ..ConcurrencyEstimate::default()
    };
    assert!(Scheduler::plan_waves(&ms1_only, &SchedulerOptions::default(), &est).is_empty());
    assert!(Scheduler::plan_waves(&[], &SchedulerOptions::default(), &est).is_empty());
}

#[test]
fn plan_waves_transition_term_and_avg_fallback() {
    // 4 spectra assumed (avg fallback, no data attached), 10 transitions per
    // SWATH at 64 bytes a point: 4 * 600 KiB + 100 MiB + 4 * 10 * 64
    let maps = make_maps(2);
    let options = SchedulerOptions {
        avg_transitions_per_swath: 10,
        ..SchedulerOptions::default()
    };
    let est = ConcurrencyEstimate {
        max_concurrent_swaths: 8,
        avg_spectra_per_swath: 4,
        ..ConcurrencyEstimate::default()
    };
    let waves = Scheduler::plan_waves(&maps, &options, &est);
    assert_eq!(waves.len(), 1);
    let per_swath = 4 * 600 * 1024 + 100 * MIB + 4 * 10 * 64;
    assert_eq!(waves[0].estimated_bytes, 2 * per_swath);
}

#[test]
fn choose_inner_batch_size_upstream_cases() {
    let options = SchedulerOptions::default();
    let choose = |total, active, threads, user| {
        Scheduler::choose_inner_batch_size(total, active, threads, user, &options)
    };
    // nothing to batch
    assert_eq!(choose(0, 1, 4, -1), 0);
    // a user size short-circuits to min(user, total)
    assert_eq!(choose(10000, 1, 4, 500), 500);
    assert_eq!(choose(1000, 1, 4, 50000), 1000);
    // target_jobs = max(active, threads * 3); ceil(total / jobs) clamped to
    // [2000, 10000], then capped at the total
    assert_eq!(choose(60000, 1, 4, -1), 5000);
    assert_eq!(choose(1_000_000, 1, 4, -1), 10000);
    assert_eq!(choose(1000, 1, 4, -1), 1000);
    assert_eq!(choose(200_000, 50, 4, -1), 4000);
    // zero threads and zero active SWATHs count as one
    assert_eq!(choose(30000, 0, 0, -1), 10000);
}

#[test]
fn available_memory_arithmetic() {
    // 16 GiB free, default reserves 2 GiB + 500 MiB + 2 GiB OSW buffer,
    // 90% of the rest
    let options = SchedulerOptions::default();
    let free_kib = 16 * 1024 * 1024;
    let available = 16 * 1024 * MIB;
    let reserved = 2 * 1024 * MIB + 500 * MIB + 2 * 1024 * MIB;
    let expected = ((available - reserved) as f64 * 0.90) as u64;
    assert_eq!(
        Scheduler::estimate_available_memory_from(Some(free_kib), &options),
        expected
    );
    // reserves above the free memory, and an unreadable free memory, give 0
    assert_eq!(
        Scheduler::estimate_available_memory_from(Some(4 * 1024 * 1024), &options),
        0
    );
    assert_eq!(Scheduler::estimate_available_memory_from(None, &options), 0);
    // the usage fraction is clamped to [0.05, 0.95]
    let greedy = SchedulerOptions {
        memory_usage_fraction: 5.0,
        ..SchedulerOptions::default()
    };
    assert_eq!(
        Scheduler::estimate_available_memory_from(Some(free_kib), &greedy),
        ((available - reserved) as f64 * 0.95) as u64
    );
    let timid = SchedulerOptions {
        memory_usage_fraction: -1.0,
        ..SchedulerOptions::default()
    };
    assert_eq!(
        Scheduler::estimate_available_memory_from(Some(free_kib), &timid),
        ((available - reserved) as f64 * 0.05) as u64
    );
}

#[test]
fn estimate_concurrency_rules() {
    let options = SchedulerOptions::default();
    // no non-MS1 map with data: the default estimate
    let maps = make_maps(3); // maps without data attached do not count
    let est = Scheduler::estimate_concurrency_with_budget(&maps, &options, || 10 * 1024 * MIB);
    assert_eq!(est, ConcurrencyEstimate::default());
    assert_eq!(est.max_concurrent_swaths, 1);
    // the budget is not even queried then
    let queried = AtomicUsize::new(0);
    let _ = Scheduler::estimate_concurrency_with_budget(&maps, &options, || {
        queried.fetch_add(1, Ordering::SeqCst);
        0
    });
    assert_eq!(queried.load(Ordering::SeqCst), 0);
}

#[test]
fn concurrency_limiter_blocks_at_capacity_and_releases() {
    // null guard: a no-op
    {
        let _slot = ScopedSlot::new(None);
    }
    // within capacity nothing blocks
    let limiter = ConcurrencyLimiter::new(2);
    {
        let _s1 = ScopedSlot::new(Some(&limiter));
        {
            let _s2 = ScopedSlot::new(Some(&limiter));
        }
        let _s3 = ScopedSlot::new(Some(&limiter));
    }
    // capacity 0 is clamped to 1
    let clamped = ConcurrencyLimiter::new(0);
    assert_eq!(clamped.capacity(), 1);
    {
        let _s = ScopedSlot::new(Some(&clamped));
    }
    // releasing without a slot is a no-op, and the limiter still works
    clamped.release_slot();
    clamped.release_slot();
    clamped.acquire_slot();
    clamped.release_slot();
}

#[test]
fn concurrency_limiter_never_exceeds_capacity_across_threads() {
    let limiter = Arc::new(ConcurrencyLimiter::new(2));
    let active = Arc::new(AtomicUsize::new(0));
    let peak = Arc::new(AtomicUsize::new(0));
    let handles: Vec<_> = (0..8)
        .map(|_| {
            let limiter = Arc::clone(&limiter);
            let active = Arc::clone(&active);
            let peak = Arc::clone(&peak);
            std::thread::spawn(move || {
                for _ in 0..50 {
                    let _slot = ScopedSlot::new(Some(&limiter));
                    let now = active.fetch_add(1, Ordering::SeqCst) + 1;
                    peak.fetch_max(now, Ordering::SeqCst);
                    std::thread::yield_now();
                    active.fetch_sub(1, Ordering::SeqCst);
                }
            })
        })
        .collect();
    for handle in handles {
        handle.join().unwrap();
    }
    assert!(peak.load(Ordering::SeqCst) <= 2);
    assert_eq!(active.load(Ordering::SeqCst), 0);
}
