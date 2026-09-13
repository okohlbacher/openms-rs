// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `SYSTEM/SysInfo.h`: byte formatting, memory readings and `MemUsage`.
//!
//! The class test has two sections. The first is a table of exact literals for
//! `bytesToHumanReadable`; the second loads a 20 MB mzML and asserts that the
//! reported working set grew by more than 10 000 KB, which is reproduced here
//! with a deliberately touched allocation instead of a file, so the test needs
//! no fixture and measures only the allocator and the kernel.

use openms::system::sys_info::{
    MAX_EVENT_LABEL_BYTES, MemUsage, bytes_to_human_readable, free_system_memory, process_id,
    process_memory_consumption, process_peak_memory_consumption,
};
use std::sync::{Mutex, MutexGuard};

/// 64 MiB — 65 536 KiB against the 10 000 KB growth the source section asserts,
/// so more than six times the margin it asks for.
const PROBE_BYTES: usize = 64 * 1024 * 1024;

/// A working-set reading is a property of the whole process, and `cargo test`
/// runs every test in this binary as a thread of one process. Any test that
/// measures memory or allocates enough to move the reading takes this lock, so
/// one test's 64 MiB block cannot be released in the middle of another's
/// measurement. `cargo nextest` gives each test its own process and does not
/// need it; the lock costs nothing there.
static MEMORY: Mutex<()> = Mutex::new(());

fn memory_guard() -> MutexGuard<'static, ()> {
    MEMORY
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}

/// Force every page of a large allocation to be resident.
fn resident_allocation() -> Vec<u8> {
    let mut block = vec![0u8; PROBE_BYTES];
    for page in (0..PROBE_BYTES).step_by(4096) {
        block[page] = 1;
    }
    std::hint::black_box(block)
}

/// Class-test section `std::string bytesToHumanReadable(UInt64 bytes)`.
///
/// All six literals are the source test's own expectations, transcribed.
#[test]
fn bytes_to_human_readable_reproduces_every_class_test_literal() {
    // The source writes these as `2ull << 00` and `2048ull << 00`.
    assert_eq!(bytes_to_human_readable(2), "2 byte");
    assert_eq!(bytes_to_human_readable(2048), "2 KiB");
    assert_eq!(bytes_to_human_readable(2048 << 10), "2 MiB");
    assert_eq!(bytes_to_human_readable(2048 << 20), "2 GiB");
    assert_eq!(bytes_to_human_readable(2048u64 << 30), "2 TiB");
    assert_eq!(bytes_to_human_readable(2048u64 << 40), "2 PiB");
}

/// Unit boundaries and the four-significant-digit format, derived from the
/// source expression `std::setprecision(4)` on a value below 1024.
#[test]
fn bytes_to_human_readable_switches_units_at_1024_with_four_digits() {
    assert_eq!(bytes_to_human_readable(0), "0 byte");
    assert_eq!(bytes_to_human_readable(1), "1 byte");
    assert_eq!(bytes_to_human_readable(1000), "1000 byte");
    assert_eq!(bytes_to_human_readable(1023), "1023 byte");
    assert_eq!(bytes_to_human_readable(1024), "1 KiB");
    assert_eq!(bytes_to_human_readable(1536), "1.5 KiB");
    // 45.34 MiB, the example in the source's own doc comment.
    assert_eq!(
        bytes_to_human_readable((45.34 * 1024.0 * 1024.0) as u64),
        "45.34 MiB"
    );
    // Six units are not enough from 2^60 bytes upwards.
    assert_eq!(
        bytes_to_human_readable(1u64 << 60),
        "Congrats. That's a lot of bytes: 1152921504606846976"
    );
    assert_eq!(
        bytes_to_human_readable(u64::MAX),
        "Congrats. That's a lot of bytes: 18446744073709551615"
    );
}

/// Class-test section `static bool getProcessMemoryConsumption(size_t& mem_virtual)`.
///
/// The source asserts that the reading grows by more than 10 000 KB after
/// loading a 20 MB file. Here a 64 MiB allocation is made resident instead.
#[cfg(target_os = "linux")]
#[test]
fn process_memory_consumption_grows_with_a_resident_allocation() {
    let _guard = memory_guard();
    let before = process_memory_consumption().expect("Linux reports VmRSS");
    assert!(before > 0);

    let block = resident_allocation();
    let after = process_memory_consumption().expect("Linux reports VmRSS");
    let growth = i128::from(after) - i128::from(before);
    assert!(
        growth > 10_000,
        "expected more than 10000 KB growth, saw {before} -> {after} KB"
    );
    drop(block);

    // The source's own comment: the memory need not come back to the kernel.
    let released = process_memory_consumption().expect("Linux reports VmRSS");
    assert!(released > 0);
}

/// The same section on a platform this port does not read memory on.
#[cfg(not(target_os = "linux"))]
#[test]
fn process_memory_consumption_is_explicitly_unavailable() {
    assert_eq!(process_memory_consumption(), None);
    assert_eq!(process_peak_memory_consumption(), None);
    assert_eq!(free_system_memory(), None);
}

/// Peak working set never falls below the current working set.
#[cfg(target_os = "linux")]
#[test]
fn peak_memory_is_at_least_the_current_memory() {
    let _guard = memory_guard();
    let block = resident_allocation();
    let current = process_memory_consumption().expect("Linux reports VmRSS");
    let peak = process_peak_memory_consumption().expect("Linux reports VmHWM");
    assert!(peak >= current, "peak {peak} KB below current {current} KB");
    drop(block);
}

/// Available physical memory is a positive figure below the machine's total.
#[cfg(target_os = "linux")]
#[test]
fn free_system_memory_is_positive() {
    let available = free_system_memory().expect("Linux publishes MemAvailable or MemFree");
    assert!(available > 0);
}

/// The process identifier is the one the operating system gave this process.
#[test]
fn process_id_is_stable_and_non_zero() {
    let first = process_id();
    assert!(first > 0);
    assert_eq!(first, process_id());
    assert_eq!(first, std::process::id());
}

/// `MemUsage`: construction records the first point, `delta` the second.
#[test]
fn mem_usage_records_both_time_points_and_reports_a_delta() {
    let _guard = memory_guard();
    let mut usage = MemUsage::new();
    assert_eq!(usage.mem_after, None);
    assert_eq!(usage.mem_after_peak, None);

    let block = resident_allocation();
    let report = usage.delta("loading").unwrap();
    assert!(report.starts_with("Memory usage (loading): "));
    assert!(report.contains("(working set delta)"));
    drop(block);

    // `delta` recorded the second point, so a second call does not resample.
    let recorded = usage.mem_after;
    let repeated = usage.delta("loading").unwrap();
    assert_eq!(usage.mem_after, recorded);
    assert_eq!(repeated, report);
}

/// `MemUsage::delta` reports growth in whole mebibytes on Linux.
#[cfg(target_os = "linux")]
#[test]
fn mem_usage_delta_counts_a_resident_allocation_in_mebibytes() {
    let _guard = memory_guard();
    let mut usage = MemUsage::new();
    let block = resident_allocation();
    usage.after();
    drop(block);

    let before = usage.mem_before.expect("Linux reports VmRSS");
    let after = usage.mem_after.expect("Linux reports VmRSS");
    let expected = (i128::from(after) - i128::from(before)) / 1024;
    assert!(expected >= 60, "expected about 64 MiB, saw {expected} MB");
    let report = usage.delta("probe").unwrap();
    assert!(
        report.starts_with(&format!(
            "Memory usage (probe): {expected} MB (working set delta)"
        )),
        "unexpected report: {report}"
    );
}

/// `MemUsage::usage` reports the absolute second reading, not a difference.
#[test]
fn mem_usage_usage_reports_the_absolute_second_reading() {
    let mut usage = MemUsage::new();
    let report = usage.usage();
    assert!(report.starts_with("Memory usage: "));
    assert!(report.contains("(working set)"));
    assert!(!report.contains("delta"));
}

/// `MemUsage::reset` forgets everything, including the first time point.
#[test]
fn mem_usage_reset_clears_all_four_readings() {
    let _guard = memory_guard();
    let mut usage = MemUsage::new();
    usage.after();
    usage.reset();
    assert_eq!(usage.mem_before, None);
    assert_eq!(usage.mem_before_peak, None);
    assert_eq!(usage.mem_after, None);
    assert_eq!(usage.mem_after_peak, None);

    // `before` records the first point again, and records it exactly where the
    // platform has a reading. The value itself is not compared against a second
    // sample: a working set is a live figure and two samples need not agree.
    usage.before();
    assert_eq!(
        usage.mem_before.is_some(),
        process_memory_consumption().is_some()
    );
    assert_eq!(usage.mem_after, None);
}

/// A reading the platform does not provide prints as `n/a`, never as a number.
#[test]
fn unavailable_readings_print_as_not_available() {
    let mut usage = MemUsage {
        mem_before: None,
        mem_before_peak: None,
        mem_after: Some(4096),
        mem_after_peak: None,
    };
    assert_eq!(
        usage.delta("io").unwrap(),
        "Memory usage (io): n/a (working set delta)"
    );
    assert_eq!(usage.usage(), "Memory usage: 4 MB (working set)");
}

/// A shrinking working set keeps its sign, which the source drops.
#[test]
fn a_negative_delta_keeps_its_sign() {
    let mut usage = MemUsage {
        mem_before: Some(4096),
        mem_before_peak: Some(8192),
        mem_after: Some(1024),
        mem_after_peak: Some(8192),
    };
    assert_eq!(
        usage.delta("release").unwrap(),
        "Memory usage (release): -3 MB (working set delta), 0 MB (peak working set delta)"
    );
}

/// A zero peak reading is treated as unsupported, as the source's `> 0` test does.
#[test]
fn a_zero_peak_reading_is_omitted_from_the_report() {
    let mut usage = MemUsage {
        mem_before: Some(0),
        mem_before_peak: Some(0),
        mem_after: Some(2048),
        mem_after_peak: Some(0),
    };
    assert_eq!(
        usage.delta("step").unwrap(),
        "Memory usage (step): 2 MB (working set delta)"
    );
    assert_eq!(usage.usage(), "Memory usage: 2 MB (working set)");
}

/// The bound on the caller-supplied event label, which the source lacks.
#[test]
fn an_oversized_event_label_is_refused() {
    let _guard = memory_guard();
    let mut usage = MemUsage::new();
    assert!(usage.delta(&"x".repeat(MAX_EVENT_LABEL_BYTES + 1)).is_err());
    assert!(usage.delta(&"x".repeat(MAX_EVENT_LABEL_BYTES)).is_ok());
}
