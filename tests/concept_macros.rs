// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `CONCEPT/Macros.h`.
//!
//! There is no `Macros_test.cpp` at core `bc9cc12`. The critical-section tests
//! check the OpenMP `critical (name)` guarantees the macro expands to (tier 4,
//! from the OpenMP specification: same-named sections exclude each other,
//! differently named ones do not); the condition tests check that the macros
//! are `debug_assert!` counterparts, compiled in only with debug assertions.

use openms::concept::macros::thread_critical;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc;
use std::time::Duration;

#[test]
fn same_named_sections_exclude_each_other() {
    static INSIDE: AtomicUsize = AtomicUsize::new(0);
    static ENTERED: AtomicUsize = AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for _ in 0..8 {
            scope.spawn(|| {
                for _ in 0..200 {
                    thread_critical("concept_macros_exclusion", || {
                        assert_eq!(INSIDE.fetch_add(1, Ordering::SeqCst), 0);
                        ENTERED.fetch_add(1, Ordering::SeqCst);
                        std::hint::spin_loop();
                        INSIDE.fetch_sub(1, Ordering::SeqCst);
                    });
                }
            });
        }
    });
    assert_eq!(ENTERED.load(Ordering::SeqCst), 8 * 200);
}

#[test]
fn differently_named_sections_run_concurrently() {
    // Thread A holds section "a" until thread B, inside section "b", signals.
    // If the two names shared one lock, B could not enter and A's wait would
    // time out.
    let (to_a, from_b) = mpsc::channel();
    let (to_b, from_a) = mpsc::channel();
    std::thread::scope(|scope| {
        let a = scope.spawn(move || {
            thread_critical("concept_macros_a", || {
                to_b.send(()).ok();
                from_b.recv_timeout(Duration::from_secs(30)).is_ok()
            })
        });
        let b = scope.spawn(move || {
            from_a.recv_timeout(Duration::from_secs(30)).ok();
            thread_critical("concept_macros_b", || to_a.send(()).is_ok())
        });
        assert!(b.join().unwrap_or(false));
        assert!(a.join().unwrap_or(false));
    });
}

#[test]
fn a_section_returns_its_body_value_and_survives_a_panic() {
    assert_eq!(thread_critical("concept_macros_value", || 6 * 7), 42);
    let panicked = std::panic::catch_unwind(|| {
        thread_critical("concept_macros_panic", || panic!("inside the section"))
    });
    assert!(panicked.is_err());
    // The section is released and can be entered again.
    assert_eq!(thread_critical("concept_macros_panic", || "again"), "again");
}

#[test]
fn conditions_that_hold_do_nothing() {
    let values = [1.0, 2.0];
    openms::openms_precondition!(!values.is_empty(), "values must not be empty");
    openms::openms_postcondition!(values.len() == 2, "two values");
}

#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "Precondition failed: two == 3 arithmetic")]
fn a_failed_precondition_panics_with_debug_assertions() {
    let two = std::hint::black_box(2);
    openms::openms_precondition!(two == 3, "arithmetic");
}

#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "Postcondition failed: done never")]
fn a_failed_postcondition_panics_with_debug_assertions() {
    let done = std::hint::black_box(false);
    openms::openms_postcondition!(done, "never");
}

#[cfg(not(debug_assertions))]
#[test]
fn failed_conditions_compile_away_without_debug_assertions() {
    let two = std::hint::black_box(2);
    openms::openms_precondition!(two == 3, "arithmetic");
    openms::openms_postcondition!(two != 2, "never");
}
