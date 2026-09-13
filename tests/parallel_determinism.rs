// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The determinism contract for [`openms::concept::parallel`].
//!
//! A parallel result must be **bit-identical** to the serial one. A scientific
//! value that moves with the thread count is not a value, and it is exactly
//! what OpenMP's `reduction(+:...)` gives up: floating-point addition is not
//! associative, so an OpenMP sum depends on how the runtime split the loop.
//!
//! These tests are the enforcement, not the documentation. They run the same
//! inputs at several thread counts and compare raw bits via `to_bits`, because
//! `==` would accept a value that drifted within a tolerance and drift is the
//! defect being excluded.

use openms::concept::parallel::{Threads, map_collect, sum_in_order};

/// Values spanning several magnitudes, so a changed summation order shows up in
/// the mantissa rather than being masked by similar-sized terms.
fn spread(n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| {
            let x = i as f64;
            match i % 4 {
                0 => x * 1.000_000_1,
                1 => x * 1e-7,
                2 => x * 1e7,
                _ => -x * 0.333_333_333_333,
            }
        })
        .collect()
}

#[test]
fn map_collect_preserves_input_order_at_every_thread_count() {
    let items: Vec<usize> = (0..5_000).collect();
    let serial = map_collect(&items, Threads::serial(), |i| i * 3 + 1);
    assert_eq!(serial.len(), items.len());
    // Order is a property of the call, not of which worker finished first.
    for &count in &[1_i64, 2, 3, 4, 8, 16, 0] {
        let parallel = map_collect(&items, Threads::from_cli(count), |i| i * 3 + 1);
        assert_eq!(
            parallel, serial,
            "thread count {count} reordered the results"
        );
    }
}

#[test]
fn sums_are_bit_identical_across_thread_counts() {
    for &n in &[0_usize, 1, 2, 63, 64, 65, 1_000, 10_000] {
        let values = spread(n);
        let baseline = sum_in_order(&values, Threads::serial(), 64);
        for &count in &[1_i64, 2, 3, 4, 8, 16, 0] {
            let got = sum_in_order(&values, Threads::from_cli(count), 64);
            assert_eq!(
                got.to_bits(),
                baseline.to_bits(),
                "n={n}, threads={count}: {got:.17e} vs {baseline:.17e}"
            );
        }
    }
}

/// The chunk size is part of the contract, not a tuning knob: it fixes the
/// summation tree. This pins that it is load-bearing, so nobody "optimises" it
/// later without realising the last bits move.
#[test]
fn the_chunk_size_is_part_of_the_contract() {
    let values = spread(4_096);
    let a = sum_in_order(&values, Threads::from_cli(8), 64);
    let b = sum_in_order(&values, Threads::from_cli(8), 512);
    // Both are correct sums; they need not agree bit-for-bit, and a caller
    // matching a published number must keep the chunk fixed.
    assert!((a - b).abs() < 1e-6 * a.abs().max(1.0), "{a} vs {b}");
    // ...but for one fixed chunk the answer never moves.
    assert_eq!(
        sum_in_order(&values, Threads::from_cli(2), 64).to_bits(),
        a.to_bits()
    );
}

#[test]
fn the_cli_thread_parameter_is_honoured_and_total() {
    assert_eq!(Threads::from_cli(1).get(), 1);
    assert_eq!(Threads::from_cli(4).get(), 4);
    // Source `-threads 0` means every available core.
    assert_eq!(Threads::from_cli(0).get(), Threads::all().get());
    // A nonsensical count runs single-threaded rather than failing a tool that
    // is otherwise ready to run, matching how the source treats it.
    assert_eq!(Threads::from_cli(-1).get(), 1);
    assert_eq!(Threads::from_cli(i64::MIN).get(), 1);
    assert!(Threads::all().get() >= 1);
    assert_eq!(Threads::default().get(), Threads::all().get());
}

/// An empty or single-element input must not take a different path.
#[test]
fn degenerate_inputs_agree_with_the_serial_path() {
    let empty: Vec<f64> = Vec::new();
    assert_eq!(sum_in_order(&empty, Threads::from_cli(8), 64), 0.0);
    assert_eq!(sum_in_order(&[42.5], Threads::from_cli(8), 64), 42.5);
    let none: Vec<u8> = Vec::new();
    assert!(map_collect(&none, Threads::from_cli(8), |b| *b).is_empty());
    // A zero chunk must not divide by zero; it is clamped to one.
    assert_eq!(sum_in_order(&[1.0, 2.0], Threads::from_cli(4), 0), 3.0);
}
