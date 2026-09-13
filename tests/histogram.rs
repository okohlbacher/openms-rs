// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Port of `Histogram_test.cpp` (core SDK bc9cc12): one test per
//! `START_SECTION`, cited by source line. The class test instantiates
//! `Histogram<float, float>` and mutates one shared object down the file; each
//! test here rebuilds the state its section sees. Transcribed literals are
//! tier 3.

use openms::Error;
use openms::math::histogram::{Histogram, MAX_BINS};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1.),
        "{actual:.17} != {expected:.17}"
    );
}

fn out_of_range(error: &Error) -> bool {
    matches!(error, Error::InvalidRange(_))
}

/// The shared `d` after `d.reset(4, 14, 2)` and the increments of L89.
fn filled() -> Histogram {
    let mut d = Histogram::with_bounds(4.0, 14.0, 2.0).unwrap();
    d.inc_by(4.0, 1.0).unwrap();
    d.inc_by(5.9, 1.0).unwrap();
    d.inc_by(8.0, 45.0).unwrap();
    d.inc_by(8.1, 1.0).unwrap();
    d.inc_by(9.9, 4.0).unwrap();
    d.inc_by(12.0, 1.0).unwrap();
    d.inc_by(13.1, 2.0).unwrap();
    d.inc_by(14.0, 3.0).unwrap();
    d
}

// L31 Histogram()
#[test]
fn section_default_constructor() {
    let d = Histogram::new();
    assert_eq!(d, Histogram::default());
    assert_eq!(d.min_bound(), 0.0);
    assert_eq!(d.max_bound(), 0.0);
    assert_eq!(d.bin_size(), 0.0);
    assert_eq!(d.len(), 0);
    assert!(d.is_empty());
    // Native: the source's minValue/maxValue dereference the end iterator of
    // the empty bin vector; this reports absence instead.
    assert_eq!(d.min_value(), None);
    assert_eq!(d.max_value(), None);
    assert!(out_of_range(&d.value_to_bin(0.0).unwrap_err()));
}

// L36 ~Histogram()
#[test]
fn section_destructor() {
    let d = Box::new(Histogram::with_bounds(0.0, 10.0, 1.0).unwrap());
    assert_eq!(d.len(), 10);
    drop(d);
}

// L42 Histogram(const Histogram& histogram)
#[test]
fn section_copy_constructor() {
    let d = Histogram::with_bounds(0.0, 10.0, 1.0).unwrap();
    let d2 = d.clone();
    assert_eq!(d, d2);
}

// L47 BinSizeType minBound() const
#[test]
fn section_min_bound() {
    let d = Histogram::with_bounds(0.0, 10.0, 1.0).unwrap();
    close(d.min_bound(), 0.0);
}

// L51 BinSizeType maxBound() const
#[test]
fn section_max_bound() {
    let d = Histogram::with_bounds(0.0, 10.0, 1.0).unwrap();
    close(d.max_bound(), 10.0);
}

// L55 BinSizeType binSize() const
#[test]
fn section_bin_size() {
    let d = Histogram::with_bounds(0.0, 10.0, 1.0).unwrap();
    close(d.bin_size(), 1.0);
}

// L59 Size size() const
#[test]
fn section_size() {
    let d = Histogram::with_bounds(0.0, 10.0, 1.0).unwrap();
    // ceil((10 - 0) / 1) == 10
    assert_eq!(d.len(), 10);
}

// L63 Histogram(BinSizeType min, BinSizeType max, BinSizeType bin_size)
#[test]
fn section_range_constructor() {
    let d3 = Histogram::with_bounds(5.5, 7.7, 0.2).unwrap();
    close(d3.min_bound(), 5.5);
    close(d3.max_bound(), 7.7);
    close(d3.bin_size(), 0.2);
    // Native: the source throws only on a non-positive bin width.
    assert!(matches!(
        Histogram::with_bounds(0.0, 1.0, 0.0).unwrap_err(),
        Error::InvalidValue(_)
    ));
    assert!(matches!(
        Histogram::with_bounds(0.0, 1.0, -1.0).unwrap_err(),
        Error::InvalidValue(_)
    ));
    assert!(matches!(
        Histogram::with_bounds(1.0, 0.0, 1.0).unwrap_err(),
        Error::InvalidValue(_)
    ));
    assert!(matches!(
        Histogram::with_bounds(0.0, f64::INFINITY, 1.0).unwrap_err(),
        Error::InvalidValue(_)
    ));
    // The bin count ceiling refuses before allocating.
    assert!(matches!(
        Histogram::with_bounds(0.0, 1.0, 1.0 / (MAX_BINS as f64 + 2.0)).unwrap_err(),
        Error::InvalidValue(_)
    ));
    // max == min is the source's explicit one-bin case.
    let single = Histogram::with_bounds(3.0, 3.0, 1.0).unwrap();
    assert_eq!(single.len(), 1);
    assert_eq!(single.value_to_bin(3.0).unwrap(), 0);
}

// L70 ValueType minValue() const
#[test]
fn section_min_value() {
    let d = Histogram::with_bounds(0.0, 10.0, 1.0).unwrap();
    assert_eq!(d.min_value(), Some(0.0));
    let filled = filled();
    assert_eq!(filled.min_value(), Some(0.0));
}

// L74 ValueType maxValue() const
#[test]
fn section_max_value() {
    let d = Histogram::with_bounds(0.0, 10.0, 1.0).unwrap();
    assert_eq!(d.max_value(), Some(0.0));
    let filled = filled();
    assert_eq!(filled.max_value(), Some(50.0));
}

// L78 ValueType operator[](Size index) const
#[test]
fn section_index() {
    let d = Histogram::with_bounds(4.0, 14.0, 2.0).unwrap();
    assert_eq!(d.len(), 5);
    for index in 0..5 {
        close(d.bin(index).unwrap(), 0.0);
    }
    // TEST_EXCEPTION(Exception::IndexOverflow, d[5])
    assert!(out_of_range(&d.bin(5).unwrap_err()));
    assert_eq!(d.bins(), &[0.0; 5]);
}

// L89 Size inc(BinSizeType val, ValueType increment = 1)
#[test]
fn section_inc() {
    let mut d = Histogram::with_bounds(4.0, 14.0, 2.0).unwrap();
    // Outside [min, max] on either side.
    assert!(out_of_range(&d.inc_by(3.9, 250.3).unwrap_err()));
    assert!(out_of_range(&d.inc_by(14.1, 250.3).unwrap_err()));

    assert_eq!(d.inc_by(4.0, 1.0).unwrap(), 0);
    assert_eq!(d.inc_by(5.9, 1.0).unwrap(), 0);
    assert_eq!(d.bins(), &[2.0, 0.0, 0.0, 0.0, 0.0]);

    assert_eq!(d.inc_by(8.0, 45.0).unwrap(), 2);
    assert_eq!(d.inc_by(8.1, 1.0).unwrap(), 2);
    assert_eq!(d.inc_by(9.9, 4.0).unwrap(), 2);
    assert_eq!(d.bins(), &[2.0, 0.0, 50.0, 0.0, 0.0]);

    assert_eq!(d.inc_by(12.0, 1.0).unwrap(), 4);
    assert_eq!(d.inc_by(13.1, 2.0).unwrap(), 4);
    // The upper bound is inclusive and always lands in the last bin.
    assert_eq!(d.inc_by(14.0, 3.0).unwrap(), 4);
    assert_eq!(d.bins(), &[2.0, 0.0, 50.0, 0.0, 6.0]);

    // The one-argument form increments by 1.
    let mut one = Histogram::with_bounds(0.0, 2.0, 1.0).unwrap();
    assert_eq!(one.inc(0.5).unwrap(), 0);
    assert_eq!(one.bins(), &[1.0, 0.0]);
    // Native: a non-finite increment or value is refused before any bin moves.
    assert!(matches!(
        one.inc_by(0.5, f64::NAN).unwrap_err(),
        Error::InvalidValue(_)
    ));
    assert!(out_of_range(&one.inc(f64::NAN).unwrap_err()));
    assert_eq!(one.bins(), &[1.0, 0.0]);
}

// L132 ConstIterator begin() const
#[test]
fn section_begin() {
    let d = filled();
    close(*d.iter().next().unwrap(), 2.0);
}

// L137 ConstIterator end() const
#[test]
fn section_end() {
    let d = filled();
    let mut it = d.iter();
    close(*it.next().unwrap(), 2.0);
    close(*it.next().unwrap(), 0.0);
    close(*it.next().unwrap(), 50.0);
    close(*it.next().unwrap(), 0.0);
    close(*it.next().unwrap(), 6.0);
    assert!(it.next().is_none());
    // The same sequence through `IntoIterator for &Histogram`.
    let mut collected: Vec<f64> = Vec::new();
    for value in &d {
        collected.push(*value);
    }
    assert_eq!(collected, vec![2.0, 0.0, 50.0, 0.0, 6.0]);
}

// L152 ValueType binValue(BinSizeType val) const
#[test]
fn section_bin_value() {
    let d = filled();
    assert!(out_of_range(&d.bin_value(3.9).unwrap_err()));
    close(d.bin_value(4.0).unwrap(), 2.0);
    close(d.bin_value(5.9).unwrap(), 2.0);
    close(d.bin_value(6.0).unwrap(), 0.0);
    close(d.bin_value(7.9).unwrap(), 0.0);
    close(d.bin_value(8.0).unwrap(), 50.0);
    close(d.bin_value(9.9).unwrap(), 50.0);
    close(d.bin_value(10.0).unwrap(), 0.0);
    close(d.bin_value(11.9).unwrap(), 0.0);
    close(d.bin_value(12.0).unwrap(), 6.0);
    close(d.bin_value(14.0).unwrap(), 6.0);
    assert!(out_of_range(&d.bin_value(14.1).unwrap_err()));
}

// L167 void reset(BinSizeType min, BinSizeType max, BinSizeType bin_size)
#[test]
fn section_reset() {
    let mut d = filled();
    d.reset(1.0, 11.0, 2.0).unwrap();
    close(d.min_bound(), 1.0);
    close(d.max_bound(), 11.0);
    assert_eq!(d.len(), 5);
    close(d.bin_size(), 2.0);
    assert_eq!(d.bins(), &[0.0; 5]);

    // Native: a refused reset leaves the histogram exactly as it was, where the
    // source clears its bins before throwing.
    let before = d.clone();
    assert!(d.reset(1.0, 11.0, 0.0).is_err());
    assert_eq!(d, before);
}

// L175 bool operator==(const Histogram& histogram) const
#[test]
fn section_equality() {
    let mut d = filled();
    d.reset(1.0, 11.0, 2.0).unwrap();
    let dist = Histogram::with_bounds(1.0, 11.0, 2.0).unwrap();
    assert_eq!(d, dist);
}

// L180 bool operator!=(const Histogram& histogram) const
#[test]
fn section_inequality() {
    let mut d = filled();
    d.reset(1.0, 11.0, 2.0).unwrap();
    // ceil((12 - 1) / 2) == 6 bins, against d's 5.
    let dist = Histogram::with_bounds(1.0, 12.0, 2.0).unwrap();
    assert_eq!(dist.len(), 6);
    assert_ne!(d, dist);
}

// L185 Histogram& operator=(const Histogram& histogram)
#[test]
fn section_assignment() {
    let mut d = filled();
    d.reset(1.0, 11.0, 2.0).unwrap();
    let mut dist = Histogram::new();
    assert_ne!(dist, d);
    dist = d.clone();
    assert_eq!(d, dist);
}

// L191 void applyLogTransformation(BinSizeType multiplier)
#[test]
fn section_apply_log_transformation() {
    let mut dist = Histogram::with_bounds(0.0, 5.0, 1.0).unwrap();
    dist.inc_by(0.5, 1.0).unwrap();
    dist.inc_by(1.5, 10.0).unwrap();
    dist.inc_by(2.5, 100.0).unwrap();
    dist.inc_by(3.5, 1000.0).unwrap();
    dist.inc_by(4.5, 10000.0).unwrap();
    dist.apply_log_transformation(1.0).unwrap();
    // TOLERANCE_ABSOLUTE(0.01); the values are ln(x + 1).
    // The first literal is ln(2); clippy flags it as an approximation of
    // `std::f64::consts::LN_2`, so it is spelled as that constant here.
    for (value, expected) in [
        (0.5, std::f64::consts::LN_2),
        (1.5, 2.3979),
        (2.5, 4.61512),
        (3.5, 6.90875),
        (4.5, 9.21044),
    ] {
        let actual = dist.bin_value(value).unwrap();
        assert!((actual - expected).abs() <= 0.01, "{actual} != {expected}");
    }
    // Derived: exactly ln(10001) for the last bin.
    close(dist.bin_value(4.5).unwrap(), 10001.0f64.ln());

    // Native: a bin at or below -1 has no logarithm; the source stores NaN.
    let mut negative = Histogram::with_bounds(0.0, 2.0, 1.0).unwrap();
    negative.inc_by(0.5, -2.0).unwrap();
    let before = negative.clone();
    assert!(matches!(
        negative.apply_log_transformation(1.0).unwrap_err(),
        Error::InvalidValue(_)
    ));
    assert_eq!(negative, before);
}

// L207 BinSizeType centerOfBin(Size bin_index) const
#[test]
fn section_center_of_bin() {
    let dist = Histogram::with_bounds(0.0, 5.0, 1.0).unwrap();
    close(dist.center_of_bin(0).unwrap(), 0.5);
    close(dist.center_of_bin(1).unwrap(), 1.5);
    close(dist.center_of_bin(2).unwrap(), 2.5);
    close(dist.center_of_bin(3).unwrap(), 3.5);
    close(dist.center_of_bin(4).unwrap(), 4.5);
    assert!(out_of_range(&dist.center_of_bin(5).unwrap_err()));
}

// L222 BinSizeType leftBorderOfBin(Size bin_index) const
#[test]
fn section_left_border_of_bin() {
    let dist = Histogram::with_bounds(0.0, 5.0, 1.0).unwrap();
    for index in 0..5 {
        assert_eq!(dist.left_border_of_bin(index).unwrap(), index as f64);
    }
    assert!(out_of_range(&dist.left_border_of_bin(5).unwrap_err()));
}

// L233 BinSizeType rightBorderOfBin(Size bin_index) const
#[test]
fn section_right_border_of_bin() {
    let dist = Histogram::with_bounds(0.0, 5.0, 1.0).unwrap();
    assert_eq!(dist.right_border_of_bin(0).unwrap(), 1.0);
    assert_eq!(dist.right_border_of_bin(1).unwrap(), 2.0);
    assert_eq!(dist.right_border_of_bin(2).unwrap(), 3.0);
    assert_eq!(dist.right_border_of_bin(3).unwrap(), 4.0);
    // The last bin is special: it holds the inclusive upper bound, so its open
    // right border is the next representable value above it. The class test
    // instantiates `Histogram<float, float>` and writes
    // `std::nextafter(5.0f, 6.0f)`; in `f64` the step is the f64 ulp.
    let next_above_five = f64::from_bits(5.0f64.to_bits() + 1);
    assert_eq!(dist.right_border_of_bin(4).unwrap(), next_above_five);
    assert!(next_above_five > 5.0);
    assert!(out_of_range(&dist.right_border_of_bin(5).unwrap_err()));
}

// Native: members with no class-test section of their own.
#[test]
fn cumulative_increments_have_no_section() {
    // incUntil: every bin below the value, optionally including it.
    let mut d = Histogram::with_bounds(0.0, 5.0, 1.0).unwrap();
    assert_eq!(d.inc_until(3.5, false, 1.0).unwrap(), 3);
    assert_eq!(d.bins(), &[1.0, 1.0, 1.0, 0.0, 0.0]);
    assert_eq!(d.inc_until(3.5, true, 1.0).unwrap(), 3);
    assert_eq!(d.bins(), &[2.0, 2.0, 2.0, 1.0, 0.0]);

    // incFrom: every bin above the value, optionally including it.
    let mut d = Histogram::with_bounds(0.0, 5.0, 1.0).unwrap();
    assert_eq!(d.inc_from(1.5, false, 1.0).unwrap(), 1);
    assert_eq!(d.bins(), &[0.0, 0.0, 1.0, 1.0, 1.0]);
    assert_eq!(d.inc_from(1.5, true, 1.0).unwrap(), 1);
    assert_eq!(d.bins(), &[0.0, 1.0, 2.0, 2.0, 2.0]);

    // getCumulativeHistogram, both directions.
    let mut up = Histogram::with_bounds(0.0, 5.0, 1.0).unwrap();
    up.add_cumulative(&[0.5, 2.5], false, true).unwrap();
    assert_eq!(up.bins(), &[1.0, 1.0, 2.0, 2.0, 2.0]);
    let mut down = Histogram::with_bounds(0.0, 5.0, 1.0).unwrap();
    down.add_cumulative(&[0.5, 2.5], true, true).unwrap();
    assert_eq!(down.bins(), &[2.0, 1.0, 1.0, 0.0, 0.0]);
    assert!(out_of_range(
        &down.add_cumulative(&[9.0], true, true).unwrap_err()
    ));

    // The data-iterator constructor.
    let filled = Histogram::from_values(&[0.5, 0.6, 4.9], 0.0, 5.0, 1.0).unwrap();
    assert_eq!(filled.bins(), &[2.0, 0.0, 0.0, 0.0, 1.0]);
    assert!(out_of_range(
        &Histogram::from_values(&[9.0], 0.0, 5.0, 1.0).unwrap_err()
    ));
}

// Native: the stream operator and the bin lookup's own guards.
#[test]
fn display_and_value_to_bin_guards() {
    let mut d = Histogram::with_bounds(0.0, 2.0, 1.0).unwrap();
    d.inc_by(0.5, 3.0).unwrap();
    assert_eq!(d.to_string(), "0.5\t3\n1.5\t0\n");

    assert_eq!(d.value_to_bin(0.0).unwrap(), 0);
    assert_eq!(d.value_to_bin(1.0).unwrap(), 1);
    // The inclusive upper bound belongs to the last bin, not to a bin the
    // arithmetic would place past the end.
    assert_eq!(d.value_to_bin(2.0).unwrap(), 1);
    assert!(out_of_range(&d.value_to_bin(-0.001).unwrap_err()));
    assert!(out_of_range(&d.value_to_bin(2.001).unwrap_err()));
    // Native: both source comparisons are false for a NaN, which then reaches
    // an unsigned conversion of floor(NaN).
    assert!(out_of_range(&d.value_to_bin(f64::NAN).unwrap_err()));
}
