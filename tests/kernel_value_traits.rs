// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::{ChromatogramPeak, Peak1D};
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

fn hash(value: &impl Hash) -> u64 {
    let mut state = DefaultHasher::new();
    value.hash(&mut state);
    state.finish()
}

#[test]
fn source_labels_and_native_numeric_formatting() {
    assert_eq!(
        Peak1D::new(123.456, 7.25).to_string(),
        "POS: 123.456 INT: 7.25"
    );
    assert_eq!(
        ChromatogramPeak::new(12.75, -2.5).to_string(),
        "POS: 12.75 INT: -2.5"
    );
    assert_eq!(Peak1D::default().to_string(), "POS: 0 INT: 0");
    assert_eq!(ChromatogramPeak::default().to_string(), "POS: 0 INT: 0");
    assert_eq!(
        format!("{:.2}", Peak1D::new(123.456, 7.25)),
        "POS: 123.46 INT: 7.25"
    );
    assert_eq!(
        format!("{:.1}", ChromatogramPeak::new(12.75, -2.5)),
        "POS: 12.8 INT: -2.5"
    );
    assert_eq!(Peak1D::new(-0.0, -0.0).to_string(), "POS: -0 INT: -0");
    assert_eq!(
        ChromatogramPeak::new(f64::INFINITY, f32::NAN).to_string(),
        "POS: inf INT: NaN"
    );
}

#[test]
fn equal_signed_zeros_have_equal_hashes() {
    let peak = Peak1D::new(0.0, 0.0);
    let chromatogram_peak = ChromatogramPeak::new(0.0, 0.0);
    for position in [0.0, -0.0] {
        for intensity in [0.0, -0.0] {
            let other = Peak1D::new(position, intensity);
            assert_eq!(peak, other);
            assert_eq!(hash(&peak), hash(&other));
            let other = ChromatogramPeak::new(position, intensity);
            assert_eq!(chromatogram_peak, other);
            assert_eq!(hash(&chromatogram_peak), hash(&other));
        }
    }
}

#[test]
fn source_equal_and_differing_value_hash_cases() {
    // Source tests assert equal values hash equally and selected changed values
    // differ. This pins no numeric digest and makes no collision-free promise.
    let peak = Peak1D::new(100.5, 1000.0);
    assert_eq!(hash(&peak), hash(&Peak1D::new(100.5, 1000.0)));
    assert_ne!(hash(&peak), hash(&Peak1D::new(200.5, 1000.0)));
    assert_ne!(hash(&peak), hash(&Peak1D::new(100.5, 1001.0)));
    let peak = ChromatogramPeak::new(10.5, 1000.0);
    assert_eq!(hash(&peak), hash(&ChromatogramPeak::new(10.5, 1000.0)));
    assert_ne!(hash(&peak), hash(&ChromatogramPeak::new(20.5, 1000.0)));
    assert_ne!(hash(&peak), hash(&ChromatogramPeak::new(10.5, 1001.0)));
    // Equality intentionally stays partial. Hashing does not make NaN equal.
    let peak = Peak1D::new(f64::NAN, 1.0);
    let copy = peak;
    assert_ne!(peak, copy);
    assert_eq!(hash(&peak), hash(&peak));
    let peak = ChromatogramPeak::new(1.0, f32::NAN);
    let copy = peak;
    assert_ne!(peak, copy);
    assert_eq!(hash(&peak), hash(&peak));
}

#[test]
fn comparator_overloads_are_native_scalar_comparisons() {
    let a = Peak1D::new(1.0, 20.0);
    let b = Peak1D::new(2.0, 10.0);
    assert!(a.mz < b.mz);
    assert!(a.mz < 2.0);
    assert!(1.0 < b.mz);
    assert!(b.intensity < a.intensity);
    assert!(b.intensity < 20.0);
    assert!(10.0 < a.intensity);
    let a = ChromatogramPeak::new(1.0, 20.0);
    let b = ChromatogramPeak::new(2.0, 10.0);
    assert!(a.rt < b.rt);
    assert!(a.rt < 2.0);
    assert!(1.0 < b.rt);
    assert!(b.intensity < a.intensity);
    assert!(b.intensity < 20.0);
    assert!(10.0 < a.intensity);
}
