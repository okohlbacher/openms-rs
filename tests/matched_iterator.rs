// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `MatchedIterator_test.cpp` (core `bc9cc12`), one test per
//! `START_SECTION`, literals transcribed (tier 3), plus the traversal rules the
//! header documents (tier 4). The spectrum-level uses are covered by
//! `tests/comparison.rs` and `tests/comparison_scorers.rs`.

use openms::data_structures::matched_iterator::{
    DaTrait, MatchTrait, MatchedIterator, MatchedIteratorMut, MzCoordinate, PpmTrait, ValueTrait,
};
use openms::{MSSpectrum, Peak1D};

type Miv<'a> = MatchedIterator<'a, f64, ValueTrait>;

const REF: [f64; 9] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 13.0];
const TARGET: [f64; 5] = [-0.01, 2.5, 3.5, 7.0, 11.0];

fn assert_at(it: &Miv<'_>, reference: f64, target: f64, ref_idx: usize, tgt_idx: usize) {
    assert_eq!(it.reference_element(), Some(&reference));
    assert_eq!(it.target_element(), Some(&target));
    assert_eq!(it.ref_idx(), ref_idx);
    assert_eq!(it.tgt_idx(), tgt_idx);
}

#[test]
fn class_test_default_constructor_is_not_the_end() {
    // MatchedIterator() and `MIV it; TEST_EQUAL(it != it.end(), true)`.
    let mut it = Miv::default();
    assert!(it != Miv::end());
    assert_eq!(it.reference_element(), None);
    assert_eq!(it.target_element(), None);
    // "assignment is the only valid thing"
    let empty: [f64; 0] = [];
    it = Miv::new(&empty, &empty, 1.0);
    assert!(it == Miv::end());
}

#[test]
fn class_test_constructor_on_containers() {
    let empty: [f64; 0] = [];
    // empty reference container
    assert!(Miv::new(&empty, &TARGET, 0.001) == Miv::end());
    // empty target container
    assert!(Miv::new(&REF, &empty, 0.001) == Miv::end());

    // actual data, tolerance 0.001: only a single hit
    let mut mi = Miv::new(&REF, &TARGET, 0.001);
    assert_at(&mi, 7.0, 7.0, 7, 3);
    assert!(mi.advance().is_ok());
    assert!(mi == Miv::end());

    // actual data, tolerance 0.5
    let mut mi = Miv::new(&REF, &TARGET, 0.5);
    assert_at(&mi, 0.0, -0.01, 0, 0);
    assert!(mi.advance().is_ok());
    assert_at(&mi, 2.0, 2.5, 2, 1);
    // post-increment: yields the current pair, then advances
    let before = mi.clone();
    assert_eq!(mi.next(), Some((2, 1)));
    assert_at(&before, 2.0, 2.5, 2, 1);
    assert_at(&mi, 3.0, 2.5, 3, 1);
    assert_eq!(mi.next(), Some((3, 1)));
    assert_at(&mi, 4.0, 3.5, 4, 2);
    assert!(mi.advance().is_ok());
    assert_at(&mi, 7.0, 7.0, 7, 3);
    assert!(mi.advance().is_ok());
    assert!(mi == Miv::end());
    assert_eq!(mi.ref_idx(), 9); // points to ref.end()
    assert_eq!(mi.tgt_idx(), 4); // points to the last element of target
}

/// A local peak type with `getMZ()`, as the source's `PpmTrait` requires.
#[derive(Clone, Copy, Debug, PartialEq)]
struct MzPeak(f64);

impl MzCoordinate for MzPeak {
    fn mz_coordinate(&self) -> f64 {
        self.0
    }
}

#[test]
fn class_test_ppm_trait_on_a_spectrum() {
    let mzs = [90.0, 100.0, 200.0, 300.0, 400.0];
    let mut s = MSSpectrum::default();
    for mz in mzs {
        s.peaks.push(Peak1D::new(mz, 0.0));
    }
    let mut s2 = s.clone();
    // Add a constant to all peaks: Math::ppmToMass(2.5, 250.0).
    let offset = (2.5 / 1e6) * 250.0;
    for peak in &mut s2.peaks {
        peak.mz += offset;
    }
    let mut it = MatchedIterator::<Peak1D, PpmTrait>::new(&s.peaks, &s2.peaks, 2.5);
    // 90, 100 and 200 do not match: 2.5 ppm of them is below the offset.
    assert_eq!(it.ref_idx(), 3);
    assert_eq!(it.tgt_idx(), 3);
    assert!(it.advance().is_ok());
    assert_eq!(it.ref_idx(), 4);
    assert_eq!(it.tgt_idx(), 4);
    assert!(it.advance().is_ok());
    assert_eq!(it.ref_idx(), 5);
    assert!(it == MatchedIterator::end());
    assert_eq!(it.tgt_idx(), 4); // target is always valid

    // The same with a local element type through MzCoordinate.
    let local: Vec<MzPeak> = mzs.iter().map(|&mz| MzPeak(mz)).collect();
    let shifted: Vec<MzPeak> = mzs.iter().map(|&mz| MzPeak(mz + offset)).collect();
    let pairs: Vec<_> = MatchedIterator::<MzPeak, PpmTrait>::new(&local, &shifted, 2.5).collect();
    assert_eq!(pairs, [(3, 3), (4, 4)]);
}

#[test]
fn class_test_equality() {
    let mut mi = Miv::new(&REF, &TARGET, 0.5);
    assert_eq!(mi.reference_element(), Some(&0.0));
    assert_eq!(mi.target_element(), Some(&-0.01));
    assert_eq!(mi.ref_idx(), 0);
    assert!(mi.advance().is_ok());
    let mi2 = mi.clone();
    assert!(mi == mi2);
    assert_eq!(mi.reference_element(), mi2.reference_element());
    assert_eq!(mi.target_element(), mi2.target_element());
    assert_eq!(mi.ref_idx(), mi2.ref_idx());
    assert_at(&mi, 2.0, 2.5, 2, 1);
}

#[test]
fn class_test_inequality() {
    let mut mi = Miv::new(&REF, &TARGET, 0.5);
    let mi2 = mi.clone();
    assert!(!(mi != mi2));
    assert!(mi.advance().is_ok());
    assert!(!(mi == mi2));
    let mi3 = mi.clone();
    assert!(!(mi != mi3));
}

#[test]
fn equality_compares_the_containers_by_identity() {
    let copy = REF;
    let a = Miv::new(&REF, &TARGET, 0.5);
    let b = Miv::new(&copy, &TARGET, 0.5);
    assert_eq!(a.ref_idx(), b.ref_idx());
    assert_eq!(a.tgt_idx(), b.tgt_idx());
    assert!(a != b);
}

/// The source test's `PairTrait` over `std::pair<double, bool>`, a
/// user-defined trait whose distance is a `float`.
struct PairTrait;

impl MatchTrait<(f64, bool)> for PairTrait {
    fn allowed_tol(tolerance: f32, _reference: &(f64, bool)) -> f32 {
        tolerance
    }
    fn diff_absolute(reference: &(f64, bool), target: &(f64, bool)) -> f64 {
        f64::from((reference.0 - target.0).abs() as f32)
    }
}

#[test]
fn class_test_arrow_operator_with_a_user_trait() {
    let pairs = [(1.0, true), (2.0, false)];
    let mut it = MatchedIterator::<(f64, bool), PairTrait>::new(&pairs, &pairs, 0.1);
    assert!(it.advance().is_ok());
    assert_eq!(it.target_element().map(|pair| pair.0), Some(2.0));
}

#[test]
fn advancing_the_end_iterator_is_an_error() {
    let mut end = Miv::end();
    assert!(end.advance().is_err());
    assert!(end == Miv::end());
    let mut exhausted = Miv::new(&REF, &TARGET, 0.001);
    assert!(exhausted.advance().is_ok());
    assert!(exhausted.advance().is_err());
    assert_eq!(exhausted.next(), None);
}

#[test]
fn iterator_yields_every_pair_and_equal_distances_prefer_the_lower_target() {
    let pairs: Vec<_> = Miv::new(&REF, &TARGET, 0.5).collect();
    assert_eq!(pairs, [(0, 0), (2, 1), (3, 1), (4, 2), (7, 3)]);
    // 1.0 lies exactly between 0.5 and 1.5: the lower one is chosen.
    let reference = [1.0];
    let target = [0.5, 1.5];
    let pairs: Vec<_> = Miv::new(&reference, &target, 1.0).collect();
    assert_eq!(pairs, [(0, 0)]);
    // A duplicate target stops the forward search before a closer one.
    let reference = [3.0];
    let target = [1.0, 1.0, 2.9];
    let pairs: Vec<_> = Miv::new(&reference, &target, 5.0).collect();
    assert_eq!(pairs, [(0, 0)]);
}

#[test]
fn value_trait_on_f32_elements_and_da_trait() {
    let reference = [1.0f32, 2.0, 3.0];
    let target = [1.05f32, 2.9];
    let pairs: Vec<_> = MatchedIterator::<f32, ValueTrait>::new(&reference, &target, 0.1).collect();
    assert_eq!(pairs, [(0, 0), (2, 1)]);

    let peaks = [Peak1D::new(100.0, 1.0), Peak1D::new(200.0, 1.0)];
    let shifted = [Peak1D::new(100.004, 1.0), Peak1D::new(200.02, 1.0)];
    let pairs: Vec<_> = MatchedIterator::<Peak1D, DaTrait>::new(&peaks, &shifted, 0.01).collect();
    assert_eq!(pairs, [(0, 0)]);
}

#[test]
fn nan_distances_do_not_panic() {
    // A NaN distance at the first target counts as "getting worse"; the
    // source would decrement before the container start, here the cursor
    // stays at the first element.
    let reference = [f64::NAN, 1.0];
    let target = [0.9, 1.0];
    let pairs: Vec<_> = Miv::new(&reference, &target, 0.5).collect();
    assert_eq!(pairs, [(1, 1)]);
    // A NaN target blocks the forward search, so nothing matches.
    let target = [f64::NAN, 1.0];
    let pairs: Vec<_> = Miv::new(&reference, &target, 0.5).collect();
    assert!(pairs.is_empty());
}

#[test]
fn mutable_cursor_edits_the_matched_targets() {
    let reference = [1.0, 5.0];
    let mut target = [0.9, 3.0, 5.2];
    let mut cursor = MatchedIteratorMut::<f64, ValueTrait>::new(&reference, &mut target, 0.5);
    assert_eq!(cursor.reference_element(), Some(&1.0));
    assert_eq!(cursor.target_element(), Some(&0.9));
    assert_eq!((cursor.ref_idx(), cursor.tgt_idx()), (0, 0));
    if let Some(element) = cursor.target_element_mut() {
        *element = 1.0;
    }
    assert!(cursor.advance().is_ok());
    assert_eq!((cursor.ref_idx(), cursor.tgt_idx()), (1, 2));
    if let Some(element) = cursor.target_element_mut() {
        *element = 5.0;
    }
    assert!(cursor.advance().is_ok());
    assert!(cursor.is_end());
    assert_eq!(cursor.target_element_mut(), None);
    assert!(cursor.advance().is_err());
    assert_eq!(target, [1.0, 3.0, 5.0]);
}
