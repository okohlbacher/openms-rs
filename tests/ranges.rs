// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `RangeManager_test.cpp` (38 sections) and
//! `SpectrumRangeManager_test.cpp` (10 sections) at core SDK `bc9cc12`, plus
//! native container and boundary tests. Expected literals are transcribed from
//! the class tests (tier 3); no C++ is executed.

use openms::Error;
use openms::kernel::ranges::{HasRangeType, MSDim, RangeBase, RangeManager, SpectrumRangeManager};
use openms::kernel::{DataArray, MobilityPeak1D, Mobilogram};
use openms::{ChromatogramPeak, MSChromatogram, MSExperiment, MSSpectrum, Peak1D};

fn rb(min: f64, max: f64) -> RangeBase {
    RangeBase::from_min_max(min, max).unwrap()
}

/// The class test's `RM::updateRanges()`: three Peak2D values
/// (RT 2 / 100 / 2, m/z 500 / 1300 / 500, intensity 1 / 47110 / 1) folded into
/// `RangeManagerContainer<RangeRT, RangeMZ, RangeIntensity, RangeMobility>`.
fn rm_update_ranges() -> RangeManager {
    let mut rm = RangeManager::experiment();
    rm.clear_ranges();
    for (rt, mz, intensity) in [
        (2.0, 500.0, 1.0f32),
        (100.0, 1300.0, 47110.0),
        (2.0, 500.0, 1.0),
    ] {
        rm.extend_rt(rt).unwrap();
        rm.extend_mz(mz).unwrap();
        rm.extend_intensity(f64::from(intensity)).unwrap();
    }
    rm
}

/// The class test's `RM::updateRanges2()`: one point only.
fn rm_update_ranges2(rm: &mut RangeManager) {
    rm.clear_ranges();
    rm.extend_rt(2.0).unwrap();
    rm.extend_mz(500.0).unwrap();
    rm.extend_intensity(1.0).unwrap();
}

fn assert_rm_reference(rm: &RangeManager) {
    assert_eq!(rm.min_rt().unwrap(), 2.0);
    assert_eq!(rm.min_mz().unwrap(), 500.0);
    assert_eq!(rm.max_rt().unwrap(), 100.0);
    assert_eq!(rm.max_mz().unwrap(), 1300.0);
    assert_eq!(rm.min_intensity().unwrap(), 1.0);
    assert_eq!(rm.max_intensity().unwrap(), 47110.0);
}

// ---------------------------------------------------------------------------
// RangeManager_test.cpp — RangeBase sections
// ---------------------------------------------------------------------------

#[test]
fn range_base_default_is_empty() {
    // START_SECTION(RangeBase())
    let b = RangeBase::new();
    assert!(b.is_empty());
    assert_eq!(b, RangeBase::default());
}

#[test]
fn range_base_min_max_constructor() {
    // START_SECTION(RangeBase(const double min, const double max))
    let b = rb(4.0, 6.0);
    assert!(!b.is_empty());
    assert_eq!(b.min().unwrap(), 4.0);
    assert_eq!(b.max().unwrap(), 6.0);
    assert!(matches!(
        RangeBase::from_min_max(6.0, 3.0),
        Err(Error::InvalidRange(_))
    ));
}

#[test]
fn range_base_copy_constructor() {
    // START_SECTION(const RangeBase& rhs)
    let b_ = rb(4.0, 6.0);
    let b = b_;
    assert!(!b.is_empty());
    assert_eq!(b.min().unwrap(), 4.0);
    assert_eq!(b.max().unwrap(), 6.0);
}

#[test]
fn range_base_assignment() {
    // START_SECTION(RangeBase& operator=(const RangeBase& rhs))
    let b_ = rb(4.0, 6.0);
    let mut b = RangeBase::new();
    assert!(b.is_empty());
    b = b_;
    assert!(!b.is_empty());
    assert_eq!(b.min().unwrap(), 4.0);
    assert_eq!(b.max().unwrap(), 6.0);
}

#[test]
fn range_base_clear() {
    // START_SECTION(void clear())
    let mut b = rb(4.0, 6.0);
    assert!(!b.is_empty());
    b.clear();
    assert!(b.is_empty());
}

#[test]
fn range_base_is_empty() {
    // START_SECTION(bool isEmpty() const) — NOT_TESTABLE in source (tested above).
    assert!(RangeBase::new().is_empty());
    assert!(!rb(1.0, 1.0).is_empty());
    assert!(matches!(
        RangeBase::new().min(),
        Err(Error::InvalidRange(_))
    ));
    assert!(matches!(
        RangeBase::new().max(),
        Err(Error::InvalidRange(_))
    ));
}

#[test]
fn range_base_contains_value() {
    // START_SECTION(bool contains(const double value) const)
    let b = rb(4.0, 6.0);
    assert!(b.contains(5.0));
    assert!(!b.contains(3.0));
    assert!(!b.contains(7.0));
    let empty = RangeBase::new();
    assert!(!empty.contains(5.0));
}

#[test]
fn range_base_contains_range() {
    // START_SECTION(bool contains(const RangeBase& inner_range) const)
    let b = rb(2.0, 6.0);
    let (inner1, inner2, inner3) = (rb(2.0, 4.0), rb(3.0, 4.0), rb(4.0, 6.0));
    let (over1, over2, outer) = (rb(1.0, 4.0), rb(3.0, 7.0), rb(1.0, 7.0));
    assert!(b.contains_range(&inner1));
    assert!(b.contains_range(&inner2));
    assert!(b.contains_range(&inner3));
    assert!(!b.contains_range(&over1));
    assert!(!b.contains_range(&over2));
    assert!(!b.contains_range(&outer));
    assert!(outer.contains_range(&b));
}

#[test]
fn range_base_set_min() {
    // START_SECTION(void setMin(const double min))
    let mut b = rb(4.0, 6.0);
    b.set_min(5.0).unwrap();
    assert_eq!(b.min().unwrap(), 5.0);
    b.set_min(7.0).unwrap(); // also increases max
    assert_eq!(b.min().unwrap(), 7.0);
    assert_eq!(b.max().unwrap(), 7.0);
}

#[test]
fn range_base_set_max() {
    // START_SECTION(void setMax(const double max))
    let mut b = rb(4.0, 6.0);
    b.set_max(5.0).unwrap();
    assert_eq!(b.max().unwrap(), 5.0);
    b.set_max(2.0).unwrap(); // also decreases min
    assert_eq!(b.min().unwrap(), 2.0);
    assert_eq!(b.max().unwrap(), 2.0);
}

#[test]
fn range_base_get_min() {
    // START_SECTION(double getMin() const) — NOT_TESTABLE in source.
    assert_eq!(rb(4.0, 6.0).min().unwrap(), 4.0);
}

#[test]
fn range_base_get_max() {
    // START_SECTION(double getMax() const) — NOT_TESTABLE in source.
    assert_eq!(rb(4.0, 6.0).max().unwrap(), 6.0);
}

#[test]
fn range_base_extend_range() {
    // START_SECTION(void extend(const RangeBase& other))
    let mut b = rb(4.0, 6.0);
    let other = rb(1.0, 8.0);
    b.extend(&other);
    assert_eq!(b.min().unwrap(), 1.0);
    assert_eq!(b.max().unwrap(), 8.0);
}

#[test]
fn range_base_extend_value() {
    // START_SECTION(void extend(const double value))
    let mut b = rb(4.0, 6.0);
    b.extend_value(1.0).unwrap();
    assert_eq!(b.min().unwrap(), 1.0);
    assert_eq!(b.max().unwrap(), 6.0);
    let mut b2 = rb(4.0, 6.0);
    b2.extend_value(8.0).unwrap();
    assert_eq!(b2.min().unwrap(), 4.0);
    assert_eq!(b2.max().unwrap(), 8.0);
    let mut b3 = rb(4.0, 6.0);
    b3.extend_value(5.0).unwrap();
    assert_eq!(b3.min().unwrap(), 4.0);
    assert_eq!(b3.max().unwrap(), 6.0);
}

#[test]
fn range_base_extend_left_right() {
    // START_SECTION(void extendLeftRight(const double by))
    let mut b = rb(4.0, 6.0);
    b.extend_left_right(1.0).unwrap();
    assert_eq!(b.min().unwrap(), 3.0);
    assert_eq!(b.max().unwrap(), 7.0);
    let mut b2 = rb(2.0, 8.0);
    b2.extend_left_right(-2.0).unwrap();
    assert_eq!(b2.min().unwrap(), 4.0);
    assert_eq!(b2.max().unwrap(), 6.0);
    b2.extend_left_right(-19.0).unwrap();
    assert!(b2.is_empty());
    let mut empty = RangeBase::new();
    empty.extend_left_right(100.0).unwrap();
    assert!(empty.is_empty());
}

#[test]
fn range_base_clamp_to() {
    // START_SECTION(void clampTo(const RangeBase& other))
    let mut b = rb(-4.0, 6.0);
    b.clamp_to(&rb(-2.0, 7.0)).unwrap();
    assert_eq!(b.min().unwrap(), -2.0);
    assert_eq!(b.max().unwrap(), 6.0);
    let mut b2 = rb(4.0, 6.0);
    b2.clamp_to(&rb(1.0, 5.0)).unwrap();
    assert_eq!(b2.min().unwrap(), 4.0);
    assert_eq!(b2.max().unwrap(), 5.0);
    b2.clamp_to(&rb(4.5, 4.5)).unwrap();
    assert_eq!(b2.min().unwrap(), 4.5);
    assert_eq!(b2.max().unwrap(), 4.5);
    let mut b3 = rb(4.0, 6.0);
    b3.clamp_to(&rb(10.0, 11.0)).unwrap();
    assert!(b3.is_empty());
    let mut b4 = rb(4.0, 6.0);
    assert!(matches!(
        b4.clamp_to(&RangeBase::new()),
        Err(Error::InvalidRange(_))
    ));
    assert_eq!(b4, rb(4.0, 6.0));
}

#[test]
fn range_base_push_into() {
    // START_SECTION(void pushInto(const RangeBase& sandbox))
    let mut b = rb(-4.0, 6.0);
    b.push_into(&rb(-2.0, 7.0)).unwrap(); // moves and clips
    assert_eq!(b.min().unwrap(), -2.0);
    assert_eq!(b.max().unwrap(), 7.0);
    let mut b2 = rb(4.0, 6.0);
    b2.push_into(&rb(1.0, 15.0)).unwrap(); // does nothing
    assert_eq!(b2.min().unwrap(), 4.0);
    assert_eq!(b2.max().unwrap(), 6.0);
    b2.push_into(&rb(4.5, 4.5)).unwrap(); // hard clip inside old range
    assert_eq!(b2.min().unwrap(), 4.5);
    assert_eq!(b2.max().unwrap(), 4.5);
    let mut b3 = rb(4.0, 6.0); // move left, no clip
    b3.push_into(&rb(-10.0, 5.0)).unwrap();
    assert_eq!(b3.min().unwrap(), 3.0);
    assert_eq!(b3.max().unwrap(), 5.0);
    b3.push_into(&rb(4.0, 10.0)).unwrap(); // move right, no clip
    assert_eq!(b3.min().unwrap(), 4.0);
    assert_eq!(b3.max().unwrap(), 6.0);
    let mut b4 = rb(4.0, 6.0); // hard clip outside old range
    b4.push_into(&rb(-10.0, -10.0)).unwrap();
    assert_eq!(b4.min().unwrap(), -10.0);
    assert_eq!(b4.max().unwrap(), -10.0);
    let mut b5 = rb(4.0, 6.0);
    assert!(matches!(
        b5.push_into(&RangeBase::new()),
        Err(Error::InvalidRange(_))
    ));
    assert_eq!(b5, rb(4.0, 6.0));
}

#[test]
fn range_base_scale_by() {
    // START_SECTION(void scaleBy(const double factor))
    let mut b = rb(4.0, 6.0);
    b.scale_by(10.0).unwrap(); // diff is 2, so extend distance to 20, by increase of 9 on each side
    assert_eq!(b.min().unwrap(), 4.0 - 9.0);
    assert_eq!(b.max().unwrap(), 6.0 + 9.0);
    // scaling empty ranges does nothing
    let (mut empty1, empty2) = (RangeBase::new(), RangeBase::new());
    empty1.scale_by(10.0).unwrap();
    assert_eq!(empty1, empty2);
}

#[test]
fn range_base_shift() {
    // START_SECTION(void shift(const double distance))
    let mut b = rb(4.0, 6.0);
    b.shift(10.0).unwrap();
    assert_eq!(b.min().unwrap(), 14.0);
    assert_eq!(b.max().unwrap(), 16.0);
    // shifting empty ranges does nothing
    let (mut empty1, empty2) = (RangeBase::new(), RangeBase::new());
    empty1.shift(10.0).unwrap();
    assert_eq!(empty1, empty2);
}

#[test]
fn range_base_center() {
    // START_SECTION(double center() const) — source returns NaN when empty; the
    // port returns None.
    assert_eq!(rb(4.0, 6.0).center(), Some(5.0));
    assert_eq!(RangeBase::new().center(), None);
}

#[test]
fn range_base_span() {
    // START_SECTION(double getSpan() const) — source returns NaN when empty; the
    // port returns None.
    assert_eq!(rb(4.0, 6.0).span(), Some(2.0));
    assert_eq!(RangeBase::new().span(), None);
}

#[test]
fn range_base_equality() {
    // START_SECTION(bool operator==(const RangeBase& rhs) const)
    let (b, b2, empty1, empty2) = (
        rb(4.0, 6.0),
        rb(4.0, 6.0),
        RangeBase::new(),
        RangeBase::new(),
    );
    assert_ne!(b, empty1);
    assert_eq!(b, b2);
    assert_eq!(empty1, empty2);
}

// ---------------------------------------------------------------------------
// RangeManager_test.cpp — RangeManager sections
// ---------------------------------------------------------------------------

#[test]
fn range_manager_constructor() {
    // START_SECTION((RangeMType())) — `ptr != nullptr`; the port checks the
    // constructed value instead.
    let rm = RangeManager::experiment();
    assert_eq!(rm.has_range(), HasRangeType::None);
    assert_eq!(rm.dims().collect::<Vec<_>>(), MSDim::ALL.to_vec());
}

#[test]
fn range_manager_destructor() {
    // START_SECTION((~RangeMType())) — `delete ptr`; dropping a value is the port.
    let rm = RangeManager::experiment();
    let _ = rm;
}

#[test]
fn range_manager_copy_constructor() {
    // START_SECTION((RangeManager(const RangeManager& rhs)))
    let rm0 = rm_update_ranges();
    let rm = rm0;
    assert_rm_reference(&rm);
}

#[test]
fn range_manager_assignment() {
    // START_SECTION((RangeManager& operator=(const RangeManager& rhs)))
    let rm0 = rm_update_ranges();
    let mut rm = RangeManager::experiment();
    assert_eq!(rm.has_range(), HasRangeType::None);
    rm = rm0;
    assert_rm_reference(&rm);
}

#[test]
fn range_manager_equality() {
    // START_SECTION((bool operator==(const RangeManager& rhs) const))
    let (mut rm0, rm) = (RangeManager::experiment(), RangeManager::experiment());
    assert!(rm == rm0);
    rm0 = rm_update_ranges();
    assert!(!(rm == rm0));
}

#[test]
fn range_manager_inequality() {
    // START_SECTION((bool operator!=(const RangeManager& rhs) const))
    let (mut rm0, rm) = (RangeManager::experiment(), RangeManager::experiment());
    assert!(!(rm != rm0));
    rm0 = rm_update_ranges();
    assert!(rm != rm0);
}

#[test]
fn range_manager_update_ranges() {
    // START_SECTION((virtual void updateRanges()=0))
    let mut rm = rm_update_ranges();
    assert_rm_reference(&rm);
    rm = rm_update_ranges(); // second time to check the initialization
    assert_eq!(rm.min_rt().unwrap(), 2.0);
    assert_eq!(rm.min_mz().unwrap(), 500.0);
    assert!(!rm.is_dim_empty(MSDim::Rt).unwrap());
    assert_eq!(rm.max_rt().unwrap(), 100.0);
    assert_eq!(rm.max_mz().unwrap(), 1300.0);
    assert!(!rm.is_dim_empty(MSDim::Mz).unwrap());
    assert_eq!(rm.min_intensity().unwrap(), 1.0);
    assert_eq!(rm.max_intensity().unwrap(), 47110.0);
    assert!(!rm.is_dim_empty(MSDim::Intensity).unwrap());
    assert!(rm.is_dim_empty(MSDim::Mobility).unwrap());

    // test with only one point
    rm_update_ranges2(&mut rm);
    assert_eq!(rm.min_rt().unwrap(), 2.0);
    assert_eq!(rm.min_mz().unwrap(), 500.0);
    assert_eq!(rm.max_rt().unwrap(), 2.0);
    assert_eq!(rm.max_mz().unwrap(), 500.0);
    assert_eq!(rm.min_intensity().unwrap(), 1.0);
    assert_eq!(rm.max_intensity().unwrap(), 1.0);
}

#[test]
fn range_manager_has_range() {
    // START_SECTION(HasRangeType hasRange() const)
    let mut rm = RangeManager::experiment();
    assert_eq!(rm.has_range(), HasRangeType::None);
    rm = rm_update_ranges();
    assert_eq!(rm.has_range(), HasRangeType::Some);
    rm.extend_mobility(56.4).unwrap();
    assert_eq!(rm.has_range(), HasRangeType::All);
}

#[test]
fn range_manager_contains_all() {
    // START_SECTION(template<typename... RangeBasesOther>
    //               bool containsAll(const RangeManager<RangeBasesOther...>& rhs) const)
    let mut rm = rm_update_ranges();
    let mut outer = rm;
    assert!(rm.contains_all(&outer).unwrap());
    assert!(outer.contains_all(&rm).unwrap());
    outer.scale_by(1.1).unwrap();
    assert!(!rm.contains_all(&outer).unwrap());
    assert!(outer.contains_all(&rm).unwrap());
    outer.scale_by(0.5).unwrap();
    assert!(rm.contains_all(&outer).unwrap());
    assert!(!outer.contains_all(&rm).unwrap());

    outer = rm;
    // empty dimensions in the rhs are considered contained
    outer.extend_mobility(56.4).unwrap(); // rm.mobility is empty
    assert!(!rm.contains_all(&outer).unwrap());
    assert!(outer.contains_all(&rm).unwrap());
    // empty dimensions do not count
    outer
        .range_for_dim_mut(MSDim::Mz)
        .unwrap()
        .scale_by(0.5)
        .unwrap(); // mz range is smaller
    rm.clear_dim(MSDim::Mz); // but now does not count anymore
    assert!(!rm.contains_all(&outer).unwrap()); // due to mobility from above
    assert!(outer.contains_all(&rm).unwrap());

    // no ranges overlap...
    let rmz = RangeManager::new(&[MSDim::Rt, MSDim::Mz]).unwrap();
    let im = RangeManager::new(&[MSDim::Intensity, MSDim::Mobility]).unwrap();
    assert!(matches!(rmz.contains_all(&im), Err(Error::InvalidRange(_))));
}

#[test]
fn range_manager_extend() {
    // START_SECTION(template<typename... RangeBasesOther>
    //               void extend(const RangeManager<RangeBasesOther...>& rhs))
    let rm = rm_update_ranges();
    let mut mid = RangeManager::new(&[MSDim::Mz, MSDim::Intensity]).unwrap();
    mid.assign(&rm).unwrap(); // assigns only overlapping dimensions
    assert_eq!(mid.min_mz().unwrap(), 500.0);
    assert_eq!(mid.max_mz().unwrap(), 1300.0);
    assert_eq!(mid.min_intensity().unwrap(), 1.0);
    assert_eq!(mid.max_intensity().unwrap(), 47110.0);

    let mut small = RangeManager::new(&[MSDim::Intensity]).unwrap();
    small.extend_intensity(123456.7).unwrap();
    mid.extend(&small).unwrap();
    assert_eq!(mid.min_mz().unwrap(), 500.0);
    assert_eq!(mid.max_mz().unwrap(), 1300.0);
    assert_eq!(mid.min_intensity().unwrap(), 1.0);
    assert_eq!(mid.max_intensity().unwrap(), 123456.7);
}

#[test]
fn range_manager_scale_by() {
    // START_SECTION(void scaleBy(const double factor))
    let mut rm = rm_update_ranges();
    rm.scale_by(2.0).unwrap();
    assert_eq!(rm.min_rt().unwrap(), 2.0 - 49.0);
    assert_eq!(rm.max_rt().unwrap(), 100.0 + 49.0);
    assert_eq!(rm.min_mz().unwrap(), 500.0 - 400.0);
    assert_eq!(rm.max_mz().unwrap(), 1300.0 + 400.0);
    assert_eq!(rm.min_intensity().unwrap(), 1.0 - (47109.0 / 2.0));
    assert_eq!(rm.max_intensity().unwrap(), 47110.0 + (47109.0 / 2.0));
    assert!(rm.is_dim_empty(MSDim::Mobility).unwrap());

    // scaling a dimension where min == max does nothing
    let mut rtmz = RangeManager::new(&[MSDim::Rt, MSDim::Mz]).unwrap();
    rtmz.extend_mz(100.0).unwrap();
    rtmz.extend_rt(50.0).unwrap();
    let copy = rtmz;
    rtmz.scale_by(2.0).unwrap();
    assert_eq!(rtmz, copy);

    // scaling empty dimensions does nothing
    let (mut rm_empty, rm_empty2) = (RangeManager::experiment(), RangeManager::experiment());
    rm_empty.scale_by(4.0).unwrap();
    assert_eq!(rm_empty, rm_empty2);
}

#[test]
fn range_manager_push_into() {
    // START_SECTION(template<typename... RangeBasesOther> void pushInto(const RangeManager<RangeBasesOther...>& sandbox))
    let mut rm = rm_update_ranges();
    let mut rmi = RangeManager::new(&[MSDim::Mz, MSDim::Intensity]).unwrap();
    rmi.extend_mz(700.0).unwrap(); // shift
    rmi.extend_mz(2000.0).unwrap();
    rmi.extend_intensity(500.0).unwrap(); // shift and clamp
    rmi.extend_intensity(600.0).unwrap();
    rm.push_into(&rmi).unwrap();
    assert_eq!(rm.min_rt().unwrap(), 2.0);
    assert_eq!(rm.max_rt().unwrap(), 100.0);
    assert_eq!(rm.min_mz().unwrap(), 500.0 + 200.0);
    assert_eq!(rm.max_mz().unwrap(), 1300.0 + 200.0);
    assert_eq!(rm.min_intensity().unwrap(), 1.0 + 499.0);
    assert_eq!(rm.max_intensity().unwrap(), 600.0); // was 47110.0
    assert!(rm.is_dim_empty(MSDim::Mobility).unwrap());

    // if no dimensions overlap...
    let rt = RangeManager::new(&[MSDim::Rt]).unwrap();
    assert!(matches!(rmi.push_into(&rt), Err(Error::InvalidRange(_))));
}

#[test]
fn range_manager_clamp_to() {
    // START_SECTION(template<typename... RangeBasesOther> void clampTo(const RangeManager<RangeBasesOther...>& rhs))
    let mut rm = rm_update_ranges();
    let mut rmi = RangeManager::new(&[MSDim::Mz, MSDim::Intensity]).unwrap();
    rmi.extend_mz(700.0).unwrap(); // clamp left
    rmi.extend_mz(2000.0).unwrap();
    rmi.extend_intensity(-10.0).unwrap(); // clamp to empty
    rmi.extend_intensity(-9.0).unwrap();
    rm.clamp_to(&rmi).unwrap();
    assert_eq!(rm.min_rt().unwrap(), 2.0); // should be untouched (since rmi.RT is empty)
    assert_eq!(rm.max_rt().unwrap(), 100.0);
    assert_eq!(rm.min_mz().unwrap(), 500.0 + 200.0);
    assert_eq!(rm.max_mz().unwrap(), 1300.0 + 0.0);
    assert!(rm.range_for_dim(MSDim::Intensity).unwrap().is_empty());
    assert!(rm.is_dim_empty(MSDim::Mobility).unwrap());

    // if no dimensions overlap...
    let rt = RangeManager::new(&[MSDim::Rt]).unwrap();
    assert!(matches!(rmi.clamp_to(&rt), Err(Error::InvalidRange(_))));
}

#[test]
fn range_manager_get_range_for_dim() {
    // START_SECTION(RangeBase& getRangeForDim(MSDim dim))
    let rm = rm_update_ranges();
    let rt = *rm.range_for_dim(MSDim::Rt).unwrap();
    let mz = *rm.range_for_dim(MSDim::Mz).unwrap();
    let int = *rm.range_for_dim(MSDim::Intensity).unwrap();
    let im = *rm.range_for_dim(MSDim::Mobility).unwrap();
    assert_eq!(rt.min().unwrap(), 2.0);
    assert_eq!(mz.min().unwrap(), 500.0);
    assert_eq!(rt.max().unwrap(), 100.0);
    assert_eq!(mz.max().unwrap(), 1300.0);
    assert_eq!(int.min().unwrap(), 1.0);
    assert_eq!(int.max().unwrap(), 47110.0);
    assert!(!rt.is_empty());
    assert!(im.is_empty());
}

#[test]
fn range_manager_clear_ranges() {
    // START_SECTION((void clearRanges()))
    let mut rm = rm_update_ranges();
    assert_rm_reference(&rm);
    assert!(!rm.is_dim_empty(MSDim::Rt).unwrap());
    assert!(!rm.is_dim_empty(MSDim::Mz).unwrap());
    assert!(!rm.is_dim_empty(MSDim::Intensity).unwrap());
    assert!(rm.is_dim_empty(MSDim::Mobility).unwrap());

    rm.clear_ranges();
    assert!(rm.is_dim_empty(MSDim::Rt).unwrap());
    assert!(rm.is_dim_empty(MSDim::Mz).unwrap());
    assert!(rm.is_dim_empty(MSDim::Intensity).unwrap());
    assert!(rm.is_dim_empty(MSDim::Mobility).unwrap());
}

#[test]
fn range_manager_print_range() {
    // START_SECTION(void printRange(std::ostream& out) const)
    let mut rm = RangeManager::experiment();
    rm.extend_rt(1.0).unwrap();
    rm.extend_mz(2.0).unwrap();
    rm.extend_intensity(3.0).unwrap();
    rm.extend_mobility(4.0).unwrap();
    assert_eq!(
        rm.to_string(),
        "rt: [1, 1]\nmz: [2, 2]\nintensity: [3, 3]\nmobility: [4, 4]\n"
    );
}

// ---------------------------------------------------------------------------
// SpectrumRangeManager_test.cpp
// ---------------------------------------------------------------------------

#[test]
fn spectrum_range_manager_constructor() {
    // START_SECTION((SpectrumRangeManager()))
    let srm = SpectrumRangeManager::new();
    assert!(srm.ms_levels().is_empty());
    assert_eq!(srm, SpectrumRangeManager::default());
}

#[test]
fn spectrum_range_manager_destructor() {
    // START_SECTION((~SpectrumRangeManager()))
    drop(SpectrumRangeManager::new());
}

#[test]
fn spectrum_range_manager_extend_rt() {
    // START_SECTION((void extendRT(double rt, UInt ms_level = 0)))
    let mut srm = SpectrumRangeManager::new();
    // ms_level 0 -> global (base) range, does not register an MS level
    srm.extend_rt(100.0, 0).unwrap();
    srm.extend_rt(200.0, 0).unwrap();
    assert_eq!(srm.global().min_rt().unwrap(), 100.0);
    assert_eq!(srm.global().max_rt().unwrap(), 200.0);
    assert!(srm.ms_levels().is_empty());

    // ms_level 2 -> level-specific range, does not touch the global range
    srm.extend_rt(10.0, 2).unwrap();
    srm.extend_rt(20.0, 2).unwrap();
    assert_eq!(srm.by_ms_level(2).unwrap().min_rt().unwrap(), 10.0);
    assert_eq!(srm.by_ms_level(2).unwrap().max_rt().unwrap(), 20.0);
    assert_eq!(srm.global().min_rt().unwrap(), 100.0); // global unchanged
    assert_eq!(srm.global().max_rt().unwrap(), 200.0);
}

#[test]
fn spectrum_range_manager_extend_mz() {
    // START_SECTION((void extendMZ(double mz, UInt ms_level = 0)))
    let mut srm = SpectrumRangeManager::new();
    srm.extend_mz(500.0, 0).unwrap();
    srm.extend_mz(600.0, 0).unwrap();
    assert_eq!(srm.global().min_mz().unwrap(), 500.0);
    assert_eq!(srm.global().max_mz().unwrap(), 600.0);

    srm.extend_mz(300.0, 2).unwrap();
    srm.extend_mz(400.0, 2).unwrap();
    assert_eq!(srm.by_ms_level(2).unwrap().min_mz().unwrap(), 300.0);
    assert_eq!(srm.by_ms_level(2).unwrap().max_mz().unwrap(), 400.0);
    assert_eq!(srm.global().min_mz().unwrap(), 500.0); // global unchanged
}

#[test]
fn spectrum_range_manager_ms_levels() {
    // START_SECTION((std::set<UInt> getMSLevels() const))
    let mut srm = SpectrumRangeManager::new();
    assert!(srm.ms_levels().is_empty());

    srm.extend_rt(1.0, 2).unwrap();
    srm.extend_rt(1.0, 3).unwrap();
    srm.extend_mz(1.0, 2).unwrap(); // level 2 already present

    let levels = srm.ms_levels();
    assert_eq!(levels.len(), 2);
    assert!(levels.contains(&2));
    assert!(levels.contains(&3));

    // global extends (ms_level 0) do not register a level
    srm.extend_rt(5.0, 0).unwrap();
    assert_eq!(srm.ms_levels().len(), 2);
}

#[test]
fn spectrum_range_manager_by_ms_level() {
    // START_SECTION((const BaseType& byMSLevel(UInt ms_level = 0) const))
    let mut srm = SpectrumRangeManager::new();
    srm.extend_rt(10.0, 2).unwrap();
    srm.extend_rt(20.0, 2).unwrap();
    assert_eq!(srm.by_ms_level(2).unwrap().min_rt().unwrap(), 10.0);

    // global ranges live in the base class, not in the per-level map:
    // byMSLevel(0) therefore has no entry and throws — the port returns None
    srm.extend_rt(100.0, 0).unwrap(); // global
    assert!(srm.by_ms_level(0).is_none());

    // unknown level throws — the port returns None
    assert!(srm.by_ms_level(99).is_none());
}

#[test]
fn spectrum_range_manager_extend_other() {
    // START_SECTION((void extend(const BaseType& other, UInt ms_level = 0)))
    let mut srm = SpectrumRangeManager::new();
    let mut other = RangeManager::spectrum_manager();
    other.extend_rt(1000.0).unwrap();
    other.extend_rt(2000.0).unwrap();

    // into a specific level
    srm.extend(&other, 5).unwrap();
    assert_eq!(srm.by_ms_level(5).unwrap().min_rt().unwrap(), 1000.0);
    assert_eq!(srm.by_ms_level(5).unwrap().max_rt().unwrap(), 2000.0);

    // into the global range (ms_level 0)
    srm.extend_rt(100.0, 0).unwrap();
    srm.extend(&other, 0).unwrap();
    assert_eq!(srm.global().min_rt().unwrap(), 100.0);
    assert_eq!(srm.global().max_rt().unwrap(), 2000.0);
}

#[test]
fn spectrum_range_manager_extend_spectrum() {
    // START_SECTION((void extendUnsafe(const MSSpectrum& spectrum, UInt ms_level = 0)))
    let mut srm = SpectrumRangeManager::new();
    let s = MSSpectrum::from_peaks(vec![Peak1D::new(700.0, 5.0), Peak1D::new(800.0, 9.0)]);

    srm.extend_spectrum(&s, 4).unwrap();
    assert_eq!(srm.by_ms_level(4).unwrap().min_mz().unwrap(), 700.0);
    assert_eq!(srm.by_ms_level(4).unwrap().max_mz().unwrap(), 800.0);
}

#[test]
fn spectrum_range_manager_clear_ranges() {
    // START_SECTION((void clearRanges()))
    let mut srm = SpectrumRangeManager::new();
    srm.extend_rt(100.0, 0).unwrap(); // global
    srm.extend_rt(10.0, 2).unwrap(); // level 2
    assert!(!srm.ms_levels().is_empty());

    srm.clear_ranges();
    assert!(srm.ms_levels().is_empty());
    assert!(srm.global().range_for_dim(MSDim::Rt).unwrap().is_empty());
}

#[test]
fn spectrum_range_manager_copy_constructor() {
    // START_SECTION((SpectrumRangeManager(const SpectrumRangeManager& source)))
    let mut srm = SpectrumRangeManager::new();
    srm.extend_rt(100.0, 0).unwrap(); // global
    srm.extend_rt(10.0, 2).unwrap();
    srm.extend_rt(20.0, 2).unwrap();

    let copy = srm.clone();
    assert_eq!(copy.global().min_rt().unwrap(), 100.0); // global copied
    assert!(copy.ms_levels().contains(&2)); // level copied
    assert_eq!(copy.by_ms_level(2).unwrap().max_rt().unwrap(), 20.0);
    assert_eq!(copy, srm);
}

// ---------------------------------------------------------------------------
// Native: dimension sets, finiteness, atomicity
// ---------------------------------------------------------------------------

#[test]
fn dimension_sets_and_presets() {
    assert!(matches!(
        RangeManager::new(&[]),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        RangeManager::new(&[MSDim::Rt, MSDim::Rt]),
        Err(Error::InvalidValue(_))
    ));
    let spectrum = RangeManager::spectrum();
    assert_eq!(
        spectrum.dims().collect::<Vec<_>>(),
        [MSDim::Mz, MSDim::Intensity, MSDim::Mobility]
    );
    assert!(!spectrum.has_dim(MSDim::Rt));
    assert!(matches!(spectrum.min_rt(), Err(Error::InvalidValue(_))));
    assert!(matches!(
        spectrum.range_for_dim(MSDim::Rt),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(
        RangeManager::chromatogram().dims().collect::<Vec<_>>(),
        [MSDim::Rt, MSDim::Intensity]
    );
    assert_eq!(
        RangeManager::mobilogram().dims().collect::<Vec<_>>(),
        [MSDim::Mobility, MSDim::Intensity]
    );
    assert_eq!(
        RangeManager::chromatogram_manager()
            .dims()
            .collect::<Vec<_>>(),
        [MSDim::Rt, MSDim::Intensity, MSDim::Mz]
    );
    assert_eq!(
        RangeManager::spectrum_manager().dims().collect::<Vec<_>>(),
        [MSDim::Mz, MSDim::Intensity, MSDim::Mobility, MSDim::Rt]
    );
    // Display follows declaration order.
    let mut cm = RangeManager::chromatogram_manager();
    cm.extend_rt(1.5).unwrap();
    assert_eq!(
        cm.to_string(),
        "rt: [1.5, 1.5]\nintensity: [, ]\nmz: [, ]\n"
    );
    // Equality ignores declaration order.
    let mut a = RangeManager::new(&[MSDim::Rt, MSDim::Mz]).unwrap();
    let mut b = RangeManager::new(&[MSDim::Mz, MSDim::Rt]).unwrap();
    a.extend_rt(1.0).unwrap();
    b.extend_rt(1.0).unwrap();
    assert_eq!(a, b);
    // A set that only partially overlaps is not equal.
    assert_ne!(a, RangeManager::new(&[MSDim::Rt]).unwrap());
    // clear_dim on an absent dimension is a no-op, as the source.
    let before = a;
    a.clear_dim(MSDim::Mobility);
    assert_eq!(a, before);
    a.clear_dim(MSDim::Rt);
    assert!(a.is_dim_empty(MSDim::Rt).unwrap());
    // Labels and the dimension list.
    assert_eq!(MSDim::Intensity.label(), "intensity");
    assert_eq!(MSDim::ALL.len(), 4);
}

#[test]
fn typed_accessor_families_delegate_to_their_dimension() {
    let mut rm = RangeManager::experiment();
    rm.set_min_rt(1.0).unwrap();
    rm.set_max_rt(2.0).unwrap();
    rm.set_min_mz(3.0).unwrap();
    rm.set_max_mz(4.0).unwrap();
    rm.set_min_intensity(5.0).unwrap();
    rm.set_max_intensity(6.0).unwrap();
    rm.set_min_mobility(7.0).unwrap();
    rm.set_max_mobility(8.0).unwrap();
    assert_eq!((rm.min_rt().unwrap(), rm.max_rt().unwrap()), (1.0, 2.0));
    assert_eq!((rm.min_mz().unwrap(), rm.max_mz().unwrap()), (3.0, 4.0));
    assert_eq!(
        (rm.min_intensity().unwrap(), rm.max_intensity().unwrap()),
        (5.0, 6.0)
    );
    assert_eq!(
        (rm.min_mobility().unwrap(), rm.max_mobility().unwrap()),
        (7.0, 8.0)
    );
    assert!(rm.contains_rt(1.5).unwrap());
    assert!(!rm.contains_mz(2.0).unwrap());
    assert!(rm.contains_intensity(6.0).unwrap());
    assert!(!rm.contains_mobility(9.0).unwrap());
    assert!(rm.contains_rt_range(&rb(1.0, 2.0)).unwrap());
    assert!(!rm.contains_mz_range(&rb(3.0, 5.0)).unwrap());
    assert!(rm.contains_intensity_range(&rb(5.5, 5.5)).unwrap());
    assert!(!rm.contains_mobility_range(&rb(0.0, 8.0)).unwrap());
    rm.extend_range(MSDim::Mz, &rb(0.0, 10.0)).unwrap();
    assert_eq!((rm.min_mz().unwrap(), rm.max_mz().unwrap()), (0.0, 10.0));
    assert_eq!(
        rm.set_min(MSDim::Mz, 20.0)
            .map(|()| rm.max_mz().unwrap())
            .unwrap(),
        20.0
    );
    assert!(rm.contains_value(MSDim::Mz, 20.0).unwrap());
    assert!(rm.contains_range(MSDim::Mz, &rb(20.0, 20.0)).unwrap());
    // min_span_if_singular widens singular dimensions only.
    rm.min_span_if_singular(4.0).unwrap();
    assert_eq!((rm.min_mz().unwrap(), rm.max_mz().unwrap()), (18.0, 22.0));
    assert_eq!((rm.min_rt().unwrap(), rm.max_rt().unwrap()), (1.0, 2.0));
    let mut single = rb(5.0, 5.0);
    single.min_span_if_singular(1.0).unwrap();
    assert_eq!(single.span(), Some(1.0));
    let mut empty = RangeBase::new();
    empty.min_span_if_singular(1.0).unwrap();
    assert!(empty.is_empty());
    // The unsafe variants report overlap without an error.
    let rt_only = RangeManager::new(&[MSDim::Rt]).unwrap();
    let mut mz_only = RangeManager::new(&[MSDim::Mz]).unwrap();
    assert!(!mz_only.assign_unsafe(&rt_only));
    assert!(!mz_only.extend_unsafe(&rt_only));
    assert!(!mz_only.push_into_unsafe(&rt_only).unwrap());
    assert!(!mz_only.clamp_to_unsafe(&rt_only));
    assert!(matches!(
        mz_only.assign(&rt_only),
        Err(Error::InvalidRange(_))
    ));
    assert!(matches!(
        mz_only.extend(&rt_only),
        Err(Error::InvalidRange(_))
    ));
    assert!(mz_only.extend_unsafe(&rm));
    assert_eq!(
        mz_only.range_for_dim(MSDim::Mz).unwrap(),
        rm.range_for_dim(MSDim::Mz).unwrap()
    );
    assert_eq!(
        rm.range_for_dim(MSDim::Rt).unwrap().non_empty_range(),
        (1.0, 2.0)
    );
    assert_eq!(
        rm.range_for_dim(MSDim::Mz).unwrap().non_empty_range(),
        (18.0, 22.0)
    );
}

#[test]
fn non_finite_input_is_rejected_and_leaves_values_unchanged() {
    assert!(matches!(
        RangeBase::singular(f64::NAN),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        RangeBase::from_min_max(f64::NEG_INFINITY, 1.0),
        Err(Error::InvalidValue(_))
    ));
    let mut b = rb(4.0, 6.0);
    for op in [
        RangeBase::set_min,
        RangeBase::set_max,
        RangeBase::extend_value,
        RangeBase::extend_left_right,
        RangeBase::min_span_if_singular,
        RangeBase::scale_by,
        RangeBase::shift,
    ] {
        assert!(matches!(op(&mut b, f64::NAN), Err(Error::InvalidValue(_))));
        assert!(matches!(
            op(&mut b, f64::INFINITY),
            Err(Error::InvalidValue(_))
        ));
        assert_eq!(b, rb(4.0, 6.0));
    }
    // Arithmetic that leaves the finite domain is refused atomically.
    let mut wide = rb(f64::MIN, f64::MAX);
    assert!(matches!(wide.scale_by(2.0), Err(Error::InvalidRange(_))));
    assert!(matches!(wide.shift(f64::MAX), Err(Error::InvalidRange(_))));
    assert_eq!(wide, rb(f64::MIN, f64::MAX));
    assert_eq!(wide.span(), Some(f64::INFINITY));
    assert_eq!(RangeBase::new().non_empty_range(), (f64::MIN, f64::MAX));
    let mut rm = rm_update_ranges();
    let before = rm;
    assert!(matches!(
        rm.extend_rt(f64::NAN),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(rm.scale_by(f64::NAN), Err(Error::InvalidValue(_))));
    assert_eq!(rm, before);
    // A manager scale that overflows one dimension changes none.
    let mut mixed = RangeManager::experiment();
    mixed
        .extend_range(MSDim::Rt, &rb(f64::MIN, f64::MAX))
        .unwrap();
    mixed.extend_range(MSDim::Mz, &rb(1.0, 2.0)).unwrap();
    let before = mixed;
    assert!(matches!(mixed.scale_by(3.0), Err(Error::InvalidRange(_))));
    assert_eq!(mixed, before);
    let mut srm = SpectrumRangeManager::new();
    assert!(matches!(
        srm.extend_rt(f64::NAN, 3),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        srm.extend_mz(f64::INFINITY, 3),
        Err(Error::InvalidValue(_))
    ));
    assert!(srm.ms_levels().is_empty());
    // Display of an empty range.
    assert_eq!(RangeBase::new().to_string(), "[, ]");
    assert_eq!(rb(-0.5, 2.25).to_string(), "[-0.5, 2.25]");
}

#[test]
fn extend_keeps_the_first_equal_endpoint_including_signed_zero() {
    let mut r = RangeBase::singular(0.0).unwrap();
    r.extend_value(-0.0).unwrap();
    assert!(r.min().unwrap().is_sign_positive());
    let mut n = RangeBase::singular(-0.0).unwrap();
    n.extend(&RangeBase::singular(0.0).unwrap());
    assert!(n.min().unwrap().is_sign_negative());
    assert!(n.max().unwrap().is_sign_negative());
}

// ---------------------------------------------------------------------------
// Native: on-demand container ranges
// ---------------------------------------------------------------------------

fn spectrum(rt: f64, level: u32, peaks: &[(f64, f32)]) -> MSSpectrum {
    MSSpectrum {
        rt,
        ms_level: level,
        peaks: peaks.iter().map(|&(mz, i)| Peak1D::new(mz, i)).collect(),
        ..Default::default()
    }
}

fn chromatogram(product_mz: f64, points: &[(f64, f32)]) -> MSChromatogram {
    let mut c = MSChromatogram::from_peaks(
        points
            .iter()
            .map(|&(rt, i)| ChromatogramPeak::new(rt, i))
            .collect(),
    );
    c.product.mz = product_mz;
    c
}

#[test]
fn spectrum_range_manager_from_peaks_and_scalar_drift_time() {
    // MSSpectrum.cpp:586-604: m/z and intensity from peaks; mobility from the
    // scalar drift time only when it is not DRIFTTIME_NOT_SET (-1).
    let mut s = spectrum(12.0, 1, &[(100.0, 5.0), (300.0, 1.0), (200.0, 9.0)]);
    let r = s.range_manager().unwrap();
    assert_eq!(
        r.dims().collect::<Vec<_>>(),
        [MSDim::Mz, MSDim::Intensity, MSDim::Mobility]
    );
    assert_eq!((r.min_mz().unwrap(), r.max_mz().unwrap()), (100.0, 300.0));
    assert_eq!(
        (r.min_intensity().unwrap(), r.max_intensity().unwrap()),
        (1.0, 9.0)
    );
    assert!(r.is_dim_empty(MSDim::Mobility).unwrap());
    assert!(!r.has_dim(MSDim::Rt)); // RT is not part of a spectrum's range
    assert_eq!(r.has_range(), HasRangeType::Some);

    s.drift_time = 0.75;
    let r = s.range_manager().unwrap();
    assert_eq!(
        (r.min_mobility().unwrap(), r.max_mobility().unwrap()),
        (0.75, 0.75)
    );
    assert_eq!(r.has_range(), HasRangeType::All);

    // The source compares against -1 exactly; other negative values count.
    s.drift_time = -2.0;
    assert_eq!(s.range_manager().unwrap().min_mobility().unwrap(), -2.0);
    s.drift_time = f64::NAN;
    assert!(matches!(s.range_manager(), Err(Error::InvalidValue(_))));

    // Empty spectrum: every dimension empty, no error.
    let empty = MSSpectrum::default();
    let r = empty.range_manager().unwrap();
    assert_eq!(r.has_range(), HasRangeType::None);
    assert!(matches!(r.min_mz(), Err(Error::InvalidRange(_))));

    // Mutation is visible immediately; there is no cache to refresh.
    s.drift_time = -1.0;
    s.peaks[0].mz = 1000.0;
    assert_eq!(s.range_manager().unwrap().max_mz().unwrap(), 1000.0);
}

#[test]
fn spectrum_range_manager_prefers_the_ion_mobility_array() {
    // MSSpectrum.cpp:592-600: an IM frame takes every value of the first
    // ion-mobility float data array and ignores the scalar drift time.
    let mut s = spectrum(1.0, 1, &[(100.0, 1.0), (200.0, 2.0), (300.0, 3.0)]);
    s.drift_time = 99.0;
    s.float_data_arrays
        .push(DataArray::new("unrelated", vec![7.0, 8.0, 9.0]));
    s.float_data_arrays.push(DataArray::new(
        "raw inverse reduced ion mobility array",
        vec![1.2, 0.8, 1.0],
    ));
    let r = s.range_manager().unwrap();
    assert_eq!(r.min_mobility().unwrap(), f64::from(0.8f32));
    assert_eq!(r.max_mobility().unwrap(), f64::from(1.2f32));

    // The UserParam fallbacks are prefixes, as in IMDataArrayUtils::getIMUnit.
    for name in [
        "Ion Mobility",
        "Ion Mobility (MS:1002815)",
        "inverse reduced ion mobility",
        "mean inverse reduced ion mobility array",
        "mean ion mobility drift time array",
        "deconvoluted ion mobility drift time array",
    ] {
        let mut t = spectrum(1.0, 1, &[(100.0, 1.0), (200.0, 2.0)]);
        t.float_data_arrays
            .push(DataArray::new(name, vec![3.0, 4.0]));
        let r = t.range_manager().unwrap();
        assert_eq!(
            (r.min_mobility().unwrap(), r.max_mobility().unwrap()),
            (3.0, 4.0),
            "{name}"
        );
    }
    // An IM array of non-finite values is rejected.
    let mut bad = spectrum(1.0, 1, &[(100.0, 1.0)]);
    bad.float_data_arrays
        .push(DataArray::new("Ion Mobility", vec![f32::NAN]));
    assert!(matches!(bad.range_manager(), Err(Error::InvalidValue(_))));
    // An empty IM array (placeholder) yields an empty mobility range.
    let mut placeholder = spectrum(1.0, 1, &[(100.0, 1.0)]);
    placeholder.drift_time = 5.0;
    placeholder
        .float_data_arrays
        .push(DataArray::new("Ion Mobility", vec![]));
    assert!(
        placeholder
            .range_manager()
            .unwrap()
            .is_dim_empty(MSDim::Mobility)
            .unwrap()
    );
}

#[test]
fn chromatogram_and_mobilogram_range_managers() {
    let c = chromatogram(505.0, &[(3.0, 2.0), (1.0, 8.0), (2.0, 4.0)]);
    let r = c.range_manager().unwrap();
    assert_eq!(r.dims().collect::<Vec<_>>(), [MSDim::Rt, MSDim::Intensity]);
    assert_eq!((r.min_rt().unwrap(), r.max_rt().unwrap()), (1.0, 3.0));
    assert_eq!(
        (r.min_intensity().unwrap(), r.max_intensity().unwrap()),
        (2.0, 8.0)
    );
    assert!(!r.has_dim(MSDim::Mz)); // the product m/z is not part of a chromatogram's range
    assert_eq!(
        MSChromatogram::default()
            .range_manager()
            .unwrap()
            .has_range(),
        HasRangeType::None
    );
    let mut nan = c.clone();
    nan.peaks[1].rt = f64::NAN;
    assert!(nan.range_manager().is_err());

    let m = Mobilogram::from_peaks(vec![
        MobilityPeak1D::new(0.9, 3.0),
        MobilityPeak1D::new(0.7, 1.0),
        MobilityPeak1D::new(1.1, 2.0),
    ]);
    let r = m.range_manager().unwrap();
    assert_eq!(
        r.dims().collect::<Vec<_>>(),
        [MSDim::Mobility, MSDim::Intensity]
    );
    assert_eq!(
        (r.min_mobility().unwrap(), r.max_mobility().unwrap()),
        (0.7, 1.1)
    );
    assert_eq!(
        (r.min_intensity().unwrap(), r.max_intensity().unwrap()),
        (1.0, 3.0)
    );
    assert_eq!(
        Mobilogram::default().range_manager().unwrap().has_range(),
        HasRangeType::None
    );
}

fn three_role_experiment() -> MSExperiment {
    let mut ms1 = spectrum(10.0, 1, &[(400.0, 100.0), (600.0, 300.0)]);
    ms1.drift_time = 0.5;
    let mut ms2 = spectrum(20.0, 2, &[(150.0, 50.0), (250.0, 900.0)]);
    ms2.float_data_arrays
        .push(DataArray::new("Ion Mobility", vec![0.2, 1.4]));
    MSExperiment {
        spectra: vec![
            ms1,
            ms2,
            spectrum(-5.0, 1, &[]), // peakless: contributes RT only
            spectrum(30.0, 3, &[(700.0, 1.0)]),
        ],
        chromatograms: vec![
            chromatogram(1200.0, &[(50.0, 5000.0), (60.0, 2.0)]),
            chromatogram(80.0, &[]), // pointless: contributes product m/z only
        ],
        ..Default::default()
    }
}

#[test]
fn experiment_spectrum_role_is_global_and_per_level() {
    // MSExperiment.cpp:689-700.
    let exp = three_role_experiment();
    let srm = exp.spectrum_range_manager().unwrap();
    let g = srm.global();
    assert_eq!(
        g.dims().collect::<Vec<_>>(),
        [MSDim::Mz, MSDim::Intensity, MSDim::Mobility, MSDim::Rt]
    );
    assert_eq!((g.min_rt().unwrap(), g.max_rt().unwrap()), (-5.0, 30.0)); // peakless RT counts
    assert_eq!((g.min_mz().unwrap(), g.max_mz().unwrap()), (150.0, 700.0));
    assert_eq!(
        (g.min_intensity().unwrap(), g.max_intensity().unwrap()),
        (1.0, 900.0)
    );
    assert_eq!(
        (g.min_mobility().unwrap(), g.max_mobility().unwrap()),
        (f64::from(0.2f32), f64::from(1.4f32))
    );
    assert_eq!(srm.ms_levels().into_iter().collect::<Vec<_>>(), [1, 2, 3]);
    let l1 = srm.by_ms_level(1).unwrap();
    assert_eq!((l1.min_rt().unwrap(), l1.max_rt().unwrap()), (-5.0, 10.0));
    assert_eq!((l1.min_mz().unwrap(), l1.max_mz().unwrap()), (400.0, 600.0));
    assert_eq!(
        (l1.min_mobility().unwrap(), l1.max_mobility().unwrap()),
        (0.5, 0.5)
    );
    let l2 = srm.by_ms_level(2).unwrap();
    assert_eq!((l2.min_rt().unwrap(), l2.max_rt().unwrap()), (20.0, 20.0));
    assert_eq!(
        (l2.min_intensity().unwrap(), l2.max_intensity().unwrap()),
        (50.0, 900.0)
    );
    let l3 = srm.by_ms_level(3).unwrap();
    assert_eq!((l3.min_mz().unwrap(), l3.max_mz().unwrap()), (700.0, 700.0));
    assert!(srm.by_ms_level(4).is_none());

    // Level 0 spectra extend the global ranges only and register no level.
    let mut optical = MSExperiment::default();
    let mut s = spectrum(3.0, 0, &[(1.0, 1.0)]);
    s.instrument_settings.scan_mode = openms::metadata::ScanMode::Emission;
    optical.spectra.push(s);
    let srm = optical.spectrum_range_manager().unwrap();
    assert!(srm.ms_levels().is_empty());
    assert_eq!(srm.global().min_rt().unwrap(), 3.0);

    // An empty experiment has no ranges and no levels.
    let srm = MSExperiment::default().spectrum_range_manager().unwrap();
    assert_eq!(srm.global().has_range(), HasRangeType::None);
    assert!(srm.ms_levels().is_empty());

    // A malformed spectrum fails the whole query.
    let mut bad = three_role_experiment();
    bad.spectra[3].peaks[0].intensity = f32::INFINITY;
    assert!(bad.spectrum_range_manager().is_err());
}

#[test]
fn experiment_chromatogram_role_uses_points_and_product_mz() {
    // MSExperiment.cpp:704-715 with MSChromatogram::getMZ() = product m/z.
    let exp = three_role_experiment();
    let c = exp.chromatogram_range_manager().unwrap();
    assert_eq!(
        c.dims().collect::<Vec<_>>(),
        [MSDim::Rt, MSDim::Intensity, MSDim::Mz]
    );
    assert_eq!((c.min_rt().unwrap(), c.max_rt().unwrap()), (50.0, 60.0));
    assert_eq!(
        (c.min_intensity().unwrap(), c.max_intensity().unwrap()),
        (2.0, 5000.0)
    );
    assert_eq!((c.min_mz().unwrap(), c.max_mz().unwrap()), (80.0, 1200.0)); // pointless chromatogram's product m/z
    assert!(!c.has_dim(MSDim::Mobility));
    // Without any product, the source contributes m/z 0.
    let mut none = MSExperiment::default();
    none.chromatograms.push(MSChromatogram::default());
    let c = none.chromatogram_range_manager().unwrap();
    assert_eq!((c.min_mz().unwrap(), c.max_mz().unwrap()), (0.0, 0.0));
    assert!(c.is_dim_empty(MSDim::Rt).unwrap());
    assert_eq!(
        MSExperiment::default()
            .chromatogram_range_manager()
            .unwrap()
            .has_range(),
        HasRangeType::None
    );
}

#[test]
fn experiment_combined_role_merges_spectra_then_chromatograms() {
    // MSExperiment.cpp:718-719.
    let exp = three_role_experiment();
    let all = exp.combined_range_manager().unwrap();
    assert_eq!(all.dims().collect::<Vec<_>>(), MSDim::ALL.to_vec());
    assert_eq!((all.min_rt().unwrap(), all.max_rt().unwrap()), (-5.0, 60.0));
    assert_eq!(
        (all.min_mz().unwrap(), all.max_mz().unwrap()),
        (80.0, 1200.0)
    );
    assert_eq!(
        (all.min_intensity().unwrap(), all.max_intensity().unwrap()),
        (1.0, 5000.0)
    );
    assert_eq!(
        (all.min_mobility().unwrap(), all.max_mobility().unwrap()),
        (f64::from(0.2f32), f64::from(1.4f32))
    );
    assert_eq!(all.has_range(), HasRangeType::All);
    // The combined result contains each role.
    assert!(
        all.contains_all(exp.spectrum_range_manager().unwrap().global())
            .unwrap()
    );
    assert!(
        all.contains_all(&exp.chromatogram_range_manager().unwrap())
            .unwrap()
    );
    // The existing spectrum-only query and the combined one agree with the
    // legacy three-dimensional summaries where they overlap.
    let legacy = exp.combined_ranges().unwrap();
    assert_eq!(legacy.rt.unwrap().min, all.min_rt().unwrap());
    assert_eq!(legacy.mz.unwrap().max, all.max_mz().unwrap());
    // Spectra first: an equal zero endpoint keeps the spectrum's sign.
    let mut zero = MSExperiment::default();
    zero.spectra.push(spectrum(0.0, 1, &[(0.0, 0.0)]));
    zero.chromatograms.push(chromatogram(-0.0, &[(-0.0, -0.0)]));
    let z = zero.combined_range_manager().unwrap();
    assert!(z.min_rt().unwrap().is_sign_positive());
    assert!(z.min_mz().unwrap().is_sign_positive());
    assert!(z.min_intensity().unwrap().is_sign_positive());
    // Empty experiment.
    assert_eq!(
        MSExperiment::default()
            .combined_range_manager()
            .unwrap()
            .has_range(),
        HasRangeType::None
    );
    // A malformed chromatogram fails the whole query.
    let mut bad = three_role_experiment();
    bad.chromatograms[0].product.mz = f64::NAN;
    assert!(bad.combined_range_manager().is_err());
    assert!(bad.chromatogram_range_manager().is_err());
    assert!(bad.spectrum_range_manager().is_ok());
}
