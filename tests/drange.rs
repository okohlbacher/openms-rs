// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Port of `DRange_test.cpp` (core SDK bc9cc12): one test per
//! `START_SECTION`, cited by source line. Literals are transcribed (tier 3).

use openms::Error;
use openms::data_structures::{
    DIntervalBase, DPosition, DPosition1, DRange, DRange1, DRange2, DRangeIntersection as Kind,
};
use std::hash::{DefaultHasher, Hash, Hasher};

fn hash(value: &impl Hash) -> u64 {
    let mut h = DefaultHasher::new();
    value.hash(&mut h);
    h.finish()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1.),
        "{actual:.17} != {expected:.17}"
    );
}
// Source globals (L27-37): p1 = (-1, -2), p2 = (3, 4), one = (1, 1), two = (2, 2).
fn p1() -> DPosition<2> {
    DPosition::xy(-1., -2.)
}
fn p2() -> DPosition<2> {
    DPosition::xy(3., 4.)
}
fn one() -> DPosition<2> {
    DPosition::xy(1., 1.)
}
fn two() -> DPosition<2> {
    DPosition::xy(2., 2.)
}
// Source global (L64): r = DRange(p1, p2).
fn r() -> DRange2 {
    DRange::new(p1(), p2())
}
fn corners(r: &DRange2) -> [f64; 4] {
    [
        r.min_position()[0],
        r.min_position()[1],
        r.max_position()[0],
        r.max_position()[1],
    ]
}

// L46 DRange(): the source comment promises zeros, the body yields the empty sentinel.
#[test]
fn default_constructor_is_the_empty_sentinel() {
    let r = DRange2::default();
    assert!(r.is_empty());
    assert_eq!(r, DRange2::empty());
    assert_ne!(r, DRange2::zero());
    assert_eq!(r.min_position()[0], f64::MAX);
    assert_eq!(r.max_position()[1], f64::MIN);
    assert_eq!(DRange2::DIMENSION, 2);
}

// L51 ~DRange()
#[test]
fn destructor_is_ordinary_drop() {
    let r = Box::new(DRange2::default());
    assert!(r.is_empty());
    drop(r);
}

// L55 DRange(lower, upper)
#[test]
fn corner_constructor() {
    assert_eq!(corners(&r()), [-1., -2., 3., 4.]);
}

// L67 copy constructor
#[test]
fn copy_constructor() {
    let r2 = r();
    assert_eq!(corners(&r2), [-1., -2., 3., 4.]);
    assert_eq!(r2, r());
}

// L75 DRange(const Base&)
#[test]
fn from_base_interval() {
    let ib: DIntervalBase<2> = r().into();
    let r2 = DRange::from(ib);
    assert_eq!(corners(&r2), [-1., -2., 3., 4.]);
    assert_eq!(*r2.base(), ib);
}

// L84 operator=(const Base&)
#[test]
fn assign_from_base_interval() {
    let ib: DIntervalBase<2> = r().into();
    let mut r2 = DRange2::default();
    assert!(r2.is_empty());
    r2 = ib.into();
    assert_eq!(corners(&r2), [-1., -2., 3., 4.]);
    *r2.base_mut() = DIntervalBase::zero();
    assert_eq!(r2, DRange2::zero());
}

// L94 operator=(const DRange&)
#[test]
fn assign_from_range() {
    let mut r2 = DRange2::default();
    assert!(r2.is_empty());
    r2 = r();
    assert_eq!(corners(&r2), [-1., -2., 3., 4.]);
}

// L103 DRange(minx, miny, maxx, maxy), including min > max normalization
#[test]
fn xy_constructor_normalizes() {
    let r2 = DRange::xy(1., 2., 3., 4.);
    assert_eq!(corners(&r2), [1., 2., 3., 4.]);
    let r = DRange::xy(2., 3., -2., -3.);
    assert_eq!(*r.min_position(), DPosition::xy(-2., -3.));
    assert_eq!(*r.max_position(), DPosition::xy(2., 3.));
}

// L114 operator==(const DRange&)
#[test]
fn equality_with_range() {
    let r = r();
    let mut r2 = r;
    assert!(r == r2);
    r2.set_min_x(0.);
    assert!(r != r2);
    r2.set_min_x(r.min_position()[0]);
    assert!(r == r2);
    r2.set_max_y(0.);
    assert!(r != r2);
    r2.set_max_y(r.max_position()[1]);
    assert!(r == r2);
}

// L127 operator==(const Base&)
#[test]
fn equality_with_base_interval() {
    let r = r();
    let mut r2: DIntervalBase<2> = r.into();
    assert!(r == r2);
    r2.set_min_x(0.);
    assert!(r != r2);
    r2.set_min_x(r.min_position()[0]);
    assert!(r == r2);
    assert!(r2 == r);
    r2.set_max_y(0.);
    assert!(r != r2);
    r2.set_max_y(r.max_position()[1]);
    assert!(r == r2);
}

// L140 encloses(position): nine literal points, half-open upper bound
#[test]
fn encloses_position() {
    let r2 = DRange::new(p1(), p2());
    let cases = [
        ((0., 0.), true),
        ((-3., -3.), false),
        ((-3., 0.), false),
        ((0., -3.), false),
        ((-3., 5.), false),
        ((0., 5.), false),
        ((5., 5.), false),
        ((5., 0.), false),
        ((5., -3.), false),
    ];
    for ((x, y), expected) in cases {
        assert_eq!(r2.encloses(&DPosition::xy(x, y)), expected, "({x},{y})");
    }
    assert!(r2.encloses(&p1()));
    assert!(!r2.encloses(&p2()));
    assert!(!r2.encloses(&DPosition::xy(3., 0.)));
}

/// The 21 `r3` configurations of L172-L264 with their `intersects` outcome;
/// `isIntersected` (L267-L359) is `!= Disjoint` for each.
fn intersection_cases() -> Vec<(DRange2, Kind)> {
    let r2 = DRange::new(p1(), p2());
    let mut r3 = r2;
    let mut cases = Vec::new();
    cases.push((r3, Kind::Inside));
    r3.set_max_x(10.);
    cases.push((r3, Kind::Intersects));
    r3.set_max(*r2.max_position() + one());
    cases.push((r3, Kind::Intersects));
    r3.set_min(*r2.max_position() + one());
    r3.set_max(*r2.max_position() + two());
    cases.push((r3, Kind::Disjoint));
    r3.set_min(*r2.min_position());
    r3.set_min_x(10.);
    r3.set_max(*r3.min_position() + one());
    cases.push((r3, Kind::Disjoint));
    r3.set_min_x(-10.);
    r3.set_min_y(-10.);
    r3.set_max(*r3.min_position() + one());
    cases.push((r3, Kind::Disjoint));
    let literal = [
        ((-10., -10., 0., -9.), Kind::Disjoint),
        ((-10., -10., 10., -9.), Kind::Disjoint),
        ((-10., 0., -9., 1.), Kind::Disjoint),
    ];
    for ((min_x, min_y, max_x, max_y), kind) in literal {
        r3.set_min_x(min_x);
        r3.set_min_y(min_y);
        r3.set_max_x(max_x);
        r3.set_max_y(max_y);
        cases.push((r3, kind));
    }
    r3.set_min_x(-10.);
    r3.set_min_y(10.);
    r3.set_max(*r3.min_position() + one());
    cases.push((r3, Kind::Disjoint));
    let literal = [
        ((-10., 0., -9., 10.), Kind::Disjoint),
        ((9., 0., 10., 10.), Kind::Disjoint),
        ((9., 0., 10., 10.), Kind::Disjoint),
        ((9., -5., 10., 0.), Kind::Disjoint),
        ((9., -5., 10., 5.), Kind::Disjoint),
        ((-5., -5., 0., 0.), Kind::Intersects),
        ((-5., -5., 5., 0.), Kind::Intersects),
        ((-5., -5., 5., 5.), Kind::Intersects),
        ((0., -5., 0., 0.), Kind::Intersects),
        ((0., -5., 5., 0.), Kind::Intersects),
        ((0., -5., 5., 5.), Kind::Intersects),
    ];
    for ((min_x, min_y, max_x, max_y), kind) in literal {
        r3.set_min_x(min_x);
        r3.set_min_y(min_y);
        r3.set_max_x(max_x);
        r3.set_max_y(max_y);
        cases.push((r3, kind));
    }
    assert_eq!(cases.len(), 21);
    cases
}

// L172 intersects(): 21 literal configurations
#[test]
fn intersects_classification() {
    let r2 = DRange::new(p1(), p2());
    for (index, (r3, kind)) in intersection_cases().into_iter().enumerate() {
        assert_eq!(r2.intersects(&r3), kind, "case {index}: {r3}");
    }
    assert_eq!(
        DRange::xy(0., 0., 1., 1.).intersects(&DRange::xy(0., 0., 1., 1.)),
        Kind::Inside
    );
    assert_eq!(
        DRange::xy(0., 0., 1., 1.).intersects(&DRange::xy(-1., -1., 2., 2.)),
        Kind::Intersects
    );
}

// L267 isIntersected(): the same 21 configurations as booleans
#[test]
fn is_intersected_is_not_disjoint() {
    let r2 = DRange::new(p1(), p2());
    let expected = [
        true, true, true, false, false, false, false, false, false, false, false, false, false,
        false, false, true, true, true, true, true, true,
    ];
    for (index, ((r3, kind), expected)) in
        intersection_cases().into_iter().zip(expected).enumerate()
    {
        assert_eq!(r2.is_intersected(&r3), expected, "case {index}");
        assert_eq!(kind != Kind::Disjoint, expected, "case {index}");
    }
}

// L362 united()
#[test]
fn united_is_the_bounding_range() {
    let r2 = DRange::new(p1(), p2());
    let mut r3 = r2;
    assert!(r2 == r2.united(&r3));
    assert!(r3 == r2.united(&r3));
    assert!(r2 == r3.united(&r2));
    assert!(r3 == r3.united(&r2));
    r3.set_min(*r2.max_position() + one());
    r3.set_max(*r2.max_position() + two());
    let mut r4 = DRange2::default();
    r4.set_min(*r2.min_position());
    r4.set_max(*r3.max_position());
    assert!(r2.united(&r3) == r4);
    assert!(r3.united(&r2) == r4);
    assert_eq!(corners(&r4), [-1., -2., 5., 6.]);
    assert_eq!(r2.united(&DRange2::empty()), r2);
    // Source quirk: the sentinel corners are min > max, so uniting two empty
    // ranges hands (MAX, MIN) to setMinMax, whose normalization swaps them
    // into the all-encompassing range instead of keeping the empty one.
    assert_eq!(
        DRange2::empty().united(&DRange2::empty()),
        DRange::new(DPosition::min_negative(), DPosition::max_positive())
    );
}

// L379 encloses(x, y)
#[test]
fn encloses_xy() {
    let r2 = DRange::new(p1(), p2());
    assert!(r2.encloses_xy(0., 0.));
    assert!(!r2.encloses_xy(-3., -3.));
    assert!(!r2.encloses_xy(-3., 0.));
    assert!(!r2.encloses_xy(0., -3.));
    assert!(!r2.encloses_xy(-3., 5.));
    assert!(!r2.encloses_xy(0., 5.));
    assert!(!r2.encloses_xy(5., 5.));
    assert!(!r2.encloses_xy(5., 0.));
    assert!(!r2.encloses_xy(5., -3.));
}

// L393 extend(double factor)
#[test]
fn extend_by_factor_keeps_the_center() {
    let mut r = r();
    assert!(matches!(
        r.extend_by_factor(-0.01),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(corners(&r), [-1., -2., 3., 4.]);
    assert!(matches!(
        r.extend_by_factor(f64::NAN),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        r.extend_by_factor(f64::INFINITY),
        Err(Error::InvalidValue(_))
    ));
    let other = *r.extend_by_factor(2.).unwrap();
    close(r.min_position()[0], -3.);
    close(r.max_position()[0], 5.);
    close(r.min_position()[1], -5.);
    close(r.max_position()[1], 7.);
    close(other.min_position()[0], -3.);
    close(other.max_position()[0], 5.);
    let mut unit = DRange1::new(DPosition1::filled(0.), DPosition1::filled(100.));
    unit.extend_by_factor(1.01).unwrap();
    close(unit.min_position()[0], -0.5);
    close(unit.max_position()[0], 100.5);
    unit.extend_by_factor(0.).unwrap();
    close(unit.min_position()[0], 50.);
    close(unit.max_position()[0], 50.);
}

// L411 extend(PositionType addition), including shrinking to a point
#[test]
fn extend_by_addition_and_shrink_to_center() {
    let mut r = r();
    let other = *r.extend_by(DPosition::xy(2., 3.));
    close(r.min_position()[0], -2.);
    close(r.max_position()[0], 4.);
    close(r.min_position()[1], -3.5);
    close(r.max_position()[1], 5.5);
    close(other.min_position()[0], -2.);
    close(other.max_position()[0], 4.);
    r.extend_by(DPosition::xy(-200., 0.));
    close(r.min_position()[0], 1.);
    close(r.max_position()[0], 1.);
    close(r.min_position()[1], -3.5);
    close(r.max_position()[1], 5.5);
}

// L435 ensureMinSpan()
#[test]
fn ensure_min_span_widens_only_narrow_dimensions() {
    let mut r = DRange::xy(-0.1, 10., 0.1, 20.);
    r.ensure_min_span(DPosition::xy(1., 3.));
    close(r.min_position()[0], -0.5);
    close(r.max_position()[0], 0.5);
    close(r.min_position()[1], 10.);
    close(r.max_position()[1], 20.);
}

// L444 swapDimensions()
#[test]
fn swap_dimensions() {
    let mut r = r();
    r.swap_dimensions();
    close(r.min_position()[0], -2.);
    close(r.max_position()[0], 4.);
    close(r.min_position()[1], -1.);
    close(r.max_position()[1], 3.);
}

// L459 pullIn(): the in/out parameter becomes the return value
#[test]
fn pull_in_clamps_to_the_closed_corners() {
    let r = DRange::new(DPosition::xy(1., 2.), DPosition::xy(3., 4.));
    let p_out_left = r.pull_in(DPosition::xy(0., 0.));
    close(p_out_left.x(), 1.);
    close(p_out_left.y(), 2.);
    let p_out_right = r.pull_in(DPosition::xy(5., 5.));
    close(p_out_right.x(), 3.);
    close(p_out_right.y(), 4.);
    let p_in = r.pull_in(DPosition::xy(2., 3.));
    close(p_in.x(), 2.);
    close(p_in.y(), 3.);
    // Transcribed std::max(min, std::min(point, max)) sends NaN to the minimum.
    assert_eq!(r.pull_in(DPosition::xy(f64::NAN, 4.)).coordinates, [1., 4.]);
}

// L480 [EXTRA] std::hash<DRange<D>>. Without Eq the unordered containers are
// represented by equality/hash agreement, as the crate's other value types do.
#[test]
fn hash_agrees_with_equality() {
    let r1 = DRange::xy(1., 2., 3., 4.);
    let r2 = DRange::xy(1., 2., 3., 4.);
    let r3 = DRange::xy(0., 0., 5., 5.);
    assert!(r1 == r2);
    assert_eq!(hash(&r1), hash(&r2));
    assert!(r1 != r3);
    assert_ne!(hash(&r1), hash(&r3));
    let mut r1d_a = DRange1::default();
    let mut r1d_b = DRange1::default();
    r1d_a.set_min(DPosition1::filled(1.));
    r1d_a.set_max(DPosition1::filled(5.));
    r1d_b.set_min(DPosition1::filled(1.));
    r1d_b.set_max(DPosition1::filled(5.));
    assert!(r1d_a == r1d_b);
    assert_eq!(hash(&r1d_a), hash(&r1d_b));
    assert_eq!(hash(&r1), hash(r1.base()));
    assert_eq!(
        hash(&DRange::xy(0., -0., 1., 1.)),
        hash(&DRange::xy(-0., 0., 1., 1.))
    );
}

// Native: delegated base accessors and stream output.
#[test]
fn delegated_base_accessors_and_display() {
    let mut range = r();
    assert_eq!(range.width(), 4.);
    assert_eq!(range.height(), 6.);
    assert_eq!(range.center(), DPosition::xy(1., 1.));
    assert_eq!(range.diagonal(), DPosition::xy(4., 6.));
    assert_eq!(
        [range.min_x(), range.min_y(), range.max_x(), range.max_y()],
        [-1., -2., 3., 4.]
    );
    range.set_min_y(-1.);
    range.set_max_x(2.);
    assert_eq!(corners(&range), [-1., -1., 2., 4.]);
    range.set_min_max(p2(), p1());
    assert_eq!(corners(&range), [-1., -2., 3., 4.]);
    range
        .set_dim_min_max(1, &DIntervalBase::<1>::zero())
        .unwrap();
    assert_eq!(corners(&range), [-1., 0., 3., 0.]);
    assert!(!range.is_empty_dim(1).unwrap());
    let mut three = DRange::<3>::default();
    three.assign(range.base());
    assert_eq!(three.min_position().coordinates, [-1., 0., f64::MAX]);
    range.clear();
    assert!(range.is_empty());
    assert!(range.is_empty_dim(0).unwrap());
    let shifted = r() + one();
    assert_eq!(corners(&shifted), [0., -1., 4., 5.]);
    let mut back = shifted - one();
    assert_eq!(back, r());
    back += two();
    back -= two();
    assert_eq!(back, r());
    assert_eq!(
        r().to_string(),
        "--DRANGE BEGIN--\nMIN --> -1 -2\nMAX --> 3 4\n--DRANGE END--\n"
    );
    assert_eq!(
        format!("{:.1}", DRange1::zero()),
        "--DRANGE BEGIN--\nMIN --> 0.0\nMAX --> 0.0\n--DRANGE END--\n"
    );
}
