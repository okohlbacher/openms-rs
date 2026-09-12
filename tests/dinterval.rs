// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Port of `DIntervalBase_test.cpp` (core SDK bc9cc12): one test per
//! `START_SECTION`, cited by source line. Literals are transcribed (tier 3).

use openms::Error;
use openms::data_structures::{DIntervalBase, DIntervalBase1, DIntervalBase2, DPosition};

type I2 = DIntervalBase2;
type I2Pos = DPosition<2>;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1.),
        "{actual:.17} != {expected:.17}"
    );
}
// Source globals (L55-60): p1 = (5, 17.5), p2 = (65, -57.5).
fn p1() -> I2Pos {
    I2Pos::xy(5., 17.5)
}
fn p2() -> I2Pos {
    I2Pos::xy(65., -57.5)
}
fn di() -> I2 {
    I2::new(I2Pos::xy(1., 2.), I2Pos::xy(3., 4.))
}

// L30 DIntervalBase() 1D
#[test]
fn default_constructor_1d_is_empty() {
    let i = DIntervalBase1::default();
    assert!(i.is_empty());
    assert_eq!(i.min_position()[0], f64::MAX);
    assert_eq!(i.max_position()[0], f64::MIN);
}

// L35 ~DIntervalBase() 1D
#[test]
fn destructor_1d_is_ordinary_drop() {
    let i = Box::new(DIntervalBase1::default());
    assert_eq!(*i, DIntervalBase1::empty());
    drop(i);
}

// L43 [EXTRA] DIntervalBase() 2D
#[test]
fn default_constructor_2d_is_empty() {
    let i = I2::default();
    assert_eq!(i, I2::empty());
    assert_eq!(I2::DIMENSION, 2);
}

// L48 [EXTRA] ~DIntervalBase() 2D
#[test]
fn destructor_2d_is_ordinary_drop() {
    let i = Box::new(I2::default());
    assert!(i.is_empty());
    drop(i);
}

// L62 operator+
#[test]
fn add_translates_both_corners() {
    let r = di() + I2Pos::xy(1., 0.5);
    close(r.min_x(), 2.);
    close(r.min_y(), 2.5);
    close(r.max_x(), 4.);
    close(r.max_y(), 4.5);
}

// L73 operator+=
#[test]
fn add_assign_translates_both_corners() {
    let mut di = di();
    di += I2Pos::xy(1., 0.5);
    let r = di;
    assert_eq!(r, di);
    close(r.min_x(), 2.);
    close(r.min_y(), 2.5);
    close(r.max_x(), 4.);
    close(r.max_y(), 4.5);
}

// L85 operator-
#[test]
fn sub_translates_both_corners() {
    let r = di() - I2Pos::xy(1., 0.5);
    close(r.min_x(), 0.);
    close(r.min_y(), 1.5);
    close(r.max_x(), 2.);
    close(r.max_y(), 3.5);
}

// L96 operator-=
#[test]
fn sub_assign_translates_both_corners() {
    let mut di = di();
    di -= I2Pos::xy(1., 0.5);
    let r = di;
    assert_eq!(r, di);
    close(r.min_x(), 0.);
    close(r.min_y(), 1.5);
    close(r.max_x(), 2.);
    close(r.max_y(), 3.5);
}

// L109 maxPosition()
#[test]
fn max_position_of_statics() {
    assert!(*I2::empty().max_position() == I2Pos::min_negative());
    assert!(*I2::zero().max_position() == I2Pos::zero());
}

// L114 minPosition()
#[test]
fn min_position_of_statics() {
    assert!(*I2::empty().min_position() == I2Pos::max_positive());
    assert!(*I2::zero().min_position() == I2Pos::zero());
}

// L119 setMinMax normalizes per dimension
#[test]
fn set_min_max_normalizes_each_dimension() {
    let mut tmp = I2::empty();
    tmp.set_min_max(p1(), p2());
    close(tmp.min_position()[0], 5.);
    close(tmp.min_position()[1], -57.5);
    close(tmp.max_position()[0], 65.);
    close(tmp.max_position()[1], 17.5);
}

// L128 setMin adjusts the maximum where needed
#[test]
fn set_min_raises_maximum() {
    let mut tmp = I2::empty();
    tmp.set_min(p1());
    assert_eq!(*tmp.min_position(), p1());
    assert_eq!(*tmp.max_position(), p1());
    tmp.set_min(p2());
    close(tmp.min_position()[0], 65.);
    close(tmp.min_position()[1], -57.5);
    close(tmp.max_position()[0], 65.);
    close(tmp.max_position()[1], 17.5);
}

// L140 setMax adjusts the minimum where needed
#[test]
fn set_max_lowers_minimum() {
    let mut tmp = I2::empty();
    tmp.set_max(p1());
    assert_eq!(*tmp.min_position(), p1());
    assert_eq!(*tmp.max_position(), p1());
    tmp.set_max(p2());
    close(tmp.min_position()[0], 5.);
    close(tmp.min_position()[1], -57.5);
    close(tmp.max_position()[0], 65.);
    close(tmp.max_position()[1], -57.5);
}

// L152 setDimMinMax
#[test]
fn set_dim_min_max_sets_one_dimension() {
    let mut tmp = I2::empty();
    let mut min_p = *tmp.min_position();
    let mut max_p = *tmp.max_position();
    tmp.set_dim_min_max(0, &DIntervalBase1::new([1.].into(), [1.1].into()))
        .unwrap();
    min_p.set_x(1.);
    max_p.set_x(1.1);
    assert_eq!(*tmp.min_position(), min_p);
    assert_eq!(*tmp.max_position(), max_p);
    assert!(matches!(
        tmp.set_dim_min_max(2, &DIntervalBase1::zero()),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(*tmp.min_position(), min_p);
}

// L164 operator==
#[test]
fn equality() {
    let mut tmp = I2::default();
    let same = tmp;
    assert!(tmp == same);
    assert!(tmp == I2::empty());
    tmp.set_max(p1());
    assert!(tmp != I2::empty());
}

// L173 operator!=
#[test]
fn inequality() {
    let mut tmp = I2::default();
    let same = tmp;
    assert!(!(tmp != same));
    assert!(!(tmp != I2::empty()));
    tmp.set_max(p1());
    assert!(tmp != I2::empty());
}

// L182 copy constructor
#[test]
fn copy_constructor() {
    let tmp = I2::new(p1(), p2());
    let tmp2 = tmp;
    assert!(tmp == tmp2);
}

// L188 DIntervalBase(minimum, maximum)
#[test]
fn corner_constructor_round_trips() {
    let tmp = I2::new(p1(), p2());
    let tmp2 = I2::new(*tmp.min_position(), *tmp.max_position());
    assert!(tmp == tmp2);
    assert_eq!(tmp.min_position().coordinates, [5., -57.5]);
}

// L194 operator=
#[test]
fn assignment() {
    let mut tmp = I2::new(p1(), p2());
    let mut tmp2 = I2::default();
    assert!(tmp != tmp2);
    tmp2 = tmp;
    assert_eq!(tmp2.min_position()[0], 5.);
    assert!(tmp == tmp2);
    tmp = I2::empty();
    tmp2 = tmp;
    assert!(tmp == tmp2);
    assert!(tmp == I2::empty());
}

// L205 clear()
#[test]
fn clear_restores_the_empty_sentinel() {
    let mut tmp = I2::default();
    assert!(tmp == I2::empty());
    tmp.set_max(p1());
    assert!(tmp != I2::empty());
    tmp.clear();
    assert!(tmp == I2::empty());
    assert!(*tmp.max_position() == I2Pos::min_negative());
    assert!(*tmp.min_position() == I2Pos::max_positive());
}

// L216 isEmpty(): min == max is not empty
#[test]
fn is_empty_only_for_the_sentinel() {
    let mut tmp = I2::default();
    assert!(tmp.is_empty());
    tmp.set_max(p1());
    assert!(!tmp.is_empty());
    tmp.clear();
    assert!(tmp.is_empty());
    tmp.set_dim_min_max(1, &DIntervalBase1::new([2.].into(), [2.].into()))
        .unwrap();
    assert!(!tmp.is_empty());
}

// L227 isEmpty(dim)
#[test]
fn is_empty_dim_checks_one_dimension() {
    let mut tmp = I2::default();
    assert!(tmp.is_empty_dim(0).unwrap());
    assert!(tmp.is_empty_dim(1).unwrap());
    tmp.set_max(p1());
    assert!(!tmp.is_empty_dim(0).unwrap());
    assert!(!tmp.is_empty_dim(1).unwrap());
    tmp.clear();
    assert!(tmp.is_empty_dim(0).unwrap());
    assert!(tmp.is_empty_dim(1).unwrap());
    tmp.set_dim_min_max(1, &DIntervalBase1::new([2.].into(), [2.].into()))
        .unwrap();
    assert!(tmp.is_empty_dim(0).unwrap());
    assert!(!tmp.is_empty_dim(1).unwrap());
    assert!(matches!(tmp.is_empty_dim(2), Err(Error::InvalidValue(_))));
}

// L242 center()
#[test]
fn center_is_the_midpoint() {
    let tmp = I2::new(p1(), p2());
    let pos = tmp.center();
    close(pos[0], 35.);
    close(pos[1], -20.);
    assert_eq!(I2::empty().center(), I2Pos::zero());
}

// L249 diagonal()
#[test]
fn diagonal_is_max_minus_min() {
    let tmp = I2::new(p1(), p2());
    let pos = tmp.diagonal();
    close(pos[0], 60.);
    close(pos[1], 75.);
}

// L256 width()
#[test]
fn width() {
    close(I2::new(p1(), p2()).width(), 60.);
    close(DIntervalBase1::new([2.].into(), [5.].into()).width(), 3.);
}

// L261 height()
#[test]
fn height() {
    close(I2::new(p1(), p2()).height(), 75.);
}

// L266 maxX()
#[test]
fn max_x() {
    close(I2::new(p1(), p2()).max_x(), 65.);
}

// L271 maxY()
#[test]
fn max_y() {
    close(I2::new(p1(), p2()).max_y(), 17.5);
}

// L276 minX()
#[test]
fn min_x() {
    close(I2::new(p1(), p2()).min_x(), 5.);
}

// L281 minY()
#[test]
fn min_y() {
    close(I2::new(p1(), p2()).min_y(), -57.5);
}

// L286 setMinX()
#[test]
fn set_min_x() {
    let mut tmp = I2::new(p1(), p2());
    tmp.set_min_x(57.67);
    close(tmp.min_x(), 57.67);
    close(tmp.max_x(), 65.);
    tmp.set_min_x(70.);
    close(tmp.max_x(), 70.);
}

// L292 setMaxX()
#[test]
fn set_max_x() {
    let mut tmp = I2::new(p1(), p2());
    tmp.set_max_x(57.67);
    close(tmp.max_x(), 57.67);
    tmp.set_max_x(1.);
    close(tmp.min_x(), 1.);
}

// L298 setMinY()
#[test]
fn set_min_y() {
    let mut tmp = I2::new(p1(), p2());
    tmp.set_min_y(57.67);
    close(tmp.min_y(), 57.67);
    close(tmp.max_y(), 57.67);
}

// L304 setMaxY()
#[test]
fn set_max_y() {
    let mut tmp = I2::new(p1(), p2());
    tmp.set_max_y(57.67);
    close(tmp.max_y(), 57.67);
    close(tmp.min_y(), -57.5);
}

// L312 assign<D2>
#[test]
fn assign_copies_the_shared_dimensions() {
    let i2 = I2::new(p1(), p2());
    let mut tmp = DIntervalBase::<3>::default();
    tmp.assign(&i2);
    close(tmp.min_position()[0], 5.);
    close(tmp.min_position()[1], -57.5);
    close(tmp.max_position()[0], 65.);
    close(tmp.max_position()[1], 17.5);
    assert_eq!(tmp.min_position()[2], f64::MAX);
    let mut tmp2 = DIntervalBase1::default();
    tmp2.assign(&i2);
    close(tmp2.min_position()[0], 5.);
    close(tmp2.max_position()[0], 65.);
}

// Native: stream output, per-dimension normalization of NaN, hashing.
#[test]
fn display_and_hash() {
    let tmp = I2::new(p1(), p2());
    assert_eq!(
        tmp.to_string(),
        "--DIntervalBase BEGIN--\nMIN --> 5 -57.5\nMAX --> 65 17.5\n--DIntervalBase END--\n"
    );
    assert_eq!(
        format!("{:.1}", DIntervalBase1::zero()),
        "--DIntervalBase BEGIN--\nMIN --> 0.0\nMAX --> 0.0\n--DIntervalBase END--\n"
    );
    let hash = |value: &I2| {
        use std::hash::{DefaultHasher, Hash, Hasher};
        let mut h = DefaultHasher::new();
        value.hash(&mut h);
        h.finish()
    };
    assert_eq!(hash(&tmp), hash(&I2::new(p2(), p1())));
    assert_ne!(hash(&tmp), hash(&I2::zero()));
}
