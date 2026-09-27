// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `DBoundingBox_test.cpp` (core `bc9cc12`): one test per
//! `START_SECTION`, literals transcribed (tier 3), plus the source
//! `std::hash<DBoundingBox<D>>` value recomputed in Python from the header's
//! `hash_combine`/`hash_float` recurrence over explicit little-endian bytes
//! (tier 4).

use openms::data_structures::{
    DBoundingBox, DBoundingBox1, DBoundingBox2, DIntervalBase2, DPosition, DPosition1, DPosition2,
};
use std::collections::{HashMap, HashSet};

type BB1 = DBoundingBox1;
type BB2 = DBoundingBox2;

#[test]
fn default_constructors_and_dimension() {
    // DBoundingBox() / ~DBoundingBox() and the [EXTRA] 2D pair.
    let one = BB1::default();
    let two = BB2::default();
    assert_eq!(BB1::DIMENSION, 1);
    assert_eq!(BB2::DIMENSION, 2);
    assert_eq!(DBoundingBox::<3>::DIMENSION, 3);
    assert!(one.is_empty());
    assert!(two.is_empty());
    assert_eq!(two.min_position()[0], f64::MAX);
    assert_eq!(two.max_position()[1], f64::MIN);
}

#[test]
fn constructor_from_two_positions() {
    let bb = BB1::new(DPosition1::filled(2.0), DPosition1::filled(5.0));
    assert_eq!(bb.min_position()[0], 2.0);
    assert_eq!(bb.max_position()[0], 5.0);
}

#[test]
fn copy_constructor_and_assignment() {
    let bb = BB2::new(DPosition2::xy(1.0, 2.0), DPosition2::xy(3.0, 4.0));
    let copy = bb;
    let mut assigned = BB2::default();
    assert_ne!(assigned, bb);
    assigned = copy;
    for other in [copy, assigned] {
        assert_eq!(other.min_position()[0], bb.min_position()[0]);
        assert_eq!(other.min_position()[1], bb.min_position()[1]);
        assert_eq!(other.max_position()[0], bb.max_position()[0]);
        assert_eq!(other.max_position()[1], bb.max_position()[1]);
    }
}

#[test]
fn assignment_from_base() {
    let base = DIntervalBase2::new(DPosition2::xy(1.0, 2.0), DPosition2::xy(3.0, 4.0));
    let copy = BB2::from(base);
    assert_eq!(copy.min_position()[0], base.min_position()[0]);
    assert_eq!(copy.min_position()[1], base.min_position()[1]);
    assert_eq!(copy.max_position()[0], base.max_position()[0]);
    assert_eq!(copy.max_position()[1], base.max_position()[1]);
    assert_eq!(DIntervalBase2::from(copy), base);
    assert_eq!(copy.as_interval(), &base);
}

#[test]
fn is_empty() {
    let mut bb = BB2::from(DIntervalBase2::empty());
    assert!(bb.is_empty());
    bb = BB2::from(DIntervalBase2::zero());
    assert!(bb.is_empty());
    bb = BB2::new(DPosition2::xy(1.0, 2.0), DPosition2::xy(3.0, 4.0));
    assert!(!bb.is_empty());
    // The derived isEmpty hides the base's sentinel test: a zero box is empty
    // for the box but not for the interval.
    let zero = BB2::from(DIntervalBase2::zero());
    assert!(zero.is_empty());
    assert!(!zero.as_interval().is_empty());
    // One degenerate dimension is enough.
    let flat = BB2::new(DPosition2::xy(1.0, 2.0), DPosition2::xy(3.0, 2.0));
    assert!(flat.is_empty());
}

#[test]
fn enlarge_with_a_position() {
    let mut bb = BB2::default();
    assert!(!bb.encloses_xy(11.0, 13.0));
    assert!(!bb.encloses_xy(10.0, 1.0));
    bb.enlarge(DPosition2::xy(11.0, 13.0));
    assert!(bb.encloses_xy(11.0, 13.0));
    assert!(!bb.encloses_xy(10.0, 1.0));
    bb.enlarge(DPosition2::xy(9.0, 0.0));
    assert!(bb.encloses_xy(11.0, 13.0));
    assert!(bb.encloses_xy(10.0, 1.0));
}

#[test]
fn enlarge_with_two_coordinates() {
    let mut bb = BB2::default();
    assert!(!bb.encloses_xy(11.0, 13.0));
    assert!(!bb.encloses_xy(10.0, 1.0));
    bb.enlarge_xy(11.0, 13.0);
    assert!(bb.encloses_xy(11.0, 13.0));
    assert!(!bb.encloses_xy(10.0, 1.0));
    bb.enlarge_xy(9.0, 0.0);
    assert!(bb.encloses_xy(11.0, 13.0));
    assert!(bb.encloses_xy(10.0, 1.0));
}

#[test]
fn equality_with_a_box_and_with_the_base() {
    let mut bb = BB2::default();
    bb.enlarge_xy(9.0, 0.0);
    let copy = bb;
    assert!(bb == copy);
    let base: DIntervalBase2 = bb.into();
    assert!(bb == base);
}

fn box_100_200_300_400() -> BB2 {
    let mut tmp = BB2::default();
    tmp.set_min_x(100.0);
    tmp.set_min_y(200.0);
    tmp.set_max_x(300.0);
    tmp.set_max_y(400.0);
    tmp
}

const ENCLOSES: [(f64, f64, bool); 20] = [
    (10.0, 200.0, false),
    (100.0, 200.0, true),
    (200.0, 200.0, true),
    (300.0, 200.0, true),
    (310.0, 200.0, false),
    (10.0, 400.0, false),
    (100.0, 400.0, true),
    (200.0, 400.0, true),
    (300.0, 400.0, true),
    (310.0, 400.0, false),
    (200.0, 190.0, false),
    (200.0, 200.0, true),
    (200.0, 300.0, true),
    (200.0, 400.0, true),
    (200.0, 410.0, false),
    (0.0, 0.0, false),
    (100.0, 200.0, true),
    (300.0, 200.0, true),
    (100.0, 400.0, true),
    (300.0, 400.0, true),
];

#[test]
fn encloses_two_coordinates() {
    let tmp = box_100_200_300_400();
    for (x, y, expected) in ENCLOSES {
        assert_eq!(tmp.encloses_xy(x, y), expected, "({x}, {y})");
    }
}

#[test]
fn encloses_a_position() {
    let tmp = box_100_200_300_400();
    for (x, y, expected) in ENCLOSES {
        assert_eq!(tmp.encloses(&DPosition2::xy(x, y)), expected, "({x}, {y})");
    }
}

#[test]
fn encloses_a_nan_coordinate_as_the_source_does() {
    // The source rejects only on `<` or `>`, both false for NaN.
    let tmp = box_100_200_300_400();
    assert!(tmp.encloses_xy(f64::NAN, 300.0));
    assert!(!tmp.encloses_xy(f64::NAN, 500.0));
}

#[test]
fn intersects() {
    let p1 = DPosition2::xy(-1.0, -2.0);
    let p2 = DPosition2::xy(3.0, 4.0);
    let one = DPosition2::xy(1.0, 1.0);
    let two = DPosition2::xy(2.0, 2.0);
    let mut r2 = BB2::default();
    r2.set_min(p1);
    r2.set_max(p2);
    let mut r3 = r2;
    assert!(r2.intersects(&r3));
    r3.set_max_x(10.0);
    assert!(r2.intersects(&r3));
    r3.set_max(*r2.max_position() + one);
    assert!(r2.intersects(&r3));
    r3.set_min(*r2.max_position() + one);
    r3.set_max(*r2.max_position() + two);
    assert!(!r2.intersects(&r3));
    r3.set_min(*r2.min_position());
    r3.set_min_x(10.0);
    r3.set_max(*r3.min_position() + one);
    assert!(!r2.intersects(&r3));
    r3.set_min_x(-10.0);
    r3.set_min_y(-10.0);
    r3.set_max(*r3.min_position() + one);
    assert!(!r2.intersects(&r3));

    // From here the source sets all four coordinates, in the order setMinX,
    // setMinY, setMaxX, setMaxY, except in the one `setMax` case below.
    let check = |r3: &mut BB2, corners: [f64; 4], expected: bool| {
        r3.set_min_x(corners[0]);
        r3.set_min_y(corners[1]);
        r3.set_max_x(corners[2]);
        r3.set_max_y(corners[3]);
        assert_eq!(r2.intersects(r3), expected, "{corners:?}");
    };
    check(&mut r3, [-10.0, -10.0, 0.0, -9.0], false);
    check(&mut r3, [-10.0, -10.0, 10.0, -9.0], false);
    check(&mut r3, [-10.0, 0.0, -9.0, 1.0], false);
    r3.set_min_x(-10.0);
    r3.set_min_y(10.0);
    r3.set_max(*r3.min_position() + one);
    assert!(!r2.intersects(&r3));
    check(&mut r3, [-10.0, 0.0, -9.0, 10.0], false);
    check(&mut r3, [9.0, 0.0, 10.0, 10.0], false);
    check(&mut r3, [9.0, 0.0, 10.0, 10.0], false);
    check(&mut r3, [9.0, -5.0, 10.0, 0.0], false);
    check(&mut r3, [9.0, -5.0, 10.0, 5.0], false);
    check(&mut r3, [-5.0, -5.0, 0.0, 0.0], true);
    check(&mut r3, [-5.0, -5.0, 5.0, 0.0], true);
    check(&mut r3, [-5.0, -5.0, 5.0, 5.0], true);
    check(&mut r3, [0.0, -5.0, 0.0, 0.0], true);
    check(&mut r3, [0.0, -5.0, 5.0, 0.0], true);
    check(&mut r3, [0.0, -5.0, 5.0, 5.0], true);
}

#[test]
fn output_operator() {
    let bb = BB1::new(DPosition1::filled(2.0), DPosition1::filled(5.0));
    assert_eq!(
        bb.to_string(),
        "--DBOUNDINGBOX BEGIN--\nMIN --> 2.0\nMAX --> 5.0\n--DBOUNDINGBOX END--\n"
    );
    // Two dimensions: coordinates separated by one space, each through
    // precisionWrapper (toStr(double, true)).
    let bb = BB2::new(DPosition2::xy(1.5, -0.25), DPosition2::xy(1.0e6, 3.0));
    assert_eq!(
        bb.to_string(),
        "--DBOUNDINGBOX BEGIN--\nMIN --> 1.5 -0.25\nMAX --> 1.0e06 3.0\n--DBOUNDINGBOX END--\n"
    );
}

#[test]
fn std_hash_contract() {
    // [EXTRA] std::hash<DBoundingBox<D>>
    let bb1 = BB2::new(DPosition2::xy(1.0, 2.0), DPosition2::xy(3.0, 4.0));
    let bb2 = BB2::new(DPosition2::xy(1.0, 2.0), DPosition2::xy(3.0, 4.0));
    let bb3 = BB2::new(DPosition2::xy(5.0, 6.0), DPosition2::xy(7.0, 8.0));
    assert!(bb1 == bb2);
    assert_eq!(bb1.source_hash64(), bb2.source_hash64());
    assert!(bb1 != bb3);
    assert_ne!(bb1.source_hash64(), bb3.source_hash64());

    let a = BB1::new(DPosition1::filled(1.0), DPosition1::filled(5.0));
    let b = BB1::new(DPosition1::filled(1.0), DPosition1::filled(5.0));
    let c = BB1::new(DPosition1::filled(2.0), DPosition1::filled(6.0));
    assert!(a == b);
    assert_eq!(a.source_hash64(), b.source_hash64());
    assert_ne!(a.source_hash64(), c.source_hash64());

    let mut set = HashSet::new();
    set.insert(bb1);
    set.insert(bb2);
    set.insert(bb3);
    assert_eq!(set.len(), 2);
    assert!(set.contains(&bb1));
    assert!(set.contains(&bb3));

    let mut map = HashMap::new();
    map.insert(bb1, "first");
    map.insert(bb2, "second");
    map.insert(bb3, "third");
    assert_eq!(map.len(), 2);
    assert_eq!(map[&bb1], "second");
    assert_eq!(map[&bb3], "third");

    let default1 = BB2::default();
    let default2 = BB2::default();
    assert!(default1 == default2);
    assert_eq!(default1.source_hash64(), default2.source_hash64());
}

#[test]
fn source_hash_value_matches_the_header_recurrence() {
    // Python: seed = 0; for v in (1.0, 5.0): seed = combine(seed, fnv(pack('<d', v))).
    let bb = BB1::new(DPosition1::filled(1.0), DPosition1::filled(5.0));
    assert_eq!(bb.source_hash64(), SOURCE_HASH_1_5);
    // Negative zero hashes as positive zero.
    let negative = BB1::new(DPosition1::filled(-0.0), DPosition1::filled(1.0));
    let positive = BB1::new(DPosition1::filled(0.0), DPosition1::filled(1.0));
    assert_eq!(negative, positive);
    assert_eq!(negative.source_hash64(), positive.source_hash64());
}

const SOURCE_HASH_1_5: u64 = 0xdbc6_4f7b_7b96_64d4;

#[test]
fn a_box_of_any_dimension() {
    let mut cube = DBoundingBox::<3>::default();
    cube.enlarge(DPosition::xyz(1.0, 2.0, 3.0));
    cube.enlarge(DPosition::xyz(-1.0, 5.0, 3.0));
    assert!(cube.is_empty()); // the third dimension has zero width
    cube.enlarge(DPosition::xyz(0.0, 0.0, 4.0));
    assert!(!cube.is_empty());
    assert!(cube.encloses(&DPosition::xyz(0.0, 3.0, 3.5)));
    assert!(!cube.encloses(&DPosition::xyz(0.0, 3.0, 4.5)));
}
