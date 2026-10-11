// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `DATASTRUCTURES/ConvexHull2D.h`: `getMapPoints` and `std::hash<ConvexHull2D>`.
//!
//! The rest of the class test is ported in `tests/features.rs`. Here: the
//! `[EXTRA] std::hash<ConvexHull2D>` section of `ConvexHull2D_test.cpp` at core
//! `bc9cc12` (tier 3), and the source hash values recomputed in Python from the
//! header's `hash_combine`/`hash_float`/`hash_int` recurrence over explicit
//! little-endian bytes (tier 4).

use openms::kernel::geometry::{ConvexHull2D, Point2D};
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};

fn outline_hull(points: &[(f64, f64)]) -> ConvexHull2D {
    let mut hull = ConvexHull2D::new();
    let points: Vec<Point2D> = points
        .iter()
        .map(|&(rt, mz)| Point2D::new(rt, mz))
        .collect();
    assert!(hull.set_hull_points(&points).is_ok());
    hull
}

fn point_hull(points: &[(f64, f64)]) -> ConvexHull2D {
    let mut hull = ConvexHull2D::new();
    for &(rt, mz) in points {
        assert!(hull.add_point(Point2D::new(rt, mz)).is_ok());
    }
    hull
}

/// `std::unordered_set` needs only `==` and the hash; Rust's `HashSet` also
/// needs `Eq`, asserted here for the NaN-free hulls of the test.
#[derive(Clone, Debug, PartialEq)]
struct Key(ConvexHull2D);

impl Eq for Key {}

impl Hash for Key {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

#[test]
fn class_test_std_hash() {
    let h1 = outline_hull(&[(1.0, 2.0), (3.0, 4.0), (5.0, 0.0)]);
    let h2 = outline_hull(&[(1.0, 2.0), (3.0, 4.0), (5.0, 0.0)]);
    assert!(h1 == h2);
    assert_eq!(h1.source_hash64(), h2.source_hash64());

    let h3 = outline_hull(&[(10.0, 20.0), (30.0, 40.0)]);
    assert!(h1 != h3);
    assert_ne!(h1.source_hash64(), h3.source_hash64());

    let mut set = HashSet::new();
    set.insert(Key(h1.clone()));
    set.insert(Key(h2.clone()));
    set.insert(Key(h3.clone()));
    assert_eq!(set.len(), 2);
    assert!(set.contains(&Key(h1.clone())));
    assert!(set.contains(&Key(h3.clone())));

    let mut map = HashMap::new();
    map.insert(Key(h1.clone()), 1);
    map.insert(Key(h3.clone()), 3);
    assert_eq!(map.len(), 2);
    assert_eq!(map[&Key(h1.clone())], 1);
    assert_eq!(map[&Key(h2.clone())], 1);
    assert_eq!(map[&Key(h3)], 3);

    let h4 = point_hull(&[(1.0, 1.0), (2.0, 2.0), (3.0, 3.0)]);
    let h5 = point_hull(&[(1.0, 1.0), (2.0, 2.0), (3.0, 3.0)]);
    assert!(h4 == h5);
    assert_eq!(h4.source_hash64(), h5.source_hash64());

    let empty1 = ConvexHull2D::new();
    let empty2 = ConvexHull2D::default();
    assert!(empty1 == empty2);
    assert_eq!(empty1.source_hash64(), empty2.source_hash64());
}

#[test]
fn source_hash_values_match_the_header_recurrence() {
    // Python over the recurrence; the point hull's outline is
    // (1,1), (2,2), (3,3), (2,2): the single-point first and last scans are
    // not listed twice, the single-point middle scan is.
    let outline = outline_hull(&[(1.0, 2.0), (3.0, 4.0), (5.0, 0.0)]);
    assert_eq!(outline.source_hash64(), 0x0dce_db1e_12e5_7ace);
    let points = point_hull(&[(1.0, 1.0), (2.0, 2.0), (3.0, 3.0)]);
    assert_eq!(points.source_hash64(), 0xe601_a5d1_f2ea_c6f6);
    assert_eq!(ConvexHull2D::new().source_hash64(), 0x5e64_b8bb_cd0f_ec0a);
}

#[test]
fn map_points_are_one_mz_box_per_scan_in_rt_order() {
    let hull = point_hull(&[(2.0, 5.0), (1.0, 3.0), (2.0, 4.0), (1.0, 7.0), (3.0, 1.0)]);
    let map = hull.map_points();
    assert_eq!(map.len(), 3);
    let expected = [(1.0, 3.0, 7.0), (2.0, 4.0, 5.0), (3.0, 1.0, 1.0)];
    for ((rt, mz), (want_rt, want_min, want_max)) in map.iter().zip(expected) {
        assert_eq!(*rt, want_rt);
        assert_eq!(mz.min_position()[0], want_min);
        assert_eq!(mz.max_position()[0], want_max);
    }
    // A single-point scan is a zero-width box.
    assert_eq!(map[2].1.width(), 0.0);
    // An outline-only hull has no map points, as setHullPoints clears them.
    assert!(
        outline_hull(&[(1.0, 2.0), (3.0, 4.0)])
            .map_points()
            .is_empty()
    );
    assert!(ConvexHull2D::new().map_points().is_empty());
}
