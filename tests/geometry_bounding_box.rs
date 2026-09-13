// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! `BoundingBox2D::intersects`, `width`, `height`, `encloses` and `union`.
//!
//! Tier 3: `DBoundingBox_test.cpp` (intersects, encloses, enlarge, equality)
//! and `DIntervalBase_test.cpp` (width, height) literals at OpenMS4-core
//! `bc9cc12`. Tier 1: the `bbox*` lines of
//! `tests/data/isotopes_source_precision/probe.tsv`, printed by the executed C++
//! SDK probe `../oracle/b2-iso-source-precision` (hashes in
//! `tests/data/isotopes_source_precision_provenance.json`). Tier 4: an
//! exhaustive symmetric grid against the interval-overlap definition.

use openms::kernel::geometry::{BoundingBox2D, Point2D};

const PROBE: &str = include_str!("data/isotopes_source_precision/probe.tsv");

fn bbox(min: (f64, f64), max: (f64, f64)) -> BoundingBox2D {
    BoundingBox2D::new(Point2D::new(min.0, min.1), Point2D::new(max.0, max.1)).unwrap()
}

fn probe_fields(kind: &str) -> Vec<Vec<&'static str>> {
    PROBE
        .lines()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|fields| fields[0] == kind)
        .collect()
}

fn hex64(text: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(text, 16).unwrap())
}

/// The source `DIntervalBase<2>` mutators, which move the opposite corner when
/// a new coordinate would invert the interval (`DIntervalBase.h:115-137` and
/// `292-317`), starting from the empty default sentinel.
#[derive(Clone, Copy, Debug)]
struct SourceBox {
    min: [f64; 2],
    max: [f64; 2],
}

impl SourceBox {
    fn empty() -> Self {
        Self {
            min: [f64::MAX; 2],
            max: [-f64::MAX; 2],
        }
    }
    fn set_min(&mut self, position: [f64; 2]) {
        self.set_min_x(position[0]);
        self.set_min_y(position[1]);
    }
    fn set_max(&mut self, position: [f64; 2]) {
        self.set_max_x(position[0]);
        self.set_max_y(position[1]);
    }
    fn set_min_x(&mut self, c: f64) {
        self.min[0] = c;
        if self.min[0] > self.max[0] {
            self.max[0] = self.min[0];
        }
    }
    fn set_min_y(&mut self, c: f64) {
        self.min[1] = c;
        if self.min[1] > self.max[1] {
            self.max[1] = self.min[1];
        }
    }
    fn set_max_x(&mut self, c: f64) {
        self.max[0] = c;
        if self.min[0] > self.max[0] {
            self.min[0] = self.max[0];
        }
    }
    fn set_max_y(&mut self, c: f64) {
        self.max[1] = c;
        if self.min[1] > self.max[1] {
            self.min[1] = self.max[1];
        }
    }
    fn plus(position: [f64; 2], offset: f64) -> [f64; 2] {
        [position[0] + offset, position[1] + offset]
    }
    fn native(self) -> BoundingBox2D {
        bbox((self.min[0], self.min[1]), (self.max[0], self.max[1]))
    }
}

#[test]
fn upstream_intersects_section_matches_literals_and_executed_cpp() {
    // DBoundingBox_test.cpp:205-312. p1 = (-1, -2), p2 = (3, 4), one, two.
    let mut r2 = SourceBox::empty();
    r2.set_min([-1.0, -2.0]);
    r2.set_max([3.0, 4.0]);
    let mut r3 = r2;
    let mut steps: Vec<(&str, SourceBox, bool)> = vec![("copy", r3, true)];
    r3.set_max_x(10.0);
    steps.push(("maxX10", r3, true));
    r3.set_max(SourceBox::plus(r2.max, 1.0));
    steps.push(("max_plus_one", r3, true));
    r3.set_min(SourceBox::plus(r2.max, 1.0));
    r3.set_max(SourceBox::plus(r2.max, 2.0));
    steps.push(("beyond_max", r3, false));
    r3.set_min(r2.min);
    r3.set_min_x(10.0);
    r3.set_max(SourceBox::plus(r3.min, 1.0));
    steps.push(("right_of", r3, false));
    r3.set_min_x(-10.0);
    r3.set_min_y(-10.0);
    r3.set_max(SourceBox::plus(r3.min, 1.0));
    steps.push(("lower_left", r3, false));
    // Each remaining step sets minX, minY, maxX, maxY in that order, except
    // `upper_left`, which sets minX and minY and then max = min + one.
    let explicit: [(&str, Option<[f64; 4]>, bool); 15] = [
        ("below_left", Some([-10.0, -10.0, 0.0, -9.0]), false),
        ("below_wide", Some([-10.0, -10.0, 10.0, -9.0]), false),
        ("left_mid", Some([-10.0, 0.0, -9.0, 1.0]), false),
        ("upper_left", None, false),
        ("left_tall", Some([-10.0, 0.0, -9.0, 10.0]), false),
        ("right_tall", Some([9.0, 0.0, 10.0, 10.0]), false),
        ("right_tall_again", Some([9.0, 0.0, 10.0, 10.0]), false),
        ("right_low", Some([9.0, -5.0, 10.0, 0.0]), false),
        ("right_span", Some([9.0, -5.0, 10.0, 5.0]), false),
        ("overlap_corner", Some([-5.0, -5.0, 0.0, 0.0]), true),
        ("overlap_wide", Some([-5.0, -5.0, 5.0, 0.0]), true),
        ("cover", Some([-5.0, -5.0, 5.0, 5.0]), true),
        ("zero_width", Some([0.0, -5.0, 0.0, 0.0]), true),
        ("half_right", Some([0.0, -5.0, 5.0, 0.0]), true),
        ("half_right_tall", Some([0.0, -5.0, 5.0, 5.0]), true),
    ];
    for (label, corners, expected) in explicit {
        match corners {
            Some(corners) => {
                r3.set_min_x(corners[0]);
                r3.set_min_y(corners[1]);
                r3.set_max_x(corners[2]);
                r3.set_max_y(corners[3]);
            }
            None => {
                r3.set_min_x(-10.0);
                r3.set_min_y(10.0);
                r3.set_max(SourceBox::plus(r3.min, 1.0));
            }
        }
        steps.push((label, r3, expected));
    }

    let executed = probe_fields("bbox");
    assert_eq!(executed.len(), steps.len());
    assert_eq!(
        steps.len(),
        21,
        "all 21 TEST_EQUAL assertions of the section"
    );
    let r2 = r2.native();
    for ((label, source, expected), cpp) in steps.into_iter().zip(executed) {
        assert_eq!(cpp[1], label);
        let corners = [cpp[2], cpp[3], cpp[4], cpp[5]].map(hex64);
        assert_eq!(
            [source.min[0], source.min[1], source.max[0], source.max[1]].map(f64::to_bits),
            corners.map(f64::to_bits),
            "{label}: transcribed setters reproduce the C++ corners"
        );
        let r3 = source.native();
        assert_eq!(r2.intersects(r3), expected, "{label}");
        assert_eq!(r3.intersects(r2), expected, "{label} (reversed)");
        assert_eq!(cpp[6] == "1", expected, "{label}: C++ r2.intersects(r3)");
        assert_eq!(cpp[7] == "1", expected, "{label}: C++ r3.intersects(r2)");
    }
}

#[test]
fn touching_borders_and_points_intersect_as_in_cpp() {
    let executed = probe_fields("bbox_touch");
    assert_eq!(executed.len(), 1);
    let a = bbox((0.0, 0.0), (1.0, 1.0));
    let touch = bbox((1.0, 1.0), (2.0, 2.0));
    let point = bbox((0.5, 0.5), (0.5, 0.5));
    let native = [
        a.intersects(touch),
        touch.intersects(a),
        a.intersects(point),
    ];
    let cpp = [executed[0][1], executed[0][2], executed[0][3]].map(|value| value == "1");
    assert_eq!(native, cpp);
    assert_eq!(native, [true; 3]);
    // A gap of one representable step in either dimension separates boxes.
    let after_one = f64::from_bits(1.0_f64.to_bits() + 1);
    let beyond_rt = bbox((after_one, 0.0), (2.0, 1.0));
    let beyond_mz = bbox((0.0, after_one), (1.0, 2.0));
    assert!(!a.intersects(beyond_rt) && !beyond_rt.intersects(a));
    assert!(!a.intersects(beyond_mz) && !beyond_mz.intersects(a));
}

#[test]
fn upstream_width_and_height_sections() {
    // DIntervalBase_test.cpp:256-264 with p1 = (5, 17.5), p2 = (65, -57.5). The
    // source constructor normalises the corners; BoundingBox2D::new requires
    // ordered corners instead.
    assert!(BoundingBox2D::new(Point2D::new(5.0, 17.5), Point2D::new(65.0, -57.5)).is_err());
    let tmp = bbox((5.0, -57.5), (65.0, 17.5));
    assert!((tmp.width() - 60.0).abs() <= 1e-5);
    assert!((tmp.height() - 75.0).abs() <= 1e-5);
    let executed = probe_fields("bbox_extent");
    assert_eq!(executed.len(), 1);
    let cpp = executed[0][1..7]
        .iter()
        .map(|text| hex64(text))
        .collect::<Vec<_>>();
    let native = [
        tmp.min().rt,
        tmp.min().mz,
        tmp.max().rt,
        tmp.max().mz,
        tmp.width(),
        tmp.height(),
    ];
    assert_eq!(
        native.map(f64::to_bits).to_vec(),
        cpp.iter().map(|value| value.to_bits()).collect::<Vec<_>>()
    );
    // Degenerate and extreme boxes: zero extent, and overflow without panic.
    let point = bbox((3.0, 4.0), (3.0, 4.0));
    assert_eq!((point.width(), point.height()), (0.0, 0.0));
    let huge = bbox((-f64::MAX, -f64::MAX), (f64::MAX, f64::MAX));
    assert_eq!(huge.width(), f64::INFINITY);
    assert_eq!(huge.height(), f64::INFINITY);
    assert_eq!(huge.rt_range().max - huge.rt_range().min, huge.width());
}

#[test]
fn upstream_encloses_sections() {
    // DBoundingBox_test.cpp:140-203: both overloads share these literals.
    let tmp = bbox((100.0, 200.0), (300.0, 400.0));
    for (x, y, expected) in [
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
    ] {
        assert_eq!(
            tmp.encloses(Point2D::new(x, y)).unwrap(),
            expected,
            "({x}, {y})"
        );
    }
}

#[test]
fn upstream_enlarge_and_equality_sections_through_union() {
    // DBoundingBox_test.cpp:102-138. The source enlarges an empty default box;
    // here emptiness is `None` and enlargement is a union with a point box.
    let enlarge = |bounds: Option<BoundingBox2D>, x: f64, y: f64| {
        let point = bbox((x, y), (x, y));
        Some(bounds.map_or(point, |bounds| bounds.union(point)))
    };
    let encloses = |bounds: Option<BoundingBox2D>, x: f64, y: f64| {
        bounds.is_some_and(|bounds| bounds.encloses(Point2D::new(x, y)).unwrap())
    };
    let mut bb2h = None;
    assert!(!encloses(bb2h, 11.0, 13.0));
    assert!(!encloses(bb2h, 10.0, 1.0));
    bb2h = enlarge(bb2h, 11.0, 13.0);
    assert!(encloses(bb2h, 11.0, 13.0));
    assert!(!encloses(bb2h, 10.0, 1.0));
    bb2h = enlarge(bb2h, 9.0, 0.0);
    assert!(encloses(bb2h, 11.0, 13.0));
    assert!(encloses(bb2h, 10.0, 1.0));
    let bb2 = enlarge(None, 9.0, 0.0).unwrap();
    let copy = bb2;
    assert_eq!(bb2, copy);
}

#[test]
fn intersects_equals_interval_overlap_on_an_exhaustive_grid() {
    let values = [0.0, 1.0, 2.0, 3.0];
    let mut boxes = Vec::new();
    for &rt_min in &values {
        for &rt_max in values.iter().filter(|&&v| v >= rt_min) {
            for &mz_min in &values {
                for &mz_max in values.iter().filter(|&&v| v >= mz_min) {
                    boxes.push(bbox((rt_min, mz_min), (rt_max, mz_max)));
                }
            }
        }
    }
    assert_eq!(boxes.len(), 100);
    for &a in &boxes {
        for &b in &boxes {
            let overlap = a.rt_range().min.max(b.rt_range().min)
                <= a.rt_range().max.min(b.rt_range().max)
                && a.mz_range().min.max(b.mz_range().min) <= a.mz_range().max.min(b.mz_range().max);
            assert_eq!(a.intersects(b), overlap, "{a:?} {b:?}");
            assert_eq!(a.intersects(b), b.intersects(a));
        }
    }
}
