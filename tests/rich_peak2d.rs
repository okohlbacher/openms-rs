// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

// The class-test transcription keeps the source's construct-then-set shape.
#![allow(clippy::field_reassign_with_default)]

//! Section-by-section port of `RichPeak2D_test.cpp`
//! (source revision bc9cc12514c768385ce121d6ca4bb710fe1983c4, tier 3:
//! transcribed class-test literals, no C++ execution).
//!
//! Source `setMetaValue("cluster_id", 4711)` stores a `DataValue` integer; the
//! Rust metadata map stores `MetaValue::from(4711i64)`. The two `[EXTRA]`
//! sections address metadata by registry index 2, which the source
//! `MetaInfoRegistry` pre-registers as `"cluster_id"`; Rust keys by name only.

use openms::concept::HasUniqueId;
use openms::kernel::{Peak2D, RichPeak2D};
use openms::metadata::MetaValue;

/// `TEST_REAL_SIMILAR` default: absolute difference within 1e-5, or ratio
/// within 1 + 1e-5.
fn real_similar(a: f64, b: f64) {
    let absdiff = (a - b).abs();
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{a} !~ {b} (absdiff {absdiff}, ratio {ratio})"
    );
}

fn cluster_id() -> MetaValue {
    MetaValue::from(4711i64)
}

/// Section `RichPeak2D()`.
#[test]
fn default_constructor() {
    let p = RichPeak2D::default();
    assert_eq!(p.peak, Peak2D::default());
    assert!(p.metadata.is_empty());
    assert_eq!(p.unique_id, 0);
    let boxed: Box<RichPeak2D> = Box::default();
    assert_eq!(*boxed, p);
}

/// Section `~RichPeak2D()`: dropping releases owned metadata.
#[test]
fn destructor() {
    let mut p = RichPeak2D::default();
    p.metadata.insert("cluster_id".into(), cluster_id());
    drop(p);
}

/// Section `RichPeak2D(const RichPeak2D& p)`.
#[test]
fn copy_constructor() {
    let mut p = RichPeak2D::default();
    p.intensity = 123.456f32;
    p.metadata.insert("cluster_id".into(), cluster_id());
    let copy_of_p = p.clone();
    real_similar(f64::from(copy_of_p.intensity), f64::from(123.456f32));
    assert_eq!(copy_of_p.metadata["cluster_id"], cluster_id());
}

/// Section `RichPeak2D(RichPeak2D&& rhs)`. The source asserts a `noexcept`
/// move constructor; a Rust move is always an infallible bitwise transfer,
/// so only the transferred values are checked.
#[test]
fn move_constructor() {
    let mut p = RichPeak2D::default();
    p.intensity = 123.456f32;
    p.metadata.insert("cluster_id".into(), cluster_id());
    p.position = [21.21, 22.22];
    let copy_of_p = p;
    real_similar(f64::from(copy_of_p.intensity), f64::from(123.456f32));
    assert_eq!(copy_of_p.metadata["cluster_id"], cluster_id());
    let i2 = copy_of_p.intensity;
    let pos2 = copy_of_p.position;
    real_similar(f64::from(i2), 123.456);
    real_similar(pos2[0], 21.21);
    real_similar(pos2[1], 22.22);
}

/// Section `RichPeak2D(const Peak2D& p)`: conversion starts with empty
/// metadata and a cleared unique ID.
#[test]
fn constructor_from_peak2d() {
    let mut p = Peak2D::default();
    p.intensity = 123.456f32;
    let copy_of_p = RichPeak2D::from(p);
    real_similar(f64::from(copy_of_p.intensity), f64::from(123.456f32));
    assert!(copy_of_p.metadata.is_empty());
    assert_eq!(copy_of_p.unique_id, 0);
    assert_eq!(RichPeak2D::from(&p), copy_of_p);
}

/// Section `explicit RichPeak2D(const PositionType& pos, const IntensityType in)`.
#[test]
fn member_constructor() {
    let p = RichPeak2D::from_position([21.21, 22.22], 123.456f32);
    let copy_of_p = p.clone();
    real_similar(f64::from(copy_of_p.intensity), 123.456);
    real_similar(copy_of_p.position[0], 21.21);
    real_similar(copy_of_p.position[1], 22.22);
    assert_eq!(RichPeak2D::new(21.21, 22.22, 123.456f32), p);
}

/// Section `RichPeak2D& operator=(const RichPeak2D& rhs)`.
#[test]
fn assignment_operator() {
    let mut p = RichPeak2D::default();
    p.intensity = 123.456f32;
    p.metadata.insert("cluster_id".into(), cluster_id());
    let mut copy_of_p = RichPeak2D::default();
    assert!(copy_of_p.metadata.is_empty());
    copy_of_p = p.clone();
    real_similar(f64::from(copy_of_p.intensity), f64::from(123.456f32));
    assert_eq!(copy_of_p.metadata["cluster_id"], cluster_id());
}

/// Section `RichPeak2D& operator=(const Peak2D& rhs)`: plain-point assignment
/// clears metadata and the unique ID. Rust returns the previous value instead
/// of destroying it in place.
#[test]
fn assignment_from_peak2d_clears_meta_info() {
    let mut p = Peak2D::default();
    p.intensity = 123.456f32;
    let mut copy_of_p = RichPeak2D::default();
    copy_of_p.metadata.insert("cluster_id".into(), cluster_id());
    copy_of_p.unique_id = 5;
    let previous = copy_of_p.replace_from_peak(p);
    real_similar(f64::from(copy_of_p.intensity), f64::from(123.456f32));
    assert!(copy_of_p.metadata.is_empty());
    assert_eq!(copy_of_p.unique_id, 0);
    assert_eq!(previous.metadata["cluster_id"], cluster_id());
    assert_eq!(previous.unique_id, 5);
}

/// Section `bool operator==(const RichPeak2D& rhs) const`.
#[test]
fn equality_operator() {
    let mut p1 = RichPeak2D::default();
    let p2 = RichPeak2D::default();
    assert!(p1 == p2);
    p1.intensity = 5.0f32;
    assert!(!(p1 == p2));
    let mut p2 = p2;
    p2.intensity = 5.0f32;
    assert!(p1 == p2);
    p1.metadata.insert("cluster_id".into(), cluster_id());
    assert!(!(p1 == p2));
    p1.metadata.remove("cluster_id");
    assert!(p1 == p2);
    // Inherited UniqueIdInterface equality participates too.
    p1.unique_id = 1;
    assert!(!(p1 == p2));
}

/// Section `bool operator!=(const RichPeak2D& rhs) const`.
#[test]
fn inequality_operator() {
    let mut p1 = RichPeak2D::default();
    let mut p2 = RichPeak2D::default();
    assert!(!(p1 != p2));
    p1.intensity = 5.0f32;
    assert!(p1 != p2);
    p2.intensity = 5.0f32;
    assert!(!(p1 != p2));
    p1.metadata.insert("cluster_id".into(), cluster_id());
    assert!(p1 != p2);
    p1.metadata.remove("cluster_id");
    assert!(!(p1 != p2));
}

/// Section `[EXTRA] meta info with copy constructor`: the copy owns its
/// metadata, so later edits to the original do not leak into it. Source
/// registry index 2 is `"cluster_id"`.
#[test]
fn meta_info_with_copy_constructor() {
    let mut p = RichPeak2D::default();
    p.metadata.insert("cluster_id".into(), "bla".into());
    let p2 = p.clone();
    assert_eq!(p.metadata["cluster_id"].as_str().unwrap(), "bla");
    assert_eq!(p2.metadata["cluster_id"].as_str().unwrap(), "bla");
    p.metadata.insert("cluster_id".into(), "bluff".into());
    assert_eq!(p.metadata["cluster_id"].as_str().unwrap(), "bluff");
    assert_eq!(p2.metadata["cluster_id"].as_str().unwrap(), "bla");
}

/// Section `[EXTRA] meta info with assignment`.
#[test]
fn meta_info_with_assignment() {
    let mut p = RichPeak2D::default();
    p.metadata.insert("cluster_id".into(), "bla".into());
    let mut p2 = RichPeak2D::default();
    assert!(p2.metadata.is_empty());
    p2 = p.clone();
    assert_eq!(p.metadata["cluster_id"].as_str().unwrap(), "bla");
    assert_eq!(p2.metadata["cluster_id"].as_str().unwrap(), "bla");
    p.metadata.insert("cluster_id".into(), "bluff".into());
    assert_eq!(p.metadata["cluster_id"].as_str().unwrap(), "bluff");
    assert_eq!(p2.metadata["cluster_id"].as_str().unwrap(), "bla");
}

/// Native check of the inherited `UniqueIdInterface` surface that the source
/// test does not exercise for this class.
#[test]
fn inherited_unique_id_interface() {
    let mut p = RichPeak2D::default();
    assert!(p.has_invalid_unique_id());
    p.set_unique_id(42);
    assert!(p.has_valid_unique_id());
    assert_eq!(p.unique_id(), 42);
    assert_eq!(p.clear_unique_id(), 1);
    p.set_unique_id_from_str("feature_17").unwrap();
    assert_eq!(p.unique_id, 17);
}
