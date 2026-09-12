// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

// The class-test transcription keeps the source's construct-then-set shape.
#![allow(clippy::field_reassign_with_default)]

//! Section-by-section port of `FeatureHandle_test.cpp`
//! (source revision bc9cc12514c768385ce121d6ca4bb710fe1983c4, tier 3:
//! transcribed class-test literals, no C++ execution).
//!
//! The source `ElementType` is `FeatureMap::value_type` = `Feature`; its
//! measured values live in the composed `BaseFeature`, which is what
//! `FeatureHandle::new` reads.

use openms::concept::HasUniqueId;
use openms::kernel::{BaseFeature, ConsensusFeature, Feature, FeatureHandle, Peak2D};
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

fn hash(value: &impl Hash) -> u64 {
    let mut hasher = DefaultHasher::new();
    value.hash(&mut hasher);
    hasher.finish()
}

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

fn element_with_id(unique_id: u64) -> Feature {
    let mut e = Feature::default();
    e.unique_id = unique_id;
    e
}

/// Section `FeatureHandle()`.
#[test]
fn default_constructor() {
    let handle = FeatureHandle::default();
    assert_eq!(handle.map_index, 0);
    assert_eq!(handle.unique_id, 0);
    assert_eq!(handle.charge, 0);
    assert_eq!(handle.width, 0.0);
    assert_eq!((handle.rt, handle.mz, handle.intensity), (0.0, 0.0, 0.0));
}

/// Section `virtual ~FeatureHandle()`.
#[test]
fn destructor() {
    let boxed: Box<FeatureHandle> = Box::default();
    drop(boxed);
}

/// Section `FeatureHandle& operator=(const FeatureHandle& rhs)`.
#[test]
fn assignment_operator() {
    let e = element_with_id(2);
    let it = FeatureHandle::new(1, &e.base);
    let mut it_copy = FeatureHandle::default();
    assert_eq!(it_copy.map_index, 0);
    it_copy = it;
    assert!(it.unique_id == it_copy.unique_id);
    assert!(it.map_index == it_copy.map_index);
    assert!(it.intensity == it_copy.intensity);
    assert!((it.rt, it.mz) == (it_copy.rt, it_copy.mz));
}

/// Section `FeatureHandle(const FeatureHandle& rhs)`.
#[test]
fn copy_constructor() {
    let e = element_with_id(2);
    let it = FeatureHandle::new(1, &e.base);
    let it_copy = it;
    assert!(it.unique_id == it_copy.unique_id);
    assert!(it.map_index == it_copy.map_index);
    assert!(it.intensity == it_copy.intensity);
    assert!((it.rt, it.mz) == (it_copy.rt, it_copy.mz));
}

/// Sections `void setCharge(ChargeType charge)` and `ChargeType getCharge()
/// const` (the getter section is `NOT_TESTABLE` in the source).
#[test]
fn set_and_get_charge() {
    let mut fh = FeatureHandle::default();
    fh.charge = -17;
    assert_eq!(fh.charge, -17);
    fh.charge = -1717;
    assert_eq!(fh.charge, -1717);
}

/// Sections `void setWidth(WidthType width)` and `WidthType getWidth() const`
/// (`NOT_TESTABLE`). The public field accepts the source's negative literal;
/// `validate` rejects it, which is a documented native strengthening.
#[test]
fn set_and_get_width() {
    let mut fh_tmp = FeatureHandle::default();
    fh_tmp.width = 10.7;
    real_similar(f64::from(fh_tmp.width), 10.7);
    assert!(fh_tmp.validate().is_ok());
    fh_tmp.width = -8.9;
    real_similar(f64::from(fh_tmp.width), -8.9);
    assert!(fh_tmp.validate().is_err());
}

/// Section `FeatureHandle(UInt64 map_index, const Peak2D& point, UInt64
/// element_index)`.
#[test]
fn constructor_from_map_index_point_and_element_index() {
    let e = Feature::default();
    let it = FeatureHandle::from_peak(1, Peak2D::new(e.rt, e.mz, e.intensity), 2);
    assert!(it.unique_id == 2);
    assert!(it.map_index == 1);
    assert!((it.rt, it.mz) == (e.rt, e.mz));
    assert_eq!(it.charge, 0);
    assert_eq!(it.width, 0.0);
    let point = Peak2D::new(44324.6, 867.4, 12.5);
    let it = FeatureHandle::from_peak(7, point, 9);
    assert_eq!((it.rt, it.mz, it.intensity), (44324.6, 867.4, 12.5));
    assert_eq!(it.key(), (7, 9));
}

/// Section `FeatureHandle(UInt64 map_index, const BaseFeature& feature)`.
#[test]
fn constructor_from_map_index_and_base_feature() {
    let mut f = Feature::default();
    f.charge = -17;
    f.rt = 44324.6;
    f.mz = 867.4;
    f.unique_id = 23;
    let f_cref: &BaseFeature = &f.base;
    let fh = FeatureHandle::new(99, f_cref);
    assert_eq!(fh.map_index, 99);
    assert_eq!(fh.unique_id, 23);
    assert_eq!(fh.rt, 44324.6);
    assert_eq!(fh.mz, 867.4);
    assert_eq!(fh.charge, -17);
}

/// Section `FeatureHandleMutable_& asMutable() const`. The source mutates a
/// handle through a const reference; Rust mutates through `&mut`, and a handle
/// stored in a consensus feature is edited by replacing the handle set, which
/// keeps the `(map_index, unique_id)` ordering invariant the source relies on.
#[test]
fn as_mutable_equivalent() {
    let mut f = ConsensusFeature::new();
    f.charge = -17;
    f.rt = 44324.6;
    f.mz = 867.4;
    f.unique_id = 23;
    let f_cref: &BaseFeature = &f.base;
    let mut fh = FeatureHandle::new(99, f_cref);
    fh.rt = -64544.3;
    assert_eq!(fh.map_index, 99);
    assert_eq!(fh.unique_id, 23);
    assert_eq!(fh.rt, -64544.3);
    assert_eq!(fh.mz, 867.4);
    assert_eq!(fh.charge, -17);

    let mut consensus = ConsensusFeature::new();
    consensus.insert(fh).unwrap();
    let mut stored = consensus.handles()[0];
    stored.rt = 1.0;
    consensus.set_handles(vec![stored]).unwrap();
    assert_eq!(consensus.handles()[0].rt, 1.0);
    assert_eq!(consensus.handles()[0].key(), (99, 23));
}

/// Section `bool operator!=(const FeatureHandle& i) const`.
#[test]
fn inequality_operator() {
    let e = element_with_id(2);
    let it1 = FeatureHandle::new(1, &e.base);
    let it2 = FeatureHandle::new(2, &e.base);
    assert!(!(it1 == it2));
    assert!(it1 != it2);
}

/// Section `bool operator==(const FeatureHandle& i) const`.
#[test]
fn equality_operator() {
    let e = element_with_id(2);
    let it1 = FeatureHandle::new(2, &e.base);
    let it2 = FeatureHandle::new(2, &e.base);
    assert!(it1 == it2);
    // Every stored value participates, as in the source operator.
    let mut changed = it2;
    changed.width = 1.0;
    assert!(it1 != changed);
    changed = it2;
    changed.charge = 1;
    assert!(it1 != changed);
}

/// Section `UInt64 getMapIndex() const`.
#[test]
fn get_map_index() {
    let e = element_with_id(2);
    let it = FeatureHandle::new(1, &e.base);
    assert!(it.map_index == 1);
}

/// Section `void setMapIndex(UInt64 i)`.
#[test]
fn set_map_index() {
    let mut it = FeatureHandle::default();
    it.map_index = 2;
    it.set_unique_id(77);
    assert!(it.map_index == 2);
    assert_eq!(it.unique_id, 77);
}

/// Section `[FeatureHandle::IndexLess] bool operator()(...)`: map index first,
/// then unique ID. The source test assigns `lhs.setUniqueId` twice (29 wins);
/// the literal sequence is reproduced, and the tuple key ordering is the
/// native comparator.
#[test]
fn index_less_is_key_ordering() {
    let mut lhs = FeatureHandle::default();
    let mut rhs = FeatureHandle::default();
    lhs.map_index = 2;
    lhs.unique_id = 77;
    rhs.map_index = 4;
    lhs.unique_id = 29;
    assert!(lhs.key() < rhs.key());
    assert!(rhs.key() >= lhs.key());
    // Equal map indices fall back to the unique ID.
    rhs.map_index = 2;
    rhs.unique_id = 30;
    assert!(lhs.key() < rhs.key());
    rhs.unique_id = 29;
    assert!(lhs.key() >= rhs.key());
}

/// Source `operator<<`: banner plus five labelled lines; charge and width are
/// not printed.
#[test]
fn stream_output_layout() {
    let mut fh = FeatureHandle::default();
    fh.rt = 44324.6;
    fh.mz = 867.4;
    fh.intensity = 12.5;
    fh.map_index = 99;
    fh.unique_id = 23;
    fh.charge = -17;
    assert_eq!(
        fh.to_string(),
        "---------- FeatureHandle -----------------\n\
         RT: 44324.6\n\
         m/z: 867.4\n\
         Intensity: 12.5\n\
         Map Index: 99\n\
         Element Id: 23\n"
    );
}

/// Source `std::hash<FeatureHandle>`: equal handles hash equally, both zero
/// signs agree, and each of the seven hashed members changes the digest.
#[test]
fn hash_covers_all_members_and_normalises_signed_zero() {
    let base = FeatureHandle::default();
    let mut negative_zero = base;
    negative_zero.rt = -0.0;
    negative_zero.mz = -0.0;
    negative_zero.intensity = -0.0;
    negative_zero.width = -0.0;
    assert_eq!(base, negative_zero);
    assert_eq!(hash(&base), hash(&negative_zero));
    let edits: [fn(&mut FeatureHandle); 7] = [
        |h| h.rt = 1.0,
        |h| h.mz = 1.0,
        |h| h.intensity = 1.0,
        |h| h.unique_id = 1,
        |h| h.map_index = 1,
        |h| h.charge = 1,
        |h| h.width = 1.0,
    ];
    for edit in edits {
        let mut changed = base;
        edit(&mut changed);
        assert_ne!(hash(&base), hash(&changed));
    }
}

/// Inherited `UniqueIdInterface` operations over the `unique_id` field.
#[test]
fn inherited_unique_id_interface() {
    let mut fh = FeatureHandle::default();
    assert!(fh.has_invalid_unique_id());
    fh.set_unique_id(23);
    assert!(fh.has_valid_unique_id());
    assert_eq!(fh.unique_id(), 23);
    let mut other = FeatureHandle::default();
    fh.swap_unique_id(&mut other);
    assert_eq!((fh.unique_id, other.unique_id), (0, 23));
    assert_eq!(other.clear_unique_id(), 1);
    fh.set_unique_id_from_str("f_23").unwrap();
    assert_eq!(fh.unique_id, 23);
}
