// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `Software_test.cpp` (core SDK bc9cc12): one test per
//! `START_SECTION`, cited by source line. Every literal is tier 3,
//! transcribed from the class test; a value derived rather than transcribed
//! says so.

use std::cmp::Ordering;
use std::collections::HashMap;
use std::hash::{BuildHasher, BuildHasherDefault, DefaultHasher};

use openms::metadata::{CVTerm, Software};

fn hash_of(software: &Software) -> u64 {
    BuildHasherDefault::<DefaultHasher>::default().hash_one(software)
}

// L29 Software(), L34 ~Software()
#[test]
fn section_constructor() {
    let software = Software::default();
    assert_eq!(software, Software::new("", ""));
}

// L38 getName, L43 setName, L49 getVersion, L54 setVersion
#[test]
fn section_accessors() {
    let tmp = Software::default();
    assert_eq!(tmp.name, "");
    assert_eq!(tmp.version, "");
    let mut tmp = Software::new("", "");
    tmp.name = "name".into();
    assert_eq!(tmp.name, "name");
    tmp.version = "0.54".into();
    assert_eq!(tmp.version, "0.54");
}

// L60 Software(const Software&)
#[test]
fn section_copy_constructor() {
    let tmp = Software::new("name", "0.54");
    let tmp2 = tmp.clone();
    assert_eq!(tmp2.version, "0.54");
    assert_eq!(tmp2.name, "name");
}

// L71 Software& operator=(const Software&)
#[test]
fn section_assignment() {
    let tmp = Software::new("name", "0.54");
    let mut tmp2 = Software::default();
    tmp2.clone_from(&tmp);
    assert_eq!(tmp2.version, "0.54");
    assert_eq!(tmp2.name, "name");
    tmp2 = Software::default();
    assert_eq!(tmp2.version, "");
    assert_eq!(tmp2.name, "");
}

// L87 operator==, L102 operator!=
#[test]
fn section_equality() {
    let empty = Software::default();
    let mut edit = Software::default();
    assert!(edit == empty);
    assert!(!(edit != empty));
    edit.version = "0.54".into();
    assert!(edit != empty);
    edit = empty.clone();
    edit.name = "name".into();
    assert!(edit != empty);

    // Derived from `Software::operator==`, which also compares the inherited
    // CVTermList.
    let mut with_term = Software::default();
    with_term
        .cv_terms
        .add(CVTerm::new(
            "MS:1000799",
            "custom unreleased software tool",
            "MS",
        ))
        .unwrap();
    assert!(with_term != empty);
}

/// `Software::operator<` (Software.cpp:37-40): `tie(name_, version_) <
/// tie(rhs.name_, rhs.version_)`. The header declares it; the class test has no
/// section for it, so these expectations are derived from that definition.
#[test]
fn name_then_version_ordering() {
    let a1 = Software::new("A", "1.0");
    let a2 = Software::new("A", "2.0");
    let b0 = Software::new("B", "0.1");
    assert_eq!(a1.cmp_name_version(&a2), Ordering::Less);
    assert_eq!(a2.cmp_name_version(&b0), Ordering::Less);
    assert_eq!(b0.cmp_name_version(&a1), Ordering::Greater);
    assert_eq!(a1.cmp_name_version(&a1.clone()), Ordering::Equal);

    let mut sorted = vec![b0.clone(), a2.clone(), a1.clone()];
    sorted.sort_by(Software::cmp_name_version);
    assert_eq!(sorted, vec![a1.clone(), a2, b0]);

    // CV terms are ignored by the ordering, not by equality.
    let mut a1_term = a1.clone();
    a1_term
        .cv_terms
        .add(CVTerm::new(
            "MS:1000799",
            "custom unreleased software tool",
            "MS",
        ))
        .unwrap();
    assert_eq!(a1.cmp_name_version(&a1_term), Ordering::Equal);
    assert_ne!(a1, a1_term);
}

// L117 [EXTRA] std::hash<Software>
#[test]
fn section_extra_hash() {
    let s1 = Software::new("TestSoftware", "1.0");
    let s2 = Software::new("TestSoftware", "1.0");
    assert!(s1 == s2);
    assert_eq!(hash_of(&s1), hash_of(&s2));

    let s3 = Software::new("OtherSoftware", "1.0");
    assert_ne!(hash_of(&s1), hash_of(&s3));
    let s4 = Software::new("TestSoftware", "2.0");
    assert_ne!(hash_of(&s1), hash_of(&s4));

    // The unordered_set / unordered_map part. `Software` is not `Eq` (its CV
    // terms may carry floating values), so it cannot key a `HashSet`; the same
    // container semantics - bucket by hash, then compare with `==` - are
    // reproduced by hand.
    let mut buckets: HashMap<u64, Vec<(Software, &str)>> = HashMap::new();
    let mut insert = |software: &Software, value: &'static str| {
        let bucket = buckets.entry(hash_of(software)).or_default();
        match bucket.iter_mut().find(|(key, _)| key == software) {
            Some(slot) => slot.1 = value,
            None => bucket.push((software.clone(), value)),
        }
    };
    insert(&s1, "value1");
    insert(&s2, "value1");
    insert(&s3, "value3");
    let size: usize = buckets.values().map(Vec::len).sum();
    assert_eq!(size, 2);
    let lookup = |software: &Software| {
        buckets
            .get(&hash_of(software))
            .and_then(|bucket| bucket.iter().find(|(key, _)| key == software))
            .map(|(_, value)| *value)
    };
    assert_eq!(lookup(&s2), Some("value1"));
    assert_eq!(lookup(&s3), Some("value3"));
}
