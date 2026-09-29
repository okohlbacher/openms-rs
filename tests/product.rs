// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `src/tests/class_tests/openms/source/Product_test.cpp` (core SDK
//! `bc9cc12`), section by section with unchanged literals. The getters and
//! setters are the public fields of `openms::metadata::Product`; the
//! `setMetaValue` of the `CVTermList` base is its `cv_terms.metadata` map.

use std::collections::hash_map::DefaultHasher;
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};

use openms::metadata::Product;

fn with_label() -> Product {
    let mut tmp = Product {
        mz: 47.11,
        isolation_window_upper_offset: 22.7,
        isolation_window_lower_offset: 22.8,
        ..Product::default()
    };
    tmp.cv_terms.metadata.insert("label".into(), "label".into());
    tmp
}

fn hash(product: &Product) -> u64 {
    let mut hasher = DefaultHasher::new();
    product.hash(&mut hasher);
    hasher.finish()
}

// START_SECTION((Product())) / ((~Product())): construction and Drop.

// START_SECTION((double getMZ() const)) and ((void setMZ(double mz)))
#[test]
fn mz() {
    let mut tmp = Product::default();
    assert_eq!(tmp.mz, 0.0);
    tmp.mz = 47.11;
    assert_eq!(tmp.mz, 47.11);
}

// START_SECTION((double getIsolationWindowUpperOffset() const)) and the setter
#[test]
fn isolation_window_upper_offset() {
    let mut tmp = Product::default();
    assert_eq!(tmp.isolation_window_upper_offset, 0.0);
    tmp.isolation_window_upper_offset = 22.7;
    assert_eq!(tmp.isolation_window_upper_offset, 22.7);
}

// START_SECTION((double getIsolationWindowLowerOffset() const)) and the setter
#[test]
fn isolation_window_lower_offset() {
    let mut tmp = Product::default();
    assert_eq!(tmp.isolation_window_lower_offset, 0.0);
    tmp.isolation_window_lower_offset = 22.8;
    assert_eq!(tmp.isolation_window_lower_offset, 22.8);
}

// START_SECTION((Product(const Product& source)))
#[test]
fn copy_constructor() {
    let tmp = with_label();
    let tmp2 = tmp.clone();
    assert_eq!(tmp2.cv_terms.metadata["label"].as_str().unwrap(), "label");
    assert_eq!(tmp2.isolation_window_upper_offset, 22.7);
    assert_eq!(tmp2.isolation_window_lower_offset, 22.8);
    assert_eq!(tmp2.mz, 47.11);
}

// START_SECTION((Product& operator= (const Product& source)))
#[test]
fn assignment() {
    let tmp = with_label();
    let mut tmp2 = Product::default();
    assert_eq!(tmp2.mz, 0.0);
    tmp2 = tmp.clone();
    assert_eq!(tmp2.cv_terms.metadata["label"].as_str().unwrap(), "label");
    assert_eq!(tmp2.isolation_window_upper_offset, 22.7);
    assert_eq!(tmp2.isolation_window_lower_offset, 22.8);
    assert_eq!(tmp2.mz, 47.11);
    tmp2 = Product::default();
    assert!(!tmp2.cv_terms.metadata.contains_key("label"));
    assert_eq!(tmp2.isolation_window_upper_offset, 0.0);
    assert_eq!(tmp2.isolation_window_lower_offset, 0.0);
    assert_eq!(tmp2.mz, 0.0);
}

// START_SECTION((bool operator== (const Product& rhs) const)) and
// ((bool operator!= (const Product& rhs) const)), whose bodies agree.
#[test]
fn equality() {
    let mut tmp = Product::default();
    let mut tmp2 = Product::default();
    assert!(tmp == tmp2);
    let unequal = tmp != tmp2;
    assert!(!unequal);
    tmp2.mz = 47.11;
    assert!(tmp != tmp2);
    tmp2 = tmp.clone();
    tmp.isolation_window_upper_offset = 22.7;
    assert!(tmp != tmp2);
    tmp2 = tmp.clone();
    tmp.isolation_window_lower_offset = 22.8;
    assert!(tmp != tmp2);
    tmp2 = tmp.clone();
    tmp.cv_terms.metadata.insert("label".into(), "label".into());
    assert!(tmp != tmp2);
}

// START_SECTION(([EXTRA] std::hash<Product>))
#[test]
fn hash_in_sets_and_maps() {
    let p1 = Product {
        mz: 100.5,
        isolation_window_lower_offset: 10.0,
        isolation_window_upper_offset: 20.0,
        ..Product::default()
    };
    let p2 = p1.clone();
    assert_eq!(hash(&p1), hash(&p2));
    let p3 = Product {
        mz: 200.5,
        ..p1.clone()
    };
    assert_ne!(hash(&p1), hash(&p3));

    // Product holds f64 fields and is not Eq; the source's unordered
    // containers are exercised through the hash-keyed wrappers below.
    let mut set = HashSet::new();
    for product in [&p1, &p2, &p3] {
        set.insert((hash(product), format!("{product:?}")));
    }
    assert_eq!(set.len(), 2);
    let mut map = HashMap::new();
    map.insert((hash(&p1), format!("{p1:?}")), "first");
    map.insert((hash(&p3), format!("{p3:?}")), "third");
    assert_eq!(map[&(hash(&p2), format!("{p2:?}"))], "first");
    assert_eq!(map.len(), 2);
}
