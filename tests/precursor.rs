// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of the `Precursor_test.cpp` sections (core SDK bc9cc12) for the
//! activation-method string lists and `std::hash<Precursor>`, cited by source
//! line. The record's accessors, copy and equality sections are covered by
//! `tests/metadata.rs` and `tests/precursor_workflow.rs`. Every literal is
//! tier 3, transcribed from the class test.

use std::collections::HashMap;
use std::hash::{BuildHasher, BuildHasherDefault, DefaultHasher};

use openms::kernel::Precursor;
use openms::metadata::{ActivationMethod, DriftTimeUnit};

fn hash_of(precursor: &Precursor) -> u64 {
    BuildHasherDefault::<DefaultHasher>::default().hash_one(precursor)
}

// L90 StringList getActivationMethodsAsString() const
#[test]
fn section_activation_methods_as_string() {
    let mut tmp = Precursor::default();
    tmp.activation_methods.insert(ActivationMethod::Cid);
    let result = tmp.activation_methods_as_names();
    assert_eq!(result.len(), 1);
    assert_eq!(result[0], "Collision-induced dissociation");
}

// L100 StringList getActivationMethodsAsShortString() const
#[test]
fn section_activation_methods_as_short_string() {
    let mut tmp = Precursor::default();
    tmp.activation_methods.insert(ActivationMethod::Cid);
    let result = tmp.activation_methods_as_short_names();
    assert_eq!(result.len(), 1);
    assert_eq!(result[0], "CID");
}

// L110 static StringList getAllNamesOfActivationMethods()
#[test]
fn section_all_names_of_activation_methods() {
    let result = Precursor::all_names_of_activation_methods();
    // SIZE_OF_ACTIVATIONMETHOD
    assert_eq!(result.len(), ActivationMethod::ALL.len());
    assert_eq!(result.len(), 19);
    // Derived: declaration order, so the first is CID's full name.
    assert_eq!(result[0], "Collision-induced dissociation");
}

// L115 static StringList getAllShortNamesOfActivationMethods()
#[test]
fn section_all_short_names_of_activation_methods() {
    let result = Precursor::all_short_names_of_activation_methods();
    assert_eq!(result.len(), 19);
    assert_eq!(result[0], "CID");
}

/// Derived from `getActivationMethodsAsString`: the set iterates in enum order,
/// whatever the insertion order.
#[test]
fn activation_method_lists_follow_the_set_order() {
    let mut tmp = Precursor::default();
    tmp.activation_methods.insert(ActivationMethod::Hcd);
    tmp.activation_methods.insert(ActivationMethod::Cid);
    tmp.activation_methods.insert(ActivationMethod::Etd);
    assert_eq!(
        tmp.activation_methods_as_short_names(),
        ["CID", "ETD", "HCD"]
    );
    assert!(
        Precursor::default()
            .activation_methods_as_names()
            .is_empty()
    );
}

fn p1() -> Precursor {
    let mut p1 = Precursor {
        mz: 500.5,
        intensity: 1000.0,
        charge: 2,
        activation_energy: 35.0,
        isolation_window_lower_offset: 1.0,
        isolation_window_upper_offset: 1.0,
        drift_time: Some(10.5),
        drift_time_unit: DriftTimeUnit::Millisecond,
        drift_window_lower_offset: 0.5,
        drift_window_upper_offset: 0.5,
        ..Precursor::default()
    };
    p1.activation_methods.insert(ActivationMethod::Cid);
    p1.possible_charge_states.push(2);
    p1.possible_charge_states.push(3);
    p1.cv_terms
        .metadata
        .insert("test_key".to_string(), "test_value".into());
    p1
}

// L453 [EXTRA] std::hash<Precursor>
#[test]
fn section_extra_hash() {
    let p1 = p1();
    let p2 = p1.clone();
    assert_eq!(hash_of(&p1), hash_of(&p2));

    let p3 = Precursor {
        mz: 600.6,
        charge: 3,
        ..Precursor::default()
    };
    assert_ne!(hash_of(&p1), hash_of(&p3));

    // unordered_set / unordered_map: `Precursor` holds floating values and is
    // not `Eq`, so the container semantics - bucket by hash, then `==` - are
    // reproduced by hand.
    let mut buckets: HashMap<u64, Vec<(Precursor, &str)>> = HashMap::new();
    let mut insert = |precursor: &Precursor, value: &'static str| {
        let bucket = buckets.entry(hash_of(precursor)).or_default();
        match bucket.iter_mut().find(|(key, _)| key == precursor) {
            Some(slot) => slot.1 = value,
            None => bucket.push((precursor.clone(), value)),
        }
    };
    insert(&p1, "first");
    insert(&p2, "first");
    insert(&p3, "third");
    assert_eq!(buckets.values().map(Vec::len).sum::<usize>(), 2);
    let lookup = |precursor: &Precursor| {
        buckets
            .get(&hash_of(precursor))
            .and_then(|bucket| bucket.iter().find(|(key, _)| key == precursor))
            .map(|(_, value)| *value)
    };
    assert_eq!(lookup(&p1), Some("first"));
    assert_eq!(lookup(&p2), Some("first"));
    assert_eq!(lookup(&p3), Some("third"));

    assert_eq!(hash_of(&p1), hash_of(&p1));
    let mut p4 = p1.clone();
    let original_hash = hash_of(&p4);
    p4.charge = 5;
    assert_ne!(original_hash, hash_of(&p4));
}

/// Derived: the hash follows equality, which identifies the two zero signs.
#[test]
fn hash_identifies_signed_zeros_as_equality_does() {
    let a = Precursor {
        activation_energy: 0.0,
        ..Precursor::default()
    };
    let b = Precursor {
        activation_energy: -0.0,
        intensity: -0.0,
        ..Precursor::default()
    };
    assert_eq!(a, b);
    assert_eq!(hash_of(&a), hash_of(&b));
}
