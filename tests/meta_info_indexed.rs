// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of the registry and `UInt`-index sections of `MetaInfo_test.cpp` and
//! `MetaInfoInterface_test.cpp` (core SDK bc9cc12), cited by source line. The
//! name-keyed sections are covered by `tests/metadata.rs` and
//! `tests/record_metadata.rs`.
//!
//! The source tests share one process-wide registry across their sections, so
//! later sections see names an earlier one registered ("testname" at 1024).
//! The registry here is caller-owned; each test builds the state its section
//! starts from, which is stated where it matters. Literals are tier 3.

use openms::identification::ProteinHit;
use openms::metadata::meta_info_indexed::IndexedMetaInfo;
use openms::metadata::meta_info_registry::MetaInfoRegistry;
use openms::metadata::{MetaInfo, MetaValue};

fn float(value: f32) -> MetaValue {
    MetaValue::try_from(f64::from(value)).unwrap()
}

/// The class tests' shared `mi` after their getValue sections and the two
/// statements that follow them (MetaInfo_test.cpp:56-79,
/// MetaInfoInterface_test.cpp:54-73): "testname" registered at 1024 by the
/// registry section, index 1024 set, "testname" set by name, then
/// `setValue("cluster_id", 4711.12f)` and `setValue(2, 4712.12f)`.
fn shared(registry: &mut MetaInfoRegistry) -> MetaInfo {
    assert_eq!(
        registry
            .register_name("testname", "testdesc", "testunit")
            .unwrap(),
        1024
    );
    let mut mi = MetaInfo::new();
    mi.set_value_at(registry, 1024, "testtesttest".into())
        .unwrap();
    mi.set_value_registered(registry, "testname", "testtesttest2".into())
        .unwrap();
    mi.set_value_registered(registry, "cluster_id", float(4711.12))
        .unwrap();
    mi.set_value_at(registry, 2, float(4712.12)).unwrap();
    mi
}

// MetaInfo_test.cpp L41 static MetaInfoRegistry& registry();
// MetaInfoInterface_test.cpp L39 static MetaInfoRegistry& metaRegistry()
#[test]
fn section_registry() {
    // The source's static registry is the caller's value here: two MetaInfo
    // objects consult the same registry because they are handed the same one.
    let mut registry = MetaInfoRegistry::default();
    registry
        .register_name("testname", "testdesc", "testunit")
        .unwrap();
    assert_eq!(registry.index_of("testname"), Some(1024));
}

// MetaInfo_test.cpp L56 getValue(UInt index, default); MetaInfoInterface_test.cpp L54 getMetaValue(UInt, default)
#[test]
fn section_get_value_by_index() {
    let mut registry = MetaInfoRegistry::default();
    registry
        .register_name("testname", "testdesc", "testunit")
        .unwrap();
    let mut mi = MetaInfo::new();
    mi.set_value_at(&registry, 1024, "testtesttest".into())
        .unwrap();
    assert_eq!(
        mi.value_at(&registry, 1024),
        Some(&MetaValue::from("testtesttest"))
    );
    // DataValue::EMPTY and the explicit default.
    assert_eq!(mi.value_at(&registry, 1025), None);
    let ten = MetaValue::from(10);
    assert_eq!(mi.value_at(&registry, 1025).unwrap_or(&ten), &ten);
    // The same through a record's MetaInfoInterface payload.
    let mut hit = ProteinHit::default();
    hit.metadata
        .set_value_at(&registry, 1024, "testtesttest".into())
        .unwrap();
    assert_eq!(
        hit.metadata.value_at(&registry, 1024),
        Some(&MetaValue::from("testtesttest"))
    );
}

// MetaInfo_test.cpp L101 getKeys(std::vector<std::string>&) and
// MetaInfoInterface_test.cpp L189 and L266 ([EXTRA]): index order, not name order.
#[test]
fn section_get_keys_names_in_index_order() {
    let mut registry = MetaInfoRegistry::default();
    let mi = shared(&mut registry);
    assert_eq!(
        mi.keys_in_index_order(&registry).unwrap(),
        ["cluster_id", "testname"]
    );
    let mut mi2 = mi.clone();
    assert_eq!(
        mi2.keys_in_index_order(&registry).unwrap(),
        ["cluster_id", "testname"]
    );
    for name in ["a", "d", "x"] {
        mi2.set_value_registered(&mut registry, name, 1.into())
            .unwrap();
    }
    assert_eq!(
        mi2.keys_in_index_order(&registry).unwrap(),
        ["cluster_id", "testname", "a", "d", "x"]
    );
    // Derived: the map itself iterates in name order.
    assert_eq!(
        mi2.keys().map(String::as_str).collect::<Vec<_>>(),
        ["a", "cluster_id", "d", "testname", "x"]
    );
    // MetaInfo_test.cpp:77-78: index 2 overwrote the value set by name.
    assert_eq!(mi.get("cluster_id"), Some(&float(4712.12)));
}

// MetaInfo_test.cpp L135 getKeys(std::vector<UInt>&); MetaInfoInterface_test.cpp L226
#[test]
fn section_get_keys_indices() {
    // "testname" holds 1024 from the registry section, as in both test files.
    let mut registry = MetaInfoRegistry::default();
    registry
        .register_name("testname", "testdesc", "testunit")
        .unwrap();
    let mut mi = MetaInfo::new();
    mi.set_value_registered(&mut registry, "label", "tag".into())
        .unwrap();
    mi.set_value_registered(&mut registry, "icon", "kreis".into())
        .unwrap();
    assert_eq!(mi.keys_as_indices(&registry).unwrap(), [3, 4]);
    for name in ["a", "d", "x"] {
        mi.set_value_registered(&mut registry, name, 1.into())
            .unwrap();
    }
    assert_eq!(
        mi.keys_as_indices(&registry).unwrap(),
        [3, 4, 1025, 1026, 1027]
    );

    // Derived: a name stored without the registry has no index.
    let mut unregistered = MetaInfo::new();
    unregistered.insert("never registered".into(), 1.into());
    assert!(unregistered.keys_as_indices(&registry).is_err());
    assert!(unregistered.keys_in_index_order(&registry).is_err());
}

// MetaInfo_test.cpp L165 exists(UInt); MetaInfoInterface_test.cpp L259 metaValueExists(UInt)
#[test]
fn section_exists_by_index() {
    let mut registry = MetaInfoRegistry::default();
    let mut mi4 = MetaInfo::new();
    assert!(!mi4.exists_at(&registry, 2));
    mi4.set_value_registered(
        &mut registry,
        "cluster_id",
        MetaValue::try_from(4712.1234).unwrap(),
    )
    .unwrap();
    assert!(mi4.exists_at(&registry, 2));
}

// MetaInfo_test.cpp L221 removeValue(UInt); MetaInfoInterface_test.cpp L307 removeMetaValue(UInt)
#[test]
fn section_remove_value_by_index() {
    let registry = MetaInfoRegistry::default();
    let mut i = MetaInfo::new();
    let i2 = MetaInfo::new();
    i.set_value_at(&registry, 1, "bla".into()).unwrap();
    assert!(i != i2);
    i.remove_value_at(&registry, 1);
    assert!(i == i2);
    // Removing a non-existing value works as well.
    i.remove_value_at(&registry, 1234);
    assert!(i == i2);
}

/// Derived: setting an index the registry does not hold is refused and stores
/// nothing, where the source stores the value under the bare index.
#[test]
fn setting_an_unregistered_index_is_refused() {
    let registry = MetaInfoRegistry::default();
    let mut mi = MetaInfo::new();
    assert!(mi.set_value_at(&registry, 4711, "x".into()).is_err());
    assert!(mi.is_empty());
}
