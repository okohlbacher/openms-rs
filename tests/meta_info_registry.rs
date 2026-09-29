// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `MetaInfoRegistry_test.cpp` (core SDK bc9cc12): one test per
//! `START_SECTION`, cited by source line. The class test shares one registry
//! `mir` across its sections; [`populated`] rebuilds that state, so each test
//! stands alone. Every literal is tier 3, transcribed from the class test.

use std::sync::Mutex;

use openms::Error;
use openms::metadata::meta_info_registry::{
    FIRST_REGISTERED_INDEX, MetaInfoRegistry, RESERVED_NAMES,
};

/// The class test's `mir` after its `registerName` section (L44-53).
fn populated() -> MetaInfoRegistry {
    let mut mir = MetaInfoRegistry::default();
    assert_eq!(
        mir.register_name("testname", "this is just a test", "")
            .unwrap(),
        1024
    );
    assert_eq!(
        mir.register_name("retention time", "this is just another test", "sec")
            .unwrap(),
        1025
    );
    assert_eq!(
        mir.register_name("another testname", "i will be set later", "me too")
            .unwrap(),
        1026
    );
    mir
}

// L32 MetaInfoRegistry(), L38 ~MetaInfoRegistry()
#[test]
fn section_constructor() {
    let registry = MetaInfoRegistry::default();
    // Derived from the source constructor: thirteen reserved names, 1..=13.
    assert_eq!(registry.len(), 13);
    for (index, name, description) in RESERVED_NAMES {
        assert_eq!(registry.index_of(name), Some(index));
        assert_eq!(registry.name_at(index).unwrap(), name);
        assert_eq!(registry.description_at(index).unwrap(), description);
        assert_eq!(registry.unit_at(index).unwrap(), "");
    }
    assert_eq!(FIRST_REGISTERED_INDEX, 1024);
}

// L44 UInt registerName(const std::string&, const std::string& = "", const std::string& = "")
#[test]
fn section_register_name() {
    let mut mir = populated();
    // Derived from the source: a second registration returns the first index
    // and leaves the stored description and unit alone.
    assert_eq!(
        mir.register_name("testname", "ignored", "ignored").unwrap(),
        1024
    );
    assert_eq!(mir.description_at(1024).unwrap(), "this is just a test");
    assert_eq!(mir.unit_at(1024).unwrap(), "");
    assert_eq!(mir.register_name("charge", "", "").unwrap(), 13);
}

// L55 setDescription(UInt, ...), L60 setDescription(const std::string&, ...)
#[test]
fn section_set_description() {
    let mut mir = populated();
    mir.set_description_at(1026, "foo").unwrap();
    assert_eq!(mir.description_at(1026).unwrap(), "foo");
    mir.set_description("another testname", "bar").unwrap();
    assert_eq!(mir.description_at(1026).unwrap(), "bar");

    // The two `Exception::InvalidValue` throws the header documents.
    assert!(matches!(
        mir.set_description_at(4711, "x"),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        mir.set_description("unregistered name", "x"),
        Err(Error::InvalidValue(_))
    ));
}

// L65 setUnit(UInt, ...), L70 setUnit(const std::string&, ...)
#[test]
fn section_set_unit() {
    let mut mir = populated();
    mir.set_unit_at(1026, "foo").unwrap();
    assert_eq!(mir.unit_at(1026).unwrap(), "foo");
    mir.set_unit("another testname", "bar").unwrap();
    assert_eq!(mir.unit_at(1026).unwrap(), "bar");

    assert!(matches!(
        mir.set_unit_at(4711, "x"),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        mir.set_unit("unregistered name", "x"),
        Err(Error::InvalidValue(_))
    ));
}

// L75 UInt getIndex(const std::string&) const
#[test]
fn section_get_index() {
    let mir = populated();
    assert_eq!(mir.index_of("testname"), Some(1024));
    assert_eq!(mir.index_of("retention time"), Some(1025));
    assert_eq!(mir.index_of("isotopic_range"), Some(1));
    assert_eq!(mir.index_of("cluster_id"), Some(2));
    // UInt(-1)
    assert_eq!(mir.index_of("unregistered name"), None);
}

// L83 std::string getName(UInt) const
#[test]
fn section_get_name() {
    let mir = populated();
    assert_eq!(mir.name_at(1).unwrap(), "isotopic_range");
    assert_eq!(mir.name_at(2).unwrap(), "cluster_id");
    assert_eq!(mir.name_at(3).unwrap(), "label");
    assert_eq!(mir.name_at(4).unwrap(), "icon");
    assert_eq!(mir.name_at(1024).unwrap(), "testname");
    assert_eq!(mir.name_at(1025).unwrap(), "retention time");
    assert!(matches!(mir.name_at(4711), Err(Error::InvalidValue(_))));
}

// L92 getDescription(UInt), L99 getDescription(const std::string&)
#[test]
fn section_get_description() {
    let mir = populated();
    assert_eq!(mir.description_at(1024).unwrap(), "this is just a test");
    assert_eq!(
        mir.description_at(1025).unwrap(),
        "this is just another test"
    );
    assert_eq!(
        mir.description_at(1).unwrap(),
        "consecutive numbering of the peaks in an isotope pattern. 0 is the monoisotopic peak"
    );
    assert_eq!(
        mir.description_at(2).unwrap(),
        "consecutive numbering of isotope clusters in a spectrum"
    );

    assert_eq!(mir.description("testname").unwrap(), "this is just a test");
    assert_eq!(
        mir.description("retention time").unwrap(),
        "this is just another test"
    );
    assert_eq!(
        mir.description("isotopic_range").unwrap(),
        "consecutive numbering of the peaks in an isotope pattern. 0 is the monoisotopic peak"
    );
    assert_eq!(
        mir.description("cluster_id").unwrap(),
        "consecutive numbering of isotope clusters in a spectrum"
    );
    assert!(matches!(mir.description_at(0), Err(Error::InvalidValue(_))));
    assert!(matches!(
        mir.description("unregistered name"),
        Err(Error::InvalidValue(_))
    ));
}

// L106 getUnit(UInt), L113 getUnit(const std::string&)
#[test]
fn section_get_unit() {
    let mir = populated();
    assert_eq!(mir.unit_at(1024).unwrap(), "");
    assert_eq!(mir.unit_at(1025).unwrap(), "sec");
    assert_eq!(mir.unit_at(1).unwrap(), "");
    assert_eq!(mir.unit_at(2).unwrap(), "");

    assert_eq!(mir.unit("testname").unwrap(), "");
    assert_eq!(mir.unit("retention time").unwrap(), "sec");
    assert_eq!(mir.unit("isotopic_range").unwrap(), "");
    assert_eq!(mir.unit("cluster_id").unwrap(), "");
    assert!(matches!(mir.unit_at(1023), Err(Error::InvalidValue(_))));
    assert!(matches!(
        mir.unit("unregistered name"),
        Err(Error::InvalidValue(_))
    ));
}

fn assert_copied(mir2: &MetaInfoRegistry) {
    assert_eq!(mir2.index_of("testname"), Some(1024));
    assert_eq!(mir2.index_of("retention time"), Some(1025));
    assert_eq!(mir2.name_at(1).unwrap(), "isotopic_range");
    assert_eq!(mir2.name_at(1024).unwrap(), "testname");
    assert_eq!(mir2.name_at(1025).unwrap(), "retention time");
    assert_eq!(mir2.description_at(1024).unwrap(), "this is just a test");
    assert_eq!(
        mir2.description_at(1025).unwrap(),
        "this is just another test"
    );
    assert_eq!(mir2.description("testname").unwrap(), "this is just a test");
    assert_eq!(
        mir2.description("retention time").unwrap(),
        "this is just another test"
    );
    assert_eq!(mir2.unit_at(1024).unwrap(), "");
    assert_eq!(mir2.unit_at(1025).unwrap(), "sec");
    assert_eq!(mir2.unit("testname").unwrap(), "");
    assert_eq!(mir2.unit("retention time").unwrap(), "sec");
}

// L120 MetaInfoRegistry(const MetaInfoRegistry&)
#[test]
fn section_copy_constructor() {
    let mir = populated();
    let mir2 = mir.clone();
    assert_copied(&mir2);
}

// L137 MetaInfoRegistry& operator=(const MetaInfoRegistry&)
#[test]
fn section_assignment() {
    let mir = populated();
    let mut mir2 = MetaInfoRegistry::default();
    assert_ne!(mir2, mir);
    mir2.clone_from(&mir);
    assert_copied(&mir2);
    // Derived: the copy carries the counter, so both assign the same next index.
    let mut original = mir;
    assert_eq!(original.register_name("next", "", "").unwrap(), 1027);
    assert_eq!(mir2.register_name("next", "", "").unwrap(), 1027);
}

// L155 [EXTRA] multithreaded example: 1e6 concurrent registrations of one new
// name each get the same index, 1027, after the three registrations above.
// The source serialises through `omp critical`; here the caller-owned registry
// sits behind a Mutex, which is the same contract.
#[test]
fn section_extra_multithreaded_example() {
    let mir = Mutex::new(populated());
    let nr_iterations: i64 = 1_000_000;
    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
        .clamp(2, 16) as i64;
    let total: i64 = std::thread::scope(|scope| {
        let workers: Vec<_> = (0..threads)
            .map(|t| {
                let mir = &mir;
                scope.spawn(move || {
                    let mut sum = 0i64;
                    let mut k = 1 + t;
                    while k <= nr_iterations {
                        let val = "newValue";
                        let mut registry = mir.lock().unwrap();
                        let index = registry.register_name(val, "", "").unwrap();
                        registry.set_description_at(index, val).unwrap();
                        drop(registry);
                        sum += i64::from(index);
                        k += threads;
                    }
                    sum
                })
            })
            .collect();
        workers.into_iter().map(|w| w.join().unwrap()).sum()
    });
    assert_eq!(total, 1027 * nr_iterations);
}
