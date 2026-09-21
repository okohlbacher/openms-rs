// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! `DATASTRUCTURES/ToolDescription.h` and `DATASTRUCTURES/ToolInfo.h`.
//!
//! `ToolDescription_test.cpp` (core `bc9cc12`) constructs and destroys a
//! `ToolDescription` and leaves its constructor and assignment sections as
//! `// TODO`; the transcription below keeps its one executable assertion and
//! adds tier-4 cases derived from `ToolDescription.cpp`: the three
//! constructors' field values, the member-wise `operator==` and the key the
//! source's `operator<` compares.

use openms::data_structures::{
    FileMapping, MappingParam, ToolDescription, ToolDescriptionInternal, ToolExternalDetails,
    ToolInfo,
};
use openms::param::{Param, ParamValue};

fn strings(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| (*value).to_owned()).collect()
}

/// `START_SECTION(ToolDescription())`: `TEST_NOT_EQUAL(ptr, null_ptr)`. The
/// default is the source's member initialisers: not internal, empty.
#[test]
fn a_default_description_is_empty_and_external() {
    let description = ToolDescription::default();
    assert!(!description.internal.is_internal);
    assert!(description.internal.name.is_empty());
    assert!(description.internal.category.is_empty());
    assert!(description.internal.types.is_empty());
    assert!(description.external_details.is_empty());
}

/// `ToolDescription(p_name, p_category, p_types)` is the constructor for
/// internal TOPP tools (`ToolDescription.cpp:56-59`).
#[test]
fn the_tool_constructor_marks_the_description_internal() {
    let description = ToolDescription::new("IDFilter", "File Filtering", &strings(&["a", "b"]));
    assert!(description.internal.is_internal);
    assert_eq!(description.internal.name, "IDFilter");
    assert_eq!(description.internal.category, "File Filtering");
    assert_eq!(description.internal.types, strings(&["a", "b"]));
    assert!(description.external_details.is_empty());
    // The types default to an empty list.
    assert!(
        ToolDescription::new("x", "y", &[])
            .internal
            .types
            .is_empty()
    );
}

/// The two `ToolDescriptionInternal` constructors (`ToolDescription.cpp:20-34`):
/// the short one leaves the description external and without a category.
#[test]
fn the_internal_constructors_set_the_given_fields() {
    let full = ToolDescriptionInternal::new(true, "name", "category", &strings(&["t"]));
    assert!(full.is_internal);
    assert_eq!(full.category, "category");
    let short = ToolDescriptionInternal::with_types("name", &strings(&["t"]));
    assert!(!short.is_internal);
    assert_eq!(short.name, "name");
    assert!(short.category.is_empty());
    assert_eq!(short.types, strings(&["t"]));
}

/// `operator==` compares all four fields (`ToolDescription.cpp:36-45`).
#[test]
fn equality_compares_every_registry_field() {
    let base = ToolDescriptionInternal::new(true, "n", "c", &strings(&["t"]));
    assert_eq!(base, base.clone());
    for other in [
        ToolDescriptionInternal::new(false, "n", "c", &strings(&["t"])),
        ToolDescriptionInternal::new(true, "m", "c", &strings(&["t"])),
        ToolDescriptionInternal::new(true, "n", "d", &strings(&["t"])),
        ToolDescriptionInternal::new(true, "n", "c", &strings(&["u"])),
    ] {
        assert_ne!(base, other);
    }
}

/// `operator<` compares `name + "." + types joined by ","`
/// (`ToolDescription.cpp:47-53`): the category and the internal flag play no
/// part, so two descriptions can be unequal and unordered.
#[test]
fn the_source_order_compares_name_and_types_only() {
    let a = ToolDescriptionInternal::new(true, "A", "x", &strings(&["b", "c"]));
    let b = ToolDescriptionInternal::new(false, "A", "y", &strings(&["b", "c"]));
    assert_eq!(a.sort_key(), "A.b,c");
    assert!(!a.source_less(&b) && !b.source_less(&a));
    assert_ne!(a, b);
    let later = ToolDescriptionInternal::with_types("A", &strings(&["d"]));
    assert!(a.source_less(&later));
    assert!(!later.source_less(&a));
    // Byte order, as std::string's operator<: '.' sorts before letters.
    let shorter = ToolDescriptionInternal::with_types("A", &[]);
    let longer = ToolDescriptionInternal::with_types("AB", &[]);
    assert!(shorter.source_less(&longer));
    assert!(!a.source_less(&a));
}

/// The external details and mappings are plain records; the mapping table is
/// ordered by id as the source's `std::map<Int, String>`.
#[test]
fn external_details_keep_their_mappings_in_id_order() {
    let mut table = MappingParam::default();
    table.mapping.insert(3, "%3".into());
    table.mapping.insert(1, "%1".into());
    table.post_moves.push(FileMapping {
        location: "%TMP/out".into(),
        target: "out".into(),
    });
    let mut param = Param::new();
    param
        .set_value("in", ParamValue::String(String::new()), "", &[])
        .unwrap();
    let details = ToolExternalDetails {
        path: "percolator".into(),
        tr_table: table,
        param,
        ..ToolExternalDetails::default()
    };
    let keys: Vec<i32> = details.tr_table.mapping.keys().copied().collect();
    assert_eq!(keys, vec![1, 3]);
    let mut description = ToolDescription::default();
    description.external_details.push(details.clone());
    assert_eq!(description.external_details[0], details);
}

/// `ToolInfo` is the six-field record `TOPPBase` hands the writers
/// (`TOPPBase.cpp:2592-2598`).
#[test]
fn tool_info_holds_the_writer_metadata() {
    let info = ToolInfo {
        version: "1.0.0".into(),
        name: "BaselineFilter".into(),
        docurl: "http://www.openms.de/doxygen/release/4.0.0/html/TOPP_BaselineFilter.html".into(),
        category: "Spectrum Processing: Peak Smoothing and Normalization".into(),
        description: "Removes the baseline from profile spectra using a top-hat filter.".into(),
        citations: strings(&["10.1038/s41592-024-02197-7"]),
    };
    assert_eq!(info.clone(), info);
    assert_eq!(ToolInfo::default().citations.len(), 0);
}
