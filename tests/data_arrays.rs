// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `src/tests/class_tests/openms/source/DataArrays_test.cpp` and
//! `src/tests/class_tests/openms/source/MetaInfoDescription_test.cpp` (core
//! SDK `bc9cc12`), section by section with unchanged literals.
//!
//! The source's `FloatDataArray`, `IntegerDataArray` and `StringDataArray`
//! are this crate's `DataArray<f32>`, `DataArray<i32>` and `DataArray<String>`;
//! their `MetaInfoDescription` base is the array's `name`, `metadata` and
//! `data_processing` fields. `operator<` and friends are `source_lt`,
//! `source_le`, `source_gt` and `source_ge`, which take the meta-info registry
//! the source reaches as a process-wide static.

use std::sync::Arc;

use openms::kernel::DataArray;
use openms::metadata::meta_info_description::{MetaInfoDescription, data_value_text};
use openms::metadata::meta_info_indexed::IndexedMetaInfo;
use openms::metadata::meta_info_registry::MetaInfoRegistry;
use openms::metadata::{DataProcessing, MetaValue, ProcessingAction};

fn array<T>(data: Vec<T>) -> DataArray<T> {
    DataArray::new("", data)
}

/// Asserts the four ordering operators against the three fixtures every
/// `operator<`..`operator>=` section of DataArrays_test.cpp uses: `a1 == a2`,
/// `a1 < a3`.
fn ordering_sections<T: PartialOrd + Clone + std::fmt::Debug>(
    a1: &DataArray<T>,
    a2: &DataArray<T>,
    a3: &DataArray<T>,
) {
    let r = MetaInfoRegistry::default();
    // operator<
    assert!(!a1.source_lt(a2, &r).unwrap());
    assert!(a1.source_lt(a3, &r).unwrap());
    assert!(!a3.source_lt(a1, &r).unwrap());
    // operator<=
    assert!(a1.source_le(a2, &r).unwrap());
    assert!(a1.source_le(a3, &r).unwrap());
    assert!(!a3.source_le(a1, &r).unwrap());
    // operator>
    assert!(!a1.source_gt(a2, &r).unwrap());
    assert!(!a1.source_gt(a3, &r).unwrap());
    assert!(a3.source_gt(a1, &r).unwrap());
    // operator>=
    assert!(a1.source_ge(a2, &r).unwrap());
    assert!(!a1.source_ge(a3, &r).unwrap());
    assert!(a3.source_ge(a1, &r).unwrap());
}

// ===========================================================================
// DataArrays_test.cpp
// ===========================================================================

// START_SECTION((FloatDataArray())) / (([EXTRA] ~FloatDataArray()))
#[test]
fn float_default() {
    let fda: DataArray<f32> = DataArray::default();
    assert_eq!(fda.data.len(), 0);
}

// START_SECTION((FloatDataArray copy constructor))
#[test]
fn float_copy() {
    let mut fda1 = array(vec![1.0_f32, 2.0, 3.0]);
    fda1.name = "test_array".into();
    let fda2 = fda1.clone();
    assert_eq!(fda2.data, [1.0, 2.0, 3.0]);
    assert_eq!(fda2.name, "test_array");
}

// START_SECTION((FloatDataArray operator==)) and ((FloatDataArray operator!=))
#[test]
fn float_equality() {
    let mut fda1 = array(vec![1.0_f32, 2.0, 3.0]);
    let mut fda2 = array(vec![1.0_f32, 2.0, 3.0]);
    let fda3 = array(vec![1.0_f32, 2.0, 4.0]);
    let fda4 = array(vec![1.0_f32, 2.0]);
    assert!(fda1 == fda2);
    assert!(fda1 != fda3);
    assert!(fda1 != fda4);
    fda1.name = "array1".into();
    assert!(fda1 != fda2);
    fda2.name = "array1".into();
    assert!(fda1 == fda2);

    let fda1 = array(vec![1.0_f32, 2.0, 3.0]);
    let fda2 = array(vec![1.0_f32, 2.0, 3.0]);
    assert!(!fda1.ne(&fda2));
    assert!(fda1 != fda3);
}

// START_SECTION((FloatDataArray operator<)) .. ((FloatDataArray operator>=))
#[test]
fn float_ordering() {
    ordering_sections(
        &array(vec![1.0_f32, 2.0, 3.0]),
        &array(vec![1.0_f32, 2.0, 3.0]),
        &array(vec![1.0_f32, 2.0, 4.0]),
    );
}

// START_SECTION((FloatDataArray with NaN))
#[test]
fn float_nan() {
    let fda1 = array(vec![1.0_f32, f32::NAN, 3.0]);
    let fda2 = array(vec![1.0_f32, f32::NAN, 3.0]);
    let _fda3 = array(vec![1.0_f32, 2.0, 3.0]);
    assert!(!fda1.eq(&fda2)); // NaN != NaN
    assert!(fda1 != fda2);
}

// START_SECTION((IntegerDataArray())) / ([EXTRA] ~IntegerDataArray())
#[test]
fn integer_default() {
    let ida: DataArray<i32> = DataArray::default();
    assert_eq!(ida.data.len(), 0);
}

// START_SECTION((IntegerDataArray copy constructor))
#[test]
fn integer_copy() {
    let mut ida1 = array(vec![1_i32, 2, 3]);
    ida1.name = "int_array".into();
    let ida2 = ida1.clone();
    assert_eq!(ida2.data, [1, 2, 3]);
    assert_eq!(ida2.name, "int_array");
}

// START_SECTION((IntegerDataArray operator==)) and ((IntegerDataArray operator!=))
#[test]
fn integer_equality() {
    let ida1 = array(vec![1_i32, 2, 3]);
    let ida2 = array(vec![1_i32, 2, 3]);
    let ida3 = array(vec![1_i32, 2, 4]);
    assert!(ida1 == ida2);
    assert!(!ida1.eq(&ida3));
    assert!(!ida1.ne(&ida2));
    assert!(ida1 != ida3);
}

// START_SECTION((IntegerDataArray operator<)) .. ((IntegerDataArray operator>=))
#[test]
fn integer_ordering() {
    ordering_sections(
        &array(vec![1_i32, 2, 3]),
        &array(vec![1_i32, 2, 3]),
        &array(vec![1_i32, 2, 4]),
    );
}

// START_SECTION((StringDataArray())) / ([EXTRA] ~StringDataArray())
#[test]
fn string_default() {
    let sda: DataArray<String> = DataArray::default();
    assert_eq!(sda.data.len(), 0);
}

fn strings(items: &[&str]) -> DataArray<String> {
    array(items.iter().map(|s| (*s).to_owned()).collect())
}

// START_SECTION((StringDataArray copy constructor))
#[test]
fn string_copy() {
    let mut sda1 = strings(&["apple", "banana", "cherry"]);
    sda1.name = "string_array".into();
    let sda2 = sda1.clone();
    assert_eq!(sda2.data, ["apple", "banana", "cherry"]);
    assert_eq!(sda2.name, "string_array");
}

// START_SECTION((StringDataArray operator==)) and ((StringDataArray operator!=))
#[test]
fn string_equality() {
    let sda1 = strings(&["apple", "banana", "cherry"]);
    let sda2 = strings(&["apple", "banana", "cherry"]);
    let sda3 = strings(&["apple", "banana", "date"]);
    assert!(sda1 == sda2);
    assert!(!sda1.eq(&sda3));
    assert!(!sda1.ne(&sda2));
    assert!(sda1 != sda3);
}

// START_SECTION((StringDataArray operator<)) .. ((StringDataArray operator>=))
#[test]
fn string_ordering() {
    ordering_sections(
        &strings(&["apple", "banana", "cherry"]),
        &strings(&["apple", "banana", "cherry"]),
        &strings(&["apple", "banana", "date"]),
    );
}

// START_SECTION((DataArrays std::vector functionality)): the vector base is
// the public `data` field.
#[test]
fn vector_functionality() {
    let mut fda: DataArray<f32> = DataArray::default();
    fda.data.push(1.0);
    fda.data.push(2.0);
    fda.data.push(3.0);
    assert_eq!(fda.data.len(), 3);
    assert_eq!(fda.data.first(), Some(&1.0));
    assert_eq!(fda.data.last(), Some(&3.0));
    fda.data.clear();
    assert_eq!(fda.data.len(), 0);
    assert!(fda.data.is_empty());
    fda.data.reserve(100);
    assert!(fda.data.capacity() >= 100);
    let ida = array(vec![1_i32, 2, 3, 4, 5]);
    let sum: i32 = ida.data.iter().sum();
    assert_eq!(sum, 15);
}

// START_SECTION((DataArrays MetaInfoDescription functionality))
#[test]
fn meta_info_description_functionality() {
    let fda: DataArray<f32> = DataArray::new("my_float_array", Vec::new());
    assert_eq!(fda.name, "my_float_array");
    let ida: DataArray<i32> = DataArray::new("my_int_array", Vec::new());
    assert_eq!(ida.name, "my_int_array");
    let sda: DataArray<String> = DataArray::new("my_string_array", Vec::new());
    assert_eq!(sda.name, "my_string_array");

    let mut fda1 = array(vec![1.0_f32, 2.0]);
    let mut fda2 = array(vec![1.0_f32, 2.0]);
    assert!(fda1 == fda2);
    fda1.name = "array1".into();
    assert!(fda1 != fda2);
    fda2.name = "array1".into();
    assert!(fda1 == fda2);
    fda1.metadata.insert("label".into(), "test_label".into());
    assert!(fda1 != fda2);
    fda2.metadata.insert("label".into(), "test_label".into());
    assert!(fda1 == fda2);
}

// START_SECTION((DataArrays comparison with mixed metadata))
#[test]
fn comparison_with_mixed_metadata() {
    let r = MetaInfoRegistry::default();
    let mut fda1 = array(vec![1.0_f32, 2.0, 3.0]);
    let mut fda2 = array(vec![1.0_f32, 2.0, 3.0]);
    fda1.name = "A".into();
    fda2.name = "B".into();
    let less_than = fda1.source_lt(&fda2, &r).unwrap();
    let greater_than = fda1.source_gt(&fda2, &r).unwrap();
    assert!(less_than != greater_than);
    assert!(!fda1.eq(&fda2));
    assert!(fda1 != fda2);
    // Derived from MetaInfoDescription.cpp:29-78: the name is consulted only
    // when the meta values differ, so neither order holds here; `>` is the
    // negation of `<=` and therefore true both ways.
    assert!(!less_than);
    assert!(greater_than);
    assert!(!fda2.source_lt(&fda1, &r).unwrap());
    assert!(fda2.source_gt(&fda1, &r).unwrap());
}

// ===========================================================================
// MetaInfoDescription_test.cpp
// ===========================================================================

// START_SECTION((MetaInfoDescription())) / ((~MetaInfoDescription()))
// START_SECTION((const std::string& getName() const))
#[test]
fn description_name_default() {
    let tmp = MetaInfoDescription::new();
    assert_eq!(tmp.name, "");
}

// START_SECTION((void setName(const std::string& name)))
#[test]
fn description_set_name() {
    let mut tmp = MetaInfoDescription::new();
    tmp.name = "name".into();
    assert_eq!(tmp.name, "name");
}

// START_SECTION((const std::vector<DataProcessing>& getDataProcessing() const))
#[test]
fn description_data_processing_default() {
    let tmp = MetaInfoDescription::new();
    assert_eq!(tmp.data_processing.len(), 0);
}

// START_SECTION((void setDataProcessing(...))) and
// ((std::vector<DataProcessing>& getDataProcessing())): the source resizes
// to one null pointer; an Arc cannot be null, so one default record stands in.
#[test]
fn description_set_data_processing() {
    let mut tmp = MetaInfoDescription::new();
    let dummy = vec![Arc::new(DataProcessing::default())];
    tmp.data_processing = dummy;
    assert_eq!(tmp.data_processing.len(), 1);
    let mut tmp = MetaInfoDescription::new();
    tmp.data_processing
        .push(Arc::new(DataProcessing::default()));
    assert_eq!(tmp.data_processing.len(), 1);
}

fn labelled() -> MetaInfoDescription {
    let mut tmp = MetaInfoDescription::new();
    tmp.name = "bla2".into();
    tmp.data_processing
        .push(Arc::new(DataProcessing::default()));
    tmp.metadata.insert("label".into(), "label".into());
    tmp
}

// START_SECTION((MetaInfoDescription(const MetaInfoDescription& source)))
#[test]
fn description_copy() {
    let tmp = labelled();
    let tmp2 = tmp.clone();
    assert_eq!(tmp2.name, "bla2");
    assert_eq!(tmp.data_processing.len(), 1);
    assert_eq!(tmp2.metadata["label"].as_str().unwrap(), "label");
}

// START_SECTION((MetaInfoDescription& operator= (const MetaInfoDescription& source)))
#[test]
fn description_assignment() {
    let tmp = labelled();
    let mut tmp2 = MetaInfoDescription::new();
    assert_eq!(tmp2.name, "");
    tmp2 = tmp.clone();
    assert_eq!(tmp2.name, "bla2");
    assert_eq!(tmp.data_processing.len(), 1);
    assert_eq!(tmp2.metadata["label"].as_str().unwrap(), "label");
    tmp2 = MetaInfoDescription::new();
    assert_eq!(tmp2.name, "");
    assert_eq!(tmp2.data_processing.len(), 0);
    assert!(!tmp2.metadata.contains_key("label"));
}

// START_SECTION((bool operator== (const MetaInfoDescription& rhs) const))
#[test]
fn description_equality() {
    let empty = MetaInfoDescription::new();
    let mut edit = empty.clone();
    assert!(edit == empty);
    edit.name = "bla2".into();
    assert!(!edit.eq(&empty));
    edit = empty.clone();
    edit.metadata.insert("label".into(), "label".into());
    assert!(!edit.eq(&empty));
}

// ---------------------------------------------------------------------------
// The ordering steps of MetaInfoDescription.cpp:29-78, each derived from the
// source, and the DataArray tie-break of DataArrays.h (vector `<` once the
// descriptions are equal).
// ---------------------------------------------------------------------------

#[test]
fn ordering_steps() {
    let mut registry = MetaInfoRegistry::default();
    // 1. fewer data-processing entries first, before anything else.
    let mut a = MetaInfoDescription::new();
    a.name = "z".into();
    let mut b = MetaInfoDescription::new();
    b.data_processing.push(Arc::new(DataProcessing::default()));
    assert!(a.source_lt(&b, &registry).unwrap());
    assert!(!b.source_lt(&a, &registry).unwrap());

    // 2a. differing meta values: a description with values sorts first.
    let mut with = MetaInfoDescription::new();
    with.name = "z".into();
    with.metadata
        .set_value_registered(&mut registry, "x", MetaValue::from(1_i32))
        .unwrap();
    let without = MetaInfoDescription::new();
    assert!(with.source_lt(&without, &registry).unwrap());
    assert!(!without.source_lt(&with, &registry).unwrap());

    // 2b. then the name.
    let mut named = with.clone();
    named.name = "a".into();
    named
        .metadata
        .set_value_registered(&mut registry, "x", MetaValue::from(2_i32))
        .unwrap();
    assert!(named.source_lt(&with, &registry).unwrap());

    // 2c. then the keys by registry index: "later" registered after "x".
    let mut k1 = MetaInfoDescription::new();
    k1.metadata
        .set_value_registered(&mut registry, "later", MetaValue::from(0_i32))
        .unwrap();
    let mut k2 = MetaInfoDescription::new();
    k2.metadata
        .set_value_registered(&mut registry, "x", MetaValue::from(0_i32))
        .unwrap();
    // {x} has the smaller index vector than {later}, although "later" < "x".
    assert!(k2.source_lt(&k1, &registry).unwrap());
    assert!(!k1.source_lt(&k2, &registry).unwrap());

    // 2d. then the values' text: "10" < "9" as strings.
    let mut v10 = MetaInfoDescription::new();
    v10.metadata
        .set_value_registered(&mut registry, "x", MetaValue::from(10_i32))
        .unwrap();
    let mut v9 = MetaInfoDescription::new();
    v9.metadata
        .set_value_registered(&mut registry, "x", MetaValue::from(9_i32))
        .unwrap();
    assert!(v10.source_lt(&v9, &registry).unwrap());

    // 3. otherwise the processing-action counts.
    let mut p1 = MetaInfoDescription::new();
    p1.data_processing.push(Arc::new(DataProcessing::default()));
    let mut p2 = p1.clone();
    let mut processing = DataProcessing::default();
    processing.actions.insert(ProcessingAction::Smoothing);
    p2.data_processing[0] = Arc::new(processing);
    assert!(p1.source_lt(&p2, &registry).unwrap());
    assert!(!p2.source_lt(&p1, &registry).unwrap());

    // A name the registry does not hold is refused when the order needs it.
    let mut unregistered = MetaInfoDescription::new();
    unregistered
        .metadata
        .insert("never registered".into(), MetaValue::from(1_i32));
    let mut other = MetaInfoDescription::new();
    other
        .metadata
        .insert("never registered".into(), MetaValue::from(2_i32));
    assert!(unregistered.source_lt(&other, &registry).is_err());
}

#[test]
fn data_value_text_matches_the_source() {
    assert_eq!(data_value_text(&MetaValue::from("abc")), "abc");
    assert_eq!(data_value_text(&MetaValue::from(-3_i32)), "-3");
    assert_eq!(data_value_text(&MetaValue::from(vec![1_i64, 2])), "[1, 2]");
    assert_eq!(
        data_value_text(&MetaValue::from(vec!["a".to_owned(), "b".to_owned()])),
        "[a, b]"
    );
    assert_eq!(data_value_text(&MetaValue::default()), "");
}

#[test]
fn data_array_ordering_follows_std_lexicographical_compare() {
    let r = MetaInfoRegistry::default();
    // A proper prefix is less.
    assert!(
        array(vec![1_i32, 2])
            .source_lt(&array(vec![1_i32, 2, 0]), &r)
            .unwrap()
    );
    // NaN is neither less nor greater, so the comparison moves on.
    assert!(
        array(vec![f32::NAN, 1.0])
            .source_lt(&array(vec![0.0_f32, 2.0]), &r)
            .unwrap()
    );
    assert!(
        !array(vec![f32::NAN, 2.0])
            .source_lt(&array(vec![0.0_f32, 1.0]), &r)
            .unwrap()
    );
    // Descriptions decide before the values.
    let mut described = array(vec![9_i32]);
    described
        .data_processing
        .push(Arc::new(DataProcessing::default()));
    assert!(array(vec![100_i32]).source_lt(&described, &r).unwrap());
    assert_eq!(described.description().data_processing.len(), 1);
}
