// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `QPXValueValidation` (`FORMAT/QPXValueValidation.h`).
//!
//! Every section of the upstream class test `QPXValueValidation_test.cpp`
//! (core `bc9cc12`) is carried over; its tables come from
//! `QPXFile::exportPSMsToQPXArrow` and `ConsensusMapArrowExport`, so here
//! they are built column by column with the same key values.
//!
//! The exact diagnostics are the executed Release build's:
//! [`diagnostics_match_the_release_build`] checks every case against the
//! `validate` records of `tests/data/qpx/cpp_release_records.tsv`, which the
//! oracle driver `F1c-qpx` (Linux x86_64 Release install
//! `openms4-release-bc9cc12-c19e494-174b576`, ibminode06, run twice) wrote
//! for the same tables built with Arrow C++.
#![cfg(feature = "arrow")]

use arrow_array::builder::{
    Float32Builder, Int32Builder, ListBuilder, StringBuilder, StructBuilder,
};
use arrow_array::{
    ArrayRef, BooleanArray, Float32Array, Int16Array, Int64Array, RecordBatch, StringArray,
    new_null_array,
};
use arrow_schema::{DataType, Field, Fields, Schema, SchemaRef};
use openms::Error;
use openms::format::arrow_io_helpers::{self as aio, ParquetWriteConfig};
use openms::format::arrow_schema_registry::{QPXFeatureSchema, QPXPSMSchema, QPXPgSchema};
use openms::format::qpx_identity;
use openms::format::qpx_value_validation::{QpxValueValidation, ValidationResult, View};
use std::collections::HashMap;
use std::sync::Arc;

// ---------------------------------------------------------------------------
// Table construction
// ---------------------------------------------------------------------------

/// A table of `schema` with the given columns; every other column is all
/// null. Panics when a non-nullable column is left out.
fn table(schema: &SchemaRef, rows: usize, columns: Vec<(&str, ArrayRef)>) -> RecordBatch {
    let mut given: HashMap<&str, ArrayRef> = columns.into_iter().collect();
    let arrays: Vec<ArrayRef> = schema
        .fields()
        .iter()
        .map(|field| match given.remove(field.name().as_str()) {
            Some(array) => array,
            None => {
                assert!(
                    field.is_nullable(),
                    "missing required column {}",
                    field.name()
                );
                new_null_array(field.data_type(), rows)
            }
        })
        .collect();
    assert!(given.is_empty(), "unknown columns {:?}", given.keys());
    RecordBatch::try_new(schema.clone(), arrays).unwrap()
}

/// Swap one column for `values`, keeping the schema (the upstream test's
/// `replaceColumn`).
fn replace_column(table: &RecordBatch, name: &str, values: ArrayRef) -> RecordBatch {
    let index = table.schema().index_of(name).unwrap();
    let mut columns = table.columns().to_vec();
    columns[index] = values;
    RecordBatch::try_new(table.schema(), columns).unwrap()
}

fn strings(values: &[&str]) -> ArrayRef {
    Arc::new(StringArray::from(values.to_vec()))
}

fn floats(values: &[Option<f32>]) -> ArrayRef {
    Arc::new(Float32Array::from(values.to_vec()))
}

fn int32_lists(field: &Field, rows: &[Option<Vec<Option<i32>>>]) -> ArrayRef {
    let DataType::List(item) = field.data_type() else {
        panic!("not a list")
    };
    let mut builder = ListBuilder::new(Int32Builder::new()).with_field(item.clone());
    for row in rows {
        match row {
            Some(values) => {
                for value in values {
                    builder.values().append_option(*value);
                }
                builder.append(true);
            }
            None => builder.append(false),
        }
    }
    Arc::new(builder.finish())
}

fn string_lists(field: &Field, rows: &[Option<Vec<Option<&str>>>]) -> ArrayRef {
    let DataType::List(item) = field.data_type() else {
        panic!("not a list")
    };
    let mut builder = ListBuilder::new(StringBuilder::new()).with_field(item.clone());
    for row in rows {
        match row {
            Some(values) => {
                for value in values {
                    builder.values().append_option(*value);
                }
                builder.append(true);
            }
            None => builder.append(false),
        }
    }
    Arc::new(builder.finish())
}

fn field_of(schema: &SchemaRef, name: &str) -> Field {
    schema.field_with_name(name).unwrap().clone()
}

fn struct_fields(field: &Field) -> Fields {
    match field.data_type() {
        DataType::List(item) => match item.data_type() {
            DataType::Struct(fields) => fields.clone(),
            other => panic!("not a list of structs: {other}"),
        },
        other => panic!("not a list: {other}"),
    }
}

/// One `intensities` entry: `(label, intensity)`; `None` is a null entry.
type Intensity<'a> = Option<(Option<&'a str>, Option<f32>)>;

fn intensities_column(field: &Field, rows: &[Option<Vec<Intensity>>]) -> ArrayRef {
    let DataType::List(item) = field.data_type() else {
        panic!("not a list")
    };
    let fields = struct_fields(field);
    let mut builder =
        ListBuilder::new(StructBuilder::from_fields(fields, 0)).with_field(item.clone());
    for row in rows {
        match row {
            Some(entries) => {
                for entry in entries {
                    let values = builder.values();
                    let (label, intensity) = entry.unwrap_or((Some(""), Some(0.0)));
                    values
                        .field_builder::<StringBuilder>(0)
                        .unwrap()
                        .append_value(label.unwrap_or(""));
                    values
                        .field_builder::<Float32Builder>(1)
                        .unwrap()
                        .append_option(intensity);
                    values.append(entry.is_some());
                }
                builder.append(true);
            }
            None => builder.append(false),
        }
    }
    Arc::new(builder.finish())
}

/// One `additional_intensities` entry: `(label, [(name, value)])`.
type Additional<'a> = (&'a str, Vec<(&'a str, Option<f32>)>);

fn additional_intensities_column(field: &Field, rows: &[Option<Vec<Additional>>]) -> ArrayRef {
    let DataType::List(item) = field.data_type() else {
        panic!("not a list")
    };
    let entry_fields = struct_fields(field);
    let inner_field = entry_fields[1].clone();
    let inner_fields = struct_fields(&inner_field);
    let DataType::List(inner_item) = inner_field.data_type() else {
        panic!("not a list")
    };
    let inner_builder = ListBuilder::new(StructBuilder::from_fields(inner_fields, 0))
        .with_field(inner_item.clone());
    let entry_builder = StructBuilder::new(
        entry_fields,
        vec![Box::new(StringBuilder::new()), Box::new(inner_builder)],
    );
    let mut builder = ListBuilder::new(entry_builder).with_field(item.clone());
    for row in rows {
        match row {
            Some(entries) => {
                for (label, values) in entry_iter(entries) {
                    let entry = builder.values();
                    entry
                        .field_builder::<StringBuilder>(0)
                        .unwrap()
                        .append_value(label);
                    let inner = entry
                        .field_builder::<ListBuilder<StructBuilder>>(1)
                        .unwrap();
                    for (name, value) in values {
                        let pair = inner.values();
                        pair.field_builder::<StringBuilder>(0)
                            .unwrap()
                            .append_value(name);
                        pair.field_builder::<Float32Builder>(1)
                            .unwrap()
                            .append_option(*value);
                        pair.append(true);
                    }
                    inner.append(true);
                    entry.append(true);
                }
                builder.append(true);
            }
            None => builder.append(false),
        }
    }
    Arc::new(builder.finish())
}

fn entry_iter<'a>(
    entries: &'a [Additional<'a>],
) -> impl Iterator<Item = (&'a str, &'a Vec<(&'a str, Option<f32>)>)> {
    entries.iter().map(|(label, values)| (*label, values))
}

/// One identification of the psm view: `(peptidoform, charge, run, scan)`;
/// a `None` scan is an empty scan list (no spectrum reference).
type Psm<'a> = (&'a str, i16, &'a str, Option<Vec<i32>>);

/// The psm view of `rows`, with `psm_id` derived as `QPXFile` derives it.
fn psm_table(rows: &[Psm]) -> RecordBatch {
    let schema = QPXPSMSchema::schema();
    let n = rows.len();
    let ids: Vec<i64> = rows
        .iter()
        .map(|(peptidoform, charge, run, scan)| {
            qpx_identity::psm_id(
                run,
                scan.as_deref().unwrap_or(&[]),
                peptidoform,
                i64::from(*charge),
            )
        })
        .collect();
    let scans: Vec<Option<Vec<Option<i32>>>> = rows
        .iter()
        .map(|(_, _, _, scan)| {
            Some(
                scan.clone()
                    .unwrap_or_default()
                    .into_iter()
                    .map(Some)
                    .collect(),
            )
        })
        .collect();
    let peptidoforms: Vec<&str> = rows.iter().map(|row| row.0).collect();
    table(
        &schema,
        n,
        vec![
            (QPXPSMSchema::PSM_ID, Arc::new(Int64Array::from(ids))),
            (QPXPSMSchema::SEQUENCE, strings(&peptidoforms)),
            (QPXPSMSchema::PEPTIDOFORM, strings(&peptidoforms)),
            (
                QPXPSMSchema::CHARGE,
                Arc::new(Int16Array::from(
                    rows.iter().map(|row| row.1).collect::<Vec<_>>(),
                )),
            ),
            (
                QPXPSMSchema::IS_DECOY,
                Arc::new(BooleanArray::from(vec![false; n])),
            ),
            (QPXPSMSchema::CALCULATED_MZ, floats(&vec![Some(0.0); n])),
            (QPXPSMSchema::OBSERVED_MZ, floats(&vec![Some(0.0); n])),
            (
                QPXPSMSchema::RUN_FILE_NAME,
                strings(&rows.iter().map(|row| row.2).collect::<Vec<_>>()),
            ),
            (
                QPXPSMSchema::SCAN,
                int32_lists(&field_of(&schema, QPXPSMSchema::SCAN), &scans),
            ),
        ],
    )
}

/// The upstream `makePSMTable(n)`: `n` identifications of run1 with scans
/// 1000, 1001, ... alternating PEPTIDER / DFPIANGER at charge 2.
fn make_psm_table(n: usize) -> RecordBatch {
    let rows: Vec<Psm> = (0..n)
        .map(|i| {
            (
                if i % 2 == 0 { "PEPTIDER" } else { "DFPIANGER" },
                2,
                "run1",
                Some(vec![1000 + i as i32]),
            )
        })
        .collect();
    psm_table(&rows)
}

/// One feature row: `(peptidoform, charge, run, rt, observed_mz)`.
type FeatureRow<'a> = (&'a str, i16, &'a str, Option<f32>, f32);

fn feature_table(
    rows: &[FeatureRow],
    intensities: Option<Vec<Option<Vec<Intensity>>>>,
) -> RecordBatch {
    let schema = QPXFeatureSchema::schema();
    let n = rows.len();
    let ids: Vec<i64> = rows
        .iter()
        .map(|(peptidoform, charge, run, rt, mz)| {
            qpx_identity::feature_id(run, peptidoform, i64::from(*charge), *rt, &[], *mz)
        })
        .collect();
    let peptidoforms: Vec<&str> = rows.iter().map(|row| row.0).collect();
    let intensities = intensities.unwrap_or_else(|| {
        (0..n)
            .map(|_| Some(vec![Some((Some("LFQ"), Some(1000.0)))]))
            .collect()
    });
    table(
        &schema,
        n,
        vec![
            (
                QPXFeatureSchema::FEATURE_ID,
                Arc::new(Int64Array::from(ids)),
            ),
            (QPXFeatureSchema::SEQUENCE, strings(&peptidoforms)),
            (QPXFeatureSchema::PEPTIDOFORM, strings(&peptidoforms)),
            (
                QPXFeatureSchema::CHARGE,
                Arc::new(Int16Array::from(
                    rows.iter().map(|row| row.1).collect::<Vec<_>>(),
                )),
            ),
            (
                QPXFeatureSchema::IS_DECOY,
                Arc::new(BooleanArray::from(vec![false; n])),
            ),
            (QPXFeatureSchema::CALCULATED_MZ, floats(&vec![Some(0.0); n])),
            (
                QPXFeatureSchema::OBSERVED_MZ,
                floats(&rows.iter().map(|row| Some(row.4)).collect::<Vec<_>>()),
            ),
            (
                QPXFeatureSchema::RUN_FILE_NAME,
                strings(&rows.iter().map(|row| row.2).collect::<Vec<_>>()),
            ),
            (
                QPXFeatureSchema::SCAN,
                int32_lists(
                    &field_of(&schema, QPXFeatureSchema::SCAN),
                    &vec![Some(vec![]); n],
                ),
            ),
            (
                QPXFeatureSchema::RT,
                floats(&rows.iter().map(|row| row.3).collect::<Vec<_>>()),
            ),
            (
                QPXFeatureSchema::INTENSITIES,
                intensities_column(
                    &field_of(&schema, QPXFeatureSchema::INTENSITIES),
                    &intensities,
                ),
            ),
        ],
    )
}

/// One pg row: `(accessions, runs, label, intensity)`; `None` lists are null.
type PgRow<'a> = (
    Option<Vec<Option<&'a str>>>,
    Option<Vec<Option<&'a str>>>,
    Option<&'a str>,
    Option<f32>,
);

fn pg_table(rows: &[PgRow], additional: Option<Vec<Option<Vec<Additional>>>>) -> RecordBatch {
    let schema = QPXPgSchema::schema();
    let n = rows.len();
    let ids: Vec<i64> = rows
        .iter()
        .enumerate()
        .map(|(i, (accessions, runs, label, _))| {
            let accessions: Vec<&str> = accessions.iter().flatten().flatten().copied().collect();
            let runs: Vec<&str> = runs.iter().flatten().flatten().copied().collect();
            // Row-unique so the identity column never masks a natural-key check.
            qpx_identity::pg_id(&accessions, &runs, *label) ^ (i as i64)
        })
        .collect();
    let anchors: Vec<&str> = rows
        .iter()
        .map(|row| {
            row.0
                .as_ref()
                .and_then(|a| a.first().copied().flatten())
                .unwrap_or("")
        })
        .collect();
    let peptides_field = field_of(&schema, QPXPgSchema::PEPTIDES);
    let DataType::List(peptides_item) = peptides_field.data_type() else {
        panic!("peptides is a list")
    };
    let mut peptides = ListBuilder::new(StructBuilder::from_fields(
        struct_fields(&peptides_field),
        0,
    ))
    .with_field(peptides_item.clone());
    for _ in 0..n {
        peptides.append(true);
    }
    let mut columns: Vec<(&str, ArrayRef)> = vec![
        (QPXPgSchema::PG_ID, Arc::new(Int64Array::from(ids))),
        (
            QPXPgSchema::PG_ACCESSIONS,
            string_lists(
                &field_of(&schema, QPXPgSchema::PG_ACCESSIONS),
                &rows.iter().map(|row| row.0.clone()).collect::<Vec<_>>(),
            ),
        ),
        (QPXPgSchema::ANCHOR_PROTEIN, strings(&anchors)),
        (
            QPXPgSchema::GROUPED_RUNS,
            string_lists(
                &field_of(&schema, QPXPgSchema::GROUPED_RUNS),
                &rows.iter().map(|row| row.1.clone()).collect::<Vec<_>>(),
            ),
        ),
        (
            QPXPgSchema::LABEL,
            Arc::new(StringArray::from(
                rows.iter().map(|row| row.2).collect::<Vec<_>>(),
            )),
        ),
        (
            QPXPgSchema::INTENSITY,
            floats(&rows.iter().map(|row| row.3).collect::<Vec<_>>()),
        ),
        (
            QPXPgSchema::IS_DECOY,
            Arc::new(BooleanArray::from(vec![false; n])),
        ),
        (QPXPgSchema::PEPTIDES, Arc::new(peptides.finish())),
    ];
    if let Some(additional) = additional {
        columns.push((
            QPXPgSchema::ADDITIONAL_INTENSITIES,
            additional_intensities_column(
                &field_of(&schema, QPXPgSchema::ADDITIONAL_INTENSITIES),
                &additional,
            ),
        ));
    }
    table(&schema, n, columns)
}

fn pg_row<'a>(
    accessions: &[&'a str],
    runs: &[&'a str],
    label: Option<&'a str>,
    intensity: Option<f32>,
) -> PgRow<'a> {
    (
        Some(accessions.iter().map(|a| Some(*a)).collect()),
        Some(runs.iter().map(|r| Some(*r)).collect()),
        label,
        intensity,
    )
}

fn validate(view: View, table: &RecordBatch) -> ValidationResult {
    QpxValueValidation::new(view).validate(table)
}

// ---------------------------------------------------------------------------
// The upstream class test, section by section
// ---------------------------------------------------------------------------

/// `QPXValueValidation(View view)` / `~QPXValueValidation()`.
#[test]
fn construction() {
    let validator = QpxValueValidation::new(View::Psm);
    assert_eq!(validator.view(), View::Psm);
}

/// `[QPXValueValidation::Result] std::string toString() const`.
#[test]
fn result_to_string_joins_diagnostics() {
    let mut result = ValidationResult::default();
    assert!(result.valid);
    assert_eq!(result.to_string(), "");
    result.valid = false;
    result.errors.push("first".into());
    assert_eq!(result.to_string(), "first");
    result.errors.push("second".into());
    assert_eq!(result.to_string(), "first; second");
}

/// `Result validate(const std::shared_ptr<arrow::Table>& table)`.
#[test]
fn a_complete_table_and_an_empty_one_are_valid() {
    let table = make_psm_table(2);
    assert_eq!(table.num_rows(), 2);
    let result = validate(View::Psm, &table);
    assert!(result.valid, "{result}");
    assert!(result.errors.is_empty());
    // A 0-row table is legal: the streaming writers emit one for empty input.
    assert!(validate(View::Psm, &make_psm_table(0)).valid);
    // (The source's null table has no counterpart: a RecordBatch reference
    // cannot be null.)
}

/// `[EXTRA] psm view: identity-bearing key components must carry a value`.
#[test]
fn psm_key_components_must_carry_a_value() {
    let table = make_psm_table(2);

    let no_scan = psm_table(&[
        ("PEPTIDER", 2, "run1", None),
        ("DFPIANGER", 2, "run1", Some(vec![1001])),
    ]);
    let result = validate(View::Psm, &no_scan);
    assert!(!result.valid);
    assert_eq!(
        result.errors,
        ["row 0 has an empty 'scan' primary-key value"]
    );

    let no_run = replace_column(&table, QPXPSMSchema::RUN_FILE_NAME, strings(&["", "run1"]));
    let result = validate(View::Psm, &no_run);
    assert!(!result.valid);
    assert_eq!(
        result.errors,
        ["row 0 has an empty or blank 'run_file_name' primary-key value"]
    );

    let blank_run = replace_column(
        &table,
        QPXPSMSchema::RUN_FILE_NAME,
        strings(&["   ", "run1"]),
    );
    assert!(!validate(View::Psm, &blank_run).valid);
    let tab_run = replace_column(
        &table,
        QPXPSMSchema::RUN_FILE_NAME,
        strings(&["\t\x0b\x0c\r\n", "run1"]),
    );
    assert!(!validate(View::Psm, &tab_run).valid);

    let empty_peptidoform = replace_column(
        &table,
        QPXPSMSchema::PEPTIDOFORM,
        strings(&["PEPTIDER", ""]),
    );
    assert_eq!(
        validate(View::Psm, &empty_peptidoform).errors,
        ["row 1 has an empty or blank 'peptidoform' primary-key value"]
    );

    // A null scan component invalidates the key.
    let schema = QPXPSMSchema::schema();
    let null_component = replace_column(
        &table,
        QPXPSMSchema::SCAN,
        int32_lists(
            &field_of(&schema, QPXPSMSchema::SCAN),
            &[Some(vec![Some(1000), None]), Some(vec![Some(1001)])],
        ),
    );
    assert_eq!(
        validate(View::Psm, &null_component).errors,
        ["row 0 has a null scan component"]
    );
}

/// `[EXTRA] psm view: the primary key must be unique`.
#[test]
fn psm_primary_key_must_be_unique() {
    // The IDMapper shape: one identification copied onto two features. The
    // copies share their psm_id too, so both checks fire.
    let duplicate = psm_table(&[
        ("PEPTIDER", 2, "run1", Some(vec![1000])),
        ("PEPTIDER", 2, "run1", Some(vec![1000])),
    ]);
    assert_eq!(duplicate.num_rows(), 2);
    let result = validate(View::Psm, &duplicate);
    assert!(!result.valid);
    assert_eq!(
        result.errors,
        [
            "row 1 repeats the QPX psm identity 'psm_id'",
            "row 1 repeats the QPX psm primary key (peptidoform, charge, run_file_name, scan)",
        ]
    );

    // Scan lists are compared as sorted lists.
    let mut reordered = psm_table(&[
        ("PEPTIDER", 2, "run1", Some(vec![2, 1])),
        ("PEPTIDER", 2, "run1", Some(vec![1, 2])),
    ]);
    // psm_id keeps scan order, so the ids differ and only the natural key fires.
    reordered = replace_column(
        &reordered,
        QPXPSMSchema::PSM_ID,
        Arc::new(Int64Array::from(vec![1, 2])),
    );
    assert_eq!(
        validate(View::Psm, &reordered).errors,
        ["row 1 repeats the QPX psm primary key (peptidoform, charge, run_file_name, scan)"]
    );
}

/// The identity column alone: a genuine collision is named separately.
#[test]
fn identity_collision_without_a_natural_key_duplicate() {
    let table = make_psm_table(2);
    let collided = replace_column(
        &table,
        QPXPSMSchema::PSM_ID,
        Arc::new(Int64Array::from(vec![7, 7])),
    );
    assert_eq!(
        validate(View::Psm, &collided).errors,
        ["row 1 repeats the QPX psm identity 'psm_id'"]
    );
}

/// `[EXTRA] feature view: unmapped features are separated by observed_mz`.
#[test]
fn feature_unmapped_rows_are_separated_by_observed_mz() {
    assert_eq!(3000.0_f64 as f32, 3000.0001_f64 as f32);
    let table = feature_table(
        &[
            ("", 2, "run1", Some(3000.0_f64 as f32), 400.0),
            ("", 2, "run1", Some(3000.0001_f64 as f32), 500.0),
        ],
        None,
    );
    assert_eq!(table.num_rows(), 2);
    let result = validate(View::Feature, &table);
    assert!(result.valid, "{result}");

    let same_mz = replace_column(
        &table,
        QPXFeatureSchema::OBSERVED_MZ,
        floats(&[Some(400.0), Some(400.0)]),
    );
    let same_mz = replace_column(
        &same_mz,
        QPXFeatureSchema::RT,
        floats(&[Some(3000.0), Some(3000.0)]),
    );
    let result = validate(View::Feature, &same_mz);
    assert!(!result.valid);
    assert_eq!(
        result.errors,
        [
            "row 1 repeats the QPX feature primary key (peptidoform, charge, run_file_name, rt, observed_mz)"
        ]
    );

    // +0 and -0 are one key value; a null rt is a value of its own.
    let signed_zero = replace_column(
        &same_mz,
        QPXFeatureSchema::RT,
        floats(&[Some(0.0), Some(-0.0)]),
    );
    assert!(!validate(View::Feature, &signed_zero).valid);
    let null_rt = replace_column(&same_mz, QPXFeatureSchema::RT, floats(&[None, Some(0.0)]));
    assert!(validate(View::Feature, &null_rt).valid);
    let both_null = replace_column(&same_mz, QPXFeatureSchema::RT, floats(&[None, None]));
    assert!(!validate(View::Feature, &both_null).valid);
}

/// `[EXTRA] feature view: non-finite key coordinates are refused`.
#[test]
fn feature_non_finite_key_coordinates_are_refused() {
    let table = feature_table(&[("", 2, "run1", Some(3000.0), 400.0)], None);

    let nan_rt = replace_column(&table, QPXFeatureSchema::RT, floats(&[Some(f32::NAN)]));
    let result = validate(View::Feature, &nan_rt);
    assert!(!result.valid);
    assert_eq!(
        result.errors,
        ["row 0 has a non-finite 'rt' primary-key value"]
    );

    let inf_mz = replace_column(
        &table,
        QPXFeatureSchema::OBSERVED_MZ,
        floats(&[Some(f32::INFINITY)]),
    );
    let result = validate(View::Feature, &inf_mz);
    assert!(!result.valid);
    assert_eq!(
        result.errors,
        ["row 0 has a non-finite 'observed_mz' primary-key value"]
    );
}

/// The feature view's `intensities` entries.
#[test]
fn feature_intensity_entries() {
    let rows = [("PEPTIDER", 2, "run1", Some(10.0), 400.0)];
    let bad = feature_table(
        &rows,
        Some(vec![Some(vec![
            Some((Some("LFQ"), Some(1.0))),
            Some((Some("LFQ"), Some(f32::NAN))),
            Some((Some("light"), Some(2.0))),
            Some((Some(" "), Some(3.0))),
            None,
        ])]),
    );
    assert_eq!(
        validate(View::Feature, &bad).errors,
        [
            "row 0 repeats intensity label 'LFQ'",
            "row 0 has a null or non-finite intensity for label 'LFQ'",
            "row 0 has non-canonical intensity label 'light'",
            "row 0 contains an empty intensities[].label",
            "row 0 contains a null intensities entry",
        ]
    );
    let tmt = feature_table(
        &rows,
        Some(vec![Some(vec![
            Some((Some("TMT126"), Some(1.0))),
            Some((Some("TMT127N"), Some(1.0))),
        ])]),
    );
    assert!(validate(View::Feature, &tmt).valid);
    let null_list = feature_table(&rows, Some(vec![None]));
    assert!(validate(View::Feature, &null_list).valid);
}

/// The pg view: membership, runs, labels and run-disjointness.
#[test]
fn pg_view_contract() {
    let good = pg_table(
        &[
            pg_row(&["P1", "P2"], &["r1", "r2"], Some("LFQ"), Some(10.0)),
            pg_row(&["P1", "P3"], &["r1"], Some("LFQ"), Some(11.0)),
            pg_row(&["P1", "P2"], &["r3"], None, None),
        ],
        None,
    );
    let result = validate(View::ProteinGroup, &good);
    assert!(result.valid, "{result}");

    // Membership is a set: order and repeats do not make rows distinct.
    let same_group = pg_table(
        &[
            pg_row(&["P2", "P1"], &["r1"], Some("LFQ"), Some(1.0)),
            pg_row(&["P1", "P2", "P1"], &["r1"], Some("LFQ"), Some(2.0)),
        ],
        None,
    );
    assert_eq!(
        validate(View::ProteinGroup, &same_group).errors,
        [
            "row 1 repeats the QPX pg identity (pg_accessions, grouped_runs, label)",
            "row 1 reuses run 'r1' in another pg row with the same (pg_accessions, label)",
        ]
    );

    // A run may feed one quantity of a group per label.
    let overlap = pg_table(
        &[
            pg_row(&["P1"], &["r1", "r2"], Some("LFQ"), Some(1.0)),
            pg_row(&["P1"], &["r2", "r3"], Some("LFQ"), Some(2.0)),
        ],
        None,
    );
    assert_eq!(
        validate(View::ProteinGroup, &overlap).errors,
        ["row 1 reuses run 'r2' in another pg row with the same (pg_accessions, label)"]
    );

    let malformed = pg_table(
        &[
            (Some(vec![]), Some(vec![Some("r1")]), Some("LFQ"), Some(1.0)),
            (
                Some(vec![Some("P1"), Some(" ")]),
                Some(vec![]),
                Some("LFQ"),
                None,
            ),
            pg_row(&["P2"], &["r1", "r1"], Some("light"), Some(f32::INFINITY)),
            pg_row(&["P3"], &["r1"], Some(""), Some(1.0)),
            (Some(vec![Some("P4")]), Some(vec![None]), None, Some(1.0)),
        ],
        None,
    );
    assert_eq!(
        validate(View::ProteinGroup, &malformed).errors,
        [
            "row 0 has an empty or blank 'anchor_protein' primary-key value",
            "row 0 has an empty 'pg_accessions' primary-key value",
            "row 1 has an empty member in pg_accessions",
            "row 1 has an empty 'grouped_runs' primary-key value",
            "row 1 must set label and intensity together, or leave both null",
            "row 2 repeats run 'r1' inside grouped_runs",
            "row 2 has non-canonical label 'light'",
            "row 2 has a non-finite intensity",
            "row 3 has an empty 'label' primary-key value",
            "row 4 has an empty member in grouped_runs",
            "row 4 must set label and intensity together, or leave both null",
        ]
    );
}

/// The pg view's `additional_intensities` must agree with the scalar label.
#[test]
fn pg_additional_intensities() {
    let rows = [
        pg_row(&["P1"], &["r1"], Some("LFQ"), Some(1.0)),
        pg_row(&["P2"], &["r1"], Some("LFQ"), Some(1.0)),
    ];
    let good = pg_table(
        &rows,
        Some(vec![
            Some(vec![(
                "LFQ",
                vec![("iBAQ", Some(1.0)), ("top3", Some(2.0))],
            )]),
            None,
        ]),
    );
    let result = validate(View::ProteinGroup, &good);
    assert!(result.valid, "{result}");

    let bad = pg_table(
        &rows,
        Some(vec![
            Some(vec![
                (
                    "TMT126",
                    vec![("iBAQ", Some(1.0)), ("iBAQ", Some(f32::NAN))],
                ),
                ("light", vec![(" ", Some(1.0))]),
            ]),
            Some(vec![("LFQ", vec![]), ("LFQ", vec![])]),
        ]),
    );
    assert_eq!(
        validate(View::ProteinGroup, &bad).errors,
        [
            "row 0 has additional_intensities label 'TMT126' that differs from its pg label",
            "row 0 repeats additional_intensities intensity name 'iBAQ'",
            "row 0 has a null or non-finite additional_intensities intensity value",
            "row 0 has non-canonical additional_intensities label 'light'",
            "row 0 has additional_intensities label 'light' that differs from its pg label",
            "row 0 contains an empty additional_intensities[].intensities[].intensity_name",
            "row 1 repeats additional_intensities label 'LFQ'",
        ]
    );
}

/// `[EXTRA] a validator spans batches, so duplicates across a batch boundary are found`.
#[test]
fn a_validator_spans_batches() {
    let table = make_psm_table(2);
    let mut validator = QpxValueValidation::new(View::Psm);
    assert!(validator.validate(&table).valid);
    let second = validator.validate(&table);
    assert!(!second.valid);
    assert!(second.to_string().contains("primary key"));
    assert_eq!(second.errors.len(), 4);
}

/// `void reset()`.
#[test]
fn reset_forgets_accumulated_keys() {
    let table = make_psm_table(2);
    let mut validator = QpxValueValidation::new(View::Psm);
    assert!(validator.validate(&table).valid);
    assert!(!validator.validate(&table).valid);
    validator.reset();
    assert!(validator.validate(&table).valid);
}

/// `[EXTRA] a rejected batch does not advance the accumulated key state`.
#[test]
fn a_rejected_batch_does_not_advance_state() {
    let good = make_psm_table(2);
    let bad = replace_column(&good, QPXPSMSchema::RUN_FILE_NAME, strings(&["", "run1"]));
    let mut validator = QpxValueValidation::new(View::Psm);
    assert!(!validator.validate(&bad).valid);
    assert!(validator.validate(&good).valid);
}

/// `[EXTRA] a table that does not carry the schema of its view is refused`.
#[test]
fn a_table_of_another_view_is_refused() {
    let psm = make_psm_table(2);
    let result = validate(View::ProteinGroup, &psm);
    assert!(!result.valid);
    assert!(!result.errors.is_empty());
    assert_eq!(
        result.errors[0],
        "schema: Field count mismatch: got 26, expected 22"
    );
    assert!(result.errors.iter().all(|e| e.starts_with("schema: ")));
}

/// `[EXTRA] a multi-chunk column is validated as one logical column`.
#[test]
fn a_concatenated_table_is_one_logical_column() {
    let one = make_psm_table(1);
    let doubled = arrow_select::concat::concat_batches(&one.schema(), [&one, &one]).unwrap();
    assert_eq!(doubled.num_rows(), 2);
    let result = validate(View::Psm, &doubled);
    assert!(!result.valid);
    assert!(result.to_string().contains("primary key"));
}

/// The generic Parquet writer refuses an invalid QPX view table before
/// touching the destination, and writes a valid one.
#[test]
fn write_table_to_parquet_applies_the_value_contract() {
    let dir = std::env::temp_dir().join(format!("qpx_value_validation_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let config = ParquetWriteConfig::default();
    let metadata = |table: &RecordBatch, file_type: &str| {
        let schema = Schema::new(table.schema().fields().clone()).with_metadata(
            [("qpx_version", "1.1"), ("file_type", file_type)]
                .into_iter()
                .map(|(k, v)| (k.to_owned(), v.to_owned()))
                .collect(),
        );
        RecordBatch::try_new(Arc::new(schema), table.columns().to_vec()).unwrap()
    };

    let good = metadata(&make_psm_table(2), "psm_file");
    let path = dir.join("good.parquet");
    aio::write_table_to_parquet(&good, &path, &config).unwrap();
    assert!(path.exists());

    let duplicate = metadata(
        &psm_table(&[
            ("PEPTIDER", 2, "run1", Some(vec![1000])),
            ("PEPTIDER", 2, "run1", Some(vec![1000])),
        ]),
        "psm_file",
    );
    let path = dir.join("duplicate.parquet");
    match aio::write_table_to_parquet(&duplicate, &path, &config) {
        Err(Error::InvalidValue(message)) => assert_eq!(
            message,
            format!(
                "ArrowIOHelpers: refusing invalid QPX psm_file table for {}: row 1 repeats the \
                 QPX psm identity 'psm_id'; row 1 repeats the QPX psm primary key (peptidoform, \
                 charge, run_file_name, scan)",
                path.display()
            )
        ),
        other => panic!("expected a refusal, got {other:?}"),
    }
    assert!(!path.exists());

    // A psm table labelled as a feature file fails the feature schema.
    let mislabelled = metadata(&make_psm_table(1), "feature_file");
    assert!(aio::write_table_to_parquet(&mislabelled, dir.join("m.parquet"), &config).is_err());
    // Other QPX file types are not value-checked.
    let run_file = metadata(&make_psm_table(1), "run_file");
    aio::write_table_to_parquet(&run_file, dir.join("run.parquet"), &config).unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
}

// ---------------------------------------------------------------------------
// Differential against the Release build
// ---------------------------------------------------------------------------

/// The `validate <tag>/...` records: validity and every diagnostic, in order.
fn release(tag: &str) -> (bool, Vec<String>) {
    let text = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/qpx/cpp_release_records.tsv"),
    )
    .unwrap();
    let mut valid = None;
    let mut errors = Vec::new();
    for line in text.lines() {
        let mut parts = line.splitn(3, '\t');
        if parts.next() != Some("validate") {
            continue;
        }
        let key = parts.next().unwrap();
        let value = parts.next().unwrap_or("").to_owned();
        let Some(rest) = key.strip_prefix(tag).and_then(|r| r.strip_prefix('/')) else {
            continue;
        };
        if rest == "valid" {
            valid = Some(value == "true");
        } else if rest.parse::<usize>().is_ok() {
            errors.push(value);
        }
    }
    (
        valid.unwrap_or_else(|| panic!("no record for {tag}")),
        errors,
    )
}

fn assert_release(tag: &str, result: &ValidationResult) {
    assert_eq!((result.valid, result.errors.clone()), release(tag), "{tag}");
}

#[test]
fn diagnostics_match_the_release_build() {
    let t2 = make_psm_table(2);
    let schema = QPXPSMSchema::schema();
    let psm_cases: Vec<(&str, RecordBatch)> = vec![
        ("psm/valid", t2.clone()),
        ("psm/empty_table", make_psm_table(0)),
        (
            "psm/no_scan",
            psm_table(&[
                ("PEPTIDER", 2, "run1", None),
                ("DFPIANGER", 2, "run1", Some(vec![1001])),
            ]),
        ),
        (
            "psm/no_run",
            replace_column(&t2, QPXPSMSchema::RUN_FILE_NAME, strings(&["", "run1"])),
        ),
        (
            "psm/blank_run",
            replace_column(&t2, QPXPSMSchema::RUN_FILE_NAME, strings(&["   ", "run1"])),
        ),
        (
            "psm/empty_peptidoform",
            replace_column(&t2, QPXPSMSchema::PEPTIDOFORM, strings(&["PEPTIDER", ""])),
        ),
        (
            "psm/null_scan_component",
            replace_column(
                &t2,
                QPXPSMSchema::SCAN,
                int32_lists(
                    &field_of(&schema, QPXPSMSchema::SCAN),
                    &[Some(vec![Some(1000), None]), Some(vec![Some(1001)])],
                ),
            ),
        ),
        (
            "psm/duplicate",
            psm_table(&[
                ("PEPTIDER", 2, "run1", Some(vec![1000])),
                ("PEPTIDER", 2, "run1", Some(vec![1000])),
            ]),
        ),
        (
            "psm/reordered_scans",
            replace_column(
                &psm_table(&[
                    ("PEPTIDER", 2, "run1", Some(vec![2, 1])),
                    ("PEPTIDER", 2, "run1", Some(vec![1, 2])),
                ]),
                QPXPSMSchema::PSM_ID,
                Arc::new(Int64Array::from(vec![1, 2])),
            ),
        ),
        (
            "psm/id_collision",
            replace_column(
                &t2,
                QPXPSMSchema::PSM_ID,
                Arc::new(Int64Array::from(vec![7, 7])),
            ),
        ),
    ];
    for (tag, table) in &psm_cases {
        assert_release(tag, &validate(View::Psm, table));
    }
    assert_release("psm/as_pg", &validate(View::ProteinGroup, &t2));
    let mut validator = QpxValueValidation::new(View::Psm);
    assert!(validator.validate(&t2).valid);
    assert_release("psm/second_batch", &validator.validate(&t2));

    let rt1 = 3000.0_f64 as f32;
    let rt2 = 3000.0001_f64 as f32;
    let unmapped = feature_table(
        &[
            ("", 2, "run1", Some(rt1), 400.0),
            ("", 2, "run1", Some(rt2), 500.0),
        ],
        None,
    );
    let same = replace_column(
        &unmapped,
        QPXFeatureSchema::OBSERVED_MZ,
        floats(&[Some(400.0), Some(400.0)]),
    );
    let same = replace_column(
        &same,
        QPXFeatureSchema::RT,
        floats(&[Some(3000.0), Some(3000.0)]),
    );
    let one = feature_table(&[("", 2, "run1", Some(3000.0), 400.0)], None);
    let single = [("PEPTIDER", 2, "run1", Some(10.0), 400.0)];
    let feature_cases: Vec<(&str, RecordBatch)> = vec![
        ("feature/unmapped", unmapped.clone()),
        ("feature/same_mz", same.clone()),
        (
            "feature/signed_zero",
            replace_column(
                &same,
                QPXFeatureSchema::RT,
                floats(&[Some(0.0), Some(-0.0)]),
            ),
        ),
        (
            "feature/null_rt",
            replace_column(&same, QPXFeatureSchema::RT, floats(&[None, Some(0.0)])),
        ),
        (
            "feature/both_null",
            replace_column(&same, QPXFeatureSchema::RT, floats(&[None, None])),
        ),
        (
            "feature/nan_rt",
            replace_column(&one, QPXFeatureSchema::RT, floats(&[Some(f32::NAN)])),
        ),
        (
            "feature/inf_mz",
            replace_column(
                &one,
                QPXFeatureSchema::OBSERVED_MZ,
                floats(&[Some(f32::INFINITY)]),
            ),
        ),
        (
            "feature/bad_intensities",
            feature_table(
                &single,
                Some(vec![Some(vec![
                    Some((Some("LFQ"), Some(1.0))),
                    Some((Some("LFQ"), Some(f32::NAN))),
                    Some((Some("light"), Some(2.0))),
                    Some((Some(" "), Some(3.0))),
                    None,
                ])]),
            ),
        ),
        (
            "feature/tmt",
            feature_table(
                &single,
                Some(vec![Some(vec![
                    Some((Some("TMT126"), Some(1.0))),
                    Some((Some("TMT127N"), Some(1.0))),
                ])]),
            ),
        ),
        (
            "feature/null_list",
            feature_table(&single, Some(vec![None])),
        ),
    ];
    for (tag, table) in &feature_cases {
        assert_release(tag, &validate(View::Feature, table));
    }

    let two = [
        pg_row(&["P1"], &["r1"], Some("LFQ"), Some(1.0)),
        pg_row(&["P2"], &["r1"], Some("LFQ"), Some(1.0)),
    ];
    let pg_cases: Vec<(&str, RecordBatch)> = vec![
        (
            "pg/good",
            pg_table(
                &[
                    pg_row(&["P1", "P2"], &["r1", "r2"], Some("LFQ"), Some(10.0)),
                    pg_row(&["P1", "P3"], &["r1"], Some("LFQ"), Some(11.0)),
                    pg_row(&["P1", "P2"], &["r3"], None, None),
                ],
                None,
            ),
        ),
        (
            "pg/same_group",
            pg_table(
                &[
                    pg_row(&["P2", "P1"], &["r1"], Some("LFQ"), Some(1.0)),
                    pg_row(&["P1", "P2", "P1"], &["r1"], Some("LFQ"), Some(2.0)),
                ],
                None,
            ),
        ),
        (
            "pg/overlap",
            pg_table(
                &[
                    pg_row(&["P1"], &["r1", "r2"], Some("LFQ"), Some(1.0)),
                    pg_row(&["P1"], &["r2", "r3"], Some("LFQ"), Some(2.0)),
                ],
                None,
            ),
        ),
        (
            "pg/malformed",
            pg_table(
                &[
                    (Some(vec![]), Some(vec![Some("r1")]), Some("LFQ"), Some(1.0)),
                    (
                        Some(vec![Some("P1"), Some(" ")]),
                        Some(vec![]),
                        Some("LFQ"),
                        None,
                    ),
                    pg_row(&["P2"], &["r1", "r1"], Some("light"), Some(f32::INFINITY)),
                    pg_row(&["P3"], &["r1"], Some(""), Some(1.0)),
                    (Some(vec![Some("P4")]), Some(vec![None]), None, Some(1.0)),
                ],
                None,
            ),
        ),
        (
            "pg/additional_good",
            pg_table(
                &two,
                Some(vec![
                    Some(vec![(
                        "LFQ",
                        vec![("iBAQ", Some(1.0)), ("top3", Some(2.0))],
                    )]),
                    None,
                ]),
            ),
        ),
        (
            "pg/additional_bad",
            pg_table(
                &two,
                Some(vec![
                    Some(vec![
                        (
                            "TMT126",
                            vec![("iBAQ", Some(1.0)), ("iBAQ", Some(f32::NAN))],
                        ),
                        ("light", vec![(" ", Some(1.0))]),
                    ]),
                    Some(vec![("LFQ", vec![]), ("LFQ", vec![])]),
                ]),
            ),
        ),
    ];
    for (tag, table) in &pg_cases {
        assert_release(tag, &validate(View::ProteinGroup, table));
    }
}
