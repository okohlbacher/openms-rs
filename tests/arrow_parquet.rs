// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

#![cfg(feature = "arrow")]

//! `ArrowSchemaRegistry`, `ParquetFile` and `ArrowIOHelpers` against the
//! executed C++ Release build.
//!
//! Every expected value here is a record of the executed Linux x86_64 Release
//! install `openms4-release-bc9cc12-c19e494-174b576` (g++ 14.4 `-O2
//! -DNDEBUG`, Arrow C++ 23.0.1), run twice on ibminode06 with identical
//! output (evidence tier 1). The records are
//! `tests/data/arrow_parquet/cpp_release_records.tsv` (sha256
//! `146a2287e1731e742f2ae82349d77a21dedbf59fd9ae3d28fb44877f0c3b408a`) and
//! `cpp_release_records_2.tsv` (sha256
//! `dab95036c4e07ff26ff3c0a4ebaf7fbfbeff09d58f0dedc63b2ff536589a4df6`); the
//! `*.parquet` files beside them were written by the executed library. See
//! `tests/data/arrow_parquet_provenance.json`. Each test names the records
//! it reads; a record a test cannot reproduce by construction says why.
//!
//! `cpp_release_concat_records.tsv` (sha256
//! `fd644358cbaeab4585ebb75a73f539d41dd2763424cfca47471663211b73a513`, run twice on
//! ibminode06, `../oracle/integ-F-concat`) records `arrow::ConcatenateTables`
//! and `concatenateAndWriteToParquet` over tables whose list element fields
//! are spelled differently.

use arrow_array::builder::{Int32Builder, ListBuilder};
use arrow_array::cast::AsArray;
use arrow_array::types::{Int64Type, TimestampMillisecondType};
use arrow_array::{
    Array, ArrayRef, BooleanArray, Float32Array, Float64Array, Int32Array, Int64Array,
    LargeStringArray, ListArray, RecordBatch, StringArray, StructArray, TimestampSecondArray,
    UInt8Array, UInt32Array, UInt64Array,
};
use arrow_schema::{DataType, Field, Fields, Schema, TimeUnit};
use openms::Error;
use openms::format::arrow_io_helpers::{self as aio, ParquetCompression, ParquetWriteConfig};
use openms::format::arrow_schema_registry::validation::{
    self, Mode, cpp_field_string, types_equal,
};
use openms::format::arrow_schema_registry::*;
use openms::format::parquet_file::ParquetFile;
use openms::kernel::features::{ColumnHeader, ConsensusMap};
use openms::metadata::{MetaInfo, MetaValue, MetaValueData};
use parquet::file::reader::{FileReader, SerializedFileReader};
use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};
use std::sync::Arc;

const RECORDS: &str = include_str!("data/arrow_parquet/cpp_release_records.tsv");
const RECORDS_2: &str = include_str!("data/arrow_parquet/cpp_release_records_2.tsv");
const RECORDS_CONCAT: &str = include_str!("data/arrow_parquet/cpp_release_concat_records.tsv");

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/arrow_parquet")
        .join(name)
}

/// Undo the driver's escaping: `\t`, `\n`, `\\` and `\xHH`.
fn unescape(text: &str) -> String {
    let mut out = String::new();
    let mut bytes = Vec::new();
    let raw = text.as_bytes();
    let mut i = 0;
    while i < raw.len() {
        if raw[i] == b'\\' && i + 1 < raw.len() {
            match raw[i + 1] {
                b't' => bytes.push(b'\t'),
                b'n' => bytes.push(b'\n'),
                b'\\' => bytes.push(b'\\'),
                b'x' => {
                    let hex = std::str::from_utf8(&raw[i + 2..i + 4]).unwrap();
                    bytes.push(u8::from_str_radix(hex, 16).unwrap());
                    i += 2;
                }
                other => panic!("unknown escape {other}"),
            }
            i += 2;
        } else {
            bytes.push(raw[i]);
            i += 1;
        }
    }
    out.push_str(&String::from_utf8(bytes).unwrap());
    out
}

fn parse(text: &str) -> Vec<(String, String, String)> {
    text.lines()
        .map(|line| {
            let mut parts = line.splitn(3, '\t');
            let kind = parts.next().unwrap().to_owned();
            let key = unescape(parts.next().unwrap());
            let value = unescape(parts.next().unwrap_or(""));
            (kind, key, value)
        })
        .collect()
}

/// Every value recorded under `kind` and `key`, in order.
fn all(text: &str, kind: &str, key: &str) -> Vec<String> {
    parse(text)
        .into_iter()
        .filter(|(k, key2, _)| k == kind && key2 == key)
        .map(|(_, _, value)| value)
        .collect()
}

/// The first record of part 1 under `kind` and `key`.
fn rec(kind: &str, key: &str) -> String {
    all(RECORDS, kind, key)
        .into_iter()
        .next()
        .unwrap_or_else(|| panic!("no record {kind}\t{key}"))
}

/// The first record of part 2 under `kind` and `key`.
fn rec2(kind: &str, key: &str) -> String {
    all(RECORDS_2, kind, key)
        .into_iter()
        .next()
        .unwrap_or_else(|| panic!("no record {kind}\t{key}"))
}

/// Assert that `result` fails the way a `throws <Name>: ...` record says:
/// the same exception kind, and the source's detail text (after the last
/// `; `) contained in the message.
fn assert_throws<T: std::fmt::Debug>(result: openms::Result<T>, record: &str) {
    let rest = record
        .strip_prefix("throws ")
        .unwrap_or_else(|| panic!("record does not throw: {record}"));
    let (name, message) = rest.split_once(": ").unwrap();
    let detail = message.rsplit("; ").next().unwrap();
    let error = result.expect_err(record);
    let text = match (&error, name) {
        (Error::InvalidValue(text), "InvalidValue") => text,
        (Error::MissingInformation(text), "MissingInformation") => text,
        _ => panic!("expected {name}, got {error:?}"),
    };
    assert!(
        text.contains(detail),
        "message {text:?} lacks the source detail {detail:?}"
    );
}

/// C `printf("%a")` as glibc prints a double.
fn hexd(v: f64) -> String {
    if v.is_nan() {
        return if v.is_sign_negative() { "-nan" } else { "nan" }.into();
    }
    let sign = if v.is_sign_negative() { "-" } else { "" };
    if v.is_infinite() {
        return format!("{sign}inf");
    }
    let bits = v.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i64;
    let fraction = bits & ((1u64 << 52) - 1);
    if exponent == 0 && fraction == 0 {
        return format!("{sign}0x0p+0");
    }
    let (lead, e) = if exponent == 0 {
        (0, -1022)
    } else {
        (1, exponent - 1023)
    };
    let digits = format!("{fraction:013x}");
    let digits = digits.trim_end_matches('0');
    let mantissa = if digits.is_empty() {
        format!("{lead}")
    } else {
        format!("{lead}.{digits}")
    };
    let plus = if e >= 0 { "+" } else { "" };
    format!("{sign}0x{mantissa}p{plus}{e}")
}

fn scratch(test: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("openms-arrow-{}-{test}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The driver's `smallTable(rows)`: `id` = 10 i, `name` = `row<i % 3>`,
/// `val` = i + 0.5.
fn small_table(rows: usize) -> RecordBatch {
    let ids: Vec<i64> = (0..rows).map(|i| i as i64 * 10).collect();
    let names: Vec<String> = (0..rows).map(|i| format!("row{}", i % 3)).collect();
    let vals: Vec<f64> = (0..rows).map(|i| i as f64 + 0.5).collect();
    let schema = Arc::new(Schema::new(vec![
        Field::new("id", DataType::Int64, true),
        Field::new("name", DataType::Utf8, true),
        Field::new("val", DataType::Float64, true),
    ]));
    RecordBatch::try_new(
        schema,
        vec![
            Arc::new(Int64Array::from(ids)),
            Arc::new(StringArray::from(names)),
            Arc::new(Float64Array::from(vals)),
        ],
    )
    .unwrap()
}

fn with_metadata(table: &RecordBatch, pairs: &[(&str, &str)]) -> RecordBatch {
    let metadata: HashMap<String, String> = pairs
        .iter()
        .map(|(k, v)| ((*k).to_owned(), (*v).to_owned()))
        .collect();
    let schema = Arc::new(Schema::new(table.schema().fields().clone()).with_metadata(metadata));
    table.clone().with_schema(schema).unwrap()
}

/// The driver's `readSchema(fn)`: row count, each field as Arrow C++ prints
/// it with its chunk count (always one here), and the schema metadata keys.
fn read_schema(path: &Path) -> String {
    let table = ParquetFile::read_table(path).unwrap();
    let mut text = format!("rows={}", table.num_rows());
    for field in table.schema().fields() {
        text.push_str(&format!(" | {} chunks=1", cpp_field_string(field)));
    }
    let schema = table.schema();
    let mut keys: Vec<&String> = schema.metadata().keys().collect();
    keys.sort();
    if !keys.is_empty() {
        let keys: Vec<&str> = keys.iter().map(|k| k.as_str()).collect();
        text.push_str(&format!(" | metadata_keys={}", keys.join(",")));
    }
    text
}

/// The driver's `describeFile(key, fn)`, read from the Parquet footer. The
/// codec is Arrow C++'s `arrow::Compression::type` number: 0 uncompressed,
/// 1 snappy, 2 gzip, 4 zstd, 5 lz4 (which Parquet stores as `LZ4_RAW`).
fn describe_file(path: &Path) -> String {
    use parquet::basic::Compression;
    let reader = SerializedFileReader::new(std::fs::File::open(path).unwrap()).unwrap();
    let metadata = reader.metadata();
    let mut text = format!(
        "rows={} row_groups={}",
        metadata.file_metadata().num_rows(),
        metadata.num_row_groups()
    );
    if metadata.num_row_groups() > 0 {
        let rows: Vec<String> = metadata
            .row_groups()
            .iter()
            .map(|g| g.num_rows().to_string())
            .collect();
        let column = metadata.row_group(0).column(0);
        let codec = match column.compression() {
            Compression::UNCOMPRESSED => "0",
            Compression::SNAPPY => "1",
            Compression::GZIP(_) => "2",
            Compression::ZSTD(_) => "4",
            Compression::LZ4_RAW => "5",
            other => panic!("codec {other:?} has no Arrow C++ LZ4 counterpart"),
        };
        text.push_str(&format!(
            " rg_rows={} codec0={codec} stats0={}",
            rows.join(","),
            column.statistics().is_some()
        ));
    }
    let arrow_schema = metadata
        .file_metadata()
        .key_value_metadata()
        .is_some_and(|kv| kv.iter().any(|entry| entry.key == "ARROW:schema"));
    text.push_str(&format!(" arrow_schema={arrow_schema}"));
    text
}

// ---------------------------------------------------------------------------
// ArrowSchemaRegistry
// ---------------------------------------------------------------------------

macro_rules! schemas {
    ($($name:ident),* $(,)?) => {
        vec![$((stringify!($name), $name::schema())),*]
    };
}

/// Records `schema <Name>` and `field <Name>/<i>`: every column of all 22
/// schemas, as Arrow C++ `Field::ToString` prints it.
#[test]
fn every_schema_matches_the_release_build_field_by_field() {
    let all_schemas = schemas!(
        ProteinSchema,
        ProteinGroupSchema,
        SearchParamsSchema,
        FeatureSchema,
        ConsensusFeatureSchema,
        PSMSchema,
        QPXPSMSchema,
        QPXFeatureSchema,
        QPXPgSchema,
        SpectraLongSchema,
        SpectraSemiWideSchema,
        ChromatogramSchema,
        ChromatogramSemiWideSchema,
        OSWPrecursorSchema,
        OSWTransitionSchema,
        OSWFeaturePrecursorSchema,
        OSWRunSchema,
        OSWFeatureSchema,
        OSWFeatureTransitionSchema,
        XICSchema,
        XIMSchema,
        XIPMSchema,
    );
    let recorded = parse(RECORDS)
        .into_iter()
        .filter(|(k, _, _)| k == "schema")
        .count();
    assert_eq!(all_schemas.len(), recorded);
    let mut fields = 0;
    for (name, schema) in all_schemas {
        assert_eq!(
            schema.fields().len().to_string(),
            rec("schema", name),
            "{name}"
        );
        assert!(schema.metadata().is_empty(), "{name}");
        for (i, field) in schema.fields().iter().enumerate() {
            assert_eq!(
                cpp_field_string(field),
                rec("field", &format!("{name}/{i}")),
                "{name}/{i}"
            );
            fields += 1;
        }
    }
    let recorded_fields = parse(RECORDS)
        .into_iter()
        .filter(|(k, _, _)| k == "field")
        .count();
    assert_eq!(fields, recorded_fields);
}

macro_rules! types {
    ($($schema:ident :: $cpp:ident => $rust:ident),* $(,)?) => {
        vec![$((concat!(stringify!($schema), "::", stringify!($cpp)), $schema::$rust())),*]
    };
}

/// Records `type <Schema>::<fn>`: every nested-type helper.
#[test]
fn every_nested_type_helper_matches_the_release_build() {
    let all_types = types!(
        ProteinSchema::modificationsType => modifications_type,
        ProteinSchema::metavaluesType => metavalues_type,
        ProteinGroupSchema::floatDataType => float_data_type,
        ProteinGroupSchema::stringDataType => string_data_type,
        ProteinGroupSchema::integerDataType => integer_data_type,
        SearchParamsSchema::metavaluesType => metavalues_type,
        FeatureSchema::convexHullType => convex_hull_type,
        FeatureSchema::metavaluesType => metavalues_type,
        ConsensusFeatureSchema::handlesType => handles_type,
        ConsensusFeatureSchema::metavaluesType => metavalues_type,
        PSMSchema::modificationsType => modifications_type,
        PSMSchema::additionalScoresType => additional_scores_type,
        PSMSchema::metavaluesType => metavalues_type,
        PSMSchema::proteinAccessionsType => protein_accessions_type,
        QPXPSMSchema::modificationsType => modifications_type,
        QPXPSMSchema::additionalScoresType => additional_scores_type,
        QPXPSMSchema::cvParamsType => cv_params_type,
        QPXPSMSchema::crossLinksType => cross_links_type,
        QPXFeatureSchema::modificationsType => modifications_type,
        QPXFeatureSchema::additionalScoresType => additional_scores_type,
        QPXFeatureSchema::cvParamsType => cv_params_type,
        QPXFeatureSchema::intensitiesType => intensities_type,
        QPXFeatureSchema::additionalIntensitiesType => additional_intensities_type,
        QPXFeatureSchema::pgAccessionsType => pg_accessions_type,
        QPXFeatureSchema::pgPositionsType => pg_positions_type,
        QPXPgSchema::groupedRunsType => grouped_runs_type,
        QPXPgSchema::additionalIntensitiesType => additional_intensities_type,
        QPXPgSchema::peptidesType => peptides_type,
        QPXPgSchema::peptideCountsType => peptide_counts_type,
        QPXPgSchema::featureCountsType => feature_counts_type,
        QPXPgSchema::additionalScoresType => additional_scores_type,
        QPXPgSchema::cvParamsType => cv_params_type,
    );
    let recorded = parse(RECORDS)
        .into_iter()
        .filter(|(k, _, _)| k == "type")
        .count();
    assert_eq!(all_types.len(), recorded);
    for (name, data_type) in all_types {
        assert_eq!(
            validation::cpp_type_string(&data_type),
            rec("type", name),
            "{name}"
        );
    }
}

/// Records `equals ...`: `DataType::Equals` with its default
/// `check_metadata = false`. The record `list item vs element
/// check_metadata` compares with metadata, which no port function offers.
#[test]
fn type_equality_follows_arrow_cpp() {
    let list = |name: &str, nullable: bool| {
        DataType::List(Arc::new(Field::new(name, DataType::Int32, nullable)))
    };
    let expect = |key: &str| rec("equals", key) == "true";
    assert_eq!(
        types_equal(&list("item", true), &list("element", true)),
        expect("list item vs element")
    );
    assert_eq!(
        types_equal(&list("item", true), &list("item", false)),
        expect("list nullable vs not null")
    );
    let strukt = |name: &str, metadata: bool| {
        let mut field = Field::new(name, DataType::Int32, true);
        if metadata {
            field = field.with_metadata(HashMap::from([("k".to_owned(), "v".to_owned())]));
        }
        DataType::Struct(Fields::from(vec![field]))
    };
    assert_eq!(
        types_equal(&strukt("a", false), &strukt("b", false)),
        expect("struct child name")
    );
    assert_eq!(
        types_equal(&strukt("a", false), &strukt("a", true)),
        expect("struct child metadata")
    );
    assert_eq!(
        types_equal(
            &DataType::Timestamp(TimeUnit::Second, None),
            &DataType::Timestamp(TimeUnit::Millisecond, None)
        ),
        expect("timestamp s vs ms")
    );
}

/// Records `validate ...`: `ArrowSchemaValidation::validate` in both modes,
/// as `valid|error count|toString()`.
#[test]
fn schema_validation_matches_the_release_build() {
    let f = |name: &str, data_type: DataType, nullable: bool| Field::new(name, data_type, nullable);
    let s = |fields: Vec<Field>| Schema::new(fields);
    let list_of = |field: Field| DataType::List(Arc::new(field));
    let ab = s(vec![
        f("a", DataType::Utf8, false),
        f("b", DataType::Float64, true),
    ]);
    let check = |key: &str, actual: Schema, expected: &Schema, mode: Mode| {
        let result = validation::validate_schema(&actual, expected, mode);
        assert_eq!(
            format!("{}|{}|{}", result.valid, result.errors.len(), result),
            rec("validate", key),
            "{key}"
        );
        // The table overload checks the same schema.
        let table = RecordBatch::new_empty(Arc::new(actual));
        assert_eq!(validation::validate(&table, expected, mode), result);
    };
    check("strict exact", ab.clone(), &ab, Mode::Strict);
    check(
        "strict missing",
        s(vec![f("a", DataType::Utf8, false)]),
        &ab,
        Mode::Strict,
    );
    check(
        "strict extra",
        s(vec![
            f("a", DataType::Utf8, false),
            f("b", DataType::Float64, true),
            f("c", DataType::Int32, true),
        ]),
        &ab,
        Mode::Strict,
    );
    check(
        "strict type",
        s(vec![
            f("a", DataType::Utf8, false),
            f("b", DataType::Int32, true),
        ]),
        &ab,
        Mode::Strict,
    );
    check(
        "strict nullability",
        s(vec![
            f("a", DataType::Utf8, true),
            f("b", DataType::Float64, false),
        ]),
        &ab,
        Mode::Strict,
    );
    check(
        "strict name",
        s(vec![
            f("x", DataType::Int32, true),
            f("b", DataType::Float64, true),
        ]),
        &ab,
        Mode::Strict,
    );
    check(
        "strict all three",
        s(vec![f("a", DataType::Int32, true)]),
        &s(vec![f("a", DataType::Utf8, false)]),
        Mode::Strict,
    );
    check(
        "strict list element name",
        s(vec![f(
            "l",
            list_of(Field::new("element", DataType::Int32, true)),
            true,
        )]),
        &s(vec![f(
            "l",
            DataType::new_list(DataType::Int32, true),
            true,
        )]),
        Mode::Strict,
    );
    check(
        "strict nested nullability",
        s(vec![f(
            "l",
            list_of(Field::new("item", DataType::Int32, false)),
            true,
        )]),
        &s(vec![f(
            "l",
            DataType::new_list(DataType::Int32, true),
            true,
        )]),
        Mode::Strict,
    );
    check(
        "strict timestamp ms vs s",
        s(vec![f(
            "t",
            DataType::Timestamp(TimeUnit::Millisecond, None),
            true,
        )]),
        &s(vec![f(
            "t",
            DataType::Timestamp(TimeUnit::Second, None),
            true,
        )]),
        Mode::Strict,
    );
    check(
        "strict timestamp tz",
        s(vec![f(
            "t",
            DataType::Timestamp(TimeUnit::Microsecond, Some("UTC".into())),
            true,
        )]),
        &s(vec![f(
            "t",
            DataType::Timestamp(TimeUnit::Second, None),
            true,
        )]),
        Mode::Strict,
    );
    let struct_list = |data_type: DataType, nullable: bool| {
        DataType::new_list(
            DataType::Struct(Fields::from(vec![Field::new("a", data_type, nullable)])),
            true,
        )
    };
    check(
        "strict nested types",
        s(vec![f("m", struct_list(DataType::Int64, false), true)]),
        &s(vec![f("m", struct_list(DataType::Int32, true), true)]),
        Mode::Strict,
    );
    check(
        "strict binary large_string uint8 float",
        s(vec![
            f("a", DataType::Binary, true),
            f("b", DataType::LargeUtf8, true),
            f("c", DataType::UInt8, true),
            f("d", DataType::Float32, true),
        ]),
        &s(vec![
            f("a", DataType::Utf8, true),
            f("b", DataType::Utf8, true),
            f("c", DataType::Int16, true),
            f("d", DataType::Boolean, true),
        ]),
        Mode::Strict,
    );
    check(
        "subset unknown",
        s(vec![
            f("a", DataType::Utf8, false),
            f("z", DataType::Float64, true),
        ]),
        &ab,
        Mode::Subset,
    );
    check(
        "subset missing",
        s(vec![f("b", DataType::Float64, true)]),
        &ab,
        Mode::Subset,
    );
    check(
        "subset type",
        s(vec![f("b", DataType::Int32, true)]),
        &ab,
        Mode::Subset,
    );
    check(
        "subset nullability",
        s(vec![f("a", DataType::Utf8, true)]),
        &ab,
        Mode::Subset,
    );
    check(
        "subset duplicate expected",
        s(vec![f("a", DataType::Int32, true)]),
        &s(vec![
            f("a", DataType::Utf8, true),
            f("a", DataType::Float64, true),
        ]),
        Mode::Subset,
    );
    check(
        "subset duplicate actual",
        s(vec![
            f("b", DataType::Int32, true),
            f("b", DataType::Float64, true),
        ]),
        &ab,
        Mode::Subset,
    );
    let empty = validation::ValidationResult::default();
    assert_eq!(
        format!("{}|{}", empty.valid, empty),
        rec("validate", "empty toString")
    );
}

// ---------------------------------------------------------------------------
// ParquetFile accessors
// ---------------------------------------------------------------------------

fn string_list(array: &dyn Array, row: usize) -> openms::Result<String> {
    let values = ParquetFile::get_string_list(Some(array), row)?;
    let mut text = values.len().to_string();
    for value in values {
        text.push('|');
        text.push_str(&value);
    }
    Ok(text)
}

/// Records `getStringList ...`: `;`-separated strings split verbatim, list
/// columns element by element, unsupported types refused.
#[test]
fn get_string_list_matches_the_release_build() {
    for s in ["a;;b;", ";", "a", "", ";;", " a ; b "] {
        let array = StringArray::from(vec![s]);
        assert_eq!(
            string_list(&array, 0).unwrap(),
            rec("getStringList", &format!("string '{s}'")),
            "{s:?}"
        );
    }
    let list =
        ListArray::from_iter_primitive::<Int64Type, _, _>(Vec::<Option<Vec<Option<i64>>>>::new());
    drop(list);
    let mut builder = ListBuilder::new(arrow_array::builder::StringBuilder::new());
    builder.values().append_value("x");
    builder.values().append_null();
    builder.values().append_value("y;z");
    builder.append(true);
    builder.append(false);
    builder.append(true);
    let strings = builder.finish();
    assert_eq!(
        string_list(&strings, 0).unwrap(),
        rec("getStringList", "list row0")
    );
    assert_eq!(
        string_list(&strings, 1).unwrap(),
        rec("getStringList", "list row1 null")
    );
    assert_eq!(
        string_list(&strings, 2).unwrap(),
        rec("getStringList", "list row2 empty")
    );
    let mut ints = ListBuilder::new(Int32Builder::new());
    ints.values().append_value(1);
    ints.append(true);
    let ints = ints.finish();
    assert_throws(
        string_list(&ints, 0),
        &rec("getStringList", "list of int32"),
    );
    assert_throws(
        string_list(&Int64Array::from(vec![1]), 0),
        &rec("getStringList", "int64 column"),
    );
    assert_eq!(
        string_list(&LargeStringArray::from(vec!["a;b"]), 0).unwrap(),
        rec("getStringList", "large_string a;b")
    );
}

/// Records `getInt64`, `getDouble`, `getBool` and `getString`: the coercing
/// scalar getters.
#[test]
fn coercing_getters_match_the_release_build() {
    let big = UInt64Array::from(vec![
        9_223_372_036_854_775_808u64,
        9_223_372_036_854_775_807,
    ]);
    assert_throws(
        ParquetFile::get_int64(Some(&big), 0, -1, false),
        &rec("getInt64", "uint64 2^63"),
    );
    assert_eq!(
        ParquetFile::get_int64(Some(&big), 1, -1, false)
            .unwrap()
            .to_string(),
        rec("getInt64", "uint64 2^63-1")
    );
    assert_throws(
        ParquetFile::get_double(Some(&big), 0, -1.0, false),
        &rec("getDouble", "uint64 2^63"),
    );
    assert_throws(
        ParquetFile::get_bool(Some(&big), 0, false, false),
        &rec("getBool", "uint64 2^63"),
    );
    let double = Float64Array::from(vec![1.5]);
    assert_throws(
        ParquetFile::get_int64(Some(&double), 0, -1, false),
        &rec("getInt64", "double"),
    );
    assert_throws(
        ParquetFile::get_bool(Some(&double), 0, false, false),
        &rec("getBool", "double"),
    );
    assert_throws(
        ParquetFile::get_string(Some(&double), 0),
        &rec("getString", "double"),
    );
    let float = Float32Array::from(vec![0.1f32]);
    assert_eq!(
        hexd(ParquetFile::get_double(Some(&float), 0, -1.0, false).unwrap()),
        rec("getDouble", "float 0.1")
    );
    assert_eq!(
        hexd(
            ParquetFile::get_double(Some(&Int64Array::from(vec![i64::MAX])), 0, -1.0, false)
                .unwrap()
        ),
        rec("getDouble", "int64 max")
    );
    assert_throws(
        ParquetFile::get_double(Some(&StringArray::from(vec!["1"])), 0, -1.0, false),
        &rec("getDouble", "string"),
    );
    let small = UInt8Array::from(vec![Some(2), None]);
    assert_eq!(
        ParquetFile::get_bool(Some(&small), 0, false, false)
            .unwrap()
            .to_string(),
        rec("getBool", "uint8 2")
    );
    assert_eq!(
        ParquetFile::get_bool(Some(&small), 1, true, true)
            .unwrap()
            .to_string(),
        rec("getBool", "uint8 null allow")
    );
    assert_throws(
        ParquetFile::get_bool(Some(&small), 1, true, false),
        &rec("getBool", "uint8 null disallow"),
    );
    assert_eq!(
        ParquetFile::get_string(Some(&small), 1).unwrap(),
        rec("getString", "uint8 null")
    );
    assert_throws(
        ParquetFile::get_int64(Some(&small), 1, 7, false),
        &rec("getInt64", "uint8 null disallow"),
    );
    assert_throws(
        ParquetFile::get_double(Some(&small), 1, 7.0, false),
        &rec("getDouble", "uint8 null disallow"),
    );
    assert_throws(
        ParquetFile::get_bool(Some(&StringArray::from(vec!["true"])), 0, false, false),
        &rec("getBool", "string"),
    );
    assert_eq!(
        ParquetFile::get_string(Some(&LargeStringArray::from(vec!["large"])), 0).unwrap(),
        rec("getString", "large_string")
    );
    let binary = arrow_array::BinaryArray::from(vec![&b"bin"[..]]);
    assert_throws(
        ParquetFile::get_string(Some(&binary), 0),
        &rec("getString", "binary"),
    );
}

/// Records `jsonEscape`, `appendOrThrow` and `finishArray`.
#[test]
fn json_escape_and_builder_helpers_match_the_release_build() {
    assert_eq!(
        ParquetFile::json_escape("a\x01\x1f\x7f\u{e9}\"\\\n\r\t\x08\x0c/z"),
        rec("jsonEscape", "controls")
    );
    let failed: Result<(), arrow_schema::ArrowError> = Err(
        arrow_schema::ArrowError::InvalidArgumentError("boom".into()),
    );
    assert_throws(
        ParquetFile::append_or_error(failed, "col"),
        &rec("appendOrThrow", "invalid").replace("'Invalid: boom'", "''"),
    );
    let mut builder = arrow_array::builder::Int64Builder::new();
    builder.append_value(77);
    let array = ParquetFile::finish_array(&mut builder, "id").unwrap();
    assert_eq!(array.len().to_string(), rec("finishArray", "int64"));
}

// ---------------------------------------------------------------------------
// ParquetFile I/O
// ---------------------------------------------------------------------------

/// The driver's five-column table of `writeTable default`.
fn default_table() -> RecordBatch {
    let mut lists = ListBuilder::new(Int32Builder::new());
    lists.values().append_value(1);
    lists.values().append_value(2);
    lists.append(true);
    lists.append(false);
    let lists = lists.finish();
    let times = TimestampSecondArray::from(vec![Some(1_700_000_000), None]);
    let schema = Arc::new(Schema::new(vec![
        Field::new("l", lists.data_type().clone(), true),
        Field::new("t", DataType::Timestamp(TimeUnit::Second, None), true),
        Field::new("u", DataType::UInt32, true),
        Field::new("ls", DataType::LargeUtf8, true),
        Field::new("nn", DataType::Int32, false),
    ]));
    RecordBatch::try_new(
        schema,
        vec![
            Arc::new(lists),
            Arc::new(times),
            Arc::new(UInt32Array::from(vec![4_000_000_000u32, 1])),
            Arc::new(LargeStringArray::from(vec![Some("big"), None])),
            Arc::new(Int32Array::from(vec![5, 6])),
        ],
    )
    .unwrap()
}

/// Records `writeTable default`, `file parquetfile_default`, `rowCount ...`
/// and `readTable garbage/missing`: `ParquetFile::writeTable` stores no Arrow
/// schema, so Arrow-only types come back as their Parquet types.
#[test]
fn parquet_file_write_table_matches_the_release_build() {
    let dir = scratch("pf-default");
    let path = dir.join("parquetfile_default.parquet");
    ParquetFile::write_table(&default_table(), &path).unwrap();
    assert_eq!(read_schema(&path), rec("writeTable", "default"));
    assert_eq!(describe_file(&path), rec("file", "parquetfile_default"));
    // The same reader on the file the executed source wrote.
    assert_eq!(
        describe_file(&data("parquetfile_default.parquet")),
        rec("file", "parquetfile_default")
    );
    assert_eq!(
        read_schema(&data("parquetfile_default.parquet")),
        rec("writeTable", "default")
    );
    let back = ParquetFile::read_table(data("parquetfile_default.parquet")).unwrap();
    let times = back.column(1).as_primitive::<TimestampMillisecondType>();
    assert_eq!(times.value(0), 1_700_000_000_000);
    assert!(times.is_null(1));
    assert_eq!(
        ParquetFile::row_count(&path).unwrap().to_string(),
        rec("rowCount", "default")
    );
    assert_eq!(
        ParquetFile::row_count(dir.join("no_such.parquet"))
            .unwrap()
            .to_string(),
        rec("rowCount", "missing")
    );
    let garbage = dir.join("garbage.bin");
    std::fs::write(&garbage, "not a Parquet file").unwrap();
    // The source lets Parquet's own std::exception escape; the port reports
    // it as an invalid value.
    assert!(rec("rowCount", "garbage").starts_with("throws std: "));
    assert!(matches!(
        ParquetFile::row_count(&garbage),
        Err(Error::InvalidValue(_))
    ));
    assert_throws(
        ParquetFile::read_table(&garbage),
        &rec("readTable", "garbage"),
    );
    let missing = rec("readTable", "missing");
    assert!(missing.contains("Failed to open parquet file"));
    assert!(matches!(
        ParquetFile::read_table(dir.join("no_such.parquet")),
        Err(Error::InvalidValue(message)) if message.contains("Failed to open parquet file")
    ));
}

/// Records `writeTable row_group_size ...`, `file parquetfile_rg2`, `file
/// parquetfile_empty*`, `readTable empty table`, `getColumn ...` and
/// `getOptionalColumn missing`.
#[test]
fn parquet_file_row_groups_and_columns_match_the_release_build() {
    let dir = scratch("pf-rg");
    let rg2 = dir.join("parquetfile_rg2.parquet");
    ParquetFile::write_table_with_row_group_size(&small_table(5), &rg2, 2).unwrap();
    assert_eq!(rec("writeTable", "row_group_size 2"), "ok");
    assert_eq!(describe_file(&rg2), rec("file", "parquetfile_rg2"));
    assert_eq!(
        describe_file(&data("parquetfile_rg2.parquet")),
        rec("file", "parquetfile_rg2")
    );
    let rg0 = dir.join("parquetfile_rg0.parquet");
    assert_throws(
        ParquetFile::write_table_with_row_group_size(&small_table(5), &rg0, 0),
        &rec("writeTable", "row_group_size 0")
            .replace("the value 'Invalid: chunk size per row_group must be greater than 0' was used but is not valid; Failed to write parquet table", "chunk size per row_group must be greater than 0"),
    );
    assert!(!rg0.exists());
    let empty = dir.join("parquetfile_empty.parquet");
    ParquetFile::write_table(&small_table(0), &empty).unwrap();
    assert_eq!(rec("writeTable", "empty table"), "ok");
    assert_eq!(describe_file(&empty), rec("file", "parquetfile_empty"));
    assert_eq!(read_schema(&empty), rec("readTable", "empty table"));
    let back = ParquetFile::read_table(&empty).unwrap();
    assert_eq!(
        ParquetFile::get_column(&back, "id")
            .unwrap()
            .len()
            .to_string(),
        rec("getColumn", "empty table id")
    );
    let empty_rg0 = dir.join("parquetfile_empty_rg0.parquet");
    ParquetFile::write_table_with_row_group_size(&small_table(0), &empty_rg0, 0).unwrap();
    assert_eq!(rec("writeTable", "empty table row_group_size 0"), "ok");
    assert_eq!(
        describe_file(&empty_rg0),
        rec("file", "parquetfile_empty_rg0")
    );
    assert_throws(
        ParquetFile::get_column(&small_table(1), "nope"),
        &rec("getColumn", "missing"),
    );
    assert_eq!(
        ParquetFile::get_optional_column(&small_table(1), "nope")
            .is_none()
            .to_string(),
        rec("getOptionalColumn", "missing")
    );
    // Source defect: ParquetFile::getColumn returns chunk(0) of a two-chunk
    // column, 2 of 5 rows. A RecordBatch column has one chunk, so the port
    // returns every row: the invariant "a column has as many values as the
    // table has rows".
    assert_eq!(rec("getColumn", "two chunks length"), "2 of 5");
    let both = arrow_select::concat::concat_batches(
        &small_table(2).schema(),
        &[small_table(2), small_table(3)],
    )
    .unwrap();
    assert_eq!(
        ParquetFile::get_column(&both, "id").unwrap().len(),
        both.num_rows()
    );
}

// ---------------------------------------------------------------------------
// ArrowIOHelpers writing
// ---------------------------------------------------------------------------

/// Records `writeTableToParquet <codec>`, `file awp_<codec>` and `readTable
/// awp_<codec>`: every codec, the stored Arrow schema, and the C++-written
/// files read back.
#[test]
fn write_table_to_parquet_codecs_match_the_release_build() {
    let dir = scratch("awp-codecs");
    for (name, compression) in [
        ("none", ParquetCompression::None),
        ("snappy", ParquetCompression::Snappy),
        ("gzip", ParquetCompression::Gzip),
        ("lz4", ParquetCompression::Lz4),
        ("zstd", ParquetCompression::Zstd),
    ] {
        let config = ParquetWriteConfig {
            compression,
            ..ParquetWriteConfig::default()
        };
        let file = format!("awp_{name}.parquet");
        let path = dir.join(&file);
        assert_eq!(rec("writeTableToParquet", name), "true");
        aio::write_table_to_parquet(&small_table(5), &path, &config).unwrap();
        let key = format!("awp_{name}");
        assert_eq!(describe_file(&path), rec("file", &key), "{name}");
        assert_eq!(describe_file(&data(&file)), rec("file", &key), "{name}");
        assert_eq!(read_schema(&path), rec("readTable", &key), "{name}");
        assert_eq!(read_schema(&data(&file)), rec("readTable", &key), "{name}");
        // The C++-written file holds the driver's table.
        assert_eq!(
            ParquetFile::read_table(data(&file)).unwrap(),
            small_table(5)
        );
        assert_eq!(ParquetFile::read_table(&path).unwrap(), small_table(5));
    }
    let config = ParquetWriteConfig {
        write_statistics: false,
        ..ParquetWriteConfig::default()
    };
    let path = dir.join("awp_nostats.parquet");
    assert_eq!(rec("writeTableToParquet", "nostats"), "true");
    aio::write_table_to_parquet(&small_table(5), &path, &config).unwrap();
    assert_eq!(describe_file(&path), rec("file", "awp_nostats"));
    assert_eq!(
        describe_file(&data("awp_nostats.parquet")),
        rec("file", "awp_nostats")
    );
}

/// Records `writeTableToParquet zstd level ...` and `gzip level ...`. Their
/// `exists=` suffix was evaluated before the write (argument evaluation
/// order), so only the success flag is evidence.
#[test]
fn write_table_to_parquet_compression_levels_match_the_release_build() {
    let dir = scratch("awp-levels");
    for (compression, name, levels) in [
        (ParquetCompression::Zstd, "zstd", [-1, 0, 1, 22, 23]),
        (ParquetCompression::Gzip, "gzip", [-1, 0, 1, 9, 10]),
    ] {
        for level in levels {
            let config = ParquetWriteConfig {
                compression,
                compression_level: level,
                ..ParquetWriteConfig::default()
            };
            let path = dir.join(format!("{name}_{level}.parquet"));
            let record = rec("writeTableToParquet", &format!("{name} level {level}"));
            let expected = record.split(' ').next().unwrap() == "true";
            let result = aio::write_table_to_parquet(&small_table(5), &path, &config);
            assert_eq!(result.is_ok(), expected, "{name} level {level}: {result:?}");
            assert_eq!(path.exists(), expected, "{name} level {level}");
            if expected {
                assert_eq!(ParquetFile::read_table(&path).unwrap(), small_table(5));
            }
        }
    }
}

/// Records `writeTableToParquet row_group_size ...`, `1100000 rows default`,
/// `data_page_size 0`, `bad directory` and `file awp_rg2`/`awp_big`.
#[test]
fn write_table_to_parquet_row_groups_and_pages_match_the_release_build() {
    let dir = scratch("awp-rg");
    let zero = ParquetWriteConfig {
        row_group_size: 0,
        ..ParquetWriteConfig::default()
    };
    assert!(rec("writeTableToParquet", "row_group_size 0").starts_with("false"));
    let path = dir.join("awp_rg0.parquet");
    assert!(aio::write_table_to_parquet(&small_table(5), &path, &zero).is_err());
    assert!(!path.exists());
    let two = ParquetWriteConfig {
        row_group_size: 2,
        ..ParquetWriteConfig::default()
    };
    let path = dir.join("awp_rg2.parquet");
    assert_eq!(rec("writeTableToParquet", "row_group_size 2"), "true");
    aio::write_table_to_parquet(&small_table(5), &path, &two).unwrap();
    assert_eq!(describe_file(&path), rec("file", "awp_rg2"));
    assert_eq!(
        describe_file(&data("awp_rg2.parquet")),
        rec("file", "awp_rg2")
    );
    let big = dir.join("awp_big.parquet");
    assert_eq!(rec("writeTableToParquet", "1100000 rows default"), "true");
    aio::write_table_to_parquet(
        &small_table(1_100_000),
        &big,
        &ParquetWriteConfig::default(),
    )
    .unwrap();
    assert_eq!(describe_file(&big), rec("file", "awp_big"));
    std::fs::remove_file(&big).unwrap();
    let page0 = ParquetWriteConfig {
        data_page_size: 0,
        ..ParquetWriteConfig::default()
    };
    let path = dir.join("awp_page0.parquet");
    assert!(rec("writeTableToParquet", "data_page_size 0").starts_with("true"));
    aio::write_table_to_parquet(&small_table(5), &path, &page0).unwrap();
    assert_eq!(ParquetFile::read_table(&path).unwrap(), small_table(5));
    assert_eq!(rec("writeTableToParquet", "bad directory"), "false");
    assert!(
        aio::write_table_to_parquet(
            &small_table(1),
            dir.join("no/such/dir/x.parquet"),
            &ParquetWriteConfig::default()
        )
        .is_err()
    );
    // The source's `null table` record has no counterpart: a RecordBatch
    // reference cannot be null.
    assert!(rec("writeTableToParquet", "null table").starts_with("false"));
}

/// Records `writeTableToParquet qpx ...`, `custom metadata` and `readTable
/// awp_metadata`.
#[test]
fn write_table_to_parquet_qpx_gate_and_metadata_match_the_release_build() {
    let dir = scratch("awp-qpx");
    let config = ParquetWriteConfig::default();
    let table = small_table(2);
    assert!(rec("writeTableToParquet", "qpx without file_type").starts_with("false"));
    let path = dir.join("awp_qpx1.parquet");
    assert!(matches!(
        aio::write_table_to_parquet(
            &with_metadata(&table, &[("qpx_version", "1.1")]),
            &path,
            &config
        ),
        Err(Error::InvalidValue(_))
    ));
    assert!(!path.exists());
    assert_eq!(rec("writeTableToParquet", "qpx other file_type"), "true");
    let path = dir.join("awp_qpx2.parquet");
    let run_file = with_metadata(&table, &[("qpx_version", "1.1"), ("file_type", "run_file")]);
    aio::write_table_to_parquet(&run_file, &path, &config).unwrap();
    assert_eq!(ParquetFile::read_table(&path).unwrap(), run_file);
    assert_eq!(
        ParquetFile::read_table(data("awp_qpx2.parquet")).unwrap(),
        run_file
    );
    // The source validates a psm_file table with QPXValueValidation and
    // refuses this one: it does not carry the QPX psm schema.
    assert!(rec("writeTableToParquet", "qpx psm_file invalid table").starts_with("false"));
    let path = dir.join("awp_qpx3.parquet");
    let psm = with_metadata(&table, &[("qpx_version", "1.1"), ("file_type", "psm_file")]);
    match aio::write_table_to_parquet(&psm, &path, &config) {
        Err(Error::InvalidValue(message)) => {
            assert!(message.contains("refusing invalid QPX psm_file table"));
            assert!(message.contains("schema: "));
        }
        other => panic!("expected a QPX value refusal, got {other:?}"),
    }
    assert!(!path.exists());
    assert_eq!(rec("writeTableToParquet", "custom metadata"), "true");
    let path = dir.join("awp_metadata.parquet");
    aio::write_table_to_parquet(&with_metadata(&table, &[("x", "y")]), &path, &config).unwrap();
    assert_eq!(read_schema(&path), rec("readTable", "awp_metadata"));
    assert_eq!(
        read_schema(&data("awp_metadata.parquet")),
        rec("readTable", "awp_metadata")
    );
    assert_eq!(
        ParquetFile::read_table(&path)
            .unwrap()
            .schema()
            .metadata()
            .get("x"),
        Some(&"y".to_owned())
    );
}

/// Records `concatenateAndWriteToParquet ...`, `readTable concat_*` and
/// `tableRowCount ...`.
#[test]
fn concatenate_and_write_to_parquet_matches_the_release_build() {
    let dir = scratch("concat");
    let config = ParquetWriteConfig::default();
    assert!(rec("concatenateAndWriteToParquet", "empty").starts_with("true"));
    let path = dir.join("concat_empty.parquet");
    aio::concatenate_and_write_to_parquet(&[], &path, &config, None).unwrap();
    assert!(!path.exists());
    let a = with_metadata(&small_table(2), &[("first", "1")]);
    let c = with_metadata(&small_table(3), &[("second", "2")]);
    let two = dir.join("concat_two.parquet");
    assert_eq!(rec("concatenateAndWriteToParquet", "two"), "true");
    aio::concatenate_and_write_to_parquet(&[a.clone(), c.clone()], &two, &config, None).unwrap();
    assert_eq!(read_schema(&two), rec("readTable", "concat_two"));
    assert_eq!(
        read_schema(&data("concat_two.parquet")),
        rec("readTable", "concat_two")
    );
    let merged = dir.join("concat_md.parquet");
    let mut metadata = aio::KeyValueMetadata::new();
    metadata.append("merged", "yes");
    assert_eq!(
        rec("concatenateAndWriteToParquet", "two with metadata"),
        "true"
    );
    aio::concatenate_and_write_to_parquet(&[a.clone(), c], &merged, &config, Some(&metadata))
        .unwrap();
    assert_eq!(read_schema(&merged), rec("readTable", "concat_md"));
    assert_eq!(
        read_schema(&data("concat_md.parquet")),
        rec("readTable", "concat_md")
    );
    // The rows are the two tables in order.
    let back = ParquetFile::read_table(&merged).unwrap();
    let ids: Vec<i64> = back.column(0).as_primitive::<Int64Type>().values().to_vec();
    assert_eq!(ids, vec![0, 10, 0, 10, 20]);
    let other = RecordBatch::new_empty(Arc::new(Schema::new(vec![Field::new(
        "id",
        DataType::Int32,
        true,
    )])));
    let path = dir.join("concat_mismatch.parquet");
    assert!(rec("concatenateAndWriteToParquet", "mismatch").starts_with("false"));
    assert!(aio::concatenate_and_write_to_parquet(&[a, other], &path, &config, None).is_err());
    assert!(!path.exists());
    let non_null = RecordBatch::try_new(
        Arc::new(Schema::new(vec![
            Field::new("id", DataType::Int64, false),
            Field::new("name", DataType::Utf8, true),
            Field::new("val", DataType::Float64, true),
        ])),
        vec![
            Arc::new(Int64Array::from(vec![1])),
            Arc::new(StringArray::from(vec!["x"])),
            Arc::new(Float64Array::from(vec![1.0])),
        ],
    )
    .unwrap();
    let path = dir.join("concat_nullability.parquet");
    assert!(rec("concatenateAndWriteToParquet", "nullability mismatch").starts_with("false"));
    assert!(
        aio::concatenate_and_write_to_parquet(&[small_table(1), non_null], &path, &config, None)
            .is_err()
    );
    assert!(!path.exists());
    assert_eq!(
        aio::table_row_count(None).to_string(),
        rec("tableRowCount", "null")
    );
    assert_eq!(
        aio::table_row_count(Some(&small_table(5))).to_string(),
        rec("tableRowCount", "five")
    );
}

/// The first record of the concatenation records under `kind` and `key`.
///
/// This driver prints Arrow's multi-line `Status` texts raw, so a line without
/// two tabs continues the previous record's value and is skipped here; the
/// single-line records this test reads are unaffected.
fn rec_concat(kind: &str, key: &str) -> String {
    RECORDS_CONCAT
        .lines()
        .filter_map(|line| {
            let mut parts = line.splitn(3, '\t');
            Some((parts.next()?, parts.next()?, parts.next()?))
        })
        .find(|(k, key2, _)| *k == kind && *key2 == key)
        .map(|(_, _, value)| value.to_owned())
        .unwrap_or_else(|| panic!("no record {kind}\t{key}"))
}

/// Records `list ...`, `struct child name` and
/// `concatenateAndWriteToParquet ...` of the concatenation oracle: tables
/// whose schemas differ in a list's element *name* are equal for Arrow C++'s
/// `Schema::Equals(check_metadata = false)` and concatenate, the result
/// keeping the first table's spelling; a different element nullability or a
/// different struct child name is not equal and does not.
#[test]
fn concatenation_follows_schema_equality_for_nested_names() {
    fn list_table(element: &str, nullable: bool, values: &[i32]) -> RecordBatch {
        let mut builder = ListBuilder::new(Int32Builder::new()).with_field(Arc::new(Field::new(
            element,
            DataType::Int32,
            nullable,
        )));
        for value in values {
            builder.values().append_value(*value);
            builder.append(true);
        }
        RecordBatch::try_from_iter([("v", Arc::new(builder.finish()) as ArrayRef)]).unwrap()
    }
    fn struct_table(child: &str, values: &[i32]) -> RecordBatch {
        let fields = Fields::from(vec![Field::new(child, DataType::Int32, true)]);
        let array = StructArray::new(
            fields,
            vec![Arc::new(Int32Array::from(values.to_vec())) as ArrayRef],
            None,
        );
        RecordBatch::try_from_iter([("s", Arc::new(array) as ArrayRef)]).unwrap()
    }
    let dir = scratch("concat_nested");
    let config = ParquetWriteConfig::default();
    // (record name, first table, second table), the driver's rows in order.
    let cases = [
        (
            "list item+element",
            list_table("item", true, &[1, 2]),
            list_table("element", true, &[3]),
        ),
        (
            "list element+item",
            list_table("element", true, &[1]),
            list_table("item", true, &[2, 3]),
        ),
        (
            "list item+item",
            list_table("item", true, &[1]),
            list_table("item", true, &[2]),
        ),
        (
            "list element nullability",
            list_table("item", true, &[1]),
            list_table("item", false, &[2]),
        ),
        (
            "list other name",
            list_table("item", true, &[1]),
            list_table("value", true, &[2]),
        ),
        (
            "struct child name",
            struct_table("a", &[1]),
            struct_table("b", &[2]),
        ),
    ];
    for (name, x, y) in cases {
        let equal = rec_concat(name, "schema_equals_no_metadata") == "true";
        assert_eq!(
            validation::schemas_equal(&x.schema(), &y.schema()),
            equal,
            "{name}: Schema::Equals"
        );
        let accepted = rec_concat(name, "ConcatenateTables_ok") == "true";
        assert_eq!(accepted, equal, "{name}: the driver's own consistency");
        let path = dir.join(format!("{}.parquet", name.replace([' ', '+'], "_")));
        let written =
            aio::concatenate_and_write_to_parquet(&[x.clone(), y.clone()], &path, &config, None);
        assert_eq!(written.is_ok(), accepted, "{name}: {written:?}");
        if accepted {
            let back = ParquetFile::read_table(&path).unwrap();
            assert_eq!(
                back.num_rows().to_string(),
                rec_concat(name, "combined_rows"),
                "{name}: rows"
            );
            assert_eq!(back.num_rows(), x.num_rows() + y.num_rows(), "{name}");
        } else {
            assert!(!path.exists(), "{name}: a refusal writes nothing");
        }
    }
    // The OpenMS writer over the item/element pair, as the driver recorded it:
    // three rows, and the Parquet file spells the element `element` because
    // `writeTableToParquet` stores no Arrow schema.
    assert_eq!(
        rec_concat("concatenateAndWriteToParquet", "item+element"),
        "true"
    );
    assert_eq!(rec_concat("concatenateAndWriteToParquet", "rowCount"), "3");
    assert_eq!(
        rec_concat("concatenateAndWriteToParquet", "read_schema"),
        "v: list<element: int32>"
    );
}

/// `[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}`,
/// the driver's `generateUuidV4 shape` pattern.
fn uuid_shape(text: &str) -> bool {
    let groups: Vec<&str> = text.split('-').collect();
    let hex = |s: &str| {
        s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    };
    groups.len() == 5
        && [8, 4, 4, 4, 12]
            .iter()
            .zip(&groups)
            .all(|(n, g)| g.len() == *n && hex(g))
        && groups[2].starts_with('4')
        && matches!(groups[3].as_bytes()[0], b'8' | b'9' | b'a' | b'b')
}

/// Records `generateUuidV4 shape` and `qpxFileMetadata ...`.
#[test]
fn qpx_file_metadata_matches_the_release_build() {
    assert_eq!(
        uuid_shape(&aio::generate_uuid_v4()).to_string(),
        rec("generateUuidV4", "shape")
    );
    assert_ne!(aio::generate_uuid_v4(), aio::generate_uuid_v4());
    let config = ParquetWriteConfig::default();
    let extra = BTreeMap::from([
        ("scan_format".to_owned(), "scan".to_owned()),
        ("aaa".to_owned(), "1".to_owned()),
    ]);
    let metadata = aio::qpx_file_metadata("psm_file", &config, &extra).unwrap();
    let records: Vec<(String, String)> = parse(RECORDS)
        .into_iter()
        .filter(|(k, key, _)| k == "qpxFileMetadata" && key.starts_with("psm_file/"))
        .map(|(_, key, value)| (key, value))
        .collect();
    assert_eq!(metadata.len(), records.len());
    for (i, (key, value)) in metadata.pairs().iter().enumerate() {
        let shown = match key.as_str() {
            "creation_date" => {
                let b = value.as_bytes();
                let shape = value.len() == 20
                    && b[4] == b'-'
                    && b[7] == b'-'
                    && b[10] == b'T'
                    && b[13] == b':'
                    && b[16] == b':'
                    && b[19] == b'Z';
                if shape {
                    "<iso8601Z>".to_owned()
                } else {
                    format!("BAD:{value}")
                }
            }
            "uuid" if uuid_shape(value) => "<uuid>".to_owned(),
            _ => value.clone(),
        };
        assert_eq!((format!("psm_file/{i}/{key}"), shown), records[i]);
    }
    for file_type in ["feature_file", "pg_file", "run_file"] {
        let metadata = aio::qpx_file_metadata(file_type, &config, &BTreeMap::new()).unwrap();
        let keys: Vec<&str> = metadata.pairs().iter().map(|(k, _)| k.as_str()).collect();
        assert_eq!(
            keys.join(","),
            rec("qpxFileMetadata", &format!("{file_type}/keys"))
        );
        for key in [
            "primary_key",
            "identity_composite",
            "compression_format",
            "creator",
            "software_provider",
            "qpx_version",
        ] {
            if let Some(value) = metadata.get(key) {
                assert_eq!(value, rec("qpxFileMetadata", &format!("{file_type}/{key}")));
            }
        }
    }
    for (name, compression) in [
        ("none", ParquetCompression::None),
        ("snappy", ParquetCompression::Snappy),
        ("gzip", ParquetCompression::Gzip),
        ("lz4", ParquetCompression::Lz4),
    ] {
        let config = ParquetWriteConfig {
            compression,
            ..ParquetWriteConfig::default()
        };
        let shown = match aio::qpx_file_metadata("psm_file", &config, &BTreeMap::new()) {
            Ok(metadata) => metadata.get("compression_format").unwrap().to_owned(),
            Err(_) => "nullptr".to_owned(),
        };
        assert_eq!(
            shown,
            rec("qpxFileMetadata", &format!("compression {name}"))
        );
    }
    let duplicate = aio::qpx_file_metadata(
        "psm_file",
        &config,
        &BTreeMap::from([("uuid".to_owned(), "mine".to_owned())]),
    )
    .unwrap();
    let uuids: String = duplicate
        .pairs()
        .iter()
        .filter(|(k, _)| k == "uuid")
        .map(|(_, v)| if v == "mine" { "mine;" } else { "<uuid>;" })
        .collect();
    assert_eq!(
        format!(
            "{}|{uuids}|FindKey={}",
            duplicate.len(),
            duplicate.find_key("uuid").unwrap()
        ),
        rec("qpxFileMetadata", "extra uuid duplicates")
    );
    // A schema map keeps the first value of a repeated key.
    assert_ne!(duplicate.to_schema_metadata()["uuid"], "mine");
}

// ---------------------------------------------------------------------------
// ArrowIOHelpers reading
// ---------------------------------------------------------------------------

/// The driver's metadata entries: name, value, value_type, value is null.
const META_ENTRIES: &[(&str, &str, &str, bool)] = &[
    ("i1", "12", "int", false),
    ("i2", " 12", "int", false),
    ("i3", "12abc", "int", false),
    ("i4", "abc", "int", false),
    ("i5", "9223372036854775807", "int", false),
    ("i6", "9223372036854775808", "int", false),
    ("i7", "-9223372036854775808", "int", false),
    ("i8", "+5", "int", false),
    ("i9", "0x10", "int", false),
    ("i10", "", "int", false),
    ("i11", "\x0b\x0c 7", "int", false),
    ("i12", "1.9", "int", false),
    ("i13", "- 1", "int", false),
    ("i14", "", "int", true),
    ("d1", "1.5", "double", false),
    ("d2", " 1.5", "double", false),
    ("d3", "1.5x", "double", false),
    ("d4", "inf", "double", false),
    ("d5", "-Infinity", "double", false),
    ("d6", "nan", "double", false),
    ("d7", "NaN", "double", false),
    ("d8", "nan(123)", "double", false),
    ("d9", "1e-310", "double", false),
    ("d10", "1e-400", "double", false),
    ("d11", "1e400", "double", false),
    ("d12", "0x1p3", "double", false),
    ("d13", "0x1p-1074", "double", false),
    ("d14", "0x1.8p-1074", "double", false),
    ("d15", "", "double", false),
    ("d16", "1,5", "double", false),
    ("d17", "-0", "double", false),
    ("d18", ".5", "double", false),
    ("d19", "5.", "double", false),
    ("d20", "1e", "double", false),
    ("d21", "2.2250738585072011e-308", "double", false),
    ("d22", "2.2250738585072014e-308", "double", false),
    ("d23", "0x1.fffffffffffffp1023", "double", false),
    ("d24", "0x1p1024", "double", false),
    ("d25", "4.9406564584124654e-324", "double", false),
    ("d26", "0e-400", "double", false),
    ("d27", "0x", "double", false),
    ("d28", "0x.8", "double", false),
    ("d29", "infinityx", "double", false),
    ("d30", "nanx", "double", false),
    ("d31", "0.1", "double", false),
    ("d32", "123456789012345678901234567890", "double", false),
    ("f1", "2.5", "float", false),
    ("il1", "[1, 2, 3]", "int_list", false),
    ("il2", "[]", "int_list", false),
    ("il3", "1,2", "int_list", false),
    ("il4", "[ 1 , 2 ]", "int_list", false),
    ("il5", "[1, x]", "int_list", false),
    ("il6", "[2147483648]", "int_list", false),
    ("il7", "[1.5]", "int_list", false),
    ("il8", "[1,]", "int_list", false),
    ("il9", "[", "int_list", false),
    ("il10", "", "int_list", false),
    ("il11", "[+3, -4]", "int_list", false),
    ("il12", "[ ]", "int_list", false),
    ("dl1", "[1.5, 2]", "double_list", false),
    ("dl2", "[]", "double_list", false),
    ("dl3", "[inf, nan, -inf]", "double_list", false),
    ("dl4", "[1e-310]", "double_list", false),
    ("dl5", "[1e400]", "double_list", false),
    ("dl6", "[0x1p3]", "double_list", false),
    ("dl7", "[1.5x]", "double_list", false),
    ("dl8", "[NaN(abc), Infinity]", "double_list", false),
    ("sl1", "[a, b ,c]", "string_list", false),
    ("sl2", "[]", "string_list", false),
    ("sl3", "[ a ]", "string_list", false),
    ("sl4", "a,b", "string_list", false),
    ("sl5", "[a, [b]]", "string_list", false),
    ("sl6", "[,]", "string_list", false),
    ("sl7", "", "string_list", false),
    ("sl8", "[\ta\n, \x0bb]", "string_list", false),
    ("s1", "x", "string", false),
    ("s2", "x", "bool", false),
    ("s3", "7", "", false),
    ("skip", "1", "int", false),
    ("dup", "1", "int", false),
    ("dup", "2.5", "double", false),
];

/// The driver's metadata column: row 0 holds every entry, row 1 is null and
/// row 2 is an empty list.
fn meta_column() -> ListArray {
    let names = StringArray::from(META_ENTRIES.iter().map(|e| Some(e.0)).collect::<Vec<_>>());
    let values = StringArray::from(
        META_ENTRIES
            .iter()
            .map(|e| if e.3 { None } else { Some(e.1) })
            .collect::<Vec<_>>(),
    );
    let types = StringArray::from(META_ENTRIES.iter().map(|e| Some(e.2)).collect::<Vec<_>>());
    let fields = Fields::from(vec![
        Field::new("name", DataType::Utf8, true),
        Field::new("value", DataType::Utf8, true),
        Field::new("value_type", DataType::Utf8, true),
    ]);
    let entries = StructArray::try_new(
        fields.clone(),
        vec![
            Arc::new(names) as ArrayRef,
            Arc::new(values),
            Arc::new(types),
        ],
        None,
    )
    .unwrap();
    let n = META_ENTRIES.len() as i32;
    ListArray::try_new(
        Arc::new(Field::new("item", DataType::Struct(fields), true)),
        arrow_buffer::OffsetBuffer::new(vec![0, n, n, n].into()),
        Arc::new(entries),
        Some(arrow_buffer::NullBuffer::from(vec![true, false, true])),
    )
    .unwrap()
}

/// The driver's `describe(DataValue)`, minus the `toString` of a floating
/// value (the port's float formatting is not under test here; the value is
/// compared exactly through its `%a` spelling).
fn describe(value: &MetaValue) -> String {
    match value.data() {
        MetaValueData::String(text) => format!("STRING_VALUE|{text}"),
        MetaValueData::Integer(number) => format!("INT_VALUE|{number}"),
        MetaValueData::Float(number) => format!("DOUBLE_VALUE|{}", hexd(*number)),
        MetaValueData::StringList(values) => format!("STRING_LIST|[{}]", values.join(", ")),
        MetaValueData::IntegerList(values) => format!(
            "INT_LIST|[{}]",
            values
                .iter()
                .map(i64::to_string)
                .collect::<Vec<_>>()
                .join(", ")
        ),
        MetaValueData::FloatList(values) => {
            let mut text = "DOUBLE_LIST".to_owned();
            for value in values {
                text.push('|');
                text.push_str(&hexd(*value));
            }
            text
        }
        MetaValueData::Empty => "EMPTY_VALUE|".to_owned(),
    }
}

/// The record without the float `toString` part `describe` leaves out.
fn comparable(record: &str) -> String {
    let parts: Vec<&str> = record.split('|').collect();
    match parts[0] {
        "DOUBLE_VALUE" => format!("DOUBLE_VALUE|{}", parts[2]),
        "DOUBLE_LIST" => {
            let mut text = "DOUBLE_LIST".to_owned();
            for part in &parts[2..] {
                text.push('|');
                text.push_str(part);
            }
            text
        }
        _ => record.to_owned(),
    }
}

/// Records `readMetaValues ...`: the typing of every entry, the exclusion,
/// a later entry replacing an earlier one, and null, empty and missing rows.
#[test]
fn read_meta_values_matches_the_release_build() {
    let column = meta_column();
    let mut target = MetaInfo::new();
    aio::read_meta_values(Some(&column), 0, &mut target, &["skip"]).unwrap();
    let records: Vec<(String, String)> = parse(RECORDS)
        .into_iter()
        .filter(|(k, key, _)| {
            k == "readMetaValues"
                && !["null row", "empty row", "null array"].contains(&key.as_str())
        })
        .map(|(_, key, value)| (key, value))
        .collect();
    assert_eq!(target.len(), records.len());
    assert!(!target.contains_key("skip"));
    for (key, record) in &records {
        let value = target.get(key).unwrap_or_else(|| panic!("{key} missing"));
        assert_eq!(describe(value), comparable(record), "{key}");
    }
    let mut null_row = MetaInfo::new();
    aio::read_meta_values(Some(&column), 1, &mut null_row, &[]).unwrap();
    assert_eq!(
        null_row.is_empty().to_string(),
        rec("readMetaValues", "null row")
    );
    let mut empty_row = MetaInfo::new();
    aio::read_meta_values(Some(&column), 2, &mut empty_row, &[]).unwrap();
    assert_eq!(
        empty_row.is_empty().to_string(),
        rec("readMetaValues", "empty row")
    );
    let mut none = MetaInfo::new();
    aio::read_meta_values(None, 0, &mut none, &[]).unwrap();
    assert_eq!(
        none.is_empty().to_string(),
        rec("readMetaValues", "null array")
    );
    // Native: a row outside the column adds nothing, and a column of another
    // type is refused (the source's unchecked cast).
    aio::read_meta_values(Some(&column), 3, &mut none, &[]).unwrap();
    assert!(none.is_empty());
    assert!(aio::read_meta_values(Some(&Int64Array::from(vec![1])), 0, &mut none, &[]).is_err());
    assert!(none.is_empty());
}

/// Records `getDoubleValue ...`, `getStringValue null array`, `isNull null
/// array` and `getColumn ...` of `ArrowIOHelpers`.
#[test]
fn arrow_io_helper_getters_match_the_release_build() {
    let doubles = Float64Array::from(vec![2.5]);
    assert_eq!(
        hexd(aio::get_double_value(Some(&doubles), 0, -1.0).unwrap()),
        rec("getDoubleValue", "value")
    );
    assert_eq!(
        hexd(aio::get_double_value(None, 0, -1.0).unwrap()),
        rec("getDoubleValue", "null array")
    );
    assert_eq!(
        aio::get_string_value(None, 0).unwrap(),
        rec("getStringValue", "null array")
    );
    assert_eq!(
        aio::is_null(None, 0).to_string(),
        rec("isNull", "null array")
    );
    assert_eq!(
        aio::get_column(&small_table(2), "nope", false)
            .is_none()
            .to_string(),
        rec("getColumn", "missing not required")
    );
    let both = arrow_select::concat::concat_batches(
        &small_table(2).schema(),
        &[small_table(2), small_table(3)],
    )
    .unwrap();
    // The second `getColumn two chunks length` record is ArrowIOHelpers's,
    // which combines chunks.
    assert_eq!(
        aio::get_column(&both, "id", true)
            .unwrap()
            .len()
            .to_string(),
        all(RECORDS, "getColumn", "two chunks length")[1]
    );
    // Typed getters: nulls and rows outside the array read as the default,
    // another type is refused.
    let ints = Int32Array::from(vec![Some(3), None]);
    assert_eq!(aio::get_int32_value(Some(&ints), 0, -1).unwrap(), 3);
    assert_eq!(aio::get_int32_value(Some(&ints), 1, -1).unwrap(), -1);
    assert_eq!(aio::get_int32_value(Some(&ints), 9, -1).unwrap(), -1);
    assert!(aio::get_int64_value(Some(&ints), 0, -1).is_err());
    assert_eq!(
        aio::get_float_value(Some(&Float32Array::from(vec![0.5f32])), 0, 0.0).unwrap(),
        0.5
    );
    assert!(aio::get_bool_value(Some(&BooleanArray::from(vec![true])), 0, false).unwrap());
    assert_eq!(
        aio::get_string_value(Some(&StringArray::from(vec!["s"])), 5).unwrap(),
        ""
    );
    assert!(aio::get_string_value(Some(&LargeStringArray::from(vec!["s"])), 0).is_err());
    assert!(
        aio::is_null(Some(&ints), 1)
            && aio::is_null(Some(&ints), 2)
            && !aio::is_null(Some(&ints), 0)
    );
}

// ---------------------------------------------------------------------------
// QPX naming and labels
// ---------------------------------------------------------------------------

/// Records `qpxRunFileName ...`, `qpxScanFormat ...`, `qpxScanComponents
/// ...`, `qpxScanFormat[] ...` and `qpxWarnOnRunNameCollisions ...`.
#[test]
fn qpx_naming_helpers_match_the_release_build() {
    for (kind, key, value) in parse(RECORDS) {
        match kind.as_str() {
            "qpxRunFileName" => assert_eq!(aio::qpx_run_file_name(&key), value, "{key:?}"),
            "qpxScanFormat" => assert_eq!(aio::qpx_scan_format(&key), value, "{key:?}"),
            "qpxScanComponents" => {
                let components = aio::qpx_scan_components(&key);
                let mut text = components.len().to_string();
                for c in components {
                    text.push_str(&format!("|{c}"));
                }
                assert_eq!(text, value, "{key:?}");
            }
            _ => {}
        }
    }
    assert_eq!(
        aio::qpx_scan_format_of_all(&["scan=1", "controllerType=0 scan=2", "foo"]).0,
        rec("qpxScanFormat[]", "scan+scan")
    );
    let (mixed, warning) = aio::qpx_scan_format_of_all(&["index=1", "scan=2"]);
    assert_eq!(mixed, rec("qpxScanFormat[]", "index+scan"));
    assert!(warning.is_some());
    assert_eq!(
        aio::qpx_scan_format_of_all(&["foo", ""]).0,
        rec("qpxScanFormat[]", "none")
    );
    assert_eq!(
        aio::qpx_scan_format_of_all::<&str>(&[]).0,
        rec("qpxScanFormat[]", "empty")
    );
    for (key, paths) in [
        ("collide", vec!["/a/run1.mzML", "/b/run1.mzML"]),
        ("same path twice", vec!["/a/x.mzML", "/a/x.mzML"]),
        ("empties", vec!["", ""]),
        ("distinct", vec!["/a/x.mzML", "/a/y.mzML"]),
    ] {
        assert_eq!(
            aio::qpx_warn_on_run_name_collisions("ctx", &paths).to_string(),
            rec("qpxWarnOnRunNameCollisions", key),
            "{key}"
        );
    }
    let lines = aio::qpx_run_name_collisions("ctx", &["/a/run1.mzML", "/b/run1.mzML"]);
    assert_eq!(lines.len(), 1);
    assert!(lines[0].contains("'run1'") && lines[0].contains("/a/run1.mzML, /b/run1.mzML"));
}

/// Records `qpxPsmRunFileName ...` of part 2, where the metadata names were
/// registered before the keys were resolved (part 1's records show the
/// source's fallback when they were not; the port has no registry).
#[test]
fn qpx_psm_run_file_name_matches_the_release_build() {
    let keys = aio::QpxRunFileNameKeys::new();
    let meta = |pairs: Vec<(&str, MetaValue)>| -> MetaInfo {
        pairs.into_iter().map(|(k, v)| (k.to_owned(), v)).collect()
    };
    let mut hit = MetaInfo::new();
    let mut id = MetaInfo::new();
    let run = |hit: &MetaInfo, id: &MetaInfo, fallback: &str| {
        aio::qpx_psm_run_file_name(hit, id, fallback, &keys)
    };
    assert_eq!(
        run(&hit, &id, "/r/fallback.mzML"),
        rec2("qpxPsmRunFileName", "fallback")
    );
    assert_eq!(
        run(&hit, &id, ""),
        rec2("qpxPsmRunFileName", "empty fallback")
    );
    id.insert("reference_file_name".into(), "/i/ident.mzML".into());
    assert_eq!(
        run(&hit, &id, "/r/fallback.mzML"),
        rec2("qpxPsmRunFileName", "identification")
    );
    id.insert("run_file_name".into(), "/i/idrun.raw".into());
    assert_eq!(
        run(&hit, &id, "/r/fallback.mzML"),
        rec2("qpxPsmRunFileName", "identification run_file_name ignored")
    );
    hit.insert("run_file_name".into(), "/h/run.raw".into());
    assert_eq!(
        run(&hit, &id, "/r/fallback.mzML"),
        rec2("qpxPsmRunFileName", "hit run_file_name")
    );
    hit.insert("reference_file_name".into(), "/h/ref.mzML".into());
    assert_eq!(
        run(&hit, &id, "/r/fallback.mzML"),
        rec2("qpxPsmRunFileName", "hit reference_file_name")
    );
    let h2 = meta(vec![
        ("reference_file_name", "".into()),
        ("run_file_name", "/h2/second.raw".into()),
    ]);
    assert_eq!(
        run(&h2, &id, "/r/fallback.mzML"),
        rec2("qpxPsmRunFileName", "empty hit reference")
    );
    let h3 = meta(vec![("reference_file_name", 42i64.into())]);
    assert_eq!(
        run(&h3, &id, "/r/fallback.mzML"),
        rec2("qpxPsmRunFileName", "numeric hit reference")
    );
    let h4 = meta(vec![(
        "reference_file_name",
        MetaValue::try_from(1.5).unwrap(),
    )]);
    assert_eq!(
        run(&h4, &id, "/r/fallback.mzML"),
        rec2("qpxPsmRunFileName", "double hit reference")
    );
    let h5 = meta(vec![(
        "reference_file_name",
        MetaValue::try_from(0.1).unwrap(),
    )]);
    assert_eq!(
        run(&h5, &id, "/r/fallback.mzML"),
        rec2("qpxPsmRunFileName", "double 0.1 hit reference")
    );
    let h6 = meta(vec![(
        "reference_file_name",
        vec!["/l/a.mzML".to_owned(), "b".to_owned()].into(),
    )]);
    assert_eq!(
        run(&h6, &id, "/r/fallback.mzML"),
        rec2("qpxPsmRunFileName", "string list hit reference")
    );
    assert_eq!(
        run(&MetaInfo::new(), &MetaInfo::new(), ""),
        rec2("qpxPsmRunFileName", "nothing at all")
    );
}

/// Records `qpxIntensityLabel ...`, `isobaric ...` and
/// `qpxIsCanonicalIntensityLabel ...`.
#[test]
fn qpx_intensity_label_matches_the_release_build() {
    let mut labels = 0;
    for (kind, key, value) in parse(RECORDS) {
        match kind.as_str() {
            "qpxIntensityLabel" => {
                let (label, channel) = key.split_once('|').unwrap();
                assert_eq!(aio::qpx_intensity_label(label, channel), value, "{key:?}");
                labels += 1;
            }
            "qpxIsCanonicalIntensityLabel" => assert_eq!(
                aio::qpx_is_canonical_intensity_label(&key).to_string(),
                value,
                "{key:?}"
            ),
            "isobaric" => {
                let prefix = if key.starts_with("tmt") {
                    "TMT"
                } else {
                    "ITRAQ"
                };
                for channel in value.split(',') {
                    assert!(
                        aio::qpx_is_canonical_intensity_label(&format!("{prefix}{channel}")),
                        "{key} {channel}"
                    );
                    // Every method name resolves, with and without a channel name.
                    assert_eq!(
                        aio::qpx_intensity_label(&format!("{key}_{channel}"), ""),
                        format!("{prefix}{channel}")
                    );
                }
            }
            _ => {}
        }
    }
    assert_eq!(labels, 34);
}

/// Records `qpxIntensityLabels ...`: labels of whole maps, grouped by source
/// file.
#[test]
fn qpx_intensity_labels_match_the_release_build() {
    let labels_of = |headers: &[(&str, &str, &str)]| {
        let mut map = ConsensusMap::default();
        for (i, (file, label, channel)) in headers.iter().enumerate() {
            let mut header = ColumnHeader {
                filename: (*file).to_owned(),
                label: (*label).to_owned(),
                ..ColumnHeader::default()
            };
            if !channel.is_empty() {
                header
                    .metadata
                    .insert("channel_name".into(), (*channel).into());
            }
            map.column_headers.insert(i as u64, header);
        }
        aio::qpx_intensity_labels(&map)
            .into_iter()
            .map(|(k, v)| format!("{k}={v};"))
            .collect::<String>()
    };
    type Header<'a> = (&'a str, &'a str, &'a str);
    let cases: &[(&str, &[Header])] = &[
        (
            "two-plex light Arg6",
            &[("r", "light", ""), ("r", "Arg6", "")],
        ),
        (
            "three-plex",
            &[("r", "", ""), ("r", "Arg6Lys4", ""), ("r", "Arg10Lys8", "")],
        ),
        ("two light", &[("r", "light", ""), ("r", "light", "")]),
        (
            "tmt",
            &[
                ("r", "tmt10plex_126", "126"),
                ("r", "tmt10plex_127N", "127N"),
            ],
        ),
        (
            "label-free two files",
            &[("a", "label-free", ""), ("b", "label-free", "")],
        ),
        (
            "silac unknown role",
            &[("r", "light", ""), ("r", "SILAC foo", "")],
        ),
        (
            "four-plex silac",
            &[
                ("r", "light", ""),
                ("r", "Arg6", ""),
                ("r", "Arg10", ""),
                ("r", "heavy", ""),
            ],
        ),
        (
            "mixed per file",
            &[
                ("a", "light", ""),
                ("a", "Lys8", ""),
                ("b", "tmt6plex_126", "126"),
            ],
        ),
        ("two heavy", &[("r", "Arg10", ""), ("r", "heavy", "")]),
        (
            "argx lysx plain",
            &[("r", "arginine", ""), ("r", "lysine", "")],
        ),
        (
            "empty channel meta",
            &[("r", "light", ""), ("r", "Arg6", "")],
        ),
    ];
    for (key, headers) in cases {
        assert_eq!(labels_of(headers), rec("qpxIntensityLabels", key), "{key}");
    }
}

// ---------------------------------------------------------------------------
// Nested columns (part 2)
// ---------------------------------------------------------------------------

/// Records `writeTableToParquet nested`, `parquet_column ...`,
/// `readTable_field ...` and `readTable_value ...` of part 2: a list column
/// is stored in the Parquet-compliant layout, and `readTable` returns a
/// `timestamp[s]` column as `timestamp[ms]` even though the Arrow schema is
/// stored, but restores `large_string`.
#[test]
fn nested_columns_round_trip_as_in_the_release_build() {
    let dir = scratch("nested");
    let mut lists = ListBuilder::new(Int32Builder::new());
    lists.values().append_value(1);
    lists.values().append_value(2);
    lists.append(true);
    lists.append(false);
    let lists = lists.finish();
    let schema = Arc::new(Schema::new(vec![
        Field::new("l", lists.data_type().clone(), true),
        Field::new("t", DataType::Timestamp(TimeUnit::Second, None), true),
        Field::new("ls", DataType::LargeUtf8, true),
    ]));
    let table = RecordBatch::try_new(
        schema,
        vec![
            Arc::new(lists),
            Arc::new(TimestampSecondArray::from(vec![Some(1_700_000_000), None])),
            Arc::new(LargeStringArray::from(vec![Some("big"), None])),
        ],
    )
    .unwrap();
    let path = dir.join("awp_nested.parquet");
    assert_eq!(rec2("writeTableToParquet", "nested"), "1");
    aio::write_table_to_parquet(&table, &path, &ParquetWriteConfig::default()).unwrap();
    for file in [path.clone(), data("awp_nested.parquet")] {
        let reader = SerializedFileReader::new(std::fs::File::open(&file).unwrap()).unwrap();
        let descr = reader.metadata().file_metadata().schema_descr();
        for i in 0..descr.num_columns() {
            assert_eq!(
                descr.column(i).path().string(),
                rec2("parquet_column", &i.to_string()),
                "{}",
                file.display()
            );
        }
        let back = ParquetFile::read_table(&file).unwrap();
        for (i, field) in back.schema().fields().iter().enumerate() {
            assert_eq!(
                cpp_field_string(field),
                rec2("readTable_field", &i.to_string()),
                "{}",
                file.display()
            );
        }
        let times = back.column(1).as_primitive::<TimestampMillisecondType>();
        assert_eq!(times.value(0).to_string(), rec2("readTable_value", "t0"));
        assert_eq!(
            u8::from(times.is_null(1)).to_string(),
            rec2("readTable_value", "t1 null")
        );
    }
}

/// Native: writing the same table twice, or two equal tables built
/// separately, gives identical bytes (the crate-decision entry's
/// determinism promise). Schema metadata lives in randomly seeded hash maps,
/// so this also checks that the writers do not encode it in hash order.
#[test]
fn parquet_output_is_byte_identical_across_writes() {
    let dir = scratch("determinism");
    let pairs = [("b", "2"), ("a", "1"), ("c", "3"), ("qpx_note", "x")];
    for (i, table) in [
        with_metadata(&small_table(50), &pairs),
        with_metadata(&small_table(50), &pairs),
    ]
    .iter()
    .enumerate()
    {
        for round in 0..2 {
            let path = dir.join(format!("awp_{i}_{round}.parquet"));
            aio::write_table_to_parquet(table, &path, &ParquetWriteConfig::default()).unwrap();
            let path = dir.join(format!("pf_{i}_{round}.parquet"));
            ParquetFile::write_table(table, &path).unwrap();
        }
    }
    for prefix in ["awp", "pf"] {
        let first = std::fs::read(dir.join(format!("{prefix}_0_0.parquet"))).unwrap();
        for (i, round) in [(0, 1), (1, 0), (1, 1)] {
            let other = std::fs::read(dir.join(format!("{prefix}_{i}_{round}.parquet"))).unwrap();
            assert_eq!(first, other, "{prefix} {i} {round}");
        }
    }
    // The metadata survives both writers.
    for prefix in ["awp", "pf"] {
        let back = ParquetFile::read_table(dir.join(format!("{prefix}_0_0.parquet"))).unwrap();
        for (key, value) in pairs {
            assert_eq!(
                back.schema().metadata().get(key).map(String::as_str),
                Some(value)
            );
        }
    }
}
