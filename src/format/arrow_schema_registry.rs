// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Arrow schema definitions of the OpenMS Parquet tables, and their validation.
//!
//! Ports `FORMAT/ArrowSchemaRegistry.h` and `ArrowSchemaRegistry.cpp`: the
//! `ArrowSchemaValidation` namespace (as [`validation`](crate::format::arrow_schema_registry::validation)) and the 22 schema
//! structs, each a set of column-name constants, the Arrow types of its nested
//! columns and the complete [`schema`](crate::format::arrow_schema_registry::ProteinSchema::schema). See
//! `docs/ARROW_PARQUET_SUPPORT.md` for the full API mapping.
//!
//! A source `std::shared_ptr<arrow::Schema>` is an
//! [`arrow_schema::SchemaRef`] and a `std::shared_ptr<arrow::DataType>` an owned
//! [`arrow_schema::DataType`]. A field built with `arrow::field(name, type)` is
//! nullable, as in Arrow C++, and a list built with `arrow::list(type)` has the
//! nullable element field `item`, which is what arrow-rs's
//! [`DataType::new_list`](arrow_schema::DataType::new_list) builds too. Every schema here therefore compares
//! equal, field by field, to the one the executed C++ library returns; the
//! rustdoc of each `schema()` lists its columns exactly as C++
//! `Schema::ToString` prints them, and `tests/arrow_parquet.rs` checks that
//! spelling against the executed Release build.
//!
//! The column constants are associated `&'static str` constants with the
//! source's names, so `ProteinSchema::ACCESSION` reads as it does in C++.
//! The schema bodies below were transcribed mechanically from the pinned
//! `ArrowSchemaRegistry.cpp` (core `bc9cc12`) by a script that refuses any
//! construct it does not recognise, so no column can be dropped silently.

use arrow_schema::{DataType, Field, Fields, Schema, SchemaRef, TimeUnit};
use std::sync::Arc;

/// `arrow::field(name, type, nullable)`.
fn field(name: &str, data_type: DataType, nullable: bool) -> Field {
    Field::new(name, data_type, nullable)
}

/// `arrow::list(type)`: a list whose element field is the nullable `item`.
fn list(element: DataType) -> DataType {
    DataType::new_list(element, true)
}

/// `arrow::struct_({...})`.
fn struct_(fields: Vec<Field>) -> DataType {
    DataType::Struct(Fields::from(fields))
}

/// `arrow::schema({...})`, without schema-level metadata.
fn schema(fields: Vec<Field>) -> SchemaRef {
    Arc::new(Schema::new(fields))
}

/// Utilities for validating Arrow table schemas against expected definitions.
///
/// Ports the `ArrowSchemaValidation` namespace.
pub mod validation {
    use arrow_array::RecordBatch;
    use arrow_schema::{DataType, Field, Schema};
    use std::fmt;

    /// Validation strictness.
    ///
    /// Ports `ArrowSchemaValidation::Mode`. The header summarises it as
    /// "Strict requires exact match, Subset allows missing and extra
    /// columns"; the implementation is more precise, and this is what it does.
    #[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
    pub enum Mode {
        /// Compare position by position: the same number of fields, the same
        /// name at each index, a compatible type and the same nullability.
        /// The source's default argument.
        #[default]
        Strict,
        /// Look every actual field up by name in the expected schema. Fields
        /// the expected schema does not know are ignored (forward
        /// compatibility), expected fields the table lacks are not reported,
        /// and a known field must have a compatible type and the same
        /// nullability. Order is irrelevant.
        Subset,
    }

    /// Result of schema validation: a validity flag and the error messages.
    ///
    /// Ports `ArrowSchemaValidation::ValidationResult`. A default result is
    /// valid with no errors, as the source's member initialisers make it.
    #[derive(Clone, Debug, PartialEq, Eq, Hash)]
    pub struct ValidationResult {
        /// Whether the schema passed; the source `valid`, default `true`.
        pub valid: bool,
        /// One message per problem, in the order they were found; the source
        /// `errors`.
        pub errors: Vec<String>,
    }

    impl Default for ValidationResult {
        fn default() -> Self {
            Self {
                valid: true,
                errors: Vec::new(),
            }
        }
    }

    /// The source `toString()`: the errors joined by `"; "`, empty when there
    /// are none.
    impl fmt::Display for ValidationResult {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            f.write_str(&self.errors.join("; "))
        }
    }

    /// Validate a table's schema against an expected schema.
    ///
    /// Ports `ArrowSchemaValidation::validate`, whose `mode` defaults to
    /// [`Mode::Strict`]. Only the table's schema is consulted, so this is
    /// [`validate_schema`] on `table.schema()`; schema-level metadata is
    /// ignored in both modes.
    ///
    /// Two types are compatible when they are equal as Arrow C++
    /// `DataType::Equals` compares them by default, or when both are
    /// timestamps, of any unit and time zone: Parquet stores no seconds, so a
    /// `timestamp[s]` column comes back as `timestamp[ms]`. See
    /// [`types_equal`] for what "equal" includes.
    ///
    /// The messages are the source's, with each type spelled as Arrow C++
    /// `DataType::ToString` spells it ([`cpp_type_string`]):
    ///
    /// ```text
    /// Field count mismatch: got 1, expected 2
    /// Field name mismatch at index 0: got 'x', expected 'a'
    /// Type mismatch for field 'b': got int32, expected double
    /// Nullability mismatch for field 'a': got nullable, expected non-null
    /// Unexpected field 'c' not in expected schema
    /// Missing field 'b'
    /// ```
    ///
    /// In strict mode a name mismatch skips that index's type and
    /// nullability checks, as the source's `continue` does. The source takes
    /// a `std::shared_ptr<arrow::Table>` and dereferences it unchecked; a
    /// [`RecordBatch`] reference cannot be null.
    pub fn validate(table: &RecordBatch, expected: &Schema, mode: Mode) -> ValidationResult {
        validate_schema(table.schema().as_ref(), expected, mode)
    }

    /// [`validate`] for a schema without a table.
    ///
    /// Native convenience: the source validates only `table->schema()`, so a
    /// reader that has a schema but no rows yet can ask the same question.
    pub fn validate_schema(actual: &Schema, expected: &Schema, mode: Mode) -> ValidationResult {
        let mut result = ValidationResult::default();
        let mut fail = |message: String| {
            result.valid = false;
            result.errors.push(message);
        };
        let actual_n = actual.fields().len();
        let expected_n = expected.fields().len();
        match mode {
            Mode::Strict => {
                if actual_n != expected_n {
                    fail(format!(
                        "Field count mismatch: got {actual_n}, expected {expected_n}"
                    ));
                }
                for (position, (a, e)) in actual
                    .fields()
                    .iter()
                    .zip(expected.fields().iter())
                    .enumerate()
                {
                    if a.name() != e.name() {
                        fail(format!(
                            "Field name mismatch at index {position}: got '{}', expected '{}'",
                            a.name(),
                            e.name()
                        ));
                        continue;
                    }
                    check_field(a, e, &mut fail);
                }
                for a in actual.fields().iter().skip(expected_n) {
                    fail(format!(
                        "Unexpected field '{}' not in expected schema",
                        a.name()
                    ));
                }
                for e in expected.fields().iter().skip(actual_n) {
                    fail(format!("Missing field '{}'", e.name()));
                }
            }
            Mode::Subset => {
                for a in actual.fields() {
                    // `Schema::GetFieldByName` answers null for an absent name
                    // and also for a name the schema carries more than once;
                    // either way the field is skipped.
                    let mut matches = expected.fields().iter().filter(|e| e.name() == a.name());
                    let (Some(e), None) = (matches.next(), matches.next()) else {
                        continue;
                    };
                    check_field(a, e, &mut fail);
                }
            }
        }
        result
    }

    fn check_field(actual: &Field, expected: &Field, fail: &mut impl FnMut(String)) {
        if !types_compatible(actual.data_type(), expected.data_type()) {
            fail(format!(
                "Type mismatch for field '{}': got {}, expected {}",
                actual.name(),
                cpp_type_string(actual.data_type()),
                cpp_type_string(expected.data_type())
            ));
        }
        if actual.is_nullable() != expected.is_nullable() {
            let spell = |nullable: bool| if nullable { "nullable" } else { "non-null" };
            fail(format!(
                "Nullability mismatch for field '{}': got {}, expected {}",
                actual.name(),
                spell(actual.is_nullable()),
                spell(expected.is_nullable())
            ));
        }
    }

    /// The source's private `areTypesCompatible`: equal, or both timestamps.
    fn types_compatible(actual: &DataType, expected: &DataType) -> bool {
        types_equal(actual, expected)
            || matches!(
                (actual, expected),
                (DataType::Timestamp(..), DataType::Timestamp(..))
            )
    }

    /// Type equality as Arrow C++ `DataType::Equals(other)` decides it with
    /// its default `check_metadata = false`.
    ///
    /// Non-nested types compare as arrow-rs compares them (a timestamp's unit
    /// and time zone included). A struct's children must agree in number,
    /// name, nullability and type. A list-like type's element field must
    /// agree in nullability and type but **not in name**: Arrow C++ compares
    /// the element name only when metadata is checked, so `list<item: int32>`
    /// equals `list<element: int32>`, the spelling another Parquet writer
    /// gives the same column. Field metadata is never compared. arrow-rs's own
    /// `==` compares both, which is why this is not simply `a == b`.
    pub fn types_equal(a: &DataType, b: &DataType) -> bool {
        let element = |x: &Field, y: &Field| {
            x.is_nullable() == y.is_nullable() && types_equal(x.data_type(), y.data_type())
        };
        match (a, b) {
            (DataType::List(x), DataType::List(y))
            | (DataType::LargeList(x), DataType::LargeList(y))
            | (DataType::ListView(x), DataType::ListView(y))
            | (DataType::LargeListView(x), DataType::LargeListView(y)) => element(x, y),
            (DataType::FixedSizeList(x, n), DataType::FixedSizeList(y, m)) => {
                n == m && element(x, y)
            }
            (DataType::Map(x, sx), DataType::Map(y, sy)) => sx == sy && element(x, y),
            (DataType::Struct(x), DataType::Struct(y)) => {
                x.len() == y.len()
                    && x.iter().zip(y.iter()).all(|(p, q)| {
                        p.name() == q.name()
                            && p.is_nullable() == q.is_nullable()
                            && types_equal(p.data_type(), q.data_type())
                    })
            }
            (DataType::Dictionary(kx, vx), DataType::Dictionary(ky, vy)) => {
                types_equal(kx, ky) && types_equal(vx, vy)
            }
            _ => a == b,
        }
    }

    /// An Arrow type spelled as Arrow C++ `DataType::ToString()` spells it.
    ///
    /// Native helper for the validation messages, which the source builds
    /// from `ToString`: `utf8` is `string`, `float64` is `double`, a list is
    /// `list<item: int32>`, a struct `struct<a: int32, b: string not null>`,
    /// a timestamp `timestamp[s]` or `timestamp[ms, tz=UTC]`. The common
    /// types, every one the registry uses and their nested forms are spelled
    /// exactly (checked against the executed Release build); a type Arrow
    /// C++ spells in a way this does not model falls back to arrow-rs's
    /// `Display`.
    pub fn cpp_type_string(data_type: &DataType) -> String {
        use arrow_schema::TimeUnit;
        let unit = |u: &TimeUnit| match u {
            TimeUnit::Second => "s",
            TimeUnit::Millisecond => "ms",
            TimeUnit::Microsecond => "us",
            TimeUnit::Nanosecond => "ns",
        };
        match data_type {
            DataType::Null => "null".to_owned(),
            DataType::Boolean => "bool".to_owned(),
            DataType::Int8 => "int8".to_owned(),
            DataType::Int16 => "int16".to_owned(),
            DataType::Int32 => "int32".to_owned(),
            DataType::Int64 => "int64".to_owned(),
            DataType::UInt8 => "uint8".to_owned(),
            DataType::UInt16 => "uint16".to_owned(),
            DataType::UInt32 => "uint32".to_owned(),
            DataType::UInt64 => "uint64".to_owned(),
            DataType::Float16 => "halffloat".to_owned(),
            DataType::Float32 => "float".to_owned(),
            DataType::Float64 => "double".to_owned(),
            DataType::Utf8 => "string".to_owned(),
            DataType::LargeUtf8 => "large_string".to_owned(),
            DataType::Utf8View => "string_view".to_owned(),
            DataType::Binary => "binary".to_owned(),
            DataType::LargeBinary => "large_binary".to_owned(),
            DataType::BinaryView => "binary_view".to_owned(),
            DataType::FixedSizeBinary(n) => format!("fixed_size_binary[{n}]"),
            DataType::Date32 => "date32[day]".to_owned(),
            DataType::Date64 => "date64[ms]".to_owned(),
            DataType::Time32(u) => format!("time32[{}]", unit(u)),
            DataType::Time64(u) => format!("time64[{}]", unit(u)),
            DataType::Duration(u) => format!("duration[{}]", unit(u)),
            DataType::Timestamp(u, None) => format!("timestamp[{}]", unit(u)),
            DataType::Timestamp(u, Some(tz)) => format!("timestamp[{}, tz={tz}]", unit(u)),
            DataType::Decimal128(p, s) => format!("decimal128({p}, {s})"),
            DataType::Decimal256(p, s) => format!("decimal256({p}, {s})"),
            DataType::List(f) => format!("list<{}>", cpp_field_string(f)),
            DataType::LargeList(f) => format!("large_list<{}>", cpp_field_string(f)),
            DataType::ListView(f) => format!("list_view<{}>", cpp_field_string(f)),
            DataType::LargeListView(f) => format!("large_list_view<{}>", cpp_field_string(f)),
            DataType::FixedSizeList(f, n) => {
                format!("fixed_size_list<{}>[{n}]", cpp_field_string(f))
            }
            DataType::Struct(fields) => {
                let children: Vec<String> = fields.iter().map(|f| cpp_field_string(f)).collect();
                format!("struct<{}>", children.join(", "))
            }
            other => other.to_string(),
        }
    }

    /// A field spelled as Arrow C++ `Field::ToString()` spells it:
    /// `name: type`, followed by ` not null` for a non-nullable field.
    pub fn cpp_field_string(field: &Field) -> String {
        let mut text = format!("{}: {}", field.name(), cpp_type_string(field.data_type()));
        if !field.is_nullable() {
            text.push_str(" not null");
        }
        text
    }
}

/// Schema for the protein identification results table (`proteins.parquet`).
///
/// Ports `ProteinSchema`.
#[derive(Clone, Copy, Debug)]
pub struct ProteinSchema;

/// Schema for the protein group (indistinguishable group) results table
/// (`protein_groups.parquet`).
///
/// Ports `ProteinGroupSchema`.
#[derive(Clone, Copy, Debug)]
pub struct ProteinGroupSchema;

/// Schema for the search-engine parameters and settings table
/// (`search_params.parquet`).
///
/// Ports `SearchParamsSchema`.
#[derive(Clone, Copy, Debug)]
pub struct SearchParamsSchema;

/// Schema for the LC-MS feature table (`FeatureMap` features,
/// `features.parquet`).
///
/// Ports `FeatureSchema`.
#[derive(Clone, Copy, Debug)]
pub struct FeatureSchema;

/// Schema for the consensus feature table (`ConsensusMap` features,
/// `consensus_features.parquet`).
///
/// Ports `ConsensusFeatureSchema`. Every field is nullable: the source builds
/// them all with the two-argument `arrow::field`.
#[derive(Clone, Copy, Debug)]
pub struct ConsensusFeatureSchema;

/// Schema for the peptide-spectrum match (PSM) results table
/// (`psms.parquet`).
///
/// Ports `PSMSchema`. Every field is nullable. The class test's section title
/// says 25 fields while its assertion expects 29; the schema has 29, the last
/// four being the fragment annotation arrays.
#[derive(Clone, Copy, Debug)]
pub struct PSMSchema;

/// Schema for the QPX PSM export (quantms Parquet eXchange format, PSM
/// table).
///
/// Ports `QPXPSMSchema`: the column names, nested Arrow types and the
/// complete schema for peptide-spectrum match data in the QPX format, which
/// `QPXFile` uses for export.
#[derive(Clone, Copy, Debug)]
pub struct QPXPSMSchema;

/// Schema for the QPX feature view (quantms Parquet eXchange format).
///
/// Ports `QPXFeatureSchema`: the column names, nested Arrow types and the
/// complete schema for consensus feature data in the QPX format, which
/// `ConsensusMapArrowExport` uses. The shared nested types delegate to
/// [`QPXPSMSchema`].
#[derive(Clone, Copy, Debug)]
pub struct QPXFeatureSchema;

/// Schema for the QPX protein group export (quantms Parquet eXchange format,
/// pg table).
///
/// Ports `QPXPgSchema`: the column names, nested Arrow types and complete
/// schema for protein group data in the QPX format. It serves both quantified
/// (`ConsensusMap`) and identification-only (search engine) output, so the
/// quantification columns are nullable.
///
/// QPX 1.1 keys this view on `anchor_protein`, `grouped_runs` and `label`
/// (the source's `@note`): one protein-group quantity applies to one
/// experimental-design fraction group and one label. The psm and feature
/// views keep their scalar `run_file_name` and nested intensity
/// representation.
#[derive(Clone, Copy, Debug)]
pub struct QPXPgSchema;

/// Schema for spectra in long format, one row per peak.
///
/// Ports `SpectraLongSchema`.
#[derive(Clone, Copy, Debug)]
pub struct SpectraLongSchema;

/// Schema for spectra in semi-wide format, one row per spectrum with list
/// columns for the peaks.
///
/// Ports `SpectraSemiWideSchema`.
#[derive(Clone, Copy, Debug)]
pub struct SpectraSemiWideSchema;

/// Schema for chromatograms in long format, one row per data point.
///
/// Ports `ChromatogramSchema`.
#[derive(Clone, Copy, Debug)]
pub struct ChromatogramSchema;

/// Schema for chromatograms in semi-wide format, one row per chromatogram
/// with list columns.
///
/// Ports `ChromatogramSemiWideSchema`.
#[derive(Clone, Copy, Debug)]
pub struct ChromatogramSemiWideSchema;

/// Schema for the OpenSWATH precursor (peptide query) table.
///
/// Ports `OSWPrecursorSchema`.
#[derive(Clone, Copy, Debug)]
pub struct OSWPrecursorSchema;

/// Schema for the OpenSWATH transition (fragment ion) table.
///
/// Ports `OSWTransitionSchema`.
#[derive(Clone, Copy, Debug)]
pub struct OSWTransitionSchema;

/// Schema for the OpenSWATH feature-level precursor intensity table.
///
/// Ports `OSWFeaturePrecursorSchema`.
#[derive(Clone, Copy, Debug)]
pub struct OSWFeaturePrecursorSchema;

/// Schema for the OpenSWATH run metadata table.
///
/// Ports `OSWRunSchema`.
#[derive(Clone, Copy, Debug)]
pub struct OSWRunSchema;

/// Schema for the OpenSWATH feature scoring results table.
///
/// Ports `OSWFeatureSchema`.
#[derive(Clone, Copy, Debug)]
pub struct OSWFeatureSchema;

/// Schema for the OpenSWATH per-transition feature scoring results table.
///
/// Ports `OSWFeatureTransitionSchema`.
#[derive(Clone, Copy, Debug)]
pub struct OSWFeatureTransitionSchema;

/// Schema for the extracted ion chromatogram (XIC) data table.
///
/// Ports `XICSchema`. Its column names are upper case, unlike every other
/// schema here.
#[derive(Clone, Copy, Debug)]
pub struct XICSchema;

/// Schema for the extracted ion mobilogram (XIM) data table.
///
/// Ports `XIMSchema`.
#[derive(Clone, Copy, Debug)]
pub struct XIMSchema;

/// Schema for the extracted ion peak-map (XIPM) data table.
///
/// Ports `XIPMSchema`: the Arrow/Parquet column contract for targeted
/// peak-map extraction, including run and transition metadata, target
/// coordinates, encoded peak payloads and per-payload compression codes.
#[derive(Clone, Copy, Debug)]
pub struct XIPMSchema;

impl ProteinSchema {
    /// Column name `accession`.
    pub const ACCESSION: &'static str = "accession";
    /// Column name `score`.
    pub const SCORE: &'static str = "score";
    /// Column name `rank`.
    pub const RANK: &'static str = "rank";
    /// Column name `coverage`.
    pub const COVERAGE: &'static str = "coverage";
    /// Column name `sequence`.
    pub const SEQUENCE: &'static str = "sequence";
    /// Column name `description`.
    pub const DESCRIPTION: &'static str = "description";
    /// Column name `is_decoy`.
    pub const IS_DECOY: &'static str = "is_decoy";
    /// Column name `run_identifier`.
    pub const RUN_IDENTIFIER: &'static str = "run_identifier";
    /// Column name `modifications`.
    pub const MODIFICATIONS: &'static str = "modifications";
    /// Column name `metavalues`.
    pub const METAVALUES: &'static str = "metavalues";
    /// The nested column type ``list<item: struct<position: int32, modification: string>>``.
    ///
    /// Source `ProteinSchema::modificationsType`, spelled as C++ `DataType::ToString` prints it.
    pub fn modifications_type() -> DataType {
        list(struct_(vec![
            field("position", DataType::Int32, true),
            field("modification", DataType::Utf8, true),
        ]))
    }
    /// The nested column type ``list<item: struct<name: string, value: string, value_type: string>>``.
    ///
    /// Source `ProteinSchema::metavaluesType`, spelled as C++ `DataType::ToString` prints it.
    pub fn metavalues_type() -> DataType {
        list(struct_(vec![
            field("name", DataType::Utf8, true),
            field("value", DataType::Utf8, true),
            field("value_type", DataType::Utf8, true),
        ]))
    }
    /// The complete Arrow schema, 10 fields.
    ///
    /// Source `ProteinSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// accession: string not null
    /// score: double not null
    /// rank: int32 not null
    /// coverage: double
    /// sequence: string
    /// description: string
    /// is_decoy: bool
    /// run_identifier: string not null
    /// modifications: list<item: struct<position: int32, modification: string>>
    /// metavalues: list<item: struct<name: string, value: string, value_type: string>> not null
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::ACCESSION, DataType::Utf8, false),
            field(Self::SCORE, DataType::Float64, false),
            field(Self::RANK, DataType::Int32, false),
            field(Self::COVERAGE, DataType::Float64, true),
            field(Self::SEQUENCE, DataType::Utf8, true),
            field(Self::DESCRIPTION, DataType::Utf8, true),
            field(Self::IS_DECOY, DataType::Boolean, true),
            field(Self::RUN_IDENTIFIER, DataType::Utf8, false),
            field(Self::MODIFICATIONS, Self::modifications_type(), true),
            field(Self::METAVALUES, Self::metavalues_type(), false),
        ])
    }
}

impl ProteinGroupSchema {
    /// Column name `group_type`.
    pub const GROUP_TYPE: &'static str = "group_type";
    /// Column name `probability`.
    pub const PROBABILITY: &'static str = "probability";
    /// Column name `accessions`.
    pub const ACCESSIONS: &'static str = "accessions";
    /// Column name `run_identifier`.
    pub const RUN_IDENTIFIER: &'static str = "run_identifier";
    /// Column name `group_index`.
    pub const GROUP_INDEX: &'static str = "group_index";
    /// Column name `float_data`.
    pub const FLOAT_DATA: &'static str = "float_data";
    /// Column name `string_data`.
    pub const STRING_DATA: &'static str = "string_data";
    /// Column name `integer_data`.
    pub const INTEGER_DATA: &'static str = "integer_data";
    /// The nested column type ``list<item: struct<name: string, values: list<item: double>>>``.
    ///
    /// Source `ProteinGroupSchema::floatDataType`, spelled as C++ `DataType::ToString` prints it.
    pub fn float_data_type() -> DataType {
        list(struct_(vec![
            field("name", DataType::Utf8, true),
            field("values", list(DataType::Float64), true),
        ]))
    }
    /// The nested column type ``list<item: struct<name: string, values: list<item: string>>>``.
    ///
    /// Source `ProteinGroupSchema::stringDataType`, spelled as C++ `DataType::ToString` prints it.
    pub fn string_data_type() -> DataType {
        list(struct_(vec![
            field("name", DataType::Utf8, true),
            field("values", list(DataType::Utf8), true),
        ]))
    }
    /// The nested column type ``list<item: struct<name: string, values: list<item: int64>>>``.
    ///
    /// Source `ProteinGroupSchema::integerDataType`, spelled as C++ `DataType::ToString` prints it.
    pub fn integer_data_type() -> DataType {
        list(struct_(vec![
            field("name", DataType::Utf8, true),
            field("values", list(DataType::Int64), true),
        ]))
    }
    /// The complete Arrow schema, 8 fields.
    ///
    /// Source `ProteinGroupSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// group_type: string not null
    /// probability: double not null
    /// accessions: list<item: string> not null
    /// run_identifier: string not null
    /// group_index: int32 not null
    /// float_data: list<item: struct<name: string, values: list<item: double>>>
    /// string_data: list<item: struct<name: string, values: list<item: string>>>
    /// integer_data: list<item: struct<name: string, values: list<item: int64>>>
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::GROUP_TYPE, DataType::Utf8, false),
            field(Self::PROBABILITY, DataType::Float64, false),
            field(Self::ACCESSIONS, list(DataType::Utf8), false),
            field(Self::RUN_IDENTIFIER, DataType::Utf8, false),
            field(Self::GROUP_INDEX, DataType::Int32, false),
            field(Self::FLOAT_DATA, Self::float_data_type(), true),
            field(Self::STRING_DATA, Self::string_data_type(), true),
            field(Self::INTEGER_DATA, Self::integer_data_type(), true),
        ])
    }
}

impl SearchParamsSchema {
    /// Column name `run_identifier`.
    pub const RUN_IDENTIFIER: &'static str = "run_identifier";
    /// Column name `search_engine`.
    pub const SEARCH_ENGINE: &'static str = "search_engine";
    /// Column name `search_engine_version`.
    pub const SEARCH_ENGINE_VERSION: &'static str = "search_engine_version";
    /// Column name `inference_engine`.
    pub const INFERENCE_ENGINE: &'static str = "inference_engine";
    /// Column name `inference_engine_version`.
    pub const INFERENCE_ENGINE_VERSION: &'static str = "inference_engine_version";
    /// Column name `date`.
    pub const DATE: &'static str = "date";
    /// Column name `score_type`.
    pub const SCORE_TYPE: &'static str = "score_type";
    /// Column name `higher_score_better`.
    pub const HIGHER_SCORE_BETTER: &'static str = "higher_score_better";
    /// Column name `significance_threshold`.
    pub const SIGNIFICANCE_THRESHOLD: &'static str = "significance_threshold";
    /// Column name `db`.
    pub const DB: &'static str = "db";
    /// Column name `db_version`.
    pub const DB_VERSION: &'static str = "db_version";
    /// Column name `taxonomy`.
    pub const TAXONOMY: &'static str = "taxonomy";
    /// Column name `charges`.
    pub const CHARGES: &'static str = "charges";
    /// Column name `mass_type`.
    pub const MASS_TYPE: &'static str = "mass_type";
    /// Column name `precursor_mass_tolerance`.
    pub const PRECURSOR_MASS_TOLERANCE: &'static str = "precursor_mass_tolerance";
    /// Column name `precursor_mass_tolerance_ppm`.
    pub const PRECURSOR_MASS_TOLERANCE_PPM: &'static str = "precursor_mass_tolerance_ppm";
    /// Column name `fragment_mass_tolerance`.
    pub const FRAGMENT_MASS_TOLERANCE: &'static str = "fragment_mass_tolerance";
    /// Column name `fragment_mass_tolerance_ppm`.
    pub const FRAGMENT_MASS_TOLERANCE_PPM: &'static str = "fragment_mass_tolerance_ppm";
    /// Column name `digestion_enzyme`.
    pub const DIGESTION_ENZYME: &'static str = "digestion_enzyme";
    /// Column name `enzyme_term_specificity`.
    pub const ENZYME_TERM_SPECIFICITY: &'static str = "enzyme_term_specificity";
    /// Column name `missed_cleavages`.
    pub const MISSED_CLEAVAGES: &'static str = "missed_cleavages";
    /// Column name `fixed_modifications`.
    pub const FIXED_MODIFICATIONS: &'static str = "fixed_modifications";
    /// Column name `variable_modifications`.
    pub const VARIABLE_MODIFICATIONS: &'static str = "variable_modifications";
    /// Column name `primary_ms_run_paths`.
    pub const PRIMARY_MS_RUN_PATHS: &'static str = "primary_ms_run_paths";
    /// Column name `metavalues`.
    pub const METAVALUES: &'static str = "metavalues";
    /// Column name `sp_metavalues`.
    pub const SP_METAVALUES: &'static str = "sp_metavalues";
    /// The nested column type ``list<item: struct<name: string, value: string, value_type: string>>``.
    ///
    /// Source `SearchParamsSchema::metavaluesType`, spelled as C++ `DataType::ToString` prints it.
    pub fn metavalues_type() -> DataType {
        list(struct_(vec![
            field("name", DataType::Utf8, true),
            field("value", DataType::Utf8, true),
            field("value_type", DataType::Utf8, true),
        ]))
    }
    /// The complete Arrow schema, 26 fields.
    ///
    /// Source `SearchParamsSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// run_identifier: string not null
    /// search_engine: string not null
    /// search_engine_version: string
    /// inference_engine: string
    /// inference_engine_version: string
    /// date: timestamp[s]
    /// score_type: string not null
    /// higher_score_better: bool not null
    /// significance_threshold: double
    /// db: string
    /// db_version: string
    /// taxonomy: string
    /// charges: string
    /// mass_type: string not null
    /// precursor_mass_tolerance: double not null
    /// precursor_mass_tolerance_ppm: bool not null
    /// fragment_mass_tolerance: double not null
    /// fragment_mass_tolerance_ppm: bool not null
    /// digestion_enzyme: string
    /// enzyme_term_specificity: string
    /// missed_cleavages: int32 not null
    /// fixed_modifications: list<item: string> not null
    /// variable_modifications: list<item: string> not null
    /// primary_ms_run_paths: list<item: string> not null
    /// metavalues: list<item: struct<name: string, value: string, value_type: string>> not null
    /// sp_metavalues: list<item: struct<name: string, value: string, value_type: string>> not null
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::RUN_IDENTIFIER, DataType::Utf8, false),
            field(Self::SEARCH_ENGINE, DataType::Utf8, false),
            field(Self::SEARCH_ENGINE_VERSION, DataType::Utf8, true),
            field(Self::INFERENCE_ENGINE, DataType::Utf8, true),
            field(Self::INFERENCE_ENGINE_VERSION, DataType::Utf8, true),
            field(
                Self::DATE,
                DataType::Timestamp(TimeUnit::Second, None),
                true,
            ),
            field(Self::SCORE_TYPE, DataType::Utf8, false),
            field(Self::HIGHER_SCORE_BETTER, DataType::Boolean, false),
            field(Self::SIGNIFICANCE_THRESHOLD, DataType::Float64, true),
            field(Self::DB, DataType::Utf8, true),
            field(Self::DB_VERSION, DataType::Utf8, true),
            field(Self::TAXONOMY, DataType::Utf8, true),
            field(Self::CHARGES, DataType::Utf8, true),
            field(Self::MASS_TYPE, DataType::Utf8, false),
            field(Self::PRECURSOR_MASS_TOLERANCE, DataType::Float64, false),
            field(Self::PRECURSOR_MASS_TOLERANCE_PPM, DataType::Boolean, false),
            field(Self::FRAGMENT_MASS_TOLERANCE, DataType::Float64, false),
            field(Self::FRAGMENT_MASS_TOLERANCE_PPM, DataType::Boolean, false),
            field(Self::DIGESTION_ENZYME, DataType::Utf8, true),
            field(Self::ENZYME_TERM_SPECIFICITY, DataType::Utf8, true),
            field(Self::MISSED_CLEAVAGES, DataType::Int32, false),
            field(Self::FIXED_MODIFICATIONS, list(DataType::Utf8), false),
            field(Self::VARIABLE_MODIFICATIONS, list(DataType::Utf8), false),
            field(Self::PRIMARY_MS_RUN_PATHS, list(DataType::Utf8), false),
            field(Self::METAVALUES, Self::metavalues_type(), false),
            field(Self::SP_METAVALUES, Self::metavalues_type(), false),
        ])
    }
}

impl FeatureSchema {
    /// Column name `unique_id`.
    pub const UNIQUE_ID: &'static str = "unique_id";
    /// Column name `parent_feature_id`.
    pub const PARENT_FEATURE_ID: &'static str = "parent_feature_id";
    /// Column name `depth`.
    pub const DEPTH: &'static str = "depth";
    /// Column name `rt`.
    pub const RT: &'static str = "rt";
    /// Column name `mz`.
    pub const MZ: &'static str = "mz";
    /// Column name `intensity`.
    pub const INTENSITY: &'static str = "intensity";
    /// Column name `charge`.
    pub const CHARGE: &'static str = "charge";
    /// Column name `quality`.
    pub const QUALITY: &'static str = "quality";
    /// Column name `quality_rt`.
    pub const QUALITY_RT: &'static str = "quality_rt";
    /// Column name `quality_mz`.
    pub const QUALITY_MZ: &'static str = "quality_mz";
    /// Column name `width`.
    pub const WIDTH: &'static str = "width";
    /// Column name `rt_bb_min`.
    pub const RT_BB_MIN: &'static str = "rt_bb_min";
    /// Column name `rt_bb_max`.
    pub const RT_BB_MAX: &'static str = "rt_bb_max";
    /// Column name `mz_bb_min`.
    pub const MZ_BB_MIN: &'static str = "mz_bb_min";
    /// Column name `mz_bb_max`.
    pub const MZ_BB_MAX: &'static str = "mz_bb_max";
    /// Column name `convex_hulls`.
    pub const CONVEX_HULLS: &'static str = "convex_hulls";
    /// Column name `metavalues`.
    pub const METAVALUES: &'static str = "metavalues";
    /// The nested column type ``list<item: struct<hull_index: int32, points: list<item: struct<x: double, y: double>>>>``.
    ///
    /// Source `FeatureSchema::convexHullType`, spelled as C++ `DataType::ToString` prints it.
    pub fn convex_hull_type() -> DataType {
        list(struct_(vec![
            field("hull_index", DataType::Int32, true),
            field(
                "points",
                list(struct_(vec![
                    field("x", DataType::Float64, true),
                    field("y", DataType::Float64, true),
                ])),
                true,
            ),
        ]))
    }
    /// The nested column type ``list<item: struct<name: string, value: string, value_type: string>>``.
    ///
    /// Source `FeatureSchema::metavaluesType`, spelled as C++ `DataType::ToString` prints it.
    pub fn metavalues_type() -> DataType {
        list(struct_(vec![
            field("name", DataType::Utf8, true),
            field("value", DataType::Utf8, true),
            field("value_type", DataType::Utf8, true),
        ]))
    }
    /// The complete Arrow schema, 17 fields.
    ///
    /// Source `FeatureSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// unique_id: int64 not null
    /// parent_feature_id: int64
    /// depth: int32 not null
    /// rt: double not null
    /// mz: double not null
    /// intensity: float not null
    /// charge: int32 not null
    /// quality: float not null
    /// quality_rt: float not null
    /// quality_mz: float not null
    /// width: float
    /// rt_bb_min: double
    /// rt_bb_max: double
    /// mz_bb_min: double
    /// mz_bb_max: double
    /// convex_hulls: list<item: struct<hull_index: int32, points: list<item: struct<x: double, y: double>>>> not null
    /// metavalues: list<item: struct<name: string, value: string, value_type: string>> not null
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::UNIQUE_ID, DataType::Int64, false),
            field(Self::PARENT_FEATURE_ID, DataType::Int64, true),
            field(Self::DEPTH, DataType::Int32, false),
            field(Self::RT, DataType::Float64, false),
            field(Self::MZ, DataType::Float64, false),
            field(Self::INTENSITY, DataType::Float32, false),
            field(Self::CHARGE, DataType::Int32, false),
            field(Self::QUALITY, DataType::Float32, false),
            field(Self::QUALITY_RT, DataType::Float32, false),
            field(Self::QUALITY_MZ, DataType::Float32, false),
            field(Self::WIDTH, DataType::Float32, true),
            field(Self::RT_BB_MIN, DataType::Float64, true),
            field(Self::RT_BB_MAX, DataType::Float64, true),
            field(Self::MZ_BB_MIN, DataType::Float64, true),
            field(Self::MZ_BB_MAX, DataType::Float64, true),
            field(Self::CONVEX_HULLS, Self::convex_hull_type(), false),
            field(Self::METAVALUES, Self::metavalues_type(), false),
        ])
    }
}

impl ConsensusFeatureSchema {
    /// Column name `unique_id`.
    pub const UNIQUE_ID: &'static str = "unique_id";
    /// Column name `rt`.
    pub const RT: &'static str = "rt";
    /// Column name `mz`.
    pub const MZ: &'static str = "mz";
    /// Column name `intensity`.
    pub const INTENSITY: &'static str = "intensity";
    /// Column name `charge`.
    pub const CHARGE: &'static str = "charge";
    /// Column name `quality`.
    pub const QUALITY: &'static str = "quality";
    /// Column name `width`.
    pub const WIDTH: &'static str = "width";
    /// Column name `handles`.
    pub const HANDLES: &'static str = "handles";
    /// Column name `metavalues`.
    pub const METAVALUES: &'static str = "metavalues";
    /// The nested column type ``list<item: struct<map_index: int64, unique_id: int64, rt: double, mz: double, intensity: float, charge: int32, width: float>>``.
    ///
    /// Source `ConsensusFeatureSchema::handlesType`, spelled as C++ `DataType::ToString` prints it.
    pub fn handles_type() -> DataType {
        list(struct_(vec![
            field("map_index", DataType::Int64, true),
            field("unique_id", DataType::Int64, true),
            field("rt", DataType::Float64, true),
            field("mz", DataType::Float64, true),
            field("intensity", DataType::Float32, true),
            field("charge", DataType::Int32, true),
            field("width", DataType::Float32, true),
        ]))
    }
    /// The nested column type ``list<item: struct<name: string, value: string, value_type: string>>``.
    ///
    /// Source `ConsensusFeatureSchema::metavaluesType`, spelled as C++ `DataType::ToString` prints it.
    pub fn metavalues_type() -> DataType {
        list(struct_(vec![
            field("name", DataType::Utf8, true),
            field("value", DataType::Utf8, true),
            field("value_type", DataType::Utf8, true),
        ]))
    }
    /// The complete Arrow schema, 9 fields.
    ///
    /// Source `ConsensusFeatureSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// unique_id: int64
    /// rt: double
    /// mz: double
    /// intensity: float
    /// charge: int32
    /// quality: float
    /// width: float
    /// handles: list<item: struct<map_index: int64, unique_id: int64, rt: double, mz: double, intensity: float, charge: int32, width: float>>
    /// metavalues: list<item: struct<name: string, value: string, value_type: string>>
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::UNIQUE_ID, DataType::Int64, true),
            field(Self::RT, DataType::Float64, true),
            field(Self::MZ, DataType::Float64, true),
            field(Self::INTENSITY, DataType::Float32, true),
            field(Self::CHARGE, DataType::Int32, true),
            field(Self::QUALITY, DataType::Float32, true),
            field(Self::WIDTH, DataType::Float32, true),
            field(Self::HANDLES, Self::handles_type(), true),
            field(Self::METAVALUES, Self::metavalues_type(), true),
        ])
    }
}

impl PSMSchema {
    /// Column name `sequence`.
    pub const SEQUENCE: &'static str = "sequence";
    /// Column name `peptidoform`.
    pub const PEPTIDOFORM: &'static str = "peptidoform";
    /// Column name `modifications`.
    pub const MODIFICATIONS: &'static str = "modifications";
    /// Column name `precursor_charge`.
    pub const PRECURSOR_CHARGE: &'static str = "precursor_charge";
    /// Column name `posterior_error_probability`.
    pub const POSTERIOR_ERROR_PROBABILITY: &'static str = "posterior_error_probability";
    /// Column name `is_decoy`.
    pub const IS_DECOY: &'static str = "is_decoy";
    /// Column name `calculated_mz`.
    pub const CALCULATED_MZ: &'static str = "calculated_mz";
    /// Column name `observed_mz`.
    pub const OBSERVED_MZ: &'static str = "observed_mz";
    /// Column name `additional_scores`.
    pub const ADDITIONAL_SCORES: &'static str = "additional_scores";
    /// Column name `protein_accessions`.
    pub const PROTEIN_ACCESSIONS: &'static str = "protein_accessions";
    /// Column name `predicted_rt`.
    pub const PREDICTED_RT: &'static str = "predicted_rt";
    /// Column name `reference_file_name`.
    pub const REFERENCE_FILE_NAME: &'static str = "reference_file_name";
    /// Column name `cv_params`.
    pub const CV_PARAMS: &'static str = "cv_params";
    /// Column name `scan`.
    pub const SCAN: &'static str = "scan";
    /// Column name `rt`.
    pub const RT: &'static str = "rt";
    /// Column name `ion_mobility`.
    pub const ION_MOBILITY: &'static str = "ion_mobility";
    /// Column name `spectrum_reference`.
    pub const SPECTRUM_REFERENCE: &'static str = "spectrum_reference";
    /// Column name `score`.
    pub const SCORE: &'static str = "score";
    /// Column name `score_type`.
    pub const SCORE_TYPE: &'static str = "score_type";
    /// Column name `higher_score_better`.
    pub const HIGHER_SCORE_BETTER: &'static str = "higher_score_better";
    /// 0-based positional column: hit position within parent identification.
    /// Distinct from `PeptideHit::getRank()`, which is round-tripped via the
    /// "rank" UserParam carried in PSM_METAVALUES.
    pub const HIT_INDEX: &'static str = "hit_index";
    /// Column name `peptide_identification_index`.
    pub const PEPTIDE_IDENTIFICATION_INDEX: &'static str = "peptide_identification_index";
    /// Column name `psm_metavalues`.
    pub const PSM_METAVALUES: &'static str = "psm_metavalues";
    /// Column name `spectrum_metavalues`.
    pub const SPECTRUM_METAVALUES: &'static str = "spectrum_metavalues";
    /// Column name `run_identifier`.
    pub const RUN_IDENTIFIER: &'static str = "run_identifier";
    /// Column name `mz_array`.
    pub const MZ_ARRAY: &'static str = "mz_array";
    /// Column name `intensity_array`.
    pub const INTENSITY_ARRAY: &'static str = "intensity_array";
    /// Column name `charge_array`.
    pub const CHARGE_ARRAY: &'static str = "charge_array";
    /// Column name `ion_type_array`.
    pub const ION_TYPE_ARRAY: &'static str = "ion_type_array";
    /// The nested column type ``list<item: struct<name: string, accession: string, positions: list<item: struct<position: string, scores: double>>>>``.
    ///
    /// Source `PSMSchema::modificationsType`, spelled as C++ `DataType::ToString` prints it.
    pub fn modifications_type() -> DataType {
        list(struct_(vec![
            field("name", DataType::Utf8, true),
            field("accession", DataType::Utf8, true),
            field(
                "positions",
                list(struct_(vec![
                    field("position", DataType::Utf8, true),
                    field("scores", DataType::Float64, true),
                ])),
                true,
            ),
        ]))
    }
    /// The nested column type ``list<item: struct<score_name: string, score_value: double, higher_better: bool>>``.
    ///
    /// Source `PSMSchema::additionalScoresType`, spelled as C++ `DataType::ToString` prints it.
    pub fn additional_scores_type() -> DataType {
        list(struct_(vec![
            field("score_name", DataType::Utf8, true),
            field("score_value", DataType::Float64, true),
            field("higher_better", DataType::Boolean, true),
        ]))
    }
    /// The nested column type ``list<item: struct<name: string, value: string, value_type: string>>``.
    ///
    /// Source `PSMSchema::metavaluesType`, spelled as C++ `DataType::ToString` prints it.
    pub fn metavalues_type() -> DataType {
        list(struct_(vec![
            field("name", DataType::Utf8, true),
            field("value", DataType::Utf8, true),
            field("value_type", DataType::Utf8, true),
        ]))
    }
    /// Arrow type for protein_accessions: `list<struct{accession, aa_before, aa_after, start, end}>`.
    ///
    /// Source `PSMSchema::proteinAccessionsType`, spelled as C++ `DataType::ToString` prints it.
    pub fn protein_accessions_type() -> DataType {
        list(struct_(vec![
            field("accession", DataType::Utf8, false),
            field("aa_before", DataType::Utf8, true),
            field("aa_after", DataType::Utf8, true),
            field("start", DataType::Int32, true),
            field("end", DataType::Int32, true),
        ]))
    }
    /// The complete Arrow schema, 29 fields.
    ///
    /// Source `PSMSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// sequence: string
    /// peptidoform: string
    /// modifications: list<item: struct<name: string, accession: string, positions: list<item: struct<position: string, scores: double>>>>
    /// precursor_charge: int32
    /// posterior_error_probability: double
    /// is_decoy: bool
    /// calculated_mz: double
    /// observed_mz: double
    /// additional_scores: list<item: struct<score_name: string, score_value: double, higher_better: bool>>
    /// protein_accessions: list<item: struct<accession: string not null, aa_before: string, aa_after: string, start: int32, end: int32>>
    /// predicted_rt: double
    /// reference_file_name: string
    /// cv_params: string
    /// scan: int32
    /// rt: double
    /// ion_mobility: double
    /// spectrum_reference: string
    /// score: double
    /// score_type: string
    /// higher_score_better: bool
    /// hit_index: int32
    /// peptide_identification_index: int32
    /// psm_metavalues: list<item: struct<name: string, value: string, value_type: string>>
    /// spectrum_metavalues: list<item: struct<name: string, value: string, value_type: string>>
    /// run_identifier: string
    /// mz_array: list<item: float>
    /// intensity_array: list<item: float>
    /// charge_array: list<item: int32>
    /// ion_type_array: list<item: string>
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::SEQUENCE, DataType::Utf8, true),
            field(Self::PEPTIDOFORM, DataType::Utf8, true),
            field(Self::MODIFICATIONS, Self::modifications_type(), true),
            field(Self::PRECURSOR_CHARGE, DataType::Int32, true),
            field(Self::POSTERIOR_ERROR_PROBABILITY, DataType::Float64, true),
            field(Self::IS_DECOY, DataType::Boolean, true),
            field(Self::CALCULATED_MZ, DataType::Float64, true),
            field(Self::OBSERVED_MZ, DataType::Float64, true),
            field(
                Self::ADDITIONAL_SCORES,
                Self::additional_scores_type(),
                true,
            ),
            field(
                Self::PROTEIN_ACCESSIONS,
                Self::protein_accessions_type(),
                true,
            ),
            field(Self::PREDICTED_RT, DataType::Float64, true),
            field(Self::REFERENCE_FILE_NAME, DataType::Utf8, true),
            field(Self::CV_PARAMS, DataType::Utf8, true),
            field(Self::SCAN, DataType::Int32, true),
            field(Self::RT, DataType::Float64, true),
            field(Self::ION_MOBILITY, DataType::Float64, true),
            field(Self::SPECTRUM_REFERENCE, DataType::Utf8, true),
            field(Self::SCORE, DataType::Float64, true),
            field(Self::SCORE_TYPE, DataType::Utf8, true),
            field(Self::HIGHER_SCORE_BETTER, DataType::Boolean, true),
            field(Self::HIT_INDEX, DataType::Int32, true),
            field(Self::PEPTIDE_IDENTIFICATION_INDEX, DataType::Int32, true),
            field(Self::PSM_METAVALUES, Self::metavalues_type(), true),
            field(Self::SPECTRUM_METAVALUES, Self::metavalues_type(), true),
            field(Self::RUN_IDENTIFIER, DataType::Utf8, true),
            field(Self::MZ_ARRAY, list(DataType::Float32), true),
            field(Self::INTENSITY_ARRAY, list(DataType::Float32), true),
            field(Self::CHARGE_ARRAY, list(DataType::Int32), true),
            field(Self::ION_TYPE_ARRAY, list(DataType::Utf8), true),
        ])
    }
}

impl QPXPSMSchema {
    /// Mandatory opaque identity and primary key; see QPXIdentity.
    pub const PSM_ID: &'static str = "psm_id";
    /// Optional cross-reference to `QPXFeatureSchema::FEATURE_ID`; null when the PSM is unlinked.
    pub const FEATURE_ID: &'static str = "feature_id";
    /// Column name `sequence`.
    pub const SEQUENCE: &'static str = "sequence";
    /// Column name `peptidoform`.
    pub const PEPTIDOFORM: &'static str = "peptidoform";
    /// Column name `modifications`.
    pub const MODIFICATIONS: &'static str = "modifications";
    /// Column name `charge`.
    pub const CHARGE: &'static str = "charge";
    /// Column name `posterior_error_probability`.
    pub const POSTERIOR_ERROR_PROBABILITY: &'static str = "posterior_error_probability";
    /// Column name `is_decoy`.
    pub const IS_DECOY: &'static str = "is_decoy";
    /// Column name `calculated_mz`.
    pub const CALCULATED_MZ: &'static str = "calculated_mz";
    /// Column name `observed_mz`.
    pub const OBSERVED_MZ: &'static str = "observed_mz";
    /// Column name `mass_error_ppm`.
    pub const MASS_ERROR_PPM: &'static str = "mass_error_ppm";
    /// Column name `additional_scores`.
    pub const ADDITIONAL_SCORES: &'static str = "additional_scores";
    /// Column name `predicted_rt`.
    pub const PREDICTED_RT: &'static str = "predicted_rt";
    /// Column name `run_file_name`.
    pub const RUN_FILE_NAME: &'static str = "run_file_name";
    /// Column name `cv_params`.
    pub const CV_PARAMS: &'static str = "cv_params";
    /// Column name `scan`.
    pub const SCAN: &'static str = "scan";
    /// Column name `rt`.
    pub const RT: &'static str = "rt";
    /// Column name `ion_mobility`.
    pub const ION_MOBILITY: &'static str = "ion_mobility";
    /// Column name `missed_cleavages`.
    pub const MISSED_CLEAVAGES: &'static str = "missed_cleavages";
    /// Column name `protein_accessions`.
    pub const PROTEIN_ACCESSIONS: &'static str = "protein_accessions";
    /// Column name `cross_links`.
    pub const CROSS_LINKS: &'static str = "cross_links";
    /// Column name `mz_array`.
    pub const MZ_ARRAY: &'static str = "mz_array";
    /// Column name `intensity_array`.
    pub const INTENSITY_ARRAY: &'static str = "intensity_array";
    /// Column name `charge_array`.
    pub const CHARGE_ARRAY: &'static str = "charge_array";
    /// Column name `ion_type_array`.
    pub const ION_TYPE_ARRAY: &'static str = "ion_type_array";
    /// Column name `ion_mobility_array`.
    pub const ION_MOBILITY_ARRAY: &'static str = "ion_mobility_array";
    /// Arrow type for modifications: `list<struct{name, accession, positions: list<struct{position, amino_acid, scores}>}>`.
    ///
    /// Source `QPXPSMSchema::modificationsType`, spelled as C++ `DataType::ToString` prints it.
    pub fn modifications_type() -> DataType {
        list(struct_(vec![
            field("name", DataType::Utf8, false),
            field("accession", DataType::Utf8, true),
            field(
                "positions",
                list(struct_(vec![
                    field("position", DataType::Int32, false),
                    field("amino_acid", DataType::Utf8, true),
                    field(
                        "scores",
                        list(struct_(vec![
                            field("score_name", DataType::Utf8, false),
                            field("score_value", DataType::Float64, false),
                            field("higher_better", DataType::Boolean, true),
                        ])),
                        true,
                    ),
                ])),
                false,
            ),
        ]))
    }
    /// Arrow type for additional scores: `list<struct{score_name, score_value, higher_better}>`.
    ///
    /// Source `QPXPSMSchema::additionalScoresType`, spelled as C++ `DataType::ToString` prints it.
    pub fn additional_scores_type() -> DataType {
        list(struct_(vec![
            field("score_name", DataType::Utf8, false),
            field("score_value", DataType::Float64, false),
            field("higher_better", DataType::Boolean, true),
        ]))
    }
    /// Arrow type for CV params: `list<struct{cv_name, cv_value}>`.
    ///
    /// Source `QPXPSMSchema::cvParamsType`, spelled as C++ `DataType::ToString` prints it.
    pub fn cv_params_type() -> DataType {
        list(struct_(vec![
            field("cv_name", DataType::Utf8, false),
            field("cv_value", DataType::Utf8, false),
        ]))
    }
    /// Arrow type for cross-links: `list<struct{xl_type, partner_sequence, ...}>`.
    ///
    /// Source `QPXPSMSchema::crossLinksType`, spelled as C++ `DataType::ToString` prints it.
    pub fn cross_links_type() -> DataType {
        list(struct_(vec![
            field("xl_type", DataType::Utf8, false),
            field("partner_sequence", DataType::Utf8, true),
            field("partner_peptidoform", DataType::Utf8, true),
            field("donor_position", DataType::Int32, false),
            field("acceptor_position", DataType::Int32, true),
            field("linker_name", DataType::Utf8, false),
            field("linker_accession", DataType::Utf8, true),
            field("linker_mass", DataType::Float64, false),
        ]))
    }
    /// Complete Arrow schema for QPX PSM table (26 fields).
    ///
    /// Source `QPXPSMSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// psm_id: int64 not null
    /// sequence: string not null
    /// peptidoform: string not null
    /// modifications: list<item: struct<name: string not null, accession: string, positions: list<item: struct<position: int32 not null, amino_acid: string, scores: list<item: struct<score_name: string not null, score_value: double not null, higher_better: bool>>>> not null>>
    /// charge: int16 not null
    /// posterior_error_probability: double
    /// is_decoy: bool not null
    /// calculated_mz: float not null
    /// observed_mz: float not null
    /// mass_error_ppm: float
    /// additional_scores: list<item: struct<score_name: string not null, score_value: double not null, higher_better: bool>>
    /// predicted_rt: float
    /// run_file_name: string not null
    /// cv_params: list<item: struct<cv_name: string not null, cv_value: string not null>>
    /// scan: list<item: int32> not null
    /// rt: float
    /// ion_mobility: float
    /// missed_cleavages: int16
    /// protein_accessions: list<item: string>
    /// cross_links: list<item: struct<xl_type: string not null, partner_sequence: string, partner_peptidoform: string, donor_position: int32 not null, acceptor_position: int32, linker_name: string not null, linker_accession: string, linker_mass: double not null>>
    /// mz_array: list<item: float>
    /// intensity_array: list<item: float>
    /// charge_array: list<item: int32>
    /// ion_type_array: list<item: string>
    /// ion_mobility_array: list<item: float>
    /// feature_id: int64
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::PSM_ID, DataType::Int64, false),
            field(Self::SEQUENCE, DataType::Utf8, false),
            field(Self::PEPTIDOFORM, DataType::Utf8, false),
            field(Self::MODIFICATIONS, Self::modifications_type(), true),
            field(Self::CHARGE, DataType::Int16, false),
            field(Self::POSTERIOR_ERROR_PROBABILITY, DataType::Float64, true),
            field(Self::IS_DECOY, DataType::Boolean, false),
            field(Self::CALCULATED_MZ, DataType::Float32, false),
            field(Self::OBSERVED_MZ, DataType::Float32, false),
            field(Self::MASS_ERROR_PPM, DataType::Float32, true),
            field(
                Self::ADDITIONAL_SCORES,
                Self::additional_scores_type(),
                true,
            ),
            field(Self::PREDICTED_RT, DataType::Float32, true),
            field(Self::RUN_FILE_NAME, DataType::Utf8, false),
            field(Self::CV_PARAMS, Self::cv_params_type(), true),
            field(Self::SCAN, list(DataType::Int32), false),
            field(Self::RT, DataType::Float32, true),
            field(Self::ION_MOBILITY, DataType::Float32, true),
            field(Self::MISSED_CLEAVAGES, DataType::Int16, true),
            field(Self::PROTEIN_ACCESSIONS, list(DataType::Utf8), true),
            field(Self::CROSS_LINKS, Self::cross_links_type(), true),
            field(Self::MZ_ARRAY, list(DataType::Float32), true),
            field(Self::INTENSITY_ARRAY, list(DataType::Float32), true),
            field(Self::CHARGE_ARRAY, list(DataType::Int32), true),
            field(Self::ION_TYPE_ARRAY, list(DataType::Utf8), true),
            field(Self::ION_MOBILITY_ARRAY, list(DataType::Float32), true),
            field(Self::FEATURE_ID, DataType::Int64, true),
        ])
    }
}

impl QPXFeatureSchema {
    /// Mandatory opaque identity and primary key; see QPXIdentity.
    pub const FEATURE_ID: &'static str = "feature_id";
    /// Optional cross-reference to the `QPXPSMSchema::PSM_ID` values this feature was quantified from.
    pub const PSM_IDS: &'static str = "psm_ids";
    /// Column name `sequence`.
    pub const SEQUENCE: &'static str = "sequence";
    /// Column name `peptidoform`.
    pub const PEPTIDOFORM: &'static str = "peptidoform";
    /// Column name `modifications`.
    pub const MODIFICATIONS: &'static str = "modifications";
    /// Column name `charge`.
    pub const CHARGE: &'static str = "charge";
    /// Column name `posterior_error_probability`.
    pub const POSTERIOR_ERROR_PROBABILITY: &'static str = "posterior_error_probability";
    /// Column name `is_decoy`.
    pub const IS_DECOY: &'static str = "is_decoy";
    /// Column name `calculated_mz`.
    pub const CALCULATED_MZ: &'static str = "calculated_mz";
    /// Column name `observed_mz`.
    pub const OBSERVED_MZ: &'static str = "observed_mz";
    /// Column name `mass_error_ppm`.
    pub const MASS_ERROR_PPM: &'static str = "mass_error_ppm";
    /// Column name `additional_scores`.
    pub const ADDITIONAL_SCORES: &'static str = "additional_scores";
    /// Column name `predicted_rt`.
    pub const PREDICTED_RT: &'static str = "predicted_rt";
    /// Column name `run_file_name`.
    pub const RUN_FILE_NAME: &'static str = "run_file_name";
    /// Column name `cv_params`.
    pub const CV_PARAMS: &'static str = "cv_params";
    /// Column name `scan`.
    pub const SCAN: &'static str = "scan";
    /// Column name `rt`.
    pub const RT: &'static str = "rt";
    /// Column name `ion_mobility`.
    pub const ION_MOBILITY: &'static str = "ion_mobility";
    /// Column name `missed_cleavages`.
    pub const MISSED_CLEAVAGES: &'static str = "missed_cleavages";
    /// Column name `intensities`.
    pub const INTENSITIES: &'static str = "intensities";
    /// Column name `additional_intensities`.
    pub const ADDITIONAL_INTENSITIES: &'static str = "additional_intensities";
    /// Column name `pg_accessions`.
    pub const PG_ACCESSIONS: &'static str = "pg_accessions";
    /// Column name `anchor_protein`.
    pub const ANCHOR_PROTEIN: &'static str = "anchor_protein";
    /// Column name `unique`.
    pub const UNIQUE: &'static str = "unique";
    /// Column name `pg_global_qvalue`.
    pub const PG_GLOBAL_QVALUE: &'static str = "pg_global_qvalue";
    /// Column name `pg_positions`.
    pub const PG_POSITIONS: &'static str = "pg_positions";
    /// Column name `ion_mobility_start`.
    pub const ION_MOBILITY_START: &'static str = "ion_mobility_start";
    /// Column name `ion_mobility_stop`.
    pub const ION_MOBILITY_STOP: &'static str = "ion_mobility_stop";
    /// Column name `gg_accessions`.
    pub const GG_ACCESSIONS: &'static str = "gg_accessions";
    /// Column name `gg_names`.
    pub const GG_NAMES: &'static str = "gg_names";
    /// Column name `id_run_file_name`.
    pub const ID_RUN_FILE_NAME: &'static str = "id_run_file_name";
    /// Column name `rt_start`.
    pub const RT_START: &'static str = "rt_start";
    /// Column name `rt_stop`.
    pub const RT_STOP: &'static str = "rt_stop";
    /// Arrow type for modifications (delegates to `QPXPSMSchema::modificationsType`).
    ///
    /// Source `QPXFeatureSchema::modificationsType`, spelled as C++ `DataType::ToString` prints it.
    pub fn modifications_type() -> DataType {
        QPXPSMSchema::modifications_type()
    }
    /// Arrow type for additional scores (delegates to `QPXPSMSchema::additionalScoresType`).
    ///
    /// Source `QPXFeatureSchema::additionalScoresType`, spelled as C++ `DataType::ToString` prints it.
    pub fn additional_scores_type() -> DataType {
        QPXPSMSchema::additional_scores_type()
    }
    /// Arrow type for CV params (delegates to `QPXPSMSchema::cvParamsType`).
    ///
    /// Source `QPXFeatureSchema::cvParamsType`, spelled as C++ `DataType::ToString` prints it.
    pub fn cv_params_type() -> DataType {
        QPXPSMSchema::cv_params_type()
    }
    /// Arrow type for intensities: `list<struct{label, intensity}>`.
    ///
    /// Source `QPXFeatureSchema::intensitiesType`, spelled as C++ `DataType::ToString` prints it.
    pub fn intensities_type() -> DataType {
        list(struct_(vec![
            field("label", DataType::Utf8, false),
            field("intensity", DataType::Float32, false),
        ]))
    }
    /// Arrow type for additional intensities: `list<struct{label, intensities: list<struct{...}>}>`.
    ///
    /// Source `QPXFeatureSchema::additionalIntensitiesType`, spelled as C++ `DataType::ToString` prints it.
    pub fn additional_intensities_type() -> DataType {
        list(struct_(vec![
            field("label", DataType::Utf8, false),
            field(
                "intensities",
                list(struct_(vec![
                    field("intensity_name", DataType::Utf8, false),
                    field("intensity_value", DataType::Float32, false),
                ])),
                false,
            ),
        ]))
    }
    /// Arrow type for protein group accessions: `list<struct{accession, start, end, pre, post}>`.
    ///
    /// Source `QPXFeatureSchema::pgAccessionsType`, spelled as C++ `DataType::ToString` prints it.
    pub fn pg_accessions_type() -> DataType {
        list(struct_(vec![
            field("accession", DataType::Utf8, false),
            field("start", DataType::Int32, true),
            field("end", DataType::Int32, true),
            field("pre", DataType::Utf8, true),
            field("post", DataType::Utf8, true),
        ]))
    }
    /// Arrow type for protein group positions: `list<struct{protein_accession, start, end}>`.
    ///
    /// Source `QPXFeatureSchema::pgPositionsType`, spelled as C++ `DataType::ToString` prints it.
    pub fn pg_positions_type() -> DataType {
        list(struct_(vec![
            field("protein_accession", DataType::Utf8, false),
            field("start", DataType::Int32, false),
            field("end", DataType::Int32, false),
        ]))
    }
    /// Complete Arrow schema for QPX feature table (33 fields).
    ///
    /// Source `QPXFeatureSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// feature_id: int64 not null
    /// sequence: string not null
    /// peptidoform: string not null
    /// modifications: list<item: struct<name: string not null, accession: string, positions: list<item: struct<position: int32 not null, amino_acid: string, scores: list<item: struct<score_name: string not null, score_value: double not null, higher_better: bool>>>> not null>>
    /// charge: int16 not null
    /// posterior_error_probability: double
    /// is_decoy: bool not null
    /// calculated_mz: float not null
    /// observed_mz: float not null
    /// mass_error_ppm: float
    /// additional_scores: list<item: struct<score_name: string not null, score_value: double not null, higher_better: bool>>
    /// predicted_rt: float
    /// run_file_name: string not null
    /// cv_params: list<item: struct<cv_name: string not null, cv_value: string not null>>
    /// scan: list<item: int32> not null
    /// rt: float
    /// ion_mobility: float
    /// missed_cleavages: int16
    /// intensities: list<item: struct<label: string not null, intensity: float not null>>
    /// additional_intensities: list<item: struct<label: string not null, intensities: list<item: struct<intensity_name: string not null, intensity_value: float not null>> not null>>
    /// pg_accessions: list<item: struct<accession: string not null, start: int32, end: int32, pre: string, post: string>>
    /// anchor_protein: string
    /// unique: bool
    /// pg_global_qvalue: double
    /// pg_positions: list<item: struct<protein_accession: string not null, start: int32 not null, end: int32 not null>>
    /// ion_mobility_start: float
    /// ion_mobility_stop: float
    /// gg_accessions: list<item: string>
    /// gg_names: list<item: string>
    /// id_run_file_name: string
    /// rt_start: float
    /// rt_stop: float
    /// psm_ids: list<item: int64>
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::FEATURE_ID, DataType::Int64, false),
            field(Self::SEQUENCE, DataType::Utf8, false),
            field(Self::PEPTIDOFORM, DataType::Utf8, false),
            field(Self::MODIFICATIONS, Self::modifications_type(), true),
            field(Self::CHARGE, DataType::Int16, false),
            field(Self::POSTERIOR_ERROR_PROBABILITY, DataType::Float64, true),
            field(Self::IS_DECOY, DataType::Boolean, false),
            field(Self::CALCULATED_MZ, DataType::Float32, false),
            field(Self::OBSERVED_MZ, DataType::Float32, false),
            field(Self::MASS_ERROR_PPM, DataType::Float32, true),
            field(
                Self::ADDITIONAL_SCORES,
                Self::additional_scores_type(),
                true,
            ),
            field(Self::PREDICTED_RT, DataType::Float32, true),
            field(Self::RUN_FILE_NAME, DataType::Utf8, false),
            field(Self::CV_PARAMS, Self::cv_params_type(), true),
            field(Self::SCAN, list(DataType::Int32), false),
            field(Self::RT, DataType::Float32, true),
            field(Self::ION_MOBILITY, DataType::Float32, true),
            field(Self::MISSED_CLEAVAGES, DataType::Int16, true),
            field(Self::INTENSITIES, Self::intensities_type(), true),
            field(
                Self::ADDITIONAL_INTENSITIES,
                Self::additional_intensities_type(),
                true,
            ),
            field(Self::PG_ACCESSIONS, Self::pg_accessions_type(), true),
            field(Self::ANCHOR_PROTEIN, DataType::Utf8, true),
            field(Self::UNIQUE, DataType::Boolean, true),
            field(Self::PG_GLOBAL_QVALUE, DataType::Float64, true),
            field(Self::PG_POSITIONS, Self::pg_positions_type(), true),
            field(Self::ION_MOBILITY_START, DataType::Float32, true),
            field(Self::ION_MOBILITY_STOP, DataType::Float32, true),
            field(Self::GG_ACCESSIONS, list(DataType::Utf8), true),
            field(Self::GG_NAMES, list(DataType::Utf8), true),
            field(Self::ID_RUN_FILE_NAME, DataType::Utf8, true),
            field(Self::RT_START, DataType::Float32, true),
            field(Self::RT_STOP, DataType::Float32, true),
            field(Self::PSM_IDS, list(DataType::Int64), true),
        ])
    }
}

impl QPXPgSchema {
    /// Mandatory opaque identity and primary key; see QPXIdentity.
    pub const PG_ID: &'static str = "pg_id";
    /// Column name `pg_accessions`.
    pub const PG_ACCESSIONS: &'static str = "pg_accessions";
    /// Column name `pg_names`.
    pub const PG_NAMES: &'static str = "pg_names";
    /// Column name `gg_accessions`.
    pub const GG_ACCESSIONS: &'static str = "gg_accessions";
    /// Column name `gg_names`.
    pub const GG_NAMES: &'static str = "gg_names";
    /// Column name `gg_qvalue`.
    pub const GG_QVALUE: &'static str = "gg_qvalue";
    /// Column name `anchor_protein`.
    pub const ANCHOR_PROTEIN: &'static str = "anchor_protein";
    /// Column name `grouped_runs`.
    pub const GROUPED_RUNS: &'static str = "grouped_runs";
    /// Column name `global_qvalue`.
    pub const GLOBAL_QVALUE: &'static str = "global_qvalue";
    /// Column name `pg_qvalue`.
    pub const PG_QVALUE: &'static str = "pg_qvalue";
    /// Column name `label`.
    pub const LABEL: &'static str = "label";
    /// Column name `intensity`.
    pub const INTENSITY: &'static str = "intensity";
    /// Column name `additional_intensities`.
    pub const ADDITIONAL_INTENSITIES: &'static str = "additional_intensities";
    /// Column name `is_decoy`.
    pub const IS_DECOY: &'static str = "is_decoy";
    /// Column name `contaminant`.
    pub const CONTAMINANT: &'static str = "contaminant";
    /// Column name `peptides`.
    pub const PEPTIDES: &'static str = "peptides";
    /// Column name `peptide_counts`.
    pub const PEPTIDE_COUNTS: &'static str = "peptide_counts";
    /// Column name `feature_counts`.
    pub const FEATURE_COUNTS: &'static str = "feature_counts";
    /// Column name `sequence_coverage`.
    pub const SEQUENCE_COVERAGE: &'static str = "sequence_coverage";
    /// Column name `molecular_weight`.
    pub const MOLECULAR_WEIGHT: &'static str = "molecular_weight";
    /// Column name `additional_scores`.
    pub const ADDITIONAL_SCORES: &'static str = "additional_scores";
    /// Column name `cv_params`.
    pub const CV_PARAMS: &'static str = "cv_params";
    /// Arrow type for grouped_runs: `list<utf8>` — the raw files of one quantification unit.
    ///
    /// Source `QPXPgSchema::groupedRunsType`, spelled as C++ `DataType::ToString` prints it.
    pub fn grouped_runs_type() -> DataType {
        list(DataType::Utf8)
    }
    /// Arrow type for additional intensities: `list<struct{label, intensities: list<struct{...}>}>`.
    ///
    /// Source `QPXPgSchema::additionalIntensitiesType`, spelled as C++ `DataType::ToString` prints it.
    pub fn additional_intensities_type() -> DataType {
        list(struct_(vec![
            field("label", DataType::Utf8, false),
            field(
                "intensities",
                list(struct_(vec![
                    field("intensity_name", DataType::Utf8, false),
                    field("intensity_value", DataType::Float32, false),
                ])),
                false,
            ),
        ]))
    }
    /// Arrow type for peptides: `list<struct{protein_name, peptide_count}>`.
    ///
    /// Source `QPXPgSchema::peptidesType`, spelled as C++ `DataType::ToString` prints it.
    pub fn peptides_type() -> DataType {
        list(struct_(vec![
            field("protein_name", DataType::Utf8, false),
            field("peptide_count", DataType::Int32, false),
        ]))
    }
    /// Arrow type for peptide_counts: struct{unique_sequences, total_sequences}.
    ///
    /// Source `QPXPgSchema::peptideCountsType`, spelled as C++ `DataType::ToString` prints it.
    pub fn peptide_counts_type() -> DataType {
        struct_(vec![
            field("unique_sequences", DataType::Int32, false),
            field("total_sequences", DataType::Int32, false),
        ])
    }
    /// Arrow type for feature_counts: struct{unique_features, total_features}.
    ///
    /// Source `QPXPgSchema::featureCountsType`, spelled as C++ `DataType::ToString` prints it.
    pub fn feature_counts_type() -> DataType {
        struct_(vec![
            field("unique_features", DataType::Int32, false),
            field("total_features", DataType::Int32, false),
        ])
    }
    /// Arrow type for additional scores (delegates to `QPXPSMSchema::additionalScoresType`).
    ///
    /// Source `QPXPgSchema::additionalScoresType`, spelled as C++ `DataType::ToString` prints it.
    pub fn additional_scores_type() -> DataType {
        QPXPSMSchema::additional_scores_type()
    }
    /// Arrow type for CV params (delegates to `QPXPSMSchema::cvParamsType`).
    ///
    /// Source `QPXPgSchema::cvParamsType`, spelled as C++ `DataType::ToString` prints it.
    pub fn cv_params_type() -> DataType {
        QPXPSMSchema::cv_params_type()
    }
    /// Complete Arrow schema for QPX pg table (22 fields).
    ///
    /// Source `QPXPgSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// pg_id: int64 not null
    /// pg_accessions: list<item: string> not null
    /// pg_names: list<item: string>
    /// gg_accessions: list<item: string>
    /// gg_names: list<item: string>
    /// gg_qvalue: double
    /// anchor_protein: string not null
    /// grouped_runs: list<item: string> not null
    /// global_qvalue: double
    /// pg_qvalue: double
    /// label: string
    /// intensity: float
    /// additional_intensities: list<item: struct<label: string not null, intensities: list<item: struct<intensity_name: string not null, intensity_value: float not null>> not null>>
    /// is_decoy: bool not null
    /// contaminant: bool
    /// peptides: list<item: struct<protein_name: string not null, peptide_count: int32 not null>> not null
    /// peptide_counts: struct<unique_sequences: int32 not null, total_sequences: int32 not null>
    /// feature_counts: struct<unique_features: int32 not null, total_features: int32 not null>
    /// sequence_coverage: float
    /// molecular_weight: float
    /// additional_scores: list<item: struct<score_name: string not null, score_value: double not null, higher_better: bool>>
    /// cv_params: list<item: struct<cv_name: string not null, cv_value: string not null>>
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::PG_ID, DataType::Int64, false),
            field(Self::PG_ACCESSIONS, list(DataType::Utf8), false),
            field(Self::PG_NAMES, list(DataType::Utf8), true),
            field(Self::GG_ACCESSIONS, list(DataType::Utf8), true),
            field(Self::GG_NAMES, list(DataType::Utf8), true),
            field(Self::GG_QVALUE, DataType::Float64, true),
            field(Self::ANCHOR_PROTEIN, DataType::Utf8, false),
            field(Self::GROUPED_RUNS, Self::grouped_runs_type(), false),
            field(Self::GLOBAL_QVALUE, DataType::Float64, true),
            field(Self::PG_QVALUE, DataType::Float64, true),
            field(Self::LABEL, DataType::Utf8, true),
            field(Self::INTENSITY, DataType::Float32, true),
            field(
                Self::ADDITIONAL_INTENSITIES,
                Self::additional_intensities_type(),
                true,
            ),
            field(Self::IS_DECOY, DataType::Boolean, false),
            field(Self::CONTAMINANT, DataType::Boolean, true),
            field(Self::PEPTIDES, Self::peptides_type(), false),
            field(Self::PEPTIDE_COUNTS, Self::peptide_counts_type(), true),
            field(Self::FEATURE_COUNTS, Self::feature_counts_type(), true),
            field(Self::SEQUENCE_COVERAGE, DataType::Float32, true),
            field(Self::MOLECULAR_WEIGHT, DataType::Float32, true),
            field(
                Self::ADDITIONAL_SCORES,
                Self::additional_scores_type(),
                true,
            ),
            field(Self::CV_PARAMS, Self::cv_params_type(), true),
        ])
    }
}

impl SpectraLongSchema {
    /// Column name `mz`.
    pub const MZ: &'static str = "mz";
    /// Column name `intensity`.
    pub const INTENSITY: &'static str = "intensity";
    /// Column name `rt`.
    pub const RT: &'static str = "rt";
    /// Column name `ion_mobility`.
    pub const ION_MOBILITY: &'static str = "ion_mobility";
    /// Column name `spectrum_index`.
    pub const SPECTRUM_INDEX: &'static str = "spectrum_index";
    /// Column name `ms_level`.
    pub const MS_LEVEL: &'static str = "ms_level";
    /// Column name `native_id`.
    pub const NATIVE_ID: &'static str = "native_id";
    /// Column name `precursor_mz`.
    pub const PRECURSOR_MZ: &'static str = "precursor_mz";
    /// Column name `precursor_charge`.
    pub const PRECURSOR_CHARGE: &'static str = "precursor_charge";
    /// Column name `precursor_intensity`.
    pub const PRECURSOR_INTENSITY: &'static str = "precursor_intensity";
    /// Column name `isolation_lower`.
    pub const ISOLATION_LOWER: &'static str = "isolation_lower";
    /// Column name `isolation_upper`.
    pub const ISOLATION_UPPER: &'static str = "isolation_upper";
    /// The complete Arrow schema, 12 fields.
    ///
    /// Source `SpectraLongSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// mz: double
    /// intensity: float
    /// rt: float
    /// ion_mobility: float
    /// spectrum_index: uint32
    /// ms_level: uint8
    /// native_id: string
    /// precursor_mz: double
    /// precursor_charge: int16
    /// precursor_intensity: float
    /// isolation_lower: double
    /// isolation_upper: double
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::MZ, DataType::Float64, true),
            field(Self::INTENSITY, DataType::Float32, true),
            field(Self::RT, DataType::Float32, true),
            field(Self::ION_MOBILITY, DataType::Float32, true),
            field(Self::SPECTRUM_INDEX, DataType::UInt32, true),
            field(Self::MS_LEVEL, DataType::UInt8, true),
            field(Self::NATIVE_ID, DataType::Utf8, true),
            field(Self::PRECURSOR_MZ, DataType::Float64, true),
            field(Self::PRECURSOR_CHARGE, DataType::Int16, true),
            field(Self::PRECURSOR_INTENSITY, DataType::Float32, true),
            field(Self::ISOLATION_LOWER, DataType::Float64, true),
            field(Self::ISOLATION_UPPER, DataType::Float64, true),
        ])
    }
}

impl SpectraSemiWideSchema {
    /// Column name `spectrum_index`.
    pub const SPECTRUM_INDEX: &'static str = "spectrum_index";
    /// Column name `rt`.
    pub const RT: &'static str = "rt";
    /// Column name `ms_level`.
    pub const MS_LEVEL: &'static str = "ms_level";
    /// Column name `native_id`.
    pub const NATIVE_ID: &'static str = "native_id";
    /// Column name `mz`.
    pub const MZ: &'static str = "mz";
    /// Column name `intensity`.
    pub const INTENSITY: &'static str = "intensity";
    /// Column name `ion_mobility`.
    pub const ION_MOBILITY: &'static str = "ion_mobility";
    /// Column name `precursor_mz`.
    pub const PRECURSOR_MZ: &'static str = "precursor_mz";
    /// Column name `precursor_charge`.
    pub const PRECURSOR_CHARGE: &'static str = "precursor_charge";
    /// Column name `precursor_intensity`.
    pub const PRECURSOR_INTENSITY: &'static str = "precursor_intensity";
    /// Column name `isolation_lower`.
    pub const ISOLATION_LOWER: &'static str = "isolation_lower";
    /// Column name `isolation_upper`.
    pub const ISOLATION_UPPER: &'static str = "isolation_upper";
    /// The complete Arrow schema, 12 fields.
    ///
    /// Source `SpectraSemiWideSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// spectrum_index: uint32
    /// rt: float
    /// ms_level: uint8
    /// native_id: string
    /// mz: list<item: double>
    /// intensity: list<item: float>
    /// ion_mobility: list<item: float>
    /// precursor_mz: double
    /// precursor_charge: int16
    /// precursor_intensity: float
    /// isolation_lower: double
    /// isolation_upper: double
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::SPECTRUM_INDEX, DataType::UInt32, true),
            field(Self::RT, DataType::Float32, true),
            field(Self::MS_LEVEL, DataType::UInt8, true),
            field(Self::NATIVE_ID, DataType::Utf8, true),
            field(Self::MZ, list(DataType::Float64), true),
            field(Self::INTENSITY, list(DataType::Float32), true),
            field(Self::ION_MOBILITY, list(DataType::Float32), true),
            field(Self::PRECURSOR_MZ, DataType::Float64, true),
            field(Self::PRECURSOR_CHARGE, DataType::Int16, true),
            field(Self::PRECURSOR_INTENSITY, DataType::Float32, true),
            field(Self::ISOLATION_LOWER, DataType::Float64, true),
            field(Self::ISOLATION_UPPER, DataType::Float64, true),
        ])
    }
}

impl ChromatogramSchema {
    /// Column name `rt`.
    pub const RT: &'static str = "rt";
    /// Column name `intensity`.
    pub const INTENSITY: &'static str = "intensity";
    /// Column name `chromatogram_index`.
    pub const CHROMATOGRAM_INDEX: &'static str = "chromatogram_index";
    /// Column name `native_id`.
    pub const NATIVE_ID: &'static str = "native_id";
    /// Column name `precursor_mz`.
    pub const PRECURSOR_MZ: &'static str = "precursor_mz";
    /// Column name `product_mz`.
    pub const PRODUCT_MZ: &'static str = "product_mz";
    /// The complete Arrow schema, 6 fields.
    ///
    /// Source `ChromatogramSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// rt: double
    /// intensity: float
    /// chromatogram_index: uint32
    /// native_id: string
    /// precursor_mz: double
    /// product_mz: double
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::RT, DataType::Float64, true),
            field(Self::INTENSITY, DataType::Float32, true),
            field(Self::CHROMATOGRAM_INDEX, DataType::UInt32, true),
            field(Self::NATIVE_ID, DataType::Utf8, true),
            field(Self::PRECURSOR_MZ, DataType::Float64, true),
            field(Self::PRODUCT_MZ, DataType::Float64, true),
        ])
    }
}

impl ChromatogramSemiWideSchema {
    /// Column name `chromatogram_index`.
    pub const CHROMATOGRAM_INDEX: &'static str = "chromatogram_index";
    /// Column name `native_id`.
    pub const NATIVE_ID: &'static str = "native_id";
    /// Column name `rt`.
    pub const RT: &'static str = "rt";
    /// Column name `intensity`.
    pub const INTENSITY: &'static str = "intensity";
    /// Column name `precursor_mz`.
    pub const PRECURSOR_MZ: &'static str = "precursor_mz";
    /// Column name `product_mz`.
    pub const PRODUCT_MZ: &'static str = "product_mz";
    /// The complete Arrow schema, 6 fields.
    ///
    /// Source `ChromatogramSemiWideSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// chromatogram_index: uint32
    /// native_id: string
    /// rt: list<item: double>
    /// intensity: list<item: float>
    /// precursor_mz: double
    /// product_mz: double
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::CHROMATOGRAM_INDEX, DataType::UInt32, true),
            field(Self::NATIVE_ID, DataType::Utf8, true),
            field(Self::RT, list(DataType::Float64), true),
            field(Self::INTENSITY, list(DataType::Float32), true),
            field(Self::PRECURSOR_MZ, DataType::Float64, true),
            field(Self::PRODUCT_MZ, DataType::Float64, true),
        ])
    }
}

impl OSWPrecursorSchema {
    /// Column name `precursor_id`.
    pub const PRECURSOR_ID: &'static str = "precursor_id";
    /// Column name `precursor_mz`.
    pub const PRECURSOR_MZ: &'static str = "precursor_mz";
    /// Column name `charge`.
    pub const CHARGE: &'static str = "charge";
    /// Column name `library_rt`.
    pub const LIBRARY_RT: &'static str = "library_rt";
    /// Column name `library_drift_time`.
    pub const LIBRARY_DRIFT_TIME: &'static str = "library_drift_time";
    /// Column name `decoy`.
    pub const DECOY: &'static str = "decoy";
    /// Column name `traml_id`.
    pub const TRAML_ID: &'static str = "traml_id";
    /// Column name `modified_sequence`.
    pub const MODIFIED_SEQUENCE: &'static str = "modified_sequence";
    /// Column name `unmodified_sequence`.
    pub const UNMODIFIED_SEQUENCE: &'static str = "unmodified_sequence";
    /// Column name `protein_accessions`.
    pub const PROTEIN_ACCESSIONS: &'static str = "protein_accessions";
    /// The complete Arrow schema, 10 fields.
    ///
    /// Source `OSWPrecursorSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// precursor_id: int64
    /// precursor_mz: double
    /// charge: int32
    /// library_rt: double
    /// library_drift_time: double
    /// decoy: bool
    /// traml_id: string
    /// modified_sequence: string
    /// unmodified_sequence: string
    /// protein_accessions: string
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::PRECURSOR_ID, DataType::Int64, true),
            field(Self::PRECURSOR_MZ, DataType::Float64, true),
            field(Self::CHARGE, DataType::Int32, true),
            field(Self::LIBRARY_RT, DataType::Float64, true),
            field(Self::LIBRARY_DRIFT_TIME, DataType::Float64, true),
            field(Self::DECOY, DataType::Boolean, true),
            field(Self::TRAML_ID, DataType::Utf8, true),
            field(Self::MODIFIED_SEQUENCE, DataType::Utf8, true),
            field(Self::UNMODIFIED_SEQUENCE, DataType::Utf8, true),
            field(Self::PROTEIN_ACCESSIONS, DataType::Utf8, true),
        ])
    }
}

impl OSWTransitionSchema {
    /// Column name `transition_id`.
    pub const TRANSITION_ID: &'static str = "transition_id";
    /// Column name `precursor_id`.
    pub const PRECURSOR_ID: &'static str = "precursor_id";
    /// Column name `traml_id`.
    pub const TRAML_ID: &'static str = "traml_id";
    /// Column name `product_mz`.
    pub const PRODUCT_MZ: &'static str = "product_mz";
    /// Column name `charge`.
    pub const CHARGE: &'static str = "charge";
    /// Column name `type`.
    pub const TYPE: &'static str = "type";
    /// Column name `annotation`.
    pub const ANNOTATION: &'static str = "annotation";
    /// Column name `ordinal`.
    pub const ORDINAL: &'static str = "ordinal";
    /// Column name `detecting`.
    pub const DETECTING: &'static str = "detecting";
    /// Column name `identifying`.
    pub const IDENTIFYING: &'static str = "identifying";
    /// Column name `quantifying`.
    pub const QUANTIFYING: &'static str = "quantifying";
    /// Column name `library_intensity`.
    pub const LIBRARY_INTENSITY: &'static str = "library_intensity";
    /// Column name `decoy`.
    pub const DECOY: &'static str = "decoy";
    /// The complete Arrow schema, 13 fields.
    ///
    /// Source `OSWTransitionSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// transition_id: int64
    /// precursor_id: int64
    /// traml_id: string
    /// product_mz: double
    /// charge: int32
    /// type: string
    /// annotation: string
    /// ordinal: int32
    /// detecting: bool
    /// identifying: bool
    /// quantifying: bool
    /// library_intensity: double
    /// decoy: bool
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::TRANSITION_ID, DataType::Int64, true),
            field(Self::PRECURSOR_ID, DataType::Int64, true),
            field(Self::TRAML_ID, DataType::Utf8, true),
            field(Self::PRODUCT_MZ, DataType::Float64, true),
            field(Self::CHARGE, DataType::Int32, true),
            field(Self::TYPE, DataType::Utf8, true),
            field(Self::ANNOTATION, DataType::Utf8, true),
            field(Self::ORDINAL, DataType::Int32, true),
            field(Self::DETECTING, DataType::Boolean, true),
            field(Self::IDENTIFYING, DataType::Boolean, true),
            field(Self::QUANTIFYING, DataType::Boolean, true),
            field(Self::LIBRARY_INTENSITY, DataType::Float64, true),
            field(Self::DECOY, DataType::Boolean, true),
        ])
    }
}

impl OSWFeaturePrecursorSchema {
    /// Column name `feature_id`.
    pub const FEATURE_ID: &'static str = "feature_id";
    /// Column name `run_id`.
    pub const RUN_ID: &'static str = "run_id";
    /// Column name `precursor_isotope`.
    pub const PRECURSOR_ISOTOPE: &'static str = "precursor_isotope";
    /// Column name `precursor_area_intensity`.
    pub const PRECURSOR_AREA_INTENSITY: &'static str = "precursor_area_intensity";
    /// Column name `precursor_apex_intensity`.
    pub const PRECURSOR_APEX_INTENSITY: &'static str = "precursor_apex_intensity";
    /// The complete Arrow schema, 5 fields.
    ///
    /// Source `OSWFeaturePrecursorSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// feature_id: int64
    /// run_id: int64
    /// precursor_isotope: int32
    /// precursor_area_intensity: double
    /// precursor_apex_intensity: double
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::FEATURE_ID, DataType::Int64, true),
            field(Self::RUN_ID, DataType::Int64, true),
            field(Self::PRECURSOR_ISOTOPE, DataType::Int32, true),
            field(Self::PRECURSOR_AREA_INTENSITY, DataType::Float64, true),
            field(Self::PRECURSOR_APEX_INTENSITY, DataType::Float64, true),
        ])
    }
}

impl OSWRunSchema {
    /// Column name `run_id`.
    pub const RUN_ID: &'static str = "run_id";
    /// Column name `filename`.
    pub const FILENAME: &'static str = "filename";
    /// The complete Arrow schema, 2 fields.
    ///
    /// Source `OSWRunSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// run_id: int64
    /// filename: string
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::RUN_ID, DataType::Int64, true),
            field(Self::FILENAME, DataType::Utf8, true),
        ])
    }
}

impl OSWFeatureSchema {
    /// Column name `feature_id`.
    pub const FEATURE_ID: &'static str = "feature_id";
    /// Column name `run_id`.
    pub const RUN_ID: &'static str = "run_id";
    /// Column name `precursor_id`.
    pub const PRECURSOR_ID: &'static str = "precursor_id";
    /// Column name `exp_rt`.
    pub const EXP_RT: &'static str = "exp_rt";
    /// Column name `exp_im`.
    pub const EXP_IM: &'static str = "exp_im";
    /// Column name `norm_rt`.
    pub const NORM_RT: &'static str = "norm_rt";
    /// Column name `delta_rt`.
    pub const DELTA_RT: &'static str = "delta_rt";
    /// Column name `left_width`.
    pub const LEFT_WIDTH: &'static str = "left_width";
    /// Column name `right_width`.
    pub const RIGHT_WIDTH: &'static str = "right_width";
    /// Column name `exp_im_leftwidth`.
    pub const EXP_IM_LEFTWIDTH: &'static str = "exp_im_leftwidth";
    /// Column name `exp_im_rightwidth`.
    pub const EXP_IM_RIGHTWIDTH: &'static str = "exp_im_rightwidth";
    /// Column name `ms1_area_intensity`.
    pub const MS1_AREA_INTENSITY: &'static str = "ms1_area_intensity";
    /// Column name `ms1_apex_intensity`.
    pub const MS1_APEX_INTENSITY: &'static str = "ms1_apex_intensity";
    /// Column name `ms1_exp_im`.
    pub const MS1_EXP_IM: &'static str = "ms1_exp_im";
    /// Column name `ms1_delta_im`.
    pub const MS1_DELTA_IM: &'static str = "ms1_delta_im";
    /// Column name `var_ms1_massdev_score`.
    pub const VAR_MS1_MASSDEV_SCORE: &'static str = "var_ms1_massdev_score";
    /// Column name `var_ms1_im_ms1_delta_score`.
    pub const VAR_MS1_IM_MS1_DELTA_SCORE: &'static str = "var_ms1_im_ms1_delta_score";
    /// Column name `var_ms1_mi_score`.
    pub const VAR_MS1_MI_SCORE: &'static str = "var_ms1_mi_score";
    /// Column name `var_ms1_mi_contrast_score`.
    pub const VAR_MS1_MI_CONTRAST_SCORE: &'static str = "var_ms1_mi_contrast_score";
    /// Column name `var_ms1_mi_combined_score`.
    pub const VAR_MS1_MI_COMBINED_SCORE: &'static str = "var_ms1_mi_combined_score";
    /// Column name `var_ms1_isotope_correlation_score`.
    pub const VAR_MS1_ISOTOPE_CORRELATION_SCORE: &'static str = "var_ms1_isotope_correlation_score";
    /// Column name `var_ms1_isotope_overlap_score`.
    pub const VAR_MS1_ISOTOPE_OVERLAP_SCORE: &'static str = "var_ms1_isotope_overlap_score";
    /// Column name `var_ms1_xcorr_coelution`.
    pub const VAR_MS1_XCORR_COELUTION: &'static str = "var_ms1_xcorr_coelution";
    /// Column name `var_ms1_xcorr_coelution_contrast`.
    pub const VAR_MS1_XCORR_COELUTION_CONTRAST: &'static str = "var_ms1_xcorr_coelution_contrast";
    /// Column name `var_ms1_xcorr_coelution_combined`.
    pub const VAR_MS1_XCORR_COELUTION_COMBINED: &'static str = "var_ms1_xcorr_coelution_combined";
    /// Column name `var_ms1_xcorr_shape`.
    pub const VAR_MS1_XCORR_SHAPE: &'static str = "var_ms1_xcorr_shape";
    /// Column name `var_ms1_xcorr_shape_contrast`.
    pub const VAR_MS1_XCORR_SHAPE_CONTRAST: &'static str = "var_ms1_xcorr_shape_contrast";
    /// Column name `var_ms1_xcorr_shape_combined`.
    pub const VAR_MS1_XCORR_SHAPE_COMBINED: &'static str = "var_ms1_xcorr_shape_combined";
    /// Column name `ms2_area_intensity`.
    pub const MS2_AREA_INTENSITY: &'static str = "ms2_area_intensity";
    /// Column name `ms2_total_area_intensity`.
    pub const MS2_TOTAL_AREA_INTENSITY: &'static str = "ms2_total_area_intensity";
    /// Column name `ms2_apex_intensity`.
    pub const MS2_APEX_INTENSITY: &'static str = "ms2_apex_intensity";
    /// Column name `ms2_exp_im`.
    pub const MS2_EXP_IM: &'static str = "ms2_exp_im";
    /// Column name `ms2_exp_im_leftwidth`.
    pub const MS2_EXP_IM_LEFTWIDTH: &'static str = "ms2_exp_im_leftwidth";
    /// Column name `ms2_exp_im_rightwidth`.
    pub const MS2_EXP_IM_RIGHTWIDTH: &'static str = "ms2_exp_im_rightwidth";
    /// Column name `ms2_delta_im`.
    pub const MS2_DELTA_IM: &'static str = "ms2_delta_im";
    /// Column name `ms2_total_mi`.
    pub const MS2_TOTAL_MI: &'static str = "ms2_total_mi";
    /// Column name `var_ms2_bseries_score`.
    pub const VAR_MS2_BSERIES_SCORE: &'static str = "var_ms2_bseries_score";
    /// Column name `var_ms2_dotprod_score`.
    pub const VAR_MS2_DOTPROD_SCORE: &'static str = "var_ms2_dotprod_score";
    /// Column name `var_ms2_intensity_score`.
    pub const VAR_MS2_INTENSITY_SCORE: &'static str = "var_ms2_intensity_score";
    /// Column name `var_ms2_isotope_correlation_score`.
    pub const VAR_MS2_ISOTOPE_CORRELATION_SCORE: &'static str = "var_ms2_isotope_correlation_score";
    /// Column name `var_ms2_isotope_overlap_score`.
    pub const VAR_MS2_ISOTOPE_OVERLAP_SCORE: &'static str = "var_ms2_isotope_overlap_score";
    /// Column name `var_ms2_library_corr`.
    pub const VAR_MS2_LIBRARY_CORR: &'static str = "var_ms2_library_corr";
    /// Column name `var_ms2_library_dotprod`.
    pub const VAR_MS2_LIBRARY_DOTPROD: &'static str = "var_ms2_library_dotprod";
    /// Column name `var_ms2_library_manhattan`.
    pub const VAR_MS2_LIBRARY_MANHATTAN: &'static str = "var_ms2_library_manhattan";
    /// Column name `var_ms2_library_rmsd`.
    pub const VAR_MS2_LIBRARY_RMSD: &'static str = "var_ms2_library_rmsd";
    /// Column name `var_ms2_library_rootmeansquare`.
    pub const VAR_MS2_LIBRARY_ROOTMEANSQUARE: &'static str = "var_ms2_library_rootmeansquare";
    /// Column name `var_ms2_library_sangle`.
    pub const VAR_MS2_LIBRARY_SANGLE: &'static str = "var_ms2_library_sangle";
    /// Column name `var_ms2_log_sn_score`.
    pub const VAR_MS2_LOG_SN_SCORE: &'static str = "var_ms2_log_sn_score";
    /// Column name `var_ms2_manhattan_score`.
    pub const VAR_MS2_MANHATTAN_SCORE: &'static str = "var_ms2_manhattan_score";
    /// Column name `var_ms2_massdev_score`.
    pub const VAR_MS2_MASSDEV_SCORE: &'static str = "var_ms2_massdev_score";
    /// Column name `var_ms2_massdev_score_weighted`.
    pub const VAR_MS2_MASSDEV_SCORE_WEIGHTED: &'static str = "var_ms2_massdev_score_weighted";
    /// Column name `var_ms2_mi_score`.
    pub const VAR_MS2_MI_SCORE: &'static str = "var_ms2_mi_score";
    /// Column name `var_ms2_mi_weighted_score`.
    pub const VAR_MS2_MI_WEIGHTED_SCORE: &'static str = "var_ms2_mi_weighted_score";
    /// Column name `var_ms2_mi_ratio_score`.
    pub const VAR_MS2_MI_RATIO_SCORE: &'static str = "var_ms2_mi_ratio_score";
    /// Column name `var_ms2_norm_rt_score`.
    pub const VAR_MS2_NORM_RT_SCORE: &'static str = "var_ms2_norm_rt_score";
    /// Column name `var_ms2_xcorr_coelution`.
    pub const VAR_MS2_XCORR_COELUTION: &'static str = "var_ms2_xcorr_coelution";
    /// Column name `var_ms2_xcorr_coelution_weighted`.
    pub const VAR_MS2_XCORR_COELUTION_WEIGHTED: &'static str = "var_ms2_xcorr_coelution_weighted";
    /// Column name `var_ms2_xcorr_shape`.
    pub const VAR_MS2_XCORR_SHAPE: &'static str = "var_ms2_xcorr_shape";
    /// Column name `var_ms2_xcorr_shape_weighted`.
    pub const VAR_MS2_XCORR_SHAPE_WEIGHTED: &'static str = "var_ms2_xcorr_shape_weighted";
    /// Column name `var_ms2_yseries_score`.
    pub const VAR_MS2_YSERIES_SCORE: &'static str = "var_ms2_yseries_score";
    /// Column name `var_ms2_elution_model_fit_score`.
    pub const VAR_MS2_ELUTION_MODEL_FIT_SCORE: &'static str = "var_ms2_elution_model_fit_score";
    /// Column name `var_ms2_im_xcorr_shape`.
    pub const VAR_MS2_IM_XCORR_SHAPE: &'static str = "var_ms2_im_xcorr_shape";
    /// Column name `var_ms2_im_xcorr_coelution`.
    pub const VAR_MS2_IM_XCORR_COELUTION: &'static str = "var_ms2_im_xcorr_coelution";
    /// Column name `var_ms2_im_delta_score`.
    pub const VAR_MS2_IM_DELTA_SCORE: &'static str = "var_ms2_im_delta_score";
    /// Column name `var_ms2_im_log_intensity`.
    pub const VAR_MS2_IM_LOG_INTENSITY: &'static str = "var_ms2_im_log_intensity";
    /// The complete Arrow schema, 65 fields.
    ///
    /// Source `OSWFeatureSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// feature_id: int64
    /// run_id: int64
    /// precursor_id: int64
    /// exp_rt: double
    /// exp_im: double
    /// norm_rt: double
    /// delta_rt: double
    /// left_width: double
    /// right_width: double
    /// exp_im_leftwidth: double
    /// exp_im_rightwidth: double
    /// ms1_area_intensity: double
    /// ms1_apex_intensity: double
    /// ms1_exp_im: double
    /// ms1_delta_im: double
    /// var_ms1_massdev_score: double
    /// var_ms1_im_ms1_delta_score: double
    /// var_ms1_mi_score: double
    /// var_ms1_mi_contrast_score: double
    /// var_ms1_mi_combined_score: double
    /// var_ms1_isotope_correlation_score: double
    /// var_ms1_isotope_overlap_score: double
    /// var_ms1_xcorr_coelution: double
    /// var_ms1_xcorr_coelution_contrast: double
    /// var_ms1_xcorr_coelution_combined: double
    /// var_ms1_xcorr_shape: double
    /// var_ms1_xcorr_shape_contrast: double
    /// var_ms1_xcorr_shape_combined: double
    /// ms2_area_intensity: double
    /// ms2_total_area_intensity: double
    /// ms2_apex_intensity: double
    /// ms2_exp_im: double
    /// ms2_exp_im_leftwidth: double
    /// ms2_exp_im_rightwidth: double
    /// ms2_delta_im: double
    /// ms2_total_mi: double
    /// var_ms2_bseries_score: double
    /// var_ms2_dotprod_score: double
    /// var_ms2_intensity_score: double
    /// var_ms2_isotope_correlation_score: double
    /// var_ms2_isotope_overlap_score: double
    /// var_ms2_library_corr: double
    /// var_ms2_library_dotprod: double
    /// var_ms2_library_manhattan: double
    /// var_ms2_library_rmsd: double
    /// var_ms2_library_rootmeansquare: double
    /// var_ms2_library_sangle: double
    /// var_ms2_log_sn_score: double
    /// var_ms2_manhattan_score: double
    /// var_ms2_massdev_score: double
    /// var_ms2_massdev_score_weighted: double
    /// var_ms2_mi_score: double
    /// var_ms2_mi_weighted_score: double
    /// var_ms2_mi_ratio_score: double
    /// var_ms2_norm_rt_score: double
    /// var_ms2_xcorr_coelution: double
    /// var_ms2_xcorr_coelution_weighted: double
    /// var_ms2_xcorr_shape: double
    /// var_ms2_xcorr_shape_weighted: double
    /// var_ms2_yseries_score: double
    /// var_ms2_elution_model_fit_score: double
    /// var_ms2_im_xcorr_shape: double
    /// var_ms2_im_xcorr_coelution: double
    /// var_ms2_im_delta_score: double
    /// var_ms2_im_log_intensity: double
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::FEATURE_ID, DataType::Int64, true),
            field(Self::RUN_ID, DataType::Int64, true),
            field(Self::PRECURSOR_ID, DataType::Int64, true),
            field(Self::EXP_RT, DataType::Float64, true),
            field(Self::EXP_IM, DataType::Float64, true),
            field(Self::NORM_RT, DataType::Float64, true),
            field(Self::DELTA_RT, DataType::Float64, true),
            field(Self::LEFT_WIDTH, DataType::Float64, true),
            field(Self::RIGHT_WIDTH, DataType::Float64, true),
            field(Self::EXP_IM_LEFTWIDTH, DataType::Float64, true),
            field(Self::EXP_IM_RIGHTWIDTH, DataType::Float64, true),
            field(Self::MS1_AREA_INTENSITY, DataType::Float64, true),
            field(Self::MS1_APEX_INTENSITY, DataType::Float64, true),
            field(Self::MS1_EXP_IM, DataType::Float64, true),
            field(Self::MS1_DELTA_IM, DataType::Float64, true),
            field(Self::VAR_MS1_MASSDEV_SCORE, DataType::Float64, true),
            field(Self::VAR_MS1_IM_MS1_DELTA_SCORE, DataType::Float64, true),
            field(Self::VAR_MS1_MI_SCORE, DataType::Float64, true),
            field(Self::VAR_MS1_MI_CONTRAST_SCORE, DataType::Float64, true),
            field(Self::VAR_MS1_MI_COMBINED_SCORE, DataType::Float64, true),
            field(
                Self::VAR_MS1_ISOTOPE_CORRELATION_SCORE,
                DataType::Float64,
                true,
            ),
            field(Self::VAR_MS1_ISOTOPE_OVERLAP_SCORE, DataType::Float64, true),
            field(Self::VAR_MS1_XCORR_COELUTION, DataType::Float64, true),
            field(
                Self::VAR_MS1_XCORR_COELUTION_CONTRAST,
                DataType::Float64,
                true,
            ),
            field(
                Self::VAR_MS1_XCORR_COELUTION_COMBINED,
                DataType::Float64,
                true,
            ),
            field(Self::VAR_MS1_XCORR_SHAPE, DataType::Float64, true),
            field(Self::VAR_MS1_XCORR_SHAPE_CONTRAST, DataType::Float64, true),
            field(Self::VAR_MS1_XCORR_SHAPE_COMBINED, DataType::Float64, true),
            field(Self::MS2_AREA_INTENSITY, DataType::Float64, true),
            field(Self::MS2_TOTAL_AREA_INTENSITY, DataType::Float64, true),
            field(Self::MS2_APEX_INTENSITY, DataType::Float64, true),
            field(Self::MS2_EXP_IM, DataType::Float64, true),
            field(Self::MS2_EXP_IM_LEFTWIDTH, DataType::Float64, true),
            field(Self::MS2_EXP_IM_RIGHTWIDTH, DataType::Float64, true),
            field(Self::MS2_DELTA_IM, DataType::Float64, true),
            field(Self::MS2_TOTAL_MI, DataType::Float64, true),
            field(Self::VAR_MS2_BSERIES_SCORE, DataType::Float64, true),
            field(Self::VAR_MS2_DOTPROD_SCORE, DataType::Float64, true),
            field(Self::VAR_MS2_INTENSITY_SCORE, DataType::Float64, true),
            field(
                Self::VAR_MS2_ISOTOPE_CORRELATION_SCORE,
                DataType::Float64,
                true,
            ),
            field(Self::VAR_MS2_ISOTOPE_OVERLAP_SCORE, DataType::Float64, true),
            field(Self::VAR_MS2_LIBRARY_CORR, DataType::Float64, true),
            field(Self::VAR_MS2_LIBRARY_DOTPROD, DataType::Float64, true),
            field(Self::VAR_MS2_LIBRARY_MANHATTAN, DataType::Float64, true),
            field(Self::VAR_MS2_LIBRARY_RMSD, DataType::Float64, true),
            field(
                Self::VAR_MS2_LIBRARY_ROOTMEANSQUARE,
                DataType::Float64,
                true,
            ),
            field(Self::VAR_MS2_LIBRARY_SANGLE, DataType::Float64, true),
            field(Self::VAR_MS2_LOG_SN_SCORE, DataType::Float64, true),
            field(Self::VAR_MS2_MANHATTAN_SCORE, DataType::Float64, true),
            field(Self::VAR_MS2_MASSDEV_SCORE, DataType::Float64, true),
            field(
                Self::VAR_MS2_MASSDEV_SCORE_WEIGHTED,
                DataType::Float64,
                true,
            ),
            field(Self::VAR_MS2_MI_SCORE, DataType::Float64, true),
            field(Self::VAR_MS2_MI_WEIGHTED_SCORE, DataType::Float64, true),
            field(Self::VAR_MS2_MI_RATIO_SCORE, DataType::Float64, true),
            field(Self::VAR_MS2_NORM_RT_SCORE, DataType::Float64, true),
            field(Self::VAR_MS2_XCORR_COELUTION, DataType::Float64, true),
            field(
                Self::VAR_MS2_XCORR_COELUTION_WEIGHTED,
                DataType::Float64,
                true,
            ),
            field(Self::VAR_MS2_XCORR_SHAPE, DataType::Float64, true),
            field(Self::VAR_MS2_XCORR_SHAPE_WEIGHTED, DataType::Float64, true),
            field(Self::VAR_MS2_YSERIES_SCORE, DataType::Float64, true),
            field(
                Self::VAR_MS2_ELUTION_MODEL_FIT_SCORE,
                DataType::Float64,
                true,
            ),
            field(Self::VAR_MS2_IM_XCORR_SHAPE, DataType::Float64, true),
            field(Self::VAR_MS2_IM_XCORR_COELUTION, DataType::Float64, true),
            field(Self::VAR_MS2_IM_DELTA_SCORE, DataType::Float64, true),
            field(Self::VAR_MS2_IM_LOG_INTENSITY, DataType::Float64, true),
        ])
    }
}

impl OSWFeatureTransitionSchema {
    /// Column name `feature_id`.
    pub const FEATURE_ID: &'static str = "feature_id";
    /// Column name `run_id`.
    pub const RUN_ID: &'static str = "run_id";
    /// Column name `transition_id`.
    pub const TRANSITION_ID: &'static str = "transition_id";
    /// Column name `area_intensity`.
    pub const AREA_INTENSITY: &'static str = "area_intensity";
    /// Column name `total_area_intensity`.
    pub const TOTAL_AREA_INTENSITY: &'static str = "total_area_intensity";
    /// Column name `apex_intensity`.
    pub const APEX_INTENSITY: &'static str = "apex_intensity";
    /// Column name `apex_rt`.
    pub const APEX_RT: &'static str = "apex_rt";
    /// Column name `rt_fwhm`.
    pub const RT_FWHM: &'static str = "rt_fwhm";
    /// Column name `masserror_ppm`.
    pub const MASSERROR_PPM: &'static str = "masserror_ppm";
    /// Column name `total_mi`.
    pub const TOTAL_MI: &'static str = "total_mi";
    /// Column name `var_intensity_score`.
    pub const VAR_INTENSITY_SCORE: &'static str = "var_intensity_score";
    /// Column name `var_intensity_ratio_score`.
    pub const VAR_INTENSITY_RATIO_SCORE: &'static str = "var_intensity_ratio_score";
    /// Column name `var_log_intensity`.
    pub const VAR_LOG_INTENSITY: &'static str = "var_log_intensity";
    /// Column name `var_xcorr_coelution`.
    pub const VAR_XCORR_COELUTION: &'static str = "var_xcorr_coelution";
    /// Column name `var_xcorr_shape`.
    pub const VAR_XCORR_SHAPE: &'static str = "var_xcorr_shape";
    /// Column name `var_log_sn_score`.
    pub const VAR_LOG_SN_SCORE: &'static str = "var_log_sn_score";
    /// Column name `var_massdev_score`.
    pub const VAR_MASSDEV_SCORE: &'static str = "var_massdev_score";
    /// Column name `var_mi_score`.
    pub const VAR_MI_SCORE: &'static str = "var_mi_score";
    /// Column name `var_mi_ratio_score`.
    pub const VAR_MI_RATIO_SCORE: &'static str = "var_mi_ratio_score";
    /// Column name `var_isotope_correlation_score`.
    pub const VAR_ISOTOPE_CORRELATION_SCORE: &'static str = "var_isotope_correlation_score";
    /// Column name `var_isotope_overlap_score`.
    pub const VAR_ISOTOPE_OVERLAP_SCORE: &'static str = "var_isotope_overlap_score";
    /// Column name `exp_im`.
    pub const EXP_IM: &'static str = "exp_im";
    /// Column name `exp_im_leftwidth`.
    pub const EXP_IM_LEFTWIDTH: &'static str = "exp_im_leftwidth";
    /// Column name `exp_im_rightwidth`.
    pub const EXP_IM_RIGHTWIDTH: &'static str = "exp_im_rightwidth";
    /// Column name `delta_im`.
    pub const DELTA_IM: &'static str = "delta_im";
    /// Column name `var_im_delta_score`.
    pub const VAR_IM_DELTA_SCORE: &'static str = "var_im_delta_score";
    /// Column name `var_im_log_intensity`.
    pub const VAR_IM_LOG_INTENSITY: &'static str = "var_im_log_intensity";
    /// Column name `var_im_xcorr_coelution_contrast`.
    pub const VAR_IM_XCORR_COELUTION_CONTRAST: &'static str = "var_im_xcorr_coelution_contrast";
    /// Column name `var_im_xcorr_shape_contrast`.
    pub const VAR_IM_XCORR_SHAPE_CONTRAST: &'static str = "var_im_xcorr_shape_contrast";
    /// Column name `var_im_xcorr_coelution_combined`.
    pub const VAR_IM_XCORR_COELUTION_COMBINED: &'static str = "var_im_xcorr_coelution_combined";
    /// Column name `var_im_xcorr_shape_combined`.
    pub const VAR_IM_XCORR_SHAPE_COMBINED: &'static str = "var_im_xcorr_shape_combined";
    /// Column name `start_position_at_5`.
    pub const START_POSITION_AT_5: &'static str = "start_position_at_5";
    /// Column name `end_position_at_5`.
    pub const END_POSITION_AT_5: &'static str = "end_position_at_5";
    /// Column name `start_position_at_10`.
    pub const START_POSITION_AT_10: &'static str = "start_position_at_10";
    /// Column name `end_position_at_10`.
    pub const END_POSITION_AT_10: &'static str = "end_position_at_10";
    /// Column name `start_position_at_50`.
    pub const START_POSITION_AT_50: &'static str = "start_position_at_50";
    /// Column name `end_position_at_50`.
    pub const END_POSITION_AT_50: &'static str = "end_position_at_50";
    /// Column name `total_width`.
    pub const TOTAL_WIDTH: &'static str = "total_width";
    /// Column name `tailing_factor`.
    pub const TAILING_FACTOR: &'static str = "tailing_factor";
    /// Column name `asymmetry_factor`.
    pub const ASYMMETRY_FACTOR: &'static str = "asymmetry_factor";
    /// Column name `slope_of_baseline`.
    pub const SLOPE_OF_BASELINE: &'static str = "slope_of_baseline";
    /// Column name `baseline_delta_2_height`.
    pub const BASELINE_DELTA_2_HEIGHT: &'static str = "baseline_delta_2_height";
    /// Column name `points_across_baseline`.
    pub const POINTS_ACROSS_BASELINE: &'static str = "points_across_baseline";
    /// Column name `points_across_half_height`.
    pub const POINTS_ACROSS_HALF_HEIGHT: &'static str = "points_across_half_height";
    /// The complete Arrow schema, 44 fields.
    ///
    /// Source `OSWFeatureTransitionSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// feature_id: int64
    /// run_id: int64
    /// transition_id: int64
    /// area_intensity: double
    /// total_area_intensity: double
    /// apex_intensity: double
    /// apex_rt: double
    /// rt_fwhm: double
    /// masserror_ppm: double
    /// total_mi: double
    /// var_intensity_score: double
    /// var_intensity_ratio_score: double
    /// var_log_intensity: double
    /// var_xcorr_coelution: double
    /// var_xcorr_shape: double
    /// var_log_sn_score: double
    /// var_massdev_score: double
    /// var_mi_score: double
    /// var_mi_ratio_score: double
    /// var_isotope_correlation_score: double
    /// var_isotope_overlap_score: double
    /// exp_im: double
    /// exp_im_leftwidth: double
    /// exp_im_rightwidth: double
    /// delta_im: double
    /// var_im_delta_score: double
    /// var_im_log_intensity: double
    /// var_im_xcorr_coelution_contrast: double
    /// var_im_xcorr_shape_contrast: double
    /// var_im_xcorr_coelution_combined: double
    /// var_im_xcorr_shape_combined: double
    /// start_position_at_5: double
    /// end_position_at_5: double
    /// start_position_at_10: double
    /// end_position_at_10: double
    /// start_position_at_50: double
    /// end_position_at_50: double
    /// total_width: double
    /// tailing_factor: double
    /// asymmetry_factor: double
    /// slope_of_baseline: double
    /// baseline_delta_2_height: double
    /// points_across_baseline: double
    /// points_across_half_height: double
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::FEATURE_ID, DataType::Int64, true),
            field(Self::RUN_ID, DataType::Int64, true),
            field(Self::TRANSITION_ID, DataType::Int64, true),
            field(Self::AREA_INTENSITY, DataType::Float64, true),
            field(Self::TOTAL_AREA_INTENSITY, DataType::Float64, true),
            field(Self::APEX_INTENSITY, DataType::Float64, true),
            field(Self::APEX_RT, DataType::Float64, true),
            field(Self::RT_FWHM, DataType::Float64, true),
            field(Self::MASSERROR_PPM, DataType::Float64, true),
            field(Self::TOTAL_MI, DataType::Float64, true),
            field(Self::VAR_INTENSITY_SCORE, DataType::Float64, true),
            field(Self::VAR_INTENSITY_RATIO_SCORE, DataType::Float64, true),
            field(Self::VAR_LOG_INTENSITY, DataType::Float64, true),
            field(Self::VAR_XCORR_COELUTION, DataType::Float64, true),
            field(Self::VAR_XCORR_SHAPE, DataType::Float64, true),
            field(Self::VAR_LOG_SN_SCORE, DataType::Float64, true),
            field(Self::VAR_MASSDEV_SCORE, DataType::Float64, true),
            field(Self::VAR_MI_SCORE, DataType::Float64, true),
            field(Self::VAR_MI_RATIO_SCORE, DataType::Float64, true),
            field(Self::VAR_ISOTOPE_CORRELATION_SCORE, DataType::Float64, true),
            field(Self::VAR_ISOTOPE_OVERLAP_SCORE, DataType::Float64, true),
            field(Self::EXP_IM, DataType::Float64, true),
            field(Self::EXP_IM_LEFTWIDTH, DataType::Float64, true),
            field(Self::EXP_IM_RIGHTWIDTH, DataType::Float64, true),
            field(Self::DELTA_IM, DataType::Float64, true),
            field(Self::VAR_IM_DELTA_SCORE, DataType::Float64, true),
            field(Self::VAR_IM_LOG_INTENSITY, DataType::Float64, true),
            field(
                Self::VAR_IM_XCORR_COELUTION_CONTRAST,
                DataType::Float64,
                true,
            ),
            field(Self::VAR_IM_XCORR_SHAPE_CONTRAST, DataType::Float64, true),
            field(
                Self::VAR_IM_XCORR_COELUTION_COMBINED,
                DataType::Float64,
                true,
            ),
            field(Self::VAR_IM_XCORR_SHAPE_COMBINED, DataType::Float64, true),
            field(Self::START_POSITION_AT_5, DataType::Float64, true),
            field(Self::END_POSITION_AT_5, DataType::Float64, true),
            field(Self::START_POSITION_AT_10, DataType::Float64, true),
            field(Self::END_POSITION_AT_10, DataType::Float64, true),
            field(Self::START_POSITION_AT_50, DataType::Float64, true),
            field(Self::END_POSITION_AT_50, DataType::Float64, true),
            field(Self::TOTAL_WIDTH, DataType::Float64, true),
            field(Self::TAILING_FACTOR, DataType::Float64, true),
            field(Self::ASYMMETRY_FACTOR, DataType::Float64, true),
            field(Self::SLOPE_OF_BASELINE, DataType::Float64, true),
            field(Self::BASELINE_DELTA_2_HEIGHT, DataType::Float64, true),
            field(Self::POINTS_ACROSS_BASELINE, DataType::Float64, true),
            field(Self::POINTS_ACROSS_HALF_HEIGHT, DataType::Float64, true),
        ])
    }
}

impl XICSchema {
    /// Column name `RUN_ID`.
    pub const RUN_ID: &'static str = "RUN_ID";
    /// Column name `SOURCE_FILE`.
    pub const SOURCE_FILE: &'static str = "SOURCE_FILE";
    /// Column name `MS_LEVEL`.
    pub const MS_LEVEL: &'static str = "MS_LEVEL";
    /// Column name `PRECURSOR_ID`.
    pub const PRECURSOR_ID: &'static str = "PRECURSOR_ID";
    /// Column name `TRANSITION_ID`.
    pub const TRANSITION_ID: &'static str = "TRANSITION_ID";
    /// Column name `MODIFIED_SEQUENCE`.
    pub const MODIFIED_SEQUENCE: &'static str = "MODIFIED_SEQUENCE";
    /// Column name `PRECURSOR_CHARGE`.
    pub const PRECURSOR_CHARGE: &'static str = "PRECURSOR_CHARGE";
    /// Column name `PRODUCT_CHARGE`.
    pub const PRODUCT_CHARGE: &'static str = "PRODUCT_CHARGE";
    /// Column name `DETECTING_TRANSITION`.
    pub const DETECTING_TRANSITION: &'static str = "DETECTING_TRANSITION";
    /// Column name `PRECURSOR_DECOY`.
    pub const PRECURSOR_DECOY: &'static str = "PRECURSOR_DECOY";
    /// Column name `PRODUCT_DECOY`.
    pub const PRODUCT_DECOY: &'static str = "PRODUCT_DECOY";
    /// Column name `TRANSITION_ORDINAL`.
    pub const TRANSITION_ORDINAL: &'static str = "TRANSITION_ORDINAL";
    /// Column name `TRANSITION_TYPE`.
    pub const TRANSITION_TYPE: &'static str = "TRANSITION_TYPE";
    /// Column name `ANNOTATION`.
    pub const ANNOTATION: &'static str = "ANNOTATION";
    /// Column name `RT_DATA`.
    pub const RT_DATA: &'static str = "RT_DATA";
    /// Column name `INTENSITY_DATA`.
    pub const INTENSITY_DATA: &'static str = "INTENSITY_DATA";
    /// Column name `RT_COMPRESSION`.
    pub const RT_COMPRESSION: &'static str = "RT_COMPRESSION";
    /// Column name `INTENSITY_COMPRESSION`.
    pub const INTENSITY_COMPRESSION: &'static str = "INTENSITY_COMPRESSION";
    /// The complete Arrow schema, 18 fields.
    ///
    /// Source `XICSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// RUN_ID: int64
    /// SOURCE_FILE: string
    /// MS_LEVEL: int64
    /// PRECURSOR_ID: int64
    /// TRANSITION_ID: int64
    /// MODIFIED_SEQUENCE: string
    /// PRECURSOR_CHARGE: int64
    /// PRODUCT_CHARGE: int64
    /// DETECTING_TRANSITION: int64
    /// PRECURSOR_DECOY: int64
    /// PRODUCT_DECOY: int64
    /// TRANSITION_ORDINAL: int64
    /// TRANSITION_TYPE: string
    /// ANNOTATION: string
    /// RT_DATA: binary
    /// INTENSITY_DATA: binary
    /// RT_COMPRESSION: int64
    /// INTENSITY_COMPRESSION: int64
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::RUN_ID, DataType::Int64, true),
            field(Self::SOURCE_FILE, DataType::Utf8, true),
            field(Self::MS_LEVEL, DataType::Int64, true),
            field(Self::PRECURSOR_ID, DataType::Int64, true),
            field(Self::TRANSITION_ID, DataType::Int64, true),
            field(Self::MODIFIED_SEQUENCE, DataType::Utf8, true),
            field(Self::PRECURSOR_CHARGE, DataType::Int64, true),
            field(Self::PRODUCT_CHARGE, DataType::Int64, true),
            field(Self::DETECTING_TRANSITION, DataType::Int64, true),
            field(Self::PRECURSOR_DECOY, DataType::Int64, true),
            field(Self::PRODUCT_DECOY, DataType::Int64, true),
            field(Self::TRANSITION_ORDINAL, DataType::Int64, true),
            field(Self::TRANSITION_TYPE, DataType::Utf8, true),
            field(Self::ANNOTATION, DataType::Utf8, true),
            field(Self::RT_DATA, DataType::Binary, true),
            field(Self::INTENSITY_DATA, DataType::Binary, true),
            field(Self::RT_COMPRESSION, DataType::Int64, true),
            field(Self::INTENSITY_COMPRESSION, DataType::Int64, true),
        ])
    }
}

impl XIMSchema {
    /// Column name `RUN_ID`.
    pub const RUN_ID: &'static str = "RUN_ID";
    /// Column name `SOURCE_FILE`.
    pub const SOURCE_FILE: &'static str = "SOURCE_FILE";
    /// Column name `MS_LEVEL`.
    pub const MS_LEVEL: &'static str = "MS_LEVEL";
    /// Column name `MOBILOGRAM_TYPE`.
    pub const MOBILOGRAM_TYPE: &'static str = "MOBILOGRAM_TYPE";
    /// Column name `PRECURSOR_ID`.
    pub const PRECURSOR_ID: &'static str = "PRECURSOR_ID";
    /// Column name `TRANSITION_ID`.
    pub const TRANSITION_ID: &'static str = "TRANSITION_ID";
    /// Column name `FEATURE_ID`.
    pub const FEATURE_ID: &'static str = "FEATURE_ID";
    /// Column name `FEATURE_RT`.
    pub const FEATURE_RT: &'static str = "FEATURE_RT";
    /// Column name `MODIFIED_SEQUENCE`.
    pub const MODIFIED_SEQUENCE: &'static str = "MODIFIED_SEQUENCE";
    /// Column name `PRECURSOR_CHARGE`.
    pub const PRECURSOR_CHARGE: &'static str = "PRECURSOR_CHARGE";
    /// Column name `PRODUCT_CHARGE`.
    pub const PRODUCT_CHARGE: &'static str = "PRODUCT_CHARGE";
    /// Column name `DETECTING_TRANSITION`.
    pub const DETECTING_TRANSITION: &'static str = "DETECTING_TRANSITION";
    /// Column name `PRECURSOR_DECOY`.
    pub const PRECURSOR_DECOY: &'static str = "PRECURSOR_DECOY";
    /// Column name `PRODUCT_DECOY`.
    pub const PRODUCT_DECOY: &'static str = "PRODUCT_DECOY";
    /// Column name `TRANSITION_ORDINAL`.
    pub const TRANSITION_ORDINAL: &'static str = "TRANSITION_ORDINAL";
    /// Column name `TRANSITION_TYPE`.
    pub const TRANSITION_TYPE: &'static str = "TRANSITION_TYPE";
    /// Column name `ANNOTATION`.
    pub const ANNOTATION: &'static str = "ANNOTATION";
    /// Column name `MOBILITY_DATA`.
    pub const MOBILITY_DATA: &'static str = "MOBILITY_DATA";
    /// Column name `INTENSITY_DATA`.
    pub const INTENSITY_DATA: &'static str = "INTENSITY_DATA";
    /// Column name `MOBILITY_COMPRESSION`.
    pub const MOBILITY_COMPRESSION: &'static str = "MOBILITY_COMPRESSION";
    /// Column name `INTENSITY_COMPRESSION`.
    pub const INTENSITY_COMPRESSION: &'static str = "INTENSITY_COMPRESSION";
    /// The complete Arrow schema, 21 fields.
    ///
    /// Source `XIMSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// RUN_ID: int64
    /// SOURCE_FILE: string
    /// MS_LEVEL: int64
    /// MOBILOGRAM_TYPE: string
    /// PRECURSOR_ID: int64
    /// TRANSITION_ID: int64
    /// FEATURE_ID: int64
    /// FEATURE_RT: double
    /// MODIFIED_SEQUENCE: string
    /// PRECURSOR_CHARGE: int64
    /// PRODUCT_CHARGE: int64
    /// DETECTING_TRANSITION: int64
    /// PRECURSOR_DECOY: int64
    /// PRODUCT_DECOY: int64
    /// TRANSITION_ORDINAL: int64
    /// TRANSITION_TYPE: string
    /// ANNOTATION: string
    /// MOBILITY_DATA: binary
    /// INTENSITY_DATA: binary
    /// MOBILITY_COMPRESSION: int64
    /// INTENSITY_COMPRESSION: int64
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::RUN_ID, DataType::Int64, true),
            field(Self::SOURCE_FILE, DataType::Utf8, true),
            field(Self::MS_LEVEL, DataType::Int64, true),
            field(Self::MOBILOGRAM_TYPE, DataType::Utf8, true),
            field(Self::PRECURSOR_ID, DataType::Int64, true),
            field(Self::TRANSITION_ID, DataType::Int64, true),
            field(Self::FEATURE_ID, DataType::Int64, true),
            field(Self::FEATURE_RT, DataType::Float64, true),
            field(Self::MODIFIED_SEQUENCE, DataType::Utf8, true),
            field(Self::PRECURSOR_CHARGE, DataType::Int64, true),
            field(Self::PRODUCT_CHARGE, DataType::Int64, true),
            field(Self::DETECTING_TRANSITION, DataType::Int64, true),
            field(Self::PRECURSOR_DECOY, DataType::Int64, true),
            field(Self::PRODUCT_DECOY, DataType::Int64, true),
            field(Self::TRANSITION_ORDINAL, DataType::Int64, true),
            field(Self::TRANSITION_TYPE, DataType::Utf8, true),
            field(Self::ANNOTATION, DataType::Utf8, true),
            field(Self::MOBILITY_DATA, DataType::Binary, true),
            field(Self::INTENSITY_DATA, DataType::Binary, true),
            field(Self::MOBILITY_COMPRESSION, DataType::Int64, true),
            field(Self::INTENSITY_COMPRESSION, DataType::Int64, true),
        ])
    }
}

impl XIPMSchema {
    /// Column name `RUN_ID`.
    pub const RUN_ID: &'static str = "RUN_ID";
    /// Column name `SOURCE_FILE`.
    pub const SOURCE_FILE: &'static str = "SOURCE_FILE";
    /// Column name `MS_LEVEL`.
    pub const MS_LEVEL: &'static str = "MS_LEVEL";
    /// Column name `PEAKMAP_TYPE`.
    pub const PEAKMAP_TYPE: &'static str = "PEAKMAP_TYPE";
    /// Column name `PRECURSOR_ID`.
    pub const PRECURSOR_ID: &'static str = "PRECURSOR_ID";
    /// Column name `TRANSITION_ID`.
    pub const TRANSITION_ID: &'static str = "TRANSITION_ID";
    /// Column name `MODIFIED_SEQUENCE`.
    pub const MODIFIED_SEQUENCE: &'static str = "MODIFIED_SEQUENCE";
    /// Column name `PRECURSOR_CHARGE`.
    pub const PRECURSOR_CHARGE: &'static str = "PRECURSOR_CHARGE";
    /// Column name `PRODUCT_CHARGE`.
    pub const PRODUCT_CHARGE: &'static str = "PRODUCT_CHARGE";
    /// Column name `DETECTING_TRANSITION`.
    pub const DETECTING_TRANSITION: &'static str = "DETECTING_TRANSITION";
    /// Column name `PRECURSOR_DECOY`.
    pub const PRECURSOR_DECOY: &'static str = "PRECURSOR_DECOY";
    /// Column name `PRODUCT_DECOY`.
    pub const PRODUCT_DECOY: &'static str = "PRODUCT_DECOY";
    /// Column name `TRANSITION_ORDINAL`.
    pub const TRANSITION_ORDINAL: &'static str = "TRANSITION_ORDINAL";
    /// Column name `TRANSITION_TYPE`.
    pub const TRANSITION_TYPE: &'static str = "TRANSITION_TYPE";
    /// Column name `ANNOTATION`.
    pub const ANNOTATION: &'static str = "ANNOTATION";
    /// Column name `TARGET_MZ`.
    pub const TARGET_MZ: &'static str = "TARGET_MZ";
    /// Column name `TARGET_RT`.
    pub const TARGET_RT: &'static str = "TARGET_RT";
    /// Column name `TARGET_ION_MOBILITY`.
    pub const TARGET_ION_MOBILITY: &'static str = "TARGET_ION_MOBILITY";
    /// Column name `RT_START`.
    pub const RT_START: &'static str = "RT_START";
    /// Column name `RT_END`.
    pub const RT_END: &'static str = "RT_END";
    /// Column name `MZ_DATA`.
    pub const MZ_DATA: &'static str = "MZ_DATA";
    /// Column name `RT_DATA`.
    pub const RT_DATA: &'static str = "RT_DATA";
    /// Column name `MOBILITY_DATA`.
    pub const MOBILITY_DATA: &'static str = "MOBILITY_DATA";
    /// Column name `INTENSITY_DATA`.
    pub const INTENSITY_DATA: &'static str = "INTENSITY_DATA";
    /// Column name `MZ_COMPRESSION`.
    pub const MZ_COMPRESSION: &'static str = "MZ_COMPRESSION";
    /// Column name `RT_COMPRESSION`.
    pub const RT_COMPRESSION: &'static str = "RT_COMPRESSION";
    /// Column name `MOBILITY_COMPRESSION`.
    pub const MOBILITY_COMPRESSION: &'static str = "MOBILITY_COMPRESSION";
    /// Column name `INTENSITY_COMPRESSION`.
    pub const INTENSITY_COMPRESSION: &'static str = "INTENSITY_COMPRESSION";
    /// The complete Arrow schema, 28 fields.
    ///
    /// Source `XIPMSchema::schema`. In column order, as C++ `Schema::ToString` prints them:
    ///
    /// ```text
    /// RUN_ID: int64
    /// SOURCE_FILE: string
    /// MS_LEVEL: int64
    /// PEAKMAP_TYPE: string
    /// PRECURSOR_ID: int64
    /// TRANSITION_ID: int64
    /// MODIFIED_SEQUENCE: string
    /// PRECURSOR_CHARGE: int64
    /// PRODUCT_CHARGE: int64
    /// DETECTING_TRANSITION: int64
    /// PRECURSOR_DECOY: int64
    /// PRODUCT_DECOY: int64
    /// TRANSITION_ORDINAL: int64
    /// TRANSITION_TYPE: string
    /// ANNOTATION: string
    /// TARGET_MZ: double
    /// TARGET_RT: double
    /// TARGET_ION_MOBILITY: double
    /// RT_START: double
    /// RT_END: double
    /// MZ_DATA: binary
    /// RT_DATA: binary
    /// MOBILITY_DATA: binary
    /// INTENSITY_DATA: binary
    /// MZ_COMPRESSION: int64
    /// RT_COMPRESSION: int64
    /// MOBILITY_COMPRESSION: int64
    /// INTENSITY_COMPRESSION: int64
    /// ```
    pub fn schema() -> SchemaRef {
        schema(vec![
            field(Self::RUN_ID, DataType::Int64, true),
            field(Self::SOURCE_FILE, DataType::Utf8, true),
            field(Self::MS_LEVEL, DataType::Int64, true),
            field(Self::PEAKMAP_TYPE, DataType::Utf8, true),
            field(Self::PRECURSOR_ID, DataType::Int64, true),
            field(Self::TRANSITION_ID, DataType::Int64, true),
            field(Self::MODIFIED_SEQUENCE, DataType::Utf8, true),
            field(Self::PRECURSOR_CHARGE, DataType::Int64, true),
            field(Self::PRODUCT_CHARGE, DataType::Int64, true),
            field(Self::DETECTING_TRANSITION, DataType::Int64, true),
            field(Self::PRECURSOR_DECOY, DataType::Int64, true),
            field(Self::PRODUCT_DECOY, DataType::Int64, true),
            field(Self::TRANSITION_ORDINAL, DataType::Int64, true),
            field(Self::TRANSITION_TYPE, DataType::Utf8, true),
            field(Self::ANNOTATION, DataType::Utf8, true),
            field(Self::TARGET_MZ, DataType::Float64, true),
            field(Self::TARGET_RT, DataType::Float64, true),
            field(Self::TARGET_ION_MOBILITY, DataType::Float64, true),
            field(Self::RT_START, DataType::Float64, true),
            field(Self::RT_END, DataType::Float64, true),
            field(Self::MZ_DATA, DataType::Binary, true),
            field(Self::RT_DATA, DataType::Binary, true),
            field(Self::MOBILITY_DATA, DataType::Binary, true),
            field(Self::INTENSITY_DATA, DataType::Binary, true),
            field(Self::MZ_COMPRESSION, DataType::Int64, true),
            field(Self::RT_COMPRESSION, DataType::Int64, true),
            field(Self::MOBILITY_COMPRESSION, DataType::Int64, true),
            field(Self::INTENSITY_COMPRESSION, DataType::Int64, true),
        ])
    }
}
