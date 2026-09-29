// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Shared Arrow/Parquet helpers for the OpenSWATH readers.
//!
//! Ports the private `FORMAT/DATAACCESS/ParquetReaderHelpers.h` /
//! `ParquetReaderHelpers.cpp` pair (namespace
//! `OpenMS::Internal::ParquetReaderHelpers`) together with the pieces the
//! three readers built on it repeat verbatim: the binary array codec
//! (`decodeBinary_`), the analyte extraction of `getAnalytes` and the run
//! listing of `getRuns`. The source keeps these out of the public headers;
//! so does this port — the module is crate-private, and only the plain data
//! types ([`Analyte`], [`RunInfo`], [`ReaderLimits`]) surface, re-exported
//! by [`xic_parquet_file`](crate::format::xic_parquet_file),
//! [`xim_parquet_file`](crate::format::xim_parquet_file) and
//! [`xipm_parquet_file`](crate::format::xipm_parquet_file).
//!
//! # Filter evaluation
//!
//! The source first tries an Arrow Dataset scan with the filter pushed down
//! (only when built `WITH_ARROW_DATASET`), then Arrow compute, and falls back
//! to its own row-by-row evaluator when Arrow cannot bind or execute the
//! expression. Every expression the parser builds is a left-to-right chain of
//! comparisons joined by AND/OR without negation, so Kleene logic with a null
//! row dropped selects exactly the rows the row-by-row evaluator selects with
//! a null comparison counted as false. This port therefore evaluates every
//! filter row by row, with the operand rules that the Arrow path and its
//! fallback produce together:
//!
//! - an integer condition compares against any signed or unsigned integer
//!   column or a boolean one (as 1/0); an unsigned value above `i64::MAX`
//!   is an error, as in the fallback;
//! - a string condition compares against a `string`, `large_string`,
//!   `binary` or `large_binary` column, byte-wise;
//! - any other column type is an error once a non-null value has to be
//!   compared ("Unsupported integer column type" / "Unsupported string column
//!   type"), which is what the fallback throws after Arrow refuses to bind.
//!
//! Native difference: an integer condition on a floating-point column, which
//! Arrow compute would evaluate after casting the literal to `double`, is
//! refused here. No OpenSWATH Parquet schema has a filterable floating-point
//! column: XIM's `FEATURE_RT` is pruned before evaluation, and XIPM exposes no
//! filter.

use crate::format::numpress_coder::{
    self as codec, NumpressCoderLimits, NumpressCompression, Work,
};
use crate::format::parquet_file::{self, ReadLimits};
use crate::format::parquet_filter::{ColumnType, Condition, FilterExpression};
use crate::{Error, Result};
use arrow_array::cast::AsArray;
use arrow_array::types::{
    Float32Type, Float64Type, Int8Type, Int16Type, Int32Type, Int64Type, UInt8Type, UInt16Type,
    UInt32Type, UInt64Type,
};
use arrow_array::{Array, ArrayRef, BooleanArray, RecordBatch};
use arrow_schema::{DataType, Schema, SchemaRef};
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

/// Resource ceilings for the OpenSWATH Parquet readers.
///
/// Native: the source reads whatever the files declare and decodes whatever
/// the binary cells hold. A reader checks every file's footer against
/// [`file`](Self::file) before it decodes any data, refuses a combined table
/// with more than [`max_total_rows`](Self::max_total_rows) rows, decodes each
/// binary cell within [`array`](Self::array) (which bounds zlib inflation as
/// well as the Numpress output) and stops once the decoded arrays of one call
/// together exceed [`max_total_values`](Self::max_total_values).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReaderLimits {
    /// Ceilings for each file's footer.
    pub file: ReadLimits,
    /// Maximum rows of all files together.
    pub max_total_rows: u64,
    /// Ceilings for decoding one binary cell.
    pub array: NumpressCoderLimits,
    /// Maximum decoded values of one call, summed over every array.
    pub max_total_values: usize,
}

impl Default for ReaderLimits {
    /// Default [`ReadLimits`] per file, 2^32 rows in total, default Numpress
    /// coder limits per cell (10,000,000 values, 128 MiB) and 2^31 decoded
    /// values per call.
    fn default() -> Self {
        Self {
            file: ReadLimits::default(),
            max_total_rows: 1 << 32,
            array: NumpressCoderLimits::default(),
            max_total_values: 1 << 31,
        }
    }
}

/// One unique run of an OpenSWATH Parquet file.
///
/// The source's `XICRunInfo`, `XIMRunInfo` and `XIPMRunInfo`, which are
/// member-for-member identical; the readers return this one type.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct RunInfo {
    /// The run identifier, the `RUN_ID` column.
    pub run_id: i64,
    /// The source file, the `SOURCE_FILE` column (empty when null or absent).
    pub source_file: String,
}

/// Analyte metadata of an OpenSWATH chromatogram or mobilogram file.
///
/// The source's `XICAnalyte` and `XIMAnalyte`, which are member-for-member
/// identical. Each `has_x` / `x` pair of the source is one `Option` here.
///
/// When the analytes are listed without nesting, every row is one unique
/// precursor–transition pair and the transition-level fields are the scalar
/// members ([`transition_id`](Self::transition_id) to
/// [`annotation`](Self::annotation)). With nesting, every row is one unique
/// precursor, the scalar transition members stay unset, and the vector
/// members ([`transition_ids`](Self::transition_ids) to
/// [`annotations`](Self::annotations)) hold one entry per unique transition
/// of that precursor, in file order; an unset integer is stored as `-1`, an
/// unset string as `""`, as in the source.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct Analyte {
    /// `PRECURSOR_ID`.
    pub precursor_id: Option<i64>,
    /// `MODIFIED_SEQUENCE` (empty when null or not requested).
    pub modified_sequence: String,
    /// `PRECURSOR_CHARGE`.
    pub precursor_charge: Option<i64>,
    /// `PRECURSOR_DECOY`.
    pub precursor_decoy: Option<i64>,
    /// `TRANSITION_ID` (unnested listings only).
    pub transition_id: Option<i64>,
    /// `PRODUCT_CHARGE` (unnested listings only).
    pub product_charge: Option<i64>,
    /// `TRANSITION_ORDINAL` (unnested listings only).
    pub transition_ordinal: Option<i64>,
    /// `DETECTING_TRANSITION` (unnested listings only).
    pub detecting_transition: Option<i64>,
    /// `PRODUCT_DECOY` (unnested listings only).
    pub product_decoy: Option<i64>,
    /// `TRANSITION_TYPE` (unnested listings only).
    pub transition_type: String,
    /// `ANNOTATION` (unnested listings only).
    pub annotation: String,
    /// `TRANSITION_ID` of each transition (nested listings only).
    pub transition_ids: Vec<i64>,
    /// `PRODUCT_CHARGE` of each transition (nested listings only).
    pub product_charges: Vec<i64>,
    /// `TRANSITION_ORDINAL` of each transition (nested listings only).
    pub transition_ordinals: Vec<i64>,
    /// `DETECTING_TRANSITION` of each transition (nested listings only).
    pub detecting_transitions: Vec<i64>,
    /// `PRODUCT_DECOY` of each transition (nested listings only).
    pub product_decoys: Vec<i64>,
    /// `TRANSITION_TYPE` of each transition (nested listings only).
    pub transition_types: Vec<String>,
    /// `ANNOTATION` of each transition (nested listings only).
    pub annotations: Vec<String>,
}

/// The source's `FilterPruningOptions`: how a filter is pruned against a
/// schema.
#[derive(Clone, Copy, Debug)]
pub(crate) struct FilterPruningOptions {
    /// Upper-case column names that exist but cannot be evaluated.
    pub(crate) unsupported_columns: &'static [&'static str],
    /// Refuse a filter that names an unsupported column.
    pub(crate) reject_unsupported_after_pruning: bool,
    /// Refuse a filter that pruning left empty.
    pub(crate) reject_empty_after_pruning: bool,
}

/// The source's `upper_`: ASCII upper case, byte by byte (`std::toupper` in
/// the C locale).
pub(crate) fn upper(input: &str) -> String {
    input.to_ascii_uppercase()
}

/// Check that every file exists, as the readers' constructors do.
///
/// # Errors
///
/// [`Error::Io`] of kind `NotFound` (the source's `Exception::FileNotFound`)
/// naming the first missing file.
pub(crate) fn require_files(filenames: &[PathBuf]) -> Result<()> {
    for filename in filenames {
        if !filename.exists() {
            return Err(Error::Io(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("file not found: {}", filename.display()),
            )));
        }
    }
    Ok(())
}

/// Two schemas are equal as Arrow C++'s `Schema::Equals(other,
/// check_metadata=false)` sees them: same field names, types and
/// nullability, in order.
fn schemas_equal(a: &Schema, b: &Schema) -> bool {
    a.fields().len() == b.fields().len()
        && a.fields().iter().zip(b.fields().iter()).all(|(x, y)| {
            x.name() == y.name()
                && x.data_type() == y.data_type()
                && x.is_nullable() == y.is_nullable()
        })
}

/// `readParquetTable(filenames)` / `readParquetTableColumns(filenames,
/// columns)`: every file read (`columns` projects), then concatenated.
///
/// # Errors
///
/// As [`parquet_file::read_file_columns`] for each file; [`Error::InvalidValue`]
/// "No parquet tables provided" for no file, "Failed to concatenate parquet
/// tables" when the files' schemas differ, and when the rows exceed
/// `limits.max_total_rows`.
pub(crate) fn read_tables(
    filenames: &[PathBuf],
    columns: Option<&[&str]>,
    limits: &ReaderLimits,
) -> Result<RecordBatch> {
    let mut tables = Vec::with_capacity(filenames.len());
    let mut rows: u64 = 0;
    for filename in filenames {
        let table = parquet_file::read_file_columns(filename, columns, &limits.file)?;
        rows = rows.saturating_add(table.num_rows() as u64);
        if rows > limits.max_total_rows {
            return Err(Error::InvalidValue(format!(
                "Parquet inputs hold more than {} rows",
                limits.max_total_rows
            )));
        }
        tables.push(table);
    }
    let first = match tables.first() {
        Some(table) => table.schema(),
        None => return Err(Error::InvalidValue("No parquet tables provided".into())),
    };
    if tables.len() == 1 {
        return Ok(tables.swap_remove(0));
    }
    for table in &tables[1..] {
        if !schemas_equal(&first, &table.schema()) {
            return Err(Error::InvalidValue(
                "Failed to concatenate parquet tables: schemas differ".into(),
            ));
        }
    }
    let batches: Vec<RecordBatch> = tables
        .into_iter()
        .map(|table| RecordBatch::try_new(first.clone(), table.columns().to_vec()))
        .collect::<std::result::Result<_, _>>()
        .map_err(|error| {
            Error::InvalidValue(format!("Failed to concatenate parquet tables: {error}"))
        })?;
    arrow_select::concat::concat_batches(&first, &batches).map_err(|error| {
        Error::InvalidValue(format!("Failed to concatenate parquet tables: {error}"))
    })
}

/// `readParquetSchemaAllFiles`: the first file's schema, provided every other
/// file's schema equals it.
///
/// # Errors
///
/// [`Error::InvalidValue`] "No parquet files provided for schema validation",
/// "Failed to read parquet schema", or "Parquet input files have incompatible
/// schemas: '<first>' vs '<other>'".
pub(crate) fn read_schema_all_files(filenames: &[PathBuf]) -> Result<SchemaRef> {
    let (first, rest) = filenames.split_first().ok_or_else(|| {
        Error::InvalidValue("No parquet files provided for schema validation".into())
    })?;
    let base = parquet_file::read_file_schema(first)?;
    for other in rest {
        let schema = parquet_file::read_file_schema(other)?;
        if !schemas_equal(&base, &schema) {
            return Err(Error::InvalidValue(format!(
                "Parquet input files have incompatible schemas: '{}' vs '{}'",
                first.display(),
                other.display()
            )));
        }
    }
    Ok(base)
}

/// `getColumn(table, name, required)`: the named column, `None` for an
/// absent optional one.
///
/// A name that occurs more than once is absent, as for Arrow C++'s
/// `Table::GetColumnByName`.
///
/// # Errors
///
/// [`Error::MissingInformation`] "Missing required column '<name>'".
pub(crate) fn column(table: &RecordBatch, name: &str, required: bool) -> Result<Option<ArrayRef>> {
    let schema = table.schema();
    let mut matches = schema
        .fields()
        .iter()
        .enumerate()
        .filter(|(_, field)| field.name() == name);
    match (matches.next(), matches.next()) {
        (Some((index, _)), None) => Ok(Some(table.column(index).clone())),
        _ if required => Err(Error::MissingInformation(format!(
            "Missing required column '{name}'"
        ))),
        _ => Ok(None),
    }
}

/// `getOptionalInt`: the integer at `row`, `None` for a null cell or an absent
/// column. Booleans read as 1/0.
///
/// # Errors
///
/// [`Error::InvalidValue`] "UINT64 value exceeds Int64 range" and
/// "Unsupported integer column type".
pub(crate) fn optional_int(array: Option<&ArrayRef>, row: usize) -> Result<Option<i64>> {
    let array = match array {
        Some(array) if !array.is_null(row) => array,
        _ => return Ok(None),
    };
    Ok(Some(match array.data_type() {
        DataType::Int64 => array.as_primitive::<Int64Type>().value(row),
        DataType::Int32 => i64::from(array.as_primitive::<Int32Type>().value(row)),
        DataType::Int16 => i64::from(array.as_primitive::<Int16Type>().value(row)),
        DataType::Int8 => i64::from(array.as_primitive::<Int8Type>().value(row)),
        DataType::UInt64 => {
            let value = array.as_primitive::<UInt64Type>().value(row);
            i64::try_from(value).map_err(|_| {
                Error::InvalidValue(format!("UINT64 value exceeds Int64 range: {value}"))
            })?
        }
        DataType::UInt32 => i64::from(array.as_primitive::<UInt32Type>().value(row)),
        DataType::UInt16 => i64::from(array.as_primitive::<UInt16Type>().value(row)),
        DataType::UInt8 => i64::from(array.as_primitive::<UInt8Type>().value(row)),
        DataType::Boolean => i64::from(array.as_boolean().value(row)),
        other => {
            return Err(Error::InvalidValue(format!(
                "Unsupported integer column type: {other}"
            )));
        }
    }))
}

/// `getOptionalDouble`: the `double` or `float` at `row`.
///
/// # Errors
///
/// [`Error::InvalidValue`] "Unsupported floating-point column type".
pub(crate) fn optional_double(array: Option<&ArrayRef>, row: usize) -> Result<Option<f64>> {
    let array = match array {
        Some(array) if !array.is_null(row) => array,
        _ => return Ok(None),
    };
    match array.data_type() {
        DataType::Float64 => Ok(Some(array.as_primitive::<Float64Type>().value(row))),
        DataType::Float32 => Ok(Some(f64::from(
            array.as_primitive::<Float32Type>().value(row),
        ))),
        other => Err(Error::InvalidValue(format!(
            "Unsupported floating-point column type: {other}"
        ))),
    }
}

/// `getOptionalString`: the `string` or `large_string` at `row`.
///
/// # Errors
///
/// [`Error::InvalidValue`] "Unsupported string column type".
pub(crate) fn optional_string(array: Option<&ArrayRef>, row: usize) -> Result<Option<String>> {
    let array = match array {
        Some(array) if !array.is_null(row) => array,
        _ => return Ok(None),
    };
    match array.data_type() {
        DataType::Utf8 => Ok(Some(array.as_string::<i32>().value(row).to_owned())),
        DataType::LargeUtf8 => Ok(Some(array.as_string::<i64>().value(row).to_owned())),
        other => Err(Error::InvalidValue(format!(
            "Unsupported string column type: {other}"
        ))),
    }
}

/// [`optional_string`] with the source's out-parameter behaviour: a null cell
/// leaves the string empty.
pub(crate) fn string_or_empty(array: Option<&ArrayRef>, row: usize) -> Result<String> {
    Ok(optional_string(array, row)?.unwrap_or_default())
}

/// `getBinaryView`: the bytes at `row`, empty for a null cell or an absent
/// column.
///
/// # Errors
///
/// [`Error::InvalidValue`] "Unsupported binary column type".
pub(crate) fn binary(array: Option<&ArrayRef>, row: usize) -> Result<&[u8]> {
    let array = match array {
        Some(array) if !array.is_null(row) => array,
        _ => return Ok(&[]),
    };
    match array.data_type() {
        DataType::Binary => Ok(array.as_binary::<i32>().value(row)),
        DataType::LargeBinary => Ok(array.as_binary::<i64>().value(row)),
        other => Err(Error::InvalidValue(format!(
            "Unsupported binary column type: {other}"
        ))),
    }
}

/// The source's `tokenize_`.
fn tokenize(expr: &str) -> Result<Vec<String>> {
    let bytes = expr.as_bytes();
    let mut tokens = Vec::new();
    let mut current: Vec<u8> = Vec::new();
    let mut quote: Option<u8> = None;
    let flush = |current: &mut Vec<u8>, tokens: &mut Vec<Vec<u8>>| {
        if !current.is_empty() {
            tokens.push(std::mem::take(current));
        }
    };
    let mut raw: Vec<Vec<u8>> = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        let c = bytes[i];
        if let Some(q) = quote {
            if c == q {
                // A closing quote always yields a token, even an empty one.
                raw.push(std::mem::take(&mut current));
                quote = None;
            } else {
                current.push(c);
            }
            i += 1;
            continue;
        }
        match c {
            b'"' | b'\'' => {
                flush(&mut current, &mut raw);
                quote = Some(c);
            }
            // std::isspace in the C locale.
            b' ' | b'\t' | b'\n' | b'\x0b' | b'\x0c' | b'\r' => flush(&mut current, &mut raw),
            b'[' | b']' | b',' | b'(' | b')' => {
                flush(&mut current, &mut raw);
                raw.push(vec![c]);
            }
            b'=' | b'!' | b'<' | b'>' => {
                flush(&mut current, &mut raw);
                let mut op = vec![c];
                if bytes.get(i + 1) == Some(&b'=') {
                    op.push(b'=');
                    i += 1;
                }
                raw.push(op);
            }
            b'&' | b'|' => {
                flush(&mut current, &mut raw);
                let mut op = vec![c];
                if bytes.get(i + 1) == Some(&c) {
                    op.push(c);
                    i += 1;
                }
                raw.push(op);
            }
            _ => current.push(c),
        }
        i += 1;
    }
    flush(&mut current, &mut raw);
    if quote.is_some() {
        return Err(Error::InvalidValue(format!(
            "Unclosed quote in filter expression: {expr}"
        )));
    }
    // Tokens split only at ASCII bytes, so each is valid UTF-8.
    for token in raw {
        tokens.push(String::from_utf8(token).map_err(|_| {
            Error::InvalidValue(format!("Invalid UTF-8 in filter expression: {expr}"))
        })?);
    }
    Ok(tokens)
}

/// `parseFilter`: a textual filter into a typed expression.
///
/// Column names are upper-cased; a column missing from `column_types` is an
/// integer column, as in the source.
///
/// # Errors
///
/// [`Error::InvalidValue`] with the source's messages: "Unclosed quote in
/// filter expression", "Missing filter operator for column", "Unsupported
/// filter operator", "IN operator expects a bracketed value list", "IN
/// operator expects at least one value", "Unclosed IN list", "Missing filter
/// value for operator" and "Dangling connector in filter expression".
pub(crate) fn parse_filter(
    filter: &str,
    column_types: &[(&str, ColumnType)],
) -> Result<FilterExpression> {
    let mut expr = FilterExpression::default();
    if filter.is_empty() {
        return Ok(expr);
    }
    let tokens = tokenize(filter)?;
    let mut i = 0;
    while i < tokens.len() {
        let column = upper(&tokens[i]);
        i += 1;
        if i >= tokens.len() {
            return Err(Error::InvalidValue(format!(
                "Missing filter operator for column: {column}"
            )));
        }
        let mut op = upper(&tokens[i]);
        i += 1;
        if op == "==" {
            op = "=".into();
        }
        if !matches!(op.as_str(), "=" | "!=" | "<" | "<=" | ">" | ">=" | "IN") {
            return Err(Error::InvalidValue(format!(
                "Unsupported filter operator: {op}"
            )));
        }
        let column_type = column_types
            .iter()
            .find(|(name, _)| *name == column)
            .map_or(ColumnType::Int, |(_, column_type)| *column_type);
        let mut values = Vec::new();
        if op == "IN" {
            let close = match tokens.get(i).map(String::as_str) {
                Some("[") => "]",
                Some("(") => ")",
                _ => {
                    return Err(Error::InvalidValue(format!(
                        "IN operator expects a bracketed value list: {filter}"
                    )));
                }
            };
            i += 1;
            while i < tokens.len() && tokens[i] != close {
                if tokens[i] != "," {
                    values.push(tokens[i].clone());
                }
                i += 1;
            }
            if values.is_empty() {
                return Err(Error::InvalidValue(format!(
                    "IN operator expects at least one value: {filter}"
                )));
            }
            if i >= tokens.len() {
                return Err(Error::InvalidValue(format!("Unclosed IN list: {filter}")));
            }
            i += 1;
        } else {
            match tokens.get(i) {
                Some(value) => values.push(value.clone()),
                None => {
                    return Err(Error::InvalidValue(format!(
                        "Missing filter value for operator: {op}"
                    )));
                }
            }
            i += 1;
        }
        expr.conditions.push(Condition {
            column,
            op,
            values,
            column_type,
        });
        if let Some(token) = tokens.get(i) {
            let connector = upper(token);
            let consumed = match connector.as_str() {
                "AND" | "&&" | "&" => {
                    expr.connectors.push("AND".into());
                    true
                }
                "OR" | "||" | "|" => {
                    expr.connectors.push("OR".into());
                    true
                }
                _ => false,
            };
            if consumed {
                i += 1;
                if i >= tokens.len() {
                    return Err(Error::InvalidValue(format!(
                        "Dangling connector in filter expression: {filter}"
                    )));
                }
            }
        }
    }
    Ok(expr)
}

/// `normalizeAndPruneColumns_` + `normalizeAndValidateFilter_`: map each
/// condition's column onto the schema's spelling, drop unsupported ones, and
/// report whether any condition is left.
fn normalize_and_validate(
    expr: &mut FilterExpression,
    schema: &Schema,
    filter_context: &str,
    options: &FilterPruningOptions,
) -> Result<bool> {
    if expr.conditions.is_empty() {
        return Ok(false);
    }
    // Later fields win when two share an upper-case name, as in the source.
    let mut name_map: HashMap<String, String> = HashMap::new();
    for field in schema.fields() {
        name_map.insert(upper(field.name()), field.name().clone());
    }
    let mut unknown = Vec::new();
    let mut unsupported = Vec::new();
    let mut kept: Vec<Condition> = Vec::with_capacity(expr.conditions.len());
    let mut connectors = Vec::new();
    for (i, condition) in expr.conditions.iter().enumerate() {
        let upper_column = upper(&condition.column);
        let Some(actual) = name_map.get(&upper_column) else {
            unknown.push(condition.column.clone());
            continue;
        };
        if options.unsupported_columns.contains(&upper_column.as_str()) {
            unsupported.push(condition.column.clone());
            continue;
        }
        if !kept.is_empty() {
            connectors.push(
                i.checked_sub(1)
                    .and_then(|previous| expr.connectors.get(previous))
                    .cloned()
                    .unwrap_or_else(|| "AND".into()),
            );
        }
        let mut condition = condition.clone();
        condition.column = actual.clone();
        kept.push(condition);
    }
    if !unknown.is_empty() {
        return Err(Error::InvalidValue(format!(
            "Unknown filter columns: {}",
            unknown.join(", ")
        )));
    }
    if !unsupported.is_empty() && options.reject_unsupported_after_pruning {
        return Err(Error::InvalidValue(format!(
            "Unsupported filter columns after pruning: {}",
            unsupported.join(", ")
        )));
    }
    expr.conditions = kept;
    expr.connectors = connectors;
    if expr.conditions.is_empty() {
        if options.reject_empty_after_pruning {
            return Err(Error::InvalidValue(format!(
                "Filter expression became empty after pruning unsupported columns: {filter_context}"
            )));
        }
        return Ok(false);
    }
    Ok(true)
}

/// A condition bound to its column, integer values parsed once.
struct BoundCondition<'a> {
    condition: &'a Condition,
    array: &'a ArrayRef,
    ints: Vec<i64>,
}

/// The value of an integer filter literal (`StringUtils::toInt64`).
fn parse_int(value: &str) -> Result<i64> {
    source_to_int64(value)
        .ok_or_else(|| Error::InvalidValue(format!("Invalid integer filter value: {value}")))
}

/// `StringUtils::toInt64`: leading and trailing spaces, tabs, line feeds and
/// carriage returns are skipped, one leading `+` is dropped, and the rest
/// must be `std::from_chars` input for an `Int64` (an optional `-`, then
/// decimal digits, in range). `None` where the source throws
/// `Exception::ConversionError`.
///
/// As in the source, `"+-5"` reads as -5: the `+` is dropped before
/// `from_chars`, which accepts the `-`.
pub(crate) fn source_to_int64(text: &str) -> Option<i64> {
    let is_ws = |c: char| matches!(c, ' ' | '\t' | '\n' | '\r');
    let trimmed = text.trim_start_matches(is_ws);
    if trimmed.is_empty() {
        return None;
    }
    let body = trimmed.strip_prefix('+').unwrap_or(trimmed);
    let digits_start = usize::from(body.starts_with('-'));
    let digits = body[digits_start..]
        .bytes()
        .take_while(u8::is_ascii_digit)
        .count();
    if digits == 0 {
        return None;
    }
    let end = digits_start + digits;
    if !body[end..].chars().all(is_ws) {
        return None;
    }
    body[..end].parse::<i64>().ok()
}

fn compare<T: PartialOrd + ?Sized>(op: &str, lhs: &T, rhs: &T) -> bool {
    match op {
        "=" => lhs == rhs,
        "!=" => lhs != rhs,
        "<" => lhs < rhs,
        "<=" => lhs <= rhs,
        ">" => lhs > rhs,
        _ => lhs >= rhs,
    }
}

/// The source's `evaluateCondition_`.
fn evaluate(bound: &BoundCondition<'_>, row: usize) -> Result<bool> {
    let condition = bound.condition;
    if condition.column_type == ColumnType::Int {
        let lhs = match bound.array.data_type() {
            DataType::Float16 | DataType::Float32 | DataType::Float64
                if !bound.array.is_null(row) =>
            {
                return Err(Error::InvalidValue(format!(
                    "Unsupported integer column type: {}",
                    bound.array.data_type()
                )));
            }
            _ => match optional_int(Some(bound.array), row)? {
                Some(value) => value,
                None => return Ok(false),
            },
        };
        if condition.op == "IN" {
            return Ok(bound.ints.contains(&lhs));
        }
        return Ok(match bound.ints.first() {
            Some(rhs) => compare(&condition.op, &lhs, rhs),
            None => false,
        });
    }
    let array = bound.array;
    if array.is_null(row) {
        return Ok(false);
    }
    let lhs: &[u8] = match array.data_type() {
        DataType::Utf8 => array.as_string::<i32>().value(row).as_bytes(),
        DataType::LargeUtf8 => array.as_string::<i64>().value(row).as_bytes(),
        DataType::Binary => array.as_binary::<i32>().value(row),
        DataType::LargeBinary => array.as_binary::<i64>().value(row),
        other => {
            return Err(Error::InvalidValue(format!(
                "Unsupported string column type: {other}"
            )));
        }
    };
    if condition.op == "IN" {
        return Ok(condition.values.iter().any(|value| value.as_bytes() == lhs));
    }
    Ok(match condition.values.first() {
        Some(rhs) => compare(&condition.op, lhs, rhs.as_bytes()),
        None => false,
    })
}

/// `applyArrowFilter`: the rows of `table` the expression selects, after
/// pruning it against the table's schema.
///
/// # Errors
///
/// [`Error::InvalidValue`] for unknown columns, pruning refusals, an operator
/// outside `=`, `!=`, `<`, `<=`, `>`, `>=`, `IN`, an integer literal that is
/// not one, and operand types the source cannot compare (see the module
/// documentation).
pub(crate) fn apply_filter(
    table: RecordBatch,
    parsed: &FilterExpression,
    filter_context: &str,
    options: &FilterPruningOptions,
) -> Result<RecordBatch> {
    if parsed.conditions.is_empty() {
        return Ok(table);
    }
    let mut pruned = parsed.clone();
    let schema = table.schema();
    if !normalize_and_validate(&mut pruned, &schema, filter_context, options)? {
        return Ok(table);
    }
    let mut bound = Vec::with_capacity(pruned.conditions.len());
    for condition in &pruned.conditions {
        if !matches!(
            condition.op.as_str(),
            "=" | "!=" | "<" | "<=" | ">" | ">=" | "IN"
        ) {
            return Err(Error::InvalidValue(format!(
                "Unsupported filter operator: {}",
                condition.op
            )));
        }
        let ints = if condition.column_type == ColumnType::Int {
            condition
                .values
                .iter()
                .map(|value| parse_int(value))
                .collect::<Result<Vec<_>>>()?
        } else {
            Vec::new()
        };
        let array = column(&table, &condition.column, false)?.ok_or_else(|| {
            Error::InvalidValue(format!(
                "Column not found in parquet schema: {}",
                condition.column
            ))
        })?;
        bound.push((condition, array, ints));
    }
    let bound: Vec<BoundCondition<'_>> = bound
        .iter()
        .map(|(condition, array, ints)| BoundCondition {
            condition,
            array,
            ints: ints.clone(),
        })
        .collect();
    let mut mask = Vec::with_capacity(table.num_rows());
    for row in 0..table.num_rows() {
        let mut result = evaluate(&bound[0], row)?;
        for (i, condition) in bound.iter().enumerate().skip(1) {
            let next = evaluate(condition, row)?;
            if pruned.connectors.get(i - 1).map(String::as_str) == Some("OR") {
                result = result || next;
            } else {
                result = result && next;
            }
        }
        mask.push(result);
    }
    arrow_select::filter::filter_record_batch(&table, &BooleanArray::from(mask)).map_err(|error| {
        Error::InvalidValue(format!(
            "Failed to build filtered table: {filter_context}: {error}"
        ))
    })
}

/// Which compression identifiers a reader accepts and how it words its
/// errors; the three readers' `decodeBinary_` differ only in these.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Codec {
    /// Whether identifier 7 (MS-Numpress PIC + zlib) is accepted.
    pub(crate) pic_zlib: bool,
    /// "RT/intensity", "mobility/intensity" or "binary".
    pub(crate) label: &'static str,
    /// "Unsupported compression type" or "Unknown compression code".
    pub(crate) unknown: &'static str,
}

/// Running total of decoded values for one call.
pub(crate) struct DecodeBudget {
    limits: ReaderLimits,
    values: usize,
}

impl DecodeBudget {
    pub(crate) fn new(limits: ReaderLimits) -> Self {
        Self { limits, values: 0 }
    }
}

fn bytes_to_doubles(bytes: &[u8], what: &str, budget: &mut DecodeBudget) -> Result<Vec<f64>> {
    if bytes.len() % 8 != 0 {
        return Err(Error::InvalidValue(format!("{what}: {}", bytes.len())));
    }
    let count = bytes.len() / 8;
    charge(budget, count)?;
    // The source copies the bytes with memcpy, so they are host order; every
    // producer the source supports is little-endian.
    Ok(bytes
        .chunks_exact(8)
        .map(|chunk| {
            let mut value = [0u8; 8];
            value.copy_from_slice(chunk);
            f64::from_le_bytes(value)
        })
        .collect())
}

fn charge(budget: &mut DecodeBudget, count: usize) -> Result<()> {
    budget.values = budget.values.saturating_add(count);
    if budget.values > budget.limits.max_total_values {
        return Err(Error::InvalidValue(format!(
            "Decoded arrays exceed {} values",
            budget.limits.max_total_values
        )));
    }
    Ok(())
}

/// `decodeBinary_`: decode one binary cell by its compression identifier.
///
/// 0 raw little-endian doubles, 1 zlib-compressed doubles, 2/3/4 MS-Numpress
/// linear/slof/PIC, 5/6/7 the same after zlib. Empty data decodes to an empty
/// array whatever the identifier says.
///
/// # Errors
///
/// [`Error::InvalidValue`] for an unknown identifier, a byte count that is
/// not a multiple of eight, a corrupt zlib or Numpress payload, and the
/// ceilings of `budget`.
pub(crate) fn decode(
    data: &[u8],
    compression: i64,
    codec_spec: Codec,
    budget: &mut DecodeBudget,
) -> Result<Vec<f64>> {
    if data.is_empty() {
        return Ok(Vec::new());
    }
    let mut work = Work::new(budget.limits.array);
    let numpress = |bytes: &[u8], mode, work: &mut Work| codec::decode_raw(bytes, mode, work);
    let values = match compression {
        0 => {
            return bytes_to_doubles(
                data,
                &format!(
                    "Invalid {} data size (not divisible by 8)",
                    codec_spec.label
                ),
                budget,
            );
        }
        1 => {
            let decoded = codec::zlib_decode(data, &mut work)?;
            return bytes_to_doubles(
                &decoded,
                &format!("Invalid decompressed {} data size", codec_spec.label),
                budget,
            );
        }
        2 => numpress(data, NumpressCompression::Linear, &mut work)?,
        3 => numpress(data, NumpressCompression::Slof, &mut work)?,
        4 => numpress(data, NumpressCompression::Pic, &mut work)?,
        5 | 6 => {
            let decoded = codec::zlib_decode(data, &mut work)?;
            let mode = if compression == 5 {
                NumpressCompression::Linear
            } else {
                NumpressCompression::Slof
            };
            numpress(&decoded, mode, &mut work)?
        }
        7 if codec_spec.pic_zlib => {
            let decoded = codec::zlib_decode(data, &mut work)?;
            numpress(&decoded, NumpressCompression::Pic, &mut work)?
        }
        other => {
            return Err(Error::InvalidValue(format!(
                "{}: {other}",
                codec_spec.unknown
            )));
        }
    };
    charge(budget, values.len())?;
    Ok(values)
}

/// The source's `appendKeyPart_` overloads, building a composite key.
fn key_int(key: &mut String, value: Option<i64>) {
    match value {
        Some(value) => key.push_str(&value.to_string()),
        None => key.push_str("NULL"),
    }
    key.push('|');
}

fn key_str(key: &mut String, value: &str) {
    key.push_str(value);
    key.push('|');
}

/// The analyte columns `getAnalytes` accepts, which are also its default.
pub(crate) const ANALYTE_COLUMNS: [&str; 11] = [
    "PRECURSOR_ID",
    "MODIFIED_SEQUENCE",
    "PRECURSOR_CHARGE",
    "PRECURSOR_DECOY",
    "TRANSITION_ID",
    "PRODUCT_CHARGE",
    "TRANSITION_ORDINAL",
    "DETECTING_TRANSITION",
    "PRODUCT_DECOY",
    "TRANSITION_TYPE",
    "ANNOTATION",
];

/// `getAnalytes` of the XIC and XIM readers.
///
/// `require_discriminator` is XIM's extra check that a nested listing
/// requests at least one precursor-level column.
///
/// # Errors
///
/// [`Error::InvalidValue`] "Unsupported analyte column" and "Column not found
/// in parquet schema" (naming the column as given), XIM's "nest_transitions=true
/// requires at least one precursor discriminator column", and the read and
/// value errors of the helpers above.
pub(crate) fn analytes<S: AsRef<str>>(
    filenames: &[PathBuf],
    columns: &[S],
    nest_transitions: bool,
    require_discriminator: bool,
    limits: &ReaderLimits,
) -> Result<Vec<Analyte>> {
    let requested: Vec<&str> = if columns.is_empty() {
        ANALYTE_COLUMNS.to_vec()
    } else {
        columns.iter().map(AsRef::as_ref).collect()
    };
    let schema = read_schema_all_files(filenames)?;
    let schema_columns: HashSet<String> = schema
        .fields()
        .iter()
        .map(|field| upper(field.name()))
        .collect();
    let mut normalized = Vec::with_capacity(requested.len());
    for name in &requested {
        let upper_name = upper(name);
        if !ANALYTE_COLUMNS.contains(&upper_name.as_str()) {
            return Err(Error::InvalidValue(format!(
                "Unsupported analyte column: {name}"
            )));
        }
        if !schema_columns.contains(&upper_name) {
            return Err(Error::InvalidValue(format!(
                "Column not found in parquet schema: {name}"
            )));
        }
        normalized.push(upper_name);
    }
    let wanted: HashSet<&str> = normalized.iter().map(String::as_str).collect();
    let want = |name: &str| wanted.contains(name);
    if nest_transitions
        && require_discriminator
        && ![
            "PRECURSOR_ID",
            "MODIFIED_SEQUENCE",
            "PRECURSOR_CHARGE",
            "PRECURSOR_DECOY",
        ]
        .iter()
        .any(|name| want(name))
    {
        return Err(Error::InvalidValue(
            "nest_transitions=true requires at least one precursor discriminator column: \
             PRECURSOR_ID, MODIFIED_SEQUENCE, PRECURSOR_CHARGE, or PRECURSOR_DECOY"
                .into(),
        ));
    }
    let projection: Vec<&str> = normalized.iter().map(String::as_str).collect();
    let table = read_tables(filenames, Some(&projection), limits)?;
    let col = |name: &str| column(&table, name, want(name));
    let precursor_id = col("PRECURSOR_ID")?;
    let modified_sequence = col("MODIFIED_SEQUENCE")?;
    let precursor_charge = col("PRECURSOR_CHARGE")?;
    let precursor_decoy = col("PRECURSOR_DECOY")?;
    let transition_id = col("TRANSITION_ID")?;
    let product_charge = col("PRODUCT_CHARGE")?;
    let transition_ordinal = col("TRANSITION_ORDINAL")?;
    let detecting_transition = col("DETECTING_TRANSITION")?;
    let product_decoy = col("PRODUCT_DECOY")?;
    let transition_type = col("TRANSITION_TYPE")?;
    let annotation = col("ANNOTATION")?;
    let int = |array: &Option<ArrayRef>, name: &str, row: usize| -> Result<Option<i64>> {
        if want(name) {
            optional_int(array.as_ref(), row)
        } else {
            Ok(None)
        }
    };
    let text = |array: &Option<ArrayRef>, name: &str, row: usize| -> Result<String> {
        if want(name) {
            string_or_empty(array.as_ref(), row)
        } else {
            Ok(String::new())
        }
    };
    let mut output: Vec<Analyte> = Vec::new();
    if !nest_transitions {
        let mut seen = HashSet::new();
        for row in 0..table.num_rows() {
            let analyte = Analyte {
                precursor_id: int(&precursor_id, "PRECURSOR_ID", row)?,
                modified_sequence: text(&modified_sequence, "MODIFIED_SEQUENCE", row)?,
                precursor_charge: int(&precursor_charge, "PRECURSOR_CHARGE", row)?,
                precursor_decoy: int(&precursor_decoy, "PRECURSOR_DECOY", row)?,
                transition_id: int(&transition_id, "TRANSITION_ID", row)?,
                product_charge: int(&product_charge, "PRODUCT_CHARGE", row)?,
                transition_ordinal: int(&transition_ordinal, "TRANSITION_ORDINAL", row)?,
                detecting_transition: int(&detecting_transition, "DETECTING_TRANSITION", row)?,
                product_decoy: int(&product_decoy, "PRODUCT_DECOY", row)?,
                transition_type: text(&transition_type, "TRANSITION_TYPE", row)?,
                annotation: text(&annotation, "ANNOTATION", row)?,
                ..Analyte::default()
            };
            let mut key = String::new();
            if want("PRECURSOR_ID") {
                key_int(&mut key, analyte.precursor_id);
            }
            if want("MODIFIED_SEQUENCE") {
                key_str(&mut key, &analyte.modified_sequence);
            }
            if want("PRECURSOR_CHARGE") {
                key_int(&mut key, analyte.precursor_charge);
            }
            if want("PRECURSOR_DECOY") {
                key_int(&mut key, analyte.precursor_decoy);
            }
            if want("TRANSITION_ID") {
                key_int(&mut key, analyte.transition_id);
            }
            if want("PRODUCT_CHARGE") {
                key_int(&mut key, analyte.product_charge);
            }
            if want("TRANSITION_ORDINAL") {
                key_int(&mut key, analyte.transition_ordinal);
            }
            if want("DETECTING_TRANSITION") {
                key_int(&mut key, analyte.detecting_transition);
            }
            if want("PRODUCT_DECOY") {
                key_int(&mut key, analyte.product_decoy);
            }
            if want("TRANSITION_TYPE") {
                key_str(&mut key, &analyte.transition_type);
            }
            if want("ANNOTATION") {
                key_str(&mut key, &analyte.annotation);
            }
            if seen.insert(key) {
                output.push(analyte);
            }
        }
        return Ok(output);
    }
    let mut precursor_index: HashMap<String, usize> = HashMap::new();
    let mut transition_seen: HashMap<String, HashSet<String>> = HashMap::new();
    for row in 0..table.num_rows() {
        // The nested listing reads the precursor columns unconditionally; an
        // unrequested one is absent from the projected table and reads unset.
        let analyte = Analyte {
            precursor_id: optional_int(precursor_id.as_ref(), row)?,
            modified_sequence: string_or_empty(modified_sequence.as_ref(), row)?,
            precursor_charge: optional_int(precursor_charge.as_ref(), row)?,
            precursor_decoy: optional_int(precursor_decoy.as_ref(), row)?,
            ..Analyte::default()
        };
        let mut precursor_key = String::new();
        if want("PRECURSOR_ID") {
            key_int(&mut precursor_key, analyte.precursor_id);
        }
        if want("MODIFIED_SEQUENCE") {
            key_str(&mut precursor_key, &analyte.modified_sequence);
        }
        if want("PRECURSOR_CHARGE") {
            key_int(&mut precursor_key, analyte.precursor_charge);
        }
        if want("PRECURSOR_DECOY") {
            key_int(&mut precursor_key, analyte.precursor_decoy);
        }
        let index = match precursor_index.get(&precursor_key) {
            Some(&index) => index,
            None => {
                output.push(analyte);
                precursor_index.insert(precursor_key.clone(), output.len() - 1);
                output.len() - 1
            }
        };
        let t_id = int(&transition_id, "TRANSITION_ID", row)?;
        let t_charge = int(&product_charge, "PRODUCT_CHARGE", row)?;
        let t_ordinal = int(&transition_ordinal, "TRANSITION_ORDINAL", row)?;
        let t_detecting = int(&detecting_transition, "DETECTING_TRANSITION", row)?;
        let t_decoy = int(&product_decoy, "PRODUCT_DECOY", row)?;
        let t_type = text(&transition_type, "TRANSITION_TYPE", row)?;
        let t_annotation = text(&annotation, "ANNOTATION", row)?;
        let mut transition_key = String::new();
        if want("TRANSITION_ID") {
            key_int(&mut transition_key, t_id);
        }
        if want("PRODUCT_CHARGE") {
            key_int(&mut transition_key, t_charge);
        }
        if want("TRANSITION_ORDINAL") {
            key_int(&mut transition_key, t_ordinal);
        }
        if want("DETECTING_TRANSITION") {
            key_int(&mut transition_key, t_detecting);
        }
        if want("PRODUCT_DECOY") {
            key_int(&mut transition_key, t_decoy);
        }
        if want("TRANSITION_TYPE") {
            key_str(&mut transition_key, &t_type);
        }
        if want("ANNOTATION") {
            key_str(&mut transition_key, &t_annotation);
        }
        if transition_seen
            .entry(precursor_key)
            .or_default()
            .insert(transition_key)
        {
            let target = &mut output[index];
            if want("TRANSITION_ID") {
                target.transition_ids.push(t_id.unwrap_or(-1));
            }
            if want("PRODUCT_CHARGE") {
                target.product_charges.push(t_charge.unwrap_or(-1));
            }
            if want("TRANSITION_ORDINAL") {
                target.transition_ordinals.push(t_ordinal.unwrap_or(-1));
            }
            if want("DETECTING_TRANSITION") {
                target.detecting_transitions.push(t_detecting.unwrap_or(-1));
            }
            if want("PRODUCT_DECOY") {
                target.product_decoys.push(t_decoy.unwrap_or(-1));
            }
            if want("TRANSITION_TYPE") {
                target.transition_types.push(t_type);
            }
            if want("ANNOTATION") {
                target.annotations.push(t_annotation);
            }
        }
    }
    Ok(output)
}

/// `getRuns` of the XIC and XIM readers: distinct `(RUN_ID, SOURCE_FILE)`
/// pairs in first-seen order.
///
/// # Errors
///
/// [`Error::MissingInformation`] when the files have neither column or no
/// `RUN_ID`; [`Error::InvalidValue`] "NULL value in required RUN_ID column".
pub(crate) fn runs(filenames: &[PathBuf], limits: &ReaderLimits) -> Result<Vec<RunInfo>> {
    let table = read_tables(filenames, Some(&["RUN_ID", "SOURCE_FILE"]), limits)?;
    let run_id = column(&table, "RUN_ID", true)?;
    let source_file = column(&table, "SOURCE_FILE", false)?;
    let mut seen = HashSet::new();
    let mut output = Vec::new();
    for row in 0..table.num_rows() {
        let info = RunInfo {
            run_id: optional_int(run_id.as_ref(), row)?.ok_or_else(|| {
                Error::InvalidValue(format!("NULL value in required RUN_ID column: {row}"))
            })?,
            source_file: string_or_empty(source_file.as_ref(), row)?,
        };
        let mut key = String::new();
        key_int(&mut key, Some(info.run_id));
        key_str(&mut key, &info.source_file);
        if seen.insert(key) {
            output.push(info);
        }
    }
    Ok(output)
}

/// `getColumns`: the column names of the files' shared schema.
///
/// # Errors
///
/// As [`read_schema_all_files`].
pub(crate) fn column_names(filenames: &[PathBuf]) -> Result<Vec<String>> {
    Ok(read_schema_all_files(filenames)?
        .fields()
        .iter()
        .map(|field| field.name().clone())
        .collect())
}

/// Collect paths into owned `PathBuf`s.
pub(crate) fn owned_paths<P: AsRef<Path>>(filenames: &[P]) -> Vec<PathBuf> {
    filenames
        .iter()
        .map(|path| path.as_ref().to_path_buf())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokenizer_and_parser_follow_the_source_grammar() {
        let types = [("ANNOTATION", ColumnType::String)];
        let expr = parse_filter(
            "precursor_id==2 || annotation = 'y3^1' AND transition_id in [1, 2]",
            &types,
        )
        .expect("valid filter");
        assert_eq!(expr.conditions.len(), 3);
        assert_eq!(expr.conditions[0].column, "PRECURSOR_ID");
        assert_eq!(expr.conditions[0].op, "=");
        assert_eq!(expr.conditions[1].values, ["y3^1"]);
        assert_eq!(expr.conditions[1].column_type, ColumnType::String);
        assert_eq!(expr.conditions[2].values, ["1", "2"]);
        assert_eq!(expr.connectors, ["OR", "AND"]);
        // Two conditions without a connector are joined by nothing: the
        // evaluator then reads AND.
        let bare = parse_filter("a=1 b=2", &[]).expect("valid filter");
        assert_eq!(bare.conditions.len(), 2);
        assert!(bare.connectors.is_empty());
        for bad in [
            "a", "a ~ 1", "a in 1", "a in []", "a in [1", "a =", "a=1 and", "a='x",
        ] {
            assert!(parse_filter(bad, &[]).is_err(), "{bad}");
        }
    }
}
