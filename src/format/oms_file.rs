// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! OMS files: SQLite databases of identification and feature data
//! (`FORMAT/OMSFile.h`, with the internal helpers `FORMAT/OMSFileLoad.h` and
//! `FORMAT/OMSFileStore.h`).
//!
//! An `.oms` file is an SQLite database whose tables mirror OpenMS's
//! `IdentificationData` (`ID_*` tables), feature and consensus maps (`FEAT_*`),
//! plus shared tables for controlled-vocabulary terms (`CVTerm`), metadata
//! values (`DataValue`, `DataValue_DataType`) and a `version` table whose
//! `OMSFile` column is the schema version.
//!
//! This module currently ports the JSON export of `OMSFileLoad`
//! ([`OmsFile::export_to_json`]), which dumps every table of any `.oms` file
//! without interpreting it. Loading into and storing from `IdentificationData`,
//! `FeatureMap` and `ConsensusMap` are recorded as open in
//! `docs/OMS_FILE_SUPPORT.md`.
//!
//! The module needs the optional `sqlite` feature, which brings in the bundled
//! SQLite of `rusqlite`.

use crate::concept::progress_logger::ProgressLogType;
use crate::format::path_io;
use crate::{Error, Result};
use rusqlite::types::ValueRef;
use rusqlite::{Connection, OpenFlags};
use std::path::Path;

/// The schema version the pinned source writes (`OMSFileStore.cpp`
/// `version_number`); readers accept any version the `version` table states.
pub const SCHEMA_VERSION: i64 = 5;

/// Largest JSON text [`OmsFile::export_to_json`] produces, in bytes.
///
/// The source builds the whole JSON document in memory without a ceiling
/// (its own `@TODO` notes it); this port refuses beyond the limit before the
/// output file is created.
pub const MAX_JSON_BYTES: usize = 1 << 30;

/// Largest number of tables an exported database may hold.
pub const MAX_TABLES: usize = 10_000;

fn sql_error(context: &str, error: rusqlite::Error) -> Error {
    // Source: raiseDBError_ throws Exception::FailedAPICall with
    // "<context>: <SQLite message>".
    Error::Io(std::io::Error::other(format!("{context}: {error}")))
}

fn limit() -> Error {
    Error::InvalidValue("OMS JSON export exceeds its size limit".into())
}

/// Row order of the tables that do not order by `id` (the source's
/// `OMSFileLoad::export_order_by_`); `version` is exported unordered.
fn order_by(table: &str) -> Option<&'static str> {
    if table.ends_with("_MetaInfo") {
        return Some("parent_id, name");
    }
    if table.ends_with("_AppliedProcessingStep") {
        return Some("parent_id, processing_step_order, score_type_id");
    }
    match table {
        "version" => None,
        "ID_IdentifiedCompound" => Some("molecule_id"),
        "ID_ParentMatch" => Some("molecule_id, parent_id, start_pos, end_pos"),
        "ID_ParentGroup_ParentSequence" => Some("group_id, parent_id"),
        "ID_ProcessingStep_InputFile" => Some("processing_step_id, input_file_id"),
        "ID_ProcessingSoftware_AssignedScore" => Some("software_id, score_type_order"),
        "ID_ObservationMatch_PeakAnnotation" => {
            Some("parent_id, processing_step_id, peak_mz, peak_annotation")
        }
        "FEAT_ConvexHull" => Some("feature_id, hull_index, point_index"),
        "FEAT_ObservationMatch" => Some("feature_id, observation_match_id"),
        "FEAT_MapMetaData" => Some("unique_id"),
        _ => Some("id"),
    }
}

/// An SQL identifier in double quotes, with embedded quotes doubled.
fn quoted(identifier: &str) -> String {
    format!("\"{}\"", identifier.replace('"', "\"\""))
}

/// Reader and writer of OMS files. Source `OMSFile`.
///
/// The source derives from `ProgressLogger` and takes its log type in the
/// constructor; the log type is kept here for the load and store operations.
/// The JSON export reports no progress, in the source as here.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct OmsFile {
    log_type: ProgressLogType,
}

impl OmsFile {
    /// An adapter with the given progress log type. Source
    /// `OMSFile(LogType log_type = LogType::NONE)`.
    pub fn new(log_type: ProgressLogType) -> Self {
        Self { log_type }
    }

    /// The progress log type.
    pub fn log_type(&self) -> ProgressLogType {
        self.log_type
    }

    /// Open `filename` read-only and read its schema version, as the source's
    /// `OMSFileLoad` constructor does.
    fn open(filename: &Path) -> Result<(Connection, i64)> {
        let connection = Connection::open_with_flags(
            filename,
            OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
        )
        .map_err(|error| sql_error("error opening OMS file", error))?;
        let version = connection
            .query_row("SELECT OMSFile FROM version", [], |row| row.get::<_, i64>(0))
            .map_err(|error| sql_error("error reading file format version number", error))?;
        Ok((connection, version))
    }

    /// The schema version an OMS file declares in its `version` table.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be opened as an SQLite database or
    /// has no readable `version` table (the source's
    /// `Exception::FailedAPICall`, "error reading file format version number").
    pub fn schema_version(&self, filename: impl AsRef<Path>) -> Result<i64> {
        Ok(Self::open(filename.as_ref())?.1)
    }

    /// Render the contents of an OMS file as JSON text.
    ///
    /// Source `OMSFileLoad::exportToJSON`: one object member per table (every
    /// table except SQLite's own), whose value is an array with one object per
    /// row; a row object maps column names to values. Integers and reals keep
    /// their SQLite storage class, text and blobs become strings, and `NULL`
    /// becomes the empty string. Rows are ordered by `id`, by
    /// `parent_id, name` for `*_MetaInfo` tables, by
    /// `parent_id, processing_step_order, score_type_id` for
    /// `*_AppliedProcessingStep` tables and by the source's fixed column lists
    /// for eight tables without an `id` column; `version` is unordered.
    ///
    /// The text is formatted as `nlohmann::json::dump(4)` formats it: members
    /// sorted by name, four-space indentation, `": "` between name and value,
    /// shortest round-trip numbers with a two-digit exponent and a `.0` on
    /// integral reals, and a final newline.
    ///
    /// # Errors
    ///
    /// As [`OmsFile::schema_version`]; [`Error::Io`] when a table cannot be
    /// read; [`Error::InvalidValue`] when a text or blob cell is not UTF-8 (the
    /// source's JSON library throws `type_error.316` there) or the text would
    /// exceed [`MAX_JSON_BYTES`], or the database holds more than
    /// [`MAX_TABLES`] tables.
    pub fn json_text(&self, filename: impl AsRef<Path>) -> Result<String> {
        let (connection, _) = Self::open(filename.as_ref())?;
        let mut tables: Vec<String> = Vec::new();
        {
            let mut query = connection
                .prepare(
                    "SELECT name FROM sqlite_master WHERE type='table' \
                     AND name NOT LIKE 'sqlite_%' ORDER BY name",
                )
                .map_err(|error| sql_error("error listing tables", error))?;
            let mut rows = query
                .query([])
                .map_err(|error| sql_error("error listing tables", error))?;
            while let Some(row) = rows
                .next()
                .map_err(|error| sql_error("error listing tables", error))?
            {
                if tables.len() >= MAX_TABLES {
                    return Err(Error::InvalidValue(
                        "OMS file holds more tables than the export accepts".into(),
                    ));
                }
                tables.push(
                    row.get::<_, String>(0)
                        .map_err(|error| sql_error("error listing tables", error))?,
                );
            }
        }
        // The source's JSON object is a std::map, which orders members
        // bytewise; SQLite's ORDER BY name uses the same BINARY collation.
        tables.sort();
        let mut out = JsonWriter::default();
        out.raw("{")?;
        for (table_index, table) in tables.iter().enumerate() {
            out.raw(if table_index == 0 { "\n" } else { ",\n" })?;
            out.indent(1)?;
            out.string(table)?;
            out.raw(": ")?;
            export_table(&connection, table, &mut out)?;
        }
        if !tables.is_empty() {
            out.raw("\n")?;
        }
        out.raw("}\n")?;
        Ok(out.text)
    }

    /// Write the contents of the OMS file `filename_in` as JSON to
    /// `filename_out`. Source `OMSFile::exportToJSON`.
    ///
    /// The text is [`OmsFile::json_text`]; it is rendered completely before the
    /// output is created and written atomically.
    ///
    /// # Errors
    ///
    /// As [`OmsFile::json_text`], plus [`Error::Io`] when the output cannot be
    /// written (the source's `Exception::FileNotWritable`).
    pub fn export_to_json(
        &self,
        filename_in: impl AsRef<Path>,
        filename_out: impl AsRef<Path>,
    ) -> Result<()> {
        let text = self.json_text(filename_in)?;
        path_io::write(filename_out.as_ref(), |writer| {
            writer.write_all(text.as_bytes())?;
            Ok(())
        })
    }
}

/// Source `exportTableToJSON_`, appending one table's array to `out`.
fn export_table(connection: &Connection, table: &str, out: &mut JsonWriter) -> Result<()> {
    let context = format!("error exporting table {table}");
    let mut sql = format!("SELECT * FROM {}", quoted(table));
    if let Some(order) = order_by(table) {
        sql.push_str(" ORDER BY ");
        sql.push_str(order);
    }
    let mut query = connection
        .prepare(&sql)
        .map_err(|error| sql_error(&context, error))?;
    // The member order of a row object is the sorted column names; a repeated
    // name keeps its last value, as assignment into the source's std::map does.
    let names: Vec<String> = query.column_names().into_iter().map(String::from).collect();
    let mut order: Vec<usize> = (0..names.len()).collect();
    order.sort_by(|a, b| names[*a].cmp(&names[*b]).then(b.cmp(a)));
    order.dedup_by(|later, earlier| names[*later] == names[*earlier]);
    let mut rows = query
        .query([])
        .map_err(|error| sql_error(&context, error))?;
    let mut count = 0usize;
    out.raw("[")?;
    while let Some(row) = rows.next().map_err(|error| sql_error(&context, error))? {
        out.raw(if count == 0 { "\n" } else { ",\n" })?;
        count += 1;
        out.indent(2)?;
        out.raw("{")?;
        for (position, column) in order.iter().enumerate() {
            out.raw(if position == 0 { "\n" } else { ",\n" })?;
            out.indent(3)?;
            out.string(&names[*column])?;
            out.raw(": ")?;
            match row
                .get_ref(*column)
                .map_err(|error| sql_error(&context, error))?
            {
                ValueRef::Null => out.string("")?,
                ValueRef::Integer(value) => out.raw(&value.to_string())?,
                ValueRef::Real(value) => out.raw(&json_number(value))?,
                ValueRef::Text(bytes) | ValueRef::Blob(bytes) => {
                    let text = std::str::from_utf8(bytes).map_err(|_| {
                        Error::InvalidValue(format!(
                            "OMS table {table} holds a cell that is not UTF-8 text"
                        ))
                    })?;
                    out.string(text)?;
                }
            }
        }
        if !order.is_empty() {
            out.raw("\n")?;
            out.indent(2)?;
        }
        out.raw("}")?;
    }
    if count > 0 {
        out.raw("\n")?;
        out.indent(1)?;
    }
    out.raw("]")
}

/// A number as `nlohmann::json` prints a double: the shortest digits that
/// round-trip, in fixed notation for decimal exponents in `(-4, 15]` (with
/// `.0` appended to an integral value) and in scientific notation otherwise,
/// the exponent signed and at least two digits long. Non-finite values print
/// `null`, as the library does.
fn json_number(value: f64) -> String {
    if !value.is_finite() {
        return "null".into();
    }
    if value == 0.0 {
        return if value.is_sign_negative() { "-0.0" } else { "0.0" }.into();
    }
    // Rust's `{:e}` is the shortest round-trip representation.
    let scientific = format!("{:e}", value.abs());
    let (mantissa, exponent) = scientific.split_once('e').unwrap_or((&scientific, "0"));
    let exponent: i32 = exponent.parse().unwrap_or(0);
    let digits: String = mantissa.chars().filter(char::is_ascii_digit).collect();
    let k = digits.len() as i32;
    // `n` is the position of the decimal point relative to the digit string.
    let n = exponent + 1;
    let sign = if value.is_sign_negative() { "-" } else { "" };
    let body = if k <= n && n <= 15 {
        format!("{digits}{}.0", "0".repeat((n - k) as usize))
    } else if 0 < n && n <= 15 {
        format!("{}.{}", &digits[..n as usize], &digits[n as usize..])
    } else if -4 < n && n <= 0 {
        format!("0.{}{digits}", "0".repeat((-n) as usize))
    } else {
        let exponent = n - 1;
        let exponent_text = format!(
            "{}{:02}",
            if exponent < 0 { '-' } else { '+' },
            exponent.abs()
        );
        if k == 1 {
            format!("{digits}e{exponent_text}")
        } else {
            format!("{}.{}e{exponent_text}", &digits[..1], &digits[1..])
        }
    };
    format!("{sign}{body}")
}

/// A bounded JSON text accumulator.
#[derive(Default)]
struct JsonWriter {
    text: String,
}

impl JsonWriter {
    fn raw(&mut self, text: &str) -> Result<()> {
        if self.text.len().saturating_add(text.len()) > MAX_JSON_BYTES {
            return Err(limit());
        }
        self.text.push_str(text);
        Ok(())
    }

    fn indent(&mut self, level: usize) -> Result<()> {
        self.raw(&" ".repeat(4 * level))
    }

    /// A JSON string as `nlohmann::json::dump` escapes it: `"` and `\`, the
    /// short escapes for backspace, form feed, newline, carriage return and
    /// tab, `\u00xx` (lower-case hex) for the other control characters, and
    /// every other character verbatim.
    fn string(&mut self, text: &str) -> Result<()> {
        let mut escaped = String::with_capacity(text.len() + 2);
        escaped.push('"');
        for c in text.chars() {
            match c {
                '"' => escaped.push_str("\\\""),
                '\\' => escaped.push_str("\\\\"),
                '\u{08}' => escaped.push_str("\\b"),
                '\u{0c}' => escaped.push_str("\\f"),
                '\n' => escaped.push_str("\\n"),
                '\r' => escaped.push_str("\\r"),
                '\t' => escaped.push_str("\\t"),
                c if (c as u32) < 0x20 => escaped.push_str(&format!("\\u{:04x}", c as u32)),
                c => escaped.push(c),
            }
        }
        escaped.push('"');
        self.raw(&escaped)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numbers_follow_the_nlohmann_format() {
        assert_eq!(json_number(0.3), "0.3");
        assert_eq!(json_number(1.0), "1.0");
        assert_eq!(json_number(-2.5), "-2.5");
        assert_eq!(json_number(1234.5678), "1234.5678");
        assert_eq!(json_number(0.0001), "0.0001");
        assert_eq!(json_number(0.00001), "1e-05");
        assert_eq!(json_number(1.5e-7), "1.5e-07");
        assert_eq!(json_number(1e14), "100000000000000.0");
        assert_eq!(json_number(1e15), "1e+15");
        assert_eq!(json_number(1.25e16), "1.25e+16");
        assert_eq!(json_number(123_456_789_012_345_680.0), "1.2345678901234568e+17");
        assert_eq!(json_number(f64::NAN), "null");
        assert_eq!(json_number(0.0), "0.0");
    }

    #[test]
    fn strings_escape_like_nlohmann() {
        let mut out = JsonWriter::default();
        out.string("a\"b\\c\nd\u{01}é/").unwrap();
        assert_eq!(out.text, "\"a\\\"b\\\\c\\nd\\u0001é/\"");
    }

    #[test]
    fn order_by_matches_the_source_table() {
        assert_eq!(order_by("version"), None);
        assert_eq!(order_by("ID_ScoreType"), Some("id"));
        assert_eq!(order_by("ID_Observation_MetaInfo"), Some("parent_id, name"));
        assert_eq!(
            order_by("ID_ObservationMatch_AppliedProcessingStep"),
            Some("parent_id, processing_step_order, score_type_id")
        );
        assert_eq!(order_by("FEAT_MapMetaData"), Some("unique_id"));
    }
}
