// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Writing of the OpenSWATH OSW (SQLite) output, the input of PyProphet.
//!
//! Ports `ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h` with its `.cpp`.
//! [`OpenSwathOSWWriter`] turns the features of a
//! [`FeatureMap`](crate::kernel::FeatureMap) into buffered table rows
//! ([`OswData`] of typed [`OswValue`] cells, [`OpenSwathOSWWriter::prepare_rows`])
//! and inserts them into the `RUN`, `FEATURE`, `FEATURE_MS1`, `FEATURE_MS2`,
//! `FEATURE_PRECURSOR` and `FEATURE_TRANSITION` tables of the file
//! ([`OpenSwathOSWWriter::write_rows`]). The tables are those of the PQP
//! transition library format; see the support document.
//!
//! The row preparation is pure Rust and always present. The members that
//! touch the database ([`write_header`](OpenSwathOSWWriter::write_header),
//! [`add_run`](OpenSwathOSWWriter::add_run),
//! [`write_rows`](OpenSwathOSWWriter::write_rows),
//! [`write_lines`](OpenSwathOSWWriter::write_lines)) need the `sqlite` feature
//! and return [`Error::Unsupported`](crate::Error::Unsupported) without it, so
//! the workflow that drives the writer compiles either way.
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_OSW_WRITER_SUPPORT.md).
//!
//! # Example
//!
//! ```
//! use openms::analysis::openswath::openswath_osw_writer::OswValue;
//!
//! // text cells are classified like the source: integer, then real, else text
//! assert_eq!(OswValue::parse("123"), OswValue::Int64(123));
//! assert_eq!(OswValue::parse("0.5"), OswValue::Double(0.5));
//! assert_eq!(OswValue::parse("NaN"), OswValue::Null);
//! assert_eq!(OswValue::parse("NIAGGLEK_2"), OswValue::Text("NIAGGLEK_2".into()));
//! ```

use super::mrm_feature_support::data_value_to_string;
use super::transition_experiment::{LightCompound, LightTransition};
use crate::format::file_info::text_format::to_str;
#[cfg(feature = "sqlite")]
use crate::format::sqlite_connector::SqliteConnector;
use crate::kernel::{Feature, FeatureMap};
use crate::metadata::{MetaValue, MetaValueData};
use crate::{Error, Result};
#[cfg(feature = "sqlite")]
use std::sync::Mutex;

/// Number of columns of a `FEATURE` row.
pub const FEATURE_COLUMN_COUNT: usize = 11;
/// Number of columns of a `FEATURE_MS1` row.
pub const FEATURE_MS1_COLUMN_COUNT: usize = 18;
/// Number of columns of a `FEATURE_MS2` row.
pub const FEATURE_MS2_COLUMN_COUNT: usize = 38;
/// Number of columns of a `FEATURE_PRECURSOR` row.
pub const FEATURE_PRECURSOR_COLUMN_COUNT: usize = 4;
/// Number of columns of a `FEATURE_TRANSITION` row.
pub const FEATURE_TRANSITION_COLUMN_COUNT: usize = 43;

/// A `FEATURE` row.
pub type FeatureRow = [OswValue; FEATURE_COLUMN_COUNT];
/// A `FEATURE_MS1` row.
pub type FeatureMS1Row = [OswValue; FEATURE_MS1_COLUMN_COUNT];
/// A `FEATURE_MS2` row.
pub type FeatureMS2Row = [OswValue; FEATURE_MS2_COLUMN_COUNT];
/// A `FEATURE_PRECURSOR` row.
pub type FeaturePrecursorRow = [OswValue; FEATURE_PRECURSOR_COLUMN_COUNT];
/// A `FEATURE_TRANSITION` row.
pub type FeatureTransitionRow = [OswValue; FEATURE_TRANSITION_COLUMN_COUNT];

/// Column names of the `FEATURE` table, in insertion order.
pub const FEATURE_COLUMNS: [&str; FEATURE_COLUMN_COUNT] = [
    "ID",
    "RUN_ID",
    "PRECURSOR_ID",
    "EXP_RT",
    "EXP_IM",
    "NORM_RT",
    "DELTA_RT",
    "LEFT_WIDTH",
    "RIGHT_WIDTH",
    "EXP_IM_LEFTWIDTH",
    "EXP_IM_RIGHTWIDTH",
];

/// Column names of the `FEATURE_MS1` table, in insertion order.
pub const FEATURE_MS1_COLUMNS: [&str; FEATURE_MS1_COLUMN_COUNT] = [
    "FEATURE_ID",
    "AREA_INTENSITY",
    "APEX_INTENSITY",
    "EXP_IM",
    "DELTA_IM",
    "VAR_MASSDEV_SCORE",
    "VAR_IM_MS1_DELTA_SCORE",
    "VAR_MI_SCORE",
    "VAR_MI_CONTRAST_SCORE",
    "VAR_MI_COMBINED_SCORE",
    "VAR_ISOTOPE_CORRELATION_SCORE",
    "VAR_ISOTOPE_OVERLAP_SCORE",
    "VAR_XCORR_COELUTION",
    "VAR_XCORR_COELUTION_CONTRAST",
    "VAR_XCORR_COELUTION_COMBINED",
    "VAR_XCORR_SHAPE",
    "VAR_XCORR_SHAPE_CONTRAST",
    "VAR_XCORR_SHAPE_COMBINED",
];

/// Column names of the `FEATURE_MS2` table, in insertion order.
pub const FEATURE_MS2_COLUMNS: [&str; FEATURE_MS2_COLUMN_COUNT] = [
    "FEATURE_ID",
    "AREA_INTENSITY",
    "TOTAL_AREA_INTENSITY",
    "APEX_INTENSITY",
    "EXP_IM",
    "EXP_IM_LEFTWIDTH",
    "EXP_IM_RIGHTWIDTH",
    "DELTA_IM",
    "TOTAL_MI",
    "VAR_BSERIES_SCORE",
    "VAR_DOTPROD_SCORE",
    "VAR_INTENSITY_SCORE",
    "VAR_ISOTOPE_CORRELATION_SCORE",
    "VAR_ISOTOPE_OVERLAP_SCORE",
    "VAR_LIBRARY_CORR",
    "VAR_LIBRARY_DOTPROD",
    "VAR_LIBRARY_MANHATTAN",
    "VAR_LIBRARY_RMSD",
    "VAR_LIBRARY_ROOTMEANSQUARE",
    "VAR_LIBRARY_SANGLE",
    "VAR_LOG_SN_SCORE",
    "VAR_MANHATTAN_SCORE",
    "VAR_MASSDEV_SCORE",
    "VAR_MASSDEV_SCORE_WEIGHTED",
    "VAR_MI_SCORE",
    "VAR_MI_WEIGHTED_SCORE",
    "VAR_MI_RATIO_SCORE",
    "VAR_NORM_RT_SCORE",
    "VAR_XCORR_COELUTION",
    "VAR_XCORR_COELUTION_WEIGHTED",
    "VAR_XCORR_SHAPE",
    "VAR_XCORR_SHAPE_WEIGHTED",
    "VAR_YSERIES_SCORE",
    "VAR_ELUTION_MODEL_FIT_SCORE",
    "VAR_IM_XCORR_SHAPE",
    "VAR_IM_XCORR_COELUTION",
    "VAR_IM_DELTA_SCORE",
    "VAR_IM_LOG_INTENSITY",
];

/// Column names of the `FEATURE_PRECURSOR` table, in insertion order.
pub const FEATURE_PRECURSOR_COLUMNS: [&str; FEATURE_PRECURSOR_COLUMN_COUNT] =
    ["FEATURE_ID", "ISOTOPE", "AREA_INTENSITY", "APEX_INTENSITY"];

/// Column names of the `FEATURE_TRANSITION` table, in insertion order.
pub const FEATURE_TRANSITION_COLUMNS: [&str; FEATURE_TRANSITION_COLUMN_COUNT] = [
    "FEATURE_ID",
    "TRANSITION_ID",
    "AREA_INTENSITY",
    "TOTAL_AREA_INTENSITY",
    "APEX_RT",
    "APEX_INTENSITY",
    "RT_FWHM",
    "MASSERROR_PPM",
    "TOTAL_MI",
    "VAR_INTENSITY_SCORE",
    "VAR_INTENSITY_RATIO_SCORE",
    "VAR_LOG_INTENSITY",
    "VAR_XCORR_COELUTION",
    "VAR_XCORR_SHAPE",
    "VAR_LOG_SN_SCORE",
    "VAR_MASSDEV_SCORE",
    "VAR_MI_SCORE",
    "VAR_MI_RATIO_SCORE",
    "VAR_ISOTOPE_CORRELATION_SCORE",
    "VAR_ISOTOPE_OVERLAP_SCORE",
    "EXP_IM",
    "EXP_IM_LEFTWIDTH",
    "EXP_IM_RIGHTWIDTH",
    "DELTA_IM",
    "VAR_IM_DELTA_SCORE",
    "VAR_IM_LOG_INTENSITY",
    "VAR_IM_XCORR_COELUTION_CONTRAST",
    "VAR_IM_XCORR_SHAPE_CONTRAST",
    "VAR_IM_XCORR_COELUTION_COMBINED",
    "VAR_IM_XCORR_SHAPE_COMBINED",
    "START_POSITION_AT_5",
    "END_POSITION_AT_5",
    "START_POSITION_AT_10",
    "END_POSITION_AT_10",
    "START_POSITION_AT_50",
    "END_POSITION_AT_50",
    "TOTAL_WIDTH",
    "TAILING_FACTOR",
    "ASYMMETRY_FACTOR",
    "SLOPE_OF_BASELINE",
    "BASELINE_DELTA_2_HEIGHT",
    "POINTS_ACROSS_BASELINE",
    "POINTS_ACROSS_HALF_HEIGHT",
];

/// The `CREATE TABLE` statements of [`OpenSwathOSWWriter::write_header`],
/// verbatim from the source (the `FEATURE_TRANSITION` columns
/// `VAR_IM_XCORR_*` carry no declared type there either).
pub const CREATE_SQL: &str = "CREATE TABLE RUN(ID INT PRIMARY KEY NOT NULL,FILENAME TEXT NOT NULL); CREATE TABLE FEATURE(ID INT PRIMARY KEY NOT NULL,RUN_ID INT NOT NULL,PRECURSOR_ID INT NOT NULL,EXP_RT REAL NOT NULL,EXP_IM REAL, NORM_RT REAL NOT NULL,DELTA_RT REAL NOT NULL,LEFT_WIDTH REAL NOT NULL,RIGHT_WIDTH REAL NOT NULL,EXP_IM_LEFTWIDTH REAL,EXP_IM_RIGHTWIDTH REAL); CREATE TABLE FEATURE_MS1(FEATURE_ID INT NOT NULL,AREA_INTENSITY REAL NOT NULL,APEX_INTENSITY REAL NOT NULL,EXP_IM REAL,DELTA_IM REAL,VAR_MASSDEV_SCORE REAL NULL,VAR_MI_SCORE REAL NULL,VAR_MI_CONTRAST_SCORE REAL NULL,VAR_MI_COMBINED_SCORE REAL NULL,VAR_ISOTOPE_CORRELATION_SCORE REAL NULL,VAR_ISOTOPE_OVERLAP_SCORE REAL NULL,VAR_IM_MS1_DELTA_SCORE REAL NULL,VAR_XCORR_COELUTION REAL NULL,VAR_XCORR_COELUTION_CONTRAST REAL NULL,VAR_XCORR_COELUTION_COMBINED REAL NULL,VAR_XCORR_SHAPE REAL NULL,VAR_XCORR_SHAPE_CONTRAST REAL NULL,VAR_XCORR_SHAPE_COMBINED REAL NULL); CREATE TABLE FEATURE_MS2(FEATURE_ID INT NOT NULL,AREA_INTENSITY REAL NOT NULL,TOTAL_AREA_INTENSITY REAL NOT NULL,APEX_INTENSITY REAL NOT NULL,EXP_IM REAL,EXP_IM_LEFTWIDTH REAL,EXP_IM_RIGHTWIDTH REAL,DELTA_IM REAL,TOTAL_MI REAL NULL,VAR_BSERIES_SCORE REAL NULL,VAR_DOTPROD_SCORE REAL NULL,VAR_INTENSITY_SCORE REAL NULL,VAR_ISOTOPE_CORRELATION_SCORE REAL NULL,VAR_ISOTOPE_OVERLAP_SCORE REAL NULL,VAR_LIBRARY_CORR REAL NULL,VAR_LIBRARY_DOTPROD REAL NULL,VAR_LIBRARY_MANHATTAN REAL NULL,VAR_LIBRARY_RMSD REAL NULL,VAR_LIBRARY_ROOTMEANSQUARE REAL NULL,VAR_LIBRARY_SANGLE REAL NULL,VAR_LOG_SN_SCORE REAL NULL,VAR_MANHATTAN_SCORE REAL NULL,VAR_MASSDEV_SCORE REAL NULL,VAR_MASSDEV_SCORE_WEIGHTED REAL NULL,VAR_MI_SCORE REAL NULL,VAR_MI_WEIGHTED_SCORE REAL NULL,VAR_MI_RATIO_SCORE REAL NULL,VAR_NORM_RT_SCORE REAL NULL,VAR_XCORR_COELUTION REAL NULL,VAR_XCORR_COELUTION_WEIGHTED REAL NULL,VAR_XCORR_SHAPE REAL NULL,VAR_XCORR_SHAPE_WEIGHTED REAL NULL,VAR_YSERIES_SCORE REAL NULL,VAR_ELUTION_MODEL_FIT_SCORE REAL NULL,VAR_IM_XCORR_SHAPE REAL NULL,VAR_IM_XCORR_COELUTION REAL NULL,VAR_IM_DELTA_SCORE REAL NULL,VAR_IM_LOG_INTENSITY REAL NULL);CREATE TABLE FEATURE_PRECURSOR(FEATURE_ID INT NOT NULL,ISOTOPE INT NOT NULL,AREA_INTENSITY REAL NOT NULL,APEX_INTENSITY REAL NOT NULL);CREATE TABLE FEATURE_TRANSITION(FEATURE_ID INT NOT NULL,TRANSITION_ID INT NOT NULL,AREA_INTENSITY REAL NOT NULL,TOTAL_AREA_INTENSITY REAL NOT NULL,APEX_RT REAL NULL,APEX_INTENSITY REAL NOT NULL,RT_FWHM REAL NOT NULL,MASSERROR_PPM REAL NULL,TOTAL_MI REAL NULL,VAR_INTENSITY_SCORE REAL NULL,VAR_INTENSITY_RATIO_SCORE REAL NULL,VAR_LOG_INTENSITY REAL NULL,VAR_XCORR_COELUTION REAL NULL,VAR_XCORR_SHAPE REAL NULL,VAR_LOG_SN_SCORE REAL NULL,VAR_MASSDEV_SCORE REAL NULL,VAR_MI_SCORE REAL NULL,VAR_MI_RATIO_SCORE REAL NULL,VAR_ISOTOPE_CORRELATION_SCORE REAL NULL,VAR_ISOTOPE_OVERLAP_SCORE REAL NULL, EXP_IM REAL NULL,EXP_IM_LEFTWIDTH REAL,EXP_IM_RIGHTWIDTH REAL,DELTA_IM REAL NULL,VAR_IM_DELTA_SCORE REAL NULL,VAR_IM_LOG_INTENSITY REAL NULL,VAR_IM_XCORR_COELUTION_CONTRAST, VAR_IM_XCORR_SHAPE_CONTRAST, VAR_IM_XCORR_COELUTION_COMBINED, VAR_IM_XCORR_SHAPE_COMBINED, START_POSITION_AT_5 REAL NULL, END_POSITION_AT_5 REAL NULL, START_POSITION_AT_10 REAL NULL, END_POSITION_AT_10 REAL NULL, START_POSITION_AT_50 REAL NULL, END_POSITION_AT_50 REAL NULL, TOTAL_WIDTH REAL NULL, TAILING_FACTOR REAL NULL, ASYMMETRY_FACTOR REAL NULL, SLOPE_OF_BASELINE REAL NULL, BASELINE_DELTA_2_HEIGHT REAL NULL, POINTS_ACROSS_BASELINE REAL NULL, POINTS_ACROSS_HALF_HEIGHT REAL NULL); ";

/// The unique ID with its top bit cleared, as a non-negative `Int64` (source
/// `SqliteHelper::clearSignBit`).
pub fn clear_sign_bit(id: u64) -> i64 {
    (id & !(1_u64 << 63)) as i64
}

/// A typed OSW table cell, bound directly to SQLite (source
/// `OpenSwathOSWWriter::OSWValue`).
///
/// A [`Double`](Self::Double) holding NaN counts as NULL ([`is_null`](Self::is_null)),
/// as in the source, whose double constructor stores a NaN as a NULL cell.
#[derive(Clone, Debug, Default, PartialEq)]
pub enum OswValue {
    /// SQL NULL.
    #[default]
    Null,
    /// An integer.
    Int64(i64),
    /// A real number.
    Double(f64),
    /// Text, stored without any numeric or NULL interpretation.
    Text(String),
}

impl OswValue {
    /// A NULL cell (source `OSWValue::null`).
    pub const fn null() -> Self {
        Self::Null
    }

    /// A text cell without numeric or NULL parsing (source `OSWValue::text`).
    pub fn text(value: impl Into<String>) -> Self {
        Self::Text(value.into())
    }

    /// A real cell; NaN becomes NULL (source `OSWValue(double)`).
    pub fn real(value: f64) -> Self {
        if value.is_nan() {
            Self::Null
        } else {
            Self::Double(value)
        }
    }

    /// Classify a text the way the source's string constructor does.
    ///
    /// The empty string and `NULL`, `nan`, `-nan`, `NaN`, `-NaN` give NULL;
    /// text that `strtoll` consumes completely (base 10, optional sign and
    /// surrounding white space) an integer; text that `strtod` consumes
    /// completely a real, hexadecimal and `inf`/`infinity` forms included (a NaN,
    /// also `nan(...)`, again NULL); any other text, and any number that over-
    /// or underflows, a text cell.
    ///
    /// Hexadecimal significands of more than 16 digits are rounded twice
    /// here, where `strtod` rounds once; no such text is expected in an
    /// identifier column.
    pub fn parse(value: &str) -> Self {
        if matches!(value, "" | "NULL" | "nan" | "-nan" | "NaN" | "-NaN") {
            return Self::Null;
        }
        let c_space = |c: char| matches!(c, ' ' | '\t' | '\n' | '\x0b' | '\x0c' | '\r');
        let trimmed = value.trim_matches(c_space);
        if let Ok(int) = trimmed.parse::<i64>() {
            return Self::Int64(int);
        }
        let Some(real) = parse_c_double(trimmed) else {
            return Self::Text(value.to_owned());
        };
        if real.is_nan() {
            return Self::Null;
        }
        let lower = trimmed.to_ascii_lowercase();
        let literal_inf = lower.trim_start_matches(['+', '-']).starts_with("inf");
        if real.is_infinite() && !literal_inf {
            return Self::Text(value.to_owned()); // ERANGE
        }
        if real != 0.0 && real.abs() < f64::MIN_POSITIVE {
            return Self::Text(value.to_owned()); // ERANGE (underflow)
        }
        if real == 0.0 && has_nonzero_significand(&lower) {
            return Self::Text(value.to_owned()); // ERANGE (underflow to zero)
        }
        Self::Double(real)
    }

    /// Whether the cell is inserted as SQL NULL: a NULL cell or a NaN.
    pub fn is_null(&self) -> bool {
        match self {
            Self::Null => true,
            Self::Double(d) => d.is_nan(),
            _ => false,
        }
    }

    /// The integer payload, if this is an integer cell (source `asInt`, which
    /// is undefined for another type).
    pub fn as_int(&self) -> Option<i64> {
        match self {
            Self::Int64(i) => Some(*i),
            _ => None,
        }
    }

    /// The real payload, if this is a real cell (source `asDouble`).
    pub fn as_double(&self) -> Option<f64> {
        match self {
            Self::Double(d) => Some(*d),
            _ => None,
        }
    }

    /// The text payload, if this is a text cell (source `asText`).
    pub fn as_text(&self) -> Option<&str> {
        match self {
            Self::Text(t) => Some(t),
            _ => None,
        }
    }

    /// The text of the cell as the source's `oswValueToString` writes it:
    /// `NULL` for a NULL or NaN cell, an integer in decimal, a real through
    /// `StringUtils::toStr(double)`, text verbatim.
    pub fn to_text(&self) -> String {
        match self {
            Self::Null => "NULL".to_owned(),
            Self::Int64(i) => i.to_string(),
            Self::Double(d) => {
                if d.is_nan() {
                    "NULL".to_owned()
                } else {
                    to_str(*d)
                }
            }
            Self::Text(t) => t.clone(),
        }
    }

    /// The cell as an SQL literal (source `oswValueToSQLLiteral`): text is
    /// quoted with doubled single quotes, everything else is [`to_text`](Self::to_text).
    pub fn to_sql_literal(&self) -> String {
        match self {
            Self::Text(t) => format!("'{}'", t.replace('\'', "''")),
            other => other.to_text(),
        }
    }
}

impl From<i64> for OswValue {
    fn from(value: i64) -> Self {
        Self::Int64(value)
    }
}

impl From<u64> for OswValue {
    /// The bits as a signed integer, as the source's `UInt64` constructor.
    fn from(value: u64) -> Self {
        Self::Int64(value as i64)
    }
}

impl From<i32> for OswValue {
    fn from(value: i32) -> Self {
        Self::Int64(i64::from(value))
    }
}

impl From<f64> for OswValue {
    /// A real cell; NaN becomes NULL.
    fn from(value: f64) -> Self {
        Self::real(value)
    }
}

#[cfg(feature = "sqlite")]
impl rusqlite::ToSql for OswValue {
    fn to_sql(&self) -> rusqlite::Result<rusqlite::types::ToSqlOutput<'_>> {
        use rusqlite::types::{ToSqlOutput, Value, ValueRef};
        Ok(match self {
            Self::Null => ToSqlOutput::Owned(Value::Null),
            Self::Double(d) if d.is_nan() => ToSqlOutput::Owned(Value::Null),
            Self::Int64(i) => ToSqlOutput::Owned(Value::Integer(*i)),
            Self::Double(d) => ToSqlOutput::Owned(Value::Real(*d)),
            Self::Text(t) => ToSqlOutput::Borrowed(ValueRef::Text(t.as_bytes())),
        })
    }
}

/// `strtod` over a text that must be consumed completely (surrounding white
/// space already removed): decimal, `inf`/`infinity`, `nan` and `nan(chars)`
/// (any case, optional sign) and C99 hexadecimal floating point.
fn parse_c_double(text: &str) -> Option<f64> {
    if text.is_empty() {
        return None;
    }
    let lower = text.to_ascii_lowercase();
    let (negative, body) = match lower.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, lower.strip_prefix('+').unwrap_or(&lower)),
    };
    if let Some(rest) = body.strip_prefix("nan") {
        let ok = rest.is_empty()
            || (rest.starts_with('(')
                && rest.ends_with(')')
                && rest.len() >= 2
                && rest[1..rest.len() - 1]
                    .chars()
                    .all(|c| c.is_ascii_alphanumeric() || c == '_'));
        return ok.then_some(f64::NAN);
    }
    if let Some(hex) = body.strip_prefix("0x") {
        let value = parse_hex_float(hex)?;
        return Some(if negative { -value } else { value });
    }
    // Rust accepts the decimal and inf/infinity grammar of strtod
    if body.starts_with(|c: char| c.is_ascii_digit() || c == '.' || c == 'i') {
        return text.parse::<f64>().ok();
    }
    None
}

/// The value of a hexadecimal floating-point literal after its `0x`:
/// `h[.h]` with at least one hex digit and an optional `p[+-]d` binary
/// exponent.
fn parse_hex_float(text: &str) -> Option<f64> {
    let (mantissa, exponent) = match text.split_once('p') {
        Some((m, e)) => (m, Some(e)),
        None => (text, None),
    };
    let (integer, fraction) = match mantissa.split_once('.') {
        Some((i, f)) => (i, f),
        None => (mantissa, ""),
    };
    if integer.is_empty() && fraction.is_empty() {
        return None;
    }
    if !integer
        .chars()
        .chain(fraction.chars())
        .all(|c| c.is_ascii_hexdigit())
    {
        return None;
    }
    let mut value = 0.0_f64;
    for c in integer.chars() {
        value = value * 16.0 + f64::from(c.to_digit(16)?);
    }
    let mut scale = 1.0_f64 / 16.0;
    for c in fraction.chars() {
        value += f64::from(c.to_digit(16)?) * scale;
        scale /= 16.0;
    }
    let exponent: i32 = match exponent {
        None => 0,
        Some(e) => {
            let digits = e.strip_prefix(['+', '-']).unwrap_or(e);
            if digits.is_empty() || !digits.chars().all(|c| c.is_ascii_digit()) {
                return None;
            }
            e.parse::<i32>().ok()?
        }
    };
    Some(value * 2.0_f64.powi(exponent))
}

/// Whether a lower-case numeric text has a nonzero digit in its significand.
fn has_nonzero_significand(lower: &str) -> bool {
    let body = lower.trim_start_matches(['+', '-']);
    let hex = body.starts_with("0x");
    let body = body.strip_prefix("0x").unwrap_or(body);
    let significand = if hex {
        body.split('p').next().unwrap_or(body)
    } else {
        body.split('e').next().unwrap_or(body)
    };
    significand
        .chars()
        .any(|c| c.is_ascii_hexdigit() && c != '0')
}

fn null_row<const N: usize>() -> [OswValue; N] {
    std::array::from_fn(|_| OswValue::Null)
}

/// The buffered rows of all OSW tables, ready for insertion (source
/// `OpenSwathOSWWriter::OSWData`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OswData {
    /// Rows of `FEATURE`.
    pub feature_rows: Vec<FeatureRow>,
    /// Rows of `FEATURE_MS1`.
    pub feature_ms1_rows: Vec<FeatureMS1Row>,
    /// Rows of `FEATURE_MS2`.
    pub feature_ms2_rows: Vec<FeatureMS2Row>,
    /// Rows of `FEATURE_PRECURSOR`.
    pub feature_precursor_rows: Vec<FeaturePrecursorRow>,
    /// Rows of `FEATURE_TRANSITION`.
    pub feature_transition_rows: Vec<FeatureTransitionRow>,
}

/// Extra bytes per cell the memory estimate adds for string storage (source
/// `string_headroom_per_value`).
const STRING_HEADROOM_PER_VALUE: u64 = 16;

#[allow(clippy::ptr_arg)] // the capacity, not the length, is the point
fn estimate_row_vector<const N: usize>(rows: &Vec<[OswValue; N]>) -> u64 {
    let row_bytes =
        std::mem::size_of::<[OswValue; N]>() as u64 + N as u64 * STRING_HEADROOM_PER_VALUE;
    rows.capacity() as u64 * row_bytes
}

impl OswData {
    /// Reserve space for about `row_count` rows per table (source
    /// `reserve(Size)`).
    pub fn reserve(&mut self, row_count: usize) {
        self.reserve_rows(row_count, row_count);
    }

    /// Reserve space for `feature_row_count` rows of each feature table and
    /// `transition_row_count` transition rows (source `reserve(Size, Size)`).
    pub fn reserve_rows(&mut self, feature_row_count: usize, transition_row_count: usize) {
        self.feature_rows.reserve(feature_row_count);
        self.feature_ms1_rows.reserve(feature_row_count);
        self.feature_ms2_rows.reserve(feature_row_count);
        self.feature_precursor_rows.reserve(feature_row_count);
        self.feature_transition_rows.reserve(transition_row_count);
    }

    /// Append the rows of `rhs`, moving them (source `append`).
    pub fn append(&mut self, mut rhs: OswData) {
        self.feature_rows.append(&mut rhs.feature_rows);
        self.feature_ms1_rows.append(&mut rhs.feature_ms1_rows);
        self.feature_ms2_rows.append(&mut rhs.feature_ms2_rows);
        self.feature_precursor_rows
            .append(&mut rhs.feature_precursor_rows);
        self.feature_transition_rows
            .append(&mut rhs.feature_transition_rows);
    }

    /// Whether no rows are buffered.
    pub fn is_empty(&self) -> bool {
        self.feature_rows.is_empty()
            && self.feature_ms1_rows.is_empty()
            && self.feature_ms2_rows.is_empty()
            && self.feature_precursor_rows.is_empty()
            && self.feature_transition_rows.is_empty()
    }

    /// The number of buffered rows over all tables (source `rowCount`).
    pub fn row_count(&self) -> usize {
        self.feature_rows.len()
            + self.feature_ms1_rows.len()
            + self.feature_ms2_rows.len()
            + self.feature_precursor_rows.len()
            + self.feature_transition_rows.len()
    }

    /// An estimate of the bytes the buffer retains (source
    /// `estimateMemoryUsage`): the capacity of each row vector times the size
    /// of a row plus 16 bytes of string headroom per cell. The row size is
    /// that of this port's cells (32 bytes), not the source's (40).
    pub fn estimate_memory_usage(&self) -> u64 {
        estimate_row_vector(&self.feature_rows)
            + estimate_row_vector(&self.feature_ms1_rows)
            + estimate_row_vector(&self.feature_ms2_rows)
            + estimate_row_vector(&self.feature_precursor_rows)
            + estimate_row_vector(&self.feature_transition_rows)
    }
}

/// The cell of a metavalue (source `oswValue(const DataValue&)`): an empty or
/// missing value is NULL, an integer an integer, a real a real (NaN as NULL),
/// anything else the classified text of its `toString`.
fn osw_value_of(value: Option<&MetaValue>) -> OswValue {
    let Some(value) = value else {
        return OswValue::Null;
    };
    match value.data() {
        MetaValueData::Empty => OswValue::Null,
        MetaValueData::Integer(i) => OswValue::Int64(*i),
        MetaValueData::Float(f) => OswValue::real(*f),
        _ => OswValue::parse(&data_value_to_string(Some(value))),
    }
}

fn meta<'a>(feature: &'a Feature, name: &str) -> Option<&'a MetaValue> {
    feature.metadata.get(name)
}

fn osw_meta(feature: &Feature, name: &str) -> OswValue {
    osw_value_of(meta(feature, name))
}

/// The values of one score metavalue: nothing, one value, or a list (source
/// private `ScoreValueList`).
enum ScoreValues {
    Empty,
    Strings(Vec<String>),
    Ints(Vec<i64>),
    Doubles(Vec<f64>),
    Single(OswValue),
}

impl ScoreValues {
    fn new(feature: &Feature, name: &str) -> Self {
        let Some(value) = meta(feature, name) else {
            return Self::Empty;
        };
        match value.data() {
            MetaValueData::Empty => Self::Empty,
            MetaValueData::StringList(v) => Self::Strings(v.clone()),
            MetaValueData::IntegerList(v) => Self::Ints(v.clone()),
            MetaValueData::FloatList(v) => Self::Doubles(v.clone()),
            _ => Self::Single(osw_value_of(Some(value))),
        }
    }

    fn len(&self) -> usize {
        match self {
            Self::Empty => 0,
            Self::Strings(v) => v.len(),
            Self::Ints(v) => v.len(),
            Self::Doubles(v) => v.len(),
            Self::Single(_) => 1,
        }
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn at(&self, index: usize) -> OswValue {
        if index >= self.len() {
            return OswValue::Null;
        }
        match self {
            Self::Empty => OswValue::Null,
            Self::Strings(v) => OswValue::parse(&v[index]),
            Self::Ints(v) => OswValue::Int64(v[index]),
            Self::Doubles(v) => OswValue::real(v[index]),
            Self::Single(value) => value.clone(),
        }
    }
}

/// `(double)DataValue` or `(int)DataValue` of a metavalue: an integer or a real
/// reads as a number; anything else, and a missing value, is an error (the
/// source throws `ConversionError` for an empty value).
fn meta_number(feature: &Feature, name: &str) -> Result<f64> {
    match meta(feature, name).map(MetaValue::data) {
        Some(MetaValueData::Integer(i)) => Ok(*i as f64),
        Some(MetaValueData::Float(f)) => Ok(*f),
        _ => Err(Error::InvalidValue(format!(
            "metavalue '{name}' of a feature is missing or not a number"
        ))),
    }
}

fn level_is(feature: &Feature, level: &str) -> bool {
    matches!(
        meta(feature, "FeatureLevel").map(MetaValue::data),
        Some(MetaValueData::String(s)) if s == level
    )
}

/// The writer of the OSW output (source `OpenSwathOSWWriter`).
///
/// A writer built with an empty file name is inactive ([`is_active`](Self::is_active));
/// the workflow then writes no OSW output.
#[derive(Debug)]
pub struct OpenSwathOSWWriter {
    output_filename: String,
    run_id: u64,
    do_write: bool,
    enable_uis_scoring: bool,
    // persistent connection, reopened on first write; the lock replaces the
    // source's connector mutex
    #[cfg(feature = "sqlite")]
    conn: Mutex<Option<SqliteConnector>>,
}

impl OpenSwathOSWWriter {
    /// A writer for `output_filename`; `uis_scores` selects the
    /// identification (UIS) transition rows. Source constructor; nothing is
    /// opened or created yet.
    pub fn new(output_filename: impl Into<String>, uis_scores: bool) -> Self {
        let output_filename = output_filename.into();
        Self {
            do_write: !output_filename.is_empty(),
            output_filename,
            run_id: 0,
            enable_uis_scoring: uis_scores,
            #[cfg(feature = "sqlite")]
            conn: Mutex::new(None),
        }
    }

    /// Whether the writer has an output file (source `isActive`).
    pub fn is_active(&self) -> bool {
        self.do_write
    }

    /// The output file name.
    pub fn output_filename(&self) -> &str {
        &self.output_filename
    }

    /// The run id written into `FEATURE.RUN_ID` (source member `run_id_`).
    pub fn run_id(&self) -> u64 {
        self.run_id
    }

    /// Set the run id used when [`prepare_rows`](Self::prepare_rows) generates
    /// `FEATURE` rows, with the sign bit cleared (source `setRunId`).
    pub fn set_run_id(&mut self, run_id: u64) {
        self.run_id = clear_sign_bit(run_id) as u64;
    }

    /// The text of one score of `feature` for SQLite insertion (source
    /// `getScore`): `NULL` for an undefined score, an integer in decimal, a real
    /// through `StringUtils::toStr(double)`, text verbatim. A list-valued score
    /// is its classified `toString`; use
    /// [`get_separate_score`](Self::get_separate_score) for the elements.
    pub fn get_score(&self, feature: &Feature, score_name: &str) -> String {
        osw_meta(feature, score_name).to_text()
    }

    /// The texts of the elements of a list-valued score of `feature` (source
    /// `getSeparateScore`); one text for a scalar, none for an undefined score.
    pub fn get_separate_score(&self, feature: &Feature, score_name: &str) -> Vec<String> {
        let values = ScoreValues::new(feature, score_name);
        (0..values.len()).map(|i| values.at(i).to_text()).collect()
    }

    /// Prepare the rows of all features of `output` for insertion (source
    /// `prepareRows`). `pep` and `transition` are unused by the source.
    ///
    /// # Errors
    ///
    /// As [`prepare_rows_into`](Self::prepare_rows_into).
    pub fn prepare_rows(
        &self,
        pep: &LightCompound,
        transition: Option<&LightTransition>,
        output: &FeatureMap,
        id: &str,
    ) -> Result<OswData> {
        let mut rows = OswData::default();
        rows.reserve(output.features.len());
        self.prepare_rows_into(&mut rows, pep, transition, output, id)?;
        Ok(rows)
    }

    /// Append the rows of all features of `output` to `rows` (source
    /// `prepareRowsInto`); `id` is the transition group identifier written to
    /// `FEATURE.PRECURSOR_ID` (classified like any text cell, so a numeric
    /// group id becomes an integer).
    ///
    /// Per feature: one `FEATURE` and one `FEATURE_MS2` row; a `FEATURE_MS1`
    /// row when the feature has the score `var_ms1_ppm_diff`; per subordinate
    /// with `FeatureLevel` `MS2` one `FEATURE_TRANSITION` row (its
    /// `MASSERROR_PPM` is the element of the feature's `masserror_ppm` list at
    /// the position of the subordinate among the MS2 subordinates), per
    /// subordinate with `FeatureLevel` `MS1` and a positive intensity one
    /// `FEATURE_PRECURSOR` row whose `ISOTOPE` is the text after
    /// `Precursor_i` in its `native_id`. With UIS scoring the transition rows
    /// are those of the `id_target_*` and `id_decoy_*` list scores (as many
    /// as `id_target_num_transitions` or `id_decoy_num_transitions` say, but
    /// no more than there are transition names), and the MS2 subordinate rows
    /// are used only when no UIS row exists at all.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `norm_RT`, `delta_rt`,
    /// `id_target_num_transitions` or `id_decoy_num_transitions` is present
    /// but not a number (the source throws a `ConversionError`).
    pub fn prepare_rows_into(
        &self,
        rows: &mut OswData,
        _pep: &LightCompound,
        _transition: Option<&LightTransition>,
        output: &FeatureMap,
        id: &str,
    ) -> Result<()> {
        let mut uis_transition_rows: Vec<FeatureTransitionRow> = Vec::new();
        let mut ms2_transition_rows: Vec<FeatureTransitionRow> = Vec::new();
        let m = osw_meta;

        for feature in &output.features {
            let feature_id = OswValue::Int64(clear_sign_bit(feature.unique_id));
            let masserror_ppm = ScoreValues::new(feature, "masserror_ppm");

            let mut ms2_subordinate_index = 0;
            for sub in &feature.subordinates {
                if level_is(sub, "MS2") {
                    let mut row: FeatureTransitionRow = null_row();
                    row[0] = feature_id.clone();
                    row[1] = m(sub, "native_id");
                    row[2] = OswValue::real(f64::from(sub.intensity));
                    row[3] = m(sub, "total_xic");
                    row[4] = m(sub, "peak_apex_position");
                    row[5] = m(sub, "peak_apex_int");
                    row[6] = m(sub, "width_at_50");
                    row[7] = masserror_ppm.at(ms2_subordinate_index);
                    row[8] = m(sub, "total_mi");
                    if meta(sub, "start_position_at_5").is_some() {
                        row[30] = m(sub, "start_position_at_5");
                        row[31] = m(sub, "end_position_at_5");
                        row[32] = m(sub, "start_position_at_10");
                        row[33] = m(sub, "end_position_at_10");
                        row[34] = m(sub, "start_position_at_50");
                        row[35] = m(sub, "end_position_at_50");
                        row[36] = m(sub, "total_width");
                        row[37] = m(sub, "tailing_factor");
                        row[38] = m(sub, "asymmetry_factor");
                        row[39] = m(sub, "slope_of_baseline");
                        row[40] = m(sub, "baseline_delta_2_height");
                        row[41] = m(sub, "points_across_baseline");
                        row[42] = m(sub, "points_across_half_height");
                    }
                    if self.enable_uis_scoring {
                        ms2_transition_rows.push(row);
                    } else {
                        rows.feature_transition_rows.push(row);
                    }
                    ms2_subordinate_index += 1;
                } else if level_is(sub, "MS1") && sub.intensity > 0.0 {
                    let native_id = m(sub, "native_id").to_text();
                    let parts: Vec<&str> = if native_id.is_empty() {
                        Vec::new()
                    } else {
                        native_id.split("Precursor_i").collect()
                    };
                    rows.feature_precursor_rows.push([
                        feature_id.clone(),
                        if parts.len() > 1 {
                            OswValue::parse(parts[1])
                        } else {
                            OswValue::Null
                        },
                        OswValue::real(f64::from(sub.intensity)),
                        m(sub, "peak_apex_int"),
                    ]);
                }
            }

            // these are missing if RT scoring is disabled
            let mut norm_rt = -1.0;
            let mut delta_rt = -1.0;
            if meta(feature, "norm_RT").is_some() {
                norm_rt = meta_number(feature, "norm_RT")?;
            }
            if meta(feature, "delta_rt").is_some() {
                delta_rt = meta_number(feature, "delta_rt")?;
            }

            rows.feature_rows.push([
                feature_id.clone(),
                OswValue::from(self.run_id),
                OswValue::parse(id),
                OswValue::real(feature.rt),
                m(feature, "im_drift"),
                OswValue::real(norm_rt),
                OswValue::real(delta_rt),
                m(feature, "leftWidth"),
                m(feature, "rightWidth"),
                m(feature, "im_drift_left"),
                m(feature, "im_drift_right"),
            ]);

            rows.feature_ms2_rows.push([
                feature_id.clone(),
                OswValue::real(f64::from(feature.intensity)),
                m(feature, "total_xic"),
                m(feature, "peak_apices_sum"),
                m(feature, "im_drift"),
                m(feature, "im_drift_left"),
                m(feature, "im_drift_right"),
                m(feature, "im_delta"),
                m(feature, "total_mi"),
                m(feature, "var_bseries_score"),
                m(feature, "var_dotprod_score"),
                m(feature, "var_intensity_score"),
                m(feature, "var_isotope_correlation_score"),
                m(feature, "var_isotope_overlap_score"),
                m(feature, "var_library_corr"),
                m(feature, "var_library_dotprod"),
                m(feature, "var_library_manhattan"),
                m(feature, "var_library_rmsd"),
                m(feature, "var_library_rootmeansquare"),
                m(feature, "var_library_sangle"),
                m(feature, "var_log_sn_score"),
                m(feature, "var_manhatt_score"),
                m(feature, "var_massdev_score"),
                m(feature, "var_massdev_score_weighted"),
                m(feature, "var_mi_score"),
                m(feature, "var_mi_weighted_score"),
                m(feature, "var_mi_ratio_score"),
                m(feature, "var_norm_rt_score"),
                m(feature, "var_xcorr_coelution"),
                m(feature, "var_xcorr_coelution_weighted"),
                m(feature, "var_xcorr_shape"),
                m(feature, "var_xcorr_shape_weighted"),
                m(feature, "var_yseries_score"),
                m(feature, "var_elution_model_fit_score"),
                m(feature, "var_im_xcorr_shape"),
                m(feature, "var_im_xcorr_coelution"),
                m(feature, "var_im_delta_score"),
                m(feature, "im_log_intensity"),
            ]);

            // only write MS1 scores if they are present
            if meta(feature, "var_ms1_ppm_diff").is_some() {
                rows.feature_ms1_rows.push([
                    feature_id.clone(),
                    m(feature, "ms1_area_intensity"),
                    m(feature, "ms1_apex_intensity"),
                    m(feature, "im_ms1_drift"),
                    m(feature, "im_ms1_delta"),
                    m(feature, "var_ms1_ppm_diff"),
                    m(feature, "var_im_ms1_delta_score"),
                    m(feature, "var_ms1_mi_score"),
                    m(feature, "var_ms1_mi_contrast_score"),
                    m(feature, "var_ms1_mi_combined_score"),
                    m(feature, "var_ms1_isotope_correlation"),
                    m(feature, "var_ms1_isotope_overlap"),
                    m(feature, "var_ms1_xcorr_coelution"),
                    m(feature, "var_ms1_xcorr_coelution_contrast"),
                    m(feature, "var_ms1_xcorr_coelution_combined"),
                    m(feature, "var_ms1_xcorr_shape"),
                    m(feature, "var_ms1_xcorr_shape_contrast"),
                    m(feature, "var_ms1_xcorr_shape_combined"),
                ]);
            }

            if self.enable_uis_scoring {
                self.prepare_uis_rows(feature, &feature_id, &mut uis_transition_rows)?;
            }
        }

        if self.enable_uis_scoring {
            let selected = if uis_transition_rows.is_empty() {
                ms2_transition_rows
            } else {
                uis_transition_rows
            };
            rows.feature_transition_rows.extend(selected);
        }
        Ok(())
    }

    /// The identification transition rows of one feature: first the target
    /// transitions, then the decoy transitions, from the list scores
    /// `id_target_*` and `id_decoy_*`.
    fn prepare_uis_rows(
        &self,
        feature: &Feature,
        feature_id: &OswValue,
        out: &mut Vec<FeatureTransitionRow>,
    ) -> Result<()> {
        // the peak shape metrics are enabled by the *target* start positions
        let start_position_at_5 = ScoreValues::new(feature, "id_target_ind_start_position_at_5");
        let enable_peak_shape_metrics =
            !start_position_at_5.is_empty() && start_position_at_5.at(0).to_text() != "0";

        // (prefix, count key, transition-name key, shape-metric prefix)
        for (kind, count_key) in [
            ("target", "id_target_num_transitions"),
            ("decoy", "id_decoy_num_transitions"),
        ] {
            let list = |suffix: &str| ScoreValues::new(feature, &format!("id_{kind}_{suffix}"));
            let names = list("transition_names");
            let area = list("area_intensity");
            let total_area = list("total_area_intensity");
            let apex_intensity = list("apex_intensity");
            let apex_position = list("peak_apex_position");
            let fwhm = list("width_at_50");
            let total_mi = list("total_mi");
            let intensity_score = list("intensity_score");
            let intensity_ratio_score = list("intensity_ratio_score");
            let log_intensity = list("ind_log_intensity");
            let xcorr_coelution = list("ind_xcorr_coelution");
            let xcorr_shape = list("ind_xcorr_shape");
            let log_sn_score = list("ind_log_sn_score");
            let massdev_score = list("ind_massdev_score");
            let mi_score = list("ind_mi_score");
            let mi_ratio_score = list("ind_mi_ratio_score");
            let isotope_correlation = list("ind_isotope_correlation");
            let isotope_overlap = list("ind_isotope_overlap");
            let im_drift = list("ind_im_drift");
            let im_drift_left = list("ind_im_drift_left");
            let im_drift_right = list("ind_im_drift_right");
            let im_delta = list("ind_im_delta");
            let im_delta_score = list("ind_im_delta_score");
            let im_log_intensity = list("ind_im_log_intensity");
            let im_contrast_coelution = list("ind_im_contrast_coelution");
            let im_contrast_shape = list("ind_im_contrast_shape");
            let im_sum_contrast_coelution = list("ind_im_sum_contrast_coelution");
            let im_sum_contrast_shape = list("ind_im_sum_contrast_shape");
            let shape_metrics = if enable_peak_shape_metrics {
                Some([
                    list("ind_start_position_at_5"),
                    list("ind_end_position_at_5"),
                    list("ind_start_position_at_10"),
                    list("ind_end_position_at_10"),
                    list("ind_start_position_at_50"),
                    list("ind_end_position_at_50"),
                    list("ind_total_width"),
                    list("ind_tailing_factor"),
                    list("ind_asymmetry_factor"),
                    list("ind_slope_of_baseline"),
                    list("ind_baseline_delta_2_height"),
                    list("ind_points_across_baseline"),
                    list("ind_points_across_half_height"),
                ])
            } else {
                None
            };

            if meta(feature, count_key).is_none() {
                continue;
            }
            let num_transitions = meta_number(feature, count_key)? as i32;
            let available = names.len() as i64;
            let count = if available > 0 {
                i64::from(num_transitions).min(available)
            } else {
                0
            };
            for i in 0..count.max(0) as usize {
                let mut row: FeatureTransitionRow = null_row();
                row[0] = feature_id.clone();
                row[1] = names.at(i);
                row[2] = area.at(i);
                row[3] = total_area.at(i);
                row[4] = apex_position.at(i);
                row[5] = apex_intensity.at(i);
                row[6] = fwhm.at(i);
                row[7] = massdev_score.at(i);
                row[8] = total_mi.at(i);
                row[9] = intensity_score.at(i);
                row[10] = intensity_ratio_score.at(i);
                row[11] = log_intensity.at(i);
                row[12] = xcorr_coelution.at(i);
                row[13] = xcorr_shape.at(i);
                row[14] = log_sn_score.at(i);
                row[15] = massdev_score.at(i);
                row[16] = mi_score.at(i);
                row[17] = mi_ratio_score.at(i);
                row[18] = isotope_correlation.at(i);
                row[19] = isotope_overlap.at(i);
                row[20] = im_drift.at(i);
                row[21] = im_drift_left.at(i);
                row[22] = im_drift_right.at(i);
                row[23] = im_delta.at(i);
                row[24] = im_delta_score.at(i);
                row[25] = im_log_intensity.at(i);
                row[26] = im_contrast_coelution.at(i);
                row[27] = im_contrast_shape.at(i);
                row[28] = im_sum_contrast_coelution.at(i);
                row[29] = im_sum_contrast_shape.at(i);
                if let Some(shape) = &shape_metrics {
                    for (offset, values) in shape.iter().enumerate() {
                        row[30 + offset] = values.at(i);
                    }
                }
                out.push(row);
            }
        }
        Ok(())
    }

    /// The SQL text of the inserts of one prepared batch (source
    /// `prepareLine`): one `INSERT` per row, tables in the order `FEATURE`,
    /// `FEATURE_MS1`, `FEATURE_PRECURSOR`, `FEATURE_MS2`,
    /// `FEATURE_TRANSITION`, for [`write_lines`](Self::write_lines).
    ///
    /// # Errors
    ///
    /// As [`prepare_rows`](Self::prepare_rows).
    pub fn prepare_line(
        &self,
        pep: &LightCompound,
        transition: Option<&LightTransition>,
        output: &FeatureMap,
        id: &str,
    ) -> Result<String> {
        let rows = self.prepare_rows(pep, transition, output, id)?;
        let mut sql = String::new();
        append_insert_sql(&mut sql, "FEATURE", &FEATURE_COLUMNS, &rows.feature_rows);
        append_insert_sql(
            &mut sql,
            "FEATURE_MS1",
            &FEATURE_MS1_COLUMNS,
            &rows.feature_ms1_rows,
        );
        append_insert_sql(
            &mut sql,
            "FEATURE_PRECURSOR",
            &FEATURE_PRECURSOR_COLUMNS,
            &rows.feature_precursor_rows,
        );
        append_insert_sql(
            &mut sql,
            "FEATURE_MS2",
            &FEATURE_MS2_COLUMNS,
            &rows.feature_ms2_rows,
        );
        append_insert_sql(
            &mut sql,
            "FEATURE_TRANSITION",
            &FEATURE_TRANSITION_COLUMNS,
            &rows.feature_transition_rows,
        );
        Ok(sql)
    }
}

fn append_insert_sql<const N: usize>(
    sql: &mut String,
    table: &str,
    columns: &[&str; N],
    rows: &[[OswValue; N]],
) {
    for row in rows {
        sql.push_str("INSERT INTO ");
        sql.push_str(table);
        sql.push_str(" (");
        sql.push_str(&columns.join(", "));
        sql.push_str(") VALUES (");
        let values: Vec<String> = row.iter().map(OswValue::to_sql_literal).collect();
        sql.push_str(&values.join(", "));
        sql.push_str("); ");
    }
}

/// The `INSERT` statement of a table with `?` placeholders (source
/// `makeInsertStatement`).
pub fn insert_statement(table: &str, columns: &[&str]) -> String {
    let placeholders = vec!["?"; columns.len()].join(", ");
    format!(
        "INSERT INTO {table} ({}) VALUES ({placeholders});",
        columns.join(", ")
    )
}

#[cfg(feature = "sqlite")]
fn sql_error(error: rusqlite::Error) -> Error {
    Error::Io(std::io::Error::other(error))
}

#[cfg(feature = "sqlite")]
fn write_table_rows<const N: usize>(
    connection: &rusqlite::Connection,
    table: &str,
    columns: &[&str; N],
    rows: &[[OswValue; N]],
) -> Result<()> {
    if rows.is_empty() {
        return Ok(());
    }
    let mut statement = connection
        .prepare(&insert_statement(table, columns))
        .map_err(sql_error)?;
    for row in rows {
        statement
            .execute(rusqlite::params_from_iter(row.iter()))
            .map_err(sql_error)?;
    }
    Ok(())
}

#[cfg(feature = "sqlite")]
impl OpenSwathOSWWriter {
    /// Create the OSW tables in the output file (source `writeHeader`): `RUN`,
    /// `FEATURE`, `FEATURE_MS1`, `FEATURE_MS2`, `FEATURE_PRECURSOR` and
    /// `FEATURE_TRANSITION`, with the statements of the source verbatim
    /// ([`CREATE_SQL`]). Unlike [`add_run`](Self::add_run) this does not check
    /// [`is_active`](Self::is_active), as in the source.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the database cannot be opened or a table exists
    /// already.
    pub fn write_header(&self) -> Result<()> {
        let conn = SqliteConnector::new(&self.output_filename)?;
        conn.execute_statement(CREATE_SQL)
    }

    /// Register a run in the `RUN` table and make it the current run (source
    /// `addRun`). Does nothing for an inactive writer.
    ///
    /// The run id is bound as an SQLite INTEGER (sign bit cleared), so
    /// `FEATURE.RUN_ID = RUN.ID` joins work; the source stored it as a BLOB
    /// before that was fixed.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] on a SQLite error, such as a duplicate run id.
    pub fn add_run(&mut self, run_id: u64, input_filename: &str) -> Result<()> {
        if !self.do_write {
            return Ok(());
        }
        let conn = SqliteConnector::new(&self.output_filename)?;
        let rid = clear_sign_bit(run_id);
        conn.connection()
            .execute(
                "INSERT INTO RUN (ID, FILENAME) VALUES (?, ?);",
                rusqlite::params![rid, input_filename],
            )
            .map_err(sql_error)?;
        self.run_id = rid as u64;
        Ok(())
    }

    fn with_connection<R>(&self, f: impl FnOnce(&SqliteConnector) -> Result<R>) -> Result<R> {
        let mut guard = self
            .conn
            .lock()
            .map_err(|_| Error::InvalidValue("the OSW connection lock is poisoned".into()))?;
        if guard.is_none() {
            let conn = SqliteConnector::new(&self.output_filename)?;
            // improve write throughput: WAL and fewer synchronous flushes; a
            // failure only costs speed, as in the source
            let _ = conn.execute_statement("PRAGMA journal_mode=WAL");
            let _ = conn.execute_statement("PRAGMA synchronous=OFF");
            *guard = Some(conn);
        }
        match guard.as_ref() {
            Some(conn) => f(conn),
            None => Err(Error::InvalidValue("no OSW connection".into())),
        }
    }

    /// Insert buffered rows in one transaction (source `writeRows`), through
    /// one prepared statement per table. Empty data is a no-op.
    ///
    /// The connection is opened at the first call and kept (WAL journal, no
    /// synchronous flush). Concurrent callers are serialized by an internal
    /// lock; the source asks for an OpenMP critical section instead.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] on a SQLite error; the transaction is rolled back, so
    /// the file holds none of the batch.
    pub fn write_rows(&self, osw_output: &OswData) -> Result<()> {
        if osw_output.is_empty() {
            return Ok(());
        }
        self.with_connection(|conn| {
            conn.execute_statement("BEGIN TRANSACTION")?;
            let inserted = (|| {
                let c = conn.connection();
                write_table_rows(c, "FEATURE", &FEATURE_COLUMNS, &osw_output.feature_rows)?;
                write_table_rows(
                    c,
                    "FEATURE_MS1",
                    &FEATURE_MS1_COLUMNS,
                    &osw_output.feature_ms1_rows,
                )?;
                write_table_rows(
                    c,
                    "FEATURE_PRECURSOR",
                    &FEATURE_PRECURSOR_COLUMNS,
                    &osw_output.feature_precursor_rows,
                )?;
                write_table_rows(
                    c,
                    "FEATURE_MS2",
                    &FEATURE_MS2_COLUMNS,
                    &osw_output.feature_ms2_rows,
                )?;
                write_table_rows(
                    c,
                    "FEATURE_TRANSITION",
                    &FEATURE_TRANSITION_COLUMNS,
                    &osw_output.feature_transition_rows,
                )
            })();
            if let Err(error) = inserted {
                let _ = conn.execute_statement("ROLLBACK TRANSACTION");
                return Err(error);
            }
            conn.execute_statement("END TRANSACTION")
        })
    }

    /// Execute prepared SQL statements in one transaction (source
    /// `writeLines`), for the output of [`prepare_line`](Self::prepare_line).
    ///
    /// # Errors
    ///
    /// [`Error::Io`] on a SQLite error, with the batch rolled back.
    pub fn write_lines(&self, to_osw_output: &[String]) -> Result<()> {
        self.with_connection(|conn| {
            conn.execute_statement("BEGIN TRANSACTION")?;
            let result = (|| {
                for line in to_osw_output {
                    conn.execute_statement(line)?;
                }
                conn.execute_statement("END TRANSACTION")
            })();
            if let Err(error) = result {
                let _ = conn.execute_statement("ROLLBACK TRANSACTION");
                return Err(error);
            }
            Ok(())
        })
    }
}

#[cfg(not(feature = "sqlite"))]
impl OpenSwathOSWWriter {
    fn needs_sqlite() -> Error {
        Error::Unsupported("the OSW writer needs the `sqlite` feature".into())
    }

    /// Needs the `sqlite` feature; always [`Error::Unsupported`] without it.
    ///
    /// # Errors
    ///
    /// [`Error::Unsupported`].
    pub fn write_header(&self) -> Result<()> {
        Err(Self::needs_sqlite())
    }

    /// Needs the `sqlite` feature; always [`Error::Unsupported`] without it,
    /// except for an inactive writer, for which it does nothing as in the
    /// source.
    ///
    /// # Errors
    ///
    /// [`Error::Unsupported`].
    pub fn add_run(&mut self, run_id: u64, _input_filename: &str) -> Result<()> {
        if !self.do_write {
            return Ok(());
        }
        let _ = run_id;
        Err(Self::needs_sqlite())
    }

    /// Needs the `sqlite` feature; empty data is a no-op, anything else
    /// [`Error::Unsupported`] without it.
    ///
    /// # Errors
    ///
    /// [`Error::Unsupported`].
    pub fn write_rows(&self, osw_output: &OswData) -> Result<()> {
        if osw_output.is_empty() {
            return Ok(());
        }
        Err(Self::needs_sqlite())
    }

    /// Needs the `sqlite` feature; always [`Error::Unsupported`] without it.
    ///
    /// # Errors
    ///
    /// [`Error::Unsupported`].
    pub fn write_lines(&self, _to_osw_output: &[String]) -> Result<()> {
        Err(Self::needs_sqlite())
    }
}
