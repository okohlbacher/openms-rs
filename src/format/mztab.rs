// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! MzTab data model: the shared cell vocabulary of `FORMAT/MzTabBase.h` and the
//! record structs plus document type of `FORMAT/MzTab.h`.
//!
//! MzTab is the PSI tab-separated reporting format. A file is a metadata
//! section of `MTD` key/value lines followed by optional protein (`PRT`),
//! peptide (`PEP`), peptide-spectrum-match (`PSM`) and small-molecule (`SML`)
//! sections — and, in the OpenMS extension, nucleic-acid (`NUC`),
//! oligonucleotide (`OLI`) and oligonucleotide-spectrum-match (`OSM`) sections.
//! Each section is a header row plus data rows with a fixed column order and
//! arbitrary trailing `opt_` columns.
//!
//! The cell vocabulary is where the format's rules live. MzTab distinguishes
//! three textual states — `null`, `NaN` and `Inf` — from an optional column
//! that is absent altogether, and every cell renders and parses itself. Those
//! states are
//! [`MzTabCellState`](crate::format::mztab::MzTabCellState) and are carried by
//! [`MzTabDouble`](crate::format::mztab::MzTabDouble) and
//! [`MzTabInteger`](crate::format::mztab::MzTabInteger); the remaining cell
//! types spell "null" as an empty payload. Every cell type implements
//! [`MzTabCell`](crate::format::mztab::MzTabCell), so a writer can render and a
//! reader can fill a column without knowing its concrete type.
//!
//! This module is the data model only. It does not read or write `.mzTab`
//! files, and it does not export a `FeatureMap`, a `ConsensusMap` or an
//! identification run to MzTab; `docs/MZTAB_SUPPORT.md` records exactly which
//! source members are and are not covered, every preserved source convention
//! and every deliberate difference.
//!
//! ```
//! use openms::format::mztab::{MzTab, MzTabCell, MzTabDouble, MzTabPSMSectionRow};
//!
//! // A default numeric cell is null, not zero.
//! let mut cell = MzTabDouble::default();
//! assert!(cell.is_null());
//! assert_eq!(cell.to_cell_string(), "null");
//!
//! // "NaN" and "Inf" are states of their own, distinct from a value.
//! cell.read_cell("NaN")?;
//! assert!(cell.is_nan() && cell.get().is_err());
//!
//! // Values render with the source's fifteen-fractional-digit convention,
//! // trailing zeros trimmed to at least one digit.
//! cell.set(34.5);
//! assert_eq!(cell.to_cell_string(), "34.5");
//! cell.set(-2.0);
//! assert_eq!(cell.to_cell_string(), "-2.0");
//!
//! let mut document = MzTab::default();
//! assert_eq!(document.meta_data.mz_tab_version.get(), "1.0.0");
//! let mut row = MzTabPSMSectionRow::default();
//! row.sequence.set("NDYKAPPQPAPGK");
//! document.psm_data.push(row);
//! assert_eq!(document.psm_section_rows().len(), 1);
//! # Ok::<(), openms::Error>(())
//! ```

use crate::chemistry::{ModificationsDB, ResidueModification, TermSpecificity};
use crate::identification::{FlankingResidue, PeptideEvidence};
use crate::metadata::{MetaInfo, MetaValue, MetaValueData};
use crate::param::value::format_float;
use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

/// Largest cell text any `read_cell`/`from_cell_string` accepts, in bytes.
///
/// The source has no ceiling: every list cell parser splits the text and pushes
/// one owned entry per field. This port refuses oversized text in a preflight,
/// before anything is allocated or the receiver is touched.
pub const MAX_CELL_BYTES: usize = 4 * 1024 * 1024;

/// Largest number of separated entries a list cell may carry.
pub const MAX_CELL_ITEMS: usize = 100_000;

/// Largest number of modification names one metadata generator accepts.
pub const MAX_MODIFICATION_NAMES: usize = 100_000;

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

fn conversion(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

fn limit() -> Error {
    bad("MzTab cell exceeds its byte or entry limit")
}

/// `StringUtils::trim`: ASCII space, tab, newline and carriage return only.
///
/// Deliberately not [`str::trim`], which also removes Unicode whitespace such
/// as U+00A0. Trimming by characters cannot split a multi-byte character.
fn trim_source(text: &str) -> &str {
    text.trim_matches(|c| c == ' ' || c == '\t' || c == '\n' || c == '\r')
}

/// The source's `toLower(s) == "null"` test, without building a lowered copy.
///
/// The source lowercases the whole cell with `std::tolower` on each `unsigned
/// char`, which is locale-dependent for bytes above 0x7F and can corrupt UTF-8.
/// An ASCII-case-insensitive comparison agrees on every input for which the
/// test can succeed and never touches non-ASCII bytes.
fn is_token(text: &str, token: &str) -> bool {
    trim_source(text).eq_ignore_ascii_case(token)
}

/// `StringUtils::split`: an empty subject yields no fields at all, so an empty
/// cell adds no entries to a list rather than one empty entry.
fn split_cells(text: &str, separator: char) -> Vec<&str> {
    if text.is_empty() {
        Vec::new()
    } else {
        text.split(separator).collect()
    }
}

fn preflight_cell(text: &str, separator: char) -> Result<()> {
    if text.len() > MAX_CELL_BYTES {
        return Err(limit());
    }
    let items = text
        .matches(separator)
        .count()
        .checked_add(1)
        .ok_or_else(limit)?;
    if items > MAX_CELL_ITEMS {
        return Err(limit());
    }
    Ok(())
}

/// `StringUtils::toDouble`. Leading `+` and surrounding source whitespace are
/// accepted; a token that overflows the double range is an error, as the
/// source's `std::from_chars` reports `result_out_of_range`.
fn parse_source_double(text: &str) -> Result<f64> {
    let token = trim_source(text);
    let value = token
        .parse::<f64>()
        .map_err(|_| conversion(format!("could not convert {token:?} to a double value")))?;
    if value.is_infinite() {
        let spelling = token.trim_start_matches(['+', '-']);
        if !(spelling.eq_ignore_ascii_case("inf") || spelling.eq_ignore_ascii_case("infinity")) {
            return Err(conversion(format!("{token:?} is outside the double range")));
        }
    }
    Ok(value)
}

/// Header cell state for `Integer` and `Double` columns.
///
/// Source `MzTabCellStateType`. The source's trailing
/// `SIZE_OF_MZTAB_CELLTYPE` counter is not ported: it exists only to size C
/// arrays, and Rust enumerates variants without a sentinel.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum MzTabCellState {
    /// The cell carries a value.
    Default,
    /// The cell is textually `null`.
    #[default]
    Null,
    /// The cell is textually `NaN`.
    NaN,
    /// The cell is textually `Inf`.
    Inf,
}

impl MzTabCellState {
    /// Textual form written for a non-[`MzTabCellState::Default`] state.
    ///
    /// `Default` has no fixed spelling — the value decides — so it renders as
    /// the empty string here and callers must format the value themselves.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Default => "",
            Self::Null => "null",
            Self::NaN => "NaN",
            Self::Inf => "Inf",
        }
    }
}

/// A cell that renders itself to, and parses itself from, MzTab text.
///
/// The source declares `isNull`, `setNull`, `toCellString` and `fromCellString`
/// separately on each of its twelve cell classes with no common base. This
/// trait states the shared contract so a reader or writer can drive a column
/// generically; each type also keeps inherent methods with the source names.
pub trait MzTabCell: Sized {
    /// Whether the cell is in its "null" state.
    fn is_null(&self) -> bool;
    /// Enter (`true`) or leave (`false`) the "null" state.
    ///
    /// The source ignores `setNull(false)` for every type whose null state is
    /// an empty payload, because there is no value to restore. Only
    /// [`MzTabDouble`], [`MzTabInteger`] and [`MzTabBoolean`] act on `false`.
    fn set_null(&mut self, null: bool);
    /// Render the cell.
    ///
    /// # Errors
    ///
    /// Only [`MzTabModification`] and [`MzTabModificationList`] can fail: the
    /// source throws `Exception::ConversionError` when a modification carries
    /// position information but no identifier.
    fn write_cell(&self) -> Result<String>;
    /// Replace the cell's content from MzTab text.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] for text the source's `fromCellString` rejects with
    /// `Exception::ConversionError`, and [`Error::InvalidValue`] when the text
    /// exceeds [`MAX_CELL_BYTES`] or [`MAX_CELL_ITEMS`].
    fn read_cell(&mut self, text: &str) -> Result<()>;
}

/// Parse one cell into a fresh value.
///
/// Native convenience: the source only offers in-place `fromCellString`.
///
/// # Errors
///
/// As [`MzTabCell::read_cell`].
pub fn parse_cell<T: Default + MzTabCell>(text: &str) -> Result<T> {
    let mut cell = T::default();
    cell.read_cell(text)?;
    Ok(cell)
}

// ---------------------------------------------------------------------------
// MzTabBase.h — numeric cells
// ---------------------------------------------------------------------------

/// A `Double` cell: a finite value, or one of `null`, `NaN` and `Inf`.
///
/// The default cell is `null` with a stored value of `0.0`, as the source's
/// default constructor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MzTabDouble {
    value: f64,
    state: MzTabCellState,
}

impl Default for MzTabDouble {
    fn default() -> Self {
        Self {
            value: 0.0,
            state: MzTabCellState::Null,
        }
    }
}

impl MzTabDouble {
    /// A cell holding `value`, as the source's explicit `MzTabDouble(double)`.
    ///
    /// A non-finite `value` is accepted and enters
    /// [`MzTabCellState::Default`], exactly as the source: it then renders as
    /// `NaN`, `inf` or `-inf` and reparses into the corresponding *state*, so
    /// the cell text is stable but the stored value is not preserved across a
    /// round trip. Use [`MzTabDouble::nan`] or [`MzTabDouble::inf`] to name the
    /// state directly.
    pub fn new(value: f64) -> Self {
        Self {
            value,
            state: MzTabCellState::Default,
        }
    }
    /// The `null` cell, same as [`Default`](Default::default).
    pub fn null() -> Self {
        Self::default()
    }
    /// The `NaN` cell.
    pub fn nan() -> Self {
        Self {
            value: 0.0,
            state: MzTabCellState::NaN,
        }
    }
    /// The `Inf` cell.
    pub fn inf() -> Self {
        Self {
            value: 0.0,
            state: MzTabCellState::Inf,
        }
    }
    /// Store `value` and enter [`MzTabCellState::Default`].
    pub fn set(&mut self, value: f64) {
        self.value = value;
        self.state = MzTabCellState::Default;
    }
    /// The stored value.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when the cell is not in
    /// [`MzTabCellState::Default`], mapping the source's
    /// `Exception::ElementNotFound` — "Did you check the cell state before
    /// querying the value?". Check with [`MzTabDouble::is_null`],
    /// [`MzTabDouble::is_nan`] and [`MzTabDouble::is_inf`], or read
    /// [`MzTabDouble::state`].
    pub fn get(&self) -> Result<f64> {
        if self.state == MzTabCellState::Default {
            Ok(self.value)
        } else {
            Err(Error::MissingInformation(format!(
                "MzTab Double cell is {}, not a value",
                self.state.as_str()
            )))
        }
    }
    /// The cell state. Native accessor; the source exposes only the three
    /// boolean predicates.
    pub fn state(&self) -> MzTabCellState {
        self.state
    }
    /// The value stored behind the state, without the state check.
    ///
    /// Native accessor for the one place the source needs it: `operator<` and
    /// `operator==` compare this raw number and ignore the state entirely. See
    /// [`MzTabDouble::source_less`].
    pub fn raw_value(&self) -> f64 {
        self.value
    }
    /// Whether the cell is `null`.
    pub fn is_null(&self) -> bool {
        self.state == MzTabCellState::Null
    }
    /// Enter `null` (`true`) or [`MzTabCellState::Default`] (`false`).
    ///
    /// `set_null(false)` leaves the stored value untouched, as the source; on a
    /// freshly constructed cell it therefore publishes `0.0`.
    pub fn set_null(&mut self, null: bool) {
        self.state = if null {
            MzTabCellState::Null
        } else {
            MzTabCellState::Default
        };
    }
    /// Whether the cell is `NaN`.
    pub fn is_nan(&self) -> bool {
        self.state == MzTabCellState::NaN
    }
    /// Enter the `NaN` state. The stored value is left untouched.
    pub fn set_nan(&mut self) {
        self.state = MzTabCellState::NaN;
    }
    /// Whether the cell is `Inf`.
    pub fn is_inf(&self) -> bool {
        self.state == MzTabCellState::Inf
    }
    /// Enter the `Inf` state. The stored value is left untouched.
    pub fn set_inf(&mut self) {
        self.state = MzTabCellState::Inf;
    }
    /// Render the cell, using the source's `StringUtils::toStr(double)`
    /// convention: fifteen *fractional* digits trimmed to at least one for
    /// `|value|` in `[1e-2, 1e4)`, and shortest-round-trip scientific notation
    /// with a `+`-free two-or-more-digit exponent outside it.
    ///
    /// Fifteen fractional digits is not fifteen significant digits, so a value
    /// such as `51.9678841193106` renders as `51.967884119310597`, spelling out
    /// the binary value rather than its shortest form. That is the source
    /// convention and the reference files depend on it. A `NaN` or infinite
    /// value in [`MzTabCellState::Default`] renders `NaN`, `inf` or `-inf`
    /// rather than the `NaN`/`Inf` state keywords.
    pub fn to_cell_string(&self) -> String {
        match self.state {
            MzTabCellState::Default => format_float(self.value, true),
            other => other.as_str().to_owned(),
        }
    }
    /// Replace the cell from MzTab text.
    ///
    /// `null`, `nan` and `inf` are recognised case-insensitively after source
    /// trimming; anything else must parse as a double.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] when the text is neither a state keyword nor a double,
    /// mapping `Exception::ConversionError`.
    pub fn from_cell_string(&mut self, text: &str) -> Result<()> {
        if text.len() > MAX_CELL_BYTES {
            return Err(limit());
        }
        if is_token(text, "null") {
            *self = Self::null();
        } else if is_token(text, "nan") {
            *self = Self::nan();
        } else if is_token(text, "inf") {
            *self = Self::inf();
        } else {
            self.set(parse_source_double(text)?);
        }
        Ok(())
    }
    /// The source's `operator<`: compares the stored values and ignores both
    /// states, so a `null` cell sorts as `0.0`.
    pub fn source_less(&self, other: &Self) -> bool {
        self.value < other.value
    }
    /// The source's `operator==`: compares the stored values and ignores both
    /// states, so `MzTabDouble::null()` and `MzTabDouble::new(0.0)` are equal
    /// even though one renders `null` and the other `0.0`.
    ///
    /// Rust `==` on this type is stricter and compares the state as well.
    pub fn source_equal(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl From<f64> for MzTabDouble {
    fn from(value: f64) -> Self {
        Self::new(value)
    }
}

impl MzTabCell for MzTabDouble {
    fn is_null(&self) -> bool {
        Self::is_null(self)
    }
    fn set_null(&mut self, null: bool) {
        Self::set_null(self, null);
    }
    fn write_cell(&self) -> Result<String> {
        Ok(self.to_cell_string())
    }
    fn read_cell(&mut self, text: &str) -> Result<()> {
        self.from_cell_string(text)
    }
}

/// An `Integer` cell: a 32-bit value, or one of `null`, `NaN` and `Inf`.
///
/// The default cell is `null` with a stored value of `0`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct MzTabInteger {
    value: i32,
    state: MzTabCellState,
}

impl MzTabInteger {
    /// A cell holding `value`, as the source's explicit `MzTabInteger(int)`.
    pub fn new(value: i32) -> Self {
        Self {
            value,
            state: MzTabCellState::Default,
        }
    }
    /// The `null` cell, same as [`Default`](Default::default).
    pub fn null() -> Self {
        Self::default()
    }
    /// The `NaN` cell.
    pub fn nan() -> Self {
        Self {
            value: 0,
            state: MzTabCellState::NaN,
        }
    }
    /// The `Inf` cell.
    pub fn inf() -> Self {
        Self {
            value: 0,
            state: MzTabCellState::Inf,
        }
    }
    /// Store `value` and enter [`MzTabCellState::Default`].
    pub fn set(&mut self, value: i32) {
        self.value = value;
        self.state = MzTabCellState::Default;
    }
    /// The stored value.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when the cell is not in
    /// [`MzTabCellState::Default`], mapping the source's
    /// `Exception::ElementNotFound`.
    pub fn get(&self) -> Result<i32> {
        if self.state == MzTabCellState::Default {
            Ok(self.value)
        } else {
            Err(Error::MissingInformation(format!(
                "MzTab Integer cell is {}, not a value",
                self.state.as_str()
            )))
        }
    }
    /// The cell state. Native accessor.
    pub fn state(&self) -> MzTabCellState {
        self.state
    }
    /// The value stored behind the state, without the state check. Native.
    pub fn raw_value(&self) -> i32 {
        self.value
    }
    /// Whether the cell is `null`.
    pub fn is_null(&self) -> bool {
        self.state == MzTabCellState::Null
    }
    /// Enter `null` (`true`) or [`MzTabCellState::Default`] (`false`).
    pub fn set_null(&mut self, null: bool) {
        self.state = if null {
            MzTabCellState::Null
        } else {
            MzTabCellState::Default
        };
    }
    /// Whether the cell is `NaN`.
    pub fn is_nan(&self) -> bool {
        self.state == MzTabCellState::NaN
    }
    /// Enter the `NaN` state. The stored value is left untouched.
    pub fn set_nan(&mut self) {
        self.state = MzTabCellState::NaN;
    }
    /// Whether the cell is `Inf`.
    pub fn is_inf(&self) -> bool {
        self.state == MzTabCellState::Inf
    }
    /// Enter the `Inf` state. The stored value is left untouched.
    pub fn set_inf(&mut self) {
        self.state = MzTabCellState::Inf;
    }
    /// Render the cell as plain decimal digits, `null`, `NaN` or `Inf`.
    pub fn to_cell_string(&self) -> String {
        match self.state {
            MzTabCellState::Default => self.value.to_string(),
            other => other.as_str().to_owned(),
        }
    }
    /// Replace the cell from MzTab text.
    ///
    /// The source parses the text as a *double* first and then requires it to
    /// be integral, because mzTab files from external sources write `4.0` in
    /// integer columns; that leniency is preserved, so `4.0` reads as `4`.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] when the text is neither a state keyword nor an
    /// integral number, mapping `Exception::ConversionError`. A value outside
    /// the `i32` range is rejected here; the source casts the double to `int`,
    /// which is undefined behaviour and on the usual platforms wraps and then
    /// trips the same equality check.
    pub fn from_cell_string(&mut self, text: &str) -> Result<()> {
        if text.len() > MAX_CELL_BYTES {
            return Err(limit());
        }
        if is_token(text, "null") {
            *self = Self::null();
        } else if is_token(text, "nan") {
            *self = Self::nan();
        } else if is_token(text, "inf") {
            *self = Self::inf();
        } else {
            let value = parse_source_double(text)?;
            if value.fract() != 0.0 || !(-2_147_483_648.0..=2_147_483_647.0).contains(&value) {
                return Err(conversion(format!(
                    "could not convert {:?} to an MzTab Integer",
                    trim_source(text)
                )));
            }
            self.set(value as i32);
        }
        Ok(())
    }
}

impl From<i32> for MzTabInteger {
    fn from(value: i32) -> Self {
        Self::new(value)
    }
}

impl MzTabCell for MzTabInteger {
    fn is_null(&self) -> bool {
        Self::is_null(self)
    }
    fn set_null(&mut self, null: bool) {
        Self::set_null(self, null);
    }
    fn write_cell(&self) -> Result<String> {
        Ok(self.to_cell_string())
    }
    fn read_cell(&mut self, text: &str) -> Result<()> {
        self.from_cell_string(text)
    }
}

/// A `Boolean` cell: `0`, `1` or `null`.
///
/// The source stores an `int` that is negative when null, and its `get()`
/// returns that raw number rather than a boolean; [`MzTabBoolean::value`] keeps
/// that signature and [`MzTabBoolean::as_bool`] is the native reading.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct MzTabBoolean {
    value: i32,
}

impl Default for MzTabBoolean {
    fn default() -> Self {
        Self { value: -1 }
    }
}

impl MzTabBoolean {
    /// A cell holding `value`, as the source's explicit `MzTabBoolean(bool)`.
    pub fn new(value: bool) -> Self {
        Self {
            value: i32::from(value),
        }
    }
    /// The `null` cell, same as [`Default`](Default::default).
    pub fn null() -> Self {
        Self::default()
    }
    /// Store `value`.
    pub fn set(&mut self, value: bool) {
        self.value = i32::from(value);
    }
    /// The source's `get()`: the raw stored integer, `-1` when the cell is null.
    pub fn value(&self) -> i32 {
        self.value
    }
    /// The cell as a boolean, `None` when null. Native accessor.
    pub fn as_bool(&self) -> Option<bool> {
        if self.value < 0 {
            None
        } else {
            Some(self.value != 0)
        }
    }
    /// Whether the cell is `null`, i.e. the stored integer is negative.
    pub fn is_null(&self) -> bool {
        self.value < 0
    }
    /// Enter `null` (`true`) or become `0` (`false`).
    ///
    /// `set_null(false)` stores `0`, not the previous value: the source has
    /// nowhere to keep it. A regression test upstream pins this polarity,
    /// because the two branches were once inverted.
    pub fn set_null(&mut self, null: bool) {
        self.value = if null { -1 } else { 0 };
    }
    /// Render the cell as `null`, `1` or `0`.
    pub fn to_cell_string(&self) -> String {
        if self.is_null() {
            "null".to_owned()
        } else if self.value != 0 {
            "1".to_owned()
        } else {
            "0".to_owned()
        }
    }
    /// Replace the cell from MzTab text.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] for anything but `null` (case-insensitive, source
    /// trimming applied) and the exact texts `0` and `1`. The source tests the
    /// two digits against the *untrimmed* text even though it trims for the
    /// `null` test, so `" 1"` is a conversion error; that asymmetry is
    /// preserved because a reader that trims would accept cells the source
    /// rejects.
    pub fn from_cell_string(&mut self, text: &str) -> Result<()> {
        if text.len() > MAX_CELL_BYTES {
            return Err(limit());
        }
        if is_token(text, "null") {
            self.set_null(true);
        } else if text == "0" {
            self.set(false);
        } else if text == "1" {
            self.set(true);
        } else {
            return Err(conversion(format!(
                "could not convert {text:?} to an MzTab Boolean"
            )));
        }
        Ok(())
    }
}

impl From<bool> for MzTabBoolean {
    fn from(value: bool) -> Self {
        Self::new(value)
    }
}

impl MzTabCell for MzTabBoolean {
    fn is_null(&self) -> bool {
        Self::is_null(self)
    }
    fn set_null(&mut self, null: bool) {
        Self::set_null(self, null);
    }
    fn write_cell(&self) -> Result<String> {
        Ok(self.to_cell_string())
    }
    fn read_cell(&mut self, text: &str) -> Result<()> {
        self.from_cell_string(text)
    }
}

/// A `String` cell. An empty payload is the `null` state.
///
/// Storing text trims source whitespace, and storing the literal `null` in any
/// letter case stores nothing, so `MzTabString::from_text("null")` is null.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabString {
    value: String,
}

impl MzTabString {
    /// A cell holding `text`, as the source's explicit
    /// `MzTabString(const std::string&)`, which forwards to `set`.
    pub fn from_text(text: &str) -> Self {
        let mut cell = Self::default();
        cell.set(text);
        cell
    }
    /// The `null` cell, same as [`Default`](Default::default).
    pub fn null() -> Self {
        Self::default()
    }
    /// Store `text`, trimmed. The literal `null` in any case stores nothing.
    pub fn set(&mut self, text: &str) {
        if is_token(text, "null") {
            self.value.clear();
        } else {
            self.value.clear();
            self.value.push_str(trim_source(text));
        }
    }
    /// The stored text, empty when the cell is null.
    pub fn get(&self) -> &str {
        &self.value
    }
    /// Whether the cell is null, i.e. the stored text is empty.
    pub fn is_null(&self) -> bool {
        self.value.is_empty()
    }
    /// Clear the cell when `null` is `true`; `false` is ignored, as the source.
    pub fn set_null(&mut self, null: bool) {
        if null {
            self.value.clear();
        }
    }
    /// Render the cell as `null` or the stored text.
    pub fn to_cell_string(&self) -> String {
        if self.is_null() {
            "null".to_owned()
        } else {
            self.value.clone()
        }
    }
    /// Replace the cell from MzTab text. Cannot fail; the source's
    /// `fromCellString` is `set`.
    pub fn from_cell_string(&mut self, text: &str) {
        self.set(text);
    }
}

impl From<&str> for MzTabString {
    fn from(text: &str) -> Self {
        Self::from_text(text)
    }
}

impl MzTabCell for MzTabString {
    fn is_null(&self) -> bool {
        Self::is_null(self)
    }
    fn set_null(&mut self, null: bool) {
        Self::set_null(self, null);
    }
    fn write_cell(&self) -> Result<String> {
        Ok(self.to_cell_string())
    }
    fn read_cell(&mut self, text: &str) -> Result<()> {
        if text.len() > MAX_CELL_BYTES {
            return Err(limit());
        }
        self.set(text);
        Ok(())
    }
}

/// One trailing `opt_…` column of a section row: a name and a nullable value.
///
/// Source `MzTabOptionalColumnEntry`, a `std::pair<std::string, MzTabString>`
/// whose comment records that the name is not nullable and the value is. A
/// named struct replaces the pair so `first`/`second` do not travel with it.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabOptionalColumnEntry {
    /// Column name, which the format requires to start with `opt_`.
    pub name: String,
    /// Column value for this row; may be null.
    pub value: MzTabString,
}

impl MzTabOptionalColumnEntry {
    /// An entry with `name` and `value`.
    pub fn new(name: impl Into<String>, value: MzTabString) -> Self {
        Self {
            name: name.into(),
            value,
        }
    }
}

// ---------------------------------------------------------------------------
// MzTabBase.h — parameters, lists and spectra references
// ---------------------------------------------------------------------------

/// A `Param` cell, rendered `[CV label, accession, name, value]`.
///
/// Null when all four parts are empty. The source's getters carry a debug-only
/// `assert(!isNull())`; a release build returns the empty strings, which is
/// what these accessors do unconditionally.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabParameter {
    cv_label: String,
    accession: String,
    name: String,
    value: String,
}

impl MzTabParameter {
    /// The null parameter, same as [`Default`](Default::default).
    pub fn null() -> Self {
        Self::default()
    }
    /// A parameter from its four parts. Native constructor; the source only
    /// offers the four setters.
    pub fn from_parts(
        cv_label: impl Into<String>,
        accession: impl Into<String>,
        name: impl Into<String>,
        value: impl Into<String>,
    ) -> Self {
        Self {
            cv_label: cv_label.into(),
            accession: accession.into(),
            name: name.into(),
            value: value.into(),
        }
    }
    /// Parse a parameter cell. Native convenience for `fromCellString`.
    ///
    /// # Errors
    ///
    /// As [`MzTabParameter::from_cell_string`].
    pub fn parse(text: &str) -> Result<Self> {
        let mut parameter = Self::default();
        parameter.from_cell_string(text)?;
        Ok(parameter)
    }
    /// Set the controlled-vocabulary label, for example `MS` or `UNIMOD`.
    pub fn set_cv_label(&mut self, cv_label: impl Into<String>) {
        self.cv_label = cv_label.into();
    }
    /// Set the term accession, for example `MS:1002453`.
    pub fn set_accession(&mut self, accession: impl Into<String>) {
        self.accession = accession.into();
    }
    /// Set the term name.
    pub fn set_name(&mut self, name: impl Into<String>) {
        self.name = name.into();
    }
    /// Set the term value.
    pub fn set_value(&mut self, value: impl Into<String>) {
        self.value = value.into();
    }
    /// The controlled-vocabulary label.
    pub fn cv_label(&self) -> &str {
        &self.cv_label
    }
    /// The term accession.
    pub fn accession(&self) -> &str {
        &self.accession
    }
    /// The term name.
    pub fn name(&self) -> &str {
        &self.name
    }
    /// The term value.
    pub fn value(&self) -> &str {
        &self.value
    }
    /// Whether all four parts are empty.
    pub fn is_null(&self) -> bool {
        self.cv_label.is_empty()
            && self.accession.is_empty()
            && self.name.is_empty()
            && self.value.is_empty()
    }
    /// Clear all four parts when `null` is `true`; `false` is ignored.
    pub fn set_null(&mut self, null: bool) {
        if null {
            self.cv_label.clear();
            self.accession.clear();
            self.name.clear();
            self.value.clear();
        }
    }
    /// Render the parameter.
    ///
    /// Name and value are wrapped in double quotes when they contain a comma
    /// *followed by a space*, which is the separator the parser splits on. A
    /// name containing a bare `,` is not quoted and produces a cell that
    /// [`MzTabParameter::from_cell_string`] rejects; this port reproduces that
    /// exactly rather than quoting more eagerly, because the reference files
    /// depend on the byte layout.
    pub fn to_cell_string(&self) -> String {
        if self.is_null() {
            return "null".to_owned();
        }
        let quoted = |text: &str| {
            if text.contains(", ") {
                format!("\"{text}\"")
            } else {
                text.to_owned()
            }
        };
        format!(
            "[{}, {}, {}, {}]",
            self.cv_label,
            self.accession,
            quoted(&self.name),
            quoted(&self.value)
        )
    }
    /// Replace the parameter from MzTab text.
    ///
    /// The source's scanner splits on unquoted commas, drops every `[` and `]`
    /// wherever they occur — inside quotes included — skips leading spaces in
    /// each field and trims each field. All of that is preserved.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] when the text does not yield exactly four fields,
    /// mapping `Exception::ConversionError`, and [`Error::InvalidValue`] when
    /// it exceeds [`MAX_CELL_BYTES`] or [`MAX_CELL_ITEMS`].
    pub fn from_cell_string(&mut self, text: &str) -> Result<()> {
        if is_token(text, "null") {
            self.set_null(true);
            return Ok(());
        }
        preflight_cell(text, ',')?;
        let mut fields: Vec<String> = Vec::new();
        let mut field = String::new();
        let mut in_quotes = false;
        for c in text.chars() {
            if c == '"' {
                in_quotes = !in_quotes;
            } else if c == ',' {
                if in_quotes {
                    field.push(',');
                } else {
                    fields.push(trim_source(&field).to_owned());
                    field.clear();
                }
            } else if c != '[' && c != ']' {
                if c == ' ' && field.is_empty() {
                    continue;
                }
                field.push(c);
            }
        }
        fields.push(trim_source(&field).to_owned());
        if fields.len() != 4 {
            return Err(conversion(format!(
                "could not convert {text:?} to an MzTab Param: {} fields instead of 4",
                fields.len()
            )));
        }
        let value = fields.pop().unwrap_or_default();
        let name = fields.pop().unwrap_or_default();
        let accession = fields.pop().unwrap_or_default();
        let cv_label = fields.pop().unwrap_or_default();
        *self = Self {
            cv_label,
            accession,
            name,
            value,
        };
        Ok(())
    }
}

impl MzTabCell for MzTabParameter {
    fn is_null(&self) -> bool {
        Self::is_null(self)
    }
    fn set_null(&mut self, null: bool) {
        Self::set_null(self, null);
    }
    fn write_cell(&self) -> Result<String> {
        Ok(self.to_cell_string())
    }
    fn read_cell(&mut self, text: &str) -> Result<()> {
        self.from_cell_string(text)
    }
}

/// A `Param[]` cell: parameters separated by `|`. Null when empty.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabParameterList {
    parameters: Vec<MzTabParameter>,
}

impl MzTabParameterList {
    /// An empty, and therefore null, list.
    pub fn null() -> Self {
        Self::default()
    }
    /// The parameters. The source returns a copy of the vector.
    pub fn get(&self) -> &[MzTabParameter] {
        &self.parameters
    }
    /// Mutable access to the parameters. Native accessor.
    pub fn get_mut(&mut self) -> &mut Vec<MzTabParameter> {
        &mut self.parameters
    }
    /// Replace the parameters.
    pub fn set(&mut self, parameters: Vec<MzTabParameter>) {
        self.parameters = parameters;
    }
    /// Whether the list is empty.
    pub fn is_null(&self) -> bool {
        self.parameters.is_empty()
    }
    /// Clear the list when `null` is `true`; `false` is ignored.
    pub fn set_null(&mut self, null: bool) {
        if null {
            self.parameters.clear();
        }
    }
    /// Render the list as `null` or `|`-separated parameters.
    pub fn to_cell_string(&self) -> String {
        if self.is_null() {
            return "null".to_owned();
        }
        self.parameters
            .iter()
            .map(MzTabParameter::to_cell_string)
            .collect::<Vec<_>>()
            .join("|")
    }
    /// Replace the list from MzTab text.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] when any field is the literal `null` — the source
    /// throws "MzTabParameter in MzTabParameterList must not be null" — or when
    /// a field is not a parameter, and [`Error::InvalidValue`] for the
    /// [`MAX_CELL_BYTES`]/[`MAX_CELL_ITEMS`] ceilings.
    ///
    /// The parsed entries *replace* the previous ones. The source pushes onto
    /// the existing vector without clearing it, so parsing twice into the same
    /// object concatenates; that is a defect, not a convention.
    pub fn from_cell_string(&mut self, text: &str) -> Result<()> {
        if is_token(text, "null") {
            self.set_null(true);
            return Ok(());
        }
        preflight_cell(text, '|')?;
        let mut next = Vec::new();
        for field in split_cells(text, '|') {
            if is_token(field, "null") {
                return Err(conversion(format!(
                    "MzTab Param in a Param[] must not be null: {text:?}"
                )));
            }
            let mut parameter = MzTabParameter::default();
            parameter.from_cell_string(field)?;
            next.push(parameter);
        }
        self.parameters = next;
        Ok(())
    }
}

impl MzTabCell for MzTabParameterList {
    fn is_null(&self) -> bool {
        Self::is_null(self)
    }
    fn set_null(&mut self, null: bool) {
        Self::set_null(self, null);
    }
    fn write_cell(&self) -> Result<String> {
        Ok(self.to_cell_string())
    }
    fn read_cell(&mut self, text: &str) -> Result<()> {
        self.from_cell_string(text)
    }
}

/// A `String[]` cell with a selectable separator. Null when empty.
///
/// The separator defaults to `|`; `,` is needed for `ambiguity_members` and GO
/// accessions, whose own values may contain `|`.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabStringList {
    entries: Vec<MzTabString>,
    separator: char,
}

impl Default for MzTabStringList {
    fn default() -> Self {
        Self {
            entries: Vec::new(),
            separator: '|',
        }
    }
}

impl MzTabStringList {
    /// An empty, and therefore null, list with the default `|` separator.
    pub fn null() -> Self {
        Self::default()
    }
    /// Choose the separator used by both rendering and parsing.
    ///
    /// Needed for `ambiguity_members` and GO accessions, which use `,` while
    /// the other string lists use `|`. The source takes a `char`, so only
    /// single-byte separators have a source equivalent; a multi-byte `char`
    /// works here but has none.
    pub fn set_separator(&mut self, separator: char) {
        self.separator = separator;
    }
    /// The current separator. Native accessor.
    pub fn separator(&self) -> char {
        self.separator
    }
    /// The entries. The source returns a copy of the vector.
    pub fn get(&self) -> &[MzTabString] {
        &self.entries
    }
    /// Mutable access to the entries. Native accessor.
    pub fn get_mut(&mut self) -> &mut Vec<MzTabString> {
        &mut self.entries
    }
    /// Replace the entries. The separator is untouched.
    pub fn set(&mut self, entries: Vec<MzTabString>) {
        self.entries = entries;
    }
    /// Whether the list is empty.
    pub fn is_null(&self) -> bool {
        self.entries.is_empty()
    }
    /// Clear the list when `null` is `true`; `false` is ignored.
    pub fn set_null(&mut self, null: bool) {
        if null {
            self.entries.clear();
        }
    }
    /// Render the list as `null` or separator-joined entries.
    pub fn to_cell_string(&self) -> String {
        if self.is_null() {
            return "null".to_owned();
        }
        self.entries
            .iter()
            .map(MzTabString::to_cell_string)
            .collect::<Vec<_>>()
            .join(&self.separator.to_string())
    }
    /// Replace the list from MzTab text, splitting on the current separator.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for the [`MAX_CELL_BYTES`]/[`MAX_CELL_ITEMS`]
    /// ceilings. Parsing a string entry cannot otherwise fail.
    ///
    /// The parsed entries *replace* the previous ones; the source appends.
    pub fn from_cell_string(&mut self, text: &str) -> Result<()> {
        if is_token(text, "null") {
            self.set_null(true);
            return Ok(());
        }
        preflight_cell(text, self.separator)?;
        let next = split_cells(text, self.separator)
            .into_iter()
            .map(MzTabString::from_text)
            .collect();
        self.entries = next;
        Ok(())
    }
}

impl MzTabCell for MzTabStringList {
    fn is_null(&self) -> bool {
        Self::is_null(self)
    }
    fn set_null(&mut self, null: bool) {
        Self::set_null(self, null);
    }
    fn write_cell(&self) -> Result<String> {
        Ok(self.to_cell_string())
    }
    fn read_cell(&mut self, text: &str) -> Result<()> {
        self.from_cell_string(text)
    }
}

/// An `Integer[]` cell: integer cells separated by `,`. Null when empty.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct MzTabIntegerList {
    entries: Vec<MzTabInteger>,
}

impl MzTabIntegerList {
    /// An empty, and therefore null, list.
    pub fn null() -> Self {
        Self::default()
    }
    /// The entries. The source returns a copy of the vector.
    pub fn get(&self) -> &[MzTabInteger] {
        &self.entries
    }
    /// Mutable access to the entries. Native accessor.
    pub fn get_mut(&mut self) -> &mut Vec<MzTabInteger> {
        &mut self.entries
    }
    /// Replace the entries.
    pub fn set(&mut self, entries: Vec<MzTabInteger>) {
        self.entries = entries;
    }
    /// Whether the list is empty.
    pub fn is_null(&self) -> bool {
        self.entries.is_empty()
    }
    /// Clear the list when `null` is `true`; `false` is ignored.
    pub fn set_null(&mut self, null: bool) {
        if null {
            self.entries.clear();
        }
    }
    /// Render the list as `null` or `,`-separated integer cells.
    pub fn to_cell_string(&self) -> String {
        if self.is_null() {
            return "null".to_owned();
        }
        self.entries
            .iter()
            .map(MzTabInteger::to_cell_string)
            .collect::<Vec<_>>()
            .join(",")
    }
    /// Replace the list from MzTab text.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] when a field is not an integer cell, and
    /// [`Error::InvalidValue`] for the [`MAX_CELL_BYTES`]/[`MAX_CELL_ITEMS`]
    /// ceilings. The parsed entries *replace* the previous ones; the source
    /// appends.
    pub fn from_cell_string(&mut self, text: &str) -> Result<()> {
        if is_token(text, "null") {
            self.set_null(true);
            return Ok(());
        }
        preflight_cell(text, ',')?;
        let mut next = Vec::new();
        for field in split_cells(text, ',') {
            let mut cell = MzTabInteger::default();
            cell.from_cell_string(field)?;
            next.push(cell);
        }
        self.entries = next;
        Ok(())
    }
}

impl MzTabCell for MzTabIntegerList {
    fn is_null(&self) -> bool {
        Self::is_null(self)
    }
    fn set_null(&mut self, null: bool) {
        Self::set_null(self, null);
    }
    fn write_cell(&self) -> Result<String> {
        Ok(self.to_cell_string())
    }
    fn read_cell(&mut self, text: &str) -> Result<()> {
        self.from_cell_string(text)
    }
}

/// A `Double[]` cell: double cells separated by `|`. Null when empty.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MzTabDoubleList {
    entries: Vec<MzTabDouble>,
}

impl MzTabDoubleList {
    /// An empty, and therefore null, list.
    pub fn null() -> Self {
        Self::default()
    }
    /// The entries. The source returns a copy of the vector.
    pub fn get(&self) -> &[MzTabDouble] {
        &self.entries
    }
    /// Mutable access to the entries. Native accessor.
    pub fn get_mut(&mut self) -> &mut Vec<MzTabDouble> {
        &mut self.entries
    }
    /// Replace the entries.
    pub fn set(&mut self, entries: Vec<MzTabDouble>) {
        self.entries = entries;
    }
    /// Whether the list is empty.
    pub fn is_null(&self) -> bool {
        self.entries.is_empty()
    }
    /// Clear the list when `null` is `true`; `false` is ignored.
    pub fn set_null(&mut self, null: bool) {
        if null {
            self.entries.clear();
        }
    }
    /// Render the list as `null` or `|`-separated double cells.
    pub fn to_cell_string(&self) -> String {
        if self.is_null() {
            return "null".to_owned();
        }
        self.entries
            .iter()
            .map(MzTabDouble::to_cell_string)
            .collect::<Vec<_>>()
            .join("|")
    }
    /// Replace the list from MzTab text.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] when a field is not a double cell, and
    /// [`Error::InvalidValue`] for the [`MAX_CELL_BYTES`]/[`MAX_CELL_ITEMS`]
    /// ceilings. The parsed entries *replace* the previous ones; the source
    /// appends.
    pub fn from_cell_string(&mut self, text: &str) -> Result<()> {
        if is_token(text, "null") {
            self.set_null(true);
            return Ok(());
        }
        preflight_cell(text, '|')?;
        let mut next = Vec::new();
        for field in split_cells(text, '|') {
            let mut cell = MzTabDouble::default();
            cell.from_cell_string(field)?;
            next.push(cell);
        }
        self.entries = next;
        Ok(())
    }
}

impl MzTabCell for MzTabDoubleList {
    fn is_null(&self) -> bool {
        Self::is_null(self)
    }
    fn set_null(&mut self, null: bool) {
        Self::set_null(self, null);
    }
    fn write_cell(&self) -> Result<String> {
        Ok(self.to_cell_string())
    }
    fn read_cell(&mut self, text: &str) -> Result<()> {
        self.from_cell_string(text)
    }
}

/// A `spectra_ref` cell: `ms_run[n]:<native spectrum identifier>`.
///
/// `n` is the one-based index of an `ms_run` entry in the metadata section. The
/// cell is null when the run index is below one or the reference text is empty.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabSpectraRef {
    ms_run: usize,
    spec_ref: String,
}

impl MzTabSpectraRef {
    /// The null reference, same as [`Default`](Default::default).
    pub fn null() -> Self {
        Self::default()
    }
    /// A reference to `spec_ref` in `ms_run[ms_file]`. Native constructor.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `ms_file` is zero or `spec_ref` is empty,
    /// the two conditions the source asserts.
    pub fn new(ms_file: usize, spec_ref: &str) -> Result<Self> {
        let mut reference = Self::default();
        reference.set_ms_file(ms_file)?;
        reference.set_spec_ref(spec_ref)?;
        Ok(reference)
    }
    /// Set the one-based `ms_run` index.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for `0`. The source asserts `index >= 1` and, in
    /// a release build where the assertion is compiled out, silently ignores
    /// the call, leaving the previous index in place.
    pub fn set_ms_file(&mut self, index: usize) -> Result<()> {
        if index < 1 {
            return Err(bad("MzTab ms_run index is one-based and must not be zero"));
        }
        self.ms_run = index;
        Ok(())
    }
    /// Set the native spectrum identifier.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for empty text. The source asserts, then in a
    /// release build logs "Spectrum reference not set." and keeps the previous
    /// value; this returns the error instead of writing to a log.
    pub fn set_spec_ref(&mut self, spec_ref: &str) -> Result<()> {
        if spec_ref.is_empty() {
            return Err(bad("MzTab spectrum reference must not be empty"));
        }
        self.spec_ref.clear();
        self.spec_ref.push_str(spec_ref);
        Ok(())
    }
    /// Identical to [`MzTabSpectraRef::set_spec_ref`].
    ///
    /// The source's `setSpecRefFile` differs from `setSpecRef` only in not
    /// logging a warning for empty input, despite the name suggesting it sets a
    /// file rather than a spectrum. Kept so the source name resolves.
    pub fn set_spec_ref_file(&mut self, spec_ref: &str) -> Result<()> {
        self.set_spec_ref(spec_ref)
    }
    /// The native spectrum identifier, empty when the reference is null.
    ///
    /// The source's `getSpecRef` asserts `!isNull()`; a release build returns
    /// the empty string, as this does.
    pub fn spec_ref(&self) -> &str {
        &self.spec_ref
    }
    /// The one-based `ms_run` index, zero when the reference is null.
    ///
    /// The source's `getMSFile` asserts `!isNull()`; a release build returns
    /// the stored number, as this does.
    pub fn ms_file(&self) -> usize {
        self.ms_run
    }
    /// Both parts, or `None` when the reference is null. Native accessor.
    pub fn resolved(&self) -> Option<(usize, &str)> {
        if self.is_null() {
            None
        } else {
            Some((self.ms_run, self.spec_ref.as_str()))
        }
    }
    /// Whether the run index is below one or the reference text is empty.
    pub fn is_null(&self) -> bool {
        self.ms_run < 1 || self.spec_ref.is_empty()
    }
    /// Reset both parts when `null` is `true`; `false` is ignored.
    pub fn set_null(&mut self, null: bool) {
        if null {
            self.ms_run = 0;
            self.spec_ref.clear();
        }
    }
    /// Render as `null` or `ms_run[n]:<reference>`.
    pub fn to_cell_string(&self) -> String {
        if self.is_null() {
            "null".to_owned()
        } else {
            format!("ms_run[{}]:{}", self.ms_run, self.spec_ref)
        }
    }
    /// Replace the reference from MzTab text.
    ///
    /// The source splits on `:` and requires exactly two fields, so a native
    /// spectrum identifier that itself contains a colon cannot be read back.
    /// That restriction is preserved.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] when the text does not split into exactly two
    /// colon-separated fields, or when the run index is not a non-negative
    /// integer. The source casts a negative index to `Size`, producing a run
    /// index near `2^64`; this rejects it.
    pub fn from_cell_string(&mut self, text: &str) -> Result<()> {
        if is_token(text, "null") {
            self.set_null(true);
            return Ok(());
        }
        preflight_cell(text, ':')?;
        let fields = split_cells(text, ':');
        if fields.len() != 2 {
            return Err(conversion(format!(
                "can not convert {text:?} to an MzTab spectra_ref"
            )));
        }
        let index_text: String = fields[0].replace("ms_run[", "").replace(']', "");
        let index = trim_source(&index_text)
            .trim_start_matches('+')
            .parse::<i64>()
            .map_err(|_| {
                conversion(format!(
                    "can not convert {text:?} to an MzTab spectra_ref: {index_text:?} is not an integer"
                ))
            })?;
        if index < 0 {
            return Err(conversion(format!(
                "MzTab ms_run index must not be negative: {index_text:?}"
            )));
        }
        let ms_run = usize::try_from(index)
            .map_err(|_| conversion("MzTab ms_run index exceeds the address range"))?;
        self.ms_run = ms_run;
        self.spec_ref.clear();
        self.spec_ref.push_str(fields[1]);
        Ok(())
    }
}

impl MzTabCell for MzTabSpectraRef {
    fn is_null(&self) -> bool {
        Self::is_null(self)
    }
    fn set_null(&mut self, null: bool) {
        Self::set_null(self, null);
    }
    fn write_cell(&self) -> Result<String> {
        Ok(self.to_cell_string())
    }
    fn read_cell(&mut self, text: &str) -> Result<()> {
        self.from_cell_string(text)
    }
}

// ---------------------------------------------------------------------------
// MzTab.h — modification cells
// ---------------------------------------------------------------------------

/// One modification or substitution of a `modifications` cell.
///
/// Rendered `pos[param]|pos[param]-identifier`, or just `identifier` when no
/// position is known. Positions are one-based residue indices; `0` denotes the
/// N-terminus and `sequence length + 1` the C-terminus per the specification.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabModification {
    positions: Vec<(usize, MzTabParameter)>,
    identifier: MzTabString,
}

impl MzTabModification {
    /// The null modification, same as [`Default`](Default::default).
    pub fn null() -> Self {
        Self::default()
    }
    /// Set the possibly ambiguous position(s) and their associated parameter,
    /// which may be null when not set.
    pub fn set_positions_and_parameters(&mut self, positions: Vec<(usize, MzTabParameter)>) {
        self.positions = positions;
    }
    /// The position/parameter pairs. The source returns a copy of the vector.
    pub fn positions_and_parameters(&self) -> &[(usize, MzTabParameter)] {
        &self.positions
    }
    /// Set the modification or substitution identifier, for example
    /// `UNIMOD:35` or `CHEMMOD:15.9949`.
    pub fn set_modification_identifier(&mut self, identifier: MzTabString) {
        self.identifier = identifier;
    }
    /// The modification or substitution identifier.
    ///
    /// The source's `getModOrSubstIdentifier` asserts `!isNull()`; a release
    /// build returns the stored cell, as this does.
    pub fn mod_or_subst_identifier(&self) -> &MzTabString {
        &self.identifier
    }
    /// Whether there are no positions and the identifier is null.
    pub fn is_null(&self) -> bool {
        self.positions.is_empty() && self.identifier.is_null()
    }
    /// Clear positions and identifier when `null` is `true`; `false` is
    /// ignored.
    pub fn set_null(&mut self, null: bool) {
        if null {
            self.positions.clear();
            self.identifier.set_null(true);
        }
    }
    /// Render the modification.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when the modification is not null but the
    /// identifier is, mapping the source's `Exception::ConversionError`
    /// "Modification or Substitution identifier MUST NOT be null or empty".
    ///
    /// # Notes
    ///
    /// The source appends each position with `std::string::operator+=` applied
    /// to a `Size`, which resolves to the `char` overload: position `3` is
    /// written as the single byte `0x03`, not as the digit `3`. This port
    /// writes decimal digits, which is what the format requires and what the
    /// source's own parser reads back. Any position-annotated cell therefore
    /// differs from the bytes the C++ produces; see `docs/MZTAB_SUPPORT.md`.
    pub fn to_cell_string(&self) -> Result<String> {
        if self.is_null() {
            return Ok("null".to_owned());
        }
        let mut positions = String::new();
        for (index, (position, parameter)) in self.positions.iter().enumerate() {
            positions.push_str(&position.to_string());
            if !parameter.is_null() {
                positions.push_str(&parameter.to_cell_string());
            }
            if index + 1 < self.positions.len() {
                positions.push('|');
            }
        }
        if self.identifier.is_null() {
            return Err(Error::MissingInformation(
                "MzTab modification or substitution identifier must not be null".into(),
            ));
        }
        if positions.is_empty() {
            Ok(self.identifier.to_cell_string())
        } else {
            Ok(format!("{positions}-{}", self.identifier.to_cell_string()))
        }
    }
    /// Replace the modification from MzTab text.
    ///
    /// Text without `-` is taken as a bare identifier. Otherwise the text is
    /// split on `-` and must yield exactly two fields, the positions and the
    /// identifier; each `|`-separated position is an integer optionally
    /// followed by a bracketed parameter.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] when the text contains more than one `-`, or a position
    /// is not an integer, or a bracketed parameter does not parse, and
    /// [`Error::InvalidValue`] for the [`MAX_CELL_BYTES`]/[`MAX_CELL_ITEMS`]
    /// ceilings.
    ///
    /// # Notes
    ///
    /// Splitting on `-` means a `CHEMMOD` identifier with a negative mass
    /// delta, which the source itself writes for any modification without a
    /// UniMod accession, cannot be read back: `CHEMMOD:-18.010565` yields two
    /// fields and is misread as position `CHEMMOD:` with identifier
    /// `18.010565`, and `8-CHEMMOD:-18.010565` yields three and is rejected.
    /// The restriction is preserved so this parser accepts exactly what the
    /// source accepts.
    ///
    /// The parsed positions *replace* the previous ones; the source appends.
    pub fn from_cell_string(&mut self, text: &str) -> Result<()> {
        if is_token(text, "null") {
            self.set_null(true);
            return Ok(());
        }
        preflight_cell(text, '-')?;
        if !text.contains('-') {
            let mut next = Self::default();
            next.identifier.from_cell_string(trim_source(text));
            *self = next;
            return Ok(());
        }
        let trimmed = trim_source(text);
        let fields = split_cells(trimmed, '-');
        if fields.len() != 2 {
            return Err(conversion(format!(
                "can not convert {text:?} to an MzTab modification"
            )));
        }
        let mut next = Self::default();
        next.identifier.from_cell_string(trim_source(fields[1]));
        preflight_cell(fields[0], '|')?;
        for position_field in split_cells(fields[0], '|') {
            match position_field.find('[') {
                None => {
                    let position = parse_position(position_field)?;
                    next.positions.push((position, MzTabParameter::default()));
                }
                Some(offset) => {
                    // `find` returns a character boundary of this string, and
                    // '[' is single-byte, so both halves are valid slices.
                    let position = parse_position(&position_field[..offset])?;
                    let mut parameter = MzTabParameter::default();
                    parameter.from_cell_string(&position_field[offset..])?;
                    next.positions.push((position, parameter));
                }
            }
        }
        *self = next;
        Ok(())
    }
}

fn parse_position(text: &str) -> Result<usize> {
    let token = trim_source(text).trim_start_matches('+');
    let value = token.parse::<i64>().map_err(|_| {
        conversion(format!(
            "could not convert {token:?} to an MzTab modification position"
        ))
    })?;
    usize::try_from(value).map_err(|_| {
        conversion(format!(
            "MzTab modification position must not be negative: {token:?}"
        ))
    })
}

impl MzTabCell for MzTabModification {
    fn is_null(&self) -> bool {
        Self::is_null(self)
    }
    fn set_null(&mut self, null: bool) {
        Self::set_null(self, null);
    }
    fn write_cell(&self) -> Result<String> {
        self.to_cell_string()
    }
    fn read_cell(&mut self, text: &str) -> Result<()> {
        self.from_cell_string(text)
    }
}

/// A `modifications` cell: modifications separated by `,`. Null when empty.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabModificationList {
    entries: Vec<MzTabModification>,
}

impl MzTabModificationList {
    /// An empty, and therefore null, list.
    pub fn null() -> Self {
        Self::default()
    }
    /// The entries. The source returns a copy of the vector.
    pub fn get(&self) -> &[MzTabModification] {
        &self.entries
    }
    /// Mutable access to the entries. Native accessor.
    pub fn get_mut(&mut self) -> &mut Vec<MzTabModification> {
        &mut self.entries
    }
    /// Replace the entries.
    pub fn set(&mut self, entries: Vec<MzTabModification>) {
        self.entries = entries;
    }
    /// Whether the list is empty.
    pub fn is_null(&self) -> bool {
        self.entries.is_empty()
    }
    /// Clear the list when `null` is `true`; `false` is ignored.
    pub fn set_null(&mut self, null: bool) {
        if null {
            self.entries.clear();
        }
    }
    /// Render the list as `null` or `,`-separated modifications.
    ///
    /// # Errors
    ///
    /// As [`MzTabModification::to_cell_string`], for the first entry that
    /// carries positions but no identifier.
    pub fn to_cell_string(&self) -> Result<String> {
        if self.is_null() {
            return Ok("null".to_owned());
        }
        let mut rendered = Vec::with_capacity(self.entries.len());
        for entry in &self.entries {
            rendered.push(entry.to_cell_string()?);
        }
        Ok(rendered.join(","))
    }
    /// Replace the list from MzTab text.
    ///
    /// Text without `[` is split on every `,`. Otherwise a comma is a separator
    /// unless it is inside an unquoted `[…]`, so parameter commas do not split
    /// the list.
    ///
    /// # Errors
    ///
    /// As [`MzTabModification::from_cell_string`] for each field, plus
    /// [`Error::InvalidValue`] for the [`MAX_CELL_BYTES`]/[`MAX_CELL_ITEMS`]
    /// ceilings.
    ///
    /// # Notes
    ///
    /// A comma inside a quoted parameter part still splits, because the
    /// source's protection test requires the scanner to be outside quotes *and*
    /// inside a bracket. Its own worked example,
    /// `3|4[a,b,,v]|8[,,"blabla, [bla]",v],1|2|3[a,b,,v]-mod:123`, is therefore
    /// split inside the quoted text, contrary to the comment above the loop.
    /// Reproduced so this parser accepts exactly what the source accepts.
    ///
    /// The parsed entries *replace* the previous ones; the source appends.
    pub fn from_cell_string(&mut self, text: &str) -> Result<()> {
        if is_token(text, "null") {
            self.set_null(true);
            return Ok(());
        }
        preflight_cell(text, ',')?;
        let mut next = Vec::new();
        for field in split_modification_entries(text) {
            let mut entry = MzTabModification::default();
            entry.from_cell_string(&field)?;
            next.push(entry);
        }
        self.entries = next;
        Ok(())
    }
}

/// The source's bell-character trick, without mutating the subject: commas
/// inside an unquoted bracket are kept, all other commas separate. Brackets and
/// quotes stay in the field text, as they do in the source.
fn split_modification_entries(text: &str) -> Vec<String> {
    if !text.contains('[') {
        return split_cells(text, ',')
            .into_iter()
            .map(str::to_owned)
            .collect();
    }
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut in_bracket = false;
    let mut in_quotes = false;
    for c in text.chars() {
        match c {
            '[' if !in_quotes => {
                in_bracket = true;
                field.push('[');
            }
            ']' if !in_quotes => {
                in_bracket = false;
                field.push(']');
            }
            '"' => {
                in_quotes = !in_quotes;
                field.push('"');
            }
            ',' if !in_quotes && in_bracket => field.push(','),
            ',' => fields.push(std::mem::take(&mut field)),
            other => field.push(other),
        }
    }
    fields.push(field);
    fields
}

impl MzTabCell for MzTabModificationList {
    fn is_null(&self) -> bool {
        Self::is_null(self)
    }
    fn set_null(&mut self, null: bool) {
        Self::set_null(self, null);
    }
    fn write_cell(&self) -> Result<String> {
        self.to_cell_string()
    }
    fn read_cell(&mut self, text: &str) -> Result<()> {
        self.from_cell_string(text)
    }
}

// ---------------------------------------------------------------------------
// Metadata section records
// ---------------------------------------------------------------------------

/// `MTD software[n]`: the software and its ordered settings.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabSoftwareMetaData {
    /// The software term.
    pub software: MzTabParameter,
    /// `software[n]-setting[m]`, keyed by the one-based `m`.
    pub setting: BTreeMap<usize, MzTabString>,
}

/// `MTD sample[n]`: description and the CV terms describing the sample.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabSampleMetaData {
    /// `sample[n]-description`.
    pub description: MzTabString,
    /// `sample[n]-species[m]`.
    pub species: BTreeMap<usize, MzTabParameter>,
    /// `sample[n]-tissue[m]`.
    pub tissue: BTreeMap<usize, MzTabParameter>,
    /// `sample[n]-cell_type[m]`.
    pub cell_type: BTreeMap<usize, MzTabParameter>,
    /// `sample[n]-disease[m]`.
    pub disease: BTreeMap<usize, MzTabParameter>,
    /// `sample[n]-custom[m]`.
    pub custom: BTreeMap<usize, MzTabParameter>,
}

/// `MTD cv[n]`: one controlled vocabulary the file references.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabCVMetaData {
    /// `cv[n]-label`.
    pub label: MzTabString,
    /// `cv[n]-full_name`.
    pub full_name: MzTabString,
    /// `cv[n]-version`.
    pub version: MzTabString,
    /// `cv[n]-url`.
    pub url: MzTabString,
}

/// `MTD instrument[n]`: name, source, analyzers and detector.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabInstrumentMetaData {
    /// `instrument[n]-name`.
    pub name: MzTabParameter,
    /// `instrument[n]-source`.
    pub source: MzTabParameter,
    /// `instrument[n]-analyzer[m]`.
    pub analyzer: BTreeMap<usize, MzTabParameter>,
    /// `instrument[n]-detector`.
    pub detector: MzTabParameter,
}

/// `MTD contact[n]`: name, affiliation and e-mail.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabContactMetaData {
    /// `contact[n]-name`.
    pub name: MzTabString,
    /// `contact[n]-affiliation`.
    pub affiliation: MzTabString,
    /// `contact[n]-email`.
    pub email: MzTabString,
}

/// `MTD fixed_mod[n]` / `variable_mod[n]`: the modification, site and position.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabModificationMetaData {
    /// The modification term, by preference a UniMod accession.
    pub modification: MzTabParameter,
    /// The modified residue, `X` for a terminal modification with no residue.
    pub site: MzTabString,
    /// `Anywhere`, `Any N-term`, `Any C-term`, `Protein N-term` or
    /// `Protein C-term`.
    pub position: MzTabString,
}

/// `MTD assay[n]`: reagent, quantification modifications and references.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabAssayMetaData {
    /// `assay[n]-quantification_reagent`.
    pub quantification_reagent: MzTabParameter,
    /// `assay[n]-quantification_mod[m]`.
    pub quantification_mod: BTreeMap<usize, MzTabModificationMetaData>,
    /// `assay[n]-sample_ref`.
    pub sample_ref: MzTabString,
    /// `assay[n]-ms_run_ref`, a list rather than a single reference; the source
    /// comment records this as an adaptation for HUPO-PSI/mzTab issue 26.
    pub ms_run_ref: Vec<i32>,
}

/// `MTD ms_run[n]`: format, location, identifier format and fragmentation.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabMSRunMetaData {
    /// `ms_run[n]-format`.
    pub format: MzTabParameter,
    /// `ms_run[n]-location`.
    pub location: MzTabString,
    /// `ms_run[n]-id_format`, the native spectrum identifier format.
    pub id_format: MzTabParameter,
    /// `ms_run[n]-fragmentation_method`.
    pub fragmentation_method: MzTabParameterList,
}

/// `MTD study_variable[n]`: the assays and samples it groups.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabStudyVariableMetaData {
    /// `study_variable[n]-assay_refs`.
    pub assay_refs: Vec<i32>,
    /// `study_variable[n]-sample_refs`.
    pub sample_refs: Vec<i32>,
    /// `study_variable[n]-description`.
    pub description: MzTabString,
}

/// The whole `MTD` section. Refer to the MzTab specification for each key.
///
/// [`Default`](Default::default) sets `mz_tab_version` to `1.0.0`, as the
/// source's constructor; every other field starts null or empty. All indexed
/// keys use one-based indices, ordered by [`BTreeMap`] so a writer emits them
/// in index order rather than hash order.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabMetaData {
    /// `mzTab-version`. Defaults to `1.0.0`.
    pub mz_tab_version: MzTabString,
    /// `mzTab-mode`: `Complete` or `Summary`.
    pub mz_tab_mode: MzTabString,
    /// `mzTab-type`: `Identification` or `Quantification`.
    pub mz_tab_type: MzTabString,
    /// `mzTab-ID`.
    pub mz_tab_id: MzTabString,
    /// `title`.
    pub title: MzTabString,
    /// `description`.
    pub description: MzTabString,
    /// `protein_search_engine_score[n]`.
    pub protein_search_engine_score: BTreeMap<usize, MzTabParameter>,
    /// `peptide_search_engine_score[n]`.
    pub peptide_search_engine_score: BTreeMap<usize, MzTabParameter>,
    /// `psm_search_engine_score[n]`.
    pub psm_search_engine_score: BTreeMap<usize, MzTabParameter>,
    /// `smallmolecule_search_engine_score[n]`.
    pub smallmolecule_search_engine_score: BTreeMap<usize, MzTabParameter>,
    /// `nucleic_acid_search_engine_score[n]`, an OpenMS extension.
    pub nucleic_acid_search_engine_score: BTreeMap<usize, MzTabParameter>,
    /// `oligonucleotide_search_engine_score[n]`, an OpenMS extension.
    pub oligonucleotide_search_engine_score: BTreeMap<usize, MzTabParameter>,
    /// `osm_search_engine_score[n]`, an OpenMS extension.
    pub osm_search_engine_score: BTreeMap<usize, MzTabParameter>,
    /// `sample_processing[n]`.
    pub sample_processing: BTreeMap<usize, MzTabParameterList>,
    /// `instrument[n]-…`.
    pub instrument: BTreeMap<usize, MzTabInstrumentMetaData>,
    /// `software[n]-…`.
    pub software: BTreeMap<usize, MzTabSoftwareMetaData>,
    /// `false_discovery_rate`.
    pub false_discovery_rate: MzTabParameterList,
    /// `publication[n]`.
    pub publication: BTreeMap<usize, MzTabString>,
    /// `contact[n]-…`.
    pub contact: BTreeMap<usize, MzTabContactMetaData>,
    /// `uri[n]`.
    pub uri: BTreeMap<usize, MzTabString>,
    /// `fixed_mod[n]-…`.
    pub fixed_mod: BTreeMap<usize, MzTabModificationMetaData>,
    /// `variable_mod[n]-…`.
    pub variable_mod: BTreeMap<usize, MzTabModificationMetaData>,
    /// `quantification_method`.
    pub quantification_method: MzTabParameter,
    /// `protein-quantification_unit`.
    pub protein_quantification_unit: MzTabParameter,
    /// `peptide-quantification_unit`.
    pub peptide_quantification_unit: MzTabParameter,
    /// `small_molecule-quantification_unit`.
    pub small_molecule_quantification_unit: MzTabParameter,
    /// `ms_run[n]-…`.
    pub ms_run: BTreeMap<usize, MzTabMSRunMetaData>,
    /// `custom[n]`.
    pub custom: BTreeMap<usize, MzTabParameter>,
    /// `sample[n]-…`.
    pub sample: BTreeMap<usize, MzTabSampleMetaData>,
    /// `assay[n]-…`.
    pub assay: BTreeMap<usize, MzTabAssayMetaData>,
    /// `study_variable[n]-…`.
    pub study_variable: BTreeMap<usize, MzTabStudyVariableMetaData>,
    /// `cv[n]-…`.
    pub cv: BTreeMap<usize, MzTabCVMetaData>,
    /// `colunit-protein` declarations, verbatim.
    pub colunit_protein: Vec<String>,
    /// `colunit-peptide` declarations, verbatim.
    pub colunit_peptide: Vec<String>,
    /// `colunit-psm` declarations, verbatim.
    pub colunit_psm: Vec<String>,
    /// `colunit-small_molecule` declarations, verbatim.
    pub colunit_small_molecule: Vec<String>,
}

impl Default for MzTabMetaData {
    fn default() -> Self {
        Self {
            mz_tab_version: MzTabString::from_text("1.0.0"),
            mz_tab_mode: MzTabString::default(),
            mz_tab_type: MzTabString::default(),
            mz_tab_id: MzTabString::default(),
            title: MzTabString::default(),
            description: MzTabString::default(),
            protein_search_engine_score: BTreeMap::new(),
            peptide_search_engine_score: BTreeMap::new(),
            psm_search_engine_score: BTreeMap::new(),
            smallmolecule_search_engine_score: BTreeMap::new(),
            nucleic_acid_search_engine_score: BTreeMap::new(),
            oligonucleotide_search_engine_score: BTreeMap::new(),
            osm_search_engine_score: BTreeMap::new(),
            sample_processing: BTreeMap::new(),
            instrument: BTreeMap::new(),
            software: BTreeMap::new(),
            false_discovery_rate: MzTabParameterList::default(),
            publication: BTreeMap::new(),
            contact: BTreeMap::new(),
            uri: BTreeMap::new(),
            fixed_mod: BTreeMap::new(),
            variable_mod: BTreeMap::new(),
            quantification_method: MzTabParameter::default(),
            protein_quantification_unit: MzTabParameter::default(),
            peptide_quantification_unit: MzTabParameter::default(),
            small_molecule_quantification_unit: MzTabParameter::default(),
            ms_run: BTreeMap::new(),
            custom: BTreeMap::new(),
            sample: BTreeMap::new(),
            assay: BTreeMap::new(),
            study_variable: BTreeMap::new(),
            cv: BTreeMap::new(),
            colunit_protein: Vec::new(),
            colunit_peptide: Vec::new(),
            colunit_psm: Vec::new(),
            colunit_small_molecule: Vec::new(),
        }
    }
}

// ---------------------------------------------------------------------------
// Section rows
// ---------------------------------------------------------------------------

/// A section row that carries optional `opt_…` columns.
///
/// Replaces the source's protected template `MzTabBase::getOptionalColumnNames_`,
/// whose only requirement on its argument is a public `opt_` member.
pub trait MzTabOptionalColumns {
    /// The row's optional columns, in the order the row declares them.
    fn optional_columns(&self) -> &[MzTabOptionalColumnEntry];
    /// Mutable access to the row's optional columns.
    fn optional_columns_mut(&mut self) -> &mut Vec<MzTabOptionalColumnEntry>;
}

/// Collect the optional column names of `rows`, first occurrence first.
///
/// The source's `getOptionalColumnNames_`: a vector rather than a set, so the
/// column order of the first row that declares a name is preserved, and a name
/// seen again is skipped.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `rows` exceeds [`MzTab::MAX_ROWS`] or the
/// distinct names exceed [`MzTab::MAX_OPTIONAL_COLUMNS`]; the row count is
/// checked before anything is allocated.
///
/// # Notes
///
/// The source deduplicates with a linear `std::find` over the names collected
/// so far, which is quadratic in the number of distinct columns. This uses an
/// auxiliary ordered set for the membership test and produces the same
/// sequence.
pub fn optional_column_names<R: MzTabOptionalColumns>(rows: &[R]) -> Result<Vec<String>> {
    if rows.len() > MzTab::MAX_ROWS {
        return Err(bad("MzTab section exceeds its row limit"));
    }
    let mut names = Vec::new();
    let mut seen = BTreeSet::new();
    for row in rows {
        for entry in row.optional_columns() {
            if seen.insert(entry.name.as_str()) {
                if names.len() >= MzTab::MAX_OPTIONAL_COLUMNS {
                    return Err(bad("MzTab section exceeds its optional-column limit"));
                }
                names.push(entry.name.clone());
            }
        }
    }
    Ok(names)
}

macro_rules! optional_columns_for {
    ($row:ty) => {
        impl MzTabOptionalColumns for $row {
            fn optional_columns(&self) -> &[MzTabOptionalColumnEntry] {
                &self.opt
            }
            fn optional_columns_mut(&mut self) -> &mut Vec<MzTabOptionalColumnEntry> {
                &mut self.opt
            }
        }
    };
}

/// `PRT` — one protein row.
///
/// [`Default`](Default::default) sets the `go_terms` and `ambiguity_members`
/// separators to `,`, as the source's constructor, because `|` occurs inside GO
/// terms and protein accessions.
#[derive(Clone, Debug, PartialEq)]
pub struct MzTabProteinSectionRow {
    /// The protein's accession.
    pub accession: MzTabString,
    /// Human readable description, i.e. the name.
    pub description: MzTabString,
    /// NEWT taxonomy for the species.
    pub taxid: MzTabInteger,
    /// Human readable name of the species.
    pub species: MzTabString,
    /// Name of the protein database.
    pub database: MzTabString,
    /// Version of the protein database.
    pub database_version: MzTabString,
    /// Search engine(s) identifying the protein.
    pub search_engine: MzTabParameterList,
    /// `best_search_engine_score[1-n]`.
    pub best_search_engine_score: BTreeMap<usize, MzTabDouble>,
    /// `search_engine_score[index1]_ms_run[index2]`.
    pub search_engine_score_ms_run: BTreeMap<usize, BTreeMap<usize, MzTabDouble>>,
    /// Identification reliability, 1-3.
    pub reliability: MzTabInteger,
    /// `num_psms_ms_run[n]`.
    pub num_psms_ms_run: BTreeMap<usize, MzTabInteger>,
    /// `num_peptides_distinct_ms_run[n]`.
    pub num_peptides_distinct_ms_run: BTreeMap<usize, MzTabInteger>,
    /// `num_peptides_unique_ms_run[n]`.
    pub num_peptides_unique_ms_run: BTreeMap<usize, MzTabInteger>,
    /// Alternative protein identifications.
    pub ambiguity_members: MzTabStringList,
    /// Modifications identified in the protein.
    pub modifications: MzTabModificationList,
    /// Location of the protein's source entry.
    pub uri: MzTabString,
    /// List of GO terms for the protein.
    pub go_terms: MzTabStringList,
    /// Amount of protein sequence identified, 0-1.
    pub coverage: MzTabDouble,
    /// `protein_abundance_assay[n]`.
    pub protein_abundance_assay: BTreeMap<usize, MzTabDouble>,
    /// `protein_abundance_study_variable[n]`.
    pub protein_abundance_study_variable: BTreeMap<usize, MzTabDouble>,
    /// `protein_abundance_stdev_study_variable[n]`.
    pub protein_abundance_stdev_study_variable: BTreeMap<usize, MzTabDouble>,
    /// `protein_abundance_std_error_study_variable[n]`.
    pub protein_abundance_std_error_study_variable: BTreeMap<usize, MzTabDouble>,
    /// Optional columns, whose names must start with `opt_`. Source `opt_`.
    pub opt: Vec<MzTabOptionalColumnEntry>,
}

impl Default for MzTabProteinSectionRow {
    fn default() -> Self {
        let mut go_terms = MzTabStringList::default();
        go_terms.set_separator(',');
        let mut ambiguity_members = MzTabStringList::default();
        ambiguity_members.set_separator(',');
        Self {
            accession: MzTabString::default(),
            description: MzTabString::default(),
            taxid: MzTabInteger::default(),
            species: MzTabString::default(),
            database: MzTabString::default(),
            database_version: MzTabString::default(),
            search_engine: MzTabParameterList::default(),
            best_search_engine_score: BTreeMap::new(),
            search_engine_score_ms_run: BTreeMap::new(),
            reliability: MzTabInteger::default(),
            num_psms_ms_run: BTreeMap::new(),
            num_peptides_distinct_ms_run: BTreeMap::new(),
            num_peptides_unique_ms_run: BTreeMap::new(),
            ambiguity_members,
            modifications: MzTabModificationList::default(),
            uri: MzTabString::default(),
            go_terms,
            coverage: MzTabDouble::default(),
            protein_abundance_assay: BTreeMap::new(),
            protein_abundance_study_variable: BTreeMap::new(),
            protein_abundance_stdev_study_variable: BTreeMap::new(),
            protein_abundance_std_error_study_variable: BTreeMap::new(),
            opt: Vec::new(),
        }
    }
}

impl MzTabProteinSectionRow {
    /// Source `RowCompare`: orders rows by accession text alone.
    pub fn row_order(&self, other: &Self) -> Ordering {
        self.accession.get().cmp(other.accession.get())
    }
}

optional_columns_for!(MzTabProteinSectionRow);

/// `PEP` — one peptide row.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MzTabPeptideSectionRow {
    /// The peptide's sequence.
    pub sequence: MzTabString,
    /// The protein's accession.
    pub accession: MzTabString,
    /// `0` false, `1` true, null otherwise: peptide is unique for the protein.
    pub unique: MzTabBoolean,
    /// Name of the sequence database.
    pub database: MzTabString,
    /// Version, optionally with the number of entries.
    pub database_version: MzTabString,
    /// Search engine(s) that identified the peptide.
    pub search_engine: MzTabParameterList,
    /// Search engine score(s) for the peptide.
    pub best_search_engine_score: BTreeMap<usize, MzTabDouble>,
    /// `search_engine_score[index1]_ms_run[index2]`.
    pub search_engine_score_ms_run: BTreeMap<usize, BTreeMap<usize, MzTabDouble>>,
    /// Identification reliability for the peptide, 1-3; `0` is null.
    pub reliability: MzTabInteger,
    /// Modifications identified in the peptide.
    pub modifications: MzTabModificationList,
    /// Time points in seconds. Semantics may vary.
    pub retention_time: MzTabDoubleList,
    /// Retention-time window, in seconds.
    pub retention_time_window: MzTabDoubleList,
    /// Precursor ion's charge.
    pub charge: MzTabInteger,
    /// Precursor ion's m/z.
    pub mass_to_charge: MzTabDouble,
    /// Location of the PSM's source entry.
    pub uri: MzTabString,
    /// Spectra identifying the peptide.
    pub spectra_ref: MzTabSpectraRef,
    /// `peptide_abundance_assay[n]`.
    pub peptide_abundance_assay: BTreeMap<usize, MzTabDouble>,
    /// `peptide_abundance_study_variable[n]`.
    pub peptide_abundance_study_variable: BTreeMap<usize, MzTabDouble>,
    /// `peptide_abundance_stdev_study_variable[n]`.
    pub peptide_abundance_stdev_study_variable: BTreeMap<usize, MzTabDouble>,
    /// `peptide_abundance_std_error_study_variable[n]`.
    pub peptide_abundance_std_error_study_variable: BTreeMap<usize, MzTabDouble>,
    /// Optional columns, whose names must start with `opt_`. Source `opt_`.
    pub opt: Vec<MzTabOptionalColumnEntry>,
}

impl MzTabPeptideSectionRow {
    /// Source `RowCompare`: orders rows by sequence, then accession.
    pub fn row_order(&self, other: &Self) -> Ordering {
        (self.sequence.get(), self.accession.get())
            .cmp(&(other.sequence.get(), other.accession.get()))
    }
}

optional_columns_for!(MzTabPeptideSectionRow);

/// `PSM` — one peptide-spectrum match row.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MzTabPSMSectionRow {
    /// The peptide's sequence.
    pub sequence: MzTabString,
    /// A unique ID of a PSM line.
    pub psm_id: MzTabInteger,
    /// List of potential parent protein accessions, as in the FASTA database.
    pub accession: MzTabString,
    /// `0` false, `1` true, null otherwise: peptide is unique for the protein.
    pub unique: MzTabBoolean,
    /// Name of the sequence database.
    pub database: MzTabString,
    /// Version, optionally with the number of entries.
    pub database_version: MzTabString,
    /// Search engine(s) that identified the peptide.
    pub search_engine: MzTabParameterList,
    /// Search engine score(s) for the peptide.
    pub search_engine_score: BTreeMap<usize, MzTabDouble>,
    /// Identification reliability for the peptide, 1-3; `0` is null.
    pub reliability: MzTabInteger,
    /// Modifications identified in the peptide.
    pub modifications: MzTabModificationList,
    /// Time points in seconds. Semantics may vary.
    pub retention_time: MzTabDoubleList,
    /// The charge of the experimental precursor ion.
    pub charge: MzTabInteger,
    /// Observed m/z of the experimental precursor ion, raw or corrected.
    pub exp_mass_to_charge: MzTabDouble,
    /// Calculated m/z of the experimental precursor ion.
    pub calc_mass_to_charge: MzTabDouble,
    /// Location of the PSM's source entry.
    pub uri: MzTabString,
    /// Spectrum for this PSM.
    pub spectra_ref: MzTabSpectraRef,
    /// Amino acid(s) in parent protein(s) before the start of this PSM.
    pub pre: MzTabString,
    /// Amino acid(s) in parent protein(s) after the end of this PSM.
    pub post: MzTabString,
    /// Start position(s) in parent protein(s).
    pub start: MzTabString,
    /// End position(s) in parent protein(s).
    pub end: MzTabString,
    /// Optional columns, whose names must start with `opt_`. Source `opt_`.
    pub opt: Vec<MzTabOptionalColumnEntry>,
}

impl MzTabPSMSectionRow {
    /// Fill `pre`, `post`, `start`, `end` and `accession` from
    /// `peptide_evidences`, overwriting whatever those five cells held.
    ///
    /// Each cell becomes the comma-joined list of the per-evidence values, in
    /// evidence order, so the five lists are positionally aligned. Per the
    /// specification an unknown flanking residue is the text `null` and a
    /// terminal one is `-`; an unknown position is the text `null` and a known
    /// one is the protein coordinate **plus one**, because MzTab counts from
    /// one while `PeptideEvidence` counts from zero.
    ///
    /// An empty `peptide_evidences` clears `pre`, `post`, `start` and `end` —
    /// and leaves `accession` alone, as the source's early return does.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `peptide_evidences` exceeds
    /// [`MAX_CELL_ITEMS`], or when a start or end coordinate is at
    /// `usize::MAX` and so cannot be incremented. Nothing is written on error.
    pub fn add_pep_evidence_to_rows(
        &mut self,
        peptide_evidences: &[PeptideEvidence],
    ) -> Result<()> {
        if peptide_evidences.is_empty() {
            self.pre = MzTabString::default();
            self.post = MzTabString::default();
            self.start = MzTabString::default();
            self.end = MzTabString::default();
            return Ok(());
        }
        if peptide_evidences.len() > MAX_CELL_ITEMS {
            return Err(bad("MzTab PSM evidence list exceeds its entry limit"));
        }
        let residue = |flank: FlankingResidue| match flank {
            FlankingResidue::Unknown => "null".to_owned(),
            FlankingResidue::NTerminus | FlankingResidue::CTerminus => "-".to_owned(),
            FlankingResidue::Residue(code) => code.to_string(),
        };
        let position = |value: Option<usize>| -> Result<String> {
            match value {
                None => Ok("null".to_owned()),
                Some(value) => value
                    .checked_add(1)
                    .map(|one_based| one_based.to_string())
                    .ok_or_else(|| bad("MzTab peptide position cannot be made one-based")),
            }
        };
        let mut pre = Vec::with_capacity(peptide_evidences.len());
        let mut post = Vec::with_capacity(peptide_evidences.len());
        let mut start = Vec::with_capacity(peptide_evidences.len());
        let mut end = Vec::with_capacity(peptide_evidences.len());
        let mut accession = Vec::with_capacity(peptide_evidences.len());
        for evidence in peptide_evidences {
            pre.push(residue(evidence.aa_before));
            post.push(residue(evidence.aa_after));
            start.push(position(evidence.start)?);
            end.push(position(evidence.end)?);
            accession.push(evidence.protein_accession.clone());
        }
        self.pre = MzTabString::from_text(&pre.join(","));
        self.post = MzTabString::from_text(&post.join(","));
        self.start = MzTabString::from_text(&start.join(","));
        self.end = MzTabString::from_text(&end.join(","));
        self.accession = MzTabString::from_text(&accession.join(","));
        Ok(())
    }
    /// Source `RowCompare`: orders rows by sequence, then the reference's run
    /// index and native identifier, then accession. The source's own `@TODO`
    /// asks whether `PSM_ID` should take part; it does not.
    pub fn row_order(&self, other: &Self) -> Ordering {
        (
            self.sequence.get(),
            self.spectra_ref.ms_file(),
            self.spectra_ref.spec_ref(),
            self.accession.get(),
        )
            .cmp(&(
                other.sequence.get(),
                other.spectra_ref.ms_file(),
                other.spectra_ref.spec_ref(),
                other.accession.get(),
            ))
    }
}

optional_columns_for!(MzTabPSMSectionRow);

/// `SML` — one small-molecule row.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MzTabSmallMoleculeSectionRow {
    /// The small molecule's identifier.
    pub identifier: MzTabStringList,
    /// Chemical formula of the identified compound.
    pub chemical_formula: MzTabString,
    /// Molecular structure in SMILES format.
    pub smiles: MzTabString,
    /// InChI key of the identified compound.
    pub inchi_key: MzTabString,
    /// Human readable description, i.e. the name.
    pub description: MzTabString,
    /// Precursor ion's m/z, as measured.
    pub exp_mass_to_charge: MzTabDouble,
    /// Precursor ion's m/z, as calculated.
    pub calc_mass_to_charge: MzTabDouble,
    /// Precursor ion's charge.
    pub charge: MzTabInteger,
    /// Time points in seconds. Semantics may vary.
    pub retention_time: MzTabDoubleList,
    /// NEWT taxonomy for the species.
    pub taxid: MzTabInteger,
    /// Human readable name of the species.
    pub species: MzTabString,
    /// Name of the used database.
    pub database: MzTabString,
    /// Version of the database, optionally with the number of compounds.
    pub database_version: MzTabString,
    /// The identification reliability, 1-3.
    pub reliability: MzTabInteger,
    /// The source entry's location.
    pub uri: MzTabString,
    /// Spectra identifying the small molecule.
    pub spectra_ref: MzTabSpectraRef,
    /// Search engine(s) identifying the small molecule.
    pub search_engine: MzTabParameterList,
    /// Search engine identification score(s).
    pub best_search_engine_score: BTreeMap<usize, MzTabDouble>,
    /// `search_engine_score[index1]_ms_run[index2]`.
    pub search_engine_score_ms_run: BTreeMap<usize, BTreeMap<usize, MzTabDouble>>,
    /// Modifications identified on the small molecule. A plain string cell in
    /// this section, not a [`MzTabModificationList`].
    pub modifications: MzTabString,
    /// `smallmolecule_abundance_assay[n]`.
    pub smallmolecule_abundance_assay: BTreeMap<usize, MzTabDouble>,
    /// `smallmolecule_abundance_study_variable[n]`.
    pub smallmolecule_abundance_study_variable: BTreeMap<usize, MzTabDouble>,
    /// `smallmolecule_abundance_stdev_study_variable[n]`.
    pub smallmolecule_abundance_stdev_study_variable: BTreeMap<usize, MzTabDouble>,
    /// `smallmolecule_abundance_std_error_study_variable[n]`.
    pub smallmolecule_abundance_std_error_study_variable: BTreeMap<usize, MzTabDouble>,
    /// Optional columns, whose names must start with `opt_`. Source `opt_`.
    pub opt: Vec<MzTabOptionalColumnEntry>,
}

optional_columns_for!(MzTabSmallMoleculeSectionRow);

/// `NUC` — one nucleic-acid row, an OpenMS extension.
///
/// Unlike [`MzTabProteinSectionRow`], the source declares no constructor here,
/// so `ambiguity_members` and `go_terms` keep the default `|` separator rather
/// than the `,` the protein section uses. Preserved.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MzTabNucleicAcidSectionRow {
    /// The nucleic acid's accession.
    pub accession: MzTabString,
    /// Human readable description, i.e. the name.
    pub description: MzTabString,
    /// NEWT taxonomy for the species.
    pub taxid: MzTabInteger,
    /// Human readable name of the species.
    pub species: MzTabString,
    /// Name of the sequence database.
    pub database: MzTabString,
    /// Version of the sequence database.
    pub database_version: MzTabString,
    /// Search engine(s) that identified the nucleic acid.
    pub search_engine: MzTabParameterList,
    /// Best search engine score(s) over all MS runs.
    pub best_search_engine_score: BTreeMap<usize, MzTabDouble>,
    /// `search_engine_score[index1]_ms_run[index2]`.
    pub search_engine_score_ms_run: BTreeMap<usize, BTreeMap<usize, MzTabDouble>>,
    /// Identification reliability, 1-3.
    pub reliability: MzTabInteger,
    /// `num_osms_ms_run[n]`.
    pub num_osms_ms_run: BTreeMap<usize, MzTabInteger>,
    /// `num_oligos_distinct_ms_run[n]`.
    pub num_oligos_distinct_ms_run: BTreeMap<usize, MzTabInteger>,
    /// `num_oligos_unique_ms_run[n]`.
    pub num_oligos_unique_ms_run: BTreeMap<usize, MzTabInteger>,
    /// Alternative nucleic acid identifications.
    pub ambiguity_members: MzTabStringList,
    /// Modifications identified in the nucleic acid.
    pub modifications: MzTabModificationList,
    /// Location of the nucleic acid's source entry.
    pub uri: MzTabString,
    /// List of GO terms. The source asks in a comment whether these make sense
    /// for nucleic acid sequences; the column exists regardless.
    pub go_terms: MzTabStringList,
    /// Fraction of nucleic acid sequence identified, 0-1.
    pub coverage: MzTabDouble,
    /// Optional columns, whose names must start with `opt_`. Source `opt_`.
    pub opt: Vec<MzTabOptionalColumnEntry>,
}

impl MzTabNucleicAcidSectionRow {
    /// Source `RowCompare`: orders rows by accession text alone.
    pub fn row_order(&self, other: &Self) -> Ordering {
        self.accession.get().cmp(other.accession.get())
    }
}

optional_columns_for!(MzTabNucleicAcidSectionRow);

/// `OLI` — one oligonucleotide row, an OpenMS extension.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MzTabOligonucleotideSectionRow {
    /// The oligonucleotide's sequence.
    pub sequence: MzTabString,
    /// The nucleic acid's accession.
    pub accession: MzTabString,
    /// `0` false, `1` true, null otherwise: the oligonucleotide maps uniquely.
    pub unique: MzTabBoolean,
    /// Search engine(s) that identified the match.
    pub search_engine: MzTabParameterList,
    /// Search engine score(s) for the match.
    pub best_search_engine_score: BTreeMap<usize, MzTabDouble>,
    /// Search engine score(s) per individual MS run.
    pub search_engine_score_ms_run: BTreeMap<usize, BTreeMap<usize, MzTabDouble>>,
    /// Identification reliability for the match, 1-3; `0` is null.
    pub reliability: MzTabInteger,
    /// Modifications identified in the oligonucleotide.
    pub modifications: MzTabModificationList,
    /// Time points in seconds. Semantics may vary.
    pub retention_time: MzTabDoubleList,
    /// Retention-time window, in seconds.
    pub retention_time_window: MzTabDoubleList,
    /// Location of the oligonucleotide's source entry.
    pub uri: MzTabString,
    /// Nucleotide in the parent sequence before the start of this match.
    pub pre: MzTabString,
    /// Nucleotide in the parent sequence after the end of this match.
    pub post: MzTabString,
    /// Start position in the parent sequence.
    pub start: MzTabInteger,
    /// End position in the parent sequence.
    pub end: MzTabInteger,
    /// Optional columns, whose names must start with `opt_`. Source `opt_`.
    pub opt: Vec<MzTabOptionalColumnEntry>,
}

impl MzTabOligonucleotideSectionRow {
    /// Source `RowCompare`: orders rows by sequence, accession, start and end.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when any of the four `start`/`end` cells
    /// involved is not a value. The source's comparator calls
    /// `MzTabInteger::get()`, which throws `Exception::ElementNotFound` for a
    /// null cell, so sorting a section with unset positions throws there too;
    /// returning the error keeps that failure visible instead of inventing an
    /// order.
    pub fn row_order(&self, other: &Self) -> Result<Ordering> {
        let (self_start, self_end) = (self.start.get()?, self.end.get()?);
        let (other_start, other_end) = (other.start.get()?, other.end.get()?);
        Ok((
            self.sequence.get(),
            self.accession.get(),
            self_start,
            self_end,
        )
            .cmp(&(
                other.sequence.get(),
                other.accession.get(),
                other_start,
                other_end,
            )))
    }
}

optional_columns_for!(MzTabOligonucleotideSectionRow);

/// `OSM` — one oligonucleotide-spectrum match row, an OpenMS extension.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MzTabOSMSectionRow {
    /// The oligonucleotide's sequence.
    pub sequence: MzTabString,
    /// Search engine(s) that identified the match.
    pub search_engine: MzTabParameterList,
    /// Search engine score(s) for the match.
    pub search_engine_score: BTreeMap<usize, MzTabDouble>,
    /// Identification reliability for the match, 1-3; `0` is null.
    pub reliability: MzTabInteger,
    /// Modifications identified in the oligonucleotide.
    pub modifications: MzTabModificationList,
    /// Time points in seconds. Semantics may vary.
    pub retention_time: MzTabDoubleList,
    /// The charge of the experimental precursor ion.
    pub charge: MzTabInteger,
    /// The m/z of the experimental precursor ion.
    pub exp_mass_to_charge: MzTabDouble,
    /// The theoretical m/z of the oligonucleotide.
    pub calc_mass_to_charge: MzTabDouble,
    /// Location of the OSM's source entry.
    pub uri: MzTabString,
    /// Reference to the spectrum underlying the match.
    pub spectra_ref: MzTabSpectraRef,
    /// Optional columns, whose names must start with `opt_`. Source `opt_`.
    pub opt: Vec<MzTabOptionalColumnEntry>,
}

impl MzTabOSMSectionRow {
    /// Source `RowCompare`: orders rows by sequence, then the reference's run
    /// index and native identifier.
    pub fn row_order(&self, other: &Self) -> Ordering {
        (
            self.sequence.get(),
            self.spectra_ref.ms_file(),
            self.spectra_ref.spec_ref(),
        )
            .cmp(&(
                other.sequence.get(),
                other.spectra_ref.ms_file(),
                other.spectra_ref.spec_ref(),
            ))
    }
}

optional_columns_for!(MzTabOSMSectionRow);

/// The `PRT` section. Source `MzTabProteinSectionRows`.
pub type MzTabProteinSectionRows = Vec<MzTabProteinSectionRow>;
/// The `PEP` section. Source `MzTabPeptideSectionRows`.
pub type MzTabPeptideSectionRows = Vec<MzTabPeptideSectionRow>;
/// The `PSM` section. Source `MzTabPSMSectionRows`.
pub type MzTabPSMSectionRows = Vec<MzTabPSMSectionRow>;
/// The `SML` section. Source `MzTabSmallMoleculeSectionRows`.
pub type MzTabSmallMoleculeSectionRows = Vec<MzTabSmallMoleculeSectionRow>;
/// The `NUC` section. Source `MzTabNucleicAcidSectionRows`.
pub type MzTabNucleicAcidSectionRows = Vec<MzTabNucleicAcidSectionRow>;
/// The `OLI` section. Source `MzTabOligonucleotideSectionRows`.
pub type MzTabOligonucleotideSectionRows = Vec<MzTabOligonucleotideSectionRow>;
/// The `OSM` section. Source `MzTabOSMSectionRows`.
pub type MzTabOSMSectionRows = Vec<MzTabOSMSectionRow>;

// ---------------------------------------------------------------------------
// The document
// ---------------------------------------------------------------------------

/// A whole MzTab document: the metadata section, the seven data sections, and
/// the comment and empty lines recorded by position.
///
/// The source keeps the sections in protected members behind getter/setter
/// pairs that do nothing but return a reference, so they are public fields
/// here; the getter names remain available as borrowing methods.
/// [`Default`](Default::default) is the source's default constructor, whose
/// metadata already declares `mzTab-version 1.0.0`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MzTab {
    /// The `MTD` section.
    pub meta_data: MzTabMetaData,
    /// The `PRT` section.
    pub protein_data: MzTabProteinSectionRows,
    /// The `PEP` section.
    pub peptide_data: MzTabPeptideSectionRows,
    /// The `PSM` section.
    pub psm_data: MzTabPSMSectionRows,
    /// The `SML` section.
    pub small_molecule_data: MzTabSmallMoleculeSectionRows,
    /// The `NUC` section.
    pub nucleic_acid_data: MzTabNucleicAcidSectionRows,
    /// The `OLI` section.
    pub oligonucleotide_data: MzTabOligonucleotideSectionRows,
    /// The `OSM` section, oligonucleotide-spectrum matches.
    pub osm_data: MzTabOSMSectionRows,
    /// Line indices of empty rows, preserved so a round trip can restore them.
    pub empty_rows: Vec<usize>,
    /// `COM` comment lines, keyed by line index.
    pub comment_rows: BTreeMap<usize, String>,
}

impl MzTab {
    /// Largest number of rows a section operation will traverse.
    pub const MAX_ROWS: usize = 10_000_000;
    /// Largest number of distinct optional columns one section may declare.
    pub const MAX_OPTIONAL_COLUMNS: usize = 100_000;

    /// The `MTD` section. Source `getMetaData`.
    pub fn meta_data(&self) -> &MzTabMetaData {
        &self.meta_data
    }
    /// Replace the `MTD` section. Source `setMetaData`.
    pub fn set_meta_data(&mut self, meta_data: MzTabMetaData) {
        self.meta_data = meta_data;
    }
    /// The `PRT` section. Source `getProteinSectionRows`.
    pub fn protein_section_rows(&self) -> &MzTabProteinSectionRows {
        &self.protein_data
    }
    /// Replace the `PRT` section. Source `setProteinSectionRows`.
    pub fn set_protein_section_rows(&mut self, rows: MzTabProteinSectionRows) {
        self.protein_data = rows;
    }
    /// The `PEP` section. Source `getPeptideSectionRows`.
    pub fn peptide_section_rows(&self) -> &MzTabPeptideSectionRows {
        &self.peptide_data
    }
    /// Replace the `PEP` section. Source `setPeptideSectionRows`.
    pub fn set_peptide_section_rows(&mut self, rows: MzTabPeptideSectionRows) {
        self.peptide_data = rows;
    }
    /// The `PSM` section. Source `getPSMSectionRows`.
    pub fn psm_section_rows(&self) -> &MzTabPSMSectionRows {
        &self.psm_data
    }
    /// Replace the `PSM` section. Source `setPSMSectionRows`.
    pub fn set_psm_section_rows(&mut self, rows: MzTabPSMSectionRows) {
        self.psm_data = rows;
    }
    /// The `SML` section. Source `getSmallMoleculeSectionRows`.
    pub fn small_molecule_section_rows(&self) -> &MzTabSmallMoleculeSectionRows {
        &self.small_molecule_data
    }
    /// Replace the `SML` section. Source `setSmallMoleculeSectionRows`.
    pub fn set_small_molecule_section_rows(&mut self, rows: MzTabSmallMoleculeSectionRows) {
        self.small_molecule_data = rows;
    }
    /// The `NUC` section. Source `getNucleicAcidSectionRows`.
    pub fn nucleic_acid_section_rows(&self) -> &MzTabNucleicAcidSectionRows {
        &self.nucleic_acid_data
    }
    /// Replace the `NUC` section. Source `setNucleicAcidSectionRows`.
    pub fn set_nucleic_acid_section_rows(&mut self, rows: MzTabNucleicAcidSectionRows) {
        self.nucleic_acid_data = rows;
    }
    /// The `OLI` section. Source `getOligonucleotideSectionRows`.
    pub fn oligonucleotide_section_rows(&self) -> &MzTabOligonucleotideSectionRows {
        &self.oligonucleotide_data
    }
    /// Replace the `OLI` section. Source `setOligonucleotideSectionRows`.
    pub fn set_oligonucleotide_section_rows(&mut self, rows: MzTabOligonucleotideSectionRows) {
        self.oligonucleotide_data = rows;
    }
    /// The `OSM` section. Source `getOSMSectionRows`.
    pub fn osm_section_rows(&self) -> &MzTabOSMSectionRows {
        &self.osm_data
    }
    /// Replace the `OSM` section. Source `setOSMSectionRows`.
    pub fn set_osm_section_rows(&mut self, rows: MzTabOSMSectionRows) {
        self.osm_data = rows;
    }
    /// Line indices of empty rows. Source `getEmptyRows`.
    pub fn empty_rows(&self) -> &[usize] {
        &self.empty_rows
    }
    /// Replace the empty-row indices. Source `setEmptyRows`.
    pub fn set_empty_rows(&mut self, rows: Vec<usize>) {
        self.empty_rows = rows;
    }
    /// Comment lines by line index. Source `getCommentRows`.
    pub fn comment_rows(&self) -> &BTreeMap<usize, String> {
        &self.comment_rows
    }
    /// Replace the comment lines. Source `setCommentRows`.
    pub fn set_comment_rows(&mut self, rows: BTreeMap<usize, String>) {
        self.comment_rows = rows;
    }

    /// Number of PSMs in the `PSM` section, which is not the number of rows:
    /// a row is duplicated per parent protein, and rows that share a `PSM_ID`
    /// are one PSM.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when any row's `PSM_ID` cell is not a
    /// value, and [`Error::InvalidValue`] when the section exceeds
    /// [`MzTab::MAX_ROWS`]. The source relies on `PSM_ID` being set for every
    /// row, as its `@note` says, and throws `Exception::ElementNotFound` from
    /// `MzTabInteger::get()` when it is not.
    pub fn number_of_psms(&self) -> Result<usize> {
        if self.psm_data.len() > Self::MAX_ROWS {
            return Err(bad("MzTab PSM section exceeds its row limit"));
        }
        let mut ids = BTreeSet::new();
        for row in &self.psm_data {
            ids.insert(row.psm_id.get()?);
        }
        Ok(ids.len())
    }

    /// Optional column names of the `PRT` section, in first-occurrence order.
    ///
    /// # Errors
    ///
    /// As [`optional_column_names`].
    pub fn protein_optional_column_names(&self) -> Result<Vec<String>> {
        optional_column_names(&self.protein_data)
    }
    /// Optional column names of the `PEP` section, in first-occurrence order.
    ///
    /// # Errors
    ///
    /// As [`optional_column_names`].
    pub fn peptide_optional_column_names(&self) -> Result<Vec<String>> {
        optional_column_names(&self.peptide_data)
    }
    /// Optional column names of the `PSM` section, in first-occurrence order.
    ///
    /// # Errors
    ///
    /// As [`optional_column_names`].
    pub fn psm_optional_column_names(&self) -> Result<Vec<String>> {
        optional_column_names(&self.psm_data)
    }
    /// Optional column names of the `SML` section, in first-occurrence order.
    ///
    /// # Errors
    ///
    /// As [`optional_column_names`].
    pub fn small_molecule_optional_column_names(&self) -> Result<Vec<String>> {
        optional_column_names(&self.small_molecule_data)
    }
    /// Optional column names of the `NUC` section, in first-occurrence order.
    ///
    /// # Errors
    ///
    /// As [`optional_column_names`].
    pub fn nucleic_acid_optional_column_names(&self) -> Result<Vec<String>> {
        optional_column_names(&self.nucleic_acid_data)
    }
    /// Optional column names of the `OLI` section, in first-occurrence order.
    ///
    /// # Errors
    ///
    /// As [`optional_column_names`].
    pub fn oligonucleotide_optional_column_names(&self) -> Result<Vec<String>> {
        optional_column_names(&self.oligonucleotide_data)
    }
    /// Optional column names of the `OSM` section, in first-occurrence order.
    ///
    /// # Errors
    ///
    /// As [`optional_column_names`].
    pub fn osm_optional_column_names(&self) -> Result<Vec<String>> {
        optional_column_names(&self.osm_data)
    }
}

/// Append one `opt_<id>_<key>` column per key, taking the value from `meta`.
///
/// Spaces in a key become underscores, because a column name may not contain
/// whitespace; the *values* are taken verbatim. A key absent from `meta` still
/// produces a column, whose value is the default — and therefore null — cell,
/// so every row of a section declares the same columns.
///
/// `keys` is ordered, so the appended columns are in key order.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `keys` would push the entry count past
/// [`MzTab::MAX_OPTIONAL_COLUMNS`]. The entries are staged and appended only on
/// success, so `opt` is unchanged on error.
///
/// # Notes
///
/// Values render with the source's `DataValue::toString(true)` convention, so a
/// float list reads `[0.5, 1.4, -2.0, 0.1]` — full precision, trailing zeros
/// trimmed to at least one fractional digit. That is not Rust's `{}` for `f64`,
/// which would write `-2`.
pub fn add_meta_info_to_optional_columns(
    keys: &BTreeSet<String>,
    opt: &mut Vec<MzTabOptionalColumnEntry>,
    id: &str,
    meta: &MetaInfo,
) -> Result<()> {
    let total = opt
        .len()
        .checked_add(keys.len())
        .ok_or_else(|| bad("MzTab optional column count overflows"))?;
    if total > MzTab::MAX_OPTIONAL_COLUMNS {
        return Err(bad("MzTab section exceeds its optional-column limit"));
    }
    let mut staged = Vec::with_capacity(keys.len());
    for key in keys {
        let name = format!("opt_{id}_{}", key.replace(' ', "_"));
        let value = match meta.get(key) {
            Some(value) => MzTabString::from_text(&meta_value_text(value)),
            None => MzTabString::default(),
        };
        staged.push(MzTabOptionalColumnEntry::new(name, value));
    }
    opt.append(&mut staged);
    Ok(())
}

/// The source's `DataValue::toString(true)`: full-precision numbers, lists
/// bracketed and separated by `", "`, and an absent value as empty text.
fn meta_value_text(value: &MetaValue) -> String {
    fn bracketed<I: Iterator<Item = String>>(parts: I) -> String {
        let joined: Vec<String> = parts.collect();
        format!("[{}]", joined.join(", "))
    }
    match value.data() {
        MetaValueData::Empty => String::new(),
        MetaValueData::String(text) => text.clone(),
        MetaValueData::Integer(number) => number.to_string(),
        MetaValueData::Float(number) => format_float(*number, true),
        MetaValueData::StringList(values) => bracketed(values.iter().cloned()),
        MetaValueData::IntegerList(values) => bracketed(values.iter().map(i64::to_string)),
        MetaValueData::FloatList(values) => {
            bracketed(values.iter().map(|number| format_float(*number, true)))
        }
    }
}

/// What [`modification_metadata`] resolved, and what it could not.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ModificationMetaDataReport {
    /// The resolved entries, keyed by the one-based position of the input name.
    ///
    /// A name the registry does not know leaves a *gap*: the source increments
    /// its index for every input, resolved or not, so the keys are not
    /// necessarily `1..=n`.
    pub metadata: BTreeMap<usize, MzTabModificationMetaData>,
    /// Input names the registry did not resolve, in input order.
    ///
    /// The source only writes "Skipping unknown residue modification" to its
    /// warning log; returning the names makes the loss visible to the caller.
    pub skipped: Vec<String>,
}

/// Build `fixed_mod`/`variable_mod` metadata for `names`, using the global
/// modification registry.
///
/// # Errors
///
/// As [`modification_metadata_with`].
pub fn modification_metadata(names: &[String]) -> Result<ModificationMetaDataReport> {
    modification_metadata_with(ModificationsDB::global(), names)
}

/// Build `fixed_mod`/`variable_mod` metadata for `names` against `db`.
///
/// Each resolved name contributes one [`MzTabModificationMetaData`]: the term
/// carries the UniMod accession upper-cased with CV label `UNIMOD` when the
/// record has one and no accession otherwise, the name is the record's short
/// identifier, `position` is the terminal specificity spelled as the MzTab
/// specification requires, and `site` is the modified residue.
///
/// Passing a caller-owned registry replaces the source's global
/// `ModificationsDB::getInstance()` singleton.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `names` exceeds [`MAX_MODIFICATION_NAMES`].
/// Unknown names are reported in [`ModificationMetaDataReport::skipped`], not
/// as an error, because the source logs and continues.
pub fn modification_metadata_with(
    db: &ModificationsDB,
    names: &[String],
) -> Result<ModificationMetaDataReport> {
    if names.len() > MAX_MODIFICATION_NAMES {
        return Err(bad("MzTab modification name list exceeds its limit"));
    }
    let mut report = ModificationMetaDataReport::default();
    for (offset, name) in names.iter().enumerate() {
        match db.get_modification(name, None, None) {
            Ok(record) => {
                report
                    .metadata
                    .insert(offset + 1, modification_entry(record));
            }
            Err(_) => report.skipped.push(name.clone()),
        }
    }
    Ok(report)
}

fn modification_entry(record: &ResidueModification) -> MzTabModificationMetaData {
    let mut term = MzTabParameter::default();
    if let Some(accession) = record.unimod_accession() {
        if !accession.is_empty() {
            term.set_cv_label("UNIMOD");
            term.set_accession(accession.to_uppercase());
        }
    }
    term.set_name(record.name());
    let position = match record.term_specificity() {
        TermSpecificity::CTerm => "Any C-term",
        TermSpecificity::NTerm => "Any N-term",
        TermSpecificity::Anywhere => "Anywhere",
        TermSpecificity::ProteinCTerm => "Protein C-term",
        TermSpecificity::ProteinNTerm => "Protein N-term",
    };
    MzTabModificationMetaData {
        modification: term,
        // The source's ResidueModification stores 'X' for "any residue"; the
        // Rust record spells that absence as None.
        site: MzTabString::from_text(&record.origin().unwrap_or('X').to_string()),
        position: MzTabString::from_text(position),
    }
}

/// Build `variable_mod` metadata, or the specified placeholder when empty.
///
/// An empty `names` yields the single entry
/// `[MS, MS:1002454, No variable modifications searched, ]` at index 1, as the
/// source.
///
/// # Errors
///
/// As [`modification_metadata_with`], plus [`Error::Parse`] only if the
/// placeholder literal itself failed to parse, which cannot happen.
pub fn variable_modification_metadata(names: &[String]) -> Result<ModificationMetaDataReport> {
    if names.is_empty() {
        return placeholder_metadata("[MS, MS:1002454, No variable modifications searched, ]");
    }
    modification_metadata(names)
}

/// Build `fixed_mod` metadata, or the specified placeholder when empty.
///
/// An empty `names` yields the single entry
/// `[MS, MS:1002453, No fixed modifications searched, ]` at index 1, as the
/// source.
///
/// # Errors
///
/// As [`modification_metadata_with`], plus [`Error::Parse`] only if the
/// placeholder literal itself failed to parse, which cannot happen.
pub fn fixed_modification_metadata(names: &[String]) -> Result<ModificationMetaDataReport> {
    if names.is_empty() {
        return placeholder_metadata("[MS, MS:1002453, No fixed modifications searched, ]");
    }
    modification_metadata(names)
}

fn placeholder_metadata(cell: &str) -> Result<ModificationMetaDataReport> {
    let mut entry = MzTabModificationMetaData::default();
    entry.modification.from_cell_string(cell)?;
    let mut report = ModificationMetaDataReport::default();
    report.metadata.insert(1, entry);
    Ok(report)
}
