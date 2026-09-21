// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Fuzzy comparison for the test suite: the library's
//! `openms::concept::fuzzy_string_comparator` (a port of
//! `OpenMS::FuzzyStringComparator`, core bc9cc12) re-exported, plus the
//! `FuzzyDiff` tool contract (topp 174b576) as **test support**.
//!
//! Integration tests include it with
//!
//! ```text
//! #[path = "support/fuzzy_string_comparator.rs"]
//! mod fuzzy;
//! ```
//!
//! The comparator itself is library code since the `FuzzyDiff` tool was
//! ported; everything this module exported before is still exported under the
//! same names, with the same behaviour. What stays here is what a test needs
//! and the library does not: [`FuzzyDiffSettings`], a bounded reader for the
//! ParamXML subset of `FuzzyDiff.ini`, and [`fuzzy_diff`], the tool's exit-code
//! contract on two files. They depend on `std` and the always-built comparator
//! only, so a test compiled without the `paramxml` feature that the TOPP
//! framework and the real tool (`openms::cli::tools::FuzzyDiff`) need can still
//! reproduce a registered `${DIFF}` comparison. `tests/topp_fuzzy_diff.rs`
//! holds [`fuzzy_diff`] to the real tool's exit code on every executed oracle
//! case, so the two cannot drift apart unnoticed.
//!
//! The upstream test suite compares tool output with
//! `FuzzyDiff -test -ini FuzzyDiff.ini [-whitelist ...]`; a Rust test reproduces a
//! registered comparison with [`FuzzyDiffSettings::load_ini`] on the pinned
//! `tests/data/fuzzy_string_comparator/FuzzyDiff.ini`, the registration's
//! whitelist ([`FuzzyDiffSettings::with_whitelist`]) and [`fuzzy_diff`] or
//! [`FuzzyDiffSettings::compare_bytes`]. Tolerances are never widened here; a test
//! that needs other values sets them explicitly and names the reason.
//!
//! The mapping to the C++ members, the preserved quirks and the native
//! differences are listed in `docs/FUZZY_STRING_COMPARATOR_SUPPORT.md`.
#![allow(dead_code, unused_imports)]

pub use openms::concept::fuzzy_string_comparator::*;

use std::fs::File;
use std::io::Read;
use std::path::Path;

/// Largest ParamXML (INI) file [`FuzzyDiffSettings::from_ini`] accepts, in bytes.
pub const MAX_INI_BYTES: usize = 16 << 20;

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    needle.is_empty()
        || haystack
            .windows(needle.len())
            .any(|window| window == needle)
}

/// Exit codes `FuzzyDiff` returns (the used subset of `TOPPBase::ExitCodes`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FuzzyDiffExit {
    /// No difference found (`EXECUTION_OK`).
    ExecutionOk = 0,
    /// An input file does not exist (`INPUT_FILE_NOT_FOUND`).
    InputFileNotFound = 1,
    /// An input file cannot be read (`INPUT_FILE_NOT_READABLE`).
    InputFileNotReadable = 2,
    /// An input file is empty (`INPUT_FILE_EMPTY`).
    InputFileEmpty = 4,
    /// A parameter is out of its registered range or the INI is invalid
    /// (`ILLEGAL_PARAMETERS`).
    IllegalParameters = 6,
    /// `in1` or `in2` is empty (`MISSING_PARAMETERS`).
    MissingParameters = 7,
    /// A malformed `matched_whitelist` entry (`IllegalArgument`, `UNKNOWN_ERROR`).
    UnknownError = 8,
    /// Differences were found (`PARSE_ERROR`; the source notes it should find a
    /// better code).
    ParseError = 10,
}

impl FuzzyDiffExit {
    /// The numeric process exit code.
    pub fn code(self) -> i32 {
        self as i32
    }
}

/// Outcome of a [`fuzzy_diff`] run.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FuzzyDiffOutcome {
    /// The exit code the C++ tool returns for the same invocation.
    pub exit: FuzzyDiffExit,
    /// The comparator log, or a one-line reason for a pre-comparison exit. The
    /// TOPPBase framework messages of the tool are not reproduced.
    pub log: Vec<u8>,
}

impl FuzzyDiffOutcome {
    /// Whether the comparison passed (exit code 0).
    pub fn passed(&self) -> bool {
        self.exit == FuzzyDiffExit::ExecutionOk
    }

    /// The log as text, with invalid UTF-8 replaced.
    pub fn log_text(&self) -> String {
        String::from_utf8_lossy(&self.log).into_owned()
    }
}

/// The parameters of a `FuzzyDiff` invocation (topp 174b576 `FuzzyDiff.cpp`).
#[derive(Clone, Debug, PartialEq)]
pub struct FuzzyDiffSettings {
    /// Acceptable relative error (`-ratio`, at least 1).
    pub ratio: f64,
    /// Acceptable absolute difference (`-absdiff`, at least 0).
    pub absdiff: f64,
    /// Lines containing one of these strings on both sides are skipped
    /// (`-whitelist`).
    pub whitelist: Vec<String>,
    /// Colon-separated pairs `first:second` (`-matched_whitelist`).
    pub matched_whitelist: Vec<String>,
    /// Verbose level 0-3 (`-verbose`).
    pub verbose: i32,
    /// Tab width for column numbers, at least 1 (`-tab_width`).
    pub tab_width: i32,
    /// Number of the first column, at least 0 (`-first_column`).
    pub first_column: i32,
    /// Sort all lines but the first before comparing (`-sort`).
    pub sort: bool,
    /// Problems found while reading an INI file; any entry makes [`fuzzy_diff`]
    /// return [`FuzzyDiffExit::IllegalParameters`], as TOPPBase rejects an invalid
    /// INI during initialisation.
    pub ini_errors: Vec<String>,
}

impl Default for FuzzyDiffSettings {
    fn default() -> Self {
        Self::registered_defaults()
    }
}

impl FuzzyDiffSettings {
    /// The defaults `FuzzyDiff` registers: ratio 1, absdiff 0, whitelist
    /// `<?xml-stylesheet`, no matched whitelist, verbose 2, tab width 8, first
    /// column 1, no sorting.
    pub fn registered_defaults() -> Self {
        Self {
            ratio: 1.0,
            absdiff: 0.0,
            whitelist: vec!["<?xml-stylesheet".to_owned()],
            matched_whitelist: Vec::new(),
            verbose: 2,
            tab_width: 8,
            first_column: 1,
            sort: false,
            ini_errors: Vec::new(),
        }
    }

    /// Path of the pinned upstream `FuzzyDiff.ini` (test-data 0cb15f2
    /// `topp/FuzzyDiff.ini`) copied into this repository's test data.
    pub fn upstream_ini_path() -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/fuzzy_string_comparator/FuzzyDiff.ini")
    }

    /// Settings of the upstream `${DIFF}` command: the registered defaults
    /// overridden by the pinned `FuzzyDiff.ini` (ratio 1.01, absdiff 0.01,
    /// whitelist `<?xml-stylesheet`, verbose 1).
    ///
    /// # Errors
    ///
    /// Returns the reason when the pinned file cannot be read.
    pub fn upstream() -> Result<Self, String> {
        Self::load_ini(&Self::upstream_ini_path())
    }

    /// Read a ParamXML (INI) file; see [`Self::from_ini`].
    ///
    /// # Errors
    ///
    /// Returns the reason when the file cannot be read or is not well-formed
    /// ParamXML.
    pub fn load_ini(path: &Path) -> Result<Self, String> {
        let file =
            File::open(path).map_err(|e| format!("cannot open '{}': {e}", path.display()))?;
        let mut bytes = Vec::new();
        file.take(MAX_INI_BYTES as u64 + 1)
            .read_to_end(&mut bytes)
            .map_err(|e| format!("cannot read '{}': {e}", path.display()))?;
        Self::from_ini(&bytes)
    }

    /// Apply the values of a ParamXML (INI) file to the registered defaults, as
    /// `FuzzyDiff -ini <file>` does.
    ///
    /// Items of the tool instance section `FuzzyDiff:1:` are read: `ratio`,
    /// `absdiff`, `whitelist`, `matched_whitelist`, `verbose`, `tab_width`,
    /// `first_column` and `sort`. The common TOPP items (`in1`, `in2`, `log`,
    /// `debug`, `threads`, `no_progress`, `force`, `test`) and the tool `version`
    /// item are accepted and ignored. Any other item in the instance section, or a
    /// value that does not parse, is recorded in [`Self::ini_errors`]. Items outside
    /// the tool's sections are ignored.
    ///
    /// # Errors
    ///
    /// Returns the reason when the input exceeds [`MAX_INI_BYTES`] or is not
    /// well-formed ParamXML.
    pub fn from_ini(bytes: &[u8]) -> Result<Self, String> {
        let items = param_xml::parse(bytes)?;
        let mut settings = Self::registered_defaults();
        for item in items {
            let in_instance = item.path == ["FuzzyDiff", "1"];
            if !in_instance {
                continue;
            }
            let problem = settings.apply_ini_item(&item);
            if let Err(problem) = problem {
                settings.ini_errors.push(problem);
            }
        }
        Ok(settings)
    }

    fn apply_ini_item(&mut self, item: &param_xml::Item) -> Result<(), String> {
        use param_xml::Value;
        let key = format!("FuzzyDiff:1:{}", item.name);
        let scalar = || match &item.value {
            Value::Scalar(text) => Ok(text.as_str()),
            Value::List(_) => Err(format!("'{key}' must be a single value")),
        };
        let list = || match &item.value {
            Value::List(values) => Ok(values.clone()),
            Value::Scalar(_) => Err(format!("'{key}' must be a list")),
        };
        let double = |text: &str| {
            text.trim()
                .parse::<f64>()
                .map_err(|_| format!("'{key}' is not a number: '{text}'"))
        };
        let int = |text: &str| {
            text.trim()
                .parse::<i32>()
                .map_err(|_| format!("'{key}' is not an integer: '{text}'"))
        };
        match item.name.as_str() {
            "ratio" => self.ratio = double(scalar()?)?,
            "absdiff" => self.absdiff = double(scalar()?)?,
            "whitelist" => self.whitelist = list()?,
            "matched_whitelist" => self.matched_whitelist = list()?,
            "verbose" => self.verbose = int(scalar()?)?,
            "tab_width" => self.tab_width = int(scalar()?)?,
            "first_column" => self.first_column = int(scalar()?)?,
            "sort" => {
                self.sort = match scalar()? {
                    "true" => true,
                    "false" => false,
                    other => return Err(format!("'{key}' is not a flag: '{other}'")),
                }
            }
            "in1" | "in2" | "log" | "debug" | "threads" | "no_progress" | "force" | "test" => {}
            other => {
                return Err(format!(
                    "Unknown (or deprecated) Parameter 'FuzzyDiff:1:{other}'"
                ));
            }
        }
        Ok(())
    }

    /// Replace the whitelist, as a registration's `-whitelist a b` does (the
    /// command-line list replaces the INI list; it is not appended).
    pub fn with_whitelist(mut self, entries: &[&str]) -> Self {
        self.whitelist = entries.iter().map(|&e| e.to_owned()).collect();
        self
    }

    /// Replace the matched whitelist, as `-matched_whitelist a:b` does.
    pub fn with_matched_whitelist(mut self, entries: &[&str]) -> Self {
        self.matched_whitelist = entries.iter().map(|&e| e.to_owned()).collect();
        self
    }

    /// A comparator configured with these settings, logging into a buffer.
    ///
    /// # Errors
    ///
    /// Returns the exit code and reason `FuzzyDiff` would fail with before
    /// comparing: an INI error, a parameter out of range, or a malformed matched
    /// whitelist entry.
    pub fn comparator(&self) -> Result<FuzzyStringComparator, (FuzzyDiffExit, String)> {
        if let Some(problem) = self.ini_errors.first() {
            return Err((FuzzyDiffExit::IllegalParameters, problem.clone()));
        }
        self.check_ranges()?;
        let matched = self.parsed_matched_whitelist()?;
        let mut comparator = FuzzyStringComparator::new();
        comparator.set_log_destination(LogDestination::Buffer);
        comparator.set_acceptable_relative(self.ratio);
        comparator.set_acceptable_absolute(self.absdiff);
        comparator.set_whitelist(self.whitelist.clone());
        comparator.set_matched_whitelist(matched);
        comparator.set_verbose_level(self.verbose);
        comparator.set_tab_width(self.tab_width);
        comparator.set_first_column(self.first_column);
        Ok(comparator)
    }

    /// Compare two in-memory texts with these settings: the `FuzzyDiff` comparison
    /// without its file checks and same-name check. Sorting applies when enabled.
    ///
    /// # Errors
    ///
    /// Returns the comparator log (or the pre-comparison reason) when the texts
    /// differ or the settings are invalid.
    pub fn compare_bytes(&self, actual: &[u8], expected: &[u8]) -> Result<(), String> {
        let mut comparator = self.comparator().map_err(|(_, reason)| reason)?;
        let passed = if self.sort {
            comparator.compare_bytes(&sorted_lines(actual), &sorted_lines(expected))
        } else {
            comparator.compare_bytes(actual, expected)
        };
        if passed {
            Ok(())
        } else {
            Err(String::from_utf8_lossy(comparator.log()).into_owned())
        }
    }

    #[allow(clippy::manual_range_contains)]
    fn check_ranges(&self) -> Result<(), (FuzzyDiffExit, String)> {
        // getDoubleOption_/getIntOption_ check a value only when it differs from
        // the registered default, with plain comparisons: NaN passes, as in C++.
        let defaults = Self::registered_defaults();
        let invalid = |name: &str, value: String| {
            Err((
                FuzzyDiffExit::IllegalParameters,
                format!("Invalid value '{value}' for parameter '{name}' given."),
            ))
        };
        if self.ratio != defaults.ratio && self.ratio < 1.0 {
            return invalid("ratio", self.ratio.to_string());
        }
        if self.absdiff != defaults.absdiff && self.absdiff < 0.0 {
            return invalid("absdiff", self.absdiff.to_string());
        }
        if self.verbose != defaults.verbose && (self.verbose < 0 || self.verbose > 3) {
            return invalid("verbose", self.verbose.to_string());
        }
        if self.tab_width != defaults.tab_width && self.tab_width < 1 {
            return invalid("tab_width", self.tab_width.to_string());
        }
        if self.first_column != defaults.first_column && self.first_column < 0 {
            return invalid("first_column", self.first_column.to_string());
        }
        Ok(())
    }

    fn parsed_matched_whitelist(&self) -> Result<Vec<(String, String)>, (FuzzyDiffExit, String)> {
        // The library's split, which the real tool uses too; its only error is
        // the source's IllegalArgument message.
        parse_matched_whitelist(&self.matched_whitelist).map_err(|error| match error {
            openms::Error::InvalidValue(message) => (FuzzyDiffExit::UnknownError, message),
            other => (FuzzyDiffExit::UnknownError, other.to_string()),
        })
    }
}

/// Run the `FuzzyDiff` tool contract on two files and return the exit code the
/// C++ tool gives for the same invocation.
///
/// Order of checks follows TOPPBase and `FuzzyDiff::main_`, as executed on the
/// oracle: during initialisation, an INI error or a parameter outside its
/// registered range, from the INI or the command line (6); then for `in1` and then
/// `in2`, an empty name (7), a missing file (1), an unreadable file (2) or an empty
/// regular file (4); then a matched-whitelist entry that does not split into
/// exactly two parts at `:` (8); then the comparison, 0 when equal and 10
/// otherwise. With `sort` the lines after the first are sorted bytewise (as
/// `std::getline` splits them, at `\n` only) and compared in memory, so the
/// same-name check does not apply, as with the source's temporary files.
pub fn fuzzy_diff(in1: &Path, in2: &Path, settings: &FuzzyDiffSettings) -> FuzzyDiffOutcome {
    let early = |exit: FuzzyDiffExit, reason: String| FuzzyDiffOutcome {
        exit,
        log: format!("{reason}\n").into_bytes(),
    };
    if let Some(problem) = settings.ini_errors.first() {
        return early(FuzzyDiffExit::IllegalParameters, problem.clone());
    }
    if let Err((exit, reason)) = settings.check_ranges() {
        return early(exit, reason);
    }
    for (name, path) in [("in1", in1), ("in2", in2)] {
        if let Err((exit, reason)) = check_input_file(name, path) {
            return early(exit, reason);
        }
    }
    let mut comparator = match settings.comparator() {
        Ok(comparator) => comparator,
        Err((exit, reason)) => return early(exit, reason),
    };
    let passed = if settings.sort {
        match (read_bounded(in1), read_bounded(in2)) {
            (Ok(a), Ok(b)) => {
                comparator.set_input_names(&in1.to_string_lossy(), &in2.to_string_lossy());
                comparator.compare_bytes(&sorted_lines(&a), &sorted_lines(&b))
            }
            (Err(reason), _) | (_, Err(reason)) => {
                return early(FuzzyDiffExit::InputFileNotReadable, reason);
            }
        }
    } else {
        comparator.compare_files(in1, in2)
    };
    FuzzyDiffOutcome {
        exit: if passed {
            FuzzyDiffExit::ExecutionOk
        } else {
            FuzzyDiffExit::ParseError
        },
        log: comparator.take_log(),
    }
}

fn check_input_file(name: &str, path: &Path) -> Result<(), (FuzzyDiffExit, String)> {
    if path.as_os_str().is_empty() {
        return Err((
            FuzzyDiffExit::MissingParameters,
            format!("Error: The required parameter '{name}' was not given or is empty!"),
        ));
    }
    let Ok(metadata) = std::fs::metadata(path) else {
        return Err((
            FuzzyDiffExit::InputFileNotFound,
            format!("Error: File not found ({})", path.display()),
        ));
    };
    if !metadata.is_dir() && File::open(path).is_err() {
        return Err((
            FuzzyDiffExit::InputFileNotReadable,
            format!("Error: File not readable ({})", path.display()),
        ));
    }
    if !metadata.is_dir() && metadata.len() == 0 {
        return Err((
            FuzzyDiffExit::InputFileEmpty,
            format!("Error: File empty ({})", path.display()),
        ));
    }
    Ok(())
}

fn read_bounded(path: &Path) -> Result<Vec<u8>, String> {
    let file = File::open(path).map_err(|e| format!("cannot open '{}': {e}", path.display()))?;
    let mut bytes = Vec::new();
    file.take(MAX_INPUT_BYTES + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| format!("cannot read '{}': {e}", path.display()))?;
    if bytes.len() as u64 > MAX_INPUT_BYTES {
        return Err(format!(
            "input '{}' exceeds the comparison limit of {MAX_INPUT_BYTES} bytes",
            path.display()
        ));
    }
    Ok(bytes)
}

/// A bounded reader for the subset of ParamXML that INI files use: `NODE`,
/// `ITEM`, `ITEMLIST` and `LISTITEM` elements with quoted attributes.
mod param_xml {
    /// Largest element nesting accepted.
    const MAX_DEPTH: usize = 64;
    /// Largest number of items accepted.
    const MAX_ITEMS: usize = 100_000;

    /// A parameter value: a scalar string or a list of strings.
    #[derive(Clone, Debug, PartialEq)]
    pub enum Value {
        /// The `value` attribute of an `ITEM`.
        Scalar(String),
        /// The `value` attributes of the `LISTITEM`s of an `ITEMLIST`.
        List(Vec<String>),
    }

    /// One `ITEM` or `ITEMLIST` with the names of its enclosing `NODE`s.
    #[derive(Clone, Debug, PartialEq)]
    pub struct Item {
        /// Names of the enclosing nodes, outermost first.
        pub path: Vec<String>,
        /// The item name.
        pub name: String,
        /// The item value.
        pub value: Value,
    }

    struct Tag {
        name: String,
        attributes: Vec<(String, String)>,
        closing: bool,
        self_closing: bool,
    }

    /// Parse the items of a ParamXML document.
    pub fn parse(bytes: &[u8]) -> Result<Vec<Item>, String> {
        if bytes.len() > super::MAX_INI_BYTES {
            return Err(format!("INI exceeds {} bytes", super::MAX_INI_BYTES));
        }
        let latin1 = declared_latin1(bytes);
        let mut items = Vec::new();
        let mut path: Vec<String> = Vec::new();
        let mut list: Option<(String, Vec<String>)> = None;
        let mut position = 0;
        while let Some(offset) = bytes
            .get(position..)
            .and_then(|rest| rest.iter().position(|&b| b == b'<'))
        {
            let start = position + offset;
            let rest = &bytes[start..];
            if rest.starts_with(b"<?") {
                position =
                    start + find(rest, b"?>").ok_or("unterminated processing instruction")? + 2;
                continue;
            }
            if rest.starts_with(b"<!--") {
                position = start + find(rest, b"-->").ok_or("unterminated comment")? + 3;
                continue;
            }
            if rest.starts_with(b"<!") {
                position = start + find(rest, b">").ok_or("unterminated declaration")? + 1;
                continue;
            }
            let (tag, length) = parse_tag(rest, latin1)?;
            position = start + length;
            let attribute = |key: &str| {
                tag.attributes
                    .iter()
                    .find(|(k, _)| k == key)
                    .map(|(_, v)| v.clone())
            };
            match (tag.name.as_str(), tag.closing) {
                ("NODE", false) => {
                    if !tag.self_closing {
                        if path.len() >= MAX_DEPTH {
                            return Err("INI nesting is too deep".into());
                        }
                        path.push(attribute("name").ok_or("NODE without name")?);
                    }
                }
                ("NODE", true) => {
                    path.pop().ok_or("unbalanced </NODE>")?;
                }
                ("ITEM", false) => {
                    let name = attribute("name").ok_or("ITEM without name")?;
                    let value = attribute("value").ok_or("ITEM without value")?;
                    items.push(Item {
                        path: path.clone(),
                        name,
                        value: Value::Scalar(value),
                    });
                }
                ("ITEMLIST", false) => {
                    let name = attribute("name").ok_or("ITEMLIST without name")?;
                    if tag.self_closing {
                        items.push(Item {
                            path: path.clone(),
                            name,
                            value: Value::List(Vec::new()),
                        });
                    } else {
                        list = Some((name, Vec::new()));
                    }
                }
                ("LISTITEM", false) => {
                    let value = attribute("value").ok_or("LISTITEM without value")?;
                    let (_, values) = list.as_mut().ok_or("LISTITEM outside ITEMLIST")?;
                    values.push(value);
                }
                ("ITEMLIST", true) => {
                    let (name, values) = list.take().ok_or("unbalanced </ITEMLIST>")?;
                    items.push(Item {
                        path: path.clone(),
                        name,
                        value: Value::List(values),
                    });
                }
                _ => {}
            }
            if items.len() > MAX_ITEMS {
                return Err("INI holds too many items".into());
            }
        }
        if !path.is_empty() || list.is_some() {
            return Err("INI ends inside an element".into());
        }
        Ok(items)
    }

    fn declared_latin1(bytes: &[u8]) -> bool {
        let head = &bytes[..bytes.len().min(200)];
        let lower: Vec<u8> = head.iter().map(u8::to_ascii_lowercase).collect();
        super::contains(&lower, b"encoding=\"iso-8859-1\"")
            || super::contains(&lower, b"encoding='iso-8859-1'")
    }

    fn find(haystack: &[u8], needle: &[u8]) -> Option<usize> {
        haystack.windows(needle.len()).position(|w| w == needle)
    }

    fn parse_tag(bytes: &[u8], latin1: bool) -> Result<(Tag, usize), String> {
        let mut i = 1;
        let closing = bytes.get(i) == Some(&b'/');
        if closing {
            i += 1;
        }
        let name_start = i;
        while bytes
            .get(i)
            .is_some_and(|b| !super::is_c_space(*b) && *b != b'>' && *b != b'/')
        {
            i += 1;
        }
        let name = String::from_utf8_lossy(&bytes[name_start..i]).into_owned();
        let mut attributes = Vec::new();
        loop {
            while bytes.get(i).is_some_and(|b| super::is_c_space(*b)) {
                i += 1;
            }
            match bytes.get(i) {
                None => return Err("unterminated tag".into()),
                Some(b'>') => {
                    return Ok((
                        Tag {
                            name,
                            attributes,
                            closing,
                            self_closing: false,
                        },
                        i + 1,
                    ));
                }
                Some(b'/') if bytes.get(i + 1) == Some(&b'>') => {
                    return Ok((
                        Tag {
                            name,
                            attributes,
                            closing,
                            self_closing: true,
                        },
                        i + 2,
                    ));
                }
                Some(_) => {}
            }
            let key_start = i;
            while bytes
                .get(i)
                .is_some_and(|b| *b != b'=' && !super::is_c_space(*b) && *b != b'>')
            {
                i += 1;
            }
            let key = String::from_utf8_lossy(&bytes[key_start..i]).into_owned();
            while bytes.get(i).is_some_and(|b| super::is_c_space(*b)) {
                i += 1;
            }
            if bytes.get(i) != Some(&b'=') {
                return Err(format!("attribute '{key}' without value"));
            }
            i += 1;
            while bytes.get(i).is_some_and(|b| super::is_c_space(*b)) {
                i += 1;
            }
            let quote = *bytes.get(i).ok_or("unterminated attribute")?;
            if quote != b'"' && quote != b'\'' {
                return Err(format!("attribute '{key}' is not quoted"));
            }
            let value_start = i + 1;
            let length = bytes
                .get(value_start..)
                .and_then(|rest| rest.iter().position(|&b| b == quote))
                .ok_or("unterminated attribute value")?;
            let raw = &bytes[value_start..value_start + length];
            attributes.push((key, decode(raw, latin1)?));
            i = value_start + length + 1;
        }
    }

    fn decode(raw: &[u8], latin1: bool) -> Result<String, String> {
        let text: String = if latin1 {
            raw.iter().map(|&b| char::from(b)).collect()
        } else {
            String::from_utf8(raw.to_vec()).map_err(|_| "attribute is not UTF-8".to_owned())?
        };
        let mut out = String::with_capacity(text.len());
        let mut rest = text.as_str();
        while let Some(amp) = rest.find('&') {
            out.push_str(&rest[..amp]);
            let after = &rest[amp + 1..];
            let semicolon = after.find(';').ok_or("unterminated entity")?;
            let entity = &after[..semicolon];
            let decoded = match entity {
                "lt" => '<',
                "gt" => '>',
                "amp" => '&',
                "quot" => '"',
                "apos" => '\'',
                _ => {
                    let code = if let Some(hex) = entity.strip_prefix("#x") {
                        u32::from_str_radix(hex, 16).ok()
                    } else if let Some(dec) = entity.strip_prefix('#') {
                        dec.parse::<u32>().ok()
                    } else {
                        None
                    };
                    code.and_then(char::from_u32)
                        .ok_or_else(|| format!("unknown entity '&{entity};'"))?
                }
            };
            out.push(decoded);
            rest = &after[semicolon + 1..];
        }
        out.push_str(rest);
        Ok(out)
    }
}
