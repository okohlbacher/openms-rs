// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Resolved tool parameters and run services: the native form of the source
//! `TOPPBase` `get*_` accessors, `parseRange_`, `inputFileReadable_`,
//! `outputFileWritable_`, `setMaxNumberOfThreads`, the progress log type and
//! the test-mode unique-id seed.
//!
//! See `docs/TOPP_CLI_SUPPORT.md` for the supported source subset.

use super::parameter::ExitCode;
use super::processing::{self, AddDataProcessing};
use crate::concept::UniqueIdGenerator;
use crate::concept::parallel::Threads;
use crate::concept::progress_logger::ProgressLogType;
use crate::metadata::{DataProcessing, ProcessingAction};
use crate::param::{Param, ParamValue};
use crate::system::file;
use crate::{Error, Result};
use std::io::Write;

/// Seed of the unique-id generator under `-test`, as `TOPPBase::main` sets it
/// (`TOPPBase.cpp:369-376`). Its first two raw draws are
/// 5233264595117471314 and 4835329514588776807.
pub const TEST_MODE_UNIQUE_ID_SEED: u64 = 19_991_231_235_959;

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

/// Parameters and services of one tool run, after defaults, INI file and
/// command line are merged and validated.
///
/// The source exposes these through protected `getStringOption_` style methods
/// and members on the tool itself. Here they are a borrowed context passed to
/// the tool body, so a tool cannot mutate its own resolved parameters mid-run,
/// and every service a later tool needs — progress logging, the thread policy,
/// the unique-id generator and processing annotations — derives from the same
/// resolved values.
#[derive(Clone, Debug)]
pub struct ToolContext {
    tool_name: String,
    version: String,
    ini_location: String,
    param: Param,
    debug_level: i64,
    test_mode: bool,
    no_progress: bool,
    force: bool,
    threads: i64,
}

impl ToolContext {
    /// A context over `param`, the resolved tree without the `<tool>:1:` prefix.
    pub(crate) fn new(tool_name: &str, version: &str, param: Param) -> Self {
        let int = |key: &str, fallback: i64| match param.value(key) {
            Ok(ParamValue::Integer(value)) => *value,
            _ => fallback,
        };
        let flag =
            |key: &str| matches!(param.value(key), Ok(ParamValue::String(text)) if text == "true");
        let debug_level = int("debug", 0);
        let threads = int("threads", 1);
        let test_mode = flag("test");
        let no_progress = flag("no_progress");
        let force = flag("force");
        Self {
            tool_name: tool_name.to_owned(),
            version: version.to_owned(),
            ini_location: format!("{tool_name}:1:"),
            param,
            debug_level,
            test_mode,
            no_progress,
            force,
            threads,
        }
    }

    /// The executable and INI section name, as `toolName_`.
    pub fn tool_name(&self) -> &str {
        &self.tool_name
    }
    /// The product version the tool reports, as the source `version_`.
    pub fn version(&self) -> &str {
        &self.version
    }
    /// The INI section this run reads, as `getIniLocation_`: always
    /// `<tool>:1:`, because `-instance` is rejected like in the source.
    pub fn ini_location(&self) -> &str {
        &self.ini_location
    }
    /// The complete resolved parameter tree, as `getParam_`.
    ///
    /// Keys carry no `<tool>:1:` prefix, so a subsection parameter reads as
    /// `algorithm:peakcount`, exactly as in the source.
    pub fn param(&self) -> &Param {
        &self.param
    }
    /// Source `-debug`, the debug level; zero disables debug output.
    pub fn debug_level(&self) -> i64 {
        self.debug_level
    }
    /// Source `-test`: outputs must not depend on the clock, the machine or
    /// absolute paths.
    pub fn test_mode(&self) -> bool {
        self.test_mode
    }
    /// Source `-no_progress`: progress logging is disabled.
    pub fn no_progress(&self) -> bool {
        self.no_progress
    }
    /// Source `-force`: the tool may override its own safety checks.
    pub fn force(&self) -> bool {
        self.force
    }
    /// Source `-threads`; 0 means every available core.
    pub fn threads(&self) -> i64 {
        self.threads
    }

    /// The progress logger type the tool hands to loaders and algorithms.
    ///
    /// Source `log_type_` (`TOPPBase.cpp:400-403`): `CMD` unless `-no_progress`
    /// is given, when it stays `NONE`.
    pub fn progress_log_type(&self) -> ProgressLogType {
        if self.no_progress {
            ProgressLogType::None
        } else {
            ProgressLogType::Cmd
        }
    }
    /// The worker-thread policy for parallel algorithms.
    ///
    /// Source `setMaxNumberOfThreads(getParamAsInt_("threads", 1))`
    /// (`TOPPBase.cpp:84-98, 408`), where zero or a negative count means every
    /// available processor. The source sets a process-wide OpenMP limit; here
    /// the policy is passed to each computation, which keeps the parallel
    /// result bit-identical to the serial one (`src/concept/parallel.rs`).
    /// A negative count yields one worker, as [`Threads::from_cli`] documents.
    pub fn thread_policy(&self) -> Threads {
        Threads::from_cli(self.threads)
    }
    /// A unique-id generator for this run's outputs.
    ///
    /// Under `-test` the generator starts from [`TEST_MODE_UNIQUE_ID_SEED`], as
    /// `UniqueIdGenerator::setSeed(19991231235959)` in `TOPPBase::main`, so the
    /// same draws are made on every machine. Otherwise it is seeded from the
    /// clock and process id. The source seeds one process-wide generator; each
    /// call here returns an independent generator starting at the first draw,
    /// so a tool should create one and reuse it.
    pub fn unique_id_generator(&self) -> UniqueIdGenerator {
        if self.test_mode {
            UniqueIdGenerator::from_seed(TEST_MODE_UNIQUE_ID_SEED)
        } else {
            UniqueIdGenerator::new()
        }
    }
    /// The processing record of this run, as `getProcessingInfo_`.
    ///
    /// See [`TEST_MODE_VERSION`](super::TEST_MODE_VERSION) for the values
    /// recorded under `-test`; otherwise the product version, the current
    /// local time and every resolved parameter as `parameter: <name>` are
    /// recorded.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a floating-point parameter cannot
    /// be represented as metadata.
    pub fn processing_info(&self, actions: &[ProcessingAction]) -> Result<DataProcessing> {
        processing::processing_info(
            &self.tool_name,
            &self.version,
            self.test_mode,
            &self.param,
            actions,
        )
    }
    /// Attach `processing` to an output map, as `addDataProcessing_`.
    ///
    /// A feature map and a consensus map gain one entry; every spectrum and
    /// chromatogram of an experiment gains one shared entry. Under `-test` a
    /// consensus map's column-header file names are reduced to base names.
    pub fn add_data_processing<T: AddDataProcessing + ?Sized>(
        &self,
        target: &mut T,
        processing: &DataProcessing,
    ) {
        target.add_data_processing(processing, self.test_mode);
    }

    fn value(&self, name: &str) -> Result<&ParamValue> {
        self.param
            .value(name)
            .map_err(|_| bad(format!("parameter '{name}' was not registered")))
    }

    /// A string, input-file, output-file or output-prefix option, as
    /// `getStringOption_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `name` is not registered or does
    /// not hold a string; the source throws `UnregisteredParameter` or
    /// `WrongParameterType`.
    pub fn string(&self, name: &str) -> Result<&str> {
        match self.value(name)? {
            ParamValue::String(text) => Ok(text),
            _ => Err(bad(format!("parameter '{name}' is not a string"))),
        }
    }
    /// An integer option, as `getIntOption_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `name` is not registered or does
    /// not hold an integer.
    pub fn int(&self, name: &str) -> Result<i64> {
        match self.value(name)? {
            ParamValue::Integer(value) => Ok(*value),
            _ => Err(bad(format!("parameter '{name}' is not an integer"))),
        }
    }
    /// A floating-point option, as `getDoubleOption_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `name` is not registered or does
    /// not hold a floating-point value. An integer is not widened: the source
    /// throws `WrongParameterType` for an integer option, and the strict update
    /// already refuses an INI value whose type differs from the registered one.
    pub fn double(&self, name: &str) -> Result<f64> {
        match self.value(name)? {
            ParamValue::Float(value) => Ok(*value),
            _ => Err(bad(format!(
                "parameter '{name}' is not a floating-point number"
            ))),
        }
    }
    /// A string, input-file or output-file list, as `getStringList_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `name` is not registered or does
    /// not hold a string list.
    pub fn string_list(&self, name: &str) -> Result<&[String]> {
        match self.value(name)? {
            ParamValue::StringList(values) => Ok(values),
            _ => Err(bad(format!("parameter '{name}' is not a string list"))),
        }
    }
    /// An integer list, as `getIntList_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `name` is not registered or does
    /// not hold an integer list.
    pub fn int_list(&self, name: &str) -> Result<&[i32]> {
        match self.value(name)? {
            ParamValue::IntegerList(values) => Ok(values),
            _ => Err(bad(format!("parameter '{name}' is not an integer list"))),
        }
    }
    /// A floating-point list, as `getDoubleList_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `name` is not registered or does
    /// not hold a floating-point list.
    pub fn double_list(&self, name: &str) -> Result<&[f64]> {
        match self.value(name)? {
            ParamValue::FloatList(values) => Ok(values),
            _ => Err(bad(format!("parameter '{name}' is not a float list"))),
        }
    }
    /// Source `getFlag_`: a registered flag is true only when it was given.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `name` is not registered, does not
    /// hold a string, or holds a string other than `true` or `false`; the source
    /// throws `WrongParameterType` and `InvalidParameter` respectively.
    pub fn flag(&self, name: &str) -> Result<bool> {
        match self.value(name)? {
            ParamValue::String(text) if text == "true" => Ok(true),
            ParamValue::String(text) if text == "false" => Ok(false),
            ParamValue::String(text) => Err(bad(format!(
                "Invalid value '{text}' for flag parameter '{name}'. Valid values are 'true' and 'false' only."
            ))),
            _ => Err(bad(format!("parameter '{name}' is not a flag"))),
        }
    }
    /// Values of a registered subsection with the subsection prefix removed,
    /// as `getParam_().copy("<name>:", true)`.
    ///
    /// # Errors
    ///
    /// Propagates parameter-tree failures; an unknown subsection yields an
    /// empty tree, as the source copy does.
    pub fn subsection(&self, name: &str) -> Result<Param> {
        self.param.copy(&format!("{name}:"), true)
    }
}

/// Source `parseRange_` for floating-point bounds (`TOPPBase.cpp:2016-2052`).
///
/// `":8"`, `"2:"`, `"2:8"` and `":"` are accepted, and an absent side leaves
/// that bound untouched. The part before the first colon sets `low` and the
/// part after the last colon sets `high`, each converted like
/// `StringUtils::toDouble`. Returns whether any bound was set. As in the
/// source, `low > high` is not an error: the caller receives an empty range.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the colon is missing, where the source
/// throws `ConversionError` rather than reading `400` as `400:400`, or when a
/// bound is not a number. Neither bound is modified on error; the source may
/// already have assigned `low` when `high` fails to convert.
pub fn parse_range(text: &str, low: &mut f64, high: &mut f64) -> Result<bool> {
    if !text.contains(':') {
        return Err(bad(format!(
            "Invalid range '{text}': expected format '[min]:[max]' (the ':' separator is missing)"
        )));
    }
    let conversion = || {
        bad(format!(
            "Could not convert string '{text}' to a range of floating point values"
        ))
    };
    let start = text.split(':').next().unwrap_or("");
    let end = text.rsplit(':').next().unwrap_or("");
    let new_low = if start.is_empty() {
        None
    } else {
        Some(to_double(start).map_err(|_| conversion())?)
    };
    let new_high = if end.is_empty() {
        None
    } else {
        Some(to_double(end).map_err(|_| conversion())?)
    };
    if let Some(value) = new_low {
        *low = value;
    }
    if let Some(value) = new_high {
        *high = value;
    }
    Ok(new_low.is_some() || new_high.is_some())
}

/// Check that `filename` can be read, as `inputFileReadable_`
/// (`TOPPBase.cpp:1968-1995`).
///
/// Returns `None` when the file exists, is readable and — unless it is a
/// directory — holds at least one byte. Otherwise writes the source's two
/// diagnostic lines to `err` and returns the exit code the source exception
/// maps to: [`ExitCode::InputFileNotFound`] (`FileNotFound`),
/// [`ExitCode::InputFileNotReadable`] (`FileNotReadable`) or
/// [`ExitCode::InputFileEmpty`] (`FileEmpty`). `param_name` names the option
/// in the first line; an empty name gives the source's generic wording.
/// Diagnostics are best effort: a failing `err` does not change the result.
pub fn input_file_readable(
    filename: &str,
    param_name: &str,
    err: &mut dyn Write,
) -> Option<ExitCode> {
    let (code, detail) = if !file::exists(filename) {
        (
            ExitCode::InputFileNotFound,
            format!("Error: File not found (the file '{filename}' does not exist)"),
        )
    } else if !file::readable(filename) {
        (
            ExitCode::InputFileNotReadable,
            format!(
                "Error: File not readable (the file '{filename}' is not readable for the current user)"
            ),
        )
    } else if !file::is_directory(filename) && file::empty(filename) {
        (
            ExitCode::InputFileEmpty,
            format!("Error: File empty (the file '{filename}' is empty)"),
        )
    } else {
        return None;
    };
    let heading = if param_name.is_empty() {
        "Cannot read input file!".to_owned()
    } else {
        format!("Cannot read input file given from parameter '-{param_name}'!")
    };
    let _ = writeln!(err, "{heading}");
    let _ = writeln!(err, "{detail}");
    Some(code)
}

/// Check that `filename` can be written, as `outputFileWritable_`
/// (`TOPPBase.cpp:1997-2013`).
///
/// Returns `None` when [`file::writable`] answers yes; that query never creates
/// or removes a file under the caller's name. Otherwise writes the source's two
/// diagnostic lines to `err` and returns [`ExitCode::CannotWriteOutputFile`],
/// the code of the source's `UnableToCreateFile`. Diagnostics are best effort.
pub fn output_file_writable(
    filename: &str,
    param_name: &str,
    err: &mut dyn Write,
) -> Option<ExitCode> {
    if file::writable(filename) {
        return None;
    }
    let heading = if param_name.is_empty() {
        "Cannot write output file!".to_owned()
    } else {
        format!("Cannot write output file given from parameter '-{param_name}'!")
    };
    let _ = writeln!(err, "{heading}");
    let _ = writeln!(
        err,
        "Error: Unable to write file (the file '{filename}' could not be created. )"
    );
    Some(ExitCode::CannotWriteOutputFile)
}

/// Whether the extension of `path` is one of `formats`, case-insensitively.
/// An empty format list accepts anything, as in source.
pub(crate) fn extension_allowed(path: &str, formats: &[String]) -> bool {
    if formats.is_empty() {
        return true;
    }
    let lower = path.to_ascii_lowercase();
    formats
        .iter()
        .any(|f| lower.ends_with(&format!(".{}", f.to_ascii_lowercase())))
}

/// The source's four whitespace characters, skipped from `index` on.
fn skip_whitespace(bytes: &[u8], mut index: usize) -> usize {
    while bytes
        .get(index)
        .is_some_and(|b| matches!(*b, b' ' | b'\t' | b'\n' | b'\r'))
    {
        index += 1;
    }
    index
}

/// Source `StringUtils::toInt32` (`StringUtils.cpp:136-166`), returning the
/// source's `ConversionError` message as the error.
///
/// Leading and trailing space, tab, newline and carriage return are skipped,
/// one `+` may precede the number, and the rest must be a complete decimal
/// `i32`. Like `std::from_chars` after the source strips `+`, a `-` may still
/// follow it.
pub(crate) fn to_int32(text: &str) -> std::result::Result<i32, String> {
    let bytes = text.as_bytes();
    let not_converted = || format!("Could not convert string '{text}' to an integer value");
    let mut cursor = skip_whitespace(bytes, 0);
    if cursor == bytes.len() {
        return Err(not_converted());
    }
    if bytes[cursor] == b'+' {
        cursor += 1;
    }
    let number_start = cursor;
    if bytes.get(cursor) == Some(&b'-') {
        cursor += 1;
    }
    let digits_start = cursor;
    while bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
        cursor += 1;
    }
    if cursor == digits_start {
        return Err(not_converted());
    }
    let value = text
        .get(number_start..cursor)
        .and_then(|digits| digits.parse::<i32>().ok())
        .ok_or_else(not_converted)?;
    let after = skip_whitespace(bytes, cursor);
    if after != bytes.len() {
        return Err(format!(
            "Prefix of string '{text}' successfully converted to an int32 value. Additional characters found at position {}",
            after + 1
        ));
    }
    Ok(value)
}

/// Length of a case-insensitive ASCII `word` at `bytes[index..]`, if present.
fn word_at(bytes: &[u8], index: usize, word: &str) -> Option<usize> {
    let candidate = bytes.get(index..index.checked_add(word.len())?)?;
    candidate
        .eq_ignore_ascii_case(word.as_bytes())
        .then_some(word.len())
}

/// Source `StringUtils::toDouble` (`StringUtils.cpp:239-276`), returning the
/// source's `ConversionError` message as the error.
///
/// Whitespace is skipped on both sides and `nan`, optionally followed by a
/// parenthesised payload, is accepted before anything else. Then one `+` may
/// precede a `std::from_chars` general-format number: an optional `-`, digits
/// with an optional fraction and exponent, or `inf`/`infinity`. A finite
/// literal that overflows is an error, as `result_out_of_range` is; underflow
/// rounds, as the oracle's libc++ fallback accepts it. Hexadecimal floats are
/// rejected, following `std::from_chars` rather than the `strtod` fallback the
/// source uses on libc++.
pub(crate) fn to_double(text: &str) -> std::result::Result<f64, String> {
    let bytes = text.as_bytes();
    let not_converted = || format!("Could not convert string '{text}' to a double value");
    let first = skip_whitespace(bytes, 0);
    if first == bytes.len() {
        return Err(not_converted());
    }
    if let Some(length) = word_at(bytes, first, "nan") {
        let mut end = first + length;
        if bytes.get(end) == Some(&b'(') {
            match bytes[end..].iter().position(|b| *b == b')') {
                Some(close) => end += close + 1,
                None => end = first,
            }
        }
        if end != first && skip_whitespace(bytes, end) == bytes.len() {
            return Ok(f64::NAN);
        }
    }
    let mut start = first;
    if bytes[start] == b'+' {
        start += 1;
    }
    let mut cursor = start;
    if bytes.get(cursor) == Some(&b'-') {
        cursor += 1;
    }
    let infinity = word_at(bytes, cursor, "infinity").or_else(|| word_at(bytes, cursor, "inf"));
    let nan = word_at(bytes, cursor, "nan");
    if let Some(length) = infinity.or(nan) {
        cursor += length;
    } else {
        let integer_start = cursor;
        while bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
            cursor += 1;
        }
        let mut digits = cursor - integer_start;
        if bytes.get(cursor) == Some(&b'.') {
            let fraction_start = cursor + 1;
            let mut fraction_end = fraction_start;
            while bytes.get(fraction_end).is_some_and(u8::is_ascii_digit) {
                fraction_end += 1;
            }
            if digits > 0 || fraction_end > fraction_start {
                digits += fraction_end - fraction_start;
                cursor = fraction_end;
            }
        }
        if digits == 0 {
            return Err(not_converted());
        }
        if matches!(bytes.get(cursor), Some(b'e' | b'E')) {
            let mut exponent = cursor + 1;
            if matches!(bytes.get(exponent), Some(b'+' | b'-')) {
                exponent += 1;
            }
            let exponent_digits = exponent;
            while bytes.get(exponent).is_some_and(u8::is_ascii_digit) {
                exponent += 1;
            }
            if exponent > exponent_digits {
                cursor = exponent;
            }
        }
    }
    let literal = text.get(start..cursor).ok_or_else(not_converted)?;
    let value: f64 = if nan.is_some() && infinity.is_none() {
        f64::NAN
    } else {
        literal.parse().map_err(|_| not_converted())?
    };
    if value.is_infinite() && infinity.is_none() {
        return Err(not_converted());
    }
    let after = skip_whitespace(bytes, cursor);
    if after != bytes.len() {
        return Err(format!(
            "Prefix of string '{text}' successfully converted to a double value. Additional characters found at position {}",
            after + 1
        ));
    }
    Ok(value)
}
