// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Parameters from a JSON document onto a known parameter tree, the load side
//! of `ParamJSONFile` (`OpenMS/FORMAT/ParamJSONFile.h`).
//!
//! A TOPP tool reads an `-ini` file whose name ends in `.json` through this
//! module (`TOPPBase.cpp:283-296` and `2619-2629` at the cli pin). The JSON
//! document carries no tree of its own: the caller hands in the tool's default
//! parameters, and the document only sets values of entries that already
//! exist. The source's store side needs the TDL library, which the Release
//! build is compiled without: there `ParamJSONFile::store` and
//! `writeToStream` throw, and so do
//! [`ParamJsonFile::store`](crate::format::param_json::ParamJsonFile::store) and
//! [`ParamJsonFile::write_to_stream`](crate::format::param_json::ParamJsonFile::write_to_stream)
//! here. The command line refuses
//! `-write_json` for the same reason (see `src/cli.rs`).
//!
//! Whether the document is JSON, and the exception text when it is not, is
//! `nlohmann::json` 3.11.2's (`src/openms/extern/nlohmann_json` at the same
//! pin), through a port of its scanner and parser; the values are read by
//! `serde_json`, and every rule after the parse is the source's
//! (`ParamJSONFile.cpp:38-167` at core pin `bc9cc12`), including the
//! conversions of the `nlohmann::json` accessors it calls and their error
//! texts. See `docs/TOPP_CLI_SUPPORT.md` ("JSON INI files").

use crate::data_structures::ToolInfo;
use crate::param::{MAX_PARAM_BYTES, Param, ParamValue};
use crate::{Error, Result};
use serde_json::Value;
use std::io::{Read, Write};
use std::path::Path;

mod nlohmann;

/// The message of the `std::runtime_error` the Release build's
/// `ParamJSONFile::writeToStream` throws, built without TDL
/// (`ParamJSONFile.cpp:323-327` at core `bc9cc12`).
pub const TDL_UNAVAILABLE: &str =
    "TDL support is not available. Rebuild with -DENABLE_TDL=ON to enable this feature.";

/// Source class `ParamJSONFile`: [`load`] reads, and the store side is the
/// Release build's, which has no TDL.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ParamJsonFile {
    /// Source `flatHierarchy`: with TDL, names would be written with their
    /// sections joined by `__` instead of nested. Kept for the API; the
    /// Release build writes nothing either way.
    pub flat_hierarchy: bool,
}

impl ParamJsonFile {
    /// Source `load` ([`load`]), a static member there.
    ///
    /// # Errors
    ///
    /// As [`load`].
    pub fn load(path: impl AsRef<Path>, param: &mut Param) -> Result<Vec<String>> {
        load(path, param)
    }

    /// Source `store` (`ParamJSONFile.cpp:170-189`): in the Release build the
    /// writer throws [`TDL_UNAVAILABLE`] after the target was opened, which
    /// leaves an empty file behind; this port refuses before it touches
    /// `filename` (CPP-350, as the command line's `-write_json`).
    /// `tool_info` is unused, as in the source.
    ///
    /// # Errors
    ///
    /// Always [`Error::Unsupported`] with [`TDL_UNAVAILABLE`].
    pub fn store(&self, filename: &str, param: &Param, tool_info: &ToolInfo) -> Result<()> {
        let _ = (filename, tool_info);
        self.write_to_stream(&mut std::io::sink(), param)
    }

    /// Source `writeToStream`, which the Release build compiles to the throw of
    /// [`TDL_UNAVAILABLE`] (`ParamJSONFile.cpp:323-327`). Nothing is written.
    ///
    /// # Errors
    ///
    /// Always [`Error::Unsupported`] with [`TDL_UNAVAILABLE`].
    pub fn write_to_stream(&self, output: &mut dyn Write, param: &Param) -> Result<()> {
        let _ = (output, param, self.flat_hierarchy);
        Err(Error::Unsupported(TDL_UNAVAILABLE.into()))
    }
}

/// The largest JSON document [`load`] reads, in bytes: the bound every
/// parameter tree of this crate already has ([`MAX_PARAM_BYTES`]).
pub const MAX_JSON_BYTES: usize = MAX_PARAM_BYTES;

/// Source `ParamJSONFile::load`: sets the values `path` names on `param`.
///
/// The document must be a JSON object. Its keys are the parameter names below
/// `<tool>:1:`, where `<tool>` is the first section of `param`
/// (`ParamJSONFile.cpp:41-43`): the instance is always `1`, whatever instance
/// `param` was built for. `__` in a key stands for `:`, and a nested object is
/// a section, unless it has a string member `class`, which makes it a file
/// value such as `{"class": "File", "path": "in.mzML"}`. A `null` value is
/// skipped.
///
/// Returns the warning lines the source writes to `OPENMS_LOG_WARN`, one for
/// each entry of `param` that holds no value (`ParamJSONFile.cpp:152-156`).
///
/// # Errors
///
/// * [`Error::Io`] when the file cannot be opened or read. The source raises
///   `FileNotFound`, `FileNotReadable` or `IOException` there
///   (`ParamJSONFile.cpp:46-60`); the caller tells them apart by the kind.
/// * [`Error::Parse`] for everything the source raises as `ParseError`
///   (`ParamJSONFile.cpp:66-69`, `86-89` and `162-165`): a document that is no
///   JSON or no object, a name `param` does not have, and a value of the wrong
///   JSON type. The message is the source's; its `what()` appends ` in: `.
/// * [`Error::InvalidValue`] for a document larger than [`MAX_JSON_BYTES`], or
///   an empty `param`, which has no first section to name the tool.
///
/// `param` is changed only when the load succeeds. The source sets values as
/// it walks the document and leaves the ones before a failure in place.
pub fn load(path: impl AsRef<Path>, param: &mut Param) -> Result<Vec<String>> {
    let bytes = read_bounded(path.as_ref())?;
    load_from_slice(&bytes, &path.as_ref().display().to_string(), param)
}

/// [`load`] for a document already in memory; `filename` is the name the
/// source's messages quote.
pub fn load_from_slice(bytes: &[u8], filename: &str, param: &mut Param) -> Result<Vec<String>> {
    if bytes.len() > MAX_JSON_BYTES {
        return Err(Error::InvalidValue(format!(
            "JSON parameter file '{filename}' exceeds {MAX_JSON_BYTES} bytes"
        )));
    }
    let tool = match param
        .iter()?
        .next()
        .and_then(|item| item.trace.first().cloned())
    {
        Some(trace) => trace.name,
        None => {
            return Err(Error::InvalidValue(
                "a JSON parameter file needs a parameter tree with a first section".into(),
            ));
        }
    };
    let document = parse(bytes)?;
    let mut loaded = param.checked_clone()?;
    let mut walk = Walk {
        filename,
        tool: &tool,
        param: &mut loaded,
        warnings: Vec::new(),
    };
    walk.object(&format!("{tool}:1:"), &document)?;
    let warnings = walk.warnings;
    *param = loaded;
    Ok(warnings)
}

fn read_bounded(path: &Path) -> Result<Vec<u8>> {
    let file = std::fs::File::open(path)?;
    if let Ok(metadata) = file.metadata() {
        if metadata.is_file() && metadata.len() > MAX_JSON_BYTES as u64 {
            return Err(Error::InvalidValue(format!(
                "JSON parameter file '{}' exceeds {MAX_JSON_BYTES} bytes",
                path.display()
            )));
        }
    }
    let mut bytes = Vec::new();
    // One byte past the bound, so a stream that is too long is refused rather
    // than cut.
    file.take(MAX_JSON_BYTES as u64 + 1)
        .read_to_end(&mut bytes)?;
    Ok(bytes)
}

/// `json::parse` (`ParamJSONFile.cpp:63`): `nlohmann::json` 3.11.2 decides
/// whether the document is JSON and words the exception when it is not
/// ([`nlohmann::check`]: `parse_error.101` with the source's line, column and
/// message, or `out_of_range.406` for a number beyond `double`); `serde_json`
/// then reads the values, with what `nlohmann::json` does differently put
/// back:
///
/// * a UTF-8 byte order mark is skipped, where `serde_json` refuses it;
/// * a NUL byte after the value ends the document, as the source's scanner
///   takes it for the end of input;
/// * `-0` is the integer 0, where `serde_json` reads the float `-0.0`
///   ([`negative_zero_as_zero`]).
///
/// A byte of the source's `last read: '...'` that is not UTF-8 is written as
/// U+FFFD, since the message is text here.
fn parse(bytes: &[u8]) -> Result<Value> {
    let end = nlohmann::check(bytes).map_err(|what| parse_error(String::from_utf8_lossy(&what)))?;
    let bytes = &bytes[..end];
    let bytes = bytes.strip_prefix(b"\xEF\xBB\xBF").unwrap_or(bytes);
    let normalised = negative_zero_as_zero(bytes);
    let bytes = normalised.as_deref().unwrap_or(bytes);
    serde_json::from_slice(bytes).map_err(|error| Error::Parse {
        line: error.line(),
        message: error.to_string(),
    })
}

/// The document with every number token `-0` written ` 0`, or `None` when it
/// has none. `nlohmann::json` reads `-0` as the integer 0 and `-0.0` as the
/// float -0.0; `serde_json` reads both as -0.0, which a floating-point
/// parameter would keep. The space keeps every later byte in place. A `-`
/// inside a string or after an exponent's `e` is left alone.
fn negative_zero_as_zero(bytes: &[u8]) -> Option<Vec<u8>> {
    let mut copy: Option<Vec<u8>> = None;
    let (mut in_string, mut escaped) = (false, false);
    for (index, &byte) in bytes.iter().enumerate() {
        if in_string {
            if escaped {
                escaped = false;
            } else if byte == b'\\' {
                escaped = true;
            } else if byte == b'"' {
                in_string = false;
            }
            continue;
        }
        if byte == b'"' {
            in_string = true;
            continue;
        }
        let token = byte == b'-'
            && bytes.get(index + 1) == Some(&b'0')
            && !matches!(bytes.get(index + 2), Some(b'0'..=b'9' | b'.' | b'e' | b'E'))
            && !matches!(index.checked_sub(1).map(|at| bytes[at]), Some(b'e' | b'E'));
        if token {
            copy.get_or_insert_with(|| bytes.to_vec())[index] = b' ';
        }
    }
    copy
}

fn parse_error(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

struct Walk<'a> {
    filename: &'a str,
    tool: &'a str,
    param: &'a mut Param,
    warnings: Vec<String>,
}

impl Walk<'_> {
    /// The source's `traverseJSONTree` (`ParamJSONFile.cpp:64-160`). Members
    /// are visited in the byte order of their names, the order of the
    /// `std::map` behind `nlohmann::json`; a name given twice keeps its last
    /// value, as in both parsers.
    fn object(&mut self, prefix: &str, node: &Value) -> Result<()> {
        let Value::Object(members) = node else {
            return Err(parse_error(format!(
                "Ignoring JSON file '{}' because of unexpected data type. Expecting a dictionary as type.",
                self.filename
            )));
        };
        let mut members: Vec<(&String, &Value)> = members.iter().collect();
        members.sort_by(|a, b| a.0.as_bytes().cmp(b.0.as_bytes()));
        for (name, value) in members {
            let key = format!("{prefix}{}", name.replace("__", ":"));
            if value.is_null() {
                continue;
            }
            if value.is_object() && !value.get("class").is_some_and(Value::is_string) {
                self.object(&format!("{key}:"), value)?;
                continue;
            }
            self.entry(&key, value)?;
        }
        Ok(())
    }

    /// One value (`ParamJSONFile.cpp:86-157`). A name too deep for a parameter
    /// tree cannot be one of `param`'s, so the depth check of
    /// [`Param::exists`] counts as absence here.
    fn entry(&mut self, key: &str, node: &Value) -> Result<()> {
        if !matches!(self.param.exists(key), Ok(true)) {
            return Err(parse_error(format!(
                "Parameter {key} passed to '{}' is invalid. To prevent usage of wrong defaults, please update/fix the parameters!",
                self.tool
            )));
        }
        let entry = self.param.entry(key)?;
        let input_file = entry.tags.contains("input file");
        let value = match &entry.value {
            ParamValue::String(_) => {
                let flag = entry.valid_strings.len() == 2
                    && ((entry.valid_strings[0] == "true" && entry.valid_strings[1] == "false")
                        || (entry.valid_strings[0] == "false" && entry.valid_strings[1] == "true"));
                if flag {
                    ParamValue::String(if boolean(node)? { "true" } else { "false" }.into())
                } else if input_file {
                    if entry.tags.contains("is_executable") && !node.is_object() {
                        ParamValue::String(string(node)?)
                    } else {
                        ParamValue::String(string(path(node)?)?)
                    }
                } else {
                    ParamValue::String(string(node)?)
                }
            }
            ParamValue::Integer(_) => ParamValue::Integer(int64(node)?),
            ParamValue::Float(_) => ParamValue::Float(double(node)?),
            ParamValue::StringList(_) => {
                let node = if input_file { path(node)? } else { node };
                ParamValue::StringList(array(node)?.iter().map(string).collect::<Result<_>>()?)
            }
            ParamValue::IntegerList(_) => {
                ParamValue::IntegerList(array(node)?.iter().map(int).collect::<Result<_>>()?)
            }
            ParamValue::FloatList(_) => {
                ParamValue::FloatList(array(node)?.iter().map(double).collect::<Result<_>>()?)
            }
            ParamValue::Empty => {
                self.warnings.push(format!(
                    "Ignoring entry '{key}' because of unknown type 'EMPTY_VALUE'."
                ));
                ParamValue::Empty
            }
        };
        self.param.set_value(key, value, "", &[])
    }
}

/// `basic_json::type_name`.
fn type_name(node: &Value) -> &'static str {
    match node {
        Value::Null => "null",
        Value::Bool(_) => "boolean",
        Value::Number(_) => "number",
        Value::String(_) => "string",
        Value::Array(_) => "array",
        Value::Object(_) => "object",
    }
}

/// `type_error` 302, which every `get<T>()` of the wrong JSON type raises.
fn type_error(expected: &str, node: &Value) -> Error {
    parse_error(format!(
        "[json.exception.type_error.302] type must be {expected}, but is {}",
        type_name(node)
    ))
}

/// `node["path"]` on the source's copy of a value: an object without the
/// member gives `null`, anything else `type_error` 305.
fn path(node: &Value) -> Result<&Value> {
    match node {
        Value::Object(members) => Ok(members.get("path").unwrap_or(&Value::Null)),
        other => Err(parse_error(format!(
            "[json.exception.type_error.305] cannot use operator[] with a string argument with {}",
            type_name(other)
        ))),
    }
}

fn boolean(node: &Value) -> Result<bool> {
    node.as_bool().ok_or_else(|| type_error("boolean", node))
}

fn string(node: &Value) -> Result<String> {
    node.as_str()
        .map(str::to_owned)
        .ok_or_else(|| type_error("string", node))
}

fn array(node: &Value) -> Result<&Vec<Value>> {
    node.as_array().ok_or_else(|| type_error("array", node))
}

/// `get<int64_t>()`, `get_arithmetic_value` (json.hpp 4610-4640): a number
/// cast with `static_cast`, and no boolean. An unsigned value above
/// `i64::MAX` wraps, as the conversion does on the reference build; a
/// fraction is truncated, and a float outside `i64` gives `i64::MIN`, the
/// x86-64 result of the conversion the language leaves undefined.
fn int64(node: &Value) -> Result<i64> {
    let Value::Number(number) = node else {
        return Err(type_error("number", node));
    };
    if let Some(value) = number.as_i64() {
        Ok(value)
    } else if let Some(value) = number.as_u64() {
        Ok(value as i64)
    } else {
        Ok(float_to_i64(number.as_f64().unwrap_or(f64::NAN)))
    }
}

/// `get<double>()`: a number, and no boolean.
fn double(node: &Value) -> Result<f64> {
    match node {
        Value::Number(number) => Ok(number.as_f64().unwrap_or(f64::NAN)),
        other => Err(type_error("number", other)),
    }
}

/// `get<int>()` for an element of an integer list, the generic arithmetic
/// `from_json` (json.hpp 4891-4934), which, unlike `get<int64_t>()`, takes a
/// boolean as 0 or 1. A value outside `i32` is cut to its low 32 bits, and a
/// float outside `i32` gives `i32::MIN`, as on the reference build.
fn int(node: &Value) -> Result<i32> {
    match node {
        Value::Bool(value) => Ok(i32::from(*value)),
        Value::Number(number) => Ok(if let Some(value) = number.as_i64() {
            value as i32
        } else if let Some(value) = number.as_u64() {
            value as i32
        } else {
            float_to_i32(number.as_f64().unwrap_or(f64::NAN))
        }),
        other => Err(type_error("number", other)),
    }
}

/// `static_cast<int64_t>(double)` as `cvttsd2si` computes it.
fn float_to_i64(value: f64) -> i64 {
    // 2^63 is exact in f64; the range is [-2^63, 2^63).
    const LIMIT: f64 = 9_223_372_036_854_775_808.0;
    if (-LIMIT..LIMIT).contains(&value) {
        value.trunc() as i64
    } else {
        i64::MIN
    }
}

/// `static_cast<int>(double)` as the 32-bit `cvttsd2si` computes it.
fn float_to_i32(value: f64) -> i32 {
    const LIMIT: f64 = 2_147_483_648.0;
    if value.is_finite() && value > -LIMIT - 1.0 && value < LIMIT {
        value.trunc() as i32
    } else {
        i32::MIN
    }
}
