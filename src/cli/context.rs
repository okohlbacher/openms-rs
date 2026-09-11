// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Resolved tool parameters, the native form of the source `get*_` accessors.

use crate::param::{Param, ParamValue};
use crate::{Error, Result};

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

/// Parameters of one tool run, after defaults, INI and command line are merged.
///
/// The source exposes these through protected `getStringOption_` style methods
/// on the tool itself. Here they are a borrowed context passed to `run`, so a
/// tool cannot mutate its own resolved parameters mid-run.
#[derive(Clone, Debug)]
pub struct ToolContext {
    param: Param,
    prefix: String,
    debug_level: i64,
    test_mode: bool,
    no_progress: bool,
    force: bool,
    threads: i64,
}

impl ToolContext {
    pub(crate) fn new(
        param: Param,
        prefix: String,
        debug_level: i64,
        test_mode: bool,
        no_progress: bool,
        force: bool,
        threads: i64,
    ) -> Self {
        Self {
            param,
            prefix,
            debug_level,
            test_mode,
            no_progress,
            force,
            threads,
        }
    }
    /// The complete resolved parameter tree, as `getParam_`.
    pub fn param(&self) -> &Param {
        &self.param
    }
    pub fn debug_level(&self) -> i64 {
        self.debug_level
    }
    pub fn test_mode(&self) -> bool {
        self.test_mode
    }
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

    fn key(&self, name: &str) -> String {
        format!("{}{}", self.prefix, name)
    }
    fn value(&self, name: &str) -> Result<&ParamValue> {
        self.param
            .value(&self.key(name))
            .map_err(|_| bad(format!("parameter '{name}' was not registered")))
    }

    pub fn string(&self, name: &str) -> Result<&str> {
        match self.value(name)? {
            ParamValue::String(text) => Ok(text),
            _ => Err(bad(format!("parameter '{name}' is not a string"))),
        }
    }
    pub fn int(&self, name: &str) -> Result<i64> {
        match self.value(name)? {
            ParamValue::Integer(value) => Ok(*value),
            _ => Err(bad(format!("parameter '{name}' is not an integer"))),
        }
    }
    pub fn double(&self, name: &str) -> Result<f64> {
        match self.value(name)? {
            ParamValue::Float(value) => Ok(*value),
            ParamValue::Integer(value) => Ok(*value as f64),
            _ => Err(bad(format!("parameter '{name}' is not a number"))),
        }
    }
    pub fn string_list(&self, name: &str) -> Result<&[String]> {
        match self.value(name)? {
            ParamValue::StringList(values) => Ok(values),
            _ => Err(bad(format!("parameter '{name}' is not a string list"))),
        }
    }
    pub fn int_list(&self, name: &str) -> Result<&[i32]> {
        match self.value(name)? {
            ParamValue::IntegerList(values) => Ok(values),
            _ => Err(bad(format!("parameter '{name}' is not an integer list"))),
        }
    }
    pub fn double_list(&self, name: &str) -> Result<&[f64]> {
        match self.value(name)? {
            ParamValue::FloatList(values) => Ok(values),
            _ => Err(bad(format!("parameter '{name}' is not a float list"))),
        }
    }
    /// Source `getFlag_`: a registered flag is true only when it was given.
    pub fn flag(&self, name: &str) -> Result<bool> {
        match self.value(name)? {
            ParamValue::String(text) => Ok(text == "true"),
            ParamValue::Integer(value) => Ok(*value != 0),
            _ => Err(bad(format!("parameter '{name}' is not a flag"))),
        }
    }
    /// Values of a registered subsection, with the subsection prefix removed.
    pub fn subsection(&self, name: &str) -> Result<Param> {
        self.param.copy(&format!("{}{name}:", self.prefix), true)
    }
}

/// Source `parseRange_`: `":8"`, `"2:"`, `"2:8"` and `""` are all accepted, and
/// an absent side leaves that bound untouched. Returns whether anything was set.
pub fn parse_range(text: &str, low: &mut f64, high: &mut f64) -> Result<bool> {
    let text = text.trim();
    if text.is_empty() {
        return Ok(false);
    }
    let Some((start, end)) = text.split_once(':') else {
        return Err(bad(format!("range '{text}' needs a colon")));
    };
    let mut set = false;
    if !start.is_empty() {
        *low = start
            .trim()
            .parse()
            .map_err(|_| bad(format!("range start '{start}' is not a number")))?;
        set = true;
    }
    if !end.is_empty() {
        *high = end
            .trim()
            .parse()
            .map_err(|_| bad(format!("range end '{end}' is not a number")))?;
        set = true;
    }
    if set && *low > *high {
        return Err(bad(format!("range '{text}' is empty")));
    }
    Ok(set)
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
