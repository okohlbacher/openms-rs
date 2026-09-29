// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Small helpers that build and read the source `defaults_` trees of the
//! feature-grouping ports, so each port states its parameters in the order and
//! with the restrictions the source's `defaults_.setValue(...)` calls give them.

use crate::error::{Error, Result};
use crate::param::{Param, ParamValue};

/// The tag list for an `{"advanced"}` parameter, or none.
fn tags(advanced: bool) -> Vec<String> {
    if advanced {
        vec!["advanced".to_string()]
    } else {
        Vec::new()
    }
}

/// `defaults_.setValue(key, value, description, tags)` for a float, with the
/// optional `setMinFloat` / `setMaxFloat` restrictions.
pub(crate) fn add_float(
    param: &mut Param,
    key: &str,
    value: f64,
    description: &str,
    advanced: bool,
    min: Option<f64>,
    max: Option<f64>,
) -> Result<()> {
    param.set_value(key, ParamValue::Float(value), description, &tags(advanced))?;
    if let Some(min) = min {
        param.set_min_float(key, min)?;
    }
    if let Some(max) = max {
        param.set_max_float(key, max)?;
    }
    Ok(())
}

/// `defaults_.setValue(key, value, description, tags)` for a string, with the
/// `setValidStrings` restriction when `valid` is not empty.
pub(crate) fn add_string(
    param: &mut Param,
    key: &str,
    value: &str,
    description: &str,
    advanced: bool,
    valid: &[&str],
) -> Result<()> {
    param.set_value(
        key,
        ParamValue::String(value.to_string()),
        description,
        &tags(advanced),
    )?;
    if !valid.is_empty() {
        let valid: Vec<String> = valid.iter().map(|s| (*s).to_string()).collect();
        param.set_valid_strings(key, &valid)?;
    }
    Ok(())
}

/// `defaults_.setValue(key, value, description, tags)` for an integer, with the
/// optional `setMinInt` / `setMaxInt` restrictions.
pub(crate) fn add_int(
    param: &mut Param,
    key: &str,
    value: i64,
    description: &str,
    advanced: bool,
    min: Option<i32>,
    max: Option<i32>,
) -> Result<()> {
    param.set_value(
        key,
        ParamValue::Integer(value),
        description,
        &tags(advanced),
    )?;
    if let Some(min) = min {
        param.set_min_int(key, min)?;
    }
    if let Some(max) = max {
        param.set_max_int(key, max)?;
    }
    Ok(())
}

/// An integer parameter.
pub(crate) fn int(param: &Param, key: &str) -> Result<i64> {
    param.value(key)?.to_i64()
}

/// A float parameter; an integer is converted as the source's `DataValue`
/// conversion to `double` does.
pub(crate) fn float(param: &Param, key: &str) -> Result<f64> {
    param.value(key)?.to_f64()
}

/// A string parameter.
pub(crate) fn string(param: &Param, key: &str) -> Result<String> {
    Ok(param.value(key)?.as_str()?.to_string())
}

/// A `"true"`/`"false"` string parameter, as the source's `toBool()` reads it.
pub(crate) fn flag(param: &Param, key: &str) -> Result<bool> {
    param.value(key)?.to_bool()
}

/// The source's `"true"`/`"false"` spelling of a flag.
pub(crate) fn flag_text(value: bool) -> &'static str {
    if value { "true" } else { "false" }
}

/// An [`Error::InvalidValue`] with `message`.
pub(crate) fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}
