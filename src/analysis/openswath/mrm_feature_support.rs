// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Metavalue access shared by the MRM feature filtering, selection and
//! confidence scoring modules.
//!
//! These are the `DataValue` conversions those headers lean on:
//! `StringUtils::toStr(const DataValue&)`, which is lenient and never throws,
//! and the `(double)` cast, which accepts an integer or a float. The source's
//! cast of any other type reads the `DataValue` union as a double, which is
//! undefined behaviour; this port reports it.

use crate::kernel::Feature;
use crate::metadata::{MetaInfo, MetaValue, MetaValueData};
use crate::{Error, Result};

/// `StringUtils::toStr(DataValue)` over an optional value: the empty string
/// for a missing or empty value, a string verbatim, an integer in decimal, a
/// float as `StringUtils::toStr(double)` and a list as `[a, b]`.
pub(crate) fn data_value_to_string(value: Option<&MetaValue>) -> String {
    fn joined<T>(values: &[T], element: impl Fn(&T) -> String) -> String {
        let mut out = String::from("[");
        for (index, value) in values.iter().enumerate() {
            if index != 0 {
                out.push_str(", ");
            }
            out.push_str(&element(value));
        }
        out.push(']');
        out
    }
    let Some(value) = value else {
        return String::new();
    };
    match value.data() {
        MetaValueData::Empty => String::new(),
        MetaValueData::String(text) => text.clone(),
        MetaValueData::Integer(number) => number.to_string(),
        MetaValueData::Float(number) => crate::param::value::format_float(*number, true),
        MetaValueData::StringList(list) => joined(list, Clone::clone),
        MetaValueData::IntegerList(list) => joined(list, ToString::to_string),
        MetaValueData::FloatList(list) => {
            joined(list, |x| crate::param::value::format_float(*x, true))
        }
    }
}

/// The lenient text of a metavalue of `metadata` (empty when absent).
pub(crate) fn meta_string(metadata: &MetaInfo, key: &str) -> String {
    data_value_to_string(metadata.get(key))
}

/// `(double)getMetaValue(key)`: an integer or a float, as a double.
///
/// # Errors
///
/// [`Error::InvalidValue`] for an empty value, a missing key and any other
/// type (the source throws `ConversionError` for an empty value and reads an
/// unrelated union member for the rest).
pub(crate) fn meta_double(metadata: &MetaInfo, key: &str) -> Result<f64> {
    match metadata.get(key).map(MetaValue::data) {
        Some(MetaValueData::Integer(number)) => Ok(*number as f64),
        Some(MetaValueData::Float(number)) => Ok(*number),
        Some(MetaValueData::Empty) | None => Err(Error::InvalidValue(format!(
            "metavalue '{key}' is missing or empty and cannot be converted to a double"
        ))),
        Some(_) => Err(Error::InvalidValue(format!(
            "metavalue '{key}' is neither an integer nor a float and cannot be converted to a double"
        ))),
    }
}

/// `native_id` of a feature as the filters read it.
pub(crate) fn native_id(feature: &Feature) -> String {
    meta_string(&feature.metadata, "native_id")
}
