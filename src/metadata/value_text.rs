// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The source text of a [`MetaValue`]: `DataValue::toString(bool)`, which
//! `StringUtils::toStr(const DataValue&, bool)` and
//! `StringUtils::appendToStr(const DataValue&, bool, std::string&)` return.
//!
//! [`MetaValue`]'s `Display` is the lenient stringification with Rust's own
//! number formatting; [`MetaValue::to_text`] is the source's text, number for
//! number. Sources at core `bc9cc12`: `DATASTRUCTURES/DataValue.cpp`
//! (`toString`), `DATASTRUCTURES/StringUtils.cpp:447-455` and
//! `DATASTRUCTURES/ListUtilsIO.h` (the vector and `VecLowPrecision` stream
//! operators). Tier 1 against `../oracle/B8a-stringutils`
//! (`tests/string_utils.rs`).

use super::{MetaValue, MetaValueData};
use crate::data_structures::string_utils::{MAX_TEXT_BYTES, ToStr, append_to_str_low_p};
use crate::error::{Error, Result};

/// Upper bound of one low-precision element's text: sign, four digits, point,
/// three fraction digits and `e-324`.
const LOW_PRECISION_MAX_BYTES: usize = 16;

impl MetaValue {
    /// `DataValue::toString(full_precision)`: the value's source text.
    ///
    /// Empty is the empty text, a string is itself, an integer its decimal
    /// digits, a float `StringUtils::toStr(value, full_precision)` (see
    /// [`to_str`](crate::data_structures::string_utils::to_str) and
    /// [`to_str_low_p`](crate::data_structures::string_utils::to_str_low_p)).
    /// A list is `[`, the elements' `toStr` texts separated by `, `, then `]`
    /// (`ListUtilsIO.h`); only a float list honours `full_precision = false`
    /// (`VecLowPrecision`). String elements are neither quoted nor escaped, so
    /// `["a", "b, c", ""]` prints `[a, b, c, ]`, as in the source. The unit is
    /// not part of the text.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a list's text could exceed
    /// [`MAX_TEXT_BYTES`] or cannot be allocated. The source has no ceiling.
    ///
    /// # Examples
    ///
    /// ```
    /// use openms::metadata::MetaValue;
    ///
    /// let value = MetaValue::try_from(12345.6789)?;
    /// assert_eq!(value.to_text(true)?, "1.23456789e04");
    /// assert_eq!(value.to_text(false)?, "1.235e04");
    /// let list = MetaValue::try_from(vec![1.5, 0.0625])?;
    /// assert_eq!(list.to_text(false)?, "[1.5, 0.062]");
    /// # Ok::<(), openms::Error>(())
    /// ```
    pub fn to_text(&self, full_precision: bool) -> Result<String> {
        match self.data() {
            MetaValueData::Empty => Ok(String::new()),
            MetaValueData::String(value) => Ok(value.clone()),
            MetaValueData::Integer(value) => Ok(value.to_str_text()),
            MetaValueData::Float(value) => Ok(if full_precision {
                value.to_str_text()
            } else {
                let mut text = String::new();
                append_to_str_low_p(*value, &mut text);
                text
            }),
            MetaValueData::StringList(values) => list_text(
                values,
                |value, text| {
                    text.push_str(value);
                },
                |value| value.len(),
            ),
            MetaValueData::IntegerList(values) => list_text(
                values,
                |value, text| value.append_to_str(text),
                ToStr::text_len_bound,
            ),
            MetaValueData::FloatList(values) if full_precision => list_text(
                values,
                |value, text| value.append_to_str(text),
                ToStr::text_len_bound,
            ),
            MetaValueData::FloatList(values) => list_text(
                values,
                |value, text| append_to_str_low_p(*value, text),
                |_| LOW_PRECISION_MAX_BYTES,
            ),
        }
    }
}

/// The `ListUtilsIO.h` vector text, its size checked before allocating.
fn list_text<T>(
    values: &[T],
    append: impl Fn(&T, &mut String),
    bound: impl Fn(&T) -> usize,
) -> Result<String> {
    let mut total = 2usize;
    for (index, value) in values.iter().enumerate() {
        let separator = if index == 0 { 0 } else { 2 };
        total = total
            .checked_add(separator)
            .and_then(|sum| sum.checked_add(bound(value)))
            .filter(|sum| *sum <= MAX_TEXT_BYTES)
            .ok_or_else(|| {
                Error::InvalidValue(format!(
                    "the text of a {}-element list may exceed {MAX_TEXT_BYTES} bytes",
                    values.len()
                ))
            })?;
    }
    let mut text = String::new();
    text.try_reserve_exact(total)
        .map_err(|_| Error::InvalidValue(format!("cannot allocate {total} bytes of list text")))?;
    text.push('[');
    for (index, value) in values.iter().enumerate() {
        if index > 0 {
            text.push_str(", ");
        }
        append(value, &mut text);
    }
    text.push(']');
    Ok(text)
}
