// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Native ListUtils operations, including its four defined conversion types.

use crate::{Error, Result};
use std::borrow::Cow;

pub const MAX_ITEMS: usize = 1_000_000;
pub const MAX_BYTES: usize = 64 * 1024 * 1024;
pub const DEFAULT_FLOAT_TOLERANCE: f64 = 0.00001;

fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
pub(crate) fn trim(text: &str) -> &str {
    text.trim_matches([' ', '\t', '\n', '\r'])
}
fn check_bytes(bytes: usize) -> Result<()> {
    if bytes > MAX_BYTES {
        Err(invalid("list byte limit exceeded"))
    } else {
        Ok(())
    }
}
fn storage<T>(count: usize, text_bytes: usize) -> Result<()> {
    if count > MAX_ITEMS {
        return Err(invalid("list item limit exceeded"));
    }
    let bytes = count
        .checked_mul(std::mem::size_of::<T>())
        .and_then(|n| n.checked_mul(2))
        .and_then(|n| {
            text_bytes
                .checked_mul(2)
                .and_then(|text| n.checked_add(text))
        })
        .ok_or_else(|| invalid("list storage size overflow"))?;
    check_bytes(bytes)
}

/// The source provides conversions for String, i32, f32 and f64. Clients may
/// implement this trait for additional owned types, like C++ specializations.
pub trait ListParse: Sized {
    fn from_list_item(text: &str) -> Result<Self>;
}
impl ListParse for String {
    fn from_list_item(text: &str) -> Result<Self> {
        check_bytes(text.len())?;
        Ok(text.to_owned())
    }
}
impl ListParse for i32 {
    fn from_list_item(text: &str) -> Result<Self> {
        check_bytes(text.len())?;
        let token = trim(text);
        // StringUtils removes one leading '+', then uses from_chars. This
        // preserves its unusual acceptance of '+-1', but never '++1'.
        let token = token.strip_prefix('+').unwrap_or(token);
        if token.starts_with('+') {
            return Err(invalid("invalid i32 list item"));
        }
        token.parse().map_err(|_| invalid("invalid i32 list item"))
    }
}
fn special_float(text: &str) -> Option<f64> {
    // The source's explicit unsigned NaN helper accepts arbitrary payload text
    // through the first closing parenthesis. The signed from_chars path permits
    // only the portable alphanumeric/underscore payload grammar.
    if text.eq_ignore_ascii_case("nan") {
        return Some(f64::NAN);
    }
    if text.len() >= 5
        && text
            .get(..3)
            .is_some_and(|prefix| prefix.eq_ignore_ascii_case("nan"))
        && text.as_bytes().get(3) == Some(&b'(')
        && text.ends_with(')')
        && !text[4..text.len() - 1].contains(')')
    {
        return Some(f64::NAN);
    }
    None
}
fn float_token(text: &str) -> Result<(&str, Option<f64>)> {
    check_bytes(text.len())?;
    let text = trim(text);
    if let Some(value) = special_float(text) {
        return Ok((text, Some(value)));
    }
    let token = text.strip_prefix('+').unwrap_or(text);
    if token.starts_with('+') {
        return Err(invalid("invalid floating-point list item"));
    }
    let unsigned = token.strip_prefix('-').unwrap_or(token);
    let sign = if token.starts_with('-') { -1.0 } else { 1.0 };
    if unsigned.eq_ignore_ascii_case("inf") || unsigned.eq_ignore_ascii_case("infinity") {
        return Ok((token, Some(sign * f64::INFINITY)));
    }
    let nan = unsigned.eq_ignore_ascii_case("nan")
        || (unsigned.len() >= 5
            && unsigned
                .get(..4)
                .is_some_and(|s| s.eq_ignore_ascii_case("nan("))
            && unsigned.ends_with(')')
            && unsigned[4..unsigned.len() - 1]
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b == b'_'));
    if nan {
        return Ok((token, Some(sign * f64::NAN)));
    }
    Ok((token, None))
}
fn lost_nonzero_mantissa(token: &str, value: f64) -> bool {
    value == 0.0
        && token
            .split(['e', 'E'])
            .next()
            .unwrap_or(token)
            .bytes()
            .any(|byte| matches!(byte, b'1'..=b'9'))
}
impl ListParse for f64 {
    fn from_list_item(text: &str) -> Result<Self> {
        let (token, special) = float_token(text)?;
        if let Some(value) = special {
            return Ok(value);
        }
        token
            .parse::<f64>()
            .ok()
            .filter(|value| value.is_finite() && !lost_nonzero_mantissa(token, *value))
            .ok_or_else(|| invalid("invalid, overflowing or underflowing f64 list item"))
    }
}
impl ListParse for f32 {
    fn from_list_item(text: &str) -> Result<Self> {
        let (token, special) = float_token(text)?;
        if let Some(value) = special {
            return Ok(value as f32);
        }
        // Direct binary32 parsing avoids an intermediate double-rounding step.
        token
            .parse::<f32>()
            .ok()
            .filter(|value| value.is_finite() && !lost_nonzero_mantissa(token, f64::from(*value)))
            .ok_or_else(|| invalid("invalid, overflowing or underflowing f32 list item"))
    }
}

/// Literal single-byte split, preserving whitespace for String values. Empty
/// text produces an empty vector; nonempty text retains empty split fields.
pub fn create<T: ListParse>(text: &str, separator: u8) -> Result<Vec<T>> {
    if text.is_empty() {
        return Ok(Vec::new());
    }
    check_bytes(text.len())?;
    let count = text.bytes().filter(|&b| b == separator).count() + 1;
    storage::<T>(count, text.len())?;
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| invalid("list allocation failed"))?;
    for bytes in text.as_bytes().split(|&b| b == separator) {
        let item = std::str::from_utf8(bytes).map_err(|_| {
            Error::Unsupported("list byte delimiter splits a UTF-8 character".into())
        })?;
        result.push(T::from_list_item(item)?);
    }
    Ok(result)
}
/// The vector overload preserves String elements verbatim; numeric elements
/// use the same checked, ASCII-trimmed conversion as the text overload.
pub fn create_from_strings<T: ListParse>(values: &[impl AsRef<str>]) -> Result<Vec<T>> {
    storage::<T>(values.len(), 0)?;
    let bytes = values
        .iter()
        .try_fold(0usize, |sum, value| sum.checked_add(value.as_ref().len()))
        .ok_or_else(|| invalid("list byte size overflow"))?;
    storage::<T>(values.len(), bytes)?;
    let mut result = Vec::new();
    result
        .try_reserve_exact(values.len())
        .map_err(|_| invalid("list allocation failed"))?;
    for value in values {
        result.push(T::from_list_item(value.as_ref())?);
    }
    Ok(result)
}

/// Source StringUtils formatting for list items. Borrowed text avoids an
/// unnecessary intermediate copy; numeric and parameter values return owned text.
pub trait ListFormat {
    fn to_list_text(&self) -> Result<Cow<'_, str>>;
}
impl ListFormat for str {
    fn to_list_text(&self) -> Result<Cow<'_, str>> {
        Ok(Cow::Borrowed(self))
    }
}
impl ListFormat for String {
    fn to_list_text(&self) -> Result<Cow<'_, str>> {
        Ok(Cow::Borrowed(self))
    }
}
impl<T: ListFormat + ?Sized> ListFormat for &T {
    fn to_list_text(&self) -> Result<Cow<'_, str>> {
        (*self).to_list_text()
    }
}
macro_rules! integer_text {
    ($($t:ty),*) => {$ (impl ListFormat for $t {
        fn to_list_text(&self) -> Result<Cow<'_, str>> { Ok(Cow::Owned(self.to_string())) }
    })*};
}
integer_text!(i8, u8, i16, u16, i32, u32, i64, u64, isize, usize);
impl ListFormat for bool {
    fn to_list_text(&self) -> Result<Cow<'_, str>> {
        Ok(Cow::Borrowed(if *self { "1" } else { "0" }))
    }
}
impl ListFormat for char {
    fn to_list_text(&self) -> Result<Cow<'_, str>> {
        if !self.is_ascii() {
            return Err(Error::Unsupported(
                "source char formatting is single-byte ASCII in UTF-8".into(),
            ));
        }
        Ok(Cow::Owned(self.to_string()))
    }
}
impl ListFormat for f64 {
    fn to_list_text(&self) -> Result<Cow<'_, str>> {
        Ok(Cow::Owned(crate::param::value::format_float(*self, true)))
    }
}
impl ListFormat for f32 {
    fn to_list_text(&self) -> Result<Cow<'_, str>> {
        Ok(Cow::Owned(crate::param::value::format_float32(*self, true)))
    }
}
impl ListFormat for crate::param::ParamValue {
    fn to_list_text(&self) -> Result<Cow<'_, str>> {
        self.to_text(true).map(Cow::Owned)
    }
}
impl ListFormat for crate::metadata::MetaValue {
    fn to_list_text(&self) -> Result<Cow<'_, str>> {
        use crate::metadata::MetaValueData;
        let content = match self.data() {
            MetaValueData::Empty => return Ok(Cow::Borrowed("")),
            MetaValueData::String(value) => return Ok(Cow::Borrowed(value)),
            MetaValueData::Integer(value) => return value.to_list_text(),
            MetaValueData::Float(value) => return value.to_list_text(),
            MetaValueData::StringList(values) => concatenate(values, ", ")?,
            MetaValueData::IntegerList(values) => concatenate(values, ", ")?,
            MetaValueData::FloatList(values) => concatenate(values, ", ")?,
        };
        check_bytes(
            content
                .len()
                .checked_add(2)
                .ok_or_else(|| invalid("list size overflow"))?,
        )?;
        Ok(Cow::Owned(format!("[{content}]")))
    }
}

/// Any finite iterator can replace the source vector/container overloads.
/// Output size and item count are checked incrementally before retained copies.
pub fn to_string_list(values: impl IntoIterator<Item = impl ListFormat>) -> Result<Vec<String>> {
    let mut output = Vec::new();
    let mut bytes = 0usize;
    for value in values {
        storage::<String>(output.len() + 1, bytes)?;
        let text = value.to_list_text()?;
        let retained = match &text {
            Cow::Borrowed(value) => value.len(),
            Cow::Owned(value) => value.capacity(),
        };
        bytes = bytes
            .checked_add(retained)
            .ok_or_else(|| invalid("list byte size overflow"))?;
        storage::<String>(output.len() + 1, bytes)?;
        // Vec's minimum amortized allocation is four String slots. Start with
        // one so the two-slots-per-item budget also bounds the first result.
        let reserved = if output.capacity() == 0 {
            output.try_reserve_exact(1)
        } else {
            output.try_reserve(1)
        };
        reserved.map_err(|_| invalid("list allocation failed"))?;
        output.push(text.into_owned());
    }
    Ok(output)
}
pub fn concatenate(
    values: impl IntoIterator<Item = impl ListFormat>,
    glue: &str,
) -> Result<String> {
    let items = to_string_list(values)?;
    if items.is_empty() {
        return Ok(String::new());
    }
    let payload: usize = items.iter().map(String::len).sum(); // already bounded by to_string_list
    let separators = glue
        .len()
        .checked_mul(items.len() - 1)
        .ok_or_else(|| invalid("list glue size overflow"))?;
    let bytes = payload
        .checked_add(separators)
        .ok_or_else(|| invalid("list join size overflow"))?;
    storage::<String>(items.len(), bytes)?;
    let mut result = String::new();
    result
        .try_reserve_exact(bytes)
        .map_err(|_| invalid("list join allocation failed"))?;
    for (i, item) in items.iter().enumerate() {
        if i != 0 {
            result.push_str(glue);
        }
        result.push_str(item);
    }
    Ok(result)
}

pub fn contains<T, E: ?Sized>(values: &[T], element: &E) -> bool
where
    T: PartialEq<E>,
{
    values.iter().any(|value| value == element)
}
/// First exact-match index. None replaces the source -1 sentinel and avoids its
/// narrowing cast to a signed 32-bit index.
pub fn get_index<T, E: ?Sized>(values: &[T], element: &E) -> Option<usize>
where
    T: PartialEq<E>,
{
    values.iter().position(|value| value == element)
}
pub fn contains_approx(values: &[f64], element: f64, tolerance: f64) -> bool {
    values
        .iter()
        .any(|&value| (value - element).abs() < tolerance)
}
pub fn contains_f64(values: &[f64], element: f64) -> bool {
    contains_approx(values, element, DEFAULT_FLOAT_TOLERANCE)
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CaseSensitivity {
    Sensitive,
    Insensitive,
}
pub fn contains_string(values: &[impl AsRef<str>], element: &str, case: CaseSensitivity) -> bool {
    values.iter().any(|value| match case {
        CaseSensitivity::Sensitive => value.as_ref() == element,
        CaseSensitivity::Insensitive => value.as_ref().eq_ignore_ascii_case(element),
    })
}
