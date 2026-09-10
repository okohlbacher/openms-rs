// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Nonmutating prefix/suffix searches and source C-locale ASCII case conversion.

use super::list::trim;
use crate::{Error, Result};
use std::ops::Range;

pub fn search_prefix(values: &[impl AsRef<str>], query: &str, trim_values: bool) -> Option<usize> {
    let query = if trim_values { trim(query) } else { query };
    values.iter().position(|value| {
        let value = value.as_ref();
        (if trim_values { trim(value) } else { value }).starts_with(query)
    })
}
pub fn search_suffix(values: &[impl AsRef<str>], query: &str, trim_values: bool) -> Option<usize> {
    let query = if trim_values { trim(query) } else { query };
    values.iter().position(|value| {
        let value = value.as_ref();
        (if trim_values { trim(value) } else { value }).ends_with(query)
    })
}
/// A checked iterator-range equivalent; the returned index is relative to the
/// original full slice. Empty ranges succeed with None.
pub fn search_prefix_in(
    values: &[impl AsRef<str>],
    range: Range<usize>,
    query: &str,
    trim_values: bool,
) -> Result<Option<usize>> {
    let slice = values
        .get(range.clone())
        .ok_or_else(|| Error::InvalidValue("invalid string-list search range".into()))?;
    Ok(search_prefix(slice, query, trim_values).map(|i| i + range.start))
}
pub fn search_suffix_in(
    values: &[impl AsRef<str>],
    range: Range<usize>,
    query: &str,
    trim_values: bool,
) -> Result<Option<usize>> {
    let slice = values
        .get(range.clone())
        .ok_or_else(|| Error::InvalidValue("invalid string-list search range".into()))?;
    Ok(search_suffix(slice, query, trim_values).map(|i| i + range.start))
}
pub fn to_upper(values: &mut [String]) {
    for value in values {
        value.make_ascii_uppercase();
    }
}
pub fn to_lower(values: &mut [String]) {
    for value in values {
        value.make_ascii_lowercase();
    }
}
