// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! OpenMS CsvFile's literal delimiter parser; deliberately not RFC CSV.

pub use super::text::Limits;
use super::text::{self, TextFile, invalid};
use crate::{Error, Result};
use std::io::{BufRead, Write};
use std::path::Path;

/// Source single-byte delimiter and enclosure options. first_n counts retained
/// lines after hardcoded # comments; zero/negative values mean no early stop.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadOptions {
    pub separator: u8,
    pub item_enclosed: bool,
    pub first_n: i32,
    pub limits: Limits,
}
impl Default for ReadOptions {
    fn default() -> Self {
        Self {
            separator: b',',
            item_enclosed: false,
            first_n: -1,
            limits: Limits::default(),
        }
    }
}
/// A separate cap is needed because a row can contain many zero-length fields.
pub const MAX_FIELDS: usize = 1_000_000;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CsvFile {
    text: TextFile,
    separator: u8,
    item_enclosed: bool,
}
impl Default for CsvFile {
    fn default() -> Self {
        Self {
            text: TextFile::default(),
            separator: b',',
            item_enclosed: false,
        }
    }
}
impl CsvFile {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn with_options(options: &ReadOptions) -> Result<Self> {
        Ok(Self {
            text: TextFile::with_limits(options.limits)?,
            separator: options.separator,
            item_enclosed: options.item_enclosed,
        })
    }
    fn text_options(options: &ReadOptions, trim_lines: bool) -> text::ReadOptions {
        text::ReadOptions {
            trim_lines,
            first_n: options.first_n,
            skip_empty_lines: false,
            comment_symbol: "#".into(),
            limits: options.limits,
        }
    }
    /// Mirrors the source filename constructor: preserve whitespace and test
    /// the hardcoded # comment prefix on the untrimmed line.
    pub fn from_reader(reader: impl BufRead, options: &ReadOptions) -> Result<Self> {
        Ok(Self {
            text: TextFile::from_reader(reader, &Self::text_options(options, false))?,
            separator: options.separator,
            item_enclosed: options.item_enclosed,
        })
    }
    pub fn from_path(path: impl AsRef<Path>, options: &ReadOptions) -> Result<Self> {
        Ok(Self {
            text: TextFile::from_path(path, &Self::text_options(options, false))?,
            separator: options.separator,
            item_enclosed: options.item_enclosed,
        })
    }
    /// Mirrors source load(): trim before the hardcoded # comment test.
    /// Replacement of both contents and settings is atomic on failure.
    pub fn load_reader(&mut self, reader: impl BufRead, options: &ReadOptions) -> Result<()> {
        let text = TextFile::from_reader(reader, &Self::text_options(options, true))?;
        *self = Self {
            text,
            separator: options.separator,
            item_enclosed: options.item_enclosed,
        };
        Ok(())
    }
    pub fn load(&mut self, path: impl AsRef<Path>, options: &ReadOptions) -> Result<()> {
        let text = TextFile::from_path(path, &Self::text_options(options, true))?;
        *self = Self {
            text,
            separator: options.separator,
            item_enclosed: options.item_enclosed,
        };
        Ok(())
    }
    pub fn row_count(&self) -> usize {
        self.text.len()
    }
    pub fn clear(&mut self) {
        self.text.clear();
    }
    pub fn separator(&self) -> u8 {
        self.separator
    }
    pub fn item_enclosed(&self) -> bool {
        self.item_enclosed
    }
    pub fn write(&self, writer: impl Write) -> Result<()> {
        self.text.write(writer)
    }
    pub fn store(&self, path: impl AsRef<Path>) -> Result<()> {
        self.text.store(path)
    }

    /// Literal join, optionally wrapping each item in double quotes without any
    /// escaping. Embedded delimiters/newlines/quotes retain source behavior.
    pub fn add_row<S: AsRef<str>>(&mut self, fields: &[S]) -> Result<()> {
        if fields.len() > MAX_FIELDS {
            return Err(invalid("CSV field count limit exceeded"));
        }
        if fields.len() > 1 && !self.separator.is_ascii() {
            return Err(Error::Unsupported(
                "non-ASCII CSV delimiter would produce invalid UTF-8".into(),
            ));
        }
        let mut size = fields.len().saturating_sub(1);
        for field in fields {
            size = size
                .checked_add(field.as_ref().len())
                .and_then(|n| n.checked_add(if self.item_enclosed { 2 } else { 0 }))
                .ok_or_else(|| invalid("CSV row size overflow"))?;
        }
        self.text.check_append(size)?;
        let mut line = String::new();
        line.try_reserve_exact(size)
            .map_err(|_| invalid("CSV row allocation failed"))?;
        for (index, field) in fields.iter().enumerate() {
            if index != 0 {
                line.push(char::from(self.separator));
            }
            if self.item_enclosed {
                line.push('"');
            }
            line.push_str(field.as_ref());
            if self.item_enclosed {
                line.push('"');
            }
        }
        self.text.push_owned(line)
    }
    /// The bool reports whether a delimiter was found, not whether the row
    /// exists. Empty rows yield (false, []); a single field is returned unchanged,
    /// even with enclosure enabled. Split fields lose their first/last byte.
    pub fn row(&self, row: usize) -> Result<(bool, Vec<String>)> {
        let line = self
            .text
            .lines()
            .get(row)
            .ok_or_else(|| invalid("CSV row index out of bounds"))?;
        if line.is_empty() {
            return Ok((false, Vec::new()));
        }
        let delimiters = line.bytes().filter(|&b| b == self.separator).count();
        let count = delimiters
            .checked_add(1)
            .filter(|&n| n <= MAX_FIELDS)
            .ok_or_else(|| invalid("CSV field count limit exceeded"))?;
        let storage = count
            .checked_mul(2 * std::mem::size_of::<String>())
            .and_then(|n| n.checked_add(line.len()))
            .filter(|&n| n <= self.text.limits().max_storage_bytes)
            .ok_or_else(|| invalid("CSV field storage limit exceeded"))?;
        let _ = storage;
        let mut fields = Vec::new();
        fields
            .try_reserve_exact(count)
            .map_err(|_| invalid("CSV field allocation failed"))?;
        for bytes in line.as_bytes().split(|&b| b == self.separator) {
            let bytes = if self.item_enclosed && delimiters != 0 {
                if bytes.is_empty() {
                    return Err(invalid("cannot remove enclosure from an empty CSV field"));
                }
                &bytes[1..bytes.len().saturating_sub(1).max(1)]
            } else {
                bytes
            };
            let value = std::str::from_utf8(bytes).map_err(|_| {
                Error::Unsupported("CSV byte delimiter/enclosure splits a UTF-8 character".into())
            })?;
            fields.push(value.to_owned());
        }
        Ok((delimiters != 0, fields))
    }
    /// Replace the output list only on success. Source's split-success boolean
    /// remains separate from checked row/enclosure/resource errors.
    pub fn get_row(&self, row: usize, fields: &mut Vec<String>) -> Result<bool> {
        let (split, result) = self.row(row)?;
        *fields = result;
        Ok(split)
    }
}
