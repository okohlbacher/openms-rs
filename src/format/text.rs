// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Source-compatible line buffers with CR/LF/CRLF reading and bounded UTF-8 storage.

use crate::{Error, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Resource limits may be lowered from these hard ceilings. Input limits include
/// skipped lines; line bytes exclude terminators, input bytes include them.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Limits {
    pub max_input_bytes: usize,
    pub max_line_bytes: usize,
    pub max_lines: usize,
    pub max_storage_bytes: usize,
    pub max_output_bytes: usize,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            max_input_bytes: 256 * 1024 * 1024,
            max_line_bytes: 16 * 1024 * 1024,
            max_lines: 1_000_000,
            max_storage_bytes: 256 * 1024 * 1024,
            max_output_bytes: 512 * 1024 * 1024,
        }
    }
}
impl Limits {
    pub(crate) fn validate(&self) -> Result<()> {
        let cap = Self::default();
        if self.max_input_bytes > cap.max_input_bytes
            || self.max_line_bytes > cap.max_line_bytes
            || self.max_lines > cap.max_lines
            || self.max_storage_bytes > cap.max_storage_bytes
            || self.max_output_bytes > cap.max_output_bytes
        {
            return Err(invalid("text limits exceed hard ceilings"));
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReadOptions {
    pub trim_lines: bool,
    /// Source semantics: only positive values stop after that many retained lines.
    /// Zero and every negative value read until EOF.
    pub first_n: i32,
    pub skip_empty_lines: bool,
    pub comment_symbol: String,
    pub limits: Limits,
}
impl Default for ReadOptions {
    fn default() -> Self {
        Self {
            trim_lines: false,
            first_n: -1,
            skip_empty_lines: false,
            comment_symbol: String::new(),
            limits: Limits::default(),
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct TextFile {
    buffer: Vec<String>,
    limits: Limits,
    storage_bytes: usize,
    dirty: bool,
}
impl PartialEq for TextFile {
    fn eq(&self, other: &Self) -> bool {
        self.buffer == other.buffer && self.limits == other.limits
    }
}
impl Eq for TextFile {}

impl TextFile {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn with_limits(limits: Limits) -> Result<Self> {
        limits.validate()?;
        Ok(Self {
            limits,
            ..Self::default()
        })
    }
    pub fn limits(&self) -> Limits {
        self.limits
    }
    pub fn lines(&self) -> &[String] {
        &self.buffer
    }
    pub fn iter(&self) -> std::slice::Iter<'_, String> {
        self.buffer.iter()
    }
    /// Like the source mutable iterator. Subsequent append/write operations
    /// recheck storage because callers can change line lengths through this view.
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, String> {
        self.dirty = true;
        self.buffer.iter_mut()
    }
    pub fn len(&self) -> usize {
        self.buffer.len()
    }
    pub fn is_empty(&self) -> bool {
        self.buffer.is_empty()
    }
    pub fn clear(&mut self) {
        self.buffer = Vec::new();
        self.storage_bytes = 0;
        self.dirty = false;
    }
    pub fn from_path(path: impl AsRef<Path>, options: &ReadOptions) -> Result<Self> {
        Self::from_reader(BufReader::new(File::open(path)?), options)
    }
    pub fn from_reader(mut reader: impl BufRead, options: &ReadOptions) -> Result<Self> {
        options.limits.validate()?;
        if options.comment_symbol.len() > options.limits.max_line_bytes {
            return Err(invalid("comment prefix exceeds text line limit"));
        }
        let mut result = Self::with_limits(options.limits)?;
        let mut budget = ReadBudget::new(options.limits);
        while let Some(bytes) = read_line(&mut reader, &mut budget)? {
            let mut line =
                String::from_utf8(bytes).map_err(|_| invalid("text input is not valid UTF-8"))?;
            if options.trim_lines {
                let start = line.len() - line.trim_start_matches([' ', '\t', '\n', '\r']).len();
                let end = start + trim(&line).len();
                line.truncate(end);
                line.drain(..start);
            }
            if (options.skip_empty_lines && line.is_empty())
                || (!options.comment_symbol.is_empty() && line.starts_with(&options.comment_symbol))
            {
                continue;
            }
            line.shrink_to_fit();
            result.push_owned(line)?;
            if options.first_n > 0 && result.len() == options.first_n as usize {
                break;
            }
        }
        Ok(result)
    }
    /// Replace the buffer only after the complete requested input succeeds.
    pub fn load(&mut self, path: impl AsRef<Path>, options: &ReadOptions) -> Result<()> {
        *self = Self::from_path(path, options)?;
        Ok(())
    }
    pub fn load_reader(&mut self, reader: impl BufRead, options: &ReadOptions) -> Result<()> {
        *self = Self::from_reader(reader, options)?;
        Ok(())
    }
    pub fn add_line(&mut self, line: impl AsRef<str>) -> Result<()> {
        let line = line.as_ref();
        self.check_append(line.len())?;
        self.push_owned(line.to_owned())
    }
    pub(crate) fn check_append(&mut self, line_bytes: usize) -> Result<()> {
        if self.dirty {
            self.storage_bytes = self.validate_storage()?;
            self.dirty = false;
        }
        if self.buffer.len() >= self.limits.max_lines || line_bytes > self.limits.max_line_bytes {
            return Err(invalid("text line count/length limit exceeded"));
        }
        let cost = line_cost(line_bytes)?;
        if cost
            > self
                .limits
                .max_storage_bytes
                .saturating_sub(self.storage_bytes)
        {
            return Err(invalid("text storage limit exceeded"));
        }
        Ok(())
    }
    pub(crate) fn push_owned(&mut self, line: String) -> Result<()> {
        self.check_append(line.len())?;
        let cost = line_cost(line.capacity())?;
        if cost
            > self
                .limits
                .max_storage_bytes
                .saturating_sub(self.storage_bytes)
        {
            return Err(invalid("text storage limit exceeded"));
        }
        // Bound the initial allocation as well as subsequent doubling: Vec's
        // default minimum is four String slots, exceeding one line's allowance.
        let reserved = if self.buffer.capacity() == 0 {
            self.buffer.try_reserve_exact(1)
        } else {
            self.buffer.try_reserve(1)
        };
        reserved.map_err(|_| invalid("text line allocation failed"))?;
        self.storage_bytes += cost;
        self.buffer.push(line);
        Ok(())
    }
    fn validate_storage(&self) -> Result<usize> {
        if self.buffer.len() > self.limits.max_lines {
            return Err(invalid("text line count limit exceeded"));
        }
        let mut bytes = 0usize;
        for line in &self.buffer {
            if line.len() > self.limits.max_line_bytes {
                return Err(invalid("text line length limit exceeded"));
            }
            bytes = bytes
                .checked_add(line_cost(line.capacity())?)
                .filter(|&n| n <= self.limits.max_storage_bytes)
                .ok_or_else(|| invalid("text storage limit exceeded"))?;
        }
        Ok(bytes)
    }
    fn output_size(&self) -> Result<()> {
        self.validate_storage()?;
        let mut remaining = self.limits.max_output_bytes;
        for line in &self.buffer {
            let (body, newline) = output_parts(line);
            let extra = if cfg!(windows) {
                body.bytes().filter(|&b| b == b'\n').count()
            } else {
                0
            };
            let size = body
                .len()
                .checked_add(extra)
                .and_then(|n| n.checked_add(if newline { NEWLINE.len() } else { 0 }))
                .ok_or_else(|| invalid("text output byte size overflow"))?;
            remaining = remaining
                .checked_sub(size)
                .ok_or_else(|| invalid("text output byte limit exceeded"))?;
        }
        Ok(())
    }
    pub fn write(&self, mut writer: impl Write) -> Result<()> {
        self.output_size()?;
        self.render(&mut writer)?;
        writer.flush()?;
        Ok(())
    }
    pub fn store(&self, path: impl AsRef<Path>) -> Result<()> {
        self.output_size()?; // precedes File::create and therefore truncation
        let mut writer = BufWriter::new(File::create(path)?);
        self.render(&mut writer)?;
        writer.flush()?;
        Ok(())
    }
    fn render(&self, writer: &mut impl Write) -> Result<()> {
        for line in &self.buffer {
            let (body, newline) = output_parts(line);
            write_platform_text(writer, body)?;
            if newline {
                writer.write_all(NEWLINE)?;
            }
        }
        Ok(())
    }
    /// Read one logical line, removing LF, CR, or CRLF. False means EOF without
    /// another line. On error the destination is unchanged; at EOF it is cleared.
    pub fn get_line(reader: &mut impl BufRead, destination: &mut String) -> Result<bool> {
        Self::get_line_with_limits(reader, destination, &Limits::default())
    }
    /// These limits apply to this single call, rather than preceding calls.
    pub fn get_line_with_limits(
        reader: &mut impl BufRead,
        destination: &mut String,
        limits: &Limits,
    ) -> Result<bool> {
        limits.validate()?;
        let mut budget = ReadBudget::new(*limits);
        let value = read_line(reader, &mut budget)?;
        match value {
            Some(bytes) => {
                let line = String::from_utf8(bytes)
                    .map_err(|_| invalid("text input is not valid UTF-8"))?;
                if line_cost(line.capacity())? > limits.max_storage_bytes {
                    return Err(invalid("text storage limit exceeded"));
                }
                *destination = line;
                Ok(true)
            }
            None => {
                destination.clear();
                Ok(false)
            }
        }
    }
}

impl<'a> IntoIterator for &'a TextFile {
    type Item = &'a String;
    type IntoIter = std::slice::Iter<'a, String>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}
impl<'a> IntoIterator for &'a mut TextFile {
    type Item = &'a mut String;
    type IntoIter = std::slice::IterMut<'a, String>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter_mut()
    }
}

pub(crate) fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
pub(crate) fn trim(line: &str) -> &str {
    line.trim_matches([' ', '\t', '\n', '\r'])
}
fn line_cost(bytes: usize) -> Result<usize> {
    // Includes at most doubled Vec capacity plus the UTF-8 payload.
    bytes
        .checked_mul(2)
        .and_then(|n| n.checked_add(2 * std::mem::size_of::<String>()))
        .ok_or_else(|| invalid("text storage byte size overflow"))
}
struct ReadBudget {
    limits: Limits,
    bytes: usize,
    lines: usize,
}
impl ReadBudget {
    fn new(limits: Limits) -> Self {
        Self {
            limits,
            bytes: 0,
            lines: 0,
        }
    }
    fn consume(&mut self, bytes: usize) -> Result<()> {
        self.bytes = self
            .bytes
            .checked_add(bytes)
            .filter(|&n| n <= self.limits.max_input_bytes)
            .ok_or_else(|| invalid("text input byte limit exceeded"))?;
        Ok(())
    }
}
fn read_line(reader: &mut impl BufRead, budget: &mut ReadBudget) -> Result<Option<Vec<u8>>> {
    if reader.fill_buf()?.is_empty() {
        return Ok(None);
    }
    if budget.lines >= budget.limits.max_lines {
        return Err(invalid("text input line count limit exceeded"));
    }
    budget.lines += 1;
    let mut line = Vec::new();
    loop {
        let buffer = reader.fill_buf()?;
        if buffer.is_empty() {
            return Ok(Some(line));
        }
        let scan = buffer
            .len()
            .min(budget.limits.max_line_bytes.saturating_sub(line.len()) + 1)
            .min(budget.limits.max_input_bytes.saturating_sub(budget.bytes) + 1);
        let window = &buffer[..scan];
        let end = window.iter().position(|&b| b == b'\r' || b == b'\n');
        let content = end.unwrap_or(scan);
        if content > budget.limits.max_line_bytes.saturating_sub(line.len()) {
            return Err(invalid("text input line length limit exceeded"));
        }
        let consumed = content + usize::from(end.is_some());
        budget.consume(consumed)?;
        line.try_reserve(content)
            .map_err(|_| invalid("text input allocation failed"))?;
        line.extend_from_slice(&window[..content]);
        let cr = end.is_some_and(|i| window[i] == b'\r');
        reader.consume(consumed);
        if cr && reader.fill_buf()?.first() == Some(&b'\n') {
            budget.consume(1)?;
            reader.consume(1);
        }
        if end.is_some() {
            return Ok(Some(line));
        }
    }
}
#[cfg(windows)]
const NEWLINE: &[u8] = b"\r\n";
#[cfg(not(windows))]
const NEWLINE: &[u8] = b"\n";
fn output_parts(line: &str) -> (&str, bool) {
    if let Some(body) = line.strip_suffix("\r\n") {
        (body, true)
    } else {
        (line, !line.ends_with('\n'))
    }
}
fn write_platform_text(writer: &mut impl Write, text: &str) -> std::io::Result<()> {
    if cfg!(windows) {
        for part in text.split_inclusive('\n') {
            if let Some(body) = part.strip_suffix('\n') {
                writer.write_all(body.as_bytes())?;
                writer.write_all(NEWLINE)?;
            } else {
                writer.write_all(part.as_bytes())?;
            }
        }
        Ok(())
    } else {
        writer.write_all(text.as_bytes())
    }
}
