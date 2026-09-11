// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::{MAX_BYTES, MAX_LABEL, Work, add, check_count, check_label, invalid, limit, mul};
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Native source parser extension point for named f64 masses and a byte reader.
/// Custom implementations own their parsing policy, budget and side effects.
/// Returned maps are validated again when consumed by an alphabet.
pub trait IMSAlphabetParser {
    fn parse(&mut self, reader: &mut dyn BufRead) -> Result<()>;
    fn elements(&self) -> &BTreeMap<String, f64>;
    /// Direct edits retain ordinary Rust semantics; later consumers revalidate them.
    fn elements_mut(&mut self) -> &mut BTreeMap<String, f64>;

    /// Opens a plain file before invoking the replaceable parser.
    fn load(&mut self, path: &Path) -> Result<()> {
        self.parse(&mut BufReader::new(std::fs::File::open(path)?))
    }
}

/// Atomic plain-text name/mass parser; first duplicate name wins.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IMSAlphabetTextParser {
    elements: BTreeMap<String, f64>,
}

impl IMSAlphabetTextParser {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn elements(&self) -> &BTreeMap<String, f64> {
        &self.elements
    }
    pub fn elements_mut(&mut self) -> &mut BTreeMap<String, f64> {
        &mut self.elements
    }

    pub fn parse(&mut self, reader: impl BufRead) -> Result<()> {
        let mut work = Work::default();
        measure_map_inner(&self.elements, &mut work, false)?;
        self.elements = parse_map(reader, &mut work)?;
        Ok(())
    }

    pub fn load(&mut self, path: impl AsRef<Path>) -> Result<()> {
        self.parse(BufReader::new(std::fs::File::open(path)?))
    }
}

impl IMSAlphabetParser for IMSAlphabetTextParser {
    fn parse(&mut self, reader: &mut dyn BufRead) -> Result<()> {
        self.parse(reader)
    }
    fn elements(&self) -> &BTreeMap<String, f64> {
        &self.elements
    }
    fn elements_mut(&mut self) -> &mut BTreeMap<String, f64> {
        &mut self.elements
    }
}

pub(super) fn parse_map(
    mut reader: impl BufRead,
    work: &mut Work,
) -> Result<BTreeMap<String, f64>> {
    let mut elements = BTreeMap::new();
    let mut line = Vec::new();
    let mut total_bytes = 0usize;
    let mut line_number = 0usize;
    while next_line(&mut reader, &mut line, &mut total_bytes, work)? {
        line_number += 1;
        work.consume(mul(line.len(), 2)?)?;
        if line.last() == Some(&b'\n') {
            line.pop();
        }
        let text = std::str::from_utf8(&line)
            .map_err(|_| syntax(line_number, "alphabet input must be UTF-8"))?;
        // The source comment precheck skips SPACE and TAB only, not all stream whitespace.
        let precheck = text.trim_start_matches([' ', '\t']);
        if precheck.is_empty() || precheck.starts_with('#') {
            continue;
        }
        let text = text.trim_start_matches(stream_space);
        let name_end = text.find(stream_space).unwrap_or(text.len());
        let name = &text[..name_end];
        if name.is_empty() {
            return Err(syntax(line_number, "missing alphabet name and mass"));
        }
        check_label(name)?;
        let mass = decimal_prefix(text[name_end..].trim_start_matches(stream_space))
            .map_err(|_| syntax(line_number, "invalid or nonfinite alphabet mass"))?;
        // A BTreeMap node compares several keys per level. This conservative
        // byte-comparison charge covers both lookup and insertion before either.
        let levels = (usize::BITS - elements.len().leading_zeros()) as usize + 1;
        work.consume(mul(mul(levels, 32)?, add(name.len(), 1)?)?)?;
        if !elements.contains_key(name) {
            check_count(add(elements.len(), 1)?)?;
            work.copy(name.len())?;
            // Logical node allowance, not an allocator-specific physical size.
            work.allocate(add(512, std::mem::size_of::<(String, f64)>())?)?;
            elements.insert(name.to_owned(), mass);
        }
    }
    Ok(elements)
}

pub(super) fn measure_map(elements: &BTreeMap<String, f64>, work: &mut Work) -> Result<()> {
    measure_map_inner(elements, work, true)
}
fn measure_map_inner(
    elements: &BTreeMap<String, f64>,
    work: &mut Work,
    check_values: bool,
) -> Result<()> {
    check_count(elements.len())?;
    work.consume(elements.len())?;
    let mut payload = mul(
        elements.len(),
        add(512, std::mem::size_of::<(String, f64)>())?,
    )?;
    for (name, mass) in elements {
        check_label(name)?;
        payload = add(payload, name.len())?;
        if payload > MAX_BYTES {
            return Err(limit());
        }
        work.consume(add(name.len(), 1)?)?;
        if check_values && !mass.is_finite() {
            return Err(invalid("custom alphabet parser returned nonfinite mass"));
        }
    }
    Ok(())
}

fn next_line(
    reader: &mut impl BufRead,
    line: &mut Vec<u8>,
    total: &mut usize,
    work: &mut Work,
) -> Result<bool> {
    line.clear();
    loop {
        let buffer = reader.fill_buf()?;
        if buffer.is_empty() {
            return Ok(!line.is_empty());
        }
        let allowance = (MAX_LABEL - line.len()).min(MAX_BYTES - *total);
        let mut take = 0;
        let mut newline = false;
        for &byte in buffer.iter().take(allowance + 1) {
            work.consume(1)?;
            take += 1;
            if take > allowance {
                return Err(invalid("IMS alphabet input byte or line limit exceeded"));
            }
            if byte == b'\n' {
                newline = true;
                break;
            }
        }
        work.copy(take)?;
        line.try_reserve_exact(take).map_err(|_| limit())?;
        line.extend_from_slice(&buffer[..take]);
        reader.consume(take);
        *total += take;
        if newline {
            return Ok(true);
        }
    }
}

// Portable native policy: decimal/scientific numeric prefix, followed by ignored
// extra text. Platform-dependent C++ hex/locale extraction is a custom-parser concern.
fn decimal_prefix(text: &str) -> Result<f64> {
    let bytes = text.as_bytes();
    let mut index = usize::from(matches!(bytes.first(), Some(b'+' | b'-')));
    if bytes
        .get(index..index + 2)
        .is_some_and(|b| b.eq_ignore_ascii_case(b"0x"))
    {
        return Err(invalid(
            "hexadecimal alphabet masses require a custom parser",
        ));
    }
    let before = index;
    while bytes.get(index).is_some_and(u8::is_ascii_digit) {
        index += 1;
    }
    let mut digits = index - before;
    if bytes.get(index) == Some(&b'.') {
        index += 1;
        let before = index;
        while bytes.get(index).is_some_and(u8::is_ascii_digit) {
            index += 1;
        }
        digits += index - before;
    }
    if digits == 0 {
        return Err(invalid("missing decimal mass"));
    }
    if matches!(bytes.get(index), Some(b'e' | b'E')) {
        index += 1;
        if matches!(bytes.get(index), Some(b'+' | b'-')) {
            index += 1;
        }
        let before = index;
        while bytes.get(index).is_some_and(u8::is_ascii_digit) {
            index += 1;
        }
        if index == before {
            return Err(invalid("incomplete mass exponent"));
        }
    }
    let mass: f64 = text[..index]
        .parse()
        .map_err(|_| invalid("invalid decimal mass"))?;
    if !mass.is_finite() {
        return Err(invalid("nonfinite decimal mass"));
    }
    Ok(mass)
}
fn syntax(line: usize, message: &str) -> Error {
    Error::Parse {
        line,
        message: message.into(),
    }
}

fn stream_space(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\r' | '\x0b' | '\x0c')
}
