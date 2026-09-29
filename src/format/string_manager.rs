// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! UTF-16 string helpers of the XML handlers (`FORMAT/HANDLERS/StringManager.h`).
//!
//! The source converts between OpenMS's narrow strings and the UTF-16 code
//! units (`XMLCh`, `char16_t`) that Xerces-C hands its SAX handlers, and
//! adds SIMD fast paths for the pure-ASCII case that dominates Base64 peak
//! data. The Rust XML layer reads UTF-8 directly, so these helpers are only
//! needed where UTF-16 data meets the port; they operate on `&[u16]` slices.
//!
//! # Not ported, and why
//!
//! - `shared_xerces_ptr<T>` and `unique_xerces_ptr<T>` are RAII guards that
//!   release Xerces-allocated memory through `XMLString::release` or the
//!   node's `release()`. Owned Rust values (`Vec<u16>`, `String`) free
//!   themselves, so no guard type exists; [`StringManager::convert_ptr`]
//!   returns the owned buffer those guards would hold.
//! - `CONST_XMLCH(s)`, a UTF-16 string literal, is `"s".encode_utf16()`.
//!
//! # Terminators
//!
//! The source takes NUL-terminated `const char16_t*` strings. Here a slice
//! carries its length, and every function that the source bounds by the
//! terminator ([`StringManager::convert_utf16`], [`StringManager::str_length`],
//! [`StringManager::parse_int`]) stops at the first `0` unit when there is
//! one, so a terminated buffer gives the source's result. Where the source
//! would read past the end of an unterminated buffer, the slice end is used.
//!
//! The SIMD code paths (`strLength`'s aligned 16-byte loads, `isASCII`'s
//! 8-unit masks, `compress64_`'s byte shuffle) are replaced by scalar loops
//! with the same results; they needed SIMDe and an AddressSanitizer
//! exemption for reading past the terminator, which safe Rust cannot do and
//! does not need.

use crate::{Error, Result};

/// Conversions between UTF-8 strings and UTF-16 code units for XML parsing.
///
/// Every function is associated; the source's default constructor and
/// destructor do nothing, so the type is an empty `Default` value.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct StringManager;

impl StringManager {
    /// The empty helper, as the source's default constructor.
    pub fn new() -> Self {
        Self
    }

    /// Transcode `text` to UTF-16 code units.
    ///
    /// The source transcodes a narrow C string with Xerces'
    /// `XMLString::transcode`, which uses the process's local code page; for
    /// ASCII that is the identity on code units. A Rust `&str` is UTF-8, so
    /// every character is transcoded, and supplementary characters become
    /// surrogate pairs. The source's `convert(const char*)` and
    /// `convert(const std::string&)` are both this function.
    pub fn convert(text: &str) -> Vec<u16> {
        text.encode_utf16().collect()
    }

    /// [`StringManager::convert`] followed by a terminating `0`, the buffer
    /// the source's `convertPtr` returns inside a `unique_xerces_ptr`.
    pub fn convert_ptr(text: &str) -> Vec<u16> {
        let mut units = Self::convert(text);
        units.push(0);
        units
    }

    /// Transcode UTF-16 code units, up to the first `0` unit, to a string.
    ///
    /// The source's `convert(const char16_t*)`: when every unit is ASCII the
    /// low bytes are copied ([`StringManager::append_ascii`]); otherwise the
    /// text is transcoded to UTF-8.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the units contain an unpaired surrogate,
    /// which has no UTF-8 encoding. Xerces' UTF-8 transcoder refuses such
    /// input too.
    pub fn convert_utf16(units: &[u16]) -> Result<String> {
        let units = &units[..Self::str_length(units)];
        let mut text = String::with_capacity(units.len());
        if Self::is_ascii(units) {
            Self::append_ascii(units, &mut text)?;
            return Ok(text);
        }
        String::from_utf16(units)
            .map_err(|_| Error::InvalidValue("UTF-16 string contains an unpaired surrogate".into()))
    }

    /// The number of code units before the first `0` unit, or the slice length
    /// when there is none.
    ///
    /// The source's `strLength` returns `0` for a null pointer; an empty slice
    /// does the same.
    pub fn str_length(units: &[u16]) -> usize {
        units
            .iter()
            .position(|&unit| unit == 0)
            .unwrap_or(units.len())
    }

    /// Parse a decimal integer from UTF-16 code units (up to the first `0`),
    /// as Xerces' `XMLString::parseInt` does, refusing values outside `i32`.
    ///
    /// The rules are Xerces': XML whitespace (space, tab, line feed, carriage
    /// return) is trimmed from both ends; the rest must be exactly what C
    /// `strtol` in base 10 consumes, which also skips a leading vertical tab
    /// or form feed and takes one optional `+` or `-` sign followed by at
    /// least one ASCII digit.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the trimmed text is empty or not entirely
    /// a number (Xerces' `NumberFormatException`), and [`Error::InvalidRange`]
    /// when the value does not fit in `i32`. Xerces raises its range error
    /// only outside the C `long` range and then truncates the `long` to `int`;
    /// [`StringManager::parse_int_source`] reproduces that.
    pub fn parse_int(units: &[u16]) -> Result<i32> {
        let value = parse_long(units)?;
        i32::try_from(value)
            .map_err(|_| Error::InvalidRange(format!("integer {value} does not fit in 32 bits")))
    }

    /// [`StringManager::parse_int`] with the source's handling of large
    /// values: anything in the 64-bit `long` range is accepted and truncated
    /// to its low 32 bits, as the source's `(int)` cast does on LP64
    /// platforms, so `"4294967296"` parses as `0`.
    ///
    /// # Errors
    ///
    /// As [`StringManager::parse_int`], with [`Error::InvalidRange`] only
    /// outside the `i64` range.
    pub fn parse_int_source(units: &[u16]) -> Result<i32> {
        Ok(parse_long(units)? as i32)
    }

    /// Whether every unit can be encoded in ASCII, that is, whether each
    /// unit is below `0x80` (no bit of `0xFF80` set). An empty slice is ASCII.
    pub fn is_ascii(units: &[u16]) -> bool {
        units.iter().all(|&unit| unit & 0xFF80 == 0)
    }

    /// Append the characters of `units`, which must all be ASCII, to
    /// `result`.
    ///
    /// The source keeps the low byte of each unit, eight units at a time,
    /// and documents the ASCII-only input as an assumption.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a unit is not ASCII; the source would
    /// append its low byte, producing a different character or invalid UTF-8.
    /// `result` is unchanged on error.
    pub fn append_ascii(units: &[u16], result: &mut String) -> Result<()> {
        if !Self::is_ascii(units) {
            return Err(Error::InvalidValue(
                "appendASCII input contains a non-ASCII code unit".into(),
            ));
        }
        let mut bytes = Vec::with_capacity(units.len());
        let mut blocks = units.chunks_exact(8);
        for block in blocks.by_ref() {
            bytes.extend_from_slice(&compress64(block));
        }
        bytes.extend(blocks.remainder().iter().map(|&unit| (unit & 0xFF) as u8));
        // Every byte is below 0x80, so this is valid UTF-8.
        result.push_str(&String::from_utf8_lossy(&bytes));
        Ok(())
    }
}

/// The source's protected `compress64_`: the low bytes of eight UTF-16 units.
fn compress64(units: &[u16]) -> [u8; 8] {
    let mut out = [0u8; 8];
    for (byte, unit) in out.iter_mut().zip(units) {
        *byte = (unit & 0xFF) as u8;
    }
    out
}

/// `XMLString::parseInt` up to its final `(int)` cast.
fn parse_long(units: &[u16]) -> Result<i64> {
    let units = &units[..StringManager::str_length(units)];
    let xml_space = |unit: &u16| matches!(*unit, 0x20 | 0x09 | 0x0A | 0x0D);
    let start = units
        .iter()
        .position(|unit| !xml_space(unit))
        .unwrap_or(units.len());
    let end = units
        .iter()
        .rposition(|unit| !xml_space(unit))
        .map_or(start, |i| i + 1);
    let trimmed = &units[start..end];
    let invalid = || {
        Error::InvalidValue(format!(
            "not an integer: {:?}",
            String::from_utf16_lossy(trimmed)
        ))
    };
    if trimmed.is_empty() {
        return Err(invalid());
    }
    // strtol: leading isspace (only \v and \f remain after the trim), an
    // optional sign, then digits.
    let mut i = 0;
    while i < trimmed.len() && matches!(trimmed[i], 0x0B | 0x0C) {
        i += 1;
    }
    let negative = match trimmed.get(i) {
        Some(0x2D) => {
            i += 1;
            true
        }
        Some(0x2B) => {
            i += 1;
            false
        }
        _ => false,
    };
    let digits = &trimmed[i..];
    if digits.is_empty() || !digits.iter().all(|&unit| (0x30..=0x39).contains(&unit)) {
        return Err(invalid());
    }
    let mut value: i64 = 0;
    for &unit in digits {
        let digit = i64::from(unit - 0x30);
        value = value
            .checked_mul(10)
            .and_then(|v| {
                if negative {
                    v.checked_sub(digit)
                } else {
                    v.checked_add(digit)
                }
            })
            .ok_or_else(|| {
                Error::InvalidRange(format!(
                    "integer {:?} is out of range",
                    String::from_utf16_lossy(trimmed)
                ))
            })?;
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn u(text: &str) -> Vec<u16> {
        text.encode_utf16().collect()
    }

    /// `XMLHandler_test.cpp` (the StringManager class test) at core
    /// `bc9cc12`, section `compress64`.
    #[test]
    fn compress64_class_test() {
        let eight_block = [
            0x0048, 0x0065, 0x006C, 0x006C, 0x006F, 0x002C, 0x0057, 0x006F,
        ];
        let eight_block_negative = [
            0x0148, 0x0165, 0x016C, 0x016C, 0x016F, 0x012C, 0x0157, 0x016F,
        ];
        let eight_block_mixed = [
            0x0042, 0x0045, 0x004C, 0x0041, 0x0142, 0x0145, 0x014C, 0x0141,
        ];
        let eight_block_kadabra = [
            0x004B, 0x0041, 0x0044, 0x0041, 0x0042, 0x0052, 0x0041, 0x0021,
        ];
        assert_eq!(&compress64(&eight_block), b"Hello,Wo");
        assert_eq!(&compress64(&eight_block_negative), b"Hello,Wo");
        assert_eq!(
            compress64(&eight_block_mixed),
            [0x42, 0x45, 0x4C, 0x41, 0x42, 0x45, 0x4C, 0x41]
        );
        let mut o4 = b"ABRA".to_vec();
        o4.extend_from_slice(&compress64(&eight_block_kadabra));
        assert_eq!(o4, b"ABRAKADABRA!");
        assert_eq!(compress64(&u("ABCDEFGH")), *b"ABCDEFGH");
    }

    #[test]
    fn parse_long_follows_strtol_after_the_xml_trim() {
        assert_eq!(parse_long(&u(" \t42\r\n")).unwrap(), 42);
        assert_eq!(parse_long(&u("\u{b}-7")).unwrap(), -7);
        assert!(parse_long(&u("7\u{b}")).is_err());
        assert!(parse_long(&u("+ 5")).is_err());
        assert!(parse_long(&u("--5")).is_err());
        assert!(parse_long(&u("0x10")).is_err());
        assert_eq!(parse_long(&u("-9223372036854775808")).unwrap(), i64::MIN);
        assert!(matches!(
            parse_long(&u("9223372036854775808")),
            Err(Error::InvalidRange(_))
        ));
    }
}
