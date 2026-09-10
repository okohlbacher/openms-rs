// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Charge suffixes and ordinals in fragment-ion names.
//!
//! Parsing and precedence follow OpenMS4-core revision `7c029e8`,
//! [`IonNaming.h`](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8/src/openms/include/OpenMS/CHEMISTRY/IonNaming.h).
//! The parsers borrow bytes, allocate nothing, and recognize ASCII letters,
//! digits and signs. They impose no input-size cap. Malformed or unrepresentable
//! numeric tokens return zero; malformed caret syntax falls through to the
//! source's trailing-sign and trailing-number checks. Only [`with_charge`]
//! copies arbitrary input, with a checked one-MiB output bound.

use crate::{Error, Result};

/// Largest charge magnitude represented as repeated signs, from the source.
pub const MAX_REPEATED_SIGNS: i32 = 8;
/// Maximum UTF-8 output bytes for [`with_charge`], including every text line.
pub const MAX_ION_NAME_BYTES: usize = 1024 * 1024;

/// Spell a charge as repeated signs through magnitude eight, then sign and
/// magnitude. Zero is empty; even `i32::MIN` needs at most eleven output bytes.
pub fn charge_suffix(charge: i32) -> String {
    if charge == 0 {
        return String::new();
    }
    let sign = if charge < 0 { "-" } else { "+" };
    let magnitude = charge.unsigned_abs();
    if magnitude <= MAX_REPEATED_SIGNS as u32 {
        return sign.repeat(magnitude as usize);
    }
    format!("{sign}{magnitude}")
}

/// Read a charge from the first line, ending at the first CR or LF.
///
/// The last caret token has priority: an optional sign and at most ten digits
/// must end at the line boundary, '/' or '*'. A valid zero token returns zero
/// immediately. Otherwise try a trailing same-sign run, then sign plus digits.
/// The numeric suffix is suppressed if '/' or '*' occurs anywhere in the line;
/// the sign-run fallback is not. Mixed signs count only the final repeated sign.
/// Unrepresentable runs return zero instead of source integer narrowing overflow.
pub fn charge_from_name(ion_name: &str) -> i32 {
    let name = first_line(ion_name);
    let Some(&last) = name.last() else {
        return 0;
    };
    if let Some(caret) = name.iter().rposition(|&byte| byte == b'^') {
        let mut begin = caret + 1;
        let negative = name.get(begin) == Some(&b'-');
        if matches!(name.get(begin), Some(b'+' | b'-')) {
            begin += 1;
        }
        let mut end = begin;
        while end < name.len() && name[end].is_ascii_digit() {
            end += 1;
        }
        if end == name.len() || matches!(name[end], b'/' | b'*') {
            if let Some(charge) = parse_magnitude(&name[begin..end], negative) {
                return charge;
            }
        }
    }
    if matches!(last, b'+' | b'-') {
        let count = name.iter().rev().take_while(|&&byte| byte == last).count();
        return i64::try_from(count)
            .ok()
            .and_then(|count| signed_charge(count, last == b'-'))
            .unwrap_or(0);
    }
    if last.is_ascii_digit() && !name.iter().any(|&byte| matches!(byte, b'/' | b'*')) {
        let mut begin = name.len();
        while begin > 0 && name[begin - 1].is_ascii_digit() {
            begin -= 1;
        }
        if begin > 0 && matches!(name[begin - 1], b'+' | b'-') {
            return parse_magnitude(&name[begin..], name[begin - 1] == b'-').unwrap_or(0);
        }
    }
    0
}

/// Add a missing nonzero charge at the end of the first line, preserving all
/// line endings and free text. An existing nonzero charge wins even if it
/// disagrees with `charge`; a zero charge adds nothing.
///
/// The complete returned string must fit [`MAX_ION_NAME_BYTES`], including when
/// no suffix is added. Oversized input/output errors before copying. Embedded
/// Unicode and NUL bytes are preserved; charge recognition remains ASCII.
pub fn with_charge(ion_name: &str, charge: i32) -> Result<String> {
    if ion_name.len() > MAX_ION_NAME_BYTES {
        return Err(output_limit());
    }
    if charge == 0 || charge_from_name(ion_name) != 0 {
        return Ok(ion_name.to_owned());
    }
    let suffix = charge_suffix(charge);
    let output_len = ion_name
        .len()
        .checked_add(suffix.len())
        .ok_or_else(output_limit)?;
    if output_len > MAX_ION_NAME_BYTES {
        return Err(output_limit());
    }
    let end = first_line(ion_name).len();
    // CR/LF are one-byte ASCII delimiters, hence valid UTF-8 slicing boundaries.
    let mut output = String::with_capacity(output_len);
    output.push_str(&ion_name[..end]);
    output.push_str(&suffix);
    output.push_str(&ion_name[end..]);
    Ok(output)
}

/// Read at most nine ASCII digits immediately after one ASCII letter.
///
/// Any ASCII letter, including lowercase, is accepted. Everything after the
/// digit run is ignored. Missing/misplaced ordinals or runs longer than nine
/// digits return zero, including overlong leading-zero runs. Callers must still
/// check the ordinal against their peptide length before indexing a sequence.
pub fn ordinal_from_name(ion_name: &str) -> u32 {
    let bytes = ion_name.as_bytes();
    if !bytes.first().is_some_and(u8::is_ascii_alphabetic) {
        return 0;
    }
    let mut value = 0_u32;
    for (index, &digit) in bytes[1..]
        .iter()
        .take_while(|&&byte| byte.is_ascii_digit())
        .enumerate()
    {
        if index == 9 {
            return 0;
        }
        value = value * 10 + u32::from(digit - b'0');
    }
    value
}

fn first_line(name: &str) -> &[u8] {
    let bytes = name.as_bytes();
    let end = bytes
        .iter()
        .position(|&byte| matches!(byte, b'\r' | b'\n'))
        .unwrap_or(bytes.len());
    &bytes[..end]
}

// Callers have already identified an ASCII digit run. Ten digits fit i64 even
// when their signed result cannot fit i32; no lossy narrowing is permitted.
fn parse_magnitude(digits: &[u8], negative: bool) -> Option<i32> {
    if digits.is_empty() || digits.len() > 10 {
        return None;
    }
    let magnitude = digits
        .iter()
        .fold(0_i64, |value, &digit| value * 10 + i64::from(digit - b'0'));
    signed_charge(magnitude, negative)
}

fn signed_charge(magnitude: i64, negative: bool) -> Option<i32> {
    i32::try_from(if negative { -magnitude } else { magnitude }).ok()
}

fn output_limit() -> Error {
    Error::InvalidValue("ion-name output exceeds 1 MiB".into())
}
