// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Helpers shared by the search-engine file adapters.
//!
//! The adapters for OMSSA, Sequest, Inspect, PepNovo and X! Tandem files read
//! peptide sequences with `AASequence::fromString(s)`, whose default is the
//! *permissive* parse (`CHEMISTRY/AASequence.cpp:1390-1397` at the core pin):
//! outside a modification bracket, a stop-codon character `*`, `#` or `+`
//! becomes the unknown residue `X` and a space is skipped. The crate's
//! [`AASequence::parse`] is the strict parse, so this module rewrites those
//! characters first.
//!
//! The Sequest, Inspect and PepNovo adapters read their files with
//! `std::getline` and index into the lines by byte, so this module also holds
//! a bounded whole-file reader that splits exactly where `getline` does and
//! the `std::string` search helpers they need with the source's `npos`
//! arithmetic. See `docs/SEARCH_ENGINE_FILES_SUPPORT.md`.

use crate::chemistry::AASequence;
use crate::{Error, Result};
use std::io::Read;

/// The largest text file the search-engine text adapters read, in bytes.
pub const MAX_TEXT_BYTES: usize = 256 * 1024 * 1024;

/// `std::string::npos`.
pub(crate) const NPOS: usize = usize::MAX;

/// Read a whole text file of at most `max_bytes` bytes.
///
/// # Errors
///
/// [`Error::Io`] when reading fails, [`Error::InvalidValue`] when the input
/// is longer than `max_bytes`, [`Error::Parse`] when it is not UTF-8.
pub(crate) fn read_text(input: impl Read, max_bytes: usize) -> Result<String> {
    let count = u64::try_from(max_bytes)
        .ok()
        .and_then(|value| value.checked_add(1))
        .ok_or_else(|| Error::InvalidValue("text byte limit overflows".into()))?;
    let mut bytes = Vec::new();
    input.take(count).read_to_end(&mut bytes)?;
    if bytes.len() > max_bytes {
        return Err(Error::InvalidValue("text byte limit exceeded".into()));
    }
    String::from_utf8(bytes).map_err(|_| Error::Parse {
        line: 0,
        message: "text input is not valid UTF-8".into(),
    })
}

/// The lines `std::getline` returns: split at LF only, so a CR stays in the
/// line, and a final LF does not start another line.
pub(crate) fn getline_lines(text: &str) -> Vec<&str> {
    let mut lines: Vec<&str> = text.split('\n').collect();
    if text.is_empty() || text.ends_with('\n') {
        lines.pop();
    }
    lines
}

/// The adapters' line clean-up: drop one final byte below 33 (a CR, for
/// instance), then trim space, tab, LF and CR.
///
/// The source tests `line[line.length() - 1] < 33` on a `char`, which is
/// signed on x86, so there it also drops the final byte of a multi-byte UTF-8
/// character. Only an ASCII byte is dropped here.
pub(crate) fn clean_line(line: &str) -> &str {
    let line = match line.as_bytes().last() {
        Some(&byte) if byte < 33 => &line[..line.len() - 1],
        _ => line,
    };
    line.trim_matches([' ', '\t', '\n', '\r'])
}

/// `std::string::find(c, from)`: the byte position of `c` at or after
/// `from`, or [`NPOS`].
pub(crate) fn find_byte(text: &str, byte: u8, from: usize) -> usize {
    text.as_bytes()
        .get(from..)
        .and_then(|rest| rest.iter().position(|&b| b == byte))
        .map_or(NPOS, |position| position + from)
}

/// `std::string::find_last_of(c)`: the last byte position of `c`, or
/// [`NPOS`].
pub(crate) fn rfind_byte(text: &str, byte: u8) -> usize {
    text.as_bytes()
        .iter()
        .rposition(|&b| b == byte)
        .unwrap_or(NPOS)
}

/// `StringUtils::substr(s, pos, n)`: `pos` clamped to the length, then at
/// most `n` bytes.
///
/// # Errors
///
/// [`Error::Parse`] when the cut splits a multi-byte UTF-8 character, which
/// the source would pass on as a broken byte sequence.
pub(crate) fn substr(text: &str, pos: usize, n: usize) -> Result<&str> {
    let begin = pos.min(text.len());
    let end = begin.saturating_add(n).min(text.len());
    text.get(begin..end).ok_or_else(|| Error::Parse {
        line: 0,
        message: format!("cutting '{text}' at bytes {begin}..{end} splits a character"),
    })
}

/// Source `AASequence::fromString(text, permissive = true)`.
///
/// Outside `(...)` and `[...]` annotations, `*`, `#` and `+` become `X` and
/// spaces are dropped; everything else, annotations included, goes to the
/// strict [`AASequence::parse`] unchanged.
///
/// # Errors
///
/// Whatever [`AASequence::parse`] reports for the rewritten text (the source's
/// `Exception::ParseError` for an unknown residue or modification).
pub fn parse_permissive_sequence(text: &str) -> Result<AASequence> {
    let mut rewritten = String::with_capacity(text.len());
    let mut depth = 0usize;
    for c in text.chars() {
        match c {
            '(' | '[' => {
                depth = depth.saturating_add(1);
                rewritten.push(c);
            }
            ')' | ']' => {
                depth = depth.saturating_sub(1);
                rewritten.push(c);
            }
            '*' | '#' | '+' if depth == 0 => rewritten.push('X'),
            ' ' if depth == 0 => {}
            _ => rewritten.push(c),
        }
    }
    AASequence::parse(&rewritten)
}

/// Attach a physical line number to a conversion or sequence error that
/// carries none.
pub(crate) fn at_line(error: Error, line: usize) -> Error {
    match error {
        Error::Parse { .. } => error,
        other => Error::Parse {
            line,
            message: other.to_string(),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn getline_splits_at_line_feeds_only() {
        assert!(getline_lines("").is_empty());
        assert_eq!(getline_lines("a\r\nb"), ["a\r", "b"]);
        assert_eq!(getline_lines("a\n\n"), ["a", ""]);
        assert_eq!(clean_line(" a b\x01"), "a b");
        assert_eq!(clean_line("a\t\r"), "a");
        assert_eq!(find_byte("a|b|c", b'|', 2), 3);
        assert_eq!(find_byte("a|b", b'|', 9), NPOS);
        assert_eq!(rfind_byte("a|b|c", b'|'), 3);
        assert_eq!(substr("abcdef", 2, NPOS).unwrap(), "cdef");
        assert_eq!(substr("abc", 7, 2).unwrap(), "");
        assert!(substr("aé", 0, 2).is_err());
    }

    #[test]
    fn permissive_parse_maps_stop_codons_and_skips_spaces() {
        assert_eq!(
            parse_permissive_sequence("PEP TIDE*").unwrap(),
            AASequence::parse("PEPTIDEX").unwrap()
        );
        assert_eq!(
            parse_permissive_sequence("PEPM(Oxidation)K").unwrap(),
            AASequence::parse("PEPM(Oxidation)K").unwrap()
        );
        assert!(parse_permissive_sequence("PEP1").is_err());
    }
}
