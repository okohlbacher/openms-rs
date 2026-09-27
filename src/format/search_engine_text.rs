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
//! characters first. See `docs/SEARCH_ENGINE_FILES_SUPPORT.md`.

use crate::chemistry::AASequence;
use crate::{Error, Result};

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
