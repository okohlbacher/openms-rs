// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! File adapter for OMSSA comma-separated result files.
//!
//! Ports `FORMAT/OMSSACSVFile.h`/`.cpp`. The file holds the results of the
//! OMSSA search engine, one peptide hit per line after a header that begins
//! with `Spectrum`:
//!
//! ```text
//! Spectrum number, Filename/id, Peptide, E-value, Mass, gi, Accession, Start,
//! Stop, Defline, Mods, Charge, Theo Mass, P-value
//! ```
//!
//! See `docs/SEARCH_ENGINE_FILES_SUPPORT.md` for the API mapping.

use crate::concept::string_conversion::{split, to_double, to_int32, trim};
use crate::format::search_engine_text::{at_line, parse_permissive_sequence};
use crate::format::text::{ReadOptions, TextFile};
use crate::identification::{PeptideHit, PeptideIdentification};
use crate::{Error, Result};
use std::io::BufRead;
use std::path::Path;

/// The number of comma-separated columns of an OMSSA CSV line whose `Defline`
/// holds no quoted commas.
pub const COLUMN_COUNT: usize = 14;

/// File adapter for OMSSA CSV files.
///
/// Each data line becomes one [`PeptideHit`] with the sequence of column 3,
/// the charge of column 12 and, as score, the *P-value* of column 14 (the
/// source's choice; the E-value of column 4 is not read). Consecutive lines
/// with the same spectrum number (column 1) share one
/// [`PeptideIdentification`] with score type `OMSSA`; every change of that
/// number starts a new one, so a spectrum that reappears later gets a second
/// identification. The source starts from the sentinel spectrum number -1, so
/// a first data line numbered -1 adds its hit to `id_data.back()`: the
/// caller's last identification, or undefined behaviour on an empty list. The
/// port always starts a new identification for the first line.
///
/// A `Defline` that contains commas must be quoted. The source counts the
/// commas between the first `"` and the next one (or the end of the line) and
/// then accepts a line with 14 fields or 14 plus that count, reading the
/// charge and P-value that many fields further right. The fields are not
/// dequoted.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct OMSSACSVFile;

impl OMSSACSVFile {
    /// An OMSSA CSV file adapter; it holds no state.
    pub fn new() -> Self {
        Self
    }

    /// Load an OMSSA CSV file into peptide identifications.
    ///
    /// The source's `load(filename, protein_identification, id_data)` never
    /// touches `protein_identification` and appends to `id_data` without
    /// clearing it; the port therefore takes no protein identification and
    /// returns the new identifications for the caller to append.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be opened or read (source
    /// `Exception::FileNotFound`, `FileNotReadable` or `IOException`); see
    /// [`OMSSACSVFile::load_reader`] for the parse errors.
    pub fn load(&self, path: impl AsRef<Path>) -> Result<Vec<PeptideIdentification>> {
        let input = TextFile::from_path(path, &ReadOptions::default())?;
        parse(&input)
    }

    /// Load OMSSA CSV content from a reader; see [`OMSSACSVFile::load`].
    ///
    /// Data lines are trimmed of space, tab, newline and carriage return
    /// before they are split, and each field is trimmed again before it is
    /// converted.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] when the first line does not begin with `Spectrum`
    /// (an empty input included; the source message "first line does not
    /// contain a description"), when a data line has neither 14 fields nor 14
    /// plus its quoted-comma count (the source message "number of columns
    /// should be 14 in line N", where N counts data lines from one), when a
    /// numeric field does not convert (the source's `ConversionError`, which it
    /// does not wrap), when the sequence does not parse, or when the P-value
    /// is not finite. The error's `line` is always the physical line.
    /// Lines are split at LF, CR or CRLF, as [`TextFile`] does; the source's
    /// `std::getline` splits at LF only, so a CR-only file is one line to the
    /// source, which then loads nothing. On error nothing is returned, whereas
    /// the source has already appended the identifications before the failing
    /// line.
    pub fn load_reader(&self, reader: impl BufRead) -> Result<Vec<PeptideIdentification>> {
        let input = TextFile::from_reader(reader, &ReadOptions::default())?;
        parse(&input)
    }
}

fn parse(input: &TextFile) -> Result<Vec<PeptideIdentification>> {
    let lines = input.lines();
    if !lines
        .first()
        .is_some_and(|line| line.starts_with("Spectrum"))
    {
        return Err(Error::Parse {
            line: 1,
            message: "first line does not contain a description".into(),
        });
    }
    let mut identifications: Vec<PeptideIdentification> = Vec::new();
    let mut current_spectrum: Option<i32> = None;
    for (data_index, raw) in lines.iter().enumerate().skip(1) {
        let line_number = data_index + 1;
        let line = trim(raw);
        let offset = quoted_comma_count(line);
        let fields = split(line, ',');
        if fields.len() != COLUMN_COUNT && fields.len() != COLUMN_COUNT + offset {
            return Err(Error::Parse {
                line: line_number,
                message: format!("number of columns should be 14 in line {data_index}"),
            });
        }
        let parse_error = |message: String| Error::Parse {
            line: line_number,
            message,
        };
        let sequence =
            parse_permissive_sequence(trim(fields[2])).map_err(|e| at_line(e, line_number))?;
        let score = to_double(trim(fields[13 + offset])).map_err(parse_error)?;
        let charge = to_int32(trim(fields[11 + offset])).map_err(parse_error)?;
        let hit =
            PeptideHit::new(score, 0, charge, sequence).map_err(|e| at_line(e, line_number))?;
        let spectrum = to_int32(trim(fields[0])).map_err(parse_error)?;
        if current_spectrum != Some(spectrum) {
            identifications.push(PeptideIdentification {
                score_type: "OMSSA".into(),
                ..PeptideIdentification::default()
            });
            current_spectrum = Some(spectrum);
        }
        if let Some(identification) = identifications.last_mut() {
            identification.hits.push(hit);
        }
    }
    Ok(identifications)
}

/// The source's count of commas between the first `"` of the line and the
/// next `"` (or the end of the line, the bounded loop of the regression fix).
fn quoted_comma_count(line: &str) -> usize {
    let Some(open) = line.find('"') else {
        return 0;
    };
    line[open + 1..]
        .chars()
        .take_while(|c| *c != '"')
        .filter(|c| *c == ',')
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quoted_commas_are_counted_up_to_the_closing_quote_or_the_end() {
        assert_eq!(quoted_comma_count("a,b"), 0);
        assert_eq!(quoted_comma_count("a,\"b,c,d\",e"), 2);
        assert_eq!(quoted_comma_count("a,\"b, c,,d"), 3);
    }
}
