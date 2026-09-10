// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! FASTA sequence records, with a streaming reader for large protein databases.
//! Sequence whitespace is removed; ASCII letters, `*`, `-`, and `.` are accepted.
//! Alphabet interpretation belongs to the chemistry layer.

use super::{parse_error, single_line};
use crate::{Error, Result};
use std::io::{BufRead, Write};

/// Identifier, description and sequence in one FASTA record.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FASTAEntry {
    /// First whitespace-delimited token following `>`.
    pub identifier: String,
    /// Remaining header text.
    pub description: String,
    /// Sequence with whitespace removed, preserving letter case.
    pub sequence: String,
}

/// Iterator over FASTA entries. Stops after the first error.
pub struct FastaReader<R> {
    reader: R,
    pending: Option<(usize, String)>,
    line: usize,
    finished: bool,
}

impl<R: BufRead> FastaReader<R> {
    /// Wrap a buffered stream without reading it yet.
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            pending: None,
            line: 0,
            finished: false,
        }
    }

    fn read_line(&mut self) -> Result<Option<String>> {
        let mut text = String::new();
        if self.reader.read_line(&mut text)? == 0 {
            return Ok(None);
        }
        self.line += 1;
        if self.line == 1 {
            text = text.trim_start_matches('\u{feff}').to_owned();
        }
        Ok(Some(text.trim().to_owned()))
    }

    fn read_entry(&mut self) -> Result<Option<FASTAEntry>> {
        let (header_line, header) = match self.pending.take() {
            Some(value) => value,
            None => loop {
                let Some(text) = self.read_line()? else {
                    return Ok(None);
                };
                if !text.is_empty() && !text.starts_with(';') {
                    break (self.line, text);
                }
            },
        };
        let header = header
            .strip_prefix('>')
            .ok_or_else(|| parse_error(header_line, "expected FASTA header beginning with >"))?
            .trim();
        let end = header.find(char::is_whitespace).unwrap_or(header.len());
        let identifier = header[..end].to_owned();
        if identifier.is_empty() {
            return Err(parse_error(header_line, "empty FASTA identifier"));
        }
        let mut entry = FASTAEntry {
            identifier,
            description: header[end..].trim().to_owned(),
            sequence: String::new(),
        };
        while let Some(text) = self.read_line()? {
            if text.starts_with('>') {
                self.pending = Some((self.line, text));
                break;
            }
            if text.starts_with(';') {
                continue;
            }
            for c in text.chars().filter(|c| !c.is_ascii_whitespace()) {
                if !valid_residue(c) {
                    return Err(parse_error(
                        self.line,
                        format!("invalid FASTA sequence character {c:?}"),
                    ));
                }
                entry.sequence.push(c);
            }
        }
        if entry.sequence.is_empty() {
            return Err(parse_error(header_line, "empty FASTA sequence"));
        }
        Ok(Some(entry))
    }
}

impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = Result<FASTAEntry>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        match self.read_entry() {
            Ok(Some(entry)) => Some(Ok(entry)),
            Ok(None) => {
                self.finished = true;
                None
            }
            Err(error) => {
                self.finished = true;
                Some(Err(error))
            }
        }
    }
}

/// Read all FASTA entries. Use [`FastaReader`] to avoid retaining the database.
pub fn read(reader: impl BufRead) -> Result<Vec<FASTAEntry>> {
    FastaReader::new(reader).collect()
}

/// Write FASTA with 80-character sequence lines.
pub fn write(mut writer: impl Write, entries: &[FASTAEntry]) -> Result<()> {
    for entry in entries {
        if entry.identifier.is_empty()
            || entry.identifier.chars().any(char::is_whitespace)
            || !single_line(&entry.identifier)
            || entry.identifier.contains('>')
            || !single_line(&entry.description)
            || entry.sequence.is_empty()
            || !entry.sequence.chars().all(valid_residue)
        {
            return Err(Error::InvalidValue(
                "invalid FASTA identifier, description or sequence".into(),
            ));
        }
    }
    for entry in entries {
        write!(writer, ">{}", entry.identifier)?;
        if !entry.description.is_empty() {
            write!(writer, " {}", entry.description)?;
        }
        writeln!(writer)?;
        for chunk in entry.sequence.as_bytes().chunks(80) {
            writer.write_all(chunk)?;
            writeln!(writer)?;
        }
    }
    writer.flush()?;
    Ok(())
}

fn valid_residue(c: char) -> bool {
    c.is_ascii_alphabetic() || matches!(c, '*' | '-' | '.')
}
