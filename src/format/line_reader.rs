// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Bounded streaming line reader with the line-end rules of the source's
//! `TextFile::getLine` (`FORMAT/TextFile.cpp`): a line ends at LF, at CR, or at
//! CRLF, and the terminator is removed. `LineReader::lf_only` instead follows
//! `std::getline`, which ends a line at LF only and keeps a CR before it.
//!
//! [`crate::format::text::TextFile`] holds a whole file under a hard ceiling of
//! one million lines, which spectral libraries and long chromatograms exceed.
//! This reader streams instead: it keeps one line in memory, and its only
//! ceilings are the line length and the total input size, both caller-chosen.

use crate::{Error, Result};
use std::io::BufRead;

/// Resource ceilings for one streamed input.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LineLimits {
    /// Longest accepted line, in bytes, excluding its terminator.
    pub max_line_bytes: usize,
    /// Largest accepted input, in bytes, including terminators.
    pub max_input_bytes: u64,
}

impl Default for LineLimits {
    /// 16 MiB per line and 4 GiB of input.
    fn default() -> Self {
        Self {
            max_line_bytes: 16 * 1024 * 1024,
            max_input_bytes: 4 * 1024 * 1024 * 1024,
        }
    }
}

/// Streaming reader of source-terminated lines.
pub(crate) struct LineReader<R> {
    reader: R,
    limits: LineLimits,
    consumed: u64,
    line_number: usize,
    carriage_return_ends_line: bool,
}

impl<R: BufRead> LineReader<R> {
    /// Lines end at LF, CR or CRLF (`TextFile::getLine`).
    pub(crate) fn new(reader: R, limits: LineLimits) -> Self {
        Self {
            reader,
            limits,
            consumed: 0,
            line_number: 0,
            carriage_return_ends_line: true,
        }
    }

    /// Lines end at LF only and keep any CR (`std::getline`).
    pub(crate) fn lf_only(reader: R, limits: LineLimits) -> Self {
        Self {
            carriage_return_ends_line: false,
            ..Self::new(reader, limits)
        }
    }

    /// One-based number of the line last returned; zero before the first.
    pub(crate) fn line_number(&self) -> usize {
        self.line_number
    }

    fn consume(&mut self, bytes: usize) -> Result<()> {
        self.consumed = self
            .consumed
            .checked_add(bytes as u64)
            .filter(|&n| n <= self.limits.max_input_bytes)
            .ok_or_else(|| Error::InvalidValue("input byte limit exceeded".into()))?;
        Ok(())
    }

    /// Reads the next line into `line`, without its terminator. Returns
    /// `false` at the end of input when no byte remains, where the source's
    /// `getLine` sets the stream's fail state; `line` is then empty.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when reading fails and [`Error::InvalidValue`]
    /// when a line or the input exceeds the limits.
    pub(crate) fn next_bytes(&mut self, line: &mut Vec<u8>) -> Result<bool> {
        line.clear();
        if self.reader.fill_buf()?.is_empty() {
            return Ok(false);
        }
        self.line_number += 1;
        loop {
            let buffer = self.reader.fill_buf()?;
            if buffer.is_empty() {
                return Ok(true);
            }
            let cr = self.carriage_return_ends_line;
            let end = buffer
                .iter()
                .position(|&b| b == b'\n' || (cr && b == b'\r'));
            let content = end.unwrap_or(buffer.len());
            if content > self.limits.max_line_bytes.saturating_sub(line.len()) {
                return Err(Error::InvalidValue(format!(
                    "line {} exceeds the line length limit of {} bytes",
                    self.line_number, self.limits.max_line_bytes
                )));
            }
            line.try_reserve(content)
                .map_err(|_| Error::InvalidValue("line allocation failed".into()))?;
            line.extend_from_slice(&buffer[..content]);
            let carriage_return = end.is_some_and(|i| buffer[i] == b'\r');
            let used = content + usize::from(end.is_some());
            self.reader.consume(used);
            self.consume(used)?;
            if carriage_return && self.reader.fill_buf()?.first() == Some(&b'\n') {
                self.reader.consume(1);
                self.consume(1)?;
            }
            if end.is_some() {
                return Ok(true);
            }
        }
    }

    /// [`Self::next_bytes`] for text that must be UTF-8.
    ///
    /// # Errors
    ///
    /// As [`Self::next_bytes`], and [`Error::Parse`] when the line is not
    /// UTF-8.
    pub(crate) fn next_line(&mut self, line: &mut String) -> Result<bool> {
        let mut bytes = std::mem::take(line).into_bytes();
        let more = self.next_bytes(&mut bytes)?;
        *line = String::from_utf8(bytes).map_err(|_| Error::Parse {
            line: self.line_number,
            message: "line is not valid UTF-8".into(),
        })?;
        Ok(more)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lines(input: &[u8]) -> Vec<Vec<u8>> {
        let mut reader = LineReader::new(input, LineLimits::default());
        let mut out = Vec::new();
        let mut line = Vec::new();
        while reader.next_bytes(&mut line).unwrap() {
            out.push(line.clone());
        }
        out
    }

    #[test]
    fn terminators_follow_text_file_get_line() {
        assert_eq!(
            lines(b"a\r\nb\rc\nd"),
            vec![b"a".to_vec(), b"b".to_vec(), b"c".to_vec(), b"d".to_vec()]
        );
        assert_eq!(lines(b"a\n\n"), vec![b"a".to_vec(), Vec::<u8>::new()]);
        assert!(lines(b"").is_empty());
        assert_eq!(lines(b"\r\r\n"), vec![Vec::<u8>::new(), Vec::new()]);
    }

    #[test]
    fn lf_only_keeps_carriage_returns() {
        let mut reader = LineReader::lf_only(&b"a\r\nb\rc\nd"[..], LineLimits::default());
        let mut line = Vec::new();
        let mut out = Vec::new();
        while reader.next_bytes(&mut line).unwrap() {
            out.push(line.clone());
        }
        assert_eq!(out, vec![b"a\r".to_vec(), b"b\rc".to_vec(), b"d".to_vec()]);
    }

    #[test]
    fn limits_are_enforced() {
        let limits = LineLimits {
            max_line_bytes: 3,
            max_input_bytes: 100,
        };
        let mut reader = LineReader::new(&b"abcd\n"[..], limits);
        let mut line = Vec::new();
        assert!(reader.next_bytes(&mut line).is_err());
        let limits = LineLimits {
            max_line_bytes: 10,
            max_input_bytes: 4,
        };
        let mut reader = LineReader::new(&b"ab\ncd\n"[..], limits);
        assert!(reader.next_bytes(&mut line).unwrap());
        assert!(reader.next_bytes(&mut line).is_err());
    }
}
