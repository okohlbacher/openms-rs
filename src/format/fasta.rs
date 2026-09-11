// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Bounded FASTA streams and the complete file lifecycle. Parsing follows the
//! source byte rules, not a protein alphabet: only space, tab, CR and LF are
//! removed from sequences. Rust records require valid UTF-8. See FASTA_SUPPORT.

use super::{parse_error, single_line};
use crate::concept::progress_logger::ProgressLogger;
use crate::{Error, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

/// Identifier, description and sequence in one FASTA record.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FASTAEntry {
    pub identifier: String,
    pub description: String,
    pub sequence: String,
}
impl FASTAEntry {
    pub fn new(
        identifier: impl Into<String>,
        description: impl Into<String>,
        sequence: impl Into<String>,
    ) -> Self {
        Self {
            identifier: identifier.into(),
            description: description.into(),
            sequence: sequence.into(),
        }
    }
    pub fn header_matches(&self, other: &Self) -> bool {
        self.identifier == other.identifier && self.description == other.description
    }
    pub fn sequence_matches(&self, other: &Self) -> bool {
        self.sequence == other.sequence
    }
}

/// Per-stream cumulative bounds, including reads repeated after a seek.
/// `max_entry_bytes` counts identifier, description and sequence bytes together.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FastaOptions {
    pub max_input_bytes: u64,
    pub max_output_bytes: u64,
    pub max_entry_bytes: usize,
    pub max_records: usize,
    pub max_work: u64,
}
impl Default for FastaOptions {
    fn default() -> Self {
        Self {
            max_input_bytes: 512 * 1024 * 1024,
            max_output_bytes: 512 * 1024 * 1024,
            max_entry_bytes: 16 * 1024 * 1024,
            max_records: 1_000_000,
            max_work: 4_000_000_000,
        }
    }
}
fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn spend(remaining: &mut u64, amount: u64) -> Result<()> {
    *remaining = remaining
        .checked_sub(amount)
        .ok_or_else(|| invalid("FASTA work limit exceeded"))?;
    Ok(())
}
fn whitespace(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\r' | b'\n')
}

/// Owning buffered reader. An error fuses iteration until an explicit seek.
/// Only the current record is retained; no line-length allocation is required.
pub struct FastaReader<R> {
    reader: R,
    options: FastaOptions,
    input_bytes: u64,
    remaining_work: u64,
    entries_read: usize,
    line: usize,
    initialized: bool,
    eof: bool,
    finished: bool,
}
impl<R: BufRead> FastaReader<R> {
    pub fn new(reader: R) -> Self {
        Self::with_options(reader, FastaOptions::default())
    }
    pub fn with_options(reader: R, options: FastaOptions) -> Self {
        Self {
            reader,
            options,
            input_bytes: 0,
            remaining_work: options.max_work,
            entries_read: 0,
            line: 1,
            initialized: false,
            eof: false,
            finished: false,
        }
    }
    pub fn into_inner(self) -> R {
        self.reader
    }
    pub fn entries_read(&self) -> usize {
        self.entries_read
    }
    fn peek(&mut self) -> Result<Option<u8>> {
        spend(&mut self.remaining_work, 1)?;
        Ok(self.reader.fill_buf()?.first().copied())
    }
    fn take(&mut self) -> Result<Option<u8>> {
        let byte = self.peek()?;
        if let Some(byte) = byte {
            if self.input_bytes >= self.options.max_input_bytes {
                return Err(invalid("FASTA input byte limit exceeded"));
            }
            self.input_bytes += 1;
            self.reader.consume(1);
            if byte == b'\n' {
                self.line = self.line.saturating_add(1);
            }
        }
        Ok(byte)
    }
    fn initialize(&mut self) -> Result<()> {
        if self.initialized {
            return Ok(());
        }
        self.initialized = true;
        // The existing native adapter accepts a leading UTF-8 BOM.
        if self.peek()? == Some(0xef) {
            for expected in [0xef, 0xbb, 0xbf] {
                if self.take()? != Some(expected) {
                    return Err(parse_error(self.line, "invalid FASTA byte-order mark"));
                }
            }
        }
        loop {
            match self.peek()? {
                Some(b'#') => loop {
                    match self.take()? {
                        Some(b'\n') => break,
                        None => {
                            self.eof = true;
                            return Ok(());
                        }
                        _ => {}
                    }
                },
                Some(byte) if whitespace(byte) => {
                    self.take()?;
                }
                _ => return Ok(()),
            }
        }
    }
    fn push(&mut self, field: &mut Vec<u8>, byte: u8, retained: &mut usize) -> Result<()> {
        if *retained >= self.options.max_entry_bytes {
            return Err(invalid("FASTA entry byte limit exceeded"));
        }
        spend(&mut self.remaining_work, 2)?; // append and final UTF-8 validation
        *retained += 1;
        field.push(byte);
        Ok(())
    }
    fn parse_entry(&mut self) -> Result<Option<FASTAEntry>> {
        self.initialize()?;
        if self.eof {
            return Ok(None);
        }
        while self.peek()?.is_some_and(whitespace) {
            self.take()?;
        }
        let header_line = self.line;
        if self.take()? != Some(b'>') {
            return Err(parse_error(
                header_line,
                "expected FASTA header beginning with >",
            ));
        }
        if self.entries_read >= self.options.max_records {
            return Err(invalid("FASTA record limit exceeded"));
        }
        let (mut id, mut description, mut sequence) = (Vec::new(), Vec::new(), Vec::new());
        let mut retained = 0;
        let has_description = loop {
            match self.take()? {
                Some(b' ' | b'\t') if !id.is_empty() => break true,
                Some(b' ' | b'\t' | b'\r') => {}
                Some(b'\n') => break false,
                None => {
                    self.eof = true;
                    return Err(parse_error(header_line, "unterminated FASTA header"));
                }
                Some(byte) => self.push(&mut id, byte, &mut retained)?,
            }
        };
        if id.is_empty() {
            return Err(parse_error(header_line, "empty FASTA identifier"));
        }
        if has_description {
            loop {
                match self.take()? {
                    Some(b'\n') => break,
                    Some(b'\r' | b'\t') => {}
                    None => {
                        self.eof = true;
                        return Err(parse_error(header_line, "unterminated FASTA description"));
                    }
                    Some(byte) => self.push(&mut description, byte, &mut retained)?,
                }
            }
        }
        loop {
            match self.take()? {
                Some(b'\n') => {
                    if self.peek()? == Some(b'>') {
                        break;
                    }
                }
                Some(b' ' | b'\t' | b'\r') => {}
                None => {
                    self.eof = true;
                    break;
                }
                Some(byte) => self.push(&mut sequence, byte, &mut retained)?,
            }
        }
        if sequence.is_empty() {
            return Err(parse_error(header_line, "empty FASTA sequence"));
        }
        let text = |bytes| {
            String::from_utf8(bytes)
                .map_err(|_| parse_error(header_line, "FASTA field is not valid UTF-8"))
        };
        let entry = FASTAEntry {
            identifier: text(id)?,
            description: text(description)?,
            sequence: text(sequence)?,
        };
        self.entries_read += 1;
        Ok(Some(entry))
    }
    pub fn next_entry(&mut self) -> Result<Option<FASTAEntry>> {
        if self.finished {
            return Ok(None);
        }
        match self.parse_entry() {
            Ok(Some(entry)) => Ok(Some(entry)),
            Ok(None) => {
                self.finished = true;
                Ok(None)
            }
            Err(error) => {
                self.finished = true;
                Err(error)
            }
        }
    }
    /// Changes the destination only when a complete entry was read.
    pub fn read_next(&mut self, output: &mut FASTAEntry) -> Result<bool> {
        if let Some(entry) = self.next_entry()? {
            *output = entry;
            Ok(true)
        } else {
            Ok(false)
        }
    }
    /// Like the source `peek`, testing EOF also makes subsequent reads return none.
    pub fn at_end(&mut self) -> Result<bool> {
        self.initialize()?;
        if self.peek()?.is_none() {
            self.eof = true;
        }
        Ok(self.eof)
    }
}
impl<R: BufRead + Seek> FastaReader<R> {
    /// Byte position of the next input byte; `None` is the source EOF sentinel -1.
    pub fn position(&mut self) -> Result<Option<u64>> {
        self.initialize()?;
        if self.eof {
            Ok(None)
        } else {
            Ok(Some(self.reader.stream_position()?))
        }
    }
    /// Seeks within the stream, clearing EOF/errors but retaining cumulative limits.
    /// PEFF/BOM initialization is not repeated. Negative positions are unrepresentable.
    pub fn set_position(&mut self, position: u64) -> Result<bool> {
        spend(&mut self.remaining_work, 3)?;
        let previous = self.reader.stream_position()?;
        let size = self.reader.seek(SeekFrom::End(0))?;
        if position > size {
            self.reader.seek(SeekFrom::Start(previous))?;
            return Ok(false);
        }
        self.reader.seek(SeekFrom::Start(position))?;
        self.initialized = true;
        self.eof = false;
        self.finished = false;
        self.line = 1; // A random byte seek has no inexpensive source line-number oracle.
        Ok(true)
    }
}
impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = Result<FASTAEntry>;
    fn next(&mut self) -> Option<Self::Item> {
        self.next_entry().transpose()
    }
}

#[derive(Clone, Copy)]
struct WriteBudget {
    records: usize,
    bytes: u64,
    work: u64,
}
impl WriteBudget {
    fn new(options: FastaOptions) -> Self {
        Self {
            records: 0,
            bytes: 0,
            work: options.max_work,
        }
    }
    fn entry(&mut self, entry: &FASTAEntry, options: FastaOptions) -> Result<()> {
        let payload = entry
            .identifier
            .len()
            .checked_add(entry.description.len())
            .and_then(|n| n.checked_add(entry.sequence.len()))
            .ok_or_else(|| invalid("FASTA entry size overflow"))?;
        if payload > options.max_entry_bytes {
            return Err(invalid("FASTA entry byte limit exceeded"));
        }
        if self.records >= options.max_records {
            return Err(invalid("FASTA record limit exceeded"));
        }
        // Charge validation and emission, including one visit for an empty record.
        spend(
            &mut self.work,
            (payload as u64)
                .checked_mul(4)
                .and_then(|n| n.checked_add(1))
                .ok_or_else(|| invalid("FASTA work size overflow"))?,
        )?;
        if entry.identifier.is_empty()
            || entry.identifier.chars().any(char::is_whitespace)
            || !single_line(&entry.identifier)
            || entry.identifier.contains('>')
            || !single_line(&entry.description)
        {
            return Err(invalid("invalid FASTA identifier or description"));
        }
        let amount = (payload as u64)
            .checked_add(3)
            .and_then(|n| n.checked_add(entry.sequence.len().div_ceil(80) as u64))
            .ok_or_else(|| invalid("FASTA output size overflow"))?;
        self.bytes = self
            .bytes
            .checked_add(amount)
            .filter(|n| *n <= options.max_output_bytes)
            .ok_or_else(|| invalid("FASTA output byte limit exceeded"))?;
        self.records += 1;
        Ok(())
    }
}
fn emit(writer: &mut impl Write, entry: &FASTAEntry) -> Result<()> {
    writeln!(writer, ">{} {}", entry.identifier, entry.description)?;
    for chunk in entry.sequence.as_bytes().chunks(80) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}
/// Streaming writer. Each entry is validated before output; an I/O failure
/// disables later writes. Call `finish` to observe flush failures.
pub struct FastaWriter<W: Write> {
    writer: Option<W>,
    options: FastaOptions,
    budget: WriteBudget,
    failed: bool,
}
impl<W: Write> FastaWriter<W> {
    pub fn new(writer: W) -> Self {
        Self::with_options(writer, FastaOptions::default())
    }
    pub fn with_options(writer: W, options: FastaOptions) -> Self {
        Self {
            writer: Some(writer),
            options,
            budget: WriteBudget::new(options),
            failed: false,
        }
    }
    pub fn write_entry(&mut self, entry: &FASTAEntry) -> Result<()> {
        if self.failed {
            return Err(invalid("FASTA writer has failed"));
        }
        let mut budget = self.budget;
        budget.entry(entry, self.options)?;
        self.budget = budget;
        let result = emit(self.writer.as_mut().expect("owned writer"), entry);
        if result.is_err() {
            self.failed = true;
        }
        result
    }
    pub fn finish(mut self) -> Result<W> {
        let mut writer = self.writer.take().expect("owned writer");
        if self.failed {
            return Err(invalid("FASTA writer has failed"));
        }
        writer.flush()?;
        Ok(writer)
    }
}
impl<W: Write> Drop for FastaWriter<W> {
    fn drop(&mut self) {
        if !self.failed {
            if let Some(writer) = &mut self.writer {
                let _ = writer.flush();
            }
        }
    }
}

pub fn read(reader: impl BufRead) -> Result<Vec<FASTAEntry>> {
    read_with_options(reader, FastaOptions::default())
}
pub fn read_with_options(reader: impl BufRead, options: FastaOptions) -> Result<Vec<FASTAEntry>> {
    FastaReader::with_options(reader, options).collect()
}
/// Source 80-byte lines, including the space after an identifier with no description.
/// All records are checked before the first write; empty sequences remain writable.
pub fn write(writer: impl Write, entries: &[FASTAEntry]) -> Result<()> {
    write_with_options(writer, entries, FastaOptions::default())
}
pub fn write_with_options(
    mut writer: impl Write,
    entries: &[FASTAEntry],
    options: FastaOptions,
) -> Result<()> {
    let mut budget = WriteBudget::new(options);
    for entry in entries {
        budget.entry(entry, options)?;
    }
    for entry in entries {
        emit(&mut writer, entry)?;
    }
    writer.flush()?;
    Ok(())
}

/// File-oriented source lifecycle, with independent simultaneous input/output.
/// Paths are plain FASTA streams, even when a compressed suffix is supplied.
/// Configure inherited source progress behavior through the public `progress`.
#[derive(Default)]
pub struct FASTAFile {
    pub options: FastaOptions,
    pub progress: ProgressLogger,
    input: Option<FastaReader<BufReader<File>>>,
    output: Option<FastaWriter<BufWriter<File>>>,
    input_size: u64,
}
impl FASTAFile {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn read_start(&mut self, path: impl AsRef<Path>) -> Result<()> {
        let file = File::open(path)?;
        let size = file.metadata()?.len();
        let mut reader = FastaReader::with_options(BufReader::new(file), self.options);
        reader.initialize()?;
        self.input = Some(reader);
        self.input_size = size;
        Ok(())
    }
    pub fn read_start_with_progress(&mut self, path: impl AsRef<Path>, label: &str) -> Result<()> {
        self.read_start(path)?;
        self.progress.start_progress(
            0,
            i64::try_from(self.input_size)
                .map_err(|_| invalid("FASTA file size exceeds progress range"))?,
            label,
        )
    }
    pub fn read_next(&mut self, entry: &mut FASTAEntry) -> Result<bool> {
        let reader = self
            .input
            .as_mut()
            .ok_or_else(|| invalid("FASTA input is not open"))?;
        let result = reader.read_next(entry)?;
        if result {
            self.progress
                .set_progress(progress_position(reader.position()?)?)?;
        }
        Ok(result)
    }
    pub fn read_next_with_progress(&mut self, entry: &mut FASTAEntry) -> Result<bool> {
        if self.read_next(entry)? {
            let position = progress_position(self.position()?)?;
            self.progress.set_progress(position)?;
            Ok(true)
        } else {
            self.progress.end_progress(0)?;
            Ok(false)
        }
    }
    pub fn position(&mut self) -> Result<Option<u64>> {
        self.input
            .as_mut()
            .ok_or_else(|| invalid("FASTA input is not open"))?
            .position()
    }
    pub fn set_position(&mut self, position: u64) -> Result<bool> {
        self.input
            .as_mut()
            .ok_or_else(|| invalid("FASTA input is not open"))?
            .set_position(position)
    }
    pub fn at_end(&mut self) -> Result<bool> {
        self.input
            .as_mut()
            .ok_or_else(|| invalid("FASTA input is not open"))?
            .at_end()
    }
    pub fn write_start(&mut self, path: impl AsRef<Path>) -> Result<()> {
        if self.output.is_some() {
            return Err(invalid("FASTA output is already open"));
        }
        check_extension(path.as_ref())?;
        self.output = Some(FastaWriter::with_options(
            BufWriter::new(File::create(path)?),
            self.options,
        ));
        Ok(())
    }
    pub fn write_next(&mut self, entry: &FASTAEntry) -> Result<()> {
        self.output
            .as_mut()
            .ok_or_else(|| invalid("FASTA output is not open"))?
            .write_entry(entry)
    }
    pub fn write_end(&mut self) -> Result<()> {
        if let Some(writer) = self.output.take() {
            writer.finish()?;
        }
        Ok(())
    }
    /// Uses a separate input stream, leaving a resident read/write session intact.
    pub fn load(&mut self, path: impl AsRef<Path>) -> Result<Vec<FASTAEntry>> {
        self.progress.start_progress(0, 1, "Loading FASTA file")?;
        let entries = read_with_options(BufReader::new(File::open(path)?), self.options)?;
        self.progress.end_progress(0)?;
        Ok(entries)
    }
    /// Replaces the destination only after the complete file succeeds.
    pub fn load_into(
        &mut self,
        path: impl AsRef<Path>,
        entries: &mut Vec<FASTAEntry>,
    ) -> Result<()> {
        *entries = self.load(path)?;
        Ok(())
    }
    /// Checks all entries before opening the destination. Transport I/O failures
    /// can leave a partial file, as in the source; caller-owned sessions are intact.
    pub fn store(&mut self, path: impl AsRef<Path>, entries: &[FASTAEntry]) -> Result<()> {
        check_extension(path.as_ref())?;
        let mut budget = WriteBudget::new(self.options);
        for entry in entries {
            budget.entry(entry, self.options)?;
        }
        self.progress.start_progress(
            0,
            i64::try_from(entries.len())
                .map_err(|_| invalid("FASTA record count exceeds progress range"))?,
            "Writing FASTA file",
        )?;
        let mut writer = BufWriter::new(File::create(path)?);
        for entry in entries {
            emit(&mut writer, entry)?;
            self.progress.next_progress()?;
        }
        writer.flush()?;
        self.progress.end_progress(0)
    }
}
fn progress_position(position: Option<u64>) -> Result<i64> {
    position.map_or(Ok(-1), |n| {
        i64::try_from(n).map_err(|_| invalid("FASTA position exceeds progress range"))
    })
}
fn check_extension(path: &Path) -> Result<()> {
    use super::file_types::{FileType, has_valid_extension};
    if !has_valid_extension(&path.to_string_lossy(), FileType::Fasta) {
        return Err(invalid("invalid FASTA output file extension"));
    }
    Ok(())
}
