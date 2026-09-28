// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Inspect result files and Inspect's trie sequence databases.
//!
//! Ports `FORMAT/InspectOutfile.h`/`.cpp`: the tab-separated result reader
//! ([`InspectOutfile::load`]), which takes precursor retention times and m/z
//! values from the searched spectrum files, and the trie database tools
//! Inspect searches with ([`InspectOutfile::generate_trie_db`],
//! [`InspectOutfile::compress_trie_db`], [`InspectOutfile::sequences`]).
//!
//! A trie database is the protein sequences separated by `*`. Its index file
//! holds one [`RECORD_LENGTH`]-byte record per protein: the protein's byte
//! position in the source database (8 bytes, little endian), its position in
//! the trie database (4 bytes, little endian) and its name (80 bytes, NUL
//! padded, `>` written as `}`).
//!
//! See `docs/SEARCH_ENGINE_FILES_SUPPORT.md` for the API mapping.

use crate::chemistry::AASequence;
use crate::data_structures::DateTime;
use crate::data_structures::string_utils::{split, to_float, to_int32, trim};
use crate::format::file_handler::type_by_content;
use crate::format::search_engine_text::{
    MAX_TEXT_BYTES, NPOS, at_line, clean_line, flanking_residue, getline_lines,
    parse_permissive_sequence, read_text,
};
use crate::format::{FileHandler, FileType, PeakFileOptions};
use crate::identification::{
    PeptideEvidence, PeptideHit, PeptideIdentification, ProteinHit, ProteinIdentification,
};
use crate::{Error, MSExperiment, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

/// Length of the source-database position in an index record.
pub const DB_POS_LENGTH: usize = 8;
/// Length of the trie-database position in an index record.
pub const TRIE_DB_POS_LENGTH: usize = 4;
/// Length of the protein name in an index record.
pub const PROTEIN_NAME_LENGTH: usize = 80;
/// Length of one index record.
pub const RECORD_LENGTH: usize = DB_POS_LENGTH + TRIE_DB_POS_LENGTH + PROTEIN_NAME_LENGTH;
/// The separator of the sequences in a trie database.
pub const TRIE_DELIMITER: u8 = b'*';
/// The score type of the peptide identifications.
pub const SCORE_TYPE: &str = "Inspect";
/// The longest line or sequence the database tools read, in bytes.
pub const MAX_RECORD_BYTES: usize = 256 * 1024 * 1024;

/// The columns of an Inspect result header (source `readOutHeader`'s output
/// parameters), zero based.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct InspectOutHeader {
    /// `#SpectrumFile`.
    pub spectrum_file_column: usize,
    /// `Scan#`.
    pub scan_column: usize,
    /// `Annotation`, the peptide.
    pub peptide_column: usize,
    /// `Protein`.
    pub protein_column: usize,
    /// `Charge`.
    pub charge_column: usize,
    /// `MQScore`.
    pub mq_score_column: usize,
    /// `p-value`.
    pub p_value_column: usize,
    /// `RecordNumber`, the protein's record in the database.
    pub record_number_column: usize,
    /// `DBFilePos`.
    pub db_file_pos_column: usize,
    /// `SpecFilePos`.
    pub spec_file_pos_column: usize,
    /// The number of columns of the header, which every result line must have.
    pub number_of_columns: usize,
}

/// The line labels of a sequence database (source `getLabels`' output
/// parameters).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DatabaseLabels {
    /// Starts the line with the accession.
    pub ac_label: String,
    /// Starts the sequence.
    pub sequence_start_label: String,
    /// Ends the sequence.
    pub sequence_end_label: String,
    /// Starts a comment line.
    pub comment_label: String,
    /// Starts the line with the species.
    pub species_label: String,
}

/// Spectrum files with the (peptide identification index, scan number) pairs
/// whose precursor they supply (the source's
/// `files_and_peptide_identification_with_scan_number`).
pub type ScanAssignments = Vec<(String, Vec<(usize, usize)>)>;

/// Inspect result and trie database reader. It holds no state, so every two
/// readers are equal, as in the source.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct InspectOutfile;

fn not_found(path: &Path) -> Error {
    Error::Io(std::io::Error::new(
        std::io::ErrorKind::NotFound,
        format!("the file '{}' could not be found", path.display()),
    ))
}

/// Source `File::exists` check before an open, then the open's own error.
fn open(path: &Path) -> Result<File> {
    if !path.exists() {
        return Err(not_found(path));
    }
    Ok(File::open(path)?)
}

/// Source `Exception::ParseError(…, expression, message)`: `what()` is
/// "`message` in: `expression`".
fn parse_error(message: impl std::fmt::Display, expression: impl std::fmt::Display) -> Error {
    Error::Parse {
        line: 0,
        message: format!("{message} in: {expression}"),
    }
}

fn check_p_value(p_value_threshold: f64) -> Result<()> {
    if (0.0..=1.0).contains(&p_value_threshold) {
        Ok(())
    } else {
        Err(Error::InvalidValue(
            "The parameters 'p_value_threshold' must be >= 0 and <=1 !".into(),
        ))
    }
}

/// Read a result file and return its lines after the header check.
fn result_lines(path: &Path) -> Result<String> {
    let text = read_text(open(path)?, MAX_TEXT_BYTES)?;
    if text.is_empty() {
        return Err(Error::Parse {
            line: 0,
            message: format!("the file '{}' is empty", path.display()),
        });
    }
    Ok(text)
}

fn field_float(text: &str, line: usize) -> Result<f64> {
    to_float(text)
        .map(f64::from)
        .map_err(|message| at_line(Error::InvalidValue(message), line))
}

/// Source `(Size) StringUtils::toInt32(text)`; a negative number, which the
/// source converts to a huge count, is refused.
fn field_count(text: &str, line: usize) -> Result<usize> {
    let value = to_int32(text).map_err(|message| at_line(Error::InvalidValue(message), line))?;
    usize::try_from(value).map_err(|_| Error::Parse {
        line,
        message: format!("negative record or scan number '{text}'"),
    })
}

fn trim_bytes(line: &[u8]) -> &[u8] {
    let is_space = |b: &u8| matches!(b, b' ' | b'\t' | b'\n' | b'\r');
    let start = line.iter().position(|b| !is_space(b)).unwrap_or(line.len());
    let end = line
        .iter()
        .rposition(|b| !is_space(b))
        .map_or(start, |p| p + 1);
    &line[start..end]
}

/// Read one line (without its line feed) into `line`, at most
/// [`MAX_RECORD_BYTES`]. Returns the bytes consumed, 0 at the end.
fn read_line(reader: &mut impl BufRead, line: &mut Vec<u8>) -> Result<usize> {
    line.clear();
    let mut consumed = 0usize;
    loop {
        let buffer = reader.fill_buf()?;
        if buffer.is_empty() {
            return Ok(consumed);
        }
        let (take, done) = match buffer.iter().position(|&b| b == b'\n') {
            Some(pos) => (pos + 1, true),
            None => (buffer.len(), false),
        };
        let content = if done {
            &buffer[..take - 1]
        } else {
            &buffer[..take]
        };
        if line.len() + content.len() > MAX_RECORD_BYTES {
            return Err(Error::InvalidValue(format!(
                "database line longer than {MAX_RECORD_BYTES} bytes"
            )));
        }
        line.extend_from_slice(content);
        reader.consume(take);
        consumed += take;
        if done {
            return Ok(consumed);
        }
    }
}

/// Source `istream::get(streambuf&, delim)`: the bytes up to (not including)
/// `delim` or the end, at most [`MAX_RECORD_BYTES`].
fn read_until_delimiter(reader: &mut impl BufRead, delimiter: u8) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    loop {
        let buffer = reader.fill_buf()?;
        if buffer.is_empty() {
            return Ok(out);
        }
        let (take, done) = match buffer.iter().position(|&b| b == delimiter) {
            Some(pos) => (pos, true),
            None => (buffer.len(), false),
        };
        if out.len() + take > MAX_RECORD_BYTES {
            return Err(Error::InvalidValue(format!(
                "trie database sequence longer than {MAX_RECORD_BYTES} bytes"
            )));
        }
        out.extend_from_slice(&buffer[..take]);
        reader.consume(take);
        if done {
            return Ok(out);
        }
    }
}

/// Source `istream::ignore(n, delim)` with an unbounded `n`: skip up to and
/// including the next `delim`. Returns false at the end of the input.
fn skip_past_delimiter(reader: &mut impl BufRead, delimiter: u8) -> Result<bool> {
    loop {
        let buffer = reader.fill_buf()?;
        if buffer.is_empty() {
            return Ok(false);
        }
        match buffer.iter().position(|&b| b == delimiter) {
            Some(pos) => {
                reader.consume(pos + 1);
                return Ok(true);
            }
            None => {
                let length = buffer.len();
                reader.consume(length);
            }
        }
    }
}

/// An output file the source opens with `ios::app` or `ios::trunc`, with its
/// current length (the source's `tellp`).
struct Output {
    writer: BufWriter<File>,
    position: u64,
}

impl Output {
    fn create(path: &Path, append: bool) -> Result<Self> {
        let file = if append {
            OpenOptions::new().create(true).append(true).open(path)
        } else {
            File::create(path)
        }
        .map_err(|error| {
            Error::Io(std::io::Error::new(
                error.kind(),
                format!(
                    "the file '{}' could not be created: {error}",
                    path.display()
                ),
            ))
        })?;
        let position = if append { file.metadata()?.len() } else { 0 };
        Ok(Self {
            writer: BufWriter::new(file),
            position,
        })
    }

    fn write(&mut self, bytes: &[u8]) -> Result<()> {
        self.writer.write_all(bytes)?;
        self.position += bytes.len() as u64;
        Ok(())
    }

    fn finish(mut self) -> Result<()> {
        self.writer.flush()?;
        Ok(())
    }
}

/// The trie-database position as the index stores it: the low four bytes,
/// little endian (the source's `memcpy` of four bytes of a `Size`).
fn trie_position_bytes(position: u64) -> [u8; TRIE_DB_POS_LENGTH] {
    let bytes = position.to_le_bytes();
    [bytes[0], bytes[1], bytes[2], bytes[3]]
}

impl InspectOutfile {
    /// A reader; it holds no state.
    pub fn new() -> Self {
        Self
    }

    /// Load an Inspect result file; see [`InspectOutfile::load_at`], with
    /// the current local date in the identifier.
    ///
    /// # Errors
    ///
    /// As [`InspectOutfile::load_at`].
    pub fn load(
        &self,
        result_file: impl AsRef<Path>,
        peptide_identifications: &mut Vec<PeptideIdentification>,
        protein_identification: &mut ProteinIdentification,
        p_value_threshold: f64,
        database: Option<&Path>,
    ) -> Result<Vec<usize>> {
        self.load_at(
            result_file,
            peptide_identifications,
            protein_identification,
            p_value_threshold,
            database,
            DateTime::now(),
        )
    }

    /// Load an Inspect result file with the run date `run_date`.
    ///
    /// Lines whose `p-value` exceeds `p_value_threshold` are skipped. One
    /// peptide identification per spectrum (a new spectrum file or scan
    /// number starts a new one) is appended to `peptide_identifications`,
    /// with score type [`SCORE_TYPE`], the threshold as significance
    /// threshold and identifier `<search engine>_<date>` (`InsPecT_<date>`
    /// when `protein_identification` names no search engine). Each hit has
    /// the `MQScore` as score (read as a `float`), rank 0, the charge, the
    /// upper-case letters between the peptide's dots as sequence, and one
    /// evidence with the flanking residues and the protein accession
    /// ([`InspectOutfile::ac_and_ac_type`]). Only the last identification's
    /// hits are sorted by score, as in the source. Each database record
    /// seen adds one protein hit with the accession to
    /// `protein_identification`, whose sequence is read from the trie
    /// `database` when one is given; its identifier and date are set.
    ///
    /// The precursor retention time and m/z of each identification come from
    /// its spectrum file ([`InspectOutfile::precursor_rt_and_mz`]). The
    /// source records the scan number for an identification when the *next*
    /// spectrum starts, and records the next spectrum's scan number (and,
    /// at a file change, the next file); the port keeps that assignment, so
    /// every identification but the last takes the precursor of the
    /// spectrum after it. `docs/SEARCH_ENGINE_FILES_SUPPORT.md` records this
    /// as a source defect.
    ///
    /// Returns the (one-based) numbers of the lines whose column count
    /// differs from the header's.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a threshold outside `[0, 1]` (source
    /// `Exception::IllegalArgument`; a NaN, which the source accepts, is
    /// refused too); [`Error::Io`] when the file is missing ("the file '…'
    /// could not be found") or unreadable; [`Error::Parse`] when it is empty
    /// ("the file '…' is empty", source `Exception::FileEmpty`), when the
    /// header lacks a column ([`InspectOutfile::read_out_header`]), for a
    /// number that does not convert (source `Exception::ConversionError`), a
    /// negative record or scan number, an invalid flanking residue or
    /// sequence; and the errors of [`InspectOutfile::sequences`] and
    /// [`InspectOutfile::precursor_rt_and_mz`]. The outputs may then hold
    /// part of the file, as in the source.
    pub fn load_at(
        &self,
        result_file: impl AsRef<Path>,
        peptide_identifications: &mut Vec<PeptideIdentification>,
        protein_identification: &mut ProteinIdentification,
        p_value_threshold: f64,
        database: Option<&Path>,
        run_date: DateTime,
    ) -> Result<Vec<usize>> {
        check_p_value(p_value_threshold)?;
        let result_file = result_file.as_ref();
        let text = result_lines(result_file)?;
        let lines = getline_lines(&text);
        let header_line = trim(clean_line(lines[0]));
        let mut line_number = 1usize;

        let identifier = if protein_identification.search_engine.is_empty() {
            format!("InsPecT_{}", run_date.date_string())
        } else {
            format!(
                "{}_{}",
                protein_identification.search_engine,
                run_date.date_string()
            )
        };
        let mut files: ScanAssignments = Vec::new();
        let mut rn_position_map: BTreeMap<usize, usize> = BTreeMap::new();
        let header = self.read_out_header(&result_file.to_string_lossy(), header_line)?;

        let mut corrupted_lines = Vec::new();
        let mut spectrum_file = String::new();
        let mut scan_number = 0usize;
        let mut peptide_identification = PeptideIdentification::default();

        for raw in &lines[1..] {
            line_number += 1;
            let line = trim(clean_line(raw));
            if line.is_empty() {
                continue;
            }
            let fields = split(line, '\t');
            if fields.len() != header.number_of_columns {
                corrupted_lines.push(line_number);
                continue;
            }
            if field_float(fields[header.p_value_column], line_number)? > p_value_threshold {
                continue;
            }
            let (accession, _) = self.ac_and_ac_type(fields[header.protein_column]);

            let record_number = field_count(fields[header.record_number_column], line_number)?;
            if let std::collections::btree_map::Entry::Vacant(entry) =
                rn_position_map.entry(record_number)
            {
                entry.insert(protein_identification.hits.len());
                protein_identification.hits.push(ProteinHit {
                    accession: accession.clone(),
                    ..ProteinHit::default()
                });
            }

            let file = fields[header.spectrum_file_column];
            let scan = field_count(fields[header.scan_column], line_number)?;
            if file != spectrum_file || scan != scan_number {
                if file != spectrum_file {
                    match files.last_mut() {
                        Some((name, scans)) if scans.is_empty() => *name = file.to_owned(),
                        _ => files.push((file.to_owned(), Vec::new())),
                    }
                }
                spectrum_file = file.to_owned();
                scan_number = scan;
                if !peptide_identification.hits.is_empty() {
                    if let Some((_, scans)) = files.last_mut() {
                        scans.push((peptide_identifications.len(), scan_number));
                    }
                    peptide_identifications.push(std::mem::take(&mut peptide_identification));
                }
                peptide_identification = PeptideIdentification {
                    identifier: identifier.clone(),
                    significance_threshold: p_value_threshold,
                    score_type: SCORE_TYPE.into(),
                    ..PeptideIdentification::default()
                };
            }

            let charge = to_int32(fields[header.charge_column])
                .map_err(|message| at_line(Error::InvalidValue(message), line_number))?;
            let score = field_float(fields[header.mq_score_column], line_number)?;
            let hit = self
                .peptide_hit(fields[header.peptide_column], &accession, charge, score)
                .map_err(|error| at_line(error, line_number))?;
            peptide_identification.hits.push(hit);
        }

        if !peptide_identification.hits.is_empty() {
            if let Some((_, scans)) = files.last_mut() {
                scans.push((peptide_identifications.len(), scan_number));
            }
            peptide_identifications.push(peptide_identification);
        }
        if files.last().is_some_and(|(_, scans)| scans.is_empty()) {
            files.pop();
        }
        if let Some(last) = peptide_identifications.last_mut() {
            last.sort()?;
        }

        if !protein_identification.hits.is_empty() {
            if let Some(database) = database {
                let mut sequences = Vec::new();
                self.sequences(database, &rn_position_map, &mut sequences)?;
                for (position, sequence) in rn_position_map.values().zip(sequences) {
                    protein_identification.hits[*position].sequence = sequence;
                }
            }
        }

        self.precursor_rt_and_mz(&files, peptide_identifications)?;
        protein_identification.date_time = Some(run_date.get());
        protein_identification.identifier = identifier;
        Ok(corrupted_lines)
    }

    fn peptide_hit(
        &self,
        annotation: &str,
        accession: &str,
        charge: i32,
        score: f64,
    ) -> Result<PeptideHit> {
        let bytes = annotation.as_bytes();
        let start = annotation.find('.').map_or(0, |dot| dot + 1);
        let end = annotation.rfind('.').unwrap_or(NPOS);
        let mut evidence = PeptideEvidence {
            protein_accession: accession.to_owned(),
            ..PeptideEvidence::default()
        };
        if start >= 2 {
            evidence.aa_before = flanking_residue(bytes[start - 2], true)?;
        }
        if end != NPOS {
            evidence.aa_after = flanking_residue(bytes.get(end + 1).copied().unwrap_or(0), false)?;
        }
        // Keep the upper-case letters between the dots; with one dot or none,
        // `substr(start, end - start)` wraps around and takes the rest.
        let stop = if end == NPOS || end < start {
            bytes.len()
        } else {
            end
        };
        let residues: String = bytes[start.min(bytes.len())..stop.min(bytes.len())]
            .iter()
            .filter(|b| b.is_ascii_uppercase())
            .map(|&b| char::from(b))
            .collect();
        let sequence: AASequence = parse_permissive_sequence(&residues)?;
        Ok(PeptideHit {
            sequence,
            score,
            rank: 0,
            charge,
            evidences: vec![evidence],
            ..PeptideHit::default()
        })
    }

    /// The database record numbers of the result lines whose `p-value` is at
    /// most `p_value_threshold`, ascending and distinct (source
    /// `getWantedRecords`).
    ///
    /// # Errors
    ///
    /// As [`InspectOutfile::load_at`] for the threshold, the file and the
    /// header, and [`Error::Parse`] for a number that does not convert or a
    /// negative record number.
    pub fn wanted_records(
        &self,
        result_file: impl AsRef<Path>,
        p_value_threshold: f64,
    ) -> Result<Vec<usize>> {
        if !(0.0..=1.0).contains(&p_value_threshold) {
            return Err(Error::InvalidValue(
                "the parameters 'p_value_threshold' must be >= 0 and <=1 !".into(),
            ));
        }
        let result_file = result_file.as_ref();
        let text = result_lines(result_file)?;
        let lines = getline_lines(&text);
        let header = self.read_out_header(&result_file.to_string_lossy(), lines[0])?;
        let mut wanted = BTreeSet::new();
        for (index, raw) in lines.iter().enumerate().skip(1) {
            let line_number = index + 1;
            let line = trim(clean_line(raw));
            if line.is_empty() {
                continue;
            }
            let fields = split(line, '\t');
            if fields.len() != header.number_of_columns {
                continue;
            }
            if field_float(fields[header.p_value_column], line_number)? > p_value_threshold {
                continue;
            }
            wanted.insert(field_count(
                fields[header.record_number_column],
                line_number,
            )?);
        }
        Ok(wanted.into_iter().collect())
    }

    /// Copy the `wanted_records` of a trie database and its index to a
    /// second database and index (source `compressTrieDB`). An empty
    /// `wanted_records` copies every record. The copied index records keep
    /// their source position and name and get the position in the second
    /// database. With `append`, the second database and index are extended
    /// and every sequence, the first included, is preceded by `*`.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] "`database` in: Same filename can not be used for
    /// original and second database!" (likewise for the indices);
    /// [`Error::Io`] when an input is missing ("the file '…' could not be
    /// found") or an output cannot be created (source
    /// `Exception::UnableToCreateFile`); [`Error::Parse`] "`index` in: index
    /// file is too short!" for a record beyond the index, after the outputs
    /// were created, as in the source.
    pub fn compress_trie_db(
        &self,
        database: impl AsRef<Path>,
        index: impl AsRef<Path>,
        wanted_records: &[usize],
        snd_database: impl AsRef<Path>,
        snd_index: impl AsRef<Path>,
        append: bool,
    ) -> Result<()> {
        let (database, index) = (database.as_ref(), index.as_ref());
        let (snd_database, snd_index) = (snd_database.as_ref(), snd_index.as_ref());
        const SAME: &str = "Same filename can not be used for original and second database!";
        if database == snd_database {
            return Err(parse_error(database.display(), SAME));
        }
        if index == snd_index {
            return Err(parse_error(index.display(), SAME));
        }
        let mut database_reader = BufReader::new(open(database)?);
        let mut index_file = open(index)?;
        let index_length = index_file.metadata()?.len();
        let all: Vec<usize>;
        let wanted_records = if wanted_records.is_empty() {
            let count = usize::try_from(index_length / RECORD_LENGTH as u64).map_err(|_| {
                Error::InvalidValue("index file has more records than addressable".into())
            })?;
            all = (0..count).collect();
            &all[..]
        } else {
            wanted_records
        };

        let mut snd_db = Output::create(snd_database, append)?;
        let mut snd_idx = Output::create(snd_index, append)?;
        let mut append = append;
        let mut failed = false;
        let mut record = [0u8; RECORD_LENGTH];
        for &wanted in wanted_records {
            let end = (wanted as u64)
                .checked_add(1)
                .and_then(|n| n.checked_mul(RECORD_LENGTH as u64));
            if end.is_none_or(|end| index_length < end) {
                return Err(parse_error(index.display(), "index file is too short!"));
            }
            index_file.seek(SeekFrom::Start(wanted as u64 * RECORD_LENGTH as u64))?;
            index_file.read_exact(&mut record)?;

            if append {
                snd_db.write(&[TRIE_DELIMITER])?;
            }
            append = true;

            let mut position = [0u8; 8];
            position[..TRIE_DB_POS_LENGTH]
                .copy_from_slice(&record[DB_POS_LENGTH..DB_POS_LENGTH + TRIE_DB_POS_LENGTH]);
            let database_position = u64::from_le_bytes(position);
            record[DB_POS_LENGTH..DB_POS_LENGTH + TRIE_DB_POS_LENGTH]
                .copy_from_slice(&trie_position_bytes(snd_db.position));
            snd_idx.write(&record)?;

            // The source's stream fails once a read yields nothing, and a
            // failed stream neither seeks nor reads again.
            if !failed {
                database_reader.seek(SeekFrom::Start(database_position))?;
                let sequence = read_until_delimiter(&mut database_reader, TRIE_DELIMITER)?;
                failed = sequence.is_empty();
                snd_db.write(&sequence)?;
            }
        }
        snd_db.finish()?;
        snd_idx.finish()
    }

    /// Generate a trie database and its index from a FASTA or SwissProt
    /// database (source `generateTrieDB`); the labels come from
    /// [`InspectOutfile::labels`]. With a non-empty `species`, only proteins
    /// whose species line contains it are taken. With `append`, the outputs
    /// are extended and every sequence, the first included, is preceded by
    /// `*`.
    ///
    /// Sequence lines are trimmed at both ends and their `*` removed; as in
    /// the source, spaces inside a line are kept (a SwissProt sequence line
    /// keeps its blocks' spaces). A protein's name is the at most
    /// [`PROTEIN_NAME_LENGTH`] bytes after the accession label and the bytes
    /// below 33 (or above 127, a negative `char`) that follow it, `>`
    /// written as `}`.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the source is missing ("the file '…' could not be
    /// found") or an output cannot be created (source
    /// `Exception::UnableToCreateFile`), the errors of
    /// [`InspectOutfile::labels`], and [`Error::InvalidValue`] for a line
    /// longer than [`MAX_RECORD_BYTES`].
    pub fn generate_trie_db(
        &self,
        source_database: impl AsRef<Path>,
        database: impl AsRef<Path>,
        index: impl AsRef<Path>,
        append: bool,
        species: &str,
    ) -> Result<()> {
        let source_database = source_database.as_ref();
        let mut source = BufReader::new(open(source_database)?);
        let labels = self.labels(source_database)?;
        let mut db = Output::create(database.as_ref(), append)?;
        let mut idx = Output::create(index.as_ref(), append)?;
        let mut append = append;

        const AC_FLAG: u8 = 1;
        const SEQUENCE_FLAG: u8 = 4;
        let species_flag: u8 = if species.is_empty() { 0 } else { 2 };
        let complete = AC_FLAG | species_flag | SEQUENCE_FLAG;
        let mut record_flags: u8 = 0;

        let ac_label = labels.ac_label.as_bytes();
        let species_label = labels.species_label.as_bytes();
        let start_label = labels.sequence_start_label.as_bytes();
        let end_label = labels.sequence_end_label.as_bytes();
        let comment_label = labels.comment_label.as_bytes();

        // The start of the current protein in the source, and the position
        // after the last line read (the source's tellg() before a getline).
        let mut source_position: u64 = 0;
        let mut position_buffer: u64 = 0;
        let mut offset: u64 = 0;
        let mut sequence: Vec<u8> = Vec::new();
        let mut record = [0u8; RECORD_LENGTH];
        let mut raw = Vec::new();

        loop {
            let consumed = read_line(&mut source, &mut raw)?;
            if consumed == 0 {
                break;
            }
            offset += consumed as u64;
            // std::getline at the end of a file without a final line feed
            // leaves tellg() at -1.
            let next_position = if raw.len() == consumed {
                u64::MAX
            } else {
                offset
            };
            let line = clean_bytes_signed(&raw).to_vec();

            if line.is_empty() || line.starts_with(comment_label) {
                position_buffer = next_position;
                continue;
            }

            if record_flags == complete {
                if !line.starts_with(end_label) {
                    sequence.extend(line.iter().filter(|&&b| b != TRIE_DELIMITER));
                } else {
                    if !sequence.is_empty() {
                        write_trie_record(
                            &mut record,
                            &sequence,
                            source_position,
                            &mut append,
                            &mut db,
                            &mut idx,
                        )?;
                        source_position = position_buffer;
                    }
                    sequence.clear();
                    record_flags = 0;
                }
            }

            if record_flags & SEQUENCE_FLAG == 0 {
                if line.starts_with(ac_label) {
                    let mut pos = ac_label.len();
                    while line.len() > pos && (line[pos] as i8) < 33 {
                        pos += 1;
                    }
                    if pos != line.len() {
                        let name = &mut record[DB_POS_LENGTH + TRIE_DB_POS_LENGTH..];
                        name.fill(0);
                        let end = (pos + PROTEIN_NAME_LENGTH).min(line.len());
                        for (target, &byte) in name.iter_mut().zip(&line[pos..end]) {
                            *target = if byte == b'>' { b'}' } else { byte };
                        }
                        record_flags |= AC_FLAG;
                    } else {
                        record_flags = 0;
                    }
                }
                if species_flag != 0 && line.starts_with(species_label) && record_flags == AC_FLAG {
                    let rest = &line[species_label.len()..];
                    if contains(rest, species.as_bytes()) {
                        record_flags |= species_flag;
                    } else {
                        record_flags = 0;
                    }
                }
                if line.starts_with(start_label)
                    && (record_flags & (AC_FLAG | species_flag)) == (AC_FLAG | species_flag)
                {
                    record_flags |= SEQUENCE_FLAG;
                }
            }
            position_buffer = next_position;
        }

        // The last record has no end label in a FASTA file.
        if record_flags == complete && !sequence.is_empty() {
            write_trie_record(
                &mut record,
                &sequence,
                source_position,
                &mut append,
                &mut db,
                &mut idx,
            )?;
        }
        db.finish()?;
        idx.finish()
    }

    /// The accession and accession type of a protein description line
    /// (source `getACAndACType`), e.g. `DQ660451.1` and `GenBank` for
    /// `gi|110174602|gb|DQ660451.1|`.
    ///
    /// Recognised: SwissProt/TrEMBL (`sp|`, `tr|`), NCBI `gi|` lines with a
    /// `gb`, `emb`, `dbj`, `ref`, `sp`, `tr` or `gnl` database (otherwise a
    /// SwissProt accession in parentheses, else the gi number with type
    /// `gi`), `ref|`, `gnl|<type>|` and `lcl|`, and otherwise a six-character
    /// SwissProt-like accession (first character one of `JLOPQUX`) in
    /// parentheses, before the first `|` or space, or at the start. Anything
    /// else is the trimmed line with type `unknown`.
    pub fn ac_and_ac_type(&self, line: &str) -> (String, String) {
        const SWISSPROT_PREFIXES: &[u8] = b"JLOPQUX";
        let swissprot_like = |accession: &[u8]| {
            accession
                .first()
                .is_some_and(|c| SWISSPROT_PREFIXES.contains(c))
        };
        let mut line = line.as_bytes().to_vec();
        if line.first() == Some(&b'>') {
            line.remove(0);
        }
        if line.last().is_some_and(|&last| (last as i8) < 33) {
            line.pop();
        }
        let mut line = trim_bytes(&line).to_vec();
        let mut accession: Vec<u8> = Vec::new();
        let mut accession_type: Vec<u8> = Vec::new();

        // A SwissProt accession in parentheses, e.g. `(P68509)`.
        let parenthesised = |line: &[u8], accession: &mut Vec<u8>, accession_type: &mut Vec<u8>| {
            let mut pos1 = find(line, b'(', 0);
            if pos1 != NPOS {
                pos1 += 1;
                let pos2 = find(line, b')', pos1);
                if pos2 != NPOS {
                    *accession = sub(line, pos1, pos2 - pos1);
                    if accession.len() == 6 && swissprot_like(accession) {
                        *accession_type = b"SwissProt".to_vec();
                    } else {
                        accession.clear();
                    }
                }
            }
        };

        if line.starts_with(b"tr") || line.starts_with(b"sp") {
            accession = sub(&line, 3, find(&line, b'|', 3).wrapping_sub(3));
            accession_type = b"SwissProt".to_vec();
        } else if line.starts_with(b"gi") {
            let mut snd = find(&line, b'|', 3);
            let mut third = 0usize;
            if snd != NPOS {
                snd += 1;
                third = find(&line, b'|', snd).wrapping_add(1);
                accession = sub(&line, third, find(&line, b'|', third).wrapping_sub(third));
                accession_type = sub(&line, snd, third.wrapping_sub(1).wrapping_sub(snd));
            }
            match accession_type.as_slice() {
                b"gb" => accession_type = b"GenBank".to_vec(),
                b"emb" => accession_type = b"EMBL".to_vec(),
                b"dbj" => accession_type = b"DDBJ".to_vec(),
                b"ref" => accession_type = b"NCBI".to_vec(),
                b"sp" | b"tr" => accession_type = b"SwissProt".to_vec(),
                b"gnl" => {
                    accession_type = accession.clone();
                    snd = find(&line, b'|', third).wrapping_add(1);
                    third = find(&line, b'|', snd);
                    if third != NPOS {
                        accession = sub(&line, snd, third.wrapping_sub(snd));
                    } else {
                        third = find(&line, b' ', snd);
                        if third != NPOS {
                            accession = sub(&line, snd, third.wrapping_sub(snd));
                        } else {
                            accession = sub(&line, snd, NPOS);
                        }
                    }
                }
                _ => {
                    parenthesised(&line, &mut accession, &mut accession_type);
                    if accession.is_empty() {
                        accession_type = b"gi".to_vec();
                        if snd != NPOS {
                            accession = sub(&line, 3, snd.wrapping_sub(4));
                        } else {
                            snd = find(&line, b' ', 3);
                            if snd != NPOS {
                                accession = sub(&line, 3, snd - 3);
                            } else {
                                accession = sub(&line, 3, NPOS);
                            }
                        }
                    }
                }
            }
        } else if line.starts_with(b"ref") {
            accession = sub(&line, 4, find(&line, b'|', 4).wrapping_sub(4));
            accession_type = b"NCBI".to_vec();
        } else if line.starts_with(b"gnl") {
            line.drain(..3);
            accession_type = sub(&line, 0, find(&line, b'|', 0));
            accession = sub(&line, accession_type.len() + 1, NPOS);
        } else if line.starts_with(b"lcl") {
            line.drain(..4);
            accession_type = b"lcl".to_vec();
            accession = line.clone();
        } else {
            parenthesised(&line, &mut accession, &mut accession_type);
            if accession.is_empty() {
                accession = sub(&line, 0, find(&line, b'|', 0));
                if accession.len() == 6 && swissprot_like(&accession) {
                    accession_type = b"SwissProt".to_vec();
                } else {
                    accession = sub(&line, 0, find(&line, b' ', 0));
                    if accession.len() == 6 && swissprot_like(&accession) {
                        accession_type = b"SwissProt".to_vec();
                    } else {
                        accession = sub(&line, 0, 6);
                        if swissprot_like(&accession) {
                            accession_type = b"SwissProt".to_vec();
                        } else {
                            accession.clear();
                        }
                    }
                }
            }
        }
        if accession.is_empty() {
            accession = trim_bytes(&line).to_vec();
            accession_type = b"unknown".to_vec();
        }
        (
            String::from_utf8_lossy(&accession).into_owned(),
            String::from_utf8_lossy(&accession_type).into_owned(),
        )
    }

    /// Set the precursor m/z and retention time of the identifications from
    /// their spectrum files (source `getPrecursorRTandMZ`): for each file,
    /// identification `ids[i]` of a pair `(i, scan)` gets the first precursor
    /// m/z and the retention time of the file's `scan`-th spectrum (one
    /// based, in file order). A file without pairs is skipped.
    ///
    /// # Errors
    ///
    /// The errors of [`InspectOutfile::experiment`]; [`Error::Parse`]
    /// "`file` in: Not enough scans in file! (`n` available, should be at
    /// least `m`)" when the file's last pair names a scan beyond it; and,
    /// where the source's behaviour is undefined, [`Error::Parse`] for a scan
    /// 0 or beyond the file, or a spectrum without precursor, and
    /// [`Error::InvalidValue`] for an index beyond `ids`.
    pub fn precursor_rt_and_mz(
        &self,
        files_and_scans: &[(String, Vec<(usize, usize)>)],
        ids: &mut [PeptideIdentification],
    ) -> Result<()> {
        for (file, scans) in files_and_scans {
            let Some(&(_, last_scan)) = scans.last() else {
                continue;
            };
            let (experiment, _) = self.experiment(file)?;
            let available = experiment.spectra.len();
            if available < last_scan {
                return Err(parse_error(
                    file,
                    format!(
                        "Not enough scans in file! ({available} available, should be at least {last_scan})"
                    ),
                ));
            }
            for &(id_index, scan) in scans {
                let spectrum = scan
                    .checked_sub(1)
                    .and_then(|i| experiment.spectra.get(i))
                    .ok_or_else(|| {
                        parse_error(
                            file,
                            format!("scan {scan} is not among the {available} scans"),
                        )
                    })?;
                let precursor = spectrum
                    .precursors
                    .first()
                    .ok_or_else(|| parse_error(file, format!("scan {scan} has no precursor")))?;
                let given = ids.len();
                let id = ids.get_mut(id_index).ok_or_else(|| {
                    Error::InvalidValue(format!(
                        "peptide identification {id_index} does not exist ({given} given)"
                    ))
                })?;
                id.mz = Some(precursor.mz);
                id.rt = Some(spectrum.rt);
            }
        }
        Ok(())
    }

    /// The line labels of a FASTA (`>`) or SwissProt (`SQ`) database, from
    /// its first non-empty line that starts with `>` or `SQ` (source
    /// `getLabels`): accession, sequence start, sequence end, comment and
    /// species labels `>`, `>`, `>`, `;`, `>` for FASTA and `AC`, `SQ`,
    /// `//`, `CC`, `OS` for SwissProt.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file is missing ("the file '…' could not be
    /// found") or unreadable; [`Error::Parse`] "`file` in: database has
    /// unknown file format (neither trie nor FASTA nor swissprot)" when no
    /// line identifies the format.
    pub fn labels(&self, source_database: impl AsRef<Path>) -> Result<DatabaseLabels> {
        let source_database = source_database.as_ref();
        let mut reader = BufReader::new(open(source_database)?);
        let mut raw = Vec::new();
        let label = |text: &str| text.to_owned();
        while read_line(&mut reader, &mut raw)? != 0 {
            let line = clean_bytes_signed(&raw);
            if line.is_empty() {
                continue;
            }
            if line.starts_with(b">") {
                return Ok(DatabaseLabels {
                    ac_label: label(">"),
                    sequence_start_label: label(">"),
                    sequence_end_label: label(">"),
                    comment_label: label(";"),
                    species_label: label(">"),
                });
            }
            if line.starts_with(b"SQ") {
                return Ok(DatabaseLabels {
                    ac_label: label("AC"),
                    sequence_start_label: label("SQ"),
                    sequence_end_label: label("//"),
                    comment_label: label("CC"),
                    species_label: label("OS"),
                });
            }
        }
        Err(parse_error(
            source_database.display(),
            "database has unknown file format (neither trie nor FASTA nor swissprot)",
        ))
    }

    /// Append the sequences of the `wanted_records` (the keys, in ascending
    /// order; the values are ignored) of a trie database to `sequences`, and
    /// return the record numbers whose sequence is empty or beyond the
    /// database (source `getSequences`).
    ///
    /// As the source's stream does after reading nothing, every record after
    /// an empty or missing one is empty too.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the database is missing ("the file '…' could not
    /// be found") or unreadable, and [`Error::InvalidValue`] for a sequence
    /// longer than [`MAX_RECORD_BYTES`].
    pub fn sequences(
        &self,
        database: impl AsRef<Path>,
        wanted_records: &BTreeMap<usize, usize>,
        sequences: &mut Vec<String>,
    ) -> Result<Vec<usize>> {
        let mut reader = BufReader::new(open(database.as_ref())?);
        let mut not_found = Vec::new();
        let mut seen = 0usize;
        let mut failed = false;
        for &wanted in wanted_records.keys() {
            let mut sequence = Vec::new();
            if !failed {
                while seen < wanted {
                    if !skip_past_delimiter(&mut reader, TRIE_DELIMITER)? {
                        break;
                    }
                    seen += 1;
                }
                if seen == wanted {
                    sequence = read_until_delimiter(&mut reader, TRIE_DELIMITER)?;
                }
                failed = sequence.is_empty();
            }
            seen = seen.max(wanted);
            if sequence.is_empty() {
                not_found.push(wanted);
            }
            sequences.push(String::from_utf8_lossy(&sequence).into_owned());
        }
        Ok(not_found)
    }

    /// Load a spectrum file whose type is recognised from its content, and
    /// return it with the type's name, e.g. `mzXML` (source
    /// `getExperiment`).
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file is missing or unreadable; [`Error::Parse`]
    /// "`file` in: Could not determine type of the file. Aborting!";
    /// [`Error::Parse`] "type is not supported for loading experiments in:
    /// `file`" for a recognised type the source's `FileHandler::loadExperiment`
    /// has no case for, such as FASTA (the pinned Release build throws exactly
    /// this `ParseError`); and the reader's errors ([`Error::Unsupported`] for
    /// a type the source loads but the crate has no native reader for).
    pub fn experiment(&self, path: impl AsRef<Path>) -> Result<(MSExperiment, String)> {
        let path = path.as_ref();
        let mut preview = Vec::new();
        open(path)?.take(65_536).read_to_end(&mut preview)?;
        let kind = type_by_content(&preview);
        let experiment = match kind {
            FileType::Unknown => {
                return Err(parse_error(
                    path.display(),
                    "Could not determine type of the file. Aborting!",
                ));
            }
            FileType::MzXml => crate::format::mzxml::load(path)?,
            FileType::MzData => crate::format::mzdata::load(path)?,
            // The cases of the source's `FileHandler::loadExperiment` switch.
            other @ (FileType::Dta
            | FileType::Dta2d
            | FileType::MzMl
            | FileType::ImzMl
            | FileType::Mgf
            | FileType::Ms2
            | FileType::SqMass
            | FileType::Xmass
            | FileType::Msp
            | FileType::BrukerTdf
            | FileType::Raw) => FileHandler::load_experiment_with_options(
                path,
                &[other],
                &PeakFileOptions::default(),
            )?,
            // The switch's `default`.
            _ => {
                return Err(parse_error(
                    "type is not supported for loading experiments",
                    path.display(),
                ));
            }
        };
        Ok((experiment, kind.name().to_owned()))
    }

    /// Set the search engine `InsPecT` and its version from the output of
    /// the Inspect executable run without parameters, which names it as
    /// `InsPecT version 20100331` (or, in older versions, `InsPecT vesrion`).
    /// Returns whether a version was found; otherwise the version is
    /// `unknown` (source `getSearchEngineAndVersion`).
    pub fn search_engine_and_version(
        &self,
        cmd_output: &str,
        protein_identification: &mut ProteinIdentification,
    ) -> bool {
        protein_identification.search_engine = "InsPecT".into();
        protein_identification.search_engine_version = "unknown".into();
        for (at, _) in cmd_output.match_indices("InsPecT ") {
            let rest = &cmd_output[at + "InsPecT ".len()..];
            let Some(rest) = rest
                .strip_prefix("version ")
                .or_else(|| rest.strip_prefix("vesrion "))
            else {
                continue;
            };
            let digits: String = rest.chars().take_while(char::is_ascii_digit).collect();
            if !digits.is_empty() {
                protein_identification.search_engine_version = digits;
                return true;
            }
        }
        false
    }

    /// Find the columns of an Inspect result header line (tab separated):
    /// `#SpectrumFile`, `Scan#`, `Annotation`, `Protein`, `Charge`,
    /// `MQScore`, `p-value`, `RecordNumber`, `DBFilePos` and `SpecFilePos`
    /// (source `readOutHeader`). A repeated name takes its last column.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] "`filename` in: at least one of the columns
    /// '#SpectrumFile', 'Scan#', 'Annotation', 'Protein', 'Charge',
    /// 'MQScore', 'p-value', 'RecordNumber', 'DBFilePos' or 'SpecFilePos' is
    /// missing!".
    pub fn read_out_header(&self, filename: &str, header_line: &str) -> Result<InspectOutHeader> {
        let fields = split(header_line, '\t');
        let mut columns: [Option<usize>; 10] = [None; 10];
        const NAMES: [&str; 10] = [
            "#SpectrumFile",
            "Scan#",
            "Annotation",
            "Protein",
            "Charge",
            "MQScore",
            "p-value",
            "RecordNumber",
            "DBFilePos",
            "SpecFilePos",
        ];
        for (position, field) in fields.iter().enumerate() {
            if let Some(slot) = NAMES.iter().position(|name| name == field) {
                columns[slot] = Some(position);
            }
        }
        let missing = || {
            parse_error(
                filename,
                "at least one of the columns '#SpectrumFile', 'Scan#', 'Annotation', 'Protein', 'Charge', 'MQScore', 'p-value', 'RecordNumber', 'DBFilePos' or 'SpecFilePos' is missing!",
            )
        };
        let mut values = [0usize; 10];
        for (value, column) in values.iter_mut().zip(columns) {
            *value = column.ok_or_else(missing)?;
        }
        Ok(InspectOutHeader {
            spectrum_file_column: values[0],
            scan_column: values[1],
            peptide_column: values[2],
            protein_column: values[3],
            charge_column: values[4],
            mq_score_column: values[5],
            p_value_column: values[6],
            record_number_column: values[7],
            db_file_pos_column: values[8],
            spec_file_pos_column: values[9],
            number_of_columns: fields.len(),
        })
    }
}

/// The source's line clean-up on `char` (signed on x86): drop one final byte
/// below 33 or above 127, then trim space, tab, line feed and carriage return.
fn clean_bytes_signed(line: &[u8]) -> &[u8] {
    let line = match line.split_last() {
        Some((&last, rest)) if (last as i8) < 33 => rest,
        _ => line,
    };
    trim_bytes(line)
}

/// Write one trie record: the sequence (after a `*` unless it is the first
/// of a new database) and its index record.
fn write_trie_record(
    record: &mut [u8; RECORD_LENGTH],
    sequence: &[u8],
    source_position: u64,
    append: &mut bool,
    db: &mut Output,
    idx: &mut Output,
) -> Result<()> {
    if *append {
        db.write(&[TRIE_DELIMITER])?;
    }
    record[..DB_POS_LENGTH].copy_from_slice(&source_position.to_le_bytes());
    record[DB_POS_LENGTH..DB_POS_LENGTH + TRIE_DB_POS_LENGTH]
        .copy_from_slice(&trie_position_bytes(db.position));
    idx.write(record)?;
    db.write(sequence)?;
    *append = true;
    Ok(())
}

/// `std::string::find(c, from)` on bytes.
fn find(text: &[u8], byte: u8, from: usize) -> usize {
    text.get(from..)
        .and_then(|rest| rest.iter().position(|&b| b == byte))
        .map_or(NPOS, |pos| pos + from)
}

/// `StringUtils::substr(s, pos, n)` on bytes: `pos` clamped, at most `n`.
fn sub(text: &[u8], pos: usize, n: usize) -> Vec<u8> {
    let start = pos.min(text.len());
    let end = start.saturating_add(n).min(text.len());
    text[start..end].to_vec()
}

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    needle.is_empty()
        || haystack
            .windows(needle.len())
            .any(|window| window == needle)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn byte_helpers_follow_std_string() {
        assert_eq!(find(b"a|b|c", b'|', 2), 3);
        assert_eq!(find(b"a|b", b'|', 9), NPOS);
        assert_eq!(sub(b"abcdef", 2, NPOS), b"cdef");
        assert_eq!(sub(b"abc", 9, 1), b"");
        assert_eq!(clean_bytes_signed(b" x \r"), b"x");
        assert_eq!(trie_position_bytes(0x1_0000_00f6), [0xf6, 0, 0, 0]);
    }
}
