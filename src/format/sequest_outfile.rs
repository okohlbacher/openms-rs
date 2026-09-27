// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Reader for Sequest `.out` result files.
//!
//! Ports `FORMAT/SequestOutfile.h`/`.cpp`. A Sequest `.out` file holds one
//! spectrum's search: a free-text header (the `.out` name, the Sequest version,
//! the run date, the precursor `(M+H)+` mass and charge, `display top N/..`),
//! a table header whose column names locate the fields, a `---` rule, and one
//! line per displayed peptide, each optionally followed by lines naming further
//! proteins. See `docs/SEARCH_ENGINE_FILES_SUPPORT.md` for the API mapping and
//! the divergences.
//!
//! The source header carries two `@todo`s, handling modifications and
//! rewriting the parser; both are still true and are recorded in the support
//! document.

use crate::chemistry::AASequence;
use crate::data_structures::DateTime;
use crate::data_structures::string_utils::{split, to_double, to_float, to_int32, trim};
use crate::format::search_engine_text::{
    MAX_TEXT_BYTES, NPOS, clean_line, find_byte, getline_lines, read_text, rfind_byte, substr,
};
use crate::identification::{
    FlankingResidue, PeptideEvidence, PeptideHit, PeptideIdentification, ProteinHit,
    ProteinIdentification,
};
use crate::metadata::MetaValue;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// Score type when the table has an `Sf` column.
pub const SCORE_TYPE: &str = "SEQUEST";

/// Score type when it has none, and `Sp` is the score.
pub const PRELIMINARY_SCORE_TYPE: &str = "SEQUEST prelim.";

/// First characters of a six-character SwissProt accession.
const SWISSPROT_PREFIXES: &str = "JLOPQUX";

/// What [`SequestOutfile::read_out_header`] finds before the result table.
///
/// Column positions are `None` where the source's are -1: the column is
/// absent from the table header.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SequestOutHeader {
    /// Run date and time, from the `MM/DD/YYYY, hh:mm AM|PM` line.
    pub datetime: DateTime,
    /// Precursor `(M+H)+` mass, read as a `float` as the source does.
    pub precursor_mz_value: f64,
    /// Precursor charge, from the `(+z)` after the mass.
    pub charge: i32,
    /// 1 when the precursor mass is monoisotopic (`MONO`), else 0.
    pub precursor_mass_type: usize,
    /// 1 when the fragment ion masses are monoisotopic (`MONO`), else 0.
    pub ion_mass_type: usize,
    /// The `N` of `display top N/..`: how many peptide lines follow at most.
    pub displayed_peptides: usize,
    /// Search engine name, e.g. `TurboSEQUEST`.
    pub sequest: String,
    /// Search engine version, e.g. `v.27 (rev. 12)`.
    pub sequest_version: String,
    /// `amino acids` or `bases`, from the `# amino acids`/`# bases` line.
    pub database_type: String,
    /// Column `#`.
    pub number_column: Option<usize>,
    /// Column `Rank/Sp`.
    pub rank_sp_column: Option<usize>,
    /// Column `Id#`.
    pub id_column: Option<usize>,
    /// Column `(M+H)+`.
    pub mh_column: Option<usize>,
    /// Column `deltCn`.
    pub delta_cn_column: Option<usize>,
    /// Column `XCorr`.
    pub xcorr_column: Option<usize>,
    /// Column `Sp`.
    pub sp_column: Option<usize>,
    /// Column `Sf`.
    pub sf_column: Option<usize>,
    /// Column `Ions`.
    pub ions_column: Option<usize>,
    /// Column `Reference`.
    pub reference_column: Option<usize>,
    /// Column `Peptide`.
    pub peptide_column: Option<usize>,
    /// The score column: `Sf` when present, else `Sp`.
    pub score_column: Option<usize>,
    /// Number of columns in the table header.
    pub number_of_columns: usize,
}

/// Sequences found by [`SequestOutfile::get_sequences`].
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SequenceLookup {
    /// The sequences found, in database order.
    pub sequences: Vec<String>,
    /// The accession and position (from the lookup map) of each sequence.
    pub found: Vec<(String, usize)>,
    /// The accessions not found, with their positions.
    pub not_found: BTreeMap<String, usize>,
}

/// Reader for Sequest `.out` files. It holds no state; the source's copy
/// constructor, assignment and `operator==` (always true) are the derived
/// ones.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SequestOutfile;

fn parse_error(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

fn at(line: usize) -> impl Fn(Error) -> Error {
    move |error| match error {
        Error::Parse { line: 0, message } => Error::Parse { line, message },
        other => other,
    }
}

fn index<'a>(fields: &'a [String], column: usize, line: usize) -> Result<&'a str> {
    fields
        .get(column)
        .map(String::as_str)
        .ok_or_else(|| Error::Parse {
            line,
            message: format!("column {column} is missing"),
        })
}

/// C `atoi` as glibc implements it (`(int) strtol(s, NULL, 10)`): leading
/// white space, an optional sign, then digits up to the first other
/// character; 0 when there are none. `strtol` saturates at the `long` range
/// and the cast then keeps the low 32 bits.
fn atoi(text: &str) -> i32 {
    let bytes = text.as_bytes();
    let mut cursor = 0;
    while bytes
        .get(cursor)
        .is_some_and(|b| matches!(b, b' ' | b'\t' | b'\n' | b'\x0b' | b'\x0c' | b'\r'))
    {
        cursor += 1;
    }
    let negative = match bytes.get(cursor) {
        Some(b'-') => {
            cursor += 1;
            true
        }
        Some(b'+') => {
            cursor += 1;
            false
        }
        _ => false,
    };
    let mut value: i64 = 0;
    while let Some(digit) = bytes.get(cursor).filter(|b| b.is_ascii_digit()) {
        let digit = i64::from(digit - b'0');
        value = if negative {
            value.saturating_mul(10).saturating_sub(digit)
        } else {
            value.saturating_mul(10).saturating_add(digit)
        };
        cursor += 1;
    }
    value as i32
}

/// A flanking residue from a Sequest peptide string. Sequest writes `-` for
/// a protein terminus, which the Rust evidence stores as the terminus marker
/// (`[` before, `]` after); the source stores the raw `-`.
fn flanking(code: u8, before: bool) -> Result<FlankingResidue> {
    match code {
        b'-' if before => Ok(FlankingResidue::NTerminus),
        b'-' => Ok(FlankingResidue::CTerminus),
        // The source reads the terminating NUL of a peptide string that ends
        // in its separator.
        0 => Ok(FlankingResidue::Unknown),
        other => FlankingResidue::from_code(char::from(other))
            .map_err(|_| parse_error(format!("invalid flanking residue '{}'", char::from(other)))),
    }
}

impl SequestOutfile {
    /// A Sequest `.out` reader; it holds no state.
    pub fn new() -> Self {
        Self
    }

    /// Load one Sequest `.out` file.
    ///
    /// Protein hits accumulate in `protein_identification` across files, as
    /// the source's TOPP caller loads a directory of `.out` files into one:
    /// an accession already among its hits is not added again. Its search
    /// engine, version, identifier (`<engine>_<run date>`) and date are
    /// replaced by this file's.
    ///
    /// `pvalues` are PeptideProphet p-values of the displayed peptides, in
    /// order; a peptide is kept when its p-value is at most
    /// `p_value_threshold`, and every peptide is kept when `pvalues` is empty.
    /// A peptide past the end of a non-empty `pvalues` is dropped (the source
    /// increments its iterator past the end, which is undefined behaviour).
    /// `ignore_proteins_per_peptide` is the source's workaround for files
    /// whose `Reference` carries a `+N` suffix without the `N` protein lines
    /// that should follow.
    ///
    /// The source's `database` parameter is only stored in a local list and
    /// never read, so the port does not take it.
    ///
    /// Returns the spectrum's identification, or `None` when no peptide was
    /// kept (the source then appends nothing).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `p_value_threshold` is outside `0..=1` or
    /// NaN (the source's `Exception::IllegalArgument`; it lets NaN through and
    /// then keeps no peptide). [`Error::Io`] when the file cannot be read, and
    /// the errors of [`SequestOutfile::read_out_header`]. [`Error::Parse`]
    /// with the source's message "Wrong number of columns in line N! (M
    /// present, should be K)" for a peptide line of the wrong width, where N
    /// counts as the source does, skipping the extra-protein lines; the
    /// error's `line` is the physical line. Also [`Error::Parse`] when a
    /// numeric field does not convert or is not finite, or a flanking residue
    /// is invalid. On error `protein_identification` is unchanged; the source
    /// has already replaced its engine, version and identifier.
    pub fn load(
        &self,
        path: impl AsRef<Path>,
        protein_identification: &mut ProteinIdentification,
        p_value_threshold: f64,
        pvalues: &[f64],
        ignore_proteins_per_peptide: bool,
    ) -> Result<Option<PeptideIdentification>> {
        if !(0.0..=1.0).contains(&p_value_threshold) {
            return Err(Error::InvalidValue(
                "the parameters 'p_value_threshold' must be >= 0 and <=1 !".into(),
            ));
        }
        let text = read_text(BufReader::new(File::open(path)?), MAX_TEXT_BYTES)?;
        let header = self.parse_out_header(&text)?;
        let lines = getline_lines(&text);
        let columns = |column: Option<usize>| {
            column.ok_or_else(|| parse_error("a required Sequest column is missing"))
        };
        let reference_column = columns(header.reference_column)?;
        let peptide_column = columns(header.peptide_column)?;
        let rank_sp_column = columns(header.rank_sp_column)?;
        let score_column = columns(header.score_column)?;

        let identifier = format!("{}_{}", header.sequest, header.datetime.date_string());
        let mut protein_hits = protein_identification.hits.clone();
        let mut ac_position_map: BTreeMap<String, usize> = BTreeMap::new();
        for hit in &protein_hits {
            let position = ac_position_map.len();
            ac_position_map
                .entry(hit.accession.clone())
                .or_insert(position);
        }

        // Skip everything up to the `---` rule.
        let mut cursor = 0usize;
        let mut line_number = 0usize;
        while cursor < lines.len() {
            let line = clean_line(lines[cursor]);
            cursor += 1;
            line_number += 1;
            if line.starts_with("---") {
                break;
            }
        }

        let score_type = if header.sf_column.is_none() {
            PRELIMINARY_SCORE_TYPE
        } else {
            SCORE_TYPE
        };
        let mut identification = PeptideIdentification {
            identifier: identifier.clone(),
            score_type: score_type.into(),
            significance_threshold: p_value_threshold,
            mz: Some(header.precursor_mz_value),
            ..PeptideIdentification::default()
        };

        let mut viewed_peptides = 0usize;
        let mut peptide_index = 0usize;
        while viewed_peptides < header.displayed_peptides {
            line_number += 1;
            let Some(raw) = lines.get(cursor) else {
                break;
            };
            cursor += 1;
            let physical = cursor;
            let line = clean_line(raw);
            if line.is_empty() {
                continue;
            }
            viewed_peptides += 1;
            let mut fields = self
                .get_columns(line, header.number_of_columns, reference_column)
                .unwrap_or_default();
            if fields.len() != header.number_of_columns {
                return Err(Error::Parse {
                    line: physical,
                    message: format!(
                        "Wrong number of columns in line {line_number}! ({} present, should be {})",
                        fields.len(),
                        header.number_of_columns
                    ),
                });
            }
            let reference = &mut fields[reference_column];
            let plus = rfind_byte(reference, b'+');
            let proteins_per_peptide = if plus == NPOS {
                0
            } else {
                let count = to_int32(&reference[plus..]).map_err(|message| Error::Parse {
                    line: physical,
                    message,
                })?;
                reference.truncate(plus);
                // The source stores the `Int` in a `Size`: a negative count
                // becomes about 2^64 lines to read, each past the end of the
                // file an empty one, which never finishes.
                usize::try_from(count).map_err(|_| Error::Parse {
                    line: physical,
                    message: format!("negative protein count {count} in the Reference column"),
                })?
            };

            let accepted = if pvalues.is_empty() {
                true
            } else {
                pvalues
                    .get(peptide_index)
                    .is_some_and(|p| *p <= p_value_threshold)
            };
            peptide_index += 1;
            if !accepted {
                if !ignore_proteins_per_peptide {
                    cursor = cursor.saturating_add(proteins_per_peptide);
                }
                continue;
            }

            let mut hit = self
                .peptide_hit(
                    &fields,
                    &header,
                    score_column,
                    peptide_column,
                    rank_sp_column,
                )
                .map_err(at(physical))?;
            let (accession, _) = self
                .get_ac_and_ac_type(&fields[reference_column])
                .map_err(at(physical))?;
            if !ac_position_map.contains_key(&accession) {
                ac_position_map.insert(accession.clone(), protein_hits.len());
                protein_hits.push(ProteinHit {
                    accession: accession.clone(),
                    ..ProteinHit::default()
                });
            }
            if let Some(evidence) = hit.evidences.first_mut() {
                evidence.protein_accession = accession;
            }
            if !ignore_proteins_per_peptide {
                // The source reads on without checking for the end of the
                // file, where `getline` leaves each line empty and an empty
                // accession is recorded per missing line.
                if proteins_per_peptide > lines.len().saturating_sub(cursor) {
                    return Err(Error::Parse {
                        line: physical,
                        message: format!(
                            "the Reference column announces {proteins_per_peptide} more proteins, but the file ends first"
                        ),
                    });
                }
                for _ in 0..proteins_per_peptide {
                    let raw = lines[cursor];
                    cursor += 1;
                    let line = clean_line(raw);
                    let rest = substr(line, 3, NPOS).map_err(at(cursor))?;
                    let (accession, _) = self.get_ac_and_ac_type(rest).map_err(at(cursor))?;
                    if !ac_position_map.contains_key(&accession) {
                        ac_position_map.insert(accession.clone(), protein_hits.len());
                        protein_hits.push(ProteinHit {
                            accession: accession.clone(),
                            ..ProteinHit::default()
                        });
                    }
                    hit.evidences.push(PeptideEvidence {
                        protein_accession: accession,
                        ..PeptideEvidence::default()
                    });
                }
            }
            identification.hits.push(hit);
        }

        let mut updated = protein_identification.clone();
        updated.search_engine = header.sequest.clone();
        updated.search_engine_version = header.sequest_version.clone();
        updated.identifier = identifier;
        updated.hits = protein_hits;
        updated.date_time = Some(header.datetime.get());
        *protein_identification = updated;
        Ok((!identification.hits.is_empty()).then_some(identification))
    }

    /// One peptide line's hit, with its first evidence's flanking residues.
    fn peptide_hit(
        &self,
        fields: &[String],
        header: &SequestOutHeader,
        score_column: usize,
        peptide_column: usize,
        rank_sp_column: usize,
    ) -> Result<PeptideHit> {
        let number = |column: usize| -> Result<f64> {
            let text = index(fields, column, 0)?;
            let value = to_double(text).map_err(parse_error)?;
            if value.is_finite() {
                Ok(value)
            } else {
                Err(parse_error(format!("non-finite value '{text}'")))
            }
        };
        let mut hit = PeptideHit {
            score: number(score_column)?,
            charge: header.charge,
            ..PeptideHit::default()
        };
        let meta = &mut hit.metadata;
        if let Some(column) = header.rank_sp_column {
            meta.insert("RankSp".into(), index(fields, column, 0)?.into());
        }
        if let Some(column) = header.id_column {
            meta.insert("SequestId".into(), atoi(index(fields, column, 0)?).into());
        }
        for (column, name) in [
            (header.mh_column, "MH"),
            (header.delta_cn_column, "DeltCn"),
            (header.xcorr_column, "XCorr"),
            (header.sp_column, "Sp"),
            (header.sf_column, "Sf"),
        ] {
            if let Some(column) = column {
                meta.insert(name.into(), MetaValue::try_from(number(column)?)?);
            }
        }
        if let Some(column) = header.ions_column {
            meta.insert("Ions".into(), index(fields, column, 0)?.into());
        }

        let with_mods = index(fields, peptide_column, 0)?;
        let bytes = with_mods.as_bytes();
        let start = find_byte(with_mods, b'.', 0).wrapping_add(1);
        let end = rfind_byte(with_mods, b'.');
        let mut evidence = PeptideEvidence::default();
        if start >= 2 {
            evidence.aa_before = flanking(bytes[start - 2], true)?;
        }
        if end < bytes.len().wrapping_add(1) {
            evidence.aa_after = flanking(bytes.get(end + 1).copied().unwrap_or(0), false)?;
        }
        // Keep the upper-case letters between the dots: modification symbols
        // such as `*` are dropped, as the source's `@todo` on modifications
        // says.
        let residues: String = substr(with_mods, start, end.wrapping_sub(start))?
            .chars()
            .filter(char::is_ascii_uppercase)
            .collect();
        hit.sequence = AASequence::parse(&residues)?;
        let rank_sp = index(fields, rank_sp_column, 0)?;
        let rank_text = substr(rank_sp, 0, find_byte(rank_sp, b'/', 0))?;
        hit.rank = to_int32(rank_text).map_err(parse_error)? as u32;
        hit.evidences.push(evidence);
        Ok(hit)
    }

    /// Split one peptide line into its columns (source `getColumns`).
    ///
    /// The line is split at spaces, empty fields dropped, and then fields
    /// merged: `a` `/` `b` into `a/b`, a field followed by one starting with
    /// `/` or ending with `/` with its successor, and a field followed by `+`
    /// and digits with that suffix. While more than `number_of_columns`
    /// fields remain, the field after `reference_column` joins it with a
    /// space, because protein names may contain spaces.
    ///
    /// Returns `None` for an empty line (the source's `false`).
    pub fn get_columns(
        &self,
        line: &str,
        number_of_columns: usize,
        reference_column: usize,
    ) -> Option<Vec<String>> {
        if line.is_empty() {
            return None;
        }
        let mut fields: Vec<String> = split(line, ' ')
            .into_iter()
            .filter(|field| !field.is_empty())
            .map(str::to_owned)
            .collect();
        let mut i = 0usize;
        while i < fields.len() {
            if i + 1 >= fields.len() {
                i += 1;
                continue;
            }
            let next = fields[i + 1].clone();
            if next == "/" && i + 2 < fields.len() {
                let third = fields.remove(i + 2);
                fields.remove(i + 1);
                fields[i].push_str(&next);
                fields[i].push_str(&third);
            } else if next.starts_with('/') || fields[i].ends_with('/') {
                fields.remove(i + 1);
                fields[i].push_str(&next);
            } else if let Some(digits) = next.strip_prefix('+') {
                if !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()) {
                    fields.remove(i + 1);
                    fields[i].push_str(&next);
                } else {
                    i += 1;
                }
            } else {
                i += 1;
            }
        }
        while fields.len() > number_of_columns && reference_column + 1 < fields.len() {
            let next = fields.remove(reference_column + 1);
            fields[reference_column].push(' ');
            fields[reference_column].push_str(&next);
        }
        Some(fields)
    }

    /// Look up the sequences of the proteins in `ac_position_map` in a FASTA
    /// database (source `getSequences`).
    ///
    /// Empty lines and lines starting with `;` are skipped; each `>` line's
    /// accession comes from [`SequestOutfile::get_ac_and_ac_type`]. Reading
    /// stops once every accession is found. The source appends to the
    /// caller's `sequences` and `found` and overwrites `not_found`; the port
    /// returns all three for this database, so a caller searching several
    /// databases appends them itself.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the database cannot be read (source
    /// `Exception::FileNotFound`), [`Error::Parse`] when it is not UTF-8 or an
    /// accession cut splits a character, [`Error::InvalidValue`] when it is
    /// larger than [`MAX_TEXT_BYTES`].
    pub fn get_sequences(
        &self,
        database: impl AsRef<Path>,
        ac_position_map: &BTreeMap<String, usize>,
    ) -> Result<SequenceLookup> {
        let text = read_text(BufReader::new(File::open(database)?), MAX_TEXT_BYTES)?;
        let mut lookup = SequenceLookup {
            not_found: ac_position_map.clone(),
            ..SequenceLookup::default()
        };
        let mut current: Option<String> = None;
        let mut sequence = String::new();
        for raw in getline_lines(&text) {
            if lookup.not_found.is_empty() {
                break;
            }
            let line = clean_line(raw);
            if line.is_empty() || line.starts_with(';') {
                continue;
            }
            if line.starts_with('>') {
                let (accession, _) = self.get_ac_and_ac_type(line)?;
                if let Some(previous) = current.take() {
                    if let Some(position) = lookup.not_found.remove(&previous) {
                        lookup.sequences.push(std::mem::take(&mut sequence));
                        lookup.found.push((previous, position));
                    }
                }
                current = lookup
                    .not_found
                    .contains_key(&accession)
                    .then_some(accession);
                sequence.clear();
            } else if current.is_some() {
                sequence.push_str(line);
            }
        }
        if let Some(previous) = current {
            if let Some(position) = lookup.not_found.remove(&previous) {
                lookup.sequences.push(sequence);
                lookup.found.push((previous, position));
            }
        }
        Ok(lookup)
    }

    /// The accession and accession type of a protein description (source
    /// `getACAndACType`), e.g. `AAD44166.1` and `GenBank` from the FASTA line
    /// `>gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]`.
    ///
    /// Recognised forms: `sp|..|`/`tr|..|` (SwissProt), `gi|number|type|ac|`
    /// with the types `gb` (GenBank), `emb` (EMBL), `dbj` (DDBJ), `ref`
    /// (NCBI), `sp`/`tr` (SwissProt) and `gnl|db|ac` (the database as type),
    /// otherwise a parenthesised six-character SwissProt accession or the gi
    /// number; `ref|..|` (NCBI), `gnl|db|ac`, `lcl|ac`; and else a
    /// six-character SwissProt accession at the start, before `|` or space,
    /// or in parentheses. Anything unrecognised is its own trimmed accession
    /// of type `unknown`. SwissProt accessions are recognised by their first
    /// letter, one of `JLOPQUX`.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] when a byte-offset cut would split a multi-byte
    /// character; the source cuts bytes.
    pub fn get_ac_and_ac_type(&self, line: &str) -> Result<(String, String)> {
        let mut accession = String::new();
        let mut accession_type = String::new();
        let line = line.strip_prefix('>').unwrap_or(line);
        let line = clean_line(line);
        let swissprot = |text: &str| {
            text.as_bytes()
                .first()
                .is_some_and(|b| SWISSPROT_PREFIXES.as_bytes().contains(b))
        };
        // The source's parenthesis rule: with no `(...)` the accession is
        // left as it is; a six-character SwissProt accession inside makes it
        // that, of type SwissProt; anything else inside clears it.
        let parenthesised =
            |line: &str, accession: &mut String, accession_type: &mut String| -> Result<()> {
                let open = find_byte(line, b'(', 0);
                if open == NPOS {
                    return Ok(());
                }
                let close = find_byte(line, b')', open + 1);
                if close == NPOS {
                    return Ok(());
                }
                let candidate = substr(line, open + 1, close - (open + 1))?;
                if candidate.len() == 6 && swissprot(candidate) {
                    *accession = candidate.to_owned();
                    *accession_type = "SwissProt".into();
                } else {
                    accession.clear();
                }
                Ok(())
            };

        if line.starts_with("tr") || line.starts_with("sp") {
            let end = find_byte(line, b'|', 3).wrapping_sub(3);
            accession = substr(line, 3, end)?.to_owned();
            accession_type = "SwissProt".into();
        } else if line.starts_with("gi") {
            let mut snd = find_byte(line, b'|', 3);
            let mut third = 0usize;
            if snd != NPOS {
                snd += 1;
                third = find_byte(line, b'|', snd).wrapping_add(1);
                let length = find_byte(line, b'|', third).wrapping_sub(third);
                accession = substr(line, third, length)?.to_owned();
                accession_type =
                    substr(line, snd, third.wrapping_sub(1).wrapping_sub(snd))?.to_owned();
            }
            match accession_type.as_str() {
                "gb" => accession_type = "GenBank".into(),
                "emb" => accession_type = "EMBL".into(),
                "dbj" => accession_type = "DDBJ".into(),
                "ref" => accession_type = "NCBI".into(),
                "sp" | "tr" => accession_type = "SwissProt".into(),
                "gnl" => {
                    accession_type = accession.clone();
                    let start = find_byte(line, b'|', third).wrapping_add(1);
                    let bar = find_byte(line, b'|', start);
                    accession = if bar != NPOS {
                        substr(line, start, bar - start)?.to_owned()
                    } else {
                        let space = find_byte(line, b' ', start);
                        if space != NPOS {
                            substr(line, start, space - start)?.to_owned()
                        } else {
                            substr(line, start, NPOS)?.to_owned()
                        }
                    };
                }
                _ => {
                    parenthesised(line, &mut accession, &mut accession_type)?;
                    if accession.is_empty() {
                        accession_type = "gi".into();
                        if snd != NPOS {
                            accession = substr(line, 3, snd.wrapping_sub(4))?.to_owned();
                        } else {
                            let space = find_byte(line, b' ', 3);
                            accession = if space != NPOS {
                                substr(line, 3, space - 3)?.to_owned()
                            } else {
                                substr(line, 3, NPOS)?.to_owned()
                            };
                        }
                    }
                }
            }
        } else if line.starts_with("ref") {
            let end = find_byte(line, b'|', 4).wrapping_sub(4);
            accession = substr(line, 4, end)?.to_owned();
            accession_type = "NCBI".into();
        } else if let Some(rest) = line.strip_prefix("gnl") {
            accession_type = substr(rest, 0, find_byte(rest, b'|', 0))?.to_owned();
            accession = substr(rest, accession_type.len() + 1, NPOS)?.to_owned();
        } else if line.starts_with("lcl") {
            accession_type = "lcl".into();
            accession = substr(line, 4, NPOS)?.to_owned();
        } else {
            parenthesised(line, &mut accession, &mut accession_type)?;
            if accession.is_empty() {
                let candidate = substr(line, 0, find_byte(line, b'|', 0))?;
                if candidate.len() == 6 && swissprot(candidate) {
                    accession = candidate.to_owned();
                    accession_type = "SwissProt".into();
                } else {
                    let candidate = substr(line, 0, find_byte(line, b' ', 0))?;
                    if candidate.len() == 6 && swissprot(candidate) {
                        accession = candidate.to_owned();
                        accession_type = "SwissProt".into();
                    } else {
                        let candidate = substr(line, 0, 6)?;
                        if swissprot(candidate) {
                            accession = candidate.to_owned();
                            accession_type = "SwissProt".into();
                        }
                    }
                }
            }
        }
        if accession.is_empty() {
            accession = trim(line).to_owned();
            accession_type = "unknown".into();
        }
        Ok((accession, accession_type))
    }

    /// Read the header of a Sequest `.out` file (source `readOutHeader`).
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be read (source
    /// `Exception::FileNotFound`), and [`Error::Parse`] with the source's
    /// messages: "No Sequest version found!" (no `.out` line followed by a
    /// line naming SEQUEST, or no version after it), "No time found!", "No
    /// precursor mass found found!", "No charge found!", "Not enough lines in
    /// file (no table header found)!", and "at least one of the columns '#',
    /// 'Rank/Sp', 'Id#', '(M+H)+', 'deltCn', 'XCorr', 'Sp', 'Ions',
    /// 'Reference' or 'Peptide' is missing!" (the source does not in fact
    /// require `Id#`). Also [`Error::Parse`] for a date, time, mass or charge
    /// that does not convert, or a `(M+H)+ mass` line too short to hold them,
    /// where the source reads past the end of its field list.
    pub fn read_out_header(&self, path: impl AsRef<Path>) -> Result<SequestOutHeader> {
        let text = read_text(BufReader::new(File::open(path)?), MAX_TEXT_BYTES)?;
        self.parse_out_header(&text)
    }

    fn parse_out_header(&self, text: &str) -> Result<SequestOutHeader> {
        let lines = getline_lines(text);
        let mut header = SequestOutHeader::default();
        let mut cursor = 0usize;
        let mut table_line = String::new();
        let mut found_table = false;
        while cursor < lines.len() {
            let line = clean_line(lines[cursor]);
            cursor += 1;
            let number = cursor;
            let fields = split(line, ',');
            if line.ends_with(".out") {
                let Some(raw) = lines.get(cursor) else {
                    return Err(parse_error("No Sequest version found!"));
                };
                cursor += 1;
                let mut version_line = clean_line(raw);
                if let Some(comma) = version_line.find(',') {
                    version_line = &version_line[..comma];
                }
                let upper = version_line.to_ascii_uppercase();
                let Some(found) = upper.find("SEQUEST") else {
                    return Err(Error::Parse {
                        line: cursor,
                        message: "No Sequest version found!".into(),
                    });
                };
                let pos1 = found + "SEQUEST".len();
                let pos = find_byte(version_line, b' ', 0);
                header.sequest = if pos == NPOS || pos >= pos1 {
                    substr(version_line, 0, pos1)?.to_owned()
                } else {
                    substr(version_line, pos, pos1 - pos)?.to_owned()
                };
                let space = find_byte(version_line, b' ', pos1);
                if space != NPOS {
                    let from = space + 1;
                    let comma = find_byte(version_line, b',', from);
                    header.sequest_version = if comma == NPOS {
                        substr(version_line, from, NPOS)?.to_owned()
                    } else {
                        substr(version_line, from, comma - from)?.to_owned()
                    };
                }
            } else if let Some(rest) = line.strip_prefix("(M+H)+ mass = ") {
                let parts = split(rest, ' ');
                let located = |message: String| Error::Parse {
                    line: number,
                    message,
                };
                let (Some(mass), Some(charge), Some(types)) =
                    (parts.first(), parts.get(3), parts.last())
                else {
                    return Err(located(format!("incomplete precursor line '{line}'")));
                };
                header.precursor_mz_value = f64::from(to_float(mass).map_err(located)?);
                header.charge =
                    to_int32(substr(charge, 1, 2).map_err(at(number))?).map_err(located)?;
                let mass_types = split(types, '/');
                let (Some(precursor), Some(ion)) = (mass_types.first(), mass_types.get(1)) else {
                    return Err(located(format!("no mass types in '{line}'")));
                };
                header.precursor_mass_type = usize::from(precursor.eq_ignore_ascii_case("MONO"));
                header.ion_mass_type = usize::from(ion.eq_ignore_ascii_case("MONO"));
            } else if fields.first().is_some_and(|f| {
                f.len() == 10 && f.as_bytes()[2] == b'/' && f.as_bytes()[5] == b'/'
            }) {
                let located = |error: Error| match error {
                    Error::Parse { message, .. } => Error::Parse {
                        line: number,
                        message,
                    },
                    other => other,
                };
                header.datetime.set_date(fields[0]).map_err(located)?;
                let clock = trim(fields.get(1).copied().unwrap_or(""));
                let parts = split(clock, ' ');
                let mut time = parts.first().copied().unwrap_or("").to_owned();
                let meridiem = parts.get(1).copied().unwrap_or("");
                let hour =
                    to_int32(substr(&time, 0, 2).map_err(at(number))?).map_err(|message| {
                        Error::Parse {
                            line: number,
                            message,
                        }
                    })?;
                if hour == 12 && meridiem == "AM" {
                    time.replace_range(..2.min(time.len()), "00");
                } else if hour != 12 && meridiem == "PM" {
                    time.replace_range(..2.min(time.len()), &(hour + 12).to_string());
                }
                time.push_str(":00");
                header.datetime.set_time(&time).map_err(located)?;
            } else if line.starts_with("# bases") {
                header.database_type = "bases".into();
            } else if line.starts_with("# amino acids") {
                header.database_type = "amino acids".into();
            } else if line.starts_with("display top")
                && fields.first().is_some_and(|f| f.starts_with("display top"))
            {
                let field = fields[0];
                let from = "display top ".len();
                let count = substr(field, from, find_byte(field, b'/', from).wrapping_sub(from))
                    .map_err(at(number))?;
                header.displayed_peptides = to_int32(count).map_err(|message| Error::Parse {
                    line: number,
                    message,
                })? as u32 as usize;
            } else if line.starts_with('#') {
                table_line = line.to_owned();
                found_table = true;
                break;
            }
        }
        if header.datetime == DateTime::default() {
            return Err(parse_error("No time found!"));
        }
        if header.sequest.is_empty() || header.sequest_version.is_empty() {
            return Err(parse_error("No Sequest version found!"));
        }
        if header.precursor_mz_value == 0.0 {
            return Err(parse_error("No precursor mass found found!"));
        }
        if header.charge == 0 {
            return Err(parse_error("No charge found!"));
        }
        if !found_table {
            return Err(parse_error(
                "Not enough lines in file (no table header found)!",
            ));
        }
        let names: Vec<&str> = split(&table_line, ' ')
            .into_iter()
            .map(trim)
            .filter(|name| !name.is_empty())
            .collect();
        header.number_of_columns = names.len();
        for (position, name) in names.iter().enumerate() {
            let slot = match *name {
                "#" => &mut header.number_column,
                "Rank/Sp" => &mut header.rank_sp_column,
                "Id#" => &mut header.id_column,
                "(M+H)+" => &mut header.mh_column,
                "deltCn" => &mut header.delta_cn_column,
                "XCorr" => &mut header.xcorr_column,
                "Sp" => &mut header.sp_column,
                "Sf" => &mut header.sf_column,
                "Ions" => &mut header.ions_column,
                "Reference" => &mut header.reference_column,
                "Peptide" => &mut header.peptide_column,
                _ => continue,
            };
            *slot = Some(position);
        }
        if header.number_column.is_none()
            || header.rank_sp_column.is_none()
            || header.mh_column.is_none()
            || header.delta_cn_column.is_none()
            || header.xcorr_column.is_none()
            || header.sp_column.is_none()
            || header.ions_column.is_none()
            || header.reference_column.is_none()
            || header.peptide_column.is_none()
        {
            return Err(parse_error(
                "at least one of the columns '#', 'Rank/Sp', 'Id#', '(M+H)+', 'deltCn', 'XCorr', 'Sp', 'Ions', 'Reference' or 'Peptide' is missing!",
            ));
        }
        header.score_column = header.sf_column.or(header.sp_column);
        Ok(header)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atoi_follows_glibc() {
        assert_eq!(atoi("  12abc"), 12);
        assert_eq!(atoi("-7"), -7);
        assert_eq!(atoi("x"), 0);
        assert_eq!(atoi("4294967297"), 1);
    }
}
