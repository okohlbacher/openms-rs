// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Reader for PepNovo de novo sequencing output.
//!
//! Ports `FORMAT/PepNovoOutfile.h`/`.cpp`. A PepNovo output file starts with a
//! free-text header (the PepNovo version, the fragment and precursor
//! tolerances, the PTMs considered), then holds one block per spectrum: a
//! `>> <file> <index> <title>` line, a tab-separated `#Index` table header and
//! one candidate sequence per line up to an empty line. See
//! `docs/SEARCH_ENGINE_FILES_SUPPORT.md`.

use crate::chemistry::ModificationsDB;
use crate::comparison::Tolerance;
use crate::concept::string_conversion::{split, to_double, to_float, to_int32, trim};
use crate::data_structures::DateTime;
use crate::format::search_engine_text::{
    MAX_TEXT_BYTES, clean_line, getline_lines, parse_permissive_sequence, read_text,
};
use crate::identification::{
    PeptideHit, PeptideIdentification, ProteinIdentification, SearchParameters,
};
use crate::metadata::MetaValue;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// Score type of every identification.
pub const SCORE_TYPE: &str = "PepNovo";

/// Maps PepNovo spectrum indices to the (RT, m/z) of their precursor (the
/// source's `IndexPosMappingType`).
pub type IndexPosMapping = BTreeMap<usize, (f64, f64)>;

/// The eight table columns the reader requires, by the key the source maps
/// them under.
const COLUMNS: [(&str, &str); 8] = [
    ("#Index", "Index"),
    ("RnkScr", "RnkScr"),
    ("PnvScr", "PnvScr"),
    ("N-Gap", "N-Gap"),
    ("C-Gap", "C-Gap"),
    ("[M+H]", "[M+H]"),
    ("Charge", "Charge"),
    ("Sequence", "Sequence"),
];

/// Everything one PepNovo output file yields.
#[derive(Clone, Debug, PartialEq)]
pub struct PepNovoResult {
    /// Search engine, version and search parameters from the header, and the
    /// identifier `<engine>_<run date>`; no protein hits.
    pub protein_identification: ProteinIdentification,
    /// One identification per spectrum block with kept candidates (see
    /// [`PepNovoOutfile::load`] for the one exception).
    pub peptide_identifications: Vec<PeptideIdentification>,
}

/// Reader for PepNovo output files. The source's copy constructor,
/// assignment and `operator==` (always true) are the derived ones.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PepNovoOutfile {
    /// Date of the run, which names the identifier. `None` takes the local
    /// clock, as the source's `DateTime::now()` does (PepNovo records no
    /// date).
    pub run_date: Option<DateTime>,
}

fn located(line: usize) -> impl Fn(String) -> Error {
    move |message| Error::Parse { line, message }
}

impl PepNovoOutfile {
    /// A reader that dates its runs by the local clock.
    pub fn new() -> Self {
        Self::default()
    }

    /// Load a PepNovo output file.
    ///
    /// Each `>> ` line starts a spectrum; its third space-separated field is
    /// the spectrum index. When `index_to_precursor` is not empty, the index
    /// must be in it and gives the RT and m/z; otherwise the fourth field is
    /// read as `<m/z>_<RT>[_...]`. Candidates whose rank score (`RnkScr`) is at
    /// least `score_threshold` become hits: charge, rank (`#Index` + 1), the
    /// rank score as score, the meta values `PnvScr`, `N-Gap`, `C-Gap` and
    /// `MZ` (from `[M+H]`), all read as `float`s, and the sequence with every
    /// PTM key of the header's `PTMs considered` replaced: keys found in
    /// `modkey_to_mod_id` become the UniMod modification (`K+42` with
    /// `Acetyl (K)` becomes `K(Acetyl)`, an N-terminal `^` key a leading
    /// `(Acetyl)`), other keys a mass tag (`K+42` becomes `K[+42]`).
    ///
    /// Every identification carries `score_threshold` as its significance
    /// threshold. An identification is reported when it has hits; after the
    /// last block, the source also reports it without hits when an earlier
    /// one was reported, and the port does the same.
    ///
    /// The source first stores the run date in the protein identification
    /// and then replaces that with a default one, so its result carries no
    /// date; the port's does not either.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be read, and [`Error::Parse`] with
    /// the source's messages: "Not enough columns (spectrum Id) in file in
    /// line N (should be 2 or more)!", "Expected an index number in line N at
    /// position 2 (line was: '...')!", "Index 'I' in line '...' not found in
    /// index table (line was: '...')!", "Precursor could not be
    /// reconstructed from title '...' in line '...' (line was: '...')!", and
    /// "Not enough columns in file in line N (should be 8)!" for a table
    /// header without the eight columns or a candidate line of another width.
    /// Also [`Error::Parse`] for a candidate field that does not convert
    /// (the source's `ConversionError`), a non-finite rank score, or a
    /// sequence that does not parse, and the lookup's error for a
    /// modification in `modkey_to_mod_id` that is not in
    /// [`ModificationsDB::global`].
    pub fn load(
        &self,
        path: impl AsRef<Path>,
        score_threshold: f64,
        index_to_precursor: &IndexPosMapping,
        modkey_to_mod_id: &BTreeMap<String, String>,
    ) -> Result<PepNovoResult> {
        let text = read_text(BufReader::new(File::open(path)?), MAX_TEXT_BYTES)?;
        let mut protein_identification = ProteinIdentification::default();
        parse_search_engine_and_version(&text, &mut protein_identification)?;
        if protein_identification.search_engine_version.is_empty() {
            protein_identification.search_engine = "PepNovo".into();
            protein_identification.search_engine_version = "unknown".into();
        }
        let date = self.run_date.unwrap_or_else(DateTime::now);
        let identifier = format!(
            "{}_{}",
            protein_identification.search_engine,
            date.date_string()
        );
        protein_identification.identifier = identifier.clone();
        let masks = modification_masks(
            &protein_identification
                .search_parameters
                .variable_modifications,
            modkey_to_mod_id,
        )?;

        let lines = getline_lines(&text);
        let mut peptide_identifications = Vec::new();
        let mut current = PeptideIdentification::default();
        let mut columns: BTreeMap<&str, usize> = BTreeMap::new();
        let mut cursor = 0usize;
        while cursor < lines.len() {
            let line = clean_line(lines[cursor]);
            cursor += 1;
            let line_number = cursor;
            if line.starts_with(">> ") {
                if !current.hits.is_empty() {
                    peptide_identifications.push(std::mem::take(&mut current));
                }
                current = self.spectrum(line, line_number, score_threshold, index_to_precursor)?;
                current.identifier = identifier.clone();
            } else if line.starts_with("#Index") {
                if columns.is_empty() {
                    for (position, name) in split(line, '\t').into_iter().enumerate() {
                        if let Some((_, key)) = COLUMNS.iter().find(|(column, _)| *column == name) {
                            columns.insert(key, position);
                        }
                    }
                    if columns.len() != COLUMNS.len() {
                        return Err(Error::Parse {
                            line: line_number,
                            message: format!(
                                "Not enough columns in file in line {line_number} (should be 8)!"
                            ),
                        });
                    }
                }
                while cursor < lines.len() {
                    let candidate = clean_line(lines[cursor]);
                    cursor += 1;
                    let candidate_number = cursor;
                    if candidate.is_empty() {
                        break;
                    }
                    let fields = split(candidate, '\t');
                    if fields.len() != COLUMNS.len() {
                        return Err(Error::Parse {
                            line: candidate_number,
                            message: format!(
                                "Not enough columns in file in line {candidate_number} (should be 8)!"
                            ),
                        });
                    }
                    if let Some(hit) =
                        candidate_hit(&fields, &columns, score_threshold, &masks, candidate_number)?
                    {
                        current.hits.push(hit);
                    }
                }
            }
        }
        if !peptide_identifications.is_empty() || !current.hits.is_empty() {
            peptide_identifications.push(current);
        }
        Ok(PepNovoResult {
            protein_identification,
            peptide_identifications,
        })
    }

    /// One `>> ` line's identification, without its hits.
    fn spectrum(
        &self,
        line: &str,
        line_number: usize,
        score_threshold: f64,
        index_to_precursor: &IndexPosMapping,
    ) -> Result<PeptideIdentification> {
        let fields = split(line, ' ');
        if fields.len() < 3 {
            return Err(Error::Parse {
                line: line_number,
                message: format!(
                    "Not enough columns (spectrum Id) in file in line {line_number} (should be 2 or more)!"
                ),
            });
        }
        let index = to_int32(trim(fields[2])).map_err(|_| Error::Parse {
            line: line_number,
            message: format!(
                "Expected an index number in line {line_number} at position 2 (line was: '{line}')!"
            ),
        })?;
        // The source stores the `Int` in a `Size`.
        let index = index as isize as usize;
        let mut identification = PeptideIdentification {
            score_type: SCORE_TYPE.into(),
            significance_threshold: score_threshold,
            ..PeptideIdentification::default()
        };
        if !index_to_precursor.is_empty() {
            let Some(&(rt, mz)) = index_to_precursor.get(&index) else {
                return Err(Error::Parse {
                    line: line_number,
                    message: format!(
                        "Index '{index}' in line '{line}' not found in index table (line was: '{line}')!"
                    ),
                });
            };
            identification.rt = Some(rt);
            identification.mz = Some(mz);
            return Ok(identification);
        }
        // The source reads the fourth field even when there is none, which
        // is undefined behaviour; the port reports an empty title.
        let title = fields.get(3).copied().unwrap_or("");
        let parts = split(title, '_');
        let reconstructed = if fields.len() >= 4 && parts.len() >= 2 {
            match to_double(parts[1]) {
                Ok(rt) => {
                    identification.rt = Some(rt);
                    to_double(parts[0]).ok().map(|mz| (rt, mz))
                }
                Err(_) => None,
            }
        } else {
            None
        };
        match reconstructed {
            Some((rt, mz)) if rt.is_finite() && mz.is_finite() => {
                identification.rt = Some(rt);
                identification.mz = Some(mz);
                Ok(identification)
            }
            _ => Err(Error::Parse {
                line: line_number,
                message: format!(
                    "Precursor could not be reconstructed from title '{title}' in line '{line}' (line was: '{line}')!"
                ),
            }),
        }
    }

    /// Read the search engine, version and search parameters from the header
    /// of a PepNovo output, e.g. the output of a PepNovo run without
    /// parameters (source `getSearchEngineAndVersion`).
    ///
    /// Up to the first `>>` line: a `PepNovo...` line gives the engine and
    /// version (`PepNovo+, Build 20081230` or `PepNovo v2.00 - ...`, split at
    /// the comma or, with three space-separated words, at the spaces); a
    /// `PM...` line the precursor and a `Fragment...` line the fragment
    /// tolerance, each its last space-separated word; a `PTM...` line the
    /// colon-separated variable modification keys unless the first is
    /// `None`. The engine and version are only replaced when found; the
    /// search parameters are always replaced by ones holding just these
    /// values.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be read; [`Error::Parse`] when a
    /// tolerance does not convert or is not finite. On error
    /// `protein_identification` is unchanged.
    pub fn get_search_engine_and_version(
        &self,
        path: impl AsRef<Path>,
        protein_identification: &mut ProteinIdentification,
    ) -> Result<()> {
        let text = read_text(BufReader::new(File::open(path)?), MAX_TEXT_BYTES)?;
        let mut updated = protein_identification.clone();
        parse_search_engine_and_version(&text, &mut updated)?;
        *protein_identification = updated;
        Ok(())
    }
}

fn tolerance(text: &str, line: usize) -> Result<Tolerance> {
    let value = f64::from(to_float(text).map_err(located(line))?);
    if value.is_finite() {
        Ok(Tolerance::Absolute(value))
    } else {
        Err(located(line)(format!("non-finite tolerance '{text}'")))
    }
}

fn parse_search_engine_and_version(
    text: &str,
    protein_identification: &mut ProteinIdentification,
) -> Result<()> {
    let mut parameters = SearchParameters::default();
    for (index, raw) in getline_lines(text).into_iter().enumerate() {
        let line = clean_line(raw);
        let number = index + 1;
        if line.is_empty() {
            continue;
        }
        if line.starts_with("PepNovo") {
            let parts = split(line, ',');
            if parts.len() == 2 {
                protein_identification.search_engine = trim(parts[0]).to_owned();
                protein_identification.search_engine_version = trim(parts[1]).to_owned();
            } else {
                let words = split(line, ' ');
                if words.len() == 3 {
                    protein_identification.search_engine = trim(words[0]).to_owned();
                    protein_identification.search_engine_version = trim(words[2]).to_owned();
                }
            }
        }
        if line.starts_with("PM") {
            let words = split(line, ' ');
            parameters.precursor_tolerance =
                tolerance(words.last().copied().unwrap_or(""), number)?;
        }
        if line.starts_with("Fragment") {
            let words = split(line, ' ');
            parameters.fragment_tolerance = tolerance(words.last().copied().unwrap_or(""), number)?;
        }
        if line.starts_with("PTM") {
            let keys: Vec<String> = split(line, ':')
                .into_iter()
                .skip(1)
                .map(|key| trim(key).to_owned())
                .collect();
            if keys.first().is_some_and(|first| first != "None") {
                parameters.variable_modifications = keys;
            }
        }
        if line.starts_with(">>") {
            break;
        }
    }
    protein_identification.search_parameters = parameters;
    Ok(())
}

/// Source `mod_mask_map`: the replacement text of each PTM key.
fn modification_masks(
    keys: &[String],
    modkey_to_mod_id: &BTreeMap<String, String>,
) -> Result<BTreeMap<String, String>> {
    let mut masks = BTreeMap::new();
    for key in keys.iter().filter(|key| !key.is_empty()) {
        let terminal = key.starts_with('^') || key.starts_with('$');
        let mask = match modkey_to_mod_id.get(key) {
            Some(id) => {
                let modification = ModificationsDB::global().get_modification(id, None, None)?;
                if terminal {
                    format!("({})", modification.name())
                } else {
                    format!(
                        "{}({})",
                        modification.origin().unwrap_or('X'),
                        modification.name()
                    )
                }
            }
            None if terminal => format!("[{key}]"),
            None => {
                let mut chars = key.chars();
                let first = chars.next().map(String::from).unwrap_or_default();
                format!("{first}[{}]", chars.as_str())
            }
        };
        masks.insert(key.clone(), mask);
    }
    Ok(masks)
}

/// One candidate line's hit, or `None` below the threshold.
fn candidate_hit(
    fields: &[&str],
    columns: &BTreeMap<&str, usize>,
    score_threshold: f64,
    masks: &BTreeMap<String, String>,
    line: usize,
) -> Result<Option<PeptideHit>> {
    let field = |key: &str| fields[columns[key]];
    let float = |key: &str| -> Result<f64> {
        let value = f64::from(to_float(field(key)).map_err(located(line))?);
        if value.is_finite() {
            Ok(value)
        } else {
            Err(located(line)(format!("non-finite {key} '{}'", field(key))))
        }
    };
    let score = float("RnkScr")?;
    if score < score_threshold {
        return Ok(None);
    }
    let mut hit = PeptideHit {
        charge: to_int32(field("Charge")).map_err(located(line))?,
        rank: to_int32(field("Index"))
            .map_err(located(line))?
            .wrapping_add(1) as u32,
        score,
        ..PeptideHit::default()
    };
    for (name, key) in [
        ("PnvScr", "PnvScr"),
        ("N-Gap", "N-Gap"),
        ("C-Gap", "C-Gap"),
        ("MZ", "[M+H]"),
    ] {
        hit.metadata
            .insert(name.into(), MetaValue::try_from(float(key)?)?);
    }
    let mut sequence = field("Sequence").to_owned();
    for (key, mask) in masks {
        if key.starts_with('^') && sequence.contains(key.as_str()) {
            sequence = sequence.replace(key.as_str(), "");
            sequence = format!("{mask}{sequence}");
        }
        sequence = sequence.replace(key.as_str(), mask);
    }
    hit.sequence = parse_permissive_sequence(&sequence).map_err(|error| match error {
        Error::Parse { message, .. } => Error::Parse { line, message },
        other => other,
    })?;
    Ok(Some(hit))
}
