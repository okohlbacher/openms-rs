// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Reader for Percolator's PSM-level tab-separated output, ported from
//! `FORMAT/PercolatorOutfile.h`.
//!
//! For PSM-level output the file extension should be `.psms`. The first row
//! must be exactly the header [`HEADER`]; every later row holds one PSM: the
//! spectrum reference `PSMId`, the three Percolator scores, the peptide with
//! its flanking residues (`K.PEPTIDE.A`) and any number of protein accessions.
//! Lines starting with `#` are comments, as for every `CsvFile`.
//!
//! Retention time, precursor m/z and charge come from the spectrum reference
//! through a [`SpectrumMetaDataLookup`]: either the reference carries them
//! itself (Mascot-Percolator references do), or it names a spectrum the caller
//! read into the lookup beforehand (MS-GF+ and X! Tandem references name an
//! index).
//!
//! Entry points: [`PercolatorOutfile::load`] and [`PercolatorOutfile::read`],
//! with [`ScoreType::from_name`] for the source's `getScoreType`. What is and is
//! not reproduced is listed in `docs/PERCOLATOR_OUTFILE_SUPPORT.md`.

use crate::chemistry::{AASequence, ModificationDefinitionsSet, ModificationsDB, TermSpecificity};
use crate::concept::boost_regex::BoostRegex;
use crate::data_structures::DateTime;
use crate::format::csv::{self, CsvFile};
use crate::format::spectrum_meta_data_lookup::{
    MetaDataFlags, SpectrumMetaData, SpectrumMetaDataLookup, to_double,
};
use crate::identification::{
    PeptideEvidence, PeptideHit, PeptideIdentification, ProteinHit, ProteinIdentification,
    SearchParameters,
};
use crate::metadata::MetaValue;
use crate::{Error, Result};
use std::collections::BTreeSet;
use std::io::BufRead;
use std::path::Path;
use std::str::FromStr;

/// The only header the source accepts, tab-separated.
pub const HEADER: &str = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds";

/// Names of the Percolator scores, in [`ScoreType`] order, the source
/// `score_type_names`. `IDFileConverter -score_type` offers exactly these.
pub const SCORE_TYPE_NAMES: [&str; 3] = ["qvalue", "PEP", "score"];

/// The three spectrum reference formats `load` registers when the lookup has
/// none, in the source's order: MS-GF+ Percolator (`INDEX1` and `CHARGE`),
/// Mascot Percolator (`INDEX0`, `RT`, `MZ` and `CHARGE`; the retention time may
/// be empty, e.g. for searches through Proteome Discoverer) and X! Tandem
/// Percolator (`INDEX0` and `CHARGE`).
///
/// The Mascot expression writes `[(scans:)(scan=)(spectrum=)]`, which Boost
/// reads as a *character class* of those characters, not as three
/// alternatives; it is kept verbatim because that is what the source matches.
pub const DEFAULT_REFERENCE_FORMATS: [&str; 3] = [
    r"_SII_(?<INDEX1>\d+)_\d+_\d+_(?<CHARGE>\d+)_\d+",
    r"spectrum:[^;]+[(scans:)(scan=)(spectrum=)](?<INDEX0>\d+)[^;]+;rt:(?<RT>\d*(\.\d+)?);mz:(?<MZ>\d+(\.\d+)?);charge:(?<CHARGE>-?\d+)",
    r"_(?<INDEX0>\d+)_(?<CHARGE>\d+)_\d+$",
];

/// Which Percolator score becomes a hit's score, the source's `ScoreType`
/// (whose `SIZE_OF_SCORETYPE` sentinel has no counterpart here).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum ScoreType {
    /// The q-value; lower is better. The source default.
    #[default]
    QValue,
    /// The posterior error probability; lower is better.
    PosteriorErrorProbability,
    /// The raw SVM score; higher is better.
    Score,
}

impl ScoreType {
    /// The score type named `name`, the source's `getScoreType`.
    ///
    /// Case is ignored (ASCII only, as `StringUtils::toLower` in the C locale).
    /// `q-value`, `qvalue` and `q value` select [`ScoreType::QValue`]; `pep`
    /// and `posterior error probability` select
    /// [`ScoreType::PosteriorErrorProbability`]; `score` selects
    /// [`ScoreType::Score`].
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for any other name, the source's
    /// `Exception::InvalidValue` "Not a valid Percolator score type".
    pub fn from_name(name: &str) -> Result<Self> {
        match name.to_ascii_lowercase().as_str() {
            "q-value" | "qvalue" | "q value" => Ok(Self::QValue),
            "pep" | "posterior error probability" => Ok(Self::PosteriorErrorProbability),
            "score" => Ok(Self::Score),
            _ => Err(Error::InvalidValue(format!(
                "Not a valid Percolator score type: {name}"
            ))),
        }
    }

    /// The entry of [`SCORE_TYPE_NAMES`] for this score type.
    pub fn name(self) -> &'static str {
        SCORE_TYPE_NAMES[self as usize]
    }

    /// The score type and higher-is-better flag the peptide identification
    /// records for this choice.
    fn identification_score(self) -> (&'static str, bool) {
        match self {
            Self::Score => ("Percolator_score", true),
            Self::QValue => ("q-value", false),
            Self::PosteriorErrorProbability => ("Posterior Error Probability", false),
        }
    }
}

impl FromStr for ScoreType {
    type Err = Error;
    fn from_str(name: &str) -> Result<Self> {
        Self::from_name(name)
    }
}

/// What one load produced.
#[derive(Clone, Debug, PartialEq)]
pub struct PercolatorDocument {
    /// The search run: identifier `id`, search engine `Percolator`, the load
    /// time, one hit per distinct accession in byte order, and the fixed and
    /// variable modifications inferred from the peptides.
    pub protein_identification: ProteinIdentification,
    /// One identification with one hit per PSM row, in file order.
    pub peptide_identifications: Vec<PeptideIdentification>,
    /// The source's log output, in order: one error per spectrum reference
    /// that could not be resolved, one warning per peptide whose `[unknown]`
    /// modifications were removed, the two "Created ..." counts, and one
    /// warning per kind of missing charge, retention time or m/z.
    pub log: Vec<String>,
}

/// The source's `PercolatorOutfile`, a stateless reader.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PercolatorOutfile;

impl PercolatorOutfile {
    /// A reader, the source's default constructor.
    pub fn new() -> Self {
        Self
    }

    /// Load a Percolator PSM file, as `load(filename, proteins, peptides,
    /// lookup, output_score)`.
    ///
    /// If `lookup` has no reference format, the three
    /// [`DEFAULT_REFERENCE_FORMATS`] are registered in it first — the lookup is
    /// modified even if the load then fails, as in the source. Each PSM's
    /// retention time, precursor m/z and precursor charge are then taken from
    /// its reference; a reference that cannot be resolved is logged, not
    /// fatal, and leaves whatever the reference itself supplied before the
    /// failure.
    ///
    /// All three scores are stored as the meta values `Percolator_score`,
    /// `Percolator_qvalue` and `Percolator_PEP`; `output_score` picks the one
    /// that becomes the hit's score, with score type `Percolator_score` (higher
    /// is better), `q-value` or `Posterior Error Probability` (lower is
    /// better).
    ///
    /// # Errors
    ///
    /// - I/O and text-limit errors of [`CsvFile`];
    /// - [`Error::InvalidValue`] when the first row is not [`HEADER`] (the
    ///   source's `Exception::ParseError`), when a score does not convert (the
    ///   source's `Exception::ConversionError`), when a peptide does not parse
    ///   as an `AASequence`, or when inferring the modification definitions
    ///   fails.
    ///
    /// Stricter than the source, each an error rather than undefined or
    /// unrepresentable behaviour: a row with fewer than five fields (the
    /// source reads `items[1]` to `items[4]` unchecked), a peptide shorter than
    /// two characters (the source reads `peptide[len - 2]`), and a score that
    /// is not finite (a hit's score and its meta values must be finite here).
    pub fn load(
        &self,
        path: impl AsRef<Path>,
        lookup: &mut SpectrumMetaDataLookup,
        output_score: ScoreType,
    ) -> Result<PercolatorDocument> {
        let csv = CsvFile::from_path(path, &read_options())?;
        self.load_csv(&csv, lookup, output_score)
    }

    /// Read Percolator PSM output from `reader`; see [`Self::load`].
    ///
    /// # Errors
    ///
    /// As for [`Self::load`].
    pub fn read(
        &self,
        reader: impl BufRead,
        lookup: &mut SpectrumMetaDataLookup,
        output_score: ScoreType,
    ) -> Result<PercolatorDocument> {
        let csv = CsvFile::from_reader(reader, &read_options())?;
        self.load_csv(&csv, lookup, output_score)
    }

    fn load_csv(
        &self,
        source: &CsvFile,
        lookup: &mut SpectrumMetaDataLookup,
        output_score: ScoreType,
    ) -> Result<PercolatorDocument> {
        if lookup.reference_formats().is_empty() {
            for format in DEFAULT_REFERENCE_FORMATS {
                lookup.add_reference_format(format)?;
            }
        }
        let lookup = &*lookup;
        let db = ModificationsDB::global();
        let rewrite = Rewrites::new()?;
        let flags =
            MetaDataFlags::RT | MetaDataFlags::PRECURSOR_MZ | MetaDataFlags::PRECURSOR_CHARGE;

        // The source reads row 0 even from an empty file, where getRow throws.
        if source.row_count() == 0 {
            return Err(Error::InvalidValue(
                "Not a valid header for Percolator (PSM level) output: the file has no rows".into(),
            ));
        }
        let (_, items) = source.row(0)?;
        let header = items.join("\t");
        if header != HEADER {
            return Err(Error::InvalidValue(format!(
                "Not a valid header for Percolator (PSM level) output: {header}"
            )));
        }

        let mut log = Vec::new();
        let mut accessions = BTreeSet::new();
        let (mut no_charge, mut no_rt, mut no_mz) = (0usize, 0usize, 0usize);
        let mut peptides = Vec::new();
        peptides
            .try_reserve_exact(source.row_count() - 1)
            .map_err(|_| Error::InvalidValue("Percolator PSM allocation failed".into()))?;
        let (score_type, higher_score_better) = output_score.identification_score();

        for row in 1..source.row_count() {
            let (_, items) = source.row(row)?;
            if items.len() < 5 {
                return Err(Error::InvalidValue(format!(
                    "Percolator output row {row} has {} fields; at least five are required",
                    items.len()
                )));
            }
            let mut meta = SpectrumMetaData::default();
            if lookup
                .meta_data_for_reference(&items[0], &mut meta, flags)
                .is_err()
            {
                log.push(format!(
                    "Error: Could not extract data for spectrum reference '{}' from row {row}",
                    items[0]
                ));
            }

            let mut hit = PeptideHit::default();
            if meta.precursor_charge != 0 {
                hit.charge = meta.precursor_charge;
            } else {
                no_charge += 1;
            }
            let mut peptide = PeptideIdentification {
                identifier: "id".into(),
                score_type: score_type.into(),
                higher_score_better,
                ..PeptideIdentification::default()
            };
            match meta.rt {
                Some(rt) => peptide.rt = Some(rt),
                None => no_rt += 1,
            }
            match meta.precursor_mz {
                Some(mz) => peptide.mz = Some(mz),
                None => no_mz += 1,
            }

            let score = finite_score(&items[1], row, "score")?;
            let qvalue = finite_score(&items[2], row, "q-value")?;
            let posterior = finite_score(&items[3], row, "posterior_error_prob")?;
            hit.metadata
                .insert("Percolator_score".into(), MetaValue::try_from(score)?);
            hit.metadata
                .insert("Percolator_qvalue".into(), MetaValue::try_from(qvalue)?);
            hit.metadata
                .insert("Percolator_PEP".into(), MetaValue::try_from(posterior)?);
            hit.score = match output_score {
                ScoreType::Score => score,
                ScoreType::QValue => qvalue,
                ScoreType::PosteriorErrorProbability => posterior,
            };

            let (text, warning) = rewrite.peptide(&items[4], db)?;
            if let Some(warning) = warning {
                log.push(warning);
            }
            hit.sequence = AASequence::parse_with_registry(&text, db)?;

            for accession in &items[5..] {
                accessions.insert(accession.clone());
                hit.evidences.push(PeptideEvidence {
                    protein_accession: accession.clone(),
                    ..PeptideEvidence::default()
                });
            }
            peptide.hits.push(hit);
            peptides.push(peptide);
        }

        let mut proteins = ProteinIdentification {
            identifier: "id".into(),
            date_time: Some(DateTime::now().iso_string()),
            search_engine: "Percolator".into(),
            ..ProteinIdentification::default()
        };
        proteins.hits = accessions
            .into_iter()
            .map(|accession| ProteinHit {
                accession,
                ..ProteinHit::default()
            })
            .collect();

        let mut definitions = ModificationDefinitionsSet::default();
        definitions.infer_from_peptides(&peptides)?;
        proteins.search_parameters = SearchParameters {
            fixed_modifications: definitions.fixed_names().into_iter().collect(),
            variable_modifications: definitions.variable_names().into_iter().collect(),
            ..SearchParameters::default()
        };

        log.push(format!(
            "Created {} protein hits.\nCreated {} peptide hits (PSMs).",
            proteins.hits.len(),
            peptides.len()
        ));
        if no_charge > 0 {
            log.push(format!(
                "{no_charge} peptide hits without charge state information."
            ));
        }
        if no_rt > 0 {
            log.push(format!(
                "{no_rt} peptide hits without retention time information."
            ));
        }
        if no_mz > 0 {
            log.push(format!(
                "{no_mz} peptide hits without mass-to-charge information."
            ));
        }
        Ok(PercolatorDocument {
            protein_identification: proteins,
            peptide_identifications: peptides,
            log,
        })
    }
}

/// `CsvFile(filename, '\t')`: tab-separated, not enclosed, all rows, `#`
/// comments dropped and lines kept untrimmed.
fn read_options() -> csv::ReadOptions {
    csv::ReadOptions {
        separator: b'\t',
        ..csv::ReadOptions::default()
    }
}

/// A score through `StringUtils::toDouble`, which must also be finite here.
fn finite_score(text: &str, row: usize, column: &str) -> Result<f64> {
    let value = to_double(text)?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::InvalidValue(format!(
            "Percolator output row {row}: {column} '{text}' is not finite"
        )))
    }
}

/// The source's peptide rewriting, `getPeptideSequence_` and
/// `resolveMisassignedNTermMods_`, with its three expressions compiled once.
struct Rewrites {
    unimod: BoostRegex,
    misassigned: BoostRegex,
    positive_shift: BoostRegex,
}

impl Rewrites {
    fn new() -> Result<Self> {
        Ok(Self {
            unimod: BoostRegex::new(r"\[UNIMOD:(\d+)\]")?,
            misassigned: BoostRegex::new(
                r"^[A-Z]\[(?<MOD1>-?\d+(\.\d+)?)\](\[(?<MOD2>-?\d+(\.\d+)?)\])?",
            )?,
            positive_shift: BoostRegex::new(r"\[(\d)")?,
        })
    }

    /// The `AASequence` text for a Percolator `peptide` column, and the
    /// warning the source logs when it removes `[unknown]` modifications.
    ///
    /// The flanking residues are dropped (it is unclear which protein they
    /// belong to): the first two characters when the second is `.`, the last
    /// two when the second-to-last is `.`. Then every `[unknown]` is removed,
    /// `[UNIMOD:n]` becomes `(UniMod:n)`, misassigned X! Tandem N-terminal mass
    /// shifts are moved to the terminus, and a positive mass shift gains its
    /// missing `+`.
    fn peptide(&self, peptide: &str, db: &ModificationsDB) -> Result<(String, Option<String>)> {
        let bytes = peptide.as_bytes();
        let length = bytes.len();
        if length < 2 {
            return Err(Error::InvalidValue(format!(
                "Percolator peptide '{peptide}' is shorter than two characters"
            )));
        }
        let start = if bytes[1] == b'.' { 2 } else { 0 };
        // The source computes `len - start - 2` in `size_t`; when that wraps
        // (e.g. `K..`) `substr` takes the rest of the string, as `None` does.
        let end = if bytes[length - 2] == b'.' {
            length.checked_sub(start + 2).map(|count| start + count)
        } else {
            None
        };
        let mut text = match end {
            Some(end) => peptide.get(start..end),
            None => peptide.get(start..),
        }
        .ok_or_else(|| {
            Error::InvalidValue(format!("Percolator peptide '{peptide}' is not valid text"))
        })?
        .to_owned();

        let mut warning = None;
        if text.contains("[unknown]") {
            warning = Some(format!(
                "Removing unknown modification(s) from peptide '{text}'"
            ));
            text = text.replace("[unknown]", "");
        }
        text = replace_all(&self.unimod, &text, |groups| {
            format!("(UniMod:{})", groups[1])
        })?;
        text = self.resolve_misassigned_n_term_mods(text, db)?;
        text = replace_all(&self.positive_shift, &text, |groups| {
            format!("[+{}", groups[1])
        })?;
        Ok((text, warning))
    }

    /// `resolveMisassignedNTermMods_`: X! Tandem may report an N-terminal mass
    /// shift on the first residue, and up to two shifts (one terminal, one on
    /// the residue) in either order.
    ///
    /// Candidates are the best registry matches within 0.01 Da for the first
    /// residue, at the N-terminus and anywhere. One shift that can be
    /// N-terminal becomes the terminal modification, unless it is
    /// carbamidomethylation on C. Of two shifts, the one that can be
    /// N-terminal is moved to the terminus; if both can, the one that cannot
    /// sit on the residue is, and otherwise the first. If neither can, the text
    /// is left for `AASequence` to reject or accept.
    fn resolve_misassigned_n_term_mods(
        &self,
        peptide: String,
        db: &ModificationsDB,
    ) -> Result<String> {
        let Some(found) = self.misassigned.search(peptide.as_bytes())? else {
            return Ok(peptide);
        };
        let group = |name: &str| found.name(name).and_then(|range| peptide.get(range));
        let Some(mod1) = group("MOD1") else {
            return Ok(peptide);
        };
        // The source takes the residue from `peptide[0]`, not from the match,
        // which differs only when `^` matched after a line separator inside
        // the text. A non-ASCII first byte has no single-character spelling
        // here, so that shape is refused rather than re-encoded.
        let first = peptide.as_bytes()[0];
        if !first.is_ascii() {
            return Err(Error::Unsupported(format!(
                "Percolator peptide '{peptide}' starts with a non-ASCII byte"
            )));
        }
        let residue = char::from(first);
        let best = |text: &str, term: TermSpecificity| -> Result<Option<String>> {
            let mass = to_double(text)?;
            Ok(db
                .best_by_mass_handle(mass, 0.01, Some(residue), Some(term))?
                .map(|modification| modification.name().to_owned()))
        };
        let n_term1 = best(mod1, TermSpecificity::NTerm)?;
        let replacement = match group("MOD2") {
            None => match n_term1 {
                // Only one shift, which may be terminal: assume it is, unless
                // it is carbamidomethylation on C.
                Some(id) if id != "Carbamidomethyl" || residue != 'C' => {
                    Some(format!(".({id}){residue}"))
                }
                _ => None,
            },
            Some(mod2) => {
                let n_term2 = best(mod2, TermSpecificity::NTerm)?;
                match (&n_term1, &n_term2) {
                    (Some(first), None) => Some(format!("({first}){residue}[{mod2}]")),
                    (None, Some(second)) => Some(format!("({second}){residue}[{mod1}]")),
                    (Some(first), Some(second)) => {
                        let on_residue1 = best(mod1, TermSpecificity::Anywhere)?.is_some();
                        let on_residue2 = best(mod2, TermSpecificity::Anywhere)?.is_some();
                        if on_residue1 && !on_residue2 {
                            Some(format!("({second}){residue}[{mod1}]"))
                        } else {
                            // The second shift must sit on the residue, or
                            // both could be either: the first is terminal.
                            Some(format!("({first}){residue}[{mod2}]"))
                        }
                    }
                    // Neither can be terminal: something is wrong, and the
                    // source lets AASequence deal with it.
                    (None, None) => None,
                }
            }
        };
        match replacement {
            Some(replacement) => replace_all(&self.misassigned, &peptide, |_| replacement.clone()),
            None => Ok(peptide),
        }
    }
}

/// `boost::regex_replace(text, regex, format)` with the default flags: every
/// non-overlapping match, left to right, is replaced by `replace` applied to
/// the match's groups (group 0 first); the text between matches is copied.
///
/// Built on the token iterator with the sub-match list `-1, 0, 1, ...`, which
/// yields for every match the text before it and then each group, and after
/// the last match the non-empty rest of the text. Match boundaries of the
/// expressions used here are ASCII bytes, so every piece is valid UTF-8.
fn replace_all(
    regex: &BoostRegex,
    text: &str,
    replace: impl Fn(&[&str]) -> String,
) -> Result<String> {
    let haystack = text.as_bytes();
    let groups = i32::try_from(regex.mark_count())
        .map_err(|_| Error::InvalidValue("too many regular-expression groups".into()))?;
    let submatches: Vec<i32> = std::iter::once(-1).chain(0..=groups).collect();
    let tokens = regex
        .tokens(haystack, &submatches)?
        .collect::<Result<Vec<_>>>()?;
    let piece = |token: &crate::concept::boost_regex::SubMatch| -> Result<&str> {
        std::str::from_utf8(token.as_bytes(haystack))
            .map_err(|_| Error::InvalidValue("regular-expression match splits a character".into()))
    };
    let mut result = String::new();
    for chunk in tokens.chunks(submatches.len()) {
        result.push_str(piece(&chunk[0])?);
        if chunk.len() == submatches.len() {
            let groups = chunk[1..].iter().map(piece).collect::<Result<Vec<_>>>()?;
            result.push_str(&replace(&groups));
        }
    }
    Ok(result)
}
