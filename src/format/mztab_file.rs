// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! MzTab file adapter: reading and writing `.mzTab` documents (`FORMAT/MzTabFile.h`).
//!
//! An MzTab file is flat tab-separated text. Every line begins with a
//! three-letter tag: `MTD` for one metadata key/value pair, `COM` for a comment,
//! and a header/data tag pair per section — `PRH`/`PRT` proteins, `PEH`/`PEP`
//! peptides, `PSH`/`PSM` peptide-spectrum matches, `SMH`/`SML` small molecules
//! and, in the OpenMS extension, `NUH`/`NUC` nucleic acids, `OLH`/`OLI`
//! oligonucleotides and `OSH`/`OSM` oligonucleotide-spectrum matches.
//!
//! The data model — every cell type, every record struct and the document
//! itself — lives in [`crate::format::mztab`]. This module is only the file
//! layer: it turns a [`MzTab`](crate::format::mztab::MzTab) into lines and lines
//! back into a [`MzTab`](crate::format::mztab::MzTab).
//!
//! # Where the work is
//!
//! Each section has a fixed prefix of required columns followed by an
//! open-ended set of optional `opt_…` columns, and several of the required
//! columns are themselves *families* indexed in brackets:
//! `best_search_engine_score[1]`, `search_engine_score[2]_ms_run[3]`,
//! `num_psms_ms_run[1]`, `protein_abundance_assay[4]`,
//! `protein_abundance_stdev_study_variable[2]` and so on. Which members of a
//! family appear is a property of the document, not of the format, so writing a
//! section means first deciding its column layout — that is
//! [`SectionLayout`](crate::format::mztab_file::SectionLayout) — and then
//! rendering every row against it. Reading is the mirror image: the header row
//! says which column holds which family member, in any order, and the reader
//! must tolerate unknown names and absent optional families.
//!
//! # Fidelity
//!
//! A document written by
//! [`MzTabFile::store`](crate::format::mztab_file::MzTabFile::store) and read
//! back with [`MzTabFile::load`](crate::format::mztab_file::MzTabFile::load)
//! compares equal, including the `null`/`NaN`/`Inf` distinction of every numeric
//! cell and including an optional column absent from a row rather than null in
//! it. Reaching that required diverging from the source in the several places
//! where its reader and writer disagree; each one is listed in
//! `docs/MZTAB_FILE_SUPPORT.md` with the source line behind it.
//!
//! ```
//! use openms::format::mztab::{MzTab, MzTabDouble, MzTabPSMSectionRow, MzTabParameter};
//! use openms::format::mztab_file::MzTabFile;
//!
//! let mut document = MzTab::default();
//! document.meta_data.mz_tab_mode.set("Summary");
//! document.meta_data.mz_tab_type.set("Identification");
//! document.meta_data.description.set("round trip");
//! document.meta_data.psm_search_engine_score.insert(
//!     1,
//!     MzTabParameter::parse("[MS, MS:1001171, Mascot:score, ]")?,
//! );
//!
//! let mut row = MzTabPSMSectionRow::default();
//! row.sequence.set("NDYKAPPQPAPGK");
//! row.psm_id.set(38);
//! row.search_engine_score
//!     .insert(1, MzTabDouble::new(51.9678841193106));
//! row.calc_mass_to_charge = MzTabDouble::nan();
//! document.psm_data.push(row);
//!
//! let adapter = MzTabFile::new();
//! let text = adapter.write_to_string(&document)?;
//! assert!(text.contains("PSH\tsequence\tPSM_ID"));
//!
//! let reloaded = adapter.load_str(&text)?;
//! assert_eq!(reloaded.psm_data, document.psm_data);
//! assert_eq!(reloaded.meta_data, document.meta_data);
//!
//! // The writer puts a blank line before each section and the reader records
//! // it, so a document built in memory gains one empty row per section on its
//! // first write. From then on the round trip is a fixed point.
//! assert_eq!(reloaded.empty_rows.len(), 1);
//! assert_eq!(adapter.write_to_string(&reloaded)?, text);
//! assert_eq!(adapter.load_str(&adapter.write_to_string(&reloaded)?)?, reloaded);
//! # Ok::<(), openms::Error>(())
//! ```

use crate::format::file_types::{self, FileType};
use crate::format::mztab::{
    MzTab, MzTabCell, MzTabDouble, MzTabMetaData, MzTabNucleicAcidSectionRow, MzTabOSMSectionRow,
    MzTabOligonucleotideSectionRow, MzTabOptionalColumnEntry, MzTabPSMSectionRow, MzTabParameter,
    MzTabParameterList, MzTabPeptideSectionRow, MzTabProteinSectionRow,
    MzTabSmallMoleculeSectionRow, MzTabString, optional_column_names,
};
use crate::format::path_io;
use crate::format::text::{Limits, TextFile};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::io::BufRead;
use std::path::Path;

/// Comment line tag.
const COMMENT: &str = "COM";
/// Metadata line tag.
const METADATA: &str = "MTD";
/// Largest number of diagnostics [`MzTabFile::load_reporting`] collects.
const MAX_DIAGNOSTICS: usize = 1024;

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

fn parse_error(line: usize, message: impl Into<String>) -> Error {
    Error::Parse {
        line,
        message: message.into(),
    }
}

/// `StringUtils::trim`: ASCII space, tab, CR and LF only, matched by characters
/// so a multi-byte character can never be split.
fn trim_source(text: &str) -> &str {
    text.trim_matches(|c| c == ' ' || c == '\t' || c == '\n' || c == '\r')
}

fn has_separator(text: &str) -> bool {
    text.contains(['\t', '\n', '\r'])
}

fn column_count(line: &str) -> usize {
    line.matches('\t').count().saturating_add(1)
}

// ---------------------------------------------------------------------------
// Bracket indices
// ---------------------------------------------------------------------------

/// Extract the bracketed index of a key such as `assay[3]`, given `assay[`.
///
/// Source: the file-static `extractBracketIndex` of `MzTabFile.cpp`, which
/// removes every occurrence of `prefix`, then every `]`, then trims and calls
/// `StringUtils::toInt32`. An empty `prefix` extracts the index from bare
/// digits.
///
/// # Errors
///
/// [`Error::Parse`] when what remains is not a decimal integer, when it is zero
/// or negative, or when it exceeds [`MzTabFile::MAX_INDEX`].
///
/// # Notes
///
/// The source returns a signed `Int` and every caller casts it to `Size`, so
/// `assay[0]` yields key `0` and `assay[-1]` yields key
/// `18446744073709551615`; both then appear in a written header as a column name
/// the format does not define. This port rejects them, because the MzTab
/// specification numbers every indexed key from one and an unbounded key would
/// make the writer's column layout unbounded too.
pub fn extract_bracket_index(key: &str, prefix: &str) -> Result<usize> {
    let stripped = if prefix.is_empty() {
        key.to_owned()
    } else {
        key.replace(prefix, "")
    };
    let digits = stripped.replace(']', "");
    let token = trim_source(&digits);
    let value: i64 = token
        .strip_prefix('+')
        .unwrap_or(token)
        .parse()
        .map_err(|_| {
            parse_error(
                0,
                format!("could not convert MzTab index {token:?} to an integer value"),
            )
        })?;
    if value < 1 {
        return Err(parse_error(
            0,
            format!("MzTab index {value} is not a positive one-based index"),
        ));
    }
    let value = usize::try_from(value).map_err(|_| parse_error(0, "MzTab index out of range"))?;
    if value > MzTabFile::MAX_INDEX {
        return Err(parse_error(
            0,
            format!("MzTab index {value} exceeds the supported maximum"),
        ));
    }
    Ok(value)
}

/// Extract the first two bracketed integers of a column name, as in
/// `search_engine_score[1]_ms_run[2]` yielding `(1, 2)`.
///
/// Source `MzTabFile::extractIndexPairsFromBrackets_`, which matches
/// `^.*?\[(\d+)\].*$` and `^.*?\[\d+\].*?\[(\d+)\].*$` and leaves either member
/// at `0` when its pattern does not match.
///
/// # Errors
///
/// [`Error::Parse`] when either bracketed group is missing, is not a positive
/// one-based index, or exceeds [`MzTabFile::MAX_INDEX`]. The source silently
/// substitutes `0`, producing the family member `…[0]`, which the format does
/// not define; refusing keeps a malformed header visible.
pub fn extract_index_pairs_from_brackets(name: &str) -> Result<(usize, usize)> {
    let mut found: Vec<&str> = Vec::new();
    let mut rest = name;
    while let Some(open) = rest.find('[') {
        // `open` indexes the ASCII '[', so the byte after it is a boundary.
        let after = &rest[open + 1..];
        let Some(close) = after.find(']') else { break };
        let digits = &after[..close];
        if !digits.is_empty() && digits.bytes().all(|byte| byte.is_ascii_digit()) {
            found.push(digits);
            if found.len() == 2 {
                break;
            }
        }
        rest = &after[close + 1..];
    }
    if found.len() < 2 {
        return Err(parse_error(
            0,
            format!("MzTab column {name:?} does not carry two bracketed indices"),
        ));
    }
    Ok((
        extract_bracket_index(found[0], "")?,
        extract_bracket_index(found[1], "")?,
    ))
}

// ---------------------------------------------------------------------------
// Column layout
// ---------------------------------------------------------------------------

type ScoreRunMap = BTreeMap<usize, BTreeMap<usize, MzTabDouble>>;

/// The bracketed column families one section's header declares, plus its
/// optional column names.
///
/// The source derives these numbers inconsistently — the protein header reads
/// them from `MzTabMetaData`, the peptide and small-molecule headers from the
/// *first* row of the section, and the row writers then emit whatever their own
/// maps happen to contain. Where those disagree the source throws
/// `Exception::Postcondition` ("Header and content differs in columns. Please
/// report this bug to the OpenMS developers.") from
/// `MzTabFile::generateMzTabSection_`. This port computes one layout per section
/// as the union of every row's keys with the metadata keys the source consults,
/// and renders every row against that layout with a `null` cell for a member the
/// row does not carry, so header and content agree by construction and no row
/// can be lost to a postcondition.
///
/// The three `count_*` sets are `num_psms_ms_run`,
/// `num_peptides_distinct_ms_run` and `num_peptides_unique_ms_run` in the
/// protein section, and `num_osms_ms_run`, `num_oligos_distinct_ms_run` and
/// `num_oligos_unique_ms_run` in the nucleic-acid section. No other section uses
/// them.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SectionLayout {
    /// Members of the `best_search_engine_score[n]` family.
    pub best_search_engine_score: BTreeSet<usize>,
    /// Score indices: the flat `search_engine_score[n]` family in the `PSM` and
    /// `OSM` sections, and the first index of `search_engine_score[s]_ms_run[r]`
    /// everywhere else.
    pub search_engine_score: BTreeSet<usize>,
    /// MS-run indices: the second index of `search_engine_score[s]_ms_run[r]`.
    pub ms_runs: BTreeSet<usize>,
    /// `num_psms_ms_run[n]` in `PRT`, `num_osms_ms_run[n]` in `NUC`.
    pub count_matches_ms_run: BTreeSet<usize>,
    /// `num_peptides_distinct_ms_run[n]` in `PRT`,
    /// `num_oligos_distinct_ms_run[n]` in `NUC`.
    pub count_distinct_ms_run: BTreeSet<usize>,
    /// `num_peptides_unique_ms_run[n]` in `PRT`, `num_oligos_unique_ms_run[n]`
    /// in `NUC`.
    pub count_unique_ms_run: BTreeSet<usize>,
    /// `<section>_abundance_assay[n]`.
    pub abundance_assay: BTreeSet<usize>,
    /// `<section>_abundance_study_variable[n]`, which also fixes the `stdev` and
    /// `std_error` columns because the three always travel as a triple.
    pub abundance_study_variable: BTreeSet<usize>,
    /// Optional column names, in the order the header declares them.
    pub optional_columns: Vec<String>,
}

fn keys_of<R, V, F>(rows: &[R], project: F) -> BTreeSet<usize>
where
    F: Fn(&R) -> &BTreeMap<usize, V>,
{
    let mut set = BTreeSet::new();
    for row in rows {
        set.extend(project(row).keys().copied());
    }
    set
}

fn score_run_keys<R, F>(rows: &[R], project: F) -> (BTreeSet<usize>, BTreeSet<usize>)
where
    F: Fn(&R) -> &ScoreRunMap,
{
    let mut scores = BTreeSet::new();
    let mut runs = BTreeSet::new();
    for row in rows {
        for (score, per_run) in project(row) {
            scores.insert(*score);
            runs.extend(per_run.keys().copied());
        }
    }
    (scores, runs)
}

impl SectionLayout {
    fn checked(self) -> Result<Self> {
        let product = self
            .search_engine_score
            .len()
            .checked_mul(self.ms_runs.len())
            .ok_or_else(|| bad("MzTab score column count overflows"))?;
        let total = [
            self.best_search_engine_score.len(),
            self.search_engine_score.len(),
            product,
            self.count_matches_ms_run.len(),
            self.count_distinct_ms_run.len(),
            self.count_unique_ms_run.len(),
            self.abundance_assay.len(),
            self.abundance_study_variable.len().saturating_mul(3),
            self.optional_columns.len(),
        ]
        .into_iter()
        .try_fold(0usize, usize::checked_add)
        .ok_or_else(|| bad("MzTab column count overflows"))?;
        if total > MzTabFile::MAX_COLUMNS {
            return Err(bad("MzTab section exceeds its column limit"));
        }
        Ok(self)
    }

    /// In `Complete` mode the source declares a score column for every
    /// `ms_run[n]` of the metadata whether or not a row carries one.
    fn with_complete_runs(mut self, meta: &MzTabMetaData) -> Self {
        if !self.search_engine_score.is_empty() && meta.mz_tab_mode.get() == "Complete" {
            self.ms_runs.extend(meta.ms_run.keys().copied());
        }
        self
    }

    /// Layout of a `PRT` section.
    ///
    /// `best_search_engine_score` is the union of the rows' keys with the
    /// `protein_search_engine_score[n]` keys of `meta`, and the abundance
    /// families are the union of the rows' keys with the `assay[n]` and
    /// `study_variable[n]` keys; both mirror what the source's
    /// `generateMzTabProteinHeader_` reads out of the metadata, except that the
    /// source uses each map's *size* and numbers the columns from one, so a
    /// metadata that declares `assay[2]` and `assay[5]` makes it write
    /// `protein_abundance_assay[1]` and `[2]`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the section exceeds
    /// [`MzTab::MAX_ROWS`](crate::format::mztab::MzTab::MAX_ROWS), its distinct
    /// optional columns exceed
    /// [`MzTab::MAX_OPTIONAL_COLUMNS`](crate::format::mztab::MzTab::MAX_OPTIONAL_COLUMNS),
    /// or the resulting column count exceeds [`MzTabFile::MAX_COLUMNS`].
    pub fn for_protein(rows: &[MzTabProteinSectionRow], meta: &MzTabMetaData) -> Result<Self> {
        let (scores, runs) = score_run_keys(rows, |row| &row.search_engine_score_ms_run);
        let mut best = keys_of(rows, |row| &row.best_search_engine_score);
        best.extend(meta.protein_search_engine_score.keys().copied());
        let mut assay = keys_of(rows, |row| &row.protein_abundance_assay);
        assay.extend(meta.assay.keys().copied());
        let mut study = keys_of(rows, |row| &row.protein_abundance_study_variable);
        study.extend(keys_of(rows, |row| {
            &row.protein_abundance_stdev_study_variable
        }));
        study.extend(keys_of(rows, |row| {
            &row.protein_abundance_std_error_study_variable
        }));
        study.extend(meta.study_variable.keys().copied());
        Self {
            best_search_engine_score: best,
            search_engine_score: scores,
            ms_runs: runs,
            count_matches_ms_run: keys_of(rows, |row| &row.num_psms_ms_run),
            count_distinct_ms_run: keys_of(rows, |row| &row.num_peptides_distinct_ms_run),
            count_unique_ms_run: keys_of(rows, |row| &row.num_peptides_unique_ms_run),
            abundance_assay: assay,
            abundance_study_variable: study,
            optional_columns: optional_column_names(rows)?,
        }
        .with_complete_runs(meta)
        .checked()
    }

    /// Layout of a `PEP` section.
    ///
    /// # Errors
    ///
    /// As [`SectionLayout::for_protein`].
    pub fn for_peptide(rows: &[MzTabPeptideSectionRow], meta: &MzTabMetaData) -> Result<Self> {
        let (scores, runs) = score_run_keys(rows, |row| &row.search_engine_score_ms_run);
        let mut best = keys_of(rows, |row| &row.best_search_engine_score);
        best.extend(meta.peptide_search_engine_score.keys().copied());
        let mut assay = keys_of(rows, |row| &row.peptide_abundance_assay);
        assay.extend(meta.assay.keys().copied());
        let mut study = keys_of(rows, |row| &row.peptide_abundance_study_variable);
        study.extend(keys_of(rows, |row| {
            &row.peptide_abundance_stdev_study_variable
        }));
        study.extend(keys_of(rows, |row| {
            &row.peptide_abundance_std_error_study_variable
        }));
        study.extend(meta.study_variable.keys().copied());
        Self {
            best_search_engine_score: best,
            search_engine_score: scores,
            ms_runs: runs,
            abundance_assay: assay,
            abundance_study_variable: study,
            optional_columns: optional_column_names(rows)?,
            ..Self::default()
        }
        .with_complete_runs(meta)
        .checked()
    }

    /// Layout of a `PSM` section, whose only indexed family is the flat
    /// `search_engine_score[n]`.
    ///
    /// The index set is the union of the rows' keys with the
    /// `psm_search_engine_score[n]` keys of `meta`. The source caps the count at
    /// one, with the comment "we currently only store one search engine score
    /// per PSM"; that limit belongs to its `ConsensusMap` and identification
    /// exporters, not to the file format, so it is not applied here.
    ///
    /// # Errors
    ///
    /// As [`SectionLayout::for_protein`].
    pub fn for_psm(rows: &[MzTabPSMSectionRow], meta: &MzTabMetaData) -> Result<Self> {
        let mut scores = keys_of(rows, |row| &row.search_engine_score);
        scores.extend(meta.psm_search_engine_score.keys().copied());
        Self {
            search_engine_score: scores,
            optional_columns: optional_column_names(rows)?,
            ..Self::default()
        }
        .checked()
    }

    /// Layout of an `SML` section.
    ///
    /// # Errors
    ///
    /// As [`SectionLayout::for_protein`].
    pub fn for_small_molecule(
        rows: &[MzTabSmallMoleculeSectionRow],
        meta: &MzTabMetaData,
    ) -> Result<Self> {
        let (scores, runs) = score_run_keys(rows, |row| &row.search_engine_score_ms_run);
        let mut best = keys_of(rows, |row| &row.best_search_engine_score);
        best.extend(meta.smallmolecule_search_engine_score.keys().copied());
        let mut assay = keys_of(rows, |row| &row.smallmolecule_abundance_assay);
        assay.extend(meta.assay.keys().copied());
        let mut study = keys_of(rows, |row| &row.smallmolecule_abundance_study_variable);
        study.extend(keys_of(rows, |row| {
            &row.smallmolecule_abundance_stdev_study_variable
        }));
        study.extend(keys_of(rows, |row| {
            &row.smallmolecule_abundance_std_error_study_variable
        }));
        study.extend(meta.study_variable.keys().copied());
        Self {
            best_search_engine_score: best,
            search_engine_score: scores,
            ms_runs: runs,
            abundance_assay: assay,
            abundance_study_variable: study,
            optional_columns: optional_column_names(rows)?,
            ..Self::default()
        }
        .with_complete_runs(meta)
        .checked()
    }

    /// Layout of a `NUC` section.
    ///
    /// # Errors
    ///
    /// As [`SectionLayout::for_protein`].
    pub fn for_nucleic_acid(
        rows: &[MzTabNucleicAcidSectionRow],
        meta: &MzTabMetaData,
    ) -> Result<Self> {
        let (scores, runs) = score_run_keys(rows, |row| &row.search_engine_score_ms_run);
        let mut best = keys_of(rows, |row| &row.best_search_engine_score);
        best.extend(meta.nucleic_acid_search_engine_score.keys().copied());
        Self {
            best_search_engine_score: best,
            search_engine_score: scores,
            ms_runs: runs,
            count_matches_ms_run: keys_of(rows, |row| &row.num_osms_ms_run),
            count_distinct_ms_run: keys_of(rows, |row| &row.num_oligos_distinct_ms_run),
            count_unique_ms_run: keys_of(rows, |row| &row.num_oligos_unique_ms_run),
            optional_columns: optional_column_names(rows)?,
            ..Self::default()
        }
        .with_complete_runs(meta)
        .checked()
    }

    /// Layout of an `OLI` section.
    ///
    /// # Errors
    ///
    /// As [`SectionLayout::for_protein`].
    pub fn for_oligonucleotide(
        rows: &[MzTabOligonucleotideSectionRow],
        meta: &MzTabMetaData,
    ) -> Result<Self> {
        let (scores, runs) = score_run_keys(rows, |row| &row.search_engine_score_ms_run);
        let mut best = keys_of(rows, |row| &row.best_search_engine_score);
        best.extend(meta.oligonucleotide_search_engine_score.keys().copied());
        Self {
            best_search_engine_score: best,
            search_engine_score: scores,
            ms_runs: runs,
            optional_columns: optional_column_names(rows)?,
            ..Self::default()
        }
        .with_complete_runs(meta)
        .checked()
    }

    /// Layout of an `OSM` section, whose only indexed family is the flat
    /// `search_engine_score[n]`.
    ///
    /// # Errors
    ///
    /// As [`SectionLayout::for_protein`].
    pub fn for_osm(rows: &[MzTabOSMSectionRow], meta: &MzTabMetaData) -> Result<Self> {
        let mut scores = keys_of(rows, |row| &row.search_engine_score);
        scores.extend(meta.osm_search_engine_score.keys().copied());
        Self {
            search_engine_score: scores,
            optional_columns: optional_column_names(rows)?,
            ..Self::default()
        }
        .checked()
    }
}

// ---------------------------------------------------------------------------
// Rendering helpers
// ---------------------------------------------------------------------------

/// Render the optional columns `names` declares, taking each value from
/// `entries` and writing `null` for a name the row does not carry.
///
/// Source `MzTabFile::addOptionalColumnsToSectionRow_`: a linear search per
/// name, first match wins, and a miss becomes `MzTabString("null")` — which is
/// the null cell, because `MzTabString::set` stores nothing for the literal text
/// `null`.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `names` exceeds
/// [`MzTab::MAX_OPTIONAL_COLUMNS`](crate::format::mztab::MzTab::MAX_OPTIONAL_COLUMNS).
pub fn optional_column_cells(
    names: &[String],
    entries: &[MzTabOptionalColumnEntry],
) -> Result<Vec<String>> {
    if names.len() > MzTab::MAX_OPTIONAL_COLUMNS {
        return Err(bad("MzTab section exceeds its optional-column limit"));
    }
    let mut cells = Vec::with_capacity(names.len());
    for name in names {
        match entries.iter().find(|entry| &entry.name == name) {
            Some(entry) => cells.push(entry.value.to_cell_string()),
            None => cells.push(MzTabString::null().to_cell_string()),
        }
    }
    Ok(cells)
}

fn push_indexed(cells: &mut Vec<String>, indices: &BTreeSet<usize>, label: &str) {
    for index in indices {
        cells.push(format!("{label}[{index}]"));
    }
}

fn push_score_run_names(cells: &mut Vec<String>, layout: &SectionLayout) {
    for score in &layout.search_engine_score {
        for run in &layout.ms_runs {
            cells.push(format!("search_engine_score[{score}]_ms_run[{run}]"));
        }
    }
}

fn push_study_variable_names(cells: &mut Vec<String>, layout: &SectionLayout, prefix: &str) {
    for index in &layout.abundance_study_variable {
        cells.push(format!("{prefix}_abundance_study_variable[{index}]"));
        cells.push(format!("{prefix}_abundance_stdev_study_variable[{index}]"));
        cells.push(format!(
            "{prefix}_abundance_std_error_study_variable[{index}]"
        ));
    }
}

fn cell_or_null<T: MzTabCell>(map: &BTreeMap<usize, T>, index: usize) -> Result<String> {
    match map.get(&index) {
        Some(value) => value.write_cell(),
        None => Ok(MzTabString::null().to_cell_string()),
    }
}

fn push_indexed_cells<T: MzTabCell>(
    cells: &mut Vec<String>,
    indices: &BTreeSet<usize>,
    map: &BTreeMap<usize, T>,
) -> Result<()> {
    for index in indices {
        cells.push(cell_or_null(map, *index)?);
    }
    Ok(())
}

fn push_score_run_cells(
    cells: &mut Vec<String>,
    layout: &SectionLayout,
    map: &ScoreRunMap,
) -> Result<()> {
    for score in &layout.search_engine_score {
        for run in &layout.ms_runs {
            match map.get(score).and_then(|per_run| per_run.get(run)) {
                Some(value) => cells.push(value.to_cell_string()),
                None => cells.push(MzTabString::null().to_cell_string()),
            }
        }
    }
    Ok(())
}

fn push_study_variable_cells(
    cells: &mut Vec<String>,
    layout: &SectionLayout,
    value: &BTreeMap<usize, MzTabDouble>,
    stdev: &BTreeMap<usize, MzTabDouble>,
    std_error: &BTreeMap<usize, MzTabDouble>,
) -> Result<()> {
    for index in &layout.abundance_study_variable {
        cells.push(cell_or_null(value, *index)?);
        cells.push(cell_or_null(stdev, *index)?);
        cells.push(cell_or_null(std_error, *index)?);
    }
    Ok(())
}

/// Join one line's cells, refusing any cell that carries a tab or a line break.
///
/// Neither can be written: a tab would add a column the header does not declare
/// and a line break would split the row. The source concatenates the cells
/// unconditionally and produces a file it cannot read back.
fn join(cells: Vec<String>) -> Result<String> {
    if cells.iter().any(|cell| has_separator(cell)) {
        return Err(bad(
            "MzTab cell contains a tab or line break and cannot be written",
        ));
    }
    Ok(cells.join("\t"))
}

// ---------------------------------------------------------------------------
// The adapter
// ---------------------------------------------------------------------------

/// File adapter for MzTab files.
///
/// Source `MzTabFile`. The nine `store…Column` setters the header declares, and
/// the seven further flags it keeps `protected` with no setter at all, are
/// public fields here, because none of them validates anything: each only
/// decides whether one optional column is written. All are `false` by default,
/// as the source's constructor leaves them. They affect writing only; reading
/// accepts a column whether or not its flag is set.
///
/// The `reliability`, `uri` and `go_terms` columns are optional in the MzTab
/// specification, so suppressing one discards the corresponding cells. A
/// document written with a flag cleared therefore does not round-trip a row
/// that carries that cell; [`MzTabFile::lossless`] enables all sixteen.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabFile {
    /// Write the `PRT` section's `reliability` column. Source
    /// `storeProteinReliabilityColumn`.
    pub store_protein_reliability: bool,
    /// Write the `PEP` section's `reliability` column. Source
    /// `storePeptideReliabilityColumn`.
    pub store_peptide_reliability: bool,
    /// Write the `PSM` section's `reliability` column. Source
    /// `storePSMReliabilityColumn`.
    pub store_psm_reliability: bool,
    /// Write the `SML` section's `reliability` column. Source
    /// `storeSmallMoleculeReliabilityColumn`.
    pub store_small_molecule_reliability: bool,
    /// Write the `PRT` section's `uri` column. Source `storeProteinUriColumn`.
    pub store_protein_uri: bool,
    /// Write the `PEP` section's `uri` column. Source `storePeptideUriColumn`.
    pub store_peptide_uri: bool,
    /// Write the `PSM` section's `uri` column. Source `storePSMUriColumn`.
    pub store_psm_uri: bool,
    /// Write the `SML` section's `uri` column. Source
    /// `storeSmallMoleculeUriColumn`.
    pub store_small_molecule_uri: bool,
    /// Write the `PRT` section's `go_terms` column. Source
    /// `storeProteinGoTerms`.
    pub store_protein_go_terms: bool,
    /// Write the `NUC` section's `reliability` column. Source
    /// `store_nucleic_acid_reliability_`, which has no public setter.
    pub store_nucleic_acid_reliability: bool,
    /// Write the `OLI` section's `reliability` column. Source
    /// `store_oligonucleotide_reliability_`, which has no public setter.
    pub store_oligonucleotide_reliability: bool,
    /// Write the `OSM` section's `reliability` column. Source
    /// `store_osm_reliability_`, which has no public setter.
    pub store_osm_reliability: bool,
    /// Write the `NUC` section's `uri` column. Source `store_nucleic_acid_uri_`,
    /// which has no public setter.
    pub store_nucleic_acid_uri: bool,
    /// Write the `OLI` section's `uri` column. Source
    /// `store_oligonucleotide_uri_`, which has no public setter.
    pub store_oligonucleotide_uri: bool,
    /// Write the `OSM` section's `uri` column. Source `store_osm_uri_`, which
    /// has no public setter.
    pub store_osm_uri: bool,
    /// Write the `NUC` section's `go_terms` column. Source
    /// `store_nucleic_acid_goterms_`, which has no public setter.
    pub store_nucleic_acid_go_terms: bool,
}

impl MzTabFile {
    /// Largest number of lines either direction will handle.
    ///
    /// The source streams through `TextFile`, which has no ceiling of its own
    /// and holds every line of the file in memory at once.
    pub const MAX_LINES: usize = 4_000_000;
    /// Largest total payload, in bytes, either direction will handle.
    pub const MAX_BYTES: usize = 512 * 1024 * 1024;
    /// Largest number of columns one row may have.
    pub const MAX_COLUMNS: usize = 200_000;
    /// Largest bracketed index any indexed key or column name may carry.
    pub const MAX_INDEX: usize = 1_000_000;

    /// A fresh adapter with every optional column disabled, as the source's
    /// default constructor.
    pub fn new() -> Self {
        Self::default()
    }

    /// Enable every optional `reliability`, `uri` and `go_terms` column.
    ///
    /// Native convenience: with all sixteen flags set, writing preserves every
    /// cell of every row, which is what a lossless round trip needs. The source
    /// offers no such shortcut and leaves seven of the flags unreachable.
    pub fn lossless() -> Self {
        Self {
            store_protein_reliability: true,
            store_peptide_reliability: true,
            store_psm_reliability: true,
            store_small_molecule_reliability: true,
            store_protein_uri: true,
            store_peptide_uri: true,
            store_psm_uri: true,
            store_small_molecule_uri: true,
            store_protein_go_terms: true,
            store_nucleic_acid_reliability: true,
            store_oligonucleotide_reliability: true,
            store_osm_reliability: true,
            store_nucleic_acid_uri: true,
            store_oligonucleotide_uri: true,
            store_osm_uri: true,
            store_nucleic_acid_go_terms: true,
        }
    }
}

/// A bounded accumulator for output lines.
#[derive(Debug, Default)]
struct LineBuffer {
    lines: Vec<String>,
    bytes: usize,
}

impl LineBuffer {
    fn push(&mut self, line: String) -> Result<()> {
        if line.contains(['\n', '\r']) {
            return Err(bad("MzTab output line contains a line break"));
        }
        if self.lines.len() >= MzTabFile::MAX_LINES {
            return Err(bad("MzTab output exceeds its line limit"));
        }
        self.bytes = self
            .bytes
            .checked_add(line.len().saturating_add(1))
            .filter(|&total| total <= MzTabFile::MAX_BYTES)
            .ok_or_else(|| bad("MzTab output exceeds its byte limit"))?;
        self.lines
            .try_reserve(1)
            .map_err(|_| bad("MzTab output allocation failed"))?;
        self.lines.push(line);
        Ok(())
    }

    fn extend(&mut self, lines: Vec<String>) -> Result<()> {
        for line in lines {
            self.push(line)?;
        }
        Ok(())
    }

    /// Push one `MTD` line, refusing a key or value that carries a tab or a
    /// line break: an embedded tab would make the value a fourth cell the
    /// reader ignores, silently truncating it on the next round trip.
    fn metadata(&mut self, key: &str, value: &str) -> Result<()> {
        if has_separator(key) || has_separator(value) {
            return Err(bad(
                "MzTab metadata key or value contains a tab or line break",
            ));
        }
        self.push(format!("{METADATA}\t{key}\t{value}"))
    }
}

fn reference_list(label: &str, values: &[i32]) -> String {
    values
        .iter()
        .map(|value| format!("{label}[{value}]"))
        .collect::<Vec<_>>()
        .join(",")
}

// ---------------------------------------------------------------------------
// Writing — metadata
// ---------------------------------------------------------------------------

impl MzTabFile {
    /// Render the `MTD` section of `meta`.
    ///
    /// Source `MzTabFile::generateMzTabMetaDataSection_`, key for key and in the
    /// same order: `mzTab-version`, `mzTab-mode` and `mzTab-type`
    /// unconditionally — a null cell becomes the text `null` — then `title` and
    /// `mzTab-ID` only when set, then `description` unconditionally, then every
    /// indexed family. A `software[n]` line is written even when the term is
    /// null, as the source does; every other optional key is skipped when null.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the section exceeds [`MzTabFile::MAX_LINES`]
    /// or [`MzTabFile::MAX_BYTES`], or when a key or value carries a tab or a
    /// line break.
    ///
    /// # Notes
    ///
    /// Two keys diverge from the source, because the source cannot read back
    /// what it writes for them.
    ///
    /// `colunit-protein`, `colunit-peptide`, `colunit-psm` and
    /// `colunit-small_molecule` are written with a tab between the key and the
    /// value. The source concatenates the two with no separator at all
    /// (`MzTabFile.cpp:1984`), so its own reader sees a two-cell line; and even
    /// a correctly separated line makes its reader throw, because it feeds the
    /// literal key `colunit` to `StringUtils::toInt32` while looking for a
    /// bracketed index that a `colunit` key never has (`MzTabFile.cpp:685`). The
    /// source also writes the PSM key as `colunit-PSM` while its reader matches
    /// the lower-case `psm`; this port writes the specification's
    /// `colunit-psm` and accepts either case on read.
    ///
    /// `nucleic_acid_search_engine_score[n]`,
    /// `oligonucleotide_search_engine_score[n]` and `osm_search_engine_score[n]`
    /// are written by the source but never read by it;
    /// [`MzTabFile::load`] reads all three.
    pub fn metadata_section_lines(&self, meta: &MzTabMetaData) -> Result<Vec<String>> {
        let mut out = LineBuffer::default();
        out.metadata("mzTab-version", &meta.mz_tab_version.to_cell_string())?;
        out.metadata("mzTab-mode", &meta.mz_tab_mode.to_cell_string())?;
        out.metadata("mzTab-type", &meta.mz_tab_type.to_cell_string())?;
        if !meta.title.is_null() {
            out.metadata("title", &meta.title.to_cell_string())?;
        }
        if !meta.mz_tab_id.is_null() {
            out.metadata("mzTab-ID", &meta.mz_tab_id.to_cell_string())?;
        }
        out.metadata("description", &meta.description.to_cell_string())?;
        for (index, value) in &meta.sample_processing {
            out.metadata(
                &format!("sample_processing[{index}]"),
                &value.to_cell_string(),
            )?;
        }
        for (label, scores) in [
            ("protein", &meta.protein_search_engine_score),
            ("peptide", &meta.peptide_search_engine_score),
            ("psm", &meta.psm_search_engine_score),
            ("smallmolecule", &meta.smallmolecule_search_engine_score),
            ("nucleic_acid", &meta.nucleic_acid_search_engine_score),
            ("oligonucleotide", &meta.oligonucleotide_search_engine_score),
            ("osm", &meta.osm_search_engine_score),
        ] {
            for (index, value) in scores {
                out.metadata(
                    &format!("{label}_search_engine_score[{index}]"),
                    &value.to_cell_string(),
                )?;
            }
        }
        for (index, instrument) in &meta.instrument {
            if !instrument.name.is_null() {
                out.metadata(
                    &format!("instrument[{index}]-name"),
                    &instrument.name.to_cell_string(),
                )?;
            }
            if !instrument.source.is_null() {
                out.metadata(
                    &format!("instrument[{index}]-source"),
                    &instrument.source.to_cell_string(),
                )?;
            }
            for (analyzer, value) in &instrument.analyzer {
                if !value.is_null() {
                    out.metadata(
                        &format!("instrument[{index}]-analyzer[{analyzer}]"),
                        &value.to_cell_string(),
                    )?;
                }
            }
            if !instrument.detector.is_null() {
                out.metadata(
                    &format!("instrument[{index}]-detector"),
                    &instrument.detector.to_cell_string(),
                )?;
            }
        }
        for (index, software) in &meta.software {
            out.metadata(
                &format!("software[{index}]"),
                &software.software.to_cell_string(),
            )?;
            for (setting, value) in &software.setting {
                out.metadata(
                    &format!("software[{index}]-setting[{setting}]"),
                    &value.to_cell_string(),
                )?;
            }
        }
        if !meta.false_discovery_rate.is_null() {
            out.metadata(
                "false_discovery_rate",
                &meta.false_discovery_rate.to_cell_string(),
            )?;
        }
        for (index, value) in &meta.publication {
            out.metadata(&format!("publication[{index}]"), &value.to_cell_string())?;
        }
        for (index, contact) in &meta.contact {
            for (suffix, value) in [
                ("name", &contact.name),
                ("affiliation", &contact.affiliation),
                ("email", &contact.email),
            ] {
                if !value.is_null() {
                    out.metadata(
                        &format!("contact[{index}]-{suffix}"),
                        &value.to_cell_string(),
                    )?;
                }
            }
        }
        for (index, value) in &meta.uri {
            out.metadata(&format!("uri[{index}]"), &value.to_cell_string())?;
        }
        for (label, mods) in [
            ("fixed_mod", &meta.fixed_mod),
            ("variable_mod", &meta.variable_mod),
        ] {
            for (index, entry) in mods {
                if !entry.modification.is_null() {
                    out.metadata(
                        &format!("{label}[{index}]"),
                        &entry.modification.to_cell_string(),
                    )?;
                }
                if !entry.site.is_null() {
                    out.metadata(
                        &format!("{label}[{index}]-site"),
                        &entry.site.to_cell_string(),
                    )?;
                }
                if !entry.position.is_null() {
                    out.metadata(
                        &format!("{label}[{index}]-position"),
                        &entry.position.to_cell_string(),
                    )?;
                }
            }
        }
        if !meta.quantification_method.is_null() {
            out.metadata(
                "quantification_method",
                &meta.quantification_method.to_cell_string(),
            )?;
        }
        for (key, value) in [
            (
                "protein-quantification_unit",
                &meta.protein_quantification_unit,
            ),
            (
                "peptide-quantification_unit",
                &meta.peptide_quantification_unit,
            ),
            (
                "small_molecule-quantification_unit",
                &meta.small_molecule_quantification_unit,
            ),
        ] {
            if !value.is_null() {
                out.metadata(key, &value.to_cell_string())?;
            }
        }
        for (index, run) in &meta.ms_run {
            if !run.format.is_null() {
                out.metadata(
                    &format!("ms_run[{index}]-format"),
                    &run.format.to_cell_string(),
                )?;
            }
            if !run.location.is_null() {
                out.metadata(
                    &format!("ms_run[{index}]-location"),
                    &run.location.to_cell_string(),
                )?;
            }
            if !run.id_format.is_null() {
                out.metadata(
                    &format!("ms_run[{index}]-id_format"),
                    &run.id_format.to_cell_string(),
                )?;
            }
            if !run.fragmentation_method.is_null() {
                out.metadata(
                    &format!("ms_run[{index}]-fragmentation_method"),
                    &run.fragmentation_method.to_cell_string(),
                )?;
            }
        }
        for (index, value) in &meta.custom {
            out.metadata(&format!("custom[{index}]"), &value.to_cell_string())?;
        }
        for (index, sample) in &meta.sample {
            for (label, terms) in [
                ("species", &sample.species),
                ("tissue", &sample.tissue),
                ("cell_type", &sample.cell_type),
                ("disease", &sample.disease),
                ("custom", &sample.custom),
            ] {
                for (inner, value) in terms {
                    out.metadata(
                        &format!("sample[{index}]-{label}[{inner}]"),
                        &value.to_cell_string(),
                    )?;
                }
            }
            if !sample.description.is_null() {
                out.metadata(
                    &format!("sample[{index}]-description"),
                    &sample.description.to_cell_string(),
                )?;
            }
        }
        for (index, assay) in &meta.assay {
            if !assay.quantification_reagent.is_null() {
                out.metadata(
                    &format!("assay[{index}]-quantification_reagent"),
                    &assay.quantification_reagent.to_cell_string(),
                )?;
            }
            for (inner, entry) in &assay.quantification_mod {
                if !entry.modification.is_null() {
                    out.metadata(
                        &format!("assay[{index}]-quantification_mod[{inner}]"),
                        &entry.modification.to_cell_string(),
                    )?;
                }
                if !entry.site.is_null() {
                    out.metadata(
                        &format!("assay[{index}]-quantification_mod[{inner}]-site"),
                        &entry.site.to_cell_string(),
                    )?;
                }
                if !entry.position.is_null() {
                    out.metadata(
                        &format!("assay[{index}]-quantification_mod[{inner}]-position"),
                        &entry.position.to_cell_string(),
                    )?;
                }
            }
            if !assay.sample_ref.is_null() {
                out.metadata(
                    &format!("assay[{index}]-sample_ref"),
                    &assay.sample_ref.to_cell_string(),
                )?;
            }
            if !assay.ms_run_ref.is_empty() {
                out.metadata(
                    &format!("assay[{index}]-ms_run_ref"),
                    &reference_list("ms_run", &assay.ms_run_ref),
                )?;
            }
        }
        for (index, variable) in &meta.study_variable {
            if !variable.assay_refs.is_empty() {
                out.metadata(
                    &format!("study_variable[{index}]-assay_refs"),
                    &reference_list("assay", &variable.assay_refs),
                )?;
            }
            if !variable.sample_refs.is_empty() {
                out.metadata(
                    &format!("study_variable[{index}]-sample_refs"),
                    &reference_list("sample", &variable.sample_refs),
                )?;
            }
            if !variable.description.is_null() {
                out.metadata(
                    &format!("study_variable[{index}]-description"),
                    &variable.description.to_cell_string(),
                )?;
            }
        }
        for (index, cv) in &meta.cv {
            for (suffix, value) in [
                ("label", &cv.label),
                ("full_name", &cv.full_name),
                ("version", &cv.version),
                ("url", &cv.url),
            ] {
                if !value.is_null() {
                    out.metadata(&format!("cv[{index}]-{suffix}"), &value.to_cell_string())?;
                }
            }
        }
        for (key, values) in [
            ("colunit-protein", &meta.colunit_protein),
            ("colunit-peptide", &meta.colunit_peptide),
            ("colunit-psm", &meta.colunit_psm),
            ("colunit-small_molecule", &meta.colunit_small_molecule),
        ] {
            for value in values {
                out.metadata(key, value)?;
            }
        }
        Ok(out.lines)
    }
}

// ---------------------------------------------------------------------------
// Writing — section headers and rows
// ---------------------------------------------------------------------------

impl MzTabFile {
    /// Render the `PRH` header line for a section with this `layout`.
    ///
    /// Source `MzTabFile::generateMzTabProteinHeader_`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when an optional column name carries a tab or a
    /// line break.
    ///
    /// # Notes
    ///
    /// The source emits the `search_engine_score[s]_ms_run[r]` columns with the
    /// runs in the outer loop and the scores in the inner one
    /// (`MzTabFile.cpp:2037`) while [`MzTabFile::protein_row`] emits the values
    /// scores-outer, runs-inner (`MzTabFile.cpp:2119`). The two orders coincide
    /// only when there is at most one score type or at most one MS run;
    /// otherwise the counts still match, the postcondition stays silent, and the
    /// values land under the wrong headers. Both sides here use scores-outer,
    /// runs-inner, which agrees with the source on every reference file and
    /// round-trips in general.
    pub fn protein_header(&self, layout: &SectionLayout) -> Result<String> {
        let mut cells = vec![
            "PRH".to_owned(),
            "accession".to_owned(),
            "description".to_owned(),
            "taxid".to_owned(),
            "species".to_owned(),
            "database".to_owned(),
            "database_version".to_owned(),
            "search_engine".to_owned(),
        ];
        push_indexed(
            &mut cells,
            &layout.best_search_engine_score,
            "best_search_engine_score",
        );
        push_score_run_names(&mut cells, layout);
        if self.store_protein_reliability {
            cells.push("reliability".to_owned());
        }
        push_indexed(&mut cells, &layout.count_matches_ms_run, "num_psms_ms_run");
        push_indexed(
            &mut cells,
            &layout.count_distinct_ms_run,
            "num_peptides_distinct_ms_run",
        );
        push_indexed(
            &mut cells,
            &layout.count_unique_ms_run,
            "num_peptides_unique_ms_run",
        );
        cells.push("ambiguity_members".to_owned());
        cells.push("modifications".to_owned());
        if self.store_protein_uri {
            cells.push("uri".to_owned());
        }
        if self.store_protein_go_terms {
            cells.push("go_terms".to_owned());
        }
        cells.push("protein_coverage".to_owned());
        push_indexed(
            &mut cells,
            &layout.abundance_assay,
            "protein_abundance_assay",
        );
        push_study_variable_names(&mut cells, layout, "protein");
        cells.extend(layout.optional_columns.iter().cloned());
        join(cells)
    }

    /// Render one `PRT` row against `layout`.
    ///
    /// Source `MzTabFile::generateMzTabSectionRow_(const MzTabProteinSectionRow&, …)`.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when the row's `modifications` cell carries
    /// positions but no identifier, and [`Error::InvalidValue`] when the layout
    /// declares more optional columns than
    /// [`MzTab::MAX_OPTIONAL_COLUMNS`](crate::format::mztab::MzTab::MAX_OPTIONAL_COLUMNS)
    /// or a cell carries a tab or a line break.
    pub fn protein_row(
        &self,
        row: &MzTabProteinSectionRow,
        layout: &SectionLayout,
    ) -> Result<String> {
        let mut cells = vec![
            "PRT".to_owned(),
            row.accession.to_cell_string(),
            row.description.to_cell_string(),
            row.taxid.to_cell_string(),
            row.species.to_cell_string(),
            row.database.to_cell_string(),
            row.database_version.to_cell_string(),
            row.search_engine.to_cell_string(),
        ];
        push_indexed_cells(
            &mut cells,
            &layout.best_search_engine_score,
            &row.best_search_engine_score,
        )?;
        push_score_run_cells(&mut cells, layout, &row.search_engine_score_ms_run)?;
        if self.store_protein_reliability {
            cells.push(row.reliability.to_cell_string());
        }
        push_indexed_cells(
            &mut cells,
            &layout.count_matches_ms_run,
            &row.num_psms_ms_run,
        )?;
        push_indexed_cells(
            &mut cells,
            &layout.count_distinct_ms_run,
            &row.num_peptides_distinct_ms_run,
        )?;
        push_indexed_cells(
            &mut cells,
            &layout.count_unique_ms_run,
            &row.num_peptides_unique_ms_run,
        )?;
        cells.push(row.ambiguity_members.to_cell_string());
        cells.push(row.modifications.to_cell_string()?);
        if self.store_protein_uri {
            cells.push(row.uri.to_cell_string());
        }
        if self.store_protein_go_terms {
            cells.push(row.go_terms.to_cell_string());
        }
        cells.push(row.coverage.to_cell_string());
        push_indexed_cells(
            &mut cells,
            &layout.abundance_assay,
            &row.protein_abundance_assay,
        )?;
        push_study_variable_cells(
            &mut cells,
            layout,
            &row.protein_abundance_study_variable,
            &row.protein_abundance_stdev_study_variable,
            &row.protein_abundance_std_error_study_variable,
        )?;
        cells.extend(optional_column_cells(&layout.optional_columns, &row.opt)?);
        join(cells)
    }

    /// Render the `PEH` header line for a section with this `layout`.
    ///
    /// Source `MzTabFile::generateMzTabPeptideHeader_`, whose score columns run
    /// runs-outer, scores-inner; see [`MzTabFile::protein_header`] for why this
    /// port orders them the other way round.
    ///
    /// # Errors
    ///
    /// As [`MzTabFile::protein_header`].
    pub fn peptide_header(&self, layout: &SectionLayout) -> Result<String> {
        let mut cells = vec![
            "PEH".to_owned(),
            "sequence".to_owned(),
            "accession".to_owned(),
            "unique".to_owned(),
            "database".to_owned(),
            "database_version".to_owned(),
            "search_engine".to_owned(),
        ];
        push_indexed(
            &mut cells,
            &layout.best_search_engine_score,
            "best_search_engine_score",
        );
        push_score_run_names(&mut cells, layout);
        if self.store_peptide_reliability {
            cells.push("reliability".to_owned());
        }
        cells.push("modifications".to_owned());
        cells.push("retention_time".to_owned());
        cells.push("retention_time_window".to_owned());
        cells.push("charge".to_owned());
        cells.push("mass_to_charge".to_owned());
        if self.store_peptide_uri {
            cells.push("uri".to_owned());
        }
        cells.push("spectra_ref".to_owned());
        push_indexed(
            &mut cells,
            &layout.abundance_assay,
            "peptide_abundance_assay",
        );
        push_study_variable_names(&mut cells, layout, "peptide");
        cells.extend(layout.optional_columns.iter().cloned());
        join(cells)
    }

    /// Render one `PEP` row against `layout`.
    ///
    /// Source `MzTabFile::generateMzTabSectionRow_(const MzTabPeptideSectionRow&, …)`.
    ///
    /// # Errors
    ///
    /// As [`MzTabFile::protein_row`].
    ///
    /// # Notes
    ///
    /// The source advances three iterators in lock-step over
    /// `peptide_abundance_study_variable`, `…_stdev_…` and `…_std_error_…` and
    /// stops as soon as any one of them is exhausted (`MzTabFile.cpp:2367`), so
    /// a row whose `stdev` map is empty writes no study-variable columns at all
    /// and trips the column-count postcondition. Here the triple is driven by
    /// the layout and a missing member becomes a `null` cell.
    pub fn peptide_row(
        &self,
        row: &MzTabPeptideSectionRow,
        layout: &SectionLayout,
    ) -> Result<String> {
        let mut cells = vec![
            "PEP".to_owned(),
            row.sequence.to_cell_string(),
            row.accession.to_cell_string(),
            row.unique.to_cell_string(),
            row.database.to_cell_string(),
            row.database_version.to_cell_string(),
            row.search_engine.to_cell_string(),
        ];
        push_indexed_cells(
            &mut cells,
            &layout.best_search_engine_score,
            &row.best_search_engine_score,
        )?;
        push_score_run_cells(&mut cells, layout, &row.search_engine_score_ms_run)?;
        if self.store_peptide_reliability {
            cells.push(row.reliability.to_cell_string());
        }
        cells.push(row.modifications.to_cell_string()?);
        cells.push(row.retention_time.to_cell_string());
        cells.push(row.retention_time_window.to_cell_string());
        cells.push(row.charge.to_cell_string());
        cells.push(row.mass_to_charge.to_cell_string());
        if self.store_peptide_uri {
            cells.push(row.uri.to_cell_string());
        }
        cells.push(row.spectra_ref.to_cell_string());
        push_indexed_cells(
            &mut cells,
            &layout.abundance_assay,
            &row.peptide_abundance_assay,
        )?;
        push_study_variable_cells(
            &mut cells,
            layout,
            &row.peptide_abundance_study_variable,
            &row.peptide_abundance_stdev_study_variable,
            &row.peptide_abundance_std_error_study_variable,
        )?;
        cells.extend(optional_column_cells(&layout.optional_columns, &row.opt)?);
        join(cells)
    }

    /// Render the `PSH` header line for a section with this `layout`.
    ///
    /// Source `MzTabFile::generateMzTabPSMHeader_`.
    ///
    /// # Errors
    ///
    /// As [`MzTabFile::protein_header`].
    pub fn psm_header(&self, layout: &SectionLayout) -> Result<String> {
        let mut cells = vec![
            "PSH".to_owned(),
            "sequence".to_owned(),
            "PSM_ID".to_owned(),
            "accession".to_owned(),
            "unique".to_owned(),
            "database".to_owned(),
            "database_version".to_owned(),
            "search_engine".to_owned(),
        ];
        push_indexed(
            &mut cells,
            &layout.search_engine_score,
            "search_engine_score",
        );
        if self.store_psm_reliability {
            cells.push("reliability".to_owned());
        }
        cells.push("modifications".to_owned());
        cells.push("retention_time".to_owned());
        cells.push("charge".to_owned());
        cells.push("exp_mass_to_charge".to_owned());
        cells.push("calc_mass_to_charge".to_owned());
        if self.store_psm_uri {
            cells.push("uri".to_owned());
        }
        cells.push("spectra_ref".to_owned());
        cells.push("pre".to_owned());
        cells.push("post".to_owned());
        cells.push("start".to_owned());
        cells.push("end".to_owned());
        cells.extend(layout.optional_columns.iter().cloned());
        join(cells)
    }

    /// Render one `PSM` row against `layout`.
    ///
    /// Source `MzTabFile::generateMzTabSectionRow_(const MzTabPSMSectionRow&, …)`.
    ///
    /// # Errors
    ///
    /// As [`MzTabFile::protein_row`].
    ///
    /// # Notes
    ///
    /// The source writes a single literal `null` score cell when the row's
    /// `search_engine_score` map is empty — a workaround, its comment says, for
    /// peptide identifications without hits in the quality-control export
    /// (`MzTabFile.cpp:2397`) — which matches the header only when the header
    /// declared exactly one score column. Driving the cells from the layout
    /// reproduces that `null` whenever the layout declares a score the row
    /// lacks, and stays consistent when it declares none or several.
    pub fn psm_row(&self, row: &MzTabPSMSectionRow, layout: &SectionLayout) -> Result<String> {
        let mut cells = vec![
            "PSM".to_owned(),
            row.sequence.to_cell_string(),
            row.psm_id.to_cell_string(),
            row.accession.to_cell_string(),
            row.unique.to_cell_string(),
            row.database.to_cell_string(),
            row.database_version.to_cell_string(),
            row.search_engine.to_cell_string(),
        ];
        push_indexed_cells(
            &mut cells,
            &layout.search_engine_score,
            &row.search_engine_score,
        )?;
        if self.store_psm_reliability {
            cells.push(row.reliability.to_cell_string());
        }
        cells.push(row.modifications.to_cell_string()?);
        cells.push(row.retention_time.to_cell_string());
        cells.push(row.charge.to_cell_string());
        cells.push(row.exp_mass_to_charge.to_cell_string());
        cells.push(row.calc_mass_to_charge.to_cell_string());
        if self.store_psm_uri {
            cells.push(row.uri.to_cell_string());
        }
        cells.push(row.spectra_ref.to_cell_string());
        cells.push(row.pre.to_cell_string());
        cells.push(row.post.to_cell_string());
        cells.push(row.start.to_cell_string());
        cells.push(row.end.to_cell_string());
        cells.extend(optional_column_cells(&layout.optional_columns, &row.opt)?);
        join(cells)
    }

    /// Render the `SMH` header line for a section with this `layout`.
    ///
    /// Source `MzTabFile::generateMzTabSmallMoleculeHeader_`.
    ///
    /// # Errors
    ///
    /// As [`MzTabFile::protein_header`].
    pub fn small_molecule_header(&self, layout: &SectionLayout) -> Result<String> {
        let mut cells = vec![
            "SMH".to_owned(),
            "identifier".to_owned(),
            "chemical_formula".to_owned(),
            "smiles".to_owned(),
            "inchi_key".to_owned(),
            "description".to_owned(),
            "exp_mass_to_charge".to_owned(),
            "calc_mass_to_charge".to_owned(),
            "charge".to_owned(),
            "retention_time".to_owned(),
            "taxid".to_owned(),
            "species".to_owned(),
            "database".to_owned(),
            "database_version".to_owned(),
        ];
        if self.store_small_molecule_reliability {
            cells.push("reliability".to_owned());
        }
        if self.store_small_molecule_uri {
            cells.push("uri".to_owned());
        }
        cells.push("spectra_ref".to_owned());
        cells.push("search_engine".to_owned());
        push_indexed(
            &mut cells,
            &layout.best_search_engine_score,
            "best_search_engine_score",
        );
        push_score_run_names(&mut cells, layout);
        cells.push("modifications".to_owned());
        push_indexed(
            &mut cells,
            &layout.abundance_assay,
            "smallmolecule_abundance_assay",
        );
        push_study_variable_names(&mut cells, layout, "smallmolecule");
        cells.extend(layout.optional_columns.iter().cloned());
        join(cells)
    }

    /// Render one `SML` row against `layout`.
    ///
    /// Source `MzTabFile::generateMzTabSectionRow_(const MzTabSmallMoleculeSectionRow&, …)`.
    ///
    /// # Errors
    ///
    /// As [`MzTabFile::protein_row`], except that this section's `modifications`
    /// cell is plain text and cannot fail.
    ///
    /// # Notes
    ///
    /// The source's small-molecule row writer never emits the
    /// `smallmolecule_abundance_assay[n]` cells even though its header declares
    /// one per assay (`MzTabFile.cpp:2548`), so a document with at least one
    /// assay fails the column-count postcondition and cannot be written at all.
    /// This port emits them, as the peptide section does.
    pub fn small_molecule_row(
        &self,
        row: &MzTabSmallMoleculeSectionRow,
        layout: &SectionLayout,
    ) -> Result<String> {
        let mut cells = vec![
            "SML".to_owned(),
            row.identifier.to_cell_string(),
            row.chemical_formula.to_cell_string(),
            row.smiles.to_cell_string(),
            row.inchi_key.to_cell_string(),
            row.description.to_cell_string(),
            row.exp_mass_to_charge.to_cell_string(),
            row.calc_mass_to_charge.to_cell_string(),
            row.charge.to_cell_string(),
            row.retention_time.to_cell_string(),
            row.taxid.to_cell_string(),
            row.species.to_cell_string(),
            row.database.to_cell_string(),
            row.database_version.to_cell_string(),
        ];
        if self.store_small_molecule_reliability {
            cells.push(row.reliability.to_cell_string());
        }
        if self.store_small_molecule_uri {
            cells.push(row.uri.to_cell_string());
        }
        cells.push(row.spectra_ref.to_cell_string());
        cells.push(row.search_engine.to_cell_string());
        push_indexed_cells(
            &mut cells,
            &layout.best_search_engine_score,
            &row.best_search_engine_score,
        )?;
        push_score_run_cells(&mut cells, layout, &row.search_engine_score_ms_run)?;
        cells.push(row.modifications.to_cell_string());
        push_indexed_cells(
            &mut cells,
            &layout.abundance_assay,
            &row.smallmolecule_abundance_assay,
        )?;
        push_study_variable_cells(
            &mut cells,
            layout,
            &row.smallmolecule_abundance_study_variable,
            &row.smallmolecule_abundance_stdev_study_variable,
            &row.smallmolecule_abundance_std_error_study_variable,
        )?;
        cells.extend(optional_column_cells(&layout.optional_columns, &row.opt)?);
        join(cells)
    }

    /// Render the `NUH` header line for a section with this `layout`.
    ///
    /// Source `MzTabFile::generateMzTabNucleicAcidHeader_`.
    ///
    /// # Errors
    ///
    /// As [`MzTabFile::protein_header`].
    ///
    /// # Notes
    ///
    /// The source numbers the `num_osms_ms_run`, `num_oligos_distinct_ms_run`
    /// and `num_oligos_unique_ms_run` columns from *zero*
    /// (`MzTabFile.cpp:2603`), alone among all bracketed families, while its row
    /// writer emits values keyed from one. This port numbers them from the keys
    /// the rows carry. The source also calls this generator with the score and
    /// best-score counts transposed (`MzTabFile.cpp:3257`), which this port does
    /// not.
    pub fn nucleic_acid_header(&self, layout: &SectionLayout) -> Result<String> {
        let mut cells = vec![
            "NUH".to_owned(),
            "accession".to_owned(),
            "description".to_owned(),
            "taxid".to_owned(),
            "species".to_owned(),
            "database".to_owned(),
            "database_version".to_owned(),
            "search_engine".to_owned(),
        ];
        push_indexed(
            &mut cells,
            &layout.best_search_engine_score,
            "best_search_engine_score",
        );
        push_score_run_names(&mut cells, layout);
        if self.store_nucleic_acid_reliability {
            cells.push("reliability".to_owned());
        }
        push_indexed(&mut cells, &layout.count_matches_ms_run, "num_osms_ms_run");
        push_indexed(
            &mut cells,
            &layout.count_distinct_ms_run,
            "num_oligos_distinct_ms_run",
        );
        push_indexed(
            &mut cells,
            &layout.count_unique_ms_run,
            "num_oligos_unique_ms_run",
        );
        cells.push("ambiguity_members".to_owned());
        cells.push("modifications".to_owned());
        if self.store_nucleic_acid_uri {
            cells.push("uri".to_owned());
        }
        if self.store_nucleic_acid_go_terms {
            cells.push("go_terms".to_owned());
        }
        cells.push("sequence_coverage".to_owned());
        cells.extend(layout.optional_columns.iter().cloned());
        join(cells)
    }

    /// Render one `NUC` row against `layout`.
    ///
    /// Source `MzTabFile::generateMzTabSectionRow_(const MzTabNucleicAcidSectionRow&, …)`.
    ///
    /// # Errors
    ///
    /// As [`MzTabFile::protein_row`].
    pub fn nucleic_acid_row(
        &self,
        row: &MzTabNucleicAcidSectionRow,
        layout: &SectionLayout,
    ) -> Result<String> {
        let mut cells = vec![
            "NUC".to_owned(),
            row.accession.to_cell_string(),
            row.description.to_cell_string(),
            row.taxid.to_cell_string(),
            row.species.to_cell_string(),
            row.database.to_cell_string(),
            row.database_version.to_cell_string(),
            row.search_engine.to_cell_string(),
        ];
        push_indexed_cells(
            &mut cells,
            &layout.best_search_engine_score,
            &row.best_search_engine_score,
        )?;
        push_score_run_cells(&mut cells, layout, &row.search_engine_score_ms_run)?;
        if self.store_nucleic_acid_reliability {
            cells.push(row.reliability.to_cell_string());
        }
        push_indexed_cells(
            &mut cells,
            &layout.count_matches_ms_run,
            &row.num_osms_ms_run,
        )?;
        push_indexed_cells(
            &mut cells,
            &layout.count_distinct_ms_run,
            &row.num_oligos_distinct_ms_run,
        )?;
        push_indexed_cells(
            &mut cells,
            &layout.count_unique_ms_run,
            &row.num_oligos_unique_ms_run,
        )?;
        cells.push(row.ambiguity_members.to_cell_string());
        cells.push(row.modifications.to_cell_string()?);
        if self.store_nucleic_acid_uri {
            cells.push(row.uri.to_cell_string());
        }
        if self.store_nucleic_acid_go_terms {
            cells.push(row.go_terms.to_cell_string());
        }
        cells.push(row.coverage.to_cell_string());
        cells.extend(optional_column_cells(&layout.optional_columns, &row.opt)?);
        join(cells)
    }

    /// Render the `OLH` header line for a section with this `layout`.
    ///
    /// Source `MzTabFile::generateMzTabOligonucleotideHeader_`.
    ///
    /// # Errors
    ///
    /// As [`MzTabFile::protein_header`].
    pub fn oligonucleotide_header(&self, layout: &SectionLayout) -> Result<String> {
        let mut cells = vec![
            "OLH".to_owned(),
            "sequence".to_owned(),
            "accession".to_owned(),
            "unique".to_owned(),
            "search_engine".to_owned(),
        ];
        push_indexed(
            &mut cells,
            &layout.best_search_engine_score,
            "best_search_engine_score",
        );
        push_score_run_names(&mut cells, layout);
        if self.store_oligonucleotide_reliability {
            cells.push("reliability".to_owned());
        }
        cells.push("modifications".to_owned());
        cells.push("retention_time".to_owned());
        cells.push("retention_time_window".to_owned());
        if self.store_oligonucleotide_uri {
            cells.push("uri".to_owned());
        }
        cells.push("pre".to_owned());
        cells.push("post".to_owned());
        cells.push("start".to_owned());
        cells.push("end".to_owned());
        cells.extend(layout.optional_columns.iter().cloned());
        join(cells)
    }

    /// Render one `OLI` row against `layout`.
    ///
    /// Source `MzTabFile::generateMzTabSectionRow_(const MzTabOligonucleotideSectionRow&, …)`.
    ///
    /// # Errors
    ///
    /// As [`MzTabFile::protein_row`].
    pub fn oligonucleotide_row(
        &self,
        row: &MzTabOligonucleotideSectionRow,
        layout: &SectionLayout,
    ) -> Result<String> {
        let mut cells = vec![
            "OLI".to_owned(),
            row.sequence.to_cell_string(),
            row.accession.to_cell_string(),
            row.unique.to_cell_string(),
            row.search_engine.to_cell_string(),
        ];
        push_indexed_cells(
            &mut cells,
            &layout.best_search_engine_score,
            &row.best_search_engine_score,
        )?;
        push_score_run_cells(&mut cells, layout, &row.search_engine_score_ms_run)?;
        if self.store_oligonucleotide_reliability {
            cells.push(row.reliability.to_cell_string());
        }
        cells.push(row.modifications.to_cell_string()?);
        cells.push(row.retention_time.to_cell_string());
        cells.push(row.retention_time_window.to_cell_string());
        if self.store_oligonucleotide_uri {
            cells.push(row.uri.to_cell_string());
        }
        cells.push(row.pre.to_cell_string());
        cells.push(row.post.to_cell_string());
        cells.push(row.start.to_cell_string());
        cells.push(row.end.to_cell_string());
        cells.extend(optional_column_cells(&layout.optional_columns, &row.opt)?);
        join(cells)
    }

    /// Render the `OSH` header line for a section with this `layout`.
    ///
    /// Source `MzTabFile::generateMzTabOSMHeader_`.
    ///
    /// # Errors
    ///
    /// As [`MzTabFile::protein_header`].
    pub fn osm_header(&self, layout: &SectionLayout) -> Result<String> {
        let mut cells = vec![
            "OSH".to_owned(),
            "sequence".to_owned(),
            "search_engine".to_owned(),
        ];
        push_indexed(
            &mut cells,
            &layout.search_engine_score,
            "search_engine_score",
        );
        if self.store_osm_reliability {
            cells.push("reliability".to_owned());
        }
        cells.push("modifications".to_owned());
        cells.push("retention_time".to_owned());
        cells.push("charge".to_owned());
        cells.push("exp_mass_to_charge".to_owned());
        cells.push("calc_mass_to_charge".to_owned());
        if self.store_osm_uri {
            cells.push("uri".to_owned());
        }
        cells.push("spectra_ref".to_owned());
        cells.extend(layout.optional_columns.iter().cloned());
        join(cells)
    }

    /// Render one `OSM` row against `layout`.
    ///
    /// Source `MzTabFile::generateMzTabSectionRow_(const MzTabOSMSectionRow&, …)`.
    ///
    /// # Errors
    ///
    /// As [`MzTabFile::protein_row`].
    pub fn osm_row(&self, row: &MzTabOSMSectionRow, layout: &SectionLayout) -> Result<String> {
        let mut cells = vec![
            "OSM".to_owned(),
            row.sequence.to_cell_string(),
            row.search_engine.to_cell_string(),
        ];
        push_indexed_cells(
            &mut cells,
            &layout.search_engine_score,
            &row.search_engine_score,
        )?;
        if self.store_osm_reliability {
            cells.push(row.reliability.to_cell_string());
        }
        cells.push(row.modifications.to_cell_string()?);
        cells.push(row.retention_time.to_cell_string());
        cells.push(row.charge.to_cell_string());
        cells.push(row.exp_mass_to_charge.to_cell_string());
        cells.push(row.calc_mass_to_charge.to_cell_string());
        if self.store_osm_uri {
            cells.push(row.uri.to_cell_string());
        }
        cells.push(row.spectra_ref.to_cell_string());
        cells.extend(optional_column_cells(&layout.optional_columns, &row.opt)?);
        join(cells)
    }
}

// ---------------------------------------------------------------------------
// Writing — the document
// ---------------------------------------------------------------------------

fn emit_section(
    out: &mut LineBuffer,
    header: String,
    rows: Vec<String>,
    label: &str,
) -> Result<()> {
    let expected = column_count(&header);
    if expected > MzTabFile::MAX_COLUMNS {
        return Err(bad("MzTab section exceeds its column limit"));
    }
    for row in &rows {
        if column_count(row) != expected {
            return Err(bad(format!(
                "MzTab {label} header and content differ in columns"
            )));
        }
    }
    out.push(String::new())?;
    out.push(header)?;
    out.extend(rows)
}

fn restore_comments_and_blanks(generated: Vec<String>, document: &MzTab) -> Result<Vec<String>> {
    let blanks: BTreeSet<usize> = document.empty_rows.iter().copied().collect();
    let comments = &document.comment_rows;
    if blanks.is_empty() && comments.is_empty() {
        return Ok(generated);
    }
    // The furthest recorded position, so the walk can reach a comment or blank
    // recorded past the point the generated lines run out. Checking it here
    // charges MzTabFile::MAX_LINES before the walk rather than during it.
    let last_recorded = blanks
        .iter()
        .next_back()
        .copied()
        .into_iter()
        .chain(comments.keys().next_back().copied())
        .max();
    if last_recorded.is_some_and(|last| last > MzTabFile::MAX_LINES) {
        return Err(bad("MzTab output exceeds its line limit"));
    }
    let mut out = LineBuffer::default();
    let mut line = 0usize;
    let mut pending = generated.into_iter().peekable();
    loop {
        if blanks.contains(&line) {
            // A recorded blank is satisfied by a generated blank when the next
            // generated line is one, so the section separators are not doubled.
            if pending.peek().is_some_and(|next| next.trim().is_empty()) {
                pending.next();
            }
            out.push(String::new())?;
        } else if let Some(comment) = comments.get(&line) {
            out.push(comment.clone())?;
        } else if let Some(next) = pending.next() {
            out.push(next)?;
        } else if last_recorded.is_none_or(|last| line >= last) {
            // Nothing generated and nothing recorded further down: done. When
            // something *is* recorded further down, the walk continues to it —
            // the source stops here and drops that tail. Nothing fills the gap,
            // so the tail moves up by as many lines as the gap holds.
            break;
        }
        line = line
            .checked_add(1)
            .filter(|&next| next <= MzTabFile::MAX_LINES)
            .ok_or_else(|| bad("MzTab output exceeds its line limit"))?;
    }
    Ok(out.lines)
}

impl MzTabFile {
    /// Render `document` as the lines of an MzTab file, without terminators.
    ///
    /// The sections follow the source's order — `PRT`, `PEP`, `PSM`, `SML`,
    /// `NUC`, `OLI`, `OSM` — each preceded by an empty line and introduced by
    /// its header row; an empty section contributes nothing. The recorded
    /// comment and empty lines are then restored at their original positions.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the document exceeds
    /// [`MzTabFile::MAX_LINES`], [`MzTabFile::MAX_BYTES`] or
    /// [`MzTabFile::MAX_COLUMNS`], when a cell carries a tab or a line break, or
    /// when a header and one of its rows disagree on the column count — which
    /// the layout makes impossible and which is therefore checked rather than
    /// assumed, in place of the source's `Exception::Postcondition`.
    /// [`Error::MissingInformation`] comes from a modification cell that carries
    /// positions but no identifier.
    ///
    /// # Notes
    ///
    /// The source restores comments and empty lines by walking the generated
    /// lines and *inserting* a blank whenever the current position was blank in
    /// the original file (`MzTabFile.cpp:3328`) — even though the generated
    /// lines already carry a blank before every section, so each one is emitted
    /// twice. Its own round-trip test cannot see this, because
    /// `FuzzyStringComparator::readNextLine_` skips blank lines outright. Here a
    /// recorded blank position consumes a generated blank when there is one, so
    /// the line count and the recorded positions survive a round trip. The
    /// source also stops as soon as the generated lines run out, dropping any
    /// comment or blank recorded past that point; those are emitted here. Their
    /// recorded positions are only reproducible when the recorded tail is
    /// contiguous with the end of the generated lines: nothing fills a gap left
    /// by a metadata key that was recorded but is not written back (a
    /// null-valued optional key, say), so a tail behind such a gap keeps its
    /// order and its content but moves up by the width of the gap.
    pub fn document_lines(&self, document: &MzTab) -> Result<Vec<String>> {
        let meta = &document.meta_data;
        let mut out = LineBuffer::default();
        out.extend(self.metadata_section_lines(meta)?)?;

        if !document.protein_data.is_empty() {
            let layout = SectionLayout::for_protein(&document.protein_data, meta)?;
            let header = self.protein_header(&layout)?;
            let mut rows = Vec::new();
            for row in &document.protein_data {
                rows.push(self.protein_row(row, &layout)?);
            }
            emit_section(&mut out, header, rows, "protein")?;
        }
        if !document.peptide_data.is_empty() {
            let layout = SectionLayout::for_peptide(&document.peptide_data, meta)?;
            let header = self.peptide_header(&layout)?;
            let mut rows = Vec::new();
            for row in &document.peptide_data {
                rows.push(self.peptide_row(row, &layout)?);
            }
            emit_section(&mut out, header, rows, "peptide")?;
        }
        if !document.psm_data.is_empty() {
            let layout = SectionLayout::for_psm(&document.psm_data, meta)?;
            let header = self.psm_header(&layout)?;
            let mut rows = Vec::new();
            for row in &document.psm_data {
                rows.push(self.psm_row(row, &layout)?);
            }
            emit_section(&mut out, header, rows, "PSM")?;
        }
        if !document.small_molecule_data.is_empty() {
            let layout = SectionLayout::for_small_molecule(&document.small_molecule_data, meta)?;
            let header = self.small_molecule_header(&layout)?;
            let mut rows = Vec::new();
            for row in &document.small_molecule_data {
                rows.push(self.small_molecule_row(row, &layout)?);
            }
            emit_section(&mut out, header, rows, "small molecule")?;
        }
        if !document.nucleic_acid_data.is_empty() {
            let layout = SectionLayout::for_nucleic_acid(&document.nucleic_acid_data, meta)?;
            let header = self.nucleic_acid_header(&layout)?;
            let mut rows = Vec::new();
            for row in &document.nucleic_acid_data {
                rows.push(self.nucleic_acid_row(row, &layout)?);
            }
            emit_section(&mut out, header, rows, "nucleic acid")?;
        }
        if !document.oligonucleotide_data.is_empty() {
            let layout = SectionLayout::for_oligonucleotide(&document.oligonucleotide_data, meta)?;
            let header = self.oligonucleotide_header(&layout)?;
            let mut rows = Vec::new();
            for row in &document.oligonucleotide_data {
                rows.push(self.oligonucleotide_row(row, &layout)?);
            }
            emit_section(&mut out, header, rows, "oligonucleotide")?;
        }
        if !document.osm_data.is_empty() {
            let layout = SectionLayout::for_osm(&document.osm_data, meta)?;
            let header = self.osm_header(&layout)?;
            let mut rows = Vec::new();
            for row in &document.osm_data {
                rows.push(self.osm_row(row, &layout)?);
            }
            emit_section(&mut out, header, rows, "OSM")?;
        }

        restore_comments_and_blanks(out.lines, document)
    }

    /// Render `document` as the complete text of an MzTab file, every line
    /// terminated with `\n`.
    ///
    /// Native convenience; the source can only write to a file.
    ///
    /// # Errors
    ///
    /// As [`MzTabFile::document_lines`].
    pub fn write_to_string(&self, document: &MzTab) -> Result<String> {
        let lines = self.document_lines(document)?;
        let mut text = String::new();
        for line in lines {
            text.push_str(&line);
            text.push('\n');
        }
        Ok(text)
    }

    /// Write `document` to `path`.
    ///
    /// Source `MzTabFile::store(const std::string&, const MzTab&) const`. The
    /// extension must name an mzTab or a tab-separated file, or be
    /// unrecognised; the source throws `Exception::UnableToCreateFile`
    /// otherwise.
    ///
    /// The whole document is rendered and checked before anything is created,
    /// and the bytes go to a sibling temporary file that is renamed into place,
    /// so a failure never leaves a truncated output where a previous version
    /// stood. The source opens the destination with `ios::trunc` first and
    /// writes as it goes.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a rejected extension and for every ceiling of
    /// [`MzTabFile::document_lines`], [`Error::MissingInformation`] from an
    /// unrenderable modification cell, and [`Error::Io`] for a filesystem
    /// failure.
    pub fn store(&self, path: impl AsRef<Path>, document: &MzTab) -> Result<()> {
        let path = path.as_ref();
        let name = path.to_string_lossy();
        if !(file_types::has_valid_extension(&name, FileType::MzTab)
            || file_types::has_valid_extension(&name, FileType::Tsv))
        {
            return Err(bad(format!(
                "invalid file extension for MzTab output, expected '{}' or '{}': {name}",
                FileType::MzTab.name(),
                FileType::Tsv.name()
            )));
        }
        let text = self.write_to_string(document)?;
        path_io::write(path, |writer| {
            writer.write_all(text.as_bytes())?;
            Ok(())
        })
    }
}

// ---------------------------------------------------------------------------
// Reading — column bookkeeping
// ---------------------------------------------------------------------------

/// Which column of a section's rows holds which member of which family.
///
/// `None` means the header did not declare the column. The source uses `0` as
/// that sentinel, which is also the index of the section tag, so for every
/// column the header omits it reads `cells[0]` — the tag itself — into the
/// field, and for a field whose type is numeric that turns into a conversion
/// error on a document the specification allows.
#[derive(Clone, Debug, Default)]
struct Columns {
    single: BTreeMap<&'static str, usize>,
    best_search_engine_score: BTreeMap<usize, usize>,
    search_engine_score: BTreeMap<usize, usize>,
    search_engine_score_ms_run: BTreeMap<usize, (usize, usize)>,
    count_matches_ms_run: BTreeMap<usize, usize>,
    count_distinct_ms_run: BTreeMap<usize, usize>,
    count_unique_ms_run: BTreeMap<usize, usize>,
    abundance_assay: BTreeMap<usize, usize>,
    abundance_study_variable: BTreeMap<usize, usize>,
    abundance_stdev_study_variable: BTreeMap<usize, usize>,
    abundance_std_error_study_variable: BTreeMap<usize, usize>,
    optional: BTreeMap<String, usize>,
    reported: bool,
}

impl Columns {
    fn get(&self, name: &str) -> Option<usize> {
        self.single.get(name).copied()
    }

    /// Highest column index the header declared, for the short-row diagnostic.
    fn last_column(&self) -> Option<usize> {
        self.single
            .values()
            .chain(self.best_search_engine_score.values())
            .chain(self.search_engine_score.values())
            .chain(self.search_engine_score_ms_run.keys())
            .chain(self.count_matches_ms_run.values())
            .chain(self.count_distinct_ms_run.values())
            .chain(self.count_unique_ms_run.values())
            .chain(self.abundance_assay.values())
            .chain(self.abundance_study_variable.values())
            .chain(self.abundance_stdev_study_variable.values())
            .chain(self.abundance_std_error_study_variable.values())
            .chain(self.optional.values())
            .copied()
            .max()
    }
}

/// One indexed column family and where its columns go.
enum Family {
    BestScore,
    Score,
    ScoreMsRun,
    CountMatches,
    CountDistinct,
    CountUnique,
    AbundanceAssay,
    AbundanceStudyVariable,
    AbundanceStdevStudyVariable,
    AbundanceStdErrorStudyVariable,
}

fn record_family(
    columns: &mut Columns,
    families: &[(&str, Family)],
    name: &str,
    index: usize,
    line: usize,
) -> Result<bool> {
    for (prefix, family) in families {
        if !name.starts_with(prefix) {
            continue;
        }
        let target = match family {
            Family::ScoreMsRun => {
                let pair = extract_index_pairs_from_brackets(name).map_err(|_| {
                    parse_error(
                        line + 1,
                        format!("MzTab column {name:?} does not carry two valid bracketed indices"),
                    )
                })?;
                columns.search_engine_score_ms_run.insert(index, pair);
                return Ok(true);
            }
            Family::BestScore => &mut columns.best_search_engine_score,
            Family::Score => &mut columns.search_engine_score,
            Family::CountMatches => &mut columns.count_matches_ms_run,
            Family::CountDistinct => &mut columns.count_distinct_ms_run,
            Family::CountUnique => &mut columns.count_unique_ms_run,
            Family::AbundanceAssay => &mut columns.abundance_assay,
            Family::AbundanceStudyVariable => &mut columns.abundance_study_variable,
            Family::AbundanceStdevStudyVariable => &mut columns.abundance_stdev_study_variable,
            Family::AbundanceStdErrorStudyVariable => {
                &mut columns.abundance_std_error_study_variable
            }
        };
        let member = extract_bracket_index(name, prefix).map_err(|_| {
            parse_error(
                line + 1,
                format!("MzTab column {name:?} has no valid bracketed index"),
            )
        })?;
        target.insert(member, index);
        return Ok(true);
    }
    Ok(false)
}

/// Parse one section header row.
///
/// A header line replaces whatever the previous one of the same section
/// declared; the source accumulates into the same maps, so a second header
/// leaves the first one's columns in place.
fn parse_header(
    singles: &[&'static str],
    families: &[(&str, Family)],
    cells: &[&str],
    line: usize,
) -> Result<Columns> {
    let mut columns = Columns::default();
    if cells.len() > MzTabFile::MAX_COLUMNS {
        return Err(bad("MzTab header exceeds its column limit"));
    }
    for (index, name) in cells.iter().enumerate().skip(1) {
        if let Some(found) = singles.iter().copied().find(|known| known == name) {
            columns.single.insert(found, index);
            continue;
        }
        if record_family(&mut columns, families, name, index, line)? {
            continue;
        }
        if name.starts_with("opt_") {
            if columns.optional.len() >= MzTab::MAX_OPTIONAL_COLUMNS {
                return Err(bad("MzTab section exceeds its optional-column limit"));
            }
            columns.optional.insert((*name).to_owned(), index);
        }
    }
    Ok(columns)
}

fn read_into<T: MzTabCell>(target: &mut T, cells: &[&str], index: Option<usize>) -> Result<()> {
    if let Some(text) = index.and_then(|column| cells.get(column)) {
        target.read_cell(text)?;
    }
    Ok(())
}

fn read_named<T: MzTabCell>(
    target: &mut T,
    cells: &[&str],
    columns: &Columns,
    name: &str,
) -> Result<()> {
    read_into(target, cells, columns.get(name))
}

fn read_family<T: Default + MzTabCell>(
    target: &mut BTreeMap<usize, T>,
    cells: &[&str],
    family: &BTreeMap<usize, usize>,
) -> Result<()> {
    for (member, column) in family {
        if let Some(text) = cells.get(*column) {
            let mut value = T::default();
            value.read_cell(text)?;
            target.insert(*member, value);
        }
    }
    Ok(())
}

fn read_score_ms_run(target: &mut ScoreRunMap, cells: &[&str], columns: &Columns) -> Result<()> {
    for (column, (score, run)) in &columns.search_engine_score_ms_run {
        if let Some(text) = cells.get(*column) {
            let mut value = MzTabDouble::default();
            value.read_cell(text)?;
            target.entry(*score).or_default().insert(*run, value);
        }
    }
    Ok(())
}

fn read_optional(
    target: &mut Vec<MzTabOptionalColumnEntry>,
    cells: &[&str],
    columns: &Columns,
) -> Result<()> {
    let mut staged = Vec::with_capacity(columns.optional.len());
    for (name, column) in &columns.optional {
        if let Some(text) = cells.get(*column) {
            staged.push(MzTabOptionalColumnEntry::new(
                name.clone(),
                MzTabString::from_text(text),
            ));
        }
    }
    target
        .try_reserve(staged.len())
        .map_err(|_| bad("MzTab optional column allocation failed"))?;
    target.append(&mut staged);
    Ok(())
}

fn note(diagnostics: &mut Vec<String>, message: String) {
    if diagnostics.len() < MAX_DIAGNOSTICS {
        diagnostics.push(message);
    }
}

fn note_short_row(
    diagnostics: &mut Vec<String>,
    columns: &Columns,
    cells: &[&str],
    label: &str,
    line: usize,
) {
    if columns
        .last_column()
        .is_some_and(|last| last >= cells.len())
    {
        note(
            diagnostics,
            format!(
                "MzTab {label} row on line {} is shorter than its header; the missing cells stay unset",
                line + 1
            ),
        );
    }
}

fn push_row<R>(rows: &mut Vec<R>, row: R) -> Result<()> {
    if rows.len() >= MzTab::MAX_ROWS {
        return Err(bad("MzTab section exceeds its row limit"));
    }
    rows.try_reserve(1)
        .map_err(|_| bad("MzTab row allocation failed"))?;
    rows.push(row);
    Ok(())
}

const PROTEIN_SINGLES: &[&str] = &[
    "accession",
    "description",
    "taxid",
    "species",
    "database",
    "database_version",
    "search_engine",
    "reliability",
    "ambiguity_members",
    "modifications",
    "uri",
    "go_terms",
    "protein_coverage",
];

const PROTEIN_FAMILIES: &[(&str, Family)] = &[
    ("best_search_engine_score[", Family::BestScore),
    ("search_engine_score[", Family::ScoreMsRun),
    ("num_psms_ms_run[", Family::CountMatches),
    ("num_peptides_distinct_ms_run[", Family::CountDistinct),
    ("num_peptides_unique_ms_run[", Family::CountUnique),
    (
        "protein_abundance_stdev_study_variable[",
        Family::AbundanceStdevStudyVariable,
    ),
    (
        "protein_abundance_std_error_study_variable[",
        Family::AbundanceStdErrorStudyVariable,
    ),
    (
        "protein_abundance_study_variable[",
        Family::AbundanceStudyVariable,
    ),
    ("protein_abundance_assay[", Family::AbundanceAssay),
];

const PEPTIDE_SINGLES: &[&str] = &[
    "sequence",
    "accession",
    "unique",
    "database",
    "database_version",
    "search_engine",
    "reliability",
    "modifications",
    "retention_time",
    "retention_time_window",
    "charge",
    "mass_to_charge",
    "uri",
    "spectra_ref",
];

const PEPTIDE_FAMILIES: &[(&str, Family)] = &[
    ("best_search_engine_score[", Family::BestScore),
    ("search_engine_score[", Family::ScoreMsRun),
    (
        "peptide_abundance_stdev_study_variable[",
        Family::AbundanceStdevStudyVariable,
    ),
    (
        "peptide_abundance_std_error_study_variable[",
        Family::AbundanceStdErrorStudyVariable,
    ),
    (
        "peptide_abundance_study_variable[",
        Family::AbundanceStudyVariable,
    ),
    ("peptide_abundance_assay[", Family::AbundanceAssay),
];

const PSM_SINGLES: &[&str] = &[
    "sequence",
    "PSM_ID",
    "accession",
    "unique",
    "database",
    "database_version",
    "search_engine",
    "reliability",
    "modifications",
    "retention_time",
    "charge",
    "exp_mass_to_charge",
    "calc_mass_to_charge",
    "uri",
    "spectra_ref",
    "pre",
    "post",
    "start",
    "end",
];

const PSM_FAMILIES: &[(&str, Family)] = &[("search_engine_score[", Family::Score)];

const SMALL_MOLECULE_SINGLES: &[&str] = &[
    "identifier",
    "chemical_formula",
    "smiles",
    "inchi_key",
    "description",
    "exp_mass_to_charge",
    "calc_mass_to_charge",
    "charge",
    "retention_time",
    "taxid",
    "species",
    "database",
    "database_version",
    "reliability",
    "uri",
    "spectra_ref",
    "search_engine",
    "modifications",
];

const SMALL_MOLECULE_FAMILIES: &[(&str, Family)] = &[
    ("best_search_engine_score[", Family::BestScore),
    ("search_engine_score[", Family::ScoreMsRun),
    (
        "smallmolecule_abundance_stdev_study_variable[",
        Family::AbundanceStdevStudyVariable,
    ),
    (
        "smallmolecule_abundance_std_error_study_variable[",
        Family::AbundanceStdErrorStudyVariable,
    ),
    (
        "smallmolecule_abundance_study_variable[",
        Family::AbundanceStudyVariable,
    ),
    ("smallmolecule_abundance_assay[", Family::AbundanceAssay),
];

const NUCLEIC_ACID_SINGLES: &[&str] = &[
    "accession",
    "description",
    "taxid",
    "species",
    "database",
    "database_version",
    "search_engine",
    "reliability",
    "ambiguity_members",
    "modifications",
    "uri",
    "go_terms",
    "sequence_coverage",
];

const NUCLEIC_ACID_FAMILIES: &[(&str, Family)] = &[
    ("best_search_engine_score[", Family::BestScore),
    ("search_engine_score[", Family::ScoreMsRun),
    ("num_osms_ms_run[", Family::CountMatches),
    ("num_oligos_distinct_ms_run[", Family::CountDistinct),
    ("num_oligos_unique_ms_run[", Family::CountUnique),
];

const OLIGONUCLEOTIDE_SINGLES: &[&str] = &[
    "sequence",
    "accession",
    "unique",
    "search_engine",
    "reliability",
    "modifications",
    "retention_time",
    "retention_time_window",
    "uri",
    "pre",
    "post",
    "start",
    "end",
];

const OLIGONUCLEOTIDE_FAMILIES: &[(&str, Family)] = &[
    ("best_search_engine_score[", Family::BestScore),
    ("search_engine_score[", Family::ScoreMsRun),
];

const OSM_SINGLES: &[&str] = &[
    "sequence",
    "search_engine",
    "reliability",
    "modifications",
    "retention_time",
    "charge",
    "exp_mass_to_charge",
    "calc_mass_to_charge",
    "uri",
    "spectra_ref",
];

const OSM_FAMILIES: &[(&str, Family)] = &[("search_engine_score[", Family::Score)];

// ---------------------------------------------------------------------------
// Reading — mandatory-column diagnostics
// ---------------------------------------------------------------------------

/// Every column the source reports as mandatory when the protein header omits
/// it, in the order it checks them.
///
/// The source writes each of these to `std::cout` once per data row, so a
/// section of a thousand rows prints the same line a thousand times, and the
/// message never reaches the caller. [`MzTabFile::load_reporting`] returns them
/// instead, once per section.
fn mandatory_protein_reports(columns: &Columns, meta: &MzTabMetaData) -> Vec<String> {
    let mut reports = Vec::new();
    for name in [
        "accession",
        "description",
        "taxid",
        "species",
        "database",
        "database_version",
        "search_engine",
    ] {
        if columns.get(name).is_none() {
            reports.push(format!("mandatory protein {name} column missing"));
        }
    }
    if columns.best_search_engine_score.is_empty() {
        reports.push("mandatory protein best_search_engine_score[1-n] column missing".to_owned());
    }
    let complete = meta.mz_tab_mode.get() == "Complete";
    let identification = meta.mz_tab_type.get() == "Identification";
    let quantification = meta.mz_tab_type.get() == "Quantification";
    if complete && columns.search_engine_score_ms_run.len() != meta.ms_run.len() {
        reports.push(format!(
            "mandatory protein search_engine_score_ms_run column(s) missing: expected {} ms_runs from the metadata section but {} provide score columns",
            meta.ms_run.len(),
            columns.search_engine_score_ms_run.len()
        ));
    }
    if complete && identification {
        for (set, name) in [
            (&columns.count_matches_ms_run, "num_psms_ms_run"),
            (
                &columns.count_distinct_ms_run,
                "num_peptides_distinct_ms_run",
            ),
            (&columns.count_unique_ms_run, "num_peptides_unique_ms_run"),
        ] {
            if set.is_empty() {
                reports.push(format!("mandatory protein {name} column(s) missing"));
            }
        }
    }
    for name in ["ambiguity_members", "modifications", "protein_coverage"] {
        if columns.get(name).is_none() {
            reports.push(format!("mandatory protein {name} column missing"));
        }
    }
    if complete && quantification && columns.abundance_assay.is_empty() {
        reports.push("mandatory protein protein_abundance_assay column(s) missing".to_owned());
    }
    if quantification {
        for (set, name) in [
            (
                &columns.abundance_study_variable,
                "protein_abundance_study_variable",
            ),
            (
                &columns.abundance_stdev_study_variable,
                "protein_abundance_stdev_study_variable",
            ),
            (
                &columns.abundance_std_error_study_variable,
                "protein_abundance_std_error_study_variable",
            ),
        ] {
            if set.is_empty() {
                reports.push(format!("mandatory {name} column(s) missing"));
            }
        }
    }
    reports
}

// ---------------------------------------------------------------------------
// Reading — section rows
// ---------------------------------------------------------------------------

fn read_protein_row(
    columns: &mut Columns,
    document: &mut MzTab,
    diagnostics: &mut Vec<String>,
    cells: &[&str],
    line: usize,
) -> Result<()> {
    if !columns.reported {
        columns.reported = true;
        for report in mandatory_protein_reports(columns, &document.meta_data) {
            note(diagnostics, report);
        }
    }
    note_short_row(diagnostics, columns, cells, "protein", line);
    let mut row = MzTabProteinSectionRow::default();
    read_named(&mut row.accession, cells, columns, "accession")?;
    read_named(&mut row.description, cells, columns, "description")?;
    read_named(&mut row.taxid, cells, columns, "taxid")?;
    read_named(&mut row.species, cells, columns, "species")?;
    read_named(&mut row.database, cells, columns, "database")?;
    read_named(
        &mut row.database_version,
        cells,
        columns,
        "database_version",
    )?;
    read_named(&mut row.search_engine, cells, columns, "search_engine")?;
    read_family(
        &mut row.best_search_engine_score,
        cells,
        &columns.best_search_engine_score,
    )?;
    read_score_ms_run(&mut row.search_engine_score_ms_run, cells, columns)?;
    read_named(&mut row.reliability, cells, columns, "reliability")?;
    read_family(
        &mut row.num_psms_ms_run,
        cells,
        &columns.count_matches_ms_run,
    )?;
    read_family(
        &mut row.num_peptides_distinct_ms_run,
        cells,
        &columns.count_distinct_ms_run,
    )?;
    read_family(
        &mut row.num_peptides_unique_ms_run,
        cells,
        &columns.count_unique_ms_run,
    )?;
    read_named(
        &mut row.ambiguity_members,
        cells,
        columns,
        "ambiguity_members",
    )?;
    read_named(&mut row.modifications, cells, columns, "modifications")?;
    read_named(&mut row.uri, cells, columns, "uri")?;
    read_named(&mut row.go_terms, cells, columns, "go_terms")?;
    read_named(&mut row.coverage, cells, columns, "protein_coverage")?;
    read_family(
        &mut row.protein_abundance_assay,
        cells,
        &columns.abundance_assay,
    )?;
    read_family(
        &mut row.protein_abundance_study_variable,
        cells,
        &columns.abundance_study_variable,
    )?;
    read_family(
        &mut row.protein_abundance_stdev_study_variable,
        cells,
        &columns.abundance_stdev_study_variable,
    )?;
    read_family(
        &mut row.protein_abundance_std_error_study_variable,
        cells,
        &columns.abundance_std_error_study_variable,
    )?;
    read_optional(&mut row.opt, cells, columns)?;
    push_row(&mut document.protein_data, row)
}

fn read_peptide_row(
    columns: &Columns,
    document: &mut MzTab,
    diagnostics: &mut Vec<String>,
    cells: &[&str],
    line: usize,
) -> Result<()> {
    note_short_row(diagnostics, columns, cells, "peptide", line);
    let mut row = MzTabPeptideSectionRow::default();
    read_named(&mut row.sequence, cells, columns, "sequence")?;
    read_named(&mut row.accession, cells, columns, "accession")?;
    read_named(&mut row.unique, cells, columns, "unique")?;
    read_named(&mut row.database, cells, columns, "database")?;
    read_named(
        &mut row.database_version,
        cells,
        columns,
        "database_version",
    )?;
    read_named(&mut row.search_engine, cells, columns, "search_engine")?;
    read_family(
        &mut row.best_search_engine_score,
        cells,
        &columns.best_search_engine_score,
    )?;
    read_score_ms_run(&mut row.search_engine_score_ms_run, cells, columns)?;
    read_named(&mut row.reliability, cells, columns, "reliability")?;
    read_named(&mut row.modifications, cells, columns, "modifications")?;
    read_named(&mut row.retention_time, cells, columns, "retention_time")?;
    read_named(
        &mut row.retention_time_window,
        cells,
        columns,
        "retention_time_window",
    )?;
    read_named(&mut row.charge, cells, columns, "charge")?;
    read_named(&mut row.mass_to_charge, cells, columns, "mass_to_charge")?;
    read_named(&mut row.uri, cells, columns, "uri")?;
    read_named(&mut row.spectra_ref, cells, columns, "spectra_ref")?;
    read_family(
        &mut row.peptide_abundance_assay,
        cells,
        &columns.abundance_assay,
    )?;
    read_family(
        &mut row.peptide_abundance_study_variable,
        cells,
        &columns.abundance_study_variable,
    )?;
    read_family(
        &mut row.peptide_abundance_stdev_study_variable,
        cells,
        &columns.abundance_stdev_study_variable,
    )?;
    read_family(
        &mut row.peptide_abundance_std_error_study_variable,
        cells,
        &columns.abundance_std_error_study_variable,
    )?;
    read_optional(&mut row.opt, cells, columns)?;
    push_row(&mut document.peptide_data, row)
}

fn read_psm_row(
    columns: &Columns,
    document: &mut MzTab,
    diagnostics: &mut Vec<String>,
    cells: &[&str],
    line: usize,
) -> Result<()> {
    note_short_row(diagnostics, columns, cells, "PSM", line);
    let mut row = MzTabPSMSectionRow::default();
    read_named(&mut row.sequence, cells, columns, "sequence")?;
    read_named(&mut row.psm_id, cells, columns, "PSM_ID")?;
    read_named(&mut row.accession, cells, columns, "accession")?;
    read_named(&mut row.unique, cells, columns, "unique")?;
    read_named(&mut row.database, cells, columns, "database")?;
    read_named(
        &mut row.database_version,
        cells,
        columns,
        "database_version",
    )?;
    read_named(&mut row.search_engine, cells, columns, "search_engine")?;
    read_family(
        &mut row.search_engine_score,
        cells,
        &columns.search_engine_score,
    )?;
    read_named(&mut row.reliability, cells, columns, "reliability")?;
    read_named(&mut row.modifications, cells, columns, "modifications")?;
    read_named(&mut row.retention_time, cells, columns, "retention_time")?;
    read_named(&mut row.charge, cells, columns, "charge")?;
    read_named(
        &mut row.exp_mass_to_charge,
        cells,
        columns,
        "exp_mass_to_charge",
    )?;
    read_named(
        &mut row.calc_mass_to_charge,
        cells,
        columns,
        "calc_mass_to_charge",
    )?;
    read_named(&mut row.uri, cells, columns, "uri")?;
    read_named(&mut row.spectra_ref, cells, columns, "spectra_ref")?;
    read_named(&mut row.pre, cells, columns, "pre")?;
    read_named(&mut row.post, cells, columns, "post")?;
    read_named(&mut row.start, cells, columns, "start")?;
    read_named(&mut row.end, cells, columns, "end")?;
    read_optional(&mut row.opt, cells, columns)?;
    push_row(&mut document.psm_data, row)
}

fn read_small_molecule_row(
    columns: &Columns,
    document: &mut MzTab,
    diagnostics: &mut Vec<String>,
    cells: &[&str],
    line: usize,
) -> Result<()> {
    note_short_row(diagnostics, columns, cells, "small molecule", line);
    let mut row = MzTabSmallMoleculeSectionRow::default();
    read_named(&mut row.identifier, cells, columns, "identifier")?;
    read_named(
        &mut row.chemical_formula,
        cells,
        columns,
        "chemical_formula",
    )?;
    read_named(&mut row.smiles, cells, columns, "smiles")?;
    read_named(&mut row.inchi_key, cells, columns, "inchi_key")?;
    read_named(&mut row.description, cells, columns, "description")?;
    read_named(
        &mut row.exp_mass_to_charge,
        cells,
        columns,
        "exp_mass_to_charge",
    )?;
    read_named(
        &mut row.calc_mass_to_charge,
        cells,
        columns,
        "calc_mass_to_charge",
    )?;
    read_named(&mut row.charge, cells, columns, "charge")?;
    read_named(&mut row.retention_time, cells, columns, "retention_time")?;
    read_named(&mut row.taxid, cells, columns, "taxid")?;
    read_named(&mut row.species, cells, columns, "species")?;
    read_named(&mut row.database, cells, columns, "database")?;
    read_named(
        &mut row.database_version,
        cells,
        columns,
        "database_version",
    )?;
    read_named(&mut row.reliability, cells, columns, "reliability")?;
    read_named(&mut row.uri, cells, columns, "uri")?;
    read_named(&mut row.spectra_ref, cells, columns, "spectra_ref")?;
    read_named(&mut row.search_engine, cells, columns, "search_engine")?;
    read_family(
        &mut row.best_search_engine_score,
        cells,
        &columns.best_search_engine_score,
    )?;
    read_score_ms_run(&mut row.search_engine_score_ms_run, cells, columns)?;
    read_named(&mut row.modifications, cells, columns, "modifications")?;
    read_family(
        &mut row.smallmolecule_abundance_assay,
        cells,
        &columns.abundance_assay,
    )?;
    read_family(
        &mut row.smallmolecule_abundance_study_variable,
        cells,
        &columns.abundance_study_variable,
    )?;
    read_family(
        &mut row.smallmolecule_abundance_stdev_study_variable,
        cells,
        &columns.abundance_stdev_study_variable,
    )?;
    read_family(
        &mut row.smallmolecule_abundance_std_error_study_variable,
        cells,
        &columns.abundance_std_error_study_variable,
    )?;
    read_optional(&mut row.opt, cells, columns)?;
    push_row(&mut document.small_molecule_data, row)
}

fn read_nucleic_acid_row(
    columns: &Columns,
    document: &mut MzTab,
    diagnostics: &mut Vec<String>,
    cells: &[&str],
    line: usize,
) -> Result<()> {
    note_short_row(diagnostics, columns, cells, "nucleic acid", line);
    let mut row = MzTabNucleicAcidSectionRow::default();
    read_named(&mut row.accession, cells, columns, "accession")?;
    read_named(&mut row.description, cells, columns, "description")?;
    read_named(&mut row.taxid, cells, columns, "taxid")?;
    read_named(&mut row.species, cells, columns, "species")?;
    read_named(&mut row.database, cells, columns, "database")?;
    read_named(
        &mut row.database_version,
        cells,
        columns,
        "database_version",
    )?;
    read_named(&mut row.search_engine, cells, columns, "search_engine")?;
    read_family(
        &mut row.best_search_engine_score,
        cells,
        &columns.best_search_engine_score,
    )?;
    read_score_ms_run(&mut row.search_engine_score_ms_run, cells, columns)?;
    read_named(&mut row.reliability, cells, columns, "reliability")?;
    read_family(
        &mut row.num_osms_ms_run,
        cells,
        &columns.count_matches_ms_run,
    )?;
    read_family(
        &mut row.num_oligos_distinct_ms_run,
        cells,
        &columns.count_distinct_ms_run,
    )?;
    read_family(
        &mut row.num_oligos_unique_ms_run,
        cells,
        &columns.count_unique_ms_run,
    )?;
    read_named(
        &mut row.ambiguity_members,
        cells,
        columns,
        "ambiguity_members",
    )?;
    read_named(&mut row.modifications, cells, columns, "modifications")?;
    read_named(&mut row.uri, cells, columns, "uri")?;
    read_named(&mut row.go_terms, cells, columns, "go_terms")?;
    read_named(&mut row.coverage, cells, columns, "sequence_coverage")?;
    read_optional(&mut row.opt, cells, columns)?;
    push_row(&mut document.nucleic_acid_data, row)
}

fn read_oligonucleotide_row(
    columns: &Columns,
    document: &mut MzTab,
    diagnostics: &mut Vec<String>,
    cells: &[&str],
    line: usize,
) -> Result<()> {
    note_short_row(diagnostics, columns, cells, "oligonucleotide", line);
    let mut row = MzTabOligonucleotideSectionRow::default();
    read_named(&mut row.sequence, cells, columns, "sequence")?;
    read_named(&mut row.accession, cells, columns, "accession")?;
    read_named(&mut row.unique, cells, columns, "unique")?;
    read_named(&mut row.search_engine, cells, columns, "search_engine")?;
    read_family(
        &mut row.best_search_engine_score,
        cells,
        &columns.best_search_engine_score,
    )?;
    read_score_ms_run(&mut row.search_engine_score_ms_run, cells, columns)?;
    read_named(&mut row.reliability, cells, columns, "reliability")?;
    read_named(&mut row.modifications, cells, columns, "modifications")?;
    read_named(&mut row.retention_time, cells, columns, "retention_time")?;
    read_named(
        &mut row.retention_time_window,
        cells,
        columns,
        "retention_time_window",
    )?;
    read_named(&mut row.uri, cells, columns, "uri")?;
    read_named(&mut row.pre, cells, columns, "pre")?;
    read_named(&mut row.post, cells, columns, "post")?;
    read_named(&mut row.start, cells, columns, "start")?;
    read_named(&mut row.end, cells, columns, "end")?;
    read_optional(&mut row.opt, cells, columns)?;
    push_row(&mut document.oligonucleotide_data, row)
}

fn read_osm_row(
    columns: &Columns,
    document: &mut MzTab,
    diagnostics: &mut Vec<String>,
    cells: &[&str],
    line: usize,
) -> Result<()> {
    note_short_row(diagnostics, columns, cells, "OSM", line);
    let mut row = MzTabOSMSectionRow::default();
    read_named(&mut row.sequence, cells, columns, "sequence")?;
    read_named(&mut row.search_engine, cells, columns, "search_engine")?;
    read_family(
        &mut row.search_engine_score,
        cells,
        &columns.search_engine_score,
    )?;
    read_named(&mut row.reliability, cells, columns, "reliability")?;
    read_named(&mut row.modifications, cells, columns, "modifications")?;
    read_named(&mut row.retention_time, cells, columns, "retention_time")?;
    read_named(&mut row.charge, cells, columns, "charge")?;
    read_named(
        &mut row.exp_mass_to_charge,
        cells,
        columns,
        "exp_mass_to_charge",
    )?;
    read_named(
        &mut row.calc_mass_to_charge,
        cells,
        columns,
        "calc_mass_to_charge",
    )?;
    read_named(&mut row.uri, cells, columns, "uri")?;
    read_named(&mut row.spectra_ref, cells, columns, "spectra_ref")?;
    read_optional(&mut row.opt, cells, columns)?;
    push_row(&mut document.osm_data, row)
}

// ---------------------------------------------------------------------------
// Reading — metadata
// ---------------------------------------------------------------------------

fn parse_reference_list(text: &str, label: &str, line: usize) -> Result<Vec<i32>> {
    // Charge the entry ceiling by counting separators, before anything is
    // collected: stripping the brackets never removes a comma, so this is the
    // field count the split below would produce. Collecting first would build
    // one `&str` per field — several million for a 16 MiB line of commas —
    // before the ceiling was consulted.
    if text.matches(',').count().saturating_add(1) > MzTabFile::MAX_COLUMNS {
        return Err(bad("MzTab reference list exceeds its entry limit"));
    }
    let stripped = text.replace(&format!("{label}["), "").replace(']', "");
    let fields: Vec<&str> = if stripped.is_empty() {
        Vec::new()
    } else {
        stripped.split(',').collect()
    };
    if fields.len() > MzTabFile::MAX_COLUMNS {
        return Err(bad("MzTab reference list exceeds its entry limit"));
    }
    let mut values = Vec::with_capacity(fields.len());
    for field in fields {
        let token = trim_source(field);
        let value: i32 = token
            .strip_prefix('+')
            .unwrap_or(token)
            .parse()
            .map_err(|_| {
                parse_error(
                    line + 1,
                    format!(
                        "could not convert MzTab {label} reference {token:?} to an integer value"
                    ),
                )
            })?;
        values.push(value);
    }
    Ok(values)
}

/// Parse one `MTD` line into `meta`.
///
/// Mirrors the source's if/else chain over the `-`-separated fields of the key,
/// branch for branch and in the same order, with two structural differences.
/// The source indexes `meta_key_fields[1]` and `meta_key_fields[2]` in more than
/// thirty branches without checking the field count, reading past the end of the
/// vector for a key such as `instrument[1]` that carries no suffix; here a
/// missing field simply fails to match. And the source reads field zero of an
/// empty vector when the key cell is empty, which this cannot do.
///
/// An unrecognised key is ignored, as in the source.
fn read_metadata_line(meta: &mut MzTabMetaData, cells: &[&str], line: usize) -> Result<()> {
    let Some((&key, &value)) = cells.get(1).zip(cells.get(2)) else {
        return Ok(());
    };
    let fields: Vec<&str> = key.split('-').collect();
    let head = fields.first().copied().unwrap_or("");
    let second = fields.get(1).copied();
    let third = fields.get(2).copied();
    let at = |message: String| parse_error(line + 1, message);
    let index_of = |text: &str, prefix: &str| {
        extract_bracket_index(text, prefix)
            .map_err(|_| at(format!("MzTab key {key:?} does not carry a valid index")))
    };

    // The file-level keys are matched against the whole key cell, as the source
    // does, because splitting on `-` cuts `mzTab-version` in two.
    if key.starts_with("mzTab-version") {
        meta.mz_tab_version.read_cell(value)?;
        return Ok(());
    }
    if key.starts_with("mzTab-mode") {
        meta.mz_tab_mode.read_cell(value)?;
        return Ok(());
    }
    if key.starts_with("mzTab-type") {
        meta.mz_tab_type.read_cell(value)?;
        return Ok(());
    }
    if key.starts_with("mzTab-ID") {
        meta.mz_tab_id.read_cell(value)?;
        return Ok(());
    }
    match head {
        "title" => {
            meta.title.set(value);
            return Ok(());
        }
        "description" => {
            meta.description.set(value);
            return Ok(());
        }
        "false_discovery_rate" => {
            meta.false_discovery_rate.read_cell(value)?;
            return Ok(());
        }
        "quantification_method" => {
            meta.quantification_method.read_cell(value)?;
            return Ok(());
        }
        _ => {}
    }
    if second == Some("quantification_unit") {
        let target = match head {
            "protein" => Some(&mut meta.protein_quantification_unit),
            "peptide" => Some(&mut meta.peptide_quantification_unit),
            "small_molecule" => Some(&mut meta.small_molecule_quantification_unit),
            _ => None,
        };
        if let Some(target) = target {
            target.read_cell(value)?;
            return Ok(());
        }
    }
    if head == "colunit" {
        let target = match second.map(str::to_ascii_lowercase).as_deref() {
            Some("protein") => Some(&mut meta.colunit_protein),
            Some("peptide") => Some(&mut meta.colunit_peptide),
            Some("psm") => Some(&mut meta.colunit_psm),
            Some("small_molecule") => Some(&mut meta.colunit_small_molecule),
            _ => None,
        };
        if let Some(target) = target {
            if target.len() >= MzTabFile::MAX_COLUMNS {
                return Err(bad("MzTab colunit list exceeds its entry limit"));
            }
            target.push(value.to_owned());
        }
        return Ok(());
    }
    if head.starts_with("sample_processing[") {
        let index = index_of(head, "sample_processing[")?;
        let mut list = MzTabParameterList::default();
        list.read_cell(value)?;
        meta.sample_processing.insert(index, list);
        return Ok(());
    }
    for (label, slot) in [
        ("protein_search_engine_score[", 0u8),
        ("peptide_search_engine_score[", 1),
        ("psm_search_engine_score[", 2),
        ("smallmolecule_search_engine_score[", 3),
        ("nucleic_acid_search_engine_score[", 4),
        ("oligonucleotide_search_engine_score[", 5),
        ("osm_search_engine_score[", 6),
    ] {
        if !head.starts_with(label) {
            continue;
        }
        let index = index_of(head, label)?;
        let mut parameter = MzTabParameter::default();
        parameter.read_cell(value)?;
        let target = match slot {
            0 => &mut meta.protein_search_engine_score,
            1 => &mut meta.peptide_search_engine_score,
            2 => &mut meta.psm_search_engine_score,
            3 => &mut meta.smallmolecule_search_engine_score,
            4 => &mut meta.nucleic_acid_search_engine_score,
            5 => &mut meta.oligonucleotide_search_engine_score,
            _ => &mut meta.osm_search_engine_score,
        };
        target.insert(index, parameter);
        return Ok(());
    }
    if head.starts_with("instrument[") {
        let index = index_of(head, "instrument[")?;
        let entry = meta.instrument.entry(index).or_default();
        match second {
            Some("name") => entry.name.read_cell(value)?,
            Some("source") => entry.source.read_cell(value)?,
            Some("detector") => entry.detector.read_cell(value)?,
            Some(field) if field.starts_with("analyzer[") && fields.len() == 2 => {
                let inner = extract_bracket_index(field, "analyzer[")
                    .map_err(|_| at(format!("MzTab key {key:?} does not carry a valid index")))?;
                let mut parameter = MzTabParameter::default();
                parameter.read_cell(value)?;
                entry.analyzer.insert(inner, parameter);
            }
            _ => {}
        }
        return Ok(());
    }
    if head.starts_with("software[") {
        let index = index_of(head, "software[")?;
        let entry = meta.software.entry(index).or_default();
        if fields.len() == 1 {
            entry.software.read_cell(value)?;
        } else if let Some(field) = second.filter(|field| field.starts_with("setting[")) {
            if fields.len() == 2 {
                let inner = extract_bracket_index(field, "setting[")
                    .map_err(|_| at(format!("MzTab key {key:?} does not carry a valid index")))?;
                entry.setting.insert(inner, MzTabString::from_text(value));
            }
        }
        return Ok(());
    }
    if head.starts_with("publication[") {
        let index = index_of(head, "publication[")?;
        meta.publication
            .insert(index, MzTabString::from_text(value));
        return Ok(());
    }
    if head.starts_with("contact[") {
        let index = index_of(head, "contact[")?;
        let entry = meta.contact.entry(index).or_default();
        match second {
            Some("name") => entry.name.set(value),
            Some("affiliation") => entry.affiliation.set(value),
            Some("email") => entry.email.set(value),
            _ => {}
        }
        return Ok(());
    }
    if head.starts_with("uri[") {
        let index = index_of(head, "uri[")?;
        meta.uri.insert(index, MzTabString::from_text(value));
        return Ok(());
    }
    for label in ["variable_mod[", "fixed_mod["] {
        if !head.starts_with(label) {
            continue;
        }
        let index = index_of(head, label)?;
        let entry = if label == "variable_mod[" {
            meta.variable_mod.entry(index).or_default()
        } else {
            meta.fixed_mod.entry(index).or_default()
        };
        match second {
            None => entry.modification.read_cell(value)?,
            Some("site") => entry.site.set(value),
            Some("position") => entry.position.set(value),
            Some(_) => {}
        }
        return Ok(());
    }
    if head.starts_with("ms_run[") {
        let index = index_of(head, "ms_run[")?;
        let entry = meta.ms_run.entry(index).or_default();
        match second {
            Some("format") => entry.format.read_cell(value)?,
            Some("location") => entry.location.set(value),
            Some("id_format") => entry.id_format.read_cell(value)?,
            Some("fragmentation_method") => entry.fragmentation_method.read_cell(value)?,
            _ => {}
        }
        return Ok(());
    }
    if head.starts_with("custom[") {
        let index = index_of(head, "custom[")?;
        let mut parameter = MzTabParameter::default();
        parameter.read_cell(value)?;
        meta.custom.insert(index, parameter);
        return Ok(());
    }
    if head.starts_with("sample[") {
        let index = index_of(head, "sample[")?;
        let entry = meta.sample.entry(index).or_default();
        if second == Some("description") {
            entry.description.set(value);
            return Ok(());
        }
        for label in ["species[", "tissue[", "cell_type[", "disease[", "custom["] {
            let Some(field) = second.filter(|field| field.starts_with(label)) else {
                continue;
            };
            let inner = extract_bracket_index(field, label)
                .map_err(|_| at(format!("MzTab key {key:?} does not carry a valid index")))?;
            let mut parameter = MzTabParameter::default();
            parameter.read_cell(value)?;
            let target = match label {
                "species[" => &mut entry.species,
                "tissue[" => &mut entry.tissue,
                "cell_type[" => &mut entry.cell_type,
                "disease[" => &mut entry.disease,
                _ => &mut entry.custom,
            };
            target.insert(inner, parameter);
            return Ok(());
        }
        return Ok(());
    }
    if head.starts_with("assay[") {
        let index = index_of(head, "assay[")?;
        if second == Some("ms_run_ref") {
            let refs = parse_reference_list(value, "ms_run", line)?;
            meta.assay.entry(index).or_default().ms_run_ref = refs;
            return Ok(());
        }
        let entry = meta.assay.entry(index).or_default();
        match second {
            Some("quantification_reagent") => entry.quantification_reagent.read_cell(value)?,
            Some("sample_ref") => entry.sample_ref.set(value),
            Some(field) if field.starts_with("quantification_mod[") => {
                let inner = extract_bracket_index(field, "quantification_mod[")
                    .map_err(|_| at(format!("MzTab key {key:?} does not carry a valid index")))?;
                let modification = entry.quantification_mod.entry(inner).or_default();
                match third {
                    None => modification.modification.read_cell(value)?,
                    Some("site") => modification.site.set(value),
                    Some("position") => modification.position.set(value),
                    Some(_) => {}
                }
            }
            _ => {}
        }
        return Ok(());
    }
    if head.starts_with("cv[") {
        let index = index_of(head, "cv[")?;
        let entry = meta.cv.entry(index).or_default();
        match second {
            Some("label") => entry.label.set(value),
            Some("full_name") => entry.full_name.set(value),
            Some("version") => entry.version.set(value),
            Some("url") => entry.url.set(value),
            _ => {}
        }
        return Ok(());
    }
    if head.starts_with("study_variable[") {
        let index = index_of(head, "study_variable[")?;
        match second {
            Some("assay_refs") => {
                let refs = parse_reference_list(value, "assay", line)?;
                meta.study_variable.entry(index).or_default().assay_refs = refs;
            }
            Some("sample_refs") => {
                let refs = parse_reference_list(value, "sample", line)?;
                meta.study_variable.entry(index).or_default().sample_refs = refs;
            }
            Some("description") => meta
                .study_variable
                .entry(index)
                .or_default()
                .description
                .set(value),
            _ => {}
        }
        return Ok(());
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Reading — the document
// ---------------------------------------------------------------------------

impl MzTabFile {
    /// Read the MzTab document at `path`.
    ///
    /// Source `MzTabFile::load(const std::string&, MzTab&)`, which fills an
    /// out-parameter; this returns the document, so a failure cannot leave a
    /// half-filled one behind. A gzip- or bzip2-compressed file is decompressed
    /// transparently, which the source does not do.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be opened or read, and everything
    /// [`MzTabFile::load_reader`] reports.
    pub fn load(&self, path: impl AsRef<Path>) -> Result<MzTab> {
        let reader = path_io::open(path.as_ref())?;
        self.load_reader(reader)
    }

    /// Read an MzTab document from `reader`.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] for a data line with fewer than three tab-separated
    /// cells, for a cell whose text its column's type rejects, and for a
    /// bracketed index that is absent, not a positive integer, or above
    /// [`MzTabFile::MAX_INDEX`]. [`Error::InvalidValue`] when the input exceeds
    /// [`MzTabFile::MAX_LINES`], [`MzTabFile::MAX_BYTES`] or
    /// [`MzTabFile::MAX_COLUMNS`], or when a cell exceeds the data model's own
    /// ceilings. [`Error::Io`] for a read failure or non-UTF-8 input.
    pub fn load_reader(&self, reader: impl BufRead) -> Result<MzTab> {
        Ok(self.load_reporting(reader)?.0)
    }

    /// Read an MzTab document from `text`.
    ///
    /// Native convenience; the source can only read from a file.
    ///
    /// # Errors
    ///
    /// As [`MzTabFile::load_reader`].
    pub fn load_str(&self, text: &str) -> Result<MzTab> {
        self.load_reader(text.as_bytes())
    }

    /// Read an MzTab document from `reader`, together with the diagnostics the
    /// source prints to `std::cout`.
    ///
    /// Each diagnostic names a column the MzTab specification makes mandatory
    /// for this document's `mzTab-mode` and `mzTab-type` and which the section's
    /// header did not declare, or a data row shorter than its header. None is an
    /// error: the source reports and carries on, and so does this. The list
    /// holds at most 1024 entries.
    ///
    /// # Errors
    ///
    /// As [`MzTabFile::load_reader`].
    pub fn load_reporting(&self, mut reader: impl BufRead) -> Result<(MzTab, Vec<String>)> {
        let limits = Limits::default();
        let mut document = MzTab::default();
        let mut diagnostics: Vec<String> = Vec::new();
        let mut protein = Columns::default();
        let mut peptide = Columns::default();
        let mut psm = Columns::default();
        let mut small_molecule = Columns::default();
        let mut nucleic_acid = Columns::default();
        let mut oligonucleotide = Columns::default();
        let mut osm = Columns::default();
        let mut empty_rows: Vec<usize> = Vec::new();
        let mut comment_rows: BTreeMap<usize, String> = BTreeMap::new();
        let mut buffer = String::new();
        let mut index = 0usize;
        let mut bytes = 0usize;

        while TextFile::get_line_with_limits(&mut reader, &mut buffer, &limits)? {
            let line = index;
            index = index
                .checked_add(1)
                .filter(|&next| next <= Self::MAX_LINES)
                .ok_or_else(|| bad("MzTab input exceeds its line limit"))?;
            bytes = bytes
                .checked_add(buffer.len())
                .filter(|&total| total <= Self::MAX_BYTES)
                .ok_or_else(|| bad("MzTab input exceeds its byte limit"))?;

            let text = trim_source(&buffer);
            // The source records every line shorter than three bytes as an empty
            // row, so a one- or two-character line is discarded rather than
            // parsed.
            if text.len() < 3 {
                if empty_rows.len() >= Self::MAX_LINES {
                    return Err(bad("MzTab input exceeds its line limit"));
                }
                empty_rows.push(line);
                continue;
            }
            // `str::get` yields None when byte three falls inside a character,
            // so a line that begins with a multi-byte character can never be
            // mistaken for a section tag. The source takes the first three bytes
            // unconditionally.
            let tag = text.get(..3);
            if tag == Some(COMMENT) {
                comment_rows.insert(line, text.to_owned());
                continue;
            }
            if column_count(text) > Self::MAX_COLUMNS {
                return Err(bad("MzTab line exceeds its column limit"));
            }
            let cells: Vec<&str> = text.split('\t').collect();
            if cells.len() < 3 {
                return Err(parse_error(
                    line + 1,
                    format!(
                        "error parsing MzTab line: {text:?}. Did you forget to use tabulator as separator?"
                    ),
                ));
            }
            match tag {
                Some(METADATA) => read_metadata_line(&mut document.meta_data, &cells, line)?,
                Some("PRH") => {
                    protein = parse_header(PROTEIN_SINGLES, PROTEIN_FAMILIES, &cells, line)?;
                }
                Some("PRT") => {
                    read_protein_row(&mut protein, &mut document, &mut diagnostics, &cells, line)?
                }
                Some("PEH") => {
                    peptide = parse_header(PEPTIDE_SINGLES, PEPTIDE_FAMILIES, &cells, line)?;
                }
                Some("PEP") => {
                    read_peptide_row(&peptide, &mut document, &mut diagnostics, &cells, line)?;
                }
                Some("PSH") => {
                    psm = parse_header(PSM_SINGLES, PSM_FAMILIES, &cells, line)?;
                }
                Some("PSM") => {
                    read_psm_row(&psm, &mut document, &mut diagnostics, &cells, line)?;
                }
                Some("SMH") => {
                    small_molecule = parse_header(
                        SMALL_MOLECULE_SINGLES,
                        SMALL_MOLECULE_FAMILIES,
                        &cells,
                        line,
                    )?;
                }
                Some("SML") => read_small_molecule_row(
                    &small_molecule,
                    &mut document,
                    &mut diagnostics,
                    &cells,
                    line,
                )?,
                Some("NUH") => {
                    nucleic_acid =
                        parse_header(NUCLEIC_ACID_SINGLES, NUCLEIC_ACID_FAMILIES, &cells, line)?;
                }
                Some("NUC") => read_nucleic_acid_row(
                    &nucleic_acid,
                    &mut document,
                    &mut diagnostics,
                    &cells,
                    line,
                )?,
                Some("OLH") => {
                    oligonucleotide = parse_header(
                        OLIGONUCLEOTIDE_SINGLES,
                        OLIGONUCLEOTIDE_FAMILIES,
                        &cells,
                        line,
                    )?;
                }
                Some("OLI") => read_oligonucleotide_row(
                    &oligonucleotide,
                    &mut document,
                    &mut diagnostics,
                    &cells,
                    line,
                )?,
                Some("OSH") => {
                    osm = parse_header(OSM_SINGLES, OSM_FAMILIES, &cells, line)?;
                }
                Some("OSM") => {
                    read_osm_row(&osm, &mut document, &mut diagnostics, &cells, line)?;
                }
                _ => {}
            }
        }

        document.empty_rows = empty_rows;
        document.comment_rows = comment_rows;
        Ok((document, diagnostics))
    }
}
