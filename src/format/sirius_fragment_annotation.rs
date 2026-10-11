// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Fragment annotations from a SIRIUS project space, as annotated spectra.
//!
//! Port of `FORMAT/DATAACCESS/SiriusFragmentAnnotation.h`
//! (`SiriusFragmentAnnotation.cpp`). See `docs/ID_MAPPER_SUPPORT.md` (section
//! "SiriusFragmentAnnotation") for the API mapping.
//!
//! One SIRIUS compound directory holds
//!
//! - `spectrum.ms`, the input OpenMS wrote for SIRIUS, whose `##n_id`, `##m_id`
//!   and `##fid` comment lines carry the native IDs, the SiriusExport m_ids and
//!   the OpenMS feature ID;
//! - `formula_candidates.tsv`, one row per candidate formula with its
//!   `formulaRank`, `molecularFormula`, `adduct` and `explainedIntensity`;
//! - `spectra/<formula>_<adduct>.tsv` (targets) and, after PASSATUTTO,
//!   `decoys/<formula>_<adduct>.tsv`, one row per annotated fragment peak.
//!
//! The functions read those files into [`MSSpectrum`](crate::kernel::MSSpectrum)s whose peaks are the
//! annotated fragments.

use super::csv::{CsvFile, ReadOptions as CsvReadOptions};
use crate::concept::constants::user_param::{
    SIRIUS_ANNOTATED_ADDUCT, SIRIUS_ANNOTATED_SUMFORMULA, SIRIUS_DECOY, SIRIUS_EXACTMASS,
    SIRIUS_EXPLANATION, SIRIUS_FEATURE_ID, SIRIUS_MZ, SIRIUS_PEAKMZ, SIRIUS_SCORE,
};
use crate::kernel::{DataArray, MSSpectrum, Peak1D};
use crate::metadata::MetaValue;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::path::Path;

/// Largest file this module reads (`spectrum.ms` or a fragment table).
/// Native ceiling; the source reads without one.
pub const MAX_FILE_BYTES: u64 = 256 * 1024 * 1024;

/// Target and decoy spectra of one SIRIUS compound directory, the source
/// nested class `SiriusFragmentAnnotation::SiriusTargetDecoySpectra`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SiriusTargetDecoySpectra {
    /// The annotated target spectrum.
    pub target: MSSpectrum,
    /// The annotated decoy spectrum, or — without decoys — an empty spectrum
    /// carrying only the concatenated native IDs and m_ids.
    pub decoy: MSSpectrum,
}

impl SiriusTargetDecoySpectra {
    /// A pair from a target and a decoy spectrum.
    pub fn new(target: MSSpectrum, decoy: MSSpectrum) -> Self {
        Self { target, decoy }
    }
}

fn read_text(path: &Path) -> Result<Option<String>> {
    let Ok(metadata) = std::fs::metadata(path) else {
        return Ok(None); // the source's failed ifstream: treated as absent
    };
    if metadata.len() > MAX_FILE_BYTES {
        return Err(Error::InvalidValue(format!(
            "SiriusFragmentAnnotation: '{}' exceeds {MAX_FILE_BYTES} bytes",
            path.display()
        )));
    }
    match std::fs::read(path) {
        Ok(bytes) => Ok(Some(String::from_utf8_lossy(&bytes).into_owned())),
        Err(_) => Ok(None),
    }
}

/// The values of every `spectrum.ms` line starting with `prefix`, as the
/// source's `getline` loop collects them (a `\r` of a CRLF file stays).
fn spectrum_ms_values(workspace: &Path, prefix: &str) -> Result<Vec<String>> {
    let Some(text) = read_text(&workspace.join("spectrum.ms"))? else {
        return Ok(Vec::new());
    };
    Ok(text
        .split('\n')
        .filter_map(|line| line.strip_prefix(prefix).map(str::to_owned))
        .collect())
}

/// The `##n_id` values of `spectrum.ms` joined by `|`, the source protected
/// `extractConcatNativeIDsFromSiriusMS_`. Empty when the file or the lines are
/// missing (the source then logs a warning).
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the file exceeds [`MAX_FILE_BYTES`].
pub fn extract_concat_native_ids(workspace: impl AsRef<Path>) -> Result<String> {
    Ok(spectrum_ms_values(workspace.as_ref(), "##n_id ")?.join("|"))
}

/// The `##m_id` values of `spectrum.ms` joined by `|`, the source protected
/// `extractConcatMIDsFromSiriusMS_`. An m_id is a native ID plus an index
/// incremented per possible identification (accurate mass search).
///
/// # Errors
///
/// As [`extract_concat_native_ids`].
pub fn extract_concat_m_ids(workspace: impl AsRef<Path>) -> Result<String> {
    Ok(spectrum_ms_values(workspace.as_ref(), "##m_id ")?.join("|"))
}

/// The first `##fid` value of `spectrum.ms` (the original OpenMS feature ID),
/// the source protected `extractFeatureIDFromSiriusMS_`. Empty when absent: an
/// fid only exists when the SIRIUS run started from an OpenMS feature rather
/// than a single MS2 spectrum.
///
/// # Errors
///
/// As [`extract_concat_native_ids`].
pub fn extract_feature_id(workspace: impl AsRef<Path>) -> Result<String> {
    Ok(spectrum_ms_values(workspace.as_ref(), "##fid ")?
        .into_iter()
        .next()
        .unwrap_or_default())
}

/// Column name to column index from the first row of `csv`, the source static
/// `extract_columnname_to_columnindex`. A repeated name keeps its first
/// index, as `std::map::insert` does.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `csv` has no rows (the source reads
/// row 0 unchecked) or the row cannot be split.
pub fn extract_columnname_to_columnindex(csv: &CsvFile) -> Result<BTreeMap<String, usize>> {
    let (_, header) = csv.row(0)?;
    let mut map = BTreeMap::new();
    for (index, name) in header.into_iter().enumerate() {
        map.entry(name).or_insert(index);
    }
    Ok(map)
}

/// Source `StringUtils::toInt32`: surrounding ASCII whitespace, one optional
/// `+`, then a complete decimal `i32`.
fn source_int32(text: &str) -> Option<i32> {
    let trimmed = text.trim_matches([' ', '\t', '\n', '\r', '\x0b', '\x0c']);
    let body = trimmed.strip_prefix('+').unwrap_or(trimmed);
    let digits = body.strip_prefix('-').unwrap_or(body);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    body.parse().ok()
}

/// Source `StringUtils::toDouble` / `toFloat` for the decimal numbers SIRIUS
/// writes: surrounding ASCII whitespace, one optional `+`, then a complete
/// number; an out-of-range finite literal is refused, as `std::from_chars`
/// refuses it. Rust's parser, not `from_chars`, reads the digits.
fn source_real<T: std::str::FromStr + Copy + Into<f64>>(text: &str) -> Option<T> {
    let trimmed = text.trim_matches([' ', '\t', '\n', '\r', '\x0b', '\x0c']);
    let body = trimmed.strip_prefix('+').unwrap_or(trimmed);
    if body.starts_with('+') {
        return None;
    }
    let value: T = body.parse().ok()?;
    let lower = body.to_ascii_lowercase();
    if value.into().is_infinite() && !lower.contains("inf") {
        return None;
    }
    Some(value)
}

fn conversion(what: &str, text: &str) -> Error {
    Error::InvalidValue(format!(
        "SiriusFragmentAnnotation: could not convert '{text}' to {what}"
    ))
}

/// The rows of `formula_candidates.tsv` as `(rank, filename, score)`.
fn formula_candidates(workspace: &Path) -> Result<Vec<(usize, String, f64)>> {
    let path = workspace.join("formula_candidates.tsv");
    if !path.is_file() {
        return Ok(Vec::new());
    }
    let csv = CsvFile::from_path(
        &path,
        &CsvReadOptions {
            separator: b'\t',
            ..CsvReadOptions::default()
        },
    )?;
    if csv.row_count() == 0 {
        return Ok(Vec::new());
    }
    let columns = extract_columnname_to_columnindex(&csv)?;
    let column = |name: &str| -> Result<usize> {
        columns.get(name).copied().ok_or_else(|| {
            Error::InvalidValue(format!(
                "SiriusFragmentAnnotation: formula_candidates.tsv has no '{name}' column"
            ))
        })
    };
    let (rank_col, formula_col, adduct_col, score_col) = (
        column("formulaRank")?,
        column("molecularFormula")?,
        column("adduct")?,
        column("explainedIntensity")?,
    );
    let mut rows = Vec::new();
    for i in 1..csv.row_count() {
        let (_, fields) = csv.row(i)?;
        let field = |index: usize| -> Result<&str> {
            fields.get(index).map(String::as_str).ok_or_else(|| {
                Error::InvalidValue(format!(
                    "SiriusFragmentAnnotation: formula_candidates.tsv row {i} is too short"
                ))
            })
        };
        let rank_text = field(rank_col)?;
        // Source `Size(toInt32(...))`: a negative rank wraps to a huge key.
        let rank = source_int32(rank_text).ok_or_else(|| conversion("an integer", rank_text))?;
        let rank = rank as i64 as u64 as usize;
        let adduct: String = field(adduct_col)?
            .chars()
            .filter(|c| !c.is_ascii_whitespace() && *c != '\x0b')
            .collect();
        let filename = format!("{}_{adduct}.tsv", field(formula_col)?);
        let score_text = field(score_col)?;
        let score =
            source_real::<f64>(score_text).ok_or_else(|| conversion("a double", score_text))?;
        rows.push((rank, filename, score));
    }
    Ok(rows)
}

/// Candidate rank to fragment-table file name (`<formula>_<adduct>.tsv`, the
/// adduct with its whitespace removed) from `formula_candidates.tsv`, the
/// source protected `extractCompoundRankingAndFilename_`. The first row wins
/// for a repeated rank (`std::map::emplace`). Empty when the file is missing.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when a required column is missing, a row
/// is too short, or a rank or score does not convert (the source throws
/// `std::out_of_range` or `Exception::ConversionError`).
pub fn extract_compound_ranking_and_filename(
    workspace: impl AsRef<Path>,
) -> Result<BTreeMap<usize, String>> {
    let mut map = BTreeMap::new();
    for (rank, filename, _) in formula_candidates(workspace.as_ref())? {
        map.entry(rank).or_insert(filename);
    }
    Ok(map)
}

/// Candidate rank to score (`explainedIntensity`) from
/// `formula_candidates.tsv`, the source protected
/// `extractCompoundRankingAndScore_`.
///
/// # Errors
///
/// As [`extract_compound_ranking_and_filename`].
pub fn extract_compound_ranking_and_score(
    workspace: impl AsRef<Path>,
) -> Result<BTreeMap<usize, f64>> {
    let mut map = BTreeMap::new();
    for (rank, _, score) in formula_candidates(workspace.as_ref())? {
        map.entry(rank).or_insert(score);
    }
    Ok(map)
}

/// Read the annotated fragment spectra of the best `max_rank` candidates of one
/// SIRIUS compound directory, the source static
/// `extractAnnotationsFromSiriusFile`.
///
/// With `decoy`, the fragment tables come from `decoys/` (SIRIUS/PASSATUTTO)
/// and `use_exact_mass` is forced off; otherwise from `spectra/`. When that
/// directory does not exist the result is empty. Otherwise one spectrum is
/// returned per rank from 1 to `min(max_rank, highest rank)`; spectrum `i`
/// carries:
///
/// - the native ID and name: the concatenated `##n_id` and `##m_id` values,
///   with `_i` appended for ranks above 1;
/// - meta values [`SIRIUS_PEAKMZ`] (`exact_mass` with `use_exact_mass`, else
///   `mz`), [`SIRIUS_ANNOTATED_SUMFORMULA`] and [`SIRIUS_ANNOTATED_ADDUCT`]
///   (the file name split at its last `_` and `.`), [`SIRIUS_DECOY`] (integer
///   1 or 0 — the source stores the `bool`, which `DataValue` holds as an
///   integer), [`SIRIUS_SCORE`] (the rank's `explainedIntensity`) and, when
///   `spectrum.ms` has one, [`SIRIUS_FEATURE_ID`].
///
/// If the rank's fragment table exists, MS level 2 and its peaks are set: one
/// [`Peak1D`] per row with the m/z from column 0 (or the exact mass from
/// column 3 with `use_exact_mass`) and the intensity from column 1; for
/// targets a float data array with the other of the two masses (named
/// [`SIRIUS_EXACTMASS`] or [`SIRIUS_MZ`]); and two string data arrays, the
/// explanation (next-to-last column, named [`SIRIUS_EXPLANATION`]) and the
/// ionization (last column, unnamed, as in the source). For decoys the float
/// array is present but empty.
///
/// | peak m/z | intensity | float array | explanation | ionization |
/// |---|---|---|---|---|
/// | 56.050855 | 20794.85 | 56.049476 | C3H5N | `[M + H]+` |
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `formula_candidates.tsv` is malformed
/// (see [`extract_compound_ranking_and_filename`]), a rank between 1 and the
/// highest rank is missing (the source's `std::map::at` throws), a fragment
/// row has fewer columns than it reads or an empty line (the source indexes
/// out of range), a number does not convert, or a file exceeds
/// [`MAX_FILE_BYTES`]. A missing or empty `formula_candidates.tsv` gives an
/// empty result; the source dereferences the end of its empty rank map first.
pub fn extract_annotations_from_sirius_file(
    workspace: impl AsRef<Path>,
    max_rank: usize,
    decoy: bool,
    use_exact_mass: bool,
) -> Result<Vec<MSSpectrum>> {
    let workspace = workspace.as_ref();
    let use_exact_mass = use_exact_mass && !decoy;
    let dir = workspace.join(if decoy { "decoys" } else { "spectra" });
    if !dir.is_dir() {
        return Ok(Vec::new());
    }
    let native_ids = extract_concat_native_ids(workspace)?;
    let m_ids = extract_concat_m_ids(workspace)?;
    let fid = extract_feature_id(workspace)?;
    let rank_filename = extract_compound_ranking_and_filename(workspace)?;
    let rank_score = extract_compound_ranking_and_score(workspace)?;
    let Some(&max_found_rank) = rank_filename.keys().next_back() else {
        return Ok(Vec::new());
    };
    if rank_score.is_empty() {
        return Ok(Vec::new());
    }
    let count = max_rank.min(max_found_rank);
    let mut result = Vec::new();
    for i in 1..=count {
        let suffix = if i > 1 {
            format!("_{i}")
        } else {
            String::new()
        };
        let filename = rank_filename.get(&i).ok_or_else(|| {
            Error::InvalidValue(format!(
                "SiriusFragmentAnnotation: formula_candidates.tsv has no rank {i}"
            ))
        })?;
        let score = *rank_score.get(&i).ok_or_else(|| {
            Error::InvalidValue(format!(
                "SiriusFragmentAnnotation: formula_candidates.tsv has no rank {i}"
            ))
        })?;
        let mut spectrum = MSSpectrum {
            native_id: format!("{native_ids}{suffix}"),
            name: format!("{m_ids}{suffix}"),
            ..MSSpectrum::default()
        };
        let meta = &mut spectrum.metadata;
        meta.insert(
            SIRIUS_PEAKMZ.into(),
            MetaValue::from(if use_exact_mass {
                SIRIUS_EXACTMASS
            } else {
                SIRIUS_MZ
            }),
        );
        // std::string::npos + 1 wraps to 0 in the source when there is no '_'.
        let underscore = filename.rfind('_');
        let sumformula = match underscore {
            Some(u) => &filename[..u],
            None => filename.as_str(),
        };
        let adduct_start = underscore.map_or(0, |u| u + 1);
        let adduct_end = filename
            .rfind('.')
            .unwrap_or(filename.len())
            .max(adduct_start);
        let adduct = &filename[adduct_start..adduct_end];
        meta.insert(
            SIRIUS_ANNOTATED_SUMFORMULA.into(),
            MetaValue::from(sumformula),
        );
        meta.insert(SIRIUS_ANNOTATED_ADDUCT.into(), MetaValue::from(adduct));
        meta.insert(SIRIUS_DECOY.into(), MetaValue::from(i64::from(decoy)));
        meta.insert(SIRIUS_SCORE.into(), MetaValue::try_from(score)?);
        if !fid.is_empty() {
            meta.insert(SIRIUS_FEATURE_ID.into(), MetaValue::from(fid.as_str()));
        }
        if let Some(text) = read_text(&dir.join(filename))? {
            read_fragments(&text, decoy, use_exact_mass, &mut spectrum)?;
        }
        result.push(spectrum);
    }
    Ok(result)
}

fn read_fragments(
    text: &str,
    decoy: bool,
    use_exact_mass: bool,
    spectrum: &mut MSSpectrum,
) -> Result<()> {
    let (main_col, extra_col) = if use_exact_mass { (3, 0) } else { (0, 3) };
    let mut peaks = Vec::new();
    let mut extra = DataArray::<f32> {
        name: if use_exact_mass {
            SIRIUS_MZ
        } else {
            SIRIUS_EXACTMASS
        }
        .into(),
        ..DataArray::default()
    };
    let mut explanations = DataArray::<String> {
        name: SIRIUS_EXPLANATION.into(),
        ..DataArray::default()
    };
    let mut ionization = DataArray::<String>::default();
    // getline: lines split at '\n'; a final newline ends the last line.
    let mut lines: Vec<&str> = text.split('\n').collect();
    if text.ends_with('\n') {
        lines.pop();
    }
    for (number, line) in lines.iter().enumerate().skip(1) {
        let fields: Vec<&str> = if line.is_empty() {
            Vec::new()
        } else {
            line.split('\t').collect()
        };
        let needed = if decoy {
            main_col.max(1) + 1
        } else {
            main_col.max(extra_col).max(1) + 1
        };
        if fields.len() < needed.max(2) {
            return Err(Error::InvalidValue(format!(
                "SiriusFragmentAnnotation: fragment table line {} has {} columns",
                number + 1,
                fields.len()
            )));
        }
        if !decoy {
            let value = source_real::<f32>(fields[extra_col])
                .ok_or_else(|| conversion("a float", fields[extra_col]))?;
            extra.data.push(value);
        }
        let mz = source_real::<f64>(fields[main_col])
            .ok_or_else(|| conversion("a double", fields[main_col]))?;
        let intensity =
            source_real::<f32>(fields[1]).ok_or_else(|| conversion("a float", fields[1]))?;
        peaks.push(Peak1D { mz, intensity });
        explanations.data.push(fields[fields.len() - 2].to_owned());
        ionization.data.push(fields[fields.len() - 1].to_owned());
    }
    spectrum.ms_level = 2;
    spectrum.peaks = peaks;
    spectrum.float_data_arrays.push(extra);
    spectrum.string_data_arrays.push(explanations);
    spectrum.string_data_arrays.push(ionization);
    Ok(())
}

/// Extract and resolve the SIRIUS target (and optionally decoy) annotations
/// of several compound directories, keyed by native ID, the source static
/// `extractAndResolveSiriusAnnotations`.
///
/// For each directory the rank-1 target spectrum is read
/// ([`extract_annotations_from_sirius_file`]); directories without one are
/// skipped. With `decoy_generation` the rank-1 decoy spectrum is read too;
/// otherwise — or when there is none — the decoy is an empty spectrum carrying
/// only the concatenated native IDs and m_ids (the source's
/// backwards-compatibility filler). A target whose [`SIRIUS_SCORE`] is at
/// least `score_threshold` (0 to 1) is kept; when several directories share
/// a (concatenated) native ID, the one with the higher score wins, a later one
/// on a tie. The result is ordered by native ID.
///
/// # Errors
///
/// As [`extract_annotations_from_sirius_file`].
pub fn extract_and_resolve_sirius_annotations(
    workspace_subdirs: &[impl AsRef<Path>],
    score_threshold: f64,
    use_exact_mass: bool,
    decoy_generation: bool,
) -> Result<Vec<SiriusTargetDecoySpectra>> {
    let mut by_native_id: BTreeMap<String, SiriusTargetDecoySpectra> = BTreeMap::new();
    let score_of = |spectrum: &MSSpectrum| -> f64 {
        spectrum
            .metadata
            .get(SIRIUS_SCORE)
            .and_then(|v| v.as_f64().ok())
            .unwrap_or(f64::NAN)
    };
    for subdir in workspace_subdirs {
        let subdir = subdir.as_ref();
        let Some(best) = extract_annotations_from_sirius_file(subdir, 1, false, use_exact_mass)?
            .into_iter()
            .next()
        else {
            continue;
        };
        let decoys = if decoy_generation {
            extract_annotations_from_sirius_file(subdir, 1, true, use_exact_mass)?
        } else {
            Vec::new()
        };
        let decoy = match decoys.into_iter().next() {
            Some(decoy) => decoy,
            None => MSSpectrum {
                native_id: extract_concat_native_ids(subdir)?,
                name: extract_concat_m_ids(subdir)?,
                ..MSSpectrum::default()
            },
        };
        let score = score_of(&best);
        if score >= score_threshold {
            match by_native_id.get_mut(&best.native_id) {
                Some(existing) => {
                    if score >= score_of(&existing.target) {
                        *existing = SiriusTargetDecoySpectra::new(best, decoy);
                    }
                }
                None => {
                    by_native_id.insert(
                        best.native_id.clone(),
                        SiriusTargetDecoySpectra::new(best, decoy),
                    );
                }
            }
        }
    }
    Ok(by_native_id.into_values().collect())
}

#[cfg(test)]
mod tests {
    use super::{source_int32, source_real};

    #[test]
    fn numeric_conversions_follow_string_utils() {
        assert_eq!(source_int32(" +3 "), Some(3));
        assert_eq!(source_int32("3x"), None);
        assert_eq!(source_real::<f64>(" 0.5\t"), Some(0.5));
        assert_eq!(source_real::<f64>("+1e3"), Some(1000.0));
        assert_eq!(source_real::<f64>("1e400"), None);
        assert_eq!(source_real::<f64>("++1"), None);
        assert_eq!(source_real::<f32>("1e39"), None);
        assert!(source_real::<f64>("inf").unwrap().is_infinite());
    }
}
