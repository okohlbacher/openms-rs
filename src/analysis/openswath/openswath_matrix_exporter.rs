// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Building and writing of OpenSWATH quantification matrices.
//!
//! Ports `ANALYSIS/OPENSWATH/OpenSwathMatrixExporter.h` with its `.cpp`.
//! [`OpenSwathMatrixExporter::build_matrix`] follows the PyProphet matrix
//! export logic of the source: the best peak group per run and precursor is
//! selected first (lowest q-value, then highest intensity, then lowest
//! feature id), then summarized to precursor, peptide, protein or gene level,
//! with an optional median or median-of-medians normalization.
//! [`OpenSwathMatrixExporter::write_matrix`] writes the result as TSV.
//!
//! As for the results exporter, the Parquet arm of the source writes an Arrow
//! table: [`OpenSwathMatrixExporter::table`] returns the typed columns and
//! [`write_matrix`](crate::analysis::openswath::openswath_matrix_exporter::OpenSwathMatrixExporter::write_matrix)
//! with [`OpenSwathExportFileFormat::Parquet`] writes them with the `arrow`
//! feature and returns [`Error::Unsupported`] without it, after the checks of
//! the source.
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_EXPORT_SUPPORT.md).
//!
//! # Example
//!
//! ```
//! use openms::analysis::openswath::openswath_export_config::OpenSwathMatrixExportConfig;
//! use openms::analysis::openswath::openswath_export_data::{
//!     OpenSwathExportRow, OpenSwathMatrixLevel,
//! };
//! use openms::analysis::openswath::openswath_matrix_exporter::OpenSwathMatrixExporter;
//!
//! let row = |run_id, name: &str, intensity, m_score, feature_id| OpenSwathExportRow {
//!     run_id,
//!     run_name: name.into(),
//!     transition_group_id: "TG1".into(),
//!     sequence: "PEPTIDEK".into(),
//!     full_peptide_name: "PEPTIDEK".into(),
//!     protein_name: "PROT1".into(),
//!     charge: 2,
//!     intensity,
//!     m_score,
//!     feature_id,
//!     ..OpenSwathExportRow::default()
//! };
//! let rows = vec![
//!     row(1, "runA", 1000.0, 0.01, 11),
//!     row(1, "runA", 9999.0, 0.05, 12), // worse q-value: loses
//!     row(2, "runB", 2000.0, 0.02, 21),
//! ];
//! let config = OpenSwathMatrixExportConfig {
//!     level: OpenSwathMatrixLevel::Precursor,
//!     ..OpenSwathMatrixExportConfig::default()
//! };
//! let matrix = OpenSwathMatrixExporter::build_matrix(&rows, &config).unwrap();
//! assert_eq!(matrix.sample_column_names, ["runA", "runB"]);
//! assert_eq!(matrix.values, [[Some(1000.0), Some(2000.0)]]);
//! ```
//!
//! [`OpenSwathMatrixExporter::build_matrix`]: crate::analysis::openswath::openswath_matrix_exporter::OpenSwathMatrixExporter::build_matrix
//! [`OpenSwathMatrixExporter::write_matrix`]: crate::analysis::openswath::openswath_matrix_exporter::OpenSwathMatrixExporter::write_matrix
//! [`OpenSwathMatrixExporter::table`]: crate::analysis::openswath::openswath_matrix_exporter::OpenSwathMatrixExporter::table
//! [`OpenSwathExportFileFormat::Parquet`]: crate::analysis::openswath::openswath_export_data::OpenSwathExportFileFormat::Parquet
//! [`Error::Unsupported`]: crate::Error::Unsupported

use super::openswath_export_config::OpenSwathMatrixExportConfig;
use super::openswath_export_data::{
    OpenSwathExportFileFormat, OpenSwathExportRow, OpenSwathMatrixLevel,
    OpenSwathMatrixNormalization, OpenSwathQuantMatrix,
};
use super::openswath_results_exporter::{ExportColumn, ExportTable};
use crate::format::file_info::text_format::ostream_g;
use crate::math::source_sort::source_sort_by;
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// The builder and writer of OpenSWATH quantification matrices (source
/// `OpenSwathMatrixExporter`, a class of two static members).
#[derive(Clone, Copy, Debug, Default)]
pub struct OpenSwathMatrixExporter;

fn precondition(message: &str) -> Error {
    Error::InvalidValue(message.to_owned())
}

/// The median of `values`: NaN for none, the middle value for an odd count and
/// the mean of the two middle ones for an even count.
fn median(mut values: Vec<f64>) -> Result<f64> {
    if values.is_empty() {
        return Ok(f64::NAN);
    }
    source_sort_by(&mut values, |a, b| a < b)?;
    let mid = values.len() / 2;
    if values.len() % 2 == 1 {
        Ok(values[mid])
    } else {
        Ok((values[mid - 1] + values[mid]) / 2.0)
    }
}

/// The sample label of every run: its name, `RUN_ID_<id>` without one, and
/// `_RUN<id>` appended when the label is already taken by a run met earlier in
/// row order. Ordered by run id.
fn build_sample_names(rows: &[OpenSwathExportRow]) -> BTreeMap<i64, String> {
    let mut names: BTreeMap<i64, String> = BTreeMap::new();
    let mut used: BTreeSet<String> = BTreeSet::new();
    for row in rows {
        if names.contains_key(&row.run_id) {
            continue;
        }
        let mut label = if row.run_name.is_empty() {
            format!("RUN_ID_{}", row.run_id)
        } else {
            row.run_name.clone()
        };
        if used.contains(&label) {
            label.push_str(&format!("_RUN{}", row.run_id));
        }
        used.insert(label.clone());
        names.insert(row.run_id, label);
    }
    names
}

/// The best peak group of every (run, transition group): lowest `m_score`,
/// then highest `intensity`, then lowest `feature_id`.
fn select_best_peak_groups(rows: &[OpenSwathExportRow]) -> Result<Vec<&OpenSwathExportRow>> {
    let mut sorted: Vec<&OpenSwathExportRow> = rows.iter().collect();
    source_sort_by(&mut sorted, |lhs, rhs| {
        if lhs.run_id != rhs.run_id {
            return lhs.run_id < rhs.run_id;
        }
        if lhs.transition_group_id != rhs.transition_group_id {
            return lhs.transition_group_id < rhs.transition_group_id;
        }
        if lhs.m_score != rhs.m_score {
            return lhs.m_score < rhs.m_score;
        }
        if lhs.intensity != rhs.intensity {
            return lhs.intensity > rhs.intensity;
        }
        lhs.feature_id < rhs.feature_id
    })?;
    let mut best: Vec<&OpenSwathExportRow> = Vec::new();
    let mut current_group: &str = "";
    let mut current_run = i64::MIN;
    for row in sorted {
        if row.run_id != current_run || row.transition_group_id != current_group {
            best.push(row);
            current_run = row.run_id;
            current_group = &row.transition_group_id;
        }
    }
    Ok(best)
}

type PeptideKey = (String, String);

struct PeptideSummaryRow {
    sequence: String,
    full_peptide_name: String,
    protein_name: String,
    gene_name: String,
    sample_values: BTreeMap<i64, f64>,
}

/// `std::accumulate(values, 0.0) / size`.
fn mean(values: &[f64]) -> f64 {
    let mut sum = 0.0;
    for v in values {
        sum += v;
    }
    sum / values.len() as f64
}

fn summarize_peptides(
    best_rows: &[&OpenSwathExportRow],
    top_n: usize,
    consistent_top: bool,
) -> Result<Vec<PeptideSummaryRow>> {
    // `emplace` keeps the first value of a key
    let mut protein_by_key: BTreeMap<PeptideKey, String> = BTreeMap::new();
    let mut gene_by_key: BTreeMap<PeptideKey, String> = BTreeMap::new();
    for row in best_rows {
        let key = (row.sequence.clone(), row.full_peptide_name.clone());
        protein_by_key
            .entry(key.clone())
            .or_insert_with(|| row.protein_name.clone());
        gene_by_key
            .entry(key)
            .or_insert_with(|| row.gene_name.clone());
    }

    let mut per_run_mean: BTreeMap<PeptideKey, BTreeMap<i64, f64>> = BTreeMap::new();
    if consistent_top {
        let mut intensities_by_precursor: BTreeMap<PeptideKey, BTreeMap<&str, Vec<f64>>> =
            BTreeMap::new();
        for row in best_rows {
            intensities_by_precursor
                .entry((row.sequence.clone(), row.full_peptide_name.clone()))
                .or_default()
                .entry(row.transition_group_id.as_str())
                .or_default()
                .push(row.intensity);
        }

        let mut selected_precursors: BTreeMap<PeptideKey, BTreeSet<&str>> = BTreeMap::new();
        for (key, precursors) in &intensities_by_precursor {
            let mut medians: Vec<(&str, f64)> = Vec::with_capacity(precursors.len());
            for (precursor, values) in precursors {
                medians.push((precursor, median(values.clone())?));
            }
            source_sort_by(&mut medians, |lhs, rhs| lhs.1 > rhs.1)?;
            let selected = selected_precursors.entry(key.clone()).or_default();
            for (precursor, _) in medians.iter().take(top_n.min(medians.len())) {
                selected.insert(*precursor);
            }
        }

        let mut grouped: BTreeMap<(PeptideKey, i64), Vec<f64>> = BTreeMap::new();
        for row in best_rows {
            let key = (row.sequence.clone(), row.full_peptide_name.clone());
            let Some(selected) = selected_precursors.get(&key) else {
                continue;
            };
            if !selected.contains(row.transition_group_id.as_str()) {
                continue;
            }
            grouped
                .entry((key, row.run_id))
                .or_default()
                .push(row.intensity);
        }
        for ((key, run_id), values) in grouped {
            per_run_mean
                .entry(key)
                .or_default()
                .insert(run_id, mean(&values));
        }
    } else {
        let mut grouped_rows: BTreeMap<(PeptideKey, i64), Vec<&OpenSwathExportRow>> =
            BTreeMap::new();
        for row in best_rows {
            grouped_rows
                .entry((
                    (row.sequence.clone(), row.full_peptide_name.clone()),
                    row.run_id,
                ))
                .or_default()
                .push(row);
        }
        for ((key, run_id), mut group_rows) in grouped_rows {
            source_sort_by(&mut group_rows, |lhs, rhs| lhs.intensity > rhs.intensity)?;
            let keep = top_n.min(group_rows.len());
            let mut sum = 0.0;
            for row in group_rows.iter().take(keep) {
                sum += row.intensity;
            }
            per_run_mean
                .entry(key)
                .or_default()
                .insert(run_id, sum / keep as f64);
        }
    }

    let mut rows = Vec::with_capacity(per_run_mean.len());
    for (key, sample_values) in per_run_mean {
        rows.push(PeptideSummaryRow {
            protein_name: protein_by_key.get(&key).cloned().unwrap_or_default(),
            gene_name: gene_by_key.get(&key).cloned().unwrap_or_default(),
            sequence: key.0,
            full_peptide_name: key.1,
            sample_values,
        });
    }
    Ok(rows)
}

fn validate_matrix_dimensions(matrix: &OpenSwathQuantMatrix) -> Result<()> {
    if matrix.identifier_rows.len() != matrix.values.len() {
        return Err(precondition(
            "Matrix identifiers and values must have identical row counts.",
        ));
    }
    for (identifiers, values) in matrix.identifier_rows.iter().zip(&matrix.values) {
        if identifiers.len() != matrix.identifier_column_names.len() {
            return Err(precondition(
                "Identifier row width does not match identifier columns.",
            ));
        }
        if values.len() != matrix.sample_column_names.len() {
            return Err(precondition(
                "Value row width does not match sample columns.",
            ));
        }
    }
    Ok(())
}

fn apply_normalization(
    matrix: &mut OpenSwathQuantMatrix,
    normalization: OpenSwathMatrixNormalization,
) -> Result<()> {
    if normalization == OpenSwathMatrixNormalization::None {
        return Ok(());
    }
    let mut sample_medians = vec![f64::NAN; matrix.sample_column_names.len()];
    for (col, sample_median) in sample_medians.iter_mut().enumerate() {
        let values: Vec<f64> = matrix
            .values
            .iter()
            .filter_map(|row| row.get(col).copied().flatten())
            .collect();
        *sample_median = median(values)?;
    }
    let mut global_median = f64::NAN;
    if normalization == OpenSwathMatrixNormalization::MedianMedian {
        let finite: Vec<f64> = sample_medians
            .iter()
            .copied()
            .filter(|m| m.is_finite())
            .collect();
        global_median = median(finite)?;
    }
    for row in &mut matrix.values {
        for (col, cell) in row.iter_mut().enumerate() {
            let Some(value) = *cell else { continue };
            let sample_median = sample_medians.get(col).copied().unwrap_or(f64::NAN);
            if !sample_median.is_finite() || sample_median == 0.0 {
                continue;
            }
            let mut value = value / sample_median;
            if normalization == OpenSwathMatrixNormalization::MedianMedian
                && global_median.is_finite()
            {
                value *= global_median;
            }
            *cell = Some(value);
        }
    }
    Ok(())
}

impl OpenSwathMatrixExporter {
    /// Build a matrix from filtered export rows (source `buildMatrix`).
    ///
    /// 1. The sample columns are the runs, one per run id in increasing id
    ///    order, labelled by the run name (`RUN_ID_<id>` for an empty one; a
    ///    label already used by another run gets `_RUN<id>` appended).
    /// 2. For each run and transition group the best peak group is kept:
    ///    lowest `m_score`, then highest `intensity`, then lowest
    ///    `feature_id`.
    /// 3. At precursor level the identifier columns are
    ///    `transition_group_id`, `Sequence`, `FullPeptideName`, `Charge` and
    ///    `ProteinName` and a cell is the intensity of the best peak group.
    ///    Higher levels first summarize each peptide (`Sequence`,
    ///    `FullPeptideName`) and run: with `consistent_top` over the `top_n`
    ///    precursors of the peptide with the highest median intensity over the
    ///    runs, the same ones in every run, otherwise over the `top_n` most
    ///    intense peak groups of that run; the summary is the mean. Peptide
    ///    level stops there. Protein and gene level then take, per protein or
    ///    gene (rows with an empty name are skipped), with `consistent_top`
    ///    the `top_n` peptides of highest median and otherwise all, and the
    ///    mean per sample of their summaries.
    /// 4. The normalization of the configuration is applied to the sample
    ///    columns.
    ///
    /// Rows come out in the order of their identifier tuples (byte-wise
    /// string comparison). A sample without a value for a row has `None`.
    /// Where the source sorts with `std::sort` on a key with ties (peak
    /// groups equal in all four keys, precursors of equal median) the
    /// permutation of the Release build's `std::sort` is reproduced.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] (source `Exception::Precondition`) for no rows,
    /// for `top_n` of zero above precursor level, when no row survives the
    /// peptide summary, and for protein or gene level when no peptide has a
    /// protein or gene annotation.
    pub fn build_matrix(
        rows: &[OpenSwathExportRow],
        config: &OpenSwathMatrixExportConfig,
    ) -> Result<OpenSwathQuantMatrix> {
        if rows.is_empty() {
            return Err(precondition(
                "No OpenSWATH rows were available for matrix export.",
            ));
        }
        if config.level != OpenSwathMatrixLevel::Precursor && config.top_n == 0 {
            return Err(precondition(
                "Matrix export requires top_n > 0 for peptide/protein/gene levels.",
            ));
        }

        let sample_names = build_sample_names(rows);
        let best_rows = select_best_peak_groups(rows)?;
        let mut matrix = OpenSwathQuantMatrix::default();

        let emit = |matrix: &mut OpenSwathQuantMatrix,
                    grouped: BTreeMap<Vec<String>, BTreeMap<i64, f64>>|
         -> Result<()> {
            matrix.identifier_rows.clear();
            matrix.values.clear();
            matrix.sample_column_names = sample_names.values().cloned().collect();
            for (identifiers, per_run) in grouped {
                matrix.identifier_rows.push(identifiers);
                let values: Vec<Option<f64>> = sample_names
                    .keys()
                    .map(|run_id| per_run.get(run_id).copied())
                    .collect();
                matrix.values.push(values);
            }
            apply_normalization(matrix, config.normalization)
        };

        if config.level == OpenSwathMatrixLevel::Precursor {
            matrix.identifier_column_names = [
                "transition_group_id",
                "Sequence",
                "FullPeptideName",
                "Charge",
                "ProteinName",
            ]
            .map(String::from)
            .to_vec();
            let mut grouped: BTreeMap<Vec<String>, BTreeMap<i64, f64>> = BTreeMap::new();
            for row in &best_rows {
                grouped
                    .entry(vec![
                        row.transition_group_id.clone(),
                        row.sequence.clone(),
                        row.full_peptide_name.clone(),
                        row.charge.to_string(),
                        row.protein_name.clone(),
                    ])
                    .or_default()
                    .insert(row.run_id, row.intensity);
            }
            emit(&mut matrix, grouped)?;
            return Ok(matrix);
        }

        let peptide_rows = summarize_peptides(&best_rows, config.top_n, config.consistent_top)?;
        if peptide_rows.is_empty() {
            return Err(precondition(
                "No OpenSWATH rows passed peptide-level matrix summarization.",
            ));
        }

        if config.level == OpenSwathMatrixLevel::Peptide {
            matrix.identifier_column_names =
                ["Sequence", "FullPeptideName"].map(String::from).to_vec();
            let mut grouped: BTreeMap<Vec<String>, BTreeMap<i64, f64>> = BTreeMap::new();
            for row in &peptide_rows {
                grouped.insert(
                    vec![row.sequence.clone(), row.full_peptide_name.clone()],
                    row.sample_values.clone(),
                );
            }
            emit(&mut matrix, grouped)?;
            return Ok(matrix);
        }

        let gene_level = config.level != OpenSwathMatrixLevel::Protein;
        matrix.identifier_column_names = vec![
            (if gene_level {
                "GeneName"
            } else {
                "ProteinName"
            })
            .to_owned(),
        ];
        let mut grouped_peptides: BTreeMap<&str, Vec<&PeptideSummaryRow>> = BTreeMap::new();
        for row in &peptide_rows {
            let key = if gene_level {
                row.gene_name.as_str()
            } else {
                row.protein_name.as_str()
            };
            if key.is_empty() {
                continue;
            }
            grouped_peptides.entry(key).or_default().push(row);
        }
        if grouped_peptides.is_empty() {
            return Err(precondition(if gene_level {
                "No gene annotations were available for matrix export."
            } else {
                "No protein annotations were available for matrix export."
            }));
        }

        let mut grouped: BTreeMap<Vec<String>, BTreeMap<i64, f64>> = BTreeMap::new();
        for (key, peptides) in grouped_peptides {
            let mut selected: Vec<&PeptideSummaryRow> = peptides.clone();
            if config.consistent_top {
                let mut medians: Vec<(&PeptideSummaryRow, f64)> =
                    Vec::with_capacity(peptides.len());
                for peptide in &peptides {
                    let values: Vec<f64> = peptide.sample_values.values().copied().collect();
                    medians.push((peptide, median(values)?));
                }
                source_sort_by(&mut medians, |lhs, rhs| lhs.1 > rhs.1)?;
                selected = medians
                    .iter()
                    .take(config.top_n.min(medians.len()))
                    .map(|(peptide, _)| *peptide)
                    .collect();
            }
            let mut per_sample: BTreeMap<i64, Vec<f64>> = BTreeMap::new();
            for peptide in selected {
                for (sample, value) in &peptide.sample_values {
                    per_sample.entry(*sample).or_default().push(*value);
                }
            }
            let sample_means: BTreeMap<i64, f64> = per_sample
                .into_iter()
                .map(|(sample, values)| (sample, mean(&values)))
                .collect();
            grouped.insert(vec![key.to_owned()], sample_means);
        }
        emit(&mut matrix, grouped)?;
        Ok(matrix)
    }

    /// The typed columns of the Parquet output of `matrix`: one non-null
    /// `utf8` column per identifier column, then one `float64` column per
    /// sample with a null for a missing value.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the matrix is empty or malformed, as for
    /// [`write_matrix`](Self::write_matrix).
    pub fn table(matrix: &OpenSwathQuantMatrix) -> Result<ExportTable> {
        Self::check_writable(matrix)?;
        let mut columns: Vec<(String, ExportColumn)> = Vec::new();
        for (col, name) in matrix.identifier_column_names.iter().enumerate() {
            let values: Vec<Option<String>> = matrix
                .identifier_rows
                .iter()
                .map(|row| Some(row[col].clone()))
                .collect();
            columns.push((name.clone(), ExportColumn::Utf8(values)));
        }
        for (col, name) in matrix.sample_column_names.iter().enumerate() {
            let values: Vec<Option<f64>> = matrix.values.iter().map(|row| row[col]).collect();
            columns.push((name.clone(), ExportColumn::Float64(values)));
        }
        Ok(ExportTable { columns })
    }

    fn check_writable(matrix: &OpenSwathQuantMatrix) -> Result<()> {
        if matrix.identifier_rows.is_empty() {
            return Err(precondition(
                "Cannot write an empty OpenSWATH quantification matrix.",
            ));
        }
        validate_matrix_dimensions(matrix)
    }

    /// Write a previously built matrix (source `writeMatrix`).
    ///
    /// The TSV has the identifier column names followed by the sample column
    /// names as header, then one line per row: the identifiers, then each value
    /// at `%g` precision 6 (`1000`, `1.23457e+06`) or the empty field for a
    /// missing value. Strings are written verbatim.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] (source `Exception::Precondition`) for a matrix
    /// without rows, with different identifier and value row counts, or a
    /// row whose width differs from its column names; nothing is written
    /// then. [`Error::Io`] when the file cannot be created or written.
    /// [`Error::Unsupported`] for Parquet in a build without the `arrow`
    /// feature; with it the file is [`table`](Self::table) written by
    /// `ExportTable::write_parquet`.
    pub fn write_matrix(
        filename: impl AsRef<Path>,
        matrix: &OpenSwathQuantMatrix,
        config: &OpenSwathMatrixExportConfig,
    ) -> Result<()> {
        Self::check_writable(matrix)?;
        match config.format {
            OpenSwathExportFileFormat::Tsv => {
                let mut file = File::create(filename.as_ref())?;
                Self::write_tsv(&mut file, matrix)?;
                file.flush()?;
                Ok(())
            }
            OpenSwathExportFileFormat::Parquet => {
                #[cfg(feature = "arrow")]
                {
                    Self::table(matrix)?.write_parquet(filename)
                }
                #[cfg(not(feature = "arrow"))]
                {
                    let _ = filename;
                    Err(Error::Unsupported(
                        "Parquet output of an OpenSWATH matrix needs the `arrow` feature".into(),
                    ))
                }
            }
        }
    }

    /// The TSV of `matrix` written to any writer; see
    /// [`write_matrix`](Self::write_matrix) for the layout. The matrix is
    /// *not* checked here.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the writer fails; [`Error::InvalidValue`] for a row
    /// shorter than its column names.
    pub fn write_tsv<W: Write>(writer: &mut W, matrix: &OpenSwathQuantMatrix) -> Result<()> {
        let mut out = String::new();
        out.push_str(&matrix.identifier_column_names.join("\t"));
        for sample in &matrix.sample_column_names {
            out.push('\t');
            out.push_str(sample);
        }
        out.push('\n');
        writer.write_all(out.as_bytes())?;
        for (identifiers, values) in matrix.identifier_rows.iter().zip(&matrix.values) {
            out.clear();
            out.push_str(&identifiers.join("\t"));
            for value in values {
                out.push('\t');
                if let Some(value) = value {
                    out.push_str(&ostream_g(*value, 6));
                }
            }
            out.push('\n');
            writer.write_all(out.as_bytes())?;
        }
        Ok(())
    }
}
