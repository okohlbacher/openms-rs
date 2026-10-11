// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Writing of filtered OpenSWATH feature results.
//!
//! Ports `ANALYSIS/OPENSWATH/OpenSwathResultsExporter.h` with its `.cpp`.
//! [`OpenSwathResultsExporter::write`] writes the rows of an
//! [`OpenSwathExportRow`] list as tab-separated text. The exporter is
//! agnostic of the input: reading the rows out of an OSW or Parquet file is the
//! business of the readers, as in the source.
//!
//! The Parquet arm of the source goes through `FORMAT/ParquetFile.h`, which is
//! not part of this module. What this module owns of it is the schema:
//! [`OpenSwathResultsExporter::table`] returns the 59 typed columns, their
//! names and their null policy as an [`ExportTable`], and a writer built on
//! the `ParquetFile` port turns that into the file. Until then
//! [`OpenSwathResultsExporter::write`] with
//! [`OpenSwathExportFileFormat::Parquet`] returns
//! [`Error::Unsupported`](crate::Error::Unsupported) and writes nothing.
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_EXPORT_SUPPORT.md).
//!
//! # Example
//!
//! ```
//! use openms::analysis::openswath::openswath_export_data::OpenSwathExportRow;
//! use openms::analysis::openswath::openswath_results_exporter::OpenSwathResultsExporter;
//!
//! let row = OpenSwathExportRow {
//!     run_id: 99,
//!     feature_id: 7777,
//!     transition_group_id: "TG1".into(),
//!     sequence: "PEPTIDEK".into(),
//!     rt: 120.5,
//!     ..OpenSwathExportRow::default()
//! };
//! let text = OpenSwathResultsExporter::tsv_string(&[row]).unwrap();
//! let mut lines = text.lines();
//! assert!(lines.next().unwrap().starts_with("run_id\tfilename\trun_name\tfeature_id"));
//! let fields: Vec<&str> = lines.next().unwrap().split('\t').collect();
//! assert_eq!(fields[3], "7777");
//! assert_eq!(fields[14], "120.5");
//! ```

use super::openswath_export_config::OpenSwathResultsExportConfig;
use super::openswath_export_data::{OpenSwathExportFileFormat, OpenSwathExportRow};
use crate::format::file_info::text_format::{ostream_g, to_str};
use crate::{Error, Result};
use std::fmt::Write as _;
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// The header line of the TSV output (without the newline), the 59 column
/// names of the source in order.
pub const RESULTS_TSV_COLUMNS: [&str; 59] = [
    "run_id",
    "filename",
    "run_name",
    "feature_id",
    "peptide_id",
    "transition_group_id",
    "precursor_id",
    "decoy",
    "Sequence",
    "FullPeptideName",
    "ProteinName",
    "GeneName",
    "Charge",
    "mz",
    "RT",
    "assay_rt",
    "delta_rt",
    "iRT",
    "assay_iRT",
    "delta_iRT",
    "Intensity",
    "aggr_prec_Peak_Area",
    "aggr_prec_Peak_Apex",
    "leftWidth",
    "rightWidth",
    "EXP_IM",
    "IM_leftWidth",
    "IM_rightWidth",
    "ms1_pep",
    "ms2_pep",
    "precursor_pep",
    "ipf_pep",
    "peak_group_rank",
    "d_score",
    "m_score",
    "pep",
    "ms2_m_score",
    "m_score_peptide_run_specific",
    "m_score_peptide_experiment_wide",
    "m_score_peptide_global",
    "m_score_protein_run_specific",
    "m_score_protein_experiment_wide",
    "m_score_protein_global",
    "m_score_gene_run_specific",
    "m_score_gene_experiment_wide",
    "m_score_gene_global",
    "alignment_group_id",
    "alignment_reference_feature_id",
    "alignment_reference_rt",
    "alignment_pep",
    "alignment_qvalue",
    "from_alignment",
    "aggr_Peak_Area",
    "aggr_Peak_Apex",
    "aggr_Fragment_Annotation",
    "ipf_FullUniModPeptideName",
    "ipf_precursor_peakgroup_pep",
    "ipf_peptidoform_pep",
    "ipf_peptidoform_m_score",
];

/// One typed column of an [`ExportTable`]; the Arrow type of the source's
/// column in parentheses.
#[derive(Clone, Debug, PartialEq)]
pub enum ExportColumn {
    /// 64-bit integers with nulls (`int64`).
    Int64(Vec<Option<i64>>),
    /// Doubles with nulls (`float64`).
    Float64(Vec<Option<f64>>),
    /// Strings with nulls (`utf8`).
    Utf8(Vec<Option<String>>),
    /// Booleans, never null (`boolean`).
    Boolean(Vec<bool>),
}

impl ExportColumn {
    /// The number of rows of the column.
    pub fn len(&self) -> usize {
        match self {
            Self::Int64(v) => v.len(),
            Self::Float64(v) => v.len(),
            Self::Utf8(v) => v.len(),
            Self::Boolean(v) => v.len(),
        }
    }

    /// Whether the column has no rows.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// A table of named, typed columns of equal length: the native stand-in for
/// the `arrow::Table` the source builds, which a Parquet writer consumes.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ExportTable {
    /// Column names and data, in output order.
    pub columns: Vec<(String, ExportColumn)>,
}

impl ExportTable {
    /// The number of rows (zero without columns).
    pub fn num_rows(&self) -> usize {
        self.columns.first().map_or(0, |(_, c)| c.len())
    }

    /// The number of columns.
    pub fn num_columns(&self) -> usize {
        self.columns.len()
    }

    /// The column called `name`.
    pub fn column(&self, name: &str) -> Option<&ExportColumn> {
        self.columns.iter().find(|(n, _)| n == name).map(|(_, c)| c)
    }
}

/// The writer of OpenSWATH feature results (source
/// `OpenSwathResultsExporter`, a class of one static member).
#[derive(Clone, Copy, Debug, Default)]
pub struct OpenSwathResultsExporter;

fn opt_f64(value: Option<f64>) -> String {
    value.map(to_str).unwrap_or_default()
}

fn opt_i64(value: Option<i64>) -> String {
    value.map(|v| v.to_string()).unwrap_or_default()
}

fn flag(value: bool) -> char {
    if value { '1' } else { '0' }
}

/// `os << double` with the default stream state: `%g` at precision 6.
fn stream_f64(value: f64) -> String {
    ostream_g(value, 6)
}

impl OpenSwathResultsExporter {
    /// Write `rows` to `filename` in the format of `config`.
    ///
    /// Empty exports are allowed: they produce a header-only TSV. The source
    /// truncates an existing file and does not remove a partial one after a
    /// failure; so does this.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be created or written (source
    /// `Exception::FileNotWritable`), and [`Error::Unsupported`] for
    /// [`OpenSwathExportFileFormat::Parquet`]: the Parquet writer needs the
    /// `FORMAT/ParquetFile.h` port. Nothing is created then.
    pub fn write(
        filename: impl AsRef<Path>,
        rows: &[OpenSwathExportRow],
        config: &OpenSwathResultsExportConfig,
    ) -> Result<()> {
        match config.format {
            OpenSwathExportFileFormat::Tsv => {
                let mut file = File::create(filename.as_ref())?;
                Self::write_tsv(&mut file, rows)?;
                file.flush()?;
                Ok(())
            }
            OpenSwathExportFileFormat::Parquet => Err(Error::Unsupported(
                "Parquet output of the OpenSWATH results needs the FORMAT/ParquetFile.h port; \
                 OpenSwathResultsExporter::table builds the typed columns it will write"
                    .into(),
            )),
        }
    }

    /// The TSV text of `rows` (header line included).
    ///
    /// # Errors
    ///
    /// Never fails in practice; the `Result` keeps the signature of
    /// [`write_tsv`](Self::write_tsv).
    pub fn tsv_string(rows: &[OpenSwathExportRow]) -> Result<String> {
        let mut bytes = Vec::new();
        Self::write_tsv(&mut bytes, rows)?;
        String::from_utf8(bytes).map_err(|e| Error::InvalidValue(e.to_string()))
    }

    /// Write `rows` as TSV to any writer: the header line, then one line per
    /// row.
    ///
    /// Numbers follow the source stream state: doubles at six significant
    /// digits (`%g`, so `120.5` and `1.23457e+06`), optional doubles through
    /// `StringUtils::toStr` at full precision (`321.5`, `1.0e-05`), optional
    /// integers in decimal, flags as `1` and `0`, an absent optional as the
    /// empty field. Strings are written verbatim: a tab or newline inside one
    /// breaks the row, exactly as in the source.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the writer fails.
    pub fn write_tsv<W: Write>(writer: &mut W, rows: &[OpenSwathExportRow]) -> Result<()> {
        let mut out = String::with_capacity(512);
        out.push_str(&RESULTS_TSV_COLUMNS.join("\t"));
        out.push('\n');
        writer.write_all(out.as_bytes())?;
        for row in rows {
            out.clear();
            let _ = write!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
                row.run_id,
                row.filename,
                row.run_name,
                row.feature_id,
                row.peptide_id,
                row.transition_group_id,
                row.precursor_id,
                flag(row.decoy),
                row.sequence,
                row.full_peptide_name,
                row.protein_name,
                row.gene_name,
                row.charge
            );
            let fields: [String; 7] = [
                stream_f64(row.mz),
                stream_f64(row.rt),
                stream_f64(row.assay_rt),
                stream_f64(row.delta_rt),
                stream_f64(row.irt),
                stream_f64(row.assay_irt),
                stream_f64(row.delta_irt),
            ];
            for field in &fields {
                out.push_str(field);
                out.push('\t');
            }
            out.push_str(&stream_f64(row.intensity));
            out.push('\t');
            out.push_str(&opt_f64(row.aggr_prec_peak_area));
            out.push('\t');
            out.push_str(&opt_f64(row.aggr_prec_peak_apex));
            out.push('\t');
            out.push_str(&stream_f64(row.left_width));
            out.push('\t');
            out.push_str(&stream_f64(row.right_width));
            for value in [
                row.exp_im,
                row.im_left_width,
                row.im_right_width,
                row.ms1_pep,
                row.ms2_pep,
                row.precursor_pep,
                row.ipf_pep,
            ] {
                out.push('\t');
                out.push_str(&opt_f64(value));
            }
            let _ = write!(out, "\t{}\t", row.peak_group_rank);
            out.push_str(&stream_f64(row.d_score));
            out.push('\t');
            out.push_str(&stream_f64(row.m_score));
            for value in [
                row.pep,
                row.ms2_m_score,
                row.peptide_run_specific_qvalue,
                row.peptide_experiment_wide_qvalue,
                row.peptide_global_qvalue,
                row.protein_run_specific_qvalue,
                row.protein_experiment_wide_qvalue,
                row.protein_global_qvalue,
                row.gene_run_specific_qvalue,
                row.gene_experiment_wide_qvalue,
                row.gene_global_qvalue,
            ] {
                out.push('\t');
                out.push_str(&opt_f64(value));
            }
            out.push('\t');
            out.push_str(&opt_i64(row.alignment_group_id));
            out.push('\t');
            out.push_str(&opt_i64(row.alignment_reference_feature_id));
            for value in [
                row.alignment_reference_rt,
                row.alignment_pep,
                row.alignment_qvalue,
            ] {
                out.push('\t');
                out.push_str(&opt_f64(value));
            }
            out.push('\t');
            out.push(flag(row.from_alignment));
            out.push('\t');
            out.push_str(&row.aggr_peak_area);
            out.push('\t');
            out.push_str(&row.aggr_peak_apex);
            out.push('\t');
            out.push_str(&row.aggr_fragment_annotation);
            out.push('\t');
            out.push_str(&row.ipf_full_peptide_name);
            for value in [
                row.ipf_precursor_peakgroup_pep,
                row.ipf_peptidoform_pep,
                row.ipf_peptidoform_m_score,
            ] {
                out.push('\t');
                out.push_str(&opt_f64(value));
            }
            out.push('\n');
            writer.write_all(out.as_bytes())?;
        }
        Ok(())
    }

    /// The typed columns of the Parquet output for `rows`, in the order and
    /// with the names of the source schema (59 columns, also for no rows).
    ///
    /// The types are `int64` for `run_id`, `feature_id`, `peptide_id`,
    /// `precursor_id`, `Charge`, `peak_group_rank` and the two alignment
    /// identifiers, `boolean` for `decoy` and `from_alignment`, `utf8` for the
    /// names, the group identifier and the aggregated text columns, `float64`
    /// for everything else. As in the source an empty string is a null, an
    /// absent optional number is a null, and the two flags and the plain
    /// numbers are never null.
    pub fn table(rows: &[OpenSwathExportRow]) -> ExportTable {
        fn int(
            rows: &[OpenSwathExportRow],
            f: impl Fn(&OpenSwathExportRow) -> i64,
        ) -> ExportColumn {
            ExportColumn::Int64(rows.iter().map(|r| Some(f(r))).collect())
        }
        fn opt_int(
            rows: &[OpenSwathExportRow],
            f: impl Fn(&OpenSwathExportRow) -> Option<i64>,
        ) -> ExportColumn {
            ExportColumn::Int64(rows.iter().map(f).collect())
        }
        fn dbl(
            rows: &[OpenSwathExportRow],
            f: impl Fn(&OpenSwathExportRow) -> f64,
        ) -> ExportColumn {
            ExportColumn::Float64(rows.iter().map(|r| Some(f(r))).collect())
        }
        fn opt_dbl(
            rows: &[OpenSwathExportRow],
            f: impl Fn(&OpenSwathExportRow) -> Option<f64>,
        ) -> ExportColumn {
            ExportColumn::Float64(rows.iter().map(f).collect())
        }
        fn text(
            rows: &[OpenSwathExportRow],
            f: impl Fn(&OpenSwathExportRow) -> &str,
        ) -> ExportColumn {
            ExportColumn::Utf8(
                rows.iter()
                    .map(|r| {
                        let s = f(r);
                        if s.is_empty() {
                            None
                        } else {
                            Some(s.to_owned())
                        }
                    })
                    .collect(),
            )
        }
        fn boolean(
            rows: &[OpenSwathExportRow],
            f: impl Fn(&OpenSwathExportRow) -> bool,
        ) -> ExportColumn {
            ExportColumn::Boolean(rows.iter().map(f).collect())
        }

        let columns: Vec<ExportColumn> = vec![
            int(rows, |r| r.run_id),
            text(rows, |r| &r.filename),
            text(rows, |r| &r.run_name),
            int(rows, |r| r.feature_id),
            int(rows, |r| r.peptide_id),
            text(rows, |r| &r.transition_group_id),
            int(rows, |r| r.precursor_id),
            boolean(rows, |r| r.decoy),
            text(rows, |r| &r.sequence),
            text(rows, |r| &r.full_peptide_name),
            text(rows, |r| &r.protein_name),
            text(rows, |r| &r.gene_name),
            int(rows, |r| i64::from(r.charge)),
            dbl(rows, |r| r.mz),
            dbl(rows, |r| r.rt),
            dbl(rows, |r| r.assay_rt),
            dbl(rows, |r| r.delta_rt),
            dbl(rows, |r| r.irt),
            dbl(rows, |r| r.assay_irt),
            dbl(rows, |r| r.delta_irt),
            dbl(rows, |r| r.intensity),
            opt_dbl(rows, |r| r.aggr_prec_peak_area),
            opt_dbl(rows, |r| r.aggr_prec_peak_apex),
            dbl(rows, |r| r.left_width),
            dbl(rows, |r| r.right_width),
            opt_dbl(rows, |r| r.exp_im),
            opt_dbl(rows, |r| r.im_left_width),
            opt_dbl(rows, |r| r.im_right_width),
            opt_dbl(rows, |r| r.ms1_pep),
            opt_dbl(rows, |r| r.ms2_pep),
            opt_dbl(rows, |r| r.precursor_pep),
            opt_dbl(rows, |r| r.ipf_pep),
            int(rows, |r| i64::from(r.peak_group_rank)),
            dbl(rows, |r| r.d_score),
            dbl(rows, |r| r.m_score),
            opt_dbl(rows, |r| r.pep),
            opt_dbl(rows, |r| r.ms2_m_score),
            opt_dbl(rows, |r| r.peptide_run_specific_qvalue),
            opt_dbl(rows, |r| r.peptide_experiment_wide_qvalue),
            opt_dbl(rows, |r| r.peptide_global_qvalue),
            opt_dbl(rows, |r| r.protein_run_specific_qvalue),
            opt_dbl(rows, |r| r.protein_experiment_wide_qvalue),
            opt_dbl(rows, |r| r.protein_global_qvalue),
            opt_dbl(rows, |r| r.gene_run_specific_qvalue),
            opt_dbl(rows, |r| r.gene_experiment_wide_qvalue),
            opt_dbl(rows, |r| r.gene_global_qvalue),
            opt_int(rows, |r| r.alignment_group_id),
            opt_int(rows, |r| r.alignment_reference_feature_id),
            opt_dbl(rows, |r| r.alignment_reference_rt),
            opt_dbl(rows, |r| r.alignment_pep),
            opt_dbl(rows, |r| r.alignment_qvalue),
            boolean(rows, |r| r.from_alignment),
            text(rows, |r| &r.aggr_peak_area),
            text(rows, |r| &r.aggr_peak_apex),
            text(rows, |r| &r.aggr_fragment_annotation),
            text(rows, |r| &r.ipf_full_peptide_name),
            opt_dbl(rows, |r| r.ipf_precursor_peakgroup_pep),
            opt_dbl(rows, |r| r.ipf_peptidoform_pep),
            opt_dbl(rows, |r| r.ipf_peptidoform_m_score),
        ];
        ExportTable {
            columns: RESULTS_TSV_COLUMNS
                .iter()
                .zip(columns)
                .map(|(name, column)| ((*name).to_owned(), column))
                .collect(),
        }
    }
}
