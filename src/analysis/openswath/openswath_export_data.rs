// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Row, table and matrix types of the OpenSWATH exports.
//!
//! Ports `ANALYSIS/OPENSWATH/OpenSwathExportData.h`, a header of plain
//! aggregates and four enumerations:
//! [`OpenSwathExportFileFormat`], [`OpenSwathIPFExportMode`],
//! [`OpenSwathMatrixLevel`] and [`OpenSwathMatrixNormalization`] with the
//! lower-case names the source's `toString` overloads return, the filtered
//! feature row [`OpenSwathExportRow`] the TSV and matrix exports read, the
//! scored feature and transition rows of the Parquet export
//! ([`OpenSwathFeatureScoreRow`], [`OpenSwathFeatureScoreTable`],
//! [`OpenSwathTransitionScoreRow`], [`OpenSwathTransitionScoreTable`]) and the
//! dense [`OpenSwathQuantMatrix`].
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_EXPORT_SUPPORT.md).
//!
//! # Example
//!
//! ```
//! use openms::analysis::openswath::openswath_export_config::OpenSwathMatrixExportConfig;
//! use openms::analysis::openswath::openswath_export_data::{
//!     OpenSwathExportRow, OpenSwathMatrixLevel,
//! };
//!
//! let row = OpenSwathExportRow::default();
//! assert_eq!(row.run_id, -1);
//! assert_eq!(row.charge, 0);
//! // The enumeration itself value-initialises to its first enumerator; the
//! // `peptide` default belongs to the export configuration
//! // (`OpenSwathExportConfig.h:66`).
//! assert_eq!(OpenSwathMatrixLevel::default().as_str(), "precursor");
//! assert_eq!(
//!     OpenSwathMatrixExportConfig::default().level,
//!     OpenSwathMatrixLevel::Peptide
//! );
//! ```
//!
//! [`OpenSwathExportFileFormat`]: crate::analysis::openswath::openswath_export_data::OpenSwathExportFileFormat
//! [`OpenSwathIPFExportMode`]: crate::analysis::openswath::openswath_export_data::OpenSwathIPFExportMode
//! [`OpenSwathMatrixLevel`]: crate::analysis::openswath::openswath_export_data::OpenSwathMatrixLevel
//! [`OpenSwathMatrixNormalization`]: crate::analysis::openswath::openswath_export_data::OpenSwathMatrixNormalization
//! [`OpenSwathExportRow`]: crate::analysis::openswath::openswath_export_data::OpenSwathExportRow
//! [`OpenSwathFeatureScoreRow`]: crate::analysis::openswath::openswath_export_data::OpenSwathFeatureScoreRow
//! [`OpenSwathFeatureScoreTable`]: crate::analysis::openswath::openswath_export_data::OpenSwathFeatureScoreTable
//! [`OpenSwathTransitionScoreRow`]: crate::analysis::openswath::openswath_export_data::OpenSwathTransitionScoreRow
//! [`OpenSwathTransitionScoreTable`]: crate::analysis::openswath::openswath_export_data::OpenSwathTransitionScoreTable
//! [`OpenSwathQuantMatrix`]: crate::analysis::openswath::openswath_export_data::OpenSwathQuantMatrix

use std::fmt;

/// Output file format of an OpenSWATH export table (source
/// `OpenSwathExportFileFormat`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum OpenSwathExportFileFormat {
    /// Tab-separated text.
    #[default]
    Tsv,
    /// Apache Parquet. Needs the `FORMAT/ParquetFile.h` port; see
    /// [`OpenSwathResultsExporter`](crate::analysis::openswath::openswath_results_exporter::OpenSwathResultsExporter).
    Parquet,
}

impl OpenSwathExportFileFormat {
    /// The source `toString(OpenSwathExportFileFormat)`: `"tsv"` or
    /// `"parquet"`.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Tsv => "tsv",
            Self::Parquet => "parquet",
        }
    }
}

impl fmt::Display for OpenSwathExportFileFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// IPF handling mode of an OpenSWATH export (source `OpenSwathIPFExportMode`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum OpenSwathIPFExportMode {
    /// Report results on peptidoform level.
    #[default]
    Peptidoform,
    /// Augment the OpenSWATH results with the IPF scores.
    Augmented,
    /// Ignore the IPF scores and report the standard OpenSWATH results only.
    Disable,
}

impl OpenSwathIPFExportMode {
    /// The source `toString(OpenSwathIPFExportMode)`: `"peptidoform"`,
    /// `"augmented"` or `"disable"`.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Peptidoform => "peptidoform",
            Self::Augmented => "augmented",
            Self::Disable => "disable",
        }
    }
}

impl fmt::Display for OpenSwathIPFExportMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Quantification matrix level of an OpenSWATH matrix export (source
/// `OpenSwathMatrixLevel`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum OpenSwathMatrixLevel {
    /// One row per precursor (transition group). The type's `Default`, as a
    /// value-initialised C++ `enum class` is its first enumerator; the
    /// *configuration's* default level is [`OpenSwathMatrixLevel::Peptide`].
    #[default]
    Precursor,
    /// One row per peptide (modified sequence).
    Peptide,
    /// One row per protein.
    Protein,
    /// One row per gene.
    Gene,
}

impl OpenSwathMatrixLevel {
    /// The source `toString(OpenSwathMatrixLevel)`: `"precursor"`,
    /// `"peptide"`, `"protein"` or `"gene"`.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Precursor => "precursor",
            Self::Peptide => "peptide",
            Self::Protein => "protein",
            Self::Gene => "gene",
        }
    }
}

impl fmt::Display for OpenSwathMatrixLevel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Normalization of an OpenSWATH matrix export (source
/// `OpenSwathMatrixNormalization`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum OpenSwathMatrixNormalization {
    /// No normalization.
    #[default]
    None,
    /// Divide every sample column by its median.
    Median,
    /// Divide every sample column by its median, then multiply by the median
    /// of the column medians.
    MedianMedian,
}

impl OpenSwathMatrixNormalization {
    /// The source `toString(OpenSwathMatrixNormalization)`: `"none"`,
    /// `"median"` or `"medianmedian"`.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Median => "median",
            Self::MedianMedian => "medianmedian",
        }
    }
}

impl fmt::Display for OpenSwathMatrixNormalization {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// A filtered feature-level row, as the results and matrix exports read it
/// (source `OpenSwathExportRow`).
///
/// The row stores only the columns the user-facing exports need; the wider
/// Parquet scoring export uses [`OpenSwathFeatureScoreRow`]. A source
/// `std::optional<T>` is an [`Option<T>`]; the id fields default to -1 as in
/// the source, everything else to zero, `false` or the empty string.
#[derive(Clone, Debug, PartialEq)]
pub struct OpenSwathExportRow {
    /// Identifier of the run (`RUN.ID`).
    pub run_id: i64,
    /// File name of the run.
    pub filename: String,
    /// Name of the run, the sample column label of the matrix exports.
    pub run_name: String,
    /// Identifier of the feature (`FEATURE.ID`).
    pub feature_id: i64,
    /// Identifier of the peptide.
    pub peptide_id: i64,
    /// Identifier of the precursor.
    pub precursor_id: i64,
    /// Identifier of the transition group.
    pub transition_group_id: String,
    /// Whether the row is a decoy.
    pub decoy: bool,
    /// Unmodified sequence.
    pub sequence: String,
    /// Full (modified) peptide name.
    pub full_peptide_name: String,
    /// Protein name(s).
    pub protein_name: String,
    /// Gene name(s).
    pub gene_name: String,
    /// Precursor charge.
    pub charge: i32,
    /// Precursor m/z.
    pub mz: f64,
    /// Experimental retention time.
    pub rt: f64,
    /// Library (assay) retention time.
    pub assay_rt: f64,
    /// Difference of the experimental and the assay retention time.
    pub delta_rt: f64,
    /// Experimental normalized retention time (iRT).
    pub irt: f64,
    /// Library normalized retention time (assay iRT).
    pub assay_irt: f64,
    /// Difference of the experimental and the assay iRT.
    pub delta_irt: f64,
    /// Intensity (area) of the peak group; the quantity of the matrices.
    pub intensity: f64,
    /// Aggregated precursor peak area.
    pub aggr_prec_peak_area: Option<f64>,
    /// Aggregated precursor peak apex.
    pub aggr_prec_peak_apex: Option<f64>,
    /// Left peak boundary (retention time).
    pub left_width: f64,
    /// Right peak boundary (retention time).
    pub right_width: f64,
    /// Experimental ion mobility.
    pub exp_im: Option<f64>,
    /// Left ion mobility boundary.
    pub im_left_width: Option<f64>,
    /// Right ion mobility boundary.
    pub im_right_width: Option<f64>,
    /// MS1 posterior error probability.
    pub ms1_pep: Option<f64>,
    /// MS2 posterior error probability.
    pub ms2_pep: Option<f64>,
    /// Precursor posterior error probability.
    pub precursor_pep: Option<f64>,
    /// IPF posterior error probability.
    pub ipf_pep: Option<f64>,
    /// Rank of the peak group within its precursor.
    pub peak_group_rank: i32,
    /// Discriminant score.
    pub d_score: f64,
    /// Run-specific peak group q-value; the criterion of the best peak group
    /// selection of the matrices.
    pub m_score: f64,
    /// Posterior error probability.
    pub pep: Option<f64>,
    /// MS2 q-value.
    pub ms2_m_score: Option<f64>,
    /// Run-specific peptide q-value.
    pub peptide_run_specific_qvalue: Option<f64>,
    /// Experiment-wide peptide q-value.
    pub peptide_experiment_wide_qvalue: Option<f64>,
    /// Global peptide q-value.
    pub peptide_global_qvalue: Option<f64>,
    /// Run-specific protein q-value.
    pub protein_run_specific_qvalue: Option<f64>,
    /// Experiment-wide protein q-value.
    pub protein_experiment_wide_qvalue: Option<f64>,
    /// Global protein q-value.
    pub protein_global_qvalue: Option<f64>,
    /// Run-specific gene q-value.
    pub gene_run_specific_qvalue: Option<f64>,
    /// Experiment-wide gene q-value.
    pub gene_experiment_wide_qvalue: Option<f64>,
    /// Global gene q-value.
    pub gene_global_qvalue: Option<f64>,
    /// Alignment group identifier.
    pub alignment_group_id: Option<i64>,
    /// Feature identifier of the alignment reference.
    pub alignment_reference_feature_id: Option<i64>,
    /// Retention time of the alignment reference.
    pub alignment_reference_rt: Option<f64>,
    /// Posterior error probability of the alignment.
    pub alignment_pep: Option<f64>,
    /// q-value of the alignment.
    pub alignment_qvalue: Option<f64>,
    /// Whether the row was recovered through the alignment.
    pub from_alignment: bool,
    /// Aggregated (multi-transition) peak areas, as text.
    pub aggr_peak_area: String,
    /// Aggregated (multi-transition) peak apexes, as text.
    pub aggr_peak_apex: String,
    /// Aggregated fragment annotations, as text.
    pub aggr_fragment_annotation: String,
    /// IPF full peptide name.
    pub ipf_full_peptide_name: String,
    /// IPF precursor peak group posterior error probability.
    pub ipf_precursor_peakgroup_pep: Option<f64>,
    /// IPF peptidoform posterior error probability.
    pub ipf_peptidoform_pep: Option<f64>,
    /// IPF peptidoform q-value.
    pub ipf_peptidoform_m_score: Option<f64>,
}

impl Default for OpenSwathExportRow {
    fn default() -> Self {
        Self {
            run_id: -1,
            filename: String::new(),
            run_name: String::new(),
            feature_id: -1,
            peptide_id: -1,
            precursor_id: -1,
            transition_group_id: String::new(),
            decoy: false,
            sequence: String::new(),
            full_peptide_name: String::new(),
            protein_name: String::new(),
            gene_name: String::new(),
            charge: 0,
            mz: 0.0,
            rt: 0.0,
            assay_rt: 0.0,
            delta_rt: 0.0,
            irt: 0.0,
            assay_irt: 0.0,
            delta_irt: 0.0,
            intensity: 0.0,
            aggr_prec_peak_area: None,
            aggr_prec_peak_apex: None,
            left_width: 0.0,
            right_width: 0.0,
            exp_im: None,
            im_left_width: None,
            im_right_width: None,
            ms1_pep: None,
            ms2_pep: None,
            precursor_pep: None,
            ipf_pep: None,
            peak_group_rank: 0,
            d_score: 0.0,
            m_score: 0.0,
            pep: None,
            ms2_m_score: None,
            peptide_run_specific_qvalue: None,
            peptide_experiment_wide_qvalue: None,
            peptide_global_qvalue: None,
            protein_run_specific_qvalue: None,
            protein_experiment_wide_qvalue: None,
            protein_global_qvalue: None,
            gene_run_specific_qvalue: None,
            gene_experiment_wide_qvalue: None,
            gene_global_qvalue: None,
            alignment_group_id: None,
            alignment_reference_feature_id: None,
            alignment_reference_rt: None,
            alignment_pep: None,
            alignment_qvalue: None,
            from_alignment: false,
            aggr_peak_area: String::new(),
            aggr_peak_apex: String::new(),
            aggr_fragment_annotation: String::new(),
            ipf_full_peptide_name: String::new(),
            ipf_precursor_peakgroup_pep: None,
            ipf_peptidoform_pep: None,
            ipf_peptidoform_m_score: None,
        }
    }
}

/// A feature-level scored row of the Parquet export (source
/// `OpenSwathFeatureScoreRow`).
#[derive(Clone, Debug, PartialEq)]
pub struct OpenSwathFeatureScoreRow {
    /// Identifier of the protein.
    pub protein_id: i64,
    /// Identifier of the peptide.
    pub peptide_id: i64,
    /// Identifier of the IPF peptidoform.
    pub ipf_peptide_id: Option<i64>,
    /// Identifier of the precursor.
    pub precursor_id: i64,
    /// Protein accession.
    pub protein_accession: String,
    /// Unmodified sequence.
    pub unmodified_sequence: String,
    /// Modified sequence.
    pub modified_sequence: String,
    /// TraML identifier of the precursor.
    pub precursor_traml_id: String,
    /// Group label of the precursor.
    pub precursor_group_label: String,
    /// Precursor m/z.
    pub precursor_mz: f64,
    /// Precursor charge.
    pub precursor_charge: i32,
    /// Library intensity of the precursor.
    pub precursor_library_intensity: Option<f64>,
    /// Library retention time of the precursor.
    pub precursor_library_rt: Option<f64>,
    /// Library drift time of the precursor.
    pub precursor_library_drift_time: Option<f64>,
    /// Identifier of the gene.
    pub gene_id: Option<i64>,
    /// Gene name.
    pub gene_name: String,
    /// Whether the gene is a decoy.
    pub gene_decoy: Option<bool>,
    /// Whether the protein is a decoy.
    pub protein_decoy: bool,
    /// Whether the peptide is a decoy.
    pub peptide_decoy: bool,
    /// Whether the precursor is a decoy.
    pub precursor_decoy: bool,
    /// Identifier of the run.
    pub run_id: i64,
    /// File name of the run.
    pub filename: String,
    /// Identifier of the feature.
    pub feature_id: i64,
    /// Experimental retention time.
    pub exp_rt: f64,
    /// Experimental ion mobility.
    pub exp_im: Option<f64>,
    /// Normalized retention time.
    pub norm_rt: f64,
    /// Retention time difference to the assay.
    pub delta_rt: f64,
    /// Left peak boundary.
    pub left_width: f64,
    /// Right peak boundary.
    pub right_width: f64,
    /// Left ion mobility boundary.
    pub im_left_width: Option<f64>,
    /// Right ion mobility boundary.
    pub im_right_width: Option<f64>,
    /// Values of the dynamically discovered `FEATURE_MS1` columns.
    pub feature_ms1_values: Vec<f64>,
    /// Values of the dynamically discovered `FEATURE_MS2` columns.
    pub feature_ms2_values: Vec<f64>,
    /// MS1 score.
    pub score_ms1_score: Option<f64>,
    /// MS1 rank.
    pub score_ms1_rank: Option<i32>,
    /// MS1 p-value.
    pub score_ms1_pvalue: Option<f64>,
    /// MS1 q-value.
    pub score_ms1_qvalue: Option<f64>,
    /// MS1 posterior error probability.
    pub score_ms1_pep: Option<f64>,
    /// MS2 score.
    pub score_ms2_score: Option<f64>,
    /// MS2 peak group rank.
    pub score_ms2_peak_group_rank: Option<i32>,
    /// MS2 p-value.
    pub score_ms2_pvalue: Option<f64>,
    /// MS2 q-value.
    pub score_ms2_qvalue: Option<f64>,
    /// MS2 posterior error probability.
    pub score_ms2_pep: Option<f64>,
    /// IPF precursor peak group posterior error probability.
    pub score_ipf_precursor_peakgroup_pep: Option<f64>,
    /// IPF posterior error probability.
    pub score_ipf_pep: Option<f64>,
    /// IPF q-value.
    pub score_ipf_qvalue: Option<f64>,
    /// Global peptide score.
    pub score_peptide_global_score: Option<f64>,
    /// Global peptide p-value.
    pub score_peptide_global_pvalue: Option<f64>,
    /// Global peptide q-value.
    pub score_peptide_global_qvalue: Option<f64>,
    /// Global peptide posterior error probability.
    pub score_peptide_global_pep: Option<f64>,
    /// Experiment-wide peptide score.
    pub score_peptide_experiment_wide_score: Option<f64>,
    /// Experiment-wide peptide p-value.
    pub score_peptide_experiment_wide_pvalue: Option<f64>,
    /// Experiment-wide peptide q-value.
    pub score_peptide_experiment_wide_qvalue: Option<f64>,
    /// Experiment-wide peptide posterior error probability.
    pub score_peptide_experiment_wide_pep: Option<f64>,
    /// Run-specific peptide score.
    pub score_peptide_run_specific_score: Option<f64>,
    /// Run-specific peptide p-value.
    pub score_peptide_run_specific_pvalue: Option<f64>,
    /// Run-specific peptide q-value.
    pub score_peptide_run_specific_qvalue: Option<f64>,
    /// Run-specific peptide posterior error probability.
    pub score_peptide_run_specific_pep: Option<f64>,
    /// Global protein score.
    pub score_protein_global_score: Option<f64>,
    /// Global protein p-value.
    pub score_protein_global_pvalue: Option<f64>,
    /// Global protein q-value.
    pub score_protein_global_qvalue: Option<f64>,
    /// Global protein posterior error probability.
    pub score_protein_global_pep: Option<f64>,
    /// Experiment-wide protein score.
    pub score_protein_experiment_wide_score: Option<f64>,
    /// Experiment-wide protein p-value.
    pub score_protein_experiment_wide_pvalue: Option<f64>,
    /// Experiment-wide protein q-value.
    pub score_protein_experiment_wide_qvalue: Option<f64>,
    /// Experiment-wide protein posterior error probability.
    pub score_protein_experiment_wide_pep: Option<f64>,
    /// Run-specific protein score.
    pub score_protein_run_specific_score: Option<f64>,
    /// Run-specific protein p-value.
    pub score_protein_run_specific_pvalue: Option<f64>,
    /// Run-specific protein q-value.
    pub score_protein_run_specific_qvalue: Option<f64>,
    /// Run-specific protein posterior error probability.
    pub score_protein_run_specific_pep: Option<f64>,
    /// Global gene score.
    pub score_gene_global_score: Option<f64>,
    /// Global gene p-value.
    pub score_gene_global_pvalue: Option<f64>,
    /// Global gene q-value.
    pub score_gene_global_qvalue: Option<f64>,
    /// Global gene posterior error probability.
    pub score_gene_global_pep: Option<f64>,
    /// Experiment-wide gene score.
    pub score_gene_experiment_wide_score: Option<f64>,
    /// Experiment-wide gene p-value.
    pub score_gene_experiment_wide_pvalue: Option<f64>,
    /// Experiment-wide gene q-value.
    pub score_gene_experiment_wide_qvalue: Option<f64>,
    /// Experiment-wide gene posterior error probability.
    pub score_gene_experiment_wide_pep: Option<f64>,
    /// Run-specific gene score.
    pub score_gene_run_specific_score: Option<f64>,
    /// Run-specific gene p-value.
    pub score_gene_run_specific_pvalue: Option<f64>,
    /// Run-specific gene q-value.
    pub score_gene_run_specific_qvalue: Option<f64>,
    /// Run-specific gene posterior error probability.
    pub score_gene_run_specific_pep: Option<f64>,
    /// Alignment group identifier.
    pub alignment_group_id: Option<i64>,
    /// Feature identifier of the alignment reference.
    pub alignment_reference_feature_id: Option<i64>,
    /// Retention time of the alignment reference.
    pub alignment_reference_rt: Option<f64>,
    /// Posterior error probability of the alignment.
    pub alignment_pep: Option<f64>,
    /// q-value of the alignment.
    pub alignment_qvalue: Option<f64>,
    /// Whether the row was recovered through the alignment.
    pub from_alignment: bool,
}

impl Default for OpenSwathFeatureScoreRow {
    fn default() -> Self {
        Self {
            protein_id: -1,
            peptide_id: -1,
            ipf_peptide_id: None,
            precursor_id: -1,
            protein_accession: String::new(),
            unmodified_sequence: String::new(),
            modified_sequence: String::new(),
            precursor_traml_id: String::new(),
            precursor_group_label: String::new(),
            precursor_mz: 0.0,
            precursor_charge: 0,
            precursor_library_intensity: None,
            precursor_library_rt: None,
            precursor_library_drift_time: None,
            gene_id: None,
            gene_name: String::new(),
            gene_decoy: None,
            protein_decoy: false,
            peptide_decoy: false,
            precursor_decoy: false,
            run_id: -1,
            filename: String::new(),
            feature_id: -1,
            exp_rt: 0.0,
            exp_im: None,
            norm_rt: 0.0,
            delta_rt: 0.0,
            left_width: 0.0,
            right_width: 0.0,
            im_left_width: None,
            im_right_width: None,
            feature_ms1_values: Vec::new(),
            feature_ms2_values: Vec::new(),
            score_ms1_score: None,
            score_ms1_rank: None,
            score_ms1_pvalue: None,
            score_ms1_qvalue: None,
            score_ms1_pep: None,
            score_ms2_score: None,
            score_ms2_peak_group_rank: None,
            score_ms2_pvalue: None,
            score_ms2_qvalue: None,
            score_ms2_pep: None,
            score_ipf_precursor_peakgroup_pep: None,
            score_ipf_pep: None,
            score_ipf_qvalue: None,
            score_peptide_global_score: None,
            score_peptide_global_pvalue: None,
            score_peptide_global_qvalue: None,
            score_peptide_global_pep: None,
            score_peptide_experiment_wide_score: None,
            score_peptide_experiment_wide_pvalue: None,
            score_peptide_experiment_wide_qvalue: None,
            score_peptide_experiment_wide_pep: None,
            score_peptide_run_specific_score: None,
            score_peptide_run_specific_pvalue: None,
            score_peptide_run_specific_qvalue: None,
            score_peptide_run_specific_pep: None,
            score_protein_global_score: None,
            score_protein_global_pvalue: None,
            score_protein_global_qvalue: None,
            score_protein_global_pep: None,
            score_protein_experiment_wide_score: None,
            score_protein_experiment_wide_pvalue: None,
            score_protein_experiment_wide_qvalue: None,
            score_protein_experiment_wide_pep: None,
            score_protein_run_specific_score: None,
            score_protein_run_specific_pvalue: None,
            score_protein_run_specific_qvalue: None,
            score_protein_run_specific_pep: None,
            score_gene_global_score: None,
            score_gene_global_pvalue: None,
            score_gene_global_qvalue: None,
            score_gene_global_pep: None,
            score_gene_experiment_wide_score: None,
            score_gene_experiment_wide_pvalue: None,
            score_gene_experiment_wide_qvalue: None,
            score_gene_experiment_wide_pep: None,
            score_gene_run_specific_score: None,
            score_gene_run_specific_pvalue: None,
            score_gene_run_specific_qvalue: None,
            score_gene_run_specific_pep: None,
            alignment_group_id: None,
            alignment_reference_feature_id: None,
            alignment_reference_rt: None,
            alignment_pep: None,
            alignment_qvalue: None,
            from_alignment: false,
        }
    }
}

/// The feature-level Parquet export table with its dynamically discovered
/// `FEATURE_MS1` / `FEATURE_MS2` columns (source `OpenSwathFeatureScoreTable`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OpenSwathFeatureScoreTable {
    /// Names of the `FEATURE_MS1` columns, in the order of every row's
    /// `feature_ms1_values`.
    pub feature_ms1_column_names: Vec<String>,
    /// Names of the `FEATURE_MS2` columns, in the order of every row's
    /// `feature_ms2_values`.
    pub feature_ms2_column_names: Vec<String>,
    /// The rows.
    pub rows: Vec<OpenSwathFeatureScoreRow>,
}

/// A transition-level scored row of the optional Parquet export (source
/// `OpenSwathTransitionScoreRow`).
#[derive(Clone, Debug, PartialEq)]
pub struct OpenSwathTransitionScoreRow {
    /// Identifier of the run.
    pub run_id: Option<i64>,
    /// Identifier of the IPF peptidoform.
    pub ipf_peptide_id: Option<i64>,
    /// Identifier of the precursor.
    pub precursor_id: i64,
    /// Identifier of the transition.
    pub transition_id: i64,
    /// TraML identifier of the transition.
    pub transition_traml_id: String,
    /// Product m/z.
    pub product_mz: f64,
    /// Charge of the fragment.
    pub transition_charge: i32,
    /// Ion type of the fragment.
    pub transition_type: String,
    /// Ordinal of the fragment ion.
    pub transition_ordinal: i32,
    /// Annotation of the fragment.
    pub annotation: String,
    /// Whether the transition is a detecting one.
    pub transition_detecting: bool,
    /// Library intensity of the transition.
    pub transition_library_intensity: Option<f64>,
    /// Whether the transition is a decoy.
    pub transition_decoy: bool,
    /// Identifier of the feature.
    pub feature_id: Option<i64>,
    /// Values of the dynamically discovered `FEATURE_TRANSITION` columns.
    pub feature_transition_values: Vec<Option<f64>>,
    /// Transition score.
    pub score_transition_score: Option<f64>,
    /// Transition rank.
    pub score_transition_rank: Option<i32>,
    /// Transition p-value.
    pub score_transition_pvalue: Option<f64>,
    /// Transition q-value.
    pub score_transition_qvalue: Option<f64>,
    /// Transition posterior error probability.
    pub score_transition_pep: Option<f64>,
}

impl Default for OpenSwathTransitionScoreRow {
    fn default() -> Self {
        Self {
            run_id: None,
            ipf_peptide_id: None,
            precursor_id: -1,
            transition_id: -1,
            transition_traml_id: String::new(),
            product_mz: 0.0,
            transition_charge: 0,
            transition_type: String::new(),
            transition_ordinal: 0,
            annotation: String::new(),
            transition_detecting: false,
            transition_library_intensity: None,
            transition_decoy: false,
            feature_id: None,
            feature_transition_values: Vec::new(),
            score_transition_score: None,
            score_transition_rank: None,
            score_transition_pvalue: None,
            score_transition_qvalue: None,
            score_transition_pep: None,
        }
    }
}

/// The transition-level Parquet export table with its dynamically discovered
/// `FEATURE_TRANSITION` columns (source `OpenSwathTransitionScoreTable`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OpenSwathTransitionScoreTable {
    /// Names of the `FEATURE_TRANSITION` columns, in the order of every row's
    /// `feature_transition_values`.
    pub feature_transition_column_names: Vec<String>,
    /// The rows.
    pub rows: Vec<OpenSwathTransitionScoreRow>,
}

/// The dense matrix of the matrix exports (source `OpenSwathQuantMatrix`).
///
/// Values are stored by row; a missing value is [`None`], so sparse run
/// coverage survives into the TSV (an empty cell) and Parquet (a null) output.
/// Row widths are *not* enforced by construction, as in the source; the
/// matrix writer checks them.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OpenSwathQuantMatrix {
    /// Names of the identifier columns.
    pub identifier_column_names: Vec<String>,
    /// One identifier tuple per row, as wide as `identifier_column_names`.
    pub identifier_rows: Vec<Vec<String>>,
    /// Names of the sample columns, one per run.
    pub sample_column_names: Vec<String>,
    /// One value row per identifier row, as wide as `sample_column_names`.
    pub values: Vec<Vec<Option<f64>>>,
}

#[cfg(test)]
mod tests {
    //! Expected values are the member initialisers and `toString` bodies of
    //! `OpenSwathExportData.h` at core `bc9cc12` (the header has no class
    //! test of its own; `OpenSwathParquetExporter_test.cpp` uses it).
    use super::*;

    #[test]
    fn to_string_matches_the_header() {
        use OpenSwathExportFileFormat as F;
        use OpenSwathIPFExportMode as I;
        use OpenSwathMatrixLevel as L;
        use OpenSwathMatrixNormalization as N;
        assert_eq!([F::Tsv.as_str(), F::Parquet.as_str()], ["tsv", "parquet"]);
        assert_eq!(
            [I::Peptidoform, I::Augmented, I::Disable].map(|x| x.as_str()),
            ["peptidoform", "augmented", "disable"]
        );
        assert_eq!(
            [L::Precursor, L::Peptide, L::Protein, L::Gene].map(|x| x.as_str()),
            ["precursor", "peptide", "protein", "gene"]
        );
        assert_eq!(
            [N::None, N::Median, N::MedianMedian].map(|x| x.as_str()),
            ["none", "median", "medianmedian"]
        );
        assert_eq!(F::Parquet.to_string(), "parquet");
        assert_eq!(N::MedianMedian.to_string(), "medianmedian");
        // A value-initialised `enum class` is its first enumerator.
        assert_eq!(F::default(), F::Tsv);
        assert_eq!(I::default(), I::Peptidoform);
        assert_eq!(L::default(), L::Precursor);
        assert_eq!(N::default(), N::None);
    }

    #[test]
    fn export_row_defaults_match_the_header() {
        let r = OpenSwathExportRow::default();
        assert_eq!(
            (r.run_id, r.feature_id, r.peptide_id, r.precursor_id),
            (-1, -1, -1, -1)
        );
        assert!(!r.decoy && !r.from_alignment);
        assert_eq!((r.charge, r.peak_group_rank), (0, 0));
        for v in [
            r.mz,
            r.rt,
            r.assay_rt,
            r.delta_rt,
            r.irt,
            r.assay_irt,
            r.delta_irt,
            r.intensity,
            r.left_width,
            r.right_width,
            r.d_score,
            r.m_score,
        ] {
            assert_eq!(v.to_bits(), 0.0f64.to_bits());
        }
        for s in [
            &r.filename,
            &r.run_name,
            &r.transition_group_id,
            &r.sequence,
            &r.full_peptide_name,
            &r.protein_name,
            &r.gene_name,
            &r.aggr_peak_area,
            &r.aggr_peak_apex,
            &r.aggr_fragment_annotation,
            &r.ipf_full_peptide_name,
        ] {
            assert!(s.is_empty());
        }
        for o in [
            r.aggr_prec_peak_area,
            r.aggr_prec_peak_apex,
            r.exp_im,
            r.im_left_width,
            r.im_right_width,
            r.ms1_pep,
            r.ms2_pep,
            r.precursor_pep,
            r.ipf_pep,
            r.pep,
            r.ms2_m_score,
            r.peptide_run_specific_qvalue,
            r.peptide_experiment_wide_qvalue,
            r.peptide_global_qvalue,
            r.protein_run_specific_qvalue,
            r.protein_experiment_wide_qvalue,
            r.protein_global_qvalue,
            r.gene_run_specific_qvalue,
            r.gene_experiment_wide_qvalue,
            r.gene_global_qvalue,
            r.alignment_reference_rt,
            r.alignment_pep,
            r.alignment_qvalue,
            r.ipf_precursor_peakgroup_pep,
            r.ipf_peptidoform_pep,
            r.ipf_peptidoform_m_score,
        ] {
            assert!(o.is_none());
        }
        assert!(r.alignment_group_id.is_none() && r.alignment_reference_feature_id.is_none());
    }

    #[test]
    fn score_row_defaults_match_the_header() {
        let r = OpenSwathFeatureScoreRow::default();
        assert_eq!(
            (
                r.protein_id,
                r.peptide_id,
                r.precursor_id,
                r.run_id,
                r.feature_id
            ),
            (-1, -1, -1, -1, -1)
        );
        assert_eq!(r.precursor_charge, 0);
        for v in [
            r.precursor_mz,
            r.exp_rt,
            r.norm_rt,
            r.delta_rt,
            r.left_width,
            r.right_width,
        ] {
            assert_eq!(v.to_bits(), 0.0f64.to_bits());
        }
        assert!(!r.protein_decoy && !r.peptide_decoy && !r.precursor_decoy && !r.from_alignment);
        assert!(r.ipf_peptide_id.is_none() && r.gene_id.is_none() && r.gene_decoy.is_none());
        assert!(r.score_ms1_rank.is_none() && r.score_ms2_peak_group_rank.is_none());
        assert!(r.score_ms2_score.is_none() && r.score_gene_run_specific_pep.is_none());
        assert!(r.alignment_group_id.is_none() && r.alignment_qvalue.is_none());
        assert!(r.feature_ms1_values.is_empty() && r.feature_ms2_values.is_empty());
        assert!(r.protein_accession.is_empty() && r.filename.is_empty());

        let t = OpenSwathTransitionScoreRow::default();
        assert_eq!((t.precursor_id, t.transition_id), (-1, -1));
        assert_eq!((t.transition_charge, t.transition_ordinal), (0, 0));
        assert_eq!(t.product_mz.to_bits(), 0.0f64.to_bits());
        assert!(!t.transition_detecting && !t.transition_decoy);
        assert!(t.run_id.is_none() && t.ipf_peptide_id.is_none() && t.feature_id.is_none());
        assert!(t.transition_library_intensity.is_none() && t.score_transition_rank.is_none());
        assert!(t.feature_transition_values.is_empty());

        let m = OpenSwathQuantMatrix::default();
        assert!(m.identifier_column_names.is_empty() && m.values.is_empty());
    }
}
