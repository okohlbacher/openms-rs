// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The row, table and option types of the OpenSWATH result exports.
//!
//! Ports `ANALYSIS/OPENSWATH/OpenSwathExportData.h` (header only): the export
//! options (file format, IPF mode, quantification-matrix level and
//! normalisation, each with its `toString`), the filtered result row the
//! TSV and matrix exports share, the feature-level and transition-level
//! scored rows and tables of the Parquet export
//! ([`OpenSwathParquetExporter`](crate::analysis::openswath::openswath_parquet_exporter::OpenSwathParquetExporter)),
//! and the dense quantification matrix. Every default is the header's member
//! initialiser; an `std::optional` is an [`Option`].
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_PARQUET_SUPPORT.md).

use std::fmt;

/// Output file format of an OpenSWATH export table (source
/// `OpenSwathExportFileFormat`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum OpenSwathExportFileFormat {
    /// Tab-separated text.
    #[default]
    Tsv,
    /// Apache Parquet.
    Parquet,
}

/// IPF handling of an OpenSWATH export (source `OpenSwathIPFExportMode`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum OpenSwathIPFExportMode {
    /// Report results on peptidoform level.
    #[default]
    Peptidoform,
    /// Augment the OpenSWATH results with IPF scores.
    Augmented,
    /// Ignore IPF scores; standard OpenSWATH results only.
    Disable,
}

/// Level of a quantification matrix export (source `OpenSwathMatrixLevel`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum OpenSwathMatrixLevel {
    /// One row per precursor.
    #[default]
    Precursor,
    /// One row per peptide.
    Peptide,
    /// One row per protein.
    Protein,
    /// One row per gene.
    Gene,
}

/// Normalisation of a quantification matrix export (source
/// `OpenSwathMatrixNormalization`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum OpenSwathMatrixNormalization {
    /// No normalisation.
    #[default]
    None,
    /// Median normalisation.
    Median,
    /// Median of medians normalisation.
    MedianMedian,
}

impl OpenSwathExportFileFormat {
    /// The source's `toString`: `tsv` or `parquet`.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Tsv => "tsv",
            Self::Parquet => "parquet",
        }
    }
}

impl OpenSwathIPFExportMode {
    /// The source's `toString`: `peptidoform`, `augmented` or `disable`.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Peptidoform => "peptidoform",
            Self::Augmented => "augmented",
            Self::Disable => "disable",
        }
    }
}

impl OpenSwathMatrixLevel {
    /// The source's `toString`: `precursor`, `peptide`, `protein` or `gene`.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Precursor => "precursor",
            Self::Peptide => "peptide",
            Self::Protein => "protein",
            Self::Gene => "gene",
        }
    }
}

impl OpenSwathMatrixNormalization {
    /// The source's `toString`: `none`, `median` or `medianmedian`.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Median => "median",
            Self::MedianMedian => "medianmedian",
        }
    }
}

macro_rules! display_as_str {
    ($($t:ty),*) => {
        $(impl fmt::Display for $t {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(self.as_str())
            }
        })*
    };
}
display_as_str!(
    OpenSwathExportFileFormat,
    OpenSwathIPFExportMode,
    OpenSwathMatrixLevel,
    OpenSwathMatrixNormalization
);

/// A filtered feature-level result row, shared by the results and matrix
/// exports (source `OpenSwathExportRow`).
///
/// It is kept separate from the Parquet row types because those also carry
/// the dynamically discovered `FEATURE_MS1` / `FEATURE_MS2` /
/// `FEATURE_TRANSITION` columns. Defaults: ids -1, numbers 0, text empty,
/// optional values absent, flags false.
#[derive(Clone, Debug, PartialEq)]
pub struct OpenSwathExportRow {
    pub run_id: i64,
    pub filename: String,
    pub run_name: String,
    pub feature_id: i64,
    pub peptide_id: i64,
    pub precursor_id: i64,
    pub transition_group_id: String,
    pub decoy: bool,
    pub sequence: String,
    pub full_peptide_name: String,
    pub protein_name: String,
    pub gene_name: String,
    pub charge: i32,
    pub mz: f64,
    pub rt: f64,
    pub assay_rt: f64,
    pub delta_rt: f64,
    pub irt: f64,
    pub assay_irt: f64,
    pub delta_irt: f64,
    pub intensity: f64,
    pub aggr_prec_peak_area: Option<f64>,
    pub aggr_prec_peak_apex: Option<f64>,
    pub left_width: f64,
    pub right_width: f64,
    pub exp_im: Option<f64>,
    pub im_left_width: Option<f64>,
    pub im_right_width: Option<f64>,
    pub ms1_pep: Option<f64>,
    pub ms2_pep: Option<f64>,
    pub precursor_pep: Option<f64>,
    pub ipf_pep: Option<f64>,
    pub peak_group_rank: i32,
    pub d_score: f64,
    pub m_score: f64,
    pub pep: Option<f64>,
    pub ms2_m_score: Option<f64>,
    pub peptide_run_specific_qvalue: Option<f64>,
    pub peptide_experiment_wide_qvalue: Option<f64>,
    pub peptide_global_qvalue: Option<f64>,
    pub protein_run_specific_qvalue: Option<f64>,
    pub protein_experiment_wide_qvalue: Option<f64>,
    pub protein_global_qvalue: Option<f64>,
    pub gene_run_specific_qvalue: Option<f64>,
    pub gene_experiment_wide_qvalue: Option<f64>,
    pub gene_global_qvalue: Option<f64>,
    pub alignment_group_id: Option<i64>,
    pub alignment_reference_feature_id: Option<i64>,
    pub alignment_reference_rt: Option<f64>,
    pub alignment_pep: Option<f64>,
    pub alignment_qvalue: Option<f64>,
    pub from_alignment: bool,
    pub aggr_peak_area: String,
    pub aggr_peak_apex: String,
    pub aggr_fragment_annotation: String,
    pub ipf_full_peptide_name: String,
    pub ipf_precursor_peakgroup_pep: Option<f64>,
    pub ipf_peptidoform_pep: Option<f64>,
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
///
/// `feature_ms1_values` and `feature_ms2_values` are aligned with the
/// table's [`OpenSwathFeatureScoreTable::feature_ms1_column_names`] and
/// [`OpenSwathFeatureScoreTable::feature_ms2_column_names`]; a shorter vector
/// leaves the remaining columns null. Defaults: ids -1, numbers 0, text
/// empty, optional values absent, flags false.
#[derive(Clone, Debug, PartialEq)]
pub struct OpenSwathFeatureScoreRow {
    pub protein_id: i64,
    pub peptide_id: i64,
    pub ipf_peptide_id: Option<i64>,
    pub precursor_id: i64,
    pub protein_accession: String,
    pub unmodified_sequence: String,
    pub modified_sequence: String,
    pub precursor_traml_id: String,
    pub precursor_group_label: String,
    pub precursor_mz: f64,
    pub precursor_charge: i32,
    pub precursor_library_intensity: Option<f64>,
    pub precursor_library_rt: Option<f64>,
    pub precursor_library_drift_time: Option<f64>,
    pub gene_id: Option<i64>,
    pub gene_name: String,
    pub gene_decoy: Option<bool>,
    pub protein_decoy: bool,
    pub peptide_decoy: bool,
    pub precursor_decoy: bool,
    pub run_id: i64,
    pub filename: String,
    pub feature_id: i64,
    pub exp_rt: f64,
    pub exp_im: Option<f64>,
    pub norm_rt: f64,
    pub delta_rt: f64,
    pub left_width: f64,
    pub right_width: f64,
    pub im_left_width: Option<f64>,
    pub im_right_width: Option<f64>,
    pub feature_ms1_values: Vec<f64>,
    pub feature_ms2_values: Vec<f64>,
    pub score_ms1_score: Option<f64>,
    pub score_ms1_rank: Option<i32>,
    pub score_ms1_pvalue: Option<f64>,
    pub score_ms1_qvalue: Option<f64>,
    pub score_ms1_pep: Option<f64>,
    pub score_ms2_score: Option<f64>,
    pub score_ms2_peak_group_rank: Option<i32>,
    pub score_ms2_pvalue: Option<f64>,
    pub score_ms2_qvalue: Option<f64>,
    pub score_ms2_pep: Option<f64>,
    pub score_ipf_precursor_peakgroup_pep: Option<f64>,
    pub score_ipf_pep: Option<f64>,
    pub score_ipf_qvalue: Option<f64>,
    pub score_peptide_global_score: Option<f64>,
    pub score_peptide_global_pvalue: Option<f64>,
    pub score_peptide_global_qvalue: Option<f64>,
    pub score_peptide_global_pep: Option<f64>,
    pub score_peptide_experiment_wide_score: Option<f64>,
    pub score_peptide_experiment_wide_pvalue: Option<f64>,
    pub score_peptide_experiment_wide_qvalue: Option<f64>,
    pub score_peptide_experiment_wide_pep: Option<f64>,
    pub score_peptide_run_specific_score: Option<f64>,
    pub score_peptide_run_specific_pvalue: Option<f64>,
    pub score_peptide_run_specific_qvalue: Option<f64>,
    pub score_peptide_run_specific_pep: Option<f64>,
    pub score_protein_global_score: Option<f64>,
    pub score_protein_global_pvalue: Option<f64>,
    pub score_protein_global_qvalue: Option<f64>,
    pub score_protein_global_pep: Option<f64>,
    pub score_protein_experiment_wide_score: Option<f64>,
    pub score_protein_experiment_wide_pvalue: Option<f64>,
    pub score_protein_experiment_wide_qvalue: Option<f64>,
    pub score_protein_experiment_wide_pep: Option<f64>,
    pub score_protein_run_specific_score: Option<f64>,
    pub score_protein_run_specific_pvalue: Option<f64>,
    pub score_protein_run_specific_qvalue: Option<f64>,
    pub score_protein_run_specific_pep: Option<f64>,
    pub score_gene_global_score: Option<f64>,
    pub score_gene_global_pvalue: Option<f64>,
    pub score_gene_global_qvalue: Option<f64>,
    pub score_gene_global_pep: Option<f64>,
    pub score_gene_experiment_wide_score: Option<f64>,
    pub score_gene_experiment_wide_pvalue: Option<f64>,
    pub score_gene_experiment_wide_qvalue: Option<f64>,
    pub score_gene_experiment_wide_pep: Option<f64>,
    pub score_gene_run_specific_score: Option<f64>,
    pub score_gene_run_specific_pvalue: Option<f64>,
    pub score_gene_run_specific_qvalue: Option<f64>,
    pub score_gene_run_specific_pep: Option<f64>,
    pub alignment_group_id: Option<i64>,
    pub alignment_reference_feature_id: Option<i64>,
    pub alignment_reference_rt: Option<f64>,
    pub alignment_pep: Option<f64>,
    pub alignment_qvalue: Option<f64>,
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
    /// Names of the `FEATURE_MS1_<name>` columns.
    pub feature_ms1_column_names: Vec<String>,
    /// Names of the `FEATURE_MS2_<name>` columns.
    pub feature_ms2_column_names: Vec<String>,
    /// The rows.
    pub rows: Vec<OpenSwathFeatureScoreRow>,
}

/// A transition-level scored row of the optional Parquet export (source
/// `OpenSwathTransitionScoreRow`).
///
/// `feature_transition_values` is aligned with
/// [`OpenSwathTransitionScoreTable::feature_transition_column_names`]; a
/// missing or absent entry is null. Defaults: ids -1, numbers 0, text
/// empty, optional values absent, flags false.
#[derive(Clone, Debug, PartialEq)]
pub struct OpenSwathTransitionScoreRow {
    pub run_id: Option<i64>,
    pub ipf_peptide_id: Option<i64>,
    pub precursor_id: i64,
    pub transition_id: i64,
    pub transition_traml_id: String,
    pub product_mz: f64,
    pub transition_charge: i32,
    pub transition_type: String,
    pub transition_ordinal: i32,
    pub annotation: String,
    pub transition_detecting: bool,
    pub transition_library_intensity: Option<f64>,
    pub transition_decoy: bool,
    pub feature_id: Option<i64>,
    pub feature_transition_values: Vec<Option<f64>>,
    pub score_transition_score: Option<f64>,
    pub score_transition_rank: Option<i32>,
    pub score_transition_pvalue: Option<f64>,
    pub score_transition_qvalue: Option<f64>,
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

/// The transition-level Parquet export table with its dynamically
/// discovered `FEATURE_TRANSITION` columns (source
/// `OpenSwathTransitionScoreTable`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OpenSwathTransitionScoreTable {
    /// Names of the `FEATURE_TRANSITION_<name>` columns.
    pub feature_transition_column_names: Vec<String>,
    /// The rows.
    pub rows: Vec<OpenSwathTransitionScoreRow>,
}

/// The dense matrix of the matrix exports (source `OpenSwathQuantMatrix`).
///
/// Values are stored by row; a missing value is `None`, which preserves
/// sparse run coverage in TSV and Parquet output.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OpenSwathQuantMatrix {
    /// Names of the identifier columns.
    pub identifier_column_names: Vec<String>,
    /// Identifier values, one vector per row.
    pub identifier_rows: Vec<Vec<String>>,
    /// Names of the sample columns.
    pub sample_column_names: Vec<String>,
    /// Sample values, one vector per row.
    pub values: Vec<Vec<Option<f64>>>,
}
