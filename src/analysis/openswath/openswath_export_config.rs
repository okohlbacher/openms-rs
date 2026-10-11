// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Filter and output configuration of the OpenSWATH exports.
//!
//! Ports `ANALYSIS/OPENSWATH/OpenSwathExportConfig.h`: the shared filter
//! configuration [`OpenSwathExportFilterConfig`], whose defaults mirror the
//! PyProphet export command line, and the three output configurations
//! [`OpenSwathResultsExportConfig`], [`OpenSwathMatrixExportConfig`] and
//! [`OpenSwathParquetExportConfig`].
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_EXPORT_SUPPORT.md).
//!
//! # Example
//!
//! ```
//! use openms::analysis::openswath::openswath_export_config::OpenSwathMatrixExportConfig;
//! use openms::analysis::openswath::openswath_export_data::OpenSwathMatrixLevel;
//!
//! let config = OpenSwathMatrixExportConfig::default();
//! assert_eq!(config.level, OpenSwathMatrixLevel::Peptide);
//! assert_eq!(config.top_n, 3);
//! assert!(config.consistent_top);
//! assert_eq!(config.filters.max_rs_peakgroup_qvalue, 0.05);
//! ```

use super::openswath_export_data::{
    OpenSwathExportFileFormat, OpenSwathIPFExportMode, OpenSwathMatrixLevel,
    OpenSwathMatrixNormalization,
};

/// The shared filter configuration of the OpenSWATH export readers (source
/// `OpenSwathExportFilterConfig`).
///
/// The defaults mirror the defaults of the corresponding PyProphet export
/// command line. The exporters of this crate do not read these values: they
/// belong to the OSW readers (`FORMAT/OSWFile.h`) that produce the rows the
/// exporters write.
#[derive(Clone, Debug, PartialEq)]
pub struct OpenSwathExportFilterConfig {
    /// IPF handling. Default [`OpenSwathIPFExportMode::Peptidoform`].
    pub ipf_mode: OpenSwathIPFExportMode,
    /// Report the aggregated transition-level quantification. Default `true`.
    pub transition_quantification: bool,
    /// Maximum posterior error probability to retain scored transitions for
    /// the quantification (requires transition-level scoring). Default 0.7.
    pub max_transition_pep: f64,
    /// Filter the results to this maximum run-specific peptidoform-level
    /// posterior error probability. Default 0.4.
    pub ipf_max_peptidoform_pep: f64,
    /// Filter the results to this maximum run-specific peak group-level
    /// q-value. Default 0.05.
    pub max_rs_peakgroup_qvalue: f64,
    /// Filter on the global peptide q-value. Default `true`.
    pub peptide: bool,
    /// Maximum global peptide q-value. Default 0.01.
    pub max_global_peptide_qvalue: f64,
    /// Filter on the global protein q-value. Default `true`.
    pub protein: bool,
    /// Maximum global protein q-value. Default 0.01.
    pub max_global_protein_qvalue: f64,
    /// Filter on the global gene q-value. Default `false`.
    pub gene: bool,
    /// Maximum global gene q-value. Default 0.01.
    pub max_global_gene_qvalue: f64,
    /// Use the alignment results. Default `false`.
    pub use_alignment: bool,
    /// Maximum alignment posterior error probability. Default 0.7.
    pub max_alignment_pep: f64,
    /// Exclude decoys. Default `true`.
    pub exclude_decoys: bool,
}

impl Default for OpenSwathExportFilterConfig {
    fn default() -> Self {
        Self {
            ipf_mode: OpenSwathIPFExportMode::Peptidoform,
            transition_quantification: true,
            max_transition_pep: 0.7,
            ipf_max_peptidoform_pep: 0.4,
            max_rs_peakgroup_qvalue: 0.05,
            peptide: true,
            max_global_peptide_qvalue: 0.01,
            protein: true,
            max_global_protein_qvalue: 0.01,
            gene: false,
            max_global_gene_qvalue: 0.01,
            use_alignment: false,
            max_alignment_pep: 0.7,
            exclude_decoys: true,
        }
    }
}

/// The configuration of a user-facing feature results export (source
/// `OpenSwathResultsExportConfig`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OpenSwathResultsExportConfig {
    /// The row filters.
    pub filters: OpenSwathExportFilterConfig,
    /// The output format. Default [`OpenSwathExportFileFormat::Tsv`].
    pub format: OpenSwathExportFileFormat,
}

/// The configuration of a matrix export (source
/// `OpenSwathMatrixExportConfig`).
#[derive(Clone, Debug, PartialEq)]
pub struct OpenSwathMatrixExportConfig {
    /// The row filters.
    pub filters: OpenSwathExportFilterConfig,
    /// The summarization level. Default [`OpenSwathMatrixLevel::Peptide`].
    pub level: OpenSwathMatrixLevel,
    /// The normalization. Default [`OpenSwathMatrixNormalization::None`].
    pub normalization: OpenSwathMatrixNormalization,
    /// How many precursors (peptide level) or peptides (protein and gene
    /// level) are summarized. Default 3; zero is refused above the precursor
    /// level.
    pub top_n: usize,
    /// Pick the same top precursors or peptides in every run. Default `true`.
    pub consistent_top: bool,
    /// The output format. Default [`OpenSwathExportFileFormat::Tsv`].
    pub format: OpenSwathExportFileFormat,
}

impl Default for OpenSwathMatrixExportConfig {
    fn default() -> Self {
        Self {
            filters: OpenSwathExportFilterConfig::default(),
            level: OpenSwathMatrixLevel::Peptide,
            normalization: OpenSwathMatrixNormalization::None,
            top_n: 3,
            consistent_top: true,
            format: OpenSwathExportFileFormat::Tsv,
        }
    }
}

/// The configuration of a scored feature Parquet export (source
/// `OpenSwathParquetExportConfig`).
#[derive(Clone, Debug, PartialEq)]
pub struct OpenSwathParquetExportConfig {
    /// The row filters.
    pub filters: OpenSwathExportFilterConfig,
    /// Also export the transition-level data. Default `true`.
    pub include_transition_data: bool,
}

impl Default for OpenSwathParquetExportConfig {
    fn default() -> Self {
        Self {
            filters: OpenSwathExportFilterConfig::default(),
            include_transition_data: true,
        }
    }
}
