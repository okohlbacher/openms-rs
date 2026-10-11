// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Configuration for OpenSWATH confidence inference: error estimation,
//! peptidoform inference (IPF) and peptide/protein/gene level-context
//! inference.
//!
//! Ports `ANALYSIS/OPENSWATH/OpenSwathInferenceConfig.h` (header-only). The
//! three configurations are plain value types in the source, not
//! `DefaultParamHandler`s, and stay plain structs here. Their defaults mirror
//! the PyProphet command-line defaults, as the source states. See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_INFERENCE_SUPPORT.md).

use crate::analysis::openswath::openswath_inference_data::{InferenceContext, InferenceLevel};
use crate::math::multiple_testing::{LfdrTransform, Pi0Method};

/// Configuration for OpenSWATH analyte-level global and local error
/// estimation, the source's `ErrorEstimationConfig`.
///
/// # Differences from the source
///
/// The source stores the pi0 method and the local-FDR transformation as
/// strings (`"bootstrap"`, `"probit"`) and converts them with
/// `MultipleTesting::toPi0Method` / `toLfdrTransform` inside every inference
/// call, which throws `std::invalid_argument` for an unknown name. The port
/// stores the enums, so an unknown name is refused once, by
/// [`Pi0Method::parse`] or [`LfdrTransform::parse`], when the configuration is
/// built.
#[derive(Clone, Debug, PartialEq)]
pub struct ErrorEstimationConfig {
    /// If true, estimate p-values parametrically (normal null,
    /// `MultipleTesting::pNorm`); otherwise empirically (`pEmp`). Default
    /// `false`.
    pub parametric: bool,
    /// If true, compute the positive false discovery rate (pFDR) instead of
    /// the FDR. Default `false`.
    pub pfdr: bool,
    /// Lambda grid for pi0 estimation. Default
    /// `0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45`.
    pub pi0_lambda: Vec<f64>,
    /// Method for choosing the pi0 tuning parameter (smoother or bootstrap).
    /// Default [`Pi0Method::Bootstrap`], the source's `"bootstrap"`.
    pub pi0_method: Pi0Method,
    /// Degrees of freedom for smoother-based pi0 estimation. Default `3`.
    pub pi0_smooth_df: i32,
    /// If true, smooth `log(pi0)` instead of pi0 directly. Default `false`.
    pub pi0_smooth_log_pi0: bool,
    /// If true, local FDR values above one are truncated to one. Default
    /// `true`.
    pub lfdr_truncate: bool,
    /// If true, enforce non-decreasing local FDR values with increasing
    /// p-values. Default `true`.
    pub lfdr_monotone: bool,
    /// Transformation applied to p-values for local FDR estimation (probit or
    /// logit). Default [`LfdrTransform::Probit`], the source's `"probit"`.
    pub lfdr_transformation: LfdrTransform,
    /// Smoothing bandwidth multiplier used during local FDR density
    /// estimation. Default `1.5`.
    pub lfdr_adj: f64,
    /// Tail clipping threshold used during local FDR calculation. Default
    /// `1e-8`.
    pub lfdr_eps: f64,
}

impl Default for ErrorEstimationConfig {
    fn default() -> Self {
        Self {
            parametric: false,
            pfdr: false,
            pi0_lambda: vec![0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45],
            pi0_method: Pi0Method::Bootstrap,
            pi0_smooth_df: 3,
            pi0_smooth_log_pi0: false,
            lfdr_truncate: true,
            lfdr_monotone: true,
            lfdr_transformation: LfdrTransform::Probit,
            lfdr_adj: 1.5,
            lfdr_eps: 1e-8,
        }
    }
}

/// Configuration for peptidoform inference, the source's
/// `PeptidoformInferenceConfig`.
///
/// Four fields — [`Self::ipf_h0`], [`Self::ipf_max_peakgroup_pep`],
/// [`Self::ipf_max_transition_pep`] and
/// [`Self::ipf_min_alignment_mapping_confidence`] — are read by the source's
/// `OSWFile` IPF readers to select rows, not by
/// [`OpenSwathPeptidoformInference`](crate::analysis::openswath::openswath_peptidoform_inference::OpenSwathPeptidoformInference)
/// itself; they are carried so a reader port can take the same value.
#[derive(Clone, Debug, PartialEq)]
pub struct PeptidoformInferenceConfig {
    /// If true, use MS1 precursor data for IPF. Default `false`.
    pub ipf_ms1_scoring: bool,
    /// If true, use MS2 precursor data for IPF. Default `false`.
    pub ipf_ms2_scoring: bool,
    /// If true, include the null hypothesis that no peptidoform explains the
    /// peak group. Default `true`.
    pub ipf_h0: bool,
    /// If true, compute the model FDR separately per number of peptidoforms
    /// (grouped) instead of pooled. Default `false`.
    pub ipf_grouped_fdr: bool,
    /// Maximum peak-group PEP to consider scored peak groups in IPF. Default
    /// `0.7`.
    pub ipf_max_peakgroup_pep: f64,
    /// Maximum precursor PEP to consider scored precursor evidence in IPF;
    /// evidence at or above it is treated as absent. Default `0.7`.
    pub ipf_max_precursor_pep: f64,
    /// Maximum precursor-peak-group PEP retained for transition-level IPF;
    /// features at or above it are dropped. Default `0.4`.
    pub ipf_max_precursor_peakgroup_pep: f64,
    /// Maximum transition PEP to consider scored transition evidence in IPF.
    /// Default `0.6`.
    pub ipf_max_transition_pep: f64,
    /// If true, propagate confident transition evidence across aligned runs.
    /// Default `false`.
    pub propagate_signal_across_runs: bool,
    /// Minimum mapping confidence retained from
    /// `FEATURE_MS2_ALIGNMENT_CANDIDATE`. Default `0.5`.
    pub ipf_min_alignment_mapping_confidence: f64,
    /// Maximum transition PEP eligible for across-run evidence propagation
    /// (inclusive). Default `0.5`.
    pub across_run_confidence_threshold: f64,
}

impl Default for PeptidoformInferenceConfig {
    fn default() -> Self {
        Self {
            ipf_ms1_scoring: false,
            ipf_ms2_scoring: false,
            ipf_h0: true,
            ipf_grouped_fdr: false,
            ipf_max_peakgroup_pep: 0.7,
            ipf_max_precursor_pep: 0.7,
            ipf_max_precursor_peakgroup_pep: 0.4,
            ipf_max_transition_pep: 0.6,
            propagate_signal_across_runs: false,
            ipf_min_alignment_mapping_confidence: 0.5,
            across_run_confidence_threshold: 0.5,
        }
    }
}

/// Configuration for peptide-, protein- and gene-level context inference,
/// the source's `LevelContextInferenceConfig`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LevelContextInferenceConfig {
    /// The analyte level to perform inference on. Default
    /// [`InferenceLevel::Peptide`].
    pub level: InferenceLevel,
    /// The context in which to perform inference (global, experiment-wide or
    /// run-specific). Default [`InferenceContext::RunSpecific`].
    pub context: InferenceContext,
    /// Error-estimation settings applied during inference.
    pub error: ErrorEstimationConfig,
}
