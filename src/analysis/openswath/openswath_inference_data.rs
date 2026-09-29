// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Compact typed rows for OpenSWATH confidence inference: the inference level
//! and context, the peptidoform-inference (IPF) evidence and result rows, and
//! the peptide/protein/gene level-context input and result rows.
//!
//! Ports `ANALYSIS/OPENSWATH/OpenSwathInferenceData.h` (header-only). The rows
//! are file-format agnostic: a reader (the source's `OSWFile`) fills them, and
//! [`LevelContextInference`](crate::analysis::openswath::level_context_inference::LevelContextInference)
//! and [`OpenSwathPeptidoformInference`](crate::analysis::openswath::openswath_peptidoform_inference::OpenSwathPeptidoformInference)
//! consume them. See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_INFERENCE_SUPPORT.md).
//!
//! References:
//!
//! > Rosenberger, G., Bludau, I., Schmitt, U. et al. Statistical control of
//! > peptide and protein error rates in large-scale targeted data-independent
//! > acquisition analyses. Nat Methods 14, 921-927 (2017).
//! > <https://doi.org/10.1038/nmeth.4398>
//!
//! > Rosenberger, G., Liu, Y., Rost, H. et al. Inference and quantification of
//! > peptidoforms in large sample cohorts by SWATH-MS. Nat Biotechnol 35,
//! > 781-788 (2017). <https://doi.org/10.1038/nbt.3908>
//!
//! Every field default below is the source's member initialiser.

use std::fmt;

/// Inference level for OpenSWATH-derived confidence inference.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum InferenceLevel {
    /// Peptidoform level (IPF); handled by
    /// [`OpenSwathPeptidoformInference`](crate::analysis::openswath::openswath_peptidoform_inference::OpenSwathPeptidoformInference),
    /// not by level-context inference.
    Peptidoform,
    /// Peptide level. The default of
    /// [`LevelContextInferenceConfig::level`](crate::analysis::openswath::openswath_inference_config::LevelContextInferenceConfig::level).
    #[default]
    Peptide,
    /// Protein level.
    Protein,
    /// Gene level.
    Gene,
}

impl InferenceLevel {
    /// The source's `toString(InferenceLevel)`: `peptidoform`, `peptide`,
    /// `protein` or `gene`.
    ///
    /// The source's trailing `return ""` after an exhaustive `switch` is
    /// unreachable for a Rust enum and has no counterpart.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Peptidoform => "peptidoform",
            Self::Peptide => "peptide",
            Self::Protein => "protein",
            Self::Gene => "gene",
        }
    }
}

impl fmt::Display for InferenceLevel {
    /// Writes [`InferenceLevel::as_str`].
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Context for peptide-, protein- and gene-level inference.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum InferenceContext {
    /// One score distribution over all rows, one best row per entity.
    Global,
    /// One shared score distribution for all rows of the call, which may carry
    /// several runs.
    ExperimentWide,
    /// Error statistics estimated independently per run id. The default.
    #[default]
    RunSpecific,
}

impl InferenceContext {
    /// The source's `toString(InferenceContext)`: `global`,
    /// `experiment-wide` or `run-specific`.
    ///
    /// As [`InferenceLevel::as_str`], the source's unreachable `return ""` has
    /// no counterpart.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Global => "global",
            Self::ExperimentWide => "experiment-wide",
            Self::RunSpecific => "run-specific",
        }
    }
}

impl fmt::Display for InferenceContext {
    /// Writes [`InferenceContext::as_str`].
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Peak-group and precursor evidence for peptidoform inference, the source's
/// `IPFPrecursorRow`.
#[derive(Clone, Debug, PartialEq)]
pub struct IPFPrecursorRow {
    /// Feature (peak group) id; `-1` by default.
    pub feature_id: i64,
    /// Posterior error probability of the MS2 peak group; `1.0` by default.
    pub ms2_peakgroup_pep: f64,
    /// Posterior error probability of the MS1 precursor signal, when scored.
    pub ms1_precursor_pep: Option<f64>,
    /// Posterior error probability of the MS2 precursor signal, when scored.
    pub ms2_precursor_pep: Option<f64>,
}

impl Default for IPFPrecursorRow {
    fn default() -> Self {
        Self {
            feature_id: -1,
            ms2_peakgroup_pep: 1.0,
            ms1_precursor_pep: None,
            ms2_precursor_pep: None,
        }
    }
}

/// Precursor-layer posterior error probability used by IPF, the source's
/// `IPFPrecursorProbabilityRow`.
#[derive(Clone, Debug, PartialEq)]
pub struct IPFPrecursorProbabilityRow {
    /// Feature (peak group) id; `-1` by default.
    pub feature_id: i64,
    /// Precursor peak-group PEP; `1.0` by default. It is the prior of the
    /// transition/peptidoform layer.
    pub precursor_peakgroup_pep: f64,
}

impl Default for IPFPrecursorProbabilityRow {
    fn default() -> Self {
        Self {
            feature_id: -1,
            precursor_peakgroup_pep: 1.0,
        }
    }
}

/// Transition-layer evidence for peptidoform inference, the source's
/// `IPFTransitionRow`.
#[derive(Clone, Debug, PartialEq)]
pub struct IPFTransitionRow {
    /// Feature (peak group) id; `-1` by default.
    pub feature_id: i64,
    /// Transition id; `-1` by default.
    pub transition_id: i64,
    /// Peptide (peptidoform) id this evidence speaks for; `-1` is the H0/null
    /// hypothesis. `-1` by default.
    pub peptide_id: i64,
    /// Transition posterior error probability; `1.0` by default.
    pub pep: f64,
    /// `1` when the transition is expected for this peptidoform (the evidence
    /// is `1 - pep`), anything else when it is not (the evidence is `pep`).
    /// `0` by default.
    pub bmask: i32,
    /// Number of peptidoforms the precursor can explain; must be positive for
    /// every non-null hypothesis. `0` by default.
    pub num_peptidoforms: i32,
    /// Alignment group for across-run propagation, when the row carries one
    /// itself; otherwise the alignment rows supply it.
    pub alignment_group_id: Option<i64>,
}

impl Default for IPFTransitionRow {
    fn default() -> Self {
        Self {
            feature_id: -1,
            transition_id: -1,
            peptide_id: -1,
            pep: 1.0,
            bmask: 0,
            num_peptidoforms: 0,
            alignment_group_id: None,
        }
    }
}

/// Alignment-group membership for cross-run IPF signal propagation, the
/// source's `IPFAlignmentRow`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IPFAlignmentRow {
    /// Alignment group id; `-1` by default. Negative ids mean "no group".
    pub alignment_group_id: i64,
    /// Feature (peak group) id; `-1` by default.
    pub feature_id: i64,
}

impl Default for IPFAlignmentRow {
    fn default() -> Self {
        Self {
            alignment_group_id: -1,
            feature_id: -1,
        }
    }
}

/// Final peptidoform inference result row, the source's `IPFResultRow`.
#[derive(Clone, Debug, PartialEq)]
pub struct IPFResultRow {
    /// Feature (peak group) id; `-1` by default.
    pub feature_id: i64,
    /// Inferred peptidoform (peptide) id; `-1` by default.
    pub peptide_id: i64,
    /// Precursor peak-group PEP of the feature; `1.0` by default.
    pub precursor_peakgroup_pep: f64,
    /// Model q-value of the peptidoform; `1.0` by default.
    pub qvalue: f64,
    /// Posterior error probability of the peptidoform; `1.0` by default.
    pub pep: f64,
}

impl Default for IPFResultRow {
    fn default() -> Self {
        Self {
            feature_id: -1,
            peptide_id: -1,
            precursor_peakgroup_pep: 1.0,
            qvalue: 1.0,
            pep: 1.0,
        }
    }
}

/// Compact level-context input row independent of file format, the source's
/// `LevelContextInputRow`.
#[derive(Clone, Debug, PartialEq)]
pub struct LevelContextInputRow {
    /// Run id; required for run-specific inference.
    pub run_id: Option<i64>,
    /// Group label the reader assigned (for example `"<run>_<entity>"`); not
    /// read by the inference itself.
    pub group_id: String,
    /// Peptide, protein or gene id; `-1` by default.
    pub entity_id: i64,
    /// Whether the row is a decoy; decoys form the null score distribution.
    pub decoy: bool,
    /// Best score of the entity in its context; `0.0` by default.
    pub score: f64,
    /// Context the row was selected for; run-specific by default.
    pub context: InferenceContext,
}

impl Default for LevelContextInputRow {
    fn default() -> Self {
        Self {
            run_id: None,
            group_id: String::new(),
            entity_id: -1,
            decoy: false,
            score: 0.0,
            context: InferenceContext::RunSpecific,
        }
    }
}

/// Final level-context result row, the source's `LevelContextResultRow`.
#[derive(Clone, Debug, PartialEq)]
pub struct LevelContextResultRow {
    /// Run id copied from the input row.
    pub run_id: Option<i64>,
    /// Entity id copied from the input row; `-1` by default.
    pub entity_id: i64,
    /// Context copied from the input row; run-specific by default.
    pub context: InferenceContext,
    /// Score copied from the input row; `0.0` by default.
    pub score: f64,
    /// Estimated p-value; `1.0` by default, `NaN` for a non-finite score.
    pub pvalue: f64,
    /// Estimated q-value; `1.0` by default, `NaN` for a non-finite score.
    pub qvalue: f64,
    /// Estimated posterior error probability (local FDR); `1.0` by default,
    /// `NaN` for a non-finite score.
    pub pep: f64,
}

impl Default for LevelContextResultRow {
    fn default() -> Self {
        Self {
            run_id: None,
            entity_id: -1,
            context: InferenceContext::RunSpecific,
            score: 0.0,
            pvalue: 1.0,
            qvalue: 1.0,
            pep: 1.0,
        }
    }
}
