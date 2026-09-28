// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Value types exchanged with the in-process Percolator.
//!
//! Port of `ANALYSIS/ID/PercolatorTypes.h` (header only, no `.cpp`) at
//! `bc9cc12`. The three source structs keep their names and fields;
//! [`RescoreInput::validate`] and [`PercolatorModel::validate`] state as
//! checks the row and column contracts that the source only documents.
//! `docs/PERCOLATOR_SUPPORT.md` maps every member.

use crate::{Error, Result};

/// Input to the domain-agnostic Percolator rescoring.
///
/// Row ordering is preserved in the output. Each row corresponds to one data
/// point; the row's semantics (PSM, transition, peak group, …) are
/// determined by the caller.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RescoreInput {
    /// `[n_rows][n_features]` scalar features per row. Rows must all have the
    /// same length.
    pub features: Vec<Vec<f64>>,
    /// Target (`false`) or decoy (`true`) label per row.
    pub is_decoy: Vec<bool>,
    /// Per-row integer key grouping rows into the same cross-validation fold.
    ///
    /// Rows sharing a key are never split across folds. Leave empty to use the
    /// row index (each row in its own group). Supply this when rows have
    /// natural duplication, such as several PSMs from one spectrum or several
    /// transitions from one precursor.
    pub cv_group_keys: Vec<i32>,
    /// Names aligned 1:1 with the feature columns; the source uses them for
    /// logging only.
    pub feature_names: Vec<String>,
    /// Optional per-row PIN-compatible scan numbers.
    ///
    /// The four PIN-compatible vectors (`scan_numbers`, `spec_file_numbers`,
    /// `exp_masses`, `calc_masses`) override, when supplied, the synthetic
    /// defaults the wrapper uses internally (`scan = row_index`,
    /// `specFileNr = 0`, `expMass = calcMass = 0.0`). Populate them when the
    /// in-process output must be bitwise-identical to running the external
    /// percolator binary on a `.pin` file derived from the same PSMs:
    /// Percolator's internal PSM sort order (`OrderScanHash` hashes
    /// `specFileNr` and `scan`) determines the cross-validation fold
    /// assignment, so the scan and spectrum-file numbers must match those
    /// that `PercolatorInfile::store` would emit. The source helper
    /// `Percolator::fillPINCompatibleFields` derives all four vectors from
    /// peptide identifications with the same conventions as
    /// `PercolatorInfile::store`.
    ///
    /// Each vector is either empty (the default is used) or has exactly one
    /// entry per feature row.
    pub scan_numbers: Vec<i32>,
    /// Optional per-row spectrum-file numbers; see [`RescoreInput::scan_numbers`].
    pub spec_file_numbers: Vec<i32>,
    /// Optional per-row experimental masses; see [`RescoreInput::scan_numbers`].
    pub exp_masses: Vec<f64>,
    /// Optional per-row calculated masses; see [`RescoreInput::scan_numbers`].
    pub calc_masses: Vec<f64>,
}

impl RescoreInput {
    /// Number of rows, taken from [`RescoreInput::features`].
    pub fn n_rows(&self) -> usize {
        self.features.len()
    }

    /// Number of feature columns: the length of the first row, or zero for
    /// an empty input.
    pub fn n_features(&self) -> usize {
        self.features.first().map_or(0, Vec::len)
    }

    /// Checks the row and column contract the source documents.
    ///
    /// Native check: the source states these as field documentation and
    /// checks a subset of them inside `Percolator::rescore`, `train` and
    /// `score`, each throwing `Exception::InvalidValue`. Here the full
    /// contract is one function every entry point calls first.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a row has a different length from the
    /// first, a feature is not finite, `is_decoy` does not have one entry per
    /// row, `feature_names` is non-empty and does not name every column,
    /// `cv_group_keys` or one of the four PIN-compatible vectors is neither
    /// empty nor one entry per row, or a mass is not finite.
    pub fn validate(&self) -> Result<()> {
        let rows = self.n_rows();
        let columns = self.n_features();
        for (index, row) in self.features.iter().enumerate() {
            if row.len() != columns {
                return Err(invalid(format!(
                    "feature row {index} has {} columns, expected {columns}",
                    row.len()
                )));
            }
            if row.iter().any(|value| !value.is_finite()) {
                return Err(invalid(format!("feature row {index} is not finite")));
            }
        }
        if self.is_decoy.len() != rows {
            return Err(invalid(format!(
                "is_decoy has {} entries for {rows} rows",
                self.is_decoy.len()
            )));
        }
        if !self.feature_names.is_empty() && self.feature_names.len() != columns {
            return Err(invalid(format!(
                "feature_names has {} entries for {columns} columns",
                self.feature_names.len()
            )));
        }
        optional_rows("cv_group_keys", self.cv_group_keys.len(), rows)?;
        optional_rows("scan_numbers", self.scan_numbers.len(), rows)?;
        optional_rows("spec_file_numbers", self.spec_file_numbers.len(), rows)?;
        optional_rows("exp_masses", self.exp_masses.len(), rows)?;
        optional_rows("calc_masses", self.calc_masses.len(), rows)?;
        if self
            .exp_masses
            .iter()
            .chain(&self.calc_masses)
            .any(|value| !value.is_finite())
        {
            return Err(invalid("PIN-compatible masses must be finite"));
        }
        Ok(())
    }
}

/// Output of Percolator rescoring, aligned 1:1 with
/// [`RescoreInput::features`].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RescoreOutput {
    /// SVM discriminant score per row.
    pub scores: Vec<f64>,
    /// q-value per row.
    pub q_values: Vec<f64>,
    /// Posterior error probability per row.
    pub peps: Vec<f64>,
}

/// Trained Percolator model: averaged SVM weights in raw feature space.
///
/// Produced by Percolator training and consumed by Percolator scoring. The
/// weights are un-normalized: they multiply raw input features directly. The
/// normalization transform learned by the SVM has already been folded into
/// the weights and bias (the source's `Normalizer::unnormalizeweight()`), so
/// callers must not normalize features before scoring; doing so would apply
/// the transform twice.
///
/// The raw SVM dot product for a row with feature vector `f` is
///
/// ```text
/// raw = sum_j(f[j] * weights[j]) + weights[n_features]     // bias last
/// ```
///
/// Scoring applies a further FDR-based rescaling on top of this raw value to
/// produce the final SVM discriminant reported in [`RescoreOutput::scores`].
#[derive(Clone, Debug, PartialEq)]
pub struct PercolatorModel {
    /// Integer schema version for the on-disk format; `1` by default.
    pub format_version: i32,
    /// Feature column names.
    ///
    /// Must be non-empty and must match [`RescoreInput::feature_names`]
    /// positionally at score time. Any string value is permitted: the bias is
    /// stored in the model file's header, so feature names carry no reserved
    /// meaning.
    pub feature_names: Vec<String>,
    /// `n_features + 1` weights; the last entry is the bias.
    pub weights: Vec<f64>,
    /// `"stdv"`, `"uni"` or `"none"`: the normalizer used during training.
    ///
    /// Informational only: all three produce raw-space weights that scoring
    /// applies directly, since the normalization transform is already folded
    /// into the weights and bias. Recorded so that reproducibility tooling
    /// can identify the learner configuration that produced the model.
    pub normalizer_type: String,
    /// Random seed used during training; informational. `0` by default.
    pub seed: i32,
}

impl Default for PercolatorModel {
    /// The source's in-class initializers: `format_version = 1`, `seed = 0`,
    /// everything else empty.
    fn default() -> Self {
        Self {
            format_version: 1,
            feature_names: Vec::new(),
            weights: Vec::new(),
            normalizer_type: String::new(),
            seed: 0,
        }
    }
}

impl PercolatorModel {
    /// Checks the model invariant `weights.len() == feature_names.len() + 1`,
    /// a non-empty feature list and finite weights.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the invariant does not hold, as the
    /// source `saveModel` and `score` throw `Exception::InvalidValue`.
    pub fn validate(&self) -> Result<()> {
        if self.feature_names.is_empty() {
            return Err(invalid("Percolator model has no feature names"));
        }
        if self.weights.len() != self.feature_names.len() + 1 {
            return Err(invalid(format!(
                "Percolator model has {} weights for {} features (expected features + 1)",
                self.weights.len(),
                self.feature_names.len()
            )));
        }
        if self.weights.iter().any(|value| !value.is_finite()) {
            return Err(invalid("Percolator model weights must be finite"));
        }
        Ok(())
    }

    /// The bias, the last weight; `None` for a model without weights.
    pub fn bias(&self) -> Option<f64> {
        self.weights.last().copied()
    }
}

fn optional_rows(name: &str, len: usize, rows: usize) -> Result<()> {
    if len == 0 || len == rows {
        Ok(())
    } else {
        Err(invalid(format!(
            "{name} has {len} entries for {rows} rows (must be empty or one per row)"
        )))
    }
}

fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}
