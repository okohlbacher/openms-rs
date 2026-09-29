// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Shared OpenSWATH context inference for peptide, protein and gene rows.
//!
//! Ports `ANALYSIS/OPENSWATH/LevelContextInference.h` and its translation
//! unit. The class is file-format agnostic and operates only on compact typed
//! rows ([`LevelContextInputRow`]). It is adapted from PyProphet's context
//! inference; see Rosenberger, G., Bludau, I., Schmitt, U. et al. Statistical
//! control of peptide and protein error rates in large-scale targeted
//! data-independent acquisition analyses. Nat Methods 14, 921-927 (2017),
//! <https://doi.org/10.1038/nmeth.4398>.
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_INFERENCE_SUPPORT.md).
//!
//! The source contains no `#pragma omp`; the port is serial.
//!
//! [`LevelContextInputRow`]: crate::analysis::openswath::openswath_inference_data::LevelContextInputRow

use crate::analysis::openswath::openswath_inference_config::{
    ErrorEstimationConfig, LevelContextInferenceConfig,
};
use crate::analysis::openswath::openswath_inference_data::{
    InferenceContext, InferenceLevel, LevelContextInputRow, LevelContextResultRow,
};
use crate::math::multiple_testing::{
    LfdrOptions, MAX_ITEMS, Pi0Smoother, lfdr, p_emp, p_norm, pi0_est, q_value,
};
use crate::processing::spline::BSplineSmoothingSpline;
use crate::{Error, Result};
use std::collections::BTreeMap;

/// Largest number of input rows one [`LevelContextInference::infer`] call
/// accepts: [`MAX_ITEMS`], the ceiling of the statistics it
/// wraps. Native; the source has no ceiling.
pub const MAX_ROWS: usize = MAX_ITEMS;

/// The smoothing spline the source's `MultipleTesting::pi0Est` fits:
/// [`BSplineSmoothingSpline`] with the smoothing parameter `-1.0` ("choose
/// automatically") and the degree `smooth_df`, read at the largest lambda.
///
/// [`pi0_est`] takes the spline through
/// [`Pi0Smoother`] because `math` may not depend on `processing`; this is the
/// adapter above that layer, so smoother-based pi0 estimation here gives the
/// source's answer rather than the lower `min(pi0_lambda)` fallback. A spline
/// that fails to fit or evaluate returns `None`, the source's `!spl.ok()`
/// branch.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct BSplinePi0Smoother;

impl Pi0Smoother for BSplinePi0Smoother {
    fn smooth_eval(&self, x: &[f64], y: &[f64], smooth_df: i32, at: f64) -> Option<f64> {
        BSplineSmoothingSpline::with_smoothing(x, y, -1.0, smooth_df)
            .ok()
            .and_then(|spline| spline.eval(at).ok())
    }
}

/// Shared OpenSWATH analysis-context inference for peptide, protein and gene
/// rows. Stateless; the source's class has only a static member.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct LevelContextInference;

/// Why [`build_error_statistics`] failed: one of the source's own
/// `Exception::Precondition`s, which the run-specific path rewraps with the
/// run id, or an error from the wrapped statistics (`std::invalid_argument` in
/// the source), which it does not catch.
enum StatsError {
    Precondition(String),
    Statistics(Error),
}

impl From<Error> for StatsError {
    fn from(error: Error) -> Self {
        Self::Statistics(error)
    }
}

impl StatsError {
    fn into_error(self) -> Error {
        match self {
            Self::Precondition(message) => Error::InvalidValue(message),
            Self::Statistics(error) => error,
        }
    }
}

/// Error statistics on the ascending finite target scores (the cutoffs), one
/// p-value, q-value and PEP per cutoff.
struct ErrorStatistics {
    cutoffs: Vec<f64>,
    pvalues: Vec<f64>,
    qvalues: Vec<f64>,
    peps: Vec<f64>,
}

impl LevelContextInference {
    /// Estimate p-values, q-values and PEPs for the given input rows.
    ///
    /// The input rows are expected to contain one best-scoring row per entity
    /// and context. Their `decoy` flag defines the target/decoy split used to
    /// build the score distributions for p-value, q-value and local-FDR/PEP
    /// estimation.
    ///
    /// For run-specific inference the target/decoy error statistics are
    /// estimated independently per `run_id`. For global and experiment-wide
    /// inference one shared score distribution is used for all rows in the
    /// call.
    ///
    /// Internally this wraps the multiple-testing utilities for
    /// target/decoy-based p-value estimation ([`p_emp`], or [`p_norm`] when
    /// [`ErrorEstimationConfig::parametric`]), pi0 estimation ([`pi0_est`]
    /// with [`BSplinePi0Smoother`]), q-value estimation ([`q_value`]) and
    /// local-FDR / posterior error probability estimation ([`lfdr`]).
    ///
    /// Each row, target or decoy, is projected onto the nearest target-score
    /// cutoff (the lowest-index one among tied cutoffs) and takes that
    /// cutoff's statistics. A row with a non-finite score keeps its run,
    /// entity, context and score and gets `NaN` p-value, q-value and PEP.
    ///
    /// # Arguments
    ///
    /// * `input` — compact input rows selected for one inference level.
    /// * `config` — inference and error-estimation configuration.
    ///
    /// # Returns
    ///
    /// One result row per input row, in input order, carrying the original
    /// entity and score plus the derived p-value, q-value and PEP. An empty
    /// input gives an empty result without checking anything else but the
    /// level.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] — the source's `Exception::Precondition`
    /// — when `config.level` is [`InferenceLevel::Peptidoform`], when
    /// run-specific inference meets a row without `run_id`, or when a score
    /// distribution (all targets or all decoys of the call, or of one run) is
    /// empty or holds only non-finite scores. In the run-specific context
    /// those distribution errors are prefixed with
    /// `Run-specific inference failed for RUN_ID <id>: ` for the first failing
    /// run in ascending run id order, as the source rewraps them. Errors from
    /// the wrapped statistics (a lambda outside `[0, 1)`, for example) are
    /// passed through unwrapped, as the source's `std::invalid_argument` is
    /// not caught by its `Exception::BaseException` handler. Returns
    /// [`Error::InvalidRange`] when `input` exceeds [`MAX_ROWS`] (native).
    pub fn infer(
        input: &[LevelContextInputRow],
        config: &LevelContextInferenceConfig,
    ) -> Result<Vec<LevelContextResultRow>> {
        if config.level == InferenceLevel::Peptidoform {
            return Err(Error::InvalidValue(
                "Peptidoform inference is not supported by LevelContextInference.".to_string(),
            ));
        }
        if input.is_empty() {
            return Ok(Vec::new());
        }
        if input.len() > MAX_ROWS {
            return Err(Error::InvalidRange(format!(
                "LevelContextInference: {} input rows exceed the ceiling of {MAX_ROWS}",
                input.len()
            )));
        }

        let mut results = vec![LevelContextResultRow::default(); input.len()];
        if config.context == InferenceContext::RunSpecific {
            let mut grouped: BTreeMap<i64, Vec<usize>> = BTreeMap::new();
            for (index, row) in input.iter().enumerate() {
                let Some(run_id) = row.run_id else {
                    return Err(Error::InvalidValue(
                        "Run-specific inference requires RUN_ID values for all rows.".to_string(),
                    ));
                };
                grouped.entry(run_id).or_default().push(index);
            }
            for (run_id, indices) in &grouped {
                assign_error_statistics(input, indices, &config.error, &mut results).map_err(
                    |error| match error {
                        StatsError::Precondition(message) => Error::InvalidValue(format!(
                            "Run-specific inference failed for RUN_ID {run_id}: {message}"
                        )),
                        StatsError::Statistics(error) => error,
                    },
                )?;
            }
        } else {
            let indices: Vec<usize> = (0..input.len()).collect();
            assign_error_statistics(input, &indices, &config.error, &mut results)
                .map_err(StatsError::into_error)?;
        }
        Ok(results)
    }
}

/// The source's file-static `findNearestMatchAscending_`: binary search for
/// the entry of the ascending `basis` nearest to `sample`, then walk left over
/// equal entries so the first of a tied run is returned.
///
/// An exact match found at either end or during the search wins; otherwise the
/// closer of the two final bracket entries, the upper one on an equal
/// distance. `basis` is never empty here (the caller refuses empty target
/// distributions), so the source's empty-basis `Exception::Precondition` is
/// unreachable and has no counterpart.
fn find_nearest_match_ascending(basis: &[f64], sample: f64) -> usize {
    if basis.len() <= 1 {
        return 0;
    }
    let mut low = 0usize;
    let mut high = basis.len() - 1;
    let mut best: Option<usize> = None;
    if basis[low] == sample {
        best = Some(low);
    } else if basis[high] == sample {
        best = Some(high);
    } else {
        while low < high - 1 {
            let mid = (low + high) / 2;
            if basis[mid] == sample {
                best = Some(mid);
            }
            if basis[mid] < sample {
                low = mid;
            } else {
                high = mid;
            }
        }
        if best.is_none() {
            let low_dist = (basis[low] - sample).abs();
            let high_dist = (basis[high] - sample).abs();
            best = Some(if low_dist < high_dist { low } else { high });
        }
    }
    let mut best = best.unwrap_or(0);
    while best > 0 && basis[best - 1] == basis[best] {
        best -= 1;
    }
    best
}

/// The source's file-static `buildErrorStatistics_`.
///
/// The source also computes PyProphet's s-value metrics (`statMetrics_`) and
/// discards them "to keep the construction order"; that computation has no
/// observable effect and is not ported.
fn build_error_statistics(
    target_scores: &[f64],
    decoy_scores: &[f64],
    config: &ErrorEstimationConfig,
) -> std::result::Result<ErrorStatistics, StatsError> {
    if target_scores.is_empty() {
        return Err(StatsError::Precondition(
            "Target score distribution is empty.".to_string(),
        ));
    }
    if decoy_scores.is_empty() {
        return Err(StatsError::Precondition(
            "Decoy score distribution is empty.".to_string(),
        ));
    }
    let mut sorted_targets: Vec<f64> = target_scores
        .iter()
        .copied()
        .filter(|s| s.is_finite())
        .collect();
    sorted_targets.sort_by(f64::total_cmp);
    let mut sorted_decoys: Vec<f64> = decoy_scores
        .iter()
        .copied()
        .filter(|s| s.is_finite())
        .collect();
    sorted_decoys.sort_by(f64::total_cmp);
    if sorted_targets.is_empty() {
        return Err(StatsError::Precondition(
            "Target score distribution only contained NaN values.".to_string(),
        ));
    }
    if sorted_decoys.is_empty() {
        return Err(StatsError::Precondition(
            "Decoy score distribution only contained NaN values.".to_string(),
        ));
    }

    let pvalues = if config.parametric {
        p_norm(&sorted_targets, &sorted_decoys)?
    } else {
        p_emp(&sorted_targets, &sorted_decoys)?
    };
    let pi0 = pi0_est(
        &pvalues,
        &config.pi0_lambda,
        config.pi0_method,
        config.pi0_smooth_df,
        config.pi0_smooth_log_pi0,
        Some(&BSplinePi0Smoother),
    )?;
    let qvalues = q_value(&pvalues, pi0.pi0, config.pfdr)?;
    let options = LfdrOptions {
        truncate: config.lfdr_truncate,
        monotone: config.lfdr_monotone,
        transform: config.lfdr_transformation,
        adj: config.lfdr_adj,
        eps: config.lfdr_eps,
        ..LfdrOptions::default()
    };
    let peps = lfdr(&pvalues, pi0.pi0, &options)?;
    Ok(ErrorStatistics {
        cutoffs: sorted_targets,
        pvalues,
        qvalues,
        peps,
    })
}

/// The source's file-static `assignErrorStatistics_` over the rows `indices`
/// of `input`.
fn assign_error_statistics(
    input: &[LevelContextInputRow],
    indices: &[usize],
    config: &ErrorEstimationConfig,
    results: &mut [LevelContextResultRow],
) -> std::result::Result<(), StatsError> {
    let mut target_scores = Vec::with_capacity(indices.len());
    let mut decoy_scores = Vec::with_capacity(indices.len());
    for &index in indices {
        let row = &input[index];
        if row.decoy {
            decoy_scores.push(row.score);
        } else {
            target_scores.push(row.score);
        }
    }
    let stats = build_error_statistics(&target_scores, &decoy_scores, config)?;
    for &index in indices {
        let row = &input[index];
        let result = &mut results[index];
        result.run_id = row.run_id;
        result.entity_id = row.entity_id;
        result.context = row.context;
        result.score = row.score;
        if !row.score.is_finite() {
            result.pvalue = f64::NAN;
            result.qvalue = f64::NAN;
            result.pep = f64::NAN;
            continue;
        }
        let matched = find_nearest_match_ascending(&stats.cutoffs, row.score);
        result.pvalue = stats.pvalues[matched];
        result.qvalue = stats.qvalues[matched];
        result.pep = stats.peps[matched];
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::find_nearest_match_ascending;

    #[test]
    fn nearest_match_prefers_first_of_ties_and_upper_on_equal_distance() {
        let basis = [1.0, 2.0, 2.0, 2.0, 5.0];
        assert_eq!(find_nearest_match_ascending(&basis, 2.0), 1);
        assert_eq!(find_nearest_match_ascending(&basis, 1.0), 0);
        assert_eq!(find_nearest_match_ascending(&basis, 5.0), 4);
        assert_eq!(find_nearest_match_ascending(&basis, 0.0), 0);
        assert_eq!(find_nearest_match_ascending(&basis, 9.0), 4);
        // 3.5 is equidistant from 2.0 and 5.0: the upper bracket wins, per the
        // source's strict `low_dist < high_dist`.
        assert_eq!(find_nearest_match_ascending(&basis, 3.5), 4);
        // 2.4 is nearer 2.0; the walk-left returns the first tied entry.
        assert_eq!(find_nearest_match_ascending(&basis, 2.4), 1);
        assert_eq!(find_nearest_match_ascending(&[7.0], -3.0), 0);
    }
}
