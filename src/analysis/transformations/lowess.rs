// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Based on OpenMS FastLowessSmoothing, whose source credits the BSD common-lisp-stat
// translation of W. S. Cleveland's NETLIB LOWESS algorithm.
// Upstream attribution: "This code below is C code, obtained from the common lisp
// stat project under the BSD licence:"
// https://raw.githubusercontent.com/blindglobe/common-lisp-stat/3bdd28c4ae3de28dce32d8b9158c1f8d1b2e3924/lib/lowess.c
// "Like much lowess code, it is derived from the initial FORTRAN code by W. S.
// Cleveland published at NETLIB." http://www.netlib.org/go/lowess.f

//! The LOWESS transformation model, `TransformationModelLowess`, and the
//! smoother it runs, `FastLowessSmoothing::lowess`.
//!
//! [`LowessModel`] sorts the anchors, smooths them with [`lowess`] and hands
//! the smoothed points to an [`InterpolatedModel`]. [`LowessModel::select_span`]
//! ports the model's optional automatic span selection by cross-validation
//! (`auto_span`), with the source's private helpers `buildSpanGrid` and
//! `scoreResiduals` as [`build_span_grid`] and [`score_residuals`]. See
//! `docs/TRANSFORMATIONS_SUPPORT.md`.

use super::{
    DataPoint, Extrapolation, InterpolatedModel, Interpolation, InterpolationOptions, bad, finite,
    validate_points,
};
use crate::Result;
use crate::math::source_sort::source_sort_by;
use crate::math::statistic_functions::{mean_absolute_deviation, quantile, root_mean_square_error};

/// Options of the LOWESS model (source `getDefaultParameters`, without the
/// automatic span selection, which is [`LowessAutoSpan`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LowessOptions {
    /// Fraction of observations per local fit, strictly positive and at most
    /// one (source `span`, default 2/3). It determines the amount of smoothing;
    /// the source recommends 0.2 to 0.8. The source states `0 < span <= 1` as
    /// a precondition of `FastLowessSmoothing::lowess` and checks it only in
    /// debug builds; it is checked here.
    pub span: f64,
    /// Robustifying passes after the initial fit (source `num_iterations`,
    /// default 3, minimum 0, no maximum). The cost is bounded by
    /// [`LowessOptions::max_work`], which every pass charges against, so an
    /// absurd count fails with the work-limit error instead of a native cap
    /// on the count (an earlier cap of 64 refused values the source accepts;
    /// the Release build's grouping with `num_iterations` 70 is reproduced
    /// by `tests/mapmatching_kd_grouping_oracle.rs`).
    pub iterations: usize,
    /// Distance within which fits are interpolated instead of computed (source
    /// `delta`). `None` uses one percent of the x range, which the source
    /// selects with a negative value; `Some` requires a finite non-negative
    /// distance.
    pub delta: Option<f64>,
    /// Interpolation between the smoothed points (source `interpolation_type`,
    /// default `cspline`).
    pub interpolation: Interpolation,
    /// Extrapolation beyond them (source `extrapolation_type`, default
    /// `four-point-linear`).
    pub extrapolation: Extrapolation,
    /// Maximum number of observations (native limit).
    pub max_points: usize,
    /// Limit on local sample operations, interpolation, and residual passes
    /// (native limit).
    pub max_work: usize,
}
impl Default for LowessOptions {
    fn default() -> Self {
        Self {
            span: 2.0 / 3.0,
            iterations: 3,
            delta: None,
            interpolation: Interpolation::CubicSpline,
            extrapolation: Extrapolation::FourPointLinear,
            max_points: 1_000_000,
            max_work: 50_000_000,
        }
    }
}
impl LowessOptions {
    fn validate(self) -> Result<()> {
        if !self.span.is_finite()
            || self.span <= 0.0
            || self.span > 1.0
            || self.max_points < 2
            || self.max_work == 0
            || self.delta.is_some_and(|v| !v.is_finite() || v < 0.0)
        {
            return Err(bad(
                "invalid LOWESS span, iterations, delta or resource limits",
            ));
        }
        Ok(())
    }
}

/// Cross-validation metric of the automatic span selection (source
/// `auto_metric` and the private `CVMetric`), in the source's
/// `names_of_cvmetric` order.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum CrossValidationMetric {
    /// `rmse`: root mean square of the absolute errors.
    Rmse,
    /// `mae`: mean absolute error. The default.
    #[default]
    Mae,
    /// `p90`: 90th percentile of the absolute errors.
    P90,
    /// `p95`: 95th percentile.
    P95,
    /// `p99`: 99th percentile.
    P99,
}
impl CrossValidationMetric {
    /// Every metric, in the source's `names_of_cvmetric` order.
    pub const ALL: [Self; 5] = [Self::Rmse, Self::Mae, Self::P90, Self::P95, Self::P99];

    /// The source parameter value naming this metric.
    pub fn name(self) -> &'static str {
        match self {
            Self::Rmse => "rmse",
            Self::Mae => "mae",
            Self::P90 => "p90",
            Self::P95 => "p95",
            Self::P99 => "p99",
        }
    }

    /// The metric a source name selects.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] for any other name. The source
    /// looks the name up with `Helpers::indexOf`, which throws
    /// `Exception::ElementNotFound`, and it does so on every fit, whether or
    /// not the automatic selection is enabled.
    pub fn from_name(name: &str) -> Result<Self> {
        Self::ALL
            .into_iter()
            .find(|metric| metric.name() == name)
            .ok_or_else(|| bad(&format!("unknown LOWESS cross-validation metric '{name}'")))
    }
}

/// Automatic span selection by cross-validation (the source's `auto_span*`
/// parameters and `auto_metric`).
///
/// When it is enabled, the source fits the model once per candidate span and
/// fold, scores the held-out errors, and keeps the best candidate; see
/// [`LowessModel::select_span`].
#[derive(Clone, Debug, PartialEq)]
pub struct LowessAutoSpan {
    /// Lower bound for the selected span (source `auto_span_min`, default
    /// 0.15, minimum 0.001); raised to at least 0.01 by the grid.
    pub span_min: f64,
    /// Upper bound (source `auto_span_max`, default 0.80, maximum 0.99);
    /// lowered to at most 0.99 by the grid.
    pub span_max: f64,
    /// Minimum number of neighbours, `span * n`, of a candidate (source
    /// `auto_min_neighbors`, default 5, minimum 3).
    pub min_neighbors: i32,
    /// Folds when there are more than 50 points (source `auto_k_folds`,
    /// default 5, minimum 2); with 50 or fewer, leave-one-out is used.
    pub k_folds: i32,
    /// How the held-out errors are scored (source `auto_metric`).
    pub metric: CrossValidationMetric,
    /// Candidate spans in (0, 1] (source `auto_span_grid`, a comma-separated
    /// list; empty uses the default grid of [`build_span_grid`]).
    pub grid: Vec<f64>,
    /// Maximum number of cross-validation fits, candidates times folds
    /// (native limit, checked before the first fit).
    pub max_fits: usize,
}
impl Default for LowessAutoSpan {
    fn default() -> Self {
        Self {
            span_min: 0.15,
            span_max: 0.80,
            min_neighbors: 5,
            k_folds: 5,
            metric: CrossValidationMetric::Mae,
            grid: Vec::new(),
            max_fits: 10_000,
        }
    }
}

/// Absolute tolerance within which two cross-validation scores tie (source
/// `kTieTol`); on a tie the larger span wins.
pub const SPAN_SCORE_TIE_TOLERANCE: f64 = 1e-12;

/// The default candidate spans of [`build_span_grid`].
pub const DEFAULT_SPAN_GRID: [f64; 13] = [
    0.01, 0.05, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90,
];

// The raise-then-lower in the body is deliberately not `f64::clamp`, which
// panics when the neighbour bound exceeds the maximum; the source ends on the
// maximum then.
#[allow(clippy::manual_clamp)]
/// The candidate span grid of the automatic selection (source private
/// `buildSpanGrid`).
///
/// The rules, in the source's order:
///
/// - an empty `candidates` uses [`DEFAULT_SPAN_GRID`];
/// - the bounds become `max(0.01, span_min)` and `min(0.99, span_max)`, and
///   are swapped if they cross;
/// - the lower bound is raised to `min(0.9, min_neighbors / n_points)` (1.0
///   when `n_points` is zero) so that no local fit has only one or two
///   neighbours;
/// - each candidate is clamped to the bounds, then the grid is sorted and
///   values closer than `1e-9` to the previously kept one are dropped;
/// - an empty result becomes `min(0.95, max(0.01, lower))`.
///
/// # Errors
///
/// Returns [`crate::Error::InvalidValue`] when a candidate or bound is not
/// finite. The source would carry a NaN through the clamps (every comparison
/// with it is false) into `std::sort`, where it is not a strict weak order.
pub fn build_span_grid(
    n_points: usize,
    candidates: &[f64],
    span_min: f64,
    span_max: f64,
    min_neighbors: i32,
) -> Result<Vec<f64>> {
    if !span_min.is_finite() || !span_max.is_finite() || candidates.iter().any(|v| !v.is_finite()) {
        return Err(bad("LOWESS span grid values must be finite"));
    }
    let mut grid = if candidates.is_empty() {
        DEFAULT_SPAN_GRID.to_vec()
    } else {
        candidates.to_vec()
    };
    let mut low = span_min.max(0.01);
    let mut high = span_max.min(0.99);
    if low > high {
        std::mem::swap(&mut low, &mut high);
    }
    let neighbors = if n_points > 0 {
        f64::from(min_neighbors) / n_points as f64
    } else {
        1.0
    };
    let lower = low.max(neighbors.min(0.9));
    for value in &mut grid {
        // Raise, then lower, in that order: when `lower` exceeds `high` (a
        // neighbour bound above the maximum) the source ends on `high`, which
        // `f64::clamp` would refuse by panicking.
        *value = value.max(lower).min(high);
    }
    grid.sort_by(f64::total_cmp);
    grid.dedup_by(|value, kept| (*value - *kept).abs() < 1e-9);
    if grid.is_empty() {
        grid.push(lower.max(0.01).min(0.95));
    }
    Ok(grid)
}

/// Score absolute held-out errors; lower is better (source private
/// `scoreResiduals`).
///
/// `rmse` is the root mean square of the errors against zero, `mae` their
/// mean, and `p90`/`p95`/`p99` the type-7 interpolated quantile of the sorted
/// errors. An empty set scores `+inf`, which the source uses to mark an
/// unusable candidate.
///
/// # Errors
///
/// Returns [`crate::Error::InvalidValue`] for a NaN error; the source only
/// ever passes finite errors, which it filters before scoring.
pub fn score_residuals(errors: &[f64], metric: CrossValidationMetric) -> Result<f64> {
    if errors.is_empty() {
        return Ok(f64::INFINITY);
    }
    if errors.iter().any(|v| v.is_nan()) {
        return Err(bad("LOWESS cross-validation errors must not be NaN"));
    }
    let percentile = |q: f64| -> Result<f64> {
        let mut sorted = errors.to_vec();
        sorted.sort_by(f64::total_cmp);
        quantile(&sorted, q)
    };
    match metric {
        CrossValidationMetric::Rmse => root_mean_square_error(errors, &vec![0.0; errors.len()]),
        CrossValidationMetric::Mae => Ok(mean_absolute_deviation(errors, 0.0)),
        CrossValidationMetric::P90 => percentile(0.90),
        CrossValidationMetric::P95 => percentile(0.95),
        CrossValidationMetric::P99 => percentile(0.99),
    }
}

/// A LOWESS smoothing model with interpolation between the smoothed points
/// (source `TransformationModelLowess`).
#[derive(Clone, Debug)]
pub struct LowessModel {
    interpolation: InterpolatedModel,
    options: LowessOptions,
}
impl LowessModel {
    /// Sort the anchors by x, smooth them, and interpolate the smoothed points
    /// (the source constructor without automatic span selection).
    ///
    /// The anchors are sorted with the permutation the C++ Release build's
    /// `std::sort` produces under the source's x-only comparison, so anchors
    /// with equal x reach the smoother in the order the source gives them.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] for invalid options, fewer than
    /// two anchors (the source's `Exception::IllegalArgument`), fewer than
    /// three distinct x values after smoothing (raised by the interpolation,
    /// as in the source), a work limit that is exceeded, or non-finite input
    /// or arithmetic.
    pub fn fit(data: &[DataPoint], options: LowessOptions) -> Result<Self> {
        options.validate()?;
        validate_points(data, options.max_points)?;
        if data.len() < 2 {
            return Err(bad("'lowess' model requires more data"));
        }
        let (x, y) = sorted_coordinates(data)?;
        let ys = lowess(&x, &y, options)?;
        let anchors: Vec<_> = x
            .iter()
            .zip(&ys)
            .map(|(&x, &y)| DataPoint::new(x, y))
            .collect();
        let interpolation = InterpolatedModel::fit(
            &anchors,
            InterpolationOptions {
                interpolation: options.interpolation,
                extrapolation: options.extrapolation,
                max_points: options.max_points,
            },
        )?;
        Ok(Self {
            interpolation,
            options,
        })
    }

    /// Select the span by cross-validation, as the source constructor does
    /// when `auto_span` is `true`; returns the selected span and its score.
    ///
    /// The anchors are sorted as in [`LowessModel::fit`]. With 50 or fewer
    /// anchors every anchor is its own fold (leave-one-out); otherwise there
    /// are `max(2, k_folds)` folds, assigned round-robin by sorted position
    /// (`CrossValidation::makeKFolds`). For each candidate of
    /// [`build_span_grid`], every fold whose complement has at least
    /// `max(3, min_neighbors)` anchors is held out in turn: a model with this
    /// candidate span and otherwise the given `options` is fitted to the rest,
    /// and the absolute errors at the held-out anchors are scored with
    /// [`score_residuals`]. A candidate with `span * n < min_neighbors` gets no
    /// errors and scores `+inf`. The first candidate starts as the best; a
    /// later one replaces it when its score is lower by more than
    /// [`SPAN_SCORE_TIE_TOLERANCE`], or ties within it and is larger
    /// (`CrossValidation::gridSearch1D` with `PreferLarger`).
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] for invalid options or grid
    /// values, fewer than two anchors, more than
    /// [`LowessAutoSpan::max_fits`] fits, and any error of a cross-validation
    /// fit: the source lets an exception from a fold's model (too few distinct
    /// x values in a training set, for instance) escape the whole selection,
    /// and so does this.
    pub fn select_span(
        data: &[DataPoint],
        options: LowessOptions,
        auto: &LowessAutoSpan,
    ) -> Result<(f64, f64)> {
        validate_points(data, options.max_points)?;
        if data.len() < 2 {
            return Err(bad("'lowess' model requires more data"));
        }
        let (x, y) = sorted_coordinates(data)?;
        let n = x.len();
        let folds = if n <= 50 {
            n
        } else {
            usize::try_from(auto.k_folds.max(2)).unwrap_or(2).min(n)
        };
        let spans = build_span_grid(
            n,
            &auto.grid,
            auto.span_min,
            auto.span_max,
            auto.min_neighbors,
        )?;
        let fits = spans
            .len()
            .checked_mul(folds)
            .ok_or_else(|| bad("LOWESS cross-validation fit count overflows"))?;
        if fits > auto.max_fits {
            return Err(bad("LOWESS cross-validation exceeds its fit limit"));
        }
        let minimum_training = usize::try_from(auto.min_neighbors.max(3)).unwrap_or(3);
        let mut best: Option<(f64, f64)> = None;
        for &span in &spans {
            let mut errors = Vec::new();
            if span * n as f64 >= f64::from(auto.min_neighbors) {
                let fold_options = LowessOptions { span, ..options };
                for fold in 0..folds {
                    // Round-robin fold membership: index i belongs to fold i % K.
                    let training: Vec<_> = (0..n)
                        .filter(|i| i % folds != fold)
                        .map(|i| DataPoint::new(x[i], y[i]))
                        .collect();
                    if training.len() < minimum_training {
                        continue;
                    }
                    let model = Self::fit(&training, fold_options)?;
                    for i in (fold..n).step_by(folds) {
                        // The source drops a non-finite error; `apply` reports
                        // a non-finite prediction as an error, which is dropped
                        // the same way.
                        if let Ok(prediction) = model.apply(x[i]) {
                            let error = (prediction - y[i]).abs();
                            if error.is_finite() {
                                errors.push(error);
                            }
                        }
                    }
                }
            }
            let score = score_residuals(&errors, auto.metric)?;
            let replace = match best {
                None => true,
                Some((best_span, best_score)) => {
                    score < best_score - SPAN_SCORE_TIE_TOLERANCE
                        || ((score - best_score).abs() <= SPAN_SCORE_TIE_TOLERANCE
                            && span > best_span)
                }
            };
            if replace {
                best = Some((span, score));
            }
        }
        // build_span_grid never returns an empty grid.
        best.ok_or_else(|| bad("LOWESS span grid is empty"))
    }

    /// Select the span by cross-validation, then fit with it (the source
    /// constructor with `auto_span` set to `true`); returns the model and the
    /// selected span's score. The model's [`LowessModel::options`] carry the
    /// selected span, which is what the source persists into its parameters
    /// (together with `auto_span` reset to `false`) so that a refit does not
    /// select again.
    ///
    /// # Errors
    ///
    /// As [`LowessModel::select_span`] and [`LowessModel::fit`].
    pub fn fit_auto(
        data: &[DataPoint],
        options: LowessOptions,
        auto: &LowessAutoSpan,
    ) -> Result<(Self, f64)> {
        let (span, score) = Self::select_span(data, options, auto)?;
        Ok((Self::fit(data, LowessOptions { span, ..options })?, score))
    }

    /// The options the model was fitted with.
    pub fn options(&self) -> LowessOptions {
        self.options
    }
    /// The smoothed points the interpolation runs through, sorted by x with
    /// repeated x averaged.
    pub fn knots(&self) -> impl Iterator<Item = (f64, f64)> + '_ {
        self.interpolation.knots()
    }
    /// Evaluate the model (source `evaluate`, which delegates to the
    /// interpolation).
    ///
    /// # Errors
    ///
    /// As [`InterpolatedModel::apply`].
    pub fn apply(&self, x: f64) -> Result<f64> {
        self.interpolation.apply(x)
    }
}

/// The anchors' coordinates in the order `std::sort` with the source's
/// x-only comparator (`cmpFirstDimension`) leaves them.
fn sorted_coordinates(data: &[DataPoint]) -> Result<(Vec<f64>, Vec<f64>)> {
    let mut sorted: Vec<(f64, f64)> = data.iter().map(|p| (p.x, p.y)).collect();
    source_sort_by(&mut sorted, |a, b| a.0 < b.0)?;
    Ok(sorted.into_iter().unzip())
}

/// Robust locally weighted linear regression (source
/// `FastLowessSmoothing::lowess`), returning the smoothed y values.
///
/// Tricube distance weights, the local linear correction, repeated-x
/// handling, robustness reweighting with the source's residual
/// pseudo-median, and linear interpolation over delta-skipped fits.
/// Inputs must be finite and sorted by non-decreasing x; signed y values are
/// valid. At least two observations are needed; the interpolated model needs
/// three distinct x coordinates after smoothing.
///
/// # Errors
///
/// Returns [`crate::Error::InvalidValue`] for invalid options, arrays of
/// different length, fewer than two or more than `max_points` observations,
/// non-finite input or arithmetic, or an exceeded work limit, and
/// [`crate::Error::UnsortedData`] when x is not sorted. The source states the
/// sortedness, length and parameter conditions as preconditions checked only
/// in debug builds.
pub fn lowess(x: &[f64], y: &[f64], options: LowessOptions) -> Result<Vec<f64>> {
    options.validate()?;
    if x.len() != y.len() || x.len() < 2 || x.len() > options.max_points {
        return Err(bad(
            "LOWESS requires equal arrays with 2..=max_points observations",
        ));
    }
    for &value in x.iter().chain(y) {
        finite(value)?;
    }
    if x.windows(2).any(|p| p[0] > p[1]) {
        return Err(crate::Error::UnsortedData);
    }
    let range = finite(x[x.len() - 1] - x[0])?;
    let delta = options.delta.unwrap_or(range * 0.01);
    let n = x.len();
    let neighbors = ((options.span * n as f64) as usize).clamp(2, n);
    let mut ys = vec![0.0; n];
    let mut weights = vec![0.0; n];
    let mut residual_weights = vec![0.0; n];
    let mut work = 0usize;
    let mut charge = |amount: usize| -> Result<()> {
        work = work
            .checked_add(amount)
            .ok_or_else(|| bad("LOWESS work count overflows"))?;
        if work > options.max_work {
            Err(bad("LOWESS exceeds configured work limit"))
        } else {
            Ok(())
        }
    };
    for iteration in 0..=options.iterations {
        let (mut left, mut right, mut i) = (0, neighbors - 1, 0);
        let mut last: Option<usize> = None;
        loop {
            // Shift the fixed-size neighborhood only when its radius decreases.
            while right < n - 1 && x[i] - x[left] > x[right + 1] - x[i] {
                charge(1)?;
                left += 1;
                right += 1;
            }
            let h = finite((x[i] - x[left]).max(x[right] - x[i]))?;
            let mut sum = 0.0;
            let mut end = left;
            for j in left..n {
                charge(1)?;
                weights[j] = 0.0;
                let distance = (x[j] - x[i]).abs();
                if distance <= 0.999 * h {
                    weights[j] = if distance > 0.001 * h {
                        let r = distance / h;
                        let c = 1.0 - r * r * r;
                        c * c * c
                    } else {
                        1.0
                    };
                    if iteration > 0 {
                        weights[j] *= residual_weights[j];
                    }
                    sum += weights[j];
                } else if x[j] > x[i] {
                    break;
                }
                end = j + 1;
            }
            finite(sum)?;
            if sum <= 0.0 {
                ys[i] = y[i];
            } else {
                charge(
                    (end - left)
                        .checked_mul(5)
                        .ok_or_else(|| bad("LOWESS work count overflows"))?,
                )?;
                for weight in &mut weights[left..end] {
                    *weight /= sum;
                }
                if h > 0.0 {
                    let center = finite((left..end).map(|j| weights[j] * x[j]).sum())?;
                    let spread = finite(
                        (left..end)
                            .map(|j| weights[j] * (x[j] - center) * (x[j] - center))
                            .sum(),
                    )?;
                    if spread.sqrt() > 0.001 * range {
                        let b = finite((x[i] - center) / spread)?;
                        for j in left..end {
                            weights[j] = finite(weights[j] * (1.0 + b * (x[j] - center)))?;
                        }
                    }
                }
                ys[i] = finite((left..end).map(|j| weights[j] * y[j]).sum())?;
            }
            if let Some(previous) = last {
                if previous + 1 < i {
                    charge(i - previous - 1)?;
                    let denominator = x[i] - x[previous];
                    for j in previous + 1..i {
                        let alpha = (x[j] - x[previous]) / denominator;
                        ys[j] = finite(alpha * ys[i] + (1.0 - alpha) * ys[previous])?;
                    }
                }
            }
            let mut fitted = i;
            let cut = finite(x[i] + delta)?;
            let mut next = i + 1;
            while next < n && x[next] <= cut {
                charge(1)?;
                if x[next] == x[fitted] {
                    ys[next] = ys[fitted];
                    fitted = next;
                }
                next += 1;
            }
            last = Some(fitted);
            if fitted == n - 1 {
                break;
            }
            i = (fitted + 1).max(next - 1);
        }
        charge(n)?;
        for i in 0..n {
            weights[i] = finite(y[i] - ys[i])?;
        }
        if iteration == options.iterations {
            break;
        }
        charge(n)?;
        let mut sorted: Vec<_> = weights.iter().map(|v| v.abs()).collect();
        sorted.sort_by(f64::total_cmp);
        // Source pseudo-median uses middle and preceding residual even for odd n.
        let cmad = finite(3.0 * (sorted[n / 2] + sorted[n / 2 - 1]))?;
        for i in 0..n {
            let r = weights[i].abs();
            residual_weights[i] = if r <= 0.001 * cmad {
                1.0
            } else if r > 0.999 * cmad {
                0.0
            } else {
                let ratio = r / cmad;
                let q = 1.0 - ratio * ratio;
                q * q
            };
        }
    }
    Ok(ys)
}
