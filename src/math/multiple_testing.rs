// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Multiple-testing correction: q-values, pi0 estimation, local FDR and
//! empirical p-values, after Storey and Tibshirani.
//!
//! Port of `src/openms/include/OpenMS/MATH/STATISTICS/MultipleTesting.h` and
//! its translation unit. See `docs/MULTIPLE_TESTING_SUPPORT.md`.
//!
//! References:
//!
//! > J. D. Storey (2002). A direct approach to false discovery rates.
//! > J. R. Statist. Soc. B 64(3):479-498.
//!
//! > J. D. Storey and R. Tibshirani (2003). Statistical significance for
//! > genome-wide experiments. PNAS 100:9440-9445.
//!
//! The source tracks PyProphet's Python implementation of the R `qvalue`
//! package rather than the papers directly; where the two differ the port
//! follows the source and says so at the item.
//!
//! Ranking is [`crate::math::rank_data::rankdata`], the kernel density estimate
//! behind [`crate::math::multiple_testing::lfdr`] is
//! [`crate::math::kernel_density::kde_fft_eval`], and the smoothing spline that
//! [`crate::math::multiple_testing::pi0_est`] needs is supplied by the caller
//! through [`crate::math::multiple_testing::Pi0Smoother`] — see that trait for
//! why, and for why passing `None` there is a different answer rather than a
//! safe one: it lowers `pi0`, and a lower `pi0` lowers every q-value and local
//! FDR derived from it.
//!
//! # Non-finite p-values
//!
//! [`crate::math::multiple_testing::q_value`] and
//! [`crate::math::multiple_testing::lfdr`] drop non-finite entries, compute on
//! what is left and write `NaN` back into the dropped positions.
//! [`crate::math::multiple_testing::compute_model_fdr`] does the opposite and
//! returns an all-`NaN` vector if *any* entry is `NaN`. Both behaviours are the
//! source's, and the header documents the asymmetry as a warning to callers.
//!
//! # Differences from the source
//!
//! * The two `std::invalid_argument` throws become [`crate::Error::InvalidValue`], and
//!   the string-keyed free functions `pi0Est(..., "smoother", ...)` and
//!   `lfdr(..., "probit", ...)` become the enum-typed calls plus
//!   [`crate::math::multiple_testing::Pi0Method::parse`] and
//!   [`crate::math::multiple_testing::LfdrTransform::parse`].
//! * `lfdr`'s six trailing parameters become
//!   [`crate::math::multiple_testing::LfdrOptions`], whose `Default` is the
//!   source's default argument list.
//! * Serial, as the source is: `MultipleTesting.cpp` carries no `#pragma omp`.

use crate::math::kernel_density::{DEFAULT_CUT, DEFAULT_GRIDSIZE, bw_nrd0, kde_fft_eval};
use crate::math::rank_data::{NanPolicy, RankMethod, rankdata_f64};
use crate::{Error, Result};
use std::f64::consts::PI;

/// Maximum number of p-values or statistics one call may consume.
///
/// Native guard; the source allocates whatever it is handed.
pub const MAX_ITEMS: usize = 50_000_000;

/// Default bandwidth inflation for [`lfdr`], `1.5`.
pub const DEFAULT_LFDR_ADJ: f64 = 1.5;

/// Default clipping constant for [`lfdr`], `1e-8`.
pub const DEFAULT_LFDR_EPS: f64 = 1e-8;

/// Default smoothing-spline degrees of freedom for [`pi0_est`], `3`.
pub const DEFAULT_SMOOTH_DF: i32 = 3;

/// `!(a > b)`: true when `a <= b` **and** when either value is `NaN`.
///
/// Spelled through `partial_cmp` because clippy rejects the negated comparison
/// on a partially ordered type. The `NaN` arm is load-bearing at every call
/// site: the source writes `if (!(x > 0.0))` precisely so that a `NaN` takes
/// the guarded branch, and `x <= b` would not.
fn not_greater(a: f64, b: f64) -> bool {
    !matches!(a.partial_cmp(&b), Some(std::cmp::Ordering::Greater))
}

fn bad(message: String) -> Error {
    Error::InvalidValue(message)
}

fn check_items(len: usize, what: &str) -> Result<()> {
    if len > MAX_ITEMS {
        return Err(Error::InvalidRange(format!(
            "{what}: {len} values exceeds the maximum {MAX_ITEMS}"
        )));
    }
    Ok(())
}

/// Stable ascending argsort, the source's file-static `argsort_asc`.
///
/// `std::stable_sort` with `a[i] < a[j]`; equal values keep their input order.
/// The values here are always finite by the time this is called.
fn argsort_asc(values: &[f64]) -> Vec<usize> {
    let mut index: Vec<usize> = (0..values.len()).collect();
    index.sort_by(|i, j| values[*i].total_cmp(&values[*j]));
    index
}

/// Which estimator [`pi0_est`] uses for the proportion of true nulls.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum Pi0Method {
    /// Fit a smoothing spline through the per-lambda estimates and read it at
    /// the largest lambda. The source's default.
    #[default]
    Smoother,
    /// Pick the lambda whose bootstrap mean squared error is smallest.
    Bootstrap,
}

impl Pi0Method {
    /// The source's `pi0MethodToString`.
    ///
    /// The source has a `default:` arm returning `"unknown"` for an enum value
    /// outside the two declared ones; a Rust enum cannot hold such a value, so
    /// that arm has no counterpart.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Smoother => "smoother",
            Self::Bootstrap => "bootstrap",
        }
    }

    /// The source's `toPi0Method`, case-insensitive.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for anything but `smoother` or
    /// `bootstrap`, where the source throws `std::invalid_argument`.
    pub fn parse(text: &str) -> Result<Self> {
        match text.to_ascii_lowercase().as_str() {
            "smoother" => Ok(Self::Smoother),
            "bootstrap" => Ok(Self::Bootstrap),
            _ => Err(bad(format!(
                "toPi0Method: invalid method '{text}', expected 'smoother' or 'bootstrap'"
            ))),
        }
    }
}

/// Which transformation [`lfdr`] applies to the p-values before estimating
/// their density.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum LfdrTransform {
    /// Inverse standard normal CDF. The source's default; the null density is
    /// then the standard normal.
    #[default]
    Probit,
    /// Log-odds. The null density is then the logistic derivative.
    Logit,
}

impl LfdrTransform {
    /// The source's `lfdrTransformToString`.
    ///
    /// As [`Pi0Method::as_str`], the source's `"unknown"` arm is unreachable
    /// here.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Probit => "probit",
            Self::Logit => "logit",
        }
    }

    /// The source's `toLfdrTransform`, case-insensitive.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for anything but `probit` or `logit`,
    /// where the source throws `std::invalid_argument`.
    pub fn parse(text: &str) -> Result<Self> {
        match text.to_ascii_lowercase().as_str() {
            "probit" => Ok(Self::Probit),
            "logit" => Ok(Self::Logit),
            _ => Err(bad(format!(
                "toLfdrTransform: invalid transform '{text}', expected 'probit' or 'logit'"
            ))),
        }
    }
}

/// The trailing arguments of [`lfdr`], which the source spells as six default
/// parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LfdrOptions {
    /// Cap the estimates at `1.0`. The source's `trunc`, default `true`.
    pub truncate: bool,
    /// Force the estimates to be non-decreasing in the p-value. The source's
    /// `monotone`, default `true`.
    pub monotone: bool,
    /// Transformation applied before the density estimate.
    pub transform: LfdrTransform,
    /// Bandwidth inflation applied to [`bw_nrd0`]; the source's `adj`, default
    /// [`DEFAULT_LFDR_ADJ`]. A larger value smooths harder.
    pub adj: f64,
    /// Clipping constant that keeps a p-value of exactly `0` or `1` off the
    /// infinities of the transform; the source's `eps`, default
    /// [`DEFAULT_LFDR_EPS`].
    pub eps: f64,
    /// Grid size handed to the kernel density estimate; default
    /// [`DEFAULT_GRIDSIZE`].
    pub gridsize: usize,
    /// Grid extension handed to the kernel density estimate; default
    /// [`DEFAULT_CUT`].
    pub cut: f64,
}

impl Default for LfdrOptions {
    /// The source's default argument list:
    /// `trunc = true, monotone = true, transf = probit, adj = 1.5,
    /// eps = 1e-8, gridsize = 512, cut = 3.0`.
    fn default() -> Self {
        Self {
            truncate: true,
            monotone: true,
            transform: LfdrTransform::Probit,
            adj: DEFAULT_LFDR_ADJ,
            eps: DEFAULT_LFDR_EPS,
            gridsize: DEFAULT_GRIDSIZE,
            cut: DEFAULT_CUT,
        }
    }
}

/// The outcome of [`pi0_est`], the source's `Math::Pi0Result`.
#[derive(Clone, Debug, PartialEq)]
pub struct Pi0Result {
    /// Estimated proportion of true null hypotheses, in `[0, 1]`.
    ///
    /// `1.0` is the most conservative value this can take — it is what
    /// `Default` carries, not what the fallbacks produce. Every fallback path
    /// returns `min(min(pi0_lambda), 1)`, which is the **smallest** of the
    /// per-lambda estimates and therefore *less* conservative than the
    /// smoothed answer, not more. See [`Pi0Smoother`].
    pub pi0: f64,
    /// The per-lambda estimates, one per entry of [`Pi0Result::lambda`].
    ///
    /// These are **not** clamped to `[0, 1]` and routinely exceed one — except
    /// on the single-lambda path, where the source stores the clamped value
    /// instead. That inconsistency is the source's and is reproduced; see
    /// `docs/MULTIPLE_TESTING_SUPPORT.md`.
    pub pi0_lambda: Vec<f64>,
    /// The lambda thresholds actually used, which is the caller's grid or the
    /// default `0.05, 0.10, ..., 0.95`.
    pub lambda: Vec<f64>,
    /// Whether the smoothing spline was fitted and read successfully. `false`
    /// means [`Pi0Result::pi0`] is the `min(min(pi0_lambda), 1)` fallback, or
    /// that the bootstrap method was used.
    ///
    /// A caller that treats this flag as cosmetic gets a different, and
    /// systematically smaller, `pi0` than the source's default call; see
    /// [`Pi0Smoother`] for the measured size of the difference and what it
    /// does to downstream q-values.
    pub pi0_smooth: bool,
}

impl Default for Pi0Result {
    /// The source's member initialisers: `pi0 = 1.0`, empty vectors,
    /// `pi0_smooth = false`.
    fn default() -> Self {
        Self {
            pi0: 1.0,
            pi0_lambda: Vec::new(),
            lambda: Vec::new(),
            pi0_smooth: false,
        }
    }
}

/// The smoothing spline [`pi0_est`] reads at the largest lambda.
///
/// The source calls `OpenMS::Math::BSplineSmoothingSpline spl(xs, ys, -1.0,
/// smooth_df)` and then `spl.eval(max_lambda)`, guarding on `spl.ok()`. That
/// class is ported, as `BSplineSmoothingSpline` in
/// `src/processing/spline/smoothing.rs`, but it lives in the `processing`
/// module, and this crate ratchets its cross-module dependency graph: `math`
/// reaches no other top-level module today and may not start. Rather than
/// duplicate a 1,800-line spline implementation, the requirement is stated as
/// this trait and satisfied above the `math` layer.
///
/// `tests/multiple_testing.rs` implements it with the ported
/// `BSplineSmoothingSpline` and reproduces the class test's pi0 literals, so
/// the seam is exercised rather than merely declared.
///
/// # Passing `None` is not a safe default
///
/// In the source the `!spl.ok()` fallback is a pathology; here it is whatever
/// a caller who omits the smoother gets, and it returns a **different number**.
/// Measured on `tests/data/test_lfdr_ref_data.csv`, the 3,170 PyProphet
/// p-values the class test ships, with the default lambda grid:
///
/// | call | `pi0` |
/// | --- | --- |
/// | `pi0_est(p, &[], Smoother, 3, false, Some(spline))` | `0.6685639` |
/// | `pi0_est(p, &[], Smoother, 3, false, None)` | `0.6403785` |
///
/// The first matches the C++ default call's literal `0.6685638`. The second is
/// `min(min(pi0_lambda), 1)`, and it is **lower**. Since
/// [`q_value`] computes `q = pi0 * m * p / rank(p)` and [`lfdr`] computes
/// `pi0 * f0 / y`, both scale linearly in `pi0`: a lower `pi0` makes every
/// q-value and every local FDR **smaller**, so more hypotheses clear any fixed
/// threshold. Omitting the smoother therefore loosens the multiple-testing
/// correction rather than tightening it, and a caller who wants the source's
/// behaviour must supply the spline. [`Pi0Result::pi0_smooth`] reports which
/// path was taken.
pub trait Pi0Smoother {
    /// Fit a smoothing spline through `(x, y)` and evaluate it at `at`.
    ///
    /// `x` is strictly increasing and has at least two entries; `y` has the
    /// same length. `smooth_df` is the spline degree the source passes as the
    /// fourth constructor argument, and the smoothing parameter it passes is
    /// `-1.0`, meaning "choose automatically".
    ///
    /// Return `None` where the source's `!spl.ok()` guard fires, so that
    /// [`pi0_est`] falls back to `min(pi0_lambda, 1)`.
    fn smooth_eval(&self, x: &[f64], y: &[f64], smooth_df: i32, at: f64) -> Option<f64>;
}

/// The source's file-static `percentile`: nearest-rank on a sorted copy.
///
/// `idx = floor(p * (n - 1) + 0.5)`, which rounds half away from zero rather
/// than interpolating, so it is neither `numpy.percentile` nor
/// [`crate::math::statistic_functions::quantile`]. Only the bootstrap branch of
/// [`pi0_est`] uses it.
fn percentile(values: &[f64], p: f64) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let position = p * (sorted.len() as f64 - 1.0);
    let index_f = (position + 0.5).floor();
    let index = if index_f >= sorted.len() as f64 || index_f.is_nan() || index_f < 0.0 {
        sorted.len() - 1
    } else {
        index_f as usize
    };
    sorted[index]
}

/// q-values, the minimum FDR at which each test is called significant.
///
/// For each finite p-value, `q = pi0 * m * p / rank(p)` with the `max` tie rank
/// over the `m` finite p-values, then a right-to-left cumulative minimum over
/// the p-ordering so that q is non-decreasing in p. The largest q is capped at
/// `1`.
///
/// # Arguments
///
/// * `pi0` — proportion of true nulls, in `[0, 1]`; typically from
///   [`pi0_est`].
/// * `pfdr` — compute the *positive* FDR, dividing additionally by
///   `1 - (1 - p)^m`, which conditions on at least one rejection.
///
/// Non-finite inputs are dropped and their positions come back `NaN`; an empty
/// input gives an empty result. A rank of zero, which cannot occur for
/// one-based ranks, would give `+inf` rather than a division by zero — the
/// source guards it and so does this.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when a finite p-value is outside `[0, 1]` or
/// `pi0` is outside `[0, 1]`, both `std::invalid_argument` in the source, and
/// [`Error::InvalidRange`] when the input exceeds [`MAX_ITEMS`].
pub fn q_value(p_values: &[f64], pi0: f64, pfdr: bool) -> Result<Vec<f64>> {
    check_items(p_values.len(), "qValue")?;
    let total = p_values.len();
    if total == 0 {
        return Ok(Vec::new());
    }
    let keep: Vec<bool> = p_values.iter().map(|v| v.is_finite()).collect();
    let p: Vec<f64> = p_values.iter().copied().filter(|v| v.is_finite()).collect();

    let mut out = vec![f64::NAN; total];
    let m = p.len();
    if m == 0 {
        return Ok(out);
    }
    if p.iter().any(|v| *v < 0.0 || *v > 1.0) {
        return Err(bad("qValue: p-values not in [0,1]".to_string()));
    }
    if !(0.0..=1.0).contains(&pi0) {
        return Err(bad("qValue: pi0 not in [0,1]".to_string()));
    }

    let order = argsort_asc(&p);
    let ranks = rankdata_f64(&p, RankMethod::Max, NanPolicy::Propagate)?;

    let mut q = vec![0.0; m];
    let m_f = m as f64;
    for i in 0..m {
        let denom = if pfdr {
            ranks[i] * (1.0 - (1.0 - p[i]).powf(m_f))
        } else {
            ranks[i]
        };
        q[i] = if denom == 0.0 {
            f64::INFINITY
        } else {
            (pi0 * m_f * p[i]) / denom
        };
    }

    // Monotonicity: cap the largest p-value's q at 1, then sweep down the
    // p-ordering taking a running minimum. The NaN guard is the source's; a
    // `pfdr` denominator of `0 * inf` is the only way to reach it.
    let last = order[m - 1];
    q[last] = q[last].min(1.0);
    if q[last].is_nan() {
        q[last] = 1.0;
    }
    for ii in (0..m - 1).rev() {
        let index = order[ii];
        let next = order[ii + 1];
        q[index] = q[index].min(q[next]);
    }

    let mut k = 0;
    for (i, slot) in out.iter_mut().enumerate() {
        if keep[i] {
            *slot = q[k];
            k += 1;
        }
    }
    Ok(out)
}

/// Estimate the proportion of true null hypotheses.
///
/// For each lambda the estimate is `#{p >= lambda} / (m (1 - lambda))`. A
/// single lambda returns that estimate capped at `1`. With several, the method
/// decides:
///
/// * [`Pi0Method::Smoother`] fits a smoothing spline through the per-lambda
///   estimates and reads it at the largest lambda, clamped to `[0, 1]`. Fewer
///   than four lambdas, fewer than two distinct ones, no `smoother`, a spline
///   that does not fit, or a `NaN` prediction all fall back to the smallest
///   per-lambda estimate, `min(min(pi0_lambda), 1)` — which is lower than the
///   smoothed value, not a conservative cap on it.
/// * [`Pi0Method::Bootstrap`] picks the lambda minimising
///   `W / (m^2 (1-lambda)^2) * (1 - W/m) + (pi0_lambda - minpi0)^2`, with
///   `minpi0` the 10th percentile of the per-lambda estimates - by the source's
///   own nearest-rank rule, `floor(0.1 (n - 1) + 0.5)`, which neither
///   interpolates nor matches
///   [`crate::math::statistic_functions::quantile`] - and `W` the
///   count at or above that lambda.
///
/// # Arguments
///
/// * `lambda` — thresholds in `[0, 1)`. An empty slice means the source's
///   default `0.05` to `0.95`, generated by the same repeated `+= 0.05` so the
///   grid carries the same accumulated rounding (`0.15000000000000002`, not
///   `0.15`).
/// * `smooth_df` — spline degree handed to `smoother`; the source's default is
///   [`DEFAULT_SMOOTH_DF`].
/// * `smooth_log_pi0` — fit the spline to `log(pi0_lambda)` and exponentiate
///   the prediction. A non-positive estimate becomes `-inf` before the fit, as
///   in the source. With this set, a non-finite prediction is *not* a fallback
///   trigger, because `exp` of `-inf` is a legitimate `0`.
/// * `smoother` — see [`Pi0Smoother`]. `None` behaves as a spline that refused
///   to fit, which returns `min(min(pi0_lambda), 1)` instead of the smoothed
///   estimate. That is a *lower* `pi0` than the source's default call — on the
///   class test's own fixture, `0.6403785` against `0.6685638` — and a lower
///   `pi0` makes downstream q-values and local FDRs smaller, not larger. It is
///   a different answer, not a safe one.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when no finite p-value is supplied, when a
/// p-value is outside `[0, 1]`, or when a lambda is outside `[0, 1)` — the
/// three `std::invalid_argument` throws — and [`Error::InvalidRange`] when the
/// input exceeds [`MAX_ITEMS`].
pub fn pi0_est(
    p_values: &[f64],
    lambda: &[f64],
    method: Pi0Method,
    smooth_df: i32,
    smooth_log_pi0: bool,
    smoother: Option<&dyn Pi0Smoother>,
) -> Result<Pi0Result> {
    check_items(p_values.len(), "pi0Est")?;
    let p: Vec<f64> = p_values.iter().copied().filter(|v| v.is_finite()).collect();
    let m = p.len();
    if m == 0 {
        return Err(bad("pi0Est: no finite p-values provided".to_string()));
    }

    let lambda_v: Vec<f64> = if lambda.is_empty() {
        let mut grid = Vec::new();
        // The source accumulates `l += 0.05`, so the grid is not exactly
        // `0.05 * k`. Reproduced literally, because the spline is fitted
        // through these abscissae.
        let mut l = 0.05;
        while l < 1.0 - 1e-12 {
            grid.push(l);
            l += 0.05;
        }
        grid
    } else {
        lambda.to_vec()
    };
    let ll = lambda_v.len();
    if ll == 0 {
        return Err(bad("pi0est: empty lambda".to_string()));
    }
    if p.iter().any(|v| *v < 0.0 || *v > 1.0) {
        return Err(bad("pi0Est: p-values not in [0,1]".to_string()));
    }
    if lambda_v.iter().any(|l| *l < 0.0 || *l >= 1.0) {
        return Err(bad("pi0Est: lambda must be in [0,1)".to_string()));
    }

    let m_f = m as f64;
    if ll == 1 {
        let l = lambda_v[0];
        let mut frac = 0.0;
        for value in &p {
            if *value >= l {
                frac += 1.0;
            }
        }
        frac /= m_f;
        let pi0 = (frac / (1.0 - l)).min(1.0);
        return Ok(Pi0Result {
            pi0,
            pi0_lambda: vec![pi0],
            lambda: lambda_v,
            pi0_smooth: false,
        });
    }

    let mut pi0s = Vec::with_capacity(ll);
    for l in &lambda_v {
        let mut count = 0.0;
        for value in &p {
            if *value >= *l {
                count += 1.0;
            }
        }
        pi0s.push((count / m_f) / (1.0 - *l));
    }

    // The source's `min(*min_element(pi0s), 1.0)` fallback. Named for what it
    // computes, not for a risk direction: it is the *smallest* per-lambda
    // estimate, so it yields a smaller pi0 than the smoother and hence smaller
    // q-values downstream. See `Pi0Smoother`.
    let min_fallback = |pi0s: &[f64], lambda_v: Vec<f64>, pi0_lambda: Vec<f64>| Pi0Result {
        pi0: pi0s.iter().copied().fold(f64::INFINITY, f64::min).min(1.0),
        pi0_lambda,
        lambda: lambda_v,
        pi0_smooth: false,
    };

    match method {
        Pi0Method::Smoother => {
            if ll < 4 {
                return Ok(min_fallback(&pi0s, lambda_v, pi0s.clone()));
            }
            let mut y = pi0s.clone();
            if smooth_log_pi0 {
                for value in y.iter_mut() {
                    *value = if *value > 0.0 {
                        value.ln()
                    } else {
                        f64::NEG_INFINITY
                    };
                }
            }
            // The source builds a `std::map<double, double>`, whose
            // `operator[]` assignment sorts by lambda and lets a repeated
            // lambda's *last* value win.
            let mut xy: Vec<(f64, f64)> = Vec::with_capacity(ll);
            for i in 0..ll {
                match xy.binary_search_by(|probe| probe.0.total_cmp(&lambda_v[i])) {
                    Ok(at) => xy[at].1 = y[i],
                    Err(at) => xy.insert(at, (lambda_v[i], y[i])),
                }
            }
            if xy.len() < 2 {
                return Ok(min_fallback(&pi0s, lambda_v, pi0s.clone()));
            }
            let xs: Vec<f64> = xy.iter().map(|pair| pair.0).collect();
            let ys: Vec<f64> = xy.iter().map(|pair| pair.1).collect();
            let max_lambda = lambda_v.iter().copied().fold(f64::NEG_INFINITY, f64::max);

            let predicted = smoother.and_then(|s| s.smooth_eval(&xs, &ys, smooth_df, max_lambda));
            let Some(mut pred) = predicted else {
                return Ok(min_fallback(&pi0s, lambda_v, pi0s.clone()));
            };
            if pred.is_nan() || (!smooth_log_pi0 && !pred.is_finite()) {
                return Ok(min_fallback(&pi0s, lambda_v, pi0s.clone()));
            }
            if smooth_log_pi0 {
                pred = pred.exp();
            }
            if pred.is_nan() {
                return Ok(min_fallback(&pi0s, lambda_v, pi0s.clone()));
            }
            Ok(Pi0Result {
                pi0: pred.clamp(0.0, 1.0),
                pi0_lambda: pi0s,
                lambda: lambda_v,
                pi0_smooth: true,
            })
        }
        Pi0Method::Bootstrap => {
            let minpi0 = percentile(&pi0s, 0.1);
            let mut counts = Vec::with_capacity(ll);
            for l in &lambda_v {
                let mut count = 0.0;
                for value in &p {
                    if *value >= *l {
                        count += 1.0;
                    }
                }
                counts.push(count);
            }
            let mut mse = vec![0.0; ll];
            for i in 0..ll {
                let l = lambda_v[i];
                let w = counts[i];
                let denom = m_f * m_f * (1.0 - l) * (1.0 - l);
                let term1 = if denom > 0.0 {
                    (w / denom) * (1.0 - w / m_f)
                } else {
                    0.0
                };
                let term2 = (pi0s[i] - minpi0) * (pi0s[i] - minpi0);
                mse[i] = term1 + term2;
            }
            let mut argmin = 0;
            let mut best = mse[0];
            for (i, value) in mse.iter().enumerate().skip(1) {
                if *value < best {
                    best = *value;
                    argmin = i;
                }
            }
            let pi0 = pi0s[argmin].min(1.0);
            Ok(Pi0Result {
                pi0,
                pi0_lambda: pi0s,
                lambda: lambda_v,
                pi0_smooth: false,
            })
        }
    }
}

/// Upper-tail probabilities of each statistic under a normal fitted to the null
/// statistics.
///
/// The null mean and its standard deviation with one degree of freedom removed
/// are taken over the finite entries of `stat0`; each entry of `stat` then
/// yields `1 - Phi((v - mu) / sigma)`, computed as the source computes it,
/// `1 - 0.5 (1 + erf(z / sqrt(2)))`.
///
/// A degenerate null — a single finite value, or several identical ones — gives
/// `sigma == 0`, and the source returns `1.0` below the mean and `0.0` at or
/// above it, treating the null as a point mass. That is preserved. Non-finite
/// entries of `stat` yield `NaN`.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `stat0` is empty or holds no finite
/// value, both `std::invalid_argument` in the source, and
/// [`Error::InvalidRange`] when either input exceeds [`MAX_ITEMS`].
pub fn p_norm(stat: &[f64], stat0: &[f64]) -> Result<Vec<f64>> {
    check_items(stat.len(), "pNorm")?;
    check_items(stat0.len(), "pNorm")?;
    if stat0.is_empty() {
        return Err(bad("pNorm: stat0 must be non-empty".to_string()));
    }
    let s0: Vec<f64> = stat0.iter().copied().filter(|v| v.is_finite()).collect();
    let m = s0.len();
    if m == 0 {
        return Err(bad("pNorm: stat0 contains no finite values".to_string()));
    }
    let mut sum = 0.0;
    for value in &s0 {
        sum += *value;
    }
    let mu = sum / m as f64;
    let mut var = 0.0;
    if m > 1 {
        for value in &s0 {
            let d = *value - mu;
            var += d * d;
        }
        var /= (m - 1) as f64;
    }
    let sigma = var.sqrt();
    let sqrt2 = 2.0_f64.sqrt();
    let degenerate = not_greater(sigma, 0.0);

    let mut out = Vec::with_capacity(stat.len());
    for value in stat.iter().copied() {
        if !value.is_finite() {
            out.push(f64::NAN);
            continue;
        }
        if degenerate {
            out.push(if value < mu { 1.0 } else { 0.0 });
            continue;
        }
        let z = (value - mu) / sigma;
        let cdf = 0.5 * (1.0 + libm::erf(z / sqrt2));
        out.push(1.0 - cdf);
    }
    Ok(out)
}

/// The standard normal quantile, `Phi^-1(p)`.
///
/// The source calls `boost::math::quantile(normal_distribution<double>(0, 1),
/// p)`. Boost is not a dependency here, so this is Acklam's rational
/// approximation followed by one Halley step against
/// `0.5 erfc(-x / sqrt(2))`, which brings the relative error to a few units in
/// the last place — well inside the `1e-2` tolerance the class test's local-FDR
/// reference vectors are compared at, and inside the `1e-4` of its `pNorm`
/// vector.
fn standard_normal_quantile(p: f64) -> f64 {
    if p <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }
    const A: [f64; 6] = [
        -3.969683028665376e+01,
        2.209460984245205e+02,
        -2.759285104469687e+02,
        1.38357751867269e+02,
        -3.066479806614716e+01,
        2.506628277459239e+00,
    ];
    const B: [f64; 5] = [
        -5.447609879822406e+01,
        1.615858368580409e+02,
        -1.556989798598866e+02,
        6.680131188771972e+01,
        -1.328068155288572e+01,
    ];
    const C: [f64; 6] = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
        4.374664141464968e+00,
        2.938163982698783e+00,
    ];
    const D: [f64; 4] = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e+00,
        3.754408661907416e+00,
    ];
    const P_LOW: f64 = 0.02425;

    let mut x = if p < P_LOW {
        let q = (-2.0 * p.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= 1.0 - P_LOW {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    };
    // One Halley refinement; the exponential overflows for |x| beyond ~38, and
    // there the approximation is already at the limit of what f64 resolves.
    if x.abs() < 37.0 {
        let e = 0.5 * libm::erfc(-x / std::f64::consts::SQRT_2) - p;
        let u = e * (2.0 * PI).sqrt() * (x * x / 2.0).exp();
        x -= u / (1.0 + x * u / 2.0);
    }
    x
}

/// Local false discovery rate, the posterior probability that a hypothesis with
/// this p-value is null.
///
/// The finite p-values are transformed — probit or logit — their density is
/// estimated with [`kde_fft_eval`] at a bandwidth of
/// `bw_nrd0(transformed) * adj`, and the estimate is divided into the null
/// density scaled by `pi0`. Probit uses the standard normal as the null;
/// logit uses the logistic derivative `e^x / (1 + e^x)^2`. A zero or negative
/// density estimate yields `+inf`, which the source produces and which
/// `truncate` then caps at `1`.
///
/// Non-finite inputs are dropped and their positions come back `NaN`.
///
/// # Arguments
///
/// * `pi0` — proportion of true nulls, in `[0, 1]`.
/// * `options` — see [`LfdrOptions`]; `LfdrOptions::default()` is the source's
///   default argument list.
///
/// # Monotonicity
///
/// With `monotone` set, the estimates are sorted by p-value, swept left to
/// right taking a running **maximum** so they are non-decreasing in p, and then
/// returned to their input positions. The source's comment and its `minrank`
/// naming say it applies a `rankdata('min')` mapping, which would give tied
/// p-values one shared value; the code indexes its `assigned` flags by the
/// original position rather than by the value, so every position is visited
/// exactly once and the mapping is the plain inverse permutation. Tied
/// p-values therefore keep the distinct values the cumulative maximum left
/// them. This port reproduces the code, not the comment, and records the
/// discrepancy in `OpenMS_CPP_ISSUES.md`.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when a finite p-value is outside `[0, 1]` or
/// `pi0` is outside `[0, 1]` — the source's two `std::invalid_argument` throws
/// — when `adj` or `eps` is not finite or positive, which the source does not
/// check, and whatever [`kde_fft_eval`] returns. Returns
/// [`Error::InvalidRange`] when the input exceeds [`MAX_ITEMS`].
pub fn lfdr(p_values: &[f64], pi0: f64, options: &LfdrOptions) -> Result<Vec<f64>> {
    check_items(p_values.len(), "lfdr")?;
    let total = p_values.len();
    let mut out = vec![f64::NAN; total];
    if total == 0 {
        return Ok(out);
    }
    let keep: Vec<bool> = p_values.iter().map(|v| v.is_finite()).collect();
    let mut p: Vec<f64> = p_values.iter().copied().filter(|v| v.is_finite()).collect();
    let m = p.len();
    if m == 0 {
        return Ok(out);
    }
    if p.iter().any(|v| *v < 0.0 || *v > 1.0) {
        return Err(bad("lfdr: p-values not in [0,1]".to_string()));
    }
    if !(0.0..=1.0).contains(&pi0) {
        return Err(bad("lfdr: pi0 not in [0,1]".to_string()));
    }
    if !options.adj.is_finite() || not_greater(options.adj, 0.0) {
        return Err(bad(
            "lfdr: the bandwidth adjustment must be finite and positive".to_string(),
        ));
    }
    if !options.eps.is_finite() || not_greater(options.eps, 0.0) || options.eps >= 0.5 {
        return Err(bad("lfdr: eps must be finite and in (0, 0.5)".to_string()));
    }

    let mut x = vec![0.0; m];
    let mut lfdr_vec;
    match options.transform {
        LfdrTransform::Probit => {
            // Clip in place, as the source does: the clipped values are what
            // the null density is evaluated at, not just the transform input.
            for value in p.iter_mut() {
                if *value < options.eps {
                    *value = options.eps;
                }
                if *value > 1.0 - options.eps {
                    *value = 1.0 - options.eps;
                }
            }
            for i in 0..m {
                x[i] = standard_normal_quantile(p[i]);
            }
            let bw = bw_nrd0(&x)? * options.adj;
            let y = kde_fft_eval(&x, bw, options.gridsize, options.cut)?;
            let norm_const = 1.0 / (2.0 * PI).sqrt();
            lfdr_vec = vec![f64::NAN; m];
            for i in 0..m {
                let f0 = norm_const * (-0.5 * x[i] * x[i]).exp();
                lfdr_vec[i] = if not_greater(y[i], 0.0) {
                    f64::INFINITY
                } else {
                    pi0 * f0 / y[i]
                };
            }
        }
        LfdrTransform::Logit => {
            for i in 0..m {
                // The source does not clip here; `eps` appears inside the
                // log-odds instead, which keeps p == 0 and p == 1 finite.
                x[i] = ((p[i] + options.eps) / (1.0 - p[i] + options.eps)).ln();
            }
            let bw = bw_nrd0(&x)? * options.adj;
            let y = kde_fft_eval(&x, bw, options.gridsize, options.cut)?;
            lfdr_vec = vec![f64::NAN; m];
            for i in 0..m {
                let ex = x[i].exp();
                let denom = 1.0 + ex;
                let dx = ex / (denom * denom);
                lfdr_vec[i] = if not_greater(y[i], 0.0) {
                    f64::INFINITY
                } else {
                    (pi0 * dx) / y[i]
                };
            }
        }
    }

    if options.truncate {
        for value in lfdr_vec.iter_mut() {
            if *value > 1.0 {
                *value = 1.0;
            }
        }
    }

    if options.monotone && m > 1 {
        let order = argsort_asc(&p);
        let mut sorted: Vec<f64> = order.iter().map(|i| lfdr_vec[*i]).collect();
        for i in 1..m {
            if sorted[i] < sorted[i - 1] {
                sorted[i] = sorted[i - 1];
            }
        }
        let mut mapped = vec![f64::NAN; m];
        for (j, index) in order.iter().enumerate() {
            mapped[*index] = sorted[j];
        }
        lfdr_vec = mapped;
    }

    let mut k = 0;
    for (i, slot) in out.iter_mut().enumerate() {
        if keep[i] {
            *slot = lfdr_vec[k];
            k += 1;
        }
    }
    Ok(out)
}

/// Model-based FDR from posterior error probabilities, the IPF q-value.
///
/// Sorts the PEPs ascending, takes the running mean of the sorted values up to
/// each element's `max` tie rank, and writes the result back to the element's
/// original position. Tied PEPs therefore share one value, the cumulative sum
/// at the last of the tie divided by its rank.
///
/// **Any `NaN` invalidates the whole result.** The returned vector is then
/// entirely `NaN` and no error is raised. The header documents this explicitly
/// as differing from [`q_value`], which localises `NaN` to the offending
/// positions, and warns that callers wanting partial results must pre-filter.
/// It is reproduced here rather than repaired, because a caller reading a
/// single position cannot tell a repaired vector from a valid one.
///
/// The source is a template over the PEP type. Only the `double` instantiation
/// is reachable from ported code; an integral instantiation cannot hold a `NaN`
/// and so can never take the propagation branch, which is why the port takes
/// `f64` and asks an integral caller to widen.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] when the input exceeds [`MAX_ITEMS`]. An
/// empty input gives an empty result, as in the source.
pub fn compute_model_fdr(data: &[f64]) -> Result<Vec<f64>> {
    check_items(data.len(), "computeModelFDR")?;
    let n = data.len();
    let mut fdr = vec![f64::NAN; n];
    if n == 0 {
        return Ok(fdr);
    }
    if data.iter().any(|v| v.is_nan()) {
        return Ok(fdr);
    }
    let order = argsort_asc(data);
    let sorted: Vec<f64> = order.iter().map(|i| data[*i]).collect();
    let ranks = rankdata_f64(&sorted, RankMethod::Max, NanPolicy::Propagate)?;
    let mut cumsum = vec![0.0; n];
    let mut acc = 0.0;
    for i in 0..n {
        acc += sorted[i];
        cumsum[i] = acc;
    }
    for i in 0..n {
        let rank = ranks[i];
        if rank.is_nan() {
            fdr[order[i]] = f64::NAN;
            continue;
        }
        // Ranks are one-based and at most n, so this index is in range; the
        // source clamps it anyway and the clamp is kept.
        let index = (rank as usize).saturating_sub(1).min(n - 1);
        fdr[order[i]] = cumsum[index] / rank;
    }
    Ok(fdr)
}

/// Empirical p-values of `stat` against the null sample `stat0`.
///
/// Both samples are pooled and sorted descending; for the `i`-th observed
/// statistic in that order, the number of null statistics ranked above it is
/// its position minus `i`, and dividing by `stat0.len()` gives the p-value.
/// Each input then takes the value at `floor(rank) - 1` of its average rank
/// among the negated statistics, and everything at or below `1 / stat0.len()`
/// is raised to that floor, so no p-value is smaller than one null observation
/// can justify.
///
/// As with [`compute_model_fdr`], the source is a template; this takes `f64`.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when either input is empty, the source's
/// `std::invalid_argument`, and [`Error::InvalidRange`] when either exceeds
/// [`MAX_ITEMS`]. The source's `std::runtime_error` on an internal count
/// mismatch cannot arise, because the permutation is a bijection.
pub fn p_emp(stat: &[f64], stat0: &[f64]) -> Result<Vec<f64>> {
    check_items(stat.len(), "pEmp")?;
    check_items(stat0.len(), "pEmp")?;
    let m = stat.len();
    let m0 = stat0.len();
    if m == 0 || m0 == 0 {
        return Err(bad("pEmp: input arrays must be non-empty".to_string()));
    }
    let mut combined = Vec::with_capacity(m + m0);
    combined.extend_from_slice(stat);
    combined.extend_from_slice(stat0);

    let n = combined.len();
    let mut perm: Vec<usize> = (0..n).collect();
    // Stable descending sort, the source's `statc[i] > statc[j]` comparator.
    perm.sort_by(|i, j| combined[*j].total_cmp(&combined[*i]));

    let mut p = Vec::with_capacity(m);
    let mut seen = 0usize;
    for (position, original) in perm.iter().enumerate() {
        if *original < m {
            p.push((position as f64 - seen as f64) / m0 as f64);
            seen += 1;
        }
    }

    let neg: Vec<f64> = stat.iter().map(|v| -*v).collect();
    let ranks = rankdata_f64(&neg, RankMethod::Average, NanPolicy::Propagate)?;

    let mut out = Vec::with_capacity(m);
    let min_p = 1.0 / m0 as f64;
    for rank in &ranks {
        let floored = rank.floor();
        let index = if floored.is_finite() && floored >= 1.0 {
            ((floored as usize) - 1).min(p.len() - 1)
        } else {
            // A NaN rank reaches an unsigned conversion in the source; here it
            // takes the same clamped last index the source's clamp would give.
            p.len() - 1
        };
        let value = p[index];
        out.push(if value <= min_p { min_p } else { value });
    }
    Ok(out)
}
