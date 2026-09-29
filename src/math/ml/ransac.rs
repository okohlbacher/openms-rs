// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! RANSAC outlier detection with pluggable models.
//!
//! Ports the four `ML/RANSAC` headers at core `bc9cc12`:
//!
//! - `RANSACModel.h`: the CRTP base `RansacModel<ModelT>` becomes the
//!   [`RansacModel`] trait, its `DPair`/`DVec`/`ModelParameters` typedefs the
//!   [`DPair`] and [`ModelParameters`] aliases;
//! - `RANSACModelLinear.h`/`.cpp`: [`RansacModelLinear`], over
//!   [`LinearRegression`];
//! - `RANSACModelQuadratic.h`/`.cpp`: [`RansacModelQuadratic`], over
//!   [`QuadraticRegression`];
//! - `RANSAC.h` (header-only): [`RansacParam`] and [`Ransac`].
//!
//! The algorithm follows the SciPy cookbook reference
//! (<http://wiki.scipy.org/Cookbook/RANSAC>) that the source was implemented
//! and tested against, with the source's own refinements: a candidate that
//! explains more points always wins, and equal counts are decided by the
//! residual sum of squares. Samples are drawn by shuffling the whole input
//! with [`RandomShuffler`], so a seed reproduces the source's result exactly.
//! See `docs/RANSAC_SUPPORT.md` for the member mapping and the evidence.

use super::linear_regression::LinearRegression;
use super::quadratic_regression::QuadraticRegression;
use crate::math::random_shuffler::RandomShuffler;
use crate::{Error, Result};
use std::fmt;

/// One data point `(dim1, dim2)`: the source's `RansacModel<>::DPair`.
pub type DPair = (f64, f64);

/// Fitted model coefficients: the source's `RansacModel<>::ModelParameters`.
/// Linear models store `[intercept, slope]`, quadratic ones `[a, b, c]`.
pub type ModelParameters = Vec<f64>;

/// The upper bound on `k * pairs.len()`, the work one [`Ransac::ransac`] call
/// may do: every iteration shuffles and scans the whole input.
///
/// Native-only guard, checked before anything is allocated; the source runs
/// whatever it is handed.
pub const MAX_RANSAC_WORK: u128 = 1 << 32;

/// The parameters of one RANSAC run.
///
/// All four numbers default to zero and `relative_d` to `false`, as the
/// source's default constructor sets them.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RansacParam {
    /// Data points: the minimum number of data points required to fit the
    /// model.
    pub n: usize,
    /// Iterations: the maximum number of iterations allowed in the algorithm.
    pub k: usize,
    /// Threshold value for determining when a data point fits a model. It is
    /// the maximal *squared* deviation, in units of the second dimension
    /// (`dim2`).
    pub t: f64,
    /// The number of close data values (according to `t`) required to assert
    /// that a model fits well to data.
    pub d: usize,
    /// Whether `d` is a percentage (0-100) of the input size rather than a
    /// count.
    pub relative_d: bool,
}

impl RansacParam {
    /// The source's full constructor.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `relative_d` is set and `d` is 100
    /// or more, where the source throws `Exception::Precondition`: use a lower
    /// value; the more outliers you expect, the lower it should be.
    pub fn new(n: usize, k: usize, t: f64, d: usize, relative_d: bool) -> Result<Self> {
        if relative_d {
            check_relative_d(d)?;
        }
        Ok(Self {
            n,
            k,
            t,
            d,
            relative_d,
        })
    }
}

impl fmt::Display for RansacParam {
    /// The source's `toString()`, character for character; `t` prints as
    /// Rust formats a double, which may differ from the source's
    /// `std::stringstream` (six significant digits) for long values.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "RANSAC param:\n  n: {}\n  k: {} iterations\n  t: {} threshold\n  d: {} inliers\n\n",
            self.n, self.k, self.t, self.d
        )
    }
}

fn check_relative_d(d: usize) -> Result<()> {
    if d >= 100 {
        Err(Error::InvalidValue(
            "RANSAC: Relative 'd' >= 100% given. Use a lower value; the more outliers you \
             expect, the lower it should be."
                .into(),
        ))
    } else {
        Ok(())
    }
}

/// A model RANSAC can fit: the source's `RansacModel<ModelT>` interface.
///
/// The source dispatches statically through CRTP to `rm_*_impl` functions; a
/// trait expresses the same contract. Implementations must be deterministic,
/// because [`Ransac`]'s reproducibility rests on it.
pub trait RansacModel {
    /// Fit the model to `data` and return its parameters.
    ///
    /// # Errors
    ///
    /// When the model cannot be fitted to these points. [`Ransac`] skips an
    /// iteration whose initial sample cannot be fitted, as the source catches
    /// every exception there.
    fn rm_fit(&self, data: &[DPair]) -> Result<ModelParameters>;

    /// A goodness-of-fit value of `data` under a freshly fitted model.
    ///
    /// For [`RansacModelLinear`] this is R-squared; for
    /// [`RansacModelQuadratic`] it is, as in the source, the chi-squared sum.
    ///
    /// # Errors
    ///
    /// When the model cannot be fitted to these points.
    fn rm_rsq(&self, data: &[DPair]) -> Result<f64>;

    /// The residual sum of squares of `data` under the model `coefficients`.
    ///
    /// # Errors
    ///
    /// When `coefficients` has fewer entries than the model needs.
    fn rm_rss(&self, data: &[DPair], coefficients: &[f64]) -> Result<f64>;

    /// The points of `data` whose squared residual under `coefficients` is
    /// strictly below `max_threshold`, in input order.
    ///
    /// The source's `@brief` says "<= max_threshold"; both of its models
    /// compare with `<`, and so do these.
    ///
    /// # Errors
    ///
    /// When `coefficients` has fewer entries than the model needs.
    fn rm_inliers(
        &self,
        data: &[DPair],
        coefficients: &[f64],
        max_threshold: f64,
    ) -> Result<Vec<DPair>>;
}

fn coefficients<const N: usize>(coefficients: &[f64], model: &str) -> Result<[f64; N]> {
    coefficients
        .get(..N)
        .and_then(|head| <[f64; N]>::try_from(head).ok())
        .ok_or_else(|| {
            Error::InvalidValue(format!(
                "{model} needs {N} coefficients, got {}",
                coefficients.len()
            ))
        })
}

fn split(data: &[DPair]) -> (Vec<f64>, Vec<f64>) {
    data.iter().copied().unzip()
}

/// The straight-line model `y = c0 + c1 x`, fitted by [`LinearRegression`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RansacModelLinear;

impl RansacModel for RansacModelLinear {
    /// `[intercept, slope]` of an unweighted least-squares line, without the
    /// goodness statistics.
    ///
    /// # Errors
    ///
    /// As [`LinearRegression::compute_regression`]: fewer than two distinct
    /// `x` values or a non-finite centroid.
    fn rm_fit(&self, data: &[DPair]) -> Result<ModelParameters> {
        let (x, y) = split(data);
        let mut regression = LinearRegression::new();
        regression.compute_regression(0.95, &x, &y, false)?;
        Ok(vec![regression.intercept(), regression.slope()])
    }

    /// The squared Pearson coefficient of `data`, computed on the fly.
    ///
    /// With exactly two points this is `0.0`: the source's
    /// `LinearRegression` computes goodness statistics only for more than two
    /// points, and reads back its initial value otherwise.
    ///
    /// # Errors
    ///
    /// As [`rm_fit`](Self::rm_fit).
    fn rm_rsq(&self, data: &[DPair]) -> Result<f64> {
        let (x, y) = split(data);
        let mut regression = LinearRegression::new();
        regression.compute_regression(0.95, &x, &y, true)?;
        Ok(regression.r_squared())
    }

    fn rm_rss(&self, data: &[DPair], coefficients_in: &[f64]) -> Result<f64> {
        let [intercept, slope] = coefficients::<2>(coefficients_in, "linear RANSAC model")?;
        let mut rss = 0.0;
        for (x, y) in data {
            let residual = y - (intercept + slope * x);
            rss += residual * residual;
        }
        Ok(rss)
    }

    fn rm_inliers(
        &self,
        data: &[DPair],
        coefficients_in: &[f64],
        max_threshold: f64,
    ) -> Result<Vec<DPair>> {
        let [intercept, slope] = coefficients::<2>(coefficients_in, "linear RANSAC model")?;
        Ok(data
            .iter()
            .copied()
            .filter(|(x, y)| {
                let residual = y - (intercept + slope * x);
                residual * residual < max_threshold
            })
            .collect())
    }
}

/// The quadratic model `y = a + b x + c x^2`, fitted by
/// [`QuadraticRegression`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RansacModelQuadratic;

impl RansacModel for RansacModelQuadratic {
    /// `[a, b, c]` of an unweighted least-squares quadratic.
    ///
    /// # Errors
    ///
    /// As [`QuadraticRegression::compute_regression`]: a singular system.
    fn rm_fit(&self, data: &[DPair]) -> Result<ModelParameters> {
        let (x, y) = split(data);
        let mut regression = QuadraticRegression::new();
        regression.compute_regression(&x, &y)?;
        Ok(vec![regression.a(), regression.b(), regression.c()])
    }

    /// The **chi-squared sum** of a fresh fit, not an R-squared: the source's
    /// `rm_rsq_impl` returns `QuadraticRegression::getChiSquared()` under this
    /// name, and callers of the quadratic model get that value.
    ///
    /// # Errors
    ///
    /// As [`rm_fit`](Self::rm_fit).
    fn rm_rsq(&self, data: &[DPair]) -> Result<f64> {
        let (x, y) = split(data);
        let mut regression = QuadraticRegression::new();
        regression.compute_regression(&x, &y)?;
        Ok(regression.chi_squared())
    }

    fn rm_rss(&self, data: &[DPair], coefficients_in: &[f64]) -> Result<f64> {
        let [a, b, c] = coefficients::<3>(coefficients_in, "quadratic RANSAC model")?;
        let mut rss = 0.0;
        for (x, y) in data {
            let diff = y - QuadraticRegression::eval_with(a, b, c, *x);
            rss += diff * diff;
        }
        Ok(rss)
    }

    fn rm_inliers(
        &self,
        data: &[DPair],
        coefficients_in: &[f64],
        max_threshold: f64,
    ) -> Result<Vec<DPair>> {
        let [a, b, c] = coefficients::<3>(coefficients_in, "quadratic RANSAC model")?;
        Ok(data
            .iter()
            .copied()
            .filter(|(x, y)| {
                let diff = y - QuadraticRegression::eval_with(a, b, c, *x);
                diff * diff < max_threshold
            })
            .collect())
    }
}

/// The RANSAC outlier detection algorithm over a model `M`.
///
/// Each iteration shuffles the whole input, fits the model to the first `n`
/// points, and collects the remaining points whose squared residual is below
/// `t`. When more than `d` points join — or every remaining point does — the
/// model is refitted to the sample plus its inliers, and that set becomes the
/// best one if it is larger than the best so far, or equally large with a
/// smaller residual sum of squares.
#[derive(Clone, Debug, PartialEq)]
pub struct Ransac<M = RansacModelLinear> {
    shuffler: RandomShuffler,
    model: M,
}

impl<M: RansacModel + Default> Ransac<M> {
    /// The source's `explicit RANSAC(uint64_t seed)`.
    ///
    /// The source passes the seed to `RandomShuffler`'s `int` constructor, so
    /// only its low 32 bits count, read as a signed value and sign-extended
    /// back to 64 bits; [`set_seed`](Self::set_seed) uses the full value. Both
    /// quirks are reproduced, so a seed selects the same permutations as in the
    /// source. The source's default seed is `time(nullptr)`; this port has no
    /// time-seeded constructor, so every run is reproducible.
    pub fn new(seed: u64) -> Self {
        Self {
            shuffler: RandomShuffler::from_int_seed(seed as i32),
            model: M::default(),
        }
    }
}

impl<M: RansacModel> Ransac<M> {
    /// A RANSAC run over an explicit model value and shuffler.
    ///
    /// Native constructor for models that carry state; the source's models
    /// are stateless and default-constructed.
    pub fn with_model(model: M, shuffler: RandomShuffler) -> Self {
        Self { shuffler, model }
    }

    /// Reseed the random shuffle with the full 64-bit `seed`.
    pub fn set_seed(&mut self, seed: u64) {
        self.shuffler.seed(seed);
    }

    /// [`ransac`](Self::ransac) with the parameters bundled in `p`.
    ///
    /// # Errors
    ///
    /// As [`ransac`](Self::ransac).
    pub fn ransac_with_param(&mut self, pairs: &[DPair], p: &RansacParam) -> Result<Vec<DPair>> {
        self.ransac(pairs, p.n, p.k, p.t, p.d, p.relative_d)
    }

    /// The points that fit the model well, as found by RANSAC.
    ///
    /// `pairs` holds `(dim1, dim2)` points. If possible, restrict `n` to the
    /// minimal number of points the model needs for a fit — 2 for linear, 3
    /// for quadratic: any higher number increases the chance of including an
    /// outlier in the sample, hence a lost iteration. At most `k` iterations
    /// run; they stop early once one model explains every point. `t` is the
    /// maximal squared deviation of an inlier, in units of `dim2`. More than
    /// `d` points besides the sample must fit for a model to count.
    ///
    /// A relative `d` (1-99 %) helps when the number of points is not known in
    /// advance but the fraction of outliers roughly is: expecting 20 %
    /// outliers, `d = 60` with `relative_d` (60 % inliers, with some margin) is
    /// a good bet for a larger input. The `n` sample points cannot also become
    /// inliers, so leave room for them. A relative `d` becomes
    /// `pairs.len() * d / 100` in integer arithmetic.
    ///
    /// Returns the sample and inliers of the best model, in shuffled order
    /// (the source says "data will be unsorted"), or an empty vector when no
    /// iteration found enough inliers.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] where the source throws
    /// `Exception::Precondition`: a relative `d` of 100 or more, or no more
    /// points than `n`. Also for a non-finite coordinate — the source has no
    /// such check and lets NaN fail every comparison — and when
    /// `k * pairs.len()` exceeds [`MAX_RANSAC_WORK`]. An error from refitting
    /// an accepted sample plus its inliers propagates, as it does in the
    /// source, where only the initial sample's fit is guarded; the shuffler
    /// has then advanced.
    pub fn ransac(
        &mut self,
        pairs: &[DPair],
        n: usize,
        k: usize,
        t: f64,
        d: usize,
        relative_d: bool,
    ) -> Result<Vec<DPair>> {
        let mut d = d;
        if relative_d {
            check_relative_d(d)?;
            d = pairs.len().saturating_mul(d) / 100;
        }
        if pairs.len() <= n {
            return Err(Error::InvalidValue(format!(
                "RANSAC: Number of total data points ({}) must be larger than number of initial \
                 points (n={n}).",
                pairs.len()
            )));
        }
        if pairs.iter().any(|(x, y)| !(x.is_finite() && y.is_finite())) {
            return Err(Error::InvalidValue("RANSAC: non-finite data point".into()));
        }
        let work = (k as u128).saturating_mul(pairs.len() as u128);
        if work > MAX_RANSAC_WORK {
            return Err(Error::InvalidValue(format!(
                "RANSAC: {k} iterations over {} points exceed the work limit of {MAX_RANSAC_WORK}",
                pairs.len()
            )));
        }

        let mut shuffled = pairs.to_vec();
        let mut best_data: Vec<DPair> = Vec::new();
        let mut best_error = f64::MAX;
        for _ in 0..k {
            // Stop once the model already includes all points.
            if best_data.len() == pairs.len() {
                break;
            }
            self.shuffler.portable_random_shuffle(&mut shuffled);
            let (sample, rest) = shuffled.split_at(n);
            // Fitting might fail if the sampled points are unfortunate.
            let Ok(coefficients) = self.model.rm_fit(sample) else {
                continue;
            };
            let also_inliers = self.model.rm_inliers(rest, &coefficients, t)?;
            // The second condition is the maximum number of inliers possible:
            // all of the remaining data.
            if also_inliers.len() > d || also_inliers.len() >= rest.len() {
                let mut better_data = Vec::with_capacity(n + also_inliers.len());
                better_data.extend_from_slice(sample);
                better_data.extend_from_slice(&also_inliers);
                let better_coefficients = self.model.rm_fit(&better_data)?;
                let better_error = self.model.rm_rss(&better_data, &better_coefficients)?;
                // More explained points win outright; they all pass the
                // threshold, while their absolute residual sum can only grow.
                // Equal counts are decided by the residual sum of squares.
                if better_data.len() > best_data.len()
                    || (better_data.len() == best_data.len() && better_error < best_error)
                {
                    best_error = better_error;
                    best_data = better_data;
                }
            }
        }
        Ok(best_data)
    }
}
