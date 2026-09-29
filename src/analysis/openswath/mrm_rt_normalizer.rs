// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Finding retention-time normalization peptides in data: outlier removal and
//! RT coverage for the linear map from experimental to normalized (iRT)
//! retention time.
//!
//! Ports `ANALYSIS/OPENSWATH/MRMRTNormalizer.h` with `MRMRTNormalizer.cpp` at
//! core `bc9cc12`, a class of static members, which become free functions.
//! The tool built on it takes a description of RT peptides and their
//! normalized retention times and writes a transformation from the measured
//! RT space into the normalized space. The principle is adapted from
//! Escher, C. et al. (2012), *Using iRT, a normalized retention time for more
//! targeted measurement of peptides*, Proteomics 12: 1111-1121.
//!
//! Points are `(experimental_rt, theoretical_rt)` pairs throughout. See
//! `docs/RANSAC_SUPPORT.md` for the member mapping and the evidence.

use crate::math::ml::linear_regression::LinearRegression;
use crate::math::ml::ransac::{Ransac, RansacModel, RansacModelLinear};
use crate::{Error, Result};
use std::str::FromStr;

/// The confidence level the source hands every regression it fits here; it
/// affects none of the statistics this module reads.
const CONFIDENCE_INTERVAL: f64 = 0.95;

/// How [`remove_outliers_iterative`] picks the next outlier candidate: the
/// source's `method` string.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum OutlierMethod {
    /// `"iter_jackknife"`: the point whose removal gives the highest R-squared
    /// of the remaining ones.
    IterJackknife,
    /// `"iter_residual"`: the point with the largest absolute residual to the
    /// fit of all points.
    IterResidual,
}

impl OutlierMethod {
    /// The source's name for the method.
    pub fn name(self) -> &'static str {
        match self {
            Self::IterJackknife => "iter_jackknife",
            Self::IterResidual => "iter_residual",
        }
    }
}

impl FromStr for OutlierMethod {
    type Err = Error;

    /// Parse the source's method name.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for any other name. The source throws
    /// `Exception::IllegalArgument` for an unknown method only when it first
    /// needs a candidate, so it accepts any name for data that already meets
    /// the R-squared limit; parsing up front refuses it every time.
    fn from_str(name: &str) -> Result<Self> {
        match name {
            "iter_jackknife" => Ok(Self::IterJackknife),
            "iter_residual" => Ok(Self::IterResidual),
            other => Err(Error::InvalidValue(format!(
                "Method {other} is not a valid method for removeOutliersIterative"
            ))),
        }
    }
}

fn unable_to_fit(message: String) -> Error {
    Error::InvalidValue(format!(
        "UnableToFit-LinearRegression-RTNormalizer: {message}"
    ))
}

/// Remove outliers from a linear RT regression with RANSAC.
///
/// Two thresholds decide acceptance: a lower R-squared limit for the
/// regression, and a lower limit on the fraction of peptides kept (the
/// coverage). RANSAC picks the candidate set — `sampling_size` points per
/// sample, at most `max_iterations` iterations, an inlier deviating from the
/// fit by less than `max_rt_threshold` in the unit of the second dimension
/// (theoretical RT) — and the set is returned only if both limits hold.
///
/// Internally RANSAC gets `n = sampling_size`, `k = max_iterations`,
/// `t = max_rt_threshold^2` and `d = floor(coverage_limit * pairs.len())`.
///
/// `seed` seeds the RANSAC shuffle, as [`Ransac::new`] does. It is a native
/// addition: the source constructs its RANSAC with the default seed
/// `time(nullptr)`, so its result can change from one second to the next.
/// Passing the same seed the source's clock would have produced gives the
/// source's result.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] starting with
/// `UnableToFit-LinearRegression-RTNormalizer` where the source throws
/// `Exception::UnableToFit`: fewer than 5 sampled peptides, fewer than 30
/// input peptides, a best set whose R-squared is below `rsq_limit`, or one
/// with fewer points than `d`. An R-squared that cannot be computed — RANSAC
/// found no set, or its points share one x value — is reported the same way.
/// Errors from [`Ransac::ransac`] propagate, and a negative or non-finite
/// `coverage_limit`, which the source converts to `size_t` with undefined
/// behaviour, is refused.
pub fn remove_outliers_ransac(
    pairs: &[(f64, f64)],
    rsq_limit: f64,
    coverage_limit: f64,
    max_iterations: usize,
    max_rt_threshold: f64,
    sampling_size: usize,
    seed: u64,
) -> Result<Vec<(f64, f64)>> {
    let n = sampling_size;
    let k = max_iterations;
    let t = max_rt_threshold * max_rt_threshold;
    let d_value = coverage_limit * pairs.len() as f64;
    if !(d_value.is_finite() && d_value >= 0.0 && d_value < usize::MAX as f64) {
        return Err(Error::InvalidValue(format!(
            "RT normalizer: coverage limit {coverage_limit} gives no valid point count"
        )));
    }
    // `(size_t)(coverage_limit * pairs.size())` truncates toward zero.
    let d = d_value as usize;

    if n < 5 {
        return Err(unable_to_fit(format!(
            "WARNING: RANSAC: {n} sampled RT peptides is below limit of 5 peptides required for \
             the RANSAC outlier detection algorithm."
        )));
    }
    if pairs.len() < 30 {
        return Err(unable_to_fit(format!(
            "WARNING: RANSAC: {} input RT peptides is below limit of 30 peptides required for \
             the RANSAC outlier detection algorithm.",
            pairs.len()
        )));
    }

    let mut ransac: Ransac<RansacModelLinear> = Ransac::new(seed);
    let new_pairs = ransac.ransac(pairs, n, k, t, d, false)?;
    let best_rsq = RansacModelLinear.rm_rsq(&new_pairs).map_err(|error| {
        unable_to_fit(format!(
            "WARNING: the {} points RANSAC kept cannot be fitted: {error}",
            new_pairs.len()
        ))
    })?;

    if best_rsq < rsq_limit {
        return Err(unable_to_fit(format!(
            "WARNING: rsq: {best_rsq} is below limit of {rsq_limit}. Validate assays for \
             RT-peptides and adjust the limit for rsq or coverage."
        )));
    }
    if new_pairs.len() < d {
        return Err(unable_to_fit(format!(
            "WARNING: number of data points: {} is below limit of {d}. Validate assays for \
             RT-peptides and adjust the limit for rsq or coverage.",
            new_pairs.len()
        )));
    }
    Ok(new_pairs)
}

/// The candidate outlier by jackknife: leaving out each point in turn, the
/// one whose removal gives the highest R-squared of a first-order linear
/// regression of the remaining points.
///
/// The source's protected `jackknifeOutlierCandidate_`, public here so it can
/// be tested directly. Ties go to the first maximum, as `std::max_element`
/// chooses; a NaN R-squared never compares greater.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] starting with
/// `UnableToFit-LinearRegression` when a leave-one-out fit fails (fewer than
/// two distinct remaining `x` values), as the source throws
/// `Exception::UnableToFit`; [`Error::InvalidRange`] when `x` and `y` differ in
/// length or are empty.
pub fn jackknife_outlier_candidate(x: &[f64], y: &[f64]) -> Result<usize> {
    check_xy(x, y)?;
    let mut rsq = Vec::with_capacity(x.len());
    for i in 0..x.len() {
        let mut x_tmp = x.to_vec();
        let mut y_tmp = y.to_vec();
        x_tmp.remove(i);
        y_tmp.remove(i);
        let mut regression = LinearRegression::new();
        regression.compute_regression(CONFIDENCE_INTERVAL, &x_tmp, &y_tmp, true)?;
        rsq.push(regression.r_squared());
    }
    Ok(max_element(&rsq))
}

/// The candidate outlier by residual: the point with the largest absolute
/// residual to the linear fit of all points.
///
/// The source's protected `residualOutlierCandidate_`, public here so it can
/// be tested directly. Ties go to the first maximum.
///
/// # Errors
///
/// As [`jackknife_outlier_candidate`], for the one fit of all points.
pub fn residual_outlier_candidate(x: &[f64], y: &[f64]) -> Result<usize> {
    check_xy(x, y)?;
    let mut regression = LinearRegression::new();
    regression.compute_regression(CONFIDENCE_INTERVAL, x, y, true)?;
    let residuals: Vec<f64> = x
        .iter()
        .zip(y)
        .map(|(xi, yi)| (yi - (regression.intercept() + regression.slope() * xi)).abs())
        .collect();
    Ok(max_element(&residuals))
}

fn check_xy(x: &[f64], y: &[f64]) -> Result<()> {
    if x.is_empty() || x.len() != y.len() {
        return Err(Error::InvalidRange(format!(
            "outlier candidate: {} x values and {} y values",
            x.len(),
            y.len()
        )));
    }
    Ok(())
}

/// `std::max_element` with `operator<`: the first position no later value
/// exceeds. The slice is non-empty.
fn max_element(values: &[f64]) -> usize {
    let mut best = 0;
    for (index, value) in values.iter().enumerate().skip(1) {
        if values[best] < *value {
            best = index;
        }
    }
    best
}

/// Remove outliers from a linear RT regression one at a time until it is
/// good enough.
///
/// Two thresholds decide acceptance: a lower R-squared limit for the
/// regression, and the lower limit of peptide coverage — the loop stops before
/// the number of points left falls below `coverage_limit` times the input
/// size. Each round fits all remaining points; if R-squared is below
/// `rsq_limit`, `method` picks a candidate outlier, and it is removed if
/// `use_chauvenet` is off or Chauvenet's criterion (see [`chauvenet`]) holds
/// for it on the assumption that the residuals are normally distributed.
/// Otherwise the loop ends.
///
/// Returns the remaining pairs, in input order, when the R-squared limit was
/// reached.
///
/// # Source quirk
///
/// Chauvenet's criterion is evaluated on the residuals of **all input pairs**
/// to the current fit, at the candidate's position **among the remaining
/// points**. Once a point has been removed, the two indices refer to different
/// peptides. The source does this, its class test's expectations depend on it,
/// and it is reproduced; `docs/RANSAC_SUPPORT.md` lists it as a defect
/// candidate.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] for fewer than three pairs (the source's
/// `Exception::IllegalArgument`), and one starting with
/// `UnableToFit-LinearRegression-RTNormalizer` when the loop ends with the
/// R-squared still below the limit (`Exception::UnableToFit`). A regression
/// that fails on the way — all remaining `x` equal, or a non-finite
/// coordinate — propagates its error, as the source's exception does.
pub fn remove_outliers_iterative(
    pairs: &[(f64, f64)],
    rsq_limit: f64,
    coverage_limit: f64,
    use_chauvenet: bool,
    method: OutlierMethod,
) -> Result<Vec<(f64, f64)>> {
    if pairs.len() < 3 {
        return Err(Error::InvalidValue(
            "Need at least 3 data points to remove outliers for the regression.".into(),
        ));
    }
    let (mut x, mut y): (Vec<f64>, Vec<f64>) = pairs.iter().copied().unzip();
    let minimum = coverage_limit * pairs.len() as f64;
    let mut rsq = 0.0;

    // Every round removes a point or ends the loop, so it runs at most
    // `pairs.len()` times.
    while x.len() as f64 >= minimum && rsq < rsq_limit {
        let mut regression = LinearRegression::new();
        regression.compute_regression(CONFIDENCE_INTERVAL, &x, &y, true)?;
        rsq = regression.r_squared();
        if !(rsq < rsq_limit) {
            break;
        }
        let intercept = regression.intercept();
        let slope = regression.slope();
        let residuals: Vec<f64> = pairs
            .iter()
            .map(|(first, second)| (second - (intercept + first * slope)).abs())
            .collect();
        let pos = match method {
            OutlierMethod::IterJackknife => jackknife_outlier_candidate(&x, &y)?,
            OutlierMethod::IterResidual => residual_outlier_candidate(&x, &y)?,
        };
        if !use_chauvenet || chauvenet(&residuals, pos)? {
            x.remove(pos);
            y.remove(pos);
        } else {
            break;
        }
    }

    if rsq < rsq_limit {
        // If the rsq is below the limit, this is an indication that something
        // went wrong.
        return Err(unable_to_fit(format!(
            "WARNING: rsq: {rsq} is below limit of {rsq_limit}. Validate assays for RT-peptides \
             and adjust the limit for rsq or coverage."
        )));
    }
    Ok(x.into_iter().zip(y).collect())
}

/// Chauvenet's criterion probability of the value at `pos`: the two-sided
/// normal tail probability `erfc(|r - mean| / (stdev * sqrt(2)))`.
///
/// The mean and the population standard deviation (dividing by `n`, computed
/// as `sqrt(mean of squares - mean^2)`) are those of all `residuals`.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty slice or a `pos` outside it;
/// the source divides by zero or reads out of bounds.
pub fn chauvenet_probability(residuals: &[f64], pos: usize) -> Result<f64> {
    if pos >= residuals.len() {
        return Err(Error::InvalidRange(format!(
            "Chauvenet: position {pos} is outside {} residuals",
            residuals.len()
        )));
    }
    let n = residuals.len() as f64;
    let mut sum = 0.0;
    let mut sum_squares = 0.0;
    for residual in residuals {
        sum += residual;
        sum_squares += residual * residual;
    }
    let mean = sum / n;
    let stdev = (sum_squares / n - mean * mean).sqrt();
    let mut d = (residuals[pos] - mean).abs() / stdev;
    d /= 2.0_f64.sqrt();
    Ok(libm::erfc(d))
}

/// Chauvenet's criterion for the value at `pos`: `true` when its
/// [`chauvenet_probability`] is below `1 / (2 n)`, meaning the outlier can be
/// removed.
///
/// # Errors
///
/// As [`chauvenet_probability`].
pub fn chauvenet(residuals: &[f64], pos: usize) -> Result<bool> {
    let prob = chauvenet_probability(residuals, pos)?;
    let criterion = 1.0 / (2 * residuals.len()) as f64;
    Ok(prob < criterion)
}

/// Whether the RT normalization peptides cover the whole RT range: at least
/// `min_bins_filled` of `nr_bins` equal bins over `rt_range` must each hold at
/// least `min_peptides_per_bin` peptides.
///
/// `rt_range` is the (estimated) full range in the normalized (theoretical)
/// RT space, and the pairs are binned by their second, theoretical, value.
/// A value at or beyond the upper end of the range falls into the last bin, as
/// in the source (which also prints a warning); values in `(-1/nr_bins, 0)`
/// of the normalized range fall into the first, because the conversion to a
/// bin truncates toward zero.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `nr_bins` is not positive, and when a
/// pair falls below the first bin or its normalized position is not finite
/// (including a zero-width range). The source writes out of bounds or has
/// undefined behaviour in all of these cases.
pub fn compute_binned_coverage(
    rt_range: (f64, f64),
    pairs: &[(f64, f64)],
    nr_bins: i32,
    min_peptides_per_bin: i32,
    min_bins_filled: i32,
) -> Result<bool> {
    let bins = usize::try_from(nr_bins)
        .ok()
        .filter(|bins| *bins > 0)
        .ok_or_else(|| {
            Error::InvalidValue(format!("RT coverage: {nr_bins} bins is not a positive count"))
        })?;
    let mut bin_counter = vec![0_i64; bins];
    for (_, theoretical) in pairs {
        // A value between [0, 1), scaled to the bin count.
        let mut norm_rt = (theoretical - rt_range.0) / (rt_range.1 - rt_range.0);
        norm_rt *= f64::from(nr_bins);
        if !norm_rt.is_finite() {
            return Err(Error::InvalidValue(format!(
                "RT coverage: theoretical RT {theoretical} has no bin in range ({}, {})",
                rt_range.0, rt_range.1
            )));
        }
        let bin = if norm_rt >= f64::from(nr_bins) {
            bins - 1
        } else if norm_rt <= -1.0 {
            return Err(Error::InvalidValue(format!(
                "RT coverage: theoretical RT {theoretical} lies below the range ({}, {})",
                rt_range.0, rt_range.1
            )));
        } else {
            // -1 < norm_rt < nr_bins: truncation toward zero gives 0..nr_bins-1.
            norm_rt as usize
        };
        bin_counter[bin] += 1;
    }
    let bins_filled = bin_counter
        .iter()
        .filter(|count| **count >= i64::from(min_peptides_per_bin))
        .count();
    Ok(bins_filled as i64 >= i64::from(min_bins_filled))
}
