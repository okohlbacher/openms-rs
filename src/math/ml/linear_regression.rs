// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Least-squares fits to a straight line, `Y(c, x) = c0 + c1 x`, with the
//! goodness-of-fit statistics the source computes alongside.
//!
//! Ports `ML/REGRESSION/LinearRegression.h` with `LinearRegression.cpp` at core
//! `bc9cc12`. See `docs/RANSAC_SUPPORT.md`, which covers the whole
//! `ML/REGRESSION` and `ML/RANSAC` family, for the member mapping and evidence.
//!
//! The source fits through the Geometric Tools Engine: the unweighted fit is
//! `gte::ApprHeightLine2`, the weighted one a 2x2 normal-equation system solved
//! by `gte::LinearSystem` (an explicit inverse). Both are transcribed here
//! operation for operation, because the order of those operations fixes the
//! last bits of every coefficient and every consumer — RANSAC, the RT
//! normalizer — compares them against thresholds. The Student's t quantile
//! comes from `statrs` in place of Boost.Math; it feeds only
//! [`LinearRegression::t_value`], [`LinearRegression::lower`] and
//! [`LinearRegression::upper`].

use crate::math::statistic_functions::{covariance, mean, variance_with_mean};
use crate::{Error, Result};
use statrs::distribution::{ContinuousCDF, StudentsT};

/// A straight-line fit and its goodness-of-fit statistics.
///
/// Next to the intercept with the y-axis and the slope of the fitted line, a
/// fit with goodness statistics records the squared Pearson coefficient, the
/// value of the t-distribution, the standard deviation and mean of the
/// residuals, the standard error of the slope, the intercept with the x-axis
/// (useful for additive series experiments), the borders of the confidence
/// interval around that x-intercept, the chi-squared value and a relative
/// standard deviation.
///
/// Every statistic reads `0.0` until a fit succeeds. After a fit **without**
/// goodness statistics — or with fewer than three points, where the source
/// skips them — only [`slope`](Self::slope), [`intercept`](Self::intercept)
/// and [`chi_squared`](Self::chi_squared) are meaningful; the others keep the
/// values of the previous fit, exactly as the source's members do.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LinearRegression {
    intercept: f64,
    slope: f64,
    x_intercept: f64,
    lower: f64,
    upper: f64,
    t_star: f64,
    r_squared: f64,
    stand_dev_residuals: f64,
    mean_residuals: f64,
    stand_error_slope: f64,
    chi_squared: f64,
    rsd: f64,
}

/// The `UnableToFit-LinearRegression` failure, as the source names it.
fn unable_to_fit(points: usize) -> Error {
    Error::InvalidValue(format!(
        "UnableToFit-LinearRegression: could not fit a linear model to the data ({points} points)"
    ))
}

fn check_lengths(x: &[f64], other: &[f64], what: &str) -> Result<()> {
    if x.len() == other.len() {
        Ok(())
    } else {
        Err(Error::InvalidRange(format!(
            "linear regression: {} x values but {} {what} values",
            x.len(),
            other.len()
        )))
    }
}

/// `gte::ApprHeightLine2<double>::Fit` over all points: the height line
/// through the centroid. Returns `(slope, intercept)`, or `None` where the
/// source's fit reports failure (a non-finite centroid, or all `x` equal).
fn height_line_fit(x: &[f64], y: &[f64]) -> Option<(f64, f64)> {
    let count = x.len();
    let mut mean_x = 0.0;
    let mut mean_y = 0.0;
    for (xi, yi) in x.iter().zip(y) {
        mean_x += xi;
        mean_y += yi;
    }
    // gte's `Vector /= scalar` multiplies by the reciprocal and zeroes the
    // vector on a zero divisor.
    if count == 0 {
        mean_x = 0.0;
        mean_y = 0.0;
    } else {
        let inverse = 1.0 / count as f64;
        mean_x *= inverse;
        mean_y *= inverse;
    }
    if !(mean_x.is_finite() && mean_y.is_finite()) {
        return None;
    }
    let mut covar00 = 0.0;
    let mut covar01 = 0.0;
    for (xi, yi) in x.iter().zip(y) {
        let dx = xi - mean_x;
        let dy = yi - mean_y;
        covar00 += dx * dx;
        covar01 += dx * dy;
    }
    // `!(covar00 > 0)` also rejects a NaN, as the source's `if (covar00 > 0)`.
    if covar00 > 0.0 {
        let slope = covar01 / covar00;
        // `intercept_ = -slope_ * mean[0] + mean[1]`.
        Some((slope, -slope * mean_x + mean_y))
    } else {
        None
    }
}

/// `sum((y - (slope * x + intercept))^2)`, in slice order.
fn chi_square(x: &[f64], y: &[f64], slope: f64, intercept: f64) -> f64 {
    let mut chi = 0.0;
    for (xi, yi) in x.iter().zip(y) {
        let residual = yi - LinearRegression::compute_point_y(*xi, slope, intercept);
        chi += residual * residual;
    }
    chi
}

/// `sum(w * (y - (slope * x + intercept))^2)`, in slice order.
fn weighted_chi_square(x: &[f64], y: &[f64], w: &[f64], slope: f64, intercept: f64) -> f64 {
    let mut chi = 0.0;
    for ((xi, yi), wi) in x.iter().zip(y).zip(w) {
        let residual = yi - LinearRegression::compute_point_y(*xi, slope, intercept);
        chi += wi * (residual * residual);
    }
    chi
}

impl LinearRegression {
    /// A regression with every coefficient and statistic at `0.0`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Fit `Y = c0 + c1 X` to the points `(x[i], y[i])` by least squares.
    ///
    /// For an "x %" confidence interval pass `confidence_interval_p = x / 100`;
    /// the 95 % interval, for example, is supposed to have a 95 % chance of
    /// containing the true value of the parameter. The interval is the one
    /// around the x-intercept, and like every other goodness statistic it is
    /// computed only when `compute_goodness` is set *and* there are more than
    /// two points; otherwise only the slope, intercept and chi-squared value
    /// change.
    ///
    /// The fit is the source's `gte::ApprHeightLine2`: centroid by summing and
    /// multiplying by `1 / n`, slope `Sxy / Sxx` about that centroid, intercept
    /// `-slope * mean_x + mean_y`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] starting with
    /// `UnableToFit-LinearRegression` where the source throws
    /// `Exception::UnableToFit`: fewer than two distinct `x` values, or a
    /// non-finite centroid (any NaN or infinite coordinate, or a sum that
    /// overflows). Returns [`Error::InvalidRange`] when `x` and `y` differ in
    /// length (the source reads `y` past its end), and [`Error::InvalidValue`]
    /// when goodness statistics are requested with a confidence level outside
    /// `(-1, 1)`, where Boost's t quantile raises a domain or overflow error.
    ///
    /// On error `self` is unchanged. The source has already overwritten the
    /// slope, intercept and chi-squared value with those of a failed fit
    /// (zeros, and the plain sum of squares of `y`) when it throws.
    pub fn compute_regression(
        &mut self,
        confidence_interval_p: f64,
        x: &[f64],
        y: &[f64],
        compute_goodness: bool,
    ) -> Result<()> {
        check_lengths(x, y, "y")?;
        let (slope, intercept) = height_line_fit(x, y).ok_or_else(|| unable_to_fit(x.len()))?;
        let mut next = self.clone();
        next.slope = slope;
        next.intercept = intercept;
        next.chi_squared = chi_square(x, y, slope, intercept);
        if compute_goodness && x.len() > 2 {
            next.compute_goodness(x, y, confidence_interval_p)?;
        }
        *self = next;
        Ok(())
    }

    /// Fit `Y = c0 + c1 X` to the points `(x[i], y[i])` weighted by `w[i]`.
    ///
    /// Solves the weighted normal equations
    /// `[Swxx Swx; Swx Sw] [c1; c0] = [Swxy; Swy]` through the explicit 2x2
    /// inverse, as `gte::LinearSystem::Solve` does. `confidence_interval_p` and
    /// `compute_goodness` are as for
    /// [`compute_regression`](Self::compute_regression).
    ///
    /// The goodness statistics are those of the **unweighted** points with the
    /// weighted line: the source hands `computeGoodness_` the bare `x` and `y`,
    /// and this port keeps that. The chi-squared value is the weighted one.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] starting with
    /// `UnableToFit-LinearRegression` when the system is singular (its
    /// determinant is exactly zero), and for any non-finite coordinate or
    /// weight. The source has no finiteness check on this path and returns
    /// NaN coefficients; this port refuses instead, so a NaN never reaches a
    /// consumer's threshold comparison. Length mismatches and the confidence
    /// level are checked as in [`compute_regression`](Self::compute_regression).
    /// On error `self` is unchanged; the source has already recomputed the
    /// chi-squared value with the previous coefficients when it throws.
    pub fn compute_regression_weighted(
        &mut self,
        confidence_interval_p: f64,
        x: &[f64],
        y: &[f64],
        w: &[f64],
        compute_goodness: bool,
    ) -> Result<()> {
        check_lengths(x, y, "y")?;
        check_lengths(x, w, "weight")?;
        if x.iter().chain(y).chain(w).any(|value| !value.is_finite()) {
            return Err(Error::InvalidValue(
                "UnableToFit-LinearRegression: non-finite coordinate or weight".into(),
            ));
        }
        let mut sum_x = 0.0;
        let mut sum_y = 0.0;
        let mut sum_xx = 0.0;
        let mut sum_xy = 0.0;
        let mut sum_w = 0.0;
        for ((xi, yi), wi) in x.iter().zip(y).zip(w) {
            sum_x += wi * xi;
            sum_y += wi * yi;
            sum_xx += wi * xi * xi;
            sum_xy += wi * xi * yi;
            sum_w += wi;
        }
        // gte::Inverse(Matrix2x2) on A = [sum_xx sum_x; sum_x sum_w].
        let det = sum_xx * sum_w - sum_x * sum_x;
        if det == 0.0 {
            return Err(unable_to_fit(x.len()));
        }
        let inv_det = 1.0 / det;
        let inv00 = sum_w * inv_det;
        let inv01 = -sum_x * inv_det;
        let inv10 = -sum_x * inv_det;
        let inv11 = sum_xx * inv_det;
        // gte's Matrix * Vector starts each row sum at zero.
        let slope = 0.0 + inv00 * sum_xy + inv01 * sum_y;
        let intercept = 0.0 + inv10 * sum_xy + inv11 * sum_y;
        if !(slope.is_finite() && intercept.is_finite()) {
            return Err(unable_to_fit(x.len()));
        }
        let mut next = self.clone();
        next.slope = slope;
        next.intercept = intercept;
        next.chi_squared = weighted_chi_square(x, y, w, slope, intercept);
        if compute_goodness && x.len() > 2 {
            next.compute_goodness(x, y, confidence_interval_p)?;
        }
        *self = next;
        Ok(())
    }

    /// `computeGoodness_`: every statistic beyond slope, intercept and
    /// chi-squared, for `x.len() > 2`.
    ///
    /// The formulas are the source's, including its non-standard relative
    /// standard deviation (after Mayr et al. 2006). Degenerate data yields the
    /// same infinities and NaNs as the source: a vertical spread of zero makes
    /// the squared Pearson coefficient `0 / 0`, a zero slope an infinite
    /// x-intercept. The source prints a diagnostic block to standard output
    /// when the relative standard deviation is negative; this port does not.
    fn compute_goodness(&mut self, x: &[f64], y: &[f64], confidence_interval_p: f64) -> Result<()> {
        if !(confidence_interval_p > -1.0 && confidence_interval_p < 1.0) {
            return Err(Error::InvalidValue(format!(
                "linear regression: confidence level {confidence_interval_p} is outside (-1, 1)"
            )));
        }
        let n = x.len();
        let n_f = n as f64;
        let x_mean = mean(x)?;
        let y_mean = mean(y)?;
        let var_x = variance_with_mean(x, x_mean)?;
        let var_y = variance_with_mean(y, y_mean)?;
        let cov_xy = covariance(x, y)?;
        let s_xx = var_x * (n - 1) as f64;

        let r_squared = (cov_xy * cov_xy) / (var_x * var_y);

        let mut sum = 0.0;
        for (xi, yi) in x.iter().zip(y) {
            sum += (yi - (self.intercept + self.slope * xi)).abs();
        }
        let mean_residuals = sum / n_f;
        let stand_dev_residuals = ((self.chi_squared - (sum * sum) / n_f) / (n - 1) as f64).sqrt();
        let stand_error_slope = stand_dev_residuals / s_xx.sqrt();
        let x_intercept = -(self.intercept / self.slope);

        let p = 1.0 - (1.0 - confidence_interval_p) / 2.0;
        let t_star = students_t_quantile((n - 2) as f64, p)?;

        let mut g = t_star / (self.slope / stand_error_slope);
        g *= g;
        let left = (x_intercept - x_mean) * g;
        let bottom = 1.0 - g;
        let d = x_intercept - x_mean;
        let right =
            t_star * (stand_dev_residuals / self.slope) * ((d * d) / s_xx + (bottom / n_f)).sqrt();
        let mut lower = x_intercept + (left + right) / bottom;
        let mut upper = x_intercept + (left - right) / bottom;
        if lower > upper {
            std::mem::swap(&mut lower, &mut upper);
        }

        let mut tmp = 0.0;
        for xi in x {
            tmp += (xi - x_mean) * (xi - x_mean);
        }
        let rsd = (100.0 / x_intercept.abs())
            * (stand_dev_residuals / self.slope)
            * ((1.0 / n_f) * (y_mean / (self.slope * self.slope * tmp))).sqrt();

        self.r_squared = r_squared;
        self.mean_residuals = mean_residuals;
        self.stand_dev_residuals = stand_dev_residuals;
        self.stand_error_slope = stand_error_slope;
        self.x_intercept = x_intercept;
        self.t_star = t_star;
        self.lower = lower;
        self.upper = upper;
        self.rsd = rsd;
        Ok(())
    }

    /// The y-intercept `c0` of the fitted line.
    pub fn intercept(&self) -> f64 {
        self.intercept
    }

    /// The slope `c1` of the fitted line.
    pub fn slope(&self) -> f64 {
        self.slope
    }

    /// The intercept of the fitted line with the x-axis, `-c0 / c1`.
    pub fn x_intercept(&self) -> f64 {
        self.x_intercept
    }

    /// The lower border of the confidence interval around the x-intercept.
    pub fn lower(&self) -> f64 {
        self.lower
    }

    /// The upper border of the confidence interval around the x-intercept.
    pub fn upper(&self) -> f64 {
        self.upper
    }

    /// The value of the t-distribution: the two-sided quantile for the
    /// requested confidence level with `n - 2` degrees of freedom.
    pub fn t_value(&self) -> f64 {
        self.t_star
    }

    /// The squared Pearson correlation coefficient of the points.
    pub fn r_squared(&self) -> f64 {
        self.r_squared
    }

    /// The standard deviation of the residuals.
    pub fn stand_dev_res(&self) -> f64 {
        self.stand_dev_residuals
    }

    /// The mean absolute residual.
    pub fn mean_res(&self) -> f64 {
        self.mean_residuals
    }

    /// The standard error of the slope.
    pub fn stand_err_slope(&self) -> f64 {
        self.stand_error_slope
    }

    /// The chi-squared value: the (weighted, for a weighted fit) residual sum
    /// of squares.
    pub fn chi_squared(&self) -> f64 {
        self.chi_squared
    }

    /// The relative standard deviation, in the non-standard form of Mayr et
    /// al. (2006) that the source uses.
    pub fn rsd(&self) -> f64 {
        self.rsd
    }

    /// Given `x`, compute `y = slope * x + intercept`.
    pub fn compute_point_y(x: f64, slope: f64, intercept: f64) -> f64 {
        slope * x + intercept
    }
}

/// Boost's `quantile(students_t(df), p)` through `statrs`.
///
/// `statrs` inverts the regularised incomplete beta function (AS 64/109)
/// rather than using Boost's Hill/Shaw inversion, so the result is not claimed
/// bit-identical; `docs/RANSAC_SUPPORT.md` records the measured agreement.
fn students_t_quantile(degrees_of_freedom: f64, p: f64) -> Result<f64> {
    let distribution = StudentsT::new(0.0, 1.0, degrees_of_freedom)
        .map_err(|error| Error::InvalidValue(format!("Student's t distribution: {error}")))?;
    // p is in (0, 1) here, which statrs requires (it asserts [0, 1]).
    Ok(distribution.inverse_cdf(p))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn failed_fits_leave_the_previous_fit_untouched() {
        let mut regression = LinearRegression::new();
        regression
            .compute_regression(0.95, &[0.0, 1.0, 2.0], &[1.0, 3.0, 5.0], true)
            .unwrap();
        let before = regression.clone();
        assert!(
            regression
                .compute_regression(0.95, &[1.0, 1.0], &[1.0, 2.0], true)
                .is_err()
        );
        assert!(
            regression
                .compute_regression(0.95, &[1.0, f64::NAN], &[1.0, 2.0], true)
                .is_err()
        );
        assert!(
            regression
                .compute_regression_weighted(0.95, &[1.0, 1.0], &[1.0, 2.0], &[1.0, 1.0], true)
                .is_err()
        );
        assert!(
            regression
                .compute_regression(1.0, &[0.0, 1.0, 2.0], &[1.0, 3.0, 5.0], true)
                .is_err()
        );
        assert_eq!(regression, before);
    }
}
