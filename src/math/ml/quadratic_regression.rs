// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Least-squares fit of a quadratic, `y = a + b x + c x^2`.
//!
//! Ports `ML/REGRESSION/QuadraticRegression.h` with `QuadraticRegression.cpp`
//! at core `bc9cc12`; see `docs/RANSAC_SUPPORT.md`. The source solves the 3x3
//! (weighted) normal equations through `gte::LinearSystem`, which inverts the
//! matrix by cofactors; that inverse and the matrix-vector product are
//! transcribed operation for operation.

use crate::{Error, Result};

/// A quadratic `y = a + b x + c x^2` fitted by (weighted) least squares.
///
/// Before the first successful fit the coefficients and the chi-squared sum
/// read back as `0.0`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct QuadraticRegression {
    a: f64,
    b: f64,
    c: f64,
    chi_squared: f64,
}

impl QuadraticRegression {
    /// A regression whose coefficients and chi-squared sum are `0.0`; it
    /// produces meaningful results only after the first successful fit.
    pub fn new() -> Self {
        Self::default()
    }

    /// Fit `y = a + b x + c x^2` to unweighted points `(x[i], y[i])`.
    ///
    /// Equivalent to [`compute_regression_weighted`](Self::compute_regression_weighted)
    /// with a weight of `1` for every point, which is how the source computes
    /// it. On success the coefficients and [`chi_squared`](Self::chi_squared)
    /// are updated.
    ///
    /// # Errors
    ///
    /// As [`compute_regression_weighted`](Self::compute_regression_weighted).
    pub fn compute_regression(&mut self, x: &[f64], y: &[f64]) -> Result<()> {
        let weights = vec![1.0; x.len()];
        self.compute_regression_weighted(x, y, &weights)
    }

    /// Fit `y = a + b x + c x^2` to points `(x[i], y[i])` weighted by `w[i]`.
    ///
    /// Solves the weighted normal equations
    /// `[Sw Swx Swxx; Swx Swxx Swxxx; Swxx Swxxx Swxxxx] [a; b; c] = [Swy; Swxy; Swxxy]`.
    /// On success the coefficients and [`chi_squared`](Self::chi_squared) are
    /// updated.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] starting with
    /// `UnableToFit-QuadraticRegression` when the system is singular — its
    /// determinant is exactly zero, as with fewer than three distinct `x`
    /// values — as the source throws `Exception::UnableToFit`. Also for any
    /// non-finite coordinate or weight, or a solution that overflows: the
    /// source checks neither and stores NaN or infinite coefficients. Returns
    /// [`Error::InvalidRange`] when the slices differ in length; the source
    /// requires `y` and `w` to be at least as long as `x` and does not check.
    /// On error `self` is unchanged, as in the source.
    pub fn compute_regression_weighted(&mut self, x: &[f64], y: &[f64], w: &[f64]) -> Result<()> {
        if x.len() != y.len() || x.len() != w.len() {
            return Err(Error::InvalidRange(format!(
                "quadratic regression: {} x values, {} y values and {} weights",
                x.len(),
                y.len(),
                w.len()
            )));
        }
        if x.iter().chain(y).chain(w).any(|value| !value.is_finite()) {
            return Err(unable_to_fit("non-finite coordinate or weight"));
        }
        let mut sum_x = 0.0;
        let mut sum_xx = 0.0;
        let mut sum_xxx = 0.0;
        let mut sum_xxxx = 0.0;
        let mut sum_y = 0.0;
        let mut sum_xy = 0.0;
        let mut sum_xxy = 0.0;
        let mut sum_w = 0.0;
        for ((xi, yi), weight) in x.iter().zip(y).zip(w) {
            sum_x += weight * xi;
            sum_xx += weight * xi * xi;
            sum_xxx += weight * xi * xi * xi;
            sum_xxxx += weight * xi * xi * xi * xi;
            sum_y += weight * yi;
            sum_xy += weight * xi * yi;
            sum_xxy += weight * xi * xi * yi;
            sum_w += weight;
        }
        let matrix = [
            [sum_w, sum_x, sum_xx],
            [sum_x, sum_xx, sum_xxx],
            [sum_xx, sum_xxx, sum_xxxx],
        ];
        let inverse = inverse_3x3(&matrix)
            .ok_or_else(|| unable_to_fit("could not fit a linear model to the data"))?;
        let rhs = [sum_y, sum_xy, sum_xxy];
        let mut solution = [0.0; 3];
        for (row, value) in inverse.iter().zip(solution.iter_mut()) {
            // gte's Matrix * Vector starts each row sum at zero.
            let mut sum = 0.0;
            for (m, b) in row.iter().zip(&rhs) {
                sum += m * b;
            }
            *value = sum;
        }
        if solution.iter().any(|value| !value.is_finite()) {
            return Err(unable_to_fit("the solution overflows"));
        }
        let [a, b, c] = solution;
        let mut chi_squared = 0.0;
        for ((xi, yi), weight) in x.iter().zip(y).zip(w) {
            let residual = yi - a - b * xi - c * xi * xi;
            chi_squared += weight * (residual * residual);
        }
        self.a = a;
        self.b = b;
        self.c = c;
        self.chi_squared = chi_squared;
        Ok(())
    }

    /// The fitted polynomial at `x`: `a + b x + c x^2` with the stored
    /// coefficients.
    pub fn eval(&self, x: f64) -> f64 {
        Self::eval_with(self.a, self.b, self.c, x)
    }

    /// A quadratic polynomial with externally supplied coefficients at `x`:
    /// `a + b x + c x x`, the constant term `a`, the linear coefficient `b` and
    /// the quadratic coefficient `c`.
    ///
    /// The source overloads `eval` with this static form.
    pub fn eval_with(a: f64, b: f64, c: f64, x: f64) -> f64 {
        a + b * x + c * x * x
    }

    /// The constant term `a` of the last successful fit, or `0.0`.
    pub fn a(&self) -> f64 {
        self.a
    }

    /// The linear coefficient `b` of the last successful fit, or `0.0`.
    pub fn b(&self) -> f64 {
        self.b
    }

    /// The quadratic coefficient `c` of the last successful fit, or `0.0`.
    pub fn c(&self) -> f64 {
        self.c
    }

    /// The weighted residual sum of squares of the last successful fit,
    /// `sum(w (y - a - b x - c x x)^2)`, or `0.0`.
    pub fn chi_squared(&self) -> f64 {
        self.chi_squared
    }
}

fn unable_to_fit(reason: &str) -> Error {
    Error::InvalidValue(format!("UnableToFit-QuadraticRegression: {reason}"))
}

/// `gte::Inverse(Matrix3x3)`: the cofactor inverse, `None` for a zero
/// determinant.
fn inverse_3x3(m: &[[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let c00 = m[1][1] * m[2][2] - m[1][2] * m[2][1];
    let c10 = m[1][2] * m[2][0] - m[1][0] * m[2][2];
    let c20 = m[1][0] * m[2][1] - m[1][1] * m[2][0];
    let det = m[0][0] * c00 + m[0][1] * c10 + m[0][2] * c20;
    if det == 0.0 {
        return None;
    }
    let inv_det = 1.0 / det;
    Some([
        [
            c00 * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            c10 * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            c20 * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ])
}
