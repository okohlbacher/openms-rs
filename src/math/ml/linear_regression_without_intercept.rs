// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Least-squares fit of a straight line through the origin.
//!
//! Native equivalent of `ML/REGRESSION/LinearRegressionWithoutIntercept.h` and
//! its `.cpp`. See `docs/LINEAR_REGRESSION_WITHOUT_INTERCEPT_SUPPORT.md`.

use crate::{Error, Result};

/// Least-squares fit to the model `y = c1 * x` (source
/// `Math::LinearRegressionWithoutIntercept`).
///
/// The accumulator keeps `sum(x * x)`, `sum(x * y)` and the number of
/// observations; the slope is their quotient once two observations are in.
/// The source's default constructor, copy and destructor map to `Default`,
/// `Clone` and ordinary drop.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LinearRegressionWithoutIntercept {
    /// Total variation in x (source `sum_xx_`).
    sum_xx: f64,
    /// Sum of products (source `sum_xy_`).
    sum_xy: f64,
    /// Number of observations (source `n_`, an `int`; unsigned 64-bit here, so
    /// it cannot overflow).
    n: u64,
}

impl LinearRegressionWithoutIntercept {
    /// An empty data set (source default constructor): the slope is NaN until
    /// two observations are added.
    pub const fn new() -> Self {
        Self {
            sum_xx: 0.0,
            sum_xy: 0.0,
            n: 0,
        }
    }

    /// Adds the observation `(x, y)` (source `addData(double, double)`).
    pub fn add_data(&mut self, x: f64, y: f64) {
        self.sum_xx += x * x;
        self.sum_xy += x * y;
        self.n += 1;
    }

    /// Adds the observations `(x[i], y[i])` for every `i` below `x.len()`, in
    /// order (source `addData(std::vector<double>&, std::vector<double>&)`).
    ///
    /// The source iterates over `x` and reads `y[i]` without a check, so a
    /// shorter `y` reads out of bounds; here it is an error and nothing is
    /// added. Entries of a longer `y` beyond `x.len()` are ignored, as in the
    /// source.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `y` has fewer entries than `x`.
    pub fn add_data_slices(&mut self, x: &[f64], y: &[f64]) -> Result<()> {
        if y.len() < x.len() {
            return Err(Error::InvalidValue(format!(
                "LinearRegressionWithoutIntercept: {} dependent values for {} independent values",
                y.len(),
                x.len()
            )));
        }
        for (&xi, &yi) in x.iter().zip(y) {
            self.add_data(xi, yi);
        }
        Ok(())
    }

    /// The slope of the estimated regression line, `sum(x * y) / sum(x * x)`
    /// (source `getSlope`); NaN with fewer than two observations ("not enough
    /// data"). A data set with `sum(x * x) = 0` divides by zero as the source
    /// does (NaN or an infinity).
    pub fn slope(&self) -> f64 {
        if self.n < 2 {
            return f64::NAN;
        }
        self.sum_xy / self.sum_xx
    }
}
