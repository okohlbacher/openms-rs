// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Least-squares fit of a line through the origin, `Y = c1 x`, accumulated one
//! observation at a time.
//!
//! Ports `ML/REGRESSION/LinearRegressionWithoutIntercept.h` with its `.cpp` at
//! core `bc9cc12`; see `docs/RANSAC_SUPPORT.md`. The source's `@brief` names
//! the model `Y(c, x) = c0 + c1 x`, but the class fits no intercept: the slope
//! is `sum(x y) / sum(x x)`, the least-squares estimate for a line through the
//! origin.

use crate::{Error, Result};

/// Running sums for a least-squares line through the origin.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LinearRegressionWithoutIntercept {
    /// Total variation in x, `sum(x * x)`.
    sum_xx: f64,
    /// Sum of products, `sum(x * y)`.
    sum_xy: f64,
    /// Number of observations.
    n: u64,
}

impl LinearRegressionWithoutIntercept {
    /// An empty data set.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add the observation `(x, y)`: `x` is the independent and `y` the
    /// dependent variable.
    pub fn add_data(&mut self, x: f64, y: f64) {
        self.sum_xx += x * x;
        self.sum_xy += x * y;
        self.n = self.n.saturating_add(1);
    }

    /// Add the observations `(x[i], y[i])` in order.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] when `y` has fewer values than `x`,
    /// and leaves the data set unchanged. The source indexes `y` by the
    /// positions of `x` without a check, reading past a shorter `y`; values of
    /// a longer `y` beyond `x.len()` are ignored, as the source's loop over
    /// `x` ignores them.
    pub fn add_data_slices(&mut self, x: &[f64], y: &[f64]) -> Result<()> {
        if y.len() < x.len() {
            return Err(Error::InvalidRange(format!(
                "regression through the origin: {} x values but {} y values",
                x.len(),
                y.len()
            )));
        }
        for (xi, yi) in x.iter().zip(y) {
            self.add_data(*xi, *yi);
        }
        Ok(())
    }

    /// The slope of the estimated regression line, `sum(x y) / sum(x x)`.
    ///
    /// Returns `None` with fewer than two observations, where the source
    /// returns a quiet NaN ("not enough data"). With two or more observations
    /// the quotient is returned as computed, so all-zero `x` values give the
    /// source's `0 / 0` NaN or an infinity.
    pub fn slope(&self) -> Option<f64> {
        if self.n < 2 {
            None
        } else {
            Some(self.sum_xy / self.sum_xx)
        }
    }

    /// The number of observations added so far.
    ///
    /// Native accessor; the source keeps the count in a private `int`, which
    /// overflows after `2^31 - 1` observations. The count here saturates at
    /// `u64::MAX` instead.
    pub fn len(&self) -> u64 {
        self.n
    }

    /// Whether no observation has been added.
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }
}
