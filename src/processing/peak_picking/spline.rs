// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use crate::{Error, Result};

fn checked(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::InvalidValue(
            "cubic spline arithmetic is not finite".into(),
        ))
    }
}

/// Natural cubic spline with zero second derivatives at both endpoints.
/// Construction follows OpenMS CubicSpline2d's tridiagonal recurrence.
#[derive(Clone, Debug)]
pub struct CubicSpline2d {
    x: Vec<f64>,
    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    d: Vec<f64>,
}
impl CubicSpline2d {
    /// Construct a spline with at most one million knots.
    pub fn new(x: &[f64], y: &[f64]) -> Result<Self> {
        Self::with_max_points(x, y, 1_000_000)
    }
    pub fn with_max_points(x: &[f64], y: &[f64], max_points: usize) -> Result<Self> {
        if x.len() != y.len() || x.len() < 2 || max_points < 2 || x.len() > max_points {
            return Err(Error::InvalidValue(
                "spline needs matching arrays of 2..=max_points knots".into(),
            ));
        }
        if x.iter().chain(y).any(|v| !v.is_finite()) {
            return Err(Error::InvalidValue("spline knots must be finite".into()));
        }
        if x.windows(2).any(|p| p[0] >= p[1]) {
            return Err(Error::InvalidValue(
                "spline coordinates must be strictly increasing".into(),
            ));
        }
        let n = x.len() - 1;
        let h: Vec<_> = x.windows(2).map(|p| p[1] - p[0]).collect();
        if h.iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidValue("spline knot spacing overflows".into()));
        }
        let mut mu = vec![0.0; n];
        let mut z = vec![0.0; n];
        for i in 1..n {
            let span = checked(x[i + 1] - x[i - 1])?;
            let l = checked(2.0 * span - h[i - 1] * mu[i - 1])?;
            mu[i] = checked(h[i] / l)?;
            z[i] = checked(
                (3.0 * (y[i + 1] * h[i - 1] - y[i] * span + y[i - 1] * h[i]) / (h[i - 1] * h[i])
                    - h[i - 1] * z[i - 1])
                    / l,
            )?;
        }
        let mut b = vec![0.0; n];
        let mut c = vec![0.0; n + 1];
        let mut d = vec![0.0; n];
        for j in (0..n).rev() {
            c[j] = checked(z[j] - mu[j] * c[j + 1])?;
            b[j] = checked((y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0)?;
            d[j] = checked((c[j + 1] - c[j]) / (3.0 * h[j]))?;
        }
        Ok(Self {
            x: x.to_vec(),
            a: y[..n].to_vec(),
            b,
            c,
            d,
        })
    }
    pub fn domain(&self) -> (f64, f64) {
        (self.x[0], self.x[self.x.len() - 1])
    }
    fn interval(&self, x: f64) -> Result<usize> {
        let (lo, hi) = self.domain();
        if !x.is_finite() || x < lo || x > hi {
            return Err(Error::InvalidValue("query outside spline domain".into()));
        }
        // At a knot use the interval starting at that knot; at the last knot use the preceding interval.
        Ok(self
            .x
            .partition_point(|&v| v <= x)
            .saturating_sub(1)
            .min(self.a.len() - 1))
    }
    pub fn eval(&self, x: f64) -> Result<f64> {
        let i = self.interval(x)?;
        let dx = x - self.x[i];
        checked(((self.d[i] * dx + self.c[i]) * dx + self.b[i]) * dx + self.a[i])
    }
    /// Evaluate derivative order one, two, or three.
    pub fn derivative(&self, x: f64, order: u8) -> Result<f64> {
        let i = self.interval(x)?;
        let dx = x - self.x[i];
        checked(match order {
            1 => self.b[i] + 2.0 * self.c[i] * dx + 3.0 * self.d[i] * dx * dx,
            2 => 2.0 * self.c[i] + 6.0 * self.d[i] * dx,
            3 => 6.0 * self.d[i],
            _ => {
                return Err(Error::InvalidValue(
                    "spline derivative order must be 1, 2, or 3".into(),
                ));
            }
        })
    }
    /// OpenMS derivative bisection for a peak bracket, assuming a positive
    /// derivative on its left side. This is not a global spline maximizer.
    pub fn peak_maximum(
        &self,
        mut left: f64,
        mut right: f64,
        tolerance: f64,
    ) -> Result<(f64, f64)> {
        self.interval(left)?;
        self.interval(right)?;
        if left > right || !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(Error::InvalidValue(
                "invalid spline maximum bracket or tolerance".into(),
            ));
        }
        for _ in 0..128 {
            let mid = left / 2.0 + right / 2.0;
            let derivative = self.derivative(mid, 1)?;
            if derivative.abs() <= f64::EPSILON || mid == left || mid == right {
                return Ok((mid, self.eval(mid)?));
            }
            if derivative < 0.0 {
                right = mid;
            } else {
                left = mid;
            }
            if right - left <= tolerance {
                let mid = left / 2.0 + right / 2.0;
                return Ok((mid, self.eval(mid)?));
            }
        }
        Err(Error::InvalidValue(
            "spline maximum did not converge within 128 iterations".into(),
        ))
    }
}
