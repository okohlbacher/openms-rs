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

use super::{
    DataPoint, Extrapolation, InterpolatedModel, Interpolation, InterpolationOptions, bad, finite,
    validate_points,
};
use crate::Result;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LowessOptions {
    /// Fraction of observations per local fit, strictly positive and at most one.
    pub span: f64,
    /// Robustifying passes after the initial fit; at most 64.
    pub iterations: usize,
    /// None uses one percent of the x range; Some requires finite nonnegative delta.
    pub delta: Option<f64>,
    pub interpolation: Interpolation,
    pub extrapolation: Extrapolation,
    pub max_points: usize,
    /// Limit on local sample operations, interpolation, and residual passes.
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
            || self.iterations > 64
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
#[derive(Clone, Debug)]
pub struct LowessModel {
    interpolation: InterpolatedModel,
    options: LowessOptions,
}
impl LowessModel {
    pub fn fit(data: &[DataPoint], options: LowessOptions) -> Result<Self> {
        options.validate()?;
        validate_points(data, options.max_points)?;
        let mut sorted: Vec<_> = data.iter().map(|p| (p.x, p.y)).collect();
        sorted.sort_by(|a, b| a.0.total_cmp(&b.0));
        let x: Vec<_> = sorted.iter().map(|p| p.0).collect();
        let y: Vec<_> = sorted.iter().map(|p| p.1).collect();
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
    pub fn options(&self) -> LowessOptions {
        self.options
    }
    pub fn knots(&self) -> impl Iterator<Item = (f64, f64)> + '_ {
        self.interpolation.knots()
    }
    pub fn apply(&self, x: f64) -> Result<f64> {
        self.interpolation.apply(x)
    }
}

/// Source FastLowessSmoothing: robust local linear regression, tricube distance
/// weights, repeated-x handling and interpolation of delta-skipped fits.
/// Inputs must be finite and sorted by nondecreasing x; signed y values are valid.
/// This helper needs at least two observations; the interpolated model needs
/// three distinct x coordinates after smoothing.
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
