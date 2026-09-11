// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni, OpenMS Rust contributors $

//! Exponentially modified Gaussian fitting with the OpenMS iRprop+ optimizer.
//!
//! Ported from `MATH/MISC/EmgGradientDescent` at revision `7c029e8`.
//! Coordinates are not rescaled: the source initialization and constraints depend
//! on their absolute units. Numerical failures return errors, including overflow
//! in the source's exponential expressions. See `docs/EMG_SUPPORT.md`.

use crate::kernel::{ChromatogramPeak, MSChromatogram, MSSpectrum, Peak1D};
use crate::{Error, Result};
use libm::erfc;
use std::f64::consts::PI;
use std::ops::Range;

/// EMG amplitude, location, Gaussian width and exponential decay parameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EmgParameters {
    pub h: f64,
    pub mu: f64,
    pub sigma: f64,
    pub tau: f64,
}

impl EmgParameters {
    /// Width and decay must be positive; every parameter must be finite.
    pub fn validate(self) -> Result<()> {
        if [self.h, self.mu, self.sigma, self.tau]
            .iter()
            .any(|v| !v.is_finite())
            || self.sigma <= 0.0
            || self.tau <= 0.0
        {
            return Err(bad(
                "EMG parameters require finite values and positive widths",
            ));
        }
        Ok(())
    }
}

/// Best training loss and the parameters at which it was observed.
#[derive(Clone, Debug, PartialEq)]
pub struct EmgEstimate {
    pub parameters: EmgParameters,
    /// Number of evaluated iterations, without the C++ exhaustion off-by-one.
    pub iterations: usize,
    /// One-based iteration which first achieved the reported best loss.
    pub best_iteration: usize,
    pub loss: f64,
    pub training_points: usize,
    /// Scalar model/gradient evaluations during parameter estimation only.
    pub evaluations: usize,
    /// The source's loss-history stopping condition was met.
    pub converged: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EmgCurve {
    pub positions: Vec<f64>,
    pub intensities: Vec<f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EmgSpectrumFit {
    pub spectrum: MSSpectrum,
    pub estimate: EmgEstimate,
    /// Input arrays have no defined aggregation/extrapolation rule.
    pub omitted_arrays: Vec<String>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EmgChromatogramFit {
    pub chromatogram: MSChromatogram,
    pub estimate: EmgEstimate,
    pub omitted_arrays: Vec<String>,
}

/// Checked source fitting options and explicit per-call resource limits.
#[derive(Clone, Debug, PartialEq)]
pub struct EmgGradientDescent {
    pub max_iterations: usize,
    pub compute_additional_points: bool,
    /// Maximum whole-input and generated point count. Must be positive.
    pub max_points: usize,
    /// Maximum scalar model/gradient evaluations across fitting and application.
    /// One optimizer iteration consumes five evaluations per training point.
    /// Input validation and training extraction are separately bounded by max_points.
    pub max_evaluations: usize,
}

impl Default for EmgGradientDescent {
    fn default() -> Self {
        Self {
            max_iterations: 100_000,
            compute_additional_points: true,
            max_points: 1_000_000,
            max_evaluations: 100_000_000,
        }
    }
}

impl EmgGradientDescent {
    pub fn validate(&self) -> Result<()> {
        if self.max_iterations == 0 || self.max_points == 0 || self.max_evaluations == 0 {
            return Err(bad(
                "EMG iteration, point and evaluation limits must be positive",
            ));
        }
        Ok(())
    }

    /// Estimate from finite intensities and at least two strictly increasing positions.
    /// Signed intensities are retained. The source's initial mean must be positive
    /// because it initializes sigma as one percent of that absolute coordinate.
    pub fn estimate_parameters(&self, xs: &[f64], ys: &[f64]) -> Result<EmgEstimate> {
        self.validate()?;
        self.estimate(xs, ys, &mut Budget::new(self.max_evaluations))
    }

    /// Evaluate supplied parameters and optionally extend the truncated side.
    /// With additional points disabled, empty or single-position inputs are valid.
    pub fn apply_parameters(&self, xs: &[f64], parameters: EmgParameters) -> Result<EmgCurve> {
        self.validate()?;
        self.apply(xs, parameters, &mut Budget::new(self.max_evaluations))
    }

    /// Fit inclusive optional m/z bounds. `Some(0.0)` is a literal bound.
    /// The returned copy preserves record metadata and omits all input data arrays.
    pub fn fit_spectrum(
        &self,
        input: &MSSpectrum,
        left: Option<f64>,
        right: Option<f64>,
    ) -> Result<EmgSpectrumFit> {
        self.preflight(input.len())?;
        // Preserve owned acquisition data while bounding the metadata clone;
        // shared processing records require only Arc handle copies here.
        input.acquisition_with_budget(&mut 50_000_000, &mut (256 * 1024 * 1024))?;
        input.validate()?;
        let xs: Vec<_> = input.peaks.iter().map(|p| p.mz).collect();
        let range = self.selected_range(&xs, left, right)?;
        let ys: Vec<_> = input.peaks[range.clone()]
            .iter()
            .map(|p| f64::from(p.intensity))
            .collect();
        let (curve, estimate) = self.fit(&xs[range], &ys)?;
        let peaks = curve
            .positions
            .iter()
            .zip(&curve.intensities)
            .map(|(&x, &y)| Ok(Peak1D::new(x, intensity(y)?)))
            .collect::<Result<Vec<_>>>()?;
        let omitted_arrays = input
            .float_data_arrays
            .iter()
            .map(|a| a.name.clone())
            .chain(input.integer_data_arrays.iter().map(|a| a.name.clone()))
            .chain(input.string_data_arrays.iter().map(|a| a.name.clone()))
            .collect();
        let mut spectrum = input.clone();
        spectrum.peaks = peaks;
        spectrum.float_data_arrays.clear();
        spectrum.integer_data_arrays.clear();
        spectrum.string_data_arrays.clear();
        spectrum.validate()?;
        Ok(EmgSpectrumFit {
            spectrum,
            estimate,
            omitted_arrays,
        })
    }

    /// Fit inclusive optional retention-time bounds, in the input coordinate units.
    pub fn fit_chromatogram(
        &self,
        input: &MSChromatogram,
        left: Option<f64>,
        right: Option<f64>,
    ) -> Result<EmgChromatogramFit> {
        self.preflight(input.len())?;
        // Preserve owned acquisition data while bounding the metadata clone;
        // shared processing records require only Arc handle copies here.
        input.acquisition_with_budget(&mut 50_000_000, &mut (256 * 1024 * 1024))?;
        input.validate()?;
        let xs: Vec<_> = input.peaks.iter().map(|p| p.rt).collect();
        let range = self.selected_range(&xs, left, right)?;
        let ys: Vec<_> = input.peaks[range.clone()]
            .iter()
            .map(|p| f64::from(p.intensity))
            .collect();
        let (curve, estimate) = self.fit(&xs[range], &ys)?;
        let peaks = curve
            .positions
            .iter()
            .zip(&curve.intensities)
            .map(|(&x, &y)| Ok(ChromatogramPeak::new(x, intensity(y)?)))
            .collect::<Result<Vec<_>>>()?;
        let omitted_arrays = input
            .float_data_arrays
            .iter()
            .map(|a| a.name.clone())
            .chain(input.integer_data_arrays.iter().map(|a| a.name.clone()))
            .chain(input.string_data_arrays.iter().map(|a| a.name.clone()))
            .collect();
        let mut chromatogram = input.clone();
        chromatogram.peaks = peaks;
        chromatogram.float_data_arrays.clear();
        chromatogram.integer_data_arrays.clear();
        chromatogram.string_data_arrays.clear();
        chromatogram.validate()?;
        Ok(EmgChromatogramFit {
            chromatogram,
            estimate,
            omitted_arrays,
        })
    }

    fn preflight(&self, count: usize) -> Result<()> {
        self.validate()?;
        if count > self.max_points {
            return Err(bad("EMG input exceeds point limit"));
        }
        Ok(())
    }

    fn selected_range(
        &self,
        xs: &[f64],
        left: Option<f64>,
        right: Option<f64>,
    ) -> Result<Range<usize>> {
        validate_positions(xs, self.max_points, 0)?;
        if left.into_iter().chain(right).any(|x| !x.is_finite())
            || left.zip(right).is_some_and(|(l, r)| l > r)
        {
            return Err(bad("EMG bounds must be finite and ordered"));
        }
        let start = left.map_or(0, |l| xs.partition_point(|&x| x < l));
        let end = right.map_or(xs.len(), |r| xs.partition_point(|&x| x <= r));
        Ok(start..end)
    }

    fn fit(&self, xs: &[f64], ys: &[f64]) -> Result<(EmgCurve, EmgEstimate)> {
        let mut budget = Budget::new(self.max_evaluations);
        let estimate = self.estimate(xs, ys, &mut budget)?;
        let curve = self.apply(xs, estimate.parameters, &mut budget)?;
        Ok((curve, estimate))
    }

    fn estimate(&self, xs: &[f64], ys: &[f64], budget: &mut Budget) -> Result<EmgEstimate> {
        validate_positions(xs, self.max_points, 2)?;
        if ys.len() != xs.len() || ys.iter().any(|y| !y.is_finite()) {
            return Err(bad("EMG intensities must be finite and aligned"));
        }
        let h_lower = ys.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let mu = initial_mean(xs, ys)?;
        let mut p = EmgParameters {
            h: h_lower,
            mu,
            sigma: mu * 0.01,
            tau: mu * 0.01 * 2.0,
        };
        p.validate()?;
        let (tx, ty) = training_set(xs, ys)?;
        let radius = finite((xs[xs.len() - 1] - xs[0]) * 0.35)?;
        let mu_left = finite(mu - radius)?;
        let mu_right = finite(mu + radius)?;
        let mut previous_gradient = [0.0; 4];
        let mut updates = [0.0; 4];
        let mut rates = [0.0125; 4];
        let mut previous_loss = f64::MAX;
        let mut best: Option<EmgEstimate> = None;
        let mut history = [0.0; 10];
        let mut history_index = 0;
        let start_evaluations = budget.used;
        for iteration in 1..=self.max_iterations {
            p.validate()?;
            budget.charge(
                tx.len()
                    .checked_mul(5)
                    .ok_or_else(|| bad("EMG evaluation count overflow"))?,
            )?;
            let loss = loss(&tx, &ty, p)?;
            if best.as_ref().is_none_or(|b| loss < b.loss) {
                best = Some(EmgEstimate {
                    parameters: p,
                    iterations: iteration,
                    best_iteration: iteration,
                    loss,
                    training_points: tx.len(),
                    evaluations: 0,
                    converged: false,
                });
            }
            let mut gradient = [
                gradient_h(&tx, &ty, p)?,
                gradient_mu(&tx, &ty, p)?,
                gradient_sigma(&tx, &ty, p)?,
                gradient_tau(&tx, &ty, p)?,
            ];
            let best = best
                .as_mut()
                .ok_or_else(|| bad("EMG has no finite estimate"))?;
            best.iterations = iteration;
            best.evaluations = budget.used - start_evaluations;
            if iteration % 50 == 0 {
                history[history_index] = loss;
                history_index = (history_index + 1) % 10;
                let mean = finite(history.iter().sum::<f64>() / 10.0)?;
                let mut squared = 0.0;
                for value in history {
                    squared = finite(squared + pow(value - mean, 2.0))?;
                }
                if sqrt(squared / 10.0) < 1.0 {
                    best.converged = true;
                    break;
                }
            }
            let mut values = [p.h, p.mu, p.sigma, p.tau];
            for i in 0..4 {
                irprop_plus(
                    previous_gradient[i],
                    &mut gradient[i],
                    &mut rates[i],
                    &mut updates[i],
                    &mut values[i],
                    loss,
                    previous_loss,
                );
                finite(values[i])?;
            }
            p = EmgParameters {
                h: h_lower.max(values[0]),
                mu: values[1].clamp(mu_left, mu_right),
                sigma: values[2].clamp(1e-4, 20.0),
                tau: values[3],
            };
            p.tau = p.tau.clamp(p.sigma, p.sigma * 15.0);
            previous_gradient = gradient;
            previous_loss = loss;
        }
        best.ok_or_else(|| bad("EMG has no finite estimate"))
    }

    fn apply(&self, xs: &[f64], p: EmgParameters, budget: &mut Budget) -> Result<EmgCurve> {
        validate_positions(
            xs,
            self.max_points,
            if self.compute_additional_points { 2 } else { 0 },
        )?;
        p.validate()?;
        budget.charge(xs.len())?;
        let mut curve = EmgCurve {
            positions: xs.to_vec(),
            intensities: xs.iter().map(|&x| model(x, p)).collect::<Result<_>>()?,
        };
        if !self.compute_additional_points {
            return Ok(curve);
        }
        let mut step = 0.0;
        for pair in xs.windows(2) {
            step = finite(step + (pair[1] - pair[0]))?;
        }
        step /= (xs.len() - 1) as f64;
        if step <= 0.0 || !step.is_finite() {
            return Err(bad("EMG extrapolation spacing is invalid"));
        }
        let mut apex = 0;
        for i in 1..xs.len() {
            if curve.intensities[i] > curve.intensities[apex] {
                apex = i;
            }
        }
        let last = xs.len() - 1;
        let left = curve.intensities[0] > curve.intensities[last];
        let limit = if left {
            xs[apex] - (xs[last] - xs[apex]) * 3.0
        } else {
            xs[apex] + (xs[apex] - xs[0]) * 3.0
        };
        finite(limit)?;
        let target = curve.intensities[if left { last } else { 0 }];
        let mut x = xs[if left { 0 } else { last }];
        let mut y = curve.intensities[if left { 0 } else { last }];
        let mut extra_x = Vec::new();
        let mut extra_y = Vec::new();
        while y > target && y > 1e-3 {
            let next = finite(if left { x - step } else { x + step })?;
            if (left && next < limit) || (!left && next > limit) {
                break;
            }
            if next == x {
                return Err(bad("EMG extrapolation cannot advance the coordinate"));
            }
            if extra_x.len() >= self.max_points - xs.len() {
                return Err(bad("EMG output exceeds point limit"));
            }
            budget.charge(1)?;
            x = next;
            y = model(x, p)?;
            extra_x.push(x);
            extra_y.push(y);
        }
        if left {
            extra_x.reverse();
            extra_y.reverse();
            extra_x.append(&mut curve.positions);
            extra_y.append(&mut curve.intensities);
            curve.positions = extra_x;
            curve.intensities = extra_y;
        } else {
            curve.positions.append(&mut extra_x);
            curve.intensities.append(&mut extra_y);
        }
        Ok(curve)
    }
}

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad("EMG numerical expression is nonfinite"))
    }
}
fn intensity(value: f64) -> Result<f32> {
    let value = value as f32;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad("EMG fitted intensity exceeds f32 range"))
    }
}
fn validate_positions(xs: &[f64], limit: usize, minimum: usize) -> Result<()> {
    if xs.len() < minimum || xs.len() > limit {
        return Err(bad("EMG coordinate count is outside configured limits"));
    }
    if xs.iter().any(|x| !x.is_finite()) {
        return Err(bad("EMG coordinates must be finite"));
    }
    if xs.windows(2).any(|v| v[0] > v[1]) {
        return Err(Error::UnsortedData);
    }
    if xs.windows(2).any(|v| v[0] == v[1]) {
        return Err(bad("EMG coordinates must be distinct"));
    }
    Ok(())
}
struct Budget {
    used: usize,
    limit: usize,
}
impl Budget {
    fn new(limit: usize) -> Self {
        Self { used: 0, limit }
    }
    fn charge(&mut self, count: usize) -> Result<()> {
        self.used = self
            .used
            .checked_add(count)
            .filter(|&n| n <= self.limit)
            .ok_or_else(|| bad("EMG evaluation limit exceeded"))?;
        Ok(())
    }
}

// Small wrappers retain the source formula grouping and floating-point powers.
fn pow(x: f64, exponent: f64) -> f64 {
    x.powf(exponent)
}
fn exp(x: f64) -> f64 {
    x.exp()
}
fn sqrt(x: f64) -> f64 {
    x.sqrt()
}
fn z(x: f64, p: EmgParameters) -> f64 {
    (1.0 / sqrt(2.0)) * (p.sigma / p.tau - (x - p.mu) / p.sigma)
}

fn initial_mean(xs: &[f64], ys: &[f64]) -> Result<f64> {
    let maximum = ys.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut i = 0;
    let mut j = xs.len() - 1;
    let mut left = xs[0];
    let mut right = xs[j];
    let mut sum = 0.0;
    for percentage in [0.6, 0.65, 0.7, 0.75, 0.8, 0.85] {
        while i < xs.len() - 1 && ys[i] <= maximum * percentage {
            left = xs[i];
            i += 1;
        }
        while j >= 1 && ys[j] <= maximum * percentage {
            right = xs[j];
            j -= 1;
        }
        sum = finite(sum + (left + right) / 2.0)?;
    }
    finite(sum / 6.0)
}

fn training_set(xs: &[f64], ys: &[f64]) -> Result<(Vec<f64>, Vec<f64>)> {
    let threshold = ys.iter().copied().fold(f64::NEG_INFINITY, f64::max) * 0.8;
    let mut indices = vec![0];
    let mut i = 1;
    while i < xs.len() - 1 && ys[i] < threshold {
        indices.push(i);
        i += 1;
    }
    indices.push(xs.len() - 1);
    let mut j = xs.len() - 2;
    while i <= j && ys[j] < threshold {
        indices.push(j);
        j -= 1;
    }
    let mut derivatives = vec![0.0; xs.len() + 1];
    derivatives[0] = 1.0;
    derivatives[xs.len()] = -1.0;
    if i > 1 {
        for k in i - 1..=(j + 1).min(xs.len() - 1) {
            derivatives[k] = finite((ys[k] - ys[k - 1]) / (xs[k] - xs[k - 1]))?;
        }
    }
    let maximum = derivatives[i..j + 2]
        .iter()
        .fold(0.0_f64, |m, d| m.max(d.abs()));
    let threshold = maximum * 0.3;
    while i < xs.len() - 1
        && i <= j
        && derivatives[i] > 0.0
        && (derivatives[i].abs() >= threshold || derivatives[i] / derivatives[i - 1] >= 0.6)
    {
        indices.push(i);
        i += 1;
    }
    while j > 0
        && i <= j
        && derivatives[j + 1] < 0.0
        && (derivatives[j + 1].abs() >= threshold || derivatives[j + 1] / derivatives[j + 2] >= 0.6)
    {
        indices.push(j);
        j -= 1;
    }
    Ok((
        indices.iter().map(|&i| xs[i]).collect(),
        indices.iter().map(|&i| ys[i]).collect(),
    ))
}

fn irprop_plus(
    previous: f64,
    gradient: &mut f64,
    rate: &mut f64,
    update: &mut f64,
    parameter: &mut f64,
    loss: f64,
    previous_loss: f64,
) {
    if previous * *gradient > 0.0 {
        *rate = (*rate * 1.2).min(2000.0);
        *update = -(*gradient / gradient.abs()) * *rate;
        *parameter += *update;
    } else if previous * *gradient < 0.0 {
        *rate = (*rate * 0.5).max(0.0);
        if loss > previous_loss {
            *parameter -= *update;
        }
        *gradient = 0.0;
    } else {
        *update = if *gradient != 0.0 {
            -(*gradient / gradient.abs()) * *rate
        } else {
            -*rate
        };
        *parameter += *update;
    }
}

fn loss(xs: &[f64], ys: &[f64], p: EmgParameters) -> Result<f64> {
    let mut sum = 0.0;
    for (&x, &y) in xs.iter().zip(ys) {
        sum = finite(sum + pow(model(x, p)? - y, 2.0) / xs.len() as f64)?;
    }
    Ok(sum)
}

// Literal source expressions: retain all three branches and per-point normalization.
// Redundant source parentheses make comparison with the C++ expressions easier.
#[allow(unused_parens)]
fn gradient_h(xs: &[f64], ys: &[f64], p: EmgParameters) -> Result<f64> {
    let EmgParameters {
        h,
        mu: u,
        sigma: s,
        tau: t,
    } = p;
    let mut sum = 0.0;
    for (&x, &y) in xs.iter().zip(ys) {
        let z = finite(z(x, p))?;
        let value = if z < 0.0 {
            ((s * exp((pow(s, 2.0) + 2.0 * t * u - 4.0 * t * x) / (2.0 * pow(t, 2.0)))
                * erfc((pow(s, 2.0) + t * (u - x)) / (sqrt(2.0) * s * t))
                * (PI
                    * h
                    * s
                    * exp((pow(s, 2.0) + 2.0 * t * u) / (2.0 * pow(t, 2.0)))
                    * erfc((pow(s, 2.0) + t * (u - x)) / (sqrt(2.0) * s * t))
                    - sqrt(2.0 * PI) * t * y * exp(x / t)))
                / pow(t, 2.0))
                / xs.len() as f64
        } else if z <= 6.71e7 {
            ((sqrt(2.0 * PI)
                * s
                * exp(1.0 / 2.0 * pow((s / t - (x - u) / s), 2.0)
                    - pow((x - u), 2.0) / (2.0 * pow(s, 2.0)))
                * erfc((s / t - (x - u) / s) / sqrt(2.0))
                * ((sqrt(PI / 2.0)
                    * h
                    * s
                    * exp(1.0 / 2.0 * pow((s / t - (x - u) / s), 2.0)
                        - pow((x - u), 2.0) / (2.0 * pow(s, 2.0)))
                    * erfc((s / t - (x - u) / s) / sqrt(2.0)))
                    / t
                    - y))
                / t)
                / xs.len() as f64
        } else {
            ((2.0
                * exp(-pow((x - u), 2.0) / (2.0 * pow(s, 2.0)))
                * ((h * exp(-pow((x - u), 2.0) / (2.0 * pow(s, 2.0))))
                    / (1.0 - (t * (x - u)) / pow(s, 2.0))
                    - y))
                / (1.0 - (t * (x - u)) / pow(s, 2.0)))
                / xs.len() as f64
        };
        sum = finite(sum + finite(value)?)?;
    }
    Ok(sum)
}
#[allow(unused_parens)]
fn gradient_mu(xs: &[f64], ys: &[f64], p: EmgParameters) -> Result<f64> {
    let EmgParameters {
        h,
        mu: u,
        sigma: s,
        tau: t,
    } = p;
    let mut sum = 0.0;
    for (&x, &y) in xs.iter().zip(ys) {
        let z = finite(z(x, p))?;
        let value = if z < 0.0 {
            (2.0 * ((sqrt(PI / 2.0)
                * h
                * s
                * exp(pow(s, 2.0) / (2.0 * pow(t, 2.0)) - (x - u) / t)
                * erfc((s / t - (x - u) / s) / sqrt(2.0)))
                / pow(t, 2.0)
                - (h * exp(pow(s, 2.0) / (2.0 * pow(t, 2.0))
                    - 1.0 / 2.0 * pow((s / t - (x - u) / s), 2.0)
                    - (x - u) / t))
                    / t)
                * ((sqrt(PI / 2.0)
                    * h
                    * s
                    * exp(pow(s, 2.0) / (2.0 * pow(t, 2.0)) - (x - u) / t)
                    * erfc((s / t - (x - u) / s) / sqrt(2.0)))
                    / t
                    - y))
                / xs.len() as f64
        } else if z <= 6.71e7 {
            (2.0 * ((sqrt(PI / 2.0)
                * h
                * s
                * exp(1.0 / 2.0 * pow((s / t - (x - u) / s), 2.0)
                    - pow((x - u), 2.0) / (2.0 * pow(s, 2.0)))
                * ((x - u) / pow(s, 2.0) + (s / t - (x - u) / s) / s)
                * erfc((s / t - (x - u) / s) / sqrt(2.0)))
                / t
                - (h * exp(-pow((x - u), 2.0) / (2.0 * pow(s, 2.0)))) / t)
                * ((sqrt(PI / 2.0)
                    * h
                    * s
                    * exp(1.0 / 2.0 * pow((s / t - (x - u) / s), 2.0)
                        - pow((x - u), 2.0) / (2.0 * pow(s, 2.0)))
                    * erfc((s / t - (x - u) / s) / sqrt(2.0)))
                    / t
                    - y))
                / xs.len() as f64
        } else {
            (2.0 * ((h * (x - u) * exp(-pow((x - u), 2.0) / (2.0 * pow(s, 2.0))))
                / (pow(s, 2.0) * (1.0 - (t * (x - u)) / pow(s, 2.0)))
                - (h * t * exp(-pow((x - u), 2.0) / (2.0 * pow(s, 2.0))))
                    / (pow(s, 2.0) * pow((1.0 - (t * (x - u)) / pow(s, 2.0)), 2.0)))
                * ((h * exp(-pow((x - u), 2.0) / (2.0 * pow(s, 2.0))))
                    / (1.0 - (t * (x - u)) / pow(s, 2.0))
                    - y))
                / xs.len() as f64
        };
        sum = finite(sum + finite(value)?)?;
    }
    Ok(sum)
}
#[allow(unused_parens)]
fn gradient_sigma(xs: &[f64], ys: &[f64], p: EmgParameters) -> Result<f64> {
    let EmgParameters {
        h,
        mu: u,
        sigma: s,
        tau: t,
    } = p;
    let mut sum = 0.0;
    for (&x, &y) in xs.iter().zip(ys) {
        let z = finite(z(x, p))?;
        let value = if z < 0.0 {
            (2.0 * ((sqrt(PI / 2.0)
                * h
                * exp(pow(s, 2.0) / (2.0 * pow(t, 2.0)) - (x - u) / t)
                * erfc((s / t - (x - u) / s) / sqrt(2.0)))
                / t
                + (sqrt(PI / 2.0)
                    * h
                    * pow(s, 2.0)
                    * exp(pow(s, 2.0) / (2.0 * pow(t, 2.0)) - (x - u) / t)
                    * erfc((s / t - (x - u) / s) / sqrt(2.0)))
                    / pow(t, 3.0)
                - (h * s
                    * exp(pow(s, 2.0) / (2.0 * pow(t, 2.0))
                        - 1.0 / 2.0 * pow((s / t - (x - u) / s), 2.0)
                        - (x - u) / t)
                    * ((x - u) / pow(s, 2.0) + 1.0 / t))
                    / t)
                * ((sqrt(PI / 2.0)
                    * h
                    * s
                    * exp(pow(s, 2.0) / (2.0 * pow(t, 2.0)) - (x - u) / t)
                    * erfc((s / t - (x - u) / s) / sqrt(2.0)))
                    / t
                    - y))
                / xs.len() as f64
        } else if z <= 6.71e7 {
            (2.0 * ((sqrt(PI / 2.0)
                * h
                * exp(1.0 / 2.0 * pow((s / t - (x - u) / s), 2.0)
                    - pow((x - u), 2.0) / (2.0 * pow(s, 2.0)))
                * erfc((s / t - (x - u) / s) / sqrt(2.0)))
                / t
                + (sqrt(PI / 2.0)
                    * h
                    * s
                    * exp(1.0 / 2.0 * pow((s / t - (x - u) / s), 2.0)
                        - pow((x - u), 2.0) / (2.0 * pow(s, 2.0)))
                    * (pow((x - u), 2.0) / pow(s, 3.0)
                        + ((x - u) / pow(s, 2.0) + 1.0 / t) * (s / t - (x - u) / s))
                    * erfc((s / t - (x - u) / s) / sqrt(2.0)))
                    / t
                - (h * s
                    * exp(-pow((x - u), 2.0) / (2.0 * pow(s, 2.0)))
                    * ((x - u) / pow(s, 2.0) + 1.0 / t))
                    / t)
                * ((sqrt(PI / 2.0)
                    * h
                    * s
                    * exp(1.0 / 2.0 * pow((s / t - (x - u) / s), 2.0)
                        - pow((x - u), 2.0) / (2.0 * pow(s, 2.0)))
                    * erfc((s / t - (x - u) / s) / sqrt(2.0)))
                    / t
                    - y))
                / xs.len() as f64
        } else {
            (2.0 * ((h * pow((x - u), 2.0) * exp(-pow((x - u), 2.0) / (2.0 * pow(s, 2.0))))
                / (pow(s, 3.0) * (1.0 - (t * (x - u)) / pow(s, 2.0)))
                - (2.0 * h * t * (x - u) * exp(-pow((x - u), 2.0) / (2.0 * pow(s, 2.0))))
                    / (pow(s, 3.0) * pow((1.0 - (t * (x - u)) / pow(s, 2.0)), 2.0)))
                * ((h * exp(-pow(x - u, 2.0) / (2.0 * pow(s, 2.0))))
                    / (1.0 - (t * (x - u)) / pow(s, 2.0))
                    - y))
                / xs.len() as f64
        };
        sum = finite(sum + finite(value)?)?;
    }
    Ok(sum)
}
#[allow(unused_parens)]
fn gradient_tau(xs: &[f64], ys: &[f64], p: EmgParameters) -> Result<f64> {
    let EmgParameters {
        h,
        mu: u,
        sigma: s,
        tau: t,
    } = p;
    let mut sum = 0.0;
    for (&x, &y) in xs.iter().zip(ys) {
        let z = finite(z(x, p))?;
        let value = if z < 0.0 {
            (2.0 * (-(sqrt(PI / 2.0)
                * h
                * s
                * exp(pow(s, 2.0) / (2.0 * pow(t, 2.0)) - (x - u) / t)
                * erfc((s / t - (x - u) / s) / sqrt(2.0)))
                / pow(t, 2.0)
                + (sqrt(PI / 2.0)
                    * h
                    * s
                    * exp(pow(s, 2.0) / (2.0 * pow(t, 2.0)) - (x - u) / t)
                    * ((x - u) / pow(t, 2.0) - pow(s, 2.0) / pow(t, 3.0))
                    * erfc((s / t - (x - u) / s) / sqrt(2.0)))
                    / t
                + (h * pow(s, 2.0)
                    * exp(pow(s, 2.0) / (2.0 * pow(t, 2.0))
                        - 1.0 / 2.0 * pow((s / t - (x - u) / s), 2.0)
                        - (x - u) / t))
                    / pow(t, 3.0))
                * ((sqrt(PI / 2.0)
                    * h
                    * s
                    * exp(pow(s, 2.0) / (2.0 * pow(t, 2.0)) - (x - u) / t)
                    * erfc((s / t - (x - u) / s) / sqrt(2.0)))
                    / t
                    - y))
                / xs.len() as f64
        } else if z <= 6.71e7 {
            (2.0 * (-(sqrt(PI / 2.0)
                * h
                * pow(s, 2.0)
                * exp(1.0 / 2.0 * pow((s / t - (x - u) / s), 2.0)
                    - pow(x - u, 2.0) / (2.0 * pow(s, 2.0)))
                * (s / t - (x - u) / s)
                * erfc((s / t - (x - u) / s) / sqrt(2.0)))
                / pow(t, 3.0)
                - (sqrt(PI / 2.0)
                    * h
                    * s
                    * exp(1.0 / 2.0 * pow((s / t - (x - u) / s), 2.0)
                        - pow(x - u, 2.0) / (2.0 * pow(s, 2.0)))
                    * erfc((s / t - (x - u) / s) / sqrt(2.0)))
                    / pow(t, 2.0)
                + (h * pow(s, 2.0) * exp(-pow(x - u, 2.0) / (2.0 * pow(s, 2.0)))) / pow(t, 3.0))
                * ((sqrt(PI / 2.0)
                    * h
                    * s
                    * exp(1.0 / 2.0 * pow((s / t - (x - u) / s), 2.0)
                        - pow(x - u, 2.0) / (2.0 * pow(s, 2.0)))
                    * erfc((s / t - (x - u) / s) / sqrt(2.0)))
                    / t
                    - y))
                / xs.len() as f64
        } else {
            ((2.0
                * h
                * (x - u)
                * exp(-pow(x - u, 2.0) / (2.0 * pow(s, 2.0)))
                * ((h * exp(-pow(x - u, 2.0) / (2.0 * pow(s, 2.0))))
                    / (1.0 - (t * (x - u)) / pow(s, 2.0))
                    - y))
                / (pow(s, 2.0) * pow((1.0 - (t * (x - u)) / pow(s, 2.0)), 2.0)))
                / xs.len() as f64
        };
        sum = finite(sum + finite(value)?)?;
    }
    Ok(sum)
}
#[allow(unused_parens)]
fn model(x: f64, p: EmgParameters) -> Result<f64> {
    let EmgParameters {
        h,
        mu: u,
        sigma: s,
        tau: t,
    } = p;
    let z = finite(z(x, p))?;
    finite(if z < 0.0 {
        ((h * s) / t)
            * sqrt(PI / 2.0)
            * exp((1.0 / 2.0) * (pow(s / t, 2.0)) - (x - u) / t)
            * erfc((1.0 / sqrt(2.0)) * (s / t - (x - u) / s))
    } else if z <= 6.71e7 {
        h * exp(-(1.0 / 2.0) * pow(((x - u) / s), 2.0))
            * (s / t)
            * sqrt(PI / 2.0)
            * exp(pow((1.0 / sqrt(2.0) * (s / t - (x - u) / s)), 2.0))
            * erfc(1.0 / sqrt(2.0) * (s / t - (x - u) / s))
    } else {
        (h * exp(-(1.0 / 2.0) * (pow(((x - u) / s), 2.0))))
            / (1.0 - (((x - u) * t) / (pow(s, 2.0))))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn training_preserves_source_collection_order_and_plateau_exclusion() {
        let xs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let ys = [0.0, 1.0, 5.0, 10.0, 10.0, 10.0, 5.0, 1.0, 0.0];
        let (x, y) = training_set(&xs, &ys).unwrap();
        assert_eq!(x, [1.0, 2.0, 3.0, 9.0, 8.0, 7.0, 4.0, 6.0]);
        assert_eq!(y, [0.0, 1.0, 5.0, 0.0, 1.0, 5.0, 10.0, 10.0]);
        assert_eq!(initial_mean(&xs, &ys).unwrap(), 5.0);
        assert_eq!(
            training_set(&[1.0, 2.0], &[1.0, 2.0]).unwrap(),
            (vec![1.0, 2.0], vec![1.0, 2.0])
        );
    }

    #[test]
    fn source_irprop_updates_include_zero_gradient_and_rollback() {
        for (
            current,
            loss,
            expected_gradient,
            expected_rate,
            expected_update,
            expected_parameter,
        ) in [
            (20.0, 13.0, 20.0, 4.8, -4.8, 855.2),
            (-20.0, 13.0, 0.0, 2.0, 0.5, 860.0),
            (-20.0, 15.0, 0.0, 2.0, 0.5, 859.5),
            (0.0, 13.0, 0.0, 4.0, -4.0, 856.0),
        ] {
            let (mut gradient, mut rate, mut update, mut value) = (current, 4.0, 0.5, 860.0);
            irprop_plus(
                10.0,
                &mut gradient,
                &mut rate,
                &mut update,
                &mut value,
                loss,
                14.0,
            );
            assert_eq!(
                (gradient, rate, update, value),
                (
                    expected_gradient,
                    expected_rate,
                    expected_update,
                    expected_parameter
                )
            );
        }
    }

    #[test]
    fn analytical_gradients_match_independent_central_differences() {
        let p = EmgParameters {
            h: 4.0,
            mu: 2.0,
            sigma: 0.8,
            tau: 1.2,
        };
        // Exercise the two nonzero model branches separately, then together.
        for xs in [
            &[0.7, 1.0, 2.0][..],
            &[3.0, 4.0, 5.0][..],
            &[0.7, 2.0, 3.0, 4.0][..],
        ] {
            let ys: Vec<_> = xs.iter().map(|&x| x * 0.3).collect();
            let analytical = [
                gradient_h(xs, &ys, p).unwrap(),
                gradient_mu(xs, &ys, p).unwrap(),
                gradient_sigma(xs, &ys, p).unwrap(),
                gradient_tau(xs, &ys, p).unwrap(),
            ];
            for i in 0..4 {
                let mut values = [p.h, p.mu, p.sigma, p.tau];
                let step = values[i].abs() * 1e-5;
                values[i] += step;
                let upper = EmgParameters {
                    h: values[0],
                    mu: values[1],
                    sigma: values[2],
                    tau: values[3],
                };
                values[i] -= 2.0 * step;
                let lower = EmgParameters {
                    h: values[0],
                    mu: values[1],
                    sigma: values[2],
                    tau: values[3],
                };
                let numerical =
                    (loss(xs, &ys, upper).unwrap() - loss(xs, &ys, lower).unwrap()) / (2.0 * step);
                assert!(
                    (analytical[i] - numerical).abs() < 1e-8 * numerical.abs().max(1.0),
                    "gradient {i}: {} vs {numerical}",
                    analytical[i]
                );
            }
        }
        let far = [-1e9];
        assert_eq!(gradient_h(&far, &[0.0], p).unwrap(), 0.0);
        assert_eq!(gradient_mu(&far, &[0.0], p).unwrap(), 0.0);
        assert_eq!(gradient_sigma(&far, &[0.0], p).unwrap(), 0.0);
        assert_eq!(gradient_tau(&far, &[0.0], p).unwrap(), 0.0);
    }
}
