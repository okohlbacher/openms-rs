// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni, OpenMS Rust contributors $

//! Exponentially modified Gaussian fitting with the OpenMS iRprop+ optimizer.
//!
//! Ported from `MATH/MISC/EmgGradientDescent` at revision `7c029e8`. Header,
//! implementation and class test are byte-identical at the current target
//! revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`, so the port carries
//! forward unchanged; `tests/data/emg_provenance.json` records both pins and
//! the verifying hashes.
//!
//! The EMG model combines a Gaussian with an exponential decay and has four
//! parameters: amplitude `h`, mean `mu`, standard deviation `sigma` and the
//! exponential relaxation time `tau` which controls tailing. Fitting uses
//! iRprop+ (Igel and Hüsken, *Improving the Rprop Learning Algorithm*, NC 2000,
//! 115-121), a resilient-backpropagation variant with an independent step size
//! per parameter: a consistent gradient sign accelerates, a sign change halves
//! the step and reverts. The model follows Kalambet, Kozmin, Mikhailova, Nagaev
//! and Tikhonov, "Reconstruction of chromatographic peaks using the
//! exponentially modified Gaussian function", *Journal of Chemometrics* 25
//! (2011) 352-356, including its three numerically distinct `z` regimes.
//!
//! Every optimizer hyper-parameter is hard coded, as in the source: the source
//! exposes only `print_debug`, `max_gd_iter` and `compute_additional_points`
//! through `DefaultParamHandler`.
//!
//! Coordinates are not rescaled: the source initialization and constraints depend
//! on their absolute units, so the same peak in minutes and in seconds fits to
//! different parameters by design. Numerical failures return errors, including
//! overflow in the source's exponential expressions. See `docs/EMG_SUPPORT.md`.

use crate::kernel::{ChromatogramPeak, MSChromatogram, MSSpectrum, Peak1D};
use crate::{Error, Result};
use libm::erfc;
use std::f64::consts::PI;
use std::ops::Range;

/// EMG amplitude, location, Gaussian width and exponential decay parameter.
///
/// The source transports these four values in the order `h`, `mu`, `sigma`,
/// `tau` through a `FloatDataArray` named `emg_parameters`; this port keeps them
/// typed and in `f64`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EmgParameters {
    /// Amplitude: peak height.
    pub h: f64,
    /// Mean: the Gaussian centre position, in the input coordinate's own units.
    pub mu: f64,
    /// Standard deviation: the Gaussian width. Must be positive.
    pub sigma: f64,
    /// Exponential relaxation time, which controls the degree of tailing. Must
    /// be positive.
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
    /// Parameters of the first iteration which achieved the best loss.
    pub parameters: EmgParameters,
    /// Number of evaluated iterations, without the C++ exhaustion off-by-one.
    /// `estimateEmgParameters` returns `iter_idx` after the `while (++iter_idx
    /// <= max_gd_iter_)` test has already incremented it past the limit.
    pub iterations: usize,
    /// One-based iteration which first achieved the reported best loss.
    pub best_iteration: usize,
    /// Mean squared error of the model over the extracted training set, at
    /// `parameters`. This is the source's `Loss_function` on the selected
    /// points, not on the whole input.
    pub loss: f64,
    /// Number of points the source's training-set extraction selected.
    pub training_points: usize,
    /// Scalar model/gradient evaluations during parameter estimation only.
    pub evaluations: usize,
    /// The source's loss-history stopping condition was met.
    pub converged: bool,
}

/// Positions and modelled intensities produced by evaluating the EMG.
///
/// The source writes these through the `out_xs` and `out_ys` out-parameters of
/// `applyEstimatedParameters`. Both vectors always have the same length, and
/// the original input positions appear unchanged; any extrapolated points are
/// prepended or appended, never interleaved.
#[derive(Clone, Debug, PartialEq)]
pub struct EmgCurve {
    /// Positions, ascending, in the input coordinate's own units.
    pub positions: Vec<f64>,
    /// Modelled intensities, one per position, in `f64`.
    pub intensities: Vec<f64>,
}

/// A fitted spectrum together with its estimate and the arrays that were dropped.
#[derive(Clone, Debug, PartialEq)]
pub struct EmgSpectrumFit {
    /// The reconstructed peak. Record metadata is preserved; the peak list is
    /// replaced by the modelled samples, which are usually more numerous than
    /// the input because the truncated side is extrapolated.
    pub spectrum: MSSpectrum,
    /// Parameters and diagnostics of the fit.
    pub estimate: EmgEstimate,
    /// Names of the input data arrays that were not carried over.
    /// Input arrays have no defined aggregation/extrapolation rule.
    pub omitted_arrays: Vec<String>,
}

/// A fitted chromatogram together with its estimate and the arrays that were dropped.
#[derive(Clone, Debug, PartialEq)]
pub struct EmgChromatogramFit {
    /// The reconstructed peak, with record metadata preserved.
    pub chromatogram: MSChromatogram,
    /// Parameters and diagnostics of the fit.
    pub estimate: EmgEstimate,
    /// Names of the input data arrays that were not carried over.
    pub omitted_arrays: Vec<String>,
}

/// Checked source fitting options and explicit per-call resource limits.
///
/// [`Default`] reproduces the source's `getDefaultParameters`: `max_gd_iter` is
/// 100,000 and `compute_additional_points` is `"true"`. The source's third
/// parameter, `print_debug` (0 to 2), has no counterpart; the information it
/// prints is returned in [`EmgEstimate`] instead.
#[derive(Clone, Debug, PartialEq)]
pub struct EmgGradientDescent {
    /// Maximum number of gradient-descent iterations. Must be positive; the
    /// source accepts zero, which makes its loop body unreachable.
    pub max_iterations: usize,
    /// Whether to extrapolate the cutoff side of the peak when applying the
    /// estimated parameters.
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
    /// Check that every configured limit is usable before any work starts.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `max_iterations`, `max_points` or
    /// `max_evaluations` is zero. The source has no equivalent: `max_gd_iter`
    /// declares a minimum of zero in its `Param`, and a zero there silently
    /// returns the `DBL_MAX` sentinels the optimizer initialises its bests with.
    pub fn validate(&self) -> Result<()> {
        if self.max_iterations == 0 || self.max_points == 0 || self.max_evaluations == 0 {
            return Err(bad(
                "EMG iteration, point and evaluation limits must be positive",
            ));
        }
        Ok(())
    }

    /// Run the gradient-descent estimation of the four EMG parameters.
    ///
    /// Source `estimateEmgParameters`, whose four out-parameters `best_h`,
    /// `best_mu`, `best_sigma` and `best_tau` and whose `UInt` return value -
    /// the number of iterations needed - are returned together in
    /// [`EmgEstimate`].
    ///
    /// # Arguments
    ///
    /// * `xs` - positions; at least two, finite and strictly increasing
    /// * `ys` - intensities, one per position and finite
    ///
    /// Estimate from finite intensities and at least two strictly increasing positions.
    /// Signed intensities are retained. The source's initial mean must be positive
    /// because it initializes sigma as one percent of that absolute coordinate.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] when `xs` descends, and
    /// [`Error::InvalidValue`] when a value is not finite, when the slice
    /// lengths differ, when positions repeat, when the point or evaluation
    /// ceiling is exceeded, or when any model, gradient or loss evaluation
    /// leaves the finite range. The source's `extractTrainingSet` throws
    /// `Exception::SizeUnderflow` for fewer than two points, and its optimizer
    /// instead breaks out of the loop on a non-finite parameter or loss and
    /// returns whatever best it had reached.
    pub fn estimate_parameters(&self, xs: &[f64], ys: &[f64]) -> Result<EmgEstimate> {
        self.validate()?;
        self.estimate(xs, ys, &mut Budget::new(self.max_evaluations))
    }

    /// Compute the EMG function on a set of points.
    ///
    /// Source `applyEstimatedParameters`. When
    /// [`compute_additional_points`](Self::compute_additional_points) is set,
    /// the algorithm detects which side of the peak is cut off and extends it.
    ///
    /// # Arguments
    ///
    /// * `xs` - positions, ascending
    /// * `parameters` - amplitude, mean, standard deviation and relaxation time
    ///
    /// Evaluate supplied parameters and optionally extend the truncated side.
    /// With additional points disabled, empty or single-position inputs are valid.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] or [`Error::InvalidValue`] as
    /// [`estimate_parameters`](Self::estimate_parameters) does, and additionally
    /// [`Error::InvalidValue`] when extrapolation is requested for fewer than
    /// two points - the source would divide by `xs.size() - 1` - or when the
    /// generated points would exceed the point ceiling.
    pub fn apply_parameters(&self, xs: &[f64], parameters: EmgParameters) -> Result<EmgCurve> {
        self.validate()?;
        self.apply(xs, parameters, &mut Budget::new(self.max_evaluations))
    }

    /// Fit a spectrum to the EMG peak model, optionally over an m/z sub-range.
    ///
    /// Source `fitEMGPeakModel<MSSpectrum>`. The method recapitulates the actual
    /// peak area of saturated or cut-off peaks and fine tunes well acquired
    /// ones. The output is a reconstruction of the input peak; additional points
    /// are often added so that the boundary intensities match.
    ///
    /// A *cutoff peak* is one whose left and right baseline intensities are not
    /// equal. A *saturated peak* is one whose maximum intensity is lower than
    /// expected because the detector saturated.
    ///
    /// # Arguments
    ///
    /// * `input` - the peak to fit; never modified, including on failure
    /// * `left` - m/z of the first point of interest, or `None` for the first peak
    /// * `right` - m/z of the last point of interest, or `None` for the last peak
    ///
    /// Fit inclusive optional m/z bounds. `Some(0.0)` is a literal bound: the
    /// source uses `0.0` as the "no bound given" sentinel through
    /// `left_pos ? PosBegin(left_pos) : begin()`, so it cannot express a
    /// boundary at zero.
    ///
    /// The returned copy preserves record metadata and omits all input data
    /// arrays. The source instead appends a four-element `FloatDataArray` named
    /// `emg_parameters` holding `h`, `mu`, `sigma` and `tau`; those values are
    /// in [`EmgSpectrumFit::estimate`] here, because an array of four entries
    /// beside a peak list of a different length violates this crate's data-array
    /// alignment invariant.
    ///
    /// # Errors
    ///
    /// As [`estimate_parameters`](Self::estimate_parameters), plus
    /// [`Error::InvalidValue`] when a bound is not finite, when `left` exceeds
    /// `right`, or when a fitted intensity is not representable as the `f32` a
    /// peak stores. The source casts to `float` unchecked.
    pub fn fit_spectrum(
        &self,
        input: &MSSpectrum,
        left: Option<f64>,
        right: Option<f64>,
    ) -> Result<EmgSpectrumFit> {
        self.preflight(input.len())?;
        // Preserve owned acquisition data while bounding the metadata clone;
        // shared processing records require only Arc handle copies here.
        let (mut copy_work, mut copy_bytes) = (50_000_000, 256 * 1024 * 1024);
        input.record_metadata_with_budget(&mut copy_work, &mut copy_bytes)?;
        input.acquisition_with_budget(&mut copy_work, &mut copy_bytes)?;
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

    /// Fit a chromatogram to the EMG peak model, optionally over an RT sub-range.
    ///
    /// Source `fitEMGPeakModel<MSChromatogram>`, the second of the source's two
    /// explicit template instantiations. Behaves exactly as
    /// [`fit_spectrum`](Self::fit_spectrum), which documents the cutoff and
    /// saturation cases, the `emg_parameters` array and the error conditions.
    ///
    /// # Arguments
    ///
    /// * `input` - the peak to fit; never modified, including on failure
    /// * `left` - RT of the first point of interest, or `None` for the first peak
    /// * `right` - RT of the last point of interest, or `None` for the last peak
    ///
    /// Fit inclusive optional retention-time bounds, in the input coordinate
    /// units: the source never converts minutes to seconds, and the two give
    /// different fits because initialization scales with the absolute mean.
    ///
    /// # Errors
    ///
    /// As [`fit_spectrum`](Self::fit_spectrum).
    pub fn fit_chromatogram(
        &self,
        input: &MSChromatogram,
        left: Option<f64>,
        right: Option<f64>,
    ) -> Result<EmgChromatogramFit> {
        self.preflight(input.len())?;
        // Preserve owned acquisition data while bounding the metadata clone;
        // shared processing records require only Arc handle copies here.
        let (mut copy_work, mut copy_bytes) = (50_000_000, 256 * 1024 * 1024);
        input.record_metadata_with_budget(&mut copy_work, &mut copy_bytes)?;
        input.acquisition_with_budget(&mut copy_work, &mut copy_bytes)?;
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
        // Source: computeMuMaxDistance(TrX), over the training positions rather
        // than the input. Training always keeps both endpoints, so for the
        // strictly increasing input this port requires the two spans are equal.
        let radius = finite(mu_max_distance(&tx))?;
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

/// Source `computeMuMaxDistance`: 35% of the span of the given positions.
///
/// Together with the initial mean this bounds how far `mu` may travel during
/// gradient descent. The source takes `std::minmax_element` over the *training*
/// positions, which are not sorted, and returns `0.0` when the container is
/// empty because both returned iterators then equal `end()`.
fn mu_max_distance(xs: &[f64]) -> f64 {
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for &x in xs {
        if x < min {
            min = x;
        }
        if x > max {
            max = x;
        }
    }
    if xs.is_empty() {
        return 0.0;
    }
    (max - min) * 0.35
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
    fn source_mu_max_distance_uses_the_unsorted_span_and_tolerates_an_empty_set() {
        // START_SECTION(double computeMuMaxDistance(const std::vector<double>& xs))
        // The source's own literal vector is not sorted; minmax_element still
        // gives (2, 9), so the answer is (9 - 2) * 0.35.
        let xs = [3.0, 2.0, 4.0, 2.0, 4.0, 5.0, 7.0, 9.0, 3.0];
        assert!((mu_max_distance(&xs) - 2.45).abs() < 1e-14);
        assert_eq!(mu_max_distance(&[]), 0.0); // empty vector case
        assert_eq!(mu_max_distance(&[7.0]), 0.0);
        // Training always keeps both input endpoints, so for the sorted input
        // this port accepts, the training span equals the input span. That is
        // why estimating over the training set is equivalent here.
        let sorted = [1.0, 2.0, 3.0, 4.0, 5.0];
        let (tx, _) = training_set(&sorted, &[1.0, 5.0, 9.0, 5.0, 1.0]).unwrap();
        assert_eq!(mu_max_distance(&tx), mu_max_distance(&sorted));
        assert_eq!(mu_max_distance(&sorted), 4.0 * 0.35);
    }

    #[test]
    fn source_compute_z_selects_the_three_documented_regimes() {
        // START_SECTION(double compute_z(x, mu, sigma, tau))
        // `compute_z` is private in C++ too and is reachable only through the
        // EmgGradientDescent_friend shim; this is that shim's equivalent.
        let p = EmgParameters {
            h: 15_515_900.0,
            mu: 14.3453,
            sigma: 0.0344277,
            tau: 0.188507,
        };
        // The three source literals, for the section's own parameters.
        let close = |a: f64, b: f64| assert!((a - b).abs() <= 1e-9 * b.abs().max(1.0), "{a} {b}");
        close(z(p.mu - 1.0 / 60.0, p), 0.471456263584609);
        close(z(p.mu + 1.0 / 60.0, p), -0.213173439809831);
        close(z(-3_333_333.0, p), 68_463_258.2588395);
        // Each literal selects a different EMG expression: the first is in
        // [0, 6.71e7], the second below zero, the third above the threshold.
        assert!((0.0..=6.71e7).contains(&z(p.mu - 1.0 / 60.0, p)));
        assert!(z(p.mu + 1.0 / 60.0, p) < 0.0);
        assert!(z(-3_333_333.0, p) > 6.71e7);
        // z changes sign at mu + sigma^2 / tau, which is where the source
        // switches from the negative-z expression to the middle one.
        let crossing = p.mu + p.sigma * p.sigma / p.tau;
        assert!(z(crossing, p).abs() < 1e-9);
        assert!(z(crossing - 1e-6, p) > 0.0);
        assert!(z(crossing + 1e-6, p) < 0.0);
        // The threshold is inclusive: z == 6.71e7 still takes the middle
        // expression, which overflows, while the next step takes the
        // asymptotic one and underflows to zero.
        let unit = EmgParameters {
            h: 1.0,
            mu: 0.0,
            sigma: 1.0,
            tau: 1.0,
        };
        assert_eq!(z(1.0, unit), 0.0);
        let at_threshold = 1.0 - 2.0_f64.sqrt() * 6.709e7;
        let past_threshold = 1.0 - 2.0_f64.sqrt() * 6.711e7;
        assert!(z(at_threshold, unit) <= 6.71e7);
        assert!(z(past_threshold, unit) > 6.71e7);
        assert!(model(at_threshold, unit).is_err());
        assert_eq!(model(past_threshold, unit).unwrap(), 0.0);
    }

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
