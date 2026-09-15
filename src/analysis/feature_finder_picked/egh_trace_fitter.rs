// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The exponential-Gaussian hybrid retention-time model of the picked feature
//! finder (`FEATUREFINDER/EGHTraceFitter.h`).
//!
//! An `EGHTraceFitter` fits the exponential-Gaussian hybrid (EGH) of Lan and
//! Jorgenson to the mass traces of one feature candidate:
//!
//! ```text
//! f(t) = H * exp(-(t - t_R)^2 / (2 sigma^2 + tau (t - t_R)))   where 2 sigma^2 + tau (t - t_R) > 0
//! f(t) = 0                                                     elsewhere
//! ```
//!
//! Lan K, Jorgenson JW. *A hybrid of exponential and gaussian functions as a
//! simple model of asymmetric chromatographic peaks.* Journal of Chromatography
//! A 915 (1-2), 1-13. `FeatureFinderAlgorithmPicked` uses it when
//! `feature:rt_shape` is `asymmetric`; the default `symmetric` uses the Gaussian
//! model of `GaussTraceFitter`. The source marks the class experimental: "Needs
//! further testing on real data", and its class test exercises the EGH only as a
//! replacement for the Gaussian, on a symmetric peak.
//!
//! The model implements the shared
//! [`TraceFitter`](crate::analysis::feature_finder_picked::trace_fitter::TraceFitter)
//! contract. The source's public nested functor `EGHTraceFunctor` becomes
//! [`EGHTraceFunctor`](crate::analysis::feature_finder_picked::egh_trace_fitter::EGHTraceFunctor),
//! and the two protected hooks the source's own tests reach through derived
//! classes, `setInitialParameters_` and `getOptimizedParameters_`, together with
//! `getAlphaBoundaries_`, are public inherent methods here so that a start point
//! can be inspected and a model rebuilt from stored parameters.
//!
//! # Preserved source arithmetic
//!
//! Every expression keeps the source's operand order, so that a single
//! evaluation (a residual, a Jacobian entry, a start point or a query) follows
//! the source to the last place in which the mathematical library agrees with
//! the oracle's. Fitted parameters are compared with the C++ fit within a
//! tolerance, not bit for bit (see the support document):
//!
//! - The residual uses the signed `sigma`; the Jacobian uses `|sigma|`.
//! - Where `2 sigma^2 + tau (t - t_R) <= 0` the model is `0`, without the
//!   baseline, and every Jacobian entry of that row is `0`.
//! - Weighting multiplies a trace's residuals and Jacobian rows by its
//!   theoretical intensity itself, not by a square root.
//! - The start point always smooths the intensity profile with the zero-padded
//!   five-point running sum and has no guard for `alpha >= 1`, unlike the
//!   Gaussian: a flat profile yields a NaN `sigma` and a non-finite `tau`, and the
//!   fit then ends at once with those values (see the support document).
//! - A start `tau` of exactly zero is replaced by `f64::EPSILON`.
//! - The retention-time bounds are the positions where the model reaches
//!   `0.043937` of its height, the FWHM uses relative height `0.5`, and the area
//!   is the Lan-Jorgenson approximation with the source's seven coefficients and
//!   the factor `0.6266571`.
//!
//! `exp`, `log`, `sqrt` and `atan` come from the `libm` crate, not from the
//! platform's C library, so the results are the same on every platform up to the
//! sign and payload of a NaN. `GaussTraceFitter` calls the platform library
//! instead; which choice both fitters share is the integrator's decision
//! (`docs/TRACE_FITTER_SUPPORT.md`, "`exp` and `log` across platforms").
//!
//! # Shared trace-fitter helpers
//!
//! The start point takes the profile, the smoothing and the half-height walks
//! from the shared
//! [`initial_shape`](crate::analysis::feature_finder_picked::trace_fitter::initial_shape)
//! with
//! [`ProfileSmoothing::Always`](crate::analysis::feature_finder_picked::trace_fitter::ProfileSmoothing::Always).
//! The fit runs the shared Levenberg-Marquardt driver
//! [`optimize`](crate::analysis::feature_finder_picked::trace_fitter::optimize),
//! the source's `TraceFitter::optimize_`, and `computeTheoretical`, the error
//! of a refused fit and the gnuplot numbers come from the same module.
//!
//! **Solver fidelity.** `optimize` calls
//! [`minimize`](crate::math::fitters::levenberg_marquardt::minimize), which is
//! not yet bit-faithful to the executed Eigen solver. Beyond the recorded
//! fixtures a fit's path can depart from Eigen's at its first trial step, so
//! fitted parameters, the status and the number of evaluations can differ from
//! the C++ well beyond the 1e-9 that the EGH fixtures meet. The gap was measured
//! for `GaussTraceFitter` in `docs/TRACE_FITTER_SUPPORT.md` ("Known gap: solver
//! fidelity beyond the fixtures"), where the departure arises inside `minimize`
//! and not in the functor; the EGH functor goes through the same driver and
//! solver, so the caveat applies to EGH fits as well, although they have not
//! been measured beyond these fixtures. The root cause is in
//! `src/math/fitters/levenberg_marquardt.rs`, under investigation in lane B3b.
//!
//! The source is serial, and so is this module. Its work is bounded by the
//! ceilings of `MassTraces::intensity_profile` and of `optimize`, both checked
//! before anything proportional to the input is allocated.
//!
//! See `docs/EGH_TRACE_FITTER_SUPPORT.md` for the API mapping, the native
//! differences and the evidence.

use crate::analysis::feature_finder_picked::helper_structs::{MassTrace, MassTraces};
use crate::analysis::feature_finder_picked::trace_fitter::{
    FEWER_RESIDUALS_THAN_PARAMETERS, ProfileSmoothing, TraceFitter, TraceFitterParams,
    compute_theoretical, initial_shape, optimize, stream_number, unable_to_fit,
};
use crate::{Error, Result};

/// Relative height at which the retention-time bounds are taken: the source
/// literal in `getOptimizedParameters_`, which its comment calls "conceptually
/// equal to 2.5 sigma" of a Gaussian.
const SIGMA_5_BOUND_ALPHA: f64 = 0.043937;

/// Relative height of the full width at half maximum.
const FWHM_ALPHA: f64 = 0.5;

/// The source's approximation of `sqrt(pi / 8)` in `getArea`.
const AREA_SIGMA_FACTOR: f64 = 0.6266571;

/// The start point of a fit, as the source's protected
/// `EGHTraceFitter::setInitialParameters_` leaves it in the fitter's members.
///
/// Returned by [`EGHTraceFitter::initial_parameters`]. The first four fields are
/// the initial parameter vector in the functor's order: see
/// [`EGHInitialParameters::to_vector`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EGHInitialParameters {
    /// Smoothed maximum of the intensity profile minus the traces' baseline:
    /// source `height_`.
    pub height: f64,
    /// Retention time of the smoothed maximum, in seconds: source `apex_rt_`.
    pub apex_rt: f64,
    /// Width estimate `sqrt(-0.5 / log(alpha) * B * A)`: source `sigma_`. NaN
    /// when `alpha >= 1` makes the radicand negative or undefined.
    pub sigma: f64,
    /// Asymmetry estimate `-1 / log(alpha) * (B - A)`, or `f64::EPSILON` when
    /// that is zero: source `tau_`.
    pub tau: f64,
    /// Retention-time span of the intensity profile, last minus first entry, in
    /// seconds: source `region_rt_span_`. The fit keeps it for
    /// [`TraceFitter::check_maximal_rt_span`].
    pub region_rt_span: f64,
}

impl EGHInitialParameters {
    /// The initial parameter vector `[height, apex_rt, sigma, tau]`, in the order
    /// of [`EGHTraceFunctor`] and of the source's `x_init`.
    pub fn to_vector(&self) -> [f64; EGHTraceFitter::NUM_PARAMS] {
        [self.height, self.apex_rt, self.sigma, self.tau]
    }
}

/// Residuals and analytic Jacobian of the EGH least-squares problem: source
/// `EGHTraceFitter::EGHTraceFunctor`.
///
/// The parameter vector is `[H, t_R, sigma, tau]`. Rows follow the traces in
/// order and, within a trace, its peaks in order, so there are
/// [`MassTraces::peak_count`] residuals. Row `k` for peak `(rt, I)` of a trace
/// with theoretical intensity `theo` is
///
/// ```text
/// t = rt - t_R,   d = 2 sigma sigma + tau t
/// residual = (baseline + theo H exp(-t^2 / d) - I) w     when d > 0
/// residual = (0 - I) w                                   otherwise
/// ```
///
/// with `w = theo` when weighted and `1` otherwise, and `I` the `f32` peak
/// intensity promoted to `f64`. The Jacobian columns, evaluated with
/// `sigma = |x[2]|` (the source comments "must be non-negative!") and all zero
/// where `d <= 0`, are
///
/// ```text
/// dH     = theo exp(-t^2/d) w
/// dt_R   = theo H exp(-t^2/d) ((4 sigma sigma + tau t) t) / (d d) w
/// dsigma = theo H exp(-t^2/d) 4 sigma t^2 / (d d) w
/// dtau   = theo H exp(-t^2/d) t t^2 / (d d) w
/// ```
///
/// For a non-negative `sigma` these are the analytic derivatives of the model
/// term. For a negative `sigma` the `sigma` column keeps the source's `|sigma|`
/// and so has the opposite sign of the residual's derivative; it is preserved,
/// not corrected, because it steers the Levenberg-Marquardt path. Where the
/// denominator is not positive the residual also drops the baseline, while
/// [`TraceFitter::gnuplot_formula`] adds it.
///
/// The source constructor takes the number of parameters (always
/// [`EGHTraceFitter::NUM_PARAMS`]) and a pointer to the protected `ModelData`
/// bundle of the traces and the weighting flag; this one borrows the traces and
/// takes the flag. The source's `operator()` and `df` return `0` and never fail;
/// [`Self::residuals`] and [`Self::jacobian`] return `Result` only because their
/// output slices are caller-sized.
#[derive(Clone, Copy, Debug)]
pub struct EGHTraceFunctor<'a> {
    traces: &'a MassTraces,
    weighted: bool,
    values: usize,
}

impl<'a> EGHTraceFunctor<'a> {
    /// A functor over `traces`, weighting each trace's rows by its theoretical
    /// intensity when `weighted` is true.
    ///
    /// The number of residuals is taken from the traces now, as the source
    /// constructor stores `getPeakCount()`.
    pub fn new(traces: &'a MassTraces, weighted: bool) -> Self {
        Self {
            traces,
            weighted,
            values: traces.peak_count(),
        }
    }

    /// Number of parameters: source `inputs()`, always
    /// [`EGHTraceFitter::NUM_PARAMS`].
    pub fn inputs(&self) -> usize {
        EGHTraceFitter::NUM_PARAMS
    }

    /// Number of residuals, one per peak: source `values()`.
    pub fn values(&self) -> usize {
        self.values
    }

    /// Fills `fvec` with the residuals at the parameter vector `x`: source
    /// `operator()(const double* x, double* fvec)`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `x` does not hold
    /// [`Self::inputs`] values or `fvec` does not hold [`Self::values`] values,
    /// leaving `fvec` unchanged. The source writes through raw pointers of the
    /// sizes it was constructed with and cannot check.
    pub fn residuals(&self, x: &[f64], fvec: &mut [f64]) -> Result<()> {
        let x = parameter_vector(x)?;
        if fvec.len() != self.values {
            return Err(Error::InvalidValue(format!(
                "EGHTraceFunctor residuals need {} slots, got {}",
                self.values,
                fvec.len()
            )));
        }
        self.fill_residuals(x, fvec);
        Ok(())
    }

    /// Fills `jacobian`, stored column-major with [`Self::values`] rows and
    /// [`Self::inputs`] columns as the source's `Eigen::Map`, with the Jacobian
    /// at the parameter vector `x`: source `df(const double* x, double* J)`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `x` does not hold
    /// [`Self::inputs`] values or `jacobian` does not hold `values * inputs`
    /// values, leaving `jacobian` unchanged.
    pub fn jacobian(&self, x: &[f64], jacobian: &mut [f64]) -> Result<()> {
        let x = parameter_vector(x)?;
        let expected = self
            .values
            .checked_mul(EGHTraceFitter::NUM_PARAMS)
            .ok_or_else(|| Error::InvalidValue("EGHTraceFunctor Jacobian size overflows".into()))?;
        if jacobian.len() != expected {
            return Err(Error::InvalidValue(format!(
                "EGHTraceFunctor Jacobian needs {expected} slots, got {}",
                jacobian.len()
            )));
        }
        let rows = self.values;
        self.fill_jacobian(x, |row, column, value| {
            if let Some(slot) = jacobian.get_mut(column * rows + row) {
                *slot = value;
            }
        });
        Ok(())
    }

    /// The residual loop of `operator()`, writing as many rows as both the
    /// traces and `fvec` hold.
    fn fill_residuals(&self, x: [f64; EGHTraceFitter::NUM_PARAMS], fvec: &mut [f64]) {
        let [height, t_r, sigma, tau] = x;
        let baseline = self.traces.baseline;
        let mut slots = fvec.iter_mut();
        for trace in self.traces.iter() {
            let weight = if self.weighted {
                trace.theoretical_int
            } else {
                1.0
            };
            for peak in &trace.peaks {
                let Some(slot) = slots.next() else {
                    return;
                };
                let t_diff = peak.rt - t_r;
                let t_diff2 = t_diff * t_diff;
                let denominator = 2.0 * sigma * sigma + tau * t_diff;
                let fegh = if denominator > 0.0 {
                    baseline + trace.theoretical_int * height * libm::exp(-t_diff2 / denominator)
                } else {
                    0.0
                };
                *slot = (fegh - f64::from(peak.intensity)) * weight;
            }
        }
    }

    /// The Jacobian loop of `df`, handing each entry to `put(row, column,
    /// value)`.
    fn fill_jacobian<P>(&self, x: [f64; EGHTraceFitter::NUM_PARAMS], mut put: P)
    where
        P: FnMut(usize, usize, f64),
    {
        let height = x[0];
        let t_r = x[1];
        let sigma = x[2].abs();
        let tau = x[3];
        let mut row = 0usize;
        for trace in self.traces.iter() {
            let weight = if self.weighted {
                trace.theoretical_int
            } else {
                1.0
            };
            let theo = trace.theoretical_int;
            for peak in &trace.peaks {
                if row >= self.values {
                    return;
                }
                let t_diff = peak.rt - t_r;
                let t_diff2 = t_diff * t_diff;
                let denominator = 2.0 * sigma * sigma + tau * t_diff;
                let (derivative_h, derivative_t_r, derivative_sigma, derivative_tau) =
                    if denominator > 0.0 {
                        let exp1 = libm::exp(-t_diff2 / denominator);
                        let denominator2 = denominator * denominator;
                        (
                            theo * exp1,
                            theo * height * exp1 * ((4.0 * sigma * sigma + tau * t_diff) * t_diff)
                                / denominator2,
                            theo * height * exp1 * 4.0 * sigma * t_diff2 / denominator2,
                            theo * height * exp1 * t_diff * t_diff2 / denominator2,
                        )
                    } else {
                        (0.0, 0.0, 0.0, 0.0)
                    };
                put(row, 0, derivative_h * weight);
                put(row, 1, derivative_t_r * weight);
                put(row, 2, derivative_sigma * weight);
                put(row, 3, derivative_tau * weight);
                row += 1;
            }
        }
    }
}

/// The exponential-Gaussian hybrid retention-time model: source
/// `EGHTraceFitter`.
///
/// Create it with [`Self::new`] (the source defaults) or
/// [`Self::with_parameters`], call [`TraceFitter::fit`], then read the fitted
/// model through the [`TraceFitter`] queries and the inherent [`Self::tau`] and
/// [`Self::sigma`], which `FeatureFinderAlgorithmPicked` stores as the feature
/// meta values `EGH_tau` and `EGH_sigma` (with `EGH_height` from
/// [`TraceFitter::height`]).
///
/// Before the first successful fit every model value is `0.0`, where the source
/// leaves its `double` members uninitialised; the queries are meaningful only
/// after a fit or [`Self::set_optimized_parameters`]. `Clone` is the source's
/// copy constructor and assignment operator, which copy the members and call
/// the no-op `updateMembers_`.
#[derive(Clone, Debug, PartialEq)]
pub struct EGHTraceFitter {
    parameters: TraceFitterParams,
    apex_rt: f64,
    height: f64,
    sigma: f64,
    tau: f64,
    sigma_5_bound: (f64, f64),
    region_rt_span: f64,
}

impl Default for EGHTraceFitter {
    /// The same as [`EGHTraceFitter::new`].
    fn default() -> Self {
        Self::new()
    }
}

impl EGHTraceFitter {
    /// Number of model parameters, `[H, t_R, sigma, tau]`: source
    /// `NUM_PARAMS_`.
    pub const NUM_PARAMS: usize = 4;

    /// Coefficients of the polynomial in `phi = atan(|tau| / |sigma|)` that
    /// scales the peak area, from table 1 of the Lan and Jorgenson paper: source
    /// `EPSILON_COEFS_`.
    pub const EPSILON_COEFS: [f64; 7] = [
        4.0, -6.293724, 9.232834, -11.342910, 9.123978, -4.173753, 0.827797,
    ];

    /// A fitter with the source defaults: `max_iteration` 500 and unweighted
    /// traces, [`TraceFitterParams::default`]. Source `EGHTraceFitter()`.
    pub fn new() -> Self {
        Self::with_parameters(TraceFitterParams::default())
    }

    /// A fitter with the given parameters: source `EGHTraceFitter()` followed
    /// by `setParameters`.
    pub fn with_parameters(parameters: TraceFitterParams) -> Self {
        Self {
            parameters,
            apex_rt: 0.0,
            height: 0.0,
            sigma: 0.0,
            tau: 0.0,
            sigma_5_bound: (0.0, 0.0),
            region_rt_span: 0.0,
        }
    }

    /// The fitted exponential time constant `tau`, in seconds: source `getTau`.
    ///
    /// Positive values describe tailing, negative values fronting peaks. The
    /// value is the optimiser's result as is; the source does not constrain it.
    pub fn tau(&self) -> f64 {
        self.tau
    }

    /// The fitted Gaussian width `sigma`, in seconds: source `getSigma`.
    ///
    /// Returned with the sign the optimiser left, as in the source: the model
    /// only uses `sigma^2` and `|sigma|`, so a negative value describes the same
    /// peak.
    pub fn sigma(&self) -> f64 {
        self.sigma
    }

    /// The start point a fit of `traces` begins from: source protected
    /// `setInitialParameters_`.
    ///
    /// Takes the shared [`initial_shape`] with [`ProfileSmoothing::Always`]:
    /// the traces' intensity profile, smoothed with a five-point running sum
    /// over zero padding whatever its length, the first strict maximum as the
    /// apex, and the walks outwards
    /// while the smoothed intensity stays above half the height. They give the
    /// left and right half-height positions `A = apex - left` and
    /// `B = right - apex` and their smoothed heights. With
    /// `alpha = (left_height + right_height) * 0.5 / height`, `tau = -1 /
    /// log(alpha) * (B - A)` (replaced by `f64::EPSILON` when exactly zero) and
    /// `sigma = sqrt(-0.5 / log(alpha) * B * A)`.
    ///
    /// Unlike the Gaussian model, nothing guards short profiles or `alpha >= 1`:
    /// a flat or single-scan profile gives `alpha = 1`, hence a NaN `sigma` and
    /// an infinite or NaN `tau`, and a baseline above the profile's edges gives
    /// `alpha > 1` and a NaN `sigma`. These values are returned as the source
    /// computes them.
    ///
    /// The source writes the members and logs every intermediate value at debug
    /// level; this returns the values and logs nothing.
    ///
    /// # Errors
    ///
    /// Propagates the errors of `initial_shape`: the peak and merge-step
    /// ceilings of [`MassTraces::intensity_profile`], a NaN retention time the
    /// merge cannot place, and [`Error::InvalidValue`] when the profile is
    /// empty because the traces hold no peak. The source then reads
    /// `smoothed[0]` of an empty vector, which is undefined behaviour.
    pub fn initial_parameters(traces: &MassTraces) -> Result<EGHInitialParameters> {
        let shape = initial_shape(traces, ProfileSmoothing::Always)?;
        let height = shape.height;
        let apex_rt = shape.apex_rt;
        let a = apex_rt - shape.left_rt;
        let b = shape.right_rt - apex_rt;
        let alpha = (shape.left_height + shape.right_height) * 0.5 / height;
        let log_alpha = libm::log(alpha);
        let mut tau = -1.0 / log_alpha * (b - a);
        if tau == 0.0 {
            tau = f64::EPSILON;
        }
        let sigma = libm::sqrt(-0.5 / log_alpha * b * a);

        Ok(EGHInitialParameters {
            height,
            apex_rt,
            sigma,
            tau,
            region_rt_span: shape.region_rt_span,
        })
    }

    /// The ordered pair of retention times where the model reaches `alpha`
    /// times its height, in seconds: source protected `getAlphaBoundaries_`.
    ///
    /// Solves equations A.2 and A.3 of the Lan and Jorgenson paper: with
    /// `L = log(alpha)` and `s = sqrt((L tau)^2 / 4 - 2 L sigma^2)`, the bounds are
    /// the apex plus the smaller and the larger of `-L tau / 2 + s` and
    /// `-L tau / 2 - s`. The smaller and larger are chosen as `std::min` and
    /// `std::max` choose them: the first operand is kept unless the second
    /// compares strictly smaller (or larger), so a NaN second operand yields the
    /// first and `+0.0` is kept against `-0.0`.
    ///
    /// `alpha` is meaningful in `(0, 1]`. Outside it the source formula runs
    /// unchanged: `alpha = 0` gives infinite or NaN bounds and a negative
    /// `alpha` NaN bounds.
    pub fn alpha_boundaries(&self, alpha: f64) -> (f64, f64) {
        let l = libm::log(alpha);
        let s =
            libm::sqrt((l * self.tau) * (l * self.tau) / 4.0 - 2.0 * l * self.sigma * self.sigma);
        // The source's `-1 * (L * tau_)`: multiplying by -1 is an exact
        // negation, so `-(l * tau)` has the same bits (a NaN's sign aside).
        let s1 = (-(l * self.tau) / 2.0) + s;
        let s2 = (-(l * self.tau) / 2.0) - s;
        let smaller = if s2 < s1 { s2 } else { s1 };
        let larger = if s1 < s2 { s2 } else { s1 };
        (self.apex_rt + smaller, self.apex_rt + larger)
    }

    /// Sets the model to the parameter vector `[H, t_R, sigma, tau]` and
    /// recomputes the retention-time bounds: source protected
    /// `getOptimizedParameters_`, which `optimize_` calls with the optimised
    /// vector.
    ///
    /// Public here so that a model can be rebuilt from stored parameters, for
    /// example a feature's `EGH_height`, `EGH_tau` and `EGH_sigma` with its
    /// retention time. The region span used by
    /// [`TraceFitter::check_maximal_rt_span`] is not part of the vector and
    /// keeps its value (`0.0` on a fitter that has not been fitted); the source
    /// leaves it uninitialised on a new fitter.
    pub fn set_optimized_parameters(&mut self, x: [f64; Self::NUM_PARAMS]) {
        self.height = x[0];
        self.apex_rt = x[1];
        self.sigma = x[2];
        self.tau = x[3];
        self.sigma_5_bound = self.alpha_boundaries(SIGMA_5_BOUND_ALPHA);
    }
}

impl TraceFitter for EGHTraceFitter {
    fn parameters(&self) -> &TraceFitterParams {
        &self.parameters
    }

    fn set_parameters(&mut self, parameters: TraceFitterParams) {
        self.parameters = parameters;
    }

    /// Fits the EGH to `traces`: source `EGHTraceFitter::fit`.
    ///
    /// Derives the start point with [`EGHTraceFitter::initial_parameters`],
    /// minimises the residuals of [`EGHTraceFunctor`] with its analytic
    /// Jacobian through the shared driver [`optimize`] under the evaluation
    /// budget [`TraceFitterParams::max_iteration`], and
    /// sets the model from the result with
    /// [`EGHTraceFitter::set_optimized_parameters`]. Every Levenberg-Marquardt
    /// status after `ImproperInputParameters` is accepted, including an
    /// exhausted budget and a start point at which the gradient already
    /// vanishes, as it does for the NaN start of a flat profile.
    ///
    /// The solver is not yet bit-faithful to the executed Eigen beyond the
    /// recorded fixtures; see "Solver fidelity" in the module documentation.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] starting with `UnableToFit-FinalSet` when
    /// the traces hold fewer than four peaks (message "Skipping feature, we
    /// always expect N>=p"), or when `max_iteration` is zero or negative, which
    /// the driver refuses as improper input ("Could not fit the gaussian to the
    /// data: Error 0", the source's wording). The start point's errors
    /// propagate, except that traces without any peak give the
    /// `UnableToFit-FinalSet` error, where the source's start point is undefined
    /// behaviour. The resource ceilings of `optimize` propagate as well: the
    /// solver's point and byte ceilings, and a fit that exhausts the work
    /// ceiling `MAX_RESIDUAL_WORK` before its configured budget.
    ///
    /// On any error the fitter is unchanged. The source has already written the
    /// start point into its members when `optimize_` throws, and keeps the
    /// previous retention-time bounds. `GaussTraceFitter` keeps its start values
    /// after a refused fit, as its source does.
    fn fit(&mut self, traces: &MassTraces) -> Result<()> {
        let values = traces.peak_count();
        if values == 0 {
            return Err(unable_to_fit(FEWER_RESIDUALS_THAN_PARAMETERS));
        }
        let initial = Self::initial_parameters(traces)?;
        let mut x = initial.to_vector();
        let functor = EGHTraceFunctor::new(traces, self.parameters.weighted);
        optimize(
            &mut x,
            values,
            |parameters, fvec| functor.fill_residuals(solver_vector(parameters), fvec),
            |parameters, jacobian| {
                functor.fill_jacobian(solver_vector(parameters), |row, column, value| {
                    jacobian.set(row, column, value);
                });
            },
            &self.parameters,
        )?;
        self.set_optimized_parameters(x);
        self.region_rt_span = initial.region_rt_span;
        Ok(())
    }

    /// The position where the model falls to `0.043937` of its height before
    /// the apex: source `getLowerRTBound`, the first of the bounds
    /// [`EGHTraceFitter::set_optimized_parameters`] stores.
    fn lower_rt_bound(&self) -> f64 {
        self.sigma_5_bound.0
    }

    /// The position where the model falls to `0.043937` of its height after
    /// the apex: source `getUpperRTBound`.
    fn upper_rt_bound(&self) -> f64 {
        self.sigma_5_bound.1
    }

    /// The fitted height `H`: source `getHeight`.
    fn height(&self) -> f64 {
        self.height
    }

    /// The fitted apex retention time `t_R`: source `getCenter`.
    fn center(&self) -> f64 {
        self.apex_rt
    }

    /// The distance between the positions at half height: source `getFWHM`,
    /// from [`EGHTraceFitter::alpha_boundaries`] at `0.5`.
    fn fwhm(&self) -> f64 {
        let (lower, upper) = self.alpha_boundaries(FWHM_ALPHA);
        upper - lower
    }

    /// Equation 12 of the Lan and Jorgenson paper at `rt`: source `getValue`.
    ///
    /// `H * exp(-(rt - t_R)^2 / (2 sigma^2 + tau (rt - t_R)))` where the
    /// denominator is positive, `0.0` elsewhere. The baseline is not added.
    fn value(&self, rt: f64) -> f64 {
        let t_diff = rt - self.apex_rt;
        let denominator = 2.0 * self.sigma * self.sigma + self.tau * t_diff;
        if denominator > 0.0 {
            self.height * libm::exp(-t_diff * t_diff / denominator)
        } else {
            0.0
        }
    }

    /// Equation 21 of the Lan and Jorgenson paper: source `getArea`.
    ///
    /// `H * (|sigma| * 0.6266571 + |tau|) * epsilon(phi)`, where
    /// `phi = atan(|tau| / |sigma|)` and `epsilon` is the polynomial with
    /// [`EGHTraceFitter::EPSILON_COEFS`], accumulated in ascending powers of
    /// `phi`. For `tau = 0` this is `4 * 0.6266571 * H |sigma|`, close to the
    /// Gaussian's `sqrt(2 pi) H sigma`. `sigma = tau = 0` gives NaN.
    fn area(&self) -> f64 {
        let abs_tau = self.tau.abs();
        let abs_sigma = self.sigma.abs();
        let phi = libm::atan(abs_tau / abs_sigma);
        let mut epsilon = Self::EPSILON_COEFS[0];
        let mut phi_pow = phi;
        for coefficient in &Self::EPSILON_COEFS[1..] {
            epsilon += phi_pow * coefficient;
            phi_pow *= phi;
        }
        self.height * (abs_sigma * AREA_SIGMA_FACTOR + abs_tau) * epsilon
    }

    /// Whether the model is narrower than `min_rt_span` times the span between
    /// its bounds requires: source `checkMinimalRTSpan`, true when
    /// `rt_bounds.1 - rt_bounds.0 < min_rt_span * (upper - lower)`.
    ///
    /// A NaN bound makes the comparison, and the result, false.
    fn check_minimal_rt_span(&self, rt_bounds: (f64, f64), min_rt_span: f64) -> bool {
        (rt_bounds.1 - rt_bounds.0) < min_rt_span * (self.sigma_5_bound.1 - self.sigma_5_bound.0)
    }

    /// Whether the model's bounds span more than `max_rt_span` times the
    /// retention-time span of the profile it was fitted to: source
    /// `checkMaximalRTSpan`, true when
    /// `upper - lower > max_rt_span * region_rt_span`.
    ///
    /// A NaN bound makes the result false.
    fn check_maximal_rt_span(&self, max_rt_span: f64) -> bool {
        (self.sigma_5_bound.1 - self.sigma_5_bound.0) > max_rt_span * self.region_rt_span
    }

    /// A gnuplot definition of the model for `trace`: source
    /// `getGnuplotFormula`.
    ///
    /// The text is
    /// `N(x)= B + (((S + T * (x - C )) > 0) ? A * exp(-1 * (x - C)**2 / ( S + T * (x - C ))) : 0)`
    /// with `N` the function name, `B` the baseline, `S = 2 * sigma * sigma`,
    /// `T = tau`, `C = rt_shift + t_R` and `A = theoretical_int * H`, each number
    /// written as a default C++ stream writes a `double` (precision 6, `%g`
    /// style), through the shared [`stream_number`].
    /// The source writes the name with `StringUtils::toStr(char)`, one byte; a
    /// non-ASCII Rust `char` is written as its UTF-8 bytes.
    fn gnuplot_formula(
        &self,
        trace: &MassTrace,
        function_name: char,
        baseline: f64,
        rt_shift: f64,
    ) -> String {
        let g = stream_number;
        let two_sigma_squared = g(2.0 * self.sigma * self.sigma);
        let tau = g(self.tau);
        let center = g(rt_shift + self.apex_rt);
        let mut formula = String::new();
        formula.push(function_name);
        formula.push_str("(x)= ");
        formula.push_str(&g(baseline));
        formula.push_str(" + (((");
        formula.push_str(&two_sigma_squared);
        formula.push_str(" + ");
        formula.push_str(&tau);
        formula.push_str(" * (x - ");
        formula.push_str(&center);
        formula.push_str(" )) > 0) ? ");
        formula.push_str(&g(trace.theoretical_int * self.height));
        formula.push_str(" * exp(-1 * (x - ");
        formula.push_str(&center);
        formula.push_str(")**2 / ( ");
        formula.push_str(&two_sigma_squared);
        formula.push_str(" + ");
        formula.push_str(&tau);
        formula.push_str(" * (x - ");
        formula.push_str(&center);
        formula.push_str(" ))) : 0)");
        formula
    }

    /// `trace.theoretical_int * value(trace.peaks[k].rt)`: source
    /// `TraceFitter::computeTheoretical`, through the shared [`compute_theoretical`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `k` is not an index into
    /// `trace.peaks`, where the source indexes without a check.
    fn compute_theoretical(&self, trace: &MassTrace, k: usize) -> Result<f64> {
        compute_theoretical(self, trace, k)
    }
}

/// The four-element parameter vector of a caller-supplied slice.
fn parameter_vector(x: &[f64]) -> Result<[f64; EGHTraceFitter::NUM_PARAMS]> {
    <[f64; EGHTraceFitter::NUM_PARAMS]>::try_from(x).map_err(|_| {
        Error::InvalidValue(format!(
            "EGHTraceFunctor needs {} parameters, got {}",
            EGHTraceFitter::NUM_PARAMS,
            x.len()
        ))
    })
}

/// The parameter vector the solver hands back, which always has the length of
/// the vector it was given; any other length yields NaN parameters instead of
/// a panic.
fn solver_vector(x: &[f64]) -> [f64; EGHTraceFitter::NUM_PARAMS] {
    parameter_vector(x).unwrap_or([f64::NAN; EGHTraceFitter::NUM_PARAMS])
}
