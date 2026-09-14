// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The Gaussian retention-time model of the picked feature finder
//! (`FEATUREFINDER/GaussTraceFitter.h`).
//!
//! `FeatureFinderAlgorithmPicked` fits this model to the mass traces of every
//! feature candidate when `feature:rt_shape` is `symmetric`, its default. The
//! model is
//!
//! ```text
//! f(rt) = height * exp(-0.5 * (rt - x0)^2 / sigma^2)
//! ```
//!
//! scaled per trace by the trace's theoretical intensity and shifted by the
//! traces' baseline. [`GaussTraceFitter`](crate::analysis::feature_finder_picked::gauss_trace_fitter::GaussTraceFitter)
//! derives start values from the summed intensity profile, fits `height`, `x0`
//! and `sigma` with the shared Levenberg-Marquardt driver
//! [`optimize`](crate::analysis::feature_finder_picked::trace_fitter::optimize),
//! and answers the queries of the
//! [`TraceFitter`](crate::analysis::feature_finder_picked::trace_fitter::TraceFitter)
//! trait from the fitted values.
//! [`GaussTraceFunctor`](crate::analysis::feature_finder_picked::gauss_trace_fitter::GaussTraceFunctor)
//! is the source's protected residual and Jacobian functor.
//!
//! This is not `crate::math::fitters::gauss`, the port of
//! `MATH/STATISTICS/GaussFitter`: that fits a one-dimensional distribution with
//! a different objective and initialisation.
//!
//! # Preserved source arithmetic
//!
//! - The sigma column of the Jacobian carries the source's factor `0.125`
//!   (`GaussTraceFitter.cpp:199`), not the true derivative's `1`. It shapes the
//!   Levenberg-Marquardt path, so the fitted parameters depend on it.
//! - The residual computes the exponent as `c_fac * d^2` with
//!   `c_fac = -0.5 / sigma^2`; [`value`] computes `(-0.5 * d^2) / sigma^2`.
//!   The two can differ in the last bit and are kept apart.
//! - Residual and Jacobian rows run trace by trace, then peak by peak.
//! - Weighting multiplies a trace's residuals and derivatives by its theoretical
//!   intensity itself, not by a function of it.
//! - The start values: the five-point running mean is skipped for profiles of
//!   three or fewer retention times; the first strict maximum wins; an `alpha`
//!   of 1 or more gives `sigma = 1.0`; otherwise
//!   `sigma = (delta_x * 0.5) / sqrt(-2 * ln(alpha))`.
//! - The constants `2.35482` (FWHM), `2.506628` (area) and `2.5` (bounds) are
//!   the source's truncated literals, not `2 * sqrt(2 ln 2)` and
//!   `sqrt(2 * pi)`.
//! - The fitted `sigma` is stored as its absolute value.
//!
//! # Platform `exp` and `log`
//!
//! `exp` and `log` are the platform C library's, as in the source, through
//! `f64::exp` and `f64::ln`. They are the only platform-dependent operations on
//! the fit path: `sqrt` is correctly rounded everywhere, and the solver calls
//! nothing else from the C library. A fit is serial and repeats bit for bit on
//! one platform, but its last bits depend on the C library, so results are not
//! bit-identical across platforms:
//!
//! - glibc 2.39 (Linux x86-64) and Apple libm (macOS arm64, where the oracle
//!   ran) return different last bits at 70 of the 37,080 distinct `exp`
//!   arguments the tests reach; the 25 `log` arguments agree. At the 3,972
//!   residual points the executed C++ recorded, both reproduce the oracle bit
//!   for bit.
//! - Fits that pass through such an argument differ between the two. For
//!   example, `start.trailing_max` deviates from the oracle by 1.72e-3 on Linux
//!   and 1.49e-3 on macOS, and of the 19,572 FeatureFinderCentroided_1 values
//!   16,699 are bit-identical on Linux and 16,696 on macOS.
//! - With Apple libm's values substituted for every `exp` and `log` on Linux,
//!   the Linux run reproduces the macOS run, so the C library is the whole
//!   difference. The acceptance criteria hold on both platforms.
//!
//! Two platform-independent choices were measured and not adopted:
//!
//! - The `libm` crate's `exp` differs from both platform libraries by one unit
//!   in the last place at 362 of the recorded residual points. That breaks the
//!   1e-14 residual criterion (up to 9.1e-12 relative) and moves fitted
//!   parameters by up to 1.13e-9.
//! - A correctly rounded `exp` and `log`, simulated by table lookup, meets every
//!   criterion and would give the same bits on every platform. Neither platform
//!   library is correctly rounded, though (glibc misrounds 31 of the arguments,
//!   Apple libm 69), so it matches the oracle in fewer last bits: 16,658 of the
//!   FeatureFinderCentroided_1 values.
//!
//! Choosing between the platform library and a correctly rounded one is left
//! to the integrator; the Gaussian and EGH fitters must make the same choice.
//!
//! # Known gap: solver fidelity
//!
//! The start values, the residuals and the Jacobian are the source's bit for
//! bit wherever they were executed. The fit that follows is not, in general.
//! The review of this package ran 79 further inputs through the product-SDK
//! `fit` and an Eigen replica of `optimize_` on macOS arm64. Against them the
//! port, on macOS arm64 (Linux x86-64 in parentheses where it differs), left
//! Eigen's residual path in 72 cases, 58 of them at the first trial step,
//! from identical start vectors, residuals and Jacobians. 21 fits missed
//! 1e-9: by 1.0e-9 to 3.2e-4 after natural termination, and by up to 4.9e-3
//! (1.2e-2) at an exhausted budget of 500. The status differed in 3 cases
//! (2), and `nfev` in 9 of the natural terminations. The agreement within
//! 1e-9 of the class-test and FeatureFinderCentroided_1 fits holds for those
//! fixtures only. The root cause is in
//! `src/math/fitters/levenberg_marquardt.rs`, under investigation in lane B3b.
//! The inputs and the executed results are the fixture
//! `tests/data/gauss_trace_fitter/solver_gap.tsv`; an ignored test prints the
//! current deviations.
//!
//! The support document `docs/TRACE_FITTER_SUPPORT.md` records the API
//! mapping, the native differences and the executed evidence.
//!
//! [`value`]: crate::analysis::feature_finder_picked::trace_fitter::TraceFitter::value

use crate::analysis::feature_finder_picked::helper_structs::{MassTrace, MassTraces};
use crate::analysis::feature_finder_picked::trace_fitter::{
    FEWER_RESIDUALS_THAN_PARAMETERS, ProfileSmoothing, TraceFitter, TraceFitterParams,
    compute_theoretical, initial_shape, optimize, stream_number, unable_to_fit,
};
use crate::math::fitters::levenberg_marquardt::DenseMatrix;
use crate::{Error, Result};

/// Fitter for retention-time profiles with a Gaussian model: source
/// `GaussTraceFitter`.
///
/// A new fitter holds the default [`TraceFitterParams`] (`max_iteration` 500,
/// unweighted) and zero model parameters; the source leaves the model members
/// uninitialised until the first fit. The queries of [`TraceFitter`] read the
/// parameters of the last [`TraceFitter::fit`] or
/// [`Self::set_initial_parameters`].
///
/// `Clone` is the source's copy constructor and assignment operator, with one
/// difference: the source copies `height_`, `x0_` and `sigma_` but not
/// `region_rt_span_`, which the copy leaves uninitialised, so
/// `checkMaximalRTSpan` on a copy reads an indeterminate value. The clone
/// copies it. `FeatureFinderAlgorithmPicked` never copies a fitter.
#[derive(Clone, Debug, PartialEq)]
pub struct GaussTraceFitter {
    parameters: TraceFitterParams,
    sigma: f64,
    x0: f64,
    height: f64,
    region_rt_span: f64,
}

impl Default for GaussTraceFitter {
    fn default() -> Self {
        Self::new()
    }
}

impl GaussTraceFitter {
    /// Number of model parameters, `height`, `x0` and `sigma`: source
    /// `NUM_PARAMS_`.
    pub const NUM_PARAMS: usize = 3;

    /// A fitter with the default parameters: source `GaussTraceFitter()`.
    pub fn new() -> Self {
        Self::with_parameters(TraceFitterParams::default())
    }

    /// A fitter with the given parameters: source `GaussTraceFitter()`
    /// followed by `setParameters`.
    pub fn with_parameters(parameters: TraceFitterParams) -> Self {
        Self {
            parameters,
            sigma: 0.0,
            x0: 0.0,
            height: 0.0,
            region_rt_span: 0.0,
        }
    }

    /// The standard deviation of the fitted Gaussian: source `getSigma`.
    ///
    /// Non-negative after a successful fit, which stores the absolute value of
    /// the optimised parameter. After [`Self::set_initial_parameters`] it is the
    /// start value, which may be NaN.
    pub fn sigma(&self) -> f64 {
        self.sigma
    }

    /// The retention-time span of the region the start values were derived
    /// from: the source's protected `region_rt_span_`, which
    /// [`TraceFitter::check_maximal_rt_span`] compares against.
    ///
    /// The source has no getter; this one exposes the member for the checks
    /// and the evidence.
    pub fn region_rt_span(&self) -> f64 {
        self.region_rt_span
    }

    /// Derive the start values from `traces`: the source's protected
    /// `setInitialParameters_`.
    ///
    /// Takes [`initial_shape`] with [`ProfileSmoothing::SkipShortProfiles`] and
    /// stores its height, apex retention time (as `x0`) and region span. From
    /// the half-maximum heights and retention times it computes
    /// `delta_x = right_rt - left_rt` and
    /// `alpha = ((left_height + right_height) * 0.5) / height`; `sigma` is
    /// `1.0` when `alpha >= 1` (the source's degenerate case, all values the
    /// same) and `(delta_x * 0.5) / sqrt(-2 * ln(alpha))` otherwise. A NaN
    /// `alpha` takes the second branch and gives a NaN `sigma`, as in the
    /// source.
    ///
    /// The source method is protected; it is public here so that the start
    /// values can be compared with the executed C++, and [`TraceFitter::fit`]
    /// calls it first.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `traces` hold no peak, and
    /// propagates the errors of [`MassTraces::intensity_profile`]. The fitter is
    /// unchanged on error.
    pub fn set_initial_parameters(&mut self, traces: &MassTraces) -> Result<()> {
        let shape = initial_shape(traces, ProfileSmoothing::SkipShortProfiles)?;
        let delta_x = shape.right_rt - shape.left_rt;
        let alpha = (shape.left_height + shape.right_height) * 0.5 / shape.height;
        let sigma = if alpha >= 1.0 {
            1.0
        } else {
            delta_x * 0.5 / (-2.0 * ln(alpha)).sqrt()
        };
        self.height = shape.height;
        self.x0 = shape.apex_rt;
        self.region_rt_span = shape.region_rt_span;
        self.sigma = sigma;
        Ok(())
    }

    /// Store an optimised parameter vector `[height, x0, sigma]`: the source's
    /// protected `getOptimizedParameters_`, which keeps `|sigma|`.
    ///
    /// Public for the same reason as [`Self::set_initial_parameters`]; the
    /// source reads its `std::vector` argument without a length check, which
    /// the array type makes unnecessary.
    pub fn set_optimized_parameters(&mut self, x: [f64; 3]) {
        self.height = x[0];
        self.x0 = x[1];
        self.sigma = x[2].abs();
    }
}

impl TraceFitter for GaussTraceFitter {
    fn parameters(&self) -> &TraceFitterParams {
        &self.parameters
    }

    fn set_parameters(&mut self, parameters: TraceFitterParams) {
        self.parameters = parameters;
    }

    /// Fit the Gaussian to `traces`: source `GaussTraceFitter::fit`.
    ///
    /// Sets the start values with [`GaussTraceFitter::set_initial_parameters`],
    /// then runs [`optimize`] on `[height, x0, sigma]` with
    /// [`GaussTraceFunctor`] over all peaks of `traces`, weighted as
    /// [`TraceFitterParams::weighted`] says, and finally stores the result with
    /// [`GaussTraceFitter::set_optimized_parameters`].
    ///
    /// # Errors
    ///
    /// - [`unable_to_fit`] with `"Skipping feature, we always expect N>=p"`
    ///   when `traces` hold fewer than three peaks. With one or two peaks the
    ///   start values are set first, as in the source, and the fitter keeps
    ///   them. Traces without any peak are refused before anything changes: the
    ///   source reads past the end of its empty profile there.
    /// - [`unable_to_fit`] with `"Could not fit the gaussian to the data: Error
    ///   0"` when `max_iteration <= 0`; the fitter keeps the start values.
    /// - The errors of [`MassTraces::intensity_profile`] and the resource
    ///   ceilings of [`optimize`], with the fitter unchanged or holding the
    ///   start values respectively.
    fn fit(&mut self, traces: &MassTraces) -> Result<()> {
        let values = traces.peak_count();
        if values == 0 {
            return Err(unable_to_fit(FEWER_RESIDUALS_THAN_PARAMETERS));
        }
        self.set_initial_parameters(traces)?;
        let mut x = [self.height, self.x0, self.sigma];
        let functor = GaussTraceFunctor::new(traces, self.parameters.weighted);
        optimize(
            &mut x,
            values,
            |point: &[f64], fvec: &mut [f64]| functor.write_residuals(point, fvec),
            |point: &[f64], jac: &mut DenseMatrix| functor.write_jacobian(point, jac),
            &self.parameters,
        )?;
        self.set_optimized_parameters(x);
        Ok(())
    }

    /// `x0 - 2.5 * sigma`: source `getLowerRTBound`.
    fn lower_rt_bound(&self) -> f64 {
        self.x0 - 2.5 * self.sigma
    }

    /// `x0 + 2.5 * sigma`: source `getUpperRTBound`.
    fn upper_rt_bound(&self) -> f64 {
        self.x0 + 2.5 * self.sigma
    }

    /// The fitted `height`: source `getHeight`.
    fn height(&self) -> f64 {
        self.height
    }

    /// The fitted `x0`: source `getCenter`.
    fn center(&self) -> f64 {
        self.x0
    }

    /// `2.35482 * sigma`, the source's truncation of `2 * sqrt(2 * ln 2)`:
    /// source `getFWHM`.
    fn fwhm(&self) -> f64 {
        2.35482 * self.sigma
    }

    /// `height * exp((-0.5 * (rt - x0)^2) / sigma^2)`: source `getValue`.
    fn value(&self, rt: f64) -> f64 {
        self.height * exp(-0.5 * pow2(rt - self.x0) / pow2(self.sigma))
    }

    /// `2.506628 * height * sigma`, the source's truncation of
    /// `sqrt(2 * pi)`: source `getArea`.
    fn area(&self) -> f64 {
        2.506628 * self.height * self.sigma
    }

    /// `(upper - lower) < min_rt_span * 5.0 * sigma`: source
    /// `checkMinimalRTSpan`, `true` when too little of the model remains.
    fn check_minimal_rt_span(&self, rt_bounds: (f64, f64), min_rt_span: f64) -> bool {
        (rt_bounds.1 - rt_bounds.0) < (min_rt_span * 5.0 * self.sigma)
    }

    /// `5.0 * sigma > max_rt_span * region_rt_span`: source
    /// `checkMaximalRTSpan`, `true` when the model is too wide.
    fn check_maximal_rt_span(&self, max_rt_span: f64) -> bool {
        5.0 * self.sigma > max_rt_span * self.region_rt_span
    }

    /// `<name>(x)= <baseline> + <theoretical_int * height> *
    /// exp(-0.5*(x-<rt_shift + x0>)**2/(<sigma>)**2)`: source
    /// `getGnuplotFormula`.
    ///
    /// The source streams `function_name` as one C++ `char`, a single byte; a
    /// `char` outside ASCII is written here as its UTF-8 encoding.
    fn gnuplot_formula(
        &self,
        trace: &MassTrace,
        function_name: char,
        baseline: f64,
        rt_shift: f64,
    ) -> String {
        format!(
            "{function_name}(x)= {} + {} * exp(-0.5*(x-{})**2/({})**2)",
            stream_number(baseline),
            stream_number(trace.theoretical_int * self.height),
            stream_number(rt_shift + self.x0),
            stream_number(self.sigma)
        )
    }

    fn compute_theoretical(&self, trace: &MassTrace, k: usize) -> Result<f64> {
        compute_theoretical(self, trace, k)
    }
}

/// `b * b`: the source's file-local `pow2`.
fn pow2(b: f64) -> f64 {
    b * b
}

/// The exponential the source calls, `exp` of the platform C library. Its last
/// bit differs between C libraries; see the module documentation.
fn exp(x: f64) -> f64 {
    x.exp()
}

/// The natural logarithm the source calls, `log` of the platform C library.
fn ln(x: f64) -> f64 {
    x.ln()
}

/// Residuals and Jacobian of the Gaussian trace model: the source's protected
/// nested class `GaussTraceFitter::GaussTraceFunctor`.
///
/// For the parameter vector `[height, x0, sigma]` and each peak of each trace,
/// in trace order and then peak order, with `d = rt - x0`,
/// `c_fac = -0.5 / sigma^2`, `e = exp(c_fac * d^2)` and `w` the trace's
/// theoretical intensity when weighted and `1.0` otherwise:
///
/// ```text
/// residual = (baseline + theoretical_int * height * e - intensity) * w
/// d/dheight = theoretical_int * e * w
/// d/dx0     = theoretical_int * height * e * d / sigma^2 * w
/// d/dsigma  = 0.125 * theoretical_int * height * e * d^2 / sigma^3 * w
/// ```
///
/// each evaluated left to right as written, with the `f32` intensity promoted
/// to `f64`. The sigma derivative keeps the source's factor `0.125`.
///
/// The source functor is protected; it is public here so that residuals and
/// Jacobians can be compared with the executed C++ at given parameter vectors.
#[derive(Clone, Copy, Debug)]
pub struct GaussTraceFunctor<'a> {
    traces: &'a MassTraces,
    weighted: bool,
    values: usize,
}

impl<'a> GaussTraceFunctor<'a> {
    /// A functor over all peaks of `traces`: source `GaussTraceFunctor(3,
    /// &data)` with `data.weighted = weighted`.
    pub fn new(traces: &'a MassTraces, weighted: bool) -> Self {
        Self {
            traces,
            weighted,
            values: traces.peak_count(),
        }
    }

    /// Number of parameters, 3: source `inputs()`.
    pub fn inputs(&self) -> usize {
        GaussTraceFitter::NUM_PARAMS
    }

    /// Number of residuals, the traces' total peak count: source `values()`.
    pub fn values(&self) -> usize {
        self.values
    }

    /// Fill `fvec` with the residuals at `x = [height, x0, sigma]`: source
    /// `operator()`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] unless `x` holds three values and `fvec`
    /// [`Self::values`] values; `fvec` is then unchanged. The source maps raw
    /// pointers without a check.
    pub fn residuals(&self, x: &[f64], fvec: &mut [f64]) -> Result<()> {
        self.check_point(x)?;
        if fvec.len() != self.values {
            return Err(Error::InvalidValue(format!(
                "a residual vector of {} values does not match {} peaks",
                fvec.len(),
                self.values
            )));
        }
        self.write_residuals(x, fvec);
        Ok(())
    }

    /// Fill `jac` with the Jacobian at `x = [height, x0, sigma]`: source `df`.
    ///
    /// Row `i`, column `j` is the derivative of residual `i` with respect to
    /// parameter `j`. The source fills a column-major Eigen map; the entries are
    /// the same.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] unless `x` holds three values and `jac`
    /// is [`Self::values`] by 3; `jac` is then unchanged.
    pub fn jacobian(&self, x: &[f64], jac: &mut DenseMatrix) -> Result<()> {
        self.check_point(x)?;
        if jac.rows() != self.values || jac.cols() != GaussTraceFitter::NUM_PARAMS {
            return Err(Error::InvalidValue(format!(
                "a {} by {} Jacobian does not match {} peaks and 3 parameters",
                jac.rows(),
                jac.cols(),
                self.values
            )));
        }
        self.write_jacobian(x, jac);
        Ok(())
    }

    fn check_point(&self, x: &[f64]) -> Result<()> {
        if x.len() == GaussTraceFitter::NUM_PARAMS {
            Ok(())
        } else {
            Err(Error::InvalidValue(format!(
                "a Gaussian trace model has 3 parameters, not {}",
                x.len()
            )))
        }
    }

    /// The residual loop; dimensions that do not match leave `fvec` unchanged.
    fn write_residuals(&self, x: &[f64], fvec: &mut [f64]) {
        let &[height, x0, sig] = x else {
            return;
        };
        if fvec.len() != self.values {
            return;
        }
        let c_fac = -0.5 / pow2(sig);
        let baseline = self.traces.baseline;
        let mut slots = fvec.iter_mut();
        for trace in self.traces {
            let weight = if self.weighted {
                trace.theoretical_int
            } else {
                1.0
            };
            for (peak, slot) in trace.peaks.iter().zip(slots.by_ref()) {
                *slot = (baseline
                    + trace.theoretical_int * height * exp(c_fac * pow2(peak.rt - x0))
                    - f64::from(peak.intensity))
                    * weight;
            }
        }
    }

    /// The Jacobian loop; dimensions that do not match leave `jac` unchanged.
    fn write_jacobian(&self, x: &[f64], jac: &mut DenseMatrix) {
        let &[height, x0, sig] = x else {
            return;
        };
        if jac.rows() != self.values || jac.cols() != GaussTraceFitter::NUM_PARAMS {
            return;
        }
        let sig_sq = pow2(sig);
        let inv_sig2 = 1.0 / sig_sq;
        let sig_3 = sig * sig_sq;
        let inv_sig3 = 1.0 / sig_3;
        let c_fac = -0.5 / sig_sq;
        let mut row = 0;
        for trace in self.traces {
            let weight = if self.weighted {
                trace.theoretical_int
            } else {
                1.0
            };
            let theo = trace.theoretical_int;
            for peak in &trace.peaks {
                let rt = peak.rt;
                let e = exp(c_fac * pow2(rt - x0));
                jac.set(row, 0, theo * e * weight);
                jac.set(row, 1, theo * height * e * (rt - x0) * inv_sig2 * weight);
                jac.set(
                    row,
                    2,
                    0.125 * theo * height * e * pow2(rt - x0) * inv_sig3 * weight,
                );
                row += 1;
            }
        }
    }
}
