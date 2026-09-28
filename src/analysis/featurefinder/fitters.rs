// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The five concrete one-dimensional model fitters.
//!
//! Ports `FEATUREFINDER/GaussFitter1D.h`, `FEATUREFINDER/BiGaussFitter1D.h`,
//! `FEATUREFINDER/EmgFitter1D.h`, `FEATUREFINDER/IsotopeFitter1D.h` and
//! `FEATUREFINDER/ExtendedIsotopeFitter1D.h` with their `.cpp` files. See
//! `docs/FEATUREFINDER_MODELS_SUPPORT.md`.
//!
//! # A source behaviour every model-building fitter keeps
//!
//! The Gaussian, bi-Gaussian, EMG and isotope fitters call the new model's
//! `setInterpolationStep(interpolation_step_)` and then `setParameters(tmp)`.
//! `setParameters` replaces the whole parameter tree with `tmp` completed from
//! the defaults, so the fitter's step is discarded and every model is sampled
//! at the model default 0.1. The isotope fitters likewise pass the
//! `isotope_model:` subsection of their parameters to the model and then
//! replace it. Both are reproduced; the fitter's `interpolation_step` still
//! sets the offset scan step of [`max_likeli_fit_offset`].

use crate::Result;
use crate::analysis::feature_finder_picked::glibc_libm;
use crate::concept::constants::PI;
use crate::error::Error;
use crate::kernel::Peak1D;
use crate::math::basic_statistics::BasicStatistics;
use crate::math::fitters::levenberg_marquardt::DenseMatrix;
use crate::param::{Param, ParamValue};

use super::base_model::{
    BaseModel, InterpolationModel, Model1D, complete_parameters, float_value, set_default,
    set_valid,
};
use super::fitter1d::{
    Fitter1D, Fitter1DCore, MAX_ITERATION_DEFAULT, QualityType, lev_marq_optimize,
    max_likeli_fit_offset, model_correlation, write_max_iteration,
};
use super::gauss_models::{BiGaussModel, EmgModel, GaussModel};
use super::isotope_model::{ExtendedIsotopeModel, IsotopeModel};

/// A Gaussian model over the enlarged bounding box, built as the source
/// builds it (see the module documentation on the discarded step).
fn build_gauss(
    step: f64,
    min_bb: f64,
    max_bb: f64,
    statistics: &BasicStatistics,
) -> Result<GaussModel> {
    let mut model = GaussModel::new()?;
    model.set_interpolation_step(step);
    let mut tmp = Param::new();
    tmp.set_value("bounding_box:min", ParamValue::Float(min_bb), "", &[])?;
    tmp.set_value("bounding_box:max", ParamValue::Float(max_bb), "", &[])?;
    tmp.set_value(
        "statistics:mean",
        ParamValue::Float(statistics.mean()),
        "",
        &[],
    )?;
    tmp.set_value(
        "statistics:variance",
        ParamValue::Float(statistics.variance()),
        "",
        &[],
    )?;
    model.set_parameters(&tmp)?;
    Ok(model)
}

/// Source `if (std::isnan(quality)) quality = -1.0;`.
fn nan_to_sentinel(quality: QualityType) -> QualityType {
    if quality.is_nan() { -1.0 } else { quality }
}

/// Gaussian distribution fitter (source `GaussFitter1D`, a
/// `MaxLikeliFitter1D`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GaussFitter1D {
    /// The shared parameters; `statistics:mean` and `statistics:variance`
    /// describe the Gaussian.
    pub core: Fitter1DCore,
}

impl GaussFitter1D {
    /// The handler name (source `setName("GaussFitter1D")`).
    pub const NAME: &'static str = "GaussFitter1D";

    /// The default parameters (source `getDefaults`): those of
    /// [`super::fitter1d::BaseFitter1D::defaults`], with `statistics:variance`
    /// described as "Variance of the model." and `statistics:mean` as "Mean
    /// value of the model.".
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// A default fitter after [`Fitter1D::set_parameters`] with `param`.
    ///
    /// # Errors
    ///
    /// As `set_parameters`.
    pub fn from_param(param: &Param) -> Result<Self> {
        let mut fitter = Self::default();
        fitter.set_parameters(param)?;
        Ok(fitter)
    }
}

impl Fitter1D for GaussFitter1D {
    fn name(&self) -> &'static str {
        Self::NAME
    }

    fn to_param(&self) -> Result<Param> {
        let mut p = Param::new();
        self.core.write_param(&mut p)?;
        set_default(
            &mut p,
            "statistics:variance",
            ParamValue::Float(self.core.statistics.variance()),
            "Variance of the model.",
            true,
        )?;
        set_default(
            &mut p,
            "statistics:mean",
            ParamValue::Float(self.core.statistics.mean()),
            "Mean value of the model.",
            true,
        )?;
        Ok(p)
    }

    fn set_parameters(&mut self, param: &Param) -> Result<Vec<String>> {
        let (complete, warnings) = complete_parameters(Self::NAME, Self::defaults()?, param)?;
        let mut core = self.core.clone();
        core.read_param(&complete)?;
        self.core = core;
        Ok(warnings)
    }

    /// Fit a [`GaussModel`] of the configured mean and variance: its bounding
    /// box is the data's range enlarged by `tolerance_stdev_bounding_box`
    /// standard deviations on each side, and its offset is fitted by
    /// [`max_likeli_fit_offset`] over one such enlargement either way, at the
    /// fitter's `interpolation_step`. A NaN quality becomes -1.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] for empty `range` (the source reads
    /// `set[0]`), and what building the model and the offset scan refuse.
    fn fit1d(&self, range: &[Peak1D]) -> Result<(QualityType, Model1D)> {
        let (mut min_bb, mut max_bb) = Fitter1DCore::data_range(range)?;
        let stdev = self.core.statistics.variance().sqrt() * self.core.tolerance_stdev_box;
        min_bb -= stdev;
        max_bb += stdev;
        let mut model = Model1D::Gauss(build_gauss(
            self.core.interpolation_step,
            min_bb,
            max_bb,
            &self.core.statistics,
        )?);
        let quality = max_likeli_fit_offset(
            &mut model,
            range,
            stdev,
            stdev,
            self.core.interpolation_step,
        )?;
        Ok((nan_to_sentinel(quality), model))
    }
}

/// Bi-Gaussian distribution fitter (source `BiGaussFitter1D`, a
/// `MaxLikeliFitter1D`).
#[derive(Clone, Debug, PartialEq)]
pub struct BiGaussFitter1D {
    /// The shared parameters; only `statistics:mean` of the statistics is
    /// used.
    pub core: Fitter1DCore,
    /// Mean `statistics:mean` and variance `statistics:variance1` (default 1):
    /// the lower half.
    pub statistics1: BasicStatistics,
    /// Mean `statistics:mean` and variance `statistics:variance2` (default 1):
    /// the upper half.
    pub statistics2: BasicStatistics,
}

impl Default for BiGaussFitter1D {
    /// The source defaults.
    fn default() -> Self {
        let core = Fitter1DCore::default();
        let mut statistics = BasicStatistics::new();
        statistics.set_mean(core.statistics.mean());
        statistics.set_variance(1.0);
        Self {
            core,
            statistics1: statistics,
            statistics2: statistics,
        }
    }
}

impl BiGaussFitter1D {
    /// The handler name (source `setName("BiGaussFitter1D")`).
    pub const NAME: &'static str = "BiGaussFitter1D";

    /// The default parameters (source `getDefaults`): those of
    /// [`super::fitter1d::BaseFitter1D::defaults`], then the advanced
    /// `statistics:variance1` 1 and `statistics:variance2` 1.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// A default fitter after [`Fitter1D::set_parameters`] with `param`.
    ///
    /// # Errors
    ///
    /// As `set_parameters`.
    pub fn from_param(param: &Param) -> Result<Self> {
        let mut fitter = Self::default();
        fitter.set_parameters(param)?;
        Ok(fitter)
    }
}

impl Fitter1D for BiGaussFitter1D {
    fn name(&self) -> &'static str {
        Self::NAME
    }

    fn to_param(&self) -> Result<Param> {
        let mut p = Param::new();
        self.core.write_param(&mut p)?;
        set_default(
            &mut p,
            "statistics:variance1",
            ParamValue::Float(self.statistics1.variance()),
            "Variance of the first gaussian, used for the lower half of the model.",
            true,
        )?;
        set_default(
            &mut p,
            "statistics:variance2",
            ParamValue::Float(self.statistics2.variance()),
            "Variance of the second gaussian, used for the upper half of the model.",
            true,
        )?;
        Ok(p)
    }

    fn set_parameters(&mut self, param: &Param) -> Result<Vec<String>> {
        let (complete, warnings) = complete_parameters(Self::NAME, Self::defaults()?, param)?;
        let mut staged = self.clone();
        staged.core.read_param(&complete)?;
        let mean = float_value(&complete, "statistics:mean")?;
        staged.statistics1.set_mean(mean);
        staged
            .statistics1
            .set_variance(float_value(&complete, "statistics:variance1")?);
        staged.statistics2.set_mean(mean);
        staged
            .statistics2
            .set_variance(float_value(&complete, "statistics:variance2")?);
        *self = staged;
        Ok(warnings)
    }

    /// Fit a [`BiGaussModel`]: the data range is enlarged below by
    /// `tolerance_stdev_bounding_box` standard deviations of the lower half
    /// and above by as many of the upper half; the offset scan spans the same
    /// two enlargements. A NaN quality becomes -1.
    ///
    /// # Errors
    ///
    /// As [`GaussFitter1D`]'s `fit1d`.
    fn fit1d(&self, range: &[Peak1D]) -> Result<(QualityType, Model1D)> {
        let (mut min_bb, mut max_bb) = Fitter1DCore::data_range(range)?;
        let stdev1 = self.statistics1.variance().sqrt() * self.core.tolerance_stdev_box;
        let stdev2 = self.statistics2.variance().sqrt() * self.core.tolerance_stdev_box;
        min_bb -= stdev1;
        max_bb += stdev2;
        let mut model = BiGaussModel::new()?;
        model.set_interpolation_step(self.core.interpolation_step);
        let mut tmp = Param::new();
        tmp.set_value("bounding_box:min", ParamValue::Float(min_bb), "", &[])?;
        tmp.set_value("bounding_box:max", ParamValue::Float(max_bb), "", &[])?;
        tmp.set_value(
            "statistics:mean",
            ParamValue::Float(self.statistics1.mean()),
            "",
            &[],
        )?;
        tmp.set_value(
            "statistics:variance1",
            ParamValue::Float(self.statistics1.variance()),
            "",
            &[],
        )?;
        tmp.set_value(
            "statistics:variance2",
            ParamValue::Float(self.statistics2.variance()),
            "",
            &[],
        )?;
        model.set_parameters(&tmp)?;
        let mut model = Model1D::BiGauss(model);
        let quality = max_likeli_fit_offset(
            &mut model,
            range,
            stdev1,
            stdev2,
            self.core.interpolation_step,
        )?;
        Ok((nan_to_sentinel(quality), model))
    }
}

/// Exponentially modified Gaussian fitter using Levenberg-Marquardt (source
/// `EmgFitter1D`, a `LevMarqFitter1D`).
#[derive(Clone, Debug, PartialEq)]
pub struct EmgFitter1D {
    /// The shared parameters; `statistics:mean` and `statistics:variance` are
    /// passed on to the model.
    pub core: Fitter1DCore,
    /// Parameter `max_iteration` (source `max_iteration_`), default 500: the
    /// Levenberg-Marquardt evaluation budget.
    pub max_iteration: i32,
    /// Parameter `init_mom`, the string `true` or `false` (default): start
    /// from method-of-moments estimates.
    pub init_mom: bool,
}

impl Default for EmgFitter1D {
    /// The source defaults.
    fn default() -> Self {
        Self {
            core: Fitter1DCore::default(),
            max_iteration: MAX_ITERATION_DEFAULT,
            init_mom: false,
        }
    }
}

/// The start parameters of an EMG fit and whether the peak is symmetric
/// (source members `height_`, `width_`, `symmetry_`, `retention_`,
/// `symmetric_`).
#[derive(Clone, Copy, Debug, PartialEq)]
struct EmgStart {
    height: f64,
    width: f64,
    symmetry: f64,
    retention: f64,
    symmetric: bool,
}

/// `sqrt(2 pi)` as the source's static member computes it.
fn sqrt2pi() -> f64 {
    (2.0 * PI).sqrt()
}
/// The source's `emg_const`.
const EMG_CONST: f64 = 2.4055;

impl EmgFitter1D {
    /// The handler name (source `setName("EmgFitter1D")`).
    pub const NAME: &'static str = "EmgFitter1D";

    /// The default parameters (source `getDefaults`): those of
    /// [`super::fitter1d::BaseFitter1D::defaults`], `max_iteration` 500, then
    /// `init_mom` `false` (`true` or `false`), with `statistics:variance`
    /// described as "Variance of the model.".
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// A default fitter after [`Fitter1D::set_parameters`] with `param`.
    ///
    /// # Errors
    ///
    /// As `set_parameters`.
    pub fn from_param(param: &Param) -> Result<Self> {
        let mut fitter = Self::default();
        fitter.set_parameters(param)?;
        Ok(fitter)
    }

    /// Build the fitted or fallback model: an [`EmgModel`] of the given
    /// parameters and bounding box, with the fitter's statistics (source
    /// lambda `build_model`). A non-finite or non-positive fitter step is
    /// replaced by 0.2 before the model's `setParameters` discards it.
    #[allow(clippy::too_many_arguments)]
    fn build_model(
        &self,
        height: f64,
        width: f64,
        symmetry: f64,
        retention: f64,
        bb_min: f64,
        bb_max: f64,
    ) -> Result<Model1D> {
        let mut model = EmgModel::new()?;
        let step = self.core.interpolation_step;
        model.set_interpolation_step(if step.is_finite() && step > 0.0 {
            step
        } else {
            0.2
        });
        let mut p = Param::new();
        let entries = [
            ("bounding_box:min", bb_min),
            ("bounding_box:max", bb_max),
            ("statistics:variance", self.core.statistics.variance()),
            ("statistics:mean", self.core.statistics.mean()),
            ("emg:height", height),
            ("emg:width", width),
            ("emg:symmetry", symmetry),
            ("emg:retention", retention),
        ];
        for (key, value) in entries {
            p.set_value(key, ParamValue::Float(value), "", &[])?;
        }
        model.set_parameters(&p)?;
        Ok(Model1D::Emg(model))
    }

    /// Start parameters from the intensity median (source
    /// `setInitialParameters_`), or from method-of-moments estimates when
    /// `init_mom` is set. `set` has at least four points.
    fn initial_parameters(&self, set: &[Peak1D]) -> EmgStart {
        if self.init_mom {
            return Self::initial_parameters_mom(set);
        }
        let mut sum = 0.0f64;
        for peak in set {
            sum += f64::from(peak.intensity);
        }
        // Source `float count`: the running sum is single precision.
        let mut median = 0usize;
        let mut count = 0.0f32;
        for (i, peak) in set.iter().enumerate() {
            count += peak.intensity;
            if f64::from(count) <= sum / 2.0 {
                median = i;
            }
        }
        let first = set[0].mz;
        let last = set[set.len() - 1].mz;
        let median_pos = set[median].mz;
        let max_peak_width = (last - median_pos).abs();
        let mut start = EmgStart {
            height: f64::from(set[median].intensity),
            width: 0.0,
            symmetry: (last - median_pos).abs() / (median_pos - first).abs(),
            retention: median_pos,
            symmetric: false,
        };
        Self::finish_symmetry(&mut start, max_peak_width);
        start
    }

    /// Method-of-moments start parameters (source `setInitialParametersMOM_`).
    fn initial_parameters_mom(set: &[Peak1D]) -> EmgStart {
        let mut weighted_sum = 0.0f64;
        let mut weight_sum = 0.0f64;
        for peak in set {
            weighted_sum += peak.mz * f64::from(peak.intensity);
            weight_sum += f64::from(peak.intensity);
        }
        let weighted_mean = weighted_sum / weight_sum;
        let mut idx = 0usize;
        let mut sum = weight_sum - f64::from(set[0].intensity);
        // Bounded as in the source (issue #6239).
        while sum > weight_sum / 2.0 && idx + 1 < set.len() {
            idx += 1;
            sum -= f64::from(set[idx].intensity);
        }
        let weighted_median = set[idx].mz;
        let mut weighted_sd = 0.0f64;
        for peak in set {
            let d = weighted_mean - peak.mz;
            weighted_sd += d * d * f64::from(peak.intensity);
        }
        weighted_sd /= weight_sum;
        weighted_sd = weighted_sd.sqrt();
        let weighted_skew = (weighted_mean - weighted_median).abs() / weighted_sd;
        let max_peak_width = (set[set.len() - 1].mz - set[idx].mz).abs();
        let cube = (weighted_skew / 2.0).powf(1.0 / 3.0);
        let mut start = EmgStart {
            height: f64::from(set[idx].intensity),
            width: 0.0,
            symmetry: weighted_sd * cube,
            retention: weighted_mean - weighted_sd * cube,
            symmetric: false,
        };
        Self::finish_symmetry(&mut start, max_peak_width);
        start
    }

    /// The shared tail of both start-parameter methods: a non-finite symmetry
    /// marks the peak symmetric with symmetry 10; a symmetry below 1 gets 5
    /// added ("the computations can lead to an overflow error at very low
    /// values of symmetry"); it is capped at the maximal peak width; and the
    /// width starts equal to it ("better for the emg function to proceed from
    /// narrow peaks").
    fn finish_symmetry(start: &mut EmgStart, max_peak_width: f64) {
        if start.symmetry.is_infinite() || start.symmetry.is_nan() {
            start.symmetric = true;
            start.symmetry = 10.0;
        }
        if start.symmetry < 1.0 {
            start.symmetry += 5.0;
        }
        // std::min(symmetry_, max_peak_width): the first unless the second is
        // strictly smaller.
        if max_peak_width < start.symmetry {
            start.symmetry = max_peak_width;
        }
        start.width = start.symmetry;
    }
}

/// The residuals of the simplified EMG (source `EgmFitterFunctor::operator()`).
fn emg_residuals(set: &[Peak1D], x: &[f64], fvec: &mut [f64]) {
    let (h, w, s, z) = (x[0], x[1], x[2], x[3]);
    let c = -EMG_CONST / 2.0f64.sqrt();
    let prefix = (h * w / s) * sqrt2pi();
    let part1 = w * w / (2.0 * s * s);
    let part2 = w / s;
    for (i, peak) in set.iter().enumerate() {
        let t = peak.mz;
        let diff = t - z;
        let yi = prefix * glibc_libm::exp(part1 - (diff / s))
            / (1.0 + glibc_libm::exp(c * ((diff / w) - part2)));
        fvec[i] = yi - f64::from(peak.intensity);
    }
}

/// The Jacobian of the simplified EMG (source `EgmFitterFunctor::df`), term
/// for term in the source's operation order.
fn emg_jacobian(set: &[Peak1D], x: &[f64], jac: &mut DenseMatrix) {
    let sqrt2pi = sqrt2pi();
    let sqrt_2 = 2.0f64.sqrt();
    let emg_const = EMG_CONST;
    let h = x[0];
    let w = x[1];
    let w2 = w * w;
    let s = x[2];
    let s2 = s * s;
    let s3 = s2 * s;
    let z = x[3];
    for (i, peak) in set.iter().enumerate() {
        let t = peak.mz;
        let diff = t - z;
        let exp1 = glibc_libm::exp((w2 / (2.0 * s2)) - (diff / s));
        let exp3 = glibc_libm::exp((-emg_const / sqrt_2) * ((diff / w) - w / s));
        let exp2 = 1.0 + exp3;
        let derivative_height = w / s * sqrt2pi * exp1 / exp2;
        let derivative_width = h / s * sqrt2pi * exp1 / exp2
            + (h * w2) / s3 * sqrt2pi * exp1 / exp2
            + (emg_const * h * w) / s * sqrt2pi * exp1 * (-diff / w2 - 1.0 / s) * exp3
                / ((exp2 * exp2) * sqrt_2);
        let derivative_symmetry = -h * w / s2 * sqrt2pi * exp1 / exp2
            + h * w / s * sqrt2pi * (-(w * w) / s3 + diff / s2) * exp1 / exp2
            + (emg_const * h * w2) / s3 * sqrt2pi * exp1 * exp3 / ((exp2 * exp2) * sqrt_2);
        let derivative_retention = h * w / s2 * sqrt2pi * exp1 / exp2
            - (emg_const * h) / s * sqrt2pi * exp1 * exp3 / ((exp2 * exp2) * sqrt_2);
        jac.set(i, 0, derivative_height);
        jac.set(i, 1, derivative_width);
        jac.set(i, 2, derivative_symmetry);
        jac.set(i, 3, derivative_retention);
    }
}

impl Fitter1D for EmgFitter1D {
    fn name(&self) -> &'static str {
        Self::NAME
    }

    fn to_param(&self) -> Result<Param> {
        let mut p = Param::new();
        self.core.write_param(&mut p)?;
        write_max_iteration(&mut p, self.max_iteration)?;
        set_default(
            &mut p,
            "init_mom",
            ParamValue::String(if self.init_mom { "true" } else { "false" }.to_owned()),
            "Initialize parameters using method of moments estimators.",
            true,
        )?;
        set_valid(&mut p, "init_mom", &["true", "false"])?;
        set_default(
            &mut p,
            "statistics:variance",
            ParamValue::Float(self.core.statistics.variance()),
            "Variance of the model.",
            true,
        )?;
        Ok(p)
    }

    fn set_parameters(&mut self, param: &Param) -> Result<Vec<String>> {
        let (complete, warnings) = complete_parameters(Self::NAME, Self::defaults()?, param)?;
        let mut staged = self.clone();
        staged.core.read_param(&complete)?;
        staged.max_iteration = complete.value("max_iteration")?.to_i32()?;
        staged.init_mom = complete.value("init_mom")?.to_bool()?;
        *self = staged;
        Ok(warnings)
    }

    /// Fit an [`EmgModel`] by Levenberg-Marquardt and return its Pearson
    /// correlation with the data (source `fit1d`, hardened for issue #6239).
    ///
    /// The fit never fails on degenerate data: fewer than four points, a
    /// non-finite position or intensity, a zero span, a non-positive or
    /// non-finite intensity sum or fitter step, a refused optimisation, a
    /// non-finite optimum or a non-finite correlation each yield quality -1
    /// and a finite fallback model centred on the most intense point (width
    /// and symmetry a sixth of the data span, at least 0.001). Otherwise the
    /// start parameters come from `init_mom`'s method, a symmetric start skips
    /// the optimisation, and the model spans the data range enlarged by
    /// `tolerance_stdev_bounding_box` standard deviations.
    ///
    /// # Errors
    ///
    /// Only what building a model refuses (a bounding box beyond
    /// [`super::base_model::MAX_SAMPLES`] samples) and the optimiser's resource
    /// ceilings; the source has neither.
    // The source's guards are NaN-sensitive (`!(span > 0.0)` is true for NaN);
    // they are transcribed literally rather than through `partial_cmp`.
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    fn fit1d(&self, set: &[Peak1D]) -> Result<(QualityType, Model1D)> {
        let n = set.len();
        if n < 4 {
            let (mut lo, mut hi, mut apex, mut height) = (0.0f64, 1.0f64, 0.0f64, 1.0f64);
            if let Some(first) = set.first() {
                lo = first.mz;
                hi = first.mz;
                apex = first.mz;
                height = f64::from(first.intensity);
                for peak in &set[1..] {
                    let pos = peak.mz;
                    if pos.is_finite() {
                        lo = lo.min(pos);
                        hi = hi.max(pos);
                    }
                    if f64::from(peak.intensity) > height {
                        height = f64::from(peak.intensity);
                        apex = pos;
                    }
                }
            }
            if !lo.is_finite() || !hi.is_finite() || hi <= lo {
                lo = 0.0;
                hi = 1.0;
            }
            if !apex.is_finite() || apex < lo || apex > hi {
                apex = lo;
            }
            let w = ((hi - lo) / 6.0).max(1e-3);
            let model = self.build_model(
                if height.is_finite() { height } else { 1.0 },
                w,
                w,
                apex,
                lo,
                hi,
            )?;
            return Ok((-1.0, model));
        }

        let mut min_pos = set[0].mz;
        let mut max_pos = set[0].mz;
        let mut apex_pos = set[0].mz;
        let mut max_int = f64::from(set[0].intensity);
        let mut weight_sum = 0.0f64;
        let mut finite_input = set[0].mz.is_finite() && set[0].intensity.is_finite();
        for peak in set {
            let pos = peak.mz;
            let intensity = f64::from(peak.intensity);
            if !pos.is_finite() || !intensity.is_finite() {
                finite_input = false;
            }
            if pos < min_pos {
                min_pos = pos;
            }
            if pos > max_pos {
                max_pos = pos;
            }
            if intensity > max_int {
                max_int = intensity;
                apex_pos = pos;
            }
            weight_sum += intensity;
        }
        if !min_pos.is_finite() || !max_pos.is_finite() || !(max_pos > min_pos) {
            let center = if apex_pos.is_finite() { apex_pos } else { 0.0 };
            min_pos = center - 0.5;
            max_pos = center + 0.5;
        }
        let span = max_pos - min_pos;
        let fb_w = (span / 6.0).max(1e-3);
        if !apex_pos.is_finite() || apex_pos < min_pos || apex_pos > max_pos {
            apex_pos = min_pos;
        }
        if !max_int.is_finite() || !(max_int > 0.0) {
            max_int = 1.0;
        }
        let step = self.core.interpolation_step;
        if !finite_input
            || !(span > 0.0)
            || !weight_sum.is_finite()
            || !(weight_sum > 0.0)
            || !step.is_finite()
            || !(step > 0.0)
        {
            let model = self.build_model(max_int, fb_w, fb_w, apex_pos, min_pos, max_pos)?;
            return Ok((-1.0, model));
        }

        let stdev = self.core.statistics.variance().sqrt() * self.core.tolerance_stdev_box;
        let min_bb = min_pos - stdev;
        let max_bb = max_pos + stdev;

        let start = self.initial_parameters(set);
        let mut x = [start.height, start.width, start.symmetry, start.retention];
        let mut fit_ok = true;
        if !start.symmetric {
            match lev_marq_optimize(
                &mut x,
                n,
                |p: &[f64], f: &mut [f64]| emg_residuals(set, p, f),
                |p: &[f64], j: &mut DenseMatrix| emg_jacobian(set, p, j),
                self.max_iteration,
            ) {
                Ok(()) => {}
                // The source catches exactly Exception::UnableToFit.
                Err(Error::InvalidValue(message))
                    if message.starts_with(
                        crate::analysis::feature_finder_picked::trace_fitter::UNABLE_TO_FIT_FINAL_SET,
                    ) =>
                {
                    fit_ok = false;
                }
                Err(other) => return Err(other),
            }
        }
        let [height, width, symmetry, retention] = x;
        if !fit_ok
            || !height.is_finite()
            || !width.is_finite()
            || !symmetry.is_finite()
            || !retention.is_finite()
        {
            let model = self.build_model(max_int, fb_w, fb_w, apex_pos, min_pos, max_pos)?;
            return Ok((-1.0, model));
        }
        let model = self.build_model(height, width, symmetry, retention, min_bb, max_bb)?;
        let correlation = model_correlation(&model, set)?;
        if !correlation.is_finite() {
            let model = self.build_model(max_int, fb_w, fb_w, apex_pos, min_pos, max_pos)?;
            return Ok((-1.0, model));
        }
        Ok((correlation, model))
    }
}

/// The `isotope_model:` subsection of a fitter's completed parameters, with
/// every root key starting with `stdev` removed (source
/// `param_.copy("isotope_model:", true)` then `removeAll("stdev")`).
fn isotope_model_subsection(sub: &Param) -> Result<Param> {
    let mut iso = sub.checked_clone()?;
    iso.remove_all("stdev")?;
    Ok(iso)
}

/// Isotope distribution fitter (source `IsotopeFitter1D`, a
/// `MaxLikeliFitter1D`).
///
/// With charge 0 it fits a [`GaussModel`] as [`GaussFitter1D`] does; otherwise
/// an [`IsotopeModel`] at `statistics:mean`, sampled from its averagine
/// formula, whose offset is then fitted.
#[derive(Clone, Debug, PartialEq)]
pub struct IsotopeFitter1D {
    /// The shared parameters, with `interpolation_step` defaulting to 0.1.
    pub core: Fitter1DCore,
    /// Parameter `charge` (Int, default 1), held as the source's `double`
    /// member.
    pub charge: f64,
    /// Parameter `isotope:stdev`, default 1: the Gaussian widening.
    pub isotope_stdev: f64,
    /// Parameter `isotope:maximum`, default 100.
    pub max_isotope: i32,
    /// The `isotope_model:` subsection of the parameters, if any. The source
    /// keeps such (unknown, warned-about) keys in its tree and passes them to
    /// the model, whose second `setParameters` then discards them; they can
    /// only make the fit fail validation.
    pub isotope_model: Param,
}

impl Default for IsotopeFitter1D {
    /// The source defaults.
    fn default() -> Self {
        Self {
            core: Fitter1DCore {
                interpolation_step: 0.1,
                ..Fitter1DCore::default()
            },
            charge: 1.0,
            isotope_stdev: 1.0,
            max_isotope: 100,
            isotope_model: Param::new(),
        }
    }
}

/// Read `charge` as the source's `CoordinateType charge_ =
/// param_.getValue("charge")`.
fn charge_value(param: &Param) -> Result<f64> {
    Ok(param.value("charge")?.to_i64()? as f64)
}

impl IsotopeFitter1D {
    /// The handler name (source `setName("IsotopeFitter1D")`).
    pub const NAME: &'static str = "IsotopeFitter1D";

    /// The default parameters (source `getDefaults`): those of
    /// [`super::fitter1d::BaseFitter1D::defaults`] with `interpolation_step`
    /// 0.1 and `statistics:variance` described as "Variance of the model.",
    /// then the advanced `charge` 1, `isotope:stdev` 1 and `isotope:maximum`
    /// 100.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().params(false)
    }

    /// A default fitter after [`Fitter1D::set_parameters`] with `param`.
    ///
    /// # Errors
    ///
    /// As `set_parameters`.
    pub fn from_param(param: &Param) -> Result<Self> {
        let mut fitter = Self::default();
        fitter.set_parameters(param)?;
        Ok(fitter)
    }

    fn params(&self, with_subsection: bool) -> Result<Param> {
        let mut p = Param::new();
        self.core.write_param(&mut p)?;
        set_default(
            &mut p,
            "statistics:variance",
            ParamValue::Float(self.core.statistics.variance()),
            "Variance of the model.",
            true,
        )?;
        set_default(
            &mut p,
            "charge",
            ParamValue::Integer(self.charge as i64),
            "Charge state of the model.",
            true,
        )?;
        set_default(
            &mut p,
            "isotope:stdev",
            ParamValue::Float(self.isotope_stdev),
            "Standard deviation of gaussian applied to the averagine isotopic pattern to simulate the inaccuracy of the mass spectrometer.",
            true,
        )?;
        set_default(
            &mut p,
            "isotope:maximum",
            ParamValue::Integer(i64::from(self.max_isotope)),
            "Maximum isotopic rank to be considered.",
            true,
        )?;
        set_default(
            &mut p,
            "interpolation_step",
            ParamValue::Float(self.core.interpolation_step),
            "Sampling rate for the interpolation of the model function.",
            true,
        )?;
        if with_subsection && !self.isotope_model.is_empty() {
            p.insert("isotope_model:", &self.isotope_model)?;
        }
        Ok(p)
    }
}

impl Fitter1D for IsotopeFitter1D {
    fn name(&self) -> &'static str {
        Self::NAME
    }

    /// The parameters, including an `isotope_model:` subsection when one was
    /// set.
    fn to_param(&self) -> Result<Param> {
        self.params(true)
    }

    fn set_parameters(&mut self, param: &Param) -> Result<Vec<String>> {
        let (complete, warnings) = complete_parameters(Self::NAME, Self::defaults()?, param)?;
        let mut staged = self.clone();
        staged.core.read_param(&complete)?;
        staged.charge = charge_value(&complete)?;
        staged.isotope_stdev = float_value(&complete, "isotope:stdev")?;
        staged.max_isotope = complete.value("isotope:maximum")?.to_i32()?;
        staged.isotope_model = complete.copy("isotope_model:", true)?;
        *self = staged;
        Ok(warnings)
    }

    /// Fit a Gaussian (charge 0) or an [`IsotopeModel`] and its offset.
    ///
    /// For the isotope model the source passes the `isotope_model:`
    /// subsection (without `stdev*` keys) to the model, sets the step, and
    /// then replaces all of that with `statistics:mean`, `charge` (truncated
    /// to Int), `isotope:mode:GaussianSD` = `isotope:stdev` and
    /// `isotope:maximum`, so the model samples at its default step 0.1 from
    /// [`IsotopeModel::formula`]. The offset scan spans the enlarged-box
    /// standard deviation either way at the fitter's step. A NaN quality
    /// becomes -1.
    ///
    /// # Errors
    ///
    /// As [`GaussFitter1D`]'s `fit1d`, and what the isotope model refuses (see
    /// [`IsotopeModel::set_samples_from_formula`]).
    fn fit1d(&self, range: &[Peak1D]) -> Result<(QualityType, Model1D)> {
        let (mut min_bb, mut max_bb) = Fitter1DCore::data_range(range)?;
        let stdev = self.core.statistics.variance().sqrt() * self.core.tolerance_stdev_box;
        min_bb -= stdev;
        max_bb += stdev;
        let mut model = if self.charge == 0.0 {
            Model1D::Gauss(build_gauss(
                self.core.interpolation_step,
                min_bb,
                max_bb,
                &self.core.statistics,
            )?)
        } else {
            let mut model = IsotopeModel::new()?;
            model.set_parameters(&isotope_model_subsection(&self.isotope_model)?)?;
            model.set_interpolation_step(self.core.interpolation_step);
            let mut tmp = Param::new();
            tmp.set_value(
                "statistics:mean",
                ParamValue::Float(self.core.statistics.mean()),
                "",
                &[],
            )?;
            tmp.set_value(
                "charge",
                ParamValue::Integer(self.charge as i32 as i64),
                "",
                &[],
            )?;
            tmp.set_value(
                "isotope:mode:GaussianSD",
                ParamValue::Float(self.isotope_stdev),
                "",
                &[],
            )?;
            tmp.set_value(
                "isotope:maximum",
                ParamValue::Integer(i64::from(self.max_isotope)),
                "",
                &[],
            )?;
            model.set_parameters(&tmp)?;
            let formula = model.formula()?;
            model.set_samples_from_formula(&formula)?;
            Model1D::Isotope(model)
        };
        let quality = max_likeli_fit_offset(
            &mut model,
            range,
            stdev,
            stdev,
            self.core.interpolation_step,
        )?;
        Ok((nan_to_sentinel(quality), model))
    }
}

/// Extended isotope distribution fitter (source `ExtendedIsotopeFitter1D`, a
/// `MaxLikeliFitter1D`).
///
/// With charge 0 it builds a [`GaussModel`] as [`GaussFitter1D`] does, but
/// without fitting its offset; otherwise an [`ExtendedIsotopeModel`] at
/// `isotope:monoisotopic_mz`. Either way the quality is the model's plain
/// Pearson correlation with the data.
#[derive(Clone, Debug, PartialEq)]
pub struct ExtendedIsotopeFitter1D {
    /// The shared parameters.
    pub core: Fitter1DCore,
    /// Parameter `charge` (Int, default 1), held as the source's `double`.
    pub charge: f64,
    /// Parameter `isotope:stdev`, default 0.1.
    pub isotope_stdev: f64,
    /// Parameter `isotope:monoisotopic_mz`, default 1.
    pub monoisotopic_mz: f64,
    /// Parameter `isotope:maximum`, default 100.
    pub max_isotope: i32,
    /// The `isotope_model:` subsection, as for [`IsotopeFitter1D`].
    pub isotope_model: Param,
}

impl Default for ExtendedIsotopeFitter1D {
    /// The source defaults.
    fn default() -> Self {
        Self {
            core: Fitter1DCore::default(),
            charge: 1.0,
            isotope_stdev: 0.1,
            monoisotopic_mz: 1.0,
            max_isotope: 100,
            isotope_model: Param::new(),
        }
    }
}

impl ExtendedIsotopeFitter1D {
    /// The handler name (source `setName("ExtendedIsotopeFitter1D")`).
    pub const NAME: &'static str = "ExtendedIsotopeFitter1D";

    /// The default parameters (source `getDefaults`): those of
    /// [`super::fitter1d::BaseFitter1D::defaults`] with `statistics:variance`
    /// described as "Variance of the model.", then the advanced `charge` 1,
    /// `isotope:stdev` 0.1, `isotope:monoisotopic_mz` 1 and `isotope:maximum`
    /// 100 (and `interpolation_step` 0.2, redeclared in place).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().params(false)
    }

    /// A default fitter after [`Fitter1D::set_parameters`] with `param`.
    ///
    /// # Errors
    ///
    /// As `set_parameters`.
    pub fn from_param(param: &Param) -> Result<Self> {
        let mut fitter = Self::default();
        fitter.set_parameters(param)?;
        Ok(fitter)
    }

    fn params(&self, with_subsection: bool) -> Result<Param> {
        let mut p = Param::new();
        self.core.write_param(&mut p)?;
        set_default(
            &mut p,
            "statistics:variance",
            ParamValue::Float(self.core.statistics.variance()),
            "Variance of the model.",
            true,
        )?;
        set_default(
            &mut p,
            "charge",
            ParamValue::Integer(self.charge as i64),
            "Charge state of the model.",
            true,
        )?;
        set_default(
            &mut p,
            "isotope:stdev",
            ParamValue::Float(self.isotope_stdev),
            "Standard deviation of gaussian applied to the averagine isotopic pattern to simulate the inaccuracy of the mass spectrometer.",
            true,
        )?;
        set_default(
            &mut p,
            "isotope:monoisotopic_mz",
            ParamValue::Float(self.monoisotopic_mz),
            "Monoisotopic m/z of the model.",
            true,
        )?;
        set_default(
            &mut p,
            "isotope:maximum",
            ParamValue::Integer(i64::from(self.max_isotope)),
            "Maximum isotopic rank to be considered.",
            true,
        )?;
        set_default(
            &mut p,
            "interpolation_step",
            ParamValue::Float(self.core.interpolation_step),
            "Sampling rate for the interpolation of the model function.",
            true,
        )?;
        if with_subsection && !self.isotope_model.is_empty() {
            p.insert("isotope_model:", &self.isotope_model)?;
        }
        Ok(p)
    }
}

impl Fitter1D for ExtendedIsotopeFitter1D {
    fn name(&self) -> &'static str {
        Self::NAME
    }

    /// The parameters, including an `isotope_model:` subsection when one was
    /// set.
    fn to_param(&self) -> Result<Param> {
        self.params(true)
    }

    fn set_parameters(&mut self, param: &Param) -> Result<Vec<String>> {
        let (complete, warnings) = complete_parameters(Self::NAME, Self::defaults()?, param)?;
        let mut staged = self.clone();
        staged.core.read_param(&complete)?;
        staged.charge = charge_value(&complete)?;
        staged.isotope_stdev = float_value(&complete, "isotope:stdev")?;
        staged.monoisotopic_mz = float_value(&complete, "isotope:monoisotopic_mz")?;
        staged.max_isotope = complete.value("isotope:maximum")?.to_i32()?;
        staged.isotope_model = complete.copy("isotope_model:", true)?;
        *self = staged;
        Ok(warnings)
    }

    /// Build the model and return its Pearson correlation with the data; a
    /// NaN correlation becomes -1. The offset is not fitted.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] for empty `range` (for charge 0 the
    /// source reads `set[0]`, otherwise its correlation throws
    /// `Exception::InvalidRange`), and what building the model refuses.
    fn fit1d(&self, range: &[Peak1D]) -> Result<(QualityType, Model1D)> {
        let model = if self.charge == 0.0 {
            let (mut min_bb, mut max_bb) = Fitter1DCore::data_range(range)?;
            let stdev = self.core.statistics.variance().sqrt() * self.core.tolerance_stdev_box;
            min_bb -= stdev;
            max_bb += stdev;
            Model1D::Gauss(build_gauss(
                self.core.interpolation_step,
                min_bb,
                max_bb,
                &self.core.statistics,
            )?)
        } else {
            let mut model = ExtendedIsotopeModel::new()?;
            model.set_parameters(&isotope_model_subsection(&self.isotope_model)?)?;
            model.set_interpolation_step(self.core.interpolation_step);
            let mut tmp = Param::new();
            tmp.set_value(
                "isotope:monoisotopic_mz",
                ParamValue::Float(self.monoisotopic_mz),
                "",
                &[],
            )?;
            tmp.set_value(
                "charge",
                ParamValue::Integer(self.charge as i32 as i64),
                "",
                &[],
            )?;
            tmp.set_value(
                "isotope:stdev",
                ParamValue::Float(self.isotope_stdev),
                "",
                &[],
            )?;
            tmp.set_value(
                "isotope:maximum",
                ParamValue::Integer(i64::from(self.max_isotope)),
                "",
                &[],
            )?;
            model.set_parameters(&tmp)?;
            Model1D::ExtendedIsotope(model)
        };
        let correlation = model_correlation(&model, range)?;
        Ok((nan_to_sentinel(correlation), model))
    }
}
