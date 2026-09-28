// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Elution-profile models sampled over a bounding box: the Gaussian, the
//! bi-Gaussian and the exponentially modified Gaussian.
//!
//! Ports `FEATUREFINDER/GaussModel.h`, `FEATUREFINDER/BiGaussModel.h` and
//! `FEATUREFINDER/EmgModel.h` with their `.cpp` files. See
//! `docs/FEATUREFINDER_MODELS_SUPPORT.md`.
//!
//! All three sample their function at `min + i * interpolation_step` for
//! `i = 0, 1, ...` up to and including the first position at or beyond `max`,
//! so the last sample lies at or past the bounding box. A bounding box with
//! `min == max` leaves the interpolation empty and keeps its previous scale and
//! offset, as the source's early return does.

use crate::Result;
use crate::analysis::feature_finder_picked::glibc_libm;
use crate::concept::constants::PI;
use crate::error::Error;
use crate::math::basic_statistics::BasicStatistics;
use crate::param::{Param, ParamValue};

use super::base_model::{
    InterpolationCore, InterpolationModel, MAX_SAMPLES, complete_parameters, float_value,
    interpolation_base_model, set_default,
};

/// Sample `function` over the bounding box `[min, max]` at the core's
/// interpolation step, in the source's loop shape.
///
/// Returns `None` for `min == max` (the source returns before sampling, after
/// clearing the data).
fn sample_box(
    core: &InterpolationCore,
    min: f64,
    max: f64,
    mut function: impl FnMut(f64) -> f64,
) -> Result<Option<Vec<f64>>> {
    if max == min {
        return Ok(None);
    }
    let step = core.interpolation_step;
    let mut data = Vec::new();
    if min < max {
        // The source loops forever on a non-positive or NaN step and samples
        // without bound on an infinite box; both are refused up front.
        let estimate = (max - min) / step;
        if !(step > 0.0 && estimate.is_finite() && estimate < MAX_SAMPLES as f64) {
            return Err(Error::InvalidValue(format!(
                "the bounding box [{min}, {max}] at interpolation step {step} needs more than \
                 {MAX_SAMPLES} samples or never ends"
            )));
        }
        data.try_reserve_exact(estimate as usize + 2)
            .map_err(|_| Error::InvalidValue("model sample allocation failed".into()))?;
    }
    let mut pos = min;
    let mut i: u32 = 0;
    while pos < max {
        if data.len() >= MAX_SAMPLES {
            return Err(Error::InvalidValue(format!(
                "the model needs more than {MAX_SAMPLES} samples"
            )));
        }
        pos = min + f64::from(i) * step;
        data.push(function(pos));
        i += 1;
    }
    Ok(Some(data))
}

/// Scale `data` so that the rectangle-rule integral over it equals `scaling`:
/// `factor = scaling / step / sum`, the source's association for the Gaussian
/// models. The sum is taken left to right, as `std::accumulate`.
fn normalize(data: &mut [f64], scaling: f64, step: f64) {
    let sum: f64 = data.iter().fold(0.0, |acc, &v| acc + v);
    let factor = scaling / step / sum;
    for value in data {
        *value *= factor;
    }
}

/// Install freshly sampled `data` with scale `step` and offset `min`, or clear
/// the data (keeping scale and offset) for an empty bounding box.
fn install(core: &mut InterpolationCore, data: Option<Vec<f64>>, min: f64) {
    match data {
        None => core.interpolation.data_mut().clear(),
        Some(data) => {
            core.interpolation.set_data(data);
            core.interpolation.set_scale(core.interpolation_step);
            core.interpolation.set_offset(min);
        }
    }
}

/// Normal distribution approximated using linear interpolation (source
/// `GaussModel`).
///
/// The samples are `normalDensity_sqrt2pi` of
/// [`BasicStatistics`] over the bounding box, scaled so that their rectangle
/// integral is the scaling factor.
#[derive(Clone, Debug, PartialEq)]
pub struct GaussModel {
    core: InterpolationCore,
    /// Parameter `bounding_box:min` (source `min_`), default 0.
    min: f64,
    /// Parameter `bounding_box:max` (source `max_`), default 1.
    max: f64,
    /// Parameters `statistics:mean` (default 0) and `statistics:variance`
    /// (default 1) (source `statistics_`).
    statistics: BasicStatistics,
}

impl GaussModel {
    /// The handler name (source `setName("GaussModel")`).
    pub const NAME: &'static str = "GaussModel";

    fn default_state() -> Self {
        let mut statistics = BasicStatistics::new();
        statistics.set_mean(0.0);
        statistics.set_variance(1.0);
        Self {
            core: InterpolationCore::default(),
            min: 0.0,
            max: 1.0,
            statistics,
        }
    }

    /// A model with the default parameters, sampled (source default
    /// constructor, which ends in `defaultsToParam_` and hence `setSamples`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure; the defaults sample 11 points.
    pub fn new() -> Result<Self> {
        let mut model = Self::default_state();
        model.set_samples()?;
        Ok(model)
    }

    /// The default parameters (source `getDefaults`): `cutoff`,
    /// `interpolation_step`, `intensity_scaling`, then the advanced
    /// `bounding_box:min` 0, `bounding_box:max` 1, `statistics:mean` 0 and
    /// `statistics:variance` 1.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default_state().to_param_impl()
    }

    /// A default model after [`crate::analysis::featurefinder::base_model::BaseModel::set_parameters`]
    /// with `param`.
    ///
    /// # Errors
    ///
    /// As `set_parameters`.
    pub fn from_param(param: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(param)?.0)
    }

    /// As [`GaussModel::from_param`], also returning the unknown-parameter
    /// warnings the source logs.
    ///
    /// # Errors
    ///
    /// As `set_parameters`.
    pub fn from_param_with_warnings(param: &Param) -> Result<(Self, Vec<String>)> {
        let mut model = Self::new()?;
        let warnings = model.set_parameters_impl(param)?;
        Ok((model, warnings))
    }

    /// The bounding box `(min, max)` the model is sampled over.
    pub fn bounding_box(&self) -> (f64, f64) {
        (self.min, self.max)
    }

    /// The mean and variance of the Gaussian.
    pub fn statistics(&self) -> &BasicStatistics {
        &self.statistics
    }

    fn to_param_impl(&self) -> Result<Param> {
        let mut p = Param::new();
        self.core.write_param(&mut p)?;
        set_default(
            &mut p,
            "bounding_box:min",
            ParamValue::Float(self.min),
            "Lower end of bounding box enclosing the data used to fit the model.",
            true,
        )?;
        set_default(
            &mut p,
            "bounding_box:max",
            ParamValue::Float(self.max),
            "Upper end of bounding box enclosing the data used to fit the model.",
            true,
        )?;
        set_default(
            &mut p,
            "statistics:mean",
            ParamValue::Float(self.statistics.mean()),
            "Centroid position of the model (Gaussian).",
            true,
        )?;
        set_default(
            &mut p,
            "statistics:variance",
            ParamValue::Float(self.statistics.variance()),
            "The variance of the Gaussian.",
            true,
        )?;
        Ok(p)
    }

    fn set_parameters_impl(&mut self, param: &Param) -> Result<Vec<String>> {
        let (complete, warnings) = complete_parameters(Self::NAME, Self::defaults()?, param)?;
        let mut staged = self.clone();
        staged.core.read_param(&complete)?;
        staged.min = float_value(&complete, "bounding_box:min")?;
        staged.max = float_value(&complete, "bounding_box:max")?;
        staged
            .statistics
            .set_mean(float_value(&complete, "statistics:mean")?);
        staged
            .statistics
            .set_variance(float_value(&complete, "statistics:variance")?);
        staged.set_samples()?;
        *self = staged;
        Ok(warnings)
    }
}

interpolation_base_model!(GaussModel);

impl InterpolationModel for GaussModel {
    fn core(&self) -> &InterpolationCore {
        &self.core
    }
    fn core_mut(&mut self) -> &mut InterpolationCore {
        &mut self.core
    }

    /// The position of the maximum: the mean.
    fn center(&self) -> Result<f64> {
        Ok(self.statistics.mean())
    }

    /// Shift the model without recomputing it: the bounding box and the mean
    /// move by `offset - interpolation().offset()`, and the parameters report
    /// the moved values.
    fn set_offset(&mut self, offset: f64) {
        let diff = offset - self.core.interpolation.offset();
        self.min += diff;
        self.max += diff;
        let mean = self.statistics.mean() + diff;
        self.statistics.set_mean(mean);
        self.core.interpolation.set_offset(offset);
    }

    /// Sample `normalDensity_sqrt2pi` over the bounding box and scale the
    /// samples to integrate to the scaling factor.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a non-empty bounding box needs more
    /// than [`MAX_SAMPLES`] samples, or the step is not positive (the source
    /// loops forever). A zero sum is not refused: as in the source, every
    /// sample becomes NaN or infinite.
    fn set_samples(&mut self) -> Result<()> {
        let statistics = self.statistics;
        let data = sample_box(&self.core, self.min, self.max, |pos| {
            statistics.normal_density_sqrt2pi(pos)
        })?
        .map(|mut data| {
            normalize(&mut data, self.core.scaling, self.core.interpolation_step);
            data
        });
        install(&mut self.core, data, self.min);
        Ok(())
    }
}

/// Asymmetric distribution of two normal distributions with different
/// variances combined at their common mean, approximated using linear
/// interpolation (source `BiGaussModel`).
#[derive(Clone, Debug, PartialEq)]
pub struct BiGaussModel {
    core: InterpolationCore,
    /// Parameter `bounding_box:min`, default 0.
    min: f64,
    /// Parameter `bounding_box:max`, default 1.
    max: f64,
    /// Mean `statistics:mean` (default 0) and variance `statistics:variance1`
    /// (default 1): the lower half.
    statistics1: BasicStatistics,
    /// Mean `statistics:mean` and variance `statistics:variance2` (default 1):
    /// the upper half.
    statistics2: BasicStatistics,
}

impl BiGaussModel {
    /// The handler name (source `setName("BiGaussModel")`).
    pub const NAME: &'static str = "BiGaussModel";

    fn default_state() -> Self {
        let mut statistics = BasicStatistics::new();
        statistics.set_mean(0.0);
        statistics.set_variance(1.0);
        Self {
            core: InterpolationCore::default(),
            min: 0.0,
            max: 1.0,
            statistics1: statistics,
            statistics2: statistics,
        }
    }

    /// A model with the default parameters, sampled (source default
    /// constructor).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn new() -> Result<Self> {
        let mut model = Self::default_state();
        model.set_samples()?;
        Ok(model)
    }

    /// The default parameters (source `getDefaults`): the three shared ones,
    /// then the advanced `bounding_box:min` 0, `bounding_box:max` 1,
    /// `statistics:mean` 0, `statistics:variance1` 1 and
    /// `statistics:variance2` 1.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default_state().to_param_impl()
    }

    /// A default model after `set_parameters` with `param`.
    ///
    /// # Errors
    ///
    /// As `set_parameters`.
    pub fn from_param(param: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(param)?.0)
    }

    /// As [`BiGaussModel::from_param`], also returning the unknown-parameter
    /// warnings.
    ///
    /// # Errors
    ///
    /// As `set_parameters`.
    pub fn from_param_with_warnings(param: &Param) -> Result<(Self, Vec<String>)> {
        let mut model = Self::new()?;
        let warnings = model.set_parameters_impl(param)?;
        Ok((model, warnings))
    }

    /// The bounding box `(min, max)` the model is sampled over.
    pub fn bounding_box(&self) -> (f64, f64) {
        (self.min, self.max)
    }

    /// The statistics of the lower and the upper half.
    pub fn statistics(&self) -> (&BasicStatistics, &BasicStatistics) {
        (&self.statistics1, &self.statistics2)
    }

    fn to_param_impl(&self) -> Result<Param> {
        let mut p = Param::new();
        self.core.write_param(&mut p)?;
        set_default(
            &mut p,
            "bounding_box:min",
            ParamValue::Float(self.min),
            "Lower end of bounding box enclosing the data used to fit the model.",
            true,
        )?;
        set_default(
            &mut p,
            "bounding_box:max",
            ParamValue::Float(self.max),
            "Upper end of bounding box enclosing the data used to fit the model.",
            true,
        )?;
        set_default(
            &mut p,
            "statistics:mean",
            ParamValue::Float(self.statistics1.mean()),
            "Centroid position of the model, this also separates both halves of the model.",
            true,
        )?;
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

    fn set_parameters_impl(&mut self, param: &Param) -> Result<Vec<String>> {
        let (complete, warnings) = complete_parameters(Self::NAME, Self::defaults()?, param)?;
        let mut staged = self.clone();
        staged.core.read_param(&complete)?;
        staged.min = float_value(&complete, "bounding_box:min")?;
        staged.max = float_value(&complete, "bounding_box:max")?;
        let mean = float_value(&complete, "statistics:mean")?;
        staged.statistics1.set_mean(mean);
        staged.statistics2.set_mean(mean);
        staged
            .statistics1
            .set_variance(float_value(&complete, "statistics:variance1")?);
        staged
            .statistics2
            .set_variance(float_value(&complete, "statistics:variance2")?);
        staged.set_samples()?;
        *self = staged;
        Ok(warnings)
    }
}

interpolation_base_model!(BiGaussModel);

impl InterpolationModel for BiGaussModel {
    fn core(&self) -> &InterpolationCore {
        &self.core
    }
    fn core_mut(&mut self) -> &mut InterpolationCore {
        &mut self.core
    }

    /// The position of the maximum: the mean of the upper half (both halves
    /// share it).
    fn center(&self) -> Result<f64> {
        Ok(self.statistics2.mean())
    }

    /// Shift the model without recomputing it: the bounding box and both
    /// means move by `offset - interpolation().offset()`.
    fn set_offset(&mut self, offset: f64) {
        let diff = offset - self.core.interpolation.offset();
        self.min += diff;
        self.max += diff;
        let mean1 = self.statistics1.mean() + diff;
        self.statistics1.set_mean(mean1);
        let mean2 = self.statistics2.mean() + diff;
        self.statistics2.set_mean(mean2);
        self.core.interpolation.set_offset(offset);
    }

    /// Sample the lower Gaussian below the mean and the upper one from the
    /// mean on, then scale the samples to integrate to the scaling factor.
    ///
    /// # Errors
    ///
    /// As [`GaussModel`]'s `set_samples`.
    fn set_samples(&mut self) -> Result<()> {
        let (lower, upper) = (self.statistics1, self.statistics2);
        let data = sample_box(&self.core, self.min, self.max, |pos| {
            if pos < lower.mean() {
                lower.normal_density_sqrt2pi(pos)
            } else {
                upper.normal_density_sqrt2pi(pos)
            }
        })?
        .map(|mut data| {
            normalize(&mut data, self.core.scaling, self.core.interpolation_step);
            data
        });
        install(&mut self.core, data, self.min);
        Ok(())
    }
}

/// Exponentially modified Gaussian model of an elution profile (source
/// `EmgModel`).
///
/// The samples are the simplified EMG of Equation 9 of doi 10.1.1.915.3568:
/// `h * w / s * sqrt(2 pi) * exp(w^2 / (2 s^2) - d / s) /
/// (1 + exp(-2.4055 / sqrt(2) * (d / w - w / s)))` with `d = pos - retention`,
/// height `h`, width `w` and symmetry `s`. They are not normalized: the
/// scaling factor has no effect on this model.
///
/// The mean and variance parameters describe the position only for
/// [`InterpolationModel::center`] and [`InterpolationModel::set_offset`]; the
/// function itself is placed by `emg:retention`.
#[derive(Clone, Debug, PartialEq)]
pub struct EmgModel {
    core: InterpolationCore,
    /// Parameter `bounding_box:min`, default 0.
    min: f64,
    /// Parameter `bounding_box:max`, default 1.
    max: f64,
    /// Parameters `statistics:mean` (default 0) and `statistics:variance`
    /// (default 1).
    statistics: BasicStatistics,
    /// Parameter `emg:height`, default 100000.
    height: f64,
    /// Parameter `emg:width`, default 5.
    width: f64,
    /// Parameter `emg:symmetry`, default 5.
    symmetry: f64,
    /// Parameter `emg:retention`, default 1200.
    retention: f64,
}

impl EmgModel {
    /// The handler name (source `setName("EmgModel")`).
    pub const NAME: &'static str = "EmgModel";

    fn default_state() -> Self {
        let mut statistics = BasicStatistics::new();
        statistics.set_mean(0.0);
        statistics.set_variance(1.0);
        Self {
            core: InterpolationCore::default(),
            min: 0.0,
            max: 1.0,
            statistics,
            height: 100_000.0,
            width: 5.0,
            symmetry: 5.0,
            retention: 1200.0,
        }
    }

    /// A model with the default parameters, sampled (source default
    /// constructor).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn new() -> Result<Self> {
        let mut model = Self::default_state();
        model.set_samples()?;
        Ok(model)
    }

    /// The default parameters (source `getDefaults`): the three shared ones,
    /// then the advanced `bounding_box:min` 0, `bounding_box:max` 1,
    /// `statistics:mean` 0, `statistics:variance` 1, `emg:height` 100000,
    /// `emg:width` 5, `emg:symmetry` 5 and `emg:retention` 1200. The source
    /// writes these as `float` literals, all exactly representable.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default_state().to_param_impl()
    }

    /// A default model after `set_parameters` with `param`.
    ///
    /// # Errors
    ///
    /// As `set_parameters`.
    pub fn from_param(param: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(param)?.0)
    }

    /// As [`EmgModel::from_param`], also returning the unknown-parameter
    /// warnings.
    ///
    /// # Errors
    ///
    /// As `set_parameters`.
    pub fn from_param_with_warnings(param: &Param) -> Result<(Self, Vec<String>)> {
        let mut model = Self::new()?;
        let warnings = model.set_parameters_impl(param)?;
        Ok((model, warnings))
    }

    /// The bounding box `(min, max)` the model is sampled over.
    pub fn bounding_box(&self) -> (f64, f64) {
        (self.min, self.max)
    }

    /// The EMG parameters `(height, width, symmetry, retention)`.
    pub fn emg_parameters(&self) -> (f64, f64, f64, f64) {
        (self.height, self.width, self.symmetry, self.retention)
    }

    fn to_param_impl(&self) -> Result<Param> {
        let mut p = Param::new();
        self.core.write_param(&mut p)?;
        let entries: [(&str, f64, &str); 8] = [
            (
                "bounding_box:min",
                self.min,
                "Lower end of bounding box enclosing the data used to fit the model.",
            ),
            (
                "bounding_box:max",
                self.max,
                "Upper end of bounding box enclosing the data used to fit the model.",
            ),
            (
                "statistics:mean",
                self.statistics.mean(),
                "Centroid position of the model.",
            ),
            (
                "statistics:variance",
                self.statistics.variance(),
                "The variance of the model.",
            ),
            (
                "emg:height",
                self.height,
                "Height of the exponentially modified Gaussian.",
            ),
            (
                "emg:width",
                self.width,
                "Width of the exponentially modified Gaussian.",
            ),
            (
                "emg:symmetry",
                self.symmetry,
                "Symmetry of the exponentially modified Gaussian.",
            ),
            (
                "emg:retention",
                self.retention,
                "Retention time of the exponentially modified Gaussian.",
            ),
        ];
        for (key, value, description) in entries {
            set_default(&mut p, key, ParamValue::Float(value), description, true)?;
        }
        Ok(p)
    }

    fn set_parameters_impl(&mut self, param: &Param) -> Result<Vec<String>> {
        let (complete, warnings) = complete_parameters(Self::NAME, Self::defaults()?, param)?;
        let mut staged = self.clone();
        staged.core.read_param(&complete)?;
        staged.min = float_value(&complete, "bounding_box:min")?;
        staged.max = float_value(&complete, "bounding_box:max")?;
        staged
            .statistics
            .set_mean(float_value(&complete, "statistics:mean")?);
        staged
            .statistics
            .set_variance(float_value(&complete, "statistics:variance")?);
        staged.height = float_value(&complete, "emg:height")?;
        staged.width = float_value(&complete, "emg:width")?;
        staged.symmetry = float_value(&complete, "emg:symmetry")?;
        staged.retention = float_value(&complete, "emg:retention")?;
        staged.set_samples()?;
        *self = staged;
        Ok(warnings)
    }
}

interpolation_base_model!(EmgModel);

/// The simplified EMG at `pos`, in the source's operation order.
///
/// `exp` is the reference build's GNU C library `exp`
/// ([`glibc_libm::exp`]), so the samples match the Linux Release build bit
/// for bit. The source's `pow(width, 2)` and `pow(symmetry, 2)` are computed
/// as products: glibc's `pow` is accurate to 0.52 ulp, so the two can differ
/// in the last bit where the exact square lies within 0.02 ulp of a rounding
/// midpoint.
pub(crate) fn emg_value(height: f64, width: f64, symmetry: f64, retention: f64, pos: f64) -> f64 {
    let sqrt_2pi = (2.0 * PI).sqrt();
    let term_sq2 = -2.4055 / 2.0f64.sqrt();
    let part1 = height * width / symmetry;
    let part2 = (width * width) / (2.0 * (symmetry * symmetry));
    let part3 = width / symmetry;
    let diff = pos - retention;
    part1 * sqrt_2pi * glibc_libm::exp(part2 - (diff / symmetry))
        / (1.0 + glibc_libm::exp(term_sq2 * ((diff / width) - part3)))
}

impl InterpolationModel for EmgModel {
    fn core(&self) -> &InterpolationCore {
        &self.core
    }
    fn core_mut(&mut self) -> &mut InterpolationCore {
        &mut self.core
    }

    /// The `statistics:mean` parameter, which [`InterpolationModel::set_offset`]
    /// moves; not the EMG's maximum.
    fn center(&self) -> Result<f64> {
        Ok(self.statistics.mean())
    }

    /// Shift the model without recomputing it: the bounding box and the mean
    /// move by `offset - interpolation().offset()`.
    ///
    /// As in the source, `emg:retention` does not move, so a model rebuilt from
    /// [`crate::analysis::featurefinder::base_model::BaseModel::to_param`]
    /// after a shift samples the EMG at its old place within the moved box.
    fn set_offset(&mut self, offset: f64) {
        let diff = offset - self.core.interpolation.offset();
        self.min += diff;
        self.max += diff;
        let mean = self.statistics.mean() + diff;
        self.statistics.set_mean(mean);
        self.core.interpolation.set_offset(offset);
    }

    /// Sample the simplified EMG over the bounding box, without normalizing.
    ///
    /// # Errors
    ///
    /// As [`GaussModel`]'s `set_samples`. Overflowing exponentials are not
    /// refused: they yield infinite or NaN samples, as in the source.
    fn set_samples(&mut self) -> Result<()> {
        let (h, w, s, r) = (self.height, self.width, self.symmetry, self.retention);
        let data = sample_box(&self.core, self.min, self.max, |pos| {
            emg_value(h, w, s, r, pos)
        })?;
        install(&mut self.core, data, self.min);
        Ok(())
    }
}
