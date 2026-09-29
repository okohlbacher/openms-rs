// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The one-dimensional model fitter interface and its two optimisation
//! strategies.
//!
//! Ports `FEATUREFINDER/Fitter1D.h`, `FEATUREFINDER/LevMarqFitter1D.h` and
//! `FEATUREFINDER/MaxLikeliFitter1D.h` with their `.cpp` files. See
//! `docs/FEATUREFINDER_MODELS_SUPPORT.md`.
//!
//! - Source `Fitter1D`, an instantiable `DefaultParamHandler` whose `fit1d`
//!   throws `Exception::NotImplemented`, is the trait [`Fitter1D`](crate::analysis::featurefinder::fitter1d::Fitter1D), its shared
//!   parameters [`Fitter1DCore`](crate::analysis::featurefinder::fitter1d::Fitter1DCore), and the base type [`BaseFitter1D`](crate::analysis::featurefinder::fitter1d::BaseFitter1D).
//! - Source `LevMarqFitter1D` contributes the parameter `max_iteration` and the
//!   protected `optimize_`, which is [`lev_marq_optimize`](crate::analysis::featurefinder::fitter1d::lev_marq_optimize). Its
//!   `GenericFunctor` (raw-pointer residual and Jacobian callbacks) becomes the
//!   two closures that function takes.
//! - Source `MaxLikeliFitter1D` contributes the protected `fitOffset_`, which
//!   is [`max_likeli_fit_offset`](crate::analysis::featurefinder::fitter1d::max_likeli_fit_offset).
//!
//! A fit's quality is the source's `Feature::QualityType`, a `float`: every
//! correlation is narrowed to `f32` where the source assigns it, and
//! comparisons between qualities happen in `f32`, as there.

use crate::Result;
use crate::analysis::feature_finder_picked::trace_fitter::{MAX_RESIDUAL_WORK, unable_to_fit};
use crate::error::Error;
use crate::kernel::Peak1D;
use crate::math::basic_statistics::BasicStatistics;
use crate::math::fitters::levenberg_marquardt::{
    DenseMatrix, LmParameters, LmStatus, minimize, preflight_points,
};
use crate::math::statistic_functions::pearson_correlation_coefficient;
use crate::param::{Param, ParamValue};

use super::base_model::Model1D;
use super::base_model::{
    InterpolationModel, MAX_SAMPLES, complete_parameters, float_value, set_default,
};

/// The quality of a fit (source `Fitter1D::QualityType`, which is
/// `Feature::QualityType`, a `float`).
pub type QualityType = f32;

/// The message `LevMarqFitter1D::optimize_` throws when there are fewer data
/// points than parameters, verbatim (it differs from `TraceFitter`'s by one
/// letter).
pub const FEWER_POINTS_THAN_PARAMETERS: &str = "Skipping feature, we always expects N>=p";

/// Most offsets times data points one [`max_likeli_fit_offset`] may evaluate.
///
/// Native-only ceiling: the source scans `[support_min - stdev1, support_min +
/// stdev2]` at `offset_step` without bound, and forever for a non-positive
/// step.
pub const MAX_OFFSET_WORK: usize = 1 << 32;

/// A one-dimensional model fitter (source abstract base `Fitter1D`, a
/// `DefaultParamHandler`).
///
/// The source's copy constructor and assignment operator are `Clone` on each
/// implementor.
pub trait Fitter1D {
    /// The handler name (source `getName`).
    fn name(&self) -> &'static str;

    /// The current parameters (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    fn to_param(&self) -> Result<Param>;

    /// Replace the parameters, taking the defaults for the missing ones
    /// (source `setParameters`). Returns the unknown-parameter warnings.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a parameter of the wrong type or
    /// outside its restriction; the fitter is then unchanged.
    fn set_parameters(&mut self, param: &Param) -> Result<Vec<String>>;

    /// Fit a model to `range` and return its quality and the model (source
    /// `fit1d(range, model)`, which returns the quality and sets `model`).
    ///
    /// # Errors
    ///
    /// As each fitter documents; [`BaseFitter1D`] returns
    /// [`Error::Unsupported`], as the source base class throws
    /// `Exception::NotImplemented`.
    fn fit1d(&self, range: &[Peak1D]) -> Result<(QualityType, Model1D)>;
}

/// The parameters every fitter shares (source `Fitter1D`'s members).
#[derive(Clone, Debug, PartialEq)]
pub struct Fitter1DCore {
    /// Parameter `interpolation_step` (source `interpolation_step_`),
    /// default 0.2: the sampling rate of the model function.
    pub interpolation_step: f64,
    /// Parameters `statistics:mean` (default 1) and `statistics:variance`
    /// (default 1) (source `statistics_`).
    pub statistics: BasicStatistics,
    /// Parameter `tolerance_stdev_bounding_box` (source `tolerance_stdev_box_`),
    /// default 3: the bounding box spans the data enlarged by this many
    /// standard deviations.
    pub tolerance_stdev_box: f64,
}

impl Default for Fitter1DCore {
    /// The source defaults.
    fn default() -> Self {
        let mut statistics = BasicStatistics::new();
        statistics.set_mean(1.0);
        statistics.set_variance(1.0);
        Self {
            interpolation_step: 0.2,
            statistics,
            tolerance_stdev_box: 3.0,
        }
    }
}

impl Fitter1DCore {
    /// Write the four shared parameters with the current values, in the
    /// source's order, all advanced: `interpolation_step`, `statistics:mean`,
    /// `statistics:variance`, `tolerance_stdev_bounding_box`.
    pub(crate) fn write_param(&self, param: &mut Param) -> Result<()> {
        set_default(
            param,
            "interpolation_step",
            ParamValue::Float(self.interpolation_step),
            "Sampling rate for the interpolation of the model function.",
            true,
        )?;
        set_default(
            param,
            "statistics:mean",
            ParamValue::Float(self.statistics.mean()),
            "Centroid position of the model.",
            true,
        )?;
        set_default(
            param,
            "statistics:variance",
            ParamValue::Float(self.statistics.variance()),
            "The variance of the model.",
            true,
        )?;
        set_default(
            param,
            "tolerance_stdev_bounding_box",
            ParamValue::Float(self.tolerance_stdev_box),
            "Bounding box has range [minimim of data, maximum of data] enlarged by tolerance_stdev_bounding_box times the standard deviation of the data.",
            true,
        )
    }

    /// Read the four shared parameters (source `Fitter1D::updateMembers_`).
    pub(crate) fn read_param(&mut self, param: &Param) -> Result<()> {
        self.tolerance_stdev_box = float_value(param, "tolerance_stdev_bounding_box")?;
        self.interpolation_step = float_value(param, "interpolation_step")?;
        self.statistics
            .set_mean(float_value(param, "statistics:mean")?);
        self.statistics
            .set_variance(float_value(param, "statistics:variance")?);
        Ok(())
    }

    /// The bounding box of `set`'s positions, as the source's loop: start at
    /// the first position and replace the minimum (maximum) only by a smaller
    /// (larger) one, so NaN positions after the first are ignored.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] for an empty set, where the source reads
    /// `set[0]` out of range.
    pub(crate) fn data_range(set: &[Peak1D]) -> Result<(f64, f64)> {
        let first = set
            .first()
            .ok_or_else(|| Error::InvalidRange("a fitter needs at least one data point".into()))?;
        let (mut min_bb, mut max_bb) = (first.mz, first.mz);
        for peak in &set[1..] {
            let tmp = peak.mz;
            if min_bb > tmp {
                min_bb = tmp;
            }
            if max_bb < tmp {
                max_bb = tmp;
            }
        }
        Ok((min_bb, max_bb))
    }
}

/// The source class `Fitter1D` itself: the shared parameters and a `fit1d`
/// that is not implemented.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BaseFitter1D {
    /// The shared parameters.
    pub core: Fitter1DCore,
}

impl BaseFitter1D {
    /// The handler name (source `DefaultParamHandler("Fitter1D")`).
    pub const NAME: &'static str = "Fitter1D";

    /// The default parameters (source `getDefaults`): `interpolation_step`
    /// 0.2, `statistics:mean` 1, `statistics:variance` 1 and
    /// `tolerance_stdev_bounding_box` 3, all advanced.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        let mut p = Param::new();
        Fitter1DCore::default().write_param(&mut p)?;
        Ok(p)
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

impl Fitter1D for BaseFitter1D {
    fn name(&self) -> &'static str {
        Self::NAME
    }
    fn to_param(&self) -> Result<Param> {
        let mut p = Param::new();
        self.core.write_param(&mut p)?;
        Ok(p)
    }
    fn set_parameters(&mut self, param: &Param) -> Result<Vec<String>> {
        let (complete, warnings) = complete_parameters(Self::NAME, Self::defaults()?, param)?;
        let mut core = self.core.clone();
        core.read_param(&complete)?;
        self.core = core;
        Ok(warnings)
    }
    /// Always [`Error::Unsupported`]: the source throws
    /// `Exception::NotImplemented`.
    fn fit1d(&self, _range: &[Peak1D]) -> Result<(QualityType, Model1D)> {
        Err(Error::Unsupported(
            "Fitter1D::fit1d is not implemented (the source throws NotImplemented)".into(),
        ))
    }
}

/// Default of `LevMarqFitter1D`'s parameter `max_iteration`.
pub const MAX_ITERATION_DEFAULT: i32 = 500;

/// Add `LevMarqFitter1D`'s parameter `max_iteration` (Int, advanced) with
/// `value` to `param`.
pub(crate) fn write_max_iteration(param: &mut Param, value: i32) -> Result<()> {
    set_default(
        param,
        "max_iteration",
        ParamValue::Integer(i64::from(value)),
        "Maximum number of iterations using by Levenberg-Marquardt algorithm.",
        true,
    )
}

/// Optimize the start parameters `x` by Levenberg-Marquardt (source
/// `LevMarqFitter1D::optimize_`).
///
/// `values` is the number of residuals, one per data point, `residual` fills
/// them and `jacobian` fills the `values x x.len()` Jacobian (the source's
/// `GenericFunctor`). The solver is
/// [`crate::math::fitters::levenberg_marquardt::minimize`] with Eigen's
/// defaults and `maxfev = max_iteration` (a budget of residual evaluations,
/// despite the name). On success `x` holds the optimum.
///
/// # Errors
///
/// Returns the stand-in for the source's `Exception::UnableToFit`
/// (`UnableToFit-FinalSet`, see [`unable_to_fit`]) when `values < x.len()`
/// ([`FEWER_POINTS_THAN_PARAMETERS`]) or the solver ends with a status up to
/// and including `ImproperInputParameters` - in particular for
/// `max_iteration <= 0` - with the message `Could not fit the gaussian to the
/// data: Error <status>`. Returns [`Error::InvalidValue`] beyond the solver's
/// point ceiling, or when the fit spends
/// [`MAX_RESIDUAL_WORK`] residual-point evaluations before its budget; the
/// source has no ceiling. `x` is unchanged on every error.
pub fn lev_marq_optimize<R, J>(
    x: &mut [f64],
    values: usize,
    residual: R,
    mut jacobian: J,
    max_iteration: i32,
) -> Result<()>
where
    R: FnMut(&[f64], &mut [f64]),
    J: FnMut(&[f64], &mut DenseMatrix),
{
    if values < x.len() {
        return Err(unable_to_fit(FEWER_POINTS_THAN_PARAMETERS));
    }
    preflight_points(values, x.len())?;
    let improper = |status: LmStatus| {
        unable_to_fit(&format!(
            "Could not fit the gaussian to the data: Error {}",
            status.code()
        ))
    };
    if x.is_empty() || max_iteration <= 0 {
        return Err(improper(LmStatus::ImproperInputParameters));
    }
    let budget = usize::try_from(max_iteration).unwrap_or(usize::MAX);
    let ceiling = (MAX_RESIDUAL_WORK / values.max(1)).max(1);
    let max_fev = budget.min(ceiling);
    let mut working = x.to_vec();
    let status = minimize(
        &mut working,
        values,
        residual,
        |point: &[f64], jac: &mut DenseMatrix| {
            jacobian(point, jac);
            0
        },
        &LmParameters {
            max_fev,
            ..LmParameters::default()
        },
    );
    if status.code() <= LmStatus::ImproperInputParameters.code() {
        return Err(improper(status));
    }
    if status == LmStatus::TooManyFunctionEvaluation && max_fev < budget {
        return Err(Error::InvalidValue(format!(
            "the fit spent its work ceiling of {max_fev} residual evaluations over {values} \
             residuals before its configured budget of {budget}"
        )));
    }
    x.copy_from_slice(&working);
    Ok(())
}

/// The Pearson correlation of `set`'s intensities with `model`'s predictions
/// at its positions, as the source computes it: both series stored as `float`
/// (`std::vector<float>`), the coefficient in `double`, the result narrowed to
/// the `float` [`QualityType`].
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty set (the source's
/// `Exception::InvalidRange`).
pub(crate) fn model_correlation<M: InterpolationModel + ?Sized>(
    model: &M,
    set: &[Peak1D],
) -> Result<QualityType> {
    let real: Vec<f64> = set.iter().map(|p| f64::from(p.intensity)).collect();
    let predicted: Vec<f64> = set
        .iter()
        .map(|p| f64::from(model.intensity(p.mz) as f32))
        .collect();
    Ok(pearson_correlation_coefficient(&real, &predicted)? as QualityType)
}

/// Fit the model's offset on the Pearson correlation with the data (source
/// `MaxLikeliFitter1D::fitOffset_`).
///
/// The correlation at the model's current offset is the first candidate; then
/// every offset from `support_min - stdev1` to `support_min + stdev2` in steps
/// of `offset_step` (accumulated by addition) is tried with
/// [`InterpolationModel::set_offset`]. The model ends at the offset of the
/// first strictly best correlation, which is returned. Correlations compare as
/// `float`, and a NaN correlation never wins, so an all-NaN scan keeps the
/// starting offset and returns NaN.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] for an empty set, and
/// [`Error::InvalidValue`] when a non-empty scan has a non-positive or NaN
/// step or would exceed [`MAX_OFFSET_WORK`] (the source loops without bound,
/// or forever). The model is then unchanged.
pub fn max_likeli_fit_offset<M: InterpolationModel + ?Sized>(
    model: &mut M,
    set: &[Peak1D],
    stdev1: f64,
    stdev2: f64,
    offset_step: f64,
) -> Result<QualityType> {
    let support_min = model.interpolation().support_min();
    let offset_min = support_min - stdev1;
    let offset_max = support_min + stdev2;
    if offset_min <= offset_max {
        let steps = (offset_max - offset_min) / offset_step;
        let work = (steps + 2.0) * set.len().max(1) as f64;
        let bounded = offset_step > 0.0
            && steps.is_finite()
            && steps < MAX_SAMPLES as f64
            && work < MAX_OFFSET_WORK as f64;
        if !bounded {
            return Err(Error::InvalidValue(format!(
                "the offset scan from {offset_min} to {offset_max} at step {offset_step} never \
                 ends or exceeds its work ceiling"
            )));
        }
    }
    let start_offset = model.interpolation().offset();
    let mut max_offset = start_offset;
    let mut max_correlation = model_correlation(model, set)?;
    let mut offset = offset_min;
    while offset <= offset_max {
        model.set_offset(offset);
        let correlation = model_correlation(model, set)?;
        if correlation > max_correlation {
            max_correlation = correlation;
            max_offset = offset;
        }
        offset += offset_step;
    }
    model.set_offset(max_offset);
    Ok(max_correlation)
}
