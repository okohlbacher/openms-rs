// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The one-dimensional model interfaces of the feature finder.
//!
//! Ports `FEATUREFINDER/BaseModel.h` (with its empty `BaseModel.cpp`),
//! `FEATUREFINDER/BaseModel_impl.h`, `FEATUREFINDER/ModelDescription.h` and
//! `FEATUREFINDER/InterpolationModel.h` (with its empty
//! `InterpolationModel.cpp`). See `docs/FEATUREFINDER_MODELS_SUPPORT.md`.
//!
//! The source's class hierarchy (`DefaultParamHandler` <- `BaseModel` <-
//! `InterpolationModel` <- five concrete models) becomes two traits,
//! [`BaseModel`] and [`InterpolationModel`], the shared state
//! [`InterpolationCore`], and the closed set [`Model1D`] that the fitters
//! return in place of `std::unique_ptr<InterpolationModel>`.
//!
//! `BaseModel_impl.h` only includes the five derived model headers so that a
//! caller can name them all through one include; [`Model1D`] is its
//! counterpart. It declares nothing else.

use crate::Result;
use crate::error::Error;
use crate::kernel::Peak1D;
use crate::math::ml::linear_interpolation::LinearInterpolation;
use crate::param::{DefaultParamHandler, Param, ParamValue};

use super::gauss_models::{BiGaussModel, EmgModel, GaussModel};
use super::isotope_model::{ExtendedIsotopeModel, IsotopeModel};

/// Most interpolation samples one model may compute.
///
/// Native-only ceiling, checked before a sampling loop allocates: the source
/// samples whatever its bounding box and step describe, and loops forever on a
/// zero or negative step.
pub const MAX_SAMPLES: usize = 10_000_000;

/// The `advanced` tag every model and fitter parameter except the three
/// [`InterpolationCore`] ones carries.
pub(crate) const ADVANCED: &str = "advanced";

/// Set one default the way the source's `defaults_.setValue` does.
pub(crate) fn set_default(
    param: &mut Param,
    key: &str,
    value: ParamValue,
    description: &str,
    advanced: bool,
) -> Result<()> {
    let tags = if advanced {
        vec![ADVANCED.to_owned()]
    } else {
        Vec::new()
    };
    param.set_value(key, value, description, &tags)
}

/// Restrict a string default to `valid` (source `setValidStrings`).
pub(crate) fn set_valid(param: &mut Param, key: &str, valid: &[&str]) -> Result<()> {
    let strings: Vec<String> = valid.iter().map(|s| (*s).to_owned()).collect();
    param.set_valid_strings(key, &strings)
}

/// A float parameter as the source's `(double)param_.getValue(key)` reads it:
/// integers convert, anything else is a `ConversionError`.
pub(crate) fn float_value(param: &Param, key: &str) -> Result<f64> {
    param.value(key)?.to_f64()
}

/// Validate `parameters` against `defaults` for the handler `name`, as the
/// source's `setParameters` does before `updateMembers_`, and return the
/// completed tree with the unknown-parameter warnings.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] for a value of the wrong type or outside its
/// restriction (the source `Exception::InvalidParameter`).
pub(crate) fn complete_parameters(
    name: &str,
    defaults: Param,
    parameters: &Param,
) -> Result<(Param, Vec<String>)> {
    let mut handler = DefaultParamHandler::new(name)?;
    handler.set_defaults(defaults)?;
    handler.set_parameters_with(parameters, |complete| complete.checked_clone())
}

/// A one-dimensional model: predicted intensity by position (source abstract
/// class `BaseModel`, a `DefaultParamHandler`).
///
/// The source's copy constructor and assignment operator are `Clone` on each
/// implementor; its `getParameters` is [`BaseModel::to_param`] and its
/// `setParameters` is [`BaseModel::set_parameters`].
pub trait BaseModel {
    /// The handler name (source `getName`): the model's class name.
    fn name(&self) -> &'static str;

    /// The current parameters, in the layout of the model's defaults (source
    /// `getParameters`).
    ///
    /// Unlike the source, a parameter the defaults do not know is not kept:
    /// the source warns about it and carries it along unused.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    fn to_param(&self) -> Result<Param>;

    /// Replace the parameters, taking the defaults for the missing ones, and
    /// recompute the model (source `setParameters` followed by
    /// `updateMembers_`). Returns the source's unknown-parameter warnings.
    ///
    /// As in the source, every earlier setter call - including
    /// [`InterpolationModel::set_interpolation_step`] - is forgotten unless
    /// `param` repeats it, because the tree is replaced, not merged.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a parameter of the wrong type or
    /// outside its restriction, and whatever the model's recomputation
    /// refuses. The model is then unchanged; the source has already replaced
    /// its parameter tree when the check throws.
    fn set_parameters(&mut self, param: &Param) -> Result<Vec<String>>;

    /// The predicted intensity at position `pos` (source `getIntensity`).
    fn intensity(&self, pos: f64) -> f64;

    /// Whether `pos` is part of the model: its intensity reaches the cut-off
    /// (source `isContained`).
    fn is_contained(&self, pos: f64) -> bool {
        self.intensity(pos) >= self.cut_off()
    }

    /// Set a peak's intensity to the prediction at its position (source
    /// `fillIntensity`). The prediction is narrowed to the peak's `f32`
    /// intensity, as the source's `setIntensity` does.
    fn fill_intensity(&self, peak: &mut Peak1D) {
        peak.intensity = self.intensity(peak.mz) as f32;
    }

    /// [`BaseModel::fill_intensity`] for every peak (source `fillIntensities`
    /// over an iterator range).
    fn fill_intensities(&self, peaks: &mut [Peak1D]) {
        for peak in peaks {
            self.fill_intensity(peak);
        }
    }

    /// The low-intensity cut-off (source `getCutOff`, parameter `cutoff`).
    fn cut_off(&self) -> f64;

    /// Set the cut-off; [`BaseModel::to_param`] reports it (source
    /// `setCutOff`).
    fn set_cut_off(&mut self, cut_off: f64);

    /// A reasonable set of samples of the model, for printing (source
    /// `getSamples(SamplesType&)`).
    fn samples(&self) -> Vec<Peak1D>;

    /// Write [`BaseModel::samples`], one per line, as the source's
    /// `getSamples(std::ostream&)` does with `Peak1D`'s stream operator:
    /// `POS: <mz> INT: <intensity>`, each number as `std::ostream` writes it by
    /// default (six significant digits).
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when writing fails.
    fn write_samples(&self, out: &mut dyn std::io::Write) -> Result<()> {
        use crate::concept::fuzzy_string_comparator::format_g;
        for sample in self.samples() {
            writeln!(
                out,
                "POS: {} INT: {}",
                format_g(sample.mz),
                format_g(f64::from(sample.intensity))
            )
            .map_err(Error::Io)?;
        }
        Ok(())
    }
}

/// The state every interpolation model shares (the members of source
/// `BaseModel` and `InterpolationModel`): the cut-off, the interpolation step,
/// the intensity scaling and the [`LinearInterpolation`] of the model
/// function.
#[derive(Clone, Debug, PartialEq)]
pub struct InterpolationCore {
    /// Parameter `cutoff` (source `cut_off_`), default 0.
    pub cut_off: f64,
    /// Parameter `interpolation_step` (source `interpolation_step_`),
    /// default 0.1: the sampling rate of the model function.
    pub interpolation_step: f64,
    /// Parameter `intensity_scaling` (source `scaling_`), default 1: the area
    /// under the model.
    pub scaling: f64,
    /// The sampled model function (source `interpolation_`).
    pub interpolation: LinearInterpolation,
}

impl Default for InterpolationCore {
    /// The source defaults with an empty, identity-mapped interpolation.
    fn default() -> Self {
        Self {
            cut_off: 0.0,
            interpolation_step: 0.1,
            scaling: 1.0,
            interpolation: LinearInterpolation::default(),
        }
    }
}

impl InterpolationCore {
    /// Write the three shared defaults, with the current values, in the
    /// source's order: `cutoff` (from `BaseModel`), then `interpolation_step`
    /// and `intensity_scaling` (from `InterpolationModel`). None is advanced.
    pub(crate) fn write_param(&self, param: &mut Param) -> Result<()> {
        set_default(
            param,
            "cutoff",
            ParamValue::Float(self.cut_off),
            "Low intensity cutoff of the model.  Peaks below this intensity are not considered part of the model.",
            false,
        )?;
        set_default(
            param,
            "interpolation_step",
            ParamValue::Float(self.interpolation_step),
            "Sampling rate for the interpolation of the model function ",
            false,
        )?;
        set_default(
            param,
            "intensity_scaling",
            ParamValue::Float(self.scaling),
            "Scaling factor used to adjust the model distribution to the intensities of the data",
            false,
        )
    }

    /// Read the three shared parameters (source `BaseModel::updateMembers_`
    /// and `InterpolationModel::updateMembers_`), keeping the interpolation.
    pub(crate) fn read_param(&mut self, param: &Param) -> Result<()> {
        self.cut_off = float_value(param, "cutoff")?;
        self.interpolation_step = float_value(param, "interpolation_step")?;
        self.scaling = float_value(param, "intensity_scaling")?;
        Ok(())
    }

    /// The interpolation's samples as peaks: position `index2key(i)`, intensity
    /// the sample narrowed to `f32` (source `InterpolationModel::getSamples`).
    pub fn samples(&self) -> Vec<Peak1D> {
        self.interpolation
            .data()
            .iter()
            .enumerate()
            .map(|(i, &value)| Peak1D::new(self.interpolation.index2key(i as f64), value as f32))
            .collect()
    }
}

/// A model approximated by linear interpolation of its sampled function
/// (source abstract class `InterpolationModel`, "derived classes have to
/// implement setSamples()").
///
/// Its predicted intensity is the interpolation's value, so
/// [`BaseModel::intensity`] of every implementor is
/// `interpolation().value(pos)`.
pub trait InterpolationModel: BaseModel {
    /// The shared state.
    fn core(&self) -> &InterpolationCore;

    /// Mutable access to the shared state, for the provided setters.
    fn core_mut(&mut self) -> &mut InterpolationCore;

    /// The model's "center", whose definition depends on the model (source
    /// `getCenter`).
    ///
    /// # Errors
    ///
    /// The source base class throws `Exception::NotImplemented`; every model
    /// ported here defines it, so this is infallible in practice.
    fn center(&self) -> Result<f64>;

    /// Shift the whole model to `offset` without recomputing it (source
    /// `setOffset`). The default moves the interpolation's offset; each model
    /// also moves the parameters that describe its position.
    fn set_offset(&mut self, offset: f64) {
        self.core_mut().interpolation.set_offset(offset);
    }

    /// Recompute the interpolation's samples from the parameters (source
    /// `setSamples`).
    ///
    /// # Errors
    ///
    /// The source base class throws `Exception::NotImplemented`, which is what
    /// [`crate::analysis::featurefinder::isotope_model::IsotopeModel`] inherits
    /// ([`Error::Unsupported`]); the other models return what their sampling
    /// refuses.
    fn set_samples(&mut self) -> Result<()>;

    /// The interpolation of the model function (source `getInterpolation`).
    fn interpolation(&self) -> &LinearInterpolation {
        &self.core().interpolation
    }

    /// The scaling: the area under the model (source `getScalingFactor`).
    fn scaling_factor(&self) -> f64 {
        self.core().scaling
    }

    /// Set the scaling; call [`InterpolationModel::set_samples`] for it to take
    /// effect (source `setScalingFactor`).
    fn set_scaling_factor(&mut self, scaling: f64) {
        self.core_mut().scaling = scaling;
    }

    /// The interpolation step (parameter `interpolation_step`).
    fn interpolation_step(&self) -> f64 {
        self.core().interpolation_step
    }

    /// Set the interpolation step; call [`InterpolationModel::set_samples`]
    /// for it to take effect (source `setInterpolationStep`).
    fn set_interpolation_step(&mut self, interpolation_step: f64) {
        self.core_mut().interpolation_step = interpolation_step;
    }
}

/// Every model the feature finder builds (source `BaseModel_impl.h`, and the
/// dynamic type behind a fitter's `std::unique_ptr<InterpolationModel>`).
#[derive(Clone, Debug, PartialEq)]
pub enum Model1D {
    /// [`GaussModel`].
    Gauss(GaussModel),
    /// [`BiGaussModel`].
    BiGauss(BiGaussModel),
    /// [`EmgModel`].
    Emg(EmgModel),
    /// [`IsotopeModel`].
    Isotope(IsotopeModel),
    /// [`ExtendedIsotopeModel`].
    ExtendedIsotope(ExtendedIsotopeModel),
}

macro_rules! dispatch {
    ($self:ident, $model:ident => $body:expr) => {
        match $self {
            Model1D::Gauss($model) => $body,
            Model1D::BiGauss($model) => $body,
            Model1D::Emg($model) => $body,
            Model1D::Isotope($model) => $body,
            Model1D::ExtendedIsotope($model) => $body,
        }
    };
}

impl BaseModel for Model1D {
    fn name(&self) -> &'static str {
        dispatch!(self, m => m.name())
    }
    fn to_param(&self) -> Result<Param> {
        dispatch!(self, m => m.to_param())
    }
    fn set_parameters(&mut self, param: &Param) -> Result<Vec<String>> {
        dispatch!(self, m => m.set_parameters(param))
    }
    fn intensity(&self, pos: f64) -> f64 {
        dispatch!(self, m => m.intensity(pos))
    }
    fn cut_off(&self) -> f64 {
        dispatch!(self, m => m.cut_off())
    }
    fn set_cut_off(&mut self, cut_off: f64) {
        dispatch!(self, m => m.set_cut_off(cut_off))
    }
    fn samples(&self) -> Vec<Peak1D> {
        dispatch!(self, m => m.samples())
    }
}

impl InterpolationModel for Model1D {
    fn core(&self) -> &InterpolationCore {
        dispatch!(self, m => m.core())
    }
    fn core_mut(&mut self) -> &mut InterpolationCore {
        dispatch!(self, m => m.core_mut())
    }
    fn center(&self) -> Result<f64> {
        dispatch!(self, m => m.center())
    }
    fn set_offset(&mut self, offset: f64) {
        dispatch!(self, m => m.set_offset(offset))
    }
    fn set_samples(&mut self) -> Result<()> {
        dispatch!(self, m => m.set_samples())
    }
}

/// Implement [`BaseModel`] for an interpolation model with a `core` field, a
/// `NAME` constant and `to_param_impl`/`set_parameters_impl` methods.
macro_rules! interpolation_base_model {
    ($ty:ty) => {
        impl $crate::analysis::featurefinder::base_model::BaseModel for $ty {
            fn name(&self) -> &'static str {
                Self::NAME
            }
            fn to_param(&self) -> $crate::Result<$crate::param::Param> {
                self.to_param_impl()
            }
            fn set_parameters(
                &mut self,
                param: &$crate::param::Param,
            ) -> $crate::Result<Vec<String>> {
                self.set_parameters_impl(param)
            }
            fn intensity(&self, pos: f64) -> f64 {
                self.core.interpolation.value(pos)
            }
            fn cut_off(&self) -> f64 {
                self.core.cut_off
            }
            fn set_cut_off(&mut self, cut_off: f64) {
                self.core.cut_off = cut_off;
            }
            fn samples(&self) -> Vec<$crate::kernel::Peak1D> {
                self.core.samples()
            }
        }
    };
}
pub(crate) use interpolation_base_model;

/// The name and parameters of a model, from which the model can be rebuilt
/// (source `ModelDescription<D>`).
///
/// The source's dimension template parameter `D` is unused by every member and
/// is not carried over.
#[derive(Clone, Debug, Default)]
pub struct ModelDescription {
    /// The model name (source `getName`/`setName`).
    pub name: String,
    /// The model parameters (source `getParam`/`setParam`).
    pub param: Param,
}

impl ModelDescription {
    /// The description of `model`: its name and current parameters (source
    /// constructor from `const BaseModel*`).
    ///
    /// # Errors
    ///
    /// As [`BaseModel::to_param`].
    pub fn from_model<M: BaseModel + ?Sized>(model: &M) -> Result<Self> {
        Ok(Self {
            name: model.name().to_owned(),
            param: model.to_param()?,
        })
    }

    /// Source `operator==`: equal names and parameters under the source's
    /// `Param` equality (names and values).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn source_equal(&self, other: &Self) -> Result<bool> {
        Ok(self.name == other.name && self.param.source_equal(&other.param)?)
    }

    /// Rebuild the described model (the reconstruction the source's class
    /// documentation promises; the source itself has no factory for it).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a name that is not one of the five
    /// interpolation models, and whatever the model's
    /// [`BaseModel::set_parameters`] refuses.
    pub fn create_model(&self) -> Result<Model1D> {
        Ok(match self.name.as_str() {
            GaussModel::NAME => Model1D::Gauss(GaussModel::from_param(&self.param)?),
            BiGaussModel::NAME => Model1D::BiGauss(BiGaussModel::from_param(&self.param)?),
            EmgModel::NAME => Model1D::Emg(EmgModel::from_param(&self.param)?),
            IsotopeModel::NAME => Model1D::Isotope(IsotopeModel::from_param(&self.param)?),
            ExtendedIsotopeModel::NAME => {
                Model1D::ExtendedIsotope(ExtendedIsotopeModel::from_param(&self.param)?)
            }
            other => {
                return Err(Error::InvalidValue(format!(
                    "'{other}' is not a feature finder model"
                )));
            }
        })
    }
}
