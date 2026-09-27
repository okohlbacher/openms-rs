// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The `DefaultParamHandler` contracts of the filters in [`crate::processing`]
//! itself: [`NLargest`], [`ThresholdMower`], [`Normalizer`], [`RankScaler`],
//! [`SqrtScaler`] and [`LinearResamplerAlign`].
//!
//! Each source class derives from `DefaultParamHandler`; here each type
//! implements [`DefaultParamHandling`]. The names, defaults, descriptions and
//! restrictions are those of the source constructors at the pinned revision
//! (`NLargest.cpp:31-35`, `ThresholdMower.cpp:15-20`, `Normalizer.cpp:15-21`,
//! `LinearResamplerAlign.h:79-84`; `RankScaler` and `SqrtScaler` declare no
//! parameters). Tier 1: each [`DefaultParamHandling::defaults`] equals the
//! `getDefaults()` that the Linux x86-64 Release build writes through
//! `ParamXMLFile` (`tests/default_param_handling.rs`,
//! `docs/DEFAULT_PARAM_HANDLING_SUPPORT.md`).

use super::{LinearResamplerAlign, NLargest, NormalizationMethod, Normalizer, RankScaler};
use super::{SqrtScaler, ThresholdMower};
use crate::param::handler::{int_parameter, unsigned_parameter};
use crate::param::{DefaultParamHandling, Param, ParamValue};
use crate::{Error, Result};

/// Source `NLargest`: the constructor's `init_()` declares `n` (default 200,
/// no restriction) and `updateMembers_` reads it into `peakcount_`, which is
/// [`NLargest::n`] here. `NLargest(UInt n)` is `NLargest { n }`.
impl DefaultParamHandling for NLargest {
    const NAME: &'static str = "NLargest";

    /// `n` = 200, "The number of peaks to keep".
    fn defaults() -> Result<Param> {
        let mut defaults = Param::new();
        defaults.set_value(
            "n",
            ParamValue::Integer(200),
            "The number of peaks to keep",
            &[],
        )?;
        Ok(defaults)
    }

    /// # Errors
    ///
    /// [`Error::InvalidValue`] when [`NLargest::n`] exceeds the source's
    /// 32-bit `int` parameter.
    fn to_param(&self) -> Result<Param> {
        let mut param = Self::defaults()?;
        param.replace_value("n", int_parameter(self.n, "n")?)?;
        Ok(param)
    }

    /// `peakcount_ = (UInt)param_.getValue("n")`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a negative `n`, where the source's
    /// conversion throws `Exception::ConversionError` (the parameter has no
    /// minimum, so `setParameters` accepts it and `updateMembers_` throws).
    fn update_members(param: &Param) -> Result<Self> {
        Ok(Self {
            n: unsigned_parameter(param, "n")? as usize,
        })
    }
}

/// Source `ThresholdMower`: the constructor declares `threshold` (default
/// 0.05, no restriction), which the source re-reads from `param_` on every
/// `filterSpectrum` call and never caches, so it has no `updateMembers_`. The
/// typed [`ThresholdMower::threshold`] is that value: changing the field is
/// what a source `setParameters` between two calls does.
impl DefaultParamHandling for ThresholdMower {
    const NAME: &'static str = "ThresholdMower";

    /// `threshold` = 0.05, "Intensity threshold, peaks below this threshold are
    /// discarded".
    fn defaults() -> Result<Param> {
        let mut defaults = Param::new();
        defaults.set_value(
            "threshold",
            ParamValue::Float(0.05),
            "Intensity threshold, peaks below this threshold are discarded",
            &[],
        )?;
        Ok(defaults)
    }

    fn to_param(&self) -> Result<Param> {
        let mut param = Self::defaults()?;
        param.replace_value("threshold", ParamValue::Float(self.threshold))?;
        Ok(param)
    }

    /// `(double)param_.getValue("threshold")`, as `filterSpectrum` reads it.
    /// A non-finite threshold is accepted here, as by the source, and refused by
    /// the filter.
    fn update_members(param: &Param) -> Result<Self> {
        Ok(Self {
            threshold: param.value("threshold")?.to_f64()?,
        })
    }
}

/// Source `Normalizer`: the constructor declares `method` (default `to_one`,
/// valid strings `to_one` and `to_TIC`) and the public `updateMembers_` copies
/// it into `method_`, which is [`Normalizer::method`] here.
impl DefaultParamHandling for Normalizer {
    const NAME: &'static str = "Normalizer";

    /// `method` = `to_one`, restricted to `to_one` and `to_TIC`.
    fn defaults() -> Result<Param> {
        let mut defaults = Param::new();
        defaults.set_value(
            "method",
            ParamValue::String("to_one".into()),
            "Normalize via dividing by TIC ('to_TIC') per spectrum (i.e. all peaks sum to 1) or normalize to max. intensity to one ('to_one') per spectrum.",
            &[],
        )?;
        defaults.set_valid_strings("method", &["to_one".into(), "to_TIC".into()])?;
        Ok(defaults)
    }

    fn to_param(&self) -> Result<Param> {
        let mut param = Self::defaults()?;
        let method = match self.method {
            NormalizationMethod::ToOne => "to_one",
            NormalizationMethod::ToTic => "to_TIC",
        };
        param.replace_value("method", ParamValue::String(method.into()))?;
        Ok(param)
    }

    /// `method_ = param_.getValue("method").toString()`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a method other than the two valid strings,
    /// which the restriction check refuses before this runs. The source stores
    /// any string and throws `Exception::InvalidValue` ("Method not known")
    /// from `filterSpectrum` instead; the enum cannot hold another method.
    fn update_members(param: &Param) -> Result<Self> {
        let method = match param.value("method")?.as_str()? {
            "to_one" => NormalizationMethod::ToOne,
            "to_TIC" => NormalizationMethod::ToTic,
            other => {
                return Err(Error::InvalidValue(format!(
                    "Normalizer: method '{other}' is not known"
                )));
            }
        };
        Ok(Self { method })
    }
}

/// Source `RankScaler`: a `DefaultParamHandler` named `RankScaler` that
/// declares no parameters and never calls `defaultsToParam_()`. Its
/// `setParameters` therefore always warns "No default parameters for
/// DefaultParameterHandler 'RankScaler' specified!", as the Release build
/// does.
impl DefaultParamHandling for RankScaler {
    const NAME: &'static str = "RankScaler";

    /// Empty.
    fn defaults() -> Result<Param> {
        Ok(Param::new())
    }

    /// Empty.
    fn to_param(&self) -> Result<Param> {
        Ok(Param::new())
    }

    fn update_members(_param: &Param) -> Result<Self> {
        Ok(Self)
    }
}

/// Source `SqrtScaler`: a `DefaultParamHandler` named `SqrtScaler` that
/// declares no parameters; see the [`RankScaler`] implementation for the
/// resulting warning.
impl DefaultParamHandling for SqrtScaler {
    const NAME: &'static str = "SqrtScaler";

    /// Empty.
    fn defaults() -> Result<Param> {
        Ok(Param::new())
    }

    /// Empty.
    fn to_param(&self) -> Result<Param> {
        Ok(Param::new())
    }

    fn update_members(_param: &Param) -> Result<Self> {
        Ok(Self)
    }
}

/// Source `LinearResamplerAlign`: the constructor declares `spacing` (default
/// 0.05) and `ppm` (the string `false`, with no valid-strings restriction), and
/// `updateMembers_` reads them into `spacing_` and `ppm_`.
impl DefaultParamHandling for LinearResamplerAlign {
    const NAME: &'static str = "LinearResamplerAlign";

    /// `spacing` = 0.05 and `ppm` = `false`, both unrestricted.
    fn defaults() -> Result<Param> {
        let mut defaults = Param::new();
        defaults.set_value(
            "spacing",
            ParamValue::Float(0.05),
            "Spacing of the resampled output peaks.",
            &[],
        )?;
        defaults.set_value(
            "ppm",
            ParamValue::String("false".into()),
            "Whether spacing is in ppm or Th",
            &[],
        )?;
        Ok(defaults)
    }

    /// `ppm` is always `false`: this port resamples on absolute grids only.
    fn to_param(&self) -> Result<Param> {
        let mut param = Self::defaults()?;
        param.replace_value("spacing", ParamValue::Float(self.spacing))?;
        Ok(param)
    }

    /// `spacing_ = param_.getValue("spacing")`, `ppm_ =
    /// param_.getValue("ppm").toBool()`.
    ///
    /// `spacing` is taken as given, as in the source; the resampling operations
    /// refuse a spacing that is not finite and positive. The native
    /// [`LinearResamplerAlign::max_points`] keeps its default.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `ppm` is neither `true` nor `false`, where
    /// the source's `toBool()` throws `Exception::ConversionError`; and
    /// [`Error::Unsupported`] for `ppm` = `true`, because the ppm branch of the
    /// source's `populate_raster_` is not ported.
    fn update_members(param: &Param) -> Result<Self> {
        let spacing = param.value("spacing")?.to_f64()?;
        if param.value("ppm")?.to_bool()? {
            return Err(Error::Unsupported(
                "LinearResamplerAlign: ppm spacing is not ported".into(),
            ));
        }
        Ok(Self {
            spacing,
            ..Self::default()
        })
    }

    /// Keeps [`LinearResamplerAlign::max_points`].
    fn keep_native_members(&self, updated: Self) -> Self {
        Self {
            max_points: self.max_points,
            ..updated
        }
    }
}
