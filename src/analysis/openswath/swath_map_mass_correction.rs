// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! m/z and ion-mobility correction of SWATH maps from calibrant peptides.
//!
//! Ports `ANALYSIS/OPENSWATH/SwathMapMassCorrection.h` with
//! `SwathMapMassCorrection.cpp`, **partially**. The source uses a set of
//! pre-determined points (the RT-normalization peptides) in a SWATH-MS map to
//! correct all maps for the m/z (and ion-mobility) shift found at those
//! points. See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_EXTRACTION_SUPPORT.md).
//!
//! Ported: the parameters (`DefaultParamHandler` becomes
//! [`SwathMapMassCorrection::defaults`], [`SwathMapMassCorrection::from_param`]
//! and [`SwathMapMassCorrection::to_param`]), the estimated-window accessors,
//! [`SwathMapMassCorrection::estimate_window`] and
//! [`SwathMapMassCorrection::find_swath_maps_pasef`].
//!
//! **Not ported:** `correctMZ` and `correctIM`, which need
//! `OpenSwathScoring::fetchSpectrumSwath`, `DIAHelpers::createMZRangePPM`
//! and `integrateWindow`, and the `LinearRegression` / `QuadraticRegression`
//! fits, none of which is ported yet; with them go the file-local helpers
//! `findBestFeature` and `findSwathMaps`, which only those two call.

use super::swath_map::SwathMap;
use super::transition_experiment::LightTransition;
use crate::kernel::MSChromatogram;
use crate::kernel::mrm::MRMTransitionGroup;
use crate::math::statistic_functions::{
    DEFAULT_R_DENSE, DEFAULT_R_SPARSE, DEFAULT_TUKEY_FACTOR, adaptive_quantile,
};
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::{Error, Result};

/// Parameter handler name of the source class.
pub const SWATH_MAP_MASS_CORRECTION_NAME: &str = "SwathMapMassCorrection";

/// The valid values of `mz_correction_function`, in source order.
pub const MZ_CORRECTION_FUNCTIONS: [&str; 8] = [
    "none",
    "regression_delta_ppm",
    "unweighted_regression",
    "weighted_regression",
    "quadratic_regression",
    "weighted_quadratic_regression",
    "weighted_quadratic_regression_delta_ppm",
    "quadratic_regression_delta_ppm",
];

/// The valid values of `im_correction_function`, in source order.
pub const IM_CORRECTION_FUNCTIONS: [&str; 2] = ["none", "linear"];

/// The transition group type the correction works on, source
/// `MRMFeatureFinderScoring::MRMTransitionGroupType`.
pub type MRMTransitionGroupType = MRMTransitionGroup<MSChromatogram, LightTransition>;

/// Correction of SWATH maps in m/z and ion mobility.
///
/// Source `SwathMapMassCorrection`, a `DefaultParamHandler`. Each public field
/// is one parameter of [`SwathMapMassCorrection::defaults`], under the same
/// name (the source's `updateMembers_` members); the estimated windows are
/// state, read and written through their accessors.
#[derive(Clone, Debug, PartialEq)]
pub struct SwathMapMassCorrection {
    /// `mz_extraction_window`: m/z extraction window width. Default -1.
    pub mz_extraction_window: f64,
    /// `mz_extraction_window_ppm`: whether the m/z window is in ppm. Default false.
    pub mz_extraction_window_ppm: bool,
    /// `ms1_im_calibration`: calibrate ion mobility from MS1 precursors instead
    /// of MS2 fragments. Default false.
    pub ms1_im_calibration: bool,
    /// `im_extraction_window`: ion-mobility extraction window width. Default -1.
    pub im_extraction_window: f64,
    /// `mz_estimation_padding_factor`: factor on the estimated m/z window
    /// (1.3 adds 30%). At least 1; default 1.3.
    pub mz_estimation_padding_factor: f64,
    /// `im_estimation_padding_factor`: factor on the estimated ion-mobility
    /// window. At least 1; default 1.3.
    pub im_estimation_padding_factor: f64,
    /// `mz_estimation_percentile`: percentile for the m/z window estimate, in
    /// `[25, 99.9]`. Default 99.
    pub mz_estimation_percentile: f64,
    /// `im_estimation_percentile`: percentile for the ion-mobility window
    /// estimate, in `[25, 99.9]`. Default 99.
    pub im_estimation_percentile: f64,
    /// `mz_correction_function`: one of [`MZ_CORRECTION_FUNCTIONS`]. Default `none`.
    pub mz_correction_function: String,
    /// `im_correction_function`: one of [`IM_CORRECTION_FUNCTIONS`]. Default `linear`.
    pub im_correction_function: String,
    /// `debug_im_file`: debug output of the ion-mobility calibration. Default empty.
    pub debug_im_file: String,
    /// `debug_mz_file`: debug output of the m/z calibration. Default empty.
    pub debug_mz_file: String,
    fragment_mz_window: f64,
    fragment_im_window: f64,
    precursor_mz_window: f64,
    precursor_im_window: f64,
}

impl Default for SwathMapMassCorrection {
    fn default() -> Self {
        Self {
            mz_extraction_window: -1.0,
            mz_extraction_window_ppm: false,
            ms1_im_calibration: false,
            im_extraction_window: -1.0,
            mz_estimation_padding_factor: 1.3,
            im_estimation_padding_factor: 1.3,
            mz_estimation_percentile: 99.0,
            im_estimation_percentile: 99.0,
            mz_correction_function: "none".into(),
            im_correction_function: "linear".into(),
            debug_im_file: String::new(),
            debug_mz_file: String::new(),
            fragment_mz_window: -1.0,
            fragment_im_window: -1.0,
            precursor_mz_window: -1.0,
            precursor_im_window: -1.0,
        }
    }
}

fn entry(name: &str, value: ParamValue, description: &str, advanced: bool) -> ParamEntry {
    let mut e = ParamEntry {
        name: name.into(),
        description: description.into(),
        value,
        ..ParamEntry::default()
    };
    if advanced {
        e.tags.insert("advanced".into());
    }
    e
}

fn float_entry(
    name: &str,
    value: f64,
    description: &str,
    min: Option<f64>,
    max: Option<f64>,
) -> ParamEntry {
    let mut e = entry(name, ParamValue::Float(value), description, false);
    if let Some(min) = min {
        e.min_float = min;
    }
    if let Some(max) = max {
        e.max_float = max;
    }
    e
}

fn string_entry(
    name: &str,
    value: &str,
    description: &str,
    advanced: bool,
    valid: &[&str],
) -> ParamEntry {
    let mut e = entry(
        name,
        ParamValue::String(value.into()),
        description,
        advanced,
    );
    e.valid_strings = valid.iter().map(|v| (*v).to_owned()).collect();
    e
}

fn float(param: &Param, key: &str) -> Result<f64> {
    match param.value(key)? {
        ParamValue::Float(value) => Ok(*value),
        _ => Err(Error::InvalidValue(format!(
            "parameter '{key}' must be a floating-point value"
        ))),
    }
}

fn text(param: &Param, key: &str) -> Result<String> {
    Ok(param.value(key)?.as_str()?.to_owned())
}

impl SwathMapMassCorrection {
    /// A corrector with the default parameters. Source constructor.
    pub fn new() -> Self {
        Self::default()
    }

    /// The default parameters, source `getDefaults` after the constructor.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure, which the constant tree cannot
    /// trigger.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters in the layout of [`SwathMapMassCorrection::defaults`],
    /// carrying the current values (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        let flag = |b: bool| if b { "true" } else { "false" };
        let entries = vec![
            float_entry(
                "mz_extraction_window",
                self.mz_extraction_window,
                "M/z extraction window width",
                None,
                None,
            ),
            string_entry(
                "mz_extraction_window_ppm",
                flag(self.mz_extraction_window_ppm),
                "Whether m/z extraction is in ppm",
                true,
                &["true", "false"],
            ),
            string_entry(
                "ms1_im_calibration",
                flag(self.ms1_im_calibration),
                "Whether to use MS1 precursor data for the ion mobility calibration (default = false, uses MS2 / fragment ions for calibration)",
                true,
                &["true", "false"],
            ),
            float_entry(
                "im_extraction_window",
                self.im_extraction_window,
                "Ion mobility extraction window width",
                None,
                None,
            ),
            float_entry(
                "mz_estimation_padding_factor",
                self.mz_estimation_padding_factor,
                "A padding factor to multiply the estimated m/z window by. For example, a factor of 1.3 will add a 30% padding to the estimated m/z window, so if the estimated m/z window is 18, then 5.4 will be added for a total estimated m/z window of 23.4. A factor of 1.0 will not add any padding to the estimated window.",
                Some(1.0),
                None,
            ),
            float_entry(
                "im_estimation_padding_factor",
                self.im_estimation_padding_factor,
                "A padding factor to multiply the estimated ion_mobility window by. For example, a factor of 1.3 will add a 30% padding to the estimated ion_mobility window, so if the estimated ion_mobility window is 0.03, then 0.009 will be added for a total estimated ion_mobility window of 0.039. A factor of 1.0 will not add any padding to the estimated window.",
                Some(1.0),
                None,
            ),
            float_entry(
                "mz_estimation_percentile",
                self.mz_estimation_percentile,
                "Percentile for m/z window estimation (25.0-99.9)",
                Some(25.0),
                Some(99.9),
            ),
            float_entry(
                "im_estimation_percentile",
                self.im_estimation_percentile,
                "Percentile for ion mobility window estimation (25.0-99.9)",
                Some(25.0),
                Some(99.9),
            ),
            string_entry(
                "mz_correction_function",
                &self.mz_correction_function,
                "Type of normalization function for m/z calibration.",
                false,
                &MZ_CORRECTION_FUNCTIONS,
            ),
            string_entry(
                "im_correction_function",
                &self.im_correction_function,
                "Type of normalization function for IM calibration.",
                false,
                &IM_CORRECTION_FUNCTIONS,
            ),
            entry(
                "debug_im_file",
                ParamValue::String(self.debug_im_file.clone()),
                "Debug file for Ion Mobility calibration.",
                false,
            ),
            entry(
                "debug_mz_file",
                ParamValue::String(self.debug_mz_file.clone()),
                "Debug file for m/z calibration.",
                false,
            ),
        ];
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries,
            nodes: Vec::new(),
        })
    }

    /// Build a corrector from parameters, as source `setParameters` followed
    /// by `updateMembers_`.
    ///
    /// Missing parameters take their defaults; unknown ones are ignored, as
    /// the source only warns (see
    /// [`SwathMapMassCorrection::from_param_with_warnings`]). The estimated
    /// windows start at -1, as after the source constructor (they are not
    /// parameters).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a parameter has the wrong value type or
    /// violates its restriction (the source `Exception::InvalidParameter`).
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`SwathMapMassCorrection::from_param`], also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`SwathMapMassCorrection::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(SWATH_MAP_MASS_CORRECTION_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, Self::from_complete_param)
    }

    fn from_complete_param(param: &Param) -> Result<Self> {
        Ok(Self {
            mz_extraction_window: float(param, "mz_extraction_window")?,
            // The source compares `== "true"`.
            mz_extraction_window_ppm: param.value("mz_extraction_window_ppm")?.as_str()? == "true",
            ms1_im_calibration: param.value("ms1_im_calibration")?.as_str()? == "true",
            im_extraction_window: float(param, "im_extraction_window")?,
            mz_estimation_padding_factor: float(param, "mz_estimation_padding_factor")?,
            im_estimation_padding_factor: float(param, "im_estimation_padding_factor")?,
            mz_estimation_percentile: float(param, "mz_estimation_percentile")?,
            im_estimation_percentile: float(param, "im_estimation_percentile")?,
            mz_correction_function: text(param, "mz_correction_function")?,
            im_correction_function: text(param, "im_correction_function")?,
            debug_im_file: text(param, "debug_im_file")?,
            debug_mz_file: text(param, "debug_mz_file")?,
            ..Self::default()
        })
    }

    /// Estimate an extraction window from absolute residuals.
    ///
    /// Source static `estimateWindow`. The residuals (for example `|delta
    /// ppm|` for m/z) are made absolute, then an adaptive half-width is taken
    /// with [`adaptive_quantile`] at `quantile` (Tukey fence `k = 1.5`; the
    /// robust, winsorized quantile is favoured while at most 1% of the
    /// residuals lie above the fence and the raw quantile once 10% or more
    /// do). The result is the half-width, or the full width (twice the
    /// half-width) when `full_width`, multiplied by `padding_factor`, in the
    /// units of the residuals. An empty input yields 0. The source's defaults
    /// are `quantile = 0.99`, `full_width = true`, `padding_factor = 1.0`.
    ///
    /// # Errors
    ///
    /// The errors of [`adaptive_quantile`]: [`Error::InvalidValue`] for a
    /// `quantile` outside `[0, 1]`, [`Error::InvalidRange`] for an input over
    /// its ceilings.
    pub fn estimate_window(
        residuals: &[f64],
        quantile: f64,
        full_width: bool,
        padding_factor: f64,
    ) -> Result<f64> {
        if residuals.is_empty() {
            return Ok(0.0);
        }
        let absolute: Vec<f64> = residuals.iter().map(|d| d.abs()).collect();
        let res = adaptive_quantile(
            &absolute,
            quantile,
            DEFAULT_TUKEY_FACTOR,
            DEFAULT_R_SPARSE,
            DEFAULT_R_DENSE,
        )?;
        let half = if full_width {
            2.0 * res.blended
        } else {
            res.blended
        };
        Ok(half * padding_factor)
    }

    /// The SWATH map of a PASEF transition group, where windows may share an
    /// m/z range and differ in ion mobility.
    ///
    /// Source `findSwathMapsPasef`. A map isolates the group when
    /// `lower < precursor m/z <= upper` and `im_lower < precursor IM <=
    /// im_upper` of the group's first transition. Of several such maps the
    /// one whose ion-mobility centre is closest to the precursor's is chosen
    /// (the first on a tie). The result holds that one map, or is empty; it
    /// is a vector for compatibility with callers (SONAR can use several).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the group has no transition, or its first
    /// transition has no ion mobility (-1); the source states the latter as a
    /// precondition checked only in debug builds, and reads out of bounds for
    /// the former.
    pub fn find_swath_maps_pasef(
        &self,
        transition_group: &MRMTransitionGroupType,
        swath_maps: &[SwathMap],
    ) -> Result<Vec<SwathMap>> {
        let first = transition_group
            .transitions()
            .first()
            .ok_or_else(|| Error::InvalidValue("the transition group has no transitions".into()))?;
        let (precursor_mz, precursor_im) = (first.precursor_mz, first.precursor_im);
        if precursor_im == -1.0 {
            return Err(Error::InvalidValue(
                "All transitions must have a valid IM value (not -1)".into(),
            ));
        }
        let mut used: Option<&SwathMap> = None;
        for m in swath_maps {
            if m.lower < precursor_mz
                && m.upper >= precursor_mz
                && m.im_lower < precursor_im
                && m.im_upper >= precursor_im
            {
                match used {
                    None => used = Some(m),
                    Some(old) => {
                        let old_diff = ((old.im_lower + old.im_upper) / 2.0 - precursor_im).abs();
                        let new_diff = ((m.im_lower + m.im_upper) / 2.0 - precursor_im).abs();
                        if old_diff > new_diff {
                            used = Some(m);
                        }
                    }
                }
            }
        }
        Ok(used.into_iter().cloned().collect())
    }

    /// The estimated fragment m/z extraction window (ppm); -1 until set.
    /// Source `getFragmentMzWindow`.
    pub fn fragment_mz_window(&self) -> f64 {
        self.fragment_mz_window
    }

    /// Set the estimated fragment m/z extraction window (ppm). Source
    /// `setFragmentMzWindow`.
    pub fn set_fragment_mz_window(&mut self, window: f64) {
        self.fragment_mz_window = window;
    }

    /// The estimated fragment ion-mobility extraction window; -1 until set.
    /// Source `getFragmentImWindow`.
    pub fn fragment_im_window(&self) -> f64 {
        self.fragment_im_window
    }

    /// Set the estimated fragment ion-mobility extraction window. Source
    /// `setFragmentImWindow`.
    pub fn set_fragment_im_window(&mut self, window: f64) {
        self.fragment_im_window = window;
    }

    /// The estimated precursor m/z extraction window (ppm); -1 until set.
    /// Source `getPrecursorMzWindow`.
    pub fn precursor_mz_window(&self) -> f64 {
        self.precursor_mz_window
    }

    /// Set the estimated precursor m/z extraction window (ppm). Source
    /// `setPrecursorMzWindow`.
    pub fn set_precursor_mz_window(&mut self, window: f64) {
        self.precursor_mz_window = window;
    }

    /// The estimated precursor ion-mobility extraction window; -1 until set.
    /// Source `getPrecursorImWindow`.
    pub fn precursor_im_window(&self) -> f64 {
        self.precursor_im_window
    }

    /// Set the estimated precursor ion-mobility extraction window. Source
    /// `setPrecursorImWindow`.
    pub fn set_precursor_im_window(&mut self, window: f64) {
        self.precursor_im_window = window;
    }
}
