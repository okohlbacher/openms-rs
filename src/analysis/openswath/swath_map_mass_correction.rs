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
//! [`SwathMapMassCorrection::defaults`](crate::analysis::openswath::swath_map_mass_correction::SwathMapMassCorrection::defaults), [`SwathMapMassCorrection::from_param`](crate::analysis::openswath::swath_map_mass_correction::SwathMapMassCorrection::from_param)
//! and [`SwathMapMassCorrection::to_param`](crate::analysis::openswath::swath_map_mass_correction::SwathMapMassCorrection::to_param)), the estimated-window accessors,
//! [`SwathMapMassCorrection::estimate_window`](crate::analysis::openswath::swath_map_mass_correction::SwathMapMassCorrection::estimate_window),
//! [`SwathMapMassCorrection::find_swath_maps_pasef`](crate::analysis::openswath::swath_map_mass_correction::SwathMapMassCorrection::find_swath_maps_pasef),
//! and the two corrections,
//! [`SwathMapMassCorrection::correct_mz`](crate::analysis::openswath::swath_map_mass_correction::SwathMapMassCorrection::correct_mz) and
//! [`SwathMapMassCorrection::correct_im`](crate::analysis::openswath::swath_map_mass_correction::SwathMapMassCorrection::correct_im), with the file-local
//! `findBestFeature` and `findSwathMaps`.
//!
//! The source runs `correctIM` under `#pragma omp parallel for`; the port is
//! serial, which is what the source does with one thread (the order in which
//! the calibration points reach the regression, and so its last bits, depend on
//! the thread schedule in the source).

use super::data_structures::SpectrumPtr;
use super::dia_helper::{create_mz_range_ppm, integrate_window};
use super::openswath_scoring::OpenSwathScoring;
use super::spectrum_access::SpectrumAccessPtr;
use super::spectrum_access_openms::SpectrumAccessQuadMZTransforming;
use super::swath_map::SwathMap;
use super::transition_experiment::{LightTargetedExperiment, LightTransition};
use crate::analysis::feature_finder_picked::glibc_libm;
use crate::analysis::transformations::{
    DataPoint, LinearOptions, ModelConfig, TransformationDescription, WindowOptions,
};
use crate::format::file_info::text_format::ostream_g;
use crate::kernel::MSChromatogram;
use crate::kernel::mrm::MRMTransitionGroup;
use crate::kernel::ranges::RangeBase;
use crate::math::ml::linear_regression::LinearRegression;
use crate::math::ml::quadratic_regression::QuadraticRegression;
use crate::math::statistic_functions::{
    DEFAULT_R_DENSE, DEFAULT_R_SPARSE, DEFAULT_TUKEY_FACTOR, adaptive_quantile,
};
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::io::Write;
use std::sync::Arc;

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

    /// Estimate and apply an m/z calibration from the best feature of every
    /// transition group (source `correctMZ`).
    ///
    /// For each group with a feature, the spectrum at the retention time of its
    /// best feature (highest overall quality, the first on ties) is fetched
    /// from the SWATH map or maps whose window holds the first transition's
    /// precursor, and every transition's fragment m/z is integrated; with MS1
    /// maps the precursor m/z is integrated too. From the fragment residuals in
    /// ppm [`fragment_mz_window`](Self::fragment_mz_window) is estimated, from
    /// the MS1 residuals [`precursor_mz_window`](Self::precursor_mz_window).
    /// Unless `mz_correction_function` is `none` or fewer than three points
    /// were found, a regression of that function is fitted and **every** map of
    /// `swath_maps` (MS1 maps included) is wrapped in a
    /// [`SpectrumAccessQuadMZTransforming`] carrying it, the `*_delta_ppm`
    /// functions in the ppm domain.
    ///
    /// `pasef` selects the window by precursor m/z and ion mobility (see
    /// [`find_swath_maps_pasef`](Self::find_swath_maps_pasef)); with
    /// `im_extraction_window` set, the fetch and the integration are
    /// restricted to the peptide's library ion mobility.
    ///
    /// With a `debug_mz_file` set, every calibration point is written to it as
    /// `mz theo_mz drift_time diff_ppm log_intensity RT`, tab separated, with
    /// 15 significant digits as the source's `writtenDigits<double>()` gives.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a group with features but no transition (the
    /// source reads `getTransitions()[0]` unchecked), for an unknown
    /// `mz_correction_function` once at least three points were found (the
    /// source's `Exception::IllegalArgument`, same message), and the errors of
    /// the spectrum access, the integration and the regression (a degenerate
    /// fit is an error, as the source's `UnableToFit`). [`Error::MissingInformation`]
    /// for a map without spectrum access when the maps are replaced. A debug
    /// file that cannot be created is an error here where the source silently
    /// writes nothing. On error `swath_maps` is unchanged.
    pub fn correct_mz(
        &mut self,
        transition_group_map: &BTreeMap<String, &MRMTransitionGroupType>,
        targeted_exp: &LightTargetedExperiment,
        swath_maps: &mut [SwathMap],
        pasef: bool,
    ) -> Result<()> {
        let ppm = self.mz_extraction_window_ppm;
        let mz_extr_window = self.mz_extraction_window;
        let corr_type = self.mz_correction_function.clone();
        let im_extraction = self.im_extraction_window;
        let padding_factor = self.mz_estimation_padding_factor;
        let is_ppm = matches!(
            corr_type.as_str(),
            "quadratic_regression_delta_ppm"
                | "weighted_quadratic_regression_delta_ppm"
                | "regression_delta_ppm"
        );
        if corr_type == "none" {
            return Ok(());
        }

        let mut debug = open_debug(
            &self.debug_mz_file,
            "mz\ttheo_mz\tdrift_time\tdiff_ppm\tlog_intensity\tRT\n",
        )?;

        let mut weights: Vec<f64> = Vec::new();
        let mut exp_mz: Vec<f64> = Vec::new();
        let mut theo_mz: Vec<f64> = Vec::new();
        let mut delta_ppm: Vec<f64> = Vec::new();
        let mut delta_ppm_ms1: Vec<f64> = Vec::new();

        let pep_im_map = pep_im_map(targeted_exp);
        let ms1_maps: Vec<SwathMap> = swath_maps.iter().filter(|m| m.ms1).cloned().collect();
        let scoring = OpenSwathScoring::default();

        for group in transition_group_map.values() {
            if group.features().is_empty() {
                continue;
            }
            let first = group.transitions().first().ok_or_else(|| {
                Error::InvalidValue("a transition group with features has no transitions".into())
            })?;
            let drift_target = pep_im_map
                .get(first.peptide_ref.as_str())
                .copied()
                .unwrap_or(0.0);

            let best_rt = find_best_feature(group);
            let used_maps = if pasef {
                self.find_swath_maps_pasef(group, swath_maps)?
            } else {
                find_swath_maps(group, swath_maps)
            };
            if used_maps.is_empty() {
                continue;
            }

            // an ion mobility window when one is configured
            let mut im_range = RangeBase::new();
            if im_extraction != -1.0 {
                im_range.set_max(drift_target)?;
                im_range.min_span_if_singular(im_extraction)?;
            }

            let spectra = scoring.fetch_spectrum_swath_maps(&used_maps, best_rt, 1, &im_range)?;
            let Some(sp) = spectra.first() else {
                continue;
            };
            for tr in group.transitions() {
                let mz_range = create_mz_range_ppm(tr.product_mz, mz_extr_window, ppm)?;
                // integrate at the position of the theoretical mass
                let Some(found) = integrate_window(sp, &mz_range, &im_range, false)? else {
                    continue;
                };
                let (mz, intensity) = (found.mz, found.intensity);
                let weight = glibc_libm::log(intensity) / glibc_libm::log(2.0);
                weights.push(weight);
                exp_mz.push(mz);
                theo_mz.push(tr.product_mz);
                let diff_ppm = (mz - tr.product_mz) * 1_000_000.0 / mz;
                delta_ppm.push(diff_ppm);
                if let Some(out) = debug.as_mut() {
                    writeln!(
                        out,
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        ostream_g(mz, 15),
                        ostream_g(tr.product_mz, 15),
                        ostream_g(drift_target, 15),
                        ostream_g(diff_ppm, 15),
                        ostream_g(weight, 15),
                        ostream_g(best_rt, 15)
                    )?;
                }
            }

            // MS1 precursor residuals for the precursor window
            if !ms1_maps.is_empty() {
                let spectra =
                    scoring.fetch_spectrum_swath_maps(&ms1_maps, best_rt, 1, &im_range)?;
                if let Some(sp_ms1) = spectra.first() {
                    let theo_prec_mz = first.precursor_mz;
                    let mz_range = create_mz_range_ppm(theo_prec_mz, mz_extr_window, ppm)?;
                    if let Some(found) = integrate_window(sp_ms1, &mz_range, &im_range, false)? {
                        let dppm = (found.mz - theo_prec_mz) / theo_prec_mz * 1e6;
                        delta_ppm_ms1.push(dppm.abs());
                    }
                }
            }
        }

        // the fragment window, then the precursor window
        let fragment_window = Self::estimate_window(
            &delta_ppm,
            self.mz_estimation_percentile / 100.0,
            true,
            padding_factor,
        )?;
        self.set_fragment_mz_window(fragment_window);
        if !delta_ppm_ms1.is_empty() {
            delta_ppm_ms1.sort_by(f64::total_cmp);
            let precursor_window = Self::estimate_window(
                &delta_ppm_ms1,
                self.mz_estimation_percentile / 100.0,
                true,
                padding_factor,
            )?;
            self.set_precursor_mz_window(precursor_window);
        }

        if exp_mz.len() < 3 {
            return Ok(());
        }
        let params: [f64; 3] = match corr_type.as_str() {
            "unweighted_regression" => {
                let mut lr = LinearRegression::new();
                lr.compute_regression(0.0, &exp_mz, &theo_mz, true)?;
                [lr.intercept(), lr.slope(), 0.0]
            }
            "weighted_regression" => {
                let mut lr = LinearRegression::new();
                lr.compute_regression_weighted(0.0, &exp_mz, &theo_mz, &weights, true)?;
                [lr.intercept(), lr.slope(), 0.0]
            }
            "quadratic_regression" => {
                let mut qr = QuadraticRegression::new();
                qr.compute_regression(&exp_mz, &theo_mz)?;
                [qr.a(), qr.b(), qr.c()]
            }
            "weighted_quadratic_regression" => {
                let mut qr = QuadraticRegression::new();
                qr.compute_regression_weighted(&exp_mz, &theo_mz, &weights)?;
                [qr.a(), qr.b(), qr.c()]
            }
            "quadratic_regression_delta_ppm" => {
                let mut qr = QuadraticRegression::new();
                qr.compute_regression(&exp_mz, &delta_ppm)?;
                [qr.a(), qr.b(), qr.c()]
            }
            "regression_delta_ppm" => {
                let mut lr = LinearRegression::new();
                lr.compute_regression(0.0, &exp_mz, &delta_ppm, true)?;
                [lr.intercept(), lr.slope(), 0.0]
            }
            "weighted_quadratic_regression_delta_ppm" => {
                let mut qr = QuadraticRegression::new();
                qr.compute_regression_weighted(&exp_mz, &delta_ppm, &weights)?;
                [qr.a(), qr.b(), qr.c()]
            }
            other => {
                return Err(Error::InvalidValue(format!(
                    "Unknown correction type {other}"
                )));
            }
        };

        // Replace the SWATH maps' access with the transforming wrapper, all or
        // none.
        let mut wrapped: Vec<SpectrumAccessPtr> = Vec::with_capacity(swath_maps.len());
        for map in swath_maps.iter() {
            let sptr = map.sptr.clone().ok_or_else(|| {
                Error::MissingInformation("SWATH map without spectrum access".into())
            })?;
            wrapped.push(Arc::new(SpectrumAccessQuadMZTransforming::quad_mz(
                sptr, params[0], params[1], params[2], is_ppm,
            )));
        }
        for (map, access) in swath_maps.iter_mut().zip(wrapped) {
            map.sptr = Some(access);
        }
        Ok(())
    }

    /// Estimate and store the ion-mobility calibration from the best feature
    /// of every transition group (source `correctIM`).
    ///
    /// Nothing happens (and `im_trafo` is left as it is) when
    /// `im_extraction_window` is negative or `im_correction_function` is
    /// `none`. Otherwise, for each group with a feature, the measured ion
    /// mobility of every fragment (or, with `ms1_im_calibration`, of the first
    /// transition's precursor in the MS1 map) is integrated at the best
    /// feature's retention time inside the peptide's library drift time
    /// widened to `im_extraction_window`. The pairs (measured, library) become
    /// the data points of `im_trafo`, which is fitted with the linear model;
    /// [`fragment_im_window`](Self::fragment_im_window) is the residual
    /// window of that fit. When MS1 maps exist their precursor mobility points
    /// are collected whichever mode is used, and
    /// [`precursor_im_window`](Self::precursor_im_window) is estimated from
    /// their raw deviations (the source assigns them to a copy whose model the
    /// assignment resets).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] starting `UnableToFit` when no calibration point
    /// was found, for a group with features but no transition, and the errors of
    /// the access, the integration, the regression and the model fit.
    pub fn correct_im(
        &mut self,
        transition_group_map: &BTreeMap<String, &MRMTransitionGroupType>,
        targeted_exp: &LightTargetedExperiment,
        swath_maps: &[SwathMap],
        pasef: bool,
        im_trafo: &mut TransformationDescription,
    ) -> Result<()> {
        let ppm = self.mz_extraction_window_ppm;
        let mz_extr_window = self.mz_extraction_window;
        let im_extraction_win = self.im_extraction_window;
        let im_estimation_padding_factor = self.im_estimation_padding_factor;
        let ms1_im = self.ms1_im_calibration;

        if im_extraction_win < 0.0 {
            return Ok(());
        }
        if self.im_correction_function == "none" {
            return Ok(());
        }
        // anything else is the linear function

        let mut debug = open_debug(&self.debug_im_file, "mz\tim\ttheo_im\tRT\tintensity\n")?;

        let pep_im_map = pep_im_map(targeted_exp);
        let ms1_maps: Vec<SwathMap> = swath_maps.iter().filter(|m| m.ms1).cloned().collect();
        let scoring = OpenSwathScoring::default();

        let mut data_im: Vec<DataPoint> = Vec::new();
        let mut exp_im: Vec<f64> = Vec::new();
        let mut theo_im: Vec<f64> = Vec::new();
        let mut exp_im_ms1_all: Vec<f64> = Vec::new();
        let mut theo_im_ms1_all: Vec<f64> = Vec::new();

        let drift_of = |peptide_ref: &str| pep_im_map.get(peptide_ref).copied().unwrap_or(0.0);

        'groups: for group in transition_group_map.values() {
            if group.features().is_empty() {
                continue;
            }
            let transitions = group.transitions();
            if transitions.is_empty() {
                return Err(Error::InvalidValue(
                    "a transition group with features has no transitions".into(),
                ));
            }
            let best_rt = find_best_feature(group);
            let used_maps = if pasef {
                self.find_swath_maps_pasef(group, swath_maps)?
            } else {
                find_swath_maps(group, swath_maps)
            };
            if used_maps.is_empty() {
                continue;
            }

            let mut exp_im_ms1_local: Vec<f64> = Vec::new();
            let mut theo_im_ms1_local: Vec<f64> = Vec::new();

            // the spectrum at the best feature's retention time; the source
            // fetches with an empty ion-mobility range here
            let empty_im = RangeBase::new();
            let fetch = |maps: &[SwathMap], im: &RangeBase| -> Result<Option<SpectrumPtr>> {
                if maps.is_empty() {
                    // Release build: no maps, no spectrum
                    return Ok(None);
                }
                Ok(scoring
                    .fetch_spectrum_swath_maps(maps, best_rt, 1, im)?
                    .into_iter()
                    .next())
            };
            let sp_ms1 = if ms1_im {
                fetch(&ms1_maps, &empty_im)?
            } else {
                None
            };
            let sp_ms2 = if ms1_im {
                None
            } else {
                fetch(&used_maps, &empty_im)?
            };

            for tr in transitions {
                if ms1_im {
                    continue;
                }
                let mz_range = create_mz_range_ppm(tr.product_mz, mz_extr_window, ppm)?;
                let drift_target = drift_of(&tr.peptide_ref);
                let mut im_range = RangeBase::new();
                if im_extraction_win != -1.0 {
                    im_range.set_max(drift_target)?;
                    im_range.min_span_if_singular(im_extraction_win)?;
                }
                // the spectrum must have a drift time array
                let Some(sp) = sp_ms2.as_ref().filter(|s| s.drift_time_array().is_some()) else {
                    continue;
                };
                let Some(found) = integrate_window(sp, &mz_range, &im_range, false)? else {
                    continue;
                };
                // skip empty windows
                if found.im <= 0.0 {
                    continue;
                }
                data_im.push(DataPoint::new(found.im, drift_target));
                exp_im.push(found.im);
                theo_im.push(drift_target);
                if let Some(out) = debug.as_mut() {
                    writeln!(
                        out,
                        "{}\t{}\t{}\t{}\t{}",
                        ostream_g(tr.precursor_mz, 15),
                        ostream_g(found.im, 15),
                        ostream_g(drift_target, 15),
                        ostream_g(best_rt, 15),
                        ostream_g(found.intensity, 15)
                    )?;
                }
            }

            // MS1 mobility points for the window estimate, whatever the mode
            if !ms1_maps.is_empty() {
                let tr0 = &transitions[0];
                let drift_target0 = drift_of(&tr0.peptide_ref);
                let mut im_range_ms1 = RangeBase::singular(drift_target0)?;
                if im_extraction_win != -1.0 {
                    im_range_ms1.min_span_if_singular(im_extraction_win)?;
                }
                let collected = fetch(&ms1_maps, &im_range_ms1)?;
                if let Some(sp) = collected
                    .as_ref()
                    .filter(|s| s.drift_time_array().is_some())
                {
                    let mz_range = create_mz_range_ppm(tr0.precursor_mz, mz_extr_window, ppm)?;
                    if let Some(found) = integrate_window(sp, &mz_range, &im_range_ms1, false)? {
                        if found.im > 0.0 {
                            exp_im_ms1_local.push(found.im);
                            theo_im_ms1_local.push(drift_target0);
                        }
                    }
                }
            }

            // MS1 extraction for the calibration itself
            if ms1_im {
                let tr = &transitions[0];
                let mz_range = create_mz_range_ppm(tr.precursor_mz, mz_extr_window, ppm)?;
                let drift_target = drift_of(&tr.peptide_ref);
                // no check of the window needed: the correction is of the IM
                let mut im_range = RangeBase::singular(drift_target)?;
                im_range.min_span_if_singular(im_extraction_win)?;
                let Some(sp) = sp_ms1.as_ref().filter(|s| s.drift_time_array().is_some()) else {
                    // the source continues with the next group here, before the
                    // MS1 points of this group are added
                    continue 'groups;
                };
                let Some(found) = integrate_window(sp, &mz_range, &im_range, false)? else {
                    continue 'groups;
                };
                if found.im <= 0.0 {
                    continue 'groups;
                }
                data_im.push(DataPoint::new(found.im, drift_target));
                exp_im.push(found.im);
                theo_im.push(drift_target);
                if let Some(out) = debug.as_mut() {
                    writeln!(
                        out,
                        "{}\t{}\t{}\t{}\t{}",
                        ostream_g(tr.precursor_mz, 15),
                        ostream_g(found.im, 15),
                        ostream_g(drift_target, 15),
                        ostream_g(best_rt, 15),
                        ostream_g(found.intensity, 15)
                    )?;
                }
            }

            exp_im_ms1_all.extend(exp_im_ms1_local);
            theo_im_ms1_all.extend(theo_im_ms1_local);
        }

        drop(debug);

        if exp_im.is_empty() {
            return Err(Error::InvalidValue(
                "UnableToFit-LinearRegression: Could not fit a linear model to the data (0 points)."
                    .into(),
            ));
        }

        // the regression of the source is computed for its debug output only;
        // its failure is the source's failure
        let mut lr = LinearRegression::new();
        lr.compute_regression(0.0, &exp_im, &theo_im, true)?;

        // store the IM transformation with the linear model
        im_trafo.set_data_points(data_im)?;
        im_trafo.fit_model(ModelConfig::Linear(LinearOptions::default()))?;

        // the MS2 window from the residual quantile of the fit
        let fragment_im_window = im_trafo.estimate_window(WindowOptions {
            quantile: self.im_estimation_percentile / 100.0,
            inverse: false,
            full_window: true,
            padding_factor: im_estimation_padding_factor,
        })?;
        self.set_fragment_im_window(fragment_im_window);

        if !exp_im_ms1_all.is_empty() {
            let ms1_points: Vec<DataPoint> = exp_im_ms1_all
                .iter()
                .zip(&theo_im_ms1_all)
                .map(|(&x, &y)| DataPoint::new(x, y))
                .collect();
            // `setDataPoints` of the source resets the model of the copy, so the
            // window is that of the raw deviations (experimental to library)
            let mut inverse_copy = im_trafo.clone();
            inverse_copy.set_data_points(ms1_points)?;
            let precursor_im_window = inverse_copy.estimate_window(WindowOptions {
                quantile: self.im_estimation_percentile / 100.0,
                inverse: false,
                full_window: true,
                padding_factor: im_estimation_padding_factor,
            })?;
            self.set_precursor_im_window(precursor_im_window);
        }
        Ok(())
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

/// The feature with the highest overall quality of a transition group, as its
/// retention time (source file-local `findBestFeature`): -1 when there is none
/// above the initial threshold of -1000; the first feature wins ties.
fn find_best_feature(group: &MRMTransitionGroupType) -> f64 {
    let mut best_rt = -1.0;
    let mut highest_score = -1000.0;
    for feature in group.features() {
        let quality = f64::from(feature.feature.base.quality);
        if quality > highest_score {
            best_rt = feature.rt;
            highest_score = quality;
        }
    }
    best_rt
}

/// The SWATH maps whose window holds the first transition's precursor,
/// `lower < precursor <= upper` (source file-local `findSwathMaps`).
fn find_swath_maps(group: &MRMTransitionGroupType, swath_maps: &[SwathMap]) -> Vec<SwathMap> {
    let Some(first) = group.transitions().first() else {
        return Vec::new();
    };
    swath_maps
        .iter()
        .filter(|m| m.lower < first.precursor_mz && m.upper >= first.precursor_mz)
        .cloned()
        .collect()
}

/// Library drift time per compound id; a later compound with the same id
/// overwrites an earlier one, as the source's map assignment does.
fn pep_im_map(targeted_exp: &LightTargetedExperiment) -> BTreeMap<&str, f64> {
    let mut map = BTreeMap::new();
    for compound in targeted_exp.compounds() {
        map.insert(compound.id.as_str(), compound.drift_time);
    }
    map
}

/// Open a debug file and write its header; `None` for an empty file name.
fn open_debug(path: &str, header: &str) -> Result<Option<std::io::BufWriter<std::fs::File>>> {
    if path.is_empty() {
        return Ok(None);
    }
    let mut out = std::io::BufWriter::new(std::fs::File::create(path)?);
    out.write_all(header.as_bytes())?;
    Ok(Some(out))
}
