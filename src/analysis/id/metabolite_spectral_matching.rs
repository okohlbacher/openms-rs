// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Spectral-library search for metabolite MS2 spectra.
//!
//! Port of `ANALYSIS/ID/MetaboliteSpectralMatching.h`
//! (`MetaboliteSpectralMatching.cpp`). See `docs/ID_MAPPER_SUPPORT.md`
//! (section "MetaboliteSpectralMatching") for the API mapping and the one
//! open gap: `merge_spectra = true` needs `SpectraMerger`, which is not ported.
//!
//! The pipeline: sort the library by precursor m/z; thin every experimental
//! spectrum with a sliding [`WindowMower`](crate::processing::window_mower::WindowMower) (20 Th windows, 5 peaks); for each
//! experimental precursor, score every library spectrum whose precursor lies
//! within the precursor tolerance and whose charge sign fits the ionization
//! mode with a hyperscore
//! ([`compute_hyper_score`](crate::analysis::id::metabolite_spectral_matching::MetaboliteSpectralMatching::compute_hyper_score)); optionally drop
//! candidates whose collision cross section (CCS) differs too much; keep the
//! best one or three per precursor; and export them as mzTab small-molecule
//! rows.

use crate::concept::constants::user_param::{
    MSM_CCS, MSM_INCHI_STRING, MSM_METABOLITE_NAME, MSM_PRECURSOR_ADDUCT, MSM_SMILES_STRING,
    MSM_SUM_FORMULA,
};
use crate::format::file_info::text_format::ostream_g;
use crate::format::mztab::{
    MzTab, MzTabDouble, MzTabDoubleList, MzTabInteger, MzTabOptionalColumnEntry,
    MzTabSmallMoleculeSectionRow, MzTabString, MzTabStringList,
};
use crate::identification::PeakAnnotation;
use crate::kernel::{MSExperiment, MSSpectrum};
use crate::metadata::{DriftTimeUnit, ImTypes};
use crate::param::value::format_float;
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::processing::SpectrumFilter;
use crate::processing::window_mower::{WindowMower, WindowMowerMethod};
use crate::{Error, Result};
use std::collections::BTreeMap;

/// Name of the source `DefaultParamHandler`, `MetaboliteSpectralMatching`.
pub const METABOLITE_SPECTRAL_MATCHING_NAME: &str = "MetaboliteSpectralMatching";

/// Largest factorial the source takes, `boost::math::max_factorial<double>`.
pub const MAX_FACTORIAL: usize = 170;

/// One spectral match, the source `SpectralMatch`.
///
/// The source's getter/setter pairs become public fields. Masses are m/z
/// values in Th despite the source's "mass" names; the retention time is in
/// seconds; CCS values are in Å² and [`ImTypes::DRIFTTIME_NOT_SET`] (-1) when
/// unavailable.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralMatch {
    /// Observed precursor m/z, `observed_precursor_mass_`.
    pub observed_precursor_mass: f64,
    /// Observed RT floored to 0.1 s, `observed_precursor_rt_`.
    pub observed_precursor_rt: f64,
    /// Library precursor m/z, `found_precursor_mass_`.
    pub found_precursor_mass: f64,
    /// Library precursor charge, `found_precursor_charge_`.
    pub found_precursor_charge: i32,
    /// Hyperscore, `matching_score_`.
    pub matching_score: f64,
    /// Index of the experimental spectrum, `observed_spectrum_idx_`.
    pub observed_spectrum_index: usize,
    /// Index of the library spectrum after sorting, `matching_spectrum_idx_`.
    pub matching_spectrum_index: usize,
    /// Native ID of the experimental spectrum, `observed_spectrum_native_id_`.
    pub observed_spectrum_native_id: String,
    /// Primary library identifier, `primary_id_`.
    pub primary_identifier: String,
    /// Secondary identifier (`HMDB_ID`), `secondary_id_`.
    pub secondary_identifier: String,
    /// Common compound name, `common_name_`.
    pub common_name: String,
    /// Sum formula, `sum_formula_`.
    pub sum_formula: String,
    /// InChI string, `inchi_string_`.
    pub inchi_string: String,
    /// SMILES string, `smiles_string_`.
    pub smiles_string: String,
    /// Precursor adduct, `precursor_adduct_`.
    pub precursor_adduct: String,
    /// Experimental CCS, `observed_ccs_`.
    pub observed_ccs: f64,
    /// Library CCS, `found_ccs_`.
    pub found_ccs: f64,
}

impl Default for SpectralMatch {
    /// The source default constructor: zeros and empty strings, both CCS
    /// values [`ImTypes::DRIFTTIME_NOT_SET`].
    fn default() -> Self {
        Self {
            observed_precursor_mass: 0.0,
            observed_precursor_rt: 0.0,
            found_precursor_mass: 0.0,
            found_precursor_charge: 0,
            matching_score: 0.0,
            observed_spectrum_index: 0,
            matching_spectrum_index: 0,
            observed_spectrum_native_id: String::new(),
            primary_identifier: String::new(),
            secondary_identifier: String::new(),
            common_name: String::new(),
            sum_formula: String::new(),
            inchi_string: String::new(),
            smiles_string: String::new(),
            precursor_adduct: String::new(),
            observed_ccs: ImTypes::DRIFTTIME_NOT_SET,
            found_ccs: ImTypes::DRIFTTIME_NOT_SET,
        }
    }
}

/// Unit of both mass tolerances (parameter `mass_error_unit`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MassErrorUnit {
    /// Parts per million (`ppm`, the default).
    #[default]
    Ppm,
    /// Thomson (`Da`).
    Da,
}

/// Which results are reported per precursor (parameter `report_mode`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ReportMode {
    /// The three best-scoring matches (`top3`, the default).
    #[default]
    Top3,
    /// The best-scoring match (`best`).
    Best,
}

/// Ionization mode (parameter `ionization_mode`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum IonizationMode {
    /// Library candidates with a negative precursor charge are skipped
    /// (`positive`, the default).
    #[default]
    Positive,
    /// Library candidates with a positive precursor charge are skipped
    /// (`negative`).
    Negative,
}

/// Spectral-library search for metabolites, the source
/// `MetaboliteSpectralMatching`.
///
/// Parameters (see [`MetaboliteSpectralMatching::defaults`]):
///
/// | parameter | field | default |
/// |---|---|---|
/// | `prec_mass_error_value` | [`precursor_mz_error`](Self::precursor_mz_error) | 100.0 |
/// | `frag_mass_error_value` | [`fragment_mz_error`](Self::fragment_mz_error) | 500.0 |
/// | `mass_error_unit` (`ppm`, `Da`) | [`mass_error_unit`](Self::mass_error_unit) | `ppm` |
/// | `report_mode` (`top3`, `best`) | [`report_mode`](Self::report_mode) | `top3` |
/// | `ionization_mode` (`positive`, `negative`) | [`ionization_mode`](Self::ionization_mode) | `positive` |
/// | `merge_spectra` (`true`, `false`) | [`merge_spectra`](Self::merge_spectra) | `true` |
/// | `ccs_error_percent` (≥ 0) | [`ccs_error_percent`](Self::ccs_error_percent) | 0.0 |
///
/// The source also derives from `ProgressLogger`; this port reports no
/// progress.
#[derive(Clone, Debug, PartialEq)]
pub struct MetaboliteSpectralMatching {
    /// Precursor m/z tolerance in [`mass_error_unit`](Self::mass_error_unit),
    /// `precursor_mz_error_`.
    pub precursor_mz_error: f64,
    /// Fragment m/z tolerance in [`mass_error_unit`](Self::mass_error_unit),
    /// `fragment_mz_error_`.
    pub fragment_mz_error: f64,
    /// Unit of both tolerances, `mz_error_unit_`.
    pub mass_error_unit: MassErrorUnit,
    /// Reported matches per precursor, `report_mode_`.
    pub report_mode: ReportMode,
    /// Ionization mode, `ion_mode_`.
    pub ionization_mode: IonizationMode,
    /// Merge MS2 spectra with the same precursor mass first, `merge_spectra_`.
    /// Requires `SpectraMerger`, which is not ported:
    /// [`MetaboliteSpectralMatching::run`] refuses `true`.
    pub merge_spectra: bool,
    /// Allowed relative CCS error in percent (0 disables filtering),
    /// `ccs_error_percent_`. The relative error
    /// `|observed - library| / library * 100` must not exceed it; typical
    /// values are 1-3 %. Candidates for which either CCS is missing are kept.
    pub ccs_error_percent: f64,
}

impl Default for MetaboliteSpectralMatching {
    fn default() -> Self {
        Self {
            precursor_mz_error: 100.0,
            fragment_mz_error: 500.0,
            mass_error_unit: MassErrorUnit::Ppm,
            report_mode: ReportMode::Top3,
            ionization_mode: IonizationMode::Positive,
            merge_spectra: true,
            ccs_error_percent: 0.0,
        }
    }
}

/// `n!` rounded to the nearest `f64` (ties to even), from exact integer
/// arithmetic. The source takes `boost::math::factorial<double>`, a table of
/// the correctly rounded values for `n <= 170`.
fn factorial_f64(n: usize) -> f64 {
    // Little-endian base-2^32 limbs of n!.
    let mut limbs: Vec<u64> = vec![1];
    for k in 2..=n as u64 {
        let mut carry = 0u64;
        for limb in &mut limbs {
            let product = *limb * k + carry;
            *limb = product & 0xffff_ffff;
            carry = product >> 32;
        }
        while carry != 0 {
            limbs.push(carry & 0xffff_ffff);
            carry >>= 32;
        }
    }
    // Bit length and correctly rounded conversion.
    let top = *limbs.last().unwrap_or(&1);
    let bits = (limbs.len() - 1) * 32 + (64 - top.leading_zeros() as usize);
    if bits <= 53 {
        let mut value = 0u64;
        for &limb in limbs.iter().rev() {
            value = (value << 32) | limb;
        }
        return value as f64;
    }
    let bit = |i: usize| -> u64 { (limbs[i / 32] >> (i % 32)) & 1 };
    let shift = bits - 53;
    let mut mantissa = 0u64;
    for i in (shift..bits).rev() {
        mantissa = (mantissa << 1) | bit(i);
    }
    let round_bit = bit(shift - 1);
    let sticky = (0..shift - 1).any(|i| bit(i) == 1);
    if round_bit == 1 && (sticky || mantissa & 1 == 1) {
        mantissa += 1;
    }
    // mantissa < 2^54 and shift < 1024 for n <= 170: exact scaling.
    (mantissa as f64) * 2f64.powi(shift as i32)
}

/// Source file-local `extractCCS`: CCS in Å² from the spectrum's drift time
/// (CCS unit, or 1/K0 converted with the Mason-Schamp relation), then the
/// first precursor's, then a numeric `CCS` meta value; else
/// [`ImTypes::DRIFTTIME_NOT_SET`].
fn extract_ccs(spectrum: &MSSpectrum, mz: f64, charge: i32) -> Result<f64> {
    let unit = spectrum.drift_time_unit;
    let time = spectrum.drift_time;
    if unit == DriftTimeUnit::CollisionCrossSection && time > 0.0 {
        return Ok(time);
    }
    if unit == DriftTimeUnit::InverseReducedMobility && time > 0.0 && mz > 0.0 && charge != 0 {
        return ImTypes::one_over_k0_to_ccs(time, mz, charge);
    }
    if let Some(precursor) = spectrum.precursors.first() {
        let time = precursor.drift_time.unwrap_or(ImTypes::DRIFTTIME_NOT_SET);
        if precursor.drift_time_unit == DriftTimeUnit::CollisionCrossSection && time > 0.0 {
            return Ok(time);
        }
        if precursor.drift_time_unit == DriftTimeUnit::InverseReducedMobility
            && time > 0.0
            && mz > 0.0
            && charge != 0
        {
            return ImTypes::one_over_k0_to_ccs(time, mz, charge);
        }
    }
    if let Some(value) = spectrum.metadata.get(MSM_CCS) {
        if let Some(ccs) = source_to_double(&value.to_string()) {
            return Ok(ccs);
        }
    }
    Ok(ImTypes::DRIFTTIME_NOT_SET)
}

/// Source `StringUtils::toDouble` for the decimal text a CCS meta value
/// holds: [`crate::concept::string_conversion::to_double`].
fn source_to_double(text: &str) -> Option<f64> {
    crate::concept::string_conversion::to_double(text).ok()
}

/// A string meta value as the source's implicit `DataValue` to `std::string`
/// conversion yields it; absent or non-string values throw there.
fn string_meta(spectrum: &MSSpectrum, key: &str) -> Result<String> {
    spectrum
        .metadata
        .get(key)
        .and_then(|value| value.as_str().ok())
        .map(str::to_owned)
        .ok_or_else(|| {
            Error::InvalidValue(format!(
                "MetaboliteSpectralMatching: library spectrum needs a string meta value '{key}'"
            ))
        })
}

fn text_meta(spectrum: &MSSpectrum, key: &str) -> Option<String> {
    spectrum.metadata.get(key).map(ToString::to_string)
}

impl MetaboliteSpectralMatching {
    /// The default parameters (source constructor), in source order.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters of this matcher in the layout of
    /// [`MetaboliteSpectralMatching::defaults`] (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        let float = |name: &str, value: f64, description: &str| ParamEntry {
            name: name.into(),
            description: description.into(),
            value: ParamValue::Float(value),
            ..ParamEntry::default()
        };
        let choice = |name: &str, value: &str, description: &str, valid: &[&str]| ParamEntry {
            name: name.into(),
            description: description.into(),
            value: ParamValue::String(value.into()),
            valid_strings: valid.iter().map(|v| (*v).to_owned()).collect(),
            ..ParamEntry::default()
        };
        let mut ccs = float(
            "ccs_error_percent",
            self.ccs_error_percent,
            "Allowed collision cross section (CCS) error in percent for filtering candidates (0 = filtering disabled). The relative error |observed - library| / library * 100 must not exceed this value. Typical values for confident identification are 1-3%. Candidates for which the experimental or the library CCS is missing are kept (not filtered). Requires CCS to be present as CCS-unit drift time, 1/K0 (converted to CCS), or a numeric 'CCS' meta value.",
        );
        ccs.min_float = 0.0;
        let entries = vec![
            float(
                "prec_mass_error_value",
                self.precursor_mz_error,
                "Error allowed for precursor ion mass.",
            ),
            float(
                "frag_mass_error_value",
                self.fragment_mz_error,
                "Error allowed for product ions.",
            ),
            choice(
                "mass_error_unit",
                match self.mass_error_unit {
                    MassErrorUnit::Ppm => "ppm",
                    MassErrorUnit::Da => "Da",
                },
                "Unit of mass error (ppm or Da)",
                &["ppm", "Da"],
            ),
            choice(
                "report_mode",
                match self.report_mode {
                    ReportMode::Top3 => "top3",
                    ReportMode::Best => "best",
                },
                "Which results shall be reported: the top-three scoring ones or the best scoring one?",
                &["top3", "best"],
            ),
            choice(
                "ionization_mode",
                match self.ionization_mode {
                    IonizationMode::Positive => "positive",
                    IonizationMode::Negative => "negative",
                },
                "Positive or negative ionization mode?",
                &["positive", "negative"],
            ),
            choice(
                "merge_spectra",
                if self.merge_spectra { "true" } else { "false" },
                "Merge MS2 spectra with the same precursor mass.",
                &["true", "false"],
            ),
            ccs,
        ];
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries,
            nodes: Vec::new(),
        })
    }

    /// Build a matcher from parameters, as source `setParameters` followed by
    /// `updateMembers_`. Missing parameters take their defaults; unknown ones
    /// are ignored (see
    /// [`MetaboliteSpectralMatching::from_param_with_warnings`]).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a wrong value type or a violated
    /// restriction (the source `Exception::InvalidParameter`).
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`MetaboliteSpectralMatching::from_param`], also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`MetaboliteSpectralMatching::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(METABOLITE_SPECTRAL_MATCHING_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, Self::from_complete_param)
    }

    fn from_complete_param(param: &Param) -> Result<Self> {
        let float = |key: &str| -> Result<f64> {
            match param.value(key)? {
                ParamValue::Float(value) => Ok(*value),
                _ => Err(Error::InvalidValue(format!(
                    "parameter '{key}' must be a floating-point value"
                ))),
            }
        };
        let text = |key: &str| -> Result<String> { Ok(param.value(key)?.as_str()?.to_owned()) };
        Ok(Self {
            precursor_mz_error: float("prec_mass_error_value")?,
            fragment_mz_error: float("frag_mass_error_value")?,
            // The source compares `mz_error_unit_ == "Da"`; anything else is ppm.
            mass_error_unit: if text("mass_error_unit")? == "Da" {
                MassErrorUnit::Da
            } else {
                MassErrorUnit::Ppm
            },
            report_mode: if text("report_mode")? == "best" {
                ReportMode::Best
            } else {
                ReportMode::Top3
            },
            ionization_mode: if text("ionization_mode")? == "negative" {
                IonizationMode::Negative
            } else {
                IonizationMode::Positive
            },
            merge_spectra: param.value("merge_spectra")?.to_bool()?,
            ccs_error_percent: float("ccs_error_percent")?,
        })
    }

    /// Hyperscore of an experimental spectrum against a library spectrum, the
    /// source static `computeHyperScore(double, bool, const MSSpectrum&, const
    /// MSSpectrum&, double)`.
    ///
    /// Library peaks between `max(mz_lower_bound, lowest experimental m/z -
    /// tolerance)` and `highest experimental m/z + tolerance` are each matched
    /// to the nearest experimental peak within the fragment tolerance
    /// (`fragment_mass_error` in ppm of the library m/z when
    /// `fragment_mass_tolerance_unit_ppm`, else Th); several library peaks may
    /// match one experimental peak. With `N` matched experimental peaks and the
    /// dot product `D` of each matched experimental intensity with the largest
    /// intensity of its library peaks, the score is `ln(D) + ln(min(N,
    /// 170)!)`, clamped at 0; fewer than three matches score 0, as does an
    /// empty experimental spectrum.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] when either spectrum is not sorted by
    /// m/z (the source assumes it) and [`Error::InvalidValue`] for a
    /// non-finite bound or a negative tolerance.
    pub fn compute_hyper_score(
        fragment_mass_error: f64,
        fragment_mass_tolerance_unit_ppm: bool,
        exp_spectrum: &MSSpectrum,
        db_spectrum: &MSSpectrum,
        mz_lower_bound: f64,
    ) -> Result<f64> {
        hyper_score(
            fragment_mass_error,
            fragment_mass_tolerance_unit_ppm,
            exp_spectrum,
            db_spectrum,
            None,
            mz_lower_bound,
        )
    }

    /// [`MetaboliteSpectralMatching::compute_hyper_score`], also appending one
    /// peak annotation per matched library peak to `annotations`, the source's
    /// overload with a `std::vector<PeptideHit::PeakAnnotation>&`.
    ///
    /// Annotations are produced only when the library spectrum has a string and
    /// an integer data array: the annotation text is the first string array's
    /// value and the charge the first integer array's value at the library
    /// peak, the m/z and intensity those of the matched experimental peak.
    ///
    /// # Errors
    ///
    /// As [`MetaboliteSpectralMatching::compute_hyper_score`], plus
    /// [`Error::InvalidValue`] when those data arrays are shorter than the
    /// library spectrum (the source's `std::vector::at` throws).
    pub fn compute_hyper_score_with_annotations(
        fragment_mass_error: f64,
        fragment_mass_tolerance_unit_ppm: bool,
        exp_spectrum: &MSSpectrum,
        db_spectrum: &MSSpectrum,
        annotations: &mut Vec<PeakAnnotation>,
        mz_lower_bound: f64,
    ) -> Result<f64> {
        hyper_score(
            fragment_mass_error,
            fragment_mass_tolerance_unit_ppm,
            exp_spectrum,
            db_spectrum,
            Some(annotations),
            mz_lower_bound,
        )
    }

    /// Search `msexp` against the spectral library `spec_db` and write the
    /// matches into the small-molecule section of `mztab_out`, the source
    /// `run(PeakMap&, PeakMap&, MzTab&, std::string& out_spectra)`.
    ///
    /// As in the source, `spec_db` is sorted by precursor m/z and every
    /// spectrum of `msexp` is thinned by a sliding window mower (20 Th, 5
    /// peaks) in place. When `out_spectra` is not empty the thinned spectra are
    /// stored there as mzML. The precursor window uses
    /// [`precursor_mz_error`](Self::precursor_mz_error) in the same unit as the
    /// fragment tolerance. Candidates whose precursor charge contradicts the
    /// ionization mode are skipped; with CCS filtering enabled, candidates
    /// whose CCS error exceeds the limit are skipped when both CCS values are
    /// known. Matches with a positive hyperscore are ranked per precursor by
    /// decreasing score (equal scores keep library order; the source's
    /// `std::sort` leaves them unspecified) and the best one or three are
    /// reported. The small-molecule rows of `mztab_out` are replaced; its other
    /// sections are untouched.
    ///
    /// Returns the source's warnings (CCS filtering enabled but no
    /// experimental CCS found).
    ///
    /// # Errors
    ///
    /// Returns [`Error::Unsupported`] when [`merge_spectra`](Self::merge_spectra)
    /// is set (`SpectraMerger` is not ported) or `out_spectra` is given without
    /// the `mzml` feature; [`Error::InvalidValue`] when a library spectrum has
    /// no precursor (the source reads `getPrecursors()[0]` unchecked) or a
    /// scored library spectrum lacks one of the string meta values `HMDB_ID`,
    /// `Sum_Formula`, `Inchi_String`, `SMILES_String`, `Precursor_Ion` (the
    /// source's implicit conversion throws); errors of the window mower and of
    /// [`MetaboliteSpectralMatching::compute_hyper_score`]; and mzML write
    /// errors. On error `mztab_out` is unchanged.
    pub fn run(
        &self,
        msexp: &mut MSExperiment,
        spec_db: &mut MSExperiment,
        mztab_out: &mut MzTab,
        out_spectra: &str,
    ) -> Result<Vec<String>> {
        if self.merge_spectra {
            return Err(Error::Unsupported(
                "MetaboliteSpectralMatching: merge_spectra needs SpectraMerger, which is not ported; set merge_spectra to false".into(),
            ));
        }
        let mut keys = Vec::with_capacity(spec_db.spectra.len());
        for spectrum in &spec_db.spectra {
            let precursor = spectrum.precursors.first().ok_or_else(|| {
                Error::InvalidValue(
                    "MetaboliteSpectralMatching: every library spectrum needs a precursor".into(),
                )
            })?;
            if precursor.mz.is_nan() {
                return Err(Error::InvalidValue(
                    "MetaboliteSpectralMatching: library precursor m/z is NaN".into(),
                ));
            }
            keys.push(precursor.mz);
        }
        // PrecursorMZLess; a stable sort keeps equal precursors in input order.
        let mut order: Vec<usize> = (0..spec_db.spectra.len()).collect();
        order.sort_by(|&a, &b| keys[a].total_cmp(&keys[b]));
        let mut sorted = Vec::with_capacity(order.len());
        let mut slots: Vec<Option<MSSpectrum>> = spec_db.spectra.drain(..).map(Some).collect();
        for index in order {
            if let Some(spectrum) = slots[index].take() {
                sorted.push(spectrum);
            }
        }
        spec_db.spectra = sorted;
        let mz_keys: Vec<f64> = spec_db.spectra.iter().map(|s| s.precursors[0].mz).collect();

        let mower = WindowMower {
            window_size: 20.0,
            peak_count: 5,
            method: WindowMowerMethod::Sliding,
            ..WindowMower::default()
        };
        mower.filter_experiment(msexp)?;

        if !out_spectra.is_empty() {
            store_spectra(out_spectra, msexp)?;
        }

        let ppm = self.mass_error_unit == MassErrorUnit::Ppm;
        let ccs_filtering = self.ccs_error_percent > 0.0;
        let mut any_observed_ccs = false;
        let mut results: Vec<SpectralMatch> = Vec::new();
        for (spec_idx, spectrum) in msexp.spectra.iter().enumerate() {
            for precursor in &spectrum.precursors {
                let precursor_mz = precursor.mz;
                let observed_ccs = extract_ccs(spectrum, precursor_mz, precursor.charge)?;
                if observed_ccs > 0.0 {
                    any_observed_ccs = true;
                }
                let offset = if ppm {
                    precursor_mz * 1e-6 * self.precursor_mz_error
                } else {
                    self.precursor_mz_error
                };
                let lower = precursor_mz - offset;
                let upper = precursor_mz + offset;
                let start = mz_keys.partition_point(|&k| k < lower);
                let end = mz_keys.partition_point(|&k| k <= upper);
                let mut partial = Vec::new();
                for search_idx in start..end.max(start) {
                    let candidate = &spec_db.spectra[search_idx];
                    let db_precursor = &candidate.precursors[0];
                    let wrong_sign = match self.ionization_mode {
                        IonizationMode::Positive => db_precursor.charge < 0,
                        IonizationMode::Negative => db_precursor.charge > 0,
                    };
                    if wrong_sign {
                        continue;
                    }
                    let library_ccs = extract_ccs(candidate, db_precursor.mz, db_precursor.charge)?;
                    if ccs_filtering && observed_ccs > 0.0 && library_ccs > 0.0 {
                        let error = (observed_ccs - library_ccs).abs() / library_ccs * 100.0;
                        if error > self.ccs_error_percent {
                            continue;
                        }
                    }
                    let score = Self::compute_hyper_score(
                        self.fragment_mz_error,
                        ppm,
                        spectrum,
                        candidate,
                        0.0,
                    )?;
                    if score > 0.0 {
                        let gnps = text_meta(candidate, "GNPS_Spectrum_ID");
                        let common_name = text_meta(candidate, "Metabolite_Name")
                            .or_else(|| gnps.clone())
                            .unwrap_or_default();
                        let primary = gnps
                            .clone()
                            .or_else(|| text_meta(candidate, "Massbank_Accession_ID"))
                            .or_else(|| text_meta(candidate, MSM_METABOLITE_NAME))
                            .unwrap_or_else(|| candidate.native_id.clone());
                        partial.push(SpectralMatch {
                            observed_precursor_mass: precursor_mz,
                            found_precursor_mass: db_precursor.mz,
                            observed_precursor_rt: (spectrum.rt * 10.0).floor() / 10.0,
                            found_precursor_charge: db_precursor.charge,
                            matching_score: score,
                            observed_spectrum_index: spec_idx,
                            matching_spectrum_index: search_idx,
                            observed_spectrum_native_id: spectrum.native_id.clone(),
                            primary_identifier: primary,
                            secondary_identifier: string_meta(candidate, "HMDB_ID")?,
                            sum_formula: string_meta(candidate, MSM_SUM_FORMULA)?,
                            common_name,
                            inchi_string: string_meta(candidate, MSM_INCHI_STRING)?,
                            smiles_string: string_meta(candidate, MSM_SMILES_STRING)?,
                            precursor_adduct: string_meta(candidate, MSM_PRECURSOR_ADDUCT)?,
                            observed_ccs,
                            found_ccs: library_ccs,
                        });
                    }
                }
                partial.sort_by(|a, b| b.matching_score.total_cmp(&a.matching_score));
                let keep = match self.report_mode {
                    ReportMode::Top3 => 3,
                    ReportMode::Best => 1,
                };
                results.extend(partial.into_iter().take(keep));
            }
        }
        let mut warnings = Vec::new();
        if ccs_filtering && !any_observed_ccs {
            warnings.push(format!(
                "Warning: 'ccs_error_percent' is set to {} but no experimental CCS could be determined for any spectrum (no CCS-unit or 1/K0 drift time and no numeric 'CCS' meta value). CCS filtering had no effect.",
                ostream_g(self.ccs_error_percent, 6)
            ));
        }
        mztab_out.small_molecule_data = export_mztab(&results);
        Ok(warnings)
    }
}

#[cfg(feature = "mzml")]
fn store_spectra(path: &str, experiment: &MSExperiment) -> Result<()> {
    crate::format::mzml::store(path, experiment)
}

#[cfg(not(feature = "mzml"))]
fn store_spectra(_path: &str, _experiment: &MSExperiment) -> Result<()> {
    Err(Error::Unsupported(
        "MetaboliteSpectralMatching: storing the spectra needs the 'mzml' feature".into(),
    ))
}

fn hyper_score(
    fragment_mass_error: f64,
    ppm: bool,
    exp_spectrum: &MSSpectrum,
    db_spectrum: &MSSpectrum,
    annotations: Option<&mut Vec<PeakAnnotation>>,
    mz_lower_bound: f64,
) -> Result<f64> {
    let (Some(first), Some(last)) = (exp_spectrum.peaks.first(), exp_spectrum.peaks.last()) else {
        return Ok(0.0);
    };
    let offset_at = |mz: f64| {
        if ppm {
            mz * fragment_mass_error * 1e-6
        } else {
            fragment_mass_error
        }
    };
    let lower = mz_lower_bound.max(first.mz - offset_at(first.mz));
    let upper = last.mz + offset_at(last.mz);
    let begin = db_spectrum.mz_begin(lower)?;
    let end = db_spectrum.mz_end(upper)?.max(begin);
    let mut matches: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for db_index in begin..end {
        let db_mz = db_spectrum.peaks[db_index].mz;
        if let Some(index) = exp_spectrum.find_nearest_with_tolerance(db_mz, offset_at(db_mz))? {
            matches.entry(index).or_default().push(db_index);
        }
    }
    let mut dot_product = 0.0;
    for (&exp_index, db_indices) in &matches {
        let db_intensity = db_indices
            .iter()
            .map(|&i| f64::from(db_spectrum.peaks[i].intensity))
            .fold(0.0, f64::max);
        dot_product += db_intensity * f64::from(exp_spectrum.peaks[exp_index].intensity);
    }
    if let Some(annotations) = annotations {
        if let (Some(texts), Some(charges)) = (
            db_spectrum.string_data_arrays.first(),
            db_spectrum.integer_data_arrays.first(),
        ) {
            let mut added = Vec::new();
            for (&exp_index, db_indices) in &matches {
                let peak = &exp_spectrum.peaks[exp_index];
                for &db_index in db_indices {
                    let (Some(text), Some(&charge)) =
                        (texts.data.get(db_index), charges.data.get(db_index))
                    else {
                        return Err(Error::InvalidValue(
                            "MetaboliteSpectralMatching: library data arrays are shorter than the spectrum".into(),
                        ));
                    };
                    added.push(PeakAnnotation {
                        annotation: text.clone(),
                        charge,
                        mz: peak.mz,
                        intensity: f64::from(peak.intensity),
                    });
                }
            }
            annotations.extend(added);
        }
    }
    let matched = matches.len();
    if matched < 3 {
        return Ok(0.0);
    }
    let term = factorial_f64(matched.min(MAX_FACTORIAL)).ln();
    let score = dot_product.ln() + term;
    Ok(if score < 0.0 { 0.0 } else { score })
}

/// Source `exportMzTab_`.
fn export_mztab(results: &[SpectralMatch]) -> Vec<MzTabSmallMoleculeSectionRow> {
    let text = |value: &str| MzTabString::from_text(value);
    let two_decimals = |value: f64| (value * 100.0).floor() / 100.0;
    results
        .iter()
        .map(|m| {
            let mut row = MzTabSmallMoleculeSectionRow::default();
            let mut identifier = MzTabStringList::default();
            identifier.set(vec![text(&m.primary_identifier)]);
            row.identifier = identifier;
            row.chemical_formula = text(&m.sum_formula);
            row.smiles = text(&m.smiles_string);
            row.inchi_key = text(&m.inchi_string);
            row.description = text(&m.common_name);
            row.exp_mass_to_charge = MzTabDouble::new(m.found_precursor_mass);
            row.charge = MzTabInteger::new(m.found_precursor_charge);
            let mut rt = MzTabDoubleList::default();
            rt.set(vec![MzTabDouble::new(m.observed_precursor_rt)]);
            row.retention_time = rt;
            row.database = text("MassBank");
            row.database_version = text("Sep 27, 2013");
            row.smallmolecule_abundance_study_variable
                .insert(1, MzTabDouble::new(0.0));
            row.smallmolecule_abundance_stdev_study_variable
                .insert(1, MzTabDouble::new(0.0));
            row.smallmolecule_abundance_std_error_study_variable
                .insert(1, MzTabDouble::new(0.0));
            let ppm_error = two_decimals(
                (m.found_precursor_mass - m.observed_precursor_mass) / m.found_precursor_mass * 1e6,
            );
            let ccs = |value: f64| {
                if value > 0.0 {
                    format_float(two_decimals(value), true)
                } else {
                    "null".to_owned()
                }
            };
            let ccs_error = if m.observed_ccs > 0.0 && m.found_ccs > 0.0 {
                format_float(
                    two_decimals((m.observed_ccs - m.found_ccs).abs() / m.found_ccs * 100.0),
                    true,
                )
            } else {
                "null".to_owned()
            };
            row.opt = vec![
                MzTabOptionalColumnEntry::new(
                    "opt_ppm_error",
                    text(&format_float(ppm_error, true)),
                ),
                MzTabOptionalColumnEntry::new("opt_adduct_ion", text(&m.precursor_adduct)),
                MzTabOptionalColumnEntry::new(
                    "opt_match_score",
                    text(&ostream_g(m.matching_score, 6)),
                ),
                MzTabOptionalColumnEntry::new("opt_sec_id", text(&m.secondary_identifier)),
                MzTabOptionalColumnEntry::new(
                    "opt_source_idx",
                    text(&m.observed_spectrum_index.to_string()),
                ),
                MzTabOptionalColumnEntry::new(
                    "opt_spec_native_id",
                    text(&m.observed_spectrum_native_id),
                ),
                MzTabOptionalColumnEntry::new("opt_observed_ccs", text(&ccs(m.observed_ccs))),
                MzTabOptionalColumnEntry::new("opt_library_ccs", text(&ccs(m.found_ccs))),
                MzTabOptionalColumnEntry::new("opt_ccs_error_percent", text(&ccs_error)),
            ];
            row
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::factorial_f64;

    #[test]
    fn factorials_are_correctly_rounded() {
        assert_eq!(factorial_f64(0), 1.0);
        assert_eq!(factorial_f64(5), 120.0);
        assert_eq!(factorial_f64(20), 2_432_902_008_176_640_000.0);
        // 25! = 15511210043330985984000000, which rounds to this double.
        assert_eq!(factorial_f64(25), 1.5511210043330986e25);
        // 170! is the largest finite factorial in f64.
        assert!(factorial_f64(170).is_finite());
        assert!((factorial_f64(170).ln() - 706.5730622457874).abs() < 1e-9);
    }
}
