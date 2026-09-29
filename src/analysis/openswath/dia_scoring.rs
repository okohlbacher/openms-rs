// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Scoring of a spectrum against a transition group in DIA (SWATH) analysis,
//! from `ANALYSIS/OPENSWATH/DIAScoring.h` (class `DIAScoring`).
//!
//! In DIA analysis a full MS2 spectrum is recorded at every chromatographic
//! point. [`DiaScoring`] computes scores from that spectrum:
//!
//! - **Isotope scores** ([`DiaScoring::dia_isotope_scores`]):
//!   - *isotope correlation*: the Pearson correlation of each fragment ion's
//!     observed isotopes with the theoretical averagine distribution, weighted
//!     by the relative intensity of the transition (more is better);
//!   - *isotope overlap*: whether a signal at `(mz - 1) / charge` exists that is
//!     stronger than the current peak, which would indicate the peak is itself
//!     an isotope of another peak. This counts how often such a peak was seen,
//!     weighted by the relative intensity of the transition (less is better).
//! - **Mass-difference score** ([`DiaScoring::dia_massdiff_score`]): the ppm
//!   difference of the experimental to the expected signal,
//!   `sum(|ppm difference|) / #observed transitions` (less is better). A list
//!   of all ppm differences is also returned, `-1` where no signal was found.
//! - **b/y ion score** ([`DiaScoring::dia_by_ion_score`]): presence of the b/y
//!   ions of the peptide.
//! - **Theoretical spectrum** ([`DiaScoring::score_with_isotopes`]): a dot
//!   product and a Manhattan score against a theoretical spectrum.
//!
//! Spectra are the light OpenSWATH spectra
//! ([`Spectrum`](crate::analysis::openswath::data_structures::Spectrum));
//! transitions are the light transition format
//! ([`LightTransition`](crate::analysis::openswath::transition_experiment::LightTransition)).
//!
//! The source memoises b/y series, prescore spectra and averagine distributions
//! in thread-local pools; this port recomputes them, which yields the same
//! values without hidden per-thread state.
//!
//! See `docs/OPENSWATH_SCORING_SUPPORT.md`.

use super::data_structures::SpectrumPtr;
use super::dia_helper::{
    WindowIntegral, by_series, create_mz_range_ppm, integrate_window_sequence,
};
use super::dia_prescoring::{DiaPrescore, PrescoreResult};
use super::itransition::MRMFeatureAccess;
use super::transition_experiment::LightTransition;
use crate::chemistry::isotopes::ProbabilityPrecision;
use crate::chemistry::{
    AASequence, CoarseIsotopePatternGenerator, CoarseMassMode, EmpiricalFormula,
    IsotopeDistribution, TheoreticalSpectrumGenerator,
};
use crate::concept::math_functions::{ppm, ppm_abs};
use crate::kernel::ranges::RangeBase;
use crate::math::statistic_functions::pearson_correlation_coefficient;
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::{Error, Result};

/// The `DefaultParamHandler` name of [`DiaScoring`] (`"DIAScoring"`).
pub const DIA_SCORING_NAME: &str = "DIAScoring";

/// The C13-C12 spacing `DIAScoring.cpp` uses.
///
/// The source file defines its own global `C13C12_MASSDIFF_U = 1.0033548`,
/// which unqualified lookup inside `namespace OpenMS` resolves to instead of
/// `Constants::C13C12_MASSDIFF_U = 1.0033548378`. `DIAPrescoring.cpp` uses the
/// latter. Both values are kept where the source uses them.
pub const DIA_SCORING_C13C12_MASSDIFF_U: f64 = 1.0033548;

/// Isotope correlation and overlap scores.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IsotopeScores {
    /// Pearson correlation against the theoretical isotope distribution
    /// (higher is better).
    pub isotope_corr: f64,
    /// Overlap with a putative preceding isotope envelope (lower is better).
    pub isotope_overlap: f64,
}

/// Fragment mass-error scores of [`DiaScoring::dia_massdiff_score`].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MassDiffScores {
    /// Mean absolute ppm error over the observed transitions (0 when none was
    /// observed).
    pub ppm_score: f64,
    /// Sum of `|ppm| * normalized_library_intensity` over the observed
    /// transitions.
    pub ppm_score_weighted: f64,
    /// Signed ppm error per transition; `-1` marks a transition without signal.
    pub diff_ppm: Vec<f64>,
}

/// The precursor mass-error score of [`DiaScoring::dia_ms1_massdiff_score`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ms1MassDiffScore {
    /// On success the absolute ppm error of the measured peak; otherwise the
    /// absolute ppm width of the extraction window, a worst-case bound (not
    /// `-1`).
    pub ppm_score: f64,
    /// Whether a signal was found in the extraction window (the source's
    /// return value).
    pub signal_found: bool,
}

/// Numbers of matched b- and y-ions of [`DiaScoring::dia_by_ion_score`].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ByIonScores {
    /// Number of matching b-ions.
    pub bseries_score: f64,
    /// Number of matching y-ions.
    pub yseries_score: f64,
}

/// DIA scoring of a spectrum against a transition group (source `DIAScoring`).
///
/// # Parameters
///
/// | name | default | meaning |
/// |---|---|---|
/// | `dia_extraction_window` | 0.05 | DIA extraction window in Th or ppm (min 0) |
/// | `dia_extraction_unit` | `Th` | window unit, `Th` or `ppm` |
/// | `dia_centroided` | `false` | use centroided DIA data |
/// | `dia_byseries_intensity_min` | 300.0 | minimum b/y intensity to consider (min 0) |
/// | `dia_byseries_ppm_diff` | 10.0 | maximum b/y ppm difference to consider (min 0) |
/// | `dia_nr_isotopes` | 4 | number of isotopes to consider (min 0) |
/// | `dia_nr_charges` | 4 | number of charges to consider (min 0) |
/// | `peak_before_mono_max_ppm_diff` | 20.0 | maximal ppm difference to count a peak at lower m/z as evidence that a peak is not monoisotopic (min 0) |
///
/// The source also constructs a `TheoreticalSpectrumGenerator` with
/// `add_metainfo = true` for the b/y series; it is not a parameter.
#[derive(Clone, Debug, PartialEq)]
pub struct DiaScoring {
    /// DIA extraction window, in Th or ppm by [`Self::dia_extraction_ppm`].
    pub dia_extract_window: f64,
    /// Whether [`Self::dia_extract_window`] is in ppm (`dia_extraction_unit`).
    pub dia_extraction_ppm: bool,
    /// Use centroided DIA data; window integration of centroided data is not
    /// implemented and fails, as in the source.
    pub dia_centroided: bool,
    /// Minimum intensity for a b/y ion to count.
    pub dia_byseries_intensity_min: f64,
    /// Maximum ppm difference for a b/y ion to count.
    pub dia_byseries_ppm_diff: f64,
    /// Number of isotopes beyond the monoisotope that are integrated.
    pub dia_nr_isotopes: u32,
    /// Number of charge states tested for a preceding isotope peak.
    pub dia_nr_charges: u32,
    /// Maximum ppm difference to count a peak at lower m/z.
    pub peak_before_mono_max_ppm_diff: f64,
    generator: TheoreticalSpectrumGenerator,
}

impl Default for DiaScoring {
    fn default() -> Self {
        Self {
            dia_extract_window: 0.05,
            dia_extraction_ppm: false,
            dia_centroided: false,
            dia_byseries_intensity_min: 300.0,
            dia_byseries_ppm_diff: 10.0,
            dia_nr_isotopes: 4,
            dia_nr_charges: 4,
            peak_before_mono_max_ppm_diff: 20.0,
            generator: TheoreticalSpectrumGenerator {
                add_metainfo: true,
                ..TheoreticalSpectrumGenerator::default()
            },
        }
    }
}

pub(crate) fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

pub(crate) fn param_entry(name: &str, value: ParamValue, description: &str) -> ParamEntry {
    ParamEntry {
        name: name.into(),
        description: description.into(),
        value,
        ..ParamEntry::default()
    }
}

pub(crate) fn float_param(
    name: &str,
    value: f64,
    description: &str,
    min: Option<f64>,
) -> ParamEntry {
    let mut e = param_entry(name, ParamValue::Float(value), description);
    if let Some(min) = min {
        e.min_float = min;
    }
    e
}

pub(crate) fn int_param(name: &str, value: i64, description: &str, min: Option<i32>) -> ParamEntry {
    let mut e = param_entry(name, ParamValue::Integer(value), description);
    if let Some(min) = min {
        e.min_int = min;
    }
    e
}

pub(crate) fn string_param(
    name: &str,
    value: &str,
    description: &str,
    valid: &[&str],
) -> ParamEntry {
    let mut e = param_entry(name, ParamValue::String(value.into()), description);
    e.valid_strings = valid.iter().map(|v| (*v).to_owned()).collect();
    e
}

pub(crate) fn flag_param(name: &str, value: bool, description: &str) -> ParamEntry {
    string_param(
        name,
        if value { "true" } else { "false" },
        description,
        &["true", "false"],
    )
}

pub(crate) fn get_float(param: &Param, key: &str) -> Result<f64> {
    param.value(key)?.to_f64()
}

pub(crate) fn get_flag(param: &Param, key: &str) -> Result<bool> {
    param.value(key)?.to_bool()
}

pub(crate) fn get_u32(param: &Param, key: &str) -> Result<u32> {
    let value = param.value(key)?.to_i64()?;
    u32::try_from(value)
        .ok()
        .filter(|&v| i32::try_from(v).is_ok())
        .ok_or_else(|| bad(format!("'{key}' must be a non-negative 32-bit integer")))
}

impl DiaScoring {
    /// A scorer with the source defaults (source `DIAScoring()`).
    pub fn new() -> Self {
        Self::default()
    }

    /// The source defaults as a [`Param`] (source `getDefaults`).
    ///
    /// # Errors
    ///
    /// Propagates a [`Param`] construction failure, which does not occur for
    /// these fixed entries.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The current parameters in the layout of [`DiaScoring::defaults`]
    /// (source `getParameters`).
    ///
    /// # Errors
    ///
    /// As [`DiaScoring::defaults`].
    pub fn to_param(&self) -> Result<Param> {
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries: vec![
                float_param(
                    "dia_extraction_window",
                    self.dia_extract_window,
                    "DIA extraction window in Th or ppm.",
                    Some(0.0),
                ),
                string_param(
                    "dia_extraction_unit",
                    if self.dia_extraction_ppm { "ppm" } else { "Th" },
                    "DIA extraction window unit",
                    &["Th", "ppm"],
                ),
                flag_param(
                    "dia_centroided",
                    self.dia_centroided,
                    "Use centroided DIA data.",
                ),
                float_param(
                    "dia_byseries_intensity_min",
                    self.dia_byseries_intensity_min,
                    "DIA b/y series minimum intensity to consider.",
                    Some(0.0),
                ),
                float_param(
                    "dia_byseries_ppm_diff",
                    self.dia_byseries_ppm_diff,
                    "DIA b/y series minimal difference in ppm to consider.",
                    Some(0.0),
                ),
                int_param(
                    "dia_nr_isotopes",
                    i64::from(self.dia_nr_isotopes),
                    "DIA number of isotopes to consider.",
                    Some(0),
                ),
                int_param(
                    "dia_nr_charges",
                    i64::from(self.dia_nr_charges),
                    "DIA number of charges to consider.",
                    Some(0),
                ),
                float_param(
                    "peak_before_mono_max_ppm_diff",
                    self.peak_before_mono_max_ppm_diff,
                    "DIA maximal difference in ppm to count a peak at lower m/z when searching for evidence that a peak might not be monoisotopic.",
                    Some(0.0),
                ),
            ],
            nodes: Vec::new(),
        })
    }

    /// Build a scorer from parameters, as source `setParameters` followed by
    /// `updateMembers_`. Missing parameters take their defaults.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a parameter has the wrong type or violates
    /// its restriction (the source `Exception::InvalidParameter`).
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`DiaScoring::from_param`], also returning the unknown-parameter
    /// warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`DiaScoring::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(DIA_SCORING_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, |param| {
            Ok(Self {
                dia_extract_window: get_float(param, "dia_extraction_window")?,
                dia_extraction_ppm: param.value("dia_extraction_unit")?.as_str()? == "ppm",
                dia_centroided: get_flag(param, "dia_centroided")?,
                dia_byseries_intensity_min: get_float(param, "dia_byseries_intensity_min")?,
                dia_byseries_ppm_diff: get_float(param, "dia_byseries_ppm_diff")?,
                dia_nr_isotopes: get_u32(param, "dia_nr_isotopes")?,
                dia_nr_charges: get_u32(param, "dia_nr_charges")?,
                peak_before_mono_max_ppm_diff: get_float(param, "peak_before_mono_max_ppm_diff")?,
                ..Self::default()
            })
        })
    }

    /// The DIA extraction window, in Th or ppm depending on
    /// [`DiaScoring::is_dia_extraction_ppm`] (source `getDIAExtractionWindow`).
    pub fn dia_extraction_window(&self) -> f64 {
        self.dia_extract_window
    }

    /// Whether the DIA extraction window is in ppm; `false` means Th (source
    /// `isDIAExtractionPPM`).
    pub fn is_dia_extraction_ppm(&self) -> bool {
        self.dia_extraction_ppm
    }

    fn integrate_at(
        &self,
        spectra: &[SpectrumPtr],
        center: f64,
        im_range: &RangeBase,
    ) -> Result<Option<WindowIntegral>> {
        let range = create_mz_range_ppm(center, self.dia_extract_window, self.dia_extraction_ppm)?;
        integrate_window_sequence(spectra, &range, im_range, self.dia_centroided)
    }

    /// Integrated intensities of the monoisotope and `dia_nr_isotopes` further
    /// isotopes at spacing `1.0033548 / |charge|` (source
    /// `getIsotopeIntysFromExpSpec_` and the loop in `diaIsotopeScoresSub_`).
    ///
    /// A non-finite isotope position (charge 0, where the source divides by
    /// zero) contributes intensity 0, as the source's NaN or infinite window
    /// finds no signal.
    fn isotope_intensities(
        &self,
        mono_mz: f64,
        spectra: &[SpectrumPtr],
        charge: i32,
        im_range: &RangeBase,
    ) -> Result<Vec<f64>> {
        let abs_charge = f64::from(charge).abs();
        let mut out = Vec::with_capacity(self.dia_nr_isotopes as usize + 1);
        for iso in 0..=self.dia_nr_isotopes {
            let center = mono_mz + f64::from(iso) * DIA_SCORING_C13C12_MASSDIFF_U / abs_charge;
            let intensity = if center.is_finite() {
                self.integrate_at(spectra, center, im_range)?
                    .map_or(0.0, |found| found.intensity)
            } else {
                0.0
            };
            out.push(intensity);
        }
        Ok(out)
    }

    /// Isotope-pattern correlation and isotope-overlap scores of a transition
    /// group (source `dia_isotope_scores`, which delegates to
    /// `diaIsotopeScoresSub_`).
    ///
    /// For every transition (fragment charge 1 when unset), the monoisotope and
    /// `dia_nr_isotopes` isotopes are integrated and correlated with the
    /// averagine distribution of weight `product_mz * |charge|`; the correlation
    /// is weighted by the transition's apex intensity relative to the peak
    /// group intensity. The overlap counts, per transition and weighted
    /// likewise, how often a peak one isotope spacing below the monoisotope is
    /// both more intense and within `peak_before_mono_max_ppm_diff`.
    ///
    /// When the peak-group intensity is not positive and finite, both scores
    /// are 0, as in the source.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when `mrmfeature` has no feature for a
    /// transition name (the source dereferences a null pointer or throws
    /// `ElementNotFound`), and the window-integration errors of
    /// [`integrate_window_sequence`].
    pub fn dia_isotope_scores(
        &self,
        transitions: &[LightTransition],
        spectra: &[SpectrumPtr],
        mrmfeature: &dyn MRMFeatureAccess,
        im_range: &RangeBase,
    ) -> Result<IsotopeScores> {
        let mut scores = IsotopeScores::default();
        let feature_intensity = f64::from(mrmfeature.intensity());
        if !feature_intensity.is_finite() || feature_intensity <= 0.0 {
            return Ok(scores);
        }
        for transition in transitions {
            let feature = mrmfeature
                .feature(&transition.transition_name)
                .ok_or_else(|| {
                    Error::MissingInformation(format!(
                        "no feature for transition '{}'",
                        transition.transition_name
                    ))
                })?;
            let rel_intensity = f64::from(feature.intensity()) / feature_intensity;
            let charge = if transition.fragment_charge != 0 {
                i32::from(transition.fragment_charge)
            } else {
                1
            };
            let isotopes =
                self.isotope_intensities(transition.product_mz, spectra, charge, im_range)?;
            let score =
                self.score_isotope_pattern_averagine(&isotopes, transition.product_mz, charge)?;
            scores.isotope_corr += score * rel_intensity;
            let (nr_occurrences, _) = self.large_peaks_before_first_isotope(
                spectra,
                transition.product_mz,
                isotopes[0],
                im_range,
            )?;
            scores.isotope_overlap += f64::from(nr_occurrences) * rel_intensity;
        }
        Ok(scores)
    }

    /// Fragment-ion mass-error scores in ppm against each transition's
    /// theoretical product m/z (source `dia_massdiff_score`).
    ///
    /// Each transition's window (width `dia_extraction_window` in the unit of
    /// `dia_extraction_unit`) is integrated and the signed ppm error of the
    /// intensity-weighted m/z recorded. Transitions without signal contribute
    /// `-1` to [`MassDiffScores::diff_ppm`] and are excluded from both scores.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `normalized_library_intensity` is shorter
    /// than `transitions` (the source reads out of bounds) or a product m/z is
    /// zero, and the window-integration errors.
    pub fn dia_massdiff_score(
        &self,
        transitions: &[LightTransition],
        spectra: &[SpectrumPtr],
        normalized_library_intensity: &[f64],
        im_range: &RangeBase,
    ) -> Result<MassDiffScores> {
        if normalized_library_intensity.len() < transitions.len() {
            return Err(bad(
                "normalized_library_intensity must have one entry per transition",
            ));
        }
        let mut out = MassDiffScores::default();
        let mut n_observed = 0usize;
        for (transition, &weight) in transitions.iter().zip(normalized_library_intensity) {
            let Some(found) = self.integrate_at(spectra, transition.product_mz, im_range)? else {
                out.diff_ppm.push(-1.0);
                continue;
            };
            n_observed += 1;
            let diff = ppm(found.mz, transition.product_mz)?;
            out.diff_ppm.push(diff);
            out.ppm_score += diff.abs();
            out.ppm_score_weighted += diff.abs() * weight;
        }
        out.ppm_score = if n_observed > 0 {
            out.ppm_score / n_observed as f64
        } else {
            0.0
        };
        Ok(out)
    }

    /// Precursor mass-error score on an MS1 spectrum (source
    /// `dia_ms1_massdiff_score`).
    ///
    /// Integrates a window around `precursor_mz` and reports the absolute ppm
    /// difference of the measured m/z. Without signal, the score is the
    /// absolute ppm width of the extraction window, a worst-case bound, and
    /// [`Ms1MassDiffScore::signal_found`] is `false`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a zero or non-finite `precursor_mz`, and the
    /// window-integration errors.
    pub fn dia_ms1_massdiff_score(
        &self,
        precursor_mz: f64,
        spectra: &[SpectrumPtr],
        im_range: &RangeBase,
    ) -> Result<Ms1MassDiffScore> {
        let range = create_mz_range_ppm(
            precursor_mz,
            self.dia_extract_window,
            self.dia_extraction_ppm,
        )?;
        match integrate_window_sequence(spectra, &range, im_range, self.dia_centroided)? {
            Some(found) => Ok(Ms1MassDiffScore {
                ppm_score: ppm_abs(found.mz, precursor_mz)?,
                signal_found: true,
            }),
            None => {
                let span = range.span().unwrap_or(0.0);
                Ok(Ms1MassDiffScore {
                    ppm_score: ppm_abs(precursor_mz + span, precursor_mz)?,
                    signal_found: false,
                })
            }
        }
    }

    /// Precursor isotope-pattern scores using an averagine model (source
    /// `dia_ms1_isotope_scores_averagine`).
    ///
    /// The isotopes at spacing `1.0033548 / |charge_state|` are integrated and
    /// correlated with the averagine distribution of weight
    /// `|precursor_mz * charge_state|` (a rough neutral-mass estimate, since the
    /// charge carrier of a negative ion is unknown). The overlap is the largest
    /// intensity ratio of a peak one isotope spacing below `precursor_mz` to the
    /// monoisotope, over charges `1..=dia_nr_charges`.
    ///
    /// The source takes `im_range` by non-`const` reference but does not change
    /// it.
    ///
    /// # Errors
    ///
    /// The window-integration and averagine-model errors, and
    /// [`Error::InvalidRange`] when the averagine model yields a different
    /// number of isotopes than were integrated (as the source's Pearson
    /// correlation throws).
    pub fn dia_ms1_isotope_scores_averagine(
        &self,
        precursor_mz: f64,
        spectra: &[SpectrumPtr],
        charge_state: i32,
        im_range: &RangeBase,
    ) -> Result<IsotopeScores> {
        let isotopes = self.isotope_intensities(precursor_mz, spectra, charge_state, im_range)?;
        let distribution = averagine_distribution(
            (precursor_mz * f64::from(charge_state)).abs(),
            self.dia_nr_isotopes as usize + 1,
        )?;
        let isotope_corr = self.score_isotope_pattern(&isotopes, &distribution)?;
        let (_, max_ratio) =
            self.large_peaks_before_first_isotope(spectra, precursor_mz, isotopes[0], im_range)?;
        Ok(IsotopeScores {
            isotope_corr,
            isotope_overlap: max_ratio,
        })
    }

    /// Precursor isotope-pattern scores using an explicit sum formula (source
    /// `dia_ms1_isotope_scores`).
    ///
    /// The theoretical pattern comes from `sum_formula` (a coarse distribution
    /// of `dia_nr_isotopes + 1` peaks); the isotope spacing uses the formula's
    /// charge, and the pattern is placed at `precursor_mz`, which may differ
    /// slightly from the formula's own m/z (for example with adducts). The
    /// overlap is as in [`DiaScoring::dia_ms1_isotope_scores_averagine`].
    ///
    /// A formula of charge 0 puts every isotope window at a non-finite m/z (the
    /// source divides by zero); those windows find no signal, as in the source.
    ///
    /// # Errors
    ///
    /// As [`DiaScoring::dia_ms1_isotope_scores_averagine`], and the isotope
    /// generator's errors for `sum_formula` (for example a negative charge).
    pub fn dia_ms1_isotope_scores(
        &self,
        precursor_mz: f64,
        spectra: &[SpectrumPtr],
        im_range: &RangeBase,
        sum_formula: &EmpiricalFormula,
    ) -> Result<IsotopeScores> {
        let isotopes =
            self.isotope_intensities(precursor_mz, spectra, sum_formula.charge(), im_range)?;
        let generator = CoarseIsotopePatternGenerator::new(
            Some(self.dia_nr_isotopes as usize + 1),
            CoarseMassMode::default(),
        )?
        .with_precision(ProbabilityPrecision::SourceSingle);
        let distribution = generator.run(sum_formula)?;
        let isotope_corr = self.score_isotope_pattern(&isotopes, &distribution)?;
        let (_, max_ratio) =
            self.large_peaks_before_first_isotope(spectra, precursor_mz, isotopes[0], im_range)?;
        Ok(IsotopeScores {
            isotope_corr,
            isotope_overlap: max_ratio,
        })
    }

    /// Count the b- and y-ions of `sequence` at fragment charge `charge` that
    /// have evidence in the spectrum (source `dia_by_ion_score`).
    ///
    /// An ion matches when a signal is found in its window, its ppm error is
    /// below `dia_byseries_ppm_diff` and its intensity above
    /// `dia_byseries_intensity_min`. The source takes `sequence` by non-`const`
    /// reference but does not change it.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a charge outside `1..=255` (the source checks
    /// `charge > 0` only in debug builds), the errors of
    /// [`by_series`](super::dia_helper::by_series) and the window-integration
    /// errors.
    pub fn dia_by_ion_score(
        &self,
        spectra: &[SpectrumPtr],
        sequence: &AASequence,
        charge: i32,
        im_range: &RangeBase,
    ) -> Result<ByIonScores> {
        let (bseries, yseries) = by_series(sequence, &self.generator, charge)?;
        let mut scores = ByIonScores::default();
        for (series, score) in [
            (&bseries, &mut scores.bseries_score),
            (&yseries, &mut scores.yseries_score),
        ] {
            for &ion_mz in series {
                if let Some(found) = self.integrate_at(spectra, ion_mz, im_range)? {
                    let diff = ppm_abs(found.mz, ion_mz)?;
                    if diff < self.dia_byseries_ppm_diff
                        && found.intensity > self.dia_byseries_intensity_min
                    {
                        *score += 1.0;
                    }
                }
            }
        }
        Ok(scores)
    }

    /// Dot-product and Manhattan scores between the spectrum and the
    /// transition group's theoretical isotope spectrum (source
    /// `score_with_isotopes`).
    ///
    /// Delegates to [`DiaPrescore`] built with `dia_extraction_window`,
    /// `dia_nr_isotopes` and `dia_nr_charges`. The window is passed as Th even
    /// when `dia_extraction_unit` is `ppm`, as in the source.
    ///
    /// # Errors
    ///
    /// As [`DiaPrescore::build_theoretical_spectrum`] and
    /// [`DiaPrescore::score_prepared`].
    pub fn score_with_isotopes(
        &self,
        spectra: &[SpectrumPtr],
        transitions: &[LightTransition],
        im_range: &RangeBase,
    ) -> Result<PrescoreResult> {
        let prescore = DiaPrescore::with_window(
            self.dia_extract_window,
            i32::try_from(self.dia_nr_isotopes).map_err(|_| bad("dia_nr_isotopes too large"))?,
            i32::try_from(self.dia_nr_charges).map_err(|_| bad("dia_nr_charges too large"))?,
        );
        let cache = prescore.build_theoretical_spectrum(transitions)?;
        DiaPrescore::score_prepared(spectra, &cache, self.dia_extract_window, im_range)
    }

    /// Search for an intense peak at lower m/z that could explain `mono_mz` as
    /// part of an isotope pattern (source `largePeaksBeforeFirstIsotope_`).
    ///
    /// For every charge `1..=dia_nr_charges`, the window at
    /// `mono_mz - 1.0033548 / charge` is integrated. Returns `(nr_occurrences,
    /// max_ratio)`: how often a found peak is more intense than `mono_int` and
    /// within `peak_before_mono_max_ppm_diff`, and the largest ratio of a found
    /// peak's intensity to `mono_int` (0 when `mono_int` is 0).
    fn large_peaks_before_first_isotope(
        &self,
        spectra: &[SpectrumPtr],
        mono_mz: f64,
        mono_int: f64,
        im_range: &RangeBase,
    ) -> Result<(u32, f64)> {
        let mut nr_occurrences = 0;
        let mut max_ratio = 0.0;
        for charge in 1..=self.dia_nr_charges {
            let center = mono_mz - DIA_SCORING_C13C12_MASSDIFF_U / f64::from(charge);
            let Some(found) = self.integrate_at(spectra, center, im_range)? else {
                continue;
            };
            let ratio = if mono_int != 0.0 {
                found.intensity / mono_int
            } else {
                0.0
            };
            if ratio > max_ratio {
                max_ratio = ratio;
            }
            let ddiff_ppm = (found.mz - center).abs() * 1e6 / center;
            if ratio > 1.0 && ddiff_ppm < self.peak_before_mono_max_ppm_diff {
                nr_occurrences += 1;
            }
        }
        Ok((nr_occurrences, max_ratio))
    }

    fn score_isotope_pattern_averagine(
        &self,
        isotopes: &[f64],
        product_mz: f64,
        charge: i32,
    ) -> Result<f64> {
        let distribution = averagine_distribution(
            (product_mz * f64::from(charge)).abs(),
            self.dia_nr_isotopes as usize + 1,
        )?;
        self.score_isotope_pattern(isotopes, &distribution)
    }

    /// Pearson correlation of observed isotope intensities with a theoretical
    /// distribution scaled to a maximum of 1, compared by order only (source
    /// `scoreIsotopePattern_(isotopes_int, IsotopeDistribution)`). A NaN
    /// correlation (a constant vector) scores 0.
    fn score_isotope_pattern(
        &self,
        isotopes: &[f64],
        distribution: &IsotopeDistribution,
    ) -> Result<f64> {
        let mut theoretical: Vec<f64> =
            distribution.peaks().iter().map(|p| p.probability).collect();
        let mut max = theoretical
            .iter()
            .fold(0.0_f64, |m, &v| if v > m { v } else { m });
        if max == 0.0 {
            max = 1.0;
        }
        theoretical.iter_mut().for_each(|v| *v /= max);
        let score = pearson_correlation_coefficient(isotopes, &theoretical)?;
        Ok(if score.is_nan() { 0.0 } else { score })
    }
}

/// The averagine distribution of `neutral_mass` limited to `nr_isotopes` peaks,
/// in the source's binary32 intensity precision (source
/// `getCachedAveragineIsotopeDistribution`, without the cache).
fn averagine_distribution(neutral_mass: f64, nr_isotopes: usize) -> Result<IsotopeDistribution> {
    CoarseIsotopePatternGenerator::new(Some(nr_isotopes), CoarseMassMode::default())?
        .with_precision(ProbabilityPrecision::SourceSingle)
        .estimate_from_peptide_weight(neutral_mass)
}
