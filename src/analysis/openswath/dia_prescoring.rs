// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Scoring of a spectrum given the library intensities of a transition group,
//! from `ANALYSIS/OPENSWATH/DIAPrescoring.h` (class `DiaPrescore`).
//!
//! In DIA (SWATH) analysis a full MS2 spectrum is recorded at every
//! chromatographic point. [`DiaPrescore`] simulates a theoretical spectrum from
//! the library intensities of a transition group (an averagine isotope pattern
//! per transition, plus two pre-isotope positions per transition) and computes
//! a Manhattan distance and a rescaled dot product between the observed
//! spectrum intensities and the simulated spectrum. See also
//! [`DiaScoring`](crate::analysis::openswath::dia_scoring::DiaScoring).
//!
//! The source keeps per-thread scratch vectors to avoid allocation; this port
//! allocates per call, which produces the same values.
//!
//! See `docs/OPENSWATH_SCORING_SUPPORT.md`.

use std::collections::BTreeMap;

use super::data_frame_writer::DataFrameWriter;
use super::dia_helper::{
    add_preisotope_weights_single, add_single_peak_isotopes_to_spec, extract_first, extract_second,
    integrate_windows_sequence, sort_by_first,
};
use super::spectrum_access::{SpectrumAccess, SpectrumPtr};
use super::stats_helpers::{dot_prod, manhattan_dist, norm};
use super::transition_experiment::{LightTargetedExperiment, LightTransition};
use super::transition_helper;
use crate::chemistry::C13C12_MASSDIFF_U;
use crate::kernel::ranges::RangeBase;
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::{Error, Result};

/// The `DefaultParamHandler` name of [`DiaPrescore`] (`"DIAPrescore"`).
pub const DIA_PRESCORE_NAME: &str = "DIAPrescore";

/// Number of pre-isotope positions added per transition (the source's
/// hard-coded `nrNegPeaks`).
const NR_NEG_PEAKS: u32 = 2;

/// The cached theoretical spectrum of one transition group (source
/// `DiaPrescore::TransitionGroupTheoreticalSpectrumCache`).
///
/// Stores precomputed theoretical m/z and intensity arrays plus rescaling
/// bounds for repeated scoring of the same transition group.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TransitionGroupTheoreticalSpectrumCache {
    /// Theoretical fragment/isotope m/z positions used for window integration.
    pub mz_theor: Vec<f64>,
    /// Non-negative normalised theoretical intensities used for the Manhattan
    /// comparison.
    pub int_theor: Vec<f64>,
    /// Signed normalised theoretical intensities including the negatively
    /// weighted pre-isotope peaks.
    pub int_theor_neg: Vec<f64>,
    /// Lower bound used to rescale the signed dot-product score.
    pub neg_val: f64,
    /// Upper bound used to rescale the signed dot-product score.
    pub pos_val: f64,
}

impl TransitionGroupTheoreticalSpectrumCache {
    /// Reset all cached spectrum data: clear the arrays and zero both scaling
    /// values (source `clear`).
    pub fn clear(&mut self) {
        *self = Self::default();
    }
}

/// The two prescores of [`DiaPrescore::score`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrescoreResult {
    /// Dot product against the signed theoretical spectrum, rescaled to
    /// `[-1, 1]` (higher is better).
    pub dotprod: f64,
    /// Manhattan distance to the non-negative theoretical spectrum (lower is
    /// better).
    pub manhattan: f64,
}

/// DIA prescoring of a spectrum against a transition group (source
/// `DiaPrescore`).
///
/// # Parameters
///
/// | name | default | meaning |
/// |---|---|---|
/// | `dia_extraction_window` | 0.1 | DIA extraction window in Th (min 0) |
/// | `nr_isotopes` | 4 | number of isotopes |
/// | `nr_charges` | 4 | number of charges; stored but never used, as in the source |
///
/// The source's `DiaPrescore(dia_extract_window, nr_isotopes, nr_charges)`
/// constructor sets the members without defining the parameter defaults, so
/// its `getParameters()` is empty; [`DiaPrescore::with_window`] builds the same
/// members, and [`DiaPrescore::to_param`] always reports them.
#[derive(Clone, Debug, PartialEq)]
pub struct DiaPrescore {
    /// DIA extraction window in Th.
    pub dia_extract_window: f64,
    /// Number of isotopes per transition in the theoretical spectrum; 0 keeps
    /// every isotope of the averagine model.
    pub nr_isotopes: i32,
    /// Number of charges (unused by the scoring, as in the source).
    pub nr_charges: i32,
}

impl Default for DiaPrescore {
    fn default() -> Self {
        Self {
            dia_extract_window: 0.1,
            nr_isotopes: 4,
            nr_charges: 4,
        }
    }
}

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

fn entry(name: &str, value: ParamValue, description: &str) -> ParamEntry {
    ParamEntry {
        name: name.into(),
        description: description.into(),
        value,
        ..ParamEntry::default()
    }
}

/// Divide `values` by `factor` in place when `factor > 0`, else leave them;
/// the source's `OpenSwath::normalize(v, factor, v)`.
fn normalize_in_place(values: &mut [f64], factor: f64) {
    if factor > 0.0 {
        for value in values {
            *value /= factor;
        }
    }
}

impl DiaPrescore {
    /// A prescorer with the source defaults (source `DiaPrescore()`).
    pub fn new() -> Self {
        Self::default()
    }

    /// A prescorer with the given members (source
    /// `DiaPrescore(double dia_extract_window, int nr_isotopes = 4, int
    /// nr_charges = 4)`).
    pub fn with_window(dia_extract_window: f64, nr_isotopes: i32, nr_charges: i32) -> Self {
        Self {
            dia_extract_window,
            nr_isotopes,
            nr_charges,
        }
    }

    /// The source defaults as a [`Param`] (source `defineDefaults`).
    ///
    /// # Errors
    ///
    /// Propagates a [`Param`] construction failure, which does not occur for
    /// these fixed entries.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The current members as a [`Param`] in the layout of
    /// [`DiaPrescore::defaults`].
    ///
    /// # Errors
    ///
    /// As [`DiaPrescore::defaults`].
    pub fn to_param(&self) -> Result<Param> {
        let mut window = entry(
            "dia_extraction_window",
            ParamValue::Float(self.dia_extract_window),
            "DIA extraction window in Th.",
        );
        window.min_float = 0.0;
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries: vec![
                window,
                entry(
                    "nr_isotopes",
                    ParamValue::Integer(i64::from(self.nr_isotopes)),
                    "nr of istopes",
                ),
                entry(
                    "nr_charges",
                    ParamValue::Integer(i64::from(self.nr_charges)),
                    "nr charges",
                ),
            ],
            nodes: Vec::new(),
        })
    }

    /// Build a prescorer from parameters, as source `setParameters` followed
    /// by `updateMembers_`. Missing parameters take their defaults.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a parameter has the wrong type, violates
    /// its restriction, or an integer does not fit the source `int`.
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`DiaPrescore::from_param`], also returning the unknown-parameter
    /// warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`DiaPrescore::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(DIA_PRESCORE_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, |param| {
            Ok(Self {
                dia_extract_window: param.value("dia_extraction_window")?.to_f64()?,
                nr_isotopes: int32(param, "nr_isotopes")?,
                nr_charges: int32(param, "nr_charges")?,
            })
        })
    }

    /// Score one or more observed spectra against a transition group (source
    /// `score`): build the theoretical spectrum from `transitions` and compare
    /// it with [`DiaPrescore::score_prepared`].
    ///
    /// The input may hold several adjacent MS2 spectra or frames around the
    /// apex, or a single pre-merged spectrum, depending on the upstream merge
    /// mode.
    ///
    /// # Errors
    ///
    /// As [`DiaPrescore::build_theoretical_spectrum`] and
    /// [`DiaPrescore::score_prepared`].
    pub fn score(
        &self,
        spectra: &[SpectrumPtr],
        transitions: &[LightTransition],
        im_range: &RangeBase,
    ) -> Result<PrescoreResult> {
        let cache = self.build_theoretical_spectrum(transitions)?;
        Self::score_prepared(spectra, &cache, self.dia_extract_window, im_range)
    }

    /// Build the theoretical DIA isotope spectrum of a transition group (source
    /// `buildTheoreticalSpectrum`).
    ///
    /// For each transition (fragment charge 1 when unset), the averagine
    /// pattern with [`DiaPrescore::nr_isotopes`] peaks scaled by the library
    /// intensity is added; then two pre-isotope positions at
    /// `product_mz - k * 1.0033548 / |charge|` are added, with weight 0 to the
    /// Manhattan spectrum and with weight `-0.5 * mean library intensity` to the
    /// dot-product spectrum. Both are sorted by m/z. The Manhattan intensities
    /// are square-rooted and normalised to sum 1; the signed intensities are
    /// normalised to unit Euclidean norm. `neg_val` and `pos_val` are the worst
    /// and best achievable dot products, used for rescaling.
    ///
    /// The source comments that the signed spectrum was built from the
    /// zero-weight one in earlier versions ("Was this a bug?"); the current
    /// behaviour is ported.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a negative [`DiaPrescore::nr_isotopes`] (the
    /// source converts it to a huge unsigned count), and the averagine model's
    /// errors (a negative or non-finite product m/z).
    pub fn build_theoretical_spectrum(
        &self,
        transitions: &[LightTransition],
    ) -> Result<TransitionGroupTheoreticalSpectrumCache> {
        let nr_isotopes = usize::try_from(self.nr_isotopes)
            .map_err(|_| bad("nr_isotopes must not be negative"))?;
        let mut cache = TransitionGroupTheoreticalSpectrumCache::default();
        let mut with_iso = Vec::with_capacity(transitions.len() * (nr_isotopes + 2));
        for transition in transitions {
            add_single_peak_isotopes_to_spec(
                transition.product_mz,
                transition.library_intensity,
                &mut with_iso,
                nr_isotopes,
                charge_of(transition),
            )?;
        }
        let mut with_neg = with_iso.clone();
        let mut avg_int = transitions
            .iter()
            .fold(0.0, |acc, t| acc + t.library_intensity);
        avg_int /= transitions.len() as f64;
        let neg_weight = 0.5 * avg_int;
        for transition in transitions {
            let charge = charge_of(transition);
            add_preisotope_weights_single(
                transition.product_mz,
                &mut with_iso,
                NR_NEG_PEAKS,
                0.0,
                C13C12_MASSDIFF_U,
                charge,
            )?;
            add_preisotope_weights_single(
                transition.product_mz,
                &mut with_neg,
                NR_NEG_PEAKS,
                -neg_weight,
                C13C12_MASSDIFF_U,
                charge,
            )?;
        }
        sort_by_first(&mut with_iso);
        sort_by_first(&mut with_neg);

        extract_first(&with_iso, &mut cache.mz_theor);
        extract_second(&with_iso, &mut cache.int_theor);
        cache.int_theor.iter_mut().for_each(|v| *v = v.sqrt());
        let total: f64 = cache.int_theor.iter().fold(0.0, |acc, v| acc + v);
        normalize_in_place(&mut cache.int_theor, total);

        extract_second(&with_neg, &mut cache.int_theor_neg);
        let neg_norm = norm(&cache.int_theor_neg);
        normalize_in_place(&mut cache.int_theor_neg, neg_norm);

        let n_neg = f64::from(NR_NEG_PEAKS) * transitions.len() as f64;
        cache.neg_val = (-neg_weight / neg_norm) * n_neg.sqrt();
        let mut best: Vec<f64> = cache
            .int_theor_neg
            .iter()
            .map(|&v| {
                if v > 0.0 {
                    v * n_neg * neg_weight / neg_norm
                } else {
                    0.0
                }
            })
            .collect();
        let best_norm = norm(&best);
        normalize_in_place(&mut best, best_norm);
        cache.pos_val = dot_prod(&best, &cache.int_theor_neg)?;
        Ok(cache)
    }

    /// Score an observed spectrum sequence against a precomputed theoretical
    /// spectrum (source static `scorePrepared`).
    ///
    /// The observed intensities in windows of width `dia_extract_window` Th
    /// around every theoretical m/z are square-rooted and normalised to sum 1;
    /// the Manhattan distance to `int_theor` is reported. They are then
    /// normalised to unit Euclidean norm, dotted with `int_theor_neg`, and the
    /// result rescaled linearly from `[neg_val, pos_val]` to `[-1, 1]`.
    ///
    /// The source's TODO notes that a combined Manhattan distance penalises a
    /// missing transition twice through the normalisation.
    ///
    /// # Errors
    ///
    /// As [`integrate_windows_sequence`](super::dia_helper::integrate_windows_sequence):
    /// in particular [`Error::MissingInformation`] for an empty theoretical
    /// spectrum. [`Error::InvalidValue`] when `int_theor` or `int_theor_neg` is
    /// shorter than `mz_theor`.
    pub fn score_prepared(
        spectra: &[SpectrumPtr],
        cache: &TransitionGroupTheoreticalSpectrumCache,
        dia_extract_window: f64,
        im_range: &RangeBase,
    ) -> Result<PrescoreResult> {
        let windows = integrate_windows_sequence(
            spectra,
            &cache.mz_theor,
            dia_extract_window,
            im_range,
            false,
        )?;
        let mut int_exp = windows.intensity;
        int_exp.iter_mut().for_each(|v| *v = v.sqrt());
        let total: f64 = int_exp.iter().fold(0.0, |acc, v| acc + v);
        normalize_in_place(&mut int_exp, total);
        let manhattan = manhattan_dist(&int_exp, &cache.int_theor)?;
        let euclid = norm(&int_exp);
        normalize_in_place(&mut int_exp, euclid);
        let raw = dot_prod(&int_exp, &cache.int_theor_neg)?;
        let dotprod = ((raw - cache.neg_val) * 2.0) / (cache.pos_val - cache.neg_val) - 1.0;
        Ok(PrescoreResult { dotprod, manhattan })
    }

    /// Prescore every transition group of `transition_exp` against every
    /// spectrum of `swath` and write the score vectors to `writer` (source
    /// `operator()`).
    ///
    /// The transitions are grouped by compound reference (in ascending order,
    /// as the source's `std::map`), which become the column names. For each
    /// spectrum `i`, two rows are stored: `score1_<RT>` with the dot products
    /// and `score2_<RT>` with the Manhattan distances, where `<RT>` is the
    /// retention time as `boost::lexical_cast<std::string>` prints a `double`
    /// (17 significant digits, `%.17g`). The source also prints a progress
    /// line per spectrum to standard output; this port does not.
    ///
    /// # Errors
    ///
    /// The errors of the spectrum access, of [`DiaPrescore::score`] and of the
    /// writer.
    pub fn score_all(
        &self,
        swath: &dyn SpectrumAccess,
        transition_exp: &LightTargetedExperiment,
        im_range: &RangeBase,
        writer: &mut dyn DataFrameWriter,
    ) -> Result<()> {
        let mut groups: BTreeMap<String, Vec<LightTransition>> = BTreeMap::new();
        transition_helper::convert(transition_exp, &mut groups);
        let names: Vec<String> = groups.keys().cloned().collect();
        writer.colnames(&names)?;
        for i in 0..swath.nr_spectra() {
            let spectrum = swath.spectrum_by_id(i)?;
            let meta = swath.spectrum_meta_by_id(i)?;
            let sequence = vec![spectrum];
            let mut dot = Vec::with_capacity(groups.len());
            let mut man = Vec::with_capacity(groups.len());
            for transitions in groups.values() {
                let result = self.score(&sequence, transitions, im_range)?;
                dot.push(result.dotprod);
                man.push(result.manhattan);
            }
            let rt = format_g17(meta.rt);
            writer.store(&format!("score1_{rt}"), &dot)?;
            writer.store(&format!("score2_{rt}"), &man)?;
        }
        Ok(())
    }
}

fn charge_of(transition: &LightTransition) -> i32 {
    if transition.fragment_charge != 0 {
        i32::from(transition.fragment_charge)
    } else {
        1
    }
}

fn int32(param: &Param, key: &str) -> Result<i32> {
    i32::try_from(param.value(key)?.to_i64()?)
        .map_err(|_| bad(format!("'{key}' does not fit the source 32-bit parameter")))
}

/// `boost::lexical_cast<std::string>(double)`: `%.17g` with trailing zeros
/// removed and the C spellings of the non-finite values.
pub(crate) fn format_g17(value: f64) -> String {
    const DIGITS: i32 = 17;
    if value.is_nan() {
        return if value.is_sign_negative() {
            "-nan"
        } else {
            "nan"
        }
        .into();
    }
    if value.is_infinite() {
        return if value < 0.0 { "-inf" } else { "inf" }.into();
    }
    if value == 0.0 {
        return if value.is_sign_negative() { "-0" } else { "0" }.into();
    }
    let scientific = format!("{:.*e}", (DIGITS - 1) as usize, value);
    let Some((mantissa, exponent)) = scientific.split_once('e') else {
        return scientific;
    };
    let exponent: i32 = exponent.parse().unwrap_or(0);
    if !(-4..DIGITS).contains(&exponent) {
        let sign = if exponent < 0 { '-' } else { '+' };
        return format!("{}e{sign}{:02}", trim_zeros(mantissa), exponent.abs());
    }
    let decimals = usize::try_from(DIGITS - 1 - exponent).unwrap_or(0);
    trim_zeros(&format!("{value:.decimals$}")).to_string()
}

fn trim_zeros(text: &str) -> &str {
    if text.contains('.') {
        text.trim_end_matches('0').trim_end_matches('.')
    } else {
        text
    }
}

#[cfg(test)]
mod tests {
    use super::format_g17;

    #[test]
    fn g17_matches_the_c_rule() {
        // Independently derived from the C `%.17g` rule.
        assert_eq!(format_g17(10.5), "10.5");
        assert_eq!(format_g17(0.1), "0.10000000000000001");
        assert_eq!(format_g17(1e20), "1e+20");
        assert_eq!(format_g17(0.0), "0");
        assert_eq!(format_g17(1234.0), "1234");
    }
}
