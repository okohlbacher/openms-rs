// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Helper functions for the DIA scoring of OpenSWATH, from
//! `ANALYSIS/OPENSWATH/DIAHelper.h` (namespace `OpenMS::DIAHelpers`).
//!
//! Window integration over light spectra ([`integrate_window`],
//! [`integrate_windows`] and their spectrum-sequence forms), extraction-window
//! arithmetic ([`adjust_extraction_window`], [`create_mz_range_ppm`]),
//! theoretical b/y masses ([`by_series`], [`theoretical_masses`]) and the
//! averagine isotope helpers used by [`DiaPrescore`](crate::analysis::openswath::dia_prescoring::DiaPrescore)
//! and [`DiaScoring`](crate::analysis::openswath::dia_scoring::DiaScoring).
//!
//! The source's out-parameters become return values: a window integration
//! returns `Ok(None)` where the source returns `false` and sets `mz = -1`,
//! `im = -1`, `intensity = 0`. The source memoises averagine distributions in a
//! thread-local hash map of up to 100000 entries; the port recomputes them, which
//! yields the same values without hidden global state.
//!
//! See `docs/OPENSWATH_SCORING_SUPPORT.md`.

use super::data_structures::Spectrum;
use super::spectrum_access::SpectrumPtr;
use crate::chemistry::isotopes::ProbabilityPrecision;
use crate::chemistry::{
    AASequence, CoarseIsotopePatternGenerator, CoarseMassMode, TheoreticalSpectrumGenerator,
};
use crate::concept::math_functions::ppm_to_mass;
use crate::kernel::ranges::RangeBase;
use crate::{Error, Result};

/// The averagine isotope spacing the source uses as its default `mannmass`,
/// in Th at charge 1 (`DEFAULT_MANN_MASS` in `DIAHelper.cpp`).
pub const DEFAULT_MANN_MASS: f64 = 1.00048;

/// The default `mannmass` of the two `addPreisotopeWeights` overloads.
///
/// The source header spells this default `1.000482`, which differs in the
/// sixth decimal from [`DEFAULT_MANN_MASS`]; both values are kept.
pub const DEFAULT_PREISOTOPE_MANN_MASS: f64 = 1.000482;

/// One integrated window: intensity-weighted m/z and ion mobility plus the
/// summed intensity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WindowIntegral {
    /// Intensity-weighted mean m/z of the signal in the window.
    pub mz: f64,
    /// Intensity-weighted mean ion mobility, or `-1` when no ion-mobility range
    /// was given (the source sentinel, kept because the scores consume it).
    pub im: f64,
    /// Summed intensity of the signal in the window; always `> 0`.
    pub intensity: f64,
}

/// The unnormalised sums [`integrate_window_sums`] accumulates.
///
/// Source `integrateWindow_` accumulates into three caller-initialised
/// doubles; this struct is those doubles.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct WindowSums {
    /// Sum of `intensity * mz`.
    pub mz: f64,
    /// Sum of `intensity * im` (only accumulated with a non-empty IM range).
    pub im: f64,
    /// Sum of intensities.
    pub intensity: f64,
}

/// The per-window results of [`integrate_windows`].
///
/// The three vectors are parallel, as the source's three out-vectors are.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IntegratedWindows {
    /// Integrated intensity per window.
    pub intensity: Vec<f64>,
    /// Intensity-weighted m/z per window; the window center for an empty window.
    pub mz: Vec<f64>,
    /// Intensity-weighted ion mobility per window; the IM range center for an
    /// empty window with an IM range, `-1` without one.
    pub im: Vec<f64>,
}

impl IntegratedWindows {
    /// Number of windows reported.
    pub fn len(&self) -> usize {
        self.intensity.len()
    }

    /// Whether no window was reported.
    pub fn is_empty(&self) -> bool {
        self.intensity.is_empty()
    }
}

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

/// Accumulate the raw sums of one spectrum in `mz_range` (and `im_range`, when
/// not empty) into `sums` (source `integrateWindow_`).
///
/// Peaks with `mz_range.min <= mz < mz_range.max` contribute: the lower bound is
/// found by binary search, the upper bound is exclusive, exactly as the source.
/// With a non-empty `im_range`, only peaks whose ion mobility lies in the closed
/// range contribute. The sums are not divided; [`integrate_window`] does that.
///
/// An empty spectrum contributes nothing; the source logs a warning.
///
/// # Notes
///
/// The m/z array must be sorted ascending. The source checks this only in debug
/// builds (`OPENMS_PRECONDITION`); this port does not check it either, because
/// the helper runs once per window and an O(n) check would dominate. Unsorted
/// input yields an unspecified but memory-safe result.
///
/// # Errors
///
/// - [`Error::MissingInformation`] when `im_range` is non-empty and the
///   spectrum has no ion-mobility array, as the source throws.
/// - [`Error::Unsupported`] for `centroided == true` on a non-empty spectrum,
///   the source `Exception::NotImplemented`.
/// - [`Error::InvalidValue`] when the intensity array (or, with an IM range, the
///   ion-mobility array) is shorter than the m/z array. The source checks these
///   lengths only in debug builds and otherwise reads out of bounds.
pub fn integrate_window_sums(
    spectrum: &Spectrum,
    mz_range: &RangeBase,
    im_range: &RangeBase,
    centroided: bool,
    sums: &mut WindowSums,
) -> Result<()> {
    let mz = &spectrum.mz_array().data;
    let intensity = &spectrum.intensity_array().data;
    if mz.is_empty() {
        return Ok(());
    }
    let drift = if im_range.is_empty() {
        None
    } else {
        match spectrum.drift_time_array() {
            Some(array) => Some(&array.data),
            None => {
                return Err(Error::MissingInformation(
                    "Cannot integrate with drift time if no drift time is available".into(),
                ));
            }
        }
    };
    if centroided {
        return Err(Error::Unsupported(
            "window integration of centroided DIA data is not implemented".into(),
        ));
    }
    if intensity.len() < mz.len() {
        return Err(bad("MZ and Intensity array need to have the same length"));
    }
    if let Some(drift) = drift {
        if drift.len() < mz.len() {
            return Err(bad("MZ and Drift Time array need to have the same length"));
        }
    }
    let (lo, hi) = range_bounds(mz_range);
    let start = mz.partition_point(|&value| value < lo);
    match drift {
        Some(drift) => {
            for i in start..mz.len() {
                if mz[i] >= hi {
                    break;
                }
                if im_range.contains(drift[i]) {
                    sums.intensity += intensity[i];
                    sums.im += intensity[i] * drift[i];
                    sums.mz += intensity[i] * mz[i];
                }
            }
        }
        None => {
            for i in start..mz.len() {
                if mz[i] >= hi {
                    break;
                }
                sums.intensity += intensity[i];
                sums.mz += intensity[i] * mz[i];
            }
        }
    }
    Ok(())
}

/// `(min, max)` of a range as the source `getMin`/`getMax` read them: an empty
/// range has `min = +max_double` and `max = -max_double`, so nothing lies in it.
fn range_bounds(range: &RangeBase) -> (f64, f64) {
    match (range.min(), range.max()) {
        (Ok(min), Ok(max)) => (min, max),
        _ => (f64::MAX, -f64::MAX),
    }
}

fn finish(sums: WindowSums, im_range: &RangeBase) -> Option<WindowIntegral> {
    if sums.intensity > 0.0 {
        Some(WindowIntegral {
            mz: sums.mz / sums.intensity,
            im: if im_range.is_empty() {
                -1.0
            } else {
                sums.im / sums.intensity
            },
            intensity: sums.intensity,
        })
    } else {
        None
    }
}

/// Integrate the intensity of one spectrum in `mz_range` (and `im_range`, if
/// not empty), returning the intensity-weighted m/z and ion mobility and the
/// total intensity (source `integrateWindow(const SpectrumPtr&, ...)`).
///
/// Returns `Ok(None)` when no signal was found (total intensity not positive);
/// the source returns `false` and sets `mz` and `im` to `-1` and the intensity
/// to `0`. Without an ion-mobility range the returned `im` is `-1`, as in the
/// source.
///
/// # Errors
///
/// As [`integrate_window_sums`].
pub fn integrate_window(
    spectrum: &Spectrum,
    mz_range: &RangeBase,
    im_range: &RangeBase,
    centroided: bool,
) -> Result<Option<WindowIntegral>> {
    let mut sums = WindowSums::default();
    integrate_window_sums(spectrum, mz_range, im_range, centroided, &mut sums)?;
    Ok(finish(sums, im_range))
}

/// Integrate a spectrum sequence in `mz_range` (and `im_range`, if not empty):
/// the sums of all spectra are pooled before dividing (source
/// `integrateWindow(const SpectrumSequence&, ...)`).
///
/// An empty sequence returns `Ok(None)`; the source logs a warning.
///
/// # Errors
///
/// As [`integrate_window_sums`], for any spectrum of the sequence.
pub fn integrate_window_sequence(
    spectra: &[SpectrumPtr],
    mz_range: &RangeBase,
    im_range: &RangeBase,
    centroided: bool,
) -> Result<Option<WindowIntegral>> {
    let mut sums = WindowSums::default();
    for spectrum in spectra {
        integrate_window_sums(spectrum, mz_range, im_range, centroided, &mut sums)?;
    }
    Ok(finish(sums, im_range))
}

fn window_range(center: f64, width: f64) -> Result<RangeBase> {
    let mut range = RangeBase::singular(center)?;
    range.min_span_if_singular(width)?;
    Ok(range)
}

fn push_windows(
    windows_center: &[f64],
    width: f64,
    im_range: &RangeBase,
    remove_zero: bool,
    mut integrate: impl FnMut(&RangeBase) -> Result<Option<WindowIntegral>>,
) -> Result<IntegratedWindows> {
    let mut out = IntegratedWindows::default();
    for &center in windows_center {
        let range = window_range(center, width)?;
        match integrate(&range)? {
            Some(found) => {
                out.intensity.push(found.intensity);
                out.mz.push(found.mz);
                out.im.push(found.im);
            }
            None if !remove_zero => {
                out.intensity.push(0.0);
                out.mz.push(center);
                out.im.push(im_range.center().unwrap_or(-1.0));
            }
            None => {}
        }
    }
    Ok(out)
}

/// Integrate one spectrum in several m/z windows of total width `width` centred
/// at `windows_center`, restricted to `im_range` when it is not empty (source
/// `integrateWindows(const SpectrumPtr&, ...)`).
///
/// A window without signal is dropped when `remove_zero` is set; otherwise it
/// reports intensity `0`, its center as m/z, and the center of `im_range` (or
/// `-1` without one) as ion mobility.
///
/// # Errors
///
/// As [`integrate_window_sums`], and [`Error::InvalidValue`] for a non-finite
/// window center or width (the source builds a NaN range).
pub fn integrate_windows(
    spectrum: &Spectrum,
    windows_center: &[f64],
    width: f64,
    im_range: &RangeBase,
    remove_zero: bool,
) -> Result<IntegratedWindows> {
    push_windows(windows_center, width, im_range, remove_zero, |range| {
        integrate_window(spectrum, range, im_range, false)
    })
}

/// Integrate a spectrum sequence in several m/z windows (source
/// `integrateWindows(const SpectrumSequence&, ...)`); see
/// [`integrate_windows`] for the window rules.
///
/// An empty sequence returns empty results; the source logs a warning.
///
/// # Errors
///
/// [`Error::MissingInformation`] when `windows_center` is empty ("No windows
/// supplied!", as the source throws; the single-spectrum form does not), and
/// the errors of [`integrate_windows`].
pub fn integrate_windows_sequence(
    spectra: &[SpectrumPtr],
    windows_center: &[f64],
    width: f64,
    im_range: &RangeBase,
    remove_zero: bool,
) -> Result<IntegratedWindows> {
    if windows_center.is_empty() {
        return Err(Error::MissingInformation("No windows supplied!".into()));
    }
    if spectra.is_empty() {
        return Ok(IntegratedWindows::default());
    }
    push_windows(windows_center, width, im_range, remove_zero, |range| {
        integrate_window_sequence(spectra, range, im_range, false)
    })
}

/// Widen `[left, right]` by an extraction window of total width
/// `mz_extract_window`, in Th or (with `mz_extraction_ppm`) in ppm of each
/// bound, and return the new `(left, right)` (source `adjustExtractionWindow`,
/// whose out-parameters are ordered `right, left`).
///
/// A left bound that becomes negative is clamped to `0`.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `mz_extract_window` is not positive and finite.
/// The source states the positivity precondition but checks it only in debug
/// builds.
pub fn adjust_extraction_window(
    left: f64,
    right: f64,
    mz_extract_window: f64,
    mz_extraction_ppm: bool,
) -> Result<(f64, f64)> {
    if !(mz_extract_window.is_finite() && mz_extract_window > 0.0) {
        return Err(bad("MZ extraction window needs to be larger than zero."));
    }
    let (mut left, mut right) = (left, right);
    if mz_extraction_ppm {
        left -= left * mz_extract_window / 2e6;
        right += right * mz_extract_window / 2e6;
    } else {
        left -= mz_extract_window / 2.0;
        right += mz_extract_window / 2.0;
    }
    if left < 0.0 {
        left = 0.0;
    }
    Ok((left, right))
}

fn fragment_charge(charge: i32) -> Result<u8> {
    u8::try_from(charge).ok().filter(|&c| c > 0).ok_or_else(|| {
        bad("For constructing b/y series we require charge being a positive integer")
    })
}

/// The theoretical b- and y-ion m/z values of `sequence` at fragment charge
/// `charge`, returned as `(b_series, y_series)` (source `getBYSeries`).
///
/// `generator` produces the spectrum between `charge` and `charge`; its peaks
/// are classified by the first character of their `IonNames` annotation, so the
/// generator must have [`add_metainfo`](TheoreticalSpectrumGenerator::add_metainfo)
/// set (the source's caller sets it too). An empty sequence yields two empty
/// vectors.
///
/// # Errors
///
/// - [`Error::InvalidValue`] for a charge outside `1..=255` (the source checks
///   `charge > 0` only in debug builds) and the generator's own errors.
/// - [`Error::MissingInformation`] when the generated spectrum has no
///   `IonNames` array of matching length; the source reads the first string
///   array unchecked.
pub fn by_series(
    sequence: &AASequence,
    generator: &TheoreticalSpectrumGenerator,
    charge: i32,
) -> Result<(Vec<f64>, Vec<f64>)> {
    let charge = fragment_charge(charge)?;
    let (mut b, mut y) = (Vec::new(), Vec::new());
    if sequence.is_empty() {
        return Ok((b, y));
    }
    let spectrum = generator.generate(sequence, charge, charge, None)?;
    let names = spectrum
        .string_data_arrays
        .first()
        .filter(|names| names.data.len() == spectrum.peaks.len())
        .ok_or_else(|| {
            Error::MissingInformation(
                "the theoretical spectrum has no ion-name annotation (set add_metainfo)".into(),
            )
        })?;
    for (peak, name) in spectrum.peaks.iter().zip(&names.data) {
        match name.as_bytes().first() {
            Some(b'y') => y.push(peak.mz),
            Some(b'b') => b.push(peak.mz),
            _ => {}
        }
    }
    Ok((b, y))
}

/// All theoretical fragment m/z values `generator` produces for `sequence` at
/// fragment charge `charge` (source `getTheorMasses`).
///
/// # Errors
///
/// [`Error::InvalidValue`] for a charge outside `1..=255` (checked by the source
/// only in debug builds) and the generator's own errors.
pub fn theoretical_masses(
    sequence: &AASequence,
    generator: &TheoreticalSpectrumGenerator,
    charge: i32,
) -> Result<Vec<f64>> {
    let charge = fragment_charge(charge)?;
    let spectrum = generator.generate(sequence, charge, charge, None)?;
    Ok(spectrum.peaks.iter().map(|peak| peak.mz).collect())
}

fn nonzero_charge(charge: i32, what: &str) -> Result<f64> {
    if charge == 0 {
        return Err(bad(format!("Charge must be non-zero when {what}.")));
    }
    Ok(f64::from(charge.unsigned_abs()))
}

/// The averagine isotope intensities for a fragment of m/z `product_mz` at
/// `|charge|`, as the source's cache computes them: a coarse generator limited
/// to `nr_isotopes` peaks (0 means unlimited) estimates the peptide pattern of
/// average weight `product_mz * |charge|`, in the source's binary32 intensity
/// precision.
fn averagine_intensities(product_mz: f64, abs_charge: f64, nr_isotopes: usize) -> Result<Vec<f64>> {
    let generator = CoarseIsotopePatternGenerator::new(
        (nr_isotopes > 0).then_some(nr_isotopes),
        CoarseMassMode::default(),
    )?
    .with_precision(ProbabilityPrecision::SourceSingle);
    let distribution = generator.estimate_from_peptide_weight(product_mz * abs_charge)?;
    Ok(distribution
        .peaks()
        .iter()
        .map(|peak| peak.probability)
        .collect())
}

/// The averagine isotope distribution of a fragment at `product_mz`, as
/// `(m/z, relative intensity)` pairs spaced `mannmass / |charge|` apart and
/// starting at `product_mz` (source `getAveragineIsotopeDistribution`, whose
/// defaults are `charge = 1`, `nr_isotopes = 4`, `mannmass = 1.00048`).
///
/// The weight passed to the averagine model is `product_mz * |charge|`; the
/// source notes this is a rough estimate that leaves the charge carriers in,
/// kept for backwards compatibility. `nr_isotopes = 0` keeps every isotope, as
/// the source's generator does.
///
/// # Errors
///
/// [`Error::InvalidValue`] for `charge == 0` (the source throws
/// `IllegalArgument`) and the averagine model's errors (negative or non-finite
/// weight).
pub fn averagine_isotope_distribution(
    product_mz: f64,
    charge: i32,
    nr_isotopes: usize,
    mannmass: f64,
) -> Result<Vec<(f64, f64)>> {
    let abs_charge = nonzero_charge(charge, "generating averagine isotope distributions")?;
    let intensities = averagine_intensities(product_mz, abs_charge, nr_isotopes)?;
    let mut mass = product_mz;
    let mut out = Vec::with_capacity(intensities.len());
    for intensity in intensities {
        out.push((mass, intensity));
        mass += mannmass / abs_charge;
    }
    Ok(out)
}

/// Simulate an isotope-resolved spectrum for `sequence`: the theoretical
/// fragment masses at `charge` and, for each, its averagine distribution with
/// the default four isotopes (source `simulateSpectrumFromAASequence`).
///
/// Returns `(first_isotope_masses, isotope_masses)`.
///
/// # Errors
///
/// As [`theoretical_masses`] and [`averagine_isotope_distribution`].
pub fn simulate_spectrum_from_sequence(
    sequence: &AASequence,
    generator: &TheoreticalSpectrumGenerator,
    charge: i32,
) -> Result<(Vec<f64>, Vec<(f64, f64)>)> {
    let first = theoretical_masses(sequence, generator, charge)?;
    let mut isotopes = Vec::new();
    for &mass in &first {
        isotopes.extend(averagine_isotope_distribution(
            mass,
            charge,
            4,
            DEFAULT_MANN_MASS,
        )?);
    }
    Ok((first, isotopes))
}

/// For every peak of an experimental spectrum `spec` of `(m/z, intensity)`,
/// append its averagine isotope pattern (`nr_isotopes` peaks, intensities
/// scaled by the peak intensity) to `isotope_masses` (source
/// `addIsotopes2Spec`).
///
/// # Errors
///
/// As [`averagine_isotope_distribution`].
pub fn add_isotopes_to_spec(
    spec: &[(f64, f64)],
    isotope_masses: &mut Vec<(f64, f64)>,
    nr_isotopes: usize,
    charge: i32,
) -> Result<()> {
    let mut staged = Vec::new();
    for &(mz, intensity) in spec {
        for (mass, relative) in
            averagine_isotope_distribution(mz, charge, nr_isotopes, DEFAULT_MANN_MASS)?
        {
            staged.push((mass, relative * intensity));
        }
    }
    isotope_masses.extend(staged);
    Ok(())
}

/// Append the averagine isotope pattern of a single peak (`mz`, intensity
/// `intensity`) with `nr_isotopes` peaks spaced [`DEFAULT_MANN_MASS`]` / |charge|`
/// to `isotope_masses` (source `addSinglePeakIsotopes2Spec`).
///
/// # Errors
///
/// As [`averagine_isotope_distribution`]; nothing is appended on error.
pub fn add_single_peak_isotopes_to_spec(
    mz: f64,
    intensity: f64,
    isotope_masses: &mut Vec<(f64, f64)>,
    nr_isotopes: usize,
    charge: i32,
) -> Result<()> {
    let abs_charge = nonzero_charge(charge, "adding isotope peaks")?;
    let intensities = averagine_intensities(mz, abs_charge, nr_isotopes)?;
    let mut mass = mz;
    isotope_masses.reserve(intensities.len());
    for relative in intensities {
        isotope_masses.push((mass, relative * intensity));
        mass += DEFAULT_MANN_MASS / abs_charge;
    }
    Ok(())
}

/// For every mass in `first_isotope_masses`, append `nr_peaks` pre-isotope
/// peaks at `mass - k * mannmass / |charge|` (`k = 1..=nr_peaks`) with weight
/// `pre_isotope_peaks_weight`, then sort `isotope_spec` by m/z (source
/// `addPreisotopeWeights(const std::vector<double>&, ...)`; defaults
/// `nr_peaks = 2`, weight `-0.5`, `mannmass = 1.000482`, `charge = 1`).
///
/// # Errors
///
/// [`Error::InvalidValue`] for `charge == 0`, as the source throws
/// `IllegalArgument`.
pub fn add_preisotope_weights(
    first_isotope_masses: &[f64],
    isotope_spec: &mut Vec<(f64, f64)>,
    nr_peaks: u32,
    pre_isotope_peaks_weight: f64,
    mannmass: f64,
    charge: i32,
) -> Result<()> {
    let abs_charge = nonzero_charge(charge, "adding pre-isotope peaks")?;
    for &mass in first_isotope_masses {
        for k in 1..=nr_peaks {
            isotope_spec.push((
                mass - (f64::from(k) * mannmass) / abs_charge,
                pre_isotope_peaks_weight,
            ));
        }
    }
    sort_by_first(isotope_spec);
    Ok(())
}

/// Append `nr_peaks` pre-isotope peaks before a single `mz`, as
/// [`add_preisotope_weights`] does per mass, but without sorting (source
/// `addPreisotopeWeights(double mz, ...)`).
///
/// # Errors
///
/// [`Error::InvalidValue`] for `charge == 0`. The source's scalar overload does
/// not check and divides by zero, producing infinite masses.
pub fn add_preisotope_weights_single(
    mz: f64,
    isotope_spec: &mut Vec<(f64, f64)>,
    nr_peaks: u32,
    pre_isotope_peaks_weight: f64,
    mannmass: f64,
    charge: i32,
) -> Result<()> {
    let abs_charge = nonzero_charge(charge, "adding pre-isotope peaks")?;
    for k in 1..=nr_peaks {
        isotope_spec.push((
            mz - (f64::from(k) * mannmass) / abs_charge,
            pre_isotope_peaks_weight,
        ));
    }
    Ok(())
}

/// Sort `(m/z, intensity)` pairs by m/z (source `sortByFirst`).
///
/// The source uses the unstable `std::sort`; this sort is stable, which is one
/// of the orders the source may produce for equal m/z values. NaN sorts last.
pub fn sort_by_first(pairs: &mut [(f64, f64)]) {
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
}

/// Append the first element of every pair to `out` (source `extractFirst`).
pub fn extract_first(pairs: &[(f64, f64)], out: &mut Vec<f64>) {
    out.extend(pairs.iter().map(|pair| pair.0));
}

/// Append the second element of every pair to `out` (source `extractSecond`).
pub fn extract_second(pairs: &[(f64, f64)], out: &mut Vec<f64>) {
    out.extend(pairs.iter().map(|pair| pair.1));
}

/// An m/z range centred at `mz_ref` of total width `dia_extraction_window`,
/// given in Th or, with `ppm`, in ppm of `mz_ref` (source `createMZRangePPM`).
///
/// # Errors
///
/// [`Error::InvalidValue`] for a non-finite center or width. A negative width
/// yields an empty range, as in the source.
pub fn create_mz_range_ppm(
    mz_ref: f64,
    dia_extraction_window: f64,
    ppm: bool,
) -> Result<RangeBase> {
    let width = if ppm {
        ppm_to_mass(dia_extraction_window, mz_ref)?
    } else {
        dia_extraction_window
    };
    window_range(mz_ref, width)
}
