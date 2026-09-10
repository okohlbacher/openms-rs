// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Precursor isolation purity and SPS fragment matching from OpenMS4-core7c029e8.
//! Scalar purity uses C13/C12 spacing and greedy, non-reusable nearest matches;
//! fuzzy scan estimates deliberately retain their separate source conventions.
//! Validation covers the numerical data used by these algorithms. Unconsumed
//! annotations and acquisition metadata are not traversed or validated; use the
//! kernel's full validation separately when that is required.

use crate::chemistry::{AASequence, C13C12_MASSDIFF_U, TheoreticalSpectrumGenerator};
use crate::comparison::Tolerance;
use crate::kernel::{MSSpectrum, Peak1D, Precursor};
use crate::{Error, Result};

#[path = "precursor_purity_batch.rs"]
mod batch;
#[path = "precursor_purity_fuzzy.rs"]
mod fuzzy;

pub const MAX_PURITY_PEAKS: usize = 1_000_000;
pub const MAX_PURITY_PRECURSORS: usize = 100_000;
pub const MAX_PURITY_ISOTOPES: usize = 1_000_000;
pub const MAX_PURITY_WORK: usize = 50_000_000;

/// Source purity metrics. The residual spectrum contains only unmatched peaks,
/// with default metadata, as in the source's newly constructed isolated window.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PurityScores {
    pub total_intensity: f64,
    pub target_intensity: f64,
    pub signal_proportion: f64,
    pub target_peak_count: usize,
    pub interfering_peak_count: usize,
    pub interfering_peaks: MSSpectrum,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct PrecursorPurity;

impl PrecursorPurity {
    /// Compute scalar isolation purity. Both isolation and tolerance boundaries
    /// are inclusive. Charge sign is ignored and zero means one; tolerance is
    /// doubled, preserving source ppm operation order. Isotopes can match even
    /// when the monoisotopic peak is absent. Empty windows return zero scores.
    /// Sorted, finite, nonnegative input peaks and finite nonnegative precursor
    /// m/z, intensity and isolation offsets are required. Unconsumed annotations
    /// and acquisition metadata are not validated. Inputs are never modified.
    pub fn compute(
        ms1: &MSSpectrum,
        precursor: &Precursor,
        tolerance: Tolerance,
    ) -> Result<PurityScores> {
        compute_with_budget(ms1, precursor, tolerance, &mut WorkBudget::new())
    }

    /// Count SPS precursor m/z windows matching default b/y fragment masses.
    /// Every precursor is counted independently, including duplicate windows.
    /// Charge zero is clamped to one. Source bounds are narrowed to f32 before
    /// matching the compact theoretical mass spectrum; offsets and precursor
    /// charges do not define these matching windows.
    pub fn count_sps_matches(
        precursors: &[Precursor],
        peptide: &AASequence,
        tolerance: Tolerance,
        max_fragment_charge: u8,
    ) -> Result<usize> {
        validate_tolerance(tolerance)?;
        if precursors.len() > MAX_PURITY_PRECURSORS {
            return Err(invalid("SPS precursor count limit exceeded"));
        }
        let mut budget = WorkBudget::new();
        budget.consume(precursors.len())?;
        for precursor in precursors {
            validate_precursor(precursor)?;
        }
        if precursors.is_empty() || peptide.is_empty() {
            return Ok(0);
        }
        let mut theoretical = Vec::new();
        TheoreticalSpectrumGenerator::default().append_mass_spectrum(
            &mut theoretical,
            peptide,
            max_fragment_charge.max(1),
        )?;
        budget.consume(theoretical.len())?;
        let mut matched = 0;
        for precursor in precursors {
            budget.consume(search_work(theoretical.len()))?;
            let width = match tolerance {
                Tolerance::Absolute(value) => value,
                Tolerance::Ppm(value) => ((value / 1e6) * precursor.mz).abs(),
            };
            finite(width, "SPS absolute tolerance")?;
            let lower = finite(precursor.mz - width, "SPS lower bound")? as f32;
            let upper = finite(precursor.mz + width, "SPS upper bound")? as f32;
            if !lower.is_finite() || !upper.is_finite() {
                return Err(invalid("SPS matching bounds exceed finite f32 range"));
            }
            let index = theoretical.partition_point(|&mz| mz < lower);
            if theoretical.get(index).is_some_and(|&mz| mz <= upper) {
                matched += 1;
            }
        }
        Ok(matched)
    }
}

pub(super) fn compute_with_budget(
    ms1: &MSSpectrum,
    precursor: &Precursor,
    tolerance: Tolerance,
    budget: &mut WorkBudget,
) -> Result<PurityScores> {
    validate_tolerance(tolerance)?;
    validate_precursor(precursor)?;
    validate_spectrum(ms1, budget)?;
    let target = precursor.mz;
    let lower = finite(
        target - precursor.isolation_window_lower_offset,
        "isolation lower bound",
    )?;
    let upper = finite(
        target + precursor.isolation_window_upper_offset,
        "isolation upper bound",
    )?;
    let charge = f64::from(precursor.charge.unsigned_abs().max(1));
    let width = finite(
        match tolerance {
            Tolerance::Absolute(value) => value * 2.0,
            Tolerance::Ppm(value) => target * value * 2.0 * 1e-6,
        },
        "doubled precursor tolerance",
    )?;
    budget.consume(search_work(ms1.len()) * 2)?;
    let start = ms1.peaks.partition_point(|peak| peak.mz < lower);
    let end = ms1.peaks.partition_point(|peak| peak.mz <= upper);
    if start == end {
        return Ok(PurityScores::default());
    }

    // Source int(lower_offset*charge), rather than dividing by isotope spacing.
    // Check the defined signed-int conversion range before truncation/use.
    let negative = finite(
        precursor.isolation_window_lower_offset * charge,
        "negative isotope estimate",
    )?
    .trunc();
    if negative > f64::from(i32::MAX) {
        return Err(invalid(
            "negative isotope estimate exceeds source integer range",
        ));
    }
    let mut isotope = -negative;
    let mut expected = finite(
        target + isotope * C13C12_MASSDIFF_U / charge,
        "isotope position",
    )?;
    if expected < lower {
        isotope += 1.0;
        expected = finite(
            target + isotope * C13C12_MASSDIFF_U / charge,
            "isotope position",
        )?;
    }
    let estimated = finite(
        (upper - expected) * charge / C13C12_MASSDIFF_U,
        "isotope iteration estimate",
    )?
    .max(0.0)
    .floor()
        + 2.0;
    if estimated > MAX_PURITY_ISOTOPES as f64 {
        return Err(invalid("precursor isotope iteration limit exceeded"));
    }
    // This conservative preflight covers every expected position. Individual
    // binary searches and physical Vec removals are charged during traversal.
    budget.consume(estimated as usize)?;
    budget.consume((end - start) * 2)?;
    let mut isolated = ms1.peaks[start..end].to_vec();
    let total_intensity = finite(
        isolated
            .iter()
            .fold(0.0, |sum, peak| sum + f64::from(peak.intensity)),
        "isolation total intensity",
    )?;
    let mut target_intensity = 0.0;
    let mut count = 0;
    let mut iterations = 0;
    let mut previous = None;
    loop {
        if expected > upper {
            break;
        }
        if previous.is_some_and(|previous| expected <= previous) {
            return Err(invalid("precursor isotope positions do not advance"));
        }
        iterations += 1;
        if iterations > MAX_PURITY_ISOTOPES {
            return Err(invalid("precursor isotope iteration limit exceeded"));
        }
        budget.consume(1 + search_work(isolated.len()))?;
        let left = finite(expected - width, "isotope matching lower bound")?;
        let right = finite(expected + width, "isotope matching upper bound")?;
        if let Some(index) = nearest(&isolated, expected) {
            if isolated[index].mz >= left && isolated[index].mz <= right {
                budget.consume(isolated.len() - index)?;
                target_intensity = finite(
                    target_intensity + f64::from(isolated[index].intensity),
                    "isolation target intensity",
                )?;
                isolated.remove(index);
                count += 1;
            }
        }
        previous = Some(expected);
        isotope += 1.0;
        expected = finite(
            target + isotope * C13C12_MASSDIFF_U / charge,
            "isotope position",
        )?;
    }
    let signal_proportion = if target_intensity > 0.0 {
        finite(
            target_intensity / total_intensity,
            "precursor signal proportion",
        )?
    } else {
        0.0
    };
    Ok(PurityScores {
        total_intensity,
        target_intensity,
        signal_proportion,
        target_peak_count: count,
        interfering_peak_count: isolated.len(),
        interfering_peaks: MSSpectrum {
            peaks: isolated,
            ..MSSpectrum::default()
        },
    })
}

fn nearest(peaks: &[Peak1D], mz: f64) -> Option<usize> {
    if peaks.is_empty() {
        return None;
    }
    let right = peaks.partition_point(|peak| peak.mz < mz);
    if right == 0 {
        Some(0)
    } else if right == peaks.len() {
        Some(right - 1)
    } else if (peaks[right].mz - mz).abs() < (peaks[right - 1].mz - mz).abs() {
        Some(right)
    } else {
        Some(right - 1)
    }
}
fn search_work(length: usize) -> usize {
    (usize::BITS - length.leading_zeros()) as usize + 1
}

pub(super) struct WorkBudget {
    remaining: usize,
}
impl WorkBudget {
    pub(super) fn new() -> Self {
        Self {
            remaining: MAX_PURITY_WORK,
        }
    }
    pub(super) fn consume(&mut self, units: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(units)
            .ok_or_else(|| invalid("precursor purity work limit exceeded"))?;
        Ok(())
    }
}
pub(super) fn validate_spectrum(spectrum: &MSSpectrum, budget: &mut WorkBudget) -> Result<()> {
    if spectrum.len() > MAX_PURITY_PEAKS || spectrum.precursors.len() > MAX_PURITY_PRECURSORS {
        return Err(invalid("precursor purity input size limit exceeded"));
    }
    budget.consume(spectrum.len().saturating_mul(3).saturating_add(1))?;
    finite(spectrum.rt, "spectrum retention time")?;
    if spectrum.ms_level == 0 {
        return Err(invalid("MS level must be positive"));
    }
    let mut previous = None;
    for peak in &spectrum.peaks {
        if !peak.mz.is_finite()
            || peak.mz < 0.0
            || !peak.intensity.is_finite()
            || peak.intensity < 0.0
        {
            return Err(invalid("purity peaks must be finite and nonnegative"));
        }
        if previous.is_some_and(|mz| peak.mz < mz) {
            return Err(Error::UnsortedData);
        }
        previous = Some(peak.mz);
    }
    Ok(())
}
pub(super) fn validate_precursor(precursor: &Precursor) -> Result<()> {
    // Constant-size numerical validation avoids repeatedly traversing unused
    // metadata when thousands of MS2 scans refer to one annotated MS1 scan.
    for value in [
        precursor.mz,
        f64::from(precursor.intensity),
        precursor.isolation_window_lower_offset,
        precursor.isolation_window_upper_offset,
        precursor.isolation_target_mz.unwrap_or(precursor.mz),
    ] {
        if !value.is_finite() || value < 0.0 {
            return Err(invalid(
                "purity precursor m/z, intensity and isolation values must be finite and nonnegative",
            ));
        }
    }
    Ok(())
}
pub(super) fn validate_tolerance(tolerance: Tolerance) -> Result<()> {
    let value = match tolerance {
        Tolerance::Absolute(value) | Tolerance::Ppm(value) => value,
    };
    if !value.is_finite() || value < 0.0 {
        return Err(invalid("purity tolerance must be finite and nonnegative"));
    }
    Ok(())
}
pub(super) fn finite(value: f64, label: &str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(format!("{label} must be finite")))
    }
}
pub(super) fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}
