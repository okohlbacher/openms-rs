// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! C13-spacing cluster detection and optional conversion to single charge.
//! Native simple and averagine/KL Deisotoper algorithms, with configurable
//! isotope sharing, preserved original indices and checked annotations.

use super::{SpectrumFilter, checked_intensity};
use crate::chemistry::{C13C12_MASSDIFF_U, CoarseIsotopePatternGenerator, PROTON_MASS_U};
use crate::comparison::Tolerance;
use crate::kernel::{DataArray, MSSpectrum, nearest};
use crate::{Error, Result};

/// One accepted isotope ladder, with indices in the original input spectrum.
/// Distinct ladders can share members when allow_shared_isotopes is enabled.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IsotopeCluster {
    pub charge: u8,
    pub peak_indices: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DeisotopingResult {
    pub spectrum: MSSpectrum,
    /// In increasing input monoisotopic m/z order, before charge conversion.
    pub clusters: Vec<IsotopeCluster>,
}

/// Options for OpenMS's simple C13-spacing/decreasing-intensity algorithm.
#[derive(Clone, Debug)]
pub struct Deisotoper {
    pub tolerance: Tolerance,
    pub min_charge: u8,
    pub max_charge: u8,
    pub keep_only_deisotoped: bool,
    pub min_isotope_peaks: usize,
    pub max_isotope_peaks: usize,
    pub make_single_charged: bool,
    pub annotate_charge: bool,
    pub annotate_isotope_peak_count: bool,
    pub annotate_features: bool,
    pub use_decreasing_model: bool,
    /// Zero/one compares M to M+1; two permits M+1 to exceed M.
    pub start_intensity_check: usize,
    pub add_up_intensity: bool,
    /// Permit different accepted ladders to reuse an isotope extension.
    /// Seeds already assigned to a ladder are always skipped.
    pub allow_shared_isotopes: bool,
    /// Bound attempted charge hypotheses and isotope lookups together.
    pub max_work: usize,
}
impl Default for Deisotoper {
    fn default() -> Self {
        Self {
            tolerance: Tolerance::Ppm(10.0),
            min_charge: 1,
            max_charge: 3,
            keep_only_deisotoped: false,
            min_isotope_peaks: 3,
            max_isotope_peaks: 10,
            make_single_charged: true,
            annotate_charge: false,
            annotate_isotope_peak_count: false,
            annotate_features: false,
            use_decreasing_model: true,
            start_intensity_check: 2,
            add_up_intensity: false,
            allow_shared_isotopes: false,
            max_work: 10_000_000,
        }
    }
}
impl Deisotoper {
    /// Resolution limits from the source; native inputs must also be finite and
    /// nonnegative (the C++ predicate only compares to the upper threshold).
    pub fn is_tolerance_supported(tolerance: Tolerance) -> bool {
        match tolerance {
            Tolerance::Absolute(v) => v.is_finite() && (0.0..=0.1).contains(&v),
            Tolerance::Ppm(v) => v.is_finite() && (0.0..=100.0).contains(&v),
        }
    }
    fn validate(&self) -> Result<()> {
        if !Self::is_tolerance_supported(self.tolerance)
            || self.min_charge == 0
            || self.max_charge < self.min_charge
            || self.min_isotope_peaks < 2
            || self.max_isotope_peaks < self.min_isotope_peaks
            || self.max_work == 0
        {
            return Err(Error::InvalidValue(
                "invalid deisotoping tolerance, charge, cluster-size or work limit".into(),
            ));
        }
        Ok(())
    }

    /// Detect clusters in a sorted, centroided spectrum. Intensities must be
    /// nonnegative. The input is untouched, including when validation fails.
    ///
    /// Charge hypotheses are tried high to low; the first accepted ladder wins.
    /// By default, ladders cannot reuse peaks from other accepted ladders.
    /// allow_shared_isotopes enables the source's shared-extension behavior.
    /// An unassigned peak has charge 0, isotope count 1 and feature number -1.
    pub fn deisotope(&self, input: &MSSpectrum) -> Result<DeisotopingResult> {
        self.validate()?;
        validate_spectrum(
            input,
            self.annotate_charge,
            self.annotate_isotope_peak_count,
            self.annotate_features,
        )?;
        if input.is_empty() {
            return Ok(DeisotopingResult {
                spectrum: input.clone(),
                clusters: Vec::new(),
            });
        }
        let mut clusters = Vec::new();
        let mut features = vec![-1_i32; input.len()];
        let mut charges = vec![0_i32; input.len()];
        let mut counts = vec![1_i32; input.len()];
        let mut output = input.clone();
        let precursor_mass = if let [p] = input.precursors.as_slice() {
            let mass = p.mz * f64::from(p.charge) - PROTON_MASS_U * f64::from(p.charge);
            if !mass.is_finite() {
                return Err(Error::InvalidValue(
                    "precursor neutral mass overflows".into(),
                ));
            }
            if p.charge != 0 && mass > 0.0 {
                Some(mass)
            } else {
                None
            }
        } else {
            None
        };
        let mut work = 0_usize;
        let mut consume_work = || -> Result<()> {
            work = work
                .checked_add(1)
                .ok_or_else(|| Error::InvalidValue("deisotoping work counter overflows".into()))?;
            if work > self.max_work {
                return Err(Error::InvalidValue(
                    "deisotoping work limit exceeded".into(),
                ));
            }
            Ok(())
        };
        for current in 0..input.len() {
            if features[current] >= 0 {
                continue;
            }
            let mz = input.peaks[current].mz;
            let tolerance = match self.tolerance {
                Tolerance::Absolute(v) => v,
                Tolerance::Ppm(v) => mz * (v / 1e6),
            };
            for charge in (self.min_charge..=self.max_charge).rev() {
                consume_work()?;
                let mass = mz * f64::from(charge) - PROTON_MASS_U * f64::from(charge);
                if !mass.is_finite() || !tolerance.is_finite() {
                    return Err(Error::InvalidValue("isotope mass/window overflows".into()));
                }
                if precursor_mass.is_some_and(|p| mass > p + tolerance) {
                    continue;
                }
                let mut indices = vec![current];
                for isotope in 1..self.max_isotope_peaks {
                    consume_work()?;
                    let expected = mz + isotope as f64 * C13C12_MASSDIFF_U / f64::from(charge);
                    if !expected.is_finite() {
                        return Err(Error::InvalidValue("expected isotope m/z overflows".into()));
                    }
                    // Input order was validated once, so each lookup stays O(log n).
                    let found = nearest(&input.peaks, expected, |p| p.mz).expect("nonempty input");
                    let previous = *indices.last().expect("cluster has its seed");
                    if input.peaks[found].mz < expected - tolerance
                        || input.peaks[found].mz > expected + tolerance
                        || found <= previous
                        || (!self.allow_shared_isotopes && features[found] >= 0)
                    {
                        break;
                    }
                    if self.use_decreasing_model
                        && isotope >= self.start_intensity_check
                        && input.peaks[found].intensity > input.peaks[previous].intensity
                    {
                        break;
                    }
                    indices.push(found);
                }
                if indices.len() < self.min_isotope_peaks {
                    continue;
                }
                let feature = i32::try_from(clusters.len())
                    .map_err(|_| Error::InvalidValue("too many isotope clusters".into()))?;
                counts[current] = i32::try_from(indices.len())
                    .map_err(|_| Error::InvalidValue("too many isotopes".into()))?;
                charges[current] = i32::from(charge);
                for &index in &indices {
                    features[index] = feature;
                }
                if self.add_up_intensity {
                    output.peaks[current].intensity = checked_intensity(
                        indices
                            .iter()
                            .map(|&i| f64::from(input.peaks[i].intensity))
                            .sum(),
                    )?;
                }
                if self.make_single_charged {
                    let converted = mz * f64::from(charge) - f64::from(charge - 1) * PROTON_MASS_U;
                    if !converted.is_finite() || converted < 0.0 {
                        return Err(Error::InvalidValue("single-charge m/z is invalid".into()));
                    }
                    output.peaks[current].mz = converted;
                }
                clusters.push(IsotopeCluster {
                    charge,
                    peak_indices: indices,
                });
                break;
            }
        }
        let keep: Vec<_> = (0..input.len())
            .filter(|&i| charges[i] > 0 || (!self.keep_only_deisotoped && features[i] < 0))
            .collect();
        if self.annotate_charge {
            output.integer_data_arrays.push(DataArray {
                name: "charge".into(),
                data: charges,
                ..DataArray::default()
            });
        }
        if self.annotate_isotope_peak_count {
            output.integer_data_arrays.push(DataArray {
                name: "iso_peak_count".into(),
                data: counts,
                ..DataArray::default()
            });
        }
        if self.annotate_features {
            output.integer_data_arrays.push(DataArray {
                name: "feature_number".into(),
                data: features,
                ..DataArray::default()
            });
        }
        output.select(&keep)?;
        output.sort_by_position()?;
        Ok(DeisotopingResult {
            spectrum: output,
            clusters,
        })
    }
}
impl SpectrumFilter for Deisotoper {
    fn filter_spectrum(&self, spectrum: &mut MSSpectrum) -> Result<()> {
        *spectrum = self.deisotope(spectrum)?.spectrum;
        Ok(())
    }
}

fn validate_spectrum(
    input: &MSSpectrum,
    annotate_charge: bool,
    annotate_isotope_peak_count: bool,
    annotate_features: bool,
) -> Result<()> {
    input.validate()?;
    if input.peaks.iter().any(|p| p.mz < 0.0 || p.intensity < 0.0) {
        return Err(Error::InvalidValue(
            "deisotoping requires nonnegative m/z and intensity".into(),
        ));
    }
    if input.peaks.windows(2).any(|p| p[0].mz > p[1].mz) {
        return Err(Error::UnsortedData);
    }
    for (enabled, name) in [
        (annotate_charge, "charge"),
        (annotate_isotope_peak_count, "iso_peak_count"),
        (annotate_features, "feature_number"),
    ] {
        if enabled
            && (input.integer_data_arrays.iter().any(|a| a.name == name)
                || input.float_data_arrays.iter().any(|a| a.name == name)
                || input.string_data_arrays.iter().any(|a| a.name == name))
        {
            return Err(Error::InvalidValue(format!(
                "deisotoping annotation {name:?} already exists"
            )));
        }
    }
    Ok(())
}

/// OpenMS Poisson-averagine isotope detection with prefix KL model checks.
/// All charge hypotheses are considered; longest accepted ladder wins, with
/// higher charge breaking ties. See `docs/AVERAGINE_DEISOTOPING_SUPPORT.md`.
#[derive(Clone, Debug)]
pub struct AveragineDeisotoper {
    pub tolerance: Tolerance,
    pub min_charge: u8,
    pub max_charge: u8,
    pub keep_only_deisotoped: bool,
    pub min_isotope_peaks: usize,
    pub max_isotope_peaks: usize,
    pub make_single_charged: bool,
    pub annotate_charge: bool,
    pub annotate_isotope_peak_count: bool,
    /// Native convenience, using the simple deisotoper's feature_number array.
    pub annotate_features: bool,
    pub add_up_intensity: bool,
    /// Permit different accepted ladders to reuse an isotope extension.
    /// Seeds already assigned to a ladder are always skipped.
    pub allow_shared_isotopes: bool,
    /// Keep the strongest N peaks before clustering, after the source's fixed
    /// 0.05 intensity threshold. None disables this selection; Some(0) is invalid.
    pub top_n: Option<usize>,
    /// Bound input visits, hypotheses, isotope lookups, Poisson vector entries
    /// and every individual KL term. Sorting is performed by the standard library.
    pub max_work: usize,
}
impl Default for AveragineDeisotoper {
    fn default() -> Self {
        Self {
            tolerance: Tolerance::Ppm(10.0),
            min_charge: 1,
            max_charge: 3,
            keep_only_deisotoped: false,
            min_isotope_peaks: 2,
            max_isotope_peaks: 10,
            make_single_charged: true,
            annotate_charge: false,
            annotate_isotope_peak_count: false,
            annotate_features: false,
            add_up_intensity: false,
            allow_shared_isotopes: true,
            top_n: Some(5000),
            max_work: 10_000_000,
        }
    }
}
impl AveragineDeisotoper {
    fn validate(&self) -> Result<()> {
        Deisotoper {
            tolerance: self.tolerance,
            min_charge: self.min_charge,
            max_charge: self.max_charge,
            min_isotope_peaks: self.min_isotope_peaks,
            max_isotope_peaks: self.max_isotope_peaks,
            max_work: self.max_work,
            ..Default::default()
        }
        .validate()?;
        if self.top_n == Some(0)
            || self.max_isotope_peaks > crate::chemistry::isotopes::MAX_ISOTOPE_PEAKS
        {
            return Err(Error::InvalidValue(
                "invalid averagine top-N or isotope-vector limit".into(),
            ));
        }
        Ok(())
    }
    /// Detect isotope ladders in sorted centroid peaks. Every candidate prefix
    /// must pass the source KL threshold for its size. Members advance strictly.
    /// By default extensions can be shared, matching C++; disabling
    /// allow_shared_isotopes makes accepted ladders disjoint. Input is untouched.
    pub fn deisotope(&self, input: &MSSpectrum) -> Result<DeisotopingResult> {
        self.validate()?;
        validate_spectrum(
            input,
            self.annotate_charge,
            self.annotate_isotope_peak_count,
            self.annotate_features,
        )?;
        if input.is_empty() {
            return Ok(DeisotopingResult {
                spectrum: input.clone(),
                clusters: Vec::new(),
            });
        }
        let mut work = 0_usize;
        let mut consume = |amount: usize| -> Result<()> {
            work = work
                .checked_add(amount)
                .filter(|v| *v <= self.max_work)
                .ok_or_else(|| {
                    Error::InvalidValue("averagine deisotoping work limit exceeded".into())
                })?;
            Ok(())
        };
        consume(input.len())?;
        // Track original indices directly; no temporary user-visible annotation
        // is needed to carry them through thresholding and intensity selection.
        let mut original_indices: Vec<_> = (0..input.len())
            .filter(|&i| f64::from(input.peaks[i].intensity) >= 0.05)
            .collect();
        if let Some(n) = self.top_n {
            if original_indices.len() > n {
                original_indices.sort_by(|&a, &b| {
                    input.peaks[b]
                        .intensity
                        .total_cmp(&input.peaks[a].intensity)
                });
                original_indices.truncate(n);
                original_indices.sort_unstable();
            }
        }
        let mut selected = input.clone();
        selected.select(&original_indices)?;
        let mut output = selected.clone();
        let mut clusters = Vec::new();
        let mut features = vec![-1_i32; selected.len()];
        let mut charges = vec![0_i32; selected.len()];
        let mut counts = vec![1_i32; selected.len()];
        let precursor_mass = if let [precursor] = selected.precursors.as_slice() {
            // Preserve this operation order: the C++ precursor constraint differs
            // from the model mass expression below at floating-point boundaries.
            let charge = f64::from(precursor.charge);
            let mass = precursor.mz * charge - PROTON_MASS_U * charge;
            if !mass.is_finite() {
                return Err(Error::InvalidValue(
                    "precursor neutral mass overflows".into(),
                ));
            }
            if precursor.charge != 0 && mass > 0.0 {
                Some(mass)
            } else {
                None
            }
        } else {
            None
        };
        for current in 0..selected.len() {
            consume(1)?;
            if features[current] >= 0 {
                continue;
            }
            let mz = selected.peaks[current].mz;
            let tolerance = match self.tolerance {
                Tolerance::Absolute(v) => v,
                Tolerance::Ppm(v) => mz * (v / 1e6),
            };
            if !tolerance.is_finite() {
                return Err(Error::InvalidValue(
                    "isotope tolerance window overflows".into(),
                ));
            }
            let mut best: Option<(u8, Vec<usize>)> = None;
            for charge in self.min_charge..=self.max_charge {
                consume(1)?;
                let q = f64::from(charge);
                let mass_for_constraint = mz * q - PROTON_MASS_U * q;
                if !mass_for_constraint.is_finite() {
                    return Err(Error::InvalidValue(
                        "candidate neutral mass overflows".into(),
                    ));
                }
                if precursor_mass.is_some_and(|p| mass_for_constraint > p + tolerance) {
                    continue;
                }
                let model_mass = q * (mz - PROTON_MASS_U);
                if !model_mass.is_finite() {
                    return Err(Error::InvalidValue("averagine model mass overflows".into()));
                }
                // Nonpositive neutral masses cannot have a Poisson isotope
                // model. Keep such low-m/z noise unassigned and continue.
                if model_mass <= 0.0 {
                    continue;
                }
                let mut indices = vec![current];
                let mut distribution: Option<Vec<f64>> = None;
                let mut observed_total = f64::from(selected.peaks[current].intensity);
                let mut model_total = 0.0;
                for isotope in 1..self.max_isotope_peaks {
                    consume(1)?;
                    let expected = mz + isotope as f64 * C13C12_MASSDIFF_U / q;
                    if !expected.is_finite() {
                        return Err(Error::InvalidValue("expected isotope m/z overflows".into()));
                    }
                    let found = nearest(&selected.peaks, expected, |peak| peak.mz)
                        .expect("selected seed exists");
                    if selected.peaks[found].mz < expected - tolerance
                        || selected.peaks[found].mz > expected + tolerance
                        || found <= *indices.last().expect("candidate has its seed")
                        || (!self.allow_shared_isotopes && features[found] >= 0)
                    {
                        break;
                    }
                    // Generate once, and only after the source's first-extension
                    // precheck. Account for all entries before its allocation.
                    if distribution.is_none() {
                        consume(self.max_isotope_peaks)?;
                        let model = CoarseIsotopePatternGenerator::approximate_intensities(
                            model_mass,
                            self.max_isotope_peaks,
                        )?;
                        model_total = model[0];
                        distribution = Some(model);
                    }
                    let model = distribution.as_ref().expect("model initialized above");
                    observed_total += f64::from(selected.peaks[found].intensity);
                    model_total += model[indices.len()];
                    let mut divergence = 0.0_f32;
                    for (position, &index) in
                        indices.iter().chain(std::iter::once(&found)).enumerate()
                    {
                        consume(1)?;
                        let observed = f64::from(selected.peaks[index].intensity) / observed_total;
                        let expected_probability = model[position] / model_total;
                        // Positive observed intensity against zero model mass has
                        // infinite divergence. Never let NaN pass the comparison.
                        if expected_probability <= 0.0 || !expected_probability.is_finite() {
                            divergence = f32::INFINITY;
                            break;
                        }
                        let term = observed * (observed / expected_probability).ln();
                        divergence = (f64::from(divergence) + term) as f32;
                    }
                    const THRESHOLDS: [f32; 7] = [0.0, 0.0, 0.05, 0.1, 0.2, 0.4, 0.6];
                    if !divergence.is_finite()
                        || divergence > THRESHOLDS[(indices.len() + 1).min(6)]
                    {
                        break;
                    }
                    indices.push(found);
                }
                if indices.len() >= self.min_isotope_peaks
                    && best.as_ref().is_none_or(|(best_charge, best_indices)| {
                        indices.len() > best_indices.len()
                            || (indices.len() == best_indices.len() && charge > *best_charge)
                    })
                {
                    best = Some((charge, indices));
                }
            }
            let Some((charge, indices)) = best else {
                continue;
            };
            let feature = i32::try_from(clusters.len())
                .map_err(|_| Error::InvalidValue("too many isotope clusters".into()))?;
            counts[current] = i32::try_from(indices.len())
                .map_err(|_| Error::InvalidValue("too many isotopes".into()))?;
            charges[current] = i32::from(charge);
            for &index in &indices {
                features[index] = feature;
            }
            if self.add_up_intensity {
                output.peaks[current].intensity = checked_intensity(
                    indices
                        .iter()
                        .map(|&i| f64::from(selected.peaks[i].intensity))
                        .sum(),
                )?;
            }
            if self.make_single_charged {
                let converted = mz * f64::from(charge) - f64::from(charge - 1) * PROTON_MASS_U;
                if !converted.is_finite() || converted < 0.0 {
                    return Err(Error::InvalidValue("single-charge m/z is invalid".into()));
                }
                output.peaks[current].mz = converted;
            }
            clusters.push(IsotopeCluster {
                charge,
                peak_indices: indices.into_iter().map(|i| original_indices[i]).collect(),
            });
        }
        let keep: Vec<_> = (0..selected.len())
            .filter(|&i| charges[i] > 0 || (!self.keep_only_deisotoped && features[i] < 0))
            .collect();
        if self.annotate_charge {
            output.integer_data_arrays.push(DataArray {
                name: "charge".into(),
                data: charges,
                ..DataArray::default()
            });
        }
        if self.annotate_isotope_peak_count {
            output.integer_data_arrays.push(DataArray {
                name: "iso_peak_count".into(),
                data: counts,
                ..DataArray::default()
            });
        }
        if self.annotate_features {
            output.integer_data_arrays.push(DataArray {
                name: "feature_number".into(),
                data: features,
                ..DataArray::default()
            });
        }
        output.select(&keep)?;
        if self.make_single_charged {
            output.sort_by_position()?;
        }
        Ok(DeisotopingResult {
            spectrum: output,
            clusters,
        })
    }
}
impl SpectrumFilter for AveragineDeisotoper {
    fn filter_spectrum(&self, spectrum: &mut MSSpectrum) -> Result<()> {
        *spectrum = self.deisotope(spectrum)?.spectrum;
        Ok(())
    }
}
