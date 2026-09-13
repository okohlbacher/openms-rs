// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Spectrum alignment and sparse binned comparisons ported from OpenMS.
//! See `docs/COMPARISON_SUPPORT.md` for source provenance, precision and limits.
//!
//! The two abstract bases of the spectrum-similarity hierarchy,
//! `COMPARISON/PeakSpectrumCompareFunctor.h` and
//! `COMPARISON/BinnedSpectrumCompareFunctor.h`, are the traits
//! [`PeakSpectrumCompareFunctor`] and [`BinnedSpectrumCompareFunctor`]; the
//! three binned scorers `COMPARISON/BinnedSharedPeakCount.h`,
//! `COMPARISON/BinnedSpectralContrastAngle.h` and
//! `COMPARISON/BinnedSumAgreeingIntensities.h` are
//! [`BinnedSharedPeakCount`], [`BinnedSpectralContrastAngle`] and
//! [`BinnedSumAgreeingIntensities`]. Their support documents are
//! `docs/PEAK_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md`,
//! `docs/BINNED_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md`,
//! `docs/BINNED_SHARED_PEAK_COUNT_SUPPORT.md`,
//! `docs/BINNED_SPECTRAL_CONTRAST_ANGLE_SUPPORT.md` and
//! `docs/BINNED_SUM_AGREEING_INTENSITIES_SUPPORT.md`.

use crate::param::DefaultParamHandler;
use crate::{Error, MSSpectrum, Precursor, Result};
use std::collections::BTreeMap;

/// Matching tolerance: absolute Th or parts per million of the reference m/z.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Tolerance {
    Absolute(f64),
    Ppm(f64),
}
impl Default for Tolerance {
    fn default() -> Self {
        Self::Absolute(0.3)
    }
}
impl Tolerance {
    fn validate(self) -> Result<()> {
        let value = match self {
            Self::Absolute(v) | Self::Ppm(v) => v,
        };
        if !value.is_finite() || value < 0.0 {
            return Err(bad("tolerance must be finite and nonnegative"));
        }
        if matches!(self, Self::Ppm(_)) && !(value as f32).is_finite() {
            return Err(bad("ppm tolerance exceeds upstream f32 range"));
        }
        Ok(())
    }
    fn at(self, mz: f64) -> f64 {
        match self {
            Self::Absolute(v) => v,
            Self::Ppm(v) => v * mz * 1e-6,
        }
    }
}

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
pub(crate) fn validate_spectrum(
    spectrum: &MSSpectrum,
    nonnegative_intensities: bool,
) -> Result<()> {
    spectrum.validate()?;
    if !spectrum.is_sorted() {
        return Err(Error::UnsortedData);
    }
    if spectrum
        .peaks
        .iter()
        .any(|p| p.mz < 0.0 || (nonnegative_intensities && p.intensity < 0.0))
    {
        return Err(bad(
            "comparison requires nonnegative m/z and nonnegative scoring intensities",
        ));
    }
    Ok(())
}
pub(crate) fn finite_score(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad("comparison arithmetic overflow"))
    }
}

/// Banded dynamic programming for absolute tolerance; directed nearest matches
/// for ppm tolerance. Intensities do not affect the alignment.
#[derive(Clone, Copy, Debug)]
pub struct SpectrumAlignment {
    pub tolerance: Tolerance,
    /// Bound on DP cells plus row/column initialization, or ppm input peaks.
    pub max_cells: usize,
}
impl Default for SpectrumAlignment {
    fn default() -> Self {
        Self {
            tolerance: Tolerance::default(),
            max_cells: 5_000_000,
        }
    }
}
#[derive(Clone, Copy)]
struct Cell {
    cost: f64,
    direction: u8,
}
struct Row {
    first: usize,
    cells: Vec<Cell>,
}
impl Row {
    fn get(&self, column: usize) -> Option<Cell> {
        column
            .checked_sub(self.first)
            .and_then(|i| self.cells.get(i))
            .copied()
    }
}
/// Shared allowance for a batch of alignments. Public single-call alignment
/// retains its caller-selected max_cells; annotation batches use this ceiling.
pub(crate) struct AlignmentWork {
    remaining: usize,
}
impl Default for AlignmentWork {
    fn default() -> Self {
        Self {
            remaining: 50_000_000,
        }
    }
}
impl AlignmentWork {
    fn consume(&mut self, cells: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(cells)
            .ok_or_else(|| bad("cumulative alignment work limit exceeded"))?;
        Ok(())
    }
}
impl SpectrumAlignment {
    /// Return ordered zero-based `(reference, target)` index pairs.
    /// Absolute matching is one-to-one. Ppm matching may reuse a target peak.
    pub fn align(
        &self,
        reference: &MSSpectrum,
        target: &MSSpectrum,
    ) -> Result<Vec<(usize, usize)>> {
        self.align_with_work(
            reference,
            target,
            &mut AlignmentWork {
                remaining: self.max_cells,
            },
        )
    }

    pub(crate) fn align_with_work(
        &self,
        reference: &MSSpectrum,
        target: &MSSpectrum,
        work: &mut AlignmentWork,
    ) -> Result<Vec<(usize, usize)>> {
        self.tolerance.validate()?;
        validate_spectrum(reference, false)?;
        validate_spectrum(target, false)?;
        if self.max_cells == 0 {
            return Err(bad("alignment cell limit must be positive"));
        }
        let mut used = reference
            .len()
            .checked_add(target.len())
            .and_then(|n| n.checked_add(1))
            .filter(|&n| n <= self.max_cells)
            .ok_or_else(|| bad("alignment exceeds configured cell limit"))?;
        work.consume(used)?;
        if reference.is_empty() || target.is_empty() {
            return Ok(Vec::new());
        }
        match self.tolerance {
            Tolerance::Ppm(_) => matched_alignment(reference, target, self.tolerance),
            Tolerance::Absolute(tolerance) => {
                finite_score(used as f64 * tolerance)?;
                let mut rows = Vec::<Row>::with_capacity(reference.len());
                let mut left = 1;
                let mut last = (0, 0);
                for i in 1..=reference.len() {
                    let pos1 = reference.peaks[i - 1].mz;
                    let mut row = Row {
                        first: left,
                        cells: Vec::new(),
                    };
                    for j in row.first..=target.len() {
                        if used == self.max_cells {
                            return Err(bad("alignment exceeds configured cell limit"));
                        }
                        work.consume(1)?;
                        used += 1;
                        let pos2 = target.peaks[j - 1].mz;
                        let distance = (pos1 - pos2).abs();
                        let off_band = pos2 > pos1
                            && distance > tolerance
                            && i < reference.len()
                            && j < target.len()
                            && reference.peaks[i].mz < pos2;
                        if pos1 > pos2 && distance > tolerance && j > left + 1 {
                            left += 1;
                        }
                        // Missing band cells have the same all-gap fallback as upstream maps.
                        let previous = |column: usize| -> Option<f64> {
                            if i == 1 {
                                Some(column as f64 * tolerance)
                            } else if column == 0 {
                                Some((i - 1) as f64 * tolerance)
                            } else {
                                rows[i - 2].get(column).map(|cell| cell.cost)
                            }
                        };
                        let diagonal =
                            distance + previous(j - 1).unwrap_or((i + j - 2) as f64 * tolerance);
                        let up = tolerance
                            + if j == 1 {
                                i as f64 * tolerance
                            } else {
                                row.get(j - 1)
                                    .map_or((i + j - 1) as f64 * tolerance, |cell| cell.cost)
                            };
                        let across =
                            tolerance + previous(j).unwrap_or((i + j - 1) as f64 * tolerance);
                        let cell = if diagonal <= up && diagonal <= across && distance <= tolerance
                        {
                            last = (i, j);
                            Cell {
                                cost: diagonal,
                                direction: 0,
                            }
                        } else if up <= across {
                            Cell {
                                cost: up,
                                direction: 1,
                            }
                        } else {
                            Cell {
                                cost: across,
                                direction: 2,
                            }
                        };
                        finite_score(cell.cost)?;
                        row.cells.push(cell);
                        if off_band {
                            break;
                        }
                    }
                    rows.push(row);
                }
                // Trace from the last chosen diagonal, not necessarily the bottom-right cell.
                let (mut i, mut j) = last;
                let mut pairs = Vec::new();
                while i > 0 && j > 0 {
                    let Some(cell) = rows[i - 1].get(j) else {
                        break;
                    };
                    match cell.direction {
                        0 => {
                            pairs.push((i - 1, j - 1));
                            i -= 1;
                            j -= 1;
                        }
                        1 => j -= 1,
                        _ => i -= 1,
                    }
                }
                pairs.reverse();
                Ok(pairs)
            }
        }
    }
}

pub(crate) fn matched_alignment(
    reference: &MSSpectrum,
    target: &MSSpectrum,
    tolerance: Tolerance,
) -> Result<Vec<(usize, usize)>> {
    // This forward traversal deliberately retains MatchedIterator's f32 distance
    // and tolerance arithmetic, including stopping at equal-distance duplicates.
    let (Tolerance::Absolute(value) | Tolerance::Ppm(value)) = tolerance;
    if !value.is_finite() || value < 0.0 || !(value as f32).is_finite() {
        return Err(bad(
            "matched-iterator tolerance must fit finite nonnegative f32",
        ));
    }
    if reference.is_empty() || target.is_empty() {
        return Ok(Vec::new());
    }
    if reference
        .peaks
        .iter()
        .chain(&target.peaks)
        .any(|p| !(p.mz as f32).is_finite())
    {
        return Err(bad("ppm alignment coordinates exceed upstream f32 range"));
    }
    let mut pairs = Vec::new();
    let mut j = 0;
    for (i, peak) in reference.peaks.iter().enumerate() {
        let allowed = match tolerance {
            Tolerance::Absolute(value) => value as f32,
            Tolerance::Ppm(value) => (value as f32 / 1e6_f32) * peak.mz as f32,
        };
        if !allowed.is_finite() {
            return Err(bad("ppm tolerance window overflow"));
        }
        let mut distance = (peak.mz - target.peaks[j].mz).abs() as f32;
        while j + 1 < target.len() {
            let next = (peak.mz - target.peaks[j + 1].mz).abs() as f32;
            if next >= distance {
                break;
            }
            j += 1;
            distance = next;
        }
        if distance <= allowed {
            pairs.push((i, j));
        }
    }
    Ok(pairs)
}

/// Distance weighting used by aligned and Zhang scores.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum DistanceWeighting {
    #[default]
    None,
    Linear,
    Gaussian,
}
impl DistanceWeighting {
    fn factor(self, distance: f64, tolerance: f64) -> Result<f64> {
        if self == Self::None {
            return Ok(1.0);
        }
        if tolerance == 0.0 && distance == 0.0 {
            return Ok(1.0);
        }
        if !tolerance.is_finite() || tolerance <= 0.0 || distance > tolerance {
            return Err(bad(
                "distance weighting requires a finite in-window tolerance",
            ));
        }
        Ok(match self {
            Self::None => 1.0,
            Self::Linear => (tolerance - distance) / tolerance,
            Self::Gaussian => erfc_small((distance / tolerance) / (3.0 * std::f64::consts::SQRT_2)),
        })
    }
}
// Convergent erf power series: inputs are [0, 1/(3 sqrt(2))], not a general erfc.
fn erfc_small(x: f64) -> f64 {
    let mut term = x;
    let mut sum = x;
    for n in 1..24 {
        term *= -x * x / f64::from(n);
        let add = term / f64::from(2 * n + 1);
        sum += add;
        if add.abs() < 1e-18 {
            break;
        }
    }
    1.0 - 2.0 / std::f64::consts::PI.sqrt() * sum
}

/// OpenMS alignment score: sum sqrt(I1 I2 factor) / sqrt(sum I1² sum I2²).
/// This is not a cosine: a self-score need not equal one and scores can exceed one.
#[derive(Clone, Copy, Debug, Default)]
pub struct SpectrumAlignmentScore {
    pub alignment: SpectrumAlignment,
    pub weighting: DistanceWeighting,
}
impl SpectrumAlignmentScore {
    /// Similarity of `reference` and `target`, as `SpectrumAlignmentScore::operator()`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] when either spectrum is not sorted by m/z,
    /// and [`Error::InvalidValue`] for a negative or non-finite intensity, an
    /// invalid tolerance, an alignment that exceeds
    /// [`SpectrumAlignment::max_cells`], or a non-finite score. A zero intensity
    /// norm on either side yields `Ok(0.0)` rather than a division by zero.
    pub fn score(&self, reference: &MSSpectrum, target: &MSSpectrum) -> Result<f64> {
        validate_spectrum(reference, true)?;
        validate_spectrum(target, true)?;
        let pairs = self.alignment.align(reference, target)?;
        let denominator = norm(reference) * norm(target);
        if denominator == 0.0 {
            return Ok(0.0);
        }
        let mut sum = 0.0;
        for (i, j) in pairs {
            let p = reference.peaks[i];
            let q = target.peaks[j];
            let factor = self
                .weighting
                .factor((p.mz - q.mz).abs(), self.alignment.tolerance.at(p.mz))?;
            sum += (f64::from(p.intensity) * f64::from(q.intensity) * factor).sqrt();
        }
        finite_score(sum / denominator)
    }
}
fn norm(spectrum: &MSSpectrum) -> f64 {
    spectrum
        .peaks
        .iter()
        .map(|p| f64::from(p.intensity).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Sparse bin-coordinate convention.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum BinUnit {
    #[default]
    Absolute,
    Ppm,
}
/// Binning layout and explicit resource limits. Bin arithmetic retains f32 source precision.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinConfig {
    pub size: f32,
    pub unit: BinUnit,
    pub spread: u32,
    pub offset: f32,
    pub max_bins: usize,
    pub max_updates: usize,
}
impl Default for BinConfig {
    fn default() -> Self {
        Self {
            size: 1.0005,
            unit: BinUnit::Absolute,
            spread: 0,
            offset: 0.4,
            max_bins: 1_000_000,
            max_updates: 10_000_000,
        }
    }
}
impl BinConfig {
    fn validate(self) -> Result<()> {
        if !self.size.is_finite()
            || self.size <= 0.0
            || !self.offset.is_finite()
            || self.max_bins == 0
            || self.max_updates == 0
        {
            return Err(bad("invalid bin size, offset, or resource limit"));
        }
        Ok(())
    }
    /// Source-compatible float-coordinate bin index; ppm bins require m/z >= 1.
    pub fn bin_index(self, mz: f64) -> Result<usize> {
        self.validate()?;
        // Validate the caller's coordinate before source-compatible narrowing:
        // tiny negatives and values just below one can round onto the boundary.
        if !mz.is_finite() || mz < 0.0 || (self.unit == BinUnit::Ppm && mz < 1.0) {
            return Err(bad("invalid m/z for binning"));
        }
        let mz = mz as f32;
        if !mz.is_finite() || mz < 0.0 || (self.unit == BinUnit::Ppm && mz < 1.0) {
            return Err(bad("invalid m/z for binning"));
        }
        let index = match self.unit {
            BinUnit::Absolute => f64::from((mz / self.size + self.offset).floor()),
            BinUnit::Ppm => (f64::from(mz.ln()) / (f64::from(self.size) * 1e-6).ln_1p()).floor(),
        };
        if !index.is_finite()
            || !(0.0..=9_007_199_254_740_990.0_f64.min((usize::MAX - 1) as f64)).contains(&index)
        {
            return Err(bad("bin index exceeds supported integer range"));
        }
        Ok(index as usize)
    }
    /// Lower m/z bound of bin `index`, as `BinnedSpectrum::getBinLowerMZ`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an invalid layout and when the bound
    /// overflows `f32`; the source returns the overflowed value.
    pub fn bin_lower_mz(self, index: usize) -> Result<f32> {
        self.validate()?;
        let result = match self.unit {
            BinUnit::Absolute => (index as f32 - self.offset) * self.size,
            BinUnit::Ppm => (1.0 + f64::from(self.size) * 1e-6).powf(index as f64) as f32,
        };
        if result.is_finite() {
            Ok(result)
        } else {
            Err(bad("bin lower m/z overflow"))
        }
    }
}

/// BinnedSpectrum from upstream KERNEL, located beside comparisons in Rust.
/// Private bins protect layout/value invariants; stored zero entries are retained.
#[derive(Clone, Debug)]
pub struct BinnedSpectrum {
    config: BinConfig,
    bins: BTreeMap<usize, f32>,
    precursors: Vec<Precursor>,
}
impl PartialEq for BinnedSpectrum {
    fn eq(&self, other: &Self) -> bool {
        self.is_compatible(other)
            && self.config.spread == other.config.spread
            && self.bins == other.bins
            && self.precursors == other.precursors
    }
}
impl BinnedSpectrum {
    /// Bin `spectrum` onto `config`, as the source's five-argument constructor.
    ///
    /// Each peak adds its intensity to its own bin and, when
    /// [`BinConfig::spread`] is nonzero, to that many neighbours on each side;
    /// bin 0 stops the downward spread. The precursor list is copied.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] for unsorted peaks, which the source only
    /// asserts in debug builds, and [`Error::InvalidValue`] for an invalid
    /// layout, an m/z below 1 under ppm binning, `f32` intensity overflow, or
    /// work beyond [`BinConfig::max_bins`] / [`BinConfig::max_updates`].
    pub fn new(spectrum: &MSSpectrum, config: BinConfig) -> Result<Self> {
        config.validate()?;
        validate_spectrum(spectrum, false)?;
        let updates = (config.spread as usize)
            .checked_mul(2)
            .and_then(|n| n.checked_add(1))
            .and_then(|n| n.checked_mul(spectrum.len()))
            .filter(|&n| n <= config.max_updates)
            .ok_or_else(|| bad("binning exceeds configured update limit"))?;
        let _ = updates;
        let mut bins = BTreeMap::new();
        for peak in &spectrum.peaks {
            if config.unit == BinUnit::Ppm && peak.mz < 1.0 {
                return Err(bad("ppm binning requires m/z >= 1"));
            }
            let index = config.bin_index(peak.mz)?;
            let last = index
                .checked_add(config.spread as usize)
                .ok_or_else(|| bad("spread bin index overflow"))?;
            for i in index.saturating_sub(config.spread as usize)..=last {
                if !bins.contains_key(&i) && bins.len() == config.max_bins {
                    return Err(bad("binning exceeds configured bin limit"));
                }
                let value = bins.entry(i).or_insert(0.0_f32);
                *value += peak.intensity;
                if !value.is_finite() {
                    return Err(bad("binned intensity overflow"));
                }
            }
        }
        Ok(Self {
            config,
            bins,
            precursors: spectrum.precursors.clone(),
        })
    }
    /// Binning layout, covering the source's `getBinSize`, `getBinSpread`,
    /// `getOffset` and the unit flag in one value.
    pub fn config(&self) -> BinConfig {
        self.config
    }
    /// Stored bins by index, as the source's const `getBins`. Explicit zeros are
    /// retained, so this is `nonZeros()` in the Eigen sense, not "nonzero values".
    pub fn bins(&self) -> &BTreeMap<usize, f32> {
        &self.bins
    }
    /// Precursors copied from the binned spectrum, as the source's const
    /// `getPrecursors`.
    pub fn precursors(&self) -> &[Precursor] {
        &self.precursors
    }
    /// Layout compatibility ignores spread, matching upstream isCompatible.
    pub fn is_compatible(&self, other: &Self) -> bool {
        self.config.unit == other.config.unit
            && self.config.size == other.config.size
            && self.config.offset == other.config.offset
    }
    /// Lookup does not insert a zero entry (unlike upstream coeffRef access).
    pub fn bin_intensity(&self, mz: f64) -> Result<f32> {
        Ok(*self.bins.get(&self.config.bin_index(mz)?).unwrap_or(&0.0))
    }
}
fn compatible(a: &BinnedSpectrum, b: &BinnedSpectrum) -> Result<()> {
    if a.is_compatible(b) {
        Ok(())
    } else {
        Err(bad("incompatible binned spectrum layouts"))
    }
}
/// OpenMS BinnedSpectralContrastAngle is cosine similarity, not acos(cosine).
/// Signed bin values are supported; a zero norm yields zero.
pub fn binned_cosine(a: &BinnedSpectrum, b: &BinnedSpectrum) -> Result<f64> {
    compatible(a, b)?;
    let norm_a = a
        .bins
        .values()
        .map(|&v| f64::from(v).powi(2))
        .sum::<f64>()
        .sqrt();
    let norm_b = b
        .bins
        .values()
        .map(|&v| f64::from(v).powi(2))
        .sum::<f64>()
        .sqrt();
    if norm_a == 0.0 || norm_b == 0.0 {
        return Ok(0.0);
    }
    let dot = a
        .bins
        .iter()
        .map(|(i, &v)| f64::from(v) * f64::from(*b.bins.get(i).unwrap_or(&0.0)))
        .sum::<f64>();
    finite_score((dot / (norm_a * norm_b)).clamp(-1.0, 1.0))
}
/// Stored-bin intersection divided by the larger stored-bin count, as in Eigen.
/// Explicit zeros still count as stored bins. Two empty spectra yield zero.
pub fn binned_shared_peak_count(a: &BinnedSpectrum, b: &BinnedSpectrum) -> Result<f64> {
    compatible(a, b)?;
    let denominator = a.bins.len().max(b.bins.len());
    if denominator == 0 {
        return Ok(0.0);
    }
    Ok(a.bins.keys().filter(|i| b.bins.contains_key(i)).count() as f64 / denominator as f64)
}
/// Sum max(0, (a+b)/2 - |a-b|), divided by the mean total intensity.
/// Negative bins are rejected; a zero denominator yields zero.
pub fn binned_sum_agreeing_intensities(a: &BinnedSpectrum, b: &BinnedSpectrum) -> Result<f64> {
    compatible(a, b)?;
    if a.bins.values().chain(b.bins.values()).any(|&v| v < 0.0) {
        return Err(bad("agreeing intensities requires nonnegative bins"));
    }
    let total = a
        .bins
        .values()
        .chain(b.bins.values())
        .map(|&v| f64::from(v))
        .sum::<f64>()
        * 0.5;
    if total == 0.0 {
        return Ok(0.0);
    }
    let sum = a
        .bins
        .iter()
        .filter_map(|(i, &v)| {
            b.bins.get(i).map(|&w| {
                let v = f64::from(v);
                let w = f64::from(w);
                ((v + w) * 0.5 - (v - w).abs()).max(0.0)
            })
        })
        .sum::<f64>();
    finite_score((sum / total).min(1.0))
}

/// Similarity in Th: max(0, window - |first precursor m/z difference|).
/// Missing precursors are treated as m/z zero, preserving the source convention.
#[derive(Clone, Copy, Debug)]
pub struct SpectrumPrecursorComparator {
    pub window: f64,
}
impl Default for SpectrumPrecursorComparator {
    fn default() -> Self {
        Self { window: 2.0 }
    }
}
impl SpectrumPrecursorComparator {
    /// Precursor similarity of `a` and `b`, as `SpectrumPrecursorComparator::operator()`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a non-finite or negative window and
    /// for an invalid spectrum.
    pub fn score(&self, a: &MSSpectrum, b: &MSSpectrum) -> Result<f64> {
        a.validate()?;
        b.validate()?;
        if !self.window.is_finite() || self.window < 0.0 {
            return Err(bad("precursor window must be finite and nonnegative"));
        }
        let a = a.precursors.first().map_or(0.0, |p| p.mz);
        let b = b.precursors.first().map_or(0.0, |p| p.mz);
        Ok((self.window - (a - b).abs()).max(0.0))
    }
}

/// Zhang many-to-many score. Uses a strict absolute tolerance boundary.
/// The per-instance Gaussian scale corrects upstream's static-first-call cache.
#[derive(Clone, Copy, Debug)]
pub struct ZhangSimilarityScore {
    pub tolerance: f64,
    pub weighting: DistanceWeighting,
    pub max_pairs: usize,
}
impl Default for ZhangSimilarityScore {
    fn default() -> Self {
        Self {
            tolerance: 0.2,
            weighting: DistanceWeighting::None,
            max_pairs: 5_000_000,
        }
    }
}
impl ZhangSimilarityScore {
    /// Similarity of `a` and `b`, as `ZhangSimilarityScore::operator()`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] for unsorted peaks and
    /// [`Error::InvalidValue`] for a negative or non-finite intensity, an
    /// invalid tolerance, more than [`ZhangSimilarityScore::max_pairs`]
    /// candidate pairs, or a non-finite score. A zero total intensity on either
    /// side yields `Ok(0.0)`.
    pub fn score(&self, a: &MSSpectrum, b: &MSSpectrum) -> Result<f64> {
        let mut sum = 0.0;
        for_each_pair(a, b, self.tolerance, false, self.max_pairs, |i, j| {
            let p = a.peaks[i];
            let q = b.peaks[j];
            let weight = self.weighting.factor((p.mz - q.mz).abs(), self.tolerance)?;
            sum += (f64::from(p.intensity) * f64::from(q.intensity) * weight).sqrt();
            Ok(())
        })?;
        let denominator = (total(a) * total(b)).sqrt();
        if denominator == 0.0 {
            Ok(0.0)
        } else {
            finite_score(sum / denominator)
        }
    }
}
/// Stein/Scott improved score, including all pairs within twice the tolerance.
#[derive(Clone, Copy, Debug)]
pub struct SteinScottImproveScore {
    pub tolerance: f64,
    pub threshold: f32,
    pub max_pairs: usize,
}
impl Default for SteinScottImproveScore {
    fn default() -> Self {
        Self {
            tolerance: 0.2,
            threshold: 0.2,
            max_pairs: 5_000_000,
        }
    }
}
impl SteinScottImproveScore {
    /// Similarity of `a` and `b`, as `SteinScottImproveScore::operator()`.
    ///
    /// Scores below [`SteinScottImproveScore::threshold`] are reported as zero.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] for unsorted peaks and
    /// [`Error::InvalidValue`] for a negative or non-finite intensity, a
    /// non-finite threshold, an invalid tolerance, more than
    /// [`SteinScottImproveScore::max_pairs`] candidate pairs, or a non-finite
    /// score. A zero intensity norm on either side yields `Ok(0.0)`.
    pub fn score(&self, a: &MSSpectrum, b: &MSSpectrum) -> Result<f64> {
        if !self.threshold.is_finite() {
            return Err(bad("score threshold must be finite"));
        }
        Tolerance::Absolute(self.tolerance).validate()?;
        let mut sum = 0.0;
        for_each_pair(a, b, 2.0 * self.tolerance, true, self.max_pairs, |i, j| {
            sum += f64::from(a.peaks[i].intensity) * f64::from(b.peaks[j].intensity);
            Ok(())
        })?;
        let denominator = norm(a) * norm(b);
        if denominator == 0.0 {
            return Ok(0.0);
        }
        let score =
            finite_score((sum - self.tolerance / 10000.0 * total(a) * total(b)) / denominator)?;
        Ok(if score < f64::from(self.threshold) {
            0.0
        } else {
            score
        })
    }
}
fn total(spectrum: &MSSpectrum) -> f64 {
    spectrum.peaks.iter().map(|p| f64::from(p.intensity)).sum()
}
fn for_each_pair(
    a: &MSSpectrum,
    b: &MSSpectrum,
    tolerance: f64,
    inclusive: bool,
    max_pairs: usize,
    mut visit: impl FnMut(usize, usize) -> Result<()>,
) -> Result<()> {
    Tolerance::Absolute(tolerance).validate()?;
    validate_spectrum(a, true)?;
    validate_spectrum(b, true)?;
    if max_pairs == 0 {
        return Err(bad("pair limit must be positive"));
    }
    let mut first = 0;
    let mut used = 0;
    for (i, p) in a.peaks.iter().enumerate() {
        while first < b.len() && b.peaks[first].mz < p.mz && p.mz - b.peaks[first].mz > tolerance {
            first += 1;
        }
        for j in first..b.len() {
            let distance = (p.mz - b.peaks[j].mz).abs();
            if b.peaks[j].mz > p.mz && distance > tolerance {
                break;
            }
            if used == max_pairs {
                return Err(bad("comparison exceeds configured pair limit"));
            }
            used += 1;
            if distance < tolerance || (inclusive && distance == tolerance) {
                visit(i, j)?;
            }
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// The spectrum-similarity functor hierarchy:
//   COMPARISON/PeakSpectrumCompareFunctor.h
//   COMPARISON/BinnedSpectrumCompareFunctor.h
//   COMPARISON/BinnedSharedPeakCount.h
//   COMPARISON/BinnedSpectralContrastAngle.h
//   COMPARISON/BinnedSumAgreeingIntensities.h
// ---------------------------------------------------------------------------

/// Largest number of stored bins a single binned comparison may examine, summed
/// over both spectra.
///
/// [`BinnedSpectrum`] already caps its own bin count at
/// [`BinConfig::max_bins`], so this is the second, comparison-side ceiling: it
/// is checked before any traversal begins, and no scorer allocates, so a
/// rejected comparison leaves both inputs untouched. The source has no such
/// ceiling; its cost is `O(bins)` per call regardless.
pub const MAX_COMPARED_BINS: usize = 4_000_000;

/// Reproduce the source's two-step construction: the abstract base names the
/// handler after itself, then the concrete functor renames it and copies its
/// (empty) defaults into the current parameters, as `defaultsToParam_` does.
fn functor_handler(base: &str, name: &str) -> Result<DefaultParamHandler> {
    let mut handler = DefaultParamHandler::new(base)?;
    handler.set_name(name)?;
    handler.defaults_to_parameters()?;
    Ok(handler)
}

/// Base for compare functors of spectra, that return a similarity value for two
/// spectra.
///
/// Implementors return a similarity value for a pair of spectra. The value
/// should be greater than or equal to zero.
///
/// The source class is an abstract `DefaultParamHandler` subclass with two pure
/// virtual `operator()` overloads. Here the two-spectrum overload is
/// [`score`](Self::score), the one-spectrum overload is
/// [`self_score`](Self::self_score), and the inherited parameter surface is
/// reached through [`handler`](Self::handler) and
/// [`handler_mut`](Self::handler_mut). The trait is object safe, so
/// `&dyn PeakSpectrumCompareFunctor` replaces the base-class pointer the source
/// hierarchy exists for.
///
/// Every derived functor in `COMPARISON/` overrides the one-spectrum overload as
/// `operator()(spec, spec)`, so the default [`self_score`](Self::self_score) is
/// that delegation and an implementor only overrides it to record a cheaper
/// closed form.
///
/// The earlier-ported peak comparators in this module -
/// [`SpectrumAlignmentScore`], [`ZhangSimilarityScore`],
/// [`SteinScottImproveScore`] and [`SpectrumPrecursorComparator`] - are source
/// descendants of this base but keep their configuration in typed `Copy` fields
/// instead of a [`DefaultParamHandler`], so they do not implement this trait
/// yet; see `docs/PEAK_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md`.
pub trait PeakSpectrumCompareFunctor {
    /// The parameter surface the source inherits from `DefaultParamHandler`,
    /// carrying the functor's registered name and its current parameters.
    fn handler(&self) -> &DefaultParamHandler;

    /// Mutable parameter surface, for the source's public `setParameters` and
    /// `setName`. The source's `updateMembers_` hook has no equivalent: a
    /// functor that derives typed state from parameters recomputes it here.
    fn handler_mut(&mut self) -> &mut DefaultParamHandler;

    /// Registered functor name, as `DefaultParamHandler::getName`.
    fn name(&self) -> &str {
        self.handler().name()
    }

    /// Similarity of `a` and `b`.
    ///
    /// # Errors
    ///
    /// Implementation defined; the source signature cannot fail but every
    /// concrete functor here reports invalid input and unbounded work instead of
    /// returning a wrong or non-finite score.
    fn score(&self, a: &MSSpectrum, b: &MSSpectrum) -> Result<f64>;

    /// Self similarity, `score(a, a)`.
    ///
    /// # Errors
    ///
    /// As [`score`](Self::score).
    fn self_score(&self, a: &MSSpectrum) -> Result<f64> {
        self.score(a, a)
    }
}

/// Base for compare functors of binned spectra.
///
/// Implementors return a value for a pair of [`BinnedSpectrum`] objects, or a
/// single one with itself. Ideally the value reflects the similarity of the
/// pair; for how each similarity is computed see the concrete functors.
///
/// The source class comment states that functors normalised to `[0, 1]` are
/// identifiable by a set `normalized` parameter. No binned functor in the
/// pinned revision registers such a parameter - only `PeakAlignment`, which
/// derives from [`PeakSpectrumCompareFunctor`], does - so that sentence does not
/// describe this hierarchy. All three binned functors here are normalised to
/// `[0, 1]` for nonnegative bins, and each says so in its own documentation.
///
/// As with [`PeakSpectrumCompareFunctor`], the two pure virtual `operator()`
/// overloads are [`score`](Self::score) and [`self_score`](Self::self_score),
/// the trait is object safe, and every source derivative implements the
/// one-spectrum overload as `operator()(spec, spec)`.
pub trait BinnedSpectrumCompareFunctor {
    /// The parameter surface the source inherits from `DefaultParamHandler`.
    fn handler(&self) -> &DefaultParamHandler;

    /// Mutable parameter surface, for the source's public `setParameters` and
    /// `setName`.
    fn handler_mut(&mut self) -> &mut DefaultParamHandler;

    /// Registered functor name, as `DefaultParamHandler::getName`.
    fn name(&self) -> &str {
        self.handler().name()
    }

    /// Similarity of `spec1` and `spec2`.
    ///
    /// # Errors
    ///
    /// Implementation defined; every functor here reports incompatible binning,
    /// work beyond [`MAX_COMPARED_BINS`] and non-finite arithmetic as
    /// [`Error::InvalidValue`].
    fn score(&self, spec1: &BinnedSpectrum, spec2: &BinnedSpectrum) -> Result<f64>;

    /// Self similarity, `score(spec, spec)`.
    ///
    /// # Errors
    ///
    /// As [`score`](Self::score).
    fn self_score(&self, spec: &BinnedSpectrum) -> Result<f64> {
        self.score(spec, spec)
    }
}

/// Refuse a comparison whose two spectra together store more than
/// [`MAX_COMPARED_BINS`] bins, before anything is traversed.
fn preflight_bins(a: &BinnedSpectrum, b: &BinnedSpectrum) -> Result<()> {
    let total = a
        .bins
        .len()
        .checked_add(b.bins.len())
        .ok_or_else(|| bad("binned comparison size overflow"))?;
    if total > MAX_COMPARED_BINS {
        return Err(bad("binned comparison exceeds the stored-bin ceiling"));
    }
    Ok(())
}

/// `Eigen::SparseVector<float>::dot`: the products of the coefficients stored at
/// indices present in both vectors, accumulated in `f32` in ascending index
/// order and widened only on return. The source assigns that `float` result to a
/// `double`, so the reduction's precision, not `double`'s, is what the score
/// inherits.
fn sparse_dot(a: &BinnedSpectrum, b: &BinnedSpectrum) -> Result<f64> {
    let mut total = 0.0_f32;
    for (index, &value) in &a.bins {
        if let Some(&other) = b.bins.get(index) {
            total += value * other;
        }
    }
    if total.is_finite() {
        Ok(f64::from(total))
    } else {
        Err(bad("binned dot product overflows f32"))
    }
}

/// `Eigen::SparseVector<float>::sum`: the stored coefficients accumulated in
/// `f32`, here in ascending index order.
fn sparse_sum(spectrum: &BinnedSpectrum) -> Result<f64> {
    let mut total = 0.0_f32;
    for &value in spectrum.bins.values() {
        total += value;
    }
    if total.is_finite() {
        Ok(f64::from(total))
    } else {
        Err(bad("binned intensity sum overflows f32"))
    }
}

/// Square root that refuses a negative or non-finite radicand instead of
/// yielding NaN. The source calls `sqrt` unguarded.
fn checked_sqrt(value: f64) -> Result<f64> {
    if value.is_finite() && value >= 0.0 {
        Ok(value.sqrt())
    } else {
        Err(bad("comparison norm is not a nonnegative finite number"))
    }
}

/// Compare functor scoring the shared peaks for similarity measurement.
///
/// The score is the number of bins occupied in both spectra divided by the
/// larger of the two occupied-bin counts, which normalises it to `[0, 1]`.
/// "Occupied" means *stored*, as `Eigen::SparseVector::nonZeros` counts stored
/// coefficients: a bin that a peak of intensity zero created, or that a spread
/// wrote a zero into, counts as occupied even though its value is zero.
///
/// The details of the score can be found in: K. Wan, I. Vidavsky, and M. Gross.
/// Comparing similar spectra: from similarity index to spectral contrast angle.
/// Journal of the American Society for Mass Spectrometry, 13(1):85-88, January
/// 2002.
///
/// The source registers no parameters, so `getParameters()` is empty and the
/// `@htmlinclude OpenMS_BinnedSharedPeakCount.parameters` block it documents is
/// empty too. Its protected `precursor_mass_tolerance_` member is never
/// initialised, never written by `updateMembers_`, never copied by the copy
/// constructor and never read; it has no counterpart here.
#[derive(Clone, Debug, PartialEq)]
pub struct BinnedSharedPeakCount {
    handler: DefaultParamHandler,
}

impl BinnedSharedPeakCount {
    /// Construct the functor with the source's name and its empty parameters.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] only if the parameter handler rejects the
    /// fixed name, which cannot happen for this literal.
    pub fn new() -> Result<Self> {
        Ok(Self {
            handler: functor_handler("BinnedSpectrumCompareFunctor", "BinnedSharedPeakCount")?,
        })
    }
}

impl BinnedSpectrumCompareFunctor for BinnedSharedPeakCount {
    fn handler(&self) -> &DefaultParamHandler {
        &self.handler
    }

    fn handler_mut(&mut self) -> &mut DefaultParamHandler {
        &mut self.handler
    }

    /// Shared occupied bins over the larger occupied-bin count.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the two spectra do not share a
    /// binning, which is the source's `Exception::IllegalArgument`, and when the
    /// two spectra together store more than [`MAX_COMPARED_BINS`] bins.
    ///
    /// Two spectra that both store no bins yield `Ok(0.0)`. The source divides
    /// by that zero denominator and returns NaN.
    fn score(&self, spec1: &BinnedSpectrum, spec2: &BinnedSpectrum) -> Result<f64> {
        compatible(spec1, spec2)?;
        preflight_bins(spec1, spec2)?;
        let denominator = spec1.bins.len().max(spec2.bins.len());
        if denominator == 0 {
            return Ok(0.0);
        }
        let shared = spec1
            .bins
            .keys()
            .filter(|index| spec2.bins.contains_key(index))
            .count();
        Ok(shared as f64 / denominator as f64)
    }
}

/// Compare functor scoring the spectral contrast angle for similarity
/// measurement.
///
/// Despite the name the source returns the *cosine* of the spectral contrast
/// angle, not the angle: `dot(a, b) / sqrt(dot(a, a) * dot(b, b))`. For
/// nonnegative bins that is in `[0, 1]`; bins may be negative in principle, in
/// which case the score is in `[-1, 1]` and is not clamped, exactly as in the
/// source.
///
/// The details of the score can be found in: K. Wan, I. Vidavsky, and M. Gross.
/// Comparing similar spectra: from similarity index to spectral contrast angle.
/// Journal of the American Society for Mass Spectrometry, 13(1):85-88, January
/// 2002.
///
/// The denominator is `sqrt(sum1 * sum2)` and not `sqrt(sum1) * sqrt(sum2)`;
/// the two differ in the last bits and the source's grouping is kept.
///
/// The source registers no parameters, and its protected
/// `precursor_mass_tolerance_` member is unused; see [`BinnedSharedPeakCount`].
///
/// ```
/// use openms::comparison::{
///     BinConfig, BinnedSpectralContrastAngle, BinnedSpectrum, BinnedSpectrumCompareFunctor,
/// };
/// use openms::{MSSpectrum, Peak1D};
///
/// let spectrum = MSSpectrum::from_peaks(vec![
///     Peak1D::new(100.0, 2.0),
///     Peak1D::new(200.0, 3.0),
/// ]);
/// let binned = BinnedSpectrum::new(&spectrum, BinConfig::default())?;
/// let functor = BinnedSpectralContrastAngle::new()?;
/// assert_eq!(functor.name(), "BinnedSpectralContrastAngle");
/// assert_eq!(functor.self_score(&binned)?, 1.0);
/// # Ok::<(), openms::Error>(())
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct BinnedSpectralContrastAngle {
    handler: DefaultParamHandler,
}

impl BinnedSpectralContrastAngle {
    /// Construct the functor with the source's name and its empty parameters.
    ///
    /// # Errors
    ///
    /// As [`BinnedSharedPeakCount::new`].
    pub fn new() -> Result<Self> {
        Ok(Self {
            handler: functor_handler(
                "BinnedSpectrumCompareFunctor",
                "BinnedSpectralContrastAngle",
            )?,
        })
    }
}

impl BinnedSpectrumCompareFunctor for BinnedSpectralContrastAngle {
    fn handler(&self) -> &DefaultParamHandler {
        &self.handler
    }

    fn handler_mut(&mut self) -> &mut DefaultParamHandler {
        &mut self.handler
    }

    /// Cosine of the spectral contrast angle between the two bin vectors.
    ///
    /// An empty or all-zero spectrum makes `sum1 * sum2` zero; the source
    /// returns a defined score of `0` there rather than `0.0 / 0.0`, and so does
    /// this.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the two spectra do not share a
    /// binning, when they together store more than [`MAX_COMPARED_BINS`] bins,
    /// and when an `f32` dot product overflows. The source only asserts
    /// compatible binning through `OPENMS_PRECONDITION`, which is a no-op in
    /// release builds, so incompatible input silently scores two different
    /// binnings against each other there; its two sibling functors throw.
    fn score(&self, spec1: &BinnedSpectrum, spec2: &BinnedSpectrum) -> Result<f64> {
        compatible(spec1, spec2)?;
        preflight_bins(spec1, spec2)?;
        let sum1 = sparse_dot(spec1, spec1)?;
        let sum2 = sparse_dot(spec2, spec2)?;
        let numerator = sparse_dot(spec1, spec2)?;
        if sum1 * sum2 == 0.0 {
            return Ok(0.0);
        }
        finite_score(numerator / checked_sqrt(sum1 * sum2)?)
    }
}

/// Sum of agreeing intensities for similarity measurement.
///
/// Per bin the score takes the mean of the two intensities minus their absolute
/// difference, `(a + b) / 2 - |a - b|`, and discards the result where it is
/// negative: bins whose intensity difference is larger than their average
/// intensity receive a weight of zero. The retained values are summed and
/// divided by the mean of the two total intensities, so perfect agreement
/// results in a similarity score of `1.0`, and the score is capped at `1.0`.
///
/// Transformation and other factors of the peptide mass spectrometry pairwise
/// peak-list comparison process. Witold E Wolski, Maciej Lalowski, Peter Martus,
/// Ralf Herwig, Patrick Giavalisco, Johan Gobom, Albert Sickmann, Hans Lehrach
/// and Knut Reinert. BMC Bioinformatics 2005, 6:285, doi:10.1186/1471-2105-6-285.
///
/// A bin stored in only one spectrum contributes `(v + 0) / 2 - |v| <= 0`, so it
/// is always discarded; the source still visits it, because the sparse sum of
/// the two bin vectors is their union, and so does this.
///
/// The source registers no parameters, and its protected
/// `precursor_mass_tolerance_` member is unused; see [`BinnedSharedPeakCount`].
#[derive(Clone, Debug, PartialEq)]
pub struct BinnedSumAgreeingIntensities {
    handler: DefaultParamHandler,
}

impl BinnedSumAgreeingIntensities {
    /// Construct the functor with the source's name and its empty parameters.
    ///
    /// # Errors
    ///
    /// As [`BinnedSharedPeakCount::new`].
    pub fn new() -> Result<Self> {
        Ok(Self {
            handler: functor_handler(
                "BinnedSpectrumCompareFunctor",
                "BinnedSumAgreeingIntensities",
            )?,
        })
    }
}

impl BinnedSpectrumCompareFunctor for BinnedSumAgreeingIntensities {
    fn handler(&self) -> &DefaultParamHandler {
        &self.handler
    }

    fn handler_mut(&mut self) -> &mut DefaultParamHandler {
        &mut self.handler
    }

    /// Agreeing intensity summed over the union of the stored bins, over the
    /// mean total intensity, capped at `1.0`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the two spectra do not share a
    /// binning, which is the source's `Exception::IllegalArgument`, when they
    /// together store more than [`MAX_COMPARED_BINS`] bins, and when an `f32`
    /// accumulation overflows.
    ///
    /// A zero mean total intensity - two empty spectra, all-zero bins, or bins
    /// that cancel - yields `Ok(0.0)`. The source divides by that zero and
    /// returns NaN when the numerator is zero too, or `1.0` when it is positive.
    fn score(&self, spec1: &BinnedSpectrum, spec2: &BinnedSpectrum) -> Result<f64> {
        compatible(spec1, spec2)?;
        preflight_bins(spec1, spec2)?;
        let sum1 = sparse_sum(spec1)?;
        let sum2 = sparse_sum(spec2)?;
        // The source builds one sparse vector over the union of both bin sets in
        // f32, truncates its negative coefficients and sums them, all before the
        // result is widened to double. The union walk below is that expression.
        let mut agreeing = 0.0_f32;
        let mut left = spec1.bins.iter();
        let mut right = spec2.bins.iter();
        let mut head_left = left.next().map(|(&index, &value)| (index, value));
        let mut head_right = right.next().map(|(&index, &value)| (index, value));
        loop {
            let (first, second) = match (head_left, head_right) {
                (Some((i, a)), Some((j, b))) if i == j => {
                    head_left = left.next().map(|(&index, &value)| (index, value));
                    head_right = right.next().map(|(&index, &value)| (index, value));
                    (a, b)
                }
                (Some((i, a)), Some((j, _))) if i < j => {
                    head_left = left.next().map(|(&index, &value)| (index, value));
                    (a, 0.0)
                }
                (Some(_), Some((_, b))) => {
                    head_right = right.next().map(|(&index, &value)| (index, value));
                    (0.0, b)
                }
                (Some((_, a)), None) => {
                    head_left = left.next().map(|(&index, &value)| (index, value));
                    (a, 0.0)
                }
                (None, Some((_, b))) => {
                    head_right = right.next().map(|(&index, &value)| (index, value));
                    (0.0, b)
                }
                (None, None) => break,
            };
            let value = (first + second) * 0.5 - (first - second).abs();
            // Eigen's cwiseMax(0) keeps the coefficient unless it is below zero.
            agreeing += if value < 0.0 { 0.0 } else { value };
        }
        if !agreeing.is_finite() {
            return Err(bad("agreeing intensity sum overflows f32"));
        }
        let denominator = (sum1 + sum2) / 2.0;
        if denominator == 0.0 {
            return Ok(0.0);
        }
        finite_score((f64::from(agreeing) / denominator).min(1.0))
    }
}

#[cfg(test)]
mod alignment_budget_tests {
    use super::*;

    #[test]
    fn repeated_alignments_charge_initialization_and_actual_cells() {
        let spectrum = MSSpectrum {
            peaks: vec![crate::Peak1D::new(100.0, 1.0)],
            ..MSSpectrum::default()
        };
        for (tolerance, per_call) in [(Tolerance::Absolute(0.0), 4), (Tolerance::Ppm(1.0), 3)] {
            let alignment = SpectrumAlignment {
                tolerance,
                max_cells: per_call,
            };
            let mut work = AlignmentWork {
                remaining: 2 * per_call,
            };
            for _ in 0..2 {
                assert_eq!(
                    alignment
                        .align_with_work(&spectrum, &spectrum, &mut work)
                        .unwrap(),
                    vec![(0, 0)]
                );
            }
            assert_eq!(work.remaining, 0);
            assert!(
                alignment
                    .align_with_work(&spectrum, &spectrum, &mut work)
                    .unwrap_err()
                    .to_string()
                    .contains("cumulative alignment")
            );
            // The public operation still starts its own configured allowance.
            assert_eq!(alignment.align(&spectrum, &spectrum).unwrap(), vec![(0, 0)]);
        }
    }
}
