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
//! [`BinnedSumAgreeingIntensities`]. The parameter-free functions
//! [`binned_shared_peak_count`], [`binned_cosine`] and
//! [`binned_sum_agreeing_intensities`] are entry points into those same three
//! implementations and cannot disagree with them. Their support documents are
//! `docs/PEAK_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md`,
//! `docs/BINNED_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md`,
//! `docs/BINNED_SHARED_PEAK_COUNT_SUPPORT.md`,
//! `docs/BINNED_SPECTRAL_CONTRAST_ANGLE_SUPPORT.md` and
//! `docs/BINNED_SUM_AGREEING_INTENSITIES_SUPPORT.md`.
//!
//! Four concrete `PeakSpectrumCompareFunctor` derivatives are
//! [`SpectrumPrecursorComparator`] (`COMPARISON/SpectrumPrecursorComparator.h`),
//! [`SpectrumCheapDPCorr`] (`COMPARISON/SpectrumCheapDPCorr.h`),
//! [`PeakAlignment`] (`COMPARISON/PeakAlignment.h`) and
//! [`SpectraSTSimilarityScore`] (`COMPARISON/SpectraSTSimilarityScore.h`), whose
//! support documents are `docs/SPECTRUM_PRECURSOR_COMPARATOR_SUPPORT.md`,
//! `docs/SPECTRUM_CHEAP_DP_CORR_SUPPORT.md`, `docs/PEAK_ALIGNMENT_SUPPORT.md`
//! and `docs/SPECTRAST_SIMILARITY_SCORE_SUPPORT.md`. The remaining three
//! derivatives - [`SpectrumAlignmentScore`], [`ZhangSimilarityScore`] and
//! [`SteinScottImproveScore`] - are still the earlier wave's typed
//! configuration structs and do not implement the trait.

use crate::param::{DefaultParamHandler, Param, ParamValue};
use crate::{Error, MSSpectrum, Peak1D, Precursor, Result};
use std::collections::BTreeMap;
use std::collections::btree_map::Entry;

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
// The three binned scores below have exactly one implementation each, shared by
// a parameter-free function here and by the functor further down that ports the
// corresponding `COMPARISON/Binned*.h` header. The functor is the authoritative
// item - it is what the header maps onto, it carries the `DefaultParamHandler`
// surface and it is reachable through `&dyn BinnedSpectrumCompareFunctor` - and
// the function is the same code path under a shorter name, so the two cannot
// disagree. Before this was collapsed the pair differed in reduction precision,
// denominator grouping, clamping and degenerate-case policy.

/// Cosine of the spectral contrast angle of two binned spectra: the
/// parameter-free entry point to [`BinnedSpectralContrastAngle`].
///
/// Despite the upstream class name this is the cosine and not the angle,
/// `dot(a, b) / sqrt(dot(a, a) * dot(b, b))`, reduced in `f32` exactly as
/// `Eigen::SparseVector<float>::dot` reduces it, with the source's denominator
/// grouping and no clamping. Signed bins are supported, and a zero
/// `dot(a, a) * dot(b, b)` yields the source's defined `0.0`.
///
/// This is the same implementation as
/// [`BinnedSpectralContrastAngle`]'s [`BinnedSpectrumCompareFunctor::score`],
/// not a second one, so the two agree bit for bit. Construct the functor when a
/// trait object or the parameter surface is wanted; call this when neither is.
///
/// # Errors
///
/// As the functor: [`Error::InvalidValue`] for incompatible binning, for more
/// than [`MAX_COMPARED_BINS`] combined stored bins, and for an `f32` dot product
/// that overflows.
pub fn binned_cosine(a: &BinnedSpectrum, b: &BinnedSpectrum) -> Result<f64> {
    spectral_contrast_angle(a, b)
}

/// Stored-bin intersection over the larger stored-bin count: the parameter-free
/// entry point to [`BinnedSharedPeakCount`].
///
/// Explicit zeros still count as stored bins, as `Eigen::SparseVector::nonZeros`
/// counts them, and two spectra that store no bins yield `0.0` where the source
/// divides by zero.
///
/// This is the same implementation as [`BinnedSharedPeakCount`]'s
/// [`BinnedSpectrumCompareFunctor::score`], not a second one.
///
/// # Errors
///
/// As the functor: [`Error::InvalidValue`] for incompatible binning and for more
/// than [`MAX_COMPARED_BINS`] combined stored bins.
pub fn binned_shared_peak_count(a: &BinnedSpectrum, b: &BinnedSpectrum) -> Result<f64> {
    shared_peak_count(a, b)
}

/// Sum of `max(0, (a + b) / 2 - |a - b|)` over the union of the stored bins,
/// divided by the mean total intensity and capped at one: the parameter-free
/// entry point to [`BinnedSumAgreeingIntensities`].
///
/// Negative bins are accepted and truncated away, as upstream; a zero mean total
/// intensity yields `0.0` where the source divides by zero.
///
/// This is the same implementation as [`BinnedSumAgreeingIntensities`]'s
/// [`BinnedSpectrumCompareFunctor::score`], not a second one.
///
/// # Errors
///
/// As the functor: [`Error::InvalidValue`] for incompatible binning, for more
/// than [`MAX_COMPARED_BINS`] combined stored bins, and for an `f32`
/// accumulation that overflows.
pub fn binned_sum_agreeing_intensities(a: &BinnedSpectrum, b: &BinnedSpectrum) -> Result<f64> {
    sum_agreeing_intensities(a, b)
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
/// The base has **seven** derivatives in the pinned tree, not the six that
/// `PeakSpectrumCompareFunctor.cpp:11-16` `#include`s for factory registration:
/// `SpectrumCheapDPCorr`, `SpectrumPrecursorComparator`, `ZhangSimilarityScore`,
/// `SpectrumAlignmentScore`, `SteinScottImproveScore`, `PeakAlignment` and -
/// absent from that include list - `SpectraSTSimilarityScore`. All seven
/// override the one-spectrum overload as `operator()(spec, spec)`, so the
/// default [`self_score`](Self::self_score) is that delegation and an
/// implementor only overrides it to record a cheaper closed form.
///
/// **Four of the seven derivatives implement this trait:**
/// [`SpectrumPrecursorComparator`], [`SpectrumCheapDPCorr`], [`PeakAlignment`]
/// and [`SpectraSTSimilarityScore`], each carrying the parameter tree its own
/// header registers. The remaining three - [`SpectrumAlignmentScore`],
/// [`ZhangSimilarityScore`] and [`SteinScottImproveScore`] - were ported in an
/// earlier wave as typed `Copy` configuration structs with no
/// [`DefaultParamHandler`]; giving them one means porting the parameter tree
/// each of them registers upstream, which is their own headers' work, and a
/// handler that did not carry those parameters would make `set_parameters`
/// silently ineffective. See
/// `docs/PEAK_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md`. Its sibling
/// [`BinnedSpectrumCompareFunctor`] has all three of its source derivatives
/// shipped here.
///
/// Only [`PeakAlignment`] never renames its handler, so a
/// `&dyn PeakSpectrumCompareFunctor` over the four reports three derived names
/// and one `"PeakSpectrumCompareFunctor"`.
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
///
/// Eigen's sparse dot is a scalar merge of the two inner iterators, so the
/// association reproduced here is the one it uses; see [`sparse_sum`] for the
/// reductions where that is not true, and note that a C++ build which contracts
/// `res += a * b` into an FMA (the default at `-ffp-contract=fast`) rounds once
/// where this rounds twice.
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
/// `f32`, here sequentially in ascending index order.
///
/// **This one is not bit-exact against every C++ build, and the fidelity claim
/// is bounded accordingly.** `SparseVector::sum()` maps the stored-value array
/// to a dense vector and calls the dense reduction, which Eigen vectorises into
/// several packet accumulators and combines at the end - an association that
/// depends on the target's packet width and on the Eigen version, and that is
/// not the sequential one. The same is true of the agreeing-intensity numerator,
/// which upstream is `s.coeffs().cwiseMax(0).sum()` over that same dense value
/// array. The port sums the same `f32` values in the same index order; what it
/// reproduces is the `f32` *precision* of the reduction, and agreement with a
/// vectorised C++ build is to `f32` reduction rounding rather than bit for bit.
/// A sequential order is chosen because it is deterministic, is what a scalar
/// build produces, and is the association a parallel implementation must
/// reproduce under `concept::parallel`.
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
///
/// [`binned_shared_peak_count`] is the parameter-free entry point to this same
/// computation.
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
        shared_peak_count(spec1, spec2)
    }
}

/// `BinnedSharedPeakCount::operator()`, the one implementation behind both that
/// functor and [`binned_shared_peak_count`].
fn shared_peak_count(spec1: &BinnedSpectrum, spec2: &BinnedSpectrum) -> Result<f64> {
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
/// [`binned_cosine`] is the parameter-free entry point to this same
/// computation. It used to be a second, `f64` implementation with a
/// `sqrt(sum1) * sqrt(sum2)` denominator and a clamp to `[-1, 1]`, which
/// disagreed with this one in the last bits; it is now this one.
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
        spectral_contrast_angle(spec1, spec2)
    }
}

/// `BinnedSpectralContrastAngle::operator()`, the one implementation behind both
/// that functor and [`binned_cosine`].
fn spectral_contrast_angle(spec1: &BinnedSpectrum, spec2: &BinnedSpectrum) -> Result<f64> {
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
///
/// [`binned_sum_agreeing_intensities`] is the parameter-free entry point to this
/// same computation. It used to be a second, `f64` implementation that rejected
/// negative bins, which the source accepts and truncates away; it is now this
/// one, so a negative bin is no longer an error from either.
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
        sum_agreeing_intensities(spec1, spec2)
    }
}

/// `BinnedSumAgreeingIntensities::operator()`, the one implementation behind both
/// that functor and [`binned_sum_agreeing_intensities`].
fn sum_agreeing_intensities(spec1: &BinnedSpectrum, spec2: &BinnedSpectrum) -> Result<f64> {
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

// ---------------------------------------------------------------------------
// The four remaining derivatives of the peak-spectrum functor hierarchy:
//   COMPARISON/SpectrumPrecursorComparator.h
//   COMPARISON/SpectrumCheapDPCorr.h
//   COMPARISON/PeakAlignment.h
//   COMPARISON/SpectraSTSimilarityScore.h
// ---------------------------------------------------------------------------

/// `boost::math::constants::root_two_pi<double>()`, the divisor in Boost's
/// normal density.
///
/// Transcribed from Boost's own decimal literal rather than computed as
/// `(2.0 * PI).sqrt()`: the two need not agree in the last bit, and
/// [`SpectrumCheapDPCorr`] inherits this constant's rounding on every matched
/// peak pair.
const ROOT_TWO_PI: f64 = 2.506628274631000502415765284811045253e0;

/// `boost::math::pdf(boost::math::normal_distribution<double>(0, sd), x)`,
/// reproduced statement by statement so that the rounding order matches.
///
/// Boost computes `exponent = x - mean`, `exponent *= -exponent`,
/// `exponent /= 2 * sd * sd`, `result = exp(exponent)` and finally
/// `result /= sd * root_two_pi`. Mean is fixed at zero here because the only
/// caller, `SpectrumCheapDPCorr::comparepeaks_`, constructs
/// `normal_distribution<double>(0., variation)`.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `sd` is not finite and positive, which
/// is Boost's `check_scale` domain error under the default policy, and when `x`
/// is NaN, which is its `check_x` domain error. An infinite `x` yields `0.0`,
/// as Boost's explicit early return does.
fn normal_pdf(sd: f64, x: f64) -> Result<f64> {
    if !sd.is_finite() || sd <= 0.0 {
        return Err(bad(
            "the Gaussian match term needs a finite positive scale; \
             a zero or negative variation makes it undefined",
        ));
    }
    if x.is_infinite() {
        return Ok(0.0);
    }
    if x.is_nan() {
        return Err(bad(
            "the Gaussian match term needs a finite position difference",
        ));
    }
    let mut exponent = x;
    exponent *= -exponent;
    exponent /= 2.0 * sd * sd;
    let mut result = exponent.exp();
    result /= sd * ROOT_TWO_PI;
    Ok(result)
}

/// A zeroed buffer of `cells` values whose allocation failure is an error
/// rather than a process abort.
///
/// Every caller checks `cells` against an explicit ceiling first, so this is
/// the second line of defence and not the bound itself.
fn zeroed<T: Copy + Default>(cells: usize) -> Result<Vec<T>> {
    let mut buffer = Vec::new();
    buffer
        .try_reserve_exact(cells)
        .map_err(|_| bad("comparison matrix allocation failed"))?;
    buffer.resize(cells, T::default());
    Ok(buffer)
}

/// Narrow a consensus coordinate computed in `f64` to the `f32` a peak stores,
/// refusing a value that does not survive the narrowing.
///
/// The source assigns the `double` expression straight into
/// `Peak1D::setIntensity`, where an out-of-range value is undefined behaviour;
/// here it is a checked error.
fn narrow(value: f64) -> Result<f32> {
    let result = value as f32;
    if result.is_finite() {
        Ok(result)
    } else {
        Err(bad("consensus intensity does not fit f32"))
    }
}

/// Build the parameter surface of a concrete `PeakSpectrumCompareFunctor`
/// derivative: the base names the handler after itself, the derivative renames
/// it, registers its defaults and copies them into the current parameters.
///
/// `name` is `None` for the one derivative that never calls `setName`.
fn derived_handler(name: Option<&str>, defaults: Param) -> Result<DefaultParamHandler> {
    let mut handler = DefaultParamHandler::new("PeakSpectrumCompareFunctor")?;
    if let Some(name) = name {
        handler.set_name(name)?;
    }
    handler.set_defaults(defaults)?;
    handler.defaults_to_parameters()?;
    Ok(handler)
}

/// Read a parameter as a plain `f64`, as the source's `(double)param_.getValue`.
fn float_parameter(handler: &DefaultParamHandler, key: &str) -> Result<f64> {
    let value = handler.parameters().value(key)?.to_f64()?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad("comparison parameters must be finite"))
    }
}

/// Read a parameter as a nonnegative integer.
///
/// The source casts these to `UInt` or `unsigned int`, so a negative value
/// silently becomes an enormous positive one; here it is refused.
fn count_parameter(handler: &DefaultParamHandler, key: &str) -> Result<u32> {
    let value = handler.parameters().value(key)?.to_i64()?;
    u32::try_from(value)
        .map_err(|_| bad("comparison count parameters must fit an unsigned 32-bit integer"))
}

/// Compare just the parent mass of two spectra.
///
/// The score is `window - |Δ precursor m/z|`, clamped to zero: the source
/// returns `0` when the distance exceeds `window` and the difference otherwise,
/// which is the same function without the redundant second subtraction. Only
/// the **first** precursor of each spectrum is read, and a spectrum with no
/// precursor contributes m/z `0`, so two spectra that both lack a precursor
/// score the full `window`. That convention is the source's and is preserved.
///
/// The single parameter `window` is registered as the integer `2` with the
/// description "Allowed deviation between precursor peaks.", exactly as
/// `SpectrumPrecursorComparator.cpp:22` registers it, and is read through the
/// handler on every call as the source reads `param_.getValue("window")`.
///
/// See `docs/SPECTRUM_PRECURSOR_COMPARATOR_SUPPORT.md`.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectrumPrecursorComparator {
    handler: DefaultParamHandler,
}

impl SpectrumPrecursorComparator {
    /// Construct the functor with the source's name and its one default.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] only if the parameter handler rejects the
    /// fixed name or the fixed one-entry default tree, which cannot happen for
    /// these literals.
    pub fn new() -> Result<Self> {
        let mut defaults = Param::new();
        defaults.set_value(
            "window",
            ParamValue::Integer(2),
            "Allowed deviation between precursor peaks.",
            &[],
        )?;
        Ok(Self {
            handler: derived_handler(Some("SpectrumPrecursorComparator"), defaults)?,
        })
    }

    /// Current `window` parameter, in Thomson.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the parameter has been replaced by a
    /// value that is not a finite number.
    pub fn window(&self) -> Result<f64> {
        float_parameter(&self.handler, "window")
    }
}

impl Default for SpectrumPrecursorComparator {
    /// The source's default construction.
    ///
    /// [`SpectrumPrecursorComparator::new`] is fallible only through the
    /// parameter handler's resource limits, which a two-word name and a
    /// one-entry default tree cannot reach, so this cannot fail in practice.
    fn default() -> Self {
        Self::new().expect("the fixed one-entry parameter tree is within the handler's limits")
    }
}

impl PeakSpectrumCompareFunctor for SpectrumPrecursorComparator {
    fn handler(&self) -> &DefaultParamHandler {
        &self.handler
    }

    fn handler_mut(&mut self) -> &mut DefaultParamHandler {
        &mut self.handler
    }

    /// Precursor similarity of `a` and `b`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an invalid spectrum and for a
    /// `window` that is not a finite nonnegative number. The source reads the
    /// parameter unchecked, so a negative window there yields a negative score
    /// for every pair; here it is refused.
    fn score(&self, a: &MSSpectrum, b: &MSSpectrum) -> Result<f64> {
        a.validate()?;
        b.validate()?;
        let window = self.window()?;
        if window < 0.0 {
            return Err(bad("precursor window must be finite and nonnegative"));
        }
        let left = a.precursors.first().map_or(0.0, |p| p.mz);
        let right = b.precursors.first().map_or(0.0, |p| p.mz);
        Ok((window - (left - right).abs()).max(0.0))
    }
}

/// Largest number of dynamic-programming cells one [`SpectrumCheapDPCorr`]
/// comparison may allocate, summed over every block the scan hands to the
/// inner alignment.
///
/// The source has no such ceiling: `dynprog_` allocates
/// `(xrun + 1) * (yrun + 1)` doubles **and** as many `int`s for every pairable
/// run, and a `variation` near its documented maximum of `1` makes one run span
/// both whole spectra. The budget is charged before either buffer is allocated,
/// so a refusal leaves both inputs and the functor untouched.
pub const MAX_DP_CORR_CELLS: usize = 1_000_000;

/// Optimal alignment of two stick spectra by dynamic programming, with a
/// Gaussian-weighted match term.
///
/// The scan walks both peak lists at once. Peaks further apart than
/// `variation` percent of their mean m/z cannot pair and are consumed one at a
/// time; where several peaks on both sides could pair, the run is handed to an
/// `O(n*m)` alignment and only pairs that could score above zero are ever
/// examined. That is the "cheap" in the class name.
///
/// Three parameters are registered, with the source's own defaults and
/// descriptions: `variation` (`0.001`), `int_cnt` (`0`) and `keeppeaks` (`0`).
/// `int_cnt` selects how the two peak heights enter the score - `0` their
/// product, `1` the square root of their product, `2` their sum and `3` their
/// agreeing intensity `max(0, (i1 + i2) / 2 - |i1 - i2|)`.
///
/// # Stateful accessors
///
/// The source's `operator()` is `const` but writes three `mutable` members: the
/// consensus spectrum, the peak map and the weighting factor. That cannot be
/// expressed behind [`PeakSpectrumCompareFunctor::score`], which really is
/// read-only here, so the port splits them:
///
/// * [`score`](PeakSpectrumCompareFunctor::score) returns the number and
///   discards the consensus. The number is unaffected, because neither
///   `factor_` nor the consensus feeds back into the score.
/// * [`compare`](Self::compare) returns the same number and records the
///   consensus and the peak map, which [`last_consensus`](Self::last_consensus)
///   and [`peak_map`](Self::peak_map) then expose - the source's
///   `lastconsensus()` and `getPeakMap()`. It also resets the factor to `0.5`,
///   as the last statement of the source's `operator()` does.
///
/// See `docs/SPECTRUM_CHEAP_DP_CORR_SUPPORT.md`.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectrumCheapDPCorr {
    handler: DefaultParamHandler,
    factor: f64,
    last_consensus: MSSpectrum,
    peak_map: BTreeMap<usize, usize>,
}

impl SpectrumCheapDPCorr {
    /// Construct the functor with the source's name, defaults and factor.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] only if the parameter handler rejects the
    /// fixed name or the fixed three-entry default tree.
    pub fn new() -> Result<Self> {
        let mut defaults = Param::new();
        defaults.set_value(
            "variation",
            ParamValue::Float(0.001),
            "Maximum difference in position (in percent of the current m/z).\n\
             Note that big values of variation ( 1 being the maximum ) result in \
             consideration of all possible pairings which has a running time of O(n*n)",
            &[],
        )?;
        defaults.set_value(
            "int_cnt",
            ParamValue::Integer(0),
            "How the peak heights are used in the score.\n\
             0 = product\n1 = sqrt(product)\n2 = sum\n3 = agreeing intensity\n",
            &[],
        )?;
        defaults.set_value(
            "keeppeaks",
            ParamValue::Integer(0),
            "Flag that states if peaks without alignment partner are kept in the consensus spectrum.",
            &[],
        )?;
        Ok(Self {
            handler: derived_handler(Some("SpectrumCheapDPCorr"), defaults)?,
            factor: 0.5,
            last_consensus: MSSpectrum::default(),
            peak_map: BTreeMap::new(),
        })
    }

    /// Weight given to the second spectrum when the next [`compare`](Self::compare)
    /// builds its consensus, as the source's `factor_`.
    pub fn factor(&self) -> f64 {
        self.factor
    }

    /// Set the weighting of the second spectrum for the next consensus, as
    /// `setFactor`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] unless `0 < factor < 1`, reproducing the
    /// source's `Exception::OutOfRange`; both bounds are exclusive there.
    pub fn set_factor(&mut self, factor: f64) -> Result<()> {
        if factor < 1.0 && factor > 0.0 {
            self.factor = factor;
            Ok(())
        } else {
            Err(Error::InvalidRange(
                "the consensus weighting factor must lie strictly between 0 and 1".into(),
            ))
        }
    }

    /// Consensus spectrum of the last [`compare`](Self::compare), as
    /// `lastconsensus()`.
    ///
    /// Before the first call this is an empty default spectrum; the source's is
    /// a default-constructed `PeakSpectrum` too. Its single precursor carries
    /// the mean of the two input precursor m/z values and the **first**
    /// spectrum's charge, which is what the source writes.
    pub fn last_consensus(&self) -> &MSSpectrum {
        &self.last_consensus
    }

    /// Peak indices of the first spectrum mapped to their partner in the second,
    /// as `getPeakMap()`.
    ///
    /// The source's `std::map<UInt, UInt>` is a [`BTreeMap`] here, so iteration
    /// order still ascends by key. Only aligned pairs appear.
    pub fn peak_map(&self) -> &BTreeMap<usize, usize> {
        &self.peak_map
    }

    /// Score `x` against `y` and record the consensus spectrum and peak map.
    ///
    /// This is the source's `operator()(x, y)` including its three `mutable`
    /// side effects, ending with the reset of the weighting factor to `0.5`.
    /// [`PeakSpectrumCompareFunctor::score`] computes the same number without
    /// them.
    ///
    /// # Errors
    ///
    /// As [`PeakSpectrumCompareFunctor::score`]. A failure leaves the recorded
    /// consensus, peak map and factor exactly as they were.
    pub fn compare(&mut self, x: &MSSpectrum, y: &MSSpectrum) -> Result<f64> {
        let outcome = self.run(x, y)?;
        self.last_consensus = outcome.consensus;
        self.peak_map = outcome.peak_map;
        self.factor = 0.5;
        Ok(outcome.score)
    }

    /// The source's `operator()(x, y)` with its outputs returned rather than
    /// written through `mutable` members.
    fn run(&self, x: &MSSpectrum, y: &MSSpectrum) -> Result<CheapDpOutcome> {
        validate_spectrum(x, true)?;
        validate_spectrum(y, true)?;
        x.len()
            .checked_add(y.len())
            .ok_or_else(|| bad("combined peak count overflows"))?;
        let variation_fraction = float_parameter(&self.handler, "variation")?;
        if variation_fraction <= 0.0 {
            return Err(bad(
                "variation must be positive; the Gaussian match term has no zero-width limit",
            ));
        }
        let int_cnt = count_parameter(&self.handler, "int_cnt")?;
        let keep_peaks = self.handler.parameters().value("keeppeaks")?.to_i64()? != 0;

        let left_precursor = x.precursors.first().cloned().unwrap_or_default();
        let right_precursor = y.precursors.first().cloned().unwrap_or_default();
        let mut consensus = MSSpectrum {
            precursors: vec![Precursor::new(
                (left_precursor.mz + right_precursor.mz) / 2.0,
                left_precursor.charge,
            )],
            ..MSSpectrum::default()
        };
        let mut peak_map = BTreeMap::new();
        let mut run = CheapDpRun {
            variation_fraction,
            int_cnt,
            keep_peaks,
            factor: self.factor,
            budget: MAX_DP_CORR_CELLS,
            consensus: &mut consensus,
            peak_map: &mut peak_map,
        };

        let (px, py) = (x.peaks.as_slice(), y.peaks.as_slice());
        let mut score = 0.0;
        let mut xi = 0;
        let mut yi = 0;
        while xi < px.len() && yi < py.len() {
            let variation = (px[xi].mz + py[yi].mz) / 2.0 * variation_fraction;
            if (px[xi].mz - py[yi].mz).abs() > variation {
                if px[xi].mz < py[yi].mz {
                    if keep_peaks {
                        let intensity = narrow(f64::from(px[xi].intensity) * (1.0 - run.factor))?;
                        run.consensus.peaks.push(Peak1D::new(px[xi].mz, intensity));
                    }
                    xi += 1;
                } else {
                    if keep_peaks {
                        let intensity = narrow(f64::from(py[yi].intensity) * run.factor)?;
                        run.consensus.peaks.push(Peak1D::new(py[yi].mz, intensity));
                    }
                    yi += 1;
                }
                continue;
            }
            let (xrun, yrun) = pairable_run(px, py, xi, yi, variation);
            if xrun > 1 && yrun > 1 {
                score += run.dynamic_program(px, py, xi, xi + xrun - 1, yi, yi + yrun - 1)?;
                xi += xrun;
                yi += yrun;
            } else {
                // The source's one-to-one consensus weights the FIRST spectrum
                // by (1 - factor); the traceback inside dynprog_ weights the
                // SECOND one by (1 - factor) instead. Both are reproduced.
                let mz = px[xi].mz * (1.0 - run.factor) + py[yi].mz * run.factor;
                let intensity = narrow(
                    f64::from(px[xi].intensity) * (1.0 - run.factor)
                        + f64::from(py[yi].intensity) * run.factor,
                )?;
                run.consensus.peaks.push(Peak1D::new(mz, intensity));
                // The source's else branch here compares the two indices rather
                // than the stored value, unlike the otherwise identical code in
                // dynprog_. It is unreachable either way: the map is cleared per
                // call and this scan visits each index of the first spectrum at
                // most once, so only the insert can run.
                run.peak_map.entry(xi).or_insert(yi);
                score += run.compare_peaks(
                    px[xi].mz,
                    py[yi].mz,
                    f64::from(px[xi].intensity),
                    f64::from(py[yi].intensity),
                )?;
                xi += 1;
                yi += 1;
            }
        }
        Ok(CheapDpOutcome {
            score: finite_score(score)?,
            consensus,
            peak_map,
        })
    }
}

impl PeakSpectrumCompareFunctor for SpectrumCheapDPCorr {
    fn handler(&self) -> &DefaultParamHandler {
        &self.handler
    }

    fn handler_mut(&mut self) -> &mut DefaultParamHandler {
        &mut self.handler
    }

    /// Optimal-alignment correlation of `a` and `b`.
    ///
    /// Identical in value to [`compare`](SpectrumCheapDPCorr::compare), which
    /// additionally records the consensus spectrum and the peak map.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] for unsorted peaks, which the source
    /// silently mis-aligns, and [`Error::InvalidValue`] for a negative
    /// intensity, a `variation` that is not positive, an `int_cnt` outside
    /// `0..=3` - where the source returns `-1` behind a `// TODO exception` -
    /// for more than [`MAX_DP_CORR_CELLS`] dynamic-programming cells, and for a
    /// non-finite score.
    fn score(&self, a: &MSSpectrum, b: &MSSpectrum) -> Result<f64> {
        Ok(self.run(a, b)?.score)
    }
}

/// What one `SpectrumCheapDPCorr::operator()` produces: the score plus the two
/// values the source publishes through `mutable` members.
struct CheapDpOutcome {
    score: f64,
    consensus: MSSpectrum,
    peak_map: BTreeMap<usize, usize>,
}

/// The parameters and accumulators `SpectrumCheapDPCorr::operator()` shares with
/// its `dynprog_` helper.
struct CheapDpRun<'a> {
    variation_fraction: f64,
    int_cnt: u32,
    keep_peaks: bool,
    factor: f64,
    budget: usize,
    consensus: &'a mut MSSpectrum,
    peak_map: &'a mut BTreeMap<usize, usize>,
}

impl CheapDpRun<'_> {
    /// `SpectrumCheapDPCorr::comparepeaks_`: a Gaussian in the position
    /// difference, whose standard deviation is `variation` percent of the mean
    /// of the two m/z values, times an intensity term chosen by `int_cnt`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the Gaussian scale is not positive,
    /// when `int_cnt` is outside `0..=3` - the source returns `-1` there, which
    /// a caller summing scores cannot distinguish from a real contribution -
    /// and when `int_cnt` is `1` and the intensity product is negative, where
    /// the source's `sqrt` yields NaN.
    fn compare_peaks(&self, posa: f64, posb: f64, inta: f64, intb: f64) -> Result<f64> {
        let variation = (posa + posb) / 2.0 * self.variation_fraction;
        let density = normal_pdf(variation, posa - posb)?;
        match self.int_cnt {
            0 => Ok(density * inta * intb),
            1 => Ok(density * checked_sqrt(inta * intb)?),
            2 => Ok(density * (inta + intb)),
            3 => Ok((density * ((inta + intb) / 2.0 - (inta - intb).abs())).max(0.0)),
            _ => Err(bad("int_cnt must be 0, 1, 2 or 3")),
        }
    }

    /// `SpectrumCheapDPCorr::dynprog_`: the optimal pairing of one run of peaks
    /// from each spectrum, plus the consensus peaks and peak-map entries its
    /// traceback emits.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the block exceeds the remaining
    /// [`MAX_DP_CORR_CELLS`] budget, when its buffers cannot be allocated, when
    /// a consensus intensity does not fit `f32`, and for anything
    /// [`compare_peaks`](Self::compare_peaks) refuses.
    fn dynamic_program(
        &mut self,
        x: &[Peak1D],
        y: &[Peak1D],
        xstart: usize,
        xend: usize,
        ystart: usize,
        yend: usize,
    ) -> Result<f64> {
        let rows = xend - xstart + 2;
        let cols = yend - ystart + 2;
        let cells = rows
            .checked_mul(cols)
            .ok_or_else(|| bad("dynamic-programming block size overflows"))?;
        self.budget = self
            .budget
            .checked_sub(cells)
            .ok_or_else(|| bad("comparison exceeds the dynamic-programming cell ceiling"))?;
        let mut dp = zeroed::<f64>(cells)?;
        let mut trace = zeroed::<i8>(cells)?;
        for i in 1..rows {
            for j in 1..cols {
                let left_peak = x[xstart + i - 1];
                let right_peak = y[ystart + j - 1];
                // The source spells this sum (y + x) here and (x + y) in the
                // scan; floating addition is commutative, so it is one value.
                let variation = (right_peak.mz + left_peak.mz) / 2.0 * self.variation_fraction;
                let align = if (left_peak.mz - right_peak.mz).abs() > variation {
                    0.0
                } else {
                    self.compare_peaks(
                        left_peak.mz,
                        right_peak.mz,
                        f64::from(left_peak.intensity),
                        f64::from(right_peak.intensity),
                    )?
                };
                let from_left = dp[i * cols + j - 1];
                let from_diagonal = dp[(i - 1) * cols + j - 1] + align;
                let from_above = dp[(i - 1) * cols + j];
                // Source: ((left > diagonal) ? left : diagonal) > above, then a
                // second strict comparison of diagonal against left. Ties
                // therefore prefer "above", and then "left" over "diagonal".
                let best_of_two = if from_left > from_diagonal {
                    from_left
                } else {
                    from_diagonal
                };
                if best_of_two > from_above {
                    if from_diagonal > from_left {
                        dp[i * cols + j] = from_diagonal;
                        trace[i * cols + j] = 5;
                    } else {
                        dp[i * cols + j] = from_left;
                        trace[i * cols + j] = -1;
                    }
                } else {
                    dp[i * cols + j] = from_above;
                    trace[i * cols + j] = 1;
                }
            }
        }

        let mut i = xend - xstart + 1;
        let mut j = yend - ystart + 1;
        loop {
            match trace[i * cols + j] {
                5 => {
                    let left_peak = x[xstart + i - 1];
                    let right_peak = y[ystart + j - 1];
                    let mz = right_peak.mz * (1.0 - self.factor) + left_peak.mz * self.factor;
                    let intensity = narrow(
                        f64::from(right_peak.intensity) * (1.0 - self.factor)
                            + f64::from(left_peak.intensity) * self.factor,
                    )?;
                    self.consensus.peaks.push(Peak1D::new(mz, intensity));
                    match self.peak_map.entry(xstart + i - 1) {
                        Entry::Vacant(slot) => {
                            slot.insert(ystart + j - 1);
                        }
                        // Unreachable: the traceback decrements i on every 5, so
                        // each key is written at most once per call.
                        Entry::Occupied(mut slot) => {
                            let previous = *slot.get();
                            slot.insert((ystart + j - 1).min(previous));
                        }
                    }
                    i -= 1;
                    j -= 1;
                }
                1 => {
                    if self.keep_peaks {
                        let peak = x[xstart + i - 1];
                        let intensity = narrow(f64::from(peak.intensity) * (1.0 - self.factor))?;
                        self.consensus.peaks.push(Peak1D::new(peak.mz, intensity));
                    }
                    i -= 1;
                }
                -1 => {
                    if self.keep_peaks {
                        let peak = y[ystart + j - 1];
                        let intensity = narrow(f64::from(peak.intensity) * self.factor)?;
                        self.consensus.peaks.push(Peak1D::new(peak.mz, intensity));
                    }
                    j -= 1;
                }
                // Every cell with i >= 1 and j >= 1 was written above, so this
                // cannot happen; the source would spin forever instead.
                _ => return Err(bad("dynamic-programming traceback reached an unset cell")),
            }
            if i == 0 || j == 0 {
                break;
            }
        }
        Ok(dp[(xend - xstart + 1) * cols + (yend - ystart + 1)])
    }
}

/// How many peaks of each spectrum, starting at `xi` and `yi`, could pair with
/// one another - the source's `xrun`/`yrun` loop.
///
/// The source spells the bounds `xit + xrun != x.end()`; both counters only ever
/// grow by one and the body breaks the moment either reaches the end, so `<` is
/// the same condition and makes the indexing obviously in range.
fn pairable_run(
    x: &[Peak1D],
    y: &[Peak1D],
    xi: usize,
    yi: usize,
    variation: f64,
) -> (usize, usize) {
    let mut xrun = 1;
    let mut yrun = 1;
    // The source writes the two disjunction terms as !(a < b); every m/z here
    // is finite, so that is exactly a >= b and is spelled so.
    while xi + xrun < x.len()
        && yi + yrun < y.len()
        && (x[xi + xrun - 1].mz + variation >= y[yi + yrun].mz
            || y[yi + yrun - 1].mz + variation >= x[xi + xrun].mz)
    {
        if y[yi + yrun - 1].mz + variation > x[xi + xrun].mz {
            xrun += 1;
        } else if x[xi + xrun - 1].mz + variation > y[yi + yrun].mz {
            yrun += 1;
        } else {
            xrun += 1;
            yrun += 1;
        }
        if xi + xrun == x.len() || yi + yrun == y.len() {
            break;
        }
    }
    (xrun, yrun)
}

/// Largest number of alignment-matrix cells one [`PeakAlignment`] comparison
/// may allocate, counting the `(n + 1) * (m + 1)` score matrix.
///
/// The source allocates that matrix plus, in `getAlignmentTraceback`, an
/// `n * m` direction matrix, with no ceiling at all: two 5000-peak spectra ask
/// it for 200 MB. The budget is checked before either buffer is allocated.
pub const MAX_ALIGNMENT_MATRIX_CELLS: usize = 4_000_000;

/// Global alignment of two peak lists with a constant gap cost.
///
/// The class comment calls this Needleman-Wunsch; the recurrence is indeed
/// global, with the first row and column pre-charged with multiples of the gap
/// cost, but the reported score is the best cell of the **last row or last
/// column** rather than the corner, so a suffix of either spectrum may be left
/// unaligned. A pair of peaks may only align when their m/z differ by at most
/// `epsilon`; otherwise the cell can only come from a gap. The gap cost is
/// `epsilon` as well - `PeakAlignment.cpp:99` sets it from the same parameter
/// under a `//TODO gapcost dependence on distance ?`.
///
/// Four parameters are registered with the source's defaults and descriptions:
/// `epsilon` (`0.2`), `normalized` (`1`), `heuristic_level` (`0`) and
/// `precursor_mass_tolerance` (`3.0`). **`normalized` is never read.** The
/// score is divided by the geometric mean of the two self-alignment scores
/// unconditionally, so clearing the flag changes nothing; it is registered here
/// so the parameter surface matches, and the port reads it no more than the
/// source does.
///
/// Two shortcuts precede the alignment, in this order: precursors further apart
/// than `precursor_mass_tolerance` score `0`, and - when `heuristic_level` is
/// nonzero - so do two spectra whose `heuristic_level` most intense peaks share
/// no m/z within `epsilon`.
///
/// This functor is the one derivative that never calls `setName`, so
/// [`name`](PeakSpectrumCompareFunctor::name) reports the base's
/// `"PeakSpectrumCompareFunctor"`. That is the source's behaviour, not an
/// omission here.
///
/// See `docs/PEAK_ALIGNMENT_SUPPORT.md`.
#[derive(Clone, Debug, PartialEq)]
pub struct PeakAlignment {
    handler: DefaultParamHandler,
}

impl PeakAlignment {
    /// Construct the functor with the source's four defaults.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] only if the parameter handler rejects the
    /// fixed four-entry default tree.
    pub fn new() -> Result<Self> {
        let mut defaults = Param::new();
        defaults.set_value(
            "epsilon",
            ParamValue::Float(0.2),
            "defines the absolute error of the mass spectrometer",
            &[],
        )?;
        defaults.set_value(
            "normalized",
            ParamValue::Integer(1),
            "is set 1 if the similarity-measurement is normalized to the range [0,1]",
            &[],
        )?;
        defaults.set_value(
            "heuristic_level",
            ParamValue::Integer(0),
            "set 0 means no heuristic is applied otherwise the given value is interpreted as \
             unsigned integer, the number of strongest peaks considered for heurisitcs - in \
             those sets of peaks has to be at least one match to conduct comparison",
            &[],
        )?;
        defaults.set_value(
            "precursor_mass_tolerance",
            ParamValue::Float(3.0),
            "Mass tolerance of the precursor peak, defines the distance of two PrecursorPeaks \
             for which they are supposed to be from different peptides",
            &[],
        )?;
        // The source never renames the handler, so the base name survives.
        Ok(Self {
            handler: derived_handler(None, defaults)?,
        })
    }

    /// Aligned `(index in spec1, index in spec2)` pairs, ascending, as
    /// `getAlignmentTraceback`.
    ///
    /// Only diagonal steps - actually aligned peaks - are reported; gap steps
    /// are traversed silently. The traceback starts at the best cell of the last
    /// row or column, preferring the earliest such cell in the row scan and
    /// overriding it only on a strict improvement in the column scan.
    ///
    /// Ties inside the matrix resolve the way the source's zero-filled
    /// direction matrix does: when no single predecessor is strictly best the
    /// cell keeps its initial `0`, which the traceback reads as "from the left",
    /// so a tie consumes a peak of the second spectrum. That is the case the
    /// source marks `// TODO the cases where all or two values are equal`, and
    /// it is reproduced rather than repaired because the reported alignment is
    /// observable.
    ///
    /// Unlike [`score`](PeakSpectrumCompareFunctor::score), this entry point
    /// applies **neither** shortcut and does **not** guard a zero variance, so
    /// it is the raw alignment of the two peak lists.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] for unsorted peaks and
    /// [`Error::InvalidValue`] for an empty spectrum - where the source divides
    /// by a zero pair count and carries NaN through the whole matrix - for a
    /// zero peak-distance variance, where the source divides by a zero sigma and
    /// returns infinities, for a negative intensity under the score's `sqrt`,
    /// for a non-finite `epsilon`, and for more than
    /// [`MAX_ALIGNMENT_MATRIX_CELLS`] matrix cells.
    pub fn alignment_traceback(
        &self,
        spec1: &MSSpectrum,
        spec2: &MSSpectrum,
    ) -> Result<Vec<(usize, usize)>> {
        validate_spectrum(spec1, true)?;
        validate_spectrum(spec2, true)?;
        let epsilon = float_parameter(&self.handler, "epsilon")?;
        // The source's sigma here has no zero-variance guard, unlike operator().
        let sigma = peak_distance_sigma(spec1, spec2, false)?;
        if sigma <= 0.0 {
            return Err(bad(
                "peak-distance variance is zero; the source divides by that sigma",
            ));
        }
        let matrix = AlignmentMatrix::fill(spec1, spec2, epsilon, sigma)?;
        Ok(matrix.traceback())
    }

    /// The source's heuristic shortcut: do the `level` most intense peaks of the
    /// two spectra share any m/z within `epsilon`?
    ///
    /// The source sorts copies of both spectra by intensity with `std::sort`,
    /// whose order among equal intensities is unspecified, so which peaks land
    /// in a tied top-`level` set is not defined there. This port sorts stably,
    /// which makes the selection deterministic without changing it whenever the
    /// intensities at the cut are distinct.
    fn heuristic_match(
        &self,
        spec1: &MSSpectrum,
        spec2: &MSSpectrum,
        epsilon: f64,
        level: usize,
    ) -> Result<bool> {
        let mut left = spec1.clone();
        let mut right = spec2.clone();
        left.sort_by_intensity(true)?;
        right.sort_by_intensity(true)?;
        for strong in &left.peaks[..level.min(left.len())] {
            for other in &right.peaks[..level.min(right.len())] {
                if (other.mz - strong.mz).abs() < epsilon {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }
}

impl PeakSpectrumCompareFunctor for PeakAlignment {
    fn handler(&self) -> &DefaultParamHandler {
        &self.handler
    }

    fn handler_mut(&mut self) -> &mut DefaultParamHandler {
        &mut self.handler
    }

    /// Normalised alignment score of `spec1` and `spec2`.
    ///
    /// The best cell of the last row or column, divided by the geometric mean of
    /// the two self-alignment scores. The best-cell search starts from
    /// `numeric_limits<double>::min()`, the smallest positive normal `f64` and
    /// **not** the most negative one, so a matrix whose last row and column are
    /// entirely negative reports that tiny positive number instead of its real
    /// maximum. This port reproduces that starting value, because the resulting
    /// score is what upstream callers have been comparing against.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] for unsorted peaks, which the source does
    /// not check, and [`Error::InvalidValue`] for a negative intensity under the
    /// score's `sqrt`, for an empty spectrum that reaches the alignment - where
    /// the source divides by a zero pair count and returns `inf` - for a
    /// non-finite parameter, for a negative `heuristic_level`, for more than
    /// [`MAX_ALIGNMENT_MATRIX_CELLS`] matrix cells, and for a zero or non-finite
    /// self-alignment product in the denominator.
    ///
    /// An empty spectrum returns `Ok(0.0)` whenever one of the two shortcuts
    /// fires first, which is what the class test observes: its empty spectrum
    /// carries no precursor, so the precursor distance to a real spectrum
    /// exceeds `precursor_mass_tolerance`.
    ///
    /// # The zero-variance guard cannot be reached usefully
    ///
    /// When every pairwise m/z distance is equal - a single peak on each side,
    /// say - the source substitutes `numeric_limits<double>::min()` for sigma.
    /// The position term is then `1 / (DBL_MIN * sqrt(2 pi))`, about `1.8e307`,
    /// so the **product** of the two self-alignment scores overflows to
    /// infinity for every nonzero `f32` intensity, down to the smallest
    /// subnormal, and the source's quotient silently becomes `0`: complete
    /// dissimilarity for a spectrum compared with itself. Here that overflow is
    /// [`Error::InvalidValue`], as is the zero denominator a spectrum of zero
    /// intensities produces.
    fn score(&self, spec1: &MSSpectrum, spec2: &MSSpectrum) -> Result<f64> {
        validate_spectrum(spec1, true)?;
        validate_spectrum(spec2, true)?;
        let epsilon = float_parameter(&self.handler, "epsilon")?;
        let precursor_tolerance = float_parameter(&self.handler, "precursor_mass_tolerance")?;
        let heuristic_level = count_parameter(&self.handler, "heuristic_level")? as usize;

        let left = spec1.precursors.first().map_or(0.0, |p| p.mz);
        let right = spec2.precursors.first().map_or(0.0, |p| p.mz);
        if (left - right).abs() > precursor_tolerance {
            return Ok(0.0);
        }
        if heuristic_level > 0 && !self.heuristic_match(spec1, spec2, epsilon, heuristic_level)? {
            return Ok(0.0);
        }

        let sigma = peak_distance_sigma(spec1, spec2, true)?;
        let matrix = AlignmentMatrix::fill(spec1, spec2, epsilon, sigma)?;
        let best = matrix.best_border_cell();
        let self1 = self_alignment_score(spec1, sigma)?;
        let self2 = self_alignment_score(spec2, sigma)?;
        let denominator = checked_sqrt(self1 * self2)?;
        if denominator == 0.0 {
            return Err(bad(
                "both self-alignment scores are zero; the source divides by that zero",
            ));
        }
        finite_score(best / denominator)
    }
}

/// `PeakAlignment::peakPairScore_`: the geometric mean of the two intensities
/// times a position term.
///
/// **The position term is not the Gaussian the formula looks like.** The source
/// writes `exp(-(fabs(pos1 - pos2)) / 2 * sigma * sigma)`, and C's precedence
/// reads that as `exp(((-|Δ|) / 2) * sigma * sigma)`: the distance enters
/// linearly and sigma **multiplies** the exponent instead of dividing it, so a
/// larger spread makes distant peaks score *less*, not more. The intended
/// `exp(-Δ² / (2σ²))` would need different parentheses. The expression is
/// transcribed exactly, because every published score from this class carries
/// it.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] for a negative intensity product, where the
/// source's `sqrt` yields NaN, and for a non-finite result.
fn peak_pair_score(pos1: f64, intens1: f64, pos2: f64, intens2: f64, sigma: f64) -> Result<f64> {
    let pi = checked_sqrt(intens1 * intens2)?;
    let pp = (1.0 / (sigma * (2.0 * crate::constants::PI).sqrt()))
        * (-((pos1 - pos2).abs()) / 2.0 * sigma * sigma).exp();
    finite_score(pi * pp)
}

/// Standard deviation of the pairwise m/z distance between two spectra, the
/// source's `mid`/`var`/`sigma` block.
///
/// Both accumulations run over every pair in row-major order, and both divide by
/// the pair count as an exact integer product widened to `f64`, exactly as the
/// source does.
///
/// `guard` selects `operator()`'s `(var == 0) ? numeric_limits<double>::min()`
/// fallback; `getAlignmentTraceback` computes the same variance without it.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when either spectrum is empty - the source
/// divides by a zero pair count and carries NaN into every cell - and when the
/// variance is not finite.
fn peak_distance_sigma(spec1: &MSSpectrum, spec2: &MSSpectrum, guard: bool) -> Result<f64> {
    let pairs = spec1
        .len()
        .checked_mul(spec2.len())
        .ok_or_else(|| bad("pairwise peak count overflows"))?;
    if pairs == 0 {
        return Err(bad(
            "peak alignment needs two non-empty spectra; the source divides by a zero pair count",
        ));
    }
    let mut mid = 0.0;
    for left in &spec1.peaks {
        for right in &spec2.peaks {
            mid += (left.mz - right.mz).abs();
        }
    }
    mid /= pairs as f64;
    let mut variance = 0.0;
    for left in &spec1.peaks {
        for right in &spec2.peaks {
            let deviation = (left.mz - right.mz).abs() - mid;
            variance += deviation * deviation;
        }
    }
    variance /= pairs as f64;
    finite_score(variance)?;
    if guard && variance == 0.0 {
        // The source's comment: "only in case of only two equal peaks in the
        // spectra sigma is 0".
        return Ok(f64::MIN_POSITIVE);
    }
    checked_sqrt(variance)
}

/// The source's `score_spec1`/`score_spec2`: every peak scored against itself.
///
/// # Errors
///
/// As [`peak_pair_score`].
fn self_alignment_score(spectrum: &MSSpectrum, sigma: f64) -> Result<f64> {
    let mut total = 0.0;
    for peak in &spectrum.peaks {
        let intensity = f64::from(peak.intensity);
        total += peak_pair_score(peak.mz, intensity, peak.mz, intensity, sigma)?;
    }
    finite_score(total)
}

/// `PeakAlignment`'s score matrix and the direction matrix its traceback reads.
///
/// One filling routine serves both entry points, so the score and the traceback
/// cannot disagree about the alignment; the source duplicates the loop and the
/// two copies already differ in their sigma.
struct AlignmentMatrix {
    rows: usize,
    cols: usize,
    cells: Vec<f64>,
    /// `1` from the diagonal, `0` from the left, `2` from above, with `0` also
    /// standing for "no strict winner", as the source's zero-filled matrix does.
    trace: Vec<u8>,
}

impl AlignmentMatrix {
    /// Fill the `(n + 1) * (m + 1)` score matrix and the `n * m` direction
    /// matrix.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a non-finite or negative `epsilon`,
    /// for more than [`MAX_ALIGNMENT_MATRIX_CELLS`] cells, for a failed
    /// allocation and for anything [`peak_pair_score`] refuses.
    fn fill(spec1: &MSSpectrum, spec2: &MSSpectrum, epsilon: f64, sigma: f64) -> Result<Self> {
        if !epsilon.is_finite() || epsilon < 0.0 {
            return Err(bad("epsilon must be finite and nonnegative"));
        }
        let rows = spec1
            .len()
            .checked_add(1)
            .ok_or_else(|| bad("alignment matrix size overflows"))?;
        let cols = spec2
            .len()
            .checked_add(1)
            .ok_or_else(|| bad("alignment matrix size overflows"))?;
        let count = rows
            .checked_mul(cols)
            .filter(|&n| n <= MAX_ALIGNMENT_MATRIX_CELLS)
            .ok_or_else(|| bad("alignment exceeds the alignment-matrix cell ceiling"))?;
        let mut cells = zeroed::<f64>(count)?;
        let mut trace = zeroed::<u8>(spec1.len() * spec2.len())?;
        // The gap cost is the same parameter as the match window.
        let gap = epsilon;
        for i in 1..rows {
            cells[i * cols] = -gap * i as f64;
        }
        for (j, cell) in cells.iter_mut().enumerate().take(cols).skip(1) {
            *cell = -gap * j as f64;
        }
        for i in 1..rows {
            for j in 1..cols {
                let left_peak = spec1.peaks[i - 1];
                let right_peak = spec2.peaks[j - 1];
                let from_left = cells[i * cols + j - 1] - gap;
                let from_above = cells[(i - 1) * cols + j] - gap;
                if (left_peak.mz - right_peak.mz).abs() <= epsilon {
                    let from_diagonal = cells[(i - 1) * cols + j - 1]
                        + peak_pair_score(
                            left_peak.mz,
                            f64::from(left_peak.intensity),
                            right_peak.mz,
                            f64::from(right_peak.intensity),
                            sigma,
                        )?;
                    cells[i * cols + j] = from_left.max(from_above.max(from_diagonal));
                    if from_diagonal > from_left && from_diagonal > from_above {
                        trace[(i - 1) * spec2.len() + j - 1] = 1;
                    } else if from_left > from_diagonal && from_left > from_above {
                        trace[(i - 1) * spec2.len() + j - 1] = 0;
                    } else if from_above > from_diagonal && from_above > from_left {
                        trace[(i - 1) * spec2.len() + j - 1] = 2;
                    }
                    // No strict winner: the cell keeps the zero it was filled
                    // with, which the traceback reads as "from the left".
                } else {
                    cells[i * cols + j] = from_left.max(from_above);
                    trace[(i - 1) * spec2.len() + j - 1] = u8::from(from_left <= from_above) * 2;
                }
            }
        }
        Ok(Self {
            rows,
            cols,
            cells,
            trace,
        })
    }

    /// The source's best-overall-score scan: the largest value in the last row
    /// or the last column, starting from `numeric_limits<double>::min()`.
    fn best_border_cell(&self) -> f64 {
        let mut best = f64::MIN_POSITIVE;
        for j in 0..self.cols {
            best = best.max(self.cells[(self.rows - 1) * self.cols + j]);
        }
        for i in 0..self.rows {
            best = best.max(self.cells[i * self.cols + self.cols - 1]);
        }
        best
    }

    /// Walk back from the best border cell, collecting the diagonal steps.
    fn traceback(&self) -> Vec<(usize, usize)> {
        let columns = self.cols - 1;
        let mut best = f64::MIN_POSITIVE;
        let mut row = 0;
        let mut column = 0;
        // Strict improvement only, so the first maximum of the last row wins and
        // the last-column scan overrides it only on a strictly larger value.
        for j in 0..self.cols {
            let value = self.cells[(self.rows - 1) * self.cols + j];
            if best < value {
                best = value;
                row = self.rows - 1;
                column = j;
            }
        }
        for i in 0..self.rows {
            let value = self.cells[i * self.cols + self.cols - 1];
            if best < value {
                best = value;
                row = i;
                column = self.cols - 1;
            }
        }
        let mut aligned = Vec::new();
        while row > 0 && column > 0 {
            match self.trace[(row - 1) * columns + column - 1] {
                1 => {
                    aligned.push((row - 1, column - 1));
                    row -= 1;
                    column -= 1;
                }
                0 => column -= 1,
                _ => row -= 1,
            }
        }
        aligned.reverse();
        aligned
    }
}

/// Peak filtering applied before a SpectraST comparison, with the source's
/// default arguments.
///
/// The source spells these as four defaulted parameters of
/// `SpectraSTSimilarityScore::preprocess`; they are gathered here because Rust
/// has no default arguments and because they always travel together.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectraStPreprocessing {
    /// Peaks at or below this intensity are dropped. Source default `2.01`,
    /// stored as `float` upstream and compared against the `float` intensity.
    pub remove_peak_intensity_threshold: f32,
    /// Peaks below `1 / cut_peaks_below` of the base peak are dropped. Source
    /// default `1000`.
    pub cut_peaks_below: u32,
    /// Fewer surviving peaks than this rejects the spectrum. Source default `5`.
    pub min_peak_number: usize,
    /// How many peaks are **examined**, not how many are kept. Source default
    /// `150`.
    pub max_peak_number: usize,
}

impl Default for SpectraStPreprocessing {
    fn default() -> Self {
        Self {
            remove_peak_intensity_threshold: 2.01,
            cut_peaks_below: 1000,
            min_peak_number: 5,
            max_peak_number: 150,
        }
    }
}

/// Dot product of SpectraST, with its dot-bias and delta-D companions.
///
/// Unlike the other peak-spectrum functors this score is meant for matching one
/// spectrum against a whole library: preprocess and transform every spectrum,
/// take the dot products, keep the best two, derive
/// [`delta_d`](Self::delta_d) from them and combine everything with
/// [`compute_f`](Self::compute_f). The method is H. Lam et al., "Development and
/// validation of a spectral library searching method for peptide identification
/// from MS/MS", Proteomics 7, 655-667, 2007.
///
/// # The scaling exponents
///
/// SpectraST scales intensity by `0.5` and m/z by `0`, and this implementation
/// applies exactly that: [`preprocess`](Self::preprocess) replaces every
/// surviving intensity with its square root and leaves m/z untouched, so no mass
/// weighting enters the dot product at all. The often-quoted `m/z^0.5` variant
/// of the published score is **not** implemented upstream, and inventing it here
/// would change every score. The only other transform is the normalisation in
/// [`transform`](Self::transform), which divides the binned vector by its own
/// Euclidean norm so that a spectrum scores exactly `1` against itself.
///
/// Binning is fixed: bin width `1`, absolute units, spread `1` and the low
/// resolution offset `0.4`, spelled out by [`Self::bin_config`]. The source's
/// own `// TODO: resolution seems rather low` sits on that line.
///
/// This is the one derivative of [`PeakSpectrumCompareFunctor`] that
/// `PeakSpectrumCompareFunctor.cpp` does not include for factory registration,
/// and the one whose constructor calls `setName` without `defaultsToParam_()`.
/// It registers no parameters, so that omission has no observable effect and is
/// reproduced as written.
///
/// See `docs/SPECTRAST_SIMILARITY_SCORE_SUPPORT.md`.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectraSTSimilarityScore {
    handler: DefaultParamHandler,
}

impl SpectraSTSimilarityScore {
    /// Construct the functor with the source's name and no parameters.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] only if the parameter handler rejects the
    /// fixed name.
    pub fn new() -> Result<Self> {
        let mut handler = DefaultParamHandler::new("PeakSpectrumCompareFunctor")?;
        // The source stops here: it never calls defaultsToParam_(), which is
        // harmless only because it registers no defaults either.
        handler.set_name("SpectraSTSimilarityScore")?;
        Ok(Self { handler })
    }

    /// The binning this score hard-codes:
    /// `BinnedSpectrum(spec, 1, false, 1, BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES)`.
    pub fn bin_config() -> BinConfig {
        BinConfig {
            size: 1.0,
            unit: BinUnit::Absolute,
            spread: 1,
            offset: 0.4,
            ..BinConfig::default()
        }
    }

    /// Dot product of two already binned spectra, the source's
    /// `operator()(const BinnedSpectrum&, const BinnedSpectrum&)`.
    ///
    /// No normalisation happens here; pass [`transform`](Self::transform)ed
    /// spectra to get a score in `[0, 1]`. The reduction is Eigen's sparse dot:
    /// `f32` products accumulated in `f32` over the shared bin indices, widened
    /// only on return.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the two spectra do not share a
    /// binning - the source hands mismatched vectors to Eigen, which asserts in
    /// a debug build and reads past the shorter one otherwise - when they
    /// together store more than [`MAX_COMPARED_BINS`] bins, and when the `f32`
    /// accumulation overflows.
    pub fn dot(&self, bin1: &BinnedSpectrum, bin2: &BinnedSpectrum) -> Result<f64> {
        compatible(bin1, bin2)?;
        preflight_bins(bin1, bin2)?;
        sparse_dot(bin1, bin2)
    }

    /// Bin `spectrum` and divide the result by its own Euclidean norm.
    ///
    /// The norm is Eigen's `SparseMatrixBase::norm()`: the `f32` square root of
    /// the `f32` sum of the squared stored coefficients. Stored zeros stay
    /// stored, because Eigen divides coefficients in place without pruning.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] for unsorted peaks and
    /// [`Error::InvalidValue`] for an m/z below zero, for a binning that
    /// overflows the configured limits, and for a zero or non-finite norm - a
    /// spectrum with no peaks or with only zero intensities - where the source
    /// divides by that zero and fills the vector with NaN.
    pub fn transform(&self, spectrum: &MSSpectrum) -> Result<BinnedSpectrum> {
        let mut binned = BinnedSpectrum::new(spectrum, Self::bin_config())?;
        let norm = sparse_norm(&binned)?;
        if !norm.is_finite() || norm <= 0.0 {
            return Err(bad(
                "cannot normalise a binned spectrum whose norm is zero or not finite",
            ));
        }
        for value in binned.bins.values_mut() {
            *value /= norm;
            if !value.is_finite() {
                return Err(bad("normalised bin is not finite"));
            }
        }
        Ok(binned)
    }

    /// How much of the dot product a few bins dominate.
    ///
    /// The numerator is the Euclidean norm of the element-wise product of the
    /// two binned vectors, again reduced in `f32`; the denominator is the dot
    /// product. `dot_product` is the source's `double dot_product = -1`
    /// sentinel: `None`, or any value that is not strictly positive, recomputes
    /// it from the two spectra. A denominator that is still not positive yields
    /// `0.0`, as the source's own guard does.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for incompatible binning, for more than
    /// [`MAX_COMPARED_BINS`] combined stored bins, for an `f32` overflow in the
    /// numerator, for a non-finite supplied `dot_product` - where the source's
    /// `<= 0.0` guard lets NaN through - and for a non-finite quotient.
    pub fn dot_bias(
        &self,
        bin1: &BinnedSpectrum,
        bin2: &BinnedSpectrum,
        dot_product: Option<f64>,
    ) -> Result<f64> {
        compatible(bin1, bin2)?;
        preflight_bins(bin1, bin2)?;
        if let Some(value) = dot_product {
            if !value.is_finite() {
                return Err(bad("a supplied dot product must be finite"));
            }
        }
        let mut squares = 0.0_f32;
        for (index, &left) in &bin1.bins {
            if let Some(&right) = bin2.bins.get(index) {
                let product = left * right;
                squares += product * product;
            }
        }
        if !squares.is_finite() {
            return Err(bad("dot-bias numerator overflows f32"));
        }
        let numerator = f64::from(squares.sqrt());
        let denominator = match dot_product {
            Some(value) if value > 0.0 => value,
            _ => self.dot(bin1, bin2)?,
        };
        if denominator <= 0.0 {
            return Ok(0.0);
        }
        finite_score(numerator / denominator)
    }

    /// Normalised distance between the best and the second best match,
    /// `(top_hit - runner_up) / top_hit`.
    ///
    /// The source notes that dot products range over `[0, 1]`; nothing checks
    /// that, and neither does this.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `top_hit` is zero, which is the
    /// source's `Exception::DivisionByZero`, and when either argument or the
    /// quotient is not finite.
    pub fn delta_d(&self, top_hit: f64, runner_up: f64) -> Result<f64> {
        if !top_hit.is_finite() || !runner_up.is_finite() {
            return Err(bad("delta_D needs two finite scores"));
        }
        if top_hit == 0.0 {
            return Err(bad("delta_D divides by a zero top hit"));
        }
        finite_score((top_hit - runner_up) / top_hit)
    }

    /// The overall SpectraST score,
    /// `0.6 * dot_product + 0.4 * delta_D - b`.
    ///
    /// The bias penalty `b` is a step function of `dot_bias`: `0.12` below
    /// `0.1` and on `(0.35, 0.4]`, `0.18` on `(0.4, 0.45]`, `0.24` above `0.45`,
    /// and zero on `[0.1, 0.35]`. The low-bias and high-bias penalties share a
    /// value; that is what `SpectraSTSimilarityScore.cpp:132` writes.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when any argument or the result is not
    /// finite. The source performs no such check and lets NaN select `b = 0`,
    /// because every comparison against NaN is false.
    pub fn compute_f(&self, dot_product: f64, delta_d: f64, dot_bias: f64) -> Result<f64> {
        if !dot_product.is_finite() || !delta_d.is_finite() || !dot_bias.is_finite() {
            return Err(bad("the SpectraST score needs three finite terms"));
        }
        let b = if dot_bias < 0.1 || (0.35 < dot_bias && dot_bias <= 0.4) {
            0.12
        } else if 0.4 < dot_bias && dot_bias <= 0.45 {
            0.18
        } else if dot_bias > 0.45 {
            0.24
        } else {
            0.0
        };
        finite_score(0.6 * dot_product + 0.4 * delta_d - b)
    }

    /// Filter `spectrum` in place and report whether it survives.
    ///
    /// Peaks are dropped unless their intensity exceeds both
    /// [`remove_peak_intensity_threshold`](SpectraStPreprocessing::remove_peak_intensity_threshold)
    /// and `1 / cut_peaks_below` of the base peak's intensity; every survivor
    /// keeps its m/z and takes the square root of its intensity, the SpectraST
    /// intensity exponent of `0.5`, computed in `f32` as upstream.
    ///
    /// # Two behaviours worth knowing before calling this
    ///
    /// * The header says the filter "cuts peaks exceeding the max_peak_number
    ///   most intense peaks". It does not. The spectrum is sorted by **m/z** and
    ///   the loop stops after examining
    ///   [`max_peak_number`](SpectraStPreprocessing::max_peak_number) peaks, so
    ///   what is kept is a prefix in m/z, not the strongest peaks, and the
    ///   result can be shorter than that bound.
    /// * The source assigns a fresh, default-constructed spectrum over the
    ///   argument, so **every piece of metadata is lost** - precursors,
    ///   retention time, native id, data arrays. That is reproduced, because a
    ///   SpectraST workflow's downstream numbers depend on the peaks it leaves
    ///   behind and silently keeping more state would be a different function.
    ///   Clone the spectrum first if the metadata matters.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an invalid spectrum, a negative m/z
    /// or intensity - the source's `sqrt` would yield NaN - a zero
    /// [`cut_peaks_below`](SpectraStPreprocessing::cut_peaks_below), where the
    /// source divides by zero and discards every peak, and a non-finite
    /// threshold. The spectrum is left untouched when any of these fires.
    pub fn preprocess(
        &self,
        spectrum: &mut MSSpectrum,
        options: SpectraStPreprocessing,
    ) -> Result<bool> {
        spectrum.validate()?;
        if options.cut_peaks_below == 0 {
            return Err(bad("cut_peaks_below must be positive"));
        }
        if !options.remove_peak_intensity_threshold.is_finite() {
            return Err(bad("the intensity threshold must be finite"));
        }
        if spectrum
            .peaks
            .iter()
            .any(|p| p.mz < 0.0 || p.intensity < 0.0)
        {
            return Err(bad(
                "SpectraST preprocessing requires nonnegative m/z and intensities",
            ));
        }
        let mut min_high_intensity = 0.0;
        if let Some(base) = spectrum.base_peak() {
            min_high_intensity =
                (1.0 / f64::from(options.cut_peaks_below)) * f64::from(base.intensity);
        }
        let mut kept = Vec::new();
        kept.try_reserve(options.max_peak_number.min(spectrum.len()))
            .map_err(|_| bad("preprocessed peak allocation failed"))?;
        let mut sorted = spectrum.clone();
        sorted.sort_by_position()?;
        for peak in sorted.peaks.iter().take(options.max_peak_number) {
            if peak.intensity > options.remove_peak_intensity_threshold
                && f64::from(peak.intensity) > min_high_intensity
            {
                kept.push(Peak1D::new(peak.mz, peak.intensity.sqrt()));
            }
        }
        let passed = kept.len() >= options.min_peak_number;
        *spectrum = MSSpectrum::from_peaks(kept);
        Ok(passed)
    }
}

impl PeakSpectrumCompareFunctor for SpectraSTSimilarityScore {
    fn handler(&self) -> &DefaultParamHandler {
        &self.handler
    }

    fn handler_mut(&mut self) -> &mut DefaultParamHandler {
        &mut self.handler
    }

    /// Normalised dot product of two peak spectra.
    ///
    /// Exactly `dot(transform(a), transform(b))`; the source spells the binning
    /// and the normalisation out a second time in this overload, with the same
    /// arguments, so the two cannot differ.
    ///
    /// # Errors
    ///
    /// As [`transform`](SpectraSTSimilarityScore::transform) and
    /// [`dot`](SpectraSTSimilarityScore::dot).
    fn score(&self, a: &MSSpectrum, b: &MSSpectrum) -> Result<f64> {
        self.dot(&self.transform(a)?, &self.transform(b)?)
    }
}

/// `Eigen::SparseMatrixBase::norm()` of a `SparseVector<float>`:
/// `sqrt(cwiseAbs2().sum())`, with the squares and their accumulation both in
/// `f32` and only the square root applied afterwards.
///
/// The accumulation order is ascending bin index, for the reason
/// [`sparse_sum`] documents: Eigen's dense redux over the stored value array is
/// vectorised, so this reproduces the reduction's `f32` precision rather than
/// its bit pattern on every build.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the `f32` accumulation overflows.
fn sparse_norm(spectrum: &BinnedSpectrum) -> Result<f32> {
    let mut total = 0.0_f32;
    for &value in spectrum.bins.values() {
        total += value * value;
    }
    if total.is_finite() {
        Ok(total.sqrt())
    } else {
        Err(bad("binned squared norm overflows f32"))
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
