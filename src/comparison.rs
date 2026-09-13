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
//! Four of the seven `PeakSpectrumCompareFunctor` derivatives and the
//! alignment primitive they share are ported as the parameterised functors
//! [`SpectrumAligner`] (`COMPARISON/SpectrumAlignment.h`),
//! [`SpectrumAlignmentScorer`] (`COMPARISON/SpectrumAlignmentScore.h`),
//! [`ZhangSimilarityScorer`] (`COMPARISON/ZhangSimilarityScore.h`) and
//! [`SteinScottImproveScorer`] (`COMPARISON/SteinScottImproveScore.h`). Their
//! support documents are `docs/SPECTRUM_ALIGNMENT_SUPPORT.md`,
//! `docs/SPECTRUM_ALIGNMENT_SCORE_SUPPORT.md`,
//! `docs/ZHANG_SIMILARITY_SCORE_SUPPORT.md` and
//! `docs/STEIN_SCOTT_IMPROVE_SCORE_SUPPORT.md`.
//!
//! Those four reproduce the source's arithmetic exactly, including the `float`
//! intensity products the C++ computes before widening to `double`, and read
//! every setting from a [`DefaultParamHandler`] on each call, as the source
//! `operator()` reads `param_`. The earlier `Copy` configuration structs
//! [`SpectrumAlignmentScore`], [`ZhangSimilarityScore`] and
//! [`SteinScottImproveScore`] remain as `f64` conveniences with no parameter
//! surface; each functor's support document states the divergence in full.

use crate::param::{DefaultParamHandler, Param, ParamValue};
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
/// The base has **seven** derivatives in the pinned tree, not the six that
/// `PeakSpectrumCompareFunctor.cpp:11-16` `#include`s for factory registration:
/// `SpectrumCheapDPCorr`, `SpectrumPrecursorComparator`, `ZhangSimilarityScore`,
/// `SpectrumAlignmentScore`, `SteinScottImproveScore`, `PeakAlignment` and -
/// absent from that include list - `SpectraSTSimilarityScore`. All seven
/// override the one-spectrum overload as `operator()(spec, spec)`, so the
/// default [`self_score`](Self::self_score) is that delegation and an
/// implementor only overrides it to record a cheaper closed form.
///
/// **This crate ships no implementor of this trait.** Four of the seven source
/// derivatives exist in this module - [`SpectrumAlignmentScore`],
/// [`ZhangSimilarityScore`], [`SteinScottImproveScore`] and
/// [`SpectrumPrecursorComparator`] - but they were ported in an earlier wave as
/// typed `Copy` configuration structs with no [`DefaultParamHandler`], and
/// giving them one means porting the parameter tree each of them registers
/// upstream, which is their own headers' work; a handler that did not carry
/// those parameters would make `set_parameters` silently ineffective. The
/// remaining three - `SpectrumCheapDPCorr`, `PeakAlignment` and
/// `SpectraSTSimilarityScore` - are not ported at all. The trait is therefore
/// the *shape* of the base, usable by a caller's own functor and exercised in
/// `tests/comparison_functors.rs`, and the port of this header is `partial`
/// until its derivatives arrive; see
/// `docs/PEAK_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md`. Its sibling
/// [`BinnedSpectrumCompareFunctor`] has all three of its source derivatives
/// shipped here.
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
// The shared alignment primitive and four of the seven
// `PeakSpectrumCompareFunctor` derivatives:
//   COMPARISON/SpectrumAlignment.h
//   COMPARISON/SpectrumAlignmentScore.h
//   COMPARISON/ZhangSimilarityScore.h
//   COMPARISON/SteinScottImproveScore.h
// ---------------------------------------------------------------------------

/// Default ceiling on the dynamic-programming cells one alignment may fill.
///
/// The source's banded alignment allocates a `std::map` row per reference peak
/// and has no ceiling at all, so two large spectra are an unbounded allocation.
/// [`SpectrumAligner::max_cells`] and [`SpectrumAlignmentScorer::max_cells`]
/// start here and are checked before the matrix is built.
pub const DEFAULT_ALIGNMENT_CELLS: usize = 5_000_000;

/// Default ceiling on the candidate peak pairs one many-to-many score may
/// examine.
///
/// [`ZhangSimilarityScore`] and [`SteinScottImproveScore`] walk a sliding window
/// whose worst case is `|s1| * |s2|` comparisons; the source has no ceiling.
/// [`ZhangSimilarityScorer::max_pairs`] and
/// [`SteinScottImproveScorer::max_pairs`] start here. Candidates are counted as
/// they are examined, including those the match predicate then rejects, because
/// examining them is the cost being bounded.
pub const DEFAULT_SCORED_PAIRS: usize = 5_000_000;

/// Read a `double` parameter, as the source's `(double)param_.getValue(key)`.
fn parameter_float(handler: &DefaultParamHandler, key: &str) -> Result<f64> {
    handler.parameters().value(key)?.to_f64()
}

/// Read a flag, as the source's `param_.getValue(key).toBool()`.
fn parameter_bool(handler: &DefaultParamHandler, key: &str) -> Result<bool> {
    handler.parameters().value(key)?.to_bool()
}

/// Register a `"true"`/`"false"` flag with the source's valid-string restriction.
fn flag_default(defaults: &mut Param, key: &str, description: &str) -> Result<()> {
    defaults.set_value(key, ParamValue::String("false".into()), description, &[])?;
    defaults.set_valid_strings(key, &["true".to_string(), "false".to_string()])
}

/// The `float` product of two peak intensities, widened afterwards.
///
/// `Peak1D::getIntensity` returns `float`, so `s1[i].getIntensity() *
/// s2[j].getIntensity()` is a single-precision multiply in every one of these
/// scorers and only the *result* is promoted when it meets a `double`. Computing
/// the same product in `f64` would be exact - an `f32` product needs at most 48
/// mantissa bits - and would therefore *not* match the source. This rounding is
/// worth about 3.5e-9 relative on the upstream `DFPIANGER` fixture.
fn source_intensity_product(a: &crate::Peak1D, b: &crate::Peak1D) -> f64 {
    f64::from(a.intensity * b.intensity)
}

/// Walk the candidate peak pairs exactly as `ZhangSimilarityScore::operator()`
/// and `SteinScottImproveScore::operator()` do.
///
/// Both bodies are the same loop over a persistent `j_left` cursor that is only
/// advanced when a target peak lies at least the window below the current
/// reference peak, and they differ solely in `matches`: `fabs(d) < tolerance`
/// for Zhang, `fabs(d) <= 2 * epsilon` for Stein/Scott. `j_left` is read when a
/// reference peak's inner loop starts and written during it, so a write takes
/// effect on the *next* reference peak, never the current one - that is the
/// source's `for (Size j = j_left; ...)` initialisation, reproduced here.
///
/// Nothing is allocated, so exceeding `max_pairs` leaves both inputs and every
/// accumulator untouched by construction.
fn source_pair_walk(
    spec1: &MSSpectrum,
    spec2: &MSSpectrum,
    max_pairs: usize,
    matches: impl Fn(f64) -> bool,
    mut visit: impl FnMut(usize, usize) -> Result<()>,
) -> Result<()> {
    if max_pairs == 0 {
        return Err(bad("pair limit must be positive"));
    }
    let mut j_left = 0;
    let mut used = 0;
    for (i, peak1) in spec1.peaks.iter().enumerate() {
        let mut j = j_left;
        while j < spec2.len() {
            let pos2 = spec2.peaks[j].mz;
            if used == max_pairs {
                return Err(bad("comparison exceeds configured pair limit"));
            }
            used += 1;
            let distance = (peak1.mz - pos2).abs();
            if matches(distance) {
                visit(i, j)?;
            } else if pos2 > peak1.mz {
                break;
            } else {
                j_left = j;
            }
            j += 1;
        }
    }
    Ok(())
}

/// Sum of squared intensities, the source's `sum1`/`sum2` accumulation.
///
/// `pow(p.getIntensity(), 2)` resolves to the `double` overload, so each `float`
/// intensity is widened first and the square is exact; the running total is
/// `double`. `powi(2)` is the same single multiplication.
fn squared_intensity_sum(spectrum: &MSSpectrum) -> f64 {
    spectrum
        .peaks
        .iter()
        .map(|p| f64::from(p.intensity).powi(2))
        .sum()
}

/// Aligns the peaks of two sorted spectra.
///
/// **Method 1**: a banded alignment - the band width comes from the `tolerance`
/// parameter - when an absolute tolerance is given. The scoring function is the
/// m/z distance between peaks; intensity plays no role.
///
/// **Method 2**: when a relative tolerance (ppm) is specified, a simple matching
/// of peaks is performed. Peaks from `s1` - usually the theoretical spectrum -
/// are assigned to the closest peak in `s2` if it lies inside the tolerance
/// window.
///
/// A peak in `s2` can be matched to none, one or several peaks in `s1`; a peak
/// in `s1` is matched to none or one peak in `s2`. Intensity is ignored. The
/// source carries a `TODO` about the `O(|s1| * log(|s2|))` complexity of this
/// second method; the port's ppm path is a single forward merge, so it is
/// `O(|s1| + |s2|)`, and the `TODO` is discharged rather than carried over.
///
/// This is the source `SpectrumAlignment` class: a `DefaultParamHandler` whose
/// only member is the parameter tree, wrapped here by composition instead of
/// inheritance. The alignment itself is [`SpectrumAlignment`], which this type
/// configures from its parameters; the two cannot disagree because there is one
/// implementation. See `docs/SPECTRUM_ALIGNMENT_SUPPORT.md`.
///
/// ```
/// use openms::comparison::SpectrumAligner;
/// use openms::{MSSpectrum, Peak1D};
///
/// let reference = MSSpectrum::from_peaks(vec![Peak1D::new(100.0, 1.0), Peak1D::new(200.0, 1.0)]);
/// let target = MSSpectrum::from_peaks(vec![Peak1D::new(100.1, 1.0)]);
/// let aligner = SpectrumAligner::new()?;
/// assert_eq!(aligner.spectrum_alignment(&reference, &target)?, [(0, 0)]);
/// # Ok::<(), openms::Error>(())
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct SpectrumAligner {
    handler: DefaultParamHandler,
    /// Ceiling on dynamic-programming cells, defaulting to
    /// [`DEFAULT_ALIGNMENT_CELLS`]. Native: the source has no ceiling.
    pub max_cells: usize,
}

impl SpectrumAligner {
    /// Construct with the source's registered name and defaults: `tolerance`
    /// `0.3` and `is_relative_tolerance` `"false"`.
    ///
    /// Reproduces `SpectrumAlignment.cpp:17-22`, which names the handler
    /// `"SpectrumAlignment"`, registers both defaults with the `"true"`/`"false"`
    /// restriction on the flag, and finishes with `defaultsToParam_()`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] only if the parameter handler rejects the
    /// fixed name or defaults, which cannot happen for these literals.
    pub fn new() -> Result<Self> {
        let mut handler = DefaultParamHandler::new("SpectrumAlignment")?;
        let mut defaults = Param::new();
        defaults.set_value(
            "tolerance",
            ParamValue::Float(0.3),
            "Defines the absolute (in Da) or relative (in ppm) tolerance",
            &[],
        )?;
        flag_default(
            &mut defaults,
            "is_relative_tolerance",
            "If true, the 'tolerance' is interpreted as ppm-value",
        )?;
        handler.set_defaults(defaults)?;
        handler.defaults_to_parameters()?;
        Ok(Self {
            handler,
            max_cells: DEFAULT_ALIGNMENT_CELLS,
        })
    }

    /// The parameter surface the source inherits from `DefaultParamHandler`.
    /// Change settings through [`DefaultParamHandler::set_parameters`] on
    /// [`handler_mut`](Self::handler_mut), as the source's `setParameters` does.
    pub fn handler(&self) -> &DefaultParamHandler {
        &self.handler
    }

    /// Mutable parameter surface, for `setParameters` and `setName`.
    pub fn handler_mut(&mut self) -> &mut DefaultParamHandler {
        &mut self.handler
    }

    /// Registered name, as `DefaultParamHandler::getName`.
    pub fn name(&self) -> &str {
        self.handler.name()
    }

    /// Current `tolerance` and `is_relative_tolerance` as the matching window
    /// [`SpectrumAlignment`] applies.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when either parameter is missing or has
    /// the wrong type, which `set_parameters` prevents.
    pub fn tolerance(&self) -> Result<Tolerance> {
        let value = parameter_float(&self.handler, "tolerance")?;
        Ok(if parameter_bool(&self.handler, "is_relative_tolerance")? {
            Tolerance::Ppm(value)
        } else {
            Tolerance::Absolute(value)
        })
    }

    /// Ordered zero-based `(reference, target)` index pairs for `s1` and `s2`.
    ///
    /// The source signature is `void getSpectrumAlignment(vector<pair<Size,
    /// Size>>& alignment, const SpectrumType1& s1, const SpectrumType2& s2)`,
    /// which clears the out-parameter first; the port returns the vector.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] when either spectrum is not sorted by
    /// m/z - the source's `Exception::IllegalArgument`, "Input to
    /// SpectrumAlignment is not sorted!" - and [`Error::InvalidValue`] for a
    /// non-finite coordinate, a negative or non-finite tolerance, a ppm
    /// coordinate outside the `f32` range the source's `MatchedIterator`
    /// narrows to, or an alignment needing more than
    /// [`max_cells`](Self::max_cells) cells. Two empty spectra, or one empty
    /// spectrum, yield an empty alignment, as upstream: the matrix is
    /// initialised, neither loop body runs, and the traceback starts at `(0, 0)`.
    pub fn spectrum_alignment(
        &self,
        s1: &MSSpectrum,
        s2: &MSSpectrum,
    ) -> Result<Vec<(usize, usize)>> {
        SpectrumAlignment {
            tolerance: self.tolerance()?,
            max_cells: self.max_cells,
        }
        .align(s1, s2)
    }
}

/// Similarity score via spectra alignment.
///
/// This class implements a simple scoring based on the alignment of spectra.
/// The alignment is implemented in [`SpectrumAligner`] and performs a dynamic
/// programming alignment of the peaks, minimising the distances between the
/// aligned peaks and maximising the number of peak pairs.
///
/// The scoring is done via the simple formula `score = sum / sqrt(sum1 * sum2)`.
/// `sum` accumulates `sqrt(I1 * I2 * factor)` over the aligned peak pairs, and
/// `sum1` and `sum2` are the sums of the squared intensities of the two spectra.
/// The class comment's "with the given exponent (default is 2)" describes a
/// parameter this class has never registered; the exponent is fixed at two by
/// the `pow(getIntensity(), 2)` in the body.
///
/// A binned version of this scoring is implemented in the
/// [`BinnedSpectralContrastAngle`] family.
///
/// **This is not a cosine**: the numerator sums square roots of intensity
/// products while the denominator sums squares, so a self-score is generally
/// greater than one - `1.4845` on the upstream `DFPIANGER` fixture.
///
/// Parameters, all registered by `SpectrumAlignmentScore.cpp:18-26`:
///
/// | Key | Default | Meaning |
/// | --- | --- | --- |
/// | `tolerance` | `0.3` | absolute (Da) or relative (ppm) tolerance |
/// | `is_relative_tolerance` | `"false"` | interpret `tolerance` as ppm |
/// | `use_linear_factor` | `"false"` | weight intensities by the relative m/z difference |
/// | `use_gaussian_factor` | `"false"` | weight them by a Gaussian of that difference |
///
/// See `docs/SPECTRUM_ALIGNMENT_SCORE_SUPPORT.md`.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectrumAlignmentScorer {
    handler: DefaultParamHandler,
    /// Ceiling on dynamic-programming cells, defaulting to
    /// [`DEFAULT_ALIGNMENT_CELLS`]. Native: the source has no ceiling.
    pub max_cells: usize,
}

impl SpectrumAlignmentScorer {
    /// Construct with the source's registered name and its four defaults.
    ///
    /// Reproduces `SpectrumAlignmentScore.cpp:15-27`: the base constructor names
    /// the handler `"PeakSpectrumCompareFunctor"`, `setName` renames it, the four
    /// defaults are registered and `defaultsToParam_()` copies them into the
    /// current parameters.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] only if the parameter handler rejects the
    /// fixed name or defaults, which cannot happen for these literals.
    pub fn new() -> Result<Self> {
        let mut handler = DefaultParamHandler::new("PeakSpectrumCompareFunctor")?;
        handler.set_name("SpectrumAlignmentScore")?;
        let mut defaults = Param::new();
        defaults.set_value(
            "tolerance",
            ParamValue::Float(0.3),
            "Defines the absolute (in Da) or relative (in ppm) tolerance",
            &[],
        )?;
        flag_default(
            &mut defaults,
            "is_relative_tolerance",
            "if true, the tolerance value is interpreted as ppm",
        )?;
        flag_default(
            &mut defaults,
            "use_linear_factor",
            "if true, the intensities are weighted with the relative m/z difference",
        )?;
        flag_default(
            &mut defaults,
            "use_gaussian_factor",
            "if true, the intensities are weighted with the relative m/z difference using a gaussian",
        )?;
        handler.set_defaults(defaults)?;
        handler.defaults_to_parameters()?;
        Ok(Self {
            handler,
            max_cells: DEFAULT_ALIGNMENT_CELLS,
        })
    }
}

impl PeakSpectrumCompareFunctor for SpectrumAlignmentScorer {
    fn handler(&self) -> &DefaultParamHandler {
        &self.handler
    }

    fn handler_mut(&mut self) -> &mut DefaultParamHandler {
        &mut self.handler
    }

    /// `sum / sqrt(sum1 * sum2)` over the alignment of `a` and `b`.
    ///
    /// The denominator grouping is the source's: one `sqrt` of the product of
    /// the two squared-intensity sums, not the product of two `sqrt`s.
    ///
    /// With `is_relative_tolerance` the per-pair window is recomputed here as
    /// `tolerance * mz1 * 1e-6` in `double`, which is **not** the window that
    /// selected the pair: the alignment's `MatchedIterator` instantiates
    /// `Math::ppmToMass` at `float` and evaluates `(tolerance / 1e6) * mz` there.
    /// A pair can therefore be matched and then carry `mz_difference >
    /// mz_tolerance`, making a linear factor negative. The source takes the
    /// square root of that negative product and returns NaN; this port reports
    /// [`Error::InvalidValue`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] for unsorted peaks and
    /// [`Error::InvalidValue`] for a non-finite coordinate or intensity, a
    /// negative intensity, a negative or non-finite tolerance, an alignment
    /// exceeding [`max_cells`](Self::max_cells), both weighting flags set at
    /// once, a zero m/z tolerance under a weighting flag, a negative weighted
    /// product, or a non-finite score.
    ///
    /// Both weighting flags set is `OPENMS_PRECONDITION(!(use_linear_factor &&
    /// use_gaussian_factor), ...)` upstream, which is compiled out of a release
    /// build and then silently lets the linear factor win; here it is an error.
    ///
    /// A zero squared-intensity sum on either side yields `Ok(0.0)`. The source
    /// computes `0.0 / sqrt(0.0)` and returns NaN. That covers two empty
    /// spectra and one empty spectrum, and the module's binned scorers make the
    /// same choice. Two spectra with no peak inside the tolerance align to no
    /// pairs and score an unremarkable `0.0` in both the source and here.
    fn score(&self, a: &MSSpectrum, b: &MSSpectrum) -> Result<f64> {
        let tolerance = parameter_float(&self.handler, "tolerance")?;
        let relative = parameter_bool(&self.handler, "is_relative_tolerance")?;
        let linear = parameter_bool(&self.handler, "use_linear_factor")?;
        let gaussian = parameter_bool(&self.handler, "use_gaussian_factor")?;
        if linear && gaussian {
            return Err(bad(
                "use either 'use_linear_factor' or 'use_gaussian_factor', not both",
            ));
        }
        validate_spectrum(a, true)?;
        validate_spectrum(b, true)?;
        let aligner = SpectrumAlignment {
            tolerance: if relative {
                Tolerance::Ppm(tolerance)
            } else {
                Tolerance::Absolute(tolerance)
            },
            max_cells: self.max_cells,
        };
        let alignment = aligner.align(a, b)?;
        let sum1 = squared_intensity_sum(a);
        let sum2 = squared_intensity_sum(b);
        let mut sum = 0.0;
        for (i, j) in alignment {
            let (p, q) = (a.peaks[i], b.peaks[j]);
            // Source: `tolerance * s1[ap.first].getMZ() * 1e-6`, in this order.
            let mz_tolerance = if relative {
                tolerance * p.mz * 1e-6
            } else {
                tolerance
            };
            let mz_difference = (p.mz - q.mz).abs();
            let factor = if linear || gaussian {
                if mz_tolerance == 0.0 {
                    return Err(bad("distance weighting divides by a zero m/z tolerance"));
                }
                if linear {
                    (mz_tolerance - mz_difference) / mz_tolerance
                } else {
                    libm::erfc(mz_difference / (3.0 * mz_tolerance * std::f64::consts::SQRT_2))
                }
            } else {
                1.0
            };
            sum += checked_sqrt(source_intensity_product(&p, &q) * factor)?;
        }
        let denominator = checked_sqrt(sum1 * sum2)?;
        if denominator == 0.0 {
            return Ok(0.0);
        }
        finite_score(sum / denominator)
    }
}

/// Similarity score of Zhang.
///
/// The details of the score can be found in: Z. Zhang, Prediction of Low-Energy
/// Collision-Induced Dissociation Spectra of Peptides, Anal. Chem., 76 (14),
/// 3908-3922, 2004.
///
/// Every peak pair closer than `tolerance` contributes `sqrt(I1 * I2 * factor)`,
/// and the total is divided by `sqrt(sum1 * sum2)` where `sum1` and `sum2` are
/// the two spectra's total intensities. Unlike [`SpectrumAlignmentScorer`] this
/// is a many-to-many comparison: no alignment restricts a peak to one partner.
///
/// Parameters, all registered by `ZhangSimilarityScore.cpp:22-31`:
///
/// | Key | Default | Meaning |
/// | --- | --- | --- |
/// | `tolerance` | `0.2` | absolute (Da) or relative (ppm) tolerance |
/// | `is_relative_tolerance` | `"false"` | **unimplemented upstream**; see [`score`](PeakSpectrumCompareFunctor::score) |
/// | `use_linear_factor` | `"false"` | weight intensities by the relative m/z difference |
/// | `use_gaussian_factor` | `"false"` | weight them by a Gaussian of that difference |
///
/// See `docs/ZHANG_SIMILARITY_SCORE_SUPPORT.md`.
#[derive(Clone, Debug, PartialEq)]
pub struct ZhangSimilarityScorer {
    handler: DefaultParamHandler,
    /// Ceiling on examined candidate pairs, defaulting to
    /// [`DEFAULT_SCORED_PAIRS`]. Native: the source has no ceiling.
    pub max_pairs: usize,
}

impl ZhangSimilarityScorer {
    /// Construct with the source's registered name and its four defaults.
    ///
    /// Reproduces `ZhangSimilarityScore.cpp:20-32`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] only if the parameter handler rejects the
    /// fixed name or defaults, which cannot happen for these literals.
    pub fn new() -> Result<Self> {
        let mut handler = DefaultParamHandler::new("PeakSpectrumCompareFunctor")?;
        handler.set_name("ZhangSimilarityScore")?;
        let mut defaults = Param::new();
        defaults.set_value(
            "tolerance",
            ParamValue::Float(0.2),
            "defines the absolute (in Da) or relative (in ppm) tolerance",
            &[],
        )?;
        flag_default(
            &mut defaults,
            "is_relative_tolerance",
            "If set to true, the tolerance is interpreted as relative",
        )?;
        flag_default(
            &mut defaults,
            "use_linear_factor",
            "if true, the intensities are weighted with the relative m/z difference",
        )?;
        flag_default(
            &mut defaults,
            "use_gaussian_factor",
            "if true, the intensities are weighted with the relative m/z difference using a gaussian",
        )?;
        handler.set_defaults(defaults)?;
        handler.defaults_to_parameters()?;
        Ok(Self {
            handler,
            max_pairs: DEFAULT_SCORED_PAIRS,
        })
    }

    /// The source's protected `getFactor_(mz_tolerance, mz_difference,
    /// is_gaussian)`.
    ///
    /// Gaussian: `erfc(mz_difference / (mz_tolerance * 3 * sqrt(2)))`. Linear:
    /// `(mz_tolerance - mz_difference) / mz_tolerance`.
    ///
    /// **The source caches the Gaussian denominator in a function-local
    /// `static const double`**, so the very first call in a process fixes
    /// `mz_tolerance * 3 * sqrt(2)` for every later call, whatever tolerance or
    /// instance it belongs to. This port evaluates it per call. Reproducing the
    /// cache would mean carrying a process-global whose value depends on which
    /// score ran first, which is not a behaviour a caller can rely on; it is
    /// recorded as an upstream defect instead.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a non-positive or non-finite
    /// `mz_tolerance`, which divides by zero upstream. The scoring loop cannot
    /// reach it: a pair exists only when `mz_difference < mz_tolerance` and
    /// `mz_difference` is a non-negative absolute value.
    pub fn factor(mz_tolerance: f64, mz_difference: f64, is_gaussian: bool) -> Result<f64> {
        if !mz_tolerance.is_finite() || mz_tolerance <= 0.0 {
            return Err(bad("Zhang weighting requires a positive m/z tolerance"));
        }
        Ok(if is_gaussian {
            libm::erfc(mz_difference / (mz_tolerance * 3.0 * 2.0_f64.sqrt()))
        } else {
            (mz_tolerance - mz_difference) / mz_tolerance
        })
    }
}

impl PeakSpectrumCompareFunctor for ZhangSimilarityScorer {
    fn handler(&self) -> &DefaultParamHandler {
        &self.handler
    }

    fn handler_mut(&mut self) -> &mut DefaultParamHandler {
        &mut self.handler
    }

    /// `sum / sqrt(sum1 * sum2)` over every peak pair closer than `tolerance`.
    ///
    /// `sum1` and `sum2` are total intensities, not squared ones - that is the
    /// difference from [`SpectrumAlignmentScorer`], together with the
    /// many-to-many pairing.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Unsupported`] when `is_relative_tolerance` is set: the
    /// source throws `Exception::NotImplemented` and carries a `TODO` to remove
    /// the parameter. The parameter is still registered here so that a `Param`
    /// tree round-trips between the two.
    ///
    /// Returns [`Error::UnsortedData`] for unsorted peaks and
    /// [`Error::InvalidValue`] for a non-finite coordinate or intensity, a
    /// negative intensity, a non-finite tolerance, more than
    /// [`max_pairs`](Self::max_pairs) examined candidates, both weighting flags
    /// set at once, or a non-finite score. The source does not check
    /// sortedness, although its sliding `j_left` cursor requires it.
    ///
    /// A zero total intensity on either side yields `Ok(0.0)`, covering two
    /// empty spectra and one empty spectrum; the source computes `0.0 /
    /// sqrt(0.0)` and returns NaN. Two spectra with no peak inside the tolerance
    /// produce no pairs and score `0.0` in both.
    ///
    /// Setting both weighting flags is accepted upstream and the Gaussian wins,
    /// because `getFactor_` takes `use_gaussian_factor` as its switch. The
    /// sibling [`SpectrumAlignmentScorer`] resolves the same clash the other way
    /// and its debug-only precondition rejects it. The port refuses it in both.
    fn score(&self, a: &MSSpectrum, b: &MSSpectrum) -> Result<f64> {
        let tolerance = parameter_float(&self.handler, "tolerance")?;
        if parameter_bool(&self.handler, "is_relative_tolerance")? {
            return Err(Error::Unsupported(
                "ZhangSimilarityScore does not implement a relative tolerance".into(),
            ));
        }
        let linear = parameter_bool(&self.handler, "use_linear_factor")?;
        let gaussian = parameter_bool(&self.handler, "use_gaussian_factor")?;
        if linear && gaussian {
            return Err(bad(
                "use either 'use_linear_factor' or 'use_gaussian_factor', not both",
            ));
        }
        if !tolerance.is_finite() {
            return Err(bad("Zhang tolerance must be finite"));
        }
        validate_spectrum(a, true)?;
        validate_spectrum(b, true)?;
        let sum1: f64 = a.peaks.iter().map(|p| f64::from(p.intensity)).sum();
        let sum2: f64 = b.peaks.iter().map(|p| f64::from(p.intensity)).sum();
        let mut sum = 0.0;
        source_pair_walk(
            a,
            b,
            self.max_pairs,
            |distance| distance < tolerance,
            |i, j| {
                let (p, q) = (a.peaks[i], b.peaks[j]);
                let factor = if linear || gaussian {
                    Self::factor(tolerance, (p.mz - q.mz).abs(), gaussian)?
                } else {
                    1.0
                };
                sum += checked_sqrt(source_intensity_product(&p, &q) * factor)?;
                Ok(())
            },
        )?;
        let denominator = checked_sqrt(sum1 * sum2)?;
        if denominator == 0.0 {
            return Ok(0.0);
        }
        finite_score(sum / denominator)
    }
}

/// Similarity score based on Stein and Scott.
///
/// This is a pairwise score function. The spectrum contains peaks, and each peak
/// is defined by two values, m/z and intensity. The score function takes the sum
/// of the products of the peak intensities from spectrum 1 and spectrum 2, but
/// only where the m/z distance between the two peaks is smaller than a given
/// window size; by default the window is the accuracy of the mass spectrometer.
/// That sum is normalised by dividing it by a distance function,
/// `sqrt(sum of squared intensities of spectrum 1 * the same for spectrum 2)`.
///
/// To distinguish close from distant spectra an additional term is subtracted.
/// It denotes the expected value of both spectra under random placement of all
/// peaks within the given mass-to-charge range. The probability that two peaks
/// with randomised intensity values lie within two epsilon of each other is a
/// constant proportional to epsilon, so the additional term is that constant
/// times the product of the two spectra's total intensities.
///
/// The details of the score can be found in: Signal Maps for Mass
/// Spectrometry-based Comparative Proteomics; Amol Prakash, Parag Mallick,
/// Jeffrey Whiteaker, Heidi Zhang, Amanda Paulovich, Mark Flory, Hookeun Lee,
/// Ruedi Aebersold and Benno Schwikowski.
///
/// Note that the window actually applied is `2 * tolerance` and that its
/// boundary is inclusive, while the constant subtracted is `tolerance / 10000`.
///
/// Parameters, both registered by `SteinScottImproveScore.cpp:21-23`:
///
/// | Key | Default | Meaning |
/// | --- | --- | --- |
/// | `tolerance` | `0.2` | the absolute error of the mass spectrometer |
/// | `threshold` | `0.2` | a score below this is reported as zero |
///
/// Neither carries a valid-string or range restriction upstream, so neither does
/// here. See `docs/STEIN_SCOTT_IMPROVE_SCORE_SUPPORT.md`.
#[derive(Clone, Debug, PartialEq)]
pub struct SteinScottImproveScorer {
    handler: DefaultParamHandler,
    /// Ceiling on examined candidate pairs, defaulting to
    /// [`DEFAULT_SCORED_PAIRS`]. Native: the source has no ceiling.
    pub max_pairs: usize,
}

impl SteinScottImproveScorer {
    /// Construct with the source's registered name and its two defaults.
    ///
    /// Reproduces `SteinScottImproveScore.cpp:17-24`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] only if the parameter handler rejects the
    /// fixed name or defaults, which cannot happen for these literals.
    pub fn new() -> Result<Self> {
        let mut handler = DefaultParamHandler::new("PeakSpectrumCompareFunctor")?;
        handler.set_name("SteinScottImproveScore")?;
        let mut defaults = Param::new();
        defaults.set_value(
            "tolerance",
            ParamValue::Float(0.2),
            "defines the absolute error of the mass spectrometer",
            &[],
        )?;
        defaults.set_value(
            "threshold",
            ParamValue::Float(0.2),
            "if the calculated score is smaller than the threshold, a zero is given back",
            &[],
        )?;
        handler.set_defaults(defaults)?;
        handler.defaults_to_parameters()?;
        Ok(Self {
            handler,
            max_pairs: DEFAULT_SCORED_PAIRS,
        })
    }
}

impl PeakSpectrumCompareFunctor for SteinScottImproveScorer {
    fn handler(&self) -> &DefaultParamHandler {
        &self.handler
    }

    fn handler_mut(&mut self) -> &mut DefaultParamHandler {
        &mut self.handler
    }

    /// `(sum - z) / sqrt(sum1 * sum2)`, reported as zero below `threshold`.
    ///
    /// `sum` adds `I1 * I2` over every peak pair no further apart than
    /// `2 * tolerance`, `sum1` and `sum2` are squared-intensity sums, and
    /// `z = tolerance / 10000 * (total1 * total2)` with that exact grouping.
    ///
    /// `threshold` is compared after a narrowing to `f32`, because the source
    /// writes `score < (float)param_.getValue("threshold")` and the `float` is
    /// then widened again for the comparison. The default `0.2` therefore
    /// compares against `0.20000000298023224`, not against `0.2`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] for unsorted peaks and
    /// [`Error::InvalidValue`] for a non-finite coordinate or intensity, a
    /// negative intensity, a non-finite tolerance or threshold, more than
    /// [`max_pairs`](Self::max_pairs) examined candidates, or a non-finite
    /// score. The source does not check sortedness, although its sliding
    /// `j_left` cursor requires it.
    ///
    /// A zero squared-intensity sum on either side yields `Ok(0.0)`, covering
    /// two empty spectra and one empty spectrum; the source divides zero by
    /// zero, gets NaN, finds `NaN < threshold` false and returns the NaN. Two
    /// spectra with no peak inside the window give `sum == 0`, so the score is
    /// `-z / sqrt(sum1 * sum2)`, a negative number that the default threshold
    /// then reports as `0.0` - in the source and here alike.
    fn score(&self, a: &MSSpectrum, b: &MSSpectrum) -> Result<f64> {
        let epsilon = parameter_float(&self.handler, "tolerance")?;
        let threshold = self.handler.parameters().value("threshold")?.to_f32()?;
        if !epsilon.is_finite() {
            return Err(bad("Stein/Scott tolerance must be finite"));
        }
        if !threshold.is_finite() {
            return Err(bad("Stein/Scott threshold must be finite"));
        }
        validate_spectrum(a, true)?;
        validate_spectrum(b, true)?;
        let constant = epsilon / 10000.0;
        let sum1 = squared_intensity_sum(a);
        let sum2 = squared_intensity_sum(b);
        let sum3: f64 = a.peaks.iter().map(|p| f64::from(p.intensity)).sum();
        let sum4: f64 = b.peaks.iter().map(|p| f64::from(p.intensity)).sum();
        let z = constant * (sum3 * sum4);
        let mut sum = 0.0;
        source_pair_walk(
            a,
            b,
            self.max_pairs,
            |distance| distance <= 2.0 * epsilon,
            |i, j| {
                sum += source_intensity_product(&a.peaks[i], &b.peaks[j]);
                Ok(())
            },
        )?;
        let denominator = checked_sqrt(sum1 * sum2)?;
        if denominator == 0.0 {
            return Ok(0.0);
        }
        let score = finite_score((sum - z) / denominator)?;
        Ok(if score < f64::from(threshold) {
            0.0
        } else {
            score
        })
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
