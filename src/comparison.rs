// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Spectrum alignment and sparse binned comparisons ported from OpenMS.
//! See `docs/COMPARISON_SUPPORT.md` for source provenance, precision and limits.

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
    pub fn config(&self) -> BinConfig {
        self.config
    }
    pub fn bins(&self) -> &BTreeMap<usize, f32> {
        &self.bins
    }
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
