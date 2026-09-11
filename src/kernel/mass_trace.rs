// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Owned trace peaks and source-compatible cached centroid/width calculations.

use super::{ConvexHull2D, Peak2D, Point2D};
use crate::{Error, Result};
use std::{
    mem::size_of,
    ops::{Index, IndexMut},
    slice,
};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MassTraceQuantMethod {
    #[default]
    Area,
    Median,
    MaxHeight,
}
impl MassTraceQuantMethod {
    pub const ALL: [Self; 3] = [Self::Area, Self::Median, Self::MaxHeight];
    pub const NAMES: [&'static str; 3] = ["area", "median", "max_height"];
    pub const fn name(self) -> &'static str {
        Self::NAMES[self as usize]
    }
    /// Exact source names; None replaces the source invalid enum sentinel.
    pub fn from_name(name: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|method| method.name() == name)
    }
}

/// One call's input count, weighted visits/comparisons and new logical payload.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MassTraceLimits {
    pub max_peaks: usize,
    pub max_work: usize,
    pub max_bytes: usize,
}
impl Default for MassTraceLimits {
    fn default() -> Self {
        Self {
            max_peaks: 1_000_000,
            max_work: 50_000_000,
            max_bytes: 256 * 1024 * 1024,
        }
    }
}

/// Source state stays cached until an explicit update. Editing a borrowed peak
/// never silently recomputes centroids, widths or previously supplied smoothing.
/// Ordinary Clone, indexing, equality and caller destruction have Rust costs.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MassTrace {
    peaks: Vec<Peak2D>,
    smoothed: Vec<f64>,
    label: String,
    centroid_mz: f64,
    centroid_rt: f64,
    centroid_sd: f64,
    centroid_im: f64,
    has_centroid_im: bool,
    fwhm: f64,
    fwhm_borders: (usize, usize),
    quant_method: MassTraceQuantMethod,
    pub fwhm_mz_avg: f64,
    pub fwhm_im_avg: f64,
    pub limits: MassTraceLimits,
}

impl MassTrace {
    pub fn new() -> Self {
        Self::default()
    }
    /// Takes a preassembled vector; peak fields are checked when consumed.
    pub fn from_peaks(peaks: Vec<Peak2D>) -> Result<Self> {
        Self::from_peaks_with_limits(peaks, MassTraceLimits::default())
    }
    pub fn from_peaks_with_limits(peaks: Vec<Peak2D>, limits: MassTraceLimits) -> Result<Self> {
        Work::new(limits, peaks.len())?;
        Ok(Self {
            peaks,
            limits,
            ..Self::default()
        })
    }
    pub fn from_slice(peaks: &[Peak2D]) -> Result<Self> {
        Self::from_slice_with_limits(peaks, MassTraceLimits::default())
    }
    pub fn from_slice_with_limits(peaks: &[Peak2D], limits: MassTraceLimits) -> Result<Self> {
        let mut work = Work::new(limits, peaks.len())?;
        let mut owned = work.vector(peaks.len())?;
        owned.extend_from_slice(peaks);
        Self::from_peaks_with_limits(owned, limits)
    }
    pub fn len(&self) -> usize {
        self.peaks.len()
    }
    pub fn is_empty(&self) -> bool {
        self.peaks.is_empty()
    }
    pub fn peaks(&self) -> &[Peak2D] {
        &self.peaks
    }
    pub fn peaks_mut(&mut self) -> &mut [Peak2D] {
        &mut self.peaks
    }
    pub fn iter(&self) -> slice::Iter<'_, Peak2D> {
        self.peaks.iter()
    }
    pub fn iter_mut(&mut self) -> slice::IterMut<'_, Peak2D> {
        self.peaks.iter_mut()
    }
    pub fn get(&self, index: usize) -> Option<&Peak2D> {
        self.peaks.get(index)
    }
    pub fn get_mut(&mut self, index: usize) -> Option<&mut Peak2D> {
        self.peaks.get_mut(index)
    }
    pub fn label(&self) -> &str {
        &self.label
    }
    pub fn set_label(&mut self, label: &str) -> Result<()> {
        let mut work = Work::new(self.limits, 0)?;
        work.allocate::<u8>(label.len())?;
        let mut owned = String::new();
        owned
            .try_reserve_exact(label.len())
            .map_err(|_| resource())?;
        owned.push_str(label);
        self.label = owned;
        Ok(())
    }
    pub const fn centroid_mz(&self) -> f64 {
        self.centroid_mz
    }
    pub const fn centroid_rt(&self) -> f64 {
        self.centroid_rt
    }
    pub const fn centroid_sd(&self) -> f64 {
        self.centroid_sd
    }
    pub const fn centroid_im(&self) -> f64 {
        self.centroid_im
    }
    pub const fn contains_im_data(&self) -> bool {
        self.has_centroid_im
    }
    pub fn set_centroid_sd(&mut self, value: f64) -> Result<()> {
        self.centroid_sd = finite(value)?;
        Ok(())
    }
    pub fn set_centroid_im(&mut self, value: f64) -> Result<()> {
        self.centroid_im = finite(value)?;
        self.has_centroid_im = true;
        Ok(())
    }
    pub const fn fwhm(&self) -> f64 {
        self.fwhm
    }
    pub const fn fwhm_borders(&self) -> (usize, usize) {
        self.fwhm_borders
    }
    pub const fn quant_method(&self) -> MassTraceQuantMethod {
        self.quant_method
    }
    pub fn set_quant_method(&mut self, method: MassTraceQuantMethod) {
        self.quant_method = method;
    }
    pub fn smoothed_intensities(&self) -> &[f64] {
        &self.smoothed
    }
    pub fn set_smoothed_intensities(&mut self, values: &[f64]) -> Result<()> {
        if values.len() != self.len() {
            return Err(bad("smoothed intensity count differs from trace size"));
        }
        let mut work = self.work()?;
        work.scan(values.len())?;
        for &value in values {
            finite(value)?;
        }
        let mut owned = work.vector(values.len())?;
        owned.extend_from_slice(values);
        self.smoothed = owned;
        Ok(())
    }
    pub fn trace_length(&self) -> Result<f64> {
        if self.len() <= 1 {
            return Ok(0.0);
        }
        self.work()?.scan(2)?;
        finite((finite(self.peaks[self.len() - 1].rt())? - finite(self.peaks[0].rt())?).abs())
    }
    pub fn average_ms1_cycle_time(&self) -> Result<f64> {
        if self.len() <= 1 {
            return Ok(0.0);
        }
        self.work()?.scan(2)?;
        finite(
            (finite(self.peaks[self.len() - 1].rt())? - finite(self.peaks[0].rt())?)
                / (self.len() - 1) as f64,
        )
    }
    pub fn compute_intensity_sum(&self) -> Result<f64> {
        self.work()?.scan(self.len())?;
        let mut sum = 0.0;
        for p in &self.peaks {
            sum = finite(sum + finite(f64::from(p.intensity))?)?;
        }
        Ok(sum)
    }
    pub fn compute_peak_area(&self) -> Result<f64> {
        let mut work = self.work()?;
        if self.is_empty() {
            return Ok(0.0);
        }
        work.scan(add(self.len(), 1)?)?;
        let mut previous_i = finite(f64::from(self.peaks[0].intensity))?;
        let mut previous_rt = finite(self.peaks[0].rt())?;
        let mut area = 0.0;
        // Include the source's first zero-width trapezoid and its operation order.
        for p in &self.peaks {
            let rt = finite(p.rt())?;
            let intensity = finite(f64::from(p.intensity))?;
            area = finite(area + (previous_i + intensity) / 2.0 * (rt - previous_rt))?;
            previous_i = intensity;
            previous_rt = rt;
        }
        Ok(area)
    }
    /// Preserves the source mixed raw/smoothed trapezoid expression exactly.
    pub fn compute_smoothed_peak_area(&self) -> Result<f64> {
        self.require_smoothed()?;
        self.work()?.scan(self.len())?;
        let mut previous_i = self.smoothed[0];
        let mut previous_rt = finite(self.peaks[0].rt())?;
        let mut area = 0.0;
        for i in 1..self.len() {
            let intensity = finite(f64::from(self.peaks[i].intensity))?;
            let rt = finite(self.peaks[i].rt())?;
            if self.smoothed[i] > 0.0 {
                area = finite(area + (previous_i + intensity) / 2.0 * (rt - previous_rt))?;
            }
            previous_i = intensity;
            previous_rt = rt;
        }
        Ok(area)
    }
    pub fn find_max_by_int_peak(&self, smoothed: bool) -> Result<usize> {
        self.find_max(smoothed, &mut self.work()?)
    }
    fn find_max(&self, smoothed: bool, work: &mut Work) -> Result<usize> {
        if smoothed {
            self.require_smoothed()?;
        }
        self.require_nonempty()?;
        work.scan(add(self.len(), 1)?)?;
        let mut best = 0;
        let mut highest = self.value(0, smoothed)?;
        for i in 0..self.len() {
            let value = self.value(i, smoothed)?;
            if value > highest {
                best = i;
                highest = value;
            }
        }
        Ok(best)
    }
    pub fn max_intensity(&self, smoothed: bool) -> Result<f64> {
        let n = if smoothed {
            self.smoothed.len()
        } else {
            self.len()
        };
        Work::new(self.limits, n)?.scan(n)?;
        let mut highest = 0.0;
        for i in 0..n {
            let value = self.value(i, smoothed)?;
            if value > highest {
                highest = value;
            }
        }
        Ok(highest)
    }
    pub fn estimate_fwhm(&mut self, smoothed: bool) -> Result<f64> {
        let mut work = self.work()?;
        work.scan(self.len())?;
        let mut previous = None;
        for p in &self.peaks {
            let rt = finite(p.rt())?;
            if previous.is_some_and(|before| before > rt) {
                return Err(bad("FWHM requires nondecreasing RT"));
            }
            previous = Some(rt);
        }
        let apex = self.find_max(smoothed, &mut work)?;
        if apex == 0 || apex == self.len() - 1 {
            // Source resets the borders but intentionally leaves cached fwhm_ alone.
            self.fwhm_borders = (0, 0);
            return Ok(0.0);
        }
        work.scan(add(self.len(), 8)?)?;
        let half = self.value(apex, smoothed)? / 2.0;
        let (mut left, mut right) = (apex, apex);
        while left > 0 && self.value(left, smoothed)? >= half {
            left -= 1;
        }
        while right + 1 < self.len() && self.value(right, smoothed)? >= half {
            right += 1;
        }
        let crossing = |bottom: usize, top: usize| -> Result<f64> {
            let low = self.value(bottom, smoothed)?;
            if low > half {
                Ok(self.peaks[bottom].rt())
            } else {
                interpolate(
                    self.peaks[bottom].rt(),
                    self.peaks[top].rt(),
                    low,
                    self.value(top, smoothed)?,
                    half,
                )
            }
        };
        let begin = crossing(left, left + 1)?;
        let end = crossing(right, right - 1)?;
        let width = finite((end - begin).abs())?;
        if width > (self.peaks[right].rt() - self.peaks[left].rt()).abs() {
            return Err(bad("interpolated FWHM exceeds its bracket"));
        }
        self.fwhm_borders = (left, right);
        self.fwhm = width;
        Ok(width)
    }
    pub fn compute_fwhm_area(&self) -> Result<f64> {
        self.fwhm_area(false, &mut self.work()?)
    }
    pub fn compute_fwhm_area_smooth(&self) -> Result<f64> {
        self.fwhm_area(true, &mut self.work()?)
    }
    fn fwhm_area(&self, smoothed: bool, work: &mut Work) -> Result<f64> {
        let (left, right) = self.fwhm_borders;
        if (left, right) == (0, 0) {
            return Ok(0.0);
        }
        if smoothed {
            self.require_smoothed()?;
        }
        work.scan(right - left + 1)?;
        let mut previous_i = self.value(left, smoothed)?;
        let mut previous_rt = finite(self.peaks[left].rt())?;
        let mut area = 0.0;
        for i in left + 1..=right {
            let intensity = self.value(i, smoothed)?;
            let rt = finite(self.peaks[i].rt())?;
            area = finite(area + (previous_i + intensity) / 2.0 * (rt - previous_rt))?;
            previous_i = intensity;
            previous_rt = rt;
        }
        Ok(area)
    }
    pub fn intensity(&self, smoothed: bool) -> Result<f64> {
        match self.quant_method {
            MassTraceQuantMethod::Area => self.fwhm_area(smoothed, &mut self.work()?),
            MassTraceQuantMethod::Median => self.median(|p| f64::from(p.intensity), true),
            MassTraceQuantMethod::MaxHeight => self.max_intensity(smoothed),
        }
    }
    pub fn convex_hull(&self) -> Result<ConvexHull2D> {
        let mut work = self.work()?;
        // Bound the existing geometry constructor: Point2D input, Scan output,
        // stable-sort scratch and conservative O(n log n) comparisons.
        work.sort(self.len())?;
        work.allocate::<[f64; 6]>(self.len())?;
        let mut points = work.vector(self.len())?;
        for p in &self.peaks {
            points.push(Point2D::new(finite(p.rt())?, finite(p.mz())?));
        }
        ConvexHull2D::from_points(&points)
    }
    pub fn update_weighted_mean_rt(&mut self) -> Result<()> {
        self.require_nonempty()?;
        self.work()?.scan(self.len())?;
        let result = if self.len() == 1 {
            finite(self.peaks[0].rt())?
        } else {
            let (mut area, mut weighted) = (0.0, 0.0);
            let mut previous_rt = finite(self.peaks[0].rt())?;
            for p in &self.peaks[1..] {
                let rt = finite(p.rt())?;
                let intensity = finite(f64::from(p.intensity))?;
                let difference = rt - previous_rt;
                weighted = finite(weighted + intensity * rt * difference)?;
                area = finite(area + intensity * difference)?;
                previous_rt = rt;
            }
            finite(weighted / area)?
        };
        self.centroid_rt = result;
        Ok(())
    }
    pub fn update_smoothed_weighted_mean_rt(&mut self) -> Result<()> {
        self.require_smoothed()?;
        self.work()?.scan(self.len())?;
        let result = if self.len() == 1 {
            finite(self.peaks[0].rt())?
        } else {
            let (mut area, mut weighted) = (0.0, 0.0);
            for (p, &intensity) in self.peaks.iter().zip(&self.smoothed) {
                if intensity > 0.0 {
                    weighted = finite(weighted + intensity * finite(p.rt())?)?;
                    area = finite(area + intensity)?;
                }
            }
            positive_weight(area)?;
            finite(weighted / area)?
        };
        self.centroid_rt = result;
        Ok(())
    }
    pub fn update_smoothed_max_rt(&mut self) -> Result<()> {
        self.require_smoothed()?;
        self.work()?.scan(self.len())?;
        let index = if self.len() == 1 {
            0
        } else {
            let (mut highest, mut index) = (-1.0, 0);
            for (i, &value) in self.smoothed.iter().enumerate() {
                if value > highest {
                    highest = value;
                    index = i;
                }
            }
            if highest <= 0.0 {
                return Err(bad("smoothed apex must be positive"));
            }
            index
        };
        self.centroid_rt = finite(self.peaks[index].rt())?;
        Ok(())
    }
    pub fn update_median_rt(&mut self) -> Result<()> {
        self.centroid_rt = self.median(Peak2D::rt, false)?;
        Ok(())
    }
    pub fn update_median_mz(&mut self) -> Result<()> {
        self.centroid_mz = self.median(Peak2D::mz, false)?;
        Ok(())
    }
    fn median(&self, value: impl Fn(&Peak2D) -> f64, intensity: bool) -> Result<f64> {
        self.require_nonempty()?;
        let mut work = self.work()?;
        work.scan(self.len())?;
        if self.len() == 1 {
            return finite(value(&self.peaks[0]));
        }
        work.sort(self.len())?;
        let mut values = work.vector(self.len())?;
        for p in &self.peaks {
            values.push(finite(value(p))?);
        }
        values.sort_unstable_by(|a, b| a.partial_cmp(b).expect("prechecked finite values"));
        let mid = values.len() / 2;
        if values.len() % 2 == 1 {
            return Ok(values[mid]);
        }
        // The private intensity median uses multiplication, coordinate medians division.
        finite(if intensity {
            0.5 * (values[mid - 1] + values[mid])
        } else {
            (values[mid - 1] + values[mid]) / 2.0
        })
    }
    pub fn update_mean_mz(&mut self) -> Result<()> {
        self.require_nonempty()?;
        self.work()?.scan(self.len())?;
        if self.len() == 1 {
            self.centroid_mz = finite(self.peaks[0].mz())?;
            return Ok(());
        }
        let mut sum = 0.0;
        for p in &self.peaks {
            sum = finite(sum + finite(p.mz())?)?;
        }
        self.centroid_mz = finite(sum / self.len() as f64)?;
        Ok(())
    }
    pub fn update_weighted_mean_mz(&mut self) -> Result<()> {
        self.require_nonempty()?;
        self.work()?.scan(self.len())?;
        let result = if self.len() == 1 {
            finite(self.peaks[0].mz())?
        } else {
            let (mut total, mut weighted) = (0.0, 0.0);
            for p in &self.peaks {
                let weight = finite(f64::from(p.intensity))?;
                total = finite(total + weight)?;
                weighted = finite(weighted + weight * finite(p.mz())?)?;
            }
            positive_weight(total)?;
            finite(weighted / total)?
        };
        self.centroid_mz = result;
        Ok(())
    }
    pub fn update_weighted_mz_sd(&mut self) -> Result<()> {
        self.require_nonempty()?;
        self.work()?.scan(self.len())?;
        let (mut total, mut weighted) = (0.0, 0.0);
        for p in &self.peaks {
            let weight = finite(f64::from(p.intensity))?;
            total = finite(total + weight)?;
            // ln(0)=-inf followed by exp(-inf)=0 is valid source arithmetic.
            let square = (2.0 * (finite(p.mz())? - self.centroid_mz).abs().ln()).exp();
            weighted = finite(weighted + weight * square)?;
        }
        positive_weight(total)?;
        self.centroid_sd = finite(weighted.sqrt() / total.sqrt())?;
        Ok(())
    }
    fn work(&self) -> Result<Work> {
        Work::new(self.limits, self.len())
    }
    fn require_nonempty(&self) -> Result<()> {
        if self.is_empty() {
            Err(bad("empty mass trace"))
        } else {
            Ok(())
        }
    }
    fn require_smoothed(&self) -> Result<()> {
        if self.smoothed.is_empty() {
            Err(bad("mass trace has no smoothed intensities"))
        } else {
            Ok(())
        }
    }
    fn value(&self, index: usize, smoothed: bool) -> Result<f64> {
        finite(if smoothed {
            self.smoothed[index]
        } else {
            f64::from(self.peaks[index].intensity)
        })
    }
}
impl Index<usize> for MassTrace {
    type Output = Peak2D;
    fn index(&self, index: usize) -> &Self::Output {
        &self.peaks[index]
    }
}
impl IndexMut<usize> for MassTrace {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.peaks[index]
    }
}
impl<'a> IntoIterator for &'a MassTrace {
    type Item = &'a Peak2D;
    type IntoIter = slice::Iter<'a, Peak2D>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}
impl<'a> IntoIterator for &'a mut MassTrace {
    type Item = &'a mut Peak2D;
    type IntoIter = slice::IterMut<'a, Peak2D>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter_mut()
    }
}

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn resource() -> Error {
    bad("mass trace resource limit exceeded")
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad("mass trace numerical value is nonfinite"))
    }
}
fn positive_weight(value: f64) -> Result<()> {
    if value < f64::EPSILON {
        Err(bad("mass trace total weight is below epsilon"))
    } else {
        Ok(())
    }
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b).ok_or_else(resource)
}
fn interpolate(xa: f64, xb: f64, ya: f64, yb: f64, y: f64) -> Result<f64> {
    if !(ya <= y && y <= yb) {
        return Err(bad("FWHM intensity does not bracket half maximum"));
    }
    if (xa - xb).abs() == 0.0 || (ya - yb).abs() == 0.0 {
        return Ok(xa);
    }
    let x = finite(xa + (y - ya) * (xb - xa) / (yb - ya))?;
    if !(xa.min(xb) <= x && x <= xa.max(xb)) {
        return Err(bad("FWHM interpolation lies outside its bracket"));
    }
    Ok(x)
}
struct Work {
    remaining: usize,
    bytes: usize,
}
impl Work {
    fn new(limits: MassTraceLimits, count: usize) -> Result<Self> {
        if count > limits.max_peaks {
            return Err(resource());
        }
        Ok(Self {
            remaining: limits.max_work,
            bytes: limits.max_bytes,
        })
    }
    fn scan(&mut self, count: usize) -> Result<()> {
        self.remaining = self.remaining.checked_sub(count).ok_or_else(resource)?;
        Ok(())
    }
    fn allocate<T>(&mut self, count: usize) -> Result<()> {
        self.scan(count)?;
        self.bytes = self
            .bytes
            .checked_sub(count.checked_mul(size_of::<T>()).ok_or_else(resource)?)
            .ok_or_else(resource)?;
        Ok(())
    }
    fn vector<T>(&mut self, count: usize) -> Result<Vec<T>> {
        self.allocate::<T>(count)?;
        let mut values = Vec::new();
        values.try_reserve_exact(count).map_err(|_| resource())?;
        Ok(values)
    }
    fn sort(&mut self, count: usize) -> Result<()> {
        if count > 1 {
            let factor = (usize::BITS - count.leading_zeros()) as usize;
            self.scan(
                count
                    .checked_mul(factor)
                    .and_then(|n| n.checked_mul(32))
                    .ok_or_else(resource)?,
            )?;
        }
        Ok(())
    }
}
