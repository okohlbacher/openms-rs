// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Borrowed isotope-group hypotheses and their source-derived summaries.

use crate::{
    Error, Result,
    kernel::{
        ChromatogramPeak, ConvexHull2D, MSChromatogram, MassTrace, MassTraceQuantMethod, Precursor,
    },
    metadata::{ChromatogramType, MetaValue},
};
use std::mem::size_of;

/// Membership and each checked operation's cumulative work/new-payload limits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FeatureHypothesisLimits {
    pub max_traces: usize,
    /// Duplicate trace references count repeatedly, as in source point totals.
    pub max_peaks: usize,
    pub max_work: usize,
    pub max_bytes: usize,
}
impl Default for FeatureHypothesisLimits {
    fn default() -> Self {
        Self {
            max_traces: 1_000_000,
            max_peaks: 10_000_000,
            max_work: 50_000_000,
            max_bytes: 256 * 1024 * 1024,
        }
    }
}

/// Ordered, possibly repeated references to existing mass traces.
///
/// Cached centroids and widths are read as stored. The first trace defines the
/// monoisotopic summary, regardless of its m/z. Clone copies only references;
/// ordinary Clone/Debug/drop have normal Rust costs. Use `checked_clone` to bound
/// the reference-vector copy. Scores retain all source f64 values, including NaN;
/// no Eq/Ord or total ordering is implied.
///
/// A trace cannot be changed or dropped while a hypothesis still borrows it:
/// ```compile_fail
/// use openms::{analysis::feature_hypothesis::FeatureHypothesis, kernel::MassTrace};
/// let mut hypothesis = FeatureHypothesis::new();
/// {
///     let trace = MassTrace::new();
///     hypothesis.add_mass_trace(&trace).unwrap();
/// }
/// assert_eq!(hypothesis.len(), 1);
/// ```
#[derive(Clone, Debug, Default)]
pub struct FeatureHypothesis<'a> {
    traces: Vec<&'a MassTrace>,
    total_peaks: usize,
    score: f64,
    charge: i64,
    pub limits: FeatureHypothesisLimits,
}

impl<'a> FeatureHypothesis<'a> {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn with_limits(limits: FeatureHypothesisLimits) -> Self {
        Self {
            limits,
            ..Self::default()
        }
    }
    pub fn len(&self) -> usize {
        self.traces.len()
    }
    pub fn is_empty(&self) -> bool {
        self.traces.is_empty()
    }
    pub fn traces(&self) -> &[&'a MassTrace] {
        &self.traces
    }
    pub fn score(&self) -> f64 {
        self.score
    }
    pub fn set_score(&mut self, score: f64) {
        self.score = score;
    }
    /// Platform-independent signed source storage. Chromatogram export checks
    /// that this value fits its narrower i32 precursor charge field.
    pub fn charge(&self) -> i64 {
        self.charge
    }
    pub fn set_charge(&mut self, charge: i64) {
        self.charge = charge;
    }
    /// Adds a borrowed reference atomically; the same trace may occur repeatedly.
    pub fn add_mass_trace(&mut self, trace: &'a MassTrace) -> Result<()> {
        let count = add(self.len(), 1)?;
        let total = add(self.total_peaks, trace.len())?;
        self.check_membership(count, total)?;
        let mut work = Work::new(self.limits);
        work.scan(1)?;
        if self.traces.len() == self.traces.capacity() {
            let capacity = self
                .traces
                .capacity()
                .saturating_mul(2)
                .max(count)
                .min(self.limits.max_traces);
            work.allocate::<&MassTrace>(capacity)?;
            work.scan(self.len())?;
            self.traces
                .try_reserve_exact(capacity - self.len())
                .map_err(|_| resource())?;
        }
        self.traces.push(trace);
        // Shared borrowing prevents trace lengths changing during membership.
        self.total_peaks = total;
        Ok(())
    }
    pub fn checked_clone(&self) -> Result<Self> {
        let mut work = self.work()?;
        let mut traces = work.vector(self.len())?;
        traces.extend_from_slice(&self.traces);
        Ok(Self { traces, ..*self })
    }
    pub fn number_of_feature_points(&self) -> usize {
        self.total_peaks
    }
    pub fn centroid_mz(&self) -> Result<f64> {
        Ok(self.first()?.centroid_mz())
    }
    pub fn centroid_rt(&self) -> Result<f64> {
        Ok(self.first()?.centroid_rt())
    }
    pub fn fwhm(&self) -> f64 {
        self.traces.first().map_or(0.0, |trace| trace.fwhm())
    }
    pub fn labels(&self) -> Result<Vec<String>> {
        let mut work = self.work()?;
        let mut labels = work.vector(self.len())?;
        for trace in &self.traces {
            labels.push(work.text(trace.label())?);
        }
        Ok(labels)
    }
    /// Underscore concatenation retains empty labels and insertion order.
    pub fn label(&self) -> Result<String> {
        let mut work = self.work()?;
        let mut bytes = self.len().saturating_sub(1);
        for trace in &self.traces {
            bytes = add(bytes, trace.label().len())?;
        }
        work.scan(self.len())?;
        let mut label = work.string(bytes)?;
        for (i, trace) in self.traces.iter().enumerate() {
            if i != 0 {
                label.push('_');
            }
            label.push_str(trace.label());
        }
        Ok(label)
    }
    pub fn all_centroid_mz(&self) -> Result<Vec<f64>> {
        self.cached_values(MassTrace::centroid_mz)
    }
    pub fn all_centroid_rt(&self) -> Result<Vec<f64>> {
        self.cached_values(MassTrace::centroid_rt)
    }
    /// Returns every cached value, including zero for traces without an IM flag.
    pub fn all_centroid_im(&self) -> Result<Vec<f64>> {
        self.cached_values(MassTrace::centroid_im)
    }
    pub fn isotope_distances(&self) -> Result<Vec<f64>> {
        let mut work = self.work()?;
        let mut distances = work.vector(self.len().saturating_sub(1))?;
        for pair in self.traces.windows(2) {
            distances.push(finite(pair[1].centroid_mz() - pair[0].centroid_mz())?);
        }
        Ok(distances)
    }
    pub fn monoisotopic_feature_intensity(&self, smoothed: bool) -> Result<f64> {
        let trace = self.first()?;
        self.check_membership(1, trace.len())?;
        let mut work = Work::new(self.limits);
        work.intensity(trace, smoothed)?;
        trace.intensity(smoothed)
    }
    pub fn all_intensities(&self, smoothed: bool) -> Result<Vec<f64>> {
        let mut work = self.work()?;
        let mut values = work.vector(self.len())?;
        for trace in &self.traces {
            work.intensity(trace, smoothed)?;
            values.push(trace.intensity(smoothed)?);
        }
        Ok(values)
    }
    pub fn summed_feature_intensity(&self, smoothed: bool) -> Result<f64> {
        let mut work = self.work()?;
        let mut sum = 0.0;
        for trace in &self.traces {
            work.intensity(trace, smoothed)?;
            sum = finite(sum + trace.intensity(smoothed)?)?;
        }
        Ok(sum)
    }
    pub fn max_intensity(&self, smoothed: bool) -> Result<f64> {
        let mut work = self.work()?;
        let mut highest = 0.0;
        for trace in &self.traces {
            work.scan(if smoothed {
                trace.smoothed_intensities().len()
            } else {
                trace.len()
            })?;
            let value = trace.max_intensity(smoothed)?;
            if value > highest {
                highest = value;
            }
        }
        Ok(highest)
    }
    /// One source scan-envelope hull per trace, using all raw RT/m/z positions.
    pub fn convex_hulls(&self) -> Result<Vec<ConvexHull2D>> {
        let mut work = self.work()?;
        let mut hulls = work.vector(self.len())?;
        for trace in &self.traces {
            // Precharge the existing MassTrace/geometry path, including scans,
            // input points, minimum Vec growth, stable-sort scratch and comparisons,
            // before calling it (128 bytes/point also covers singleton Scan capacity).
            work.sort(trace.len())?;
            work.scan(mul(trace.len(), 4)?)?;
            work.allocate::<[f64; 16]>(trace.len())?;
            hulls.push(trace.convex_hull()?);
        }
        Ok(hulls)
    }
    /// Raw, RT-sorted BasePeak chromatograms. Each precursor uses the first
    /// trace's cached m/z, the hypothesis charge and `peptide_sequence` ID text.
    /// Equal RTs keep their input order. Empty hypotheses and charges outside
    /// i32 are checked errors replacing unsafe/implementation-defined source cases.
    pub fn chromatograms(&self, feature_id: u64) -> Result<Vec<MSChromatogram>> {
        let mz = finite(self.first()?.centroid_mz())?;
        let charge = i32::try_from(self.charge).map_err(|_| bad("precursor charge exceeds i32"))?;
        let mut work = self.work()?;
        let mut chromatograms = work.vector(self.len())?;
        // u64 and usize decimal strings plus underscore fit this fixed allowance.
        work.allocate::<u8>(64)?;
        let id = feature_id.to_string();
        for (i, trace) in self.traces.iter().enumerate() {
            work.scan(mul(trace.len(), 2)?)?;
            work.sort(trace.len())?;
            work.allocate::<ChromatogramPeak>(trace.len())?; // stable-sort scratch
            let mut peaks = work.vector(trace.len())?;
            for peak in trace.peaks() {
                peaks.push(ChromatogramPeak::new(
                    finite(peak.rt())?,
                    finite(f64::from(peak.intensity))? as f32,
                ));
            }
            peaks.sort_by(|a, b| a.rt.partial_cmp(&b.rt).expect("finite RT preflight"));
            work.allocate::<u8>(64)?;
            let name = format!("{feature_id}_{i}");
            let native_id = work.text(&name)?;
            // One sparse BTree node plus owned key/value strings; account before
            // constructing it. No source/acquisition/processing payload is copied.
            work.allocate::<u8>(add(512, mul(24, size_of::<(String, MetaValue)>())?)?)?;
            let key = work.text("peptide_sequence")?;
            let value = work.text(&id)?;
            let mut precursor = Precursor::new(mz, charge);
            precursor
                .cv_terms
                .metadata
                .insert(key, MetaValue::from(value));
            chromatograms.push(MSChromatogram {
                peaks,
                native_id,
                name,
                precursor,
                chromatogram_type: ChromatogramType::BasePeak,
                ..MSChromatogram::default()
            });
        }
        Ok(chromatograms)
    }
    fn cached_values(&self, value: impl Fn(&MassTrace) -> f64) -> Result<Vec<f64>> {
        let mut work = self.work()?;
        let mut values = work.vector(self.len())?;
        values.extend(self.traces.iter().map(|trace| value(trace)));
        Ok(values)
    }
    fn first(&self) -> Result<&MassTrace> {
        self.traces
            .first()
            .copied()
            .ok_or_else(|| bad("hypothesis contains no mass traces"))
    }
    fn check_membership(&self, count: usize, peaks: usize) -> Result<()> {
        if count > self.limits.max_traces || peaks > self.limits.max_peaks {
            Err(resource())
        } else {
            Ok(())
        }
    }
    fn work(&self) -> Result<Work> {
        self.check_membership(self.len(), self.total_peaks)?;
        let mut work = Work::new(self.limits);
        work.scan(self.len())?;
        Ok(work)
    }
}

/// Source CmpMassTraceByMZ predicate; ordinary f64 `<`, not a total ordering.
pub fn mass_trace_mz_less(left: &MassTrace, right: &MassTrace) -> bool {
    left.centroid_mz() < right.centroid_mz()
}
/// Source CmpHypothesesByScore predicate; descending score, unordered NaNs false.
pub fn hypothesis_score_greater(
    left: &FeatureHypothesis<'_>,
    right: &FeatureHypothesis<'_>,
) -> bool {
    left.score() > right.score()
}
/// Native name for the source's public two-field FeatureFindingMetabo Range.
/// It is an unvalidated value record; no generic range algebra is implied.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MetaboIsotopeMassWindow {
    pub left_boundary: f64,
    pub right_boundary: f64,
}

struct Work {
    remaining: usize,
    bytes: usize,
}
impl Work {
    fn new(limits: FeatureHypothesisLimits) -> Self {
        Self {
            remaining: limits.max_work,
            bytes: limits.max_bytes,
        }
    }
    fn scan(&mut self, count: usize) -> Result<()> {
        self.remaining = self.remaining.checked_sub(count).ok_or_else(resource)?;
        Ok(())
    }
    fn allocate<T>(&mut self, count: usize) -> Result<()> {
        self.scan(count)?;
        self.bytes = self
            .bytes
            .checked_sub(mul(count, size_of::<T>())?)
            .ok_or_else(resource)?;
        Ok(())
    }
    fn vector<T>(&mut self, count: usize) -> Result<Vec<T>> {
        self.allocate::<T>(count)?;
        let mut values = Vec::new();
        values.try_reserve_exact(count).map_err(|_| resource())?;
        Ok(values)
    }
    fn string(&mut self, count: usize) -> Result<String> {
        self.allocate::<u8>(count)?;
        let mut value = String::new();
        value.try_reserve_exact(count).map_err(|_| resource())?;
        Ok(value)
    }
    fn text(&mut self, value: &str) -> Result<String> {
        let mut copy = self.string(value.len())?;
        copy.push_str(value);
        Ok(copy)
    }
    fn sort(&mut self, count: usize) -> Result<()> {
        if count > 1 {
            self.scan(mul(
                mul(count, (usize::BITS - count.leading_zeros()) as usize)?,
                32,
            )?)?;
        }
        Ok(())
    }
    fn intensity(&mut self, trace: &MassTrace, smoothed: bool) -> Result<()> {
        match trace.quant_method() {
            MassTraceQuantMethod::Median => {
                self.scan(trace.len())?;
                if trace.len() > 1 {
                    self.sort(trace.len())?;
                    self.allocate::<f64>(trace.len())?;
                }
            }
            MassTraceQuantMethod::Area => {
                let (left, right) = trace.fwhm_borders();
                if (left, right) != (0, 0) {
                    self.scan(right - left + 1)?;
                }
            }
            MassTraceQuantMethod::MaxHeight => self.scan(if smoothed {
                trace.smoothed_intensities().len()
            } else {
                trace.len()
            })?,
        }
        Ok(())
    }
}
fn add(left: usize, right: usize) -> Result<usize> {
    left.checked_add(right).ok_or_else(resource)
}
fn mul(left: usize, right: usize) -> Result<usize> {
    left.checked_mul(right).ok_or_else(resource)
}
fn resource() -> Error {
    bad("feature hypothesis resource limit or allocation failure")
}
fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad(
            "feature hypothesis consumed nonfinite data or arithmetic",
        ))
    }
}
