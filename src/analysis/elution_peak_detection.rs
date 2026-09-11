// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Source elution splitting, smoothing, extrema and explicit width filtering.

use crate::concept::progress_logger::{ProgressLogType, ProgressLogger};
use crate::kernel::{MassTrace, Peak2D};
use crate::processing::smoothing::SavitzkyGolayFilter;
use crate::{Error, Result};
use std::mem::size_of;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ElutionPeakWidthFiltering {
    Off,
    #[default]
    Fixed,
    Auto,
}

/// All six source parameters. `min_fwhm` also controls valley separation when
/// width filtering is Off/Auto; Auto requires an explicit width-filter call.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElutionPeakDetectionOptions {
    pub chrom_fwhm: f64,
    pub chrom_peak_snr: f64,
    pub width_filtering: ElutionPeakWidthFiltering,
    pub min_fwhm: f64,
    pub max_fwhm: f64,
    pub masstrace_snr_filtering: bool,
}
impl Default for ElutionPeakDetectionOptions {
    fn default() -> Self {
        Self {
            chrom_fwhm: 5.0,
            chrom_peak_snr: 3.0,
            width_filtering: ElutionPeakWidthFiltering::Fixed,
            min_fwhm: 1.0,
            max_fwhm: 60.0,
            masstrace_snr_filtering: false,
        }
    }
}

/// One operation's aggregate input/output counts, work and allocated payload.
/// Input and generated peak counts are each capped; all copies share bytes/work.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ElutionPeakDetectionLimits {
    pub max_traces: usize,
    pub max_peaks: usize,
    pub max_work: usize,
    pub max_bytes: usize,
}
impl Default for ElutionPeakDetectionLimits {
    fn default() -> Self {
        Self {
            max_traces: 1_000_000,
            max_peaks: 10_000_000,
            max_work: 50_000_000,
            max_bytes: 256 * 1024 * 1024,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ElutionExtrema {
    pub maxima: Vec<usize>,
    pub minima: Vec<usize>,
}

/// Input traces and returned output change only after the whole call succeeds.
/// Successful rejection still publishes source smoothing/FWHM changes. Progress
/// backend output and state cannot be rolled back with scientific results.
#[derive(Clone)]
pub struct ElutionPeakDetection {
    pub options: ElutionPeakDetectionOptions,
    pub limits: ElutionPeakDetectionLimits,
    pub logger: ProgressLogger,
}
impl Default for ElutionPeakDetection {
    fn default() -> Self {
        let mut logger = ProgressLogger::new();
        logger.set_log_type(ProgressLogType::Cmd);
        Self {
            options: ElutionPeakDetectionOptions::default(),
            limits: ElutionPeakDetectionLimits::default(),
            logger,
        }
    }
}
impl ElutionPeakDetection {
    pub fn new() -> Self {
        Self::default()
    }
    /// Source single-trace overload: no progress events.
    pub fn detect_peaks(&mut self, trace: &mut MassTrace) -> Result<Vec<MassTrace>> {
        let mut work = Work::new(self.limits);
        work.inputs(std::slice::from_ref(trace))?;
        let mut staged = work.copy_trace(trace)?;
        let mut output = Vec::new();
        self.detect(&mut staged, &mut output, &mut work)?;
        *trace = staged;
        Ok(output)
    }
    /// Returns previous output ownership without inspecting or destroying it.
    pub fn detect_peaks_into(
        &mut self,
        trace: &mut MassTrace,
        output: &mut Vec<MassTrace>,
    ) -> Result<Vec<MassTrace>> {
        Ok(std::mem::replace(output, self.detect_peaks(trace)?))
    }
    /// Serial input/segment order replaces the source's optional OpenMP order.
    pub fn detect_peaks_many(&mut self, traces: &mut [MassTrace]) -> Result<Vec<MassTrace>> {
        let mut work = Work::new(self.limits);
        let mut staged = work.copy_inputs(traces)?;
        let count = i64::try_from(traces.len()).map_err(|_| resource())?;
        self.logger
            .start_progress(0, count, "elution peak detection")?;
        let mut output = Vec::new();
        let result: Result<()> = (|| {
            for (i, trace) in staged.iter_mut().enumerate() {
                self.logger.set_progress(i as i64)?;
                self.detect(trace, &mut output, &mut work)?;
            }
            Ok(())
        })();
        let end = self.logger.end_progress(0);
        result?;
        end?;
        for (target, source) in traces.iter_mut().zip(staged) {
            *target = source;
        }
        Ok(output)
    }
    pub fn detect_peaks_many_into(
        &mut self,
        traces: &mut [MassTrace],
        output: &mut Vec<MassTrace>,
    ) -> Result<Vec<MassTrace>> {
        Ok(std::mem::replace(output, self.detect_peaks_many(traces)?))
    }
    /// Recomputes every input FWHM, then retains inclusive sorted ranks
    /// floor(n*0.05)..=floor(n*0.95), independently of width_filtering.
    pub fn filter_by_peak_width(&self, traces: &mut [MassTrace]) -> Result<Vec<MassTrace>> {
        let mut work = Work::new(self.limits);
        let mut staged = work.copy_inputs(traces)?;
        let mut widths = work.vector(staged.len())?;
        for (index, trace) in staged.iter_mut().enumerate() {
            work.summary(trace.len())?;
            widths.push((trace.estimate_fwhm(true)?, index));
        }
        work.sort(widths.len())?;
        widths.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.cmp(&b.1)));
        let lower = (widths.len() as f64 * 0.05).floor() as usize;
        let upper = (widths.len() as f64 * 0.95).floor() as usize;
        let mut output = Vec::new();
        for (rank, &(_, index)) in widths.iter().enumerate() {
            work.consume(1)?;
            if rank >= lower && rank <= upper {
                let copy = work.copy_trace(&staged[index])?;
                work.push(&mut output, copy)?;
            }
        }
        for (target, source) in traces.iter_mut().zip(staged) {
            *target = source;
        }
        Ok(output)
    }
    pub fn filter_by_peak_width_into(
        &self,
        traces: &mut [MassTrace],
        output: &mut Vec<MassTrace>,
    ) -> Result<Vec<MassTrace>> {
        Ok(std::mem::replace(
            output,
            self.filter_by_peak_width(traces)?,
        ))
    }
    /// Degree two, even windows incremented, with the inherited odd-frame cap
    /// 1023. Fewer than three peaks ignore the window and copy raw intensities.
    pub fn smooth_data(&self, trace: &mut MassTrace, window: i32) -> Result<()> {
        let mut work = Work::new(self.limits);
        work.input(trace)?;
        smooth(trace, window, &mut work)
    }
    pub fn find_local_extrema(
        &self,
        trace: &MassTrace,
        neighbors: usize,
    ) -> Result<ElutionExtrema> {
        let mut work = Work::new(self.limits);
        work.input(trace)?;
        extrema(trace, neighbors, self.options.min_fwhm, &mut work)
    }
    pub fn compute_mass_trace_noise(&self, trace: &MassTrace) -> Result<f64> {
        let mut work = Work::new(self.limits);
        work.input(trace)?;
        noise(trace, &mut work)
    }
    /// Nonempty zero-noise/zero-span divisions error when their result is not
    /// finite. Empty traces return zero, preserving the separate source branch.
    pub fn compute_mass_trace_snr(&self, trace: &MassTrace) -> Result<f64> {
        let mut work = Work::new(self.limits);
        work.input(trace)?;
        if trace.is_empty() {
            return Ok(0.0);
        }
        let noise = noise(trace, &mut work)?;
        work.summary(trace.len())?;
        let noise_area = finite(noise * trace.trace_length()?)?;
        finite(trace.compute_peak_area()? / noise_area)
    }
    pub fn compute_apex_snr(&self, trace: &MassTrace) -> Result<f64> {
        let mut work = Work::new(self.limits);
        work.input(trace)?;
        apex_snr(trace, &mut work)
    }
    fn detect(
        &self,
        trace: &mut MassTrace,
        output: &mut Vec<MassTrace>,
        work: &mut Work,
    ) -> Result<()> {
        work.consume(8)?;
        let window =
            finite(finite(self.options.chrom_fwhm)? / trace.average_ms1_cycle_time()?)?.ceil();
        // Source converts double -> Size -> Int. Reject undefined/narrowing
        // boundaries rather than platform-dependent saturation or wrapping.
        if window < 0.0 || window > f64::from(i32::MAX) {
            return Err(bad(
                "elution window cannot be represented as a nonnegative i32",
            ));
        }
        let window = window as i32;
        smooth(trace, window, work)?;
        let found = extrema(trace, window as usize / 2, self.options.min_fwhm, work)?;
        if found.maxima.is_empty() {
            return Ok(());
        }
        if found.maxima.len() == 1 {
            let width_ok = self.width_ok(trace, work)?;
            let snr_ok = self.snr_ok(trace, work)?;
            if width_ok && snr_ok {
                work.summary(trace.len())?;
                trace.update_smoothed_max_rt()?;
                if self.options.width_filtering != ElutionPeakWidthFiltering::Fixed {
                    trace.estimate_fwhm(true)?;
                }
                let copy = work.copy_trace(trace)?;
                work.push(output, copy)?;
            }
            return Ok(());
        }
        let mut first = 0;
        for (ordinal, last) in found
            .minima
            .into_iter()
            .chain(std::iter::once(trace.len() - 1))
            .enumerate()
        {
            work.consume(1)?;
            let end = add(last, 1)?;
            let source = trace
                .peaks()
                .get(first..end)
                .ok_or_else(|| bad("invalid source split interval"))?;
            let mut peaks = work.vector(source.len())?;
            peaks.extend_from_slice(source);
            work.summary(source.len())?;
            let mut split = MassTrace::from_peaks_with_limits(peaks, trace.limits)?;
            work.allocate(mul(source.len(), size_of::<f64>())?)?;
            split.set_smoothed_intensities(&trace.smoothed_intensities()[first..end])?;
            first = end;
            if trace.contains_im_data() {
                split.set_centroid_im(trace.centroid_im())?;
            }
            split.fwhm_mz_avg = trace.fwhm_mz_avg;
            split.fwhm_im_avg = trace.fwhm_im_avg;
            let width_ok = self.width_ok(&mut split, work)?;
            // Source intentionally evaluates the original trace for every split.
            let snr_ok = self.snr_ok(trace, work)?;
            if width_ok && snr_ok {
                let label_len = add(trace.label().len(), 32)?;
                work.consume(mul(label_len, 2)?)?;
                work.allocate(mul(label_len, 2)?)?;
                split.set_label(&format!("{}.{}", trace.label(), ordinal + 1))?;
                work.summary(split.len())?;
                split.update_smoothed_max_rt()?;
                split.update_weighted_mean_mz()?;
                split.update_weighted_mz_sd()?;
                split.set_quant_method(trace.quant_method());
                if self.options.width_filtering != ElutionPeakWidthFiltering::Fixed {
                    split.estimate_fwhm(true)?;
                }
                work.push(output, split)?;
            }
        }
        Ok(())
    }
    fn width_ok(&self, trace: &mut MassTrace, work: &mut Work) -> Result<bool> {
        if self.options.width_filtering != ElutionPeakWidthFiltering::Fixed {
            return Ok(true);
        }
        let minimum = finite(self.options.min_fwhm)?;
        let maximum = finite(self.options.max_fwhm)?;
        work.summary(trace.len())?;
        let width = trace.estimate_fwhm(true)?;
        Ok(!(width < minimum || width > maximum))
    }
    fn snr_ok(&self, original: &MassTrace, work: &mut Work) -> Result<bool> {
        if !self.options.masstrace_snr_filtering {
            return Ok(true);
        }
        let threshold = finite(self.options.chrom_peak_snr)?;
        Ok(apex_snr(original, work)? >= threshold)
    }
}

fn smooth(trace: &mut MassTrace, window: i32, work: &mut Work) -> Result<()> {
    let n = trace.len();
    let mut raw = work.vector(n)?;
    for peak in trace.peaks() {
        raw.push(finite(f64::from(peak.intensity))?);
    }
    let values = if n < 3 {
        raw
    } else {
        let frame = window.max(3) as usize;
        let frame = add(frame, usize::from(frame % 2 == 0))?;
        if frame > 1023 {
            return Err(bad("Savitzky–Golay odd frame exceeds 1023"));
        }
        // Covers degree-two QR, coefficient rows, their Vec slots, and scratch
        // before the existing smoother allocates. No independent per-trace reset.
        work.consume(add(mul(mul(frame, frame)?, 16)?, mul(frame, 256)?)?)?;
        work.allocate(add(mul(mul(frame, frame / 2 + 1)?, 8)?, mul(frame, 128)?)?)?;
        let mut rts = work.vector(n)?;
        for peak in trace.peaks() {
            rts.push(peak.rt());
        }
        work.consume(mul(n, if frame > n { 4 } else { add(frame, 4)? })?)?;
        work.allocate(mul(n, 8)?)?;
        let mut filtered = SavitzkyGolayFilter::new(frame, 2)?.filter(&rts, &raw)?;
        work.consume(n)?;
        for value in &mut filtered {
            // Source spectrum intensity storage rounds to f32 before promotion.
            *value = finite(f64::from(*value as f32))?;
        }
        filtered
    };
    work.consume(mul(n, 3)?)?;
    work.allocate(mul(n, 8)?)?;
    trace.set_smoothed_intensities(&values)
}

fn extrema(
    trace: &MassTrace,
    neighbors: usize,
    min_fwhm: f64,
    work: &mut Work,
) -> Result<ElutionExtrema> {
    let values = trace.smoothed_intensities();
    let n = values.len();
    if n != trace.len() {
        return Err(bad("mass trace must be smoothed before extrema detection"));
    }
    if n == 0 {
        return Ok(ElutionExtrema::default());
    }
    work.consume(n)?;
    for &value in values {
        finite(value)?;
    }
    let mut maxima = work.vector(n)?;
    let mut minima = work.vector(n)?;
    if n < 3 {
        maxima.push(usize::from(n == 2 && values[1] > values[0]));
        return Ok(ElutionExtrema { maxima, minima });
    }
    // Precheck all source endpoint additions before they can wrap.
    add(n - 1, neighbors)?;
    let mut indices = work.vector(n)?;
    indices.extend(0..n);
    work.sort(n)?;
    indices.sort_unstable_by(|&a, &b| values[a].partial_cmp(&values[b]).unwrap().then(a.cmp(&b)));
    let mut used = work.vector(n)?;
    used.resize(n, false);
    for index in indices {
        work.consume(1)?;
        let height = values[index];
        if used[index] || height <= 0.0 {
            continue;
        }
        let start = index.saturating_sub(neighbors);
        let end = (index + neighbors).min(n);
        let mut maximum = true;
        for &value in &values[start..end] {
            work.consume(1)?;
            if value > height {
                maximum = false;
                break;
            }
        }
        if maximum {
            maxima.push(index);
            work.consume(end - start)?;
            used[start..end].fill(true);
        }
    }
    work.sort(maxima.len())?;
    maxima.sort_unstable();
    if maxima.len() > 1 {
        let distance = finite(min_fwhm)? / 2.0;
        let (mut left, mut right) = (0, 1);
        while right < maxima.len() {
            work.consume(16)?;
            let mut low = maxima[left] + 1;
            let mut high = maxima[right] - 1;
            while low + 1 < high {
                work.consume(2)?;
                let mid = low + (high - low) / 2;
                if values[mid] <= values[mid + 1] {
                    high = mid;
                } else {
                    low = mid;
                }
            }
            let valley = if values[low] < values[high] {
                low
            } else {
                high
            };
            let valley_height = values[valley].max(1.0);
            let left_height = values[maxima[left]];
            let right_height = values[maxima[right]];
            let left_rt = finite(trace[maxima[left]].rt())?;
            let middle_rt = finite(trace[valley].rt())?;
            let right_rt = finite(trace[maxima[right]].rt())?;
            if left_height / valley_height >= 2.0
                && right_height / valley_height >= 2.0
                && finite((middle_rt - left_rt).abs())? >= distance
                && finite((right_rt - middle_rt).abs())? >= distance
            {
                minima.push(valley);
                left = right;
            } else if left_height <= right_height {
                left = right;
            }
            right += 1;
        }
    }
    Ok(ElutionExtrema { maxima, minima })
}
fn noise(trace: &MassTrace, work: &mut Work) -> Result<f64> {
    let values = trace.smoothed_intensities();
    work.consume(mul(values.len(), 4)?)?;
    let mut sum = 0.0;
    for (peak, &value) in trace.peaks().iter().zip(values) {
        let difference = finite(f64::from(peak.intensity) - value)?;
        sum = finite(sum + difference * difference)?;
    }
    if values.is_empty() {
        Ok(0.0)
    } else {
        finite((sum / values.len() as f64).sqrt())
    }
}
fn apex_snr(trace: &MassTrace, work: &mut Work) -> Result<f64> {
    let noise = noise(trace, work)?;
    if noise > 0.0 {
        work.summary(trace.len())?;
        finite(trace.max_intensity(true)? / noise)
    } else {
        Ok(0.0)
    }
}

struct Work {
    limits: ElutionPeakDetectionLimits,
    remaining: usize,
    bytes: usize,
    output_peaks: usize,
}
impl Work {
    fn new(limits: ElutionPeakDetectionLimits) -> Self {
        Self {
            limits,
            remaining: limits.max_work,
            bytes: limits.max_bytes,
            output_peaks: 0,
        }
    }
    fn consume(&mut self, n: usize) -> Result<()> {
        self.remaining = self.remaining.checked_sub(n).ok_or_else(resource)?;
        Ok(())
    }
    fn allocate(&mut self, n: usize) -> Result<()> {
        self.bytes = self.bytes.checked_sub(n).ok_or_else(resource)?;
        Ok(())
    }
    fn vector<T>(&mut self, n: usize) -> Result<Vec<T>> {
        self.consume(n)?;
        self.allocate(mul(n, size_of::<T>())?)?;
        let mut result = Vec::new();
        result.try_reserve_exact(n).map_err(|_| resource())?;
        Ok(result)
    }
    fn sort(&mut self, n: usize) -> Result<()> {
        let logarithm = usize::BITS as usize - n.leading_zeros() as usize;
        self.consume(mul(mul(n, logarithm)?, 32)?)
    }
    fn summary(&mut self, n: usize) -> Result<()> {
        self.consume(add(mul(n, 64)?, 256)?)
    }
    fn input(&mut self, trace: &MassTrace) -> Result<()> {
        limit(trace.len(), self.limits.max_peaks)?;
        self.consume(1)
    }
    fn inputs(&mut self, traces: &[MassTrace]) -> Result<()> {
        limit(traces.len(), self.limits.max_traces)?;
        self.consume(traces.len())?;
        let mut total = 0;
        for trace in traces {
            total = add(total, trace.len())?;
            limit(total, self.limits.max_peaks)?;
        }
        Ok(())
    }
    fn copy_trace(&mut self, trace: &MassTrace) -> Result<MassTrace> {
        let n = add(
            mul(trace.len(), size_of::<Peak2D>())?,
            mul(trace.smoothed_intensities().len(), 8)?,
        )?;
        let n = add(add(n, trace.label().len())?, size_of::<MassTrace>())?;
        self.consume(n)?;
        self.allocate(n)?;
        Ok(trace.clone())
    }
    fn copy_inputs(&mut self, traces: &[MassTrace]) -> Result<Vec<MassTrace>> {
        self.inputs(traces)?;
        let mut result = self.vector(traces.len())?;
        for trace in traces {
            result.push(self.copy_trace(trace)?);
        }
        Ok(result)
    }
    fn push(&mut self, output: &mut Vec<MassTrace>, trace: MassTrace) -> Result<()> {
        limit(add(output.len(), 1)?, self.limits.max_traces)?;
        let count = add(self.output_peaks, trace.len())?;
        limit(count, self.limits.max_peaks)?;
        if output.len() == output.capacity() {
            let capacity = output
                .capacity()
                .saturating_mul(2)
                .clamp(1, self.limits.max_traces);
            self.consume(capacity)?;
            self.allocate(mul(capacity, size_of::<MassTrace>())?)?;
            output
                .try_reserve_exact(capacity - output.len())
                .map_err(|_| resource())?;
        }
        self.output_peaks = count;
        output.push(trace);
        Ok(())
    }
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b).ok_or_else(resource)
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b).ok_or_else(resource)
}
fn limit(n: usize, maximum: usize) -> Result<()> {
    if n > maximum { Err(resource()) } else { Ok(()) }
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad("nonfinite elution calculation"))
    }
}
fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn resource() -> Error {
    bad("elution peak detection resource limit exceeded")
}
