// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Apex-first mass-trace extension, including source IM and FWHM branches.

use crate::concept::{
    constants::user_param,
    progress_logger::{ProgressLogType, ProgressLogger},
};
use crate::kernel::{
    AreaIter, MSExperiment, MSSpectrum, MassTrace, MassTraceLimits, MassTraceQuantMethod, Peak1D,
    Peak2D,
};
use crate::{Error, Result};
use std::mem::size_of;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum TraceTerminationCriterion {
    #[default]
    Outlier,
    SampleRate,
}

/// All eleven source parameters. Finite negative scalar settings are retained;
/// a negative maximum length disables that check, just as in the source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MassTraceDetectionOptions {
    pub mass_error_ppm: f64,
    pub noise_threshold_int: f64,
    pub chrom_peak_snr: f64,
    pub ion_mobility_tolerance: f64,
    pub reestimate_mt_sd: bool,
    pub quant_method: MassTraceQuantMethod,
    pub trace_termination_criterion: TraceTerminationCriterion,
    pub trace_termination_outliers: usize,
    pub min_sample_rate: f64,
    pub min_trace_length: f64,
    pub max_trace_length: f64,
}
impl Default for MassTraceDetectionOptions {
    fn default() -> Self {
        Self {
            mass_error_ppm: 20.0,
            noise_threshold_int: 10.0,
            chrom_peak_snr: 3.0,
            ion_mobility_tolerance: 0.01,
            reestimate_mt_sd: true,
            quant_method: MassTraceQuantMethod::Area,
            trace_termination_criterion: TraceTerminationCriterion::Outlier,
            trace_termination_outliers: 5,
            min_sample_rate: 0.5,
            min_trace_length: 5.0,
            max_trace_length: -1.0,
        }
    }
}
impl MassTraceDetectionOptions {
    fn validate(self) -> Result<()> {
        for value in [
            self.mass_error_ppm,
            self.noise_threshold_int,
            self.chrom_peak_snr,
            self.ion_mobility_tolerance,
            self.min_sample_rate,
            self.min_trace_length,
            self.max_trace_length,
        ] {
            finite(value)?;
        }
        Ok(())
    }
}

/// Cumulative limits for one operation, including area reconstruction, failed
/// candidate traces, sorting, and accepted-trace summary calculations.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MassTraceDetectionLimits {
    pub max_spectra: usize,
    pub max_peaks: usize,
    pub max_traces: usize,
    pub max_work: usize,
    pub max_bytes: usize,
}
impl Default for MassTraceDetectionLimits {
    fn default() -> Self {
        Self {
            max_spectra: 1_000_000,
            max_peaks: 10_000_000,
            max_traces: 1_000_000,
            max_work: 50_000_000,
            max_bytes: 256 * 1024 * 1024,
        }
    }
}

/// Configuration and last successful run's availability flags. Progress output
/// is caller-controlled and cannot be rolled back with scientific output.
#[derive(Clone)]
pub struct MassTraceDetection {
    pub options: MassTraceDetectionOptions,
    pub limits: MassTraceDetectionLimits,
    pub logger: ProgressLogger,
    arrays: Arrays,
    ccs_warning: bool,
}
impl Default for MassTraceDetection {
    fn default() -> Self {
        let mut logger = ProgressLogger::new();
        logger.set_log_type(ProgressLogType::Cmd);
        Self {
            options: MassTraceDetectionOptions::default(),
            limits: MassTraceDetectionLimits::default(),
            logger,
            arrays: Arrays::default(),
            ccs_warning: false,
        }
    }
}
impl MassTraceDetection {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn has_fwhm_mz(&self) -> bool {
        self.arrays.mz.is_some()
    }
    pub fn has_fwhm_im(&self) -> bool {
        self.arrays.im_width.is_some()
    }
    pub fn has_centroid_im(&self) -> bool {
        self.arrays.im.is_some()
    }
    /// First source-recognized IM array has CCS units and tolerance is below 1.
    /// This diagnostic does not change the exact-name scientific IM selection.
    pub fn ccs_tolerance_warning(&self) -> bool {
        self.ccs_warning
    }
    /// Zero means unlimited traces subject to the operation's resource limits.
    pub fn run(&mut self, input: &MSExperiment, max_traces: usize) -> Result<Vec<MassTrace>> {
        let mut work = Work::new(self.limits);
        let (output, arrays, warning) = self.execute(input, max_traces, &mut work)?;
        self.arrays = arrays;
        self.ccs_warning = warning;
        Ok(output)
    }
    /// Replaces output on success and returns its prior ownership. The previous
    /// output is never scanned, cloned, validated or destroyed by this method.
    pub fn run_into(
        &mut self,
        input: &MSExperiment,
        output: &mut Vec<MassTrace>,
        max_traces: usize,
    ) -> Result<Vec<MassTrace>> {
        Ok(std::mem::replace(output, self.run(input, max_traces)?))
    }
    /// Source area overload: equal endpoints are a complete no-op. Other calls
    /// stage the consumed cursor and replacement output until all work succeeds.
    /// A nonreachable endpoint errors; source dereferencing end was undefined.
    pub fn run_area<'a>(
        &mut self,
        begin: &mut AreaIter<'a>,
        end: &AreaIter<'a>,
        output: &mut Vec<MassTrace>,
    ) -> Result<Option<Vec<MassTrace>>> {
        if begin == end {
            return Ok(None);
        }
        let mut work = Work::new(self.limits);
        let mut cursor = begin.clone();
        let mut input = MSExperiment::default();
        let mut spectrum = MSSpectrum::default();
        let mut count = 0usize;
        while &cursor != end {
            work.consume(1)?;
            let item = cursor
                .next()
                .ok_or_else(|| bad("area endpoint is not reachable"))?;
            count = add(count, 1)?;
            limit(count, self.limits.max_peaks, "area peak count")?;
            if item.spectrum.rt != spectrum.rt {
                if spectrum.rt != -1.0 {
                    limit(
                        add(input.spectra.len(), 1)?,
                        self.limits.max_spectra,
                        "area scan count",
                    )?;
                    work.push(&mut input.spectra, spectrum)?;
                    spectrum = MSSpectrum::default();
                } else {
                    // Source clear(false) discards a preceding sentinel-RT group.
                    work.consume(spectrum.peaks.len())?;
                    spectrum.peaks.clear();
                }
                spectrum.rt = item.spectrum.rt;
            }
            work.push(&mut spectrum.peaks, *item.peak)?;
        }
        limit(
            add(input.spectra.len(), 1)?,
            self.limits.max_spectra,
            "area scan count",
        )?;
        work.push(&mut input.spectra, spectrum)?;
        let (replacement, arrays, warning) = self.execute(&input, 0, &mut work)?;
        self.arrays = arrays;
        self.ccs_warning = warning;
        *begin = cursor;
        Ok(Some(std::mem::replace(output, replacement)))
    }

    fn execute(
        &mut self,
        input: &MSExperiment,
        max_traces: usize,
        work: &mut Work,
    ) -> Result<(Vec<MassTrace>, Arrays, bool)> {
        self.options.validate()?;
        limit(
            input.spectra.len(),
            self.limits.max_spectra,
            "spectrum count",
        )?;
        work.consume(input.spectra.len())?;
        let warning = ccs_warning(input, self.options.ion_mobility_tolerance, work)?;
        let mut scans = Vec::new();
        let mut apices = Vec::new();
        let mut total = 0usize;
        let mut input_peaks = 0usize;
        let apex_threshold =
            finite(self.options.chrom_peak_snr * self.options.noise_threshold_int)?;
        for spectrum in &input.spectra {
            if spectrum.ms_level != 1 {
                continue;
            }
            finite(spectrum.rt)?;
            input_peaks = add(input_peaks, spectrum.peaks.len())?;
            limit(input_peaks, self.limits.max_peaks, "input peak count")?;
            work.consume(spectrum.peaks.len())?;
            check_array_lengths(spectrum, work)?;
            let mut indices = Vec::new();
            let mut last_mz = None;
            for (index, peak) in spectrum.peaks.iter().enumerate() {
                let intensity = finite(f64::from(peak.intensity))?;
                if intensity > self.options.noise_threshold_int {
                    finite(peak.mz)?;
                    if last_mz.is_some_and(|mz| mz > peak.mz) {
                        return Err(Error::UnsortedData);
                    }
                    last_mz = Some(peak.mz);
                    if intensity > apex_threshold {
                        work.push(
                            &mut apices,
                            Apex {
                                intensity,
                                scan: scans.len(),
                                peak: indices.len(),
                            },
                        )?;
                    }
                    work.push(&mut indices, index)?;
                }
            }
            let offset = total;
            total = add(total, indices.len())?;
            work.push(
                &mut scans,
                Scan {
                    spectrum,
                    indices,
                    offset,
                },
            )?;
        }
        if scans.len() < 3 {
            return Err(bad("at least three MS1 spectra are required"));
        }
        work.sort(apices.len())?;
        // Total order by original encounter index avoids sort scratch and exactly
        // reproduces stable intensity sort followed by reverse iteration.
        apices.sort_unstable_by(|a, b| {
            a.intensity
                .partial_cmp(&b.intensity)
                .unwrap()
                .then(a.scan.cmp(&b.scan))
                .then(a.peak.cmp(&b.peak))
        });
        let arrays = Arrays::discover(&scans, work)?;
        let mut visited = work.vector::<bool>(total)?;
        visited.resize(total, false);
        let end = i64::try_from(total).map_err(|_| bad("progress range overflow"))?;
        self.logger.start_progress(0, end, "mass trace detection")?;
        let result = self.extend(&scans, &apices, arrays, &mut visited, max_traces, work);
        // Balance owned progress nesting even when the checked computation fails.
        let ended = self.logger.end_progress(0);
        let output = result?;
        ended?;
        Ok((output, arrays, warning))
    }

    fn extend(
        &mut self,
        scans: &[Scan<'_>],
        apices: &[Apex],
        arrays: Arrays,
        visited: &mut [bool],
        maximum: usize,
        work: &mut Work,
    ) -> Result<Vec<MassTrace>> {
        let mut output = Vec::new();
        let mut detected = 0usize;
        for apex in apices.iter().rev() {
            work.consume(1)?;
            if visited[scans[apex.scan].offset + apex.peak] {
                continue;
            }
            let peak = scans[apex.scan].peak(apex.peak);
            let mut growing =
                Growing::new(*apex, peak, &scans[apex.scan], arrays, self.options, work)?;
            let mut directions = [Direction::new(apex.scan), Direction::new(apex.scan)];
            loop {
                let mut moved = false;
                // Down precedes up, including centroid/variance and shared hit-rate updates.
                for side in 0..2 {
                    if !directions[side].active {
                        continue;
                    }
                    let next = if side == 0 {
                        directions[side].index.checked_sub(1)
                    } else {
                        directions[side]
                            .index
                            .checked_add(1)
                            .filter(|i| *i < scans.len())
                    };
                    let Some(index) = next else {
                        continue;
                    };
                    moved = true;
                    work.consume(1)?;
                    let scan = &scans[index];
                    if !scan.indices.is_empty() {
                        let candidate = find_candidate(scan, &growing, arrays, self.options, work)?;
                        if let Some(candidate) =
                            candidate.filter(|candidate| !visited[scan.offset + candidate])
                        {
                            growing.accept(index, candidate, scan, arrays, self.options, work)?;
                            directions[side].hits += 1;
                            directions[side].missed = 0;
                        } else {
                            directions[side].missed += 1;
                        }
                    }
                    directions[side].index = index;
                    directions[side].scans += 1;
                    let hits = directions[0].hits + directions[1].hits + 1;
                    let count = directions[0].scans + directions[1].scans + 1;
                    directions[side].active = match self.options.trace_termination_criterion {
                        TraceTerminationCriterion::Outlier => {
                            directions[side].missed <= self.options.trace_termination_outliers
                        }
                        TraceTerminationCriterion::SampleRate => {
                            !(directions[side].scans > 5
                                && (hits as f64 / count as f64) < self.options.min_sample_rate)
                        }
                    };
                }
                if !moved {
                    break;
                }
            }
            let count = directions[0].scans + directions[1].scans + 1;
            let adjusted = count - directions[0].missed - directions[1].missed;
            let range = finite(
                (scans[growing.last_scan].spectrum.rt - scans[growing.first_scan].spectrum.rt)
                    .abs(),
            )?;
            if range < self.options.min_trace_length
                || (self.options.max_trace_length >= 0.0 && range > self.options.max_trace_length)
                || (growing.indices.len() as f64 / adjusted as f64) < self.options.min_sample_rate
            {
                continue;
            }
            limit(
                add(output.len(), 1)?,
                self.limits.max_traces,
                "output trace count",
            )?;
            work.consume(growing.indices.len())?;
            for &(scan, peak) in &growing.indices {
                visited[scans[scan].offset + peak] = true;
            }
            let trace = growing.finish(
                scans,
                arrays,
                self.options,
                output.len() + 1,
                self.limits,
                work,
            )?;
            detected = add(detected, trace.len())?;
            work.push(&mut output, trace)?;
            self.logger.set_progress(
                i64::try_from(detected).map_err(|_| bad("progress count overflow"))?,
            )?;
            if maximum != 0 && output.len() == maximum {
                break;
            }
        }
        Ok(output)
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct Arrays {
    mz: Option<usize>,
    im: Option<usize>,
    im_width: Option<usize>,
}
impl Arrays {
    fn discover(scans: &[Scan<'_>], work: &mut Work) -> Result<Self> {
        let mut result = Self::default();
        work.consume(scans.len())?;
        if let Some(scan) = scans
            .iter()
            .find(|scan| !scan.spectrum.float_data_arrays.is_empty())
        {
            for (index, array) in scan.spectrum.float_data_arrays.iter().enumerate() {
                work.consume(add(1, array.name.len())?)?;
                match array.name.as_str() {
                    user_param::FWHM_MZ_ppm if result.mz.is_none() => result.mz = Some(index),
                    user_param::ION_MOBILITY if result.im.is_none() => result.im = Some(index),
                    user_param::FWHM_IM if result.im_width.is_none() => {
                        result.im_width = Some(index)
                    }
                    _ => (),
                }
            }
        }
        for (index, name) in [
            (result.mz, user_param::FWHM_MZ_ppm),
            (result.im, user_param::ION_MOBILITY),
            (result.im_width, user_param::FWHM_IM),
        ] {
            if let Some(index) = index {
                work.consume(scans.len())?;
                for scan in scans {
                    let array = scan
                        .spectrum
                        .float_data_arrays
                        .get(index)
                        .ok_or_else(|| bad("inconsistent trace metadata arrays"))?;
                    work.consume(array.name.len())?;
                    if array.name != name {
                        return Err(bad("inconsistent trace metadata array positions"));
                    }
                    if array.data.is_empty() && !scan.indices.is_empty() {
                        return Err(bad("trace metadata array length mismatch"));
                    }
                }
            }
        }
        Ok(result)
    }
}
struct Scan<'a> {
    spectrum: &'a MSSpectrum,
    indices: Vec<usize>,
    offset: usize,
}
impl Scan<'_> {
    fn peak(&self, index: usize) -> Peak1D {
        self.spectrum.peaks[self.indices[index]]
    }
    fn value(&self, array: usize, index: usize) -> Result<f64> {
        finite(f64::from(
            self.spectrum.float_data_arrays[array].data[self.indices[index]],
        ))
    }
}
#[derive(Clone, Copy)]
struct Apex {
    intensity: f64,
    scan: usize,
    peak: usize,
}
struct Direction {
    index: usize,
    active: bool,
    scans: usize,
    hits: usize,
    missed: usize,
}
impl Direction {
    fn new(index: usize) -> Self {
        Self {
            index,
            active: true,
            scans: 0,
            hits: 0,
            missed: 0,
        }
    }
}

struct Weighted {
    centroid: f64,
    counter: f64,
    denominator: f64,
}
impl Weighted {
    fn new(value: f64, intensity: f64) -> Result<Self> {
        let mut result = Self {
            centroid: value,
            counter: finite(intensity * value)?,
            denominator: intensity,
        };
        // The source deliberately includes the apex twice in incremental means.
        result.add(value, intensity)?;
        Ok(result)
    }
    fn add(&mut self, value: f64, intensity: f64) -> Result<()> {
        let counter = finite(1.0 + intensity * value / self.counter)?;
        let denominator = finite(1.0 + intensity / self.denominator)?;
        self.centroid = finite(self.centroid * (counter / denominator))?;
        self.counter = finite(self.counter * counter)?;
        self.denominator = finite(self.denominator * denominator)?;
        Ok(())
    }
}
struct Growing {
    indices: Vec<(usize, usize)>,
    mz_widths: Vec<f64>,
    im_widths: Vec<f64>,
    mz: Weighted,
    im: Option<Weighted>,
    sd: f64,
    intensity: f64,
    first_scan: usize,
    last_scan: usize,
}
impl Growing {
    fn new(
        apex: Apex,
        peak: Peak1D,
        scan: &Scan<'_>,
        arrays: Arrays,
        options: MassTraceDetectionOptions,
        work: &mut Work,
    ) -> Result<Self> {
        work.consume(32)?;
        let mz = Weighted::new(peak.mz, apex.intensity)?;
        let im = arrays
            .im
            .map(|array| Weighted::new(scan.value(array, apex.peak)?, apex.intensity))
            .transpose()?;
        let sd = finite((mz.centroid / 1e6) * options.mass_error_ppm)?;
        let mut result = Self {
            indices: Vec::new(),
            mz_widths: Vec::new(),
            im_widths: Vec::new(),
            mz,
            im,
            sd,
            intensity: apex.intensity,
            first_scan: apex.scan,
            last_scan: apex.scan,
        };
        result.record(apex.scan, apex.peak, scan, arrays, work)?;
        Ok(result)
    }
    fn record(
        &mut self,
        scan_index: usize,
        peak: usize,
        scan: &Scan<'_>,
        arrays: Arrays,
        work: &mut Work,
    ) -> Result<()> {
        work.consume(4)?;
        work.push(&mut self.indices, (scan_index, peak))?;
        if let Some(array) = arrays.mz {
            work.push(&mut self.mz_widths, scan.value(array, peak)?)?;
        }
        if let Some(array) = arrays.im_width {
            work.push(&mut self.im_widths, scan.value(array, peak)?)?;
        }
        self.first_scan = self.first_scan.min(scan_index);
        self.last_scan = self.last_scan.max(scan_index);
        Ok(())
    }
    fn accept(
        &mut self,
        index: usize,
        candidate: usize,
        scan: &Scan<'_>,
        arrays: Arrays,
        options: MassTraceDetectionOptions,
        work: &mut Work,
    ) -> Result<()> {
        work.consume(48)?;
        let peak = scan.peak(candidate);
        let intensity = f64::from(peak.intensity);
        self.mz.add(peak.mz, intensity)?;
        if let (Some(im), Some(array)) = (&mut self.im, arrays.im) {
            im.add(scan.value(array, candidate)?, intensity)?;
        }
        self.record(index, candidate, scan, arrays, work)?;
        if options.reestimate_mt_sd {
            let denom1 = self.intensity.ln() + 2.0 * self.sd.ln();
            let denom2 = intensity.ln() + 2.0 * (peak.mz - self.mz.centroid).abs().ln();
            let denom = (denom1.exp() + denom2.exp()).sqrt();
            self.intensity = finite(self.intensity + intensity)?;
            let sd = finite(denom / self.intensity.sqrt())?;
            if sd > f64::EPSILON {
                self.sd = sd;
            }
        }
        Ok(())
    }
    fn finish(
        mut self,
        scans: &[Scan<'_>],
        arrays: Arrays,
        options: MassTraceDetectionOptions,
        number: usize,
        limits: MassTraceDetectionLimits,
        work: &mut Work,
    ) -> Result<MassTrace> {
        work.sort(self.indices.len())?;
        self.indices.sort_unstable_by_key(|&(scan, _)| scan);
        // Conservative shared precharge covers constructor, all three source
        // summary traversals, label formatting/copy, and vector destruction.
        work.consume(add(mul(self.indices.len(), 64)?, 64)?)?;
        let mut peaks = work.vector(self.indices.len())?;
        for (scan, peak) in self.indices {
            let p = scans[scan].peak(peak);
            peaks.push(Peak2D::new(scans[scan].spectrum.rt, p.mz, p.intensity));
        }
        let mut trace = MassTrace::from_peaks_with_limits(
            peaks,
            MassTraceLimits {
                max_peaks: limits.max_peaks,
                max_work: limits.max_work,
                max_bytes: limits.max_bytes,
            },
        )?;
        trace.update_weighted_mean_rt()?;
        trace.update_weighted_mean_mz()?;
        if arrays.mz.is_some() {
            trace.fwhm_mz_avg = median(&mut self.mz_widths, work)?;
        }
        if arrays.im_width.is_some() {
            trace.fwhm_im_avg = median(&mut self.im_widths, work)?;
        }
        if let Some(im) = self.im {
            trace.set_centroid_im(im.centroid)?;
        }
        trace.set_quant_method(options.quant_method);
        trace.update_weighted_mz_sd()?;
        work.allocate(64)?;
        trace.set_label(&format!("T{number}"))?;
        Ok(trace)
    }
}
fn find_candidate(
    scan: &Scan<'_>,
    growing: &Growing,
    arrays: Arrays,
    options: MassTraceDetectionOptions,
    work: &mut Work,
) -> Result<Option<usize>> {
    let lower = finite(growing.mz.centroid - 3.0 * growing.sd)?;
    let upper = finite(growing.mz.centroid + 3.0 * growing.sd)?;
    let candidate = if let (Some(array), Some(im)) = (arrays.im, &growing.im) {
        let lower_im = finite(im.centroid - options.ion_mobility_tolerance)?;
        let upper_im = finite(im.centroid + options.ion_mobility_tolerance)?;
        let first = bound(scan, lower, false, work)?;
        let last = bound(scan, upper, true, work)?;
        let mut best = None;
        for index in first..last {
            work.consume(4)?;
            let value = scan.value(array, index)?;
            if value >= lower_im
                && value <= upper_im
                && best.is_none_or(|old| {
                    (scan.peak(index).mz - growing.mz.centroid).abs()
                        < (scan.peak(old).mz - growing.mz.centroid).abs()
                })
            {
                best = Some(index);
            }
        }
        best
    } else {
        let next = bound(scan, growing.mz.centroid, false, work)?;
        Some(if next == 0 {
            0
        } else if next == scan.indices.len()
            || growing.mz.centroid - scan.peak(next - 1).mz
                <= scan.peak(next).mz - growing.mz.centroid
        {
            next - 1
        } else {
            next
        })
    };
    Ok(candidate.filter(|&index| scan.peak(index).mz >= lower && scan.peak(index).mz <= upper))
}
fn bound(scan: &Scan<'_>, mz: f64, upper: bool, work: &mut Work) -> Result<usize> {
    let (mut first, mut last) = (0, scan.indices.len());
    while first < last {
        work.consume(1)?;
        let middle = first + (last - first) / 2;
        let before = if upper {
            scan.peak(middle).mz <= mz
        } else {
            scan.peak(middle).mz < mz
        };
        if before {
            first = middle + 1;
        } else {
            last = middle;
        }
    }
    Ok(first)
}
fn median(values: &mut [f64], work: &mut Work) -> Result<f64> {
    work.sort(values.len())?;
    values.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    let middle = values.len() / 2;
    if values.len() % 2 == 1 {
        Ok(values[middle])
    } else {
        finite((values[middle - 1] + values[middle]) / 2.0)
    }
}
fn check_array_lengths(spectrum: &MSSpectrum, work: &mut Work) -> Result<()> {
    let count = add(
        add(
            spectrum.float_data_arrays.len(),
            spectrum.integer_data_arrays.len(),
        )?,
        spectrum.string_data_arrays.len(),
    )?;
    work.consume(count)?;
    for len in spectrum
        .float_data_arrays
        .iter()
        .map(|a| a.data.len())
        .chain(spectrum.integer_data_arrays.iter().map(|a| a.data.len()))
        .chain(spectrum.string_data_arrays.iter().map(|a| a.data.len()))
    {
        if len != 0 && len != spectrum.peaks.len() {
            return Err(bad("spectrum data array length mismatch"));
        }
    }
    Ok(())
}
fn ccs_warning(input: &MSExperiment, tolerance: f64, work: &mut Work) -> Result<bool> {
    work.consume(input.spectra.len())?;
    for spectrum in &input.spectra {
        work.consume(spectrum.float_data_arrays.len())?;
        for array in &spectrum.float_data_arrays {
            work.consume(mul(add(array.name.len(), 1)?, 32)?)?;
            if let Some(ccs) = im_unit_is_ccs(&array.name) {
                return Ok(ccs && tolerance < 1.0);
            }
        }
    }
    Ok(false)
}
fn im_unit_is_ccs(name: &str) -> Option<bool> {
    // All nine proper descendants in pinned PSI-MS have millisecond or VSSC
    // units. The generic parent term is not a child of itself; none has CCS.
    match name {
        "mean ion mobility array"
        | "mean inverse reduced ion mobility array"
        | "raw ion mobility array"
        | "raw inverse reduced ion mobility array"
        | "deconvoluted ion mobility array"
        | "deconvoluted inverse reduced ion mobility array" => return Some(false),
        "mean ion mobility drift time array" => return Some(false),
        "raw ion mobility drift time array" => return Some(false),
        "deconvoluted ion mobility drift time array" => return Some(false),
        _ => (),
    }
    if name.starts_with(user_param::MEAN_INVERSE_REDUCED_ION_MOBILITY_ARRAY)
        || name.starts_with(user_param::INVERSE_REDUCED_ION_MOBILITY)
    {
        return Some(false);
    }
    if name.starts_with(user_param::ION_MOBILITY) {
        return Some(
            !(name.contains("MS:1002815") || name.contains("MS:1003006"))
                && name.contains("MS:1002954"),
        );
    }
    None
}

struct Work {
    remaining: usize,
    bytes: usize,
}
impl Work {
    fn new(limits: MassTraceDetectionLimits) -> Self {
        Self {
            remaining: limits.max_work,
            bytes: limits.max_bytes,
        }
    }
    fn consume(&mut self, count: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(count)
            .ok_or_else(|| bad("mass trace detection work limit exceeded"))?;
        Ok(())
    }
    fn allocate(&mut self, bytes: usize) -> Result<()> {
        self.bytes = self
            .bytes
            .checked_sub(bytes)
            .ok_or_else(|| bad("mass trace detection allocation limit exceeded"))?;
        Ok(())
    }
    fn vector<T>(&mut self, len: usize) -> Result<Vec<T>> {
        self.consume(len)?;
        self.allocate(mul(len, size_of::<T>())?)?;
        let mut result = Vec::new();
        result
            .try_reserve_exact(len)
            .map_err(|_| bad("mass trace detection allocation failed"))?;
        Ok(result)
    }
    fn push<T>(&mut self, values: &mut Vec<T>, value: T) -> Result<()> {
        self.consume(1)?;
        if values.len() == values.capacity() {
            let next = values.len().saturating_mul(2).max(1);
            self.consume(values.len())?;
            self.allocate(mul(next, size_of::<T>())?)?;
            values
                .try_reserve_exact(next - values.len())
                .map_err(|_| bad("mass trace detection allocation failed"))?;
        }
        values.push(value);
        Ok(())
    }
    fn sort(&mut self, len: usize) -> Result<()> {
        let levels = usize::BITS as usize - len.leading_zeros() as usize;
        self.consume(mul(len, mul(levels + 1, 4)?)?)
    }
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b)
        .ok_or_else(|| bad("mass trace detection size overflow"))
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b)
        .ok_or_else(|| bad("mass trace detection size overflow"))
}
fn limit(value: usize, maximum: usize, what: &str) -> Result<()> {
    if value > maximum {
        Err(bad(&format!("{what} limit exceeded")))
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad("nonfinite mass trace detection value"))
    }
}
fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn vector_initialization_is_precharged_before_allocation() {
        let mut work = Work {
            remaining: 2,
            bytes: 10,
        };
        assert!(work.vector::<bool>(3).is_err());
        assert_eq!(work.bytes, 10);
        let mut work = Work {
            remaining: 3,
            bytes: 3,
        };
        let mut values = work.vector::<bool>(3).unwrap();
        values.resize(3, false);
        assert_eq!((work.remaining, work.bytes), (0, 0));
    }
    #[test]
    fn source_iterative_mean_literals_and_shared_trace_finalization_budget() {
        let (value, intensity) = (150.22, 25_000_000.0);
        let mut mean = Weighted {
            centroid: value,
            counter: value * intensity,
            denominator: intensity,
        };
        mean.add(150.34, 23_043_030.).unwrap();
        let expected = (value * intensity + 150.34 * 23_043_030.) / (intensity + 23_043_030.);
        assert!((mean.centroid - expected).abs() < 1e-12);
        mean.add(150.11, 1_932_392.).unwrap();
        let expected = (value * intensity + 150.34 * 23_043_030. + 150.11 * 1_932_392.)
            / (intensity + 23_043_030. + 1_932_392.);
        assert!((mean.centroid - expected).abs() < 1e-12);
        let mut d = MassTraceDetection::new();
        d.logger.set_log_type(ProgressLogType::None);
        d.options.min_trace_length = 0.;
        let input = MSExperiment {
            spectra: (0..3)
                .map(|i| MSSpectrum {
                    rt: i as f64,
                    peaks: vec![Peak1D::new(100., 100.), Peak1D::new(200., 100.)],
                    ..Default::default()
                })
                .collect(),
            ..Default::default()
        };
        let mut work = Work::new(d.limits);
        d.execute(&input, 1, &mut work).unwrap();
        let one_cost = d.limits.max_work - work.remaining;
        d.limits.max_work = one_cost;
        assert!(d.run(&input, 1).is_ok());
        assert!(d.run(&input, 0).is_err());
    }
}
