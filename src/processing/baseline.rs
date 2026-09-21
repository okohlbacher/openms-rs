// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Baseline filtering with methods from mathematical morphology.
//!
//! Ports the header-only `OpenMS/PROCESSING/BASELINE/MorphologicalFilter.h`
//! (its `.cpp` is empty) at core `bc9cc12`. The API mapping, the source
//! conventions this module keeps and the executed C++ evidence are in
//! `docs/MORPHOLOGICAL_FILTER_SUPPORT.md`.
//!
//! The fundamental operations are erosion and dilation with respect to a
//! structuring element, which here is a straight line of `struc_size` samples.
//! For an input `x_0, x_1, ...` the erosion at `i` is the minimum of the
//! window `x_{i - struc_size/2} ..= x_{i + struc_size/2}` and the dilation is
//! its maximum, with the window clipped at both ends of the signal. Baseline
//! filtering uses the top-hat transform: the signal minus its opening, where
//! the opening is the dilation of the erosion.
//!
//! The source's `@note`s carry over: the filter is designed for uniformly
//! spaced profile data, and the data must be sorted by ascending m/z. Sorting
//! is checked here (the source does not check it).
//!
//! # The source's single-sample element
//!
//! The source's van Herk erosion and dilation never write their last output
//! sample when the element is one sample long and the signal has more than
//! five samples. That sample keeps the value the output buffer already held.
//! This module reproduces that result, because the `BaselineFilter` tool
//! reaches it whenever an element in Thomson is narrower than the peak
//! spacing, which is the common case for centroided MS2 spectra: there the
//! source top-hat keeps the last peak and zeroes every other one. Every other
//! element length and signal length gives the clipped-window result above;
//! the executed sweep behind this statement is described in the support
//! document.

use super::{AcquisitionCopies, SpectrumFilter};
use crate::concept::progress_logger::{ProgressLogger, ProgressReporter, progress_value};
use crate::kernel::{MSChromatogram, MSExperiment, MSSpectrum, SpectrumType};
use crate::{Error, Result};
use std::collections::VecDeque;

/// The progress label of source `filterExperiment`
/// (`MorphologicalFilter.h:306`).
pub const BASELINE_PROGRESS_LABEL: &str = "filtering baseline";

/// Morphological operation, the source `method` parameter.
///
/// The source validates `method` against a string list and silently does
/// nothing for a string outside it; an enum makes that case unrepresentable,
/// so the source's `Exception::IllegalArgument` note for `filterRange` has no
/// Rust counterpart.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MorphologicalMethod {
    /// `identity`: copy the input unchanged.
    Identity,
    /// `erosion`: sliding-window minimum, with the source's single-sample
    /// element behaviour described in the module documentation.
    Erosion,
    /// `dilation`: sliding-window maximum, with the source's single-sample
    /// element behaviour described in the module documentation.
    Dilation,
    /// `opening`: the dilation of the erosion.
    Opening,
    /// `closing`: the erosion of the dilation.
    Closing,
    /// `gradient`: the dilation minus the erosion.
    ///
    /// With a single-sample element on more than five samples the source
    /// subtracts, at the last sample, a value its erosion left in a buffer
    /// shared by every call. See [`MorphologicalFilter`] for how far that
    /// history is reproduced.
    Gradient,
    /// `tophat`: the input minus the opening. The source default and the
    /// operation used for baseline removal.
    #[default]
    TopHat,
    /// `bothat`: OpenMS's signed bottom-hat, the input minus the closing, so
    /// the result is generally non-positive.
    BottomHat,
    /// `erosion_simple`: the source's direct-window erosion, which gives the
    /// clipped-window minimum at every sample for every element length.
    ErosionSimple,
    /// `dilation_simple`: the source's direct-window dilation, which gives the
    /// clipped-window maximum at every sample for every element length.
    DilationSimple,
}

/// Length of the structuring element, the source `struc_elem_length` and
/// `struc_elem_unit` parameters.
///
/// The source default is `Thomson(3.0)`; its parameter documentation asks for
/// an element wider than the expected peak width. Whatever the unit, the
/// spectrum and chromatogram methods round the resulting sample count up to an
/// odd number, as source `filter(MSSpectrum&)` does.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum StructuringElement {
    /// `DataPoints`: a sample count. Zero is rejected.
    ///
    /// The source truncates its double parameter to an unsigned count and
    /// `filter(MSSpectrum&)` turns a count of zero into one; this port refuses
    /// zero instead of choosing the width for the caller.
    DataPoints(usize),
    /// `Thomson`: a width in m/z, converted to samples with the spectrum's
    /// global average spacing (see [`MorphologicalFilter::effective_window`]).
    ///
    /// Must be finite and positive. The source accepts zero, which its
    /// conversion turns into one sample, and a negative width, whose
    /// conversion to an unsigned count is undefined behaviour in C++.
    Thomson(f64),
    /// Chromatogram extension with no source counterpart: a width in
    /// seconds, converted with the chromatogram's global average spacing.
    Seconds(f64),
}

/// Morphological filter over spectra, chromatograms and plain intensity
/// ranges (source `MorphologicalFilter`).
///
/// The source class is a `DefaultParamHandler` and `ProgressLogger`; here the
/// three parameters are typed public fields, and the progress of
/// `filterExperiment` goes to a caller's logger through
/// [`MorphologicalFilter::filter_experiment_with_progress`].
/// The source declares, but never defines, a copy constructor, so it cannot
/// be copied; this filter has no state between calls and is `Copy`.
///
/// # Source state between calls
///
/// Source `filterRange` keeps a function-local `static` buffer that is shared
/// by every call in a process and is only ever grown. The gradient's last
/// sample under the single-sample element reads a stale value from it. A
/// fresh TOPP process starts with an empty buffer, so this port models the
/// buffer as follows:
///
/// - [`SpectrumFilter::filter_experiment`] carries one buffer through the
///   spectra in order, as source `filterExperiment` does inside the
///   `BaselineFilter` tool. Spectra with fewer than two peaks never reach the
///   source's `filterRange` and leave the buffer untouched.
/// - [`SpectrumFilter::filter_spectrum`], [`MorphologicalFilter::filter_range`]
///   and [`MorphologicalFilter::filter_chromatogram`] start from an empty
///   buffer, as the first source call in a process does.
///
/// A C++ program that calls the source filter several times sees the history
/// of all earlier calls instead; that history is not reproduced across
/// separate calls here. Only [`MorphologicalMethod::Gradient`] can observe it.
#[derive(Clone, Copy, Debug)]
pub struct MorphologicalFilter {
    /// Operation to apply (source `method`, default `tophat`).
    pub method: MorphologicalMethod,
    /// Element length and unit (source `struc_elem_length` and
    /// `struc_elem_unit`, default 3 Thomson).
    pub structuring_element: StructuringElement,
}

impl Default for MorphologicalFilter {
    /// The source defaults: top-hat with a 3 Thomson element.
    fn default() -> Self {
        Self {
            method: MorphologicalMethod::TopHat,
            structuring_element: StructuringElement::Thomson(3.0),
        }
    }
}

/// The source erosion and dilation fall back to the direct-window method when
/// the signal has at most this many samples (`size <= 5`).
const SOURCE_SIMPLE_MAX_SAMPLES: usize = 5;

/// The source's function-local `static` buffer in `filterRange`: grown to the
/// signal length, never shrunk or cleared.
#[derive(Debug, Default)]
struct RangeBuffer {
    values: Vec<f32>,
}

impl RangeBuffer {
    /// The first `len` samples, growing with zeros as `std::vector::resize` does.
    fn prefix(&mut self, len: usize) -> &mut [f32] {
        if self.values.len() < len {
            self.values.resize(len, 0.0);
        }
        &mut self.values[..len]
    }
}

impl MorphologicalFilter {
    /// [`SpectrumFilter::filter_experiment`], reporting progress to `progress`
    /// as the source's `ProgressLogger` base does.
    ///
    /// Source `filterExperiment` calls `startProgress(0, exp.size(),
    /// "filtering baseline")` (`MorphologicalFilter.h:306`), `setProgress(i)`
    /// after filtering spectrum `i` (`:310`), so the values run from `0` to
    /// `n - 1`, and `endProgress()` after the last (`:312`). This makes the
    /// same calls; the filtered experiment is the one
    /// [`SpectrumFilter::filter_experiment`] produces. The metadata-copy
    /// preflight happens before the section starts, so an experiment refused
    /// there prints nothing.
    ///
    /// # Errors
    ///
    /// As [`SpectrumFilter::filter_experiment`], and the errors of `progress`.
    /// An error inside the section still ends it, which the source does not
    /// do; see [`ProgressReporter::section`]. The experiment is unchanged on
    /// error.
    pub fn filter_experiment_with_progress(
        &self,
        experiment: &mut MSExperiment,
        progress: &mut ProgressLogger,
    ) -> Result<()> {
        self.filter_experiment_reporting(experiment, &mut ProgressReporter::new(Some(progress)))
    }

    /// Source `filterExperiment` into a copy, inside one progress section.
    fn filter_experiment_reporting(
        &self,
        experiment: &mut MSExperiment,
        reporter: &mut ProgressReporter<'_>,
    ) -> Result<()> {
        for spectrum in &experiment.spectra {
            AcquisitionCopies::default().spectrum(spectrum)?;
        }
        let records = progress_value(experiment.spectra.len())?;
        let mut spectra = experiment.spectra.clone();
        let mut buffer = RangeBuffer::default();
        reporter.section(0, records, BASELINE_PROGRESS_LABEL, |reporter| {
            for (index, spectrum) in spectra.iter_mut().enumerate() {
                self.filter_spectrum_with(spectrum, &mut buffer)?;
                // :310, `setProgress(i)`.
                reporter.set_count(index)?;
            }
            Ok(())
        })?;
        experiment.spectra = spectra;
        Ok(())
    }

    /// A filter with checked options.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for `DataPoints(0)` and for a
    /// `Thomson` or `Seconds` width that is not finite and positive.
    pub fn new(
        method: MorphologicalMethod,
        structuring_element: StructuringElement,
    ) -> Result<Self> {
        let filter = Self {
            method,
            structuring_element,
        };
        filter.validate_options()?;
        Ok(filter)
    }

    fn validate_options(&self) -> Result<()> {
        match self.structuring_element {
            StructuringElement::DataPoints(0) => Err(Error::InvalidValue(
                "structuring element must contain at least one point".into(),
            )),
            StructuringElement::Thomson(width) | StructuringElement::Seconds(width)
                if !width.is_finite() || width <= 0.0 =>
            {
                Err(Error::InvalidValue(
                    "structuring element width must be finite and positive".into(),
                ))
            }
            _ => Ok(()),
        }
    }

    /// The odd element length in samples for a record with these positions.
    ///
    /// A `DataPoints` count is used as given. A `Thomson` or `Seconds` width is
    /// converted as source `filter(MSSpectrum&)` does: the samples are assumed
    /// to be uniformly spaced, the average spacing comes from the first and
    /// last position and the sample count, and the width becomes
    /// `ceil(width * (N - 1) / (last - first))` samples. Either way an even
    /// count is then rounded up to the next odd one. A record with fewer than
    /// two positions yields one sample for a width.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for invalid options, non-finite
    /// positions, a width conversion over a zero or non-finite span (where the
    /// source divides by zero and converts the result to an unsigned count,
    /// which is undefined behaviour), or a count that overflows `usize`;
    /// [`Error::UnsortedData`] when the positions are not sorted ascending.
    pub fn effective_window(&self, positions: &[f64]) -> Result<usize> {
        self.validate_options()?;
        validate_positions(positions)?;
        let points = match self.structuring_element {
            StructuringElement::DataPoints(points) => points,
            StructuringElement::Thomson(width) | StructuringElement::Seconds(width) => {
                if positions.len() < 2 {
                    return Ok(1);
                }
                let span = positions[positions.len() - 1] - positions[0];
                if !span.is_finite() || span <= 0.0 {
                    return Err(Error::InvalidValue(
                        "coordinate-based morphology needs a positive finite signal span".into(),
                    ));
                }
                let points = (width * (positions.len() - 1) as f64 / span).ceil();
                if !points.is_finite() || points >= usize::MAX as f64 {
                    return Err(Error::InvalidValue(
                        "structuring element length overflows".into(),
                    ));
                }
                points as usize
            }
        };
        if points % 2 == 0 {
            points
                .checked_add(1)
                .ok_or_else(|| Error::InvalidValue("structuring element length overflows".into()))
        } else {
            Ok(points)
        }
    }

    /// Apply the operation to a plain intensity range (source `filterRange`).
    ///
    /// The source writes into a caller-allocated output range; this returns a
    /// new vector whose samples start at zero, which is what a sample the
    /// source leaves unwritten holds in the source class test. The element
    /// must be given in `DataPoints`; source `filterRange` reads the
    /// `struc_elem_length` number whatever the unit says.
    ///
    /// An even count is rounded up to odd, as the spectrum method does. Source
    /// `filterRange` uses an even count as given, which its van Herk erosion
    /// and dilation turn into asymmetric windows; this port does not reproduce
    /// those. A singleton range is filtered like any other, so a singleton
    /// top-hat range is zero while [`SpectrumFilter::filter_spectrum`] leaves a
    /// singleton spectrum unchanged, as the source wrapper does.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for invalid options, a unit other than
    /// `DataPoints`, a non-finite intensity, or a subtraction that overflows
    /// `f32` (the source would store an infinity).
    pub fn filter_range(&self, intensities: &[f32]) -> Result<Vec<f32>> {
        self.validate_options()?;
        let StructuringElement::DataPoints(points) = self.structuring_element else {
            return Err(Error::InvalidValue(
                "intensity-only morphology requires DataPoints units".into(),
            ));
        };
        let points = if points % 2 == 0 {
            points.saturating_add(1)
        } else {
            points
        };
        self.apply(intensities, points, &mut RangeBuffer::default())
    }

    /// Source `filterRange` for an odd element of `points` samples, reading
    /// and updating the source's shared buffer.
    fn apply(&self, input: &[f32], points: usize, buffer: &mut RangeBuffer) -> Result<Vec<f32>> {
        if input.iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidValue(
                "morphology intensities must be finite".into(),
            ));
        }
        use MorphologicalMethod::*;
        let len = input.len();
        let mut output = vec![0.0; len];
        match self.method {
            Identity => output.copy_from_slice(input),
            Erosion => source_extrema(input, points, false, &mut output),
            Dilation => source_extrema(input, points, true, &mut output),
            ErosionSimple => clipped_extrema(input, points, false, &mut output, len),
            DilationSimple => clipped_extrema(input, points, true, &mut output, len),
            Opening => {
                let scratch = buffer.prefix(len);
                source_extrema(input, points, false, scratch);
                source_extrema(scratch, points, true, &mut output);
            }
            Closing => {
                let scratch = buffer.prefix(len);
                source_extrema(input, points, true, scratch);
                source_extrema(scratch, points, false, &mut output);
            }
            Gradient => {
                let scratch = buffer.prefix(len);
                source_extrema(input, points, false, scratch);
                source_extrema(input, points, true, &mut output);
                for (out, &eroded) in output.iter_mut().zip(scratch.iter()) {
                    *out = checked_difference(*out, eroded)?;
                }
            }
            TopHat | BottomHat => {
                let scratch = buffer.prefix(len);
                let first_is_max = self.method == BottomHat;
                source_extrema(input, points, first_is_max, scratch);
                source_extrema(scratch, points, !first_is_max, &mut output);
                for (out, &value) in output.iter_mut().zip(input) {
                    *out = checked_difference(value, *out)?;
                }
            }
        }
        Ok(output)
    }

    /// Apply the operation to a chromatogram's intensities.
    ///
    /// A convenience extension: the source class filters spectra and
    /// experiments only. The element is given in `Seconds` or `DataPoints`
    /// and converted like the spectrum element; records with fewer than two
    /// peaks are left unchanged. The source's shared buffer starts empty (see
    /// [`MorphologicalFilter`]).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for invalid options, a `Thomson`
    /// element, a failed element conversion, non-finite intensities or an
    /// overflowing subtraction; [`Error::UnsortedData`] for unsorted retention
    /// times; and the chromatogram's own validation errors. The chromatogram
    /// is unchanged on error.
    pub fn filter_chromatogram(&self, chromatogram: &mut MSChromatogram) -> Result<()> {
        self.validate_options()?;
        chromatogram.validate()?;
        if matches!(self.structuring_element, StructuringElement::Thomson(_)) {
            return Err(Error::InvalidValue(
                "chromatogram morphology uses Seconds or DataPoints".into(),
            ));
        }
        let positions: Vec<_> = chromatogram.peaks.iter().map(|p| p.rt).collect();
        let points = self.effective_window(&positions)?;
        if chromatogram.len() <= 1 {
            return Ok(());
        }
        let values: Vec<_> = chromatogram.peaks.iter().map(|p| p.intensity).collect();
        let output = self.apply(&values, points, &mut RangeBuffer::default())?;
        for (peak, value) in chromatogram.peaks.iter_mut().zip(output) {
            peak.intensity = value;
        }
        Ok(())
    }

    /// Source `filter(MSSpectrum&)` with the source's shared buffer passed in.
    fn filter_spectrum_with(
        &self,
        spectrum: &mut MSSpectrum,
        buffer: &mut RangeBuffer,
    ) -> Result<()> {
        self.validate_options()?;
        spectrum.validate()?;
        if matches!(self.structuring_element, StructuringElement::Seconds(_)) {
            return Err(Error::InvalidValue(
                "spectrum morphology uses Thomson or DataPoints".into(),
            ));
        }
        let positions: Vec<_> = spectrum.peaks.iter().map(|p| p.mz).collect();
        let points = self.effective_window(&positions)?;
        if spectrum.len() > 1 {
            let values: Vec<_> = spectrum.peaks.iter().map(|p| p.intensity).collect();
            let output = self.apply(&values, points, buffer)?;
            for (peak, value) in spectrum.peaks.iter_mut().zip(output) {
                peak.intensity = value;
            }
        }
        spectrum.spectrum_type = SpectrumType::Profile;
        Ok(())
    }
}

impl SpectrumFilter for MorphologicalFilter {
    /// Filter one spectrum's intensities (source `filter(MSSpectrum&)`).
    ///
    /// The element length comes from [`MorphologicalFilter::effective_window`]
    /// for this spectrum. As in the source, the spectrum is marked `Profile`
    /// and a spectrum with fewer than two peaks keeps its intensities. The
    /// source's shared buffer starts empty (see [`MorphologicalFilter`]).
    /// m/z values, metadata and data arrays are untouched.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for invalid options, a `Seconds`
    /// element, a failed element conversion, non-finite intensities or an
    /// overflowing subtraction; [`Error::UnsortedData`] when the peaks are not
    /// sorted by m/z, which the source documents as a precondition without
    /// checking it; and the spectrum's own validation errors. The spectrum is
    /// unchanged on error.
    fn filter_spectrum(&self, spectrum: &mut MSSpectrum) -> Result<()> {
        self.filter_spectrum_with(spectrum, &mut RangeBuffer::default())
    }

    /// Filter every spectrum in order (source `filterExperiment`).
    ///
    /// The element length is computed for each spectrum separately, and one
    /// source buffer is carried through the spectra in order (see
    /// [`MorphologicalFilter`]), so the work stays serial: spectrum order is
    /// observable through [`MorphologicalMethod::Gradient`]. Chromatograms and
    /// experiment metadata are unchanged, as in the source.
    ///
    /// The spectra are copied before they are filtered, so an error leaves the
    /// experiment untouched, and the copy is metered **per spectrum** rather
    /// than once for the whole experiment, as the trait's default does: one
    /// ledger over a whole run rejects real data, because the ceiling would
    /// then shrink as the run grows. The benchmark's 40,856-spectrum UK222
    /// file exhausts the shared ledger, while every one of its spectra is far
    /// inside the per-record allowance. What is copied is metadata this
    /// experiment already holds in memory, so the record is the level the
    /// bound belongs at.
    ///
    /// # Errors
    ///
    /// Returns the first spectrum's error from
    /// [`SpectrumFilter::filter_spectrum`], or [`Error::InvalidValue`] when one
    /// spectrum's metadata exceeds the processing copy budget. The experiment
    /// is unchanged on error.
    ///
    /// Reports no progress; see
    /// [`MorphologicalFilter::filter_experiment_with_progress`].
    fn filter_experiment(&self, experiment: &mut MSExperiment) -> Result<()> {
        self.filter_experiment_reporting(experiment, &mut ProgressReporter::silent())
    }
}

fn validate_positions(positions: &[f64]) -> Result<()> {
    if positions.iter().any(|x| !x.is_finite()) {
        return Err(Error::InvalidValue("coordinates must be finite".into()));
    }
    if positions.windows(2).any(|p| p[0] > p[1]) {
        return Err(Error::UnsortedData);
    }
    Ok(())
}

/// `a - b` in `f32`, as the source computes it, refusing an overflow.
fn checked_difference(a: f32, b: f32) -> Result<f32> {
    let difference = a - b;
    if difference.is_finite() {
        Ok(difference)
    } else {
        Err(Error::InvalidValue("intensity overflow".into()))
    }
}

/// Source `applyErosion_` / `applyDilation_` into `target`.
///
/// The source takes the direct-window method for `size <= struc_size` or
/// `size <= 5` and van Herk's block method otherwise. Both give clipped-window
/// extrema, except that van Herk's method with a one-sample element writes
/// every output sample but the last. That sample keeps whatever `target`
/// held. The windows themselves are computed with a monotonic deque rather
/// than van Herk's prefix and suffix blocks; both do linear work.
///
/// The source passes the element length as a signed `Int`, so a length of
/// 2^31 or more would take the van Herk branch with a negative length, which
/// is undefined behaviour; `usize` lengths here always take the defined path.
fn source_extrema(values: &[f32], points: usize, maximum: bool, target: &mut [f32]) {
    let len = values.len();
    let written = if points == 1 && len > SOURCE_SIMPLE_MAX_SAMPLES {
        len - 1
    } else {
        len
    };
    clipped_extrema(values, points, maximum, target, written);
}

/// Clipped-window minima or maxima for the first `written` centres into
/// `target`, in linear time with a monotonic deque of indices.
fn clipped_extrema(
    values: &[f32],
    points: usize,
    maximum: bool,
    target: &mut [f32],
    written: usize,
) {
    let radius = points / 2;
    let mut queue: VecDeque<usize> = VecDeque::new();
    let mut next = 0;
    for (center, slot) in target
        .iter_mut()
        .enumerate()
        .take(written.min(values.len()))
    {
        let start = center.saturating_sub(radius);
        let end = center
            .saturating_add(radius)
            .saturating_add(1)
            .min(values.len());
        while next < end {
            while queue.back().is_some_and(|&i| {
                if maximum {
                    values[i] < values[next]
                } else {
                    values[i] > values[next]
                }
            }) {
                queue.pop_back();
            }
            queue.push_back(next);
            next += 1;
        }
        while queue.front().is_some_and(|&i| i < start) {
            queue.pop_front();
        }
        // The window always contains its centre, so the queue is never empty.
        *slot = values[queue.front().copied().unwrap_or(center)];
    }
}
