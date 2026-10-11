// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Extraction of chromatograms (XICs) from spectral data.
//!
//! Ports `ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h` with
//! `ChromatogramExtractorAlgorithm.cpp`. The algorithm takes a set of
//! extraction coordinates and extracts the signal of a spectral map at the
//! product-ion m/z (and optionally ion mobility) within each coordinate's
//! retention-time range. The coordinates are up to the caller;
//! [`prepare_coordinates`](crate::analysis::openswath::chromatogram_extractor::prepare_coordinates)
//! builds them for the common MS1 and MS2 cases. For MS2 extraction the map is
//! assumed to come from a SWATH (data-independent acquisition) experiment. See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_EXTRACTION_SUPPORT.md).
//!
//! # Source walk, preserved
//!
//! The top-hat integrator walks outwards from the first peak at or above the
//! target m/z and sums every peak strictly inside the window. Two properties of
//! the source walk are reproduced exactly, because the upstream class test
//! encodes the first of them in its expected values:
//!
//! - When the target lies at least two peaks to the right of the first peak of
//!   the spectrum, the leftward walk stops *on* the first peak without adding
//!   it, even when it lies inside the window
//!   (`ChromatogramExtractorAlgorithm.cpp:77-103`). The class test's
//!   `400.05 -> 8400` omits the intensity 8 of the first peak for this reason.
//! - When the target lies above every peak, the last peak is added twice: once
//!   as the "current" peak (`:63-75`) and again as the first step of the
//!   leftward walk (`:77-103`).
//!
//! Both are recorded as C++ defect candidates in the support document.
//!
//! A third property decides the result of a mixed extraction. The source
//! advances its ion-mobility iterator only in the ion-mobility overload of the
//! integrator, so after a coordinate without ion mobility (a negative
//! `ion_mobility`) the iterator lags behind the m/z iterator for the rest of
//! the spectrum, and every later ion-mobility coordinate filters on the ion
//! mobility of an earlier peak.
//! [`extract_chromatograms`](crate::analysis::openswath::chromatogram_extractor_algorithm::ChromatogramExtractorAlgorithm::extract_chromatograms)
//! reproduces the lag (see
//! [`extract_value_tophat_im_lagged`](crate::analysis::openswath::chromatogram_extractor_algorithm::ChromatogramExtractorAlgorithm::extract_value_tophat_im_lagged));
//! an extraction in which every coordinate carries ion mobility, or none does,
//! is not affected.
//!
//! # Ceilings
//!
//! One extraction appends one point per matching spectrum and coordinate. The
//! product of the spectrum count and the coordinate count is checked against
//! [`MAX_EXTRACTED_POINTS`](crate::analysis::openswath::chromatogram_extractor_algorithm::MAX_EXTRACTED_POINTS) before anything is allocated.

use super::data_structures::Chromatogram;
use super::spectrum_access::SpectrumAccess;
use crate::concept::progress_logger::{ProgressLogger, ProgressReporter, progress_value};
use crate::{Error, Result};
use std::cmp::Ordering;

/// Largest number of (spectrum, coordinate) pairs one extraction may visit.
///
/// Native ceiling; the source has none. `2^30` pairs are 16 GiB of output
/// arrays at worst, well above a full SWATH run.
pub const MAX_EXTRACTED_POINTS: usize = 1 << 30;

/// The progress label of [`ChromatogramExtractorAlgorithm::extract_chromatograms`],
/// source `ChromatogramExtractorAlgorithm.cpp:275`.
pub const EXTRACTION_PROGRESS_LABEL: &str = "Extracting chromatograms";

/// Where and what to extract: m/z, ion mobility and an RT range, with an id.
///
/// Source `ChromatogramExtractorAlgorithm::ExtractionCoordinates`. All numeric
/// fields default to 0, as in the source.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ExtractionCoordinates {
    /// m/z value around which should be extracted.
    pub mz: f64,
    /// Ion mobility value around which should be extracted; a negative value
    /// disables ion-mobility filtering for this coordinate.
    pub ion_mobility: f64,
    /// Precursor m/z value; the algorithm currently ignores it, as the source
    /// does.
    pub mz_precursor: f64,
    /// RT start of the extraction, in seconds.
    pub rt_start: f64,
    /// RT end of the extraction, in seconds. When `rt_end - rt_start` is not
    /// positive the whole chromatogram is extracted.
    pub rt_end: f64,
    /// Identifier (usually the native id of the output chromatogram).
    pub id: String,
}

impl ExtractionCoordinates {
    /// Whether `left` sorts before `right` by m/z (`left.mz < right.mz`).
    ///
    /// Source `SortExtractionCoordinatesByMZ`, a strict-weak-order predicate;
    /// NaN compares false either way, as in the source.
    pub fn sort_by_mz(left: &Self, right: &Self) -> bool {
        left.mz < right.mz
    }

    /// Whether `left` sorts after `right` by m/z (`left.mz > right.mz`).
    ///
    /// Source `SortExtractionCoordinatesReverseByMZ`.
    pub fn sort_reverse_by_mz(left: &Self, right: &Self) -> bool {
        left.mz > right.mz
    }

    /// An [`Ordering`] by m/z for `sort_by`, equal to
    /// [`ExtractionCoordinates::sort_by_mz`] on non-NaN values.
    ///
    /// Native helper: a NaN m/z compares equal to everything, which is what
    /// the source predicate yields for NaN under `std::stable_sort`.
    pub fn cmp_by_mz(left: &Self, right: &Self) -> Ordering {
        left.mz.partial_cmp(&right.mz).unwrap_or(Ordering::Equal)
    }
}

/// The m/z filter profile, source `getFilterNr_` (1 = top-hat, 2 = Bartlett).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExtractionFilter {
    /// Sum every peak strictly inside the window.
    Tophat,
    /// A triangular (Bartlett) weighting; declared by the source, never
    /// implemented there.
    Bartlett,
}

impl ExtractionFilter {
    /// Parse a filter name, source `getFilterNr_`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for anything but `"tophat"` or `"bartlett"`,
    /// the source's `Exception::IllegalArgument`.
    pub fn from_name(filter: &str) -> Result<Self> {
        match filter {
            "tophat" => Ok(Self::Tophat),
            "bartlett" => Ok(Self::Bartlett),
            _ => Err(Error::InvalidValue(
                "Filter either needs to be tophat or bartlett".into(),
            )),
        }
    }
}

/// Chromatogram extraction at given coordinates.
///
/// Source `ChromatogramExtractorAlgorithm`, a `ProgressLogger` with no state of
/// its own; the progress logger is passed to
/// [`ChromatogramExtractorAlgorithm::extract_chromatograms_with_progress`]
/// instead of being inherited.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ChromatogramExtractorAlgorithm;

/// The `[left, right]` m/z window around `mz`, exclusive at both ends in use.
fn mz_window(mz: f64, mz_extraction_window: f64, ppm: bool) -> (f64, f64) {
    if ppm {
        (
            mz - mz * mz_extraction_window / 2.0 * 1.0e-6,
            mz + mz * mz_extraction_window / 2.0 * 1.0e-6,
        )
    } else {
        (
            mz - mz_extraction_window / 2.0,
            mz + mz_extraction_window / 2.0,
        )
    }
}

/// The shared walk of both source `extract_value_tophat` overloads.
///
/// `accept(i)` is the ion-mobility test of peak `i` (always true for the m/z
/// only overload). The walk bounds depend on m/z alone, as in the source.
fn tophat_walk(
    mz: &[f64],
    intensity: &[f64],
    pos: &mut usize,
    target: f64,
    left: f64,
    right: f64,
    accept: impl Fn(usize) -> bool,
) -> f64 {
    let n = mz.len();
    let mut integrated = 0.0;
    if n == 0 {
        return integrated;
    }
    let inside = |i: usize| mz[i] > left && mz[i] < right;

    // :52-58 advance until the first peak at or above the target.
    while *pos < n && mz[*pos] < target {
        *pos += 1;
    }
    let it = *pos;

    // :63-75 the current peak, or the last one when the walk passed the end.
    let current = if it == n { n - 1 } else { it };
    if inside(current) && accept(current) {
        integrated += intensity[current];
    }

    // :77-103 walk left. The first step may land on the first peak, which is
    // then added by the special case; the loop itself never adds index 0.
    let mut walker = it;
    if it != 0 {
        walker -= 1;
        if walker == 0 && inside(walker) && accept(walker) {
            integrated += intensity[walker];
        }
    }
    while walker != 0 && inside(walker) {
        if accept(walker) {
            integrated += intensity[walker];
        }
        walker -= 1;
    }

    // :105-119 walk right.
    let mut walker = it;
    if it != n {
        walker += 1;
    }
    while walker != n && inside(walker) {
        if accept(walker) {
            integrated += intensity[walker];
        }
        walker += 1;
    }
    integrated
}

fn check_walk_arrays(mz: usize, other: usize, what: &str, pos: usize) -> Result<()> {
    if other < mz {
        return Err(Error::InvalidValue(format!(
            "{what} array ({other}) is shorter than the m/z array ({mz})"
        )));
    }
    if pos > mz {
        return Err(Error::InvalidValue(format!(
            "walk position {pos} is past the end of the m/z array ({mz})"
        )));
    }
    Ok(())
}

impl ChromatogramExtractorAlgorithm {
    /// A new extractor. Source default constructor.
    pub fn new() -> Self {
        Self
    }

    /// Integrate the intensity within `mz ± mz_extraction_window / 2` and
    /// advance the walk position.
    ///
    /// Source `extract_value_tophat` (m/z only). The source's iterator triple
    /// `mz_start`, `mz_it`, `mz_end` becomes the slice `mz` with the index
    /// `pos` (the source's `mz_it`, which moves `int_it` along with it); the
    /// returned value is the source's `integrated_intensity`, which the source
    /// overwrites (not adds to) on every call.
    ///
    /// `mz_extraction_window` is the full window in Th, or in ppm when `ppm`
    /// is set (a 50 ppm window extracts 25 ppm on either side). The window is
    /// open: a peak exactly on either edge is not summed.
    ///
    /// # Notes
    ///
    /// The call moves `pos` forward to the first peak at or above `mz` and
    /// can **not** extract anything to the left of where `pos` already is, so
    /// all targets must be extracted in ascending m/z order. A descending
    /// target is not an error; it silently yields the wrong sum, as in the
    /// source. The source walk's two quirks (the first peak skipped by a long
    /// leftward walk, the last peak counted twice when the target is past the
    /// end) are reproduced; see the module documentation.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `intensity` is shorter than `mz` or `pos`
    /// is past the end of `mz`; the source reads out of bounds there.
    pub fn extract_value_tophat(
        &self,
        mz: &[f64],
        intensity: &[f64],
        pos: &mut usize,
        target_mz: f64,
        mz_extraction_window: f64,
        ppm: bool,
    ) -> Result<f64> {
        check_walk_arrays(mz.len(), intensity.len(), "intensity", *pos)?;
        let (left, right) = mz_window(target_mz, mz_extraction_window, ppm);
        Ok(tophat_walk(
            mz,
            intensity,
            pos,
            target_mz,
            left,
            right,
            |_| true,
        ))
    }

    /// Integrate the intensity within `mz ± mz_extraction_window / 2` and
    /// `im ± im_extraction_window / 2` and advance the walk position.
    ///
    /// Source `extract_value_tophat` (m/z and ion mobility). The spectrum is
    /// sorted by m/z with one ion-mobility value per peak (`ion_mobility`, the
    /// source's `im_it`); peaks outside the ion-mobility window are skipped,
    /// but the walk is bounded by m/z alone. Both windows are open. The rules
    /// and notes of [`ChromatogramExtractorAlgorithm::extract_value_tophat`]
    /// apply, including the ascending-order requirement.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `intensity` or `ion_mobility` is shorter
    /// than `mz`, or `pos` is past the end of `mz`.
    #[allow(clippy::too_many_arguments)]
    pub fn extract_value_tophat_im(
        &self,
        mz: &[f64],
        intensity: &[f64],
        ion_mobility: &[f64],
        pos: &mut usize,
        target_mz: f64,
        target_im: f64,
        mz_extraction_window: f64,
        im_extraction_window: f64,
        ppm: bool,
    ) -> Result<f64> {
        let mut im_pos = *pos;
        self.extract_value_tophat_im_lagged(
            mz,
            intensity,
            ion_mobility,
            pos,
            &mut im_pos,
            target_mz,
            target_im,
            mz_extraction_window,
            im_extraction_window,
            ppm,
        )
    }

    /// As [`Self::extract_value_tophat_im`], with the ion-mobility walk
    /// position kept apart from the m/z position, as the source keeps it.
    ///
    /// The source advances the ion-mobility iterator `im_it` only inside the
    /// ion-mobility overload of `extract_value_tophat`; the m/z-only overload
    /// advances `mz_it` and `int_it` and leaves `im_it` behind. A spectrum
    /// that mixes coordinates with and without ion mobility (a coordinate
    /// whose `ion_mobility` is negative takes the m/z-only overload) therefore
    /// reads the ion mobility of an earlier peak for every ion-mobility
    /// coordinate after the first m/z-only one, and the loss grows with each
    /// further m/z-only advance. `FeatureFinderIdentification` produces such
    /// spectra whenever a peptide has no ion-mobility annotation, and the
    /// features it reports depend on the lag, so it is reproduced here;
    /// [`ChromatogramExtractorAlgorithm::extract_chromatograms`] keeps one
    /// `im_pos` per spectrum.
    ///
    /// `im_pos` is the source's `im_it - im_start`; it never exceeds `pos`
    /// and advances by as many peaks as `pos` does in this call. The source
    /// reads before the start of the array where `im_pos` lags behind a walk
    /// that reaches the first peaks; those reads are undefined behaviour in
    /// the source and count here as outside the window.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `intensity` or `ion_mobility` is shorter
    /// than `mz`, `pos` is past the end of `mz`, or `im_pos` is past `pos`.
    #[allow(clippy::too_many_arguments)]
    pub fn extract_value_tophat_im_lagged(
        &self,
        mz: &[f64],
        intensity: &[f64],
        ion_mobility: &[f64],
        pos: &mut usize,
        im_pos: &mut usize,
        target_mz: f64,
        target_im: f64,
        mz_extraction_window: f64,
        im_extraction_window: f64,
        ppm: bool,
    ) -> Result<f64> {
        check_walk_arrays(mz.len(), intensity.len(), "intensity", *pos)?;
        check_walk_arrays(mz.len(), ion_mobility.len(), "ion mobility", *pos)?;
        if *im_pos > *pos {
            return Err(Error::InvalidValue(format!(
                "ion mobility walk position {} is past the m/z walk position {}",
                *im_pos, *pos
            )));
        }
        let (left, right) = mz_window(target_mz, mz_extraction_window, ppm);
        let left_im = target_im - im_extraction_window / 2.0;
        let right_im = target_im + im_extraction_window / 2.0;
        let lag = *pos - *im_pos;
        let before = *pos;
        let integrated = tophat_walk(mz, intensity, pos, target_mz, left, right, |i| {
            i >= lag && ion_mobility[i - lag] > left_im && ion_mobility[i - lag] < right_im
        });
        *im_pos += *pos - before;
        Ok(integrated)
    }

    /// Extract chromatograms at the m/z, ion mobility and RT of each
    /// coordinate, appending one point per matching spectrum.
    ///
    /// Source `extractChromatograms`. `output[k]` receives the chromatogram of
    /// `extraction_coordinates[k]`: for every spectrum whose RT lies in
    /// `[rt_start, rt_end]` (or every spectrum, when `rt_end - rt_start` is
    /// not positive) the spectrum RT is appended to its time array and the
    /// integrated intensity to its intensity array. Existing points are kept,
    /// as in the source, which appends too.
    ///
    /// `mz_extraction_window` is the full window in Th, or ppm when `ppm` is
    /// set. `im_extraction_window` is the full ion-mobility window; a positive
    /// value requests ion-mobility extraction for every coordinate with a
    /// non-negative `ion_mobility`, and then every non-empty spectrum must
    /// carry an ion-mobility array. `filter` is `"tophat"` (`"bartlett"` is
    /// declared but not implemented, as in the source).
    ///
    /// The coordinates must be sorted by ascending m/z, since one walk position
    /// is shared by all coordinates within a spectrum.
    ///
    /// # Errors
    ///
    /// All errors leave `output` unchanged; the source may have appended to
    /// some chromatograms before it throws.
    ///
    /// - [`Error::InvalidValue`] when `output` and `extraction_coordinates`
    ///   differ in length, when `filter` is unknown, when the coordinates are
    ///   not sorted by m/z (all three the source's `IllegalArgument`), when ion
    ///   mobility is requested and a non-empty spectrum has no ion-mobility
    ///   array (`IllegalArgument` too), or when a spectrum's intensity or
    ///   ion-mobility array is shorter than its m/z array (undefined in the
    ///   source).
    /// - [`Error::Unsupported`] for the Bartlett filter, when the first point
    ///   would be extracted with it (the source's `NotImplemented`).
    /// - [`Error::InvalidRange`] when the spectrum count times the coordinate
    ///   count exceeds [`MAX_EXTRACTED_POINTS`].
    /// - The errors of the spectrum access.
    ///
    /// As in the source, an input without spectra returns `Ok(())` before any
    /// argument is checked.
    #[allow(clippy::too_many_arguments)]
    pub fn extract_chromatograms(
        &self,
        input: &dyn SpectrumAccess,
        output: &mut [Chromatogram],
        extraction_coordinates: &[ExtractionCoordinates],
        mz_extraction_window: f64,
        ppm: bool,
        im_extraction_window: f64,
        filter: &str,
    ) -> Result<()> {
        self.extract_reporting(
            input,
            output,
            extraction_coordinates,
            mz_extraction_window,
            ppm,
            im_extraction_window,
            filter,
            &mut ProgressReporter::silent(),
        )
    }

    /// [`ChromatogramExtractorAlgorithm::extract_chromatograms`], reporting to
    /// `progress` as the source's `ProgressLogger` base does:
    /// `startProgress(0, n, "Extracting chromatograms")`, `setProgress(i)`
    /// before spectrum `i`, `endProgress()` (`:275-278`, `:352`).
    ///
    /// # Errors
    ///
    /// As [`ChromatogramExtractorAlgorithm::extract_chromatograms`], and the
    /// errors of `progress`.
    #[allow(clippy::too_many_arguments)]
    pub fn extract_chromatograms_with_progress(
        &self,
        input: &dyn SpectrumAccess,
        output: &mut [Chromatogram],
        extraction_coordinates: &[ExtractionCoordinates],
        mz_extraction_window: f64,
        ppm: bool,
        im_extraction_window: f64,
        filter: &str,
        progress: &mut ProgressLogger,
    ) -> Result<()> {
        self.extract_reporting(
            input,
            output,
            extraction_coordinates,
            mz_extraction_window,
            ppm,
            im_extraction_window,
            filter,
            &mut ProgressReporter::new(Some(progress)),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn extract_reporting(
        &self,
        input: &dyn SpectrumAccess,
        output: &mut [Chromatogram],
        coordinates: &[ExtractionCoordinates],
        mz_extraction_window: f64,
        ppm: bool,
        im_extraction_window: f64,
        filter: &str,
        reporter: &mut ProgressReporter<'_>,
    ) -> Result<()> {
        let input_size = input.nr_spectra();
        if input_size < 1 {
            return Ok(());
        }
        if output.len() != coordinates.len() {
            return Err(Error::InvalidValue(format!(
                "Output and extraction coordinates need to have the same size: {} != {}",
                output.len(),
                coordinates.len()
            )));
        }
        let used_filter = ExtractionFilter::from_name(filter)?;
        if coordinates
            .windows(2)
            .any(|pair| ExtractionCoordinates::sort_reverse_by_mz(&pair[0], &pair[1]))
        {
            return Err(Error::InvalidValue(
                "Input to extractChromatogram needs to be sorted by m/z".into(),
            ));
        }
        input_size
            .checked_mul(coordinates.len())
            .filter(|&pairs| pairs <= MAX_EXTRACTED_POINTS)
            .ok_or_else(|| {
                Error::InvalidRange(format!(
                    "{input_size} spectra x {} coordinates exceed the extraction ceiling of {MAX_EXTRACTED_POINTS} points",
                    coordinates.len()
                ))
            })?;

        // Build the appended points aside, commit only on success.
        let mut times: Vec<Vec<f64>> = vec![Vec::new(); coordinates.len()];
        let mut intensities: Vec<Vec<f64>> = vec![Vec::new(); coordinates.len()];
        let has_im = im_extraction_window > 0.0;
        let records = progress_value(input_size)?;
        reporter.section(0, records, EXTRACTION_PROGRESS_LABEL, |reporter| {
            for scan_idx in 0..input_size {
                reporter.set_count(scan_idx)?;
                let spectrum = input.spectrum_by_id(scan_idx)?;
                let meta = input.spectrum_meta_by_id(scan_idx)?;
                let mz = &spectrum.mz_array().data;
                let intensity = &spectrum.intensity_array().data;
                if mz.is_empty() {
                    continue;
                }
                if intensity.len() < mz.len() {
                    return Err(Error::InvalidValue(format!(
                        "spectrum {scan_idx}: intensity array ({}) is shorter than the m/z array ({})",
                        intensity.len(),
                        mz.len()
                    )));
                }
                let im: &[f64] = if has_im {
                    match spectrum.drift_time_array() {
                        Some(array) => &array.data,
                        None => {
                            return Err(Error::InvalidValue(
                                "Requested ion mobility extraction but no ion mobility array found."
                                    .into(),
                            ));
                        }
                    }
                } else {
                    &[]
                };
                if has_im && im.len() < mz.len() {
                    return Err(Error::InvalidValue(format!(
                        "spectrum {scan_idx}: ion mobility array ({}) is shorter than the m/z array ({})",
                        im.len(),
                        mz.len()
                    )));
                }
                let current_rt = meta.rt;
                let mut pos = 0usize;
                // the source's `im_it`: advanced by the ion-mobility walk only
                let mut im_pos = 0usize;
                for (k, coord) in coordinates.iter().enumerate() {
                    if coord.rt_end - coord.rt_start > 0.0
                        && (current_rt < coord.rt_start || current_rt > coord.rt_end)
                    {
                        continue;
                    }
                    let use_im = coord.ion_mobility >= 0.0 && has_im;
                    let integrated = match used_filter {
                        ExtractionFilter::Tophat if use_im => self.extract_value_tophat_im_lagged(
                            mz,
                            intensity,
                            im,
                            &mut pos,
                            &mut im_pos,
                            coord.mz,
                            coord.ion_mobility,
                            mz_extraction_window,
                            im_extraction_window,
                            ppm,
                        )?,
                        ExtractionFilter::Tophat => self.extract_value_tophat(
                            mz,
                            intensity,
                            &mut pos,
                            coord.mz,
                            mz_extraction_window,
                            ppm,
                        )?,
                        ExtractionFilter::Bartlett => {
                            return Err(Error::Unsupported(
                                "the Bartlett extraction filter is not implemented".into(),
                            ));
                        }
                    };
                    times[k].push(current_rt);
                    intensities[k].push(integrated);
                }
            }
            Ok(())
        })?;

        for ((chromatogram, time), intensity) in output.iter_mut().zip(times).zip(intensities) {
            chromatogram.time_array_mut().data.extend(time);
            chromatogram.intensity_array_mut().data.extend(intensity);
        }
        Ok(())
    }
}
