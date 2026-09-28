// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Adding up individual spectra.
//!
//! Ports `ANALYSIS/OPENSWATH/SpectrumAddition.h` with `SpectrumAddition.cpp`
//! (a class of static members, which become free functions). The spectra are
//! resampled in m/z with a given sampling rate and then added up, which may
//! lose accuracy, especially with an inappropriate sampling rate. See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_DATAACCESS_SUPPORT.md).
//!
//! # Resampling
//!
//! Resampling is the source's `LinearResamplerAlign::raster` over a
//! prepopulated grid, reproduced here for the two element types the source
//! instantiates it with (`double` arrays and `Peak1D`): every input intensity
//! is split between the two neighbouring grid points in proportion to the
//! distance, so the total intensity is conserved; intensity below the first
//! grid point goes to the first, and intensity at or beyond the last grid point
//! to the last. The crate's
//! [`LinearResamplerAlign`](crate::processing::LinearResamplerAlign) builds its
//! own grid and does not accumulate into an existing one, so it cannot be
//! reused here. The source's `verifySpacing_` only logs a warning, and it
//! compares the resampler's own default spacing (not the sampling rate used
//! here); the port does not log it.
//!
//! # Ceilings
//!
//! The grid has `(max - min) / sampling_rate + 1` points (truncated). The source
//! converts that to `int` unchecked, so a zero or tiny sampling rate is
//! undefined behaviour; the port refuses a non-finite or non-positive rate and
//! a grid above [`MAX_RESAMPLED_POINTS`], before allocating.

use super::data_structures::{BinaryDataArray, Spectrum, SpectrumPtr};
use super::spectrum_access::{SpectrumSequence, filter_by_drift};
use crate::kernel::ranges::RangeBase;
use crate::kernel::{MSSpectrum, Peak1D};
use crate::{Error, Result};
use std::cmp::Ordering;
use std::sync::Arc;

/// Upper bound on the number of points of a resampling grid (ten million).
pub const MAX_RESAMPLED_POINTS: usize = 10_000_000;

fn grid_size(min: f64, max: f64, sampling_rate: f64) -> Result<usize> {
    if !sampling_rate.is_finite() || sampling_rate <= 0.0 {
        return Err(Error::InvalidValue(format!(
            "the sampling rate must be positive and finite, got {sampling_rate} (a single-point spectrum cannot estimate one)"
        )));
    }
    let points = (max - min) / sampling_rate + 1.0;
    if !points.is_finite() || points < 1.0 || points > MAX_RESAMPLED_POINTS as f64 {
        return Err(Error::InvalidValue(format!(
            "a resampling grid of {points} points is outside [1, {MAX_RESAMPLED_POINTS}]"
        )));
    }
    // Truncation toward zero, as the source's conversion to int.
    Ok(points as usize)
}

/// `LinearResamplerAlign::raster` over separate `double` arrays, adding into
/// `grid_int`.
fn raster_arrays(raw_mz: &[f64], raw_int: &[f64], grid_mz: &[f64], grid_int: &mut [f64]) {
    let n = grid_mz.len();
    if n == 0 {
        return;
    }
    let mut r = 0;
    while r < raw_mz.len() && raw_mz[r] < grid_mz[0] {
        grid_int[0] += raw_int[r];
        r += 1;
    }
    let mut j = 0;
    while r < raw_mz.len() {
        while j < n && grid_mz[j] < raw_mz[r] {
            j += 1;
        }
        if j != 0 {
            j -= 1;
        }
        if j + 1 == n {
            break;
        }
        let dist_left = (raw_mz[r] - grid_mz[j]).abs();
        let dist_right = (raw_mz[r] - grid_mz[j + 1]).abs();
        grid_int[j] += raw_int[r] * dist_right / (dist_left + dist_right);
        grid_int[j + 1] += raw_int[r] * dist_left / (dist_left + dist_right);
        r += 1;
    }
    while r < raw_mz.len() {
        grid_int[j] += raw_int[r];
        r += 1;
    }
}

/// `LinearResamplerAlign::raster` over `Peak1D`, adding into `grid`.
///
/// Intensities are `f32`: an addition of two intensities is an `f32` addition
/// and a weighted share is computed in `f64` and stored as `f32`, as the
/// source's `setIntensity(getIntensity() + ...)` does.
fn raster_peaks(raw: &[Peak1D], grid: &mut [Peak1D]) {
    let n = grid.len();
    if n == 0 {
        return;
    }
    let mut r = 0;
    while r < raw.len() && raw[r].mz < grid[0].mz {
        grid[0].intensity += raw[r].intensity;
        r += 1;
    }
    let mut j = 0;
    while r < raw.len() {
        while j < n && grid[j].mz < raw[r].mz {
            j += 1;
        }
        if j != 0 {
            j -= 1;
        }
        if j + 1 == n {
            break;
        }
        let dist_left = (raw[r].mz - grid[j].mz).abs();
        let dist_right = (raw[r].mz - grid[j + 1].mz).abs();
        let intensity = f64::from(raw[r].intensity);
        grid[j].intensity = (f64::from(grid[j].intensity)
            + intensity * dist_right / (dist_left + dist_right)) as f32;
        grid[j + 1].intensity = (f64::from(grid[j + 1].intensity)
            + intensity * dist_left / (dist_left + dist_right))
            as f32;
        r += 1;
    }
    while r < raw.len() {
        grid[j].intensity += raw[r].intensity;
        r += 1;
    }
}

// Sort key consistent with `<` on every non-NaN value (it merges -0 and +0,
// which `<` treats as equal) and total, so a NaN cannot break the sort.
fn mz_key(a: f64, b: f64) -> Ordering {
    let norm = |x: f64| if x == 0.0 { 0.0 } else { x };
    norm(a).total_cmp(&norm(b))
}

/// Sort every data array of `spec` by the m/z array, stably.
///
/// Equal m/z values keep their order. Empty arrays are skipped, as in the
/// source, so an empty ion-mobility placeholder does not block the sort.
///
/// The source sorts with `<`, which is undefined for a NaN m/z; the port sorts
/// NaN (with a clear sign bit) after every number.
///
/// # Errors
///
/// [`Error::InvalidValue`] when a non-empty array is not as long as the m/z
/// array (the source then reads out of bounds); `spec` is then unchanged.
pub fn sort_spectrum_by_mz(spec: &mut Spectrum) -> Result<()> {
    let n = spec.mz_array().data.len();
    if let Some(bad) = spec
        .data_arrays()
        .iter()
        .find(|a| !a.data.is_empty() && a.data.len() != n)
    {
        return Err(Error::InvalidValue(format!(
            "data array '{}' has {} entries for {n} m/z values",
            bad.description,
            bad.data.len()
        )));
    }
    let mz = &spec.mz_array().data;
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| mz_key(mz[a], mz[b]));
    for array in spec.data_arrays_mut() {
        if array.data.is_empty() {
            continue;
        }
        array.data = order.iter().map(|&i| array.data[i]).collect();
    }
    Ok(())
}

fn check_spectrum_pair(spectrum: &Spectrum) -> Result<()> {
    if spectrum.mz_array().data.len() != spectrum.intensity_array().data.len() {
        return Err(Error::InvalidValue(format!(
            "the m/z array has {} entries but the intensity array {}",
            spectrum.mz_array().data.len(),
            spectrum.intensity_array().data.len()
        )));
    }
    Ok(())
}

/// Add up a sequence of spectra by resampling them onto a common grid and
/// adding their intensities.
///
/// - An empty sequence yields an empty spectrum; a single spectrum is returned
///   as it is (the same handle, neither resampled nor filtered).
/// - When the **first** spectrum is empty the result is empty, whatever the
///   other spectra hold, as in the source.
/// - The grid runs from the smallest first m/z to the largest last m/z of the
///   non-empty spectra (each is assumed sorted), in steps of `sampling_rate`.
///   A negative `sampling_rate` asks for an estimate: the smallest distance
///   between neighbouring m/z values in any spectrum, starting from the first
///   spectrum's m/z span.
/// - With `filter_zeros`, only grid points with a positive intensity are kept.
///
/// Only the m/z and intensity arrays are used; further arrays are dropped.
///
/// # Errors
///
/// [`Error::InvalidValue`] when a spectrum's intensity array is not as long as
/// its m/z array, or when the sampling rate (given or estimated) or the grid
/// size is refused (see the module notes).
pub fn add_up_spectra(
    all_spectra: &[SpectrumPtr],
    sampling_rate: f64,
    filter_zeros: bool,
) -> Result<SpectrumPtr> {
    if all_spectra.len() == 1 {
        return Ok(Arc::clone(&all_spectra[0]));
    }
    let Some(first) = all_spectra.first() else {
        return Ok(Arc::new(Spectrum::new()));
    };
    if first.mz_array().data.is_empty() {
        return Ok(Arc::new(Spectrum::new()));
    }
    for spectrum in all_spectra {
        check_spectrum_pair(spectrum)?;
    }
    let first_mz = &first.mz_array().data;
    let mut min = first_mz[0];
    let mut max = first_mz[first_mz.len() - 1];
    let mut min_spacing = max - min;
    for spectrum in all_spectra {
        let mz = &spectrum.mz_array().data;
        if mz.is_empty() {
            continue;
        }
        if sampling_rate < 0.0 {
            for pair in mz.windows(2) {
                if min_spacing > pair[1] - pair[0] {
                    min_spacing = pair[1] - pair[0];
                }
            }
        }
        if mz[0] < min {
            min = mz[0];
        }
        if mz[mz.len() - 1] > max {
            max = mz[mz.len() - 1];
        }
    }
    let rate = if sampling_rate < 0.0 {
        min_spacing
    } else {
        sampling_rate
    };
    let points = grid_size(min, max, rate)?;
    let grid_mz: Vec<f64> = (0..points).map(|i| min + i as f64 * rate).collect();
    let mut grid_int = vec![0.0; points];
    for spectrum in all_spectra {
        raster_arrays(
            &spectrum.mz_array().data,
            &spectrum.intensity_array().data,
            &grid_mz,
            &mut grid_int,
        );
    }
    if !filter_zeros {
        return Ok(Arc::new(Spectrum::from_arrays(grid_mz, grid_int)));
    }
    let (mz, intensity): (Vec<f64>, Vec<f64>) = grid_mz
        .into_iter()
        .zip(grid_int)
        .filter(|&(_, i)| i > 0.0)
        .unzip();
    Ok(Arc::new(Spectrum::from_arrays(mz, intensity)))
}

/// Add up ion-mobility spectra: filter each to the drift-time window
/// `im_range` (inclusive), then [`add_up_spectra`].
///
/// Resampling is not supported on three-dimensional data, hence the filter
/// first. An empty sequence yields an empty spectrum; an empty `im_range`
/// skips the filter. Note that a single spectrum is filtered (a copy) before
/// being returned, unlike in [`add_up_spectra`].
///
/// # Errors
///
/// Those of [`filter_by_drift`] (a spectrum without an ion-mobility array) and
/// of [`add_up_spectra`].
pub fn add_up_spectra_in_mobility_range(
    all_spectra: &[SpectrumPtr],
    im_range: &RangeBase,
    sampling_rate: f64,
    filter_zeros: bool,
) -> Result<SpectrumPtr> {
    if all_spectra.is_empty() {
        return Ok(Arc::new(Spectrum::new()));
    }
    if im_range.is_empty() {
        return add_up_spectra(all_spectra, sampling_rate, filter_zeros);
    }
    let (lo, hi) = (im_range.min()?, im_range.max()?);
    let filtered = all_spectra
        .iter()
        .map(|s| filter_by_drift(s, lo, hi).map(Arc::new))
        .collect::<Result<SpectrumSequence>>()?;
    add_up_spectra(&filtered, sampling_rate, filter_zeros)
}

/// Concatenate a spectrum sequence into a single spectrum, sorted by m/z.
///
/// The result has as many data arrays as the first spectrum, the extra ones
/// under the first spectrum's descriptions; array `k` of every spectrum is
/// appended to array `k` of the result. No resampling takes place.
///
/// # Errors
///
/// [`Error::InvalidValue`] when a spectrum has more data arrays than the first
/// (the source writes out of bounds), and the errors of
/// [`sort_spectrum_by_mz`] (arrays whose lengths do not line up).
pub fn concatenate_spectra(all_spectra: &[SpectrumPtr]) -> Result<SpectrumPtr> {
    let mut added = Spectrum::new();
    if let Some(first) = all_spectra.first() {
        for extra in first.data_arrays().iter().skip(2) {
            added.push_data_array(BinaryDataArray::with_description(
                Vec::new(),
                extra.description.clone(),
            ));
        }
    }
    let width = added.data_arrays().len();
    for spectrum in all_spectra {
        if spectrum.data_arrays().len() > width {
            return Err(Error::InvalidValue(format!(
                "a spectrum has {} data arrays, more than the first spectrum's {width}",
                spectrum.data_arrays().len()
            )));
        }
        for (target, source) in added
            .data_arrays_mut()
            .iter_mut()
            .zip(spectrum.data_arrays())
        {
            target.data.extend_from_slice(&source.data);
        }
    }
    sort_spectrum_by_mz(&mut added)?;
    Ok(Arc::new(added))
}

/// Add up OpenMS spectra by resampling them and adding their intensities.
///
/// The rules of [`add_up_spectra`] apply, with two source differences kept:
///
/// - with a negative `sampling_rate` (estimation), the m/z range also takes
///   every peak into account, not only the first and last of each spectrum;
/// - intensities are `f32`: each spectrum is resampled onto its own zeroed grid
///   and then added to the running total in `f32`.
///
/// A single spectrum is returned as a clone; the result of a resampling carries
/// no metadata (a default [`MSSpectrum`]). The source's debug-only
/// precondition that the spectra carry no float data arrays (no ion mobility)
/// is not checked; the arrays are ignored, as in a release build.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the sampling rate (given or estimated) or the
/// grid size is refused.
pub fn add_up_ms_spectra(
    all_spectra: &[MSSpectrum],
    sampling_rate: f64,
    filter_zeros: bool,
) -> Result<MSSpectrum> {
    if all_spectra.len() == 1 {
        return Ok(all_spectra[0].clone());
    }
    let Some(first) = all_spectra.first() else {
        return Ok(MSSpectrum::default());
    };
    let (Some(first_peak), Some(last_peak)) = (first.peaks.first(), first.peaks.last()) else {
        return Ok(MSSpectrum::default());
    };
    let mut min = first_peak.mz;
    let mut max = last_peak.mz;
    let mut min_spacing = max - min;
    for spectrum in all_spectra {
        let peaks = &spectrum.peaks;
        let (Some(head), Some(tail)) = (peaks.first(), peaks.last()) else {
            continue;
        };
        if sampling_rate < 0.0 {
            for (k, peak) in peaks.iter().enumerate() {
                if k > 0 && min_spacing > peak.mz - peaks[k - 1].mz {
                    min_spacing = peak.mz - peaks[k - 1].mz;
                }
                if peak.mz < min {
                    min = peak.mz;
                }
                if peak.mz > max {
                    max = peak.mz;
                }
            }
        }
        if head.mz < min {
            min = head.mz;
        }
        if tail.mz > max {
            max = tail.mz;
        }
    }
    let rate = if sampling_rate < 0.0 {
        min_spacing
    } else {
        sampling_rate
    };
    let points = grid_size(min, max, rate)?;
    let grid: Vec<Peak1D> = (0..points)
        .map(|i| Peak1D::new(min + i as f64 * rate, 0.0))
        .collect();
    let mut master = grid.clone();
    for spectrum in all_spectra {
        let mut output = grid.clone();
        raster_peaks(&spectrum.peaks, &mut output);
        for (m, o) in master.iter_mut().zip(&output) {
            m.intensity += o.intensity;
        }
    }
    let peaks = if filter_zeros {
        master.into_iter().filter(|p| p.intensity > 0.0).collect()
    } else {
        master
    };
    Ok(MSSpectrum {
        peaks,
        ..MSSpectrum::default()
    })
}
