// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Biosaur2 hills: linking peaks across consecutive scans into m/z traces
//! (`detectHills_`, `linkScanToHills_`, `computeHillMzStep_`), the summary
//! statistics and length filter (`processHills_`), the valley split
//! (`splitHills_`), and the small numerical helpers these stages share
//! (`meanFilter_`, `calculateMedian_`, `calculatePPM_`, `calibrateMass_`).
//!
//! The helpers are private members of the source class. They are public here,
//! as free functions, so that the differential tests can compare each one with
//! the Release build on its own.

use super::{Biosaur2Algorithm, Hill, source_invalid_value};
use crate::kernel::MSSpectrum;
use crate::math::libstdcxx::std_max;
use crate::math::source_sort::source_sort_by;
use crate::math::statistic_functions::median;
use crate::metadata::{ImTypes, IonMobilityFormat};
use crate::processing::noise_estimation::x86::cvttsd2si32;
use crate::{Error, Result};
use std::collections::BTreeMap;

/// Mass difference of `mz1` from the reference `mz2`, in ppm, the source
/// `calculatePPM_`: `(mz1 - mz2) / mz2 * 1e6`, and 0 when `|mz2| < 1e-10`.
pub fn calculate_ppm(mz1: f64, mz2: f64) -> f64 {
    if mz2.abs() < 1e-10 {
        return 0.0;
    }
    (mz1 - mz2) / mz2 * 1e6
}

/// Median of `values`, 0 for an empty input, the source `calculateMedian_`
/// (`Math::median` on a sorted copy; an even count averages the two middle
/// values).
///
/// # Errors
///
/// Only when the sort permutation cannot be allocated.
pub fn calculate_median(values: &[f64]) -> Result<f64> {
    if values.is_empty() {
        return Ok(0.0);
    }
    let mut sorted = values.to_vec();
    median(&mut sorted)
}

/// Moving-average filter with zero padding, the source `meanFilter_`.
///
/// `window` is the half-width of the kernel (kernel length `2 * window + 1`),
/// mirroring the Python `meanfilt` with NumPy's `same` convolution: every
/// output is the sum of the in-range inputs of its window divided by the full
/// kernel length.
pub fn mean_filter(data: &[f64], window: usize) -> Vec<f64> {
    let mut result = vec![0.0; data.len()];
    if data.is_empty() {
        return result;
    }
    let kernel_len = window.saturating_mul(2).saturating_add(1);
    for (i, out) in result.iter_mut().enumerate() {
        let mut sum = 0.0;
        // j runs from i - window to i + window; zeros outside the data.
        let first = i.saturating_sub(window);
        let last = (i.saturating_add(window)).min(data.len() - 1);
        for value in &data[first..=last] {
            sum += value;
        }
        *out = sum / kernel_len as f64;
    }
    result
}

/// Native ceiling on the histogram of [`calibrate_mass`]: 2^24 bins (128 MiB of
/// counters).
pub const MAX_CALIBRATION_BINS: usize = 1 << 24;

/// Estimate a systematic shift and a spread from a distribution of mass
/// errors (ppm), the source `calibrateMass_`.
///
/// The errors are binned (`bin_width`, the source default 0.05 ppm) between
/// their minimum and maximum, and the mean and standard deviation of the bin
/// centres weighted by the bin counts are returned, the variance floored at
/// 0.01. `(0, 10)` is returned for no errors, fewer than 5 bins, fewer than 10
/// binned errors, or a non-finite spread; a mean at or beyond the range's
/// larger bound retries with 0.25 ppm bins.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] where the source recurses without end: the
/// 0.25 ppm retry itself asks for a retry. No finite input reaches that. Also,
/// natively, when the histogram would need more than
/// [`MAX_CALIBRATION_BINS`] bins (an error range above 800,000 ppm at the
/// default width; the errors come from linking within `htol` or `itol`) or
/// cannot be allocated; the source allocates one counter per bin, up to
/// `i32::MAX` of them.
pub fn calibrate_mass(mass_errors: &[f64], bin_width: f64) -> Result<(f64, f64)> {
    if mass_errors.is_empty() {
        return Ok((0.0, 10.0));
    }
    let mut min_error = mass_errors[0];
    let mut max_error = mass_errors[0];
    for &error in &mass_errors[1..] {
        if error < min_error {
            min_error = error;
        }
        if max_error < error {
            max_error = error;
        }
    }
    let mass_left = -min_error;
    let mass_right = max_error;
    let n_bins = cvttsd2si32((mass_left + mass_right) / bin_width);
    if n_bins < 5 {
        return Ok((0.0, 10.0));
    }
    let bins = n_bins as usize;
    if bins > MAX_CALIBRATION_BINS {
        return Err(Error::InvalidValue(format!(
            "calibrateMass_ would need {bins} bins, more than the {MAX_CALIBRATION_BINS} allowed"
        )));
    }
    let mut bin_counts: Vec<u64> = Vec::new();
    bin_counts
        .try_reserve_exact(bins)
        .map_err(|_| Error::InvalidValue(format!("cannot allocate {bins} calibration bins")))?;
    bin_counts.resize(bins, 0);
    for &error in mass_errors {
        let bin = cvttsd2si32((error + mass_left) / bin_width);
        if bin >= 0 && bin < n_bins {
            bin_counts[bin as usize] += 1;
        }
    }
    let mut sum_x = 0.0;
    let mut sum_x2 = 0.0;
    let mut sum_w = 0.0;
    for (i, &count) in bin_counts.iter().enumerate() {
        let x = -mass_left + (i as f64 + 0.5) * bin_width;
        let w = count as f64;
        sum_x += w * x;
        sum_x2 += w * x * x;
        sum_w += w;
    }
    if sum_w < 10.0 {
        return Ok((0.0, 10.0));
    }
    let mean = sum_x / sum_w;
    let variance = (sum_x2 / sum_w) - (mean * mean);
    let sigma = (if variance < 0.01 { 0.01 } else { variance }).sqrt();
    if mean.abs() >= std_max(mass_left, mass_right) {
        if bin_width == 0.25 {
            return Err(Error::InvalidValue(
                "calibrateMass_ would retry with 0.25 ppm bins without end".into(),
            ));
        }
        return calibrate_mass(mass_errors, 0.25);
    }
    if sigma.is_infinite() || sigma.is_nan() {
        return Ok((0.0, 10.0));
    }
    Ok((mean, sigma))
}

/// Whether a peak passes the intensity and m/z window.
fn in_window(algorithm: &Biosaur2Algorithm, mz: f64, intensity: f64) -> bool {
    !(intensity < algorithm.mini || mz < algorithm.minmz || mz > algorithm.maxmz)
}

/// The m/z bin width for hill linking, the source `computeHillMzStep_`:
/// `htol_ppm * 1e-6` times the largest m/z of a peak inside the window, or 0
/// when there is none or `htol_ppm` is not positive.
pub fn compute_hill_mz_step(
    algorithm: &Biosaur2Algorithm,
    spectra: &[&MSSpectrum],
    htol_ppm: f64,
) -> f64 {
    let mut max_mz_value = 0.0;
    for spectrum in spectra {
        for peak in &spectrum.peaks {
            if !in_window(algorithm, peak.mz, f64::from(peak.intensity)) {
                continue;
            }
            if peak.mz > max_mz_value {
                max_mz_value = peak.mz;
            }
        }
    }
    if max_mz_value <= 0.0 || htol_ppm <= 0.0 {
        return 0.0;
    }
    htol_ppm * 1e-6 * max_mz_value
}

/// The linking state carried from one scan to the next.
#[derive(Default)]
struct LinkState<'a> {
    /// Hill position of each peak of the previous scan (`usize::MAX`: none).
    prev_peak_to_hill: Vec<usize>,
    prev_spectrum: Option<&'a MSSpectrum>,
    /// m/z bin -> peaks of the previous scan, in descending intensity.
    prev_fast_dict: BTreeMap<i32, Vec<usize>>,
    /// Ion-mobility bin of each peak of the previous scan.
    prev_im_bins: Vec<i32>,
}

/// Detect hills (continuous m/z traces) in `spectra`, the source
/// `detectHills_`.
///
/// Scans through the spectra and groups peaks with similar m/z values across
/// consecutive scans into hills (see the source's `linkScanToHills_` for the rule);
/// with `mass_diffs` the signed ppm difference of every link is collected for
/// the hill calibration. `use_im` enables ion-mobility gating (2-D m/z-IM
/// hills) for spectra with per-peak ion mobility. The hills are returned
/// unprocessed: only the per-point vectors, the length, the retention-time
/// bounds and apex and the intensity sum are set.
///
/// # Errors
///
/// As the source's `linkScanToHills_`.
pub fn detect_hills(
    algorithm: &Biosaur2Algorithm,
    spectra: &[&MSSpectrum],
    htol_ppm: f64,
    use_im: bool,
    mut mass_diffs: Option<&mut Vec<f64>>,
) -> Result<Vec<Hill>> {
    let mz_step = compute_hill_mz_step(algorithm, spectra, htol_ppm);
    if let Some(diffs) = mass_diffs.as_deref_mut() {
        diffs.clear();
    }
    let mut hills = Vec::new();
    let mut state = LinkState::default();
    for (scan_idx, spectrum) in spectra.iter().enumerate() {
        link_scan_to_hills(
            algorithm,
            spectrum,
            scan_idx,
            htol_ppm,
            mz_step,
            use_im,
            &mut hills,
            &mut state,
            mass_diffs.as_deref_mut(),
        )?;
    }
    Ok(hills)
}

/// The per-point ion mobility of a spectrum's peaks, when it has a per-peak
/// ion-mobility array.
fn im_values(spectrum: &MSSpectrum) -> Result<Option<&[f32]>> {
    if ImTypes::determine_im_format(spectrum) != IonMobilityFormat::PerPeak {
        return Ok(None);
    }
    let (index, _) = spectrum.im_data()?;
    Ok(Some(&spectrum.float_data_arrays[index].data))
}

fn im_too_short(peak_idx: usize) -> Error {
    source_invalid_value(
        "Ion mobility array shorter than peak list.",
        &peak_idx.to_string(),
    )
}

/// Link the peaks of one scan to the hills of the previous scan or start new
/// hills, the source `linkScanToHills_`.
///
/// The peaks inside the intensity and m/z window are visited in descending
/// intensity (the source's `std::sort`, ties as the Release build orders
/// them). Each looks up the previous scan's peaks in its own and the two
/// neighbouring m/z bins (`mz_step` wide), keeps, when ion-mobility gating is
/// active for this scan, those at most one ion-mobility bin (`paseftol` wide)
/// away, and links to the most intense of them within `htol_ppm` (the last
/// one wins a tie) that no stronger peak of this scan has claimed. A peak
/// without a partner starts a new hill. A scan without peaks in the window
/// breaks every hill.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] with the source's text `Ion mobility array
/// shorter than peak list.` when a peak in the window has no entry in the
/// spectrum's ion-mobility array; and natively when a peak in the window has
/// a non-finite m/z or intensity, which the source sorts and bins with
/// undefined behaviour.
#[allow(clippy::too_many_arguments)]
fn link_scan_to_hills<'a>(
    algorithm: &Biosaur2Algorithm,
    spectrum: &'a MSSpectrum,
    scan_idx: usize,
    htol_ppm: f64,
    mz_step: f64,
    use_im_global: bool,
    hills: &mut Vec<Hill>,
    state: &mut LinkState<'a>,
    mut mass_diffs: Option<&mut Vec<f64>>,
) -> Result<()> {
    let rt = spectrum.rt;
    let mut valid = Vec::new();
    for (index, peak) in spectrum.peaks.iter().enumerate() {
        let intensity = f64::from(peak.intensity);
        if !in_window(algorithm, peak.mz, intensity) {
            continue;
        }
        if !peak.mz.is_finite() || !intensity.is_finite() {
            return Err(Error::InvalidValue(format!(
                "Biosaur2 needs finite peaks: spectrum at RT {rt} has a non-finite m/z or intensity at peak {index}"
            )));
        }
        valid.push(index);
    }
    if valid.is_empty() {
        state.prev_fast_dict.clear();
        state.prev_peak_to_hill.clear();
        state.prev_spectrum = None;
        return Ok(());
    }

    let n_peaks = spectrum.peaks.len();
    let intensities: Vec<f64> = valid
        .iter()
        .map(|&i| f64::from(spectrum.peaks[i].intensity))
        .collect();
    let im = im_values(spectrum)?;
    let use_im_current = use_im_global && im.is_some();
    let mut im_bin_per_peak = vec![0i32; n_peaks];
    if use_im_current {
        if let Some(values) = im {
            for &peak_idx in &valid {
                let value = values.get(peak_idx).ok_or_else(|| im_too_short(peak_idx))?;
                im_bin_per_peak[peak_idx] = if algorithm.paseftol > 0.0 {
                    cvttsd2si32(f64::from(*value) / algorithm.paseftol)
                } else {
                    0
                };
            }
        }
    }

    let mut order: Vec<usize> = (0..valid.len()).collect();
    source_sort_by(&mut order, |&a, &b| intensities[a] > intensities[b])?;
    let basic_id_sorted: Vec<usize> = order.iter().map(|&local| valid[local]).collect();
    let fast_array: Vec<i32> = basic_id_sorted
        .iter()
        .map(|&peak_idx| {
            if mz_step > 0.0 {
                cvttsd2si32(spectrum.peaks[peak_idx].mz / mz_step)
            } else {
                0
            }
        })
        .collect();
    let mut fast_dict: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
    for (&bin, &peak_idx) in fast_array.iter().zip(&basic_id_sorted) {
        fast_dict.entry(bin).or_default().push(peak_idx);
    }

    let mut curr_peak_to_hill = vec![usize::MAX; n_peaks];
    let prev_len = state.prev_spectrum.map_or(0, MSSpectrum::len);
    let mut banned_prev = vec![false; prev_len];

    for (&idx, &fm) in basic_id_sorted.iter().zip(&fast_array) {
        let fi = if use_im_current {
            im_bin_per_peak[idx]
        } else {
            0
        };
        let mut assigned_hill = usize::MAX;
        let keys = [fm, fm.wrapping_sub(1), fm.wrapping_add(1)];
        if keys
            .iter()
            .any(|key| state.prev_fast_dict.contains_key(key))
        {
            let mut all_idx: Vec<usize> = Vec::new();
            for key in &keys {
                if let Some(peaks) = state.prev_fast_dict.get(key) {
                    all_idx.extend_from_slice(peaks);
                }
            }
            if use_im_current && !state.prev_im_bins.is_empty() {
                let prev_im_bins = &state.prev_im_bins;
                all_idx.retain(|&prev| {
                    prev < prev_im_bins.len()
                        && (i64::from(prev_im_bins[prev]) - i64::from(fi)).abs() <= 1
                });
            }

            let mut best_intensity = 0.0;
            let mut best_prev: Option<usize> = None;
            let mut best_mass_diff_with_sign = 0.0;
            let mz_cur = spectrum.peaks[idx].mz;
            if let Some(prev_spectrum) = state.prev_spectrum {
                for &prev in &all_idx {
                    if prev >= prev_spectrum.len() || banned_prev[prev] {
                        continue;
                    }
                    let prev_peak = prev_spectrum.peaks[prev];
                    let cur_intensity = f64::from(prev_peak.intensity);
                    let mass_diff_with_sign = (mz_cur - prev_peak.mz) / mz_cur * 1e6;
                    if mass_diff_with_sign.abs() <= htol_ppm && cur_intensity >= best_intensity {
                        best_intensity = cur_intensity;
                        best_prev = Some(prev);
                        best_mass_diff_with_sign = mass_diff_with_sign;
                    }
                }
            }
            if let Some(prev) = best_prev {
                banned_prev[prev] = true;
                if let Some(diffs) = mass_diffs.as_deref_mut() {
                    diffs.push(best_mass_diff_with_sign);
                }
                if prev < state.prev_peak_to_hill.len() {
                    assigned_hill = state.prev_peak_to_hill[prev];
                }
            }
        }

        let hill_id = if assigned_hill != usize::MAX && assigned_hill < hills.len() {
            assigned_hill
        } else {
            hills.push(Hill {
                hill_idx: hills.len(),
                ..Hill::default()
            });
            hills.len() - 1
        };
        let hill = &mut hills[hill_id];
        let intensity = f64::from(spectrum.peaks[idx].intensity);
        hill.scan_indices.push(scan_idx);
        hill.peak_indices.push(idx);
        hill.mz_values.push(spectrum.peaks[idx].mz);
        hill.intensities.push(intensity);
        hill.rt_values.push(rt);
        hill.drift_times.push(spectrum.drift_time);
        // The per-peak ion mobility is recorded whether or not gating is on;
        // a spectrum with only a drift time records -1 (the source's
        // `shouldThrowForMissingIM_` never throws).
        let ion_mobility = match im {
            Some(values) => f64::from(*values.get(idx).ok_or_else(|| im_too_short(idx))?),
            None => {
                let _ = should_throw_for_missing_im(spectrum);
                -1.0
            }
        };
        hill.ion_mobilities.push(ion_mobility);
        hill.length = hill.scan_indices.len();
        hill.rt_end = rt;
        if hill.length == 1 {
            hill.rt_start = rt;
            hill.rt_apex = rt;
            hill.intensity_apex = intensity;
        } else if intensity > hill.intensity_apex {
            hill.intensity_apex = intensity;
            hill.rt_apex = rt;
        }
        hill.intensity_sum += intensity;
        curr_peak_to_hill[idx] = hill_id;
    }

    state.prev_fast_dict = fast_dict;
    state.prev_peak_to_hill = curr_peak_to_hill;
    state.prev_im_bins = im_bin_per_peak;
    state.prev_spectrum = Some(spectrum);
    Ok(())
}

/// Whether a spectrum with a drift time but no ion-mobility array is an error,
/// the source `shouldThrowForMissingIM_`: never. *The Python reference
/// implementation gracefully degrades when ion mobility arrays are missing.*
pub fn should_throw_for_missing_im(spectrum: &MSSpectrum) -> bool {
    let _ = spectrum;
    false
}

/// Filter hills by length and compute their summary statistics, the source
/// `processHills_`.
///
/// Hills with at least `min_length` points are kept, with the median drift
/// time, the intensity-weighted mean m/z (left at its previous value when the
/// intensities sum to zero) and the intensity-weighted mean ion mobility
/// (the median when the intensities sum to zero).
///
/// # Errors
///
/// Only when a median's sort permutation cannot be allocated.
pub fn process_hills(hills: &[Hill], min_length: usize) -> Result<Vec<Hill>> {
    let mut processed = Vec::new();
    for hill in hills {
        if hill.length < min_length {
            continue;
        }
        let mut out = hill.clone();
        out.drift_time_median = calculate_median(&hill.drift_times)?;
        if !hill.mz_values.is_empty() && hill.mz_values.len() == hill.intensities.len() {
            let mut weighted_mz_sum = 0.0;
            let mut intensity_sum = 0.0;
            for (&mz, &intensity) in hill.mz_values.iter().zip(&hill.intensities) {
                weighted_mz_sum += mz * intensity;
                intensity_sum += intensity;
            }
            if intensity_sum > 0.0 {
                out.mz_weighted_mean = weighted_mz_sum / intensity_sum;
            }
        }
        if !hill.ion_mobilities.is_empty() && hill.ion_mobilities.len() == hill.intensities.len() {
            let mut weighted_im_sum = 0.0;
            let mut intensity_sum_im = 0.0;
            for (&im, &intensity) in hill.ion_mobilities.iter().zip(&hill.intensities) {
                weighted_im_sum += im * intensity;
                intensity_sum_im += intensity;
            }
            out.ion_mobility_median = if intensity_sum_im > 0.0 {
                weighted_im_sum / intensity_sum_im
            } else {
                calculate_median(&hill.ion_mobilities)?
            };
        } else {
            out.ion_mobility_median = calculate_median(&hill.ion_mobilities)?;
        }
        processed.push(out);
    }
    Ok(processed)
}

/// Split hills at valley positions to separate co-eluting species, the source
/// `splitHills_` (the Python `split_peaks`).
///
/// Only hills of at least `2 * min_length` points are examined. Their
/// intensities are smoothed ([`mean_filter`], half-width 3); a first pass
/// collects candidate valleys whose left and right maxima both exceed the
/// valley by the factor `hvf`, at least `min_length` points apart and from the
/// ends (keeping, among close candidates, the one with the larger product of
/// the two ratios); a second pass keeps those whose right maximum up to the next
/// candidate still exceeds the valley by `hvf`. The hill is cut there, and
/// segments of at least `min_length` points become hills with recomputed
/// bounds, apex, intensity sum and (for a positive sum) weighted m/z. The
/// first segment keeps the hill's index; later segments take fresh indices
/// above every input index, in order. The drift-time and ion-mobility medians
/// are copied from the original hill, not recomputed.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] with the source's text `Split hill meta data
/// arrays are inconsistent.` when a hill's drift-time or ion-mobility vector
/// is shorter than its scan list, and natively when a hill's per-point vectors
/// are shorter than its `length` (the source's `.at()` throws
/// `std::out_of_range` or indexes out of bounds), or when `min_length` is 0,
/// which the `minlh` parameter refuses and for which the source reads before
/// the smoothed trace.
pub fn split_hills(hills: &[Hill], hvf: f64, min_length: usize) -> Result<Vec<Hill>> {
    if min_length == 0 {
        return Err(Error::InvalidValue(
            "splitHills_ needs a minimum hill length of at least 1".into(),
        ));
    }
    let mut split = Vec::new();
    let mut next_hill_idx = 0usize;
    for hill in hills {
        next_hill_idx = next_hill_idx.max(hill.hill_idx.wrapping_add(1));
    }

    for hill in hills {
        if hill.length < min_length.saturating_mul(2) {
            split.push(hill.clone());
            continue;
        }
        let hill_len = hill.length;
        if hill.intensities.len() < hill_len
            || hill.scan_indices.len() < hill_len
            || hill.peak_indices.len() < hill_len
            || hill.mz_values.len() < hill_len
            || hill.rt_values.len() < hill_len
        {
            return Err(Error::InvalidValue(format!(
                "hill {} holds fewer points than its length {hill_len}",
                hill.hill_idx
            )));
        }
        let smoothed = mean_filter(&hill.intensities[..hill_len], 3);

        let mut min_idx_list: Vec<usize> = Vec::new();
        let mut recheck_positions: Vec<usize> = Vec::new();
        let c_len = hill_len - min_length;
        let mut idx = min_length - 1;
        let mut l_idx = 0usize;
        let mut min_val = 0.0;
        while idx <= c_len {
            if let Some(&last) = min_idx_list.last() {
                if idx >= last + min_length {
                    l_idx = last;
                }
            }
            let valley_intensity = smoothed[idx];
            if valley_intensity <= 0.0 {
                idx += 1;
                continue;
            }
            let mut left_max = 0.0;
            for &value in &smoothed[l_idx..idx] {
                left_max = std_max(left_max, value);
            }
            if left_max == 0.0 {
                idx += 1;
                continue;
            }
            let l_r = left_max / valley_intensity;
            if l_r >= hvf {
                let mut right_max = 0.0;
                for &value in &smoothed[idx + 1..hill_len] {
                    right_max = std_max(right_max, value);
                }
                if right_max > 0.0 {
                    let r_r = right_max / valley_intensity;
                    if r_r >= hvf {
                        let mult_val = l_r * r_r;
                        let candidate_pos = idx + usize::from(l_r > r_r);
                        if min_length <= candidate_pos && candidate_pos <= c_len {
                            let far_enough = min_idx_list
                                .last()
                                .is_none_or(|&last| candidate_pos >= last + min_length);
                            if far_enough {
                                min_idx_list.push(candidate_pos);
                                recheck_positions.push(idx);
                                min_val = mult_val;
                            } else if mult_val > min_val {
                                if let Some(last) = min_idx_list.last_mut() {
                                    *last = candidate_pos;
                                }
                                if let Some(last) = recheck_positions.last_mut() {
                                    *last = idx;
                                }
                                min_val = mult_val;
                            }
                        }
                    }
                }
            }
            idx += 1;
        }

        let mut final_splits = Vec::new();
        for (k, (&min_idx, &recheck_idx)) in min_idx_list.iter().zip(&recheck_positions).enumerate()
        {
            let end_idx = min_idx_list.get(k + 1).copied().unwrap_or(hill_len);
            if recheck_idx + 1 >= end_idx || recheck_idx >= hill_len {
                continue;
            }
            let recheck_val = smoothed[recheck_idx];
            if recheck_val <= 0.0 {
                continue;
            }
            let mut right_max = 0.0;
            for &value in &smoothed[recheck_idx + 1..end_idx] {
                right_max = std_max(right_max, value);
            }
            if right_max / recheck_val >= hvf {
                final_splits.push(min_idx);
            }
        }
        if final_splits.is_empty() {
            split.push(hill.clone());
            continue;
        }

        let mut boundaries = Vec::with_capacity(final_splits.len() + 2);
        boundaries.push(0);
        boundaries.extend_from_slice(&final_splits);
        boundaries.push(hill_len);
        for (s, bounds) in boundaries.windows(2).enumerate() {
            let (seg_start, seg_end) = (bounds[0], bounds[1]);
            if seg_end <= seg_start || seg_end - seg_start < min_length {
                continue;
            }
            let range = seg_start..seg_end;
            let drift_times = hill.drift_times.get(range.clone());
            let ion_mobilities = hill.ion_mobilities.get(range.clone());
            let (Some(drift_times), Some(ion_mobilities)) = (drift_times, ion_mobilities) else {
                return Err(source_invalid_value(
                    "Split hill meta data arrays are inconsistent.",
                    &(seg_end - seg_start).to_string(),
                ));
            };
            let mut segment = Hill {
                scan_indices: hill.scan_indices[range.clone()].to_vec(),
                peak_indices: hill.peak_indices[range.clone()].to_vec(),
                mz_values: hill.mz_values[range.clone()].to_vec(),
                intensities: hill.intensities[range.clone()].to_vec(),
                rt_values: hill.rt_values[range].to_vec(),
                drift_times: drift_times.to_vec(),
                ion_mobilities: ion_mobilities.to_vec(),
                ..hill.clone_summary()
            };
            segment.length = segment.scan_indices.len();
            segment.rt_start = segment.rt_values[0];
            segment.rt_end = segment.rt_values[segment.length - 1];
            let mut weighted_mz_sum = 0.0;
            let mut intensity_sum = 0.0;
            for (&mz, &intensity) in segment.mz_values.iter().zip(&segment.intensities) {
                weighted_mz_sum += mz * intensity;
                intensity_sum += intensity;
            }
            if intensity_sum > 0.0 {
                segment.mz_weighted_mean = weighted_mz_sum / intensity_sum;
            }
            segment.intensity_sum = intensity_sum;
            // `std::max_element`: the first maximum.
            let mut apex_idx = 0;
            for (i, &value) in segment.intensities.iter().enumerate() {
                if segment.intensities[apex_idx] < value {
                    apex_idx = i;
                }
            }
            segment.intensity_apex = segment.intensities[apex_idx];
            segment.rt_apex = segment.rt_values[apex_idx];
            if s == 0 {
                segment.hill_idx = hill.hill_idx;
            } else {
                segment.hill_idx = next_hill_idx;
                next_hill_idx = next_hill_idx.wrapping_add(1);
            }
            split.push(segment);
        }
    }
    Ok(split)
}

impl Hill {
    /// A copy of the summary fields with empty per-point vectors.
    fn clone_summary(&self) -> Self {
        Self {
            scan_indices: Vec::new(),
            peak_indices: Vec::new(),
            mz_values: Vec::new(),
            intensities: Vec::new(),
            rt_values: Vec::new(),
            drift_times: Vec::new(),
            ion_mobilities: Vec::new(),
            ..*self
        }
    }
}
