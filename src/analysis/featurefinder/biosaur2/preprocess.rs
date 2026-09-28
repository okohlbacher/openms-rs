// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Biosaur2 preprocessing: profile centroiding (`centroidProfileSpectra_`),
//! TOF intensity filtering (`processTOF_`) and PASEF/TIMS centroiding in joint
//! m/z-ion-mobility space (`centroidPASEFData_`).
//!
//! Each step reads the working spectra and writes the ones it changes into a
//! parallel `Option<MSSpectrum>` slot, so that the stored experiment is only
//! replaced once a run has succeeded (see
//! [`Biosaur2Algorithm::run_with_details`]). The source modifies the
//! experiment in place; the resulting spectra are the same.

use super::{Biosaur2Algorithm, Biosaur2Message, info};
use crate::Result;
use crate::concept::constants::user_param;
use crate::kernel::{DataArray, MSSpectrum, Peak1D, SpectrumType};
use crate::math::source_sort::source_sort_by;
use crate::metadata::{ImTypes, IonMobilityFormat};
use crate::processing::noise_estimation::x86::cvttsd2si32;
use crate::processing::peak_picking::{PeakPickerHiRes, PickingCompatibility};
use std::collections::BTreeMap;

/// The spectrum at working position `k`: its replacement if a step made one,
/// otherwise the stored spectrum.
fn current<'a>(
    spectra: &'a [MSSpectrum],
    ms1: &[usize],
    replaced: &'a [Option<MSSpectrum>],
    k: usize,
) -> &'a MSSpectrum {
    replaced[k].as_ref().unwrap_or(&spectra[ms1[k]])
}

/// Centroid the profile spectra with `PeakPickerHiRes`, the source
/// `centroidProfileSpectra_`.
///
/// The picker runs with its defaults and `signal_to_noise` 0, and with every
/// source compatibility behaviour ([`PickingCompatibility::source`]), because
/// the source picker accepts what the native default refuses. A spectrum stored
/// as centroided is kept as it is; every other one is replaced by the picked
/// spectrum, with the input's retention time, MS level and (when set) drift
/// time, marked centroided. The source rebuilds the experiment from the
/// spectra only, so chromatograms and settings survive, as they do here.
///
/// # Errors
///
/// Whatever [`PeakPickerHiRes::pick_spectrum`] refuses.
pub(crate) fn centroid_profile_spectra(
    spectra: &[MSSpectrum],
    ms1: &[usize],
    replaced: &mut [Option<MSSpectrum>],
    messages: &mut Vec<Biosaur2Message>,
) -> Result<()> {
    info(
        messages,
        "Centroiding profile spectra using PeakPickerHiRes...".to_owned(),
    );
    let picker = PeakPickerHiRes {
        signal_to_noise: 0.0,
        compatibility: PickingCompatibility::source(),
        ..PeakPickerHiRes::default()
    };
    let mut before = 0usize;
    let mut after = 0usize;
    let mut picked = Vec::with_capacity(ms1.len());
    for k in 0..ms1.len() {
        let input = current(spectra, ms1, replaced, k);
        before += input.len();
        if input.spectrum_type == SpectrumType::Centroid {
            after += input.len();
            picked.push(None);
            continue;
        }
        let mut output = picker.pick_spectrum(input)?.spectrum;
        output.rt = input.rt;
        output.ms_level = input.ms_level;
        output.spectrum_type = SpectrumType::Centroid;
        if input.drift_time >= 0.0 {
            output.drift_time = input.drift_time;
        }
        after += output.len();
        picked.push(Some(output));
    }
    for (slot, output) in replaced.iter_mut().zip(picked) {
        if output.is_some() {
            *slot = output;
        }
    }
    info(
        messages,
        format!("Centroiding: {before} profile points -> {after} centroided peaks"),
    );
    Ok(())
}

/// `mean` and population standard deviation of `values`, summed left to right
/// as `std::accumulate` and the source's loop do.
fn mean_and_std(values: &[f64]) -> (f64, f64) {
    let sum = values.iter().fold(0.0, |acc, &v| acc + v);
    let mean = sum / values.len() as f64;
    let mut sq_sum = 0.0;
    for &value in values {
        sq_sum += (value - mean) * (value - mean);
    }
    (mean, (sq_sum / values.len() as f64).sqrt())
}

/// The 50 Th bin of an m/z value: `static_cast<int>(mz / 50.0)`.
fn tof_bin(mz: f64) -> i32 {
    cvttsd2si32(mz / 50.0)
}

/// Apply TOF-specific intensity filtering, the source `processTOF_`.
///
/// Per 50 Th m/z bin, a threshold is learned from the log10 intensities of the
/// first 25 spectra (mean + 2 standard deviations over bins with at least 150
/// values, the first pass of the Python reference) and then refined on each of
/// the first 50 spectra (mean + 3 standard deviations where that spectrum has
/// at least 150 values in the bin; the last such spectrum wins). Every peak
/// whose intensity is not strictly above its bin's threshold (150 for a bin
/// without one) is removed, with its data-array entries. Unlike the original
/// C++ code, the training range is not restricted to `[minmz, maxmz]`, to match
/// the Python implementation.
///
/// `log10` and `pow` are the `libm` crate's, so the thresholds are the same on
/// every platform; the Release build calls glibc's. The differential test
/// compares the filtered spectra bit for bit on its TOF input.
///
/// # Errors
///
/// Whatever [`MSSpectrum::select`] refuses: a non-empty data array whose length
/// differs from the peak count (the source's `checkDataArraySizes_`).
pub(crate) fn process_tof(
    spectra: &[MSSpectrum],
    ms1: &[usize],
    replaced: &mut [Option<MSSpectrum>],
    messages: &mut Vec<Biosaur2Message>,
) -> Result<()> {
    info(
        messages,
        "Applying TOF-specific intensity filtering...".to_owned(),
    );
    let count = ms1.len();
    let any_im = (0..count).any(|k| {
        ImTypes::determine_im_format(current(spectra, ms1, replaced, k)) != IonMobilityFormat::None
    });
    if any_im {
        info(
            messages,
            "TOF filtering: input data contains ion-mobility arrays.".to_owned(),
        );
    }

    let log_bins = |spectrum: &MSSpectrum| {
        let mut bins: BTreeMap<i32, Vec<f64>> = BTreeMap::new();
        for peak in &spectrum.peaks {
            let bin = tof_bin(peak.mz);
            let intensity = f64::from(peak.intensity);
            if intensity <= 0.0 {
                continue;
            }
            let log_intensity = libm::log10(intensity);
            if !log_intensity.is_finite() {
                continue;
            }
            bins.entry(bin).or_default().push(log_intensity);
        }
        bins
    };

    // Phase 1: initial thresholds from the first 25 spectra, mean + 2 sigma.
    let mut intensity_bins: BTreeMap<i32, Vec<f64>> = BTreeMap::new();
    for k in 0..count.min(25) {
        for (bin, mut values) in log_bins(current(spectra, ms1, replaced, k)) {
            intensity_bins.entry(bin).or_default().append(&mut values);
        }
    }
    let mut thresholds: BTreeMap<i32, f64> = BTreeMap::new();
    for (&bin, values) in &intensity_bins {
        // Every stored value is finite already; the source filters again.
        if values.len() >= 150 {
            let (mean, std_dev) = mean_and_std(values);
            if !mean.is_finite() || !std_dev.is_finite() {
                continue;
            }
            thresholds.insert(bin, libm::pow(10.0, mean + 2.0 * std_dev));
        }
    }
    // Phase 2: refinement on each of the first 50 spectra, mean + 3 sigma.
    for k in 0..count.min(50) {
        for (bin, values) in log_bins(current(spectra, ms1, replaced, k)) {
            if values.len() >= 150 {
                let (mean, std_dev) = mean_and_std(&values);
                if !mean.is_finite() || !std_dev.is_finite() {
                    continue;
                }
                thresholds.insert(bin, libm::pow(10.0, mean + 3.0 * std_dev));
            }
        }
    }

    let mut before = 0usize;
    let mut after = 0usize;
    let mut filtered = Vec::with_capacity(count);
    for k in 0..count {
        let spectrum = current(spectra, ms1, replaced, k);
        before += spectrum.len();
        let keep: Vec<usize> = spectrum
            .peaks
            .iter()
            .enumerate()
            .filter(|(_, peak)| {
                let threshold = thresholds.get(&tof_bin(peak.mz)).copied().unwrap_or(150.0);
                f64::from(peak.intensity) > threshold
            })
            .map(|(index, _)| index)
            .collect();
        let mut output = spectrum.clone();
        output.select(&keep)?;
        after += output.len();
        filtered.push(output);
    }
    for (slot, output) in replaced.iter_mut().zip(filtered) {
        *slot = Some(output);
    }
    info(
        messages,
        format!("TOF filtering: {before} peaks -> {after} peaks"),
    );
    Ok(())
}

/// Centroid PASEF/TIMS spectra in joint m/z-ion mobility space, the source
/// `centroidPASEFData_`.
///
/// The source documentation: *Performs 2D clustering of peaks in the m/z and
/// ion mobility dimensions to reduce data complexity for ion mobility-enabled
/// instruments.* It mirrors the Python `centroid_pasef_scan`: the peaks inside
/// the intensity and m/z window are sorted by their coarse m/z bin (`mz_step`
/// wide, the source's `std::sort`), and each unclaimed peak starts a cluster
/// that chains every later peak within `htol` ppm of the current one, at most
/// one bin further in m/z and one `max(IM) * pasef_tolerance` bin in ion
/// mobility, and within `pasef_tolerance` in ion mobility. A cluster of at
/// least `pasefminlh` points whose summed intensity reaches `pasefmini`
/// becomes one peak at the intensity-weighted m/z and ion mobility; the
/// spectrum's data arrays are replaced by the one ion-mobility array
/// (`Ion Mobility`).
///
/// Spectra of another MS level, without a per-peak ion-mobility array, with an
/// array of another length than the peak list, without peaks, or whose largest
/// ion mobility is not positive are left unchanged. A spectrum none of whose
/// peaks lies in the window loses its peaks and data arrays. Nothing happens
/// when `mz_step` or `pasef_tolerance` is not positive.
///
/// `spectra` are one group's working spectra; a changed spectrum is written to
/// the matching slot of `owned`.
///
/// # Errors
///
/// Only an introsort refusal of [`source_sort_by`], which `<` on integers
/// cannot provoke.
pub fn centroid_pasef_data(
    algorithm: &Biosaur2Algorithm,
    spectra: &[&MSSpectrum],
    owned: &mut [Option<MSSpectrum>],
    mz_step: f64,
    pasef_tolerance: f64,
    messages: &mut Vec<Biosaur2Message>,
) -> Result<()> {
    if mz_step <= 0.0 || pasef_tolerance <= 0.0 {
        return Ok(());
    }
    let hill_mz_accuracy = algorithm.htol;
    let ion_mobility_accuracy = pasef_tolerance;
    let mut before = 0usize;
    let mut after = 0usize;

    for (k, spectrum) in spectra.iter().enumerate() {
        if spectrum.ms_level != 1 {
            continue;
        }
        before += spectrum.len();
        if ImTypes::determine_im_format(spectrum) != IonMobilityFormat::PerPeak {
            continue;
        }
        let (im_index, _) = spectrum.im_data()?;
        let im_array = &spectrum.float_data_arrays[im_index].data;
        if im_array.len() != spectrum.len() || spectrum.is_empty() {
            continue;
        }

        let mut mz_ar = Vec::new();
        let mut intensity_ar = Vec::new();
        let mut im_ar = Vec::new();
        for (peak, &im) in spectrum.peaks.iter().zip(im_array) {
            let intensity = f64::from(peak.intensity);
            if intensity < algorithm.mini || peak.mz < algorithm.minmz || peak.mz > algorithm.maxmz
            {
                continue;
            }
            mz_ar.push(peak.mz);
            intensity_ar.push(intensity);
            im_ar.push(f64::from(im));
        }
        if mz_ar.is_empty() {
            let mut cleared = (*spectrum).clone();
            cleared.clear(false);
            owned[k] = Some(cleared);
            continue;
        }
        // `std::max_element`: the first maximum under `<`.
        let mut max_im = im_ar[0];
        for &im in &im_ar[1..] {
            if max_im < im {
                max_im = im;
            }
        }
        if max_im <= 0.0 {
            continue;
        }
        let ion_mobility_step = max_im * ion_mobility_accuracy;
        if ion_mobility_step <= 0.0 {
            continue;
        }

        let n = mz_ar.len();
        let mz_fast: Vec<i32> = mz_ar.iter().map(|&mz| cvttsd2si32(mz / mz_step)).collect();
        let im_fast: Vec<i32> = im_ar
            .iter()
            .map(|&im| cvttsd2si32(im / ion_mobility_step))
            .collect();
        let mut order: Vec<usize> = (0..n).collect();
        source_sort_by(&mut order, |&a, &b| mz_fast[a] < mz_fast[b])?;
        let mz_sorted: Vec<f64> = order.iter().map(|&i| mz_ar[i]).collect();
        let intensity_sorted: Vec<f64> = order.iter().map(|&i| intensity_ar[i]).collect();
        let im_sorted: Vec<f64> = order.iter().map(|&i| im_ar[i]).collect();
        let mz_fast_sorted: Vec<i32> = order.iter().map(|&i| mz_fast[i]).collect();
        let im_fast_sorted: Vec<i32> = order.iter().map(|&i| im_fast[i]).collect();

        let mut banned = vec![false; n];
        let mut mz_new = Vec::new();
        let mut intensity_new = Vec::new();
        let mut im_new = Vec::new();
        let mut peak_idx = 0usize;
        while peak_idx < n {
            let mut tmp = Vec::new();
            if !banned[peak_idx] {
                let mass_accuracy_cur = mz_sorted[peak_idx] * 1e-6 * hill_mz_accuracy;
                let mz_val_int = mz_fast_sorted[peak_idx];
                let im_val_int = im_fast_sorted[peak_idx];
                tmp.push(peak_idx);
                let mut peak_idx_2 = peak_idx + 1;
                while peak_idx_2 < n {
                    if !banned[peak_idx_2] {
                        let mz_val_int_2 = mz_fast_sorted[peak_idx_2];
                        if i64::from(mz_val_int_2) - i64::from(mz_val_int) > 1 {
                            break;
                        }
                        if (mz_sorted[peak_idx] - mz_sorted[peak_idx_2]).abs() <= mass_accuracy_cur
                        {
                            let im_val_int_2 = im_fast_sorted[peak_idx_2];
                            if (i64::from(im_val_int) - i64::from(im_val_int_2)).abs() <= 1
                                && (im_sorted[peak_idx] - im_sorted[peak_idx_2]).abs()
                                    <= ion_mobility_accuracy
                            {
                                // The chain continues from the peak just
                                // added, as the Python reference does.
                                tmp.push(peak_idx_2);
                                peak_idx = peak_idx_2;
                            }
                        }
                    }
                    peak_idx_2 += 1;
                }
            }

            let l_new = tmp.len();
            if l_new >= algorithm.pasefminlh {
                if l_new == 1 {
                    let i_val_new = intensity_sorted[peak_idx];
                    if i_val_new >= algorithm.pasefmini {
                        mz_new.push(mz_sorted[peak_idx]);
                        intensity_new.push(i_val_new);
                        im_new.push(im_sorted[peak_idx]);
                        banned[peak_idx] = true;
                    }
                } else {
                    let mut i_val_new = 0.0;
                    for &idx in &tmp {
                        i_val_new += intensity_sorted[idx];
                    }
                    if i_val_new >= algorithm.pasefmini {
                        let mut mz_weighted_sum = 0.0;
                        let mut im_weighted_sum = 0.0;
                        for &idx in &tmp {
                            mz_weighted_sum += mz_sorted[idx] * intensity_sorted[idx];
                            im_weighted_sum += im_sorted[idx] * intensity_sorted[idx];
                        }
                        mz_new.push(mz_weighted_sum / i_val_new);
                        intensity_new.push(i_val_new);
                        im_new.push(im_weighted_sum / i_val_new);
                        for &idx in &tmp {
                            banned[idx] = true;
                        }
                    }
                }
            }
            peak_idx += 1;
        }

        let mut centroided = (*spectrum).clone();
        centroided.clear(false);
        centroided.peaks = mz_new
            .iter()
            .zip(&intensity_new)
            .map(|(&mz, &intensity)| Peak1D::new(mz, intensity as f32))
            .collect();
        centroided.float_data_arrays = vec![DataArray::new(
            user_param::ION_MOBILITY,
            im_new.iter().map(|&im| im as f32).collect(),
        )];
        after += centroided.len();
        owned[k] = Some(centroided);
    }
    info(
        messages,
        format!("PASEF centroiding: {before} peaks -> {after} centroided clusters"),
    );
    Ok(())
}
