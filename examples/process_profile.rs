// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Estimate profile noise, pick iterative centroids, and retain local top peaks.

use openms::Result;
use openms::kernel::{DataArray, MSSpectrum, Peak1D, SpectrumType};
use openms::processing::SpectrumFilter;
use openms::processing::iterative::PeakPickerIterative;
use openms::processing::mean_noise::SignalToNoiseEstimatorMeanIterative;
use openms::processing::window_mower::WindowMower;
use std::io::{self, Write};

fn main() -> Result<()> {
    let mut input = MSSpectrum::from_peaks(
        (0..=400)
            .map(|i| {
                let mz = 100.0 + f64::from(i) * 0.01;
                let intensity = [(100.5, 100.0), (101.5, 80.0), (102.5, 60.0), (103.5, 20.0)]
                    .iter()
                    .fold(1.0, |sum, &(center, height)| {
                        sum + height * (-0.5 * ((mz - center) / 0.05).powi(2)).exp()
                    });
                Peak1D::new(mz, intensity as f32)
            })
            .collect(),
    );
    input.spectrum_type = SpectrumType::Profile;
    let noise = SignalToNoiseEstimatorMeanIterative {
        window_length: 2.0,
        ..Default::default()
    }
    .estimate_spectrum(&input)?;
    // Noise describes observed samples. The picker reports this array as omitted,
    // since centroiding has no source aggregation rule for such annotations.
    input.float_data_arrays.push(DataArray::new(
        "mean_signal_to_noise",
        noise.signal_to_noise.iter().map(|&v| v as f32).collect(),
    ));
    let result = PeakPickerIterative::default().pick_spectrum(&input)?;
    let count = result.picked.spectrum.len();
    let mut selected = result.picked.spectrum;
    WindowMower::default().filter_spectrum(&mut selected)?;
    let mut out = io::BufWriter::new(io::stdout().lock());
    writeln!(
        out,
        "profile_samples\titerative_centroids\tretained_centroids\tsparse_noise_windows_percent"
    )?;
    writeln!(
        out,
        "{}\t{count}\t{}\t{:.6}",
        input.len(),
        selected.len(),
        noise.sparse_window_percent
    )?;
    writeln!(out, "mz\tintegrated_intensity\tleft_mz\tright_mz")?;
    for (i, p) in selected.peaks.iter().enumerate() {
        writeln!(
            out,
            "{:.6}\t{:.6}\t{:.6}\t{:.6}",
            p.mz,
            p.intensity,
            selected.float_data_arrays[1].data[i],
            selected.float_data_arrays[2].data[i]
        )?;
    }
    out.flush()?;
    Ok(())
}
