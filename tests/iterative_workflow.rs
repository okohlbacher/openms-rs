// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::kernel::{DataArray, MSSpectrum, Peak1D, SpectrumType};
use openms::processing::SpectrumFilter;
use openms::processing::iterative::PeakPickerIterative;
use openms::processing::mean_noise::SignalToNoiseEstimatorMeanIterative;
use openms::processing::window_mower::WindowMower;

fn trace() -> MSSpectrum {
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
    input.native_id = "scan=1".into();
    input
        .metadata
        .insert("sample".into(), "iterative workflow".into());
    let noise = SignalToNoiseEstimatorMeanIterative {
        window_length: 2.0,
        ..Default::default()
    }
    .estimate_spectrum(&input)
    .unwrap();
    assert_eq!(noise.sparse_window_percent, 0.0);
    assert_eq!(noise.noise.len(), input.len());
    assert!(noise.noise.iter().all(|v| v.is_finite() && *v >= 1.0));
    input.float_data_arrays.push(DataArray::new(
        "mean_signal_to_noise",
        noise.signal_to_noise.iter().map(|&v| v as f32).collect(),
    ));
    input
}

#[test]
fn mean_noise_iterative_picking_and_window_selection_preserve_scientific_associations() {
    let input = trace();
    let original = input.clone();
    let result = PeakPickerIterative::default()
        .pick_spectrum(&input)
        .unwrap();
    assert_eq!(result.regions.len(), 4);
    assert_eq!(result.picked.omitted_arrays, ["mean_signal_to_noise"]);
    assert_eq!(result.picked.spectrum.spectrum_type, SpectrumType::Centroid);
    assert_eq!(input, original);
    for ((peak, region), center) in result
        .picked
        .spectrum
        .peaks
        .iter()
        .zip(&result.regions)
        .zip([100.5, 101.5, 102.5, 103.5])
    {
        assert!((peak.mz - center).abs() < 0.00002);
        // Check the reported abundance independently against original samples.
        let raw_sum: f64 = input.peaks[region.left_index..=region.right_index]
            .iter()
            .map(|p| f64::from(p.intensity))
            .sum();
        assert_eq!(peak.intensity, raw_sum as f32);
    }
    let before = result.picked.spectrum;
    let mut selected = before.clone();
    WindowMower::default()
        .filter_spectrum(&mut selected)
        .unwrap();
    assert_eq!(selected.peaks, &before.peaks[..2]);
    assert_eq!(selected.native_id, input.native_id);
    assert_eq!(selected.metadata, input.metadata);
    assert_eq!(selected.float_data_arrays.len(), 3);
    for (retained, prior) in selected
        .float_data_arrays
        .iter()
        .zip(&before.float_data_arrays)
    {
        assert_eq!(retained.name, prior.name);
        assert_eq!(retained.data, &prior.data[..2]);
    }
}

#[cfg(feature = "mzml")]
#[test]
fn selected_iterative_centroids_and_width_arrays_round_trip_through_mzml() {
    let result = PeakPickerIterative::default()
        .pick_spectrum(&trace())
        .unwrap();
    let mut selected = result.picked.spectrum;
    WindowMower::default()
        .filter_spectrum(&mut selected)
        .unwrap();
    let experiment = openms::MSExperiment {
        spectra: vec![selected.clone()],
        ..Default::default()
    };
    let mut xml = Vec::new();
    openms::format::mzml::write(&mut xml, &experiment).unwrap();
    let read = openms::format::mzml::read(xml.as_slice()).unwrap();
    assert_eq!(read.spectra.len(), 1);
    assert_eq!(read.spectra[0].peaks, selected.peaks);
    assert_eq!(
        read.spectra[0].float_data_arrays,
        selected.float_data_arrays
    );
    assert_eq!(read.spectra[0].spectrum_type, SpectrumType::Centroid);
    assert_eq!(read.spectra[0].native_id, selected.native_id);
    assert_eq!(read.spectra[0].metadata, selected.metadata);
}
