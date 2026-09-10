// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::analysis::emg::{EmgGradientDescent, EmgParameters};
use openms::analysis::peak_integrator::{IntegrationMethod, PeakIntegrator};
use openms::kernel::{ChromatogramPeak, MSChromatogram};

fn parameters() -> EmgParameters {
    EmgParameters {
        h: 1000.0,
        mu: 100.0,
        sigma: 1.0,
        tau: 2.0,
    }
}

fn sampled_trace(positions: &[f64]) -> MSChromatogram {
    let curve = EmgGradientDescent {
        compute_additional_points: false,
        ..Default::default()
    }
    .apply_parameters(positions, parameters())
    .unwrap();
    MSChromatogram::from_peaks(
        curve
            .positions
            .iter()
            .zip(&curve.intensities)
            .map(|(&x, &y)| ChromatogramPeak::new(x, y as f32))
            .collect(),
    )
}

#[test]
fn integrated_emg_obeys_independent_continuous_area_and_centroid_identities() {
    // Convolution by a normalized exponential preserves Gaussian area and
    // shifts its centroid by tau. Wide bounds leave less than 1e-9 of the area
    // outside the grid; f32 intensities and quadrature set the test tolerance.
    let positions: Vec<_> = (0..=6600).map(|i| 87.0 + f64::from(i) * 0.01).collect();
    let trace = sampled_trace(&positions);
    let integrated = PeakIntegrator {
        integration_method: IntegrationMethod::Trapezoid,
        ..Default::default()
    }
    .integrate_chromatogram(&trace, 87.0, 153.0)
    .unwrap();
    let p = parameters();
    let total = p.h * p.sigma * (2.0 * std::f64::consts::PI).sqrt();
    assert!((integrated.area / total - 1.0).abs() < 2e-7);
    let moment: f64 = trace
        .peaks
        .windows(2)
        .map(|pair| {
            (pair[1].rt - pair[0].rt)
                * (pair[0].rt * f64::from(pair[0].intensity)
                    + pair[1].rt * f64::from(pair[1].intensity))
                / 2.0
        })
        .sum();
    assert!((moment / integrated.area - (p.mu + p.tau)).abs() < 2e-5);
}

#[cfg(feature = "mzml")]
#[test]
fn cropped_fit_can_be_exchanged_and_integrated_without_unaligned_parameters() {
    use openms::kernel::DataArray;
    let positions: Vec<_> = (0..=35).map(|i| 97.0 + f64::from(i) * 0.2).collect();
    let mut input = sampled_trace(&positions);
    input.native_id = "synthetic_emg".into();
    input.metadata.insert("sample".into(), "cropped".into());
    input
        .integer_data_arrays
        .push(DataArray::new("sample_index", (0..36).collect()));
    let original = input.clone();
    let result = EmgGradientDescent::default()
        .fit_chromatogram(&input, None, None)
        .unwrap();
    assert_eq!(input, original);
    assert!(result.chromatogram.len() > input.len());
    assert_eq!(result.omitted_arrays, ["sample_index"]);
    assert!(result.chromatogram.float_data_arrays.is_empty());
    assert!(result.estimate.loss < 1.0);
    let experiment = openms::MSExperiment {
        chromatograms: vec![result.chromatogram.clone()],
        ..Default::default()
    };
    let mut xml = Vec::new();
    openms::format::mzml::write(&mut xml, &experiment).unwrap();
    let decoded = openms::format::mzml::read(xml.as_slice()).unwrap();
    let trace = &decoded.chromatograms[0];
    assert_eq!(trace.peaks, result.chromatogram.peaks);
    assert_eq!(trace.metadata, input.metadata);
    assert_eq!(trace.native_id, input.native_id);
    let integrator = PeakIntegrator {
        integration_method: IntegrationMethod::Trapezoid,
        ..Default::default()
    };
    let left = trace.peaks.first().unwrap().rt;
    let right = trace.peaks.last().unwrap().rt;
    let reconstructed = integrator
        .integrate_chromatogram(trace, left, right)
        .unwrap();
    let measured = integrator
        .integrate_chromatogram(&input, 97.0, 104.0)
        .unwrap();
    let p = parameters();
    let expected = p.h * p.sigma * (2.0 * std::f64::consts::PI).sqrt();
    assert!((reconstructed.area - expected).abs() < (measured.area - expected).abs());
    let automatic = PeakIntegrator {
        emg: Some(EmgGradientDescent::default()),
        ..integrator
    }
    .integrate_chromatogram(&input, 97.0, 104.0)
    .unwrap();
    assert_eq!(automatic, reconstructed);
}
