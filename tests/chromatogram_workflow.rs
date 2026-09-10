// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::analysis::peak_integrator::{IntegrationMethod, PeakIntegrator};
use openms::kernel::{ChromatogramPeak, DataArray, MSChromatogram};
use openms::processing::chromatogram::{
    ChromatogramPickingMethod, ChromatogramSmoothing, PeakPickerChromatogram,
};

fn trace() -> MSChromatogram {
    // Exact binary positions intentionally cannot be represented in the source
    // f32 boundary annotations. Literal intensities sum to 286; at spacing 1/8,
    // trapezoidal integration between zero endpoints is exactly 35.75.
    let signal = [
        0., 0., 0., 1., 2., 4., 8., 16., 32., 48., 64., 48., 32., 16., 8., 4., 2., 1., 0., 0., 0.,
    ];
    let mut input = MSChromatogram::from_peaks(
        signal
            .iter()
            .enumerate()
            .map(|(i, &y)| ChromatogramPeak::new(100_000_000.0 + i as f64 / 8.0, y))
            .collect(),
    );
    input.native_id = "transition=1".into();
    input.metadata.insert("sample".into(), "workflow".into());
    input.integer_data_arrays.push(DataArray::new(
        "original_sample",
        (0..signal.len() as i32).collect(),
    ));
    input
}

fn picker(method: ChromatogramPickingMethod) -> PeakPickerChromatogram {
    PeakPickerChromatogram {
        method,
        smoothing: ChromatogramSmoothing::Gaussian { width: 1.0 },
        signal_to_noise: 0.0,
        seed_signal_to_noise: 0.0,
        ..Default::default()
    }
}

#[test]
fn picked_regions_integrate_raw_samples_without_rounding_their_boundaries() {
    let input = trace();
    let original = input.clone();
    for method in [
        ChromatogramPickingMethod::Legacy,
        ChromatogramPickingMethod::Corrected,
    ] {
        let result = picker(method).pick_chromatogram(&input).unwrap();
        assert_eq!(result.regions.len(), 1);
        let region = &result.regions[0];
        assert_eq!(region.apex_index, 10);
        let left = input.peaks[region.left_index].rt;
        let right = input.peaks[region.right_index].rt;
        assert!(left < right);
        assert_eq!(left as f32, right as f32);
        let sum = PeakIntegrator::default()
            .integrate_chromatogram(&input, left, right)
            .unwrap();
        assert_eq!(sum.area, 286.0);
        assert_eq!(sum.height, 64.0);
        assert_eq!(sum.apex_pos, 100_000_001.25);
        assert_eq!(
            sum.hull_points.len(),
            region.right_index - region.left_index + 1
        );
        assert_eq!(sum.hull_points.first().unwrap()[0], left);
        assert_eq!(sum.hull_points.last().unwrap()[0], right);
        let integrator = PeakIntegrator {
            integration_method: IntegrationMethod::Trapezoid,
            ..Default::default()
        };
        let area = integrator
            .integrate_chromatogram(&input, left, right)
            .unwrap();
        assert_eq!(area.area, 35.75);
        let baseline = integrator
            .estimate_background_chromatogram(&input, left, right, sum.apex_pos)
            .unwrap();
        assert_eq!(baseline.area, 0.0);
        assert_eq!(baseline.height, 0.0);
        let shape = integrator
            .calculate_shape_metrics_chromatogram(&input, left, right, sum.height, sum.apex_pos)
            .unwrap();
        assert_eq!(shape.width_at_50, 0.5);
        assert_eq!(shape.points_across_half_height, 5);
        assert_eq!(shape.tailing_factor, 1.0);
        assert_eq!(shape.asymmetry_factor, 1.0);
        let picked = &result.picked.chromatogram;
        let abundance = picked
            .float_data_arrays
            .iter()
            .find(|array| array.name == "IntegratedIntensity")
            .unwrap();
        assert_eq!(abundance.data, [286.0]);
        assert_eq!(picked.metadata, input.metadata);
        assert_eq!(picked.native_id, input.native_id);
        assert_eq!(
            result.smoothed.integer_data_arrays,
            input.integer_data_arrays
        );
        assert_eq!(result.picked.omitted_arrays, ["original_sample"]);
        // Integrating the smoothed signal would change the reported raw sum.
        assert_ne!(result.smoothed.peaks, input.peaks);
    }
    assert_eq!(input, original);
}

#[cfg(feature = "mzml")]
#[test]
fn picked_chromatogram_annotations_survive_native_mzml_interchange() {
    let input = trace();
    let result = picker(ChromatogramPickingMethod::Corrected)
        .pick_chromatogram(&input)
        .unwrap();
    let experiment = openms::MSExperiment {
        chromatograms: vec![result.picked.chromatogram.clone()],
        ..Default::default()
    };
    let mut xml = Vec::new();
    openms::format::mzml::write(&mut xml, &experiment).unwrap();
    let read = openms::format::mzml::read(xml.as_slice()).unwrap();
    assert_eq!(read.chromatograms.len(), 1);
    assert_eq!(
        read.chromatograms[0].peaks,
        experiment.chromatograms[0].peaks
    );
    assert_eq!(
        read.chromatograms[0].float_data_arrays,
        experiment.chromatograms[0].float_data_arrays
    );
    assert_eq!(read.chromatograms[0].metadata, input.metadata);
}
