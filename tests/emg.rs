// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::analysis::emg::{EmgGradientDescent, EmgParameters};
use openms::kernel::{ChromatogramPeak, DataArray, MSChromatogram, MSSpectrum, Peak1D};

// Literal cutoff example from pinned EmgGradientDescent_test.cpp. The source
// fits f32 container intensities; its slice-level operations retain f64 values.
const X: [f64; 12] = [
    15.34253311,
    15.35624981,
    15.36995029,
    15.38366699,
    15.39736652,
    15.41156673,
    15.42574978,
    15.44018364,
    15.45436668,
    15.46856689,
    15.48274994,
    15.49695015,
];
const Y: [f64; 12] = [
    3.48297429,
    15.54384613,
    50.31319046,
    151.8971405,
    411.25631714,
    946.44311523,
    1642.56152344,
    2118.89526367,
    2055.13647461,
    1665.13232422,
    1275.53015137,
    1009.70056152,
];
fn chromatogram() -> MSChromatogram {
    MSChromatogram {
        peaks: X
            .iter()
            .zip(Y)
            .map(|(&x, y)| ChromatogramPeak::new(x, y as f32))
            .collect(),
        ..Default::default()
    }
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}
fn parameters() -> EmgParameters {
    EmgParameters {
        h: 15_515_900.0,
        mu: 14.3453,
        sigma: 0.0344277,
        tau: 0.188507,
    }
}

#[test]
fn source_cutoff_fit_and_container_equivalence() {
    let input = chromatogram();
    let fitter = EmgGradientDescent::default();
    let fit = fitter.fit_chromatogram(&input, None, None).unwrap();
    assert_eq!(fit.chromatogram.len(), 28);
    let p = fit.estimate.parameters;
    close(p.h, 3791.07, 1e-5);
    close(p.mu, 15.4227, 1e-5);
    close(p.sigma, 0.0210588, 1e-5);
    close(p.tau, 0.0476741, 1e-5);
    assert!(fit.estimate.best_iteration <= fit.estimate.iterations);
    assert_eq!(fit.estimate.training_points, 12);
    assert_eq!(fit.estimate.evaluations, fit.estimate.iterations * 12 * 5);
    let spectrum = MSSpectrum {
        peaks: input
            .peaks
            .iter()
            .map(|p| Peak1D::new(p.rt, p.intensity))
            .collect(),
        ..Default::default()
    };
    let spec_fit = fitter.fit_spectrum(&spectrum, None, None).unwrap();
    assert_eq!(spec_fit.estimate, fit.estimate);
    assert!(
        spec_fit
            .spectrum
            .peaks
            .iter()
            .zip(fit.chromatogram.peaks)
            .all(|(a, b)| a.mz == b.rt && a.intensity == b.intensity)
    );
}

#[test]
fn source_model_values_and_exact_sample_application() {
    let fitter = EmgGradientDescent {
        compute_additional_points: false,
        ..Default::default()
    };
    let p = parameters();
    let xs = [-3_333_333.0, p.mu - 1.0 / 60.0, p.mu + 1.0 / 60.0];
    let curve = fitter.apply_parameters(&xs, p).unwrap();
    assert_eq!(curve.positions, xs);
    assert_eq!(curve.intensities[0], 0.0);
    close(curve.intensities[1], 1992032.65711041, 1e-13);
    close(curve.intensities[2], 4088964.97520213, 1e-13);
    let p = EmgParameters {
        mu: 860.719,
        sigma: 2.06566,
        tau: 11.3104,
        ..p
    };
    let curve = fitter
        .apply_parameters(&[p.mu - 1.0, p.mu + 1.0], p)
        .unwrap();
    close(curve.intensities[0], 1992033.06584247, 1e-13);
    close(curve.intensities[1], 4088968.52957875, 1e-13);
    assert!(
        fitter
            .apply_parameters(&[], p)
            .unwrap()
            .positions
            .is_empty()
    );
    assert_eq!(
        fitter.apply_parameters(&[p.mu], p).unwrap().positions.len(),
        1
    );
}

#[test]
fn additional_points_keep_original_samples_and_last_threshold_crossing() {
    let p = EmgParameters {
        h: 100.0,
        mu: 10.0,
        sigma: 1.0,
        tau: 1.0,
    };
    let original = [7.0, 8.0, 9.0, 10.0, 11.0, 12.0];
    let fitter = EmgGradientDescent::default();
    let curve = fitter.apply_parameters(&original, p).unwrap();
    assert_eq!(&curve.positions[..original.len()], &original);
    assert!(curve.positions.len() > original.len());
    let n = curve.intensities.len();
    assert!(curve.intensities[n - 1] <= curve.intensities[0]);
    assert!(curve.intensities[n - 2] > curve.intensities[0]);
    assert!(curve.positions.windows(2).all(|w| w[1] - w[0] == 1.0));

    let low = fitter
        .apply_parameters(&original, EmgParameters { h: 0.01, ..p })
        .unwrap();
    let last = low.intensities.len() - 1;
    assert!(low.intensities[last] <= 1e-3);
    assert!(low.intensities[last - 1] > 1e-3);
    assert!(low.intensities[last] > low.intensities[0]);
}

#[test]
fn convergence_keeps_first_best_parameters_and_counts_evaluated_iterations() {
    let xs = [0.01, 0.0101, 0.0102];
    let fitter = EmgGradientDescent {
        max_iterations: 50,
        ..Default::default()
    };
    let fit = fitter.estimate_parameters(&xs, &[0.0; 3]).unwrap();
    assert!(fit.converged);
    assert_eq!(fit.iterations, 50);
    assert_eq!(fit.best_iteration, 1);
    assert_eq!(fit.training_points, 2);
    assert_eq!(fit.evaluations, 500);
    assert_eq!(fit.loss, 0.0);
    close(fit.parameters.mu, 0.0101, 1e-15);
    close(fit.parameters.sigma, 0.000101, 1e-15);
    let fit = EmgGradientDescent {
        max_iterations: 49,
        ..fitter
    }
    .estimate_parameters(&xs, &[0.0; 3])
    .unwrap();
    assert_eq!(fit.iterations, 49);
    assert!(!fit.converged);
}

#[test]
fn fitted_copy_preserves_metadata_and_reports_all_omitted_arrays() {
    let mut input = chromatogram();
    input.name = "cutoff reference".into();
    input.metadata.insert("origin".into(), "EMG fixture".into());
    input
        .float_data_arrays
        .push(DataArray::new("profile quality", vec![1.0; 12]));
    input
        .integer_data_arrays
        .push(DataArray::new("source index", (0..12).collect()));
    input
        .string_data_arrays
        .push(DataArray::new("label", vec!["raw".to_owned(); 12]));
    let before = input.clone();
    let output = EmgGradientDescent::default()
        .fit_chromatogram(&input, None, None)
        .unwrap();
    assert_eq!(input, before);
    assert_eq!(output.chromatogram.name, input.name);
    assert_eq!(output.chromatogram.metadata, input.metadata);
    assert_eq!(
        output.omitted_arrays,
        ["profile quality", "source index", "label"]
    );
    assert!(output.chromatogram.float_data_arrays.is_empty());
    assert!(output.chromatogram.integer_data_arrays.is_empty());
    assert!(output.chromatogram.string_data_arrays.is_empty());
    output.chromatogram.validate().unwrap();
}

#[test]
fn optional_bounds_are_inclusive_and_zero_is_not_a_sentinel() {
    let fitter = EmgGradientDescent {
        max_iterations: 1,
        compute_additional_points: false,
        ..Default::default()
    };
    let input = chromatogram();
    let selected = fitter
        .fit_chromatogram(&input, Some(X[2]), Some(X[9]))
        .unwrap();
    assert_eq!(
        selected
            .chromatogram
            .peaks
            .iter()
            .map(|p| p.rt)
            .collect::<Vec<_>>(),
        X[2..=9]
    );
    assert_eq!(selected.estimate.iterations, 1);
    assert!(!selected.estimate.converged);
    assert!(fitter.fit_chromatogram(&input, None, Some(0.0)).is_err());
    assert!(
        fitter
            .fit_chromatogram(&input, Some(X[4]), Some(X[3]))
            .is_err()
    );
    assert!(
        fitter
            .fit_chromatogram(&input, Some(f64::NAN), None)
            .is_err()
    );
}

#[test]
fn evaluation_limit_covers_estimation_and_application_together() {
    let input = chromatogram();
    let fitter = EmgGradientDescent {
        max_iterations: 1,
        compute_additional_points: false,
        max_evaluations: 60,
        ..Default::default()
    };
    assert_eq!(fitter.estimate_parameters(&X, &Y).unwrap().evaluations, 60);
    assert!(fitter.fit_chromatogram(&input, None, None).is_err());
    let fitter = EmgGradientDescent {
        max_evaluations: 72,
        ..fitter
    };
    assert!(fitter.fit_chromatogram(&input, None, None).is_ok());
    let fitter = EmgGradientDescent {
        max_evaluations: 11,
        ..fitter
    };
    assert!(fitter.apply_parameters(&X, parameters()).is_err());
    let fitter = EmgGradientDescent {
        max_points: 12,
        ..Default::default()
    };
    assert!(fitter.fit_chromatogram(&input, None, None).is_err());
    assert_eq!(input, chromatogram());
}

#[test]
fn invalid_and_nonfinite_inputs_fail_without_partial_results() {
    let fitter = EmgGradientDescent::default();
    for (x, y) in [
        (&[][..], &[][..]),
        (&[1.0][..], &[1.0][..]),
        (&[1.0, 2.0][..], &[1.0][..]),
        (&[1.0, 1.0][..], &[1.0, 2.0][..]),
        (&[2.0, 1.0][..], &[1.0, 2.0][..]),
        (&[1.0, f64::INFINITY][..], &[1.0, 2.0][..]),
        (&[1.0, 2.0][..], &[1.0, f64::NAN][..]),
        (&[-2.0, 2.0][..], &[1.0, 1.0][..]),
    ] {
        assert!(fitter.estimate_parameters(x, y).is_err());
    }
    for fitter in [
        EmgGradientDescent {
            max_iterations: 0,
            ..Default::default()
        },
        EmgGradientDescent {
            max_points: 0,
            ..Default::default()
        },
        EmgGradientDescent {
            max_evaluations: 0,
            ..Default::default()
        },
    ] {
        assert!(fitter.estimate_parameters(&X, &Y).is_err());
    }
    for p in [
        EmgParameters {
            sigma: 0.0,
            ..parameters()
        },
        EmgParameters {
            tau: -1.0,
            ..parameters()
        },
        EmgParameters {
            h: f64::INFINITY,
            ..parameters()
        },
    ] {
        assert!(fitter.apply_parameters(&X, p).is_err());
    }
    let fitter = EmgGradientDescent {
        compute_additional_points: false,
        ..fitter
    };
    // Raw exp(z²)*erfc(z) overflows in the source at ordinary tail positions.
    assert!(
        fitter
            .apply_parameters(
                &[-40.0],
                EmgParameters {
                    h: 1.0,
                    mu: 0.0,
                    sigma: 1.0,
                    tau: 1.0
                }
            )
            .is_err()
    );
    let mut input = chromatogram();
    input
        .float_data_arrays
        .push(DataArray::new("misaligned", vec![1.0]));
    let before = input.clone();
    assert!(fitter.fit_chromatogram(&input, None, None).is_err());
    assert_eq!(input, before);
}
