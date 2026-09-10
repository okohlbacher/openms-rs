// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Independent literal C++ test oracles; see data/emg_provenance.json.

use openms::analysis::emg::{EmgGradientDescent, EmgParameters};
use openms::kernel::{ChromatogramPeak, MSChromatogram, MSSpectrum, Peak1D};

struct Trace {
    positions: Vec<f64>,
    raw: Vec<f64>,
    stored: Vec<f32>,
}
fn trace(name: &str) -> Trace {
    let mut result = Trace {
        positions: vec![],
        raw: vec![],
        stored: vec![],
    };
    for line in include_str!("data/emg_traces.tsv").lines().skip(1) {
        let v: Vec<_> = line.split('\t').collect();
        if v[0] != name {
            continue;
        }
        assert_eq!(v[1].parse::<usize>().unwrap(), result.positions.len());
        let x: f64 = v[2].parse().unwrap();
        let y: f64 = v[3].parse().unwrap();
        assert_eq!(x.to_bits(), u64::from_str_radix(v[4], 16).unwrap());
        assert_eq!(y.to_bits(), u64::from_str_radix(v[5], 16).unwrap());
        let stored = f32::from_bits(u32::from_str_radix(v[6], 16).unwrap());
        assert_eq!(stored, y as f32);
        result.positions.push(x);
        result.raw.push(y);
        result.stored.push(stored);
    }
    assert!(!result.positions.is_empty());
    result
}
fn close(got: f64, expected: f64, relative: f64) {
    assert!(
        (got - expected).abs() <= relative * expected.abs(),
        "got {got:.17e}, expected {expected:.17e}"
    );
}
fn direct() -> EmgGradientDescent {
    EmgGradientDescent {
        compute_additional_points: false,
        ..Default::default()
    }
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
fn seven_literal_fits_match_source_parameters_counts_and_container_overloads() {
    let fitter = EmgGradientDescent::default();
    for line in include_str!("data/emg_fits.tsv").lines().skip(1) {
        let fields: Vec<_> = line.split('\t').collect();
        let t = trace(fields[0]);
        assert_eq!(t.positions.len(), fields[1].parse::<usize>().unwrap());
        let chromatogram = MSChromatogram {
            peaks: t
                .positions
                .iter()
                .zip(&t.stored)
                .map(|(&x, &y)| ChromatogramPeak::new(x, y))
                .collect(),
            ..Default::default()
        };
        let fit = fitter
            .fit_chromatogram(&chromatogram, None, None)
            .unwrap_or_else(|e| panic!("{}: {e}", fields[0]));
        assert_eq!(
            fit.chromatogram.len(),
            fields[2].parse::<usize>().unwrap(),
            "{}",
            fields[0]
        );
        assert_eq!(
            fit.estimate.training_points,
            fields[7].parse::<usize>().unwrap()
        );
        let p = fit.estimate.parameters;
        for (got, expected) in [p.h, p.mu, p.sigma, p.tau].into_iter().zip(&fields[3..7]) {
            // Source compares parameters after storing them in a float data array.
            close(f64::from(got as f32), expected.parse().unwrap(), 1e-5);
        }
        let spectrum = MSSpectrum {
            peaks: chromatogram
                .peaks
                .iter()
                .map(|p| Peak1D::new(p.rt, p.intensity))
                .collect(),
            ..Default::default()
        };
        let other = fitter.fit_spectrum(&spectrum, None, None).unwrap();
        assert_eq!(other.estimate, fit.estimate);
        assert_eq!(other.spectrum.len(), fit.chromatogram.len());
        for (a, b) in other.spectrum.peaks.iter().zip(&fit.chromatogram.peaks) {
            assert_eq!(a.mz, b.rt);
            assert_eq!(a.intensity, b.intensity);
        }
        if !fields[8].is_empty() {
            // This is the full raw-f64 loss with f32-rounded fitted parameters,
            // not the fit's selected-training-point loss.
            let rounded = EmgParameters {
                h: f64::from(p.h as f32),
                mu: f64::from(p.mu as f32),
                sigma: f64::from(p.sigma as f32),
                tau: f64::from(p.tau as f32),
            };
            let values = direct().apply_parameters(&t.positions, rounded).unwrap();
            let loss = values
                .intensities
                .iter()
                .zip(&t.raw)
                .fold(0.0, |sum, (&y, &raw)| {
                    sum + (y - raw).powf(2.0) / t.raw.len() as f64
                });
            close(loss, fields[8].parse().unwrap(), 1e-5);
        }
    }
}

#[test]
fn source_initial_mean_and_training_count_are_observable_before_updates() {
    let fitter = EmgGradientDescent {
        max_iterations: 1,
        ..direct()
    };
    for (name, mean, count) in [
        ("glutamate", 2.69743333333333, 107),
        ("saturated_min", 2.69516110583333, 77),
        ("saturated_cutoff_sec", 865.1314205, 61),
        ("cutoff_sec", 926.90050115, 12),
    ] {
        let t = trace(name);
        let e = fitter.estimate_parameters(&t.positions, &t.raw).unwrap();
        close(e.parameters.mu, mean, 1e-14);
        assert_eq!(
            e.parameters.h,
            t.raw.iter().copied().fold(f64::NEG_INFINITY, f64::max)
        );
        assert_eq!(e.parameters.sigma, e.parameters.mu * 0.01);
        assert_eq!(e.parameters.tau, e.parameters.sigma * 2.0);
        assert_eq!(e.training_points, count);
        assert_eq!(e.iterations, 1);
        assert_eq!(e.best_iteration, 1);
        assert!(!e.converged);
    }
}

#[test]
fn source_scalar_model_goldens_cover_all_three_branches_and_units() {
    let p = parameters();
    for (x, expected) in [
        (p.mu - 1.0 / 60.0, 1_992_032.657_110_41),
        (p.mu + 1.0 / 60.0, 4_088_964.975_202_13),
        (-3_333_333.0, 0.0),
    ] {
        close(
            direct().apply_parameters(&[x], p).unwrap().intensities[0],
            expected,
            1e-13,
        );
    }
    let p = EmgParameters {
        mu: 860.719,
        sigma: 2.06566,
        tau: 11.3104,
        ..p
    };
    for (x, expected) in [
        (p.mu - 1.0, 1_992_033.065_842_47),
        (p.mu + 1.0, 4_088_968.529_578_75),
        (-200_000_000.0, 0.0),
    ] {
        close(
            direct().apply_parameters(&[x], p).unwrap().intensities[0],
            expected,
            1e-13,
        );
    }
}

#[test]
fn source_fixed_parameter_left_extension_has_exact_count_and_literal_endpoint() {
    let t = trace("saturated_cutoff_min");
    let no_extension = direct()
        .apply_parameters(&t.positions, parameters())
        .unwrap();
    assert_eq!(no_extension.positions, t.positions);
    close(no_extension.intensities[0], 2_144_281.147_222_8, 1e-13);
    let extended = EmgGradientDescent::default()
        .apply_parameters(&t.positions, parameters())
        .unwrap();
    assert_eq!(extended.positions.len(), 71);
    close(extended.positions[0], 14.2717555076923, 1e-14);
    close(extended.intensities[0], 108_845.941_990_663, 1e-12);
    assert_eq!(&extended.positions[5..], &t.positions);
    assert_eq!(&extended.intensities[5..], &no_extension.intensities);
}

#[test]
fn full_call_evaluation_and_generated_point_limits_are_inclusive() {
    let t = trace("saturated_cutoff_min");
    let mut fit = EmgGradientDescent {
        max_points: 71,
        max_evaluations: 71,
        ..Default::default()
    };
    assert_eq!(
        fit.apply_parameters(&t.positions, parameters())
            .unwrap()
            .positions
            .len(),
        71
    );
    fit.max_points = 70;
    assert!(fit.apply_parameters(&t.positions, parameters()).is_err());
    fit.max_points = 71;
    fit.max_evaluations = 70;
    assert!(fit.apply_parameters(&t.positions, parameters()).is_err());

    let t = trace("cutoff_min");
    let input = MSChromatogram {
        peaks: t
            .positions
            .iter()
            .zip(&t.stored)
            .map(|(&x, &y)| ChromatogramPeak::new(x, y))
            .collect(),
        ..Default::default()
    };
    let mut fit = EmgGradientDescent {
        max_iterations: 1,
        compute_additional_points: false,
        max_evaluations: 72,
        ..Default::default()
    };
    assert!(fit.fit_chromatogram(&input, None, None).is_ok());
    fit.max_evaluations = 71;
    assert!(fit.fit_chromatogram(&input, None, None).is_err());
}

#[test]
fn ordered_training_selection_matches_independent_source_translation() {
    let fitter = EmgGradientDescent {
        max_iterations: 1,
        ..direct()
    };
    for row in include_str!("data/emg_training.tsv").lines().skip(1) {
        let fields: Vec<_> = row.split('\t').collect();
        let t = trace(fields[0]);
        let expected_indices: Vec<usize> =
            fields[3].split(',').map(|v| v.parse().unwrap()).collect();
        let fit = fitter.estimate_parameters(&t.positions, &t.raw).unwrap();
        assert_eq!(fit.training_points, expected_indices.len());
        close(fit.parameters.mu, fields[1].parse().unwrap(), 1e-14);
        close(fit.loss, fields[2].parse().unwrap(), 5e-14);
        let predicted = direct()
            .apply_parameters(&t.positions, fit.parameters)
            .unwrap();
        let expected = expected_indices.iter().fold(0.0, |sum, &i| {
            sum + (predicted.intensities[i] - t.raw[i]).powf(2.0) / expected_indices.len() as f64
        });
        // Exactly the same scalar results are accumulated in the source-selected
        // order; sorting that set can alter the last bits of the reported loss.
        assert_eq!(fit.loss.to_bits(), expected.to_bits(), "{}", fields[0]);
    }
    // If the second point already exceeds the 80% threshold, source i==1
    // suppresses derivative calculation. All interior plateau samples are skipped.
    let fit = fitter
        .estimate_parameters(&[100., 101., 102., 103., 104.], &[1., 10., 10., 10., 1.])
        .unwrap();
    assert_eq!(fit.training_points, 2);
}

#[test]
#[allow(clippy::excessive_precision)] // Preserve the independent decimal oracle output.
fn erfc_tails_and_zero_branch_boundary_match_independent_scalar_oracles() {
    // Python stdlib math.erfc applied to the literal source scalar expression;
    // these are independent of the Rust libm implementation and the optimizer.
    let p = EmgParameters {
        h: 1.,
        mu: 0.,
        sigma: 1.,
        tau: 1.,
    };
    for (x, y) in [
        (-26., 5.9773356228680293e-149),
        (-10., 1.7392632041272148e-23),
        (-1., 0.25557335662268754),
        (0., 0.65567954241879833),
        (1., 0.76017345053314034),
        (10., 0.0001876257132043801),
        (100., 1.5374074625819143e-43),
        (700., 4.0747394393902098e-304),
    ] {
        close(
            direct().apply_parameters(&[x], p).unwrap().intensities[0],
            y,
            2e-13,
        );
    }
    let tail = direct().apply_parameters(&[740., 746.], p).unwrap();
    assert!(tail.intensities[0] > 0.0 && tail.intensities[0] < f64::MIN_POSITIVE);
    assert!(
        tail.intensities[0]
            .to_bits()
            .abs_diff(1.7292297604443629e-321_f64.to_bits())
            <= 2
    );
    assert_eq!(tail.intensities[1], 0.0);
    let around = direct()
        .apply_parameters(
            &[
                f64::from_bits(1_f64.to_bits() - 1),
                1.,
                f64::from_bits(1_f64.to_bits() + 1),
            ],
            p,
        )
        .unwrap();
    for y in around.intensities {
        close(y, 0.76017345053314034, 2e-15);
    }
}

#[test]
fn source_exponential_overflow_and_asymptotic_boundary_have_checked_results() {
    let p = EmgParameters {
        h: 1.,
        mu: 0.,
        sigma: 1.,
        tau: 1.,
    };
    assert!(direct().apply_parameters(&[-40.], p).is_err());
    // The source can return a prior finite best after a later gradient failure.
    // The checked native API deliberately returns an error for that trajectory.
    let xs = [100., 101., 102., 10_000.];
    let initial = EmgParameters {
        h: 10.,
        mu: 101.,
        sigma: 1.01,
        tau: 2.02,
    };
    assert!(direct().apply_parameters(&xs, initial).is_ok());
    assert!(
        direct()
            .estimate_parameters(&xs, &[1., 10., 1., 1.])
            .is_err()
    );
    let below_boundary = 1.0 - 2.0_f64.sqrt() * 6.70e7;
    let above_boundary = 1.0 - 2.0_f64.sqrt() * 6.72e7;
    assert!(direct().apply_parameters(&[below_boundary], p).is_err());
    assert_eq!(
        direct()
            .apply_parameters(&[above_boundary], p)
            .unwrap()
            .intensities,
        [0.0]
    );
    for p in [
        EmgParameters { sigma: 0., ..p },
        EmgParameters { tau: 0., ..p },
        EmgParameters {
            h: f64::INFINITY,
            ..p
        },
        EmgParameters { mu: f64::NAN, ..p },
    ] {
        assert!(direct().apply_parameters(&[1.], p).is_err());
    }
}
