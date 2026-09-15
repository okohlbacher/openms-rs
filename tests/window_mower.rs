// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::format::dta;
use openms::kernel::DataArray;
use openms::processing::SpectrumFilter;
use openms::processing::window_mower::{WindowMower, WindowMowerMethod};
use openms::{ChromatogramPeak, MSChromatogram, MSExperiment, MSSpectrum, Peak1D};
use std::io::Cursor;

fn spectrum(points: &[(f64, f32)]) -> MSSpectrum {
    let mut input =
        MSSpectrum::from_peaks(points.iter().map(|&(x, y)| Peak1D::new(x, y)).collect());
    input.name = "original spectrum".into();
    input.metadata.insert("note".into(), "preserve this".into());
    input.float_data_arrays.push(DataArray::new(
        "float marker",
        (0..points.len()).map(|i| i as f32 + 0.5).collect(),
    ));
    input
        .float_data_arrays
        .push(DataArray::new("empty placeholder", vec![]));
    input.integer_data_arrays.push(DataArray::new(
        "original index",
        (0..points.len() as i32).collect(),
    ));
    input.string_data_arrays.push(DataArray::new(
        "string marker",
        (0..points.len()).map(|i| format!("p{i}")).collect(),
    ));
    input
}
fn options(method: WindowMowerMethod, window_size: f64, peak_count: usize) -> WindowMower {
    WindowMower {
        method,
        window_size,
        peak_count,
        ..Default::default()
    }
}
fn selected(input: &MSSpectrum, expected: &[usize]) -> MSSpectrum {
    let mut result = input.clone();
    result.select(expected).unwrap();
    result
}

#[test]
fn pinned_transformers_fixture_has_source_56_sliding_and_30_jumping_peaks() {
    // Literal WindowMower_test.cpp assertions against its unchanged DTA fixture.
    let input = dta::read(Cursor::new(include_bytes!("data/Transformers_tests.dta"))).unwrap();
    assert_eq!(input.len(), 121);
    let default = WindowMower::default();
    assert_eq!(
        (default.window_size, default.peak_count, default.method),
        (50.0, 2, WindowMowerMethod::Sliding)
    );
    assert_eq!(default.filtered_spectrum(&input).unwrap().len(), 56);
    let mut jumping = input.clone();
    default.filter_jumping(&mut jumping).unwrap();
    assert_eq!(jumping.len(), 30);
    let mut sliding = input.clone();
    options(WindowMowerMethod::Jumping, 50.0, 2)
        .filter_sliding(&mut sliding)
        .unwrap();
    assert_eq!(sliding.len(), 56);
    assert_eq!(
        "slide".parse::<WindowMowerMethod>().unwrap(),
        WindowMowerMethod::Sliding
    );
    assert_eq!(
        "jump".parse::<WindowMowerMethod>().unwrap(),
        WindowMowerMethod::Jumping
    );
    assert!("moving".parse::<WindowMowerMethod>().is_err());
}

#[test]
fn source_triangle_retains_exact_peak_and_annotation_indices() {
    // Independent source class-test triangle; the comment claiming one final
    // retained peak for width10 contradicts round(0.9*2)=2 and the asserted20.
    let points: Vec<_> = (0..100)
        .map(|i| {
            (
                i as f64,
                if i < 50 {
                    (i as f64 + 0.1) as f32
                } else {
                    ((100 - i) as f64 + 0.2) as f32
                },
            )
        })
        .collect();
    let input = spectrum(&points);
    for method in [WindowMowerMethod::Sliding, WindowMowerMethod::Jumping] {
        assert_eq!(
            options(method, 50.0, 2).filtered_spectrum(&input).unwrap(),
            selected(&input, &[48, 49, 50, 51])
        );
    }
    let expected = [
        8, 9, 18, 19, 28, 29, 38, 39, 48, 49, 50, 51, 60, 61, 70, 71, 80, 81, 90, 91,
    ];
    assert_eq!(
        options(WindowMowerMethod::Jumping, 10.0, 2)
            .filtered_spectrum(&input)
            .unwrap(),
        selected(&input, &expected)
    );
}

#[test]
fn sliding_stops_at_first_end_reaching_window_and_preserves_original_order() {
    let input = spectrum(&[(2.0, 2.0), (0.0, 100.0), (1.0, 1.0)]);
    let config = options(WindowMowerMethod::Sliding, 50.0, 1);
    assert_eq!(config.retained_indices(&input).unwrap(), vec![1]);
    let input = spectrum(&[(2.0, 10.0), (0.0, 9.0), (1.0, 1.0)]);
    assert_eq!(
        options(WindowMowerMethod::Sliding, 2.0, 1)
            .retained_indices(&input)
            .unwrap(),
        vec![0, 1]
    );
    assert_eq!(
        options(WindowMowerMethod::Jumping, 2.0, 1)
            .retained_indices(&input)
            .unwrap(),
        vec![1]
    );
}

#[test]
fn strict_windows_reset_after_gaps_and_scale_the_last_quota() {
    let input = spectrum(&[(0.0, 1.0), (9.0, 10.0), (10.0, 100.0), (100.0, 1000.0)]);
    assert_eq!(
        options(WindowMowerMethod::Jumping, 10.0, 1)
            .retained_indices(&input)
            .unwrap(),
        vec![1, 2]
    );
    let input = spectrum(&[(100.0, 1.0)]);
    assert!(
        options(WindowMowerMethod::Jumping, 10.0, 2)
            .retained_indices(&input)
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        options(WindowMowerMethod::Sliding, 10.0, 2)
            .retained_indices(&input)
            .unwrap(),
        vec![0]
    );
    for (span, expected) in [
        (f64::from_bits(2.5_f64.to_bits() - 1), 0),
        (2.5, 1),
        (f64::from_bits(2.5_f64.to_bits() + 1), 1),
    ] {
        let input = spectrum(&[(0.0, 1.0), (span, 2.0)]);
        assert_eq!(
            options(WindowMowerMethod::Jumping, 10.0, 2)
                .retained_indices(&input)
                .unwrap()
                .len(),
            expected
        );
    }
}

#[test]
fn duplicate_membership_differs_between_methods_and_treats_signed_zeros_as_equal() {
    let input = spectrum(&[(10.0, 1.0), (0.0, 100.0), (-0.0, 100.0), (0.0, 1.0)]);
    assert_eq!(
        options(WindowMowerMethod::Sliding, 50.0, 1)
            .retained_indices(&input)
            .unwrap(),
        vec![1, 2, 3]
    );
    assert_eq!(
        options(WindowMowerMethod::Jumping, 10.0, 1)
            .retained_indices(&input)
            .unwrap(),
        vec![1, 2]
    );
    let zeros = spectrum(&[(0.0, -0.0), (0.0, 0.0), (10.0, 1.0)]);
    assert_eq!(
        options(WindowMowerMethod::Jumping, 10.0, 1)
            .retained_indices(&zeros)
            .unwrap(),
        vec![0, 1]
    );
    let before = input.clone();
    for method in [WindowMowerMethod::Sliding, WindowMowerMethod::Jumping] {
        let result = options(method, 10.0, 1).filtered_spectrum(&input).unwrap();
        result.validate().unwrap();
        assert_eq!(input, before);
    }
}

#[test]
fn stable_ties_and_signed_intensities_are_deterministic() {
    let input = spectrum(&[(2.0, -1.0), (0.0, -1.0), (1.0, -2.0), (10.0, -3.0)]);
    assert_eq!(
        options(WindowMowerMethod::Sliding, 50.0, 1)
            .retained_indices(&input)
            .unwrap(),
        vec![1]
    );
    assert_eq!(
        options(WindowMowerMethod::Jumping, 10.0, 1)
            .retained_indices(&input)
            .unwrap(),
        vec![1]
    );
    let all = spectrum(&[(2.0, 0.0), (0.0, 0.0), (1.0, 0.0)]);
    assert_eq!(
        options(WindowMowerMethod::Sliding, 50.0, 2)
            .retained_indices(&all)
            .unwrap(),
        vec![1, 2]
    );
}

// Deliberately simple source interpreter: clone every window, fully stable-sort
// by intensity, and use linear equality membership. Tests keep these cases tiny.
fn reference(
    input: &MSSpectrum,
    method: WindowMowerMethod,
    width: f64,
    quota: usize,
) -> Vec<usize> {
    let mut order: Vec<_> = (0..input.len()).collect();
    order.sort_by(|&a, &b| input.peaks[a].mz.partial_cmp(&input.peaks[b].mz).unwrap());
    let mut chosen = Vec::new();
    let mut start = 0;
    while start < order.len() {
        let mut end = start;
        while end < order.len() && input.peaks[order[end]].mz - input.peaks[order[start]].mz < width
        {
            end += 1;
        }
        let count = if method == WindowMowerMethod::Jumping && end == order.len() {
            ((input.peaks[order[end - 1]].mz - input.peaks[order[start]].mz) / width * quota as f64)
                .round() as usize
        } else {
            quota
        };
        let mut current = order[start..end].to_vec();
        current.sort_by(|&a, &b| {
            input.peaks[b]
                .intensity
                .partial_cmp(&input.peaks[a].intensity)
                .unwrap()
        });
        chosen.extend(current.into_iter().take(count));
        if end == order.len() {
            break;
        }
        start = if method == WindowMowerMethod::Sliding {
            start + 1
        } else {
            end
        };
    }
    let output_order: Vec<_> = if method == WindowMowerMethod::Sliding {
        (0..input.len()).collect()
    } else {
        order
    };
    output_order
        .into_iter()
        .filter(|&i| {
            chosen.iter().any(|&j| {
                if method == WindowMowerMethod::Sliding {
                    input.peaks[i].mz == input.peaks[j].mz
                } else {
                    input.peaks[i] == input.peaks[j]
                }
            })
        })
        .collect()
}

#[test]
fn optimized_membership_matches_independent_small_source_interpreter() {
    for seed in 0..15 {
        let points: Vec<_> = (0..19)
            .map(|i| {
                (
                    ((i * 7 + seed * 3) % 13) as f64 - 3.0,
                    ((i * 11 + seed) % 7) as f32 - 3.0,
                )
            })
            .collect();
        let input = spectrum(&points);
        for method in [WindowMowerMethod::Sliding, WindowMowerMethod::Jumping] {
            for width in [0.5, 2.0, 5.5, 50.0] {
                for quota in [0, 1, 3, usize::MAX] {
                    assert_eq!(
                        options(method, width, quota)
                            .retained_indices(&input)
                            .unwrap(),
                        reference(&input, method, width, quota),
                        "seed{seed}, {method:?}, width{width}, quota{quota}"
                    );
                }
            }
        }
    }
}

#[test]
fn empty_zero_quota_and_huge_quota_have_defined_behavior() {
    for method in [WindowMowerMethod::Sliding, WindowMowerMethod::Jumping] {
        assert!(
            options(method, 10.0, 2)
                .filtered_spectrum(&MSSpectrum::default())
                .unwrap()
                .is_empty()
        );
        let input = spectrum(&[(0.0, 1.0), (1.0, 2.0)]);
        assert_eq!(
            options(method, 10.0, 0).filtered_spectrum(&input).unwrap(),
            selected(&input, &[])
        );
        assert_eq!(
            options(method, 10.0, usize::MAX)
                .retained_indices(&input)
                .unwrap(),
            vec![0, 1]
        );
    }
}

#[test]
fn invalid_data_parameters_and_work_fail_without_mutating_spectrum() {
    let input = spectrum(&[(0.0, 1.0), (1.0, 2.0), (2.0, 3.0)]);
    for config in [
        WindowMower {
            window_size: 0.0,
            ..Default::default()
        },
        WindowMower {
            window_size: -1.0,
            ..Default::default()
        },
        WindowMower {
            window_size: f64::NAN,
            ..Default::default()
        },
        WindowMower {
            window_size: f64::INFINITY,
            ..Default::default()
        },
        WindowMower {
            max_points: 0,
            ..Default::default()
        },
        WindowMower {
            max_points: 2,
            ..Default::default()
        },
        // The work ceiling is `max_work + work_per_point * points`, so a work
        // refusal needs both terms small: the last case leaves the rate
        // positive and still lands under what these three peaks cost.
        WindowMower {
            max_work: 0,
            work_per_point: 0,
            ..Default::default()
        },
        WindowMower {
            max_work: 1,
            work_per_point: 0,
            ..Default::default()
        },
        WindowMower {
            max_work: 20,
            work_per_point: 0,
            ..Default::default()
        },
        WindowMower {
            max_work: 0,
            work_per_point: 6,
            ..Default::default()
        },
    ] {
        let mut value = input.clone();
        assert!(config.filter_spectrum(&mut value).is_err());
        assert_eq!(value, input);
    }
    let mut bad = input.clone();
    bad.integer_data_arrays[0].data.pop();
    let before = bad.clone();
    assert!(WindowMower::default().filter_spectrum(&mut bad).is_err());
    assert_eq!(bad, before);
    let mut span = spectrum(&[(-f64::MAX, 1.0), (f64::MAX, 2.0)]);
    let before = span.clone();
    assert!(WindowMower::default().filter_spectrum(&mut span).is_err());
    assert_eq!(span, before);
    let nonfinite = spectrum(&[(0.0, f32::INFINITY)]);
    assert!(WindowMower::default().retained_indices(&nonfinite).is_err());
    let nonfinite = spectrum(&[(f64::NAN, 1.0)]);
    assert!(WindowMower::default().retained_indices(&nonfinite).is_err());
}

#[test]
fn experiment_limits_are_per_record_and_errors_do_not_partially_commit() {
    let valid = spectrum(&[(2.0, 2.0), (0.0, 100.0), (1.0, 1.0)]);
    let mut experiment = MSExperiment {
        spectra: vec![valid.clone(), valid.clone()],
        chromatograms: vec![MSChromatogram::from_peaks(vec![ChromatogramPeak::new(
            1.0, 5.0,
        )])],
        ..Default::default()
    };
    experiment
        .settings
        .metadata
        .insert("experiment".into(), "keep".into());
    let before = experiment.clone();
    // Both ceilings are per record, so a configuration that accepts one
    // spectrum accepts an experiment of any number of copies of it. A run-wide
    // ledger shrank as the run grew and refused real data: the benchmark's
    // 1.2 GB Velos run holds 88,434,492 peaks in 43,745 spectra, 88 times the
    // default million, while its largest spectrum holds 16,766.
    let tight = WindowMower {
        max_work: 60,
        work_per_point: 0,
        max_points: 3,
        ..Default::default()
    };
    assert!(tight.filtered_spectrum(&valid).is_ok());
    let mut many = before.clone();
    many.spectra = vec![valid.clone(); 64];
    tight.filter_experiment(&mut many).unwrap();
    let expected = tight.filtered_spectrum(&valid).unwrap();
    assert_eq!(many.spectra, vec![expected; 64]);
    // A record that exceeds either ceiling on its own is still refused, and the
    // experiment is left untouched.
    for config in [
        WindowMower {
            max_points: 2,
            ..Default::default()
        },
        WindowMower {
            max_work: 20,
            work_per_point: 0,
            ..Default::default()
        },
        WindowMower {
            max_points: 0,
            ..Default::default()
        },
    ] {
        assert!(config.filter_experiment(&mut experiment).is_err());
        assert_eq!(experiment, before);
    }
    let mut invalid = before.clone();
    invalid.spectra[1].string_data_arrays[0].data.pop();
    let saved = invalid.clone();
    assert!(
        WindowMower::default()
            .filter_experiment(&mut invalid)
            .is_err()
    );
    assert_eq!(invalid, saved);
    let expected = WindowMower::default().filtered_spectrum(&valid).unwrap();
    WindowMower::default()
        .filter_experiment(&mut experiment)
        .unwrap();
    assert_eq!(experiment.spectra, vec![expected.clone(), expected]);
    assert_eq!(experiment.chromatograms, before.chromatograms);
    assert_eq!(experiment.settings.metadata, before.settings.metadata);
}

/// The work a spectrum can cost is bounded by its own point count, and the
/// defaults admit the largest spectrum of a real run.
///
/// Sliding cost is `Θ(n · w)` for mean window occupancy `w`, and `w` is set by
/// the data, not by `n`: peaks spread over about twice the window width put
/// `n/2` points in each of `n/2` windows, which is quadratic. The ceiling is
/// `max_work + work_per_point * n`, so such a spectrum is refused after work
/// linear in its own size instead of running to completion, while a spectrum of
/// the same point count at a realistic occupancy is accepted. The accepted case
/// under the defaults is the largest spectrum of the benchmark's 1.2 GB LTQ
/// Orbitrap Velos run: 16,766 peaks over about 1,800 Th.
#[test]
fn work_is_bounded_by_the_point_count_and_the_defaults_admit_a_real_spectrum() {
    let over = |count: u32, span: f64| -> Vec<(f64, f32)> {
        (0..count)
            .map(|i| {
                (
                    200.0 + f64::from(i) * span / f64::from(count),
                    i as f32 + 1.0,
                )
            })
            .collect()
    };
    assert!(
        WindowMower::default()
            .retained_indices(&spectrum(&over(16_766, 1_800.0)))
            .is_ok()
    );
    // A record over the per-spectrum point ceiling is refused before any work.
    let huge = MSSpectrum::from_peaks(
        (0..1_000_001)
            .map(|i| Peak1D::new(200.0 + f64::from(i) * 1e-3, 1.0))
            .collect(),
    );
    assert!(WindowMower::default().retained_indices(&huge).is_err());

    // Same point count, same ceiling: occupancy alone decides. At 1,024 units
    // per point a 2,000-peak spectrum may spend 2,048,000; spread over 1,800 Th
    // it needs about a quarter of that, packed into 100 Th about twice it.
    let ceiling = WindowMower {
        max_work: 0,
        work_per_point: 1_024,
        ..Default::default()
    };
    assert!(
        ceiling
            .retained_indices(&spectrum(&over(2_000, 1_800.0)))
            .is_ok()
    );
    assert!(
        ceiling
            .retained_indices(&spectrum(&over(2_000, 100.0)))
            .is_err()
    );
}
