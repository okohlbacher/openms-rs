// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::comparison::*;
use openms::format::dta;
use openms::processing::{Normalizer, SpectrumFilter};
use openms::{Error, MSSpectrum, Peak1D, Precursor};

fn spectrum(mzs: &[f64], intensities: &[f32]) -> MSSpectrum {
    MSSpectrum::from_peaks(
        mzs.iter()
            .zip(intensities)
            .map(|(&mz, &i)| Peak1D::new(mz, i))
            .collect(),
    )
}
fn unit(mzs: &[f64]) -> MSSpectrum {
    spectrum(mzs, &vec![1.0; mzs.len()])
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}
fn golden() -> MSSpectrum {
    dta::read(include_bytes!("data/comparison_dfpianger.dta").as_slice()).unwrap()
}
fn config() -> BinConfig {
    BinConfig {
        size: 1.5,
        spread: 2,
        offset: 0.0,
        ..Default::default()
    }
}

#[test]
fn upstream_absolute_and_ppm_alignment_golden() {
    let a = dta::read(include_bytes!("data/comparison_alignment_1.dta").as_slice()).unwrap();
    let b = dta::read(include_bytes!("data/comparison_alignment_2.dta").as_slice()).unwrap();
    let mut aligner = SpectrumAlignment {
        tolerance: Tolerance::Absolute(1.01),
        ..Default::default()
    };
    assert_eq!(
        aligner.align(&a, &b).unwrap(),
        [(0, 0), (1, 1), (3, 3), (4, 5), (6, 6)]
    );
    aligner.tolerance = Tolerance::Ppm(10.0);
    assert_eq!(aligner.align(&a, &b).unwrap(), [(6, 6)]);
    aligner.tolerance = Tolerance::Ppm(10000.0);
    assert_eq!(
        aligner.align(&a, &b).unwrap(),
        [(0, 0), (1, 1), (2, 2), (3, 3), (4, 5), (5, 5), (6, 6)]
    );
    let a = golden();
    assert_eq!(
        SpectrumAlignment::default().align(&a, &a).unwrap().len(),
        a.len()
    );
    let mut b = a.clone();
    b.peaks.truncate(100);
    assert_eq!(
        SpectrumAlignment::default().align(&a, &b).unwrap().len(),
        100
    );
}

#[test]
fn absolute_dp_enforces_one_to_one_and_ppm_is_directed() {
    let a = unit(&[100.0, 100.2]);
    let b = unit(&[100.11]);
    let absolute = SpectrumAlignment {
        tolerance: Tolerance::Absolute(0.15),
        ..Default::default()
    };
    assert_eq!(absolute.align(&a, &b).unwrap(), [(1, 0)]);
    let ppm = SpectrumAlignment {
        tolerance: Tolerance::Ppm(2000.0),
        ..Default::default()
    };
    assert_eq!(ppm.align(&a, &b).unwrap(), [(0, 0), (1, 0)]);
    assert_eq!(ppm.align(&b, &a).unwrap(), [(0, 1)]);
    let tie = SpectrumAlignment {
        tolerance: Tolerance::Ppm(20000.0),
        ..Default::default()
    };
    assert_eq!(
        tie.align(&unit(&[100.0]), &unit(&[99.0, 101.0])).unwrap(),
        [(0, 0)]
    );
}

#[test]
fn source_ppm_float_ties_and_duplicate_stopping_are_preserved() {
    let aligner = SpectrumAlignment {
        tolerance: Tolerance::Ppm(10000.0),
        ..Default::default()
    };
    // MatchedIterator stops on equal adjacent distances, even when a later
    // target peak beyond a duplicate would be closer. Preserve this source quirk.
    assert!(
        aligner
            .align(&unit(&[200.0]), &unit(&[100.0, 100.0, 200.0]))
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        aligner
            .align(&unit(&[100.0, 100.0]), &unit(&[100.0, 100.0]))
            .unwrap(),
        [(0, 0), (1, 0)]
    );
}

#[test]
fn weighted_ppm_rejects_matches_rounded_outside_exact_tolerance() {
    let a = unit(&[1000.0]);
    let b = unit(&[1000.01000000002]);
    let alignment = SpectrumAlignment {
        tolerance: Tolerance::Ppm(10.0),
        ..Default::default()
    };
    assert_eq!(alignment.align(&a, &b).unwrap(), [(0, 0)]);
    let mut score = SpectrumAlignmentScore {
        alignment,
        weighting: DistanceWeighting::None,
    };
    assert!(score.score(&a, &b).unwrap().is_finite());
    for weighting in [DistanceWeighting::Linear, DistanceWeighting::Gaussian] {
        score.weighting = weighting;
        assert!(score.score(&a, &b).is_err());
    }
    let absolute = BinnedSpectrum::new(&unit(&[100.0]), config()).unwrap();
    assert!(absolute.bin_intensity(-1e-100).is_err());
    let ppm = BinnedSpectrum::new(
        &unit(&[100.0]),
        BinConfig {
            unit: BinUnit::Ppm,
            ..config()
        },
    )
    .unwrap();
    assert!(ppm.bin_intensity(0.999999999).is_err());
}

#[test]
fn alignment_empty_exact_boundary_and_validation() {
    let exact = SpectrumAlignment {
        tolerance: Tolerance::Absolute(0.0),
        ..Default::default()
    };
    assert_eq!(
        exact.align(&unit(&[1.0, 2.0]), &unit(&[1.0, 2.1])).unwrap(),
        [(0, 0)]
    );
    assert!(
        exact
            .align(&MSSpectrum::default(), &unit(&[1.0]))
            .unwrap()
            .is_empty()
    );
    assert!(matches!(
        exact.align(&unit(&[2.0, 1.0]), &unit(&[1.0])),
        Err(Error::UnsortedData)
    ));
    for tolerance in [
        Tolerance::Absolute(-1.0),
        Tolerance::Ppm(f64::NAN),
        Tolerance::Ppm(f64::MAX),
    ] {
        assert!(
            SpectrumAlignment {
                tolerance,
                ..Default::default()
            }
            .align(&unit(&[1.0]), &unit(&[1.0]))
            .is_err()
        );
    }
    assert!(exact.align(&unit(&[f64::INFINITY]), &unit(&[1.0])).is_err());
    assert!(
        SpectrumAlignment {
            max_cells: 4,
            ..Default::default()
        }
        .align(&unit(&[1.0, 2.0]), &unit(&[1.0, 2.0]))
        .is_err()
    );
    let two = unit(&[1.0, 2.0]);
    assert!(
        SpectrumAlignment {
            max_cells: 5,
            ..Default::default()
        }
        .align(&two, &two)
        .is_err()
    );
}

#[test]
fn upstream_alignment_score_retains_noncosine_normalization() {
    let mut a = golden();
    Normalizer::default().filter_spectrum(&mut a).unwrap();
    let score = SpectrumAlignmentScore::default();
    // The historical source literal uses absolute tolerance 0.01. Independently
    // summing this fixture's normalized intensities gives the tight expectation.
    close(score.score(&a, &a).unwrap(), 1.48268, 0.01);
    close(score.score(&a, &a).unwrap(), 1.4845010143546342, 1e-14);
    let mut b = a.clone();
    b.peaks.truncate(100);
    close(score.score(&a, &b).unwrap(), 3.82472, 1e-5);
    close(
        score
            .score(&unit(&[1.0]), &spectrum(&[1.0], &[4.0]))
            .unwrap(),
        0.5,
        1e-15,
    );
}

#[test]
fn alignment_weighting_linear_and_gaussian_reference_values() {
    let a = unit(&[10.0]);
    let b = unit(&[10.5]);
    let mut score = SpectrumAlignmentScore {
        alignment: SpectrumAlignment {
            tolerance: Tolerance::Absolute(1.0),
            ..Default::default()
        },
        weighting: DistanceWeighting::Linear,
    };
    close(score.score(&a, &b).unwrap(), 0.5_f64.sqrt(), 1e-15);
    score.weighting = DistanceWeighting::Gaussian;
    // Python math.erfc(0.5/(3*sqrt(2))) = 0.8676323347781927.
    close(
        score.score(&a, &b).unwrap(),
        0.8676323347781927_f64.sqrt(),
        1e-15,
    );
    score.alignment.tolerance = Tolerance::Absolute(0.0);
    close(score.score(&a, &a).unwrap(), 1.0, 1e-15);
    assert_eq!(score.score(&a, &MSSpectrum::default()).unwrap(), 0.0);
    assert_eq!(score.score(&spectrum(&[10.0], &[0.0]), &a).unwrap(), 0.0);
    assert!(score.score(&spectrum(&[10.0], &[-1.0]), &a).is_err());
}

#[test]
fn upstream_bin_indices_and_lower_boundaries() {
    let cfg = BinConfig {
        size: 10.0,
        unit: BinUnit::Ppm,
        spread: 0,
        offset: 0.0,
        ..Default::default()
    };
    for (mz, index) in [(1.0, 0), (10.0, 230259), (100.0, 460519), (1000.0, 690778)] {
        assert_eq!(cfg.bin_index(mz).unwrap(), index);
    }
    close(
        f64::from(cfg.bin_lower_mz(1000).unwrap()),
        (1.0_f64 + 10e-6).powi(1000),
        1e-7,
    );
    let cfg = BinConfig {
        size: 1.0,
        spread: 0,
        offset: 0.5,
        ..Default::default()
    };
    for mz in [1.0, 10.0, 100.0, 1000.0] {
        assert_eq!(
            cfg.bin_index(mz - 0.01).unwrap(),
            cfg.bin_index(mz + 0.01).unwrap()
        );
    }
    assert_eq!(cfg.bin_lower_mz(0).unwrap(), -0.5);
    assert_eq!(cfg.bin_lower_mz(1000).unwrap(), 999.5);
}

#[test]
fn upstream_binned_storage_and_comparison_golden() {
    let a = golden();
    let mut b = a.clone();
    b.peaks.pop();
    let ba = BinnedSpectrum::new(&a, config()).unwrap();
    let bb = BinnedSpectrum::new(&b, config()).unwrap();
    assert_eq!(ba.bins().len(), 347);
    assert_eq!(ba.bins()[&658], 501645.0);
    assert_eq!(ba.precursors(), a.precursors);
    close(binned_shared_peak_count(&ba, &bb).unwrap(), 0.997118, 1e-6);
    close(
        binned_sum_agreeing_intensities(&ba, &bb).unwrap(),
        0.99707,
        1e-5,
    );
    let cfg = BinConfig {
        offset: 0.4,
        ..config()
    };
    let ba = BinnedSpectrum::new(&a, cfg).unwrap();
    let bb = BinnedSpectrum::new(&b, cfg).unwrap();
    close(binned_cosine(&ba, &bb).unwrap(), 0.999985, 1e-5);
    close(binned_cosine(&ba, &ba).unwrap(), 1.0, 1e-15);
    close(
        binned_sum_agreeing_intensities(&ba, &ba).unwrap(),
        1.0,
        1e-15,
    );
}

#[test]
fn bin_spread_boundary_overlap_zeros_and_readonly_lookup() {
    let cfg = BinConfig {
        size: 1.0,
        spread: 1,
        offset: 0.0,
        ..Default::default()
    };
    let bins = BinnedSpectrum::new(&spectrum(&[0.0, 1.0], &[2.0, 3.0]), cfg).unwrap();
    assert_eq!(
        bins.bins()
            .iter()
            .map(|(&i, &v)| (i, v))
            .collect::<Vec<_>>(),
        [(0, 5.0), (1, 5.0), (2, 3.0)]
    );
    assert_eq!(bins.bin_intensity(99.0).unwrap(), 0.0);
    assert_eq!(bins.bins().len(), 3);
    let zero = BinnedSpectrum::new(&spectrum(&[1.0], &[0.0]), cfg).unwrap();
    assert_eq!(binned_cosine(&zero, &zero).unwrap(), 0.0);
    assert_eq!(binned_shared_peak_count(&zero, &zero).unwrap(), 1.0);
    assert_eq!(binned_sum_agreeing_intensities(&zero, &zero).unwrap(), 0.0);
    let empty = BinnedSpectrum::new(&MSSpectrum::default(), cfg).unwrap();
    assert_eq!(binned_shared_peak_count(&empty, &empty).unwrap(), 0.0);
}

#[test]
fn binning_rejects_bad_config_ranges_overflow_and_limits() {
    for cfg in [
        BinConfig {
            size: 0.0,
            ..config()
        },
        BinConfig {
            size: f32::NAN,
            ..config()
        },
        BinConfig {
            offset: f32::INFINITY,
            ..config()
        },
        BinConfig {
            max_bins: 0,
            ..config()
        },
    ] {
        assert!(BinnedSpectrum::new(&unit(&[1.0]), cfg).is_err());
    }
    let ppm = BinConfig {
        unit: BinUnit::Ppm,
        ..config()
    };
    assert!(BinnedSpectrum::new(&unit(&[0.99]), ppm).is_err());
    assert!(config().bin_index(f64::MAX).is_err());
    assert!(config().bin_index(-1e-100).is_err());
    assert!(
        BinConfig {
            unit: BinUnit::Ppm,
            ..config()
        }
        .bin_index(0.999999999)
        .is_err()
    );
    assert!(BinnedSpectrum::new(&unit(&[2.0, 1.0]), config()).is_err());
    assert!(
        BinnedSpectrum::new(
            &unit(&[1.0]),
            BinConfig {
                max_bins: 1,
                ..config()
            }
        )
        .is_err()
    );
    assert!(
        BinnedSpectrum::new(
            &unit(&[1.0]),
            BinConfig {
                spread: u32::MAX,
                ..config()
            }
        )
        .is_err()
    );
    assert!(BinnedSpectrum::new(&spectrum(&[1.0, 1.0], &[f32::MAX, f32::MAX]), config()).is_err());
    assert!(ppm.bin_lower_mz(usize::MAX).is_err());
}

#[test]
fn binned_compatibility_signed_cosine_and_symmetric_scores() {
    let cfg = BinConfig {
        size: 1.0,
        spread: 0,
        offset: 0.0,
        ..Default::default()
    };
    let a = BinnedSpectrum::new(&spectrum(&[1.0, 2.0], &[1.0, 2.0]), cfg).unwrap();
    let neg = BinnedSpectrum::new(&spectrum(&[1.0, 2.0], &[-1.0, -2.0]), cfg).unwrap();
    close(binned_cosine(&a, &neg).unwrap(), -1.0, 1e-15);
    // binned_sum_agreeing_intensities is now the BinnedSumAgreeingIntensities
    // implementation itself rather than a second one, and the source accepts
    // negative bins: it truncates every coefficient of
    // (a + b)/2 - |a - b| that falls below zero. Here that is every bin, and the
    // mean total intensity (3 + -3)/2 is zero, so the score is the guarded 0.0
    // instead of the error the separate f64 implementation used to return.
    assert_eq!(binned_sum_agreeing_intensities(&a, &neg).unwrap(), 0.0);
    let incompatible =
        BinnedSpectrum::new(&unit(&[1.0]), BinConfig { offset: 0.5, ..cfg }).unwrap();
    assert!(binned_cosine(&a, &incompatible).is_err());
    assert!(binned_shared_peak_count(&a, &incompatible).is_err());
    assert!(binned_sum_agreeing_intensities(&a, &incompatible).is_err());
    let spread = BinnedSpectrum::new(&unit(&[1.0]), BinConfig { spread: 1, ..cfg }).unwrap();
    assert!(a.is_compatible(&spread));
    close(
        binned_cosine(&a, &spread).unwrap(),
        binned_cosine(&spread, &a).unwrap(),
        1e-15,
    );
    close(
        binned_sum_agreeing_intensities(&a, &spread).unwrap(),
        binned_sum_agreeing_intensities(&spread, &a).unwrap(),
        1e-15,
    );
}

#[test]
fn upstream_precursor_comparison_and_missing_convention() {
    let a = dta::read(include_bytes!("data/Transformers_tests.dta").as_slice()).unwrap();
    let b = dta::read(include_bytes!("data/comparison_transformers_2.dta").as_slice()).unwrap();
    let compare = SpectrumPrecursorComparator::default();
    close(compare.score(&a, &b).unwrap(), 1.7685, 1e-10);
    assert_eq!(compare.score(&a, &a).unwrap(), 2.0);
    assert_eq!(
        compare
            .score(&MSSpectrum::default(), &MSSpectrum::default())
            .unwrap(),
        2.0
    );
    let missing = MSSpectrum::default();
    let distant = MSSpectrum {
        precursors: vec![Precursor::new(5.0, 1)],
        ..Default::default()
    };
    assert_eq!(compare.score(&missing, &distant).unwrap(), 0.0);
}

#[test]
fn upstream_zhang_golden_and_stein_scott_formula() {
    let mut a = golden();
    Normalizer::default().filter_spectrum(&mut a).unwrap();
    let zhang = ZhangSimilarityScore::default();
    close(zhang.score(&a, &a).unwrap(), 1.82682, 1e-5);
    let mut b = a.clone();
    b.peaks.truncate(100);
    close(zhang.score(&a, &b).unwrap(), 0.328749, 1e-6);
    let a = spectrum(
        &[500.0, 600.0, 700.0, 800.0, 900.0],
        &[500.0, 600.0, 700.0, 800.0, 900.0],
    );
    let expected = 1.0 - (0.2 / 10000.0) * 3500.0_f64.powi(2) / 2550000.0;
    close(
        SteinScottImproveScore::default().score(&a, &a).unwrap(),
        expected,
        1e-15,
    );
}

#[test]
fn many_to_many_strict_vs_inclusive_boundaries_and_limits() {
    let a = unit(&[1.0]);
    let b = unit(&[2.0]);
    let z = ZhangSimilarityScore {
        tolerance: 1.0,
        ..Default::default()
    };
    assert_eq!(z.score(&a, &b).unwrap(), 0.0);
    let s = SteinScottImproveScore {
        tolerance: 0.5,
        threshold: 0.0,
        ..Default::default()
    };
    close(s.score(&a, &b).unwrap(), 0.99995, 1e-15);
    let repeated = unit(&[1.0, 1.0, 1.0]);
    assert!(
        ZhangSimilarityScore {
            max_pairs: 2,
            ..Default::default()
        }
        .score(&repeated, &repeated)
        .is_err()
    );
    assert!(
        SteinScottImproveScore {
            max_pairs: 2,
            ..Default::default()
        }
        .score(&repeated, &repeated)
        .is_err()
    );
    assert_eq!(z.score(&MSSpectrum::default(), &a).unwrap(), 0.0);
    assert!(z.score(&spectrum(&[1.0], &[-1.0]), &a).is_err());
}

#[test]
fn gaussian_scale_is_instance_specific_instead_of_upstream_static_cache() {
    let a = unit(&[1.0]);
    let b = unit(&[1.1]);
    let narrow = ZhangSimilarityScore {
        tolerance: 0.2,
        weighting: DistanceWeighting::Gaussian,
        ..Default::default()
    };
    let broad = ZhangSimilarityScore {
        tolerance: 2.0,
        ..narrow
    };
    let n = narrow.score(&a, &b).unwrap();
    let wide = broad.score(&a, &b).unwrap();
    assert!(n < wide && wide <= 1.0);
    assert_eq!(n, narrow.score(&a, &b).unwrap());
}

// Independent map-shaped transcription of upstream SpectrumAlignment.h. The
// production implementation uses compact contiguous rows instead of these maps.
fn source_dp(a: &MSSpectrum, b: &MSSpectrum, tolerance: f64) -> Vec<(usize, usize)> {
    use std::collections::BTreeMap;
    let mut scores = BTreeMap::new();
    let mut trace = BTreeMap::new();
    for i in 0..=a.len() {
        scores.insert((i, 0), i as f64 * tolerance);
    }
    for j in 0..=b.len() {
        scores.insert((0, j), j as f64 * tolerance);
    }
    let mut left = 1;
    let mut last = (0, 0);
    for i in 1..=a.len() {
        let mut j = left;
        while j <= b.len() {
            let p = a.peaks[i - 1].mz;
            let q = b.peaks[j - 1].mz;
            let d = (p - q).abs();
            let stop = q > p && d > tolerance && i < a.len() && j < b.len() && a.peaks[i].mz < q;
            if p > q && d > tolerance && j > left + 1 {
                left += 1;
            }
            let get = |r, c| *scores.get(&(r, c)).unwrap_or(&((r + c) as f64 * tolerance));
            let diag = d + get(i - 1, j - 1);
            let horizontal = tolerance + get(i, j - 1);
            let vertical = tolerance + get(i - 1, j);
            let (cost, previous) = if diag <= horizontal && diag <= vertical && d <= tolerance {
                last = (i, j);
                (diag, (i - 1, j - 1))
            } else if horizontal <= vertical {
                (horizontal, (i, j - 1))
            } else {
                (vertical, (i - 1, j))
            };
            scores.insert((i, j), cost);
            trace.insert((i, j), previous);
            if stop {
                break;
            }
            j += 1;
        }
    }
    let mut out = Vec::new();
    let (mut i, mut j) = last;
    while i > 0 && j > 0 {
        let previous = trace.get(&(i, j)).copied().unwrap_or((0, 0));
        if previous == (i - 1, j - 1) {
            out.push((i - 1, j - 1));
        }
        (i, j) = previous;
    }
    out.reverse();
    out
}

#[test]
fn compact_banded_dp_matches_source_map_oracle_across_boundaries() {
    let mut seed = 77_u64;
    for case in 0..200 {
        let mut next = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((seed >> 32) % 31) as f64 * 0.25
        };
        let mut a: Vec<_> = (0..case % 8).map(|_| next()).collect();
        let mut b: Vec<_> = (0..(case / 3) % 9).map(|_| next()).collect();
        a.sort_by(f64::total_cmp);
        b.sort_by(f64::total_cmp);
        let a = unit(&a);
        let b = unit(&b);
        for tolerance in [0.0, 0.25, 0.75, 2.0, 10.0] {
            let config = SpectrumAlignment {
                tolerance: Tolerance::Absolute(tolerance),
                ..Default::default()
            };
            let found = config.align(&a, &b).unwrap();
            assert_eq!(
                found,
                source_dp(&a, &b, tolerance),
                "case {case} tol {tolerance}"
            );
            assert!(found.windows(2).all(|p| p[0].0 < p[1].0 && p[0].1 < p[1].1));
            assert!(
                found
                    .iter()
                    .all(|&(i, j)| (a.peaks[i].mz - b.peaks[j].mz).abs() <= tolerance)
            );
        }
    }
}

#[test]
fn pair_resource_bound_includes_rejected_strict_boundary_candidates() {
    let a = unit(&[1.0, 1.0, 1.0]);
    let comparison = ZhangSimilarityScore {
        tolerance: 0.0,
        max_pairs: 2,
        ..Default::default()
    };
    assert!(comparison.score(&a, &a).is_err());
}

#[test]
fn binned_equality_compares_data_and_layout_not_resource_budgets() {
    let a = BinnedSpectrum::new(&unit(&[1.0]), config()).unwrap();
    let b = BinnedSpectrum::new(
        &unit(&[1.0]),
        BinConfig {
            max_bins: 100,
            max_updates: 100,
            ..config()
        },
    )
    .unwrap();
    assert_eq!(a, b);
    let shifted = BinnedSpectrum::new(
        &MSSpectrum::default(),
        BinConfig {
            offset: 0.5,
            ..config()
        },
    )
    .unwrap();
    let unshifted = BinnedSpectrum::new(&MSSpectrum::default(), config()).unwrap();
    assert_ne!(shifted, unshifted);
}
