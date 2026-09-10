// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::analysis::psm_scoring::*;
use openms::chemistry::{AASequence, TheoreticalSpectrumGenerator};
use openms::comparison::Tolerance;
use openms::kernel::DataArray;
use openms::{MSSpectrum, Peak1D};

fn theo(sequence: &str, max_charge: u8) -> MSSpectrum {
    TheoreticalSpectrumGenerator {
        add_metainfo: true,
        ..Default::default()
    }
    .generate(&AASequence::parse(sequence).unwrap(), 1, max_charge, None)
    .unwrap()
}
fn spectrum(peaks: &[(f64, f32)], names: &[&str]) -> MSSpectrum {
    MSSpectrum {
        peaks: peaks.iter().map(|&(mz, i)| Peak1D::new(mz, i)).collect(),
        string_data_arrays: if names.is_empty() {
            vec![]
        } else {
            vec![DataArray::new(
                "IonNames",
                names.iter().map(|s| (*s).into()).collect(),
            )]
        },
        ..Default::default()
    }
}
fn near(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t, "{a} != {b}");
}

#[test]
fn source_peptide_psm_goldens_and_empty_results() {
    let mut theoretical = theo("PEPTIDE", 1);
    assert_eq!(theoretical.len(), 11);
    for tolerance in [Tolerance::Absolute(0.1), Tolerance::Ppm(10.0)] {
        let hyper = HyperScore { tolerance };
        let morpheus = MorpheusScore { tolerance };
        near(
            hyper.compute(&theoretical, &theoretical).unwrap(),
            13.8516496,
            1e-7,
        );
        assert_eq!(
            morpheus.compute(&theoretical, &theoretical).unwrap().score,
            12.0
        );
        assert_eq!(
            hyper.compute(&MSSpectrum::default(), &theoretical).unwrap(),
            0.0
        );
        assert_eq!(
            morpheus
                .compute(&MSSpectrum::default(), &theoretical)
                .unwrap(),
            MorpheusResult::default()
        );
    }
    theoretical = theo("PEPTIDE", 3);
    for tolerance in [Tolerance::Absolute(0.1), Tolerance::Ppm(10.0)] {
        near(
            HyperScore { tolerance }
                .compute(&theoretical, &theoretical)
                .unwrap(),
            67.8210771,
            1e-7,
        );
        assert_eq!(
            MorpheusScore { tolerance }
                .compute(&theoretical, &theoretical)
                .unwrap()
                .score,
            34.0
        );
    }
    let unmatched = theo("YYYYYY", 3);
    assert_eq!(
        HyperScore {
            tolerance: Tolerance::Absolute(1e-5)
        }
        .compute(&unmatched, &theoretical)
        .unwrap(),
        0.0
    );
    let no_match = MorpheusScore {
        tolerance: Tolerance::Absolute(1e-5),
    }
    .compute(&theo("PEPTIDE", 1), &theo("EDITPEP", 1))
    .unwrap();
    assert_eq!(
        (no_match.matches, no_match.score, no_match.mean_error_da),
        (0, 0.0, 1e10)
    );
}

#[test]
fn source_squared_mass_example_distinguishes_ppm_and_da_matching() {
    let mut theoretical = theo("PEPTIDE", 3);
    let mut experimental = theoretical.clone();
    for (t, e) in theoretical.peaks.iter_mut().zip(&mut experimental.peaks) {
        let mz = t.mz.powi(2);
        e.mz = mz;
        t.mz = mz + 9e-6 * mz;
    }
    let abs = MorpheusScore {
        tolerance: Tolerance::Absolute(0.1),
    }
    .compute(&experimental, &theoretical)
    .unwrap();
    assert_eq!(abs.matches, 4);
    near(f64::from(abs.score), 4.1212, 1e-4);
    assert_eq!(
        MorpheusScore {
            tolerance: Tolerance::Ppm(10.0)
        }
        .compute(&experimental, &theoretical)
        .unwrap()
        .score,
        34.0
    );
    near(
        HyperScore {
            tolerance: Tolerance::Absolute(0.1),
        }
        .compute(&experimental, &theoretical)
        .unwrap(),
        3.401197,
        1e-6,
    );
    near(
        HyperScore {
            tolerance: Tolerance::Ppm(10.0),
        }
        .compute(&experimental, &theoretical)
        .unwrap(),
        67.8210771,
        1e-7,
    );
}

#[test]
fn morpheus_uses_two_reuse_passes_and_theoretical_count_for_errors() {
    let exp = spectrum(&[(100.0, 1.0), (100.05, 3.0), (110.0, 4.0)], &[]);
    let theoretical = spectrum(&[(100.0, 1.0), (100.02, 1.0)], &[]);
    let result = MorpheusScore {
        tolerance: Tolerance::Absolute(0.1),
    }
    .compute(&exp, &theoretical)
    .unwrap();
    assert_eq!(result.matches, 2);
    assert_eq!(result.n_peaks, 2);
    assert_eq!(result.total_intensity, 8.0);
    assert_eq!(result.matched_intensity, 4.0);
    assert_eq!(result.score, 2.5);
    near(f64::from(result.mean_error_da), 0.025, 1e-8);
    near(f64::from(result.mean_error_ppm), 250.0, 1e-4);
}

#[test]
fn morpheus_charge_mismatch_advances_source_pointer_without_retry() {
    let exp = spectrum(&[(100.0, 1.0), (100.01, 1.0)], &[]);
    let theoretical = spectrum(&[(100.0, 1.0)], &[]);
    let result = MorpheusScore {
        tolerance: Tolerance::Absolute(0.1),
    }
    .compute_with_charges(&exp, &[2, 1], &theoretical, &[1])
    .unwrap();
    // First pass cannot retry the theoretical peak after the first wrong charge;
    // second pass can still count a later experimental peak's ion current.
    assert_eq!(
        (result.matches, result.matched_intensity, result.score),
        (0, 1.0, 0.5)
    );
    assert_eq!(result.mean_error_da, 1e10);
}

#[test]
fn detailed_hyperscore_counts_all_terminal_series_and_source_error_denominator() {
    let theoretical = spectrum(
        &[
            (100.0, 1.0),
            (200.0, 1.0),
            (300.0, 1.0),
            (400.0, 1.0),
            (500.0, 1.0),
        ],
        &["a2+", "c2+", "x2+", "xl$z2+", "precursor"],
    );
    let exp = spectrum(
        &[
            (100.01, 2.0),
            (200.01, 2.0),
            (300.01, 2.0),
            (400.01, 2.0),
            (500.01, 2.0),
        ],
        &[],
    );
    let scorer = HyperScore {
        tolerance: Tolerance::Absolute(0.1),
    };
    near(
        scorer.compute(&exp, &theoretical).unwrap(),
        11.0_f64.ln(),
        1e-12,
    );
    let detail = scorer.compute_with_detail(&exp, &theoretical).unwrap();
    assert_eq!(
        (detail.matched_prefix_ions, detail.matched_suffix_ions),
        (2, 2)
    );
    near(detail.score, 11.0_f64.ln() + 2.0 * 2.0_f64.ln(), 1e-12);
    near(detail.mean_error, 0.05 / 4.0, 1e-12); // Error includes the matched unclassified precursor.
}

#[test]
fn hyperscore_boundaries_charge_awareness_and_crosslink_names() {
    let theoretical = spectrum(
        &[(100.0, 2.0), (200.0, 3.0), (300.0, 4.0)],
        &["xl$b2+", "xl$y2+", "a2+"],
    );
    let exp = spectrum(&[(100.5, 1.0), (200.5, 1.0), (300.5, 1.0)], &[]);
    let scorer = HyperScore {
        tolerance: Tolerance::Absolute(0.5),
    };
    near(
        scorer.compute(&exp, &theoretical).unwrap(),
        10.0_f64.ln(),
        1e-12,
    ); // Inclusive basic boundary, all peaks in dot.
    assert_eq!(
        scorer
            .compute_with_charges(&exp, &[1, 1, 1], &theoretical, &[1, 1, 1])
            .unwrap(),
        0.0
    ); // Strict charge-aware boundary.
    let scorer = HyperScore {
        tolerance: Tolerance::Absolute(0.6),
    };
    near(
        scorer
            .compute_with_charges(&exp, &[1, 1, 1], &theoretical, &[1, 1, 1])
            .unwrap(),
        6.0_f64.ln(),
        1e-12,
    ); // b/y only.
    near(
        scorer
            .compute_with_charges(&exp, &[2, 1, 1], &theoretical, &[1, 1, 1])
            .unwrap(),
        4.0_f64.ln(),
        1e-12,
    );
}

#[test]
fn charge_aware_intensity_support_collapses_ordinals_and_is_atomic() {
    let theoretical = spectrum(
        &[(100.0, 1.0), (200.0, 1.0), (300.0, 1.0), (400.0, 1.0)],
        &["b1+", "b1++", "b2+", "y1+"],
    );
    let exp = spectrum(
        &[(100.0, 2.0), (200.0, 3.0), (300.0, 4.0), (400.0, 5.0)],
        &[],
    );
    let mut sums = [10.0; 4];
    let score = HyperScore::default()
        .compute_with_intensity_sum(&exp, &[1; 4], &theoretical, &[1; 4], &mut sums)
        .unwrap();
    assert_eq!(sums, [15.0, 14.0, 10.0, 15.0]);
    near(score, 15.0_f64.ln() + 2.0_f64.ln(), 1e-12); // Two distinct prefixes, one suffix.
    let mut bad_theoretical = theoretical.clone();
    bad_theoretical.string_data_arrays[0].data[3] = "y99+".into();
    let before = sums;
    assert!(
        HyperScore::default()
            .compute_with_intensity_sum(&exp, &[1; 4], &bad_theoretical, &[1; 4], &mut sums)
            .is_err()
    );
    assert_eq!(sums, before);
}

#[test]
fn scoring_rejects_bad_annotations_arrays_and_nonfinite_or_zero_tic_results() {
    let exp = spectrum(&[(100.0, 1.0)], &[]);
    let theoretical = spectrum(&[(100.0, 1.0)], &["b1+"]);
    assert!(HyperScore::default().compute(&exp, &exp).is_err());
    assert!(
        HyperScore::default()
            .compute(&exp, &spectrum(&[(100.0, 1.0)], &[""]))
            .is_err()
    );
    assert!(
        HyperScore::default()
            .compute_with_charges(&exp, &[], &theoretical, &[1])
            .is_err()
    );
    assert!(
        MorpheusScore::default()
            .compute_with_charges(&exp, &[1], &theoretical, &[])
            .is_err()
    );
    assert!(
        MorpheusScore::default()
            .compute(&spectrum(&[(100.0, 0.0)], &[]), &theoretical)
            .is_err()
    );
    let huge = spectrum(&[(100.0, f32::MAX)], &["b1+"]);
    assert!(HyperScore::default().compute(&huge, &huge).is_err()); // f32 product overflow in the basic source overload.
    assert!(
        HyperScore::default()
            .compute_with_charges(&huge, &[1], &huge, &[1])
            .unwrap()
            .is_finite()
    );
    assert!(
        HyperScore {
            tolerance: Tolerance::Absolute(-1.0)
        }
        .compute(&exp, &theoretical)
        .is_err()
    );
    assert!(
        MorpheusScore::default()
            .compute(&spectrum(&[(0.0, 1.0)], &[]), &spectrum(&[(0.0, 1.0)], &[]))
            .is_err()
    );
}
