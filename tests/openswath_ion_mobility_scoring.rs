// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! IonMobilityScoring.
//!
//! Expected values are transcribed from the upstream class test at the core pin
//! `bc9cc12` (evidence tier 3), `src/tests/class_tests/openms/source/
//! IonMobilityScoring_test.cpp`, cited per test by section name. Tests marked
//! "derived" check an independently derived invariant or a native guard
//! (tier 4).

use openms::Error;
use openms::analysis::openswath::data_structures::{BinaryDataArray, Spectrum, SpectrumPtr};
use openms::analysis::openswath::ion_mobility_scoring::{
    DriftExtraction, IonMobilityScoring, sum_aligned_mobilograms,
};
use openms::analysis::openswath::transition_experiment::LightTransition;
use openms::kernel::mrm::OpenSwathScores;
use openms::kernel::ranges::RangeBase;
use openms::kernel::{MobilityPeak1D, Mobilogram};
use std::sync::Arc;

fn real_similar(a: f64, b: f64) {
    let absdiff = (a - b).abs();
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{a} !~ {b} (absdiff {absdiff}, ratio {ratio})"
    );
}

fn range(center: f64, width: f64) -> RangeBase {
    let mut r = RangeBase::singular(center).unwrap();
    r.min_span_if_singular(width).unwrap();
    r
}

fn transitions() -> Vec<LightTransition> {
    [(500.2, "group1"), (600.5, "group2")]
        .iter()
        .map(|(mz, name)| LightTransition {
            product_mz: *mz,
            precursor_mz: 700.2,
            fragment_charge: 1,
            transition_name: (*name).into(),
            ..LightTransition::default()
        })
        .collect()
}

fn with_im(mz: Vec<f64>, intensity: Vec<f64>, im: Vec<f64>, described: bool) -> SpectrumPtr {
    let mut s = Spectrum::from_arrays(mz, intensity);
    let desc = if described { "Ion Mobility" } else { "" };
    s.push_data_array(BinaryDataArray::with_description(im, desc));
    Arc::new(s)
}

fn spec() -> SpectrumPtr {
    let mut mz = vec![500.2; 8];
    mz.extend(vec![500.3; 8]);
    mz.extend([600.2, 600.3, 600.4, 600.5, 600.6, 600.7, 600.8, 600.9]);
    let mut int = vec![10., 20., 30., 40., 40., 30., 20., 10.];
    int.extend([10., 20., 30., 40., 40., 30., 20., 10.]);
    int.extend(vec![20.; 8]);
    let im = vec![
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.3, 0.1, 0.2,
        0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
    ];
    with_im(mz, int, im, true)
}

fn ms1spec() -> SpectrumPtr {
    with_im(
        vec![700.2, 700.3, 700.4, 700.5, 700.6, 700.7, 700.8, 700.9],
        vec![10., 20., 30., 40., 40., 30., 20., 10.],
        vec![0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
        true,
    )
}

fn extraction(window: f64) -> DriftExtraction {
    DriftExtraction {
        dia_extraction_window: window,
        dia_extraction_ppm: false,
        drift_extra: 0.25,
    }
}

fn drift(scores: &mut OpenSwathScores, target: f64, im: RangeBase, window: f64) {
    IonMobilityScoring::drift_scoring(
        &[spec()],
        &transitions(),
        scores,
        target,
        im,
        &extraction(window),
        false,
        None,
    )
    .unwrap();
}

/// IonMobilityScoring_test.cpp:[EXTRA] driftScoring, tests #1-#6.
#[test]
fn drift_scoring() {
    let mut scores = OpenSwathScores::default();
    // #1: an ion-mobility array without description is not a drift array
    let empty = with_im(Vec::new(), Vec::new(), Vec::new(), false);
    IonMobilityScoring::drift_scoring(
        &[empty],
        &transitions(),
        &mut scores,
        1.0,
        range(1.0, 1.0),
        &extraction(0.3),
        false,
        None,
    )
    .unwrap();
    real_similar(scores.im_drift, 0.0);
    real_similar(scores.im_drift_weighted, 0.0);
    real_similar(scores.im_delta_score, 0.0);
    real_similar(scores.im_xcorr_shape_score, 0.0);
    real_similar(scores.im_xcorr_coelution_score, 0.0);
    // #2
    let s = spec();
    assert_eq!(s.mz_array().data.len(), 24);
    assert_eq!(s.drift_time_array().unwrap().data.len(), 24);
    drift(&mut scores, 1.0, range(1.0, 1.0), 0.3);
    real_similar(scores.im_drift, (0.705405 + 0.4) / 2.0);
    real_similar(scores.im_drift_weighted, 0.662790697674419);
    real_similar(scores.im_delta_score, (0.294595 + 0.6) / 2.0);
    real_similar(scores.im_xcorr_shape_score, 0.892124778448826);
    real_similar(scores.im_xcorr_coelution_score, 2.73205080756888);
    // #3
    drift(&mut scores, 1.0, range(1.0, 1.0), 0.1);
    real_similar(scores.im_drift, (0.5 + 0.4) / 2.0);
    real_similar(scores.im_drift_weighted, 0.489473684210526);
    real_similar(scores.im_delta_score, (0.5 + 0.6) / 2.0);
    real_similar(scores.im_xcorr_shape_score, 0.833513903989399);
    real_similar(scores.im_xcorr_coelution_score, 0.910683602522959);
    // #4: exactly one entry in the mobilogram
    drift(&mut scores, 1.05, range(1.05, 0.1), 0.3);
    real_similar(scores.im_drift, 1.1);
    real_similar(scores.im_drift_weighted, 1.1);
    real_similar(scores.im_delta_score, 0.05);
    real_similar(scores.im_xcorr_shape_score, 0.0);
    real_similar(scores.im_xcorr_coelution_score, 1.0);
    // #5: one zero transition
    drift(
        &mut scores,
        1.1,
        RangeBase::from_min_max(1.0, 1.3).unwrap(),
        0.3,
    );
    real_similar(scores.im_drift, 1.16666666666667);
    real_similar(scores.im_drift_weighted, 1.16666666666667);
    real_similar(scores.im_delta_score, 0.0666666666666667);
    real_similar(scores.im_xcorr_shape_score, 1.0 / 3.0);
    real_similar(scores.im_xcorr_coelution_score, 3.73205080756888);
    // #6: all-zero transitions
    drift(&mut scores, 1.1, range(3.0, 1.0), 0.3);
    real_similar(scores.im_drift, -1.0);
    real_similar(scores.im_drift_weighted, -1.0);
    real_similar(scores.im_delta_score, -1.0);
    assert!(scores.im_xcorr_shape_score.is_nan());
    real_similar(scores.im_xcorr_coelution_score, 0.0);
}

/// Derived: peak picking needs the unported PeakPickerMobilogram and fails
/// without touching the scores.
#[test]
fn drift_scoring_peak_picking_is_unsupported() {
    let mut scores = OpenSwathScores::default();
    let before = scores;
    let r = IonMobilityScoring::drift_scoring(
        &[spec()],
        &transitions(),
        &mut scores,
        1.0,
        range(1.0, 1.0),
        &extraction(0.3),
        true,
        None,
    );
    assert!(matches!(r, Err(Error::Unsupported(_))));
    assert_eq!(scores, before);
}

/// IonMobilityScoring_test.cpp:[EXTRA] driftScoringMS1.
#[test]
fn drift_scoring_ms1() {
    let mut scores = OpenSwathScores::default();
    let run = |spectra: &[SpectrumPtr], scores: &mut OpenSwathScores| {
        IonMobilityScoring::drift_scoring_ms1(
            spectra,
            &transitions(),
            scores,
            1.0,
            range(1.0, 1.0),
            &extraction(0.3),
            None,
        )
        .unwrap();
    };
    // spectra without a described drift array are skipped
    run(
        &[with_im(Vec::new(), Vec::new(), Vec::new(), false)],
        &mut scores,
    );
    let s = ms1spec();
    assert_eq!(s.mz_array().data.len(), 8);
    assert_eq!(s.drift_time_array().unwrap().data.len(), 8);
    run(&[s], &mut scores);
    real_similar(scores.im_ms1_delta_score, 0.7);
}

fn contrast(scores: &mut OpenSwathScores, im: RangeBase, window: f64) {
    IonMobilityScoring::drift_scoring_ms1_contrast(
        &[spec()],
        &[ms1spec()],
        &transitions(),
        scores,
        im,
        &extraction(window),
        None,
    )
    .unwrap();
}

/// IonMobilityScoring_test.cpp:[EXTRA] driftScoringMS1Contrast.
#[test]
fn drift_scoring_ms1_contrast() {
    let mut scores = OpenSwathScores::default();
    contrast(&mut scores, range(1.0, 1.0), 0.3);
    real_similar(scores.im_ms1_contrast_coelution, 5.62132034355964);
    real_similar(scores.im_ms1_contrast_shape, 0.50991093654836);
    real_similar(scores.im_ms1_sum_contrast_coelution, 2.0);
    real_similar(scores.im_ms1_sum_contrast_shape, 0.56486260935015);
    contrast(&mut scores, range(1.0, 1.0), 0.1);
    real_similar(scores.im_ms1_contrast_coelution, 6.0);
    real_similar(scores.im_ms1_contrast_shape, 0.0);
    real_similar(scores.im_ms1_sum_contrast_coelution, 6.0);
    real_similar(scores.im_ms1_sum_contrast_shape, 0.0);
    contrast(&mut scores, range(1.05, 0.1), 0.3);
    real_similar(scores.im_ms1_contrast_coelution, 1.0);
    real_similar(scores.im_ms1_contrast_shape, 0.0);
    real_similar(scores.im_ms1_sum_contrast_coelution, 1.0);
    real_similar(scores.im_ms1_sum_contrast_shape, 0.0);
    contrast(&mut scores, range(1.15, 0.3), 0.3);
    real_similar(scores.im_ms1_contrast_coelution, 3.0);
    real_similar(scores.im_ms1_contrast_shape, 0.0);
    real_similar(scores.im_ms1_sum_contrast_coelution, 3.0);
    real_similar(scores.im_ms1_sum_contrast_shape, 0.0);
    contrast(&mut scores, range(3.0, 1.0), 0.3);
    real_similar(scores.im_ms1_contrast_coelution, 0.0);
    assert!(scores.im_ms1_contrast_shape.is_nan());
    real_similar(scores.im_ms1_sum_contrast_coelution, 0.0);
    assert!(scores.im_ms1_sum_contrast_shape.is_nan());
}

/// Derived: the grid merges values closer than `eps`, alignment inserts zeros,
/// and summing adds point by point (IonMobilityScoring.cpp `computeGrid_`,
/// `alignToGrid_`, `sumAlignedMobilograms`).
#[test]
fn grid_align_and_sum() {
    let m = |pts: &[(f64, f32)]| {
        Mobilogram::from_peaks(
            pts.iter()
                .map(|&(mobility, intensity)| MobilityPeak1D {
                    mobility,
                    intensity,
                })
                .collect(),
        )
    };
    let a = m(&[(0.1, 1.0), (0.3, 3.0)]);
    let b = m(&[(0.1000000001, 2.0), (0.2, 5.0)]);
    let grid = IonMobilityScoring::compute_grid(&[a.clone(), b.clone()], 1e-6);
    assert_eq!(grid, vec![0.1, 0.2, 0.3]);
    let (aligned_a, _) = IonMobilityScoring::align_to_grid(&a, &grid, 1e-6).unwrap();
    let ints: Vec<f32> = aligned_a.peaks.iter().map(|p| p.intensity).collect();
    assert_eq!(ints, vec![1.0, 0.0, 3.0]);
    let (aligned_b, _) = IonMobilityScoring::align_to_grid(&b, &grid, 1e-6).unwrap();
    let sum = sum_aligned_mobilograms(&[aligned_a, aligned_b.clone()]).unwrap();
    let ints: Vec<f32> = sum.peaks.iter().map(|p| p.intensity).collect();
    assert_eq!(ints, vec![3.0, 5.0, 3.0]);
    assert!(sum_aligned_mobilograms(&[aligned_b, a]).is_err());
    assert!(sum_aligned_mobilograms(&[]).unwrap().peaks.is_empty());
}

/// Derived: drift_id_scoring with an identification transition that has no
/// signal sets the four contrast scores to 0 (IonMobilityScoring.cpp
/// `driftIdScoring`, the `identification_im == -1` branch).
#[test]
fn drift_id_scoring_without_signal() {
    let mut scores = OpenSwathScores {
        im_ind_contrast_coelution: 9.0,
        ..OpenSwathScores::default()
    };
    let id = [LightTransition {
        product_mz: 900.0,
        transition_name: "id".into(),
        ..LightTransition::default()
    }];
    IonMobilityScoring::drift_id_scoring(
        &[spec()],
        &id,
        &transitions(),
        &mut scores,
        1.0,
        range(1.0, 1.0),
        &extraction(0.3),
        true,
        None,
    )
    .unwrap();
    assert_eq!(scores.im_drift, -1.0);
    real_similar(scores.im_delta, 2.0);
    assert_eq!(scores.im_ind_contrast_coelution, 0.0);
    assert_eq!(scores.im_log_intensity, 0.0);
}
