// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! MRMScoring.
//!
//! Expected values are transcribed from the upstream class test at the core pin
//! `bc9cc12` (evidence tier 3), `src/tests/class_tests/openms/source/
//! MRMScoring_test.cpp`, cited per test by section name. Tests marked
//! "derived" check an independently derived invariant (tier 4).

use openms::analysis::openswath::itransition::SignalToNoisePtr;
use openms::analysis::openswath::mock_objects::{MockFeature, MockMRMFeature, MockSignalToNoise};
use openms::analysis::openswath::mrm_scoring::MRMScoring;
use openms::analysis::openswath::transition_experiment::{LightCompound, LightTransition};
use std::collections::BTreeMap;
use std::sync::Arc;

fn real_similar(a: f64, b: f64) {
    let absdiff = (a - b).abs();
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{a} !~ {b} (absdiff {absdiff}, ratio {ratio})"
    );
}

const INTENSITY1: [f64; 11] = [
    5.97543668746948,
    4.2749171257019,
    3.3301842212677,
    4.08597040176392,
    5.50307035446167,
    5.24326848983765,
    8.40812492370605,
    2.83419919013977,
    6.94378805160522,
    7.69957494735718,
    4.08597040176392,
];
const INTENSITY2: [f64; 11] = [
    15.8951349258423,
    41.5446395874023,
    76.0746307373047,
    109.069435119629,
    111.90364074707,
    169.79216003418,
    121.043930053711,
    63.0136985778809,
    44.6150207519531,
    21.4926776885986,
    7.93575811386108,
];
const MS1: [f64; 11] = [
    0.0, 110.0, 200.0, 270.0, 320.0, 350.0, 360.0, 350.0, 320.0, 270.0, 200.0,
];

fn trace(values: &[f64]) -> Arc<MockFeature> {
    Arc::new(MockFeature {
        intensity_vec: values.to_vec(),
        ..MockFeature::default()
    })
}

fn ids(names: &[&str]) -> Vec<String> {
    names.iter().map(|s| (*s).to_string()).collect()
}

/// MRMScoring_test.cpp `fill_mock_objects`.
fn mock1() -> (MockMRMFeature, Vec<String>) {
    let mut features = BTreeMap::new();
    features.insert("group1".to_string(), trace(&INTENSITY1));
    features.insert("group2".to_string(), trace(&INTENSITY2));
    let mut precursors = BTreeMap::new();
    precursors.insert("ms1trace".to_string(), trace(&MS1));
    let f = MockMRMFeature {
        features,
        precursor_features: precursors,
        ..MockMRMFeature::default()
    };
    (f, ids(&["group1", "group2"]))
}

/// MRMScoring_test.cpp `fill_mock_objects2`.
fn mock2() -> (MockMRMFeature, Vec<String>, Vec<String>) {
    let (mut f, native) = mock1();
    let mut precursors = BTreeMap::new();
    precursors.insert("ms1trace1".to_string(), trace(&MS1));
    precursors.insert(
        "ms1trace2".to_string(),
        trace(&[
            10.0, 115.0, 180.0, 280.0, 330.0, 340.0, 390.0, 320.0, 300.0, 250.0, 100.0,
        ]),
    );
    precursors.insert("ms1trace3".to_string(), trace(&[0.0; 11]));
    f.precursor_features = precursors;
    (f, ids(&["ms1trace1", "ms1trace2", "ms1trace3"]), native)
}

fn abs_sums(m: &openms::analysis::openswath::mrm_scoring::XCorrMatrixType) -> Vec<f64> {
    // Row-major order, as the class test's nested loop (r outer, c inner).
    let mut out = Vec::new();
    for r in 0..m.rows() {
        for c in 0..m.cols() {
            out.push(m.get(r, c).unwrap().data.iter().map(|e| e.1.abs()).sum());
        }
    }
    out
}

fn check_auto_and_cross(m: &openms::analysis::openswath::mrm_scoring::XCorrMatrixType) {
    assert_eq!(m.rows(), 2);
    assert_eq!(m.cols(), 2);
    let auto = m.get(0, 0).unwrap();
    assert_eq!(auto.data.len(), 21);
    assert_eq!(auto.data[10].0, 0);
    assert_eq!(auto.data[11].0, 1);
    assert_eq!(auto.data[9].0, -1);
    assert_eq!(auto.data[12].0, 2);
    assert_eq!(auto.data[8].0, -2);
    real_similar(auto.data[10].1, 1.0);
    real_similar(auto.data[11].1, -0.227352707759245);
    real_similar(auto.data[9].1, -0.227352707759245);
    real_similar(auto.data[12].1, -0.07501116);
    real_similar(auto.data[8].1, -0.07501116);
    let cross = m.get(0, 1).unwrap();
    real_similar(cross.data[12].1, -0.31165141);
    real_similar(cross.data[11].1, -0.35036919);
    real_similar(cross.data[10].1, 0.03129565);
    real_similar(cross.data[9].1, 0.30204049);
    real_similar(cross.data[8].1, 0.13012441);
    real_similar(cross.data[7].1, 0.39698322);
    real_similar(cross.data[6].1, 0.16608774);
}

/// MRMScoring_test.cpp:initializeXCorrMatrix.
#[test]
fn initialize_xcorr_matrix() {
    let (f, native) = mock1();
    let mut s = MRMScoring::new();
    s.initialize_xcorr_matrix(&f, &native).unwrap();
    check_auto_and_cross(s.xcorr_matrix());
    // derived: the raw-data overload gives the same matrix
    let mut t = MRMScoring::new();
    t.initialize_xcorr_matrix_from_data(&[INTENSITY1.to_vec(), INTENSITY2.to_vec()])
        .unwrap();
    assert_eq!(t.xcorr_matrix(), s.xcorr_matrix());
}

/// MRMScoring_test.cpp:initializeXCorrContrastMatrix.
#[test]
fn initialize_xcorr_contrast_matrix() {
    let (f, native) = mock1();
    let mut s = MRMScoring::new();
    s.initialize_xcorr_contrast_matrix(&f, &native, &native)
        .unwrap();
    check_auto_and_cross(s.xcorr_contrast_matrix());
}

/// MRMScoring_test.cpp:initializeXCorrPrecursorContrastMatrix.
#[test]
fn initialize_xcorr_precursor_contrast_matrix() {
    let (f, pre, native) = mock2();
    let mut s = MRMScoring::new();
    s.initialize_xcorr_precursor_contrast_matrix(&f, &pre, &native)
        .unwrap();
    let m = s.xcorr_precursor_contrast_matrix();
    assert_eq!(m.rows(), 3);
    assert_eq!(m.cols(), 2);
    let sums = abs_sums(m);
    for (got, exp) in sums
        .iter()
        .zip([3.40949220, 6.19794611, 3.68912454, 6.60757921, 0.0, 0.0])
    {
        real_similar(*got, exp);
    }
}

/// MRMScoring_test.cpp:initializeXCorrPrecursorCombinedMatrix.
#[test]
fn initialize_xcorr_precursor_combined_matrix() {
    let (f, pre, native) = mock2();
    let mut s = MRMScoring::new();
    s.initialize_xcorr_precursor_combined_matrix(&f, &pre, &native)
        .unwrap();
    let m = s.xcorr_precursor_combined_matrix();
    assert_eq!(m.rows(), 5);
    assert_eq!(m.cols(), 5);
    let sums = abs_sums(m);
    for (i, exp) in [
        (0, 5.86440677),
        (1, 6.05410398),
        (2, 0.0),
        (3, 3.40949220),
        (4, 6.19794611),
        (6, 6.30751744),
        (7, 0.0),
        (8, 3.68912454),
        (9, 6.60757921),
        (12, 0.0),
        (13, 0.0),
        (14, 0.0),
        (18, 3.13711983),
        (19, 3.57832717),
        (24, 6.78303987),
    ] {
        real_similar(sums[i], exp);
    }
}

/// MRMScoring_test.cpp:test_calcXcorrCoelutionScore, test_calcXcorrShapeScore
/// and the two weighted variants.
#[test]
fn xcorr_coelution_and_shape() {
    let (f, native) = mock1();
    let mut s = MRMScoring::new();
    s.initialize_xcorr_matrix(&f, &native).unwrap();
    real_similar(s.calc_xcorr_coelution_score(), 1.0 + 3.0_f64.sqrt());
    real_similar(
        s.calc_xcorr_coelution_weighted_score(&[0.5, 0.5]).unwrap(),
        1.5,
    );
    real_similar(s.calc_xcorr_shape_score(), (1.0 + 0.3969832 + 1.0) / 3.0);
    real_similar(
        s.calc_xcorr_shape_weighted_score(&[0.5, 0.5]).unwrap(),
        0.6984916,
    );
    // derived: mismatched weights are refused
    assert!(s.calc_xcorr_shape_weighted_score(&[1.0]).is_err());
}

/// MRMScoring_test.cpp:test_calcSeparateXcorrContrastCoelutionScore and
/// test_calcSeparateXcorrContrastShapeScore.
#[test]
fn separate_xcorr_contrast() {
    let (f, native) = mock1();
    let mut s = MRMScoring::new();
    s.initialize_xcorr_contrast_matrix(&f, &native, &native)
        .unwrap();
    let c = s.calc_separate_xcorr_contrast_coelution_score();
    real_similar(c[0], 1.5);
    real_similar(c[1], 1.5);
    let sh = s.calc_separate_xcorr_contrast_shape_score();
    real_similar(sh[0], 0.698492);
    real_similar(sh[1], 0.698492);
}

/// MRMScoring_test.cpp:test_single_transition_xcorr_scores.
#[test]
fn single_transition_xcorr() {
    let (f, mut native) = mock1();
    native.truncate(1);
    let mut s = MRMScoring::new();
    s.initialize_xcorr_matrix(&f, &native).unwrap();
    assert_eq!(s.xcorr_matrix().rows(), 1);
    assert_eq!(s.xcorr_matrix().cols(), 1);
    real_similar(s.calc_xcorr_coelution_score(), 0.0);
    real_similar(s.calc_xcorr_coelution_weighted_score(&[1.0]).unwrap(), 0.0);
    real_similar(s.calc_xcorr_shape_score(), 1.0);
    real_similar(s.calc_xcorr_shape_weighted_score(&[1.0]).unwrap(), 1.0);
}

/// MRMScoring_test.cpp:calcXcorrPrecursorContrastCoelutionScore,
/// calcXcorrPrecursorCombinedCoelutionScore, calcXcorrPrecursorContrastShapeScore
/// and calcXcorrPrecursorCombinedShapeScore.
#[test]
fn xcorr_precursor_scores() {
    let (f, pre, native) = mock2();
    let mut s = MRMScoring::new();
    s.initialize_xcorr_precursor_contrast_matrix(&f, &pre, &native)
        .unwrap();
    real_similar(s.calc_xcorr_precursor_contrast_coelution_score(), 8.7328638);
    real_similar(s.calc_xcorr_precursor_contrast_shape_score(), 0.3772868);
    // derived: the SumFrag variants compute the same reductions
    assert_eq!(
        s.calc_xcorr_precursor_contrast_sum_frag_coelution_score(),
        s.calc_xcorr_precursor_contrast_coelution_score()
    );
    assert_eq!(
        s.calc_xcorr_precursor_contrast_sum_frag_shape_score(),
        s.calc_xcorr_precursor_contrast_shape_score()
    );
    s.initialize_xcorr_precursor_combined_matrix(&f, &pre, &native)
        .unwrap();
    real_similar(s.calc_xcorr_precursor_combined_coelution_score(), 8.4321855);
    real_similar(s.calc_xcorr_precursor_combined_shape_score(), 0.5079334);
}

/// MRMScoring_test.cpp:test_Library_score.
#[test]
fn library_score() {
    let f1 = Arc::new(MockFeature {
        intensity: 782.38073_f32,
        ..MockFeature::default()
    });
    let f2 = Arc::new(MockFeature {
        intensity: 58.384506_f32,
        ..MockFeature::default()
    });
    let mut features = BTreeMap::new();
    features.insert("group1".to_string(), f1);
    features.insert("group2".to_string(), f2.clone());
    features.insert("group3".to_string(), f2);
    let f = MockMRMFeature {
        features,
        ..MockMRMFeature::default()
    };
    let transitions: Vec<LightTransition> =
        [(1.0, "group1"), (10000.0, "group2"), (2000.0, "group3")]
            .iter()
            .map(|(lib, name)| LightTransition {
                library_intensity: *lib,
                transition_name: (*name).into(),
                ..LightTransition::default()
            })
            .collect();
    let r = MRMScoring::calc_library_score(&f, &transitions).unwrap();
    real_similar(r.correlation, -0.654591316);
    real_similar(r.norm_manhattan, 0.5800337593);
    real_similar(r.manhattan, 1.279644714);
    real_similar(r.dotprod, 0.34514801);
    real_similar(r.spectral_angle, 1.483262);
    real_similar(r.rmsd, 0.6727226674);
}

/// MRMScoring_test.cpp:test_RT_score.
#[test]
fn rt_score() {
    let pep = LightCompound {
        rt: 100.0,
        ..LightCompound::default()
    };
    real_similar(MRMScoring::calc_rt_score(&pep, 100.0), 0.0);
    real_similar(MRMScoring::calc_rt_score(&pep, 0.0), 100.0);
    // derived: the source's "unset" sentinel scores 0
    let unset = LightCompound {
        rt: -1000.0,
        ..LightCompound::default()
    };
    assert_eq!(MRMScoring::calc_rt_score(&unset, 5.0), 0.0);
}

/// MRMScoring_test.cpp:test_SN_score.
#[test]
fn sn_score() {
    let estimators: Vec<SignalToNoisePtr> = vec![
        Arc::new(MockSignalToNoise { sn_value: 500.0 }),
        Arc::new(MockSignalToNoise { sn_value: 1500.0 }),
    ];
    let (f, _) = mock1();
    real_similar(MRMScoring::calc_sn_score(&f, &estimators), 1000.0);
    let sep = MRMScoring::calc_separate_sn_score(&f, &estimators);
    real_similar(sep[0], 6.21461);
    real_similar(sep[1], 7.31322);
    assert_eq!(MRMScoring::calc_sn_score(&f, &[]), 0.0);
}

/// MRMScoring_test.cpp:initializeMIMatrix, test_calcMIScore and
/// test_calcMIWeightedScore.
#[test]
fn mi_matrix_and_scores() {
    let (f, native) = mock1();
    let mut s = MRMScoring::new();
    s.initialize_mi_matrix(&f, &native).unwrap();
    let m = s.mi_matrix();
    assert_eq!(m.rows(), 2);
    assert_eq!(m.cols(), 2);
    real_similar(*m.get(0, 0).unwrap(), 3.2776);
    real_similar(*m.get(0, 1).unwrap(), 3.2776);
    real_similar(*m.get(1, 1).unwrap(), 3.4594);
    real_similar(*m.get(1, 0).unwrap(), 0.0);
    real_similar(s.calc_mi_score(), 3.3382);
    real_similar(s.calc_mi_weighted_score(&[0.5, 0.5]).unwrap(), 3.3231);
}

/// MRMScoring_test.cpp:initializeMIContrastMatrix and
/// test_calcSeparateMIContrastScore.
#[test]
fn mi_contrast() {
    let (f, native1) = mock1();
    let native2: Vec<String> = native1.iter().rev().cloned().collect();
    let mut s = MRMScoring::new();
    s.initialize_mi_contrast_matrix(&f, &native1, &native2)
        .unwrap();
    let m = s.mi_contrast_matrix();
    real_similar(*m.get(0, 0).unwrap(), 3.2776);
    real_similar(*m.get(0, 1).unwrap(), 3.2776);
    real_similar(*m.get(1, 1).unwrap(), 3.2776);
    real_similar(*m.get(1, 0).unwrap(), 3.4594);
    let sep = s.calc_separate_mi_contrast_score();
    real_similar(sep[0], 3.27761);
    real_similar(sep[1], 3.36852);
}

/// MRMScoring_test.cpp:test_single_transition_mi_scores.
#[test]
fn single_transition_mi() {
    let (f, mut native) = mock1();
    native.truncate(1);
    let mut s = MRMScoring::new();
    s.initialize_mi_matrix(&f, &native).unwrap();
    assert_eq!(s.mi_matrix().rows(), 1);
    real_similar(*s.mi_matrix().get(0, 0).unwrap(), 3.2776);
    real_similar(s.calc_mi_score(), 3.2776);
    real_similar(s.calc_mi_weighted_score(&[1.0]).unwrap(), 3.2776);
}

/// MRMScoring_test.cpp:initializeMIPrecursorContrastMatrix,
/// initializeMIPrecursorCombinedMatrix, test_calcMIPrecursorContrastScore and
/// test_calcMIPrecursorCombinedScore.
#[test]
fn mi_precursor() {
    let (f, pre, native) = mock2();
    let mut s = MRMScoring::new();
    s.initialize_mi_precursor_contrast_matrix(&f, &pre, &native)
        .unwrap();
    let m = s.mi_precursor_contrast_matrix();
    assert_eq!(m.rows(), 3);
    assert_eq!(m.cols(), 2);
    real_similar(m.data().iter().sum::<f64>(), 12.01954465);
    real_similar(s.calc_mi_precursor_contrast_score(), 2.003257);
    s.initialize_mi_precursor_combined_matrix(&f, &pre, &native)
        .unwrap();
    let m = s.mi_precursor_combined_matrix();
    assert_eq!(m.rows(), 5);
    assert_eq!(m.cols(), 5);
    real_similar(m.data().iter().sum::<f64>(), 48.98726953);
    real_similar(s.calc_mi_precursor_combined_score(), 1.959490);
}

/// Derived: an unknown native ID is an error, not a null dereference.
#[test]
fn unknown_native_id() {
    let (f, _) = mock1();
    let mut s = MRMScoring::new();
    assert!(s.initialize_xcorr_matrix(&f, &ids(&["nope"])).is_err());
}

// ---------------------------------------------------------------------------
// OpenSwathScores (OpenSwathScores_test.cpp)
// ---------------------------------------------------------------------------

use openms::analysis::openswath::openswath_scores::OpenSwathScoresUsage;
use openms::kernel::mrm::OpenSwathScores;

/// OpenSwathScores_test.cpp:EXTRA OpenSwath_Scores().
#[test]
fn scores_defaults() {
    let s = OpenSwathScores::default();
    real_similar(s.ms1_xcorr_coelution_score, -1.0);
    real_similar(s.ms1_xcorr_shape_score, -1.0);
    real_similar(s.ms1_mi_score, -1.0);
    // derived: OpenSwathScores.h `OpenSwath_Scores_Usage` member initialisers
    let u = OpenSwathScoresUsage::default();
    assert!(u.use_coelution_score && u.use_ms2_isotope_scores && u.use_ionseries_scores);
    assert!(!u.use_peak_shape_metrics);
}

/// OpenSwathScores_test.cpp:get_quick_lda_score.
#[test]
fn scores_quick_lda() {
    real_similar(
        OpenSwathScores::get_quick_lda_score(1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
        -0.5319046 + 2.1643962 + 8.0353047 + 0.1458914 + -1.6901925 + -0.8002824,
    );
}

fn ones() -> OpenSwathScores {
    OpenSwathScores {
        library_corr: 1.0,
        library_norm_manhattan: 1.0,
        norm_rt_score: 1.0,
        xcorr_coelution_score: 1.0,
        xcorr_shape_score: 1.0,
        log_sn_score: 1.0,
        elution_model_fit_score: 1.0,
        ..OpenSwathScores::default()
    }
}

/// OpenSwathScores_test.cpp:calculate_lda_prescore.
#[test]
fn scores_lda_prescore() {
    real_similar(OpenSwathScores::default().calculate_lda_prescore(), 0.0);
    real_similar(
        ones().calculate_lda_prescore(),
        -0.34664267 + 2.98700722 + 7.05496384 + 0.09445371 + -5.71823862 + -0.72989582 + 1.88443209,
    );
}

/// OpenSwathScores_test.cpp:calculate_lda_single_transition.
#[test]
fn scores_lda_single_transition() {
    real_similar(
        OpenSwathScores::default().calculate_lda_single_transition(),
        0.0,
    );
    real_similar(
        ones().calculate_lda_single_transition(),
        7.05496384 + -0.72989582 + -1.08443209,
    );
}

/// OpenSwathScores_test.cpp:calculate_swath_lda_prescore.
#[test]
fn scores_swath_lda_prescore() {
    real_similar(
        OpenSwathScores::default().calculate_swath_lda_prescore(),
        0.0,
    );
    let s = OpenSwathScores {
        library_corr: 1.0,
        library_norm_manhattan: 1.0,
        norm_rt_score: 1.0,
        isotope_correlation: 1.0,
        isotope_overlap: 1.0,
        massdev_score: 1.0,
        xcorr_coelution_score: 1.0,
        xcorr_shape_score: 1.0,
        yseries_score: 1.0,
        log_sn_score: 1.0,
        ..OpenSwathScores::default()
    };
    real_similar(
        s.calculate_swath_lda_prescore(),
        -0.19011762
            + 2.47298914
            + 5.63906731
            + -0.62640133
            + 0.36006925
            + 0.08814003
            + 0.13978311
            + -1.16475032
            + -0.19267813
            + -0.61712054,
    );
}
