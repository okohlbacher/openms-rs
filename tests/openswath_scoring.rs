// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! OpenSwathScoring.
//!
//! Expected values are transcribed from the upstream class test at the core pin
//! `bc9cc12` (evidence tier 3), `src/tests/class_tests/openms/source/
//! OpenSwathScoring_test.cpp`, cited per test by section name. The source marks
//! the four `calculate*` sections NOT_TESTABLE (covered by
//! `MRMFeatureFinderScoring_test.cpp`); the tests here marked "derived" check
//! the driver against the scores it composes (tier 4).

use openms::Error;
use openms::analysis::openswath::data_structures::SpectrumPtr;
use openms::analysis::openswath::dia_scoring::DiaScoring;
use openms::analysis::openswath::mock_objects::{MockFeature, MockMRMFeature};
use openms::analysis::openswath::mrm_scoring::MRMScoring;
use openms::analysis::openswath::openswath_scores::OpenSwathScoresUsage;
use openms::analysis::openswath::openswath_scoring::OpenSwathScoring;
use openms::analysis::openswath::spectrum_access::SpectrumAccessPtr;
use openms::analysis::openswath::spectrum_access_openms::SimpleOpenMSSpectraFactory;
use openms::analysis::openswath::swath_map::SwathMap;
use openms::analysis::openswath::transition_experiment::{LightCompound, LightTransition};
use openms::kernel::mrm::{OpenSwathIndScores, OpenSwathScores};
use openms::kernel::ranges::RangeBase;
use openms::kernel::{DataArray, MSExperiment, MSSpectrum, Peak1D};
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

fn access(spectra: Vec<MSSpectrum>) -> SpectrumAccessPtr {
    let mut exp = MSExperiment::default();
    exp.spectra = spectra;
    SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(Arc::new(exp)).unwrap()
}

fn spectrum(peaks: &[(f64, f32)], rt: f64) -> MSSpectrum {
    MSSpectrum {
        peaks: peaks.iter().map(|&(mz, i)| Peak1D::new(mz, i)).collect(),
        rt,
        ..MSSpectrum::default()
    }
}

/// OpenSwathScoring_test.cpp `generateImSpec`.
fn im_spec(k_min: i32, k_max: i32, rt: f64) -> MSSpectrum {
    let mut s = spectrum(
        &(k_min..k_max)
            .map(|k| (100.0 + f64::from(k), 1.0))
            .collect::<Vec<_>>(),
        rt,
    );
    s.float_data_arrays.push(DataArray::new(
        "Ion Mobility",
        (k_min..k_max).map(|k| k as f32).collect(),
    ));
    s
}

fn scorer(method: &str) -> OpenSwathScoring {
    let mut sc = OpenSwathScoring::new();
    sc.initialize(
        1.0,
        1,
        0.005,
        0.20,
        0.0,
        &OpenSwathScoresUsage::default(),
        method,
        "fixed",
        true,
        true,
    )
    .unwrap();
    sc
}

fn mz(s: &SpectrumPtr) -> &[f64] {
    &s.mz_array().data
}
fn int(s: &SpectrumPtr) -> &[f64] {
    &s.intensity_array().data
}
fn im(s: &SpectrumPtr) -> &[f64] {
    &s.drift_time_array().unwrap().data
}

fn check(values: &[f64], expected: &[f64]) {
    assert_eq!(values.len(), expected.len(), "{values:?} vs {expected:?}");
    for (v, e) in values.iter().zip(expected) {
        real_similar(*v, *e);
    }
}

/// OpenSwathScoring_test.cpp:initialize.
#[test]
fn initialize() {
    let mut sc = OpenSwathScoring::new();
    let su = OpenSwathScoresUsage::default();
    sc.initialize(
        100.0, 1, 0.01, 0.20, 0.0, &su, "simple", "fixed", true, true,
    )
    .unwrap();
    // derived: the method strings are validated, as the source throws
    let before = sc;
    assert!(matches!(
        sc.initialize(1.0, 1, 0.01, 0.2, 0.0, &su, "sum", "fixed", true, true),
        Err(Error::InvalidValue(_))
    ));
    assert!(
        sc.initialize(1.0, 1, 0.01, 0.2, 0.0, &su, "simple", "static", true, true)
            .is_err()
    );
    assert_eq!(sc, before);
}

/// OpenSwathScoring_test.cpp:fetchSpectrumSwath, without ion mobility.
#[test]
fn fetch_spectrum_swath() {
    let none = RangeBase::new();
    // empty map, default (simple) addition
    let sp = OpenSwathScoring::new()
        .fetch_spectrum_swath(access(Vec::new()).as_ref(), 20.0, 1, &none)
        .unwrap();
    assert!(sp.is_empty());

    let one = access(vec![spectrum(&[(20.0, 200.0)], 20.0)]);
    assert_eq!(one.nr_spectra(), 1);
    for method in ["resample", "simple"] {
        let sp = scorer(method)
            .fetch_spectrum_swath(one.as_ref(), 20.0, 1, &none)
            .unwrap();
        assert_eq!(sp.len(), 1);
        check(mz(&sp[0]), &[20.0]);
        check(int(&sp[0]), &[200.0]);
    }

    let three = access(vec![
        spectrum(&[(20.0, 200.0)], 10.0),
        spectrum(&[(20.0, 200.0)], 20.0),
        spectrum(&[(20.0, 200.0)], 30.0),
    ]);
    assert_eq!(three.nr_spectra(), 3);
    let sp = scorer("resample")
        .fetch_spectrum_swath(three.as_ref(), 20.0, 3, &none)
        .unwrap();
    assert_eq!(sp.len(), 1);
    check(mz(&sp[0]), &[20.0]);
    check(int(&sp[0]), &[600.0]);
    let sp = scorer("simple")
        .fetch_spectrum_swath(three.as_ref(), 20.0, 3, &none)
        .unwrap();
    assert_eq!(sp.len(), 3);
    for s in &sp {
        check(mz(s), &[20.0]);
        check(int(s), &[200.0]);
    }

    let four = access(vec![
        spectrum(&[(20.0, 200.0)], 10.0),
        spectrum(&[(20.001, 200.0)], 20.0),
        spectrum(&[(250.001, 300.0)], 50.0),
        spectrum(&[(250.002, 500.0)], 60.0),
    ]);
    assert_eq!(four.nr_spectra(), 4);
    let sp = scorer("resample")
        .fetch_spectrum_swath(four.as_ref(), 20.0, 3, &none)
        .unwrap();
    assert_eq!(sp.len(), 1);
    check(mz(&sp[0]), &[20.0, 20.005, 250.0]);
    check(int(&sp[0]), &[360.0, 40.0, 300.0]);
    let sp = scorer("simple")
        .fetch_spectrum_swath(four.as_ref(), 20.0, 3, &none)
        .unwrap();
    assert_eq!(sp.len(), 3);
    check(mz(&sp[0]), &[20.001]);
    check(int(&sp[0]), &[200.0]);
    check(mz(&sp[1]), &[20.0]);
    check(int(&sp[1]), &[200.0]);
    check(mz(&sp[2]), &[250.0]);
    check(int(&sp[2]), &[300.0]);
}

/// OpenSwathScoring_test.cpp:fetchSpectrumSwath - extra (ion mobility).
#[test]
fn fetch_spectrum_swath_im() {
    let mut im_range = RangeBase::singular(3.0).unwrap();
    im_range.min_span_if_singular(2.0).unwrap();

    let one = access(vec![im_spec(1, 6, 20.0)]);
    assert_eq!(one.nr_spectra(), 1);
    let sp = scorer("resample")
        .fetch_spectrum_swath(one.as_ref(), 20.0, 1, &im_range)
        .unwrap();
    assert_eq!(sp.len(), 1);
    check(mz(&sp[0]), &[102.0, 103.0, 104.0]);
    check(int(&sp[0]), &[1.0, 1.0, 1.0]);
    check(im(&sp[0]), &[2.0, 3.0, 4.0]);
    let sp = scorer("simple")
        .fetch_spectrum_swath(one.as_ref(), 20.0, 1, &im_range)
        .unwrap();
    assert_eq!(sp.len(), 1);
    check(mz(&sp[0]), &[101.0, 102.0, 103.0, 104.0, 105.0]);
    check(int(&sp[0]), &[1.0; 5]);
    check(im(&sp[0]), &[1.0, 2.0, 3.0, 4.0, 5.0]);

    let spectra = || {
        vec![
            im_spec(1, 3, 19.0),
            im_spec(1, 6, 20.0),
            im_spec(3, 6, 21.0),
        ]
    };
    let with_tail = {
        let mut v = spectra();
        v.push(im_spec(1, 6, 250.0));
        v
    };
    for map in [access(spectra()), access(with_tail)] {
        let sp = scorer("resample")
            .fetch_spectrum_swath(map.as_ref(), 20.0, 3, &im_range)
            .unwrap();
        assert_eq!(sp.len(), 1);
        check(mz(&sp[0]), &[102.0, 102.0, 103.0, 103.0, 104.0, 104.0]);
        check(int(&sp[0]), &[1.0; 6]);
        check(im(&sp[0]), &[2.0, 2.0, 3.0, 3.0, 4.0, 4.0]);
        let sp = scorer("simple")
            .fetch_spectrum_swath(map.as_ref(), 20.0, 3, &im_range)
            .unwrap();
        assert_eq!(sp.len(), 3);
        check(mz(&sp[0]), &[101.0, 102.0, 103.0, 104.0, 105.0]);
        check(int(&sp[0]), &[1.0; 5]);
        check(im(&sp[0]), &[1.0, 2.0, 3.0, 4.0, 5.0]);
        check(mz(&sp[1]), &[101.0, 102.0]);
        check(int(&sp[1]), &[1.0, 1.0]);
        check(im(&sp[1]), &[1.0, 2.0]);
        check(mz(&sp[2]), &[103.0, 104.0, 105.0]);
        check(im(&sp[2]), &[3.0, 4.0, 5.0]);
    }

    // resampling with drift data but no IM range concatenates everything
    let map = access(spectra());
    let sp = scorer("resample")
        .fetch_spectrum_swath(map.as_ref(), 20.0, 3, &RangeBase::new())
        .unwrap();
    assert_eq!(sp.len(), 1);
    assert_eq!(mz(&sp[0]).len(), 10);
    assert_eq!(int(&sp[0]).len(), 10);
    assert_eq!(im(&sp[0]).len(), 10);
    real_similar(mz(&sp[0])[0], 101.0);
    real_similar(mz(&sp[0])[1], 101.0);
    real_similar(mz(&sp[0])[8], 105.0);
    real_similar(mz(&sp[0])[9], 105.0);
}

/// Derived: the scoring form of the fetch returns nothing for an empty map in
/// both modes (OpenSwathScoring.cpp `fetchSpectrumSwath(..., out)` returns
/// early), and SONAR maps give one merged spectrum.
#[test]
fn fetch_spectrum_swath_into_and_maps() {
    let none = RangeBase::new();
    let empty = access(Vec::new());
    assert!(
        scorer("resample")
            .fetch_spectrum_swath_into(empty.as_ref(), 20.0, 1, &none)
            .unwrap()
            .is_empty()
    );
    let map = |mz: f64| SwathMap {
        sptr: Some(access(vec![spectrum(&[(mz, 100.0)], 20.0)])),
        ..SwathMap::default()
    };
    let sp = scorer("simple")
        .fetch_spectrum_swath_maps(&[map(20.0), map(30.0)], 20.0, 1, &none)
        .unwrap();
    assert_eq!(sp.len(), 1);
    check(mz(&sp[0]), &[20.0, 30.0]);
    check(int(&sp[0]), &[100.0, 100.0]);
    assert!(
        scorer("simple")
            .fetch_spectrum_swath_maps(&[], 20.0, 1, &none)
            .is_err()
    );
}

fn feature() -> MockMRMFeature {
    let trace = |v: &[f64], apex: f32| {
        Arc::new(MockFeature {
            intensity_vec: v.to_vec(),
            intensity: apex,
            rt: 20.0,
            ..MockFeature::default()
        })
    };
    let mut features = BTreeMap::new();
    features.insert("t1".to_string(), trace(&[1.0, 3.0, 6.0, 3.0, 1.0], 6.0));
    features.insert("t2".to_string(), trace(&[2.0, 5.0, 9.0, 4.0, 1.0], 9.0));
    MockMRMFeature {
        features,
        intensity: 15.0,
        rt: 20.0,
        ..MockMRMFeature::default()
    }
}

/// Derived: the chromatographic and ID scores equal the MRMScoring values they
/// are built from, and switched-off scores are left alone.
#[test]
fn chromatographic_scores_compose_mrm_scoring() {
    let f = feature();
    let ids = vec!["t1".to_string(), "t2".to_string()];
    let weights = [0.4, 0.6];
    let sc = OpenSwathScoring::new();
    let mut scores = OpenSwathScores::default();
    sc.calculate_chromatographic_scores(&f, &ids, &[], &weights, &[], &mut scores)
        .unwrap();
    let mut mrm = MRMScoring::new();
    mrm.initialize_xcorr_matrix(&f, &ids).unwrap();
    assert_eq!(
        scores.xcorr_coelution_score,
        mrm.calc_xcorr_coelution_score()
    );
    assert_eq!(scores.xcorr_shape_score, mrm.calc_xcorr_shape_score());
    assert_eq!(
        scores.weighted_xcorr_shape,
        mrm.calc_xcorr_shape_weighted_score(&weights).unwrap()
    );
    mrm.initialize_mi_matrix(&f, &ids).unwrap();
    assert_eq!(scores.mi_score, mrm.calc_mi_score());
    assert_eq!(scores.nr_peaks, 2.0);
    assert_eq!(scores.sn_ratio, 0.0);
    assert_eq!(scores.log_sn_score, 0.0);
    // no precursors: the MS1 scores keep their -1 "not computed" defaults
    assert_eq!(scores.ms1_xcorr_coelution_score, -1.0);

    let mut ind = OpenSwathIndScores::default();
    sc.calculate_chromatographic_id_scores(&f, &ids[..1], &ids, &[], &mut ind)
        .unwrap();
    let mut c = MRMScoring::new();
    c.initialize_xcorr_contrast_matrix(&f, &ids[..1], &ids)
        .unwrap();
    assert_eq!(
        ind.ind_xcorr_coelution_score,
        c.calc_separate_xcorr_contrast_coelution_score()
    );
    assert_eq!(ind.ind_log_sn_score, Vec::<f64>::new());
}

/// Derived: the library scores and the normalised RT score
/// (OpenSwathScoring.cpp `calculateLibraryScores`).
#[test]
fn library_scores() {
    let f = feature();
    let transitions: Vec<LightTransition> = [("t1", 2.0), ("t2", 3.0)]
        .iter()
        .map(|(n, lib)| LightTransition {
            transition_name: (*n).into(),
            library_intensity: *lib,
            ..LightTransition::default()
        })
        .collect();
    let compound = LightCompound {
        rt: 30.0,
        ..LightCompound::default()
    };
    let mut sc = OpenSwathScoring::new();
    sc.rt_normalization_factor = 100.0;
    let mut scores = OpenSwathScores::default();
    sc.calculate_library_scores(&f, &transitions, &compound, 20.0, &mut scores)
        .unwrap();
    let lib = MRMScoring::calc_library_score(&f, &transitions).unwrap();
    assert_eq!(scores.library_corr, lib.correlation);
    assert_eq!(scores.library_rootmeansquare, lib.rmsd);
    real_similar(scores.raw_rt_score, 10.0);
    real_similar(scores.norm_rt_score, 0.1);
    real_similar(scores.normalized_experimental_rt, 20.0);
}

/// Derived: the DIA scores on one SWATH map equal the DiaScoring values they
/// are built from; the per-transition mass errors are returned.
#[test]
fn dia_scores_compose_dia_scoring() {
    let map = SwathMap {
        sptr: Some(access(vec![spectrum(
            &[
                (499.99, 100.0),
                (500.0, 200.0),
                (600.0, 50.0),
                (601.0, 20.0),
            ],
            20.0,
        )])),
        ..SwathMap::default()
    };
    let transitions: Vec<LightTransition> = [("t1", 500.0), ("t2", 600.0)]
        .iter()
        .map(|(n, mz)| LightTransition {
            transition_name: (*n).into(),
            product_mz: *mz,
            precursor_mz: 700.0,
            library_intensity: 1.0,
            fragment_charge: 1,
            ..LightTransition::default()
        })
        .collect();
    let f = feature();
    let dia = DiaScoring::new();
    let mut sc = scorer("simple");
    sc.su.use_im_scores = false;
    let mut scores = OpenSwathScores::default();
    let ppm = sc
        .calculate_dia_scores(
            &f,
            &transitions,
            &[0.5, 0.5],
            std::slice::from_ref(&map),
            None,
            &dia,
            &LightCompound::default(),
            &mut scores,
            -1.0,
            &RangeBase::new(),
            None,
            -1,
        )
        .unwrap();
    let spectra = sc
        .fetch_spectrum_swath_into(
            map.sptr.as_ref().unwrap().as_ref(),
            20.0,
            1,
            &RangeBase::new(),
        )
        .unwrap();
    let md = dia
        .dia_massdiff_score(&transitions, &spectra, &[0.5, 0.5], &RangeBase::new())
        .unwrap();
    assert_eq!(ppm, md.diff_ppm);
    assert_eq!(scores.massdev_score, md.ppm_score);
    let iso = dia
        .dia_isotope_scores(&transitions, &spectra, &f, &RangeBase::new())
        .unwrap();
    assert_eq!(scores.isotope_correlation, iso.isotope_corr);
    // an empty compound sequence skips the b/y scores
    assert_eq!(scores.bseries_score, 0.0);
}
