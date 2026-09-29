// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! DIAHelper, DIAPrescoring and DIAScoring.
//!
//! Expected values are transcribed from the upstream class tests at the core
//! pin `bc9cc12` (evidence tier 3), cited per test as `<file>:<section>`:
//! `DIAHelper_test.cpp`, `DIAPrescoring_test.cpp` and `DIAScoring_test.cpp`
//! under `src/tests/class_tests/openms/source/`. Tests marked "derived" check an
//! independently derived invariant or a native guard (tier 4).

use openms::Error;
use openms::analysis::openswath::data_structures::SpectrumPtr;
use openms::analysis::openswath::data_structures::{BinaryDataArray, Spectrum};
use openms::analysis::openswath::dia_helper::{
    self, IntegratedWindows, adjust_extraction_window, integrate_window, integrate_window_sequence,
    integrate_windows, integrate_windows_sequence,
};
use openms::analysis::openswath::dia_prescoring::DiaPrescore;
use openms::analysis::openswath::dia_scoring::DiaScoring;
use openms::analysis::openswath::mock_objects::{MockFeature, MockMRMFeature};
use openms::analysis::openswath::transition_experiment::LightTransition;
use openms::chemistry::{AASequence, TheoreticalSpectrumGenerator};
use openms::kernel::ranges::RangeBase;
use openms::param::ParamValue;
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

fn range(center: f64, width: f64) -> RangeBase {
    let mut r = RangeBase::singular(center).unwrap();
    r.min_span_if_singular(width).unwrap();
    r
}

fn empty() -> RangeBase {
    RangeBase::new()
}

fn with_im(mz: Vec<f64>, intensity: Vec<f64>, im: Vec<f64>) -> Spectrum {
    let mut s = Spectrum::from_arrays(mz, intensity);
    s.push_data_array(BinaryDataArray::with_description(im, "Ion Mobility"));
    s
}

// ---------------------------------------------------------------------------
// DIAHelper
// ---------------------------------------------------------------------------

fn helper_im_spec() -> Spectrum {
    with_im(
        vec![100., 101., 102., 103., 104., 105., 106.],
        vec![1.; 7],
        vec![1., 2., 3., 4., 5., 6., 7.],
    )
}

fn helper_spec() -> Spectrum {
    Spectrum::from_arrays(vec![100., 101., 102., 103., 104., 105., 106.], vec![1.; 7])
}

/// DIAHelper_test.cpp:integrateWindow(SpectrumPtr, ...).
#[test]
fn helper_integrate_window_single() {
    let im = range(3.5, 3.0);
    let mz = range(102.0, 2.0);
    let mz2 = range(105.0, 8.0);
    // empty spectrum: no signal
    assert_eq!(
        integrate_window(&Spectrum::new(), &mz, &empty(), false).unwrap(),
        None
    );
    // IM requested but absent
    assert!(matches!(
        integrate_window(&helper_spec(), &mz, &im, false),
        Err(Error::MissingInformation(_))
    ));
    let r = integrate_window(&helper_im_spec(), &mz, &empty(), false)
        .unwrap()
        .unwrap();
    real_similar(r.mz, 101.5);
    real_similar(r.intensity, 2.0);
    real_similar(r.im, -1.0);
    let r = integrate_window(&helper_im_spec(), &mz2, &im, false)
        .unwrap()
        .unwrap();
    real_similar(r.im, 3.5);
    real_similar(r.intensity, 4.0);
    let r = integrate_window(&helper_im_spec(), &mz, &im, false)
        .unwrap()
        .unwrap();
    real_similar(r.im, 2.5);
    real_similar(r.intensity, 2.0);
}

/// DIAHelper_test.cpp:integrateWindow(SpectrumSequence, ...).
#[test]
fn helper_integrate_window_sequence() {
    let im = range(3.5, 3.0);
    let mz = range(102.0, 2.0);
    let mz2 = range(105.0, 8.0);
    assert_eq!(
        integrate_window_sequence(&[], &mz, &empty(), false).unwrap(),
        None
    );
    let empty_seq: Vec<SpectrumPtr> = vec![Arc::new(Spectrum::new())];
    assert_eq!(
        integrate_window_sequence(&empty_seq, &mz, &empty(), false).unwrap(),
        None
    );
    let seq: Vec<SpectrumPtr> = vec![Arc::new(helper_im_spec())];
    let r = integrate_window_sequence(&seq, &mz, &empty(), false)
        .unwrap()
        .unwrap();
    real_similar(r.mz, 101.5);
    real_similar(r.intensity, 2.0);
    real_similar(r.im, -1.0);
    let r = integrate_window_sequence(&seq, &mz2, &im, false)
        .unwrap()
        .unwrap();
    real_similar(r.im, 3.5);
    real_similar(r.intensity, 4.0);
    let r = integrate_window_sequence(&seq, &mz, &im, false)
        .unwrap()
        .unwrap();
    real_similar(r.im, 2.5);
    real_similar(r.intensity, 2.0);
}

fn check_zero_windows(w: &IntegratedWindows) {
    assert_eq!(w.len(), 3);
    assert_eq!(w.mz.len(), 3);
    assert_eq!(w.im.len(), 3);
    for (i, center) in [101.0, 103.0, 105.0].iter().enumerate() {
        real_similar(w.intensity[i], 0.0);
        real_similar(w.mz[i], *center);
        real_similar(w.im[i], -1.0);
    }
}

/// DIAHelper_test.cpp:integrateWindows(SpectrumPtr, ...).
#[test]
fn helper_integrate_windows_single() {
    let im = range(3.5, 3.0);
    let windows = [101.0, 103.0, 105.0];
    let w = integrate_windows(&Spectrum::new(), &windows, 2.0, &empty(), true).unwrap();
    assert!(w.is_empty() && w.mz.is_empty() && w.im.is_empty());
    let w = integrate_windows(&Spectrum::new(), &windows, 2.0, &empty(), false).unwrap();
    check_zero_windows(&w);
    let w = integrate_windows(&helper_spec(), &windows, 2.0, &empty(), false).unwrap();
    assert_eq!(w.len(), 3);
    assert_eq!(w.im.len(), 3);
    real_similar(w.intensity[0], 2.0);
    real_similar(w.mz[0], 100.5);
    real_similar(w.im[0], -1.0);
    let w = integrate_windows(&helper_im_spec(), &[102.0], 2.0, &im, false).unwrap();
    assert_eq!(w.len(), 1);
    real_similar(w.intensity[0], 2.0);
    real_similar(w.mz[0], 101.5);
    real_similar(w.im[0], 2.5);
}

/// DIAHelper_test.cpp:integrateWindows(SpectrumSequence, ...).
#[test]
fn helper_integrate_windows_sequence() {
    let im = range(3.5, 3.0);
    let windows = [101.0, 103.0, 105.0];
    let empty_seq: Vec<SpectrumPtr> = vec![Arc::new(Spectrum::new())];
    assert!(matches!(
        integrate_windows_sequence(&empty_seq, &[], 2.0, &empty(), false),
        Err(Error::MissingInformation(_))
    ));
    let w = integrate_windows_sequence(&empty_seq, &windows, 2.0, &im, true).unwrap();
    assert!(w.is_empty() && w.mz.is_empty() && w.im.is_empty());
    let w = integrate_windows_sequence(&empty_seq, &windows, 2.0, &empty(), false).unwrap();
    check_zero_windows(&w);
    let seq: Vec<SpectrumPtr> = vec![Arc::new(helper_spec())];
    let w = integrate_windows_sequence(&seq, &windows, 2.0, &empty(), false).unwrap();
    assert_eq!(w.len(), 3);
    real_similar(w.intensity[0], 2.0);
    real_similar(w.mz[0], 100.5);
    real_similar(w.im[0], -1.0);
    let seq: Vec<SpectrumPtr> = vec![Arc::new(helper_im_spec())];
    let w = integrate_windows_sequence(&seq, &[102.0], 2.0, &im, false).unwrap();
    assert_eq!(w.len(), 1);
    real_similar(w.intensity[0], 2.0);
    real_similar(w.mz[0], 101.5);
    real_similar(w.im[0], 2.5);
}

/// Derived: an empty window with an IM range reports the IM range center, as
/// `DIAHelper.cpp` pushes `range_im.center()`.
#[test]
fn helper_empty_window_reports_im_center() {
    let im = range(3.5, 3.0);
    let w = integrate_windows(&helper_im_spec(), &[200.0], 2.0, &im, false).unwrap();
    real_similar(w.im[0], 3.5);
    real_similar(w.mz[0], 200.0);
}

/// DIAHelper_test.cpp:[EXTRA] adjustExtractionWindow.
#[test]
fn helper_adjust_extraction_window() {
    let (l, r) = adjust_extraction_window(500.0, 500.0, 0.5, false).unwrap();
    real_similar(l, 500.0 - 0.25);
    real_similar(r, 500.0 + 0.25);
    let (l, r) = adjust_extraction_window(500.0, 500.0, 10.0, true).unwrap();
    real_similar(l, 500.0 - 500.0 * 5.0 / 1e6);
    real_similar(r, 500.0 + 500.0 * 5.0 / 1e6);
    let (l, r) = adjust_extraction_window(500.0, 500.0, 1000.0, false).unwrap();
    real_similar(l, 0.0);
    real_similar(r, 500.0 + 500.0);
    // derived: the source's positivity precondition is checked
    assert!(adjust_extraction_window(500.0, 500.0, 0.0, false).is_err());
}

fn by_generator() -> TheoreticalSpectrumGenerator {
    TheoreticalSpectrumGenerator {
        add_metainfo: true,
        ..TheoreticalSpectrumGenerator::default()
    }
}

/// DIAScoring_test.cpp:[EXTRA] getBYSeries (the DIAHelper_test section only
/// calls the functions without checks).
#[test]
fn helper_by_series() {
    let generator = by_generator();
    let mut a = AASequence::parse("SYVAWDR").unwrap();
    let (b, y) = dia_helper::by_series(&a, &generator, 1).unwrap();
    let b_exp = [251.10323, 350.17164, 421.20875, 607.28807, 722.31501];
    let y_exp = [
        175.11955, 290.14649, 476.22580, 547.26291, 646.33133, 809.39466,
    ];
    assert_eq!(b.len(), 5);
    assert_eq!(y.len(), 6);
    for (got, exp) in b.iter().zip(b_exp) {
        real_similar(*got, exp);
    }
    for (got, exp) in y.iter().zip(y_exp) {
        real_similar(*got, exp);
    }
    a.set_modification(1, "Phospho").unwrap();
    let (b, y) = dia_helper::by_series(&a, &generator, 1).unwrap();
    assert_eq!(b.len(), 5);
    assert_eq!(y.len(), 6);
    for (got, exp) in b.iter().zip(b_exp) {
        real_similar(*got, exp + 79.9657);
    }
    for (i, (got, exp)) in y.iter().zip(y_exp).enumerate() {
        real_similar(*got, if i == 5 { exp + 79.9657 } else { exp });
    }
    // DIAHelper_test.cpp:[EXTRA] getBYSeries_test calls getTheorMasses too.
    let masses = dia_helper::theoretical_masses(&a, &generator, 1).unwrap();
    assert_eq!(masses.len(), 11);
}

/// Derived: a generator without ion names is refused rather than read blind.
#[test]
fn helper_by_series_needs_ion_names() {
    let a = AASequence::parse("SYVAWDR").unwrap();
    assert!(matches!(
        dia_helper::by_series(&a, &TheoreticalSpectrumGenerator::default(), 1),
        Err(Error::MissingInformation(_))
    ));
    assert!(dia_helper::by_series(&a, &by_generator(), 0).is_err());
}

/// Derived: averagine helpers space isotopes by `mannmass / |charge|` and
/// reject charge 0, as `DIAHelper.cpp` throws `IllegalArgument`.
#[test]
fn helper_averagine_spacing_and_charge() {
    let d = dia_helper::averagine_isotope_distribution(500.0, -2, 4, 1.00048).unwrap();
    assert_eq!(d.len(), 4);
    for (i, (mass, _)) in d.iter().enumerate() {
        real_similar(*mass, 500.0 + i as f64 * 1.00048 / 2.0);
    }
    assert!(dia_helper::averagine_isotope_distribution(500.0, 0, 4, 1.00048).is_err());
    let mut out = Vec::new();
    assert!(dia_helper::add_single_peak_isotopes_to_spec(500.0, 1.0, &mut out, 4, 0).is_err());
    assert!(out.is_empty());
    let mut spec = vec![(10.0, 1.0)];
    dia_helper::add_preisotope_weights(&[12.0], &mut spec, 2, -0.5, 1.0, 1).unwrap();
    assert_eq!(spec, vec![(10.0, 1.0), (10.0, -0.5), (11.0, -0.5)]);
}

// ---------------------------------------------------------------------------
// DIAPrescoring
// ---------------------------------------------------------------------------

fn tr(product_mz: f64, name: &str, library_intensity: f64) -> LightTransition {
    LightTransition {
        product_mz,
        fragment_charge: 1,
        transition_name: name.into(),
        library_intensity,
        ..LightTransition::default()
    }
}

fn peak_block(first: &[f64]) -> Vec<f64> {
    first.to_vec()
}

const PRESCORE_INT: [f64; 48] = [
    10., 20., 50., 100., 50., 20., 10., 3., 7., 15., 30., 15., 7., 3., 1., 3., 9., 15., 9., 3., 1.,
    3., 9., 3., 10., 20., 50., 100., 50., 20., 10., 3., 7., 15., 30., 15., 7., 3., 1., 3., 9., 15.,
    9., 3., 1., 3., 9., 3.,
];

fn prescore_transitions() -> Vec<LightTransition> {
    vec![tr(500.0, "group1", 5.0), tr(600.0, "group2", 5.0)]
}

/// DIAPrescoring_test.cpp:test score function with perfect first transition.
#[test]
fn prescore_perfect_first_transition() {
    let mz = vec![
        499.97, 499.98, 499.99, 500.0, 500.01, 500.02, 500.03, 500.97, 500.98, 500.99, 501.0,
        501.01, 501.02, 501.03, 501.97, 501.98, 501.99, 502.0, 502.01, 502.02, 502.03, 502.99,
        503.0, 503.01, 599.97, 599.98, 599.99, 600.0, 600.01, 600.02, 600.03, 600.97, 600.98,
        600.99, 601.0, 601.01, 601.02, 601.03, 601.97, 601.98, 601.99, 602.0, 602.01, 602.02,
        602.03, 602.99, 603.0, 603.01,
    ];
    let seq: Vec<SpectrumPtr> = vec![Arc::new(Spectrum::from_arrays(
        mz,
        peak_block(&PRESCORE_INT),
    ))];
    let r = DiaPrescore::with_window(0.05, 4, 4)
        .score(&seq, &prescore_transitions(), &empty())
        .unwrap();
    real_similar(r.dotprod, 0.916131286812994);
    real_similar(r.manhattan, 0.23670593984202);
}

fn second_half_mz() -> Vec<f64> {
    vec![
        599.97, 599.98, 599.99, 600.0, 600.01, 600.02, 600.03, 600.97, 600.98, 600.99, 601.0,
        601.01, 601.02, 601.03, 601.97, 601.98, 601.99, 602.0, 602.01, 602.02, 602.03, 602.99,
        603.0, 603.01,
    ]
}

fn shifted_mz() -> Vec<f64> {
    let mut mz = vec![
        498.97, 498.98, 498.99, 499.0, 499.01, 499.02, 499.03, 499.97, 499.98, 499.99, 500.0,
        500.01, 500.02, 500.03, 500.97, 500.98, 500.99, 501.0, 501.01, 501.02, 501.03, 501.99,
        502.0, 502.01,
    ];
    mz.extend(second_half_mz());
    mz
}

/// DIAPrescoring_test.cpp:test score function missing first transition.
#[test]
fn prescore_missing_first_transition() {
    let seq: Vec<SpectrumPtr> = vec![Arc::new(Spectrum::from_arrays(
        second_half_mz(),
        PRESCORE_INT[24..].to_vec(),
    ))];
    let r = DiaPrescore::with_window(0.05, 4, 4)
        .score(&seq, &prescore_transitions(), &empty())
        .unwrap();
    real_similar(r.dotprod, 0.627263258948172);
    real_similar(r.manhattan, 0.984211129641047);
}

/// DIAPrescoring_test.cpp:test score function with shifted first transition.
#[test]
fn prescore_shifted_first_transition() {
    let seq: Vec<SpectrumPtr> = vec![Arc::new(Spectrum::from_arrays(
        shifted_mz(),
        PRESCORE_INT.to_vec(),
    ))];
    let r = DiaPrescore::with_window(0.05, 4, 4)
        .score(&seq, &prescore_transitions(), &empty())
        .unwrap();
    real_similar(r.dotprod, 0.43738312515644);
    real_similar(r.manhattan, 0.557433222328531);
}

/// DIAPrescoring_test.cpp:test score function missing first transition due to
/// different ion mobility.
#[test]
fn prescore_missing_first_transition_by_im() {
    let im = vec![
        1., 1., 3., 1., 1., 1., 1., 2., 2., 3., 1., 2., 1., 2., 1., 2., 1., 2., 3., 2., 2., 2., 2.,
        2., 7., 6., 7., 6., 8., 6., 8., 7., 7., 7., 8., 7., 7., 8., 8., 6., 8., 8., 6., 6., 6., 6.,
        8., 6.,
    ];
    let seq: Vec<SpectrumPtr> = vec![Arc::new(with_im(shifted_mz(), PRESCORE_INT.to_vec(), im))];
    let r = DiaPrescore::with_window(0.05, 4, 4)
        .score(&seq, &prescore_transitions(), &range(7.0, 2.0))
        .unwrap();
    real_similar(r.dotprod, 0.627263258948172);
    real_similar(r.manhattan, 0.984211129641047);
}

/// DIAPrescoring_test.cpp:DiaPrescore() — and derived: the defaults round-trip
/// through a `Param` (DIAPrescoring.cpp `defineDefaults`).
#[test]
fn prescore_defaults() {
    let p = DiaPrescore::defaults().unwrap();
    assert_eq!(
        p.value("dia_extraction_window").unwrap(),
        &ParamValue::Float(0.1)
    );
    assert_eq!(p.value("nr_isotopes").unwrap(), &ParamValue::Integer(4));
    assert_eq!(p.value("nr_charges").unwrap(), &ParamValue::Integer(4));
    assert_eq!(DiaPrescore::from_param(&p).unwrap(), DiaPrescore::new());
}

// ---------------------------------------------------------------------------
// DIAScoring
// ---------------------------------------------------------------------------

const SCORE_MZ: [f64; 48] = [
    498.97, 498.98, 498.99, 499.0, 499.01, 499.02, 499.03, 499.97, 499.98, 499.99, 500.0, 500.01,
    500.02, 500.03, 500.97, 500.98, 500.99, 501.0, 501.01, 501.02, 501.03, 501.99, 502.0, 502.01,
    599.97, 599.98, 599.99, 600.0, 600.01, 600.02, 600.03, 600.97, 600.98, 600.99, 601.0, 601.01,
    601.02, 601.03, 601.97, 601.98, 601.99, 602.0, 602.01, 602.02, 602.03, 602.99, 603.0, 603.01,
];

fn prepare_spectrum() -> Spectrum {
    Spectrum::from_arrays(SCORE_MZ.to_vec(), PRESCORE_INT.to_vec())
}

fn prepare_im_spectrum() -> Spectrum {
    let intensity = vec![
        10., 20., 50., 100., 50., 20., 10., 3., 7., 15., 30., 200000., 15., 7., 3., 1., 3., 9.,
        15., 9., 3., 1., 3., 9., 3., 10., 20., 50., 100., 50., 20., 10., 3., 7., 15., 30., 15., 7.,
        3., 1., 3., 9., 15., 9., 3., 1., 3., 9., 200000., 3.,
    ];
    let mz = vec![
        498.97, 498.98, 498.99, 499.0, 499.01, 499.02, 499.03, 499.97, 499.98, 499.99, 500.0,
        500.005, 500.01, 500.02, 500.03, 500.97, 500.98, 500.99, 501.0, 501.01, 501.02, 501.03,
        501.99, 502.0, 502.01, 599.97, 599.98, 599.99, 600.0, 600.01, 600.02, 600.03, 600.97,
        600.98, 600.99, 601.0, 601.01, 601.02, 601.03, 601.97, 601.98, 601.99, 602.0, 602.01,
        602.02, 602.03, 602.99, 603.0, 603.01, 603.01,
    ];
    let im = vec![
        2., 3., 2., 3., 4., 2., 3., 2., 3., 2., 3., 6., 2., 3., 3., 2., 3., 2., 3., 4., 2., 3., 2.,
        3., 3., 2., 3., 2., 3., 4., 2., 3., 2., 3., 2., 3., 4., 2., 3., 2., 3., 2., 3., 4., 2., 3.,
        2., 3., 6., 3.,
    ];
    assert_eq!(im.len(), mz.len());
    assert_eq!(im.len(), intensity.len());
    with_im(mz, intensity, im)
}

fn shift(mut s: Spectrum) -> Spectrum {
    let data = &mut s.mz_array_mut().data;
    let n = data.len();
    let half = n as f64 / 2.0;
    for (i, v) in data.iter_mut().enumerate() {
        // DIAScoring_test.cpp prepareShiftedSpectrum: the loop bound is `size / 2.0`
        let ppm = if (i as f64) < half { 15.0 } else { 10.0 };
        *v += *v / 1_000_000.0 * ppm;
    }
    s
}

fn seq(s: Spectrum) -> Vec<SpectrumPtr> {
    vec![Arc::new(s)]
}

fn p_dia(large: bool) -> DiaScoring {
    let mut p = DiaScoring::defaults().unwrap();
    let set = |p: &mut openms::param::Param, key: &str, value: ParamValue| {
        p.set_value(key, value, "", &[]).unwrap();
    };
    set(
        &mut p,
        "dia_extraction_window",
        ParamValue::Float(if large { 0.5 } else { 0.05 }),
    );
    set(
        &mut p,
        "dia_extraction_unit",
        ParamValue::String("Th".into()),
    );
    set(&mut p, "dia_centroided", ParamValue::String("false".into()));
    set(
        &mut p,
        "dia_byseries_intensity_min",
        ParamValue::Float(30.0),
    );
    set(&mut p, "dia_byseries_ppm_diff", ParamValue::Float(50.0));
    set(&mut p, "dia_nr_isotopes", ParamValue::Integer(4));
    set(&mut p, "dia_nr_charges", ParamValue::Integer(4));
    DiaScoring::from_param(&p).unwrap()
}

fn mrm_feature(total: f32) -> MockMRMFeature {
    let mut features = BTreeMap::new();
    features.insert(
        "group1".to_string(),
        Arc::new(MockFeature {
            intensity: 0.3,
            ..MockFeature::default()
        }),
    );
    features.insert(
        "group2".to_string(),
        Arc::new(MockFeature {
            intensity: 0.7,
            ..MockFeature::default()
        }),
    );
    MockMRMFeature {
        features,
        intensity: total,
        ..MockMRMFeature::default()
    }
}

fn plain_tr(product_mz: f64, charge: i8, name: &str) -> LightTransition {
    LightTransition {
        product_mz,
        fragment_charge: charge,
        transition_name: name.into(),
        ..LightTransition::default()
    }
}

/// DIAScoring_test.cpp:[EXTRA] forward / forward with IM / forward negative
/// charge (-1 and -2) dia_isotope_scores.
#[test]
fn scoring_isotope_forward() {
    let s = p_dia(false);
    let f = mrm_feature(0.7);
    for (spectra, im) in [
        (seq(prepare_spectrum()), empty()),
        (seq(prepare_im_spectrum()), range(3.0, 2.0)),
    ] {
        let r = s
            .dia_isotope_scores(&[plain_tr(600.0, 1, "group2")], &spectra, &f, &im)
            .unwrap();
        real_similar(r.isotope_corr, 0.995335798317618);
        real_similar(r.isotope_overlap, 0.0);
    }
    let spectra = seq(prepare_spectrum());
    let r = s
        .dia_isotope_scores(&[plain_tr(600.0, -1, "group2")], &spectra, &f, &empty())
        .unwrap();
    real_similar(r.isotope_corr, 0.995335798317618);
    real_similar(r.isotope_overlap, 0.0);
    let r = s
        .dia_isotope_scores(&[plain_tr(600.0, -2, "group2")], &spectra, &f, &empty())
        .unwrap();
    real_similar(r.isotope_corr, 0.717092518007138);
    real_similar(r.isotope_overlap, 0.0);
}

/// DIAScoring_test.cpp:[EXTRA] backward / backward with IM dia_isotope_scores.
#[test]
fn scoring_isotope_backward() {
    let s = p_dia(false);
    let f = mrm_feature(0.3);
    for (spectra, im) in [
        (seq(prepare_spectrum()), empty()),
        (seq(prepare_im_spectrum()), range(3.0, 2.0)),
    ] {
        let r = s
            .dia_isotope_scores(&[plain_tr(500.0, 1, "group1")], &spectra, &f, &im)
            .unwrap();
        real_similar(r.isotope_corr, 0.959692139694113);
        real_similar(r.isotope_overlap, 1.0);
    }
}

/// DIAScoring_test.cpp:dia_isotope_scores and [EXTRA] with IM.
#[test]
fn scoring_isotope_both() {
    let s = p_dia(false);
    let f = mrm_feature(1.0);
    let transitions = [plain_tr(500.0, 1, "group1"), plain_tr(600.0, 1, "group2")];
    for (spectra, im) in [
        (seq(prepare_spectrum()), empty()),
        (seq(prepare_im_spectrum()), range(3.0, 2.0)),
    ] {
        let r = s
            .dia_isotope_scores(&transitions, &spectra, &f, &im)
            .unwrap();
        real_similar(
            r.isotope_corr,
            0.995335798317618 * 0.7 + 0.959692139694113 * 0.3,
        );
        real_similar(r.isotope_overlap, 0.0 * 0.7 + 1.0 * 0.3);
    }
}

/// DIAScoring_test.cpp:dia_ms1_isotope_scores and [EXTRA] with IM (the
/// averagine form).
#[test]
fn scoring_ms1_isotope_averagine() {
    let s = p_dia(false);
    for (spectra, im) in [
        (seq(prepare_spectrum()), empty()),
        (seq(prepare_im_spectrum()), range(3.0, 2.0)),
    ] {
        let r = s
            .dia_ms1_isotope_scores_averagine(500.0, &spectra, 1, &im)
            .unwrap();
        real_similar(r.isotope_corr, 0.959692139694113);
        real_similar(r.isotope_overlap, 240.0 / 74.0);
        let r = s
            .dia_ms1_isotope_scores_averagine(500.0, &spectra, 2, &im)
            .unwrap();
        real_similar(r.isotope_corr, 0.680028918385546);
        real_similar(r.isotope_overlap, 240.0 / 74.0);
        let r = s
            .dia_ms1_isotope_scores_averagine(499.0, &spectra, 1, &im)
            .unwrap();
        real_similar(r.isotope_corr, 0.995485552148335);
        real_similar(r.isotope_overlap, 0.0);
    }
}

/// DIAScoring_test.cpp:dia_massdiff_score and [EXTRA] with IM.
#[test]
fn scoring_massdiff() {
    let s = p_dia(true);
    let transitions = [plain_tr(500.0, 1, "group1"), plain_tr(600.0, 1, "group2")];
    for (spectra, im) in [
        (seq(shift(prepare_spectrum())), empty()),
        (seq(shift(prepare_im_spectrum())), range(3.0, 2.0)),
    ] {
        let r = s
            .dia_massdiff_score(&transitions, &spectra, &[0.7, 0.3], &im)
            .unwrap();
        real_similar(r.ppm_score, (15.0 + 10.0) / 2.0);
        real_similar(r.ppm_score_weighted, 15.0 * 0.7 + 10.0 * 0.3);
        assert_eq!(r.diff_ppm.len(), 2);
    }
    // derived: the library intensities must cover every transition
    let spectra = seq(prepare_spectrum());
    assert!(
        s.dia_massdiff_score(&transitions, &spectra, &[0.7], &empty())
            .is_err()
    );
}

/// DIAScoring_test.cpp:dia_ms1_massdiff_score and [EXTRA] with IM.
#[test]
fn scoring_ms1_massdiff() {
    let s = p_dia(true);
    for (spectra, im) in [
        (seq(shift(prepare_spectrum())), empty()),
        (seq(shift(prepare_im_spectrum())), range(3.0, 2.0)),
    ] {
        let r = s.dia_ms1_massdiff_score(500.0, &spectra, &im).unwrap();
        assert!(r.signal_found);
        real_similar(r.ppm_score, 15.0);
        let r = s.dia_ms1_massdiff_score(600.0, &spectra, &im).unwrap();
        assert!(r.signal_found);
        real_similar(r.ppm_score, 10.0);
        let r = s.dia_ms1_massdiff_score(100.0, &spectra, &im).unwrap();
        assert!(!r.signal_found);
        real_similar(r.ppm_score, 0.5 * 1_000_000.0 / 100.0);
    }
}

/// DIAScoring_test.cpp:dia_by_ion_score.
#[test]
fn scoring_by_ion() {
    let s = p_dia(false);
    let mz = vec![
        350.17164,
        421.20875,
        421.20875 + 79.9657,
        547.26291,
        646.33133,
        809.39466 + 79.9657,
    ];
    let spectra = seq(Spectrum::from_arrays(mz, vec![100.0; 6]));
    let mut a = AASequence::parse("SYVAWDR").unwrap();
    let r = s.dia_by_ion_score(&spectra, &a, 1, &empty()).unwrap();
    real_similar(r.bseries_score, 2.0);
    real_similar(r.yseries_score, 2.0);
    a.set_modification(1, "Phospho").unwrap();
    let r = s.dia_by_ion_score(&spectra, &a, 1, &empty()).unwrap();
    real_similar(r.bseries_score, 1.0);
    real_similar(r.yseries_score, 3.0);
}

/// DIAScoring_test.cpp:score_with_isotopes and [EXTRA] with IM.
#[test]
fn scoring_with_isotopes() {
    let s = p_dia(false);
    let transitions = prescore_transitions();
    for (spectra, im) in [
        (seq(prepare_spectrum()), empty()),
        (seq(prepare_im_spectrum()), range(3.0, 2.0)),
    ] {
        let r = s.score_with_isotopes(&spectra, &transitions, &im).unwrap();
        real_similar(r.dotprod, 0.43738312458795);
        real_similar(r.manhattan, 0.55743322213368);
    }
}

/// Derived: DIAScoring.cpp constructor defaults, and the parameter round trip.
#[test]
fn scoring_defaults() {
    let p = DiaScoring::defaults().unwrap();
    assert_eq!(
        p.value("dia_extraction_window").unwrap(),
        &ParamValue::Float(0.05)
    );
    assert_eq!(
        p.value("dia_extraction_unit").unwrap(),
        &ParamValue::String("Th".into())
    );
    assert_eq!(
        p.value("dia_byseries_intensity_min").unwrap(),
        &ParamValue::Float(300.0)
    );
    assert_eq!(
        p.value("dia_byseries_ppm_diff").unwrap(),
        &ParamValue::Float(10.0)
    );
    assert_eq!(p.value("dia_nr_isotopes").unwrap(), &ParamValue::Integer(4));
    assert_eq!(p.value("dia_nr_charges").unwrap(), &ParamValue::Integer(4));
    assert_eq!(
        p.value("peak_before_mono_max_ppm_diff").unwrap(),
        &ParamValue::Float(20.0)
    );
    let s = DiaScoring::from_param(&p).unwrap();
    assert_eq!(s, DiaScoring::new());
    assert!(!s.is_dia_extraction_ppm());
    real_similar(s.dia_extraction_window(), 0.05);
}

/// Derived: a peak-group intensity of zero yields zero isotope scores, as the
/// source returns early.
#[test]
fn scoring_isotope_zero_feature_intensity() {
    let s = p_dia(false);
    let spectra = seq(prepare_spectrum());
    let r = s
        .dia_isotope_scores(
            &[plain_tr(600.0, 1, "group2")],
            &spectra,
            &mrm_feature(0.0),
            &empty(),
        )
        .unwrap();
    assert_eq!(r.isotope_corr, 0.0);
    assert_eq!(r.isotope_overlap, 0.0);
}
