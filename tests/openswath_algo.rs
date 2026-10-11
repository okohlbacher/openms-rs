// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! OPENSWATHALGO: data structures, SwathMap, transitions, mocks, StatsHelpers,
//! Scoring and DataFrameWriter.
//!
//! Expected values are transcribed from the upstream class tests at the core
//! pin `bc9cc12` (evidence tier 3), cited per test as
//! `<file>:<section>`: `openswathalgo/Datastructures_test.cpp`,
//! `openswathalgo/SwathMap_test.cpp`, `openswathalgo/Scoring_test.cpp` and
//! `openms/source/StatsHelpers_test.cpp`. Tests marked "derived" check an
//! independently derived invariant or a native guard (tier 4).

use openms::Error;
use openms::analysis::openswath::data_frame_writer::{CsvWriter, DataFrameWriter, DataMatrix};
use openms::analysis::openswath::data_structures::{
    BinaryDataArray, Chromatogram, Spectrum, SpectrumMeta, is_ion_mobility_description,
};
use openms::analysis::openswath::itransition::{
    MRMFeatureAccess, SignalToNoiseAccess, TransitionGroupAccess,
};
use openms::analysis::openswath::mock_objects::{
    MockFeature, MockMRMFeature, MockSignalToNoise, MockTransitionGroup,
};
use openms::analysis::openswath::scoring::{self, XCorrArrayType};
use openms::analysis::openswath::stats_helpers::{self, MeanAndStddev};
use openms::analysis::openswath::swath_map::SwathMap;
use openms::analysis::openswath::transition_experiment::{
    FragmentIonType, LightCompound, LightModification, LightTargetedExperiment, LightTransition,
    TransitionFlags,
};
use openms::analysis::openswath::transition_helper;
use openms::system::file::TempDir;
use std::collections::{BTreeMap, HashSet};
use std::sync::Arc;

fn real_similar(a: f64, b: f64) {
    let absdiff = (a - b).abs();
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{a} !~ {b} (absdiff {absdiff}, ratio {ratio})"
    );
}

// ---------------------------------------------------------------------------
// DataStructures
// ---------------------------------------------------------------------------

/// Datastructures_test.cpp:OSSpectrum_empty.
#[test]
fn spectrum_empty() {
    let s = Spectrum::new();
    assert!(s.drift_time_array().is_none());
    assert_eq!(s.mz_array().data.len(), 0);
    assert_eq!(s.intensity_array().data.len(), 0);
    assert_eq!(s.data_arrays().len(), 2);
}

/// Datastructures_test.cpp:OSSpectrum_data.
#[test]
fn spectrum_data_legacy_im_name() {
    let mut s = Spectrum::new();
    s.set_mz_array(BinaryDataArray::new(vec![1.5]));
    s.set_intensity_array(BinaryDataArray::new(vec![100.1]));
    s.push_data_array(BinaryDataArray::with_description(
        vec![300.1],
        "Ion Mobility",
    ));
    let im = s.drift_time_array().expect("ion-mobility array");
    assert_eq!(s.mz_array().data.len(), 1);
    assert_eq!(s.intensity_array().data.len(), 1);
    assert_eq!(im.data.len(), 1);
    real_similar(s.mz_array().data[0], 1.5);
    real_similar(s.intensity_array().data[0], 100.1);
    real_similar(im.data[0], 300.1);
}

/// Datastructures_test.cpp:OSSpectrum_data_2.
#[test]
fn spectrum_data_im_name_variants() {
    let mut s = Spectrum::new();
    s.set_mz_array(BinaryDataArray::new(vec![1.5]));
    s.set_intensity_array(BinaryDataArray::new(vec![100.1]));
    s.push_data_array(BinaryDataArray::with_description(
        vec![300.1],
        "Ion Mobility (MS:1002476)",
    ));
    real_similar(s.drift_time_array().unwrap().data[0], 300.1);
    let last = s.data_arrays().len() - 1;
    let cases = [
        ("", false),
        ("Ion Mobility (blah)", true),
        ("Ion mOBILITY (blah)", false),
        ("raw inverse reduced ion mobility array", true),
        ("mean ion mobility array", true),
        ("inverse reduced ion mobility", true),
        ("mean inverse reduced ion mobility array", true),
        ("Ion Mobility Centroid", true),
    ];
    for (description, found) in cases {
        s.data_arrays_mut()[last].description = description.to_string();
        assert_eq!(s.drift_time_array().is_some(), found, "{description:?}");
        assert_eq!(is_ion_mobility_description(description), found);
    }
}

/// Derived: `setDriftTimeArray` appends under the description "Ion Mobility".
#[test]
fn set_drift_time_array_appends_and_names() {
    let mut s = Spectrum::from_arrays(vec![1.0], vec![2.0]);
    s.set_drift_time_array(BinaryDataArray::with_description(vec![0.5], "anything"));
    assert_eq!(s.data_arrays().len(), 3);
    assert_eq!(s.drift_time_array().unwrap().description, "Ion Mobility");
    assert_eq!(s.drift_time_array_index(), Some(2));
}

/// Derived: the two leading arrays cannot be removed (the source's null arrays).
#[test]
fn set_data_arrays_refuses_fewer_than_two() {
    let mut s = Spectrum::from_arrays(vec![1.0], vec![2.0]);
    let err = s.set_data_arrays(vec![BinaryDataArray::default()]);
    assert!(matches!(err, Err(Error::InvalidValue(_))));
    assert_eq!(s.mz_array().data, vec![1.0]);
    let mut c = Chromatogram::new();
    assert!(c.set_data_arrays(Vec::new()).is_err());
    c.set_data_arrays(vec![
        BinaryDataArray::new(vec![1.0]),
        BinaryDataArray::new(vec![5.0]),
        BinaryDataArray::with_description(vec![7.0], "extra"),
    ])
    .unwrap();
    assert_eq!(c.time_array().data, vec![1.0]);
    assert_eq!(c.intensity_array().data, vec![5.0]);
    assert_eq!(c.data_arrays().len(), 3);
}

/// Derived: `SpectrumMeta::RTLess`.
#[test]
fn spectrum_meta_rt_less() {
    let a = SpectrumMeta {
        rt: 1.0,
        ..SpectrumMeta::default()
    };
    let b = SpectrumMeta {
        rt: 2.0,
        ..SpectrumMeta::default()
    };
    assert!(SpectrumMeta::rt_less(&a, &b));
    assert!(!SpectrumMeta::rt_less(&b, &a));
    assert!(!SpectrumMeta::rt_less(&a, &a));
}

// ---------------------------------------------------------------------------
// SwathMap
// ---------------------------------------------------------------------------

/// SwathMap_test.cpp:testIsEqual.
#[test]
fn swath_map_is_equal() {
    let map1 = SwathMap::default();
    let map2 = SwathMap::default();
    assert!(map1.is_equal(&map2));
    let map3 = SwathMap::new(1.0, 2.0, 1.5, false);
    let map4 = SwathMap::new(1.0, 2.0, 1.5, false);
    let map5 = SwathMap::new(1.0, 2.0, 1.5, true);
    assert!(map3.is_equal(&map4));
    assert!(!map3.is_equal(&map5));
    let map8 = SwathMap::with_im(1.0, 2.0, 1.5, -1.0, -1.0, false);
    assert!(map3.is_equal(&map8));
    let map9 = SwathMap::with_im(1.0, 2.0, 1.5, 1.0, 1.1, false);
    let map10 = SwathMap::with_im(1.0, 2.0, 1.5, 1.0, 1.1, false);
    assert!(map9.is_equal(&map10));
    let map11 = SwathMap::with_im(1.0, 2.0, 1.5, 1.3, 1.4, false);
    let map12 = SwathMap::with_im(1.0, 2.0, 1.5, 1.0, 1.2, false);
    assert!(!map9.is_equal(&map11));
    assert!(!map9.is_equal(&map12));
    // Derived: the six-map section's maps 6 and 7 differ by m/z bounds.
    let map6 = SwathMap::new(1.0, 3.0, 2.0, false);
    let map7 = SwathMap::new(2.0, 3.0, 2.5, false);
    assert!(!map3.is_equal(&map6));
    assert!(!map3.is_equal(&map7));
    // Derived: the tolerance is strict.
    let near = SwathMap::new(1.0 + 1e-3, 2.0, 1.5, false);
    assert!(!map3.is_equal_with_tolerance(&near, 1e-3 / 2.0));
    assert!(map3.is_equal_with_tolerance(&near, 1e-2));
    assert!(map1.sptr.is_none());
    assert_eq!((map1.im_lower, map1.im_upper), (-1.0, -1.0));
}

// ---------------------------------------------------------------------------
// TransitionExperiment, Transitions, TransitionHelper
// ---------------------------------------------------------------------------

/// Derived from the header's `stringToFragmentIonType`/`fragmentIonTypeToString`.
#[test]
fn fragment_ion_type_round_trip() {
    let all = [
        FragmentIonType::AIon,
        FragmentIonType::BIon,
        FragmentIonType::CIon,
        FragmentIonType::XIon,
        FragmentIonType::YIon,
        FragmentIonType::ZIon,
        FragmentIonType::ZPrimeIon,
        FragmentIonType::ZDotIon,
        FragmentIonType::Precursor,
        FragmentIonType::BMinusH2O,
        FragmentIonType::YMinusH2O,
        FragmentIonType::BMinusNH3,
        FragmentIonType::YMinusNH3,
        FragmentIonType::Empty,
    ];
    for t in all {
        assert_eq!(FragmentIonType::from_type_str(t.as_str()), t);
    }
    assert_eq!(
        FragmentIonType::from_type_str("precursor"),
        FragmentIonType::Precursor
    );
    assert_eq!(
        FragmentIonType::from_type_str("b-H20"),
        FragmentIonType::BMinusH2O
    );
    assert_eq!(
        FragmentIonType::from_type_str("y-H20"),
        FragmentIonType::YMinusH2O
    );
    assert_eq!(
        FragmentIonType::from_type_str("Y"),
        FragmentIonType::Unknown
    );
    assert_eq!(FragmentIonType::Unknown.as_str(), "unknown");
    assert_eq!(FragmentIonType::Empty as u8, 255);
    assert_eq!(FragmentIonType::Unknown as u8, 0);
}

/// Derived from the header: flags default, annotation, sentinels, identity.
#[test]
fn light_transition_accessors() {
    let flags = TransitionFlags::default();
    assert!(!flags.decoy() && flags.detecting() && flags.quantifying() && !flags.identifying());
    let mut t = LightTransition::default();
    assert_eq!(t.annotation(), "");
    assert!(!t.is_precursor_im_set());
    assert!(!t.is_product_charge_state_set());
    assert_eq!(t.fragment_nr, -1);
    t.set_fragment_type("y");
    assert_eq!(t.annotation(), "y");
    t.fragment_nr = 6;
    assert_eq!(t.annotation(), "y6");
    t.fragment_charge = 2;
    assert_eq!(t.annotation(), "y6^2");
    assert_eq!(t.product_charge_state(), 2);
    t.set_fragment_type("prec");
    assert_eq!(t.annotation(), "prec");
    t.set_fragment_type("q");
    assert_eq!(t.annotation(), "unknown6^2");
    t.set_decoy(true);
    t.set_detecting_transition(false);
    t.set_identifying_transition(true);
    t.set_quantifying_transition(false);
    assert!(t.decoy() && !t.is_detecting_transition());
    assert!(t.is_identifying_transition() && !t.is_quantifying_transition());
    t.precursor_im = 0.9;
    assert!(t.is_precursor_im_set());
    t.transition_name = "tr".into();
    t.peptide_ref = "pep".into();
    assert_eq!(t.native_id(), "tr");
    assert_eq!(t.compound_ref(), "pep");
    let other = LightTransition {
        transition_name: "tr".into(),
        product_mz: 5.0,
        ..LightTransition::default()
    };
    assert_eq!(t, other);
    let set: HashSet<LightTransition> = [t, other].into_iter().collect();
    assert_eq!(set.len(), 1);
    let m1 = LightModification {
        location: 1,
        unimod_id: 35,
    };
    assert_ne!(
        m1,
        LightModification {
            location: 1,
            unimod_id: 4
        }
    );
}

/// Derived from the header: compound defaults and lookups.
#[test]
fn light_compound_and_lookup() {
    let c = LightCompound::default();
    assert_eq!(c.drift_time(), -1.0);
    assert!(c.rt.is_nan() && c.rt_start.is_nan() && c.rt_end.is_nan());
    assert!(c.is_peptide());
    let mut exp = LightTargetedExperiment::new();
    let mk = |id: &str, seq: &str| LightCompound {
        id: id.into(),
        sequence: seq.into(),
        ..LightCompound::default()
    };
    exp.compounds = vec![mk("a", "FIRST"), mk("b", "B"), mk("a", "LAST")];
    // getCompoundByRef: the map is filled in order, so the last duplicate wins.
    assert_eq!(exp.compound_by_ref("a").unwrap().sequence, "LAST");
    assert_eq!(exp.peptide_by_ref("b").unwrap().sequence, "B");
    assert!(exp.compound_by_ref("missing").is_none());
    assert_eq!(exp.compound_reference_map()["a"].sequence, "LAST");
    // findPeptide: linear search, the first duplicate wins.
    assert_eq!(
        transition_helper::find_peptide(&exp, "a").unwrap().sequence,
        "FIRST"
    );
    assert!(transition_helper::find_peptide(&exp, "zzz").is_none());
}

/// Derived from TransitionHelper.cpp: `convert` groups by compound reference
/// and appends to the caller's map.
#[test]
fn transition_helper_convert() {
    let mut exp = LightTargetedExperiment::new();
    for (name, pep) in [("t1", "p1"), ("t2", "p2"), ("t3", "p1")] {
        exp.transitions.push(LightTransition {
            transition_name: name.into(),
            peptide_ref: pep.into(),
            ..LightTransition::default()
        });
    }
    let mut map: BTreeMap<String, Vec<LightTransition>> = BTreeMap::new();
    map.insert(
        "p1".into(),
        vec![LightTransition {
            transition_name: "t0".into(),
            ..LightTransition::default()
        }],
    );
    transition_helper::convert(&exp, &mut map);
    let names =
        |k: &str| -> Vec<String> { map[k].iter().map(|t| t.transition_name.clone()).collect() };
    assert_eq!(names("p1"), ["t0", "t1", "t3"]);
    assert_eq!(names("p2"), ["t2"]);
}

// ---------------------------------------------------------------------------
// MockObjects
// ---------------------------------------------------------------------------

/// Derived from MockObjects.cpp: each mock returns what its fields hold.
#[test]
fn mock_objects() {
    let feature = MockFeature {
        rt_vec: vec![1.0, 2.0],
        intensity_vec: vec![3.0, 4.0],
        intensity: 5.0,
        rt: 6.0,
    };
    let mut mrm = MockMRMFeature {
        intensity: 7.0,
        rt: 8.0,
        metavalue: 9.0,
        ..MockMRMFeature::default()
    };
    mrm.features.insert("b".into(), Arc::new(feature.clone()));
    mrm.features
        .insert("a".into(), Arc::new(MockFeature::default()));
    mrm.precursor_features
        .insert("p".into(), Arc::new(feature.clone()));
    assert_eq!(mrm.native_ids(), ["a", "b"]);
    assert_eq!(mrm.precursor_ids(), ["p"]);
    assert_eq!(mrm.len(), 2);
    assert_eq!(mrm.meta_value("anything").unwrap(), 9.0);
    let f = mrm.feature("b").unwrap();
    assert_eq!(f.rt_values().unwrap(), vec![1.0, 2.0]);
    assert_eq!(f.intensity_values().unwrap(), vec![3.0, 4.0]);
    assert_eq!((f.intensity(), f.rt()), (5.0, 6.0));
    assert!(mrm.feature("missing").is_none());
    assert_eq!(mrm.native_ids().len(), 2, "a failed lookup inserts nothing");
    assert_eq!(mrm.precursor_feature("p").unwrap().rt(), 6.0);
    let group = MockTransitionGroup {
        size: 3,
        native_ids: vec!["x".into()],
        library_intensities: vec![0.5],
    };
    assert_eq!(group.len(), 3);
    assert_eq!(group.native_ids(), ["x"]);
    assert_eq!(group.library_intensities().unwrap(), vec![0.5]);
    let sn = MockSignalToNoise { sn_value: 4.5 };
    assert_eq!(sn.value_at_rt(-100.0), 4.5);
}

// ---------------------------------------------------------------------------
// StatsHelpers
// ---------------------------------------------------------------------------

/// StatsHelpers_test.cpp:testscorefunction.
#[test]
fn stats_helpers_norm_dot_manhattan() {
    let mz = [
        599.97, 599.98, 599.99, 600.0, 600.01, 600.02, 600.03, 600.97, 600.98, 600.99, 601.0,
        601.01, 601.02, 601.03, 601.97, 601.98, 601.99, 602.0, 602.01, 602.02, 602.03, 602.99,
        603.0, 603.01,
    ];
    let norm = stats_helpers::norm(&mz);
    let normalized = stats_helpers::normalize(&mz, norm);
    real_similar(stats_helpers::norm(&normalized), 1.0);
    real_similar(
        stats_helpers::dot_prod(&normalized, &normalized).unwrap(),
        1.0,
    );
    real_similar(
        stats_helpers::manhattan_dist(&normalized, &normalized).unwrap(),
        0.0,
    );
}

/// Derived from StatsHelpers.cpp: normalisation rules, scoring and guards.
#[test]
fn stats_helpers_derived() {
    assert_eq!(stats_helpers::normalize(&[2.0, 4.0], 0.0), vec![0.0, 0.0]);
    let mut out = vec![7.0, 8.0, 9.0];
    stats_helpers::normalize_into(&[2.0, 4.0], -1.0, &mut out);
    assert_eq!(
        out,
        vec![7.0, 8.0],
        "a non-positive factor leaves the reused output"
    );
    stats_helpers::normalize_into(&[2.0, 4.0], 2.0, &mut out);
    assert_eq!(out, vec![1.0, 2.0]);
    assert_eq!(stats_helpers::norm(&[]), 0.0);
    assert_eq!(stats_helpers::norm(&[3.0, 4.0]), 5.0);
    assert_eq!(
        stats_helpers::dot_prod_i32(&[1, 2], &[3, 4, 99]).unwrap(),
        11.0
    );
    assert!(stats_helpers::dot_prod_i32(&[i32::MAX], &[2]).is_err());
    assert_eq!(stats_helpers::dot_prod_f32(&[1.5], &[2.0]).unwrap(), 3.0);
    assert!(stats_helpers::dot_prod(&[1.0, 2.0], &[1.0]).is_err());
    // Identical vectors: sqrt-normalised dot product 1, Manhattan distance 0.
    let v = vec![1.0, 4.0, 9.0];
    real_similar(
        stats_helpers::dotprod_scoring(v.clone(), v.clone()).unwrap(),
        1.0,
    );
    real_similar(stats_helpers::manhattan_scoring(v.clone(), v).unwrap(), 0.0);
    // [1,4] vs [4,1]: sqrt -> [1,2],[2,1]; unit vectors dot = 4/5.
    real_similar(
        stats_helpers::dotprod_scoring(vec![1.0, 4.0], vec![4.0, 1.0]).unwrap(),
        0.8,
    );
    // sqrt -> [1,2]/3 and [2,1]/3: |1/3-2/3| + |2/3-1/3| = 2/3.
    real_similar(
        stats_helpers::manhattan_scoring(vec![1.0, 4.0], vec![4.0, 1.0]).unwrap(),
        2.0 / 3.0,
    );
    let acc: MeanAndStddev = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]
        .into_iter()
        .collect();
    assert_eq!(acc.count(), 8);
    real_similar(acc.mean(), 5.0);
    real_similar(acc.standard_stddev(), 2.0);
    real_similar(acc.sample_variance(), 32.0 / 7.0);
    real_similar(acc.stddev(), (32.0f64 / 7.0).sqrt());
    assert_eq!(MeanAndStddev::new().variance(), 0.0);
}

/// Derived: the deprecated `cor_pearson` agrees with the crate's Pearson
/// correlation it points to.
#[test]
#[allow(deprecated)]
fn cor_pearson_matches_math() {
    let x = [1.0, 2.0, 3.0, 5.0];
    let y = [2.0, 1.0, 4.0, 3.0];
    let expected =
        openms::math::statistic_functions::pearson_correlation_coefficient(&x, &y).unwrap();
    real_similar(stats_helpers::cor_pearson(&x, &y).unwrap(), expected);
    assert_eq!(
        stats_helpers::cor_pearson(&[1.0, 1.0], &[2.0, 3.0]).unwrap(),
        0.0
    );
}

// ---------------------------------------------------------------------------
// Scoring
// ---------------------------------------------------------------------------

const ARR1: [f64; 6] = [0.0, 1.0, 3.0, 5.0, 2.0, 0.0];
const ARR2: [f64; 6] = [1.0, 3.0, 5.0, 2.0, 0.0, 0.0];

/// Scoring_test.cpp:double_NormalizedManhattanDist_test.
#[test]
fn normalized_manhattan_dist() {
    let (mut a, mut b) = (ARR1.to_vec(), ARR2.to_vec());
    real_similar(
        scoring::normalized_manhattan_dist(&mut a, &mut b).unwrap(),
        0.15151515,
    );
    assert!(scoring::normalized_manhattan_dist(&mut [], &mut []).is_err());
}

/// Scoring_test.cpp:double_RootMeanSquareDeviation_test.
#[test]
fn root_mean_square_deviation() {
    real_similar(
        scoring::root_mean_square_deviation(&ARR1, &ARR2).unwrap(),
        1.91485421551,
    );
    assert!(scoring::root_mean_square_deviation(&[1.0], &[1.0, 2.0]).is_err());
}

/// Scoring_test.cpp:double_SpectralAngle_test.
#[test]
fn spectral_angle() {
    let pi = std::f64::consts::PI;
    let sa = |a: &[f64], b: &[f64]| scoring::spectral_angle(a, b).unwrap();
    real_similar(sa(&ARR1, &ARR2), 0.7699453419277419);
    real_similar(sa(&[0.0, 0.0, 0.0], &[0.0, 0.0, 0.0]), pi / 2.0);
    real_similar(sa(&[1.0, 0.0, 0.0], &[1.0, 0.0, 0.0]), 0.0);
    real_similar(sa(&[1.0, 0.0, 0.0], &[-1.0, 0.0, 0.0]), pi);
    real_similar(sa(&[1.0, 0.0, 0.0], &[0.0, 1.0, 0.0]), pi / 2.0);
    real_similar(
        sa(
            &[0.03174064, 0.11582065, 0.63258941],
            &[0.71882213, 0.00087569, 0.36516896],
        ),
        1.0597217204768459,
    );
    real_similar(
        sa(
            &[0.6608937, 0.0726909, 0.40912141],
            &[0.52081914, 0.71088, 0.0175557],
        ),
        0.9449782659258582,
    );
    real_similar(
        sa(
            &[0.58858475, 0.08963515, 0.08578046],
            &[0.76180969, 0.72763536, 0.50090751],
        ),
        0.6547156284689354,
    );
    real_similar(
        sa(
            &[0.08653022, 0.11595108, 0.74268632],
            &[0.55176333, 0.16783033, 0.70364679],
        ),
        0.5418305329889055,
    );
}

/// Scoring_test.cpp:void_normalize_sum_test.
#[test]
fn normalize_sum() {
    let mut d = ARR1.to_vec();
    scoring::normalize_sum(&mut d);
    for (got, want) in d
        .iter()
        .zip([0.0, 0.09090909, 0.27272727, 0.45454545, 0.18181818, 0.0])
    {
        real_similar(*got, want);
    }
    let mut zeros = vec![0.0, 0.0];
    scoring::normalize_sum(&mut zeros);
    assert_eq!(zeros, vec![0.0, 0.0]);
}

/// Scoring_test.cpp:standardize_data_test.
#[test]
fn standardize_data() {
    let (mut a, mut b) = (ARR1.to_vec(), ARR2.to_vec());
    scoring::standardize_data(&mut a);
    scoring::standardize_data(&mut b);
    let want_a = [
        -1.03479296,
        -0.47036043,
        0.65850461,
        1.78736965,
        0.09407209,
        -1.03479296,
    ];
    let want_b = [
        -0.47036043,
        0.65850461,
        1.78736965,
        0.09407209,
        -1.03479296,
        -1.03479296,
    ];
    for (g, w) in a.iter().zip(want_a).chain(b.iter().zip(want_b)) {
        real_similar(*g, w);
    }
}

fn check_xcorr_five(result: &XCorrArrayType, offset: usize) {
    let want = [0.15634218, 0.8215339, 0.4159292, -0.567846, -0.7374631];
    for (k, w) in want.iter().enumerate() {
        assert_eq!(result.data[offset + k].0, k as i32 - 2);
        real_similar(result.data[offset + k].1, *w);
    }
}

/// Scoring_test.cpp:test_calculateCrossCorrelation.
#[test]
fn calculate_cross_correlation() {
    let (mut a, mut b) = (ARR1.to_vec(), ARR2.to_vec());
    scoring::standardize_data(&mut a);
    scoring::standardize_data(&mut b);
    let mut result = scoring::calculate_cross_correlation(&a, &b, 2, 1).unwrap();
    for entry in result.iter_mut() {
        entry.1 /= 6.0;
    }
    check_xcorr_five(&result, 0);
}

/// Scoring_test.cpp:test_MRMFeatureScoring_normalizedCrossCorrelation and
/// ..._full_range.
#[test]
fn normalized_cross_correlation() {
    let (mut a, mut b) = (ARR1.to_vec(), ARR2.to_vec());
    let result = scoring::normalized_cross_correlation(&mut a, &mut b, 2, 1).unwrap();
    check_xcorr_five(&result, 0);
    let (mut a, mut b) = (ARR1.to_vec(), ARR2.to_vec());
    let result = scoring::normalized_cross_correlation(&mut a, &mut b, 6, 1).unwrap();
    assert_eq!(result.data.len(), 13);
    check_xcorr_five(&result, 4);
}

fn generate_test_signal(n: usize, freq1: f64, freq2: f64) -> Vec<f64> {
    (0..n)
        .map(|i| (i as f64 * freq1).sin() + 0.5 * (i as f64 * freq2).cos())
        .collect()
}

fn xcorr_reference(d1: &[f64], d2: &[f64], maxdelay: i32, lag: i32) -> Vec<(i32, f64)> {
    let n = d1.len() as i32;
    let mut out = Vec::new();
    let mut delay = -maxdelay;
    while delay <= maxdelay {
        let (start, end) = (0.max(-delay), n.min(n - delay));
        let mut sxy = 0.0;
        for i in start..end {
            sxy += d1[i as usize] * d2[(i + delay) as usize];
        }
        out.push((delay, sxy));
        delay += lag;
    }
    out
}

/// Scoring_test.cpp:test_calculateCrossCorrelation_equivalence (the source
/// pins its Eigen version against this left-to-right reference loop).
#[test]
fn cross_correlation_equivalence() {
    let (mut a, mut b) = (ARR1.to_vec(), ARR2.to_vec());
    scoring::standardize_data(&mut a);
    scoring::standardize_data(&mut b);
    let cases = [
        (a, b, 2),
        (
            generate_test_signal(50, 0.1, 0.3),
            generate_test_signal(50, 0.2, 0.5),
            10,
        ),
        (
            generate_test_signal(200, 0.05, 0.15),
            generate_test_signal(200, 0.07, 0.23),
            20,
        ),
    ];
    for (d1, d2, maxdelay) in cases {
        let got = scoring::calculate_cross_correlation(&d1, &d2, maxdelay, 1).unwrap();
        let want = xcorr_reference(&d1, &d2, maxdelay, 1);
        assert_eq!(got.data.len(), want.len());
        for (g, w) in got.data.iter().zip(&want) {
            assert_eq!(g.0, w.0);
            real_similar(g.1, w.1);
        }
    }
}

/// Scoring_test.cpp:test_xcorr_mathematical_properties.
#[test]
fn xcorr_mathematical_properties() {
    let sig = generate_test_signal(50, 0.2, 0.5);
    let (mut s1, mut s2) = (sig.clone(), sig);
    let auto = scoring::normalized_cross_correlation(&mut s1, &mut s2, 10, 1).unwrap();
    let best = scoring::xcorr_array_get_max_peak(&auto).unwrap();
    assert_eq!(auto.data[best].0, 0);
    real_similar(auto.data.iter().find(|e| e.0 == 0).unwrap().1, 1.0);
    let sig1 = generate_test_signal(50, 0.2, 0.5);
    let sig2 = generate_test_signal(50, 0.3, 0.7);
    let ab =
        scoring::normalized_cross_correlation(&mut sig1.clone(), &mut sig2.clone(), 5, 1).unwrap();
    let ba =
        scoring::normalized_cross_correlation(&mut sig2.clone(), &mut sig1.clone(), 5, 1).unwrap();
    for e in &ab {
        let mirror = ba.data.iter().find(|m| m.0 == -e.0).unwrap();
        real_similar(e.1, mirror.1);
    }
}

/// Scoring_test.cpp:test_xcorr_edge_cases.
#[test]
fn xcorr_edge_cases() {
    let result =
        scoring::normalized_cross_correlation(&mut [0.0; 10], &mut [0.0; 10], 10, 1).unwrap();
    assert!(result.data.iter().all(|e| e.1 == 0.0));
    let result =
        scoring::normalized_cross_correlation(&mut [5.0; 10], &mut [5.0; 10], 10, 1).unwrap();
    assert!(result.data.iter().all(|e| e.1.abs() <= 1e-5));
    let d1 = [1.0, 2.0];
    let d2 = [3.0, 4.0];
    let got = scoring::calculate_cross_correlation(&d1, &d2, 1, 1).unwrap();
    let want = xcorr_reference(&d1, &d2, 1, 1);
    for (g, w) in got.data.iter().zip(&want) {
        assert_eq!(g.0, w.0);
        real_similar(g.1, w.1);
    }
    let single = scoring::calculate_cross_correlation(&[42.0], &[7.0], 0, 1).unwrap();
    assert_eq!(single.data, vec![(0, 42.0 * 7.0)]);
}

/// Derived: the preconditions the source compiles out.
#[test]
fn scoring_guards() {
    assert!(scoring::calculate_cross_correlation(&[1.0], &[1.0], 1, 0).is_err());
    assert!(scoring::calculate_cross_correlation(&[1.0], &[1.0, 2.0], 1, 1).is_err());
    assert!(scoring::calculate_cross_correlation(&[1.0], &[1.0], i32::MAX, 1).is_err());
    assert!(
        scoring::calculate_cross_correlation(&[1.0], &[1.0], -1, 1)
            .unwrap()
            .data
            .is_empty()
    );
    let mut empty: [f64; 0] = [];
    assert!(scoring::normalized_cross_correlation(&mut empty, &mut [], 1, 1).is_err());
    assert!(scoring::xcorr_array_get_max_peak(&XCorrArrayType::default()).is_err());
    // Ties keep the first maximum.
    let tie = XCorrArrayType {
        data: vec![(-1, 0.5), (0, 0.9), (1, 0.9)],
    };
    assert_eq!(scoring::xcorr_array_get_max_peak(&tie).unwrap(), 1);
    assert!(scoring::ranked_mutual_information(&[3], &[0], 2, 0).is_err());
    assert!(scoring::ranked_mutual_information(&[], &[], 0, 0).is_err());
}

const RANK_DATA1: [f64; 11] = [
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
const RANK_DATA2: [f64; 11] = [
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

/// Scoring_test.cpp:test_computeAndAppendRank.
#[test]
fn compute_and_append_rank() {
    let mut result = vec![99, 99];
    scoring::compute_and_append_rank(&RANK_DATA1, &mut result).unwrap();
    assert_eq!(result, vec![7, 4, 1, 2, 6, 5, 10, 0, 8, 9, 2]);
}

/// Scoring_test.cpp:test_rankedMutualInformation.
#[test]
fn ranked_mutual_information() {
    let (mut r1, mut r2) = (Vec::new(), Vec::new());
    let max1 = scoring::compute_and_append_rank(&RANK_DATA1, &mut r1).unwrap();
    let max2 = scoring::compute_and_append_rank(&RANK_DATA2, &mut r2).unwrap();
    assert_eq!(max1, *r1.iter().max().unwrap());
    assert_eq!(max2, *r2.iter().max().unwrap());
    let mi = |a: &[u32], b: &[u32], m1: u32, m2: u32| {
        scoring::ranked_mutual_information(a, b, m1, m2).unwrap()
    };
    // The source's TEST_REAL_SIMILAR against 3.2776 (a 5-digit literal).
    real_similar(mi(&r1, &r2, max1, max2), 3.2776);
    real_similar(mi(&[0], &[0], 0, 0), 0.0);
    let a = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    let b = [0, 1, 5, 4, 4, 2, 3, 1, 0, 2];
    real_similar(mi(&a, &b, 0, 5), 0.0);
    real_similar(mi(&a, &b, 0, 5), mi(&b, &a, 5, 0));
    let a = [0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7];
    let b = [0, 1, 5, 4, 4, 2, 3, 1, 0, 2, 6, 7, 7, 6, 5, 3];
    real_similar(mi(&a, &b, 7, 7), 2.0);
    let a = [0, 1, 2, 3, 4, 4, 5, 6, 5, 1];
    let b = [6, 7, 8, 4, 5, 1, 2, 0, 3, 0];
    real_similar(mi(&a, &b, 6, 8), 2.52193);
    real_similar(mi(&a, &b, 6, 8), mi(&b, &a, 8, 6));
}

/// Derived from Scoring.cpp: `computeRankVector` appends one rank vector per row.
#[test]
fn compute_rank_vector_appends() {
    let mut ranks = vec![vec![42u32]];
    let max =
        scoring::compute_rank_vector(&[vec![5.0, 3.0, 3.0, 8.0], vec![1.0]], &mut ranks).unwrap();
    assert_eq!(max, vec![3, 0]);
    assert_eq!(ranks, vec![vec![42], vec![2, 0, 0, 3], vec![0]]);
}

// ---------------------------------------------------------------------------
// DataFrameWriter
// ---------------------------------------------------------------------------

/// Derived from DataFrameWriter.cpp: the in-memory table and the TSV layout.
#[test]
fn data_frame_writers() {
    let mut m = DataMatrix::new();
    m.colnames(&["a".into(), "b".into()]).unwrap();
    m.store("row1", &[1.0, 2.0]).unwrap();
    m.store("row2", &[3.0]).unwrap();
    assert_eq!(m.column_names(), ["a", "b"]);
    assert_eq!(m.row_names(), ["row1", "row2"]);
    assert_eq!(m.rows(), [vec![1.0, 2.0], vec![3.0]]);

    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("out.tsv");
    let mut w = CsvWriter::create(&path).unwrap();
    w.colnames(&["c1".into(), "c2".into(), "c3".into()])
        .unwrap();
    w.store("r", &[1234.5678, 0.000012345, 3.0]).unwrap();
    w.store("empty", &[]).unwrap();
    w.finish().unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    assert_eq!(text, "c1\tc2\tc3\nr\t1234.6\t1.2345e-05\t3\nempty\t\n");
    assert!(CsvWriter::create(dir.path().join("missing/x.tsv")).is_err());
}
