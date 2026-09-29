// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `LabeledPairFinder_test.cpp` and
//! `FeatureGroupingAlgorithmLabeled_test.cpp` at source revision
//! bc9cc12514c768385ce121d6ca4bb710fe1983c4
//! (`src/tests/class_tests/openms/source/`).
//!
//! Evidence tier 3: every expected value is a transcribed class-test literal,
//! except the sections marked "tier 4", whose expected values are derived by
//! hand from the source as shown at the assertion. No expected value was
//! derived from this port's output. The tier-1 differential against the
//! retained `TOPP_FeatureLinkerLabeled_1`/`_2` outputs is
//! `tests/mapmatching_labeled_grouping_topp.rs`.

use openms::analysis::mapmatching::base_group_finder::BaseGroupFinder;
use openms::analysis::mapmatching::feature_grouping_algorithm::FeatureGroupingAlgorithm;
use openms::analysis::mapmatching::feature_grouping_algorithm_labeled::FeatureGroupingAlgorithmLabeled;
use openms::analysis::mapmatching::labeled_pair_finder::LabeledPairFinder;
use openms::concept::UniqueIdGenerator;
use openms::kernel::conversion_helper::MapConversion;
use openms::kernel::features::{ColumnHeader, ConsensusMap, Feature, FeatureMap};
use openms::param::{Param, ParamValue};

/// `TEST_REAL_SIMILAR` with the ClassTest default relative tolerance of 1e-5.
fn similar(actual: f64, expected: f64) {
    let tolerance = 1e-5 * expected.abs().max(actual.abs()).max(1e-5);
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} is not similar to {expected}"
    );
}

/// The ten features both class tests build: `(RT, m/z, intensity)`, all of
/// charge 1 and overall quality 1. Positions are `float` literals in the
/// source (`setRT(1.0f)`), so they are stored as the `f32` values widened.
fn ten_features() -> FeatureMap {
    let values: [(f32, f32, f32); 10] = [
        // start
        (1.0, 1.0, 4.0),
        // best
        (1.5, 5.0, 2.0),
        // inside (down, up, left, right)
        (1.0, 5.0, 0.0),
        (3.0, 5.0, 0.0),
        (1.5, 4.8, 0.0),
        (1.5, 5.2, 0.0),
        // outside (down, up, left, right)
        (0.0, 5.0, 0.0),
        (4.0, 5.0, 0.0),
        (1.5, 4.0, 0.0),
        (1.5, 6.0, 0.0),
    ];
    let features = values
        .iter()
        .map(|&(rt, mz, intensity)| {
            let mut feature = Feature::new(f64::from(rt), f64::from(mz), intensity);
            feature.charge = 1;
            feature.quality = 1.0;
            feature
        })
        .collect();
    FeatureMap::from_features(features)
}

/// The class tests' parameters: no estimation, RT distance 0.4 with
/// deviations 1 (below) and 2 (above), m/z distance 4 with deviation 0.6.
fn manual_param() -> Param {
    let mut p = LabeledPairFinder::defaults().unwrap();
    set(&mut p, "rt_estimate", ParamValue::String("false".into()));
    set(&mut p, "rt_pair_dist", ParamValue::Float(0.4));
    set(&mut p, "rt_dev_low", ParamValue::Float(1.0));
    set(&mut p, "rt_dev_high", ParamValue::Float(2.0));
    set(&mut p, "mz_pair_dists", ParamValue::FloatList(vec![4.0]));
    set(&mut p, "mz_dev", ParamValue::Float(0.6));
    p
}

fn set(param: &mut Param, key: &str, value: ParamValue) {
    let description = param.description(key).unwrap().to_string();
    let tags: Vec<String> = param.tags(key).unwrap().iter().cloned().collect();
    param.set_value(key, value, &description, &tags).unwrap();
}

/// Column headers 5 (light) and 8 (heavy), both of file `filename`.
fn light_heavy_headers(out: &mut ConsensusMap) {
    out.column_headers.insert(
        5,
        ColumnHeader {
            label: "light".into(),
            filename: "filename".into(),
            ..ColumnHeader::default()
        },
    );
    let mut heavy = out.column_headers[&5].clone();
    heavy.label = "heavy".into();
    out.column_headers.insert(8, heavy);
}

// ---------------------------------------------------------------------------
// LabeledPairFinder_test.cpp
// ---------------------------------------------------------------------------

#[test]
fn labeled_pair_finder_defaults() {
    // Tier 4: the source constructor's `defaults_.setValue` calls.
    let finder = LabeledPairFinder::default();
    assert!(finder.rt_estimate);
    assert_eq!(finder.rt_pair_dist, -20.0);
    assert_eq!(finder.rt_dev_low, 15.0);
    assert_eq!(finder.rt_dev_high, 15.0);
    assert_eq!(finder.mz_pair_dists, vec![4.0]);
    assert_eq!(finder.mz_dev, 0.05);
    assert!(!finder.mrm);
    let p = LabeledPairFinder::defaults().unwrap();
    assert_eq!(LabeledPairFinder::from_param(&p).unwrap(), finder);
    assert!(p.tags("mrm").unwrap().iter().any(|t| t == "advanced"));
}

#[test]
fn labeled_pair_finder_run() {
    let mut pm = LabeledPairFinder::from_param(&manual_param()).unwrap();

    let mut output = ConsensusMap::new();
    assert!(pm.run(&[], &mut output).is_err());

    let input = vec![MapConversion::feature_map_to_consensus(5, &ten_features(), None).unwrap()];
    light_heavy_headers(&mut output);
    pm.run(&input, &mut output).unwrap();

    assert_eq!(output.len(), 1);
    let handles = output.features[0].handles();
    similar(handles[0].mz, 1.0);
    similar(handles[0].rt, 1.0);
    let last = handles.last().unwrap();
    similar(last.mz, 5.0);
    similar(last.rt, 1.5);
    similar(f64::from(output.features[0].quality), 0.959346);
    assert_eq!(output.features[0].charge, 1);
}

#[cfg(feature = "featurexml")]
#[test]
fn labeled_pair_finder_rt_estimation() {
    // Test automated RT parameter estimation.
    let mut p2 = LabeledPairFinder::defaults().unwrap();
    set(&mut p2, "rt_estimate", ParamValue::String("true".into()));
    set(&mut p2, "mz_pair_dists", ParamValue::FloatList(vec![4.0]));
    set(&mut p2, "mz_dev", ParamValue::Float(0.2));
    let pm2 = LabeledPairFinder::from_param(&p2).unwrap();

    let features2 =
        openms::format::featurexml::load("tests/data/labeled_pair_finder.featureXML").unwrap();
    let input2 = vec![MapConversion::feature_map_to_consensus(5, &features2, None).unwrap()];
    let mut output2 = ConsensusMap::new();
    light_heavy_headers(&mut output2);
    let mut generator = UniqueIdGenerator::from_seed(19_991_231_235_959);
    let estimate = pm2
        .run_with_generator(&input2, &mut output2, &mut generator)
        .unwrap();
    assert!(estimate.estimated);
    assert_eq!(output2.len(), 250);
}

#[test]
fn labeled_pair_finder_rejects_invalid_headers() {
    // Tier 4: the source's four `Exception::IllegalArgument` checks.
    let pm = LabeledPairFinder::from_param(&manual_param()).unwrap();
    let input = vec![MapConversion::feature_map_to_consensus(5, &ten_features(), None).unwrap()];
    let mut generator = UniqueIdGenerator::from_seed(1);

    // Two input maps.
    let mut out = ConsensusMap::new();
    light_heavy_headers(&mut out);
    let two = vec![input[0].clone(), ConsensusMap::new()];
    assert!(
        pm.run_with_generator(&two, &mut out, &mut generator)
            .is_err()
    );

    // One column header.
    let mut out = ConsensusMap::new();
    out.column_headers.insert(5, ColumnHeader::default());
    assert!(
        pm.run_with_generator(&input, &mut out, &mut generator)
            .is_err()
    );

    // Different file names.
    let mut out = ConsensusMap::new();
    light_heavy_headers(&mut out);
    out.column_headers.get_mut(&8).unwrap().filename = "other".into();
    assert!(
        pm.run_with_generator(&input, &mut out, &mut generator)
            .is_err()
    );

    // No heavy label.
    let mut out = ConsensusMap::new();
    light_heavy_headers(&mut out);
    out.column_headers.get_mut(&8).unwrap().label = "medium".into();
    let before = out.clone();
    assert!(
        pm.run_with_generator(&input, &mut out, &mut generator)
            .is_err()
    );
    assert_eq!(out, before);
}

#[test]
fn labeled_pair_finder_greedy_selection_and_identifications() {
    // Tier 4, derived from the source: for the light feature at RT 1, m/z 1
    // the window is RT [0.4, 3.4] and m/z [4.4, 5.6], which holds five
    // features (m/z 4.8, 5.0, 5.2 at RT 1.5; m/z 5.0 at RT 1.0 and 3.0); no
    // other feature has a partner, so five candidate pairs exist, and the
    // greedy step keeps only the best one, because every candidate shares the
    // light feature. The
    // appended identifications of the input map and the protein
    // identifications are carried to the result.
    let mut features = ten_features();
    for (i, feature) in features.features.iter_mut().enumerate() {
        feature.unique_id = 100 + i as u64;
    }
    let pm = LabeledPairFinder::from_param(&manual_param()).unwrap();
    let mut input = MapConversion::feature_map_to_consensus(5, &features, None).unwrap();
    input
        .unassigned_peptide_identifications
        .push(openms::identification::PeptideIdentification::default());
    let mut out = ConsensusMap::new();
    light_heavy_headers(&mut out);
    let mut generator = UniqueIdGenerator::from_seed(7);
    let estimate = pm
        .run_with_generator(&[input], &mut out, &mut generator)
        .unwrap();
    assert!(!estimate.estimated);
    assert_eq!(estimate.rt_pair_dist, 0.4);
    assert_eq!(out.len(), 1);
    let keys: Vec<_> = out.features[0].handles().iter().map(|h| h.key()).collect();
    assert_eq!(keys, vec![(5, 100), (8, 101)]);
    assert_eq!(out.unassigned_peptide_identifications.len(), 1);
    // Monoisotopic consensus: mean RT and intensity, minimum m/z.
    similar(out.features[0].rt, 1.25);
    similar(out.features[0].mz, 1.0);
    similar(f64::from(out.features[0].intensity), 3.0);
    assert_ne!(out.features[0].unique_id, 0);
}

// ---------------------------------------------------------------------------
// FeatureGroupingAlgorithmLabeled_test.cpp
// ---------------------------------------------------------------------------

#[test]
fn feature_grouping_algorithm_labeled_group() {
    let mut fga = FeatureGroupingAlgorithmLabeled::default();
    let mut out = ConsensusMap::new();

    // Test exception (no input).
    assert!(fga.group(&[], &mut out).is_err());

    // Real test.
    let input = vec![ten_features()];
    fga = FeatureGroupingAlgorithmLabeled::from_param(&manual_param()).unwrap();

    // Test exception (no file name set in out).
    assert!(fga.group(&input, &mut out).is_err());

    light_heavy_headers(&mut out);
    fga.group(&input, &mut out).unwrap();

    // TOLERANCE_ABSOLUTE(0.001) in the source.
    assert_eq!(out.len(), 1);
    assert!((f64::from(out.features[0].quality) - 0.959346).abs() <= 0.001);
    assert_eq!(out.features[0].len(), 2);
    let handles = out.features[0].handles();
    assert!((handles[0].mz - 1.0).abs() <= 0.001);
    assert!((handles[0].rt - 1.0).abs() <= 0.001);
    assert!((handles[1].mz - 5.0).abs() <= 0.001);
    assert!((handles[1].rt - 1.5).abs() <= 0.001);
}

#[test]
fn feature_grouping_algorithm_labeled_parameters() {
    // Tier 4: the source inserts `LabeledPairFinder().getParameters()` at the
    // root of its defaults.
    assert_eq!(
        FeatureGroupingAlgorithmLabeled::defaults().unwrap(),
        LabeledPairFinder::defaults().unwrap()
    );
    let algorithm = FeatureGroupingAlgorithmLabeled::from_param(&manual_param()).unwrap();
    assert_eq!(algorithm.to_param().unwrap(), manual_param());
}
