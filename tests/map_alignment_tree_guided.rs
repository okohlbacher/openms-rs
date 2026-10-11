// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `MapAlignmentAlgorithmTreeGuided` against the upstream class tests at core
//! pin `bc9cc12` (evidence tier 3): `MapAlignmentAlgorithmTreeGuided_test.cpp`
//! with its three featureXML inputs (`tests/data/map_alignment_tree_guided_in*.featureXML`)
//! and the second section of `TransformationModelDefaults_test.cpp`. The
//! guide tree is the class test's literal, as the clustering is not ported;
//! its distances check the similarity. The remaining cases are derived from
//! the pinned source (tier 4).

#![cfg(feature = "featurexml")]

use openms::Error;
use openms::analysis::mapmatching::alignment_tree_guided::{
    MapAlignmentAlgorithmTreeGuided, SeqAndRtList, TREE_GUIDED_MODEL_TYPES,
    peptide_identifications_pearson_similarity,
};
use openms::analysis::transformations::transformation_model_defaults;
use openms::concept::unique_id::UniqueIdGenerator;
use openms::format::featurexml;
use openms::kernel::{BaseFeature, Feature, FeatureMap};
use openms::param::{Param, ParamValue};
use std::collections::BTreeMap;

/// ClassTest `isRealSimilar`: absolute difference at most 1e-5, or ratio at
/// most 1 + 1e-5 (`src/testframework/source/CONCEPT/ClassTest.cpp` at the
/// pin).
fn similar(actual: f64, expected: f64) {
    let ratio = if actual.abs() > expected.abs() {
        actual / expected
    } else {
        expected / actual
    };
    assert!(
        (actual - expected).abs() <= 1e-5 || (ratio - 1.0).abs() <= 1e-5,
        "{actual} !~ {expected}"
    );
}

fn maps() -> Vec<FeatureMap> {
    (0..3)
        .map(|i| {
            featurexml::load(format!(
                "tests/data/map_alignment_tree_guided_in{i}.featureXML"
            ))
            .unwrap()
        })
        .collect()
}

/// The class test's tree: `(0, 2, 1.84834e-04)` then `(0, 1, 0.505752)`.
const TREE: [(usize, usize); 2] = [(0, 2), (0, 1)];

#[test]
fn defaults_keep_the_model_subtree() {
    // TransformationModelDefaults_test.cpp, second section.
    let defaults = MapAlignmentAlgorithmTreeGuided::defaults().unwrap();
    let model = defaults.copy("model:", true).unwrap();
    assert!(
        model
            .source_equal(&transformation_model_defaults("b_spline").unwrap())
            .unwrap()
    );
    assert_eq!(
        defaults.value("model_type").unwrap().as_str().unwrap(),
        "b_spline"
    );
    assert_eq!(
        defaults.entry("model_type").unwrap().valid_strings,
        TREE_GUIDED_MODEL_TYPES.map(String::from).to_vec()
    );
    assert_eq!(
        defaults
            .value("align_algorithm:use_feature_rt")
            .unwrap()
            .as_str()
            .unwrap(),
        "true"
    );
    // The rest of the identification alignment's defaults are its own.
    assert_eq!(
        defaults
            .value("align_algorithm:min_run_occur")
            .unwrap()
            .to_i64()
            .unwrap(),
        2
    );
    // The class test sets the defaults back.
    let mut aligner = MapAlignmentAlgorithmTreeGuided::new().unwrap();
    assert!(
        aligner
            .set_parameters(&aligner.to_param().unwrap())
            .unwrap()
            .is_empty()
    );
    assert_eq!(aligner, MapAlignmentAlgorithmTreeGuided::new().unwrap());
}

#[test]
fn parameters_select_the_model() {
    let mut p = Param::new();
    p.set_value("model_type", ParamValue::String("lowess".into()), "", &[])
        .unwrap();
    p.set_value("model:lowess:span", ParamValue::Float(0.5), "", &[])
        .unwrap();
    let aligner = MapAlignmentAlgorithmTreeGuided::from_param(&p).unwrap();
    assert_eq!(aligner.model_type, "lowess");
    assert_eq!(
        aligner
            .model_param()
            .unwrap()
            .value("span")
            .unwrap()
            .to_f64()
            .unwrap(),
        0.5
    );
    let mut bad = Param::new();
    bad.set_value("model_type", ParamValue::String("spline".into()), "", &[])
        .unwrap();
    assert!(MapAlignmentAlgorithmTreeGuided::from_param(&bad).is_err());
}

#[test]
fn build_tree_inputs_and_similarity() {
    let maps = maps();
    let (lists, ranges) = MapAlignmentAlgorithmTreeGuided::extract_seq_and_rt(&maps).unwrap();
    // Identification counts: the indirect test of extractSeqAndRt_.
    assert_eq!(ranges.len(), 3);
    assert_eq!(ranges[0].len(), 6);
    assert_eq!(ranges[1].len(), 5);
    assert_eq!(ranges[2].len(), 5);
    assert!(ranges.iter().all(|r| r.windows(2).all(|w| w[0] <= w[1])));
    // ClusterHierarchical's distance is `1 - similarity` in float; the first
    // average-linkage merge is the closest pair itself.
    let distance = |i: usize, j: usize| {
        1.0f32 - peptide_identifications_pearson_similarity(&lists[i], &lists[j]).unwrap()
    };
    similar(f64::from(distance(2, 0)), 1.84834e-04);
    // Derived: the second merge joins map 1 to {0, 2} at the Lance-Williams
    // average 0.5 d(0, 1) + 0.5 d(2, 1), as AverageLinkage computes it.
    let merged = 0.5f32 * distance(1, 0) + 0.5f32 * distance(2, 1);
    similar(f64::from(merged), 0.505752);
    assert!(distance(1, 0) > distance(2, 0) && distance(2, 1) > distance(2, 0));
    // The similarity is symmetric.
    assert_eq!(distance(0, 1), distance(1, 0));
}

#[test]
fn similarity_follows_the_source_arithmetic() {
    let list = |entries: &[(&str, &[f64])]| -> SeqAndRtList {
        entries
            .iter()
            .map(|(s, rts)| ((*s).to_owned(), rts.to_vec()))
            .collect()
    };
    let empty = SeqAndRtList::new();
    assert_eq!(
        peptide_identifications_pearson_similarity(&empty, &empty).unwrap(),
        0.0
    );
    // No shared sequence: the source's Pearson coefficient of empty ranges
    // throws InvalidRange.
    let a = list(&[("A", &[1.0])]);
    let b = list(&[("B", &[1.0])]);
    assert!(matches!(
        peptide_identifications_pearson_similarity(&a, &b),
        Err(Error::InvalidRange(_))
    ));
    assert!(peptide_identifications_pearson_similarity(&a, &empty).is_err());
    // "Median" by position of an unsorted list, as Math::median(..., true):
    // [9, 1, 5] gives 1, not 5; perfectly correlated medians give Pearson 1.
    let a = list(&[("A", &[9.0, 1.0, 5.0]), ("B", &[10.0]), ("C", &[20.0])]);
    let b = list(&[("A", &[2.0]), ("B", &[11.0]), ("C", &[21.0])]);
    assert_eq!(
        peptide_identifications_pearson_similarity(&a, &b).unwrap(),
        1.0
    );
    let sorted = list(&[("A", &[1.0, 5.0, 9.0]), ("B", &[10.0]), ("C", &[20.0])]);
    assert!(peptide_identifications_pearson_similarity(&sorted, &b).unwrap() < 1.0);
    // The walk stops at the end of the shorter list: {A, B, D, E} against
    // {A, B, C} takes three steps (A, B, then C < D), so 2 shared of 3.
    let a = list(&[("A", &[1.0]), ("B", &[2.0]), ("D", &[3.0]), ("E", &[4.0])]);
    let b = list(&[("A", &[1.0]), ("B", &[3.0]), ("C", &[5.0])]);
    assert_eq!(
        peptide_identifications_pearson_similarity(&a, &b).unwrap(),
        1.0f32 * 2.0 / 3.0
    );
}

/// The class test's `median_feature_rt_per_seq`: per feature the sequence
/// of the identification closest in RT, the feature RT under it, duplicates
/// removed, and the median.
fn median_feature_rt_per_seq(map: &FeatureMap) -> BTreeMap<String, f64> {
    let mut rts: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    for feature in &map.features {
        let mut best = None;
        let mut best_distance = f64::MAX;
        for id in &feature.peptide_identifications {
            if let Some(hit) = id.hits.first() {
                let d = (id.rt.unwrap() - feature.rt).abs();
                if d < best_distance {
                    best_distance = d;
                    best = Some(hit.sequence.to_string());
                }
            }
        }
        if let Some(sequence) = best {
            rts.entry(sequence).or_default().push(feature.rt);
        }
    }
    rts.into_iter()
        .map(|(s, mut v)| {
            v.sort_by(f64::total_cmp);
            v.dedup();
            let n = v.len();
            let median = if n % 2 == 0 {
                (v[n / 2 - 1] + v[n / 2]) / 2.0
            } else {
                v[n / 2]
            };
            (s, median)
        })
        .collect()
}

#[test]
fn class_test_sections() {
    let mut maps = maps();
    let maps_orig = maps.clone();
    let mut aligner = MapAlignmentAlgorithmTreeGuided::new().unwrap();
    let mut generator = UniqueIdGenerator::from_seed(1);
    let (_, ranges) = MapAlignmentAlgorithmTreeGuided::extract_seq_and_rt(&maps).unwrap();

    // treeGuidedAlignment
    let merged = aligner
        .tree_guided_alignment(&TREE, &mut maps, &ranges, &mut generator)
        .unwrap();
    let map_transformed = &merged.map_transformed;
    assert_eq!(map_transformed.len(), 15);
    assert_eq!(merged.trafo_order[0], 1);
    assert_eq!(merged.trafo_order[2], 0);
    assert_eq!(merged.trafo_order, vec![1, 2, 0]);
    assert_eq!(map_transformed.features[0].unique_id, 20);
    assert_eq!(map_transformed.features[14].unique_id, 14);
    for (i, feature) in map_transformed.features.iter().enumerate() {
        assert_eq!(feature.metadata.contains_key("original_RT"), i >= 5, "{i}");
    }
    // The combined maps are stored at the smaller index; the rest is cleared.
    assert_eq!(maps[0].len(), 15);
    assert!(maps[1].is_empty() && maps[2].is_empty());

    // computeTrafosByOriginalRT
    let result = aligner
        .compute_trafos_by_original_rt(&maps_orig, map_transformed)
        .unwrap();
    let trafos = &result.transformations;
    assert_eq!(trafos.len(), 3);
    let mut has_non_identity_shift = false;
    for (map, trafo) in maps_orig.iter().zip(trafos) {
        assert!(!trafo.data_points().is_empty());
        assert_eq!(trafo.model_type(), "b_spline");
        let medians = median_feature_rt_per_seq(map);
        for point in trafo.data_points() {
            let median = medians.get(&point.note).expect("sequence of the map");
            similar(point.x, *median);
            if (point.x - point.y).abs() > 1e-6 {
                has_non_identity_shift = true;
            }
        }
    }
    assert!(has_non_identity_shift);

    // computeTransformedFeatureMaps
    let mut transformed = maps_orig.clone();
    MapAlignmentAlgorithmTreeGuided::compute_transformed_feature_maps(&mut transformed, trafos)
        .unwrap();
    for map in &transformed {
        for feature in &map.features {
            assert!(feature.metadata.contains_key("original_RT"));
        }
    }
    assert!(
        MapAlignmentAlgorithmTreeGuided::compute_transformed_feature_maps(
            &mut transformed,
            &trafos[..2]
        )
        .is_err()
    );
}

#[test]
fn align_with_tree_runs_the_three_steps() {
    let mut maps = maps();
    let mut aligner = MapAlignmentAlgorithmTreeGuided::new().unwrap();
    let mut generator = UniqueIdGenerator::from_seed(1);
    let result = aligner
        .align_with_tree(&mut maps, &TREE, &mut generator)
        .unwrap();
    assert_eq!(result.transformations.len(), 3);
    for (map, trafo) in maps.iter().zip(&result.transformations) {
        for feature in &map.features {
            let original = feature
                .metadata
                .get("original_RT")
                .unwrap()
                .as_f64()
                .unwrap();
            assert_eq!(feature.rt, trafo.apply(original).unwrap());
        }
    }
    assert!(result.info.iter().any(|l| l == "Alignment based on:"));
}

fn identified(rt: f64, sequence: &str, unique_id: u64) -> Feature {
    use openms::identification::{PeptideHit, PeptideIdentification};
    let mut feature = Feature::from(BaseFeature::new(rt, 500.0, 100.0));
    feature.base.unique_id = unique_id;
    feature
        .base
        .peptide_identifications
        .push(PeptideIdentification {
            rt: Some(rt),
            hits: vec![PeptideHit {
                sequence: sequence.parse().unwrap(),
                ..PeptideHit::default()
            }],
            ..PeptideIdentification::default()
        });
    feature
}

#[test]
fn refusals_and_unique_id_conflicts() {
    let mut aligner = MapAlignmentAlgorithmTreeGuided::new().unwrap();
    let mut generator = UniqueIdGenerator::from_seed(7);
    let mut maps = maps();
    let (_, ranges) = MapAlignmentAlgorithmTreeGuided::extract_seq_and_rt(&maps).unwrap();
    for tree in [[(0, 3)], [(1, 1)]] {
        assert!(
            aligner
                .tree_guided_alignment(&tree, &mut maps, &ranges, &mut generator)
                .is_err()
        );
    }
    assert!(
        aligner
            .tree_guided_alignment(&TREE, &mut maps, &ranges[..2], &mut generator)
            .is_err()
    );
    assert!(
        aligner
            .tree_guided_alignment(&[], &mut [], &[], &mut generator)
            .is_err()
    );
    // A map without identifications.
    let mut unidentified = maps.clone();
    unidentified[1] =
        FeatureMap::from_features(vec![Feature::from(BaseFeature::new(1.0, 500.0, 1.0))]);
    let (_, ranges) = MapAlignmentAlgorithmTreeGuided::extract_seq_and_rt(&unidentified).unwrap();
    let error = aligner
        .tree_guided_alignment(&TREE, &mut unidentified, &ranges, &mut generator)
        .unwrap_err();
    assert_eq!(
        error.to_string(),
        "missing information: FeatureMap originating from 'UNKNOWN' contains no Peptide Identifications. Cannot align!"
    );

    // Two maps whose features share unique IDs: FeatureMap::operator+=
    // replaces the repeated ones and logs how many.
    let sequences = ["AAA", "CCC", "DDD", "EEE", "FFF", "GGG"];
    let mut pair: Vec<FeatureMap> = (0..2)
        .map(|m| {
            FeatureMap::from_features(
                sequences
                    .iter()
                    .enumerate()
                    .map(|(i, s)| identified(10.0 * i as f64 + m as f64, s, i as u64 + 1))
                    .collect(),
            )
        })
        .collect();
    let (_, ranges) = MapAlignmentAlgorithmTreeGuided::extract_seq_and_rt(&pair).unwrap();
    let mut linear = Param::new();
    linear
        .set_value("model_type", ParamValue::String("linear".into()), "", &[])
        .unwrap();
    let mut aligner = MapAlignmentAlgorithmTreeGuided::from_param(&linear).unwrap();
    let merged = aligner
        .tree_guided_alignment(&[(0, 1)], &mut pair, &ranges, &mut generator)
        .unwrap();
    assert!(
        merged
            .info
            .iter()
            .any(|l| l == "Replaced 6 invalid uniqueID's")
    );
    let mut ids: Vec<u64> = merged
        .map_transformed
        .features
        .iter()
        .map(|f| f.unique_id)
        .collect();
    assert_eq!(&ids[..6], &[1, 2, 3, 4, 5, 6]);
    ids.sort_unstable();
    ids.dedup();
    assert_eq!(ids.len(), 12);
    merged.map_transformed.validate().unwrap();
}
