// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Identification-based map alignment against the upstream class test at core
//! pin `bc9cc12` (evidence tier 3): `MapAlignmentAlgorithmIdentification_test.cpp`
//! with its two idXML inputs copied byte for byte (hashes in
//! `tests/data/map_alignment_provenance.json`), plus independently derived
//! checks of the state, the parameters and the feature-map path.

#![cfg(feature = "idxml")]

use openms::Error;
use openms::analysis::mapmatching::alignment_identification::{
    MapAlignmentAlgorithmIdentification, ScoreFilter,
};
use openms::format::idxml;
use openms::identification::{PeptideHit, PeptideIdentification};
use openms::kernel::{Feature, FeatureMap};
use openms::param::{Param, ParamValue};

fn peptides() -> Vec<Vec<PeptideIdentification>> {
    [
        include_bytes!("data/map_alignment_identification_1.idXML").as_slice(),
        include_bytes!("data/map_alignment_identification_2.idXML").as_slice(),
    ]
    .into_iter()
    .map(|bytes| idxml::read(bytes).unwrap().peptide_identifications)
    .collect()
}

/// ClassTest `isRealSimilar`: absolute difference at most 1e-5, or ratio at
/// most 1 + 1e-5 (`src/testframework/source/CONCEPT/ClassTest.cpp:364` at the
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

fn aligner() -> MapAlignmentAlgorithmIdentification {
    let mut aligner = MapAlignmentAlgorithmIdentification::default();
    // The class test sets an unknown parameter, which the source only warns
    // about.
    let mut params = MapAlignmentAlgorithmIdentification::defaults().unwrap();
    params
        .set_value("peptide_score_threshold", ParamValue::Float(0.0), "", &[])
        .unwrap();
    let warnings = aligner.set_parameters(&params).unwrap();
    assert_eq!(warnings.len(), 1, "{warnings:?}");
    aligner
}

#[test]
fn align_without_with_internal_and_with_external_reference() {
    let mut peptides = peptides();
    let mut aligner = aligner();
    // Section `align`: alignment without reference.
    let result = aligner.align(&peptides, None).unwrap();
    let transforms = result.transformations;
    assert_eq!(transforms.len(), 2);
    assert_eq!(transforms[0].data_points().len(), 10);
    assert_eq!(transforms[1].data_points().len(), 10);
    let mut reference_rts = Vec::new();
    for (a, b) in transforms[0]
        .data_points()
        .iter()
        .zip(transforms[1].data_points())
    {
        // Both RT transforms map to a common RT scale.
        similar(a.y, b.y);
        reference_rts.push(a.x);
    }
    assert_eq!(
        result.info,
        [
            "",
            "Alignment based on:",
            "- 10 data points for sample 1",
            "- 10 data points for sample 2",
            ""
        ]
    );
    // Alignment with internal reference.
    let result = aligner.align(&peptides, Some(0)).unwrap();
    let transforms = result.transformations;
    assert_eq!(transforms.len(), 2);
    assert_eq!(transforms[0].model_type(), "identity");
    assert_eq!(transforms[1].data_points().len(), 10);
    for (point, reference) in transforms[1].data_points().iter().zip(&reference_rts) {
        similar(point.y, *reference);
    }
    assert_eq!(result.info[2], "- 10 data points for sample 1 (reference)");

    // Section `setReference`: alignment with external reference.
    aligner.set_reference(&peptides[0]).unwrap();
    peptides.remove(0);
    let transforms = aligner.align(&peptides, None).unwrap().transformations;
    assert_eq!(transforms.len(), 1);
    assert_eq!(transforms[0].data_points().len(), 10);
    for (point, reference) in transforms[0].data_points().iter().zip(&reference_rts) {
        similar(point.y, *reference);
    }
}

#[test]
fn repeated_align_with_internal_reference_does_not_leak_stale_state() {
    // Section "[EXTRA] repeated align() with internal reference does not leak
    // stale reference state": min_run_occur 3 exceeds the two runs.
    let peptides = peptides();
    let mut aligner = MapAlignmentAlgorithmIdentification::default();
    let mut params = MapAlignmentAlgorithmIdentification::defaults().unwrap();
    params
        .set_value("min_run_occur", ParamValue::Integer(3), "", &[])
        .unwrap();
    aligner.set_parameters(&params).unwrap();
    let first = aligner.align(&peptides, Some(0)).unwrap();
    let second = aligner.align(&peptides, Some(0)).unwrap();
    for result in [&first, &second] {
        assert_eq!(result.transformations.len(), 2);
        assert_eq!(result.transformations[0].model_type(), "identity");
        assert_eq!(result.warnings.len(), 1, "{:?}", result.warnings);
        assert!(result.warnings[0].contains("Using 2 instead."));
    }
    let a = first.transformations[1].data_points();
    let b = second.transformations[1].data_points();
    assert!(!a.is_empty());
    assert_eq!(a, b);
}

#[test]
fn parameters_follow_the_source_layout() {
    let defaults = MapAlignmentAlgorithmIdentification::defaults().unwrap();
    assert_eq!(
        defaults.value("min_score").unwrap(),
        &ParamValue::Float(0.05)
    );
    assert_eq!(defaults.entry("min_run_occur").unwrap().min_int, 2);
    assert_eq!(defaults.entry("max_rt_shift").unwrap().min_float, 0.0);
    assert_eq!(
        defaults.valid_strings("use_adducts").unwrap(),
        ["true".to_owned(), "false".to_owned()]
    );
    let mut aligner = MapAlignmentAlgorithmIdentification::default();
    let mut params = Param::new();
    params
        .set_value("min_run_occur", ParamValue::Integer(1), "", &[])
        .unwrap();
    assert!(matches!(
        aligner.set_parameters(&params),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        aligner.align(&peptides(), Some(2)),
        Err(Error::InvalidValue(_))
    ));
}

fn feature(rt: f64, rts_and_hits: &[(f64, &str, f64)]) -> Feature {
    let mut feature = Feature::new(rt, 500.0, 100.0);
    for &(id_rt, sequence, score) in rts_and_hits {
        let mut id = PeptideIdentification::new();
        id.rt = Some(id_rt);
        id.hits
            .push(PeptideHit::new(score, 1, 2, sequence.parse().unwrap()).unwrap());
        feature.base.peptide_identifications.push(id);
    }
    feature
}

#[test]
fn feature_maps_use_the_closest_identification_and_set_the_score_filter() {
    // Independently derived from the source's feature-map path: with
    // use_feature_rt, the feature RT is recorded under the sequence of its
    // identification closest in RT that passes the score filter; a feature
    // whose identifications all fail it is recorded under the empty
    // sequence, as the source does.
    let map = |shift: f64| {
        FeatureMap::from_features(vec![
            feature(
                10.0 + shift,
                &[(9.0, "PEPTIDE", 5.0), (10.5, "PEPTIDER", 1.0)],
            ),
            feature(20.0 + shift, &[(20.0, "ELVISK", 0.5)]),
            feature(30.0 + shift, &[(30.0, "SAMPLER", 5.0)]),
            feature(40.0 + shift, &[(40.0, "ANOTHERK", 5.0)]),
        ])
    };
    let maps = [map(0.0), map(2.0)];
    let mut aligner = MapAlignmentAlgorithmIdentification::default();
    let mut params = MapAlignmentAlgorithmIdentification::defaults().unwrap();
    for (key, value) in [
        ("use_feature_rt", ParamValue::String("true".into())),
        ("score_cutoff", ParamValue::String("true".into())),
        ("min_score", ParamValue::Float(2.0)),
        ("max_rt_shift", ParamValue::Float(0.0)),
    ] {
        params.set_value(key, value, "", &[]).unwrap();
    }
    aligner.set_parameters(&params).unwrap();
    let result = aligner.align(&maps, None).unwrap();
    assert_eq!(aligner.score_filter(), ScoreFilter::AtLeast);
    let notes: Vec<_> = result.transformations[0]
        .data_points()
        .iter()
        .map(|p| p.note.as_str())
        .collect();
    // PEPTIDER (1.0) fails the cut-off, so the first feature is PEPTIDE even
    // though PEPTIDER is closer; ELVISK fails, so its feature goes under "".
    assert_eq!(notes, ["", "ANOTHERK", "PEPTIDE", "SAMPLER"]);
    let first = &result.transformations[0].data_points()[2];
    assert_eq!((first.x, first.y), (10.0, 11.0));
}
