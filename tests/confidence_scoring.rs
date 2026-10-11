// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! ConfidenceScoring.
//!
//! Every section of the upstream class test at the core pin `bc9cc12`
//! (`ConfidenceScoring_test.cpp`, evidence tier 3) is transcribed; tests marked
//! "native" cover behaviour that test does not reach (tier 4).

use openms::Error;
use openms::analysis::openswath::confidence_scoring::{ConfidenceScoring, Glm, RtNorm};
use openms::analysis::openswath::reaction_monitoring_transition::ReactionMonitoringTransition;
use openms::analysis::openswath::targeted_experiment::TargetedExperiment;
use openms::analysis::openswath::targeted_experiment_helper::{Peptide, RetentionTime};
use openms::analysis::transformations::TransformationDescription;
use openms::kernel::{Feature, FeatureMap};

fn real_similar(a: f64, b: f64) {
    if a == b {
        return;
    }
    let absdiff = (a - b).abs();
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{a} !~ {b} (absdiff {absdiff}, ratio {ratio})"
    );
}

fn peptide(id: &str, rt: f64) -> Peptide {
    let mut p = Peptide::default();
    p.id = id.into();
    let mut retention_time = RetentionTime::default();
    retention_time.set_rt(rt);
    p.rts.push(retention_time);
    p
}

fn transition(
    native_id: &str,
    peptide_ref: &str,
    product_mz: f64,
    intensity: f64,
) -> ReactionMonitoringTransition {
    let mut t = ReactionMonitoringTransition::new();
    t.set_native_id(native_id);
    t.set_precursor_mz(400.0);
    t.set_product_mz(product_mz);
    t.set_peptide_ref(peptide_ref);
    t.set_library_intensity(intensity);
    t
}

/// The two-assay library of the upstream test: pep_1 (RT 50) with tr_10 and
/// pep_2 (RT 60) with tr_20, one transition each.
fn library() -> TargetedExperiment {
    let mut library = TargetedExperiment::new();
    library.add_peptide(peptide("pep_1", 50.0));
    library.add_transition(transition("tr_10", "pep_1", 500.0, 500.0));
    library.add_peptide(peptide("pep_2", 60.0));
    library.add_transition(transition("tr_20", "pep_2", 500.0, 500.0));
    library
}

fn feature(peptide_ref: Option<&str>, sub_native_id: Option<&str>) -> Feature {
    let mut f = Feature::default();
    f.rt = 60.0;
    f.quality = -1.0;
    if let Some(r) = peptide_ref {
        f.metadata.insert("PeptideRef".into(), r.into());
    }
    if let Some(id) = sub_native_id {
        f.subordinates.push(subordinate(Some(id)));
    }
    f
}

fn subordinate(native_id: Option<&str>) -> Feature {
    let mut sub = Feature::default();
    sub.intensity = 1.0;
    sub.mz = 500.0;
    if let Some(id) = native_id {
        sub.metadata.insert("native_id".into(), id.into());
    }
    sub
}

fn message(err: &Error) -> &str {
    match err {
        Error::InvalidValue(m) => m,
        other => panic!("unexpected error {other}"),
    }
}

#[test]
fn initialize_and_initialize_glm() {
    // START_SECTION(void initialize(...)) and (void initializeGlm(...))
    let mut scoring = ConfidenceScoring::new(false);
    scoring.initialize(
        TargetedExperiment::new(),
        0,
        0,
        TransformationDescription::default(),
    );
    scoring.initialize_glm(0.0, -1.0, -1.0);
    assert_eq!(
        scoring.glm(),
        Glm {
            intercept: 0.0,
            rt_coef: -1.0,
            int_coef: -1.0
        }
    );
}

#[test]
fn score_map() {
    // START_SECTION((void scoreMap(FeatureMap & features)))
    let mut scoring = ConfidenceScoring::new(true); // test mode
    scoring.initialize(
        TargetedExperiment::new(),
        0,
        0,
        TransformationDescription::default(),
    );
    scoring.initialize_glm(0.0, -1.0, -1.0);
    let mut features = FeatureMap::default();
    // fewer than two assays in the library
    assert!(scoring.score_map(&mut features).is_err());

    features
        .features
        .push(feature(Some("pep_1"), Some("tr_10")));
    features
        .features
        .push(feature(Some("pep_2"), Some("tr_20")));
    scoring.initialize(library(), 0, 0, TransformationDescription::default());
    scoring.score_map(&mut features).unwrap();

    real_similar(f64::from(features.features[0].quality), 0.0);
    real_similar(f64::from(features.features[1].quality), 1.0);
    let metadata = |i: usize, key: &str| features.features[i].metadata[key].as_f64().unwrap();
    // the absolute computed score for each feature
    real_similar(metadata(0, "GLM_score"), 0.0);
    real_similar(metadata(1, "GLM_score"), 0.5);
    // the local fdr score (1 - quality)
    real_similar(metadata(0, "local_FDR"), 1.0);
    real_similar(metadata(1, "local_FDR"), 0.0);
    // the cached retention time range of the library
    let norm = scoring.rt_norm();
    assert_eq!(
        norm,
        RtNorm {
            min_rt: 50.0,
            max_rt: 60.0
        }
    );
}

#[test]
fn extra_test_exceptions() {
    // START_SECTION(([EXTRA] test exceptions))
    let mut scoring = ConfidenceScoring::new(true);
    scoring.initialize(library(), 0, 0, TransformationDescription::default());
    scoring.initialize_glm(0.0, -1.0, -1.0);
    let mut features = FeatureMap::default();
    // no PeptideRef: the feature cannot be mapped to the assay
    features.features.push(feature(None, None));
    features.features.push(feature(None, None));
    let err = scoring.score_map(&mut features).unwrap_err();
    assert_eq!(
        message(&err),
        "Feature does not contain meta value 'PeptideRef' (reference to assay)"
    );
    // after adding the meta value there is still an exception: no subordinates
    features.features[0]
        .metadata
        .insert("PeptideRef".into(), "pep_1".into());
    features.features[1]
        .metadata
        .insert("PeptideRef".into(), "pep_2".into());
    let err = scoring.score_map(&mut features).unwrap_err();
    assert_eq!(
        message(&err),
        "Feature intensities were empty - please provide feature subordinate with intensities"
    );
    // sub-features that cannot be mapped to the transitions (no native_id)
    features.features[0].subordinates = vec![subordinate(None)];
    features.features[1].subordinates = vec![subordinate(None)];
    let err = scoring.score_map(&mut features).unwrap_err();
    assert_eq!(
        message(&err),
        "Did not find a feature for each assay provided - each feature needs to have n subordinates with the meta-value 'native_id' set to the corresponding transition."
    );
    features.features[0].subordinates = vec![subordinate(Some("tr_10"))];
    features.features[1].subordinates = vec![subordinate(Some("tr_20"))];
    scoring.score_map(&mut features).unwrap();
    real_similar(f64::from(features.features[0].quality), 0.0);
    real_similar(f64::from(features.features[1].quality), 1.0);
}

#[test]
fn native_glm_and_rt_norm_formulas() {
    // tier 4: sigmoid(intercept + rt_coef * d^2 + int_coef * dist), derived
    // from the formula in ConfidenceScoring.h.
    let glm = Glm {
        intercept: 1.0,
        rt_coef: -0.5,
        int_coef: 2.0,
    };
    let lm: f64 = 1.0 - 0.5 * 4.0 + 2.0 * 0.25;
    assert_eq!(glm.evaluate(2.0, 0.25), 1.0 / (1.0 + (-lm).exp()));
    assert_eq!(glm.evaluate(0.0, 0.0), 1.0 / (1.0 + (-1.0_f64).exp()));
    let norm = RtNorm {
        min_rt: 10.0,
        max_rt: 30.0,
    };
    assert_eq!(norm.apply(10.0), 0.0);
    assert_eq!(norm.apply(30.0), 100.0);
    assert_eq!(norm.apply(20.0), 50.0);
}

#[test]
fn native_n_transitions_keeps_the_most_intense_peaks() {
    // tier 4: extractIntensities_ erases the lowest intensities until
    // n_transitions remain, then lists the survivors by m/z. One assay with
    // three transitions, feature with the same three at intensities
    // (m/z 500: 5, 600: 1, 700: 3), n_transitions = 2 keeps m/z 500 and 700:
    // normalised (5/8, 3/8); the assay's library intensities (500: 4, 600: 9,
    // 700: 2) are cut to the same two: (4/6, 2/6). Distance
    // |5/8 - 4/6| + |3/8 - 2/6|.
    let mut lib = TargetedExperiment::new();
    lib.add_peptide(peptide("true", 50.0));
    lib.add_peptide(peptide("other", 60.0));
    for (id, mz, intensity) in [("a", 500.0, 4.0), ("b", 600.0, 9.0), ("c", 700.0, 2.0)] {
        lib.add_transition(transition(id, "true", mz, intensity));
    }
    lib.add_transition(transition("o", "other", 500.0, 1.0));
    let mut f = Feature::default();
    f.rt = 50.0;
    f.metadata.insert("PeptideRef".into(), "true".into());
    for (id, mz, intensity) in [("a", 500.0, 5.0_f32), ("b", 600.0, 1.0), ("c", 700.0, 3.0)] {
        let mut sub = Feature::default();
        sub.mz = mz;
        sub.intensity = intensity;
        sub.metadata.insert("native_id".into(), id.into());
        f.subordinates.push(sub);
    }
    let mut features = FeatureMap::from_features(vec![f]);
    let mut scoring = ConfidenceScoring::new(false);
    scoring.initialize(lib, 0, 2, TransformationDescription::default());
    scoring.initialize_glm(0.0, 0.0, -1.0);
    scoring.score_map(&mut features).unwrap();
    let dist = (5.0 / 8.0 - 4.0 / 6.0_f64).abs() + (3.0 / 8.0 - 2.0 / 6.0_f64).abs();
    let expected = 1.0 / (1.0 + dist.exp());
    let got = features.features[0].metadata["GLM_score"].as_f64().unwrap();
    assert!((got - expected).abs() < 1e-12, "{got} vs {expected}");
    // the "other" assay has one transition, fewer than the feature's two, and
    // is skipped as a decoy: no decoys, the local FDR is 0 / 0
    assert!(
        features.features[0].metadata["local_FDR"]
            .as_f64()
            .unwrap()
            .is_nan()
    );
}

#[test]
fn native_decoy_draws_are_reproducible() {
    // tier 4: with n_decoys > 0 the decoy list is shuffled by the
    // seed-0 engine; two scorers give the same result and a repeated call on
    // one scorer continues the stream rather than restarting it.
    let mut lib = TargetedExperiment::new();
    for i in 0..6 {
        lib.add_peptide(peptide(&format!("p{i}"), 10.0 * f64::from(i)));
        lib.add_transition(transition(
            &format!("t{i}"),
            &format!("p{i}"),
            500.0,
            100.0 + f64::from(i),
        ));
    }
    let make_features = || {
        let mut f = Feature::default();
        f.rt = 20.0;
        f.metadata.insert("PeptideRef".into(), "p2".into());
        let mut sub = Feature::default();
        sub.mz = 500.0;
        sub.intensity = 7.0;
        sub.metadata.insert("native_id".into(), "t2".into());
        f.subordinates.push(sub);
        FeatureMap::from_features(vec![f])
    };
    let run = || {
        let mut s = ConfidenceScoring::new(false);
        s.initialize(lib.clone(), 3, 0, TransformationDescription::default());
        s.initialize_glm(0.0, -0.01, -1.0);
        let mut map = make_features();
        s.score_map(&mut map).unwrap();
        map.features[0].metadata["local_FDR"].as_f64().unwrap()
    };
    assert_eq!(run(), run());
    // n_decoys larger than the unrelated assays is clamped to all of them
    let mut s = ConfidenceScoring::new(false);
    s.initialize(lib.clone(), 99, 0, TransformationDescription::default());
    s.initialize_glm(0.0, -0.01, -1.0);
    let mut map = make_features();
    s.score_map(&mut map).unwrap();
    assert_eq!(s.n_decoys(), 0);
}
