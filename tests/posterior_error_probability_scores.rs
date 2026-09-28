// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tests for the identification-facing members of
//! `MATH/STATISTICS/PosteriorErrorProbabilityModel.h` (core SDK `bc9cc12`):
//! `getScore_`, `transformScore_`, `extractAndTransformScores` and
//! `updateScores`. The class test has no section for them; every expected
//! value here is derived from `PosteriorErrorProbabilityModel.cpp:874-1164`
//! (cited per test), or is the invariant that an updated score is the
//! model's `computeProbability` of the transformed score.

use openms::analysis::id::posterior_error_probability_scores::{
    extract_and_transform_scores, get_score, transform_score, update_scores,
};
use openms::identification::{PeptideHit, PeptideIdentification, ProteinIdentification};
use openms::math::posterior_error_probability::{
    IncorrectComponent, OutlierHandling, PepParameters, PosteriorErrorProbabilityModel,
};
use openms::metadata::MetaValue;

fn hit(score: f64, charge: i32) -> PeptideHit {
    PeptideHit {
        score,
        charge,
        ..PeptideHit::default()
    }
}

fn run(identifier: &str, engine: &str) -> ProteinIdentification {
    ProteinIdentification {
        identifier: identifier.into(),
        search_engine: engine.into(),
        ..ProteinIdentification::default()
    }
}

/// An E-value identification: lower is better, score type "E-Value".
fn evalue_id(identifier: &str, hits: Vec<PeptideHit>) -> PeptideIdentification {
    PeptideIdentification {
        identifier: identifier.into(),
        score_type: "E-Value".into(),
        higher_score_better: false,
        hits,
        ..PeptideIdentification::default()
    }
}

// getScore_ (cpp:874-895)
#[test]
fn get_score_prefers_the_current_type_then_meta_values() {
    let mut h = hit(0.25, 2);
    assert_eq!(get_score(&["E-Value"], &h, "E-Value").unwrap(), 0.25);
    h.metadata
        .insert("expect".into(), MetaValue::try_from(0.5).unwrap());
    h.metadata
        .insert("mvh_score".into(), MetaValue::from(7_i32));
    assert_eq!(get_score(&["EValue", "expect"], &h, "Mascot").unwrap(), 0.5);
    // `<type>_score` is consulted after the plain name, per requested type.
    assert_eq!(get_score(&["mvh"], &h, "x").unwrap(), 7.0);
    // The first requested type that matches wins, even over a later exact one.
    assert_eq!(get_score(&["expect", "x"], &h, "x").unwrap(), 0.5);
    assert!(get_score(&["f-val"], &h, "x").is_err());
    h.metadata.insert("f-val".into(), "text".into());
    assert!(get_score(&["f-val"], &h, "x").is_err());
}

// transformScore_ (cpp:900-956)
#[test]
fn transform_score_per_engine() {
    let with = |name: &str, value: f64| {
        let mut h = hit(1.0, 2);
        h.metadata
            .insert(name.into(), MetaValue::try_from(value).unwrap());
        h
    };
    let close = |a: f64, b: f64| assert!((a - b).abs() < 1e-12, "{a} vs {b}");
    close(
        transform_score("XTANDEM", &with("E-Value", 1e-3), "x").unwrap(),
        3.0,
    );
    close(
        transform_score("SIMTANDEM", &with("E-Value", 1e-2), "x").unwrap(),
        2.0,
    );
    close(
        transform_score("OMSSA", &with("OMSSA", 1e-4), "x").unwrap(),
        4.0,
    );
    close(
        transform_score("MSGFPLUS", &with("MS:1002053", 1e-5), "x").unwrap(),
        5.0,
    );
    close(
        transform_score("MS-GF+", &with("expect", 1e-1), "x").unwrap(),
        1.0,
    );
    close(
        transform_score("COMET", &with("MS:1002257", 1e-6), "x").unwrap(),
        6.0,
    );
    close(
        transform_score("MSFRAGGER", &with("expect", 1e-2), "x").unwrap(),
        2.0,
    );
    close(
        transform_score("SPECTRAST", &with("f-val", 0.25), "x").unwrap(),
        25.0,
    );
    close(
        transform_score("MYRIMATCH", &with("mvh", 12.5), "x").unwrap(),
        12.5,
    );
    close(
        transform_score("SAGE", &with("ln(hyperscore)", 3.5), "x").unwrap(),
        3.5,
    );
    close(
        transform_score("PROSE", &with("hyperscore", 2.5), "x").unwrap(),
        2.5,
    );
    close(
        transform_score("MASCOT", &with("EValue", 1e-3), "x").unwrap(),
        3.0,
    );
    // Issue #740: a Mascot hit scored 0 gives NaN.
    let mut zero = with("EValue", 1e-3);
    zero.score = 0.0;
    assert!(transform_score("MASCOT", &zero, "x").unwrap().is_nan());
    // Listed engines without a transformation are refused, as in the source.
    for engine in ["TIDE-SEARCH", "OPENMS/CONSENSUSID_BEST", "XTandem"] {
        assert!(transform_score(engine, &with("E-Value", 1e-3), "x").is_err());
    }
}

// extractAndTransformScores (cpp:958-1085)
#[test]
fn extract_top_hits_with_target_decoy_split() {
    let proteins = vec![run("r1", "XTandem"), run("r2", "Comet")];
    let peptides = vec![
        // Unsorted hits: sort() puts the lowest E-value first.
        evalue_id("r1", vec![hit(1e-1, 2), hit(1e-4, 2)]),
        evalue_id("r1", vec![hit(1e-2, 3)]),
        evalue_id("r1", vec![hit(0.5, 2)]),
        evalue_id("r1", vec![]),
        // Another run's identifications never reach the XTandem scores.
        PeptideIdentification {
            score_type: "expect".into(),
            ..evalue_id("r2", vec![hit(1e-9, 2)])
        },
    ];
    let scores =
        extract_and_transform_scores(&proteins, &peptides, false, true, true, 0.05).unwrap();
    // COMET has one score (not more than two), so only XTANDEM enters.
    assert_eq!(scores.keys().collect::<Vec<_>>(), ["XTANDEM"]);
    let x = &scores["XTANDEM"];
    let expected: Vec<f64> = [1e-4_f64, 1e-2, 0.5].iter().map(|v| -v.log10()).collect();
    assert_eq!(x[0], expected);
    // Raw top-hit score < 0.05 is a target, otherwise a decoy.
    assert_eq!(x[1], expected[..2]);
    assert_eq!(x[2], expected[2..]);
}

#[test]
fn extract_all_hits_and_split_by_charge() {
    let proteins = vec![run("r1", "XTandem")];
    let peptides = vec![
        evalue_id("r1", vec![hit(1e-1, 2), hit(1e-2, 3)]),
        evalue_id("r1", vec![hit(1e-3, 2), hit(1e-4, 3)]),
        evalue_id("r1", vec![hit(1e-5, 2), hit(1e-6, 3)]),
    ];
    let all = extract_and_transform_scores(&proteins, &peptides, false, false, true, 0.05).unwrap();
    assert_eq!(all["XTANDEM"][0].len(), 6);
    // Without top_hits_only the target and decoy lists stay empty.
    assert!(all["XTANDEM"][1].is_empty() && all["XTANDEM"][2].is_empty());

    let split =
        extract_and_transform_scores(&proteins, &peptides, true, false, false, 0.05).unwrap();
    assert_eq!(split.keys().collect::<Vec<_>>(), ["XTANDEM,2", "XTANDEM,3"]);
    let charge2: Vec<f64> = [1e-1_f64, 1e-3, 1e-5].iter().map(|v| -v.log10()).collect();
    assert_eq!(split["XTANDEM,2"][0], charge2);

    // split_charge without any hit is a precondition violation.
    assert!(extract_and_transform_scores(&proteins, &[], true, true, false, 0.05).is_err());
}

#[test]
fn extract_resolves_consensus_runs_and_refuses_untransformable_engines() {
    let mut consensus = run("c", "OpenMS/ConsensusID_best");
    consensus
        .metadata
        .insert("ConsensusIDBaseSearch".into(), "XTandem:2017.2.1".into());
    let peptides: Vec<_> = (1..=3)
        .map(|i| evalue_id("c", vec![hit(10f64.powi(-i), 2)]))
        .collect();
    let scores =
        extract_and_transform_scores(&[consensus], &peptides, false, true, false, 0.05).unwrap();
    assert_eq!(scores.keys().collect::<Vec<_>>(), ["XTANDEM"]);
    let expected: Vec<f64> = (1..=3).map(|i| -10f64.powi(-i).log10()).collect();
    assert_eq!(scores["XTANDEM"][0], expected);

    // tide-search is listed but transformScore_ has no branch for it.
    let tide = run("t", "tide-search");
    let peptides = vec![evalue_id("t", vec![hit(0.1, 2)])];
    assert!(extract_and_transform_scores(&[tide], &peptides, false, true, false, 0.05).is_err());
}

// updateScores (cpp:1088-1162)
#[test]
fn update_scores_writes_probabilities() {
    let mut fit_scores = vec![
        -0.39, 0.06, 0.12, 0.48, 0.94, 1.01, 1.67, 1.68, 1.76, 1.80, 2.44, 3.25, 3.72, 4.12, 4.28,
        4.60, 4.92, 5.28, 5.53, 6.22,
    ];
    let mut model = PosteriorErrorProbabilityModel::with_parameters(PepParameters {
        number_of_bins: 10,
        incorrectly_assigned: IncorrectComponent::Gumbel,
        ..PepParameters::default()
    });
    assert!(model.fit(&mut fit_scores, OutlierHandling::None).unwrap());

    let proteins = vec![run("r1", "XTandem"), run("r2", "Mascot")];
    let mut peptides = vec![
        evalue_id("r1", vec![hit(1e-3, 2), hit(1e-1, 3)]),
        evalue_id("r2", vec![hit(1e-3, 2)]),
    ];
    let outcome =
        update_scores(&model, "xtandem", 2, false, true, &proteins, &mut peptides).unwrap();
    assert!(!outcome.unable_to_fit_data && !outcome.data_might_not_be_well_fit);
    // Only the charge-2 hit is updated; the old score is kept as a meta value.
    let expected = model.compute_probability(3.0).unwrap();
    assert!((peptides[0].hits[0].score - expected).abs() < 1e-12);
    assert_eq!(
        peptides[0].hits[0].metadata["E-Value_score"]
            .as_f64()
            .unwrap(),
        1e-3
    );
    assert_eq!(peptides[0].hits[1].score, 1e-1);
    assert!(!peptides[0].hits[1].metadata.contains_key("E-Value_score"));
    // Another engine's run is untouched.
    assert_eq!(peptides[1].hits[0].score, 1e-3);

    // prob_correct stores one minus the error probability.
    let mut peptides = vec![evalue_id("r1", vec![hit(1e-3, 2)])];
    update_scores(&model, "XTandem", 0, true, false, &proteins, &mut peptides).unwrap();
    assert!((peptides[0].hits[0].score - (1.0 - expected)).abs() < 1e-12);

    // A Mascot score of 0 transforms to NaN and becomes probability 1.
    let mut mascot = vec![PeptideIdentification {
        identifier: "r2".into(),
        score_type: "Mascot".into(),
        hits: vec![{
            let mut h = hit(0.0, 2);
            h.metadata
                .insert("EValue".into(), MetaValue::try_from(1e-3).unwrap());
            h
        }],
        ..PeptideIdentification::default()
    }];
    update_scores(&model, "Mascot", 0, false, false, &proteins, &mut mascot).unwrap();
    assert_eq!(mascot[0].hits[0].score, 1.0);
}
