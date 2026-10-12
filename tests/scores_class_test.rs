// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

//! Every assertion of `Scores_test.cpp` at core `bc9cc12`, section by section.

use openms::analysis::scores::{
    ScoreType, all_score_names, is_known_score_type, normalize_score_name,
};
use std::collections::BTreeSet;

#[test]
fn parse_id_type() {
    // type-name strings -> enum (case-insensitive; '-', '_', ' ' ignored)
    for (text, expected) in [
        ("raw", ScoreType::Raw),
        ("RAW", ScoreType::Raw),
        ("raw evalue", ScoreType::RawEValue),
        ("q-value", ScoreType::QValue),
        ("qvalue", ScoreType::QValue),
        ("FDR", ScoreType::Fdr),
        ("false discovery rate", ScoreType::Fdr),
        (
            "Posterior Error Probability",
            ScoreType::PosteriorErrorProbability,
        ),
        ("pep", ScoreType::PosteriorErrorProbability),
        ("Posterior Probability", ScoreType::PosteriorProbability),
        ("pp", ScoreType::PosteriorProbability),
        // a trailing "_score" suffix is stripped before matching
        ("q-value_score", ScoreType::QValue),
    ] {
        assert_eq!(ScoreType::parse(text).unwrap(), expected, "{text}");
    }
    assert!(ScoreType::parse("not-a-score-type").is_err());
}

#[test]
fn is_score_type() {
    assert!(ScoreType::Raw.matches("hyperscore"));
    assert!(ScoreType::Raw.matches("ln(hyperscore)"));
    assert!(ScoreType::QValue.matches("q-value"));
    assert!(ScoreType::QValue.matches("q-value_score"));
    assert!(ScoreType::RawEValue.matches("expect"));
    assert!(ScoreType::Fdr.matches("FDR"));
    assert!(!ScoreType::QValue.matches("hyperscore"));
    assert!(!ScoreType::Raw.matches("unknown_name"));
}

#[test]
fn is_higher_better() {
    assert!(ScoreType::Raw.higher_is_better());
    assert!(ScoreType::PosteriorProbability.higher_is_better());
    assert!(!ScoreType::RawEValue.higher_is_better());
    assert!(!ScoreType::PosteriorErrorProbability.higher_is_better());
    assert!(!ScoreType::Fdr.higher_is_better());
    assert!(!ScoreType::QValue.higher_is_better());
}

#[test]
fn get_all_id_score_names() {
    let names: Vec<_> = all_score_names().collect();
    // 9 RAW + 6 RAW_EVAL + 1 PP + 5 PEP + 3 FDR + 5 QVAL
    assert_eq!(names.len(), 29);
    let set: BTreeSet<_> = names.into_iter().collect();
    for name in [
        "q-value",
        "expect",
        "hyperscore",
        "FDR",
        "Posterior Probability",
    ] {
        assert!(set.contains(name), "{name}");
    }
}

#[test]
fn get_id_names_for_type() {
    assert_eq!(ScoreType::Raw.names().len(), 9);
    assert!(ScoreType::Raw.names().contains(&"hyperscore"));
    assert_eq!(ScoreType::RawEValue.names().len(), 6);
    assert!(ScoreType::RawEValue.names().contains(&"expect"));
    assert_eq!(ScoreType::PosteriorProbability.names().len(), 1);
    assert_eq!(ScoreType::PosteriorErrorProbability.names().len(), 5);
    assert_eq!(ScoreType::Fdr.names().len(), 3);
    assert!(ScoreType::Fdr.names().contains(&"fdr"));
    assert_eq!(ScoreType::QValue.names().len(), 5);
    assert!(ScoreType::QValue.names().contains(&"q-value"));
}

#[test]
fn find_id_type_by_name() {
    assert_eq!(ScoreType::from_name("hyperscore"), Some(ScoreType::Raw));
    assert_eq!(ScoreType::from_name("expect"), Some(ScoreType::RawEValue));
    assert_eq!(ScoreType::from_name("q-value"), Some(ScoreType::QValue));
    assert_eq!(
        ScoreType::from_name("Posterior Error Probability"),
        Some(ScoreType::PosteriorErrorProbability)
    );
    assert_eq!(ScoreType::from_name("definitely_not_a_score"), None);
}

#[test]
fn normalize_score_name_strips_only_the_suffix() {
    assert_eq!(normalize_score_name("q-value_score"), "q-value");
    assert_eq!(normalize_score_name("q-value"), "q-value");
    assert_eq!(normalize_score_name("hyperscore"), "hyperscore");
}

#[test]
fn is_known_score_type_normalizes_first() {
    assert!(is_known_score_type("q-value"));
    assert!(is_known_score_type("q-value_score"));
    assert!(is_known_score_type("hyperscore"));
    assert!(is_known_score_type("expect"));
    assert!(!is_known_score_type("definitely_not_a_score"));
}
