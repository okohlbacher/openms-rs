// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

//! `IDScoreSwitcherAlgorithm_test.cpp` at core `bc9cc12`, the
//! `IDScoreSwitcher` TOPP tests 1 and 2 (which call `switchScores` on every
//! peptide or protein identification), and invariants derived from the source
//! `switchScores` arithmetic.

use openms::analysis::id::id_score_switcher::{IDScoreSwitcherAlgorithm, IDSwitchResult};
use openms::analysis::scores::ScoreType;
use openms::identification::{PeptideHit, PeptideIdentification};
use openms::metadata::MetaValue;
use openms::param::{Param, ParamValue};

fn hit(score: f64, meta: &[(&str, f64)]) -> PeptideHit {
    let mut hit = PeptideHit {
        score,
        ..PeptideHit::default()
    };
    for (key, value) in meta {
        hit.metadata
            .insert((*key).into(), MetaValue::try_from(*value).unwrap());
    }
    hit
}

fn id(score_type: &str, higher: bool, hits: Vec<PeptideHit>) -> PeptideIdentification {
    PeptideIdentification {
        score_type: score_type.into(),
        higher_score_better: higher,
        hits,
        ..PeptideIdentification::default()
    }
}

#[test]
fn find_score_type_class_test() {
    use IDScoreSwitcherAlgorithm as S;
    let r = S::find_score_type(
        &id("Posterior Error Probability", false, vec![hit(0.05, &[])]),
        ScoreType::PosteriorErrorProbability,
    );
    assert!(r.is_main_score_type);
    assert_eq!(r.score_name, "Posterior Error Probability");

    let r = S::find_score_type(
        &id(
            "XTandem",
            true,
            vec![hit(100.0, &[("pep", 0.01), ("other_score", 0.5)])],
        ),
        ScoreType::PosteriorErrorProbability,
    );
    assert!(!r.is_main_score_type);
    assert_eq!(r.score_name, "pep");

    let r = S::find_score_type(
        &id("Mascot", true, vec![hit(50.0, &[("e_value", 0.001)])]),
        ScoreType::PosteriorErrorProbability,
    );
    assert!(!r.is_main_score_type);
    assert!(r.score_name.is_empty());

    let r = S::find_score_type(
        &id("Mascot", true, vec![hit(75.0, &[("PEP", 0.02)])]),
        ScoreType::PosteriorErrorProbability,
    );
    assert!(!r.is_main_score_type);
    assert_eq!(r.score_name, "PEP");

    let r = S::find_score_type(
        &id(
            "SEQUEST:xcorr",
            true,
            vec![hit(2.5, &[("pep_score", 0.03)])],
        ),
        ScoreType::PosteriorErrorProbability,
    );
    assert!(!r.is_main_score_type);
    assert_eq!(r.score_name, "pep_score");

    let r = S::find_score_type(
        &id("q-value", false, vec![hit(0.02, &[])]),
        ScoreType::QValue,
    );
    assert!(r.is_main_score_type);
    assert_eq!(r.score_name, "q-value");

    let r = S::find_score_type(
        &id("Mascot", true, vec![hit(60.0, &[("qvalue", 0.05)])]),
        ScoreType::QValue,
    );
    assert!(!r.is_main_score_type);
    assert_eq!(r.score_name, "qvalue");

    let r = S::find_score_type(
        &id("XTandem", true, vec![hit(120.0, &[("FDR", 0.01)])]),
        ScoreType::Fdr,
    );
    assert!(!r.is_main_score_type);
    assert_eq!(r.score_name, "FDR");
}

#[test]
fn defaults_match_the_source_declarations() {
    let d = IDScoreSwitcherAlgorithm::defaults().unwrap();
    for key in [
        "new_score",
        "new_score_orientation",
        "new_score_type",
        "old_score",
    ] {
        assert_eq!(
            d.value(key).unwrap(),
            &ParamValue::String(String::new()),
            "{key}"
        );
    }
    assert_eq!(
        d.value("proteins").unwrap(),
        &ParamValue::String("false".into())
    );
    // The source default orientation "" is not a valid string: setting
    // parameters without choosing one fails validation, as in the source.
    assert!(IDScoreSwitcherAlgorithm::from_param(&Param::new()).is_err());
    let mut p = Param::new();
    p.set_value(
        "new_score_orientation",
        ParamValue::String("higher_better".into()),
        "",
        &[],
    )
    .unwrap();
    p.set_value("new_score", ParamValue::String("x".into()), "", &[])
        .unwrap();
    let s = IDScoreSwitcherAlgorithm::from_param(&p).unwrap();
    assert!(s.higher_better());
    assert_eq!(s.effective_new_score_type(), "x");
    assert_eq!(
        IDScoreSwitcherAlgorithm::from_param(&s.to_param().unwrap()).unwrap(),
        s
    );
}

#[test]
fn switch_scores_backs_up_differing_old_scores_with_a_tilde() {
    // Derived from the source switchScores: an existing old-score meta value
    // that differs by more than 1e-6 (relative) diverts the backup to "~".
    let switcher = IDScoreSwitcherAlgorithm {
        new_score: "Percolator_PEP".into(),
        new_score_orientation: Some(false),
        ..IDScoreSwitcherAlgorithm::default()
    };
    let mut pep = id(
        "q-value",
        false,
        vec![
            hit(0.01, &[("Percolator_PEP", 0.2)]),
            hit(0.02, &[("Percolator_PEP", 0.3), ("q-value", 0.02)]),
            hit(0.03, &[("Percolator_PEP", 0.4), ("q-value", 0.5)]),
        ],
    );
    assert_eq!(switcher.switch_scores(&mut pep).unwrap(), 3);
    assert_eq!(pep.score_type, "Percolator_PEP");
    assert!(!pep.higher_score_better);
    let f = |h: &PeptideHit, k: &str| h.metadata.get(k).map(|v| v.as_f64().unwrap());
    assert_eq!(pep.hits[0].score, 0.2);
    assert_eq!(f(&pep.hits[0], "q-value"), Some(0.01));
    assert_eq!(f(&pep.hits[1], "q-value"), Some(0.02));
    assert_eq!(f(&pep.hits[1], "q-value~"), None);
    assert_eq!(f(&pep.hits[2], "q-value"), Some(0.5));
    assert_eq!(f(&pep.hits[2], "q-value~"), Some(0.03));

    let mut missing = id("q-value", false, vec![hit(0.01, &[])]);
    let before = missing.clone();
    assert!(switcher.switch_scores(&mut missing).is_err());
    assert_eq!(missing, before);
}

#[test]
fn switch_to_score_type_and_back() {
    let mut ids = vec![id(
        "XTandem",
        true,
        vec![hit(12.0, &[("Posterior Error Probability_score", 0.001)])],
    )];
    let result = IDScoreSwitcherAlgorithm::switch_to_score_type(&mut ids, "pep").unwrap();
    assert_eq!(
        result,
        IDSwitchResult {
            original_score_name: "XTandem".into(),
            original_score_higher_better: true,
            original_score_type: ScoreType::Raw,
            requested_score_higher_better: false,
            requested_score_type: ScoreType::PosteriorErrorProbability,
            requested_score_name: "Posterior Error Probability".into(),
            score_switched: true,
        }
    );
    assert_eq!(ids[0].hits[0].score, 0.001);
    assert_eq!(
        IDScoreSwitcherAlgorithm::switch_back_score_type(&mut ids, &result).unwrap(),
        1
    );
    assert_eq!(ids[0].score_type, "XTandem");
    assert!(ids[0].higher_score_better);
    assert_eq!(ids[0].hits[0].score, 12.0);
}

/// `IDScoreSwitcherAlgorithm_test.cpp` section `switchToGeneralScoreType`.
///
/// The class test loads `IDScoreSwitcherAlgorithm_test_input.idXML`, whose
/// only peptide identification is transcribed here: the upstream file is
/// line-wrapped inside attribute values (the `xmlns:xsi` URI and a
/// `SearchParameters` attribute name are split), which the C++ reader
/// tolerates and the native idXML reader refuses.
#[test]
fn switch_to_general_score_type_class_test() {
    let mut first = hit(
        0.0,
        &[
            ("consensus_support", 1.0),
            ("Posterior Error Probability_score", 1.544466758592018e-16),
        ],
    );
    first.sequence = openms::chemistry::AASequence::parse("LVTDLTK").unwrap();
    first.charge = 2;
    first
        .metadata
        .insert("target_decoy".into(), MetaValue::from("target"));
    let mut ids = vec![id("q-value", false, vec![first])];
    let mut switcher = IDScoreSwitcherAlgorithm::default();
    let count = switcher
        .switch_to_general_score_type(&mut ids, ScoreType::PosteriorErrorProbability)
        .unwrap();
    assert_eq!(ids[0].score_type, "Posterior Error Probability");
    // Derived from the source switchScores: the new score replaces the main
    // score and the old one is kept under the old score type. The default
    // orientation ("", lower is better) already matches PEP's: no warning.
    assert_eq!(count.hits, 1);
    assert_eq!(count.warnings.len(), 0);
    assert!(!ids[0].higher_score_better);
    assert_eq!(ids[0].hits[0].score, 1.544466758592018e-16);
    assert_eq!(ids[0].hits[0].metadata["q-value"].as_f64().unwrap(), 0.0);
    assert_eq!(switcher.new_score, "Posterior Error Probability_score");
    assert_eq!(switcher.new_score_type, "Posterior Error Probability");
}

#[cfg(feature = "idxml")]
mod files {
    use super::*;
    use openms::format::idxml;
    use std::path::Path;

    fn load(name: &str) -> idxml::IdXmlDocument {
        idxml::load(
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/data")
                .join(name),
        )
        .unwrap()
    }

    fn assert_same_meta(got: &openms::metadata::MetaInfo, want: &openms::metadata::MetaInfo) {
        assert_eq!(
            got.keys().collect::<Vec<_>>(),
            want.keys().collect::<Vec<_>>()
        );
        for (key, value) in want {
            match value.as_f64() {
                Ok(number) => assert_eq!(got[key].as_f64().unwrap(), number, "{key}"),
                Err(_) => assert_eq!(got[key].to_string(), value.to_string(), "{key}"),
            }
        }
    }

    /// TOPP_IDScoreSwitcher_2: `-new_score "Posterior Probability_score"
    /// -new_score_type "Posterior Probability" -new_score_orientation
    /// higher_better -proteins`.
    #[test]
    fn topp_idscoreswitcher_2_proteins() {
        let switcher = IDScoreSwitcherAlgorithm {
            new_score: "Posterior Probability_score".into(),
            new_score_type: "Posterior Probability".into(),
            new_score_orientation: Some(true),
            old_score: String::new(),
            proteins: true,
        };
        let mut proteins = load("IDScoreSwitcher_2_input.idXML").protein_identifications;
        for id in &mut proteins {
            switcher.switch_scores(id).unwrap();
        }
        let want = load("IDScoreSwitcher_2_output.idXML").protein_identifications;
        assert_eq!(proteins.len(), want.len());
        for (g, w) in proteins.iter().zip(&want) {
            assert_eq!(g.score_type, w.score_type);
            assert_eq!(g.higher_score_better, w.higher_score_better);
            assert_eq!(g.hits.len(), w.hits.len());
            for (gh, wh) in g.hits.iter().zip(&w.hits) {
                assert_eq!(gh.accession, wh.accession);
                assert_eq!(gh.score, wh.score);
                assert_same_meta(&gh.metadata, &wh.metadata);
            }
        }
    }
}
