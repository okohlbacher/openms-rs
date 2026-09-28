// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

//! Class tests of `ConsensusIDAlgorithm{Average,Best,Worst,Ranks}` at core
//! `bc9cc12`, and the `ConsensusID` TOPP tests 1, 4, 5 and 6 (the only
//! upstream evidence for PEPMatrix and PEPIons, whose class tests are
//! `NOT_TESTABLE`), run on the single spectrum group of
//! `ConsensusID_1_input.idXML`.

use openms::analysis::id::consensus_id::{
    ConsensusIDAlgorithm, ConsensusIDAlgorithmAverage, ConsensusIDAlgorithmBest,
    ConsensusIDAlgorithmPEPIons, ConsensusIDAlgorithmPEPMatrix, ConsensusIDAlgorithmRanks,
    ConsensusIDAlgorithmWorst, ConsensusIDFilter,
};
use openms::analysis::id::needleman_wunsch::ScoringMatrix;
use openms::chemistry::AASequence;
use openms::identification::{PeptideHit, PeptideIdentification};
use openms::param::{Param, ParamValue};

fn hit(sequence: &str, score: f64) -> PeptideHit {
    PeptideHit {
        sequence: AASequence::parse(sequence).unwrap(),
        score,
        ..PeptideHit::default()
    }
}

/// The three runs of the class tests (`ConsensusIDAlgorithmBest_test.cpp`
/// lines 40-90; identical in the Average, Worst and Ranks tests).
fn class_test_ids() -> Vec<PeptideIdentification> {
    let run = |hits: &[(&str, f64)]| PeptideIdentification {
        score_type: "Posterior Error Probability".into(),
        higher_score_better: false,
        hits: hits.iter().map(|(s, v)| hit(s, *v)).collect(),
        ..PeptideIdentification::default()
    };
    vec![
        run(&[("A", 0.1), ("B", 0.2), ("C", 0.3), ("D", 0.4), ("E", 0.5)]),
        run(&[("C", 0.2), ("A", 0.4), ("B", 0.6)]),
        run(&[
            ("F", 0.0),
            ("C", 0.1),
            ("G", 0.2),
            ("D", 0.3),
            ("B", 0.4),
            ("E", 0.5),
            ("H", 0.6),
            ("I", 0.7),
            ("J", 0.8),
            ("K", 0.9),
        ]),
    ]
}

fn with_considered_hits<T, F: Fn(&Param) -> T>(from: F, hits: i64) -> T {
    let mut param = Param::new();
    param
        .set_value("filter:considered_hits", ParamValue::Integer(hits), "", &[])
        .unwrap();
    from(&param)
}

/// Sequence and score of every output hit; the class tests use an absolute
/// tolerance of 0.01 (`TOLERANCE_ABSOLUTE(0.01)`).
fn check(ids: &[PeptideIdentification], expected: &[(&str, f64)]) {
    assert_eq!(ids.len(), 1);
    let hits = &ids[0].hits;
    assert_eq!(hits.len(), expected.len());
    for (hit, (sequence, score)) in hits.iter().zip(expected) {
        assert_eq!(hit.sequence.to_string(), *sequence);
        assert!(
            (hit.score - score).abs() <= 0.01,
            "{sequence}: {} vs {score}",
            hit.score
        );
    }
}

fn support(hit: &PeptideHit) -> f64 {
    hit.metadata["consensus_support"].as_f64().unwrap()
}

#[test]
fn best_class_test() {
    let consensus = with_considered_hits(|p| ConsensusIDAlgorithmBest::from_param(p).unwrap(), 0);
    let mut f = class_test_ids();
    consensus.apply(&mut f, 0).unwrap();
    check(
        &f,
        &[
            ("F", 0.0),
            ("A", 0.1),
            ("C", 0.1),
            ("B", 0.2),
            ("G", 0.2),
            ("D", 0.3),
            ("E", 0.5),
            ("H", 0.6),
            ("I", 0.7),
            ("J", 0.8),
            ("K", 0.9),
        ],
    );
    let mut ids = class_test_ids();
    ids[2].higher_score_better = true;
    let before = ids.clone();
    assert!(consensus.apply(&mut ids, 0).is_err());
    assert_eq!(ids, before, "an error leaves the input unchanged");
}

#[test]
fn best_support_edge_cases_issue_2020() {
    // ConsensusIDAlgorithmBest_test.cpp lines 149-179.
    let mut consensus = ConsensusIDAlgorithmBest::default();
    let mut id = vec![PeptideIdentification::default()];
    id[0].hits = vec![hit("PEPTIDE", 0.0), hit("PEPTIDE", 0.0)];
    consensus.apply(&mut id, 2).unwrap(); // two runs, but one produced no hits
    assert_eq!(id.len(), 1);
    assert_eq!(id[0].hits.len(), 1);
    assert_eq!(id[0].hits[0].sequence.to_string(), "PEPTIDE");
    assert!((support(&id[0].hits[0]) - 1.0).abs() <= 1e-5);

    let mut param = Param::new();
    param
        .set_value(
            "filter:count_empty",
            ParamValue::String("true".into()),
            "",
            &[],
        )
        .unwrap();
    consensus = ConsensusIDAlgorithmBest::from_param(&param).unwrap();
    let again = id[0].hits[0].clone();
    id[0].hits.push(again);
    consensus.apply(&mut id, 2).unwrap();
    assert_eq!(id.len(), 1);
    assert_eq!(id[0].hits.len(), 1);
    assert_eq!(id[0].hits[0].sequence.to_string(), "PEPTIDE");
    assert!(support(&id[0].hits[0]).abs() <= 1e-5);
}

#[test]
fn worst_class_test() {
    let consensus = with_considered_hits(|p| ConsensusIDAlgorithmWorst::from_param(p).unwrap(), 0);
    let mut f = class_test_ids();
    consensus.apply(&mut f, 0).unwrap();
    check(
        &f,
        &[
            ("F", 0.0),
            ("G", 0.2),
            ("C", 0.3),
            ("A", 0.4),
            ("D", 0.4),
            ("E", 0.5),
            ("B", 0.6),
            ("H", 0.6),
            ("I", 0.7),
            ("J", 0.8),
            ("K", 0.9),
        ],
    );
    let mut ids = class_test_ids();
    ids[2].higher_score_better = true;
    assert!(consensus.apply(&mut ids, 0).is_err());
}

#[test]
fn average_class_test() {
    let consensus =
        with_considered_hits(|p| ConsensusIDAlgorithmAverage::from_param(p).unwrap(), 5);
    let mut f = class_test_ids();
    consensus.apply(&mut f, 0).unwrap();
    check(
        &f,
        &[
            ("F", 0.0),
            ("C", 0.2),
            ("G", 0.2),
            ("A", 0.25),
            ("D", 0.35),
            ("B", 0.4),
            ("E", 0.5),
        ],
    );
    // The two "0.2" scores are not equal (floating-point effects), so the
    // ranks of the hits differ.
    assert!(f[0].hits[1].score < f[0].hits[2].score);
    let mut ids = class_test_ids();
    ids[2].higher_score_better = true;
    assert!(consensus.apply(&mut ids, 0).is_err());
}

#[test]
fn ranks_class_test() {
    let consensus = with_considered_hits(|p| ConsensusIDAlgorithmRanks::from_param(p).unwrap(), 5);
    let mut f = class_test_ids();
    consensus.apply(&mut f, 0).unwrap();
    check(
        &f,
        &[
            ("C", 0.8),
            ("A", 0.6),
            ("B", 0.5333),
            ("F", 0.33333),
            ("D", 0.26666),
            ("G", 0.2),
            ("E", 0.06666),
        ],
    );
}

#[test]
fn defaults_match_the_source_declarations() {
    let d = ConsensusIDAlgorithmPEPMatrix::defaults().unwrap();
    assert_eq!(
        d.value("filter:considered_hits").unwrap(),
        &ParamValue::Integer(0)
    );
    assert_eq!(
        d.value("filter:min_support").unwrap(),
        &ParamValue::Float(0.0)
    );
    assert_eq!(
        d.value("filter:count_empty").unwrap(),
        &ParamValue::String("false".into())
    );
    assert_eq!(
        d.value("filter:keep_old_scores").unwrap(),
        &ParamValue::String("false".into())
    );
    assert_eq!(
        d.value("matrix").unwrap(),
        &ParamValue::String("PAM30MS".into())
    );
    assert_eq!(d.value("penalty").unwrap(), &ParamValue::Integer(5));
    let d = ConsensusIDAlgorithmPEPIons::defaults().unwrap();
    assert_eq!(d.value("mass_tolerance").unwrap(), &ParamValue::Float(0.5));
    assert_eq!(d.value("min_shared").unwrap(), &ParamValue::Integer(2));

    let mut bad = Param::new();
    bad.set_value("filter:min_support", ParamValue::Float(1.5), "", &[])
        .unwrap();
    assert!(ConsensusIDAlgorithmBest::from_param(&bad).is_err());
    let mut bad = Param::new();
    bad.set_value("matrix", ParamValue::String("Identity".into()), "", &[])
        .unwrap();
    assert!(ConsensusIDAlgorithmPEPMatrix::from_param(&bad).is_err());

    let round = ConsensusIDAlgorithmPEPMatrix {
        filter: ConsensusIDFilter {
            considered_hits: 6,
            min_support: 0.25,
            count_empty: true,
            keep_old_scores: true,
        },
        matrix: ScoringMatrix::Identity,
        penalty: 3,
    };
    assert_eq!(
        ConsensusIDAlgorithmPEPMatrix::from_param(&round.to_param().unwrap()).unwrap(),
        round
    );
}

#[test]
fn empty_input_is_left_alone() {
    let mut ids: Vec<PeptideIdentification> = Vec::new();
    ConsensusIDAlgorithmBest::default()
        .apply(&mut ids, 0)
        .unwrap();
    assert!(ids.is_empty());
}

#[cfg(feature = "idxml")]
mod topp {
    use super::*;
    use openms::format::idxml;
    use std::path::Path;

    fn load(name: &str) -> Vec<PeptideIdentification> {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data")
            .join(name);
        idxml::load(path).unwrap().peptide_identifications
    }

    /// `ConsensusID -in ConsensusID_1_input.idXML -algorithm ...`: the three
    /// runs have one identification each, all at the same RT and m/z, so the
    /// tool's grouping yields one group of three and calls
    /// `apply(ids, {}, 3)` (topp `ConsensusID.cpp`, idXML branch).
    ///
    /// The expected files store 15 significant digits, so scores and support
    /// are compared to 1e-13.
    fn topp_case(algorithm: &dyn ConsensusIDAlgorithm, expected: &str) {
        let mut ids = load("ConsensusID_1_input.idXML");
        assert_eq!(ids.len(), 3);
        algorithm.apply(&mut ids, 3).unwrap();
        let expected = load(expected);
        assert_eq!(expected.len(), 1);
        assert_eq!(ids[0].score_type, expected[0].score_type);
        assert_eq!(ids[0].higher_score_better, expected[0].higher_score_better);
        let got = &ids[0].hits;
        let want = &expected[0].hits;
        assert_eq!(got.len(), want.len());
        for (g, w) in got.iter().zip(want) {
            assert_eq!(g.sequence, w.sequence);
            assert_eq!(g.charge, w.charge);
            assert!(
                (g.score - w.score).abs() <= 1e-13,
                "{}: score {} vs {}",
                w.sequence,
                g.score,
                w.score
            );
            assert!(
                (support(g) - support(w)).abs() <= 1e-13,
                "{}: support {} vs {}",
                w.sequence,
                support(g),
                support(w)
            );
        }
    }

    #[test]
    fn topp_consensusid_1_pepmatrix_pam30ms() {
        topp_case(
            &ConsensusIDAlgorithmPEPMatrix::default(),
            "ConsensusID_1_output.idXML",
        );
    }

    #[test]
    fn topp_consensusid_4_pepmatrix_identity_considered_6() {
        let mut algorithm = ConsensusIDAlgorithmPEPMatrix {
            matrix: ScoringMatrix::Identity,
            ..ConsensusIDAlgorithmPEPMatrix::default()
        };
        algorithm.filter.considered_hits = 6;
        topp_case(&algorithm, "ConsensusID_4_output.idXML");
    }

    #[test]
    fn topp_consensusid_5_pepions() {
        topp_case(
            &ConsensusIDAlgorithmPEPIons::default(),
            "ConsensusID_5_output.idXML",
        );
    }

    #[test]
    fn topp_consensusid_6_best_min_support() {
        let mut algorithm = ConsensusIDAlgorithmBest::default();
        algorithm.filter.min_support = 0.5;
        topp_case(&algorithm, "ConsensusID_6_output.idXML");
    }
}
