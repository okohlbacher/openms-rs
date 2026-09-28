// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

//! `IDDecoyProbability` against the executed C++ Release build.
//!
//! The upstream class test has no assertions (every section is a TODO) and no
//! TOPP test runs the tool, so the expected values come from an oracle:
//! `oracle/I1-iddecoy/driver.cpp`, built against the Linux x86_64 Release
//! install `openms4-release-bc9cc12-c19e494-174b576` and run twice on
//! ibminode06 with identical output (`tests/data/IDDecoyProbability_oracle_ibminode06.tsv`,
//! sha256 `804014f9…56ce1`). The inputs below are the driver's, built from
//! exact integer arithmetic so both languages see the same doubles.

use openms::analysis::id::id_decoy_probability::IDDecoyProbability;
use openms::chemistry::AASequence;
use openms::identification::{PeptideHit, PeptideIdentification};
use openms::metadata::MetaValue;
use openms::param::ParamValue;

fn pep(kind: &str, higher: bool, score: f64, td: Option<&str>) -> PeptideIdentification {
    let mut hit = PeptideHit {
        score,
        sequence: AASequence::parse("PEPTIDE").unwrap(),
        ..PeptideHit::default()
    };
    if let Some(td) = td {
        hit.metadata
            .insert("target_decoy".into(), MetaValue::from(td));
    }
    PeptideIdentification {
        score_type: kind.into(),
        higher_score_better: higher,
        hits: vec![hit],
        ..PeptideIdentification::default()
    }
}

fn target_raw(i: i32) -> f64 {
    10.0 + f64::from((i * 37) % 101) / 5.0
}
fn decoy_raw(j: i32) -> f64 {
    f64::from((j * 53) % 97) / 8.0
}
fn target_e(i: i32) -> f64 {
    f64::from((i * 37) % 101 + 1) / 100000.0
}
fn decoy_e(j: i32) -> f64 {
    f64::from((j * 53) % 97 + 1) / 1000.0
}

/// (score, score type, meta values) per output hit of one oracle case.
fn oracle(case: &str) -> Vec<(f64, String, Vec<(String, String)>)> {
    include_str!("data/IDDecoyProbability_oracle_ibminode06.tsv")
        .lines()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|f| f[0] == case)
        .map(|f| {
            let meta = f[5..]
                .iter()
                .map(|kv| {
                    let (k, v) = kv.split_once('=').unwrap();
                    (k.to_owned(), v.to_owned())
                })
                .collect();
            (f[2].parse().unwrap(), f[3].to_owned(), meta)
        })
        .collect()
}

/// Largest absolute score difference; meta values and score types must match.
fn compare(got: &[PeptideIdentification], case: &str) -> f64 {
    let want = oracle(case);
    assert_eq!(got.len(), want.len(), "{case}");
    let mut worst = 0.0_f64;
    for (id, (score, kind, meta)) in got.iter().zip(&want) {
        assert!(id.higher_score_better);
        assert_eq!(&id.score_type, kind);
        let hit = &id.hits[0];
        worst = worst.max((hit.score - score).abs());
        assert_eq!(hit.metadata.len(), meta.len());
        for (key, value) in meta {
            let v = &hit.metadata[key];
            match value.parse::<f64>() {
                Ok(number) if key != "target_decoy" => assert_eq!(v.as_f64().unwrap(), number),
                _ => assert_eq!(&v.to_string(), value),
            }
        }
    }
    eprintln!("{case}: max |difference| to the Release oracle = {worst:e}");
    worst
}

/// Probabilities agree with the Release build to this absolute difference.
/// The fits reproduce the Release build's Levenberg-Marquardt path; the
/// remaining freedom is the platform `pow`/`exp`/`log10`/`tgamma` rounding.
const TOLERANCE: f64 = 1e-12;

#[test]
fn combined_higher_is_better_matches_release_oracle() {
    let mut ids: Vec<_> = (0..200)
        .map(|i| pep("XTandem", true, target_raw(i), Some("target")))
        .chain((0..150).map(|j| pep("XTandem", true, decoy_raw(j), Some("decoy"))))
        .collect();
    IDDecoyProbability::default().apply(&mut ids).unwrap();
    let worst = compare(&ids, "A");
    assert!(worst <= TOLERANCE, "A: max difference {worst:e}");
}

#[test]
fn combined_e_values_match_release_oracle() {
    let mut ids: Vec<_> = (0..200)
        .map(|i| pep("E-Value", false, target_e(i), Some("target")))
        .chain((0..150).map(|j| pep("E-Value", false, decoy_e(j), Some("decoy"))))
        .collect();
    IDDecoyProbability::default().apply(&mut ids).unwrap();
    let worst = compare(&ids, "B");
    assert!(worst <= TOLERANCE, "B: max difference {worst:e}");
}

#[test]
fn zero_e_value_is_refused_where_the_source_writes_nan() {
    // Oracle case B0: the Release build scores the zero E-value as -nan.
    assert!(oracle("B0").iter().any(|(score, _, _)| score.is_nan()));
    let mut ids: Vec<_> = (0..200)
        .map(|i| pep("E-Value", false, target_e(i), Some("target")))
        .chain(std::iter::once(pep("E-Value", false, 0.0, Some("target"))))
        .chain((0..150).map(|j| pep("E-Value", false, decoy_e(j), Some("decoy"))))
        .collect();
    let before = ids.clone();
    assert!(IDDecoyProbability::default().apply(&mut ids).is_err());
    assert_eq!(ids, before);
}

#[test]
fn forward_reverse_with_twenty_bins_matches_release_oracle() {
    let fwd: Vec<_> = (0..200)
        .map(|i| pep("XTandem", true, target_raw(i), None))
        .collect();
    let rev: Vec<_> = (0..150)
        .map(|j| pep("XTandem", true, decoy_raw(j), None))
        .collect();
    let algorithm = IDDecoyProbability {
        number_of_bins: 20,
        ..IDDecoyProbability::default()
    };
    let out = algorithm.apply_forward_reverse(&fwd, &rev).unwrap();
    let worst = compare(&out, "C");
    assert!(worst <= TOLERANCE, "C: max difference {worst:e}");
}

#[test]
fn defaults_match_the_source_declarations() {
    let d = IDDecoyProbability::defaults().unwrap();
    assert_eq!(d.value("number_of_bins").unwrap(), &ParamValue::Integer(40));
    assert_eq!(
        d.value("lower_score_better_default_value_if_zero").unwrap(),
        &ParamValue::Float(50.0)
    );
    assert!(d.has_tag("number_of_bins", "advanced").unwrap());
    let p = IDDecoyProbability {
        number_of_bins: 12,
        lower_score_better_default_value_if_zero: 30.0,
    };
    assert_eq!(
        IDDecoyProbability::from_param(&p.to_param().unwrap()).unwrap(),
        p
    );
}
