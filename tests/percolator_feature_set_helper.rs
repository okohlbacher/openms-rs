// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `ANALYSIS/ID/PercolatorFeatureSetHelper.h`.
//!
//! Expectations come from `PercolatorFeatureSetHelper_test.cpp` at `bc9cc12`
//! and its six retained idXML fixtures: each `*.topperc.idXML` input is
//! processed and compared with the `*.topperc_check.idXML` file the C++
//! produced. The class test compares only meta-value key sets and the
//! registered feature list (stored in the check file as the
//! `extra_features` search-parameter meta value); this test also compares
//! every meta value, since the check files retain them. See
//! `tests/data/percolator_provenance.json`.
#![cfg(feature = "idxml")]

use openms::analysis::id::percolator_feature_set_helper::{
    accession_order, add_andes_features, add_comet_features, add_mascot_features,
    add_mascot_features_source, add_msfragger_features, add_msgf_features, add_xtandem_features,
    assign_delta_score, check_extra_features, merge_multise_protein_ids,
    merge_multise_protein_ids_at, rescale_fragment_feature,
};
use openms::chemistry::AASequence;
use openms::format::idxml::{self, IdXmlDocument};
use openms::identification::{
    PeptideHit, PeptideIdentification, ProteinHit, ProteinIdentification,
};
use openms::metadata::{MetaValue, MetaValueData};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/percolator")
        .join(name)
}

fn load(name: &str) -> IdXmlDocument {
    idxml::load(data(name)).unwrap()
}

fn close(a: f64, b: f64) -> bool {
    // The check files print doubles with 15 significant digits.
    a == b || (a - b).abs() <= 1e-12 * a.abs().max(b.abs()).max(1.0)
}

fn same_value(a: &MetaValue, b: &MetaValue) -> bool {
    match (a.data(), b.data()) {
        (MetaValueData::Float(x), MetaValueData::Float(y)) => close(*x, *y),
        _ => a == b,
    }
}

/// The class test's `check_pepids`, strengthened to compare values.
fn check_pepids(check: &[PeptideIdentification], against: &[PeptideIdentification]) {
    assert_eq!(check.len(), against.len());
    for (c, a) in check.iter().zip(against) {
        assert_eq!(c.hits.len(), a.hits.len());
        for (ch, ah) in c.hits.iter().zip(&a.hits) {
            let ck: Vec<_> = ch.metadata.keys().collect();
            let ak: Vec<_> = ah.metadata.keys().collect();
            assert_eq!(ck, ak, "meta keys of {}", ch.sequence.as_str());
            for (key, value) in &ch.metadata {
                assert!(
                    same_value(value, &ah.metadata[key]),
                    "{key} of {}: expected {value:?}, got {:?}",
                    ch.sequence.as_str(),
                    ah.metadata[key]
                );
            }
        }
    }
}

/// The class test's `check_proids`.
fn check_proids(check: &[ProteinIdentification], against: &[ProteinIdentification], fs: &[String]) {
    assert_eq!(check.len(), against.len());
    for (c, a) in check.iter().zip(against) {
        assert_eq!(c.hits.len(), a.hits.len());
    }
    let efc = check[0].search_parameters.metadata["extra_features"]
        .as_str()
        .unwrap();
    assert_eq!(efc, fs.join(","));
}

#[test]
fn merge_multise_protein_ids_unions_by_accession() {
    let mut comet = load("comet.topperc.idXML").protein_identifications;
    let mut msgf = load("msgf.topperc.idXML").protein_identifications;
    let mut expected: BTreeSet<String> = msgf[0].hits.iter().map(|h| h.accession.clone()).collect();
    expected.extend(comet[0].hits.iter().map(|h| h.accession.clone()));
    let expected_db = msgf[0].search_parameters.database.clone();

    let mut merged = Vec::new();
    merge_multise_protein_ids(&mut merged, &mut msgf).unwrap();
    assert_eq!(merged.len(), 1);
    assert_eq!(merged[0].search_engine, "MS-GF+");

    merge_multise_protein_ids(&mut merged, &mut comet).unwrap();
    assert_eq!(merged.len(), 1);
    assert_eq!(merged[0].search_engine, "multiple");
    assert_eq!(merged[0].search_parameters.database, expected_db);
    assert_eq!(merged[0].hits.len(), expected.len());
    let got: BTreeSet<String> = merged[0].hits.iter().map(|h| h.accession.clone()).collect();
    assert_eq!(got, expected);
}

#[test]
fn merge_records_engine_settings_and_run_metadata() {
    // Rust-only: the class test checks neither the identifier, the per-engine
    // settings nor the swap; the expectations follow the source text.
    let mut msgf = load("msgf.topperc.idXML").protein_identifications;
    let incoming = msgf[0].clone();
    let now = chrono::NaiveDate::from_ymd_opt(2026, 9, 28)
        .unwrap()
        .and_hms_opt(10, 11, 12)
        .unwrap();
    let mut merged = Vec::new();
    merge_multise_protein_ids_at(&mut merged, &mut msgf, now).unwrap();
    let run = &merged[0];
    assert_eq!(run.identifier, "TopPerc_2026-09-28");
    assert_eq!(run.date_time.as_deref(), Some("2026-09-28T10:11:12"));
    let meta = &run.search_parameters.metadata;
    assert_eq!(
        meta["SE:MS-GF+"].as_str().unwrap(),
        incoming.search_engine_version
    );
    assert_eq!(
        meta["MS-GF+:db"].as_str().unwrap(),
        incoming.search_parameters.database
    );
    assert!(meta.contains_key("MS-GF+:enzyme_term_specificity"));
    assert!(meta.contains_key("MS-GF+:fragment_mass_tolerance_unit"));
    // The swap leaves the incoming run with the target's former, empty list.
    assert!(msgf[0].hits.is_empty());
    // Sorted by accession.
    assert!(
        run.hits
            .windows(2)
            .all(|w| accession_order(&w[0], &w[1]).is_le())
    );
    // A second run of the same engine is not re-recorded and keeps the name.
    let mut again = vec![incoming];
    let before = merged[0].search_parameters.clone();
    merge_multise_protein_ids_at(&mut merged, &mut again, now).unwrap();
    assert_eq!(merged[0].search_engine, "MS-GF+");
    assert_eq!(merged[0].search_parameters, before);
    // Empty input is an error and changes nothing.
    let snapshot = merged.clone();
    assert!(merge_multise_protein_ids_at(&mut merged, &mut [], now).is_err());
    assert_eq!(merged, snapshot);
}

#[test]
fn set_union_keeps_target_hit_for_equal_accessions() {
    let hit = |accession: &str, score: f64| ProteinHit::new(score, 0, accession, "").unwrap();
    let mut all = vec![ProteinIdentification {
        search_engine: "A".into(),
        hits: vec![hit("P1", 1.0), hit("P3", 1.0)],
        ..Default::default()
    }];
    let mut new = vec![ProteinIdentification {
        search_engine: "A".into(),
        hits: vec![hit("P3", 2.0), hit("P2", 2.0)],
        ..Default::default()
    }];
    merge_multise_protein_ids(&mut all, &mut new).unwrap();
    let got: Vec<_> = all[0]
        .hits
        .iter()
        .map(|h| (h.accession.as_str(), h.score))
        .collect();
    assert_eq!(got, [("P1", 1.0), ("P2", 2.0), ("P3", 1.0)]);
    // The incoming run keeps its hits, sorted.
    let sorted: Vec<_> = new[0].hits.iter().map(|h| h.accession.as_str()).collect();
    assert_eq!(sorted, ["P2", "P3"]);
}

#[test]
fn msgf_features_match_retained_output() {
    let check = load("msgf.topperc_check.idXML");
    let mut doc = load("msgf.topperc.idXML");
    let mut fs = Vec::new();
    add_msgf_features(&mut doc.peptide_identifications, &mut fs);
    check_pepids(&check.peptide_identifications, &doc.peptide_identifications);
    check_proids(
        &check.protein_identifications,
        &doc.protein_identifications,
        &fs,
    );
}

#[test]
fn xtandem_features_match_retained_output() {
    let check = load("xtandem.topperc_check.idXML");
    let mut doc = load("xtandem.topperc.idXML");
    let mut fs = Vec::new();
    add_xtandem_features(&mut doc.peptide_identifications, &mut fs).unwrap();
    check_pepids(&check.peptide_identifications, &doc.peptide_identifications);
    check_proids(
        &check.protein_identifications,
        &doc.protein_identifications,
        &fs,
    );
}

#[test]
fn comet_features_match_retained_output() {
    let check = load("comet.topperc_check.idXML");
    let mut doc = load("comet.topperc.idXML");
    let mut fs = Vec::new();
    add_comet_features(&mut doc.peptide_identifications, &mut fs).unwrap();
    check_pepids(&check.peptide_identifications, &doc.peptide_identifications);
    check_proids(
        &check.protein_identifications,
        &doc.protein_identifications,
        &fs,
    );
}

#[test]
fn feature_errors_leave_input_unchanged() {
    // Rust-only: the source throws ConversionError midway (or dereferences an
    // empty vector); the port refuses atomically.
    let mut doc = load("comet.topperc.idXML");
    doc.peptide_identifications[2].hits[0]
        .metadata
        .remove("MS:1002257");
    let before = doc.peptide_identifications.clone();
    let mut fs = Vec::new();
    assert!(add_comet_features(&mut doc.peptide_identifications, &mut fs).is_err());
    assert_eq!(doc.peptide_identifications, before);
    assert!(fs.is_empty());

    let mut empty: Vec<PeptideIdentification> = Vec::new();
    assert!(add_xtandem_features(&mut empty, &mut fs).is_err());
    let mut doc = load("xtandem.topperc.idXML");
    doc.peptide_identifications[1].hits[0]
        .metadata
        .remove("nextscore");
    let before = doc.peptide_identifications.clone();
    assert!(add_xtandem_features(&mut doc.peptide_identifications, &mut fs).is_err());
    assert_eq!(doc.peptide_identifications, before);
    assert!(fs.is_empty());
}

#[test]
fn andes_features_collect_numeric_prefixed_values() {
    let mut hit = PeptideHit::default();
    let meta = &mut hit.metadata;
    meta.insert("andes:RankScore".into(), MetaValue::try_from(12.5).unwrap());
    meta.insert("andes:NumMatchedMainIons".into(), MetaValue::from(7_i64));
    meta.insert(
        "andes:DeltaRankScore".into(),
        MetaValue::try_from(3.25).unwrap(),
    );
    meta.insert("andes:flag".into(), MetaValue::from("yes"));
    meta.insert("target_decoy".into(), MetaValue::from("target"));
    let pid = PeptideIdentification {
        hits: vec![hit],
        ..Default::default()
    };
    let mut fs = Vec::new();
    add_andes_features(&[pid], &mut fs);
    assert_eq!(fs.len(), 3);
    for name in [
        "andes:RankScore",
        "andes:NumMatchedMainIons",
        "andes:DeltaRankScore",
    ] {
        assert!(fs.iter().any(|f| f == name));
    }
    assert!(!fs.iter().any(|f| f == "andes:flag"));
    assert!(!fs.iter().any(|f| f == "target_decoy"));
    // Byte order, as the source's std::set.
    assert_eq!(
        fs,
        [
            "andes:DeltaRankScore",
            "andes:NumMatchedMainIons",
            "andes:RankScore"
        ]
    );
}

#[test]
fn check_extra_features_keeps_features_on_every_psm() {
    let mut hit_a = PeptideHit::default();
    hit_a
        .metadata
        .insert("on_all".into(), MetaValue::try_from(1.0).unwrap());
    hit_a
        .metadata
        .insert("on_first_only".into(), MetaValue::try_from(2.0).unwrap());
    let mut hit_b = PeptideHit::default();
    hit_b
        .metadata
        .insert("on_all".into(), MetaValue::try_from(3.0).unwrap());
    let psms = [hit_a, hit_b];

    let mut extra: Vec<String> = ["on_all", "on_first_only", "on_none"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let removed = check_extra_features(&psms, &mut extra);
    assert_eq!(extra, ["on_all"]);
    assert_eq!(removed, ["on_first_only", "on_none"]);

    let mut none: Vec<String> = Vec::new();
    check_extra_features(&psms, &mut none);
    assert!(none.is_empty());
}

#[test]
fn msfragger_features_are_registered_in_source_order() {
    let mut fs = vec!["existing".to_string()];
    add_msfragger_features(&mut fs);
    assert_eq!(
        fs,
        [
            "existing",
            "MS:1001330",
            "hyperscore",
            "nextscore",
            "isotope_error"
        ]
    );
}

fn mascot_ids() -> Vec<PeptideIdentification> {
    let hit = |score: f64, seq: &str| {
        let mut h = PeptideHit::new(score, 0, 2, AASequence::parse(seq).unwrap()).unwrap();
        h.metadata
            .insert("MS:1001171".into(), MetaValue::try_from(score).unwrap());
        h
    };
    vec![PeptideIdentification {
        hits: vec![
            hit(20.0, "PEPTIDE"),
            hit(50.0, "PEPM(Oxidation)IDE"),
            hit(35.0, "PEPTIDER"),
        ],
        higher_score_better: true,
        ..Default::default()
    }]
}

#[test]
fn mascot_features_sort_and_annotate() {
    // Independently derived: the source text sorts by score, computes
    // consecutive differences with 0 for the last hit, and flags
    // modifications (bool promoted to an int DataValue).
    let mut ids = mascot_ids();
    let mut fs = Vec::new();
    add_mascot_features(&mut ids, &mut fs).unwrap();
    assert_eq!(fs, ["MS:1001171", "MASCOT:delta_score", "MASCOT:hasMod"]);
    let hits = &ids[0].hits;
    let scores: Vec<f64> = hits.iter().map(|h| h.score).collect();
    assert_eq!(scores, [50.0, 35.0, 20.0]);
    let deltas: Vec<f64> = hits
        .iter()
        .map(|h| h.metadata["MASCOT:delta_score"].as_f64().unwrap())
        .collect();
    assert_eq!(deltas, [15.0, 15.0, 0.0]);
    let mods: Vec<i64> = hits
        .iter()
        .map(|h| h.metadata["MASCOT:hasMod"].as_i64().unwrap())
        .collect();
    assert_eq!(mods, [1, 0, 0]);

    // The source variant: only the sort survives.
    let mut ids = mascot_ids();
    let mut fs = Vec::new();
    add_mascot_features_source(&mut ids, &mut fs).unwrap();
    assert_eq!(fs, ["MS:1001171", "MASCOT:delta_score", "MASCOT:hasMod"]);
    let hits = &ids[0].hits;
    assert_eq!(hits[0].score, 50.0);
    assert!(
        hits.iter()
            .all(|h| !h.metadata.contains_key("MASCOT:delta_score"))
    );
    assert!(
        hits.iter()
            .all(|h| !h.metadata.contains_key("MASCOT:hasMod"))
    );
}

#[test]
fn protected_helpers() {
    // rescaleFragmentFeature_: (1+7)^2 / (1+min(n,7))^2.
    assert_eq!(rescale_fragment_feature(1.0, 7), 1.0);
    assert_eq!(rescale_fragment_feature(1.0, 100), 1.0);
    assert_eq!(rescale_fragment_feature(1.0, 3), 4.0);
    assert_eq!(rescale_fragment_feature(2.0, 0), 128.0);
    // assignDeltaScore_: empty is a no-op; one hit gets 0.
    let mut none: Vec<PeptideHit> = Vec::new();
    assign_delta_score(&mut none, "s", "d").unwrap();
    let mut one = vec![PeptideHit::default()];
    one[0].metadata.insert("s".into(), MetaValue::from(5_i64));
    assign_delta_score(&mut one, "s", "d").unwrap();
    assert_eq!(one[0].metadata["d"].as_f64().unwrap(), 0.0);
    // A missing score is refused without mutation.
    let mut two = vec![one[0].clone(), PeptideHit::default()];
    let before = two.clone();
    assert!(assign_delta_score(&mut two, "s", "x").is_err());
    assert_eq!(two, before);
}
