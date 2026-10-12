// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `ANALYSIS/ID/IonIdentityMolecularNetworking` against the upstream class
//! test (`IonIdentityMolecularNetworking_test.cpp` at core `bc9cc12`), whose
//! literals are transcribed here, plus native error and edge cases.

use openms::analysis::id::ion_identity_molecular_networking::{
    annotate_consensus_map, store_supplementary_pair_table, write_supplementary_pair_table,
};
use openms::concept::constants::user_param::{
    IIMN_ADDUCT_PARTNERS, IIMN_ANNOTATION_NETWORK_NUMBER, IIMN_BEST_ION, IIMN_LINKED_GROUPS,
    IIMN_ROW_ID,
};
use openms::kernel::{ConsensusFeature, ConsensusMap};
use openms::metadata::MetaValue;
use openms::{Error, Result};

/// The class test's `mkCF` helper.
fn cf(mz: f64, groups: &[&str], best_ion: &str) -> ConsensusFeature {
    let mut feature = ConsensusFeature::new();
    feature.mz = mz;
    feature.rt = 100.0;
    if !groups.is_empty() {
        let list: Vec<String> = groups.iter().map(|g| (*g).to_owned()).collect();
        feature
            .metadata
            .insert(IIMN_LINKED_GROUPS.into(), MetaValue::from(list));
    }
    if !best_ion.is_empty() {
        feature
            .metadata
            .insert(IIMN_BEST_ION.into(), MetaValue::from(best_ion));
    }
    feature
}

fn int(map: &ConsensusMap, i: usize, key: &str) -> i64 {
    map.features[i].metadata[key].as_i64().unwrap()
}

fn text(map: &ConsensusMap, i: usize, key: &str) -> String {
    map.features[i].metadata[key].to_string()
}

fn has(map: &ConsensusMap, i: usize, key: &str) -> bool {
    map.features[i].metadata.contains_key(key)
}

#[test]
fn annotate_documented_scenario() -> Result<()> {
    let mut cm = ConsensusMap::from_features(vec![
        cf(100.0, &["1", "2"], ""),
        cf(122.0, &["2", "3"], ""),
        cf(200.0, &["4", "5"], ""),
    ]);
    annotate_consensus_map(&mut cm)?;
    assert_eq!(int(&cm, 0, IIMN_ROW_ID), 1);
    assert_eq!(int(&cm, 1, IIMN_ROW_ID), 2);
    assert_eq!(int(&cm, 2, IIMN_ROW_ID), 3);
    assert_eq!(int(&cm, 0, IIMN_ANNOTATION_NETWORK_NUMBER), 1);
    assert_eq!(int(&cm, 1, IIMN_ANNOTATION_NETWORK_NUMBER), 1);
    assert_eq!(int(&cm, 2, IIMN_ANNOTATION_NETWORK_NUMBER), 2);
    assert_eq!(text(&cm, 0, IIMN_ADDUCT_PARTNERS), "2");
    assert_eq!(text(&cm, 1, IIMN_ADDUCT_PARTNERS), "1");
    assert!(!has(&cm, 2, IIMN_ADDUCT_PARTNERS));
    for i in 0..3 {
        assert!(!has(&cm, i, IIMN_LINKED_GROUPS));
    }
    Ok(())
}

#[test]
fn annotate_chain() -> Result<()> {
    let mut cm = ConsensusMap::from_features(vec![
        cf(100.0, &["1", "2"], ""),
        cf(120.0, &["2", "3"], ""),
        cf(140.0, &["3", "4"], ""),
    ]);
    annotate_consensus_map(&mut cm)?;
    for i in 0..3 {
        assert_eq!(int(&cm, i, IIMN_ANNOTATION_NETWORK_NUMBER), 1);
    }
    assert_eq!(text(&cm, 0, IIMN_ADDUCT_PARTNERS), "2");
    assert_eq!(text(&cm, 1, IIMN_ADDUCT_PARTNERS), "1;3");
    assert_eq!(text(&cm, 2, IIMN_ADDUCT_PARTNERS), "2");
    Ok(())
}

#[test]
fn feature_without_groups_gets_only_a_row_id() -> Result<()> {
    let mut cm = ConsensusMap::from_features(vec![
        cf(100.0, &["1", "2"], ""),
        cf(122.0, &["2", "3"], ""),
        cf(300.0, &[], ""),
    ]);
    annotate_consensus_map(&mut cm)?;
    assert_eq!(int(&cm, 2, IIMN_ROW_ID), 3);
    assert!(!has(&cm, 2, IIMN_ANNOTATION_NETWORK_NUMBER));
    assert!(!has(&cm, 2, IIMN_ADDUCT_PARTNERS));
    Ok(())
}

/// Component numbers follow the first vertex of each component: a feature
/// without groups between two networks does not take a number, and the
/// network of the later features is numbered after the earlier ones.
/// Independently derived from Boost's `connected_components` over a `vecS`
/// vertex list.
#[test]
fn network_numbers_follow_vertex_order() -> Result<()> {
    let mut cm = ConsensusMap::from_features(vec![
        cf(100.0, &["7"], ""),
        cf(110.0, &[], ""),
        cf(120.0, &["8"], ""),
        cf(130.0, &["7"], ""),
    ]);
    annotate_consensus_map(&mut cm)?;
    assert_eq!(int(&cm, 0, IIMN_ANNOTATION_NETWORK_NUMBER), 1);
    assert_eq!(int(&cm, 2, IIMN_ANNOTATION_NETWORK_NUMBER), 2);
    assert_eq!(int(&cm, 3, IIMN_ANNOTATION_NETWORK_NUMBER), 1);
    assert_eq!(text(&cm, 0, IIMN_ADDUCT_PARTNERS), "4");
    assert_eq!(text(&cm, 3, IIMN_ADDUCT_PARTNERS), "1");
    Ok(())
}

#[test]
fn write_pair_table_matches_class_test() -> Result<()> {
    let mut cm = ConsensusMap::from_features(vec![
        cf(100.0, &["1", "2"], "[M+H]+"),
        cf(122.0, &["2", "3"], "[M+Na]+"),
        cf(200.0, &["4", "5"], "[M+H]+"),
    ]);
    annotate_consensus_map(&mut cm)?;
    let mut out = Vec::new();
    assert!(write_supplementary_pair_table(&cm, &mut out)?);
    let text = String::from_utf8(out).unwrap();
    let lines: Vec<_> = text.lines().filter(|l| !l.is_empty()).collect();
    assert_eq!(lines.len(), 2);
    assert_eq!(lines[0], "ID1,ID2,EdgeType,Score,Annotation");
    assert_eq!(lines[1], "1,2,MS1 annotation,0,[M+H]+ [M+Na]+ dm/z=22.0");
    Ok(())
}

/// A three-feature chain writes each unordered pair once. Scores are the
/// partner counts of both sides minus 2; best ions default to `default`.
#[test]
fn write_pair_table_deduplicates_reverse_edges() -> Result<()> {
    let mut cm = ConsensusMap::from_features(vec![
        cf(100.0, &["1", "2"], ""),
        cf(120.5, &["2", "3"], ""),
        cf(140.0, &["3", "4"], "x,y"),
    ]);
    annotate_consensus_map(&mut cm)?;
    let mut out = Vec::new();
    write_supplementary_pair_table(&cm, &mut out)?;
    let text = String::from_utf8(out).unwrap();
    assert_eq!(
        text,
        "ID1,ID2,EdgeType,Score,Annotation\n\
         1,2,MS1 annotation,1,default default dm/z=20.5\n\
         2,3,MS1 annotation,1,default x_y dm/z=19.5\n"
    );
    Ok(())
}

#[test]
fn unannotated_or_empty_maps_write_nothing() -> Result<()> {
    let mut out = Vec::new();
    assert!(!write_supplementary_pair_table(
        &ConsensusMap::new(),
        &mut out
    )?);
    let cm = ConsensusMap::from_features(vec![cf(1.0, &["1"], "")]);
    assert!(!write_supplementary_pair_table(&cm, &mut out)?);
    assert!(out.is_empty());
    let dir = std::env::temp_dir().join(format!("iimn-{}", std::process::id()));
    std::fs::create_dir_all(&dir)?;
    let path = dir.join("pairs.csv");
    assert!(!store_supplementary_pair_table(&cm, &path)?);
    assert!(!path.exists());
    std::fs::remove_dir_all(&dir)?;
    Ok(())
}

#[test]
fn bad_meta_values_are_errors_and_leave_the_map_unchanged() -> Result<()> {
    let mut cm = ConsensusMap::from_features(vec![cf(100.0, &["1", "abc"], "")]);
    let before = cm.clone();
    assert!(matches!(
        annotate_consensus_map(&mut cm),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(cm, before);
    let mut feature = cf(100.0, &[], "");
    feature
        .metadata
        .insert(IIMN_LINKED_GROUPS.into(), MetaValue::from("1"));
    let mut cm = ConsensusMap::from_features(vec![feature]);
    assert!(matches!(
        annotate_consensus_map(&mut cm),
        Err(Error::InvalidValue(_))
    ));

    let mut cm = ConsensusMap::from_features(vec![cf(100.0, &["1"], ""), cf(110.0, &["1"], "")]);
    annotate_consensus_map(&mut cm)?;
    cm.features[0]
        .metadata
        .insert(IIMN_ADDUCT_PARTNERS.into(), MetaValue::from("9"));
    assert!(matches!(
        write_supplementary_pair_table(&cm, Vec::new()),
        Err(Error::InvalidValue(_))
    ));
    Ok(())
}
