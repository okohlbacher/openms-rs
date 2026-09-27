// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//
// PeptideIdentification::buildUIDsFromAllPepIDs / buildUIDFromPepID and
// buildUSI(mapping) over the class-test consensus map MQEvidence_2.consensusXML
// (copied from src/tests/class_tests/openms/data at bc9cc12).
//
// Expected values are EXECUTED C++ output (tier 1): case `uids` of
// ../oracle/B8b-peptide-id/driver.cpp, built against the Linux x86_64 Release
// install openms4-release-bc9cc12-c19e494-174b576 and run twice on ibminode06
// (byte-identical); see tests/peptide_identification.rs for the fixture layout.
// The source class test asserts the first three multimap entries and three
// buildUIDFromPepID values (tier 3); the oracle covers all 109 entries.
#![cfg(feature = "consensusxml")]

use openms::format::consensusxml;
use openms::identification::{IdentifierMSRunMapper, PeptideIdentification};

const ORACLE: &str = include_str!("data/peptide_identification_oracle.tsv");

fn rows() -> Vec<Vec<&'static str>> {
    ORACLE
        .lines()
        .filter_map(|line| line.strip_prefix("uids\t"))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/data/MQEvidence_2.consensusXML"
);

/// The class-test map is loaded with a raised work ceiling. At the pinned
/// revision the reader charges every modified peptide sequence a conservative
/// registry-lookup allowance (`AASequence::parse_with_budget`, about 3.8
/// million units against the bundled registry), so this 186 KB file's 32
/// modified hits exceed the default `max_work` of 50 million, which the C++
/// reader has no counterpart for. That reader limit is not what this file tests.
fn load() -> openms::kernel::ConsensusMap {
    let options = consensusxml::ReadOptions {
        max_work: 1_000_000_000,
        ..consensusxml::ReadOptions::default()
    };
    consensusxml::load_with_options(FIXTURE, &options).unwrap()
}

#[test]
fn every_uid_and_location_matches_the_release_multimap_in_order() {
    let map = load();
    let rows = rows();
    let count = rows.iter().find(|row| row[0] == "count").unwrap();
    assert_eq!(count[2], map.features.len().to_string());
    assert_eq!(
        count[3],
        map.unassigned_peptide_identifications.len().to_string()
    );

    let uids = PeptideIdentification::build_uids_from_all(&map).unwrap();
    // Flatten in multimap order: ascending key, then insertion order.
    let flat: Vec<(String, String, String)> = uids
        .iter()
        .flat_map(|(uid, locations)| {
            locations.iter().map(move |(feature, position)| {
                (
                    uid.clone(),
                    // The source's Size(-1) marks an unassigned identification.
                    feature.map_or_else(|| u64::MAX.to_string(), |index| index.to_string()),
                    position.to_string(),
                )
            })
        })
        .collect();
    let expected: Vec<(String, String, String)> = rows
        .iter()
        .filter(|row| row[0] == "uid")
        .map(|row| (row[1].to_owned(), row[2].to_owned(), row[3].to_owned()))
        .collect();
    assert_eq!(count[1], expected.len().to_string());
    assert_eq!(expected.len(), 109);
    assert_eq!(flat, expected);

    // Source literals (PeptideIdentification_test.cpp): the first entry is the
    // unassigned identification 4, the next two are features 11 and 4.
    let first: Vec<_> = uids.iter().take(3).collect();
    assert!(first[0].0.ends_with("ES-0014b_2.mzML|spectrum=112"));
    assert_eq!(*first[0].1, vec![(None, 4)]);
    assert_eq!(*first[1].1, vec![(Some(11), 0)]);
    assert_eq!(*first[2].1, vec![(Some(4), 0)]);
}

#[test]
fn single_uids_and_mapped_usis_match_the_release_build() {
    let map = load();
    let rows = rows();
    let mapping = IdentifierMSRunMapper::from_runs(&map.protein_identifications).unwrap();
    let feature0: Vec<_> = rows.iter().filter(|row| row[0] == "cf0").collect();
    assert_eq!(
        feature0.len(),
        map.features[0].base.peptide_identifications.len()
    );
    for row in feature0 {
        let index: usize = row[1].parse().unwrap();
        let id = &map.features[0].base.peptide_identifications[index];
        assert_eq!(id.build_uid(&mapping).unwrap(), row[2]);
    }
    for (index, id) in map.unassigned_peptide_identifications.iter().enumerate() {
        let label = format!("usi_unassigned_{index}");
        let row = rows.iter().find(|row| row[0] == label).unwrap();
        let usi = id
            .build_usi_with_mapping(&mapping, "PXD000001", true)
            .unwrap();
        assert_eq!(
            row[1..],
            [
                format!("valid={}", u8::from(usi.is_valid())),
                format!("collection={}", usi.collection),
                format!("ms_run={}", usi.ms_run),
                format!("type={}", usi.index_type),
                format!("index={}", usi.index),
                format!("interpretation={}", usi.interpretation),
                format!("text={usi}"),
            ],
            "{label}"
        );
    }
    assert_eq!(map.unassigned_peptide_identifications.len(), 34);
}
