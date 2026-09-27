// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Class tests for the X! Tandem and OMSSA XML readers.
//!
//! Every expected value is transcribed from the upstream class test named at
//! each case, at the core pin (tier 2), and run against the unmodified
//! upstream fixtures described in `tests/data/search_engine_files_provenance.json`.
//! Cases the upstream tests do not cover are derived from the source's control
//! flow and say so.
#![cfg(feature = "idxml")]

use openms::Error;
use openms::chemistry::{AASequence, ModificationDefinitionsSet};
use openms::data_structures::DateTime;
use openms::format::xtandem_xml_file::XTandemXMLFile;
use openms::metadata::MetaValueData;
use std::io::Cursor;
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn reader() -> XTandemXMLFile {
    XTandemXMLFile {
        run_date: Some("2026-09-27 12:00:00".parse::<DateTime>().unwrap()),
        ..XTandemXMLFile::new()
    }
}

#[test]
fn xtandem_load_matches_the_class_test() {
    // XTandemXMLFile_test.cpp, section load.
    let mods = ModificationDefinitionsSet::from_comma_separated(
        "",
        "Carbamidomethyl (C),Oxidation (M),Carboxymethyl (C)",
    )
    .unwrap();
    let result = reader()
        .load(data("XTandemXMLFile_test.xml"), &mods)
        .unwrap();
    assert_eq!(result.peptide_identifications.len(), 303);
    assert_eq!(result.protein_identification.hits.len(), 497);
    // "should have picked up the default N-terminal modifications":
    assert_eq!(result.modification_definitions.variable_names().len(), 6);
    assert_eq!(result.modification_definitions.fixed_names().len(), 0);

    let result = reader()
        .load(data("XTandemXMLFile_test_2.xml"), &mods)
        .unwrap();
    assert_eq!(result.peptide_identifications.len(), 2);
    assert_eq!(result.protein_identification.hits.len(), 21);
    // "no additional modifications in this case":
    assert_eq!(result.modification_definitions.variable_names().len(), 3);
    assert_eq!(result.modification_definitions.fixed_names().len(), 0);
}

#[test]
fn xtandem_first_domain_of_the_second_fixture() {
    // Derived from XTandemXMLFile.cpp and the fixture. The fixture holds the
    // model groups 39 (20 proteins) and 2 (one protein), every domain with the
    // sequence YEDEINKR. The source's std::map orders the identifications by
    // X! Tandem id, so group 2 comes first. In group 39 (charge 2) the first
    // domain "39.1.1" (start 270, end 277, pre "YRNK", post "TNAE",
    // hyperscore 27.9) lies in protein uid 4003, and the 19 domains that
    // repeat its sequence add evidences to the same hit.
    let mods = ModificationDefinitionsSet::from_comma_separated("", "").unwrap();
    let result = reader()
        .load(data("XTandemXMLFile_test_2.xml"), &mods)
        .unwrap();
    let ids = &result.peptide_identifications;
    assert_eq!(ids.len(), 2);
    assert_eq!(
        ids[0].spectrum_reference(),
        "2013-01-18_093003_proteo227_3148_1_tandem_input_file.mzData scan 2 (charge 2)"
    );
    assert_eq!(
        ids[1].spectrum_reference(),
        "2013-01-18_093003_proteo227_3148_1_tandem_input_file.mzData scan 39 (charge 2)"
    );
    assert_eq!(ids[0].hits.len(), 1);
    assert_eq!(ids[0].hits[0].evidences.len(), 1);
    let id = &ids[1];
    assert_eq!(id.hits.len(), 1);
    let hit = &id.hits[0];
    assert_eq!(hit.sequence, AASequence::parse("YEDEINKR").unwrap());
    assert_eq!(hit.evidences.len(), 20);
    assert_eq!(id.score_type, "XTandem");
    assert!(id.higher_score_better);
    assert_eq!(id.identifier, "XTandem_2026-09-27");
    assert_eq!(hit.charge, 2);
    assert_eq!(hit.score, 27.9);
    let meta = |name: &str| hit.metadata.get(name).map(|value| value.data().clone());
    assert_eq!(meta("nextscore"), Some(MetaValueData::Float(22.6)));
    assert_eq!(meta("delta"), Some(MetaValueData::Float(0.0029)));
    assert_eq!(meta("mass"), Some(MetaValueData::Float(1066.5164)));
    assert_eq!(meta("E-Value"), Some(MetaValueData::Float(0.12)));
    assert_eq!(meta("hyperscore"), Some(MetaValueData::Float(27.9)));
    assert_eq!(meta("y_score"), Some(MetaValueData::Float(12.8)));
    assert_eq!(meta("y_ions"), Some(MetaValueData::Integer(7)));
    assert_eq!(meta("b_score"), Some(MetaValueData::Float(10.0)));
    assert_eq!(meta("b_ions"), Some(MetaValueData::Integer(1)));
    assert_eq!(meta("a_score"), None);
    let first = &hit.evidences[0];
    assert_eq!(
        first.protein_accession,
        "P04264 SWISS-PROT:P04264 Tax_Id=9606 Gene_Symbol=KRT1 Keratin, type II cytoskeletal 1"
    );
    assert_eq!(first.start, Some(269));
    assert_eq!(first.end, Some(276));
    assert_eq!(first.aa_before.code(), 'K');
    assert_eq!(first.aa_after.code(), 'T');
    let protein = &result.protein_identification;
    assert!(!protein.higher_score_better);
    assert_eq!(protein.search_engine, "XTandem");
    assert_eq!(protein.date_time.as_deref(), Some("2026-09-27 12:00:00"));
    let scores: Vec<f64> = protein.hits.iter().map(|hit| hit.score).collect();
    assert!(
        scores.windows(2).all(|pair| pair[0] <= pair[1]),
        "{scores:?}"
    );
}

#[test]
fn xtandem_aa_outside_a_group_is_refused() {
    // Derived from XTandemXMLFile.cpp: the source logs a non-fatal error and
    // then reads the top of an empty stack.
    let text = "<bioml><aa type=\"M\" at=\"1\" modified=\"15.9949\"/></bioml>";
    let mods = ModificationDefinitionsSet::from_comma_separated("", "").unwrap();
    match reader().load_reader(Cursor::new(text), &mods) {
        Err(Error::Parse { line, message }) => {
            assert_eq!(line, 1);
            assert!(message.contains("outside of a 'group'"), "{message}");
        }
        other => panic!("expected a parse error, got {other:?}"),
    }
}

#[test]
fn xtandem_missing_required_attribute_is_fatal() {
    // Derived from XMLHandler::attributeAsString_: fatalError on a missing
    // required attribute.
    let text = "<bioml>\n<group type=\"model\" z=\"2\"><protein uid=\"1\" label=\"P\">\
                <peptide><domain id=\"1.1.1\" pre=\"K\" post=\"R\" start=\"2\" end=\"4\"/>\
                </peptide></protein></group></bioml>";
    let mods = ModificationDefinitionsSet::from_comma_separated("", "").unwrap();
    match reader().load_reader(Cursor::new(text), &mods) {
        Err(Error::Parse { line, message }) => {
            assert_eq!(line, 2);
            assert_eq!(message, "Required attribute 'seq' not present!");
        }
        other => panic!("expected a parse error, got {other:?}"),
    }
}
