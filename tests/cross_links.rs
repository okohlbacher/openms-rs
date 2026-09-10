// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Source expectations: pinned CrossLinksDB_test.cpp and OBODataProvider.cpp.

use openms::chemistry::{
    CrossLinksDB, ModificationRecord, ModificationsDB, OboReadOptions, ResidueModification,
    TermSpecificity,
};
use std::collections::BTreeSet;
use std::sync::Arc;

const SYNTHETIC: &str = r#"format-version: 1.2

[Term]
id: XLMOD:95000
name: Two
synonym: "Alias Two" EXACT []
property_value: reactionSites: "2" xsd:nonNegativeInteger
property_value: monoisotopicMass: "12.5" xsd:double
property_value: specificities: "(K,Protein N-term)&(D,Protein C-term)" xsd:string

[Term]
id: XLMOD:95001
name: Mono
property_value: reactionSites: "1" xsd:nonNegativeInteger
property_value: monoisotopicMass: "7" xsd:double
property_value: specificities: "(K)" xsd:string

[Term]
id: XLMOD:95002
name: NoFlag
property_value: monoisotopicMass: "7" xsd:double
property_value: specificities: "(C)" xsd:string

[Term]
id: XLMOD:95003
name: ZeroTerminal
property_value: reactionSites: "2" xsd:nonNegativeInteger
property_value: monoisotopicMass: "0" xsd:double
property_value: specificities: "(Protein N-term)" xsd:string

[Term]
id: XLMOD:95004
name: Ambiguous
property_value: reactionSites: "2" xsd:nonNegativeInteger
property_value: monoisotopicMass: "6" xsd:double
property_value: specificities: "(B,J,Z,X,R)" xsd:string
"#;

#[test]
fn source_cross_link_names_specificities_and_accessions() {
    let global = CrossLinksDB::global();
    assert!(std::ptr::eq(global, CrossLinksDB::global()));
    let db = global.database();
    assert!(db.len() > 10);
    let dss = db
        .get_modification("DSS", Some('K'), Some(TermSpecificity::Anywhere))
        .unwrap();
    assert_eq!(dss.full_id(), "DSS (K)");
    assert_eq!(dss.obo_accession(), Some("XLMOD:02001"));
    assert_eq!(dss.record_id(), None);
    assert_eq!(dss.accession(), "XLMOD:02001");
    assert_eq!(
        db.get_modification("DSS-d0", Some('K'), Some(TermSpecificity::Anywhere))
            .unwrap(),
        dss
    );
    assert_eq!(
        db.get_modification("XLMOD:02001", Some('K'), Some(TermSpecificity::Anywhere))
            .unwrap(),
        dss
    );
    assert_eq!(
        db.get_modification("DSS", None, Some(TermSpecificity::NTerm))
            .unwrap()
            .full_id(),
        "DSS (N-term)"
    );
    assert_eq!(
        db.get_modification("EDC", None, Some(TermSpecificity::CTerm))
            .unwrap()
            .full_id(),
        "EDC (C-term)"
    );
    let edc: BTreeSet<_> = db
        .find("EDC", None, None)
        .iter()
        .map(|m| (m.origin(), m.term_specificity()))
        .collect();
    assert_eq!(
        edc,
        BTreeSet::from([
            (Some('D'), TermSpecificity::Anywhere),
            (Some('E'), TermSpecificity::Anywhere),
            (Some('K'), TermSpecificity::Anywhere),
            (Some('S'), TermSpecificity::Anywhere),
            (Some('T'), TermSpecificity::Anywhere),
            (Some('Y'), TermSpecificity::Anywhere),
            (None, TermSpecificity::NTerm),
            (None, TermSpecificity::CTerm),
        ])
    );
    assert!(
        db.find("EDC", Some('R'), Some(TermSpecificity::Anywhere))
            .is_empty()
    );
    assert!(db.entries().iter().any(|m| m.full_id() == "EDC (T)"));
}

#[test]
fn source_mass_queries_include_isobaric_and_negative_cross_links() {
    let db = CrossLinksDB::global().database();
    let names: BTreeSet<_> = db
        .search_by_mass(138.06807961, 0.00001, Some('K'), None)
        .unwrap()
        .iter()
        .map(|m| m.full_id().to_owned())
        .collect();
    assert!(names.contains("DSS (K)"));
    assert!(names.contains("BS3 (K)"));
    let terminal = db
        .search_by_mass(138.068, 0.01, None, Some(TermSpecificity::NTerm))
        .unwrap();
    assert_eq!(terminal.len(), 2);
    assert_eq!(
        terminal
            .iter()
            .map(|m| m.full_id())
            .collect::<BTreeSet<_>>(),
        BTreeSet::from(["BS3 (N-term)", "DSS (N-term)"])
    );
    assert!(
        db.search_by_mass(800000000.0, 0.1, Some('S'), None)
            .unwrap()
            .is_empty()
    );
    let edc = db.get_modification("EDC (E)", None, None).unwrap();
    assert_eq!(edc.diff_mono_mass(), -18.01056027);
    assert!(
        db.search_by_mass(
            -18.01056027,
            0.0,
            Some('E'),
            Some(TermSpecificity::Anywhere)
        )
        .unwrap()
        .iter()
        .any(|m| m.name() == "EDC")
    );
    let dss = db.get_modification("DSS (K)", None, None).unwrap();
    assert!(
        dss.diff_formula().is_empty(),
        "the source XLMOD provider does not infer a formula from one mass"
    );
    assert!(db.search_by_mass(f64::NAN, 0.01, None, None).is_err());
}

#[test]
fn source_search_list_uses_ontology_accessions_and_is_sorted() {
    let global = CrossLinksDB::global();
    let names = global.all_search_modifications();
    assert_eq!(names, global.all_search_modifications());
    assert_eq!(names.len(), global.database().len());
    assert!(names.windows(2).all(|pair| pair[0] <= pair[1]));
    for required in ["EDC (S)", "EDC (E)", "DSS (K)", "BS3 (N-term)"] {
        assert!(names.iter().any(|name| name == required));
    }
    assert!(!names.iter().any(|name| name == "DSS"));
}

#[test]
fn cross_link_loader_forces_source_filter_and_preserves_terminal_mapping() {
    let options = OboReadOptions::default();
    assert!(!options.cross_links_only);
    let links = CrossLinksDB::from_obo(SYNTHETIC.as_bytes(), &options).unwrap();
    assert!(
        !options.cross_links_only,
        "the caller's options are unchanged"
    );
    let db = links.database();
    assert_eq!(db.len(), 6); // Two:4; NoFlag:1; Ambiguous:R only
    assert!(db.find("Mono", None, None).is_empty());
    assert!(db.find("ZeroTerminal", None, None).is_empty());
    assert_eq!(
        db.find("NoFlag", None, None).len(),
        1,
        "source rejects reactionSites1, rather than requiring2"
    );
    assert_eq!(
        db.find("Ambiguous", None, None)[0].full_id(),
        "Ambiguous (R)"
    );
    assert_eq!(db.find("Alias Two", None, None).len(), 4);
    assert_eq!(
        db.get_modification("Two (N-term)", None, None)
            .unwrap()
            .term_specificity(),
        TermSpecificity::NTerm
    );
    assert_eq!(
        db.get_modification("Two (C-term)", None, None)
            .unwrap()
            .term_specificity(),
        TermSpecificity::CTerm
    );
    let general = ModificationsDB::from_obo(SYNTHETIC.as_bytes(), &options).unwrap();
    assert!(general.find("Two", None, None).is_empty());
    assert_eq!(general.find("Mono", None, None).len(), 1);
}

#[test]
fn owned_additions_and_handles_do_not_mutate_the_global_cross_link_database() {
    let global = CrossLinksDB::global();
    let initial_names = global.all_search_modifications();
    let owned_handle = {
        let mut owned =
            CrossLinksDB::from_obo(SYNTHETIC.as_bytes(), &OboReadOptions::default()).unwrap();
        let first = owned
            .database()
            .get_modification_handle("Two (K)", None, None)
            .unwrap();
        let custom = ResidueModification::from_record(ModificationRecord {
            full_id: "Custom cross-link (K)".into(),
            name: "Custom cross-link".into(),
            origin: Some('K'),
            diff_mono_mass: 25.0,
            ..Default::default()
        })
        .unwrap();
        owned.database_mut().extend_records(vec![custom]).unwrap();
        assert_eq!(first.diff_mono_mass(), 12.5);
        assert!(
            owned
                .database()
                .get_modification("Custom cross-link (K)", None, None)
                .is_ok()
        );
        assert!(
            !owned
                .all_search_modifications()
                .iter()
                .any(|name| name == "Custom cross-link (K)")
        );
        let second = owned
            .database()
            .get_modification_handle("Two (K)", None, None)
            .unwrap();
        assert!(Arc::ptr_eq(&first, &second));
        first
    };
    assert_eq!(owned_handle.full_id(), "Two (K)");
    assert_eq!(owned_handle.diff_mono_mass(), 12.5);
    assert_eq!(global.all_search_modifications(), initial_names);
    assert!(
        global
            .database()
            .find("Custom cross-link", None, None)
            .is_empty()
    );
}

#[test]
fn reader_limits_are_preserved_by_the_cross_link_wrapper() {
    for options in [
        OboReadOptions {
            max_input_bytes: 20,
            ..Default::default()
        },
        OboReadOptions {
            max_line_bytes: 10,
            ..Default::default()
        },
        OboReadOptions {
            max_terms: 1,
            ..Default::default()
        },
        OboReadOptions {
            max_records: 1,
            ..Default::default()
        },
        OboReadOptions {
            max_aliases: 1,
            ..Default::default()
        },
    ] {
        assert!(CrossLinksDB::from_obo(SYNTHETIC.as_bytes(), &options).is_err());
    }
    let malformed = SYNTHETIC.replace("\"12.5\"", "\"NaN\"");
    assert!(CrossLinksDB::from_obo(malformed.as_bytes(), &OboReadOptions::default()).is_err());
    assert!(
        CrossLinksDB::from_obo(&b""[..], &OboReadOptions::default())
            .unwrap()
            .database()
            .is_empty()
    );
}
