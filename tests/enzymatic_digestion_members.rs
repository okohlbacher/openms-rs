// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `EnzymaticDigestion.h` (the `Specificity` enum, its name table, the enzyme
//! name constants, `filterByMissedCleavages`) and `DigestionEnzymeDB.h`
//! (`hasEnzyme(const DigestionEnzymeType*)`, the lowercase name index).
//!
//! Expected values: `ProteaseDB_test.cpp` (hasEnzyme by pointer) at core
//! `bc9cc12` (tier 2), and the executed Release oracle
//! `../oracle/A5a-chem` (ibminode06, Linux x86_64 Release install
//! bc9cc12/c19e494/174b576; tier 1) for the name lookups, the specificity
//! names and the missed-cleavage filter.

use openms::chemistry::digestion::{DigestionSpecificity, ProteaseDB, ProteaseDigestion};
use openms::chemistry::enzymatic_digestion::{
    NAMES_OF_SPECIFICITY, NO_CLEAVAGE, SIZE_OF_SPECIFICITY, Specificity, UNSPECIFIC_CLEAVAGE,
};
use openms::chemistry::rnase::RNaseDB;

#[test]
fn specificity_values_and_names() {
    assert_eq!(SIZE_OF_SPECIFICITY, NAMES_OF_SPECIFICITY.len());
    let expected = [
        ("none", Specificity::None, 0),
        ("semi", Specificity::Semi, 1),
        ("full", Specificity::Full, 2),
        ("unknown", Specificity::Unknown, 3),
        ("no-cterm", Specificity::NoCTerm, 8),
        ("no-nterm", Specificity::NoNTerm, 9),
        ("Full", Specificity::Unknown, 3),
        ("", Specificity::Unknown, 3),
    ];
    for (name, spec, value) in expected {
        assert_eq!(Specificity::by_name(name), spec, "{name:?}");
        assert_eq!(Specificity::by_name(name).value(), value, "{name:?}");
    }
    assert_eq!(Specificity::default(), Specificity::Full);
    assert_eq!(Specificity::NoCTerm.to_string(), "no-cterm");
    assert_eq!(NO_CLEAVAGE, "no cleavage");
    assert_eq!(UNSPECIFIC_CLEAVAGE, "unspecific cleavage");
    assert!(ProteaseDB::global().has_enzyme(NO_CLEAVAGE));
    assert!(ProteaseDB::global().has_enzyme(UNSPECIFIC_CLEAVAGE));
}

#[test]
fn source_specificity_round_trip_and_refusals() {
    let mut digestion = ProteaseDigestion::default();
    assert_eq!(digestion.source_specificity(), Specificity::Full);
    digestion.set_source_specificity(Specificity::Semi).unwrap();
    assert_eq!(digestion.specificity, DigestionSpecificity::Semi);
    for unsupported in [
        Specificity::Unknown,
        Specificity::NoCTerm,
        Specificity::NoNTerm,
    ] {
        assert!(digestion.set_source_specificity(unsupported).is_err());
        assert_eq!(digestion.specificity, DigestionSpecificity::Semi);
    }
    for supported in [
        DigestionSpecificity::None,
        DigestionSpecificity::Semi,
        DigestionSpecificity::Full,
    ] {
        assert_eq!(
            DigestionSpecificity::try_from(Specificity::from(supported)).unwrap(),
            supported
        );
    }
}

#[test]
fn filter_by_missed_cleavages_passes_the_internal_site_count() {
    let trypsin = ProteaseDigestion::default();
    for (sequence, count) in [("ACKDERFGK", 2), ("ACKDEFGK", 1), ("ACKPDEFGK", 0)] {
        let mut seen = -1;
        let filtered = trypsin
            .filter_by_missed_cleavages(sequence, |mc| {
                seen = mc;
                mc > 1
            })
            .unwrap();
        assert_eq!(seen, count, "{sequence}");
        assert_eq!(filtered, count > 1, "{sequence}");
    }
}

#[test]
fn registry_lookup_by_lowercase_name_and_by_record() {
    let db = ProteaseDB::global();
    for (name, present) in [
        ("Trypsin", true),
        ("trypsin", true),
        ("TRYPSIN", false),
        ("Try", false),
        ("lys-c", true),
        ("arg-c", true),
        ("Clostripain", true),
        ("clostripain", false),
    ] {
        assert_eq!(db.has_enzyme(name), present, "{name}");
    }
    assert_eq!(db.get_enzyme("trypsin").unwrap().name(), "Trypsin");

    let trypsin = db.get_enzyme("Trypsin").unwrap();
    assert!(db.has_enzyme_record(trypsin));
    let copy = *trypsin;
    assert!(!db.has_enzyme_record(&copy));

    let rnases = RNaseDB::global();
    let t1 = rnases.get_enzyme("RNase_T1").unwrap();
    assert!(rnases.has_enzyme_record(&t1));
    let detached = (*t1).clone();
    assert!(!rnases.has_enzyme_record(&detached));
}
