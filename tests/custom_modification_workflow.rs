// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Caller-owned registry lifetimes, source formula precedence and interchange.

use openms::chemistry::{
    AASequence, EmpiricalFormula, ModificationRecord, ModificationsDB, ModifiedPeptideGenerator,
    ResidueModification, SequenceModification, TermSpecificity,
};
use std::sync::Arc;

fn formula(text: &str) -> EmpiricalFormula {
    text.parse().unwrap()
}
fn record(
    name: &str,
    origin: char,
    delta: f64,
    diff: &str,
    absolute: Option<&str>,
) -> ResidueModification {
    ResidueModification::from_record(ModificationRecord {
        name: name.into(),
        full_name: format!("Laboratory {name}"),
        origin: Some(origin),
        diff_mono_mass: delta,
        diff_formula: formula(diff),
        absolute_formula: absolute.map(formula),
        ..Default::default()
    })
    .unwrap()
}
fn known(sequence: &AASequence, index: usize) -> &Arc<ResidueModification> {
    match sequence.residue_modification(index).unwrap().unwrap() {
        SequenceModification::Known(record) => record,
        _ => panic!("expected resolved registry chemistry"),
    }
}
fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 1e-10, "{actual} vs {expected}");
}

#[test]
fn caller_registry_resolves_names_and_masses_without_static_storage() {
    let database =
        ModificationsDB::from_records(vec![record("LabO", 'M', 15.994915, "O", None)]).unwrap();
    let handle = database
        .get_modification_handle("LabO", Some('M'), Some(TermSpecificity::Anywhere))
        .unwrap();
    let parsed = AASequence::parse_with_registry("AM(LabO)K", &database).unwrap();
    let numeric = AASequence::parse_with_registry("AM[+15.994915]K", &database).unwrap();
    let mut edited = AASequence::parse("AMK").unwrap();
    edited
        .set_modification_with_registry(1, "LabO", &database)
        .unwrap();
    assert_eq!(parsed, numeric);
    assert_eq!(parsed, edited);
    assert!(Arc::ptr_eq(known(&parsed, 1), &handle));
    let weak = Arc::downgrade(&handle);
    drop(handle);
    drop(database);
    assert!(weak.upgrade().is_some());
    assert_eq!(
        parsed.formula().unwrap(),
        AASequence::parse("AM(Oxidation)K")
            .unwrap()
            .formula()
            .unwrap()
    );
    assert_eq!(parsed.suffix(2).unwrap().as_str(), "MK");
    assert_eq!(parsed.fragment_ions(2).unwrap().len(), 8);
    assert!(AASequence::parse(&parsed.to_string()).is_err());
    drop(parsed);
    drop(numeric);
    drop(edited);
    assert!(
        weak.upgrade().is_none(),
        "last owned peptide releases its custom chemistry"
    );
}

#[test]
fn generated_variants_share_records_and_retain_custom_chemistry() {
    let (variants, weak) = {
        let database =
            ModificationsDB::from_records(vec![record("MassOnly", 'M', 12.5, "", None)]).unwrap();
        let handles =
            ModifiedPeptideGenerator::get_modifications_with_registry(&["MassOnly (M)"], &database)
                .unwrap();
        let weak = Arc::downgrade(&handles[0]);
        let output = ModifiedPeptideGenerator::default()
            .variable_modifications(&handles, &"AMMK".parse().unwrap(), 2, false)
            .unwrap();
        (output, weak)
    };
    assert_eq!(variants.len(), 3);
    assert!(Arc::ptr_eq(known(&variants[0], 2), known(&variants[2], 2)));
    let base: AASequence = "AMMK".parse().unwrap();
    for (sequence, shift) in variants.iter().zip([12.5, 12.5, 25.]) {
        near(
            sequence.mono_mass().unwrap(),
            base.mono_mass().unwrap() + shift,
        );
        assert!(sequence.formula().is_err());
        assert!(sequence.fragment_ions(2).is_ok());
    }
    assert!(weak.upgrade().is_some());
    drop(variants);
    assert!(weak.upgrade().is_none());
}

#[test]
fn absolute_formula_replaces_a_free_residue_and_can_resolve_unknown_chemistry() {
    let records = vec![
        record("Restore", 'X', 1., "", Some("C3H7NO2"))
            .with_absolute_masses(999., 998.)
            .unwrap(),
        record("Replace", 'M', 1., "", Some("C3H7NO2"))
            .with_absolute_masses(999., 998.)
            .unwrap(),
        record("DeltaWins", 'M', 1., "O", Some("C3H7NO2"))
            .with_absolute_masses(999., 998.)
            .unwrap(),
    ];
    let db = ModificationsDB::from_records(records).unwrap();
    for input in ["AX(Restore)", "AM(Replace)"] {
        let actual = AASequence::parse_with_registry(input, &db).unwrap();
        let expected: AASequence = "AA".parse().unwrap();
        assert_eq!(actual.formula().unwrap(), expected.formula().unwrap());
        near(actual.mono_mass().unwrap(), expected.mono_mass().unwrap());
        assert_eq!(
            actual.average_mass().unwrap(),
            expected.average_mass().unwrap()
        );
        for (a, e) in actual
            .fragment_ions(2)
            .unwrap()
            .iter()
            .zip(expected.fragment_ions(2).unwrap())
        {
            near(a.mz, e.mz);
        }
        // An absolute replacement makes the total formula known, but is not a
        // supplied delta formula, which would require a known starting residue.
        assert!(
            actual
                .residue_modification(1)
                .unwrap()
                .unwrap()
                .diff_formula()
                .is_err()
        );
    }
    let actual = AASequence::parse_with_registry("AM(DeltaWins)", &db).unwrap();
    let expected: AASequence = "AM(Oxidation)".parse().unwrap();
    assert_eq!(actual.formula().unwrap(), expected.formula().unwrap());
    assert_eq!(actual.mono_mass().unwrap(), expected.mono_mass().unwrap());
}

#[test]
fn no_change_and_terminal_formula_rules_remain_distinct() {
    let mut terminal = ModificationRecord {
        name: "Terminal".into(),
        full_name: "Laboratory terminal".into(),
        term_specificity: TermSpecificity::NTerm,
        diff_mono_mass: 12.5,
        absolute_formula: Some(formula("C3H7NO2")),
        mono_mass: 999.,
        ..Default::default()
    };
    let db = ModificationsDB::from_records(vec![
        record("NoChange", 'T', 0., "", Some("C4H7NO2"))
            .with_absolute_masses(101., 101.)
            .unwrap(),
        record("UnknownNoChange", 'X', 0., "", Some("C3H7NO2")),
        ResidueModification::from_record(terminal.clone()).unwrap(),
    ])
    .unwrap();
    let actual = AASequence::parse_with_registry("T(NoChange)A", &db).unwrap();
    let expected: AASequence = "TA".parse().unwrap();
    assert_eq!(actual.formula().unwrap(), expected.formula().unwrap());
    assert_eq!(actual.mono_mass().unwrap(), expected.mono_mass().unwrap());
    let unresolved = AASequence::parse_with_registry("X(UnknownNoChange)", &db).unwrap();
    assert!(unresolved.formula().is_err());
    assert!(unresolved.mono_mass().is_err());
    assert!(unresolved.to_unimod_string().is_err());
    let actual = AASequence::parse_with_registry(".(Terminal)A", &db).unwrap();
    near(
        actual.mono_mass().unwrap(),
        AASequence::parse("A").unwrap().mono_mass().unwrap() + 12.5,
    );
    assert!(
        actual.formula().is_err(),
        "terminal formula uses its delta, not absolute replacement"
    );
    terminal.diff_mono_mass = -1e6;
    let bad =
        ModificationsDB::from_records(vec![ResidueModification::from_record(terminal).unwrap()])
            .unwrap();
    let mut unchanged = expected.clone();
    assert!(
        unchanged
            .set_n_terminal_modification_with_registry("Terminal", &bad)
            .is_err()
    );
    assert_eq!(unchanged, expected);
}

#[test]
fn vocabulary_accessions_and_unimod_mass_export_never_invent_ids() {
    let mut r = ModificationRecord {
        name: "MOD:99999".into(),
        full_name: "Laboratory oxygen".into(),
        origin: Some('M'),
        obo_accession: Some("MOD:99999".into()),
        diff_mono_mass: 15.994915,
        diff_formula: formula("O"),
        ..Default::default()
    };
    let db =
        ModificationsDB::from_records(vec![ResidueModification::from_record(r.clone()).unwrap()])
            .unwrap();
    let p = AASequence::parse_with_registry("AM(MOD:99999)", &db).unwrap();
    assert_eq!(
        p.residue_modification(1).unwrap().unwrap().record_id(),
        None
    );
    assert_eq!(p.to_accession_string(), "AM(MOD:99999)");
    assert_eq!(
        AASequence::parse_with_registry(&p.to_accession_string(), &db).unwrap(),
        p
    );
    let mass_text = p.to_unimod_string().unwrap();
    assert!(mass_text.starts_with("AM["));
    assert!(!mass_text.contains("UniMod"));
    near(
        AASequence::parse(&mass_text).unwrap().mono_mass().unwrap(),
        p.mono_mass().unwrap(),
    );
    r.origin = None;
    r.term_specificity = TermSpecificity::NTerm;
    r.diff_formula = formula("");
    r.diff_mono_mass = 12.5;
    let db =
        ModificationsDB::from_records(vec![ResidueModification::from_record(r).unwrap()]).unwrap();
    let p = AASequence::parse_with_registry(".(MOD:99999)A", &db).unwrap();
    let mass_text = p.to_unimod_string().unwrap();
    assert!(mass_text.starts_with(".["));
    near(
        AASequence::parse(&mass_text).unwrap().mono_mass().unwrap(),
        p.mono_mass().unwrap(),
    );
}

#[cfg(feature = "idxml")]
#[test]
fn custom_registry_idxml_roundtrip_checks_exact_chemistry_before_output() {
    use openms::format::idxml::{self, IdXmlDocument};
    use openms::identification::{PeptideHit, PeptideIdentification, ProteinIdentification};
    let db =
        ModificationsDB::from_records(vec![record("LabO", 'M', 15.994915, "O", None)]).unwrap();
    let p = AASequence::parse_with_registry("AM(LabO)K", &db).unwrap();
    let document = IdXmlDocument {
        protein_identifications: vec![ProteinIdentification {
            identifier: "custom".into(),
            date_time: Some("2026-09-10T12:00:00".into()),
            ..Default::default()
        }],
        peptide_identifications: vec![PeptideIdentification {
            identifier: "custom".into(),
            hits: vec![PeptideHit {
                sequence: p,
                ..Default::default()
            }],
            ..Default::default()
        }],
        ..Default::default()
    };
    let mut bytes = Vec::new();
    idxml::write(&mut bytes, &document).unwrap();
    let mut restored =
        idxml::read_with_registry(bytes.as_slice(), &Default::default(), &db).unwrap();
    assert_eq!(restored, idxml::read(bytes.as_slice()).unwrap());
    assert!(
        restored.protein_identifications[0]
            .search_parameters
            .metadata
            .remove("modification_definitions")
            .is_some()
    );
    assert_eq!(restored, document);
    let conflicting =
        ModificationsDB::from_records(vec![record("LabO", 'M', 31.98983, "O2", None)]).unwrap();
    let mut untouched = vec![99];
    assert!(
        idxml::write_with_registry(&mut untouched, &document, &Default::default(), &conflicting)
            .is_err()
    );
    assert_eq!(untouched, [99]);
    drop(db);
    assert!(
        restored.peptide_identifications[0].hits[0]
            .sequence
            .formula()
            .is_ok()
    );
}

#[test]
fn full_id_only_caller_record_has_reconstructable_sequence_text() {
    let record = ResidueModification::from_record(ModificationRecord {
        full_id: "Custom-only-ID".into(),
        origin: Some('M'),
        diff_mono_mass: 12.5,
        ..Default::default()
    })
    .unwrap();
    let db = ModificationsDB::from_records(vec![record]).unwrap();
    let sequence = AASequence::parse_with_registry("AM(Custom-only-ID)K", &db).unwrap();
    assert_eq!(sequence.to_string(), "AM(Custom-only-ID)K");
    assert_eq!(
        AASequence::parse_with_registry(&sequence.to_string(), &db).unwrap(),
        sequence
    );
    assert_eq!(
        AASequence::parse_with_registry(&sequence.to_accession_string(), &db).unwrap(),
        sequence
    );
}
