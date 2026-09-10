// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::{
    EmpiricalFormula, ModificationRecord, ModificationsDB, NeutralLoss, OboReadOptions,
    ResidueModification, TermSpecificity,
};
use std::collections::BTreeSet;
use std::io::{self, BufReader, Cursor, Read};
use std::sync::Arc;

fn record() -> ResidueModification {
    ResidueModification::from_record(ModificationRecord {
        name: "Laboratory delta".into(),
        full_name: "Caller supplied modification".into(),
        origin: Some('M'),
        diff_mono_mass: 12.5,
        diff_average_mass: 12.6,
        obo_accession: Some("LAB:1".into()),
        synonyms: BTreeSet::from(["My mass".into()]),
        ..Default::default()
    })
    .unwrap()
}
const OBO: &str = "[Term]\nid: MOD:900\nname: Native record\nsynonym: \"Second name\" EXACT []\nproperty_value: Origin: \"M,K,M,B,J,Z\" xsd:string\nproperty_value: DiffMono: \"12.5\" xsd:float\n";

#[test]
fn owned_records_and_handles_share_identity_and_release_without_leaks() {
    let db = ModificationsDB::from_records(vec![record()]).unwrap();
    let first = db
        .get_modification_handle("LAB:1", Some('M'), None)
        .unwrap();
    assert_eq!(first.record_id(), None);
    assert_eq!(first.accession(), "LAB:1");
    assert_eq!(first.unimod_accession(), None);
    assert_eq!(first.obo_accession(), Some("LAB:1"));
    let second = db.find_handles("My mass", None, None).pop().unwrap();
    assert!(Arc::ptr_eq(&first, &second));
    assert!(std::ptr::eq(
        first.as_ref(),
        db.get_modification("My mass", None, None).unwrap()
    ));
    assert!(Arc::ptr_eq(
        &first,
        &db.search_by_mass_handles(12.5, 0.0, None, None).unwrap()[0]
    ));
    assert!(
        db.best_by_mass_handle(12.5, 0.0, None, None)
            .unwrap()
            .is_none()
    );
    assert!(Arc::ptr_eq(
        &first,
        &db.best_by_mass_handle(12.5, 0.1, None, None)
            .unwrap()
            .unwrap()
    ));
    let weak = Arc::downgrade(&first);
    drop(db);
    drop(second);
    assert_eq!(first.full_id(), "Laboratory delta (M)");
    drop(first);
    assert!(weak.upgrade().is_none());
}

#[test]
fn descriptor_validation_and_absolute_composition_are_independent() {
    let formula = EmpiricalFormula::parse("C4H9NO3").unwrap();
    let m = ResidueModification::from_record(ModificationRecord {
        name: "Replacement".into(),
        origin: Some('X'),
        term_specificity: TermSpecificity::NTerm,
        absolute_formula: Some(formula.clone()),
        record_id: Some(42),
        neutral_losses: vec![
            NeutralLoss::new(EmpiricalFormula::parse("H2O").unwrap(), 18.0, 18.1).unwrap(),
        ],
        ..Default::default()
    })
    .unwrap();
    assert_eq!(m.full_id(), "Replacement (N-term)");
    assert_eq!(m.origin(), None);
    assert_eq!(m.absolute_formula(), Some(&formula));
    assert!(m.diff_formula().is_empty());
    assert_eq!(m.record_id(), Some(42));
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            ResidueModification::from_record(ModificationRecord {
                name: "Bad".into(),
                mono_mass: bad,
                ..Default::default()
            })
            .is_err()
        );
    }
    for name in ["", "bad\nname", "bad\0name"] {
        assert!(
            ResidueModification::from_record(ModificationRecord {
                name: name.into(),
                ..Default::default()
            })
            .is_err()
        );
    }
    assert!(NeutralLoss::new(formula, f64::NAN, 1.0).is_err());
}

#[test]
fn obo_stream_limits_cover_expansion_existing_registry_and_atomic_extension() {
    let options = OboReadOptions::default();
    let db = ModificationsDB::from_obo(Cursor::new(OBO), &options).unwrap();
    assert_eq!(db.len(), 2);
    assert_eq!(
        db.entries().iter().map(|m| m.origin()).collect::<Vec<_>>(),
        [Some('K'), Some('M')]
    );
    assert_eq!(db.find("Second name", None, None).len(), 2);
    for limited in [
        OboReadOptions {
            max_input_bytes: OBO.len() - 1,
            ..options.clone()
        },
        OboReadOptions {
            max_line_bytes: 10,
            ..options.clone()
        },
        OboReadOptions {
            max_terms: 0,
            ..options.clone()
        },
        OboReadOptions {
            max_records: 1,
            ..options.clone()
        },
        OboReadOptions {
            max_aliases: 1,
            ..options.clone()
        },
        OboReadOptions {
            max_registry_bytes: 1,
            ..options.clone()
        },
    ] {
        let mut target = ModificationsDB::from_records(vec![record()]).unwrap();
        let before = target.entries()[0].clone();
        assert!(target.extend_obo(Cursor::new(OBO), &limited).is_err());
        assert_eq!(target.len(), 1);
        assert!(Arc::ptr_eq(&before, &target.entries()[0]));
        assert!(target.find("MOD:900", None, None).is_empty());
    }
    let mut target = db.clone();
    assert!(
        target
            .extend_obo(
                Cursor::new(""),
                &OboReadOptions {
                    max_records: 1,
                    ..options
                }
            )
            .is_err()
    );
    assert_eq!(target.len(), 2);
}

#[test]
fn malformed_and_interrupted_streams_never_publish_partial_records() {
    struct Broken {
        data: Cursor<Vec<u8>>,
    }
    impl Read for Broken {
        fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
            if self.data.position() == self.data.get_ref().len() as u64 {
                Err(io::Error::other("injected"))
            } else {
                self.data.read(out)
            }
        }
    }
    let mut db = ModificationsDB::from_records(vec![record()]).unwrap();
    assert!(
        db.extend_obo(
            BufReader::with_capacity(
                7,
                Broken {
                    data: Cursor::new(OBO.as_bytes().to_vec())
                }
            ),
            &Default::default()
        )
        .is_err()
    );
    assert_eq!(db.len(), 1);
    for tail in [
        "property_value: DiffMono: \"NaN\" xsd:float",
        "property_value: DiffFormula: \"Invalid\" xsd:string",
        "property_value: Something: unquoted",
        "synonym: unquoted",
        "property_value: TermSpec: \"sideways\" xsd:string",
    ] {
        let text = format!("{OBO}{tail}\n");
        assert!(
            db.extend_obo(text.as_bytes(), &Default::default()).is_err(),
            "{tail}"
        );
        assert_eq!(db.len(), 1);
    }
    assert!(
        db.extend_obo(&b"[Term]\nid: MOD:\xff\n"[..], &Default::default())
            .is_err()
    );
}

#[test]
fn full_table_schema_stays_unchanged_while_global_appends_monolinks() {
    let table = include_str!("../resources/modifications/openms-rust-modifications.tsv");
    let db = ModificationsDB::from_tsv(table).unwrap();
    assert_eq!(db.len(), 3035);
    assert_eq!(ModificationsDB::global().len(), 3127);
    assert!(
        ModificationsDB::global().entries()[3035..]
            .iter()
            .all(|m| m.record_id().is_none() && m.obo_accession().is_some())
    );
}

#[test]
fn record_order_distinguishes_every_stored_identity_and_chemical_field() {
    use std::cmp::Ordering;
    let baseline = ModificationRecord {
        name: "Shared ID".into(),
        full_name: "Shared description".into(),
        full_id: "Shared ID (M)".into(),
        origin: Some('M'),
        record_id: Some(1),
        obo_accession: Some("LAB:1".into()),
        ..Default::default()
    };
    let edits: &[fn(&mut ModificationRecord)] = &[
        |r| r.record_id = Some(2),
        |r| r.obo_accession = Some("LAB:2".into()),
        |r| {
            r.synonyms.insert("Another alias".into());
        },
        |r| r.absolute_formula = Some(EmpiricalFormula::default()),
        |r| r.absolute_formula = Some(EmpiricalFormula::parse("C").unwrap()),
        |r| r.absolute_formula = Some(EmpiricalFormula::parse("C+").unwrap()),
        |r| r.name = "Another short ID".into(),
        |r| r.full_name = "Another description".into(),
        |r| r.full_id = "Another full ID (M)".into(),
        |r| r.origin = Some('K'),
        |r| r.term_specificity = TermSpecificity::NTerm,
        |r| r.diff_formula = EmpiricalFormula::parse("C").unwrap(),
        |r| r.diff_formula = EmpiricalFormula::parse("(13)C").unwrap(),
        |r| r.diff_formula = EmpiricalFormula::parse("+").unwrap(),
        |r| r.diff_mono_mass = 1.0,
        |r| r.diff_average_mass = 1.0,
        |r| r.mono_mass = 1.0,
        |r| r.average_mass = 1.0,
        |r| r.hidden = true,
        |r| r.classification = "Custom".into(),
        |r| {
            r.neutral_losses
                .push(NeutralLoss::new(EmpiricalFormula::default(), 0., 0.).unwrap())
        },
        |r| {
            r.neutral_losses
                .push(NeutralLoss::new(EmpiricalFormula::parse("H2O").unwrap(), 0., 0.).unwrap())
        },
        |r| {
            r.neutral_losses
                .push(NeutralLoss::new(EmpiricalFormula::parse("+").unwrap(), 0., 0.).unwrap())
        },
        |r| {
            r.neutral_losses
                .push(NeutralLoss::new(EmpiricalFormula::default(), 1., 0.).unwrap())
        },
        |r| {
            r.neutral_losses
                .push(NeutralLoss::new(EmpiricalFormula::default(), 0., 1.).unwrap())
        },
    ];
    let original = ResidueModification::from_record(baseline.clone()).unwrap();
    let mut values = vec![original.clone()];
    for edit in edits {
        let mut changed = baseline.clone();
        edit(&mut changed);
        let changed = ResidueModification::from_record(changed).unwrap();
        assert_ne!(original, changed);
        assert_ne!(original.cmp(&changed), Ordering::Equal);
        values.push(changed);
    }
    for left in &values {
        for right in &values {
            assert_eq!(left.cmp(right) == Ordering::Equal, left == right);
            assert_eq!(left.partial_cmp(right), Some(left.cmp(right)));
            assert_eq!(left.cmp(right), right.cmp(left).reverse());
        }
    }
    let all: BTreeSet<_> = values.into_iter().collect();
    assert_eq!(all.len(), edits.len() + 1);
    assert_eq!(original.cmp(&original.clone()), Ordering::Equal);
}

#[test]
fn signed_zero_order_agrees_with_eq_in_all_masses_and_neutral_losses() {
    use std::cmp::Ordering;
    let make = |mass| {
        ResidueModification::from_record(ModificationRecord {
            name: "Zero".into(),
            diff_mono_mass: mass,
            diff_average_mass: mass,
            mono_mass: mass,
            average_mass: mass,
            neutral_losses: vec![
                NeutralLoss::new(EmpiricalFormula::default(), mass, mass).unwrap(),
            ],
            ..Default::default()
        })
        .unwrap()
    };
    let negative = make(-0.0);
    let positive = make(0.0);
    assert_eq!(negative, positive);
    assert_eq!(negative.cmp(&positive), Ordering::Equal);
    assert_eq!(negative.partial_cmp(&positive), Some(Ordering::Equal));
    assert_eq!(
        negative.neutral_losses()[0].cmp(&positive.neutral_losses()[0]),
        Ordering::Equal
    );
    assert_eq!(
        BTreeSet::from([negative.clone(), positive.clone()]).len(),
        1
    );
    assert!(make(-1.0) < negative);
    assert!(make(1.0) > positive);

    let first = NeutralLoss::new(EmpiricalFormula::parse("H2O").unwrap(), 18., 18.).unwrap();
    let second = NeutralLoss::new(EmpiricalFormula::parse("NH3").unwrap(), 17., 17.).unwrap();
    let left = ResidueModification::from_record(ModificationRecord {
        name: "Loss order".into(),
        neutral_losses: vec![first.clone(), second.clone()],
        ..Default::default()
    })
    .unwrap();
    let right = ResidueModification::from_record(ModificationRecord {
        name: "Loss order".into(),
        neutral_losses: vec![second, first],
        ..Default::default()
    })
    .unwrap();
    assert_ne!(left, right);
    assert_ne!(left.cmp(&right), Ordering::Equal);
}
