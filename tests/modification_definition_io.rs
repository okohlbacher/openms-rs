// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::chemistry::{
    AASequence, EmpiricalFormula, ModificationProvenance, ModificationRecord, ModificationsDB,
    NeutralLoss, ResidueModification, TermSpecificity,
};
use openms::format::modification_definitions as io;
use openms::identification::{
    PeptideHit, PeptideIdentification, ProteinIdentification, SearchParameters,
};
use std::sync::Arc;

fn formula(text: &str) -> EmpiricalFormula {
    text.parse().unwrap()
}
fn record(name: &str, origin: char, diff: &str) -> ResidueModification {
    let diff = formula(diff);
    // Source ModificationDefinitionIO_test helper sets only mono, not average.
    ResidueModification::from_record(ModificationRecord {
        name: name.into(),
        origin: Some(origin),
        diff_mono_mass: diff.mono_mass(),
        diff_formula: diff,
        ..Default::default()
    })
    .unwrap()
}
#[test]
fn source_collect_search_only_hit_cv_exclusion_and_all_specificities() {
    let mut db = ModificationsDB::global().clone();
    db.extend_records(vec![
        record("TestIO:Adduct", 'K', "C9H11N2O8P"),
        record("TestIO:VarOnly", 'S', "HPO3"),
        record("TestIO:TwoSites", 'K', "C2H2O"),
        record("TestIO:TwoSites", 'Y', "C2H2O"),
    ])
    .unwrap();
    let run = ProteinIdentification {
        identifier: "run1".into(),
        search_parameters: SearchParameters {
            variable_modifications: vec![
                "Oxidation (M)".into(),
                "TestIO:VarOnly (S)".into(),
                "TestIO:TwoSites".into(),
            ],
            ..Default::default()
        },
        ..Default::default()
    };
    let peptide = PeptideIdentification {
        identifier: "run1".into(),
        hits: vec![PeptideHit {
            sequence: AASequence::parse_with_registry("PEPTK(TestIO:Adduct)M(Oxidation)IDE", &db)
                .unwrap(),
            ..Default::default()
        }],
        ..Default::default()
    };
    let definitions = io::collect(
        &[
            run,
            ProteinIdentification {
                identifier: "run2".into(),
                ..Default::default()
            },
        ],
        &[peptide],
        &db,
    )
    .unwrap();
    assert_eq!(definitions.len(), 1);
    let names: Vec<_> = definitions["run1"].iter().map(|r| r.full_id()).collect();
    assert_eq!(
        names,
        [
            "TestIO:Adduct (K)",
            "TestIO:TwoSites (K)",
            "TestIO:TwoSites (Y)",
            "TestIO:VarOnly (S)"
        ]
    );
}
#[test]
fn source_registration_recovery_count_and_literal_peptide_chemistry() {
    // Source ModificationDefinitionIO_test.cpp registerFrom section: same name,
    // formula, peptide, formula and rounded mono-mass literals; no C++ execution.
    let definition = record("TestIO:FromBlob", 'K', "C9H11N2O8P");
    let text = definition.to_definition_string().unwrap();
    let mut db = ModificationsDB::default();
    let report = io::register_from(&format!("{text};not a record;{text}"), &mut db).unwrap();
    assert_eq!(report.registered, 2);
    assert_eq!(report.diagnostics.len(), 1);
    assert_eq!(report.diagnostics[0].record, 2);
    assert_eq!(db.len(), 1);
    let peptide = AASequence::parse_with_registry("AEADNLDDK(TestIO:FromBlob)K", &db).unwrap();
    assert_eq!(peptide.formula().unwrap(), formula("C54H86N15O28P1"));
    assert!((peptide.mono_mass().unwrap() - 1423.5504442334).abs() < 1e-8);
    assert_eq!(
        io::register_search_parameters(&SearchParameters::default(), &mut db)
            .unwrap()
            .registered,
        0
    );
}
#[test]
fn escaped_fields_explicit_masses_losses_and_nine_field_compatibility() {
    let water = formula("H2O");
    let original = ResidueModification::from_record(ModificationRecord {
        name: "lab|name;back\\slash".into(),
        full_name: "explicit | ; \\ chemical name".into(),
        full_id: "kept unusual ID".into(),
        origin: None,
        term_specificity: TermSpecificity::ProteinCTerm,
        diff_formula: formula("H-1N-1O"),
        diff_mono_mass: 7.123456789012345,
        diff_average_mass: -0.0,
        neutral_losses: vec![
            NeutralLoss::new(water.clone(), water.mono_mass(), water.average_mass()).unwrap(),
        ],
        ..Default::default()
    })
    .unwrap();
    let encoded = original.to_definition_string().unwrap();
    assert!(encoded.contains("lab\\|name\\;back\\\\slash"));
    let decoded = ResidueModification::from_definition_string(&encoded).unwrap();
    assert_eq!(decoded, original);
    assert_eq!(
        decoded.diff_mono_mass().to_bits(),
        original.diff_mono_mass().to_bits()
    );
    assert_eq!(decoded.diff_average_mass().to_bits(), (-0.0f64).to_bits());
    let joined = format!(";;{encoded};;;{encoded};");
    assert_eq!(
        ResidueModification::split_definition_records(&joined).unwrap(),
        [&encoded, &encoded]
    );
    let nine = ResidueModification::from_definition_string("1|Named|||K|none|||").unwrap();
    assert_eq!(nine.full_id(), "Named (K)");
    assert_eq!(nine.diff_mono_mass(), 0.0);
    assert_eq!(nine.provenance(), ModificationProvenance::Defined);
    assert!(ResidueModification::from_definition_string("1|Named|||K|none||0|0||extra").is_err());
}
#[test]
fn deterministic_union_noop_existing_empty_and_duplicate_run_last_wins() {
    let a = Arc::new(record("TestIO:AttA", 'K', "C2H2O"));
    let b = Arc::new(record("TestIO:AttB", 'R', "C2H2O"));
    assert_eq!(
        io::encode(&[b.clone(), a.clone(), a.clone()]).unwrap(),
        io::encode(&[a.clone(), b.clone()]).unwrap()
    );
    let mut parameters = SearchParameters::default();
    io::attach(&mut parameters, &[]).unwrap();
    assert!(!parameters.metadata.contains_key(io::METADATA_KEY));
    parameters
        .metadata
        .insert(io::METADATA_KEY.into(), "".into());
    io::attach(&mut parameters, &[]).unwrap();
    assert_eq!(parameters.metadata[io::METADATA_KEY].as_str().unwrap(), "");
    io::attach(&mut parameters, std::slice::from_ref(&a)).unwrap();
    let once = parameters.clone();
    io::attach(&mut parameters, &[a]).unwrap();
    assert_eq!(parameters, once);
    io::attach(&mut parameters, &[b]).unwrap();
    assert_eq!(
        ResidueModification::split_definition_records(
            parameters.metadata[io::METADATA_KEY].as_str().unwrap()
        )
        .unwrap()
        .len(),
        2
    );
    let first = ProteinIdentification {
        identifier: "run".into(),
        search_parameters: parameters,
        ..Default::default()
    };
    let last = ProteinIdentification {
        identifier: "run".into(),
        search_parameters: once.clone(),
        ..Default::default()
    };
    assert_eq!(
        io::encode_by_run(&[first, last], &Default::default()).unwrap()["run"],
        once.metadata[io::METADATA_KEY].as_str().unwrap()
    );
}
#[test]
fn provenance_is_not_chemistry_and_portable_writer_rejects_unrepresented_fields() {
    let r = record("lab", 'K', "O");
    assert!(io::is_definition(&r));
    let cv = r.clone().with_provenance(ModificationProvenance::Cv);
    assert_eq!(r, cv);
    assert_eq!(r.cmp(&cv), std::cmp::Ordering::Equal);
    assert!(!io::is_definition(&cv));
    assert!(
        ModificationsDB::global()
            .entries()
            .iter()
            .all(|r| r.provenance() == ModificationProvenance::Cv)
    );
    let absolute = r.clone().with_absolute_masses(123.5, 124.0).unwrap();
    assert!(io::encode(&[Arc::new(absolute)]).is_err());
    let odd_loss = ResidueModification::from_record(ModificationRecord {
        name: "odd-loss".into(),
        neutral_losses: vec![NeutralLoss::new(formula("H2O"), 1.0, 2.0).unwrap()],
        ..Default::default()
    })
    .unwrap();
    assert!(io::encode(&[Arc::new(odd_loss)]).is_err());
    let charged_empty = ResidueModification::from_record(ModificationRecord {
        name: "charged-empty".into(),
        diff_formula: formula("+"),
        ..Default::default()
    })
    .unwrap();
    assert!(io::encode(&[Arc::new(charged_empty)]).is_err());
}
#[test]
fn conflict_malformed_strict_and_exhausted_shared_budgets_are_atomic() {
    let first = record("collision", 'K', "O");
    let second = record("collision", 'K', "H2");
    let new = record("new", 'S', "O");
    let mut db = ModificationsDB::from_records(vec![first.clone()]).unwrap();
    let before = db.entries().to_vec();
    assert!(
        io::register_from(
            &format!(
                "{};{}",
                new.to_definition_string().unwrap(),
                second.to_definition_string().unwrap()
            ),
            &mut db
        )
        .is_err()
    );
    assert_eq!(db.entries(), before);
    let mut params = SearchParameters::default();
    params.metadata.insert(
        io::METADATA_KEY.into(),
        format!("{};broken", new.to_definition_string().unwrap()).into(),
    );
    assert!(io::register_search_parameters(&params, &mut db).is_err());
    assert_eq!(db.entries(), before);
    params.metadata.insert(
        io::METADATA_KEY.into(),
        first.to_definition_string().unwrap().into(),
    );
    let mut work = 0;
    let mut bytes = usize::MAX;
    assert!(
        io::register_search_parameters_with_budget(&params, &mut db, &mut work, &mut bytes)
            .is_err()
    );
    assert_eq!(db.entries(), before);
    let clone = params.clone();
    assert!(io::attach_with_budget(&mut params, &[Arc::new(new)], &mut 100, &mut 0).is_err());
    assert_eq!(params, clone);
}
#[test]
fn anonymous_mass_tags_and_defined_named_terminal_hits_are_distinct() {
    let terminal = ResidueModification::from_record(ModificationRecord {
        name: "lab-end".into(),
        origin: None,
        term_specificity: TermSpecificity::NTerm,
        diff_mono_mass: 12.5,
        diff_average_mass: 13.5,
        ..Default::default()
    })
    .unwrap();
    let db = ModificationsDB::from_records(vec![terminal]).unwrap();
    let peptide = PeptideIdentification {
        identifier: "run".into(),
        hits: vec![PeptideHit {
            sequence: AASequence::parse_with_registry("(lab-end)AK[999]", &db).unwrap(),
            ..Default::default()
        }],
        ..Default::default()
    };
    let defs = io::collect(&[], &[peptide], &db).unwrap();
    assert_eq!(defs["run"].len(), 1);
    assert_eq!(defs["run"][0].full_id(), "lab-end (N-term)");
}

#[test]
fn typed_registration_retains_fields_without_claiming_portable_projection() {
    let definition = record("absolute", 'M', "")
        .with_absolute_masses(500.0, 501.0)
        .unwrap();
    let mut db = ModificationsDB::default();
    assert!(!db.has_defined_modification("absolute"));
    let first = db.register_definition(&definition).unwrap();
    assert!(db.has_defined_modification("absolute"));
    let again = db.register_definition(&definition).unwrap();
    assert!(Arc::ptr_eq(&first, &again));
    assert_eq!(first.mono_mass(), 500.0);
    assert!(io::encode(&[first]).is_err());
    let tiny = ResidueModification::from_record(ModificationRecord {
        name: "tiny".into(),
        diff_mono_mass: 1e-20,
        diff_average_mass: 1e20,
        ..Default::default()
    })
    .unwrap();
    assert!(
        tiny.to_definition_string()
            .unwrap()
            .contains("|1e-20|1e+20|")
    );
}

#[test]
#[allow(clippy::excessive_precision)] // Literal source input intentionally retained.
fn literal_residue_modification_codec_records_use_none_for_anywhere() {
    // ResidueModification_test.cpp lines780/711: literal record and independent
    // exact-roundtrip explicit mass, rather than a record built by our encoder.
    let definition =
        ResidueModification::from_definition_string("1|TestDef:Der|||K|none|O|15.994915|15.9994|")
            .unwrap();
    assert_eq!(definition.full_id(), "TestDef:Der (K)");
    assert_eq!(definition.term_specificity(), TermSpecificity::Anywhere);
    assert!(
        definition
            .to_definition_string()
            .unwrap()
            .contains("|K|none|O1|")
    );
    let literal = "1|TestDef:RoundTrip|TestDef:RoundTrip (K)|a\\|b\\;c\\\\d|K|none|C9H11N2O8P1|306.025304840900048|306.16|H3O4P1,H2O1";
    let decoded = ResidueModification::from_definition_string(literal).unwrap();
    assert_eq!(decoded.full_name(), "a|b;c\\d");
    assert_eq!(
        decoded.diff_mono_mass().to_bits(),
        306.025304840900048f64.to_bits()
    );
    assert_eq!(decoded.diff_average_mass().to_bits(), 306.16f64.to_bits());
    assert_eq!(decoded.neutral_losses().len(), 2);
    assert_eq!(decoded.neutral_losses()[0].formula(), &formula("H3PO4"));
    assert!(ResidueModification::from_definition_string("1|bad|||K|Anywhere|O|1|1|").is_err());
}
