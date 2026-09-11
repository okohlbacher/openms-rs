use openms::chemistry::proforma::*;
use openms::chemistry::{
    AASequence, EmpiricalFormula, ModificationRecord, ModificationsDB, PROTON_MASS_U,
    ResidueModification,
};
use std::sync::Arc;

fn empty_db() -> ModificationsDB {
    ModificationsDB::from_records(vec![]).unwrap()
}
fn base(sequence: &str) -> f64 {
    AASequence::parse(sequence).unwrap().mono_mass().unwrap()
}
fn close(left: f64, right: f64, tolerance: f64) {
    assert!(
        (left - right).abs() <= tolerance,
        "{left:.15} != {right:.15}"
    );
}
fn mono(text: &str, db: &mut ModificationsDB) -> f64 {
    Peptidoform::parse(text)
        .unwrap()
        .mono_mass(db)
        .unwrap()
        .value
}
fn mass(value: f64) -> Modification {
    modification(ModificationTag::MassDelta(MassDelta {
        mass: value,
        ..Default::default()
    }))
}
fn modification(tag: ModificationTag) -> Modification {
    Modification {
        alternatives: vec![(tag, None)],
        resolved_mod: None,
    }
}
fn formula(value: &str, charge: Option<i32>) -> Modification {
    modification(ModificationTag::FormulaTag(FormulaTag {
        formula_string: value.into(),
        charge,
    }))
}
fn named(value: &str) -> Modification {
    modification(ModificationTag::NamedMod(NamedMod {
        name: value.into(),
        cv_hint: None,
    }))
}
fn element(aa: char, modifications: Vec<Modification>) -> SequenceElement {
    SequenceElement {
        amino_acid: aa,
        modifications,
    }
}
fn chain(elements: Vec<SequenceElement>) -> Peptidoform {
    Peptidoform {
        sequence: elements.into_iter().map(SequenceSection::Element).collect(),
        ..Default::default()
    }
}

#[test]
fn all_source_scalar_mass_examples_and_reference_literals() {
    let mut db = ModificationsDB::global().clone();
    // Source ProFormaParser_test.cpp2582ff comparisons use AASequence, with
    // separate rounded literals for DFPIANGER and PEPTIDE. Registry stored
    // deltas are rounded independently of AASequence formula masses; 2e-6
    // tolerance preserves that inherited source-data distinction.
    for (proforma, peptide) in [
        ("PEPTIDE", "PEPTIDE"),
        ("PEM[UNIMOD:35]TIDE", "PEM(Oxidation)TIDE"),
        ("[UNIMOD:1]-PEPTIDE-[UNIMOD:2]", "(Acetyl)PEPTIDE(Amidated)"),
        ("DFPIANGER", "DFPIANGER"),
        ("PEPTS[UNIMOD:21]IDE", "PEPTS(Phospho)IDE"),
        ("PEPTC[UNIMOD:4]IDE", "PEPTC(Carbamidomethyl)IDE"),
        (
            "[UNIMOD:1]-PEM[UNIMOD:35]PTIDES[UNIMOD:21]K",
            "(Acetyl)PEM(Oxidation)PTIDES(Phospho)K",
        ),
    ] {
        close(mono(proforma, &mut db), base(peptide), 2e-6);
    }
    close(mono("DFPIANGER", &mut db), 1017.48796, 1e-4);
    close(mono("PEPTIDE", &mut db), 799.3599, 1e-4);
    close(
        mono("PEM[+15.9949]TIDE", &mut db),
        base("PEMTIDE") + 15.9949,
        2e-5,
    );
    close(
        mono("PEPTS[+79.966331]IDE", &mut db),
        mono("PEPTS[UNIMOD:21]IDE", &mut db),
        1e-8,
    );
    close(
        mono("[Formula:C2H2O]-PEPTIDE", &mut db),
        mono("[UNIMOD:1]-PEPTIDE", &mut db),
        1e-6,
    );
    close(
        mono("[UNIMOD:1]-PEPTIDE-[UNIMOD:2]", &mut db) - base("PEPTIDE"),
        42.010565 - 0.984016,
        1e-6,
    );
}

#[test]
fn source_can_issue_try_and_charge_overloads() {
    let mut db = ModificationsDB::global().clone();
    for text in [
        "PEPTIDE",
        "PEM[UNIMOD:35]TIDE",
        "PEM[+15.9949]TIDE",
        "PEM[Formula:O]TIDE",
        "PEM[UNIMOD:35]PTIDES[UNIMOD:21]K",
    ] {
        let pf = Peptidoform::parse(text).unwrap();
        assert!(pf.can_calculate_mass(&mut db).unwrap().value);
        assert!(
            pf.mass_calculation_issues(&mut db)
                .unwrap()
                .value
                .is_empty()
        );
        let attempt = pf.try_mono_mass(&mut db).unwrap();
        assert!(attempt.issues.is_empty());
        close(
            attempt.value.unwrap(),
            pf.mono_mass(&mut db).unwrap().value,
            1e-10,
        );
        close(
            pf.try_mz(2, &mut db).unwrap().value.unwrap(),
            pf.mz(2, &mut db).unwrap().value,
            1e-10,
        );
    }
    let pf = Peptidoform::parse("PEM[UnknownMod12345]TIDE").unwrap();
    assert!(!pf.can_calculate_mass(&mut db).unwrap().value);
    let attempt = pf.try_mono_mass(&mut db).unwrap();
    assert!(attempt.value.is_none());
    assert_eq!(
        attempt.issues[0].issue_type,
        ConversionIssueType::UnresolvedMod
    );
    assert_eq!(attempt.issues[0].position, Some(2));
    assert!(!attempt.warnings.is_empty());
    assert!(pf.mono_mass(&mut db).is_err());
    let pf = Peptidoform::parse("PEPTIDE").unwrap();
    assert!(pf.mz(0, &mut db).is_err());
    assert!(pf.try_mz(0, &mut db).unwrap().value.is_none());
    let ion = PeptidoformIon::parse("PEPTIDE/2").unwrap();
    close(
        ion.mz(&mut db).unwrap().value,
        (base("PEPTIDE") + 2. * PROTON_MASS_U) / 2.,
        1e-10,
    );
    close(
        ion.try_mz(&mut db).unwrap().value.unwrap(),
        pf.mz(2, &mut db).unwrap().value,
        1e-10,
    );
    let ion = PeptidoformIon::parse("PEPTIDE").unwrap();
    assert!(ion.mz(&mut db).is_err());
    assert_eq!(
        ion.try_mz(&mut db).unwrap().issues[0].description,
        "No charge state specified"
    );
}

#[test]
fn source_unlocalised_labile_global_and_crosslink_examples() {
    let mut db = ModificationsDB::global().clone();
    for (text, sequence, delta) in [
        ("[+79.966331]?PEPTIDE", "PEPTIDE", 79.966331),
        ("{+162.0528}PEPTIDE", "PEPTIDE", 162.0528),
        ("<[+57.0215]@C>PEPTCIDECK", "PEPTCIDECK", 2. * 57.0215),
        (
            "EVTSEKC[-2.0156#XL1]LEMSC[#XL1]EFD",
            "EVTSEKCLEMSCEFD",
            -2.0156,
        ),
        ("PEPTC[-2.01565#XL1]IDEC[#XL1]K", "PEPTCIDECK", -2.01565),
    ] {
        close(mono(text, &mut db), base(sequence) + delta, 1e-4);
    }
    close(
        mono("<[+57.0215]@C>PEPTCIDECK", &mut db),
        mono("PEPTC[+57.0215]IDEC[+57.0215]K", &mut db),
        1e-4,
    );
    for value in ["138.06807961", "138.068"] {
        let ion =
            PeptidoformIon::parse(&format!("PEPTIDEK[+{value}#XL1]//SEQUENCEK[#XL1]")).unwrap();
        assert!(ion.can_calculate_mass(&mut db).unwrap().value);
        let expected = base("PEPTIDEK") + base("SEQUENCEK") + value.parse::<f64>().unwrap();
        let actual = ion.mono_mass(&mut db).unwrap().value;
        if value == "138.06807961" {
            close(actual, expected, 1e-6);
        }
        // The source's shorter 138.068 spelling compares try against strict,
        // not the literal delta: nearest-mass resolution selects DSS 138.06808.
        close(
            ion.try_mono_mass(&mut db).unwrap().value.unwrap(),
            actual,
            1e-10,
        );
    }
    let chim = PeptidoformIon::parse("PEPTIDE/2+EDITPEP/3").unwrap();
    assert!(chim.can_calculate_mass(&mut db).unwrap().value);
    assert!(chim.mono_mass(&mut db).is_err());
    assert!(chim.try_mono_mass(&mut db).unwrap().value.is_none());
    for pf in &chim.chains {
        close(pf.mono_mass(&mut db).unwrap().value, 799.3599, 1e-4);
    }
}

#[test]
fn independent_free_residue_order_empty_and_unknown_codes() {
    let mut db = empty_db();
    // Exact independent formulas: source adds each free residue minus water,
    // then one water, rather than one combined empirical-formula evaluation.
    let water = EmpiricalFormula::parse("H2O").unwrap().mono_mass();
    let free = [("A", "C3H7NO2"), ("G", "C2H5NO2"), ("M", "C5H11NO2S")];
    let mut expected = 0.;
    for (_, formula) in free {
        expected += EmpiricalFormula::parse(formula).unwrap().mono_mass() - water;
    }
    expected += water;
    assert_eq!(mono("AGM", &mut db), expected);
    assert_eq!(
        Peptidoform::default().mono_mass(&mut db).unwrap().value,
        water
    );
    for aa in ['B', 'Z', 'X'] {
        assert_eq!(
            chain(vec![element(aa, vec![])])
                .mono_mass(&mut db)
                .unwrap()
                .value,
            0.
        );
    }
    close(mono("BZ", &mut db), -water, 1e-12);
    for aa in ['J', 'U', 'O'] {
        assert!(
            chain(vec![element(aa, vec![])])
                .mono_mass(&mut db)
                .unwrap()
                .value
                .is_finite()
        );
    }
    for aa in ['a', 'é', '\0'] {
        let report = chain(vec![element(aa, vec![])])
            .mass_calculation_issues(&mut db)
            .unwrap();
        assert_eq!(
            report.value[0].issue_type,
            ConversionIssueType::UnsupportedFeature
        );
    }
}

#[test]
fn formula_mass_fallback_differs_from_resolution_charge_rules() {
    let mut db = empty_db();
    let water = EmpiricalFormula::parse("H2O").unwrap().mono_mass();
    for (text, charge) in [("", None), ("O", Some(2)), ("H1+", None)] {
        let pf = chain(vec![element('A', vec![formula(text, charge)])]);
        close(
            pf.mono_mass(&mut db).unwrap().value,
            base("A") + EmpiricalFormula::parse(text).unwrap().mono_mass(),
            1e-10,
        );
        assert!(
            pf.mass_calculation_issues(&mut db)
                .unwrap()
                .value
                .is_empty()
        );
    }
    let pf = chain(vec![element('A', vec![formula("not_a_formula", None)])]);
    assert!(pf.try_mono_mass(&mut db).unwrap().value.is_none());
    assert!(db.is_empty());
    let a = Peptidoform::parse("PEPT[INFO:x|Formula:Zn1:z+2]IDE").unwrap();
    let b = Peptidoform::parse("PEPT[Formula:Zn1:z+2|INFO:x]IDE").unwrap();
    assert_eq!(
        a.mono_mass(&mut db).unwrap().value,
        b.mono_mass(&mut db).unwrap().value
    );
    assert!(water > 18.);
}

#[test]
fn preserved_range_and_ambiguous_source_defects() {
    let mut db = empty_db();
    // CPP-008, CPP-009: these assertions record current upstream defects.
    close(mono("(M[+16]A)[+1]", &mut db), base("MA") + 1., 1e-10);
    close(mono("M[+16]A[+1]", &mut db), base("MA") + 17., 1e-10);
    let pf = Peptidoform::parse("(M[UnknownMod999]A)").unwrap();
    assert!(pf.can_calculate_mass(&mut db).unwrap().value);
    let pf = Peptidoform::parse("(?I[UnknownMod999]L)").unwrap();
    assert!(pf.can_calculate_mass(&mut db).unwrap().value);
    close(pf.mono_mass(&mut db).unwrap().value, base("I"), 1e-10);
    close(mono("(?I[+10]L)", &mut db), base("I") + 10., 1e-10);
    close(mono("(?LI[+10])", &mut db), base("L"), 1e-10);
    let different = Peptidoform::parse("(?AG)")
        .unwrap()
        .mass_calculation_issues(&mut db)
        .unwrap();
    assert_eq!(
        different.value[0].issue_type,
        ConversionIssueType::AmbiguousRegion
    );
    let empty = Peptidoform {
        sequence: vec![
            SequenceSection::AmbiguousRegion(AmbiguousRegion::default()),
            SequenceSection::Element(element('a', vec![])),
        ],
        ..Default::default()
    };
    assert_eq!(
        empty.mass_calculation_issues(&mut db).unwrap().value[0].position,
        Some(1)
    );
}

#[test]
fn source_label_only_order_and_dedup_group_boundaries() {
    let mut db = empty_db();
    let forward = PeptidoformIon::parse("K[#XL1]//K[+138.068#XL1]").unwrap();
    let reverse = PeptidoformIon::parse("K[+138.068#XL1]//K[#XL1]").unwrap();
    close(
        forward.mono_mass(&mut db).unwrap().value,
        2. * base("K"),
        1e-10,
    );
    close(
        reverse.mono_mass(&mut db).unwrap().value,
        2. * base("K") + 138.068,
        1e-10,
    );
    // Only the first alternative's label participates in deduplication.
    close(
        mono("K[INFO:x|+10#XL1]K[+10#XL1]", &mut db),
        base("KK") + 20.,
        1e-10,
    );

    let mut pf = Peptidoform::parse("K[+10#XL1]K[+10#XL1]").unwrap();
    let mut same = mass(10.);
    same.alternatives[0].1 = Some(Label {
        label_type: LabelType::Crosslink,
        identifier: "XL1".into(),
        score: None,
    });
    pf.unlocalised_mods.push(UnlocalisedMod {
        modifications: vec![same.clone()],
        occurrence: Some(2),
    });
    pf.labile_mods.push(LabileModification {
        modification: same.clone(),
    });
    pf.global_mods
        .push(GlobalModEntry::GlobalModification(GlobalModification {
            modification: same,
            locations: vec!["K".into()],
        }));
    close(
        pf.mono_mass(&mut db).unwrap().value,
        base("KK") + 60.,
        1e-10,
    );
    let mut pf = chain(vec![
        element('K', vec![mass(10.)]),
        element('K', vec![mass(10.)]),
    ]);
    for section in &mut pf.sequence {
        let SequenceSection::Element(e) = section else {
            panic!()
        };
        e.modifications[0].alternatives[0].1 = Some(Label {
            label_type: LabelType::Branch,
            identifier: "same".into(),
            score: Some(f64::NAN),
        });
    }
    close(
        pf.mono_mass(&mut db).unwrap().value,
        base("KK") + 20.,
        1e-10,
    );
}

#[test]
fn signed_occurrences_and_global_location_scope() {
    let mut db = empty_db();
    let mut pf = Peptidoform::parse("K(K)(?KL)").unwrap();
    pf.unlocalised_mods.push(UnlocalisedMod {
        modifications: vec![mass(2.)],
        occurrence: Some(-3),
    });
    // Only ordinary K is counted, once despite duplicate locations; "N-term"
    // and multi-byte locations are ignored. Isotope replacement has no mass effect.
    pf.global_mods
        .push(GlobalModEntry::GlobalModification(GlobalModification {
            modification: mass(3.),
            locations: vec!["N-term".into(), "é".into(), "K".into(), "K".into()],
        }));
    pf.global_mods
        .push(GlobalModEntry::IsotopeReplacement(IsotopeReplacement {
            isotope: "13C".into(),
        }));
    // K and L differ, so use equal I/L alternatives for a calculable region.
    if let SequenceSection::AmbiguousRegion(r) = &mut pf.sequence[2] {
        r.elements[0].amino_acid = 'I';
    }
    close(
        pf.mono_mass(&mut db).unwrap().value,
        base("KKI") - 3.,
        1e-10,
    );
}

#[test]
fn issue_positions_order_prefixes_and_error_precedence() {
    let mut db = empty_db();
    let mut pf = chain(vec![element('a', vec![named("missing")])]);
    pf.n_term_mods.push(named("n"));
    pf.c_term_mods.push(named("c"));
    pf.unlocalised_mods.push(UnlocalisedMod {
        modifications: vec![named("u")],
        occurrence: None,
    });
    pf.labile_mods.push(LabileModification {
        modification: named("l"),
    });
    let ion = PeptidoformIon {
        chains: vec![pf],
        is_chimeric: true,
        ..Default::default()
    };
    let issues = ion.mass_calculation_issues(&mut db).unwrap().value;
    assert_eq!(issues.len(), 6);
    assert!(
        issues
            .iter()
            .all(|i| i.description.starts_with("Chain 0: "))
    );
    assert_eq!(
        issues.iter().map(|i| i.position).collect::<Vec<_>>(),
        [Some(0), Some(0), Some(0), Some(0), None, None]
    );
    assert!(
        ion.mono_mass(&mut db)
            .unwrap_err()
            .to_string()
            .contains("Unknown amino acid")
    );
    assert!(
        ion.try_mono_mass(&mut db).unwrap().issues[0]
            .description
            .contains("chimeric")
    );
    assert_eq!(
        ion.try_mz(&mut db).unwrap().issues[0].description,
        "No charge state specified"
    );
    let empty = PeptidoformIon {
        is_chimeric: true,
        ..Default::default()
    };
    assert_eq!(empty.mono_mass(&mut db).unwrap().value, 0.);
    assert_eq!(empty.try_mono_mass(&mut db).unwrap().value, Some(0.));
    assert!(empty.can_calculate_mass(&mut db).unwrap().value);
}

#[test]
fn charges_use_only_integer_fields_and_check_overflow() {
    let mut db = empty_db();
    let pf = Peptidoform::parse("BZ").unwrap();
    let mass = pf.mono_mass(&mut db).unwrap().value;
    for charge in [-3, -1, 1, 2, i32::MIN] {
        close(
            pf.mz(charge, &mut db).unwrap().value,
            (mass + f64::from(charge) * PROTON_MASS_U) / f64::from(charge.unsigned_abs()),
            1e-12,
        );
    }
    let mut ion = PeptidoformIon {
        chains: vec![pf],
        charge: Some(ChargeState::Adducts(vec![
            AdductIon {
                formula: "not chemistry".into(),
                charge: 2,
                occurrence: Some(2),
            },
            AdductIon {
                formula: "Na".into(),
                charge: -1,
                occurrence: Some(2),
            },
        ])),
        ..Default::default()
    };
    close(
        ion.mz(&mut db).unwrap().value,
        (mass + 2. * PROTON_MASS_U) / 2.,
        1e-12,
    );
    ion.charge = Some(ChargeState::Adducts(vec![AdductIon {
        formula: "".into(),
        charge: i32::MAX,
        occurrence: Some(2),
    }]));
    assert!(ion.try_mz(&mut db).is_err());
    ion.charge = Some(ChargeState::Simple(0));
    assert_eq!(
        ion.try_mz(&mut db).unwrap().issues[0].description,
        "Charge state is zero"
    );
}

#[test]
fn existing_handles_empty_annotations_and_first_chemistry_selection() {
    let mut db = empty_db();
    let record = Arc::new(
        ResidueModification::from_record(ModificationRecord {
            name: "custom".into(),
            diff_mono_mass: 12.,
            ..Default::default()
        })
        .unwrap(),
    );
    let pf = chain(vec![element(
        'A',
        vec![Modification {
            alternatives: vec![],
            resolved_mod: Some(record.clone()),
        }],
    )]);
    close(pf.mono_mass(&mut db).unwrap().value, base("A") + 12., 1e-10);
    let mut empty = pf.clone();
    let SequenceSection::Element(e) = &mut empty.sequence[0] else {
        panic!()
    };
    e.modifications[0].resolved_mod = None;
    assert!(empty.try_mono_mass(&mut db).unwrap().value.is_none());
    let mut annotations = mass(20.);
    annotations.alternatives.insert(
        0,
        (
            ModificationTag::InfoTag(InfoTag {
                text: "note".into(),
            }),
            None,
        ),
    );
    let pf = chain(vec![element('A', vec![annotations])]);
    close(pf.mono_mass(&mut db).unwrap().value, base("A") + 20., 1e-10);
    let pf = chain(vec![element(
        'A',
        vec![modification(ModificationTag::PositionConstraint(
            PositionConstraint::default(),
        ))],
    )]);
    close(pf.mono_mass(&mut db).unwrap().value, base("A"), 1e-10);
    assert_eq!(Arc::strong_count(&record), 2);
}

#[test]
fn source_copy_passes_make_formula_interning_observable() {
    // CPP-020: valid public AST, not claimed to arise from text parsing of this NamedMod.
    let pf = chain(vec![
        element('M', vec![named("M[Formula:Cl101]")]),
        element('M', vec![formula("Cl101", None)]),
    ]);
    let before = pf.clone();
    let delta = EmpiricalFormula::parse("Cl101").unwrap().mono_mass();
    let mut strict = empty_db();
    assert!(pf.mono_mass(&mut strict).is_err());
    assert!(strict.is_empty());
    let mut db = empty_db();
    let attempt = pf.try_mono_mass(&mut db).unwrap();
    assert!(attempt.issues.is_empty());
    close(attempt.value.unwrap(), base("MM") + delta, 1e-8);
    // Validation resolved its second copy, but source calculates its first copy.
    close(
        pf.mono_mass(&mut db).unwrap().value,
        base("MM") + 2. * delta,
        1e-8,
    );
    assert_eq!(db.len(), 1);
    assert_eq!(pf, before);
    let ion = PeptidoformIon {
        chains: vec![pf.clone()],
        ..Default::default()
    };
    close(
        ion.try_mono_mass(&mut empty_db()).unwrap().value.unwrap(),
        base("MM") + 2. * delta,
        1e-8,
    );
    let mut db = empty_db();
    let report = pf.mass_calculation_issues(&mut db).unwrap();
    assert_eq!(report.value.len(), 1);
    assert_eq!(db.len(), 1);
}

#[test]
fn late_numerical_and_copy_budget_errors_do_not_publish_registry_or_ast() {
    let pf = chain(vec![
        element('M', vec![formula("Cl101", None)]),
        element('A', vec![mass(f64::MAX), mass(f64::MAX)]),
    ]);
    let before = pf.clone();
    let mut db = empty_db();
    assert!(pf.try_mono_mass(&mut db).is_err());
    assert!(db.is_empty());
    assert_eq!(pf, before);
    let huge = Peptidoform {
        name: Some("x".repeat(MAX_PROFORMA_MASS_TEXT_BYTES + 1)),
        ..Default::default()
    };
    assert!(huge.mono_mass(&mut db).is_err());
    assert!(db.is_empty());
    let mut ignored = Peptidoform::parse("A").unwrap();
    ignored.charge = Some(ChargeState::Simple(0));
    let before = ignored.clone();
    assert!(ignored.mono_mass(&mut db).is_ok());
    assert_eq!(ignored, before);
}
