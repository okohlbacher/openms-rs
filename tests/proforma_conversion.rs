use openms::chemistry::proforma::*;
use openms::chemistry::{
    AASequence, EmpiricalFormula, ModificationProvenance, ModificationRecord, ModificationsDB,
    ResidueModification, SequenceModification,
};
use std::sync::Arc;

fn parse(text: &str) -> Peptidoform {
    Peptidoform::parse(text).unwrap()
}
fn permissive(text: &str) -> AASequence {
    parse(text)
        .to_aa_sequence_with_policy(
            ConversionPolicy::BestEffort,
            &mut ModificationsDB::global().clone(),
        )
        .unwrap()
        .value
}
fn record(name: &str, mass: f64, formula: &str) -> Arc<ResidueModification> {
    Arc::new(
        ResidueModification::from_record(ModificationRecord {
            name: name.into(),
            full_id: format!("{name} (M)"),
            origin: Some('M'),
            diff_mono_mass: mass,
            provenance: ModificationProvenance::MassOnly,
            diff_formula: EmpiricalFormula::parse(formula).unwrap(),
            ..ModificationRecord::default()
        })
        .unwrap(),
    )
}
fn attached(records: &[Arc<ResidueModification>]) -> Peptidoform {
    Peptidoform {
        sequence: vec![SequenceSection::Element(SequenceElement {
            amino_acid: 'M',
            modifications: records
                .iter()
                .map(|record| Modification {
                    alternatives: vec![],
                    resolved_mod: Some(Arc::clone(record)),
                })
                .collect(),
        })],
        ..Peptidoform::default()
    }
}
fn first_record(pf: &Peptidoform) -> &Arc<ResidueModification> {
    let SequenceSection::Element(element) = &pf.sequence[0] else {
        panic!()
    };
    element.modifications[0].resolved_mod.as_ref().unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

#[test]
fn source_basic_and_all_policy_operations() {
    let mut db = ModificationsDB::global().clone();
    for text in ["PEPTIDE", "PEPM[UNIMOD:35]TIDE", "[UNIMOD:1]-PEPTIDE"] {
        let pf = parse(text);
        assert!(pf.is_representable_as_aa_sequence(&mut db).unwrap().value);
        let a = pf.to_aa_sequence(&mut db).unwrap().value;
        for policy in [
            ConversionPolicy::FailOnLoss,
            ConversionPolicy::DropUnlocalised,
            ConversionPolicy::BestEffort,
        ] {
            assert_eq!(
                pf.to_aa_sequence_with_policy(policy, &mut db)
                    .unwrap()
                    .value,
                a
            );
        }
        assert_eq!(
            a.as_str(),
            if text.contains("PEPM") {
                "PEPMTIDE"
            } else {
                "PEPTIDE"
            }
        );
    }
    assert_eq!(
        permissive("PEPM[UNIMOD:35]TIDE")
            .residue_modification(3)
            .unwrap()
            .unwrap()
            .record_id(),
        Some(35)
    );
    let pf = parse("[+79.966331]?PEPTIDE");
    assert!(!pf.is_representable_as_aa_sequence(&mut db).unwrap().value);
    assert!(pf.to_aa_sequence(&mut db).is_err());
    for policy in [
        ConversionPolicy::DropUnlocalised,
        ConversionPolicy::BestEffort,
    ] {
        assert_eq!(
            pf.to_aa_sequence_with_policy(policy, &mut db)
                .unwrap()
                .value,
            AASequence::parse("PEPTIDE").unwrap()
        );
    }
    let unknown = parse("PEPM[+123456.789]TIDE");
    assert!(!permissive("PEPM[+123456.789]TIDE").is_modified());
    assert!(unknown.to_aa_sequence(&mut db).is_err());
    for text in [
        "PEPM[+15.9949]TIDE",
        "PEPC[+57.02146]TIDE",
        "PEPS[+79.96633]TIDE",
    ] {
        assert!(permissive(text).is_modified());
    }
}

#[test]
fn source_formula_interning_exact_identity_and_unusable_tags() {
    let mut db = ModificationsDB::global().clone();
    let mut previous = None;
    for text in [
        "PEPA[Formula:C2H2O]TIDE",
        "PEPA[Formula:C2H2O1]TIDE",
        "PEPA[Formula:H2C2O1]TIDE",
    ] {
        let seq = parse(text).to_aa_sequence(&mut db).unwrap().value;
        let SequenceModification::Known(m) = seq.residue_modification(3).unwrap().unwrap() else {
            panic!()
        };
        assert_eq!(m.full_id(), "A[Formula:C2H2O1]");
        assert_eq!(m.provenance(), ModificationProvenance::MassOnly);
        assert!(m.name().is_empty());
        if let Some(ref old) = previous {
            assert!(Arc::ptr_eq(old, m));
        }
        previous = Some(Arc::clone(m));
    }
    let charged = parse("PEPM[Formula:O:z+1]TIDE");
    assert!(charged.to_aa_sequence(&mut db).is_err());
    assert!(
        !charged
            .to_aa_sequence_with_policy(ConversionPolicy::BestEffort, &mut db)
            .unwrap()
            .value
            .is_modified()
    );
    let formula = permissive("[Formula:C2H2O]-PEPTIDE");
    close(
        formula.mono_mass().unwrap(),
        permissive("[UNIMOD:1]-PEPTIDE").mono_mass().unwrap(),
        2e-6,
    );
    assert!(
        !parse("pepm[Formula:O]tide")
            .aa_sequence_conversion_issues(&mut db)
            .unwrap()
            .value
            .is_empty()
    );
    assert!(
        parse("peptide")
            .is_representable_as_aa_sequence(&mut db)
            .unwrap()
            .value
    );
    assert!(parse("peptide").to_aa_sequence(&mut db).is_err());
}

#[test]
fn source_combination_repetition_cancellation_and_losses() {
    let oxygen = EmpiricalFormula::parse("O").unwrap();
    let base = AASequence::parse("PEPMTIDE").unwrap();
    for (text, count) in [
        ("PEPM[Oxidation][Formula:O]TIDE", 2),
        ("PEPM[Formula:O][Formula:O][Formula:O]TIDE", 3),
        ("PEPM[Oxidation][Oxidation]TIDE", 2),
    ] {
        let result = permissive(text);
        let expected = base
            .formula()
            .unwrap()
            .checked_add(&EmpiricalFormula::parse(&format!("O{count}")).unwrap())
            .unwrap();
        assert_eq!(result.formula().unwrap(), expected);
        close(
            result.mono_mass().unwrap(),
            base.mono_mass().unwrap() + f64::from(count) * oxygen.mono_mass(),
            1e-10,
        );
        assert!(
            parse(text)
                .to_aa_sequence(&mut ModificationsDB::global().clone())
                .is_err()
        );
        assert_eq!(
            result
                .residue_modification(3)
                .unwrap()
                .unwrap()
                .known()
                .unwrap()
                .provenance(),
            ModificationProvenance::MassOnly
        );
    }
    let cancel = permissive("PEPM[Formula:H2O][Formula:H-2O-1]TIDE");
    assert_eq!(cancel, base);
    let loss = permissive("PEPM[Formula:O][Formula:H-2]TIDE");
    assert_eq!(
        loss.formula().unwrap(),
        base.formula()
            .unwrap()
            .checked_add(&EmpiricalFormula::parse("H-2O").unwrap())
            .unwrap()
    );
    let oxygen_record = record("same", oxygen.mono_mass(), "O");
    let seq = attached(&[Arc::clone(&oxygen_record), oxygen_record])
        .to_aa_sequence(&mut ModificationsDB::default())
        .unwrap()
        .value;
    assert_eq!(
        seq.formula().unwrap(),
        AASequence::parse("M")
            .unwrap()
            .formula()
            .unwrap()
            .checked_add(&EmpiricalFormula::parse("O2").unwrap())
            .unwrap()
    );
}

#[test]
fn declared_mass_combination_thresholds_and_source_record_collisions() {
    let oxygen = EmpiricalFormula::parse("O").unwrap().mono_mass();
    let cases = [
        (
            vec![record("a", 100., "O"), record("b", 200., "O")],
            300.,
            true,
        ),
        (
            vec![record("a", 100., "O"), record("b", 200., "")],
            300.,
            false,
        ),
    ];
    for (records, mass, warning) in cases {
        let mut db = ModificationsDB::default();
        let result = attached(&records).to_aa_sequence(&mut db).unwrap();
        let m = result
            .value
            .residue_modification(0)
            .unwrap()
            .unwrap()
            .known()
            .unwrap();
        assert!(m.diff_formula().is_empty());
        assert_eq!(m.diff_mono_mass(), mass);
        assert_eq!(
            result
                .warnings
                .iter()
                .any(|w| matches!(w, ConversionWarning::FormulaMassDisagreement { .. })),
            warning
        );
        assert!(m.name().is_empty());
        assert_eq!(m.full_id(), "M[+300.0]");
    }
    for (sum, modified) in [
        (1e-6, false),
        (f64::from_bits(1e-6f64.to_bits() + 1), true),
        (-1e-6, false),
    ] {
        let records = [record("a", 0., ""), record("b", sum, "")];
        assert_eq!(
            attached(&records)
                .to_aa_sequence(&mut ModificationsDB::default())
                .unwrap()
                .value
                .is_modified(),
            modified
        );
    }
    let records = [record("a", oxygen, "O"), record("b", oxygen + 0.0005, "O")];
    assert!(
        !attached(&records)
            .to_aa_sequence(&mut ModificationsDB::default())
            .unwrap()
            .value
            .residue_modification(0)
            .unwrap()
            .unwrap()
            .known()
            .unwrap()
            .diff_formula()
            .is_empty()
    );
    let collision = ResidueModification::from_record(ModificationRecord {
        full_id: "M[+300.0]".into(),
        name: "collision".into(),
        origin: Some('M'),
        diff_mono_mass: 17.,
        ..ModificationRecord::default()
    })
    .unwrap();
    let mut db = ModificationsDB::from_records(vec![collision]).unwrap();
    let result = attached(&[record("a", 100., ""), record("b", 200., "")])
        .to_aa_sequence(&mut db)
        .unwrap();
    assert_eq!(
        result
            .value
            .residue_modification(0)
            .unwrap()
            .unwrap()
            .name(),
        "collision"
    );
}

#[test]
fn source_annotation_alternatives_definition_and_terminal_rules() {
    let mut db = ModificationsDB::global().clone();
    for text in [
        "PEPM[Formula:O|INFO:my note]TIDE",
        "PEPM[INFO:my note|Formula:O]TIDE",
    ] {
        assert!(
            parse(text)
                .aa_sequence_conversion_issues(&mut db)
                .unwrap()
                .value
                .is_empty()
        );
        assert!(
            parse(text)
                .to_aa_sequence(&mut db)
                .unwrap()
                .value
                .is_modified()
        );
    }
    assert!(
        parse("PEPM[Oxidation|Phospho]TIDE")
            .aa_sequence_conversion_issues(&mut db)
            .unwrap()
            .value
            .iter()
            .any(|i| i.issue_type == ConversionIssueType::AlternativeMods)
    );
    assert!(
        parse("PEPM[INFO:|UNIMOD:999999]TIDE")
            .aa_sequence_conversion_issues(&mut db)
            .unwrap()
            .value
            .iter()
            .any(|i| i.issue_type == ConversionIssueType::UnresolvedMod)
    );
    assert!(
        parse("PEPM[UNIMOD:447]TIDE")
            .to_aa_sequence(&mut db)
            .is_err()
    );
    let first = permissive("[#XL1][UNIMOD:1]-PEPTIDE");
    assert_eq!(
        first.n_terminal_modification().unwrap().record_id(),
        Some(1)
    );
    assert!(
        parse("[UNIMOD:1][Formula:O]-PEPTIDE")
            .to_aa_sequence(&mut db)
            .is_err()
    );
    let seq = permissive("PEPM[Formula:C2H2O|INFO:Oxidation]TIDE");
    assert_eq!(
        seq.residue_modification(3)
            .unwrap()
            .unwrap()
            .known()
            .unwrap()
            .diff_formula(),
        &EmpiricalFormula::parse("C2H2O").unwrap()
    );
}

#[test]
fn source_full_issue_order_and_lossy_flattening() {
    let mut pf = parse("<[Formula:O]@M>[+1]?{Glycan:Hex}M[unknown|other#XL1](?IL)(MA)[+1]");
    pf.n_term_mods = parse("[unknown|other]-M").n_term_mods;
    pf.c_term_mods = parse("M-[unknown|other]").c_term_mods;
    let issues = pf
        .aa_sequence_conversion_issues(&mut ModificationsDB::global().clone())
        .unwrap()
        .value;
    let kinds: Vec<_> = issues.iter().map(|i| i.issue_type).collect();
    assert_eq!(
        kinds,
        vec![
            ConversionIssueType::UnlocalisedMod,
            ConversionIssueType::LabileMod,
            ConversionIssueType::GlobalMod,
            ConversionIssueType::UnresolvedMod,
            ConversionIssueType::AlternativeMods,
            ConversionIssueType::CrossLink,
            ConversionIssueType::AmbiguousRegion,
            ConversionIssueType::ModifiedRange,
            ConversionIssueType::UnresolvedMod,
            ConversionIssueType::AlternativeMods,
            ConversionIssueType::UnresolvedMod,
            ConversionIssueType::AlternativeMods
        ]
    );
    assert_eq!(issues[6].position, Some(1));
    assert_eq!(issues[7].position, Some(2));
    assert_eq!(issues[8].position, None);
    let a = pf
        .to_aa_sequence_with_policy(
            ConversionPolicy::BestEffort,
            &mut ModificationsDB::global().clone(),
        )
        .unwrap()
        .value;
    let b = pf
        .to_aa_sequence_with_policy(
            ConversionPolicy::DropUnlocalised,
            &mut ModificationsDB::global().clone(),
        )
        .unwrap()
        .value;
    assert_eq!(a, b);
    assert_eq!(a.as_str(), "MIMA");
    assert!(!a.is_modified());
    let iso = parse("<13C>M");
    assert!(
        iso.is_representable_as_aa_sequence(&mut ModificationsDB::default())
            .unwrap()
            .value
    );
    assert_eq!(
        iso.to_aa_sequence(&mut ModificationsDB::default())
            .unwrap()
            .value,
        AASequence::parse("M").unwrap()
    );
}

#[test]
fn cpp_029_030_031_predicate_and_cursor_defects_are_preserved() {
    // CPP-029: the issue collector accepts annotations-only brackets, but strict attachment rejects null.
    let pf = parse("M[INFO:note]");
    let mut db = ModificationsDB::default();
    assert!(pf.is_representable_as_aa_sequence(&mut db).unwrap().value);
    assert!(
        pf.to_aa_sequence(&mut db)
            .unwrap_err()
            .to_string()
            .contains("Unresolved modification at position 0")
    );
    assert_eq!(
        pf.to_aa_sequence_with_policy(ConversionPolicy::BestEffort, &mut db)
            .unwrap()
            .value
            .as_str(),
        "M"
    );
    // CPP-030: terminal crosslink labels are absent from diagnostics and then disappear.
    let mut terminal = parse("M");
    terminal.n_term_mods = vec![Modification {
        alternatives: vec![(
            ModificationTag::InfoTag(InfoTag::default()),
            Some(Label {
                label_type: LabelType::Crosslink,
                identifier: "link".into(),
                score: None,
            }),
        )],
        resolved_mod: None,
    }];
    assert!(
        terminal
            .is_representable_as_aa_sequence(&mut db)
            .unwrap()
            .value
    );
    assert!(
        !terminal
            .to_aa_sequence(&mut db)
            .unwrap()
            .value
            .is_modified()
    );
    // CPP-031: source increments the attachment cursor for an empty ambiguous section.
    let mut empty = attached(&[record("custom", 1., "")]);
    empty.sequence.insert(
        0,
        SequenceSection::AmbiguousRegion(AmbiguousRegion { elements: vec![] }),
    );
    assert!(
        empty
            .to_aa_sequence_with_policy(ConversionPolicy::BestEffort, &mut db)
            .unwrap_err()
            .to_string()
            .contains("index out of bounds")
    );
}

#[test]
fn reverse_source_precedence_and_named_definition_identity() {
    let mut db = ModificationsDB::global().clone();
    let base = AASequence::parse("PEPM(Oxidation)TIDE").unwrap();
    let pf = Peptidoform::from_aa_sequence(&base).unwrap();
    assert_eq!(
        pf.to_text(WriteMode::Canonical).unwrap(),
        "PEPM[UNIMOD:35]TIDE"
    );
    let round = parse(&pf.to_text(WriteMode::Lossless).unwrap())
        .to_aa_sequence(&mut db)
        .unwrap()
        .value;
    assert_eq!(round, base);
    let definition = ResidueModification::from_record(ModificationRecord {
        name: "Test:Adduct".into(),
        origin: Some('M'),
        diff_formula: EmpiricalFormula::parse("C2H2O").unwrap(),
        diff_mono_mass: 42.010565,
        provenance: ModificationProvenance::Defined,
        ..ModificationRecord::default()
    })
    .unwrap();
    let m = db.register_definition(&definition).unwrap();
    let sequence = attached(&[Arc::clone(&m)])
        .to_aa_sequence(&mut db)
        .unwrap()
        .value;
    let written = Peptidoform::from_aa_sequence(&sequence).unwrap();
    assert_eq!(
        written.to_text(WriteMode::Lossless).unwrap(),
        "M[Formula:C2H2O1|INFO:Test:Adduct]"
    );
    assert!(Arc::ptr_eq(first_record(&written), &m));
    let back = parse(&written.to_text(WriteMode::Lossless).unwrap())
        .to_aa_sequence(&mut db)
        .unwrap()
        .value;
    let SequenceModification::Known(back_m) = back.residue_modification(0).unwrap().unwrap() else {
        panic!()
    };
    assert!(Arc::ptr_eq(back_m, &m));
    let named = attached(&[record("ordinary", 123.45, "")])
        .to_aa_sequence(&mut db)
        .unwrap()
        .value;
    assert_eq!(
        Peptidoform::from_aa_sequence(&named)
            .unwrap()
            .to_text(WriteMode::Lossless)
            .unwrap(),
        "M[ordinary]"
    );
}

#[test]
fn source_formula_complete_and_chemical_roundtrip_examples() {
    let source = "AEADNLDDK[Formula:C9H11N2O8P|INFO:NuXL:U-H2O]K";
    let mut db = ModificationsDB::global().clone();
    let sequence = parse(source).to_aa_sequence(&mut db).unwrap().value;
    assert_eq!(sequence.as_str(), "AEADNLDDKK");
    close(sequence.mono_mass().unwrap(), 1423.5504442334, 2e-6);
    assert_eq!(
        sequence.formula().unwrap(),
        EmpiricalFormula::parse("C54H86N15O28P").unwrap()
    );
    let out = Peptidoform::from_aa_sequence(&sequence).unwrap();
    assert!(
        out.to_text(WriteMode::Canonical)
            .unwrap()
            .contains("Formula:")
    );
    let back = parse(&out.to_text(WriteMode::Canonical).unwrap())
        .to_aa_sequence(&mut db)
        .unwrap()
        .value;
    assert_eq!(back.formula().unwrap(), sequence.formula().unwrap());
    for text in [
        "PEPTIDE",
        "(Acetyl)PEPM(Oxidation)TIDE",
        "PEPM[+0.00335]TIDE",
        "PEPM[-0.00335]TIDE",
        "PEPM[+12345.6789]TIDE",
        "AEADNLDDK[+306.025304840900048]K",
    ] {
        let sequence = AASequence::parse(text).unwrap();
        let pf = Peptidoform::from_aa_sequence(&sequence).unwrap();
        let written = pf.to_text(WriteMode::Lossless).unwrap();
        assert!(!written.contains("[]"));
        let _ = parse(&written);
        close(
            pf.mono_mass(&mut ModificationsDB::global().clone())
                .unwrap()
                .value,
            sequence.mono_mass().unwrap(),
            2e-6,
        );
    }
}

#[test]
fn native_mass_tags_exact_fixed_text_and_zero_base_placeholders() {
    for (mass, expected) in [
        (0., "+0"),
        (-0., "+0"),
        (0.00335, "+0.00335"),
        (-0.00335, "-0.00335"),
        (12345.6789, "+12345.6789"),
    ] {
        let m = Arc::new(
            ResidueModification::from_record(ModificationRecord {
                full_id: "M[anonymous]".into(),
                origin: Some('M'),
                diff_mono_mass: mass,
                ..ModificationRecord::default()
            })
            .unwrap(),
        );
        let seq = attached(&[m])
            .to_aa_sequence(&mut ModificationsDB::default())
            .unwrap()
            .value;
        let pf = Peptidoform::from_aa_sequence(&seq).unwrap();
        let SequenceSection::Element(elem) = &pf.sequence[0] else {
            panic!()
        };
        let ModificationTag::MassDelta(tag) = &elem.modifications[0].alternatives[0].0 else {
            panic!()
        };
        assert_eq!(tag.original_text, expected);
    }
    for text in [
        "M[+0.00000001]",
        "M[+123456789]",
        "M[147.035405]",
        "X[999]",
        "B[999]",
        "Z[999]",
        ".[12.345]M",
        "M.[12.345]",
    ] {
        let seq = AASequence::parse_with_registry(text, &ModificationsDB::default()).unwrap();
        let pf = Peptidoform::from_aa_sequence(&seq).unwrap();
        close(
            pf.mono_mass(&mut ModificationsDB::default()).unwrap().value,
            seq.mono_mass().unwrap(),
            1e-8,
        );
        if text.starts_with(['B', 'Z', 'X']) {
            close(
                first_record(&pf).diff_mono_mass(),
                999. + EmpiricalFormula::parse("H2O").unwrap().mono_mass(),
                1e-12,
            );
            assert!(
                seq.residue_modification(0)
                    .unwrap()
                    .unwrap()
                    .diff_mono_mass()
                    .is_err()
            );
        }
    }
}

#[test]
fn immutable_inputs_registry_transaction_and_detached_handles() {
    let pf = parse("M[Formula:Cl101]A[not-found]");
    let snapshot = pf.clone();
    let mut db = ModificationsDB::default();
    assert!(pf.to_aa_sequence(&mut db).is_err());
    assert!(db.is_empty());
    assert_eq!(pf, snapshot);
    let issues = pf.aa_sequence_conversion_issues(&mut db).unwrap();
    assert!(!issues.value.is_empty());
    assert_eq!(db.len(), 1);
    assert_eq!(pf, snapshot);
    let sequence = parse("M[Formula:Cl101]")
        .to_aa_sequence(&mut db)
        .unwrap()
        .value;
    let detached = Peptidoform::from_aa_sequence(&sequence).unwrap();
    drop(db);
    drop(sequence);
    assert_eq!(first_record(&detached).full_id(), "M[Formula:Cl101]");
    let empty = Peptidoform {
        n_term_mods: parse("[Formula:O]-M").n_term_mods,
        ..Peptidoform::default()
    };
    let result = empty
        .to_aa_sequence(&mut ModificationsDB::default())
        .unwrap()
        .value;
    assert!(result.is_empty());
    assert!(result.n_terminal_modification().is_some());
    assert_eq!(result.mono_mass().unwrap(), 0.);
    assert_eq!(
        Peptidoform::from_aa_sequence(&result)
            .unwrap()
            .n_term_mods
            .len(),
        1
    );
}

#[test]
fn executed_source_formatter_fixture_through_complete_reverse_api() {
    let mut rows = 0;
    for line in include_str!("data/proforma_conversion_mass_text.tsv")
        .lines()
        .filter(|line| !line.starts_with('#'))
    {
        let (bits, expected) = line.split_once('\t').unwrap();
        let mass = f64::from_bits(u64::from_str_radix(bits, 16).unwrap());
        let m = Arc::new(
            ResidueModification::from_record(ModificationRecord {
                full_id: "M[probe]".into(),
                origin: Some('M'),
                // Independent source fields: valid cached absolute mass while
                // the reverse conversion reads any finite declared delta.
                mono_mass: 1000.,
                diff_mono_mass: mass,
                provenance: ModificationProvenance::MassOnly,
                ..ModificationRecord::default()
            })
            .unwrap(),
        );
        let sequence = attached(&[m])
            .to_aa_sequence(&mut ModificationsDB::default())
            .unwrap()
            .value;
        let pf = Peptidoform::from_aa_sequence(&sequence).unwrap();
        let SequenceSection::Element(element) = &pf.sequence[0] else {
            panic!()
        };
        let ModificationTag::MassDelta(delta) = &element.modifications[0].alternatives[0].0 else {
            panic!()
        };
        assert_eq!(delta.original_text, expected, "bits {bits}");
        assert_eq!(
            pf.to_text(WriteMode::Lossless).unwrap(),
            format!("M[{expected}]"),
            "bits {bits}"
        );
        assert_eq!(
            parse(&format!("M[{expected}]"))
                .to_text(WriteMode::Lossless)
                .unwrap(),
            format!("M[{expected}]")
        );
        rows += 1;
    }
    assert_eq!(rows, 1078);
}

#[test]
fn cpp_037_charged_formula_and_reverse_defined_mass_precedence() {
    let charged = EmpiricalFormula::parse("O").unwrap().with_charge(1);
    let neutral = EmpiricalFormula::parse("O").unwrap();
    let declared_sum = charged.mono_mass() + neutral.mono_mass();
    let a = Arc::new(
        ResidueModification::from_record(ModificationRecord {
            full_id: "a".into(),
            origin: Some('M'),
            diff_mono_mass: charged.mono_mass(),
            diff_formula: charged,
            ..ModificationRecord::default()
        })
        .unwrap(),
    );
    let b = Arc::new(
        ResidueModification::from_record(ModificationRecord {
            full_id: "b".into(),
            origin: Some('M'),
            diff_mono_mass: neutral.mono_mass(),
            diff_formula: neutral,
            ..ModificationRecord::default()
        })
        .unwrap(),
    );
    let sequence = attached(&[a, b])
        .to_aa_sequence(&mut ModificationsDB::default())
        .unwrap()
        .value;
    let m = sequence
        .residue_modification(0)
        .unwrap()
        .unwrap()
        .known()
        .unwrap();
    // CPP-037: source formula_sum.toString drops charge before neutral interning.
    // This is source-reviewed conversion behavior, not an executed C++ conversion.
    assert_eq!(m.diff_formula(), &EmpiricalFormula::parse("O2").unwrap());
    close(
        declared_sum - m.diff_mono_mass(),
        openms::chemistry::PROTON_MASS_U,
        1e-12,
    );
    let def = Arc::new(
        ResidueModification::from_record(ModificationRecord {
            name: "defined-mass".into(),
            origin: Some('M'),
            diff_mono_mass: 123.4567,
            ..ModificationRecord::default()
        })
        .unwrap(),
    );
    let sequence = attached(&[def])
        .to_aa_sequence(&mut ModificationsDB::default())
        .unwrap()
        .value;
    assert_eq!(
        Peptidoform::from_aa_sequence(&sequence)
            .unwrap()
            .to_text(WriteMode::Lossless)
            .unwrap(),
        "M[+123.4567|INFO:defined-mass]"
    );
    let cv = Arc::new(
        ResidueModification::from_record(ModificationRecord {
            name: "defined-cv".into(),
            record_id: Some(35),
            origin: Some('M'),
            diff_mono_mass: 123.4567,
            ..ModificationRecord::default()
        })
        .unwrap(),
    );
    let sequence = attached(&[cv])
        .to_aa_sequence(&mut ModificationsDB::default())
        .unwrap()
        .value;
    assert_eq!(
        Peptidoform::from_aa_sequence(&sequence)
            .unwrap()
            .to_text(WriteMode::Lossless)
            .unwrap(),
        "M[UNIMOD:35]"
    );
}
