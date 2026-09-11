use openms::chemistry::proforma::*;
use openms::chemistry::{
    ModificationRecord, ModificationsDB, ResidueModification, TermSpecificity,
};
use std::sync::Arc;

fn element(amino_acid: char) -> SequenceElement {
    SequenceElement {
        amino_acid,
        modifications: vec![],
    }
}
fn peptide(sequence: &str) -> Peptidoform {
    Peptidoform {
        sequence: sequence
            .chars()
            .map(|c| SequenceSection::Element(element(c)))
            .collect(),
        ..Default::default()
    }
}
fn modification(tag: ModificationTag) -> Modification {
    Modification {
        alternatives: vec![(tag, None)],
        ..Default::default()
    }
}
fn named(name: &str) -> Modification {
    modification(ModificationTag::NamedMod(NamedMod {
        name: name.into(),
        ..Default::default()
    }))
}
fn attach(pf: &mut Peptidoform, index: usize, modification: Modification) {
    let SequenceSection::Element(element) = &mut pf.sequence[index] else {
        panic!("ordinary test element")
    };
    element.modifications.push(modification);
}
fn label(score: Option<f64>) -> Label {
    Label {
        label_type: LabelType::Ambiguous,
        identifier: "g1".into(),
        score,
    }
}
fn scored(tag: ModificationTag, score: f64) -> Modification {
    Modification {
        alternatives: vec![(tag, Some(label(Some(score))))],
        ..Default::default()
    }
}
fn delta(value: f64, spelling: &str) -> ModificationTag {
    ModificationTag::MassDelta(MassDelta {
        mass: value,
        original_text: spelling.into(),
        ..Default::default()
    })
}
fn info(text: &str) -> ModificationTag {
    ModificationTag::InfoTag(InfoTag { text: text.into() })
}

#[test]
fn source_literal_writer_and_precision_assertions() {
    // ProFormaParser_test.cpp:1146–1167 and 1337–1407. These construct the
    // source AST directly; no source/native parser execution is claimed.
    assert_eq!(
        peptide("PEPTIDE").to_text(WriteMode::default()).unwrap(),
        "PEPTIDE"
    );
    let mut pf = peptide("EMK");
    attach(
        &mut pf,
        1,
        modification(ModificationTag::CvAccession(CvAccession {
            database: CvDatabase::Unimod,
            accession: "35".into(),
        })),
    );
    for mode in [WriteMode::Lossless, WriteMode::Canonical] {
        assert_eq!(pf.to_text(mode).unwrap(), "EM[UNIMOD:35]K");
    }
    let mut termini = peptide("PEPTIDE");
    termini.n_term_mods.push(named("Acetyl"));
    termini.c_term_mods.push(named("Amidated"));
    assert_eq!(
        termini.to_text(WriteMode::Lossless).unwrap(),
        "[Acetyl]-PEPTIDE-[Amidated]"
    );
    for (mass, spelling, canonical) in [
        (15.99, "+15.99", "EM[+15.9900]K"),
        (15.99491234, "+15.99491234", "EM[+15.9949]K"),
    ] {
        let mut pf = peptide("EMK");
        attach(&mut pf, 1, modification(delta(mass, spelling)));
        assert_eq!(
            pf.to_text(WriteMode::Lossless).unwrap(),
            format!("EM[{spelling}]K")
        );
        assert_eq!(pf.to_text(WriteMode::Canonical).unwrap(), canonical);
    }
    let mut pf = peptide("EMK");
    attach(&mut pf, 1, named("Oxidation"));
    assert_eq!(pf.to_text(WriteMode::Canonical).unwrap(), "EM[Oxidation]K");
}

#[test]
fn source_literal_crosslinks_chimeric_charges_and_names() {
    // Exact source assertions at 298, 548, 571, 635 and 649.
    let mut pf = peptide("EMEVTKSESPEK");
    let xl = Label {
        label_type: LabelType::Crosslink,
        identifier: "XL1".into(),
        score: None,
    };
    attach(
        &mut pf,
        5,
        Modification {
            alternatives: vec![(
                ModificationTag::CvAccession(CvAccession {
                    database: CvDatabase::Xlmod,
                    accession: "02001".into(),
                }),
                Some(xl.clone()),
            )],
            ..Default::default()
        },
    );
    attach(
        &mut pf,
        11,
        Modification {
            alternatives: vec![(info(""), Some(xl))],
            ..Default::default()
        },
    );
    assert_eq!(
        pf.to_text(WriteMode::Lossless).unwrap(),
        "EMEVTK[XLMOD:02001#XL1]SESPEK[#XL1]"
    );
    let mut ion = PeptidoformIon {
        chains: vec![peptide("EMEVEESPEK"), peptide("ELVISLIVER")],
        is_chimeric: true,
        ..Default::default()
    };
    assert_eq!(
        ion.to_text(WriteMode::Lossless).unwrap(),
        "EMEVEESPEK+ELVISLIVER"
    );
    ion.chains[0].charge = Some(ChargeState::Simple(2));
    ion.chains[1].charge = Some(ChargeState::Simple(3));
    assert_eq!(
        ion.to_text(WriteMode::Lossless).unwrap(),
        "EMEVEESPEK/2+ELVISLIVER/3"
    );
    ion.chains = vec![peptide("AANSIPYQVSLNS"), peptide("AKEQFERQTA")];
    ion.chains[0].name = Some("Trypsin".into());
    ion.chains[1].name = Some("Keratin".into());
    assert_eq!(
        ion.to_text(WriteMode::Lossless).unwrap(),
        "(>Trypsin)AANSIPYQVSLNS+(>Keratin)AKEQFERQTA"
    );
    let mut pf = peptide("PEPTIDE");
    pf.name = Some("sp|P12345|PROT_HUMAN".into());
    assert_eq!(
        pf.to_text(WriteMode::Lossless).unwrap(),
        "(>sp|P12345|PROT_HUMAN)PEPTIDE"
    );
}

#[test]
fn independent_complete_structure_preserves_order_and_annotations() {
    let formula = FormulaTag {
        formula_string: "[13C]2H-1".into(),
        charge: Some(-2),
    };
    let mut pf = peptide("A");
    pf.name = Some("α peptide".into());
    pf.global_mods = vec![
        GlobalModEntry::IsotopeReplacement(IsotopeReplacement {
            isotope: "15N".into(),
        }),
        GlobalModEntry::GlobalModification(GlobalModification {
            modification: named("z-last"),
            locations: vec!["K".into(), "N-term".into(), "K".into()],
        }),
        GlobalModEntry::GlobalModification(GlobalModification {
            modification: named("a-first"),
            locations: vec![],
        }),
    ];
    pf.unlocalised_mods.push(UnlocalisedMod {
        modifications: vec![named("Phospho"), modification(info("β"))],
        occurrence: Some(-2),
    });
    pf.labile_mods.push(LabileModification {
        modification: modification(ModificationTag::GlycanComposition(GlycanComposition {
            components: vec![
                (GlycanComponent::Name("Hex".into()), 1),
                (GlycanComponent::Formula(formula.clone()), 0),
                (GlycanComponent::Name("HexNAc".into()), -3),
            ],
        })),
    });
    pf.n_term_mods = vec![named("Acetyl"), named("Other")];
    let mut alternative = element('Q');
    alternative
        .modifications
        .push(modification(ModificationTag::FormulaTag(formula)));
    pf.sequence
        .push(SequenceSection::AmbiguousRegion(AmbiguousRegion {
            elements: vec![element('D'), alternative],
        }));
    pf.sequence
        .push(SequenceSection::ModifiedRange(ModifiedRange {
            elements: vec![element('P'), element('E')],
            modifications: vec![
                modification(ModificationTag::PositionConstraint(PositionConstraint {
                    residues: vec!['M', 'K', 'M'],
                    n_term: true,
                    c_term: true,
                })),
                named("range"),
            ],
        }));
    pf.c_term_mods = vec![named("Amidated"), Modification::default()];
    let expected = "(>α peptide)<15N><[z-last]@K,N-term,K><[a-first]>[Phospho][INFO:β]^-2?{Glycan:Hex1Formula:[13C]2H-1:z-20HexNAc-3}[Acetyl][Other]-A(?DQ[Formula:[13C]2H-1:z-2])(PE)[Position:N-term,C-term,MKM][range]-[Amidated][]";
    for mode in [WriteMode::Lossless, WriteMode::Canonical] {
        assert_eq!(pf.to_text(mode).unwrap(), expected);
    }
}

#[test]
fn all_database_and_mass_hints_and_signed_counts() {
    let mut pf = peptide("X");
    for (db, prefix, hint) in [
        (CvDatabase::Unimod, "UNIMOD", "U"),
        (CvDatabase::Mod, "MOD", "M"),
        (CvDatabase::Resid, "RESID", "R"),
        (CvDatabase::Xlmod, "XLMOD", "X"),
        (CvDatabase::Gno, "GNO", "G"),
    ] {
        let mut one = peptide("X");
        attach(
            &mut one,
            0,
            Modification {
                alternatives: vec![
                    (
                        ModificationTag::CvAccession(CvAccession {
                            database: db,
                            accession: "00abc".into(),
                        }),
                        None,
                    ),
                    (
                        ModificationTag::NamedMod(NamedMod {
                            cv_hint: Some(db),
                            name: "ñ".into(),
                        }),
                        None,
                    ),
                ],
                ..Default::default()
            },
        );
        assert_eq!(
            one.to_text(WriteMode::Canonical).unwrap(),
            format!("X[{prefix}:00abc|{hint}:ñ]")
        );
    }
    for (source, prefix) in [
        (MassDeltaSource::None, ""),
        (MassDeltaSource::Obs, "Obs:"),
        (MassDeltaSource::U, "U:"),
        (MassDeltaSource::M, "M:"),
        (MassDeltaSource::R, "R:"),
        (MassDeltaSource::X, "X:"),
        (MassDeltaSource::G, "G:"),
    ] {
        attach(
            &mut pf,
            0,
            modification(ModificationTag::MassDelta(MassDelta {
                source,
                mass: -0.0,
                original_text: "-0".into(),
            })),
        );
        let expected = format!("[{prefix}+-0.0000]");
        assert!(
            pf.to_text(WriteMode::Canonical)
                .unwrap()
                .ends_with(&expected)
        );
    }
    pf.n_term_mods
        .push(modification(ModificationTag::FormulaTag(FormulaTag {
            formula_string: "H".into(),
            charge: Some(i32::MIN),
        })));
    pf.unlocalised_mods.push(UnlocalisedMod {
        modifications: vec![],
        occurrence: Some(i32::MAX),
    });
    assert!(
        pf.to_text(WriteMode::Lossless)
            .unwrap()
            .starts_with("^2147483647?[Formula:H:z-2147483648]-X")
    );
}

#[test]
fn stream_state_is_preserved_inside_chains_and_reset_between_them() {
    let mut a = peptide("A");
    attach(&mut a, 0, scored(info(""), 0.12345678));
    attach(&mut a, 0, scored(delta(1.25, ""), 0.12345678));
    attach(&mut a, 0, scored(info(""), 0.12345678));
    assert_eq!(
        a.to_text(WriteMode::Lossless).unwrap(),
        "A[#g1(0.123457)][+1.2500#g1(0.1235)][#g1(0.1235)]"
    );
    assert_eq!(
        a.to_text(WriteMode::Canonical).unwrap(),
        "A[#g1(0.12)][+1.2500#g1(0.12)][#g1(0.12)]"
    );
    let mut b = peptide("B");
    attach(&mut b, 0, scored(delta(1.25, "+1.250000"), 0.12345678));
    assert_eq!(
        b.to_text(WriteMode::Lossless).unwrap(),
        "B[+1.250000#g1(0.123457)]"
    );
    let ion = PeptidoformIon {
        chains: vec![a, b],
        ..Default::default()
    };
    assert!(
        ion.to_text(WriteMode::Lossless)
            .unwrap()
            .ends_with("//B[+1.250000#g1(0.123457)]")
    );
}

#[test]
fn default_float_six_significant_digits_rounding_and_extremes() {
    // Independent classic-locale defaultfloat expectations, not upstream tests.
    for (score, expected) in [
        (0.0, "0"),
        (-0.0, "-0"),
        (0.125, "0.125"),
        (1.23456789, "1.23457"),
        (123456.7, "123457"),
        (999999.9, "1e+06"),
        (0.00009999999, "0.0001"),
        (0.00001, "1e-05"),
        (f64::MAX, "1.79769e+308"),
        (f64::MIN_POSITIVE, "2.22507e-308"),
        (f64::from_bits(1), "4.94066e-324"),
    ] {
        let mut pf = peptide("A");
        attach(&mut pf, 0, scored(info(""), score));
        assert_eq!(
            pf.to_text(WriteMode::Lossless).unwrap(),
            format!("A[#g1({expected})]")
        );
    }
    let mut pf = peptide("A");
    attach(&mut pf, 0, modification(delta(f64::MAX, "")));
    let text = pf.to_text(WriteMode::Canonical).unwrap();
    assert_eq!(text.len(), 318);
    assert!(text.ends_with(".0000]"));
}

#[test]
fn source_label_only_labile_and_empty_container_quirks() {
    let mut pf = Peptidoform {
        name: Some(String::new()),
        ..Default::default()
    };
    let labelled = Modification {
        alternatives: vec![(info(""), Some(label(None))), (info(""), None)],
        ..Default::default()
    };
    pf.labile_mods.push(LabileModification {
        modification: labelled.clone(),
    });
    pf.n_term_mods.push(labelled);
    pf.unlocalised_mods.push(UnlocalisedMod::default());
    pf.sequence
        .push(SequenceSection::AmbiguousRegion(AmbiguousRegion::default()));
    pf.sequence
        .push(SequenceSection::ModifiedRange(ModifiedRange::default()));
    pf.c_term_mods.push(Modification::default());
    assert_eq!(
        pf.to_text(WriteMode::Lossless).unwrap(),
        "(>)?{INFO:#g1|INFO:}[#g1|INFO:]-(?)()-[]"
    );
    assert_eq!(
        Peptidoform::default()
            .to_text(WriteMode::Canonical)
            .unwrap(),
        ""
    );
    assert_eq!(
        PeptidoformIon::default()
            .to_text(WriteMode::Canonical)
            .unwrap(),
        ""
    );
}

#[test]
fn source_charge_placement_adduct_sum_and_name_omissions() {
    let mut a = peptide("A");
    a.charge = Some(ChargeState::Simple(2));
    assert_eq!(a.to_text(WriteMode::Lossless).unwrap(), "A");
    let mut ion = PeptidoformIon {
        name: Some("not emitted".into()),
        chains: vec![a],
        charge: Some(ChargeState::Simple(-3)),
        ..Default::default()
    };
    assert_eq!(ion.to_text(WriteMode::Lossless).unwrap(), "A/-3");
    ion.is_chimeric = true;
    assert_eq!(ion.to_text(WriteMode::Lossless).unwrap(), "A/2/-3");
    ion.charge = Some(ChargeState::Adducts(vec![
        AdductIon {
            formula: "Na".into(),
            charge: 1,
            occurrence: Some(2),
        },
        AdductIon {
            formula: "H".into(),
            charge: -1,
            occurrence: None,
        },
    ]));
    assert_eq!(
        ion.to_text(WriteMode::Lossless).unwrap(),
        "A/2/[Na:z+1^2,H:z-1]1+"
    );
    ion.charge = Some(ChargeState::Adducts(vec![]));
    assert_eq!(ion.to_text(WriteMode::Lossless).unwrap(), "A/2/[]0+");
    ion.charge = Some(ChargeState::Adducts(vec![AdductIon {
        charge: i32::MIN,
        ..Default::default()
    }]));
    assert_eq!(
        ion.to_text(WriteMode::Lossless).unwrap(),
        "A/2/[:z-2147483648]2147483648-"
    );
    ion.charge = Some(ChargeState::Simple(i32::MAX));
    assert_eq!(ion.to_text(WriteMode::Lossless).unwrap(), "A/2/2147483647");
}

#[test]
fn checked_consumed_numeric_and_character_state_and_atomic_errors() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut pf = peptide("A");
        attach(&mut pf, 0, modification(delta(value, "literal-δ")));
        // Lossless mode does not consume the independent numeric value.
        assert_eq!(pf.to_text(WriteMode::Lossless).unwrap(), "A[literal-δ]");
        assert!(pf.to_text(WriteMode::Canonical).is_err());
        let mut pf = peptide("B");
        attach(&mut pf, 0, scored(info(""), value));
        assert!(pf.to_text(WriteMode::Lossless).is_err());
    }
    assert!(peptide("é").to_text(WriteMode::Lossless).is_err());
    let mut pf = peptide("A");
    attach(
        &mut pf,
        0,
        modification(ModificationTag::PositionConstraint(PositionConstraint {
            residues: vec!['α'],
            ..Default::default()
        })),
    );
    assert!(pf.to_text(WriteMode::Canonical).is_err());
    // NUL is a valid source char byte: no invented grammar validation.
    assert_eq!(peptide("\0").to_text(WriteMode::Lossless).unwrap(), "\0");
    let overflow = ChargeState::Adducts(vec![
        AdductIon {
            charge: i32::MAX,
            ..Default::default()
        },
        AdductIon {
            charge: 1,
            ..Default::default()
        },
    ]);
    let mut pf = peptide("A");
    pf.charge = Some(overflow.clone());
    assert_eq!(pf.to_text(WriteMode::Lossless).unwrap(), "A");
    let mut ion = PeptidoformIon {
        chains: vec![pf],
        charge: Some(overflow),
        ..Default::default()
    };
    let saved = ion.clone();
    assert!(ion.to_text(WriteMode::Lossless).is_err());
    assert_eq!(ion, saved);
    ion.charge = Some(ChargeState::Adducts(vec![AdductIon {
        charge: i32::MIN,
        occurrence: Some(-1),
        ..Default::default()
    }]));
    assert!(ion.to_text(WriteMode::Canonical).is_err());
}

#[test]
fn shared_output_nodes_and_work_are_bounded_before_publication() {
    let mut pf = peptide("A");
    pf.name = Some("x".repeat(MAX_PROFORMA_TEXT_BYTES - 4));
    assert_eq!(
        pf.to_text(WriteMode::Lossless).unwrap().len(),
        MAX_PROFORMA_TEXT_BYTES
    );
    pf.name.as_mut().unwrap().push('x');
    assert!(pf.to_text(WriteMode::Lossless).is_err());
    let mut ion = PeptidoformIon {
        chains: vec![
            Peptidoform {
                name: Some("x".repeat(MAX_PROFORMA_TEXT_BYTES / 2)),
                ..Default::default()
            };
            2
        ],
        ..Default::default()
    };
    assert!(ion.to_text(WriteMode::Canonical).is_err());
    ion.chains.clear();
    // A one-item AST can hold a list beyond the shared node cap, even when
    // every empty string would emit only a delimiter.
    ion.charge = Some(ChargeState::Adducts(vec![
        AdductIon::default();
        MAX_PROFORMA_NODES + 1
    ]));
    assert!(
        ion.to_text(WriteMode::Canonical)
            .unwrap_err()
            .to_string()
            .contains("node")
    );
    let mut pf = peptide("A");
    let many = scored(info(""), 0.5);
    if let SequenceSection::Element(element) = &mut pf.sequence[0] {
        element.modifications = vec![many; 25_000];
    }
    assert!(
        pf.to_text(WriteMode::Lossless)
            .unwrap_err()
            .to_string()
            .contains("work")
    );
}

#[test]
fn owned_resolved_chemistry_is_shared_but_not_serialized() {
    let registry = ModificationsDB::from_records(vec![
        ResidueModification::from_record(ModificationRecord {
            name: "CallerO".into(),
            origin: Some('M'),
            diff_mono_mass: 15.994915,
            diff_formula: "O".parse().unwrap(),
            ..Default::default()
        })
        .unwrap(),
    ])
    .unwrap();
    let handle = registry
        .get_modification_handle("CallerO", Some('M'), Some(TermSpecificity::Anywhere))
        .unwrap();
    let weak = Arc::downgrade(&handle);
    let mut modification = named("arbitrary source spelling");
    modification.resolved_mod = Some(handle);
    let cloned = modification.clone();
    assert!(Arc::ptr_eq(
        modification.resolved_mod.as_ref().unwrap(),
        cloned.resolved_mod.as_ref().unwrap()
    ));
    let mut pf = peptide("M");
    attach(&mut pf, 0, modification);
    assert_eq!(
        pf.to_text(WriteMode::Canonical).unwrap(),
        "M[arbitrary source spelling]"
    );
    drop(registry);
    drop(cloned);
    assert!(weak.upgrade().is_some());
    drop(pf);
    assert!(weak.upgrade().is_none());
    let issue = ConversionIssue {
        issue_type: ConversionIssueType::UnresolvedMod,
        description: "caller diagnostic".into(),
        position: None,
    };
    assert_eq!(issue.clone(), issue);
    let group = CrossLinkGroup {
        label: "XL1".into(),
        sites: vec![(0, usize::MAX), (1, 0)],
    };
    assert_eq!(group.clone(), group);
}
