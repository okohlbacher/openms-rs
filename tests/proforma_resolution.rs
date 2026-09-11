use openms::chemistry::proforma::*;
use openms::chemistry::{
    EmpiricalFormula, ModificationProvenance, ModificationRecord, ModificationsDB,
    ResidueModification, TermSpecificity,
};
use std::sync::Arc;

fn record(
    name: &str,
    origin: Option<char>,
    term: TermSpecificity,
    formula: &str,
    mass: f64,
) -> ResidueModification {
    ResidueModification::from_record(ModificationRecord {
        name: name.into(),
        origin,
        term_specificity: term,
        diff_formula: EmpiricalFormula::parse(formula).unwrap(),
        diff_mono_mass: mass,
        ..ModificationRecord::default()
    })
    .unwrap()
}
fn modification(tag: ModificationTag) -> Modification {
    Modification {
        alternatives: vec![(tag, None)],
        resolved_mod: None,
    }
}
fn named(name: &str) -> ModificationTag {
    ModificationTag::NamedMod(NamedMod {
        name: name.into(),
        cv_hint: None,
    })
}
fn formula(value: &str) -> ModificationTag {
    ModificationTag::FormulaTag(FormulaTag {
        formula_string: value.into(),
        charge: None,
    })
}
fn info(name: &str) -> ModificationTag {
    ModificationTag::InfoTag(InfoTag { text: name.into() })
}
fn mass(value: f64) -> ModificationTag {
    ModificationTag::MassDelta(MassDelta {
        mass: value,
        ..MassDelta::default()
    })
}
fn chain(residue: char, tags: Vec<ModificationTag>) -> Peptidoform {
    Peptidoform {
        sequence: vec![SequenceSection::Element(SequenceElement {
            amino_acid: residue,
            modifications: tags.into_iter().map(modification).collect(),
        })],
        ..Peptidoform::default()
    }
}
fn mods(pf: &Peptidoform) -> &[Modification] {
    match &pf.sequence[0] {
        SequenceSection::Element(e) => &e.modifications,
        _ => panic!(),
    }
}
fn handle(pf: &Peptidoform, i: usize) -> &Arc<ResidueModification> {
    mods(pf)[i].resolved_mod.as_ref().unwrap()
}

#[test]
fn source_literal_accession_name_formula_and_residue_incompatibility() {
    // Direct source resolution assertions at ProFormaParser_test.cpp1415ff,
    // plus resolution-only portions of the formula/INFO conversion cases.
    let mut db = ModificationsDB::global().clone();
    for text in ["EM[UNIMOD:35]K", "EM[Oxidation]K"] {
        let mut pf = Peptidoform::parse(text).unwrap();
        pf.resolve_modifications(&mut db).unwrap();
        let SequenceSection::Element(m) = &pf.sequence[1] else {
            panic!()
        };
        assert_eq!(
            m.modifications[0].resolved_mod.as_ref().unwrap().full_id(),
            "Oxidation (M)"
        );
    }
    for text in [
        "PEPM[UNIMOD:447]TIDE",
        "PEPM[INFO:x|UNIMOD:447]TIDE",
        "PEPT[INFO:|Glycan:Hex]IDE",
        "pepm[Formula:O]tide",
    ] {
        let mut pf = Peptidoform::parse(text).unwrap();
        pf.resolve_modifications(&mut db).unwrap();
        let SequenceSection::Element(e) = &pf.sequence[3] else {
            panic!()
        };
        assert!(e.modifications[0].resolved_mod.is_none());
    }
    let mut pf = Peptidoform::parse("PEPM[Formula:O]TIDE").unwrap();
    pf.resolve_modifications(&mut db).unwrap();
    let SequenceSection::Element(e) = &pf.sequence[3] else {
        panic!()
    };
    let m = e.modifications[0].resolved_mod.as_ref().unwrap();
    assert_eq!(m.full_id(), "M[Formula:O1]");
    assert!(m.name().is_empty());
    assert_eq!(m.provenance(), ModificationProvenance::MassOnly);
    assert_eq!(
        m.diff_mono_mass(),
        EmpiricalFormula::parse("O").unwrap().mono_mass()
    );
}

#[test]
fn equivalent_formulas_share_handle_and_new_records_are_visible_immediately() {
    let mut db = ModificationsDB::default();
    let mut pf = chain(
        'A',
        vec![
            formula("C2H2O"),
            formula("C2H2O1"),
            formula("H2C2O1"),
            named("A[Formula:C2H2O1]"),
        ],
    );
    assert!(pf.resolve_modifications(&mut db).unwrap().is_empty());
    assert_eq!(db.len(), 1);
    for i in 1..4 {
        assert!(Arc::ptr_eq(handle(&pf, 0), handle(&pf, i)));
    }
    let retained = Arc::clone(handle(&pf, 0));
    assert!(pf.resolve_modifications(&mut db).unwrap().is_empty());
    assert_eq!(db.len(), 1);
    assert!(Arc::ptr_eq(&retained, handle(&pf, 0)));
    let mut second = chain(
        'M',
        vec![
            formula("O"),
            mass(EmpiricalFormula::parse("O").unwrap().mono_mass()),
        ],
    );
    second.resolve_modifications(&mut db).unwrap();
    assert!(Arc::ptr_eq(handle(&second, 0), handle(&second, 1)));
}

#[test]
fn isotope_canonical_order_and_site_identity_do_not_use_mass_keys() {
    let mut db = ModificationsDB::default();
    let mut pf = chain('M', vec![formula("C(13)CH"), formula("H(13)CC")]);
    pf.n_term_mods.push(modification(formula("C(13)CH")));
    pf.c_term_mods.push(modification(formula("C(13)CH")));
    pf.resolve_modifications(&mut db).unwrap();
    assert_eq!(handle(&pf, 0).full_id(), "M[Formula:(13)C1C1H1]");
    assert!(Arc::ptr_eq(handle(&pf, 0), handle(&pf, 1)));
    assert_eq!(
        pf.n_term_mods[0].resolved_mod.as_ref().unwrap().full_id(),
        ".n[Formula:(13)C1C1H1]"
    );
    assert_eq!(
        pf.c_term_mods[0].resolved_mod.as_ref().unwrap().full_id(),
        ".c[Formula:(13)C1C1H1]"
    );
    assert_eq!(db.len(), 3);
}

#[test]
fn formula_absolute_masses_preserve_free_residue_zero_and_terminal_defaults() {
    let delta = EmpiricalFormula::parse("H-1").unwrap();
    let mut db = ModificationsDB::default();
    for aa in ['B', 'Z', 'X'] {
        let mut pf = chain(aa, vec![formula("H-1")]);
        pf.resolve_modifications(&mut db).unwrap();
        assert_eq!(handle(&pf, 0).mono_mass(), delta.mono_mass());
        assert_eq!(handle(&pf, 0).average_mass(), delta.average_mass());
        assert!(handle(&pf, 0).full_name().starts_with("[-"));
    }
    let mut pf = chain('A', vec![formula("O")]);
    pf.n_term_mods.push(modification(formula("O")));
    pf.c_term_mods.push(modification(formula("O")));
    pf.resolve_modifications(&mut db).unwrap();
    let o = EmpiricalFormula::parse("O").unwrap();
    let a = EmpiricalFormula::parse("C3H7NO2").unwrap();
    assert_eq!(handle(&pf, 0).mono_mass(), o.mono_mass() + a.mono_mass());
    assert_eq!(
        handle(&pf, 0).average_mass(),
        o.average_mass() + a.average_mass()
    );
    for (modification, base) in [(&pf.n_term_mods[0], "H"), (&pf.c_term_mods[0], "OH")] {
        let m = modification.resolved_mod.as_ref().unwrap();
        assert_eq!(
            m.mono_mass(),
            o.mono_mass() + EmpiricalFormula::parse(base).unwrap().mono_mass()
        );
        assert_eq!(m.average_mass(), 0.);
        assert_eq!(m.diff_average_mass(), o.average_mass());
    }
}

#[test]
fn unusable_formula_states_are_unresolved_without_swallowing_resource_errors() {
    let mut db = ModificationsDB::default();
    let mut pf = chain(
        'M',
        vec![
            formula(""),
            formula("C0"),
            formula("[13C2]"),
            formula("bad"),
            formula("H1+"),
            formula("C2147483648"),
            ModificationTag::FormulaTag(FormulaTag {
                formula_string: "Zn1".into(),
                charge: Some(2),
            }),
        ],
    );
    pf.resolve_modifications(&mut db).unwrap();
    assert!(mods(&pf).iter().all(|m| m.resolved_mod.is_none()));
    assert_eq!(db.len(), 0);
    let mut no_residue = Peptidoform {
        n_term_mods: vec![modification(formula("O"))],
        labile_mods: vec![LabileModification {
            modification: modification(formula("O")),
        }],
        ..Peptidoform::default()
    };
    no_residue.resolve_modifications(&mut db).unwrap();
    assert!(no_residue.n_term_mods[0].resolved_mod.is_some());
    assert!(
        no_residue.labile_mods[0]
            .modification
            .resolved_mod
            .is_none()
    );
}

#[test]
fn source_accession_anywhere_preference_differs_from_unrestricted_named_search() {
    let build = |full: &str, term| {
        ResidueModification::from_record(ModificationRecord {
            name: "Shared".into(),
            full_id: full.into(),
            record_id: Some(17),
            origin: Some('M'),
            term_specificity: term,
            ..ModificationRecord::default()
        })
        .unwrap()
    };
    let mut db = ModificationsDB::from_records(vec![
        build("first terminal", TermSpecificity::NTerm),
        build("first anywhere", TermSpecificity::Anywhere),
        build("second anywhere", TermSpecificity::Anywhere),
    ])
    .unwrap();
    assert!(
        db.get_modification_handle("Shared", Some('M'), None)
            .is_err()
    );
    let mut pf = chain(
        'M',
        vec![
            named("Shared"),
            ModificationTag::CvAccession(CvAccession {
                database: CvDatabase::Unimod,
                accession: "17".into(),
            }),
        ],
    );
    let warnings = pf.resolve_modifications(&mut db).unwrap();
    assert_eq!(handle(&pf, 0).full_id(), "first terminal");
    assert_eq!(handle(&pf, 1).full_id(), "first anywhere");
    assert!(matches!(
        &warnings[..],
        [ResolutionWarning::AmbiguousAccession { .. }]
    ));
}

#[test]
fn every_cv_prefix_exact_alias_precedence_and_ignored_named_hint() {
    let kinds = [
        (CvDatabase::Unimod, "UNIMOD"),
        (CvDatabase::Mod, "MOD"),
        (CvDatabase::Resid, "RESID"),
        (CvDatabase::Xlmod, "XLMOD"),
        (CvDatabase::Gno, "GNO"),
    ];
    let mut entries: Vec<_> = kinds
        .iter()
        .map(|(_, name)| {
            record(
                &format!("{name}:7"),
                Some('M'),
                TermSpecificity::Anywhere,
                "",
                1.,
            )
        })
        .collect();
    entries.push(
        ResidueModification::from_record(ModificationRecord {
            name: "normalized".into(),
            record_id: Some(7),
            origin: Some('M'),
            ..ModificationRecord::default()
        })
        .unwrap(),
    );
    let mut db = ModificationsDB::from_records(entries).unwrap();
    let mut tags: Vec<_> = kinds
        .iter()
        .map(|(database, _)| {
            ModificationTag::CvAccession(CvAccession {
                database: *database,
                accession: "7".into(),
            })
        })
        .collect();
    tags.push(named("uNiMoD:7"));
    tags.push(ModificationTag::NamedMod(NamedMod {
        cv_hint: Some(CvDatabase::Gno),
        name: "MOD:7".into(),
    }));
    let mut pf = chain('M', tags);
    pf.resolve_modifications(&mut db).unwrap();
    for (i, (_, prefix)) in kinds.iter().enumerate() {
        assert_eq!(handle(&pf, i).name(), format!("{prefix}:7"));
    }
    assert_eq!(handle(&pf, 5).name(), "normalized");
    assert_eq!(handle(&pf, 6).name(), "MOD:7");
}

#[test]
fn source_wildcard_queries_anonymous_x_and_empty_aliases() {
    let a = record("firstA", Some('A'), TermSpecificity::Anywhere, "", 1.);
    let anonymous = ResidueModification::from_record(ModificationRecord {
        full_id: "X[123]".into(),
        origin: Some('X'),
        diff_mono_mass: 123.,
        ..ModificationRecord::default()
    })
    .unwrap();
    let wildcard = record("namedX", Some('X'), TermSpecificity::Anywhere, "", 2.);
    let mut db = ModificationsDB::from_records(vec![a, anonymous, wildcard]).unwrap();
    for aa in ['X', '.', '?', '\0'] {
        let mut pf = chain(aa, vec![named("firstA")]);
        pf.resolve_modifications(&mut db).unwrap();
        assert_eq!(handle(&pf, 0).name(), "firstA");
    }
    let mut pf = chain('M', vec![named("X[123]"), named("namedX"), named("")]);
    pf.resolve_modifications(&mut db).unwrap();
    assert!(mods(&pf)[0].resolved_mod.is_none());
    assert_eq!(handle(&pf, 1).name(), "namedX");
    assert_eq!(handle(&pf, 2).name(), "namedX");
    let mut x = chain('X', vec![named("X[123]")]);
    x.resolve_modifications(&mut db).unwrap();
    assert_eq!(handle(&x, 0).full_id(), "X[123]");
}

#[test]
fn strict_mass_window_zero_filter_and_provider_ties() {
    let mut db = ModificationsDB::from_records(vec![
        record("zero", Some('M'), TermSpecificity::Anywhere, "", 0.),
        record(
            "low",
            Some('M'),
            TermSpecificity::Anywhere,
            "",
            2. - 0.00390625,
        ),
        record(
            "high",
            Some('M'),
            TermSpecificity::Anywhere,
            "",
            2. + 0.00390625,
        ),
    ])
    .unwrap();
    let mut pf = chain('M', vec![mass(0.), mass(0.001), mass(2.), mass(1e50)]);
    pf.resolve_modifications(&mut db).unwrap();
    assert_eq!(handle(&pf, 0).name(), "zero");
    assert!(mods(&pf)[1].resolved_mod.is_none());
    assert_eq!(handle(&pf, 2).name(), "low");
    assert!(mods(&pf)[3].resolved_mod.is_none());
    let mut db = ModificationsDB::from_records(vec![record(
        "edge",
        Some('M'),
        TermSpecificity::Anywhere,
        "",
        0.01,
    )])
    .unwrap();
    let mut pf = chain('M', vec![mass(0.)]);
    pf.resolve_modifications(&mut db).unwrap();
    assert!(mods(&pf)[0].resolved_mod.is_none());
}

#[test]
fn first_chemistry_only_and_empty_alternatives_keep_source_behavior() {
    let mut db = ModificationsDB::from_records(vec![record(
        "known",
        Some('M'),
        TermSpecificity::Anywhere,
        "O",
        16.,
    )])
    .unwrap();
    let retained = Arc::clone(&db.entries()[0]);
    let mut pf = chain('M', vec![]);
    let SequenceSection::Element(e) = &mut pf.sequence[0] else {
        panic!()
    };
    e.modifications.push(Modification {
        alternatives: vec![
            (info("note"), None),
            (named("absent"), None),
            (named("known"), None),
        ],
        resolved_mod: Some(Arc::clone(&retained)),
    });
    e.modifications.push(Modification {
        alternatives: vec![],
        resolved_mod: Some(Arc::clone(&retained)),
    });
    e.modifications.push(Modification {
        alternatives: vec![
            (info(""), None),
            (
                ModificationTag::GlycanComposition(GlycanComposition::default()),
                None,
            ),
            (named("known"), None),
        ],
        resolved_mod: None,
    });
    let warnings = pf.resolve_modifications(&mut db).unwrap();
    assert!(
        matches!(&warnings[..],[ResolutionWarning::ModificationNotFound {name}] if name=="absent")
    );
    assert!(mods(&pf)[0].resolved_mod.is_none());
    assert!(Arc::ptr_eq(handle(&pf, 1), &retained));
    assert!(mods(&pf)[2].resolved_mod.is_none());
}

#[test]
fn defined_info_identity_wins_and_inline_interning_still_occurs() {
    let definition = record(
        "Test4a:Disagree",
        Some('K'),
        TermSpecificity::Anywhere,
        "C2H2O",
        42.,
    );
    let mut db = ModificationsDB::from_records(vec![definition]).unwrap();
    let original = Arc::clone(&db.entries()[0]);
    let mut pf = Peptidoform::parse("PEPK[Formula:O|INFO:Test4a:Disagree]TIDE").unwrap();
    let warnings = pf.resolve_modifications(&mut db).unwrap();
    assert_eq!(db.len(), 2);
    let SequenceSection::Element(k) = &pf.sequence[3] else {
        panic!()
    };
    assert!(Arc::ptr_eq(
        k.modifications[0].resolved_mod.as_ref().unwrap(),
        &original
    ));
    assert!(
        matches!(&warnings[..],[ResolutionWarning::DefinitionDisagreement {definition,inline}] if definition.name()=="Test4a:Disagree" && inline.full_id()=="K[Formula:O1]")
    );
    let mut same = Peptidoform::parse("PEPK[INFO:Test4a:Disagree|Formula:C2H2O]TIDE").unwrap();
    assert!(same.resolve_modifications(&mut db).unwrap().is_empty());
    // Formula identity takes precedence even over a deliberately different mass.
    assert_ne!(
        original.diff_mono_mass(),
        EmpiricalFormula::parse("C2H2O").unwrap().mono_mass()
    );
}

#[test]
fn vocabulary_info_does_not_override_and_mass_agreement_uses_one_microdalton() {
    let cv = record("CvName", Some('M'), TermSpecificity::Anywhere, "O", 16.)
        .with_provenance(ModificationProvenance::Cv);
    let mut db = ModificationsDB::from_records(vec![cv]).unwrap();
    let mut pf = Peptidoform::parse("M[INFO:CvName|Formula:C2H2O]").unwrap();
    assert!(pf.resolve_modifications(&mut db).unwrap().is_empty());
    assert_eq!(handle(&pf, 0).full_id(), "M[Formula:C2H2O1]");
    for (offset, warnings) in [(1e-6, 0), (1.000001e-6, 1)] {
        let mut db = ModificationsDB::from_records(vec![
            record("definition", Some('M'), TermSpecificity::Anywhere, "", 0.),
            record("inline", Some('M'), TermSpecificity::Anywhere, "", offset),
        ])
        .unwrap();
        let mut pf = Peptidoform::parse("M[inline|INFO:definition]").unwrap();
        assert_eq!(pf.resolve_modifications(&mut db).unwrap().len(), warnings);
        assert_eq!(handle(&pf, 0).name(), "definition");
    }
}

#[test]
fn every_source_group_is_visited_but_range_element_annotations_are_untouched() {
    let mut db = ModificationsDB::from_records(vec![
        record("all", Some('M'), TermSpecificity::Anywhere, "", 1.),
        record("all", None, TermSpecificity::NTerm, "", 2.),
        record("all", None, TermSpecificity::CTerm, "", 3.),
    ])
    .unwrap();
    let keep = Arc::clone(&db.entries()[0]);
    let element = || SequenceElement {
        amino_acid: 'M',
        modifications: vec![modification(named("all"))],
    };
    let mut skipped = element();
    skipped.modifications[0].alternatives = vec![(mass(f64::NAN), None)];
    skipped.modifications[0].resolved_mod = Some(Arc::clone(&keep));
    let mut pf = Peptidoform {
        sequence: vec![
            SequenceSection::Element(element()),
            SequenceSection::AmbiguousRegion(AmbiguousRegion {
                elements: vec![element()],
            }),
            SequenceSection::ModifiedRange(ModifiedRange {
                elements: vec![skipped],
                modifications: vec![modification(named("all"))],
            }),
        ],
        n_term_mods: vec![modification(named("all"))],
        c_term_mods: vec![modification(named("all"))],
        unlocalised_mods: vec![UnlocalisedMod {
            modifications: vec![modification(named("all"))],
            occurrence: Some(-42),
        }],
        labile_mods: vec![LabileModification {
            modification: modification(named("all")),
        }],
        global_mods: vec![
            GlobalModEntry::GlobalModification(GlobalModification {
                modification: modification(named("all")),
                locations: vec!["uninterpreted".into()],
            }),
            GlobalModEntry::IsotopeReplacement(IsotopeReplacement {
                isotope: "uninterpreted".into(),
            }),
        ],
        ..Peptidoform::default()
    };
    assert!(pf.resolve_modifications(&mut db).unwrap().is_empty());
    assert!(Arc::ptr_eq(handle(&pf, 0), &keep));
    let SequenceSection::AmbiguousRegion(region) = &pf.sequence[1] else {
        panic!()
    };
    assert!(region.elements[0].modifications[0].resolved_mod.is_some());
    let SequenceSection::ModifiedRange(range) = &pf.sequence[2] else {
        panic!()
    };
    assert!(range.modifications[0].resolved_mod.is_some());
    assert!(Arc::ptr_eq(
        range.elements[0].modifications[0]
            .resolved_mod
            .as_ref()
            .unwrap(),
        &keep
    ));
    assert_eq!(
        pf.n_term_mods[0]
            .resolved_mod
            .as_ref()
            .unwrap()
            .diff_mono_mass(),
        2.
    );
    assert_eq!(
        pf.c_term_mods[0]
            .resolved_mod
            .as_ref()
            .unwrap()
            .diff_mono_mass(),
        3.
    );
    assert!(
        pf.unlocalised_mods[0].modifications[0]
            .resolved_mod
            .is_some()
    );
    assert!(pf.labile_mods[0].modification.resolved_mod.is_some());
    let GlobalModEntry::GlobalModification(g) = &pf.global_mods[0] else {
        panic!()
    };
    assert!(g.modification.resolved_mod.is_some());
}

#[test]
fn late_invalid_mass_and_large_lookup_roll_back_registry_and_all_handles() {
    for late in [mass(f64::INFINITY), named(&"x".repeat(60_000))] {
        let mut db = ModificationsDB::default();
        let mut pf = chain('M', vec![formula("O"), late]);
        let text = pf.clone();
        assert!(pf.resolve_modifications(&mut db).is_err());
        assert_eq!(pf, text);
        assert!(db.is_empty());
    }
}

#[test]
fn existing_formula_key_wins_without_chemistry_or_site_revalidation() {
    let prior = ResidueModification::from_record(ModificationRecord {
        name: "prior".into(),
        full_id: "M[Formula:O1]".into(),
        origin: Some('K'),
        diff_formula: EmpiricalFormula::parse("H2").unwrap(),
        diff_mono_mass: 123.,
        ..ModificationRecord::default()
    })
    .unwrap();
    let mut db = ModificationsDB::from_records(vec![prior]).unwrap();
    let original = Arc::clone(&db.entries()[0]);
    let mut pf = chain('M', vec![formula("O")]);
    pf.resolve_modifications(&mut db).unwrap();
    assert!(Arc::ptr_eq(handle(&pf, 0), &original));
    assert_eq!(db.len(), 1);
    assert_eq!(handle(&pf, 0).origin(), Some('K'));
}

#[test]
fn defined_gate_handles_anonymous_records_unicode_and_exact_names() {
    let anonymous = ResidueModification::from_record(ModificationRecord {
        full_id: "defined anonymous".into(),
        origin: Some('M'),
        diff_mono_mass: 100.,
        ..ModificationRecord::default()
    })
    .unwrap();
    let unicode = record("化学", Some('M'), TermSpecificity::Anywhere, "O", 16.);
    let mut db = ModificationsDB::from_records(vec![anonymous, unicode]).unwrap();
    let mut pf =
        Peptidoform::parse("M[Formula:O|INFO:defined anonymous][化学][неизвестно]").unwrap();
    let warnings = pf.resolve_modifications(&mut db).unwrap();
    assert_eq!(handle(&pf, 0).full_id(), "defined anonymous");
    assert_eq!(handle(&pf, 1).name(), "化学");
    assert!(
        matches!(&warnings[..], [ResolutionWarning::DefinitionDisagreement {..}, ResolutionWarning::ModificationNotFound {name}] if name == "неизвестно")
    );
    let mut scalar = chain('🧪', vec![formula("O")]);
    scalar.resolve_modifications(&mut db).unwrap();
    assert!(mods(&scalar)[0].resolved_mod.is_none());
}
