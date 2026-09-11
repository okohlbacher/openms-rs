#![allow(clippy::field_reassign_with_default)]
use openms::chemistry::proforma::*;
use openms::chemistry::{
    AASequence, EmpiricalFormula, ModificationRecord, ModificationsDB, PROTON_MASS_U,
    ProteinProteinCrossLink, ResidueModification, TheoreticalIonSeries,
    TheoreticalSpectrumGenerator, TheoreticalSpectrumGeneratorXLMS,
};
use openms::{MSSpectrum, kernel::SpectrumType};
use std::sync::Arc;
fn empty_db() -> ModificationsDB {
    ModificationsDB::from_records(vec![]).unwrap()
}
fn global() -> ModificationsDB {
    ModificationsDB::global().clone()
}
fn pf(text: &str) -> Peptidoform {
    Peptidoform::parse(text).unwrap()
}
fn ion(text: &str) -> PeptidoformIon {
    PeptidoformIon::parse(text).unwrap()
}
fn options(types: &str) -> SpectrumGenerationOptions {
    SpectrumGenerationOptions {
        ion_types: types.into(),
        ..Default::default()
    }
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-8, "{a} != {b}")
}
fn named(name: &str) -> ModificationTag {
    ModificationTag::NamedMod(NamedMod {
        name: name.into(),
        cv_hint: None,
    })
}
fn delta(mass: f64) -> ModificationTag {
    ModificationTag::MassDelta(MassDelta {
        mass,
        ..Default::default()
    })
}
fn linked(
    tag: ModificationTag,
    id: &str,
    record: Option<Arc<ResidueModification>>,
) -> Modification {
    Modification {
        alternatives: vec![(
            tag,
            Some(Label {
                label_type: LabelType::Crosslink,
                identifier: id.into(),
                score: None,
            }),
        )],
        resolved_mod: record,
    }
}
fn record(name: &str, mass: f64) -> Arc<ResidueModification> {
    Arc::new(
        ResidueModification::from_record(ModificationRecord {
            name: name.into(),
            full_id: format!("{name} (M)"),
            origin: Some('M'),
            diff_mono_mass: mass,
            ..Default::default()
        })
        .unwrap(),
    )
}
fn two(a: Modification, b: Modification) -> PeptidoformIon {
    let mut a_pf = pf("AM");
    let mut b_pf = pf("MA");
    if let SequenceSection::Element(e) = &mut a_pf.sequence[1] {
        e.modifications.push(a)
    }
    if let SequenceSection::Element(e) = &mut b_pf.sequence[0] {
        e.modifications.push(b)
    }
    PeptidoformIon {
        chains: vec![a_pf, b_pf],
        ..Default::default()
    }
}
fn label() -> Modification {
    linked(ModificationTag::InfoTag(InfoTag::default()), "XL1", None)
}
fn precursor(s: &MSSpectrum) -> f64 {
    let i = s.string_data_arrays[0]
        .data
        .iter()
        .position(|n| n == "[M+H]")
        .unwrap();
    s.peaks[i].mz
}
fn source_full(text: &str) -> f64 {
    let mut mass = 0.;
    for c in text.chars() {
        let formula = match c {
            'A' => "C3H7NO2",
            'M' => "C5H11NO2S",
            'G' => "C2H5NO2",
            _ => panic!(),
        };
        mass += EmpiricalFormula::parse(formula).unwrap().mono_mass()
            - EmpiricalFormula::parse("H2O").unwrap().mono_mass();
    }
    mass + EmpiricalFormula::parse("H2O").unwrap().mono_mass()
}

#[test]
fn source_ten_sections_keep_their_original_assertion_strength() {
    let mut db = global();
    let mut p = pf("PEPTIDE");
    p.resolve_modifications(&mut db).unwrap();
    assert!(p.can_generate_spectrum(&mut db).unwrap().value);
    let mut bad = pf("PEM[UnknownMod999]PTIDE");
    bad.resolve_modifications(&mut db).unwrap();
    assert!(!bad.can_generate_spectrum(&mut db).unwrap().value);
    assert!(
        !bad.spectrum_generation_issues(&mut db)
            .unwrap()
            .value
            .is_empty()
    );
    let mut one = ion("PEPTIDE/2");
    one.chains[0].resolve_modifications(&mut db).unwrap();
    assert!(one.can_generate_spectrum(&mut db).unwrap().value);
    assert!(
        one.spectrum_generation_issues(&mut db)
            .unwrap()
            .value
            .is_empty()
    );
    let mut modified = pf("PEM[UNIMOD:35]PTIDE");
    modified.resolve_modifications(&mut db).unwrap();
    assert!(
        modified
            .spectrum_generation_issues(&mut db)
            .unwrap()
            .value
            .is_empty()
    );
    let s = p
        .generate_spectrum(&SpectrumGenerationOptions::default(), &mut db)
        .unwrap()
        .value;
    assert!(s.len() >= 10);
    let bym = p.generate_spectrum(&options("byM"), &mut db).unwrap().value;
    assert!(bym.len() >= s.len());
    assert!(
        !p.generate_spectrum(&options("abcxyz"), &mut db)
            .unwrap()
            .value
            .is_empty()
    );
    assert!(
        !one.generate_spectrum(&options("by"), &mut db)
            .unwrap()
            .value
            .is_empty()
    );
    let mut xl = ion("PEPK[+138.068#XL1]IDE//ANOK[#XL1]THER");
    for chain in &mut xl.chains {
        chain.resolve_modifications(&mut db).unwrap();
    }
    assert!(xl.can_generate_spectrum(&mut db).unwrap().value);
    assert!(
        xl.spectrum_generation_issues(&mut db)
            .unwrap()
            .value
            .is_empty()
    );
    assert!(
        !xl.generate_spectrum(&options("by"), &mut db)
            .unwrap()
            .value
            .is_empty()
    );
    let missing = ion("PEPTIDE//ANOTHER");
    assert!(!missing.can_generate_spectrum(&mut db).unwrap().value);
    assert!(
        !missing
            .spectrum_generation_issues(&mut db)
            .unwrap()
            .value
            .is_empty()
    );
    let chim = ion("PEPTIDE+ANOTHER");
    assert!(!chim.can_generate_spectrum(&mut db).unwrap().value);
}
#[test]
fn every_ordinary_flag_matches_the_real_backend_and_source_metadata() {
    let p = pf("PEPTIDE");
    let sequence = AASequence::parse("PEPTIDE").unwrap();
    for types in [
        "a", "b", "c", "x", "y", "z", "M", "I", "abcxyzMI", "", "BY", " bby?! ",
    ] {
        let mut g = TheoreticalSpectrumGenerator::default();
        g.ion_series = [
            ('a', TheoreticalIonSeries::A),
            ('b', TheoreticalIonSeries::B),
            ('c', TheoreticalIonSeries::C),
            ('x', TheoreticalIonSeries::X),
            ('y', TheoreticalIonSeries::Y),
            ('z', TheoreticalIonSeries::Z),
        ]
        .into_iter()
        .filter_map(|(c, s)| types.contains(c).then_some(s))
        .collect();
        g.add_metainfo = true;
        g.add_precursor_peaks = types.contains('M');
        g.add_abundant_immonium_ions = types.contains('I');
        let actual = p
            .generate_spectrum(&options(types), &mut empty_db())
            .unwrap()
            .value;
        assert_eq!(
            actual,
            g.generate(&sequence, 1, 1, None).unwrap(),
            "{types}"
        );
        assert_eq!(actual.ms_level, 2);
        assert_eq!(actual.spectrum_type, SpectrumType::Centroid);
        assert_eq!(actual.precursors[0].charge, 2);
    }
    let s = p
        .generate_spectrum(&options("by"), &mut empty_db())
        .unwrap()
        .value;
    assert_eq!(s.len(), 11);
    assert!(s.string_data_arrays[0].data.contains(&"b2+".into()));
    assert!(!s.string_data_arrays[0].data.contains(&"b1+".into()));
}
#[test]
fn ordinary_losses_metadata_and_modified_chemistry_are_forwarded() {
    for text in [
        "PEPTIDE",
        "[UNIMOD:1]-PEM[UNIMOD:35]PTIDE-[UNIMOD:2]",
        "AM[Formula:O]A",
    ] {
        for loss in [false, true] {
            for meta in [false, true] {
                let mut db = global();
                let p = pf(text);
                let seq = p.to_aa_sequence(&mut db).unwrap().value;
                let opt = SpectrumGenerationOptions {
                    min_charge: 1,
                    max_charge: 2,
                    add_losses: loss,
                    add_metainfo: meta,
                    ..Default::default()
                };
                let g = TheoreticalSpectrumGenerator {
                    add_losses: loss,
                    add_metainfo: meta,
                    ..Default::default()
                };
                assert_eq!(
                    p.generate_spectrum(&opt, &mut db).unwrap().value,
                    g.generate(&seq, 1, 2, None).unwrap()
                );
            }
        }
    }
}
#[test]
fn all_ion_issue_branches_have_exact_source_order_and_position() {
    let mut cases = vec![
        (PeptidoformIon::default(), "No peptide chains to fragment"),
        (
            PeptidoformIon {
                chains: vec![pf("A")],
                is_chimeric: true,
                ..Default::default()
            },
            "Theoretical spectrum generation not supported for chimeric spectra.",
        ),
        (
            PeptidoformIon {
                chains: vec![pf("A"); 3],
                ..Default::default()
            },
            "Only two-chain cross-links are currently supported for spectrum generation",
        ),
        (ion("AA//AA"), "Cross-link label not found in both chains"),
        (
            ion("AM[#XL1]//M[#XL2]A"),
            "Cross-link labels don't match between chains",
        ),
    ];
    cases[0].0.is_chimeric = true;
    for (p, text) in cases {
        let result = p.spectrum_generation_issues(&mut empty_db()).unwrap().value;
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].position, Some(0));
        assert_eq!(
            result[0].issue_type,
            ConversionIssueType::UnsupportedFeature
        );
        assert_eq!(result[0].description, text);
        let error = p
            .generate_spectrum(&options("by"), &mut empty_db())
            .unwrap_err()
            .to_string();
        assert!(error.contains(&format!("Spectrum generation failed: {text}; ")));
    }
}
#[test]
fn advisory_success_does_not_hide_strict_or_generator_failures() {
    let p = pf("M[INFO:note]");
    assert!(p.can_generate_spectrum(&mut empty_db()).unwrap().value);
    assert!(
        p.generate_spectrum(&options("by"), &mut empty_db())
            .unwrap_err()
            .to_string()
            .contains("Unresolved modification at position 0")
    );
    let p = pf("AX");
    assert!(p.can_generate_spectrum(&mut empty_db()).unwrap().value);
    assert!(
        p.generate_spectrum(&options("by"), &mut empty_db())
            .is_err()
    );
    let p = Peptidoform::default();
    assert_eq!(
        p.generate_spectrum(&options("by"), &mut empty_db())
            .unwrap()
            .value,
        MSSpectrum::default()
    );
    assert!(
        pf("A")
            .generate_spectrum(&options("cx"), &mut empty_db())
            .is_err()
    );
    for charge in [0, -1, 256, i32::MAX] {
        let mut o = options("by");
        o.min_charge = charge;
        o.max_charge = charge;
        assert!(pf("AA").generate_spectrum(&o, &mut empty_db()).is_err());
    }
    let mut o = options("y");
    o.min_charge = 255;
    o.max_charge = 255;
    let s = pf("AA")
        .generate_spectrum(&o, &mut empty_db())
        .unwrap()
        .value;
    assert_eq!(s.precursors[0].charge, 256);
}
#[test]
fn ignored_chain_and_ion_context_cannot_change_spectra() {
    let mut p = pf("PEPTIDE");
    let expected = p
        .generate_spectrum(&options("by"), &mut empty_db())
        .unwrap()
        .value;
    p.name = Some("name".into());
    p.charge = Some(ChargeState::Adducts(vec![AdductIon {
        formula: "not even a formula".into(),
        charge: i32::MIN,
        occurrence: Some(i32::MAX),
    }]));
    assert_eq!(
        p.generate_spectrum(&options("by"), &mut empty_db())
            .unwrap()
            .value,
        expected
    );
    let value = PeptidoformIon {
        name: Some("ion".into()),
        chains: vec![p],
        charge: Some(ChargeState::Simple(-99)),
        is_chimeric: false,
    };
    assert_eq!(
        value
            .generate_spectrum(&options("by"), &mut empty_db())
            .unwrap()
            .value,
        expected
    );
}
#[test]
fn source_repeated_issue_pass_can_use_a_new_formula_definition() {
    let mut p = pf("MM[Formula:Cl101]");
    if let SequenceSection::Element(e) = &mut p.sequence[0] {
        e.modifications.push(Modification {
            alternatives: vec![(named("M[Formula:Cl101]"), None)],
            resolved_mod: None,
        });
    }
    let before = p.clone();
    let mut db = empty_db();
    let report = p.spectrum_generation_issues(&mut db).unwrap();
    assert!(report.value.is_empty());
    assert_eq!(db.len(), 1);
    assert_eq!(report.warnings.len(), 1);
    assert_eq!(p, before);
    let mut db = empty_db();
    let s = PeptidoformIon {
        chains: vec![p.clone()],
        ..Default::default()
    }
    .generate_spectrum(&options("y"), &mut db)
    .unwrap();
    assert_eq!(db.len(), 1);
    assert_eq!(s.warnings.len(), 1);
    assert_eq!(p, before);
}
#[test]
fn cpp038_linker_mass_is_counted_twice_with_one_attached_endpoint() {
    let d = record("D", 100.);
    let mut db = ModificationsDB::from_records(vec![(*d).clone()]).unwrap();
    let p = two(linked(delta(100.), "XL1", Some(d)), label());
    let before = p.clone();
    let o = SpectrumGenerationOptions {
        min_charge: 2,
        max_charge: 2,
        ..options("M")
    };
    let s = p.generate_spectrum(&o, &mut db).unwrap().value;
    let base = source_full("AM") + source_full("MA");
    near(precursor(&s), (base + 200. + 2. * PROTON_MASS_U) / 2.);
    let chemical = (base + 100. + 2. * PROTON_MASS_U) / 2.;
    near(precursor(&s) - chemical, 50.);
    near(p.mono_mass(&mut db).unwrap().value, base + 100.);
    assert_eq!(p, before);
    assert_eq!(
        s.string_data_arrays[0]
            .data
            .iter()
            .filter(|n| n.as_str() == "[M+H]")
            .count(),
        2
    );
    assert_eq!(s.ms_level, MSSpectrum::default().ms_level);
    assert!(s.precursors.is_empty());
}
#[test]
fn cpp038_two_attached_endpoints_make_three_source_copies() {
    let d = record("D", 100.);
    let mut db = ModificationsDB::from_records(vec![(*d).clone()]).unwrap();
    let p = two(
        linked(delta(100.), "XL1", Some(d.clone())),
        linked(delta(100.), "XL1", Some(d)),
    );
    let s = p.generate_spectrum(&options("M"), &mut db).unwrap().value;
    let base = source_full("AM") + source_full("MA");
    near(precursor(&s), base + 300. + PROTON_MASS_U);
    near(p.mono_mass(&mut db).unwrap().value, base + 100.);
}
#[test]
fn original_resolution_timing_changes_separate_named_linker_mass() {
    let d = record("D", 100.);
    let mut db = ModificationsDB::from_records(vec![(*d).clone()]).unwrap();
    let unresolved = two(linked(named("D"), "XL1", None), label());
    let resolved = two(linked(named("D"), "XL1", Some(d)), label());
    assert!(unresolved.can_generate_spectrum(&mut db).unwrap().value);
    let a = unresolved
        .generate_spectrum(&options("M"), &mut db)
        .unwrap()
        .value;
    let b = resolved
        .generate_spectrum(&options("M"), &mut db)
        .unwrap()
        .value;
    near(precursor(&b) - precursor(&a), 100.);
}
#[test]
fn raw_mass_wins_over_handle_and_strict_threshold_chooses_beta() {
    for a in [-0.01, 0., 0.0009, 0.001, 0.0011] {
        let ar = record("A", a + 0.00001);
        let br = record("B", 7.);
        let mut db = ModificationsDB::from_records(vec![(*ar).clone(), (*br).clone()]).unwrap();
        let p = two(
            linked(delta(a), "XL1", Some(ar)),
            linked(delta(7.), "XL1", Some(br)),
        );
        let s = p.generate_spectrum(&options("M"), &mut db).unwrap().value;
        let link = if a > 0.001 { a } else { 7. };
        near(
            precursor(&s),
            source_full("AM") + source_full("MA") + (a + 0.00001) + 7. + link + PROTON_MASS_U,
        );
    }
}
#[test]
fn only_first_alternative_and_ordinary_residue_labels_are_found() {
    let mut p = two(label(), label());
    if let SequenceSection::Element(e) = &mut p.chains[0].sequence[1] {
        e.modifications[0]
            .alternatives
            .insert(0, (ModificationTag::InfoTag(InfoTag::default()), None));
    }
    assert!(!p.can_generate_spectrum(&mut empty_db()).unwrap().value);
    let mut p = two(label(), label());
    if let SequenceSection::Element(e) = &mut p.chains[0].sequence[1] {
        p.chains[0].n_term_mods = std::mem::take(&mut e.modifications);
    }
    assert!(!p.can_generate_spectrum(&mut empty_db()).unwrap().value);
    let mut p = two(label(), label());
    if let SequenceSection::Element(e) = &mut p.chains[0].sequence[1] {
        e.modifications.push(linked(delta(999.), "XL-other", None));
    }
    assert!(p.can_generate_spectrum(&mut empty_db()).unwrap().value);
}
#[test]
fn cpp047_range_cursor_uses_source_position_instead_of_flattened_position() {
    let mut p = two(label(), label());
    p.chains[0].sequence[0] = SequenceSection::ModifiedRange(ModifiedRange {
        elements: vec![
            SequenceElement {
                amino_acid: 'A',
                modifications: vec![],
            },
            SequenceElement {
                amino_acid: 'G',
                modifications: vec![],
            },
        ],
        ..Default::default()
    });
    let actual = p
        .generate_spectrum(&options("by"), &mut empty_db())
        .unwrap()
        .value;
    let mut link = ProteinProteinCrossLink::default();
    link.alpha = Some(Arc::new(AASequence::parse("AGM").unwrap()));
    link.beta = Some(Arc::new(AASequence::parse("MA").unwrap()));
    link.cross_link_position = (0, 0);
    let mut g = TheoreticalSpectrumGeneratorXLMS::default();
    g.options.add_a_ions = false;
    g.options.add_precursor_peaks = false;
    let mut expected = MSSpectrum::default();
    for alpha in [true, false] {
        g.get_crosslink_ion_spectrum(&mut expected, &link, alpha, 1, 1)
            .unwrap();
    }
    assert_eq!(actual, expected);
    link.cross_link_position.0 = 2;
    let mut corrected = MSSpectrum::default();
    for alpha in [true, false] {
        g.get_crosslink_ion_spectrum(&mut corrected, &link, alpha, 1, 1)
            .unwrap();
    }
    assert_ne!(actual, corrected);
}
#[test]
fn xlms_empty_ion_flags_keep_k_linked_and_charges_while_i_is_ignored() {
    let p = two(label(), label());
    let mut o = options("");
    o.add_metainfo = false;
    let a = p.generate_spectrum(&o, &mut empty_db()).unwrap().value;
    assert_eq!(a.len(), 1);
    assert!(a.string_data_arrays.is_empty());
    assert_eq!(a.integer_data_arrays[0].data, vec![1]);
    o.ion_types = "I ??? BY".into();
    assert_eq!(p.generate_spectrum(&o, &mut empty_db()).unwrap().value, a);
    o.min_charge = -2;
    o.max_charge = -2;
    let negative = p.generate_spectrum(&o, &mut empty_db()).unwrap().value;
    assert!(negative.peaks[0].mz < 0.);
    o.min_charge = 2;
    o.max_charge = 1;
    assert!(
        p.generate_spectrum(&o, &mut empty_db())
            .unwrap()
            .value
            .is_empty()
    );
    o.min_charge = 0;
    o.max_charge = 0;
    assert!(p.generate_spectrum(&o, &mut empty_db()).is_err());
}
#[test]
fn late_conversion_and_backend_failures_roll_back_registry_and_ast() {
    let p = pf("M[Formula:Cl101]A");
    let before = p.clone();
    let mut db = empty_db();
    let mut o = options("by");
    o.min_charge = 0;
    assert!(p.generate_spectrum(&o, &mut db).is_err());
    assert!(db.is_empty());
    assert_eq!(p, before);
    let mut p = two(
        linked(
            ModificationTag::FormulaTag(FormulaTag {
                formula_string: "Cl101".into(),
                charge: None,
            }),
            "XL1",
            None,
        ),
        label(),
    );
    if let SequenceSection::Element(e) = &mut p.chains[1].sequence[1] {
        e.amino_acid = '?';
    }
    let before = p.clone();
    assert!(p.generate_spectrum(&options("by"), &mut db).is_err());
    assert!(db.is_empty());
    assert_eq!(p, before);
}

#[test]
fn warning_order_follows_alpha_conversion_before_beta_resolution() {
    let mut p = two(label(), linked(named("missing_beta"), "XL1", None));
    if let SequenceSection::Element(e) = &mut p.chains[0].sequence[1] {
        for (name, mass) in [("one", 100.), ("two", 200.)] {
            e.modifications.push(Modification {
                alternatives: vec![],
                resolved_mod: Some(Arc::new(
                    ResidueModification::from_record(ModificationRecord {
                        name: name.into(),
                        origin: Some('M'),
                        diff_mono_mass: mass,
                        diff_formula: EmpiricalFormula::parse("O").unwrap(),
                        ..Default::default()
                    })
                    .unwrap(),
                )),
            });
        }
    }
    let result = p
        .generate_spectrum(&options("by"), &mut empty_db())
        .unwrap();
    assert!(!result.value.is_empty());
    assert!(matches!(&result.warnings[..], [
        ConversionWarning::FormulaMassDisagreement { summed_mass: 300., .. },
        ConversionWarning::Resolution(ResolutionWarning::ModificationNotFound { name })
    ] if name == "missing_beta"));
}
