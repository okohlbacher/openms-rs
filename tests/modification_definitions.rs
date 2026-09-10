// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Source: pinned ModificationDefinition/ModificationDefinitionsSet class tests.

use openms::Error;
use openms::chemistry::{
    AASequence, EmpiricalFormula, ModificationDefinition, ModificationDefinitionsSet,
    ModificationMassMode, ModificationMatchOptions, ModificationsDB, ModifiedPeptideGenerator,
    ResidueModification, TermSpecificity,
};
use openms::identification::{PeptideHit, PeptideIdentification};
use std::collections::{BTreeSet, HashSet};
use std::sync::Arc;

fn sequence(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}
fn names(values: &[&str]) -> BTreeSet<String> {
    values.iter().map(|name| (*name).to_owned()).collect()
}
fn id(sequences: &[&str]) -> PeptideIdentification {
    PeptideIdentification {
        hits: sequences
            .iter()
            .enumerate()
            .map(|(index, text)| {
                PeptideHit::new(index as f64, index as u32, 2, sequence(text)).unwrap()
            })
            .collect(),
        ..Default::default()
    }
}
fn options(residue: &str, mode: ModificationMassMode) -> ModificationMatchOptions {
    ModificationMatchOptions {
        residue: residue.into(),
        mass_mode: mode,
        ..Default::default()
    }
}

#[test]
fn source_defaults_mutation_and_safe_hash_equality() {
    let mut definition = ModificationDefinition::default();
    assert!(definition.fixed);
    assert_eq!(definition.max_occurrences, 0);
    assert_eq!(definition.modification_name(), "");
    assert!(definition.modification().is_err());
    definition.set_modification("Acetyl (N-term)").unwrap();
    let saved = definition.clone();
    assert!(definition.set_modification("missing modification").is_err());
    assert_eq!(definition, saved);
    assert_eq!(
        definition,
        ModificationDefinition::new("Acetyl (N-term)").unwrap()
    );
    definition.fixed = false;
    assert_ne!(definition, saved);
    definition.fixed = true;
    definition.max_occurrences = 2;
    assert_ne!(definition, saved);
    let carboxymethyl =
        ModificationDefinition::with_options("Carboxymethyl (C)", false, 2).unwrap();
    assert!(!carboxymethyl.fixed);
    assert_eq!(carboxymethyl.max_occurrences, 2);
    let hashes = HashSet::from([
        ModificationDefinition::default(),
        saved.clone(),
        saved,
        definition,
    ]);
    assert_eq!(hashes.len(), 3);
}

#[test]
fn source_partition_identity_first_wins_and_fixed_wins_merged_view() {
    let mut definitions = ModificationDefinitionsSet::from_names(
        &["Phospho (Y)", "Phospho (T)", "Phospho (S)"],
        &["Carbamidomethyl (C)", "Phospho (S)"],
    )
    .unwrap();
    assert_eq!(definitions.len(), 5);
    assert_eq!(definitions.fixed_modifications().len(), 3);
    assert_eq!(definitions.variable_modifications().len(), 2);
    assert_eq!(definitions.modification_names().len(), 4);
    assert_eq!(
        definitions
            .modifications()
            .iter()
            .map(|m| m.modification_name())
            .collect::<Vec<_>>(),
        [
            "Carbamidomethyl (C)",
            "Phospho (S)",
            "Phospho (T)",
            "Phospho (Y)"
        ]
    );
    assert!(definitions.modifications()[1].fixed);
    definitions
        .add_modification(ModificationDefinition::with_options("Phospho (S)", true, 10).unwrap())
        .unwrap();
    assert_eq!(
        definitions
            .fixed_modifications()
            .next()
            .unwrap()
            .max_occurrences,
        0
    );
    let variable = ModificationDefinition::with_options("Phospho (S)", false, 3).unwrap();
    let fixed = ModificationDefinition::new("Phospho (S)").unwrap();
    definitions
        .set_modifications(&[variable.clone(), fixed])
        .unwrap();
    assert_eq!(definitions.len(), 1);
    assert_eq!(definitions.variable_modifications().next(), Some(&variable));
}

#[test]
fn replacement_is_atomic_and_csv_preserves_untrimmed_source_tokens() {
    let mut definitions = ModificationDefinitionsSet::from_comma_separated(
        "Phospho (S),Phospho (T),Phospho (Y)",
        "Carbamidomethyl (C)",
    )
    .unwrap();
    definitions.max_modifications = 7;
    let saved = definitions.clone();
    for (fixed, variable) in [
        ("Phospho (S)", "unknown"),
        ("Phospho (S), Phospho (T)", ""),
        ("Phospho (S),", ""),
    ] {
        assert!(definitions.set_comma_separated(fixed, variable).is_err());
        assert_eq!(definitions, saved);
    }
    assert!(
        definitions
            .add_modification(ModificationDefinition::default())
            .is_err()
    );
    assert!(
        definitions
            .set_modifications(&[ModificationDefinition::default()])
            .is_err()
    );
    assert_eq!(definitions, saved);
    definitions.set_names(&[], &[]).unwrap();
    assert!(definitions.is_empty());
    assert_eq!(definitions.max_modifications, 7);
}

#[test]
fn compatibility_uses_source_fixed_residue_rule_and_ignores_stored_counts() {
    let mut definitions = ModificationDefinitionsSet::from_names(
        &["Carbamidomethyl (C)"],
        &["Phospho (S)", "Phospho (T)", "Phospho (Y)"],
    )
    .unwrap();
    for (peptide, expected) in [
        ("CCTKPESER", false),
        ("C(Carbamidomethyl)CTKPESER", false),
        ("C(Carbamidomethyl)C(Carbamidomethyl)TKPESER", true),
        (
            "C(Carbamidomethyl)C(Carbamidomethyl)T(Phospho)TKPESER",
            true,
        ),
        ("(Acetyl)CCTKPESER", false),
        (
            "(Acetyl)C(Carbamidomethyl)C(Carbamidomethyl)TKPES(Phospho)ER",
            false,
        ),
        (
            "(Acetyl)C(Carbamidomethyl)C(Carbamidomethyl)T(Phospho)KPES(Phospho)ER",
            false,
        ),
    ] {
        assert_eq!(
            definitions.is_compatible(&sequence(peptide)).unwrap(),
            expected,
            "{peptide}"
        );
    }
    definitions.max_modifications = 1;
    definitions
        .set_modifications(&[
            ModificationDefinition::with_options("Carbamidomethyl (C)", true, 1).unwrap(),
            ModificationDefinition::with_options("Phospho (S)", false, 1).unwrap(),
        ])
        .unwrap();
    assert!(
        definitions
            .is_compatible(&sequence(
                "C(Carbamidomethyl)C(Carbamidomethyl)S(Phospho)S(Phospho)"
            ))
            .unwrap()
    );
}

#[test]
fn terminal_fixed_predicate_retains_source_quirks() {
    let acetyl = ModificationDefinitionsSet::from_names(&["Acetyl (N-term)"], &[]).unwrap();
    assert!(acetyl.is_compatible(&sequence("PEPTIDE")).unwrap());
    assert!(acetyl.is_compatible(&sequence("(Acetyl)PEPTIDE")).unwrap());
    assert!(!acetyl.is_compatible(&sequence("XPEPTIDE")).unwrap());
    let pyro = ModificationDefinitionsSet::from_names(&["Gln->pyro-Glu (N-term Q)"], &[]).unwrap();
    assert!(
        !pyro
            .is_compatible(&sequence("(Gln->pyro-Glu)QPEPTIDE"))
            .unwrap()
    );
    assert!(
        ModificationDefinitionsSet::default()
            .is_compatible(&sequence("BZX"))
            .unwrap()
    );
}

#[test]
fn source_delta_match_goldens_and_fixed_variable_filters() {
    let definitions = ModificationDefinitionsSet::from_names(
        &["Gln->pyro-Glu (N-term Q)"],
        &["Glu->pyro-Glu (N-term E)", "Oxidation (M)"],
    )
    .unwrap();
    let mut query = ModificationMatchOptions {
        residue: "E".into(),
        term_specificity: Some(TermSpecificity::NTerm),
        tolerance: 0.1,
        ..Default::default()
    };
    let matches = definitions.find_matches(-18.0, &query).unwrap();
    assert_eq!(matches.len(), 1);
    assert_eq!(
        matches[0].definition.modification_name(),
        "Glu->pyro-Glu (N-term E)"
    );
    query.term_specificity = Some(TermSpecificity::Anywhere);
    assert!(definitions.find_matches(-18.0, &query).unwrap().is_empty());
    query.term_specificity = Some(TermSpecificity::NTerm);
    query.residue = "Q".into();
    assert!(definitions.find_matches(-18.0, &query).unwrap().is_empty());
    query.residue.clear();
    query.consider_variable = false;
    assert!(definitions.find_matches(-18.0, &query).unwrap().is_empty());
    query.consider_variable = true;
    query.tolerance = 2.0;
    assert_eq!(
        definitions
            .find_matches(-18.0, &query)
            .unwrap()
            .iter()
            .map(|m| m.definition.modification_name())
            .collect::<Vec<_>>(),
        ["Glu->pyro-Glu (N-term E)", "Gln->pyro-Glu (N-term Q)"]
    );
}

#[test]
fn equal_mass_errors_keep_partition_and_full_id_order_inclusively() {
    let definitions = ModificationDefinitionsSet::from_names(
        &["Phospho (Y)", "Phospho (S)"],
        &["Phospho (T)", "Phospho (S)"],
    )
    .unwrap();
    let delta = definitions
        .fixed_modifications()
        .next()
        .unwrap()
        .modification()
        .unwrap()
        .diff_mono_mass();
    let query = ModificationMatchOptions {
        tolerance: 0.0,
        ..Default::default()
    };
    let matches = definitions.find_matches(delta, &query).unwrap();
    assert_eq!(
        matches
            .iter()
            .map(|m| (m.definition.modification_name(), m.definition.fixed))
            .collect::<Vec<_>>(),
        [
            ("Phospho (S)", true),
            ("Phospho (Y)", true),
            ("Phospho (S)", false),
            ("Phospho (T)", false)
        ]
    );
    assert!(matches.iter().all(|m| m.mass_error == 0.0));
    let distance = (delta - 80.0).abs();
    let at_boundary = ModificationMatchOptions {
        tolerance: distance,
        ..Default::default()
    };
    assert_eq!(
        definitions.find_matches(80.0, &at_boundary).unwrap().len(),
        4
    );
}

#[test]
fn absolute_fallback_uses_full_formula_minus_water_and_source_alias_filter() {
    let definitions = ModificationDefinitionsSet::from_names(&[], &["Oxidation (M)"]).unwrap();
    let delta = definitions
        .variable_modifications()
        .next()
        .unwrap()
        .modification()
        .unwrap()
        .diff_mono_mass();
    let absolute = delta
        + (EmpiricalFormula::parse("C5H11NO2S").unwrap().mono_mass()
            - EmpiricalFormula::parse("H2O").unwrap().mono_mass());
    for residue in ["M", "Methionine", "Met", "MET"] {
        let query = ModificationMatchOptions {
            tolerance: 0.0,
            ..options(residue, ModificationMassMode::Absolute)
        };
        let matches = definitions.find_matches(absolute, &query).unwrap();
        assert_eq!(matches.len(), 1, "{residue}");
        assert_eq!(matches[0].mass_error.to_bits(), 0.0_f64.to_bits());
    }
    let empty = options("", ModificationMassMode::Absolute);
    assert!(
        definitions
            .find_matches(absolute, &empty)
            .unwrap()
            .is_empty()
    );
    assert_eq!(definitions.find_matches(0.0, &empty).unwrap().len(), 1);
    assert!(
        definitions
            .find_matches(
                absolute,
                &options("methionine", ModificationMassMode::Absolute)
            )
            .unwrap()
            .is_empty()
    );
    for residue in [".", "X", "Missing"] {
        assert!(
            definitions
                .find_matches(absolute, &options(residue, ModificationMassMode::Absolute))
                .is_err()
        );
    }
}

#[test]
fn explicit_absolute_masses_override_fallback_without_changing_registry_or_delta() {
    let db = ModificationsDB::global();
    let original = db.get_modification("Acetyl (N-term)", None, None).unwrap();
    assert_eq!(original.mono_mass(), 0.0);
    assert_eq!(original.average_mass(), 0.0);
    let owned = original.clone().with_absolute_masses(123.5, 124.0).unwrap();
    assert_eq!(owned.mono_mass(), 123.5);
    assert_eq!(owned.average_mass(), 124.0);
    assert_eq!(owned.diff_mono_mass(), original.diff_mono_mass());
    assert_ne!(&owned, original);
    let mut definitions = ModificationDefinitionsSet::default();
    definitions
        .add_modification(ModificationDefinition::from_modification(&owned, true, 0))
        .unwrap();
    for residue in ["", "M", "B", "Z", "X", ".", "Missing"] {
        assert_eq!(
            definitions
                .find_matches(123.5, &options(residue, ModificationMassMode::Absolute))
                .unwrap()
                .len(),
            1
        );
    }
    assert_eq!(
        definitions
            .find_matches(
                original.diff_mono_mass(),
                &ModificationMatchOptions::default()
            )
            .unwrap()
            .len(),
        1
    );
    for (mono, average) in [(f64::NAN, 1.0), (1.0, f64::INFINITY)] {
        assert!(
            original
                .clone()
                .with_absolute_masses(mono, average)
                .is_err()
        );
    }
    let negative = original.clone().with_absolute_masses(-2.0, -3.0).unwrap();
    definitions
        .set_modifications(&[ModificationDefinition::from_modification(
            &negative, true, 0,
        )])
        .unwrap();
    assert_eq!(
        definitions
            .find_matches(-2.0, &options("", ModificationMassMode::Absolute))
            .unwrap()
            .len(),
        1
    );
    assert!(matches!(
        definitions.find_matches(123.5, &options("X", ModificationMassMode::Absolute)),
        Err(Error::Unsupported(_))
    ));
    let huge = original
        .clone()
        .with_absolute_masses(f64::MAX, f64::MAX)
        .unwrap();
    definitions
        .set_modifications(&[ModificationDefinition::from_modification(&huge, true, 0)])
        .unwrap();
    assert!(
        definitions
            .find_matches(-f64::MAX, &options("", ModificationMassMode::Absolute))
            .is_err()
    );
}

#[test]
fn source_inference_pools_all_hits_and_keeps_absent_sites_out_of_null_counts() {
    let mut definitions = ModificationDefinitionsSet::default();
    definitions.max_modifications = 1;
    let peptides = [
        id(&["AC(Carbamidomethyl)M", "(Acetyl)AEM"]),
        id(&["AC(Carbamidomethyl)M(Oxidation)"]),
    ];
    definitions.infer_from_peptides(&peptides).unwrap();
    assert_eq!(definitions.fixed_names(), names(&["Carbamidomethyl (C)"]));
    assert_eq!(
        definitions.variable_names(),
        names(&["Acetyl (N-term)", "Oxidation (M)"])
    );
    assert_eq!(definitions.max_modifications, 1);
    assert!(
        definitions
            .modifications()
            .iter()
            .all(|m| m.max_occurrences == 0)
    );
    for hit in peptides.iter().flat_map(|id| &id.hits) {
        assert!(definitions.is_compatible(&hit.sequence).unwrap());
    }
    definitions
        .infer_from_peptides(&[id(&["M(Oxidation)M"])])
        .unwrap();
    assert!(definitions.fixed_names().is_empty());
    assert_eq!(definitions.variable_names(), names(&["Oxidation (M)"]));
    definitions
        .infer_from_peptides(&[id(&["(Acetyl)A", ""])])
        .unwrap();
    assert_eq!(definitions.variable_names(), names(&["Acetyl (N-term)"]));
    definitions.infer_from_peptides(&[]).unwrap();
    assert!(definitions.is_empty());
    assert_eq!(definitions.max_modifications, 1);
}

#[test]
fn anonymous_known_mass_is_not_an_invented_registry_definition() {
    let mut definitions =
        ModificationDefinitionsSet::from_names(&["Carbamidomethyl (C)"], &[]).unwrap();
    for text in ["K[+12.3456789]", "X[147.035399]"] {
        let peptides = [id(&["C(Carbamidomethyl)", text])];
        assert!(
            peptides[0].hits[1]
                .sequence
                .residue_modification(0)
                .unwrap()
                .unwrap()
                .mass_tag()
                .is_some()
        );
        definitions.infer_from_peptides(&peptides).unwrap();
        let anonymous = definitions
            .fixed_modifications()
            .find(|definition| definition.mass_tag().is_some())
            .unwrap();
        assert_eq!(anonymous.modification_name(), text);
        assert!(matches!(
            anonymous.modification(),
            Err(Error::Unsupported(_))
        ));
        assert_eq!(definitions.len(), 2);
        assert!(
            definitions
                .is_compatible(&peptides[0].hits[1].sequence)
                .unwrap()
        );
    }
}

#[test]
fn numerical_errors_and_resource_limits_are_checked_atomically() {
    let mut definitions =
        ModificationDefinitionsSet::from_names(&["Acetyl (N-term)"], &["Oxidation (M)"]).unwrap();
    for (mass, tolerance) in [(f64::NAN, 1.0), (1.0, f64::INFINITY), (1.0, -1.0)] {
        let query = ModificationMatchOptions {
            tolerance,
            ..Default::default()
        };
        assert!(definitions.find_matches(mass, &query).is_err());
    }
    assert!(
        definitions
            .find_matches(
                0.0,
                &ModificationMatchOptions {
                    consider_fixed: false,
                    consider_variable: false,
                    ..Default::default()
                }
            )
            .is_err()
    );
    assert!(
        definitions
            .find_matches(0.0, &options("é", ModificationMassMode::Delta))
            .is_err()
    );
    let limitless = definitions.clone();
    definitions.max_work = 1;
    assert_eq!(
        definitions, limitless,
        "resource settings are not chemical equality"
    );
    assert!(
        definitions
            .infer_from_peptides(&[id(&["C(Carbamidomethyl)"])])
            .is_err()
    );
    assert!(definitions.is_compatible(&sequence("PEPTIDE")).is_err());
    assert!(
        definitions
            .find_matches(0.0, &ModificationMatchOptions::default())
            .is_err()
    );
    assert!(definitions.set_names(&[], &["Oxidation (M)"]).is_err());
    assert_eq!(definitions, limitless);
    definitions.max_work = 0;
    assert!(definitions.set_names(&[], &[]).is_err());
    assert_eq!(definitions, limitless);
}

#[test]
fn inferred_full_id_cannot_silently_choose_conflicting_custom_chemistry() {
    let original = ModificationsDB::global()
        .get_modification("Oxidation (M)", None, None)
        .unwrap();
    // Owned registry handles outlive the caller-owned record/registry safely.
    let changed = Arc::new(original.clone().with_absolute_masses(123.0, 124.0).unwrap());
    let identical_copy = Arc::new(original.clone());
    let modified = |modification: Arc<ResidueModification>| {
        let mut peptide = sequence("M");
        ModifiedPeptideGenerator::default()
            .apply_fixed_modifications(&[modification], &mut peptide)
            .unwrap();
        PeptideHit {
            sequence: peptide,
            ..Default::default()
        }
    };
    let original_hit = modified(Arc::new(original.clone()));
    let changed_hit = modified(changed);
    let mut definitions =
        ModificationDefinitionsSet::from_names(&[], &["Carbamidomethyl (C)"]).unwrap();
    let before = definitions.clone();
    for hits in [
        vec![original_hit.clone(), changed_hit.clone()],
        vec![changed_hit, original_hit.clone()],
    ] {
        let ids = [PeptideIdentification {
            hits,
            ..Default::default()
        }];
        assert!(matches!(
            definitions.infer_from_peptides(&ids),
            Err(Error::InvalidValue(_))
        ));
        assert_eq!(definitions, before);
    }
    definitions
        .infer_from_peptides(&[PeptideIdentification {
            hits: vec![original_hit, modified(identical_copy)],
            ..Default::default()
        }])
        .unwrap();
    assert_eq!(definitions.fixed_names(), names(&["Oxidation (M)"]));
    assert!(definitions.variable_names().is_empty());
    assert_eq!(
        definitions
            .fixed_modifications()
            .next()
            .unwrap()
            .modification()
            .unwrap(),
        original
    );
}
