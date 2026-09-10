// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Anonymous definition behavior derived from ResidueModification::
//! createUnknownFromMassString, AASequence::parseModSquareBrackets_, and
//! ModificationDefinitionsSet at the pinned source revision.

use openms::Error;
use openms::chemistry::{
    AASequence, EmpiricalFormula, ModificationDefinition, ModificationDefinitionsSet,
    ModificationMassMode, ModificationMatchOptions, SequenceModification, TermSpecificity,
};
use openms::identification::{PeptideHit, PeptideIdentification};
use std::collections::BTreeSet;

fn sequence(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}
fn residue_annotation(peptide: &AASequence) -> &SequenceModification {
    let annotation = peptide.residue_modification(0).unwrap().unwrap();
    assert!(annotation.mass_tag().is_some());
    annotation
}
fn definition_set(annotation: &SequenceModification) -> ModificationDefinitionsSet {
    let mut definitions = ModificationDefinitionsSet::default();
    definitions
        .add_modification(ModificationDefinition::from_annotation(annotation, true, 0))
        .unwrap();
    definitions
}
fn mass(formula: &str) -> f64 {
    EmpiricalFormula::parse(formula).unwrap().mono_mass()
}
fn assert_matches(
    definitions: &ModificationDefinitionsSet,
    value: f64,
    mode: ModificationMassMode,
) {
    let found = definitions
        .find_matches(
            value,
            &ModificationMatchOptions {
                mass_mode: mode,
                tolerance: 0.0,
                ..Default::default()
            },
        )
        .unwrap();
    assert_eq!(found.len(), 1);
    assert_eq!(found[0].mass_error, 0.0);
    assert!(found[0].definition.mass_tag().is_some());
}
fn id(peptides: &[&str]) -> PeptideIdentification {
    PeptideIdentification {
        hits: peptides
            .iter()
            .map(|text| PeptideHit {
                sequence: sequence(text),
                ..Default::default()
            })
            .collect(),
        ..Default::default()
    }
}

#[test]
fn owned_annotation_survives_source_drop_and_preserves_settings_and_spelling() {
    let mut definition = {
        let peptide = sequence("K[+12.345678901]");
        ModificationDefinition::from_annotation(residue_annotation(&peptide), false, 7)
    };
    assert_eq!(definition.modification_name(), "K[+12.345678901]");
    let tag = definition.mass_tag().unwrap();
    assert_eq!(tag.input(), "+12.345678901");
    assert_eq!(tag.origin(), Some('K'));
    assert_eq!(tag.term_specificity(), TermSpecificity::Anywhere);
    assert_eq!(tag.full_id(), definition.modification_name());
    assert!(!definition.fixed);
    assert_eq!(definition.max_occurrences, 7);
    assert!(matches!(
        definition.modification(),
        Err(Error::Unsupported(_))
    ));
    // Constructing a definition never registers a global name as a side effect.
    assert!(ModificationDefinition::new(definition.modification_name()).is_err());
    let named = sequence("M(Oxidation)");
    definition.set_annotation(named.residue_modification(0).unwrap().unwrap());
    assert!(definition.mass_tag().is_none());
    assert_eq!(definition.modification().unwrap().name(), "Oxidation");
    assert!(!definition.fixed);
    assert_eq!(definition.max_occurrences, 7);
    let peptide = sequence("K[+12.345678901]");
    definition.set_annotation(residue_annotation(&peptide));
    let before = definition.clone();
    assert!(definition.set_modification("not-a-modification").is_err());
    assert_eq!(definition, before);
}

#[test]
fn source_empty_short_id_rule_precedes_exact_full_id_membership() {
    let fixed = sequence("K[+12.345678901]");
    let alternate = sequence("K[+13.345678901]");
    let mut definitions = definition_set(residue_annotation(&fixed));
    definitions
        .add_modification(ModificationDefinition::from_annotation(
            residue_annotation(&alternate),
            false,
            0,
        ))
        .unwrap();
    definitions
        .add_modification(ModificationDefinition::with_options("Acetyl (K)", false, 0).unwrap())
        .unwrap();
    // Source user-defined modifications all have empty short IDs. Either listed
    // anonymous tag passes the fixed short-name check, then full-ID membership.
    assert!(definitions.is_compatible(&fixed).unwrap());
    assert!(definitions.is_compatible(&alternate).unwrap());
    assert!(!definitions.is_compatible(&sequence("K")).unwrap());
    assert!(
        !definitions
            .is_compatible(&sequence("K[+14.345678901]"))
            .unwrap()
    );
    assert!(!definitions.is_compatible(&sequence("K(Acetyl)")).unwrap());
}

#[test]
fn source_residue_mass_search_uses_full_anchor_and_original_operation_order() {
    let full = mass("C6H14N2O2");
    let internal = full - mass("H2O");
    let delta_peptide = sequence("K[+12.345678901]");
    let definitions = definition_set(residue_annotation(&delta_peptide));
    assert_matches(&definitions, 12.345678901, ModificationMassMode::Delta);
    assert_matches(
        &definitions,
        12.345678901 + full,
        ModificationMassMode::Absolute,
    );

    let absolute_peptide = sequence("K[211.123456789]");
    let definitions = definition_set(residue_annotation(&absolute_peptide));
    let source_delta = 211.123456789 - internal;
    assert_matches(&definitions, source_delta, ModificationMassMode::Delta);
    assert_matches(
        &definitions,
        source_delta + full,
        ModificationMassMode::Absolute,
    );
    assert!(
        definitions
            .find_matches(
                211.123456789,
                &ModificationMatchOptions {
                    mass_mode: ModificationMassMode::Absolute,
                    tolerance: 0.0,
                    ..Default::default()
                }
            )
            .unwrap()
            .is_empty(),
        "internal tag value is not the source full-residue anchor"
    );
}

#[test]
fn terminal_tags_use_h_and_oh_anchors_for_delta_and_absolute_spelling() {
    for (text, input, delta, absolute, term) in [
        (
            ".[+11.123456789]A",
            "+11.123456789",
            11.123456789,
            11.123456789 + mass("H"),
            TermSpecificity::NTerm,
        ),
        (
            "A.[+12.987654321]",
            "+12.987654321",
            12.987654321,
            12.987654321 + mass("HO"),
            TermSpecificity::CTerm,
        ),
        (
            ".[31.222222222]A",
            "31.222222222",
            31.222222222 - mass("H"),
            31.222222222,
            TermSpecificity::NTerm,
        ),
        (
            "A.[21.111111111]",
            "21.111111111",
            21.111111111 - mass("HO"),
            21.111111111,
            TermSpecificity::CTerm,
        ),
    ] {
        let peptide = sequence(text);
        let annotation = if term == TermSpecificity::NTerm {
            peptide.n_terminal_modification()
        } else {
            peptide.c_terminal_modification()
        }
        .unwrap();
        let definitions = definition_set(annotation);
        assert_eq!(
            definitions
                .fixed_modifications()
                .next()
                .unwrap()
                .mass_tag()
                .unwrap()
                .input(),
            input
        );
        assert_matches(&definitions, delta, ModificationMassMode::Delta);
        assert_matches(&definitions, absolute, ModificationMassMode::Absolute);
        assert!(definitions.is_compatible(&peptide).unwrap());
        // Source generic terminal origin X does not force ordinary terminal slots.
        assert!(definitions.is_compatible(&sequence("A")).unwrap());
    }
}

#[test]
fn unresolved_absolute_tags_are_retained_without_inventing_unmodified_deltas() {
    for residue in ['B', 'Z', 'X'] {
        let peptide = sequence(&format!("{residue}[201.123456789]"));
        let definitions = definition_set(residue_annotation(&peptide));
        assert_matches(
            &definitions,
            201.123456789 + mass("H2O"),
            ModificationMassMode::Absolute,
        );
        let tag = definitions
            .fixed_modifications()
            .next()
            .unwrap()
            .mass_tag()
            .unwrap();
        assert_eq!(tag.residue_mono_mass(), Some(201.123456789));
        assert_eq!(tag.delta_mono_mass(), None);
        assert!(matches!(
            definitions.find_matches(1.0, &ModificationMatchOptions::default()),
            Err(Error::Unsupported(_))
        ));
        assert!(
            definitions
                .find_matches(
                    1.0,
                    &ModificationMatchOptions {
                        term_specificity: Some(TermSpecificity::NTerm),
                        ..Default::default()
                    }
                )
                .unwrap()
                .is_empty(),
            "unmatched unknown chemistry need not be evaluated"
        );
        assert!(definitions.is_compatible(&peptide).unwrap());
    }
}

#[test]
fn anonymous_inference_preserves_lexeme_identity_and_checks_payload_before_mutation() {
    let mut definitions = ModificationDefinitionsSet::default();
    definitions
        .infer_from_peptides(&[id(&["K[+12.345678901]", "K[+12.3456789010]", "K"])])
        .unwrap();
    assert!(definitions.fixed_names().is_empty());
    assert_eq!(
        definitions.variable_names(),
        BTreeSet::from([
            "K[+12.345678901]".to_string(),
            "K[+12.3456789010]".to_string()
        ])
    );
    assert!(
        definitions
            .modifications()
            .iter()
            .all(|definition| definition.mass_tag().is_some())
    );
    let ids = [id(&["K[+12.345678901]"])];
    let full_id_bytes = ids[0].hits[0]
        .sequence
        .residue_modification(0)
        .unwrap()
        .unwrap()
        .full_id()
        .len();
    let needed = 1 + 1 + 1 + 2 + full_id_bytes; // identification, hit, residue, termini, spelling
    definitions.max_work = needed - 1;
    let before = definitions.clone();
    assert!(definitions.infer_from_peptides(&ids).is_err());
    assert_eq!(definitions, before);
    definitions.max_work = needed;
    definitions.infer_from_peptides(&ids).unwrap();
    assert_eq!(
        definitions.fixed_names(),
        BTreeSet::from(["K[+12.345678901]".to_string()])
    );
    assert!(definitions.variable_names().is_empty());
}
