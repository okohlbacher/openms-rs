// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Literal source compatibility/inference/matching assertions plus independent
//! boundary checks. See data/modification_generation_provenance.json.

use openms::chemistry::{
    AASequence, ModificationDefinition, ModificationDefinitionsSet, ModificationMassMode,
    ModificationMatchOptions, ModificationsDB, TermSpecificity,
};
use openms::identification::{PeptideHit, PeptideIdentification};
use std::collections::BTreeSet;

fn peptide(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}
fn identifiers(texts: &[&str]) -> PeptideIdentification {
    PeptideIdentification {
        hits: texts
            .iter()
            .map(|s| PeptideHit {
                sequence: peptide(s),
                ..Default::default()
            })
            .collect(),
        ..Default::default()
    }
}
fn names(values: &[&str]) -> BTreeSet<String> {
    values.iter().map(|s| (*s).to_owned()).collect()
}

#[test]
fn literal_source_compatibility_cases_and_unenforced_counts() {
    for row in include_str!("data/modification_generation_compatibility.tsv")
        .lines()
        .skip(1)
    {
        let fields: Vec<_> = row.split('\t').collect();
        let fixed: Vec<_> = fields[1].split(';').collect();
        let variable: Vec<_> = fields[2].split(';').collect();
        let defs = ModificationDefinitionsSet::from_names(&fixed, &variable).unwrap();
        assert_eq!(
            defs.is_compatible(&peptide(fields[3])).unwrap(),
            fields[4].parse::<bool>().unwrap()
        );
    }
    let mut defs = ModificationDefinitionsSet::default();
    defs.add_modification(
        ModificationDefinition::with_options("Carbamidomethyl (C)", true, 1).unwrap(),
    )
    .unwrap();
    defs.add_modification(ModificationDefinition::with_options("Phospho (S)", false, 1).unwrap())
        .unwrap();
    defs.max_modifications = 1;
    assert!(
        defs.is_compatible(&peptide(
            "C(Carbamidomethyl)C(Carbamidomethyl)S(Phospho)S(Phospho)"
        ))
        .unwrap()
    );
    let terminal = ModificationDefinitionsSet::from_names(&["Acetyl (N-term)"], &[]).unwrap();
    assert!(terminal.is_compatible(&peptide("AAA")).unwrap());
}

#[test]
fn source_definition_set_identity_and_merged_fixed_precedence() {
    let mut defs = ModificationDefinitionsSet::default();
    defs.add_modification(ModificationDefinition::with_options("Oxidation (M)", true, 1).unwrap())
        .unwrap();
    defs.add_modification(ModificationDefinition::with_options("Oxidation (M)", true, 5).unwrap())
        .unwrap();
    defs.add_modification(ModificationDefinition::with_options("Oxidation (M)", false, 2).unwrap())
        .unwrap();
    assert_eq!(defs.len(), 2);
    assert_eq!(
        defs.fixed_modifications().next().unwrap().max_occurrences,
        1
    );
    assert_eq!(
        defs.variable_modifications()
            .next()
            .unwrap()
            .max_occurrences,
        2
    );
    assert_eq!(defs.modifications().len(), 1);
    assert!(defs.modifications()[0].fixed);
    assert_eq!(defs.modifications()[0].max_occurrences, 1);
    let before = defs.clone();
    assert!(
        defs.set_names(&["Carbamidomethyl (C)"], &["not-a-modification"])
            .is_err()
    );
    assert_eq!(defs, before);
    assert!(
        defs.add_modification(ModificationDefinition::default())
            .is_err()
    );
    assert_eq!(defs, before);
}

#[test]
fn source_negative_delta_matching_and_inclusive_error_boundary() {
    let defs = ModificationDefinitionsSet::from_names(
        &["Gln->pyro-Glu (N-term Q)"],
        &["Glu->pyro-Glu (N-term E)", "Oxidation (M)"],
    )
    .unwrap();
    let mut options = ModificationMatchOptions {
        residue: "E".into(),
        term_specificity: Some(TermSpecificity::NTerm),
        tolerance: 0.1,
        ..Default::default()
    };
    assert_eq!(
        defs.find_matches(-18., &options).unwrap()[0]
            .definition
            .modification_name(),
        "Glu->pyro-Glu (N-term E)"
    );
    options.term_specificity = Some(TermSpecificity::Anywhere);
    assert!(defs.find_matches(-18., &options).unwrap().is_empty());
    options.term_specificity = Some(TermSpecificity::NTerm);
    options.residue = "Q".into();
    assert!(defs.find_matches(-18., &options).unwrap().is_empty());
    options.residue.clear();
    options.tolerance = 2.;
    assert_eq!(
        defs.find_matches(-18., &options)
            .unwrap()
            .iter()
            .map(|m| m.definition.modification_name())
            .collect::<Vec<_>>(),
        ["Glu->pyro-Glu (N-term E)", "Gln->pyro-Glu (N-term Q)"]
    );
    options.consider_fixed = false;
    options.consider_variable = false;
    assert!(defs.find_matches(-18., &options).is_err());

    let ox =
        ModificationDefinitionsSet::from_names(&["Oxidation (M)"], &["Oxidation (M)"]).unwrap();
    let delta = ox
        .fixed_modifications()
        .next()
        .unwrap()
        .modification()
        .unwrap()
        .diff_mono_mass();
    let options = ModificationMatchOptions {
        tolerance: delta,
        ..Default::default()
    };
    let matches = ox.find_matches(0., &options).unwrap();
    assert_eq!(matches.len(), 2);
    assert!(matches[0].definition.fixed && !matches[1].definition.fixed);
    let options = ModificationMatchOptions {
        tolerance: f64::from_bits(delta.to_bits() - 1),
        ..options
    };
    assert!(ox.find_matches(0., &options).unwrap().is_empty());
}

#[test]
fn absolute_mass_matching_preserves_stored_mass_and_source_residue_fallback() {
    let db = ModificationsDB::global();
    let oxidation = db.get_modification("Oxidation (M)", None, None).unwrap();
    let mut defs = ModificationDefinitionsSet::from_names(&["Oxidation (M)"], &[]).unwrap();
    // Source fallback computes full-residue mass, subtracts water, then adds delta.
    let mass = peptide("M").formula().unwrap().mono_mass()
        - "H2O"
            .parse::<openms::chemistry::EmpiricalFormula>()
            .unwrap()
            .mono_mass()
        + oxidation.diff_mono_mass();
    let mut options = ModificationMatchOptions {
        residue: "M".into(),
        mass_mode: ModificationMassMode::Absolute,
        tolerance: 0.,
        ..Default::default()
    };
    assert_eq!(defs.find_matches(mass, &options).unwrap().len(), 1);
    options.residue = "Methionine".into();
    assert_eq!(defs.find_matches(mass, &options).unwrap().len(), 1);
    options.residue.clear();
    assert_eq!(defs.find_matches(0., &options).unwrap().len(), 1);
    options.residue = "X".into();
    assert!(defs.find_matches(mass, &options).is_err());

    let absolute = db
        .get_modification("Acetyl (N-term)", None, None)
        .unwrap()
        .clone()
        .with_absolute_masses(123.456, 124.5)
        .unwrap();
    assert_eq!(absolute.mono_mass(), 123.456);
    assert_eq!(absolute.average_mass(), 124.5);
    defs.set_modifications(&[ModificationDefinition::from_modification(
        &absolute, true, 0,
    )])
    .unwrap();
    // Positive stored mass bypasses residue lookup, even for ambiguous or unknown
    // selectors accepted by the source's origin-X filter.
    for residue in ["", "B", "Z", "X", "not-a-residue"] {
        options.residue = residue.into();
        let matched = defs.find_matches(123.456, &options).unwrap();
        assert_eq!(matched.len(), 1);
        assert_eq!(matched[0].mass_error, 0.);
    }
    assert!(
        absolute
            .clone()
            .with_absolute_masses(f64::INFINITY, 0.)
            .is_err()
    );
    assert!(absolute.with_absolute_masses(0., f64::NAN).is_err());
}

#[test]
fn source_inference_pools_all_hits_and_retains_terminal_predicate_quirk() {
    let ids = [
        identifiers(&["AC(Carbamidomethyl)M", "(Acetyl)AEM"]),
        identifiers(&["AC(Carbamidomethyl)M(Oxidation)"]),
    ];
    let mut defs = ModificationDefinitionsSet::default();
    defs.max_modifications = 7;
    defs.infer_from_peptides(&ids).unwrap();
    assert_eq!(defs.fixed_names(), names(&["Carbamidomethyl (C)"]));
    assert_eq!(
        defs.variable_names(),
        names(&["Acetyl (N-term)", "Oxidation (M)"])
    );
    assert_eq!(defs.max_modifications, 7);
    let pyro = peptide(".(Gln->pyro-Glu)QAA");
    defs.infer_from_peptides(&[identifiers(&[".(Gln->pyro-Glu)QAA"])])
        .unwrap();
    assert_eq!(defs.fixed_names(), names(&["Gln->pyro-Glu (N-term Q)"]));
    assert!(!defs.is_compatible(&pyro).unwrap());
    defs.infer_from_peptides(&[identifiers(&["C[999]"])])
        .unwrap();
    let anonymous = peptide("C[999]");
    let annotation = anonymous.residue_modification(0).unwrap().unwrap();
    assert_eq!(defs.fixed_names(), names(&[annotation.full_id()]));
    assert_eq!(
        defs.fixed_modifications().next().unwrap().mass_tag(),
        annotation.mass_tag()
    );
}
