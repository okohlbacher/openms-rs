// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::{
    AASequence, EmpiricalFormula, ModificationsDB, ModifiedPeptideGenerator as Generator,
    ResidueModification,
};
use std::sync::Arc;

fn seq(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}
fn mods(names: &[&str]) -> Vec<Arc<ResidueModification>> {
    Generator::get_modifications(names).unwrap()
}
fn strings(peptides: &[AASequence]) -> Vec<String> {
    peptides.iter().map(ToString::to_string).collect()
}

#[test]
fn pinned_fixed_carbamidomethyl_and_terminal_examples() {
    // Literal cases from ModifiedPeptideGenerator_test.cpp at revision 7c029e8.
    let fixed = mods(&["Carbamidomethyl (C)"]);
    for (input, expected) in [
        ("AAAACAAAA", "AAAAC(Carbamidomethyl)AAAA"),
        ("AAAAAAAAA", "AAAAAAAAA"),
        (
            "AAAACAAC(Carbamidomethyl)AAA",
            "AAAAC(Carbamidomethyl)AAC(Carbamidomethyl)AAA",
        ),
        (
            "AAAACAAC(Oxidation)AAA",
            "AAAAC(Carbamidomethyl)AAC(Oxidation)AAA",
        ),
    ] {
        let mut peptide = seq(input);
        Generator::default()
            .apply_fixed_modifications(&fixed, &mut peptide)
            .unwrap();
        assert_eq!(peptide, seq(expected));
    }
    let mut peptide = seq("K(Carbamyl)AAAAAAAA");
    Generator::default()
        .apply_fixed_modifications(&mods(&["Carbamyl (N-term)"]), &mut peptide)
        .unwrap();
    assert_eq!(peptide, seq(".(Carbamyl)K(Carbamyl)AAAAAAAA"));
}

#[test]
fn reverse_sites_and_weighted_26_71_and_four_variant_source_goldens() {
    let generator = Generator::default();
    let ox = mods(&["Oxidation (M)"]);
    let input = seq("AAMAAAMAA");
    let expected = [
        "AAMAAAM(Oxidation)AA",
        "AAM(Oxidation)AAAMAA",
        "AAM(Oxidation)AAAM(Oxidation)AA",
    ];
    assert_eq!(
        strings(
            &generator
                .variable_modifications(&ox, &input, 2, false)
                .unwrap()
        ),
        expected
    );
    assert_eq!(
        strings(
            &generator
                .variable_modifications(&ox, &input, 1, false)
                .unwrap()
        ),
        expected[..2]
    );
    let alternatives = mods(&["Glutathione (C)", "Carbamidomethyl (C)"]);
    let full = generator
        .variable_modifications(&alternatives, &seq("ACAACAACA"), 3, false)
        .unwrap();
    assert_eq!(full.len(), 26);
    let distinct: std::collections::BTreeSet<_> = strings(&full).into_iter().collect();
    assert_eq!(distinct.len(), 26);
    let alternatives = mods(&["Glutathione (C)", "Carbamidomethyl (C)", "Oxidation (M)"]);
    assert_eq!(
        generator
            .variable_modifications(&alternatives, &seq("ACMACMACA"), 3, false)
            .unwrap()
            .len(),
        71
    );
    let alternatives = mods(&["Carbamyl (N-term)", "Oxidation (M)"]);
    let actual = generator
        .variable_modifications(&alternatives, &seq("KAAAAAAAMA"), 2, true)
        .unwrap();
    assert_eq!(
        actual,
        [
            seq("KAAAAAAAMA"),
            seq("KAAAAAAAM(Oxidation)A"),
            seq(".(Carbamyl)KAAAAAAAMA"),
            seq(".(Carbamyl)KAAAAAAAM(Oxidation)A")
        ]
    );
}

#[test]
fn resolution_deduplicates_and_conflicting_fixed_records_have_stable_precedence() {
    let resolved = mods(&["Oxidation (C)", "Carbamidomethyl (C)", "Oxidation (C)"]);
    assert_eq!(resolved.len(), 2);
    assert_eq!(resolved[0].full_id(), "Carbamidomethyl (C)");
    assert!(Generator::get_modifications(&["not-a-modification"]).is_err());
    let mut peptide = seq("CC(Oxidation)");
    Generator::default()
        .apply_fixed_modifications(&resolved, &mut peptide)
        .unwrap();
    assert_eq!(peptide, seq("C(Oxidation)C(Oxidation)"));
    let mut reverse = resolved;
    reverse.reverse();
    let mut again = seq("CC(Oxidation)");
    Generator::default()
        .apply_fixed_modifications(&reverse, &mut again)
        .unwrap();
    assert_eq!(peptide, again);
    let mut peptide = seq("A");
    Generator::default()
        .apply_fixed_modifications(
            &mods(&["Carbamyl (N-term)", "Acetyl (N-term)"]),
            &mut peptide,
        )
        .unwrap();
    assert_eq!(peptide.n_terminal_modification().unwrap().name(), "Acetyl");
}

#[test]
fn source_terminal_fast_path_general_duplicates_and_existing_overwrite_are_retained() {
    let generator = Generator::default();
    let pyro = mods(&["Gln->pyro-Glu (N-term Q)"]);
    let input = seq(".(Acetyl)QA");
    let one = generator
        .variable_modifications(&pyro, &input, 1, false)
        .unwrap();
    assert_eq!(one.len(), 1);
    assert_eq!(one[0].n_terminal_modification().unwrap().name(), "Acetyl");
    assert_eq!(
        one[0].residue_modification(0).unwrap().unwrap().full_id(),
        pyro[0].full_id()
    );
    let general = generator
        .variable_modifications(&pyro, &input, 2, false)
        .unwrap();
    assert_eq!(general.len(), 1);
    assert_eq!(
        general[0].n_terminal_modification().unwrap().full_id(),
        pyro[0].full_id()
    );
    assert!(general[0].residue_modification(0).unwrap().is_none());
    let duplicate = generator
        .variable_modifications(&pyro, &seq("QA"), 2, false)
        .unwrap();
    assert_eq!(duplicate.len(), 2);
    assert_eq!(duplicate[0], duplicate[1]);
    let mut wrong_origin = seq("AA");
    generator
        .apply_fixed_modifications(&pyro, &mut wrong_origin)
        .unwrap();
    assert_eq!(
        wrong_origin.n_terminal_modification().unwrap().full_id(),
        pyro[0].full_id()
    );
    assert!(AASequence::parse(&wrong_origin.to_string()).is_err());
    // Generic terminal origin X has no cached residue in source's fast path.
    assert!(
        generator
            .variable_modifications(&mods(&["Carbamyl (N-term)"]), &seq("X"), 1, false)
            .is_err()
    );
}

#[test]
fn occupied_anonymous_residues_are_protected_but_matching_terminal_rules_still_apply() {
    let generator = Generator::default();
    let input = seq("C[+0.001]MC(Carbamidomethyl)X[123.456789]");
    let original = input.clone();
    let alternatives = mods(&["Carbamidomethyl (C)", "Oxidation (M)"]);
    let generated = generator
        .variable_modifications(&alternatives, &input, 3, false)
        .unwrap();
    assert_eq!(generated.len(), 1);
    assert_eq!(
        generated[0].residue_modification(0).unwrap(),
        original.residue_modification(0).unwrap()
    );
    assert_eq!(
        generated[0].residue_modification(3).unwrap(),
        original.residue_modification(3).unwrap()
    );
    assert!(generated[0].mono_mass().is_ok());
    assert!(generated[0].formula().is_err());
    assert_eq!(input, original);
    let mut input = seq(".[+0.001]QA");
    generator
        .apply_fixed_modifications(&mods(&["Gln->pyro-Glu (N-term Q)"]), &mut input)
        .unwrap();
    assert_eq!(
        input.n_terminal_modification().unwrap().name(),
        "Gln->pyro-Glu"
    );
}

#[test]
fn protein_termini_are_ignored_and_empty_terminal_states_have_source_zero_chemistry() {
    let generator = Generator::default();
    let protein = mods(&["Acetyl (Protein N-term)"]);
    let mut input = seq("MA");
    generator
        .apply_fixed_modifications(&protein, &mut input)
        .unwrap();
    assert!(!input.is_modified());
    assert!(
        generator
            .variable_modifications(&protein, &input, 2, false)
            .unwrap()
            .is_empty()
    );
    let terminal = mods(&["Carbamyl (N-term)"]);
    let empty = AASequence::default();
    assert!(
        generator
            .variable_modifications(&terminal, &empty, 1, false)
            .unwrap()
            .is_empty()
    );
    let generated = generator
        .variable_modifications(&terminal, &empty, 2, true)
        .unwrap();
    assert_eq!(generated.len(), 2);
    assert_eq!(generated[0], empty);
    assert!(generated[1].is_empty());
    assert!(generated[1].is_modified());
    assert_eq!(generated[1].mono_mass().unwrap(), 0.0);
    assert_eq!(generated[1].formula().unwrap().to_string(), "");
    let mut fixed = empty;
    generator
        .apply_fixed_modifications(&terminal, &mut fixed)
        .unwrap();
    assert_eq!(fixed, generated[1]);
}

#[test]
fn append_and_fixed_chemistry_errors_leave_all_original_data_unchanged() {
    let generator = Generator::default();
    let invalid = mods(&["Met->Hse (C-term M)"]);
    // The unchecked source terminal pass removes sulfur from an A-only peptide.
    // Native typed placement is retained where finite valid chemistry exists;
    // negative final atom counts are still errors.
    let mut peptide = seq("AAA");
    let original = peptide.clone();
    assert!(
        generator
            .apply_fixed_modifications(&invalid, &mut peptide)
            .is_err()
    );
    assert_eq!(peptide, original);
    let mut output = vec![seq("PEPTIDE")];
    let saved = output.clone();
    assert!(
        generator
            .apply_variable_modifications(&invalid, &peptide, 2, &mut output, true)
            .is_err()
    );
    assert_eq!(output, saved);
    generator
        .apply_variable_modifications(&mods(&["Oxidation (M)"]), &seq("MM"), 1, &mut output, false)
        .unwrap();
    assert_eq!(output[0], saved[0]);
    assert_eq!(output.len(), 3);
}

#[test]
fn point_site_work_count_and_overflow_guards_are_atomic() {
    let mods = mods(&["Oxidation (M)"]);
    let peptide = seq("MM");
    for generator in [
        Generator {
            max_residues: 1,
            ..Default::default()
        },
        Generator {
            max_sites: 1,
            ..Default::default()
        },
        Generator {
            max_work: 1,
            ..Default::default()
        },
        Generator {
            max_outputs: 3,
            ..Default::default()
        },
        Generator {
            max_output_bytes: 1,
            ..Default::default()
        },
    ] {
        let mut output = vec![seq("A")];
        let original = output.clone();
        assert!(
            generator
                .apply_variable_modifications(&mods, &peptide, 2, &mut output, false)
                .is_err()
        );
        assert_eq!(output, original);
    }
    let mut output = vec![seq("A")];
    Generator {
        max_outputs: 4,
        ..Default::default()
    }
    .apply_variable_modifications(&mods, &peptide, 2, &mut output, false)
    .unwrap();
    assert_eq!(output.len(), 4);
    let generator = Generator {
        max_outputs: usize::MAX,
        max_output_bytes: usize::MAX,
        max_work: usize::MAX,
        ..Default::default()
    };
    assert!(
        generator
            .variable_modifications(&mods, &seq(&"M".repeat(128)), usize::MAX, false)
            .is_err()
    );
}

#[test]
fn payload_limits_cover_owned_anonymous_strings_and_existing_noop_output() {
    let generator = Generator {
        max_output_bytes: 2_000,
        ..Default::default()
    };
    let short = seq("AX[123.456789]K");
    assert_eq!(
        generator
            .variable_modifications(&[], &short, 0, true)
            .unwrap(),
        [short]
    );
    let long = seq(&format!("AX[{}123.456789]K", "0".repeat(900)));
    assert!(
        generator
            .variable_modifications(&[], &long, 0, true)
            .is_err()
    );
    let mut existing = vec![long];
    let saved = existing.clone();
    assert!(
        generator
            .apply_variable_modifications(&[], &seq("A"), 0, &mut existing, false)
            .is_err()
    );
    assert_eq!(existing, saved);
    let mut existing = vec![seq("A")];
    assert!(
        Generator {
            max_outputs: 0,
            ..Default::default()
        }
        .apply_variable_modifications(&[], &seq("A"), 0, &mut existing, false)
        .is_err()
    );
}

fn custom(
    site: &str,
    term: &str,
    delta: f64,
    average_delta: f64,
    formula: &str,
    absolute: f64,
) -> Arc<ResidueModification> {
    let table = format!(
        "90000\tLabMass\tCustom test mass\t{site}\t{term}\t{delta}\t{average_delta}\t{formula}\t0\tOther\t\n"
    );
    let database = ModificationsDB::from_tsv(&table).unwrap();
    let record = database
        .get_modification("LabMass", None, None)
        .unwrap()
        .clone()
        .with_absolute_masses(absolute, absolute)
        .unwrap();
    // The returned shared record outlives this temporary caller registry.
    Arc::new(record)
}

#[test]
fn formula_free_known_records_use_source_free_residue_arithmetic() {
    let water = EmpiricalFormula::parse("H2O").unwrap().mono_mass();
    let free_m = EmpiricalFormula::parse("C5H11NO2S").unwrap().mono_mass();
    for (record, expected_internal) in [
        (
            custom("M", "anywhere", 12.5, 12.6, "", 0.0),
            (free_m + 12.5) - water,
        ),
        (
            custom("M", "anywhere", 12.5, 12.6, "", 250.0),
            250.0 - water,
        ),
    ] {
        let mut peptide = seq("MA");
        Generator::default()
            .apply_fixed_modifications(&[record], &mut peptide)
            .unwrap();
        assert!(peptide.formula().is_err());
        assert!(peptide.average_mass().is_err());
        assert!(
            peptide
                .residue_modification(0)
                .unwrap()
                .unwrap()
                .diff_formula()
                .is_err()
        );
        let first = peptide.fragment_ions(1).unwrap().remove(0);
        assert_eq!(
            first.mz,
            expected_internal + openms::chemistry::PROTON_MASS_U
        );
        assert!(peptide.mz(2).unwrap().is_finite());
        assert!(peptide.prefix(1).unwrap().formula().is_err());
        assert_eq!(peptide.suffix(1).unwrap(), seq("A"));
    }
    // Only a declared absolute supplies an otherwise unknown base mass.
    let mut peptide = seq("X");
    Generator::default()
        .apply_fixed_modifications(
            &[custom("X", "anywhere", 1.0, 0.0, "", 300.0)],
            &mut peptide,
        )
        .unwrap();
    assert_eq!(peptide.mono_mass().unwrap(), 300.0);
    assert!(peptide.formula().is_err());
    let mut unresolved = seq("X");
    Generator::default()
        .apply_fixed_modifications(
            &[custom("X", "anywhere", 1.0, 0.0, "", 0.0)],
            &mut unresolved,
        )
        .unwrap();
    assert!(unresolved.mono_mass().is_err());
}

#[test]
fn formula_precedence_and_absolute_only_noop_follow_source_changes_residue_gate() {
    let generator = Generator::default();
    let original = seq("MA");
    let mut absolute_only = original.clone();
    generator
        .apply_fixed_modifications(
            &[custom("M", "anywhere", 0.0, 0.0, "", 900.0)],
            &mut absolute_only,
        )
        .unwrap();
    assert!(absolute_only.is_modified());
    assert_eq!(
        absolute_only.formula().unwrap(),
        original.formula().unwrap()
    );
    assert_eq!(
        absolute_only.mono_mass().unwrap(),
        original.mono_mass().unwrap()
    );
    assert!(
        absolute_only
            .residue_modification(0)
            .unwrap()
            .unwrap()
            .diff_formula()
            .unwrap()
            .is_empty()
    );
    let mut formula_first = original.clone();
    generator
        .apply_fixed_modifications(
            &[custom("M", "anywhere", 12.5, 12.6, "O", 900.0)],
            &mut formula_first,
        )
        .unwrap();
    assert_eq!(
        formula_first.mono_mass().unwrap(),
        seq("M(Oxidation)A").mono_mass().unwrap()
    );
    assert_eq!(
        formula_first.formula().unwrap(),
        seq("M(Oxidation)A").formula().unwrap()
    );
    let mut charge_only = original.clone();
    generator
        .apply_fixed_modifications(
            &[custom("M", "anywhere", 0.0, 0.0, "+", 900.0)],
            &mut charge_only,
        )
        .unwrap();
    assert_eq!(
        charge_only.mono_mass().unwrap(),
        original.mono_mass().unwrap()
    );
    assert_eq!(charge_only.formula().unwrap(), original.formula().unwrap());
    // A nonzero average difference also activates the source changes-residue
    // gate, so an explicit monoisotopic absolute mass must then take effect.
    let mut average_activates = seq("M");
    generator
        .apply_fixed_modifications(
            &[custom("M", "anywhere", 0.0, 1.0, "", 250.0)],
            &mut average_activates,
        )
        .unwrap();
    assert_eq!(average_activates.mono_mass().unwrap(), 250.0);
    assert!(average_activates.average_mass().is_err());
}

#[test]
fn formula_free_terminal_uses_declared_delta_and_custom_errors_are_atomic() {
    let generator = Generator::default();
    let original = seq("MA");
    let mut terminal = original.clone();
    generator
        .apply_fixed_modifications(
            &[custom("N-term", "n-term", 12.5, 12.6, "", 900.0)],
            &mut terminal,
        )
        .unwrap();
    assert!((terminal.mono_mass().unwrap() - original.mono_mass().unwrap() - 12.5).abs() < 1e-12);
    assert!(terminal.formula().is_err());
    assert!(terminal.average_mass().is_err());
    assert!(
        terminal
            .n_terminal_modification()
            .unwrap()
            .diff_formula()
            .is_err()
    );
    assert_eq!(terminal.suffix(1).unwrap(), seq("A"));
    for record in [
        custom("M", "anywhere", 1.0, 0.0, "", 1.0),
        custom("M", "anywhere", -1000.0, 0.0, "", 0.0),
        custom("N-term", "n-term", -1000.0, 0.0, "", 0.0),
    ] {
        let mut input = original.clone();
        assert!(
            generator
                .apply_fixed_modifications(std::slice::from_ref(&record), &mut input)
                .is_err()
        );
        assert_eq!(input, original);
        let mut output = vec![seq("PEPTIDE")];
        let before = output.clone();
        assert!(
            generator
                .apply_variable_modifications(&[record], &input, 2, &mut output, true)
                .is_err()
        );
        assert_eq!(output, before);
    }
}
