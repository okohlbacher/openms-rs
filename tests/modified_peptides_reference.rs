// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Independent source assertions and branch-specific reference cases.
//! Source hashes and the distinction between literal/derived cases are recorded
//! in data/modification_generation_provenance.json and the reference review.

use openms::chemistry::{
    AASequence, IonSeries, ModificationsDB, ModifiedPeptideGenerator, ResidueModification,
};
use std::sync::Arc;

fn peptide(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}
fn mods(text: &str) -> Vec<Arc<ResidueModification>> {
    let names: Vec<_> = text.split(';').filter(|s| !s.is_empty()).collect();
    ModifiedPeptideGenerator::get_modifications(&names).unwrap()
}
fn terminal_name(peptide: &AASequence) -> Option<&str> {
    peptide.n_terminal_modification().map(|m| m.full_id())
}

#[test]
fn seven_fixed_modification_outputs_match_literal_upstream_assertions() {
    let generator = ModifiedPeptideGenerator::default();
    for row in include_str!("data/modification_generation_fixed.tsv")
        .lines()
        .skip(1)
    {
        let fields: Vec<_> = row.split('\t').collect();
        let mut input = peptide(fields[2]);
        generator
            .apply_fixed_modifications(&mods(fields[1]), &mut input)
            .unwrap();
        assert_eq!(input, peptide(fields[3]), "{}", fields[0]);
    }
}

#[test]
fn nineteen_variable_cases_match_source_counts_and_specified_order() {
    let generator = ModifiedPeptideGenerator::default();
    for row in include_str!("data/modification_generation_variable.tsv")
        .lines()
        .skip(1)
    {
        let fields: Vec<_> = row.split('\t').collect();
        let input = peptide(fields[2]);
        let output = generator
            .variable_modifications(
                &mods(fields[1]),
                &input,
                fields[3].parse().unwrap(),
                fields[4].parse().unwrap(),
            )
            .unwrap();
        assert_eq!(
            output.len(),
            fields[5].parse::<usize>().unwrap(),
            "{}",
            fields[0]
        );
        let expected: Vec<_> = include_str!("data/modification_generation_expected.tsv")
            .lines()
            .skip(1)
            .filter_map(|line| {
                let v: Vec<_> = line.split('\t').collect();
                (v[0] == fields[0]).then(|| peptide(v[2]))
            })
            .collect();
        match fields[6] {
            "ordered" => assert_eq!(output, expected, "{}", fields[0]),
            "set" => {
                let mut actual: Vec<_> = output.iter().map(ToString::to_string).collect();
                let mut expected: Vec<_> = expected.iter().map(ToString::to_string).collect();
                actual.sort();
                expected.sort();
                assert_eq!(actual, expected, "{}", fields[0]);
            }
            "count" => {}
            _ => panic!("unknown fixture comparison"),
        }
    }
}

#[test]
fn terminal_fast_path_general_path_and_duplicate_placements_remain_distinct() {
    let generator = ModifiedPeptideGenerator::default();
    let pyro = mods("Gln->pyro-Glu (N-term Q)");
    let input = peptide("QAA");
    let one = generator
        .variable_modifications(&pyro, &input, 1, false)
        .unwrap();
    assert_eq!(one.len(), 1);
    assert!(one[0].n_terminal_modification().is_none());
    assert_eq!(
        one[0].residue_modification(0).unwrap().unwrap().full_id(),
        "Gln->pyro-Glu (N-term Q)"
    );
    let general = generator
        .variable_modifications(&pyro, &input, 2, false)
        .unwrap();
    // The general path inserts this terminal alternative twice: once in its
    // terminal prepass and once for the matching first Q residue.
    assert_eq!(general.len(), 2);
    assert_eq!(general[0], general[1]);
    assert_eq!(terminal_name(&general[0]), Some("Gln->pyro-Glu (N-term Q)"));
    assert!(general[0].residue_modification(0).unwrap().is_none());
    assert_eq!(one[0].formula().unwrap(), general[0].formula().unwrap());
    assert!(one[0].mono_mass().unwrap().is_finite());
    let original_ions = input.fragment_ions(2).unwrap();
    let residue_ions = one[0].fragment_ions(2).unwrap();
    let terminal_ions = general[0].fragment_ions(2).unwrap();
    for ((original, residue), terminal) in
        original_ions.iter().zip(&residue_ions).zip(&terminal_ions)
    {
        // Residue::setFormula recomputes mass from elemental composition. The
        // true terminal slot instead uses the declared rounded mass delta.
        // Complementary y ions lacking the Q change in neither path.
        let (residue_delta, terminal_delta) = if original.series == IonSeries::B {
            (
                pyro[0].diff_formula().mono_mass() / f64::from(original.charge),
                pyro[0].diff_mono_mass() / f64::from(original.charge),
            )
        } else {
            (0., 0.)
        };
        assert!((residue.mz - (original.mz + residue_delta)).abs() < 1e-9);
        assert!((terminal.mz - (original.mz + terminal_delta)).abs() < 1e-9);
    }

    let generic = mods("Carbamyl (N-term)");
    assert!(
        generator
            .variable_modifications(&generic, &peptide("AAA"), 1, false)
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        generator
            .variable_modifications(&generic, &peptide("AAA"), 2, false)
            .unwrap()
            .len(),
        1
    );
    let protein_terminal = mods("Deamidated (Protein N-term F)");
    assert!(
        generator
            .variable_modifications(&protein_terminal, &peptide("FAA"), 2, false)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn existing_and_foreign_origin_terminal_assignments_follow_source_pass_order() {
    let generator = ModifiedPeptideGenerator::default();
    let pyro = mods("Gln->pyro-Glu (N-term Q)");
    let input = peptide(".(Acetyl)QAA");
    let one = generator
        .variable_modifications(&pyro, &input, 1, false)
        .unwrap();
    assert_eq!(one.len(), 1);
    assert_eq!(terminal_name(&one[0]), Some("Acetyl (N-term)"));
    assert!(one[0].residue_modification(0).unwrap().is_some());
    let general = generator
        .variable_modifications(&pyro, &input, 2, false)
        .unwrap();
    assert_eq!(general.len(), 1);
    assert_eq!(terminal_name(&general[0]), Some("Gln->pyro-Glu (N-term Q)"));
    assert!(general[0].residue_modification(0).unwrap().is_none());

    // Source's terminal prepass does not inspect the first residue's origin.
    let mut foreign = peptide("AAA");
    generator
        .apply_fixed_modifications(&pyro, &mut foreign)
        .unwrap();
    assert_eq!(terminal_name(&foreign), Some("Gln->pyro-Glu (N-term Q)"));
    assert!(foreign.formula().unwrap().mono_mass().is_finite());
    assert!(AASequence::parse(&foreign.to_string()).is_err());
}

#[test]
fn fixed_conflicts_use_full_id_order_and_skip_preexisting_residue_annotations() {
    let generator = ModifiedPeptideGenerator::default();
    let alternatives = mods("Glutathione (C);Carbamidomethyl (C);Glutathione (C)");
    assert_eq!(
        alternatives.iter().map(|m| m.full_id()).collect::<Vec<_>>(),
        ["Carbamidomethyl (C)", "Glutathione (C)"]
    );
    let mut sequence = peptide("CC(Oxidation)C[999]");
    let original = sequence.clone();
    generator
        .apply_fixed_modifications(&alternatives, &mut sequence)
        .unwrap();
    assert_eq!(
        sequence.residue_modification(0).unwrap().unwrap().full_id(),
        "Glutathione (C)"
    );
    assert_eq!(
        sequence.residue_modification(1).unwrap(),
        original.residue_modification(1).unwrap()
    );
    assert_eq!(
        sequence.residue_modification(2).unwrap(),
        original.residue_modification(2).unwrap()
    );
}

#[test]
fn weighted_sites_have_all_eight_distinct_variants_in_source_group_order() {
    // Independent two-site Cartesian oracle: right site, left site, both sites.
    // Full-ID order replaces unspecified source pointer order within each site;
    // the earlier (right) site varies fastest in the two-site group.
    let result = ModifiedPeptideGenerator::default()
        .variable_modifications(
            &mods("Glutathione (C);Carbamidomethyl (C)"),
            &peptide("CC"),
            2,
            false,
        )
        .unwrap();
    let expected = [
        "CC(Carbamidomethyl)",
        "CC(Glutathione)",
        "C(Carbamidomethyl)C",
        "C(Glutathione)C",
        "C(Carbamidomethyl)C(Carbamidomethyl)",
        "C(Carbamidomethyl)C(Glutathione)",
        "C(Glutathione)C(Carbamidomethyl)",
        "C(Glutathione)C(Glutathione)",
    ]
    .map(peptide);
    assert_eq!(result, expected);
}

#[test]
fn custom_mass_only_records_follow_residue_override_and_terminal_delta_branches() {
    // Shared test records outlive their caller registries without leaked storage.
    fn custom(row: &str, absolute: f64) -> Arc<ResidueModification> {
        let db = ModificationsDB::from_tsv(row).unwrap();
        Arc::new(
            db.entries()[0]
                .as_ref()
                .clone()
                .with_absolute_masses(absolute, 0.)
                .unwrap(),
        )
    }
    fn apply(input: &str, modification: Arc<ResidueModification>) -> AASequence {
        let mut input = peptide(input);
        ModifiedPeptideGenerator::default()
            .apply_fixed_modifications(&[modification], &mut input)
            .unwrap();
        input
    }
    fn mass_only(sequence: &AASequence, expected_mass: f64) {
        assert!((sequence.mono_mass().unwrap() - expected_mass).abs() < 1e-9);
        assert!(sequence.formula().is_err());
        assert!(sequence.average_mass().is_err());
    }

    let base = peptide("M");
    let water_mass = "H2O"
        .parse::<openms::chemistry::EmpiricalFormula>()
        .unwrap()
        .mono_mass();
    let delta = custom(
        "1\tDeltaOnly\tDelta-only reference\tM\tanywhere\t10.25\t11\t\t0\tOther\t",
        0.,
    );
    let changed = apply("M", delta);
    let source_delta_mass = water_mass + (base.formula().unwrap().mono_mass() + 10.25 - water_mass);
    mass_only(&changed, source_delta_mass);
    assert_eq!(
        changed.mono_mass().unwrap().to_bits(),
        source_delta_mass.to_bits()
    );
    assert!(
        changed
            .residue_modification(0)
            .unwrap()
            .unwrap()
            .diff_formula()
            .is_err()
    );

    // Source adopts a nonzero absolute FREE-residue mass when no formula is
    // available, then subtracts water for the internal residue representation.
    let absolute = custom(
        "2\tAbsolute\tAbsolute reference\tM\tanywhere\t10.25\t11\t\t0\tOther\t",
        300.,
    );
    let changed = apply("MA", absolute);
    let alanine_internal = "C3H5NO"
        .parse::<openms::chemistry::EmpiricalFormula>()
        .unwrap()
        .mono_mass();
    mass_only(&changed, 300. + alanine_internal);
    let ions = changed.fragment_ions(1).unwrap();
    let expected_b1 = 300. - water_mass + openms::chemistry::PROTON_MASS_U;
    assert!((ions[0].mz - expected_b1).abs() < 1e-9);
    assert_eq!(ions[1], peptide("MA").fragment_ions(1).unwrap()[1]);

    let unknown = custom(
        "3\tUnknownAbsolute\tUnknown absolute reference\tX\tanywhere\t1\t0\t\t0\tOther\t",
        300.,
    );
    mass_only(&apply("X", unknown), 300.);

    // Residue.cpp's changes_residue predicate ignores absolute-only records;
    // applying an annotation does not adopt its unrelated absolute mass.
    let no_op = custom(
        "4\tNoChange\tAbsolute-only reference\tM\tanywhere\t0\t0\t\t0\tOther\t",
        300.,
    );
    let unchanged = apply("M", no_op);
    assert!(unchanged.is_modified());
    assert_eq!(unchanged.mono_mass().unwrap(), base.mono_mass().unwrap());
    assert_eq!(unchanged.formula().unwrap(), base.formula().unwrap());
    assert_eq!(
        unchanged.average_mass().unwrap(),
        base.average_mass().unwrap()
    );

    // A supplied delta formula wins over both absolute and declared delta mass.
    let formula = custom(
        "5\tFormulaWins\tFormula precedence reference\tM\tanywhere\t10.25\t11\tO1\t0\tOther\t",
        300.,
    );
    let oxidized = apply("M", formula);
    assert_eq!(
        oxidized.formula().unwrap(),
        peptide("M(Oxidation)").formula().unwrap()
    );
    assert!(
        (oxidized.mono_mass().unwrap() - peptide("M(Oxidation)").mono_mass().unwrap()).abs() < 1e-9
    );

    let average_only = custom(
        "6\tAverageOnly\tAverage-only reference\tM\tanywhere\t0\t2\t\t0\tOther\t",
        0.,
    );
    mass_only(&apply("M", average_only), base.mono_mass().unwrap());

    // True terminal slots use declared differences, never the absolute override.
    let terminal = custom(
        "7\tTerminalMass\tTerminal reference\tN-term\tn-term\t12.5\t13\t\t0\tOther\t",
        300.,
    );
    let changed = apply("M", terminal);
    assert!(changed.n_terminal_modification().is_some());
    mass_only(&changed, base.mono_mass().unwrap() + 12.5);

    let invalid = custom(
        "8\tNegativeAbsolute\tNegative absolute reference\tM\tanywhere\t1\t0\t\t0\tOther\t",
        -1.,
    );
    let mut unchanged = base.clone();
    assert!(
        ModifiedPeptideGenerator::default()
            .apply_fixed_modifications(&[invalid], &mut unchanged)
            .is_err()
    );
    assert_eq!(unchanged, base);
}

#[cfg(feature = "idxml")]
#[test]
fn source_empty_annotated_sequence_is_rejected_before_xml_output() {
    use openms::format::idxml::{self, IdXmlDocument};
    use openms::identification::{PeptideHit, PeptideIdentification, ProteinIdentification};

    let mut empty = peptide("");
    ModifiedPeptideGenerator::default()
        .apply_fixed_modifications(&mods("Acetyl (N-term)"), &mut empty)
        .unwrap();
    assert!(empty.is_empty() && empty.is_modified());
    assert_eq!(terminal_name(&empty), Some("Acetyl (N-term)"));
    // AASequence.cpp enters formula/mass accumulation only for nonempty data.
    assert_eq!(empty.mono_mass().unwrap(), 0.);
    assert_eq!(empty.formula().unwrap().to_string(), "");
    let doc = IdXmlDocument {
        protein_identifications: vec![ProteinIdentification {
            identifier: "empty-sequence-reference".into(),
            date_time: Some("2026-09-10T12:00:00".into()),
            ..Default::default()
        }],
        peptide_identifications: vec![PeptideIdentification {
            identifier: "empty-sequence-reference".into(),
            hits: vec![PeptideHit {
                sequence: empty,
                ..Default::default()
            }],
            ..Default::default()
        }],
        ..Default::default()
    };
    let before = doc.clone();
    let mut output = b"existing caller output".to_vec();
    let original = output.clone();
    let error = idxml::write(&mut output, &doc).unwrap_err();
    assert!(
        error.to_string().contains("modification placement"),
        "{error}"
    );
    assert_eq!(output, original);
    assert_eq!(doc, before);
}

#[test]
fn append_limits_account_for_existing_results_including_noop_generation() {
    let generator = ModifiedPeptideGenerator::default();
    let input = peptide("MM");
    let ox = mods("Oxidation (M)");
    let mut output = vec![peptide("A")];
    generator
        .apply_variable_modifications(&ox, &input, 1, &mut output, false)
        .unwrap();
    assert_eq!(
        output,
        [
            peptide("A"),
            peptide("MM(Oxidation)"),
            peptide("M(Oxidation)M")
        ]
    );

    let limited = ModifiedPeptideGenerator {
        max_outputs: 2,
        ..Default::default()
    };
    let mut output = vec![peptide("A")];
    let before = output.clone();
    assert!(
        limited
            .apply_variable_modifications(&ox, &input, 1, &mut output, false)
            .is_err()
    );
    assert_eq!(output, before);
    output.push(peptide("AA"));
    output.push(peptide("AAA"));
    let before = output.clone();
    assert!(
        limited
            .apply_variable_modifications(&[], &input, 0, &mut output, false)
            .is_err()
    );
    assert_eq!(output, before);

    let limited = ModifiedPeptideGenerator {
        max_output_bytes: 1,
        ..Default::default()
    };
    let mut output = vec![peptide("AAAA")];
    let before = output.clone();
    assert!(
        limited
            .apply_variable_modifications(&[], &input, 0, &mut output, false)
            .is_err()
    );
    assert_eq!(output, before);
}
