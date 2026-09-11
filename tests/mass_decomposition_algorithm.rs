// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::{
    AASequence, DecompositionResidueSet, MassDecomposition, MassDecompositionAlgorithm,
    MassDecompositionOptions, ModificationRecord, ModificationsDB, ResidueModification,
};

#[test]
fn source_default_and_two_literal_peptide_counts() {
    let defaults = MassDecompositionOptions::default();
    assert_eq!(defaults.decomp_weights_precision, 0.01);
    assert_eq!(defaults.tolerance, 0.3);
    assert_eq!(defaults.residue_set.name(), "Natural19WithoutI");
    let mass = AASequence::parse("DFPIANGER")
        .unwrap()
        .mono_mass_for(openms::chemistry::PeptideFragmentType::Internal, 0)
        .unwrap();
    for (tolerance, count) in [(0.0001, 842), (0.001, 911)] {
        let algorithm = MassDecompositionAlgorithm::from_options(MassDecompositionOptions {
            tolerance,
            ..defaults.clone()
        })
        .unwrap();
        let result = algorithm.decompositions(mass).unwrap();
        assert_eq!(result.len(), count);
        assert!(
            result
                .iter()
                .any(|value| value.equals_text("A1 D1 E1 F1 G1 L1 N1 P1 R1").unwrap())
        );
    }
}

fn record(name: &str, origin: char, absolute: f64, delta: f64) -> ResidueModification {
    ResidueModification::from_record(ModificationRecord {
        name: name.into(),
        origin: Some(origin),
        mono_mass: absolute,
        diff_mono_mass: delta,
        ..Default::default()
    })
    .unwrap()
}
fn mass(algorithm: &MassDecompositionAlgorithm, symbol: char) -> f64 {
    algorithm
        .alphabet()
        .iter()
        .find(|&&(letter, _)| letter == symbol)
        .unwrap()
        .1
}

#[test]
fn all_residue_sets_keep_source_membership_and_ambiguous_limits() {
    for set in DecompositionResidueSet::ALL {
        assert_eq!(set.name().parse::<DecompositionResidueSet>().unwrap(), set);
        let result = MassDecompositionAlgorithm::from_options(MassDecompositionOptions {
            residue_set: set,
            ..Default::default()
        });
        let count = match set {
            DecompositionResidueSet::All
            | DecompositionResidueSet::Ambiguous
            | DecompositionResidueSet::AmbiguousWithoutX => {
                assert!(result.is_err());
                continue;
            }
            DecompositionResidueSet::AllNatural => 22,
            DecompositionResidueSet::Natural19WithoutI
            | DecompositionResidueSet::Natural19WithoutL => 19,
            _ => 20,
        };
        assert_eq!(result.unwrap().alphabet().len(), count);
    }
    assert!("natural20".parse::<DecompositionResidueSet>().is_err());
    let registry = ModificationsDB::from_records(vec![
        record("repairB", 'B', 100.0, 0.0),
        record("repairZ", 'Z', 101.0, 0.0),
    ])
    .unwrap();
    let algorithm = MassDecompositionAlgorithm::with_registry(
        MassDecompositionOptions {
            residue_set: DecompositionResidueSet::AmbiguousWithoutX,
            fixed_modifications: vec!["repairB".into(), "repairZ".into()],
            ..Default::default()
        },
        &registry,
    )
    .unwrap();
    assert_eq!(algorithm.alphabet().len(), 25);
    assert_eq!(mass(&algorithm, 'B'), 100.0);
    assert_eq!(mass(&algorithm, 'J'), mass(&algorithm, 'I'));
}

#[test]
fn custom_absolute_delta_ignored_labels_and_fixed_adjusted_variable_masses() {
    let base = MassDecompositionAlgorithm::new().unwrap();
    let registry = ModificationsDB::from_records(vec![
        record("a-absolute", 'A', 200.0, 5.0),
        record("b-delta", 'A', 0.0, 3.0),
        record("c-variable", 'A', 0.0, 10.0),
        record("a-ignore", 'X', 12.0, 5.0),
        record("b-empty", 'A', 0.0, 0.0),
    ])
    .unwrap();
    let algorithm = MassDecompositionAlgorithm::with_registry(
        MassDecompositionOptions {
            fixed_modifications: vec!["b-delta".into(), "a-absolute".into(), "b-delta".into()],
            variable_modifications: vec!["c-variable".into(), "b-empty".into(), "a-ignore".into()],
            ..Default::default()
        },
        &registry,
    )
    .unwrap();
    assert_eq!(mass(&algorithm, 'A'), 203.0); // Full-ID order; direct absolute, no H2O subtraction.
    assert_eq!(mass(&algorithm, 'c'), 213.0); // Ignored definitions consume a and b.
    assert!(
        !algorithm
            .alphabet()
            .iter()
            .any(|&(c, _)| c == 'a' || c == 'b')
    );
    assert_eq!(algorithm.diagnostics().len(), 2);
    assert_eq!(mass(&algorithm, 'G'), mass(&base, 'G'));
    drop(registry);
    assert_eq!(mass(&algorithm, 'c'), 213.0);
}

#[test]
fn source_absolute_override_needs_no_delta_and_terminal_specificity_is_ignored() {
    let mut record = ModificationRecord {
        name: "end".into(),
        origin: Some('A'),
        mono_mass: 123.0,
        term_specificity: openms::chemistry::TermSpecificity::NTerm,
        ..Default::default()
    };
    let first = ResidueModification::from_record(record.clone()).unwrap();
    record.name = "missingI".into();
    record.origin = Some('I');
    record.mono_mass = 0.0;
    record.diff_mono_mass = 77.0;
    let registry = ModificationsDB::from_records(vec![
        first,
        ResidueModification::from_record(record).unwrap(),
    ])
    .unwrap();
    let algorithm = MassDecompositionAlgorithm::with_registry(
        MassDecompositionOptions {
            fixed_modifications: vec!["end".into(), "missingI".into()],
            ..Default::default()
        },
        &registry,
    )
    .unwrap();
    assert_eq!(mass(&algorithm, 'A'), 123.0);
    assert_eq!(mass(&algorithm, 'I'), 77.0); // Missing map key starts at zero.
}

#[test]
fn failed_reconfiguration_and_append_are_atomic() {
    let mut algorithm = MassDecompositionAlgorithm::new().unwrap();
    let before = algorithm.alphabet().to_vec();
    for precision in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e-12, 1e20] {
        assert!(
            algorithm
                .set_options(MassDecompositionOptions {
                    decomp_weights_precision: precision,
                    ..Default::default()
                })
                .is_err()
        );
        assert_eq!(algorithm.alphabet(), before);
    }
    let mut output = vec![MassDecomposition::parse("Z1").unwrap()];
    let before_output = output.clone();
    let pointer = output.as_ptr();
    for query in [-1.0, f64::NAN, f64::INFINITY, 0.0, 1e30, 1e8] {
        assert!(algorithm.append_decompositions(&mut output, query).is_err());
        assert_eq!(output, before_output);
        assert_eq!(output.as_ptr(), pointer);
    }
    algorithm
        .append_decompositions(&mut output, mass(&algorithm, 'G'))
        .unwrap();
    assert_eq!(output[0], before_output[0]);
    assert!(
        output[1..]
            .iter()
            .any(|value| value.equals_text("G1").unwrap())
    );
}

#[test]
fn modification_limits_and_nonpositive_resolved_weights_are_checked() {
    let records: Vec<_> = (0..27)
        .map(|i| record(&format!("ignored{i:02}"), 'X', 1.0, 0.0))
        .collect();
    let registry = ModificationsDB::from_records(records).unwrap();
    let names = registry
        .entries()
        .iter()
        .map(|r| r.full_id().to_owned())
        .collect();
    assert!(
        MassDecompositionAlgorithm::with_registry(
            MassDecompositionOptions {
                variable_modifications: names,
                ..Default::default()
            },
            &registry
        )
        .is_err()
    );
    for absolute in [-1.0, 1e-9, f64::MAX] {
        let registry =
            ModificationsDB::from_records(vec![record("invalidA", 'A', absolute, 0.0)]).unwrap();
        assert!(
            MassDecompositionAlgorithm::with_registry(
                MassDecompositionOptions {
                    fixed_modifications: vec!["invalidA".into()],
                    ..Default::default()
                },
                &registry
            )
            .is_err()
        );
    }
    for tolerance in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            MassDecompositionAlgorithm::from_options(MassDecompositionOptions {
                tolerance,
                ..Default::default()
            })
            .is_err()
        );
    }
    let options = MassDecompositionOptions {
        fixed_modifications: vec!["missing".into(); 10_001],
        ..Default::default()
    };
    assert!(MassDecompositionAlgorithm::from_options(options).is_err());
    let options = MassDecompositionOptions {
        fixed_modifications: vec!["n".repeat(1024 * 1024 + 1)],
        ..Default::default()
    };
    assert!(MassDecompositionAlgorithm::from_options(options).is_err());
}

#[test]
fn ambiguous_lookup_uses_first_provider_and_full_id_set_identity() {
    let records = vec![
        record("same", 'A', 200.0, 0.0),
        record("same", 'C', 250.0, 0.0),
    ];
    let registry = ModificationsDB::from_records(records).unwrap();
    assert!(registry.get_modification("same", None, None).is_err()); // Global lookup policy stays unchanged.
    let options = MassDecompositionOptions {
        fixed_modifications: vec!["same".into(), "same (A)".into()],
        ..Default::default()
    };
    let algorithm = MassDecompositionAlgorithm::with_registry(options, &registry).unwrap();
    assert_eq!(mass(&algorithm, 'A'), 200.0);
    assert_ne!(mass(&algorithm, 'C'), 250.0);
    let registry = ModificationsDB::from_records(vec![
        record("same", 'A', 200.0, 0.0),
        record("same", 'A', 250.0, 0.0),
    ])
    .unwrap();
    let algorithm = MassDecompositionAlgorithm::with_registry(
        MassDecompositionOptions {
            fixed_modifications: vec!["same".into()],
            ..Default::default()
        },
        &registry,
    )
    .unwrap();
    assert_eq!(mass(&algorithm, 'A'), 200.0);
}
