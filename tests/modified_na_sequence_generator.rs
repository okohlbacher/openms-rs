// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::{
    ModifiedNASequenceGenerator, NAFragmentType, NASequence, Ribonucleotide, RibonucleotideDB,
    RibonucleotideRecord, RibonucleotideTermSpecificity as Term,
};
use std::sync::Arc;

fn r(code: &str, origin: char, term: Term, formula: &str) -> Arc<Ribonucleotide> {
    Arc::new(
        Ribonucleotide::from_record(RibonucleotideRecord {
            code: code.into(),
            origin,
            term_specificity: term,
            formula: formula.parse().unwrap(),
            ..Default::default()
        })
        .unwrap(),
    )
}
fn mods(names: &[&str]) -> Vec<Arc<Ribonucleotide>> {
    names
        .iter()
        .map(|name| RibonucleotideDB::global().get(name).unwrap())
        .collect()
}
fn strings(sequences: &[NASequence]) -> Vec<String> {
    sequences.iter().map(ToString::to_string).collect()
}

#[test]
fn source_fixed_literal_and_existing_modification_protection() {
    // Pinned ModifiedNASequenceGenerator_test.cpp:28–44.
    let generator = ModifiedNASequenceGenerator::default();
    let mut sequence: NASequence = "AUAUAUA".parse().unwrap();
    generator
        .apply_fixed_modifications(&mods(&["s4U"]), &mut sequence)
        .unwrap();
    assert_eq!(sequence.to_string(), "A[s4U]A[s4U]A[s4U]A");
    assert_eq!(sequence, "A[s4U]A[s4U]A[s4U]A".parse().unwrap());
    let original = sequence.clone();
    generator
        .apply_fixed_modifications(&mods(&["m3U"]), &mut sequence)
        .unwrap();
    assert_eq!(sequence, original);
}

#[test]
fn source_variable_counts_and_caller_order() {
    // The four literal count cases at source test lines 50–91.
    let generator = ModifiedNASequenceGenerator::default();
    let sequence = "AUAUAUA".parse().unwrap();
    let modifications = mods(&["s4U", "m3U"]);
    let variants = generator
        .variable_modifications(&modifications, &sequence, 1, true)
        .unwrap();
    assert_eq!(
        strings(&variants),
        [
            "AUAUAUA",
            "AUAUA[s4U]A",
            "AUAUA[m3U]A",
            "AUA[s4U]AUA",
            "AUA[m3U]AUA",
            "A[s4U]AUAUA",
            "A[m3U]AUAUA",
        ]
    );
    assert_eq!(
        generator
            .variable_modifications(&modifications, &sequence, 1, false)
            .unwrap()
            .len(),
        6
    );
    assert_eq!(
        generator
            .variable_modifications(&modifications, &sequence, 3, true)
            .unwrap()
            .len(),
        27
    );
    assert_eq!(
        generator
            .variable_modifications(&mods(&["s4U", "m3U", "m1A"]), &sequence, 7, true)
            .unwrap()
            .len(),
        432
    );
    // Derived coefficients of (1+2x)^3, including the unchanged sequence.
    assert_eq!(
        generator
            .variable_modifications(&modifications, &sequence, 2, true)
            .unwrap()
            .len(),
        19
    );
    assert_eq!(
        generator
            .variable_modifications(&modifications, &sequence, usize::MAX, true)
            .unwrap()
            .len(),
        27
    );
}

#[test]
fn exact_subset_order_and_increasing_placement_count() {
    let generator = ModifiedNASequenceGenerator::default();
    let sequence = "AAAA".parse().unwrap();
    let modifications = [r("a", 'A', Term::Anywhere, "")];
    let variants = generator
        .variable_modifications(&modifications, &sequence, 2, false)
        .unwrap();
    assert_eq!(
        strings(&variants),
        [
            "AAAa", "AAaA", "AaAA", "aAAA", "AAaa", "AaAa", "AaaA", "aAAa", "aAaA", "aaAA"
        ]
    );
}

#[test]
fn fixed_first_terminal_last_residue_and_original_slot_protection() {
    let generator = ModifiedNASequenceGenerator::default();
    let t1 = r("first3", 'X', Term::ThreePrime, "H");
    let t2 = r("second3", 'U', Term::ThreePrime, "O");
    let f1 = r("first5", 'X', Term::FivePrime, "H");
    let f2 = r("second5", 'U', Term::FivePrime, "O");
    let a = r("firstU", 'U', Term::Anywhere, "H");
    let b = r("lastU", 'U', Term::Anywhere, "O");
    let mut sequence = "UU".parse().unwrap();
    generator
        .apply_fixed_modifications(
            &[t1.clone(), t2, f1.clone(), f2, a, b.clone()],
            &mut sequence,
        )
        .unwrap();
    assert!(Arc::ptr_eq(sequence.three_prime_mod().unwrap(), &t1));
    assert!(Arc::ptr_eq(sequence.five_prime_mod().unwrap(), &f1));
    assert!(
        sequence
            .residues()
            .iter()
            .all(|record| Arc::ptr_eq(record, &b))
    );
    let before = sequence.clone();
    generator
        .apply_fixed_modifications(&mods(&["s4U", "5'-p", "3'-p"]), &mut sequence)
        .unwrap();
    assert_eq!(sequence, before);
}

#[test]
fn maximum_one_terminal_specificity_quirk_is_a_residue_replacement() {
    let generator = ModifiedNASequenceGenerator::default();
    let terminal = r("terminalU", 'U', Term::FivePrime, "H");
    let registry = RibonucleotideDB::from_entries(vec![
        openms::chemistry::RibonucleotideEntry {
            ribonucleotide: terminal.clone(),
            alternatives: None,
        },
        openms::chemistry::RibonucleotideEntry {
            ribonucleotide: mods(&["U"]).remove(0),
            alternatives: None,
        },
    ])
    .unwrap();
    let sequence = NASequence::parse_with_registry("U", &registry).unwrap();
    let single = generator
        .variable_modifications(std::slice::from_ref(&terminal), &sequence, 1, false)
        .unwrap();
    assert_eq!(single.len(), 1);
    assert!(Arc::ptr_eq(&single[0].residues()[0], &terminal));
    assert!(single[0].five_prime_mod().is_none());
    assert!(single[0].checked_string_with_registry(&registry).is_err());
    let general = generator
        .variable_modifications(std::slice::from_ref(&terminal), &sequence, 2, false)
        .unwrap();
    assert_eq!(general.len(), 1);
    assert!(Arc::ptr_eq(general[0].five_prime_mod().unwrap(), &terminal));
    assert_eq!(general[0].residues(), sequence.residues());
    let mut occupied = sequence;
    occupied
        .set_five_prime_mod(Some(r("existing", '.', Term::Anywhere, "O")))
        .unwrap();
    assert_eq!(
        generator
            .variable_modifications(std::slice::from_ref(&terminal), &occupied, 1, false)
            .unwrap()
            .len(),
        1
    );
    assert!(
        generator
            .variable_modifications(&[terminal], &occupied, 2, false)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn empty_sequences_terminal_sites_and_source_iteration_order() {
    let generator = ModifiedNASequenceGenerator::default();
    let modifications = [
        r("five", 'Q', Term::FivePrime, "H"),
        r("three", 'X', Term::ThreePrime, "O"),
    ];
    let empty = NASequence::new();
    assert!(
        generator
            .variable_modifications(&modifications, &empty, 1, false)
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        generator
            .variable_modifications(&modifications, &empty, 1, true)
            .unwrap()
            .as_slice(),
        std::slice::from_ref(&empty)
    );
    let variants = generator
        .variable_modifications(&modifications, &empty, 2, true)
        .unwrap();
    assert_eq!(
        strings(&variants),
        ["", "[five]", "[three]", "[five][three]"]
    );
    assert!(variants.iter().all(NASequence::is_empty));
    for variant in &variants {
        assert!(variant.formula(NAFragmentType::Full, 0).unwrap().is_empty());
        assert_eq!(variant.mono_mass(NAFragmentType::Full, 0).unwrap(), 0.0);
    }
    let mut fixed = empty;
    generator
        .apply_fixed_modifications(&modifications, &mut fixed)
        .unwrap();
    assert_eq!(fixed, variants[3]);
}

#[test]
fn repeated_arc_dedup_preserves_distinct_allocations_and_chemistry() {
    let generator = ModifiedNASequenceGenerator::default();
    let sequence = "U".parse().unwrap();
    let first = r("same", 'U', Term::Anywhere, "H");
    let equal_copy = Arc::new(first.as_ref().clone());
    let different = r("same", 'U', Term::Anywhere, "O");
    let variants = generator
        .variable_modifications(
            &[
                first.clone(),
                first.clone(),
                equal_copy.clone(),
                different.clone(),
                equal_copy.clone(),
            ],
            &sequence,
            2,
            false,
        )
        .unwrap();
    assert_eq!(variants.len(), 3);
    assert!(Arc::ptr_eq(&variants[0].residues()[0], &first));
    assert!(Arc::ptr_eq(&variants[1].residues()[0], &equal_copy));
    assert!(Arc::ptr_eq(&variants[2].residues()[0], &different));
    assert_eq!(variants[0], variants[1]);
    assert_ne!(variants[0], variants[2]);
    assert_ne!(
        variants[0].formula(NAFragmentType::Full, 0).unwrap(),
        variants[2].formula(NAFragmentType::Full, 0).unwrap()
    );
}

#[test]
fn no_op_candidates_and_ambiguous_records_are_not_removed_or_expanded() {
    let generator = ModifiedNASequenceGenerator::default();
    let sequence: NASequence = "UU".parse().unwrap();
    let variants = generator
        .variable_modifications(&mods(&["U"]), &sequence, 2, true)
        .unwrap();
    assert_eq!(variants.len(), 4);
    assert!(variants.iter().all(|variant| variant == &sequence));
    let ambiguous = r("choice?", 'U', Term::Anywhere, "H-1");
    let variants = generator
        .variable_modifications(&[ambiguous], &sequence, 1, false)
        .unwrap();
    assert_eq!(variants.len(), 2);
    assert_eq!(strings(&variants), ["U[choice?]", "[choice?]U"]);
}

#[test]
fn generator_accepts_signed_records_without_recalculating_chemistry() {
    let record = Arc::new(
        Ribonucleotide::from_record(RibonucleotideRecord {
            code: "negative".into(),
            origin: 'U',
            mono_mass: -123.0,
            average_mass: -456.0,
            formula: "H-10".parse().unwrap(),
            ..Default::default()
        })
        .unwrap(),
    );
    let mut sequence = "U".parse().unwrap();
    ModifiedNASequenceGenerator::default()
        .apply_fixed_modifications(std::slice::from_ref(&record), &mut sequence)
        .unwrap();
    assert!(Arc::ptr_eq(&sequence.residues()[0], &record));
    assert_eq!(sequence.residues()[0].mono_mass(), -123.0);
    assert!(sequence.mono_mass(NAFragmentType::Full, 0).unwrap() < 0.0);
}

#[test]
fn owned_sulfur_context_survives_and_missing_context_is_not_replaced() {
    let generator = ModifiedNASequenceGenerator::default();
    let custom_end = r("5'-p*", 'X', Term::FivePrime, "H3");
    let weak = Arc::downgrade(&custom_end);
    let registry = RibonucleotideDB::from_records(vec![custom_end.as_ref().clone()]).unwrap();
    let context = registry.get("5'-p*").unwrap();
    let sequence = NASequence::from_records_with_registry(mods(&["A", "U"]), &registry).unwrap();
    let variant = generator
        .variable_modifications(&[r("A*", 'A', Term::Anywhere, "H")], &sequence, 1, false)
        .unwrap()
        .remove(0);
    drop(registry);
    drop(sequence);
    assert!(Arc::ptr_eq(
        variant.suffix(1).unwrap().five_prime_mod().unwrap(),
        &context
    ));
    drop(custom_end);
    assert!(weak.upgrade().is_none());
    let absent = RibonucleotideDB::default();
    let mut no_context =
        NASequence::from_records_with_registry(mods(&["A", "U"]), &absent).unwrap();
    generator
        .apply_fixed_modifications(&[r("A*", 'A', Term::Anywhere, "H")], &mut no_context)
        .unwrap();
    assert!(no_context.suffix(1).is_err());
}

#[test]
fn no_op_append_still_checks_existing_count_and_owned_payload() {
    let sequence: NASequence = "U".parse().unwrap();
    let mut output = vec![sequence.clone()];
    let before = output.clone();
    let generator = ModifiedNASequenceGenerator {
        max_outputs: 0,
        ..Default::default()
    };
    assert!(
        generator
            .apply_variable_modifications(&[], &sequence, 0, &mut output, false)
            .is_err()
    );
    assert_eq!(output, before);
    let large = Arc::new(
        Ribonucleotide::from_record(RibonucleotideRecord {
            code: "U".into(),
            origin: 'U',
            name: "x".repeat(65_536),
            ..Default::default()
        })
        .unwrap(),
    );
    let mut output = vec![NASequence::from_records(vec![large]).unwrap()];
    let before = output.clone();
    let generator = ModifiedNASequenceGenerator {
        max_output_bytes: 1000,
        ..Default::default()
    };
    assert!(
        generator
            .apply_variable_modifications(&[], &sequence, 0, &mut output, false)
            .is_err()
    );
    assert_eq!(output, before);
    let generator = ModifiedNASequenceGenerator {
        max_outputs: 0,
        max_output_bytes: 0,
        ..Default::default()
    };
    assert!(
        generator
            .variable_modifications(&[], &sequence, 0, false)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn count_work_sites_and_bytes_fail_atomically() {
    let sequence = "UUU".parse().unwrap();
    let modifications = mods(&["m3U", "s4U"]);
    for generator in [
        ModifiedNASequenceGenerator {
            max_outputs: 6,
            ..Default::default()
        },
        ModifiedNASequenceGenerator {
            max_work: 20,
            ..Default::default()
        },
        ModifiedNASequenceGenerator {
            max_sites: 2,
            ..Default::default()
        },
        ModifiedNASequenceGenerator {
            max_residues: 2,
            ..Default::default()
        },
        ModifiedNASequenceGenerator {
            max_output_bytes: 1,
            ..Default::default()
        },
    ] {
        let mut output = vec!["A".parse().unwrap()];
        let before = output.clone();
        assert!(
            generator
                .apply_variable_modifications(&modifications, &sequence, 2, &mut output, true)
                .is_err()
        );
        assert_eq!(output, before);
        let mut fixed = sequence.clone();
        // Count alone permits a single fixed output; all other tight limits fail.
        if generator.max_outputs != 6 {
            assert!(
                generator
                    .apply_fixed_modifications(&modifications, &mut fixed)
                    .is_err()
            );
            assert_eq!(fixed, sequence);
        }
    }
    let long = NASequence::parse(&"U".repeat(64)).unwrap();
    let overflow = ModifiedNASequenceGenerator {
        max_outputs: usize::MAX,
        max_output_bytes: usize::MAX,
        ..Default::default()
    };
    assert!(
        overflow
            .variable_modifications(&modifications, &long, 64, false)
            .is_err()
    );
}

#[test]
fn later_sequence_validation_failure_never_commits_original_or_variants() {
    // Valid input just below the sequence's 16 MiB rendered-text allowance.
    // Replacing its final A by a maximum-length code crosses that limit.
    let existing = r(&"x".repeat(4094), 'X', Term::Anywhere, "");
    let mut residues = vec![existing; 4095];
    residues.push(mods(&["A"]).remove(0));
    let sequence = NASequence::from_records(residues).unwrap();
    let replacement = [r(&"y".repeat(4096), 'A', Term::Anywhere, "")];
    let generator = ModifiedNASequenceGenerator::default();
    let mut output = vec![NASequence::new()];
    let before = output.clone();
    assert!(
        generator
            .apply_variable_modifications(&replacement, &sequence, 1, &mut output, true)
            .is_err()
    );
    assert_eq!(output, before);
    let mut fixed = sequence.clone();
    assert!(
        generator
            .apply_fixed_modifications(&replacement, &mut fixed)
            .is_err()
    );
    assert_eq!(fixed, sequence);
}

#[test]
fn source_no_compatible_and_zero_paths_append_exactly_one_original() {
    let generator = ModifiedNASequenceGenerator::default();
    let sequence = "A".parse().unwrap();
    for modifications in [Vec::new(), mods(&["s4U"])] {
        for maximum in [0, 1, 2] {
            let mut output = vec!["G".parse().unwrap()];
            generator
                .apply_variable_modifications(&modifications, &sequence, maximum, &mut output, true)
                .unwrap();
            assert_eq!(strings(&output), ["G", "A"]);
            assert!(
                generator
                    .variable_modifications(&modifications, &sequence, maximum, false)
                    .unwrap()
                    .is_empty()
            );
        }
    }
}
