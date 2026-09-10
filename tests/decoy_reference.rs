// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Literal source strings and separately identified hand-derived state/ordering
// checks. See data/decoy_provenance.json. No C++ execution produced these data.

use openms::chemistry::decoy_generator::DecoyGenerator;
use openms::chemistry::{AASequence, Protease};

fn sequence(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}

fn enzyme(name: &str) -> Protease {
    Protease::from_name(name).unwrap()
}

fn text(values: &[AASequence]) -> Vec<&str> {
    values.iter().map(AASequence::as_str).collect()
}

#[test]
fn all_thirteen_literal_source_strings_preserve_the_shared_call_history() {
    let mut generator = DecoyGenerator::with_seed(4711);
    let mut stateful_order = 0;
    let mut rows = 0;
    for line in include_str!("data/decoy_source_goldens.tsv")
        .lines()
        .skip(1)
    {
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields.len(), 8);
        assert_eq!(fields[3], "4711");
        let input = sequence(fields[1]);
        let actual = match fields[0] {
            "reverse_protein" => generator.reverse_protein(&input).unwrap(),
            "reverse_peptides" => generator
                .reverse_peptides(&input, enzyme(fields[2]))
                .unwrap(),
            "shuffle_peptides" => {
                stateful_order += 1;
                assert_eq!(fields[4].parse::<usize>().unwrap(), stateful_order);
                generator
                    .shuffle_peptides(&input, enzyme(fields[2]), 100)
                    .unwrap()
            }
            "shuffle" => {
                let variant = fields[5].parse::<usize>().unwrap();
                let before = generator.clone();
                let variants = generator
                    .shuffle(&input, enzyme(fields[2]), variant + 1)
                    .unwrap();
                assert_eq!(
                    generator, before,
                    "outer shuffle must ignore receiver history"
                );
                variants[variant].clone()
            }
            operation => panic!("unknown fixture operation {operation}"),
        };
        assert_eq!(actual.as_str(), fields[6], "source lines {}", fields[7]);
        rows += 1;
    }
    assert_eq!(rows, 13);
    assert_eq!(stateful_order, 3);
}

#[test]
fn zero_attempts_cache_the_original_and_reseeding_keeps_cached_choices() {
    let input = sequence("TESTPEPTIDE");
    let trypsin = enzyme("Trypsin");
    let mut generator = DecoyGenerator::with_seed(4711);
    assert_eq!(
        generator.shuffle_peptides(&input, trypsin, 0).unwrap(),
        input
    );
    let cached = generator.clone();
    assert_eq!(
        generator.shuffle_peptides(&input, trypsin, 100).unwrap(),
        input
    );
    assert_eq!(generator, cached, "a cache hit must consume no draws");
    generator.set_seed(u64::MAX);
    let reseeded = generator.clone();
    assert_eq!(
        generator.shuffle_peptides(&input, trypsin, 100).unwrap(),
        input
    );
    assert_eq!(generator, reseeded);

    let mut fresh = DecoyGenerator::with_seed(4711);
    assert_eq!(
        fresh
            .shuffle_peptides(&input, trypsin, 100)
            .unwrap()
            .as_str(),
        "DIESETEPTTP"
    );
    fresh.set_seed(0);
    let cached_decoy = fresh.clone();
    assert_eq!(
        fresh.shuffle_peptides(&input, trypsin, 0).unwrap().as_str(),
        "DIESETEPTTP"
    );
    assert_eq!(fresh, cached_decoy);
}

#[test]
fn a_cached_final_product_is_reused_across_enzyme_and_nonfinal_context() {
    let mut generator = DecoyGenerator::with_seed(4711);
    // The source outer Trypsin golden ETEPTSRRTP|DIE uses a fresh seed 4711
    // for TESTRPEPTR. Both enzymes leave that isolated product undivided.
    assert_eq!(
        generator
            .shuffle_peptides(&sequence("TESTRPEPTR"), enzyme("no cleavage"), 100)
            .unwrap()
            .as_str(),
        "ETEPTSRRTP"
    );
    // That cached final-product choice now moves the nonfinal cleavage residue.
    // The uncached final IDE has zero attempts and therefore remains unchanged.
    let input = sequence("TESTRPEPTRIDE");
    assert_eq!(
        generator
            .shuffle_peptides(&input, enzyme("Trypsin"), 0)
            .unwrap()
            .as_str(),
        "ETEPTSRRTPIDE"
    );
    let before = generator.clone();
    assert_eq!(
        generator
            .shuffle_peptides(&input, enzyme("Trypsin"), 100)
            .unwrap()
            .as_str(),
        "ETEPTSRRTPIDE"
    );
    assert_eq!(generator, before);
}

#[test]
fn short_final_ranges_consume_draws_even_when_no_better_decoy_exists() {
    let input = sequence("AG");
    let no_cleavage = enzyme("no cleavage");
    let mut attempted = DecoyGenerator::with_seed(4711);
    let mut skipped = attempted.clone();
    assert_eq!(
        attempted
            .shuffle_peptides(&input, no_cleavage, 100)
            .unwrap(),
        input
    );
    assert_eq!(
        skipped.shuffle_peptides(&input, no_cleavage, 0).unwrap(),
        input
    );
    // Every two-residue permutation has forward or reverse identity one.
    // Cache contents agree, but 100 attempts advance the random stream.
    assert_ne!(attempted, skipped);
    attempted.set_seed(4711);
    assert_eq!(attempted, skipped);

    let input = sequence("AKR");
    let mut attempted = DecoyGenerator::with_seed(4711);
    let mut skipped = attempted.clone();
    // Nonfinal AK shuffles only A; final R also has a one-element range.
    assert_eq!(
        attempted
            .shuffle_peptides(&input, enzyme("Trypsin"), 100)
            .unwrap(),
        input
    );
    assert_eq!(
        skipped
            .shuffle_peptides(&input, enzyme("Trypsin"), 0)
            .unwrap(),
        input
    );
    assert_eq!(
        attempted, skipped,
        "one-element shuffle ranges consume no draws"
    );
}

#[test]
fn positional_anchors_and_unspecific_length_then_start_order_match_source_loops() {
    let mut generator = DecoyGenerator::with_seed(4711);
    let before = generator.clone();
    assert_eq!(
        generator
            .reverse_peptides(&sequence("AK"), enzyme("Trypsin"))
            .unwrap()
            .as_str(),
        "KA"
    );
    // Asp-N cuts AC | DPEF | DGK. The source still anchors each nonfinal
    // product's final letter: AC | EPDF | KGD.
    assert_eq!(
        generator
            .reverse_peptides(&sequence("ACDPEFDGK"), enzyme("Asp-N"))
            .unwrap()
            .as_str(),
        "ACEPDFKGD"
    );
    // Source AASequence digestion emits A, B, C, AB, BC, ABC in this order.
    assert_eq!(
        generator
            .reverse_peptides(&sequence("ABC"), enzyme("unspecific cleavage"))
            .unwrap()
            .as_str(),
        "ABCABBCCBA"
    );
    assert_eq!(generator, before);
    assert_eq!(
        generator
            .shuffle_peptides(&sequence("ABC"), enzyme("unspecific cleavage"), 0)
            .unwrap()
            .as_str(),
        "ABCABBCABC"
    );
}

#[test]
fn outer_variants_ignore_receiver_state_keep_short_duplicates_and_redigest_products() {
    let mut generator = DecoyGenerator::with_seed(u64::MAX);
    let input = sequence("TESTPEPTIDE");
    generator
        .shuffle_peptides(&input, enzyme("Trypsin"), 0)
        .unwrap();
    let before = generator.clone();
    let variants = generator.shuffle(&input, enzyme("Trypsin"), 2).unwrap();
    assert_eq!(text(&variants), ["DIESETEPTTP", "PTEPIDEETTS"]);
    assert_eq!(generator, before);
    assert_eq!(
        text(
            &generator
                .shuffle(&sequence("AKR"), enzyme("Trypsin"), 3)
                .unwrap()
        ),
        ["AKR", "AKR", "AKR"]
    );
    assert_eq!(
        text(
            &generator
                .shuffle(&sequence("AG"), enzyme("no cleavage"), 2)
                .unwrap()
        ),
        ["AG", "AG"]
    );

    let variants = generator
        .shuffle(&sequence("ABC"), enzyme("unspecific cleavage"), 1)
        .unwrap();
    let result = variants[0].as_str();
    // The outer A,B,C,AB,BC contribute seven unchanged letters. The long ABC
    // is redigested, contributing the same seven-letter prefix plus a shuffle
    // of ABC. This checks source expansion without inventing an RNG golden.
    assert_eq!(result.len(), 17);
    assert!(result.starts_with("ABCABBCABCABBC"));
    let mut tail = result.as_bytes()[14..].to_vec();
    tail.sort_unstable();
    assert_eq!(tail, b"ABC");
    assert_eq!(generator, before);
}

#[test]
fn empty_entrypoints_differ_and_all_modified_forms_are_rejected_atomically() {
    let mut generator = DecoyGenerator::with_seed(4711);
    let empty = AASequence::default();
    let trypsin = enzyme("Trypsin");
    let original = generator.clone();
    assert_eq!(generator.reverse_protein(&empty).unwrap(), empty);
    assert!(generator.reverse_peptides(&empty, trypsin).is_err());
    for attempts in [0, 100] {
        assert!(
            generator
                .shuffle_peptides(&empty, trypsin, attempts)
                .is_err()
        );
        assert_eq!(generator, original);
    }
    assert_eq!(
        text(&generator.shuffle(&empty, trypsin, 2).unwrap()),
        ["", ""]
    );
    assert!(generator.shuffle(&empty, trypsin, 0).unwrap().is_empty());
    assert!(
        generator
            .shuffle(&sequence("AMK"), trypsin, 0)
            .unwrap()
            .is_empty()
    );

    let mut residue_tag = sequence("AMK");
    residue_tag.set_mass_tag(1, "+0.123456789").unwrap();
    let mut n_term = sequence("AMK");
    n_term.set_n_terminal_mass_tag("+0.123456789").unwrap();
    let mut c_term = sequence("AMK");
    c_term.set_c_terminal_mass_tag("+0.123456789").unwrap();
    for modified in [sequence("AM(Oxidation)K"), residue_tag, n_term, c_term] {
        assert!(modified.is_modified());
        assert!(generator.reverse_protein(&modified).is_err());
        assert!(generator.reverse_peptides(&modified, trypsin).is_err());
        assert!(generator.shuffle_peptides(&modified, trypsin, 0).is_err());
        assert!(generator.shuffle(&modified, trypsin, 0).is_err());
        assert_eq!(generator, original);
    }
}
