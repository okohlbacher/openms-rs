// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::decoy_generator::{
    DecoyGenerator, MAX_DECOY_INPUT_RESIDUES, MAX_DECOY_VARIANTS,
};
use openms::chemistry::{AASequence, Protease};

fn sequence(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}

fn composition(text: &str) -> [usize; 26] {
    let mut counts = [0; 26];
    for letter in text.bytes() {
        counts[usize::from(letter - b'A')] += 1;
    }
    counts
}

#[test]
fn all_registered_enzymes_support_both_peptide_transformations() {
    let original = sequence("MACDKPERWQ");
    for enzyme in Protease::ALL {
        let mut generator = DecoyGenerator::with_seed(4711);
        let reversed = generator.reverse_peptides(&original, enzyme).unwrap();
        let shuffled = generator.shuffle_peptides(&original, enzyme, 5).unwrap();
        if enzyme != Protease::UnspecificCleavage {
            assert_eq!(
                composition(reversed.as_str()),
                composition(original.as_str())
            );
            assert_eq!(
                composition(shuffled.as_str()),
                composition(original.as_str())
            );
            assert_eq!(reversed.formula().unwrap(), original.formula().unwrap());
            assert_eq!(shuffled.formula().unwrap(), original.formula().unwrap());
        } else {
            let n = original.len();
            assert_eq!(reversed.len(), n * (n + 1) * (n + 2) / 6);
            assert_eq!(shuffled.len(), reversed.len());
        }
    }
}

#[test]
fn noncanonical_parent_letters_need_no_mass_or_formula() {
    let original = sequence("BZXJUO");
    let mut generator = DecoyGenerator::with_seed(u64::MAX);
    assert_eq!(
        generator.reverse_protein(&original).unwrap().as_str(),
        "OUJXZB"
    );
    assert!(original.formula().is_err());
    let shuffled = generator
        .shuffle_peptides(&original, Protease::NoCleavage, 100)
        .unwrap();
    assert_eq!(
        composition(shuffled.as_str()),
        composition(original.as_str())
    );
    assert!(shuffled.formula().is_err());
}

#[test]
fn reversals_and_outer_variants_do_not_change_receiver_history() {
    let target = sequence("TESTRPEPTRIDE");
    let mut generator = DecoyGenerator::with_seed(123);
    generator
        .shuffle_peptides(&target, Protease::Trypsin, 2)
        .unwrap();
    let saved = generator.clone();
    generator.reverse_protein(&target).unwrap();
    generator.reverse_peptides(&target, Protease::AspN).unwrap();
    let first = generator.shuffle(&target, Protease::TrypsinP, 3).unwrap();
    let fresh = DecoyGenerator::with_seed(987654321)
        .shuffle(&target, Protease::TrypsinP, 3)
        .unwrap();
    assert_eq!(first, fresh);
    assert_eq!(generator, saved);
}

#[test]
fn zero_attempt_cache_hit_avoids_even_an_extreme_future_attempt_limit() {
    let target = sequence("ACDEFGHIK");
    let mut generator = DecoyGenerator::with_seed(1);
    assert_eq!(
        generator
            .shuffle_peptides(&target, Protease::NoCleavage, 0)
            .unwrap(),
        target
    );
    generator.set_seed(987);
    let saved = generator.clone();
    assert_eq!(
        generator
            .shuffle_peptides(&target, Protease::NoCleavage, usize::MAX)
            .unwrap(),
        target
    );
    assert_eq!(generator, saved);
}

#[test]
fn every_modification_kind_errors_without_changing_state_even_for_zero_variants() {
    let mut generator = DecoyGenerator::with_seed(4711);
    for text in [
        "M(Oxidation)AC",
        "(Acetyl)AC",
        "AC(Amidated)",
        "AX[123.4567]",
        ".[+1.23456789]AC",
    ] {
        let protein = sequence(text);
        assert!(protein.is_modified(), "{text}");
        let saved = generator.clone();
        assert!(generator.reverse_protein(&protein).is_err());
        assert!(
            generator
                .reverse_peptides(&protein, Protease::Trypsin)
                .is_err()
        );
        assert!(
            generator
                .shuffle_peptides(&protein, Protease::Trypsin, 0)
                .is_err()
        );
        assert!(generator.shuffle(&protein, Protease::Trypsin, 0).is_err());
        assert_eq!(generator, saved);
    }
}

#[test]
fn unspecific_nested_products_and_aggregate_variant_lengths_are_bounded() {
    let generator = DecoyGenerator::with_seed(4711);
    // ABC -> A|B|C|AB|BC|ABC; only the length-three product is redigested
    // inside outer shuffle, making 1+1+1+2+2+10 = 17 residues per variant.
    let results = generator
        .shuffle(&sequence("ABC"), Protease::UnspecificCleavage, 2)
        .unwrap();
    assert_eq!(
        results.iter().map(AASequence::len).collect::<Vec<_>>(),
        [17, 17]
    );
    let crowded = sequence(&"A".repeat(446));
    // 446*447/2 products fit the product cap, but their overlapping residue
    // total exceeds the output cap before any reversal/candidate allocation.
    assert!(
        generator
            .reverse_peptides(&crowded, Protease::UnspecificCleavage)
            .unwrap_err()
            .to_string()
            .contains("output residue limit")
    );
    assert!(
        generator
            .reverse_peptides(&sequence(&"A".repeat(447)), Protease::UnspecificCleavage)
            .unwrap_err()
            .to_string()
            .contains("peptide product limit")
    );
    assert!(
        generator
            .shuffle(
                &sequence(&"A".repeat(1001)),
                Protease::NoCleavage,
                MAX_DECOY_VARIANTS
            )
            .unwrap_err()
            .to_string()
            .contains("output residue limit")
    );
}

#[test]
fn input_and_variant_limits_are_checked_and_empty_variants_remain_defined() {
    let generator = DecoyGenerator::with_seed(1);
    assert_eq!(
        generator.reverse_protein(&AASequence::default()).unwrap(),
        AASequence::default()
    );
    assert!(
        generator
            .reverse_peptides(&AASequence::default(), Protease::NoCleavage)
            .is_err()
    );
    assert_eq!(
        generator
            .shuffle(&AASequence::default(), Protease::NoCleavage, 3)
            .unwrap(),
        vec![AASequence::default(); 3]
    );
    assert!(
        generator
            .shuffle(
                &AASequence::default(),
                Protease::NoCleavage,
                MAX_DECOY_VARIANTS + 1
            )
            .unwrap_err()
            .to_string()
            .contains("variant limit")
    );
    // X avoids expensive known-residue formula construction in this input-cap test.
    let oversized = sequence(&"X".repeat(MAX_DECOY_INPUT_RESIDUES + 1));
    assert!(
        generator
            .reverse_protein(&oversized)
            .unwrap_err()
            .to_string()
            .contains("input residue limit")
    );
}
