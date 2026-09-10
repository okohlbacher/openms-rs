// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::{AAIndex, AAIndexScale, AASequence};

fn sequence(input: &str) -> AASequence {
    AASequence::parse(input).unwrap()
}

#[test]
fn all_accessions_are_stable_and_have_checked_numeric_lookups() {
    let names = [
        "KHAG800101",
        "VASM830103",
        "NADH010106",
        "NADH010107",
        "WILM950102",
        "ROBB760107",
        "OOBM850104",
        "FAUJ880111",
        "FINA770101",
        "ARGP820102",
    ];
    // Source AAIndex.h A entries; the independent reference fixture covers all
    // 200 values without borrowing production table storage.
    let alanine = [49.1, 0.159, 5.0, -2.0, 2.62, 0.0, -2.49, 0.0, 1.08, 1.18];
    for ((scale, name), expected) in AAIndexScale::ALL.into_iter().zip(names).zip(alanine) {
        assert_eq!(scale.accession(), name);
        assert_eq!(scale.value('A').unwrap(), expected);
        for residue in ['B', 'J', 'O', 'U', 'X', 'Z', 'a', '*', 'é', '\0'] {
            assert!(scale.value(residue).is_err(), "{name}, {residue:?}");
        }
    }
}

#[test]
fn source_indicator_memberships_do_not_use_textbook_or_numeric_index_categories() {
    for c in 'A'..='Z' {
        assert_eq!(
            AAIndex::aliphatic(c),
            if "AGFIMLPV".contains(c) { 1.0 } else { 0.0 }
        );
        assert_eq!(AAIndex::acidic(c), if "DE".contains(c) { 1.0 } else { 0.0 });
        assert_eq!(
            AAIndex::basic(c),
            if "KRHW".contains(c) { 1.0 } else { 0.0 }
        );
        assert_eq!(
            AAIndex::polar(c),
            if "STYHCNQW".contains(c) { 1.0 } else { 0.0 }
        );
    }
    assert_eq!(AAIndex::basic('W'), 1.0);
    assert_eq!(AAIndexScale::Fauj880111.value('W').unwrap(), 0.0);
    for c in ['a', 'k', '?', '\0', 'é', '🧬'] {
        assert_eq!(
            [
                AAIndex::aliphatic(c),
                AAIndex::acidic(c),
                AAIndex::basic(c),
                AAIndex::polar(c)
            ],
            [0.0; 4]
        );
    }
}

#[test]
#[allow(clippy::approx_constant)] // Source arginine energy delta 6.28, not tau.
fn source_gb500_goldens_and_direct_split_order_are_preserved() {
    // Literal AAIndex_test.cpp values and its 0.01 absolute tolerance, rev7c029e8.
    for (peptide, expected) in [
        ("ALEGDEK", 1337.53),
        ("GTVVTGR", 1442.70),
        ("EHVLLAR", 1442.70),
    ] {
        assert!(
            (AAIndex::calculate_gb(&sequence(peptide), 500.0).unwrap() - expected).abs() <= 0.01
        );
    }
    // Source-derived arithmetic for ARR, not captured C++ output. Split pairing
    // and ordinary direct exp/sum order differ by one ULP from unconditional LSE.
    let rt: f64 = ((6.022_136_7e23 * 1.380_657e-23) / 1000.0) * 500.0;
    let mut association = 0.0;
    association += (916.84 / rt).exp();
    association += ((881.82 + 6.28) / rt).exp() + (1000.0 / rt).exp();
    association += ((882.98 + 6.28) / rt).exp() + (1000.0 / rt).exp();
    association += ((882.98 - 95.82) / rt).exp();
    let direct = rt * association.ln() / 2.0_f64.ln();
    assert_eq!(
        AAIndex::calculate_gb(&sequence("ARR"), 500.0)
            .unwrap()
            .to_bits(),
        direct.to_bits()
    );
    let ak = AAIndex::calculate_gb(&sequence("AK"), 500.0).unwrap();
    let ka = AAIndex::calculate_gb(&sequence("KA"), 500.0).unwrap();
    assert!((ak - 1337.5334766479912).abs() < 3e-12);
    assert!((ka - 1321.6972307017786).abs() < 3e-12);
    assert!(ak > ka + 15.0); // The first sidechain is omitted.
}

#[test]
fn stable_gb_at_100_kelvin_retains_ties_and_matches_decimal_derivations() {
    // Independently evaluated at 90 Decimal digits using source table/constant
    // values; source's 100 K tests only assert inequality, not these numbers.
    for (peptide, expected) in [
        ("ALEGDEK", 1337.0032102826751),
        ("GTVVTGR", 1442.6950408889634),
        ("EHVLLAR", 1442.6950408889634),
        ("ARR", 1443.5264914079446),
    ] {
        let value = AAIndex::calculate_gb(&sequence(peptide), 100.0).unwrap();
        assert!((value - expected).abs() < 5e-12, "{peptide}: {value}");
    }
    let single = AAIndex::calculate_gb(&sequence("AR"), 100.0).unwrap();
    let tied = AAIndex::calculate_gb(&sequence("ARR"), 100.0).unwrap();
    assert!((tied - single - 0.8314505189811898).abs() < 5e-12);
}

#[test]
fn positive_extreme_temperatures_and_empty_identity_remain_finite() {
    let empty = sequence("");
    let expected = (916.84_f64 - 95.82) / 2.0_f64.ln();
    for temperature in [
        f64::from_bits(1),
        f64::MIN_POSITIVE,
        1e-200,
        100.0,
        500.0,
        1e20,
        1e100,
        f64::MAX,
    ] {
        assert_eq!(
            AAIndex::calculate_gb(&empty, temperature).unwrap(),
            expected
        );
        let value = AAIndex::calculate_gb(&sequence("AK"), temperature).unwrap();
        assert!(value.is_finite() && value > 0.0);
    }
    let limiting = 926.74 / 2.0_f64.ln();
    for temperature in [
        f64::from_bits(1),
        f64::from_bits(61),
        f64::MIN_POSITIVE,
        1e-200,
    ] {
        assert_eq!(
            AAIndex::calculate_gb(&sequence("AK"), temperature).unwrap(),
            limiting
        );
    }
}

#[test]
fn all_annotations_are_ignored_without_mutation_or_formula_requirements() {
    let plain = sequence("ACMK");
    let annotated = sequence("(Acetyl)AC(Carbamidomethyl)M[+12.3456789]K.[+0.123456789]");
    assert!(annotated.formula().is_err());
    let before = annotated.clone();
    for temperature in [100.0, 500.0, 1e100] {
        assert_eq!(
            AAIndex::calculate_gb(&plain, temperature).unwrap(),
            AAIndex::calculate_gb(&annotated, temperature).unwrap()
        );
    }
    assert_eq!(annotated, before);
}

#[test]
fn invalid_temperatures_and_late_unsupported_residues_are_checked() {
    for temperature in [0.0, -0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for peptide in ["", "A"] {
            assert!(AAIndex::calculate_gb(&sequence(peptide), temperature).is_err());
        }
    }
    for residue in ['B', 'J', 'O', 'U', 'X', 'Z'] {
        let peptide = sequence(&format!("AK{residue}"));
        let before = peptide.clone();
        for temperature in [f64::from_bits(1), 100.0, 500.0] {
            assert!(AAIndex::calculate_gb(&peptide, temperature).is_err());
        }
        assert_eq!(peptide, before);
    }
}
