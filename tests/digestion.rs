// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Source goldens: pinned7c029e8 ProteaseDigestion/EnzymaticDigestion/ProteaseDB tests.

use openms::chemistry::digestion::{MAX_DIGESTED_RESIDUES, MAX_DIGESTION_SEQUENCE_LENGTH};
use openms::chemistry::{
    AASequence, DigestionSpecificity as Specificity, EmpiricalFormula, ProductValidation, Protease,
    ProteaseDB, ProteaseDigestion,
};
use std::collections::BTreeSet;

fn sequence(text: &str) -> AASequence {
    text.parse().unwrap()
}
fn formula(text: &str) -> EmpiricalFormula {
    text.parse().unwrap()
}
fn checked_missed() -> ProductValidation {
    ProductValidation {
        ignore_missed_cleavages: false,
        ..Default::default()
    }
}

#[test]
fn complete_registry_metadata_and_source_synonyms() {
    let db = ProteaseDB::global();
    assert_eq!(db.enzymes().len(), 33);
    assert_eq!(Protease::ALL.len(), 33);
    let names: BTreeSet<_> = db.names().into_iter().collect();
    assert_eq!(names.len(), 33);
    for &enzyme in &Protease::ALL {
        assert_eq!(enzyme.to_string(), enzyme.metadata().name());
        assert_eq!(Protease::from_name(enzyme.name()).unwrap(), enzyme);
        assert_eq!(db.get_enzyme(enzyme.name()).unwrap(), enzyme.metadata());
        for &synonym in enzyme.metadata().synonyms() {
            assert_eq!(Protease::from_name(synonym).unwrap(), enzyme);
        }
    }
    assert_eq!(Protease::from_name("Clostripain").unwrap(), Protease::ArgC);
    assert_eq!(Protease::from_name("no_cut").unwrap(), Protease::NoCleavage);
    assert_eq!(
        Protease::from_name("nonspecific").unwrap(),
        Protease::UnspecificCleavage
    );
    assert_eq!(
        Protease::from_name("Glu-C").unwrap(),
        Protease::GlutamylEndopeptidase
    );
    assert!(!db.has_enzyme("Try"));
    assert!(!db.has_enzyme("trypsin"));
    let trypsin = db.get_enzyme("Trypsin").unwrap();
    assert_eq!(trypsin.regex(), "(?<=[KRX])(?!P)");
    assert_eq!(trypsin.n_term_gain(), formula("H"));
    assert_eq!(trypsin.c_term_gain(), formula("OH"));
    assert_eq!(trypsin.psi_id(), "MS:1001251");
    assert_eq!(trypsin.xtandem_id(), "[KR]|{P}");
    assert_eq!(trypsin.comet_id(), Some(1));
    assert_eq!(trypsin.omssa_id(), Some(0));
    assert_eq!(trypsin.msgf_id(), None);
    assert!(Protease::Iodosobenzoate.metadata().n_term_gain().is_empty());
    assert!(Protease::Iodosobenzoate.metadata().c_term_gain().is_empty());
    assert_eq!(
        db.enzymes_by_regex("(?<=[RX])(?!P)").unwrap()[0].name(),
        "Arg-C"
    );
    assert_eq!(db.enzymes_by_regex("(?<=[MX])").unwrap().len(), 2);
    assert!(!db.has_regex("(?<=[P])(?!P)"));
    assert!(db.enzymes_by_regex("(?<=[P])(?!P)").is_err());
    assert!(db.xtandem_names().contains(&"Trypsin"));
    assert!(db.xtandem_names().contains(&"no cleavage"));
    assert!(db.omssa_names().contains(&"Trypsin"));
    assert!(!db.omssa_names().contains(&"leukocyte elastase"));
    assert!(db.comet_names().contains(&"unspecific cleavage"));
    assert!(db.msgf_names().contains(&"unspecific cleavage"));
    for name in ["full", "semi", "none"] {
        assert_eq!(name.parse::<Specificity>().unwrap().to_string(), name);
    }
    assert!("no-cterm".parse::<Specificity>().is_err());
    assert!("unknown".parse::<Specificity>().is_err());
}

#[test]
fn every_compiled_rule_matches_independent_regex_context_goldens() {
    // Fixtures were generated with Python's actual regex engine from the pinned
    // XML, independently of the native predicates. All 26^2+26^3 contexts cover
    // every possible internal bond for the registry's maximum lookaround width.
    for line in include_str!("data/enzyme_cleavage_golden.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let fields: Vec<_> = line.split('\t').collect();
        let enzyme = Protease::from_name(fields[0]).unwrap();
        let expected_count: usize = fields[1].parse().unwrap();
        let expected_hash = u64::from_str_radix(fields[2], 16).unwrap();
        let mut count = 0;
        let mut hash = 14695981039346656037_u64;
        let mut add = |text: &[u8], position| {
            let text = std::str::from_utf8(text).unwrap();
            let cleaves = enzyme.is_cleavage_site(text, position).unwrap();
            count += usize::from(cleaves);
            hash = (hash ^ u64::from(cleaves)).wrapping_mul(1099511628211);
        };
        for a in b'A'..=b'Z' {
            for b in b'A'..=b'Z' {
                add(&[a, b], 1);
            }
        }
        for a in b'A'..=b'Z' {
            for b in b'A'..=b'Z' {
                for c in b'A'..=b'Z' {
                    add(&[a, b, c], 2);
                }
            }
        }
        assert_eq!(count, expected_count, "{} context count", enzyme.name());
        assert_eq!(hash, expected_hash, "{} context pattern", enzyme.name());
    }
}

#[test]
fn cleavage_context_boundaries_and_ambiguous_codes_are_exact() {
    assert_eq!(
        Protease::Trypsin.cleavage_sites("WKPAMRPAKRA").unwrap(),
        [0, 9, 10, 11]
    );
    assert_eq!(
        Protease::TrypsinP.cleavage_sites("WKPAMRPAKRA").unwrap(),
        [0, 2, 6, 9, 10, 11]
    );
    assert_eq!(
        Protease::AspN.cleavage_sites("DADBD").unwrap(),
        [0, 2, 3, 4, 5]
    );
    assert_eq!(
        Protease::AspNB.cleavage_sites("DADBD").unwrap(),
        [0, 2, 4, 5]
    );
    assert_eq!(
        Protease::FormicAcid.cleavage_sites("ADDBXA").unwrap(),
        [0, 2, 3, 4, 6]
    );
    assert_eq!(
        Protease::ProlineEndopeptidase
            .cleavage_sites("PHPAHPPRPA")
            .unwrap(),
        [0, 3, 9, 10]
    );
    assert_eq!(
        Protease::ProlineEndopeptidaseHKR
            .cleavage_sites("PHPAHPPRPA")
            .unwrap(),
        [0, 1, 3, 6, 7, 9, 10]
    );
    // X is explicitly in nearly every cleavage class, but not in iodosobenzoate.
    assert_eq!(
        Protease::Iodobenzoate.cleavage_sites("AXWA").unwrap(),
        [0, 2, 3, 4]
    );
    assert_eq!(
        Protease::Iodosobenzoate.cleavage_sites("AXWA").unwrap(),
        [0, 3, 4]
    );
    for &enzyme in &Protease::ALL {
        assert_eq!(enzyme.cleavage_sites("").unwrap(), [0]);
        assert_eq!(enzyme.cleavage_sites("A").unwrap(), [0, 1]);
        assert!(enzyme.is_cleavage_site("A", 0).unwrap());
        assert!(enzyme.is_cleavage_site("A", 1).unwrap());
        assert!(!enzyme.is_cleavage_site("A", 2).unwrap());
        assert!(!enzyme.is_cleavage_site("A", usize::MAX).unwrap());
    }
    for text in ["a", "A*", "A K", "A\nK", "Ä", "M(Oxidation)"] {
        assert!(Protease::Trypsin.cleavage_sites(text).is_err());
    }
}

#[test]
fn source_fully_specific_products_and_count_golden() {
    let mut d = ProteaseDigestion::default();
    for (mc, expected) in [(0, 4), (1, 7), (2, 9), (3, 10), (usize::MAX, 10)] {
        d.missed_cleavages = mc;
        assert_eq!(d.peptide_count(&sequence("ARCRDRE")).unwrap(), expected);
    }
    d.missed_cleavages = 1;
    assert_eq!(
        d.digest_unmodified("ARCDRE").unwrap(),
        ["AR", "CDR", "E", "ARCDR", "CDRE"]
    );
    d.min_length = 3;
    d.max_length = Some(4);
    assert_eq!(d.digest_unmodified("ARCDRE").unwrap(), ["CDR", "CDRE"]);
    assert_eq!(d.peptide_count_unmodified("ARCDRE").unwrap(), 2);
    d.set_enzyme("Asp-N").unwrap();
    d.min_length = 1;
    d.max_length = None;
    d.missed_cleavages = 0;
    assert_eq!(d.digest_unmodified("DADBD").unwrap(), ["DA", "D", "B", "D"]);
    assert!(d.digest(&sequence("")).unwrap().is_empty());
    assert_eq!(d.count_internal_cleavage_sites("").unwrap(), 0);
}

#[test]
fn source_semispecific_order_missed_cleavages_and_length_filters() {
    let mut d = ProteaseDigestion {
        enzyme: Protease::TrypsinP,
        specificity: Specificity::Semi,
        ..Default::default()
    };
    assert_eq!(
        d.digest_unmodified("AGLRMHP").unwrap(),
        [
            "AGLR", "MHP", "GLR", "MH", "LR", "M", "R", "AGL", "HP", "AG", "P", "A"
        ]
    );
    assert_eq!(d.peptide_count_unmodified("AGLRMHP").unwrap(), 12);
    let products = d.digest_ranges("AGLRMHP").unwrap();
    assert!(products.iter().all(|p| p.missed_cleavages == 0));
    d.missed_cleavages = 2;
    let products = d.digest_unmodified("AGLRMHPQGHKWYV").unwrap();
    assert!(products.contains(&"AGLRMHPQGHKW"));
    assert!(products.contains(&"GHKWYV"));
    assert!(products.contains(&"A"));
    d.missed_cleavages = 1;
    assert_eq!(
        d.digest_unmodified("AHLW").unwrap(),
        ["AHLW", "HLW", "AHL", "LW", "AH", "W", "A"]
    );
    d.min_length = 2;
    assert!(
        d.digest_unmodified("AGLRMHP")
            .unwrap()
            .iter()
            .all(|p| p.len() >= 2)
    );
    d.min_length = 3;
    d.max_length = Some(10);
    let products = d.digest_unmodified("AGLRMHPQGHKWYV").unwrap();
    assert!(!products.contains(&"AG"));
    assert!(!products.contains(&"AGLRMHPQGHKW"));
}

#[test]
fn unrestricted_substrings_match_source_and_ignore_missed_ceiling() {
    for specificity in [Specificity::Full, Specificity::Semi, Specificity::None] {
        let mut d = ProteaseDigestion {
            enzyme: Protease::UnspecificCleavage,
            specificity,
            ..Default::default()
        };
        assert_eq!(d.peptide_count_unmodified("ABCDEFGHIJ").unwrap(), 55);
        assert_eq!(
            d.digest_unmodified("ABC").unwrap(),
            ["A", "AB", "ABC", "B", "BC", "C"]
        );
        assert!(
            d.digest_ranges("ABC")
                .unwrap()
                .iter()
                .all(|p| p.missed_cleavages == p.len() - 1)
        );
        d.min_length = 8;
        d.max_length = Some(10);
        assert_eq!(
            d.digest_unmodified("ABCDEFGHIJ").unwrap(),
            [
                "ABCDEFGH",
                "ABCDEFGHI",
                "ABCDEFGHIJ",
                "BCDEFGHI",
                "BCDEFGHIJ",
                "CDEFGHIJ"
            ]
        );
        assert!(d.digest_ranges("ACD").unwrap().is_empty());
        assert!(d.digest_ranges("").unwrap().is_empty());
    }
    let d = ProteaseDigestion {
        specificity: Specificity::None,
        ..Default::default()
    };
    assert_eq!(
        d.digest_unmodified("RKR").unwrap(),
        ["R", "RK", "RKR", "K", "KR", "R"]
    );
    let p = sequence("RKR");
    assert_eq!(d.digest(&p).unwrap().len(), 6);
    // The source allows explicit missed-cleavage checking during validity even
    // though SPEC_NONE enumeration ignores that ceiling.
    assert!(!d.is_valid_product(&p, 0..3, checked_missed()).unwrap());
    assert!(
        d.is_valid_product(&p, 0..3, ProductValidation::default())
            .unwrap()
    );
}

#[test]
fn generation_matches_exhaustive_validity_oracle_and_ranges_are_unique() {
    for enzyme in [
        Protease::Trypsin,
        Protease::TrypsinP,
        Protease::AspN,
        Protease::ProlineEndopeptidase,
        Protease::NoCleavage,
        Protease::UnspecificCleavage,
    ] {
        for text in ["MACKPDERK", "HPRPDDA", "DAD", "AAA", "RK", "R", ""] {
            let cuts = enzyme.cleavage_sites(text).unwrap();
            for specificity in [Specificity::Full, Specificity::Semi, Specificity::None] {
                for missed in [0, 1, usize::MAX] {
                    let d = ProteaseDigestion {
                        enzyme,
                        specificity,
                        missed_cleavages: missed,
                        min_length: 2,
                        max_length: Some(5),
                        ..Default::default()
                    };
                    let products = d.digest_ranges(text).unwrap();
                    let actual: BTreeSet<_> = products
                        .iter()
                        .map(|p| (p.start, p.end, p.missed_cleavages))
                        .collect();
                    assert_eq!(actual.len(), products.len());
                    assert_eq!(d.peptide_count_unmodified(text).unwrap(), products.len());
                    let mut expected = BTreeSet::new();
                    for start in 0..text.len() {
                        for end in start + 1..=text.len() {
                            let mc = cuts.iter().filter(|&&p| p > start && p < end).count();
                            let n = cuts.contains(&start);
                            let c = cuts.contains(&end);
                            let unrestricted = enzyme == Protease::UnspecificCleavage
                                || specificity == Specificity::None;
                            let valid = unrestricted
                                || (mc <= missed
                                    && match specificity {
                                        Specificity::Full => n && c,
                                        Specificity::Semi => n || c,
                                        Specificity::None => true,
                                    });
                            if valid && (2..=5).contains(&(end - start)) {
                                expected.insert((start, end, mc));
                            }
                        }
                    }
                    assert_eq!(actual, expected, "{enzyme} {specificity} {missed} {text}");
                }
            }
        }
    }
}

#[test]
fn source_validity_termini_and_missed_cleavages() {
    let text = "ABCDEFGKABCRAAAKAARPBBBB";
    let mut d = ProteaseDigestion::default();
    let defaults = ProductValidation::default();
    for (range, full, semi) in [
        (0..3, false, true),
        (0..8, true, true),
        (8..12, true, true),
        (8..16, true, true),
        (0..19, false, true),
        (8..11, false, true),
        (3..9, false, false),
        (1..8, false, true),
        (0..text.len(), true, true),
    ] {
        d.specificity = Specificity::Full;
        assert_eq!(
            d.is_valid_product_unmodified(text, range.clone(), defaults)
                .unwrap(),
            full
        );
        d.specificity = Specificity::Semi;
        assert_eq!(
            d.is_valid_product_unmodified(text, range.clone(), defaults)
                .unwrap(),
            semi
        );
        d.specificity = Specificity::None;
        assert!(
            d.is_valid_product_unmodified(text, range, defaults)
                .unwrap()
        );
    }
    d.specificity = Specificity::Semi;
    assert!(
        d.is_valid_product_unmodified(text, 8..11, checked_missed())
            .unwrap()
    );
    assert!(
        !d.is_valid_product_unmodified(text, 8..13, checked_missed())
            .unwrap()
    );
    d.missed_cleavages = 1;
    assert!(
        d.is_valid_product_unmodified(text, 8..13, checked_missed())
            .unwrap()
    );
    assert!(
        !d.is_valid_product_unmodified(text, 8..18, checked_missed())
            .unwrap()
    );
    d.missed_cleavages = 2;
    assert!(
        d.is_valid_product_unmodified(text, 8..18, checked_missed())
            .unwrap()
    );
    assert_eq!(d.count_missed_cleavages(text, 8..18).unwrap(), 2);
    assert_eq!(d.count_internal_cleavage_sites(text).unwrap(), 3);
    assert!(
        !d.is_valid_product_unmodified(text, 0..text.len(), checked_missed())
            .unwrap()
    );
    d.missed_cleavages = 3;
    assert!(
        d.is_valid_product_unmodified(text, 0..text.len(), checked_missed())
            .unwrap()
    );
}

#[test]
fn methionine_and_random_asp_pro_are_validation_only_allowances() {
    let mut d = ProteaseDigestion::default();
    let nterm = ProductValidation {
        allow_nterm_protein_cleavage: true,
        ..Default::default()
    };
    for start in [1, 2] {
        assert!(
            d.is_valid_product_unmodified("MBCDEFGKABCRAAAKAA", start..8, nterm)
                .unwrap()
        );
        assert!(
            !d.is_valid_product_unmodified(
                "MBCDEFGKABCRAAAKAA",
                start..8,
                ProductValidation::default()
            )
            .unwrap()
        );
    }
    let text = "MCABCDPEFGKACDPBCRAAAKAARPBBDPBBCDP";
    let random = ProductValidation {
        allow_random_asp_pro_cleavage: true,
        ..Default::default()
    };
    d.specificity = Specificity::Semi;
    assert!(d.is_valid_product_unmodified(text, 6..9, random).unwrap());
    assert!(
        !d.is_valid_product_unmodified(text, 6..9, ProductValidation::default())
            .unwrap()
    );
    let random_checked = ProductValidation {
        ignore_missed_cleavages: false,
        ..random
    };
    assert!(
        !d.is_valid_product_unmodified(text, 6..16, random_checked)
            .unwrap()
    );
    assert!(
        d.is_valid_product_unmodified(text, 11..18, random_checked)
            .unwrap()
    );
    d.missed_cleavages = 1;
    assert!(
        d.is_valid_product_unmodified(text, 6..16, random_checked)
            .unwrap()
    );
    d.enzyme = Protease::NoCleavage;
    d.specificity = Specificity::Full;
    assert!(
        d.is_valid_product_unmodified("MADPDE", 1..6, nterm)
            .unwrap()
    );
    assert!(
        !d.is_valid_product_unmodified("MADPDE", 1..6, ProductValidation::default())
            .unwrap()
    );
    assert!(
        d.is_valid_product_unmodified("MADPDE", 0..3, random)
            .unwrap()
    );
    assert!(
        !d.is_valid_product_unmodified("MADPDE", 0..3, ProductValidation::default())
            .unwrap()
    );
    // Restoring the N-terminus includes any genuine skipped site in that prefix.
    d.enzyme = Protease::Trypsin;
    d.missed_cleavages = 0;
    let checked_nterm = ProductValidation {
        ignore_missed_cleavages: false,
        ..nterm
    };
    assert!(
        !d.is_valid_product_unmodified("MKAR", 2..4, checked_nterm)
            .unwrap()
    );
    d.missed_cleavages = 1;
    assert!(
        d.is_valid_product_unmodified("MKAR", 2..4, checked_nterm)
            .unwrap()
    );
}

#[test]
fn validity_keeps_multiresidue_enzyme_context_at_range_edges() {
    let d = ProteaseDigestion {
        enzyme: Protease::ProlineEndopeptidase,
        specificity: Specificity::None,
        ..Default::default()
    };
    // HP|A is a real cleavage even when a candidate starts at P. The source's
    // SPEC_NONE substring tokenization loses H and misses this internal site.
    assert_eq!(d.count_missed_cleavages("HPA", 1..3).unwrap(), 1);
    assert!(
        !d.is_valid_product_unmodified("HPA", 1..3, checked_missed())
            .unwrap()
    );
    assert!(
        d.is_valid_product_unmodified("HPA", 1..3, ProductValidation::default())
            .unwrap()
    );
    assert_eq!(d.count_internal_cleavage_sites("PA").unwrap(), 0);
}

#[test]
fn modified_subsequences_keep_only_original_terminal_modifications() {
    let protein = sequence(".(Acetyl)AC(Carbamidomethyl)KR.(Amidated)");
    for specificity in [Specificity::Full, Specificity::Semi, Specificity::None] {
        let d = ProteaseDigestion {
            specificity,
            missed_cleavages: 1,
            ..Default::default()
        };
        let products = d.digest(&protein).unwrap();
        assert_eq!(
            products.len(),
            d.digest_ranges(protein.as_str()).unwrap().len()
        );
        for p in products {
            assert_eq!(p.sequence, protein.subsequence(p.start..p.end).unwrap());
            assert_eq!(p.sequence.n_terminal_modification().is_some(), p.start == 0);
            assert_eq!(
                p.sequence.c_terminal_modification().is_some(),
                p.end == protein.len()
            );
            let has_c = p.start <= 1 && p.end > 1;
            assert_eq!(p.sequence.to_string().contains("Carbamidomethyl"), has_c);
        }
    }
}

#[test]
fn empty_ranges_overflows_configuration_and_resource_errors_are_safe() {
    let mut d = ProteaseDigestion::default();
    for range in [
        0..0,
        3..3,
        std::ops::Range { start: 2, end: 1 },
        0..4,
        usize::MAX..usize::MAX,
    ] {
        assert!(
            !d.is_valid_product_unmodified("AAA", range.clone(), ProductValidation::default())
                .unwrap()
        );
        assert!(d.count_missed_cleavages("AAA", range).is_err());
    }
    assert!(
        !d.is_valid_product_unmodified(
            "",
            0..0,
            ProductValidation {
                allow_nterm_protein_cleavage: true,
                ..Default::default()
            }
        )
        .unwrap()
    );
    assert!(d.set_enzyme("Unknown").is_err());
    assert_eq!(d.enzyme, Protease::Trypsin);
    for bad in [
        ProteaseDigestion {
            min_length: 0,
            ..d.clone()
        },
        ProteaseDigestion {
            min_length: 3,
            max_length: Some(2),
            ..d.clone()
        },
        ProteaseDigestion {
            max_products: 0,
            ..d.clone()
        },
        ProteaseDigestion {
            max_work: 0,
            ..d.clone()
        },
        ProteaseDigestion {
            max_work: 2,
            ..d.clone()
        },
        ProteaseDigestion {
            enzyme: Protease::UnspecificCleavage,
            max_products: 2,
            ..d.clone()
        },
    ] {
        assert!(bad.digest_ranges("AAA").is_err());
    }
    assert!(
        Protease::Trypsin
            .cleavage_sites(&"A".repeat(MAX_DIGESTION_SEQUENCE_LENGTH + 1))
            .is_err()
    );
    let huge = sequence(&"A".repeat(2000));
    let d = ProteaseDigestion {
        specificity: Specificity::None,
        min_length: 1890,
        max_length: Some(2000),
        ..d
    };
    let ranges = d.digest_ranges(huge.as_str()).unwrap();
    assert!(ranges.iter().map(|p| p.len()).sum::<usize>() > MAX_DIGESTED_RESIDUES);
    assert!(d.digest(&huge).is_err());
}
