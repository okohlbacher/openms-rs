// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::na_sequence::{MAX_NA_SEQUENCE_RESIDUES, MAX_NA_SEQUENCE_TEXT_BYTES};
use openms::chemistry::{
    ELECTRON_MASS_U, EmpiricalFormula, NAFragmentType as F, NASequence, Ribonucleotide,
    RibonucleotideDB, RibonucleotideRecord, RibonucleotideTermSpecificity as Term,
};
use std::collections::{BTreeSet, HashSet};
use std::sync::Arc;

fn formula(text: &str) -> EmpiricalFormula {
    text.parse().unwrap()
}
fn record(code: &str, composition: &str, term: Term) -> Ribonucleotide {
    Ribonucleotide::from_record(RibonucleotideRecord {
        code: code.into(),
        formula: formula(composition),
        term_specificity: term,
        mono_mass: 1234.5,
        average_mass: -42.0,
        ..Default::default()
    })
    .unwrap()
}
fn sequence(text: &str) -> NASequence {
    NASequence::parse(text).unwrap()
}

#[test]
fn source_parser_end_overwrites_and_space_rules() {
    let rna = sequence("[5'-p]A[5'-p*]C[3'-p][3'-c]");
    assert_eq!(rna.to_string(), "*ACc");
    assert_eq!(rna.len(), 2);
    assert_eq!(sequence("pA U C Gp").to_string(), "pAUCGp");
    assert_eq!(sequence("p").len(), 0);
    assert!(sequence("p").five_prime_mod().is_some());
    assert!(sequence("pp").three_prime_mod().is_some());
    assert_eq!(sequence("*c").to_string(), "*c");
    for malformed in ["[A", "A]", "[]", "A\tG", "A\nG", "[A]p]", "é"] {
        assert!(NASequence::parse(malformed).is_err(), "{malformed:?}");
    }
}

#[test]
fn source_slice_limits_zero_lengths_and_end_retention() {
    let rna = sequence("pAUCGp");
    assert_eq!(rna.subsequence(0, None).unwrap(), rna);
    assert_eq!(rna.subsequence(1, None).unwrap().to_string(), "UCGp");
    assert_eq!(rna.subsequence(0, Some(2)).unwrap().to_string(), "pAU");
    assert_eq!(rna.subsequence(2, Some(1)).unwrap().to_string(), "C");
    assert_eq!(
        rna.subsequence(2, Some(usize::MAX)).unwrap().to_string(),
        "CGp"
    );
    assert_eq!(rna.prefix(0).unwrap().to_string(), "p");
    assert_eq!(rna.suffix(0).unwrap().to_string(), "p");
    assert_eq!(rna.subsequence(0, Some(0)).unwrap().to_string(), "p");
    assert_eq!(rna.subsequence(1, Some(0)).unwrap().to_string(), "");
    assert!(rna.prefix(4).is_err());
    assert!(rna.suffix(4).is_err());
    assert!(rna.subsequence(4, Some(0)).is_err());
    assert!(NASequence::new().prefix(0).is_err());
    assert!(NASequence::new().suffix(0).is_err());
    assert!(NASequence::new().subsequence(0, None).is_err());
}

#[test]
fn all_fragment_branches_follow_elemental_algebra_and_electron_correction() {
    let rna = sequence("GG");
    let base = formula("C20H25N10O12P");
    assert_eq!(rna.formula(F::Full, 0).unwrap(), base);
    assert_eq!(rna.formula(F::Full, -2).unwrap(), formula("C20H23N10O12P"));
    for (fragment, correction) in [
        (F::AIon, "H-2O-1"),
        (F::BIon, ""),
        (F::CIon, "H-1PO2"),
        (F::DIon, "HPO3"),
        (F::WIon, "HPO3"),
        (F::XIon, "H-1PO2"),
        (F::YIon, ""),
        (F::ZIon, "H-2O-1"),
    ] {
        let expected = base
            .checked_add(&formula(correction))
            .unwrap()
            .checked_sub(&formula("H"))
            .unwrap();
        assert_eq!(rna.formula(fragment, -1).unwrap(), expected);
        assert_eq!(
            rna.mono_mass(fragment, -1).unwrap(),
            expected.mono_mass() + ELECTRON_MASS_U
        );
        assert_eq!(
            rna.average_mass(fragment, -1).unwrap(),
            expected.average_mass() + ELECTRON_MASS_U
        );
    }
    assert_eq!(
        rna.formula(F::AminusB, -1).unwrap(),
        formula("C15H17N5O10P")
    );
    let ended = sequence("pGGp");
    let with_ends = base.checked_add(&formula("H2P2O6")).unwrap();
    assert_eq!(ended.formula(F::Full, 0).unwrap(), with_ends);
    for fragment in [
        F::Internal,
        F::FivePrime,
        F::ThreePrime,
        F::Precursor,
        F::BIonMinusH2O,
        F::YIonMinusH2O,
        F::BIonMinusNH3,
        F::YIonMinusNH3,
        F::NonIdentified,
        F::Unannotated,
    ] {
        assert_eq!(ended.formula(fragment, -3).unwrap(), base);
        assert_eq!(
            ended.mono_mass(fragment, -3).unwrap(),
            base.mono_mass() + 3.0 * ELECTRON_MASS_U
        );
    }
}

#[test]
fn empty_and_placeholder_chemistry_is_source_defined() {
    let mut empty = NASequence::new();
    empty
        .set_five_prime_mod(Some(Arc::new(record("end", "C100", Term::Anywhere))))
        .unwrap();
    for fragment in [F::Full, F::AminusB, F::WIon, F::Internal] {
        assert_eq!(
            empty.formula(fragment, -2).unwrap(),
            EmpiricalFormula::default()
        );
        assert_eq!(
            empty.mono_mass(fragment, -2).unwrap(),
            2.0 * ELECTRON_MASS_U
        );
        assert_eq!(
            empty.average_mass(fragment, 3).unwrap(),
            -3.0 * ELECTRON_MASS_U
        );
    }
    let n = Arc::new(record("N", "", Term::Anywhere));
    let missing = NASequence::from_records(vec![n.clone(), n]).unwrap();
    assert_eq!(missing.formula(F::Full, 0).unwrap(), formula("H-1PO2"));
    assert!(
        missing
            .formula(F::Full, i32::MIN)
            .unwrap_err()
            .to_string()
            .contains("overflows")
    );
    empty.clear();
    assert_eq!(empty, NASequence::new());
}

#[test]
fn owned_records_ignore_declared_masses_and_keep_complete_value_identity() {
    let a = record("Q", "H2", Term::Anywhere);
    let b = record("Q", "H3", Term::Anywhere);
    let one = NASequence::from_records(vec![Arc::new(a.clone())]).unwrap();
    let same = NASequence::from_records(vec![Arc::new(a.clone())]).unwrap();
    let different = NASequence::from_records(vec![Arc::new(b.clone())]).unwrap();
    assert_eq!(one, same);
    assert_ne!(one, different);
    assert_eq!(one.to_string(), different.to_string());
    assert_eq!(
        one.mono_mass(F::Full, 0).unwrap(),
        formula("H2").mono_mass()
    );
    assert_eq!(
        one.average_mass(F::Full, 0).unwrap(),
        formula("H2").average_mass()
    );
    assert_eq!(
        BTreeSet::from([one.clone(), same.clone(), different.clone()]).len(),
        2
    );
    assert_eq!(HashSet::from([one.clone(), same, different]).len(), 2);
    let db = RibonucleotideDB::from_records(vec![b]).unwrap();
    assert!(one.checked_string_with_registry(&db).is_err());
    let owned = {
        let db = RibonucleotideDB::from_records(vec![a]).unwrap();
        let parsed = NASequence::parse_with_registry("Q", &db).unwrap();
        assert_eq!(parsed.checked_string_with_registry(&db).unwrap(), "Q");
        parsed
    };
    assert_eq!(owned, one);
    assert_eq!(owned.get_residue(0).unwrap().formula(), &formula("H2"));
}

#[test]
fn custom_sulfur_context_survives_registry_drop_and_changes_only_relevant_identity() {
    fn database(end_formula: &str) -> RibonucleotideDB {
        RibonucleotideDB::from_records(vec![
            record("C*", "C", Term::Anywhere),
            record("A", "H", Term::Anywhere),
            record("5'-p*", end_formula, Term::FivePrime),
        ])
        .unwrap()
    }
    let canonical_end = database("H2PO2S");
    let other_end = database("H2PO3S");
    let sulfur = NASequence::parse_with_registry("[C*]A", &canonical_end).unwrap();
    let other = NASequence::parse_with_registry("[C*]A", &other_end).unwrap();
    assert_ne!(sulfur, other);
    assert_eq!(
        sulfur.formula(F::Full, 0).unwrap(),
        other.formula(F::Full, 0).unwrap()
    );
    assert!(sulfur.checked_string_with_registry(&other_end).is_err());
    let plain = NASequence::parse_with_registry("A", &canonical_end).unwrap();
    assert_eq!(
        plain,
        NASequence::parse_with_registry("A", &other_end).unwrap()
    );
    assert_eq!(
        NASequence::new(),
        NASequence::parse_with_registry("", &other_end).unwrap()
    );
    drop(canonical_end);
    let suffix = sulfur.suffix(1).unwrap();
    assert_eq!(suffix.to_string(), "*A");
    assert_eq!(suffix.formula(F::WIon, 0).unwrap(), formula("H2PO2S"));
    assert_eq!(sulfur.subsequence(1, None).unwrap(), suffix);
    let missing = sulfur.clone().with_phosphorothioate_end(None).unwrap();
    assert!(missing.suffix(1).is_err());
    assert!(missing.prefix(1).is_ok());
    let last_star =
        NASequence::from_records(vec![Arc::new(record("C*", "C", Term::Anywhere))]).unwrap();
    assert_eq!(last_star.suffix(0).unwrap().to_string(), "*");
    assert_eq!(sulfur.formula(F::CIon, 0).unwrap().count("S").unwrap(), 1);
    assert_eq!(
        last_star.formula(F::CIon, 0).unwrap().count("S").unwrap(),
        1
    );
}

#[test]
fn explicit_placement_unicode_codes_and_nonroundtrippable_text_are_checked() {
    let q = record("Q", "H", Term::FivePrime);
    let unicode = record("α", "C", Term::Anywhere);
    let db = RibonucleotideDB::from_records(vec![q, unicode]).unwrap();
    let bare = NASequence::parse_with_registry("Q", &db).unwrap();
    assert_eq!(bare.len(), 1);
    assert!(bare.five_prime_mod().is_none());
    let bracketed = NASequence::parse_with_registry("[Q]", &db).unwrap();
    assert!(bracketed.is_empty());
    assert_eq!(bracketed.to_string(), "[Q]");
    assert_eq!(
        NASequence::parse_with_registry("[α]", &db)
            .unwrap()
            .to_string(),
        "[α]"
    );
    assert!(NASequence::parse_with_registry("α", &db).is_err());
    let mut placed = NASequence::new();
    placed
        .set_five_prime_mod(Some(db.get("α").unwrap()))
        .unwrap();
    assert!(placed.checked_string_with_registry(&db).is_err());
    placed.set_sequence(vec![db.get("Q").unwrap()]).unwrap();
    assert_eq!(placed.to_string(), "[α]Q");
    let before = placed.clone();
    assert!(placed.set_residue(1, db.get("Q").unwrap()).is_err());
    assert_eq!(placed, before);
    assert!(placed.get_residue(usize::MAX).is_err());
    placed.set_residue(0, db.get("α").unwrap()).unwrap();
    assert_eq!(placed.to_string(), "[α][α]");
}

#[test]
fn input_and_replacement_limits_and_formula_overflow_are_atomic() {
    assert!(NASequence::parse(&" ".repeat(MAX_NA_SEQUENCE_TEXT_BYTES + 1)).is_err());
    let r = Arc::new(Ribonucleotide::default());
    let mut seq = NASequence::from_records(vec![r.clone()]).unwrap();
    let before = seq.clone();
    assert!(
        seq.set_sequence(vec![r; MAX_NA_SEQUENCE_RESIDUES + 1])
            .is_err()
    );
    assert_eq!(seq, before);
    let max = Arc::new(record("Q", "C2147483647", Term::Anywhere));
    let plus = Arc::new(record("R", "C", Term::Anywhere));
    let overflow = NASequence::from_records(vec![max, plus]).unwrap();
    let preserved = overflow.clone();
    assert!(overflow.formula(F::Full, 0).is_err());
    assert!(overflow.mono_mass(F::Full, 0).is_err());
    assert_eq!(overflow, preserved);
}

#[test]
fn large_formula_uses_one_cumulative_work_budget() {
    let a = RibonucleotideDB::global().get("A").unwrap();
    let rna = NASequence::from_records(vec![a; 100_000]).unwrap();
    let error = rna.formula(F::Full, 0).unwrap_err();
    assert!(
        error.to_string().contains("RNA operation work limit"),
        "{error}"
    );
    assert_eq!(rna.len(), 100_000);
    assert_eq!(rna.get_residue(99_999).unwrap().code(), "A");
}
