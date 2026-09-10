// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Derived from pinned AASequence_test.cpp 7c029e8; expected mass values below
// are upstream test constants, not output captured from a C++ build.
use openms::Error;
use openms::chemistry::{
    AASequence, EmpiricalFormula, ModificationsDB, PROTON_MASS_U, SequenceModification,
    TermSpecificity,
};
fn seq(s: &str) -> AASequence {
    AASequence::parse(s).unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:.14} != {b:.14}");
}
fn unavailable<T: std::fmt::Debug>(value: openms::Result<T>) {
    assert!(matches!(value, Err(Error::Unsupported(_))), "{value:?}");
}

#[test]
fn anonymous_spellings_are_shared_by_independently_owned_fragment_slices() {
    let text = format!("+12.3456789{}", "0".repeat(900));
    let peptide = seq(&format!("AM[{text}]KA"));
    let prefix = peptide.prefix(3).unwrap();
    let suffix = peptide.suffix(3).unwrap();
    let tag = peptide
        .residue_modification(1)
        .unwrap()
        .unwrap()
        .mass_tag()
        .unwrap();
    for (fragment, index) in [(&prefix, 1), (&suffix, 0)] {
        let other = fragment
            .residue_modification(index)
            .unwrap()
            .unwrap()
            .mass_tag()
            .unwrap();
        // Fragment creation must not recopy long immutable decimal spellings.
        assert!(std::ptr::eq(tag.input(), other.input()));
        assert!(std::ptr::eq(tag.full_id(), other.full_id()));
        assert_eq!(tag, other);
    }
    drop(peptide);
    let original = prefix
        .residue_modification(1)
        .unwrap()
        .unwrap()
        .mass_tag()
        .unwrap()
        .clone();
    let expected_mass = prefix.mono_mass().unwrap();
    assert_eq!(original.input(), text);
    let mut changed = prefix.clone();
    changed.set_mass_tag(1, "+1.23456789").unwrap();
    assert_eq!(prefix.mono_mass().unwrap(), expected_mass);
    assert_eq!(
        prefix
            .residue_modification(1)
            .unwrap()
            .unwrap()
            .mass_tag()
            .unwrap(),
        &original
    );
    assert_eq!(
        suffix
            .residue_modification(0)
            .unwrap()
            .unwrap()
            .mass_tag()
            .unwrap(),
        &original
    );
    assert_ne!(
        changed
            .residue_modification(1)
            .unwrap()
            .unwrap()
            .mass_tag()
            .unwrap(),
        &original
    );
}

#[test]
fn ambiguous_residues_are_represented_without_placeholder_chemistry() {
    for residue in ["B", "Z", "X"] {
        for text in [
            residue.to_owned(),
            format!("PEP{residue}TIDE"),
            format!("(Acetyl){residue}(Amidated)"),
        ] {
            let sequence = seq(&text);
            unavailable(sequence.formula());
            unavailable(sequence.mono_mass());
            unavailable(sequence.average_mass());
            unavailable(sequence.mz(2));
            assert_eq!(seq(&sequence.to_string()), sequence);
        }
        let mut sequence = seq(residue);
        let before = sequence.clone();
        assert!(sequence.set_mass_tag(0, "+0").is_err());
        assert_eq!(sequence, before);
    }
    assert_eq!(seq("BZXXZB").as_str(), "BZXXZB");
    assert!(!seq("BZXXZB").is_modified());
    assert!(seq("J").formula().is_ok());
    assert_eq!(seq("J").mono_mass().unwrap(), seq("I").mono_mass().unwrap());
}

#[test]
fn source_integer_absolute_mass_lookup_cases() {
    let cases = [
        ("PEPTIDEK[136]", 7, "Label:13C(6)15N(2)"),
        ("PEPS[167]TIDEK", 3, "Phospho"),
        ("PEPC[160]TIDEK", 3, "Carbamidomethyl"),
        ("PEPM[147]TIDEK", 3, "Oxidation"),
        ("PEPT[181]TIDEK", 3, "Phospho"),
        ("PEPY[243]TIDEK", 3, "Phospho"),
        ("PEPR[166]TIDEK", 3, "Label:13C(6)15N(4)"),
    ];
    for (text, index, name) in cases {
        let sequence = seq(text);
        let modification = sequence.residue_modification(index).unwrap().unwrap();
        assert_eq!(modification.name(), name);
        assert!(modification.known().is_some());
        assert!(modification.record_id().is_some());
        assert!(sequence.formula().is_ok());
        assert!(sequence.average_mass().is_ok());
        assert_eq!(seq(&sequence.to_unimod_string().unwrap()), sequence);
    }
}

#[test]
fn source_explicit_numeric_terminal_lookup() {
    for (text, name) in [
        ("n[+40]CPEPTIDEK", "Pyro-carbamidomethyl"),
        ("n[-17]QPEPTIDEK", "Gln->pyro-Glu"),
        ("[+42].MVLVQDLLHPTAASEAR", "Acetyl"),
        ("[+304.207].ETC[+57.0215]RQLGLGTNIYNAER", "TMTpro"),
        ("n[+34]TGSESSQTGTSTTSSR", "Dimethyl:2H(4)13C(2)"),
        ("n[35]TGSESSQTGTSTTSSR", "Dimethyl:2H(4)13C(2)"),
    ] {
        let sequence = seq(text);
        assert_eq!(sequence.n_terminal_modification().unwrap().name(), name);
        assert!(
            sequence
                .n_terminal_modification()
                .unwrap()
                .known()
                .is_some()
        );
        assert_eq!(seq(&sequence.to_string()), sequence);
    }
    let sequence = seq("PEPTIDEc[-1]");
    assert_eq!(
        sequence.c_terminal_modification().unwrap().name(),
        "Amidated"
    );
}

#[test]
fn source_decimal_precision_separates_known_and_anonymous_masses() {
    let coarse = seq("PEPM[147.0354]TIDEK");
    assert_eq!(
        coarse.residue_modification(3).unwrap().unwrap().name(),
        "Oxidation"
    );
    let precise = seq("PEPM[147.035405]TIDEK");
    let annotation = precise.residue_modification(3).unwrap().unwrap();
    assert!(annotation.mass_tag().is_some());
    assert_eq!(annotation.mass_tag().unwrap().input(), "147.035405");
    assert_eq!(annotation.full_id(), "M[147.035405]");
    assert_eq!(annotation.origin(), Some('M'));
    assert_eq!(annotation.term_specificity(), TermSpecificity::Anywhere);
    assert_eq!(annotation.record_id(), None);
    near(
        precise.mono_mass().unwrap(),
        seq("PEPTIDEK").mono_mass().unwrap() + 147.035405,
        1e-10,
    );
    unavailable(precise.formula());
    unavailable(precise.average_mass());
    unavailable(annotation.diff_formula());
    unavailable(annotation.diff_average_mass());
}

#[test]
fn source_fractional_delta_regressions_preserve_exact_mass_and_spelling() {
    let unmodified = seq("PEPTIDE").mono_mass().unwrap();
    for (text, delta) in [
        ("-0.5", -0.5),
        ("+0.001", 0.001),
        ("-0.001", -0.001),
        ("+0.000000001", 1e-9),
        ("+0.00335", 0.00335),
    ] {
        let input = format!("PEPT[{text}]IDE");
        let sequence = seq(&input);
        let annotation = sequence.residue_modification(3).unwrap().unwrap();
        let tag = annotation.mass_tag().unwrap();
        assert!(tag.is_delta());
        assert_eq!(tag.mass(), delta);
        assert_eq!(tag.delta_mono_mass(), Some(delta));
        assert_eq!(annotation.diff_mono_mass().unwrap(), delta);
        near(sequence.mono_mass().unwrap() - unmodified, delta, 2e-13);
        assert_eq!(sequence.to_string(), input);
        assert_eq!(sequence.to_unimod_string().unwrap(), input);
        assert_eq!(seq(&sequence.to_string()), sequence);
    }
}

#[test]
fn source_annotated_mass_goldens() {
    for text in [
        "PEPTC[+57.02]IDE",
        "PEPTC(Carbamidomethyl)IDE",
        "PEPTC[160.030654]IDE",
        "PEPTX[160.030654]IDE",
    ] {
        near(seq(text).mono_mass().unwrap(), 959.39066, 1e-4);
    }
    near(
        seq("PEPTM[-30]IDE").mono_mass().unwrap(),
        930.4004 - 29.992806,
        1e-4,
    );
    near(seq("PEPTM[-30.4004]IDE").mono_mass().unwrap(), 900.0, 1e-4);
    near(seq("PEPTIDE").mono_mass().unwrap(), 799.36001, 1e-4);
    near(seq("IDE").mono_mass().unwrap(), 375.1641677975, 1e-10);
    near(
        seq("DFPANGERX[113.0840643509]").mono_mass().unwrap(),
        904.4038997864 + 113.0840643509,
        1e-9,
    );
}

#[test]
fn absolute_tags_on_ambiguous_residues_have_only_mono_mass() {
    for residue in ['B', 'Z', 'X'] {
        let sequence = seq(&format!("PEPT{residue}[999.000]IDE"));
        near(
            sequence.mono_mass().unwrap(),
            seq("PEPTIDE").mono_mass().unwrap() + 999.,
            1e-10,
        );
        let annotation = sequence.residue_modification(4).unwrap().unwrap();
        let tag = annotation.mass_tag().unwrap();
        assert_eq!(tag.mass(), 999.0);
        assert!(!tag.is_delta());
        assert_eq!(tag.delta_mono_mass(), None);
        assert_eq!(tag.residue_mono_mass(), Some(999.));
        unavailable(annotation.diff_mono_mass());
        unavailable(sequence.formula());
        unavailable(sequence.average_mass());
        assert_eq!(seq(&sequence.to_string()), sequence);
        near(
            sequence.mz(2).unwrap(),
            sequence.mono_mass().unwrap() / 2. + PROTON_MASS_U,
            1e-12,
        );
    }
    // X[0] is an explicit known zero INTERNAL mass, unlike an unresolved bare X.
    near(
        seq("X[0]").mono_mass().unwrap(),
        EmpiricalFormula::parse("H2O").unwrap().mono_mass(),
        1e-12,
    );
}

#[test]
fn source_anonymous_terminal_absolute_h_and_oh_conventions() {
    let base = seq("IDE").mono_mass().unwrap();
    for (text, delta, full_id, term) in [
        (
            ".[1601.2384790319]IDE",
            1600.230654,
            ".n[1601.2384790319]",
            TermSpecificity::NTerm,
        ),
        (
            "n[+1600.230654]IDE",
            1600.230654,
            ".n[+1600.230654]",
            TermSpecificity::NTerm,
        ),
        (
            "IDE.[1617.2333940319]",
            1600.230654,
            ".c[1617.2333940319]",
            TermSpecificity::CTerm,
        ),
        (
            "IDEc[+1600.230654]",
            1600.230654,
            ".c[+1600.230654]",
            TermSpecificity::CTerm,
        ),
    ] {
        let sequence = seq(text);
        near(sequence.mono_mass().unwrap(), base + delta, 1e-9);
        let annotation = if term == TermSpecificity::NTerm {
            sequence.n_terminal_modification()
        } else {
            sequence.c_terminal_modification()
        }
        .unwrap();
        assert_eq!(annotation.full_id(), full_id);
        assert_eq!(annotation.term_specificity(), term);
        assert_eq!(annotation.origin(), None);
        assert_eq!(annotation.mass_tag().unwrap().residue_mono_mass(), None);
        near(annotation.diff_mono_mass().unwrap(), delta, 1e-9);
        unavailable(sequence.formula());
        unavailable(sequence.average_mass());
        assert_eq!(seq(&sequence.to_string()), sequence);
    }
}

#[test]
fn slicing_recovers_known_chemistry_and_preserves_owned_annotation_identity() {
    let source = seq("(Acetyl)ABX[999.000]CDZ");
    unavailable(source.mono_mass());
    unavailable(source.formula());
    let sliced = source.subsequence(2..5).unwrap();
    assert_eq!(sliced.to_string(), "X[999.000]CD");
    assert!(sliced.mono_mass().is_ok());
    unavailable(sliced.formula());
    let known = sliced.suffix(2).unwrap();
    assert_eq!(known, seq("CD"));
    assert!(known.formula().is_ok());
    assert!(known.average_mass().is_ok());
    let first = source.prefix(1).unwrap();
    assert_eq!(first, seq("(Acetyl)A"));
    let empty = source.subsequence(1..1).unwrap();
    assert_eq!(empty, AASequence::default());
    assert_eq!(empty.mono_mass().unwrap(), 0.0);
    assert!(source.subsequence(0..99).is_err());
    assert!(
        source
            .subsequence(std::ops::Range { start: 4, end: 2 })
            .is_err()
    );
}

#[test]
fn mass_only_b_y_fragments_follow_the_retained_annotation() {
    let base = seq("ACMK");
    let modified = seq("AC[+12.3456789]MK");
    let original = base.fragment_ions(2).unwrap();
    let ions = modified.fragment_ions(2).unwrap();
    for (before, after) in original.iter().zip(&ions) {
        let retained = match after.series {
            openms::chemistry::IonSeries::B => after.ordinal >= 2,
            openms::chemistry::IonSeries::Y => after.ordinal >= 3,
        };
        near(
            after.mz - before.mz,
            if retained {
                12.3456789 / f64::from(after.charge)
            } else {
                0.
            },
            1e-12,
        );
    }
    unavailable(seq("ABK").fragment_ions(1));
    assert!(modified.fragment_ions(0).is_err());
}

#[test]
fn tag_setters_are_atomic_and_do_not_mutate_the_registry() {
    let db = ModificationsDB::global();
    let records = db.len();
    let mut sequence = seq("ACMK");
    sequence.set_mass_tag(1, "+12.3456789").unwrap();
    let before = sequence.clone();
    for invalid in [
        "", "NaN", "+inf", "1e3", "--1", "1.2.3", " 12", "+-1", "-999",
    ] {
        assert!(sequence.set_mass_tag(1, invalid).is_err(), "{invalid}");
        assert_eq!(sequence, before);
    }
    assert!(sequence.set_mass_tag(99, "+1").is_err());
    assert_eq!(sequence, before);
    sequence.set_n_terminal_mass_tag("+123.456789").unwrap();
    sequence.set_c_terminal_mass_tag("234.567891").unwrap();
    assert_eq!(seq(&sequence.to_string()), sequence);
    let before = sequence.clone();
    assert!(sequence.set_n_terminal_mass_tag("-9999").is_err());
    assert_eq!(sequence, before);
    sequence.set_n_terminal_modification("").unwrap();
    sequence.set_c_terminal_modification("").unwrap();
    sequence.set_modification(1, "").unwrap();
    assert_eq!(sequence, seq("ACMK"));
    sequence.set_mass_tag(2, "+15.99").unwrap();
    assert!(matches!(
        sequence.residue_modification(2).unwrap(),
        Some(SequenceModification::Known(_))
    ));
    assert_eq!(db.len(), records);
}

#[test]
fn malformed_and_nonfinite_tags_are_rejected() {
    for text in [
        "PEP[]TIDE",
        "PEP[+]TIDE",
        "PEP[.]TIDE",
        "PEP[1e3]TIDE",
        "PEP[1.0e3]TIDE",
        "PEP[NaN]TIDE",
        "PEP[inf]TIDE",
        "PEP[+1",
        "PEP[ 1]TIDE",
        "PEP[1\n]TIDE",
        "PEP[1[2]]TIDE",
        "[1]",
        "n[1]",
        "[+1].",
        "PEP.[+1]TIDE",
        "X[+0]",
        "X[-0]",
        "B[+1]",
        "Z[-1]",
        "A[-999]",
        "A[+1][+2]",
        "A(Oxidation)[+1]",
    ] {
        assert!(AASequence::parse(text).is_err(), "accepted {text}");
    }
    let overflowing = format!("A[{}]", "9".repeat(400));
    assert!(AASequence::parse(&overflowing).is_err());
    let underflowing = format!("A[0.{}1]", "0".repeat(400));
    assert!(AASequence::parse(&underflowing).is_err());
    assert!(AASequence::parse(&format!("A[{}]", "0".repeat(1025))).is_err());
    for text in ["b", "z", "x", "A B", "AB*", "α"] {
        assert!(AASequence::parse(text).is_err(), "{text}");
    }
}

#[test]
fn numeric_residue_tags_do_not_move_to_termini_after_slicing_or_setters() {
    // Native correction: C++ reclassifies a boundary Q[111] as N-terminal
    // pyro-Glu, making immutable anonymous annotations context dependent.
    let internal = seq("AQ[111]AR");
    let beginning = internal.subsequence(1..4).unwrap();
    let single = internal.subsequence(1..2).unwrap();
    let ending = internal.subsequence(0..2).unwrap();
    for sequence in [beginning, single, ending] {
        assert!(sequence.n_terminal_modification().is_none());
        assert!(sequence.c_terminal_modification().is_none());
        assert_eq!(seq(&sequence.to_string()), sequence);
        assert_eq!(seq(&sequence.to_unimod_string().unwrap()), sequence);
        unavailable(sequence.formula());
    }
    let mut sequence = seq("A");
    sequence.set_mass_tag(0, "+42").unwrap();
    // Integer +42 resolves the source's first Anywhere hit, Ala->Xle.
    assert!(sequence.residue_modification(0).unwrap().is_some());
    assert_eq!(seq(&sequence.to_string()), sequence);
    sequence.set_mass_tag(0, "+42.010565").unwrap();
    assert_eq!(sequence.to_string(), "A[+42.010565]");
    assert!(
        sequence
            .residue_modification(0)
            .unwrap()
            .unwrap()
            .mass_tag()
            .is_some()
    );
    assert!(sequence.n_terminal_modification().is_none());
    assert_eq!(seq(&sequence.to_string()), sequence);
    assert!(
        seq("n[+42]A")
            .n_terminal_modification()
            .unwrap()
            .known()
            .is_some()
    );
    assert!(
        seq("A.[-1]")
            .c_terminal_modification()
            .unwrap()
            .known()
            .is_some()
    );
    assert!(
        seq("C[143]PEPTIDEK")
            .residue_modification(0)
            .unwrap()
            .unwrap()
            .mass_tag()
            .is_some()
    );
    assert_eq!(
        seq("n[-17]QPEPTIDEK")
            .n_terminal_modification()
            .unwrap()
            .name(),
        "Gln->pyro-Glu"
    );
}
