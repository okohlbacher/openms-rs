// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Source expectations: pinned AASequence/ModificationsDB/Residue class tests.

use openms::chemistry::{
    AASequence, EmpiricalFormula, ModificationsDB, PROTON_MASS_U, ProteaseDigestion,
    TermSpecificity,
};

fn peptide(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}
fn formula(text: &str) -> EmpiricalFormula {
    EmpiricalFormula::parse(text).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.12} differs from {b:.12}");
}

#[test]
fn entire_pinned_registry_loads_and_resolves_specificities() {
    let db = ModificationsDB::global();
    assert_eq!(db.len(), 3035 + 92); // UniMod/custom, then pinned XLMOD monolinks.
    assert!(!db.is_empty());
    let oxidation = db
        .get_modification("Oxidation", Some('M'), Some(TermSpecificity::Anywhere))
        .unwrap();
    assert_eq!(oxidation.full_id(), "Oxidation (M)");
    assert_eq!(oxidation.record_id(), Some(35));
    assert_eq!(oxidation.accession(), "UniMod:35");
    assert_eq!(oxidation.diff_formula(), &formula("O"));
    near(oxidation.diff_mono_mass(), 15.994915, 1e-12);
    near(oxidation.diff_average_mass(), 15.9994, 1e-12);
    for name in [
        "UniMod:35",
        "unimod:35",
        "UNIMOD:35",
        "Oxidation (M)",
        "Oxidation or Hydroxylation",
    ] {
        assert_eq!(
            db.get_modification(name, Some('M'), Some(TermSpecificity::Anywhere))
                .unwrap(),
            oxidation
        );
    }
    assert!(db.get_modification("Oxidation", None, None).is_err());
    assert!(
        db.get_modification("Acetyl", Some('A'), Some(TermSpecificity::Anywhere))
            .is_err()
    );
    assert!(db.find("not-a-modification", None, None).is_empty());
    let acetyl = db
        .get_modification("Acetyl (N-term)", Some('P'), Some(TermSpecificity::NTerm))
        .unwrap();
    assert_eq!(acetyl.origin(), None);
    assert_eq!(acetyl.diff_formula(), &formula("C2H2O"));
}

#[test]
fn neutral_losses_and_mass_search_preserve_declared_data() {
    let db = ModificationsDB::global();
    let phospho = db
        .get_modification("Phospho", Some('S'), Some(TermSpecificity::Anywhere))
        .unwrap();
    assert!(
        phospho
            .neutral_losses()
            .iter()
            .any(|loss| loss.formula() == &formula("H3PO4")
                && (loss.mono_mass() - 97.976896).abs() < 1e-12)
    );
    let found = db
        .search_by_mass(57.021464, 1e-6, Some('C'), Some(TermSpecificity::Anywhere))
        .unwrap();
    assert!(found.iter().any(|m| m.name() == "Carbamidomethyl"));
    assert_eq!(
        db.best_by_mass(15.9949, 0.0001, Some('M'), Some(TermSpecificity::Anywhere))
            .unwrap()
            .unwrap()
            .name(),
        "Oxidation"
    );
    assert!(db.best_by_mass(1e9, 0.01, None, None).unwrap().is_none());
    for (mass, tolerance) in [(f64::NAN, 0.1), (0.0, f64::INFINITY), (0.0, -1.0)] {
        assert!(db.search_by_mass(mass, tolerance, None, None).is_err());
    }
}

#[test]
fn custom_table_parser_is_checked_and_preserves_source_terminal_wildcard() {
    let row = "1\tTest\tTest name\tC-term\tanywhere\t1\t1\tH1\t0\tOther\t";
    let db = ModificationsDB::from_tsv(row).unwrap();
    assert_eq!(db.entries()[0].full_id(), "Test (X)");
    assert!(
        db.get_modification("Test", Some('A'), Some(TermSpecificity::Anywhere))
            .is_ok()
    );
    assert_eq!(db.search_by_mass(0.0, 1.0, None, None).unwrap().len(), 1);
    assert!(db.best_by_mass(0.0, 1.0, None, None).unwrap().is_none());
    assert!(db.best_by_mass(1.0, 0.0, None, None).unwrap().is_none());
    for text in [
        "short",
        &row.replace("\t1\t1\tH1", "\tNaN\t1\tH1"),
        &row.replace("anywhere", "unknown"),
        &row.replace("H1", "Nonexistent"),
        &row.replace("\t0\tOther", "\t2\tOther"),
    ] {
        assert!(ModificationsDB::from_tsv(text).is_err(), "accepted {text}");
    }
}

#[test]
fn named_residue_and_isotope_modifications_match_source_formulas() {
    let oxidized = peptide("PEPM(Oxidation)TIDEK");
    assert_eq!(
        oxidized.formula().unwrap(),
        peptide("PEPMTIDEK")
            .formula()
            .unwrap()
            .checked_add(&formula("O"))
            .unwrap()
    );
    near(
        oxidized.mono_mass().unwrap(),
        oxidized.formula().unwrap().mono_mass(),
        1e-12,
    );
    assert_eq!(oxidized, peptide("PEPM(UniMod:35)TIDEK"));
    let cam = peptide("PEPC(Carbamidomethyl)TIDEK");
    near(
        cam.mono_mass().unwrap() - peptide("PEPCTIDEK").mono_mass().unwrap(),
        formula("C2H3NO").mono_mass(),
        1e-10,
    );
    let label = peptide("CNARCK(Label:13C(6)15N(2))NCNCN");
    assert_eq!(
        label.residue_modification(5).unwrap().unwrap().name(),
        "Label:13C(6)15N(2)"
    );
    assert_eq!(label.formula().unwrap().count("(13)C").unwrap(), 6);
    assert_eq!(label.formula().unwrap().count("(15)N").unwrap(), 2);
    // Pinned ElementDB: 6*(13.003355-12) + 2*(15.000109-14.003074).
    near(
        label.mono_mass().unwrap() - peptide(label.as_str()).mono_mass().unwrap(),
        8.0142,
        1e-10,
    );
    assert_eq!(label, peptide(&label.to_string()));
}

#[test]
fn terminal_grammar_and_accession_strings_match_source_examples() {
    let explicit = peptide(".(UniMod:1)PEPC(UniMod:4)PEPM(UniMod:35)PEPR.(UniMod:2)");
    assert_eq!(
        explicit,
        peptide("n(Acetyl)PEPC(Carbamidomethyl)PEPM(Oxidation)PEPRc(Amidated)")
    );
    assert_eq!(
        explicit.to_unimod_string().unwrap(),
        ".(UniMod:1)PEPC(UniMod:4)PEPM(UniMod:35)PEPR.(UniMod:2)"
    );
    assert_eq!(
        explicit.to_string(),
        ".(Acetyl)PEPC(Carbamidomethyl)PEPM(Oxidation)PEPR.(Amidated)"
    );
    assert_eq!(
        peptide("Q(Gln->pyro-Glu)PEPTIDEK"),
        peptide("(UniMod:28)QPEPTIDEK")
    );
    assert_eq!(
        peptide("PEPTIDEM(UniMod:10)"),
        peptide("PEPTIDEM.(UniMod:10)")
    );
    assert_eq!(peptide("M(UniMod:10)"), peptide("M.(UniMod:10)"));
    assert_eq!(peptide("A(Amidated)"), peptide("A.(Amidated)"));
    assert_eq!(
        peptide("ANLVFK(Label:13C(6)15N(2))EIEK(Label:2H(4))(Amidated)")
            .c_terminal_modification()
            .unwrap()
            .name(),
        "Amidated"
    );
    assert_eq!(
        peptide("(UniMod:51)CPEPTIDE")
            .n_terminal_modification()
            .unwrap()
            .term_specificity(),
        TermSpecificity::ProteinNTerm
    );
}

#[test]
fn terminal_masses_use_declared_deltas_and_average_uses_formula() {
    let base = peptide("DFPIANGER");
    let modified = peptide("(Acetyl)DFPIANGER(Amidated)");
    near(
        modified.mono_mass().unwrap(),
        base.mono_mass().unwrap() + 42.010565 - 0.984016,
        1e-10,
    );
    assert_eq!(
        modified.formula().unwrap(),
        base.formula()
            .unwrap()
            .checked_add(&formula("C2H3N"))
            .unwrap()
    );
    near(
        modified.average_mass().unwrap(),
        modified.formula().unwrap().average_mass(),
        1e-10,
    );
    near(
        modified.mz(2).unwrap(),
        (modified.mono_mass().unwrap() + 2.0 * PROTON_MASS_U) / 2.0,
        1e-10,
    );
    near(
        peptide("(UniMod:51)CPEPTIDE").mono_mass().unwrap(),
        902.3691545801998 + 788.725777,
        1e-6,
    );
    near(
        peptide("(dNIC)DFPIANGER").mono_mass().unwrap(),
        base.mono_mass().unwrap() + 109.048119,
        1e-6,
    );
}

#[test]
fn mutation_replaces_removes_and_is_atomic_on_errors() {
    let mut seq = peptide("CMK");
    seq.set_modification(0, "Carbamidomethyl").unwrap();
    seq.set_modification(1, "Oxidation").unwrap();
    seq.set_n_terminal_modification("Acetyl").unwrap();
    seq.set_c_terminal_modification("Amidated").unwrap();
    assert_eq!(
        seq,
        peptide("(Acetyl)C(Carbamidomethyl)M(Oxidation)K(Amidated)")
    );
    let before = seq.clone();
    assert!(seq.set_modification(50, "Oxidation").is_err());
    assert!(seq.set_modification(1, "nonexistent").is_err());
    assert!(seq.set_n_terminal_modification("Met->Hse").is_err());
    assert_eq!(seq, before);
    seq.set_modification(0, "").unwrap();
    seq.set_modification(1, "").unwrap();
    seq.set_n_terminal_modification("").unwrap();
    seq.set_c_terminal_modification("").unwrap();
    assert!(!seq.is_modified());
    assert_eq!(seq, peptide("CMK"));
    assert!(
        AASequence::default()
            .set_n_terminal_modification("Acetyl")
            .is_err()
    );
}

#[test]
fn display_retains_explicit_protein_terminal_specificity() {
    let mut seq = peptide("PEPTIDE");
    seq.set_n_terminal_modification("Acetyl (Protein N-term)")
        .unwrap();
    seq.set_c_terminal_modification("Amidated (Protein C-term)")
        .unwrap();
    assert_eq!(peptide(&seq.to_string()), seq);
    assert!(seq.to_string().contains("Protein N-term"));
}

#[test]
fn subsequences_and_digestion_retain_only_original_terminal_modifications() {
    let seq = peptide("(Acetyl)AC(Carbamidomethyl)KRPM(Oxidation)K(Amidated)");
    assert_eq!(
        seq.prefix(3).unwrap(),
        peptide("(Acetyl)AC(Carbamidomethyl)K")
    );
    assert_eq!(seq.suffix(4).unwrap(), peptide("RPM(Oxidation)K(Amidated)"));
    assert_eq!(
        seq.subsequence(1..6).unwrap(),
        peptide("C(Carbamidomethyl)KRPM(Oxidation)")
    );
    assert_eq!(seq.prefix(0).unwrap(), AASequence::default());
    assert!(seq.suffix(8).is_err());
    let products = ProteaseDigestion::default().digest(&seq).unwrap();
    assert_eq!(products.len(), 2);
    assert_eq!(products[0].sequence, seq.prefix(3).unwrap());
    assert_eq!(products[1].sequence, seq.suffix(4).unwrap());
}

#[test]
fn modified_fragment_series_follow_retained_residues_and_termini() {
    let base = peptide("ACMK");
    let modified = peptide("(Acetyl)AC(Carbamidomethyl)M(Oxidation)K(Amidated)");
    let unmodified = base.fragment_ions(2).unwrap();
    let actual = modified.fragment_ions(2).unwrap();
    for (ion, original) in actual.iter().zip(&unmodified) {
        let delta = match ion.series {
            openms::chemistry::IonSeries::B => {
                42.010565
                    + if ion.ordinal >= 2 {
                        formula("C2H3NO").mono_mass()
                    } else {
                        0.0
                    }
                    + if ion.ordinal >= 3 {
                        formula("O").mono_mass()
                    } else {
                        0.0
                    }
            }
            openms::chemistry::IonSeries::Y => {
                -0.984016
                    + if ion.ordinal >= 2 {
                        formula("O").mono_mass()
                    } else {
                        0.0
                    }
                    + if ion.ordinal >= 3 {
                        formula("C2H3NO").mono_mass()
                    } else {
                        0.0
                    }
            }
        };
        near(ion.mz - original.mz, delta / f64::from(ion.charge), 1e-10);
    }
}

#[test]
fn invalid_modification_syntax_and_specificity_are_errors() {
    for text in [
        "M()",
        "M(Oxidation",
        "M(Oxidation))",
        "M(Oxidation)(Oxidation)",
        "(Acetyl)(Acetyl)PEP",
        "PEP(Amidated)(Amidated)",
        "(Amidated)PEP",
        "PEPTIDEM(UniMod:10)K",
        "PQ(UniMod:28)EPTIDEK",
        "PC(UniMod:26)EPTIDEK",
        "PEP.(Amidated)TIDE",
        "(Acetyl)",
        "A(unknown)",
    ] {
        assert!(AASequence::parse(text).is_err(), "accepted {text}");
    }
    let deep = format!("M({}Oxidation{})", "(".repeat(130), ")".repeat(130));
    assert!(AASequence::parse(&deep).is_err());
}
