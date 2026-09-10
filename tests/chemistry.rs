// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS contributors $
// Golden expectations derive from OpenMS4-core revision 7c029e8 class tests:
// EmpiricalFormula_test.cpp, AASequence_test.cpp, and ProteaseDigestion_test.cpp.

use openms::chemistry::{
    AASequence, EmpiricalFormula, IonSeries, PROTON_MASS_U, Protease, ProteaseDigestion, element,
    element_table,
};

fn formula(text: &str) -> EmpiricalFormula {
    EmpiricalFormula::parse(text).unwrap()
}

fn sequence(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.12}, expected {expected:.12}, tolerance {tolerance}"
    );
}

#[test]
fn formula_upstream_isotope_mass_and_aliases() {
    // EmpiricalFormula_test.cpp: C5(13)C4 equals (12)C5(13)C4 by mass.
    let labeled = formula("C5(13)C4\n\n ");
    near(labeled.mono_mass(), 112.013419, 0.000002);
    near(
        labeled.mono_mass(),
        formula("(12)C5(13)C4").mono_mass(),
        1e-12,
    );
    assert_eq!(labeled.count("C").unwrap(), 5);
    assert_eq!(labeled.count("(13)C").unwrap(), 4);
    assert_eq!(labeled.count("N").unwrap(), 0);
    assert_eq!(labeled.atom_count(), 9);
    assert_eq!(formula("D2O"), formula("(2)H2O"));
    assert_eq!(formula("T2O"), formula("(3)H2O"));
    assert!(formula("C5H2").mono_mass() < formula("C5D2").mono_mass());
    assert!(formula("C5D2").mono_mass() < formula("C5T2").mono_mass());
    near(formula("(15)N").average_mass(), 15.000109, 1e-12);
}

#[test]
fn formula_charge_is_protonation_not_electron_ionization() {
    for (text, charge) in [("C2+", 1), ("C2+3", 3), ("C2-2", -2)] {
        let charged = formula(text);
        assert_eq!(charged.charge(), charge);
        near(
            charged.mono_mass(),
            24.0 + f64::from(charge) * PROTON_MASS_U,
            1e-12,
        );
        near(
            charged.average_mass(),
            formula("C2").average_mass() + f64::from(charge) * PROTON_MASS_U,
            1e-12,
        );
    }
    assert_eq!(formula("H-2").count("H").unwrap(), -2);
    assert_eq!(formula("H-2").charge(), 0);
    assert_eq!(formula("H1-2").count("H").unwrap(), 1);
    assert_eq!(formula("H1-2").charge(), -2);
    assert_eq!(formula("H-2-2").count("H").unwrap(), -2);
    assert_eq!(formula("H-2-2").charge(), -2);
    for charge in [-2, -1, 0, 1, 2, i32::MIN, i32::MAX] {
        let ion = formula("C1H2").with_charge(charge);
        assert_eq!(formula(&ion.to_string()), ion);
    }
    assert_eq!(formula("+2").charge(), 2);
    assert_eq!(formula("-").charge(), -1);
    assert!(formula("+2").is_empty());
    assert_eq!(formula("").mono_mass(), 0.0);
    assert_eq!(formula("  \n\t "), EmpiricalFormula::default());
    near(
        formula("C6H12O6-2").mz().unwrap(),
        (180.0633903828 - 2.0 * PROTON_MASS_U) / 2.0,
        1e-10,
    );
    assert!(formula("C").mz().is_err());
    assert!(formula("C-10+1").mz().is_err());
}

#[test]
fn formula_upstream_algebra_and_contains() {
    assert_eq!(formula("C3H8").checked_scale(3).unwrap(), formula("C9H24"));
    assert_eq!(
        formula("C6").checked_add(&formula("C-6H2")).unwrap(),
        formula("H2")
    );
    assert_eq!(
        formula("C5H12").checked_sub(&formula("CH12")).unwrap(),
        formula("C4")
    );
    assert_eq!(
        formula("C4").checked_sub(&formula("C4H-2")).unwrap(),
        formula("H2")
    );
    assert!(formula("C6C-6").is_empty());
    assert!(formula("C0H0").is_empty());
    assert_eq!(
        formula("C2H4+2").checked_scale(-2).unwrap(),
        formula("C-4H-8-4")
    );
    assert_eq!(
        formula("C2+2").checked_sub(&formula("C1+1")).unwrap(),
        formula("C1+1")
    );
    assert_eq!(
        formula("C2H4+2").checked_scale(0).unwrap(),
        EmpiricalFormula::default()
    );
    let metabolite = formula("C12H36N2");
    for part in [
        "C12H36N2",
        "C-12H36N2",
        "C11H36N2",
        "N2",
        "H3",
        "P-1",
        "",
        "K-1H2",
    ] {
        assert!(metabolite.contains(&formula(part)), "{part}");
    }
    assert!(!metabolite.contains(&formula("P1")));
    assert!(!metabolite.contains(&formula("KH-2")));
}

#[test]
fn formula_rejects_malformed_input_and_checks_overflow() {
    for text in [
        "2H",
        "C6 H12",
        "C6\nH12",
        "h2o",
        "Xx",
        "(14)C",
        "(13C",
        "(C)",
        "(13)",
        "(65536)C",
        "C(OH)2",
        "[13C]",
        "C++",
        "C+2H",
        "C--2",
        "C+-2",
        "C+2147483648",
        "C2147483648",
        "C-2147483649",
        "Cé",
        "💧",
    ] {
        assert!(EmpiricalFormula::parse(text).is_err(), "accepted {text:?}");
    }
    assert!(formula("C2147483647").checked_add(&formula("C1")).is_err());
    assert!(formula("C-2147483648").checked_sub(&formula("C1")).is_err());
    assert!(formula("C-2147483648").checked_scale(-1).is_err());
    assert!(formula("+2147483647").checked_add(&formula("+1")).is_err());
    assert!(EmpiricalFormula::parse("C2147483647C1").is_err());
    let small = formula("C1");
    for text in ["", "C2", "CC", "(13)C2", "Xx"] {
        assert!(small.count(text).is_err());
    }
}

#[test]
fn whole_element_table_is_available_and_most_abundant_isotope_is_used() {
    assert_eq!(element_table().len(), 84);
    assert_eq!(element("C"), element("Carbon"));
    let carbon = element("C").unwrap();
    assert_eq!(carbon.name(), "Carbon");
    assert_eq!(carbon.symbol(), "C");
    assert_eq!(carbon.atomic_number(), 6);
    assert!(element("Unknown").is_none());
    for entry in element_table() {
        assert!(entry.mono_mass().is_finite() && entry.mono_mass() > 0.0);
        assert!(entry.average_mass().is_finite() && entry.average_mass() > 0.0);
        near(
            formula(entry.symbol()).mono_mass(),
            entry.mono_mass(),
            1e-12,
        );
        for isotope in entry.isotopes() {
            let label = format!("({}){}", isotope.mass_number, entry.symbol());
            near(formula(&label).mono_mass(), isotope.mass, 1e-12);
            near(formula(&label).average_mass(), isotope.mass, 1e-12);
        }
    }
    // AASequence_test.cpp explicitly checks Se uses isotope 80, not lightest 74.
    near(formula("Se").mono_mass(), 79.9165213, 1e-10);
    // Intentional correction: declared Ir tables are used, not upstream's
    // accidental call to buildElement_("Iridium", "Ir", ..., rhenium_...).
    near(formula("Ir").mono_mass(), 192.962924, 1e-10);
    near(formula("(191)Ir").mono_mass(), 190.960591, 1e-10);
    assert!(EmpiricalFormula::parse("(187)Ir").is_err());
}

#[test]
fn peptide_upstream_formula_mono_and_average_masses() {
    assert_eq!(
        sequence("ACDEF").formula().unwrap(),
        formula("O10SH33N5C24")
    );
    let peptide = sequence("DFPIANGER");
    // The historical upstream literal predates the rounded ElementDB table;
    // it differs by 5.6 microdaltons and passes the upstream relative tolerance.
    near(peptide.mono_mass().unwrap(), 1017.487958568, 6e-6);
    assert_eq!(peptide.formula().unwrap(), formula("C44H67N13O15"));
    near(peptide.mono_mass().unwrap(), 1017.4879641373, 1e-9);
    near(peptide.average_mass().unwrap(), 1018.08088, 0.01);
    near(peptide.mz(1).unwrap(), 1018.495240604071, 1e-9);
    near(peptide.mz(2).unwrap(), 509.751258535421, 1e-9);
    near(peptide.mz(-2).unwrap(), 507.736705601879, 1e-9);
    assert!(peptide.mz(0).is_err());
    assert!(peptide.mz(i32::MIN).is_err());
    // Selenium-containing golden formula/masses in AASequence_test.cpp.
    let seleno = formula("C75H122N20O32S2Se1");
    near(seleno.mono_mass(), 1958.7140766518, 1e-9);
    near(seleno.average_mass(), 1958.981404189803, 1e-9);
    // The upstream sequence contains one oxidation. Unmodified has one fewer O.
    let unmodified = sequence("PEPTIDESEKUEMCER");
    assert_eq!(
        unmodified
            .formula()
            .unwrap()
            .checked_add(&formula("O"))
            .unwrap(),
        seleno
    );
    assert_eq!(sequence("O").formula().unwrap(), formula("C12H21N3O3"));
    assert_eq!(
        sequence("I").formula().unwrap(),
        sequence("L").formula().unwrap()
    );
    assert_eq!(
        sequence("J").formula().unwrap(),
        sequence("I").formula().unwrap()
    );
    assert_eq!(sequence("U").formula().unwrap(), formula("C3H7NO2Se"));
}

#[test]
fn every_supported_residue_matches_its_upstream_full_formula() {
    // Independently transcribed full amino-acid formulas from ResidueDB.cpp.
    for (residue, expected) in [
        ("A", "C3H7NO2"),
        ("C", "C3H7NO2S"),
        ("D", "C4H7NO4"),
        ("E", "C5H9NO4"),
        ("F", "C9H11NO2"),
        ("G", "C2H5NO2"),
        ("H", "C6H9N3O2"),
        ("I", "C6H13NO2"),
        ("K", "C6H14N2O2"),
        ("L", "C6H13NO2"),
        ("M", "C5H11NO2S"),
        ("N", "C4H8N2O3"),
        ("P", "C5H9NO2"),
        ("Q", "C5H10N2O3"),
        ("R", "C6H14N4O2"),
        ("S", "C3H7NO3"),
        ("T", "C4H9NO3"),
        ("V", "C5H11NO2"),
        ("W", "C11H12N2O2"),
        ("Y", "C9H11NO3"),
        ("U", "C3H7NO2Se"),
        ("O", "C12H21N3O3"),
        ("J", "C6H13NO2"),
    ] {
        assert_eq!(
            sequence(residue).formula().unwrap(),
            formula(expected),
            "{residue}"
        );
    }
}

#[test]
fn peptide_parsing_subsequence_and_empty_boundaries() {
    let peptide = sequence("DFPIANGER");
    assert_eq!(peptide.len(), 9);
    assert_eq!(peptide.as_str(), "DFPIANGER");
    assert_eq!(peptide.subsequence(3..6).unwrap(), sequence("IAN"));
    assert_eq!(peptide.subsequence(0..3).unwrap(), sequence("DFP"));
    assert_eq!(peptide.subsequence(6..9).unwrap(), sequence("GER"));
    assert!(peptide.subsequence(0..10).is_err());
    let (start, end) = (4, 3);
    assert!(peptide.subsequence(start..end).is_err());
    assert_eq!(peptide.subsequence(9..9).unwrap(), AASequence::default());
    let empty = sequence("");
    assert!(empty.is_empty());
    assert_eq!(empty.mono_mass().unwrap(), 0.0);
    assert!(empty.formula().unwrap().is_empty());
    assert!(empty.mz(1).is_err());
    for text in ["peptide", "PEP TIDE", "PEP\nTIDE", "PEPé", "*", "-", "123"] {
        assert!(AASequence::parse(text).is_err(), "accepted {text:?}");
    }
}

#[test]
fn fragment_ions_use_upstream_b_y_formulas_and_charge() {
    // AASequence_test.cpp single-residue b/y alanine formulas are exercised
    // as b1/y1 of AA, giving a real internal peptide bond.
    let ions = sequence("AA").fragment_ions(2).unwrap();
    assert_eq!(ions.len(), 4);
    assert_eq!(
        (ions[0].series, ions[0].ordinal, ions[0].charge),
        (IonSeries::B, 1, 1)
    );
    assert_eq!(
        (ions[1].series, ions[1].ordinal, ions[1].charge),
        (IonSeries::Y, 1, 1)
    );
    near(ions[0].mz, 72.044390626271, 1e-10);
    near(ions[1].mz, 90.054955690071, 1e-10);
    near(ions[2].mz, (ions[0].mz + PROTON_MASS_U) / 2.0, 1e-12);
    near(ions[3].mz, (ions[1].mz + PROTON_MASS_U) / 2.0, 1e-12);
    let peptide = sequence("PEPTIDE");
    let fragments = peptide.fragment_ions(3).unwrap();
    assert_eq!(fragments.len(), 36);
    for pair in fragments.chunks_exact(2) {
        let charge = f64::from(pair[0].charge);
        assert_eq!(pair[0].ordinal + pair[1].ordinal, peptide.len());
        near(
            (pair[0].mz + pair[1].mz) * charge - 2.0 * charge * PROTON_MASS_U,
            peptide.mono_mass().unwrap(),
            1e-9,
        );
    }
    assert!(sequence("").fragment_ions(1).unwrap().is_empty());
    assert!(sequence("A").fragment_ions(1).unwrap().is_empty());
    assert!(sequence("AA").fragment_ions(0).is_err());
}

fn digest_strings(digest: &ProteaseDigestion, input: &str) -> Vec<String> {
    digest
        .digest(&sequence(input))
        .unwrap()
        .into_iter()
        .map(|p| p.sequence.to_string())
        .collect()
}

#[test]
fn trypsin_matches_upstream_zero_and_missed_cleavage_products() {
    let mut digest = ProteaseDigestion::default();
    for (input, expected) in [
        ("ACDE", vec!["ACDE"]),
        ("ACKDE", vec!["ACK", "DE"]),
        ("ACRDE", vec!["ACR", "DE"]),
        ("ARCRDRE", vec!["AR", "CR", "DR", "E"]),
        ("RKR", vec!["R", "K", "R"]),
    ] {
        assert_eq!(digest_strings(&digest, input), expected);
    }
    digest.missed_cleavages = 1;
    assert_eq!(
        digest_strings(&digest, "ARCDRE"),
        ["AR", "CDR", "E", "ARCDR", "CDRE"]
    );
    assert_eq!(digest_strings(&digest, "RKR"), ["R", "K", "R", "RK", "KR"]);
    let products = digest.digest(&sequence("ARCDRE")).unwrap();
    let spans: Vec<_> = products
        .iter()
        .map(|p| (p.start, p.end, p.missed_cleavages))
        .collect();
    assert_eq!(
        spans,
        [(0, 2, 0), (2, 5, 0), (5, 6, 0), (0, 5, 1), (2, 6, 1)]
    );
    digest.missed_cleavages = usize::MAX;
    assert_eq!(
        digest_strings(&digest, "RKR"),
        ["R", "K", "R", "RK", "KR", "RKR"]
    );
    assert_eq!(digest_strings(&digest, "ACDE"), ["ACDE"]);
}

#[test]
fn trypsin_proline_boundary_and_length_filters() {
    let mut digest = ProteaseDigestion::default();
    for input in ["ACKPDE", "ACRPDE", "WKP", "MRP"] {
        assert_eq!(digest_strings(&digest, input), [input]);
    }
    assert_eq!(digest_strings(&digest, "ACKPDERA"), ["ACKPDER", "A"]);
    assert_eq!(digest_strings(&digest, "ACRPDEKA"), ["ACRPDEK", "A"]);
    digest.enzyme = Protease::TrypsinP;
    assert_eq!(digest_strings(&digest, "ACKPDE"), ["ACK", "PDE"]);
    assert_eq!(digest_strings(&digest, "ACRPDE"), ["ACR", "PDE"]);
    digest.enzyme = Protease::Trypsin;
    digest.missed_cleavages = 1;
    digest.min_length = 3;
    digest.max_length = Some(4);
    // Exact upstream filtering golden (three candidates discarded).
    assert_eq!(digest_strings(&digest, "ARCDRE"), ["CDR", "CDRE"]);
    digest.min_length = 4;
    assert_eq!(digest_strings(&digest, "ARCDRE"), ["CDRE"]);
    digest.min_length = 5;
    assert!(digest.digest(&sequence("ARCDRE")).is_err());
    digest.min_length = 0;
    assert!(digest.digest(&sequence("ARCDRE")).is_err());
}

#[test]
fn digestion_empty_terminal_sites_no_cleavage_and_duplicate_sequences() {
    let mut digest = ProteaseDigestion::default();
    assert!(digest.digest(&sequence("")).unwrap().is_empty());
    assert_eq!(digest_strings(&digest, "R"), ["R"]);
    assert_eq!(digest_strings(&digest, "RR"), ["R", "R"]);
    assert_eq!(digest_strings(&digest, "RK"), ["R", "K"]);
    assert_eq!(digest_strings(&digest, "PRP"), ["PRP"]);
    digest.enzyme = Protease::NoCleavage;
    digest.missed_cleavages = usize::MAX;
    assert_eq!(digest_strings(&digest, "ARCDRE"), ["ARCDRE"]);
    digest.max_length = Some(5);
    assert!(digest.digest(&sequence("ARCDRE")).unwrap().is_empty());
}
