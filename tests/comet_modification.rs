// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

//! `ANALYSIS/ID/CometModification.h`.
//!
//! Every section of the upstream class test
//! `src/tests/class_tests/openms/source/CometModification_test.cpp` at the core
//! pin `bc9cc12`, with its literals. Tests marked "native" check port-only
//! behaviour derived from the source code.

use openms::analysis::id::comet_modification::CometModification;
use openms::chemistry::{ModificationsDB, TermSpecificity};

fn real_similar(a: f64, b: f64) {
    let absdiff = (a - b).abs();
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{a} !~ {b} (absdiff {absdiff}, ratio {ratio})"
    );
}

/// An entry of the class test's shape: mass 43.0058 unless given.
fn entry(residues: &str, term_distance: i32, nc_term: i32) -> CometModification {
    CometModification {
        mass: 43.0058,
        residues: residues.into(),
        term_distance,
        nc_term,
        ..CometModification::default()
    }
}

#[test]
fn default_entry() {
    // CometModification(); the struct's default member initialisers.
    let d = CometModification::default();
    assert_eq!(d.mass, 0.0);
    assert!(d.residues.is_empty());
    assert_eq!(
        (
            d.binary_group,
            d.max_mods_per_peptide,
            d.term_distance,
            d.nc_term,
            d.required
        ),
        (0, 5, -1, 0, false)
    );
    assert_eq!(CometModification::MASS_TOLERANCE, 1e-6);
}

#[test]
fn construction_from_modifications() {
    let db = ModificationsDB::global();
    let ox = db
        .get_modification("Oxidation", Some('M'), Some(TermSpecificity::Anywhere))
        .unwrap();
    let cm = CometModification::from_modification(ox, 0, 3);
    real_similar(cm.mass, 15.994915);
    assert_eq!(cm.residues, "M");
    assert_eq!(cm.binary_group, 0);
    assert_eq!(cm.max_mods_per_peptide, 3);
    assert_eq!(cm.term_distance, -1);
    assert_eq!(cm.nc_term, 0);
    assert!(!cm.required);

    let ac = db
        .get_modification("Acetyl", None, Some(TermSpecificity::NTerm))
        .unwrap();
    let cm = CometModification::from_modification(ac, 0, 5);
    real_similar(cm.mass, 42.010565);
    assert_eq!(cm.residues, "n");
    assert_eq!((cm.term_distance, cm.nc_term), (0, 2));

    let am = db
        .get_modification("Amidated", None, Some(TermSpecificity::CTerm))
        .unwrap();
    let cm = CometModification::from_modification(am, 0, 5);
    assert_eq!(cm.residues, "c");
    assert_eq!((cm.term_distance, cm.nc_term), (0, 3));

    assert_eq!(
        CometModification::from_modification(ox, 2, 5).binary_group,
        2
    );
}

#[test]
fn protein_terminal_construction() {
    // native, from the constructor's protein-terminal branches: a protein
    // N-term modification on a specific residue keeps that residue.
    let db = ModificationsDB::global();
    let ac = db
        .get_modification("Acetyl", None, Some(TermSpecificity::ProteinNTerm))
        .unwrap();
    let cm = CometModification::from_modification(ac, 0, 1);
    assert_eq!((cm.term_distance, cm.nc_term), (0, 0));
    let expected = match ac.origin() {
        None | Some('X') => "n".to_owned(),
        Some(residue) => residue.to_string(),
    };
    assert_eq!(cm.residues, expected);
}

#[test]
fn mergeability() {
    assert!(entry("K", -1, 0).is_mergeable_with(&entry("R", -1, 0)));

    let mut m = entry("M", -1, 0);
    m.mass = 15.9949;
    assert!(!entry("K", -1, 0).is_mergeable_with(&m));

    let mut a = entry("K", -1, 0);
    a.binary_group = 1;
    let mut b = entry("R", -1, 0);
    b.binary_group = 2;
    assert!(!a.is_mergeable_with(&b));

    assert!(entry("n", 0, 2).is_mergeable_with(&entry("K", -1, 0)));
    assert!(entry("c", 0, 3).is_mergeable_with(&entry("K", -1, 0)));
    assert!(!entry("n", 0, 2).is_mergeable_with(&entry("c", 0, 3)));
    assert!(!entry("n", 0, 0).is_mergeable_with(&entry("K", -1, 0)));
    assert!(!entry("c", 0, 1).is_mergeable_with(&entry("K", -1, 0)));
    assert!(entry("n", 0, 0).is_mergeable_with(&entry("K", 0, 0)));
    assert!(!entry("n", 0, 0).is_mergeable_with(&entry("c", 0, 1)));

    let mut near = entry("R", -1, 0);
    near.mass = 43.0058 + 5e-7;
    assert!(entry("K", -1, 0).is_mergeable_with(&near));
    let mut far = entry("R", -1, 0);
    far.mass = 43.0058 + 2e-6;
    assert!(!entry("K", -1, 0).is_mergeable_with(&far));
}

#[test]
fn merging() {
    let mut a = entry("K", -1, 0);
    a.max_mods_per_peptide = 3;
    let mut b = entry("R", -1, 0);
    b.max_mods_per_peptide = 5;
    a.merge(&b);
    assert_eq!(a.residues, "KR");
    assert_eq!((a.term_distance, a.nc_term), (-1, 0));
    assert_eq!(a.max_mods_per_peptide, 5);

    let mut a = entry("K", -1, 0);
    a.merge(&entry("R", -1, 0));
    a.merge(&entry("S", -1, 0));
    assert_eq!(a.residues, "KRS");

    let mut nterm = entry("n", 0, 2);
    nterm.merge(&entry("K", -1, 0));
    assert_eq!(nterm.residues, "nK");
    assert_eq!((nterm.term_distance, nterm.nc_term), (-1, 0));

    let mut aa = entry("K", -1, 0);
    aa.merge(&entry("n", 0, 2));
    assert_eq!(aa.residues, "Kn");
    assert_eq!((aa.term_distance, aa.nc_term), (-1, 0));

    let mut nterm = entry("n", 0, 2);
    for residue in ["R", "S", "T"] {
        nterm.merge(&entry(residue, -1, 0));
    }
    assert_eq!(nterm.residues, "nRST");
    assert_eq!((nterm.term_distance, nterm.nc_term), (-1, 0));

    let mut cterm = entry("c", 0, 3);
    cterm.merge(&entry("K", -1, 0));
    assert_eq!(cterm.residues, "cK");
    assert_eq!((cterm.term_distance, cterm.nc_term), (-1, 0));

    let mut a = entry("KR", -1, 0);
    a.merge(&entry("K", -1, 0));
    assert_eq!(a.residues, "KR");

    let mut a = entry("K", -1, 0);
    let mut b = entry("R", -1, 0);
    b.required = true;
    a.merge(&b);
    assert!(a.required);
}

#[test]
fn comet_strings() {
    let mut m = CometModification {
        mass: 15.9949,
        residues: "M".into(),
        max_mods_per_peptide: 3,
        ..CometModification::default()
    };
    assert_eq!(
        m.to_comet_string(1),
        "variable_mod01 = 15.9949 M 0 3 -1 0 0 0.0"
    );
    let sty = CometModification {
        mass: 79.966331,
        residues: "STY".into(),
        max_mods_per_peptide: 3,
        ..CometModification::default()
    };
    assert!(sty.to_comet_string(9).starts_with("variable_mod09"));
    assert!(sty.to_comet_string(10).starts_with("variable_mod10"));
    let nk = CometModification {
        mass: 42.010565,
        residues: "nK".into(),
        max_mods_per_peptide: 3,
        ..CometModification::default()
    };
    assert_eq!(
        nk.to_comet_string(1),
        "variable_mod01 = 42.0106 nK 0 3 -1 0 0 0.0"
    );
    m.required = true;
    assert_eq!(
        m.to_comet_string(1),
        "variable_mod01 = 15.9949 M 0 3 -1 0 1 0.0"
    );
}

#[test]
fn comet_string_number_formatting() {
    // native: std::ostream defaults (six significant digits, no padding
    // beyond setw(2) on the index).
    let mut m = CometModification {
        mass: -17.026549,
        residues: "Q".into(),
        ..CometModification::default()
    };
    assert_eq!(
        m.to_comet_string(123),
        "variable_mod123 = -17.0265 Q 0 5 -1 0 0 0.0"
    );
    m.mass = 0.0;
    assert_eq!(m.to_comet_string(0), "variable_mod00 = 0 Q 0 5 -1 0 0 0.0");
}

#[test]
fn merging_lists() {
    let merged = CometModification::merge_modifications(&[
        entry("K", -1, 0),
        entry("R", -1, 0),
        entry("S", -1, 0),
    ]);
    assert_eq!(merged.len(), 1);
    assert_eq!(merged[0].residues, "KRS");

    let mut m = entry("M", -1, 0);
    m.mass = 15.9949;
    assert_eq!(
        CometModification::merge_modifications(&[entry("K", -1, 0), m.clone()]).len(),
        2
    );

    let merged = CometModification::merge_modifications(&[
        entry("n", 0, 0),
        entry("n", 0, 2),
        entry("K", -1, 0),
        entry("R", -1, 0),
        entry("S", -1, 0),
    ]);
    assert_eq!(merged.len(), 2);
    assert_eq!(merged[0].residues, "n");
    assert_eq!((merged[0].term_distance, merged[0].nc_term), (0, 0));
    for residue in ['n', 'K', 'R', 'S'] {
        assert!(merged[1].residues.contains(residue));
    }
    assert_eq!((merged[1].term_distance, merged[1].nc_term), (-1, 0));

    assert!(CometModification::merge_modifications(&[]).is_empty());
    let single = CometModification::merge_modifications(&[m]);
    assert_eq!(single.len(), 1);
    assert_eq!(single[0].residues, "M");
}

#[test]
fn a_merged_terminal_entry_accepts_the_opposite_terminal() {
    // native: documents a source behaviour. After a peptide N-term entry
    // merges with an amino-acid entry it has term_distance -1, so it is no
    // longer "peptide N-term" and a later peptide C-term entry of the same
    // mass merges into it, giving a residue string with both 'n' and 'c'.
    let merged = CometModification::merge_modifications(&[
        entry("n", 0, 2),
        entry("K", -1, 0),
        entry("c", 0, 3),
    ]);
    assert_eq!(merged.len(), 1);
    assert_eq!(merged[0].residues, "nKc");
}
