// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `CHEMISTRY/ElementDB.h` and the `Element.h` value members.
//!
//! The lookup and ordering cases are `ElementDB_test.cpp` and
//! `Element_test.cpp` at core `bc9cc12` (tier 2). The isotope elements, the
//! `D`/`T` aliases and the stream text are derived from `ElementDB.cpp:565-660`
//! (`storeIsotopes_`, the alias block) and `Element.cpp:147-163` at the same
//! pin (tier 4).

use openms::chemistry::element_db::{
    atomic_numbers, get_element, get_element_by_atomic_number, has_atomic_number, has_element,
    names, symbols,
};
use openms::chemistry::{EmpiricalFormula, element, element_table};
use std::collections::HashSet;
use std::ptr;

#[test]
fn lookups_return_the_same_element() {
    let carbon = get_element("Carbon").unwrap();
    assert!(ptr::eq(carbon, names()["Carbon"]));
    assert!(ptr::eq(carbon, symbols()["C"]));
    assert!(ptr::eq(carbon, atomic_numbers()[&6]));
    assert!(ptr::eq(carbon, get_element_by_atomic_number(6).unwrap()));
    assert!(ptr::eq(
        get_element("Hydrogen").unwrap(),
        get_element("H").unwrap()
    ));
    assert!(has_element("Carbon"));
    assert!(has_atomic_number(6));
    assert!(!has_element("Carbonium"));
    assert!(!has_atomic_number(0));
    assert!(get_element_by_atomic_number(200).is_none());
    // The natural table and the database agree.
    assert!(ptr::eq(element("C").unwrap(), carbon));
    assert_eq!(atomic_numbers().len(), element_table().len());
}

#[test]
fn isotope_elements_and_aliases() {
    let c13 = get_element("(13)C").unwrap();
    assert!(ptr::eq(c13, get_element("(13)Carbon").unwrap()));
    assert_eq!(c13.name(), "(13)Carbon");
    assert_eq!(c13.atomic_number(), 6);
    assert!(c13.is_isotope());
    assert!(!get_element("C").unwrap().is_isotope());
    assert_eq!(c13.mono_mass(), 13.003355000000001);
    assert_eq!(c13.average_mass(), c13.mono_mass());
    let distribution = c13.isotope_distribution();
    assert_eq!(distribution.len(), 1);
    assert_eq!(distribution.peaks()[0].probability, 1.0);
    // Isotope elements are not in the atomic-number map.
    assert!(ptr::eq(atomic_numbers()[&6], get_element("C").unwrap()));

    // The source adds D and T as symbols of (2)H and (3)H, including 3H with
    // zero natural abundance.
    assert!(ptr::eq(
        get_element("D").unwrap(),
        get_element("(2)H").unwrap()
    ));
    assert!(ptr::eq(
        get_element("T").unwrap(),
        get_element("(3)Hydrogen").unwrap()
    ));
    assert!(names().get("D").is_none());

    // One isotope element per tabulated isotope.
    let tabulated: usize = element_table().iter().map(|e| e.isotopes().len()).sum();
    assert_eq!(names().len(), element_table().len() + tabulated);
    assert_eq!(symbols().len(), element_table().len() + tabulated + 2);
}

#[test]
fn formulas_key_isotope_elements_as_labelled_atoms() {
    let c13 = get_element("(13)C").unwrap();
    let labelled = EmpiricalFormula::from_element(2, c13, 0);
    assert_eq!(labelled, EmpiricalFormula::parse("(13)C2").unwrap());
    assert_eq!(labelled.number_of(c13), 2);
    assert!(labelled.has_element(c13));
    assert!(!labelled.has_element(get_element("C").unwrap()));
    let heavy_water = EmpiricalFormula::parse("D2O").unwrap();
    assert_eq!(heavy_water.number_of(get_element("D").unwrap()), 2);
}

/// `bool operator < (const Element& element) const`.
#[test]
fn ordering_matches_the_source_cases() {
    let h = get_element("H").unwrap();
    let c = get_element("Carbon").unwrap();
    let o = get_element("O").unwrap();
    let s = get_element("S").unwrap();
    assert!(h < c);
    assert!(c < o);
    assert!(c >= c);
    assert!(s >= c);
    // Same atomic number: the monoisotopic weight decides.
    assert!(c < get_element("(13)C").unwrap());
}

/// `[EXTRA] std::hash<Element>`: equal elements hash equally, different ones
/// differently.
#[test]
fn hashes_follow_equality() {
    let carbon = *get_element("C").unwrap();
    let copy = carbon;
    assert_eq!(carbon, copy);
    assert_eq!(carbon.source_hash64(), copy.source_hash64());
    let nitrogen = *get_element("N").unwrap();
    assert_ne!(carbon, nitrogen);
    assert_ne!(carbon.source_hash64(), nitrogen.source_hash64());
    let mut set = HashSet::new();
    set.insert(carbon);
    set.insert(copy);
    set.insert(nitrogen);
    assert_eq!(set.len(), 2);
    let all: HashSet<u64> = element_table().iter().map(|e| e.source_hash64()).collect();
    assert_eq!(all.len(), element_table().len());
}

#[test]
fn stream_text_follows_the_source_operator() {
    assert_eq!(
        get_element("C").unwrap().source_stream_text(),
        "Carbon C 6 12.0107 12 12.0=98.93% 13.003355000000001=1.07%"
    );
    // Zero-abundance isotopes are omitted: 3H.
    assert_eq!(
        get_element("H").unwrap().source_stream_text(),
        "Hydrogen H 1 1.00794 1.00783 1.0078250319=99.9885% 2.01410178=0.0115%"
    );
    assert_eq!(
        get_element("(13)C").unwrap().source_stream_text(),
        "(13)Carbon (13)C 6 13.0034 13.0034 13.003355000000001=100%"
    );
}
