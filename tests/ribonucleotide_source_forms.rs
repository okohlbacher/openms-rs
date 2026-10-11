// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `CHEMISTRY/Ribonucleotide.h`: the source stream text and the source hash.
//!
//! The hash cases are the `[EXTRA] std::hash<Ribonucleotide>` section of
//! `Ribonucleotide_test.cpp` at core `bc9cc12` (tier 2). The stream text is
//! derived from `operator<<` in `Ribonucleotide.cpp:21-28` at the same pin and
//! the `EmpiricalFormula` `operator<<` it calls (tier 4).

use openms::chemistry::{
    EmpiricalFormula, Ribonucleotide, RibonucleotideRecord, RibonucleotideTermSpecificity,
};
use std::collections::HashMap;

#[allow(clippy::too_many_arguments)]
fn ribo(
    name: &str,
    code: &str,
    formula: &str,
    origin: char,
    mono: f64,
    average: f64,
    term: RibonucleotideTermSpecificity,
) -> Ribonucleotide {
    Ribonucleotide::from_record(RibonucleotideRecord {
        name: name.into(),
        code: code.into(),
        new_code: code.into(),
        html_code: code.into(),
        formula: EmpiricalFormula::parse(formula).unwrap(),
        origin,
        mono_mass: mono,
        average_mass: average,
        term_specificity: term,
        baseloss_formula: EmpiricalFormula::parse("C5H10O5").unwrap(),
    })
    .unwrap()
}

fn adenosine(term: RibonucleotideTermSpecificity) -> Ribonucleotide {
    ribo("adenosine", "A", "C10H13N5O4", 'A', 267.0968, 267.24, term)
}

#[test]
fn source_hash_follows_equality() {
    let anywhere = RibonucleotideTermSpecificity::Anywhere;
    let ribo1 = adenosine(anywhere);
    let ribo2 = adenosine(anywhere);
    assert_eq!(ribo1, ribo2);
    assert_eq!(ribo1.source_hash64(), ribo2.source_hash64());

    let ribo3 = ribo(
        "guanosine",
        "G",
        "C10H13N5O5",
        'G',
        283.0917,
        283.24,
        anywhere,
    );
    assert_ne!(ribo1, ribo3);
    assert_ne!(ribo1.source_hash64(), ribo3.source_hash64());

    let copy = ribo1.clone();
    assert_eq!(ribo1.source_hash64(), copy.source_hash64());

    // The source keys an unordered_map by the hash; the native Hash must agree
    // with equality in the same way.
    let mut map = HashMap::new();
    map.insert(ribo1.clone(), "first");
    map.insert(ribo3.clone(), "third");
    assert_eq!(map.len(), 2);
    map.insert(ribo2, "second");
    assert_eq!(map.len(), 2);
    assert_eq!(map[&ribo1], "second");
    assert_eq!(map[&ribo3], "third");

    assert_eq!(
        Ribonucleotide::default().source_hash64(),
        Ribonucleotide::default().source_hash64()
    );
    assert_ne!(
        adenosine(RibonucleotideTermSpecificity::FivePrime).source_hash64(),
        adenosine(RibonucleotideTermSpecificity::ThreePrime).source_hash64()
    );
}

#[test]
fn source_stream_text_uses_the_source_formula_form() {
    let anywhere = RibonucleotideTermSpecificity::Anywhere;
    assert_eq!(
        adenosine(anywhere).source_stream_text(),
        "Ribonucleotide 'A' (adenosine, C10H13N5O4)"
    );
    // A count of one prints the bare symbol in the source's formula operator<<.
    let amp = ribo(
        "adenosine monophosphate",
        "pA",
        "C10H14N5O7P",
        'A',
        347.063,
        347.22,
        anywhere,
    );
    assert_eq!(
        amp.source_stream_text(),
        "Ribonucleotide 'pA' (adenosine monophosphate, C10H14N5O7P)"
    );
    // The default record: code ".", name "unknown ribonucleotide", empty formula.
    assert_eq!(
        Ribonucleotide::default().source_stream_text(),
        "Ribonucleotide '.' (unknown ribonucleotide, )"
    );
}
