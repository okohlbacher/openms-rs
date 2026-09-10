// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::rnase::matches_code_pattern;
use openms::chemistry::{
    DigestionEnzymeRNA, DigestionEnzymeRNARecord, EmpiricalFormula, NASequence, RNaseDB,
    RNaseDigestion, Ribonucleotide, RibonucleotideDB, RibonucleotideRecord,
};
use std::collections::{BTreeSet, HashSet};
use std::sync::Arc;

fn sequence(text: &str) -> NASequence {
    NASequence::parse(text).unwrap()
}
fn strings(values: &[NASequence]) -> Vec<String> {
    values.iter().map(ToString::to_string).collect()
}
fn enzyme(name: &str, after: &str, before: &str, five: &str, three: &str) -> DigestionEnzymeRNA {
    DigestionEnzymeRNA::from_record(DigestionEnzymeRNARecord {
        name: name.into(),
        cuts_after: after.into(),
        cuts_before: before.into(),
        five_prime_gain: five.into(),
        three_prime_gain: three.into(),
        ..Default::default()
    })
    .unwrap()
}
fn record(code: &str, formula: &str) -> Ribonucleotide {
    Ribonucleotide::from_record(RibonucleotideRecord {
        name: code.into(),
        code: code.into(),
        formula: EmpiricalFormula::parse(formula).unwrap(),
        ..Default::default()
    })
    .unwrap()
}

#[test]
fn registry_defaults_fields_aliases_and_complete_value_identity() {
    let registry = RNaseDB::global();
    assert!(std::ptr::eq(registry, RNaseDB::global()));
    assert_eq!(registry.enzymes().len(), 14);
    assert_eq!(registry.names().len(), 14);
    for name in registry.names() {
        assert!(registry.has_enzyme(name));
    }
    let t1 = registry.get_enzyme("RNase_T1").unwrap();
    assert!(Arc::ptr_eq(&t1, &registry.get_enzyme("rnase_t1").unwrap()));
    assert!(!registry.has_enzyme("RNASE_T1"));
    assert!(registry.get_enzyme("missing").is_err());
    assert!(!registry.has_regex(""));
    assert!(registry.enzyme_by_regex("G(?!m)").is_err());
    assert_eq!(t1.cuts_after(), "G(?!m)");
    assert_eq!(t1.cuts_before(), ".*");
    assert_eq!(t1.three_prime_gain(), "p");
    assert_eq!(DigestionEnzymeRNA::default().name(), "unknown_enzyme");
    let mut changed = t1.record().clone();
    changed.three_prime_gain = "c".into();
    let changed = DigestionEnzymeRNA::from_record(changed).unwrap();
    assert_ne!(&changed, t1.as_ref());
    assert_eq!(HashSet::from([changed, t1.as_ref().clone()]).len(), 2);
    let huge = DigestionEnzymeRNARecord {
        name: "x".repeat(65_537),
        ..Default::default()
    };
    assert!(DigestionEnzymeRNA::from_record(huge).is_err());
}

#[test]
fn later_provider_replacement_preserves_literal_source_index_erasure() {
    let mut first = enzyme("Example", "G", "", "", "").record().clone();
    first.regex = "A".into();
    first.synonyms = BTreeSet::from(["shared".into()]);
    let mut second = first.clone();
    second.name = "other".into();
    let mut replacement = first.clone();
    replacement.regex = "U".into();
    replacement.synonyms = BTreeSet::from(["fresh".into()]);
    let registry = RNaseDB::from_records(vec![
        DigestionEnzymeRNA::from_record(first).unwrap(),
        DigestionEnzymeRNA::from_record(second).unwrap(),
        DigestionEnzymeRNA::from_record(replacement).unwrap(),
    ])
    .unwrap();
    assert_eq!(registry.names(), ["other", "Example"]);
    assert!(registry.has_enzyme("other"));
    assert!(registry.has_enzyme("example"));
    assert!(registry.has_enzyme("fresh"));
    // Removing the replaced record erases shared indices even if a different
    // retained record most recently supplied that alias or regex (source quirk).
    assert!(!registry.has_enzyme("shared"));
    assert!(!registry.has_regex("A"));
    let owned = registry.enzyme_by_regex("U").unwrap();
    drop(registry);
    assert_eq!(owned.name(), "Example");
}

#[test]
fn raw_code_predicates_preserve_regex_search_and_lookaround_precedence() {
    for (pattern, code, expected) in [
        ("G(?!m)", "m1G", true),
        ("G(?!m)", "Gm", false),
        ("G(?!m)", "GmG", true),
        ("G|A(?!m)", "Gm", true),
        ("G|A(?!m)", "Am", false),
        ("G|Q(?!m)", "Qm", false),
        ("G|Q(?!m)", "Gm", true),
        (".*(?!m)$", "Um", true),
        ("^[^C]+$", "m5C", false),
        ("^[^C]+$", "", false),
        ("^[^C]+$", "ä\nU", true),
        ("(?<!m6)A", "m6A", false),
        ("(?<!m6)A", "m6AA", true),
        ("(?<!m5)C", "m5C", false),
        ("U|P|\\]|D|5", "m1Y", false),
        ("U|P|m1Y|D|5", "m1Y", true),
        ("U|P|\\]|D|5", "a]", true),
        ("A", "m6A", true),
        ("", "", true),
        ("m6A", "[m6A]", true),
    ] {
        assert_eq!(
            matches_code_pattern(pattern, code).unwrap(),
            expected,
            "{pattern} {code}"
        );
    }
    assert!(matches_code_pattern("^A$", "A").is_err());
    assert!(matches_code_pattern("A", &"A".repeat(65_537)).is_err());
}

#[test]
fn positions_use_source_start_then_missed_order_and_half_open_coordinates() {
    let mut digest = RNaseDigestion::default();
    digest.missed_cleavages = usize::MAX;
    let products = digest
        .digest_with_positions(&sequence("pAUGUCGCAG"))
        .unwrap();
    assert_eq!(
        products
            .iter()
            .map(|p| (p.start, p.end, p.missed_cleavages))
            .collect::<Vec<_>>(),
        [
            (0, 3, 0),
            (0, 6, 1),
            (0, 9, 2),
            (3, 6, 0),
            (3, 9, 1),
            (6, 9, 0)
        ]
    );
    assert_eq!(
        products
            .iter()
            .map(|p| p.sequence.to_string())
            .collect::<Vec<_>>(),
        ["pAUGp", "pAUGUCGp", "pAUGUCGCAG", "UCGp", "UCGCAG", "CAG"]
    );
}

#[test]
fn length_bounds_empty_input_and_unspecific_order_are_checked_before_subtraction() {
    let mut digest = RNaseDigestion::new("unspecific cleavage").unwrap();
    assert!(digest.digest(&sequence("p")).unwrap().is_empty());
    let rna = sequence("ACG");
    digest.min_length = 4;
    assert!(digest.digest(&rna).unwrap().is_empty());
    digest.min_length = 3;
    digest.max_length = 2;
    assert!(digest.digest(&rna).unwrap().is_empty());
    digest.min_length = 2;
    digest.max_length = usize::MAX;
    digest.missed_cleavages = usize::MAX;
    assert_eq!(strings(&digest.digest(&rna).unwrap()), ["ACp", "ACG", "CG"]);
    digest.set_enzyme("no cleavage").unwrap();
    assert_eq!(strings(&digest.digest(&rna).unwrap()), ["ACG"]);
    digest.max_length = 2;
    assert!(digest.digest(&rna).unwrap().is_empty());
    let mut output = vec![rna];
    digest.digest_into(&sequence(""), &mut output).unwrap();
    assert!(output.is_empty());
}

#[test]
fn enzyme_gains_replace_internal_sulfur_ends_and_preserve_original_boundaries() {
    let digest = RNaseDigestion::default();
    let products = digest.digest(&sequence("pA[G*]Up")).unwrap();
    assert_eq!(strings(&products), ["pA[G*]p", "Up"]);
    assert!(products[1].five_prime_mod().is_none());
    let digest = RNaseDigestion::new("RNase_H").unwrap();
    assert_eq!(
        strings(&digest.digest(&sequence("pACGp")).unwrap()),
        ["pA", "pC", "pGp"]
    );
    let mut all = digest.clone();
    all.missed_cleavages = 1;
    assert_eq!(
        strings(&all.digest(&sequence("ACG")).unwrap()),
        ["A", "AC", "pC", "pCG", "pG"]
    );
}

#[test]
fn source_multi_position_patterns_match_raw_records_in_each_direction() {
    let custom = Arc::new(enzyme("compound", "A,G", "U", "", "p"));
    let digest = RNaseDigestion::with_enzyme(custom, RibonucleotideDB::global()).unwrap();
    assert_eq!(
        strings(&digest.digest(&sequence("CA[m1G]UAC")).unwrap()),
        ["CA[m1G]p", "UAC"]
    );
    let custom = Arc::new(enzyme("unconstrained", "", "", "", ""));
    let digest = RNaseDigestion::with_enzyme(custom, RibonucleotideDB::global()).unwrap();
    assert_eq!(
        strings(&digest.digest(&sequence("ACG")).unwrap()),
        ["A", "C", "G"]
    );
    let custom = Arc::new(enzyme("two upstream", ",", "", "", ""));
    let digest = RNaseDigestion::with_enzyme(custom, RibonucleotideDB::global()).unwrap();
    assert_eq!(
        strings(&digest.digest(&sequence("ACG")).unwrap()),
        ["AC", "G"]
    );
    let digest = RNaseDigestion::new("mazF").unwrap();
    assert_eq!(
        strings(&digest.digest(&sequence("ACAA")).unwrap()),
        ["ACAA"]
    );
}

#[test]
fn custom_gains_use_source_code_rewriting_and_keep_owned_chemistry() {
    let registry =
        RibonucleotideDB::from_records(vec![record("fp", "H2"), record("[cap]", "H2O")]).unwrap();
    let enzyme = Arc::new(enzyme("custom", "C", "", "fp", "cap"));
    let digest = RNaseDigestion::with_enzyme(enzyme.clone(), &registry).unwrap();
    drop(registry);
    let result = digest.digest(&sequence("ACA")).unwrap();
    assert_eq!(result[0].three_prime_mod().unwrap().code(), "[cap]");
    assert_eq!(result[1].five_prime_mod().unwrap().code(), "fp");
    assert_eq!(
        result[0].three_prime_mod().unwrap().formula(),
        &EmpiricalFormula::parse("H2O").unwrap()
    );
    // Merely registering cap does not satisfy source lookup of [cap].
    let wrong =
        RibonucleotideDB::from_records(vec![record("fp", "H2"), record("cap", "O")]).unwrap();
    assert!(RNaseDigestion::with_enzyme(enzyme, &wrong).is_err());
}

#[test]
fn failed_configuration_preserves_enzyme_gains_and_existing_options() {
    let mut digest = RNaseDigestion::default();
    digest.missed_cleavages = 1;
    let before = digest.digest(&sequence("AGUC")).unwrap();
    assert!(digest.set_enzyme("not_an_enzyme").is_err());
    let invalid = Arc::new(enzyme("unsupported", "(?=C)", "", "", ""));
    assert!(
        digest
            .set_enzyme_with_registry(invalid, RibonucleotideDB::global())
            .is_err()
    );
    assert_eq!(digest.enzyme().name(), "RNase_T1");
    assert_eq!(digest.missed_cleavages, 1);
    assert_eq!(digest.digest(&sequence("AGUC")).unwrap(), before);
}

#[test]
fn work_count_bytes_and_missing_sulfur_context_fail_atomically() {
    let rna = sequence("ACGU");
    let original = vec![sequence("pAGp")];
    for (products, work, bytes, residues) in [
        (1, 50_000_000, 256 * 1024 * 1024, 100),
        (100, 0, 256 * 1024 * 1024, 100),
        (100, 50_000_000, 200, 100),
        (100, 50_000_000, 256 * 1024 * 1024, 1),
    ] {
        let mut digest = RNaseDigestion::new("unspecific cleavage").unwrap();
        digest.max_products = products;
        digest.max_work = work;
        digest.max_output_bytes = bytes;
        digest.max_residues = residues;
        let mut output = original.clone();
        assert!(digest.digest_into(&rna, &mut output).is_err());
        assert_eq!(output, original);
    }
    let empty_registry = RibonucleotideDB::default();
    let rna = NASequence::from_records_with_registry(
        vec![
            RibonucleotideDB::global().get("G*").unwrap(),
            RibonucleotideDB::global().get("A").unwrap(),
        ],
        &empty_registry,
    )
    .unwrap();
    let mut output = original.clone();
    assert!(
        RNaseDigestion::default()
            .digest_into(&rna, &mut output)
            .is_err()
    );
    assert_eq!(output, original);
}

#[test]
fn many_short_fragments_are_budgeted_by_their_own_residues() {
    let input = sequence(&"G".repeat(1_000));
    let mut digest = RNaseDigestion::default();
    digest.max_output_bytes = 4_000_000;
    digest.max_length = 1;
    digest.missed_cleavages = usize::MAX;
    let output = digest.digest_with_positions(&input).unwrap();
    assert_eq!(output.len(), 1_000);
    assert_eq!(output[999].sequence.to_string(), "G");
    assert_eq!(output[500].start, 500);
}

#[test]
fn long_common_prefix_enzyme_keys_consume_the_registry_work_budget() {
    let prefix = "A".repeat(1_024);
    let records = (0..2_048)
        .map(|i| enzyme(&format!("{prefix}{i:04}"), "G", "", "", ""))
        .collect();
    let error = RNaseDB::from_records(records).unwrap_err();
    assert!(error.to_string().contains("work limit"));
}
