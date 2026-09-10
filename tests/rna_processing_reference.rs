// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Pinned literal references and independent scientific checks. No C++ execution.
// See data/rna_processing_provenance.json for source/derived distinctions.

use openms::chemistry::{
    EmpiricalFormula, ModifiedNASequenceGenerator, NAFragmentType, NASequence,
    NucleicAcidSpectrumGenerator, RNaseDB, RNaseDigestion, Ribonucleotide, RibonucleotideDB,
    RibonucleotideRecord, RibonucleotideTermSpecificity, rnase::matches_code_pattern,
};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

fn sequence(text: &str) -> NASequence {
    NASequence::parse(text).unwrap()
}
fn rows(text: &str) -> impl Iterator<Item = Vec<&str>> {
    text.lines().skip(1).map(|line| line.split('\t').collect())
}
fn codes(sequences: &[NASequence]) -> Vec<String> {
    sequences.iter().map(ToString::to_string).collect()
}
fn modification(
    code: &str,
    origin: char,
    term: RibonucleotideTermSpecificity,
) -> Arc<Ribonucleotide> {
    Arc::new(
        Ribonucleotide::from_record(RibonucleotideRecord {
            code: code.into(),
            origin,
            term_specificity: term,
            ..RibonucleotideRecord::default()
        })
        .unwrap(),
    )
}

#[test]
fn twelve_literal_rnase_cases_and_thirty_eight_ordered_products() {
    let expected: Vec<_> =
        rows(include_str!("data/rna_processing_digestion_outputs.tsv")).collect();
    let mut count = 0;
    let mut products = 0;
    for row in rows(include_str!("data/rna_processing_digestion_cases.tsv")) {
        let mut digestion = RNaseDigestion::new(row[1]).unwrap();
        digestion.missed_cleavages = row[3].parse().unwrap();
        digestion.min_length = row[4].parse().unwrap();
        digestion.max_length = row[5].parse().unwrap();
        let output = digestion.digest(&sequence(row[2])).unwrap();
        let literal: Vec<_> = expected
            .iter()
            .filter(|r| r[0] == row[0])
            .map(|r| r[2])
            .collect();
        assert_eq!(output.len(), row[7].parse::<usize>().unwrap());
        assert_eq!(codes(&output), literal, "source line {}", row[6]);
        count += 1;
        products += output.len();
    }
    assert_eq!((count, products), (12, 38));
}

#[test]
fn all_registered_enzyme_fields_and_6048_independent_code_predicates() {
    let registry = RNaseDB::global();
    assert_eq!(registry.enzymes().len(), 14);
    let mut enzymes = 0;
    for r in rows(include_str!("data/rna_processing_enzymes.tsv")) {
        let enzyme = registry.get_enzyme(r[0]).unwrap();
        assert_eq!(enzyme.name(), r[0]);
        assert_eq!(enzyme.regex_description(), r[1]);
        assert_eq!(enzyme.cuts_after(), r[2]);
        assert_eq!(enzyme.cuts_before(), r[3]);
        assert_eq!(enzyme.three_prime_gain(), r[4]);
        assert_eq!(enzyme.five_prime_gain(), r[5]);
        enzymes += 1;
    }
    let mut count = 0;
    for r in rows(include_str!("data/rna_processing_code_predicates.tsv")) {
        assert_eq!(
            matches_code_pattern(r[0], r[2]).unwrap(),
            r[3] == "1",
            "{r:?}"
        );
        count += 1;
    }
    assert_eq!((enzymes, count), (14, 6048));
    // Search may succeed at a later occurrence after an earlier forbidden one.
    for (pattern, code, expected) in [
        ("G(?!m)", "GmG", true),
        ("G(?!m)", "Gm", false),
        ("G|A(?!m)", "Gm", true),
        ("G|Q(?!m)", "QmQ", true),
        ("(?<!m6)A", "m6AmA", true),
        ("(?<!m6)A", "m6A", false),
        ("(?<!m5)C", "m5CC", true),
        (".*(?!m)$", "Um", true),
        ("U|P|\\]|D|5", "]", true),
        ("U|P|m1Y|D|5", "]", false),
        ("^[^C]+$", "", false),
        ("^[^C]+$", "Cm", false),
    ] {
        assert_eq!(matches_code_pattern(pattern, code).unwrap(), expected);
    }
}

#[test]
fn rnase_coordinates_outer_ends_and_sulfur_gain_override() {
    let mut enzyme = RNaseDigestion::new("RNase_T1").unwrap();
    enzyme.missed_cleavages = 2;
    let output = enzyme
        .digest_with_positions(&sequence("pAUGUCGCAG"))
        .unwrap();
    assert_eq!(
        output
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
    // Slicing after G* supplies 5'-p*, then the enzyme's null five-prime gain
    // overwrites it. The original outer phosphate ends remain attached.
    let input = sequence("p[G*]Ap");
    assert!(input.suffix(1).unwrap().five_prime_mod().is_some());
    let output = RNaseDigestion::new("RNase_T1")
        .unwrap()
        .digest(&input)
        .unwrap();
    assert_eq!(codes(&output), ["p[G*]p", "Ap"]);
    assert!(output[1].five_prime_mod().is_none());
    let output = RNaseDigestion::new("RNase_H")
        .unwrap()
        .digest(&sequence("ACG"))
        .unwrap();
    assert_eq!(codes(&output), ["A", "pC", "pG"]);
    // mazF tests three separate record codes; it never considers boundary 0.
    let output = RNaseDigestion::new("mazF")
        .unwrap()
        .digest(&sequence("ACAACA"))
        .unwrap();
    assert_eq!(codes(&output), ["ACAp", "ACA"]);
}

fn cartesian_strings(choices: &[&[&str]]) -> BTreeSet<String> {
    // Independent full product: no generation limits, compatible-site masks or
    // production traversal. This also verifies identities, not just counts.
    choices
        .iter()
        .fold(BTreeSet::from([String::new()]), |prefixes, options| {
            prefixes
                .into_iter()
                .flat_map(|prefix| options.iter().map(move |part| format!("{prefix}{part}")))
                .collect()
        })
}

#[test]
fn literal_fixed_and_variable_modifications_with_complete_cartesian_identities() {
    let db = RibonucleotideDB::global();
    let generator = ModifiedNASequenceGenerator::default();
    let input = sequence("AUAUAUA");
    let mut fixed = input.clone();
    generator
        .apply_fixed_modifications(&[db.get("s4U").unwrap()], &mut fixed)
        .unwrap();
    let literals: Vec<_> =
        rows(include_str!("data/rna_processing_modification_outputs.tsv")).collect();
    assert_eq!(fixed, sequence(literals[0][2]));
    let mut count = 0;
    for r in rows(include_str!("data/rna_processing_modification_counts.tsv")) {
        let mods: Vec<_> = r[0].split('|').map(|name| db.get(name).unwrap()).collect();
        let output = generator
            .variable_modifications(&mods, &input, r[1].parse().unwrap(), r[2] == "true")
            .unwrap();
        assert_eq!(
            output.len(),
            r[4].parse::<usize>().unwrap(),
            "source line {}",
            r[5]
        );
        if r[1] == "1" {
            let actual: BTreeSet<_> = codes(&output).into_iter().collect();
            let expected: BTreeSet<_> = literals
                .iter()
                .filter(|s| s[0] == "single" && (r[2] == "true" || s[1] != "0"))
                .map(|s| s[2].to_string())
                .collect();
            assert_eq!(actual, expected);
        } else {
            let a: &[&str] = if r[1] == "7" { &["A", "[m1A]"] } else { &["A"] };
            let u: &[&str] = &["U", "[m3U]", "[s4U]"];
            assert_eq!(
                codes(&output).into_iter().collect::<BTreeSet<_>>(),
                cartesian_strings(&[a, u, a, u, a, u, a])
            );
        }
        if r[2] == "true" {
            assert_eq!(output[0], input);
        }
        count += 1;
    }
    assert_eq!(count, 4);
}

fn typed_state(seq: &NASequence) -> String {
    format!(
        "{}|{}|{}",
        seq.three_prime_mod().map_or("-", |r| r.code()),
        seq.five_prime_mod().map_or("-", |r| r.code()),
        seq.residues()
            .iter()
            .map(|r| r.code())
            .collect::<Vec<_>>()
            .join(",")
    )
}

#[test]
fn modified_rna_terminal_subset_order_and_nested_alternatives() {
    use RibonucleotideTermSpecificity::{Anywhere, FivePrime, ThreePrime};
    let generator = ModifiedNASequenceGenerator::default();
    let mods = [
        modification("T", 'U', ThreePrime),
        modification("F", 'U', FivePrime),
        modification("R1", 'A', Anywhere),
        modification("R2", 'A', Anywhere),
    ];
    let input = sequence("A");
    let output = generator
        .variable_modifications(&mods, &input, 3, true)
        .unwrap();
    // Ascending source site map T,F,R. One-site order R,F,T; two-site order
    // FR,TR,TF. Selected-site alternatives nest in ascending site order.
    let expected = [
        "-|-|A", "-|-|R1", "-|-|R2", "-|F|A", "T|-|A", "-|F|R1", "-|F|R2", "T|-|R1", "T|-|R2",
        "T|F|A", "T|F|R1", "T|F|R2",
    ];
    assert_eq!(output.iter().map(typed_state).collect::<Vec<_>>(), expected);
    let f2 = modification("F2", 'U', FivePrime);
    let mods = [mods[1].clone(), f2, mods[2].clone(), mods[3].clone()];
    let output = generator
        .variable_modifications(&mods, &input, 2, false)
        .unwrap();
    assert_eq!(
        output.iter().map(typed_state).collect::<Vec<_>>(),
        [
            "-|-|R1", "-|-|R2", "-|F|A", "-|F2|A", "-|F|R1", "-|F|R2", "-|F2|R1", "-|F2|R2"
        ]
    );
    // Fixed terminals choose first, while all candidates matching the original
    // unmodified residue remain eligible and the last ordinary candidate wins.
    let mut fixed = input;
    generator
        .apply_fixed_modifications(&mods, &mut fixed)
        .unwrap();
    assert_eq!(typed_state(&fixed), "-|F|R2");
}

#[test]
fn modified_rna_fast_path_and_allocation_identity_are_observable() {
    use RibonucleotideTermSpecificity::{Anywhere, FivePrime};
    let generator = ModifiedNASequenceGenerator::default();
    let terminal = modification("terminal-A", 'A', FivePrime);
    let input = sequence("A");
    let one = generator
        .variable_modifications(std::slice::from_ref(&terminal), &input, 1, false)
        .unwrap();
    assert_eq!(typed_state(&one[0]), "-|-|terminal-A");
    let two = generator
        .variable_modifications(std::slice::from_ref(&terminal), &input, 2, false)
        .unwrap();
    assert_eq!(typed_state(&two[0]), "-|terminal-A|A");
    let db = RibonucleotideDB::from_records(vec![
        (*terminal).clone(),
        (*RibonucleotideDB::global().get("A").unwrap()).clone(),
    ])
    .unwrap();
    assert!(one[0].checked_string_with_registry(&db).is_err());
    assert!(
        generator
            .variable_modifications(
                std::slice::from_ref(&terminal),
                &NASequence::default(),
                1,
                false
            )
            .unwrap()
            .is_empty()
    );
    let empty = generator
        .variable_modifications(&[terminal], &NASequence::default(), 2, false)
        .unwrap();
    assert_eq!(empty.len(), 1);
    assert!(empty[0].is_empty());
    assert!(empty[0].five_prime_mod().is_some());
    let first = modification("same", 'A', Anywhere);
    let equal_copy = Arc::new((*first).clone());
    let mods = [first.clone(), first.clone(), equal_copy.clone()];
    let output = generator
        .variable_modifications(&mods, &input, 1, false)
        .unwrap();
    assert_eq!(output.len(), 2); // repeated handle removed; equal allocation retained
    assert!(Arc::ptr_eq(&output[0].residues()[0], &first));
    assert!(Arc::ptr_eq(&output[1].residues()[0], &equal_copy));
    assert_eq!(output[0], output[1]);
}

fn single_series(series: &str) -> NucleicAcidSpectrumGenerator {
    NucleicAcidSpectrumGenerator {
        add_metainfo: true,
        add_first_prefix_ion: true,
        add_a_ions: series == "a",
        add_b_ions: series == "b",
        add_c_ions: series == "c",
        add_d_ions: series == "d",
        add_w_ions: series == "w",
        add_x_ions: series == "x",
        add_y_ions: series == "y",
        add_z_ions: series == "z",
        add_a_minus_b_ions: series == "a-B",
        ..NucleicAcidSpectrumGenerator::default()
    }
}

#[test]
fn all_132_literal_ion_encodings_and_126_source_compared_values() {
    let mut groups: BTreeMap<(&str, &str), Vec<Vec<&str>>> = BTreeMap::new();
    let mut encoded = 0;
    for r in rows(include_str!("data/rna_processing_source_ions.tsv")) {
        let literal: f64 = r[5].parse().unwrap();
        assert_eq!(literal.to_bits(), u64::from_str_radix(r[6], 16).unwrap());
        groups.entry((r[1], r[2])).or_default().push(r);
        encoded += 1;
    }
    assert_eq!((encoded, groups.len()), (132, 18));
    let mut compared = 0;
    for ((input, series), references) in groups {
        let output = single_series(series)
            .generate(&sequence(input), -1, -1)
            .unwrap();
        assert_eq!(output.peaks.len(), 7);
        assert_eq!(output.integer_data_arrays[0].data, [-1; 7]);
        let expected: Vec<_> = references.iter().filter(|r| r[8] == "true").collect();
        assert_eq!(expected.len(), output.peaks.len());
        for (index, (peak, reference)) in output.peaks.iter().zip(expected).enumerate() {
            let expected: f64 = reference[5].parse().unwrap();
            // Rounded Ariadne values embedded in source tests; source macro
            // tolerance is unavailable. This explicit 0.001 Da native check is
            // not a claim of a C++ runtime value or exact macro tolerance.
            assert!(
                (peak.mz - expected).abs() <= 0.001,
                "{input} {series}{} actual {}, literal {} at source line {}",
                index + 1,
                peak.mz,
                expected,
                reference[7]
            );
            assert_eq!(peak.intensity, 1.0);
            let label = if series == "a-B" {
                format!("a{}-B", index + 1)
            } else {
                format!("{series}{}", index + 1)
            };
            assert_eq!(output.string_data_arrays[0].data[index], label);
            compared += 1;
        }
    }
    assert_eq!(compared, 126);
}

#[test]
fn actual_default_multiple_spectrum_source_case() {
    // Source constructs but never applies its all-series/metadata Param.
    let generator = NucleicAcidSpectrumGenerator::default();
    let input = sequence("[m1A]UCCACAGp");
    let charges = BTreeSet::from([-1, -3, -5]);
    let multiple = generator.generate_multiple(&input, &charges, -1).unwrap();
    assert_eq!(multiple.len(), 3);
    for charge in charges {
        let direct = generator.generate(&input, -1, charge).unwrap();
        assert_eq!(multiple[&charge], direct);
        assert_eq!(direct.peaks.len(), 13 * charge.unsigned_abs() as usize);
        assert!(direct.string_data_arrays.is_empty());
        assert!(direct.integer_data_arrays.is_empty());
    }
}

#[test]
fn independent_ambiguous_base_loss_and_precursor_linkage_conventions() {
    let generator = single_series("a-B");
    let output = generator.generate(&sequence("[mA?]C"), -1, -1).unwrap();
    assert_eq!(output.peaks.len(), 2);
    assert_eq!(output.string_data_arrays[0].data, ["a1-B", "a1-B"]);
    assert!(output.peaks.iter().all(|p| p.intensity == 0.5));
    let methyl = EmpiricalFormula::parse("CH2").unwrap().mono_mass();
    assert!((output.peaks[1].mz - output.peaks[0].mz - methyl).abs() < 1e-12);
    let generator = NucleicAcidSpectrumGenerator {
        add_metainfo: true,
        add_precursor_peaks: true,
        ..NucleicAcidSpectrumGenerator::default()
    };
    let plain = generator.generate(&sequence("ACG"), -1, -1).unwrap();
    let sulfur = generator.generate(&sequence("[A*]CG"), -1, -1).unwrap();
    let mass_of_m = |s: &openms::kernel::MSSpectrum| {
        let i = s.string_data_arrays[0]
            .data
            .iter()
            .position(|n| n == "M")
            .unwrap();
        s.peaks[i].mz
    };
    // The optimized both-fragment precursor omits the first linkage sulfur.
    assert_eq!(mass_of_m(&plain).to_bits(), mass_of_m(&sulfur).to_bits());
    let formula_only = NucleicAcidSpectrumGenerator {
        add_b_ions: false,
        add_y_ions: false,
        ..generator
    };
    let sulfur = formula_only.generate(&sequence("[A*]CG"), -1, -1).unwrap();
    let plain = formula_only.generate(&sequence("ACG"), -1, -1).unwrap();
    let shift = EmpiricalFormula::parse("SO-1").unwrap().mono_mass();
    assert!((mass_of_m(&sulfur) - mass_of_m(&plain) - shift).abs() < 1e-10);
    // Published neutral-mass guards are formula-based and explicitly rounded.
    for (text, literal) in [
        ("[m1A]UCCACAGp", 2585.3800),
        ("[m1A]UCCACA[G*]p", 2585.3800),
        ("[m1A]UC[C*]AC[A*]Gp", 2617.334342),
    ] {
        assert!(
            (sequence(text).mono_mass(NAFragmentType::Full, 0).unwrap() - literal).abs() < 0.01
        );
    }
}
