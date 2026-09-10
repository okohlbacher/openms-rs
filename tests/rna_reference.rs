// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Literal C++ class-test fixtures and independent source-derived RNA checks.
// See data/rna_provenance.json; no C++ execution generated these references.

use openms::chemistry::{
    EmpiricalFormula, NAFragmentType, NASequence, Ribonucleotide, RibonucleotideDB,
    RibonucleotideRecord, RibonucleotideTermSpecificity,
};
use std::collections::{BTreeSet, HashSet};
use std::sync::Arc;

fn fragment(name: &str) -> NAFragmentType {
    match name {
        "Full" => NAFragmentType::Full,
        "AIon" => NAFragmentType::AIon,
        "BIon" => NAFragmentType::BIon,
        "CIon" => NAFragmentType::CIon,
        "DIon" => NAFragmentType::DIon,
        "WIon" => NAFragmentType::WIon,
        "XIon" => NAFragmentType::XIon,
        "YIon" => NAFragmentType::YIon,
        "ZIon" => NAFragmentType::ZIon,
        "AminusB" => NAFragmentType::AminusB,
        _ => panic!("unknown source fragment {name}"),
    }
}
fn formula(text: &str) -> EmpiricalFormula {
    EmpiricalFormula::parse(text).unwrap()
}
fn sequence(text: &str) -> NASequence {
    NASequence::parse(text).unwrap()
}

#[test]
fn all_fifteen_literal_source_formula_assertions() {
    let mut count = 0;
    for line in include_str!("data/rna_source_formulas.tsv").lines().skip(1) {
        let fields: Vec<_> = line.split('\t').collect();
        let expected = fields[3]
            .split('|')
            .try_fold(EmpiricalFormula::default(), |sum, term| {
                sum.checked_add(&formula(term))
            })
            .unwrap();
        assert_eq!(
            sequence(fields[0])
                .formula(fragment(fields[1]), fields[2].parse().unwrap())
                .unwrap(),
            expected,
            "source line {}: {line}",
            fields[4]
        );
        count += 1;
    }
    assert_eq!(count, 15);
}

#[test]
fn all_twenty_three_literal_source_mono_and_average_assertions() {
    let mut counts = [0, 0];
    for line in include_str!("data/rna_source_masses.tsv").lines().skip(1) {
        let fields: Vec<_> = line.split('\t').collect();
        let peptide = sequence(fields[1]);
        let kind = fragment(fields[2]);
        let charge = fields[3].parse().unwrap();
        let average = fields[0] == "average";
        let actual = if average {
            peptide.average_mass(kind, charge).unwrap()
        } else {
            peptide.mono_mass(kind, charge).unwrap()
        } / fields[4].parse::<f64>().unwrap();
        let expected: f64 = fields[6].parse().unwrap();
        assert_eq!(
            expected.to_bits(),
            u64::from_str_radix(fields[7], 16).unwrap()
        );
        // These are rounded external-tool values reproduced in C++ tests, not
        // high-precision implementation outputs. The 1201.3 average reference
        // has only one decimal; the GGG charge -2 reference is also approximate.
        // This explicit native tolerance does not claim a missing ClassTest
        // macro's exact settings. Exact elemental references are tested above.
        let tolerance = if fields[5] == "1201.3" {
            0.05
        } else {
            expected.abs() * 1e-5
        };
        assert!(
            (actual - expected).abs() <= tolerance,
            "source line {}: actual {actual}, reference {expected}, tolerance {tolerance}",
            fields[8]
        );
        counts[usize::from(average)] += 1;
    }
    assert_eq!(counts, [18, 5]);
}

#[test]
fn all_ten_source_slices_and_full_length_rejections() {
    let mut count = 0;
    for line in include_str!("data/rna_source_slices.tsv").lines().skip(1) {
        let fields: Vec<_> = line.split('\t').collect();
        let original = sequence(fields[0]);
        let argument = fields[2].parse().unwrap();
        let result = match fields[1] {
            "prefix" => original.prefix(argument),
            "suffix" => original.suffix(argument),
            "subsequence" => original.subsequence(
                argument,
                (fields[3] != "-").then(|| fields[3].parse().unwrap()),
            ),
            _ => panic!("unknown slice"),
        }
        .unwrap();
        assert_eq!(result.to_string(), fields[4], "source line {}", fields[5]);
        count += 1;
    }
    assert_eq!(count, 10);
    let original = sequence("pACp");
    for length in [2, 3, usize::MAX] {
        assert!(original.prefix(length).is_err());
        assert!(original.suffix(length).is_err());
    }
    assert!(original.subsequence(2, None).is_err());
    assert_eq!(original.subsequence(0, Some(usize::MAX)).unwrap(), original);
}

fn expected_record_mass(kind: &str, bits: &str, form: &EmpiricalFormula, average: bool) -> f64 {
    if kind == "formula" {
        if average {
            form.average_mass()
        } else {
            form.mono_mass()
        }
    } else {
        assert_eq!(kind, "declared");
        f64::from_bits(u64::from_str_radix(bits, 16).unwrap())
    }
}

#[test]
fn complete_independent_378_record_projection_and_last_code_wins() {
    let registry = RibonucleotideDB::global();
    let rows: Vec<_> = include_str!("data/rna_registry_reference.tsv")
        .lines()
        .skip(1)
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .collect();
    assert_eq!(rows.len(), 378);
    assert_eq!(registry.entries().len(), rows.len());
    for (record, row) in registry.entries().iter().zip(&rows) {
        assert_eq!(record.code(), row[2], "{} {}", row[0], row[1]);
        let expected_formula = formula(row[3]);
        assert_eq!(record.formula(), &expected_formula, "{}", row[2]);
        assert_eq!(record.origin(), row[4].chars().next().unwrap());
        let term = match row[5] {
            "Anywhere" => RibonucleotideTermSpecificity::Anywhere,
            "FivePrime" => RibonucleotideTermSpecificity::FivePrime,
            "ThreePrime" => RibonucleotideTermSpecificity::ThreePrime,
            _ => panic!("invalid fixture term"),
        };
        assert_eq!(record.term_specificity(), term);
        assert_eq!(record.baseloss_formula(), &formula(row[6]));
        assert_eq!(
            record.mono_mass().to_bits(),
            expected_record_mass(row[7], row[8], &expected_formula, false).to_bits(),
            "mono {} {}",
            row[0],
            row[1]
        );
        assert_eq!(
            record.average_mass().to_bits(),
            expected_record_mass(row[9], row[10], &expected_formula, true).to_bits(),
            "average {} {}",
            row[0],
            row[1]
        );
        if !row[11].is_empty() {
            let [left, right] = registry.alternatives(row[2]).unwrap();
            assert_eq!(left.code(), row[11]);
            assert_eq!(right.code(), row[12]);
        }
    }
    assert_eq!(
        rows.iter().map(|row| row[2]).collect::<BTreeSet<_>>().len(),
        375
    );
    assert_eq!(registry.entries()[0].code(), "io6A");
    for (code, expected) in [("pm1acp3Y", "C13H18N3O11P"), ("pm2,7Gm", "C13H18N5O8P")] {
        assert_eq!(registry.get(code).unwrap().formula(), &formula(expected));
    }
    assert_eq!(registry.get("Am").unwrap().name(), "2'-O-methyladenosine");
    assert_eq!(registry.get_prefix("m1AmCGU").unwrap().code(), "m1Am");
    assert!(registry.get("bla").is_err());
    assert!(registry.get_prefix("blam1A").is_err());
}

#[test]
fn real_all_record_canonicalization_and_literal_terminal_overwrite_semantics() {
    let registry = RibonucleotideDB::global();
    // NASequence_test.cpp's intended all-record loop cannot reach its assertion
    // (a code cannot equal both c and p). This loop actually visits every row.
    for record in registry.entries() {
        let parsed =
            NASequence::parse_with_registry(&format!("[{}]", record.code()), registry).unwrap();
        if matches!(record.code(), "3'-p" | "3'-c") {
            // Source Display emits a solitary p/c for an empty three-prime
            // end. That spelling is a five-prime end or invalid on reparsing.
            assert!(parsed.checked_string_with_registry(registry).is_err());
        } else {
            let rendered = parsed
                .checked_string_with_registry(registry)
                .unwrap_or_else(|error| panic!("{}: {error}", record.code()));
            assert_eq!(
                NASequence::parse_with_registry(&rendered, registry).unwrap(),
                parsed
            );
        }
        let contextual =
            NASequence::parse_with_registry(&format!("A[{}]A", record.code()), registry).unwrap();
        let rendered = contextual.checked_string_with_registry(registry).unwrap();
        assert_eq!(
            NASequence::parse_with_registry(&rendered, registry).unwrap(),
            contextual
        );
    }
    for (text, expected) in [
        ("A[5'-p]C[3'-p][5'-p*]", "*ACp"),
        ("[3'-p]AC[3'-c]", "ACc"),
        ("p A C p", "pACp"),
        ("[A][C]", "AC"),
        ("p", "p"),
        ("*p", "*p"),
        ("pp", "pp"),
    ] {
        assert_eq!(sequence(text).to_string(), expected, "{text}");
    }
    for text in ["c", "Am", " pAC", "ACp ", "A\tC", "[missing]", "A[Am"] {
        assert!(NASequence::parse(text).is_err(), "{text}");
    }
    let original = sequence("pA[C*]p");
    assert_eq!(original.prefix(0).unwrap().to_string(), "p");
    assert_eq!(original.suffix(0).unwrap().to_string(), "*p");
    assert_eq!(
        sequence("[A*]Cp")
            .subsequence(1, Some(0))
            .unwrap()
            .to_string(),
        "*"
    );
}

#[test]
fn source_empty_charge_unknown_composition_and_legacy_fragment_fallbacks() {
    use NAFragmentType::*;
    let electron = 1.0 / 1_822.888_502_047_7;
    for text in ["", "p", "pp", "*p", "[pN]"] {
        let empty = sequence(text);
        assert!(empty.is_empty());
        for charge in [-3, -1, 0, 1, 4] {
            assert_eq!(
                empty.formula(Full, charge).unwrap(),
                EmpiricalFormula::default()
            );
            assert_eq!(
                empty.mono_mass(Full, charge).unwrap(),
                -f64::from(charge) * electron
            );
            assert_eq!(
                empty.average_mass(AminusB, charge).unwrap(),
                -f64::from(charge) * electron
            );
        }
    }
    // Source N has explicitly empty formula and zero mass. Linkages still add
    // H-1PO2 per bond; this is a source computation, not an invented N formula.
    assert_eq!(
        sequence("NNN").formula(Full, 0).unwrap(),
        formula("H-2O4P2")
    );
    let base = formula("C19H25N8O11P"); // A + C + H-1PO2, no end groups.
    for kind in [
        Internal,
        FivePrime,
        ThreePrime,
        Precursor,
        BIonMinusH2O,
        YIonMinusH2O,
        BIonMinusNH3,
        YIonMinusNH3,
        NonIdentified,
        Unannotated,
    ] {
        assert_eq!(sequence("pACp").formula(kind, 7).unwrap(), base);
        assert_eq!(
            sequence("pACp").mono_mass(kind, 7).unwrap(),
            base.mono_mass() - 7.0 * electron
        );
    }
}

fn custom(code: &str, form: &str, mono_mass: f64) -> Ribonucleotide {
    Ribonucleotide::from_record(RibonucleotideRecord {
        code: code.into(),
        formula: formula(form),
        mono_mass,
        average_mass: -17.0,
        ..Default::default()
    })
    .unwrap()
}

#[test]
fn sequence_custom_chemistry_survives_registry_drop_and_uses_complete_identity() {
    let left_record = Arc::new(custom("custom", "C2H4", 999.0));
    let right_record = Arc::new(custom("custom", "C2H6", 999.0));
    let declared_difference = Arc::new(custom("custom", "C2H4", 1000.0));
    let left = NASequence::from_records(vec![left_record.clone()]).unwrap();
    let equal_copy = NASequence::from_records(vec![Arc::new((*left_record).clone())]).unwrap();
    let right = NASequence::from_records(vec![right_record]).unwrap();
    let different_declared = NASequence::from_records(vec![declared_difference]).unwrap();
    assert_eq!(left, equal_copy);
    assert_ne!(left, right);
    assert_ne!(left, different_declared);
    assert_eq!(left.to_string(), right.to_string());
    assert_eq!(
        left.mono_mass(NAFragmentType::Full, 0).unwrap(),
        formula("C2H4").mono_mass()
    );
    assert_eq!(
        left.average_mass(NAFragmentType::Full, 0).unwrap(),
        formula("C2H4").average_mass()
    );
    let hash_values: HashSet<_> = [
        left.clone(),
        equal_copy.clone(),
        right.clone(),
        different_declared.clone(),
    ]
    .into_iter()
    .collect();
    let ordered_values: BTreeSet<_> = [left.clone(), equal_copy, right.clone(), different_declared]
        .into_iter()
        .collect();
    assert_eq!(hash_values.len(), 3);
    assert_eq!(ordered_values.len(), 3);
    let parsed = {
        let registry = RibonucleotideDB::from_records(vec![(*left_record).clone()]).unwrap();
        let parsed = NASequence::parse_with_registry("[custom]", &registry).unwrap();
        assert_eq!(
            parsed.checked_string_with_registry(&registry).unwrap(),
            "[custom]"
        );
        assert!(right.checked_string_with_registry(&registry).is_err());
        parsed
    };
    assert_eq!(parsed, left);
    assert_eq!(
        parsed.formula(NAFragmentType::Full, 0).unwrap(),
        formula("C2H4")
    );
}

#[test]
fn original_provider_rows_equal_embedded_records_and_source_branch_precedence() {
    use openms::chemistry::ribonucleotide_db::read_tsv;
    let source_tsv = include_str!("../resources/rna/Custom_RNA_modifications.tsv");
    let report = read_tsv(source_tsv).unwrap();
    assert_eq!(report.skipped, 0);
    assert_eq!(report.entries.len(), 45);
    let registry = RibonucleotideDB::global();
    for (entry, embedded) in report.entries.iter().zip(&registry.entries()[333..]) {
        assert_eq!(entry.ribonucleotide.as_ref(), embedded.as_ref());
    }
    let header = "name\tshort_name\tnew_nomenclature\toriginating_base\trnamods_abbrev\thtml_abbrev\tformula\tmonoisotopic_mass\taverage_mass\talternatives\n";
    let rows = concat!(
        "zero\tz\t\tA\t\t\tH2\t0\t0\t\n",
        "missing\tv\t\tA\t\t\tH2\tNone\t\t\n",
        "terminal\tdAm\tN\tA\t\t\tH2\t0\t0\t\n",
        "exception\tdXm\tGN\tG\t\t\tH2\t0\t0\t\n",
        "deoxy\tdX?\t\tA\t\t\tH2\t0\t0\tz v\n",
        "ambiguous\tX?\t\tA\t\t\tH2\t0\t0\tz ignored v"
    );
    // Deliberately no newline on the last row: source getline still emits it.
    let report = read_tsv(&format!("{header}{rows}")).unwrap();
    assert_eq!(report.skipped, 0);
    assert_eq!(report.entries.len(), 6);
    let records: Vec<_> = report
        .entries
        .iter()
        .map(|entry| &entry.ribonucleotide)
        .collect();
    assert_eq!(records[0].mono_mass(), formula("H2").mono_mass());
    assert_eq!(records[0].average_mass(), formula("H2").average_mass());
    assert_eq!(records[1].mono_mass(), 0.0);
    assert_eq!(records[1].average_mass(), 0.0);
    assert_eq!(
        records[2].term_specificity(),
        RibonucleotideTermSpecificity::FivePrime
    );
    for record in &records[2..4] {
        assert_eq!(record.baseloss_formula(), &formula("C5H10O5"));
    }
    assert_eq!(
        records[3].term_specificity(),
        RibonucleotideTermSpecificity::Anywhere
    );
    assert_eq!(records[4].baseloss_formula(), &formula("C5H10O4"));
    assert!(
        !report.entries[4].is_ambiguous(),
        "deoxy branch precedes ambiguity"
    );
    assert_eq!(
        report.entries[5].alternatives.as_ref().unwrap(),
        &["z", "v"]
    );

    #[cfg(feature = "rna-json")]
    {
        use openms::chemistry::ribonucleotide_db::read_modomics_json;
        let report = read_modomics_json(include_str!("../resources/rna/Modomics.json")).unwrap();
        assert_eq!(report.entries.len(), 333);
        for (entry, embedded) in report.entries.iter().zip(&registry.entries()[..333]) {
            assert_eq!(
                entry.ribonucleotide.as_ref(),
                embedded.as_ref(),
                "{}",
                embedded.code()
            );
        }
        // JSON uses lexical object keys, keeps explicit zero mono mass, ignores
        // the content of four reference moieties, and distinguishes ?* at the
        // provider level from Ribonucleotide::is_ambiguous().
        let source = r#"{
            "2":{"name":"fallback","short_name":"tail","reference_moiety":["A"],"formula":"H2"},
            "10":{"name":"zero","short_name":"ApN","reference_moiety":["A"],"formula":"H2","mass_monoiso":0},
            "11":{"name":"wild","short_name":"pN","reference_moiety":[null,7,{},[]],"formula":"H2"},
            "12":{"name":"ambiguous","short_name":"X?*","reference_moiety":["A"],"formula":"H2","alternatives":["ApN","tail","unused"]}
        }"#;
        let report = read_modomics_json(source).unwrap();
        assert_eq!(report.skipped, 0);
        let entries = report.entries;
        assert_eq!(
            entries
                .iter()
                .map(|e| e.ribonucleotide.code())
                .collect::<Vec<_>>(),
            ["ApN", "pN", "X?*", "tail"]
        );
        assert_eq!(entries[0].ribonucleotide.mono_mass(), 0.0);
        assert_eq!(
            entries[0].ribonucleotide.term_specificity(),
            RibonucleotideTermSpecificity::Anywhere
        );
        assert_eq!(entries[1].ribonucleotide.origin(), 'X');
        assert_eq!(
            entries[1].ribonucleotide.term_specificity(),
            RibonucleotideTermSpecificity::FivePrime
        );
        assert!(entries[2].is_ambiguous());
        assert!(!entries[2].ribonucleotide.is_ambiguous());
        assert_eq!(entries[2].alternatives.as_ref().unwrap(), &["ApN", "tail"]);
        assert_eq!(
            entries[3].ribonucleotide.mono_mass(),
            formula("H2").mono_mass()
        );
        assert_eq!(entries[3].ribonucleotide.average_mass(), 0.0);
    }
}
