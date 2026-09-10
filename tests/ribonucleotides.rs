// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

#[cfg(feature = "rna-json")]
use openms::chemistry::ribonucleotide_db::read_modomics_json;
use openms::chemistry::ribonucleotide_db::read_tsv;
use openms::chemistry::{
    Ribonucleotide, RibonucleotideDB, RibonucleotideEntry, RibonucleotideRecord,
    RibonucleotideTermSpecificity,
};
use std::collections::{BTreeSet, HashSet};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::sync::Arc;

const HEADER: &str = "name\tshort_name\tnew_nomenclature\toriginating_base\trnamods_abbrev\thtml_abbrev\tformula\tmonoisotopic_mass\taverage_mass";
fn nucleotide(code: &str, formula: &str, mass: f64) -> Ribonucleotide {
    Ribonucleotide::from_record(RibonucleotideRecord {
        code: code.into(),
        formula: formula.parse().unwrap(),
        mono_mass: mass,
        ..Default::default()
    })
    .unwrap()
}
fn hash(record: &Ribonucleotide) -> u64 {
    let mut hasher = DefaultHasher::new();
    record.hash(&mut hasher);
    hasher.finish()
}
fn entry(code: &str, mass: f64, alternatives: Option<[&str; 2]>) -> RibonucleotideEntry {
    RibonucleotideEntry {
        ribonucleotide: Arc::new(nucleotide(code, "C", mass)),
        alternatives: alternatives.map(|[a, b]| [a.into(), b.into()]),
    }
}

#[test]
fn source_defaults_and_independent_declared_fields() {
    let default = Ribonucleotide::default();
    assert_eq!(default.name(), "unknown ribonucleotide");
    assert_eq!(default.code(), ".");
    assert_eq!(default.new_code(), "");
    assert_eq!(default.html_code(), ".");
    assert_eq!(default.origin(), '.');
    assert_eq!(default.mono_mass(), 0.0);
    assert_eq!(default.average_mass(), 0.0);
    assert!(default.formula().is_empty());
    assert_eq!(default.baseloss_formula(), &"C5H10O5".parse().unwrap());
    assert!(!default.is_modified());
    assert!(!default.is_ambiguous());
    assert_eq!(
        default.to_string(),
        "Ribonucleotide '.' (unknown ribonucleotide, )"
    );
    let record = Ribonucleotide::from_record(RibonucleotideRecord {
        code: "A".into(),
        origin: 'A',
        formula: "C10H13N5O4".parse().unwrap(),
        mono_mass: -12.5,
        average_mass: 500.0,
        ..Default::default()
    })
    .unwrap();
    assert_eq!(record.mono_mass(), -12.5);
    assert_eq!(record.average_mass(), 500.0);
    assert_ne!(record.formula().mono_mass(), record.mono_mass());
    assert!(!record.is_modified());
    assert!(nucleotide("mA?", "", 0.0).is_ambiguous());
    assert!(!nucleotide("mA?*", "", 0.0).is_ambiguous());
}

#[test]
fn every_stored_field_participates_in_value_identity() {
    let base = RibonucleotideRecord::default();
    let mut records = vec![Ribonucleotide::from_record(base.clone()).unwrap()];
    for field in 0..10 {
        let mut changed = base.clone();
        match field {
            0 => changed.name = "other".into(),
            1 => changed.code = "other".into(),
            2 => changed.new_code = "other".into(),
            3 => changed.html_code = "other".into(),
            4 => changed.formula = "C+".parse().unwrap(),
            5 => changed.origin = 'X',
            6 => changed.mono_mass = 1.0,
            7 => changed.average_mass = 1.0,
            8 => changed.term_specificity = RibonucleotideTermSpecificity::FivePrime,
            _ => changed.baseloss_formula = "C5H10O5+".parse().unwrap(),
        }
        records.push(Ribonucleotide::from_record(changed).unwrap());
    }
    assert_eq!(records.iter().collect::<BTreeSet<_>>().len(), 11);
    assert_eq!(records.iter().collect::<HashSet<_>>().len(), 11);
    for left in &records {
        for right in &records {
            assert_eq!(left == right, left.cmp(right).is_eq());
        }
    }
    let positive = Ribonucleotide::default();
    let negative = Ribonucleotide::from_record(RibonucleotideRecord {
        mono_mass: -0.0,
        average_mass: -0.0,
        ..Default::default()
    })
    .unwrap();
    assert_eq!(positive, negative);
    assert_eq!(positive.cmp(&negative), std::cmp::Ordering::Equal);
    assert_eq!(hash(&positive), hash(&negative));
}

#[test]
fn invalid_record_values_and_text_limits_are_checked() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            Ribonucleotide::from_record(RibonucleotideRecord {
                mono_mass: value,
                ..Default::default()
            })
            .is_err()
        );
        assert!(
            Ribonucleotide::from_record(RibonucleotideRecord {
                average_mass: value,
                ..Default::default()
            })
            .is_err()
        );
    }
    for code in [String::new(), "x".repeat(4097)] {
        assert!(
            Ribonucleotide::from_record(RibonucleotideRecord {
                code,
                ..Default::default()
            })
            .is_err()
        );
    }
    assert!(
        Ribonucleotide::from_record(RibonucleotideRecord {
            name: "x".repeat(65537),
            ..Default::default()
        })
        .is_err()
    );
    assert!(nucleotide("α", "", 0.0).is_modified());
}

#[test]
fn duplicate_codes_and_ambiguities_follow_final_code_map() {
    let db = RibonucleotideDB::from_entries(vec![
        entry("A", 1.0, None),
        entry("B", 2.0, None),
        entry("X?", 10.0, Some(["A", "B"])),
        entry("A", 3.0, None),
        entry("X?", 20.0, Some(["A", "missing"])),
    ])
    .unwrap();
    assert_eq!(db.len(), 5);
    assert_eq!(db.entries()[0].mono_mass(), 1.0);
    assert_eq!(db.get("A").unwrap().mono_mass(), 3.0);
    assert_eq!(db.get("X?").unwrap().mono_mass(), 20.0);
    let [first, second] = db.alternatives("X?").unwrap();
    assert_eq!(first.mono_mass(), 3.0);
    assert_eq!(second.mono_mass(), 2.0);
    assert_eq!(db.diagnostics().len(), 1);
    assert_eq!(db.diagnostics()[0].index, 5);
    let copy = db.clone();
    assert!(Arc::ptr_eq(&copy.get("A").unwrap(), &first));
    let ignored = entry("NoAmbiguity", 0.0, Some(["", "A"]));
    assert!(!ignored.is_ambiguous());
    assert!(
        RibonucleotideDB::from_entries(vec![ignored])
            .unwrap()
            .diagnostics()
            .is_empty()
    );
}

#[test]
fn successful_later_ambiguity_replaces_mapping_and_handles_outlive_registry() {
    let db = RibonucleotideDB::from_entries(vec![
        entry("A", 1.0, None),
        entry("B", 2.0, None),
        entry("X?", 0.0, Some(["A", "B"])),
        entry("X?", 0.0, Some(["B", "A"])),
    ])
    .unwrap();
    let [first, second] = db.alternatives("X?").unwrap();
    assert_eq!(first.code(), "B");
    assert_eq!(second.code(), "A");
    let weak = Arc::downgrade(&first);
    drop(db);
    assert_eq!(first.mono_mass(), 2.0);
    drop(first);
    assert!(weak.upgrade().is_none());
}

#[test]
fn longest_prefix_is_utf8_safe_and_is_independent_from_exact_lookup() {
    let db = RibonucleotideDB::from_records(vec![
        nucleotide("α", "", 1.0),
        nucleotide("αβ", "", 2.0),
        nucleotide("abc", "", 3.0),
    ])
    .unwrap();
    assert_eq!(db.get_prefix("αβγ").unwrap().code(), "αβ");
    assert_eq!(db.get_prefix("αx").unwrap().code(), "α");
    assert_eq!(db.get_prefix("abcdef").unwrap().code(), "abc");
    assert!(db.get("αβγ").is_err());
    assert!(db.get_prefix("").is_err());
    assert!(db.get_prefix("ab").is_err());
    assert!(db.get(&"a".repeat(4097)).is_err());
    assert!(db.alternatives(&"a".repeat(4097)).is_err());
}

#[test]
fn registry_limits_include_existing_record_and_alternative_payload() {
    let ribo = Arc::new(Ribonucleotide::default());
    let records = vec![
        RibonucleotideEntry {
            ribonucleotide: ribo,
            alternatives: None
        };
        100001
    ];
    assert!(RibonucleotideDB::from_entries(records).is_err());
    let mut record = entry("A", 0.0, None);
    record.alternatives = Some(["x".repeat(4097), "A".into()]);
    assert!(RibonucleotideDB::from_entries(vec![record]).is_err());
    // Shared Arcs still count their logical record payload for each input entry.
    let large = Arc::new(
        Ribonucleotide::from_record(RibonucleotideRecord {
            name: "x".repeat(65536),
            ..Default::default()
        })
        .unwrap(),
    );
    let records = vec![
        RibonucleotideEntry {
            ribonucleotide: large,
            alternatives: None
        };
        2100
    ];
    assert!(RibonucleotideDB::from_entries(records).is_err());
}

#[test]
fn tsv_reports_bad_rows_keeps_valid_eof_and_preserves_prime_and_mass_semantics() {
    let text = format!(
        "# source comment\n{HEADER}\textra\nwrong\nbad\tA\t\tA\t\tA\tC\tnan\t0\nname\t5′-p\tN\tX\tignored\thtml′\tH3PO4\tNone\t0"
    );
    let report = read_tsv(&text).unwrap();
    assert_eq!(report.skipped, 2);
    assert_eq!(
        report
            .diagnostics
            .iter()
            .map(|d| d.index)
            .collect::<Vec<_>>(),
        [3, 4]
    );
    assert_eq!(report.entries.len(), 1);
    let record = &report.entries[0].ribonucleotide;
    assert_eq!(record.code(), "5'-p");
    assert_eq!(record.html_code(), "html'");
    assert_eq!(
        record.term_specificity(),
        RibonucleotideTermSpecificity::FivePrime
    );
    assert_eq!(record.mono_mass(), 0.0);
    assert_eq!(record.average_mass(), record.formula().average_mass());
    assert!(read_tsv("no header").is_err());
    assert!(read_tsv("# only comment\n").is_err());
}

#[test]
fn tsv_ambiguity_first_last_fields_and_resource_errors() {
    let text = format!(
        "{HEADER}\nname\tX?\t\tX\t\t\t-\t0\t0\tA ignored B\nname\tX?\t\tX\t\t\t-\t0\t0\tmissing\n"
    );
    let report = read_tsv(&text).unwrap();
    assert_eq!(
        report.entries[0].alternatives,
        Some(["A".into(), "B".into()])
    );
    assert_eq!(report.skipped, 1);
    let oversize = format!("{HEADER}\nname\t{}\t\tA\t\t\t\t0\t0", "x".repeat(4097));
    assert!(read_tsv(&oversize).is_err());
    assert!(read_tsv(&" ".repeat(16 * 1024 * 1024 + 1)).is_err());
}

#[test]
fn global_registry_is_complete_without_optional_json_runtime() {
    let db = RibonucleotideDB::global();
    assert_eq!(db.len(), 378);
    assert_eq!(db.entries()[0].code(), "io6A");
    assert_eq!(db.get("A").unwrap().mono_mass(), 267.0968);
    assert_eq!(db.get("A").unwrap().average_mass(), 267.241);
    assert!(db.diagnostics().is_empty());
}

#[cfg(feature = "rna-json")]
#[test]
fn json_rows_report_type_errors_empty_codes_and_declared_zero() {
    let text = r#"{
      "2":{"name":"null abbrev","short_name":"B","reference_moiety":["A"],"formula":"C","abbrev":null},
      "10":{"name":"retained","short_name":"A","reference_moiety":["A"],"formula":"C","mass_monoiso":0,"mass_avg":-3},
      "3":{"name":"empty","short_name":"","reference_moiety":["A"],"formula":""},
      "4":{"name":"calculated","short_name":"C","reference_moiety":["A"],"formula":"C","mass_monoiso":null}
    }"#;
    let report = read_modomics_json(text).unwrap();
    assert_eq!(
        report
            .entries
            .iter()
            .map(|e| e.ribonucleotide.code())
            .collect::<Vec<_>>(),
        ["A", "C"]
    );
    assert_eq!(report.entries[0].ribonucleotide.mono_mass(), 0.0);
    assert_eq!(report.entries[0].ribonucleotide.average_mass(), -3.0);
    assert_eq!(report.entries[1].ribonucleotide.mono_mass(), 12.0);
    assert_eq!(report.entries[1].ribonucleotide.average_mass(), 0.0);
    assert_eq!(report.skipped, 2);
    assert_eq!(
        report
            .diagnostics
            .iter()
            .map(|d| d.index)
            .collect::<Vec<_>>(),
        [2, 3]
    );
}

#[cfg(feature = "rna-json")]
#[test]
fn json_document_bounds_and_shape_errors_are_checked_before_record_output() {
    assert!(read_modomics_json("{").is_err());
    assert!(read_modomics_json(&format!("{}0{}", "[".repeat(65), "]".repeat(65))).is_err());
    assert!(read_modomics_json(&format!("[{}0]", "0,".repeat(250001))).is_err());
    assert!(read_modomics_json(&format!("\"{}\"", "x".repeat(6 * 65536 + 1))).is_err());
    assert!(read_modomics_json("null").unwrap().entries.is_empty());
    assert_eq!(read_modomics_json("3").unwrap().skipped, 1);
    let malformed = r#"[{"name":"bad","short_name":"A?","reference_moiety":["A"],"formula":"C","alternatives":[1,2]}, {"name":"bad","short_name":"C","reference_moiety":["A","C"],"formula":"C"}]"#;
    assert_eq!(read_modomics_json(malformed).unwrap().skipped, 2);
}

#[cfg(feature = "rna-json")]
#[test]
fn json_and_tsv_preserve_their_distinct_empty_first_alternative_behavior() {
    let json = r#"[{"name":"ambiguous","short_name":"X?","reference_moiety":["X"],"formula":"","alternatives":["","G"]}]"#;
    let report = read_modomics_json(json).unwrap();
    assert_eq!(report.entries.len(), 1);
    assert!(report.entries[0].alternatives.is_none());
    let tsv = format!("{HEADER}\nambiguous\tX?\t\tX\t\t\t\t0\t0\t G");
    let report = read_tsv(&tsv).unwrap();
    assert_eq!(
        report.entries[0].alternatives,
        Some(["".into(), "G".into()])
    );
    assert!(!report.entries[0].is_ambiguous());

    // Source ignores the alternatives field for codes outside its ambiguity
    // branch. A field that is not consumed must not trigger its code-size cap.
    let ignored = "x".repeat(5000);
    let json = format!(
        r#"[{{"name":"plain","short_name":"A","reference_moiety":["A"],"formula":"","alternatives":["{ignored}","G"]}}]"#
    );
    assert!(
        read_modomics_json(&json).unwrap().entries[0]
            .alternatives
            .is_none()
    );
    let tsv = format!("{HEADER}\ndeoxy\tdX?\t\tX\t\t\t\t0\t0\t{ignored} G");
    assert!(read_tsv(&tsv).unwrap().entries[0].alternatives.is_none());
}
