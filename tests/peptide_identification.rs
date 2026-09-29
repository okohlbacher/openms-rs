// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//
// PeptideIdentification members beyond the record (UIDs, USIs, the score
// comparator, the source empty(), label and reference text, hashing) and the
// USI class they return.
//
// Expected values are EXECUTED C++ output (tier 1) unless a test says
// otherwise: ../oracle/B8b-peptide-id/driver.cpp, compiled with libOpenMS's
// Release flags against the Linux x86_64 Release install
// openms4-release-bc9cc12-c19e494-174b576 and run twice on ibminode06, both
// runs byte-identical. tests/data/peptide_identification_oracle.tsv holds its
// stdout, every line prefixed with the case name; the manifest
// tests/data/peptide_identification_provenance.json records the hashes. The
// consensus-map case is replayed by tests/peptide_identification_consensus.rs.
// Tests marked "source literal" transcribe PeptideIdentification_test.cpp at
// bc9cc12 (tier 3).

use openms::Error;
use openms::chemistry::AASequence;
use openms::concept::boost_regex::BoostRegex;
use openms::identification::usi::{Usi, UsiIndexType};
use openms::identification::{
    IdentifierMSRunMapper, PeptideEvidence, PeptideHit, PeptideIdentification,
    ProteinIdentification, SpectrumIdentification,
};
use openms::metadata::MetaValue;
use std::collections::hash_map::DefaultHasher;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::hash::{Hash, Hasher};

const ORACLE: &str = include_str!("data/peptide_identification_oracle.tsv");

/// The oracle lines of one case, split into tab-separated fields (the case
/// name removed).
fn case(name: &str) -> Vec<Vec<&'static str>> {
    let rows: Vec<Vec<&str>> = ORACLE
        .lines()
        .filter_map(|line| line.strip_prefix(name)?.strip_prefix('\t'))
        .map(|rest| rest.split('\t').collect())
        .collect();
    assert!(!rows.is_empty(), "oracle case {name} missing");
    rows
}

/// The value inside `prefix[...]` of an oracle label.
fn bracketed<'a>(label: &'a str, prefix: &str) -> Option<&'a str> {
    label
        .strip_prefix(prefix)?
        .strip_prefix('[')?
        .strip_suffix(']')
}

fn flag(value: bool) -> &'static str {
    if value { "1" } else { "0" }
}

/// The driver's `print_usi` record for `usi`, without its label.
fn usi_fields(usi: &Usi) -> Vec<String> {
    vec![
        format!("valid={}", flag(usi.is_valid())),
        format!("collection={}", usi.collection),
        format!("ms_run={}", usi.ms_run),
        format!("type={}", usi.index_type),
        format!("index={}", usi.index),
        format!("interpretation={}", usi.interpretation),
        format!("text={usi}"),
    ]
}

#[track_caller]
fn assert_usi(row: &[&str], usi: &Usi) {
    assert_eq!(&row[1..], usi_fields(usi).as_slice(), "{}", row[0]);
}

fn hit(sequence: &str, charge: i32, score: f64) -> PeptideHit {
    PeptideHit::new(score, 0, charge, AASequence::parse(sequence).unwrap()).unwrap()
}

fn hash_of<T: Hash>(value: &T) -> u64 {
    let mut state = DefaultHasher::new();
    value.hash(&mut state);
    state.finish()
}

// ---------------------------------------------------------------------------
// USI class (case usi_class).
// ---------------------------------------------------------------------------

#[test]
fn usi_constructors_validity_and_text_match_the_release_build() {
    let rows = case("usi_class");
    let find = |label: &str| rows.iter().find(|row| row[0] == label).unwrap();
    assert_usi(find("default"), &Usi::default());
    assert_usi(
        find("ctor"),
        &Usi::new("PXD1", "run.mzML", UsiIndexType::Index, "7", "PEPTIDE/2"),
    );
    assert_usi(
        find("ctor_missing_index"),
        &Usi::new("PXD1", "run.mzML", UsiIndexType::Scan, "", ""),
    );
    assert_usi(
        find("from_scan"),
        &Usi::from_scan_number("PXD1", "f.mzML", -3, "PEPTIDE/2"),
    );
    assert_usi(
        find("from_native_scan"),
        &Usi::from_native_id(
            "PXD1",
            "f.mzML",
            "controllerType=0 controllerNumber=1 scan=5",
        ),
    );
    assert_usi(
        find("from_native_other"),
        &Usi::from_native_id("PXD1", "f.mzML", "custom:thing"),
    );
    assert_eq!(find("cv")[1..], [Usi::CV_ACCESSION, Usi::CV_NAME]);
    assert_eq!(
        find("stream")[1],
        Usi::new("PXD1", "r.mzML", UsiIndexType::Scan, "3", "").to_string()
    );
}

#[test]
fn usi_parsing_keeps_the_source_partial_state_and_error_text() {
    let rows = case("usi_class");
    let mut checked = 0;
    for row in &rows {
        if let Some(text) = bracketed(row[0], "parse") {
            let mut parsed = Usi::default();
            let accepted = parsed.from_string(text);
            assert_eq!(
                row[1..],
                [
                    format!("from_string={}", flag(accepted)),
                    format!("is_valid_usi={}", flag(Usi::is_valid_usi(text))),
                    format!("try_parse={}", flag(Usi::try_parse(text).is_some())),
                ],
                "{text}"
            );
            let state = rows
                .iter()
                .find(|other| bracketed(other[0], "parsed") == Some(text))
                .unwrap();
            assert_usi(state, &parsed);
            let constructed = rows
                .iter()
                .find(|other| bracketed(other[0], "ctor") == Some(text))
                .unwrap();
            match Usi::parse(text) {
                Ok(usi) => {
                    assert_eq!(constructed[1], "ok", "{text}");
                    assert_eq!(usi, parsed);
                    assert_eq!(text.parse::<Usi>().unwrap(), usi);
                }
                Err(Error::Parse { message, .. }) => {
                    assert_eq!(constructed[1..], ["throw", "Parse Error", &message]);
                }
                Err(other) => panic!("{text}: {other}"),
            }
            checked += 1;
        }
    }
    assert_eq!(checked, 20);
}

#[test]
fn usi_index_type_names_are_case_insensitive_and_refuse_others() {
    let rows = case("usi_class");
    let mut checked = 0;
    for row in &rows {
        let Some(text) = bracketed(row[0], "index_type") else {
            continue;
        };
        match UsiIndexType::parse(text) {
            Ok(index_type) => {
                assert_eq!(row[1], index_type.as_str(), "{text}");
                assert_eq!(text.parse::<UsiIndexType>().unwrap(), index_type);
            }
            // The driver printed its label before the call threw.
            Err(Error::InvalidValue(message)) => {
                assert_eq!(row[2..], ["throw", "InvalidValue", &message], "{text}")
            }
            Err(other) => panic!("{text}: {other}"),
        }
        checked += 1;
    }
    assert_eq!(checked, 9);
}

#[test]
fn usi_basenames_follow_the_source_separator_and_uri_rules() {
    let rows = case("usi_class");
    let mut checked = 0;
    for row in &rows {
        if let Some(path) = bracketed(row[0], "basename") {
            assert_eq!(
                format!("[{}]", Usi::extract_basename(path)),
                row[1],
                "{path}"
            );
            checked += 1;
        }
    }
    assert_eq!(checked, 15);
}

#[test]
fn usi_scan_numbers_follow_the_native_id_patterns_and_to_int32() {
    let rows = case("usi_class");
    let mut checked = 0;
    for row in &rows {
        if let Some(native_id) = bracketed(row[0], "scan_number") {
            let expected = match row[1] {
                "none" => None,
                value => Some(value.parse::<i32>().unwrap()),
            };
            assert_eq!(
                Usi::extract_scan_number_from_native_id(native_id),
                expected,
                "{native_id:?}"
            );
            checked += 1;
        }
    }
    assert_eq!(checked, 32);
}

/// The port matches the source's `prefix(?<GROUP>\d+)` patterns without a
/// regular-expression engine. This compiles the source patterns with the
/// crate's Boost-compatible engine and checks that its last group-1 token
/// converts to the same scan number on every oracle input (independently
/// derived, tier 4, on top of the tier-1 values above).
#[test]
fn usi_scan_prefix_matching_agrees_with_the_boost_patterns() {
    fn pattern(id: &str) -> &'static str {
        if ["scan=", "controllerType=", "function=", "frame="]
            .iter()
            .any(|p| id.starts_with(p))
        {
            r"scan=(?<GROUP>\d+)"
        } else if id.starts_with("index=") {
            r"index=(?<GROUP>\d+)"
        } else if id.starts_with("scanId=") {
            r"scanId=(?<GROUP>\d+)"
        } else if id.starts_with("scanID=") {
            r"scanID=(?<GROUP>\d+)"
        } else if id.starts_with("spectrum=") {
            r"spectrum=(?<GROUP>\d+)"
        } else if id.starts_with("file=") {
            r"file=(?<GROUP>\d+)"
        } else {
            r"(?<GROUP>\d+)"
        }
    }
    let native = [
        "scan=",
        "scanId=",
        "scanID=",
        "controllerType=",
        "function=",
        "sample=",
        "index=",
        "spectrum=",
        "file=",
        "frame=",
    ];
    let rows = case("usi_class");
    for row in &rows {
        let Some(id) = bracketed(row[0], "scan_number") else {
            continue;
        };
        let regex = BoostRegex::new(pattern(id)).unwrap();
        let last = regex
            .tokens(id.as_bytes(), &[1])
            .unwrap()
            .map(|token| token.unwrap())
            .filter(|token| token.matched)
            .last()
            .map(|token| String::from_utf8(token.as_bytes(id.as_bytes()).to_vec()).unwrap());
        let from_regex = native
            .iter()
            .any(|p| id.starts_with(p))
            .then_some(last)
            .flatten()
            .and_then(|digits| digits.parse::<i32>().ok());
        let direct = Usi::extract_scan_number_from_native_id(id);
        if let Some(scan) = from_regex {
            assert_eq!(direct, Some(scan), "{id:?}");
        }
    }
}

#[test]
fn usi_ordering_is_the_source_operator_less() {
    let rows = case("usi_class");
    let order = [
        Usi::new("b", "r", UsiIndexType::Scan, "1", ""),
        Usi::new("a", "r", UsiIndexType::NativeId, "1", ""),
        Usi::new("a", "r", UsiIndexType::Scan, "2", ""),
        Usi::new("a", "r", UsiIndexType::Index, "1", ""),
        Usi::new("a", "r", UsiIndexType::Scan, "10", ""),
        Usi::new("a", "r", UsiIndexType::Scan, "2", "X"),
        Usi::new("a", "q", UsiIndexType::Scan, "9", ""),
    ];
    let less: Vec<&Vec<&str>> = rows.iter().filter(|row| row[0] == "less").collect();
    assert_eq!(less.len(), order.len());
    for (i, row) in less.iter().enumerate() {
        assert_eq!(row[1], i.to_string());
        for (j, other) in order.iter().enumerate() {
            assert_eq!(row[2 + j], flag(order[i] < *other), "{i} < {j}");
        }
        let next = &order[(i + 1) % order.len()];
        assert_eq!(
            row[2 + order.len()],
            format!("{}{}", flag(order[i] == order[i]), flag(order[i] != *next))
        );
    }
}

// ---------------------------------------------------------------------------
// buildUSI (cases usi_build and usi_mapping).
// ---------------------------------------------------------------------------

#[test]
fn build_usi_reads_scan_numbers_and_the_first_hit_as_proforma() {
    let rows = case("usi_build");
    let mut checked = 0;
    for row in &rows {
        let label = row[0];
        let usi = if let Some(reference) = bracketed(label, "ref") {
            let mut id = PeptideIdentification::new();
            id.set_spectrum_reference(reference);
            id.hits = vec![hit("EM(Oxidation)K", 2, 100.0), hit("PEPTIDEK", 3, 50.0)];
            id.build_usi("sample.mzML", "PXD000561", false).unwrap()
        } else if let Some(reference) = bracketed(label, "ref_interp") {
            let mut id = PeptideIdentification::new();
            id.set_spectrum_reference(reference);
            // Not sorted: the first hit is the interpretation, whatever its score.
            id.hits = vec![hit("EM(Oxidation)K", 2, 100.0), hit("PEPTIDEK", 3, 50.0)];
            id.build_usi("sample.mzML", "PXD000561", true).unwrap()
        } else if let Some(spec) = bracketed(label, "interp") {
            let (sequence, charge) = spec.rsplit_once('/').unwrap();
            let mut id = PeptideIdentification::new();
            id.set_spectrum_reference("scan=1");
            id.hits = vec![hit(sequence, charge.parse().unwrap(), 1.0)];
            id.build_usi("s.mzML", Usi::LOCAL_COLLECTION, true).unwrap()
        } else {
            let mut id = PeptideIdentification::new();
            id.set_spectrum_reference("scan=3");
            let interpretation = match label {
                "interp_nohits" => true,
                "default_dataset" => false,
                other => panic!("unexpected oracle label {other}"),
            };
            id.build_usi("s.mzML", Usi::LOCAL_COLLECTION, interpretation)
                .unwrap()
        };
        assert_usi(row, &usi);
        checked += 1;
    }
    assert_eq!(checked, 65);
}

#[test]
fn build_usi_with_mapping_resolves_merged_runs_and_legacy_base_names() {
    let mut run1 = ProteinIdentification::new();
    run1.identifier = "run1".into();
    run1.primary_ms_run_paths = vec!["/data/file1.mzML".into()];
    let mut run2 = ProteinIdentification::new();
    run2.identifier = "run2".into();
    run2.primary_ms_run_paths = vec![
        "/data/fileA.mzML".into(),
        "file:///C:/data/fileB.mzML".into(),
    ];
    let mut run3 = ProteinIdentification::new();
    run3.identifier = "run3".into();
    let mapping = IdentifierMSRunMapper::from_runs(&[run1, run2, run3]).unwrap();

    let rows = case("usi_mapping");
    let expect = |label: &str, id: &PeptideIdentification| {
        let row = rows.iter().find(|row| row[0] == label).unwrap();
        match id.build_usi_with_mapping(&mapping, "PXD000561", true) {
            Ok(usi) => assert_usi(row, &usi),
            // Source: ConversionError for a negative id_merge_index.
            Err(error) => {
                assert_eq!(row[1..3], ["throw", "ConversionError"], "{label}: {error}")
            }
        }
    };
    let mut id = PeptideIdentification::new();
    id.identifier = "run1".into();
    id.set_spectrum_reference("scan=100");
    expect("single", &id);
    id.identifier = "run2".into();
    expect("merged_default", &id);
    id.metadata.insert("id_merge_index".into(), 1_i64.into());
    expect("merged_1", &id);
    id.metadata.insert("id_merge_index".into(), 99_i64.into());
    expect("merged_99", &id);
    id.metadata.insert("id_merge_index".into(), (-1_i64).into());
    expect("merged_neg", &id);
    id.metadata.remove("id_merge_index");
    id.identifier = "unknown".into();
    expect("unknown", &id);
    id.metadata
        .insert("base_name".into(), "C:\\raw\\legacy.mzML".into());
    expect("unknown_base_name", &id);
    id.identifier = "run3".into();
    expect("no_paths_base_name", &id);
    let mut bare = PeptideIdentification::new();
    bare.identifier = "run1".into();
    expect("no_reference", &bare);
    assert_eq!(rows.len(), 9);
}

// Source literals: PeptideIdentification_test.cpp, buildUSI sections.
#[test]
fn build_usi_source_literals() {
    let mut id = PeptideIdentification::new();
    id.set_spectrum_reference("controllerType=0 controllerNumber=1 scan=12345");
    id.rt = Some(100.0);
    id.mz = Some(500.0);
    let usi = id.build_usi("sample.mzML", "PXD000561", false).unwrap();
    assert!(usi.is_valid());
    assert_eq!(usi.index, "12345");
    assert_eq!(usi.index_type, UsiIndexType::Scan);
    let local = id
        .build_usi("sample.mzML", Usi::LOCAL_COLLECTION, false)
        .unwrap();
    assert_eq!(local.collection, "local");
    id.hits.push(hit("PEPTIDEK", 2, 100.0));
    id.sort().unwrap();
    assert_eq!(
        id.build_usi("sample.mzML", "PXD000561", true)
            .unwrap()
            .interpretation,
        "PEPTIDEK/2"
    );
    let empty = PeptideIdentification::new()
        .build_usi("test.mzML", "PXD000564", false)
        .unwrap();
    assert!(!empty.is_valid());
    assert_eq!(empty.to_string(), "");
    assert!(!Usi::is_valid_usi(&empty.to_string()));
    let uri = "file:///C:/Users/bielow/data/ES-0014b_2.mzML";
    let mut scan = PeptideIdentification::new();
    scan.set_spectrum_reference("scan=12345");
    assert_eq!(
        scan.build_usi(Usi::extract_basename(uri), "PXD000561", false)
            .unwrap()
            .to_string(),
        "mzspec:PXD000561:ES-0014b_2.mzML:scan:12345"
    );
}

// ---------------------------------------------------------------------------
// buildUIDFromPepID error paths (case uid_errors).
// ---------------------------------------------------------------------------

#[test]
fn build_uid_uses_the_single_run_path_or_map_index_and_reports_missing_information() {
    let mut one = ProteinIdentification::new();
    one.identifier = "one".into();
    one.primary_ms_run_paths = vec!["/data/a.mzML".into()];
    let mut two = ProteinIdentification::new();
    two.identifier = "two".into();
    two.primary_ms_run_paths = vec!["/data/b.mzML".into(), "/data/c.mzML".into()];
    let mapping = IdentifierMSRunMapper::from_runs(&[one, two]).unwrap();

    let rows = case("uid_errors");
    let expect = |label: &str, id: &PeptideIdentification| {
        let row = rows.iter().find(|row| row[0] == label).unwrap();
        match id.build_uid(&mapping) {
            Ok(uid) => assert_eq!(row[1..], ["ok", uid.as_str()], "{label}"),
            // The driver printed "<label>\tok\t" before the call threw.
            Err(Error::MissingInformation(message)) => assert_eq!(
                row[1..],
                ["ok", label, "throw", "MissingInformation", &message],
                "{label}"
            ),
            Err(other) => panic!("{label}: {other}"),
        }
    };
    let mut id = PeptideIdentification::new();
    id.identifier = "one".into();
    expect("missing_reference", &id);
    id.set_spectrum_reference("scan=5");
    expect("single_run", &id);
    id.identifier = "two".into();
    expect("merged_no_map_index", &id);
    id.metadata.insert("map_index".into(), 3_i64.into());
    expect("merged_map_index_int", &id);
    id.metadata.insert("map_index".into(), "x7".into());
    expect("merged_map_index_string", &id);
    id.metadata
        .insert("map_index".into(), MetaValue::try_from(2.5).unwrap());
    expect("merged_map_index_double", &id);
    id.identifier = "unknown".into();
    expect("unknown_run_with_map_index", &id);
    id.metadata.remove("map_index");
    expect("unknown_run_without_map_index", &id);
    let mut int_reference = PeptideIdentification::new();
    int_reference.identifier = "one".into();
    int_reference
        .metadata
        .insert("spectrum_reference".into(), 17_i64.into());
    expect("int_reference", &int_reference);
    int_reference.set_spectrum_reference("");
    expect("empty_reference", &int_reference);
    assert_eq!(rows.len(), 10);
}

// ---------------------------------------------------------------------------
// empty(), labels, references and sorting (cases empty, labels, sort).
// ---------------------------------------------------------------------------

#[test]
fn source_empty_consults_run_hits_threshold_score_type_and_orientation_only() {
    let rows = case("empty");
    let expect = |label: &str, id: &PeptideIdentification| {
        let row = rows.iter().find(|row| row[0] == label).unwrap();
        assert_eq!(row[1], flag(id.source_empty()), "{label}");
    };
    let mut id = PeptideIdentification::new();
    expect("default", &id);
    id.identifier = "x".into();
    expect("identifier", &id);
    id = PeptideIdentification::new();
    id.score_type = "q".into();
    expect("score_type", &id);
    id = PeptideIdentification::new();
    id.higher_score_better = false;
    expect("lower_better", &id);
    id = PeptideIdentification::new();
    id.significance_threshold = 1.0;
    expect("threshold1", &id);
    id.significance_threshold = 0.0;
    expect("threshold0", &id);
    id = PeptideIdentification::new();
    id.rt = Some(5.0);
    id.mz = Some(6.0);
    id.metadata.insert("foo".into(), 1_i64.into());
    expect("rt_mz_meta", &id);
    id.hits.push(PeptideHit::default());
    expect("hit", &id);
    assert!(!id.is_empty());
    assert_eq!(rows.len(), 8);
}

#[test]
fn labels_and_references_are_stringified_as_the_source_does() {
    let rows = case("labels");
    let expect = |label: &str, value: String| {
        let row = rows.iter().find(|row| row[0] == label).unwrap();
        assert_eq!(row[1], format!("[{value}]"), "{label}");
    };
    let mut id = PeptideIdentification::new();
    expect("label_default", id.experiment_label());
    id.set_experiment_label("A");
    expect("label_A", id.experiment_label());
    // An empty label is ignored: the earlier label and its entry stay.
    id.set_experiment_label("");
    expect("label_after_empty", id.experiment_label());
    assert!(id.metadata.contains_key("experiment_label"));
    id.metadata.insert("experiment_label".into(), 12_i64.into());
    expect("label_int", id.experiment_label());

    let mut reference = PeptideIdentification::new();
    expect("ref_default", reference.spectrum_reference());
    for (label, value) in [
        ("ref_int0", MetaValue::from(0_i64)),
        ("ref_int_neg", MetaValue::from(-42_i64)),
        ("ref_double", MetaValue::try_from(1.5).unwrap()),
        (
            "ref_double_big",
            MetaValue::try_from(123_456_789.125).unwrap(),
        ),
        (
            "ref_string_list",
            MetaValue::from(vec!["a".to_owned(), "b".to_owned()]),
        ),
    ] {
        reference
            .metadata
            .insert("spectrum_reference".into(), value);
        expect(label, reference.spectrum_reference());
    }
    reference.set_spectrum_reference("");
    expect("ref_set_empty", reference.spectrum_reference());
    assert!(reference.metadata.contains_key("spectrum_reference"));
}

#[test]
fn sort_is_stable_and_the_comparator_is_strict() {
    let rows = case("sort");
    let mut id = PeptideIdentification::new();
    id.hits = vec![
        hit("AAA", 1, 5.0),
        hit("CCC", 1, 7.0),
        hit("DDD", 1, 5.0),
        hit("EEE", 1, 9.0),
        hit("FFF", 1, 7.0),
    ];
    let order = |id: &PeptideIdentification| -> Vec<String> {
        id.hits.iter().map(|hit| hit.sequence.to_string()).collect()
    };
    id.sort().unwrap();
    assert_eq!(order(&id), rows[0][1..]);
    assert_eq!(rows[0][0], "higher");
    id.higher_score_better = false;
    id.sort().unwrap();
    assert_eq!(order(&id), rows[1][1..]);
    assert_eq!(rows[1][0], "lower");

    let (a, b, a2) = (hit("A", 1, 1.0), hit("B", 1, 2.0), hit("C", 1, 1.0));
    for (row, higher) in [(&rows[2], true), (&rows[3], false)] {
        let compare = PeptideIdentification::score_comparator(higher);
        assert_eq!(
            row[1],
            format!(
                "{}{}{}",
                flag(compare(&a, &b)),
                flag(compare(&b, &a)),
                flag(compare(&a, &a2))
            )
        );
    }
}

// Source literals: PeptideIdentification_test.cpp, getReferencingHits.
#[test]
fn referencing_hits_is_available_over_any_hit_slice() {
    let mut hits = Vec::new();
    for (sequence, score, accession) in [
        ("FIRSTPROTEIN", 23.0, "TEST_PROTEIN1"),
        ("SECONDPROTEIN", 10.0, "TEST_PROTEIN2"),
        ("THIRDPROTEIN", 11.0, "TEST_PROTEIN3"),
    ] {
        let mut candidate = hit(sequence, 1, score);
        candidate.evidences.push(PeptideEvidence {
            protein_accession: accession.into(),
            ..PeptideEvidence::default()
        });
        hits.push(candidate);
    }
    let query = BTreeSet::from(["TEST_PROTEIN2".to_owned(), "TEST_PROTEIN3".to_owned()]);
    let found: Vec<PeptideHit> = PeptideIdentification::hits_referencing(&hits, &query)
        .into_iter()
        .cloned()
        .collect();
    assert_eq!(found.len(), 2);
    assert_eq!(found[0].sequence.to_string(), "SECONDPROTEIN");
    assert_eq!(found[1].sequence.to_string(), "THIRDPROTEIN");
    let id = PeptideIdentification {
        hits,
        ..PeptideIdentification::default()
    };
    let own: Vec<PeptideHit> = id.referencing_hits(&query).into_iter().cloned().collect();
    assert_eq!(own, found);
}

// Source literals: PeptideIdentification_test.cpp, std::hash<PeptideIdentification>.
#[test]
fn hash_is_consistent_with_equality_and_keys_hash_containers() {
    let make = || {
        let mut id = SpectrumIdentification::new();
        id.rt = Some(123.456);
        id.mz = Some(789.012);
        id.identifier = "test_id".into();
        id.score_type = "Mascot".into();
        id.higher_score_better = true;
        id.significance_threshold = 0.05;
        id
    };
    let (id1, id2) = (make(), make());
    assert_eq!(id1, id2);
    assert_eq!(hash_of(&id1), hash_of(&id2));
    let mut id3 = PeptideIdentification::new();
    id3.rt = Some(999.999);
    id3.mz = Some(111.111);
    id3.identifier = "different_id".into();
    assert_ne!(id1, id3);
    assert_ne!(hash_of(&id1), hash_of(&id3));
    let set: HashSet<PeptideIdentification> = [id1.clone(), id2.clone(), id3.clone()].into();
    assert_eq!(set.len(), 2);
    let mut map = HashMap::new();
    map.insert(id1.clone(), "first");
    map.insert(id3.clone(), "third");
    assert_eq!(map[&id2], "first");
    assert_eq!(map[&id3], "third");

    let with_hit = |sequence: &str, charge: i32, score: f64| PeptideIdentification {
        hits: vec![hit(sequence, charge, score)],
        ..PeptideIdentification::default()
    };
    let (id4, id5) = (with_hit("PEPTIDE", 2, 50.0), with_hit("PEPTIDE", 2, 50.0));
    assert_eq!(id4, id5);
    assert_eq!(hash_of(&id4), hash_of(&id5));
    let id6 = with_hit("DIFFERENT", 3, 99.0);
    assert_ne!(id4, id6);
    assert_ne!(hash_of(&id4), hash_of(&id6));
    // Unset coordinates are equal and hash equally (the source's NaN rule).
    let (id7, id8) = (PeptideIdentification::new(), PeptideIdentification::new());
    assert!(id7.rt.is_none() && id7.mz.is_none());
    assert_eq!(id7, id8);
    assert_eq!(hash_of(&id7), hash_of(&id8));
    // Independently derived: the two zeros are equal, so they hash equally.
    let mut zero = PeptideIdentification::new();
    zero.rt = Some(0.0);
    let mut negative_zero = PeptideIdentification::new();
    negative_zero.rt = Some(-0.0);
    assert_eq!(zero, negative_zero);
    assert_eq!(hash_of(&zero), hash_of(&negative_zero));
    // A present zero differs from an absent coordinate.
    assert_ne!(zero, PeptideIdentification::new());
}

// Source literals: USI_test.cpp at bc9cc12 (tier 3), one assertion per value.
#[test]
fn usi_class_test_literals() {
    let usi = Usi::parse("mzspec:PXD000561:sample.mzML:scan:12345:PEPTIDEK/2").unwrap();
    assert!(usi.is_valid());
    assert_eq!(
        (
            usi.collection.as_str(),
            usi.ms_run.as_str(),
            usi.index.as_str()
        ),
        ("PXD000561", "sample.mzML", "12345")
    );
    assert_eq!(usi.interpretation, "PEPTIDEK/2");
    let colon = Usi::parse("mzspec:PXD000561:sample.mzML:scan:12345:EM[UNIMOD:35]K/2").unwrap();
    assert_eq!(colon.interpretation, "EM[UNIMOD:35]K/2");
    assert_eq!(
        colon.to_string(),
        "mzspec:PXD000561:sample.mzML:scan:12345:EM[UNIMOD:35]K/2"
    );
    assert!(matches!(
        Usi::parse("invalid_usi_string"),
        Err(Error::Parse { .. })
    ));
    assert!(matches!(
        Usi::parse("mzspec:PXD000561"),
        Err(Error::Parse { .. })
    ));
    assert_eq!(
        Usi::new("PXD000561", "sample.mzML", UsiIndexType::Index, "0", "").to_string(),
        "mzspec:PXD000561:sample.mzML:index:0"
    );
    let mut parsed = Usi::default();
    assert!(parsed.from_string("mzspec:PXD000561:sample.mzML:index:0"));
    assert_eq!(parsed.index_type, UsiIndexType::Index);
    assert!(!parsed.from_string("invalid"));
    assert!(!parsed.from_string(""));
    assert_eq!(
        Usi::from_scan_number("PXD000561", "sample.mzML", 12345, "EM[UNIMOD:35]K/2").to_string(),
        "mzspec:PXD000561:sample.mzML:scan:12345:EM[UNIMOD:35]K/2"
    );
    assert_eq!(
        Usi::from_native_id("PXD000561", "sample.mzML", "scan=12345").to_string(),
        "mzspec:PXD000561:sample.mzML:scan:12345"
    );
    let fallback = Usi::from_native_id("PXD000561", "sample.mzML", "custom_id_format");
    assert_eq!(fallback.index_type, UsiIndexType::NativeId);
    assert_eq!(fallback.index, "custom_id_format");
    for (native_id, scan) in [
        ("scan=12345", Some(12345)),
        (
            "controllerType=0 controllerNumber=1 scan=54321",
            Some(54321),
        ),
        ("scanId=99999", Some(99999)),
        ("index=100", Some(100)),
        ("12345", Some(12345)),
        ("custom_format_no_number", None),
    ] {
        assert_eq!(Usi::extract_scan_number_from_native_id(native_id), scan);
    }
    assert_eq!(UsiIndexType::parse("INDEX").unwrap(), UsiIndexType::Index);
    assert!(UsiIndexType::parse("invalid_type").is_err());
    for (path, base) in [
        ("/data/experiments/run1.mzML", "run1.mzML"),
        ("/sample.mzML", "sample.mzML"),
        ("C:\\Users\\test\\data\\run1.mzML", "run1.mzML"),
        ("file:///path/to/sample.mzML", "sample.mzML"),
        (
            "file:///C:/data/experiments/ES-0014b_2.mzML",
            "ES-0014b_2.mzML",
        ),
        ("run1.raw", "run1.raw"),
        ("/data/sub\\folder/sample.mzML", "sample.mzML"),
    ] {
        assert_eq!(Usi::extract_basename(path), base);
    }
    assert_eq!(Usi::CV_ACCESSION, "MS:1003063");
    assert_eq!(Usi::CV_NAME, "universal spectrum identifier");
}
