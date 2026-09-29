// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//
// ProteaseDB (CHEMISTRY/ProteaseDB.h and its DigestionEnzymeDB base).
//
// Expected values are EXECUTED C++ output (tier 1): ../oracle/B8b-dbs/driver.cpp,
// compiled with libOpenMS's Release flags against the Linux x86_64 Release
// install openms4-release-bc9cc12-c19e494-174b576 and run twice on ibminode06,
// both runs byte-identical. tests/data/protease_db_oracle.tsv is the "protease"
// case's stdout, tests/data/protease_db_oracle_writetsv.tsv the file
// ProteaseDB::writeTSV wrote; tests/data/protease_db_provenance.json records
// the hashes. Tests marked "source literal" transcribe ProteaseDB_test.cpp at
// bc9cc12 (tier 3).

use openms::chemistry::{EmpiricalFormula, Protease, ProteaseDB};
use std::collections::{BTreeMap, BTreeSet};

const ORACLE: &str = include_str!("data/protease_db_oracle.tsv");
const ORACLE_TSV: &str = include_str!("data/protease_db_oracle_writetsv.tsv");

/// Oracle rows starting with `tag`, split on tabs, the tag removed.
fn rows(tag: &str) -> Vec<Vec<&'static str>> {
    let rows: Vec<Vec<&str>> = ORACLE
        .lines()
        .filter_map(|line| line.strip_prefix(tag)?.strip_prefix('\t'))
        .map(|rest| rest.split('\t').collect())
        .collect();
    assert!(!rows.is_empty(), "oracle has no {tag} rows");
    rows
}

fn single(tag: &str) -> Vec<&'static str> {
    let rows = rows(tag);
    assert_eq!(rows.len(), 1, "{tag}");
    rows[0][0].split('|').collect()
}

fn cpp_id(text: &str) -> Option<i32> {
    match text.parse::<i32>().unwrap() {
        -1 => None,
        id => Some(id),
    }
}

fn formula(text: &str) -> EmpiricalFormula {
    if text.is_empty() {
        EmpiricalFormula::default()
    } else {
        text.parse().unwrap()
    }
}

#[test]
fn every_enzyme_definition_matches_the_source_singleton() {
    let db = ProteaseDB::global();
    let count: usize = rows("count")[0][0].parse().unwrap();
    assert_eq!(count, 33);
    assert_eq!(db.names().len(), count);
    let mut seen = BTreeSet::new();
    for row in rows("iter") {
        let enzyme = db.get_enzyme(row[1]).unwrap();
        assert_eq!(enzyme.name(), row[1]);
        assert_eq!(enzyme.regex(), row[2], "{}", row[1]);
        let synonyms: Vec<&str> = if row[3].is_empty() {
            Vec::new()
        } else {
            row[3].split('|').collect()
        };
        assert_eq!(enzyme.synonyms(), synonyms.as_slice(), "{}", row[1]);
        assert_eq!(enzyme.xtandem_id(), row[4], "{}", row[1]);
        assert_eq!(enzyme.comet_id(), cpp_id(row[5]), "{}", row[1]);
        assert_eq!(enzyme.msgf_id(), cpp_id(row[6]), "{}", row[1]);
        assert_eq!(enzyme.omssa_id(), cpp_id(row[7]), "{}", row[1]);
        assert_eq!(enzyme.psi_id(), row[8], "{}", row[1]);
        assert_eq!(enzyme.n_term_gain(), formula(row[9]), "{}", row[1]);
        assert_eq!(enzyme.c_term_gain(), formula(row[10]), "{}", row[1]);
        assert_eq!(enzyme.description(), row[11], "{}", row[1]);
        assert!(seen.insert(row[1]));
    }
    assert_eq!(seen, db.names().into_iter().collect());
}

#[test]
fn name_lookup_uses_names_lower_cased_names_and_synonyms() {
    // Every name, lower- and upper-cased name, synonym (in all three
    // spellings) and a few unknown keys, through hasEnzyme and getEnzyme.
    let db = ProteaseDB::global();
    let keys = rows("key");
    assert_eq!(keys.len(), 125);
    for row in keys {
        let (key, has, got) = (row[0], row[1] == "1", row[2]);
        assert_eq!(db.has_enzyme(key), has, "hasEnzyme({key:?})");
        match db.get_enzyme(key) {
            Ok(enzyme) => assert_eq!(enzyme.name(), got, "getEnzyme({key:?})"),
            Err(_) => assert_eq!(got, "!ElementNotFound", "getEnzyme({key:?})"),
        }
    }
    // Protease::from_name and ProteaseDigestion::set_enzyme resolve through
    // the same index, as source ProteaseDigestion::setEnzyme does.
    assert_eq!(Protease::from_name("trypsin").unwrap(), Protease::Trypsin);
    assert!(Protease::from_name("TRYPSIN").is_err());
}

#[test]
fn regex_lookup_returns_the_enzyme_the_source_index_keeps() {
    let db = ProteaseDB::global();
    for row in rows("regex") {
        let (regex, has, got) = (row[0], row[1] == "1", row[2]);
        assert_eq!(db.has_regex(regex), has, "hasRegEx({regex:?})");
        match db.get_enzyme_by_regex(regex) {
            Ok(enzyme) => assert_eq!(enzyme.name(), got, "getEnzymeByRegEx({regex:?})"),
            Err(error) => {
                let message = got.strip_prefix("!IllegalArgument:").unwrap();
                assert!(error.to_string().contains(message), "{error} vs {message}");
            }
        }
    }
    // Two enzymes share "(?<=[MX])"; the index keeps the later one, and the
    // native all-matches lookup returns both in registration order.
    let both: Vec<_> = db
        .enzymes_by_regex("(?<=[MX])")
        .unwrap()
        .into_iter()
        .map(|e| e.name())
        .collect();
    assert_eq!(both, ["CNBr", "cyanogen-bromide"]);
    // source literal: getEnzymeByRegEx(RKP) is Arg-C
    assert_eq!(
        db.get_enzyme_by_regex("(?<=[RX])(?!P)").unwrap().name(),
        "Arg-C"
    );
}

#[test]
fn membership_of_an_enzyme_is_identity_not_equality() {
    let db = ProteaseDB::global();
    let expected: BTreeMap<&str, bool> = rows("has_ptr")
        .into_iter()
        .map(|row| (row[0], row[1] == "1"))
        .collect();
    let trypsin = db.get_enzyme("Trypsin").unwrap();
    assert_eq!(db.contains(trypsin), expected["registered"]);
    let copy = *trypsin;
    assert_eq!(copy, *trypsin);
    assert_eq!(db.contains(&copy), expected["copy"]);
    assert!(!expected["other"]);
}

/// The C++ name listings walk a set of pointers, so their order is heap
/// address order: the oracle's Linux order starts with Arg-C/P, not Trypsin.
/// The port lists in registration order (the pinned Enzymes.xml order,
/// `Protease::ALL`); the members must agree exactly.
#[test]
fn name_listings_have_the_source_members_in_registration_order() {
    let db = ProteaseDB::global();
    let registration: Vec<&str> = Protease::ALL.iter().map(|p| p.name()).collect();
    assert_eq!(db.names(), registration);
    let cpp_order = single("all_names");
    assert_eq!(cpp_order[0], "Arg-C/P");
    assert_ne!(cpp_order, registration);
    for (tag, native) in [
        ("all_names", db.names()),
        ("xtandem_names", db.xtandem_names()),
        ("comet_names", db.comet_names()),
        ("omssa_names", db.omssa_names()),
        ("msgf_names", db.msgf_names()),
    ] {
        let cpp: BTreeSet<&str> = single(tag).into_iter().collect();
        let rust: BTreeSet<&str> = native.iter().copied().collect();
        assert_eq!(rust, cpp, "{tag}");
        assert_eq!(rust.len(), native.len(), "{tag} has duplicates");
        // Registration order is kept within each engine listing.
        let positions: Vec<usize> = native
            .iter()
            .map(|n| registration.iter().position(|r| r == n).unwrap())
            .collect();
        assert!(positions.windows(2).all(|w| w[0] < w[1]), "{tag}");
    }
    // source literals
    assert!(db.xtandem_names().contains(&"no cleavage"));
    assert!(!db.omssa_names().contains(&"leukocyte elastase"));
}

#[test]
fn write_tsv_has_the_source_header_and_lines() {
    let db = ProteaseDB::global();
    let mut out = Vec::new();
    db.write_tsv(&mut out).unwrap();
    let text = String::from_utf8(out).unwrap();
    assert!(text.ends_with('\n') && ORACLE_TSV.ends_with('\n'));
    let rust: Vec<&str> = text.lines().collect();
    let cpp: Vec<&str> = ORACLE_TSV.lines().collect();
    assert_eq!(rust[0], "OpenMS_AllowedEnzymes");
    assert_eq!(rust[0], cpp[0]);
    // Same lines; order differs only by the source's address order.
    let mut rust_sorted = rust[1..].to_vec();
    let mut cpp_sorted = cpp[1..].to_vec();
    rust_sorted.sort_unstable();
    cpp_sorted.sort_unstable();
    assert_eq!(rust_sorted, cpp_sorted);
    assert_eq!(&rust[1..], db.names().as_slice());
    assert_eq!(&cpp[1..], single("all_names").as_slice());

    let dir = std::env::temp_dir().join(format!("openms-protease-db-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("enzymes.tsv");
    db.write_tsv_file(&path).unwrap();
    assert_eq!(std::fs::read_to_string(&path).unwrap(), text);
    // The source's unchecked ofstream reports nothing; the port reports it.
    assert!(
        db.write_tsv_file(dir.join("missing").join("x.tsv"))
            .is_err()
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn source_literals_of_the_class_test() {
    let db = ProteaseDB::global();
    assert!(!db.has_enzyme("Try"));
    assert!(db.has_enzyme("Trypsin"));
    assert_eq!(db.get_enzyme("Trypsin").unwrap().name(), "Trypsin");
    assert_eq!(db.get_enzyme("Clostripain").unwrap().name(), "Arg-C");
    assert!(db.get_enzyme("DOESNOTEXIST").is_err());
    assert!(!db.has_regex("(?<=[P])(?!P)"));
    assert!(db.has_regex("(?<=[RX])(?!P)"));
    assert!(db.iter().count() >= 10);
    let names = db.names();
    assert!(names.contains(&"Trypsin"));
    assert!(!names.contains(&"Tryptryp"));
    assert_eq!(db.names().len(), names.len());
}
