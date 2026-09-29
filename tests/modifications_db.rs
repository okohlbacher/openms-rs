// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//
// ModificationsDB (CHEMISTRY/ModificationsDB.h) against the executed source.
//
// Expected values are EXECUTED C++ output (tier 1): ../oracle/B8b-dbs/driver.cpp,
// compiled with libOpenMS's Release flags against the Linux x86_64 Release
// install openms4-release-bc9cc12-c19e494-174b576 and run twice on ibminode06,
// both runs byte-identical. The driver builds ModificationsDB from the three
// providers the Rust global registry is built from (unimod.xml,
// custom_mods.xml, XLMOD.obo, byte-identical to resources/modifications), so
// the registries can be compared record by record.
// tests/data/modifications_db_oracle.tsv is its "mods" stdout,
// tests/data/modifications_db_oracle_writetsv.tsv the file writeTSV wrote and
// tests/data/modifications_db_oracle_singleton.tsv the getInstance() case;
// tests/data/modifications_db_provenance.json records the hashes.

use openms::chemistry::{
    EmpiricalFormula, ModificationProvenance, ModificationRecord, ModificationsDB,
    ResidueModification, TermSpecificity,
};
use std::collections::{BTreeSet, HashMap};
use std::sync::{Arc, OnceLock};

const ORACLE: &str = include_str!("data/modifications_db_oracle.tsv");
const ORACLE_TSV: &str = include_str!("data/modifications_db_oracle_writetsv.tsv");
const ORACLE_SINGLETON: &str = include_str!("data/modifications_db_oracle_singleton.tsv");

fn rows(tag: &str) -> Vec<Vec<&'static str>> {
    let rows: Vec<Vec<&str>> = ORACLE
        .lines()
        .filter_map(|line| line.strip_prefix(tag)?.strip_prefix('\t'))
        .map(|rest| rest.split('\t').collect())
        .collect();
    assert!(!rows.is_empty(), "oracle has no {tag} rows");
    rows
}

fn db() -> &'static ModificationsDB {
    ModificationsDB::global()
}

/// Position of a record in the global registry, by identity.
fn index_of(record: &ResidueModification) -> usize {
    static POSITIONS: OnceLock<HashMap<usize, usize>> = OnceLock::new();
    let positions = POSITIONS.get_or_init(|| {
        db().entries()
            .iter()
            .enumerate()
            .map(|(index, entry)| (Arc::as_ptr(entry) as usize, index))
            .collect()
    });
    positions[&(std::ptr::from_ref(record) as usize)]
}

fn indices(text: &str) -> Vec<usize> {
    if text.is_empty() {
        Vec::new()
    } else {
        text.split(',').map(|i| i.parse().unwrap()).collect()
    }
}

fn optional_index(text: &str) -> Option<usize> {
    (text != "null").then(|| text.parse().unwrap())
}

fn residue(text: &str) -> Option<char> {
    text.chars().next()
}

fn term(text: &str) -> Option<TermSpecificity> {
    match text {
        "ANYWHERE" => Some(TermSpecificity::Anywhere),
        "C_TERM" => Some(TermSpecificity::CTerm),
        "N_TERM" => Some(TermSpecificity::NTerm),
        "PROTEIN_C_TERM" => Some(TermSpecificity::ProteinCTerm),
        "PROTEIN_N_TERM" => Some(TermSpecificity::ProteinNTerm),
        "NUMBER" => None,
        other => panic!("unknown specificity {other}"),
    }
}

fn exact(text: &str) -> f64 {
    text.parse().unwrap()
}

/// Keys whose C++ name set is not in load order: the set holds pointers, and
/// here the allocator placed a later record at a lower address. The source's
/// "last match" is then an allocator artefact; the port's is the last record
/// in registry order.
fn set_in_load_order(set: &[usize]) -> bool {
    set.windows(2).all(|w| w[0] < w[1])
}

#[test]
fn global_registry_is_the_source_registry_record_for_record() {
    let count: usize = rows("count")[0][0].parse().unwrap();
    assert_eq!(count, 3127);
    assert_eq!(db().len(), count);
    for row in rows("mod") {
        let i: usize = row[0].parse().unwrap();
        let m = db().entries()[i].as_ref();
        assert_eq!(m.full_id(), row[1], "{i}");
        assert_eq!(m.name(), row[2], "{i}");
        assert_eq!(m.full_name(), row[3], "{i}");
        assert_eq!(m.unimod_accession().unwrap_or_default(), row[4], "{i}");
        let record_id: i64 = row[5].parse().unwrap();
        assert_eq!(m.record_id().map_or(-1, i64::from), record_id, "{i}");
        assert_eq!(m.obo_accession().unwrap_or(""), row[6], "{i}");
        // The source stores a terminal record's origin as 'X'.
        assert_eq!(m.origin().unwrap_or('X').to_string(), row[7], "{i}");
        assert_eq!(Some(m.term_specificity()), term(row[8]), "{i}");
        assert_eq!(m.diff_mono_mass().to_bits(), exact(row[9]).to_bits(), "{i}");
        assert_eq!(
            m.diff_average_mass().to_bits(),
            exact(row[10]).to_bits(),
            "{i}"
        );
        assert_eq!(m.mono_mass().to_bits(), exact(row[11]).to_bits(), "{i}");
        assert_eq!(m.average_mass().to_bits(), exact(row[12]).to_bits(), "{i}");
        let formula = if row[13].is_empty() {
            EmpiricalFormula::default()
        } else {
            row[13].parse().unwrap()
        };
        assert_eq!(m.diff_formula(), &formula, "{i}");
        let synonyms: BTreeSet<String> = if row[14].is_empty() {
            BTreeSet::new()
        } else {
            row[14].split('|').map(str::to_owned).collect()
        };
        assert_eq!(m.synonyms(), &synonyms, "{i}");
        // isUserDefined: empty Id, non-empty FullId.
        let user_defined = m.name().is_empty() && !m.full_id().is_empty();
        assert_eq!(user_defined, row[15] == "1", "{i}");
        // Provenance 1 is CV.
        assert_eq!(row[16], "1");
        assert_eq!(m.provenance(), ModificationProvenance::Cv, "{i}");
    }
}

#[test]
fn every_lookup_key_has_the_source_index_entry() {
    let keys = rows("key");
    assert!(keys.len() > 7000);
    let mut out_of_order = 0;
    for row in keys {
        let key = row[0];
        assert_eq!(db().has(key), row[1] == "1", "has({key:?})");
        match db().find_modification_index(key) {
            Ok(index) => assert_eq!(index.to_string(), row[2], "findModificationIndex({key:?})"),
            Err(error) => {
                let message = row[2]
                    .strip_prefix("!ElementNotFound:the element '")
                    .and_then(|m| m.strip_suffix("' could not be found"))
                    .unwrap_or_else(|| panic!("{key:?}: C++ gave {}, Rust {error}", row[2]));
                assert!(error.to_string().contains(message), "{key:?}: {error}");
            }
        }
        // searchModifications without filters: the same records.
        let cpp_set = indices(row[3]);
        let rust_set: BTreeSet<usize> = db()
            .find(key, None, None)
            .into_iter()
            .map(index_of)
            .collect();
        assert_eq!(
            rust_set,
            cpp_set.iter().copied().collect::<BTreeSet<_>>(),
            "searchModifications({key:?})"
        );
        // searchModificationsFast: the last match and the multiple flag.
        let (fast, multiple) = db().search_modifications_fast(key, None, None);
        assert_eq!(multiple, row[5] == "1", "fast({key:?}) multiple");
        let fast = fast.map(index_of);
        if set_in_load_order(&cpp_set) {
            assert_eq!(fast, optional_index(row[4]), "fast({key:?})");
        } else {
            out_of_order += 1;
            assert_eq!(fast, cpp_set.iter().max().copied(), "fast({key:?})");
        }
    }
    // Pinned so a regenerated oracle cannot silently widen the exception.
    assert_eq!(out_of_order, 6);
}

#[test]
fn filtered_fast_search_follows_the_source_residue_rule() {
    let queries = rows("fast");
    assert_eq!(queries.len(), 31 * 13 * 6);
    for row in queries {
        let (name, res, spec) = (row[0], residue(row[1]), term(row[2]));
        let cpp_set = indices(row[5]);
        let (fast, multiple) = db().search_modifications_fast(name, res, spec);
        assert_eq!(multiple, row[4] == "1", "{row:?}");
        assert!(set_in_load_order(&cpp_set), "{row:?}");
        assert_eq!(fast.map(index_of), optional_index(row[3]), "{row:?}");
        // The native `find` treats the wildcard residues X, '.' and '?' as
        // literal origins, where the source matches every origin: a
        // documented gap. For every other residue it agrees with the source.
        let rust_set: BTreeSet<usize> = db()
            .find(name, res, spec)
            .into_iter()
            .map(index_of)
            .collect();
        let cpp: BTreeSet<usize> = cpp_set.iter().copied().collect();
        if matches!(res, Some('X' | '.' | '?')) {
            assert!(rust_set.is_subset(&cpp), "{row:?}");
        } else {
            assert_eq!(rust_set, cpp, "{row:?}");
        }
    }
}

#[test]
fn mass_searches_match_the_source() {
    let queries = rows("mass");
    assert_eq!(queries.len(), 38);
    for row in queries {
        let (mass, tolerance) = (exact(row[0]), exact(row[1]));
        let (res, spec) = (residue(row[2]), term(row[3]));
        let unsorted = indices(row[4]);
        let sorted = indices(row[5]);
        let best = optional_index(row[6]);
        let unsorted_ids: Vec<&str> = if row[7].is_empty() {
            Vec::new()
        } else {
            row[7].split('|').collect()
        };
        let sorted_ids: Vec<&str> = if row[8].is_empty() {
            Vec::new()
        } else {
            row[8].split('|').collect()
        };
        // Both sorted overloads: records, and their full IDs.
        let rust_sorted = db()
            .search_by_mass_sorted(mass, tolerance, res, spec)
            .unwrap();
        assert_eq!(
            rust_sorted.iter().map(|m| index_of(m)).collect::<Vec<_>>(),
            sorted,
            "{row:?}"
        );
        assert_eq!(
            rust_sorted.iter().map(|m| m.full_id()).collect::<Vec<_>>(),
            sorted_ids,
            "{row:?}"
        );
        let rust_unsorted: Vec<usize> = db()
            .search_by_mass(mass, tolerance, res, spec)
            .unwrap()
            .into_iter()
            .map(index_of)
            .collect();
        let rust_best = db()
            .best_by_mass(mass, tolerance, res, spec)
            .unwrap()
            .map(index_of);
        if matches!(res, Some('X' | '.' | '?')) {
            // The documented wildcard-residue gap of the native searches.
            assert!(
                rust_unsorted.iter().all(|i| unsorted.contains(i)),
                "{row:?}"
            );
        } else {
            assert_eq!(rust_unsorted, unsorted, "{row:?}");
            assert_eq!(
                rust_unsorted
                    .iter()
                    .map(|&i| db().entries()[i].full_id())
                    .collect::<Vec<_>>(),
                unsorted_ids,
                "{row:?}"
            );
            assert_eq!(rust_best, best, "{row:?}");
        }
    }
}

#[test]
fn search_modification_finds_the_equal_stored_record() {
    for row in rows("exact") {
        if row[0] == "changed0" {
            assert_eq!(row[1], "null");
            let first = db().entries()[0].as_ref().clone();
            let changed = first.with_absolute_masses(1.0, 1.0).unwrap();
            assert!(db().search_modification(&changed).is_none());
            continue;
        }
        let i: usize = row[0].parse().unwrap();
        let copy = db().entries()[i].as_ref().clone();
        let found = db().search_modification(&copy).map(index_of);
        assert_eq!(found, optional_index(row[1]), "{i}");
    }
}

#[test]
fn search_modifications_are_sorted_as_the_source_sorts_them() {
    let count: usize = rows("search_mods")[0][0].parse().unwrap();
    let cpp: Vec<&str> = rows("search_mod").into_iter().map(|row| row[0]).collect();
    assert_eq!(cpp.len(), count);
    assert_eq!(db().all_search_modifications(), cpp);
    // source literals (ModificationsDB_test.cpp)
    for present in [
        "Phospho (S)",
        "Sulfo (S)",
        "NIC (N-term)",
        "Dehydrated (N-term C)",
    ] {
        assert!(cpp.contains(&present));
    }
    assert!(!cpp.contains(&"Phospho"));
}

#[test]
fn write_tsv_reproduces_the_source_file() {
    let mut out = Vec::new();
    db().write_tsv(&mut out).unwrap();
    assert_eq!(String::from_utf8(out).unwrap(), ORACLE_TSV);
}

#[test]
fn add_modification_keeps_the_source_duplicate_and_key_rules() {
    let adds = rows("add");
    let row = |name: &str| adds.iter().find(|r| r[0] == name).unwrap().clone();
    let mut registry = db().clone();
    assert_eq!(registry.has(""), row("has_empty_before")[1] == "1");

    // A record whose full ID is registered returns the registered record.
    let duplicate = registry
        .add_modification(db().entries()[0].as_ref().clone())
        .unwrap();
    let dup = row("dup");
    assert!(Arc::ptr_eq(
        &duplicate,
        &registry.entries()[dup[1].parse::<usize>().unwrap()]
    ));
    assert_eq!(registry.len().to_string(), dup[2]);

    // A user-defined record known only by its full ID.
    let new = row("new");
    let record = ResidueModification::from_record(ModificationRecord {
        full_id: "Phospho (A)".into(),
        ..Default::default()
    })
    .unwrap();
    let added = registry.add_modification(record).unwrap();
    assert!(Arc::ptr_eq(&added, registry.entries().last().unwrap()) == (new[1] == "1"));
    assert_eq!(registry.len().to_string(), new[2]);
    assert_eq!(registry.has("Phospho (A)"), new[3] == "1");
    // The source indexes the empty Id, FullName and UniMod accession too.
    assert_eq!(registry.has(""), new[4] == "1");
    assert_eq!(
        registry
            .find_modification_index("Phospho (A)")
            .unwrap()
            .to_string(),
        new[5]
    );
    assert_eq!(added.name(), new[6]);
    assert_eq!(added.full_name(), new[7]);
    assert_eq!(added.unimod_accession().unwrap_or_default(), new[8]);
    assert_eq!(added.origin().unwrap_or('X').to_string(), new[9]);
    assert_eq!(Some(added.term_specificity()), term(new[10]));

    // A named record: indexed under its short ID and derived full ID.
    let second = row("second");
    let record = ResidueModification::from_record(ModificationRecord {
        name: "B8bTest".into(),
        origin: Some('K'),
        diff_mono_mass: 1.5,
        ..Default::default()
    })
    .unwrap();
    let added = registry.add_modification(record.clone()).unwrap();
    assert_eq!(added.full_id(), second[1]);
    assert_eq!(registry.len().to_string(), second[2]);
    assert_eq!(registry.has("B8bTest"), second[3] == "1");
    assert_eq!(registry.has("B8bTest (K)"), second[4] == "1");
    let again = registry.add_modification(record).unwrap();
    let second_again = row("second_again");
    assert_eq!(Arc::ptr_eq(&again, &added), second_again[1] == "1");
    assert_eq!(registry.len().to_string(), second_again[2]);

    // The global registry itself is untouched.
    assert_eq!(db().len(), 3127);
    assert!(!db().has("Phospho (A)"));
}

#[test]
fn instantiation_is_reported_once_the_global_registry_exists() {
    // The source singleton additionally loads PSI-MOD.obo, 3,610 records; the
    // port does not bundle PSI-MOD (docs/MODIFICATION_REGISTRY_SUPPORT.md).
    let lines: Vec<&str> = ORACLE_SINGLETON.lines().collect();
    assert_eq!(lines[0], "instantiated_before\t0");
    assert_eq!(lines[1], "instantiated_after\t1");
    assert_eq!(lines[2], "count\t3610");
    let _ = ModificationsDB::global();
    assert!(ModificationsDB::is_instantiated());
}

#[test]
fn unimod_prefix_retry_follows_the_six_letter_source_rule() {
    // Source: size > 6 and lower-cased first six letters "unimod" retries as
    // "UniMod" + rest.
    for name in ["unimod:21", "UNIMOD:21", "uNiMoD:21"] {
        let (found, _) = db().search_modifications_fast(name, Some('S'), None);
        assert_eq!(found.unwrap().full_id(), "Phospho (S)", "{name}");
    }
    assert!(
        db().search_modifications_fast("unimod", None, None)
            .0
            .is_none()
    );
    assert!(
        db().search_modifications_fast("unimod21", None, None)
            .0
            .is_none()
    );
}

#[test]
fn checked_boundaries() {
    assert!(
        db().search_by_mass_sorted(f64::NAN, 0.1, None, None)
            .is_err()
    );
    assert!(db().search_by_mass_sorted(1.0, -0.1, None, None).is_err());
    assert!(
        db().search_by_mass_sorted(1.0, f64::INFINITY, None, None)
            .is_err()
    );
    assert!(db().find_modification_index("").is_err());
    assert!(!db().has(""));
    let dir = std::env::temp_dir().join(format!("openms-mods-db-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("mods.tsv");
    db().write_tsv_file(&path).unwrap();
    assert_eq!(std::fs::read_to_string(&path).unwrap(), ORACLE_TSV);
    assert!(
        db().write_tsv_file(dir.join("missing").join("x.tsv"))
            .is_err()
    );
    std::fs::remove_dir_all(&dir).unwrap();
}
