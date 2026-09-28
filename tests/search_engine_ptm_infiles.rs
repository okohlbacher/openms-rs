// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! PTMXML tables and the search-engine input files that read them
//! (Inspect, Sequest): the upstream class-test sections
//! (`PTMXMLFile_test.cpp`, `InspectInfile_test.cpp`, `SequestInfile_test.cpp`
//! at core `bc9cc12`, fixtures copied from its `src/tests/class_tests/openms/data/`)
//! and a byte-exact comparison with the files the executed C++ Release build
//! wrote (`tests/data/search_engine_oracle/{inspect,sequest}_infile.out`,
//! ibminode06, `driver.cpp.txt`).
#![cfg(feature = "idxml")]

use openms::format::inspect_infile::InspectInfile;
use openms::format::ptm_xml_file::PTMXMLFile;
use openms::Error;
use std::collections::BTreeMap;
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn oracle(name: &str) -> String {
    std::fs::read_to_string(data(&format!("search_engine_oracle/{name}"))).unwrap()
}

fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("f3b_{tag}_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The driver's `dumpFile` framing of a written file.
fn dump(text: &str) -> String {
    assert!(text
        .bytes()
        .all(|b| b == b'\n' || b == b'\t' || (32..127).contains(&b)));
    format!("file_bytes {}\n{text}\n<<EOF\n", text.len())
}

fn parse_message(error: Error) -> String {
    match error {
        Error::Parse { message, .. } => message,
        other => panic!("expected a parse error, got {other:?}"),
    }
}

// ---------------------------------------------------------------- PTMXMLFile

#[test]
fn ptmxml_load() {
    let table = PTMXMLFile::new().load(data("PTMs.xml")).unwrap();
    assert_eq!(table["TEST"].0, "N2O2-CH3");
    assert_eq!(table["TEST"].1, "KLR");
}

#[test]
fn ptmxml_store_round_trips() {
    let file = PTMXMLFile::new();
    let table = file.load(data("PTMs.xml")).unwrap();
    let dir = temp_dir("ptmxml");
    let path = dir.join("PTMs.xml");
    file.store(&path, &table).unwrap();
    let back = file.load(&path).unwrap();
    assert_eq!(back["TEST"].0, "N2O2-CH3");
    assert_eq!(back["TEST"].1, "KLR");
    // The source layout, byte for byte: it is the fixture's.
    assert_eq!(
        std::fs::read_to_string(&path).unwrap(),
        std::fs::read_to_string(data("PTMs.xml")).unwrap()
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn ptmxml_missing_file_is_an_io_error() {
    let error = PTMXMLFile::new()
        .load(data("does_not_exist_PTMs.xml"))
        .unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error:?}");
}

// ------------------------------------------------------------- InspectInfile

const LINE: &str = "10.3+,KRLNH,fix:+16,C:16-,cterm:-16,nterm";

fn expected_inspect_mods() -> BTreeMap<String, Vec<String>> {
    [
        ("PTM_0", ["KRLNH", "10.3", "FIX"]),
        ("PTM_1", ["C", "16", "OPT"]),
        ("PTM_2", ["CTERM", "-16", "OPT"]),
        ("PTM_3", ["NTERM", "-16", "OPT"]),
    ]
    .into_iter()
    .map(|(k, v)| (k.to_owned(), v.iter().map(|s| (*s).to_owned()).collect()))
    .collect()
}

#[test]
fn inspect_infile_copy_assignment_and_equality() {
    let mut first = InspectInfile::new();
    first.set_spectra("dummy");
    let second = first.clone();
    let mut third = InspectInfile::new();
    third.set_spectra("dummy");
    first = InspectInfile::new();
    assert!(second == third);
    assert!(first == InspectInfile::new());

    // The source equality ignores the database.
    let mut with_db = third.clone();
    with_db.set_db("other");
    assert!(with_db == third);
    assert_eq!(with_db.db(), "other");
}

#[test]
fn inspect_infile_handle_ptms_errors() {
    let ptms = data("Inspect_PTMs.xml");
    let mut file = InspectInfile::new();

    let error = file
        .handle_ptms("Phosphorylation", data("a"), true)
        .unwrap_err();
    match error {
        Error::Io(e) => assert!(e.to_string().contains("could not be found"), "{e}"),
        other => panic!("{other:?}"),
    }
    assert_eq!(
        parse_message(file.handle_ptms("2H20,KRLNH,fix", &ptms, true).unwrap_err()),
        "There's something wrong with this modification. Aborting! in: 2H20,KRLNH,fix"
    );
    assert_eq!(
        parse_message(file.handle_ptms("10.3+", &ptms, true).unwrap_err()),
        "No residues for modification given. Aborting! in: 10.3+"
    );
    assert_eq!(
        parse_message(
            file.handle_ptms("10.3+,KRLNH,stat,PTM_0", &ptms, true)
                .unwrap_err()
        ),
        "There's something wrong with the type of this modification. Aborting! in: 10.3+,KRLNH,stat,PTM_0"
    );
    assert_eq!(
        parse_message(
            file.handle_ptms("Phosphorylation:Phosphorylation", &ptms, true)
                .unwrap_err()
        ),
        "There's already a modification with this name. Aborting! in: Phosphorylation"
    );
    assert!(file.modifications().is_empty());
}

#[test]
fn inspect_infile_handle_ptms() {
    let ptms = data("Inspect_PTMs.xml");
    let mut file = InspectInfile::new();
    // average masses
    file.handle_ptms(LINE, &ptms, false).unwrap();
    assert_eq!(*file.modifications(), expected_inspect_mods());
    // monoisotopic masses
    file.handle_ptms(LINE, &ptms, true).unwrap();
    assert_eq!(*file.modifications(), expected_inspect_mods());
}

#[test]
fn inspect_infile_named_modification_from_the_table() {
    // Independently derived: HO3P = H 1.00782503207 + 3 O 15.99491461956 +
    // P 30.97376163 = 79.96633 (mono, cut after five decimals).
    let mut file = InspectInfile::new();
    file.handle_ptms("Phosphorylation", data("Inspect_PTMs.xml"), true)
        .unwrap();
    let entry = &file.modifications()["Phosphorylation"];
    assert_eq!(entry[0], "STYDHCR");
    assert_eq!(entry[1], "79.96633");
    assert_eq!(entry[2], "OPT");
}

fn class_test_inspect_infile() -> InspectInfile {
    let mut file = InspectInfile::new();
    file.set_spectra("dummy4712");
    assert_eq!(file.spectra(), "dummy4712");
    file.set_db("dummy4711");
    assert_eq!(file.db(), "dummy4711");
    file.set_enzyme("Trypsin");
    assert_eq!(file.enzyme(), "Trypsin");
    let ptms = data("Inspect_PTMs.xml");
    file.handle_ptms(LINE, &ptms, false).unwrap();
    file.handle_ptms(LINE, &ptms, true).unwrap();
    file.handle_ptms(LINE, &ptms, false).unwrap();
    assert_eq!(*file.modifications(), expected_inspect_mods());
    file.set_modifications_per_peptide(2);
    assert_eq!(file.modifications_per_peptide(), 2);
    file.set_blind(1);
    assert_eq!(file.blind(), 1);
    file.set_max_ptm_size(250.0);
    assert_eq!(file.max_ptm_size(), 250.0);
    file.set_precursor_mass_tolerance(1.3);
    assert_eq!(file.precursor_mass_tolerance(), 1.3f32);
    file.set_peak_mass_tolerance(0.3);
    assert_eq!(file.peak_mass_tolerance(), 0.3f32);
    file.set_multicharge(1);
    assert_eq!(file.multicharge(), 1);
    file.set_instrument("ESI-ION-TRAP");
    assert_eq!(file.instrument(), "ESI-ION-TRAP");
    file.set_tag_count(1);
    assert_eq!(file.tag_count(), 1);
    file
}

#[test]
fn inspect_infile_store_matches_the_class_test_template() {
    let mut file = class_test_inspect_infile();
    let dir = temp_dir("inspect");
    let path = dir.join("inspect.tsv");
    file.store(&path).unwrap();
    assert_eq!(
        std::fs::read(&path).unwrap(),
        std::fs::read(data("InspectInfile_test_template1.txt")).unwrap()
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn inspect_infile_matches_the_release_build() {
    let mut file = class_test_inspect_infile();
    let dir = temp_dir("inspect_oracle");
    let mut out = String::new();
    // An unknown extension passes the source's extension check.
    file.store(dir.join("inspect_infile.tmp")).unwrap();
    out.push_str("store .tmp succeeded\n");
    let path = dir.join("inspect_infile.tsv");
    file.store(&path).unwrap();
    out.push_str(&dump(&std::fs::read_to_string(&path).unwrap()));
    for (name, entry) in file.modifications() {
        out.push_str(&format!("mod {name}"));
        for value in entry {
            out.push_str(&format!(" [{value}]"));
        }
        out.push('\n');
    }
    assert_eq!(out, oracle("inspect_infile.out"));
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn inspect_infile_store_refuses_a_known_other_extension() {
    let mut file = InspectInfile::new();
    let dir = temp_dir("inspect_ext");
    assert!(matches!(
        file.store(dir.join("inspect.mzML")),
        Err(Error::InvalidValue(_))
    ));
    assert!(!dir.join("inspect.mzML").exists());
    std::fs::remove_dir_all(&dir).unwrap();
}
