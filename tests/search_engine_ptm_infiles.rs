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

use openms::Error;
use openms::format::inspect_infile::InspectInfile;
use openms::format::ptm_xml_file::PTMXMLFile;
use openms::format::sequest_infile::SequestInfile;
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
    assert!(
        text.bytes()
            .all(|b| b == b'\n' || b == b'\t' || (32..127).contains(&b))
    );
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

// ------------------------------------------------------------- SequestInfile

const SEQUEST_LINE: &str =
    "10.3+,KRLNH,fix:+16,C:16-,cterm,opt:-16,nterm,fix:17-,cterm_prot:-17,nterm_prot,fix";

fn expected_sequest_mods() -> BTreeMap<String, Vec<String>> {
    [
        ("PTM_0", ["KRLNH", "10.3", "FIX"]),
        ("PTM_1", ["C", "16", "OPT"]),
        ("PTM_2", ["CTERM", "-16", "OPT"]),
        ("PTM_3", ["NTERM", "-16", "FIX"]),
        ("PTM_4", ["CTERM_PROT", "-17", "OPT"]),
        ("PTM_5", ["NTERM_PROT", "-17", "FIX"]),
    ]
    .into_iter()
    .map(|(k, v)| (k.to_owned(), v.iter().map(|s| (*s).to_owned()).collect()))
    .collect()
}

const STANDARD_ENZYME_TABLE: &str = "[SEQUEST_ENZYME_INFO]
0.  AspN                    0     D             -
1.  AspN_DE                 0     DE            -
2.  Chymotrypsin            1     FWYL          -
3.  Chymotrypsin_WYF        1     FWY           -
4.  Clostripain             1     R             -
5.  Cyanogen_Bromide        1     M             -
6.  Elastase                1     ALIV          P
7.  Elastase/Tryp/Chymo     1     ALIVKRWFY     P
8.  GluC                    1     E             -
9.  GluC_ED                 1     ED            -
10.  IodosoBenzoate          1     W             -
11.  LysC                    1     K             -
12.  No_Enzyme               0     -             -
13.  Proline_Endopept        1     P             -
14.  Trypsin                 1     KRLNH         -
15.  Trypsin/Chymo           1     KRLFWYN       -
16.  Trypsin_Strict          1     KR            -
";

#[test]
fn sequest_infile_copy_and_equality() {
    let mut first = SequestInfile::new();
    first.set_database("dummy");
    let second = first.clone();
    assert!(first == second);
    assert!(second != SequestInfile::new());
    let mut third = SequestInfile::new();
    third.set_database("dummy");
    assert!(first == third);
    third.set_database("another dummy");
    assert!(first != third);
}

#[test]
fn sequest_infile_handle_ptms_errors() {
    let ptms = data("Sequest_PTMs.xml");
    let mut file = SequestInfile::new();
    assert!(matches!(
        file.handle_ptms("Phosphorylation", data("a"), true),
        Err(Error::Io(_))
    ));
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
}

/// The class test's `file` after its sections, in their order; `add_z`
/// adds the class test's `Z_TestEnzyme` (the oracle driver does not).
fn class_test_sequest_infile(add_z: bool) -> SequestInfile {
    let mut file = SequestInfile::new();
    assert_eq!(file.enzyme_info_as_string(), STANDARD_ENZYME_TABLE);
    if add_z {
        file.add_enzyme_info("Z_TestEnzyme", "1", "RMW", "-");
        assert_eq!(
            file.enzyme_info_as_string(),
            format!("{STANDARD_ENZYME_TABLE}17.  Z_TestEnzyme            1     RMW           -\n")
        );
    }
    let ptms = data("Sequest_PTMs.xml");
    file.handle_ptms(SEQUEST_LINE, &ptms, false).unwrap();
    assert_eq!(*file.modifications(), expected_sequest_mods());
    file.handle_ptms(SEQUEST_LINE, &ptms, true).unwrap();
    assert_eq!(*file.modifications(), expected_sequest_mods());
    file.handle_ptms(SEQUEST_LINE, &ptms, false).unwrap();
    assert_eq!(*file.modifications(), expected_sequest_mods());

    let database = r"\\bude\langwisc\sequest_test\Analysis.mzXML.fasta";
    file.set_database(database);
    assert_eq!(file.database(), database);
    file.set_neutral_losses_for_ions("0 1 1");
    assert_eq!(file.neutral_losses_for_ions(), "0 1 1");
    file.set_ion_series_weights("0 1.0 0 0 0 0 0 1.0 0");
    assert_eq!(file.ion_series_weights(), "0 1.0 0 0 0 0 0 1.0 0");
    file.set_partial_sequence("SEQVEST TEST");
    assert_eq!(file.partial_sequence(), "SEQVEST TEST");
    file.set_sequence_header_filter("homo~sapiens !mus musculus");
    assert_eq!(file.sequence_header_filter(), "homo~sapiens !mus musculus");
    file.set_precursor_mass_tolerance(1.3);
    assert_eq!(file.precursor_mass_tolerance(), 1.3f32);
    file.set_peak_mass_tolerance(0.3);
    assert_eq!(file.peak_mass_tolerance(), 0.3f32);
    file.set_match_peak_tolerance(1.2);
    assert_eq!(file.match_peak_tolerance(), 1.2f32);
    file.set_ion_cutoff_percentage(0.3);
    assert_eq!(file.ion_cutoff_percentage(), 0.3f32);
    file.set_protein_mass_filter("30.2 0");
    assert_eq!(file.protein_mass_filter(), "30.2 0");
    file.set_peptide_mass_unit(0);
    assert_eq!(file.peptide_mass_unit(), 0);
    file.set_output_lines(10);
    assert_eq!(file.output_lines(), 10);
    // The class test expects 18: it added Z_TestEnzyme to the 17 standard
    // enzymes; the oracle driver did not and printed 17.
    let expected_bad = if add_z { 18 } else { 17 };
    assert_eq!(file.set_enzyme("i_dont_exist_enzyme"), expected_bad);
    assert!(file.enzyme_name().is_err());
    assert_eq!(file.set_enzyme("Trypsin"), 0);
    assert_eq!(file.enzyme_number(), 14);
    assert_eq!(file.enzyme_name().unwrap(), "Trypsin");
    file.set_max_aa_per_mod_per_peptide(4);
    assert_eq!(file.max_aa_per_mod_per_peptide(), 4);
    file.set_max_mods_per_peptide(3);
    assert_eq!(file.max_mods_per_peptide(), 3);
    file.set_nucleotide_reading_frame(0);
    assert_eq!(file.nucleotide_reading_frame(), 0);
    file.set_max_internal_cleavage_sites(2);
    assert_eq!(file.max_internal_cleavage_sites(), 2);
    file.set_match_peak_count(5);
    assert_eq!(file.match_peak_count(), 5);
    file.set_match_peak_allowed_error(4);
    assert_eq!(file.match_peak_allowed_error(), 4);
    file.set_show_fragment_ions(true);
    assert!(file.show_fragment_ions());
    file.set_print_duplicate_references(true);
    assert!(file.print_duplicate_references());
    file.set_remove_precursor_near_peaks(true);
    assert!(file.remove_precursor_near_peaks());
    file.set_mass_type_parent(true);
    assert!(file.mass_type_parent());
    file.set_mass_type_fragment(true);
    assert!(file.mass_type_fragment());
    file.set_normalize_xcorr(true);
    assert!(file.normalize_xcorr());
    file.set_residues_in_upper_case(true);
    assert!(file.residues_in_upper_case());
    file
}

#[test]
fn sequest_infile_store_matches_the_class_test_template() {
    let file = class_test_sequest_infile(true);
    let dir = temp_dir("sequest");
    let path = dir.join("sequest.params");
    file.store(&path).unwrap();
    assert_eq!(
        std::fs::read(&path).unwrap(),
        std::fs::read(data("SequestInfile_test_template1.txt")).unwrap()
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn sequest_infile_matches_the_release_build() {
    let file = class_test_sequest_infile(false);
    let mut out = String::from("enzyme_bad 17\nenzyme_trypsin 0\n");
    out.push_str(&dump(&file.to_params_string().unwrap()));
    for (name, entry) in file.modifications() {
        out.push_str(&format!("mod {name}"));
        for value in entry {
            out.push_str(&format!(" [{value}]"));
        }
        out.push('\n');
    }
    out.push_str(&format!(
        "enzyme_info\n{}<<EOF\n",
        file.enzyme_info_as_string()
    ));
    assert_eq!(out, oracle("sequest_infile.out"));
}

#[test]
fn sequest_infile_add_enzyme_info_deduplicates_repeated_residues() {
    let mut file = SequestInfile::new();
    file.add_enzyme_info("Test", "1", "RKR", "P");
    assert_eq!(file.enzyme_info()["Test"], ["1", "KR", "P"]);
    assert_eq!(file.enzyme_name().unwrap(), "Test");
    file.add_enzyme_info("Test", "1", "RK", "P");
    assert_eq!(file.enzyme_info()["Test"], ["1", "RK", "P"]);
}

#[test]
fn sequest_infile_groups_differential_modifications_by_mass() {
    // Independently derived from the source's store(): optional residue
    // masses are summed per residue, then residues with equal sums share one
    // "mass residues" entry, in ascending mass order, masses written by
    // StringUtils::toStr(float) ("16.0").
    let mut file = SequestInfile::new();
    file.handle_ptms("16,MW:16,C:-1,K", data("Sequest_PTMs.xml"), true)
        .unwrap();
    let text = file.to_params_string().unwrap();
    assert!(
        text.contains("\ndiff_search_options = -1.0 K 16.0 CMW\n"),
        "{text}"
    );
    assert!(text.contains("\nterm_diff_search_options = 0 0\n"));
    assert!(
        SequestInfile::new()
            .to_params_string()
            .unwrap()
            .contains("\ndiff_search_options = 0 X\n")
    );
}
