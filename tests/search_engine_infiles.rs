// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Search-engine input file writers: the upstream class-test sections
//! (`XTandemInfile_test.cpp` at core `bc9cc12`) and a byte-exact comparison
//! with the files the executed C++ Release build wrote
//! (`tests/data/search_engine_oracle/*_infile.out`, written on ibminode06 by
//! `driver.cpp.txt`; see `tests/search_engine_oracle.rs`).

use openms::chemistry::ModificationDefinitionsSet;
use openms::format::xtandem_infile::{ErrorUnit, MassType, XTandemInfile};
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn oracle(name: &str) -> String {
    std::fs::read_to_string(data(&format!("search_engine_oracle/{name}"))).unwrap()
}

/// The driver's `dumpFile` framing of a written file (printable ASCII,
/// tab and line feed only, which is all these writers produce here).
fn dump(text: &str) -> String {
    assert!(text
        .bytes()
        .all(|b| b == b'\n' || b == b'\t' || (32..127).contains(&b)));
    format!("file_bytes {}\n{text}\n<<EOF\n", text.len())
}

/// The class test's `xml_file` after its setter sections, in their order.
fn xtandem_class_test_infile() -> XTandemInfile {
    let mut x = XTandemInfile::new();
    x.set_fragment_mass_tolerance(13.0);
    assert_eq!(x.fragment_mass_tolerance(), 13.0);
    x.set_precursor_mass_tolerance_plus(14.0);
    assert_eq!(x.precursor_mass_tolerance_plus(), 14.0);
    x.set_precursor_mass_tolerance_minus(15.0);
    assert_eq!(x.precursor_mass_tolerance_minus(), 15.0);
    x.set_precursor_mass_error_unit(ErrorUnit::Daltons);
    assert_eq!(x.precursor_mass_error_unit(), ErrorUnit::Daltons);
    x.set_precursor_mass_error_unit(ErrorUnit::Ppm);
    assert_eq!(x.precursor_mass_error_unit(), ErrorUnit::Ppm);
    x.set_number_of_threads(16);
    assert_eq!(x.number_of_threads(), 16);
    let sets =
        ModificationDefinitionsSet::from_comma_separated("Oxidation (M)", "Carboxymethyl (C)")
            .unwrap();
    x.set_modifications(sets.clone());
    assert!(*x.modifications() == sets);
    x.set_output_filename("blubb_new_outputfilename");
    assert_eq!(x.output_filename(), "blubb_new_outputfilename");
    x.set_input_filename("blubb_new_inputfilename");
    assert_eq!(x.input_filename(), "blubb_new_inputfilename");
    x.set_taxonomy_filename("blubb_new_taxonomy_file");
    assert_eq!(x.taxonomy_filename(), "blubb_new_taxonomy_file");
    x.set_default_parameters_filename("blubb_new_default_parameters_file");
    assert_eq!(
        x.default_parameters_filename(),
        "blubb_new_default_parameters_file"
    );
    x.set_taxon("blubb_taxon");
    assert_eq!(x.taxon(), "blubb_taxon");
    x.set_max_precursor_charge(17);
    assert_eq!(x.max_precursor_charge(), 17);
    x.set_number_of_missed_cleavages(18);
    assert_eq!(x.number_of_missed_cleavages(), 18);
    x.set_max_valid_evalue(19.0);
    assert_eq!(x.max_valid_evalue(), 19.0);
    x.set_precursor_error_type(MassType::Monoisotopic);
    assert_eq!(x.precursor_error_type(), MassType::Monoisotopic);
    x.set_precursor_error_type(MassType::Average);
    assert_eq!(x.precursor_error_type(), MassType::Average);
    x.set_fragment_mass_error_unit(ErrorUnit::Daltons);
    assert_eq!(x.fragment_mass_error_unit(), ErrorUnit::Daltons);
    x.set_fragment_mass_error_unit(ErrorUnit::Ppm);
    assert_eq!(x.fragment_mass_error_unit(), ErrorUnit::Ppm);
    x
}

#[test]
fn xtandem_infile_defaults_are_the_source_constructor() {
    let x = XTandemInfile::new();
    assert_eq!(x.fragment_mass_tolerance(), 0.3);
    assert_eq!(x.precursor_mass_tolerance_plus(), 2.0);
    assert_eq!(x.precursor_mass_tolerance_minus(), 2.0);
    assert_eq!(x.fragment_mass_error_unit(), ErrorUnit::Daltons);
    assert_eq!(x.precursor_mass_error_unit(), ErrorUnit::Daltons);
    assert_eq!(x.precursor_error_type(), MassType::Monoisotopic);
    assert_eq!(x.max_precursor_charge(), 4);
    assert_eq!(x.number_of_threads(), 1);
    assert_eq!(x.cleavage_site(), "[KR]|{P}");
    assert!(!x.semi_cleavage());
    assert!(!x.allow_isotope_error());
    assert_eq!(x.number_of_missed_cleavages(), 1);
    assert_eq!(x.output_results(), "valid");
    assert_eq!(x.max_valid_evalue(), 0.01);
    assert!(x.modifications().is_empty());
}

#[test]
fn xtandem_infile_rejects_an_unknown_output_result_type() {
    let mut x = XTandemInfile::new();
    assert!(x.set_output_results("bogus").is_err());
    assert_eq!(x.output_results(), "valid");
    for result in ["all", "stochastic", "valid"] {
        x.set_output_results(result).unwrap();
        assert_eq!(x.output_results(), result);
    }
}

#[test]
fn xtandem_infile_matches_the_release_build() {
    // The class test's write section, then the driver's two extra writes
    // (ignore_member_parameters, then force_default_mods).
    let mut x = xtandem_class_test_infile();
    let sets = ModificationDefinitionsSet::from_comma_separated(
        "Oxidation (M),Dimethyl (N-term),Carboxymethyl (C)",
        "Ammonium (C-term),ICDID (C)",
    )
    .unwrap();
    x.set_modifications(sets);
    x.set_allow_isotope_error(true);
    x.set_semi_cleavage(true);
    x.set_output_results("all").unwrap();

    let mut out = String::new();
    let mut messages = Vec::new();
    for (ignore, force) in [(false, false), (true, false), (false, true)] {
        let mut written = Vec::new();
        messages.extend(x.write_to(&mut written, ignore, force).unwrap());
        out.push_str(&dump(&String::from_utf8(written).unwrap()));
    }
    let cpp = oracle("xtandem_infile.out");
    // The C++ OPENMS_LOG_INFO lines, deduplicated by its logger, follow the
    // three file dumps.
    assert!(cpp.len() > out.len());
    let (files, log) = cpp.split_at(out.len());
    assert_eq!(out, files);
    let corrected = "X! Tandem config file: Fixed modification and variable modification to origin 'C' detected. Using corrected mass of 80.0626 instead of 138.068.";
    assert_eq!(log.lines().next(), Some(corrected), "{log}");
    assert!(log.contains("occurred 2 times"));
    assert_eq!(messages.iter().filter(|m| *m == corrected).count(), 2);
}

#[test]
fn xtandem_infile_write_creates_the_file() {
    let dir = std::env::temp_dir().join(format!("xtandem_infile_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("input.xml");
    let x = XTandemInfile::new();
    x.write(&path, true, false).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    assert!(text.starts_with("<?xml version=\"1.0\"?>\n"));
    assert!(text.ends_with("</bioml>\n"));
    assert!(!text.contains("default parameters"));
    assert!(x
        .write(dir.join("missing/input.xml"), false, false)
        .is_err());
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn xtandem_infile_quick_options_cover_default_n_terminal_modifications() {
    let mut x = XTandemInfile::new();
    x.set_modifications(
        ModificationDefinitionsSet::from_comma_separated(
            "",
            "Acetyl (N-term),Gln->pyro-Glu (N-term Q),Glu->pyro-Glu (N-term E)",
        )
        .unwrap(),
    );
    let mut written = Vec::new();
    x.write_to(&mut written, false, false).unwrap();
    let text = String::from_utf8(written).unwrap();
    assert!(text.contains("label=\"protein, quick pyrolidone\">yes<"));
    assert!(text.contains("label=\"protein, quick acetyl\">yes<"));
    assert!(text.contains("label=\"residue, potential modification mass\"></note>"));

    // Forced: no quick options, every modification written explicitly.
    let mut forced = Vec::new();
    x.write_to(&mut forced, false, true).unwrap();
    let forced = String::from_utf8(forced).unwrap();
    assert!(!forced.contains("quick"));
    let line = forced
        .lines()
        .find(|l| l.contains("potential modification mass"))
        .unwrap();
    assert_eq!(line.matches('@').count(), 3);

    // Another peptide N-terminal variable modification forces them too.
    let mut x = XTandemInfile::new();
    x.set_modifications(
        ModificationDefinitionsSet::from_comma_separated("", "Acetyl (N-term),Dimethyl (N-term)")
            .unwrap(),
    );
    let mut written = Vec::new();
    x.write_to(&mut written, false, false).unwrap();
    assert!(!String::from_utf8(written).unwrap().contains("quick acetyl"));
}
