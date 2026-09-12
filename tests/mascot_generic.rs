// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//
//! Tests for `src/format/mascot_generic.rs`, the port of
//! `FORMAT/MascotGenericFile.h`.
//!
//! Every `START_SECTION` of the upstream `MascotGenericFile_test.cpp` is
//! reproduced; the mapping and the asserted literals are recorded in
//! `docs/MASCOT_GENERIC_SUPPORT.md` and
//! `tests/data/mascot_generic_provenance.json`.

use openms::concept::constants::user_param::{MSM_METABOLITE_NAME, MSM_SMILES_STRING};
use openms::format::mascot_generic::{
    self as mgf, CarryOver, ExperimentCollector, MascotGenericFile, ReadOptions,
};
use openms::interfaces::MSDataConsumer;
use openms::kernel::SpectrumType;
use openms::metadata::{ExperimentalSettings, SourceFile};
use openms::param::ParamValue;
use openms::{Error, MSExperiment, MSSpectrum, Peak1D, Precursor};
use std::ops::ControlFlow;

const INFILE: &str = include_str!("data/MascotInfile_test.mascot_in");
const GNPS: &str = include_str!("data/MascotGenericFile_GNPS.mgf");

fn writer() -> MascotGenericFile {
    MascotGenericFile::new().expect("default Mascot generic writer")
}

fn strings(values: &[&str]) -> ParamValue {
    ParamValue::StringList(values.iter().map(|v| (*v).to_owned()).collect())
}

fn store_to_string(
    file: &mut MascotGenericFile,
    filename: &str,
    experiment: &MSExperiment,
    compact: bool,
) -> String {
    let mut bytes = Vec::new();
    file.store_to(&mut bytes, filename, experiment, compact)
        .expect("store to a stream");
    String::from_utf8(bytes).expect("MGF output is UTF-8")
}

// ---------------------------------------------------------------------------
// START_SECTION(MascotGenericFile())
// ---------------------------------------------------------------------------

#[test]
fn default_construction_registers_the_source_parameters() {
    let file = writer();
    // The upstream section only asserts that construction yields a non-null
    // pointer; the Rust equivalent is that the source defaults are present.
    assert_eq!(
        file.parameters().value("database").unwrap(),
        &ParamValue::String("MSDB".into())
    );
    assert_eq!(
        file.parameters().value("search_type").unwrap(),
        &ParamValue::String("MIS".into())
    );
    assert_eq!(
        file.parameters().value("missed_cleavages").unwrap(),
        &ParamValue::Integer(1)
    );
    assert_eq!(
        file.parameters().value("precursor_mass_tolerance").unwrap(),
        &ParamValue::Float(3.0)
    );
    assert_eq!(
        file.parameters().value("internal:boundary").unwrap(),
        &ParamValue::String("GZWgAaYKjHFeUaLOLEIOMq".into())
    );
    assert!(!file.store_compact());
    // updateMembers_ expanded the five default specificity groups into ten
    // per-residue keys.
    let groups = file.special_modification_groups();
    assert_eq!(groups.len(), 10);
    assert_eq!(groups["Phospho (S)"], "Phospho (ST)");
    assert_eq!(groups["Deamidated (Q)"], "Deamidated (NQ)");
    assert_eq!(groups["Cation:Na (E)"], "Cation:Na (DE)");
}

// ---------------------------------------------------------------------------
// START_SECTION(virtual ~MascotGenericFile())
// ---------------------------------------------------------------------------

#[test]
fn dropping_the_writer_releases_it() {
    // The upstream destructor section has no assertion: it only proves `delete`
    // on a heap instance is well formed. The Rust analogue is that the value is
    // owned, droppable and clonable without shared state.
    let file = writer();
    let copy = file.clone();
    drop(file);
    assert_eq!(copy.special_modification_groups().len(), 10);
    drop(copy);
}

// ---------------------------------------------------------------------------
// START_SECTION(template <typename MapType> void load(const std::string&, MapType&))
// ---------------------------------------------------------------------------

#[test]
fn load_upstream_infile_fixture() {
    let experiment = mgf::read(INFILE.as_bytes()).unwrap();
    assert_eq!(experiment.spectra.len(), 1);
    assert_eq!(experiment.spectra[0].peaks.len(), 9);
    let spectrum = &experiment.spectra[0];
    assert_eq!(spectrum.ms_level, 2);
    assert_eq!(spectrum.spectrum_type, SpectrumType::Centroid);
    assert_eq!(spectrum.native_id, "index=0");
    assert_eq!(spectrum.rt, 25.379);
    assert_eq!(spectrum.precursors[0].mz, 1998.0);
    assert_eq!(spectrum.precursors[0].charge, 0);
    // The title gains the native ID so repeated titles stay distinguishable.
    assert_eq!(
        spectrum.metadata["TITLE"].as_str().unwrap(),
        "Testtitle_index=0"
    );
    assert_eq!(spectrum.peaks[0], Peak1D::new(1.0, 1.0));
    assert_eq!(spectrum.peaks[8], Peak1D::new(9.0, 81.0));
    // The fixture's fifth and seventh peak lines carry a trailing comment in a
    // third whitespace field, which the source parses as the optional per-peak
    // charge and discards.
    assert_eq!(spectrum.peaks[4], Peak1D::new(5.0, 25.0));
    assert_eq!(spectrum.peaks[6], Peak1D::new(7.0, 49.0));
    // The same content through the path-taking entry point.
    let via_path = mgf::load("tests/data/MascotInfile_test.mascot_in").unwrap();
    assert_eq!(via_path, experiment);
    let mut destination = MSExperiment::default();
    mgf::read_into(INFILE.as_bytes(), &mut destination).unwrap();
    assert_eq!(destination, experiment);
}

// ---------------------------------------------------------------------------
// START_SECTION(void store(std::ostream&, const std::string&, const PeakMap&, bool))
// ---------------------------------------------------------------------------

#[test]
fn store_to_stream_reproduces_the_upstream_expectations() {
    let mut experiment = mgf::read(INFILE.as_bytes()).unwrap();
    let mut file = writer();
    let mut parameters = file.parameters().clone();
    parameters
        .set_value(
            "fixed_modifications",
            strings(&["Carbamidomethyl (C)", "Phospho (S)"]),
            "",
            &[],
        )
        .unwrap();
    parameters
        .set_value(
            "variable_modifications",
            strings(&["Oxidation (M)", "Deamidated (N)", "Deamidated (Q)"]),
            "",
            &[],
        )
        .unwrap();
    file.set_parameters(&parameters).unwrap();

    let text = store_to_string(&mut file, "test", &experiment, false);
    for expected in [
        "BEGIN IONS\nTITLE=Testtitle_index=0\nPEPMASS=1998.0\nRTINSECONDS=25.379000000000001\nSCANS=0",
        "1.0 1.0\n2.0 4.0\n3.0 9.0\n4.0 16.0\n5.0 25.0\n6.0 36.0\n7.0 49.0\n8.0 64.0\n9.0 81.0\nEND IONS\n",
        "MODS=Carbamidomethyl (C)\n",
        "MODS=Phospho (ST)\n",
        "IT_MODS=Deamidated (NQ)",
        "IT_MODS=Oxidation (M)",
    ] {
        assert!(text.contains(expected), "missing {expected:?} in\n{text}");
    }

    // Without a stored TITLE the writer composes one from m/z, RT, native ID
    // and file name.
    experiment.spectra[0].metadata.remove("TITLE");
    let text = store_to_string(&mut file, "test", &experiment, false);
    for expected in [
        "BEGIN IONS\nTITLE=1998.0_25.379000000000001_index=0_test\nPEPMASS=1998.0\nRTINSECONDS=25.379000000000001\nSCANS=0",
        "1.0 1.0\n2.0 4.0\n3.0 9.0\n4.0 16.0\n5.0 25.0\n6.0 36.0\n7.0 49.0\n8.0 64.0\n9.0 81.0\nEND IONS\n",
        "MODS=Carbamidomethyl (C)\n",
        "MODS=Phospho (ST)\n",
        "IT_MODS=Deamidated (NQ)",
        "IT_MODS=Oxidation (M)",
    ] {
        assert!(text.contains(expected), "missing {expected:?} in\n{text}");
    }

    // Reset the parameters and check the compact form.
    let mut file = writer();
    let spectrum = MSSpectrum {
        native_id: "index=250".into(),
        ms_level: 2,
        rt: 234.567_890_1,
        precursors: vec![Precursor::new(901.234_567_8, 0)],
        // Intensity zero is omitted from the compact form.
        peaks: vec![
            Peak1D::new(567.890_123_4, 0.0),
            // The upstream literal is 2345.678901, whose nearest f32 is 2345.679.
            Peak1D::new(890.123_456_7, 2_345.679),
        ],
        ..Default::default()
    };
    let compact = MSExperiment {
        spectra: vec![spectrum],
        ..Default::default()
    };
    let text = store_to_string(&mut file, "test", &compact, true);
    assert!(
        text.contains(
            "BEGIN IONS\nTITLE=901.23457_234.568_index=250_test\nPEPMASS=901.23457\nRTINSECONDS=234.568\nSCANS=250\n890.12346 2345.679\nEND IONS"
        ),
        "compact output was\n{text}"
    );
    assert!(file.store_compact());
}

/// A compact spectrum that already carries a `TITLE` gets *significant* digits,
/// not fixed decimals, until something sets the stream's `fixed` flag.
///
/// `MascotGenericFile.cpp` streams `fixed` only in the branch that generates a
/// title, and again on every compact peak line; with a `TITLE` meta value
/// present the stream is still in its default float format, so
/// `setprecision(5)`/`setprecision(3)` mean five and three significant digits:
/// `901.23` and `235`. Because the flag lives in the ostream it stays set for
/// the rest of the file, so the *second* such spectrum is written with fixed
/// decimals. Found by a second-model review; the port used to write
/// `901.23457`/`234.568` for both, which is more precision than the source
/// keeps.
#[test]
fn compact_output_with_a_stored_title_uses_significant_digits() {
    let titled = |index: u32| MSSpectrum {
        native_id: format!("index={index}"),
        ms_level: 2,
        rt: 234.567_890_1,
        precursors: vec![Precursor::new(901.234_567_8, 0)],
        peaks: vec![Peak1D::new(890.123_456_7, 2_345.679)],
        metadata: [("TITLE".to_owned(), openms::metadata::MetaValue::from("x"))]
            .into_iter()
            .collect(),
        ..Default::default()
    };
    let mut file = writer();
    let experiment = MSExperiment {
        spectra: vec![titled(0)],
        ..Default::default()
    };
    let text = store_to_string(&mut file, "test", &experiment, true);
    assert!(
        text.contains("TITLE=x\nPEPMASS=901.23\nRTINSECONDS=235\n"),
        "compact output was\n{text}"
    );
    // The peak line is fixed in both cases.
    assert!(text.contains("\n890.12346 2345.679\n"), "{text}");
    // The flag the peak line set is sticky, so the next spectrum is fixed.
    let experiment = MSExperiment {
        spectra: vec![titled(0), titled(1)],
        ..Default::default()
    };
    let text = store_to_string(&mut file, "test", &experiment, true);
    assert!(
        text.contains("PEPMASS=901.23\nRTINSECONDS=235\n"),
        "compact output was\n{text}"
    );
    assert!(
        text.contains("PEPMASS=901.23457\nRTINSECONDS=234.568\n"),
        "compact output was\n{text}"
    );
    // Without a stored title nothing changes: the generated title sets `fixed`
    // before the precursor m/z is written.
    let mut plain = titled(0);
    plain.metadata.clear();
    let experiment = MSExperiment {
        spectra: vec![plain],
        ..Default::default()
    };
    let text = store_to_string(&mut file, "test", &experiment, true);
    assert!(
        text.contains("TITLE=901.23457_234.568_index=0_test\nPEPMASS=901.23457\n"),
        "compact output was\n{text}"
    );
}

// ---------------------------------------------------------------------------
// START_SECTION(void store(const std::string&, const PeakMap&, bool))
// ---------------------------------------------------------------------------

#[test]
fn store_to_file_round_trips() {
    let directory = std::env::temp_dir().join("openms_mascot_generic_store");
    std::fs::create_dir_all(&directory).unwrap();
    let path = directory.join("MascotGenericFile_1.mgf");
    let experiment = mgf::read(INFILE.as_bytes()).unwrap();
    let mut file = writer();
    let report = file.store(&path, &experiment, false).unwrap();
    assert_eq!(report.written, 1);
    assert_eq!(report.skipped_ms_level, 0);
    let copy = mgf::load(&path).unwrap();
    assert_eq!(experiment.spectra.len(), copy.spectra.len());
    assert_eq!(
        experiment.spectra[0].peaks.len(),
        copy.spectra[0].peaks.len()
    );
    assert!((experiment.spectra[0].rt - copy.spectra[0].rt).abs() < 1e-9);
    assert!(
        (experiment.spectra[0].precursors[0].mz - copy.spectra[0].precursors[0].mz).abs() < 1e-9
    );
    // The source refuses a destination whose suffix is not .mgf.
    let bad = directory.join("MascotGenericFile_1.txt");
    assert!(matches!(
        file.store(&bad, &experiment, false),
        Err(Error::InvalidValue(_))
    ));
    assert!(!bad.exists());
    std::fs::remove_file(&path).ok();
    std::fs::remove_dir(&directory).ok();
}

// ---------------------------------------------------------------------------
// START_SECTION(COMPOUND_NAME to Metabolite_Name mapping)
// ---------------------------------------------------------------------------

#[test]
fn compound_name_maps_to_the_metabolite_name_user_param() {
    let text = "BEGIN IONS\n\
                TITLE=Test spectrum\n\
                PEPMASS=500.0\n\
                CHARGE=1\n\
                COMPOUND_NAME=Caffeine\n\
                100.0 1000.0\n\
                200.0 2000.0\n\
                END IONS\n";
    let experiment = mgf::read(text.as_bytes()).unwrap();
    assert_eq!(experiment.spectra.len(), 1);
    let spectrum = &experiment.spectra[0];
    assert!(spectrum.metadata.contains_key(MSM_METABOLITE_NAME));
    assert_eq!(
        spectrum.metadata[MSM_METABOLITE_NAME].as_str().unwrap(),
        "Caffeine"
    );
    assert_eq!(spectrum.peaks.len(), 2);
    assert!((spectrum.precursors[0].mz - 500.0).abs() < 1e-9);
    assert_eq!(spectrum.precursors[0].charge, 1);
}

// ---------------------------------------------------------------------------
// START_SECTION(GNPS MGF file - 3-Des-Microcystein_LR)
// ---------------------------------------------------------------------------

#[test]
fn gnps_library_spectrum() {
    let experiment = mgf::read(GNPS.as_bytes()).unwrap();
    assert_eq!(experiment.spectra.len(), 1);
    let spectrum = &experiment.spectra[0];
    assert!(spectrum.metadata.contains_key("GNPS_Spectrum_ID"));
    assert_eq!(
        spectrum.metadata["GNPS_Spectrum_ID"].as_str().unwrap(),
        "CCMSLIB00000001547"
    );
    assert!(spectrum.metadata.contains_key(MSM_METABOLITE_NAME));
    assert_eq!(
        spectrum.metadata[MSM_METABOLITE_NAME].as_str().unwrap(),
        "3-Des-Microcystein_LR"
    );
    assert!((spectrum.precursors[0].mz - 981.54).abs() < 1e-9);
    assert_eq!(spectrum.precursors[0].charge, 1);
    assert_eq!(spectrum.peaks.len(), 43);
    let base = spectrum
        .peaks
        .iter()
        .find(|p| (p.mz - 599.352783).abs() < 0.001)
        .expect("the base peak is present");
    assert!((f64::from(base.intensity) - 764_523.0).abs() < 1e-3);
    // The fixture's ADDUCT= and ION_MODE= lines match no source prefix and are
    // therefore silently dropped; SMILES= and SCANS= are kept.
    assert!(!spectrum.metadata.contains_key("ADDUCT"));
    assert!(!spectrum.metadata.contains_key("ION_MODE"));
    assert!(!spectrum.metadata.contains_key("IONMODE"));
    assert!(spectrum.metadata.contains_key(MSM_SMILES_STRING));
    assert_eq!(spectrum.metadata["Scan_ID"].as_str().unwrap(), "1");
    assert_eq!(spectrum.ms_level, 2);
}

// ---------------------------------------------------------------------------
// START_SECTION(SEQ sequence query field - single and multiple)
// ---------------------------------------------------------------------------

#[test]
fn seq_single_is_a_one_element_string_list_and_round_trips() {
    let text = "BEGIN IONS\n\
                TITLE=seq_single\n\
                PEPMASS=500.0\n\
                CHARGE=2+\n\
                SEQ=PEPTIDER\n\
                100.0 1000.0\n\
                200.0 2000.0\n\
                END IONS\n";
    let experiment = mgf::read(text.as_bytes()).unwrap();
    assert_eq!(experiment.spectra.len(), 1);
    let sequences = experiment.spectra[0].metadata["SEQ"]
        .as_string_list()
        .unwrap();
    assert_eq!(sequences.len(), 1);
    assert_eq!(sequences[0], "PEPTIDER");

    let mut file = writer();
    let out = store_to_string(&mut file, "test", &experiment, false);
    let copy = mgf::read(out.as_bytes()).unwrap();
    assert_eq!(copy.spectra.len(), 1);
    let round_trip = copy.spectra[0].metadata["SEQ"].as_string_list().unwrap();
    assert_eq!(round_trip.len(), 1);
    assert_eq!(round_trip[0], "PEPTIDER");
}

#[test]
fn seq_multiple_accumulates_in_order_and_round_trips() {
    let text = "BEGIN IONS\n\
                TITLE=seq_multi\n\
                PEPMASS=600.0\n\
                CHARGE=2+\n\
                SEQ=PEPTIDEA\n\
                SEQ=PEPTIDEB\n\
                SEQ=PEPTIDEC\n\
                100.0 1000.0\n\
                END IONS\n";
    let experiment = mgf::read(text.as_bytes()).unwrap();
    let sequences = experiment.spectra[0].metadata["SEQ"]
        .as_string_list()
        .unwrap();
    assert_eq!(sequences.len(), 3);
    assert_eq!(sequences[0], "PEPTIDEA");
    assert_eq!(sequences[1], "PEPTIDEB");
    assert_eq!(sequences[2], "PEPTIDEC");

    let mut file = writer();
    let out = store_to_string(&mut file, "test", &experiment, false);
    let copy = mgf::read(out.as_bytes()).unwrap();
    let round_trip = copy.spectra[0].metadata["SEQ"].as_string_list().unwrap();
    assert_eq!(round_trip.len(), 3);
    assert_eq!(round_trip[0], "PEPTIDEA");
    assert_eq!(round_trip[1], "PEPTIDEB");
    assert_eq!(round_trip[2], "PEPTIDEC");
}

#[test]
fn seq_set_programmatically_is_written_in_both_store_modes() {
    let spectrum = MSSpectrum {
        native_id: "index=0".into(),
        ms_level: 2,
        rt: 100.0,
        precursors: vec![Precursor::new(500.0, 2)],
        peaks: vec![Peak1D::new(100.0, 1000.0)],
        metadata: [(
            "SEQ".to_owned(),
            openms::metadata::MetaValue::from(vec!["PEPTIDER".to_owned()]),
        )]
        .into_iter()
        .collect(),
        ..Default::default()
    };
    let experiment = MSExperiment {
        spectra: vec![spectrum],
        ..Default::default()
    };
    let mut file = writer();
    assert!(store_to_string(&mut file, "test", &experiment, false).contains("SEQ=PEPTIDER"));
    let mut file = writer();
    assert!(store_to_string(&mut file, "test", &experiment, true).contains("SEQ=PEPTIDER"));
}

/// A query with the maximum number of `SEQ=` lines is linear, not quadratic.
///
/// The source reads the accumulated list out of the meta value, appends one
/// entry and writes the whole list back on *every* `SEQ=` line, so a 600 kB
/// block with 100,000 of them copies about five billion strings. Found by a
/// second-model review; the list is accumulated in the reader and stored once
/// when the block ends. Measured on the Linux gate node in release mode:
/// 433.9 s before, 0.04 s after.
#[test]
fn many_seq_lines_in_one_query_are_linear() {
    let mut text = String::from("BEGIN IONS\nPEPMASS=500.0\n");
    for index in 0..100_000 {
        text.push_str("SEQ=PEPTIDER");
        text.push_str(&(index % 10).to_string());
        text.push('\n');
    }
    text.push_str("100.0 5.0\nEND IONS\n");
    let experiment = mgf::read(text.as_bytes()).unwrap();
    let sequences = experiment.spectra[0].metadata["SEQ"]
        .as_string_list()
        .unwrap()
        .to_vec();
    assert_eq!(sequences.len(), 100_000);
    assert_eq!(sequences[0], "PEPTIDER0");
    assert_eq!(sequences[99_999], "PEPTIDER9");
    // One line past the ceiling is refused, and the ceiling is the source's
    // only protection against an unbounded list.
    text.insert_str(
        text.len() - "100.0 5.0\nEND IONS\n".len(),
        "SEQ=ONE_TOO_MANY\n",
    );
    assert!(matches!(
        mgf::read(text.as_bytes()),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn seq_does_not_bleed_across_blocks_in_either_carry_over_mode() {
    let text = "BEGIN IONS\n\
                TITLE=first\n\
                PEPMASS=500.0\n\
                SEQ=FIRSTONE\n\
                100.0 1000.0\n\
                END IONS\n\
                BEGIN IONS\n\
                TITLE=second\n\
                PEPMASS=600.0\n\
                200.0 2000.0\n\
                END IONS\n";
    for carry_over in [CarryOver::Reset, CarryOver::Source] {
        let options = ReadOptions {
            carry_over,
            ..Default::default()
        };
        let experiment = mgf::read_with_options(text.as_bytes(), &options).unwrap();
        assert_eq!(experiment.spectra.len(), 2);
        let first = experiment.spectra[0].metadata["SEQ"]
            .as_string_list()
            .unwrap();
        assert_eq!(first.len(), 1);
        assert_eq!(first[0], "FIRSTONE");
        assert!(!experiment.spectra[1].metadata.contains_key("SEQ"));
        assert!(!experiment.spectra[1].metadata.contains_key("TITLE_STALE"));
        assert_eq!(
            experiment.spectra[1].metadata["TITLE"].as_str().unwrap(),
            "second_index=1"
        );
    }
}

// ---------------------------------------------------------------------------
// Native tests: source quirks, no-panic guarantees and resource bounds
// ---------------------------------------------------------------------------

#[test]
fn carry_over_reproduces_the_source_bleed_and_reset_prevents_it() {
    let text = "BEGIN IONS\n\
                PEPMASS=500.0 42\n\
                CHARGE=3+\n\
                RTINSECONDS=11.5\n\
                MSLEVEL=3\n\
                ORGANISM=E. coli\n\
                100.0 1.0\n\
                END IONS\n\
                BEGIN IONS\n\
                PEPMASS=600.0\n\
                200.0 2.0\n\
                END IONS\n";
    let source = mgf::read_with_options(
        text.as_bytes(),
        &ReadOptions {
            carry_over: CarryOver::Source,
            ..Default::default()
        },
    )
    .unwrap();
    let second = &source.spectra[1];
    // Only peaks, native ID, TITLE and SEQ are reset by getNextSpectrum_.
    assert_eq!(second.precursors[0].charge, 3);
    assert_eq!(second.precursors[0].intensity, 42.0);
    assert_eq!(second.rt, 11.5);
    assert_eq!(second.ms_level, 3);
    assert_eq!(second.metadata["ORGANISM"].as_str().unwrap(), "E. coli");
    assert_eq!(second.native_id, "index=1");

    let reset = mgf::read(text.as_bytes()).unwrap();
    let second = &reset.spectra[1];
    assert_eq!(second.precursors[0].charge, 0);
    assert_eq!(second.precursors[0].intensity, 0.0);
    assert_eq!(second.rt, -1.0);
    assert_eq!(second.ms_level, 2);
    assert!(!second.metadata.contains_key("ORGANISM"));
}

#[test]
fn lines_outside_a_block_are_ignored_including_the_parameter_header() {
    let mut file = writer();
    let experiment = mgf::read(INFILE.as_bytes()).unwrap();
    let text = store_to_string(&mut file, "test", &experiment, false);
    assert!(text.contains("CHARGE=1,2,3\n"));
    // The global CHARGE default and every other header line lie outside a
    // block, so re-reading the file yields one spectrum with charge 0.
    let copy = mgf::read(text.as_bytes()).unwrap();
    assert_eq!(copy.spectra.len(), 1);
    assert_eq!(copy.spectra[0].precursors[0].charge, 0);
}

#[test]
fn a_block_with_no_peak_line_merges_into_the_next() {
    let text = "BEGIN IONS\n\
                PEPMASS=1.0\n\
                END IONS\n\
                BEGIN IONS\n\
                PEPMASS=2.0\n\
                100.0 5.0\n\
                END IONS\n";
    let experiment = mgf::read(text.as_bytes()).unwrap();
    assert_eq!(experiment.spectra.len(), 1);
    assert_eq!(experiment.spectra[0].precursors[0].mz, 2.0);
    assert_eq!(experiment.spectra[0].native_id, "index=0");
}

#[test]
fn truncation_before_the_peak_list_is_silent_and_after_it_is_an_error() {
    let before = "BEGIN IONS\nPEPMASS=1.0\n";
    assert_eq!(mgf::read(before.as_bytes()).unwrap().spectra.len(), 0);
    let after = "BEGIN IONS\nPEPMASS=1.0\n100.0 5.0\n";
    assert!(matches!(
        mgf::read(after.as_bytes()),
        Err(Error::Parse { .. })
    ));
}

#[test]
fn a_header_line_after_the_peak_list_is_a_parse_error() {
    let text = "BEGIN IONS\nPEPMASS=1.0\n100.0 5.0\nCHARGE=2+\nEND IONS\n";
    let error = mgf::read(text.as_bytes()).unwrap_err();
    match error {
        Error::Parse { message, .. } => {
            assert!(
                message.contains("does not contain m/z and intensity"),
                "{message}"
            );
        }
        other => panic!("unexpected error {other:?}"),
    }
}

#[test]
fn charge_accepts_a_sign_suffix_but_rejects_lists_and_negatives() {
    let read_charge = |value: &str| {
        let text = format!("BEGIN IONS\nPEPMASS=1.0\nCHARGE={value}\n100.0 5.0\nEND IONS\n");
        mgf::read(text.as_bytes()).map(|e| e.spectra[0].precursors[0].charge)
    };
    assert_eq!(read_charge("2+").unwrap(), 2);
    assert_eq!(read_charge("+2").unwrap(), 2);
    assert_eq!(read_charge("2").unwrap(), 2);
    // The source strips '+' only, so a trailing '-' leaves "2-" unparsable and
    // a charge list keeps its commas.
    assert!(matches!(read_charge("2-"), Err(Error::Parse { .. })));
    assert!(matches!(read_charge("1,2,3"), Err(Error::Parse { .. })));
    // A leading minus is a complete integer and is accepted.
    assert_eq!(read_charge("-2").unwrap(), -2);
}

#[test]
fn pepmass_accepts_one_or_two_fields_and_rejects_three() {
    let text = "BEGIN IONS\nPEPMASS=500.25 1200.5\n100.0 5.0\nEND IONS\n";
    let experiment = mgf::read(text.as_bytes()).unwrap();
    assert_eq!(experiment.spectra[0].precursors[0].mz, 500.25);
    assert_eq!(experiment.spectra[0].precursors[0].intensity, 1200.5);
    let text = "BEGIN IONS\nPEPMASS=500.25 1200.5 3\n100.0 5.0\nEND IONS\n";
    let error = mgf::read(text.as_bytes()).unwrap_err();
    match error {
        Error::Parse { message, .. } => {
            assert!(message.contains("Cannot parse PEPMASS"), "{message}")
        }
        other => panic!("unexpected error {other:?}"),
    }
}

/// `PEPMASS=` is tab-substituted but *not* whitespace-collapsed, so a double
/// space yields an empty middle field and the three-entry error.
///
/// `MascotGenericFile.h` calls `simplify` only on the peak lines; `PEPMASS`
/// gets `substitute('\t', ' ')` and then a plain `split`, and `StringUtils`
/// keeps empty chunks. Found by a second-model review: the port collapsed the
/// run and accepted m/z 500 with intensity 10.
#[test]
fn pepmass_is_not_whitespace_collapsed() {
    let read = |line: &str| {
        let text = format!("BEGIN IONS\n{line}\n100.0 5.0\nEND IONS\n");
        mgf::read(text.as_bytes())
    };
    match read("PEPMASS=500  10").unwrap_err() {
        Error::Parse { message, .. } => {
            assert!(message.contains("but 3 were present"), "{message}")
        }
        other => panic!("unexpected error {other:?}"),
    }
    // A tab is still an accepted separator.
    let experiment = read("PEPMASS=500\t10").unwrap();
    assert_eq!(experiment.spectra[0].precursors[0].mz, 500.0);
    assert_eq!(experiment.spectra[0].precursors[0].intensity, 10.0);
    // And a single space remains the ordinary two-field form.
    let experiment = read("PEPMASS=500 10").unwrap();
    assert_eq!(experiment.spectra[0].precursors[0].intensity, 10.0);
}

/// `StringUtils::substr` clamps its start position, so a header line shorter
/// than its own key stores an empty value instead of failing.
///
/// Found by a second-model review: `.get(offset..)` treated the out-of-range
/// offset as an error, so a bare `NAME` line was refused where the source
/// stores an empty metabolite name and a bare `MSLEVEL` takes its MS2 fallback.
#[test]
fn a_header_line_shorter_than_its_key_yields_an_empty_value() {
    let text = "BEGIN IONS\nNAME\nPEPMASS=500\n100.0 5.0\nEND IONS\n";
    let experiment = mgf::read(text.as_bytes()).unwrap();
    assert_eq!(
        experiment.spectra[0].metadata[MSM_METABOLITE_NAME]
            .as_str()
            .unwrap(),
        ""
    );
    let text = "BEGIN IONS\nMSLEVEL\nPEPMASS=500\n100.0 5.0\nEND IONS\n";
    let experiment = mgf::read(text.as_bytes()).unwrap();
    assert_eq!(experiment.spectra[0].ms_level, 2);
    assert_eq!(
        experiment.spectra[0].metadata["MSLEVEL"].as_str().unwrap(),
        "2"
    );
    // An offset inside a multi-byte character is still refused, because the
    // source's byte slice would be an ill-formed string.
    let text = "BEGIN IONS\nPEPMASSé=1.0\n100.0 5.0\nEND IONS\n";
    assert!(matches!(
        mgf::read(text.as_bytes()),
        Err(Error::Parse { .. })
    ));
    // A key with nothing after it is still the empty-value conversion error.
    let text = "BEGIN IONS\nPEPMASS\n100.0 5.0\nEND IONS\n";
    assert!(matches!(
        mgf::read(text.as_bytes()),
        Err(Error::Parse { .. })
    ));
}

/// `toDouble`/`toInt32` reject a second `+` and an overflowing decimal.
///
/// `StringUtils` strips one leading `+` and hands the rest to
/// `std::from_chars`, which refuses a `+` of its own and reports an
/// overflowing literal as `result_out_of_range`. Rust's parser accepts both, so
/// the port had to refuse them explicitly. Found by a second-model review.
#[test]
fn numeric_conversion_refuses_a_second_plus_and_an_overflowing_literal() {
    let read = |line: &str| {
        let text = format!("BEGIN IONS\n{line}\n100.0 5.0\nEND IONS\n");
        mgf::read(text.as_bytes())
    };
    assert!(matches!(read("PEPMASS=++5"), Err(Error::Parse { .. })));
    assert_eq!(read("PEPMASS=+5").unwrap().spectra[0].precursors[0].mz, 5.0);
    // `CHARGE=` is the exception: it removes *every* '+' before converting, so
    // a second one never reaches the conversion.
    assert_eq!(
        read("PEPMASS=1.0\nCHARGE=++2").unwrap().spectra[0].precursors[0].charge,
        2
    );
    // An overflowing retention time is a conversion error, not an infinity.
    assert!(matches!(
        read("PEPMASS=1.0\nRTINSECONDS=1e999"),
        Err(Error::Parse { .. })
    ));
    // In the TITLE branch that conversion error takes the fallback path, so the
    // title is stored and the block loads.
    let text = "BEGIN IONS\nTITLE=run, 1e999 min\nPEPMASS=1.0\n100.0 5.0\nEND IONS\n";
    let experiment = mgf::read(text.as_bytes()).unwrap();
    assert_eq!(
        experiment.spectra[0].metadata["TITLE"].as_str().unwrap(),
        "run, 1e999 min"
    );
}

#[test]
fn title_with_minutes_sets_the_retention_time_and_stores_no_title() {
    let text = "BEGIN IONS\n\
                TITLE= Cmpd 1, +MSn(595.3), 10.9 min\n\
                PEPMASS=595.3\n\
                100.0 5.0\n\
                END IONS\n";
    let experiment = mgf::read(text.as_bytes()).unwrap();
    let spectrum = &experiment.spectra[0];
    assert!((spectrum.rt - 654.0).abs() < 1e-9);
    assert!(!spectrum.metadata.contains_key("TITLE"));
}

#[test]
fn title_with_minutes_but_unparsable_falls_back_to_the_first_equals_split() {
    // The chunk "TITLE=a min b" contains "min", so its first token "TITLE=a" is
    // converted and fails, and the source stores only the text between the
    // first and second '='.
    let text = "BEGIN IONS\nTITLE=a min b=c\nPEPMASS=1.0\n100.0 5.0\nEND IONS\n";
    let experiment = mgf::read(text.as_bytes()).unwrap();
    let spectrum = &experiment.spectra[0];
    assert_eq!(spectrum.metadata["TITLE"].as_str().unwrap(), "a min b");
    assert_eq!(spectrum.rt, -1.0);
}

/// A retention time already committed by an earlier chunk survives a later
/// chunk's conversion failure.
///
/// `MascotGenericFile.h` calls `spectrum.setRT(...)` *inside* the chunk loop,
/// immediately after each successful conversion, and the single enclosing
/// `catch` does not undo it — it only adds the fallback `TITLE`. Found by a
/// second-model review: the port staged the retention time and discarded it.
#[test]
fn a_committed_title_retention_time_survives_a_later_failure() {
    let text = "BEGIN IONS\n\
                RTINSECONDS=7\n\
                TITLE=run, 2 min, bad min, 3 min\n\
                PEPMASS=1.0\n\
                100.0 5.0\n\
                END IONS\n";
    let experiment = mgf::read(text.as_bytes()).unwrap();
    let spectrum = &experiment.spectra[0];
    // The second chunk converted and set 120 s; the third failed and aborted
    // the loop, so the fourth chunk is never visited.
    assert!((spectrum.rt - 120.0).abs() < 1e-9, "{}", spectrum.rt);
    assert_eq!(
        spectrum.metadata["TITLE"].as_str().unwrap(),
        "run, 2 min, bad min, 3 min"
    );
    // Two TITLE lines in one block: the first commits 60 s, the second commits
    // 120 s before failing.
    let text = "BEGIN IONS\n\
                TITLE=first, 1 min\n\
                TITLE=second, 2 min, bad min\n\
                PEPMASS=1.0\n\
                100.0 5.0\n\
                END IONS\n";
    let experiment = mgf::read(text.as_bytes()).unwrap();
    assert!((experiment.spectra[0].rt - 120.0).abs() < 1e-9);
}

#[test]
fn a_second_title_line_replaces_the_first_without_the_native_id_suffix() {
    let text = "BEGIN IONS\nTITLE=one\nTITLE=two\nPEPMASS=1.0\n100.0 5.0\nEND IONS\n";
    let experiment = mgf::read(text.as_bytes()).unwrap();
    assert_eq!(
        experiment.spectra[0].metadata["TITLE"].as_str().unwrap(),
        "two"
    );
}

#[test]
fn mslevel_falls_back_to_two_when_unparsable() {
    let unparsable = "BEGIN IONS\nPEPMASS=1.0\nMSLEVEL=abc\n100.0 5.0\nEND IONS\n";
    let experiment = mgf::read(unparsable.as_bytes()).unwrap();
    assert_eq!(experiment.spectra[0].ms_level, 2);
    assert_eq!(
        experiment.spectra[0].metadata["MSLEVEL"].as_str().unwrap(),
        "2"
    );
    // `std::stoi` consumes a numeric prefix and ignores the rest.
    let prefix = "BEGIN IONS\nPEPMASS=1.0\nMSLEVEL=3 (HCD)\n100.0 5.0\nEND IONS\n";
    assert_eq!(mgf::read(prefix.as_bytes()).unwrap().spectra[0].ms_level, 3);
    // Out of `int` range: the source falls back to MS2 silently.
    let huge = "BEGIN IONS\nPEPMASS=1.0\nMSLEVEL=99999999999999\n100.0 5.0\nEND IONS\n";
    let experiment = mgf::read(huge.as_bytes()).unwrap();
    assert_eq!(experiment.spectra[0].ms_level, 2);
    assert!(!experiment.spectra[0].metadata.contains_key("MSLEVEL"));
    // A non-positive level is refused rather than stored.
    let zero = "BEGIN IONS\nPEPMASS=1.0\nMSLEVEL=0\n100.0 5.0\nEND IONS\n";
    assert!(matches!(
        mgf::read(zero.as_bytes()),
        Err(Error::Parse { .. })
    ));
}

#[test]
fn non_ascii_input_is_handled_without_panicking() {
    // A title, a metabolite name and a whole line of non-ASCII text: the source
    // hands `isdigit` a negative char here, which is undefined behaviour.
    let text = "BEGIN IONS\n\
                TITLE=日本語のタイトル\n\
                COMPOUND_NAME=Caffeïne — ☕\n\
                日本語\n\
                PEPMASS=500.0\n\
                100.0 5.0\n\
                END IONS\n";
    let experiment = mgf::read(text.as_bytes()).unwrap();
    let spectrum = &experiment.spectra[0];
    assert_eq!(
        spectrum.metadata["TITLE"].as_str().unwrap(),
        "日本語のタイトル_index=0"
    );
    assert_eq!(
        spectrum.metadata[MSM_METABOLITE_NAME].as_str().unwrap(),
        "Caffeïne — ☕"
    );
    assert_eq!(spectrum.peaks.len(), 1);
    // A key whose value offset lands inside a multi-byte character is refused
    // instead of producing an ill-formed string.
    let split = "BEGIN IONS\nPEPMASSé=1.0\n100.0 5.0\nEND IONS\n";
    let error = mgf::read(split.as_bytes()).unwrap_err();
    match error {
        Error::Parse { message, .. } => {
            assert!(message.contains("multi-byte"), "{message}");
        }
        other => panic!("unexpected error {other:?}"),
    }
    // A non-ASCII file name reduces to the empty stem without slicing bytes.
    let mut file = writer();
    let mut experiment = experiment;
    experiment.spectra[0].metadata.remove("TITLE");
    let out = store_to_string(&mut file, "日本語.mgf", &experiment, false);
    assert!(out.contains("_index=0_\n"), "{out}");
}

#[test]
fn peak_lines_reject_missing_whitespace_and_bad_numbers() {
    let no_space = "BEGIN IONS\nPEPMASS=1.0\n100.0\nEND IONS\n";
    match mgf::read(no_space.as_bytes()).unwrap_err() {
        Error::Parse { message, .. } => {
            assert!(
                message.contains("does not contain m/z and intensity"),
                "{message}"
            );
        }
        other => panic!("unexpected error {other:?}"),
    }
    let bad = "BEGIN IONS\nPEPMASS=1.0\n100.0 abc\nEND IONS\n";
    match mgf::read(bad.as_bytes()).unwrap_err() {
        Error::Parse { message, .. } => {
            assert!(
                message.contains("could not be converted to a number"),
                "{message}"
            );
        }
        other => panic!("unexpected error {other:?}"),
    }
    // A NaN or an infinity parses in the source and is stored silently; this
    // port refuses it, because no consumer of an MSSpectrum accepts one.
    for value in ["1.0 nan", "1.0 inf", "1e400 1.0", "1.0 1e400"] {
        let text = format!("BEGIN IONS\nPEPMASS=1.0\n{value}\nEND IONS\n");
        assert!(
            matches!(mgf::read(text.as_bytes()), Err(Error::Parse { .. })),
            "{value} was accepted"
        );
    }
    // A line starting with 'i' or '-' is not a peak line at all: only an ASCII
    // digit opens the peak list, so "inf 1.0" is read as an unknown header key
    // and the block ends up truncated and silently dropped.
    let leading = "BEGIN IONS\nPEPMASS=1.0\ninf 1.0\nEND IONS\n";
    assert_eq!(mgf::read(leading.as_bytes()).unwrap().spectra.len(), 0);
}

#[test]
fn read_options_filter_by_retention_time_m_z_and_intensity() {
    let text = "BEGIN IONS\n\
                PEPMASS=1.0\n\
                RTINSECONDS=10\n\
                100.0 5.0\n\
                200.0 50.0\n\
                300.0 500.0\n\
                END IONS\n\
                BEGIN IONS\n\
                PEPMASS=2.0\n\
                RTINSECONDS=100\n\
                150.0 7.0\n\
                END IONS\n";
    let options = ReadOptions {
        rt_range: Some(0.0..50.0),
        mz_range: Some(150.0..250.0),
        intensity_range: Some(10.0..100.0),
        ..Default::default()
    };
    let experiment = mgf::read_with_options(text.as_bytes(), &options).unwrap();
    assert_eq!(experiment.spectra.len(), 1);
    assert_eq!(experiment.spectra[0].peaks.len(), 1);
    assert_eq!(experiment.spectra[0].peaks[0], Peak1D::new(200.0, 50.0));
    // Inverted and non-finite ranges are refused before anything is read.
    let broken = ReadOptions {
        rt_range: Some(10.0..0.0),
        ..Default::default()
    };
    assert!(matches!(
        mgf::read_with_options(text.as_bytes(), &broken),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn resource_ceilings_are_enforced_before_allocation() {
    let text = "BEGIN IONS\nPEPMASS=1.0\n100.0 5.0\n200.0 6.0\nEND IONS\n";
    let options = ReadOptions {
        limits: openms::format::mascot_generic::Limits {
            max_peaks: 1,
            ..Default::default()
        },
        ..Default::default()
    };
    assert!(matches!(
        mgf::read_with_options(text.as_bytes(), &options),
        Err(Error::InvalidValue(_))
    ));
    let options = ReadOptions {
        limits: openms::format::mascot_generic::Limits {
            max_spectra: 0,
            ..Default::default()
        },
        ..Default::default()
    };
    assert!(matches!(
        mgf::read_with_options(text.as_bytes(), &options),
        Err(Error::InvalidValue(_))
    ));
    let options = ReadOptions {
        limits: openms::format::mascot_generic::Limits {
            max_bytes: 8,
            ..Default::default()
        },
        ..Default::default()
    };
    assert!(matches!(
        mgf::read_with_options(text.as_bytes(), &options),
        Err(Error::InvalidValue(_))
    ));
    // Limits above the hard ceilings are refused.
    let options = ReadOptions {
        limits: openms::format::mascot_generic::Limits {
            max_peaks: usize::MAX,
            ..Default::default()
        },
        ..Default::default()
    };
    assert!(matches!(
        mgf::read_with_options(text.as_bytes(), &options),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn ten_thousand_peaks_are_refused_by_the_writer() {
    let spectrum = MSSpectrum {
        native_id: "index=0".into(),
        ms_level: 2,
        rt: 1.0,
        precursors: vec![Precursor::new(500.0, 2)],
        peaks: (0..10_000)
            .map(|i| Peak1D::new(f64::from(i) + 1.0, 1.0))
            .collect(),
        ..Default::default()
    };
    let experiment = MSExperiment {
        spectra: vec![spectrum],
        ..Default::default()
    };
    let mut file = writer();
    let mut bytes = Vec::new();
    match file.store_to(&mut bytes, "test", &experiment, false) {
        Err(Error::InvalidValue(message)) => {
            assert!(message.contains("the upper limit is 10,000"), "{message}");
        }
        other => panic!("unexpected result {other:?}"),
    }
}

#[test]
fn a_zero_precursor_and_non_ms2_levels_are_skipped_with_a_report() {
    let experiment = MSExperiment {
        spectra: vec![
            MSSpectrum {
                ms_level: 2,
                rt: 5.0,
                precursors: vec![Precursor::new(0.0, 0)],
                ..Default::default()
            },
            MSSpectrum {
                ms_level: 1,
                rt: 6.0,
                ..Default::default()
            },
            MSSpectrum {
                ms_level: 0,
                rt: 7.0,
                instrument_settings: openms::metadata::InstrumentSettings {
                    scan_mode: openms::metadata::ScanMode::Emission,
                    ..Default::default()
                },
                ..Default::default()
            },
            MSSpectrum {
                native_id: "index=3".into(),
                ms_level: 2,
                rt: 8.0,
                precursors: vec![Precursor::new(700.0, 0)],
                peaks: vec![Peak1D::new(1.0, 1.0)],
                ..Default::default()
            },
        ],
        ..Default::default()
    };
    let mut file = writer();
    let mut bytes = Vec::new();
    let report = file
        .store_to(&mut bytes, "run.mgf", &experiment, false)
        .unwrap();
    assert_eq!(report.written, 1);
    assert_eq!(report.skipped_without_precursor, 1);
    assert_eq!(report.skipped_ms_level, 2);
    assert!(
        report
            .warnings
            .iter()
            .any(|w| w.contains("No precursor m/z information for spectrum with rt 5")),
        "{:?}",
        report.warnings
    );
    assert!(
        report
            .warnings
            .iter()
            .any(|w| w.contains("MSLevel is set to 0")),
        "{:?}",
        report.warnings
    );
    assert!(
        report
            .warnings
            .iter()
            .any(|w| w.contains("no native ID accession")),
        "{:?}",
        report.warnings
    );
}

#[test]
fn the_parameter_header_matches_the_source_line_order_and_formatting() {
    let mut file = writer();
    let mut parameters = file.parameters().clone();
    parameters
        .set_value("decoy", ParamValue::String("true".into()), "", &[])
        .unwrap();
    parameters
        .set_value("number_of_hits", ParamValue::Integer(5), "", &[])
        .unwrap();
    parameters
        .set_value(
            "email",
            ParamValue::String("user@example.org".into()),
            "",
            &[],
        )
        .unwrap();
    file.set_parameters(&parameters).unwrap();
    let empty = MSExperiment::default();
    let text = store_to_string(&mut file, "test", &empty, false);
    let expected = "COM=OpenMS_search\n\
                    USERNAME=OpenMS\n\
                    USEREMAIL=user@example.org\n\
                    FORMAT=Mascot generic\n\
                    TOLU=Da\n\
                    ITOLU=Da\n\
                    FORMVER=1.01\n\
                    DB=MSDB\n\
                    DECOY=1\n\
                    SEARCH=MIS\n\
                    REPORT=5\n\
                    CLE=Trypsin\n\
                    MASS=monoisotopic\n\
                    INSTRUMENT=Default\n\
                    PFA=1\n\
                    TOL=3\n\
                    ITOL=0.3\n\
                    TAXONOMY=All entries\n\
                    CHARGE=1,2,3\n";
    assert_eq!(text, expected);
    // FORMAT must stay inside the first five lines: OpenMS recognises its own
    // MGF files by it when the suffix is not .mgf.
    assert!(
        text.lines().take(5).any(|line| line.starts_with("FORMAT=")),
        "{text}"
    );
}

#[test]
fn report_is_auto_and_optional_lines_are_omitted_by_default() {
    let mut file = writer();
    let mut parameters = file.parameters().clone();
    parameters
        .set_value("search_title", ParamValue::String(String::new()), "", &[])
        .unwrap();
    file.set_parameters(&parameters).unwrap();
    let text = store_to_string(&mut file, "test", &MSExperiment::default(), false);
    assert!(!text.contains("COM="));
    assert!(!text.contains("USEREMAIL="));
    assert!(!text.contains("DECOY="));
    assert!(text.contains("REPORT=AUTO\n"));
}

#[test]
fn http_format_writes_mime_parts_and_the_peak_list_enclosure() {
    let mut file = writer();
    let mut parameters = file.parameters().clone();
    parameters
        .set_value(
            "internal:HTTP_format",
            ParamValue::String("true".into()),
            "",
            &[],
        )
        .unwrap();
    file.set_parameters(&parameters).unwrap();
    let enclosure = file.http_peak_list_enclosure("run.mgf").unwrap();
    assert_eq!(
        enclosure.0,
        "--GZWgAaYKjHFeUaLOLEIOMq\nContent-Disposition: form-data; name=\"FILE\"; filename=\"run.mgf\"\n\n"
    );
    assert_eq!(enclosure.1, "\n\n--GZWgAaYKjHFeUaLOLEIOMq--\n");
    let text = store_to_string(&mut file, "run.mgf", &MSExperiment::default(), false);
    assert!(text.starts_with(
        "--GZWgAaYKjHFeUaLOLEIOMq\nContent-Disposition: form-data; name=\"COM\"\n\nOpenMS_search\n"
    ));
    assert!(text.ends_with("\n\n--GZWgAaYKjHFeUaLOLEIOMq--\n"));
    assert!(text.contains(&enclosure.0));
}

#[test]
fn internal_content_selects_header_or_peak_list_only() {
    let experiment = mgf::read(INFILE.as_bytes()).unwrap();
    for (content, header, block) in [
        ("all", true, true),
        ("header_only", true, false),
        ("peaklist_only", false, true),
    ] {
        let mut file = writer();
        let mut parameters = file.parameters().clone();
        parameters
            .set_value(
                "internal:content",
                ParamValue::String(content.into()),
                "",
                &[],
            )
            .unwrap();
        file.set_parameters(&parameters).unwrap();
        let text = store_to_string(&mut file, "test", &experiment, false);
        assert_eq!(text.contains("FORMVER=1.01"), header, "{content}");
        assert_eq!(text.contains("BEGIN IONS"), block, "{content}");
    }
    // write_header_to alone produces the header a caller embedding its own peak
    // list needs.
    let file = writer();
    let mut bytes = Vec::new();
    file.write_header_to(&mut bytes).unwrap();
    assert!(String::from_utf8(bytes).unwrap().contains("FORMVER=1.01"));
}

#[test]
fn skip_spectrum_charges_suppresses_the_charge_line() {
    let experiment = MSExperiment {
        spectra: vec![MSSpectrum {
            native_id: "index=0".into(),
            ms_level: 2,
            rt: 1.0,
            precursors: vec![Precursor::new(500.0, -3)],
            peaks: vec![Peak1D::new(1.0, 1.0)],
            ..Default::default()
        }],
        ..Default::default()
    };
    let mut file = writer();
    let text = store_to_string(&mut file, "test", &experiment, false);
    // The source writes magnitude and sign separately, so a negative charge
    // comes out as "-3-" and its own reader cannot read it back.
    assert!(text.contains("CHARGE=-3-\n"), "{text}");
    assert!(matches!(
        mgf::read(text.as_bytes()),
        Err(Error::Parse { .. })
    ));

    let mut file = writer();
    let mut parameters = file.parameters().clone();
    parameters
        .set_value(
            "skip_spectrum_charges",
            ParamValue::String("true".into()),
            "",
            &[],
        )
        .unwrap();
    file.set_parameters(&parameters).unwrap();
    let text = store_to_string(&mut file, "test", &experiment, false);
    assert!(!text.contains("CHARGE=-3"), "{text}");
    assert!(text.contains("CHARGE=1,2,3\n"), "{text}");
}

#[test]
fn scans_follows_the_native_id_type_accession_table() {
    let spectrum = |native_id: &str| MSSpectrum {
        native_id: native_id.into(),
        ms_level: 2,
        rt: 1.0,
        precursors: vec![Precursor::new(500.0, 0)],
        peaks: vec![Peak1D::new(1.0, 1.0)],
        ..Default::default()
    };
    let file = writer();
    let scans = |native_id: &str, accession: &str| {
        let mut bytes = Vec::new();
        file.write_spectrum(&mut bytes, &spectrum(native_id), "run", accession)
            .unwrap();
        let text = String::from_utf8(bytes).unwrap();
        text.lines()
            .find_map(|line| line.strip_prefix("SCANS="))
            .expect("a SCANS line")
            .to_owned()
    };
    // UNKNOWN: the text after the native ID's last '='.
    assert_eq!(scans("index=250", "UNKNOWN"), "250");
    assert_eq!(scans("no equals here", "UNKNOWN"), "no equals here");
    // Thermo: scan=NUMBER.
    assert_eq!(
        scans(
            "controllerType=0 controllerNumber=1 scan=11515",
            "MS:1000768"
        ),
        "11515"
    );
    // Bruker TDF: the trailing tokens are ignored and the last scan= wins.
    assert_eq!(scans("frame=3 scan=9 precursor=2", "MS:1002818"), "9");
    // An index= native ID is one less than the scan number consumers expect.
    assert_eq!(scans("index=250", "MS:1000774"), "251");
    // WIFF: cycle * 1000 + experiment.
    assert_eq!(
        scans("sample=1 period=1 cycle=96 experiment=1", "MS:1000770"),
        "96001"
    );
    // A bare number.
    assert_eq!(scans("abc 4711", "MS:1001530"), "4711");
    // Failures write the source's -1 sentinel rather than aborting the file.
    assert_eq!(scans("index=250", "MS:1000768"), "-1");
    assert_eq!(scans("index=250", "MS:9999999"), "-1");
}

/// Only the *last* match is converted, and the WIFF pair is validated only for
/// the final match — both as `SpectrumNativeIDParser` does.
///
/// The token iterator collects every match and then takes `matches.back()`
/// before calling `toInt32`, so an earlier convertible match is not a fallback,
/// and the WIFF branch inspects only `matches[size-2]`/`matches[size-1]`. The
/// conversion is `toInt32`, so the width is 32 bits.
///
/// `experiment >= 1000` raises `Exception::InvalidValue`, which
/// `extractScanNumber` does not catch — only `ConversionError` — so it aborts
/// the store instead of writing the sentinel. Found by a second-model review;
/// this test previously asserted `SCANS=-1` for that input.
#[test]
fn the_last_scan_number_match_wins_before_conversion() {
    let spectrum = |native_id: &str| MSSpectrum {
        native_id: native_id.into(),
        ms_level: 2,
        rt: 1.0,
        precursors: vec![Precursor::new(500.0, 0)],
        peaks: vec![Peak1D::new(1.0, 1.0)],
        ..Default::default()
    };
    let file = writer();
    let written = |native_id: &str, accession: &str| -> Result<String, Error> {
        let mut bytes = Vec::new();
        file.write_spectrum(&mut bytes, &spectrum(native_id), "run", accession)?;
        let text = String::from_utf8(bytes).unwrap();
        Ok(text
            .lines()
            .find_map(|line| line.strip_prefix("SCANS="))
            .expect("a SCANS line")
            .to_owned())
    };
    let scans = |native_id: &str, accession: &str| written(native_id, accession).unwrap();
    // The last match overflows i32, so the conversion fails and the sentinel is
    // written; the earlier, convertible match is not used.
    assert_eq!(scans("scan=7 scan=2147483648", "MS:1000768"), "-1");
    assert_eq!(scans("scan=2147483647", "MS:1000768"), "2147483647");
    // index= adds one in int arithmetic, which overflows at INT_MAX; the port
    // writes the sentinel where the source has signed overflow.
    assert_eq!(scans("index=2147483646", "MS:1000774"), "2147483647");
    assert_eq!(scans("index=2147483647", "MS:1000774"), "-1");
    // WIFF: only the final cycle/experiment pair is examined, so an earlier
    // out-of-range experiment does not spoil a valid later pair.
    assert_eq!(
        scans("cycle=1 experiment=1000 cycle=2 experiment=3", "MS:1000770"),
        "2003"
    );
    // Boost's \s covers vertical tab and form feed.
    assert_eq!(scans("cycle=96\u{b}experiment=1", "MS:1000770"), "96001");
    assert_eq!(scans("cycle=96\u{c}experiment=1", "MS:1000770"), "96001");
    // The final pair's experiment of 1000 or more is an uncaught InvalidValue.
    assert!(matches!(
        written("sample=1 cycle=1 experiment=1000", "MS:1000770"),
        Err(Error::InvalidValue(_))
    ));
    // The same input through a whole store aborts it rather than writing a file.
    let experiment = MSExperiment {
        spectra: vec![spectrum("cycle=1 experiment=1000")],
        settings: ExperimentalSettings {
            source_files: vec![SourceFile {
                native_id_type_accession: "MS:1000770".into(),
                ..Default::default()
            }],
            ..Default::default()
        },
        ..Default::default()
    };
    let mut file = writer();
    let mut bytes = Vec::new();
    assert!(matches!(
        file.store_to(&mut bytes, "test", &experiment, false),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn the_consumer_interface_streams_one_spectrum_at_a_time() {
    let text = "BEGIN IONS\nPEPMASS=1.0\n100.0 1.0\nEND IONS\n\
                BEGIN IONS\nPEPMASS=2.0\n200.0 2.0\nEND IONS\n\
                BEGIN IONS\nPEPMASS=3.0\n300.0 3.0\nEND IONS\n";
    let mut collector = ExperimentCollector::default();
    mgf::consume(text.as_bytes(), &mut collector).unwrap();
    assert_eq!(collector.experiment.spectra.len(), 3);
    assert_eq!(collector.experiment.spectra[2].precursors[0].mz, 3.0);

    struct Stopper(usize);
    impl MSDataConsumer for Stopper {
        fn set_expected_size(&mut self, _s: usize, _c: usize) -> openms::Result<()> {
            Ok(())
        }
        fn set_experimental_settings(
            &mut self,
            _settings: &openms::metadata::ExperimentalSettings,
        ) -> openms::Result<()> {
            Ok(())
        }
        fn consume_spectrum(
            &mut self,
            _spectrum: &mut MSSpectrum,
        ) -> openms::Result<ControlFlow<()>> {
            self.0 += 1;
            Ok(if self.0 == 2 {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            })
        }
        fn consume_chromatogram(
            &mut self,
            _chromatogram: &mut openms::MSChromatogram,
        ) -> openms::Result<ControlFlow<()>> {
            Ok(ControlFlow::Continue(()))
        }
    }
    let mut stopper = Stopper(0);
    mgf::consume(text.as_bytes(), &mut stopper).unwrap();
    assert_eq!(stopper.0, 2);
}

#[test]
fn special_modification_groups_can_be_reconfigured() {
    let mut file = writer();
    let mut parameters = file.parameters().clone();
    parameters
        .set_value(
            "special_modifications",
            ParamValue::String("Phospho (STY)".into()),
            "",
            &[],
        )
        .unwrap();
    file.set_parameters(&parameters).unwrap();
    let groups = file.special_modification_groups();
    assert_eq!(groups.len(), 3);
    assert_eq!(groups["Phospho (Y)"], "Phospho (STY)");
    // An empty list clears the map.
    let mut parameters = file.parameters().clone();
    parameters
        .set_value(
            "special_modifications",
            ParamValue::String(String::new()),
            "",
            &[],
        )
        .unwrap();
    file.set_parameters(&parameters).unwrap();
    assert!(file.special_modification_groups().is_empty());
}

#[test]
fn duplicate_modification_groups_collapse_to_one_line() {
    let mut file = writer();
    let mut parameters = file.parameters().clone();
    parameters
        .set_value(
            "variable_modifications",
            strings(&["Deamidated (N)", "Deamidated (Q)"]),
            "",
            &[],
        )
        .unwrap();
    file.set_parameters(&parameters).unwrap();
    let text = store_to_string(&mut file, "test", &MSExperiment::default(), false);
    assert_eq!(text.matches("IT_MODS=Deamidated (NQ)\n").count(), 1);
}

#[test]
fn an_unknown_parameter_is_reported_as_a_warning_and_a_bad_value_as_an_error() {
    let mut file = writer();
    let mut parameters = file.parameters().clone();
    parameters
        .set_value("not_a_mascot_parameter", ParamValue::Integer(1), "", &[])
        .unwrap();
    let warnings = file.set_parameters(&parameters).unwrap();
    assert!(
        warnings
            .iter()
            .any(|w| w.contains("not_a_mascot_parameter")),
        "{warnings:?}"
    );

    let mut file = writer();
    let mut parameters = file.parameters().clone();
    parameters
        .set_value("search_type", ParamValue::String("NOPE".into()), "", &[])
        .unwrap();
    assert!(file.set_parameters(&parameters).is_err());
    // The refused update left the previous value in place.
    assert_eq!(
        file.parameters().value("search_type").unwrap(),
        &ParamValue::String("MIS".into())
    );
}

#[test]
fn the_output_byte_ceiling_is_checked_before_the_file_is_created() {
    let directory = std::env::temp_dir().join("openms_mascot_generic_ceiling");
    std::fs::create_dir_all(&directory).unwrap();
    let path = directory.join("bounded.mgf");
    let experiment = mgf::read(INFILE.as_bytes()).unwrap();
    let mut file = writer();
    let options = openms::format::mascot_generic::WriteOptions {
        max_output_bytes: 16,
        ..Default::default()
    };
    assert!(
        file.store_with_options(&path, &experiment, false, &options)
            .is_err()
    );
    assert!(!path.exists());
    std::fs::remove_dir_all(&directory).ok();
}
