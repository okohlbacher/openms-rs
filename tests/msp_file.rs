// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tests for `src/format/msp.rs`, the port of `FORMAT/MSPFile.h`.
//!
//! Every `START_SECTION` of the upstream class test `MSPFile_test.cpp` (core
//! `bc9cc12`) is reproduced with its literals and its fixture
//! `MSPFile_test.msp`. The source's `AnnotatedMSRun` overloads are exercised
//! through the `(experiment, identifications)` pair that stands in for it.
//! The writer is also checked byte for byte against the executed Release
//! build (tier 1, `MSPFile_release_store_*.msp`). The remaining tests are
//! independently derived from the source text (tier 4). Fixture hashes are recorded in
//! `tests/data/msp_file_provenance.json`.

use openms::Error;
use openms::format::FileType;
use openms::format::msp::{InstrumentFilter, MSPFile, ReadLimits};
use openms::kernel::DataArray;
use openms::metadata::MetaValue;
use openms::param::ParamValue;
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/msp_file")
        .join(name)
}

fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("msp_file_{tag}_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

// START_SECTION(MSPFile(const MSPFile &rhs)) and
// START_SECTION(MSPFile& operator=(const MSPFile &rhs)).
#[test]
fn copies_carry_the_parameters() {
    let f2 = MSPFile::new();
    let mut p = f2.to_param().unwrap();
    p.set_value("instrument", ParamValue::String("it".into()), "", &[])
        .unwrap();
    let f1 = MSPFile::from_param(&p).unwrap();
    assert_ne!(f1.to_param().unwrap(), f2.to_param().unwrap());
    let f3 = f1;
    assert_eq!(f1.to_param().unwrap(), f3.to_param().unwrap());
    assert_eq!(f3.instrument, InstrumentFilter::It);
}

// START_SECTION(void load(const std::string &filename,
// std::vector< PeptideIdentification > &ids, PeakMap &exp))
#[test]
fn load_matches_the_class_test() {
    let msp_file = MSPFile::new();
    let (ids, exp) = msp_file.load(data("MSPFile_test.msp")).unwrap();
    assert_eq!(exp.len(), 7);
    assert_eq!(ids.len(), 7);
    // Loading again returns fresh containers (the source's regression check
    // that `ids` is replaced, not appended to).
    let (ids, exp) = msp_file.load(data("MSPFile_test.msp")).unwrap();
    assert_eq!(exp.len(), 7);
    assert_eq!(ids.len(), 7);

    // DocumentIdentifier.
    assert_eq!(
        exp.settings.document.loaded_file_path,
        data("MSPFile_test.msp").to_str().unwrap()
    );
    assert_eq!(exp.settings.document.loaded_file_type, FileType::Msp);
    for (i, spectrum) in exp.spectra.iter().enumerate() {
        assert_eq!(spectrum.native_id, format!("index={i}"));
    }
    let hit = &ids[5].hits[0];
    assert_eq!(
        hit.sequence.to_string(),
        ".(Acetyl)AAAAAAGAGPEM(Oxidation)VR"
    );
    assert_eq!(hit.peak_annotations[0].annotation, "a3");
    assert_eq!(hit.peak_annotations[1].annotation, "b3");
    assert_eq!(hit.peak_annotations[1].charge, 1);
    assert_eq!(hit.peak_annotations[2].annotation, "?");
    assert_eq!(hit.peak_annotations[2].charge, 0);
    assert_eq!(hit.peak_annotations[3].annotation, "y4");
    assert_eq!(hit.peak_annotations[3].charge, 2);
    assert_eq!(
        ids[6].hits[0].sequence.to_string(),
        ".(Acetyl)AAAAAAVGPGAGGAGSAVPGGAGPC(Carbamidomethyl)ATVSVFPGAR"
    );

    let qtof = MSPFile {
        instrument: InstrumentFilter::Qtof,
        ..msp_file
    };
    let (ids, exp) = qtof.load(data("MSPFile_test.msp")).unwrap();
    assert_eq!(exp.len(), 2);
    assert_eq!(ids.len(), 2);
    assert_eq!(exp.spectra[0].native_id, "index=0");
    assert_eq!(exp.spectra[1].native_id, "index=1");

    let it = MSPFile {
        instrument: InstrumentFilter::It,
        ..msp_file
    };
    let (ids, exp) = it.load(data("MSPFile_test.msp")).unwrap();
    assert_eq!(exp.len(), 5);
    assert_eq!(ids.len(), 5);
    assert_eq!(exp.spectra[0].native_id, "index=2");
    assert_eq!(exp.spectra[1].native_id, "index=3");
    assert_eq!(exp.spectra[2].native_id, "index=4");
}

// START_SECTION(void store(const std::string& filename, const AnnotatedMSRun&
// annot_exp) const)
#[test]
fn store_round_trip_matches_the_class_test() {
    let msp_file = MSPFile::new();
    let (ids, exp) = msp_file.load(data("MSPFile_test.msp")).unwrap();
    let dir = temp_dir("store");
    let filename = dir.join("out.msp");
    msp_file.store(&filename, &exp, &ids).unwrap();
    let (ids, exp) = msp_file.load(&filename).unwrap();
    assert_eq!(ids.len(), 7);
    assert_eq!(exp.len(), 7);
    for id in &ids[..5] {
        assert_eq!(id.hits.len(), 1);
    }
    let modified: Vec<bool> = ids[..5]
        .iter()
        .map(|id| id.hits[0].sequence.is_modified())
        .collect();
    assert_eq!(modified, [false, false, false, true, false]);
    let charges: Vec<i32> = ids[..5].iter().map(|id| id.hits[0].charge).collect();
    assert_eq!(charges, [2, 2, 2, 2, 3]);
    std::fs::remove_dir_all(&dir).unwrap();
}

// MSPFile.cpp:25-39.
#[test]
fn defaults_match_the_source_constructor() {
    let p = MSPFile::defaults().unwrap();
    let text = |key: &str| p.value(key).unwrap().as_str().unwrap().to_owned();
    assert_eq!(text("parse_headers"), "false");
    assert_eq!(text("parse_peakinfo"), "true");
    assert_eq!(text("parse_firstpeakinfo_only"), "true");
    assert_eq!(text("instrument"), "");
    assert_eq!(MSPFile::from_param(&p).unwrap(), MSPFile::new());
    let mut bad = p;
    bad.set_value("instrument", ParamValue::String("orbitrap".into()), "", &[])
        .unwrap();
    assert!(MSPFile::from_param(&bad).is_err());
}

// Peak m/z and intensity go through StringUtils::toFloat (MSPFile.cpp:262-275);
// MW sets the precursor to (MW + z * PROTON_MASS_U) / z (MSPFile.cpp:124-133).
#[test]
fn peaks_are_floats_and_mw_sets_the_precursor() {
    let (_, exp) = MSPFile::new().load(data("MSPFile_test.msp")).unwrap();
    let first = &exp.spectra[0];
    assert_eq!(first.peaks.len(), 14);
    assert_eq!(first.peaks[0].mz, f64::from(608.7_f32));
    assert_eq!(first.peaks[0].intensity, 974.0);
    assert_eq!(first.precursors.len(), 1);
    let expected = (1857.918 + 2.0 * openms::concept::constants::PROTON_MASS_U) / 2.0;
    assert_eq!(first.precursors[0].mz, expected);
    // Headers are not parsed by default.
    assert!(first.metadata.is_empty());
}

// parseHeader_ (MSPFile.cpp:353-371): trimmed key=value tokens with one '='.
#[test]
fn parse_headers_stores_comment_tokens() {
    let reader = MSPFile {
        parse_headers: true,
        ..MSPFile::new()
    };
    let (_, exp) = reader.load(data("MSPFile_test.msp")).unwrap();
    let meta = &exp.spectra[0].metadata;
    assert_eq!(meta["Spec"], MetaValue::from("Consensus"));
    assert_eq!(meta["Inst"], MetaValue::from("qtof"));
    assert_eq!(meta["Protein"], MetaValue::from("\"TNNC1_BOVIN\""));
    assert!(!meta.contains_key("Comment:"));
}

// parse_firstpeakinfo_only = false keeps every comma-separated annotation; the
// class test's commented-out expectation for ids[5] (y2-H2O, charge 1).
#[test]
fn all_peak_annotations_can_be_kept() {
    let reader = MSPFile {
        parse_firstpeakinfo_only: false,
        ..MSPFile::new()
    };
    let (ids, _) = reader.load(data("MSPFile_test.msp")).unwrap();
    let annotations = &ids[5].hits[0].peak_annotations;
    assert_eq!(annotations[2].annotation, "y2-H2O");
    assert_eq!(annotations[2].charge, 1);
    assert_eq!(annotations.len(), 5);
}

#[test]
fn peak_info_can_be_skipped() {
    let reader = MSPFile {
        parse_peakinfo: false,
        ..MSPFile::new()
    };
    let (ids, exp) = reader.load(data("MSPFile_test.msp")).unwrap();
    assert!(ids.iter().all(|id| id.hits[0].peak_annotations.is_empty()));
    assert_eq!(exp.spectra[0].peaks.len(), 14);
}

fn read(text: &str) -> openms::Result<usize> {
    MSPFile::new()
        .read(text.as_bytes(), ReadLimits::default())
        .map(|(ids, _)| ids.len())
}

// Defensive refusals where the source indexes without checks.
#[test]
fn malformed_entries_are_refused() {
    let parse_error = |text: &str| {
        let err = read(text).unwrap_err();
        assert!(matches!(err, Error::Parse { .. }), "{text:?}: {err:?}");
    };
    parse_error("Name:PEPTIDE/2\n");
    parse_error("Name: PEPTIDE\n");
    parse_error("Name: PEPTIDE/x\n");
    parse_error("Num peaks: 1\n1.0\t2.0\t\"a\"\n");
    parse_error("Name: PEPTIDE/2\nNum peaks: 1\n1.0\n");
    // parse_peakinfo requires the annotation field, quoted.
    parse_error("Name: PEPTIDE/2\nNum peaks: 1\n1.0\t2.0\n");
    parse_error("Name: PEPTIDE/2\nNum peaks: 1\n1.0\t2.0\ta\n");
    parse_error("Name: PEPTIDE/2\nNum peaks: 1\n1.0\t2.0\t\"/x\"\n");
    parse_error("Name: PEPTMDE/2\nComment: Mods=2/4,M,Oxidation\n");
    parse_error("Name: PEPTIDE/2\nComment: Mods=1/1,E\n");
    // A charge of 0 would divide by zero.
    assert!(matches!(
        read("Name: PEPTIDE/0\nMW: 100.0\n").unwrap_err(),
        Error::InvalidValue(_)
    ));
}

// The line that ends a peak block is consumed (MSPFile.cpp:244).
#[test]
fn the_line_after_a_peak_block_is_consumed() {
    let text = "Name: PEPTIDE/2\nNum peaks: 1\n1.0\t2.0\t\"a\"\nName: PEPTIDEK/2\n";
    assert_eq!(read(text).unwrap(), 1);
    let text = "Name: PEPTIDE/2\nNum peaks: 1\n1.0\t2.0\t\"a\"\n\nName: PEPTIDEK/2\n";
    assert_eq!(read(text).unwrap(), 2);
}

#[test]
fn limits_are_enforced() {
    let limits = ReadLimits {
        max_peaks: 13,
        ..ReadLimits::default()
    };
    let file = std::fs::File::open(data("MSPFile_test.msp")).unwrap();
    let err = MSPFile::new()
        .read(std::io::BufReader::new(file), limits)
        .unwrap_err();
    assert!(matches!(err, Error::InvalidValue(_)), "{err:?}");
    let limits = ReadLimits {
        max_spectra: 6,
        ..ReadLimits::default()
    };
    let err = MSPFile::new()
        .load_with_limits(data("MSPFile_test.msp"), limits)
        .unwrap_err();
    assert!(matches!(err, Error::InvalidValue(_)), "{err:?}");
}

// store (MSPFile.cpp:374-489): layout of one entry, intensities scaled to a
// maximum of 10000, ion names from the IonNames array.
#[test]
fn write_layout_follows_the_source() {
    let (ids, mut exp) = MSPFile::new().load(data("MSPFile_test.msp")).unwrap();
    let mut out = Vec::new();
    MSPFile::new().write(&exp, &ids, &mut out).unwrap();
    let text = String::from_utf8(out).unwrap();
    let blocks: Vec<&str> = text.split("\n\n").collect();
    // Entry 3: AM(O)FDIF..., Mods=1/1,M,Oxidation.
    let lines: Vec<&str> = blocks[3].lines().collect();
    assert_eq!(lines[0], "Name: AM(O)FDIFVLGAEDGCISTK/2");
    assert!(lines[1].starts_with("MW: "));
    assert_eq!(lines[2], "Comment: Mods=1/1,M,Oxidation Inst=it");
    assert_eq!(lines[3], "Num peaks: 14");
    // 608.7f32 prints as 608.7 with six significant digits; 974 of a maximum
    // of 10000 stays 974.
    assert_eq!(lines[4], "608.7\t974\t\"?\"");
    // Entry 5: N-terminal acetyl is position 0, residue M11 oxidised.
    let lines: Vec<&str> = blocks[5].lines().collect();
    assert_eq!(lines[0], "Name: AAAAAAGAGPEM(O)VR/2");
    assert_eq!(
        lines[2],
        "Comment: Mods=2/0,A,Acetyl/11,M,Oxidation Inst=it"
    );
    // 566 of a maximum of 566 scales to 10000.
    assert_eq!(lines[6], "400.3\t10000\t\"?\"");

    // IonNames entries are written quoted; a short array is refused.
    let names: Vec<String> = (0..14).map(|i| format!("ion{i}")).collect();
    exp.spectra[0].string_data_arrays.push(DataArray {
        name: "IonNames".into(),
        data: names,
        ..DataArray::default()
    });
    let mut out = Vec::new();
    MSPFile::new().write(&exp, &ids, &mut out).unwrap();
    let text = String::from_utf8(out).unwrap();
    assert!(text.contains("608.7\t974\t\"ion0\"\n"));
    exp.spectra[0].string_data_arrays[0].data.truncate(3);
    assert!(MSPFile::new().write(&exp, &ids, Vec::new()).is_err());
    // Mismatched pairs are refused.
    assert!(MSPFile::new().write(&exp, &ids[..3], Vec::new()).is_err());
    // A known foreign extension is refused (FileHandler::hasValidExtension).
    let dir = temp_dir("ext");
    assert!(
        MSPFile::new()
            .store(dir.join("x.mzML"), &exp, &ids)
            .is_err()
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

// Tier 1: the OpenMS4 Release build (ibminode06, install
// openms4-release-bc9cc12-c19e494-174b576, ../oracle/F4-msp, two identical
// runs) loaded MSPFile_test.msp into an AnnotatedMSRun and stored it with the
// default parameters, then again after attaching an IonNames array
// ("ion0".."ion13") to the first spectrum. The port writes the same bytes,
// the `MW:` digits included.
#[test]
fn write_matches_the_release_build_byte_for_byte() {
    let (ids, mut exp) = MSPFile::new().load(data("MSPFile_test.msp")).unwrap();
    let mut out = Vec::new();
    MSPFile::new().write(&exp, &ids, &mut out).unwrap();
    let expected = std::fs::read_to_string(data("MSPFile_release_store_default.msp")).unwrap();
    assert_eq!(String::from_utf8(out).unwrap(), expected);

    let names: Vec<String> = (0..exp.spectra[0].peaks.len())
        .map(|i| format!("ion{i}"))
        .collect();
    exp.spectra[0].string_data_arrays.push(DataArray {
        name: "IonNames".into(),
        data: names,
        ..DataArray::default()
    });
    let dir = temp_dir("release");
    let filename = dir.join("out.msp");
    MSPFile::new().store(&filename, &exp, &ids).unwrap();
    let written = std::fs::read_to_string(&filename).unwrap();
    let expected = std::fs::read_to_string(data("MSPFile_release_store_ion_names.msp")).unwrap();
    assert_eq!(written, expected);
    std::fs::remove_dir_all(&dir).unwrap();
}
