// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tests for `src/format/gnps_mgf.rs`, the port of `FORMAT/GNPSMGFFile.h`.
//!
//! Upstream has no class test for `GNPSMGFFile` at core `bc9cc12`. The TOPP
//! differentials reproduce the three MGF-producing `GNPSExport` TOPP tests of
//! test-data `0cb15f2` (`topp/CMakeLists.txt:989-1001`): the same
//! consensusXML and mzML inputs, the parameters of `GNPSExport_1_mostint.ini`,
//! `GNPSExport_2_merged.ini` and `GNPSExport_3_binsize.ini`, and the retained
//! C++ outputs `GNPSExport_{1,2,3}_out.mgf`. Upstream compares with FuzzyDiff
//! (`FuzzyDiff.ini`: ratio 1.01, absdiff 0.01); this port is held to line
//! equality except for 14 header values per file that the retained files
//! spell inexactly, which must read back as the identical `f64` (see
//! `RESPELLED_PER_FILE`). Fixtures and hashes are recorded in
//! `tests/data/gnps_mgf_provenance.json`.

#![cfg(all(feature = "consensusxml", feature = "mzml"))]

use openms::Error;
use openms::format::gnps_mgf::{
    DEF_COSINE_SIMILARITY, DEF_MERGE_BIN_SIZE, DEF_PEPT_CUTOFF, GNPSMGFFile, OutputType,
};
use openms::identification::PeptideIdentification;
use openms::kernel::{ConsensusFeature, ConsensusMap, FeatureHandle};
use openms::metadata::MetaValue;
use openms::param::ParamValue;
use std::path::{Path, PathBuf};

fn data(dir: &str, name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(dir)
        .join(name)
}

fn mzmls() -> [PathBuf; 2] {
    [
        data("gnps_mgf", "GNPSExport_mzml1.mzML"),
        data("gnps_mgf", "GNPSExport_mzml2.mzML"),
    ]
}

/// The retained outputs were written by an older OpenMS number formatter that
/// spells 14 `PEPMASS=`/`RTINSECONDS=` values per file differently from their
/// exactly rounded 15-fraction-digit text in the last digit (for example
/// `RTINSECONDS=242.627651999999984` for the double whose expansion is
/// `242.62765199999998344...`, which `StringUtils::toStr` at the pin spells
/// `...983`); each such pair reads back as the same double. Upstream FuzzyDiff
/// reads numbers as numbers, so every line must agree exactly, or be the same
/// `KEY=` with a value that parses to the identical `f64`, which is the
/// FuzzyDiff criterion with no tolerance. The count 14 (8 `PEPMASS`, 6
/// `RTINSECONDS`, in each of the three files) was derived independently of
/// this crate, by formatting every header value of the retained files with
/// exact 15-digit rounding in Python's `decimal`.
const RESPELLED_PER_FILE: usize = 14;

fn assert_same_text(actual: &str, expected: &str) {
    let a: Vec<&str> = actual.split('\n').collect();
    let e: Vec<&str> = expected.split('\n').collect();
    assert_eq!(a.len(), e.len(), "line count");
    let mut respelled = 0;
    for (line, (x, y)) in a.iter().zip(&e).enumerate() {
        if x == y {
            continue;
        }
        let same_number = match (x.split_once('='), y.split_once('=')) {
            (Some((kx, vx)), Some((ky, vy))) if kx == ky => {
                matches!((vx.parse::<f64>(), vy.parse::<f64>()), (Ok(p), Ok(q)) if p.to_bits() == q.to_bits())
            }
            _ => false,
        };
        assert!(same_number, "line {}: {x:?} != {y:?}", line + 1);
        respelled += 1;
    }
    assert_eq!(respelled, RESPELLED_PER_FILE);
}

fn run_case(writer: &GNPSMGFFile, expected: &str) {
    let cm =
        openms::format::consensusxml::load(data("gnps_tables", "GNPSExport_cons1.consensusXML"))
            .unwrap();
    let mut out = Vec::new();
    writer.write(&cm, &mzmls(), &mut out).unwrap();
    let expected = std::fs::read_to_string(data("gnps_mgf", expected)).unwrap();
    assert_same_text(&String::from_utf8(out).unwrap(), &expected);
}

// TOPP_GNPSExport_1 (topp/CMakeLists.txt:989-990), GNPSExport_1_mostint.ini:
// output_type most_intense, peptide_cutoff 10, ms2_bin_size 2,
// merged_spectra:cos_similarity 0.85.
#[test]
fn topp_gnps_export_1_most_intense() {
    let writer = GNPSMGFFile {
        output_type: OutputType::MostIntense,
        peptide_cutoff: 10,
        ms2_bin_size: 2.0,
        cos_similarity: 0.85,
    };
    run_case(&writer, "GNPSExport_1_out.mgf");
}

// TOPP_GNPSExport_2 (topp/CMakeLists.txt:995-996), GNPSExport_2_merged.ini:
// output_type merged_spectra, peptide_cutoff 10, ms2_bin_size 2,
// merged_spectra:cos_similarity 0.85.
#[test]
fn topp_gnps_export_2_merged_spectra() {
    let writer = GNPSMGFFile {
        output_type: OutputType::MergedSpectra,
        peptide_cutoff: 10,
        ms2_bin_size: 2.0,
        cos_similarity: 0.85,
    };
    run_case(&writer, "GNPSExport_2_out.mgf");
}

// TOPP_GNPSExport_3 (topp/CMakeLists.txt:999-1000), GNPSExport_3_binsize.ini:
// output_type most_intense, peptide_cutoff 10, ms2_bin_size 0.02,
// merged_spectra:cos_similarity 0.85.
#[test]
fn topp_gnps_export_3_bin_size() {
    let writer = GNPSMGFFile {
        output_type: OutputType::MostIntense,
        peptide_cutoff: 10,
        ms2_bin_size: 0.02,
        cos_similarity: 0.85,
    };
    run_case(&writer, "GNPSExport_3_out.mgf");
}

// The same case through `store`, the source's file entry point.
#[test]
fn store_writes_the_file() {
    let dir = std::env::temp_dir().join(format!("gnps_mgf_store_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let out = dir.join("out.mgf");
    let writer = GNPSMGFFile {
        peptide_cutoff: 10,
        ms2_bin_size: 2.0,
        cos_similarity: 0.85,
        ..GNPSMGFFile::new()
    };
    writer
        .store(
            data("gnps_tables", "GNPSExport_cons1.consensusXML"),
            &mzmls(),
            &out,
        )
        .unwrap();
    let expected = std::fs::read_to_string(data("gnps_mgf", "GNPSExport_1_out.mgf")).unwrap();
    assert_same_text(&std::fs::read_to_string(&out).unwrap(), &expected);
    std::fs::remove_dir_all(&dir).unwrap();
}

// Defaults from GNPSMGFFile.h:61-72 and the constructor, GNPSMGFFile.cpp:31-43.
#[test]
fn defaults_match_the_source_constructor() {
    assert_eq!(DEF_COSINE_SIMILARITY, 0.9);
    // BinnedSpectrum::DEFAULT_BIN_WIDTH_HIRES is the float 0.02f.
    assert_eq!(DEF_MERGE_BIN_SIZE, f64::from(0.02_f32));
    assert_eq!(DEF_PEPT_CUTOFF, 5);
    let p = GNPSMGFFile::defaults().unwrap();
    assert_eq!(
        p.value("output_type").unwrap(),
        &ParamValue::String("most_intense".into())
    );
    assert_eq!(p.value("peptide_cutoff").unwrap(), &ParamValue::Integer(5));
    assert_eq!(
        p.value("ms2_bin_size").unwrap(),
        &ParamValue::Float(f64::from(0.02_f32))
    );
    assert_eq!(
        p.value("merged_spectra:cos_similarity").unwrap(),
        &ParamValue::Float(0.9)
    );
    assert_eq!(GNPSMGFFile::from_param(&p).unwrap(), GNPSMGFFile::new());
}

#[test]
fn parameters_round_trip_and_restrictions_hold() {
    let writer = GNPSMGFFile {
        output_type: OutputType::MergedSpectra,
        peptide_cutoff: -1,
        ms2_bin_size: 2.0,
        cos_similarity: 0.85,
    };
    let p = writer.to_param().unwrap();
    assert_eq!(GNPSMGFFile::from_param(&p).unwrap(), writer);

    // setValidStrings("output_type", {"merged_spectra","most_intense"}).
    let mut bad = p.clone();
    bad.set_value(
        "output_type",
        ParamValue::String("full_spectra".into()),
        "",
        &[],
    )
    .unwrap();
    assert!(GNPSMGFFile::from_param(&bad).is_err());
    // setMinInt("peptide_cutoff", -1).
    let mut bad = p.clone();
    bad.set_value("peptide_cutoff", ParamValue::Integer(-2), "", &[])
        .unwrap();
    assert!(GNPSMGFFile::from_param(&bad).is_err());
    // setMinFloat("ms2_bin_size", 0).
    let mut bad = p;
    bad.set_value("ms2_bin_size", ParamValue::Float(-0.5), "", &[])
        .unwrap();
    assert!(GNPSMGFFile::from_param(&bad).is_err());
}

fn annotated_feature(map_index: i64, spectrum_index: i64) -> ConsensusFeature {
    let mut feature = ConsensusFeature::default();
    feature.base.mz = 300.0;
    feature
        .insert(FeatureHandle {
            map_index: 0,
            intensity: 10.0,
            ..FeatureHandle::default()
        })
        .unwrap();
    let mut id = PeptideIdentification::default();
    id.metadata
        .insert("map_index".into(), MetaValue::from(map_index));
    id.metadata
        .insert("spectrum_index".into(), MetaValue::from(spectrum_index));
    feature.base.peptide_identifications.push(id);
    feature
}

// The source reads `pepts[0]` of an empty vector; this port refuses.
#[test]
fn a_feature_without_annotation_is_refused() {
    let mut cm = ConsensusMap::new();
    cm.features.push(ConsensusFeature::default());
    let err = GNPSMGFFile::new()
        .write(&cm, &mzmls(), Vec::new())
        .unwrap_err();
    assert!(matches!(err, Error::MissingInformation(_)), "{err:?}");
}

// The source indexes `mzml_file_paths[map_index]` and the experiment without
// checks; this port refuses out-of-range references.
#[test]
fn out_of_range_references_are_refused() {
    let mut cm = ConsensusMap::new();
    cm.features.push(annotated_feature(0, 1_000_000));
    let err = GNPSMGFFile::new()
        .write(&cm, &mzmls(), Vec::new())
        .unwrap_err();
    assert!(matches!(err, Error::InvalidValue(_)), "{err:?}");

    let mut cm = ConsensusMap::new();
    cm.features.push(annotated_feature(0, 0));
    let none: [PathBuf; 0] = [];
    let err = GNPSMGFFile::new()
        .write(&cm, &none, Vec::new())
        .unwrap_err();
    assert!(matches!(err, Error::InvalidValue(_)), "{err:?}");
}

// `(int)getMetaValue("map_index")` throws Exception::ConversionError for a
// non-integer DataValue (DataValue.cpp:519-527).
#[test]
fn a_non_integer_annotation_is_refused() {
    let mut cm = ConsensusMap::new();
    let mut feature = annotated_feature(0, 0);
    feature.base.peptide_identifications[0]
        .metadata
        .insert("map_index".into(), MetaValue::from("0"));
    cm.features.push(feature);
    let err = GNPSMGFFile::new()
        .write(&cm, &mzmls(), Vec::new())
        .unwrap_err();
    assert!(matches!(err, Error::InvalidValue(_)), "{err:?}");
}
