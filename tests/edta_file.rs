// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tests for `src/format/edta.rs`, the port of `FORMAT/EDTAFile.h`.
//!
//! Every `START_SECTION` of the upstream `EDTAFile_test.cpp` (core `bc9cc12`)
//! is reproduced with its literals; the fixtures are the upstream class-test
//! data files `src/tests/class_tests/openms/data/EDTAFile_*` at the same
//! revision, recorded with their sha256 in
//! `tests/data/edta_file_provenance.json`. The remaining tests are native
//! regressions for behaviour read from `EDTAFile.cpp`.

use openms::Error;
use openms::concept::UniqueIdGenerator;
use openms::format::edta::EDTAFile;
use openms::kernel::{
    BaseFeature, ConsensusFeature, ConsensusMap, Feature, FeatureHandle, FeatureMap,
};
use openms::system::file::TempDir;
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/edta_file")
        .join(name)
}

fn ids() -> UniqueIdGenerator {
    UniqueIdGenerator::from_seed(19_991_231_235_959)
}

fn read(text: &str) -> openms::Result<ConsensusMap> {
    EDTAFile::new().read(text.as_bytes(), "in.edta", &mut ids())
}

/// The source's `TEST_REAL_SIMILAR` with its default tolerances: an absolute
/// difference of at most 1e-5 or a ratio within 1 +- 1e-5.
fn similar(a: f64, b: f64) -> bool {
    if (a - b).abs() <= 1e-5 {
        return true;
    }
    if a == 0.0 || b == 0.0 {
        return false;
    }
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    (1.0 / (1.0 + 1e-5)..=1.0 + 1e-5).contains(&ratio)
}

// START_SECTION(EDTAFile()) / START_SECTION(virtual ~EDTAFile())
#[test]
fn construction_and_destruction() {
    let file = EDTAFile::new();
    assert_eq!(file, EDTAFile);
    let _ = file;
}

// START_SECTION(void load(const std::string &filename, ConsensusMap &consensus_map))
#[test]
fn load_matches_the_source_section() {
    let file = EDTAFile::new();
    let mut generator = ids();
    let fm = file
        .load(data("EDTAFile_test_1.edta"), &mut generator)
        .unwrap();
    assert_eq!(fm.len(), 2);
    assert_eq!(fm.features[0].base.rt, 321.0);
    assert_eq!(fm.features[0].base.mz, 405.233);
    assert_eq!(fm.features[0].base.intensity, 24_543_534.0);
    assert_eq!(fm.features[0].base.charge, 2);
    assert_eq!(fm.features[0].base.metadata["mymeta"].to_string(), "lala");
    assert_eq!(fm.features[1].base.rt, 322.0);
    assert_eq!(fm.features[1].base.mz, 406.207);
    assert_eq!(fm.features[1].base.intensity, 4_343_344.0);
    assert_eq!(fm.features[1].base.charge, 3);
    assert_eq!(fm.features[1].base.metadata["mymeta"].to_string(), "blubb");

    let fm = file
        .load(data("EDTAFile_test_3.edta"), &mut generator)
        .unwrap();
    assert_eq!(fm.len(), 3);

    assert!(matches!(
        file.load(data("EDTAFile_test_2.edta"), &mut generator),
        Err(Error::Parse { .. })
    ));
    match file.load(data("EDTAFile_test_3_doesnotexist.edta"), &mut generator) {
        Err(Error::Io(error)) => assert_eq!(error.kind(), std::io::ErrorKind::NotFound),
        other => panic!("expected FileNotFound, got {other:?}"),
    }
}

// START_SECTION((void store(const std::string& filename, const ConsensusMap& map) const))
#[test]
fn consensus_store_round_trips_as_the_source_section() {
    let file = EDTAFile::new();
    let mut generator = ids();
    let cm = file
        .load(data("EDTAFile_test_4.edta"), &mut generator)
        .unwrap();
    let dir = TempDir::new(false).unwrap();
    let outfile = dir.path().join("out.edta");
    file.store_consensus(&outfile, &cm).unwrap();
    let cm2 = file.load(&outfile, &mut generator).unwrap();

    assert_eq!(cm.len(), cm2.len());
    for (a, b) in cm.features.iter().zip(&cm2.features) {
        assert!(similar(a.base.rt, b.base.rt));
        assert!(similar(a.base.mz, b.base.mz));
        assert!(similar(
            f64::from(a.base.intensity),
            f64::from(b.base.intensity)
        ));
        assert_eq!(a.base.charge, b.base.charge);
        assert_eq!(a.len(), b.len());
        // cannot test for meta values, since they are not written to EDTA
    }
}

// START_SECTION((void store(const std::string& filename, const FeatureMap& map) const))
#[cfg(feature = "featurexml")]
#[test]
fn feature_store_round_trips_as_the_source_section() {
    let fm = openms::format::featurexml::load(data("EDTAFile_test_out_1.featureXML")).unwrap();
    let file = EDTAFile::new();
    let dir = TempDir::new(false).unwrap();
    let outfile = dir.path().join("out.edta");
    file.store_features(&outfile, &fm).unwrap();
    let cm = file.load(&outfile, &mut ids()).unwrap();

    assert_eq!(fm.len(), cm.len());
    for (a, b) in fm.features.iter().zip(&cm.features) {
        assert!(similar(a.base.rt, b.base.rt));
        assert!(similar(a.base.mz, b.base.mz));
        assert!(similar(
            f64::from(a.base.intensity),
            f64::from(b.base.intensity)
        ));
    }
}

// ----------------------------------------------------------------------------
// Native regressions derived from EDTAFile.cpp
// ----------------------------------------------------------------------------

#[test]
fn consensus_variant_reads_sub_features_and_skips_na_quadruplets() {
    let cm = EDTAFile::new()
        .load(data("EDTAFile_test_4.edta"), &mut ids())
        .unwrap();
    assert_eq!(cm.len(), 3);
    // Line 2: the second sub-feature is `NA NA NA NA` and is skipped.
    assert_eq!(cm.features[0].len(), 1);
    let handle = cm.features[0].handles()[0];
    assert_eq!(handle.map_index, 0);
    assert_eq!(handle.rt, 331.0);
    assert_eq!(handle.mz, 405.111);
    assert_eq!(handle.intensity, 23_423.0);
    assert_eq!(handle.charge, 2);
    assert_eq!(cm.features[1].len(), 2);
    assert_eq!(cm.features[1].handles()[1].map_index, 1);
    assert_eq!(cm.features[1].handles()[1].mz, 207.2323);
    // Twelve header columns: three quadruplets, two of them sub-features, so
    // two column headers.
    assert_eq!(cm.column_headers.len(), 2);
    for (index, header) in &cm.column_headers {
        assert_eq!(header.label, format!("EDTA_Map {index}"));
        assert_eq!(header.size, 3);
        assert!(header.filename.ends_with("EDTAFile_test_4.edta"));
    }
}

#[test]
fn unique_ids_are_drawn_in_source_order() {
    let text = "RT MZ INT CHARGE RT1 MZ1 INT1 CHARGE1 RT2 MZ2 INT2 CHARGE2\n\
                1 2 3 1 1 2 3 1 NA NA NA NA\n\
                4 5 6 1 4 5 6 1 4 6 6 1\n";
    let map = read(text).unwrap();
    let mut expected = ids();
    let stream: Vec<u64> = (0..7).map(|_| expected.get_unique_id()).collect();
    // map, consensus 1, sub 1, sub 2 (skipped, but its ID is drawn),
    // consensus 2, sub 1, sub 2.
    assert_eq!(map.unique_id, stream[0]);
    assert_eq!(map.features[0].base.unique_id, stream[1]);
    assert_eq!(map.features[0].handles()[0].unique_id, stream[2]);
    assert_eq!(map.features[1].base.unique_id, stream[4]);
    assert_eq!(map.features[1].handles()[0].unique_id, stream[5]);
    assert_eq!(map.features[1].handles()[1].unique_id, stream[6]);
}

#[test]
fn separator_is_chosen_from_the_first_line() {
    let comma = read("1,2,3\n4,5,6\n").unwrap();
    assert_eq!(comma.len(), 2);
    assert_eq!(comma.features[1].base.mz, 5.0);
    // Without a charge column the charge is never set.
    assert_eq!(comma.features[1].base.charge, 0);
    // A tab in the first line wins over spaces.
    let tab = read("RT\tm z\tint\tcharge\n1\t2\t3\t4\n").unwrap();
    assert_eq!(tab.features[0].base.charge, 4);
    // A single column without separators splits on space: too few columns.
    assert!(matches!(read("1\n"), Err(Error::Parse { line: 1, .. })));
    assert!(matches!(read(""), Err(Error::Parse { line: 1, .. })));
}

#[test]
fn header_is_mandatory_beyond_three_columns() {
    let error = read("1 2 3 4\n").unwrap_err();
    match error {
        Error::Parse { line, message } => {
            assert_eq!(line, 1);
            assert!(message.contains("No HEADER provided"), "{message}");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn na_and_missing_columns_become_minus_one() {
    let map = read("RT mz int charge\nNA 2 3\n").unwrap();
    assert_eq!(map.features[0].base.rt, -1.0);
    assert_eq!(map.features[0].base.charge, -1);
}

#[test]
fn a_fourth_cell_of_a_three_column_file_is_still_converted() {
    // The source converts column four to an integer even when the charge is
    // never used, so a non-number there is a parse error.
    assert!(matches!(
        read("1 2 3\n4 5 6 x\n"),
        Err(Error::Parse { line: 2, .. })
    ));
}

#[test]
fn meta_data_needs_a_column_name() {
    let error = read("RT mz int charge\n1 2 3 4 extra\n").unwrap_err();
    match error {
        Error::Parse { line, message } => {
            assert_eq!(line, 2);
            assert!(
                message.contains("Missing meta data header for column 5!"),
                "{message}"
            );
        }
        other => panic!("{other:?}"),
    }
    // Empty meta cells are ignored rather than refused.
    let map = read("RT\tmz\tint\tcharge\tnote\n1\t2\t3\t4\t \n").unwrap();
    assert!(map.features[0].base.metadata.is_empty());
}

#[test]
fn empty_lines_are_skipped() {
    let map = read("1 2 3\n\n   \n4 5 6\n").unwrap();
    assert_eq!(map.len(), 2);
}

#[test]
fn non_finite_values_are_refused() {
    assert!(matches!(
        read("1 2 inf\n"),
        Err(Error::Parse { line: 1, .. })
    ));
    assert!(matches!(
        read("1 2 1e300\n"),
        Err(Error::Parse { line: 1, .. })
    ));
    assert!(matches!(
        read("RT mz int\n1 nan 3\n"),
        Err(Error::Parse { line: 2, .. })
    ));
}

fn consensus_with_handles(map_indices: &[u64]) -> ConsensusMap {
    let mut consensus = ConsensusFeature::new();
    consensus.base = BaseFeature::new(10.0, 500.25, 1000.5);
    consensus.base.charge = 2;
    for (n, &index) in map_indices.iter().enumerate() {
        let mut feature = BaseFeature::new(10.0 + n as f64, 500.25, 100.0);
        feature.unique_id = n as u64 + 1;
        feature.charge = 2;
        consensus
            .insert(FeatureHandle::new(index, &feature))
            .unwrap();
    }
    ConsensusMap::from_features(vec![consensus, ConsensusFeature::new()])
}

#[test]
fn consensus_writer_layout() {
    let map = consensus_with_handles(&[0, 2]);
    let mut out = Vec::new();
    EDTAFile::new().write_consensus(&mut out, &map).unwrap();
    let text = String::from_utf8(out).unwrap();
    assert_eq!(
        text,
        "RT\tm/z\tintensity\tcharge\tRT1\tm/z1\tintensity1\tcharge1\tRT2\tm/z2\tintensity2\tcharge2\n\
         10.0\t500.25\t1000.5\t2\t10.0\t500.25\t100.0\t2\t11.0\t500.25\t100.0\t2\n\
         0.0\t0.0\t0.0\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"
    );
    // Reading it back renumbers the second handle's map from 2 to 1, as the
    // source does: handles are written in consecutive columns.
    let back = read(&text).unwrap();
    assert_eq!(back.features[0].handles()[1].map_index, 1);
}

#[test]
fn consensus_store_checks_the_extension() {
    let dir = TempDir::new(false).unwrap();
    let map = consensus_with_handles(&[0]);
    let file = EDTAFile::new();
    assert!(matches!(
        file.store_consensus(dir.path().join("out.featureXML"), &map),
        Err(Error::InvalidValue(_))
    ));
    assert!(!dir.path().join("out.featureXML").exists());
    // An unknown extension is permitted, as the source's check permits it.
    file.store_consensus(dir.path().join("out.unknownext"), &map)
        .unwrap();
}

#[test]
fn feature_writer_layout() {
    let mut feature = Feature::new(1459.875, 343.25, 2.5e9);
    feature.base.charge = 1;
    let map = FeatureMap::from_features(vec![feature]);
    let mut out = Vec::new();
    EDTAFile::new().write_features(&mut out, &map).unwrap();
    assert_eq!(
        String::from_utf8(out).unwrap(),
        "RT\tm/z\tintensity\tcharge\n1459.875\t343.25\t2.5e09\t1\n"
    );
}
