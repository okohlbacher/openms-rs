// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tests for `src/format/gnps_quantification.rs` and
//! `src/format/gnps_meta_value.rs`, the ports of
//! `FORMAT/GNPSQuantificationFile.h` and `FORMAT/GNPSMetaValueFile.h`.
//!
//! The single `START_SECTION` of each upstream class test (core `bc9cc12`) is
//! reproduced with its literals. The TOPP differentials compare this writer's
//! bytes with the expected outputs of the upstream `GNPSExport` TOPP tests,
//! `GNPSExport_1_out_quant.txt` and `GNPSExport_4_out_meta_values.tsv`
//! (test-data `0cb15f2`), for the same consensusXML inputs; they need the
//! `consensusxml` feature. Fixtures and hashes are recorded in
//! `tests/data/gnps_tables_provenance.json`.

use openms::format::gnps_meta_value::GNPSMetaValueFile;
use openms::format::gnps_quantification::GNPSQuantificationFile;
use openms::kernel::{BaseFeature, ColumnHeader, ConsensusFeature, ConsensusMap, FeatureHandle};
use openms::metadata::MetaValue;
use std::path::{Path, PathBuf};

// Only the consensusxml-gated tests read a fixture; without the feature this would be dead code.
#[cfg_attr(not(feature = "consensusxml"), allow(dead_code))]
fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/gnps_tables")
        .join(name)
}

fn header(filename: &str, label: &str, size: usize) -> ColumnHeader {
    ColumnHeader {
        filename: filename.into(),
        label: label.into(),
        size,
        ..ColumnHeader::default()
    }
}

/// The class tests' map: two input maps.
fn two_maps() -> ConsensusMap {
    let mut cm = ConsensusMap::new();
    cm.column_headers
        .insert(0, header("/data/runA.mzML", "lf", 10));
    cm.column_headers
        .insert(1, header("/data/runB.mzML", "lf", 12));
    cm
}

fn handle(map_index: u64, unique_id: u64, rt: f64, mz: f64, intensity: f32) -> FeatureHandle {
    let mut feature = BaseFeature::new(rt, mz, intensity);
    feature.unique_id = unique_id;
    FeatureHandle::new(map_index, &feature)
}

fn quant_text(cm: &ConsensusMap) -> String {
    let mut out = Vec::new();
    GNPSQuantificationFile::write(cm, &mut out).unwrap();
    String::from_utf8(out).unwrap()
}

// GNPSQuantificationFile_test.cpp:
// START_SECTION((static void store(const ConsensusMap& consensus_map, const std::string& output_file)))
#[test]
fn quantification_store_matches_the_source_section() {
    let mut cm = two_maps();
    let mut cf = ConsensusFeature::new();
    cf.base.rt = 100.0;
    cf.base.mz = 500.0;
    cf.base.intensity = 2000.0;
    cf.base.charge = 2;
    cf.base.unique_id = 42;
    cf.insert(handle(0, 111, 100.1, 500.0, 1200.0)).unwrap();
    cf.insert(handle(1, 222, 99.9, 500.0, 800.0)).unwrap();
    cm.features.push(cf);

    let dir = openms::system::file::TempDir::new(false).unwrap();
    let fn_ = dir.path().join("quant.txt");
    GNPSQuantificationFile::store(&cm, &fn_).unwrap();
    let text = std::fs::read_to_string(&fn_).unwrap();
    let lines: Vec<&str> = text.lines().collect();

    assert_eq!(lines.len(), 5);
    assert_eq!(lines[0], "#MAP\tid\tfilename\tlabel\tsize");
    assert_eq!(
        lines[1],
        "#CONSENSUS\trt_cf\tmz_cf\tintensity_cf\tcharge_cf\twidth_cf\tquality_cf\trt_0\tmz_0\tintensity_0\tcharge_0\twidth_0\trt_1\tmz_1\tintensity_1\tcharge_1\twidth_1"
    );
    assert_eq!(lines[2], "MAP\t0\t/data/runA.mzML\tlf\t10");
    assert_eq!(lines[3], "MAP\t1\t/data/runB.mzML\tlf\t12");
    let c: Vec<&str> = lines[4].split('\t').collect();
    assert_eq!(c.len(), 17);
    assert_eq!(c[0], "CONSENSUS");
    // rt_cf is intentionally not asserted as exact text in the source test.
    assert_eq!(c[2], "500.0"); // mz_cf
    assert_eq!(c[3], "2000.0"); // intensity_cf
    assert_eq!(c[4], "2"); // charge_cf
    assert_eq!(c[9], "1200.0"); // intensity of the map-0 sub-feature
    assert_eq!(c[14], "800.0"); // intensity of the map-1 sub-feature
}

// GNPSMetaValueFile_test.cpp:
// START_SECTION((static void store(const ConsensusMap& consensus_map, const std::string& output_file)))
#[test]
fn meta_value_store_matches_the_source_section() {
    let cm = two_maps();
    let dir = openms::system::file::TempDir::new(false).unwrap();
    let fn_ = dir.path().join("meta.tsv");
    GNPSMetaValueFile::store(&cm, &fn_).unwrap();
    let text = std::fs::read_to_string(&fn_).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), 3);
    assert_eq!(lines[0], "\tfilename\tATTRIBUTE_MAPID");
    assert_eq!(lines[1], "0\trunA.mzML\tMAP0");
    assert_eq!(lines[2], "1\trunB.mzML\tMAP1");
}

// Upstream TOPP test GNPSExport_1: -in_cm GNPSExport_cons1.consensusXML
// -out_quantification, compared with GNPSExport_1_out_quant.txt. The input
// has no IIMN annotation, so the tool writes the loaded map unchanged.
//
// The upstream TOPP test compares with FuzzyDiff, which reads numbers as
// numbers. The retained file spells 52 values (on 44 lines) differently from
// their exactly rounded 15-fraction-digit text in the last digits (for example
// `266.169825393919894` for the double whose expansion is
// `266.16982539391989348...`, and `1177.885337635699671` for
// `1177.88533763569967049...`); every such pair reads back as the same double.
// So every line and field must agree, a text field exactly and a number as the
// identical `f64`, which is the FuzzyDiff criterion with no tolerance. The
// count 52 was derived independently of this crate, by formatting every
// fixed-notation value of the retained file with exact 15-digit rounding.
#[cfg(feature = "consensusxml")]
#[test]
fn quantification_matches_the_gnps_export_1_topp_output() {
    let cm = openms::format::consensusxml::load(data("GNPSExport_cons1.consensusXML")).unwrap();
    let expected = std::fs::read_to_string(data("GNPSExport_1_out_quant.txt")).unwrap();
    let actual = quant_text(&cm);
    let actual_lines: Vec<&str> = actual.split('\n').collect();
    let expected_lines: Vec<&str> = expected.split('\n').collect();
    assert_eq!(actual_lines.len(), expected_lines.len());
    let mut respelled = 0;
    for (line, (a, e)) in actual_lines.iter().zip(&expected_lines).enumerate() {
        let a_fields: Vec<&str> = a.split('\t').collect();
        let e_fields: Vec<&str> = e.split('\t').collect();
        assert_eq!(a_fields.len(), e_fields.len(), "line {line}");
        for (x, y) in a_fields.iter().zip(&e_fields) {
            if x == y {
                continue;
            }
            let (xv, yv) = (x.parse::<f64>(), y.parse::<f64>());
            assert!(
                matches!((xv, yv), (Ok(p), Ok(q)) if p.to_bits() == q.to_bits()),
                "line {line}: {x} != {y}"
            );
            respelled += 1;
        }
    }
    assert_eq!(respelled, 52);
}

// Upstream TOPP test GNPSExport_4: -in_cm GNPSExport_cons2.consensusXML
// -out_meta_values, compared with GNPSExport_4_out_meta_values.tsv. The IIMN
// annotation the tool applies first changes features only, not the column
// headers this table is built from.
#[cfg(feature = "consensusxml")]
#[test]
fn meta_values_match_the_gnps_export_4_topp_output() {
    let cm = openms::format::consensusxml::load(data("GNPSExport_cons2.consensusXML")).unwrap();
    let expected = std::fs::read_to_string(data("GNPSExport_4_out_meta_values.tsv")).unwrap();
    let mut out = Vec::new();
    GNPSMetaValueFile::write(&cm, &mut out).unwrap();
    assert_eq!(String::from_utf8(out).unwrap(), expected);
}

// ----------------------------------------------------------------------------
// Native regressions derived from the two .cpp files
// ----------------------------------------------------------------------------

#[test]
fn iimn_columns_follow_the_first_feature() {
    let mut cm = two_maps();
    let mut first = ConsensusFeature::new();
    first
        .base
        .metadata
        .insert("row ID".into(), MetaValue::from(7i64));
    first
        .base
        .metadata
        .insert("best ion".into(), MetaValue::from("[M+H]+"));
    first.insert(handle(1, 5, 1.0, 2.0, 3.0)).unwrap();
    let second = ConsensusFeature::new();
    cm.features = vec![first, second];
    let text = quant_text(&cm);
    let lines: Vec<&str> = text.lines().collect();
    assert!(lines[1].starts_with(
        "#CONSENSUS\trt_cf\tmz_cf\tintensity_cf\tcharge_cf\twidth_cf\tquality_cf\trow ID\tbest ion\tpartners\tannotation network number\trt_0"
    ));
    let c: Vec<&str> = lines[4].split('\t').collect();
    assert_eq!(&c[7..11], &["7", "[M+H]+", "", ""]);
    // Map 0 has no handle: five empty fields, then map 1.
    assert_eq!(&c[11..16], &["", "", "", "", ""]);
    assert_eq!(&c[16..], &["1.0", "2.0", "3.0", "0", "0.0"]);
    // The second feature lacks every IIMN value.
    let c: Vec<&str> = lines[5].split('\t').collect();
    assert_eq!(&c[7..11], &["", "", "", ""]);
    // Without `row ID` on the first feature there are no IIMN columns, even
    // when a later feature has one.
    cm.features.swap(0, 1);
    assert!(!quant_text(&cm).contains("row ID"));
}

#[test]
fn handles_beyond_the_column_count_are_not_written() {
    let mut cm = ConsensusMap::new();
    cm.column_headers.insert(0, header("a.mzML", "", 1));
    cm.column_headers.insert(5, header("b.mzML", "", 1));
    let mut cf = ConsensusFeature::new();
    cf.insert(handle(5, 1, 1.0, 2.0, 3.0)).unwrap();
    cm.features.push(cf);
    let text = quant_text(&cm);
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines[3], "MAP\t5\tb.mzML\t\t1");
    // Two maps give columns 0 and 1 only; the handle of map 5 is dropped.
    assert_eq!(lines[4].split('\t').count(), 7 + 10);
    assert!(lines[4].ends_with("\t\t\t\t\t\t\t\t\t\t"));
}

#[test]
fn text_fields_substitute_tabs_and_refuse_newlines() {
    let mut cm = ConsensusMap::new();
    cm.column_headers
        .insert(0, header("dir\\with\ttab.mzML", "l\tx", 1));
    let text = quant_text(&cm);
    assert!(text.contains("MAP\t0\tdir\\with_tab.mzML\tl_x\t1\n"));
    let mut meta = Vec::new();
    GNPSMetaValueFile::write(&cm, &mut meta).unwrap();
    assert_eq!(
        String::from_utf8(meta).unwrap(),
        "\tfilename\tATTRIBUTE_MAPID\n0\twith_tab.mzML\tMAP0\n"
    );
    cm.column_headers.insert(1, header("bad\nname", "", 1));
    assert!(GNPSQuantificationFile::write(&cm, Vec::new()).is_err());
}

#[test]
fn empty_map_writes_headers_only() {
    let cm = ConsensusMap::new();
    assert_eq!(
        quant_text(&cm),
        "#MAP\tid\tfilename\tlabel\tsize\n#CONSENSUS\trt_cf\tmz_cf\tintensity_cf\tcharge_cf\twidth_cf\tquality_cf\n"
    );
    let mut meta = Vec::new();
    GNPSMetaValueFile::write(&cm, &mut meta).unwrap();
    assert_eq!(
        String::from_utf8(meta).unwrap(),
        "\tfilename\tATTRIBUTE_MAPID\n"
    );
}
