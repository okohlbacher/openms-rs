// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

#![cfg(feature = "arrow")]

//! The features table of `FeatureMapArrowIO` against the executed C++
//! Release build.
//!
//! Expected values are the records of the executed Linux x86_64 Release
//! install `openms4-release-bc9cc12-c19e494-174b576` (Arrow C++ 23.0.1), run
//! twice on ibminode06 with identical records and files (evidence tier 1):
//! `tests/data/arrow_parquet/cpp_release_records_4.tsv`. The driver imported
//! the features table of the retained bundle
//! `FileConverter_9_output.featureparquet` (test-data `0cb15f2`), exported a
//! constructed FeatureMap to `fm_features.parquet`, read that back, and
//! imported a slice of it that orphans two subordinates.

use arrow_array::RecordBatch;
use openms::format::arrow_io_helpers::{self as aio, ParquetWriteConfig};
use openms::format::feature_map_arrow_io::FeatureMapArrowIO as Fma;
use openms::format::parquet_file::ParquetFile;
use openms::kernel::features::{Feature, FeatureMap};
use openms::kernel::geometry::{ConvexHull2D, Point2D};
use openms::metadata::{MetaInfo, MetaValue, MetaValueData};
use std::path::{Path, PathBuf};

const RECORDS: &str = include_str!("data/arrow_parquet/cpp_release_records_4.tsv");

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/arrow_parquet")
        .join(name)
}

fn records() -> Vec<(String, String, String)> {
    RECORDS
        .lines()
        .map(|line| {
            let mut parts = line.splitn(3, '\t');
            (
                parts.next().unwrap().to_owned(),
                parts.next().unwrap().to_owned(),
                parts.next().unwrap_or("").to_owned(),
            )
        })
        .collect()
}

fn rec(kind: &str, key: &str) -> String {
    records()
        .into_iter()
        .find(|(k, kk, _)| k == kind && kk == key)
        .map(|(_, _, v)| v)
        .unwrap_or_else(|| panic!("no record {kind}\t{key}"))
}

/// C `printf("%a")` as glibc prints a double.
fn hexd(v: f64) -> String {
    let sign = if v.is_sign_negative() { "-" } else { "" };
    let bits = v.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i64;
    let fraction = bits & ((1u64 << 52) - 1);
    if exponent == 0 && fraction == 0 {
        return format!("{sign}0x0p+0");
    }
    let (lead, e) = if exponent == 0 {
        (0, -1022)
    } else {
        (1, exponent - 1023)
    };
    let digits = format!("{fraction:013x}");
    let digits = digits.trim_end_matches('0');
    let mantissa = if digits.is_empty() {
        format!("{lead}")
    } else {
        format!("{lead}.{digits}")
    };
    let plus = if e >= 0 { "+" } else { "" };
    format!("{sign}0x{mantissa}p{plus}{e}")
}

fn describe(value: &MetaValue) -> String {
    let list = |items: Vec<String>| format!("[{}]", items.join(", "));
    match value.data() {
        MetaValueData::String(text) => format!("STRING_VALUE|{text}"),
        MetaValueData::Integer(n) => format!("INT_VALUE|{n}"),
        MetaValueData::Float(x) => format!("DOUBLE_VALUE|{}", hexd(*x)),
        MetaValueData::StringList(v) => format!("STRING_LIST|{}", list(v.clone())),
        MetaValueData::IntegerList(v) => {
            format!("INT_LIST|{}", list(v.iter().map(i64::to_string).collect()))
        }
        MetaValueData::FloatList(v) => {
            let mut text = "DOUBLE_LIST".to_owned();
            for x in v {
                text.push('|');
                text.push_str(&hexd(*x));
            }
            text
        }
        MetaValueData::Empty => "EMPTY_VALUE|".to_owned(),
    }
}

fn metas(meta: &MetaInfo) -> String {
    meta.iter()
        .map(|(k, v)| format!("{k}={};", describe(v)))
        .collect()
}

/// The driver's `dumpFeature`, one line per feature.
fn dump(key: &str, feature: &Feature, lines: &mut Vec<(String, String)>) {
    let hulls: String = feature
        .convex_hulls
        .iter()
        .map(|hull| {
            let points: String = hull
                .hull_points()
                .iter()
                .map(|p| format!("{},{};", hexd(p.rt), hexd(p.mz)))
                .collect();
            format!("[{points}]")
        })
        .collect();
    lines.push((
        key.to_owned(),
        format!(
            "{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
            feature.unique_id,
            hexd(feature.rt),
            hexd(feature.mz),
            hexd(f64::from(feature.intensity)),
            feature.charge,
            hexd(f64::from(feature.quality)),
            hexd(f64::from(feature.quality_rt)),
            hexd(f64::from(feature.quality_mz)),
            hexd(f64::from(feature.width)),
            hulls,
            metas(&feature.metadata),
            feature.subordinates.len()
        ),
    ));
    for (i, sub) in feature.subordinates.iter().enumerate() {
        dump(&format!("{key}/{i}"), sub, lines);
    }
}

/// Compare a map with the driver's `dumpMap(tag, ...)` records, all of them.
fn assert_map_matches(tag: &str, map: &FeatureMap) {
    assert_eq!(
        map.features.len().to_string(),
        rec("features", tag),
        "{tag}"
    );
    let mut ours = Vec::new();
    for (i, feature) in map.features.iter().enumerate() {
        dump(&format!("{tag}/{i}"), feature, &mut ours);
    }
    let theirs: Vec<(String, String)> = records()
        .into_iter()
        .filter(|(k, key, _)| k == "feature" && key.starts_with(&format!("{tag}/")))
        .map(|(_, key, value)| (key, value))
        .collect();
    assert_eq!(ours, theirs, "{tag}");
}

fn hull(points: &[(f64, f64)]) -> ConvexHull2D {
    let points: Vec<Point2D> = points.iter().map(|&(x, y)| Point2D::new(x, y)).collect();
    let mut hull = ConvexHull2D::new();
    hull.set_hull_points(&points).unwrap();
    hull
}

fn feature(id: u64, rt: f64, mz: f64, intensity: f32) -> Feature {
    let mut f = Feature::new(rt, mz, intensity);
    f.unique_id = id;
    f
}

/// The driver's `constructed()` map.
fn constructed() -> FeatureMap {
    let mut f1 = feature(11, 100.5, 500.25, 1234.5);
    f1.charge = 2;
    f1.quality = 0.75;
    f1.quality_rt = 0.5;
    f1.quality_mz = 0.25;
    f1.set_width(3.5).unwrap();
    f1.convex_hulls = vec![
        hull(&[(99.0, 500.0), (102.0, 500.0), (102.0, 500.5), (99.0, 500.5)]),
        hull(&[(98.5, 501.0), (103.0, 501.2)]),
    ];
    f1.metadata
        .insert("FWHM".into(), MetaValue::try_from(2.5).unwrap());
    f1.metadata.insert("label".into(), "top".into());
    f1.metadata.insert(
        "score_list".into(),
        MetaValue::try_from(vec![0.5, 1.25]).unwrap(),
    );
    let mut s1 = feature(12, 100.0, 500.0, 10.0);
    s1.charge = 2;
    s1.metadata.insert("isotope".into(), 0i64.into());
    s1.subordinates = vec![feature(13, 100.1, 500.1, 5.0)];
    let mut s3 = feature(14, 101.0, 501.0, 7.0);
    s3.convex_hulls = vec![hull(&[(100.0, 501.0), (101.5, 501.0)])];
    f1.subordinates = vec![s1, s3];
    let mut f2 = feature(21, 200.0, 600.0, 0.0);
    f2.charge = -1;
    let mut f3 = feature(9_223_372_036_854_775_813, -1.0, 0.0, 1.0);
    f3.metadata.insert("ints".into(), vec![1i64, 2].into());
    FeatureMap {
        features: vec![f1, f2, f3],
        ..FeatureMap::default()
    }
}

/// Records `import FileConverter_9` and `feature FileConverter_9/...`: the
/// retained C++-written features table imports as the executed source
/// imports it.
#[test]
fn retained_features_table_imports_as_in_the_release_build() {
    let table = ParquetFile::read_table(data(
        "fixtures/topp/FileConverter_9_output.featureparquet/features.parquet",
    ))
    .unwrap();
    let mut map = FeatureMap::default();
    let warnings = Fma::import_features_from_arrow(&table, &mut map).unwrap();
    assert_eq!(rec("import", "FileConverter_9"), "1");
    assert!(warnings.is_empty());
    assert_map_matches("FileConverter_9", &map);
}

/// Records `export ...` and `roundtrip ...`: the port's features table
/// equals the one the executed source wrote for the same map, and both read
/// back to the source's reconstruction (with `FWHM` taken from `width`).
#[test]
fn exported_features_table_matches_the_release_build() {
    let map = constructed();
    assert_map_matches("constructed", &map);
    let table = Fma::export_features_to_arrow(&map).unwrap();
    assert_eq!(table.num_rows().to_string(), rec("export", "rows"));
    let dir = std::env::temp_dir().join(format!("openms-fma-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("fm_features.parquet");
    assert_eq!(rec("export", "written"), "1");
    aio::write_table_to_parquet(&table, &path, &ParquetWriteConfig::default()).unwrap();
    // The same file content as the executed source's, read the same way.
    assert_eq!(
        ParquetFile::read_table(&path).unwrap(),
        ParquetFile::read_table(data("fm_features.parquet")).unwrap()
    );
    assert_eq!(rec("roundtrip", "ok"), "1");
    for file in [path, data("fm_features.parquet")] {
        let mut back = FeatureMap::default();
        let warnings =
            Fma::import_features_from_arrow(&ParquetFile::read_table(&file).unwrap(), &mut back)
                .unwrap();
        assert!(warnings.is_empty());
        assert_map_matches("roundtrip", &back);
    }
}

/// Records `orphans ...` and `empty ok`: subordinates whose parent row is
/// missing become top-level features in row order, with the source's
/// warning; an empty table adds nothing.
#[test]
fn orphaned_subordinates_become_top_level_as_in_the_release_build() {
    let table = Fma::export_features_to_arrow(&constructed()).unwrap();
    let sliced: RecordBatch = table.slice(1, table.num_rows() - 1);
    let mut map = FeatureMap::default();
    let warnings = Fma::import_features_from_arrow(&sliced, &mut map).unwrap();
    assert_eq!(rec("orphans", "ok"), "1");
    assert_map_matches("orphans", &map);
    assert_eq!(
        warnings,
        vec![
            "FeatureMapArrowIO: Could not find parent feature with id 11 for feature 12. Will \
             be added as top-level."
                .to_owned(),
            "FeatureMapArrowIO: Could not find parent feature with id 11 for feature 14. Will \
             be added as top-level."
                .to_owned(),
        ]
    );
    let mut empty = FeatureMap::default();
    Fma::import_features_from_arrow(&table.slice(0, 0), &mut empty).unwrap();
    assert_eq!(format!("1|{}", empty.features.len()), rec("empty", "ok"));
}

/// Native: a table of the wrong shape is refused and leaves the map alone.
#[test]
fn a_wrong_table_is_refused_without_touching_the_map() {
    let mut map = constructed();
    let before = map.clone();
    let wrong = ParquetFile::read_table(data(
        "fixtures/topp/FileConverter_9_output.featureparquet/proteins.parquet",
    ))
    .unwrap();
    assert!(Fma::import_features_from_arrow(&wrong, &mut map).is_err());
    assert_eq!(map, before);
}
