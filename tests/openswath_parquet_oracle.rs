// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

#![cfg(feature = "arrow")]

//! The OpenSWATH Parquet readers and writers and `MSExperimentArrowExport`
//! against the executed C++ Release build (evidence tier 1).
//!
//! Every expected value is a record of `../oracle/F1b-openswath/driver.cpp`,
//! run on ibminode06 against the Linux x86_64 Release install
//! `openms4-release-bc9cc12-c19e494-174b576` (g++ 14.4 `-O2 -DNDEBUG`, Arrow
//! C++ 23.0.1), twice, with identical records and byte-identical output
//! files. The records are `tests/data/openswath_parquet/cpp_release/records.tsv`
//! (sha256 `91a8ca259376ba1573c2842d3494828a94dfae37bf39811a48ed5a8a0894e5ee`);
//! the Parquet files beside it were written by the executed library. Doubles
//! are compared through C's `%a`, so bit for bit. See
//! `tests/data/openswath_parquet_provenance.json`.

use arrow_array::cast::AsArray;
use arrow_array::{Array, RecordBatch};
use arrow_schema::DataType;
use openms::analysis::openswath::openswath_helper::compute_precursor_id;
use openms::analysis::openswath::peak_map_extractor::ExtractedPeakMap;
use openms::analysis::openswath::transition_experiment::{
    LightCompound, LightTargetedExperiment, LightTransition,
};
use openms::format::mobilogram_parquet_consumer::{
    MobilogramAnnotation, MobilogramParquetConsumer,
};
use openms::format::ms_chromatogram_parquet_consumer::MSChromatogramParquetConsumer;
use openms::format::msexperiment_arrow_export::{
    ArrowChromatogramExportConfig, ArrowExportFormat, ArrowSpectraExportConfig,
    MSExperimentArrowExport,
};
use openms::format::parquet_file::ParquetFile;
use openms::format::parquet_filter::ParquetFilter;
use openms::format::xic_parquet_file::{XicAnalyte, XicChromatogram, XicParquetFile, XicQuery};
use openms::format::xim_parquet_file::{XimMobilogram, XimParquetFile, XimQuery};
use openms::format::xipm_parquet_consumer::XIPMParquetConsumer;
use openms::format::xipm_parquet_file::{XipmParquetFile, XipmPeakMap};
use openms::{
    ChromatogramPeak, Error, MSChromatogram, MSExperiment, MSSpectrum, MobilityPeak1D, Mobilogram,
    Peak1D, Precursor,
};
use std::collections::BTreeMap;
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/openswath_parquet")
        .join(name)
}

fn release(name: &str) -> PathBuf {
    data("cpp_release").join(name)
}

fn scratch(test: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("openms-oswpqo-{}-{test}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The records, keyed by their leading field(s).
fn records() -> BTreeMap<String, String> {
    let text = std::fs::read_to_string(release("records.tsv")).unwrap();
    let mut map = BTreeMap::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        let (key, value) = match parts.len() {
            2 => (parts[0].to_owned(), parts[1].to_owned()),
            3 => (format!("{}\t{}", parts[0], parts[1]), parts[2].to_owned()),
            _ => panic!("bad record {line}"),
        };
        map.insert(key, value);
    }
    map
}

fn record(records: &BTreeMap<String, String>, key: &str) -> String {
    records
        .get(key)
        .unwrap_or_else(|| panic!("no record {key}"))
        .clone()
}

/// C `printf("%a")` as glibc prints a double.
fn hexd(v: f64) -> String {
    if v.is_nan() {
        return if v.is_sign_negative() { "-nan" } else { "nan" }.into();
    }
    let sign = if v.is_sign_negative() { "-" } else { "" };
    if v.is_infinite() {
        return format!("{sign}inf");
    }
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

fn opt(v: Option<i64>) -> String {
    v.map_or_else(|| "NULL".into(), |v| v.to_string())
}

fn optd(v: Option<f64>) -> String {
    v.map_or_else(|| "NULL".into(), hexd)
}

fn join<T: ToString>(v: &[T]) -> String {
    v.iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn dbls(v: &[f64]) -> String {
    v.iter().map(|&x| hexd(x)).collect::<Vec<_>>().join(",")
}

/// The driver's `xicRow`.
fn xic_row(c: &XicChromatogram) -> String {
    format!(
        "{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
        c.run_id,
        c.source_file,
        c.ms_level,
        opt(c.precursor_id),
        opt(c.transition_id),
        c.modified_sequence,
        opt(c.precursor_charge),
        opt(c.product_charge),
        opt(c.detecting_transition),
        opt(c.precursor_decoy),
        opt(c.product_decoy),
        opt(c.transition_ordinal),
        c.transition_type,
        c.annotation,
        dbls(&c.rt),
        dbls(&c.intensity)
    )
}

/// The driver's `ximRow`.
fn xim_row(m: &XimMobilogram) -> String {
    format!(
        "{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
        m.run_id,
        m.source_file,
        m.ms_level,
        m.mobilogram_type,
        opt(m.precursor_id),
        opt(m.transition_id),
        opt(m.feature_id),
        optd(m.feature_rt),
        m.modified_sequence,
        opt(m.precursor_charge),
        opt(m.product_charge),
        opt(m.detecting_transition),
        opt(m.precursor_decoy),
        opt(m.product_decoy),
        opt(m.transition_ordinal),
        m.transition_type,
        m.annotation,
        dbls(&m.mobility),
        dbls(&m.intensity)
    )
}

/// The driver's `xipmRow`.
fn xipm_row(p: &XipmPeakMap) -> String {
    format!(
        "{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
        p.run_id,
        p.source_file,
        p.ms_level,
        p.peakmap_type,
        opt(p.precursor_id),
        opt(p.transition_id),
        p.modified_sequence,
        opt(p.precursor_charge),
        opt(p.product_charge),
        opt(p.detecting_transition),
        opt(p.precursor_decoy),
        opt(p.product_decoy),
        opt(p.transition_ordinal),
        p.transition_type,
        p.annotation,
        hexd(p.target_mz),
        optd(p.target_rt),
        optd(p.target_ion_mobility),
        optd(p.rt_start),
        optd(p.rt_end),
        dbls(&p.mz),
        dbls(&p.rt),
        dbls(&p.ion_mobility),
        dbls(&p.intensity)
    )
}

/// The driver's `analytes`.
fn analytes(v: &[XicAnalyte]) -> String {
    let mut out = v.len().to_string();
    for a in v {
        out.push_str(&format!(
            " {{{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}}}",
            opt(a.precursor_id),
            a.modified_sequence,
            opt(a.precursor_charge),
            opt(a.precursor_decoy),
            opt(a.transition_id),
            opt(a.product_charge),
            opt(a.transition_ordinal),
            opt(a.detecting_transition),
            opt(a.product_decoy),
            a.transition_type,
            a.annotation,
            join(&a.transition_ids),
            join(&a.product_charges),
            join(&a.transition_ordinals),
            join(&a.detecting_transitions),
            join(&a.product_decoys),
            a.transition_types.join(","),
            a.annotations.join(",")
        ));
    }
    out
}

/// Compare a record with a Rust result: the record's text, or any error for
/// an `EXC` record (`InvalidValue` for the source's `InvalidValue` and
/// `ConversionError`).
fn check(key: &str, expected: &str, actual: openms::Result<String>) {
    if expected.starts_with("EXC ") {
        match actual {
            Err(Error::InvalidValue(_)) => {}
            other => panic!("{key}: expected {expected}, got {other:?}"),
        }
    } else {
        assert_eq!(
            actual.unwrap_or_else(|e| panic!("{key}: {e}")),
            expected,
            "{key}"
        );
    }
}

fn ids(rows: &[XicChromatogram]) -> String {
    format!(
        "{} {}",
        rows.len(),
        rows.iter()
            .map(|c| c
                .transition_id
                .map_or_else(|| "0".into(), |v| v.to_string()))
            .collect::<Vec<_>>()
            .join(",")
    )
}

#[test]
fn xic_filters_match_release() {
    let records = records();
    let xic = XicParquetFile::new(data("XICParquetFile_1_input.xic")).unwrap();
    let mut seen = 0;
    for (key, expected) in &records {
        let Some(filter) = key.strip_prefix("xic_filter\t") else {
            continue;
        };
        let query = XicQuery {
            filter: filter.into(),
            ..XicQuery::default()
        };
        check(
            key,
            expected.trim_end(),
            xic.chromatograms(&query)
                .map(|rows| ids(&rows).trim_end().to_owned()),
        );
        seen += 1;
    }
    assert_eq!(seen, 28);
    let run = |query: XicQuery| xic.chromatograms(&query).map(|r| r.len().to_string());
    check(
        "xic_args\tprecursor_id=3,ms_level=2",
        &record(&records, "xic_args\tprecursor_id=3,ms_level=2"),
        run(XicQuery {
            precursor_id: Some(3),
            ms_level: Some(2),
            ..XicQuery::default()
        }),
    );
    check(
        "xic_args\tprecursor_charge=3",
        &record(&records, "xic_args\tprecursor_charge=3"),
        run(XicQuery {
            precursor_charge: Some(3),
            ..XicQuery::default()
        }),
    );
    check(
        "xic_args\tmodified_sequence=PEPTIDEF,filter",
        &record(&records, "xic_args\tmodified_sequence=PEPTIDEF,filter"),
        xic.chromatograms(&XicQuery {
            modified_sequence: "PEPTIDEF".into(),
            filter: "transition_id>16".into(),
            ..XicQuery::default()
        })
        .map(|r| ids(&r)),
    );
    let typed = |filter: ParquetFilter| {
        xic.chromatograms_with_filter(&filter)
            .map(|r| r.len().to_string())
    };
    let mut f = ParquetFilter::new();
    f.eq_str("PRECURSOR_ID", "2");
    check(
        "xic_typed\teq_str PRECURSOR_ID 2",
        &record(&records, "xic_typed\teq_str PRECURSOR_ID 2"),
        typed(f),
    );
    let mut f = ParquetFilter::new();
    f.eq_int("MODIFIED_SEQUENCE", 5);
    check(
        "xic_typed\teq_int MODIFIED_SEQUENCE 5",
        &record(&records, "xic_typed\teq_int MODIFIED_SEQUENCE 5"),
        typed(f),
    );
    let mut f = ParquetFilter::new();
    f.in_int("TRANSITION_ID", &[]);
    check(
        "xic_typed\tin_int TRANSITION_ID empty",
        &record(&records, "xic_typed\tin_int TRANSITION_ID empty"),
        typed(f),
    );
    let mut f = ParquetFilter::new();
    f.eq_int("precursor_id", 4)
        .or_next()
        .eq_int("transition_id", 1);
    check(
        "xic_typed\tlowercase column precursor_id",
        &record(&records, "xic_typed\tlowercase column precursor_id"),
        xic.chromatograms_with_filter(&f).map(|r| ids(&r)),
    );
}

#[test]
fn xic_metadata_matches_release() {
    let records = records();
    let xic = XicParquetFile::new(data("XICParquetFile_1_input.xic")).unwrap();
    let mut seen = 0;
    for (key, expected) in &records {
        let Some(spec) = key.strip_prefix("xic_analytes\t") else {
            continue;
        };
        let (columns, mode) = spec.rsplit_once(' ').unwrap();
        let columns: Vec<&str> = columns.split(',').collect();
        check(
            key,
            expected,
            xic.analytes(&columns, mode == "nested")
                .map(|a| analytes(&a)),
        );
        seen += 1;
    }
    assert_eq!(seen, 12);
    assert_eq!(
        xic.columns().unwrap().join(","),
        record(&records, "xic_columns")
    );
    let runs = xic.runs().unwrap();
    assert_eq!(
        format!("{} {} {}", runs.len(), runs[0].run_id, runs[0].source_file),
        record(&records, "xic_runs")
    );
    assert_eq!(
        xic_row(&xic.load().unwrap()[0]),
        record(&records, "xic_load_row0")
    );
}

#[test]
fn xim_matches_release() {
    let records = records();
    let xim = XimParquetFile::new(data("XIMParquetFile_23_input.xim")).unwrap();
    let mut seen = 0;
    for (key, expected) in &records {
        let Some(filter) = key.strip_prefix("xim_filter\t") else {
            continue;
        };
        let query = XimQuery {
            filter: filter.into(),
            ..XimQuery::default()
        };
        check(
            key,
            expected,
            xim.mobilograms(&query).map(|r| r.len().to_string()),
        );
        seen += 1;
    }
    assert_eq!(seen, 6);
    check(
        "xim_args\tfeature_rt=25.2",
        &record(&records, "xim_args\tfeature_rt=25.2"),
        xim.mobilograms(&XimQuery {
            feature_rt: Some(25.200000476837157),
            ..XimQuery::default()
        })
        .map(|r| r.len().to_string()),
    );
    check(
        "xim_args\ttype=ms2,precursor=2",
        &record(&records, "xim_args\ttype=ms2,precursor=2"),
        xim.mobilograms(&XimQuery {
            precursor_id: Some(2),
            mobilogram_type: "ms2".into(),
            ..XimQuery::default()
        })
        .map(|r| r.len().to_string()),
    );
    let none: [&str; 0] = [];
    for nest in [false, true] {
        let mode = if nest { "nested" } else { "flat" };
        check(
            &format!("xim_analytes\tTRANSITION_ID {mode}"),
            &record(&records, &format!("xim_analytes\tTRANSITION_ID {mode}")),
            xim.analytes(&["TRANSITION_ID"], nest)
                .map(|a| a.len().to_string()),
        );
        assert_eq!(
            analytes(&xim.analytes(&none, nest).unwrap()),
            record(&records, &format!("xim_analytes\tdefault {mode}")),
            "{mode}"
        );
    }
    let rows = xim.load().unwrap();
    assert_eq!(
        format!("{} || {}", xim_row(&rows[5]), xim_row(&rows[0])),
        record(&records, "xim_load_row5")
    );
}

/// The driver's `library()`.
fn library() -> LightTargetedExperiment {
    let mut exp = LightTargetedExperiment::default();
    exp.compounds.push(LightCompound {
        id: "pep1".into(),
        sequence: "PEPTIDE".into(),
        charge: 2,
        ..LightCompound::default()
    });
    exp.compounds.push(LightCompound {
        id: "42".into(),
        sequence: "PEPTIDEK".into(),
        charge: 0,
        ..LightCompound::default()
    });
    let mut t = LightTransition {
        transition_name: "tr1".into(),
        peptide_ref: "pep1".into(),
        fragment_charge: 1,
        fragment_nr: 7,
        ..LightTransition::default()
    };
    t.set_fragment_type("y");
    t.set_detecting_transition(true);
    exp.transitions.push(t.clone());
    let mut t2 = t.clone();
    t2.transition_name = "7".into();
    t2.peptide_ref = "42".into();
    t2.fragment_charge = 0;
    t2.fragment_nr = -1;
    t2.set_fragment_type("");
    t2.set_decoy(true);
    exp.transitions.push(t2);
    let mut t3 = t;
    t3.transition_name = "-3".into();
    t3.peptide_ref = "orphan".into();
    t3.fragment_charge = 2;
    t3.fragment_nr = 4;
    t3.set_fragment_type("b");
    exp.transitions.push(t3);
    exp
}

fn chrom(id: &str, points: &[(f64, f32)]) -> MSChromatogram {
    MSChromatogram {
        native_id: id.into(),
        peaks: points
            .iter()
            .map(|&(rt, i)| ChromatogramPeak::new(rt, i))
            .collect(),
        ..MSChromatogram::default()
    }
}

/// The driver's `partXicConsumer` chromatograms.
fn consumer_chromatograms() -> Vec<MSChromatogram> {
    vec![
        chrom("tr1", &[(100.0, 1000.0)]),
        chrom("tr1", &[(100.0, 1000.0), (101.0, 900.0)]),
        chrom("tr1", &[(100.0, 1000.0), (101.0, 900.0), (102.5, 1234.5)]),
        chrom(
            "pep1_Precursor_i0",
            &[(100.0, 10.0), (101.0, 20.0), (102.0, 30.0)],
        ),
        chrom("7", &[(5.0, 1.0), (6.0, 2.0), (7.0, 3.0)]),
        chrom("-3", &[(-5.0, 1.0), (-4.0, 2.0), (-3.0, 3.0)]),
        chrom("tr1", &[]),
        chrom("42_Precursor_i1", &[(1.0, 0.0), (2.0, 0.0), (3.0, 0.0)]),
        chrom("orphan_Precursor_i0", &[(1.0, 5.0), (2.0, 6.0), (3.0, 7.0)]),
    ]
}

fn rows_of(record: &str) -> Vec<String> {
    record.split(" || ").skip(1).map(str::to_owned).collect()
}

/// The source writer's `.xic`, read here, is the source reader's rows.
#[test]
fn cpp_written_xic_reads_as_release() {
    let records = records();
    let expected = record(&records, "xicw_read");
    let rows = XicParquetFile::new(release("consumer.xic"))
        .unwrap()
        .load()
        .unwrap();
    assert_eq!(
        rows.len().to_string(),
        expected.split(" || ").next().unwrap()
    );
    let actual: Vec<String> = rows.iter().map(xic_row).collect();
    assert_eq!(actual, rows_of(&expected));
    // One and two points: the source stored RT with a fixed point of 0.
    assert!(rows[0].rt[0].is_nan() && rows[1].rt.iter().all(|v| v.is_nan()));
    // Negative RTs wrapped through the unsigned 32-bit Numpress integers.
    assert!(rows[5].rt[0] > 4.0e8);
    // The chromatogram consumed after finalize() was dropped.
    assert_eq!(
        XicParquetFile::new(release("consumer_after.xic"))
            .unwrap()
            .load()
            .unwrap()
            .len(),
        1
    );
    assert_eq!(record(&records, "xicw_after_finalize_read"), "1");
    // One unknown native id cost the source every chromatogram.
    assert!(record(&records, "xicw_badid").contains("expected length 2 but got length 1"));
    assert_eq!(
        XicParquetFile::new(release("consumer_badid.xic"))
            .unwrap()
            .load()
            .unwrap()
            .len(),
        0
    );
}

/// The Rust writer on the driver's input: every row the source's, except the
/// arrays the source cannot represent (native: stored losslessly).
#[test]
fn rust_written_xic_matches_release() {
    let records = records();
    let expected = rows_of(&record(&records, "xicw_read"));
    let out = scratch("xicw").join("consumer.xic");
    let mut consumer = MSChromatogramParquetConsumer::new(&out, 11, "src.mzML", &library());
    for mut c in consumer_chromatograms() {
        consumer.consume(&mut c).unwrap();
    }
    consumer.finalize().unwrap();
    let rows = XicParquetFile::new(&out).unwrap().load().unwrap();
    assert_eq!(rows.len(), expected.len());
    for (i, (row, expected)) in rows.iter().zip(&expected).enumerate() {
        let native = [0usize, 1, 5].contains(&i);
        if native {
            let input: Vec<f64> = consumer_chromatograms()[i]
                .peaks
                .iter()
                .map(|p| p.rt)
                .collect();
            assert_eq!(row.rt, input, "row {i}");
            let mut patched = row.clone();
            patched.rt = Vec::new();
            let expected_fields: Vec<&str> = expected.split('|').collect();
            let actual = xic_row(&patched);
            let actual_fields: Vec<&str> = actual.split('|').collect();
            for (k, (a, e)) in actual_fields.iter().zip(&expected_fields).enumerate() {
                if k != 14 {
                    assert_eq!(a, e, "row {i} field {k}");
                }
            }
        } else {
            assert_eq!(&xic_row(row), expected, "row {i}");
        }
    }
    let mut after = chrom("tr1", &[(1.0, 2.0), (2.0, 3.0), (3.0, 4.0)]);
    assert!(consumer.consume(&mut after).is_err());
}

fn mobilogram(peaks: &[(f64, f32)]) -> Mobilogram {
    Mobilogram::from_peaks(
        peaks
            .iter()
            .map(|&(m, i)| MobilityPeak1D::new(m, i))
            .collect(),
    )
}

/// The source writer's `.xim` read here, and the Rust writer on the same
/// input, both give the source reader's rows.
#[test]
fn xim_writer_and_reader_match_release() {
    let records = records();
    let expected = rows_of(&record(&records, "ximw_read"));
    let cpp: Vec<String> = XimParquetFile::new(release("consumer.xim"))
        .unwrap()
        .load()
        .unwrap()
        .iter()
        .map(xim_row)
        .collect();
    assert_eq!(cpp, expected);
    let out = scratch("ximw").join("consumer.xim");
    let mut consumer =
        MobilogramParquetConsumer::new(&out, (1u64 << 63) | 5, "src.mzML", &library());
    let one = mobilogram(&[(0.9, 100.0)]);
    let three = mobilogram(&[(0.8, 1000.0), (0.9, 2000.0), (1.05, 3000.0)]);
    let empty = Mobilogram::default();
    let annotation =
        |kind: &str, level: i64, transition: i64, native: &str, rt: f64, feature: i64| {
            MobilogramAnnotation {
                mobilogram_type: kind.into(),
                ms_level: Some(level),
                transition_id: Some(transition),
                transition_native_id: native.into(),
                feature_rt: Some(rt),
                feature_id: Some(feature),
            }
        };
    consumer
        .consume_mobilogram(&one, &MobilogramAnnotation::default())
        .unwrap();
    consumer
        .consume_mobilogram(&three, &annotation("ms2", 2, -1, "tr1", 50.0, 7))
        .unwrap();
    consumer
        .consume_mobilogram(&three, &annotation("ms2", 2, 99, "tr1", f64::NAN, -1))
        .unwrap();
    consumer
        .consume_mobilogram(&three, &annotation("ms1", 1, 5, "unknown", 12.5, 3))
        .unwrap();
    consumer
        .consume_mobilogram(&empty, &annotation("", -1, -1, "7", -1.0, -1))
        .unwrap();
    consumer
        .consume_mobilogram(&three, &annotation("ms2", 2, -1, "-3", 1.0, 1))
        .unwrap();
    consumer.finalize().unwrap();
    assert!(record(&records, "ximw_write").contains("EXC IllegalArgument"));
    assert!(
        consumer
            .consume_mobilogram(&one, &MobilogramAnnotation::default())
            .is_err()
    );
    let rust: Vec<String> = XimParquetFile::new(&out)
        .unwrap()
        .load()
        .unwrap()
        .iter()
        .map(xim_row)
        .collect();
    assert_eq!(rust, expected);
}

fn peak_map(
    id: &str,
    target: (f64, f64, f64),
    window: (f64, f64),
    arrays: [&[f64]; 4],
) -> ExtractedPeakMap {
    ExtractedPeakMap {
        native_id: id.into(),
        target_mz: target.0,
        target_rt: target.1,
        target_ion_mobility: target.2,
        rt_start: window.0,
        rt_end: window.1,
        mz: arrays[0].to_vec(),
        rt: arrays[1].to_vec(),
        ion_mobility: arrays[2].to_vec(),
        intensity: arrays[3].to_vec(),
    }
}

/// The source writer's `.xipm` read here, and the Rust writer on the same
/// input, both give the source reader's rows.
#[test]
fn xipm_writer_and_reader_match_release() {
    let records = records();
    let expected = rows_of(&record(&records, "xipmw_read"));
    let cpp: Vec<String> = XipmParquetFile::new(release("consumer.xipm"))
        .unwrap()
        .load()
        .unwrap()
        .iter()
        .map(xipm_row)
        .collect();
    assert_eq!(cpp, expected);
    let out = scratch("xipmw").join("consumer.xipm");
    let mut consumer = XIPMParquetConsumer::new(&out, &library());
    let nan = f64::NAN;
    consumer
        .consume_peak_map(
            &peak_map(
                "tr1",
                (500.2, 100.0, 1.1),
                (95.0, 105.0),
                [
                    &[500.19, 500.20],
                    &[100.0, 101.0],
                    &[1.05, 1.08],
                    &[1000.0, 900.0],
                ],
            ),
            7,
            "run1.mzML",
            2,
        )
        .unwrap();
    consumer
        .consume_peak_map(
            &peak_map(
                &compute_precursor_id("pep1", 0),
                (600.2, 100.0, 1.1),
                (95.0, 105.0),
                [&[600.19], &[100.0], &[1.06], &[500.0]],
            ),
            7,
            "run1.mzML",
            1,
        )
        .unwrap();
    consumer
        .consume_peak_map(
            &peak_map(
                "7",
                (300.0, nan, -1.0),
                (10.0, 20.0),
                [
                    &[299.9, 300.0, 300.1],
                    &[10.0, 15.0, 20.0],
                    &[0.7, 0.8, 0.9],
                    &[1.0, 2.0, 3.0],
                ],
            ),
            u64::MAX,
            "",
            2,
        )
        .unwrap();
    consumer
        .consume_peak_map(
            &peak_map("-3", (250.0, nan, 0.0), (0.0, -1.0), [&[], &[], &[], &[]]),
            1,
            "x",
            2,
        )
        .unwrap();
    assert!(record(&records, "xipmw_write").contains("ragged-EXC IllegalArgument"));
    assert!(
        consumer
            .consume_peak_map(
                &peak_map(
                    "tr1",
                    (1.0, 1.0, 1.0),
                    (0.0, 1.0),
                    [&[1.0], &[1.0, 2.0], &[1.0], &[1.0]]
                ),
                1,
                "x",
                2
            )
            .is_err()
    );
    consumer.finalize().unwrap();
    let rust: Vec<String> = XipmParquetFile::new(&out)
        .unwrap()
        .load()
        .unwrap()
        .iter()
        .map(xipm_row)
        .collect();
    assert_eq!(rust, expected);
}

/// The driver's `testExperiment()`.
fn arrow_experiment() -> MSExperiment {
    let spectrum = |rt: f64, level: u32, id: &str, peaks: &[(f64, f32)]| MSSpectrum {
        rt,
        ms_level: level,
        native_id: id.into(),
        peaks: peaks.iter().map(|&(mz, i)| Peak1D::new(mz, i)).collect(),
        ..MSSpectrum::default()
    };
    let mut exp = MSExperiment::default();
    exp.spectra.push(spectrum(
        100.0,
        1,
        "spectrum=0",
        &[(100.0, 1000.0), (200.0, 2000.0), (300.0, 3000.0)],
    ));
    let mut ms2 = spectrum(
        110.123456789,
        2,
        "spectrum=1",
        &[(150.0, 500.0), (250.0, 750.0)],
    );
    let mut precursor = Precursor::new(500.0, 2);
    precursor.intensity = 50000.0;
    precursor.isolation_window_lower_offset = 1.5;
    precursor.isolation_window_upper_offset = 2.5;
    ms2.precursors = vec![precursor];
    exp.spectra.push(ms2);
    exp.spectra.push(spectrum(
        200.0,
        1,
        "spectrum=2",
        &[(150.0, 1500.0), (250.0, 2500.0)],
    ));
    exp.spectra.push(spectrum(250.0, 3, "spectrum=3", &[]));
    let chromatogram = |id: &str, precursor: f64, product: f64, points: &[(f64, f32)]| {
        let mut c = chrom(id, points);
        c.precursor.mz = precursor;
        c.product.mz = product;
        c
    };
    exp.chromatograms.push(chromatogram(
        "TIC",
        500.0,
        200.0,
        &[(10.0, 100.0), (20.0, 200.0), (30.0, 150.0)],
    ));
    exp.chromatograms.push(chromatogram(
        "XIC_1",
        600.0,
        300.0,
        &[(15.0, 50.0), (25.0, 75.0)],
    ));
    exp
}

/// Column-by-column equality of two tables; list columns compared value by
/// value, so the list element field's name (`item` or Parquet's `element`)
/// does not matter.
fn assert_tables_equal(name: &str, cpp: &RecordBatch, rust: &RecordBatch) {
    let names = |t: &RecordBatch| {
        t.schema()
            .fields()
            .iter()
            .map(|f| f.name().clone())
            .collect::<Vec<_>>()
    };
    assert_eq!(names(cpp), names(rust), "{name}");
    assert_eq!(cpp.num_rows(), rust.num_rows(), "{name}");
    for (i, (a, b)) in cpp.columns().iter().zip(rust.columns()).enumerate() {
        match (a.data_type(), b.data_type()) {
            (DataType::List(_), DataType::List(_)) => {
                let (a, b) = (a.as_list::<i32>(), b.as_list::<i32>());
                for row in 0..a.len() {
                    assert_eq!(
                        a.is_null(row),
                        b.is_null(row),
                        "{name} column {i} row {row}"
                    );
                    assert_eq!(
                        a.value(row).as_ref(),
                        b.value(row).as_ref(),
                        "{name} column {i} row {row}"
                    );
                }
            }
            (x, y) => {
                assert_eq!(x, y, "{name} column {i}");
                assert_eq!(a.as_ref(), b.as_ref(), "{name} column {i}");
            }
        }
    }
}

#[test]
fn arrow_export_matches_release() {
    let records = records();
    let exp = arrow_experiment();
    let long = ArrowSpectraExportConfig::default();
    let wide = ArrowSpectraExportConfig {
        format: ArrowExportFormat::SemiWide,
        ..ArrowSpectraExportConfig::default()
    };
    let cases = [
        ("spectra_long", long.clone()),
        ("spectra_wide", wide.clone()),
        (
            "spectra_long_windows",
            ArrowSpectraExportConfig {
                min_mz: 150.0,
                max_mz: 250.0,
                min_rt: 105.0,
                ..long.clone()
            },
        ),
        (
            "spectra_wide_windows",
            ArrowSpectraExportConfig {
                min_mz: 150.0,
                max_mz: 250.0,
                ms_levels: vec![1, 3],
                ..wide
            },
        ),
        (
            "spectra_long_columns",
            ArrowSpectraExportConfig {
                columns: vec![
                    "native_id".into(),
                    "precursor_charge".into(),
                    "bogus".into(),
                    "mz".into(),
                ],
                ..long
            },
        ),
    ];
    for (name, config) in &cases {
        assert_eq!(
            MSExperimentArrowExport::spectra_arrow_column_names(&exp, config)
                .unwrap()
                .join(","),
            record(&records, &format!("arrow_names\t{name}")),
            "{name}"
        );
        assert_eq!(record(&records, &format!("arrow_parquet\t{name}")), "true");
        let cpp = ParquetFile::read_table(release(&format!("{name}.parquet"))).unwrap();
        let rust = MSExperimentArrowExport::spectra_to_arrow(&exp, config).unwrap();
        assert_tables_equal(name, &cpp, &rust);
    }
    let chrom_cases = [
        ("chroms_long", ArrowChromatogramExportConfig::default()),
        (
            "chroms_wide_rt",
            ArrowChromatogramExportConfig {
                format: ArrowExportFormat::SemiWide,
                min_rt: 12.0,
                max_rt: 22.0,
                ..ArrowChromatogramExportConfig::default()
            },
        ),
    ];
    for (name, config) in &chrom_cases {
        assert_eq!(
            MSExperimentArrowExport::chromatogram_arrow_column_names(&exp, config).join(","),
            record(&records, &format!("arrow_names\t{name}")),
            "{name}"
        );
        let cpp = ParquetFile::read_table(release(&format!("{name}.parquet"))).unwrap();
        let rust = MSExperimentArrowExport::chromatograms_to_arrow(&exp, config).unwrap();
        assert_tables_equal(name, &cpp, &rust);
    }
}
