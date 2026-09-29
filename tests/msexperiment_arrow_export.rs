// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

#![cfg(feature = "arrow")]

//! `MSExperimentArrowExport` against `MSExperimentArrowExport_test.cpp` at
//! core `bc9cc12`: every section, with its literals (evidence tier 3). The
//! class test reaches the tables through the C Data Interface; here the
//! `*_to_arrow` functions return them directly. Tests marked "invariant"
//! check values that follow from the input by the header's column
//! definitions (independently derived), not from Rust output.

use arrow_array::cast::AsArray;
use arrow_array::types::{Float32Type, Float64Type, Int16Type, UInt8Type, UInt32Type};
use arrow_array::{Array, RecordBatch};
use arrow_schema::DataType;
use openms::format::msexperiment_arrow_export::{
    ArrowChromatogramExportConfig, ArrowExportFormat, ArrowSpectraExportConfig,
    MSExperimentArrowExport, ParquetCompression, ParquetWriteConfig,
};
use openms::format::parquet_file::ParquetFile;
use openms::{ChromatogramPeak, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, Precursor};
use std::path::PathBuf;

fn scratch(test: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("openms-msexparrow-{}-{test}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn spectrum(rt: f64, level: u32, native_id: &str, peaks: &[(f64, f32)]) -> MSSpectrum {
    MSSpectrum {
        rt,
        ms_level: level,
        native_id: native_id.into(),
        peaks: peaks.iter().map(|&(mz, i)| Peak1D::new(mz, i)).collect(),
        ..MSSpectrum::default()
    }
}

/// `createTestExperiment()`.
fn test_experiment() -> MSExperiment {
    let mut exp = MSExperiment::default();
    exp.spectra.push(spectrum(
        100.0,
        1,
        "spectrum=0",
        &[(100.0, 1000.0), (200.0, 2000.0), (300.0, 3000.0)],
    ));
    let mut ms2 = spectrum(110.0, 2, "spectrum=1", &[(150.0, 500.0), (250.0, 750.0)]);
    let mut precursor = Precursor::new(500.0, 2);
    precursor.intensity = 50000.0;
    precursor.isolation_window_lower_offset = 1.5;
    precursor.isolation_window_upper_offset = 1.5;
    ms2.precursors = vec![precursor];
    exp.spectra.push(ms2);
    exp.spectra.push(spectrum(
        200.0,
        1,
        "spectrum=2",
        &[(150.0, 1500.0), (250.0, 2500.0)],
    ));
    exp
}

fn chromatogram(
    native_id: &str,
    precursor: f64,
    product: f64,
    points: &[(f64, f32)],
) -> MSChromatogram {
    let mut chrom = MSChromatogram {
        native_id: native_id.into(),
        peaks: points
            .iter()
            .map(|&(rt, i)| ChromatogramPeak::new(rt, i))
            .collect(),
        ..MSChromatogram::default()
    };
    chrom.precursor.mz = precursor;
    chrom.product.mz = product;
    chrom
}

fn long() -> ArrowSpectraExportConfig {
    ArrowSpectraExportConfig::default()
}

fn has(table: &RecordBatch, name: &str) -> bool {
    table.schema().index_of(name).is_ok()
}

fn field_type(table: &RecordBatch, name: &str) -> DataType {
    table
        .schema()
        .field_with_name(name)
        .unwrap()
        .data_type()
        .clone()
}

/// `exportSpectraToArrowCDataInterface - empty experiment`.
#[test]
fn spectra_empty_experiment() {
    let table =
        MSExperimentArrowExport::spectra_to_arrow(&MSExperiment::default(), &long()).unwrap();
    assert_eq!(table.num_rows(), 0);
    assert!(table.num_columns() > 0);
}

/// `exportSpectraToArrowCDataInterface - long format basic`.
#[test]
fn spectra_long_basic() {
    let table = MSExperimentArrowExport::spectra_to_arrow(&test_experiment(), &long()).unwrap();
    assert_eq!(table.num_rows(), 7);
    for name in [
        "mz",
        "intensity",
        "rt",
        "spectrum_index",
        "ms_level",
        "native_id",
        "precursor_mz",
        "precursor_charge",
    ] {
        assert!(has(&table, name), "{name}");
    }
    assert_eq!(field_type(&table, "mz"), DataType::Float64);
    assert_eq!(field_type(&table, "intensity"), DataType::Float32);
    assert_eq!(field_type(&table, "ms_level"), DataType::UInt8);
}

/// Invariant: the long table's values are the peaks with their spectrum's
/// metadata; precursor columns are null for the MS1 spectra.
#[test]
fn spectra_long_values() {
    let table = MSExperimentArrowExport::spectra_to_arrow(&test_experiment(), &long()).unwrap();
    let names: Vec<String> = table
        .schema()
        .fields()
        .iter()
        .map(|f| f.name().clone())
        .collect();
    assert_eq!(
        names,
        [
            "mz",
            "intensity",
            "rt",
            "spectrum_index",
            "ms_level",
            "native_id",
            "precursor_mz",
            "precursor_charge",
            "precursor_intensity",
            "isolation_lower",
            "isolation_upper"
        ]
    );
    let col = |name: &str| table.column(table.schema().index_of(name).unwrap()).clone();
    let mz = col("mz");
    let mz = mz.as_primitive::<Float64Type>();
    assert_eq!(
        mz.values().to_vec(),
        [100.0, 200.0, 300.0, 150.0, 250.0, 150.0, 250.0]
    );
    let index = col("spectrum_index");
    assert_eq!(
        index.as_primitive::<UInt32Type>().values().to_vec(),
        [0, 0, 0, 1, 1, 2, 2]
    );
    let level = col("ms_level");
    assert_eq!(
        level.as_primitive::<UInt8Type>().values().to_vec(),
        [1, 1, 1, 2, 2, 1, 1]
    );
    let rt = col("rt");
    assert_eq!(rt.as_primitive::<Float32Type>().value(3), 110.0f32);
    let prec = col("precursor_mz");
    let prec = prec.as_primitive::<Float64Type>();
    assert_eq!(prec.null_count(), 5);
    assert_eq!(prec.value(3), 500.0);
    let charge = col("precursor_charge");
    assert_eq!(charge.as_primitive::<Int16Type>().value(4), 2);
    let pint = col("precursor_intensity");
    assert_eq!(pint.as_primitive::<Float32Type>().value(3), 50000.0);
    let lower = col("isolation_lower");
    assert_eq!(lower.as_primitive::<Float64Type>().value(3), 1.5);
    let native = col("native_id");
    assert_eq!(native.as_string::<i32>().value(6), "spectrum=2");
}

/// `exportSpectraToArrowCDataInterface - long format with MS level filter`.
#[test]
fn spectra_long_ms_level_filter() {
    let config = ArrowSpectraExportConfig {
        ms_levels: vec![1],
        ..long()
    };
    assert_eq!(
        MSExperimentArrowExport::spectra_to_arrow(&test_experiment(), &config)
            .unwrap()
            .num_rows(),
        5
    );
}

/// `exportSpectraToArrowCDataInterface - long format with RT filter`.
#[test]
fn spectra_long_rt_filter() {
    let config = ArrowSpectraExportConfig {
        min_rt: 100.0,
        max_rt: 110.0,
        ..long()
    };
    assert_eq!(
        MSExperimentArrowExport::spectra_to_arrow(&test_experiment(), &config)
            .unwrap()
            .num_rows(),
        5
    );
}

/// `exportSpectraToArrowCDataInterface - long format with m/z filter`.
#[test]
fn spectra_long_mz_filter() {
    let config = ArrowSpectraExportConfig {
        min_mz: 150.0,
        max_mz: 250.0,
        ..long()
    };
    assert_eq!(
        MSExperimentArrowExport::spectra_to_arrow(&test_experiment(), &config)
            .unwrap()
            .num_rows(),
        5
    );
}

/// `exportSpectraToArrowCDataInterface - long format with column selection`.
#[test]
fn spectra_long_column_selection() {
    let config = ArrowSpectraExportConfig {
        columns: vec!["mz".into(), "intensity".into(), "rt".into()],
        ..long()
    };
    let table = MSExperimentArrowExport::spectra_to_arrow(&test_experiment(), &config).unwrap();
    assert_eq!(table.num_columns(), 3);
    assert!(has(&table, "mz") && has(&table, "intensity") && has(&table, "rt"));
    assert!(!has(&table, "spectrum_index"));
}

/// `exportSpectraToArrowCDataInterface - long format no precursor info`.
#[test]
fn spectra_long_no_precursor_info() {
    let config = ArrowSpectraExportConfig {
        include_precursor_info: false,
        ..long()
    };
    let table = MSExperimentArrowExport::spectra_to_arrow(&test_experiment(), &config).unwrap();
    assert!(!has(&table, "precursor_mz"));
    assert!(!has(&table, "precursor_charge"));
    assert!(!has(&table, "isolation_lower"));
}

/// `exportSpectraToArrowCDataInterface - semi-wide format basic`.
#[test]
fn spectra_semi_wide_basic() {
    let config = ArrowSpectraExportConfig {
        format: ArrowExportFormat::SemiWide,
        ..long()
    };
    let table = MSExperimentArrowExport::spectra_to_arrow(&test_experiment(), &config).unwrap();
    assert_eq!(table.num_rows(), 3);
    assert!(matches!(field_type(&table, "mz"), DataType::List(_)));
    assert!(matches!(field_type(&table, "intensity"), DataType::List(_)));
    assert_eq!(field_type(&table, "rt"), DataType::Float32);
    // Invariant: the second row's list holds the MS2 peaks.
    let mz = table
        .column(table.schema().index_of("mz").unwrap())
        .as_list::<i32>()
        .value(1);
    assert_eq!(
        mz.as_primitive::<Float64Type>().values().to_vec(),
        [150.0, 250.0]
    );
}

/// `exportSpectraToArrowCDataInterface - semi-wide format with MS level filter`.
#[test]
fn spectra_semi_wide_ms_level_filter() {
    let config = ArrowSpectraExportConfig {
        format: ArrowExportFormat::SemiWide,
        ms_levels: vec![2],
        ..long()
    };
    assert_eq!(
        MSExperimentArrowExport::spectra_to_arrow(&test_experiment(), &config)
            .unwrap()
            .num_rows(),
        1
    );
}

/// `getSpectraArrowColumnNames - long format` and `- with column filter`.
#[test]
fn spectra_column_names() {
    let columns =
        MSExperimentArrowExport::spectra_arrow_column_names(&test_experiment(), &long()).unwrap();
    for name in [
        "mz",
        "intensity",
        "rt",
        "spectrum_index",
        "ms_level",
        "precursor_mz",
    ] {
        assert!(columns.iter().any(|c| c == name), "{name}");
    }
    let config = ArrowSpectraExportConfig {
        columns: vec!["mz".into(), "intensity".into()],
        ..long()
    };
    let columns =
        MSExperimentArrowExport::spectra_arrow_column_names(&test_experiment(), &config).unwrap();
    assert_eq!(columns, ["mz", "intensity"]);
}

/// `exportChromatogramsToArrowCDataInterface - empty chromatograms`.
#[test]
fn chromatograms_empty() {
    let table = MSExperimentArrowExport::chromatograms_to_arrow(
        &MSExperiment::default(),
        &ArrowChromatogramExportConfig::default(),
    )
    .unwrap();
    assert_eq!(table.num_rows(), 0);
}

fn two_chromatograms() -> MSExperiment {
    let mut exp = MSExperiment::default();
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

/// `exportChromatogramsToArrowCDataInterface - long format`.
#[test]
fn chromatograms_long() {
    let table = MSExperimentArrowExport::chromatograms_to_arrow(
        &two_chromatograms(),
        &ArrowChromatogramExportConfig::default(),
    )
    .unwrap();
    assert_eq!(table.num_rows(), 5);
    for name in [
        "rt",
        "intensity",
        "chromatogram_index",
        "precursor_mz",
        "product_mz",
    ] {
        assert!(has(&table, name), "{name}");
    }
    // Invariant: per-point chromatogram metadata.
    let product = table.column(table.schema().index_of("product_mz").unwrap());
    assert_eq!(
        product.as_primitive::<Float64Type>().values().to_vec(),
        [200.0, 200.0, 200.0, 300.0, 300.0]
    );
}

/// `exportChromatogramsToArrowCDataInterface - semi-wide format`.
#[test]
fn chromatograms_semi_wide() {
    let mut exp = MSExperiment::default();
    exp.chromatograms.push(chromatogram(
        "TIC",
        0.0,
        0.0,
        &[(10.0, 100.0), (20.0, 200.0)],
    ));
    let config = ArrowChromatogramExportConfig {
        format: ArrowExportFormat::SemiWide,
        ..ArrowChromatogramExportConfig::default()
    };
    let table = MSExperimentArrowExport::chromatograms_to_arrow(&exp, &config).unwrap();
    assert_eq!(table.num_rows(), 1);
    assert!(matches!(field_type(&table, "rt"), DataType::List(_)));
}

/// `exportChromatogramsToArrowCDataInterface - with RT filter`.
#[test]
fn chromatograms_rt_filter() {
    let mut exp = MSExperiment::default();
    exp.chromatograms.push(chromatogram(
        "TIC",
        0.0,
        0.0,
        &[(10.0, 100.0), (20.0, 200.0), (30.0, 150.0)],
    ));
    let config = ArrowChromatogramExportConfig {
        min_rt: 15.0,
        max_rt: 25.0,
        ..ArrowChromatogramExportConfig::default()
    };
    assert_eq!(
        MSExperimentArrowExport::chromatograms_to_arrow(&exp, &config)
            .unwrap()
            .num_rows(),
        1
    );
}

/// `[EXTRA] a failed write leaves no partial .parquet behind`.
#[test]
fn failed_write_leaves_no_file() {
    let mut exp = MSExperiment::default();
    exp.spectra.push(spectrum(1.0, 1, "", &[(100.0, 200.0)]));
    let dir = scratch("failed");
    let good = dir.join("good.parquet");
    MSExperimentArrowExport::export_spectra_to_parquet(
        &exp,
        &good,
        &long(),
        &ParquetWriteConfig::default(),
    )
    .unwrap();
    assert!(good.exists());
    let failed = dir.join("failed.parquet");
    let no_row_group = ParquetWriteConfig {
        row_group_size: 0,
        ..ParquetWriteConfig::default()
    };
    assert!(
        MSExperimentArrowExport::export_spectra_to_parquet(&exp, &failed, &long(), &no_row_group)
            .is_err()
    );
    assert!(!failed.exists());
}

/// `ArrowSpectraExportConfig`, `ArrowChromatogramExportConfig` and
/// `ParquetWriteConfig - default values`.
#[test]
fn config_defaults() {
    let config = ArrowSpectraExportConfig::default();
    assert_eq!(config.format, ArrowExportFormat::Long);
    assert!(config.ms_levels.is_empty());
    assert_eq!(
        (config.min_rt, config.max_rt, config.min_mz, config.max_mz),
        (0.0, 0.0, 0.0, 0.0)
    );
    assert!(config.columns.is_empty());
    assert!(config.include_precursor_info);
    assert!(config.include_ion_mobility);
    let config = ArrowChromatogramExportConfig::default();
    assert_eq!(config.format, ArrowExportFormat::Long);
    assert_eq!((config.min_rt, config.max_rt), (0.0, 0.0));
    assert!(config.columns.is_empty());
    let config = ParquetWriteConfig::default();
    assert_eq!(config.compression, ParquetCompression::Zstd);
    assert_eq!(config.compression_level, 3);
    assert_eq!(config.row_group_size, 128 * 1024 * 1024);
    assert!(config.write_statistics);
    assert_eq!(config.data_page_size, 1024 * 1024);
}

/// `exportSpectraToParquet - basic export`, `- with compression options`,
/// `- with filtering`, `- empty experiment`.
#[test]
fn spectra_to_parquet() {
    let dir = scratch("spectra");
    let basic = dir.join("basic.parquet");
    MSExperimentArrowExport::export_spectra_to_parquet(
        &test_experiment(),
        &basic,
        &long(),
        &ParquetWriteConfig::default(),
    )
    .unwrap();
    assert!(std::fs::metadata(&basic).unwrap().len() > 0);
    // Invariant: the file holds the table.
    assert_eq!(ParquetFile::read_table(&basic).unwrap().num_rows(), 7);
    let compressed = dir.join("zstd9.parquet");
    let config = ParquetWriteConfig {
        compression: ParquetCompression::Zstd,
        compression_level: 9,
        row_group_size: 1024 * 1024,
        ..ParquetWriteConfig::default()
    };
    MSExperimentArrowExport::export_spectra_to_parquet(
        &test_experiment(),
        &compressed,
        &long(),
        &config,
    )
    .unwrap();
    assert!(compressed.exists());
    let ms2 = dir.join("ms2.parquet");
    let filter = ArrowSpectraExportConfig {
        ms_levels: vec![2],
        ..long()
    };
    MSExperimentArrowExport::export_spectra_to_parquet(
        &test_experiment(),
        &ms2,
        &filter,
        &ParquetWriteConfig::default(),
    )
    .unwrap();
    assert!(ms2.exists());
    let empty = dir.join("empty.parquet");
    MSExperimentArrowExport::export_spectra_to_parquet(
        &MSExperiment::default(),
        &empty,
        &long(),
        &ParquetWriteConfig::default(),
    )
    .unwrap();
    assert!(empty.exists());
}

/// `exportChromatogramsToParquet - basic export` and `- with compression options`.
#[test]
fn chromatograms_to_parquet() {
    let dir = scratch("chroms");
    let basic = dir.join("basic.parquet");
    MSExperimentArrowExport::export_chromatograms_to_parquet(
        &two_chromatograms(),
        &basic,
        &ArrowChromatogramExportConfig::default(),
        &ParquetWriteConfig::default(),
    )
    .unwrap();
    assert!(std::fs::metadata(&basic).unwrap().len() > 0);
    let mut exp = MSExperiment::default();
    exp.chromatograms.push(chromatogram(
        "SRM_transition_1",
        450.0,
        175.0,
        &[(5.0, 500.0), (10.0, 1000.0), (15.0, 750.0)],
    ));
    let snappy = dir.join("snappy.parquet");
    let config = ParquetWriteConfig {
        compression: ParquetCompression::Snappy,
        ..ParquetWriteConfig::default()
    };
    MSExperimentArrowExport::export_chromatograms_to_parquet(
        &exp,
        &snappy,
        &ArrowChromatogramExportConfig::default(),
        &config,
    )
    .unwrap();
    assert!(snappy.exists());
}

/// Native: an RT window on unsorted spectra is refused (the source's
/// binary search assumes sorted spectra), and an inverted window selects
/// nothing (the source's iterator loop runs past the end).
#[test]
fn native_rt_window_edges() {
    let mut exp = test_experiment();
    exp.spectra.swap(0, 2);
    let config = ArrowSpectraExportConfig {
        min_rt: 100.0,
        ..long()
    };
    assert!(matches!(
        MSExperimentArrowExport::spectra_to_arrow(&exp, &config),
        Err(openms::Error::UnsortedData)
    ));
    let inverted = ArrowSpectraExportConfig {
        min_rt: 150.0,
        max_rt: 105.0,
        ..long()
    };
    assert_eq!(
        MSExperimentArrowExport::spectra_to_arrow(&test_experiment(), &inverted)
            .unwrap()
            .num_rows(),
        0
    );
}
