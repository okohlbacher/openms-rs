// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

#![cfg(feature = "arrow")]

//! The OpenSWATH Parquet writers (`.xic`, `.xim`) against the upstream class
//! tests and a retained C++-written file.
//!
//! Sections of `MSChromatogramParquetConsumer_test.cpp` and
//! `MobilogramParquetConsumer_test.cpp` at core `bc9cc12` are reproduced with
//! their literals (evidence tier 3). The chromatogram sections that build
//! their library with `TraMLFile` + `OpenSwathDataAccessHelper::
//! convertTargetedExp` (unported here) run where the library does not matter
//! (no chromatogram; destructor without output directory) with an empty one;
//! `MSChromatogramParquetConsumer_basic` runs with the library written out
//! from `MSChromatogramParquetConsumer_1_input.TraML` by hand (compound id,
//! sequence and charge; transition name, reference and product charge).
//!
//! `chrom_mzml_matches_the_retained_cpp_xic` is a differential check: the
//! 18 chromatograms of `MSChromatogramParquetConsumer_1_output.chrom.mzML` are
//! the ones the upstream fixture `XICParquetFile_1_input.xic` holds, written
//! there by the C++ consumer. Every row's ids, sequence, charges and flags,
//! and every decoded RT and intensity (C++ Numpress encoding decoded here
//! against Rust Numpress encoding decoded here), must be identical.

use openms::analysis::openswath::transition_experiment::{
    LightCompound, LightTargetedExperiment, LightTransition,
};
use openms::format::mobilogram_parquet_consumer::{
    MobilogramAnnotation, MobilogramParquetConsumer,
};
use openms::format::ms_chromatogram_parquet_consumer::MSChromatogramParquetConsumer;
use openms::format::xic_parquet_file::XicParquetFile;
use openms::format::xim_parquet_file::XimParquetFile;
use openms::{ChromatogramPeak, Error, MSChromatogram, MobilityPeak1D, Mobilogram};
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/openswath_parquet")
        .join(name)
}

fn scratch(test: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("openms-oswpqw-{}-{test}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn compound(id: &str, sequence: &str, charge: i32) -> LightCompound {
    LightCompound {
        id: id.into(),
        sequence: sequence.into(),
        charge,
        ..LightCompound::default()
    }
}

fn transition(name: &str, peptide_ref: &str) -> LightTransition {
    LightTransition {
        transition_name: name.into(),
        peptide_ref: peptide_ref.into(),
        ..LightTransition::default()
    }
}

/// The class tests' `pep1` / `tr1` library (`trailing_empty_chromatograms`).
fn pep1_library() -> LightTargetedExperiment {
    let mut library = LightTargetedExperiment::default();
    library.compounds.push(compound("pep1", "PEPTIDE", 2));
    let mut tr1 = transition("tr1", "pep1");
    tr1.precursor_mz = 600.2;
    tr1.product_mz = 500.2;
    tr1.fragment_charge = 1;
    tr1.fragment_nr = 7;
    tr1.set_fragment_type("y");
    tr1.set_detecting_transition(true);
    library.transitions.push(tr1);
    library
}

/// The IPF library of the `mixed_target_and_decoy_identifying_transitions`
/// sections of both class tests.
fn ipf_library() -> LightTargetedExperiment {
    let mut library = LightTargetedExperiment::default();
    library
        .compounds
        .push(compound("2407", "TSHVENDYIDNPSLALT(UniMod:21)TGPK", 3));
    let mut target = transition("tr_target", "2407");
    target.precursor_mz = 784.6967;
    target.product_mz = 326.1459;
    target.fragment_charge = 1;
    target.fragment_nr = 3;
    target.set_fragment_type("b");
    target.set_detecting_transition(true);
    target.set_identifying_transition(false);
    target.set_decoy(false);
    let mut decoy = target.clone();
    decoy.transition_name = "tr_ident_decoy".into();
    decoy.product_mz = 302.1710;
    decoy.set_detecting_transition(false);
    decoy.set_identifying_transition(true);
    decoy.set_decoy(true);
    library.transitions.push(target);
    library.transitions.push(decoy);
    library
}

fn chromatogram(native_id: &str, peaks: &[(f64, f32)]) -> MSChromatogram {
    MSChromatogram {
        native_id: native_id.into(),
        peaks: peaks
            .iter()
            .map(|&(rt, intensity)| ChromatogramPeak::new(rt, intensity))
            .collect(),
        ..MSChromatogram::default()
    }
}

// ------------------------------------------- MSChromatogramParquetConsumer ----

/// `MSChromatogramParquetConsumer_1_input.TraML` as the light library the
/// class test converts it into (ids, sequences, charges only).
#[cfg(feature = "mzml")]
fn traml_library() -> LightTargetedExperiment {
    let mut library = LightTargetedExperiment::default();
    library
        .compounds
        .push(compound("PEPTIDEA_Extra", "PEPTIDEAAAA", 2));
    for letter in ['A', 'B', 'C', 'D', 'E', 'F'] {
        let id = format!("PEPTIDE{letter}");
        library.compounds.push(compound(&id, &id, 2));
    }
    for (index, letter) in ['A', 'B', 'C', 'D', 'E', 'F'].iter().enumerate() {
        for k in 1..=3 {
            let mut t = transition(&(index * 3 + k).to_string(), &format!("PEPTIDE{letter}"));
            t.fragment_charge = 2;
            t.fragment_nr = -1;
            library.transitions.push(t);
        }
    }
    library
}

/// `START_SECTION(MSChromatogramParquetConsumer_basic)`.
#[cfg(feature = "mzml")]
#[test]
fn chromatogram_consumer_basic() {
    let experiment =
        openms::format::mzml::load(data("MSChromatogramParquetConsumer_1_output.chrom.mzML"))
            .unwrap();
    assert!(!experiment.chromatograms.is_empty());
    let out = scratch("xic-basic").join("out.xic");
    {
        let mut consumer =
            MSChromatogramParquetConsumer::new(&out, 1, "test_source", &traml_library());
        for chrom in &experiment.chromatograms {
            consumer.consume(&mut chrom.clone()).unwrap();
        }
    }
    assert!(out.exists());
    let chroms = XicParquetFile::new(&out).unwrap().load().unwrap();
    assert_eq!(chroms.len(), experiment.chromatograms.len());
    assert!(!chroms[0].rt.is_empty());
    assert!(!chroms[0].intensity.is_empty());
}

/// Differential: the same chromatograms, written here and by the C++
/// consumer (`XICParquetFile_1_input.xic`), agree row for row.
#[cfg(feature = "mzml")]
#[test]
fn chrom_mzml_matches_the_retained_cpp_xic() {
    let experiment =
        openms::format::mzml::load(data("MSChromatogramParquetConsumer_1_output.chrom.mzML"))
            .unwrap();
    let out = scratch("xic-diff").join("out.xic");
    let mut consumer = MSChromatogramParquetConsumer::new(&out, 5, "src", &traml_library());
    for chrom in &experiment.chromatograms {
        consumer.consume(&mut chrom.clone()).unwrap();
    }
    consumer.finalize().unwrap();
    let ours = XicParquetFile::new(&out).unwrap().load().unwrap();
    let theirs = XicParquetFile::new(data("XICParquetFile_1_input.xic"))
        .unwrap()
        .load()
        .unwrap();
    assert_eq!(ours.len(), theirs.len());
    for (a, b) in ours.iter().zip(&theirs) {
        assert_eq!(a.ms_level, b.ms_level);
        assert_eq!(a.precursor_id, b.precursor_id);
        assert_eq!(a.transition_id, b.transition_id);
        assert_eq!(a.modified_sequence, b.modified_sequence);
        assert_eq!(a.precursor_charge, b.precursor_charge);
        assert_eq!(a.product_charge, b.product_charge);
        assert_eq!(a.detecting_transition, b.detecting_transition);
        assert_eq!(a.precursor_decoy, b.precursor_decoy);
        assert_eq!(a.product_decoy, b.product_decoy);
        assert_eq!(a.transition_ordinal, b.transition_ordinal);
        assert_eq!(a.transition_type, b.transition_type);
        assert_eq!(a.annotation, b.annotation);
        assert_eq!(a.rt.len(), b.rt.len());
        for (x, y) in a.rt.iter().zip(&b.rt) {
            assert_eq!(
                x.to_bits(),
                y.to_bits(),
                "RT of transition {:?}",
                a.transition_id
            );
        }
        for (x, y) in a.intensity.iter().zip(&b.intensity) {
            assert_eq!(
                x.to_bits(),
                y.to_bits(),
                "intensity of transition {:?}",
                a.transition_id
            );
        }
    }
}

/// `START_SECTION(MSChromatogramParquetConsumer_empty_chromatograms)`.
#[test]
fn chromatogram_consumer_empty() {
    let out = scratch("xic-empty").join("out.xic");
    {
        let _consumer = MSChromatogramParquetConsumer::new(
            &out,
            1,
            "test_source",
            &LightTargetedExperiment::default(),
        );
    }
    let chroms = XicParquetFile::new(&out).unwrap().load().unwrap();
    assert_eq!(chroms.len(), 0);
}

/// `START_SECTION(MSChromatogramParquetConsumer_trailing_empty_chromatograms_are_written)`.
#[test]
fn chromatogram_consumer_trailing_empty() {
    let out = scratch("xic-trailing").join("out.xic");
    {
        let mut consumer =
            MSChromatogramParquetConsumer::new(&out, 1, "test_source", &pep1_library());
        consumer
            .consume(&mut chromatogram("tr1", &[(100.0, 1000.0)]))
            .unwrap();
        consumer.consume(&mut chromatogram("tr1", &[])).unwrap();
    }
    let chroms = XicParquetFile::new(&out).unwrap().load().unwrap();
    assert_eq!(chroms.len(), 2);
    assert!(!chroms[0].rt.is_empty());
    assert!(chroms[1].rt.is_empty());
    assert!(chroms[1].intensity.is_empty());
    // Native: a one-point chromatogram keeps its RT (the source stores it
    // with a Numpress fixed point of 0).
    assert_eq!(chroms[0].rt, [100.0]);
    assert_eq!(chroms[0].annotation, "y7^1");
}

/// `START_SECTION(MSChromatogramParquetConsumer_destructor_no_throw)`.
#[test]
fn chromatogram_consumer_destructor_no_throw() {
    let out = std::env::temp_dir()
        .join("openms_missing_dir")
        .join("xic_out.xic");
    let _ = std::fs::remove_dir_all(out.parent().unwrap());
    {
        let _consumer = MSChromatogramParquetConsumer::new(
            &out,
            1,
            "test_source",
            &LightTargetedExperiment::default(),
        );
    }
    assert!(!out.exists());
    // An explicit finalize reports what the destructor swallows.
    let mut consumer = MSChromatogramParquetConsumer::new(
        &out,
        1,
        "test_source",
        &LightTargetedExperiment::default(),
    );
    assert!(matches!(consumer.finalize(), Err(Error::Io(_))));
}

/// `START_SECTION(MSChromatogramParquetConsumer_mixed_target_and_decoy_identifying_transitions_keep_single_target_analyte)`.
#[test]
fn chromatogram_consumer_ipf_precursor_stays_target() {
    let out = scratch("xic-ipf").join("out.xic");
    {
        let mut consumer =
            MSChromatogramParquetConsumer::new(&out, 1, "test_source", &ipf_library());
        consumer
            .consume(&mut chromatogram("tr_target", &[(100.0, 1000.0)]))
            .unwrap();
        consumer
            .consume(&mut chromatogram("tr_ident_decoy", &[(100.0, 1000.0)]))
            .unwrap();
    }
    let xic = XicParquetFile::new(&out).unwrap();
    let chroms = xic.load().unwrap();
    assert_eq!(chroms.len(), 2);
    for chrom in &chroms {
        assert_eq!(chrom.precursor_decoy, Some(0));
    }
    let none: [&str; 0] = [];
    assert_eq!(xic.analytes(&none, true).unwrap().len(), 1);
}

/// A native id the library does not know is refused and nothing is written
/// for it (the source appends RUN_ID, SOURCE_FILE and MS_LEVEL before it
/// throws).
#[test]
fn chromatogram_consumer_unknown_native_id() {
    let out = scratch("xic-unknown").join("out.xic");
    let mut consumer = MSChromatogramParquetConsumer::new(&out, 1, "s", &pep1_library());
    let mut unknown = chromatogram("nope", &[(1.0, 2.0)]);
    assert!(matches!(
        consumer.consume(&mut unknown),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(unknown.peaks.len(), 1);
    let mut precursor = chromatogram("missing_Precursor_i0", &[(1.0, 2.0)]);
    assert!(matches!(
        consumer.consume(&mut precursor),
        Err(Error::InvalidValue(_))
    ));
    let mut known = chromatogram("pep1_Precursor_i0", &[(1.0, 2.0), (2.0, 3.0), (3.0, 4.0)]);
    consumer.consume(&mut known).unwrap();
    assert!(known.peaks.is_empty());
    consumer.finalize().unwrap();
    let chroms = XicParquetFile::new(&out).unwrap().load().unwrap();
    assert_eq!(chroms.len(), 1);
    assert_eq!(chroms[0].ms_level, 1);
    assert_eq!(chroms[0].annotation, "Precursor_i0");
    assert_eq!(chroms[0].transition_id, None);
    assert!(consumer.consume(&mut chromatogram("tr1", &[])).is_err());
}

// ----------------------------------------------- MobilogramParquetConsumer ----

fn mobilogram(peaks: &[(f64, f32)]) -> Mobilogram {
    Mobilogram::from_peaks(
        peaks
            .iter()
            .map(|&(mobility, intensity)| MobilityPeak1D::new(mobility, intensity))
            .collect(),
    )
}

/// `START_SECTION(MobilogramParquetConsumer_basic)`.
#[test]
fn mobilogram_consumer_basic() {
    let out = scratch("xim-basic").join("out.xim");
    let mut m = mobilogram(&[(1.0, 10.0), (2.0, 20.0)]);
    m.rt = 123.45;
    {
        let mut consumer = MobilogramParquetConsumer::new(
            &out,
            1,
            "test_source",
            &LightTargetedExperiment::default(),
        );
        consumer
            .consume_mobilogram(&m, &MobilogramAnnotation::default())
            .unwrap();
    }
    assert!(out.exists());
}

/// `START_SECTION(MobilogramParquetConsumer_const_mobilogram_not_modified)`.
#[test]
fn mobilogram_consumer_input_unchanged() {
    let out = scratch("xim-const").join("out.xim");
    let mut m = mobilogram(&[(0.5, 100.0), (1.5, 200.0)]);
    m.rt = 1.0;
    {
        let mut consumer =
            MobilogramParquetConsumer::new(&out, 1, "src", &LightTargetedExperiment::default());
        consumer
            .consume_mobilogram(&m, &MobilogramAnnotation::default())
            .unwrap();
    }
    assert_eq!(m.peaks.len(), 2);
    assert_eq!(m.peaks[0].mobility, 0.5);
    assert_eq!(m.peaks[1].mobility, 1.5);
}

/// `START_SECTION(MobilogramParquetConsumer_round_trip_file_written)`.
#[test]
fn mobilogram_consumer_round_trip_file_written() {
    let out = scratch("xim-rt").join("out.xim");
    let mut m = mobilogram(&[(0.80, 1000.0), (0.90, 2000.0), (1.00, 3000.0)]);
    m.rt = 50.0;
    {
        let mut consumer = MobilogramParquetConsumer::new(
            &out,
            42,
            "run_src",
            &LightTargetedExperiment::default(),
        );
        for native in ["tr1", "tr2"] {
            let annotation = MobilogramAnnotation {
                mobilogram_type: "ms2".into(),
                ms_level: Some(2),
                transition_id: None,
                transition_native_id: native.into(),
                feature_rt: Some(50.0),
                feature_id: Some(7),
            };
            consumer.consume_mobilogram(&m, &annotation).unwrap();
        }
    }
    assert!(out.exists());
    assert!(std::fs::metadata(&out).unwrap().len() > 0);
}

/// `START_SECTION(MobilogramParquetConsumer_empty_mobilograms)`.
#[test]
fn mobilogram_consumer_empty() {
    let out = scratch("xim-empty").join("out.xim");
    {
        let _consumer = MobilogramParquetConsumer::new(
            &out,
            1,
            "test_source",
            &LightTargetedExperiment::default(),
        );
    }
    assert!(out.exists());
    assert!(
        XimParquetFile::new(&out)
            .unwrap()
            .load()
            .unwrap()
            .is_empty()
    );
}

/// `START_SECTION(MobilogramParquetConsumer_destructor_flushes)` and
/// `(..._explicit_finalize_writes_file)`.
#[test]
fn mobilogram_consumer_destructor_and_finalize_write() {
    let m = mobilogram(&[(1.0, 100.0)]);
    for explicit in [false, true] {
        let out = scratch(&format!("xim-final-{explicit}")).join("out.xim");
        {
            let mut consumer = MobilogramParquetConsumer::new(
                &out,
                1,
                "test_source",
                &LightTargetedExperiment::default(),
            );
            consumer
                .consume_mobilogram(&m, &MobilogramAnnotation::default())
                .unwrap();
            if explicit {
                consumer.finalize().unwrap();
                assert!(
                    consumer
                        .consume_mobilogram(&m, &MobilogramAnnotation::default())
                        .is_err()
                );
            }
        }
        assert!(out.exists());
        assert!(std::fs::metadata(&out).unwrap().len() > 0);
    }
}

/// `START_SECTION(MobilogramParquetConsumer_mixed_target_and_decoy_identifying_transitions_keep_single_target_analyte)`.
#[test]
fn mobilogram_consumer_ipf_precursor_stays_target() {
    let out = scratch("xim-ipf").join("out.xim");
    let mut m = mobilogram(&[(1.0, 1000.0)]);
    m.rt = 100.0;
    {
        let mut consumer = MobilogramParquetConsumer::new(&out, 1, "test_source", &ipf_library());
        for native in ["tr_target", "tr_ident_decoy"] {
            let annotation = MobilogramAnnotation {
                mobilogram_type: "transition".into(),
                ms_level: Some(2),
                transition_id: None,
                transition_native_id: native.into(),
                feature_rt: Some(100.0),
                feature_id: Some(7),
            };
            consumer.consume_mobilogram(&m, &annotation).unwrap();
        }
    }
    let xim = XimParquetFile::new(&out).unwrap();
    let mobilograms = xim.load().unwrap();
    assert_eq!(mobilograms.len(), 2);
    for m in &mobilograms {
        assert_eq!(m.precursor_decoy, Some(0));
    }
    let none: [&str; 0] = [];
    assert_eq!(xim.analytes(&none, true).unwrap().len(), 1);
}

// ----------------------------------------------------- XIPMParquetConsumer ----

use openms::analysis::openswath::openswath_helper::compute_precursor_id;
use openms::analysis::openswath::peak_map_extractor::ExtractedPeakMap;
use openms::format::xipm_parquet_consumer::XIPMParquetConsumer;
use openms::format::xipm_parquet_file::{XipmParquetFile, XipmQuery};

/// The class tests' `makeExperiment_`.
fn xipm_library() -> LightTargetedExperiment {
    let mut library = pep1_library();
    library.compounds[0].rt = 100.0;
    library.compounds[0].drift_time = 1.1;
    library
}

/// `makeMixedIPFExperiment_`.
fn xipm_ipf_library() -> LightTargetedExperiment {
    let mut library = xipm_library();
    let mut decoy = transition("tr1_ident_decoy", "pep1");
    decoy.precursor_mz = 600.2;
    decoy.product_mz = 501.2;
    decoy.fragment_charge = 1;
    decoy.fragment_nr = 3;
    decoy.set_fragment_type("b");
    decoy.set_detecting_transition(false);
    decoy.set_identifying_transition(true);
    decoy.set_decoy(true);
    library.transitions.push(decoy);
    library
}

/// `makeTransitionPeakMap_(native_id, target_mz, mz_values)`.
fn transition_peak_map(native_id: &str, target_mz: f64, mz: &[f64]) -> ExtractedPeakMap {
    ExtractedPeakMap {
        native_id: native_id.into(),
        target_mz,
        target_rt: 100.0,
        target_ion_mobility: 1.1,
        rt_start: 95.0,
        rt_end: 105.0,
        mz: mz.to_vec(),
        rt: vec![100.0; mz.len()],
        ion_mobility: vec![1.05; mz.len()],
        intensity: vec![1000.0; mz.len()],
    }
}

/// `makeTransitionPeakMap_()`.
fn tr1_peak_map() -> ExtractedPeakMap {
    ExtractedPeakMap {
        native_id: "tr1".into(),
        target_mz: 500.2,
        target_rt: 100.0,
        target_ion_mobility: 1.1,
        rt_start: 95.0,
        rt_end: 105.0,
        mz: vec![500.19, 500.20],
        rt: vec![100.0, 101.0],
        ion_mobility: vec![1.05, 1.08],
        intensity: vec![1000.0, 900.0],
    }
}

/// `makePrecursorPeakMap_()`.
fn precursor_peak_map() -> ExtractedPeakMap {
    ExtractedPeakMap {
        native_id: compute_precursor_id("pep1", 0),
        target_mz: 600.2,
        target_rt: 100.0,
        target_ion_mobility: 1.1,
        rt_start: 95.0,
        rt_end: 105.0,
        mz: vec![600.19],
        rt: vec![100.0],
        ion_mobility: vec![1.06],
        intensity: vec![500.0],
    }
}

/// `writeTestFile_(run_id, source_file)` of `XIPMParquetFile_test.cpp`.
fn write_xipm(test: &str) -> PathBuf {
    let out = scratch(test).join("out.xipm");
    let mut consumer = XIPMParquetConsumer::new(&out, &xipm_library());
    consumer
        .consume_peak_map(&tr1_peak_map(), 7, "run1.mzML", 2)
        .unwrap();
    consumer
        .consume_peak_map(&precursor_peak_map(), 7, "run1.mzML", 1)
        .unwrap();
    consumer.finalize().unwrap();
    out
}

/// `START_SECTION(XIPMParquetConsumer_basic_roundtrip)`.
#[test]
fn xipm_consumer_basic_roundtrip() {
    let out = scratch("xipm-basic").join("out.xipm");
    {
        let mut consumer = XIPMParquetConsumer::new(&out, &xipm_library());
        consumer
            .consume_peak_map(&tr1_peak_map(), 7, "run1.mzML", 2)
            .unwrap();
        consumer
            .consume_peak_map(&precursor_peak_map(), 7, "run1.mzML", 1)
            .unwrap();
    }
    assert!(out.exists());
    assert!(std::fs::metadata(&out).unwrap().len() > 0);
    let peak_maps = XipmParquetFile::new(&out).unwrap().load().unwrap();
    assert_eq!(peak_maps.len(), 2);
    assert_eq!(peak_maps[0].target_rt, Some(100.0));
    assert_eq!(peak_maps[1].target_rt, Some(100.0));
}

/// `START_SECTION(XIPMParquetConsumer_empty_file)`.
#[test]
fn xipm_consumer_empty_file() {
    let out = scratch("xipm-empty").join("out.xipm");
    {
        let _consumer = XIPMParquetConsumer::new(&out, &xipm_library());
    }
    assert!(out.exists());
}

/// `START_SECTION(XIPMParquetConsumer_explicit_finalize_writes_file)`.
#[test]
fn xipm_consumer_explicit_finalize() {
    let out = scratch("xipm-final").join("out.xipm");
    {
        let mut consumer = XIPMParquetConsumer::new(&out, &xipm_library());
        consumer
            .consume_peak_map(&tr1_peak_map(), 7, "run1.mzML", 2)
            .unwrap();
        consumer.finalize().unwrap();
    }
    assert!(out.exists());
    assert!(std::fs::metadata(&out).unwrap().len() > 0);
}

/// `START_SECTION(XIPMParquetConsumer_mixed_target_and_decoy_identifying_transitions_keep_target_precursor_decoy)`.
#[test]
fn xipm_consumer_ipf_precursor_stays_target() {
    let out = scratch("xipm-ipf").join("out.xipm");
    {
        let mut consumer = XIPMParquetConsumer::new(&out, &xipm_ipf_library());
        consumer
            .consume_peak_map(
                &transition_peak_map("tr1", 500.2, &[500.19, 500.20]),
                7,
                "run1.mzML",
                2,
            )
            .unwrap();
        consumer
            .consume_peak_map(
                &transition_peak_map("tr1_ident_decoy", 501.2, &[501.19, 501.20]),
                7,
                "run1.mzML",
                2,
            )
            .unwrap();
    }
    let peak_maps = XipmParquetFile::new(&out).unwrap().load().unwrap();
    assert_eq!(peak_maps.len(), 2);
    let mut product_decoy_zero = 0;
    let mut product_decoy_one = 0;
    for peak_map in &peak_maps {
        assert_eq!(peak_map.precursor_decoy, Some(0));
        match peak_map.product_decoy {
            Some(0) => product_decoy_zero += 1,
            Some(1) => product_decoy_one += 1,
            other => panic!("product decoy {other:?}"),
        }
    }
    assert_eq!(product_decoy_zero, 1);
    assert_eq!(product_decoy_one, 1);
}

/// `XIPMParquetFile_test.cpp`, every section, on the file its
/// `writeTestFile_` / `writeEmptyTestFile_` write.
#[test]
fn xipm_file_class_test_on_written_file() {
    let file = write_xipm("xipm-file");
    let xipm = XipmParquetFile::new(&file).unwrap();
    // load
    let peak_maps = xipm.load().unwrap();
    assert_eq!(peak_maps.len(), 2);
    assert_eq!(peak_maps[0].mz.len(), peak_maps[0].rt.len());
    assert_eq!(peak_maps[0].mz.len(), peak_maps[0].ion_mobility.len());
    assert_eq!(peak_maps[0].mz.len(), peak_maps[0].intensity.len());
    assert_eq!(peak_maps[0].target_rt, Some(100.0));
    // getPeakMaps(...filters...)
    let precursor = xipm
        .peak_maps(&XipmQuery {
            ms_level: Some(1),
            run_id: Some(7),
            peakmap_type: "precursor".into(),
            ..XipmQuery::default()
        })
        .unwrap();
    assert_eq!(precursor.len(), 1);
    assert_eq!(precursor[0].peakmap_type, "precursor");
    let transition = xipm
        .peak_maps(&XipmQuery {
            transition_id: Some(1),
            ms_level: Some(2),
            run_id: Some(7),
            peakmap_type: "transition".into(),
            ..XipmQuery::default()
        })
        .unwrap();
    assert_eq!(transition.len(), 1);
    assert_eq!(transition[0].annotation, "y7^1");
    assert_eq!(transition[0].target_rt, Some(100.0));
    // getPeakMaps_multi_file
    let multi = XipmParquetFile::from_files(&[&file, &file]).unwrap();
    assert_eq!(multi.load().unwrap().len(), 4);
    // getRuns
    let runs = xipm.runs().unwrap();
    assert_eq!(runs.len(), 1);
    assert_eq!(runs[0].run_id, 7);
    // getColumns
    let columns = xipm.columns().unwrap();
    for wanted in ["MZ_DATA", "RT_DATA", "MOBILITY_DATA", "TARGET_RT"] {
        assert!(columns.iter().any(|c| c == wanted), "{wanted}");
    }
    // getPeakMaps_empty_file, getRuns_empty_file
    let empty = scratch("xipm-file-empty").join("out.xipm");
    XIPMParquetConsumer::new(&empty, &xipm_library())
        .finalize()
        .unwrap();
    let empty = XipmParquetFile::new(&empty).unwrap();
    assert!(empty.load().unwrap().is_empty());
    assert!(empty.runs().unwrap().is_empty());
}

/// The file written here and the retained upstream `XIPMParquetFile_1_input.xipm`
/// (the same two peak maps, written by the C++ consumer before TARGET_RT
/// existed) decode to the same rows.
#[test]
fn xipm_written_matches_the_retained_cpp_xipm() {
    let ours = XipmParquetFile::new(write_xipm("xipm-diff"))
        .unwrap()
        .load()
        .unwrap();
    let theirs = XipmParquetFile::new(data("XIPMParquetFile_1_input.xipm"))
        .unwrap()
        .load()
        .unwrap();
    assert_eq!(ours.len(), theirs.len());
    for (a, b) in ours.iter().zip(&theirs) {
        assert_eq!(a.run_id, b.run_id);
        assert_eq!(a.source_file, b.source_file);
        assert_eq!(a.ms_level, b.ms_level);
        assert_eq!(a.peakmap_type, b.peakmap_type);
        assert_eq!(a.precursor_id, b.precursor_id);
        assert_eq!(a.transition_id, b.transition_id);
        assert_eq!(a.modified_sequence, b.modified_sequence);
        assert_eq!(a.precursor_charge, b.precursor_charge);
        assert_eq!(a.product_charge, b.product_charge);
        assert_eq!(a.detecting_transition, b.detecting_transition);
        assert_eq!(a.precursor_decoy, b.precursor_decoy);
        assert_eq!(a.product_decoy, b.product_decoy);
        assert_eq!(a.transition_ordinal, b.transition_ordinal);
        assert_eq!(a.transition_type, b.transition_type);
        assert_eq!(a.annotation, b.annotation);
        assert_eq!(a.target_mz.to_bits(), b.target_mz.to_bits());
        assert_eq!(a.target_rt, b.target_rt);
        assert_eq!(a.target_ion_mobility, b.target_ion_mobility);
        assert_eq!(a.rt_start, b.rt_start);
        assert_eq!(a.rt_end, b.rt_end);
        for (x, y) in [
            (&a.mz, &b.mz),
            (&a.rt, &b.rt),
            (&a.ion_mobility, &b.ion_mobility),
            (&a.intensity, &b.intensity),
        ] {
            assert_eq!(x.len(), y.len());
            for (p, q) in x.iter().zip(y.iter()) {
                assert_eq!(p.to_bits(), q.to_bits());
            }
        }
    }
}
