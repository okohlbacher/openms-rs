// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! ANALYSIS/OPENSWATH extraction: `ChromatogramExtractorAlgorithm`,
//! `ChromatogramExtractor`, `PeakMapExtractor` and `SwathWindowLoader`.
//!
//! Expected values are transcribed from the upstream class tests at the core
//! pin `bc9cc12` (evidence tier 3), cited per test as `<file>:<section>`, all
//! under `src/tests/class_tests/openms/source/`. Tests marked "derived" check
//! an independently derived invariant or a native guard (tier 4). Tests marked
//! "oracle" replay values printed by the C++ Release build on ibminode06
//! (tier 1; see `tests/data/openswath_extraction_provenance.json`).

use openms::Error;
#[cfg(feature = "mzml")]
use openms::analysis::openswath::chromatogram_extractor::ChromatogramExtractor;
use openms::analysis::openswath::chromatogram_extractor::{
    prepare_coordinates, return_chromatogram,
};
use openms::analysis::openswath::chromatogram_extractor_algorithm::{
    ChromatogramExtractorAlgorithm, ExtractionCoordinates, ExtractionFilter,
};
use openms::analysis::openswath::data_structures::Chromatogram;
use openms::analysis::openswath::openswath_helper::compute_transition_group_id;
use openms::analysis::openswath::peak_map_extractor::{ExtractedPeakMap, PeakMapExtractor};
use openms::analysis::openswath::spectrum_access_openms::SpectrumAccessOpenMS;
use openms::analysis::openswath::swath_map::SwathMap;
use openms::analysis::openswath::swath_map_mass_correction::{
    MRMTransitionGroupType, SwathMapMassCorrection,
};
#[cfg(feature = "mzml")]
use openms::analysis::openswath::swath_qc::ChargeDistribution;
use openms::analysis::openswath::swath_qc::SwathQC;
use openms::analysis::openswath::swath_window_loader::{
    annotate_swath_maps, annotate_swath_maps_from_file, parse_swath_windows, read_swath_windows,
};
use openms::analysis::openswath::transition_experiment::{
    LightCompound, LightTargetedExperiment, LightTransition,
};
use openms::constants::C13C12_MASSDIFF_U;
use openms::kernel::{DataArray, MSExperiment, MSSpectrum, Peak1D};
use openms::metadata::{ChromatogramType, ExperimentalSettings, MetaValue};
use openms::param::ParamValue;
use std::path::PathBuf;
use std::sync::Arc;

fn real_similar(a: f64, b: f64) {
    let absdiff = (a - b).abs();
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{a} !~ {b} (absdiff {absdiff}, ratio {ratio})"
    );
}

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn coord(mz: f64, rt_start: f64, rt_end: f64, id: &str, im: f64) -> ExtractionCoordinates {
    ExtractionCoordinates {
        mz,
        ion_mobility: im,
        rt_start,
        rt_end,
        id: id.into(),
        ..ExtractionCoordinates::default()
    }
}

/// `find_max_helper` of ChromatogramExtractorAlgorithm_test.cpp: the first
/// maximum intensity and its RT.
fn find_max(chrom: &Chromatogram) -> (f64, f64) {
    let mut max_value = -1.0;
    let mut found_at = -1.0;
    for (&rt, &intensity) in chrom
        .time_array()
        .data
        .iter()
        .zip(&chrom.intensity_array().data)
    {
        if intensity > max_value {
            max_value = intensity;
            found_at = rt;
        }
    }
    (max_value, found_at)
}

fn access(exp: MSExperiment) -> SpectrumAccessOpenMS {
    SpectrumAccessOpenMS::new(Arc::new(exp))
}

fn spectrum(rt: f64, peaks: &[(f64, f32)], im: Option<&[f32]>) -> MSSpectrum {
    let mut s = MSSpectrum {
        rt,
        ms_level: 2,
        peaks: peaks.iter().map(|&(mz, i)| Peak1D::new(mz, i)).collect(),
        ..MSSpectrum::default()
    };
    if let Some(im) = im {
        s.float_data_arrays
            .push(DataArray::new("Ion Mobility", im.to_vec()));
    }
    s
}

// ---------------------------------------------------------------------------
// ChromatogramExtractorAlgorithm
// ---------------------------------------------------------------------------

/// The 22-point arrays shared by the `extract_value_tophat` sections.
const MZ_ARR: [f64; 22] = [
    400.0, 400.01, 400.02, 400.03, 400.04, 400.05, 400.06, 400.07, 400.08, 400.09, 400.1, 400.11,
    400.12, 400.13, 400.14, 400.15, 400.16, 400.17, 400.18, 400.19, 450.0, 500.0,
];
const INT_ARR: [f64; 22] = [
    8.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 900.0, 800.0, 700.0, 600.0,
    500.0, 400.0, 300.0, 200.0, 100.0, 0.0, 10.0, 10.0,
];
const IM_ARR: [f64; 22] = [
    100.0, 200.01, 100.02, 200.03, 100.04, 200.05, 100.06, 200.07, 100.08, 200.09, 100.1, 200.11,
    100.12, 200.13, 100.14, 200.15, 100.16, 200.17, 100.18, 200.19, 300.1, 300.2,
];

/// ChromatogramExtractorAlgorithm_test.cpp: `extract_value_tophat`, Th and
/// ppm windows, including the documented non-monotonic misuse.
#[test]
fn tophat_class_test_sequence() {
    let x = ChromatogramExtractorAlgorithm::new();
    let mut pos = 0usize;
    let f = |mz: f64, w: f64, ppm: bool, pos: &mut usize| {
        x.extract_value_tophat(&MZ_ARR, &INT_ARR, pos, mz, w, ppm)
            .unwrap()
    };
    real_similar(f(399.805, 0.2, false, &mut pos), 0.0);
    real_similar(f(399.91, 0.2, false, &mut pos), 108.0);
    real_similar(f(400.0, 0.2, false, &mut pos), 4508.0);
    real_similar(f(400.05, 0.2, false, &mut pos), 8400.0);
    let v = f(400.1, 0.2, false, &mut pos);
    real_similar(v, 9000.0);
    assert_eq!(v as i64, 9000);
    real_similar(f(400.28, 0.2, false, &mut pos), 100.0);
    real_similar(f(500.0, 0.2, false, &mut pos), 10.0);
    // Not monotonically increasing: the class test only asserts "not 9000".
    assert_ne!(f(400.1, 0.2, false, &mut pos) as i64, 9000);

    let mut pos = 0usize;
    real_similar(f(399.89, 500.0, true, &mut pos), 0.0);
    real_similar(f(399.91, 500.0, true, &mut pos), 8.0);
    real_similar(f(399.92, 500.0, true, &mut pos), 108.0);
    real_similar(f(400.0, 500.0, true, &mut pos), 4508.0);
    real_similar(f(400.05, 500.0, true, &mut pos), 8400.0);
    real_similar(f(400.1, 500.0, true, &mut pos), 9000.0);
}

/// ChromatogramExtractorAlgorithm_test.cpp: `[EXTRA IM] extract_value_tophat`.
#[test]
fn tophat_im_class_test_sequence() {
    let x = ChromatogramExtractorAlgorithm::new();
    let mut pos = 0usize;
    let f = |mz: f64, im: f64, pos: &mut usize| {
        x.extract_value_tophat_im(&MZ_ARR, &INT_ARR, &IM_ARR, pos, mz, im, 0.2, 0.3, false)
            .unwrap()
    };
    real_similar(f(399.805, 100.0, &mut pos), 0.0);
    real_similar(f(399.91, 100.0, &mut pos), 8.0);
    real_similar(f(400.0, 100.0, &mut pos), 2008.0);
    real_similar(f(400.05, 100.0, &mut pos), 4100.0);
    let v = f(400.1, 100.0, &mut pos);
    real_similar(v, 4100.0);
    assert_eq!(v as i64, 4100);
    real_similar(f(400.28, 100.0, &mut pos), 0.0);
    real_similar(f(400.28, 200.0, &mut pos), 0.0);
    real_similar(f(400.28, 200.1, &mut pos), 0.0);
    real_similar(f(500.0, 300.0, &mut pos), 0.0);
    real_similar(f(500.0, 300.1, &mut pos), 10.0);
    assert_ne!(f(400.1, 100.0, &mut pos) as i64, 9000);
}

/// ChromatogramExtractorAlgorithm_test.cpp: the two `[EXTRA]` sections with a
/// 0.001 Th window around 400.0001.
#[test]
fn tophat_narrow_window_first_point() {
    let x = ChromatogramExtractorAlgorithm::new();
    let mut pos = 0usize;
    real_similar(
        x.extract_value_tophat(&MZ_ARR, &INT_ARR, &mut pos, 399.805, 0.2, false)
            .unwrap(),
        0.0,
    );
    real_similar(
        x.extract_value_tophat(&MZ_ARR, &INT_ARR, &mut pos, 400.0001, 0.001, false)
            .unwrap(),
        8.0,
    );
    let mut pos = 0usize;
    real_similar(
        x.extract_value_tophat_im(
            &MZ_ARR, &INT_ARR, &IM_ARR, &mut pos, 399.805, 100.0, 0.2, 0.3, false,
        )
        .unwrap(),
        0.0,
    );
    real_similar(
        x.extract_value_tophat_im(
            &MZ_ARR, &INT_ARR, &IM_ARR, &mut pos, 400.0001, 100.0, 0.001, 0.3, false,
        )
        .unwrap(),
        8.0,
    );
}

/// Oracle: the source walk's two quirks, measured on the C++ Release build.
/// A target above every peak counts the last peak twice; a leftward walk of
/// two or more steps never adds the first peak.
#[test]
fn tophat_source_walk_quirks() {
    let x = ChromatogramExtractorAlgorithm::new();
    // Last peak counted twice: 100 + 7 + 7.
    let mz = [100.0, 200.0, 200.02];
    let int = [1.0, 100.0, 7.0];
    let mut pos = 0usize;
    let v = x
        .extract_value_tophat(&mz, &int, &mut pos, 200.05, 0.2, false)
        .unwrap();
    assert_eq!(v, 114.0);
    assert_eq!(pos, 3);
    // First peak skipped: window (99.8, 100.2) holds all three peaks.
    let mz = [100.0, 100.05, 100.1];
    let int = [1.0, 10.0, 100.0];
    let mut pos = 0usize;
    let v = x
        .extract_value_tophat(&mz, &int, &mut pos, 100.1, 0.4, false)
        .unwrap();
    assert_eq!(v, 110.0);
}

/// Derived: the walk refuses a short intensity array and a position past the
/// end instead of reading out of bounds, and an empty spectrum yields zero.
#[test]
fn tophat_guards() {
    let x = ChromatogramExtractorAlgorithm::new();
    let mut pos = 0usize;
    assert!(matches!(
        x.extract_value_tophat(&[1.0, 2.0], &[1.0], &mut pos, 1.0, 1.0, false),
        Err(Error::InvalidValue(_))
    ));
    let mut pos = 3usize;
    assert!(matches!(
        x.extract_value_tophat(&[1.0, 2.0], &[1.0, 2.0], &mut pos, 1.0, 1.0, false),
        Err(Error::InvalidValue(_))
    ));
    let mut pos = 0usize;
    assert!(matches!(
        x.extract_value_tophat_im(
            &[1.0, 2.0],
            &[1.0, 2.0],
            &[1.0],
            &mut pos,
            1.0,
            1.0,
            1.0,
            1.0,
            false
        ),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(
        x.extract_value_tophat(&[], &[], &mut pos, 1.0, 1.0, false)
            .unwrap(),
        0.0
    );
}

/// ChromatogramExtractorAlgorithm_test.cpp: `extractChromatograms` over
/// ChromatogramExtractor_input.mzML, and the missing-ion-mobility exception.
#[cfg(feature = "mzml")]
#[test]
fn extract_chromatograms_class_test_file() {
    let bytes = std::fs::read(data(
        "openswath_extraction_ChromatogramExtractor_input.mzML",
    ))
    .unwrap();
    let exp = openms::format::mzml::read(bytes.as_slice()).unwrap();
    let input = access(exp);
    let coordinates = vec![
        coord(618.31, 0.0, -1.0, "tr1", 0.0),
        coord(628.45, 0.0, -1.0, "tr2", 0.0),
        coord(654.38, 0.0, -1.0, "tr3", 0.0),
    ];
    let mut out = vec![Chromatogram::new(); 3];
    ChromatogramExtractorAlgorithm::new()
        .extract_chromatograms(&input, &mut out, &coordinates, 0.05, false, -1.0, "tophat")
        .unwrap();
    assert_eq!(out[0].time_array().data.len(), 59);
    assert_eq!(out[0].intensity_array().data.len(), 59);
    let (max_value, found_at) = find_max(&out[1]);
    real_similar(max_value, 169.792);
    real_similar(found_at, 3120.26);
    let (max_value, found_at) = find_max(&out[2]);
    real_similar(max_value, 577.33);
    real_similar(found_at, 3120.26);
    let (max_value, found_at) = find_max(&out[0]);
    real_similar(max_value, 35.593);
    real_similar(found_at, 3055.16);

    // There is no ion mobility, so this should not work; output unchanged.
    let before = out.clone();
    assert!(matches!(
        ChromatogramExtractorAlgorithm::new().extract_chromatograms(
            &input,
            &mut out,
            &coordinates,
            0.05,
            false,
            1.0,
            "tophat"
        ),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(out, before);
}

/// The four-spectrum ion-mobility experiment of the `[EXTRA]
/// extractChromatograms` section.
fn im_experiment() -> SpectrumAccessOpenMS {
    let mut exp = MSExperiment::new();
    for i in 0..4 {
        let mut peaks = Vec::new();
        let mut im = Vec::new();
        for base in [618.3, 628.4] {
            for k in 0..10 {
                peaks.push((base + f64::from(k) * 0.01, (100 * i + k * 2) as f32));
                im.push((100 + k * 10) as f32);
            }
        }
        exp.spectra.push(spectrum(f64::from(i), &peaks, Some(&im)));
    }
    access(exp)
}

/// ChromatogramExtractorAlgorithm_test.cpp: `[EXTRA] extractChromatograms`,
/// without, with a small and with a larger ion-mobility window.
#[test]
fn extract_chromatograms_ion_mobility() {
    let input = im_experiment();
    let coordinates = vec![
        coord(618.31, 0.0, -1.0, "tr1", 120.0),
        coord(628.45, 0.0, -1.0, "tr2", 170.0),
    ];
    for (im_window, first, second) in [
        (-1.0, 1830.0, 2790.0),
        (15.0, 304.0, 314.0),
        (30.0, 303.0 + 304.0 + 305.0, 313.0 + 314.0 + 315.0),
    ] {
        let mut out = vec![Chromatogram::new(); 2];
        ChromatogramExtractorAlgorithm::new()
            .extract_chromatograms(
                &input,
                &mut out,
                &coordinates,
                0.10,
                false,
                im_window,
                "tophat",
            )
            .unwrap();
        assert_eq!(out[0].time_array().data.len(), 4);
        assert_eq!(out[0].intensity_array().data.len(), 4);
        let (max_value, found_at) = find_max(&out[0]);
        real_similar(max_value, first);
        real_similar(found_at, 3.0);
        let (max_value, found_at) = find_max(&out[1]);
        real_similar(max_value, second);
        real_similar(found_at, 3.0);
    }
}

/// Derived: argument checks, their source order, RT windows, appending and
/// atomicity.
#[test]
fn extract_chromatograms_contract() {
    let x = ChromatogramExtractorAlgorithm::new();
    let input = im_experiment();
    let sorted = vec![
        coord(618.31, 0.0, -1.0, "a", -1.0),
        coord(628.45, 1.5, 2.5, "b", -1.0),
    ];
    // Size mismatch.
    let mut one = vec![Chromatogram::new()];
    assert!(matches!(
        x.extract_chromatograms(&input, &mut one, &sorted, 0.1, false, -1.0, "tophat"),
        Err(Error::InvalidValue(_))
    ));
    // Unknown filter.
    let mut out = vec![Chromatogram::new(); 2];
    assert!(matches!(
        x.extract_chromatograms(&input, &mut out, &sorted, 0.1, false, -1.0, "gauss"),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(
        ExtractionFilter::from_name("bartlett").unwrap(),
        ExtractionFilter::Bartlett
    );
    // Bartlett is declared, not implemented; the output stays unchanged.
    assert!(matches!(
        x.extract_chromatograms(&input, &mut out, &sorted, 0.1, false, -1.0, "bartlett"),
        Err(Error::Unsupported(_))
    ));
    assert_eq!(out, vec![Chromatogram::new(); 2]);
    // Unsorted coordinates.
    let unsorted = vec![sorted[1].clone(), sorted[0].clone()];
    assert!(matches!(
        x.extract_chromatograms(&input, &mut out, &unsorted, 0.1, false, -1.0, "tophat"),
        Err(Error::InvalidValue(_))
    ));
    // RT window [1.5, 2.5] keeps only the spectrum at RT 2; results append.
    x.extract_chromatograms(&input, &mut out, &sorted, 0.1, false, -1.0, "tophat")
        .unwrap();
    x.extract_chromatograms(&input, &mut out, &sorted, 0.1, false, -1.0, "tophat")
        .unwrap();
    assert_eq!(out[0].time_array().data.len(), 8);
    assert_eq!(out[1].time_array().data, vec![2.0, 2.0]);
    // An empty input returns before any check, as in the source.
    let empty = access(MSExperiment::new());
    let mut none: Vec<Chromatogram> = Vec::new();
    x.extract_chromatograms(&empty, &mut none, &sorted, 0.1, false, 1.0, "nonsense")
        .unwrap();
}

/// Derived: the sort predicates and the m/z ordering helper.
#[test]
fn extraction_coordinates_ordering() {
    let a = coord(100.0, 0.0, 0.0, "a", 0.0);
    let b = coord(200.0, 0.0, 0.0, "b", 0.0);
    assert!(ExtractionCoordinates::sort_by_mz(&a, &b));
    assert!(!ExtractionCoordinates::sort_by_mz(&b, &a));
    assert!(ExtractionCoordinates::sort_reverse_by_mz(&b, &a));
    assert_eq!(
        ExtractionCoordinates::cmp_by_mz(&a, &b),
        std::cmp::Ordering::Less
    );
    let d = ExtractionCoordinates::default();
    assert_eq!(
        (d.mz, d.ion_mobility, d.rt_start, d.rt_end),
        (0.0, 0.0, 0.0, 0.0)
    );
}

// ---------------------------------------------------------------------------
// PeakMapExtractor
// ---------------------------------------------------------------------------

const PM_MZ: [f64; 5] = [499.99, 500.00, 500.02, 500.03, 500.10];
const PM_INT: [f32; 5] = [10.0, 11.0, 12.0, 13.0, 14.0];
const PM_IM: [f32; 5] = [1.00, 1.05, 1.10, 1.15, 1.30];

fn pm_spectrum(rt: f64, offset: f32) -> MSSpectrum {
    let peaks: Vec<(f64, f32)> = PM_MZ
        .iter()
        .zip(PM_INT)
        .map(|(&mz, i)| (mz, i + offset))
        .collect();
    spectrum(rt, &peaks, Some(&PM_IM))
}

/// PeakMapExtractor_test.cpp: `extractPeakMaps`.
#[test]
fn peak_maps_class_test() {
    let mut exp = MSExperiment::new();
    exp.spectra.push(pm_spectrum(10.0, 0.0));
    exp.spectra.push(pm_spectrum(20.0, 100.0));
    let input = access(exp);
    let coordinates = vec![
        coord(500.02, 0.0, -1.0, "tr1", 1.10),
        coord(500.02, 15.0, 25.0, "tr2", 1.10),
    ];
    let out = PeakMapExtractor::new()
        .extract_peak_maps(&input, &coordinates, 0.05, false, 0.20, "tophat")
        .unwrap();
    assert_eq!(out.len(), 2);
    assert_eq!(out[0].native_id, "tr1");
    assert_eq!(out[0].mz.len(), 6);
    assert_eq!(out[0].mz.len(), out[0].rt.len());
    assert_eq!(out[0].mz.len(), out[0].ion_mobility.len());
    assert_eq!(out[0].mz.len(), out[0].intensity.len());
    real_similar(out[0].rt[0], 10.0);
    real_similar(out[0].rt[3], 20.0);
    real_similar(out[0].mz[0], 500.00);
    real_similar(out[0].mz[2], 500.03);
    real_similar(out[0].intensity[0], 11.0);
    real_similar(out[0].intensity[5], 113.0);
    assert_eq!(out[1].native_id, "tr2");
    assert_eq!(out[1].mz.len(), 3);
    real_similar(out[1].rt[0], 20.0);
    real_similar(out[1].intensity[0], 111.0);
    // Derived: the target fields of the source's pre-extraction fill.
    assert!(out[0].target_rt.is_nan());
    assert_eq!(out[1].target_rt, 20.0);
    assert_eq!(out[1].target_ion_mobility, 1.10);
}

/// PeakMapExtractor_test.cpp: `[EXTRA] extractPeakMaps throws without ion
/// mobility data`.
#[test]
fn peak_maps_require_ion_mobility() {
    let mut exp = MSExperiment::new();
    exp.spectra.push(spectrum(10.0, &[(500.0, 100.0)], None));
    let input = access(exp);
    let coordinates = vec![coord(500.0, 0.0, 0.0, "tr1", 1.0)];
    assert!(matches!(
        PeakMapExtractor::new().extract_peak_maps(
            &input,
            &coordinates,
            0.05,
            false,
            0.20,
            "tophat"
        ),
        Err(Error::InvalidValue(_))
    ));
}

/// PeakMapExtractor_test.cpp: `[EXTRA] extractPeakMaps includes peaks on mz and
/// IM window boundaries`.
#[test]
fn peak_maps_closed_windows() {
    let mut exp = MSExperiment::new();
    exp.spectra.push(pm_spectrum(10.0, 0.0));
    let input = access(exp);
    let coordinates = vec![coord(500.02, 0.0, -1.0, "tr1", 1.10)];
    let out = PeakMapExtractor::new()
        .extract_peak_maps(&input, &coordinates, 0.06, false, 0.20, "tophat")
        .unwrap();
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].mz.len(), 4);
    real_similar(*out[0].mz.first().unwrap(), 499.99);
    real_similar(*out[0].mz.last().unwrap(), 500.03);
    real_similar(*out[0].ion_mobility.first().unwrap(), 1.00);
    real_similar(*out[0].ion_mobility.last().unwrap(), 1.15);
}

/// Derived: the filter is checked only with spectra and coordinates present,
/// a negative IM window records the raw ion mobility without filtering, and
/// mismatched arrays are refused.
#[test]
fn peak_maps_contract() {
    let x = PeakMapExtractor::new();
    let empty = access(MSExperiment::new());
    let coordinates = vec![coord(500.02, 1.0, 3.0, "tr1", 1.10)];
    let out = x
        .extract_peak_maps(&empty, &coordinates, 0.05, false, 0.2, "bartlett")
        .unwrap();
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].target_rt, 2.0);
    assert!(out[0].is_empty());

    let mut exp = MSExperiment::new();
    exp.spectra.push(pm_spectrum(2.0, 0.0));
    let input = access(exp);
    assert!(matches!(
        x.extract_peak_maps(&input, &coordinates, 0.05, false, 0.2, "bartlett"),
        Err(Error::InvalidValue(_))
    ));
    let out = x
        .extract_peak_maps(&input, &coordinates, 0.05, false, -1.0, "tophat")
        .unwrap();
    assert_eq!(out[0].len(), 3);
    assert_eq!(out[0].ion_mobility[0], f64::from(1.05f32));

    let mut exp = MSExperiment::new();
    exp.spectra
        .push(spectrum(2.0, &[(500.0, 1.0), (500.01, 2.0)], Some(&[1.0])));
    let input = access(exp);
    assert!(matches!(
        x.extract_peak_maps(&input, &coordinates, 0.05, false, 0.2, "tophat"),
        Err(Error::InvalidValue(_))
    ));
    assert!(ExtractedPeakMap::default().target_rt.is_nan());
}

// ---------------------------------------------------------------------------
// SwathWindowLoader
// ---------------------------------------------------------------------------

fn swath_maps() -> Vec<SwathMap> {
    vec![
        SwathMap::new(399.0, 426.0, 1.0, false),
        SwathMap::new(424.0, 451.0, 2.0, false),
        SwathMap::new(450.0, 475.0, 3.0, false),
        SwathMap::new(474.0, 501.0, 4.0, false),
    ]
}

fn check_annotated(maps: &[SwathMap]) {
    for (map, (lower, upper)) in maps.iter().zip([
        (400.0, 425.0),
        (425.0, 450.0),
        (450.0, 475.0),
        (475.0, 500.0),
    ]) {
        real_similar(map.lower, lower);
        real_similar(map.upper, upper);
    }
}

/// SwathWindowLoader_test.cpp: `readSwathWindows`, with and without a
/// header, and the missing-file exception.
#[test]
fn read_swath_windows_class_test() {
    let (lower, upper) =
        read_swath_windows(data("openswath_extraction_SwathWindowFile.txt")).unwrap();
    assert!(matches!(
        read_swath_windows("/does/not/exist_swath_windows.txt"),
        Err(Error::Io(ref e)) if e.kind() == std::io::ErrorKind::NotFound
    ));
    assert_eq!(lower.len(), upper.len());
    for (i, v) in [400.0, 425.0, 450.0, 475.0].into_iter().enumerate() {
        real_similar(lower[i], v);
        real_similar(upper[i], v + 25.0);
    }
    let (lower_nh, upper_nh) =
        read_swath_windows(data("openswath_extraction_SwathWindowFileNoHeader.txt")).unwrap();
    assert_eq!(lower_nh.len(), upper.len());
    assert_eq!(lower_nh, lower);
    assert_eq!(upper_nh, upper);
}

/// SwathWindowLoader_test.cpp: `annotateSwathMapsFromFile`, in order, sorted,
/// the count and order exceptions, and `force`.
#[test]
fn annotate_swath_maps_class_test() {
    let file = data("openswath_extraction_SwathWindowFile.txt");
    let mut maps = swath_maps();
    let warnings = annotate_swath_maps_from_file(&file, &mut maps, false, false).unwrap();
    assert!(warnings.is_empty());
    check_annotated(&maps);

    let mut inverted = swath_maps();
    inverted.reverse();
    let mut maps = inverted.clone();
    assert_eq!(maps[0].lower, 474.0);
    assert_eq!(maps[0].center, 4.0);
    assert_eq!(maps[3].lower, 399.0);
    assert_eq!(maps[3].center, 1.0);
    annotate_swath_maps_from_file(&file, &mut maps, true, false).unwrap();
    check_annotated(&maps);
    for (i, map) in maps.iter().enumerate() {
        real_similar(map.center, (i + 1) as f64);
    }

    let mut too_large = vec![SwathMap::default(); 5];
    assert!(matches!(
        annotate_swath_maps_from_file(&file, &mut too_large, true, false),
        Err(Error::InvalidValue(_))
    ));
    let mut too_small = vec![SwathMap::default(); 3];
    assert!(matches!(
        annotate_swath_maps_from_file(&file, &mut too_small, true, false),
        Err(Error::InvalidValue(_))
    ));

    // Wrong order and no sorting fails; the maps are left unchanged.
    let mut maps = inverted.clone();
    assert!(matches!(
        annotate_swath_maps_from_file(&file, &mut maps, false, false),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(maps[0].lower, 474.0);
    assert_eq!(maps[0].center, 4.0);

    // Wrong order, no sorting, force: overwritten, still in reverse order.
    let mut maps = inverted;
    let warnings = annotate_swath_maps_from_file(&file, &mut maps, false, true).unwrap();
    assert!(!warnings.is_empty());
    check_annotated(&maps);
    for (i, map) in maps.iter().enumerate() {
        real_similar(map.center, (4 - i) as f64);
    }
}

/// Derived: MS1 maps are skipped, both ways a count mismatch is detected, and
/// malformed rows are refused rather than repeating the previous row.
#[test]
fn swath_windows_native_contract() {
    let mut maps = vec![SwathMap::new(0.0, 2000.0, 0.0, true)];
    maps.extend(swath_maps());
    let (lower, upper) =
        parse_swath_windows("400 425\n425\t450\n450 475\n475 500 extra\n").unwrap();
    annotate_swath_maps(&lower, &upper, &mut maps, false, false).unwrap();
    assert_eq!((maps[0].lower, maps[0].upper), (0.0, 2000.0));
    check_annotated(&maps[1..]);
    assert!(matches!(
        annotate_swath_maps(&[1.0], &[], &mut maps, false, false),
        Err(Error::InvalidValue(_))
    ));

    assert!(matches!(
        parse_swath_windows("400 425\n\n425 450\n"),
        Err(Error::Parse { line: 2, .. })
    ));
    assert!(matches!(
        parse_swath_windows("400 425\n425 450\n\n"),
        Err(Error::Parse { line: 3, .. })
    ));
    assert!(matches!(parse_swath_windows(""), Err(Error::Parse { .. })));
    assert!(matches!(
        parse_swath_windows("lower upper\n"),
        Err(Error::Parse { line: 2, .. })
    ));
    assert!(matches!(
        parse_swath_windows("400\n"),
        Err(Error::Parse { line: 1, .. })
    ));
    assert!(matches!(
        parse_swath_windows("425 400\n"),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        parse_swath_windows("400 425\n400 abc\n"),
        Err(Error::Parse { line: 2, .. })
    ));
    let (lower, upper) = parse_swath_windows("  4e2\t425\r\n425 450\r\n").unwrap();
    assert_eq!((lower, upper), (vec![400.0, 425.0], vec![425.0, 450.0]));
}

// ---------------------------------------------------------------------------
// ChromatogramExtractor
// ---------------------------------------------------------------------------

/// ChromatogramExtractor_input.TraML as the source's `convertTargetedExp`
/// turns it into a light experiment. TraMLFile and `convertTargetedExp` are
/// not ported yet, so the conversion is transcribed; the oracle prints the
/// converted compounds and transitions (see the provenance manifest).
fn traml_light_experiment() -> LightTargetedExperiment {
    let mut exp = LightTargetedExperiment::new();
    for (id, seq, rt) in [("tr_gr1", "PEPTIDEA", 44.0), ("tr_gr2", "PEPTIDEB", 2.0)] {
        exp.compounds.push(LightCompound {
            id: id.into(),
            sequence: seq.into(),
            rt,
            ..LightCompound::default()
        });
    }
    for (name, pep, prec, prod, lib) in [
        ("tr1", "tr_gr1", 500.0, 628.45, 1.0),
        ("tr2", "tr_gr1", 500.0, 654.38, 2.0),
        ("tr3", "tr_gr2", 501.0, 618.31, 10000.0),
    ] {
        exp.transitions.push(LightTransition {
            transition_name: name.into(),
            peptide_ref: pep.into(),
            precursor_mz: prec,
            product_mz: prod,
            library_intensity: lib,
            ..LightTransition::default()
        });
    }
    exp
}

/// ChromatogramExtractor_test.cpp: `prepare_coordinates`, transitions.
#[test]
fn prepare_coordinates_transitions() {
    let exp = traml_light_experiment();
    let mut chroms = Vec::new();
    let mut coords = Vec::new();
    prepare_coordinates(&mut chroms, &mut coords, &exp, 1.0, false, 0).unwrap();
    assert_eq!(chroms.len(), coords.len());
    assert_eq!(coords.len(), 3);
    assert_eq!(coords[0].mz, 618.31);
    assert_eq!(coords[1].mz, 628.45);
    assert_eq!(coords[2].mz, 654.38);
    real_similar(coords[0].rt_start, 1.5);
    real_similar(coords[1].rt_start, 43.5);
    real_similar(coords[2].rt_start, 43.5);
    real_similar(coords[0].rt_end, 2.5);
    real_similar(coords[1].rt_end, 44.5);
    real_similar(coords[2].rt_end, 44.5);
    assert_eq!(coords[0].id, "tr3");
    assert_eq!(coords[1].id, "tr1");
    assert_eq!(coords[2].id, "tr2");
    // Derived: MS2 fields.
    assert_eq!(coords[0].mz_precursor, 501.0);
    assert_eq!(coords[0].ion_mobility, -1.0);
}

/// ChromatogramExtractor_test.cpp: `prepare_coordinates`, peptides, and the
/// issue #7284 isotope regression.
#[test]
fn prepare_coordinates_peptides_and_isotopes() {
    let exp = traml_light_experiment();
    let mut chroms = Vec::new();
    let mut coords = Vec::new();
    prepare_coordinates(&mut chroms, &mut coords, &exp, 1.0, true, 0).unwrap();
    assert_eq!(chroms.len(), coords.len());
    assert_eq!(coords.len(), 2);
    assert_eq!(coords[0].mz, 500.0);
    assert_eq!(coords[1].mz, 501.0);
    real_similar(coords[0].rt_start, 43.5);
    real_similar(coords[1].rt_start, 1.5);
    real_similar(coords[0].rt_end, 44.5);
    real_similar(coords[1].rt_end, 2.5);
    assert_eq!(coords[0].id, "tr_gr1_Precursor_i0");
    assert_eq!(coords[1].id, "tr_gr2_Precursor_i0");
    assert_eq!(compute_transition_group_id(&coords[0].id), "tr_gr1");
    assert_eq!(compute_transition_group_id(&coords[1].id), "tr_gr2");

    let mut chroms = Vec::new();
    let mut coords = Vec::new();
    prepare_coordinates(&mut chroms, &mut coords, &exp, 1.0, true, 2).unwrap();
    assert_eq!(chroms.len(), coords.len());
    assert_eq!(coords.len(), 6);
    let d = C13C12_MASSDIFF_U;
    for (c, (mz, id)) in coords.iter().zip([
        (500.0, "tr_gr1_Precursor_i0"),
        (501.0, "tr_gr2_Precursor_i0"),
        (500.0 + d, "tr_gr1_Precursor_i1"),
        (501.0 + d, "tr_gr2_Precursor_i1"),
        (500.0 + 2.0 * d, "tr_gr1_Precursor_i2"),
        (501.0 + 2.0 * d, "tr_gr2_Precursor_i2"),
    ]) {
        real_similar(c.mz, mz);
        assert_eq!(c.id, id);
    }
}

/// Derived: the three RT-window modes, a compound without transitions, a
/// dangling compound reference, and atomicity.
#[test]
fn prepare_coordinates_contract() {
    let mut exp = traml_light_experiment();
    let mut chroms = Vec::new();
    let mut coords = Vec::new();
    // Negative window: whole RT range.
    prepare_coordinates(&mut chroms, &mut coords, &exp, -1.0, false, 0).unwrap();
    assert!(coords.iter().all(|c| c.rt_start == 0.0 && c.rt_end == -1.0));
    // NaN window without RT ranges: refused, nothing appended.
    assert!(matches!(
        prepare_coordinates(&mut chroms, &mut coords, &exp, f64::NAN, false, 0),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!((chroms.len(), coords.len()), (3, 3));
    // NaN window with RT ranges: copied.
    for c in &mut exp.compounds {
        c.rt_start = 10.0;
        c.rt_end = 20.0;
    }
    let mut chroms = Vec::new();
    let mut coords = Vec::new();
    prepare_coordinates(&mut chroms, &mut coords, &exp, f64::NAN, false, 0).unwrap();
    assert!(
        coords
            .iter()
            .all(|c| c.rt_start == 10.0 && c.rt_end == 20.0)
    );
    // RT-limited without RT: refused.
    exp.compounds[0].rt = f64::NAN;
    assert!(matches!(
        prepare_coordinates(&mut chroms, &mut coords, &exp, 5.0, true, 0),
        Err(Error::InvalidValue(_))
    ));
    // A compound without transitions gets one coordinate at m/z 0, no isotopes.
    let mut exp = traml_light_experiment();
    exp.compounds.push(LightCompound {
        id: "lonely".into(),
        rt: 7.0,
        drift_time: 3.0,
        ..LightCompound::default()
    });
    let mut chroms = Vec::new();
    let mut coords = Vec::new();
    prepare_coordinates(&mut chroms, &mut coords, &exp, 1.0, true, 1).unwrap();
    assert_eq!(coords.len(), 5);
    assert_eq!(coords[0].id, "lonely_Precursor_i0");
    assert_eq!(
        (
            coords[0].mz,
            coords[0].rt_start,
            coords[0].rt_end,
            coords[0].ion_mobility
        ),
        (0.0, 0.0, -1.0, 0.0)
    );
    // A transition referencing an unknown compound.
    exp.transitions[0].peptide_ref = "missing".into();
    assert!(matches!(
        prepare_coordinates(&mut chroms, &mut coords, &exp, 1.0, false, 0),
        Err(Error::MissingInformation(_))
    ));
}

/// ChromatogramExtractor_test.cpp: `return_chromatogram` after extraction
/// over ChromatogramExtractor_input.mzML.
#[cfg(feature = "mzml")]
#[test]
fn return_chromatogram_class_test() {
    let bytes = std::fs::read(data(
        "openswath_extraction_ChromatogramExtractor_input.mzML",
    ))
    .unwrap();
    let exp = Arc::new(openms::format::mzml::read(bytes.as_slice()).unwrap());
    let input = SpectrumAccessOpenMS::new(Arc::clone(&exp));
    let light = traml_light_experiment();
    let mut output = Vec::new();
    let mut coords = Vec::new();
    prepare_coordinates(&mut output, &mut coords, &light, -1.0, false, 0).unwrap();
    ChromatogramExtractor::new()
        .extract_chromatograms(&input, &mut output, &coords, 0.05, false, "tophat")
        .unwrap();
    let mut chromatograms = Vec::new();
    return_chromatogram(
        &output,
        &coords,
        &light,
        &exp.spectra[0],
        &mut chromatograms,
        false,
        0.0,
    )
    .unwrap();
    assert_eq!(chromatograms.len(), 3);
    assert_eq!(
        chromatograms[0].chromatogram_type,
        ChromatogramType::SelectedReactionMonitoring
    );
    real_similar(chromatograms[0].product.mz, 618.31);
    assert!(
        chromatograms[0]
            .precursor
            .cv_terms
            .metadata
            .contains_key("peptide_sequence")
    );
    // Derived: the rest of the annotation.
    assert_eq!(chromatograms[0].native_id, "tr3");
    assert_eq!(chromatograms[0].precursor.mz, 501.0);
    assert_eq!(chromatograms[0].peaks.len(), 59);
    assert_eq!(
        chromatograms[0].precursor.cv_terms.metadata["peptide_sequence"].to_string(),
        "PEPTIDEB"
    );
    assert_eq!(
        chromatograms[0].precursor.isolation_window_lower_offset,
        exp.spectra[0].precursors[0].isolation_window_lower_offset
    );
}

/// Derived: MS1 conversion, ion-mobility annotation, data-processing marks and
/// the error paths.
#[test]
fn return_chromatogram_contract() {
    let light = traml_light_experiment();
    let mut output = Vec::new();
    let mut coords = Vec::new();
    prepare_coordinates(&mut output, &mut coords, &light, -1.0, true, 1).unwrap();
    coords[0].ion_mobility = 0.9;
    let mut settings = MSSpectrum::default();
    settings
        .data_processing
        .push(Arc::new(openms::metadata::DataProcessing::default()));
    let mut chromatograms = Vec::new();
    return_chromatogram(
        &output,
        &coords,
        &light,
        &settings,
        &mut chromatograms,
        true,
        0.2,
    )
    .unwrap();
    assert_eq!(chromatograms.len(), 4);
    let c = &chromatograms[0];
    assert_eq!(c.chromatogram_type, ChromatogramType::BasePeak);
    assert_eq!(c.native_id, "tr_gr1_Precursor_i0");
    assert_eq!(c.precursor.mz, 500.0);
    assert_eq!(c.precursor.drift_time, Some(0.9));
    assert_eq!(c.precursor.drift_window_lower_offset, 0.1);
    assert_eq!(
        c.precursor.cv_terms.metadata["peptide_sequence"].to_string(),
        "PEPTIDEA"
    );
    assert!(
        c.data_processing[0]
            .metadata
            .contains_key("performed_on_spectra")
    );
    assert!(settings.data_processing[0].metadata.is_empty());
    // Coordinates with no drift time (-1) get none.
    assert_eq!(chromatograms[1].precursor.drift_time, None);

    // MS2 with an unknown transition id, and too few coordinates.
    let mut bad = coords.clone();
    bad[0].id = "nope".into();
    let before = chromatograms.len();
    assert!(matches!(
        return_chromatogram(
            &output,
            &bad,
            &light,
            &settings,
            &mut chromatograms,
            false,
            0.0
        ),
        Err(Error::MissingInformation(_))
    ));
    assert!(matches!(
        return_chromatogram(
            &output,
            &coords[..1],
            &light,
            &settings,
            &mut chromatograms,
            true,
            0.0
        ),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(chromatograms.len(), before);
}

// ---------------------------------------------------------------------------
// SwathQC
// ---------------------------------------------------------------------------

/// SwathQC_test.cpp: `isSubsampledSpectrum_`.
#[test]
fn swath_qc_subsample_class_test() {
    let f = SwathQC::is_subsampled_spectrum;
    assert!(f(0, 100, 4));
    assert!(f(10, 100, 4));
    assert!(!f(10, 4, 10));
    assert!(!f(10, 4, 11));
    let r = [1, 0, 0, 1, 0, 1, 0, 0, 1, 0];
    for (i, &expected) in r.iter().enumerate() {
        assert_eq!(f(10, 4, i), expected == 1, "10/4 at {i}");
    }
    for i in 0..10 {
        assert!(!f(10, 0, i));
    }
    for i in 0..4 {
        assert!(f(4, 4, i));
    }
    let r5 = [1, 0, 0, 1, 0];
    for (i, &expected) in r5.iter().enumerate() {
        assert_eq!(f(5, 2, i), expected == 1, "5/2 at {i}");
    }
}

#[cfg(feature = "mzml")]
fn orbitrap_maps() -> Vec<SwathMap> {
    let bytes = std::fs::read(data(
        "openswath_extraction_PeakPickerHiRes_orbitrap_sn1_out.mzML",
    ))
    .unwrap();
    let exp = openms::format::mzml::read(bytes.as_slice()).unwrap();
    vec![SwathMap {
        sptr: Some(Arc::new(access(exp))),
        ms1: true,
        ..SwathMap::default()
    }]
}

/// SwathQC_test.cpp: static `getChargeDistribution` over
/// PeakPickerHiRes_orbitrap_sn1_out.mzML, and `storeJSON` after the
/// processing callback.
#[cfg(feature = "mzml")]
#[test]
fn swath_qc_class_test() {
    let maps = orbitrap_maps();
    let cd = SwathQC::charge_distribution_of_maps(&maps, 10, 0.04).unwrap();
    let expected: ChargeDistribution = [(1, 12), (6, 1), (10, 4)].into_iter().collect();
    assert_eq!(cd, expected);

    let bytes = std::fs::read(data(
        "openswath_extraction_PeakPickerHiRes_orbitrap_sn1_out.mzML",
    ))
    .unwrap();
    let exp = openms::format::mzml::read(bytes.as_slice()).unwrap();
    let mut qc = SwathQC::new(10, 0.04);
    let count = exp.spectra.iter().filter(|s| s.ms_level == 1).count();
    qc.set_nr_ms1_spectra(count);
    {
        let mut f = qc.spectra_processing_func();
        for s in exp.spectra.iter().filter(|s| s.ms_level == 1) {
            f(s).unwrap();
        }
    }
    let dir = std::env::temp_dir().join(format!("openswath_qc_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("qc.json");
    qc.store_json(&path).unwrap();
    let written = std::fs::read_to_string(&path).unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
    // The class test compares against this document with FuzzyStringComparator.
    let expected_json = "{\n  \"ChargeDistributionMS1\": [\n    [\n      1,\n      12\n    ],\n    [\n      6,\n      1\n    ],\n    [\n      10,\n      4\n    ]\n  ]\n}\n";
    assert_eq!(written, expected_json);
}

/// Derived: experimental settings, the MS-level filter and the empty JSON.
#[test]
fn swath_qc_contract() {
    let mut qc = SwathQC::new(5, 0.04);
    let mut es = ExperimentalSettings::default();
    qc.apply_experimental_settings(&es).unwrap();
    es.metadata
        .insert("nr_ms1_spectra".into(), MetaValue::from(-3_i64));
    assert!(matches!(
        qc.apply_experimental_settings(&es),
        Err(Error::InvalidValue(_))
    ));
    es.metadata
        .insert("nr_ms1_spectra".into(), MetaValue::from("7"));
    assert!(matches!(
        qc.exp_settings_func()(&es),
        Err(Error::InvalidValue(_))
    ));
    es.metadata
        .insert("nr_ms1_spectra".into(), MetaValue::from(7_i64));
    qc.apply_experimental_settings(&es).unwrap();
    // Once set, later settings are ignored.
    es.metadata
        .insert("nr_ms1_spectra".into(), MetaValue::from("x"));
    qc.apply_experimental_settings(&es).unwrap();
    let ms2 = MSSpectrum {
        ms_level: 2,
        ..MSSpectrum::default()
    };
    qc.process_spectrum(&ms2).unwrap();
    assert!(qc.charge_distribution().is_empty());
    assert_eq!(qc.to_json(), "{\n  \"ChargeDistributionMS1\": []\n}\n");
    let ms1 = SwathMap {
        ms1: true,
        ..SwathMap::default()
    };
    assert!(matches!(
        SwathQC::charge_distribution_of_maps(&[ms1], 1, 0.04),
        Err(Error::MissingInformation(_))
    ));
}

// ---------------------------------------------------------------------------
// SwathMapMassCorrection
// ---------------------------------------------------------------------------

/// Derived: parameters round-trip, restrictions are enforced, accessors.
#[test]
fn swath_map_mass_correction_params() {
    let defaults = SwathMapMassCorrection::defaults().unwrap();
    let smc = SwathMapMassCorrection::from_param(&defaults).unwrap();
    assert_eq!(smc, SwathMapMassCorrection::new());
    assert_eq!(smc.to_param().unwrap(), defaults);
    let mut p = defaults.clone();
    p.set_value(
        "mz_correction_function",
        ParamValue::String("quadratic_regression".into()),
        "",
        &[],
    )
    .unwrap();
    p.set_value(
        "mz_extraction_window_ppm",
        ParamValue::String("true".into()),
        "",
        &[],
    )
    .unwrap();
    let smc = SwathMapMassCorrection::from_param(&p).unwrap();
    assert_eq!(smc.mz_correction_function, "quadratic_regression");
    assert!(smc.mz_extraction_window_ppm);
    let mut bad = defaults.clone();
    bad.set_value(
        "mz_correction_function",
        ParamValue::String("cubic".into()),
        "",
        &[],
    )
    .unwrap();
    assert!(SwathMapMassCorrection::from_param(&bad).is_err());
    let mut bad = defaults;
    bad.set_value(
        "mz_estimation_percentile",
        ParamValue::Float(100.0),
        "",
        &[],
    )
    .unwrap();
    assert!(SwathMapMassCorrection::from_param(&bad).is_err());

    let mut smc = SwathMapMassCorrection::new();
    assert_eq!(
        (
            smc.fragment_mz_window(),
            smc.fragment_im_window(),
            smc.precursor_mz_window(),
            smc.precursor_im_window()
        ),
        (-1.0, -1.0, -1.0, -1.0)
    );
    smc.set_fragment_mz_window(1.0);
    smc.set_fragment_im_window(2.0);
    smc.set_precursor_mz_window(3.0);
    smc.set_precursor_im_window(4.0);
    assert_eq!(
        (
            smc.fragment_mz_window(),
            smc.fragment_im_window(),
            smc.precursor_mz_window(),
            smc.precursor_im_window()
        ),
        (1.0, 2.0, 3.0, 4.0)
    );
}

fn pasef_map(lower: f64, upper: f64, im_lower: f64, im_upper: f64) -> SwathMap {
    SwathMap::with_im(
        lower,
        upper,
        (lower + upper) / 2.0,
        im_lower,
        im_upper,
        false,
    )
}

fn pasef_group(im: f64) -> MRMTransitionGroupType {
    let mut tg = MRMTransitionGroupType::new();
    tg.add_transition(
        LightTransition {
            transition_name: "t".into(),
            precursor_mz: 412.0,
            precursor_im: im,
            product_mz: 500.0,
            ..LightTransition::default()
        },
        "t",
    )
    .unwrap();
    tg
}

/// Oracle: `findSwathMapsPasef` picks the map whose IM centre is closest,
/// the first on a tie; guards derived.
#[test]
fn swath_map_mass_correction_pasef() {
    let smc = SwathMapMassCorrection::new();
    let maps = vec![
        pasef_map(400.0, 425.0, 5.0, 12.0),
        pasef_map(400.0, 425.0, 9.0, 14.0),
        pasef_map(400.0, 425.0, 10.0, 12.5),
        pasef_map(412.0, 425.0, 0.0, 20.0),
        pasef_map(400.0, 412.0, 0.0, 20.0),
        pasef_map(400.0, 425.0, 11.0, 20.0),
    ];
    let used = smc
        .find_swath_maps_pasef(&pasef_group(11.0), &maps)
        .unwrap();
    assert_eq!(used.len(), 1);
    assert_eq!(
        (
            used[0].lower,
            used[0].upper,
            used[0].im_lower,
            used[0].im_upper
        ),
        (400.0, 425.0, 10.0, 12.5)
    );
    let tie = vec![
        pasef_map(400.0, 425.0, 9.0, 13.0),
        pasef_map(400.0, 425.0, 10.0, 12.0),
    ];
    let used = smc.find_swath_maps_pasef(&pasef_group(11.0), &tie).unwrap();
    assert_eq!((used[0].im_lower, used[0].im_upper), (9.0, 13.0));
    assert!(
        smc.find_swath_maps_pasef(&pasef_group(11.0), &[])
            .unwrap()
            .is_empty()
    );
    assert!(matches!(
        smc.find_swath_maps_pasef(&pasef_group(-1.0), &maps),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        smc.find_swath_maps_pasef(&MRMTransitionGroupType::new(), &maps),
        Err(Error::InvalidValue(_))
    ));
}

// ---------------------------------------------------------------------------
// Tier 1: replay of the C++ Release oracle (ibminode06)
// ---------------------------------------------------------------------------

/// `(section, key, rest)` records of tests/data/openswath_extraction_oracle.tsv.
fn oracle() -> Vec<(String, String, String)> {
    std::fs::read_to_string(data("openswath_extraction_oracle.tsv"))
        .unwrap()
        .lines()
        .filter_map(|line| {
            let mut parts = line.splitn(3, '\t');
            Some((
                parts.next()?.to_string(),
                parts.next()?.to_string(),
                parts.next().unwrap_or("").to_string(),
            ))
        })
        .collect()
}

fn oracle_get(records: &[(String, String, String)], section: &str, key: &str) -> String {
    records
        .iter()
        .find(|(s, k, _)| s == section && k == key)
        .unwrap_or_else(|| panic!("no oracle record {section}/{key}"))
        .2
        .clone()
}

/// The C++ `%.17g` of a double; the parse is exact.
fn g17(text: &str) -> f64 {
    match text {
        "nan" => f64::NAN,
        "-nan" => -f64::NAN,
        _ => text.parse().unwrap(),
    }
}

fn g17_list(text: &str) -> Vec<f64> {
    if text.is_empty() {
        Vec::new()
    } else {
        text.split(',').map(g17).collect()
    }
}

fn same(a: f64, b: f64) -> bool {
    a == b || (a.is_nan() && b.is_nan())
}

/// Oracle: every `extract_value_tophat` call of the class-test sequences and
/// the quirk cases, value and walk position, bit for bit.
#[test]
fn oracle_tophat() {
    let records = oracle();
    let x = ChromatogramExtractorAlgorithm::new();
    type Call = (f64, f64, bool, f64, f64);
    let run = |name: &str, mz: &[f64], int: &[f64], im: Option<&[f64]>, calls: &[Call]| {
        let mut pos = 0usize;
        for (c, &(target, w, ppm, tim, imw)) in calls.iter().enumerate() {
            if target < 0.0 {
                pos = 0;
                continue;
            }
            let v = match im {
                Some(im) => x
                    .extract_value_tophat_im(mz, int, im, &mut pos, target, tim, w, imw, ppm)
                    .unwrap(),
                None => x
                    .extract_value_tophat(mz, int, &mut pos, target, w, ppm)
                    .unwrap(),
            };
            let rec = oracle_get(&records, "tophat", &format!("{name}#{c}"));
            let (value, p) = rec.split_once("\tpos=").unwrap();
            assert!(same(v, g17(value)), "{name}#{c}: {v} vs {value}");
            assert_eq!(pos.to_string(), p, "{name}#{c} position");
        }
    };
    let r = (-1.0, 0.0, false, 0.0, 0.0);
    let th: Vec<Call> = [
        (399.805, 0.2, false),
        (399.91, 0.2, false),
        (400.0, 0.2, false),
        (400.05, 0.2, false),
        (400.1, 0.2, false),
        (400.28, 0.2, false),
        (500.0, 0.2, false),
        (400.1, 0.2, false),
        (-1.0, 0.0, false),
        (399.89, 500.0, true),
        (399.91, 500.0, true),
        (399.92, 500.0, true),
        (400.0, 500.0, true),
        (400.05, 500.0, true),
        (400.1, 500.0, true),
        (-1.0, 0.0, false),
        (399.805, 0.2, false),
        (400.0001, 0.001, false),
    ]
    .iter()
    .map(|&(t, w, p)| (t, w, p, 0.0, 0.0))
    .collect();
    run("th", &MZ_ARR, &INT_ARR, None, &th);
    let im = |t: f64, w: f64, tim: f64| (t, w, false, tim, 0.3);
    run(
        "th_im",
        &MZ_ARR,
        &INT_ARR,
        Some(&IM_ARR),
        &[
            im(399.805, 0.2, 100.0),
            im(399.91, 0.2, 100.0),
            im(400.0, 0.2, 100.0),
            im(400.05, 0.2, 100.0),
            im(400.1, 0.2, 100.0),
            im(400.28, 0.2, 100.0),
            im(400.28, 0.2, 200.0),
            im(400.28, 0.2, 200.1),
            im(500.0, 0.2, 300.0),
            im(500.0, 0.2, 300.1),
            im(400.1, 0.2, 100.0),
            r,
            im(399.805, 0.2, 100.0),
            im(400.0001, 0.001, 100.0),
        ],
    );
    run(
        "quirk_end",
        &[100.0, 200.0, 200.02],
        &[1.0, 100.0, 7.0],
        None,
        &[(200.05, 0.2, false, 0.0, 0.0)],
    );
    run(
        "quirk_first",
        &[100.0, 100.05, 100.1],
        &[1.0, 10.0, 100.0],
        None,
        &[(100.1, 0.4, false, 0.0, 0.0)],
    );
    run(
        "quirk_end_im",
        &[100.0, 200.0, 200.02],
        &[1.0, 100.0, 7.0],
        Some(&[5.0, 5.0, 5.0]),
        &[(200.05, 0.2, false, 5.0, 1.0)],
    );
    run(
        "single",
        &[100.0],
        &[3.0],
        None,
        &[
            (100.0, 0.2, false, 0.0, 0.0),
            (100.05, 0.2, false, 0.0, 0.0),
        ],
    );
}

#[cfg(feature = "mzml")]
fn check_chroms(records: &[(String, String, String)], name: &str, out: &[Chromatogram]) {
    for (k, c) in out.iter().enumerate() {
        let rt = g17_list(&oracle_get(records, "xic", &format!("{name}#{k}.rt")));
        let int = g17_list(&oracle_get(records, "xic", &format!("{name}#{k}.int")));
        assert_eq!(c.time_array().data, rt, "{name}#{k} RT");
        assert_eq!(c.intensity_array().data, int, "{name}#{k} intensity");
    }
}

/// Oracle: every extracted point over ChromatogramExtractor_input.mzML and
/// over the ion-mobility experiment at three IM windows, bit for bit.
#[cfg(feature = "mzml")]
#[test]
fn oracle_extract_chromatograms() {
    let records = oracle();
    let bytes = std::fs::read(data(
        "openswath_extraction_ChromatogramExtractor_input.mzML",
    ))
    .unwrap();
    let input = access(openms::format::mzml::read(bytes.as_slice()).unwrap());
    let coordinates = vec![
        coord(618.31, 0.0, -1.0, "tr1", 0.0),
        coord(628.45, 0.0, -1.0, "tr2", 0.0),
        coord(654.38, 0.0, -1.0, "tr3", 0.0),
    ];
    let mut out = vec![Chromatogram::new(); 3];
    ChromatogramExtractorAlgorithm::new()
        .extract_chromatograms(&input, &mut out, &coordinates, 0.05, false, -1.0, "tophat")
        .unwrap();
    check_chroms(&records, "file", &out);

    let input = im_experiment();
    let coordinates = vec![
        coord(618.31, 0.0, -1.0, "tr1", 120.0),
        coord(628.45, 0.0, -1.0, "tr2", 170.0),
    ];
    for (w, name) in [(-1.0, "im-1"), (15.0, "im15"), (30.0, "im30")] {
        let mut out = vec![Chromatogram::new(); 2];
        ChromatogramExtractorAlgorithm::new()
            .extract_chromatograms(&input, &mut out, &coordinates, 0.10, false, w, "tophat")
            .unwrap();
        check_chroms(&records, name, &out);
    }
    // The C++ Bartlett run throws before appending anything (partial size 0).
    assert_eq!(oracle_get(&records, "xic", "bartlett_partial"), "0");
    assert_eq!(oracle_get(&records, "xic", "empty_input"), "no exception");
}

/// Oracle: every extracted peak-map point and target field.
#[test]
fn oracle_peak_maps() {
    let records = oracle();
    let mut exp = MSExperiment::new();
    exp.spectra.push(pm_spectrum(10.0, 0.0));
    exp.spectra.push(pm_spectrum(20.0, 100.0));
    let input = access(exp);
    let coordinates = vec![
        coord(500.02, 0.0, -1.0, "tr1", 1.10),
        coord(500.02, 15.0, 25.0, "tr2", 1.10),
    ];
    for (w, ws) in [
        (0.05, "0.050000000000000003"),
        (0.06, "0.059999999999999998"),
    ] {
        for (iw, iws) in [(0.20, "0.20000000000000001"), (-1.0, "-1")] {
            let out = PeakMapExtractor::new()
                .extract_peak_maps(&input, &coordinates, w, false, iw, "tophat")
                .unwrap();
            for (k, m) in out.iter().enumerate() {
                let key = format!("w{ws}_iw{iws}#{k}");
                let meta = oracle_get(&records, "peakmap", &format!("{key}.meta"));
                let f: Vec<&str> = meta.split('\t').collect();
                assert_eq!(f[0], m.native_id);
                assert!(same(g17(f[1]), m.target_mz));
                assert!(same(g17(f[2]), m.target_rt));
                assert!(same(g17(f[3]), m.target_ion_mobility));
                assert!(same(g17(f[4]), m.rt_start));
                assert!(same(g17(f[5]), m.rt_end));
                for (field, values) in [
                    ("mz", &m.mz),
                    ("rt", &m.rt),
                    ("im", &m.ion_mobility),
                    ("int", &m.intensity),
                ] {
                    let expected =
                        g17_list(&oracle_get(&records, "peakmap", &format!("{key}.{field}")));
                    assert_eq!(values, &expected, "{key}.{field}");
                }
            }
        }
    }
}

/// Oracle: the light experiment `convertTargetedExp` builds from the TraML
/// (checks the transcription in `traml_light_experiment`) and every
/// coordinate `prepare_coordinates` produces.
#[test]
fn oracle_prepare_coordinates() {
    let records = oracle();
    let exp = traml_light_experiment();
    let compounds: Vec<&String> = records
        .iter()
        .filter(|(s, k, _)| s == "light" && k == "compound")
        .map(|r| &r.2)
        .collect();
    assert_eq!(compounds.len(), exp.compounds.len());
    for (c, rec) in exp.compounds.iter().zip(compounds) {
        let text = format!(
            "{}\tseq={}\tname={}\trt={}\trt_start=nan\trt_end=nan\tdrift={}\tcharge={}",
            c.id, c.sequence, c.compound_name, c.rt, c.drift_time, c.charge
        );
        assert_eq!(&text, rec);
        assert!(c.rt_start.is_nan() && c.rt_end.is_nan());
    }
    let transitions: Vec<&String> = records
        .iter()
        .filter(|(s, k, _)| s == "light" && k == "transition")
        .map(|r| &r.2)
        .collect();
    for (t, rec) in exp.transitions.iter().zip(transitions) {
        let f: Vec<&str> = rec.split('\t').collect();
        assert_eq!(f[0], t.transition_name);
        assert_eq!(f[1], format!("pep={}", t.peptide_ref));
        assert_eq!(g17(&f[2][5..]), t.precursor_mz);
        assert_eq!(g17(&f[3][5..]), t.product_mz);
        assert_eq!(g17(&f[4][4..]), t.library_intensity);
        assert_eq!(g17(&f[5][8..]), t.precursor_im);
    }
    let check = |name: &str, exp: &LightTargetedExperiment, w: f64, ms1: bool, iso: i32| {
        let mut chroms = Vec::new();
        let mut coords = Vec::new();
        prepare_coordinates(&mut chroms, &mut coords, exp, w, ms1, iso).unwrap();
        assert_eq!(
            oracle_get(&records, "coord", &format!("{name}.n")),
            format!("{}\t{}", coords.len(), chroms.len())
        );
        for (i, c) in coords.iter().enumerate() {
            let rec = oracle_get(&records, "coord", &format!("{name}#{i}"));
            let f: Vec<&str> = rec.split('\t').collect();
            assert_eq!(f[0], c.id, "{name}#{i}");
            for (field, v) in
                f[1..]
                    .iter()
                    .zip([c.mz, c.mz_precursor, c.ion_mobility, c.rt_start, c.rt_end])
            {
                let (_, text) = field.split_once('=').unwrap();
                assert!(same(g17(text), v), "{name}#{i} {field} vs {v}");
            }
        }
    };
    check("ms2_w1", &exp, 1.0, false, 0);
    check("ms1_w1", &exp, 1.0, true, 0);
    check("ms1_w1_iso2", &exp, 1.0, true, 2);
    check("ms2_wneg", &exp, -1.0, false, 0);
    check("ms1_wneg_isoneg", &exp, -1.0, true, -3);
    assert!(oracle_get(&records, "coord", "ms2_wnan").starts_with("EXCEPTION IllegalArgument"));
    let mut lonely = exp.clone();
    lonely.compounds.push(LightCompound {
        id: "lonely".into(),
        rt: 7.0,
        drift_time: 3.0,
        ..LightCompound::default()
    });
    check("lonely", &lonely, 1.0, true, 1);
}

/// Oracle: the annotation `return_chromatogram` gives, MS2 and MS1.
#[cfg(feature = "mzml")]
#[test]
fn oracle_return_chromatogram() {
    let records = oracle();
    let bytes = std::fs::read(data(
        "openswath_extraction_ChromatogramExtractor_input.mzML",
    ))
    .unwrap();
    let exp = Arc::new(openms::format::mzml::read(bytes.as_slice()).unwrap());
    let input = SpectrumAccessOpenMS::new(Arc::clone(&exp));
    let light = traml_light_experiment();
    let seq = |c: &openms::kernel::MSChromatogram| {
        c.precursor
            .cv_terms
            .metadata
            .get("peptide_sequence")
            .map_or("<none>".to_string(), ToString::to_string)
    };
    let drift = |c: &openms::kernel::MSChromatogram| c.precursor.drift_time.unwrap_or(-1.0);

    let mut output = Vec::new();
    let mut coords = Vec::new();
    prepare_coordinates(&mut output, &mut coords, &light, -1.0, false, 0).unwrap();
    ChromatogramExtractor::new()
        .extract_chromatograms(&input, &mut output, &coords, 0.05, false, "tophat")
        .unwrap();
    let mut chroms = Vec::new();
    return_chromatogram(
        &output,
        &coords,
        &light,
        &exp.spectra[0],
        &mut chroms,
        false,
        0.0,
    )
    .unwrap();
    for (i, c) in chroms.iter().enumerate() {
        let p = &c.precursor;
        let text = format!(
            "{}\ttype=5\tprec={}\tcharge={}\tseq={}\tprod={}\tiso_lo={}\tiso_hi={}\tdrift={}\tnpeaks={}\tndp={}",
            c.native_id,
            p.mz,
            p.charge,
            seq(c),
            c.product.mz,
            p.isolation_window_lower_offset,
            p.isolation_window_upper_offset,
            drift(c),
            c.peaks.len(),
            c.data_processing.len()
        );
        assert_eq!(
            c.chromatogram_type,
            ChromatogramType::SelectedReactionMonitoring
        );
        let expected = oracle_get(&records, "retchrom", &format!("ms2#{i}"));
        // `%.17g` of the product m/z; Rust prints the shortest round trip.
        let f: Vec<&str> = expected.split('\t').collect();
        assert!(same(g17(&f[5][5..]), c.product.mz));
        let mut ours: Vec<String> = text.split('\t').map(str::to_string).collect();
        ours[5] = f[5].to_string();
        assert_eq!(ours.join("\t"), expected);
    }

    let mut output = Vec::new();
    let mut coords = Vec::new();
    prepare_coordinates(&mut output, &mut coords, &light, -1.0, true, 1).unwrap();
    coords[0].ion_mobility = 0.9;
    let mut chroms = Vec::new();
    return_chromatogram(
        &output,
        &coords,
        &light,
        &exp.spectra[0],
        &mut chroms,
        true,
        0.2,
    )
    .unwrap();
    for (i, c) in chroms.iter().enumerate() {
        assert_eq!(c.chromatogram_type, ChromatogramType::BasePeak);
        let expected = oracle_get(&records, "retchrom", &format!("ms1#{i}"));
        let f: Vec<&str> = expected.split('\t').collect();
        assert_eq!(f[0], c.native_id);
        assert!(same(g17(&f[2][5..]), c.precursor.mz));
        assert_eq!(f[3], format!("charge={}", c.precursor.charge));
        assert_eq!(f[4], format!("seq={}", seq(c)));
        assert!(same(g17(&f[5][5..]), c.product.mz));
        assert!(same(g17(&f[6][6..]), drift(c)));
        assert!(same(g17(&f[7][4..]), c.precursor.drift_window_lower_offset));
        assert!(same(g17(&f[8][4..]), c.precursor.drift_window_upper_offset));
        assert_eq!(f[9], format!("npeaks={}", c.peaks.len()));
    }
}

/// Oracle: the SWATH-window reader on the class-test files and on edge
/// files. Files the C++ reads with a repeated row (blank lines) or reads into
/// a range error are refused here; see the module documentation.
#[test]
fn oracle_swath_windows() {
    let records = oracle();
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let accepted = [
        (
            "SwathWindowFile.txt",
            "openswath_extraction_SwathWindowFile.txt".to_string(),
        ),
        (
            "SwathWindowFileNoHeader.txt",
            "openswath_extraction_SwathWindowFileNoHeader.txt".to_string(),
        ),
    ];
    for (key, file) in accepted {
        let (lower, upper) = read_swath_windows(dir.join(file)).unwrap();
        let rec = oracle_get(&records, "swathwin", key);
        let (lo, up) = rec.split_once('\t').unwrap();
        assert_eq!(lower, g17_list(&lo[6..]));
        assert_eq!(upper, g17_list(&up[6..]));
    }
    for (key, content) in [
        ("swath_crlf.txt", "lower upper\r\n400 425\r\n425 450\r\n"),
        ("swath_extra_cols.txt", "400 425 999\n425 450 x\n"),
    ] {
        let (lower, upper) = parse_swath_windows(content).unwrap();
        let rec = oracle_get(&records, "swathwin", key);
        let (lo, up) = rec.split_once('\t').unwrap();
        assert_eq!(lower, g17_list(&lo[6..]));
        assert_eq!(upper, g17_list(&up[6..]));
    }
    // C++ repeats the previous row for a blank line; refused here.
    assert_eq!(
        oracle_get(&records, "swathwin", "swath_blank_trailing.txt"),
        "lower=400,425,425\tupper=425,450,450"
    );
    assert!(parse_swath_windows("400 425\n425 450\n\n").is_err());
    assert_eq!(
        oracle_get(&records, "swathwin", "swath_blank_middle.txt"),
        "lower=400,400,425\tupper=425,425,450"
    );
    assert!(parse_swath_windows("400 425\n\n425 450\n").is_err());
    // Rows the C++ refuses as illegal ranges are refused here too.
    for (key, content) in [
        ("swath_one_col.txt", "400 425\n430\n"),
        ("swath_inverted.txt", "400 425\n425 400\n"),
        ("swath_garbage.txt", "400 425\n425 abc\n"),
    ] {
        assert!(oracle_get(&records, "swathwin", key).starts_with("EXCEPTION InvalidValue"));
        assert!(parse_swath_windows(content).is_err(), "{key}");
    }
    // Annotation.
    let file = data("openswath_extraction_SwathWindowFile.txt");
    let render = |maps: &[SwathMap]| {
        maps.iter()
            .map(|m| {
                format!(
                    "[{},{},{},{}]",
                    m.lower,
                    m.upper,
                    m.center,
                    if m.ms1 { 1 } else { 0 }
                )
            })
            .collect::<String>()
    };
    let mut maps = swath_maps();
    annotate_swath_maps_from_file(&file, &mut maps, false, false).unwrap();
    assert_eq!(render(&maps), oracle_get(&records, "annotate", "inorder"));
    let mut inv = swath_maps();
    inv.reverse();
    let mut maps = inv.clone();
    annotate_swath_maps_from_file(&file, &mut maps, true, false).unwrap();
    assert_eq!(render(&maps), oracle_get(&records, "annotate", "sorted"));
    let mut maps = inv.clone();
    assert!(annotate_swath_maps_from_file(&file, &mut maps, false, false).is_err());
    assert_eq!(
        render(&maps),
        oracle_get(&records, "annotate", "unsorted_noforce.after")
    );
    let mut maps = inv;
    annotate_swath_maps_from_file(&file, &mut maps, false, true).unwrap();
    assert_eq!(
        render(&maps),
        oracle_get(&records, "annotate", "unsorted_force")
    );
    let mut maps = vec![SwathMap::new(0.0, 2000.0, 0.0, true)];
    maps.extend(swath_maps());
    annotate_swath_maps_from_file(&file, &mut maps, false, false).unwrap();
    assert_eq!(
        render(&maps),
        oracle_get(&records, "annotate", "ms1_skipped")
    );
}

/// Oracle: the whole `isSubsampledSpectrum_` grid (totals and sample counts
/// 0..=12, indices 0..=13), 7 of 37 and `size_t(-1)` samples.
#[test]
fn oracle_swath_qc_subsample_grid() {
    let records = oracle();
    let bits = |t: usize, s: usize, n: usize| {
        (0..n)
            .map(|i| {
                if SwathQC::is_subsampled_spectrum(t, s, i) {
                    '1'
                } else {
                    '0'
                }
            })
            .collect::<String>()
    };
    for t in 0..=12 {
        for s in 0..=12 {
            assert_eq!(
                bits(t, s, 14),
                oracle_get(&records, "subsample", &format!("{t}_{s}")),
                "{t}/{s}"
            );
        }
    }
    assert_eq!(bits(37, 7, 40), oracle_get(&records, "subsample", "37_7"));
    assert_eq!(
        bits(10, usize::MAX, 12),
        oracle_get(&records, "subsample", "10_max")
    );
}

/// Oracle: charge distributions for several sample counts and tolerances.
#[cfg(feature = "mzml")]
#[test]
fn oracle_swath_qc_charge_distribution() {
    let records = oracle();
    let maps = orbitrap_maps();
    for (n, ns) in [
        (10usize, "10"),
        (1, "1"),
        (2, "2"),
        (0, "0"),
        (usize::MAX, "-1"),
    ] {
        for (tol, ts) in [(0.04, "0.040000000000000001"), (0.01, "0.01")] {
            let cd = SwathQC::charge_distribution_of_maps(&maps, n, tol).unwrap();
            let text: String = cd.iter().map(|(c, k)| format!("{c}:{k},")).collect();
            assert_eq!(
                text,
                oracle_get(&records, "chargedist", &format!("{ns}_{ts}")),
                "{ns}/{ts}"
            );
        }
    }
    assert_eq!(
        oracle_get(&records, "qcjson", "store").replace("\\n", "\n"),
        "{\n  \"ChargeDistributionMS1\": [\n    [\n      1,\n      12\n    ],\n    [\n      6,\n      1\n    ],\n    [\n      10,\n      4\n    ]\n  ]\n}\n"
    );
}

/// Oracle: `estimateWindow` over eight residual sets, three quantiles, both
/// widths and two padding factors.
#[test]
fn oracle_estimate_window() {
    let records = oracle();
    let mut n = 0;
    for (_, key, input) in records.iter().filter(|(s, _, _)| s == "estimate_input") {
        let values = g17_list(input);
        for (q, qs) in [
            (0.99, "0.98999999999999999"),
            (0.5, "0.5"),
            (0.9, "0.90000000000000002"),
        ] {
            for full in [true, false] {
                for (pad, ps) in [(1.0, "1"), (1.3, "1.3")] {
                    let k = format!("{key}_{qs}_{}_{ps}", if full { 1 } else { 0 });
                    let expected = g17(&oracle_get(&records, "estimate", &k));
                    let got =
                        SwathMapMassCorrection::estimate_window(&values, q, full, pad).unwrap();
                    assert!(same(got, expected), "{k}: {got} vs {expected}");
                    n += 1;
                }
            }
        }
    }
    assert_eq!(n, 8 * 12);
}

/// Oracle: the default parameters, names, values, restrictions, tags and
/// descriptions, in source order.
#[test]
fn oracle_swath_map_mass_correction_defaults() {
    let records = oracle();
    let defaults = SwathMapMassCorrection::defaults().unwrap();
    let expected: Vec<&(String, String, String)> = records
        .iter()
        .filter(|(s, _, _)| s == "smc_param")
        .collect();
    let entries = &defaults.root().entries;
    assert_eq!(entries.len(), expected.len());
    for (e, (_, name, rest)) in entries.iter().zip(expected) {
        assert_eq!(&e.name, name);
        let f: Vec<&str> = rest.splitn(7, '\t').collect();
        match &e.value {
            ParamValue::Float(v) => {
                assert_eq!(f[1], "type=2");
                assert_eq!(g17(f[0]), *v, "{name}");
            }
            ParamValue::String(s) => {
                assert_eq!(f[1], "type=0");
                assert_eq!(f[0], s);
            }
            other => panic!("unexpected value {other:?}"),
        }
        let tags: String = e.tags.iter().map(|t| format!("{t};")).collect();
        assert_eq!(f[2], format!("tags={tags}"));
        let valid: String = e.valid_strings.iter().map(|t| format!("{t};")).collect();
        assert_eq!(f[3], format!("valid={valid}"));
        assert_eq!(g17(&f[4][4..]), e.min_float, "{name} min");
        assert_eq!(g17(&f[5][4..]), e.max_float, "{name} max");
        assert_eq!(f[6], format!("desc={}", e.description));
    }
    let smc = SwathMapMassCorrection::new();
    let w = oracle_get(&records, "smc_get", "windows");
    assert_eq!(
        w,
        format!(
            "{}\t{}\t{}\t{}",
            smc.fragment_mz_window(),
            smc.fragment_im_window(),
            smc.precursor_mz_window(),
            smc.precursor_im_window()
        )
    );
}

// ---------------------------------------------------------------------------
// SwathMapMassCorrection::correctMZ and correctIM (group O integration)
// ---------------------------------------------------------------------------
//
// Transcribed from SwathMapMassCorrection_test.cpp at the core pin (tier 3):
// the mock spectra, the transition groups and every expected literal.

use openms::analysis::openswath::spectrum_access::SpectrumAccessPtr;
use openms::analysis::transformations::TransformationDescription;
use openms::kernel::mrm::MRMFeature;
use std::collections::BTreeMap;

fn smmc_transition(mz: f64, prec: f64, im: f64, pep: &str, name: &str) -> LightTransition {
    LightTransition {
        product_mz: mz,
        precursor_mz: prec,
        precursor_im: im,
        peptide_ref: pep.into(),
        transition_name: name.into(),
        ..LightTransition::default()
    }
}

/// `addTransitions`: tr1..tr4 of `pep1` (drift time 11) into `group` and `exp`.
fn smmc_add_transitions(group: &mut MRMTransitionGroupType) -> LightTargetedExperiment {
    let mut exp = LightTargetedExperiment::default();
    for (i, mz) in [500.0, 600.0, 700.0, 800.0].into_iter().enumerate() {
        let name = format!("tr{}", i + 1);
        let tr = smmc_transition(mz, 412.0, 11.0, "pep1", &name);
        group.add_transition(tr.clone(), &name).unwrap();
        exp.transitions.push(tr);
    }
    exp.compounds.push(LightCompound {
        id: "pep1".into(),
        drift_time: 11.0,
        ..LightCompound::default()
    });
    exp
}

/// `addTransitionsPep2` and `addTransitionsPep3`.
fn smmc_add_pep(
    group: &mut MRMTransitionGroupType,
    exp: &mut LightTargetedExperiment,
    (mz, prec, im, pep, name): (f64, f64, f64, &str, &str),
) {
    let tr = smmc_transition(mz, prec, im, pep, name);
    group.add_transition(tr.clone(), name).unwrap();
    exp.transitions.push(tr);
    exp.compounds.push(LightCompound {
        id: pep.into(),
        drift_time: im,
        ..LightCompound::default()
    });
}

fn smmc_feature() -> MRMFeature {
    let mut feature = MRMFeature::default();
    feature.rt = 3120.0;
    feature
}

fn smmc_map(peaks: &[(f64, f32)], rt: f64, im: Option<&[f32]>) -> SpectrumAccessPtr {
    let mut exp = MSExperiment::default();
    exp.spectra.push(spectrum(rt, peaks, im));
    Arc::new(access(exp))
}

fn smmc_swath(sptr: SpectrumAccessPtr, lower: f64, upper: f64, im: (f64, f64)) -> SwathMap {
    SwathMap {
        sptr: Some(sptr),
        lower,
        upper,
        center: 412.5,
        im_lower: im.0,
        im_upper: im.1,
        ms1: false,
        ..SwathMap::default()
    }
}

fn smmc_first_mz(map: &SwathMap) -> Vec<f64> {
    map.sptr
        .as_ref()
        .unwrap()
        .spectrum_by_id(0)
        .unwrap()
        .mz_array()
        .data
        .clone()
}

fn smmc_groups<'a>(
    groups: &[&'a MRMTransitionGroupType],
) -> BTreeMap<String, &'a MRMTransitionGroupType> {
    groups
        .iter()
        .enumerate()
        .map(|(i, g)| (format!("group{}", i + 1), *g))
        .collect()
}

/// SwathMapMassCorrection_test.cpp: `correctMZ`, every section.
#[test]
fn swath_map_mass_correction_correct_mz_class_test() {
    // targets for correction are: 500.00, 600.00, 700.00, 800.00
    // "measured data" as input:   500.02, 600.00, 699.97, 800.02
    let mut group = MRMTransitionGroupType::default();
    group.add_feature(smmc_feature()).unwrap();
    let targ_exp = smmc_add_transitions(&mut group);
    let map_of_groups = smmc_groups(&[&group, &group, &group]);

    let sptr = smmc_map(
        &[
            (500.02, 150.0),
            (600.00, 150.0),
            (699.97, 22500.01), // double the weight of all other data
            (800.02, 150.0),
        ],
        3121.0,
        None,
    );
    // never used: the second PASEF map
    let sptr2 = smmc_map(
        &[
            (499.97, 1_500_000.0),
            (600.00, 150.0),
            (699.97, 150.0),
            (800.02, 150.0),
        ],
        3122.0,
        None,
    );
    let map = smmc_swath(sptr, 400.0, 425.0, (10.0, 20.0));
    // the class test never sets `mapPASEF.upper`, which stays 0: the window can
    // never contain a precursor
    let map_pasef = smmc_swath(sptr2, 400.0, 0.0, (20.0, 30.0));

    let pasef = false;
    let run = |function: &str, window: f64, maps: &[SwathMap], pasef: bool| -> Vec<SwathMap> {
        let mut mc = SwathMapMassCorrection::new();
        mc.mz_correction_function = function.into();
        mc.mz_extraction_window = window;
        let mut maps = maps.to_vec();
        mc.correct_mz(&map_of_groups, &targ_exp, &mut maps, pasef)
            .unwrap();
        maps
    };

    // should work with empty maps, with the default (none) and a regression
    let mut mc = SwathMapMassCorrection::new();
    let mut empty: Vec<SwathMap> = Vec::new();
    mc.correct_mz(&map_of_groups, &targ_exp, &mut empty, pasef)
        .unwrap();
    mc.mz_correction_function = "unweighted_regression".into();
    mc.correct_mz(&map_of_groups, &targ_exp, &mut empty, pasef)
        .unwrap();
    // nothing was found: the window estimates are 0 (empty residuals)
    assert_eq!(mc.fragment_mz_window(), 0.0);

    // none: the spectrum is untouched
    let maps = run("none", 1.0, std::slice::from_ref(&map), pasef);
    let data = smmc_first_mz(&maps[0]);
    for (got, want) in data.iter().zip([500.02, 600.00, 699.97, 800.02]) {
        real_similar(*got, want);
    }

    // unweighted regression, window 0.05
    let linear = |a: f64, b: f64| [500.02, 600.00, 699.97, 800.02].map(|x| a + b * x);
    let maps = run(
        "unweighted_regression",
        0.05,
        std::slice::from_ref(&map),
        pasef,
    );
    for (got, want) in smmc_first_mz(&maps[0])
        .iter()
        .zip(linear(-0.00428216, 0.999986))
    {
        real_similar(*got, want);
    }
    // window 1.0
    let maps = run(
        "unweighted_regression",
        1.0,
        std::slice::from_ref(&map),
        pasef,
    );
    for (got, want) in smmc_first_mz(&maps[0])
        .iter()
        .zip(linear(-0.0219795, 1.00003))
    {
        real_similar(*got, want);
    }
    // the same call again: the source test lists a third set of literals for
    // the identical parameters (-0.0315101, 1.00005); it is the same call and
    // the tolerance of TEST_REAL_SIMILAR admits both
    for (got, want) in smmc_first_mz(&maps[0])
        .iter()
        .zip(linear(-0.0315101, 1.00005))
    {
        real_similar(*got, want);
    }

    // quadratic regression
    let quad =
        |a: f64, b: f64, c: f64| [500.02, 600.00, 699.97, 800.02].map(|x| a + b * x + c * x * x);
    let maps = run(
        "quadratic_regression",
        1.0,
        std::slice::from_ref(&map),
        pasef,
    );
    for (got, want) in smmc_first_mz(&maps[0]).iter().zip(quad(
        -0.7395987927448004,
        1.002305255194642,
        -1.750157412772069e-06,
    )) {
        real_similar(*got, want);
    }
    // weighted quadratic regression
    let maps = run(
        "weighted_quadratic_regression",
        1.0,
        std::slice::from_ref(&map),
        pasef,
    );
    for (got, want) in smmc_first_mz(&maps[0]).iter().zip(quad(
        -0.8323316718451679,
        1.002596944948891,
        -1.967834556637627e-06,
    )) {
        real_similar(*got, want);
    }
    // the delta-ppm functions
    let maps = run(
        "quadratic_regression_delta_ppm",
        1.0,
        std::slice::from_ref(&map),
        pasef,
    );
    for (got, want) in smmc_first_mz(&maps[0]).iter().zip([
        499.997160932778,
        600.010219722383,
        699.988081672119,
        800.004537672719,
    ]) {
        real_similar(*got, want);
    }
    let maps = run(
        "weighted_quadratic_regression_delta_ppm",
        1.0,
        std::slice::from_ref(&map),
        pasef,
    );
    for (got, want) in smmc_first_mz(&maps[0]).iter().zip([
        499.996336995751,
        600.013185628794,
        699.992311403648,
        800.005854568825,
    ]) {
        real_similar(*got, want);
    }

    // PASEF on: only the window containing the IM of the transitions is used,
    // so the result is that of the first map alone
    let maps = run(
        "unweighted_regression",
        0.05,
        &[map.clone(), map_pasef],
        true,
    );
    for (got, want) in smmc_first_mz(&maps[0])
        .iter()
        .zip(linear(-0.00428216, 0.999986))
    {
        real_similar(*got, want);
    }
}

/// SwathMapMassCorrection_test.cpp: `correctIM`, every section.
#[test]
fn swath_map_mass_correction_correct_im_class_test() {
    // m/z targets:   500, 600, 700, 800, 900, 950
    // mobility:       11,  11,  11,  11,  15,  20
    // measured input: 22, 21.5, 20.5, 21, 24, 31
    let mut gr1 = MRMTransitionGroupType::default();
    let mut gr2 = MRMTransitionGroupType::default();
    let mut gr3 = MRMTransitionGroupType::default();
    for g in [&mut gr1, &mut gr2, &mut gr3] {
        g.add_feature(smmc_feature()).unwrap();
    }
    let mut targ_exp = smmc_add_transitions(&mut gr1);
    smmc_add_pep(&mut gr2, &mut targ_exp, (900.0, 500.0, 15.0, "pep2", "tr5"));
    smmc_add_pep(&mut gr3, &mut targ_exp, (950.0, 600.0, 20.0, "pep3", "tr6"));
    let groups = smmc_groups(&[&gr1, &gr2, &gr3]);

    let peaks = [
        (500.02, 150.0_f32),
        (600.00, 150.0),
        (699.97, 22500.01),
        (800.02, 150.0),
        (900.02, 150.0),
        (950.02, 150.0),
    ];
    let im = [22.0_f32, 21.5, 20.5, 21.0, 24.0, 31.0];
    // the class test adds the second mock spectrum to `exp` (not `exp2`), so
    // the first map holds two spectra at RT 3121
    let im_second = [22.3_f32, 21.5, 20.5, 21.5, 24.4, 31.0];
    let mut exp = MSExperiment::default();
    exp.spectra.push(spectrum(3121.0, &peaks, Some(&im)));
    exp.spectra.push(spectrum(3121.0, &peaks, Some(&im_second)));
    let sptr: SpectrumAccessPtr = Arc::new(access(exp));
    let sptr2: SpectrumAccessPtr = Arc::new(access(MSExperiment::default()));
    let mut exp_ms1 = MSExperiment::default();
    exp_ms1.spectra.push(spectrum(
        3121.0,
        &[(412.02, 150.0), (500.02, 150.0), (600.01, 150.0)],
        Some(&[22.0, 24.0, 31.0]),
    ));
    let sptr_ms1: SpectrumAccessPtr = Arc::new(access(exp_ms1));

    let map = smmc_swath(sptr, 400.0, 800.0, (10.0, 30.0));
    let map_pasef = smmc_swath(sptr2, 400.0, 800.0, (100.0, 200.0));
    let ms1_map = SwathMap {
        sptr: Some(sptr_ms1),
        ms1: true,
        ..SwathMap::default()
    };

    let correct = |configure: &dyn Fn(&mut SwathMapMassCorrection),
                   maps: &[SwathMap],
                   pasef: bool|
     -> Result<TransformationDescription, Error> {
        let mut mc = SwathMapMassCorrection::new();
        configure(&mut mc);
        let mut trafo = TransformationDescription::default();
        mc.correct_im(&groups, &targ_exp, maps, pasef, &mut trafo)?;
        Ok(trafo)
    };
    let ms2_config = |mc: &mut SwathMapMassCorrection| {
        mc.mz_correction_function = "none".into();
        mc.mz_extraction_window = 1.0;
        mc.im_extraction_window = 100.0;
    };

    // empty maps: the default im_extraction_window of -1 returns at once
    let trafo = correct(&|_| {}, &[], false).unwrap();
    real_similar(trafo.apply(10.0).unwrap(), 10.0);
    real_similar(trafo.apply(100.0).unwrap(), 100.0);
    let trafo = correct(
        &|mc| mc.mz_correction_function = "unweighted_regression".into(),
        &[],
        false,
    )
    .unwrap();
    real_similar(trafo.apply(10.0).unwrap(), 10.0);
    real_similar(trafo.apply(100.0).unwrap(), 100.0);

    // MS2-based alignment, with and without the MS1 map
    for maps in [vec![map.clone(), ms1_map.clone()], vec![map.clone()]] {
        let trafo = correct(&ms2_config, &maps, false).unwrap();
        real_similar(trafo.apply(10.0).unwrap(), 0.889721627408994);
        real_similar(trafo.apply(20.0).unwrap(), 10.0974304068522);
        real_similar(trafo.apply(30.0).unwrap(), 19.3051391862955);
    }

    // a single peptide: the regression is only an intercept
    let mut map_single = map.clone();
    map_single.upper = 425.0;
    let trafo = correct(&ms2_config, &[map_single], false).unwrap();
    for x in [10.0, 20.0, 30.0] {
        real_similar(trafo.apply(x).unwrap(), 11.0);
    }

    // PASEF data: the second window is excluded
    let trafo = correct(
        &ms2_config,
        &[map.clone(), map_pasef, ms1_map.clone()],
        true,
    )
    .unwrap();
    real_similar(trafo.apply(10.0).unwrap(), 0.889721627408994);
    real_similar(trafo.apply(20.0).unwrap(), 10.0974304068522);
    real_similar(trafo.apply(30.0).unwrap(), 19.3051391862955);

    // MS1 calibration with an MS1 map but no MS2 map: no point, UnableToFit
    let ms1_config = |mc: &mut SwathMapMassCorrection| {
        ms2_config(mc);
        mc.ms1_im_calibration = true;
    };
    match correct(&ms1_config, &[ms1_map.clone()], false) {
        Err(Error::InvalidValue(message)) => assert!(message.starts_with("UnableToFit")),
        other => panic!("expected UnableToFit, got {other:?}"),
    }

    // MS1 ion mobility alignment
    let trafo = correct(&ms1_config, &[map, ms1_map], false).unwrap();
    real_similar(trafo.apply(10.0).unwrap(), 0.835820895522389);
    real_similar(trafo.apply(20.0).unwrap(), 10.089552238806);
    real_similar(trafo.apply(30.0).unwrap(), 19.3432835820896);
}

// ---------------------------------------------------------------------------
// return_chromatogram over the heavyweight TargetedExperiment (group O
// integration)
// ---------------------------------------------------------------------------

/// `ChromatogramExtractor_input.TraML` through the ported `TraMLFile` and
/// `convertTargetedExp` is the experiment this file transcribed before either
/// existed, and `return_chromatogram<TargetedExperiment>` annotates exactly
/// what `return_chromatogram<LightTargetedExperiment>` does for it, MS2 and MS1.
#[cfg(all(feature = "mzml", feature = "traml"))]
#[test]
fn return_chromatogram_for_the_targeted_experiment_matches_the_light_instantiation() {
    use openms::analysis::openswath::data_access_helper::convert_targeted_exp;
    let heavy = openms::format::traml::load(data(
        "openswath_extraction_ChromatogramExtractor_input.TraML",
    ))
    .unwrap();
    let converted = convert_targeted_exp(&heavy).unwrap();
    let transcribed = traml_light_experiment();
    assert_eq!(converted.compounds().len(), transcribed.compounds().len());
    for (got, want) in converted.compounds().iter().zip(transcribed.compounds()) {
        assert_eq!(got.id, want.id);
        assert_eq!(got.sequence, want.sequence);
        assert_eq!(got.rt, want.rt);
    }
    assert_eq!(
        converted.transitions().len(),
        transcribed.transitions().len()
    );
    for (got, want) in converted
        .transitions()
        .iter()
        .zip(transcribed.transitions())
    {
        assert_eq!(got.transition_name, want.transition_name);
        assert_eq!(got.peptide_ref, want.peptide_ref);
        assert_eq!(got.precursor_mz, want.precursor_mz);
        assert_eq!(got.product_mz, want.product_mz);
        assert_eq!(got.library_intensity, want.library_intensity);
    }

    let bytes = std::fs::read(data(
        "openswath_extraction_ChromatogramExtractor_input.mzML",
    ))
    .unwrap();
    let exp = Arc::new(openms::format::mzml::read(bytes.as_slice()).unwrap());
    let input = SpectrumAccessOpenMS::new(Arc::clone(&exp));
    for (ms1, isotopes, im_width) in [(false, 0, 0.0), (true, 1, 0.2)] {
        let mut output = Vec::new();
        let mut coords = Vec::new();
        prepare_coordinates(&mut output, &mut coords, &converted, -1.0, ms1, isotopes).unwrap();
        if ms1 {
            coords[0].ion_mobility = 0.9;
        } else {
            ChromatogramExtractor::new()
                .extract_chromatograms(&input, &mut output, &coords, 0.05, false, "tophat")
                .unwrap();
        }
        let mut from_light = Vec::new();
        return_chromatogram(
            &output,
            &coords,
            &converted,
            &exp.spectra[0],
            &mut from_light,
            ms1,
            im_width,
        )
        .unwrap();
        let mut from_heavy = Vec::new();
        return_chromatogram(
            &output,
            &coords,
            &heavy,
            &exp.spectra[0],
            &mut from_heavy,
            ms1,
            im_width,
        )
        .unwrap();
        assert!(!from_heavy.is_empty());
        assert_eq!(from_heavy, from_light, "ms1 = {ms1}");
    }

    // Derived from `extract_id_<TargetedExperiment>`: an id that is neither a
    // peptide nor a compound annotates the empty sequence and charge 0, where
    // the light instantiation dereferences a missing compound.
    let mut coord = coord(500.0, 0.0, 1.0, "tr1_Precursor_i0", -1.0);
    coord.mz_precursor = 500.0;
    let mut stray = heavy.clone();
    stray.set_peptides(Vec::new());
    let mut out = Vec::new();
    return_chromatogram(
        &[Chromatogram::new()],
        std::slice::from_ref(&coord),
        &stray,
        &exp.spectra[0],
        &mut out,
        true,
        0.0,
    )
    .unwrap();
    assert_eq!(
        out[0].precursor.cv_terms.metadata["peptide_sequence"],
        MetaValue::from("")
    );
    assert_eq!(out[0].precursor.charge, 0);
}
