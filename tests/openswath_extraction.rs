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
use openms::analysis::openswath::chromatogram_extractor_algorithm::{
    ChromatogramExtractorAlgorithm, ExtractionCoordinates, ExtractionFilter,
};
use openms::analysis::openswath::data_structures::Chromatogram;
use openms::analysis::openswath::peak_map_extractor::{ExtractedPeakMap, PeakMapExtractor};
use openms::analysis::openswath::spectrum_access::SpectrumAccess;
use openms::analysis::openswath::spectrum_access_openms::SpectrumAccessOpenMS;
use openms::analysis::openswath::swath_map::SwathMap;
use openms::analysis::openswath::swath_window_loader::{
    annotate_swath_maps, annotate_swath_maps_from_file, parse_swath_windows, read_swath_windows,
};
use openms::kernel::{DataArray, MSExperiment, MSSpectrum, Peak1D};
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
    let mut s = MSSpectrum::default();
    s.rt = rt;
    s.ms_level = 2;
    s.peaks = peaks.iter().map(|&(mz, i)| Peak1D::new(mz, i)).collect();
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
    let mut f = |mz: f64, w: f64, ppm: bool, pos: &mut usize| {
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
    let mut f = |mz: f64, im: f64, pos: &mut usize| {
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
