// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Scan mobility, ion-mobility bulk exports and the RT/m/z raster of
//! `KERNEL/MSExperiment.h`, plus the mobility half of `KERNEL/AreaIterator.h`.
//!
//! Every literal named "source" below is transcribed from
//! `MSExperiment_test.cpp` or `AreaIterator_test.cpp` at Core SDK `bc9cc12`
//! (tier 3, source review). The remaining tests are native invariants,
//! resource bounds and atomicity (tier 4).

use openms::Error;
use openms::kernel::experiment_mobility::{
    ExperimentMobilityLimits, FlatPeakDataIm, RtMzRaster, SpectrumPeakDataIm,
};
use openms::kernel::ranges::{MSDim, RangeBase, RangeManager};
use openms::kernel::spectrum_mobility::RasterAggregation;
use openms::kernel::{AreaBounds, AreaOptions, DataArray, MSExperiment, MSSpectrum, Peak1D};

fn scan(rt: f64, drift_time: f64, level: u32, mzs: &[f64]) -> MSSpectrum {
    MSSpectrum {
        rt,
        drift_time,
        ms_level: level,
        peaks: mzs.iter().map(|&mz| Peak1D::new(mz, mz as f32)).collect(),
        ..Default::default()
    }
}

/// AreaIterator_test.cpp's five-scan fixture, including its drift times and the
/// empty third scan.
fn area_fixture() -> MSExperiment {
    MSExperiment {
        spectra: vec![
            scan(2., 1., 1, &[502., 510.]),
            scan(4., 1.4, 1, &[504., 506.]),
            scan(6., 1.6, 1, &[]),
            scan(8., 1.8, 1, &[504.1, 506.1]),
            scan(10., 1.99, 1, &[502.1, 510.1]),
        ],
        ..Default::default()
    }
}

/// MSExperiment_test.cpp's `exp_area`: `set2DData` of six RT/m/z points,
/// grouped into three scans, with drift times 0, 1 and 2 assigned afterwards.
fn experiment_area_fixture() -> MSExperiment {
    MSExperiment {
        spectra: vec![
            scan(-1., 0., 1, &[2.]),
            scan(1., 1., 1, &[2., 3.]),
            scan(2., 2., 1, &[10., 11., 12.]),
        ],
        ..Default::default()
    }
}

/// A spectrum whose peaks each carry an ion mobility value, as the source class
/// test builds them: a `UserParam` "Ion Mobility" float data array.
fn im_scan(rt: f64, level: u32, peaks: &[(f64, f32)], mobilities: &[f32]) -> MSSpectrum {
    MSSpectrum {
        rt,
        ms_level: level,
        peaks: peaks
            .iter()
            .map(|&(mz, intensity)| Peak1D::new(mz, intensity))
            .collect(),
        float_data_arrays: vec![DataArray::new("Ion Mobility", mobilities.to_vec())],
        ..Default::default()
    }
}

/// MSExperiment_test.cpp's `get2DPeakDataIMPerSpectrum` fixture.
fn im_export_fixture() -> MSExperiment {
    MSExperiment {
        spectra: vec![
            im_scan(1., 1, &[(100., 1000.), (200., 2000.)], &[0.8, 1.2]),
            im_scan(2., 1, &[(150., 1500.), (250., 2500.)], &[1.5, 2.0]),
            im_scan(3., 2, &[(175., 1750.)], &[1.8]),
        ],
        ..Default::default()
    }
}

fn bounds(rt: (f64, f64), mz: (f64, f64)) -> AreaBounds {
    AreaBounds::new(rt.0, rt.1, mz.0, mz.1).unwrap()
}

// ---------------------------------------------------------------------------
// AreaIterator.h: getDriftTime, lowIM/highIM
// ---------------------------------------------------------------------------

#[test]
fn source_area_iterator_drift_time_and_scan_mobility_window() {
    let exp = area_fixture();
    // START_SECTION((CoordinateType getDriftTime() const)): RTBegin(3),
    // RTEnd(9), lowMZ(503), highMZ(509), lowIM(1), highIM(1.5). Only scan 1
    // (drift time 1.4) is inside the mobility window.
    let options = AreaOptions::new(bounds((3., 9.), (503., 509.)), 1)
        .with_mobility(1., 1.5)
        .unwrap();
    let points: Vec<_> = exp.area_iter(options).unwrap().collect();
    assert_eq!(
        points.iter().map(|p| p.peak.mz).collect::<Vec<_>>(),
        [504., 506.]
    );
    assert_eq!(
        points.iter().map(|p| p.spectrum.rt).collect::<Vec<_>>(),
        [4., 4.]
    );
    assert_eq!(
        points.iter().map(|p| p.drift_time()).collect::<Vec<_>>(),
        [1.4, 1.4]
    );

    // [EXTRA] Overall test: a mobility window covering every scan selects the
    // same eight peaks as the unrestricted traversals do.
    let all = [502., 510., 504., 506., 504.1, 506.1, 502.1, 510.1];
    for options in [
        AreaOptions::new(bounds((0., 15.), (500., 520.)), 1),
        AreaOptions::new(AreaBounds::default(), 1),
        AreaOptions::new(AreaBounds::default(), 1)
            .with_mobility(0., 2.)
            .unwrap(),
        AreaOptions::new(bounds((0., 15.), (500., 520.)), 1)
            .with_mobility(0., 2.)
            .unwrap(),
    ] {
        assert_eq!(
            exp.area_iter(options)
                .unwrap()
                .map(|p| p.peak.mz)
                .collect::<Vec<_>>(),
            all
        );
    }
    // [EXTRA] Overall test: "Test with empty IM range" lowIM(0) highIM(0.9).
    assert_eq!(
        exp.area_iter(
            AreaOptions::new(AreaBounds::default(), 1)
                .with_mobility(0., 0.9)
                .unwrap()
        )
        .unwrap()
        .count(),
        0
    );
}

#[test]
fn mutable_area_iteration_snapshots_the_scan_drift_time() {
    let mut exp = area_fixture();
    let options = AreaOptions::new(bounds((3., 9.), (503., 509.)), 1)
        .with_mobility(1., 1.5)
        .unwrap();
    let observed: Vec<_> = exp
        .area_iter_mut(options)
        .unwrap()
        .map(|p| (p.peak_index, p.rt, p.drift_time))
        .collect();
    assert_eq!(observed, [(0, 4., 1.4), (1, 4., 1.4)]);
}

// ---------------------------------------------------------------------------
// MSExperiment.h: areaBegin / areaBeginConst from a RangeManager
// ---------------------------------------------------------------------------

fn source_range_manager(mobility: Option<(f64, f64)>) -> RangeManager {
    let mut manager = RangeManager::experiment();
    manager
        .extend_range(MSDim::Rt, &RangeBase::from_min_max(0., 2.).unwrap())
        .unwrap();
    manager
        .extend_range(MSDim::Mz, &RangeBase::from_min_max(3., 11.).unwrap())
        .unwrap();
    if let Some((min, max)) = mobility {
        manager
            .extend_range(MSDim::Mobility, &RangeBase::from_min_max(min, max).unwrap())
            .unwrap();
    }
    manager
}

#[test]
fn source_range_manager_area_overloads() {
    let mut exp = experiment_area_fixture();
    // START_SECTION(areaBeginConst(const RangeManagerType& range)): RangeRT
    // (0,2) and RangeMZ (3,11), intensity and mobility left empty.
    let manager = source_range_manager(None);
    assert_eq!(
        exp.area_begin_from_ranges(&manager, 1)
            .unwrap()
            .map(|p| p.peak.mz)
            .collect::<Vec<_>>(),
        [3., 10., 11.]
    );
    // "add mobility as constraint": RangeMobility (2,2) keeps only the scan
    // whose drift time is 2.
    let restricted = source_range_manager(Some((2., 2.)));
    assert_eq!(
        exp.area_begin_from_ranges(&restricted, 1)
            .unwrap()
            .map(|p| p.peak.mz)
            .collect::<Vec<_>>(),
        [10., 11.]
    );
    // START_SECTION(areaBegin(const RangeManagerType& range)): the mutable
    // overload selects the same peaks.
    assert_eq!(
        exp.area_begin_mut_from_ranges(&manager, 1)
            .unwrap()
            .map(|p| p.peak.mz)
            .collect::<Vec<_>>(),
        [3., 10., 11.]
    );
    assert_eq!(
        exp.area_begin_mut_from_ranges(&restricted, 1)
            .unwrap()
            .map(|p| p.peak.mz)
            .collect::<Vec<_>>(),
        [10., 11.]
    );
}

#[test]
fn range_manager_area_options_treat_empty_and_missing_dimensions_as_unrestricted() {
    let exp = experiment_area_fixture();
    // An all-empty experiment manager restricts nothing, so every MS1 peak of
    // the fixture is visited, including the scan at RT -1.
    let empty = RangeManager::experiment();
    assert_eq!(exp.area_begin_from_ranges(&empty, 1).unwrap().count(), 6);
    // A spectrum manager carries only m/z and intensity; the absent RT and
    // mobility dimensions behave exactly like empty ones.
    let mut spectrum = RangeManager::spectrum();
    spectrum
        .extend_range(MSDim::Mz, &RangeBase::from_min_max(3., 11.).unwrap())
        .unwrap();
    assert!(!spectrum.has_dim(MSDim::Rt));
    assert_eq!(
        exp.area_begin_from_ranges(&spectrum, 1)
            .unwrap()
            .map(|p| p.peak.mz)
            .collect::<Vec<_>>(),
        [3., 10., 11.]
    );
    // The intensity dimension is ignored: the source area iterator has no
    // intensity filter, and the fixture's intensities are all outside (0, 1).
    let mut with_intensity = RangeManager::experiment();
    with_intensity
        .extend_range(MSDim::Intensity, &RangeBase::from_min_max(0., 1.).unwrap())
        .unwrap();
    assert_eq!(
        exp.area_begin_from_ranges(&with_intensity, 1)
            .unwrap()
            .count(),
        6
    );
    // Source byte narrowing applies to the level, exactly as for the scalar
    // wrapper: 256 selects level zero, which no fixture scan has.
    assert_eq!(exp.area_begin_from_ranges(&empty, 256).unwrap().count(), 0);
    let options = AreaOptions::from_range_manager(&empty, 256).unwrap();
    assert_eq!(options.ms_level, 256);
    assert_eq!(options.mobility, None);
}

// ---------------------------------------------------------------------------
// MSExperiment.h: IMBegin / IMEnd
// ---------------------------------------------------------------------------

fn drift_time_fixture() -> MSExperiment {
    MSExperiment {
        spectra: [30., 40., 45., 50.]
            .iter()
            .map(|&dt| scan(-1., dt, 1, &[]))
            .collect(),
        ..Default::default()
    }
}

#[test]
fn source_im_begin_and_im_end() {
    let exp = drift_time_fixture();
    // START_SECTION((Iterator IMBegin(CoordinateType im)))
    assert_eq!(exp.im_begin(20.).unwrap(), 0);
    assert_eq!(exp.spectra[exp.im_begin(20.).unwrap()].drift_time, 30.);
    assert_eq!(exp.spectra[exp.im_begin(30.).unwrap()].drift_time, 30.);
    assert_eq!(exp.spectra[exp.im_begin(31.).unwrap()].drift_time, 40.);
    assert_eq!(exp.im_begin(55.).unwrap(), exp.spectra.len());
    // START_SECTION((Iterator IMEnd(CoordinateType rt)))
    assert_eq!(exp.spectra[exp.im_end(20.).unwrap()].drift_time, 30.);
    assert_eq!(exp.spectra[exp.im_end(30.).unwrap()].drift_time, 40.);
    assert_eq!(exp.spectra[exp.im_end(31.).unwrap()].drift_time, 40.);
    assert_eq!(exp.im_end(55.).unwrap(), exp.spectra.len());
}

#[test]
fn scan_mobility_search_checks_order_and_finiteness() {
    let exp = MSExperiment::default();
    assert_eq!(exp.im_begin(1.).unwrap(), 0);
    assert_eq!(exp.im_end(1.).unwrap(), 0);
    let mut unsorted = drift_time_fixture();
    unsorted.spectra[1].drift_time = 90.;
    assert!(matches!(unsorted.im_begin(40.), Err(Error::UnsortedData)));
    assert!(matches!(unsorted.im_end(40.), Err(Error::UnsortedData)));
    let mut nonfinite = drift_time_fixture();
    nonfinite.spectra[2].drift_time = f64::NAN;
    assert!(nonfinite.im_begin(40.).is_err());
    assert!(matches!(
        drift_time_fixture().im_begin(f64::NAN),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        drift_time_fixture().im_end(f64::INFINITY),
        Err(Error::InvalidValue(_))
    ));
    // Unset drift times all equal the source sentinel, so they are sorted and
    // every query lands at one end of the run.
    let unset = MSExperiment {
        spectra: vec![MSSpectrum::default(), MSSpectrum::default()],
        ..Default::default()
    };
    assert_eq!(unset.im_begin(-1.).unwrap(), 0);
    assert_eq!(unset.im_end(-1.).unwrap(), 2);
}

// ---------------------------------------------------------------------------
// MSExperiment.h: isIMFrame
// ---------------------------------------------------------------------------

#[test]
fn source_is_im_frame() {
    // START_SECTION(bool MSExperiment::isIMFrame() const): one RT for all four
    // scans, drift times 30, 40, 45, 50.
    let mut exp = MSExperiment {
        spectra: [30., 40., 45., 50.]
            .iter()
            .map(|&dt| scan(1., dt, 1, &[]))
            .collect(),
        ..Default::default()
    };
    assert!(exp.is_im_frame().unwrap());
    exp.spectra[3].rt = 2.; // changing RT ...
    assert!(!exp.is_im_frame().unwrap());
    exp.spectra[3].rt = 1.; // undo
    exp.spectra[2].drift_time = exp.spectra[1].drift_time; // duplicate drift time
    assert!(!exp.is_im_frame().unwrap());
    exp.spectra[3].rt = 2.;
    assert!(!exp.is_im_frame().unwrap());
}

#[test]
fn im_frame_predicate_edge_cases_follow_the_source() {
    // An empty run is not a frame.
    assert!(!MSExperiment::default().is_im_frame().unwrap());
    // A single scan is reported as a frame even with no ion mobility at all,
    // because the source only compares against the previous drift time and
    // there is none. Recorded as a source defect candidate.
    let lone = MSExperiment {
        spectra: vec![scan(1., -1., 1, &[])],
        ..Default::default()
    };
    assert!(lone.is_im_frame().unwrap());
    // Two scans without drift times are not a frame: both hold the sentinel.
    let pair = MSExperiment {
        spectra: vec![scan(1., -1., 1, &[]), scan(1., -1., 1, &[])],
        ..Default::default()
    };
    assert!(!pair.is_im_frame().unwrap());
    // Only the immediately preceding scan is compared, so a repeated but
    // non-adjacent drift time still reports a frame.
    let revisited = MSExperiment {
        spectra: vec![
            scan(1., 1., 1, &[]),
            scan(1., 2., 1, &[]),
            scan(1., 1., 1, &[]),
        ],
        ..Default::default()
    };
    assert!(revisited.is_im_frame().unwrap());
    // Non-finite coordinates are refused rather than silently answered.
    let bad_rt = MSExperiment {
        spectra: vec![scan(f64::NAN, 1., 1, &[])],
        ..Default::default()
    };
    assert!(bad_rt.is_im_frame().is_err());
    let bad_drift = MSExperiment {
        spectra: vec![scan(1., f64::INFINITY, 1, &[])],
        ..Default::default()
    };
    assert!(bad_drift.is_im_frame().is_err());
}

// ---------------------------------------------------------------------------
// MSExperiment.h: get2DPeakDataIMPerSpectrum
// ---------------------------------------------------------------------------

#[test]
fn source_get_2d_peak_data_im_per_spectrum() {
    let exp = im_export_fixture();
    // Test 1: full range, MS level 1.
    let out = exp
        .get_2d_peak_data_im_per_spectrum(bounds((0., 4.), (0., 300.)), 1)
        .unwrap();
    assert_eq!(out.rt.len(), 2);
    assert_eq!(out.mz.len(), 2);
    assert_eq!(out.intensity.len(), 2);
    assert_eq!(out.ion_mobility.len(), 2);
    assert_eq!(out.rt, [1.0, 2.0]);
    assert_eq!(out.mz[0], [100.0, 200.0]);
    assert_eq!(out.intensity[0], [1000.0, 2000.0]);
    assert_eq!(out.ion_mobility[0], [0.8, 1.2]);
    assert_eq!(out.mz[1], [150.0, 250.0]);
    assert_eq!(out.intensity[1], [1500.0, 2500.0]);
    assert_eq!(out.ion_mobility[1], [1.5, 2.0]);

    // Test 2: limited RT range.
    let out = exp
        .get_2d_peak_data_im_per_spectrum(bounds((1.5, 2.5), (0., 300.)), 1)
        .unwrap();
    assert_eq!(out.rt, [2.0]);
    assert_eq!(out.ion_mobility, [vec![1.5, 2.0]]);

    // Test 3: limited m/z range. The second peak of the first spectrum and the
    // second of the second are outside, so only one row with one peak remains.
    let out = exp
        .get_2d_peak_data_im_per_spectrum(bounds((0., 4.), (120., 180.)), 1)
        .unwrap();
    assert_eq!(out.rt, [2.0]);
    assert_eq!(out.mz, [vec![150.0]]);
    assert_eq!(out.ion_mobility, [vec![1.5]]);

    // Test 4: MS level 2.
    let out = exp
        .get_2d_peak_data_im_per_spectrum(bounds((0., 4.), (0., 300.)), 2)
        .unwrap();
    assert_eq!(out.rt, [3.0]);
    assert_eq!(out.mz, [vec![175.0]]);
    assert_eq!(out.intensity, [vec![1750.0]]);
    assert_eq!(out.ion_mobility, [vec![1.8]]);

    // Test 5: empty range.
    let out = exp
        .get_2d_peak_data_im_per_spectrum(bounds((5., 6.), (0., 300.)), 1)
        .unwrap();
    assert_eq!(out, SpectrumPeakDataIm::default());

    // Test 6: a spectrum without ion mobility data yields the -1 sentinel.
    let plain = MSExperiment {
        spectra: vec![MSSpectrum {
            rt: 1.,
            ms_level: 1,
            peaks: vec![Peak1D::new(100., 1000.), Peak1D::new(200., 2000.)],
            ..Default::default()
        }],
        ..Default::default()
    };
    let out = plain
        .get_2d_peak_data_im_per_spectrum(bounds((0., 4.), (0., 300.)), 1)
        .unwrap();
    assert_eq!(out.rt.len(), 1);
    assert_eq!(out.ion_mobility, [vec![-1.0, -1.0]]);
}

// ---------------------------------------------------------------------------
// MSExperiment.h: get2DPeakDataIM
// ---------------------------------------------------------------------------

#[test]
fn source_get_2d_peak_data_im() {
    let exp = MSExperiment {
        spectra: im_export_fixture().spectra[..2].to_vec(),
        ..Default::default()
    };
    // Test 1: full range, MS level 1.
    let out = exp
        .get_2d_peak_data_im(bounds((0., 4.), (0., 300.)), 1)
        .unwrap();
    assert_eq!(out.rt, [1.0, 1.0, 2.0, 2.0]);
    assert_eq!(out.mz, [100.0, 200.0, 150.0, 250.0]);
    assert_eq!(out.intensity, [1000.0, 2000.0, 1500.0, 2500.0]);
    assert_eq!(out.ion_mobility, [0.8, 1.2, 1.5, 2.0]);

    // Test 2: limited RT range.
    let out = exp
        .get_2d_peak_data_im(bounds((1.5, 2.5), (0., 300.)), 1)
        .unwrap();
    assert_eq!(out.rt, [2.0, 2.0]);
    assert_eq!(out.mz, [150.0, 250.0]);
    assert_eq!(out.intensity, [1500.0, 2500.0]);
    assert_eq!(out.ion_mobility, [1.5, 2.0]);

    // Test 3: a spectrum without ion mobility data yields the -1 sentinel.
    let plain = MSExperiment {
        spectra: vec![MSSpectrum {
            rt: 1.,
            ms_level: 1,
            peaks: vec![Peak1D::new(100., 1000.), Peak1D::new(200., 2000.)],
            ..Default::default()
        }],
        ..Default::default()
    };
    let out = plain
        .get_2d_peak_data_im(bounds((0., 4.), (0., 300.)), 1)
        .unwrap();
    assert_eq!(out.rt.len(), 2);
    assert_eq!(out.ion_mobility, [-1.0, -1.0]);
}

#[test]
fn flat_mobility_export_reads_each_peak_own_spectrum_at_rt_minus_one() {
    // The source declares its row cursor inside the per-peak loop, so a
    // spectrum at RT exactly -1 reports -1 mobility even though it carries an
    // ion mobility array. This port reads the array.
    let exp = MSExperiment {
        spectra: vec![im_scan(-1., 1, &[(100., 10.), (200., 20.)], &[0.4, 0.5])],
        ..Default::default()
    };
    let flat = exp.get_2d_peak_data_im(AreaBounds::default(), 1).unwrap();
    assert_eq!(flat.rt, [-1.0, -1.0]);
    assert_eq!(flat.ion_mobility, [0.4, 0.5]);
    // The row-grouping export keeps the source's cursor semantics, so the same
    // run cannot open a row and is refused rather than writing past the end.
    assert!(matches!(
        exp.get_2d_peak_data_im_per_spectrum(AreaBounds::default(), 1),
        Err(Error::InvalidValue(ref message)) if message.contains("no existing output row")
    ));
}

#[test]
fn mobility_exports_append_and_continue_the_callers_last_row() {
    let exp = MSExperiment {
        spectra: vec![im_scan(-1., 1, &[(100., 10.)], &[0.4])],
        ..Default::default()
    };
    // Appending a selection whose first RT is -1 continues the caller's final
    // row, and the source fetches no array for it, so the value is the -1
    // sentinel rather than 0.4.
    let mut rows = SpectrumPeakDataIm {
        rt: vec![7.0],
        mz: vec![vec![1.0]],
        intensity: vec![vec![2.0]],
        ion_mobility: vec![vec![3.0]],
    };
    exp.append_2d_peak_data_im_per_spectrum(AreaBounds::default(), 1, &mut rows)
        .unwrap();
    assert_eq!(rows.rt, [7.0]);
    assert_eq!(rows.mz, [vec![1.0, 100.0]]);
    assert_eq!(rows.intensity, [vec![2.0, 10.0]]);
    assert_eq!(rows.ion_mobility, [vec![3.0, -1.0]]);

    // The flat export appends to all four columns and keeps their contents.
    let mut flat = FlatPeakDataIm {
        rt: vec![7.0],
        mz: vec![1.0],
        intensity: vec![2.0],
        ion_mobility: vec![3.0],
    };
    exp.append_2d_peak_data_im(AreaBounds::default(), 1, &mut flat)
        .unwrap();
    assert_eq!(flat.rt, [7.0, -1.0]);
    assert_eq!(flat.mz, [1.0, 100.0]);
    assert_eq!(flat.intensity, [2.0, 10.0]);
    assert_eq!(flat.ion_mobility, [3.0, 0.4]);

    // An empty selection leaves the caller's output untouched.
    let before = flat.clone();
    exp.append_2d_peak_data_im(bounds((100., 200.), (0., 1.)), 1, &mut flat)
        .unwrap();
    assert_eq!(flat, before);
}

#[test]
fn merged_rows_use_the_first_spectrum_ion_mobility_array() {
    // Two spectra sharing one retention time merge into a single row, and the
    // source fetches the ion mobility array only when the row starts, so the
    // second spectrum's peaks are looked up in the first spectrum's array.
    let exp = MSExperiment {
        spectra: vec![
            im_scan(1., 1, &[(100., 10.), (200., 20.)], &[0.1, 0.2]),
            im_scan(1., 1, &[(300., 30.), (400., 40.)], &[0.3, 0.4]),
        ],
        ..Default::default()
    };
    let out = exp
        .get_2d_peak_data_im_per_spectrum(AreaBounds::default(), 1)
        .unwrap();
    assert_eq!(out.rt, [1.0]);
    assert_eq!(out.mz, [vec![100.0, 200.0, 300.0, 400.0]]);
    assert_eq!(out.ion_mobility, [vec![0.1, 0.2, 0.1, 0.2]]);
    // The flat export is unaffected: every peak reads its own array.
    let flat = exp.get_2d_peak_data_im(AreaBounds::default(), 1).unwrap();
    assert_eq!(flat.ion_mobility, [0.1, 0.2, 0.3, 0.4]);
}

#[test]
fn short_or_nonfinite_ion_mobility_arrays_are_refused_atomically() {
    // The source indexes the array with the peak's index inside its spectrum
    // and never checks the length, reading out of bounds.
    let short = MSExperiment {
        spectra: vec![im_scan(1., 1, &[(100., 10.), (200., 20.)], &[0.4])],
        ..Default::default()
    };
    assert!(matches!(
        short.get_2d_peak_data_im(AreaBounds::default(), 1),
        Err(Error::InvalidValue(ref message)) if message.contains("shorter than")
    ));
    let mut rows = SpectrumPeakDataIm::default();
    assert!(
        short
            .append_2d_peak_data_im_per_spectrum(AreaBounds::default(), 1, &mut rows)
            .is_err()
    );
    assert_eq!(rows, SpectrumPeakDataIm::default());
    // An empty but present array is the same case; the source's unit test only
    // reaches it through a name match.
    let empty = MSExperiment {
        spectra: vec![im_scan(1., 1, &[(100., 10.)], &[])],
        ..Default::default()
    };
    assert!(empty.get_2d_peak_data_im(AreaBounds::default(), 1).is_err());
    let nonfinite = MSExperiment {
        spectra: vec![im_scan(1., 1, &[(100., 10.)], &[f32::NAN])],
        ..Default::default()
    };
    assert!(
        nonfinite
            .get_2d_peak_data_im(AreaBounds::default(), 1)
            .is_err()
    );
    // An array whose name does not denote ion mobility is not consulted, so the
    // sentinel is written and the length never matters.
    let unrelated = MSExperiment {
        spectra: vec![MSSpectrum {
            rt: 1.,
            ms_level: 1,
            peaks: vec![Peak1D::new(100., 10.), Peak1D::new(200., 20.)],
            float_data_arrays: vec![DataArray::new("signal to noise", vec![1.0])],
            ..Default::default()
        }],
        ..Default::default()
    };
    assert_eq!(
        unrelated
            .get_2d_peak_data_im(AreaBounds::default(), 1)
            .unwrap()
            .ion_mobility,
        [-1.0, -1.0]
    );
}

#[test]
fn mobility_export_alignment_and_ceilings_are_checked_before_writing() {
    let exp = im_export_fixture();
    let mut misaligned = FlatPeakDataIm {
        rt: vec![1.0],
        mz: vec![1.0],
        intensity: vec![1.0],
        ion_mobility: Vec::new(),
    };
    assert!(matches!(
        exp.append_2d_peak_data_im(AreaBounds::default(), 1, &mut misaligned),
        Err(Error::InvalidValue(ref message)) if message.contains("not aligned")
    ));
    let mut misaligned_rows = SpectrumPeakDataIm {
        rt: vec![1.0],
        mz: vec![vec![1.0]],
        intensity: vec![vec![1.0]],
        ion_mobility: vec![Vec::new()],
    };
    assert!(matches!(
        exp.append_2d_peak_data_im_per_spectrum(AreaBounds::default(), 1, &mut misaligned_rows),
        Err(Error::InvalidValue(ref message)) if message.contains("not aligned")
    ));
    let zero = ExperimentMobilityLimits {
        max_spectra: 0,
        max_peaks: 0,
        max_output_points: 0,
        max_output_rows: 0,
        max_work: 0,
        max_bytes: 0,
    };
    for limits in [
        zero,
        ExperimentMobilityLimits {
            max_output_points: 1,
            ..Default::default()
        },
        ExperimentMobilityLimits {
            max_output_rows: 1,
            ..Default::default()
        },
        ExperimentMobilityLimits {
            max_bytes: 0,
            ..Default::default()
        },
        ExperimentMobilityLimits {
            max_work: 8,
            ..Default::default()
        },
    ] {
        assert!(
            matches!(
                exp.get_2d_peak_data_im_with_limits(AreaBounds::default(), 1, limits),
                Err(Error::InvalidValue(ref message)) if message.contains("limit exceeded")
            ) || matches!(
                exp.get_2d_peak_data_im_per_spectrum_with_limits(
                    AreaBounds::default(),
                    1,
                    limits
                ),
                Err(Error::InvalidValue(ref message)) if message.contains("limit exceeded")
            )
        );
    }
    // An empty run can use all-zero ceilings.
    assert_eq!(
        MSExperiment::default()
            .get_2d_peak_data_im_with_limits(AreaBounds::default(), 1, zero)
            .unwrap(),
        FlatPeakDataIm::default()
    );
}

// ---------------------------------------------------------------------------
// MSExperiment.h: rasterizeRTMZ
// ---------------------------------------------------------------------------

fn raster_fixture() -> MSExperiment {
    MSExperiment {
        spectra: vec![
            MSSpectrum {
                rt: 1.,
                ms_level: 1,
                peaks: vec![
                    Peak1D::new(100., 10.),
                    Peak1D::new(110., 30.),
                    Peak1D::new(250., 20.),
                ],
                ..Default::default()
            },
            MSSpectrum {
                rt: 2.,
                ms_level: 1,
                peaks: vec![Peak1D::new(100., 5.)],
                ..Default::default()
            },
            MSSpectrum {
                rt: 3.,
                ms_level: 1,
                peaks: vec![Peak1D::new(300., 7.)],
                ..Default::default()
            },
            MSSpectrum {
                rt: 3.,
                ms_level: 2,
                peaks: vec![Peak1D::new(100., 1000.)],
                ..Default::default()
            },
        ],
        ..Default::default()
    }
}

#[test]
fn rasterize_rt_mz_bins_clamps_edges_and_aggregates() {
    let exp = raster_fixture();
    let raster = RtMzRaster::new(2, 2, 1., 3., 100., 300., 1);
    // Row-major, m/z slowest: [ (mz0,rt0), (mz0,rt1), (mz1,rt0), (mz1,rt1) ].
    // The RT-3 and m/z-300 peak lands in the last bin of both axes rather than
    // past the end, and the MS2 spectrum does not contribute.
    assert_eq!(
        exp.rasterize_rt_mz(&raster).unwrap(),
        [40.0, 5.0, 20.0, 7.0]
    );
    assert_eq!(
        exp.rasterize_rt_mz(&raster.with_aggregation(RasterAggregation::Max))
            .unwrap(),
        [30.0, 5.0, 20.0, 7.0]
    );
    // Peaks outside the m/z window are excluded, not clamped into it.
    let narrow = RtMzRaster::new(1, 1, 1., 3., 100., 150., 1);
    assert_eq!(exp.rasterize_rt_mz(&narrow).unwrap(), [45.0]);
    // An MS level is matched exactly.
    assert_eq!(
        exp.rasterize_rt_mz(&RtMzRaster::new(1, 1, 1., 3., 100., 300., 2))
            .unwrap(),
        [1000.0]
    );
    assert_eq!(
        exp.rasterize_rt_mz(&RtMzRaster::new(1, 1, 1., 3., 100., 300., 3))
            .unwrap(),
        [0.0]
    );
    // An empty run and an empty window both produce a zeroed image.
    assert_eq!(
        MSExperiment::default().rasterize_rt_mz(&raster).unwrap(),
        [0.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(
        exp.rasterize_rt_mz(&RtMzRaster::new(2, 2, 20., 25., 100., 300., 1))
            .unwrap(),
        [0.0, 0.0, 0.0, 0.0]
    );
}

#[test]
fn rasterize_rt_mz_rejects_invalid_grids_and_unsorted_input() {
    let exp = raster_fixture();
    assert!(matches!(
        exp.rasterize_rt_mz(&RtMzRaster::new(0, 2, 1., 3., 100., 300., 1)),
        Err(Error::InvalidValue(ref message)) if message.contains("RT bins")
    ));
    assert!(matches!(
        exp.rasterize_rt_mz(&RtMzRaster::new(2, 0, 1., 3., 100., 300., 1)),
        Err(Error::InvalidValue(ref message)) if message.contains("m/z bins")
    ));
    assert!(matches!(
        exp.rasterize_rt_mz(&RtMzRaster::new(2, 2, 3., 3., 100., 300., 1)),
        Err(Error::InvalidRange(_))
    ));
    assert!(matches!(
        exp.rasterize_rt_mz(&RtMzRaster::new(2, 2, 1., 3., 300., 100., 1)),
        Err(Error::InvalidRange(_))
    ));
    assert!(
        exp.rasterize_rt_mz(&RtMzRaster::new(2, 2, f64::NAN, 3., 100., 300., 1))
            .is_err()
    );
    // The pixel ceiling and the overflowing product are both refused, and
    // nothing is allocated for either.
    assert!(matches!(
        exp.rasterize_rt_mz(&RtMzRaster::new(5000, 5000, 1., 3., 100., 300., 1)),
        Err(Error::InvalidValue(ref message)) if message.contains("MAX_RASTER_PIXELS")
    ));
    let overflow = RtMzRaster::new(usize::MAX, 2, 1., 3., 100., 300., 1);
    assert!(overflow.pixels().is_err());
    assert!(matches!(
        exp.rasterize_rt_mz(&overflow),
        Err(Error::InvalidValue(ref message)) if message.contains("overflows")
    ));
    assert_eq!(
        RtMzRaster::new(3, 4, 1., 3., 100., 300., 1)
            .pixels()
            .unwrap(),
        12
    );
    // Sortedness is a checked precondition, unlike the source's @note.
    let mut unsorted_rt = raster_fixture();
    unsorted_rt.spectra[0].rt = 9.;
    assert!(matches!(
        unsorted_rt.rasterize_rt_mz(&RtMzRaster::new(2, 2, 1., 3., 100., 300., 1)),
        Err(Error::UnsortedData)
    ));
    let mut unsorted_mz = raster_fixture();
    unsorted_mz.spectra[0].peaks.swap(0, 2);
    assert!(matches!(
        unsorted_mz.rasterize_rt_mz(&RtMzRaster::new(2, 2, 1., 3., 100., 300., 1)),
        Err(Error::UnsortedData)
    ));
    let mut nonfinite = raster_fixture();
    nonfinite.spectra[1].peaks[0].intensity = f32::NAN;
    assert!(
        nonfinite
            .rasterize_rt_mz(&RtMzRaster::new(2, 2, 1., 3., 100., 300., 1))
            .is_err()
    );
    assert_eq!(MSExperiment::MAX_RASTER_PIXELS, 16_777_216);
    assert_eq!(MSExperiment::MAX_MOBILITY_ITEMS, 100_000_000);
}

// ---------------------------------------------------------------------------
// Scan-mobility area selection: native validation
// ---------------------------------------------------------------------------

#[test]
fn scan_mobility_area_bounds_are_validated_and_drift_times_must_be_finite() {
    let exp = area_fixture();
    let base = AreaOptions::new(AreaBounds::default(), 1);
    assert!(base.with_mobility(2., 1.).is_err());
    assert!(base.with_mobility(f64::NAN, 1.).is_err());
    assert!(base.with_mobility(0., f64::INFINITY).is_err());
    // A zero-width mobility window keeps the scans exactly on it.
    assert_eq!(
        exp.area_iter(base.with_mobility(1.4, 1.4).unwrap())
            .unwrap()
            .map(|p| p.peak.mz)
            .collect::<Vec<_>>(),
        [504., 506.]
    );
    // Without a mobility bound the drift time is never read, so a NaN drift
    // time on an unrelated scan does not affect the traversal.
    let mut nan_drift = area_fixture();
    nan_drift.spectra[4].drift_time = f64::NAN;
    assert_eq!(exp.area_iter(base).unwrap().count(), 8);
    assert_eq!(nan_drift.area_iter(base).unwrap().count(), 8);
    // With a mobility bound it is refused rather than silently excluded, as
    // RangeBase::contains would do.
    assert!(
        nan_drift
            .area_iter(base.with_mobility(0., 2.).unwrap())
            .is_err()
    );
    // An IM frame carries per-peak mobility but no scalar drift time, so the
    // source sentinel -1 keeps it out of any window above -1.
    let frame = MSExperiment {
        spectra: vec![im_scan(1., 1, &[(100., 10.)], &[0.5])],
        ..Default::default()
    };
    assert_eq!(
        frame
            .area_iter(base.with_mobility(0., 1.).unwrap())
            .unwrap()
            .count(),
        0
    );
    assert_eq!(
        frame
            .area_iter(base.with_mobility(-1., 1.).unwrap())
            .unwrap()
            .count(),
        1
    );
    assert_eq!(frame.area_iter(base).unwrap().count(), 1);
}

// ---------------------------------------------------------------------------
// AreaIterator.h: getRT, and the MS-level cases of the overall test
// ---------------------------------------------------------------------------

#[test]
fn source_area_iterator_rt_and_ms_level_cases() {
    let exp = area_fixture();
    // START_SECTION((CoordinateType getRT() const)): RTBegin(3), RTEnd(9),
    // lowMZ(503), highMZ(509). Four peaks across two scans, with the m/z and
    // the scan RT checked at every step.
    let observed: Vec<_> = exp
        .area_iter(AreaOptions::new(bounds((3., 9.), (503., 509.)), 1))
        .unwrap()
        .map(|p| (p.peak.mz, p.spectrum.rt))
        .collect();
    assert_eq!(observed, [(504., 4.), (506., 4.), (504.1, 8.), (506.1, 8.)]);

    // [EXTRA] Overall test: an explicit `msLevel(1)` selects the same eight
    // peaks as the default level, and level 3 selects nothing.
    assert_eq!(
        exp.area_iter(
            AreaOptions::new(bounds((0., 15.), (500., 520.)), 1)
                .with_mobility(0., 2.)
                .unwrap()
        )
        .unwrap()
        .count(),
        8
    );
    assert_eq!(
        exp.area_iter(AreaOptions::new(bounds((0., 15.), (500., 520.)), 3))
            .unwrap()
            .count(),
        0
    );
    // "Test with empty (no MS level 1) experiment": every scan at level 2.
    let mut level_two = area_fixture();
    for spectrum in &mut level_two.spectra {
        spectrum.ms_level = 2;
    }
    assert_eq!(
        level_two
            .area_iter(AreaOptions::new(bounds((0., 15.), (500., 520.)), 1))
            .unwrap()
            .count(),
        0
    );
    assert_eq!(
        level_two
            .area_iter(AreaOptions::new(bounds((0., 15.), (500., 520.)), 2))
            .unwrap()
            .count(),
        8
    );
}

// ---------------------------------------------------------------------------
// MSExperiment.h: the backward-compatible combined-range delegates
// ---------------------------------------------------------------------------

#[test]
fn source_backward_compatible_range_delegates_include_scan_mobility() {
    // START_SECTION((Backward compatibility tests)): one MS1 spectrum with a
    // scalar drift time and one chromatogram with two points. getMinRT() and
    // its siblings read the combined range manager, so the port answers them
    // through combined_range_manager().
    let mut exp = MSExperiment::default();
    exp.spectra.push(MSSpectrum {
        ms_level: 1,
        rt: 30.,
        drift_time: 50.,
        peaks: vec![Peak1D::new(100., 1000.)],
        ..Default::default()
    });
    exp.chromatograms.push(openms::kernel::MSChromatogram {
        peaks: vec![
            openms::kernel::ChromatogramPeak::new(10., 500.),
            openms::kernel::ChromatogramPeak::new(20., 1500.),
        ],
        // The source also stores a "product_mz" meta value here; neither
        // implementation reads it for the m/z range, which is why the minimum
        // m/z below is the unset product m/z of zero.
        ..Default::default()
    });
    let combined = exp.combined_range_manager().unwrap();
    assert_eq!(combined.min_rt().unwrap(), 10.0);
    assert_eq!(combined.max_rt().unwrap(), 30.0);
    assert_eq!(combined.min_mz().unwrap(), 0.0);
    assert_eq!(combined.max_mz().unwrap(), 100.0);
    assert_eq!(combined.min_intensity().unwrap(), 500.0);
    assert_eq!(combined.max_intensity().unwrap(), 1500.0);
    assert_eq!(combined.min_mobility().unwrap(), 50.0);
    assert_eq!(combined.max_mobility().unwrap(), 50.0);
}
