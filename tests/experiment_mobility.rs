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
use openms::kernel::ranges::{HasRangeType, MSDim, RangeBase, RangeManager};
use openms::kernel::spectrum_mobility::RasterAggregation;
use openms::kernel::{
    AreaBounds, AreaOptions, ChromatogramPeak, DataArray, MSChromatogram, MSExperiment, MSSpectrum,
    MzAggregation, MzRtRegion, Peak1D,
};
use openms::metadata::ContactPerson;

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

// ---------------------------------------------------------------------------
// MSExperiment_test.cpp sections whose literals no earlier package asserted.
//
// Every section below carries more than five assertion macros, so the mapping
// rule requires a port rather than a citation. The fixtures and the asserted
// values are transcribed from the pinned `MSExperiment_test.cpp`; the source
// line range appears at each test.
// ---------------------------------------------------------------------------

/// One contact named "Name" plus one spectrum holding peaks at m/z 5 and 10 —
/// the fixture both assignment sections build (`MSExperiment_test.cpp:117-126`
/// and `:144-153`).
fn assignment_fixture() -> MSExperiment {
    let mut settings = openms::metadata::ExperimentalSettings::new();
    settings.contacts = vec![ContactPerson {
        first_name: "Name".into(),
        ..Default::default()
    }];
    MSExperiment {
        spectra: vec![MSSpectrum {
            peaks: vec![Peak1D::new(5., 0.), Peak1D::new(10., 0.)],
            ..Default::default()
        }],
        settings,
        ..Default::default()
    }
}

#[test]
fn source_copy_assignment_carries_contacts_spectra_and_mz_range() {
    // START_SECTION((MSExperiment& operator= (const MSExperiment& source)))
    // MSExperiment_test.cpp:115-140. `Clone` is the port's copy assignment.
    let source = assignment_fixture();
    let mut copy = MSExperiment::default();
    assert_eq!(copy.len(), 0); // `PeakMap tmp2;` starts empty
    copy = source.clone();
    assert_eq!(copy.settings.contacts.len(), 1);
    assert_eq!(copy.settings.contacts[0].first_name, "Name");
    assert_eq!(copy.len(), 1);
    let ranges = copy.combined_range_manager().unwrap();
    assert_eq!(ranges.min_mz().unwrap(), 5.0);
    assert_eq!(ranges.max_mz().unwrap(), 10.0);
    // `tmp2 = PeakMap();` — assignment from a fresh temporary empties the target.
    copy = MSExperiment::default();
    assert_eq!(copy.settings.contacts.len(), 0);
    assert_eq!(copy.len(), 0);
    // The source of a copy is untouched, which is what separates this section
    // from the move-assignment one below.
    assert_eq!(source.len(), 1);
    assert_eq!(source.settings.contacts[0].first_name, "Name");
    // `operator!=` (MSExperiment_test.cpp:191-203): a run that differs in
    // contacts or in spectrum count is not equal to an empty one.
    assert_ne!(copy, source);
    assert_eq!(copy, MSExperiment::default());
}

#[test]
fn source_move_assignment_transfers_everything_and_empties_the_origin() {
    // START_SECTION((MSExperiment& operator= (const MSExperiment&& source)))
    // MSExperiment_test.cpp:142-174. Rust moves by value, so the section's
    // `PeakMap tmp2 = std::move(tmp); TEST_EQUAL(tmp.size(), 0)` pair needs a
    // move that leaves the origin observable: `std::mem::take` is exactly the
    // source's move-assign-then-default-the-origin behaviour.
    let mut origin = assignment_fixture();
    let original = origin.clone();
    let moved = std::mem::take(&mut origin);
    assert_eq!(moved, original); // should be equal to the original
    assert_eq!(moved.settings.contacts.len(), 1);
    assert_eq!(moved.settings.contacts[0].first_name, "Name");
    assert_eq!(moved.len(), 1);
    let ranges = moved.combined_range_manager().unwrap();
    assert_eq!(ranges.min_mz().unwrap(), 5.0);
    assert_eq!(ranges.max_mz().unwrap(), 10.0);
    // test move
    assert_eq!(origin.len(), 0);
    // `tmp2 = PeakMap();` — rvalue assignment over a populated run.
    let mut target = moved;
    assert_eq!(target.len(), 1);
    target = MSExperiment::default();
    assert_eq!(target.settings.contacts.len(), 0);
    assert_eq!(target.len(), 0);
    // A plain Rust move of the whole value keeps the peaks too, so the
    // section's primary claim does not rest on `mem::take` alone.
    let relocated = original;
    assert_eq!(relocated.spectra[0].peaks[1].mz, 10.0);
}

/// One spectrum with one peak and a scalar drift time — the shape the
/// `updateRanges` fixture of `MSExperiment_test.cpp:431-468` repeats.
fn im_peak_scan(rt: f64, drift_time: f64, level: u32, mz: f64, intensity: f32) -> MSSpectrum {
    MSSpectrum {
        rt,
        drift_time,
        ms_level: level,
        peaks: vec![Peak1D::new(mz, intensity)],
        ..Default::default()
    }
}

/// One spectrum with one peak and no scalar drift time.
fn plain_peak_scan(rt: f64, level: u32, mz: f64, intensity: f32) -> MSSpectrum {
    MSSpectrum {
        rt,
        ms_level: level,
        peaks: vec![Peak1D::new(mz, intensity)],
        ..Default::default()
    }
}

/// A chromatogram whose `Product` m/z is set, as `chrom1.setProduct(prod1)` does.
fn product_chromatogram(mz: f64, points: &[(f64, f32)]) -> MSChromatogram {
    let mut chromatogram = MSChromatogram {
        peaks: points
            .iter()
            .map(|&(rt, intensity)| ChromatogramPeak::new(rt, intensity))
            .collect(),
        ..Default::default()
    };
    chromatogram.product.mz = mz;
    chromatogram
}

#[test]
fn source_update_ranges_combined_per_level_and_chromatogram_literals() {
    // START_SECTION((virtual void updateRanges())) MSExperiment_test.cpp:429-603,
    // the largest section of the class test (64 assertion macros). The port has
    // no `updateRanges()`; the three managers are computed on demand, so the
    // section's repeated "second time to check the initialization" becomes a
    // second query.
    let mut exp = MSExperiment {
        spectra: vec![
            im_peak_scan(30., 99., 1, 5.0, -5.0),
            im_peak_scan(40., 99., 1, 7.0, -7.0),
            im_peak_scan(45., 199., 3, 9.0, -10.0),
            im_peak_scan(50., 66., 3, 10.0, -9.0),
        ],
        ..Default::default()
    };
    for _ in 0..2 {
        let combined = exp.combined_range_manager().unwrap();
        assert_eq!(combined.min_mz().unwrap(), 5.0);
        assert_eq!(combined.max_mz().unwrap(), 10.0);
        assert_eq!(combined.min_intensity().unwrap(), -10.0);
        assert_eq!(combined.max_intensity().unwrap(), -5.0);
        assert_eq!(combined.min_rt().unwrap(), 30.0);
        assert_eq!(combined.max_rt().unwrap(), 50.0);
        assert_eq!(combined.min_mobility().unwrap(), 66.0);
        assert_eq!(combined.max_mobility().unwrap(), 199.0);
        assert_eq!(exp.ms_levels(), [1, 3]);
        assert_eq!(exp.total_peak_count().unwrap(), 4);
    }
    // The MS1 slice of the per-level manager, asserted twice by the source's
    // own `for (int l = 0; l < 2; ++l)` loop.
    let initial_ms_levels = exp.ms_levels();
    for _ in 0..2 {
        let spectra = exp.spectrum_range_manager().unwrap();
        let level1 = spectra.by_ms_level(1).unwrap();
        assert_eq!(level1.min_mz().unwrap(), 5.0);
        assert_eq!(level1.max_mz().unwrap(), 7.0);
        assert_eq!(level1.min_intensity().unwrap(), -7.0);
        assert_eq!(level1.max_intensity().unwrap(), -5.0);
        assert_eq!(level1.min_rt().unwrap(), 30.0);
        assert_eq!(level1.max_rt().unwrap(), 40.0);
        assert_eq!(level1.min_mobility().unwrap(), 99.0);
        assert_eq!(level1.max_mobility().unwrap(), 99.0);
        assert_eq!(exp.ms_levels(), initial_ms_levels);
        assert_eq!(exp.total_peak_count().unwrap(), 4);
    }

    // "test with only one peak": MSExperiment_test.cpp:524-555.
    exp = MSExperiment {
        spectra: vec![im_peak_scan(30., 99., 1, 5.0, -5.0)],
        ..Default::default()
    };
    let combined = exp.combined_range_manager().unwrap();
    assert_eq!(combined.min_mz().unwrap(), 5.0);
    assert_eq!(combined.max_mz().unwrap(), 5.0);
    assert_eq!(combined.min_intensity().unwrap(), -5.0);
    assert_eq!(combined.max_intensity().unwrap(), -5.0);
    assert_eq!(combined.min_rt().unwrap(), 30.0);
    assert_eq!(combined.max_rt().unwrap(), 30.0);
    assert_eq!(combined.min_mobility().unwrap(), 99.0);
    assert_eq!(combined.max_mobility().unwrap(), 99.0);
    let spectra = exp.spectrum_range_manager().unwrap();
    let global = spectra.global();
    assert_eq!(global.min_mz().unwrap(), 5.0);
    assert_eq!(global.max_mz().unwrap(), 5.0);
    assert_eq!(global.min_intensity().unwrap(), -5.0);
    assert_eq!(global.max_intensity().unwrap(), -5.0);
    assert_eq!(global.min_rt().unwrap(), 30.0);
    assert_eq!(global.max_rt().unwrap(), 30.0);
    assert_eq!(global.min_mobility().unwrap(), 99.0);
    assert_eq!(global.max_mobility().unwrap(), 99.0);

    // "test ranges with a chromatogram": MSExperiment_test.cpp:557-600. These
    // two chromatograms carry real `Product` m/z values, 100 and 80.
    exp.chromatograms = vec![
        product_chromatogram(100.0, &[(0.3, 10.0), (0.2, 10.2)]),
        product_chromatogram(80.0, &[(0.2, 10.2), (0.1, 10.4)]),
    ];
    let combined = exp.combined_range_manager().unwrap();
    assert_eq!(combined.min_mz().unwrap(), 5.0);
    assert_eq!(combined.max_mz().unwrap(), 100.0);
    assert_eq!(combined.min_intensity().unwrap(), -5.0);
    assert_eq!(combined.max_intensity().unwrap(), f64::from(10.4_f32));
    assert_eq!(combined.min_rt().unwrap(), 0.1);
    assert_eq!(combined.max_rt().unwrap(), 30.0); // overall range still 30
    let chromatograms = exp.chromatogram_range_manager().unwrap();
    assert_eq!(chromatograms.min_mz().unwrap(), 80.0);
    assert_eq!(chromatograms.max_mz().unwrap(), 100.0);
    assert_eq!(chromatograms.min_intensity().unwrap(), 10.0);
    assert_eq!(chromatograms.max_intensity().unwrap(), f64::from(10.4_f32));
    assert_eq!(chromatograms.min_rt().unwrap(), 0.1);
    assert_eq!(chromatograms.max_rt().unwrap(), 0.3); // chromatogram range 0.1-0.3
}

/// The source's `createPeakMapWithRTs`, optionally followed by its
/// `setMSLevel`: peakless spectra at the given retention times and levels.
fn rt_only_experiment(rts: &[f64], levels: &[u32]) -> MSExperiment {
    MSExperiment {
        spectra: rts
            .iter()
            .zip(levels)
            .map(|(&rt, &ms_level)| MSSpectrum {
                rt,
                ms_level,
                ..Default::default()
            })
            .collect(),
        ..Default::default()
    }
}

#[test]
fn source_closest_spectrum_in_rt_literals_any_level_and_per_level() {
    // START_SECTION(ConstIterator getClosestSpectrumInRT(const double RT) const)
    // MSExperiment_test.cpp:768-805. One Rust function answers both source
    // overloads, with `ms_level == 0` as the "any level" wildcard the
    // single-argument overload implies.
    let exp = rt_only_experiment(&[30., 40., 45., 50.], &[1, 1, 1, 1]);
    for (rt, expected) in [
        (-200.0, 0),
        (20.0, 0),
        (31.0, 0),
        (34.9, 0),
        (39.0, 1),
        (41.0, 1),
        (42.4, 1),
        (44.0, 2),
        (47.0, 2),
        (47.6, 3),
        (51.0, 3),
        (5_100_000.0, 3),
    ] {
        assert_eq!(exp.closest_spectrum_in_rt(rt, 0).unwrap(), Some(expected));
    }
    assert_eq!(
        MSExperiment::default()
            .closest_spectrum_in_rt(47.6, 0)
            .unwrap(),
        None
    );
    // The same `{30, 40, 45, 50}` fixture carries the four `RTBegin`/`RTEnd`
    // sections — the mutable pair at MSExperiment_test.cpp:738-765 and the
    // const pair at :958-1007, which assert the same literals. One Rust index
    // serves all four.
    assert_eq!(exp.spectra[exp.rt_begin(20.).unwrap()].rt, 30.0);
    assert_eq!(exp.spectra[exp.rt_begin(30.).unwrap()].rt, 30.0);
    assert_eq!(exp.spectra[exp.rt_begin(31.).unwrap()].rt, 40.0);
    assert_eq!(exp.rt_begin(55.).unwrap(), exp.len());
    assert_eq!(exp.spectra[exp.rt_end(20.).unwrap()].rt, 30.0);
    assert_eq!(exp.spectra[exp.rt_end(30.).unwrap()].rt, 40.0);
    assert_eq!(exp.spectra[exp.rt_end(31.).unwrap()].rt, 40.0);
    assert_eq!(exp.rt_end(55.).unwrap(), exp.len());

    // START_SECTION(ConstIterator getClosestSpectrumInRT(const double RT, UInt
    // ms_level) const) MSExperiment_test.cpp:824-883.
    let exp = rt_only_experiment(
        &[30., 31., 32., 40., 41., 50., 60., 61.],
        &[1, 2, 2, 1, 2, 1, 1, 2],
    );
    for (rt, level, expected) in [
        (-200.0, 1, 0),
        (-200.0, 2, 1),
        (20.0, 1, 0),
        (31.0, 1, 0),
        (34.9, 1, 0),
        (20.0, 2, 1),
        (31.0, 2, 1),
        (31.4, 2, 1),
        (39.0, 1, 3),
        (41.0, 1, 3),
        (42.4, 1, 3),
        (45.5, 1, 5),
        (49.0, 1, 5),
        (54.5, 1, 5),
        (55.1, 1, 6),
        (59.1, 1, 6),
        (5_100_000.0, 1, 6),
        (58.0, 2, 7),
        (63.0, 2, 7),
        (5_100_000.0, 2, 7),
    ] {
        assert_eq!(
            exp.closest_spectrum_in_rt(rt, level).unwrap(),
            Some(expected)
        );
    }
    assert_eq!(
        MSExperiment::default()
            .closest_spectrum_in_rt(47.6, 1)
            .unwrap(),
        None
    );
    // The one assertion of this section the port does not reproduce: the
    // source's two-argument overload answers `cend()` for `ms_level == 0`
    // because no scan carries that level, while the port reserves `0` as the
    // wildcard standing in for the single-argument overload, so it answers the
    // first scan of any level. Recorded under *Native differences* in
    // docs/EXPERIMENT_MOBILITY_SUPPORT.md.
    assert_eq!(exp.closest_spectrum_in_rt(-200.0, 0).unwrap(), Some(0));
}

#[test]
fn source_is_sorted_literals_separate_rt_and_mz_checks() {
    // START_SECTION(bool isSorted(bool check_mz = true ) const)
    // MSExperiment_test.cpp:1046-1090.
    let peaks = || {
        vec![
            Peak1D::new(1000., 1.),
            Peak1D::new(1001., 1.),
            Peak1D::new(1002., 1.),
        ]
    };
    let mut exp = MSExperiment {
        spectra: vec![
            MSSpectrum {
                rt: 1.,
                peaks: peaks(),
                ..Default::default()
            },
            MSSpectrum {
                rt: 2.,
                peaks: peaks(),
                ..Default::default()
            },
        ],
        ..Default::default()
    };
    // The source labels its first block "test with identical RTs" but assigns
    // 1.0 and 2.0, so it is the same ascending case its next block repeats.
    for _ in 0..2 {
        assert!(exp.is_sorted(false));
        assert!(exp.is_sorted(true));
    }
    // "test with a reversed spectrum": the retention-time order still holds.
    exp.spectra[0].peaks.reverse();
    assert!(exp.is_sorted(false));
    assert!(!exp.is_sorted(true));
    // "test with reversed RTs".
    exp.spectra.reverse();
    assert!(!exp.is_sorted(false));
    assert!(!exp.is_sorted(true));
}

#[test]
fn source_precursor_spectrum_literals_for_both_overloads() {
    // START_SECTION((ConstIterator getPrecursorSpectrum(ConstIterator) const))
    // MSExperiment_test.cpp:1149-1178, and the `int` overload at :1181-1208.
    // One Rust function serves both: `None` replaces the past-the-end iterator
    // and the `-1`.
    let mut exp = MSExperiment {
        spectra: vec![MSSpectrum::default(); 10],
        ..Default::default()
    };
    for (index, level) in [1_u32, 2, 1, 2, 2].into_iter().enumerate() {
        exp.spectra[index].ms_level = level;
    }
    for (index, expected) in [
        (0, None),
        (1, Some(0)),
        (2, None),
        (3, Some(2)),
        (4, Some(2)),
    ] {
        assert_eq!(exp.precursor_spectrum_index(index).unwrap(), expected);
    }
    // `getPrecursorSpectrum(exp.end()) == exp.end()`: the source hands the
    // past-the-end iterator straight back, while the port rejects an
    // out-of-range index instead of inventing an answer for it.
    assert!(exp.precursor_spectrum_index(exp.len()).is_err());

    for (index, level) in [2_u32, 1, 1, 1, 1].into_iter().enumerate() {
        exp.spectra[index].ms_level = level;
    }
    for index in 0..5 {
        assert_eq!(exp.precursor_spectrum_index(index).unwrap(), None);
    }
    assert!(exp.precursor_spectrum_index(exp.len()).is_err());
}

#[test]
fn source_swap_exchanges_comment_spectra_levels_and_ranges() {
    // START_SECTION((void swap(MSExperiment &from)))
    // MSExperiment_test.cpp:1356-1380. `std::mem::swap` is the counterpart.
    let mut first = MSExperiment {
        spectra: vec![MSSpectrum {
            ms_level: 2,
            peaks: vec![Peak1D::new(0., 0.5), Peak1D::new(0., 1.7)],
            ..Default::default()
        }],
        ..Default::default()
    };
    first.settings.comment = "stupid comment".into();
    let mut second = MSExperiment::default();
    std::mem::swap(&mut first, &mut second);

    assert_eq!(first.settings.comment, "");
    assert_eq!(first.len(), 0);
    assert_eq!(
        first.combined_range_manager().unwrap().has_range(),
        HasRangeType::None
    );
    assert_eq!(first.ms_levels().len(), 0);
    assert_eq!(first.total_peak_count().unwrap(), 0);

    assert_eq!(second.settings.comment, "stupid comment");
    assert_eq!(second.len(), 1);
    assert_eq!(
        second
            .combined_range_manager()
            .unwrap()
            .min_intensity()
            .unwrap(),
        0.5
    );
    assert_eq!(second.ms_levels().len(), 1);
    assert_eq!(second.total_peak_count().unwrap(), 2);
}

#[test]
fn source_add_chromatogram_appends_without_touching_earlier_entries() {
    // START_SECTION((void addChromatogram(const MSChromatogram&)))
    // MSExperiment_test.cpp:1510-1534. The port's counterpart is
    // `chromatograms.push`, so the section pins that the push appends and
    // leaves the existing entry byte-identical.
    let first = product_chromatogram(0., &[(0.1, 10.0), (0.2, 10.2)]);
    let second = product_chromatogram(0., &[(0.2, 10.2), (0.3, 10.4)]);
    let mut exp = MSExperiment::default();
    assert_eq!(exp.chromatograms.len(), 0);
    exp.chromatograms.push(first.clone());
    assert_eq!(exp.chromatograms.len(), 1);
    assert_eq!(exp.chromatograms[0], first);
    exp.chromatograms.push(second.clone());
    assert_eq!(exp.chromatograms.len(), 2);
    assert_eq!(exp.chromatograms[0], first);
    assert_eq!(exp.chromatograms[1], second);

    // START_SECTION((void setChromatograms(const std::vector<MSChromatogram>&)))
    // MSExperiment_test.cpp:1485-1507: assigning the container replaces it
    // wholesale and preserves each entry. The target starts with an unrelated
    // chromatogram, so the assertion below also shows the replacement.
    let mut assigned = MSExperiment {
        chromatograms: vec![product_chromatogram(999., &[(9., 9.)])],
        ..Default::default()
    };
    assigned.chromatograms = vec![first.clone(), second.clone()];
    assert_eq!(assigned.chromatograms.len(), 2);
    assert_eq!(assigned.chromatograms[0], first);
    assert_eq!(assigned.chromatograms[1], second);

    // START_SECTION((std::vector<MSChromatogram>& getChromatograms()))
    // MSExperiment_test.cpp:1543-1553: the non-const accessor exists so the
    // caller can swap the container in and out. The port's public field does
    // the same through `std::mem::swap`.
    let mut detached = Vec::new();
    std::mem::swap(&mut assigned.chromatograms, &mut detached);
    assert_eq!(assigned.chromatograms.len(), 0);
    assert_eq!(detached.len(), 2);
    std::mem::swap(&mut assigned.chromatograms, &mut detached);
    assert_eq!(assigned.chromatograms.len(), 2);
    assert_eq!(detached.len(), 0);
}

#[test]
fn source_sort_spectra_reset_and_settings_assignment_literals() {
    // START_SECTION((void sortSpectra(bool sort_mz = true)))
    // MSExperiment_test.cpp:1010-1043. The `set2DData` fixture groups four
    // points into two scans whose peaks arrive out of m/z order.
    let mut exp = MSExperiment {
        spectra: vec![
            MSSpectrum {
                rt: 1.,
                peaks: vec![Peak1D::new(5., 0.), Peak1D::new(3., 0.)],
                ..Default::default()
            },
            MSSpectrum {
                rt: 2.,
                peaks: vec![Peak1D::new(14., 0.), Peak1D::new(11., 0.)],
                ..Default::default()
            },
        ],
        ..Default::default()
    };
    exp.sort_spectra(true).unwrap();
    assert_eq!(exp.spectra[0].peaks[0].mz, 3.0);
    assert_eq!(exp.spectra[0].peaks[1].mz, 5.0);
    assert_eq!(exp.spectra[1].peaks[0].mz, 11.0);
    assert_eq!(exp.spectra[1].peaks[1].mz, 14.0);

    // START_SECTION((MSExperiment& operator=(const ExperimentalSettings&)))
    // MSExperiment_test.cpp:1140-1146, and the two `getExperimentalSettings`
    // accessors at :1124-1138: the port's public `settings` field serves all
    // three.
    let mut labelled = MSExperiment::default();
    labelled.settings.comment = "test".into();
    assert_eq!(labelled.settings.comment, "test");
    // The receiver already holds a spectrum, so the assertion also shows that
    // assigning the settings leaves the run's data alone, as the source's
    // `operator=(const ExperimentalSettings&)` does.
    let mut receiver = MSExperiment {
        spectra: vec![plain_peak_scan(7., 1, 70.0, 700.0)],
        ..Default::default()
    };
    receiver.settings = labelled.settings.clone();
    assert_eq!(receiver.settings.comment, "test");
    assert_eq!(receiver.len(), 1);

    // START_SECTION(void clear(bool clear_meta_data))
    // MSExperiment_test.cpp:1383-1401, and `reset()` at :1092-1121, which the
    // port answers with `clear(true)` because it keeps no range cache for
    // `reset()` to drop in addition.
    let mut edit = MSExperiment {
        spectra: vec![MSSpectrum::default(); 5],
        chromatograms: vec![MSChromatogram::default(); 5],
        ..Default::default()
    };
    edit.settings.sample.name = "bla".into();
    edit.settings.metadata.insert("label".into(), "bla".into());
    edit.clear(false);
    assert_eq!(edit.len(), 0);
    assert!(edit.chromatograms.is_empty());
    assert_ne!(edit, MSExperiment::default()); // the metadata survived
    edit.clear(true);
    assert!(edit.is_empty());
    assert_eq!(edit, MSExperiment::default());
}

/// One spectrum with several peaks.
fn multi_peak_scan(rt: f64, level: u32, points: &[(f64, f32)]) -> MSSpectrum {
    MSSpectrum {
        rt,
        ms_level: level,
        peaks: points
            .iter()
            .map(|&(mz, intensity)| Peak1D::new(mz, intensity))
            .collect(),
        ..Default::default()
    }
}

/// The four-spectrum aggregation fixture of `MSExperiment_test.cpp`, rebuilt
/// here for the two extraction cases below.
fn aggregation_fixture() -> MSExperiment {
    MSExperiment {
        spectra: vec![
            multi_peak_scan(1.0, 1, &[(100., 1000.), (200., 2000.), (300., 3000.)]),
            multi_peak_scan(2.0, 2, &[(150., 1500.), (250., 2500.)]),
            multi_peak_scan(3.0, 1, &[(100., 1100.), (200., 2100.), (300., 3100.)]),
            multi_peak_scan(4.0, 1, &[(100., 1200.), (200., 2200.), (300., 3200.)]),
        ],
        ..Default::default()
    }
}

#[test]
fn source_matrix_xic_extraction_over_two_ranges() {
    // START_SECTION((std::vector<MSChromatogram> extractXICsFromMatrix(...)))
    // MSExperiment_test.cpp:2316-2389. The section's single case is the
    // two-row matrix; `tests/experiment_aggregation.rs::source_matrix_all_reducers`
    // covers the one-row form for every reducer but never this one.
    let exp = aggregation_fixture();
    let xics = exp
        .extract_xics_from_matrix(
            &[[90.0, 110.0, 0.0, 3.5], [190.0, 210.0, 0.0, 5.0]],
            1,
            MzAggregation::Sum,
        )
        .unwrap();
    assert_eq!(xics.len(), 2);
    assert_eq!(
        xics[0]
            .peaks
            .iter()
            .map(|p| p.intensity)
            .collect::<Vec<_>>(),
        [1000.0, 1100.0]
    );
    assert_eq!(
        xics[1]
            .peaks
            .iter()
            .map(|p| p.intensity)
            .collect::<Vec<_>>(),
        [2000.0, 2100.0, 2200.0]
    );
}

#[test]
fn source_custom_reduction_xic_averages_the_whole_mz_window() {
    // START_SECTION((... extractXICs(..., MzReductionFunctionType))) "Test 3:
    // Custom reduction function (average intensity)",
    // MSExperiment_test.cpp:2610-2643. Tests 1 and 2 of that section are
    // asserted by
    // `tests/experiment_aggregation.rs::source_xic_full_rt_product_and_default_metadata`;
    // this third case, whose product m/z is the (90 + 310) / 2 midpoint, is not.
    let exp = aggregation_fixture();
    let xics = exp
        .extract_xics_with(
            &[MzRtRegion::new(90.0, 310.0, 0.0, 5.0).unwrap()],
            1,
            |peaks| MzAggregation::Mean.reduce(peaks),
        )
        .unwrap();
    assert_eq!(xics.len(), 1);
    assert_eq!(
        xics[0]
            .peaks
            .iter()
            .map(|p| (p.rt, p.intensity))
            .collect::<Vec<_>>(),
        [(1.0, 2000.0), (3.0, 2100.0), (4.0, 2200.0)]
    );
    assert_eq!(xics[0].product.mz, (90.0 + 310.0) / 2.0);
}

/// The chromatogram every `updateRanges` sub-case of the dual-range block adds:
/// two points and a `"product_mz"` *meta value* that is not the `Product` m/z,
/// so the m/z the ranges see stays the unset `0`.
fn meta_only_chromatogram() -> MSChromatogram {
    let mut chromatogram = MSChromatogram {
        peaks: vec![
            ChromatogramPeak::new(10., 500.),
            ChromatogramPeak::new(20., 1500.),
        ],
        ..Default::default()
    };
    chromatogram.metadata.insert(
        "product_mz".into(),
        openms::metadata::MetaValue::try_from(305.0).unwrap(),
    );
    chromatogram
}

#[test]
fn source_spectrum_ranges_global_and_per_level_literals() {
    // START_SECTION((const SpectrumRangeManagerType& spectrumRanges() const))
    // MSExperiment_test.cpp:2650-2693.
    let exp = MSExperiment {
        spectra: vec![
            plain_peak_scan(30., 1, 100.0, 1000.0),
            plain_peak_scan(35., 2, 200.0, 2000.0),
        ],
        ..Default::default()
    };
    let spectra = exp.spectrum_range_manager().unwrap();
    let global = spectra.global();
    assert_eq!(global.min_mz().unwrap(), 100.0);
    assert_eq!(global.max_mz().unwrap(), 200.0);
    assert_eq!(global.min_intensity().unwrap(), 1000.0);
    assert_eq!(global.max_intensity().unwrap(), 2000.0);
    let level1 = spectra.by_ms_level(1).unwrap();
    assert_eq!(level1.min_mz().unwrap(), 100.0);
    assert_eq!(level1.max_mz().unwrap(), 100.0);
    assert_eq!(level1.min_intensity().unwrap(), 1000.0);
    assert_eq!(level1.max_intensity().unwrap(), 1000.0);
    let level2 = spectra.by_ms_level(2).unwrap();
    assert_eq!(level2.min_mz().unwrap(), 200.0);
    assert_eq!(level2.max_mz().unwrap(), 200.0);
    assert_eq!(level2.min_intensity().unwrap(), 2000.0);
    assert_eq!(level2.max_intensity().unwrap(), 2000.0);
}

#[test]
fn source_chromatogram_ranges_span_every_chromatogram() {
    // START_SECTION((const ChromatogramRangeManagerType& chromatogramRanges() const))
    // MSExperiment_test.cpp:2696-2732. Two chromatograms whose retention-time
    // and intensity windows interleave, so the manager has to span both.
    let exp = MSExperiment {
        chromatograms: vec![
            product_chromatogram(0., &[(10., 500.), (20., 1500.)]),
            product_chromatogram(0., &[(15., 800.), (25., 1800.)]),
        ],
        ..Default::default()
    };
    let chromatograms = exp.chromatogram_range_manager().unwrap();
    assert_eq!(chromatograms.min_rt().unwrap(), 10.0);
    assert_eq!(chromatograms.max_rt().unwrap(), 25.0);
    assert_eq!(chromatograms.min_intensity().unwrap(), 500.0);
    assert_eq!(chromatograms.max_intensity().unwrap(), 1800.0);
}

#[test]
fn source_dual_range_update_ranges_four_cases_and_three_levels() {
    // START_SECTION((void updateRanges())) MSExperiment_test.cpp:2735-2914, the
    // second section of that name (50 assertion macros). The port computes the
    // three managers on demand, so each sub-case queries them directly.

    // Test case 1: Empty experiment.
    let empty = MSExperiment::default();
    assert_eq!(
        empty.spectrum_range_manager().unwrap().global().has_range(),
        HasRangeType::None
    );
    assert_eq!(
        empty.chromatogram_range_manager().unwrap().has_range(),
        HasRangeType::None
    );
    assert_eq!(
        empty.combined_range_manager().unwrap().has_range(),
        HasRangeType::None
    );

    // Test case 2: Experiment with only spectra.
    let spectra_only = MSExperiment {
        spectra: vec![plain_peak_scan(30., 1, 100.0, 1000.0)],
        ..Default::default()
    };
    let spectra = spectra_only.spectrum_range_manager().unwrap();
    assert_eq!(spectra.global().has_range(), HasRangeType::Some);
    assert_eq!(spectra.global().min_mz().unwrap(), 100.0);
    assert_eq!(spectra.global().max_mz().unwrap(), 100.0);
    assert_eq!(
        spectra_only
            .chromatogram_range_manager()
            .unwrap()
            .has_range(),
        HasRangeType::None
    );
    let combined = spectra_only.combined_range_manager().unwrap();
    assert_eq!(combined.has_range(), HasRangeType::Some);
    assert_eq!(combined.min_mz().unwrap(), 100.0);
    assert_eq!(combined.max_mz().unwrap(), 100.0);
    assert_eq!(combined.min_rt().unwrap(), 30.0);
    assert_eq!(combined.max_rt().unwrap(), 30.0);
    assert_eq!(combined.min_intensity().unwrap(), 1000.0);
    assert_eq!(combined.max_intensity().unwrap(), 1000.0);

    // Test case 3: Experiment with only chromatograms.
    let chromatograms_only = MSExperiment {
        chromatograms: vec![meta_only_chromatogram()],
        ..Default::default()
    };
    assert_eq!(
        chromatograms_only
            .spectrum_range_manager()
            .unwrap()
            .global()
            .has_range(),
        HasRangeType::None
    );
    let chromatograms = chromatograms_only.chromatogram_range_manager().unwrap();
    assert_eq!(chromatograms.has_range(), HasRangeType::All);
    assert_eq!(chromatograms.min_rt().unwrap(), 10.0);
    assert_eq!(chromatograms.max_rt().unwrap(), 20.0);
    let combined = chromatograms_only.combined_range_manager().unwrap();
    assert_eq!(combined.has_range(), HasRangeType::Some);
    assert_eq!(combined.min_rt().unwrap(), 10.0);
    assert_eq!(combined.max_rt().unwrap(), 20.0);
    assert_eq!(combined.min_intensity().unwrap(), 500.0);
    assert_eq!(combined.max_intensity().unwrap(), 1500.0);

    // Test case 4: Experiment with both spectra and chromatograms.
    let both = MSExperiment {
        spectra: vec![plain_peak_scan(30., 1, 100.0, 1000.0)],
        chromatograms: vec![meta_only_chromatogram()],
        ..Default::default()
    };
    let spectra = both.spectrum_range_manager().unwrap();
    assert_eq!(spectra.global().has_range(), HasRangeType::Some);
    assert_eq!(spectra.global().min_mz().unwrap(), 100.0);
    assert_eq!(spectra.global().max_mz().unwrap(), 100.0);
    assert_eq!(spectra.global().min_intensity().unwrap(), 1000.0);
    assert_eq!(spectra.global().max_intensity().unwrap(), 1000.0);
    let chromatograms = both.chromatogram_range_manager().unwrap();
    assert_eq!(chromatograms.has_range(), HasRangeType::All);
    assert_eq!(chromatograms.min_rt().unwrap(), 10.0);
    assert_eq!(chromatograms.max_rt().unwrap(), 20.0);
    assert_eq!(chromatograms.min_intensity().unwrap(), 500.0);
    assert_eq!(chromatograms.max_intensity().unwrap(), 1500.0);
    let combined = both.combined_range_manager().unwrap();
    assert_eq!(combined.has_range(), HasRangeType::Some);
    // The source annotates this `0` with "TODO: Why 0? precursor m/z not set?":
    // the chromatogram's unset `Product` m/z of zero extends the m/z dimension.
    assert_eq!(combined.min_mz().unwrap(), 0.0);
    assert_eq!(combined.max_mz().unwrap(), 100.0);
    assert_eq!(combined.min_rt().unwrap(), 10.0);
    assert_eq!(combined.max_rt().unwrap(), 30.0);
    assert_eq!(combined.min_intensity().unwrap(), 500.0);
    assert_eq!(combined.max_intensity().unwrap(), 1500.0);

    // The section's tail: three MS levels.
    let levels = MSExperiment {
        spectra: vec![
            plain_peak_scan(30., 1, 100.0, 1000.0),
            plain_peak_scan(35., 2, 200.0, 2000.0),
            plain_peak_scan(40., 3, 300.0, 3000.0),
        ],
        ..Default::default()
    };
    let spectra = levels.spectrum_range_manager().unwrap();
    for (level, mz) in [(1, 100.0), (2, 200.0), (3, 300.0)] {
        let manager = spectra.by_ms_level(level).unwrap();
        assert_eq!(manager.min_mz().unwrap(), mz);
        assert_eq!(manager.max_mz().unwrap(), mz);
    }
    assert_eq!(spectra.global().min_mz().unwrap(), 100.0);
    assert_eq!(spectra.global().max_mz().unwrap(), 300.0);
    assert_eq!(spectra.global().min_intensity().unwrap(), 1000.0);
    assert_eq!(spectra.global().max_intensity().unwrap(), 3000.0);
}
