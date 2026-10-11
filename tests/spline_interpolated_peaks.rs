// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Section-by-section ports of `SplinePackage_test.cpp` and
//! `SplineInterpolatedPeaks_test.cpp` at source revision
//! bc9cc12514c768385ce121d6ca4bb710fe1983c4
//! (`src/tests/class_tests/openms/source/`).
//!
//! Evidence tier 3: every expected value is a transcribed class-test literal,
//! compared as the class test compares it (`TEST_EQUAL` exactly,
//! `TEST_REAL_SIMILAR` at relative 1e-5). No expected value was derived from
//! this port's output.

use openms::kernel::{ChromatogramPeak, MSChromatogram, MSSpectrum, Peak1D};
use openms::processing::spline::spline_interpolated_peaks::{
    DEFAULT_NAVIGATOR_SCALING, SplineInterpolatedPeaks,
};
use openms::processing::spline::spline_package::SplinePackage;

fn similar(actual: f64, expected: f64) {
    let tolerance = 1e-5 * expected.abs().max(actual.abs());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} is not similar to {expected}"
    );
}

#[test]
fn spline_package_class_test() {
    let mz = [413.8, 413.9, 414.0, 414.1, 414.2];
    let intensity = [0.0, 100.2, 20.3, 2000.4, 4.3];
    let sp1 = SplinePackage::new(&mz, &intensity).unwrap();
    assert_eq!(sp1.pos_min(), 413.8);
    assert_eq!(sp1.pos_max(), 414.2);
    similar(sp1.pos_step_width(), 0.1);
    assert!(sp1.is_in_package(414.05));
    similar(sp1.eval(414.05).unwrap(), 1134.08593750018);
    // a single point is refused (Exception::IllegalArgument)
    assert!(SplinePackage::new(&[413.9], &[100.2]).is_err());
    // two points interpolate linearly
    let sp4 = SplinePackage::new(&[413.8, 413.9], &[0.0, 100.2]).unwrap();
    similar(sp4.eval(413.85).unwrap(), 50.1);
}

fn gauss(x: f64, centre: f64) -> f64 {
    (-(x - centre).powi(2) / (2.0 * 0.15 * 0.15)).exp()
}

fn two_gaussians() -> (Vec<f64>, Vec<f64>) {
    let mut pos = Vec::new();
    let mut intensity = Vec::new();
    for i in 0..11 {
        let x = 416.3 + 0.1 * f64::from(i);
        pos.push(x);
        intensity.push(gauss(x, 416.8));
    }
    for i in 0..11 {
        let x = 418.2 + 0.1 * f64::from(i);
        pos.push(x);
        intensity.push(gauss(x, 418.7));
    }
    (pos, intensity)
}

#[test]
fn spline_interpolated_peaks_class_test() {
    let (pos, intensity) = two_gaussians();

    // constructors from vectors, spectrum and chromatogram
    similar(
        SplineInterpolatedPeaks::new(&pos, &intensity)
            .unwrap()
            .pos_min(),
        416.3,
    );
    let mut spectrum = MSSpectrum {
        rt: 1789.0714,
        ..MSSpectrum::default()
    };
    let mut chromatogram = MSChromatogram::default();
    for (&p, &i) in pos.iter().zip(&intensity) {
        spectrum.peaks.push(Peak1D::new(p, i as f32));
        chromatogram.peaks.push(ChromatogramPeak::new(p, i as f32));
    }
    similar(
        SplineInterpolatedPeaks::from_spectrum(&spectrum)
            .unwrap()
            .pos_min(),
        416.3,
    );
    similar(
        SplineInterpolatedPeaks::from_chromatogram(&chromatogram)
            .unwrap()
            .pos_min(),
        416.3,
    );

    let spectrum2 = SplineInterpolatedPeaks::new(&pos, &intensity).unwrap();
    assert_eq!(spectrum2.pos_min(), 416.3);
    assert_eq!(spectrum2.pos_max(), 419.2);
    assert_eq!(spectrum2.size(), 2);
    let nav = || spectrum2.navigator(DEFAULT_NAVIGATOR_SCALING).unwrap();

    // Navigator::eval outside the Gaussians
    assert_eq!(nav().eval(400.0).unwrap(), 0.0);
    assert_eq!(nav().eval(417.8).unwrap(), 0.0);
    assert_eq!(nav().eval(500.0).unwrap(), 0.0);
    // near the edge and near the maximum
    similar(nav().eval(416.33).unwrap(), 0.007848195698809);
    similar(nav().eval(416.81).unwrap(), 0.997572728799559);
    // first package, then a search in the last
    let mut n = nav();
    similar(n.eval(416.81).unwrap(), 0.997572728799559);
    similar(n.eval(418.75).unwrap(), 0.944147611428987);
    // last package, then a search in the first
    let mut n2 = nav();
    similar(n2.eval(418.75).unwrap(), 0.944147611428987);
    similar(n2.eval(416.81).unwrap(), 0.997572728799559);

    // Navigator::getNextPos: within a package, to the next one, beyond range
    assert_eq!(nav().next_pos(417.0).unwrap(), 417.07);
    assert_eq!(nav().next_pos(417.29).unwrap(), 418.2);
    similar(nav().next_pos(500.0).unwrap(), 419.2);
}

#[test]
fn spline_interpolated_peaks_single_point_joins_the_next_package() {
    let mut pos3 = Vec::new();
    let mut intensity3 = Vec::new();
    for i in 0..4 {
        pos3.push(400.0 + f64::from(i) * 0.5);
        intensity3.push(10.0);
    }
    pos3.push(407.5);
    intensity3.push(10.0);
    for i in 0..4 {
        pos3.push(410.0 + f64::from(i) * 0.5);
        intensity3.push(10.0);
    }
    let spectrum3 = SplineInterpolatedPeaks::new(&pos3, &intensity3).unwrap();
    assert_eq!(spectrum3.size(), 2);
    let nav = || spectrum3.navigator(DEFAULT_NAVIGATOR_SCALING).unwrap();
    assert_eq!(nav().eval(405.0).unwrap(), 0.0);
    assert_eq!(nav().eval(408.0).unwrap(), 10.0);

    // a single point is refused (Exception::IllegalArgument)
    assert!(SplineInterpolatedPeaks::new(&[407.5], &[10.0]).is_err());
}

/// Tier 4: `getNavigator` refuses a spectrum without packages
/// (`Exception::InvalidSize`); an all-zero spectrum has none here, where the
/// source reads past its work vector.
#[test]
fn spline_interpolated_peaks_without_packages() {
    let zeros = SplineInterpolatedPeaks::new(&[1.0, 2.0, 3.0], &[0.0, 0.0, 0.0]).unwrap();
    assert_eq!(zeros.size(), 0);
    assert!(zeros.navigator(DEFAULT_NAVIGATOR_SCALING).is_err());
}
