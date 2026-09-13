// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs` against its class
//! test and an executed product-SDK oracle.
//!
//! Every `START_SECTION` of `FeatureFinderAlgorithmPickedHelperStructs_test.cpp`
//! (core bc9cc12) is transcribed with its literals and its comparison macros.
//! `TEST_EQUAL` is exact equality after converting the expected value to the
//! actual value's type. `TEST_REAL_SIMILAR` is ClassTest's `isRealSimilar`
//! with its default tolerances, absolute 1e-5 and ratio 1 + 1e-5. The source
//! test builds shared state section by section; each test here rebuilds that
//! state in the same order.
//!
//! The oracle replay rebuilds every input of
//! `../oracle/feature-finder-picked-helper-structs/driver.cpp`, formats each
//! result the way the driver prints it, and requires the C++ run's output row
//! for row, with floating-point values compared bit for bit on every platform.
//! That is deliberately stricter than a tolerance; the support document's
//! evidence section gives the reasons and says when it must be re-measured.

// The class-test literals are transcribed verbatim, including digits beyond
// the precision of their type, so the f32 values match the C++ literals exactly.
#![allow(clippy::excessive_precision)]

use openms::Error;
use openms::analysis::feature_finder_picked::helper_structs::{
    IsotopePattern, MassTrace, MassTraces, PatternPeak, Seed, TheoreticalIsotopePattern, TracePeak,
};
use openms::kernel::Point2D;

/// Output of `run.sh` in `../oracle/feature-finder-picked-helper-structs/`.
const ORACLE: &str = include_str!("data/feature_finder_picked_helper_structs_oracle.tsv");

/// `mt1` of the class test: retention time and `f32` intensity; m/z is 1000.
const MT1: [(f64, f32); 10] = [
    (677.1, 1.08268226589),
    (677.4, 1.58318959267),
    (677.7, 2.22429840363),
    (678.0, 3.00248879081),
    (678.3, 3.89401804768),
    (678.6, 4.8522452777),
    (678.9, 5.80919229659),
    (679.2, 6.68216169129),
    (679.5, 7.38493077109),
    (679.8, 7.84158938645),
];

/// `mt2` of the class test (`p2_4` to `p2_6`); m/z is 1001.
const MT2: [(f64, f32); 3] = [
    (678.0, 0.750622197703),
    (678.3, 0.97350451192),
    (678.6, 1.21306131943),
];

/// ClassTest `isRealSimilar` with the default `absdiff_max_allowed` 1e-5 and
/// `ratio_max_allowed` 1 + 1e-5 (ClassTest.cpp:35-38 and 364-470). The source
/// evaluates in `long double`; the transcribed cases are far from both limits.
fn is_real_similar(number_1: f64, number_2: f64) -> bool {
    const ABSDIFF_MAX_ALLOWED: f64 = 1e-5;
    const RATIO_MAX_ALLOWED: f64 = 1.0 + 1e-5;
    if number_1.is_nan() || number_2.is_nan() {
        return false;
    }
    let is_absdiff_small = (number_1 - number_2).abs() <= ABSDIFF_MAX_ALLOWED;
    if number_1 == 0.0 || number_2 == 0.0 {
        return (number_1 == 0.0 && number_2 == 0.0) || is_absdiff_small;
    }
    let mut ratio = number_1 / number_2;
    if ratio < 0.0 {
        return is_absdiff_small;
    }
    if ratio < 1.0 {
        ratio = 1.0 / ratio;
    }
    ratio <= RATIO_MAX_ALLOWED || is_absdiff_small
}

fn assert_real_similar(actual: f64, expected: f64, what: &str) {
    assert!(
        is_real_similar(actual, expected),
        "{what}: {actual} is not similar to {expected}"
    );
}

/// The class test's `mt1`; `p1_k` is the peak with identity `(1, k)`.
fn class_mt1() -> MassTrace {
    let mut trace = MassTrace {
        theoretical_int: 0.8,
        ..MassTrace::default()
    };
    for (k, (rt, intensity)) in MT1.into_iter().enumerate() {
        trace
            .peaks
            .push(TracePeak::new(1, k + 1, rt, 1000.0, intensity));
    }
    trace
}

/// The class test's `mt2`; `p2_k` is the peak with identity `(2, k)`.
fn class_mt2() -> MassTrace {
    let mut trace = MassTrace {
        theoretical_int: 0.2,
        ..MassTrace::default()
    };
    for (k, (rt, intensity)) in MT2.into_iter().enumerate() {
        trace
            .peaks
            .push(TracePeak::new(2, k + 4, rt, 1001.0, intensity));
    }
    trace
}

/// `mt` holding only `mt1`, as it is before `mt2` is added: the `mt1` section
/// `updateMaximum` has already run when the class test copies it in.
fn class_mt_with_mt1() -> MassTraces {
    let mut mt1 = class_mt1();
    mt1.update_maximum();
    let mut mt = MassTraces::new();
    mt.push(mt1);
    mt
}

/// `mt` holding `mt1` and `mt2`.
fn class_mt() -> MassTraces {
    let mut mt = class_mt_with_mt1();
    mt.push(class_mt2());
    mt
}

/// `mt` after the border cases that precede the `computeIntensityProfile`
/// section: a leading, a gap and a trailing peak on the second trace.
fn class_mt_with_border_cases() -> MassTraces {
    let mut mt = class_mt();
    mt[1]
        .peaks
        .insert(0, TracePeak::new(2, 0, 676.8, 1001.0, 0.286529652));
    mt[1]
        .peaks
        .push(TracePeak::new(2, 7, 679.2, 1001.0, 0.72952935));
    mt[1]
        .peaks
        .push(TracePeak::new(2, 8, 680.1, 1001.0, 0.672624672));
    mt
}

// ---------------------------------------------------------------------------
// The 14 START_SECTIONs, in source order.
// ---------------------------------------------------------------------------

/// Section `[IsotopePattern] IsotopePattern(Size size)`.
#[test]
fn isotope_pattern_constructor_sizes_every_vector() {
    let expected_size = 10;
    let pattern = IsotopePattern::new(expected_size).unwrap();
    assert_eq!(pattern.intensity.len(), expected_size);
    assert_eq!(pattern.mz_score.len(), expected_size);
    assert_eq!(pattern.peak.len(), expected_size);
    assert_eq!(pattern.spectrum.len(), expected_size);
    assert_eq!(pattern.theoretical_mz.len(), expected_size);
}

/// Section `[MassTrace] ConvexHull2D getConvexhull() const`.
#[test]
fn mass_trace_convex_hull_encloses_the_trace() {
    let mt1 = class_mt1();
    let ch = mt1.convex_hull().unwrap();
    let p1_10_mz = mt1.peaks[9].mz;
    assert!(ch.encloses(Point2D::new(679.8, p1_10_mz)).unwrap());
    assert!(!ch.encloses(Point2D::new(679.8, p1_10_mz + 1.0)).unwrap());
    assert!(!ch.encloses(Point2D::new(679.9, p1_10_mz)).unwrap());
}

/// Section `[MassTrace] void updateMaximum()`.
#[test]
fn mass_trace_update_maximum_finds_p1_10() {
    let mut mt1 = class_mt1();
    mt1.update_maximum();
    let max_peak = mt1.max_peak.unwrap();
    assert_eq!((max_peak.spectrum, max_peak.peak), (1, 10));
    assert_eq!(mt1.max_rt, 679.8);
}

/// Section `[MassTrace] double getAvgMZ() const`.
#[test]
fn mass_trace_avg_mz_is_intensity_weighted() {
    assert_eq!(class_mt1().avg_mz(), 1000.0);

    let mut mt_avg = MassTrace::default();
    mt_avg.peaks.push(TracePeak::new(0, 1, 100.0, 10.5, 1000.0));
    mt_avg.peaks.push(TracePeak::new(0, 2, 100.0, 10.0, 100.0));
    mt_avg.peaks.push(TracePeak::new(0, 3, 100.0, 9.5, 10.0));
    assert_real_similar(mt_avg.avg_mz(), 10.4459, "mt_avg.getAvgMZ()");
}

/// Section `[MassTrace] bool isValid() const`.
#[test]
fn mass_trace_is_valid_needs_three_peaks() {
    let mt1 = class_mt1();
    assert!(mt1.is_valid());
    let mut mt_non_valid = MassTrace::default();
    mt_non_valid.peaks.push(mt1.peaks[9]);
    assert!(!mt_non_valid.is_valid());
    mt_non_valid.peaks.push(mt1.peaks[8]);
    assert!(!mt_non_valid.is_valid());
    mt_non_valid.peaks.push(mt1.peaks[7]);
    assert!(mt_non_valid.is_valid());
}

/// Section `[MassTraces] MassTraces()`.
#[test]
fn mass_traces_constructor_sets_max_trace_zero() {
    assert_eq!(class_mt_with_mt1().max_trace, 0);
}

/// Section `[MassTraces] Size getPeakCount() const`.
#[test]
fn mass_traces_peak_count() {
    assert_eq!(class_mt_with_mt1().peak_count(), 10);
    assert_eq!(MassTraces::new().peak_count(), 0);
}

/// Section `[MassTraces] bool isValid(double seed_mz, double trace_tolerance)`.
#[test]
fn mass_traces_is_valid_needs_two_traces_and_the_seed() {
    let mut invalid_traces = MassTraces::new();
    invalid_traces.push(class_mt1());
    assert!(!invalid_traces.is_valid(600.0, 0.03));

    let mt = class_mt();
    assert!(mt.is_valid(1000.0, 0.00));
    assert!(mt.is_valid(1001.003, 0.03));
    assert!(!mt.is_valid(1002.0, 0.003));
}

/// Section `[MassTraces] Size getTheoreticalmaxPosition() const`.
#[test]
fn mass_traces_theoretical_max_position() {
    assert!(matches!(
        MassTraces::new().theoretical_max_position(),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(class_mt().theoretical_max_position().unwrap(), 0);
}

/// Section `[MassTraces] void updateBaseline()`.
#[test]
fn mass_traces_update_baseline_is_the_lowest_intensity() {
    let mut empty_traces = MassTraces::new();
    empty_traces.update_baseline();
    assert_eq!(empty_traces.baseline, 0.0);

    let mut mt = class_mt();
    mt.update_baseline();
    let p2_4 = MT2[0].1;
    assert_eq!(mt.baseline, f64::from(p2_4));
}

/// Section `[MassTraces] std::pair<double,double> getRTBounds() const`.
#[test]
fn mass_traces_rt_bounds() {
    assert!(matches!(
        MassTraces::new().rt_bounds(),
        Err(Error::InvalidValue(_))
    ));
    let bounds = class_mt().rt_bounds().unwrap();
    assert_eq!(bounds.0, 677.1);
    assert_eq!(bounds.1, 679.8);
}

/// Section `[MassTraces] void computeIntensityProfile(std::list<...>) const`.
#[test]
fn mass_traces_intensity_profile_merges_the_traces() {
    let profile = class_mt_with_border_cases().intensity_profile().unwrap();
    assert_eq!(profile.len(), 12);

    // (rt, expected intensity) with the source's f32 literal sums.
    let expected: [(f64, f32); 12] = [
        (676.8, 0.286529652),
        (677.1, 1.08268226589),
        (677.4, 1.58318959267),
        (677.7, 2.22429840363),
        (678.0, 3.00248879081 + 0.750622197703),
        (678.3, 3.89401804768 + 0.97350451192),
        (678.6, 4.8522452777 + 1.21306131943),
        (678.9, 5.80919229659),
        (679.2, 6.68216169129 + 0.72952935),
        (679.5, 7.38493077109),
        (679.8, 7.84158938645),
        (680.1, 0.672624672),
    ];
    for (index, (entry, (rt, intensity))) in profile.iter().zip(expected).enumerate() {
        assert_real_similar(entry.0, rt, &format!("profile[{index}].first"));
        assert_real_similar(
            entry.1,
            f64::from(intensity),
            &format!("profile[{index}].second"),
        );
    }
}

/// Section `[Seed] bool operator<(const Seed &rhs) const`.
#[test]
fn seed_orders_by_intensity() {
    let s1 = Seed {
        intensity: 100.0,
        ..Seed::default()
    };
    let s2 = Seed {
        intensity: 200.0,
        ..Seed::default()
    };
    let s3 = Seed {
        intensity: 300.0,
        ..Seed::default()
    };

    assert!(s1.is_less_intense_than(&s2));
    assert!(s1.is_less_intense_than(&s3));
    assert!(s2.is_less_intense_than(&s3));

    assert!(!s2.is_less_intense_than(&s1));
    assert!(!s3.is_less_intense_than(&s1));
    assert!(!s3.is_less_intense_than(&s2));
}

/// Section `[TheoreticalIsotopePattern] Size size() const`.
#[test]
fn theoretical_isotope_pattern_size() {
    let mut theo_pattern = TheoreticalIsotopePattern::default();
    assert_eq!(theo_pattern.len(), 0);
    theo_pattern.intensity.push(0.7);
    theo_pattern.intensity.push(0.2);
    theo_pattern.intensity.push(0.1);
    assert_eq!(theo_pattern.len(), 3);
}

// ---------------------------------------------------------------------------
// Tier 1: replay of the executed product-SDK driver.
// ---------------------------------------------------------------------------

/// Formats results as `driver.cpp` prints them.
#[derive(Default)]
struct Replay {
    lines: Vec<String>,
    next_identity: usize,
}

impl Replay {
    fn row(&mut self, case: &str, quantity: &str, kind: &str, value: String) {
        self.lines
            .push(format!("{case}\t{quantity}\t{kind}\t{value}"));
    }

    fn f64(&mut self, case: &str, quantity: &str, value: f64) {
        self.row(case, quantity, "f64", format!("{:016x}", value.to_bits()));
    }

    fn usize(&mut self, case: &str, quantity: &str, value: usize) {
        self.row(case, quantity, "usize", value.to_string());
    }

    fn isize(&mut self, case: &str, quantity: &str, value: i64) {
        self.row(case, quantity, "isize", value.to_string());
    }

    fn boolean(&mut self, case: &str, quantity: &str, value: bool) {
        self.row(case, quantity, "bool", value.to_string());
    }

    /// The driver prints `Precondition` for a caught `Exception::Precondition`,
    /// which this port reports as `Error::InvalidValue`.
    fn error(&mut self, case: &str, quantity: &str, error: &Error) {
        let name = match error {
            Error::InvalidValue(_) => "Precondition".to_owned(),
            other => format!("unexpected {other:?}"),
        };
        self.row(case, quantity, "error", name);
    }

    /// A peak with a fresh identity, as each `Peak1D` of the driver has a fresh
    /// address.
    fn peak(&mut self, rt: f64, mz: f64, intensity: f32) -> TracePeak {
        self.next_identity += 1;
        TracePeak::new(self.next_identity, 0, rt, mz, intensity)
    }

    fn trace(&mut self, peaks: &[(f64, f64, f32)]) -> MassTrace {
        let mut trace = MassTrace::default();
        for &(rt, mz, intensity) in peaks {
            let peak = self.peak(rt, mz, intensity);
            trace.peaks.push(peak);
        }
        trace
    }

    fn profile(&mut self, case: &str, traces: &MassTraces) {
        let profile = traces.intensity_profile().unwrap();
        self.usize(case, "profile.len", profile.len());
        for (index, (rt, intensity)) in profile.into_iter().enumerate() {
            self.f64(case, &format!("profile[{index}].rt"), rt);
            self.f64(case, &format!("profile[{index}].intensity"), intensity);
        }
    }

    fn hull(&mut self, case: &str, trace: &MassTrace) {
        let points = trace.convex_hull().unwrap().hull_points();
        self.usize(case, "hull.len", points.len());
        for (index, point) in points.into_iter().enumerate() {
            self.f64(case, &format!("hull[{index}].rt"), point.rt);
            self.f64(case, &format!("hull[{index}].mz"), point.mz);
        }
    }

    fn maximum(&mut self, case: &str, trace: &mut MassTrace) {
        trace.update_maximum();
        let index = trace.max_peak.and_then(|max| {
            trace
                .peaks
                .iter()
                .position(|peak| (peak.spectrum, peak.peak) == (max.spectrum, max.peak))
        });
        self.isize(
            case,
            "max_peak.index",
            index.map_or(-1, |index| index as i64),
        );
        if trace.max_peak.is_some() {
            self.f64(case, "max_rt", trace.max_rt);
        }
    }
}

/// `main` of `driver.cpp`, statement for statement.
fn replay() -> Vec<String> {
    let mut r = Replay::default();

    {
        let c = "isotope_pattern_new";
        let pattern = IsotopePattern::new(10).unwrap();
        r.usize(c, "peak.len", pattern.peak.len());
        r.usize(c, "spectrum.len", pattern.spectrum.len());
        r.usize(c, "intensity.len", pattern.intensity.len());
        r.usize(c, "mz_score.len", pattern.mz_score.len());
        r.usize(c, "theoretical_mz.len", pattern.theoretical_mz.len());
        let peaks_not_found = pattern.peak.iter().all(|p| *p == PatternPeak::NotFound);
        let zeros = (0..10).all(|i| {
            pattern.spectrum[i] == 0
                && pattern.intensity[i] == 0.0
                && pattern.mz_score[i] == 0.0
                && pattern.theoretical_mz[i] == 0.0
        });
        r.boolean(c, "peak.all_minus_one", peaks_not_found);
        r.boolean(c, "vectors.all_zero", zeros);
        r.usize(
            c,
            "theoretical_pattern.intensity.len",
            pattern.theoretical_pattern.intensity.len(),
        );
    }

    let mut mt1 = MassTrace {
        theoretical_int: 0.8,
        ..MassTrace::default()
    };
    let mut p1 = Vec::new();
    for (rt, intensity) in MT1 {
        let peak = r.peak(rt, 1000.0, intensity);
        p1.push(peak);
        mt1.peaks.push(peak);
    }

    {
        let c = "mass_trace_convex_hull";
        let ch = mt1.convex_hull().unwrap();
        let mz = p1[9].mz;
        r.boolean(
            c,
            "encloses(679.8,1000)",
            ch.encloses(Point2D::new(679.8, mz)).unwrap(),
        );
        r.boolean(
            c,
            "encloses(679.8,1001)",
            ch.encloses(Point2D::new(679.8, mz + 1.0)).unwrap(),
        );
        r.boolean(
            c,
            "encloses(679.9,1000)",
            ch.encloses(Point2D::new(679.9, mz)).unwrap(),
        );
        r.hull(c, &mt1);
    }

    r.maximum("mass_trace_update_maximum", &mut mt1);

    {
        let c = "mass_trace_avg_mz";
        r.f64(c, "mt1", mt1.avg_mz());
        let mt_avg = r.trace(&[
            (100.0, 10.5, 1000.0),
            (100.0, 10.0, 100.0),
            (100.0, 9.5, 10.0),
        ]);
        r.f64(c, "mt_avg", mt_avg.avg_mz());
    }

    {
        let c = "mass_trace_is_valid";
        r.boolean(c, "mt1", mt1.is_valid());
        let mut mt_non_valid = MassTrace::default();
        mt_non_valid.peaks.push(p1[9]);
        r.boolean(c, "one_peak", mt_non_valid.is_valid());
        mt_non_valid.peaks.push(p1[8]);
        r.boolean(c, "two_peaks", mt_non_valid.is_valid());
        mt_non_valid.peaks.push(p1[7]);
        r.boolean(c, "three_peaks", mt_non_valid.is_valid());
    }

    let mut mt = MassTraces::new();
    let mut empty_traces = MassTraces::new();
    mt.push(mt1.clone());

    r.usize("mass_traces_new", "max_trace", mt.max_trace);
    r.usize("mass_traces_peak_count", "mt", mt.peak_count());
    r.usize("mass_traces_peak_count", "empty", empty_traces.peak_count());

    let mut mt2 = MassTrace {
        theoretical_int: 0.2,
        ..MassTrace::default()
    };
    for (rt, intensity) in MT2 {
        let peak = r.peak(rt, 1001.0, intensity);
        mt2.peaks.push(peak);
    }
    mt.push(mt2);

    {
        let c = "mass_traces_is_valid";
        let mut invalid_traces = MassTraces::new();
        invalid_traces.push(mt1.clone());
        r.boolean(
            c,
            "one_trace(600,0.03)",
            invalid_traces.is_valid(600.0, 0.03),
        );
        r.boolean(c, "mt(1000,0)", mt.is_valid(1000.0, 0.00));
        r.boolean(c, "mt(1001.003,0.03)", mt.is_valid(1001.003, 0.03));
        r.boolean(c, "mt(1002,0.003)", mt.is_valid(1002.0, 0.003));
    }

    {
        let c = "mass_traces_theoretical_max_position";
        match empty_traces.theoretical_max_position() {
            Ok(position) => r.usize(c, "empty", position),
            Err(error) => r.error(c, "empty", &error),
        }
        r.usize(c, "mt", mt.theoretical_max_position().unwrap());
    }

    {
        let c = "mass_traces_update_baseline";
        empty_traces.update_baseline();
        r.f64(c, "empty", empty_traces.baseline);
        mt.update_baseline();
        r.f64(c, "mt", mt.baseline);
    }

    {
        let c = "mass_traces_rt_bounds";
        match empty_traces.rt_bounds() {
            Ok((min, max)) => {
                r.f64(c, "empty.min", min);
                r.f64(c, "empty.max", max);
            }
            Err(error) => r.error(c, "empty", &error),
        }
        let (min, max) = mt.rt_bounds().unwrap();
        r.f64(c, "mt.min", min);
        r.f64(c, "mt.max", max);
    }

    let p2_0 = r.peak(676.8, 1001.0, 0.286529652);
    mt[1].peaks.insert(0, p2_0);
    let p2_7 = r.peak(679.2, 1001.0, 0.72952935);
    mt[1].peaks.push(p2_7);
    let p2_8 = r.peak(680.1, 1001.0, 0.672624672);
    mt[1].peaks.push(p2_8);

    r.profile("mass_traces_intensity_profile", &mt);

    {
        let c = "seed_less";
        let s1 = Seed::new(0, 0, 100.0);
        let s2 = Seed::new(0, 0, 200.0);
        let s3 = Seed::new(0, 0, 300.0);
        r.boolean(c, "s1<s2", s1.is_less_intense_than(&s2));
        r.boolean(c, "s1<s3", s1.is_less_intense_than(&s3));
        r.boolean(c, "s2<s3", s2.is_less_intense_than(&s3));
        r.boolean(c, "s2<s1", s2.is_less_intense_than(&s1));
        r.boolean(c, "s3<s1", s3.is_less_intense_than(&s1));
        r.boolean(c, "s3<s2", s3.is_less_intense_than(&s2));
        r.boolean(c, "s1<s1", s1.is_less_intense_than(&s1));
        let n = Seed::new(0, 0, f32::NAN);
        r.boolean(c, "nan<s1", n.is_less_intense_than(&s1));
        r.boolean(c, "s1<nan", s1.is_less_intense_than(&n));
        let neg0 = Seed::new(0, 0, -0.0);
        let pos0 = Seed::new(0, 0, 0.0);
        r.boolean(c, "-0<+0", neg0.is_less_intense_than(&pos0));
    }

    {
        let c = "theoretical_isotope_pattern_size";
        let mut theo_pattern = TheoreticalIsotopePattern::default();
        r.usize(c, "empty", theo_pattern.len());
        theo_pattern.intensity.extend([0.7, 0.2, 0.1]);
        r.usize(c, "three", theo_pattern.len());
    }

    // ======== Boundary cases beyond the class test ========

    {
        let mut traces = MassTraces::new();
        let t = r.trace(&[(1.0, 500.0, 1.0), (2.0, 500.0, 2.0), (3.0, 500.0, 4.0)]);
        traces.push(t);
        let t = r.trace(&[
            (2.0, 501.0, 8.0),
            (1.0, 501.0, 16.0),
            (2.0, 501.0, 32.0),
            (2.0, 501.0, 64.0),
        ]);
        traces.push(t);
        let t = r.trace(&[
            (0.5, 502.0, 128.0),
            (3.0, 502.0, 256.0),
            (3.0, 502.0, 512.0),
            (4.0, 502.0, 1024.0),
        ]);
        traces.push(t);
        r.profile("boundary_profile_unsorted_duplicates", &traces);
    }

    {
        let mut traces = MassTraces::new();
        traces.push(MassTrace::default());
        let t = r.trace(&[(1.0, 500.0, 1.5), (2.0, 500.0, 2.5)]);
        traces.push(t);
        r.profile("boundary_profile_empty_first_trace", &traces);
    }

    {
        let mut traces = MassTraces::new();
        for i in 0..10 {
            let t = r.trace(&[(7.0, 500.0 + f64::from(i), 0.1)]);
            traces.push(t);
        }
        r.profile("boundary_profile_double_accumulation", &traces);
    }

    {
        let mut tie = r.trace(&[
            (1.0, 500.0, 5.0),
            (2.0, 500.0, 7.0),
            (3.0, 500.0, 7.0),
            (4.0, 500.0, 3.0),
        ]);
        r.maximum("boundary_update_maximum_tie", &mut tie);
        let mut first_nan = r.trace(&[(1.0, 500.0, f32::NAN), (2.0, 500.0, 1.0)]);
        r.maximum("boundary_update_maximum_first_nan", &mut first_nan);
        let mut later_nan =
            r.trace(&[(1.0, 500.0, 1.0), (2.0, 500.0, f32::NAN), (3.0, 500.0, 2.0)]);
        r.maximum("boundary_update_maximum_later_nan", &mut later_nan);
        let mut empty = MassTrace::default();
        empty.update_maximum();
        r.boolean(
            "boundary_update_maximum_empty",
            "max_peak.is_null",
            empty.max_peak.is_none(),
        );
    }

    {
        let mut ties = MassTraces::new();
        for value in [0.5, 0.7, 0.7] {
            ties.push(MassTrace {
                theoretical_int: value,
                ..MassTrace::default()
            });
        }
        r.usize(
            "boundary_theoretical_max_position",
            "tie",
            ties.theoretical_max_position().unwrap(),
        );
        let mut leading_nan = MassTraces::new();
        for value in [f64::NAN, 1.0] {
            leading_nan.push(MassTrace {
                theoretical_int: value,
                ..MassTrace::default()
            });
        }
        r.usize(
            "boundary_theoretical_max_position",
            "leading_nan",
            leading_nan.theoretical_max_position().unwrap(),
        );
    }

    {
        let mut later_nan = MassTraces::new();
        let t = r.trace(&[(1.0, 500.0, 3.0), (2.0, 500.0, f32::NAN)]);
        later_nan.push(t);
        let t = r.trace(&[(1.0, 501.0, 1.0)]);
        later_nan.push(t);
        later_nan.update_baseline();
        r.f64("boundary_update_baseline", "later_nan", later_nan.baseline);
        let mut leading_nan = MassTraces::new();
        let t = r.trace(&[(1.0, 500.0, f32::NAN), (2.0, 500.0, 1.0)]);
        leading_nan.push(t);
        leading_nan.update_baseline();
        r.boolean(
            "boundary_update_baseline",
            "leading_nan.is_nan",
            leading_nan.baseline.is_nan(),
        );
        let mut skips_empty = MassTraces::new();
        skips_empty.push(MassTrace::default());
        let t = r.trace(&[(1.0, 500.0, 0.25), (2.0, 500.0, 0.125)]);
        skips_empty.push(t);
        skips_empty.update_baseline();
        r.f64(
            "boundary_update_baseline",
            "skips_empty_trace",
            skips_empty.baseline,
        );
    }

    {
        let mut no_peaks = MassTraces::new();
        no_peaks.push(MassTrace::default());
        let (min, max) = no_peaks.rt_bounds().unwrap();
        r.f64("boundary_rt_bounds", "no_peaks.min", min);
        r.f64("boundary_rt_bounds", "no_peaks.max", max);
        let mut with_nan = MassTraces::new();
        let t = r.trace(&[(f64::NAN, 500.0, 1.0), (5.0, 500.0, 1.0), (2.0, 500.0, 1.0)]);
        with_nan.push(t);
        let (min, max) = with_nan.rt_bounds().unwrap();
        r.f64("boundary_rt_bounds", "with_nan.min", min);
        r.f64("boundary_rt_bounds", "with_nan.max", max);
    }

    {
        let empty = MassTrace::default();
        r.boolean("boundary_avg_mz", "empty.is_nan", empty.avg_mz().is_nan());
        let zero = r.trace(&[(1.0, 10.0, 0.0), (2.0, 20.0, 0.0)]);
        r.boolean(
            "boundary_avg_mz",
            "zero_intensity.is_nan",
            zero.avg_mz().is_nan(),
        );
        let inexact = r.trace(&[
            (1.0, 400.123456789, 0.3),
            (2.0, 400.123556789, 0.7),
            (3.0, 400.123656789, 0.1),
        ]);
        r.f64("boundary_avg_mz", "inexact", inexact.avg_mz());
    }

    {
        let mut traces = MassTraces::new();
        traces.push(MassTrace::default());
        let t = r.trace(&[(1.0, 500.25, 1.0), (2.0, 500.25, 1.0), (3.0, 500.25, 1.0)]);
        traces.push(t);
        r.boolean(
            "boundary_traces_is_valid",
            "nan_then_match(500,0.25)",
            traces.is_valid(500.0, 0.25),
        );
        r.boolean(
            "boundary_traces_is_valid",
            "nan_then_miss(500,0.2499)",
            traces.is_valid(500.0, 0.2499),
        );
    }

    {
        let t = r.trace(&[
            (2.0, 10.0, 1.0),
            (1.0, 11.0, 1.0),
            (2.0, 12.0, 1.0),
            (3.0, 9.5, 1.0),
        ]);
        r.hull("boundary_convex_hull_unsorted", &t);
    }

    r.lines
}

#[test]
fn oracle_replay_matches_the_product_sdk_output_row_for_row() {
    let expected: Vec<&str> = ORACLE.lines().collect();
    let actual = replay();
    for (index, (row, want)) in actual.iter().zip(&expected).enumerate() {
        assert_eq!(row, want, "oracle row {}", index + 1);
    }
    assert_eq!(actual.len(), expected.len(), "oracle row count");
}

// ---------------------------------------------------------------------------
// Tier 4: native boundaries with no C++ analogue.
// ---------------------------------------------------------------------------

#[test]
fn profile_of_no_traces_is_empty() {
    assert!(MassTraces::new().intensity_profile().unwrap().is_empty());
}

#[test]
fn profile_refuses_a_nan_retention_time_that_the_source_loops_on() {
    // A NaN in a later trace meets a profile entry: the source loop never ends.
    let mut traces = MassTraces::new();
    let mut first = MassTrace::default();
    first.peaks.push(TracePeak::new(0, 0, 1.0, 500.0, 1.0));
    traces.push(first);
    let mut second = MassTrace::default();
    second
        .peaks
        .push(TracePeak::new(1, 0, f64::NAN, 500.0, 1.0));
    traces.push(second);
    assert!(matches!(
        traces.intensity_profile(),
        Err(Error::InvalidValue(_))
    ));

    // A NaN copied from the first trace, never compared, passes through.
    let mut traces = MassTraces::new();
    let mut only = MassTrace::default();
    only.peaks.push(TracePeak::new(0, 0, f64::NAN, 500.0, 1.0));
    traces.push(only);
    let profile = traces.intensity_profile().unwrap();
    assert_eq!(profile.len(), 1);
    assert!(profile[0].0.is_nan());

    // A NaN appended past the profile's end, never compared, passes through.
    let mut traces = MassTraces::new();
    traces.push(MassTrace::default());
    let mut second = MassTrace::default();
    second
        .peaks
        .push(TracePeak::new(1, 0, f64::NAN, 500.0, 2.0));
    traces.push(second);
    let profile = traces.intensity_profile().unwrap();
    assert_eq!(profile.len(), 1);
    assert!(profile[0].0.is_nan());
    assert_eq!(profile[0].1, 2.0);
}

#[test]
fn profile_ceilings_are_checked_before_allocating() {
    // One trace one peak over the peak ceiling.
    let mut traces = MassTraces::new();
    traces.push(MassTrace {
        peaks: vec![TracePeak::default(); MassTraces::MAX_PEAKS + 1],
        ..MassTrace::default()
    });
    assert!(matches!(
        traces.intensity_profile(),
        Err(Error::InvalidValue(_))
    ));

    // Exactly the peak ceiling, spread over 250 traces of 4000 peaks: within
    // MAX_PEAKS, but the merge bound 4000 * (249 * 250 / 2 + 249) exceeds
    // MAX_PROFILE_STEPS.
    let mut traces = MassTraces::new();
    for _ in 0..250 {
        traces.push(MassTrace {
            peaks: vec![TracePeak::default(); 4000],
            ..MassTrace::default()
        });
    }
    assert_eq!(traces.peak_count(), MassTraces::MAX_PEAKS);
    assert!(matches!(
        traces.intensity_profile(),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn convex_hull_refuses_non_finite_coordinates_and_oversized_traces() {
    let mut trace = MassTrace::default();
    trace
        .peaks
        .push(TracePeak::new(0, 0, f64::INFINITY, 500.0, 1.0));
    assert!(matches!(trace.convex_hull(), Err(Error::InvalidValue(_))));

    let mut trace = MassTrace::default();
    trace.peaks.push(TracePeak::new(0, 0, 1.0, f64::NAN, 1.0));
    assert!(matches!(trace.convex_hull(), Err(Error::InvalidValue(_))));

    let trace = MassTrace {
        peaks: vec![TracePeak::new(0, 0, 1.0, 500.0, 1.0); MassTraces::MAX_PEAKS + 1],
        ..MassTrace::default()
    };
    assert!(matches!(trace.convex_hull(), Err(Error::InvalidValue(_))));

    assert!(MassTrace::default().convex_hull().unwrap().is_empty());
}

#[test]
fn isotope_pattern_size_ceiling_and_initial_values() {
    assert!(matches!(
        IsotopePattern::new(IsotopePattern::MAX_SIZE + 1),
        Err(Error::InvalidValue(_))
    ));
    let pattern = IsotopePattern::new(3).unwrap();
    assert_eq!(pattern.peak, vec![PatternPeak::NotFound; 3]);
    assert_eq!(pattern.spectrum, vec![0; 3]);
    assert_eq!(pattern.intensity, vec![0.0; 3]);
    assert_eq!(pattern.mz_score, vec![0.0; 3]);
    assert_eq!(pattern.theoretical_mz, vec![0.0; 3]);
    assert_eq!(
        pattern.theoretical_pattern,
        TheoreticalIsotopePattern::default()
    );
    assert_eq!(IsotopePattern::new(0).unwrap(), IsotopePattern::default());

    assert_eq!(PatternPeak::Found(7).index(), Some(7));
    assert_eq!(PatternPeak::NotFound.index(), None);
    assert_eq!(PatternPeak::Removed.index(), None);
}

#[test]
fn reserve_is_bounded_and_leaves_the_collection_unchanged() {
    let mut traces = class_mt();
    traces.max_trace = 1;
    traces.baseline = 0.5;
    let before = traces.clone();
    assert!(matches!(
        traces.reserve(MassTraces::MAX_TRACES + 1),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(traces, before);
    traces.reserve(1).unwrap();
    traces.reserve(64).unwrap();
    assert_eq!(traces, before);
}

#[test]
fn clear_keeps_max_trace_and_baseline() {
    let mut traces = class_mt();
    traces.max_trace = 1;
    traces.baseline = 0.25;
    traces.clear();
    assert!(traces.is_empty());
    assert_eq!(traces.len(), 0);
    assert_eq!(traces.max_trace, 1);
    assert_eq!(traces.baseline, 0.25);
}

#[test]
fn collection_accessors_mirror_the_exported_vector_methods() {
    let mut traces = class_mt();
    assert_eq!(traces.len(), 2);
    assert!(traces.get(2).is_none());
    assert_eq!(traces.get(1).unwrap().theoretical_int, 0.2);
    assert_eq!(traces.last().unwrap().theoretical_int, 0.2);
    traces.last_mut().unwrap().theoretical_int = 0.3;
    traces.get_mut(0).unwrap().theoretical_int = 0.7;
    assert_eq!(traces[1].theoretical_int, 0.3);
    assert_eq!(traces.as_slice()[0].theoretical_int, 0.7);
    for trace in &mut traces {
        trace.theoretical_int *= 2.0;
    }
    let values: Vec<f64> = (&traces).into_iter().map(|t| t.theoretical_int).collect();
    assert_eq!(values, vec![1.4, 0.6]);
    assert_eq!(traces.iter_mut().count(), 2);
    assert!(MassTraces::new().last().is_none());
}

#[test]
fn seed_ordering_ignores_the_position() {
    let early = Seed::new(0, 0, 10.0);
    let late = Seed::new(99, 7, 10.0);
    assert!(!early.is_less_intense_than(&late));
    assert!(!late.is_less_intense_than(&early));
    assert_ne!(early, late);
}

#[test]
fn update_maximum_and_avg_mz_on_an_empty_trace() {
    let mut trace = MassTrace {
        max_rt: 3.0,
        ..MassTrace::default()
    };
    trace.update_maximum();
    assert!(trace.max_peak.is_none());
    assert_eq!(trace.max_rt, 3.0);
    assert!(!trace.is_valid());
}

#[test]
fn update_baseline_keeps_its_value_when_no_trace_has_peaks() {
    let mut traces = MassTraces::new();
    traces.push(MassTrace::default());
    traces.baseline = 42.0;
    traces.update_baseline();
    assert_eq!(traces.baseline, 42.0);
}
