// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Class-test coverage for the spectrum-similarity functor hierarchy:
//! `COMPARISON/PeakSpectrumCompareFunctor.h`,
//! `COMPARISON/BinnedSpectrumCompareFunctor.h`,
//! `COMPARISON/BinnedSharedPeakCount.h`,
//! `COMPARISON/BinnedSpectralContrastAngle.h` and
//! `COMPARISON/BinnedSumAgreeingIntensities.h`.

use openms::comparison::{
    BinConfig, BinUnit, BinnedSharedPeakCount, BinnedSpectralContrastAngle, BinnedSpectrum,
    BinnedSpectrumCompareFunctor, BinnedSumAgreeingIntensities, PeakSpectrumCompareFunctor,
};
use openms::format::dta;
use openms::param::DefaultParamHandler;
use openms::{Error, MSSpectrum, Peak1D, Result};

/// `PILISSequenceDB_DFPIANGER_1.dta`, the fixture every one of the three scorer
/// class tests loads, byte-identical to the upstream copy.
fn golden() -> MSSpectrum {
    dta::read(include_bytes!("data/comparison_dfpianger.dta").as_slice()).unwrap()
}

/// `BinnedSpectrum bs1(s1, 1.5, false, 2, offset)`, the upstream construction.
fn binned(spectrum: &MSSpectrum, offset: f32) -> BinnedSpectrum {
    BinnedSpectrum::new(
        spectrum,
        BinConfig {
            size: 1.5,
            unit: BinUnit::Absolute,
            spread: 2,
            offset,
            ..Default::default()
        },
    )
    .unwrap()
}

fn peaks(mzs: &[f64], intensities: &[f32]) -> MSSpectrum {
    MSSpectrum::from_peaks(
        mzs.iter()
            .zip(intensities)
            .map(|(&mz, &intensity)| Peak1D::new(mz, intensity))
            .collect(),
    )
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-5,
        "{actual} is not TEST_REAL_SIMILAR to {expected}"
    );
}

// ---------------------------------------------------------------------------
// COMPARISON/PeakSpectrumCompareFunctor.h
//
// Every one of its six class-test sections is NOT_TESTABLE upstream, because
// the class is a pure interface. The Rust analogue of "a derived class can be
// written and called through the base" is a trait implementation exercised
// through a trait object, which is what these tests do.
// ---------------------------------------------------------------------------

/// Minimal derived functor: the sum over matched m/z of the intensity products.
/// It exists only to exercise the base trait, as a derived class would.
struct SharedIntensityProduct {
    handler: DefaultParamHandler,
}

impl SharedIntensityProduct {
    fn new() -> Result<Self> {
        let mut handler = DefaultParamHandler::new("PeakSpectrumCompareFunctor")?;
        handler.set_name("SharedIntensityProduct")?;
        handler.defaults_to_parameters()?;
        Ok(Self { handler })
    }
}

impl PeakSpectrumCompareFunctor for SharedIntensityProduct {
    fn handler(&self) -> &DefaultParamHandler {
        &self.handler
    }
    fn handler_mut(&mut self) -> &mut DefaultParamHandler {
        &mut self.handler
    }
    fn score(&self, a: &MSSpectrum, b: &MSSpectrum) -> Result<f64> {
        let mut total = 0.0;
        for left in &a.peaks {
            for right in &b.peaks {
                if left.mz == right.mz {
                    total += f64::from(left.intensity) * f64::from(right.intensity);
                }
            }
        }
        Ok(total)
    }
}

#[test]
fn peak_spectrum_compare_functor_base_is_a_named_handler_and_a_pair_of_operators() {
    // START_SECTION(PeakSpectrumCompareFunctor())
    // START_SECTION(~PeakSpectrumCompareFunctor())
    // The base names its handler after itself and the derived functor renames
    // it, so the observable result of construction is the derived name with an
    // empty parameter tree. Destruction is drop glue with no observable effect.
    let functor = SharedIntensityProduct::new().unwrap();
    assert_eq!(functor.name(), "SharedIntensityProduct");
    assert_eq!(functor.handler().name(), "SharedIntensityProduct");
    assert!(functor.handler().parameters().is_empty());
    assert!(functor.handler().defaults().is_empty());

    // START_SECTION(PeakSpectrumCompareFunctor(const PeakSpectrumCompareFunctor& source))
    // START_SECTION(PeakSpectrumCompareFunctor& operator=(const PeakSpectrumCompareFunctor& source))
    // The C++ copy constructor and assignment operator both forward to
    // DefaultParamHandler, so name and parameters survive; here that is Clone.
    let copy = SharedIntensityProduct {
        handler: functor.handler().clone(),
    };
    assert_eq!(copy.name(), functor.name());
    assert_eq!(copy.handler().parameters(), functor.handler().parameters());

    // START_SECTION(double operator () (const PeakSpectrum& a, const PeakSpectrum& b) const)
    let a = peaks(&[100.0, 200.0, 300.0], &[1.0, 2.0, 4.0]);
    let b = peaks(&[200.0, 300.0], &[3.0, 5.0]);
    assert_eq!(functor.score(&a, &b).unwrap(), 2.0 * 3.0 + 4.0 * 5.0);

    // START_SECTION(double operator () (const PeakSpectrum& a) const)
    // The default self_score is the delegation every source derivative writes.
    assert_eq!(
        functor.self_score(&a).unwrap(),
        functor.score(&a, &a).unwrap()
    );
    assert_eq!(functor.self_score(&a).unwrap(), 1.0 + 4.0 + 16.0);

    // The base exists for polymorphic use; the trait is object safe.
    let erased: &dyn PeakSpectrumCompareFunctor = &functor;
    assert_eq!(erased.name(), "SharedIntensityProduct");
    assert_eq!(erased.self_score(&a).unwrap(), 21.0);
}

#[test]
fn peak_spectrum_compare_functor_handler_is_mutable_as_in_the_source() {
    // DefaultParamHandler::setName is public on the source base, so a functor's
    // name is not frozen at construction; handler_mut is that access.
    let mut functor = SharedIntensityProduct::new().unwrap();
    functor.handler_mut().set_name("renamed").unwrap();
    assert_eq!(functor.name(), "renamed");
}

// ---------------------------------------------------------------------------
// COMPARISON/BinnedSpectrumCompareFunctor.h
//
// All eight of its class-test sections are NOT_TESTABLE upstream. Two of them
// name a nested BinnedSpectrumCompareFunctor::IncompatibleBinning exception
// class that no longer exists in the pinned header; its role is taken by
// Exception::IllegalArgument thrown from the derived functors, which is what
// the incompatible-binning test below covers.
// ---------------------------------------------------------------------------

#[test]
fn binned_spectrum_compare_functor_base_is_a_named_handler_and_a_pair_of_operators() {
    // START_SECTION(BinnedSpectrumCompareFunctor())
    // START_SECTION(~BinnedSpectrumCompareFunctor())
    // START_SECTION((BinnedSpectrumCompareFunctor(const BinnedSpectrumCompareFunctor &source)))
    // START_SECTION((BinnedSpectrumCompareFunctor& operator=(const BinnedSpectrumCompareFunctor &source)))
    // Construction, copy and assignment are the base's DefaultParamHandler
    // behaviour, observable through every derivative.
    let shared = BinnedSharedPeakCount::new().unwrap();
    let angle = BinnedSpectralContrastAngle::new().unwrap();
    let agreeing = BinnedSumAgreeingIntensities::new().unwrap();
    let named: [(&dyn BinnedSpectrumCompareFunctor, &str); 3] = [
        (&shared, "BinnedSharedPeakCount"),
        (&angle, "BinnedSpectralContrastAngle"),
        (&agreeing, "BinnedSumAgreeingIntensities"),
    ];
    for (functor, name) in named {
        assert_eq!(functor.name(), name);
        assert!(functor.handler().parameters().is_empty());
        assert!(functor.handler().defaults().is_empty());
    }

    // START_SECTION((virtual double operator()(const BinnedSpectrum &spec1, const BinnedSpectrum &spec2) const =0))
    // START_SECTION((virtual double operator()(const BinnedSpectrum &spec) const =0))
    // Both pure virtuals are reachable through the base, and the one-spectrum
    // overload is the delegation every derivative implements.
    let spectrum = binned(&golden(), 0.4);
    let functors: [&dyn BinnedSpectrumCompareFunctor; 3] = [&shared, &angle, &agreeing];
    for functor in functors {
        let pairwise = functor.score(&spectrum, &spectrum).unwrap();
        assert_eq!(functor.self_score(&spectrum).unwrap(), pairwise);
        // Each of the three is normalised, and self-similarity is exactly one:
        // shared bins over the larger bin count is len/len; the contrast angle
        // is sum1 / sqrt(sum1 * sum1); and every bin agrees with itself, so the
        // agreeing sum equals the mean total intensity. None of these depends on
        // a transcribed literal.
        assert_eq!(pairwise, 1.0);
    }
}

#[test]
fn binned_functors_reject_incompatible_binning() {
    // Replaces the two class-test sections naming the removed nested
    // BinnedSpectrumCompareFunctor::IncompatibleBinning exception. Upstream,
    // BinnedSharedPeakCount and BinnedSumAgreeingIntensities throw
    // Exception::IllegalArgument here while BinnedSpectralContrastAngle only
    // asserts through OPENMS_PRECONDITION, which is a no-op in release builds;
    // this port makes all three an Error::InvalidValue.
    let spectrum = golden();
    let reference = binned(&spectrum, 0.0);
    let different_size = BinnedSpectrum::new(
        &spectrum,
        BinConfig {
            size: 2.0,
            unit: BinUnit::Absolute,
            spread: 2,
            offset: 0.0,
            ..Default::default()
        },
    )
    .unwrap();
    let different_offset = binned(&spectrum, 0.4);

    let shared = BinnedSharedPeakCount::new().unwrap();
    let angle = BinnedSpectralContrastAngle::new().unwrap();
    let agreeing = BinnedSumAgreeingIntensities::new().unwrap();
    let functors: [&dyn BinnedSpectrumCompareFunctor; 3] = [&shared, &angle, &agreeing];
    for functor in functors {
        for other in [&different_size, &different_offset] {
            assert!(matches!(
                functor.score(&reference, other),
                Err(Error::InvalidValue(_))
            ));
            assert!(matches!(
                functor.score(other, &reference),
                Err(Error::InvalidValue(_))
            ));
        }
    }

    // Spread is deliberately not part of compatibility, in the source and here.
    let wider_spread = BinnedSpectrum::new(
        &spectrum,
        BinConfig {
            size: 1.5,
            unit: BinUnit::Absolute,
            spread: 3,
            offset: 0.0,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(shared.score(&reference, &wider_spread).is_ok());
}

// ---------------------------------------------------------------------------
// COMPARISON/BinnedSharedPeakCount.h
// ---------------------------------------------------------------------------

#[test]
fn binned_shared_peak_count_constructs_copies_and_assigns() {
    // START_SECTION(BinnedSharedPeakCount())
    // START_SECTION(~BinnedSharedPeakCount())
    let functor = BinnedSharedPeakCount::new().unwrap();
    assert_eq!(functor.name(), "BinnedSharedPeakCount");
    assert!(functor.handler().parameters().is_empty());

    // START_SECTION((BinnedSharedPeakCount(const BinnedSharedPeakCount &source)))
    let copy = functor.clone();
    assert_eq!(copy.name(), functor.name());
    assert_eq!(copy.handler().parameters(), functor.handler().parameters());

    // START_SECTION((BinnedSharedPeakCount& operator=(const BinnedSharedPeakCount &source)))
    let mut assigned = BinnedSharedPeakCount::new().unwrap();
    assigned.handler_mut().set_name("scratch").unwrap();
    assert_eq!(assigned.name(), "scratch");
    assigned = functor.clone();
    assert_eq!(assigned.name(), functor.name());
    assert_eq!(
        assigned.handler().parameters(),
        functor.handler().parameters()
    );
}

#[test]
fn binned_shared_peak_count_scores_the_upstream_fixture() {
    // START_SECTION((double operator()(const BinnedSpectrum &spec1, const BinnedSpectrum &spec2) const))
    let s1 = golden();
    let mut s2 = s1.clone();
    s2.peaks.pop(); // s2.pop_back()
    let bs1 = binned(&s1, 0.0);
    let bs2 = binned(&s2, 0.0);
    close(
        BinnedSharedPeakCount::new()
            .unwrap()
            .score(&bs1, &bs2)
            .unwrap(),
        0.997118,
    );
    assert_eq!(
        BinnedSharedPeakCount::new()
            .unwrap()
            .score(&bs1, &bs1)
            .unwrap(),
        1.0
    );

    // Independently derived: dropping the last peak removes one of the 347
    // stored bins and leaves the rest shared, so the score is 346/347. The
    // upstream literal 0.997118 is that ratio to six places.
    assert_eq!(bs1.bins().len(), 347);
    assert_eq!(bs2.bins().len(), 346);
    assert_eq!(
        BinnedSharedPeakCount::new()
            .unwrap()
            .score(&bs1, &bs2)
            .unwrap(),
        346.0 / 347.0
    );
}

#[test]
fn binned_shared_peak_count_self_similarity_is_one() {
    // START_SECTION((double operator()(const BinnedSpectrum &spec) const ))
    let bs1 = binned(&golden(), 0.4);
    close(
        BinnedSharedPeakCount::new()
            .unwrap()
            .self_score(&bs1)
            .unwrap(),
        1.0,
    );
}

#[test]
fn binned_shared_peak_count_counts_stored_bins_not_nonzero_values() {
    // Eigen counts stored coefficients, so a bin holding an explicit zero is
    // shared. Two spectra that store nothing divide by zero upstream and return
    // NaN; here the defined score is zero.
    let layout = BinConfig {
        size: 1.0,
        unit: BinUnit::Absolute,
        spread: 0,
        offset: 0.0,
        ..Default::default()
    };
    let zero = BinnedSpectrum::new(&peaks(&[1.0], &[0.0]), layout).unwrap();
    let empty = BinnedSpectrum::new(&MSSpectrum::default(), layout).unwrap();
    let functor = BinnedSharedPeakCount::new().unwrap();
    assert_eq!(zero.bins().len(), 1);
    assert_eq!(functor.score(&zero, &zero).unwrap(), 1.0);
    assert_eq!(functor.score(&empty, &empty).unwrap(), 0.0);
    assert_eq!(functor.score(&zero, &empty).unwrap(), 0.0);

    // Disjoint bins share nothing; the denominator is the larger bin count.
    let left = BinnedSpectrum::new(&peaks(&[1.0, 2.0, 3.0], &[1.0, 1.0, 1.0]), layout).unwrap();
    let right = BinnedSpectrum::new(&peaks(&[10.0, 11.0], &[1.0, 1.0]), layout).unwrap();
    assert_eq!(functor.score(&left, &right).unwrap(), 0.0);
    let overlap = BinnedSpectrum::new(&peaks(&[2.0, 10.0], &[1.0, 1.0]), layout).unwrap();
    assert_eq!(functor.score(&left, &overlap).unwrap(), 1.0 / 3.0);
    assert_eq!(functor.score(&overlap, &left).unwrap(), 1.0 / 3.0);
}

// ---------------------------------------------------------------------------
// COMPARISON/BinnedSpectralContrastAngle.h
// ---------------------------------------------------------------------------

#[test]
fn binned_spectral_contrast_angle_constructs_copies_and_assigns() {
    // START_SECTION(BinnedSpectralContrastAngle())
    // START_SECTION(~BinnedSpectralContrastAngle())
    let functor = BinnedSpectralContrastAngle::new().unwrap();
    assert_eq!(functor.name(), "BinnedSpectralContrastAngle");
    assert!(functor.handler().parameters().is_empty());

    // START_SECTION((BinnedSpectralContrastAngle(const BinnedSpectralContrastAngle &source)))
    let copy = functor.clone();
    assert_eq!(copy.name(), functor.name());
    assert_eq!(copy.handler().parameters(), functor.handler().parameters());

    // START_SECTION((BinnedSpectralContrastAngle& operator=(const BinnedSpectralContrastAngle &source)))
    let mut assigned = BinnedSpectralContrastAngle::new().unwrap();
    assigned.handler_mut().set_name("scratch").unwrap();
    assert_eq!(assigned.name(), "scratch");
    assigned = functor.clone();
    assert_eq!(assigned.name(), functor.name());
    assert_eq!(
        assigned.handler().parameters(),
        functor.handler().parameters()
    );
}

#[test]
fn binned_spectral_contrast_angle_scores_the_upstream_fixture() {
    // START_SECTION((double operator()(const BinnedSpectrum &spec1, const BinnedSpectrum &spec2) const))
    let s1 = golden();
    let mut s2 = s1.clone();
    s2.peaks.pop();
    let bs1 = binned(&s1, 0.4);
    let bs2 = binned(&s2, 0.4);
    let functor = BinnedSpectralContrastAngle::new().unwrap();
    close(functor.score(&bs1, &bs2).unwrap(), 0.999985);

    // The upstream "empty / all-zero spectrum must yield a defined score of 0"
    // regression case.
    let empty = binned(&MSSpectrum::default(), 0.4);
    assert_eq!(functor.score(&empty, &bs1).unwrap(), 0.0);
    assert_eq!(functor.score(&bs1, &empty).unwrap(), 0.0);
    assert_eq!(functor.score(&empty, &empty).unwrap(), 0.0);
}

#[test]
fn binned_spectral_contrast_angle_self_similarity_is_one() {
    // START_SECTION((double operator()(const BinnedSpectrum &spec) const ))
    let bs1 = binned(&golden(), 0.4);
    let functor = BinnedSpectralContrastAngle::new().unwrap();
    close(functor.self_score(&bs1).unwrap(), 1.0);
    // Derived, not transcribed: sum1 / sqrt(sum1 * sum1) is exactly one because
    // the square of an f32-valued f64 is exact in f64.
    assert_eq!(functor.self_score(&bs1).unwrap(), 1.0);
}

#[test]
fn binned_spectral_contrast_angle_is_a_cosine_over_the_bin_vectors() {
    // Derived values: orthogonal bin vectors score 0, antiparallel ones -1, and
    // the source neither clamps nor converts the cosine into an angle.
    let layout = BinConfig {
        size: 1.0,
        unit: BinUnit::Absolute,
        spread: 0,
        offset: 0.0,
        ..Default::default()
    };
    let functor = BinnedSpectralContrastAngle::new().unwrap();
    let a = BinnedSpectrum::new(&peaks(&[1.0, 2.0], &[3.0, 4.0]), layout).unwrap();
    let orthogonal = BinnedSpectrum::new(&peaks(&[10.0], &[5.0]), layout).unwrap();
    assert_eq!(functor.score(&a, &orthogonal).unwrap(), 0.0);

    // cos = (3*4 + 4*3) / sqrt(25 * 25) = 24/25.
    let swapped = BinnedSpectrum::new(&peaks(&[1.0, 2.0], &[4.0, 3.0]), layout).unwrap();
    assert_eq!(functor.score(&a, &swapped).unwrap(), 24.0 / 25.0);
    assert_eq!(functor.score(&swapped, &a).unwrap(), 24.0 / 25.0);

    // Negative bins are permitted upstream and the score is not clamped.
    let opposed = BinnedSpectrum::new(&peaks(&[1.0, 2.0], &[-3.0, -4.0]), layout).unwrap();
    assert_eq!(functor.score(&a, &opposed).unwrap(), -1.0);
}

// ---------------------------------------------------------------------------
// COMPARISON/BinnedSumAgreeingIntensities.h
// ---------------------------------------------------------------------------

#[test]
fn binned_sum_agreeing_intensities_constructs_copies_and_assigns() {
    // START_SECTION(BinnedSumAgreeingIntensities())
    // START_SECTION(~BinnedSumAgreeingIntensities())
    let functor = BinnedSumAgreeingIntensities::new().unwrap();
    assert_eq!(functor.name(), "BinnedSumAgreeingIntensities");
    assert!(functor.handler().parameters().is_empty());

    // START_SECTION((BinnedSumAgreeingIntensities(const BinnedSumAgreeingIntensities &source)))
    let copy = functor.clone();
    assert_eq!(copy.name(), functor.name());
    assert_eq!(copy.handler().parameters(), functor.handler().parameters());

    // START_SECTION((BinnedSumAgreeingIntensities& operator=(const BinnedSumAgreeingIntensities &source)))
    let mut assigned = BinnedSumAgreeingIntensities::new().unwrap();
    assigned.handler_mut().set_name("scratch").unwrap();
    assert_eq!(assigned.name(), "scratch");
    assigned = functor.clone();
    assert_eq!(assigned.name(), functor.name());
    assert_eq!(
        assigned.handler().parameters(),
        functor.handler().parameters()
    );
}

#[test]
fn binned_sum_agreeing_intensities_scores_the_upstream_fixture() {
    // START_SECTION((double operator()(const BinnedSpectrum &spec1, const BinnedSpectrum &spec2) const))
    let s1 = golden();
    let mut s2 = s1.clone();
    s2.peaks.pop();
    let bs1 = binned(&s1, 0.0);
    let bs2 = binned(&s2, 0.0);
    let functor = BinnedSumAgreeingIntensities::new().unwrap();
    close(functor.score(&bs1, &bs2).unwrap(), 0.99707);
    assert_eq!(functor.score(&bs1, &bs1).unwrap(), 1.0);
}

#[test]
fn binned_sum_agreeing_intensities_self_similarity_is_one() {
    // START_SECTION((double operator()(const BinnedSpectrum &spec) const ))
    let bs1 = binned(&golden(), 0.4);
    let functor = BinnedSumAgreeingIntensities::new().unwrap();
    close(functor.self_score(&bs1).unwrap(), 1.0);
    // Derived: every bin contributes (v + v)/2 - 0 = v exactly in f32, so the
    // agreeing sum is the bin sum and the denominator is the same value.
    assert_eq!(functor.self_score(&bs1).unwrap(), 1.0);
}

#[test]
fn binned_sum_agreeing_intensities_discards_bins_that_disagree() {
    // Derived from the documented rule: a bin whose intensity difference exceeds
    // its average intensity receives weight zero, and a bin present in only one
    // spectrum always does.
    let layout = BinConfig {
        size: 1.0,
        unit: BinUnit::Absolute,
        spread: 0,
        offset: 0.0,
        ..Default::default()
    };
    let functor = BinnedSumAgreeingIntensities::new().unwrap();
    let a = BinnedSpectrum::new(&peaks(&[1.0], &[4.0]), layout).unwrap();

    // Equal intensities: (4+4)/2 - 0 = 4 over (4+4)/2 = 4, so exactly one.
    assert_eq!(functor.score(&a, &a).unwrap(), 1.0);

    // |4 - 2| = 2 <= (4 + 2)/2 = 3, so the bin keeps weight 3 - 2 = 1 and the
    // score is 1 / 3.
    let smaller = BinnedSpectrum::new(&peaks(&[1.0], &[2.0]), layout).unwrap();
    assert_eq!(functor.score(&a, &smaller).unwrap(), 1.0 / 3.0);
    assert_eq!(functor.score(&smaller, &a).unwrap(), 1.0 / 3.0);

    // |4 - 1| = 3 > (4 + 1)/2 = 2.5, so the bin is discarded entirely.
    let far = BinnedSpectrum::new(&peaks(&[1.0], &[1.0]), layout).unwrap();
    assert_eq!(functor.score(&a, &far).unwrap(), 0.0);

    // A bin in only one spectrum is always discarded.
    let elsewhere = BinnedSpectrum::new(&peaks(&[9.0], &[4.0]), layout).unwrap();
    assert_eq!(functor.score(&a, &elsewhere).unwrap(), 0.0);

    // A zero mean total intensity yields a defined zero; upstream divides by it.
    let empty = BinnedSpectrum::new(&MSSpectrum::default(), layout).unwrap();
    assert_eq!(functor.score(&empty, &empty).unwrap(), 0.0);
    let zero = BinnedSpectrum::new(&peaks(&[1.0], &[0.0]), layout).unwrap();
    assert_eq!(functor.score(&zero, &zero).unwrap(), 0.0);

    // Unlike the free function binned_sum_agreeing_intensities, the functor
    // accepts negative bins, as the source does - and then the truncation makes
    // a wholly negative spectrum score zero against itself, because
    // (v + v) / 2 - 0 = v is itself below zero. That is source behaviour, not a
    // port choice, and it is why "perfect agreement scores 1.0" holds only for
    // nonnegative bins.
    let negative = BinnedSpectrum::new(&peaks(&[1.0], &[-4.0]), layout).unwrap();
    assert_eq!(functor.score(&negative, &negative).unwrap(), 0.0);
}

#[test]
fn binned_comparisons_are_bounded_and_never_return_a_nonfinite_score() {
    // Every scorer streams over the stored bins and allocates nothing, so the
    // ceiling is a pure refusal that leaves both inputs untouched.
    assert_eq!(openms::comparison::MAX_COMPARED_BINS, 4_000_000);

    let layout = BinConfig {
        size: 1.0,
        unit: BinUnit::Absolute,
        spread: 0,
        offset: 0.0,
        ..Default::default()
    };
    // f32 accumulation of the squares overflows before f64 would, which the
    // contrast angle reports instead of propagating an infinity.
    let huge = BinnedSpectrum::new(&peaks(&[1.0, 2.0], &[3.0e38, 3.0e38]), layout).unwrap();
    assert!(matches!(
        BinnedSpectralContrastAngle::new()
            .unwrap()
            .score(&huge, &huge),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        BinnedSumAgreeingIntensities::new()
            .unwrap()
            .score(&huge, &huge),
        Err(Error::InvalidValue(_))
    ));
    // Counting bins cannot overflow, so the shared-peak count still answers.
    assert_eq!(
        BinnedSharedPeakCount::new()
            .unwrap()
            .score(&huge, &huge)
            .unwrap(),
        1.0
    );
}
