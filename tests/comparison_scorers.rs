// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Class-test coverage for the four concrete `PeakSpectrumCompareFunctor`
//! derivatives: `COMPARISON/SpectrumPrecursorComparator.h`,
//! `COMPARISON/SpectrumCheapDPCorr.h`, `COMPARISON/PeakAlignment.h` and
//! `COMPARISON/SpectraSTSimilarityScore.h`.
//!
//! Every `START_SECTION` of the four upstream class tests has a counterpart
//! here, named after it in a comment. Values marked "class-test literal" are
//! transcribed from the upstream test and are evidence tier 3; values marked
//! "derived" are computed from the algorithm's definition and are stronger.

use openms::comparison::{
    BinConfig, BinnedSpectrum, MAX_ALIGNMENT_MATRIX_CELLS, PeakAlignment,
    PeakSpectrumCompareFunctor, SpectraSTSimilarityScore, SpectraStPreprocessing,
    SpectrumCheapDPCorr, SpectrumPrecursorComparator,
};
use openms::format::dta;
use openms::param::{Param, ParamValue};
use openms::{Error, MSSpectrum, Peak1D, Precursor};

/// `Transformers_tests.dta`, 121 peaks, precursor MH+ 739.771 at charge 2.
fn transformers_1() -> MSSpectrum {
    dta::read(include_bytes!("data/Transformers_tests.dta").as_slice()).unwrap()
}

/// `Transformers_tests_2.dta`, 93 peaks, precursor MH+ 739.308 at charge 2.
fn transformers_2() -> MSSpectrum {
    dta::read(include_bytes!("data/comparison_transformers_2.dta").as_slice()).unwrap()
}

/// `PILISSequenceDB_DFPIANGER_1.dta`, 127 peaks, precursor 1019.74 at charge 1.
fn dfpianger() -> MSSpectrum {
    dta::read(include_bytes!("data/comparison_dfpianger.dta").as_slice()).unwrap()
}

/// The peak list of the `index`-th spectrum of the retained upstream
/// `SpectraSTSimilarityScore_1.msp`.
///
/// Only `Num peaks:` and the tab-separated `m/z intensity annotation` lines
/// that follow it are read; this is not an MSP reader and does not pretend to
/// be one. The first two spectra are identical and the third is the same
/// intensity pattern at well-separated m/z.
fn msp(index: usize) -> MSSpectrum {
    let text = include_str!("data/comparison_spectrast_1.msp");
    let mut spectra: Vec<Vec<Peak1D>> = Vec::new();
    let mut remaining = 0usize;
    for line in text.lines() {
        if let Some(count) = line.strip_prefix("Num peaks:") {
            remaining = count.trim().parse().unwrap();
            spectra.push(Vec::new());
            continue;
        }
        if remaining == 0 {
            continue;
        }
        let mut fields = line.split('\t');
        let mz: f64 = fields.next().unwrap().trim().parse().unwrap();
        let intensity: f32 = fields.next().unwrap().trim().parse().unwrap();
        spectra.last_mut().unwrap().push(Peak1D::new(mz, intensity));
        remaining -= 1;
    }
    MSSpectrum::from_peaks(spectra[index].clone())
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} is not within {tolerance} of {expected}"
    );
}

fn spectrum(mzs: &[f64], intensities: &[f32]) -> MSSpectrum {
    MSSpectrum::from_peaks(
        mzs.iter()
            .zip(intensities)
            .map(|(&mz, &intensity)| Peak1D::new(mz, intensity))
            .collect(),
    )
}

// ---------------------------------------------------------------------------
// COMPARISON/SpectrumPrecursorComparator.h - six sections
// ---------------------------------------------------------------------------

#[test]
fn precursor_comparator_construction_copy_and_assignment() {
    // START_SECTION(SpectrumPrecursorComparator())
    // START_SECTION(~SpectrumPrecursorComparator())
    // Construction registers the derived name and the single "window" default;
    // destruction is drop glue with no observable effect.
    let functor = SpectrumPrecursorComparator::new().unwrap();
    assert_eq!(functor.name(), "SpectrumPrecursorComparator");
    assert_eq!(functor.handler().parameters().size(), 1);
    // The source registers the integer 2, not the float 2.0.
    assert_eq!(
        *functor.handler().parameters().value("window").unwrap(),
        ParamValue::Integer(2)
    );
    assert_eq!(
        functor
            .handler()
            .parameters()
            .description("window")
            .unwrap(),
        "Allowed deviation between precursor peaks."
    );
    assert_eq!(functor.window().unwrap(), 2.0);
    assert_eq!(SpectrumPrecursorComparator::default(), functor);

    // START_SECTION(SpectrumPrecursorComparator(const SpectrumPrecursorComparator& source))
    // START_SECTION(SpectrumPrecursorComparator& operator=(const ... & source))
    // Both forward to DefaultParamHandler, so name and parameters survive.
    let copy = functor.clone();
    assert_eq!(copy.name(), functor.name());
    assert_eq!(copy.handler().parameters(), functor.handler().parameters());
    let mut assigned = SpectrumPrecursorComparator::new().unwrap();
    assigned.clone_from(&functor);
    assert_eq!(assigned.name(), functor.name());
    assert_eq!(
        assigned.handler().parameters(),
        functor.handler().parameters()
    );
}

#[test]
fn precursor_comparator_scores_the_parent_mass_distance() {
    // START_SECTION(double operator () (const PeakSpectrum& a, const PeakSpectrum& b) const)
    let a = transformers_1();
    let b = transformers_2();
    let functor = SpectrumPrecursorComparator::new().unwrap();
    // Class-test literal 1.7685. Derived: DTA turns MH+ 739.771 and 739.308 at
    // charge 2 into 370.3891382333855 and 370.1576382333855, whose distance is
    // 0.2315, and 2 - 0.2315 = 1.7685.
    close(functor.score(&a, &b).unwrap(), 1.7685, 1e-9);
    close(
        functor.score(&a, &b).unwrap(),
        2.0 - (a.precursors[0].mz - b.precursors[0].mz).abs(),
        0.0,
    );
    // Class-test literal 2. Derived: a zero distance leaves the whole window.
    assert_eq!(functor.score(&a, &a).unwrap(), 2.0);

    // START_SECTION(double operator () (const PeakSpectrum& a) const)
    // The source's one-argument overload is operator()(spec, spec).
    assert_eq!(functor.self_score(&a).unwrap(), 2.0);
}

#[test]
fn precursor_comparator_missing_precursor_convention_and_guards() {
    let functor = SpectrumPrecursorComparator::new().unwrap();
    // A spectrum with no precursor contributes m/z zero, so two such spectra
    // score the full window. This is the source convention, not an accident.
    let empty = MSSpectrum::default();
    assert_eq!(functor.score(&empty, &empty).unwrap(), 2.0);
    let distant = MSSpectrum {
        precursors: vec![Precursor::new(5.0, 1)],
        ..MSSpectrum::default()
    };
    assert_eq!(functor.score(&empty, &distant).unwrap(), 0.0);

    // Only the first precursor is read.
    let two = MSSpectrum {
        precursors: vec![Precursor::new(100.0, 1), Precursor::new(500.0, 1)],
        ..MSSpectrum::default()
    };
    let one = MSSpectrum {
        precursors: vec![Precursor::new(101.0, 1)],
        ..MSSpectrum::default()
    };
    assert_eq!(functor.score(&two, &one).unwrap(), 1.0);

    // A negative window is refused; the source would return negative scores.
    let mut negative = SpectrumPrecursorComparator::new().unwrap();
    let mut parameters = Param::new();
    parameters
        .set_value("window", ParamValue::Integer(-1), "", &[])
        .unwrap();
    negative.handler_mut().set_parameters(&parameters).unwrap();
    assert!(matches!(
        negative.score(&empty, &empty),
        Err(Error::InvalidValue(_))
    ));
}

// ---------------------------------------------------------------------------
// COMPARISON/SpectrumCheapDPCorr.h - nine sections
// ---------------------------------------------------------------------------

#[test]
fn cheap_dp_corr_construction_copy_and_assignment() {
    // START_SECTION(SpectrumCheapDPCorr())
    // START_SECTION(~SpectrumCheapDPCorr())
    let functor = SpectrumCheapDPCorr::new().unwrap();
    assert_eq!(functor.name(), "SpectrumCheapDPCorr");
    assert_eq!(functor.handler().parameters().size(), 3);
    assert_eq!(
        *functor.handler().parameters().value("variation").unwrap(),
        ParamValue::Float(0.001)
    );
    assert_eq!(
        *functor.handler().parameters().value("int_cnt").unwrap(),
        ParamValue::Integer(0)
    );
    assert_eq!(
        *functor.handler().parameters().value("keeppeaks").unwrap(),
        ParamValue::Integer(0)
    );
    // The source's constructor sets factor_ = 0.5 before defaultsToParam_().
    assert_eq!(functor.factor(), 0.5);
    assert!(functor.last_consensus().is_empty());
    assert!(functor.peak_map().is_empty());

    // START_SECTION(SpectrumCheapDPCorr(const SpectrumCheapDPCorr& source))
    // START_SECTION(SpectrumCheapDPCorr& operator=(const SpectrumCheapDPCorr& source))
    let copy = functor.clone();
    assert_eq!(copy.handler().parameters(), functor.handler().parameters());
    assert_eq!(copy.name(), functor.name());
    let assigned = functor.clone();
    assert_eq!(
        assigned.handler().parameters(),
        functor.handler().parameters()
    );
    assert_eq!(assigned.name(), functor.name());
}

#[test]
fn cheap_dp_corr_reproduces_the_upstream_golden_scores() {
    // START_SECTION(double operator () (const PeakSpectrum& a, const PeakSpectrum& b) const)
    let a = transformers_1();
    let b = transformers_2();
    let mut functor = SpectrumCheapDPCorr::new().unwrap();
    // Class-test literals 10145.4 and 12295.5 at TOLERANCE_ABSOLUTE(0.1). The
    // tighter expectations are an independent reimplementation of the scan,
    // dynprog_ and Boost's normal pdf in Python, which reproduced both literals
    // as well as the consensus and peak-map sizes asserted below.
    let cross = functor.compare(&a, &b).unwrap();
    close(cross, 10145.4, 0.1);
    close(cross, 10145.449278148695, 1e-6);
    let self_score = functor.compare(&a, &a).unwrap();
    close(self_score, 12295.5, 0.1);
    close(self_score, 12295.522100159595, 1e-6);

    // The class test repeats both on a second, freshly constructed object.
    let mut fresh = SpectrumCheapDPCorr::new().unwrap();
    close(fresh.compare(&a, &b).unwrap(), 10145.4, 0.1);
    close(fresh.compare(&a, &a).unwrap(), 12295.5, 0.1);

    // START_SECTION(const PeakSpectrum& lastconsensus() const)
    // Class-test literal 121 after the (spec1, spec1) call. Derived: a spectrum
    // aligned against itself pairs all 121 of its peaks, and with keeppeaks
    // cleared only paired peaks enter the consensus.
    assert_eq!(functor.last_consensus().len(), 121);
    assert_eq!(functor.last_consensus().len(), a.len());
    // The consensus carries one precursor at the mean of the two inputs.
    assert_eq!(functor.last_consensus().precursors.len(), 1);
    close(
        functor.last_consensus().precursors[0].mz,
        a.precursors[0].mz,
        1e-12,
    );
    assert_eq!(
        functor.last_consensus().precursors[0].charge,
        a.precursors[0].charge
    );

    // START_SECTION((Map<UInt, UInt> getPeakMap() const))
    // Class-test literal 121. Derived: every peak of the self-alignment is
    // mapped, and to itself.
    assert_eq!(functor.peak_map().len(), 121);
    for (&from, &to) in functor.peak_map() {
        assert_eq!(from, to);
    }

    // START_SECTION(double operator () (const PeakSpectrum& a) const)
    close(functor.self_score(&a).unwrap(), 12295.5, 0.1);
    // The read-only trait entry point agrees with the recording one exactly.
    assert_eq!(functor.score(&a, &b).unwrap(), cross);
}

#[test]
fn cheap_dp_corr_factor_is_range_checked_and_reset_after_a_comparison() {
    // START_SECTION(void setFactor(double f))
    let mut functor = SpectrumCheapDPCorr::new().unwrap();
    functor.set_factor(0.3).unwrap();
    assert_eq!(functor.factor(), 0.3);
    // Both bounds are exclusive upstream.
    assert!(matches!(
        functor.set_factor(1.1),
        Err(Error::InvalidRange(_))
    ));
    assert!(matches!(
        functor.set_factor(1.0),
        Err(Error::InvalidRange(_))
    ));
    assert!(matches!(
        functor.set_factor(0.0),
        Err(Error::InvalidRange(_))
    ));
    assert_eq!(functor.factor(), 0.3);

    // The factor weights the consensus only, so a comparison run at 0.3 scores
    // exactly as one run at the default, and the source resets it afterwards.
    let a = spectrum(&[100.0, 200.0], &[10.0, 20.0]);
    let b = spectrum(&[100.0, 200.0], &[30.0, 40.0]);
    let weighted = functor.compare(&a, &b).unwrap();
    assert_eq!(functor.factor(), 0.5);
    let balanced = functor.compare(&a, &b).unwrap();
    assert_eq!(weighted, balanced);

    // The first consensus used factor 0.3, the second 0.5.
    let mut at_three_tenths = SpectrumCheapDPCorr::new().unwrap();
    at_three_tenths.set_factor(0.3).unwrap();
    at_three_tenths.compare(&a, &b).unwrap();
    close(
        f64::from(at_three_tenths.last_consensus().peaks[0].intensity),
        10.0 * 0.7 + 30.0 * 0.3,
        1e-5,
    );
    close(
        f64::from(functor.last_consensus().peaks[0].intensity),
        10.0 * 0.5 + 30.0 * 0.5,
        1e-5,
    );
}

#[test]
fn cheap_dp_corr_intensity_terms_and_kept_peaks() {
    // int_cnt selects the intensity term; two peaks at the same m/z isolate it
    // from the alignment entirely, so each branch has a closed form.
    let a = spectrum(&[100.0], &[9.0]);
    let b = spectrum(&[100.0], &[4.0]);
    // The Gaussian factor is pdf(N(0, 100 * 0.001), 0) = 1 / (0.1 * sqrt(2 pi)).
    let density = 1.0 / (0.1 * (2.0 * std::f64::consts::PI).sqrt());
    for (int_cnt, expected) in [
        (0, density * 36.0),
        (1, density * 6.0),
        (2, density * 13.0),
        (3, density * (6.5 - 5.0)),
    ] {
        let mut functor = SpectrumCheapDPCorr::new().unwrap();
        let mut parameters = functor.handler().parameters().clone();
        parameters
            .set_value("int_cnt", ParamValue::Integer(int_cnt), "", &[])
            .unwrap();
        functor.handler_mut().set_parameters(&parameters).unwrap();
        close(functor.score(&a, &b).unwrap(), expected, 1e-9);
    }

    // The source returns -1 for any other int_cnt, behind a "// TODO exception";
    // a caller summing scores cannot tell that from a real contribution.
    let mut broken = SpectrumCheapDPCorr::new().unwrap();
    let mut parameters = broken.handler().parameters().clone();
    parameters
        .set_value("int_cnt", ParamValue::Integer(4), "", &[])
        .unwrap();
    broken.handler_mut().set_parameters(&parameters).unwrap();
    assert!(matches!(broken.score(&a, &b), Err(Error::InvalidValue(_))));

    // keeppeaks adds the unpaired peaks to the consensus, weighted by the
    // factor of the spectrum they came from. Only peaks the scan reaches are
    // considered: the loop stops as soon as either list is exhausted, so a
    // trailing tail is dropped whether or not the flag is set.
    let left = spectrum(&[100.0, 300.0], &[10.0, 30.0]);
    let right = spectrum(&[200.0, 300.0], &[20.0, 40.0]);
    let mut keeping = SpectrumCheapDPCorr::new().unwrap();
    keeping.compare(&left, &right).unwrap();
    assert_eq!(keeping.last_consensus().len(), 1);
    let mut parameters = keeping.handler().parameters().clone();
    parameters
        .set_value("keeppeaks", ParamValue::Integer(1), "", &[])
        .unwrap();
    keeping.handler_mut().set_parameters(&parameters).unwrap();
    keeping.compare(&left, &right).unwrap();
    // 100 Th and 200 Th have no partner and survive at (1 - factor) and factor
    // of their intensity respectively; 300 Th pairs. The source would read an
    // uninitialised member for the same decision inside dynprog_; this port
    // reads the parameter in both places.
    assert_eq!(keeping.last_consensus().len(), 3);
    close(
        f64::from(keeping.last_consensus().peaks[0].intensity),
        5.0,
        1e-5,
    );
    close(
        f64::from(keeping.last_consensus().peaks[1].intensity),
        10.0,
        1e-5,
    );
    close(
        f64::from(keeping.last_consensus().peaks[2].intensity),
        35.0,
        1e-5,
    );
}

#[test]
fn cheap_dp_corr_refuses_undefined_input() {
    let functor = SpectrumCheapDPCorr::new().unwrap();
    // Unsorted peaks mis-align silently upstream.
    let unsorted = spectrum(&[200.0, 100.0], &[1.0, 1.0]);
    let sorted = spectrum(&[100.0, 200.0], &[1.0, 1.0]);
    assert!(matches!(
        functor.score(&unsorted, &sorted),
        Err(Error::UnsortedData)
    ));
    // Two peaks at m/z zero make the Gaussian scale zero, where Boost raises a
    // domain error rather than dividing by it.
    let at_zero = spectrum(&[0.0], &[1.0]);
    assert!(matches!(
        functor.score(&at_zero, &at_zero),
        Err(Error::InvalidValue(_))
    ));
    // A non-positive variation is the same degenerate Gaussian.
    let mut zero_variation = SpectrumCheapDPCorr::new().unwrap();
    let mut parameters = zero_variation.handler().parameters().clone();
    parameters
        .set_value("variation", ParamValue::Float(0.0), "", &[])
        .unwrap();
    zero_variation
        .handler_mut()
        .set_parameters(&parameters)
        .unwrap();
    assert!(matches!(
        zero_variation.score(&sorted, &sorted),
        Err(Error::InvalidValue(_))
    ));
    // Two empty spectra never enter the scan and score zero.
    assert_eq!(
        functor
            .score(&MSSpectrum::default(), &MSSpectrum::default())
            .unwrap(),
        0.0
    );
}

#[test]
fn cheap_dp_corr_bounds_its_dynamic_programming_block() {
    // A variation of 1 is the source's documented maximum and makes one run
    // span both spectra, which is the O(n*n) case its own description warns
    // about. 2000 peaks a side is just over the millionth cell.
    let wide: Vec<Peak1D> = (0..2000)
        .map(|i| Peak1D::new(1000.0 + f64::from(i), 1.0))
        .collect();
    let wide = MSSpectrum::from_peaks(wide);
    let mut functor = SpectrumCheapDPCorr::new().unwrap();
    let mut parameters = functor.handler().parameters().clone();
    parameters
        .set_value("variation", ParamValue::Float(1.0), "", &[])
        .unwrap();
    functor.handler_mut().set_parameters(&parameters).unwrap();
    let refused = functor.compare(&wide, &wide).unwrap_err();
    assert!(
        refused.to_string().contains("dynamic-programming cell"),
        "{refused}"
    );
    // The refusal left the recorded state untouched.
    assert!(functor.last_consensus().is_empty());
    assert!(functor.peak_map().is_empty());
}

// ---------------------------------------------------------------------------
// COMPARISON/PeakAlignment.h - seven sections
// ---------------------------------------------------------------------------

#[test]
fn peak_alignment_construction_copy_and_assignment() {
    // START_SECTION(PeakAlignment())
    // START_SECTION(~PeakAlignment())
    let functor = PeakAlignment::new().unwrap();
    // The one derivative that never calls setName, so the base name survives.
    assert_eq!(functor.name(), "PeakSpectrumCompareFunctor");
    assert_eq!(functor.handler().parameters().size(), 4);
    assert_eq!(
        *functor.handler().parameters().value("epsilon").unwrap(),
        ParamValue::Float(0.2)
    );
    assert_eq!(
        *functor.handler().parameters().value("normalized").unwrap(),
        ParamValue::Integer(1)
    );
    assert_eq!(
        *functor
            .handler()
            .parameters()
            .value("heuristic_level")
            .unwrap(),
        ParamValue::Integer(0)
    );
    assert_eq!(
        *functor
            .handler()
            .parameters()
            .value("precursor_mass_tolerance")
            .unwrap(),
        ParamValue::Float(3.0)
    );

    // START_SECTION((PeakAlignment(const PeakAlignment &source)))
    // START_SECTION((PeakAlignment& operator=(const PeakAlignment &source)))
    let copy = functor.clone();
    assert_eq!(copy.name(), functor.name());
    assert_eq!(copy.handler().parameters(), functor.handler().parameters());
    let assigned = functor.clone();
    assert_eq!(assigned.name(), functor.name());
    assert_eq!(
        assigned.handler().parameters(),
        functor.handler().parameters()
    );
}

#[test]
fn peak_alignment_reproduces_the_upstream_golden_score() {
    // START_SECTION((double operator()(const PeakSpectrum &spec1, const PeakSpectrum &spec2) const))
    let s1 = dfpianger();
    let mut s2 = dfpianger();
    s2.peaks.pop();
    let functor = PeakAlignment::new().unwrap();
    // Class-test literal 0.997477. The tighter expectation is an independent
    // reimplementation of the matrix fill, the DBL_MIN-seeded border scan and
    // the source's own exp(-|d|/2 * sigma * sigma) position term.
    close(functor.score(&s1, &s2).unwrap(), 0.997477, 1e-6);
    close(functor.score(&s1, &s2).unwrap(), 0.9974770186204426, 1e-12);

    // Empty spectra: the class test asserts zero, and it is the precursor
    // shortcut that produces it - the empty spectrum has no precursor, so its
    // m/z reads as zero and the distance to 1019.74 exceeds the 3 Th default.
    let empty = MSSpectrum::default();
    assert_eq!(functor.score(&empty, &s2).unwrap(), 0.0);
    assert_eq!(functor.score(&s1, &empty).unwrap(), 0.0);
    // With that shortcut disarmed the source divides by a zero pair count and
    // returns infinity; here the undefined case is an error.
    assert!(matches!(
        functor.score(&empty, &MSSpectrum::default()),
        Err(Error::InvalidValue(_))
    ));

    // START_SECTION((double operator()(const PeakSpectrum &spec) const))
    // Class-test literal 1. Derived: the diagonal of a self-alignment scores
    // exactly the self-alignment sum, so the quotient is exactly one.
    assert_eq!(functor.self_score(&s1).unwrap(), 1.0);
}

#[test]
fn peak_alignment_traceback_is_the_diagonal_for_a_self_alignment() {
    // START_SECTION((vector< pair<Size,Size> > getAlignmentTraceback(const PeakSpectrum &spec1, const PeakSpectrum &spec2) const))
    let s1 = dfpianger();
    let functor = PeakAlignment::new().unwrap();
    let traceback = functor.alignment_traceback(&s1, &s1).unwrap();
    // Class-test expectation: 127 pairs (i, i). Derived: aligning a spectrum
    // with itself, every diagonal step is strictly better than either gap.
    assert_eq!(traceback.len(), 127);
    assert_eq!(traceback.len(), s1.len());
    for (index, &(row, column)) in traceback.iter().enumerate() {
        assert_eq!(row, index);
        assert_eq!(column, index);
    }

    // The traceback is ascending and never revisits an index on either side.
    let mut previous = None;
    for &(row, column) in &traceback {
        if let Some((last_row, last_column)) = previous {
            assert!(row > last_row && column > last_column);
        }
        previous = Some((row, column));
    }

    // Unlike the score, the traceback has no zero-variance guard, so two
    // one-peak spectra at the same m/z hit the source's division by a zero
    // sigma; that is refused rather than returning infinities.
    let single = spectrum(&[500.0], &[1.0]);
    assert!(matches!(
        functor.alignment_traceback(&single, &single),
        Err(Error::InvalidValue(_))
    ));
    // The score does guard the zero variance, but its DBL_MIN substitute is a
    // cure worse than the disease: the position term becomes
    // 1 / (DBL_MIN * sqrt(2 pi)) ~ 1.8e307, so the product of the two
    // self-alignment scores overflows to infinity for every nonzero f32
    // intensity - the smallest subnormal, 1.4e-45, still leaves the square
    // above 1.8e308 - and the source silently reports 0, complete
    // dissimilarity, for a spectrum compared with itself. The guard branch is
    // therefore unusable upstream, and every way into it is an error here.
    for intensity in [f32::MIN_POSITIVE, 1.0, 1.0e30] {
        assert!(matches!(
            functor.self_score(&spectrum(&[500.0], &[intensity])),
            Err(Error::InvalidValue(_))
        ));
    }
    // A zero intensity instead makes both self-alignment scores zero, which is
    // the source's other division by zero.
    assert!(matches!(
        functor.self_score(&spectrum(&[500.0], &[0.0])),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn peak_alignment_shortcuts_and_unread_normalized_flag() {
    let functor = PeakAlignment::new().unwrap();
    let near = MSSpectrum {
        precursors: vec![Precursor::new(500.0, 2)],
        ..spectrum(&[100.0, 200.0], &[1.0, 2.0])
    };
    let far = MSSpectrum {
        precursors: vec![Precursor::new(504.0, 2)],
        ..spectrum(&[100.0, 200.0], &[1.0, 2.0])
    };
    // 4 Th apart is beyond the 3 Th default, so the precursor shortcut fires.
    assert_eq!(functor.score(&near, &far).unwrap(), 0.0);
    let mut wide = PeakAlignment::new().unwrap();
    let mut parameters = wide.handler().parameters().clone();
    parameters
        .set_value("precursor_mass_tolerance", ParamValue::Float(5.0), "", &[])
        .unwrap();
    wide.handler_mut().set_parameters(&parameters).unwrap();
    assert_eq!(wide.score(&near, &far).unwrap(), 1.0);

    // The heuristic shortcut: with level 1 only the strongest peak of each
    // spectrum is considered, and 200 Th against 900 Th shares nothing.
    let strong_low = spectrum(&[200.0, 900.0], &[10.0, 1.0]);
    let strong_high = spectrum(&[200.0, 900.0], &[1.0, 10.0]);
    let mut heuristic = PeakAlignment::new().unwrap();
    let mut parameters = heuristic.handler().parameters().clone();
    parameters
        .set_value("heuristic_level", ParamValue::Integer(1), "", &[])
        .unwrap();
    heuristic.handler_mut().set_parameters(&parameters).unwrap();
    assert_eq!(heuristic.score(&strong_low, &strong_high).unwrap(), 0.0);
    // Level 2 sees both peaks of each and finds the match.
    parameters
        .set_value("heuristic_level", ParamValue::Integer(2), "", &[])
        .unwrap();
    heuristic.handler_mut().set_parameters(&parameters).unwrap();
    assert!(heuristic.score(&strong_low, &strong_high).unwrap() > 0.0);
    // A negative level is a huge unsigned value upstream; here it is refused.
    parameters
        .set_value("heuristic_level", ParamValue::Integer(-1), "", &[])
        .unwrap();
    heuristic.handler_mut().set_parameters(&parameters).unwrap();
    assert!(matches!(
        heuristic.score(&strong_low, &strong_high),
        Err(Error::InvalidValue(_))
    ));

    // "normalized" is registered but never read: clearing it changes nothing.
    let plain = functor.score(&near, &near).unwrap();
    let mut cleared = PeakAlignment::new().unwrap();
    let mut parameters = cleared.handler().parameters().clone();
    parameters
        .set_value("normalized", ParamValue::Integer(0), "", &[])
        .unwrap();
    cleared.handler_mut().set_parameters(&parameters).unwrap();
    assert_eq!(cleared.score(&near, &near).unwrap(), plain);
}

#[test]
fn peak_alignment_bounds_its_matrix_and_rejects_unsorted_peaks() {
    let functor = PeakAlignment::new().unwrap();
    let unsorted = spectrum(&[200.0, 100.0], &[1.0, 1.0]);
    let sorted = spectrum(&[100.0, 200.0], &[1.0, 1.0]);
    assert!(matches!(
        functor.score(&unsorted, &sorted),
        Err(Error::UnsortedData)
    ));

    // 2100 x 2100 peaks ask for 4_414_201 cells, past the ceiling; the source
    // allocates that plus a second matrix without checking anything.
    let side = 2100;
    assert!((side + 1) * (side + 1) > MAX_ALIGNMENT_MATRIX_CELLS);
    let wide = MSSpectrum::from_peaks(
        (0..side)
            .map(|i| Peak1D::new(100.0 + i as f64, 1.0))
            .collect(),
    );
    let refused = functor.score(&wide, &wide).unwrap_err();
    assert!(
        refused.to_string().contains("alignment-matrix cell"),
        "{refused}"
    );
}

// ---------------------------------------------------------------------------
// COMPARISON/SpectraSTSimilarityScore.h - twelve sections
// ---------------------------------------------------------------------------

#[test]
fn spectrast_construction_copy_and_assignment() {
    // START_SECTION(SpectraSTSimilarityScore())
    // START_SECTION(~SpectraSTSimilarityScore())
    let functor = SpectraSTSimilarityScore::new().unwrap();
    assert_eq!(functor.name(), "SpectraSTSimilarityScore");
    // The constructor calls setName and stops: no defaults, and no
    // defaultsToParam_() either.
    assert!(functor.handler().parameters().is_empty());
    assert!(functor.handler().defaults().is_empty());

    // START_SECTION(SpectraSTSimilarityScore(const SpectraSTSimilarityScore& source))
    // START_SECTION(SpectraSTSimilarityScore& operator = (const ... & source))
    let copy = functor.clone();
    assert_eq!(copy.name(), functor.name());
    assert_eq!(copy.handler().parameters(), functor.handler().parameters());
    let assigned = functor.clone();
    assert_eq!(assigned.name(), functor.name());
    assert_eq!(
        assigned.handler().parameters(),
        functor.handler().parameters()
    );
}

#[test]
fn spectrast_peak_spectrum_dot_products() {
    // START_SECTION(double operator () (const PeakSpectrum& spec) const)
    // START_SECTION(double operator () (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const)
    let functor = SpectraSTSimilarityScore::new().unwrap();
    let s1 = msp(0);
    let s2 = msp(1);
    let s3 = msp(2);
    assert_eq!(s1.peaks, s2.peaks);
    // Class-test literal 1 at TOLERANCE_ABSOLUTE(0.01). Derived: a normalised
    // vector dotted with itself is one up to the f32 reduction's rounding.
    close(functor.self_score(&s1).unwrap(), 1.0, 1e-6);
    close(functor.score(&s1, &s2).unwrap(), 1.0, 1e-6);
    // Class-test literal 0. Derived: with bin width 1 and spread 1 the two peak
    // sets occupy disjoint bins, so the intersection the sparse dot walks is
    // empty and the score is exactly zero.
    assert_eq!(functor.score(&s1, &s3).unwrap(), 0.0);

    // START_SECTION((double operator()(const BinnedSpectrum &bin1, const BinnedSpectrum &bin2) const))
    close(
        functor
            .dot(
                &functor.transform(&s1).unwrap(),
                &functor.transform(&s2).unwrap(),
            )
            .unwrap(),
        1.0,
        1e-6,
    );
    assert_eq!(
        functor
            .dot(
                &functor.transform(&s1).unwrap(),
                &functor.transform(&s3).unwrap()
            )
            .unwrap(),
        0.0
    );
    // The peak-spectrum overload is exactly dot(transform, transform).
    assert_eq!(
        functor.score(&s1, &s3).unwrap(),
        functor
            .dot(
                &functor.transform(&s1).unwrap(),
                &functor.transform(&s3).unwrap()
            )
            .unwrap()
    );
    // Incompatible binning is refused; upstream hands Eigen mismatched vectors.
    let other = BinnedSpectrum::new(
        &s1,
        BinConfig {
            size: 2.0,
            ..SpectraSTSimilarityScore::bin_config()
        },
    )
    .unwrap();
    assert!(matches!(
        functor.dot(&functor.transform(&s1).unwrap(), &other),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn spectrast_transform_normalises_the_binned_vector() {
    // START_SECTION(BinnedSpectrum transform(const PeakSpectrum& spec))
    let functor = SpectraSTSimilarityScore::new().unwrap();
    let s1 = spectrum(&[0.5, 1.5, 2.5, 3.5], &[1.0, 0.0, 2.0, 3.0]);
    let binned = functor.transform(&s1).unwrap();
    // Derived: floor(mz + 0.4) puts the four peaks in bins 0, 1, 2 and 3, and
    // spread 1 adds each intensity to its neighbours, with bin 0 unable to
    // spread downwards. That gives stored bins 1, 3, 5, 5, 3 at indices 0..4,
    // whose norm is sqrt(69).
    let norm = 69.0_f32.sqrt();
    let expected = [1.0_f32, 3.0, 5.0, 5.0, 3.0];
    assert_eq!(binned.bins().len(), expected.len());
    for ((&index, &value), (position, &raw)) in
        binned.bins().iter().zip(expected.iter().enumerate())
    {
        assert_eq!(index, position);
        assert_eq!(value, raw / norm);
    }
    // Class-test literals 0.1205, 0.3614, 0.602 and 0.602 at
    // TOLERANCE_ABSOLUTE(0.01), read off the first four stored coefficients.
    let values: Vec<f64> = binned.bins().values().map(|&v| f64::from(v)).collect();
    close(values[0], 0.1205, 0.01);
    close(values[1], 0.3614, 0.01);
    close(values[2], 0.602, 0.01);
    close(values[3], 0.602, 0.01);

    // A spectrum with nothing to normalise by is refused; the source divides by
    // that zero and fills the vector with NaN.
    assert!(matches!(
        functor.transform(&MSSpectrum::default()),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        functor.transform(&spectrum(&[100.0], &[0.0])),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn spectrast_dot_bias_measures_domination_by_few_bins() {
    // START_SECTION(double dot_bias(const BinnedSpectrum &bin1, const BinnedSpectrum &bin2, double dot_product=-1) const)
    let functor = SpectraSTSimilarityScore::new().unwrap();
    let s1 = spectrum(&[1.0, 2.0, 3.0, 4.0], &[1.0, 0.0, 2.0, 3.0]);
    let s2 = spectrum(&[1.0, 2.0, 3.0, 4.0, 5.0], &[0.0, 4.0, 5.0, 6.0, 0.0]);
    let config = SpectraSTSimilarityScore::bin_config();
    let bin1 = BinnedSpectrum::new(&s1, config).unwrap();
    let bin2 = BinnedSpectrum::new(&s2, config).unwrap();
    // Class-test literal 98.585 at TOLERANCE_ABSOLUTE(0.01). Derived: the two
    // raw binned vectors are (1, 1, 3, 5, 5, 3) and (0, 4, 9, 15, 11, 6, 0),
    // whose element-wise products over the shared indices are 0, 4, 27, 75, 55
    // and 18; the sum of their squares is 9719.
    let expected = f64::from(9719.0_f32.sqrt());
    close(
        functor.dot_bias(&bin1, &bin2, Some(1.0)).unwrap(),
        98.585,
        0.01,
    );
    assert_eq!(functor.dot_bias(&bin1, &bin2, Some(1.0)).unwrap(), expected);
    // Symmetric in its two arguments.
    assert_eq!(
        functor.dot_bias(&bin2, &bin1, Some(1.0)).unwrap(),
        functor.dot_bias(&bin1, &bin2, Some(1.0)).unwrap()
    );

    // Two spectra sharing no bin: the recomputed dot product is zero, which the
    // source's own guard turns into a defined zero rather than a NaN.
    let orthogonal1 = BinnedSpectrum::new(&spectrum(&[10.0], &[1.0]), config).unwrap();
    let orthogonal2 = BinnedSpectrum::new(&spectrum(&[20.0], &[1.0]), config).unwrap();
    assert_eq!(
        functor.dot_bias(&orthogonal1, &orthogonal2, None).unwrap(),
        0.0
    );
    // -1 and 0 are the source's sentinels for "recompute the dot product".
    assert_eq!(
        functor
            .dot_bias(&orthogonal1, &orthogonal2, Some(-1.0))
            .unwrap(),
        0.0
    );
    assert_eq!(
        functor
            .dot_bias(&orthogonal1, &orthogonal2, Some(0.0))
            .unwrap(),
        0.0
    );
    // A supplied denominator of 2 halves the ratio.
    close(
        functor.dot_bias(&bin1, &bin2, Some(2.0)).unwrap(),
        expected / 2.0,
        1e-12,
    );
    // NaN slips past the source's `<= 0.0` guard; here it is refused.
    assert!(matches!(
        functor.dot_bias(&bin1, &bin2, Some(f64::NAN)),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn spectrast_preprocess_filters_and_square_roots_the_intensities() {
    // START_SECTION(bool preprocess(PeakSpectrum &spec, float remove_peak_intensity_threshold=2.01, UInt cut_peaks_below=1000, Size min_peak_number=5, Size max_peak_number=150))
    let functor = SpectraSTSimilarityScore::new().unwrap();
    let mut s1 = msp(0);
    // Class-test call preprocess(s1, 2, 10000) expecting six survivors.
    // Derived: the ten intensities are 2, 2, 5, 2, 3, 4, 5, 5, 2, 3, the base
    // peak is 5 so the relative floor is 5/10000, and exactly the six peaks
    // above intensity 2 survive.
    assert!(
        functor
            .preprocess(
                &mut s1,
                SpectraStPreprocessing {
                    remove_peak_intensity_threshold: 2.0,
                    cut_peaks_below: 10000,
                    ..SpectraStPreprocessing::default()
                },
            )
            .unwrap()
    );
    assert_eq!(s1.len(), 6);
    // The SpectraST intensity exponent is 0.5, applied in f32.
    assert_eq!(s1.peaks[0].intensity, 5.0_f32.sqrt());
    assert_eq!(s1.peaks[0].mz, 430.3);

    // min_peak_number 12 rejects the same six survivors.
    let mut s2 = msp(1);
    assert!(
        !functor
            .preprocess(
                &mut s2,
                SpectraStPreprocessing {
                    remove_peak_intensity_threshold: 2.0,
                    cut_peaks_below: 1000,
                    min_peak_number: 12,
                    ..SpectraStPreprocessing::default()
                },
            )
            .unwrap()
    );
    assert_eq!(s2.len(), 6);

    // max_peak_number 8 stops after examining eight peaks, all of which pass.
    let mut s3 = msp(2);
    assert!(
        functor
            .preprocess(
                &mut s3,
                SpectraStPreprocessing {
                    remove_peak_intensity_threshold: 1.0,
                    cut_peaks_below: 10000,
                    min_peak_number: 5,
                    max_peak_number: 8,
                },
            )
            .unwrap()
    );
    assert_eq!(s3.len(), 8);
    // It is a prefix in m/z, not the eight strongest peaks: the two dropped
    // peaks are the highest m/z, and one of them is stronger than the survivor
    // at 200 Th.
    assert_eq!(s3.peaks[0].mz, 200.0);
    assert_eq!(s3.peaks[7].mz, 900.0);

    // The whole spectrum is replaced, so metadata does not survive; that is
    // what `spec = tmp` does upstream.
    let mut with_metadata = MSSpectrum {
        precursors: vec![Precursor::new(500.0, 2)],
        native_id: "scan=1".into(),
        ..msp(0)
    };
    functor
        .preprocess(&mut with_metadata, SpectraStPreprocessing::default())
        .unwrap();
    assert!(with_metadata.precursors.is_empty());
    assert!(with_metadata.native_id.is_empty());

    // A zero cut_peaks_below divides by zero upstream and discards everything.
    let mut untouched = msp(0);
    assert!(matches!(
        functor.preprocess(
            &mut untouched,
            SpectraStPreprocessing {
                cut_peaks_below: 0,
                ..SpectraStPreprocessing::default()
            },
        ),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(untouched.len(), 10);
}

#[test]
fn spectrast_delta_d_and_compute_f() {
    // START_SECTION(double delta_D(double top_hit, double runner_up))
    let functor = SpectraSTSimilarityScore::new().unwrap();
    // Class-test literals 0.2 and 0.96, and the DivisionByZero exception.
    close(functor.delta_d(5.0, 4.0).unwrap(), 0.2, 1e-12);
    close(functor.delta_d(25.0, 1.0).unwrap(), 0.96, 1e-12);
    assert!(matches!(
        functor.delta_d(0.0, 5.0),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        functor.delta_d(f64::NAN, 5.0),
        Err(Error::InvalidValue(_))
    ));

    // START_SECTION((double compute_F(double dot_product, double delta_D, double dot_bias)))
    // Upstream: NOT_TESTABLE, "pretty straightforward function". Derived from
    // the step function at SpectraSTSimilarityScore.cpp:132.
    close(functor.compute_f(1.0, 0.2, 0.2).unwrap(), 0.68, 1e-12);
    // Below 0.1 and on (0.35, 0.4] the penalty is 0.12.
    close(functor.compute_f(1.0, 0.2, 0.05).unwrap(), 0.56, 1e-12);
    close(functor.compute_f(1.0, 0.2, 0.4).unwrap(), 0.56, 1e-12);
    // The band [0.1, 0.35] is unpenalised, both bounds included.
    close(functor.compute_f(1.0, 0.2, 0.1).unwrap(), 0.68, 1e-12);
    close(functor.compute_f(1.0, 0.2, 0.35).unwrap(), 0.68, 1e-12);
    // (0.4, 0.45] costs 0.18 and anything above 0.45 costs 0.24.
    close(functor.compute_f(1.0, 0.2, 0.45).unwrap(), 0.5, 1e-12);
    close(functor.compute_f(1.0, 0.2, 0.46).unwrap(), 0.44, 1e-12);
    // The source lets NaN select b = 0 because every comparison is false.
    assert!(matches!(
        functor.compute_f(f64::NAN, 0.2, 0.2),
        Err(Error::InvalidValue(_))
    ));
}

// ---------------------------------------------------------------------------
// The base trait over all four derivatives
// ---------------------------------------------------------------------------

#[test]
fn all_four_derivatives_are_usable_through_the_base_trait() {
    let s1 = dfpianger();
    let functors: Vec<Box<dyn PeakSpectrumCompareFunctor>> = vec![
        Box::new(SpectrumPrecursorComparator::new().unwrap()),
        Box::new(SpectrumCheapDPCorr::new().unwrap()),
        Box::new(PeakAlignment::new().unwrap()),
        Box::new(SpectraSTSimilarityScore::new().unwrap()),
    ];
    let names: Vec<&str> = functors.iter().map(|f| f.name()).collect();
    assert_eq!(
        names,
        [
            "SpectrumPrecursorComparator",
            // PeakAlignment never renames its handler.
            "SpectrumCheapDPCorr",
            "PeakSpectrumCompareFunctor",
            "SpectraSTSimilarityScore",
        ]
    );
    for functor in &functors {
        // Every source derivative implements the one-spectrum overload as
        // operator()(spec, spec), which is the trait's default.
        assert_eq!(
            functor.self_score(&s1).unwrap(),
            functor.score(&s1, &s1).unwrap()
        );
        assert!(functor.self_score(&s1).unwrap() > 0.0);
    }
}
