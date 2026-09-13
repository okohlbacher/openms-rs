// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `COMPARISON/SpectrumAlignment.h`, `COMPARISON/SpectrumAlignmentScore.h`,
//! `COMPARISON/ZhangSimilarityScore.h` and
//! `COMPARISON/SteinScottImproveScore.h`: every `START_SECTION` of the four
//! upstream class tests, plus the boundaries those tests leave open.

use openms::comparison::{
    PeakSpectrumCompareFunctor, SpectrumAligner, SpectrumAlignmentScorer, SteinScottImproveScorer,
    Tolerance, ZhangSimilarityScorer,
};
use openms::format::dta;
use openms::param::{Param, ParamValue};
use openms::processing::{Normalizer, SpectrumFilter};
use openms::{Error, MSSpectrum, Peak1D};

fn spectrum(mzs: &[f64], intensities: &[f32]) -> MSSpectrum {
    MSSpectrum::from_peaks(
        mzs.iter()
            .zip(intensities)
            .map(|(&mz, &i)| Peak1D::new(mz, i))
            .collect(),
    )
}

fn unit(mzs: &[f64]) -> MSSpectrum {
    spectrum(mzs, &vec![1.0; mzs.len()])
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

/// `PILISSequenceDB_DFPIANGER_1.dta`, retained byte-identical.
fn dfpianger() -> MSSpectrum {
    dta::read(include_bytes!("data/comparison_dfpianger.dta").as_slice()).unwrap()
}

/// The upstream fixture after `Normalizer` with `method = to_one`, which is what
/// all three score class tests feed their functors.
fn normalized_dfpianger() -> MSSpectrum {
    let mut spectrum = dfpianger();
    Normalizer::default()
        .filter_spectrum(&mut spectrum)
        .unwrap();
    spectrum
}

/// `s2.resize(100)` in the class tests: the truncation happens *before* the
/// second normalisation, so the retained peaks keep the full spectrum's scale.
fn truncated(spectrum: &MSSpectrum, len: usize) -> MSSpectrum {
    let mut copy = spectrum.clone();
    copy.peaks.truncate(len);
    copy
}

fn with_float(handler_parameters: &Param, key: &str, value: f64) -> Param {
    let mut parameters = handler_parameters.clone();
    parameters
        .set_value(key, ParamValue::Float(value), "", &[])
        .unwrap();
    parameters
}

fn with_flag(handler_parameters: &Param, key: &str, value: bool) -> Param {
    let mut parameters = handler_parameters.clone();
    parameters
        .set_value(
            key,
            ParamValue::String(if value { "true" } else { "false" }.into()),
            "",
            &[],
        )
        .unwrap();
    parameters
}

// ---------------------------------------------------------------------------
// SpectrumAlignment_test.cpp - five sections
// ---------------------------------------------------------------------------

/// `START_SECTION(SpectrumAlignment())` and
/// `START_SECTION(virtual ~SpectrumAlignment())`.
#[test]
fn spectrum_aligner_constructs_with_source_defaults_and_drops() {
    let aligner = SpectrumAligner::new().unwrap();
    assert_eq!(aligner.name(), "SpectrumAlignment");
    let parameters = aligner.handler().parameters();
    assert_eq!(
        parameters.value("tolerance").unwrap().to_f64().unwrap(),
        0.3
    );
    assert!(
        !parameters
            .value("is_relative_tolerance")
            .unwrap()
            .to_bool()
            .unwrap()
    );
    assert_eq!(aligner.tolerance().unwrap(), Tolerance::Absolute(0.3));
    // The defaults tree keeps the source's restriction on the flag.
    assert_eq!(
        aligner
            .handler()
            .defaults()
            .valid_strings("is_relative_tolerance")
            .unwrap(),
        ["true".to_string(), "false".to_string()]
    );
    drop(aligner);
}

/// `START_SECTION(SpectrumAlignment(const SpectrumAlignment &source))` and
/// `START_SECTION(SpectrumAlignment& operator=(const SpectrumAlignment &source))`.
#[test]
fn spectrum_aligner_copy_and_assignment_carry_name_and_parameters() {
    let mut first = SpectrumAligner::new().unwrap();
    let adjusted = with_float(first.handler().parameters(), "tolerance", 0.2);
    first.handler_mut().set_parameters(&adjusted).unwrap();

    let copy = first.clone();
    assert_eq!(first.name(), copy.name());
    assert_eq!(first.handler().parameters(), copy.handler().parameters());
    assert!(first.handler().source_equal(copy.handler()).unwrap());
    assert_eq!(copy.tolerance().unwrap(), Tolerance::Absolute(0.2));

    let mut second = SpectrumAligner::new().unwrap();
    assert_eq!(second.tolerance().unwrap(), Tolerance::Absolute(0.3));
    second = first.clone();
    assert_eq!(first.name(), second.name());
    assert_eq!(&adjusted, second.handler().parameters());
}

/// `START_SECTION(template <typename SpectrumType> void getSpectrumAlignment(...))`,
/// the only section of this class test with assertions on data: 14 of them.
#[test]
fn spectrum_aligner_reproduces_the_upstream_alignment_golden() {
    let full = dfpianger();
    let aligner = SpectrumAligner::new().unwrap();
    // TEST_EQUAL(alignment.size(), s1.size())
    assert_eq!(
        aligner.spectrum_alignment(&full, &full).unwrap().len(),
        full.len()
    );
    assert_eq!(full.len(), 127);
    // TEST_EQUAL(alignment.size(), 100) after s2.resize(100)
    let head = truncated(&full, 100);
    assert_eq!(aligner.spectrum_alignment(&full, &head).unwrap().len(), 100);

    let s3 = dta::read(include_bytes!("data/comparison_alignment_1.dta").as_slice()).unwrap();
    let s4 = dta::read(include_bytes!("data/comparison_alignment_2.dta").as_slice()).unwrap();

    let mut banded = SpectrumAligner::new().unwrap();
    let parameters = with_float(banded.handler().parameters(), "tolerance", 1.01);
    banded.handler_mut().set_parameters(&parameters).unwrap();
    assert_eq!(
        banded.spectrum_alignment(&s3, &s4).unwrap(),
        [(0, 0), (1, 1), (3, 3), (4, 5), (6, 6)]
    );

    // p.setValue("is_relative_tolerance", "true"); p.setValue("tolerance", 10.0)
    let mut relative = SpectrumAligner::new().unwrap();
    let parameters = with_flag(
        relative.handler().parameters(),
        "is_relative_tolerance",
        true,
    );
    let parameters = with_float(&parameters, "tolerance", 10.0);
    relative.handler_mut().set_parameters(&parameters).unwrap();
    assert_eq!(relative.tolerance().unwrap(), Tolerance::Ppm(10.0));
    assert_eq!(relative.spectrum_alignment(&s3, &s4).unwrap(), [(6, 6)]);

    // one percent tolerance
    let parameters = with_float(relative.handler().parameters(), "tolerance", 1e4);
    relative.handler_mut().set_parameters(&parameters).unwrap();
    assert_eq!(
        relative.spectrum_alignment(&s3, &s4).unwrap(),
        [(0, 0), (1, 1), (2, 2), (3, 3), (4, 5), (5, 5), (6, 6)]
    );
}

// ---------------------------------------------------------------------------
// SpectrumAlignmentScore_test.cpp - six sections
// ---------------------------------------------------------------------------

/// `START_SECTION(SpectrumAlignmentScore())` and
/// `START_SECTION(virtual ~SpectrumAlignmentScore())`.
#[test]
fn spectrum_alignment_scorer_constructs_with_source_defaults_and_drops() {
    let scorer = SpectrumAlignmentScorer::new().unwrap();
    assert_eq!(scorer.name(), "SpectrumAlignmentScore");
    let parameters = scorer.handler().parameters();
    assert_eq!(
        parameters.value("tolerance").unwrap().to_f64().unwrap(),
        0.3
    );
    for flag in [
        "is_relative_tolerance",
        "use_linear_factor",
        "use_gaussian_factor",
    ] {
        assert!(!parameters.value(flag).unwrap().to_bool().unwrap());
    }
    assert_eq!(parameters.size(), 4);
    drop(scorer);
}

/// `START_SECTION(SpectrumAlignmentScore(const SpectrumAlignmentScore &source))`
/// and `START_SECTION(SpectrumAlignmentScore& operator=(...))`.
#[test]
fn spectrum_alignment_scorer_copy_and_assignment_carry_name_and_parameters() {
    let mut first = SpectrumAlignmentScorer::new().unwrap();
    let adjusted = with_float(first.handler().parameters(), "tolerance", 0.2);
    first.handler_mut().set_parameters(&adjusted).unwrap();

    let copy = first.clone();
    assert_eq!(first.name(), copy.name());
    assert_eq!(first.handler().parameters(), copy.handler().parameters());

    let mut second = SpectrumAlignmentScorer::new().unwrap();
    assert_eq!(
        second
            .handler()
            .parameters()
            .value("tolerance")
            .unwrap()
            .to_f64()
            .unwrap(),
        0.3
    );
    second = first.clone();
    assert_eq!(second.name(), "SpectrumAlignmentScore");
    assert_eq!(&adjusted, second.handler().parameters());
}

/// `START_SECTION(double operator () (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const)`,
/// the pairwise section: `TEST_REAL_SIMILAR(score, 1.48268)` and
/// `TEST_REAL_SIMILAR(score, 3.82472)` under `TOLERANCE_ABSOLUTE(0.01)`.
#[test]
fn spectrum_alignment_scorer_reproduces_the_upstream_pairwise_scores() {
    let normalized = normalized_dfpianger();
    let scorer = SpectrumAlignmentScorer::new().unwrap();
    let self_score = scorer.score(&normalized, &normalized).unwrap();
    close(self_score, 1.48268, 0.01);
    // Independent full-precision model of the source expression, including the
    // `float` intensity product and the `sqrt(sum1 * sum2)` denominator.
    assert_eq!(self_score, 1.484_501_010_820_065);

    let head = truncated(&normalized, 100);
    let truncated_score = scorer.score(&normalized, &head).unwrap();
    close(truncated_score, 3.82472, 1e-5);
    assert_eq!(truncated_score, 3.8247227487373805);
}

/// `START_SECTION(double operator()(const PeakSpectrum &spec) const)`:
/// `TEST_REAL_SIMILAR(score, 1.48268)`.
#[test]
fn spectrum_alignment_scorer_self_score_is_the_pairwise_score() {
    let normalized = normalized_dfpianger();
    let scorer = SpectrumAlignmentScorer::new().unwrap();
    let self_score = scorer.self_score(&normalized).unwrap();
    close(self_score, 1.48268, 0.01);
    assert_eq!(self_score, 1.484_501_010_820_065);
    // The source override is literally `operator()(spec, spec)`.
    assert_eq!(self_score, scorer.score(&normalized, &normalized).unwrap());
    // Self similarity above one: this normalisation is not a cosine.
    assert!(self_score > 1.0);
}

// ---------------------------------------------------------------------------
// ZhangSimilarityScore_test.cpp - six sections
// ---------------------------------------------------------------------------

/// `START_SECTION(ZhangSimilarityScore())` and
/// `START_SECTION(~ZhangSimilarityScore())`.
#[test]
fn zhang_scorer_constructs_with_source_defaults_and_drops() {
    let scorer = ZhangSimilarityScorer::new().unwrap();
    assert_eq!(scorer.name(), "ZhangSimilarityScore");
    let parameters = scorer.handler().parameters();
    assert_eq!(
        parameters.value("tolerance").unwrap().to_f64().unwrap(),
        0.2
    );
    for flag in [
        "is_relative_tolerance",
        "use_linear_factor",
        "use_gaussian_factor",
    ] {
        assert!(!parameters.value(flag).unwrap().to_bool().unwrap());
    }
    assert_eq!(parameters.size(), 4);
    drop(scorer);
}

/// `START_SECTION(ZhangSimilarityScore(const ZhangSimilarityScore& source))` and
/// `START_SECTION(ZhangSimilarityScore& operator = (const ZhangSimilarityScore& source))`.
#[test]
fn zhang_scorer_copy_and_assignment_carry_name_and_parameters() {
    let original = ZhangSimilarityScorer::new().unwrap();
    let copy = original.clone();
    assert_eq!(copy.name(), original.name());
    assert_eq!(copy.handler().parameters(), original.handler().parameters());

    let mut assigned = ZhangSimilarityScorer::new().unwrap();
    let scratch = with_float(assigned.handler().parameters(), "tolerance", 0.9);
    assigned.handler_mut().set_parameters(&scratch).unwrap();
    assigned = original.clone();
    assert_eq!(assigned.name(), original.name());
    assert_eq!(
        assigned
            .handler()
            .parameters()
            .value("tolerance")
            .unwrap()
            .to_f64()
            .unwrap(),
        0.2
    );
}

/// `START_SECTION(double operator () (const PeakSpectrum& spec) const)`:
/// `TEST_REAL_SIMILAR(score, 1.82682)`.
#[test]
fn zhang_scorer_reproduces_the_upstream_self_score() {
    let normalized = normalized_dfpianger();
    let scorer = ZhangSimilarityScorer::new().unwrap();
    let score = scorer.self_score(&normalized).unwrap();
    close(score, 1.82682, 1e-5);
    assert_eq!(score, 1.8268153570547176);
    assert_eq!(score, scorer.score(&normalized, &normalized).unwrap());
}

/// `START_SECTION(double operator () (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const)`:
/// `TEST_REAL_SIMILAR(score, 1.82682)` and `TEST_REAL_SIMILAR(score, 0.328749)`.
#[test]
fn zhang_scorer_reproduces_the_upstream_pairwise_scores() {
    let normalized = normalized_dfpianger();
    let scorer = ZhangSimilarityScorer::new().unwrap();
    close(
        scorer.score(&normalized, &normalized).unwrap(),
        1.82682,
        1e-5,
    );

    let head = truncated(&normalized, 100);
    let truncated_score = scorer.score(&normalized, &head).unwrap();
    close(truncated_score, 0.328749, 1e-6);
    assert_eq!(truncated_score, 0.32874865683513527);
}

// ---------------------------------------------------------------------------
// SteinScottImproveScore_test.cpp - six sections
// ---------------------------------------------------------------------------

/// The class test's fixture: five peaks at 500..900 whose intensity equals their
/// m/z, so every product is exactly representable and the result has a closed
/// form.
fn stein_scott_fixture() -> MSSpectrum {
    spectrum(
        &[500.0, 600.0, 700.0, 800.0, 900.0],
        &[500.0, 600.0, 700.0, 800.0, 900.0],
    )
}

/// `START_SECTION(SteinScottImproveScore())` and
/// `START_SECTION(virtual ~SteinScottImproveScore())`.
#[test]
fn stein_scott_scorer_constructs_with_source_defaults_and_drops() {
    let scorer = SteinScottImproveScorer::new().unwrap();
    assert_eq!(scorer.name(), "SteinScottImproveScore");
    let parameters = scorer.handler().parameters();
    assert_eq!(
        parameters.value("tolerance").unwrap().to_f64().unwrap(),
        0.2
    );
    assert_eq!(
        parameters.value("threshold").unwrap().to_f64().unwrap(),
        0.2
    );
    assert_eq!(parameters.size(), 2);
    // Neither parameter carries a restriction upstream.
    assert!(parameters.valid_strings("tolerance").is_err());
    drop(scorer);
}

/// `START_SECTION(SteinScottImproveScore(const SteinScottImproveScore& source))`
/// and `START_SECTION(SteinScottImproveScore& operator = (...))`.
#[test]
fn stein_scott_scorer_copy_and_assignment_carry_name_and_parameters() {
    let original = SteinScottImproveScorer::new().unwrap();
    let copy = original.clone();
    assert_eq!(copy.name(), original.name());
    assert_eq!(copy.handler().parameters(), original.handler().parameters());

    let mut assigned = SteinScottImproveScorer::new().unwrap();
    let scratch = with_float(assigned.handler().parameters(), "threshold", 0.9);
    assigned.handler_mut().set_parameters(&scratch).unwrap();
    assigned = original.clone();
    assert_eq!(assigned.name(), original.name());
    assert_eq!(
        assigned.handler().parameters(),
        original.handler().parameters()
    );
}

/// `START_SECTION(double operator () (const PeakSpectrum& spec) const)`:
/// `if (score > 0.99) score = 1; TEST_REAL_SIMILAR(score, 1)`.
#[test]
fn stein_scott_scorer_reproduces_the_upstream_self_score() {
    let fixture = stein_scott_fixture();
    let scorer = SteinScottImproveScorer::new().unwrap();
    let score = scorer.self_score(&fixture).unwrap();
    assert!(score > 0.99);
    // Derived rather than transcribed: every pair inside the +-0.4 window is a
    // peak with itself, so `sum` is the exact sum of squares 2550000, `sum1` and
    // `sum2` are the same, and `z = 0.2 / 10000 * 3500^2`.
    let z = (0.2 / 10000.0) * (3500.0 * 3500.0);
    assert_eq!(score, (2550000.0 - z) / 2550000.0);
    assert_eq!(score, 0.9999039215686274);
    assert_eq!(score, scorer.score(&fixture, &fixture).unwrap());
}

/// `START_SECTION(double operator () (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const)`:
/// the same fixture built twice.
#[test]
fn stein_scott_scorer_reproduces_the_upstream_pairwise_score() {
    let first = stein_scott_fixture();
    let second = stein_scott_fixture();
    let scorer = SteinScottImproveScorer::new().unwrap();
    let score = scorer.score(&first, &second).unwrap();
    assert!(score > 0.99);
    assert_eq!(score, 0.9999039215686274);
}

// ---------------------------------------------------------------------------
// Boundaries the class tests leave open
// ---------------------------------------------------------------------------

/// Every degenerate denominator the source divides by zero on.
#[test]
fn empty_and_zero_intensity_spectra_score_a_defined_zero() {
    let empty = MSSpectrum::default();
    let populated = unit(&[100.0, 200.0]);
    let silent = spectrum(&[100.0, 200.0], &[0.0, 0.0]);

    let alignment = SpectrumAlignmentScorer::new().unwrap();
    let zhang = ZhangSimilarityScorer::new().unwrap();
    let stein = SteinScottImproveScorer::new().unwrap();

    for (left, right) in [
        (&empty, &empty),
        (&empty, &populated),
        (&populated, &empty),
        (&silent, &populated),
    ] {
        assert_eq!(alignment.score(left, right).unwrap(), 0.0);
        assert_eq!(zhang.score(left, right).unwrap(), 0.0);
        assert_eq!(stein.score(left, right).unwrap(), 0.0);
    }
    assert_eq!(alignment.self_score(&empty).unwrap(), 0.0);
    assert_eq!(zhang.self_score(&empty).unwrap(), 0.0);
    assert_eq!(stein.self_score(&empty).unwrap(), 0.0);

    // An empty alignment between two populated spectra is not degenerate: the
    // source scores it too, and only Stein/Scott's expected-overlap term is
    // nonzero there.
    let low = unit(&[100.0]);
    let high = unit(&[500.0]);
    assert_eq!(
        SpectrumAligner::new()
            .unwrap()
            .spectrum_alignment(&low, &high)
            .unwrap(),
        []
    );
    assert_eq!(alignment.score(&low, &high).unwrap(), 0.0);
    assert_eq!(zhang.score(&low, &high).unwrap(), 0.0);
    // (0 - 0.2/10000 * 1 * 1) / 1 is negative, so the 0.2 threshold zeroes it.
    assert_eq!(stein.score(&low, &high).unwrap(), 0.0);
    let mut permissive = SteinScottImproveScorer::new().unwrap();
    let parameters = with_float(permissive.handler().parameters(), "threshold", -1.0);
    permissive
        .handler_mut()
        .set_parameters(&parameters)
        .unwrap();
    assert_eq!(
        permissive.score(&low, &high).unwrap(),
        -((0.2 / 10000.0) * (1.0 * 1.0))
    );
}

/// The upstream `Exception::NotImplemented` on `is_relative_tolerance`.
#[test]
fn zhang_refuses_the_unimplemented_relative_tolerance() {
    let mut scorer = ZhangSimilarityScorer::new().unwrap();
    let parameters = with_flag(scorer.handler().parameters(), "is_relative_tolerance", true);
    scorer.handler_mut().set_parameters(&parameters).unwrap();
    let peaks = unit(&[100.0]);
    assert!(matches!(
        scorer.score(&peaks, &peaks),
        Err(Error::Unsupported(_))
    ));
    // The parameter is still registered, so a Param tree round-trips.
    assert!(
        scorer
            .handler()
            .parameters()
            .value("is_relative_tolerance")
            .unwrap()
            .to_bool()
            .unwrap()
    );
}

/// `use_linear_factor` and `use_gaussian_factor`, in both scorers that have them.
#[test]
fn weighting_factors_reproduce_the_source_expressions() {
    let left = unit(&[10.0]);
    let right = unit(&[10.5]);

    let mut scorer = SpectrumAlignmentScorer::new().unwrap();
    let wide = with_float(scorer.handler().parameters(), "tolerance", 1.0);
    scorer.handler_mut().set_parameters(&wide).unwrap();
    assert_eq!(scorer.score(&left, &right).unwrap(), 1.0_f64);

    let linear = with_flag(scorer.handler().parameters(), "use_linear_factor", true);
    scorer.handler_mut().set_parameters(&linear).unwrap();
    // factor = (1.0 - 0.5) / 1.0, then sqrt(1 * 1 * factor) / sqrt(1 * 1).
    assert_eq!(scorer.score(&left, &right).unwrap(), 0.5_f64.sqrt());

    let gaussian = with_flag(
        &with_flag(scorer.handler().parameters(), "use_linear_factor", false),
        "use_gaussian_factor",
        true,
    );
    scorer.handler_mut().set_parameters(&gaussian).unwrap();
    close(
        scorer.score(&left, &right).unwrap(),
        0.8676323347781927_f64.sqrt(),
        1e-15,
    );

    // ZhangSimilarityScore::getFactor_ is public here and has the same two arms.
    assert_eq!(
        ZhangSimilarityScorer::factor(1.0, 0.5, false).unwrap(),
        0.5_f64
    );
    close(
        ZhangSimilarityScorer::factor(1.0, 0.5, true).unwrap(),
        0.8676323347781927,
        1e-15,
    );
    assert!(ZhangSimilarityScorer::factor(0.0, 0.0, false).is_err());
    assert!(ZhangSimilarityScorer::factor(f64::NAN, 0.0, true).is_err());

    let mut zhang = ZhangSimilarityScorer::new().unwrap();
    let wide = with_flag(
        &with_float(zhang.handler().parameters(), "tolerance", 1.0),
        "use_linear_factor",
        true,
    );
    zhang.handler_mut().set_parameters(&wide).unwrap();
    assert_eq!(zhang.score(&left, &right).unwrap(), 0.5_f64.sqrt());
}

/// Both weighting flags at once: an upstream debug-only precondition in one
/// scorer and silently accepted in the other; refused in both here.
#[test]
fn both_weighting_flags_at_once_are_refused() {
    let peaks = unit(&[10.0]);
    let mut scorer = SpectrumAlignmentScorer::new().unwrap();
    let both = with_flag(
        &with_flag(scorer.handler().parameters(), "use_linear_factor", true),
        "use_gaussian_factor",
        true,
    );
    scorer.handler_mut().set_parameters(&both).unwrap();
    assert!(matches!(
        scorer.score(&peaks, &peaks),
        Err(Error::InvalidValue(_))
    ));

    let mut zhang = ZhangSimilarityScorer::new().unwrap();
    let both = with_flag(
        &with_flag(zhang.handler().parameters(), "use_linear_factor", true),
        "use_gaussian_factor",
        true,
    );
    zhang.handler_mut().set_parameters(&both).unwrap();
    assert!(matches!(
        zhang.score(&peaks, &peaks),
        Err(Error::InvalidValue(_))
    ));
}

/// `getFactor_` caches its Gaussian denominator in a function-local `static`,
/// so upstream the first call in a process fixes the scale for every later one.
#[test]
fn zhang_gaussian_scale_follows_the_instance_not_the_first_call() {
    let left = unit(&[1.0]);
    let right = unit(&[1.1]);
    let build = |tolerance: f64| {
        let mut scorer = ZhangSimilarityScorer::new().unwrap();
        let parameters = with_flag(
            &with_float(scorer.handler().parameters(), "tolerance", tolerance),
            "use_gaussian_factor",
            true,
        );
        scorer.handler_mut().set_parameters(&parameters).unwrap();
        scorer
    };
    let narrow = build(0.2);
    let broad = build(2.0);
    let narrow_score = narrow.score(&left, &right).unwrap();
    let broad_score = broad.score(&left, &right).unwrap();
    assert!(narrow_score < broad_score && broad_score <= 1.0);
    // Repeating in the other order changes nothing; a cached denominator would.
    assert_eq!(broad_score, broad.score(&left, &right).unwrap());
    assert_eq!(narrow_score, narrow.score(&left, &right).unwrap());
    assert_eq!(
        narrow_score,
        ZhangSimilarityScorer::factor(0.2, 0.1_f64, true)
            .unwrap()
            .sqrt()
    );
}

/// The ppm window that selects a pair is not the ppm window that weights it.
#[test]
fn relative_tolerance_weighting_reports_the_source_nan_instead() {
    let left = unit(&[1000.0]);
    let right = unit(&[1000.01000000002]);
    let mut scorer = SpectrumAlignmentScorer::new().unwrap();
    let relative = with_float(
        &with_flag(scorer.handler().parameters(), "is_relative_tolerance", true),
        "tolerance",
        10.0,
    );
    scorer.handler_mut().set_parameters(&relative).unwrap();
    // MatchedIterator narrows to f32 and matches; the score's own f64 window is
    // 10.0 * 1000.0 * 1e-6 = 0.01, which the 0.01000000002 difference exceeds.
    let mut aligner = SpectrumAligner::new().unwrap();
    let alignment_parameters = with_float(
        &with_flag(
            aligner.handler().parameters(),
            "is_relative_tolerance",
            true,
        ),
        "tolerance",
        10.0,
    );
    aligner
        .handler_mut()
        .set_parameters(&alignment_parameters)
        .unwrap();
    assert_eq!(aligner.spectrum_alignment(&left, &right).unwrap(), [(0, 0)]);
    assert_eq!(scorer.score(&left, &right).unwrap(), 1.0);
    for flag in ["use_linear_factor", "use_gaussian_factor"] {
        let mut weighted = scorer.clone();
        let parameters = with_flag(weighted.handler().parameters(), flag, true);
        weighted.handler_mut().set_parameters(&parameters).unwrap();
        // Linear: a negative radicand. Gaussian: erfc stays finite, so only the
        // linear arm can fail here.
        let outcome = weighted.score(&left, &right);
        if flag == "use_linear_factor" {
            assert!(matches!(outcome, Err(Error::InvalidValue(_))));
        } else {
            assert!(outcome.unwrap() < 1.0);
        }
    }
}

/// The many-to-many walk, its window boundaries and its cursor.
#[test]
fn pair_walk_is_many_to_many_with_source_boundary_rules() {
    let mut zhang = ZhangSimilarityScorer::new().unwrap();
    let one = with_float(zhang.handler().parameters(), "tolerance", 1.0);
    zhang.handler_mut().set_parameters(&one).unwrap();
    // Zhang's window is strict: |1.0 - 2.0| < 1.0 is false.
    assert_eq!(zhang.score(&unit(&[1.0]), &unit(&[2.0])).unwrap(), 0.0);
    // ... and everything strictly inside counts in full, because the default
    // factor is 1.0 and the distance never weights the product.
    assert_eq!(
        zhang.score(&unit(&[1.0]), &unit(&[1.999999999])).unwrap(),
        1.0
    );
    // With the linear factor the same pair is weighted by (1.0 - d) / 1.0.
    let weighted = with_flag(zhang.handler().parameters(), "use_linear_factor", true);
    let mut linear = zhang.clone();
    linear.handler_mut().set_parameters(&weighted).unwrap();
    // The difference is formed exactly as the source forms it, so the expected
    // value has to carry the same cancellation rather than a decimal literal.
    let difference = (1.0_f64 - 1.999999999).abs();
    assert_eq!(
        linear.score(&unit(&[1.0]), &unit(&[1.999999999])).unwrap(),
        ((1.0 - difference) / 1.0_f64).sqrt()
    );

    // Stein/Scott's window is 2 * tolerance and its boundary is inclusive.
    let mut stein = SteinScottImproveScorer::new().unwrap();
    let parameters = with_float(
        &with_float(stein.handler().parameters(), "tolerance", 0.5),
        "threshold",
        0.0,
    );
    stein.handler_mut().set_parameters(&parameters).unwrap();
    close(
        stein.score(&unit(&[1.0]), &unit(&[2.0])).unwrap(),
        0.99995,
        1e-15,
    );

    // One reference peak can take several targets, which an alignment forbids.
    let repeated = unit(&[1.0, 1.0, 1.0]);
    let single = unit(&[1.0]);
    close(
        zhang.score(&single, &repeated).unwrap(),
        3.0 / 3.0_f64.sqrt(),
        1e-15,
    );
    assert_eq!(
        SpectrumAligner::new()
            .unwrap()
            .spectrum_alignment(&single, &repeated)
            .unwrap()
            .len(),
        1
    );
}

/// Resource ceilings, which the source does not have.
#[test]
fn resource_ceilings_refuse_unbounded_work() {
    let repeated = unit(&[1.0, 1.0, 1.0]);
    let mut zhang = ZhangSimilarityScorer::new().unwrap();
    zhang.max_pairs = 2;
    assert!(zhang.score(&repeated, &repeated).is_err());
    zhang.max_pairs = 0;
    assert!(zhang.score(&repeated, &repeated).is_err());

    let mut stein = SteinScottImproveScorer::new().unwrap();
    stein.max_pairs = 2;
    assert!(stein.score(&repeated, &repeated).is_err());

    let two = unit(&[1.0, 2.0]);
    let mut aligner = SpectrumAligner::new().unwrap();
    aligner.max_cells = 4;
    assert!(aligner.spectrum_alignment(&two, &two).is_err());
    aligner.max_cells = 5_000_000;
    assert_eq!(aligner.spectrum_alignment(&two, &two).unwrap().len(), 2);

    let mut scorer = SpectrumAlignmentScorer::new().unwrap();
    scorer.max_cells = 4;
    assert!(scorer.score(&two, &two).is_err());
}

/// Input the source accepts silently and the port refuses.
#[test]
fn invalid_input_is_refused_rather_than_scored() {
    let descending = MSSpectrum::from_peaks(vec![Peak1D::new(2.0, 1.0), Peak1D::new(1.0, 1.0)]);
    let ascending = unit(&[1.0, 2.0]);
    let negative = spectrum(&[1.0], &[-1.0]);

    let alignment = SpectrumAlignmentScorer::new().unwrap();
    let zhang = ZhangSimilarityScorer::new().unwrap();
    let stein = SteinScottImproveScorer::new().unwrap();
    let aligner = SpectrumAligner::new().unwrap();

    assert!(matches!(
        aligner.spectrum_alignment(&descending, &ascending),
        Err(Error::UnsortedData)
    ));
    for outcome in [
        alignment.score(&descending, &ascending),
        zhang.score(&descending, &ascending),
        stein.score(&descending, &ascending),
    ] {
        assert!(matches!(outcome, Err(Error::UnsortedData)));
    }
    for outcome in [
        alignment.score(&negative, &ascending),
        zhang.score(&negative, &ascending),
        stein.score(&negative, &ascending),
    ] {
        assert!(matches!(outcome, Err(Error::InvalidValue(_))));
    }

    // A negative tolerance produces an empty alignment upstream; here the
    // alignment primitive rejects it outright.
    let mut backwards = SpectrumAligner::new().unwrap();
    let parameters = with_float(backwards.handler().parameters(), "tolerance", -1.0);
    backwards.handler_mut().set_parameters(&parameters).unwrap();
    assert!(
        backwards
            .spectrum_alignment(&ascending, &ascending)
            .is_err()
    );

    // The valid-string restriction survives into set_parameters.
    let mut flagged = ZhangSimilarityScorer::new().unwrap();
    let mut wrong = flagged.handler().parameters().clone();
    wrong
        .set_value(
            "use_linear_factor",
            ParamValue::String("yes".into()),
            "",
            &[],
        )
        .unwrap();
    assert!(flagged.handler_mut().set_parameters(&wrong).is_err());
}

/// The base class exists to be held by pointer; the trait replaces that.
#[test]
fn the_three_scorers_are_usable_as_trait_objects() {
    let normalized = normalized_dfpianger();
    let functors: Vec<Box<dyn PeakSpectrumCompareFunctor>> = vec![
        Box::new(SpectrumAlignmentScorer::new().unwrap()),
        Box::new(ZhangSimilarityScorer::new().unwrap()),
        Box::new(SteinScottImproveScorer::new().unwrap()),
    ];
    let names: Vec<&str> = functors.iter().map(|f| f.name()).collect();
    assert_eq!(
        names,
        [
            "SpectrumAlignmentScore",
            "ZhangSimilarityScore",
            "SteinScottImproveScore"
        ]
    );
    for functor in &functors {
        let self_score = functor.self_score(&normalized).unwrap();
        assert!(
            self_score >= 0.0,
            "{} returned {self_score}",
            functor.name()
        );
        assert_eq!(
            self_score,
            functor.score(&normalized, &normalized).unwrap(),
            "{}",
            functor.name()
        );
    }
}

/// The score is symmetric for these three functors on the upstream fixture, an
/// invariant no class test checks.
#[test]
fn scores_are_symmetric_on_the_upstream_fixture() {
    let normalized = normalized_dfpianger();
    let head = truncated(&normalized, 100);
    let alignment = SpectrumAlignmentScorer::new().unwrap();
    let zhang = ZhangSimilarityScorer::new().unwrap();
    let stein = SteinScottImproveScorer::new().unwrap();
    close(
        alignment.score(&normalized, &head).unwrap(),
        alignment.score(&head, &normalized).unwrap(),
        1e-12,
    );
    assert_eq!(
        zhang.score(&normalized, &head).unwrap(),
        zhang.score(&head, &normalized).unwrap()
    );
    assert_eq!(
        stein.score(&normalized, &head).unwrap(),
        stein.score(&head, &normalized).unwrap()
    );
}
