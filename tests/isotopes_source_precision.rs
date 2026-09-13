// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Source-precision coarse isotope patterns and the source `trimLeft`.
//!
//! Tier 1: `tests/data/isotopes_source_precision/probe.tsv` was printed by the
//! executed C++ SDK probe `../oracle/b2-iso-source-precision` (product SDK,
//! core 4fdec46; hashes in `tests/data/isotopes_source_precision_provenance.json`).
//! Masses are binary64 and `Peak1D` intensities binary32; both are compared bit
//! for bit. Tier 3: `CoarseIsotopeDistribution_test.cpp` and
//! `IsotopeDistribution_test.cpp` literals at `bc9cc12`, compared with the
//! ClassTest `TEST_REAL_SIMILAR` rule. Tier 4: native errors and limits.

use openms::chemistry::isotopes::{
    AveragineComposition, CoarseIsotopePatternGenerator as Generator, CoarseMassMode,
    IsotopeDistribution, IsotopePeak, ProbabilityPrecision,
};
use openms::chemistry::{EmpiricalFormula, element};

const PROBE: &str = include_str!("data/isotopes_source_precision/probe.tsv");
const SINGLE: ProbabilityPrecision = ProbabilityPrecision::SourceSingle;
const APPROXIMATE: CoarseMassMode = CoarseMassMode::Approximate;
const NOMINAL: CoarseMassMode = CoarseMassMode::Nominal;
/// FeatureFinderAlgorithmPicked `mass_window_width` and the FFC_1
/// `intensity_percentage_optional` (0.1 %).
const MASS_WINDOW_WIDTH: f64 = 100.0;
const OPTIONAL_CUTOFF: f64 = 0.1 / 100.0;

fn formula(text: &str) -> EmpiricalFormula {
    text.parse().unwrap()
}

fn source(max: Option<usize>, mode: CoarseMassMode) -> Generator {
    Generator::new(max, mode).unwrap().with_precision(SINGLE)
}

fn distribution(peaks: &[(f64, f64)]) -> IsotopeDistribution {
    IsotopeDistribution::from_peaks(
        peaks
            .iter()
            .map(|&(mass, probability)| IsotopePeak { mass, probability })
            .collect(),
    )
    .unwrap()
}

fn probe_line(kind: &str, label: &str) -> Vec<&'static str> {
    PROBE
        .lines()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .find(|fields| fields[0] == kind && fields[1] == label)
        .unwrap_or_else(|| panic!("probe line {kind} {label}"))
}

/// `(mass bits, Peak1D intensity bits)` of a probe distribution line.
fn probe_peaks(kind: &str, label: &str) -> Vec<(u64, u32)> {
    let fields = probe_line(kind, label);
    let count: usize = fields[2].parse().unwrap();
    let peaks: Vec<(u64, u32)> = fields[3..]
        .iter()
        .map(|peak| {
            let (mass, intensity) = peak.split_once(':').unwrap();
            (
                u64::from_str_radix(mass, 16).unwrap(),
                u32::from_str_radix(intensity, 16).unwrap(),
            )
        })
        .collect();
    assert_eq!(peaks.len(), count, "{kind} {label}");
    peaks
}

/// Bits of a source-precision result; every weight must be a binary32 value.
fn bits(distribution: &IsotopeDistribution) -> Vec<(u64, u32)> {
    distribution
        .peaks()
        .iter()
        .map(|peak| {
            let single = peak.probability as f32;
            assert_eq!(
                f64::from(single).to_bits(),
                peak.probability.to_bits(),
                "source precision stores exact binary32 weights"
            );
            (peak.mass.to_bits(), single.to_bits())
        })
        .collect()
}

/// ClassTest `isRealSimilar` (`ClassTest.cpp:364-490`): absolute difference
/// at most 1e-5, or ratio at most 1 + 1e-5, with the source's zero cases.
fn real_similar(actual: f64, expected: f64) -> bool {
    if actual.is_nan() || expected.is_nan() {
        return false;
    }
    let small = (actual - expected).abs() <= 1e-5;
    if actual == 0.0 {
        return expected == 0.0 || small;
    }
    if expected == 0.0 {
        return small;
    }
    let mut ratio = actual / expected;
    if ratio < 0.0 {
        return small;
    }
    if ratio < 1.0 {
        ratio = 1.0 / ratio;
    }
    ratio <= 1.0 + 1e-5 || small
}

fn similar(actual: f64, expected: f64, what: &str) {
    assert!(
        real_similar(actual, expected),
        "{what}: got {actual:.17}, expected {expected}"
    );
}

#[test]
fn element_tables_narrow_to_the_source_peak1d_values() {
    for symbol in ["H", "C", "N", "O", "S", "P", "Br", "Se"] {
        let native: Vec<(u64, u32)> = element(symbol)
            .unwrap()
            .isotopes()
            .iter()
            .map(|isotope| (isotope.mass.to_bits(), (isotope.abundance as f32).to_bits()))
            .collect();
        assert_eq!(native, probe_peaks("element", symbol), "{symbol}");
    }
}

#[test]
fn run_matches_executed_cpp_bits() {
    for (label, text, max, mode) in [
        ("C6H12O6_max3", "C6H12O6", Some(3), APPROXIMATE),
        ("C6H12O6+2_max3", "C6H12O6+2", Some(3), APPROXIMATE),
        ("C6H14O6_max3", "C6H14O6", Some(3), APPROXIMATE),
        ("C222N190O110_max11", "C222N190O110", Some(11), APPROXIMATE),
        ("Br2_max5", "Br2", Some(5), APPROXIMATE),
        ("CBr2_max7", "CBr2", Some(7), APPROXIMATE),
        ("C160_max10", "C160", Some(10), APPROXIMATE),
        ("H2_max11", "H2", Some(11), APPROXIMATE),
        ("C100_max11_round", "C100", Some(11), NOMINAL),
        ("C2_max0", "C2", None, APPROXIMATE),
    ] {
        let pattern = source(max, mode).run(&formula(text)).unwrap();
        assert_eq!(bits(&pattern), probe_peaks("run", label), "{label}");
    }
}

#[test]
fn averagine_estimates_match_executed_cpp_bits() {
    let generator = source(Some(3), APPROXIMATE);
    for mass in [100.0, 1000.0, 10000.0] {
        let cases = [
            ("peptide", generator.estimate_from_peptide_weight(mass)),
            ("rna", generator.estimate_from_rna_weight(mass)),
            ("dna", generator.estimate_from_dna_weight(mass)),
        ];
        for (family, pattern) in cases {
            let label = format!("{family}_max3_{mass}");
            assert_eq!(
                bits(&pattern.unwrap()),
                probe_peaks("estimate", &label),
                "{label}"
            );
        }
    }
    let unbounded = source(None, APPROXIMATE)
        .estimate_from_peptide_weight(1234.2)
        .unwrap();
    assert_eq!(unbounded.len(), 317); // CoarseIsotopeDistribution_test.cpp:101
    assert_eq!(
        bits(&unbounded),
        probe_peaks("estimate", "peptide_max0_1234.2")
    );
}

#[test]
fn feature_finder_picked_windows_match_executed_cpp_bits_and_trims() {
    // FeatureFinderAlgorithmPicked.cpp:364-374: max_isotopes 20, window centre
    // 0.5 * width + index * width, source trimLeft then trimRight.
    for index in 0..=80_usize {
        let label = format!("ffap20_w{index}");
        let mass = 0.5 * MASS_WINDOW_WIDTH + index as f64 * MASS_WINDOW_WIDTH;
        let mut pattern = source(Some(20), APPROXIMATE)
            .estimate_from_peptide_weight(mass)
            .unwrap();
        assert_eq!(bits(&pattern), probe_peaks("window", &label), "{label}");
        let trims = probe_line("window_trim", &label);
        let before = pattern.len();
        pattern.trim_left_source(OPTIONAL_CUTOFF).unwrap();
        let after_left = pattern.len();
        pattern.trim_right(OPTIONAL_CUTOFF).unwrap();
        assert_eq!(
            [before, after_left, pattern.len()].map(|n| n.to_string()),
            [trims[2], trims[3], trims[4]].map(str::to_string),
            "{label}"
        );
    }
}

#[test]
fn override_weights_narrow_like_source_insertion_and_match_cpp_windows() {
    // FeatureFinderAlgorithmPicked.cpp:163-170 with abundance_12C = 90, as the
    // intended two-isotope distribution (`set`), and max_isotopes 20 + 1000.
    let abundance_12c = 90.0;
    let carbon = distribution(&[
        (12.0, abundance_12c / 100.0),
        (13.0, 1.0 - (abundance_12c / 100.0)),
    ]);
    let narrowed: Vec<(u64, u32)> = carbon
        .peaks()
        .iter()
        .map(|peak| (peak.mass.to_bits(), (peak.probability as f32).to_bits()))
        .collect();
    assert_eq!(narrowed, probe_peaks("override", "set_12C_90"));
    let mut generator = source(Some(1020), APPROXIMATE);
    generator.set_isotope_override("C", carbon).unwrap();
    for index in [0_usize, 1, 5, 10] {
        let mass = 0.5 * 100.0 + index as f64 * 100.0;
        let label = format!("set_12C_90_max1020_w{index}");
        let pattern = generator.estimate_from_peptide_weight(mass).unwrap();
        assert_eq!(
            bits(&pattern),
            probe_peaks("override_window", &label),
            "{label}"
        );
    }
    // The source builds that override by inserting into a default-constructed
    // IsotopeDistribution, which already holds (0, 1). The executed C++ keeps
    // the stray peak and yields different, longer patterns; the native
    // override validation rejects the input instead of reproducing it.
    assert_eq!(
        probe_peaks("override", "ffap_insert_12C_90")[0],
        (0.0_f64.to_bits(), 1.0_f32.to_bits())
    );
    for index in [0_usize, 1, 5, 10] {
        let stray = probe_peaks(
            "override_window",
            &format!("ffap_insert_12C_90_max1020_w{index}"),
        );
        let intended = probe_peaks("override_window", &format!("set_12C_90_max1020_w{index}"));
        assert!(stray.len() > intended.len());
    }
    let stray = distribution(&[(0.0, 1.0), (12.0, 0.9), (13.0, 0.1)]);
    assert!(
        source(Some(1020), APPROXIMATE)
            .set_isotope_override("C", stray)
            .is_err()
    );
}

#[test]
fn raw_convolution_matches_executed_cpp_bits() {
    let carbon = IsotopeDistribution::from_peaks(
        element("C")
            .unwrap()
            .isotopes()
            .iter()
            .map(|isotope| IsotopePeak {
                mass: isotope.mass,
                probability: isotope.abundance,
            })
            .collect(),
    )
    .unwrap();
    let squared = source(None, APPROXIMATE)
        .convolve(&carbon, &carbon)
        .unwrap();
    assert_eq!(bits(&squared), probe_peaks("convolve", "C_C_max0"));
    assert_eq!(
        source(None, APPROXIMATE)
            .convolve_power(&carbon, 2)
            .unwrap(),
        squared
    );
    // CoarseIsotopeDistribution_test.cpp:111-120.
    let identity = source(Some(1), APPROXIMATE)
        .convolve(
            &IsotopeDistribution::default(),
            &IsotopeDistribution::default(),
        )
        .unwrap();
    assert_eq!(
        bits(&identity),
        probe_peaks("convolve", "identity_identity_max1")
    );
    assert_eq!(identity.len(), 1);
    assert_eq!(identity.peaks()[0].mass, 0.0);
    assert_eq!(identity.peaks()[0].probability, 1.0);
}

#[test]
fn source_trim_left_keeps_every_peak_when_none_reaches_the_cutoff() {
    let c160 = source(Some(10), APPROXIMATE).run(&formula("C160")).unwrap();
    let executed = probe_line("trim", "C160_max10");
    let mut all_below_left = c160.clone();
    all_below_left.trim_left_source(0.9).unwrap();
    let mut all_below_right = c160.clone();
    all_below_right.trim_right(0.9).unwrap();
    let mut both = c160.clone();
    both.trim_right(0.2).unwrap();
    both.trim_left_source(0.2).unwrap();
    assert_eq!(
        [
            c160.len(),
            all_below_left.len(),
            all_below_right.len(),
            both.len()
        ]
        .map(|n| n.to_string()),
        [executed[2], executed[3], executed[4], executed[5]].map(str::to_string)
    );
    assert_eq!(all_below_left, c160, "source trimLeft erases nothing");
    let mut native = c160.clone();
    native.trim_left(0.9).unwrap();
    assert!(
        native.is_empty(),
        "the native default still removes every peak"
    );

    let mut empty = IsotopeDistribution::empty();
    empty.trim_left_source(0.5).unwrap();
    assert_eq!(empty.len().to_string(), probe_line("trim", "empty")[2]);

    let mut renormalized = both.clone();
    renormalized.renormalize_with(SINGLE).unwrap();
    assert_eq!(
        bits(&renormalized),
        probe_peaks("renormalize", "C160_max10_trim0.2")
    );

    // Equality with the cutoff is retained; invalid cutoffs leave data intact.
    let mut gapped = distribution(&[(0.0, 0.01), (1.0, 0.2), (2.0, 0.0), (3.0, 0.1)]);
    gapped.trim_left_source(0.2).unwrap();
    assert_eq!(gapped.peaks()[0].mass, 1.0);
    let old = gapped.clone();
    assert!(gapped.trim_left_source(f64::NAN).is_err());
    assert!(gapped.trim_left_source(-1.0).is_err());
    assert_eq!(gapped, old);
}

#[test]
fn fragment_weights_match_executed_cpp_bits() {
    let generator = source(None, APPROXIMATE);
    for (label, precursor, fragment) in [
        ("peptide_200_100_iso01", 200.0, 100.0),
        ("peptide_2000_1000_iso01", 2000.0, 1000.0),
    ] {
        let joint = generator
            .estimate_fragment_from_weights(
                precursor,
                fragment,
                &[0, 1],
                AveragineComposition::PEPTIDE,
            )
            .unwrap();
        assert_eq!(bits(&joint), probe_peaks("fragment", label), "{label}");
    }
    let nominal = source(Some(11), NOMINAL);
    let c1 = nominal.run(&formula("C1")).unwrap();
    let c2 = nominal.run(&formula("C2")).unwrap();
    let joint = generator
        .calc_fragment_isotope_dist(&c1, &c2, &[0, 1, 2], formula("C1").mono_mass())
        .unwrap();
    assert_eq!(bits(&joint), probe_peaks("fragment", "calc_C1_C2_iso012"));
}

#[test]
fn coarse_class_test_literals_hold_in_source_precision() {
    // CoarseIsotopeDistribution_test.cpp:122-180 (run).
    let glucose = source(Some(3), APPROXIMATE)
        .run(&formula("C6H12O6"))
        .unwrap();
    assert_eq!(glucose.len(), 3);
    similar(glucose.peaks()[0].mass, 180.063, "glucose m0");
    similar(glucose.peaks()[0].probability, 0.923456, "glucose p0");
    similar(glucose.peaks()[2].mass, 182.0701, "glucose m2");
    similar(glucose.peaks()[2].probability, 0.013232, "glucose p2");
    let charged = source(Some(3), APPROXIMATE)
        .run(&formula("C6H12O6+2"))
        .unwrap();
    assert_eq!(charged.len(), 3);
    similar(charged.peaks()[0].mass, 182.077943, "charged m0");
    similar(charged.peaks()[0].probability, 0.923246, "charged p0");
    similar(charged.peaks()[2].mass, 184.0846529, "charged m2");
    similar(charged.peaks()[2].probability, 0.0132435, "charged p2");
    let explicit = source(Some(3), APPROXIMATE)
        .run(&formula("C6H14O6"))
        .unwrap();
    assert_eq!(explicit.len(), 3);
    assert!((explicit.peaks()[0].mass - 182.077943).abs() < 0.005);
    assert!(
        source(Some(3), APPROXIMATE)
            .run(&formula("C6H12O6-2"))
            .is_err()
    );

    // :182-251 (convolvePow_).
    let heavy = source(Some(11), APPROXIMATE)
        .run(&formula("C222N190O110"))
        .unwrap();
    let expected = [
        0.0349429, 0.109888, 0.180185, 0.204395, 0.179765, 0.130358, 0.0809864, 0.0442441,
        0.0216593, 0.00963707, 0.0039406,
    ];
    assert_eq!(heavy.len(), expected.len());
    for (index, (peak, value)) in heavy.peaks().iter().zip(expected).enumerate() {
        assert_eq!(peak.mass.round(), 7084.0 + index as f64);
        similar(peak.probability, value, "C222N190O110");
    }
    for (text, max, first, values) in [
        (
            "Br2",
            5,
            158.0,
            vec![0.2569476, 0.0, 0.49990478, 0.0, 0.24314761],
        ),
        (
            "CBr2",
            7,
            170.0,
            vec![
                0.254198270573,
                0.002749339427,
                0.494555798854,
                0.005348981146,
                0.240545930573,
                0.002601679427,
            ],
        ),
    ] {
        let pattern = source(Some(max), APPROXIMATE).run(&formula(text)).unwrap();
        assert_eq!(pattern.len(), values.len(), "{text}");
        for (index, (peak, value)) in pattern.peaks().iter().zip(values).enumerate() {
            assert_eq!(peak.mass.round(), first + index as f64, "{text}");
            similar(peak.probability, value, text);
        }
    }

    // :253-295 (estimateFromWeightAndComp, estimateFromPeptideWeight).
    let generator = source(Some(3), APPROXIMATE);
    let composed = generator
        .estimate_from_weight_and_comp(1000.0, AveragineComposition::PEPTIDE)
        .unwrap();
    assert_eq!(
        composed,
        generator.estimate_from_peptide_weight(1000.0).unwrap()
    );
    for (mass, probability, first_mass, rounded) in [
        (100.0, 0.949735, 100.170, 100.0),
        (1000.0, 0.586906, 999.714, 1000.0),
        (10000.0, 0.046495, 9994.041, 9994.0),
    ] {
        let pattern = generator.estimate_from_peptide_weight(mass).unwrap();
        similar(pattern.peaks()[0].probability, probability, "peptide p0");
        similar(pattern.peaks()[0].mass, first_mass, "peptide m0");
        let nominal = source(Some(3), NOMINAL)
            .estimate_from_peptide_weight(mass)
            .unwrap();
        similar(nominal.peaks()[0].mass, rounded, "rounded peptide m0");
    }

    // :297-353 (Poisson approximations against the source-precision truth).
    for mass in [20.0, 300.0, 1000.0, 2500.0] {
        let approximation = Generator::approximate_from_peptide_weight(mass, 20, 1).unwrap();
        let intensities = Generator::approximate_intensities(mass, 20).unwrap();
        let truth = source(Some(approximation.len()), APPROXIMATE)
            .estimate_from_peptide_weight(mass)
            .unwrap();
        for q in [
            approximation
                .peaks()
                .iter()
                .map(|p| p.probability)
                .collect::<Vec<_>>(),
            intensities,
        ] {
            let mut kl = 0.0;
            let mut sum = 0.0;
            for (p, q) in truth.peaks().iter().map(|p| p.probability).zip(q) {
                sum += q;
                if p != 0.0 {
                    kl += p * (p / q).ln();
                }
            }
            similar(sum, 1.0, "Poisson sum");
            assert!(kl > 0.0 && kl < 0.05, "mass {mass}: KL {kl}");
        }
    }

    // :400-439 (estimateFromPeptideWeightAndS, renormalized tail).
    for (sulfur, tail) in [
        (0, 0.00290370998965918),
        (1, 0.0439547771832361),
        (2, 0.0804989104418586),
        (3, 0.117023432503842),
    ] {
        let mut pattern = generator
            .estimate_from_peptide_weight_and_sulfur(100.0, sulfur)
            .unwrap();
        pattern.renormalize_with(SINGLE).unwrap();
        similar(
            pattern.peaks().last().unwrap().probability,
            tail,
            "sulfur tail",
        );
    }

    // :441-472 (RNA and DNA averagine).
    for (mass, rna, dna) in [
        (100.0, 0.958166, 0.958166),
        (1000.0, 0.668538, 0.657083),
        (10000.0, 0.080505, 0.075138),
    ] {
        let r = generator.estimate_from_rna_weight(mass).unwrap();
        let d = generator.estimate_from_dna_weight(mass).unwrap();
        similar(r.peaks()[0].probability, rna, "RNA p0");
        similar(d.peaks()[0].probability, dna, "DNA p0");
    }

    // :474-652 (fragment estimates with isolated M0 and M+1).
    let unbounded = source(None, APPROXIMATE);
    for (composition, cases) in [
        (
            AveragineComposition::PEPTIDE,
            [
                (200.0, 100.0, 0.954654801320083),
                (2000.0, 100.0, 0.975984866212216),
                (20000.0, 100.0, 0.995783521351781),
                (2000.0, 1000.0, 0.741290977639283),
                (20000.0, 1000.0, 0.95467154987681),
                (20000.0, 10000.0, 0.542260764523188),
            ],
        ),
        (
            AveragineComposition::DNA,
            [
                (200.0, 100.0, 0.963845242419331),
                (2000.0, 100.0, 0.978300783455351),
                (20000.0, 100.0, 0.995652529512413),
                (2000.0, 1000.0, 0.776727852910751),
                (20000.0, 1000.0, 0.95504592203456),
                (20000.0, 10000.0, 0.555730613643729),
            ],
        ),
        (
            AveragineComposition::RNA,
            [
                (200.0, 100.0, 0.963845242419331),
                (2000.0, 100.0, 0.977854088814216),
                (20000.0, 100.0, 0.995465661923629),
                (2000.0, 1000.0, 0.784037437107401),
                (20000.0, 1000.0, 0.955768644474843),
                (20000.0, 10000.0, 0.558201381343203),
            ],
        ),
    ] {
        for (precursor, fragment, expected) in cases {
            let mut joint = unbounded
                .estimate_fragment_from_weights(precursor, fragment, &[0, 1], composition)
                .unwrap();
            joint.renormalize_with(SINGLE).unwrap();
            similar(joint.peaks()[0].probability, expected, "fragment M0");
        }
        let whole = unbounded
            .estimate_fragment_from_weights(200.0, 200.0, &[0, 1], composition)
            .unwrap();
        let precursor = if composition == AveragineComposition::PEPTIDE {
            unbounded.estimate_from_peptide_weight(200.0).unwrap()
        } else {
            source(Some(2), APPROXIMATE)
                .estimate_from_weight_and_comp(200.0, composition)
                .unwrap()
        };
        for (a, b) in whole.peaks().iter().zip(precursor.peaks()) {
            assert_eq!(a.mass, b.mass);
        }
    }
    for (precursor, fragment, mass) in [(200.0, 100.0, 100.170), (2000.0, 100.0, 100.170)] {
        let joint = unbounded
            .estimate_fragment_from_weights(
                precursor,
                fragment,
                &[0, 1],
                AveragineComposition::PEPTIDE,
            )
            .unwrap();
        similar(joint.peaks()[0].mass, mass, "fragment m0");
    }
    let rounded = source(None, NOMINAL)
        .estimate_fragment_from_weights(200.0, 100.0, &[0, 1], AveragineComposition::PEPTIDE)
        .unwrap();
    assert_eq!(rounded.peaks()[0].mass, 100.0);

    // :654-717 (calcFragmentIsotopeDist).
    let nominal = source(Some(11), NOMINAL);
    let iso1 = nominal.run(&formula("C1")).unwrap();
    let iso2 = nominal.run(&formula("C2")).unwrap();
    let mono = formula("C1").mono_mass();
    let mut iso3 = unbounded
        .calc_fragment_isotope_dist(&iso1, &iso2, &[0, 1, 2], mono)
        .unwrap();
    iso3.renormalize_with(SINGLE).unwrap();
    let calc_mass = source(Some(11), APPROXIMATE).run(&formula("C1")).unwrap();
    for (a, b) in calc_mass.peaks().iter().zip(iso3.peaks()) {
        assert_eq!(a.mass, b.mass);
        similar(a.probability, b.probability, "all precursors isolated");
    }
    let mut iso4 = unbounded
        .calc_fragment_isotope_dist(&iso1, &iso2, &[0, 1], mono)
        .unwrap();
    iso4.renormalize_with(SINGLE).unwrap();
    assert_eq!(calc_mass.peaks()[0].mass, iso4.peaks()[0].mass);
    assert_eq!(calc_mass.peaks()[1].mass, iso4.peaks()[1].mass);
    similar(iso1.peaks()[0].probability, 0.989300, "C1 p0");
    similar(iso1.peaks()[1].probability, 0.010700, "C1 p1");
    similar(iso4.peaks()[0].probability, 0.989524, "iso4 p0");
    similar(iso4.peaks()[1].probability, 0.010479, "iso4 p1");
    let iso5 = source(None, NOMINAL)
        .calc_fragment_isotope_dist(&iso1, &iso2, &[0, 1], mono)
        .unwrap();
    for ((a, b), (mass, rounded)) in iso3
        .peaks()
        .iter()
        .zip(iso5.peaks())
        .zip([(12.0, 12.0), (13.0033548378, 13.0)])
    {
        similar(a.mass, mass, "calculated mass");
        assert_eq!(b.mass, rounded);
    }
}

#[test]
fn distribution_class_test_literals_hold_in_source_precision() {
    // IsotopeDistribution_test.cpp at bc9cc12.
    assert_eq!(IsotopeDistribution::default().len(), 1);
    let eleven = source(Some(11), APPROXIMATE);
    let eleven_round = source(Some(11), NOMINAL);
    let c4 = eleven.run(&formula("C4")).unwrap();
    assert_eq!(c4, eleven.run(&formula("C4")).unwrap());
    assert_ne!(eleven_round.run(&formula("C4")).unwrap(), c4);
    assert_ne!(IsotopeDistribution::default(), c4);
    assert_eq!(c4.len(), 5);
    let mut cleared = c4.clone();
    cleared.clear();
    assert_eq!(cleared.len(), 0);

    let mut h2 = eleven.run(&formula("H2")).unwrap();
    similar(h2.max_mass().unwrap(), 6.02907, "H2 max");
    similar(h2.min_mass().unwrap(), 2.01565, "H2 min");
    let h2_round = eleven_round.run(&formula("H2")).unwrap();
    assert_eq!(h2_round.max_mass(), Some(6.0));
    assert_eq!(h2_round.min_mass(), Some(2.0));
    similar(c4.min_mass().unwrap(), 48.0, "C4 min");
    assert_eq!(
        eleven_round.run(&formula("C4")).unwrap().min_mass(),
        Some(48.0)
    );
    let mut h2_min = h2.clone();
    for mass in [11.2, 10.2] {
        h2.insert(IsotopePeak {
            mass,
            probability: 2.0,
        })
        .unwrap();
    }
    similar(h2.max_mass().unwrap(), 11.2, "inserted max");
    for mass in [1.2, 10.2] {
        h2_min
            .insert(IsotopePeak {
                mass,
                probability: 2.0,
            })
            .unwrap();
    }
    similar(h2_min.min_mass().unwrap(), 1.2, "inserted min");

    let c1 = eleven_round.run(&formula("C1")).unwrap();
    assert_eq!(c1.most_abundant().unwrap().mass, 12.0);
    let mut c100 = eleven_round.run(&formula("C100")).unwrap();
    assert_eq!(c100.most_abundant().unwrap().mass, 1201.0);
    c100.clear();
    assert_eq!(c100.most_abundant(), None); // source: Peak1D(0, 1)

    let ten = source(Some(10), APPROXIMATE);
    let mut c160 = ten.run(&formula("C160")).unwrap();
    assert_ne!(c160.len(), 3);
    c160.trim_right(0.2).unwrap();
    assert_eq!(c160.len(), 3);
    c160.trim_left_source(0.2).unwrap();
    assert_eq!(c160.len(), 2);
    c160.renormalize_with(SINGLE).unwrap();
    let sum: f64 = c160.peaks().iter().map(|p| p.probability).sum();
    similar(sum, 1.0, "renormalized sum");
}

#[test]
fn source_precision_is_explicit_and_the_default_is_unchanged() {
    // CoarseIsotopeDistribution_test.cpp:44-91 (constructors, mass labels).
    let unbounded = Generator::default();
    assert_eq!(unbounded.max_peaks(), None, "source getMaxIsotope() == 0");
    assert_eq!(unbounded.mass_mode(), APPROXIMATE);
    let bounded = Generator::new(Some(117), APPROXIMATE).unwrap();
    assert_eq!(
        (bounded.max_peaks(), bounded.mass_mode()),
        (Some(117), APPROXIMATE)
    );
    let rounded = Generator::new(Some(117), NOMINAL).unwrap();
    assert_eq!(
        (rounded.max_peaks(), rounded.mass_mode()),
        (Some(117), NOMINAL)
    );
    assert!(Generator::new(Some(0), APPROXIMATE).is_err());

    let default = Generator::new(Some(11), APPROXIMATE).unwrap();
    assert_eq!(default.precision(), ProbabilityPrecision::Double);
    assert_eq!(
        Generator::default().precision(),
        ProbabilityPrecision::Double
    );
    assert_eq!(
        ProbabilityPrecision::default(),
        ProbabilityPrecision::Double
    );
    let explicit_double = default.clone().with_precision(ProbabilityPrecision::Double);
    let heavy = formula("C222N190O110");
    let native = default.run(&heavy).unwrap();
    assert_eq!(explicit_double.run(&heavy).unwrap(), native);
    assert!(
        native
            .peaks()
            .iter()
            .any(|p| f64::from(p.probability as f32) != p.probability),
        "the default does not narrow to binary32"
    );
    let mut toggled = default.clone();
    toggled.set_precision(SINGLE);
    assert_eq!(toggled.precision(), SINGLE);
    let single = toggled.run(&heavy).unwrap();
    assert_ne!(single, native);
    for (a, b) in single.peaks().iter().zip(native.peaks()) {
        assert_eq!(a.mass.to_bits(), b.mass.to_bits(), "same mass labels here");
        similar(a.probability, b.probability, "binary32 against f64");
    }
}

#[test]
fn source_precision_rejects_values_outside_binary32_and_keeps_limits() {
    let mut huge = source(None, APPROXIMATE);
    huge.set_isotope_override("C", distribution(&[(12.0, 1e39), (13.0, 1.0)]))
        .unwrap();
    assert!(huge.run(&formula("C1")).is_err());
    let mut double_huge = Generator::new(None, APPROXIMATE).unwrap();
    double_huge
        .set_isotope_override("C", distribution(&[(12.0, 1e39), (13.0, 1.0)]))
        .unwrap();
    assert!(double_huge.run(&formula("C1")).is_ok());

    let mut overflow = source(None, APPROXIMATE);
    overflow
        .set_isotope_override("C", distribution(&[(12.0, 3e38), (13.0, 3e38)]))
        .unwrap();
    assert!(overflow.run(&formula("C2")).is_err());

    let mut wide_weight = distribution(&[(0.0, 1e39)]);
    let old = wide_weight.clone();
    assert!(wide_weight.renormalize_with(SINGLE).is_err());
    assert_eq!(wide_weight, old);
    let mut zero = distribution(&[(0.0, 0.0)]);
    assert!(zero.renormalize_with(SINGLE).is_err());
    let mut tiny = distribution(&[(0.0, 1e-50)]);
    assert!(
        tiny.renormalize_with(SINGLE).is_err(),
        "narrows to zero weight"
    );
    let mut empty = IsotopeDistribution::empty();
    empty.renormalize_with(SINGLE).unwrap();
    assert!(empty.is_empty());

    let identity = IsotopeDistribution::default();
    let wide = distribution(&[(0.0, 0.5), (9999.0, 0.5)]);
    assert!(source(None, APPROXIMATE).convolve(&wide, &wide).is_err());
    assert_eq!(
        source(Some(3), APPROXIMATE)
            .convolve(&wide, &wide)
            .unwrap()
            .peaks()
            .iter()
            .map(|p| p.probability)
            .collect::<Vec<_>>(),
        [0.25, 0.0, 0.0]
    );
    let huge_gap = distribution(&[(0.0, 0.5), (1e12, 0.5)]);
    assert!(
        source(None, APPROXIMATE)
            .convolve(&huge_gap, &identity)
            .is_err()
    );
    assert!(
        source(Some(3), APPROXIMATE)
            .convolve(&huge_gap, &identity)
            .is_ok()
    );
    let carbon = distribution(&[(12.0, 0.9893), (13.003355, 0.0107)]);
    assert_eq!(
        source(None, APPROXIMATE)
            .convolve_power(&carbon, 0)
            .unwrap(),
        identity
    );
    assert!(
        source(None, APPROXIMATE)
            .convolve_power(&IsotopeDistribution::empty(), 3)
            .unwrap()
            .is_empty()
    );
    assert!(source(None, APPROXIMATE).run(&formula("C-1H2")).is_err());
    assert!(
        source(None, APPROXIMATE)
            .calc_fragment_isotope_dist(&carbon, &carbon, &[], 12.0)
            .is_err()
    );
}
