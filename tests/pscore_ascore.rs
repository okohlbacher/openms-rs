// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

//! `ANALYSIS/ID/PScore.h` and `ANALYSIS/ID/AScore.h`.
//!
//! Expected values are the literals of the upstream class tests at the core
//! pin `bc9cc12`: `src/tests/class_tests/openms/source/PScore_test.cpp` and
//! `AScore_test.cpp`, with their DTA inputs copied unchanged from
//! `src/tests/class_tests/openms/data/` (hashes in
//! `tests/data/pscore_ascore_provenance.json`). `real_similar` is the class
//! tests' `TEST_REAL_SIMILAR` at its default tolerance (absolute or relative
//! 1e-5). Tests marked "native" check port-only behaviour: refusals where the
//! source is undefined, and invariants derived from the definitions.

use openms::MSSpectrum;
use openms::Peak1D;
use openms::analysis::id::ascore::{AScore, MAX_ASCORE_PERMUTATIONS, ProbablePhosphoSites};
use openms::analysis::id::pscore::{DEFAULT_MZ_WINDOW, PScore};
use openms::chemistry::AASequence;
use openms::chemistry::theoretical::TheoreticalSpectrumGenerator;
use openms::format::dta;
use openms::identification::PeptideHit;
use openms::kernel::MSExperiment;
use openms::param::{Param, ParamValue};
use std::cmp::Ordering;

fn real_similar(a: f64, b: f64) {
    let absdiff = (a - b).abs();
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{a} !~ {b} (absdiff {absdiff}, ratio {ratio})"
    );
}

fn load(bytes: &[u8]) -> MSSpectrum {
    dta::read(bytes).unwrap()
}

fn pscore_dta() -> MSSpectrum {
    load(include_bytes!("data/PScore_test.dta"))
}

fn ascore_input(n: usize) -> MSSpectrum {
    match n {
        1 => load(include_bytes!("data/Ascore_test_input1.dta")),
        2 => load(include_bytes!("data/Ascore_test_input2.dta")),
        3 => load(include_bytes!("data/Ascore_test_input3.dta")),
        _ => load(include_bytes!("data/Ascore_test_input4.dta")),
    }
}

fn mz_and_intensities(spec: &MSSpectrum) -> (Vec<f64>, Vec<f64>) {
    (
        spec.peaks.iter().map(|p| p.mz).collect(),
        spec.peaks.iter().map(|p| f64::from(p.intensity)).collect(),
    )
}

fn with(values: &[(&str, ParamValue)]) -> AScore {
    let mut param = Param::new();
    for (key, value) in values {
        param.set_value(key, value.clone(), "", &[]).unwrap();
    }
    AScore::from_param(&param).unwrap()
}

fn seq(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}

// ---------------------------------------------------------------- PScore

#[test]
fn rank_in_mz_window_counts_higher_neighbours() {
    // PScore_test.cpp: calculateIntensityRankInMZWindow.
    let values: Vec<f64> = (0..100).map(f64::from).collect();
    let ranks = PScore::calculate_intensity_rank_in_mz_window(&values, &values, 9.9).unwrap();
    assert_eq!(ranks.len(), values.len());
    for &rank in &ranks[..ranks.len() - 4] {
        assert_eq!(rank, 4);
    }
    let ranks = PScore::calculate_intensity_rank_in_mz_window(&values, &values, 10.1).unwrap();
    assert_eq!(ranks.len(), values.len());
    for &rank in &ranks[..ranks.len() - 5] {
        assert_eq!(rank, 5);
    }
}

#[test]
fn peak_level_spectra_keep_top_peaks_per_window() {
    // PScore_test.cpp: calculatePeakLevelSpectra.
    let spec = pscore_dta();
    let (mz, intensities) = mz_and_intensities(&spec);
    let ranks = PScore::calculate_intensity_rank_in_mz_window(&mz, &intensities, 100.0).unwrap();
    assert_eq!(ranks.len(), spec.len());
    let pls = PScore::calculate_peak_level_spectra(&spec, &ranks, 0, 1).unwrap();
    assert_eq!(pls.len(), 2);
    let level0 = [169.65, 231.51, 362.22, 508.47, 579.61, 629.66, 712.18];
    for (i, &expected) in level0.iter().enumerate() {
        real_similar(pls[&0].peaks[i].mz, expected);
    }
    let level1 = [
        149.93, 169.65, 231.51, 263.88, 318.38, 362.22, 389.84, 489.86, 508.47, 562.72, 579.61,
        629.66, 712.18,
    ];
    for (i, &expected) in level1.iter().enumerate() {
        real_similar(pls[&1].peaks[i].mz, expected);
    }
}

#[test]
fn pscore_of_self_matching_spectra() {
    // PScore_test.cpp: computePScore(..., const PeakSpectrum&, ...).
    let spec = pscore_dta();
    let (mz, intensities) = mz_and_intensities(&spec);
    let theo_spec = MSSpectrum {
        peaks: spec.peaks.clone(),
        ..MSSpectrum::default()
    };
    let ranks = PScore::calculate_intensity_rank_in_mz_window(&mz, &intensities, 100.0).unwrap();
    let pls = PScore::calculate_peak_level_spectra(&spec, &ranks, 0, 0).unwrap();
    let top_1 = PScore::compute_pscore(0.1, true, &pls, &theo_spec, DEFAULT_MZ_WINDOW).unwrap();
    let pls = PScore::calculate_peak_level_spectra(&spec, &ranks, 0, 1).unwrap();
    let top_2 = PScore::compute_pscore(0.1, true, &pls, &theo_spec, DEFAULT_MZ_WINDOW).unwrap();
    real_similar(top_1, 83.867454);
    real_similar(top_2, 154.682242);

    let generator = TheoreticalSpectrumGenerator {
        add_first_prefix_ion: true,
        ..TheoreticalSpectrumGenerator::default()
    };
    let spec = generator.generate(&seq("IFSQVGK"), 1, 1, None).unwrap();
    assert_eq!(spec.len(), 12);
    let (mz, intensities) = mz_and_intensities(&spec);
    let ranks = PScore::calculate_intensity_rank_in_mz_window(&mz, &intensities, 100.0).unwrap();
    let pls = PScore::calculate_peak_level_spectra(&spec, &ranks, 0, 0).unwrap();
    let all_match = PScore::compute_pscore(0.1, true, &pls, &spec, DEFAULT_MZ_WINDOW).unwrap();
    real_similar(all_match, 240.0);
}

#[test]
fn rank_map_is_the_per_spectrum_rank() {
    // native: the source test leaves calculateRankMap untested ("convenience").
    let spec = pscore_dta();
    let (mz, intensities) = mz_and_intensities(&spec);
    let experiment = MSExperiment {
        spectra: vec![spec.clone(), MSSpectrum::default(), spec],
        ..MSExperiment::default()
    };
    let map = PScore::calculate_rank_map(&experiment, 100.0).unwrap();
    let ranks = PScore::calculate_intensity_rank_in_mz_window(&mz, &intensities, 100.0).unwrap();
    assert_eq!(map, vec![ranks.clone(), Vec::new(), ranks]);
}

#[test]
fn best_of_overload_keeps_the_source_level_probability() {
    // native: the vector overload uses p = level / mz_window. Level 1 of the
    // DTA spectrum matched against itself gives N = n = its peak count and
    // p = 0.01, so the score is -10 log10(0.01^N) = 20 N; level 0 has p = 0
    // and any match there scores +inf.
    let spec = pscore_dta();
    let (mz, intensities) = mz_and_intensities(&spec);
    let ranks = PScore::calculate_intensity_rank_in_mz_window(&mz, &intensities, 100.0).unwrap();
    let pls = PScore::calculate_peak_level_spectra(&spec, &ranks, 1, 1).unwrap();
    let level1 = pls[&1].clone();
    let score = PScore::compute_pscore_best_of(
        0.1,
        true,
        &pls,
        std::slice::from_ref(&level1),
        DEFAULT_MZ_WINDOW,
    )
    .unwrap();
    real_similar(score, 20.0 * level1.len() as f64);
    let pls = PScore::calculate_peak_level_spectra(&spec, &ranks, 0, 0).unwrap();
    let score = PScore::compute_pscore_best_of(
        0.1,
        true,
        &pls,
        std::slice::from_ref(&spec),
        DEFAULT_MZ_WINDOW,
    )
    .unwrap();
    assert_eq!(score, f64::INFINITY);
    // An empty level spectrum is the source's findNearest precondition.
    let mut empty = pls.clone();
    empty.insert(3, MSSpectrum::default());
    assert!(PScore::compute_pscore_best_of(0.1, true, &empty, &[spec], 100.0).is_err());
}

#[test]
fn correction_terms_follow_the_source_tables() {
    // native: the source test leaves these untested; values transcribed from
    // PScore.cpp, where they are free functions the declared statics never reach.
    real_similar(PScore::mass_correction_term(1600.0), 24.0);
    assert_eq!(PScore::cleavage_correction_term(0, false), 53.2);
    assert_eq!(PScore::cleavage_correction_term(1, true), 42.1);
    assert_eq!(PScore::cleavage_correction_term(1, false), 31.1);
    assert_eq!(PScore::cleavage_correction_term(2, true), 17.0);
    assert_eq!(PScore::cleavage_correction_term(3, true), 0.0);
    let table = [42.0, 28.0, 22.0, 16.0, 9.0, 5.0, 2.0, 0.0, 0.0];
    for (modifications, &expected) in table.iter().enumerate() {
        assert_eq!(
            PScore::modification_correction_term(modifications),
            expected
        );
    }
}

#[test]
fn pscore_refuses_where_the_source_is_undefined() {
    // native.
    assert!(PScore::calculate_intensity_rank_in_mz_window(&[1.0, 2.0], &[1.0], 10.0).is_err());
    assert!(PScore::calculate_intensity_rank_in_mz_window(&[1.0], &[f64::NAN], 10.0).is_err());
    assert!(
        PScore::calculate_intensity_rank_in_mz_window(&[], &[], 10.0)
            .unwrap()
            .is_empty()
    );
    let spec = pscore_dta();
    let too_many = vec![0; spec.len() + 1];
    assert!(PScore::calculate_peak_level_spectra(&spec, &too_many, 0, 1).is_err());
    assert!(
        PScore::calculate_peak_level_spectra(&MSSpectrum::default(), &[], 0, 9)
            .unwrap()
            .is_empty()
    );
    assert!(
        PScore::calculate_peak_level_spectra(&spec, &[0; 3], 5, 2)
            .unwrap()
            .is_empty()
    );
    let mut unsorted = spec.clone();
    unsorted.peaks.reverse();
    let pls = PScore::calculate_peak_level_spectra(&spec, &vec![0; spec.len()], 0, 0).unwrap();
    assert!(PScore::compute_pscore(0.1, true, &pls, &unsorted, 100.0).is_err());
}

// ---------------------------------------------------------------- AScore

#[test]
fn cumulative_score() {
    // AScore_test.cpp: computeCumulativeScore_.
    assert!(AScore::compute_cumulative_score(1, 5, 0.1).is_err());
    real_similar(AScore::compute_cumulative_score(1, 1, 0.1).unwrap(), 0.1);
    real_similar(AScore::compute_cumulative_score(3, 1, 0.1).unwrap(), 0.271);
}

fn rows(values: &[f64]) -> Vec<Vec<f64>> {
    values.iter().map(|&v| vec![v; 10]).collect()
}

fn check(
    site: &ProbablePhosphoSites,
    first: usize,
    second: usize,
    seq_1: usize,
    seq_2: usize,
    depth: usize,
) {
    assert_eq!(
        (
            site.first,
            site.second,
            site.seq_1,
            site.seq_2,
            site.peak_depth
        ),
        (first, second, seq_1, seq_2, depth)
    );
}

#[test]
fn highest_scoring_permutations() {
    // AScore_test.cpp: determineHighestScoringPermutations_.
    let permutations = vec![vec![1, 3, 5], vec![3, 5, 6], vec![1, 3, 6], vec![1, 5, 6]];

    let scores_1 = rows(&[0.1, 0.2, 0.3, 0.4]);
    let ranking = AScore::rank_weighted_permutation_peptide_scores(&scores_1).unwrap();
    real_similar(ranking.last().unwrap().0, 0.4);
    let sites =
        AScore::determine_highest_scoring_permutations(&scores_1, &permutations, &ranking).unwrap();
    assert_eq!(sites.len(), 3);
    check(&sites[0], 1, 3, 3, 1, 1);
    check(&sites[1], 5, 3, 3, 2, 1);
    check(&sites[2], 6, 3, 3, 0, 1);

    let scores_3 = rows(&[0.1, 0.4, 0.3, 0.2]);
    let ranking = AScore::rank_weighted_permutation_peptide_scores(&scores_3).unwrap();
    real_similar(ranking.last().unwrap().0, 0.4);
    let sites =
        AScore::determine_highest_scoring_permutations(&scores_3, &permutations, &ranking).unwrap();
    assert_eq!(sites.len(), 3);
    check(&sites[0], 3, 1, 1, 3, 1);
    check(&sites[1], 5, 1, 1, 2, 1);
    check(&sites[2], 6, 1, 1, 0, 1);

    let scores_2 = rows(&[0.2, 0.3, 0.4, 0.1]);
    let ranking = AScore::rank_weighted_permutation_peptide_scores(&scores_2).unwrap();
    real_similar(ranking.last().unwrap().0, 0.4);
    let sites =
        AScore::determine_highest_scoring_permutations(&scores_2, &permutations, &ranking).unwrap();
    assert_eq!(sites.len(), 3);
    check(&sites[0], 1, 5, 2, 1, 1);
    check(&sites[1], 3, 5, 2, 3, 1);
    check(&sites[2], 6, 5, 2, 0, 1);

    let best = vec![
        55.0, 60.0, 75.0, 100.0, 90.0, 120.0, 125.0, 120.0, 100.0, 90.0,
    ];
    let other = vec![40.0, 50.0, 53.0, 60.0, 50.0, 53.0, 59.0, 53.0, 50.0, 40.0];
    let mut scores = vec![best, other.clone()];
    let permutations = vec![vec![3], vec![6]];
    let ranking = AScore::rank_weighted_permutation_peptide_scores(&scores).unwrap();
    real_similar(ranking.last().unwrap().0, 94.10714285714287);
    let sites =
        AScore::determine_highest_scoring_permutations(&scores, &permutations, &ranking).unwrap();
    assert_eq!(sites.len(), 1);
    check(&sites[0], 3, 6, 0, 1, 6);

    let permutations = vec![
        vec![3, 5],
        vec![5, 6],
        vec![3, 7],
        vec![3, 6],
        vec![5, 7],
        vec![6, 7],
    ];
    for _ in 0..4 {
        scores.push(other.clone());
    }
    let ranking = AScore::rank_weighted_permutation_peptide_scores(&scores).unwrap();
    let sites =
        AScore::determine_highest_scoring_permutations(&scores, &permutations, &ranking).unwrap();
    assert_eq!(sites.len(), 2);
    check(&sites[0], 3, 7, 0, 4, 6);
    check(&sites[1], 5, 6, 0, 3, 6);
}

fn candidates(seq_1: usize, seq_2: usize, first: usize, second: usize) -> ProbablePhosphoSites {
    ProbablePhosphoSites {
        first,
        second,
        seq_1,
        seq_2,
        ..ProbablePhosphoSites::default()
    }
}

fn first_last(spectrum: &MSSpectrum) -> (f64, f64) {
    (spectrum.peaks[0].mz, spectrum.peaks[spectrum.len() - 1].mz)
}

fn qsp_permutations() -> Vec<Vec<usize>> {
    vec![vec![1], vec![2], vec![4], vec![7], vec![10]]
}

#[test]
fn site_determining_ions() {
    // AScore_test.cpp: computeSiteDeterminingIons_, default tolerance 0.05 Da.
    let ascore = AScore::default();
    let th = ascore
        .create_theoretical_spectra(&qsp_permutations(), &seq("QSSVTQVTEQSPK"))
        .unwrap();
    let ions = ascore
        .compute_site_determining_ions(&th, &candidates(3, 4, 10, 7))
        .unwrap();
    assert_eq!((ions[0].len(), ions[1].len()), (6, 6));

    // (sequence, positions, sizes, first/last of seq_1-side, first/last of seq_2-side)
    type Case = (&'static str, [usize; 2], [usize; 2], [f64; 2], [f64; 2]);
    let cases: [Case; 4] = [
        (
            "VTEQSP",
            [1, 4],
            [6, 6],
            [203.102, 538.19],
            [201.123, 540.17],
        ),
        (
            "TYQYS",
            [0, 4],
            [7, 7],
            [106.05, 636.206],
            [186.016, 640.201],
        ),
        (
            "TSTYQYSYPP",
            [2, 6],
            [8, 8],
            [370.101, 917.403],
            [290.135, 997.37],
        ),
        (
            "ATPGNLGSSVLMYK",
            [12, 8],
            [8, 4],
            [390.142, 1128.57],
            [310.176, 1208.54],
        ),
    ];
    for (sequence, [a, b], [size_a, size_b], ends_a, ends_b) in cases {
        let permutations = vec![vec![a], vec![b]];
        let th = ascore
            .create_theoretical_spectra(&permutations, &seq(sequence))
            .unwrap();
        let ions = ascore
            .compute_site_determining_ions(&th, &candidates(0, 1, a, b))
            .unwrap();
        assert_eq!(
            (ions[0].len(), ions[1].len()),
            (size_a, size_b),
            "{sequence}"
        );
        let (lo, hi) = first_last(&ions[0]);
        real_similar(lo, ends_a[0]);
        real_similar(hi, ends_a[1]);
        let (lo, hi) = first_last(&ions[1]);
        real_similar(lo, ends_b[0]);
        real_similar(hi, ends_b[1]);

        let ions = ascore
            .compute_site_determining_ions(&th, &candidates(1, 0, b, a))
            .unwrap();
        assert_eq!(
            (ions[1].len(), ions[0].len()),
            (size_a, size_b),
            "{sequence}"
        );
        let (lo, hi) = first_last(&ions[1]);
        real_similar(lo, ends_a[0]);
        real_similar(hi, ends_a[1]);
        let (lo, hi) = first_last(&ions[0]);
        real_similar(lo, ends_b[0]);
        real_similar(hi, ends_b[1]);
    }
}

#[test]
fn sites_and_permutations() {
    // AScore_test.cpp: getSites_ (on toUniModString of an unmodified peptide,
    // which is its residue string) and computePermutations_.
    assert_eq!(
        AScore::get_sites(seq("VTQSPSSP").as_str()),
        vec![1, 3, 5, 6]
    );
    let tupel = [1, 2, 3, 4];
    assert_eq!(
        AScore::compute_permutations(&tupel, 1).unwrap(),
        vec![vec![1], vec![2], vec![3], vec![4]]
    );
    assert_eq!(
        AScore::compute_permutations(&tupel, 2).unwrap(),
        vec![
            vec![1, 2],
            vec![1, 3],
            vec![1, 4],
            vec![2, 3],
            vec![2, 4],
            vec![3, 4]
        ]
    );
    assert_eq!(
        AScore::compute_permutations(&tupel, 3).unwrap(),
        vec![vec![1, 2, 3], vec![1, 2, 4], vec![1, 3, 4], vec![2, 3, 4]]
    );
    assert_eq!(
        AScore::compute_permutations(&tupel, 4).unwrap(),
        vec![vec![1, 2, 3, 4]]
    );
    assert!(AScore::compute_permutations(&[], 0).unwrap().is_empty());
}

#[test]
fn remove_phosphosites() {
    // AScore_test.cpp: removePhosphositesFromSequence_.
    let without = AScore::remove_phosphosites_from_sequence("QSSVTQVTEQS(Phospho)PK").unwrap();
    assert_eq!(without.to_string(), "QSSVTQVTEQSPK");
}

#[test]
fn windows_and_matched_ions() {
    // AScore_test.cpp: peakPickingPerWindowsInSpectrum_ and numberOfMatchedIons_.
    let spectrum = ascore_input(3);
    let windows = AScore::peak_picking_per_windows_in_spectrum(&spectrum).unwrap();
    assert_eq!(windows.len(), 8);
    assert_eq!(windows[0].len(), 1);
    assert_eq!(windows[1].len(), 1);
    assert_eq!(windows[4].len(), 0);
    assert_eq!(windows[7].len(), 1);

    let ascore = with(&[("fragment_mass_tolerance", ParamValue::Float(0.5))]);
    let th = ascore
        .create_theoretical_spectra(&qsp_permutations(), &seq("QSSVTQVTEQSPK"))
        .unwrap();
    let last = &th[th.len() - 1]; // QSSVTQVTEQS(Phospho)PK
    assert_eq!(
        ascore.number_of_matched_ions(last, &windows[0], 1).unwrap(),
        1
    );
    assert_eq!(
        ascore.number_of_matched_ions(last, &windows[1], 1).unwrap(),
        1
    );
}

#[test]
fn cumulative_binomial_probability_of_the_best_peptide() {
    // AScore_test.cpp: calculateCumulativeBinominalProbabilityScore, with the
    // 0.5 Da tolerance the preceding section left on the shared instance.
    let ascore = with(&[("fragment_mass_tolerance", ParamValue::Float(0.5))]);
    let site = ProbablePhosphoSites {
        seq_1: 4,
        seq_2: 3,
        peak_depth: 6,
        first: 10,
        second: 7,
        ..ProbablePhosphoSites::default()
    };
    let windows = AScore::peak_picking_per_windows_in_spectrum(&ascore_input(3)).unwrap();
    let th = ascore
        .create_theoretical_spectra(&qsp_permutations(), &seq("QSSVTQVTEQSPK"))
        .unwrap();
    let ions = ascore.compute_site_determining_ions(&th, &site).unwrap();
    let n_trials = ions[0].len();
    let p = site.peak_depth as f64 / 100.0;
    let mut n_first = 0;
    for window in &windows {
        n_first += ascore
            .number_of_matched_ions(&ions[0], window, site.peak_depth)
            .unwrap();
    }
    let p_first = AScore::compute_cumulative_score(n_trials, n_first, p).unwrap();
    real_similar(-10.0 * p_first.log10(), 53.5336889240929);
}

#[test]
fn theoretical_spectra_per_placement() {
    // AScore_test.cpp: createTheoreticalSpectra_.
    let th = AScore::default()
        .create_theoretical_spectra(&qsp_permutations(), &seq("QSSVTQVTEQSPK"))
        .unwrap();
    assert_eq!(th.len(), 5);
    assert_eq!(th[0].name, "QS(Phospho)SVTQVTEQSPK");
    assert_eq!(th[4].name, "QSSVTQVTEQS(Phospho)PK");
    real_similar(th[4].peaks[0].mz, 147.11340);
    real_similar(th[4].peaks[2].mz, 244.166);
    real_similar(th[4].peaks[21].mz, 1352.57723);
}

fn ascore_1(hit: &PeptideHit) -> f64 {
    hit.metadata["AScore_1"].as_f64().unwrap()
}

fn run(ascore: &AScore, input: usize, sequence: &str) -> PeptideHit {
    let hit = PeptideHit::new(1.0, 1, 1, seq(sequence)).unwrap();
    let mut spectrum = ascore_input(input);
    ascore.compute(&hit, &mut spectrum).unwrap()
}

#[test]
fn compute_in_da() {
    // AScore_test.cpp: compute, Da cases.
    let ascore = with(&[("fragment_mass_tolerance", ParamValue::Float(0.6))]);
    let hit1 = run(&ascore, 1, "QSSVT(Phospho)QSK");
    real_similar(ascore_1(&hit1), 8.65157151899052);
    assert_eq!(hit1.sequence.to_string(), "QSS(Phospho)VTQSK");
    assert_eq!(
        hit1.metadata["search_engine_sequence"].as_str().unwrap(),
        "QSSVT(Phospho)QSK"
    );

    let hit2 = run(&ascore, 2, "RIRLT(Phospho)ATTR");
    real_similar(ascore_1(&hit2), 18.8755623850511);
    assert_eq!(hit2.sequence.to_string(), "RIRLT(Phospho)ATTR");

    let hit3 = run(&ascore, 3, "QSSVTQVTEQS(Phospho)PK");
    real_similar(ascore_1(&hit3), 88.3030731386678);
    assert_eq!(hit3.sequence.to_string(), "QSSVTQVTEQS(Phospho)PK");

    let ascore = with(&[("fragment_mass_tolerance", ParamValue::Float(0.05))]);
    let hit4 = run(&ascore, 4, "ATPGNLGSSVLHS(Phospho)K");
    real_similar(ascore_1(&hit4), 49.2714597801023);
    assert_eq!(hit4.sequence.to_string(), "ATPGNLGSSVLHS(Phospho)K");
}

#[test]
fn compute_in_ppm() {
    // AScore_test.cpp: compute, "PPM UNIT TEST" and the unambiguous case.
    let ascore = with(&[
        ("fragment_mass_tolerance", ParamValue::Float(700.0)),
        ("fragment_mass_unit", ParamValue::String("ppm".into())),
    ]);
    let hit5 = run(&ascore, 1, "QSSVT(Phospho)QSK");
    real_similar(ascore_1(&hit5), 6.53833235677545);
    assert_eq!(hit5.sequence.to_string(), "QSS(Phospho)VTQSK");

    let ascore = with(&[
        ("fragment_mass_tolerance", ParamValue::Float(70.0)),
        ("fragment_mass_unit", ParamValue::String("ppm".into())),
    ]);
    let hit6 = run(&ascore, 4, "ATPGNLGSSVLHS(Phospho)K");
    real_similar(ascore_1(&hit6), 40.6506162613816);
    assert_eq!(hit6.sequence.to_string(), "ATPGNLGSSVLHS(Phospho)K");

    let hit7 = run(&ascore, 4, "PEPT(Phospho)IDE");
    let unambiguous = ascore.to_param().unwrap();
    real_similar(
        hit7.score,
        unambiguous
            .value("unambiguous_score")
            .unwrap()
            .to_f64()
            .unwrap(),
    );
}

#[test]
fn compute_score_is_the_lowest_site_score() {
    // native: the returned score is the minimum of the AScore_<n> values and
    // AScore_pep_score is the best weighted placement score.
    let ascore = with(&[("fragment_mass_tolerance", ParamValue::Float(0.6))]);
    let hit = run(&ascore, 3, "QSSVTQVTEQS(Phospho)PK");
    assert_eq!(hit.score, ascore_1(&hit));
    assert!(hit.metadata.contains_key("AScore_pep_score"));
    assert!(!hit.metadata.contains_key("AScore_2"));
}

#[test]
fn compute_aborted_cases_reset_the_score() {
    // native, from the branches of AScore::compute.
    let ascore = AScore::default();
    let hit = PeptideHit::new(5.0, 1, 1, seq("QSSVT(Phospho)QSK")).unwrap();
    let mut empty = MSSpectrum::default();
    assert_eq!(ascore.compute(&hit, &mut empty).unwrap().score, -1.0);
    let mut spectrum = ascore_input(1);
    let unphosphorylated = PeptideHit::new(5.0, 1, 1, seq("QSSVTQSK")).unwrap();
    assert_eq!(
        ascore
            .compute(&unphosphorylated, &mut spectrum)
            .unwrap()
            .score,
        -1.0
    );
    let short = AScore {
        max_peptide_length: 3,
        ..AScore::default()
    };
    assert_eq!(short.compute(&hit, &mut spectrum).unwrap().score, -1.0);
    let few = AScore {
        max_permutations: 3,
        ..AScore::default()
    };
    let out = few.compute(&hit, &mut spectrum).unwrap();
    assert_eq!(out.score, -1.0);
    assert!(out.metadata.is_empty());
}

#[test]
fn parameters_round_trip_and_validate() {
    // native: the source defaults (AScore.cpp constructor) and restrictions.
    let defaults = AScore::defaults().unwrap();
    assert_eq!(
        defaults.value("fragment_mass_tolerance").unwrap(),
        &ParamValue::Float(0.05)
    );
    assert_eq!(
        defaults.value("fragment_mass_unit").unwrap(),
        &ParamValue::String("Da".into())
    );
    assert_eq!(
        defaults.value("max_peptide_length").unwrap(),
        &ParamValue::Integer(40)
    );
    assert_eq!(
        defaults.value("max_num_perm").unwrap(),
        &ParamValue::Integer(16384)
    );
    assert_eq!(
        defaults.value("unambiguous_score").unwrap(),
        &ParamValue::Integer(1000)
    );
    assert!(defaults.has_tag("max_num_perm", "advanced").unwrap());
    assert_eq!(AScore::from_param(&defaults).unwrap(), AScore::default());
    let custom = AScore {
        fragment_mass_tolerance: 12.5,
        fragment_tolerance_ppm: true,
        max_peptide_length: 0,
        max_permutations: 7,
        unambiguous_score: -3,
    };
    assert_eq!(
        AScore::from_param(&custom.to_param().unwrap()).unwrap(),
        custom
    );

    let mut negative = Param::new();
    negative
        .set_value("fragment_mass_tolerance", ParamValue::Float(-1.0), "", &[])
        .unwrap();
    assert!(AScore::from_param(&negative).is_err());
    let mut unit = Param::new();
    unit.set_value(
        "fragment_mass_unit",
        ParamValue::String("Th".into()),
        "",
        &[],
    )
    .unwrap();
    assert!(AScore::from_param(&unit).is_err());
    let mut unknown = Param::new();
    unknown
        .set_value("no_such_parameter", ParamValue::Integer(1), "", &[])
        .unwrap();
    let (_, warnings) = AScore::from_param_with_warnings(&unknown).unwrap();
    assert!(!warnings.is_empty());
}

#[test]
fn ascore_refuses_where_the_source_is_undefined() {
    // native.
    assert!(AScore::peptide_score(&[1.0; 9]).is_err());
    assert!(AScore::determine_highest_scoring_permutations(&[], &[], &[]).is_err());
    // No runner-up: the only other placement keeps the best site.
    let scores = rows(&[0.1, 0.2]);
    let ranking = AScore::rank_weighted_permutation_peptide_scores(&scores).unwrap();
    assert!(
        AScore::determine_highest_scoring_permutations(&scores, &[vec![1], vec![1]], &ranking)
            .is_err()
    );
    let ascore = AScore::default();
    assert!(
        ascore
            .compute_site_determining_ions(&[], &ProbablePhosphoSites::default())
            .is_err()
    );
    // More events than sites: no placement.
    assert!(AScore::compute_permutations(&[1], 2).unwrap().is_empty());
    // The placement ceiling is checked before allocation.
    let sites: Vec<usize> = (0..40).collect();
    assert!(AScore::compute_permutations(&sites, 20).is_err());
    assert!(
        AScore::compute_permutations(&sites[..21], 1).unwrap().len() <= MAX_ASCORE_PERMUTATIONS
    );
    let negative = MSSpectrum {
        peaks: vec![Peak1D {
            mz: -5.0,
            intensity: 1.0,
        }],
        ..MSSpectrum::default()
    };
    assert!(AScore::peak_picking_per_windows_in_spectrum(&negative).is_err());
    let wide = MSSpectrum {
        peaks: vec![
            Peak1D {
                mz: 1.0,
                intensity: 1.0,
            },
            Peak1D {
                mz: 1e12,
                intensity: 1.0,
            },
        ],
        ..MSSpectrum::default()
    };
    assert!(AScore::peak_picking_per_windows_in_spectrum(&wide).is_err());
    let bad_tolerance = AScore {
        fragment_mass_tolerance: f64::NAN,
        ..AScore::default()
    };
    let hit = PeptideHit::new(1.0, 1, 1, seq("QSSVT(Phospho)QSK")).unwrap();
    assert!(bad_tolerance.compute(&hit, &mut ascore_input(1)).is_err());
}

#[test]
fn spectrum_difference_skips_ions_within_one_tolerance() {
    // native, from getSpectrumDifference_: a match removes both peaks and
    // every further peak of either side within tolerance of the matched pair.
    let ascore = AScore::default(); // 0.05 Da
    let peaks = |mz: &[f64]| -> Vec<Peak1D> {
        mz.iter().map(|&mz| Peak1D { mz, intensity: 1.0 }).collect()
    };
    let first = peaks(&[100.0, 200.0, 200.04, 300.0, 400.0]);
    let second = peaks(&[200.01, 200.03, 350.0]);
    let difference = ascore.spectrum_difference(&first, &second);
    let mz: Vec<f64> = difference.iter().map(|p| p.mz).collect();
    assert_eq!(mz, vec![100.0, 300.0, 400.0]);
    assert_eq!(ascore.compare_mz(1.0, 1.04), Ordering::Equal);
    assert_eq!(ascore.compare_mz(1.0, 1.06), Ordering::Less);
    assert_eq!(ascore.compare_mz(1.06, 1.0), Ordering::Greater);
    let ppm = AScore {
        fragment_mass_tolerance: 10.0,
        fragment_tolerance_ppm: true,
        ..AScore::default()
    };
    assert_eq!(ppm.compare_mz(1000.0, 1000.009), Ordering::Equal);
    assert_eq!(ppm.compare_mz(1000.0, 1000.011), Ordering::Less);
    real_similar(
        ppm.compute_base_probability(1000.0),
        2.0 * 10.0 / 100.0 * 1000.0 * 1e-6,
    );
}
