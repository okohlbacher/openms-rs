// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `ANALYSIS/ID/Percolator.h` and `PercolatorTypes.h`.
//!
//! Each test ports one section of `Percolator_test.cpp` at `bc9cc12`; the
//! section is named in the test's comment. The source generates its
//! synthetic data with glibc `rand()` after `srand(seed)`, reproduced here by
//! [`GlibcRand`] so the inputs are the source's inputs (up to the last bit of
//! the platform `log`/`cos` in the Box-Muller helper, which the properties
//! checked do not depend on). The one literal-valued section, the TRON
//! characterization, is asserted at the source's own relative tolerance.
//! See `tests/data/percolator_provenance.json`.

use openms::Error;
use openms::analysis::id::percolator::{
    FeatureNormalizer, MainScore, PepMethod, Percolator, PercolatorParameters, load_model,
    read_model, save_model, write_model,
};
use openms::analysis::id::percolator_types::{PercolatorModel, RescoreInput, RescoreOutput};
use openms::chemistry::AASequence;
use openms::identification::{PeptideHit, PeptideIdentification};
use openms::metadata::MetaValue;
use std::f64::consts::PI;

/// glibc's `rand()` (the `TYPE_3` additive feedback generator `srand`
/// selects), so the class test's `std::srand(seed)` data is reproduced.
struct GlibcRand {
    r: Vec<u32>,
    k: usize,
}

impl GlibcRand {
    fn new(seed: u32) -> Self {
        let mut r = vec![0u32; 344];
        r[0] = if seed == 0 { 1 } else { seed };
        for i in 1..31 {
            let hi = (r[i - 1] as i32) / 127_773;
            let lo = (r[i - 1] as i32) % 127_773;
            let mut word = 16_807 * lo - 2_836 * hi;
            if word < 0 {
                word += 2_147_483_647;
            }
            r[i] = word as u32;
        }
        for i in 31..34 {
            r[i] = r[i - 31];
        }
        for i in 34..344 {
            r[i] = r[i - 31].wrapping_add(r[i - 3]);
        }
        Self { r, k: 344 }
    }

    fn next(&mut self) -> u32 {
        let value = self.r[self.k - 31].wrapping_add(self.r[self.k - 3]);
        self.r.push(value);
        self.k += 1;
        value >> 1
    }

    /// `static_cast<double>(std::rand()) / RAND_MAX`.
    fn rand01(&mut self) -> f64 {
        f64::from(self.next()) / 2_147_483_647.0
    }

    /// The test's Box-Muller `randn` (two uniform draws).
    fn randn(&mut self) -> f64 {
        let u1 = self.rand01().max(1e-9);
        let u2 = self.rand01();
        (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
    }
}

/// `makeModeratelySeparableInput_`.
fn moderately_separable(n_per_class: usize, seed: u32) -> RescoreInput {
    let mut rng = GlibcRand::new(seed);
    let mut input = RescoreInput {
        feature_names: vec!["signal".into(), "noise".into()],
        ..Default::default()
    };
    for _ in 0..n_per_class {
        let signal = 1.0 + 0.6 * rng.randn();
        input.features.push(vec![signal, rng.rand01()]);
        input.is_decoy.push(false);
    }
    for _ in 0..n_per_class {
        let signal = 0.6 * rng.randn();
        input.features.push(vec![signal, rng.rand01()]);
        input.is_decoy.push(true);
    }
    input
}

fn median_target_beats_decoy(out: &RescoreOutput, input: &RescoreInput) -> bool {
    let (mut t, mut d): (Vec<f64>, Vec<f64>) = (Vec::new(), Vec::new());
    for (i, &score) in out.scores.iter().enumerate() {
        if input.is_decoy[i] {
            d.push(score);
        } else {
            t.push(score);
        }
    }
    t.sort_by(f64::total_cmp);
    d.sort_by(f64::total_cmp);
    t[t.len() / 2] > d[d.len() / 2]
}

fn params(seed: i32) -> PercolatorParameters {
    PercolatorParameters {
        seed,
        ..PercolatorParameters::default()
    }
}

fn percolator(seed: i32) -> Percolator {
    Percolator::with_parameters(params(seed))
}

fn all_in_unit(values: &[f64]) -> bool {
    values.iter().all(|&v| (0.0..=1.0 + 1e-9).contains(&v))
}

fn descending_order(scores: &[f64]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..scores.len()).collect();
    idx.sort_by(|&a, &b| scores[b].total_cmp(&scores[a]));
    idx
}

#[test]
fn glibc_rand_matches_the_reference_sequence() {
    // glibc 2.x: srand(1); rand() x3.
    let mut rng = GlibcRand::new(1);
    assert_eq!(rng.next(), 1_804_289_383);
    assert_eq!(rng.next(), 846_930_886);
    assert_eq!(rng.next(), 1_681_692_777);
}

#[test]
fn rescore_trivially_separable() {
    // Section `RescoreOutput rescore(const RescoreInput& input)`.
    let mut input = RescoreInput {
        feature_names: vec!["separator".into(), "noise".into()],
        ..Default::default()
    };
    for i in 0..200 {
        input
            .features
            .push(vec![1.0 + 0.01 * f64::from(i), f64::from(i % 7)]);
        input.is_decoy.push(false);
    }
    for i in 0..200 {
        input
            .features
            .push(vec![-1.0 + 0.01 * f64::from(i), f64::from(i % 11)]);
        input.is_decoy.push(true);
    }
    let out = percolator(42).rescore(&input).unwrap();
    assert_eq!(out.scores.len(), 400);
    assert_eq!(out.q_values.len(), 400);
    assert_eq!(out.peps.len(), 400);
    assert!(median_target_beats_decoy(&out, &input));
}

#[test]
fn reproducibility_with_fixed_seed_is_bit_identical() {
    // Section `[EXTRA] reproducibility with fixed seed is bit-identical`.
    let mut input = RescoreInput {
        feature_names: vec!["f".into()],
        ..Default::default()
    };
    for i in 0..400 {
        let decoy = i % 2 == 1;
        input
            .features
            .push(vec![if decoy { -1.0 } else { 1.0 } + 0.01 * f64::from(i)]);
        input.is_decoy.push(decoy);
    }
    let mut p1 = Percolator::with_parameters(PercolatorParameters {
        num_threads: 1,
        ..params(17)
    });
    let mut p2 = Percolator::with_parameters(PercolatorParameters {
        num_threads: 1,
        ..params(17)
    });
    let a = p1.rescore(&input).unwrap();
    let b = p2.rescore(&input).unwrap();
    assert_eq!(a, b);
    // Native: the thread count does not change the bits.
    let mut p3 = Percolator::with_parameters(PercolatorParameters {
        num_threads: 3,
        ..params(17)
    });
    assert_eq!(p3.rescore(&input).unwrap(), a);
}

#[test]
fn q_values_in_expected_range() {
    // Section `[EXTRA] q-values in expected range`.
    let mut rng = GlibcRand::new(42);
    let mut input = RescoreInput {
        feature_names: vec!["f".into(), "noise".into()],
        ..Default::default()
    };
    for _ in 0..500 {
        let t = 1.0 + 0.6 * (rng.rand01() - 0.5) * 2.0;
        let d = 0.6 * (rng.rand01() - 0.5) * 2.0;
        let noise_t = rng.rand01();
        input.features.push(vec![t, noise_t]);
        input.is_decoy.push(false);
        let noise_d = rng.rand01();
        input.features.push(vec![d, noise_d]);
        input.is_decoy.push(true);
    }
    let out = Percolator::default().rescore(&input).unwrap();
    assert!(all_in_unit(&out.q_values));
    assert!(all_in_unit(&out.peps));
    let idx = descending_order(&out.scores);
    let n = idx.len();
    let top: f64 = idx[..50].iter().map(|&i| out.q_values[i]).sum::<f64>() / 50.0;
    let bottom: f64 = idx[n - 50..].iter().map(|&i| out.q_values[i]).sum::<f64>() / 50.0;
    assert!(top < bottom);
}

fn psm_ids(seed: u32, n: usize, td_key: &str) -> Vec<PeptideIdentification> {
    let mut rng = GlibcRand::new(seed);
    (0..n)
        .map(|i| {
            let decoy = i % 2 == 1;
            let mut hit =
                PeptideHit::new(42.0, 0, 2, AASequence::parse("PEPTIDER").unwrap()).unwrap();
            hit.metadata.insert(
                td_key.into(),
                MetaValue::from(if decoy { "decoy" } else { "target" }),
            );
            let f = if decoy { 0.0 } else { 1.5 } + 0.9 * rng.randn();
            hit.metadata
                .insert("f".into(), MetaValue::try_from(f).unwrap());
            PeptideIdentification {
                identifier: "run1".into(),
                rt: Some(0.1 * i as f64),
                hits: vec![hit],
                ..Default::default()
            }
        })
        .collect()
}

#[test]
fn rescore_peptide_identifications() {
    // Section `void rescore(std::vector<PeptideIdentification>&, const StringList&)`.
    let mut rng = GlibcRand::new(7);
    let mut peps: Vec<PeptideIdentification> = (0..800)
        .map(|i| {
            let decoy = i % 2 == 1;
            let mut hit = PeptideHit::default();
            hit.metadata.insert(
                "target_decoy".into(),
                MetaValue::from(if decoy { "decoy" } else { "target" }),
            );
            let f = if decoy { 0.0 } else { 1.0 } + 0.6 * (rng.rand01() - 0.5) * 2.0;
            hit.metadata
                .insert("feat_sep".into(), MetaValue::try_from(f).unwrap());
            hit.metadata.insert(
                "feat_noise".into(),
                MetaValue::try_from(rng.rand01()).unwrap(),
            );
            PeptideIdentification {
                identifier: "run1".into(),
                rt: Some(i as f64 * 0.1),
                hits: vec![hit],
                ..Default::default()
            }
        })
        .collect();
    Percolator::default()
        .rescore_peptide_ids(&mut peps, &["feat_sep".into(), "feat_noise".into()])
        .unwrap();
    let (mut ts, mut ds, mut nt, mut nd) = (0.0, 0.0, 0, 0);
    for pid in &peps {
        assert_eq!(pid.hits.len(), 1);
        let hit = &pid.hits[0];
        assert!(hit.metadata.contains_key("percolator_q_value"));
        assert!(hit.metadata.contains_key("percolator_pep"));
        let s = hit.metadata["percolator_score"].as_f64().unwrap();
        if hit.metadata["target_decoy"].as_str().unwrap() == "decoy" {
            ds += s;
            nd += 1;
        } else {
            ts += s;
            nt += 1;
        }
    }
    assert!(nt > 0 && nd > 0);
    assert!(ts / f64::from(nt) > ds / f64::from(nd));
}

#[test]
fn missing_target_decoy_meta_value_is_refused() {
    // Section `[EXTRA] missing target_decoy meta throws`.
    let mut peps: Vec<PeptideIdentification> = (0..20)
        .map(|_| {
            let mut hit = PeptideHit::default();
            hit.metadata
                .insert("feat_x".into(), MetaValue::try_from(1.0).unwrap());
            PeptideIdentification {
                hits: vec![hit],
                ..Default::default()
            }
        })
        .collect();
    let before = peps.clone();
    let result = Percolator::default().rescore_peptide_ids(&mut peps, &["feat_x".into()]);
    assert!(matches!(result, Err(Error::InvalidValue(_))));
    assert_eq!(peps, before);
}

#[test]
fn realistic_scale_dataset() {
    // Section `[EXTRA] output makes sense on a realistic-scale dataset`.
    let mut rng = GlibcRand::new(2026);
    let mut input = RescoreInput {
        feature_names: vec!["f0".into(), "f1".into(), "noise".into()],
        ..Default::default()
    };
    let mut easy = Vec::new();
    for _ in 0..500 {
        let row = vec![
            2.0 + 0.5 * rng.randn(),
            1.0 + 0.5 * rng.randn(),
            rng.randn(),
        ];
        input.features.push(row);
        input.is_decoy.push(false);
        easy.push(true);
    }
    for _ in 0..500 {
        let row = vec![0.5 + 0.5 * rng.randn(), 0.5 * rng.randn(), rng.randn()];
        input.features.push(row);
        input.is_decoy.push(false);
        easy.push(false);
    }
    for _ in 0..1000 {
        let row = vec![0.5 * rng.randn(), 0.0 * rng.randn(), rng.randn()];
        input.features.push(row);
        input.is_decoy.push(true);
        easy.push(false);
    }
    let out = percolator(2026).rescore(&input).unwrap();
    let n = input.features.len();
    assert_eq!(out.scores.len(), n);
    assert!(all_in_unit(&out.q_values));
    assert!(all_in_unit(&out.peps));
    let idx = descending_order(&out.scores);
    let top_targets = idx[..100].iter().filter(|&&i| !input.is_decoy[i]).count();
    let top_easy = idx[..100].iter().filter(|&&i| easy[i]).count();
    assert!(top_targets >= 80, "{top_targets}");
    assert!(top_easy >= 60, "{top_easy}");
    let bottom_non_easy = idx[n - 100..].iter().filter(|&&i| !easy[i]).count();
    assert!(bottom_non_easy >= 80);
    let confident = (0..n)
        .filter(|&i| !input.is_decoy[i] && out.q_values[i] < 0.01)
        .count();
    assert!(confident >= 50, "{confident}");
    let top_pep: f64 = idx[..100].iter().map(|&i| out.peps[i]).sum();
    let bottom_pep: f64 = idx[n - 100..].iter().map(|&i| out.peps[i]).sum();
    assert!(top_pep / 100.0 < bottom_pep / 100.0);
}

#[test]
fn report_as_main_score_promotes_an_output() {
    // Section `[EXTRA] report_as_main_score promotes a Percolator output to hit.getScore()`.
    let source = psm_ids(3, 800, "target_decoy");
    let base = PercolatorParameters {
        test_fdr: 0.1,
        train_fdr: 0.1,
        ..params(3)
    };
    let cases = [
        (MainScore::QValue, "percolator_q_value", "q-value", false),
        (
            MainScore::Pep,
            "percolator_pep",
            "Posterior Error Probability",
            false,
        ),
        (MainScore::Svm, "percolator_score", "svm", true),
    ];
    for (mode, key, score_type, higher) in cases {
        let mut peps = source.clone();
        let mut p = Percolator::with_parameters(PercolatorParameters {
            report_as_main_score: mode,
            ..base.clone()
        });
        p.rescore_peptide_ids(&mut peps, &["f".into()]).unwrap();
        for pid in &peps {
            assert_eq!(pid.score_type, score_type);
            assert_eq!(pid.higher_score_better, higher);
            let hit = &pid.hits[0];
            let promoted = hit.metadata[key].as_f64().unwrap();
            assert!((hit.score - promoted).abs() <= 1e-12);
        }
    }
    let mut peps = source.clone();
    Percolator::with_parameters(base.clone())
        .rescore_peptide_ids(&mut peps, &["f".into()])
        .unwrap();
    assert_eq!(peps[0].hits[0].score, 42.0);
    assert!(peps[0].hits[0].metadata.contains_key("percolator_q_value"));
}

#[test]
fn use_pi0_false_disables_pi0_correction() {
    // Section `[EXTRA] use_pi0=false disables pi0 correction`.
    let mut rng = GlibcRand::new(11);
    let mut input = RescoreInput {
        feature_names: vec!["f0".into()],
        ..Default::default()
    };
    for _ in 0..600 {
        input.features.push(vec![1.5 + rng.randn()]);
        input.is_decoy.push(false);
    }
    for _ in 0..600 {
        input.features.push(vec![rng.randn()]);
        input.is_decoy.push(true);
    }
    let config = PercolatorParameters {
        test_fdr: 0.05,
        train_fdr: 0.05,
        ..params(17)
    };
    let mut on = Percolator::with_parameters(config.clone());
    let mut off = Percolator::with_parameters(PercolatorParameters {
        use_pi0: false,
        ..config
    });
    assert_eq!(on.pi0(), -1.0);
    let out_on = on.rescore(&input).unwrap();
    let out_off = off.rescore(&input).unwrap();
    assert_eq!(off.pi0(), 1.0);
    assert!(on.pi0() < 1.0);
    let idx = descending_order(&out_on.scores);
    let (mut mean_on, mut mean_off, mut counted) = (0.0, 0.0, 0);
    for &i in &idx {
        if counted == 100 {
            break;
        }
        if input.is_decoy[i] {
            continue;
        }
        mean_on += out_on.q_values[i];
        mean_off += out_off.q_values[i];
        counted += 1;
    }
    assert_eq!(counted, 100);
    assert!(mean_on < mean_off);
}

#[test]
fn malformed_input_is_refused() {
    // Section `[EXTRA] malformed input throws InvalidValue`.
    let mut p = Percolator::default();
    assert!(matches!(
        p.rescore(&RescoreInput::default()),
        Err(Error::InvalidValue(_))
    ));
    let mismatch = RescoreInput {
        features: vec![vec![1.0], vec![1.0]],
        is_decoy: vec![false],
        ..Default::default()
    };
    assert!(matches!(p.rescore(&mismatch), Err(Error::InvalidValue(_))));
    let ragged = RescoreInput {
        features: vec![vec![1.0, 2.0], vec![1.0]],
        is_decoy: vec![false, true],
        ..Default::default()
    };
    assert!(matches!(p.rescore(&ragged), Err(Error::InvalidValue(_))));
    // Native: non-finite features and an unported PEP method.
    let nan = RescoreInput {
        features: vec![vec![f64::NAN], vec![1.0]],
        is_decoy: vec![false, true],
        ..Default::default()
    };
    assert!(matches!(p.rescore(&nan), Err(Error::InvalidValue(_))));
    let mut nonparametric = Percolator::with_parameters(PercolatorParameters {
        pep_method: PepMethod::Nonparametric,
        ..PercolatorParameters::default()
    });
    assert!(matches!(
        nonparametric.rescore(&moderately_separable(50, 1)),
        Err(Error::Unsupported(_))
    ));
    assert_eq!(p.pi0(), -1.0);
}

#[test]
fn train_records_the_model() {
    // Section `PercolatorModel train(const RescoreInput& input)`.
    let input = moderately_separable(300, 42);
    let mut p = percolator(17);
    let model = p.train(&input).unwrap();
    assert_eq!(model.format_version, 1);
    assert_eq!(model.feature_names, ["signal", "noise"]);
    assert_eq!(model.weights.len(), 3);
    assert_eq!(model.normalizer_type, "stdv");
    assert_eq!(model.seed, 17);
    assert!(model.weights[0] > 0.0);
    assert_eq!(p.svm_weights(), std::slice::from_ref(&model.weights));
    assert_eq!(p.pi0(), -1.0);
}

#[test]
fn score_with_a_trained_model() {
    // Sections `RescoreOutput score(...)` and `[EXTRA] train on A score disjoint B`.
    let a = moderately_separable(300, 42);
    let b = moderately_separable(300, 99);
    let mut p = percolator(17);
    let model = p.train(&a).unwrap();
    let weights = p.svm_weights().to_vec();
    let out_a = p.score(&a, &model).unwrap();
    assert_eq!(out_a.scores.len(), a.features.len());
    assert!(median_target_beats_decoy(&out_a, &a));
    let out_b = p.score(&b, &model).unwrap();
    assert_eq!(out_b.scores.len(), b.features.len());
    assert!(median_target_beats_decoy(&out_b, &b));
    // score() leaves the training weights alone.
    assert_eq!(p.svm_weights(), weights.as_slice());
}

#[test]
fn train_then_score_on_overlapping_input() {
    // Section `[EXTRA] train then score on same overlapping input`.
    let input = moderately_separable(400, 42);
    let mut p_tr = percolator(17);
    let mut p_re = percolator(17);
    let model = p_tr.train(&input).unwrap();
    let out_ts = p_tr.score(&input, &model).unwrap();
    let out_re = p_re.rescore(&input).unwrap();
    assert!(all_in_unit(&out_ts.q_values));
    assert!(all_in_unit(&out_ts.peps));
    assert!(median_target_beats_decoy(&out_ts, &input));
    let concordant = |o: &RescoreOutput| {
        let idx = descending_order(&o.scores);
        let top = idx.len() / 5;
        idx[..top].iter().filter(|&&i| !input.is_decoy[i]).count() as f64 / top as f64
    };
    let (ts, re) = (concordant(&out_ts), concordant(&out_re));
    assert!(ts > 0.70 && re > 0.70, "{ts} {re}");
    assert!((ts - re).abs() < 0.10);
}

#[test]
fn score_rejects_mismatched_models() {
    // Sections `[EXTRA] score rejects feature-count mismatch`,
    // `... feature-name mismatch` and `... empty feature_names on either side`.
    let input = moderately_separable(300, 42);
    let mut p = Percolator::default();
    let model = p.train(&input).unwrap();
    let mut fewer = input.clone();
    for row in &mut fewer.features {
        row.pop();
    }
    fewer.feature_names.pop();
    assert!(matches!(
        p.score(&fewer, &model),
        Err(Error::InvalidValue(_))
    ));
    let mut renamed = input.clone();
    renamed.feature_names = vec!["signal".into(), "different_name".into()];
    assert!(matches!(
        p.score(&renamed, &model),
        Err(Error::InvalidValue(_))
    ));
    let mut unnamed = input.clone();
    unnamed.feature_names.clear();
    assert!(matches!(
        p.score(&unnamed, &model),
        Err(Error::InvalidValue(_))
    ));
    let nameless_model = PercolatorModel {
        weights: vec![0.1; 3],
        normalizer_type: "stdv".into(),
        ..PercolatorModel::default()
    };
    assert!(matches!(
        p.score(&input, &nameless_model),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn save_and_load_round_trip() {
    // Sections `static void saveModel(...)` and `[EXTRA] feature named 'm0'
    // round-trips through save/load`.
    let dir = std::env::temp_dir().join(format!("openms_percolator_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    for first in ["signal", "m0"] {
        let mut input = moderately_separable(300, 42);
        input.feature_names[0] = first.into();
        let mut p = percolator(17);
        let model = p.train(&input).unwrap();
        let path = dir.join(format!("{first}.model"));
        save_model(&model, &path).unwrap();
        let loaded = load_model(&path).unwrap();
        // max_digits10 makes the round trip exact.
        assert_eq!(loaded, model);
        let out_mem = p.score(&input, &model).unwrap();
        let mut q = percolator(17);
        q.train(&input).unwrap();
        let out_disk = q.score(&input, &loaded).unwrap();
        assert_eq!(out_mem, out_disk);
    }
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn model_file_layout() {
    // Rust-only, derived from the saveModel source text.
    let model = PercolatorModel {
        format_version: 1,
        feature_names: vec!["f0".into(), "f1".into()],
        weights: vec![0.1, -2.0, 0.5],
        normalizer_type: "uni".into(),
        seed: 7,
    };
    let mut bytes = Vec::new();
    write_model(&model, &mut bytes).unwrap();
    assert_eq!(
        String::from_utf8(bytes).unwrap(),
        "# OpenMS Percolator model file\nformat_version: 1\nnormalizer: uni\nseed: 7\nn_features: 2\nbias: 0.5\n# features \u{2014} one 'name<TAB>weight' per line\nf0\t0.10000000000000001\nf1\t-2\n"
    );
    let broken = PercolatorModel {
        weights: vec![1.0],
        ..model.clone()
    };
    assert!(matches!(
        write_model(&broken, Vec::new()),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn load_model_is_strict() {
    // Section `static PercolatorModel loadModel(const String& filename)`.
    let missing = std::env::temp_dir().join("openms_percolator_no_such_model_file");
    assert!(matches!(load_model(&missing), Err(Error::Io(_))));
    let cases = [
        "format_version: 99\nnormalizer: stdv\nseed: 1\nn_features: 1\nbias: 0.1\nx\t0.5\n",
        "format_version: 1\nnormalizer: stdv\nseed: 1\nn_features: 1\nx\t0.5\n",
        "format_version: 1\nnormalize: stdv\nnormalizer: stdv\nseed: 1\nn_features: 1\nbias: 0.1\nx\t0.5\n",
        "format_version: 1\nnormalizer: bogus\nseed: 1\nn_features: 1\nbias: 0.1\nx\t0.5\n",
        "format_version: 1\nformat_version: 1\nnormalizer: stdv\nseed: 1\nn_features: 1\nbias: 0.1\nx\t0.5\n",
        "format_version: 1\nnormalizer: stdv\nseed: 1\nn_features: 2\nbias: 0.1\nx\t0.5\n",
    ];
    for body in cases {
        assert!(
            matches!(read_model(body.as_bytes()), Err(Error::Parse { .. })),
            "{body}"
        );
    }
    let ok = read_model(
        "# c\n\nformat_version: 1\nnormalizer: none\nseed: -3\nn_features: 1\nbias: 0.25\nx\t0.5\n"
            .as_bytes(),
    )
    .unwrap();
    assert_eq!(ok.weights, [0.5, 0.25]);
    assert_eq!(ok.seed, -3);
}

#[test]
fn same_instance_consecutive_calls_are_bit_identical() {
    // Section `[EXTRA] same-instance consecutive calls are bit-identical`.
    let input = moderately_separable(500, 123);
    let mut p = Percolator::with_parameters(PercolatorParameters {
        num_threads: 1,
        ..params(123)
    });
    let a = p.rescore(&input).unwrap();
    let b = p.rescore(&input).unwrap();
    assert_eq!(a, b);
}

#[test]
fn input_order_invariance() {
    // Section `[EXTRA] input-order invariance Jaccard characterization`.
    // The source shuffles with std::mt19937(456) + std::shuffle; any
    // permutation exercises the property, and a fixed stride permutation is
    // used here (independently derived).
    let mut base = moderately_separable(800, 456);
    base.scan_numbers = (1..=base.features.len() as i32).collect();
    let n = base.features.len();
    let perm: Vec<usize> = (0..n).map(|i| (i * 7919) % n).collect();
    let permuted = RescoreInput {
        feature_names: base.feature_names.clone(),
        features: perm.iter().map(|&i| base.features[i].clone()).collect(),
        is_decoy: perm.iter().map(|&i| base.is_decoy[i]).collect(),
        scan_numbers: perm.iter().map(|&i| base.scan_numbers[i]).collect(),
        ..Default::default()
    };
    let out1 = percolator(456).rescore(&base).unwrap();
    let out2 = percolator(456).rescore(&permuted).unwrap();
    let mut scores2 = vec![0.0; n];
    let mut q2 = vec![0.0; n];
    for (k, &i) in perm.iter().enumerate() {
        scores2[i] = out2.scores[k];
        q2[i] = out2.q_values[k];
    }
    let r = pearson(&out1.scores, &scores2);
    assert!(r >= 0.99, "{r}");
    let a: Vec<usize> = (0..n)
        .filter(|&i| !base.is_decoy[i] && out1.q_values[i] <= 0.01)
        .collect();
    let b: Vec<usize> = (0..n)
        .filter(|&i| !base.is_decoy[i] && q2[i] <= 0.01)
        .collect();
    let inter = a.iter().filter(|i| b.contains(i)).count();
    let union = a.len() + b.len() - inter;
    let jaccard = if union == 0 {
        1.0
    } else {
        inter as f64 / union as f64
    };
    assert!(jaccard >= 0.95, "{jaccard}");
}

fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let (mut sx, mut sy, mut sxx, mut syy, mut sxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for (a, b) in x.iter().zip(y) {
        sx += a;
        sy += b;
        sxx += a * a;
        syy += b * b;
        sxy += a * b;
    }
    let den = ((n * sxx - sx * sx) * (n * syy - sy * sy)).sqrt();
    if den > 1e-15 {
        (n * sxy - sx * sy) / den
    } else {
        0.0
    }
}

#[test]
fn seed_change_is_statistically_equivalent() {
    // Section `[EXTRA] seed change produces different but statistically equivalent scores`.
    let input = moderately_separable(800, 789);
    let out17 = percolator(17).rescore(&input).unwrap();
    let out42 = percolator(42).rescore(&input).unwrap();
    assert!(out17.scores.iter().zip(&out42.scores).any(|(a, b)| a != b));
    assert!(pearson(&out17.scores, &out42.scores) >= 0.9);
    let count = |o: &RescoreOutput| {
        (0..o.q_values.len())
            .filter(|&i| !input.is_decoy[i] && o.q_values[i] <= 0.05)
            .count() as i64
    };
    let (c17, c42) = (count(&out17), count(&out42));
    let maxc = c17.max(c42);
    let ratio = if maxc > 0 {
        (c17 - c42).abs() as f64 / maxc as f64
    } else {
        0.0
    };
    assert!(ratio <= 0.20);
}

fn three_feature_input(seed: u32) -> RescoreInput {
    let mut rng = GlibcRand::new(seed);
    let mut input = RescoreInput {
        feature_names: vec!["f0".into(), "f1".into(), "noise".into()],
        ..Default::default()
    };
    for _ in 0..600 {
        let row = vec![
            2.0 + 0.5 * rng.randn(),
            1.0 + 0.5 * rng.randn(),
            rng.randn(),
        ];
        input.features.push(row);
        input.is_decoy.push(false);
    }
    for _ in 0..600 {
        let row = vec![0.5 * rng.randn(), 0.5 * rng.randn(), rng.randn()];
        input.features.push(row);
        input.is_decoy.push(true);
    }
    input
}

fn separation(out: &RescoreOutput, input: &RescoreInput) -> f64 {
    let (mut st, mut sd, mut nt, mut nd) = (0.0, 0.0, 0.0, 0.0);
    for (i, &s) in out.scores.iter().enumerate() {
        if input.is_decoy[i] {
            sd += s;
            nd += 1.0;
        } else {
            st += s;
            nt += 1.0;
        }
    }
    st / nt - sd / nd
}

#[test]
fn normalizer_none_trains() {
    // Section `[EXTRA] normalizer="none" trains and produces sensible output`.
    let input = three_feature_input(2027);
    let out_none = Percolator::with_parameters(PercolatorParameters {
        normalizer: FeatureNormalizer::None,
        ..params(2027)
    })
    .rescore(&input)
    .unwrap();
    assert!(all_in_unit(&out_none.q_values));
    assert!(all_in_unit(&out_none.peps));
    assert!(separation(&out_none, &input) > 0.0);
    let out_stdv = percolator(2027).rescore(&input).unwrap();
    assert!(pearson(&out_none.scores, &out_stdv.scores) >= 0.95);
    let differing = out_none
        .scores
        .iter()
        .zip(&out_stdv.scores)
        .filter(|(a, b)| a != b)
        .count();
    assert!(differing >= input.features.len() / 10);
}

#[test]
fn initial_direction_with_and_without_inversion() {
    // Section `[EXTRA] initial_direction explicit feature with and without "-" inversion`.
    let input = three_feature_input(2028);
    let run = |direction: &str| {
        Percolator::with_parameters(PercolatorParameters {
            initial_direction: direction.into(),
            ..params(2028)
        })
        .rescore(&input)
    };
    let auto = run("").unwrap();
    let pos = run("f0").unwrap();
    let neg = run("-f0").unwrap();
    assert!(separation(&auto, &input) > 0.0);
    let pos_sep = separation(&pos, &input);
    assert!(pos_sep > 0.5);
    assert!(separation(&neg, &input) <= 0.5 * pos_sep);
    // Case-insensitive name match, as FeatureNames::getFeatureNumber.
    assert!(run("F0").is_ok());
    assert!(matches!(
        run("no_such_feature"),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn target_decoy_metavalue_overrides_the_key() {
    // Section `[EXTRA] target_decoy_metavalue overrides the meta-key for high-level rescore`.
    let peps = psm_ids(2029, 800, "my_td");
    let base = PercolatorParameters {
        test_fdr: 0.1,
        train_fdr: 0.1,
        ..params(2029)
    };
    let mut ok = peps.clone();
    Percolator::with_parameters(PercolatorParameters {
        target_decoy_metavalue: "my_td".into(),
        ..base.clone()
    })
    .rescore_peptide_ids(&mut ok, &["f".into()])
    .unwrap();
    assert!(ok[0].hits[0].metadata.contains_key("percolator_q_value"));
    let mut bad = peps;
    assert!(matches!(
        Percolator::with_parameters(base.clone()).rescore_peptide_ids(&mut bad, &["f".into()]),
        Err(Error::InvalidValue(_))
    ));
}

/// The TRON characterization fixture: 200 rows, 4 features, closed form.
fn tron_fixture() -> RescoreInput {
    let mut input = RescoreInput {
        feature_names: vec!["f0".into(), "f1".into(), "f2".into(), "f3".into()],
        ..Default::default()
    };
    for i in 0..200i32 {
        let target = i < 100;
        let sign = if target { 1.0 } else { -1.0 };
        let n0 = f64::from(((i * 7) % 17) - 8) / 30.0;
        let n1 = f64::from(((i * 11) % 13) - 6) / 30.0;
        let n2 = f64::from(((i * 13) % 11) - 5) / 30.0;
        let n3 = f64::from(((i * 17) % 19) - 9) / 30.0;
        input
            .features
            .push(vec![sign + n0, sign * 0.5 + n1, n2, n3]);
        input.is_decoy.push(!target);
    }
    input
}

fn tron_percolator() -> Percolator {
    Percolator::with_parameters(PercolatorParameters {
        seed: 1,
        c_pos: 1.0,
        c_neg: 1.0,
        num_iterations: 1,
        nested_xval_bins: 1,
        test_fdr: 0.05,
        train_fdr: 0.05,
        num_threads: 1,
        post_processing_tdc: false,
        train_best_positive: false,
        normalizer: FeatureNormalizer::Stdv,
        initial_direction: "f0".into(),
        ..PercolatorParameters::default()
    })
}

fn assert_relative(actual: f64, expected: f64, tolerance: f64, what: &str) {
    let scale = actual.abs().max(expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{what}: {actual:.17} vs {expected:.17}"
    );
}

#[test]
fn tron_solver_capture() {
    // Section `[EXTRA] TRON solver capture for BLAS-to-Eigen migration`:
    // literals and the 1e-8 relative tolerance are the source's.
    let input = tron_fixture();
    let mut p = tron_percolator();
    let model = p.train(&input).unwrap();
    assert_eq!(model.weights.len(), 5);
    let expected = [
        0.322043788444676,
        0.539442429779238,
        -0.0742443430244037,
        0.105641051671916,
        -0.408004346094467,
    ];
    for (i, (&actual, &want)) in model.weights.iter().zip(&expected).enumerate() {
        assert_relative(actual, want, 1e-8, &format!("weight {i}"));
    }
    let out = p.score(&input, &model).unwrap();
    assert_eq!(out.scores.len(), 200);
    for (row, want) in [
        (0, 5.27386027410112),
        (50, 6.7241872137246),
        (100, -1.01713704846917),
        (199, -0.648701128332891),
    ] {
        assert_relative(out.scores[row], want, 1e-8, &format!("score row {row}"));
    }
}

#[test]
fn subsampled_training_scores_every_row() {
    // Rust-only: subset_max_train below the row count trains on a sample and
    // still scores all rows (the source's second pass).
    let input = moderately_separable(400, 5);
    let mut p = Percolator::with_parameters(PercolatorParameters {
        subset_max_train: 300,
        ..params(5)
    });
    let out = p.rescore(&input).unwrap();
    assert_eq!(out.scores.len(), 800);
    assert!(median_target_beats_decoy(&out, &input));
    assert!(out.q_values.iter().any(|&q| q < 1.0));
}

#[test]
fn isotonic_and_nested_configurations_run() {
    // Rust-only: the PAVA PEP path, the uni normalizer, the Cpos/Cneg grid
    // search and nested cross-validation all produce bounded output.
    let input = moderately_separable(300, 8);
    for config in [
        Percolator::with_parameters(PercolatorParameters {
            pep_method: PepMethod::Isotonic,
            ..params(8)
        }),
        Percolator::with_parameters(PercolatorParameters {
            normalizer: FeatureNormalizer::Uni,
            c_pos: 0.0,
            ..params(8)
        }),
        Percolator::with_parameters(PercolatorParameters {
            nested_xval_bins: 2,
            num_iterations: 3,
            ..params(8)
        }),
    ] {
        let mut p = config;
        let out = p.rescore(&input).unwrap();
        assert!(all_in_unit(&out.q_values));
        assert!(all_in_unit(&out.peps));
        assert!(median_target_beats_decoy(&out, &input));
    }
}

#[test]
fn parameters_round_trip() {
    // Rust-only: defaults() lists the 17 source parameters and from_param
    // reads them back.
    let defaults = Percolator::defaults().unwrap();
    let names = [
        "c_pos",
        "c_neg",
        "test_fdr",
        "train_fdr",
        "num_iterations",
        "initial_direction",
        "pep_method",
        "seed",
        "target_decoy_metavalue",
        "normalizer",
        "train_best_positive",
        "post_processing_tdc",
        "nested_xval_bins",
        "subset_max_train",
        "report_as_main_score",
        "num_threads",
        "use_pi0",
    ];
    for name in names {
        assert!(defaults.value(name).is_ok(), "{name}");
    }
    let custom = Percolator::with_parameters(PercolatorParameters {
        seed: -4,
        pep_method: PepMethod::Isotonic,
        normalizer: FeatureNormalizer::Uni,
        report_as_main_score: MainScore::Pep,
        use_pi0: false,
        subset_max_train: 0,
        ..PercolatorParameters::default()
    });
    let back = Percolator::from_param(&custom.to_param().unwrap()).unwrap();
    assert_eq!(back.parameters.seed, -4);
    assert_eq!(back.parameters.pep_method, PepMethod::Isotonic);
    assert_eq!(back.parameters.normalizer, FeatureNormalizer::Uni);
    assert_eq!(back.parameters.report_as_main_score, MainScore::Pep);
    assert!(!back.parameters.use_pi0);
    assert_eq!(back.parameters.subset_max_train, 0);
}

#[test]
fn fill_pin_compatible_fields_follows_the_pin_writer() {
    // Rust-only, derived from fillPINCompatibleFields' source text.
    let hit = |charge: i32| {
        PeptideHit::new(1.0, 0, charge, AASequence::parse("PEPTIDE").unwrap()).unwrap()
    };
    let mut a = PeptideIdentification {
        mz: Some(500.5),
        hits: vec![hit(2), hit(3)],
        ..Default::default()
    };
    a.set_spectrum_reference("scan=17");
    a.metadata
        .insert("file_origin".into(), MetaValue::from("a.mzML"));
    let mut b = PeptideIdentification {
        hits: vec![hit(2)],
        ..Default::default()
    };
    b.set_spectrum_reference("scan=42");
    b.metadata
        .insert("file_origin".into(), MetaValue::from("b.mzML"));
    b.hits[0]
        .metadata
        .insert("CalcMass".into(), MetaValue::try_from(1234.5).unwrap());
    let empty = PeptideIdentification::default();
    let ids = [a.clone(), empty, b];
    let mut input = RescoreInput::default();
    Percolator::fill_pin_compatible_fields(&ids, true, &mut input).unwrap();
    assert_eq!(input.scan_numbers, [17, 17, 42]);
    assert_eq!(input.spec_file_numbers, [0, 0, 1]);
    assert_eq!(input.exp_masses, [500.5, 500.5, 0.0]);
    assert_eq!(input.calc_masses[2], 1234.5);
    let mz2 = a.hits[0].sequence.mz(2).unwrap();
    assert_eq!(input.calc_masses[0], mz2);
    let mut first_only = RescoreInput::default();
    Percolator::fill_pin_compatible_fields(&ids, false, &mut first_only).unwrap();
    assert_eq!(first_only.scan_numbers, [17, 42]);
}
