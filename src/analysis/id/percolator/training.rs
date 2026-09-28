// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//
// Translated from Percolator (https://github.com/percolator/percolator,
// commit eb157f7, as vendored by OpenMS), Copyright 2006-2012 Lukas Käll and
// the Percolator authors, licensed under the Apache License, Version 2.0; see
// docs/PERCOLATOR_SUPPORT.md for the notice.

//! Semi-supervised training: the feature normalizers (`StdvNormalizer`,
//! `UniNormalizer`, `NoNormalizer`), the initial-direction logic of
//! `SanityCheck`, and `CrossValidation` (three folds, optional nested
//! grid search over the soft-margin costs, fold-averaged raw weights).
//!
//! Only the configuration the OpenMS wrapper constructs is ported:
//! `quickValidation`, per-iteration test reporting, RESET decoy
//! subsampling (`decoyFractionTraining < 1`), weight-file input and
//! `pin.xml` default weights are unreachable from `Percolator.h` and are
//! absent.

use super::engine::Psm;
use super::engine::{AlgIn, PseudoRandom, Scores};
use super::svm::l2_svm;
use crate::concept::parallel::{Threads, map_collect};
use crate::{Error, Result};

/// `Normalizer::subclass_type`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum NormalizerKind {
    Stdv,
    Uni,
    None,
}

/// A fitted normalizer: `(x - sub) / div` per feature.
#[derive(Clone, Debug)]
pub(super) struct Normalizer {
    kind: NormalizerKind,
    sub: Vec<f64>,
    div: Vec<f64>,
}

impl Normalizer {
    /// `setSet` over the rows (targets first, then decoys, as the source's
    /// data sets are visited).
    pub(super) fn fit(kind: NormalizerKind, rows: &[&[f64]], n_features: usize) -> Self {
        let mut sub = vec![0.0; n_features];
        let mut div = vec![0.0; n_features];
        match kind {
            NormalizerKind::Stdv => {
                let mut n = 0.0;
                for row in rows {
                    n += 1.0;
                    for ix in 0..n_features {
                        sub[ix] += row[ix];
                    }
                }
                for value in &mut sub {
                    if n > 0.0 {
                        *value /= n;
                    }
                }
                for row in rows {
                    for ix in 0..n_features {
                        let d = row[ix] - sub[ix];
                        div[ix] += d * d;
                    }
                }
                for value in &mut div {
                    if *value <= 0.0 || n == 0.0 {
                        *value = 1.0;
                    } else {
                        *value = (*value / n).sqrt();
                    }
                }
            }
            NormalizerKind::Uni => {
                let mut mins = vec![1e100; n_features];
                let mut maxs = vec![-1e100; n_features];
                for row in rows {
                    for ix in 0..n_features {
                        // std::min(features[ix], mins[ix]) keeps the
                        // feature unless the minimum is strictly smaller;
                        // std::max(features[ix], maxs[ix]) likewise.
                        if !(mins[ix] < row[ix]) {
                            mins[ix] = row[ix];
                        }
                        if !(row[ix] < maxs[ix]) {
                            maxs[ix] = row[ix];
                        }
                    }
                }
                for ix in 0..n_features {
                    sub[ix] = mins[ix];
                    div[ix] = maxs[ix] - mins[ix];
                    if div[ix] <= 0.0 {
                        div[ix] = 1.0;
                    }
                }
            }
            NormalizerKind::None => {
                div.iter_mut().for_each(|v| *v = 1.0);
            }
        }
        Self { kind, sub, div }
    }

    /// `normalize`: in place, `(x - sub) / div`.
    pub(super) fn normalize(&self, row: &mut [f64]) {
        for ix in 0..self.sub.len() {
            row[ix] = (row[ix] - self.sub[ix]) / self.div[ix];
        }
    }

    /// `unnormalizeweight`: weights for normalized features to weights for
    /// raw features, bias last.
    pub(super) fn unnormalize_weights(&self, input: &[f64]) -> Vec<f64> {
        let n = self.sub.len();
        let mut out = vec![0.0; n + 1];
        match self.kind {
            NormalizerKind::None => out.copy_from_slice(&input[..n + 1]),
            NormalizerKind::Stdv | NormalizerKind::Uni => {
                let mut sum = 0.0;
                for i in 0..n {
                    out[i] = input[i] / self.div[i];
                    sum += self.sub[i] * input[i] / self.div[i];
                }
                out[n] = input[n] - sum;
            }
        }
        out
    }
}

/// The initial-direction state of `SanityCheck`.
pub(super) struct SanityCheck {
    /// `initDefaultDir`: 0, or a one-based feature index, negative when lower
    /// values are better.
    init_default_dir: i64,
    init_positives: i32,
    test_fdr: f64,
    initial_train_fdr: f64,
}

impl SanityCheck {
    /// `setInitDefaultDir(0)`, `setInitDefaultDirName(name)`,
    /// `checkAndSetDefaultDir()`.
    ///
    /// A leading `-` means lower is better. Names match the feature names
    /// case-insensitively (ASCII), first match wins.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a non-empty name matches no feature (the
    /// source's `MyException`, rethrown by the wrapper as
    /// `Exception::InvalidValue`).
    pub(super) fn new(initial_direction: &str, feature_names: &[String]) -> Result<Self> {
        let mut init_default_dir = 0;
        if !initial_direction.is_empty() {
            let (sign, name) = match initial_direction.strip_prefix('-') {
                Some(rest) => (-1, rest),
                None => (1, initial_direction),
            };
            let number = feature_names
                .iter()
                .position(|f| f.len() == name.len() && f.eq_ignore_ascii_case(name))
                .map_or(0, |i| i as i64 + 1);
            init_default_dir = sign * number;
            if init_default_dir == 0 {
                return Err(Error::InvalidValue(format!(
                    "initial_direction not resolved: Initial direction feature name \"{name}\" not found"
                )));
            }
        }
        Ok(Self {
            init_default_dir,
            init_positives: 0,
            test_fdr: 0.01,
            initial_train_fdr: 0.01,
        })
    }

    /// `getDefaultDirection`: the named feature's unit vector, or per fold
    /// the best single feature of that fold's training set.
    fn default_direction(
        &self,
        w: &mut [Vec<f64>],
        train: &mut [Scores],
        n_features: usize,
        psms: &[Psm],
    ) -> Result<()> {
        if self.init_default_dir == 0 {
            for (set, weights) in w.iter_mut().enumerate() {
                train[set].init_direction(self.initial_train_fdr, weights, n_features, psms)?;
            }
        } else {
            let index = (self.init_default_dir.unsigned_abs() - 1) as usize;
            let value = if self.init_default_dir < 0 { -1.0 } else { 1.0 };
            for weights in w.iter_mut() {
                weights.iter_mut().for_each(|v| *v = 0.0);
                weights[index] = value;
            }
        }
        Ok(())
    }

    /// `getInitDirection`: set the starting weights and count the test-set
    /// targets they separate.
    #[allow(clippy::too_many_arguments)]
    fn init_direction(
        &mut self,
        test: &mut [Scores],
        train: &mut [Scores],
        w: &mut [Vec<f64>],
        test_fdr: f64,
        initial_train_fdr: f64,
        n_features: usize,
        psms: &[Psm],
    ) -> Result<i32> {
        self.test_fdr = test_fdr;
        self.initial_train_fdr = initial_train_fdr;
        self.default_direction(w, train, n_features, psms)?;
        self.init_positives = 0;
        for (set, weights) in w.iter().enumerate() {
            self.init_positives +=
                test[set].calc_scores_and_qvals(weights, test_fdr, false, psms)?;
        }
        Ok(self.init_positives)
    }

    /// `validateDirection`: reset to the default direction when training
    /// found no targets, fewer than the initial direction, or flipped the
    /// named feature's sign.
    fn validate_direction(
        &self,
        test: &mut [Scores],
        train: &mut [Scores],
        w: &mut [Vec<f64>],
        n_features: usize,
        psms: &[Psm],
    ) -> Result<bool> {
        let mut over_fdr = 0;
        for (set, weights) in w.iter().enumerate() {
            over_fdr += test[set].calc_scores_and_qvals(weights, self.test_fdr, false, psms)?;
        }
        if over_fdr <= 0 || self.init_positives > over_fdr {
            self.default_direction(w, train, n_features, psms)?;
            return Ok(false);
        }
        if self.init_default_dir != 0 {
            let index = (self.init_default_dir.unsigned_abs() - 1) as usize;
            for weights in w.iter() {
                if weights[index] * self.init_default_dir as f64 <= 0.0 {
                    self.default_direction(w, train, n_features, psms)?;
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }
}

/// `CandidateCposCfrac`.
#[derive(Clone, Debug)]
struct Candidate {
    cpos: f64,
    cfrac: f64,
    set: usize,
    nested_set: usize,
    ww: Vec<f64>,
}

/// `CrossValidation` as the wrapper configures it.
pub(super) struct CrossValidation {
    test_fdr: f64,
    selection_fdr: f64,
    initial_selection_fdr: f64,
    selected_cpos: f64,
    selected_cneg: f64,
    niter: u32,
    use_pi0: bool,
    nested_xval_bins: usize,
    train_best_positive: bool,
    threads: Threads,
    num_folds: usize,
    n_features: usize,
    svm_inputs: Vec<AlgIn>,
    weights: Vec<Vec<f64>>,
    candidates: Vec<Candidate>,
    train_scores: Vec<Scores>,
    test_scores: Vec<Scores>,
    candidates_cpos: Vec<f64>,
    candidates_cfrac: Vec<f64>,
}

/// The parameters of one cross-validation run.
pub(super) struct CvConfig {
    pub(super) test_fdr: f64,
    pub(super) train_fdr: f64,
    pub(super) c_pos: f64,
    pub(super) c_neg: f64,
    pub(super) niter: u32,
    pub(super) use_pi0: bool,
    pub(super) nested_xval_bins: usize,
    pub(super) train_best_positive: bool,
    pub(super) threads: Threads,
}

impl CrossValidation {
    pub(super) fn new(config: &CvConfig, n_features: usize) -> Self {
        let mut selection_fdr = config.train_fdr;
        let mut initial_selection_fdr = config.train_fdr;
        if selection_fdr <= 0.0 {
            selection_fdr = config.test_fdr;
            if initial_selection_fdr <= 0.0 {
                initial_selection_fdr = config.test_fdr;
            }
        }
        Self {
            test_fdr: config.test_fdr,
            selection_fdr,
            initial_selection_fdr,
            selected_cpos: config.c_pos,
            selected_cneg: config.c_neg,
            niter: config.niter,
            use_pi0: config.use_pi0,
            nested_xval_bins: config.nested_xval_bins,
            train_best_positive: config.train_best_positive,
            threads: config.threads,
            num_folds: 3,
            n_features,
            svm_inputs: Vec::new(),
            weights: Vec::new(),
            candidates: Vec::new(),
            train_scores: Vec::new(),
            test_scores: Vec::new(),
            candidates_cpos: Vec::new(),
            candidates_cfrac: Vec::new(),
        }
    }

    /// `preIterationSetup`: folds, initial direction, grid.
    pub(super) fn pre_iteration_setup(
        &mut self,
        fullset: &mut Scores,
        sanity: &mut SanityCheck,
        psms: &[Psm],
        rng: &mut PseudoRandom,
    ) -> Result<i32> {
        self.svm_inputs = (0..self.num_folds * self.nested_xval_bins)
            .map(|_| AlgIn::default())
            .collect();
        fullset.create_xval_sets_by_spectrum(
            &mut self.train_scores,
            &mut self.test_scores,
            self.num_folds,
            psms,
            rng,
        )?;
        self.weights = vec![vec![0.0; self.n_features + 1]; self.num_folds];
        let positives = sanity.init_direction(
            &mut self.test_scores,
            &mut self.train_scores,
            &mut self.weights,
            self.test_fdr,
            self.initial_selection_fdr,
            self.n_features,
            psms,
        )?;
        self.initialize_grid_search(fullset.target_decoy_size_ratio());
        Ok(positives)
    }

    /// `initializeGridSearch`.
    fn initialize_grid_search(&mut self, target_decoy_size_ratio: f64) {
        if self.selected_cpos > 0.0 {
            self.candidates_cpos.push(self.selected_cpos);
        } else {
            self.candidates_cpos.extend([10.0, 1.0, 0.1]);
        }
        if self.selected_cpos > 0.0 && self.selected_cneg > 0.0 {
            self.candidates_cfrac
                .push(self.selected_cneg / self.selected_cpos);
        } else {
            self.candidates_cfrac.extend([
                1.0 * target_decoy_size_ratio,
                3.0 * target_decoy_size_ratio,
                10.0 * target_decoy_size_ratio,
            ]);
        }
        for set in 0..self.num_folds {
            for nested_set in 0..self.nested_xval_bins {
                for &cpos in &self.candidates_cpos {
                    for &cfrac in &self.candidates_cfrac {
                        self.candidates.push(Candidate {
                            cpos,
                            cfrac,
                            set,
                            nested_set,
                            ww: vec![0.0; self.n_features + 1],
                        });
                    }
                }
            }
        }
    }

    /// `train`: `niter` cross-validation steps, then score each test fold.
    pub(super) fn train(&mut self, psms: &[Psm], rng: &mut PseudoRandom) -> Result<()> {
        for i in 0..self.niter {
            let fdr = if i == 0 {
                self.initial_selection_fdr
            } else {
                self.selection_fdr
            };
            self.do_step(fdr, psms, rng)?;
        }
        for set in 0..self.num_folds {
            self.test_scores[set].calc_scores_and_qvals(
                &self.weights[set],
                self.test_fdr,
                false,
                psms,
            )?;
        }
        Ok(())
    }

    /// `doStep`.
    fn do_step(&mut self, selection_fdr: f64, psms: &[Psm], rng: &mut PseudoRandom) -> Result<i32> {
        for set in 0..self.num_folds {
            self.train_scores[set].calc_scores_and_qvals(
                &self.weights[set],
                selection_fdr,
                true,
                psms,
            )?;
        }
        let mut nested_test_scores: Vec<Vec<Scores>> = Vec::with_capacity(self.num_folds);
        for set in 0..self.num_folds {
            let mut nested_train: Vec<Scores> = (0..self.nested_xval_bins)
                .map(|_| Scores::new(self.use_pi0))
                .collect();
            let mut nested_test: Vec<Scores> = (0..self.nested_xval_bins)
                .map(|_| Scores::new(self.use_pi0))
                .collect();
            if self.nested_xval_bins > 1 {
                self.train_scores[set].create_xval_sets_by_spectrum(
                    &mut nested_train,
                    &mut nested_test,
                    self.nested_xval_bins,
                    psms,
                    rng,
                )?;
            } else {
                nested_train[0] = self.train_scores[set].clone();
                nested_test[0] = self.train_scores[set].clone();
            }
            for (nested_fold, train) in nested_train.iter_mut().enumerate() {
                let input = &mut self.svm_inputs[set * self.nested_xval_bins + nested_fold];
                train.generate_negative_training_set(input);
                train.generate_positive_training_set(
                    input,
                    selection_fdr,
                    self.train_best_positive,
                    psms,
                )?;
            }
            nested_test_scores.push(nested_test);
        }
        let inputs = &self.svm_inputs;
        let nested = self.nested_xval_bins;
        let n_features = self.n_features;
        let trained: Vec<Result<Vec<f64>>> =
            map_collect(&self.candidates, self.threads, |candidate| {
                train_candidate(
                    &inputs[candidate.set * nested + candidate.nested_set],
                    candidate,
                    n_features,
                    psms,
                )
            });
        for (candidate, weights) in self.candidates.iter_mut().zip(trained) {
            candidate.ww = weights?;
        }
        self.merge_candidates(selection_fdr, &mut nested_test_scores, psms)
    }

    /// `mergeCpCnPairs`: pick each fold's best (cpos, cfrac) by test-set
    /// positives (`>=`, so later candidates win ties), retraining on the
    /// whole fold when nested cross-validation chose it.
    fn merge_candidates(
        &mut self,
        selection_fdr: f64,
        nested_test_scores: &mut [Vec<Scores>],
        psms: &[Psm],
    ) -> Result<i32> {
        let folds = self.num_folds;
        let mut best_true_poses = vec![0i32; folds];
        let mut best_cposes = vec![1.0; folds];
        let mut best_cfracs = vec![1.0; folds];
        let per_set = self.candidates.len() / folds;
        for set in 0..folds {
            let mut intermediate: Vec<((f64, f64), i32)> = Vec::new();
            for index in set * per_set..(set + 1) * per_set {
                let candidate = &self.candidates[index];
                let tp = nested_test_scores[set][candidate.nested_set].calc_scores_and_qvals(
                    &candidate.ww,
                    self.test_fdr,
                    true,
                    psms,
                )?;
                let key = (candidate.cpos, candidate.cfrac);
                match intermediate.iter_mut().find(|(k, _)| *k == key) {
                    Some((_, total)) => *total += tp,
                    None => intermediate.push((key, tp)),
                }
                if self.nested_xval_bins <= 1 && tp >= best_true_poses[set] {
                    best_true_poses[set] = tp;
                    self.weights[set] = candidate.ww.clone();
                    best_cposes[set] = candidate.cpos;
                    best_cfracs[set] = candidate.cfrac;
                }
            }
            if self.nested_xval_bins > 1 {
                for &cpos in &self.candidates_cpos {
                    for &cfrac in &self.candidates_cfrac {
                        let tp = intermediate
                            .iter()
                            .find(|(k, _)| *k == (cpos, cfrac))
                            .map_or(0, |(_, v)| *v);
                        if tp >= best_true_poses[set] {
                            best_true_poses[set] = tp;
                            best_cposes[set] = cpos;
                            best_cfracs[set] = cfrac;
                        }
                    }
                }
            }
        }
        if self.nested_xval_bins > 1 {
            for set in 0..folds {
                let index = set * self.nested_xval_bins;
                let mut input = std::mem::take(&mut self.svm_inputs[index]);
                self.train_scores[set].generate_negative_training_set(&mut input);
                self.train_scores[set].generate_positive_training_set(
                    &mut input,
                    selection_fdr,
                    self.train_best_positive,
                    psms,
                )?;
                let rows: Vec<&[f64]> = input
                    .vals
                    .iter()
                    .map(|&i| psms[i].features.as_slice())
                    .collect();
                self.weights[set] = l2_svm(
                    &rows,
                    &input.y,
                    self.n_features,
                    best_cposes[set],
                    best_cposes[set] * best_cfracs[set],
                )?;
                self.svm_inputs[index] = input;
            }
        }
        let mut best_true_pos = 0.0;
        for set in 0..folds {
            best_true_pos += f64::from(self.train_scores[set].calc_scores_and_qvals(
                &self.weights[set],
                self.test_fdr,
                false,
                psms,
            )?);
        }
        Ok((best_true_pos / (folds - 1) as f64) as i32)
    }

    /// `postIterationProcessing`: validate the direction (resetting it and
    /// rescoring every test fold with fold 0's weights when it fails), then
    /// merge the test folds into `fullset`.
    pub(super) fn post_iteration_processing(
        &mut self,
        fullset: &mut Scores,
        sanity: &SanityCheck,
        psms: &[Psm],
        rng: &mut PseudoRandom,
    ) -> Result<()> {
        let valid = sanity.validate_direction(
            &mut self.test_scores,
            &mut self.train_scores,
            &mut self.weights,
            self.n_features,
            psms,
        )?;
        if !valid {
            for set in 0..self.num_folds {
                self.test_scores[set].calc_scores_and_qvals(
                    &self.weights[0],
                    self.selection_fdr,
                    false,
                    psms,
                )?;
            }
        }
        fullset.merge(
            &mut self.test_scores,
            self.selection_fdr,
            false,
            &mut self.weights,
            psms,
            rng,
        )
    }

    /// `getAvgWeights`: the mean of the folds' raw-space weights.
    pub(super) fn average_weights(&self, normalizer: &Normalizer) -> Vec<f64> {
        let n = self.n_features + 1;
        let mut weights = vec![0.0; n];
        let folds = self.weights.len() as f64;
        for set_weights in &self.weights {
            let raw = normalizer.unnormalize_weights(set_weights);
            for ix in 0..n {
                weights[ix] += raw[ix] / folds;
            }
        }
        weights
    }
}

/// `trainCpCnPair`: one SVM for one (cpos, cfrac) candidate.
fn train_candidate(
    input: &AlgIn,
    candidate: &Candidate,
    n_features: usize,
    psms: &[Psm],
) -> Result<Vec<f64>> {
    let rows: Vec<&[f64]> = input
        .vals
        .iter()
        .map(|&i| psms[i].features.as_slice())
        .collect();
    l2_svm(
        &rows,
        &input.y,
        n_features,
        candidate.cpos,
        candidate.cfrac * candidate.cpos,
    )
}
