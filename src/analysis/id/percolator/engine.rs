// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//
// Translated from Percolator (https://github.com/percolator/percolator,
// commit eb157f7, as vendored by OpenMS under
// src/openms/thirdparty/percolator), Copyright 2006-2012 Lukas Käll,
// licensed under the Apache License, Version 2.0; see
// docs/PERCOLATOR_SUPPORT.md for the notice.

//! The score bookkeeping of the vendored Percolator: `PseudoRandom`,
//! `PSMDescription`, `ScoreHolder` and its orderings, `Scores`, the
//! `PosteriorEstimator` q-value, p-value and pi0 routines, and the
//! OpenMS-local `InProcessSampler`.
//!
//! Every `std::sort` of the source goes through
//! [`source_sort_by`], because the order of equal elements it leaves is
//! observable (it decides cross-validation folds, TDC survivors and tie
//! breaks), and every `std::unique` through [`std_unique`], which keeps the
//! libstdc++ tail. The vendored code's `NO_TERMINATE` flag is always set by
//! the OpenMS wrapper, so its `MyException` paths are the warnings they
//! become there; only `std::vector::at` range errors remain errors.

use super::pep::InferPep;
use crate::math::libstdcxx::lower_bound;
use crate::math::source_sort::source_sort_by;
use crate::{Error, Result};

/// `PseudoRandom::kRandMax`, the modulus of the generator.
pub(super) const K_RAND_MAX: u64 = 4_294_967_291;

/// Percolator's linear congruential generator (`PseudoRandom.cpp`).
///
/// The source keeps the state in a process-wide static; here the state is
/// owned by one training or scoring call.
#[derive(Clone, Debug)]
pub(super) struct PseudoRandom {
    seed: u64,
}

impl PseudoRandom {
    /// `PseudoRandom::setSeed(static_cast<unsigned long>(seed))`: a negative
    /// `int` seed sign-extends to 64 bits, as the source's conversion does.
    pub(super) fn new(seed: i32) -> Self {
        Self {
            seed: i64::from(seed) as u64,
        }
    }

    /// `lcg_rand`: `seed = (seed * 279470273) % 4294967291` in 64-bit
    /// unsigned arithmetic, wrapping as `uint64_t` does.
    pub(super) fn lcg_rand(&mut self) -> u64 {
        self.seed = self.seed.wrapping_mul(279_470_273) % K_RAND_MAX;
        self.seed
    }

    /// `lcg_uniform_rand`: `lcg_rand() / (kRandMax + 1.0)`.
    pub(super) fn lcg_uniform_rand(&mut self) -> f64 {
        self.lcg_rand() as f64 / (K_RAND_MAX as f64 + 1.0)
    }
}

/// `LabelType`, ordered by its integer value.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum Label {
    Decoy,
    Target,
}

/// The fields of `PSMDescription` the in-process path reads.
#[derive(Clone, Debug)]
pub(super) struct Psm {
    /// Feature row; normalized in place before training, as the source's
    /// `FeatureMemoryPool` rows are.
    pub(super) features: Vec<f64>,
    pub(super) exp_mass: f64,
    pub(super) scan: u32,
    pub(super) spec_file: u32,
    /// The input row this PSM was built from (the source's `row_%08zu` id).
    pub(super) row: usize,
}

/// `ScoreHolder`: a score, its q-value and PEP, and the PSM it belongs to
/// (an index into the caller's PSM arena).
#[derive(Clone, Debug)]
pub(super) struct ScoreHolder {
    pub(super) score: f64,
    pub(super) q: f64,
    pub(super) pep: f64,
    pub(super) psm: usize,
    pub(super) label: Label,
}

impl ScoreHolder {
    pub(super) fn new(psm: usize, label: Label) -> Self {
        Self {
            score: 0.0,
            q: 0.0,
            pep: 0.0,
            psm,
            label,
        }
    }
    pub(super) fn is_target(&self) -> bool {
        self.label == Label::Target
    }
    pub(super) fn is_decoy(&self) -> bool {
        self.label == Label::Decoy
    }
}

/// `operator>(ScoreHolder, ScoreHolder)`: score, then scan, then experimental
/// mass, then label.
fn greater(a: &ScoreHolder, b: &ScoreHolder, psms: &[Psm]) -> bool {
    let (pa, pb) = (&psms[a.psm], &psms[b.psm]);
    (a.score > b.score)
        || (a.score == b.score && pa.scan > pb.scan)
        || (a.score == b.score && pa.scan == pb.scan && pa.exp_mass > pb.exp_mass)
        || (a.score == b.score
            && pa.scan == pb.scan
            && pa.exp_mass == pb.exp_mass
            && a.label > b.label)
}

/// `operator<(ScoreHolder, ScoreHolder)`.
fn less(a: &ScoreHolder, b: &ScoreHolder, psms: &[Psm]) -> bool {
    let (pa, pb) = (&psms[a.psm], &psms[b.psm]);
    (a.score < b.score)
        || (a.score == b.score && pa.scan < pb.scan)
        || (a.score == b.score && pa.scan == pb.scan && pa.exp_mass < pb.exp_mass)
        || (a.score == b.score
            && pa.scan == pb.scan
            && pa.exp_mass == pb.exp_mass
            && a.label < b.label)
}

/// `fast_uint_hash` on `unsigned int`.
fn fast_uint_hash(x: u32) -> u32 {
    let mut x = ((x >> 16) ^ x).wrapping_mul(0x045d_9f3b);
    x = ((x >> 16) ^ x).wrapping_mul(0x045d_9f3b);
    (x >> 16) ^ x
}

/// `OrderScanHash`'s key: `fast_uint_hash(specFileNr) ^ (fast_uint_hash(scan)
/// << 1)`, shifted and combined in 32 bits before the widening to `size_t`.
fn scan_hash(psm: &Psm) -> u64 {
    u64::from(fast_uint_hash(psm.spec_file) ^ (fast_uint_hash(psm.scan) << 1))
}

/// `OrderScanMassCharge`.
fn order_scan_mass_charge(a: &ScoreHolder, b: &ScoreHolder, psms: &[Psm]) -> bool {
    let (pa, pb) = (&psms[a.psm], &psms[b.psm]);
    if pa.spec_file != pb.spec_file {
        return pa.spec_file < pb.spec_file;
    }
    if pa.scan != pb.scan {
        return pa.scan < pb.scan;
    }
    if pa.exp_mass != pb.exp_mass {
        return pa.exp_mass < pb.exp_mass;
    }
    a.score > b.score
}

/// `OrderScanLabel`.
fn order_scan_label(a: &ScoreHolder, b: &ScoreHolder, psms: &[Psm]) -> bool {
    let (pa, pb) = (&psms[a.psm], &psms[b.psm]);
    if pa.spec_file != pb.spec_file {
        return pa.spec_file < pb.spec_file;
    }
    if pa.scan != pb.scan {
        return pa.scan < pb.scan;
    }
    a.label > b.label
}

/// `std::unique` as libstdc++ implements it (`__unique`): the kept elements
/// are compacted to the front and the new end is returned; the tail keeps
/// whatever the compaction did not overwrite, which the source observes where
/// it does not erase it (`generatePositiveTrainingSet`).
pub(super) fn std_unique<T: Clone>(items: &mut [T], equal: impl Fn(&T, &T) -> bool) -> usize {
    let n = items.len();
    let mut found = None;
    for i in 0..n.saturating_sub(1) {
        if equal(&items[i], &items[i + 1]) {
            found = Some(i);
            break;
        }
    }
    let Some(mut dest) = found else {
        return n;
    };
    let mut first = dest + 1;
    loop {
        first += 1;
        if first >= n {
            break;
        }
        if !equal(&items[dest], &items[first]) {
            dest += 1;
            items[dest] = items[first].clone();
        }
    }
    dest + 1
}

/// `std::max(a, b)`: `a < b ? b : a`.
pub(super) fn std_max(a: f64, b: f64) -> f64 {
    if a < b { b } else { a }
}

/// `std::min(a, b)`: `b < a ? b : a`.
pub(super) fn std_min(a: f64, b: f64) -> f64 {
    if b < a { b } else { a }
}

/// `mymin(a, b)`: `a > b ? b : a`, the `partial_sum` operator of the source's
/// q-value conversions.
fn mymin(a: f64, b: f64) -> f64 {
    if a > b { b } else { a }
}

/// `std::partial_sum(v.rbegin(), v.rend(), v.rbegin(), mymin)`.
fn reverse_running_min(values: &mut [f64]) {
    let n = values.len();
    if n == 0 {
        return;
    }
    let mut acc = values[n - 1];
    for index in (0..n - 1).rev() {
        acc = mymin(acc, values[index]);
        values[index] = acc;
    }
}

/// `std::vector::at`: an out-of-range index is the source's
/// `std::out_of_range`, which the wrapper reports as a failed run.
fn at(values: &[f64], index: i64) -> Result<f64> {
    usize::try_from(index)
        .ok()
        .and_then(|i| values.get(i).copied())
        .ok_or_else(|| Error::InvalidValue("vector::_M_range_check: index out of range".into()))
}

/// Mix-max counts (`PosteriorEstimator::getMixMaxCounts`): walking from the
/// worst score, one entry per decoy of the running target and decoy counts,
/// ties flushed together.
fn mix_max_counts(combined: &[(f64, bool)]) -> (Vec<f64>, Vec<f64>) {
    let (mut h_w, mut h_z) = (Vec::new(), Vec::new());
    let (mut cnt_z, mut cnt_w, mut queue) = (0i64, 0i64, 0i64);
    let n = combined.len();
    for k in (0..n).rev() {
        if combined[k].1 {
            cnt_w += 1;
        } else {
            cnt_z += 1;
            queue += 1;
        }
        if k == 0 || combined[k].0 != combined[k - 1].0 {
            for _ in 0..queue {
                h_w.push(cnt_w as f64);
                h_z.push(cnt_z as f64);
            }
            queue = 0;
        }
    }
    (h_w, h_z)
}

/// `PosteriorEstimator::getQValues` with `includeNegativesInResult = true`
/// (set by every `Scores::calcQvals`): one q-value per entry of `combined`,
/// which is sorted best first.
pub(super) fn q_values(
    pi0: f64,
    combined: &[(f64, bool)],
    skip_decoys_plus_one: bool,
    null_target_win_prob: f64,
) -> Result<Vec<f64>> {
    let decoy_factor = null_target_win_prob / (1.0 - null_target_win_prob);
    let (h_w, h_z) = if pi0 < 1.0 {
        mix_max_counts(combined)
    } else {
        (Vec::new(), Vec::new())
    };
    let mut e_f1_mod_run_tot = 0.0;
    let mut n_z_ge_w: i64 = if skip_decoys_plus_one { 0 } else { 1 };
    let mut n_w_ge_w: i64 = 0;
    let (mut decoy_queue, mut target_queue) = (0i64, 0i64);
    let mut q = Vec::with_capacity(combined.len());
    let n = combined.len();
    for k in 0..n {
        if combined[k].1 {
            n_w_ge_w += 1;
            target_queue += 1;
        } else {
            n_z_ge_w += 1;
            decoy_queue += 1;
        }
        if k + 1 == n || combined[k].0 != combined[k + 1].0 {
            if pi0 < 1.0 && decoy_queue > 0 {
                let j = h_w.len() as i64 - (n_z_ge_w - 1);
                let cnt_w = at(&h_w, j)? as i64;
                let cnt_z = at(&h_z, j)? as i64;
                let mut est = (cnt_w as f64 - pi0 * cnt_z as f64) / ((1.0 - pi0) * cnt_z as f64);
                est = if est > 1.0 { 1.0 } else { est };
                est = if est < 0.0 { 0.0 } else { est };
                e_f1_mod_run_tot += decoy_queue as f64 * est * (1.0 - pi0);
            }
            target_queue += decoy_queue;
            let fdr = (n_z_ge_w as f64 * pi0 + e_f1_mod_run_tot) / (n_w_ge_w.max(1) as f64)
                * decoy_factor;
            for _ in 0..target_queue {
                q.push(std_min(fdr, 1.0));
            }
            decoy_queue = 0;
            target_queue = 0;
        }
    }
    reverse_running_min(&mut q);
    Ok(q)
}

/// `PosteriorEstimator::getPValues`: empirical p-values of the targets
/// against the decoys, ascending.
pub(super) fn p_values(combined: &[(f64, bool)]) -> Vec<f64> {
    let mut p = Vec::new();
    let (mut n_decoys, mut pos_same, mut neg_same) = (1usize, 0usize, 0usize);
    let n = combined.len();
    for k in 0..n {
        if combined[k].1 {
            pos_same += 1;
        } else {
            neg_same += 1;
        }
        if k + 1 == n || combined[k].0 != combined[k + 1].0 {
            for ix in 0..pos_same {
                p.push(
                    n_decoys as f64
                        + (neg_same.wrapping_mul(ix + 1)) as f64 / (pos_same + 1) as f64,
                );
            }
            n_decoys += neg_same;
            neg_same = 0;
            pos_same = 0;
        }
    }
    let total = n_decoys as f64;
    for value in &mut p {
        *value /= total;
    }
    p
}

/// `numLambda` and `maxLambda` of `PosteriorEstimator.cpp`.
const NUM_LAMBDA: u32 = 100;
const MAX_LAMBDA: f64 = 0.5;

/// `PosteriorEstimator::checkSeparation`: true when no p-value reaches
/// `maxLambda / numLambda`.
fn check_separation(p: &[f64]) -> bool {
    let min_lambda = (1.0 / f64::from(NUM_LAMBDA)) * MAX_LAMBDA;
    lower_bound(p, |v| *v < min_lambda) == p.len()
}

/// `PosteriorEstimator::estimatePi0` (Storey 2002) with 100 bootstrap
/// rounds drawn from `rng`.
fn estimate_pi0(p: &[f64], rng: &mut PseudoRandom) -> Result<f64> {
    let (mut lambdas, mut pi0s) = (Vec::new(), Vec::new());
    let n = p.len();
    for ix in 0..=NUM_LAMBDA {
        let lambda = (f64::from(ix + 1) / f64::from(NUM_LAMBDA)) * MAX_LAMBDA;
        let start = lower_bound(p, |v| *v < lambda);
        let wl = (n - start) as f64;
        let pi0 = wl / n as f64 / (1.0 - lambda);
        if pi0 > 0.0 {
            lambdas.push(lambda);
            pi0s.push(pi0);
        }
    }
    if pi0s.is_empty() {
        // "too good separation": with NO_TERMINATE and the static usePi0_
        // (true) the source returns 1.
        return Ok(1.0);
    }
    let mut min_pi0 = pi0s[0];
    for &value in &pi0s[1..] {
        if value < min_pi0 {
            min_pi0 = value;
        }
    }
    let mut mse = vec![0.0; pi0s.len()];
    let mut boot = Vec::new();
    for _ in 0..100 {
        bootstrap(p, &mut boot, rng)?;
        let nb = boot.len();
        for (ix, &lambda) in lambdas.iter().enumerate() {
            let start = lower_bound(&boot, |v| *v < lambda);
            let wl = (nb - start) as f64;
            let pi0_boot = wl / nb as f64 / (1.0 - lambda);
            mse[ix] += (pi0_boot - min_pi0) * (pi0_boot - min_pi0);
        }
    }
    let mut min_ix = 0;
    for ix in 1..mse.len() {
        if mse[ix] < mse[min_ix] {
            min_ix = ix;
        }
    }
    Ok(std_max(std_min(pi0s[min_ix], 1.0), 0.0))
}

/// The source's `bootstrap<double>`: at most 1000 draws with replacement,
/// sorted ascending.
fn bootstrap(input: &[f64], out: &mut Vec<f64>, rng: &mut PseudoRandom) -> Result<()> {
    out.clear();
    let n = input.len() as f64;
    let draws = input.len().min(1000);
    for _ in 0..draws {
        let draw = (rng.lcg_uniform_rand() * n) as usize;
        let value = *input
            .get(draw)
            .ok_or_else(|| Error::InvalidValue("bootstrap draw out of range".into()))?;
        out.push(value);
    }
    source_sort_by(out, |a, b| a < b)
}

/// One cross-validation training input (`AlgIn`): feature rows by PSM index,
/// labels, and the positive and negative counts.
#[derive(Clone, Debug, Default)]
pub(super) struct AlgIn {
    pub(super) vals: Vec<usize>,
    pub(super) y: Vec<f64>,
    pub(super) positives: usize,
    pub(super) negatives: usize,
}

/// `Scores`: a list of score holders with their target/decoy totals and the
/// pi0 and null-win-probability that go into their q-values.
#[derive(Clone, Debug)]
pub(super) struct Scores {
    use_pi0: bool,
    pi0: f64,
    target_decoy_size_ratio: f64,
    null_target_win_prob: f64,
    total_decoys: u32,
    total_targets: u32,
    pub(super) scores: Vec<ScoreHolder>,
}

impl Scores {
    pub(super) fn new(use_pi0: bool) -> Self {
        Self {
            use_pi0,
            pi0: 1.0,
            target_decoy_size_ratio: 1.0,
            null_target_win_prob: 0.5,
            total_decoys: 0,
            total_targets: 0,
            scores: Vec::new(),
        }
    }

    pub(super) fn pi0(&self) -> f64 {
        self.pi0
    }

    pub(super) fn target_decoy_size_ratio(&self) -> f64 {
        self.target_decoy_size_ratio
    }

    fn add(&mut self, holder: ScoreHolder) {
        self.scores.push(holder);
    }

    /// `populateWithPSMs`: the targets, then the decoys, in registration
    /// order, all scored 0.
    pub(super) fn populate(&mut self, targets: &[usize], decoys: &[usize]) {
        self.scores.clear();
        self.scores.extend(
            targets
                .iter()
                .map(|&psm| ScoreHolder::new(psm, Label::Target)),
        );
        self.scores.extend(
            decoys
                .iter()
                .map(|&psm| ScoreHolder::new(psm, Label::Decoy)),
        );
        self.total_targets = targets.len() as u32;
        self.total_decoys = decoys.len() as u32;
        self.target_decoy_size_ratio = f64::from(self.total_targets) / f64::from(self.total_decoys);
    }

    /// `scoreAndAddPSM`: `score = sum_j f[j] * w[j] + w[n]`, in feature order
    /// (unlike [`Scores::calc_scores`], which sums from the back), then append.
    pub(super) fn score_and_add(
        &mut self,
        mut holder: ScoreHolder,
        raw_weights: &[f64],
        psms: &[Psm],
    ) {
        let features = &psms[holder.psm].features;
        let n = features.len();
        for j in 0..n {
            holder.score += features[j] * raw_weights[j];
        }
        holder.score += raw_weights[n];
        if holder.is_target() {
            self.total_targets += 1;
        } else {
            self.total_decoys += 1;
        }
        self.scores.push(holder);
    }

    /// `merge`: per fold, sort best first, set pi0, compute q-values,
    /// optionally rescale the scores (and the fold's weights), then
    /// concatenate and [`Scores::post_merge_step`].
    pub(super) fn merge(
        &mut self,
        folds: &mut [Scores],
        fdr: f64,
        skip_normalize_scores: bool,
        weights: &mut [Vec<f64>],
        psms: &[Psm],
        rng: &mut PseudoRandom,
    ) -> Result<()> {
        self.scores.clear();
        for (fold, w) in folds.iter_mut().zip(weights.iter_mut()) {
            source_sort_by(&mut fold.scores, |a, b| greater(a, b, psms))?;
            fold.check_separation_and_set_pi0(rng)?;
            fold.calc_qvals(fdr, false)?;
            if !skip_normalize_scores {
                fold.normalize_scores(fdr, w);
            }
            self.scores.extend(fold.scores.iter().cloned());
        }
        self.post_merge_step(psms, rng)
    }

    /// `postMergeStep`: sort best first, recount, and set pi0.
    pub(super) fn post_merge_step(&mut self, psms: &[Psm], rng: &mut PseudoRandom) -> Result<()> {
        source_sort_by(&mut self.scores, |a, b| greater(a, b, psms))?;
        self.total_decoys = self.scores.iter().filter(|s| s.is_decoy()).count() as u32;
        self.total_targets = self.scores.iter().filter(|s| s.is_target()).count() as u32;
        self.target_decoy_size_ratio =
            f64::from(self.total_targets) / std_max(1.0, f64::from(self.total_decoys));
        self.check_separation_and_set_pi0(rng)
    }

    /// `calcScore`: the bias, then the features added from the last to the
    /// first.
    fn calc_score(features: &[f64], w: &[f64]) -> f64 {
        let mut ix = features.len();
        let mut score = w[ix];
        while ix > 0 {
            ix -= 1;
            score += features[ix] * w[ix];
        }
        score
    }

    /// `calcScores`: score every PSM with `w`, then sort best first.
    pub(super) fn calc_scores(&mut self, w: &[f64], psms: &[Psm]) -> Result<()> {
        for holder in &mut self.scores {
            holder.score = Self::calc_score(&psms[holder.psm].features, w);
        }
        source_sort_by(&mut self.scores, |a, b| greater(a, b, psms))
    }

    /// `calcScoresAndQvals`: [`Scores::calc_scores`] then
    /// [`Scores::calc_qvals`]; the number of targets below `fdr`.
    pub(super) fn calc_scores_and_qvals(
        &mut self,
        w: &[f64],
        fdr: f64,
        skip_decoys_plus_one: bool,
        psms: &[Psm],
    ) -> Result<i32> {
        self.calc_scores(w, psms)?;
        self.calc_qvals(fdr, skip_decoys_plus_one)
    }

    fn score_label_pairs(&self) -> Vec<(f64, bool)> {
        self.scores
            .iter()
            .map(|s| (s.score, s.is_target()))
            .collect()
    }

    /// `calcQvals`: q-values in the current order; the number of targets
    /// with `q < fdr`.
    pub(super) fn calc_qvals(&mut self, fdr: f64, skip_decoys_plus_one: bool) -> Result<i32> {
        let combined = self.score_label_pairs();
        let q = q_values(
            self.pi0,
            &combined,
            skip_decoys_plus_one,
            self.null_target_win_prob,
        )?;
        let mut positives = 0;
        for (holder, value) in self.scores.iter_mut().zip(q) {
            holder.q = value;
            if holder.q < fdr && holder.is_target() {
                positives += 1;
            }
        }
        Ok(positives)
    }

    /// `generateNegativeTrainingSet`: every decoy, label -1.
    pub(super) fn generate_negative_training_set(&self, data: &mut AlgIn) {
        data.vals.clear();
        data.y.clear();
        for holder in self.scores.iter().filter(|s| s.is_decoy()) {
            data.vals.push(holder.psm);
            data.y.push(-1.0);
        }
        data.negatives = data.vals.len();
    }

    /// `generatePositiveTrainingSet`: the targets with `q <= fdr`, label 1,
    /// after the negatives. With `train_best_positive` the holders are first
    /// sorted by spectrum and label and de-duplicated per spectrum and label
    /// **without erasing the tail**, and only the kept prefix is re-sorted
    /// and read, exactly as the source leaves this `Scores`.
    pub(super) fn generate_positive_training_set(
        &mut self,
        data: &mut AlgIn,
        fdr: f64,
        train_best_positive: bool,
        psms: &[Psm],
    ) -> Result<()> {
        data.vals.truncate(data.negatives);
        data.y.truncate(data.negatives);
        let mut end = self.scores.len();
        if train_best_positive {
            source_sort_by(&mut self.scores, |a, b| order_scan_label(a, b, psms))?;
            end = std_unique(&mut self.scores, |x, y| {
                let (px, py) = (&psms[x.psm], &psms[y.psm]);
                px.spec_file == py.spec_file && px.scan == py.scan && x.label == y.label
            });
            source_sort_by(&mut self.scores[..end], |a, b| greater(a, b, psms))?;
        }
        let mut positives = 0;
        for holder in &self.scores[..end] {
            if holder.is_target() && holder.q <= fdr {
                data.vals.push(holder.psm);
                data.y.push(1.0);
                positives += 1;
            }
        }
        data.positives = positives;
        Ok(())
    }

    /// `weedOutRedundantTDC`: keep the best-scoring holder per (spectrum
    /// file, scan, experimental mass), then [`Scores::post_merge_step`].
    pub(super) fn weed_out_redundant_tdc(
        &mut self,
        psms: &[Psm],
        rng: &mut PseudoRandom,
    ) -> Result<()> {
        source_sort_by(&mut self.scores, |a, b| order_scan_mass_charge(a, b, psms))?;
        let end = std_unique(&mut self.scores, |x, y| {
            let (px, py) = (&psms[x.psm], &psms[y.psm]);
            px.spec_file == py.spec_file && px.scan == py.scan && px.exp_mass == py.exp_mass
        });
        self.scores.truncate(end);
        self.post_merge_step(psms, rng)
    }

    /// `recalculateSizes`: anything not a target counts as a decoy.
    fn recalculate_sizes(&mut self) {
        self.total_targets = self.scores.iter().filter(|s| s.is_target()).count() as u32;
        self.total_decoys = self.scores.len() as u32 - self.total_targets;
        self.target_decoy_size_ratio = f64::from(self.total_targets) / f64::from(self.total_decoys);
    }

    /// `getInitDirection`: the single feature (and sign) that separates the
    /// most targets at `initial_selection_fdr`, with decoys+1 skipped.
    ///
    /// Only the first `n_features` entries of `direction` are reset; the
    /// bias entry keeps its previous value, as in the source.
    pub(super) fn init_direction(
        &mut self,
        initial_selection_fdr: f64,
        direction: &mut [f64],
        n_features: usize,
        psms: &[Psm],
    ) -> Result<i32> {
        let mut best_positives = -1;
        let mut best_feature: i64 = -1;
        let mut low_best = false;
        for feature in 0..n_features {
            for holder in &mut self.scores {
                holder.score = psms[holder.psm].features[feature];
            }
            source_sort_by(&mut self.scores, |a, b| less(a, b, psms))?;
            for i in 0..2 {
                if i == 1 {
                    self.scores.reverse();
                }
                let positives = self.calc_qvals(initial_selection_fdr, true)?;
                if positives > best_positives {
                    best_positives = positives;
                    best_feature = feature as i64;
                    low_best = i == 0;
                }
            }
        }
        for value in direction.iter_mut().take(n_features) {
            *value = 0.0;
        }
        if best_positives <= 0 {
            // "cannot find an initial direction": NO_TERMINATE falls back to
            // the first feature.
            best_feature = 0;
        }
        if best_feature >= 0 {
            let index = best_feature as usize;
            if let Some(slot) = direction.get_mut(index) {
                *slot = if low_best { -1.0 } else { 1.0 };
            }
        }
        Ok(best_positives)
    }

    /// `checkSeparationAndSetPi0`.
    pub(super) fn check_separation_and_set_pi0(&mut self, rng: &mut PseudoRandom) -> Result<()> {
        let combined = self.score_label_pairs();
        let p = p_values(&combined);
        self.pi0 = 1.0;
        if check_separation(&p) {
            // "too good separation": NO_TERMINATE keeps pi0 = 1.
        } else if self.use_pi0 {
            self.pi0 = estimate_pi0(&p, rng)?;
        }
        Ok(())
    }

    /// `createXvalSetsBySpectrum` with `decoyFractionTraining = 1` and one
    /// decoy per target, the only configuration the wrapper reaches.
    ///
    /// Sorts this set by the spectrum hash, then deals whole spectra to
    /// folds chosen by `rng`, never to a fold that is already full when a
    /// new spectrum starts.
    ///
    /// # Errors
    ///
    /// Native guard: [`Error::InvalidValue`] if ten million consecutive
    /// draws find no fold with room, which the source would loop on forever.
    pub(super) fn create_xval_sets_by_spectrum(
        &mut self,
        train: &mut Vec<Scores>,
        test: &mut Vec<Scores>,
        folds: usize,
        psms: &[Psm],
        rng: &mut PseudoRandom,
    ) -> Result<()> {
        let use_pi0 = self.use_pi0;
        train.resize_with(folds, || Scores::new(use_pi0));
        test.resize_with(folds, || Scores::new(use_pi0));
        let mut remain = vec![0i64; folds];
        let mut ix = self.scores.len() as i64;
        let mut fold = folds as i64;
        while fold > 0 {
            fold -= 1;
            remain[fold as usize] = ix / (fold + 1);
            ix -= remain[fold as usize];
        }
        source_sort_by(&mut self.scores, |a, b| {
            scan_hash(&psms[a.psm]) < scan_hash(&psms[b.psm])
        })?;
        if self.scores.is_empty() {
            // "no scored PSMs": NO_TERMINATE returns here.
            return Ok(());
        }
        let modulus = folds as u64;
        let mut previous = psms[self.scores[0].psm].scan;
        let mut rand_index = (rng.lcg_rand() % modulus) as usize;
        for holder in &self.scores {
            let current = psms[holder.psm].scan;
            if previous != current {
                rand_index = (rng.lcg_rand() % modulus) as usize;
                let mut draws = 0u32;
                while remain[rand_index] <= 0 {
                    draws += 1;
                    if draws > 10_000_000 {
                        return Err(Error::InvalidValue(
                            "cross-validation fold assignment does not terminate".into(),
                        ));
                    }
                    rand_index = (rng.lcg_rand() % modulus) as usize;
                }
            }
            for i in 0..folds {
                if i == rand_index {
                    test[i].add(holder.clone());
                }
                if i != rand_index || folds == 1 {
                    train[i].add(holder.clone());
                }
            }
            remain[rand_index] -= 1;
            previous = current;
        }
        // decoysPerTarget = 1, decoyFractionTraining = 1: factor 1, win
        // probability 1/2.
        let factor = 2.0 - 1.0 * 1.0;
        let null_target_win_prob = factor / (1.0 + 1.0);
        for i in 0..folds {
            train[i].recalculate_sizes();
            test[i].recalculate_sizes();
            train[i].null_target_win_prob = null_target_win_prob;
            test[i].null_target_win_prob = null_target_win_prob;
        }
        Ok(())
    }

    /// `normalizeScores`: rescale so the score at `fdr` maps to 0 and the
    /// median decoy to -1, and fold the same transform into `weights`
    /// (`Normalizer::endScoreNormalizeWeights`).
    ///
    /// When the median decoy does not score below the FDR threshold the
    /// source (with `NO_TERMINATE`) shifts by one instead of rescaling.
    pub(super) fn normalize_scores(&mut self, fdr: f64, weights: &mut [f64]) {
        let median_index = self.total_decoys / 2;
        if self.scores.is_empty() {
            return;
        }
        let mut fdr_score = self.scores[0].score;
        let mut median_decoy_score = fdr_score + 1.0;
        let mut decoys = 0u32;
        for holder in &self.scores {
            if holder.q < fdr {
                fdr_score = holder.score;
            }
            if holder.is_decoy() {
                decoys += 1;
                if decoys == median_index {
                    median_decoy_score = holder.score;
                    break;
                }
            }
        }
        let mut diff = fdr_score - median_decoy_score;
        if diff <= 0.0 {
            diff = 1.0;
            fdr_score += 1.0;
        }
        for holder in &mut self.scores {
            holder.score -= fdr_score;
            holder.score /= diff;
        }
        let n = weights.len() - 1;
        for w in weights.iter_mut().take(n) {
            *w /= diff;
        }
        weights[n] = (weights[n] - fdr_score) / diff;
    }

    /// `calcPep(spline = false, interp = false, pava)`: target PEPs from
    /// their q-values by monotone regression (PAVA or the I-spline fit),
    /// decoy PEPs by interpolation between the neighbouring targets.
    pub(super) fn calc_pep(&mut self, pava: bool) -> Result<()> {
        let target_q: Vec<f64> = self
            .scores
            .iter()
            .filter(|s| s.is_target())
            .map(|s| s.q)
            .collect();
        let regressor = InferPep::new(!pava);
        let mut target_pep = regressor.q_to_pep(&target_q)?;
        let mut target_q = target_q;
        target_pep.push(1.0);
        target_q.push(1.0);
        let (mut index, mut l_q, mut l_pep) = (0usize, 0.0, 0.0);
        for holder in &mut self.scores {
            if holder.is_target() {
                holder.pep = target_pep[index];
                l_pep = target_pep[index];
                l_q = target_q[index];
                index += 1;
            } else {
                let pep =
                    InferPep::interpolate(holder.q, l_q, target_q[index], l_pep, target_pep[index]);
                holder.pep = std_max(0.0, std_min(1.0, pep));
            }
        }
        Ok(())
    }
}

/// The OpenMS-local `sampleTrainingRowIndices` (`InProcessSampler.cpp`):
/// reservoir sampling of at most `max_train` rows, keyed by (scan,
/// experimental mass) so rows of one spectrum travel together, returned in
/// ascending row order. `max_train == 0` or `>= rows` selects every row.
pub(super) fn sample_training_rows(
    scans: &[i32],
    exp_masses: &[f64],
    max_train: usize,
    rng: &mut PseudoRandom,
) -> Vec<usize> {
    let n = scans.len();
    if max_train == 0 || max_train >= n {
        return (0..n).collect();
    }
    // std::map<std::pair<int, double>, size_t>: keys compare by scan, then
    // by mass under `<`; the masses are finite, so total_cmp agrees.
    let mut priorities: std::collections::BTreeMap<(i32, OrderedMass), u64> =
        std::collections::BTreeMap::new();
    let mut heap: std::collections::BinaryHeap<(u64, usize)> = std::collections::BinaryHeap::new();
    let mut upper_limit = u64::MAX;
    for i in 0..n {
        let key = (scans[i], OrderedMass(exp_masses[i]));
        let priority = match priorities.get(&key) {
            Some(&p) => p,
            None => {
                let p = rng.lcg_rand();
                priorities.insert(key, p);
                p
            }
        };
        if heap.len() < max_train || priority < upper_limit {
            heap.push((priority, i));
            if heap.len() > max_train {
                if let Some((top, _)) = heap.pop() {
                    upper_limit = top;
                }
            }
        }
    }
    let mut out: Vec<usize> = heap.into_iter().map(|(_, i)| i).collect();
    out.sort_unstable();
    out
}

/// A finite mass with the total order `std::pair`'s `<` gives it.
#[derive(Clone, Copy, Debug)]
struct OrderedMass(f64);
impl PartialEq for OrderedMass {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}
impl Eq for OrderedMass {}
impl PartialOrd for OrderedMass {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for OrderedMass {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // -0.0 and 0.0 are one key under `<`, as in the source map.
        if self.0 == other.0 {
            std::cmp::Ordering::Equal
        } else {
            self.0.total_cmp(&other.0)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lcg_matches_the_recurrence() {
        let mut rng = PseudoRandom::new(1);
        assert_eq!(rng.lcg_rand(), 279_470_273);
        assert_eq!(rng.lcg_rand(), (279_470_273u64 * 279_470_273) % K_RAND_MAX);
    }

    #[test]
    fn unique_keeps_the_libstdcxx_tail() {
        let mut v = vec![1, 1, 2, 3, 3, 4];
        let end = std_unique(&mut v, |a, b| a == b);
        assert_eq!(end, 4);
        assert_eq!(v, [1, 2, 3, 4, 3, 4]);
    }

    #[test]
    fn q_values_of_a_clean_split() {
        // targets 3, then 1 decoy: decoys+1 counting.
        let combined = [(4.0, true), (3.0, true), (2.0, true), (1.0, false)];
        let q = q_values(1.0, &combined, false, 0.5).unwrap();
        assert_eq!(q, [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0]);
    }

    #[test]
    fn p_values_count_decoys_above() {
        let combined = [(4.0, true), (3.0, false), (2.0, true), (1.0, false)];
        // nDecoys starts at 1 and ends at 3.
        assert_eq!(p_values(&combined), [1.0 / 3.0, 2.0 / 3.0]);
    }
}
