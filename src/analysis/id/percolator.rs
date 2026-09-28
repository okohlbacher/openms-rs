// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! In-process Percolator: semi-supervised target/decoy rescoring with
//! q-values and posterior error probabilities.
//!
//! Port of `ANALYSIS/ID/Percolator.h` and `Percolator.cpp` at `bc9cc12`,
//! together with the vendored Percolator subset it wraps
//! (`src/openms/thirdparty/percolator`, upstream commit `eb157f7`, Apache
//! License 2.0, with liblinear's BSD-3 TRON solver). No crate implements
//! Percolator, so the vendored algorithm is translated; its parts live in
//! the private submodules. `docs/PERCOLATOR_SUPPORT.md` maps every member,
//! carries the third-party notices, and lists the divergences.
//!
//! The public API has two groups of entry points.
//!
//! Combined training and scoring (original Percolator semantics, per-fold
//! score normalization on the merged set):
//! [`Percolator::rescore_peptide_ids`] for PSM-shaped data, which writes the
//! results back to each hit as meta values, and [`Percolator::rescore`],
//! domain-agnostic over a feature matrix, labels, optional cross-validation
//! group keys and feature names. The latter applies to PSM, transition or
//! peak-group rescoring, or any other semi-supervised target/decoy setting.
//!
//! Separate training and scoring (train on A, predict on B, with
//! fold-averaged weights): [`Percolator::train`] returns a
//! [`PercolatorModel`], [`Percolator::score`] applies one, and
//! [`save_model`] / [`load_model`] persist it as plain text.
//! `rescore` and `train` + `score` produce scores on the same scale but are
//! not bitwise equal on the same input; see [`Percolator::score`].
//!
//! # Preconditions
//!
//! At least ~100 decoys and enough discriminable targets are needed for
//! meaningful results. Features are continuous-valued; categorical features
//! must be one-hot encoded by the caller. Rows that share structure (several
//! transitions of one precursor) must share a cross-validation key, or
//! q-values will be optimistic.
//!
//! # Thread safety and reproducibility
//!
//! The source relies on process-wide statics (the random generator, the
//! feature count, the normalizer, the sanity-check direction) and is not
//! safe to call concurrently even on different instances. Here all of that
//! state belongs to one [`Percolator`] value, so instances are independent
//! and `Send`. Results are reproducible given the seed and the input order;
//! unlike the source, they do not depend on the thread count
//! ([`PercolatorParameters::num_threads`] only bounds parallelism).

mod engine;
mod pep;
mod svm;
mod training;

use self::engine::{Label, PseudoRandom, Psm, ScoreHolder, Scores, sample_training_rows};
use self::training::{CrossValidation, CvConfig, Normalizer, NormalizerKind, SanityCheck};
use super::percolator_types::{PercolatorModel, RescoreInput, RescoreOutput};
use crate::concept::parallel::Threads;
use crate::format::percolator_infile::{ScanNumberPattern, scan_identifier};
use crate::identification::{PeptideHit, PeptideIdentification};
use crate::math::source_sort::{TemporaryBuffer, source_stable_sort_permutation};
use crate::metadata::{MetaValue, MetaValueData};
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::{Error, Result};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

/// Parameter-handler name of the source class.
pub const PERCOLATOR_NAME: &str = "Percolator";

/// Native ceiling on `rows * features` of one call (2^30 values, 8 GiB of
/// `f64`). The source has none.
pub const MAX_INPUT_VALUES: usize = 1 << 30;

/// Native ceiling on the number of rows: the source numbers rows with `int`
/// (`static_cast<int>(i)` scan defaults).
pub const MAX_ROWS: usize = i32::MAX as usize;

/// Posterior error probability estimator (`pep_method`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PepMethod {
    /// `logistic_regression`, the default: despite the legacy name, PEPs
    /// are derived from the target q-values by a monotone I-spline fit
    /// (the vendored `InferPEP` with `ISPLINE_TRR`).
    #[default]
    LogisticRegression,
    /// `isotonic`: the same q-value route, fitted by PAVA.
    Isotonic,
    /// `nonparametric`: kernel smoothing (`PosteriorEstimator::estimatePEP`,
    /// a roughness-penalised logistic spline). **Not ported**: selecting it
    /// is refused with [`Error::Unsupported`].
    Nonparametric,
}

impl PepMethod {
    fn name(self) -> &'static str {
        match self {
            Self::LogisticRegression => "logistic_regression",
            Self::Isotonic => "isotonic",
            Self::Nonparametric => "nonparametric",
        }
    }
}

/// Feature normalizer used during training (`normalizer`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum FeatureNormalizer {
    /// `stdv`: zero mean, unit standard deviation (default).
    #[default]
    Stdv,
    /// `uni`: min/max scaling to `[0, 1]`.
    Uni,
    /// `none`: no normalization.
    None,
}

impl FeatureNormalizer {
    /// The source's name of this normalizer, as stored in a model file.
    pub fn name(self) -> &'static str {
        match self {
            Self::Stdv => "stdv",
            Self::Uni => "uni",
            Self::None => "none",
        }
    }
    fn kind(self) -> NormalizerKind {
        match self {
            Self::Stdv => NormalizerKind::Stdv,
            Self::Uni => NormalizerKind::Uni,
            Self::None => NormalizerKind::None,
        }
    }
}

/// Which Percolator output becomes a hit's main score after
/// [`Percolator::rescore_peptide_ids`] (`report_as_main_score`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MainScore {
    /// `none` (default): the main score and score type stay untouched.
    #[default]
    None,
    /// `q-value`: score type `q-value`, lower is better.
    QValue,
    /// `pep`: score type `Posterior Error Probability`, lower is better.
    Pep,
    /// `svm`: score type `svm`, higher is better.
    Svm,
}

impl MainScore {
    fn name(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::QValue => "q-value",
            Self::Pep => "pep",
            Self::Svm => "svm",
        }
    }
}

/// The parameters of [`Percolator`], one field per source parameter (see
/// [`Percolator::defaults`] for their descriptions and restrictions).
#[derive(Clone, Debug, PartialEq)]
pub struct PercolatorParameters {
    /// `c_pos`: positive regularization constant C+ (`>= 0`; 0 selects it by
    /// cross-validation from 10, 1 and 0.1).
    pub c_pos: f64,
    /// `c_neg`: negative regularization constant C- (`>= 0`; 0 selects the
    /// ratio by cross-validation).
    pub c_neg: f64,
    /// `test_fdr`: FDR threshold for SVM direction selection and score
    /// normalization, in `[0, 1]`.
    pub test_fdr: f64,
    /// `train_fdr`: training FDR for the SVM iterations, in `[0, 1]`.
    pub train_fdr: f64,
    /// `num_iterations`: SVM training iterations, `>= 1`.
    pub num_iterations: u32,
    /// `initial_direction`: feature name to start from, `-` prefix when
    /// lower is better; empty for automatic selection.
    pub initial_direction: String,
    /// `pep_method`.
    pub pep_method: PepMethod,
    /// `seed`: random seed for cross-validation splitting.
    pub seed: i32,
    /// `target_decoy_metavalue`: hit meta value holding
    /// `target`/`decoy` for [`Percolator::rescore_peptide_ids`].
    pub target_decoy_metavalue: String,
    /// `normalizer`.
    pub normalizer: FeatureNormalizer,
    /// `train_best_positive`: train on the best PSM per spectrum only.
    pub train_best_positive: bool,
    /// `post_processing_tdc`: target-decoy competition before q-values and
    /// PEPs.
    pub post_processing_tdc: bool,
    /// `nested_xval_bins`: nested cross-validation bins for the C grid
    /// search, `>= 1` (1 disables it).
    pub nested_xval_bins: u32,
    /// `subset_max_train`: cap on the training set size (0 = all rows).
    pub subset_max_train: u32,
    /// `report_as_main_score`.
    pub report_as_main_score: MainScore,
    /// `num_threads`: workers for cross-validation, in `1..=3`.
    pub num_threads: u32,
    /// `use_pi0`: pi0 correction of q-values and PEPs.
    pub use_pi0: bool,
}

impl Default for PercolatorParameters {
    /// The source defaults.
    fn default() -> Self {
        Self {
            c_pos: 0.1,
            c_neg: 0.0,
            test_fdr: 0.01,
            train_fdr: 0.01,
            num_iterations: 10,
            initial_direction: String::new(),
            pep_method: PepMethod::LogisticRegression,
            seed: 1,
            target_decoy_metavalue: "target_decoy".into(),
            normalizer: FeatureNormalizer::Stdv,
            train_best_positive: true,
            post_processing_tdc: true,
            nested_xval_bins: 1,
            subset_max_train: 100_000,
            report_as_main_score: MainScore::None,
            num_threads: 3,
            use_pi0: true,
        }
    }
}

/// The in-process Percolator: its [`PercolatorParameters`] and the two
/// results of the last call, read with [`Percolator::svm_weights`] and
/// [`Percolator::pi0`].
#[derive(Clone, Debug)]
pub struct Percolator {
    /// The parameters (`DefaultParamHandler` state of the source).
    pub parameters: PercolatorParameters,
    svm_weights: Vec<Vec<f64>>,
    last_pi0: f64,
    /// The generator state `score` continues from (the source's
    /// process-wide `PseudoRandom::seed_`, initially 1).
    rng: PseudoRandom,
}

impl Default for Percolator {
    fn default() -> Self {
        Self::with_parameters(PercolatorParameters::default())
    }
}

/// One input row's PIN-compatible keys (`RowKeys`).
struct RowKeys {
    scan: i32,
    spec_file: i32,
    exp_mass: f64,
}

/// `makeRowKeys`: explicit scan numbers, else the cross-validation key, else
/// the row index; spectrum file and mass default to 0.
fn row_keys(input: &RescoreInput, i: usize) -> RowKeys {
    let scan = if !input.scan_numbers.is_empty() {
        input.scan_numbers[i]
    } else if !input.cv_group_keys.is_empty() {
        input.cv_group_keys[i]
    } else {
        i as i32
    };
    RowKeys {
        scan,
        spec_file: input.spec_file_numbers.get(i).copied().unwrap_or(0),
        exp_mass: input.exp_masses.get(i).copied().unwrap_or(0.0),
    }
}

/// `makePSMFromRow`: the keys convert to `unsigned int` as in the source.
fn make_psm(input: &RescoreInput, i: usize) -> Psm {
    let keys = row_keys(input, i);
    Psm {
        features: input.features[i].clone(),
        exp_mass: keys.exp_mass,
        scan: keys.scan as u32,
        spec_file: keys.spec_file as u32,
        row: i,
    }
}

fn label_of(is_decoy: bool) -> Label {
    if is_decoy {
        Label::Decoy
    } else {
        Label::Target
    }
}

impl Percolator {
    /// A Percolator with the source defaults.
    pub fn new() -> Self {
        Self::default()
    }

    /// A Percolator with the given parameters.
    pub fn with_parameters(parameters: PercolatorParameters) -> Self {
        Self {
            parameters,
            svm_weights: Vec::new(),
            last_pi0: -1.0,
            rng: PseudoRandom::new(1),
        }
    }

    /// SVM weights of the last [`Percolator::rescore`] or
    /// [`Percolator::train`]: one fold-averaged weight vector in raw feature
    /// space with the bias last, wrapped in a one-element outer vector to
    /// keep the shape of older per-fold consumers. Empty before the first
    /// training; [`Percolator::score`] leaves it unchanged.
    pub fn svm_weights(&self) -> &[Vec<f64>] {
        &self.svm_weights
    }

    /// Pi0 (null fraction) of the last [`Percolator::rescore`] or
    /// [`Percolator::score`], estimated on the merged, post-TDC score set
    /// right before the PEP calculation; `1.0` when [`PercolatorParameters::use_pi0`]
    /// is off. `-1.0` before any such call, after [`Percolator::train`] (which
    /// does not score), and after a call that failed.
    pub fn pi0(&self) -> f64 {
        self.last_pi0
    }

    fn validate(&self) -> Result<()> {
        let bad = |message: &str| Err(Error::InvalidValue(message.into()));
        if self.parameters.c_pos < 0.0 || !self.parameters.c_pos.is_finite() {
            return bad("c_pos must be a finite value >= 0");
        }
        if self.parameters.c_neg < 0.0 || !self.parameters.c_neg.is_finite() {
            return bad("c_neg must be a finite value >= 0");
        }
        if !(0.0..=1.0).contains(&self.parameters.test_fdr)
            || !(0.0..=1.0).contains(&self.parameters.train_fdr)
        {
            return bad("test_fdr and train_fdr must lie in [0, 1]");
        }
        if self.parameters.num_iterations < 1 || self.parameters.nested_xval_bins < 1 {
            return bad("num_iterations and nested_xval_bins must be at least 1");
        }
        if !(1..=3).contains(&self.parameters.num_threads) {
            return bad("num_threads must lie in 1..=3");
        }
        if self.parameters.pep_method == PepMethod::Nonparametric {
            return Err(Error::Unsupported(
                "pep_method 'nonparametric' (PosteriorEstimator::estimatePEP) is not ported".into(),
            ));
        }
        Ok(())
    }

    /// `validateRescoreInput_` plus the native checks of
    /// [`RescoreInput::validate`] and the resource ceilings.
    fn check_input(input: &RescoreInput) -> Result<()> {
        if input.features.is_empty() {
            return Err(Error::InvalidValue(
                "RescoreInput has empty features matrix".into(),
            ));
        }
        if input.features.len() > MAX_ROWS
            || input
                .features
                .len()
                .checked_mul(input.n_features())
                .is_none_or(|cells| cells > MAX_INPUT_VALUES)
        {
            return Err(Error::InvalidValue(
                "RescoreInput exceeds the Percolator input ceiling".into(),
            ));
        }
        if input.n_features() == 0 {
            return Err(Error::InvalidValue(
                "RescoreInput has no feature columns".into(),
            ));
        }
        input.validate()
    }

    fn cv_config(&self) -> CvConfig {
        CvConfig {
            test_fdr: self.parameters.test_fdr,
            train_fdr: self.parameters.train_fdr,
            c_pos: self.parameters.c_pos,
            c_neg: self.parameters.c_neg,
            niter: self.parameters.num_iterations,
            use_pi0: self.parameters.use_pi0,
            nested_xval_bins: self.parameters.nested_xval_bins as usize,
            train_best_positive: self.parameters.train_best_positive,
            threads: Threads::from_cli(i64::from(self.parameters.num_threads)),
        }
    }

    /// Feature names as the source registers them: the caller's, else
    /// `feat_<j>`.
    fn registered_names(input: &RescoreInput) -> Vec<String> {
        (0..input.n_features())
            .map(|j| {
                input
                    .feature_names
                    .get(j)
                    .cloned()
                    .unwrap_or_else(|| format!("feat_{j}"))
            })
            .collect()
    }

    /// Normalize the training PSMs and run cross-validated training;
    /// returns the fold-merged score set and the averaged raw weights.
    fn run_training(
        &self,
        input: &RescoreInput,
        psms: &mut [Psm],
        rng: &mut PseudoRandom,
    ) -> Result<(Scores, Vec<f64>)> {
        let n_features = input.n_features();
        let names = Self::registered_names(input);
        let mut sanity = SanityCheck::new(&self.parameters.initial_direction, &names)?;
        let (targets, decoys): (Vec<usize>, Vec<usize>) =
            (0..psms.len()).partition(|&i| !input.is_decoy[psms[i].row]);
        let normalizer = {
            let rows: Vec<&[f64]> = targets
                .iter()
                .chain(&decoys)
                .map(|&i| psms[i].features.as_slice())
                .collect();
            Normalizer::fit(self.parameters.normalizer.kind(), &rows, n_features)
        };
        for &i in targets.iter().chain(&decoys) {
            normalizer.normalize(&mut psms[i].features);
        }
        let psms: &[Psm] = psms;
        let mut scores = Scores::new(self.parameters.use_pi0);
        scores.populate(&targets, &decoys);
        let mut cv = CrossValidation::new(&self.cv_config(), n_features);
        let trained: Result<Vec<f64>> = (|| {
            cv.pre_iteration_setup(&mut scores, &mut sanity, psms, rng)?;
            cv.train(psms, rng)?;
            cv.post_iteration_processing(&mut scores, &sanity, psms, rng)?;
            Ok(cv.average_weights(&normalizer))
        })();
        match trained {
            Ok(weights) => Ok((scores, weights)),
            Err(error) => Err(Error::InvalidValue(format!(
                "Percolator training failed: {error}"
            ))),
        }
    }

    /// The scoring tail shared by rescore and score: optional TDC,
    /// q-values, pi0, PEPs, and extraction by input row.
    fn finish_scoring(
        &mut self,
        mut scores: Scores,
        psms: &[Psm],
        n_rows: usize,
        rng: &mut PseudoRandom,
    ) -> Result<RescoreOutput> {
        let result = (|| {
            if self.parameters.post_processing_tdc {
                scores.weed_out_redundant_tdc(psms, rng)?;
            }
            scores.calc_qvals(self.parameters.test_fdr, false)?;
            let pi0 = scores.pi0();
            scores.calc_pep(self.parameters.pep_method == PepMethod::Isotonic)?;
            Ok(pi0)
        })();
        let pi0 = result.map_err(|error: Error| {
            Error::InvalidValue(format!("Percolator scoring failed: {error}"))
        })?;
        self.last_pi0 = pi0;
        Ok(extract(&scores, psms, n_rows))
    }

    /// Rescore a feature matrix domain-agnostically.
    ///
    /// Returns per-row SVM scores, q-values and PEPs aligned 1:1 with
    /// `input.features` (no reordering). The rows are not interpreted: for
    /// transition rescoring the q-values and PEPs are transition-level, and
    /// no roll-up to unique entities is done.
    ///
    /// When [`PercolatorParameters::subset_max_train`] is below the row count, a
    /// reservoir sample keyed by (scan, experimental mass) is trained on and
    /// the trained model is then applied to every row, as upstream
    /// `Caller.cpp` does.
    ///
    /// Target-decoy competition ([`PercolatorParameters::post_processing_tdc`], on by
    /// default) keeps one row per (spectrum file, scan, experimental mass);
    /// rows that lose are returned with `score = 0`, `q_value = 1`,
    /// `pep = 1`. To tell them from a row that scored zero, disable TDC or
    /// populate the PIN-compatible fields so the key is informative.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a malformed input (no rows, ragged rows,
    /// mismatched label or optional-field lengths, non-finite values, no
    /// feature column, above [`MAX_INPUT_VALUES`]), an unresolvable
    /// [`PercolatorParameters::initial_direction`], invalid parameters, or a failed
    /// training or scoring step; [`Error::Unsupported`] for
    /// [`PepMethod::Nonparametric`]. The source's own sanity failures (too
    /// few decoys, no discriminative feature) are warnings there, because
    /// the wrapper runs the vendored code with `NO_TERMINATE`; they are
    /// warnings (silently) here too.
    pub fn rescore(&mut self, input: &RescoreInput) -> Result<RescoreOutput> {
        self.last_pi0 = -1.0;
        Self::check_input(input)?;
        self.validate()?;
        let n_rows = input.n_rows();
        let mut rng = PseudoRandom::new(self.parameters.seed);
        let (scans, masses): (Vec<i32>, Vec<f64>) = (0..n_rows)
            .map(|i| {
                let keys = row_keys(input, i);
                (keys.scan, keys.exp_mass)
            })
            .unzip();
        let training_rows = sample_training_rows(
            &scans,
            &masses,
            self.parameters.subset_max_train as usize,
            &mut rng,
        );
        let did_subsample = training_rows.len() < n_rows;
        let mut psms: Vec<Psm> = training_rows.iter().map(|&i| make_psm(input, i)).collect();
        let (scores, mut weights) = self.run_training(input, &mut psms, &mut rng)?;
        self.svm_weights = if weights.is_empty() {
            Vec::new()
        } else {
            vec![weights.clone()]
        };
        if did_subsample {
            let all: Vec<Psm> = (0..n_rows).map(|i| make_psm(input, i)).collect();
            let prepared = (|| {
                let mut scores = Scores::new(self.parameters.use_pi0);
                for (index, psm) in all.iter().enumerate() {
                    scores.score_and_add(
                        ScoreHolder::new(index, label_of(input.is_decoy[psm.row])),
                        &weights,
                        &all,
                    );
                }
                scores.post_merge_step(&all, &mut rng)?;
                scores.calc_qvals(self.parameters.test_fdr, false)?;
                scores.normalize_scores(self.parameters.test_fdr, &mut weights);
                Ok(scores)
            })();
            let scores = prepared.map_err(|error: Error| {
                Error::InvalidValue(format!("Percolator scoring failed: {error}"))
            })?;
            let output = self.finish_scoring(scores, &all, n_rows, &mut rng);
            self.rng = rng;
            return output;
        }
        let output = self.finish_scoring(scores, &psms, n_rows, &mut rng);
        self.rng = rng;
        output
    }

    /// Train a Percolator model on feature rows; no scoring.
    ///
    /// Runs the full semi-supervised cross-validation training on every row
    /// (subsampling, if wanted, is the caller's). The returned model holds
    /// fold-averaged SVM weights in raw feature space and can be applied to
    /// a different input with [`Percolator::score`]. The weights are also
    /// mirrored into [`Percolator::svm_weights`]; [`Percolator::pi0`] becomes
    /// `-1`. Missing feature names are recorded as `feat_<j>`.
    ///
    /// # Errors
    ///
    /// As [`Percolator::rescore`], for malformed input or a failed training.
    pub fn train(&mut self, input: &RescoreInput) -> Result<PercolatorModel> {
        self.last_pi0 = -1.0;
        Self::check_input(input)?;
        self.validate()?;
        let mut rng = PseudoRandom::new(self.parameters.seed);
        let mut psms: Vec<Psm> = (0..input.n_rows()).map(|i| make_psm(input, i)).collect();
        let (_, weights) = self.run_training(input, &mut psms, &mut rng)?;
        self.rng = rng;
        self.svm_weights = if weights.is_empty() {
            Vec::new()
        } else {
            vec![weights.clone()]
        };
        Ok(PercolatorModel {
            format_version: 1,
            feature_names: Self::registered_names(input),
            weights,
            normalizer_type: self.parameters.normalizer.name().into(),
            seed: self.parameters.seed,
        })
    }

    /// Score feature rows with a trained model; no training.
    ///
    /// Applies `model.weights` to the raw features: `raw(i) = sum_j
    /// features[i][j] * weights[j] + weights[n]`. Callers must not
    /// pre-normalize: the normalization is folded into the weights. The
    /// post-scoring pipeline (q-values, optional TDC, PEPs) runs on this
    /// input's own target/decoy distribution, and the scores are then
    /// rescaled so the median decoy maps to -1 and the score at
    /// [`PercolatorParameters::test_fdr`] to 0:
    /// `score(i) = (raw(i) - fdrScore) / (fdrScore - medianDecoyScore)`.
    /// Rescoring rescales per fold instead, which is why `rescore(X)` and
    /// `score(X, train(X))` share a scale but are not bitwise equal; for raw
    /// dot products use the model weights directly. TDC losers get the
    /// defaults as in [`Percolator::rescore`].
    ///
    /// The source does not reseed its process-wide generator here, so the
    /// pi0 bootstrap continues from whatever state earlier calls left.
    /// This port continues from the state this instance's previous call
    /// left (seed 1 for a fresh instance), which equals the source for one
    /// instance used alone.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for malformed input, a model whose weight
    /// count is not `n_features + 1`, empty feature names on either side,
    /// feature names that differ from the input's, or a failed scoring step.
    pub fn score(
        &mut self,
        input: &RescoreInput,
        model: &PercolatorModel,
    ) -> Result<RescoreOutput> {
        self.last_pi0 = -1.0;
        Self::check_input(input)?;
        let n_features = input.n_features();
        if model.weights.len() != n_features + 1 {
            return Err(Error::InvalidValue(format!(
                "PercolatorModel weight count does not match input feature count (expected n_features + 1 for bias): {} vs {}",
                model.weights.len(),
                n_features + 1
            )));
        }
        if model.feature_names.is_empty() {
            return Err(Error::InvalidValue(
                "PercolatorModel has no feature_names — cannot validate input alignment".into(),
            ));
        }
        if input.feature_names.is_empty() {
            return Err(Error::InvalidValue(
                "RescoreInput.feature_names required when scoring against a trained model".into(),
            ));
        }
        if model.feature_names.len() != input.feature_names.len() {
            return Err(Error::InvalidValue(format!(
                "PercolatorModel feature_names size mismatch with input: {} vs {}",
                model.feature_names.len(),
                input.feature_names.len()
            )));
        }
        if let Some(column) = model
            .feature_names
            .iter()
            .zip(&input.feature_names)
            .position(|(a, b)| a != b)
        {
            return Err(Error::InvalidValue(format!(
                "PercolatorModel feature_names differ from input at column {column}: {} != {}",
                model.feature_names[column], input.feature_names[column]
            )));
        }
        self.validate()?;
        let n_rows = input.n_rows();
        let psms: Vec<Psm> = (0..n_rows).map(|i| make_psm(input, i)).collect();
        let mut rng = self.rng.clone();
        let prepared = (|| {
            let mut scores = Scores::new(self.parameters.use_pi0);
            for i in 0..n_rows {
                scores.score_and_add(
                    ScoreHolder::new(i, label_of(input.is_decoy[i])),
                    &model.weights,
                    &psms,
                );
            }
            scores.post_merge_step(&psms, &mut rng)?;
            scores.calc_qvals(self.parameters.test_fdr, false)?;
            let mut weights = model.weights.clone();
            scores.normalize_scores(self.parameters.test_fdr, &mut weights);
            Ok(scores)
        })();
        let scores = prepared.map_err(|error: Error| {
            Error::InvalidValue(format!("Percolator scoring failed: {error}"))
        })?;
        let output = self.finish_scoring(scores, &psms, n_rows, &mut rng);
        self.rng = rng;
        output
    }

    /// Rescore PSMs in place, a PSM-shaped wrapper over
    /// [`Percolator::rescore`].
    ///
    /// Every hit becomes one row (identification `i`'s hits share
    /// cross-validation key `i`) and gets the meta values
    /// `percolator_score` (SVM discriminant), `percolator_q_value` and
    /// `percolator_pep`. `feature_names` name numeric hit meta values; when
    /// empty, they are discovered from the first hit: every integer or float
    /// meta value except the target/decoy key, `percolator_score`,
    /// `percolator_q_value`, `percolator_pep`, `q-value`, `q-value_score`,
    /// `score`, `RT`, `MZ`, `target_decoy` and `isDecoy`, sorted by name.
    /// A hit is a decoy when its [`PercolatorParameters::target_decoy_metavalue`] text
    /// is `decoy`.
    ///
    /// With [`PercolatorParameters::report_as_main_score`] other than `none`, the
    /// chosen output also replaces each hit's score, the identification's
    /// score type and direction are set, and its hits are re-sorted (stable,
    /// as the source's `PeptideIdentification::sort`).
    ///
    /// An empty `peptide_ids` is left alone.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the first identification has no hits and
    /// no feature names are given, no numeric feature is found, a hit lacks
    /// the target/decoy meta value or a feature, a feature is not an integer
    /// or float (the source's `ConversionError`), or [`Percolator::rescore`]
    /// fails. On error `peptide_ids` is unchanged.
    pub fn rescore_peptide_ids(
        &mut self,
        peptide_ids: &mut [PeptideIdentification],
        feature_names: &[String],
    ) -> Result<()> {
        if peptide_ids.is_empty() {
            return Ok(());
        }
        let td_meta = self.parameters.target_decoy_metavalue.clone();
        let features = if feature_names.is_empty() {
            let first = peptide_ids[0].hits.first().ok_or_else(|| {
                Error::InvalidValue("first PeptideIdentification has no hits".into())
            })?;
            let discovered = discover_feature_names(first, &td_meta);
            if discovered.is_empty() {
                return Err(Error::InvalidValue(
                    "no numeric features found on first PeptideHit".into(),
                ));
            }
            discovered
        } else {
            feature_names.to_vec()
        };
        let mut input = RescoreInput {
            feature_names: features.clone(),
            ..RescoreInput::default()
        };
        let mut locations = Vec::new();
        for (i, identification) in peptide_ids.iter().enumerate() {
            for (j, hit) in identification.hits.iter().enumerate() {
                let label = hit.metadata.get(&td_meta).ok_or_else(|| {
                    Error::InvalidValue(format!(
                        "PeptideHit missing target/decoy meta value: {td_meta}"
                    ))
                })?;
                let mut row = Vec::with_capacity(features.len());
                for name in &features {
                    let value = hit.metadata.get(name).ok_or_else(|| {
                        Error::InvalidValue(format!(
                            "PeptideHit missing feature meta value: {name}"
                        ))
                    })?;
                    row.push(match value.data() {
                        MetaValueData::Integer(v) => *v as f64,
                        MetaValueData::Float(v) => *v,
                        _ => {
                            return Err(Error::InvalidValue(format!(
                                "feature meta value '{name}' is not numeric"
                            )));
                        }
                    });
                }
                input.features.push(row);
                input.is_decoy.push(label.to_string() == "decoy");
                input.cv_group_keys.push(
                    i32::try_from(i).map_err(|_| {
                        Error::InvalidValue("too many peptide identifications".into())
                    })?,
                );
                locations.push((i, j));
            }
        }
        let output = self.rescore(&input)?;
        let mut updated = peptide_ids.to_vec();
        for (row, &(i, j)) in locations.iter().enumerate() {
            let hit = &mut updated[i].hits[j];
            hit.metadata.insert(
                "percolator_score".into(),
                MetaValue::source_float(output.scores[row]),
            );
            hit.metadata.insert(
                "percolator_q_value".into(),
                MetaValue::source_float(output.q_values[row]),
            );
            hit.metadata.insert(
                "percolator_pep".into(),
                MetaValue::source_float(output.peps[row]),
            );
        }
        if self.parameters.report_as_main_score != MainScore::None {
            let mut touched = vec![false; updated.len()];
            for (row, &(i, j)) in locations.iter().enumerate() {
                let identification = &mut updated[i];
                identification.hits[j].score = match self.parameters.report_as_main_score {
                    MainScore::QValue => output.q_values[row],
                    MainScore::Pep => output.peps[row],
                    MainScore::Svm | MainScore::None => output.scores[row],
                };
                if !touched[i] {
                    touched[i] = true;
                    let (score_type, higher) = match self.parameters.report_as_main_score {
                        MainScore::Svm | MainScore::None => ("svm", true),
                        MainScore::Pep => ("Posterior Error Probability", false),
                        MainScore::QValue => ("q-value", false),
                    };
                    identification.score_type = score_type.into();
                    identification.higher_score_better = higher;
                }
            }
            for (identification, touched) in updated.iter_mut().zip(touched) {
                if touched {
                    stable_sort_hits(identification)?;
                }
            }
        }
        peptide_ids.clone_from_slice(&updated);
        Ok(())
    }

    /// Fill the four PIN-compatible optional fields of `input` from
    /// peptide identifications, with the conventions of the `.pin` writer
    /// ([`crate::format::percolator_infile`]).
    ///
    /// - scan: from the scan identifier (spectrum reference, `spectrum_id`,
    ///   or the one-based identification index) with the pattern the first
    ///   identification's identifier selects; `-1` when none is found.
    /// - spectrum file: a dense index over distinct `file_origin` +
    ///   `|` + `id_merge_index` texts, in first-seen order.
    /// - experimental mass: the precursor m/z (not converted to a neutral
    ///   mass), 0 without one.
    /// - calculated mass: the hit's numeric `CalcMass`, else the sequence
    ///   m/z at the hit's charge, 0 when that fails.
    ///
    /// Call this before [`Percolator::rescore`] when the in-process output
    /// must match the external percolator on the same `.pin` data: without
    /// it the row index stands in for the scan number, which changes the
    /// cross-validation split. `peptide_ids` must parallel `input.features`:
    /// one row per hit (`flatten_hits`) or per identification (its first
    /// hit); identifications without hits contribute no row. An empty
    /// `peptide_ids` leaves `input` untouched.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a `CalcMass` meta value is not an
    /// integer or float (the source's `ConversionError`); `input` is then
    /// unchanged.
    pub fn fill_pin_compatible_fields(
        peptide_ids: &[PeptideIdentification],
        flatten_hits: bool,
        input: &mut RescoreInput,
    ) -> Result<()> {
        if peptide_ids.is_empty() {
            return Ok(());
        }
        let pattern = ScanNumberPattern::for_native_id(&scan_identifier(&peptide_ids[0], 0));
        let mut spec_files: HashMap<String, i32> = HashMap::new();
        let (mut scans, mut files, mut exp, mut calc) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for (index, identification) in peptide_ids.iter().enumerate() {
            if identification.hits.is_empty() {
                continue;
            }
            let identifier = scan_identifier(identification, index + 1);
            let scan = pattern.extract(&identifier).unwrap_or(-1);
            let text = |key: &str| {
                identification
                    .metadata
                    .get(key)
                    .map(ToString::to_string)
                    .unwrap_or_default()
            };
            let key = format!("{}|{}", text("file_origin"), text("id_merge_index"));
            let next = spec_files.len() as i32;
            let spec_file = *spec_files.entry(key).or_insert(next);
            let exp_mass = identification.mz.unwrap_or(0.0);
            let count = if flatten_hits {
                identification.hits.len()
            } else {
                1
            };
            for hit in &identification.hits[..count] {
                let calc_mass = match hit.metadata.get("CalcMass") {
                    Some(value) => match value.data() {
                        MetaValueData::Integer(v) => *v as f64,
                        MetaValueData::Float(v) => *v,
                        _ => {
                            return Err(Error::InvalidValue(
                                "CalcMass meta value is not numeric".into(),
                            ));
                        }
                    },
                    None => hit.sequence.mz(hit.charge).unwrap_or(0.0),
                };
                scans.push(scan);
                files.push(spec_file);
                exp.push(exp_mass);
                calc.push(calc_mass);
            }
        }
        input.scan_numbers = scans;
        input.spec_file_numbers = files;
        input.exp_masses = exp;
        input.calc_masses = calc;
        Ok(())
    }

    /// The source defaults as a [`Param`], in declaration order.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure, which the constant tree
    /// cannot trigger.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The current parameters in the layout of [`Percolator::defaults`].
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a count does not fit the source's 32-bit
    /// integer parameter.
    pub fn to_param(&self) -> Result<Param> {
        let int = |name: &str, value: u32| -> Result<i64> {
            i32::try_from(value).map(i64::from).map_err(|_| {
                Error::InvalidValue(format!("'{name}' does not fit the source 32-bit parameter"))
            })
        };
        let flag = |value: bool| if value { "true" } else { "false" };
        let entries = vec![
            float_entry(
                "c_pos",
                self.parameters.c_pos,
                "Positive regularization constant (C+).",
                Some(0.0),
                None,
            ),
            float_entry(
                "c_neg",
                self.parameters.c_neg,
                "Negative regularization constant (C-). 0 = auto-select.",
                Some(0.0),
                None,
            ),
            float_entry(
                "test_fdr",
                self.parameters.test_fdr,
                "Test FDR threshold for SVM direction selection.",
                Some(0.0),
                Some(1.0),
            ),
            float_entry(
                "train_fdr",
                self.parameters.train_fdr,
                "Training FDR for SVM iterations.",
                Some(0.0),
                Some(1.0),
            ),
            int_entry(
                "num_iterations",
                int("num_iterations", self.parameters.num_iterations)?,
                "Number of SVM training iterations.",
                Some(1),
                None,
            ),
            string_entry(
                "initial_direction",
                &self.parameters.initial_direction,
                "Feature name to use as the initial scoring direction. Prefix with '-' if LOWER values are more target-like (e.g. 'COMET:lnExpect' means higher is better; '-COMET:lnExpect' means lower is better). Leave empty (default) to let Percolator auto-select by testing every feature in both directions and picking the one that separates the most targets at train_fdr.",
                &[],
            ),
            string_entry(
                "pep_method",
                self.parameters.pep_method.name(),
                "Posterior error probability estimator: 'logistic_regression' (target-decoy competition + isotonic fit; legacy name, not actually logistic regression), 'isotonic' (target-q → pep monotone regression), or 'nonparametric' (kernel smoothing).",
                &["logistic_regression", "isotonic", "nonparametric"],
            ),
            int_entry(
                "seed",
                i64::from(self.parameters.seed),
                "Random seed for CV splitting.",
                None,
                None,
            ),
            string_entry(
                "target_decoy_metavalue",
                &self.parameters.target_decoy_metavalue,
                "Meta-value on PeptideHit indicating target/decoy ('target'/'decoy').",
                &[],
            ),
            string_entry(
                "normalizer",
                self.parameters.normalizer.name(),
                "Feature normalizer: 'stdv' (zero-mean, unit-stdv, default), 'uni' (min/max), or 'none' (no normalization).",
                &["stdv", "uni", "none"],
            ),
            string_entry(
                "train_best_positive",
                flag(self.parameters.train_best_positive),
                "During TRAINING, use only the best PSM per spectrum as a positive example. Affects the SVM (weights change). Orthogonal to post_processing_tdc. Default is 'true' because typical OpenMS inputs are multi-hit-per-spectrum (Comet/MSGF+/Sage top-N, OSW peak groups). Inert on one-hit-per-spectrum inputs. Note: the upstream percolator binary defaults this to off — PercolatorAdapter preserves the CLI contract by mirroring its flag state when writing this param.",
                &["true", "false"],
            ),
            string_entry(
                "post_processing_tdc",
                flag(self.parameters.post_processing_tdc),
                "AFTER training, before PEP/q-value calculation, apply Target-Decoy Competition: keep only the highest-scoring PSM per (scan, mass, charge) and drop the rest from the score set. Does NOT touch the SVM — only the score distribution used for FDR. Default is 'true' because without TDC on multi-hit input, q-values count every candidate as an independent 'chance' and drift optimistic. Inert on one-hit-per-spectrum inputs. Note: the upstream percolator binary defaults this to off — PercolatorAdapter preserves the CLI contract by mirroring its flag state when writing this param.",
                &["true", "false"],
            ),
            int_entry(
                "nested_xval_bins",
                int("nested_xval_bins", self.parameters.nested_xval_bins)?,
                "Number of nested CV bins for grid search of Cpos/Cneg. Upstream default for cross validation is 3; 1 disables nested search.",
                Some(1),
                None,
            ),
            int_entry(
                "subset_max_train",
                int("subset_max_train", self.parameters.subset_max_train)?,
                "Cap on the training set size (0 = use all). Default 100000 bounds training time on very large inputs (>1M PSMs) at negligible quality loss. Sampling is reservoir-based, keyed by (scan, exp_mass) so all PSMs sharing a spectrum travel together; the trained model is then applied to ALL input rows via a second scoring pass (matches upstream Caller.cpp behavior). Inert when the input size is at or below the cap. Note: upstream percolator binary defaults this to 0 — PercolatorAdapter preserves the CLI contract by passing its own int option (CLI default 0) through to this param.",
                Some(0),
                None,
            ),
            string_entry(
                "report_as_main_score",
                self.parameters.report_as_main_score.name(),
                "Which Percolator output to stamp as the main PeptideHit score after rescoring. 'none' (default) leaves hit.getScore() and the PeptideIdentification score_type untouched — callers who want the original search score preserved (or who run their own post-processing like PercolatorAdapter) rely on this. 'q-value' sets main score to the PSM-level q-value and marks higher_score_better=false. 'pep' sets main score to the PEP and marks higher_score_better=false. 'svm' sets main score to the SVM discriminant and marks higher_score_better=true. All three Percolator outputs remain accessible as percolator_{score,q_value,pep} meta values regardless of this setting.",
                &["none", "q-value", "pep", "svm"],
            ),
            int_entry(
                "num_threads",
                int("num_threads", self.parameters.num_threads)?,
                "Number of OpenMP threads for cross-validation (one per CV fold; max useful = 3). Matches upstream's default.",
                Some(1),
                Some(3),
            ),
            string_entry(
                "use_pi0",
                flag(self.parameters.use_pi0),
                "Enable pi0 correction when computing q-values and PEPs (default). When true, pi0 is estimated per CV fold via bootstrap spline fit on the target score distribution using decoys as the null. When false, pi0 is fixed at 1.0 (pure target-decoy FDR, no bootstrap estimation). Affects both Scores::calcQvals and Scores::calcPep. Upstream percolator has no CLI flag for this; the in-library Scores(bool usePi0) constructor is the only toggle, matched here.",
                &["true", "false"],
            ),
        ];
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries,
            nodes: Vec::new(),
        })
    }

    /// Build a Percolator from parameters, as source `setParameters`
    /// followed by `updateMembers_`. Missing parameters take their defaults;
    /// unknown ones are ignored (see
    /// [`Percolator::from_param_with_warnings`]).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a parameter has the wrong type or
    /// violates its restriction (the source `Exception::InvalidParameter`).
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`Percolator::from_param`], also returning the unknown-parameter
    /// warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`Percolator::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(PERCOLATOR_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, Self::from_complete_param)
    }

    fn from_complete_param(param: &Param) -> Result<Self> {
        let float = |key: &str| match param.value(key)? {
            ParamValue::Float(value) => Ok(*value),
            _ => Err(Error::InvalidValue(format!(
                "parameter '{key}' must be a floating-point value"
            ))),
        };
        let int = |key: &str| param.value(key)?.to_i64();
        let text = |key: &str| Ok::<String, Error>(param.value(key)?.as_str()?.to_owned());
        let unsigned = |key: &str| {
            u32::try_from(int(key)?)
                .map_err(|_| Error::InvalidValue(format!("parameter '{key}' must not be negative")))
        };
        let seed = i32::try_from(int("seed")?)
            .map_err(|_| Error::InvalidValue("parameter 'seed' does not fit 32 bits".into()))?;
        Ok(Self::with_parameters(PercolatorParameters {
            c_pos: float("c_pos")?,
            c_neg: float("c_neg")?,
            test_fdr: float("test_fdr")?,
            train_fdr: float("train_fdr")?,
            num_iterations: unsigned("num_iterations")?,
            seed,
            pep_method: match text("pep_method")?.as_str() {
                "isotonic" => PepMethod::Isotonic,
                "nonparametric" => PepMethod::Nonparametric,
                _ => PepMethod::LogisticRegression,
            },
            normalizer: match text("normalizer")?.as_str() {
                "uni" => FeatureNormalizer::Uni,
                "none" => FeatureNormalizer::None,
                _ => FeatureNormalizer::Stdv,
            },
            train_best_positive: text("train_best_positive")? == "true",
            post_processing_tdc: text("post_processing_tdc")? == "true",
            nested_xval_bins: unsigned("nested_xval_bins")?,
            subset_max_train: unsigned("subset_max_train")?,
            initial_direction: text("initial_direction")?,
            use_pi0: text("use_pi0")? == "true",
            num_threads: unsigned("num_threads")?,
            report_as_main_score: match text("report_as_main_score")?.as_str() {
                "q-value" => MainScore::QValue,
                "pep" => MainScore::Pep,
                "svm" => MainScore::Svm,
                _ => MainScore::None,
            },
            target_decoy_metavalue: text("target_decoy_metavalue")?,
        }))
    }
}

fn entry(name: &str, value: ParamValue, description: &str) -> ParamEntry {
    ParamEntry {
        name: name.into(),
        description: description.into(),
        value,
        ..ParamEntry::default()
    }
}

fn float_entry(
    name: &str,
    value: f64,
    description: &str,
    min: Option<f64>,
    max: Option<f64>,
) -> ParamEntry {
    let mut e = entry(name, ParamValue::Float(value), description);
    if let Some(min) = min {
        e.min_float = min;
    }
    if let Some(max) = max {
        e.max_float = max;
    }
    e
}

fn int_entry(
    name: &str,
    value: i64,
    description: &str,
    min: Option<i32>,
    max: Option<i32>,
) -> ParamEntry {
    let mut e = entry(name, ParamValue::Integer(value), description);
    if let Some(min) = min {
        e.min_int = min;
    }
    if let Some(max) = max {
        e.max_int = max;
    }
    e
}

fn string_entry(name: &str, value: &str, description: &str, valid: &[&str]) -> ParamEntry {
    let mut e = entry(name, ParamValue::String(value.into()), description);
    e.valid_strings = valid.iter().map(|v| (*v).to_owned()).collect();
    e
}

/// `discoverFeatureNames_`: integer and float meta values of `hit`, minus
/// the target/decoy key and the blocklist, sorted.
fn discover_feature_names(hit: &PeptideHit, td_meta: &str) -> Vec<String> {
    const BLOCKLIST: [&str; 10] = [
        "percolator_score",
        "percolator_q_value",
        "percolator_pep",
        "q-value",
        "q-value_score",
        "score",
        "RT",
        "MZ",
        "target_decoy",
        "isDecoy",
    ];
    let mut out: Vec<String> = hit
        .metadata
        .iter()
        .filter(|(key, _)| key.as_str() != td_meta && !BLOCKLIST.contains(&key.as_str()))
        .filter(|(_, value)| {
            matches!(
                value.data(),
                MetaValueData::Integer(_) | MetaValueData::Float(_)
            )
        })
        .map(|(key, _)| key.clone())
        .collect();
    out.sort();
    out
}

/// `PeptideIdentification::sort`: `std::stable_sort` by score in the
/// identification's direction, without the Rust record's validation.
fn stable_sort_hits(identification: &mut PeptideIdentification) -> Result<()> {
    let higher = identification.higher_score_better;
    let hits = &identification.hits;
    let order = source_stable_sort_permutation(
        hits.len(),
        |a, b| {
            if higher {
                hits[a].score > hits[b].score
            } else {
                hits[a].score < hits[b].score
            }
        },
        TemporaryBuffer::Allocate,
    )?;
    let sorted: Vec<PeptideHit> = order
        .iter()
        .map(|&i| identification.hits[i].clone())
        .collect();
    identification.hits = sorted;
    Ok(())
}

/// The output extraction of `rescore`/`score`: defaults 0/1/1, then every
/// surviving score holder writes its row.
fn extract(scores: &Scores, psms: &[Psm], n_rows: usize) -> RescoreOutput {
    let mut out = RescoreOutput {
        scores: vec![0.0; n_rows],
        q_values: vec![1.0; n_rows],
        peps: vec![1.0; n_rows],
    };
    for holder in &scores.scores {
        let row = psms[holder.psm].row;
        if row < n_rows {
            out.scores[row] = holder.score;
            out.q_values[row] = holder.q;
            out.peps[row] = holder.pep;
        }
    }
    out
}

/// The model file's second comment line, with its em dash.
const FEATURES_COMMENT: &str = "# features \u{2014} one 'name<TAB>weight' per line";

/// Serialize a model as plain text (`Percolator::saveModel`).
///
/// A comment line, the header keys `format_version`, `normalizer`, `seed`,
/// `n_features` and `bias` as `key: value`, a second comment, then one
/// `name<TAB>weight` row per feature. The bias lives in the header, so
/// feature names carry no reserved values. Doubles are written as the
/// source's stream at `max_digits10` (17 significant digits, `%.17g`), so
/// [`load_model`] round-trips losslessly. The format is not the external
/// percolator's multi-column `.weights` format.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `weights.len() != feature_names.len() + 1`
/// (the source's `InvalidValue`), and, as native checks the source lacks,
/// when a feature name contains a tab or a line break or starts with `#`,
/// or the normalizer is not `stdv`, `uni` or `none`: the source writes such
/// a model, but its own reader then misreads or rejects it.
/// [`Error::Io`] when the file cannot be written (the source's
/// `UnableToCreateFile`).
pub fn save_model(model: &PercolatorModel, path: impl AsRef<Path>) -> Result<()> {
    let text = model_text(model)?;
    let mut file = std::fs::File::create(path).map_err(Error::Io)?;
    file.write_all(text.as_bytes()).map_err(Error::Io)
}

/// [`save_model`] to a writer.
///
/// # Errors
///
/// As [`save_model`].
pub fn write_model(model: &PercolatorModel, mut writer: impl Write) -> Result<()> {
    writer
        .write_all(model_text(model)?.as_bytes())
        .map_err(Error::Io)
}

fn model_text(model: &PercolatorModel) -> Result<String> {
    if model.weights.len() != model.feature_names.len() + 1 {
        return Err(Error::InvalidValue(format!(
            "PercolatorModel invariant broken: weights.size() != feature_names.size()+1: {} vs {}",
            model.weights.len(),
            model.feature_names.len() + 1
        )));
    }
    if !matches!(model.normalizer_type.as_str(), "stdv" | "uni" | "none") {
        return Err(Error::InvalidValue(
            "PercolatorModel normalizer must be stdv, uni or none".into(),
        ));
    }
    if let Some(name) = model
        .feature_names
        .iter()
        .find(|n| n.contains(['\t', '\n', '\r']) || n.starts_with('#'))
    {
        return Err(Error::InvalidValue(format!(
            "feature name {name:?} cannot be stored in a Percolator model file"
        )));
    }
    let bias = model.weights[model.feature_names.len()];
    let mut text = String::new();
    text.push_str("# OpenMS Percolator model file\n");
    text.push_str(&format!("format_version: {}\n", model.format_version));
    text.push_str(&format!("normalizer: {}\n", model.normalizer_type));
    text.push_str(&format!("seed: {}\n", model.seed));
    text.push_str(&format!("n_features: {}\n", model.feature_names.len()));
    text.push_str(&format!("bias: {}\n", format_g17(bias)));
    text.push_str(FEATURES_COMMENT);
    text.push('\n');
    for (name, weight) in model.feature_names.iter().zip(&model.weights) {
        text.push_str(&format!("{name}\t{}\n", format_g17(*weight)));
    }
    Ok(text)
}

/// `std::ostream << double` at precision 17 (`%.17g`).
fn format_g17(value: f64) -> String {
    if value.is_nan() {
        return if value.is_sign_negative() {
            "-nan"
        } else {
            "nan"
        }
        .into();
    }
    if value.is_infinite() {
        return if value < 0.0 { "-inf" } else { "inf" }.into();
    }
    if value == 0.0 {
        return if value.is_sign_negative() { "-0" } else { "0" }.into();
    }
    let scientific = format!("{value:.16e}");
    let (mantissa, exponent) = scientific.split_once('e').unwrap_or((&scientific, "0"));
    let exponent: i32 = exponent.parse().unwrap_or(0);
    if !(-4..17).contains(&exponent) {
        let mantissa = trim_fraction(mantissa);
        let sign = if exponent < 0 { '-' } else { '+' };
        format!("{mantissa}e{sign}{:02}", exponent.abs())
    } else {
        let decimals = (16 - exponent).max(0) as usize;
        trim_fraction(&format!("{value:.decimals$}")).to_owned()
    }
}

fn trim_fraction(text: &str) -> &str {
    if text.contains('.') {
        text.trim_end_matches('0').trim_end_matches('.')
    } else {
        text
    }
}

/// Deserialize a model written by [`save_model`] (`Percolator::loadModel`).
///
/// # Errors
///
/// [`Error::Io`] when the file cannot be opened (the source's
/// `FileNotFound`), and [`Error::Parse`] as [`read_model`].
pub fn load_model(path: impl AsRef<Path>) -> Result<PercolatorModel> {
    let file = std::fs::File::open(path).map_err(Error::Io)?;
    read_model(BufReader::new(file))
}

/// [`load_model`] from a reader.
///
/// The reader is strict, as the source's: unknown or duplicate header
/// keys, a missing required key (`format_version`, `normalizer`, `seed`,
/// `n_features`, `bias`), a `format_version` other than 1, a normalizer
/// other than `stdv`/`uni`/`none`, an unparsable value, and a declared
/// `n_features` that differs from the number of feature rows are rejected.
/// Blank lines and lines starting with `#` are skipped; a line with a tab
/// is a feature row, any other line a `key: value` header. Numbers are read
/// as `std::stoi`/`std::stol`/`std::stod` read them: leading whitespace and
/// trailing text after the number are ignored. The bias is appended as the
/// last weight.
///
/// # Errors
///
/// [`Error::Parse`] with the one-based line number for every violation
/// above (line 0 for the whole-file checks); [`Error::Io`] when reading
/// fails.
pub fn read_model(reader: impl BufRead) -> Result<PercolatorModel> {
    const KNOWN: [&str; 5] = ["format_version", "normalizer", "seed", "n_features", "bias"];
    let parse = |line: usize, message: String| Error::Parse { line, message };
    let mut model = PercolatorModel {
        format_version: 0,
        ..PercolatorModel::default()
    };
    let mut seen: Vec<&str> = Vec::new();
    let mut declared: i64 = -1;
    let mut bias = None;
    for (index, line) in reader.lines().enumerate() {
        let number = index + 1;
        let line = line.map_err(Error::Io)?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some(tab) = line.find('\t') {
            let name = &line[..tab];
            let weight = stod(&line[tab + 1..])
                .ok_or_else(|| parse(number, format!("cannot parse weight on line '{line}'")))?;
            model.feature_names.push(name.to_owned());
            model.weights.push(weight);
            continue;
        }
        let colon = line.find(':').ok_or_else(|| {
            parse(
                number,
                format!("expected 'key: value' or 'name<TAB>weight', got '{line}'"),
            )
        })?;
        let key = trim_ws(&line[..colon]);
        let value = trim_ws(&line[colon + 1..]);
        let Some(known) = KNOWN.iter().copied().find(|k| *k == key) else {
            return Err(parse(number, format!("unknown header key '{key}'")));
        };
        if seen.contains(&known) {
            return Err(parse(number, format!("duplicate header key '{key}'")));
        }
        seen.push(known);
        let cannot = || parse(number, format!("cannot parse value for '{key}'"));
        match known {
            "format_version" => model.format_version = stoi(value).ok_or_else(cannot)?,
            "seed" => model.seed = stoi(value).ok_or_else(cannot)?,
            "n_features" => declared = stol(value).ok_or_else(cannot)?,
            "bias" => bias = Some(stod(value).ok_or_else(cannot)?),
            _ => {
                if !matches!(value, "stdv" | "uni" | "none") {
                    return Err(parse(
                        number,
                        "invalid normalizer value (expected stdv|uni|none)".into(),
                    ));
                }
                model.normalizer_type = value.to_owned();
            }
        }
    }
    for required in KNOWN {
        if !seen.contains(&required) {
            return Err(parse(
                0,
                format!("missing required header key '{required}'"),
            ));
        }
    }
    if model.format_version != 1 {
        return Err(parse(0, "unsupported format_version (expected 1)".into()));
    }
    if declared < 0 || declared as u64 != model.feature_names.len() as u64 {
        return Err(parse(
            0,
            format!(
                "declared n_features does not match feature row count: {declared} vs {}",
                model.feature_names.len()
            ),
        ));
    }
    let bias = bias.ok_or_else(|| parse(0, "'bias' header key not found".into()))?;
    model.weights.push(bias);
    Ok(model)
}

/// The loader's `trim`: strip spaces, tabs, carriage returns and line feeds.
fn trim_ws(text: &str) -> &str {
    text.trim_matches([' ', '\t', '\r', '\n'])
}

/// The integer prefix `std::stol` reads: C whitespace, an optional sign,
/// then decimal digits; `None` without digits or on overflow.
fn stol(text: &str) -> Option<i64> {
    let text = text.trim_start_matches([' ', '\t', '\n', '\x0b', '\x0c', '\r']);
    let bytes = text.as_bytes();
    let mut end = usize::from(matches!(bytes.first(), Some(b'+' | b'-')));
    let digits_start = end;
    while end < bytes.len() && bytes[end].is_ascii_digit() {
        end += 1;
    }
    if end == digits_start {
        return None;
    }
    text[..end].parse::<i64>().ok()
}

/// `std::stoi`: [`stol`] limited to `int`.
fn stoi(text: &str) -> Option<i32> {
    stol(text).and_then(|v| i32::try_from(v).ok())
}

/// The floating prefix `std::stod` reads: C whitespace, an optional sign,
/// then `inf`, `infinity` or `nan` (any case), or decimal digits with an
/// optional fraction and exponent; `None` when no number starts the text.
/// Hexadecimal floats, which `strtod` also reads, are not accepted.
fn stod(text: &str) -> Option<f64> {
    let text = text.trim_start_matches([' ', '\t', '\n', '\x0b', '\x0c', '\r']);
    let bytes = text.as_bytes();
    let mut i = usize::from(matches!(bytes.first(), Some(b'+' | b'-')));
    let negative = bytes.first() == Some(&b'-');
    let rest = text[i..].to_ascii_lowercase();
    if rest.starts_with("inf") {
        return Some(if negative {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        });
    }
    if rest.starts_with("nan") {
        return Some(if negative { -f64::NAN } else { f64::NAN });
    }
    let digit = |i: usize| bytes.get(i).is_some_and(u8::is_ascii_digit);
    let mut count = 0;
    while digit(i) {
        i += 1;
        count += 1;
    }
    if bytes.get(i) == Some(&b'.') {
        i += 1;
        while digit(i) {
            i += 1;
            count += 1;
        }
    }
    if count == 0 {
        return None;
    }
    if matches!(bytes.get(i), Some(b'e' | b'E')) {
        let mut j = i + 1;
        if matches!(bytes.get(j), Some(b'+' | b'-')) {
            j += 1;
        }
        if digit(j) {
            while digit(j) {
                j += 1;
            }
            i = j;
        }
    }
    text[..i].parse::<f64>().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g17_matches_printf() {
        assert_eq!(format_g17(0.1), "0.10000000000000001");
        assert_eq!(format_g17(1.0), "1");
        assert_eq!(format_g17(-0.408004346094467), "-0.40800434609446701");
        assert_eq!(format_g17(1e-5), "1.0000000000000001e-05");
        assert_eq!(format_g17(1e20), "1e+20");
        assert_eq!(format_g17(123456.0), "123456");
    }

    #[test]
    fn number_prefixes_follow_the_c_library() {
        assert_eq!(stoi(" 12abc"), Some(12));
        assert_eq!(stoi("abc"), None);
        assert_eq!(stod("1.5\t2"), Some(1.5));
        assert_eq!(stod("1e"), Some(1.0));
        assert_eq!(stod("-inf"), Some(f64::NEG_INFINITY));
        assert_eq!(stod(""), None);
    }
}
