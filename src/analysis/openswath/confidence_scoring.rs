// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Confidence scoring of SRM/MRM/PRM features against a targeted assay library.
//!
//! Ports `ANALYSIS/OPENSWATH/ConfidenceScoring.h` and its `.cpp`. For each
//! detected feature [`ConfidenceScoring::score_map`] scores the feature
//! against the matching assay (its `PeptideRef` metavalue) and against decoy
//! assays drawn from the same library, with a fitted binomial GLM over the
//! normalised retention time difference and the Manhattan distance of the
//! transition intensities, and writes `GLM_score`, `local_FDR` and an overall
//! quality of `1 - local_FDR` back to the feature.
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/MRM_FEATURE_FINDER_SUPPORT.md).
//!
//! # Example
//!
//! ```
//! use openms::analysis::openswath::confidence_scoring::ConfidenceScoring;
//! use openms::analysis::openswath::reaction_monitoring_transition::ReactionMonitoringTransition;
//! use openms::analysis::openswath::targeted_experiment::TargetedExperiment;
//! use openms::analysis::openswath::targeted_experiment_helper::Peptide;
//! use openms::analysis::transformations::TransformationDescription;
//! use openms::kernel::{Feature, FeatureMap};
//!
//! let mut library = TargetedExperiment::new();
//! for (id, rt, transition) in [("pep_1", 50.0, "tr_10"), ("pep_2", 60.0, "tr_20")] {
//!     let mut peptide = Peptide::default();
//!     peptide.id = id.into();
//!     let mut retention_time = openms::analysis::openswath::targeted_experiment_helper::RetentionTime::default();
//!     retention_time.set_rt(rt);
//!     peptide.rts.push(retention_time);
//!     library.add_peptide(peptide);
//!     let mut t = ReactionMonitoringTransition::new();
//!     t.set_native_id(transition);
//!     t.set_peptide_ref(id);
//!     t.set_product_mz(500.0);
//!     t.set_library_intensity(500.0);
//!     library.add_transition(t);
//! }
//! let mut features = FeatureMap::default();
//! for (peptide, transition) in [("pep_1", "tr_10"), ("pep_2", "tr_20")] {
//!     let mut f = Feature::default();
//!     f.rt = 60.0;
//!     f.metadata.insert("PeptideRef".into(), peptide.into());
//!     let mut sub = Feature::default();
//!     sub.intensity = 1.0;
//!     sub.mz = 500.0;
//!     sub.metadata.insert("native_id".into(), transition.into());
//!     f.subordinates.push(sub);
//!     features.features.push(f);
//! }
//! let mut scoring = ConfidenceScoring::new(false);
//! scoring.initialize(library, 0, 0, TransformationDescription::default());
//! scoring.initialize_glm(0.0, -1.0, -1.0);
//! scoring.score_map(&mut features).unwrap();
//! // pep_2 has its expected RT at the feature's RT: GLM score 0.5, FDR 0
//! assert_eq!(features.features[1].quality, 1.0);
//! ```
//!
//! [`ConfidenceScoring::score_map`]: crate::analysis::openswath::confidence_scoring::ConfidenceScoring::score_map

use crate::analysis::openswath::mrm_feature_support::{meta_string, native_id};
use crate::analysis::openswath::scoring::normalize_sum;
use crate::analysis::openswath::targeted_experiment::TargetedExperiment;
use crate::analysis::openswath::targeted_experiment_helper::Peptide;
use crate::analysis::transformations::TransformationDescription;
use crate::concept::math_functions::RandomShuffler;
use crate::kernel::{Feature, FeatureMap};
use crate::metadata::MetaValue;
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};

/// The binomial GLM `sigmoid(intercept + rt_coef * diff_rt^2 + int_coef * dist_int)`
/// (source `ConfidenceScoring::GLM_`).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Glm {
    /// The intercept.
    pub intercept: f64,
    /// The coefficient of the squared normalised retention time difference.
    pub rt_coef: f64,
    /// The coefficient of the Manhattan intensity distance.
    pub int_coef: f64,
}

impl Glm {
    /// Evaluate at (`diff_rt`, `dist_int`); a probability in `[0, 1]`.
    pub fn evaluate(&self, diff_rt: f64, dist_int: f64) -> f64 {
        let lm = self.intercept + self.rt_coef * diff_rt * diff_rt + self.int_coef * dist_int;
        1.0 / (1.0 + (-lm).exp())
    }
}

/// Mapping of retention times into `[0, 100]` by the library's retention time
/// range (source `ConfidenceScoring::RTNorm_`).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RtNorm {
    /// The smallest assay retention time of the library.
    pub min_rt: f64,
    /// The largest assay retention time of the library.
    pub max_rt: f64,
}

impl RtNorm {
    /// `(rt - min_rt) / (max_rt - min_rt) * 100`.
    pub fn apply(&self, rt: f64) -> f64 {
        (rt - self.min_rt) / (self.max_rt - self.min_rt) * 100.0
    }
}

/// The `boost::bimap<double, multiset_of<double>>` the source collects
/// intensities in: unique keys on the left (a repeated key is ignored) and
/// the right view ordered by intensity, ties in insertion order.
#[derive(Default)]
struct IntensityMap {
    // (key, value, insertion sequence)
    entries: Vec<(f64, f64, usize)>,
    next: usize,
}

impl IntensityMap {
    fn insert(&mut self, key: f64, value: f64) {
        if self.entries.iter().any(|e| e.0 == key) {
            return;
        }
        self.entries.push((key, value, self.next));
        self.next += 1;
    }

    fn len(&self) -> usize {
        self.entries.len()
    }

    /// `intensity_map.right.erase(intensity_map.right.begin())`.
    fn erase_lowest(&mut self) {
        let lowest = self
            .entries
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.1.total_cmp(&b.1).then(a.2.cmp(&b.2)))
            .map(|(i, _)| i);
        if let Some(i) = lowest {
            self.entries.remove(i);
        }
    }

    /// The entries by ascending key: the left view.
    fn by_key(&self) -> Vec<(f64, f64)> {
        let mut sorted: Vec<(f64, f64)> = self.entries.iter().map(|e| (e.0, e.1)).collect();
        sorted.sort_by(|a, b| a.0.total_cmp(&b.0));
        sorted
    }
}

/// `extractIntensities_`: keep the `n_transitions` highest intensities (all
/// for `0`) and list them ordered by key, negative values (missing, `-1`)
/// raised to zero.
fn extract_intensities(intensity_map: &mut IntensityMap, n_transitions: usize) -> Vec<f64> {
    if n_transitions > 0 {
        let diff = intensity_map.len() as i64 - n_transitions as i64;
        for _ in 0..diff.max(0) {
            intensity_map.erase_lowest();
        }
    }
    intensity_map
        .by_key()
        .into_iter()
        .map(|(_, intensity)| 0.0_f64.max(intensity))
        .collect()
}

fn illegal(message: &str) -> Error {
    Error::InvalidValue(message.to_owned())
}

/// Scores features against a targeted assay library (source
/// `ConfidenceScoring`).
///
/// Call [`initialize`](Self::initialize) and
/// [`initialize_glm`](Self::initialize_glm), then
/// [`score_map`](Self::score_map). Like the source, a scorer keeps state
/// between `score_map` calls: the transition lookup grows by the library's
/// transitions on every call, and `n_decoys` is clamped for good.
#[derive(Clone, Debug)]
pub struct ConfidenceScoring {
    glm: Glm,
    rt_norm: RtNorm,
    library: TargetedExperiment,
    decoy_index: Vec<i32>,
    n_decoys: usize,
    transition_map: BTreeMap<String, Vec<i32>>,
    n_transitions: usize,
    rt_trafo: TransformationDescription,
    shuffler: RandomShuffler,
}

impl ConfidenceScoring {
    /// An empty scorer.
    ///
    /// With `test_mode` false the decoy shuffler is seeded with `0`, so the
    /// decoy draws are reproducible across runs and platforms; with `test_mode`
    /// true it is seeded with the current time in seconds, as the source's
    /// `time(nullptr)`. The members the source leaves uninitialised start at
    /// zero.
    pub fn new(test_mode: bool) -> Self {
        let seed = if test_mode {
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map_or(0, |d| d.as_secs())
        } else {
            0
        };
        Self {
            glm: Glm::default(),
            rt_norm: RtNorm::default(),
            library: TargetedExperiment::new(),
            decoy_index: Vec::new(),
            n_decoys: 0,
            transition_map: BTreeMap::new(),
            n_transitions: 0,
            rt_trafo: TransformationDescription::default(),
            shuffler: RandomShuffler::with_seed(seed),
        }
    }

    /// The configured GLM.
    pub fn glm(&self) -> Glm {
        self.glm
    }

    /// The retention time normalisation of the last [`score_map`](Self::score_map).
    pub fn rt_norm(&self) -> RtNorm {
        self.rt_norm
    }

    /// The number of decoys per feature, `0` meaning all unrelated assays.
    pub fn n_decoys(&self) -> usize {
        self.n_decoys
    }

    /// Install the library, the number of decoys, the number of transitions
    /// kept for the intensity distance (`0` keeps all) and the retention time
    /// transformation applied to the features.
    pub fn initialize(
        &mut self,
        library: TargetedExperiment,
        n_decoys: usize,
        n_transitions: usize,
        rt_trafo: TransformationDescription,
    ) {
        self.library = library;
        self.n_decoys = n_decoys;
        self.n_transitions = n_transitions;
        self.rt_trafo = rt_trafo;
    }

    /// Install the GLM coefficients fitted on a training set.
    pub fn initialize_glm(&mut self, intercept: f64, rt_coef: f64, int_coef: f64) {
        self.glm = Glm {
            intercept,
            rt_coef,
            int_coef,
        };
    }

    /// The retention time of an assay (source `getAssayRT_`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the assay has no retention time (the
    /// source's `Exception::IllegalArgument`; its debug-build precondition
    /// does not exist in a Release build).
    pub fn get_assay_rt(&self, assay: &Peptide) -> Result<f64> {
        assay.retention_time()
    }

    /// The Manhattan (L1) distance of two vectors (source `manhattanDist_`),
    /// over the length of `x`.
    ///
    /// # Panics
    ///
    /// Never; a `y` shorter than `x` counts the missing entries as zero (the
    /// source reads past the end).
    pub fn manhattan_dist(&self, x: &[f64], y: &[f64]) -> f64 {
        let mut sum = 0.0;
        for (i, &xi) in x.iter().enumerate() {
            sum += (xi - y.get(i).copied().unwrap_or(0.0)).abs();
        }
        sum
    }

    fn choose_decoys(&mut self) -> Result<()> {
        if self.n_decoys == 0 {
            return Ok(()); // the list is already initialised
        }
        let mut budget = self.decoy_index.len().saturating_mul(64).saturating_add(64);
        self.shuffler
            .shuffle_limited(&mut self.decoy_index, &mut budget)
    }

    /// Score one feature against one assay (source `scoreAssay_`); returns
    /// the GLM output.
    ///
    /// The assay's transitions are reduced to those named in `transition_ids`
    /// (all when empty), matched to the feature's intensities by product m/z,
    /// and the Manhattan distance of the sum-normalised intensity vectors and
    /// the normalised retention time difference go through the GLM.
    /// `feature_intensities` is normalised in place.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the feature has no intensities, when the
    /// assay and the feature differ in the number of usable transitions, and
    /// when the assay has no retention time.
    pub fn score_assay(
        &self,
        assay: &Peptide,
        feature_rt: f64,
        feature_intensities: &mut [f64],
        transition_ids: &BTreeSet<String>,
    ) -> Result<f64> {
        let assay_rt = self.rt_norm.apply(self.get_assay_rt(assay)?);
        let diff_rt = assay_rt - feature_rt;
        let mut intensity_map = IntensityMap::default();
        if let Some(indices) = self.transition_map.get(&assay.id) {
            for &index in indices {
                let transition = usize::try_from(index)
                    .ok()
                    .and_then(|i| self.library.transitions().get(i))
                    .ok_or_else(|| illegal("transition index out of range"))?;
                // for the "true" assay the same transitions as the feature's
                if !transition_ids.is_empty() && !transition_ids.contains(transition.native_id()) {
                    continue;
                }
                intensity_map.insert(transition.product_mz(), transition.library_intensity());
            }
        }
        let mut assay_intensities =
            extract_intensities(&mut intensity_map, feature_intensities.len());
        if feature_intensities.is_empty() {
            return Err(illegal(
                "Feature intensities were empty - please provide feature subordinate with intensities",
            ));
        }
        if feature_intensities.len() != assay_intensities.len() {
            return Err(illegal(
                "Did not find a feature for each assay provided - each feature needs to have n subordinates with the meta-value 'native_id' set to the corresponding transition.",
            ));
        }
        normalize_sum(feature_intensities);
        normalize_sum(&mut assay_intensities);
        let dist_int = self.manhattan_dist(feature_intensities, &assay_intensities);
        Ok(self.glm.evaluate(diff_rt, dist_int))
    }

    /// Score one feature against its true assay and a sample of decoys
    /// (source `scoreFeature_`), writing `GLM_score`, `local_FDR` and the
    /// overall quality.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the feature has no `PeptideRef`, names an
    /// assay absent from the library or carries no usable intensities, and the
    /// errors of [`score_assay`](Self::score_assay).
    pub fn score_feature(&mut self, feature: &mut Feature) -> Result<()> {
        // predictors from the feature
        let feature_rt = self.rt_norm.apply(self.rt_trafo.apply(feature.rt)?);
        let mut intensity_map = IntensityMap::default();
        // for the true assay the same transitions as the feature's are used
        let mut trans_id_map: BTreeMap<u64, (f64, String)> = BTreeMap::new();
        for sub in &feature.subordinates {
            intensity_map.insert(sub.mz, f64::from(sub.intensity));
            // std::map::operator[]: the last subordinate of an m/z wins
            trans_id_map.insert(sub.mz.to_bits(), (sub.mz, native_id(sub)));
        }
        let mut feature_intensities = extract_intensities(&mut intensity_map, self.n_transitions);
        // `intensity_map` now holds only the transitions needed later
        let mut transition_ids = BTreeSet::new();
        for (mz, _) in intensity_map.by_key() {
            transition_ids.insert(
                trans_id_map
                    .get(&mz.to_bits())
                    .map(|(_, id)| id.clone())
                    .unwrap_or_default(),
            );
        }

        let mut scores: Vec<f64> = Vec::new(); // the true score first
        let true_id = meta_string(&feature.metadata, "PeptideRef");
        if true_id.is_empty() {
            return Err(illegal(
                "Feature does not contain meta value 'PeptideRef' (reference to assay)",
            ));
        }
        let true_assay = self.library.peptide_by_ref(&true_id)?.clone();
        scores.push(self.score_assay(
            &true_assay,
            feature_rt,
            &mut feature_intensities,
            &transition_ids,
        )?);

        // decoy assays
        self.choose_decoys()?;
        let mut counter = 0usize;
        let decoys: Vec<i32> = self.decoy_index.clone();
        for decoy_index in decoys {
            let decoy_assay = usize::try_from(decoy_index)
                .ok()
                .and_then(|i| self.library.peptides().get(i))
                .ok_or_else(|| illegal("decoy index out of range"))?
                .clone();
            // skip the true assay and assays with too few transitions
            let n_transitions = self.transition_map.get(&decoy_assay.id).map_or(0, Vec::len);
            if decoy_assay.id == true_id || n_transitions < feature_intensities.len() {
                continue;
            }
            scores.push(self.score_assay(
                &decoy_assay,
                feature_rt,
                &mut feature_intensities,
                &BTreeSet::new(),
            )?);
            counter += 1;
            if self.n_decoys > 0 && counter >= self.n_decoys {
                break; // enough decoys
            }
        }

        // count decoy scores greater than the true score
        let greater = scores[1..].iter().filter(|&&s| s > scores[0]).count();
        let n_scores = scores.len();
        feature
            .metadata
            .insert("GLM_score".into(), MetaValue::source_float(scores[0]));
        // 0 / 0 is NaN when no decoy qualified, as in the source
        let local_fdr = greater as f64 / (n_scores as f64 - 1.0);
        feature
            .metadata
            .insert("local_FDR".into(), MetaValue::source_float(local_fdr));
        feature.quality = (1.0 - local_fdr) as f32;
        Ok(())
    }

    /// Score every feature of `features` in place.
    ///
    /// Each feature must carry the assay ID in `PeptideRef` and one
    /// subordinate per transition of the assay with `native_id` set. Writes
    /// `GLM_score`, `local_FDR` and sets the overall quality to
    /// `1 - local_FDR`; the library's retention time range is computed and
    /// cached on the way. When `n_decoys` exceeds the number of unrelated
    /// assays it is clamped to `0` (all of them).
    ///
    /// The source reports progress through a `ProgressLogger`; no progress is
    /// reported here.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the library has fewer than two peptides
    /// ("There need to be at least 2 assays in the library for
    /// ConfidenceScoring."), and the errors of
    /// [`score_feature`](Self::score_feature). Features scored before an
    /// error keep their results, as in the source.
    pub fn score_map(&mut self, features: &mut FeatureMap) -> Result<()> {
        // are there enough assays in the library?
        let n_assays = self.library.peptides().len();
        if n_assays < 2 {
            return Err(illegal(
                "There need to be at least 2 assays in the library for ConfidenceScoring.",
            ));
        }
        // (n_assays - 1 < n_decoys: the source logs a warning here)
        if n_assays - 1 <= self.n_decoys {
            self.n_decoys = 0; // use all available assays
        }
        self.decoy_index = (0..n_assays)
            .map(|i| {
                i32::try_from(i).map_err(|_| illegal("the library has more assays than an Int"))
            })
            .collect::<Result<Vec<_>>>()?;

        // mapping between assays and transitions
        for (i, transition) in self.library.transitions().iter().enumerate() {
            let index = i32::try_from(i)
                .map_err(|_| illegal("the library has more transitions than an Int"))?;
            self.transition_map
                .entry(transition.peptide_ref().to_owned())
                .or_default()
                .push(index);
        }

        // min./max. RT in the library
        let mut min_rt = f64::INFINITY;
        let mut max_rt = f64::NEG_INFINITY;
        for peptide in self.library.peptides() {
            let current_rt = self.get_assay_rt(peptide)?;
            if current_rt == -1.0 {
                continue; // indicates a missing value
            }
            min_rt = min_rt.min(current_rt);
            max_rt = max_rt.max(current_rt);
        }
        self.rt_norm = RtNorm { min_rt, max_rt };

        for feature in &mut features.features {
            self.score_feature(feature)?;
        }
        Ok(())
    }
}

impl Default for ConfidenceScoring {
    fn default() -> Self {
        Self::new(false)
    }
}
