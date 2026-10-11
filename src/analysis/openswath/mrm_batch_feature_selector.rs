// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Iterative application of an [`MRMFeatureSelector`].
//!
//! Ports `ANALYSIS/OPENSWATH/MRMBatchFeatureSelector.h` and its `.cpp`: the
//! class with deleted constructor and destructor and three static functions.
//! [`batch_mrm_features`] calls a selector once per parameter set, feeding
//! the selection of each cycle to the next; [`batch_mrm_features_score`] and
//! [`batch_mrm_features_qmip`] do so with the two selectors of
//! [`mrm_feature_selector`](crate::analysis::openswath::mrm_feature_selector).
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/MRM_FEATURE_SELECTOR_SUPPORT.md).
//!
//! # Example
//!
//! ```
//! use openms::analysis::openswath::mrm_batch_feature_selector::batch_mrm_features_score;
//! use openms::analysis::openswath::mrm_feature_selector::{PartitionSolver, SelectorParameters};
//! use openms::kernel::FeatureMap;
//!
//! // no parameter sets: the selection is empty, as the source clears it first
//! let selected = batch_mrm_features_score(&FeatureMap::default(), &[], PartitionSolver).unwrap();
//! assert!(selected.is_empty());
//! ```
//!
//! [`MRMFeatureSelector`]: crate::analysis::openswath::mrm_feature_selector::MRMFeatureSelector
//! [`batch_mrm_features`]: crate::analysis::openswath::mrm_batch_feature_selector::batch_mrm_features
//! [`batch_mrm_features_score`]: crate::analysis::openswath::mrm_batch_feature_selector::batch_mrm_features_score
//! [`batch_mrm_features_qmip`]: crate::analysis::openswath::mrm_batch_feature_selector::batch_mrm_features_qmip

use crate::Result;
use crate::analysis::openswath::mrm_feature_selector::{
    LpSolver, MRMFeatureSelector, MRMFeatureSelectorQMIP, MRMFeatureSelectorScore,
    SelectorParameters,
};
use crate::kernel::FeatureMap;

/// Call `feature_selector.select_mrm_feature` once per entry of `parameters`,
/// the result of each cycle being the input of the next.
///
/// With no parameter sets the result is an empty map (the source clears its
/// output first and never fills it).
///
/// # Errors
///
/// The first error of a selector call; nothing is returned in that case (the
/// source leaves its out-parameter at the last completed cycle).
pub fn batch_mrm_features(
    feature_selector: &dyn MRMFeatureSelector,
    features: &FeatureMap,
    parameters: &[SelectorParameters],
) -> Result<FeatureMap> {
    let mut input_features = features.clone();
    let mut selected_features = FeatureMap::default();
    for params in parameters {
        selected_features = feature_selector.select_mrm_feature(&input_features, params)?;
        input_features = selected_features.clone();
    }
    Ok(selected_features)
}

/// [`batch_mrm_features`] with a [`MRMFeatureSelectorScore`] solving with `solver`.
///
/// # Errors
///
/// As [`batch_mrm_features`].
pub fn batch_mrm_features_score<S: LpSolver>(
    features: &FeatureMap,
    parameters: &[SelectorParameters],
    solver: S,
) -> Result<FeatureMap> {
    batch_mrm_features(&MRMFeatureSelectorScore::new(solver), features, parameters)
}

/// [`batch_mrm_features`] with a [`MRMFeatureSelectorQMIP`] solving with `solver`.
///
/// # Errors
///
/// As [`batch_mrm_features`].
pub fn batch_mrm_features_qmip<S: LpSolver>(
    features: &FeatureMap,
    parameters: &[SelectorParameters],
    solver: S,
) -> Result<FeatureMap> {
    batch_mrm_features(&MRMFeatureSelectorQMIP::new(solver), features, parameters)
}
