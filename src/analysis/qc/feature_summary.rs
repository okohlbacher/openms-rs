// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Detected compounds as a metabolomics QC metric (source
//! `QC/FeatureSummary.h`): a summary of the features of a featureXML file.
//! See `docs/QC_SUPPORT.md`.

use super::qc_base::{QCBase, Requires, Status};
use crate::kernel::FeatureMap;
use crate::{Error, Result};

/// The summary values [`FeatureSummary::compute`] calculates.
///
/// Source `FeatureSummary::Result`; the default is a count of 0 and a mean of
/// 0, and equality compares both fields, source `operator==`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FeatureSummaryResult {
    /// The number of (top-level) features, the detected compounds.
    pub feature_count: u32,
    /// The mean of the features' `rt_deviation` meta values, in seconds; 0
    /// when no feature carries one.
    pub rt_shift_mean: f32,
}

/// Detected compounds as a metabolomics QC metric.
///
/// Source `FeatureSummary`, a stateless metric that summarises the features
/// of a featureXML file.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FeatureSummary;

impl FeatureSummary {
    /// Source `getName`.
    pub const NAME: &'static str = "Summary of features from featureXML file";

    /// A new metric, source `FeatureSummary()`.
    pub fn new() -> Self {
        Self
    }

    /// Summarise `feature_map`, source `compute`: the number of features and
    /// the mean of the `rt_deviation` meta value over the features that carry
    /// it.
    ///
    /// The mean is accumulated in `f32`, each value narrowed to `f32` first,
    /// and divided in `f32`, exactly as the source's `float` arithmetic does.
    /// Subordinate features are neither counted nor averaged.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a feature's `rt_deviation` is not
    /// numeric (a string, a list or empty), where the source's `DataValue`
    /// conversion to `float` throws `Exception::ConversionError`, and when the
    /// map holds more than `u32::MAX` features, which the source's `UInt`
    /// count would truncate.
    pub fn compute(&self, feature_map: &FeatureMap) -> Result<FeatureSummaryResult> {
        let feature_count = u32::try_from(feature_map.features.len()).map_err(|_| {
            Error::InvalidValue("feature count exceeds the source's 32-bit counter".into())
        })?;
        let mut sum_rt_deviations = 0.0_f32;
        let mut rt_count = 0_u32;
        for feature in &feature_map.features {
            if let Some(value) = feature.metadata.get("rt_deviation") {
                sum_rt_deviations += value.as_f64()? as f32;
                rt_count += 1;
            }
        }
        let rt_shift_mean = if rt_count != 0 {
            sum_rt_deviations / rt_count as f32
        } else {
            0.0
        };
        Ok(FeatureSummaryResult {
            feature_count,
            rt_shift_mean,
        })
    }
}

impl QCBase for FeatureSummary {
    fn name(&self) -> &str {
        Self::NAME
    }
    /// Features with unfiltered peptide identifications, source
    /// `requirements`.
    fn requirements(&self) -> Status {
        Status::from_flag(Requires::PreFdrFeat)
    }
}
