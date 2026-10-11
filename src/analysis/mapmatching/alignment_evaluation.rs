// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Evaluation of an alignment against a ground truth ("CAAP" evaluation).
//!
//! Ports `ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithm.h`,
//! `MapAlignmentEvaluationAlgorithmPrecision.h` and
//! `MapAlignmentEvaluationAlgorithmRecall.h` with their `.cpp` files. The
//! source's abstract base class becomes the [`MapAlignmentEvaluationAlgorithm`](crate::analysis::mapmatching::alignment_evaluation::MapAlignmentEvaluationAlgorithm)
//! trait, whose provided [`is_same_handle`](crate::analysis::mapmatching::alignment_evaluation::MapAlignmentEvaluationAlgorithm::is_same_handle)
//! is the base class's non-virtual `isSameHandle`; the two subclasses become
//! [`MapAlignmentEvaluationPrecision`](crate::analysis::mapmatching::alignment_evaluation::MapAlignmentEvaluationPrecision) and [`MapAlignmentEvaluationRecall`](crate::analysis::mapmatching::alignment_evaluation::MapAlignmentEvaluationRecall).
//! See `docs/MAP_ALIGNMENT_SUPPORT.md`.

use crate::kernel::{ConsensusFeature, ConsensusMap, FeatureHandle};
use crate::{Error, Result};

/// An algorithm that scores an input consensus map against a ground truth
/// (source `MapAlignmentEvaluationAlgorithm`).
pub trait MapAlignmentEvaluationAlgorithm {
    /// Score `input` against `ground_truth` (source `evaluate`, whose `out`
    /// parameter is the return value).
    ///
    /// Two handles are the same feature when [`is_same_handle`] says so with
    /// the given deviations: `rt_dev` and `mz_dev` in the maps' units,
    /// `int_dev` in intensity, and `use_charge` to require equal charges.
    ///
    /// [`is_same_handle`]: MapAlignmentEvaluationAlgorithm::is_same_handle
    ///
    /// # Errors
    ///
    /// As the implementation documents.
    fn evaluate(
        &self,
        input: &ConsensusMap,
        ground_truth: &ConsensusMap,
        rt_dev: f64,
        mz_dev: f64,
        int_dev: f32,
        use_charge: bool,
    ) -> Result<f64>;

    /// Whether two handles are the same feature (source `isSameHandle`),
    /// judged by retention time, m/z and intensity rather than by map and
    /// element index: they are, unless the RT difference exceeds `rt_dev`,
    /// the m/z difference exceeds `mz_dev`, the intensity difference exceeds
    /// `int_dev`, or `use_charge` is set and the charges differ.
    ///
    /// The intensity difference is taken in single precision, as the source
    /// computes it on `Peak2D::IntensityType`. A NaN difference exceeds
    /// nothing, as in the source.
    fn is_same_handle(
        &self,
        lhs: &FeatureHandle,
        rhs: &FeatureHandle,
        rt_dev: f64,
        mz_dev: f64,
        int_dev: f32,
        use_charge: bool,
    ) -> bool {
        is_same_handle(lhs, rhs, rt_dev, mz_dev, int_dev, use_charge)
    }
}

/// [`MapAlignmentEvaluationAlgorithm::is_same_handle`] as a free function.
pub fn is_same_handle(
    lhs: &FeatureHandle,
    rhs: &FeatureHandle,
    rt_dev: f64,
    mz_dev: f64,
    int_dev: f32,
    use_charge: bool,
) -> bool {
    if (lhs.rt - rhs.rt).abs() > rt_dev {
        return false;
    }
    if (lhs.mz - rhs.mz).abs() > mz_dev {
        return false;
    }
    if (lhs.intensity - rhs.intensity).abs() > int_dev {
        return false;
    }
    !(use_charge && lhs.charge != rhs.charge)
}

/// The ground-truth consensus features that group at least two elements;
/// both measures ignore the others.
fn usable_ground_truth(ground_truth: &ConsensusMap) -> Vec<&ConsensusFeature> {
    ground_truth
        .features
        .iter()
        .filter(|f| f.len() >= 2)
        .collect()
}

/// Number of handles of `truth` that match at least one handle of `tool`.
fn shared_handles(
    truth: &ConsensusFeature,
    tool: &ConsensusFeature,
    rt_dev: f64,
    mz_dev: f64,
    int_dev: f32,
    use_charge: bool,
) -> usize {
    truth
        .handles()
        .iter()
        .filter(|gt| {
            tool.handles()
                .iter()
                .any(|t| is_same_handle(t, gt, rt_dev, mz_dev, int_dev, use_charge))
        })
        .count()
}

/// Alignment precision (source `MapAlignmentEvaluationAlgorithmPrecision`).
///
/// With `N` ground-truth consensus features of at least two elements,
/// `precision = 1/N * sum_i (shared_i / tool_i)`: for each of them, `shared_i`
/// sums, over the input consensus features of at least two elements that
/// share any element with it, the number of shared elements, and `tool_i` the
/// sizes of those input features. A ground-truth feature with nothing shared
/// contributes zero.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MapAlignmentEvaluationPrecision;

impl MapAlignmentEvaluationAlgorithm for MapAlignmentEvaluationPrecision {
    /// # Errors
    ///
    /// None; always `Ok`. A ground truth with no consensus feature of at
    /// least two elements yields NaN, which is what the source computes
    /// (`1.0 / 0 * 0`); unlike recall, the source does not guard it.
    fn evaluate(
        &self,
        input: &ConsensusMap,
        ground_truth: &ConsensusMap,
        rt_dev: f64,
        mz_dev: f64,
        int_dev: f32,
        use_charge: bool,
    ) -> Result<f64> {
        let truth = usable_ground_truth(ground_truth);
        let mut sum = 0.0;
        for gt in &truth {
            let (mut shared_total, mut tool_total) = (0usize, 0usize);
            for tool in &input.features {
                let shared = shared_handles(gt, tool, rt_dev, mz_dev, int_dev, use_charge);
                if tool.len() >= 2 && shared > 0 {
                    shared_total += shared;
                    tool_total += tool.len();
                }
            }
            if shared_total != 0 {
                sum += shared_total as f64 / tool_total as f64;
            }
        }
        Ok((1.0 / truth.len() as f64) * sum)
    }
}

/// Alignment recall (source `MapAlignmentEvaluationAlgorithmRecall`).
///
/// With `N` ground-truth consensus features of at least two elements,
/// `recall = 1/N * sum_i (shared_i / (m_i * gt_i))`: for each of them,
/// `shared_i` sums the shared elements over the input consensus features of at
/// least two elements that share any, `m_i` counts those input features, and
/// `gt_i` is the size of the ground-truth feature.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MapAlignmentEvaluationRecall;

impl MapAlignmentEvaluationAlgorithm for MapAlignmentEvaluationRecall {
    /// An empty input map recalls nothing and scores 0.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the ground truth has no consensus
    /// feature of at least two elements (the source's
    /// `Exception::InvalidParameter`: recall is undefined then, also when both
    /// maps are empty).
    fn evaluate(
        &self,
        input: &ConsensusMap,
        ground_truth: &ConsensusMap,
        rt_dev: f64,
        mz_dev: f64,
        int_dev: f32,
        use_charge: bool,
    ) -> Result<f64> {
        let truth = usable_ground_truth(ground_truth);
        if truth.is_empty() {
            return Err(Error::InvalidValue(
                "Recall is undefined: the ground-truth consensus map contains no consensus feature with at least two elements.".into(),
            ));
        }
        if input.features.is_empty() {
            return Ok(0.0);
        }
        let mut sum = 0.0;
        for gt in &truth {
            let (mut shared_total, mut matched) = (0usize, 0usize);
            for tool in &input.features {
                let shared = shared_handles(gt, tool, rt_dev, mz_dev, int_dev, use_charge);
                if tool.len() >= 2 && shared > 0 {
                    shared_total += shared;
                    matched += 1;
                }
            }
            // The source counts every ground-truth element once per input
            // feature and divides by the input size again, which is exactly
            // the ground-truth feature's size.
            let gt_size = gt.len();
            if shared_total != 0 {
                sum += shared_total as f64 / (matched * gt_size) as f64;
            }
        }
        Ok((1.0 / truth.len() as f64) * sum)
    }
}
