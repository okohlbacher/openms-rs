// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Hierarchical clustering of peak spectra through binned similarity.
//!
//! The binned-spectrum overload of `ClusterHierarchical::cluster`
//! (`ML/CLUSTERING/ClusterHierarchical.h`). The generic overload lives in
//! [`crate::math::ml::clustering::cluster_hierarchical`]; this one needs
//! [`BinnedSpectrum`] and [`BinnedSpectrumCompareFunctor`], which `math` may
//! not depend on. See `docs/CLUSTERING_SUPPORT.md`.

use crate::comparison::{BinConfig, BinUnit, BinnedSpectrum, BinnedSpectrumCompareFunctor};
use crate::data_structures::binary_tree_node::BinaryTreeNode;
use crate::data_structures::distance_matrix::DistanceMatrix;
use crate::math::ml::clustering::cluster_functor::ClusterFunctor;
use crate::math::ml::clustering::cluster_hierarchical::ClusterHierarchical;
use crate::{MSSpectrum, Result};

/// Clusters peak spectra by a binned similarity (source
/// `ClusterHierarchical::cluster(std::vector<PeakSpectrum>&, const
/// BinnedSpectrumCompareFunctor&, double sz, UInt sp, float offset, const
/// ClusterFunctor&, std::vector<BinaryTreeNode>&, DistanceMatrix<float>&)`).
///
/// Each spectrum is binned with bin size `size`, bin spread `spread` and bin
/// offset `offset` in absolute (Th) units, so that `comparator` can be
/// applied. The source passes its `double` bin size to the `float` parameter
/// of the `BinnedSpectrum` constructor; `size` is converted to `f32` the same
/// way. The source's low-resolution defaults are a bin size of `1.0005`
/// (`BinnedSpectrum::DEFAULT_BIN_WIDTH_LOWRES`) and an offset of `0.4`
/// (`DEFAULT_BIN_OFFSET_LOWRES`).
///
/// `original_distance` is always re-created, whatever its size: filled with 1
/// and then set, below the diagonal, to `1 - comparator(b_i, b_j)` through
/// [`DistanceMatrix::set`], which maintains the tracked minimum. The
/// similarity must be normalized to `[0, 1]`. Clustering stops at `clusterer`
/// reaching the threshold of `hierarchical`, converted to `f32`. The result is
/// the tree for further investigation with
/// [`ClusterAnalyzer`](crate::math::ml::clustering::cluster_analyzer::ClusterAnalyzer).
///
/// # Errors
///
/// The errors of [`BinnedSpectrum::new`] (unsorted peaks, an invalid layout,
/// work beyond the [`BinConfig`] limits), of `comparator` and of `clusterer`
/// (fewer than two spectra), and of [`DistanceMatrix::resize`].
#[allow(clippy::too_many_arguments, reason = "source signature")]
pub fn cluster_binned_spectra(
    hierarchical: &ClusterHierarchical,
    data: &[MSSpectrum],
    comparator: &dyn BinnedSpectrumCompareFunctor,
    size: f64,
    spread: u32,
    offset: f32,
    clusterer: &dyn ClusterFunctor,
    original_distance: &mut DistanceMatrix<f32>,
) -> Result<Vec<BinaryTreeNode>> {
    let config = BinConfig {
        size: size as f32,
        unit: BinUnit::Absolute,
        spread,
        offset,
        ..BinConfig::default()
    };
    let mut binned_data = Vec::with_capacity(data.len());
    for spectrum in data {
        binned_data.push(BinnedSpectrum::new(spectrum, config)?);
    }

    original_distance.clear();
    original_distance.resize(data.len(), 1.0)?;
    for i in 0..binned_data.len() {
        for j in 0..i {
            let similarity = comparator.score(&binned_data[i], &binned_data[j])?;
            original_distance.set(i, j, (1.0 - similarity) as f32)?;
        }
    }

    clusterer.cluster(original_distance, hierarchical.threshold() as f32)
}
