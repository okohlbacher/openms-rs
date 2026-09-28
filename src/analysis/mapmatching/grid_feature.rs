// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! A feature as QT feature grouping stores it in its grid.
//!
//! Native equivalent of `DATASTRUCTURES/GridFeature.h`. The header lives under
//! `DATASTRUCTURES`, but it needs `BaseFeature`, `AASequence` and
//! `PeptideIdentification`; `data_structures` may not depend on the kernel
//! (the kernel already reaches it), so the port lives beside its consumer,
//! `QTClusterFinder`. See `docs/CLUSTERING_SUPPORT.md`.

use std::cmp::Ordering;
use std::collections::BTreeSet;

use crate::chemistry::AASequence;
use crate::kernel::features::BaseFeature;

/// A peptide sequence ordered as the source's `AASequence::operator<` orders
/// it, as far as that order is reproducible.
///
/// The source's `std::set<AASequence>` orders by residue count first (a
/// shorter sequence is smaller), then puts a sequence without an N-terminal
/// modification before one with, then compares residue by residue. This
/// wrapper orders by residue count, then N-terminal modification presence,
/// then the unmodified one-letter string, then C-terminal modification
/// presence. Where the source then compares modification *pointers* (whose
/// order is the registry's allocation order) or modification ids, this falls
/// back to the crate's [`AASequence`] order, which is deterministic but not
/// the source's. Two keys are equal exactly when the sequences are equal.
///
/// The order is observable: it is the iteration order of
/// [`GridFeature::annotations`], and `QTCluster` breaks ties between equally
/// good annotations by it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SourceOrderedSequence(pub AASequence);

impl Ord for SourceOrderedSequence {
    fn cmp(&self, other: &Self) -> Ordering {
        let (a, b) = (&self.0, &other.0);
        a.len()
            .cmp(&b.len())
            .then_with(|| {
                a.n_terminal_modification()
                    .is_some()
                    .cmp(&b.n_terminal_modification().is_some())
            })
            .then_with(|| a.as_str().cmp(b.as_str()))
            .then_with(|| {
                a.c_terminal_modification()
                    .is_some()
                    .cmp(&b.c_terminal_modification().is_some())
            })
            .then_with(|| a.cmp(b))
    }
}

impl PartialOrd for SourceOrderedSequence {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Representation of a feature in a hash grid (source `GridFeature`).
///
/// A grid feature can be stored in a
/// [`HashGrid`](crate::math::ml::clustering::hash_grid::HashGrid) and refers to
/// a [`BaseFeature`] (a feature or a consensus feature). It is used for QT
/// feature grouping (`QTClusterFinder`). The source holds a reference to the
/// feature; so does this port, for the lifetime `'a`.
#[derive(Clone, Debug)]
pub struct GridFeature<'a> {
    /// The contained feature (source `feature_`).
    feature: &'a BaseFeature,
    /// Index of the feature map or consensus map (source `map_index_`).
    map_index: usize,
    /// Index of the feature in the map (source `feature_index_`).
    feature_index: usize,
    /// Peptide sequences annotated to the feature (source `annotations_`).
    annotations: BTreeSet<SourceOrderedSequence>,
}

impl<'a> GridFeature<'a> {
    /// A grid feature for `feature`, the `feature_index`-th feature of map
    /// `map_index` (source constructor).
    ///
    /// The annotations are the sequences of the first hit of each of the
    /// feature's peptide identifications; identifications without hits are
    /// skipped (the source comments that this "shouldn't be the case").
    pub fn new(feature: &'a BaseFeature, map_index: usize, feature_index: usize) -> Self {
        let annotations = feature
            .peptide_identifications
            .iter()
            .filter_map(|identification| identification.hits.first())
            .map(|hit| SourceOrderedSequence(hit.sequence.clone()))
            .collect();
        Self {
            feature,
            map_index,
            feature_index,
            annotations,
        }
    }

    /// The feature (source `getFeature`).
    pub fn feature(&self) -> &'a BaseFeature {
        self.feature
    }

    /// The map index (source `getMapIndex`).
    pub fn map_index(&self) -> usize {
        self.map_index
    }

    /// The feature index (source `getFeatureIndex`).
    pub fn feature_index(&self) -> usize {
        self.feature_index
    }

    /// The ID of the grid feature, the feature index (source `getID`). The
    /// source converts the `Size` index to `Int`; an index above `i32::MAX`
    /// wraps here as it does in the source's two's-complement conversion.
    pub fn id(&self) -> i32 {
        self.feature_index as i32
    }

    /// The peptide sequences annotated to the feature (source
    /// `getAnnotations`), in the source's set order; see
    /// [`SourceOrderedSequence`].
    pub fn annotations(&self) -> &BTreeSet<SourceOrderedSequence> {
        &self.annotations
    }

    /// The feature's retention time (source `getRT`).
    pub fn rt(&self) -> f64 {
        self.feature.rt
    }

    /// The feature's m/z (source `getMZ`).
    pub fn mz(&self) -> f64 {
        self.feature.mz
    }
}
