// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! One cluster of 2D points, as grid-based clustering builds it.
//!
//! Native equivalent of `ML/CLUSTERING/GridBasedCluster.h`. See
//! `docs/CLUSTERING_SUPPORT.md`.

use std::cmp::Ordering;

use crate::data_structures::dinterval::DIntervalBase;
use crate::data_structures::dposition::DPosition;

/// Centre of a cluster (source `GridBasedCluster::Point`, `DPosition<2>`).
pub type Point = DPosition<2>;

/// Bounding box of a cluster (source `GridBasedCluster::Rectangle`,
/// `DBoundingBox<2>`). `DBoundingBox` is `DIntervalBase` plus `enlarge`,
/// `encloses`, `intersects` and `isEmpty`, none of which this header's
/// interface exposes; the grid-based clustering enlarges boxes through
/// [`enlarge`].
pub type Rectangle = DIntervalBase<2>;

/// Source `DBoundingBox::enlarge(p)`: per dimension, lowers the minimum to
/// `p[i]` when `p[i]` is smaller and raises the maximum when it is larger.
pub fn enlarge(bounding_box: &mut Rectangle, p: &Point) {
    let mut min = *bounding_box.min_position();
    let mut max = *bounding_box.max_position();
    for i in 0..2 {
        if p[i] < min[i] {
            min[i] = p[i];
        }
        if p[i] > max[i] {
            max[i] = p[i];
        }
    }
    bounding_box.set_min_max(min, max);
}

/// Basic data structure for clustering (source `GridBasedCluster`).
///
/// Each point of a cluster can (optionally) carry two properties A and B. For
/// two points to be in the same cluster they need the same property A (for
/// example the same charge) and different properties B (for example
/// originating from different maps). `-1` means a property is not set.
///
/// **Comparison is by centre y only**, as in the source: `<`, `>` and `==`
/// compare `centre().y()` and ignore everything else, so two different
/// clusters with the same centre y are equal. `extendClustersY` relies on
/// this ordering.
#[derive(Clone, Debug)]
pub struct GridBasedCluster {
    /// Centre of the cluster (source `centre_`).
    centre: Point,
    /// Bounding box of the cluster, `(min, max)` in x and y (source
    /// `bounding_box_`).
    bounding_box: Rectangle,
    /// Indices of the points in the cluster (source `point_indices_`).
    point_indices: Vec<i32>,
    /// Property A shared by all points (source `property_A_`).
    property_a: i32,
    /// Property B of each point (source `properties_B_`).
    properties_b: Vec<i32>,
}

impl GridBasedCluster {
    /// A cluster with properties (source five-argument constructor). The
    /// arguments are stored as given; `properties_b` is not checked against
    /// `point_indices`, as in the source.
    pub fn new(
        centre: Point,
        bounding_box: Rectangle,
        point_indices: Vec<i32>,
        property_a: i32,
        properties_b: Vec<i32>,
    ) -> Self {
        Self {
            centre,
            bounding_box,
            point_indices,
            property_a,
            properties_b,
        }
    }

    /// A cluster without properties (source three-argument constructor):
    /// property A is `-1` and every point's property B is `-1`.
    pub fn without_properties(
        centre: Point,
        bounding_box: Rectangle,
        point_indices: Vec<i32>,
    ) -> Self {
        let properties_b = vec![-1; point_indices.len()];
        Self {
            centre,
            bounding_box,
            point_indices,
            property_a: -1,
            properties_b,
        }
    }

    /// The cluster centre (source `getCentre`).
    pub fn centre(&self) -> &Point {
        &self.centre
    }

    /// The bounding box (source `getBoundingBox`).
    pub fn bounding_box(&self) -> &Rectangle {
        &self.bounding_box
    }

    /// The indices of the points in the cluster (source `getPoints`).
    pub fn points(&self) -> &[i32] {
        &self.point_indices
    }

    /// Property A (source `getPropertyA`).
    pub fn property_a(&self) -> i32 {
        self.property_a
    }

    /// Property B of every point (source `getPropertiesB`).
    pub fn properties_b(&self) -> &[i32] {
        &self.properties_b
    }
}

impl PartialEq for GridBasedCluster {
    /// Source `operator==`: equal centre y. A NaN centre y is equal to nothing.
    fn eq(&self, other: &Self) -> bool {
        self.centre.y() == other.centre.y()
    }
}

impl PartialOrd for GridBasedCluster {
    /// Source `operator<` / `operator>`: by centre y.
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.centre.y().partial_cmp(&other.centre.y())
    }
}
