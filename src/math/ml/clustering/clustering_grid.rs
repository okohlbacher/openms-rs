// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! A rectilinear grid on which clusters of 2D points are registered.
//!
//! Native equivalent of `ML/CLUSTERING/ClusteringGrid.h`. See
//! `docs/CLUSTERING_SUPPORT.md`.

use std::collections::BTreeMap;

use crate::data_structures::dposition::DPosition;
use crate::math::libstdcxx::lower_bound;
use crate::{Error, Result};

/// Coordinates `(i, j)` of a grid cell (source `ClusteringGrid::CellIndex`,
/// `std::pair<int,int>`).
pub type CellIndex = (i32, i32);

/// Coordinates in the x-y plane (source `ClusteringGrid::Point`,
/// `DPosition<2>`).
pub type Point = DPosition<2>;

/// Data structure to store 2D data to be clustered, for example (m/z,
/// retention time) coordinates from multiplex filtering (source
/// `ClusteringGrid`).
///
/// The grid is given by its lines in x and in y; cell `(i, j)` lies between
/// lines `i - 1` and `i` in x and `j - 1` and `j` in y. Each non-empty cell
/// holds the indices of the clusters centred in it, in insertion order.
#[derive(Clone, Debug, PartialEq)]
pub struct ClusteringGrid {
    /// Grid lines in x (source `grid_spacing_x_`).
    grid_spacing_x: Vec<f64>,
    /// Grid lines in y (source `grid_spacing_y_`).
    grid_spacing_y: Vec<f64>,
    /// `[min, max]` of the grid in x (source `range_x_`): the first and last
    /// grid line.
    range_x: (f64, f64),
    /// `[min, max]` of the grid in y (source `range_y_`).
    range_y: (f64, f64),
    /// Grid cell index mapped to the clusters in it (source `cells_`, a
    /// `std::map` of `std::list<int>`).
    cells: BTreeMap<CellIndex, Vec<i32>>,
}

impl ClusteringGrid {
    /// A grid with lines `grid_spacing_x` and `grid_spacing_y` (source
    /// constructor).
    ///
    /// Source `@note`: the vectors are assumed to be sorted. As in the source
    /// this is not checked; [`index`](Self::index) on unsorted lines returns
    /// what the source's `std::lower_bound` returns on them.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when either vector is empty, where the source
    /// reads `front()` of an empty vector, or has more lines than an `int`
    /// cell index can address.
    pub fn new(grid_spacing_x: Vec<f64>, grid_spacing_y: Vec<f64>) -> Result<Self> {
        let range = |lines: &[f64], axis: &str| -> Result<(f64, f64)> {
            if i32::try_from(lines.len()).is_err() {
                return Err(Error::InvalidValue(format!(
                    "ClusteringGrid: too many grid lines in {axis}"
                )));
            }
            match (lines.first(), lines.last()) {
                (Some(&first), Some(&last)) => Ok((first, last)),
                _ => Err(Error::InvalidValue(format!(
                    "ClusteringGrid: no grid lines in {axis}"
                ))),
            }
        };
        let range_x = range(&grid_spacing_x, "x")?;
        let range_y = range(&grid_spacing_y, "y")?;
        Ok(Self {
            grid_spacing_x,
            grid_spacing_y,
            range_x,
            range_y,
            cells: BTreeMap::new(),
        })
    }

    /// The grid lines in x (source `getGridSpacingX`, which returns a copy).
    pub fn grid_spacing_x(&self) -> &[f64] {
        &self.grid_spacing_x
    }

    /// The grid lines in y (source `getGridSpacingY`, which returns a copy).
    pub fn grid_spacing_y(&self) -> &[f64] {
        &self.grid_spacing_y
    }

    /// Adds the cluster with index `cluster_index` (in the cluster list) to the
    /// cell `cell_index` (source `addCluster`). The cell is created if needed.
    pub fn add_cluster(&mut self, cell_index: CellIndex, cluster_index: i32) {
        self.cells
            .entry(cell_index)
            .or_default()
            .push(cluster_index);
    }

    /// Removes the cluster with index `cluster_index` from the cell
    /// `cell_index`, and removes the cell if no other cluster is left (source
    /// `removeCluster`). As `std::list::remove`, every occurrence of the index
    /// is removed. An absent cell is ignored.
    pub fn remove_cluster(&mut self, cell_index: CellIndex, cluster_index: i32) {
        if let Some(clusters) = self.cells.get_mut(&cell_index) {
            clusters.retain(|&index| index != cluster_index);
            if clusters.is_empty() {
                self.cells.remove(&cell_index);
            }
        }
    }

    /// Removes all clusters from the grid, and hence all cells (source
    /// `removeAllClusters`).
    pub fn remove_all_clusters(&mut self) {
        self.cells.clear();
    }

    /// The indices of the clusters centred in cell `cell_index`, in insertion
    /// order (source `getClusters`); empty for a cell that holds none, as the
    /// pinned source returns after its `map::end()` fix.
    pub fn clusters(&self, cell_index: CellIndex) -> Vec<i32> {
        self.cells.get(&cell_index).cloned().unwrap_or_default()
    }

    /// The index `(i, j)` of the cell in which `position` lies (source
    /// `getIndex`).
    ///
    /// `i` is the number of x grid lines not greater than `x`, as the source's
    /// `std::lower_bound` with `std::less_equal` computes it, so a point on a
    /// line belongs to the cell above it and a point on the first line to
    /// cell 0; likewise `j` in y.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `position` lies outside the grid, that is,
    /// outside `[first, last]` line in x or in y (source
    /// `Exception::IllegalArgument`, with the message "This position (x,y)=(x,y)
    /// is outside the range of the grid. (min <= x <= max, min <= y <= max)").
    /// A NaN coordinate compares false against both bounds and is accepted, as
    /// in the source, and lands in cell 0 of its axis.
    pub fn index(&self, position: &Point) -> Result<CellIndex> {
        let (x, y) = (position.x(), position.y());
        if x < self.range_x.0 || x > self.range_x.1 || y < self.range_y.0 || y > self.range_y.1 {
            return Err(Error::InvalidValue(format!(
                "This position (x,y)=({x},{y}) is outside the range of the grid. ({} <= x <= {}, {} <= y <= {})",
                self.range_x.0, self.range_x.1, self.range_y.0, self.range_y.1
            )));
        }
        // Both lengths fit `i32`, checked at construction.
        let i = lower_bound(&self.grid_spacing_x, |line| *line <= x) as i32;
        let j = lower_bound(&self.grid_spacing_y, |line| *line <= y) as i32;
        Ok((i, j))
    }

    /// Whether any cluster is registered in cell `cell_index` (source
    /// `isNonEmptyCell`).
    ///
    /// The source documents `Exception::IllegalArgument` and
    /// `Exception::InvalidValue` for this function, but its implementation is
    /// a plain lookup that never throws; so is this one.
    pub fn is_non_empty_cell(&self, cell_index: CellIndex) -> bool {
        self.cells.contains_key(&cell_index)
    }

    /// The number of grid cells occupied by one or more clusters (source
    /// `getCellCount`).
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }
}
