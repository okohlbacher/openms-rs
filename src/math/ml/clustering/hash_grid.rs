// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! A map from 2D coordinates to values, bucketed into grid cells.
//!
//! Native equivalent of `ML/CLUSTERING/HashGrid.h` (header-only template). The
//! source stores cells in a `boost::unordered_map` and each cell's pairs in a
//! `boost::unordered_multimap`, so its iteration order is whatever Boost's
//! hashing and bucket count produce. That order is not reproduced: cells
//! iterate in ascending cell index and pairs within a cell in insertion order.
//! Everything else (the cell of a coordinate, sizes, erasure counts, the grid
//! dimension) is the source's. See `docs/CLUSTERING_SUPPORT.md`.

use std::collections::BTreeMap;

use crate::data_structures::dposition::DPosition;
use crate::{Error, Result};

/// Coordinate of a stored pair (source `HashGrid::ClusterCenter`,
/// `DPosition<2, double>`).
pub type ClusterCenter = DPosition<2>;

/// Index of a cell (source `HashGrid::CellIndex`, `DPosition<2, Int64>`).
pub type CellIndex = [i64; 2];

/// Contents of one cell (source `HashGrid::CellContent`): its (coordinate,
/// value) pairs in insertion order.
pub type CellContent<C> = Vec<(ClusterCenter, C)>;

/// Container for (2-dimensional coordinate, value) pairs (source
/// `HashGrid<Cluster>`).
///
/// A hash grid consists of cells. The key of each cell is a pair of integers;
/// each pair is assigned to the cell holding its coordinate. The container
/// implements most of a standard multimap interface. `C` is the stored type
/// (source template parameter `Cluster`).
///
/// The source's element iterators are replaced by [`iter`](Self::iter) /
/// [`iter_mut`](Self::iter_mut), which yield each pair with its cell index
/// (the source iterator's `index()`), and by the `(cell, position)` location
/// [`insert`](Self::insert) returns and [`erase_at`](Self::erase_at) takes.
#[derive(Clone, Debug)]
pub struct HashGrid<C> {
    /// The non-empty and emptied cells (source `cells_`).
    cells: BTreeMap<CellIndex, CellContent<C>>,
    /// Upper-right corner of the key space for cells (source `grid_dimension_`).
    grid_dimension: CellIndex,
    /// Dimension of the cells (source public `cell_dimension`).
    cell_dimension: ClusterCenter,
}

impl<C> HashGrid<C> {
    /// An empty grid with cells of size `cell_dimension` (source constructor).
    /// The grid dimension starts at `(0, 0)`.
    pub fn new(cell_dimension: ClusterCenter) -> Self {
        Self {
            cells: BTreeMap::new(),
            grid_dimension: [0, 0],
            cell_dimension,
        }
    }

    /// Dimension of the cells (source public member `cell_dimension`).
    pub fn cell_dimension(&self) -> &ClusterCenter {
        &self.cell_dimension
    }

    /// Upper-right corner of the key space for cells (source public member
    /// `grid_dimension`): per dimension the largest cell index inserted so
    /// far, and at least 0. [`clear`](Self::clear) does not reset it, as in
    /// the source.
    pub fn grid_dimension(&self) -> &CellIndex {
        &self.grid_dimension
    }

    /// Inserts a (coordinate, value) pair (source `insert`) and returns its
    /// location: the cell index and the position within the cell.
    ///
    /// # Errors
    ///
    /// As [`cell_index_at_cluster_center`](Self::cell_index_at_cluster_center);
    /// the grid is unchanged on error.
    pub fn insert(&mut self, key: ClusterCenter, value: C) -> Result<(CellIndex, usize)> {
        let cell_key = self.cell_index_at_cluster_center(&key)?;
        let cell = self.cells.entry(cell_key).or_default();
        for (current, new) in self.grid_dimension.iter_mut().zip(cell_key) {
            if *current < new {
                *current = new;
            }
        }
        cell.push((key, value));
        Ok((cell_key, cell.len() - 1))
    }

    /// Erases the pair at `position` in cell `cell` (source `erase(iterator)`)
    /// and returns it. Later pairs of the cell move down by one position. The
    /// cell itself stays, empty, as in the source.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when no such pair exists (an invalid iterator
    /// in the source).
    pub fn erase_at(&mut self, cell: CellIndex, position: usize) -> Result<(ClusterCenter, C)> {
        let content = self
            .cells
            .get_mut(&cell)
            .filter(|content| position < content.len())
            .ok_or_else(|| Error::InvalidValue("HashGrid: no element at this location".into()))?;
        Ok(content.remove(position))
    }

    /// Erases every pair whose coordinate equals `key` and returns how many
    /// were erased (source `erase(const key_type&)`). Coordinates compare with
    /// `==`, so a NaN coordinate matches nothing.
    ///
    /// # Errors
    ///
    /// As [`cell_index_at_cluster_center`](Self::cell_index_at_cluster_center).
    pub fn erase(&mut self, key: &ClusterCenter) -> Result<usize> {
        let cell_key = self.cell_index_at_cluster_center(key)?;
        let Some(content) = self.cells.get_mut(&cell_key) else {
            return Ok(0);
        };
        let before = content.len();
        content.retain(|(center, _)| center != key);
        Ok(before - content.len())
    }

    /// Removes all cells (source `clear`).
    pub fn clear(&mut self) {
        self.cells.clear();
    }

    /// All pairs with their cell index: cells in ascending index, pairs in
    /// insertion order (source `begin()` / `end()`, whose order is Boost's
    /// hash order instead).
    pub fn iter(&self) -> impl Iterator<Item = (CellIndex, &(ClusterCenter, C))> + '_ {
        self.cells
            .iter()
            .flat_map(|(index, content)| content.iter().map(move |pair| (*index, pair)))
    }

    /// As [`iter`](Self::iter), with mutable values (source non-const
    /// iterator). The coordinates stay read-only, as the source's `value_type`
    /// has a `const` key.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = (CellIndex, &ClusterCenter, &mut C)> + '_ {
        self.cells.iter_mut().flat_map(|(index, content)| {
            content
                .iter_mut()
                .map(move |(center, value)| (*index, &*center, value))
        })
    }

    /// Whether the grid holds no pair (source `empty`).
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The number of pairs (source `size`).
    pub fn len(&self) -> usize {
        self.cells.values().map(Vec::len).sum()
    }

    /// All cells, including emptied ones, in ascending cell index (source
    /// `grid_begin()` / `grid_end()`, const and non-const, in Boost's order).
    pub fn grid_iter(&self) -> impl Iterator<Item = (&CellIndex, &CellContent<C>)> + '_ {
        self.cells.iter()
    }

    /// As [`grid_iter`](Self::grid_iter), with mutable cell contents (source
    /// non-const `grid_begin()`, "currently needed non-const by
    /// HierarchicalClustering").
    pub fn grid_iter_mut(
        &mut self,
    ) -> impl Iterator<Item = (&CellIndex, &mut CellContent<C>)> + '_ {
        self.cells.iter_mut()
    }

    /// The cell at index `x` (source const `grid_at`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when there is no such cell (source
    /// `std::out_of_range` from `unordered_map::at`).
    pub fn grid_at(&self, x: &CellIndex) -> Result<&CellContent<C>> {
        self.cells.get(x).ok_or_else(|| no_cell(x))
    }

    /// The cell at index `x`, mutable (source non-const `grid_at`).
    ///
    /// # Errors
    ///
    /// As [`grid_at`](Self::grid_at).
    pub fn grid_at_mut(&mut self, x: &CellIndex) -> Result<&mut CellContent<C>> {
        self.cells.get_mut(x).ok_or_else(|| no_cell(x))
    }

    /// The cell at index `x` if present (source const `grid_find`, which
    /// returns `grid_end()` otherwise).
    pub fn grid_find(&self, x: &CellIndex) -> Option<&CellContent<C>> {
        self.cells.get(x)
    }

    /// The cell at index `x` if present, mutable (source non-const
    /// `grid_find`).
    pub fn grid_find_mut(&mut self, x: &CellIndex) -> Option<&mut CellContent<C>> {
        self.cells.get_mut(x)
    }

    /// The cell index of the coordinate `key` (source
    /// `cellIndexAtClusterCenter`): per dimension `floor(key / cell_dimension)`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when a computed index is outside `i64` (source
    /// `Exception::OutOfRange`). The source compares against
    /// `numeric_limits<Int64>::max()` converted to `double`, which is `2^63`,
    /// so it lets `2^63` itself through to an undefined conversion; here
    /// `2^63` is refused too, and so is a NaN index (a NaN coordinate or a
    /// zero-by-zero division), which the source also converts undefinedly.
    pub fn cell_index_at_cluster_center(&self, key: &ClusterCenter) -> Result<CellIndex> {
        let mut index = [0_i64; 2];
        for (dimension, slot) in index.iter_mut().enumerate() {
            let t = (key[dimension] / self.cell_dimension[dimension]).floor();
            // -2^63 is exact in f64; 2^63 is the first value above i64::MAX.
            if t.is_nan() || t < -9_223_372_036_854_775_808.0 || t >= 9_223_372_036_854_775_808.0 {
                return Err(Error::InvalidRange(
                    "HashGrid: cell index exceeds the Int64 range".into(),
                ));
            }
            *slot = t as i64;
        }
        Ok(index)
    }
}

fn no_cell(x: &CellIndex) -> Error {
    Error::InvalidRange(format!("HashGrid: no cell at ({}, {})", x[0], x[1]))
}

/// Hash value of an integer position (source free function template
/// `hash_value(const DPosition<N, T>&)` for an integral `T`): the XOR of the
/// coordinates, since `boost::hash` of an integer is the integer itself.
///
/// As the source's class test notes, this hash is weak: `(1, 1)` and `(2, 2)`
/// both hash to 0. The floating-point instantiation, which Boost's container
/// uses for [`ClusterCenter`] keys, hashes through `boost::hash<double>` and is
/// not ported, because this port's grid does not hash.
pub fn hash_value(coordinates: &[u64]) -> u64 {
    coordinates
        .iter()
        .fold(0, |hash, coordinate| hash ^ coordinate)
}
