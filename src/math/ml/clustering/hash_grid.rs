// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! A map from 2D coordinates to values, bucketed into grid cells.
//!
//! Native equivalent of `ML/CLUSTERING/HashGrid.h` (header-only template). The
//! source stores cells in a `boost::unordered_map` and each cell's pairs in a
//! `boost::unordered_multimap`, so its iteration order is whatever Boost's
//! hashing, bucket counts and bucket groups produce. That order is
//! reproduced (through [`BoostTable`](crate::math::hash_order::BoostTable),
//! validated against the executed C++ by `tests/hash_order_oracle.rs`):
//! cells iterate in Boost 1.89's order for the hash
//! `boost::hash<Int64>(i) ^ boost::hash<Int64>(j)`, and the pairs within a
//! cell in its multimap's order for the hash of the two `double` bit patterns,
//! with equal keys adjacent. Everything else (the cell of a coordinate, sizes,
//! erasure counts, the grid dimension) is the source's. See
//! `docs/CLUSTERING_SUPPORT.md`.
//!
//! The stored pairs sit in vectors kept in iteration order: `insert` and the
//! erasures keep each cell's vector in step with the emulated table. The
//! `*_mut` accessors that hand out a cell's vector let a caller change values
//! in place; adding, removing or reordering pairs through them is not
//! supported.

use std::cell::OnceCell;
use std::collections::{BTreeMap, HashMap};

use crate::data_structures::dposition::DPosition;
use crate::math::hash_order::{BoostTable, boost_hash_f64, boost_hash_i64};
use crate::{Error, Result};

/// Coordinate of a stored pair (source `HashGrid::ClusterCenter`,
/// `DPosition<2, double>`).
pub type ClusterCenter = DPosition<2>;

/// Index of a cell (source `HashGrid::CellIndex`, `DPosition<2, Int64>`).
pub type CellIndex = [i64; 2];

/// Contents of one cell (source `HashGrid::CellContent`): its (coordinate,
/// value) pairs in the iteration order of Boost's `unordered_multimap`.
pub type CellContent<C> = Vec<(ClusterCenter, C)>;

/// `hash_value(const DPosition<2, Int64>&)`: the XOR of `boost::hash<Int64>`.
fn cell_hash(index: &CellIndex) -> u64 {
    boost_hash_i64(index[0]) ^ boost_hash_i64(index[1])
}

/// `hash_value(const DPosition<2, double>&)`: the XOR of `boost::hash<double>`.
fn center_hash(center: &ClusterCenter) -> u64 {
    hash_value_f64(&[center[0], center[1]])
}

/// One cell: its pairs in iteration order and the emulated multimap that
/// decides that order.
#[derive(Clone, Debug)]
struct Cell<C> {
    content: CellContent<C>,
    /// Table node of each entry of `content`.
    nodes: Vec<usize>,
    table: BoostTable<ClusterCenter>,
}

impl<C> Cell<C> {
    fn new() -> Self {
        Self {
            content: Vec::new(),
            nodes: Vec::new(),
            table: BoostTable::new(),
        }
    }
}

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
    /// The non-empty and emptied cells (source `cells_`), by cell index; the
    /// iteration order is `outer`'s.
    cells: BTreeMap<CellIndex, Cell<C>>,
    /// The emulated `boost::unordered_map` of the cells; node `k` is the
    /// `k`-th cell created since the last `clear`.
    outer: BoostTable<CellIndex>,
    /// The cell index of each node of `outer`.
    outer_cells: Vec<CellIndex>,
    /// The cell indices in iteration order, computed on first use after a
    /// change.
    order: OnceCell<Vec<CellIndex>>,
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
            outer: BoostTable::new(),
            outer_cells: Vec::new(),
            order: OnceCell::new(),
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
    /// location: the cell index and the position within the cell, which is
    /// its position in the cell's iteration order after the insertion (an
    /// insertion may reorder the cell's other pairs).
    ///
    /// # Errors
    ///
    /// As [`cell_index_at_cluster_center`](Self::cell_index_at_cluster_center);
    /// the grid is unchanged on error.
    pub fn insert(&mut self, key: ClusterCenter, value: C) -> Result<(CellIndex, usize)> {
        let cell_key = self.cell_index_at_cluster_center(&key)?;
        if !self.cells.contains_key(&cell_key) {
            let (node, _) = self.outer.insert_unique(cell_key, cell_hash(&cell_key));
            debug_assert_eq!(node, self.outer_cells.len());
            self.outer_cells.push(cell_key);
            self.cells.insert(cell_key, Cell::new());
            self.order = OnceCell::new();
        }
        for (current, new) in self.grid_dimension.iter_mut().zip(cell_key) {
            if *current < new {
                *current = new;
            }
        }
        let cell = self
            .cells
            .get_mut(&cell_key)
            .expect("the cell was created above");
        let node = cell.table.insert_equiv(key, center_hash(&key));
        // The new node, then the iteration order of the whole cell: a rehash
        // may permute the existing pairs.
        let order = cell.table.order();
        let mut by_node: HashMap<usize, (ClusterCenter, C)> = HashMap::with_capacity(order.len());
        for (id, pair) in cell.nodes.drain(..).zip(cell.content.drain(..)) {
            by_node.insert(id, pair);
        }
        by_node.insert(node, (key, value));
        let mut position = 0;
        for (index, id) in order.iter().enumerate() {
            if *id == node {
                position = index;
            }
            let pair = by_node.remove(id).expect("every live node holds a pair");
            cell.content.push(pair);
            cell.nodes.push(*id);
        }
        Ok((cell_key, position))
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
        let target = self
            .cells
            .get_mut(&cell)
            .filter(|target| position < target.content.len())
            .ok_or_else(|| Error::InvalidValue("HashGrid: no element at this location".into()))?;
        let node = target.nodes.remove(position);
        target.table.remove(node);
        Ok(target.content.remove(position))
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
        let Some(target) = self.cells.get_mut(&cell_key) else {
            return Ok(0);
        };
        let mut erased = 0;
        let mut index = 0;
        while index < target.content.len() {
            if target.content[index].0 == *key {
                let node = target.nodes.remove(index);
                target.table.remove(node);
                target.content.remove(index);
                erased += 1;
            } else {
                index += 1;
            }
        }
        Ok(erased)
    }

    /// Removes all cells (source `clear`). Like the source's
    /// `unordered_map::clear`, the cell table keeps its bucket count.
    pub fn clear(&mut self) {
        self.cells.clear();
        self.outer.clear();
        self.outer_cells.clear();
        self.order = OnceCell::new();
    }

    /// The cell indices in the iteration order of the cell map.
    fn cell_order(&self) -> &[CellIndex] {
        self.order.get_or_init(|| {
            self.outer
                .order()
                .into_iter()
                .map(|node| self.outer_cells[node])
                .collect()
        })
    }

    /// All pairs with their cell index: cells in the order of Boost's
    /// `unordered_map`, pairs within a cell in its `unordered_multimap`'s
    /// order (source `begin()` / `end()`).
    pub fn iter(&self) -> impl Iterator<Item = (CellIndex, &(ClusterCenter, C))> + '_ {
        self.cell_order().iter().flat_map(move |index| {
            self.cells[index]
                .content
                .iter()
                .map(move |pair| (*index, pair))
        })
    }

    /// As [`iter`](Self::iter), with mutable values (source non-const
    /// iterator). The coordinates stay read-only, as the source's `value_type`
    /// has a `const` key.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = (CellIndex, &ClusterCenter, &mut C)> + '_ {
        let order: Vec<CellIndex> = self.cell_order().to_vec();
        let position: HashMap<CellIndex, usize> = order
            .iter()
            .enumerate()
            .map(|(index, cell)| (*cell, index))
            .collect();
        let mut cells: Vec<(CellIndex, &mut Cell<C>)> = self
            .cells
            .iter_mut()
            .map(|(index, cell)| (*index, cell))
            .collect();
        cells.sort_by_key(|(index, _)| position[index]);
        cells.into_iter().flat_map(|(index, cell)| {
            cell.content
                .iter_mut()
                .map(move |(center, value)| (index, &*center, value))
        })
    }

    /// Whether the grid holds no pair (source `empty`).
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The number of pairs (source `size`).
    pub fn len(&self) -> usize {
        self.cells.values().map(|cell| cell.content.len()).sum()
    }

    /// All cells, including emptied ones, in the iteration order of the cell
    /// map (source `grid_begin()` / `grid_end()`).
    pub fn grid_iter(&self) -> impl Iterator<Item = (&CellIndex, &CellContent<C>)> + '_ {
        self.cell_order().iter().map(move |index| {
            let (key, cell) = self
                .cells
                .get_key_value(index)
                .expect("every ordered cell exists");
            (key, &cell.content)
        })
    }

    /// As [`grid_iter`](Self::grid_iter), with mutable cell contents (source
    /// non-const `grid_begin()`, "currently needed non-const by
    /// HierarchicalClustering").
    pub fn grid_iter_mut(
        &mut self,
    ) -> impl Iterator<Item = (&CellIndex, &mut CellContent<C>)> + '_ {
        let position: HashMap<CellIndex, usize> = self
            .cell_order()
            .iter()
            .enumerate()
            .map(|(index, cell)| (*cell, index))
            .collect();
        let mut cells: Vec<(&CellIndex, &mut CellContent<C>)> = self
            .cells
            .iter_mut()
            .map(|(index, cell)| (index, &mut cell.content))
            .collect();
        cells.sort_by_key(|(index, _)| position[*index]);
        cells.into_iter()
    }

    /// The cell at index `x` (source const `grid_at`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when there is no such cell (source
    /// `std::out_of_range` from `unordered_map::at`).
    pub fn grid_at(&self, x: &CellIndex) -> Result<&CellContent<C>> {
        self.cells
            .get(x)
            .map(|cell| &cell.content)
            .ok_or_else(|| no_cell(x))
    }

    /// The cell at index `x`, mutable (source non-const `grid_at`).
    ///
    /// # Errors
    ///
    /// As [`grid_at`](Self::grid_at).
    pub fn grid_at_mut(&mut self, x: &CellIndex) -> Result<&mut CellContent<C>> {
        self.cells
            .get_mut(x)
            .map(|cell| &mut cell.content)
            .ok_or_else(|| no_cell(x))
    }

    /// The cell at index `x` if present (source const `grid_find`, which
    /// returns `grid_end()` otherwise).
    pub fn grid_find(&self, x: &CellIndex) -> Option<&CellContent<C>> {
        self.cells.get(x).map(|cell| &cell.content)
    }

    /// The cell at index `x` if present, mutable (source non-const
    /// `grid_find`).
    pub fn grid_find_mut(&mut self, x: &CellIndex) -> Option<&mut CellContent<C>> {
        self.cells.get_mut(x).map(|cell| &mut cell.content)
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
            // [-2^63, 2^63) is exactly the i64 range in f64; NaN is outside it.
            if !(-9_223_372_036_854_775_808.0..9_223_372_036_854_775_808.0).contains(&t) {
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
/// both hash to 0. The grid's cell index (`DPosition<2, Int64>`) uses this
/// instantiation; its [`ClusterCenter`] keys use [`hash_value_f64`].
pub fn hash_value(coordinates: &[u64]) -> u64 {
    coordinates
        .iter()
        .fold(0, |hash, coordinate| hash ^ coordinate)
}

/// Hash value of a floating-point position (the same template for
/// `T = double`): the XOR of `boost::hash<double>` of the coordinates, which
/// is the bit pattern of `value + 0` (`-0.0` hashes as `+0.0`).
pub fn hash_value_f64(coordinates: &[f64]) -> u64 {
    coordinates
        .iter()
        .fold(0, |hash, coordinate| hash ^ boost_hash_f64(*coordinate))
}
