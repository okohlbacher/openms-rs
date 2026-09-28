// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! A symmetric distance matrix that stores only the cells below the diagonal.
//!
//! Native equivalent of `DATASTRUCTURES/DistanceMatrix.h` (header-only
//! template). The source keeps an array of row pointers, row `i` holding the
//! `i` cells `(i, 0) .. (i, i - 1)`; this port keeps the same triangular
//! layout as one `Vec` per row, so every index the source computes addresses
//! the same cell here. See `docs/CLUSTERING_SUPPORT.md`.

use std::fmt;

use crate::concept::fuzzy_string_comparator::format_g;
use crate::{Error, Result};

/// Largest number of stored (below-diagonal) cells a matrix may hold.
///
/// `dimension * (dimension - 1) / 2` cells are allocated up front; the source
/// allocates without a bound and throws `Exception::OutOfMemory` only when the
/// allocator fails. 2^28 cells are 1 GiB of `f32`, a 23,170-element matrix.
pub const MAX_DISTANCE_MATRIX_CELLS: usize = 1 << 28;

/// A two-dimensional distance matrix, similar to a general matrix.
///
/// Only the elements below the main diagonal are stored. Access to an element
/// above the diagonal is translated to its mirror below the diagonal, and the
/// main diagonal reads as zero because a self-distance is assumed to be zero.
///
/// The coordinates of the minimal element are tracked for as long as values are
/// written only through [`set`](Self::set). The other mutators
/// ([`set_quick`](Self::set_quick), [`reduce`](Self::reduce)) may leave the
/// tracked minimum stale; call [`update_min_element`](Self::update_min_element)
/// before relying on it again, exactly as the source documents.
///
/// `V` is the cell type (source template parameter `Value`); the clustering
/// code uses `f32`. `V::default()` is the value read on the diagonal and must
/// therefore be zero, which it is for every numeric type.
#[derive(Clone, Debug)]
pub struct DistanceMatrix<V> {
    /// Row `i` holds cells `(i, 0) .. (i, i - 1)`; row 0 is empty (source
    /// `matrix_`, whose row 0 is a null pointer).
    rows: Vec<Vec<V>>,
    /// Number of accessible rows, which is also the number of columns (source
    /// `dimensionsize_`).
    dimension: usize,
    /// Row and column of the minimal element, row greater than column (source
    /// `min_element_`).
    min_element: (usize, usize),
}

impl<V> Default for DistanceMatrix<V> {
    /// An empty matrix with no rows (source default constructor); the tracked
    /// minimum is `(0, 0)`.
    fn default() -> Self {
        Self {
            rows: Vec::new(),
            dimension: 0,
            min_element: (0, 0),
        }
    }
}

fn out_of_range(what: &str) -> Error {
    Error::InvalidRange(format!("DistanceMatrix: {what}"))
}

impl<V: Copy + PartialOrd + Default> DistanceMatrix<V> {
    /// An empty matrix (source default constructor).
    pub fn new() -> Self {
        Self::default()
    }

    /// A `dimension` x `dimension` matrix whose every off-diagonal cell holds
    /// `value` (source detailed constructor; the source's default `value` is
    /// `Value()`, pass `V::default()` for it). The main diagonal still reads as
    /// zero. The tracked minimum is `(1, 0)`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the matrix would hold more than
    /// [`MAX_DISTANCE_MATRIX_CELLS`] cells or the allocation fails; the source
    /// throws `Exception::OutOfMemory` for a failed allocation.
    pub fn with_value(dimension: usize, value: V) -> Result<Self> {
        let rows = allocate_rows(dimension, value)?;
        Ok(Self {
            rows,
            dimension,
            min_element: (1, 0),
        })
    }

    /// The value at row `i`, column `j` (source `getValue` and `operator()`,
    /// all four overloads).
    ///
    /// A cell above the diagonal reads its mirror, a diagonal cell reads
    /// `V::default()` (zero).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when `i` or `j` is not below
    /// [`dimension_size`](Self::dimension_size) (source
    /// `Exception::OutOfRange`).
    pub fn get(&self, i: usize, j: usize) -> Result<V> {
        if i >= self.dimension || j >= self.dimension {
            return Err(out_of_range("coordinates out of range"));
        }
        if i == j {
            return Ok(V::default());
        }
        let (row, column) = if i < j { (j, i) } else { (i, j) };
        self.cell(row, column)
    }

    /// Stores `value` at `(i, j)` and keeps the tracked minimum up to date
    /// (source `setValue`).
    ///
    /// Writing to a diagonal cell does nothing. When the written cell is not
    /// the tracked minimum, a smaller value becomes the new minimum. When it is
    /// the tracked minimum and the value grows, the minimum is recomputed with
    /// [`update_min_element`](Self::update_min_element). This is the source
    /// after the fix of issue #9488 (DATA-24), which the pinned class test
    /// checks.
    ///
    /// The source reads the tracked minimum's cell without checking it; after
    /// [`reduce`](Self::reduce) that cell may no longer exist, and the source
    /// then reads freed memory. Here a tracked minimum outside the current
    /// triangle is recomputed first.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when `i` or `j` is out of range (source
    /// `Exception::OutOfRange`).
    pub fn set(&mut self, i: usize, j: usize, value: V) -> Result<()> {
        if i >= self.dimension || j >= self.dimension {
            return Err(out_of_range("coordinates out of range"));
        }
        if i == j {
            return Ok(());
        }
        let (row, column) = if i < j { (j, i) } else { (i, j) };
        if !self.is_stored(self.min_element) {
            self.update_min_element()?;
        }
        let (min_row, min_column) = self.min_element;
        if row != min_row || column != min_column {
            *self.cell_mut(row, column)? = value;
            if value < self.cell(min_row, min_column)? {
                self.min_element = (row, column);
            }
        } else if value <= self.cell(min_row, min_column)? {
            *self.cell_mut(row, column)? = value;
        } else {
            *self.cell_mut(row, column)? = value;
            self.update_min_element()?;
        }
        Ok(())
    }

    /// Stores `value` at `(i, j)` without maintaining the tracked minimum
    /// (source `setValueQuick`).
    ///
    /// Source `@note`: this may invalidate the tracked minimum; call
    /// [`update_min_element`](Self::update_min_element) before using it
    /// again. Writing to a diagonal cell does nothing.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when `i` or `j` is out of range (source
    /// `Exception::OutOfRange`).
    pub fn set_quick(&mut self, i: usize, j: usize, value: V) -> Result<()> {
        if i >= self.dimension || j >= self.dimension {
            return Err(out_of_range("coordinates out of range"));
        }
        if i == j {
            return Ok(());
        }
        let (row, column) = if i < j { (j, i) } else { (i, j) };
        *self.cell_mut(row, column)? = value;
        Ok(())
    }

    /// Resets to an empty matrix with the tracked minimum at `(0, 0)` (source
    /// `clear`).
    pub fn clear(&mut self) {
        self.rows = Vec::new();
        self.min_element = (0, 0);
        self.dimension = 0;
    }

    /// Replaces the content with a `dimension` x `dimension` matrix filled with
    /// `value` (source `resize`, whose default `value` is `Value()`).
    ///
    /// Source `@note`: all previous content is invalidated. The tracked
    /// minimum becomes `(1, 0)`.
    ///
    /// # Errors
    ///
    /// As [`with_value`](Self::with_value). The matrix is unchanged on error;
    /// the source leaves an empty matrix behind when an allocation fails.
    pub fn resize(&mut self, dimension: usize, value: V) -> Result<()> {
        let rows = allocate_rows(dimension, value)?;
        self.rows = rows;
        self.dimension = dimension;
        self.min_element = (1, 0);
        Ok(())
    }

    /// Removes row `j` and then column `j`, shrinking the matrix by one
    /// dimension (source `reduce`).
    ///
    /// Source `@note`: this may invalidate the tracked minimum; call
    /// [`update_min_element`](Self::update_min_element) before using it.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when `j` is not below
    /// [`dimension_size`](Self::dimension_size) (source
    /// `Exception::OutOfRange`).
    pub fn reduce(&mut self, j: usize) -> Result<()> {
        if j >= self.dimension {
            return Err(out_of_range("row to remove out of range"));
        }
        // Rows below `j` are untouched; every later row loses its column `j`
        // and moves up by one, which is what the source's copy loop does. The
        // former row `j + 1` keeps `j` cells, so row 0 stays empty.
        self.rows.remove(j);
        for row in self.rows.iter_mut().skip(j) {
            if j < row.len() {
                row.remove(j);
            }
        }
        self.dimension -= 1;
        Ok(())
    }

    /// Number of rows, equal to the number of columns (source
    /// `dimensionsize`).
    pub fn dimension_size(&self) -> usize {
        self.dimension
    }

    /// Recomputes the tracked minimum (source `updateMinElement`).
    ///
    /// The scan starts at `(1, 0)` and replaces the candidate only with a
    /// strictly smaller cell, row by row with `std::min_element` within each
    /// row, so ties keep the earliest row and, within a row, the earliest
    /// column. A matrix of dimension 1 keeps `(1, 0)`, a cell that does not
    /// exist, as the source does.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] for an empty matrix (source
    /// `Exception::OutOfRange`).
    pub fn update_min_element(&mut self) -> Result<()> {
        self.min_element = (1, 0);
        if self.dimension < 1 {
            return Err(out_of_range("no element to access"));
        }
        if self.dimension != 1 {
            for r in 2..self.dimension {
                let row = &self.rows[r];
                let mut best = 0;
                for (c, value) in row.iter().enumerate().skip(1) {
                    if *value < row[best] {
                        best = c;
                    }
                }
                let current = self.cell(self.min_element.0, self.min_element.1)?;
                if row[best] < current {
                    self.min_element = (r, best);
                }
            }
        }
        Ok(())
    }

    /// Row and column of the minimal element, the row always greater than the
    /// column (source `getMinElementCoordinates`).
    ///
    /// The value is only as current as the last [`set`](Self::set),
    /// [`resize`](Self::resize) or
    /// [`update_min_element`](Self::update_min_element); see the type
    /// documentation.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] for an empty matrix (source
    /// `Exception::OutOfRange`).
    pub fn min_element_coordinates(&self) -> Result<(usize, usize)> {
        if self.dimension == 0 {
            return Err(out_of_range("no element to access"));
        }
        Ok(self.min_element)
    }

    fn is_stored(&self, (row, column): (usize, usize)) -> bool {
        row < self.dimension && column < row
    }

    fn cell(&self, row: usize, column: usize) -> Result<V> {
        self.rows
            .get(row)
            .filter(|_| row < self.dimension)
            .and_then(|cells| cells.get(column))
            .copied()
            .ok_or_else(|| out_of_range("cell is not stored"))
    }

    fn cell_mut(&mut self, row: usize, column: usize) -> Result<&mut V> {
        let dimension = self.dimension;
        self.rows
            .get_mut(row)
            .filter(|_| row < dimension)
            .and_then(|cells| cells.get_mut(column))
            .ok_or_else(|| out_of_range("cell is not stored"))
    }
}

fn allocate_rows<V: Copy>(dimension: usize, value: V) -> Result<Vec<Vec<V>>> {
    let cells = dimension
        .checked_mul(dimension.saturating_sub(1))
        .map(|product| product / 2)
        .filter(|&cells| cells <= MAX_DISTANCE_MATRIX_CELLS)
        .ok_or_else(|| {
            Error::InvalidValue(format!(
                "DistanceMatrix: dimension {dimension} exceeds {MAX_DISTANCE_MATRIX_CELLS} cells"
            ))
        })?;
    let _ = cells;
    let oom = || Error::InvalidValue("DistanceMatrix: out of memory".into());
    let mut rows = Vec::new();
    rows.try_reserve_exact(dimension).map_err(|_| oom())?;
    for i in 0..dimension {
        let mut row = Vec::new();
        row.try_reserve_exact(i).map_err(|_| oom())?;
        row.resize(i, value);
        rows.push(row);
    }
    Ok(rows)
}

impl<V: PartialEq> PartialEq for DistanceMatrix<V> {
    /// Cell-wise equality (source `operator==`).
    ///
    /// Matrices of different dimension always compare unequal, as the pinned
    /// source does after its out-of-bounds fix. The tracked minimum is not
    /// compared. A NaN cell makes two matrices unequal.
    fn eq(&self, other: &Self) -> bool {
        if self.dimension != other.dimension {
            return false;
        }
        (1..self.dimension).all(|i| self.rows[i][..i] == other.rows[i][..i])
    }
}

impl<V: Copy + PartialOrd + Default + Into<f64>> fmt::Display for DistanceMatrix<V> {
    /// Prints the full square matrix (source `operator<<`).
    ///
    /// Every cell is written as `std::ostream` writes it at precision 6
    /// (`%g`) and followed by a tab; each row ends with a newline. As in the
    /// source, the stream width of 8 applies to the first value only (the
    /// width is consumed by the first insertion), so cell `(0, 0)` is padded to
    /// eight characters. The source colours the diagonal red only when the
    /// stream is a terminal `std::cout`/`std::cerr`; a formatter is never one,
    /// so no escape codes are written.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.dimension {
            for j in 0..self.dimension {
                let value: f64 = self.get(i, j).map_err(|_| fmt::Error)?.into();
                let text = format_g(value);
                if i == 0 && j == 0 {
                    write!(f, "{text:>8}\t")?;
                } else {
                    write!(f, "{text}\t")?;
                }
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reduce_first_row_keeps_triangle() {
        let mut dm = DistanceMatrix::with_value(3, 1.0_f32).unwrap();
        dm.set(1, 0, 2.0).unwrap();
        dm.set(2, 0, 3.0).unwrap();
        dm.set(2, 1, 4.0).unwrap();
        dm.reduce(0).unwrap();
        assert_eq!(dm.dimension_size(), 2);
        assert_eq!(dm.get(1, 0).unwrap(), 4.0);
    }
}
