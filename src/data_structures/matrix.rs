// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A dense two-dimensional matrix stored in column-major order.
//!
//! Port of `src/openms/include/OpenMS/DATASTRUCTURES/Matrix.h` (core SDK
//! `bc9cc12`; `Matrix.cpp` only instantiates it for `int`, `double` and
//! `float`). See `docs/MATRIX_SUPPORT.md`.
//!
//! Elements are stored column after column (Fortran order, as Eigen and NumPy's
//! `order='F'` expect): element `(row, col)` lives at `col * rows + row`, so
//! [`Matrix::iter`] visits `(0,0), (1,0), …, (rows-1,0), (0,1), …`.
//!
//! # Example
//!
//! ```
//! use openms::data_structures::matrix::Matrix;
//!
//! let mut m = Matrix::filled(2, 3, 6_i32);
//! m[(1, 2)] = 112;
//! m[(0, 0)] = 100;
//! assert_eq!(m.as_slice(), &[100, 6, 6, 6, 6, 112]);
//! assert_eq!(m.to_string(), "   100      6      6 \n     6      6    112 \n");
//! ```

use std::fmt;
use std::ops::{Index, IndexMut};

use crate::error::{Error, Result};

/// A dense matrix with column-major storage, as the source's
/// `Matrix<Value>`.
///
/// | C++ | Rust |
/// |---|---|
/// | `Matrix()` | [`Matrix::new`], [`Default`] |
/// | `Matrix(rows, cols, value = Value())` | [`Matrix::filled`]; [`Matrix::with_size`] for `Value()` |
/// | copy/move construction and assignment, destructor | [`Clone`], Rust moves, [`Drop`] of the owned buffer |
/// | `rows()`, `cols()`, `size()`, `empty()` | [`Matrix::rows`], [`Matrix::cols`], [`Matrix::size`], [`Matrix::is_empty`] |
/// | `getValue(row, col)` (const and mutable) | [`Matrix::get`], [`Matrix::get_mut`] (checked); [`Index`]/[`IndexMut`] |
/// | `setValue(row, col, value)` | [`Matrix::set_value`] (checked) |
/// | `operator()(row, col)` (const and mutable) | `m[(row, col)]` via [`Index`]/[`IndexMut`] |
/// | `data()` (const and mutable) | [`Matrix::as_slice`], [`Matrix::as_mut_slice`] |
/// | `innerStride()`, `outerStride()`, `rowMajor()` | [`Matrix::inner_stride`], [`Matrix::outer_stride`], [`Matrix::ROW_MAJOR`] |
/// | `resize(rows, cols)` | [`Matrix::resize`] |
/// | `fill(value)`, `clear()` | [`Matrix::fill`], [`Matrix::clear`] |
/// | `setMatrix<T, ROWS, COLS>(array)` | [`Matrix::set_matrix`] |
/// | `begin`/`end`/`cbegin`/`cend` | [`Matrix::iter`], [`Matrix::iter_mut`], `IntoIterator for &Matrix` / `&mut Matrix` |
/// | `maxValue()`, `minValue()` | [`Matrix::max_value`], [`Matrix::min_value`] |
/// | `operator==`, `operator!=` | [`PartialEq`]; [`Matrix::checked_eq`] for the precondition |
/// | `operator<<(std::ostream&, const Matrix&)` | [`fmt::Display`] |
#[derive(Clone, Debug, Default)]
pub struct Matrix<T> {
    data: Vec<T>,
    rows: usize,
    cols: usize,
}

impl<T> Matrix<T> {
    /// Whether the storage is row-major: never, as the source's
    /// `rowMajor()`.
    pub const ROW_MAJOR: bool = false;

    /// An empty 0 × 0 matrix, as the source's default constructor.
    pub fn new() -> Self {
        Self {
            data: Vec::new(),
            rows: 0,
            cols: 0,
        }
    }

    /// Number of rows.
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Number of columns.
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Number of stored elements, the source's `size()`.
    pub fn size(&self) -> usize {
        self.data.len()
    }

    /// Whether no element is stored, the source's `empty()`.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Distance between two consecutive elements of one column: always 1.
    pub fn inner_stride(&self) -> usize {
        1
    }

    /// Distance between the starts of two consecutive columns: the row count.
    ///
    /// The source returns it as `int` via `static_cast<int>(rows_)`; this port
    /// keeps the full `usize`.
    pub fn outer_stride(&self) -> usize {
        self.rows
    }

    fn offset(&self, row: usize, col: usize) -> Option<usize> {
        if row < self.rows && col < self.cols {
            // Both factors are bounded by the storage length, which fits usize.
            Some(col * self.rows + row)
        } else {
            None
        }
    }

    /// The element at `(row, col)`, or `None` outside the matrix.
    ///
    /// The source's `getValue` does not check; an out-of-range column reads
    /// past the buffer and an out-of-range row aliases the next column.
    pub fn get(&self, row: usize, col: usize) -> Option<&T> {
        self.offset(row, col).and_then(|i| self.data.get(i))
    }

    /// The mutable element at `(row, col)`, or `None` outside the matrix.
    pub fn get_mut(&mut self, row: usize, col: usize) -> Option<&mut T> {
        match self.offset(row, col) {
            Some(i) => self.data.get_mut(i),
            None => None,
        }
    }

    /// Store `value` at `(row, col)`, the source's `setValue`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when `(row, col)` lies outside the matrix (the
    /// source does not check).
    pub fn set_value(&mut self, row: usize, col: usize, value: T) -> Result<()> {
        let (rows, cols) = (self.rows, self.cols);
        let slot = self.get_mut(row, col).ok_or_else(|| {
            Error::InvalidRange(format!(
                "element ({row}, {col}) is outside a {rows} x {cols} matrix"
            ))
        })?;
        *slot = value;
        Ok(())
    }

    /// The column-major storage, the source's `data()`.
    pub fn as_slice(&self) -> &[T] {
        &self.data
    }

    /// The mutable column-major storage, the source's non-const `data()`.
    pub fn as_mut_slice(&mut self) -> &mut [T] {
        &mut self.data
    }

    /// Iterate the elements in storage (column-major) order, the source's
    /// `begin()`/`end()` and `cbegin()`/`cend()`.
    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.data.iter()
    }

    /// Iterate the elements mutably in storage order.
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, T> {
        self.data.iter_mut()
    }

    /// Empty the matrix and set both dimensions to zero.
    pub fn clear(&mut self) {
        self.data.clear();
        self.rows = 0;
        self.cols = 0;
    }
}

impl<T: Clone> Matrix<T> {
    /// A `rows` × `cols` matrix with every element set to `value`, the
    /// source's `Matrix(rows, cols, value)`.
    ///
    /// # Panics
    ///
    /// When `rows * cols` overflows `usize` (the source wraps and allocates the
    /// wrapped size); use [`Matrix::try_filled`] for untrusted dimensions.
    pub fn filled(rows: usize, cols: usize, value: T) -> Self {
        match Self::try_filled(rows, cols, value) {
            Ok(matrix) => matrix,
            Err(error) => panic!("{error}"),
        }
    }

    /// As [`Matrix::filled`], refusing dimensions whose product overflows.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when `rows * cols` overflows `usize`.
    pub fn try_filled(rows: usize, cols: usize, value: T) -> Result<Self> {
        let len = element_count(rows, cols)?;
        Ok(Self {
            data: vec![value; len],
            rows,
            cols,
        })
    }

    /// Set every element to `value`.
    pub fn fill(&mut self, value: T) {
        self.data.fill(value);
    }
}

impl<T: Clone + Default> Matrix<T> {
    /// A `rows` × `cols` matrix of `T::default()`, the source's
    /// `Matrix(rows, cols)` with `Value()`.
    ///
    /// # Panics
    ///
    /// When `rows * cols` overflows `usize`.
    pub fn with_size(rows: usize, cols: usize) -> Self {
        Self::filled(rows, cols, T::default())
    }

    /// Set the dimensions to `rows` × `cols`, the source's `resize`.
    ///
    /// As in the source this resizes the flat buffer only: existing elements
    /// keep their storage positions (so they move to other `(row, col)`
    /// positions when the row count changes), the buffer is truncated or
    /// extended with `T::default()`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when `rows * cols` overflows `usize`; the
    /// matrix is unchanged.
    pub fn resize(&mut self, rows: usize, cols: usize) -> Result<()> {
        let len = element_count(rows, cols)?;
        self.rows = rows;
        self.cols = cols;
        self.data.resize(len, T::default());
        Ok(())
    }

    /// Resize to `ROWS` × `COLS` and copy `array[i][j]` to `(i, j)`, the
    /// source's `setMatrix<T, ROWS, COLS>`.
    pub fn set_matrix<const ROWS: usize, const COLS: usize>(&mut self, array: &[[T; COLS]; ROWS]) {
        // ROWS * COLS elements exist in `array`, so the product cannot overflow.
        self.rows = ROWS;
        self.cols = COLS;
        self.data.clear();
        self.data.reserve(ROWS * COLS);
        for col in 0..COLS {
            for row in array.iter() {
                self.data.push(row[col].clone());
            }
        }
    }
}

impl<T: PartialOrd + Clone + Default> Matrix<T> {
    /// The largest element, or `T::default()` for an empty matrix, the
    /// source's `maxValue()`.
    ///
    /// Reproduces `std::max_element`: the first element wins ties, and an
    /// element is taken only when the current best compares `<` it, so a NaN
    /// is never taken after the first position and a leading NaN is kept.
    pub fn max_value(&self) -> T {
        let mut iter = self.data.iter();
        let Some(mut best) = iter.next() else {
            return T::default();
        };
        for value in iter {
            if best < value {
                best = value;
            }
        }
        best.clone()
    }

    /// The smallest element, or `T::default()` for an empty matrix, the
    /// source's `minValue()` (`std::min_element`: first wins ties, an element
    /// is taken only when it compares `<` the current best).
    pub fn min_value(&self) -> T {
        let mut iter = self.data.iter();
        let Some(mut best) = iter.next() else {
            return T::default();
        };
        for value in iter {
            if value < best {
                best = value;
            }
        }
        best.clone()
    }
}

impl<T: PartialEq> Matrix<T> {
    /// The source's `operator==` with its precondition made a checked error.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the dimensions differ; the source raises
    /// `Exception::Precondition` there in a debug build and, in a release
    /// build, compares the flat buffers only (so a 4 × 5 and a 5 × 4 matrix of
    /// equal elements compare equal).
    pub fn checked_eq(&self, rhs: &Self) -> Result<bool> {
        if self.rows != rhs.rows || self.cols != rhs.cols {
            return Err(Error::InvalidValue(format!(
                "matrices must have the same dimensions for comparison ({} x {} vs {} x {})",
                self.rows, self.cols, rhs.rows, rhs.cols
            )));
        }
        Ok(self.data == rhs.data)
    }
}

/// Equal dimensions and equal elements. Matrices of different dimensions are
/// unequal, where the source's release build compares only the flat buffers
/// and its debug build raises a precondition violation (see
/// [`Matrix::checked_eq`]).
impl<T: PartialEq> PartialEq for Matrix<T> {
    fn eq(&self, rhs: &Self) -> bool {
        self.rows == rhs.rows && self.cols == rhs.cols && self.data == rhs.data
    }
}

impl<T: Eq> Eq for Matrix<T> {}

fn element_count(rows: usize, cols: usize) -> Result<usize> {
    rows.checked_mul(cols).ok_or_else(|| {
        Error::InvalidRange(format!("a {rows} x {cols} matrix has too many elements"))
    })
}

impl<T> Index<(usize, usize)> for Matrix<T> {
    type Output = T;

    /// The source's `operator()(row, col)`.
    ///
    /// # Panics
    ///
    /// When `(row, col)` lies outside the matrix; use [`Matrix::get`] for a
    /// checked read.
    fn index(&self, (row, col): (usize, usize)) -> &T {
        let (rows, cols) = (self.rows, self.cols);
        match self.get(row, col) {
            Some(value) => value,
            None => panic!("element ({row}, {col}) is outside a {rows} x {cols} matrix"),
        }
    }
}

impl<T> IndexMut<(usize, usize)> for Matrix<T> {
    /// The source's mutable `operator()(row, col)`.
    ///
    /// # Panics
    ///
    /// When `(row, col)` lies outside the matrix; use [`Matrix::get_mut`].
    fn index_mut(&mut self, (row, col): (usize, usize)) -> &mut T {
        let (rows, cols) = (self.rows, self.cols);
        match self.get_mut(row, col) {
            Some(value) => value,
            None => panic!("element ({row}, {col}) is outside a {rows} x {cols} matrix"),
        }
    }
}

impl<'a, T> IntoIterator for &'a Matrix<T> {
    type Item = &'a T;
    type IntoIter = std::slice::Iter<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.iter()
    }
}

impl<'a, T> IntoIterator for &'a mut Matrix<T> {
    type Item = &'a mut T;
    type IntoIter = std::slice::IterMut<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.iter_mut()
    }
}

/// An element type the source's stream operator can print: its text under
/// `std::setprecision(6)` before `std::setw(6)` pads it.
pub trait MatrixCell {
    /// The element as `os << std::setprecision(6) << value` writes it.
    fn stream_text(&self) -> String;
}

macro_rules! integer_cell {
    ($($t:ty),*) => {$(
        impl MatrixCell for $t {
            fn stream_text(&self) -> String {
                self.to_string()
            }
        }
    )*};
}
integer_cell!(
    i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize
);

impl MatrixCell for f64 {
    fn stream_text(&self) -> String {
        crate::concept::fuzzy_string_comparator::format_g(*self)
    }
}

impl MatrixCell for f32 {
    fn stream_text(&self) -> String {
        // `ostream << float` promotes to double for `%g` formatting.
        crate::concept::fuzzy_string_comparator::format_g(f64::from(*self))
    }
}

/// The source's `operator<<`: one line per row, each element right-aligned
/// in a field of six (`std::setw(6)`) after `std::setprecision(6)`, followed
/// by one space, every line ending in `'\n'`.
impl<T: MatrixCell> fmt::Display for Matrix<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for row in 0..self.rows {
            for col in 0..self.cols {
                if let Some(value) = self.get(row, col) {
                    write!(f, "{:>6} ", value.stream_text())?;
                }
            }
            writeln!(f)?;
        }
        Ok(())
    }
}
