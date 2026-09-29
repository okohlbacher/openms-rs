# Dense matrix

`data_structures::matrix` covers `DATASTRUCTURES/Matrix.h` and `DATASTRUCTURES/Matrix.cpp`
at SDK `bc9cc12`. `Matrix.cpp` only instantiates the template for `int`, `double` and
`float`; the Rust type is generic, and its stream output is implemented for every primitive
integer, `f32` and `f64`.

Tests: [`tests/matrix.rs`](../tests/matrix.rs), every `START_SECTION` of `Matrix_test.cpp`
with unchanged literals, and the members the class test does not reach (`data`, the strides,
`setValue`, `clear`, `minValue`/`maxValue`, `resize` shrinking, floating-point output), each
checked against a value that follows from the header's definition.

`MatrixEigen.h` (the internal zero-copy Eigen view) is a separate header and is not covered
here; `as_slice` with `outer_stride` is the same column-major view.

## API mapping

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `Matrix()` | `Matrix::new()`, `Default` | 0 × 0. |
| `Matrix(rows, cols, value = Value())` | `Matrix::filled(rows, cols, value)`, `Matrix::with_size(rows, cols)`, `Matrix::try_filled` | `filled` panics when `rows * cols` overflows (the source wraps); `try_filled` returns `InvalidRange`. |
| copy/move constructors and assignments, `~Matrix()` | `Clone`, moves, `Drop` | |
| `rows()`, `cols()`, `size()`, `empty()` | `rows()`, `cols()`, `size()`, `is_empty()` | |
| `getValue(row, col)` const / mutable | `get(row, col) -> Option<&T>`, `get_mut(row, col) -> Option<&mut T>` | Checked; the source reads out of bounds and lets an out-of-range row alias the next column. |
| `setValue(row, col, value)` | `set_value(row, col, value) -> Result<()>` | `InvalidRange` outside the matrix. |
| `operator()(row, col)` const / mutable | `m[(row, col)]` (`Index`/`IndexMut`) | Panics outside the matrix, like slice indexing. |
| `data()` const / mutable | `as_slice()`, `as_mut_slice()` | Column-major. |
| `innerStride()`, `outerStride()`, `rowMajor()` | `inner_stride()`, `outer_stride()`, `Matrix::ROW_MAJOR` | `outer_stride` is `usize`, not the source's `static_cast<int>`. |
| `resize(rows, cols)` | `resize(rows, cols) -> Result<()>` | Resizes the flat buffer exactly as `std::vector::resize`: elements keep their storage positions, new slots are `T::default()`. `InvalidRange` on overflow, matrix unchanged. |
| `fill(value)`, `clear()` | `fill(value)`, `clear()` | |
| `setMatrix<T, ROWS, COLS>(const T (&)[ROWS][COLS])` | `set_matrix(&[[T; COLS]; ROWS])` | Const generics. |
| `begin`/`end`, `cbegin`/`cend` | `iter()`, `iter_mut()`, `IntoIterator for &Matrix` and `&mut Matrix` | Storage order. |
| `maxValue()`, `minValue()` | `max_value()`, `min_value()` | `T::default()` when empty; `std::max_element`/`min_element` semantics, including NaN handling (a leading NaN is kept, a later one never taken). |
| `operator==`, `operator!=` | `PartialEq`, `Eq`; `checked_eq(&rhs) -> Result<bool>` | See below. |
| `operator<<(std::ostream&, const Matrix&)` | `Display` | `setprecision(6)` then `setw(6)` per element, one trailing space per element, `'\n'` per row; floats via `%g` (`concept::fuzzy_string_comparator::format_g`). |

## Divergence

`operator==` carries `OPENMS_PRECONDITION(rows_ == rhs.rows_ && cols_ == rhs.cols_)`, which
throws `Exception::Precondition` in a debug build and is compiled out of a release build, where
only the flat buffers are compared: a 4 × 5 and a 5 × 4 matrix of equal elements compare equal.
`PartialEq` compares the dimensions too, so matrices of different shape are unequal;
`checked_eq` reproduces the precondition as an `InvalidValue` error. The class test's three
`TEST_PRECONDITION_VIOLATED` cases are checked against `checked_eq`.

## Evidence

Tier 3: `Matrix_test.cpp` transcribed with unchanged literals, including the whitespace of the
stream-output dump. No C++ was executed.
