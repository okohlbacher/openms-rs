// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `src/tests/class_tests/openms/source/Matrix_test.cpp` (core SDK
//! `bc9cc12`), section by section with unchanged literals, plus the members
//! the class test does not reach, each checked against a value derived from
//! `Matrix.h` itself.

use openms::data_structures::matrix::Matrix;

/// The fixture `mi` the class test builds up across sections: resized to
/// 2 x 3 and filled with 7.
fn fixture_mi() -> Matrix<i32> {
    let mut mi = Matrix::new();
    mi.resize(2, 2).unwrap();
    mi.fill(3);
    mi.resize(2, 3).unwrap();
    mi.fill(7);
    mi
}

// START_SECTION((Matrix()))
#[test]
fn default_constructor() {
    let mi1: Matrix<i32> = Matrix::new();
    assert_eq!(mi1.size(), 0);
    assert_eq!(mi1.cols(), 0);
    assert_eq!(mi1.rows(), 0);
    assert_eq!(mi1.iter().count(), 0);
    assert_eq!((&mi1).into_iter().count(), 0);
    assert_eq!(mi1.to_string(), "");
    assert_eq!(Matrix::<i32>::default(), mi1);
}

// START_SECTION((~Matrix())) is Drop; nothing to test.

// START_SECTION((void resize(size_type i, size_type j)))
#[test]
fn resize() {
    let mi = fixture_mi();
    for (r, c) in [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)] {
        assert_eq!(mi[(r, c)], 7);
    }
}

// START_SECTION((Matrix(const Matrix & source)))
#[test]
fn copy_constructor_and_column_major_iteration() {
    let mi = fixture_mi();
    let mi2 = mi.clone();
    assert_eq!(mi2.cols(), 3);
    assert_eq!(mi2.rows(), 2);
    for (r, c) in [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)] {
        assert_eq!(mi2[(r, c)], 7);
    }
    let (mut row, mut col) = (0, 0);
    for value in mi2.iter() {
        assert_eq!(*value, *mi.get(row, col).unwrap());
        row += 1;
        if row == mi2.rows() {
            row = 0;
            col += 1;
        }
    }
    let (mut row, mut col) = (0, 0);
    for value in &mi2 {
        assert_eq!(*value, *mi.get(row, col).unwrap());
        row += 1;
        if row == mi2.rows() {
            row = 0;
            col += 1;
        }
    }
}

// START_SECTION((Matrix& operator = (const Matrix & rhs)))
#[test]
fn assignment() {
    let mi = fixture_mi();
    let mut mi3: Matrix<i32> = Matrix::new();
    assert_eq!(mi3.size(), 0);
    mi3 = mi.clone();
    assert_eq!(mi3.cols(), 3);
    assert_eq!(mi3.rows(), 2);
    for (r, c) in [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)] {
        assert_eq!(mi3[(r, c)], 7);
    }
}

// The class test's `mi(1,1)=17;` between sections, then the three
// operator() sections and the second assignment section, in order.
#[test]
fn element_access_sections() {
    let mut mi = fixture_mi();
    mi[(1, 1)] = 17;
    // START_SECTION((const_reference operator()(...) const))
    {
        let micr = &mi;
        assert_eq!(micr[(1, 1)], 17);
    }
    // START_SECTION((reference operator()(size_type const i, size_type const j)))
    mi[(1, 2)] = 33;
    {
        let micr = &mi;
        assert_eq!(micr[(1, 2)], 33);
    }
    // START_SECTION((reference operator() (size_type const i, size_type const j)))
    mi[(1, 0)] = 44;
    {
        let micr = &mi;
        assert_eq!(micr[(1, 0)], 44);
    }
    // START_SECTION((void operator()(size_type const i, size_type const j) = value_type value))
    mi[(1, 1)] = 18;
    assert_eq!(mi[(1, 1)], 18);
}

// START_SECTION((Value& getValue(size_t const i, size_t const j)))
#[test]
fn get_value_mutable() {
    let mut test_matrix = Matrix::filled(2, 2, 0_i32);
    *test_matrix.get_mut(0, 0).unwrap() = 100;
    *test_matrix.get_mut(0, 1).unwrap() = 200;
    *test_matrix.get_mut(1, 0).unwrap() = 300;
    *test_matrix.get_mut(1, 1).unwrap() = 400;
    assert_eq!(test_matrix[(0, 0)], 100);
    assert_eq!(test_matrix[(0, 1)], 200);
    assert_eq!(test_matrix[(1, 0)], 300);
    assert_eq!(test_matrix[(1, 1)], 400);
    *test_matrix.get_mut(0, 0).unwrap() += 5;
    assert_eq!(test_matrix[(0, 0)], 105);
}

// START_SECTION((Matrix(const SizeType rows, const SizeType cols, ValueType value = ValueType())))
// START_SECTION((SizeType cols() const)) / ((SizeType rows() const)): the
// source's section titles are swapped relative to their bodies; the bodies
// are kept.
#[test]
fn fill_constructor_rows_cols() {
    let mi5 = Matrix::filled(4, 5, 6_i32);
    assert_eq!(mi5.size(), 20);
    assert_eq!(mi5.rows(), 4);
    assert_eq!(mi5.cols(), 5);
    // `Matrix<float> mf(6,7,8);` is constructed at file scope and unused.
    let mf = Matrix::filled(6, 7, 8.0_f32);
    assert_eq!(mf.size(), 42);
}

// START_SECTION((bool operator==(Matrix const &rhs) const))
#[test]
fn equality() {
    let mut mi1 = Matrix::filled(4, 5, 6_i32);
    mi1[(2, 3)] = 17;
    let mi2 = Matrix::filled(4, 5, 6_i32);
    assert_ne!(mi1, mi2);
    mi1[(2, 3)] = 6;
    assert_eq!(mi1, mi2);
    let mi3 = Matrix::filled(5, 4, 6_i32);
    let mi4 = Matrix::filled(4, 4, 6_i32);
    let mi5 = Matrix::filled(5, 5, 6_i32);
    // TEST_PRECONDITION_VIOLATED for each: the checked comparison refuses.
    assert!(mi1.checked_eq(&mi3).is_err());
    assert!(mi1.checked_eq(&mi4).is_err());
    assert!(mi1.checked_eq(&mi5).is_err());
    assert!(mi1.checked_eq(&mi2).unwrap());
    // Operator form: different dimensions are unequal (documented divergence
    // from the source's release build, which compares the flat buffers).
    assert_ne!(mi1, mi3);
    assert_ne!(mi1, mi4);
    assert_ne!(mi1, mi5);
}

// START_SECTION((template <int ROWS, int COLS> void setMatrix(const ValueType matrix[ROWS][COLS])))
#[test]
fn set_matrix() {
    let test_matrix: [[f64; 4]; 4] = [
        [0.0, 2.5, 3.0, 0.1],
        [0.0, 1.0, 5.9, 0.2],
        [0.0, 2.0, 5.6, 0.1],
        [0.0, 2.0, 3.0, 0.1],
    ];
    let mut my_matrix: Matrix<f64> = Matrix::new();
    my_matrix.set_matrix(&test_matrix);
    for (i, row) in test_matrix.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            assert_eq!(my_matrix[(i, j)], *value);
        }
    }
}

// START_SECTION((template <typename Value> std::ostream & operator<<(...)))
#[test]
fn stream_output() {
    let mut mi = Matrix::filled(2, 3, 6_i32);
    mi[(1, 2)] = 112;
    mi[(0, 0)] = 100;
    mi[(1, 1)] = 111;
    mi[(0, 2)] = 103;
    let matrix_dump = "   100      6    103 \n     6    111    112 \n";
    assert_eq!(mi.to_string(), matrix_dump);
}

// Members the class test does not reach; expected values follow from the
// definitions in Matrix.h.

#[test]
fn storage_layout_and_strides() {
    let mut m = Matrix::filled(2, 3, 0_i32);
    m.set_value(1, 2, 5).unwrap();
    // data_[col * rows_ + row] = data_[2 * 2 + 1]
    assert_eq!(m.as_slice()[5], 5);
    m.as_mut_slice()[0] = 9;
    assert_eq!(m[(0, 0)], 9);
    assert_eq!(m.inner_stride(), 1);
    assert_eq!(m.outer_stride(), 2);
    const { assert!(!Matrix::<i32>::ROW_MAJOR) };
    assert!(m.set_value(2, 0, 1).is_err());
    assert!(m.get(0, 3).is_none());
    assert!(
        m.get(2, 0).is_none(),
        "an out-of-range row must not alias the next column"
    );
    for value in &mut m {
        *value += 1;
    }
    assert_eq!(m.as_slice(), &[10, 1, 1, 1, 1, 6]);
    assert!(!m.is_empty());
    m.clear();
    assert!(m.is_empty());
    assert_eq!((m.rows(), m.cols(), m.size()), (0, 0, 0));
}

#[test]
fn resize_keeps_flat_positions() {
    // std::vector::resize: the buffer keeps its prefix, new slots are Value().
    let mut m: Matrix<i32> = Matrix::new();
    m.set_matrix(&[[1, 2], [3, 4]]); // storage [1, 3, 2, 4]
    m.resize(3, 2).unwrap();
    assert_eq!(m.as_slice(), &[1, 3, 2, 4, 0, 0]);
    assert_eq!(m[(0, 1)], 4);
    m.resize(1, 1).unwrap();
    assert_eq!(m.as_slice(), &[1]);
    assert!(m.resize(usize::MAX, 2).is_err());
    assert_eq!((m.rows(), m.cols()), (1, 1));
    assert!(Matrix::try_filled(usize::MAX, 2, 0_u8).is_err());
    assert_eq!(Matrix::<f64>::with_size(2, 2).as_slice(), &[0.0; 4]);
}

#[test]
fn min_and_max_values() {
    let empty: Matrix<f64> = Matrix::new();
    assert_eq!(empty.max_value(), 0.0);
    assert_eq!(empty.min_value(), 0.0);
    let mut m: Matrix<f64> = Matrix::new();
    m.set_matrix(&[[0.0, 2.5, 3.0, 0.1], [0.0, 1.0, 5.9, 0.2]]);
    assert_eq!(m.max_value(), 5.9);
    assert_eq!(m.min_value(), 0.0);
    // std::max_element / min_element with a NaN: a leading NaN is kept, a
    // later one is never taken (every comparison with it is false).
    let mut leading = Matrix::filled(1, 3, 1.0_f64);
    leading[(0, 0)] = f64::NAN;
    assert!(leading.max_value().is_nan());
    assert!(leading.min_value().is_nan());
    let mut trailing = Matrix::filled(1, 3, 1.0_f64);
    trailing[(0, 2)] = f64::NAN;
    assert_eq!(trailing.max_value(), 1.0);
    assert_eq!(trailing.min_value(), 1.0);
}

#[test]
fn stream_output_of_floating_point() {
    // setprecision(6) is %g with six significant digits; setw(6) pads.
    let mut m: Matrix<f64> = Matrix::new();
    m.set_matrix(&[[0.1, 1234567.0], [-2.5, 1.0 / 3.0]]);
    assert_eq!(m.to_string(), "   0.1 1.23457e+06 \n  -2.5 0.333333 \n");
    let f = Matrix::filled(1, 1, 8.0_f32);
    assert_eq!(f.to_string(), "     8 \n");
}

#[test]
#[should_panic(expected = "outside a 2 x 2 matrix")]
fn index_outside_panics() {
    let m = Matrix::filled(2, 2, 0_i32);
    let _ = m[(2, 0)];
}
