// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Port of `DPosition_test.cpp` (core SDK bc9cc12): one test per
//! `START_SECTION`, cited by source line. Literals are transcribed (tier 3).

use openms::data_structures::{DPosition, DPosition1, DPosition2};
use std::cmp::Ordering;
use std::hash::{DefaultHasher, Hash, Hasher};

fn hash(value: &impl Hash) -> u64 {
    let mut h = DefaultHasher::new();
    value.hash(&mut h);
    h.finish()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1.),
        "{actual:.17} != {expected:.17}"
    );
}

// L38 DPosition()
#[test]
fn default_constructor_is_all_zero() {
    let d10 = DPosition::<10>::default();
    assert_eq!(d10.coordinates, [0.; 10]);
    assert_eq!(DPosition::<10>::new(), d10);
}

// L43 ~DPosition()
#[test]
fn destructor_is_ordinary_drop() {
    let d10 = Box::new(DPosition::<10>::new());
    assert_eq!(d10[9], 0.);
    drop(d10);
}

// L47 swap
#[test]
fn swap_exchanges_all_coordinates() {
    let mut i = DPosition::xyz(1., 2., 3.);
    let mut j = DPosition::xyz(4., 5., 6.);
    std::mem::swap(&mut i, &mut j);
    assert_eq!(i.coordinates, [4., 5., 6.]);
    assert_eq!(j.coordinates, [1., 2., 3.]);
}

// L61 abs. The Int64 half of the section exercises the integer std::abs
// overload; the coordinate type is fixed to f64 here, so only the double
// half is ported, with the same exact (not similar) comparisons.
#[test]
fn abs_uses_double_precision() {
    let small_negative_double = -f64::EPSILON;
    let j = DPosition::xyz(-1.4444, -small_negative_double, small_negative_double).abs();
    assert_eq!(j[0], 1.4444);
    assert_eq!(j[1], -small_negative_double);
    assert_eq!(j[2], -small_negative_double);
    assert_eq!(DPosition::xy(-0., -3.).abs().coordinates, [0., 3.]);
}

// L83 operator[] const
#[test]
fn const_index_reads_zero_and_get_is_checked() {
    let i = DPosition::<3>::new();
    assert_eq!(i[0], 0.);
    assert_eq!(i[1], 0.);
    assert_eq!(i[2], 0.);
    assert_eq!(i.get(2), Some(0.));
    assert_eq!(i.get(3), None);
}

// L83 TEST_PRECONDITION_VIOLATED(i[3])
#[test]
#[should_panic]
fn const_index_out_of_range_panics() {
    let i = DPosition::<3>::new();
    let _ = i[3];
}

// L91 operator[] mutable
#[test]
fn mutable_index_writes_each_coordinate() {
    let mut i = DPosition::<3>::new();
    i[0] = 1.;
    assert_eq!(i.coordinates, [1., 0., 0.]);
    i[1] = 2.;
    assert_eq!(i.coordinates, [1., 2., 0.]);
    i[2] = 3.;
    assert_eq!(i.coordinates, [1., 2., 3.]);
    assert_eq!(i.get_mut(3), None);
    *i.get_mut(0).unwrap() = 7.;
    assert_eq!(i[0], 7.);
}

// L91 TEST_PRECONDITION_VIOLATED(i[3] = 4.0)
#[test]
#[should_panic]
fn mutable_index_out_of_range_panics() {
    let mut i = DPosition::<3>::new();
    i[3] = 4.;
}

// L109 copy constructor
#[test]
fn copy_constructor_copies_all_coordinates() {
    let mut p = DPosition::<3>::new();
    p[0] = 12.3;
    p[1] = 23.4;
    p[2] = 34.5;
    let copy_of_p = p;
    assert_eq!(copy_of_p[0], p[0]);
    assert_eq!(copy_of_p[1], p[1]);
    assert_eq!(copy_of_p[2], p[2]);
    assert_eq!(DPosition::<3>::size(), DPosition::<3>::size());
    assert_eq!(copy_of_p, p);
}

// L121 operator=
#[test]
fn assignment_copies_all_coordinates() {
    let mut p = DPosition::<3>::new();
    p[0] = 12.3;
    p[1] = 23.4;
    p[2] = 34.5;
    let mut copy_of_p = DPosition::<3>::filled(9.);
    assert_ne!(copy_of_p, p);
    copy_of_p = p;
    assert_eq!(copy_of_p.coordinates, [12.3, 23.4, 34.5]);
    assert_eq!(copy_of_p[0], p[0]);
}

// L134 DPosition(CoordinateType x)
#[test]
fn filled_constructor_sets_all_dimensions() {
    let p = DPosition::<3>::filled(12.34);
    close(p[0], 12.34);
    close(p[1], 12.34);
    close(p[2], 12.34);
}

// L141 DPosition(x, y)
#[test]
fn xy_constructor() {
    let p = DPosition::xy(1., 2.);
    close(p[0], 1.);
    close(p[1], 2.);
}

// L147 DPosition(x, y, z)
#[test]
fn xyz_constructor() {
    let p = DPosition::xyz(1., 2., 3.);
    close(p[0], 1.);
    close(p[1], 2.);
    close(p[2], 3.);
}

// L154 inner product operator*
#[test]
fn inner_product() {
    let i = DPosition::xyz(2., 3., 4.);
    let j = DPosition::xyz(3., 4., 5.);
    close(i * j, 6. + 12. + 20.);
    close(i.dot(&j), 38.);
}

fn one_to_ten() -> DPosition<10> {
    DPosition::from([1., 2., 3., 4., 5., 6., 7., 8., 9., 10.])
}

// L177 const begin()
#[test]
fn const_begin_yields_first_coordinate() {
    let i = one_to_ten();
    assert_eq!(i.iter().next(), Some(&1.));
    assert_eq!((&i).into_iter().next(), Some(&1.));
}

// L182 const end(): copying begin..end yields all ten values
#[test]
fn const_iteration_copies_all_ten_values() {
    let i = one_to_ten();
    let v: Vec<f64> = i.iter().copied().collect();
    assert_eq!(v.len(), 10);
    for (k, value) in v.iter().enumerate() {
        close(*value, (k + 1) as f64);
    }
    assert_eq!(i.as_slice(), &v[..]);
    assert_eq!(i.into_iter().collect::<Vec<_>>(), v);
}

// L201 mutable begin(): *begin = 11
#[test]
fn mutable_begin_writes_first_coordinate() {
    let mut i = one_to_ten();
    assert_eq!(*i.iter_mut().next().unwrap(), 1.);
    *i.iter_mut().next().unwrap() = 11.;
    assert_eq!(i[0], 11.);
    for value in &mut i {
        *value += 0.;
    }
    assert_eq!(i[0], 11.);
}

// L208 mutable end(): copy after the modification
#[test]
fn mutable_iteration_sees_the_modification() {
    let mut i = one_to_ten();
    *i.iter_mut().next().unwrap() = 11.;
    let v: Vec<f64> = i.iter().copied().collect();
    assert_eq!(v.len(), 10);
    close(v[0], 11.);
    close(v[1], 2.);
    close(v[9], 10.);
}

// L227 static size()
#[test]
fn size_is_the_dimension() {
    assert_eq!(DPosition::<777>::size(), 777);
    assert_eq!(DPosition::<3>::size(), 3);
    assert_eq!(DPosition::<1>::size(), 1);
    assert_eq!(DPosition::<123>::size(), 123);
    assert_eq!(DPosition::<123>::DIMENSION, 123);
}

// L237 clear()
#[test]
fn clear_sets_all_dimensions_to_zero() {
    let mut p = DPosition::xyz(1.2, 2.3, 3.4);
    close(p[0], 1.2);
    close(p[1], 2.3);
    close(p[2], 3.4);
    p.clear();
    assert_eq!(p.coordinates, [0., 0., 0.]);
}

// L251 operator==
#[test]
fn equality_compares_every_dimension() {
    let mut p1 = DPosition::<3>::new();
    let mut p2 = DPosition::<3>::new();
    assert!(p1 == p2);
    p1[0] = 1.234;
    assert!(p1 != p2);
    p2[0] = 1.234;
    assert!(p1 == p2);
    p1[1] = 1.345;
    assert!(p1 != p2);
    p2[1] = 1.345;
    assert!(p1 == p2);
    p1[2] = 1.456;
    assert!(p1 != p2);
    p2[2] = 1.456;
    assert!(p1 == p2);
}

// L271 operator!=
#[test]
fn inequality_is_the_negation_of_equality() {
    let mut p1 = DPosition::<3>::new();
    let mut p2 = DPosition::<3>::new();
    assert!(!(p1 != p2));
    p1[0] = 1.234;
    assert!(!(p1 == p2));
    p2[0] = 1.234;
    assert!(!(p1 != p2));
    p1[1] = 1.345;
    assert!(!(p1 == p2));
    p2[1] = 1.345;
    assert!(!(p1 != p2));
    p1[2] = 1.456;
    assert!(!(p1 == p2));
    p2[2] = 1.456;
    assert!(!(p1 != p2));
}

// L291 operator< is lexicographic over dimensions 0..D-1
#[test]
fn less_than_is_lexicographic() {
    let mut p1 = DPosition::<3>::new();
    let mut p2 = DPosition::<3>::new();
    assert_ne!(p1.partial_cmp(&p2), Some(Ordering::Less));
    for dim in 0..3 {
        p1[dim] = p2[dim] - 0.1;
        assert!(p1 < p2);
        p2[dim] = p1[dim] - 0.1;
        assert_ne!(p1.partial_cmp(&p2), Some(Ordering::Less));
        p2[dim] = p1[dim];
    }
    // A later dimension does not outweigh an earlier one.
    assert!(DPosition::xy(0., 100.) < DPosition::xy(1., 0.));
    assert_eq!(
        DPosition::xy(f64::NAN, 0.).partial_cmp(&DPosition::xy(1., 0.)),
        None
    );
}

// L314 operator>
#[test]
fn greater_than_is_lexicographic() {
    let mut p1 = DPosition::<3>::new();
    let mut p2 = DPosition::<3>::new();
    assert_ne!(p1.partial_cmp(&p2), Some(Ordering::Greater));
    p1[0] = p2[0] - 0.1;
    assert_ne!(p1.partial_cmp(&p2), Some(Ordering::Greater));
    p2[0] = p1[0] - 0.1;
    assert!(p1 > p2);
}

// L325 operator>=
#[test]
fn greater_equal_is_lexicographic() {
    let mut p1 = DPosition::<3>::new();
    let mut p2 = DPosition::<3>::new();
    assert!(p1 >= p2);
    p1[0] = p2[0] - 0.1;
    assert_eq!(p1.partial_cmp(&p2), Some(Ordering::Less));
    p2[0] = p1[0] - 0.1;
    assert!(p1 >= p2);
}

// L336 operator<=
#[test]
fn less_equal_is_lexicographic() {
    let mut p1 = DPosition::<3>::new();
    let mut p2 = DPosition::<3>::new();
    assert!(p1 <= p2);
    p1[0] = p2[0] - 0.1;
    assert!(p1 <= p2);
    p2[0] = p1[0] - 0.1;
    assert_eq!(p1.partial_cmp(&p2), Some(Ordering::Greater));
}

// L346 unary operator-
#[test]
fn negation_flips_signs_and_is_an_involution() {
    let mut p1 = DPosition::<3>::new();
    p1[0] = 5.;
    let mut p2 = -p1;
    assert!(p1 != p2);
    assert_eq!(p2[0], -5.);
    p2 = -p2;
    assert!(p1 == p2);
}

// L355 binary operator-
#[test]
fn subtraction_is_component_wise() {
    let p1 = DPosition::xyz(1.234, 2.234, 3.234);
    let p2 = DPosition::xyz(0.234, 0.234, 0.234);
    let p3 = DPosition::xyz(1., 2., 3.);
    let d = p1 - p2;
    close(d[0], p3[0]);
    close(d[1], p3[1]);
    let e = p2 - p1;
    close(e[0], -p3[0]);
    close(e[1], -p3[1]);
}

// L372 binary operator+
#[test]
fn addition_is_component_wise() {
    let p1 = DPosition::xyz(-1., -2., -3.);
    let p2 = DPosition::xyz(1., 2., 3.);
    let p3 = DPosition::<3>::new();
    assert!((p1 + p2) == p3);
}

// L383 DPosition(x, y) with float literals
#[test]
fn xy_constructor_from_float_literals() {
    let p1 = DPosition::xy(f64::from(11.0f32), f64::from(12.1f32));
    close(p1[0], 11.);
    close(p1[1], f64::from(12.1f32));
    let p = DPosition::xy(12.34, 56.78);
    close(p[0], 12.34);
    close(p[1], 56.78);
}

// L392 getX
#[test]
fn x_reads_the_first_dimension() {
    let p1 = DPosition::xy(11., f64::from(12.1f32));
    close(p1.x(), 11.);
}

// L397 getY
#[test]
fn y_reads_the_second_dimension() {
    let p1 = DPosition::xy(11., f64::from(12.1f32));
    close(p1.y(), f64::from(12.1f32));
}

// L402 setX
#[test]
fn set_x_writes_only_the_first_dimension() {
    let mut p1 = DPosition::xy(11., f64::from(12.1f32));
    p1.set_x(5.);
    close(p1[0], 5.);
    close(p1[1], f64::from(12.1f32));
}

// L409 setY
#[test]
fn set_y_writes_only_the_second_dimension() {
    let mut p1 = DPosition::xy(11., f64::from(12.1f32));
    p1.set_y(5.);
    close(p1[0], 11.);
    close(p1[1], 5.);
}

// L416 operator*=
#[test]
fn scalar_multiply_assign() {
    let mut p1 = DPosition::xy(3., 4.);
    p1 *= 5.;
    let p2 = DPosition::xy(15., 20.);
    close(p1[0], p2[0]);
    close(p1[1], p2[1]);
}

// L424 operator+=
#[test]
fn add_assign() {
    let mut p1 = DPosition::xy(3., 4.);
    let p2 = DPosition::xy(15., 20.);
    p1 += p2;
    let p3 = DPosition::xy(18., 24.);
    close(p1[0], p3[0]);
    close(p1[1], p3[1]);
}

// L433 operator-=
#[test]
fn sub_assign() {
    let mut p1 = DPosition::xy(3., 4.);
    let p2 = DPosition::xy(18., 24.);
    p1 -= p2;
    let p3 = DPosition::xy(-15., -20.);
    close(p1[0], p3[0]);
    close(p1[1], p3[1]);
}

// L442 operator/=
#[test]
fn scalar_divide_assign() {
    let mut p1 = DPosition::xy(15., 20.);
    p1 /= 5.;
    let p2 = DPosition::xy(3., 4.);
    close(p1[0], p2[0]);
    close(p1[1], p2[1]);
}

fn unit_square() -> [DPosition2; 4] {
    [
        DPosition::xy(0., 0.),
        DPosition::xy(0., 1.),
        DPosition::xy(1., 0.),
        DPosition::xy(1., 1.),
    ]
}

// L450 spatiallyGreaterEqual: 16 literal combinations
#[test]
fn spatially_greater_equal_truth_table() {
    let [p00, p01, p10, p11] = unit_square();
    let expected = [
        [true, false, false, false],
        [true, true, false, false],
        [true, false, true, false],
        [true, true, true, true],
    ];
    for (row, a) in [p00, p01, p10, p11].iter().enumerate() {
        for (col, b) in [p00, p01, p10, p11].iter().enumerate() {
            assert_eq!(
                a.spatially_greater_equal(b),
                expected[row][col],
                "{row},{col}"
            );
        }
    }
}

// L473 spatiallyLessEqual: 16 literal combinations
#[test]
fn spatially_less_equal_truth_table() {
    let [p00, p01, p10, p11] = unit_square();
    let expected = [
        [true, true, true, true],
        [false, true, false, true],
        [false, false, true, true],
        [false, false, false, true],
    ];
    for (row, a) in [p00, p01, p10, p11].iter().enumerate() {
        for (col, b) in [p00, p01, p10, p11].iter().enumerate() {
            assert_eq!(a.spatially_less_equal(b), expected[row][col], "{row},{col}");
        }
    }
}

// L496 zero()
#[test]
fn zero_static() {
    assert_eq!(DPosition1::zero()[0], 0.);
    assert_eq!(DPosition::<4>::zero().coordinates, [0.; 4]);
}

// L501 minPositive() == numeric_limits<double>::min()
#[test]
fn min_positive_static() {
    assert_eq!(DPosition1::min_positive()[0], f64::MIN_POSITIVE);
    assert_eq!(DPosition1::min_positive()[0], 2.2250738585072014e-308);
}

// L506 minNegative() == -numeric_limits<double>::max()
#[test]
fn min_negative_static() {
    assert_eq!(DPosition1::min_negative()[0], -f64::MAX);
    assert_eq!(DPosition1::min_negative()[0], f64::MIN);
    assert!(DPosition1::min_negative()[0].is_finite());
}

// L511 maxPositive() == numeric_limits<double>::max()
#[test]
fn max_positive_static() {
    assert_eq!(DPosition1::max_positive()[0], f64::MAX);
    assert!(DPosition1::max_positive()[0].is_finite());
}

// L516 [EXTRA] int DPosition: the coordinate type is fixed to f64, so the
// integer instantiation is exercised with the same integral values.
#[test]
fn integral_values_spatially_greater_equal() {
    let [p00, p01, p10, p11] = unit_square();
    assert!(p00.spatially_greater_equal(&p00));
    assert!(!p00.spatially_greater_equal(&p01));
    assert!(!p00.spatially_greater_equal(&p10));
    assert!(!p00.spatially_greater_equal(&p11));
    assert!(p01.spatially_greater_equal(&p00));
    assert!(p01.spatially_greater_equal(&p01));
    assert!(!p01.spatially_greater_equal(&p10));
    assert!(!p01.spatially_greater_equal(&p11));
    assert!(p10.spatially_greater_equal(&p00));
    assert!(!p10.spatially_greater_equal(&p01));
    assert!(p10.spatially_greater_equal(&p10));
    assert!(!p10.spatially_greater_equal(&p11));
    assert!(p11.spatially_greater_equal(&p00));
    assert!(p11.spatially_greater_equal(&p01));
    assert!(p11.spatially_greater_equal(&p10));
    assert!(p11.spatially_greater_equal(&p11));
}

// L541 [EXTRA] char DPosition: 'a' = 97 and 'b' = 98 as f64.
#[test]
fn character_code_values_negate_and_order() {
    let mut pa1 = DPosition::<3>::new();
    pa1[0] = f64::from(b'a');
    let mut pb2 = -pa1;
    assert!(pa1 != pb2);
    pb2 = -pb2;
    assert!(pa1 == pb2);
    let pa = DPosition1::filled(f64::from(b'a'));
    let pb = DPosition1::filled(f64::from(b'b'));
    assert!(pa < pb);
}

// L558 [EXTRA] scalar multiplication, both operand orders
#[test]
fn scalar_multiplication_both_orders() {
    let p1 = DPosition::xy(3., 4.);
    let p2 = p1 * 5.;
    let p3 = 5. * p1;
    let expected = DPosition::xy(15., 20.);
    close(p2[0], expected[0]);
    close(p2[1], expected[1]);
    close(p3[0], expected[0]);
    close(p3[1], expected[1]);
}

// L574 [EXTRA] scalar division
#[test]
fn scalar_division() {
    let p1 = DPosition::xy(15., 20.);
    let p2 = p1 / 5.;
    let expected = DPosition::xy(3., 4.);
    close(p2[0], expected[0]);
    close(p2[1], expected[1]);
}

// L589 [EXTRA] std::hash<DPosition<2>>. Unordered containers need Eq, which
// the f64 coordinates do not provide; the container assertions are
// represented by equality/hash agreement, as the crate's other value types do.
#[test]
fn hash_of_two_dimensional_positions() {
    let p1 = DPosition::xy(1.5, 2.5);
    let p2 = DPosition::xy(1.5, 2.5);
    let p3 = DPosition::xy(3.5, 2.5);
    assert_eq!(hash(&p1), hash(&p2));
    assert_ne!(hash(&p1), hash(&p3));
    assert_eq!(p1, p2);
    assert_ne!(p1, p3);
    assert_eq!(hash(&DPosition::xy(0., -0.)), hash(&DPosition::xy(-0., 0.)));
    assert_eq!(DPosition::xy(0., -0.), DPosition::xy(-0., 0.));
}

// L622 [EXTRA] std::hash<DPosition<1>>
#[test]
fn hash_of_one_dimensional_positions() {
    let p1 = DPosition1::filled(1.5);
    let p2 = DPosition1::filled(1.5);
    let p3 = DPosition1::filled(2.5);
    assert_eq!(hash(&p1), hash(&p2));
    assert_ne!(hash(&p1), hash(&p3));
    assert_eq!(p1, p2);
}

// Native: stream output is space separated without a newline.
#[test]
fn display_is_space_separated() {
    assert_eq!(DPosition::xyz(1., 2.5, -3.).to_string(), "1 2.5 -3");
    assert_eq!(format!("{:.2}", DPosition::xy(1., 2.555)), "1.00 2.56");
    assert_eq!(DPosition1::filled(0.1).to_string(), "0.1");
    let array: [f64; 2] = DPosition::xy(1., 2.).into();
    assert_eq!(array, [1., 2.]);
}
