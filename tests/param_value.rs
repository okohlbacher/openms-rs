// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::param::{ParamValue as V, ParamValueType as T};
use std::cmp::Ordering;
use std::hash::{Hash, Hasher};

#[test]
fn source_storage_types_copy_move_assignment_and_empty_distinction() {
    let cases = [
        (V::default(), T::Empty),
        (V::from("test char"), T::String),
        (V::from(String::from("test string")), T::String),
        (V::from(-3000i16), T::Integer),
        (V::from(3000u16), T::Integer),
        (V::from(-3000i32), T::Integer),
        (V::from(3000u32), T::Integer),
        (V::from(-3000i64), T::Integer),
        (V::try_from(3000u64).unwrap(), T::Integer),
        (V::try_from(3000usize).unwrap(), T::Integer),
        (V::try_from(-3000isize).unwrap(), T::Integer),
        (V::from(3.0f32), T::Float),
        (V::from(-3.4f64), T::Float),
        (V::from(vec!["a".to_owned(), "b".to_owned()]), T::StringList),
        (V::from(vec![1i32, 2, 3]), T::IntegerList),
        (V::from(vec![1.2f64, 2.3, 3.4]), T::FloatList),
    ];
    for (value, kind) in cases {
        assert_eq!(value.value_type(), kind);
        let mut moved = value.checked_clone().unwrap();
        assert_eq!(moved, value);
        let saved = std::mem::take(&mut moved);
        assert_eq!(saved, value);
        assert_eq!(moved, V::EMPTY);
        moved = saved;
        assert_eq!(moved, value);
    }
    assert!(V::EMPTY.is_empty());
    for nonempty in [
        V::from(""),
        V::StringList(vec![]),
        V::IntegerList(vec![]),
        V::FloatList(vec![]),
    ] {
        assert!(!nonempty.is_empty());
    }
    assert_eq!(
        V::from(1.23f32).to_f64().unwrap().to_bits(),
        f64::from(1.23f32).to_bits()
    );
}

#[test]
fn source_typed_conversions_and_boolean_rules() {
    assert_eq!(V::from(-55).to_i16().unwrap(), -55);
    assert_eq!(V::from(-55).to_i32().unwrap(), -55);
    assert_eq!(V::from(-55).to_i64().unwrap(), -55);
    assert_eq!(V::from(55).to_u16().unwrap(), 55);
    assert_eq!(V::from(55).to_u32().unwrap(), 55);
    assert_eq!(V::from(55).to_u64().unwrap(), 55);
    assert_eq!(V::from(55).to_f32().unwrap(), 55.0);
    assert_eq!(V::from(55).to_f64().unwrap(), 55.0);
    assert_eq!(V::from(5.4).to_f64().unwrap(), 5.4);
    assert_eq!(V::from(5.4).to_f32().unwrap(), 5.4f32);
    assert_eq!(
        String::try_from(&V::from("test string")).unwrap(),
        "test string"
    );
    assert_eq!(V::from(vec![1, 2, 3]).to_int_vector().unwrap(), [1, 2, 3]);
    assert_eq!(
        V::from(vec![1.2, 2.3]).to_double_vector().unwrap(),
        [1.2, 2.3]
    );
    assert_eq!(
        V::from(vec![String::from("a")]).to_string_vector().unwrap(),
        ["a"]
    );
    assert!(V::from("true").to_bool().unwrap());
    assert!(!V::from("false").to_bool().unwrap());
    for wrong in [
        V::EMPTY,
        V::from("TRUE"),
        V::from(" false"),
        V::from("bla"),
        V::from(12),
        V::from(34.45),
    ] {
        assert!(wrong.to_bool().is_err());
    }
    assert_eq!(V::EMPTY.to_char().unwrap(), None);
    assert_eq!(
        V::from("hello\0world").to_char().unwrap(),
        Some("hello\0world")
    );
    assert!(V::from(12).to_char().is_err());
    for wrong in [
        V::EMPTY,
        V::from("5.4"),
        V::IntegerList(vec![]),
        V::FloatList(vec![5.4]),
    ] {
        assert!(wrong.to_f64().is_err());
        assert!(wrong.to_i64().is_err());
    }
    assert!(V::from(55.4).to_i64().is_err());
    assert!(V::from(-55).to_u64().is_err());
    assert!(V::from("a,b").to_string_vector().is_err());
    assert!(V::from(vec![1.2]).to_int_vector().is_err());
    assert!(V::from(vec![1]).to_double_vector().is_err());
    assert!(String::try_from(&V::from(12)).is_err());
}

#[test]
fn native_checked_integer_storage_and_narrowing_boundaries() {
    assert_eq!(
        V::try_from(i64::MAX as u64).unwrap().to_i64().unwrap(),
        i64::MAX
    );
    assert!(V::try_from(i64::MAX as u64 + 1).is_err());
    assert!(V::try_from(u64::MAX).is_err());
    assert_eq!(V::from(i64::MIN).to_i64().unwrap(), i64::MIN);
    assert_eq!(V::from(i16::MIN).to_i16().unwrap(), i16::MIN);
    assert_eq!(V::from(i16::MAX).to_i16().unwrap(), i16::MAX);
    assert!(V::from(i32::from(i16::MIN) - 1).to_i16().is_err());
    assert!(V::from(i32::from(i16::MAX) + 1).to_i16().is_err());
    assert!(V::from(i64::from(i32::MAX) + 1).to_i32().is_err());
    assert!(V::from(i64::from(u32::MAX) + 1).to_u32().is_err());
    assert!(V::from(-1).to_usize().is_err());
    assert!(V::from(f64::MAX).to_f32().is_err());
    assert_eq!(V::from(f64::from(f32::MAX)).to_f32().unwrap(), f32::MAX);
    assert_eq!(V::from(f64::from_bits(1)).to_f32().unwrap(), 0.0);
    // Source integers convert directly to f32, not via an intermediate f64.
    let n = (1i64 << 62) + (1i64 << 38) + 1;
    assert_eq!(V::from(n).to_f32().unwrap().to_bits(), (n as f32).to_bits());
}

#[test]
fn all_source_to_string_literal_goldens_and_stream_distinction() {
    let cases = [
        (V::EMPTY, "", ""),
        (V::from("hello"), "hello", "hello"),
        (V::from(5), "5", "5"),
        (V::from(47.11), "47.109999999999999", "47.11"),
        (V::from(-23456.78), "-2.345678e04", "-2.346e04"),
        (
            V::from(vec![
                "test string".to_owned(),
                "string2".to_owned(),
                "last string".to_owned(),
            ]),
            "[test string, string2, last string]",
            "[test string, string2, last string]",
        ),
        (
            V::from(vec![1, 2, 3, 4, 5]),
            "[1, 2, 3, 4, 5]",
            "[1, 2, 3, 4, 5]",
        ),
        (
            V::from(vec![1.2, 47.11, 1.2345678e05]),
            "[1.2, 47.109999999999999, 1.2345678e05]",
            "[1.2, 47.11, 1.235e05]",
        ),
    ];
    for (value, full, low) in cases {
        assert_eq!(value.to_text(true).unwrap(), full);
        assert_eq!(value.to_text(false).unwrap(), low);
    }
    let stream: String = [
        V::from(5),
        V::from(100u32),
        V::from(1.111),
        V::from(1.1),
        V::from("hello "),
        V::from("world"),
        V::EMPTY,
    ]
    .iter()
    .map(|v| v.to_stream_text().unwrap())
    .collect();
    assert_eq!(stream, "51001.1111.1hello world");
    assert_eq!(V::from(5.0).to_text(true).unwrap(), "5.0");
    assert_eq!(V::from(5.0).to_stream_text().unwrap(), "5");
    assert_eq!(V::from(1e6).to_stream_text().unwrap(), "1e+06");
    assert_eq!(V::from(1e-5).to_stream_text().unwrap(), "1e-05");
    assert_eq!(
        V::from(vec!["a,b".to_owned(), "[x]".to_owned(), "".to_owned()])
            .to_text(true)
            .unwrap(),
        "[a,b, [x], ]"
    );
}

#[test]
fn float_boundaries_precision_and_special_values() {
    for (value, expected) in [
        (0.0, "0.0"),
        (-0.0, "-0.0"),
        (0.01, "0.01"),
        (10000.0, "1.0e04"),
        (1e-4, "1.0e-04"),
        (1e100, "1.0e100"),
        (f64::from_bits(1), "5.0e-324"),
        (f64::INFINITY, "inf"),
        (f64::NEG_INFINITY, "-inf"),
        (f64::NAN, "NaN"),
    ] {
        assert_eq!(V::from(value).to_text(true).unwrap(), expected, "{value:?}");
    }
    assert_eq!(
        V::from(f64::from_bits(0.01f64.to_bits() - 1))
            .to_text(true)
            .unwrap(),
        "9.999999999999998e-03"
    );
    assert_eq!(
        V::from(f64::from_bits(10000.0f64.to_bits() - 1))
            .to_text(true)
            .unwrap(),
        "9999.999999999998181"
    );
    assert_eq!(V::from(f64::INFINITY).to_f32().unwrap(), f32::INFINITY);
    assert!(V::from(f64::NAN).to_f32().unwrap().is_nan());
    let negative_zero = V::from(-0.0);
    assert_eq!(
        negative_zero.to_f64().unwrap().to_bits(),
        (-0.0f64).to_bits()
    );
    assert_ne!(V::from(f64::NAN), V::from(f64::NAN));
    assert_ne!(V::from(vec![f64::NAN]), V::from(vec![f64::NAN]));
    assert_eq!(
        V::from(vec![f64::INFINITY, f64::NEG_INFINITY, f64::NAN])
            .to_text(false)
            .unwrap(),
        "[inf, -inf, NaN]"
    );
}

#[test]
fn source_comparisons_are_type_strict_and_lists_use_size_only() {
    assert_ne!(V::from(5), V::from(5.0));
    assert_ne!(V::from("5"), V::from(5));
    assert_eq!(V::from(-0.0), V::from(0.0));
    assert!(V::from("a").source_less(&V::from("b")));
    assert!(V::from(5).source_greater(&V::from(4)));
    assert!(V::from(1.1).source_less(&V::from(2.2)));
    assert!(!V::EMPTY.source_less(&V::EMPTY));
    assert!(!V::from(f64::NAN).source_less(&V::from(1.0)));
    assert!(!V::from(f64::NAN).source_greater(&V::from(1.0)));
    assert_eq!(V::from(5).partial_cmp(&V::from(5.0)), None);
    for (a, b, c) in [
        (V::from(vec![9]), V::from(vec![-9]), V::from(vec![0, 0])),
        (
            V::from(vec![9.0]),
            V::from(vec![-9.0]),
            V::from(vec![0.0, 0.0]),
        ),
        (
            V::from(vec![String::from("z")]),
            V::from(vec![String::from("a")]),
            V::from(vec![String::new(), String::new()]),
        ),
    ] {
        assert_ne!(a, b);
        assert!(!a.source_less(&b));
        assert!(!a.source_greater(&b));
        assert_eq!(a.partial_cmp(&b), None);
        assert!(a.source_less(&c));
        assert!(c.source_greater(&a));
        assert_eq!(a.partial_cmp(&a), Some(Ordering::Equal));
    }
    assert!(V::from(vec![f64::NAN]).source_less(&V::from(vec![0.0, 0.0])));
}

#[test]
#[allow(clippy::approx_constant)] // Literal 3.14 is a pinned ParamValue test input.
fn source_hash_recurrence_literals_and_signed_zero_contract() {
    // Independently evaluated HashUtils.h FNV-1a/golden-ratio recurrence using
    // explicit little-endian bytes, not a hash generated by this implementation.
    for (value, expected) in [
        (V::EMPTY, 0xaf63bb4c8601b479),
        (V::from("test"), 0x33851e51bcb4bb2e),
        (V::from(42), 0x0d59e515dec3f9a3),
        (V::from(3.14), 0x0a4ed4905e3cb625),
        (V::from(0.0), 0xe4ab8baccf524eae),
        (
            V::from(vec!["a".to_owned(), "b".to_owned()]),
            0x4e8c47b9f79c4af8,
        ),
        (V::from(vec![1, -2, 3]), 0x48b4091067caadb8),
        (V::from(vec![1.1, 2.2]), 0x56df3d4e28346bec),
    ] {
        assert_eq!(value.source_hash64().unwrap(), expected);
    }
    fn rust_hash(value: &V) -> u64 {
        let mut h = std::collections::hash_map::DefaultHasher::new();
        value.hash(&mut h);
        h.finish()
    }
    assert_eq!(rust_hash(&V::from(0.0)), rust_hash(&V::from(-0.0)));
    assert_eq!(
        rust_hash(&V::from(vec![0.0])),
        rust_hash(&V::from(vec![-0.0]))
    );
    assert_eq!(V::from(-0.0).source_hash64().unwrap(), 0xe4ab8baccf524eae);
    assert_ne!(rust_hash(&V::from(5)), rust_hash(&V::from(5.0)));
}

#[test]
fn bounded_conversions_fail_without_changing_values() {
    let value = V::StringList(vec![String::new(); 1_000_001]);
    assert!(value.to_text(true).is_err());
    assert!(value.to_string_vector().is_err());
    assert!(value.checked_clone().is_err());
    assert!(value.source_hash64().is_err());
    assert_eq!(value.as_string_list().unwrap().len(), 1_000_001);
    let floats = V::FloatList(vec![1.0; 1_000_000]);
    assert!(floats.to_text(true).is_err()); // formatting work is precharged in full
    assert_eq!(floats.as_float_list().unwrap()[999_999], 1.0);
    let original = V::from(vec![String::from("owned")]);
    let mut copy = original.to_string_vector().unwrap();
    copy[0].push('!');
    assert_eq!(original.as_string_list().unwrap()[0], "owned");
}
