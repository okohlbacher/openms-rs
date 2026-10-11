// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `CONCEPT/CommonEnums.h`, `CONCEPT/Types.h` and `CONCEPT/HashUtils.h`.
//!
//! None of the three headers has a class test at core `bc9cc12`. Expected
//! values are the header's own literals (tier 3: the name tables, enum
//! discriminants, ASCII codes and `numeric_limits<T>::digits10` values) or are
//! derived independently of the port (tier 4: the published FNV-1a 64-bit test
//! vectors, and FNV-1a and Boost `hash_combine` recomputed in Python from the
//! header's formulas over explicit little-endian bytes). See
//! `docs/CONCEPT_BASICS_SUPPORT.md`.

use openms::chemistry::hydrophobicity::HydrophobicityScale;
use openms::concept::common_enums::{DIM_NAMES, DIM_NAMES_SHORT, DimUnit, MZ_UNIT_NAMES, MzUnits};
use openms::concept::hash_utils::{
    FNV_OFFSET_BASIS, FNV_PRIME, GOLDEN_RATIO, fnv1a_hash_bytes, fnv1a_hash_string, hash_char,
    hash_combine, hash_float, hash_int,
};
use openms::concept::types::{self, LONG_DOUBLE_WRITTEN_DIGITS, ascii, written_digits};

// ---------------------------------------------------------------- CommonEnums

#[test]
fn dim_unit_discriminants_and_names_follow_the_header() {
    // CommonEnums.h: RT = 0 counting up, SIZE_OF_DIM_UNITS = 7.
    assert_eq!(DimUnit::SIZE_OF_DIM_UNITS, 7);
    let expected = [
        (DimUnit::Rt, 0, "RT [s]", "RT"),
        (DimUnit::Mz, 1, "m/z [Th]", "m/z"),
        (DimUnit::Int, 2, "intensity", "int"),
        (DimUnit::ImMs, 3, "IM [milliseconds]", "IM"),
        (DimUnit::ImVssc, 4, "IM [vs / cm2]", "IM"),
        (DimUnit::FaimsCv, 5, "FAIMS CV", "FCV"),
        (DimUnit::ImCcs, 6, "CCS [angstrom^2]", "CCS"),
    ];
    for (unit, index, name, short) in expected {
        assert_eq!(unit.index(), index);
        assert_eq!(unit as usize, index);
        assert_eq!(DimUnit::from_index(index), Some(unit));
        assert_eq!(DimUnit::ALL[index], unit);
        assert_eq!(unit.name(), name);
        assert_eq!(DIM_NAMES[index], name);
        assert_eq!(unit.short_name(), short);
        assert_eq!(DIM_NAMES_SHORT[index], short);
        assert_eq!(unit.to_string(), name);
    }
    assert_eq!(DimUnit::from_index(DimUnit::SIZE_OF_DIM_UNITS), None);
    assert_eq!(DimUnit::from_index(usize::MAX), None);
}

#[test]
fn mz_units_discriminants_names_and_parsing() {
    assert_eq!(MzUnits::SIZE_OF_MZ_UNITS, 2);
    assert_eq!(MZ_UNIT_NAMES, ["Da", "ppm"]);
    assert_eq!(MzUnits::Da.index(), 0);
    assert_eq!(MzUnits::Ppm.index(), 1);
    assert_eq!(MzUnits::from_index(0), Some(MzUnits::Da));
    assert_eq!(MzUnits::from_index(1), Some(MzUnits::Ppm));
    assert_eq!(MzUnits::from_index(2), None);
    assert_eq!(MzUnits::Da.name(), "Da");
    assert_eq!(MzUnits::Ppm.to_string(), "ppm");
    assert_eq!("Da".parse::<MzUnits>().ok(), Some(MzUnits::Da));
    assert_eq!("ppm".parse::<MzUnits>().ok(), Some(MzUnits::Ppm));
    // The lookup in MZ_UNIT_NAMES is exact, so case and spacing matter.
    for bad in ["da", "PPM", " ppm", "", "Th"] {
        assert!(bad.parse::<MzUnits>().is_err(), "{bad:?}");
    }
}

#[test]
fn hydrophobicity_scale_keeps_the_source_discriminants() {
    // HydrophobicityScaleMethod: KYTE_DOOLITTLE = 0 .. EISENBERG_CONSENSUS = 6.
    let expected = [
        HydrophobicityScale::KyteDoolittle,
        HydrophobicityScale::Eisenberg,
        HydrophobicityScale::HoppWoods,
        HydrophobicityScale::BullBreese,
        HydrophobicityScale::BlackMould,
        HydrophobicityScale::Guy,
        HydrophobicityScale::EisenbergConsensus,
    ];
    for (index, scale) in expected.into_iter().enumerate() {
        assert_eq!(scale as usize, index);
        assert_eq!(HydrophobicityScale::ALL[index], scale);
    }
}

// ---------------------------------------------------------------- Types

#[test]
fn type_aliases_have_the_source_widths() {
    assert_eq!(std::mem::size_of::<types::Int32>(), 4);
    assert_eq!(std::mem::size_of::<types::UInt32>(), 4);
    assert_eq!(std::mem::size_of::<types::Int64>(), 8);
    assert_eq!(std::mem::size_of::<types::UInt64>(), 8);
    assert_eq!(std::mem::size_of::<types::Int>(), 4);
    assert_eq!(std::mem::size_of::<types::UInt>(), 4);
    assert_eq!(std::mem::size_of::<types::Byte>(), 1);
    assert_eq!(std::mem::size_of::<types::Uid>(), 8);
    assert_eq!(std::mem::size_of::<types::Time>(), 8);
    assert_eq!(
        std::mem::size_of::<types::Size>(),
        std::mem::size_of::<usize>()
    );
    assert_eq!(
        std::mem::size_of::<types::SignedSize>(),
        std::mem::size_of::<isize>()
    );
    // Signedness: Int32/Int64/Int/Time/SignedSize are signed; the rest are not.
    // Checked at compile time, so they cannot silently stop running.
    const _: () = assert!(types::Int32::MIN < 0 && types::Int64::MIN < 0 && types::Int::MIN < 0);
    const _: () = assert!(types::Time::MIN < 0 && types::SignedSize::MIN < 0);
    assert_eq!(types::UInt32::MIN, 0);
    assert_eq!(types::Byte::MIN, 0);
}

#[test]
fn ascii_codes_follow_the_header() {
    let expected: [(u8, u8); 15] = [
        (ascii::BACKSPACE, 8),
        (ascii::BELL, 7),
        (ascii::CARRIAGE_RETURN, 13),
        (ascii::HORIZONTAL_TAB, 9),
        (ascii::NEWLINE, 10),
        (ascii::RETURN, 10),
        (ascii::SPACE, 32),
        (ascii::TAB, 9),
        (ascii::VERTICAL_TAB, 11),
        (ascii::COLON, b':'),
        (ascii::COMMA, b','),
        (ascii::EXCLAMATION_MARK, b'!'),
        (ascii::POINT, b'.'),
        (ascii::QUESTION_MARK, b'?'),
        (ascii::SEMICOLON, b';'),
    ];
    for (actual, code) in expected {
        assert_eq!(actual, code);
    }
    // ASCII__RETURN is defined as ASCII__NEWLINE, not as the carriage return.
    assert_eq!(ascii::RETURN, ascii::NEWLINE);
    assert_ne!(ascii::RETURN, ascii::CARRIAGE_RETURN);
    assert_eq!(ascii::TAB, ascii::HORIZONTAL_TAB);
}

#[test]
fn written_digits_follow_each_specialisation() {
    // float and double: numeric_limits<T>::digits10.
    assert_eq!(written_digits::<f32>(), 6);
    assert_eq!(written_digits::<f64>(), 15);
    // int and unsigned int: numeric_limits<int / unsigned int>::digits10 = 9.
    assert_eq!(written_digits::<i32>(), 9);
    assert_eq!(written_digits::<u32>(), 9);
    // long and unsigned long deliberately return the int precision, not 18/19.
    assert_eq!(written_digits::<i64>(), 9);
    assert_eq!(written_digits::<u64>(), 9);
    assert_eq!(written_digits::<isize>(), 9);
    assert_eq!(written_digits::<usize>(), 9);
    // Every other type takes the generic template's default stream precision 6.
    assert_eq!(written_digits::<i8>(), 6);
    assert_eq!(written_digits::<u8>(), 6);
    assert_eq!(written_digits::<i16>(), 6);
    assert_eq!(written_digits::<u16>(), 6);
    assert_eq!(written_digits::<i128>(), 6);
    assert_eq!(written_digits::<u128>(), 6);
    // x86_64 Linux long double (80-bit extended): digits10 = 18.
    assert_eq!(LONG_DOUBLE_WRITTEN_DIGITS, 18);
}

#[test]
fn written_digits_are_digits10_of_the_mantissa() {
    // numeric_limits<T>::digits10 is floor((p - 1) * log10(2)) for a p-bit
    // significand: the largest d with 10^d <= 2^(p - 1). The header's rationale
    // (2^24 / 10^6 = 16.777216, 2^53 / 10^15 = 9.007199254740992) rests on it.
    let f32_digits = u32::try_from(written_digits::<f32>()).unwrap_or(0);
    let f64_digits = u32::try_from(written_digits::<f64>()).unwrap_or(0);
    let f32_half = 1u64 << (f32::MANTISSA_DIGITS - 1);
    let f64_half = 1u64 << (f64::MANTISSA_DIGITS - 1);
    assert!(10u64.pow(f32_digits) <= f32_half && f32_half < 10u64.pow(f32_digits + 1));
    assert!(10u64.pow(f64_digits) <= f64_half && f64_half < 10u64.pow(f64_digits + 1));
    assert_eq!(written_digits::<f32>(), f32::DIGITS as i32);
    assert_eq!(written_digits::<f64>(), f64::DIGITS as i32);
}

// ---------------------------------------------------------------- HashUtils

#[test]
fn fnv_constants_are_the_published_64_bit_parameters() {
    assert_eq!(FNV_OFFSET_BASIS, 14_695_981_039_346_656_037);
    assert_eq!(FNV_OFFSET_BASIS, 0xcbf2_9ce4_8422_2325);
    assert_eq!(FNV_PRIME, 1_099_511_628_211);
    assert_eq!(FNV_PRIME, 0x0000_0100_0000_01b3);
    assert_eq!(GOLDEN_RATIO, 0x9e37_79b9_7f4a_7c15);
}

#[test]
fn fnv1a_matches_the_published_test_vectors() {
    // FNV-1a 64-bit reference vectors (Fowler/Noll/Vo test suite).
    assert_eq!(fnv1a_hash_string(""), 0xcbf2_9ce4_8422_2325);
    assert_eq!(fnv1a_hash_string("a"), 0xaf63_dc4c_8601_ec8c);
    assert_eq!(fnv1a_hash_string("foobar"), 0x8594_4171_f739_67e8);
    assert_eq!(fnv1a_hash_bytes(b"foobar"), fnv1a_hash_string("foobar"));
    assert_eq!(fnv1a_hash_bytes(&[]), FNV_OFFSET_BASIS);
}

#[test]
fn hash_char_is_one_fnv_step_over_the_unsigned_byte() {
    assert_eq!(hash_char(b'a'), 0xaf63_dc4c_8601_ec8c);
    for byte in 0..=u8::MAX {
        assert_eq!(hash_char(byte), fnv1a_hash_bytes(&[byte]));
    }
    // A signed C++ char of -1 is the unsigned byte 0xff.
    assert_eq!(hash_char(0xff), 0xaf64_724c_8602_eb6e);
}

#[test]
fn hash_int_hashes_the_little_endian_bytes_of_its_width() {
    // Python: FNV-1a over struct.pack('<i', 1).
    assert_eq!(hash_int(1i32), 0xad2a_ca77_4798_5764);
    assert_eq!(hash_int(1u32), hash_int(1i32));
    assert_eq!(hash_int(1i32), fnv1a_hash_bytes(&[1, 0, 0, 0]));
    assert_eq!(hash_int(1i64), fnv1a_hash_bytes(&[1, 0, 0, 0, 0, 0, 0, 0]));
    assert_eq!(hash_int(-1i16), fnv1a_hash_bytes(&[0xff, 0xff]));
    assert_eq!(hash_int(0x1234u16), fnv1a_hash_bytes(&[0x34, 0x12]));
    assert_eq!(hash_int(255u8), 0xaf64_724c_8602_eb6e);
    assert_eq!(hash_int(true), fnv1a_hash_bytes(&[1]));
    assert_eq!(hash_int(false), fnv1a_hash_bytes(&[0]));
    assert_eq!(hash_int('a'), fnv1a_hash_bytes(&[0x61, 0, 0, 0]));
    assert_eq!(
        hash_int(usize::MAX),
        fnv1a_hash_bytes(&usize::MAX.to_le_bytes())
    );
    // The width is part of the value: the same number hashes differently.
    assert_ne!(hash_int(1u8), hash_int(1u32));
}

#[test]
fn hash_float_normalises_negative_zero_only() {
    // Python: FNV-1a over struct.pack('<d', 0.0) and struct.pack('<f', 0.0).
    assert_eq!(hash_float(0.0f64), 0xa8c7_f832_281a_39c5);
    assert_eq!(hash_float(-0.0f64), 0xa8c7_f832_281a_39c5);
    assert_eq!(hash_float(0.0f32), 0x4d25_767f_9dce_13f5);
    assert_eq!(hash_float(-0.0f32), 0x4d25_767f_9dce_13f5);
    // Python: FNV-1a over struct.pack('<d', 1.5).
    assert_eq!(hash_float(1.5f64), 0xaa95_e932_29a2_7c80);
    assert_eq!(
        hash_float(-1.5f64),
        fnv1a_hash_bytes(&(-1.5f64).to_le_bytes())
    );
    assert_ne!(hash_float(1.5f64), hash_float(-1.5f64));
    assert_eq!(hash_float(2.5f32), fnv1a_hash_bytes(&2.5f32.to_le_bytes()));
    // NaN is hashed by its bits; distinct payloads differ, as in the source.
    let quiet = f64::from_bits(0x7ff8_0000_0000_0000);
    let payload = f64::from_bits(0x7ff8_0000_0000_0001);
    assert_eq!(
        hash_float(quiet),
        fnv1a_hash_bytes(&quiet.to_bits().to_le_bytes())
    );
    assert_ne!(hash_float(quiet), hash_float(payload));
}

#[test]
fn hash_combine_is_the_boost_golden_ratio_recurrence() {
    let mut seed = 0;
    hash_combine(&mut seed, 0);
    assert_eq!(seed, GOLDEN_RATIO);

    // 1 ^ (2 + GR + (1 << 6) + (1 >> 2))
    let mut seed = 1;
    hash_combine(&mut seed, 2);
    assert_eq!(seed, 0x9e37_79b9_7f4a_7c56);

    // Python, modulo 2^64.
    let mut seed = 0xdead_beef;
    hash_combine(&mut seed, 0x1234_5678);
    assert_eq!(seed, 0x9e37_79f1_aa34_40e7);

    // Every addition and the left shift wrap.
    let mut seed = u64::MAX;
    hash_combine(&mut seed, u64::MAX);
    assert_eq!(seed, 0x21c8_8646_80b5_842c);
}

#[test]
fn equal_values_hash_equally() {
    // The contract the header states: if a == b then hash(a) == hash(b).
    let values = [
        0.0f64,
        -0.0,
        1.0,
        -1.0,
        f64::MIN_POSITIVE,
        f64::MAX,
        f64::INFINITY,
    ];
    for a in values {
        for b in values {
            if a == b {
                assert_eq!(hash_float(a), hash_float(b), "{a} {b}");
            }
        }
    }
}
