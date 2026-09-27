// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Basic type aliases, ASCII codes and output precisions: `CONCEPT/Types.h`.
//!
//! The crate writes Rust's own primitive types throughout (`i32`, `usize`,
//! `u8`, ...); the aliases here exist so that a source signature can be read
//! off one-to-one and so that the whole header has a counterpart. The source
//! comments' descriptions are carried at each item.
//!
//! # API mapping
//!
//! | Source (`Types.h` at core `bc9cc12`) | Rust |
//! |---|---|
//! | `Int32`, `UInt32`, `Int64`, `UInt64` | [`Int32`], [`UInt32`], [`Int64`], [`UInt64`] |
//! | `Time` (`time_t`) | [`Time`], `i64` |
//! | `UInt` (`unsigned int`), `Int` (`int`) | [`UInt`], [`Int`] |
//! | `Byte` (`uint8_t`), `UID` (`uint64_t`) | [`Byte`], [`Uid`] |
//! | `Size` (`size_t`), `SignedSize` (`ptrdiff_t`) | [`Size`], [`SignedSize`] |
//! | `enum ASCII` and its 15 `ASCII__*` enumerators | the [`ascii`] constants |
//! | `writtenDigits<T>()`, its generic fallback and the `float`, `double`, `int`, `unsigned int`, `long int`, `unsigned long int` and `long double` specialisations | [`written_digits`] over [`WrittenDigits`] |
//!
//! # `time_t`
//!
//! [`Time`] is `i64`. `time_t` is a signed 64-bit count of seconds on every
//! 64-bit target and on 32-bit Windows; only 32-bit glibc and musl builds
//! without 64-bit time keep a 32-bit `time_t`. The port has no C interface
//! through which a `time_t` passes, so there is no ABI to match, and the wider
//! type removes the year-2038 limit rather than adding one.
//!
//! # `long double`
//!
//! Rust has no `long double`. The source's specialisation returns
//! `std::numeric_limits<long double>::digits10` except on Windows, where it
//! returns `double`'s 15 because the Microsoft I/O libraries treat `long
//! double` as `double`. [`LONG_DOUBLE_WRITTEN_DIGITS`] records the value on
//! the x86_64 Linux reference build (18, the 80-bit extended format); no Rust
//! type carries it.

/// Signed integer type (32 bit): source `Int32`.
pub type Int32 = i32;
/// Unsigned integer type (32 bit): source `UInt32`.
pub type UInt32 = u32;
/// Signed integer type (64 bit): source `Int64`.
pub type Int64 = i64;
/// Unsigned integer type (64 bit): source `UInt64`.
pub type UInt64 = u64;
/// A point in time as seconds since the epoch: source `Time`, "a synonym for
/// `time_t`". See the module notes on its width.
pub type Time = i64;
/// Unsigned integer type: source `UInt`, an `unsigned int` (32 bit on every
/// supported target).
pub type UInt = u32;
/// Signed integer type: source `Int`, an `int` (32 bit on every supported
/// target).
pub type Int = i32;
/// Byte data, 8 bits long and always unsigned: source `Byte`.
pub type Byte = u8;
/// A unique object ID as an unsigned 64-bit value: source `UID`. The owned
/// identifier type of the port is `concept::unique_id::UniqueId`.
pub type Uid = u64;
/// Size type, which can hold the result of `size()`: source `Size`.
pub type Size = usize;
/// Signed size type, e.g. for a pointer difference: source `SignedSize`
/// (`ptrdiff_t`).
pub type SignedSize = isize;

/// The source's `enum ASCII` character codes, as bytes.
///
/// The enumerators are named `ASCII__<NAME>` in the source; here they are
/// `ascii::<NAME>`. `RETURN` equals `NEWLINE` and `TAB` equals
/// `HORIZONTAL_TAB`, as in the source.
pub mod ascii {
    /// `ASCII__BACKSPACE`, `'\b'`.
    pub const BACKSPACE: u8 = 0x08;
    /// `ASCII__BELL`, `'\a'`.
    pub const BELL: u8 = 0x07;
    /// `ASCII__CARRIAGE_RETURN`, `'\r'`.
    pub const CARRIAGE_RETURN: u8 = b'\r';
    /// `ASCII__HORIZONTAL_TAB`, `'\t'`.
    pub const HORIZONTAL_TAB: u8 = b'\t';
    /// `ASCII__NEWLINE`, `'\n'`.
    pub const NEWLINE: u8 = b'\n';
    /// `ASCII__RETURN`, defined as `ASCII__NEWLINE` (`'\n'`, not `'\r'`).
    pub const RETURN: u8 = NEWLINE;
    /// `ASCII__SPACE`, `' '`.
    pub const SPACE: u8 = b' ';
    /// `ASCII__TAB`, defined as `ASCII__HORIZONTAL_TAB`.
    pub const TAB: u8 = HORIZONTAL_TAB;
    /// `ASCII__VERTICAL_TAB`, `'\v'`.
    pub const VERTICAL_TAB: u8 = 0x0b;
    /// `ASCII__COLON`, `':'`.
    pub const COLON: u8 = b':';
    /// `ASCII__COMMA`, `','`.
    pub const COMMA: u8 = b',';
    /// `ASCII__EXCLAMATION_MARK`, `'!'`.
    pub const EXCLAMATION_MARK: u8 = b'!';
    /// `ASCII__POINT`, `'.'`.
    pub const POINT: u8 = b'.';
    /// `ASCII__QUESTION_MARK`, `'?'`.
    pub const QUESTION_MARK: u8 = b'?';
    /// `ASCII__SEMICOLON`, `';'`.
    pub const SEMICOLON: u8 = b';';
}

/// `std::numeric_limits<long double>::digits10` on the x86_64 Linux reference
/// build: what the source's `writtenDigits<long double>()` returns there. On
/// Windows the source returns 15 instead. See the module notes.
pub const LONG_DOUBLE_WRITTEN_DIGITS: Int = 18;

/// The number of digits commonly used for writing a value of a type (a.k.a.
/// precision): source `writtenDigits<T>()`.
///
/// The source unifies precision handling with these functions ("so please
/// don't use ad-hoc numbers"). A type decimal digit count is what the type
/// represents without loss: 6 for single and 15 for double precision, since
/// `2^24 / 10^6 = 16.777216` and `2^53 / 10^15 = 9.007199254740992`, so
/// rounding removes the remaining difference.
///
/// The default, [`WrittenDigits::WRITTEN_DIGITS`] = 6, is the source's generic
/// template, "the default precision of 6 according to 27.4.4.1
/// basic_ios constructors (C++ Standard)"; any type can opt in with an empty
/// `impl`. The specialisations are:
///
/// | Type | Source specialisation | Value |
/// |---|---|---|
/// | `f32` | `float`: `numeric_limits<float>::digits10` | 6 |
/// | `f64` | `double`: `numeric_limits<double>::digits10` | 15 |
/// | `i32` | `int`: `numeric_limits<int>::digits10` | 9 |
/// | `u32` | `unsigned int`: `numeric_limits<unsigned int>::digits10` | 9 |
/// | `i64`, `isize` | `long int`: deliberately `numeric_limits<int>::digits10` | 9 |
/// | `u64`, `usize` | `unsigned long int`: deliberately `numeric_limits<unsigned int>::digits10` | 9 |
/// | `i8`, `u8`, `i16`, `u16`, `i128`, `u128` | none, generic | 6 |
///
/// The source specialises the integer types only "not to bother people who
/// unintentionally provide an int argument"; the `long` specialisations return
/// the `int` precision, not `long`'s 18, and that quirk is kept. `i64`, `u64`,
/// `isize` and `usize` map to `long` and `unsigned long` as on the LP64 Linux
/// reference build, where `int64_t` and `size_t` are those types. (On Windows,
/// where `long` is 32 bits, `int64_t` is `long long` and falls to the generic
/// 6.)
pub trait WrittenDigits {
    /// The precision for this type.
    const WRITTEN_DIGITS: Int = 6;
}

impl WrittenDigits for f32 {
    const WRITTEN_DIGITS: Int = f32::DIGITS as Int;
}

impl WrittenDigits for f64 {
    const WRITTEN_DIGITS: Int = f64::DIGITS as Int;
}

/// `std::numeric_limits<int>::digits10` and `<unsigned int>::digits10`.
const INT_DIGITS10: Int = 9;

impl WrittenDigits for i32 {
    const WRITTEN_DIGITS: Int = INT_DIGITS10;
}

impl WrittenDigits for u32 {
    const WRITTEN_DIGITS: Int = INT_DIGITS10;
}

impl WrittenDigits for i64 {
    const WRITTEN_DIGITS: Int = INT_DIGITS10;
}

impl WrittenDigits for u64 {
    const WRITTEN_DIGITS: Int = INT_DIGITS10;
}

impl WrittenDigits for isize {
    const WRITTEN_DIGITS: Int = INT_DIGITS10;
}

impl WrittenDigits for usize {
    const WRITTEN_DIGITS: Int = INT_DIGITS10;
}

impl WrittenDigits for i8 {}
impl WrittenDigits for u8 {}
impl WrittenDigits for i16 {}
impl WrittenDigits for u16 {}
impl WrittenDigits for i128 {}
impl WrittenDigits for u128 {}

/// The number of digits commonly used for writing a `T`: source
/// `writtenDigits<T>(const T& = T())`.
///
/// The source's argument is unused and only selects the type; here the type
/// is the generic parameter. See [`WrittenDigits`] for the values.
///
/// # Examples
///
/// ```
/// use openms::concept::types::written_digits;
///
/// assert_eq!(written_digits::<f32>(), 6);
/// assert_eq!(written_digits::<f64>(), 15);
/// assert_eq!(written_digits::<i64>(), 9);
/// assert_eq!(written_digits::<u16>(), 6);
/// ```
pub const fn written_digits<T: WrittenDigits>() -> Int {
    T::WRITTEN_DIGITS
}
