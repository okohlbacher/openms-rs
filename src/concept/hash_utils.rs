// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Hash utilities for OpenMS classes: `CONCEPT/HashUtils.h`.
//!
//! The source header provides hash functions that are consistent (equal inputs
//! give equal hashes), fast (the FNV-1a algorithm) and low in collisions
//! (combining with golden-ratio mixing). They satisfy the requirement a hash
//! container places on a hash: if `a == b` then `hash(a) == hash(b)`. Over
//! stable identifiers (strings, not pointers) the values are reproducible
//! across process runs.
//!
//! # API mapping
//!
//! | Source (`HashUtils.h` at core `bc9cc12`) | Rust |
//! |---|---|
//! | `fnv1a_hash_bytes(const void*, std::size_t)` | [`fnv1a_hash_bytes`] over a byte slice |
//! | `fnv1a_hash_string(const std::string&)` | [`fnv1a_hash_string`] over `&str` |
//! | `hash_combine(std::size_t&, std::size_t)` | [`hash_combine`] |
//! | `template<T> hash_int(T)` (integral `T`) | [`hash_int`] over [`HashInt`] |
//! | `hash_char(char)` | [`hash_char`] |
//! | `template<T> hash_float(T)` (floating-point `T`) | [`hash_float`] over [`HashFloat`] |
//!
//! # Differences from the source
//!
//! - **Width.** The source returns `std::size_t` and starts from 64-bit
//!   constants, so on a 32-bit host every step is truncated to 32 bits and the
//!   values differ. These functions are 64-bit on every target and return
//!   `u64`; on the 64-bit hosts the source is built for, the values are equal.
//! - **Byte order.** The source's `@note` warns that `hash_int` and
//!   `hash_float` hash "the platform's native byte representation" and differ
//!   between big- and little-endian systems, which is why it calls them
//!   unsuitable for persistent storage. [`hash_int`] and [`hash_float`] hash
//!   the little-endian bytes on every target: the same value as the source on
//!   the little-endian x86_64 and arm64 hosts, and the same value everywhere,
//!   so the endianness caveat does not apply. The 32-bit width caveat above
//!   is likewise gone.
//! - **`long double`.** `hash_float` also accepts `long double`, whose width and
//!   padding depend on the platform. Rust has no such type, so [`HashFloat`]
//!   is implemented for `f32` and `f64` only.
//! - **Integral types.** C++ `std::is_integral` admits `bool`, the character
//!   types and every integer width. [`HashInt`] covers `bool`, `char` (4 bytes,
//!   as C++ `char32_t`), and `i8` to `i128`, `u8` to `u128`, `isize`, `usize`;
//!   a byte-sized C++ `char` is [`hash_char`] or `hash_int(u8)`/`hash_int(i8)`.
//!
//! These are the functions the ported `ParamValue::source_hash64` and the
//! controlled-vocabulary FNV hash compute; other native types hash through
//! `std::hash::Hash` with the caller's hasher and do not use them.

/// FNV-1a 64-bit offset basis, `14695981039346656037`.
pub const FNV_OFFSET_BASIS: u64 = 14_695_981_039_346_656_037;

/// FNV-1a 64-bit prime, `1099511628211`.
pub const FNV_PRIME: u64 = 1_099_511_628_211;

/// The 64-bit golden-ratio constant `0x9e3779b97f4a7c15` (`2^64 / phi`) that
/// [`hash_combine`] adds.
pub const GOLDEN_RATIO: u64 = 0x9e37_79b9_7f4a_7c15;

/// FNV-1a hash of a byte sequence: source `fnv1a_hash_bytes(data, size)`.
///
/// FNV-1a is a non-cryptographic hash: for each byte the hash is XORed with
/// the byte and multiplied by [`FNV_PRIME`] modulo 2^64, starting from
/// [`FNV_OFFSET_BASIS`]. It is fast, simple to verify and distributes well for
/// hash tables. The empty sequence hashes to the offset basis.
///
/// # Examples
///
/// ```
/// use openms::concept::hash_utils::{fnv1a_hash_bytes, FNV_OFFSET_BASIS};
///
/// assert_eq!(fnv1a_hash_bytes(b""), FNV_OFFSET_BASIS);
/// assert_eq!(fnv1a_hash_bytes(b"a"), 0xaf63_dc4c_8601_ec8c);
/// ```
pub const fn fnv1a_hash_bytes(data: &[u8]) -> u64 {
    let mut hash = FNV_OFFSET_BASIS;
    let mut index = 0;
    while index < data.len() {
        hash ^= data[index] as u64;
        hash = hash.wrapping_mul(FNV_PRIME);
        index += 1;
    }
    hash
}

/// FNV-1a hash of a string's UTF-8 bytes: source `fnv1a_hash_string(s)`.
///
/// The source offers it as "a consistent alternative to `std::hash<std::string>`,
/// which may vary between standard library implementations"; the same holds
/// against Rust's randomly seeded `HashMap` hasher.
///
/// # Examples
///
/// ```
/// use openms::concept::hash_utils::fnv1a_hash_string;
///
/// assert_eq!(fnv1a_hash_string("foobar"), 0x8594_4171_f739_67e8);
/// ```
pub const fn fnv1a_hash_string(text: &str) -> u64 {
    fnv1a_hash_bytes(text.as_bytes())
}

/// Combines `value` into the hash `seed` with golden-ratio mixing: source
/// `hash_combine(seed, value)`.
///
/// This is Boost's formula, `seed ^= value + 0x9e3779b97f4a7c15 + (seed << 6) +
/// (seed >> 2)`, with every addition and the left shift wrapping modulo 2^64 as
/// the source's unsigned arithmetic does. The constant is [`GOLDEN_RATIO`]. The
/// source updates `seed` in place; so does this function.
///
/// # Examples
///
/// ```
/// use openms::concept::hash_utils::{hash_combine, GOLDEN_RATIO};
///
/// let mut seed = 0;
/// hash_combine(&mut seed, 0);
/// assert_eq!(seed, GOLDEN_RATIO);
/// ```
pub fn hash_combine(seed: &mut u64, value: u64) {
    let mixed = value
        .wrapping_add(GOLDEN_RATIO)
        .wrapping_add(*seed << 6)
        .wrapping_add(*seed >> 2);
    *seed ^= mixed;
}

/// An integral type [`hash_int`] accepts: source `static_assert(std::is_integral_v<T>)`.
pub trait HashInt: Copy {
    /// The value's bytes, least significant first.
    fn hash_bytes(self) -> Vec<u8>;
}

macro_rules! hash_int_impl {
    ($($integer:ty),*) => {$(
        impl HashInt for $integer {
            fn hash_bytes(self) -> Vec<u8> {
                self.to_le_bytes().to_vec()
            }
        }
    )*};
}
hash_int_impl!(
    i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize
);

impl HashInt for bool {
    /// One byte, 0 or 1: C++ `sizeof(bool)` is 1 on every supported target.
    fn hash_bytes(self) -> Vec<u8> {
        vec![u8::from(self)]
    }
}

impl HashInt for char {
    /// The four little-endian bytes of the code point, as C++ `char32_t`.
    fn hash_bytes(self) -> Vec<u8> {
        u32::from(self).to_le_bytes().to_vec()
    }
}

/// Hash of an integer: FNV-1a over its bytes, source `hash_int(value)`.
///
/// The source hashes the native byte representation and `sizeof(T)` bytes;
/// this hashes the little-endian representation of the same width, which is
/// the native one on the hosts the source is built for (see the module notes).
///
/// # Examples
///
/// ```
/// use openms::concept::hash_utils::{fnv1a_hash_bytes, hash_int};
///
/// assert_eq!(hash_int(1u32), fnv1a_hash_bytes(&[1, 0, 0, 0]));
/// assert_eq!(hash_int(-1i16), fnv1a_hash_bytes(&[0xff, 0xff]));
/// ```
pub fn hash_int<T: HashInt>(value: T) -> u64 {
    fnv1a_hash_bytes(&value.hash_bytes())
}

/// Hash of one byte-sized character: source `hash_char(c)`.
///
/// The source casts its `char` to `unsigned char` and applies one FNV-1a step,
/// the same value as [`fnv1a_hash_bytes`] of that single byte. A signed C++
/// `char` of `-1` is the byte `0xff`; pass `c as u8`.
///
/// # Examples
///
/// ```
/// use openms::concept::hash_utils::{fnv1a_hash_string, hash_char};
///
/// assert_eq!(hash_char(b'a'), fnv1a_hash_string("a"));
/// ```
pub const fn hash_char(character: u8) -> u64 {
    (FNV_OFFSET_BASIS ^ character as u64).wrapping_mul(FNV_PRIME)
}

/// A floating-point type [`hash_float`] accepts: source
/// `static_assert(std::is_floating_point_v<T>)`, for `float` and `double`.
pub trait HashFloat: Copy {
    /// The value's bytes, least significant first, after `-0.0` became `+0.0`.
    fn normalized_hash_bytes(self) -> Vec<u8>;
}

impl HashFloat for f32 {
    fn normalized_hash_bytes(self) -> Vec<u8> {
        let value = if self == 0.0 { 0.0f32 } else { self };
        value.to_le_bytes().to_vec()
    }
}

impl HashFloat for f64 {
    fn normalized_hash_bytes(self) -> Vec<u8> {
        let value = if self == 0.0 { 0.0f64 } else { self };
        value.to_le_bytes().to_vec()
    }
}

/// Hash of a floating-point value: FNV-1a over its raw bytes, source
/// `hash_float(value)`.
///
/// Hashing the raw bytes is consistent with `==` on floating-point values, with
/// one repair the source makes and this keeps: `-0.0` and `+0.0` compare equal
/// but differ in their sign bit, so `-0.0` is normalised to `+0.0` first and
/// both hash alike. NaN compares unequal to itself, so the hash contract places
/// no demand on it; distinct NaN payloads hash differently, as in the source.
///
/// # Examples
///
/// ```
/// use openms::concept::hash_utils::hash_float;
///
/// assert_eq!(hash_float(-0.0f64), hash_float(0.0f64));
/// assert_eq!(hash_float(-0.0f32), hash_float(0.0f32));
/// ```
pub fn hash_float<T: HashFloat>(value: T) -> u64 {
    fnv1a_hash_bytes(&value.normalized_hash_bytes())
}
