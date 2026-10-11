// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Owned parameter values, distinct from finite-only identification metadata.

use super::ParamWork;
use crate::{Error, Result};
use std::cmp::Ordering;
use std::hash::{Hash, Hasher};

const MAX_VALUE_BYTES: usize = 64 * 1024 * 1024;
const MAX_VALUE_ELEMENTS: usize = 1_000_000;

/// Source value discriminants, including empty as a distinct seventh type.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum ParamValueType {
    String = 0,
    Integer = 1,
    Float = 2,
    StringList = 3,
    IntegerList = 4,
    FloatList = 5,
    Empty = 6,
}

/// Source parameter storage. NaN and infinities are valid floating-point values.
/// Equality is type-strict; NaN remains unequal to itself, so there is no `Eq`.
#[derive(Clone, Debug, Default, PartialEq)]
pub enum ParamValue {
    #[default]
    Empty,
    String(String),
    Integer(i64),
    Float(f64),
    StringList(Vec<String>),
    IntegerList(Vec<i32>),
    FloatList(Vec<f64>),
}

fn conversion(target: &str) -> Error {
    Error::InvalidValue(format!("parameter value cannot be converted to {target}"))
}
fn limit() -> Error {
    Error::InvalidValue("parameter value exceeds its resource limit".into())
}

macro_rules! integer_conversion {
    ($($name:ident: $type:ty),* $(,)?) => {$ (
        /// Convert an integer parameter, checking the target integer range.
        pub fn $name(&self) -> Result<$type> {
            <$type>::try_from(self.to_i64()?).map_err(|_| conversion(stringify!($type)))
        }
    )*};
}

impl ParamValue {
    pub const EMPTY: Self = Self::Empty;

    pub fn value_type(&self) -> ParamValueType {
        match self {
            Self::Empty => ParamValueType::Empty,
            Self::String(_) => ParamValueType::String,
            Self::Integer(_) => ParamValueType::Integer,
            Self::Float(_) => ParamValueType::Float,
            Self::StringList(_) => ParamValueType::StringList,
            Self::IntegerList(_) => ParamValueType::IntegerList,
            Self::FloatList(_) => ParamValueType::FloatList,
        }
    }

    pub fn is_empty(&self) -> bool {
        matches!(self, Self::Empty)
    }

    pub fn as_str(&self) -> Result<&str> {
        match self {
            Self::String(value) => Ok(value),
            _ => Err(conversion("string")),
        }
    }
    pub fn as_string_list(&self) -> Result<&[String]> {
        match self {
            Self::StringList(value) => Ok(value),
            _ => Err(conversion("string list")),
        }
    }
    pub fn as_integer_list(&self) -> Result<&[i32]> {
        match self {
            Self::IntegerList(value) => Ok(value),
            _ => Err(conversion("integer list")),
        }
    }
    pub fn as_float_list(&self) -> Result<&[f64]> {
        match self {
            Self::FloatList(value) => Ok(value),
            _ => Err(conversion("floating-point list")),
        }
    }

    pub fn to_i64(&self) -> Result<i64> {
        match self {
            Self::Integer(value) => Ok(*value),
            _ => Err(conversion("integer")),
        }
    }
    integer_conversion! {
        to_i8: i8, to_u8: u8, to_i16: i16, to_u16: u16,
        to_i32: i32, to_u32: u32, to_u64: u64,
        to_isize: isize, to_usize: usize,
    }

    /// Integer-to-float conversion may round, just as in the source.
    pub fn to_f64(&self) -> Result<f64> {
        match self {
            Self::Integer(value) => Ok(*value as f64),
            Self::Float(value) => Ok(*value),
            _ => Err(conversion("floating-point number")),
        }
    }
    /// Preserve explicit NaN/infinity; reject overflow of a finite input.
    pub fn to_f32(&self) -> Result<f32> {
        let value = self.to_f64()?;
        let result = match self {
            Self::Integer(value) => *value as f32,
            _ => value as f32,
        };
        if value.is_finite() && !result.is_finite() {
            Err(conversion("finite f32"))
        } else {
            Ok(result)
        }
    }
    /// Only the exact, untrimmed strings `true` and `false` are booleans.
    pub fn to_bool(&self) -> Result<bool> {
        match self.as_str()? {
            "true" => Ok(true),
            "false" => Ok(false),
            _ => Err(conversion("boolean")),
        }
    }
    /// Safe counterpart of the source nullable character pointer.
    /// The complete string is borrowed, including any embedded NUL character.
    pub fn to_char(&self) -> Result<Option<&str>> {
        if self.is_empty() {
            Ok(None)
        } else {
            self.as_str().map(Some)
        }
    }

    pub fn to_string_vector(&self) -> Result<Vec<String>> {
        self.as_string_list()?;
        match self.checked_clone()? {
            Self::StringList(value) => Ok(value),
            _ => unreachable!(),
        }
    }
    pub fn to_int_vector(&self) -> Result<Vec<i32>> {
        self.as_integer_list()?;
        match self.checked_clone()? {
            Self::IntegerList(value) => Ok(value),
            _ => unreachable!(),
        }
    }
    pub fn to_double_vector(&self) -> Result<Vec<f64>> {
        self.as_float_list()?;
        match self.checked_clone()? {
            Self::FloatList(value) => Ok(value),
            _ => unreachable!(),
        }
    }
    /// Bounded owned copy; ordinary Rust moves/assignments need no conversion.
    pub fn checked_clone(&self) -> Result<Self> {
        let mut work = ParamWork::default();
        let bytes = self.measure(&mut work)?;
        work.copy(bytes)?;
        Ok(self.clone())
    }

    /// Source `toString`: full precision by default in C++; pass `true` here.
    /// Lists use unquoted `[a, b]` notation, which is deliberately not a codec.
    pub fn to_text(&self, full_precision: bool) -> Result<String> {
        self.to_text_with_work(full_precision, &mut ParamWork::default())
    }

    pub(crate) fn to_text_with_work(
        &self,
        full_precision: bool,
        work: &mut ParamWork,
    ) -> Result<String> {
        self.render(|x| format_float(x, full_precision), work)
    }

    /// Source output-stream operator with the default classic-locale stream
    /// settings (six significant digits). Unlike `to_text`, `5.0` becomes `5`.
    pub fn to_stream_text(&self) -> Result<String> {
        self.to_stream_text_with_work(&mut ParamWork::default())
    }

    pub(crate) fn to_stream_text_with_work(&self, work: &mut ParamWork) -> Result<String> {
        self.render(format_stream_float, work)
    }

    fn render(&self, float: impl Fn(f64) -> String, work: &mut ParamWork) -> Result<String> {
        self.measure(work)?;
        let float_count = match self {
            Self::Float(_) => 1,
            Self::FloatList(values) => values.len(),
            _ => 0,
        };
        // Includes short-lived numeric formatting buffers, not just the final
        // output. Charge the whole list before formatting its first element.
        work.copy(float_count.checked_mul(256).ok_or_else(limit)?)?;
        let capacity = match self {
            Self::Empty => 0,
            Self::String(s) => s.len(),
            Self::Integer(_) | Self::Float(_) => 64,
            Self::StringList(values) => values.iter().try_fold(2usize, |n, s| {
                n.checked_add(s.len())
                    .and_then(|n| n.checked_add(2))
                    .ok_or_else(limit)
            })?,
            Self::IntegerList(values) => values
                .len()
                .checked_mul(16)
                .and_then(|n| n.checked_add(2))
                .ok_or_else(limit)?,
            Self::FloatList(values) => values
                .len()
                .checked_mul(66)
                .and_then(|n| n.checked_add(2))
                .ok_or_else(limit)?,
        };
        work.copy(capacity)?;
        let mut output = String::new();
        output.try_reserve_exact(capacity).map_err(|_| limit())?;
        match self {
            Self::Empty => {}
            Self::String(s) => output.push_str(s),
            Self::Integer(value) => output.push_str(&value.to_string()),
            Self::Float(value) => output.push_str(&float(*value)),
            Self::StringList(values) => {
                output.push('[');
                for (i, value) in values.iter().enumerate() {
                    if i != 0 {
                        output.push_str(", ");
                    }
                    output.push_str(value);
                }
                output.push(']');
            }
            Self::IntegerList(values) => {
                output.push('[');
                for (i, value) in values.iter().enumerate() {
                    if i != 0 {
                        output.push_str(", ");
                    }
                    output.push_str(&value.to_string());
                }
                output.push(']');
            }
            Self::FloatList(values) => {
                output.push('[');
                for (i, value) in values.iter().enumerate() {
                    if i != 0 {
                        output.push_str(", ");
                    }
                    output.push_str(&float(*value));
                }
                output.push(']');
            }
        }
        Ok(output)
    }

    /// Source `<`: unlike types are false; lists compare only their lengths.
    pub fn source_less(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::String(a), Self::String(b)) => a < b,
            (Self::Integer(a), Self::Integer(b)) => a < b,
            (Self::Float(a), Self::Float(b)) => a < b,
            (Self::StringList(a), Self::StringList(b)) => a.len() < b.len(),
            (Self::IntegerList(a), Self::IntegerList(b)) => a.len() < b.len(),
            (Self::FloatList(a), Self::FloatList(b)) => a.len() < b.len(),
            _ => false,
        }
    }
    /// Source `>`: unlike types are false; lists compare only their lengths.
    pub fn source_greater(&self, other: &Self) -> bool {
        other.source_less(self)
    }

    /// Source FNV/combine hash on a 64-bit little-endian SDK host. This explicit
    /// method is bounded; the standard `Hash` trait follows normal Rust semantics.
    pub fn source_hash64(&self) -> Result<u64> {
        self.measure(&mut ParamWork::default())?;
        Ok(self.hash64())
    }

    fn hash64(&self) -> u64 {
        let mut seed = fnv(&[self.value_type() as u8]);
        match self {
            Self::Empty => {}
            Self::String(s) => combine(&mut seed, fnv(s.as_bytes())),
            Self::Integer(x) => combine(&mut seed, fnv(&x.to_le_bytes())),
            Self::Float(x) => combine(&mut seed, hash_float(*x)),
            Self::StringList(values) => {
                combine(&mut seed, fnv(&(values.len() as u64).to_le_bytes()));
                for s in values {
                    combine(&mut seed, fnv(s.as_bytes()));
                }
            }
            Self::IntegerList(values) => {
                combine(&mut seed, fnv(&(values.len() as u64).to_le_bytes()));
                for x in values {
                    combine(&mut seed, fnv(&x.to_le_bytes()));
                }
            }
            Self::FloatList(values) => {
                combine(&mut seed, fnv(&(values.len() as u64).to_le_bytes()));
                for x in values {
                    combine(&mut seed, hash_float(*x));
                }
            }
        }
        seed
    }

    /// Charge complete logical contents before comparison, hashing or copying.
    pub(crate) fn measure(&self, work: &mut ParamWork) -> Result<usize> {
        let mut bytes = std::mem::size_of::<Self>();
        work.consume(bytes)?;
        let mut add = |n: usize| -> Result<()> {
            bytes = bytes
                .checked_add(n)
                .filter(|n| *n <= MAX_VALUE_BYTES)
                .ok_or_else(limit)?;
            work.consume(n)
        };
        match self {
            Self::String(s) => add(s.len())?,
            Self::StringList(values) => {
                if values.len() > MAX_VALUE_ELEMENTS {
                    return Err(limit());
                }
                add(values
                    .len()
                    .checked_mul(std::mem::size_of::<String>())
                    .ok_or_else(limit)?)?;
                for s in values {
                    add(s.len())?;
                }
            }
            Self::IntegerList(values) => {
                if values.len() > MAX_VALUE_ELEMENTS {
                    return Err(limit());
                }
                add(values
                    .len()
                    .checked_mul(std::mem::size_of::<i32>())
                    .ok_or_else(limit)?)?;
            }
            Self::FloatList(values) => {
                if values.len() > MAX_VALUE_ELEMENTS {
                    return Err(limit());
                }
                add(values
                    .len()
                    .checked_mul(std::mem::size_of::<f64>())
                    .ok_or_else(limit)?)?;
            }
            _ => {}
        }
        Ok(bytes)
    }
}

impl PartialOrd for ParamValue {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let order = match (self, other) {
            (Self::Empty, Self::Empty) => return Some(Ordering::Equal),
            (Self::String(a), Self::String(b)) => return Some(a.cmp(b)),
            (Self::Integer(a), Self::Integer(b)) => return Some(a.cmp(b)),
            (Self::Float(a), Self::Float(b)) => return a.partial_cmp(b),
            (Self::StringList(a), Self::StringList(b)) => a.len().cmp(&b.len()),
            (Self::IntegerList(a), Self::IntegerList(b)) => a.len().cmp(&b.len()),
            (Self::FloatList(a), Self::FloatList(b)) => a.len().cmp(&b.len()),
            _ => return None,
        };
        if order == Ordering::Equal && self != other {
            None
        } else {
            Some(order)
        }
    }
}
impl Hash for ParamValue {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u64(self.hash64());
    }
}

macro_rules! from_integer {
    ($($type:ty),* $(,)?) => {$ (
        impl From<$type> for ParamValue {
            fn from(value: $type) -> Self { Self::Integer(i64::from(value)) }
        }
    )*};
}
from_integer!(i8, u8, i16, u16, i32, u32);
impl From<i64> for ParamValue {
    fn from(value: i64) -> Self {
        Self::Integer(value)
    }
}
impl TryFrom<u64> for ParamValue {
    type Error = Error;
    fn try_from(value: u64) -> Result<Self> {
        i64::try_from(value)
            .map(Self::Integer)
            .map_err(|_| conversion("signed i64 storage"))
    }
}
impl TryFrom<usize> for ParamValue {
    type Error = Error;
    fn try_from(value: usize) -> Result<Self> {
        i64::try_from(value)
            .map(Self::Integer)
            .map_err(|_| conversion("signed i64 storage"))
    }
}
impl TryFrom<isize> for ParamValue {
    type Error = Error;
    fn try_from(value: isize) -> Result<Self> {
        i64::try_from(value)
            .map(Self::Integer)
            .map_err(|_| conversion("signed i64 storage"))
    }
}
impl From<f64> for ParamValue {
    fn from(value: f64) -> Self {
        Self::Float(value)
    }
}
impl From<f32> for ParamValue {
    fn from(value: f32) -> Self {
        Self::Float(f64::from(value))
    }
}
impl From<String> for ParamValue {
    fn from(value: String) -> Self {
        Self::String(value)
    }
}
impl From<&str> for ParamValue {
    fn from(value: &str) -> Self {
        Self::String(value.into())
    }
}
impl From<Vec<String>> for ParamValue {
    fn from(value: Vec<String>) -> Self {
        Self::StringList(value)
    }
}
impl From<Vec<i32>> for ParamValue {
    fn from(value: Vec<i32>) -> Self {
        Self::IntegerList(value)
    }
}
impl From<Vec<f64>> for ParamValue {
    fn from(value: Vec<f64>) -> Self {
        Self::FloatList(value)
    }
}

macro_rules! try_into_scalar {
    ($($type:ty => $method:ident),* $(,)?) => {$ (
        impl TryFrom<&ParamValue> for $type {
            type Error = Error;
            fn try_from(value: &ParamValue) -> Result<Self> { value.$method() }
        }
    )*};
}
try_into_scalar! {
    i8 => to_i8, u8 => to_u8, i16 => to_i16, u16 => to_u16,
    i32 => to_i32, u32 => to_u32, i64 => to_i64, u64 => to_u64,
    isize => to_isize, usize => to_usize, f32 => to_f32, f64 => to_f64,
    bool => to_bool, Vec<String> => to_string_vector,
    Vec<i32> => to_int_vector, Vec<f64> => to_double_vector,
}
impl TryFrom<&ParamValue> for String {
    type Error = Error;
    fn try_from(value: &ParamValue) -> Result<Self> {
        value.as_str()?;
        match value.checked_clone()? {
            ParamValue::String(value) => Ok(value),
            _ => unreachable!(),
        }
    }
}

fn fnv(bytes: &[u8]) -> u64 {
    bytes.iter().fold(14_695_981_039_346_656_037, |h, b| {
        (h ^ u64::from(*b)).wrapping_mul(1_099_511_628_211)
    })
}
fn combine(seed: &mut u64, value: u64) {
    *seed ^= value
        .wrapping_add(0x9e3779b97f4a7c15)
        .wrapping_add(*seed << 6)
        .wrapping_add(*seed >> 2);
}
fn hash_float(value: f64) -> u64 {
    fnv(&(if value == 0.0 { 0.0 } else { value }).to_le_bytes())
}

/// `StringUtils::toStr(double, full_precision)`: the source's one
/// `NumericFormatting::appendNumeric<double>`, which
/// [`crate::concept::numeric_formatting`] ports (15 fraction digits, or 3 in
/// the reduced form with its three-digit scientific mantissa).
pub(crate) fn format_float(value: f64, full_precision: bool) -> String {
    let mut text = String::new();
    crate::concept::numeric_formatting::append_numeric_f64(
        value,
        &mut text,
        if full_precision { 15 } else { 3 },
        !full_precision,
    );
    text
}

/// `StringUtils::toStr(float, full_precision)`, `appendNumeric<float>`: six
/// fraction digits in the full form and the `float` threshold literals; see
/// [`format_float`].
pub(crate) fn format_float32(value: f32, full_precision: bool) -> String {
    let mut text = String::new();
    crate::concept::numeric_formatting::append_numeric_f32(
        value,
        &mut text,
        if full_precision { 6 } else { 3 },
        !full_precision,
    );
    text
}

/// A `double` inserted into a default `std::ostream`: `%g` at the stream's
/// default precision of six, with glibc's `nan`/`-nan` and `-0`.
fn format_stream_float(value: f64) -> String {
    crate::concept::numeric_formatting::ostream_g(value, 6)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formatting_shares_work_and_allocation_across_values() {
        let value = ParamValue::from(47.11);
        let mut work = ParamWork {
            remaining_work: 500,
            remaining_bytes: 1024,
        };
        assert_eq!(
            value.to_text_with_work(true, &mut work).unwrap(),
            "47.109999999999999"
        );
        assert!(value.to_text_with_work(true, &mut work).is_err());
        assert_eq!(value.to_f64().unwrap(), 47.11);
        let mut work = ParamWork {
            remaining_work: 10_000,
            remaining_bytes: 300,
        };
        assert!(value.to_text_with_work(true, &mut work).is_err());
    }

    #[test]
    fn float_text_is_the_one_append_numeric_with_its_even_digit_ties() {
        // The exact shortest-round-trip ties of the executed differential in
        // tests/file_info_text_format.rs (oracle rows by bit pattern): where
        // Rust's `{:e}` alone rounds up (`...813e05`), `std::to_chars` takes
        // the even digit. The text of a ParamValue is that routine's.
        assert_eq!(
            format_float(f64::from_bits(0xc100_800d_a680_0000), true),
            "-1.3516970629882812e05"
        );
        assert_eq!(
            format_float(f64::from_bits(0x4315_50f7_dca7_0001), true),
            "1.5000000000000002e15"
        );
        assert_eq!(
            format_float32(f32::from_bits(0x4a06_7a05), true),
            "2.2032652e06"
        );
        assert_eq!(
            format_float32(f32::from_bits(0x4a18_9681), true),
            "2.5000002e06"
        );
        assert_eq!(
            ParamValue::from(f64::from_bits(0xc100_800d_a680_0000))
                .to_text(true)
                .unwrap(),
            "-1.3516970629882812e05"
        );
        // The stream form: six significant digits, glibc's spelling of -0.
        assert_eq!(format_stream_float(-0.0), "-0");
        assert_eq!(format_stream_float(5.0), "5");
        assert_eq!(format_stream_float(1.0e20), "1e+20");
        assert_eq!(format_stream_float(0.0001234567), "0.000123457");
    }

    #[test]
    fn f32_formatter_retains_source_precision_and_binary32_thresholds() {
        assert_eq!(format_float32(47.11, true), "47.110001");
        assert_eq!(format_float32(47.11, false), "47.11");
        assert_eq!(format_float32(0.01, true), "0.01");
        assert_eq!(format_float32(10000.0, true), "1.0e04");
        assert_eq!(format_float32(f32::from_bits(1), true), "1.0e-45");
        assert_eq!(format_float32(f32::INFINITY, true), "inf");
        assert_eq!(format_float32(-0.0, true), "-0.0");
    }
}
