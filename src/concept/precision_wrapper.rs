// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Output of a floating-point value with appropriate precision:
//! `CONCEPT/PrecisionWrapper.h`.
//!
//! # API mapping
//!
//! | Source (`PrecisionWrapper.h` at core `bc9cc12`) | Rust |
//! |---|---|
//! | `template <typename FloatingPointType> struct PrecisionWrapper` | [`PrecisionWrapper`](crate::concept::precision_wrapper::PrecisionWrapper) |
//! | `explicit PrecisionWrapper(const FloatingPointType rhs)` | [`PrecisionWrapper::new`](crate::concept::precision_wrapper::PrecisionWrapper::new) |
//! | `PrecisionWrapper(const PrecisionWrapper& rhs)` | `Clone` and `Copy` |
//! | `PrecisionWrapper()`, private and "intentionally not implemented" | no `Default` |
//! | `FloatingPointType const ref_` | [`PrecisionWrapper::value`](crate::concept::precision_wrapper::PrecisionWrapper::value) |
//! | `precisionWrapper(const FloatingPointType rhs)` | [`precision_wrapper`] |
//! | `operator<<(std::ostream&, const PrecisionWrapper<FloatingPointType>&)` | `Display` for `PrecisionWrapper<f64>` and `PrecisionWrapper<f32>` |
//!
//! The source writes the value through `StringUtils::toStr(value, true)`,
//! which is `NumericFormatting::appendNumeric(value, target,
//! writtenDigits<T>(), false)`: [`crate::concept::numeric_formatting`] with 15
//! fraction digits for `double` and 6 for `float`. The same text is
//! `format::file_info::text_format::to_str` and `to_str_f32`.
//!
//! The source's purpose is to "set the appropriate precision for output
//! temporarily", restoring the stream's precision afterwards "so that no side
//! effects remain". Rust formatting keeps no precision between items, so there
//! is nothing to restore; a formatter's own precision or width (`{:.2}`,
//! `{:>10}`) is ignored, as the source ignores the stream's. The source's note
//! that the make-function cannot return a `const&` concerns C++ value
//! categories and has no Rust counterpart. `long double` is not covered:
//! Rust has none.

use std::fmt;

use super::numeric_formatting::{append_numeric_f32, append_numeric_f64};
use super::types::written_digits;

/// A floating-point value that writes itself at full precision: source
/// `PrecisionWrapper<FloatingPointType>`.
///
/// # Examples
///
/// The upstream class test's literals:
///
/// ```
/// use openms::concept::precision_wrapper::{precision_wrapper, PrecisionWrapper};
///
/// assert_eq!(precision_wrapper(0.123_456_789_012_345_68).to_string(), "0.123456789012346");
/// assert_eq!(precision_wrapper(0.123_456_79_f32).to_string(), "0.123457");
/// assert_eq!(precision_wrapper(1.0).to_string(), "1.0");
/// assert_eq!(PrecisionWrapper::new(100.5).to_string(), "100.5");
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrecisionWrapper<T> {
    /// The wrapped value: source member `ref_`, a `const` copy of the argument.
    pub value: T,
}

impl<T> PrecisionWrapper<T> {
    /// Wraps `value`: source `explicit PrecisionWrapper(const FloatingPointType rhs)`.
    ///
    /// The source notes that the make-function [`precision_wrapper`] is
    /// normally preferred.
    pub const fn new(value: T) -> Self {
        Self { value }
    }
}

/// Wraps `value` for output with appropriate precision: source make-function
/// `precisionWrapper(rhs)`, which deduces the type from its argument.
pub const fn precision_wrapper<T>(value: T) -> PrecisionWrapper<T> {
    PrecisionWrapper::new(value)
}

impl fmt::Display for PrecisionWrapper<f64> {
    /// Writes `StringUtils::toStr(value, true)` for a `double`.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut text = String::new();
        append_numeric_f64(self.value, &mut text, digits::<f64>(), false);
        formatter.write_str(&text)
    }
}

impl fmt::Display for PrecisionWrapper<f32> {
    /// Writes `StringUtils::toStr(value, true)` for a `float`.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut text = String::new();
        append_numeric_f32(self.value, &mut text, digits::<f32>(), false);
        formatter.write_str(&text)
    }
}

/// `writtenDigits<T>()` as a digit count.
fn digits<T: super::types::WrittenDigits>() -> usize {
    usize::try_from(written_digits::<T>()).unwrap_or(0)
}
