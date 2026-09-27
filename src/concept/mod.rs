// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Shared constants, runtime logging, progress reporting, and owned unique IDs.

/// Boost.Regex-compatible regular expressions over `fancy-regex`.
pub mod boost_regex;
pub mod constants;
/// Fuzzy comparison of text that tolerates numeric differences, from
/// `CONCEPT/FuzzyStringComparator.h`: the comparator behind `FuzzyDiff`.
pub mod fuzzy_string_comparator;
pub mod log_stream;
/// General numeric helpers: ppm and Dalton tolerances, rounding, binning,
/// interval transforms and binomial statistics, from `MATH/MathFunctions.h`.
pub mod math_functions;
/// Data parallelism with a determinism contract: a parallel result must be
/// bit-identical to the serial one.
pub mod parallel;
pub mod progress_logger;
pub mod unique_id;
pub use unique_id::{HasUniqueId, UniqueId, UniqueIdGenerator};
/// Enums shared across the library, from `CONCEPT/CommonEnums.h`: plotting-axis
/// units and m/z tolerance units with their name tables.
pub mod common_enums;
/// FNV-1a hashing and golden-ratio hash combining, from `CONCEPT/HashUtils.h`.
pub mod hash_utils;
/// Named critical sections and the condition macros, from `CONCEPT/Macros.h`.
pub mod macros;
/// Locale-independent decimal text for floating-point values, from
/// `CONCEPT/Detail/NumericFormatting.h`: the rules behind every
/// `StringUtils::toStr` of a float.
pub mod numeric_formatting;
/// Output of a floating-point value with appropriate precision, from
/// `CONCEPT/PrecisionWrapper.h`.
pub mod precision_wrapper;
/// Basic type aliases, ASCII codes and output precisions, from `CONCEPT/Types.h`.
pub mod types;
