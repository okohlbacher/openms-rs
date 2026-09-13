// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Descriptive statistics: summary measures, ranks, binned counts.
//!
//! This is the port of the C++ `OpenMS::Math` statistics namespace, which the
//! rest of the library builds on:
//!
//! - `MATH/StatisticFunctions.h` — the free-function library
//!   ([`crate::math::statistic_functions`]).
//! - `MATH/STATISTICS/BasicStatistics.h` — an accumulating weighted
//!   distribution ([`crate::math::basic_statistics::BasicStatistics`]).
//! - `MATH/STATISTICS/RankData.h` — SciPy-compatible ranking
//!   ([`crate::math::rank_data::rankdata`]).
//! - `MATH/STATISTICS/Histogram.h` — a binned counter
//!   ([`crate::math::histogram::Histogram`]).
//!
//! Every function computes in `f64`, as the source does, and reproduces the
//! source's accumulation order rather than a mathematically equivalent
//! rearrangement: floating-point addition is not associative, so a "cleaner"
//! fold would silently change results that other ported code compares against.
//!
//! Degrees of freedom are **not** uniform in the source and are preserved
//! exactly as written: [`crate::math::statistic_functions::variance`] and
//! [`crate::math::statistic_functions::covariance`] divide by `n - 1`,
//! [`crate::math::statistic_functions::mean_square_error`] and
//! [`crate::math::statistic_functions::mean_absolute_deviation`] divide by `n`,
//! and [`crate::math::basic_statistics::BasicStatistics`] divides by the
//! probability mass. See `docs/STATISTIC_FUNCTIONS_SUPPORT.md`,
//! `docs/BASIC_STATISTICS_SUPPORT.md`, `docs/RANK_DATA_SUPPORT.md` and
//! `docs/HISTOGRAM_SUPPORT.md`.

/// An accumulating weighted distribution with a normal approximation.
pub mod basic_statistics;
/// A binned counter over a closed value range.
pub mod histogram;
/// SciPy-compatible ranking with selectable tie and NaN handling.
pub mod rank_data;
/// Means, medians, quantiles, deviations and correlation coefficients.
pub mod statistic_functions;
