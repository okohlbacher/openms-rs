// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Numerical routines ported from the OpenMS `MATH` domain.
//!
//! - `MATH/StatisticFunctions.h` — the free-function library
//!   ([`crate::math::statistic_functions`]).
//! - `MATH/STATISTICS/BasicStatistics.h` — an accumulating weighted
//!   distribution ([`crate::math::basic_statistics::BasicStatistics`]).
//! - `MATH/STATISTICS/RankData.h` — SciPy-compatible ranking
//!   ([`crate::math::rank_data::rankdata`]).
//! - `MATH/STATISTICS/Histogram.h` — a binned counter
//!   ([`crate::math::histogram::Histogram`]).
//! - `MATH/STATISTICS/{Gauss,GammaDistribution,GumbelDistribution,GumbelMaxLikelihood}Fitter.h`
//!   — parameter estimation for named distributions ([`crate::math::fitters`]).
//! - `MATH/STATISTICS/KernelDensityEstimation.h` — FFT-based kernel density
//!   estimation ([`crate::math::kernel_density`]), on the radix-2 transform in
//!   [`crate::math::fft`] that replaces the vendored evergreen FFT.
//! - `MATH/STATISTICS/MultipleTesting.h` — q-values, pi0 and local FDR
//!   ([`crate::math::multiple_testing`]).
//! - `MATH/STATISTICS/PosteriorErrorProbabilityModel.h` — the EM-fitted score
//!   mixture ([`crate::math::posterior_error_probability`]).
//!
//! Everything here computes in `f64`, as the source does — the fitter headers
//! use `double` throughout and never `float` — and reproduces the source's
//! accumulation order rather than a mathematically equivalent rearrangement:
//! floating-point addition is not associative, so a "cleaner" fold would
//! silently change results that other ported code compares against.
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
/// A radix-2 decimation-in-frequency FFT and its packed real transform.
pub mod fft;
/// Parameter estimation for named distributions.
pub mod fitters;
/// A binned counter over a closed value range.
pub mod histogram;
/// FFT-based Gaussian kernel density estimation.
pub mod kernel_density;
/// q-values, pi0 estimation and local false discovery rates.
pub mod multiple_testing;
/// An EM-fitted two-component mixture of search-engine scores.
pub mod posterior_error_probability;
/// SciPy-compatible ranking with selectable tie and NaN handling.
pub mod rank_data;
/// Means, medians, quantiles, deviations and correlation coefficients.
pub mod statistic_functions;
