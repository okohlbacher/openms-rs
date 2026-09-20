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
//! Three of the modules here are not formulas but *substrate*: they are the
//! arithmetic and the ordering the source's own results depend on, and every
//! port that has to reproduce the Release build's exact bits or its exact
//! permutation needs them, not only the one that first measured them.
//!
//! - `crate::math::x86_64` (crate-private) — the SSE2 instructions behind `double`
//!   arithmetic: IEEE 754 fixes every finite result, but not which NaN bit
//!   pattern an operation yields, and the Release build's answer is the
//!   instruction set's, not the host's. A value printed as `-nan` is
//!   `0xfff8000000000000`, and only this module produces it on an arm64 host.
//! - `crate::math::libstdcxx` (crate-private) — `std::__lower_bound` and
//!   `std::__upper_bound`, which decide *positions* rather than values, on keys
//!   the standard's own precondition need not hold for.
//! - [`crate::math::source_sort`] — the permutations `std::sort` and
//!   `std::stable_sort` leave, comparison by comparison and move by move. An
//!   order statistic read positionally out of a sorted range is only defined
//!   once that permutation is, so this is arithmetic's ordering counterpart and
//!   belongs beside it.
//!
//! All three were promoted out of `analysis::feature_finder_picked`, which
//! measured them first; the move changed paths and module documentation only.
//! Decision D16 of `docs/EARLY_TOPP_WORK_PACKAGES.md` records why reproducing
//! an unspecified `std::sort` permutation is in scope at all.
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
/// `std::__lower_bound` and `std::__upper_bound` as the Release build's
/// libstdc++ implements them.
pub(crate) mod libstdcxx;
/// q-values, pi0 estimation and local false discovery rates.
pub mod multiple_testing;
/// An EM-fitted two-component mixture of search-engine scores.
pub mod posterior_error_probability;
/// SciPy-compatible ranking with selectable tie and NaN handling.
pub mod rank_data;
/// The C++ Release build's `std::sort` and `std::stable_sort` permutations.
pub mod source_sort;
/// Means, medians, quantiles, deviations and correlation coefficients.
pub mod statistic_functions;
/// The SSE2 instruction behaviour behind the Release build's `double`
/// arithmetic, where IEEE 754 alone does not fix the result.
pub(crate) mod x86_64;
