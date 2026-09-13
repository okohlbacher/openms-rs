// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Fitters that estimate the parameters of a named distribution from data.
//!
//! One module per `MATH/STATISTICS` header, plus the Levenberg-Marquardt
//! solver all four of them share. Each fitter keeps its own result type,
//! because the C++ nests a differently-defaulted result struct inside each
//! class, and each keeps the source's initial guess, termination rule and
//! parameter transform - a different starting point on a non-convex surface
//! converges somewhere else.
//!
//! See `docs/DISTRIBUTION_FITTERS_SUPPORT.md`.

/// Least-squares fit of a Gamma distribution to (x, y) points.
pub mod gamma;
/// Least-squares fit of a Gaussian to (x, y) points.
pub mod gauss;
/// Least-squares fit of a Gumbel distribution to (x, y) points.
pub mod gumbel;
/// Weighted maximum-likelihood fit of a Gumbel distribution to raw samples.
pub mod gumbel_max_likelihood;
/// The Levenberg-Marquardt solver the four fitters share.
pub mod levenberg_marquardt;
