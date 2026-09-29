// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of the `ML` domain that need nothing above `math`: clustering,
//! RANSAC, regression, interpolation and the other generic learners.
//!
//! Submodules are registered below. An `ML` port that needs kernel or format
//! types belongs in the analysis domain that uses it instead, because `math`
//! is a sink in the module graph and must stay one.

/// `LinearRegression`: least-squares line fits with goodness-of-fit statistics.
pub mod linear_regression;
/// `LinearRegressionWithoutIntercept`: least-squares line through the origin.
pub mod linear_regression_without_intercept;
/// `QuadraticRegression`: (weighted) least-squares quadratic fits.
pub mod quadratic_regression;
/// `RANSAC`, `RANSACModel`, `RANSACModelLinear`, `RANSACModelQuadratic`: outlier detection.
pub mod ransac;
