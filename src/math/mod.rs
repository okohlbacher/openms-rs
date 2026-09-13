// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Numerical routines ported from the OpenMS `MATH` domain.
//!
//! The first group here is the distribution fitters of `MATH/STATISTICS`; see
//! [`fitters`](crate::math::fitters).
//!
//! Everything in this module computes in `f64`, as the C++ does: the four
//! fitter headers use `double` throughout and never `float`.

/// Parameter estimation for named distributions.
pub mod fitters;
