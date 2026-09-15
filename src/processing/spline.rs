// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Spline interpolation, least-squares B-spline smoothing and spline extremum
//! search, ported from the OpenMS `MATH/MISC` splines.
//!
//! | C++ header | Rust module | Support document |
//! |---|---|---|
//! | `CubicSpline2d.h` | [`cubic`](crate::processing::spline::cubic) | `docs/CUBIC_SPLINE2D_SUPPORT.md` |
//! | `BSpline2d.h` | [`b_spline`](crate::processing::spline::b_spline) | `docs/BSPLINE2D_SUPPORT.md` |
//! | `BSplineSmoothingSpline.h` | [`smoothing`](crate::processing::spline::smoothing) | `docs/BSPLINE_SMOOTHING_SPLINE_SUPPORT.md` |
//! | `SplineBisection.h` | [`bisection`](crate::processing::spline::bisection) | `docs/SPLINE_BISECTION_SUPPORT.md` |
//!
//! The two spline families answer an out-of-domain query differently, and the
//! difference is inherited from the source rather than chosen here:
//!
//! * [`CubicSpline2d`](crate::processing::spline::CubicSpline2d) interpolates on
//!   the closed knot interval only and rejects anything outside it, as the C++
//!   throws `Exception::IllegalArgument`.
//! * [`BSpline2d`](crate::processing::spline::BSpline2d) extrapolates: its basis
//!   functions vanish more than two node intervals from their node, so beyond
//!   that distance the evaluation reduces to the fitted mean of the ordinates.
//!
//! Both are exercised against an executed probe of the pinned C++ sources; the
//! support documents record the evidence tier for every claim.

/// Cubic B-spline least-squares smoother over a uniform node grid (`BSpline2d.h`).
pub mod b_spline;
/// Bisection search for the maximum of a spline (`SplineBisection.h`).
pub mod bisection;
/// Natural cubic-spline interpolation of a 2D data set (`CubicSpline2d.h`).
pub mod cubic;
/// Smoothing spline balancing residual sum of squares against model size
/// (`BSplineSmoothingSpline.h`).
pub mod smoothing;

pub use b_spline::{BSpline2d, BoundaryCondition};
pub use bisection::{SplineFunction, spline_bisection};
pub use cubic::{CubicSpline2d, CubicSpline2dFitter};
pub use smoothing::BSplineSmoothingSpline;
