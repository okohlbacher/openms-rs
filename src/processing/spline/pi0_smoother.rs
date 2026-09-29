// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The smoothing spline `MultipleTesting::pi0Est` constructs, supplied to
//! [`crate::math::multiple_testing::pi0_est`].
//!
//! Completes the port of `src/openms/include/OpenMS/MATH/STATISTICS/MultipleTesting.h`
//! (`MultipleTesting.cpp`, `pi0Est`'s `Pi0Method::Smoother` branch). The source
//! builds `BSplineSmoothingSpline spl(xs, ys, -1.0, smooth_df)` and reads
//! `spl.eval(max_lambda)` guarded by `spl.ok()`. `math` may not import
//! `processing` (`tools/check_module_cycles.py`), so the math function states the
//! requirement as the [`Pi0Smoother`](crate::math::multiple_testing::Pi0Smoother) trait and this module satisfies it with the
//! crate's port of that very class, [`BSplineSmoothingSpline`].
//!
//! [`pi0_est_smoothed`](crate::processing::spline::pi0_smoother::pi0_est_smoothed) is the source's `pi0Est` signature exactly — the same
//! five parameters with the same defaults available as constants — so a caller
//! that wants the source's default call no longer has to know about the seam.
//! See `docs/MULTIPLE_TESTING_SUPPORT.md`.

use crate::Result;
use crate::math::multiple_testing::{Pi0Method, Pi0Result, Pi0Smoother, pi0_est};
use crate::processing::spline::smoothing::BSplineSmoothingSpline;

/// The source's smoothing parameter for `pi0Est`: `-1.0`, meaning "choose the
/// residual budget automatically" (`MultipleTesting.cpp`, `pi0Est`).
pub const PI0_SMOOTHING_PARAMETER: f64 = -1.0;

/// [`Pi0Smoother`] backed by the ported `BSplineSmoothingSpline`.
///
/// `smooth_eval(x, y, smooth_df, at)` is the source's
/// `BSplineSmoothingSpline spl(xs, ys, -1.0, smooth_df); spl.eval(at)`. Where the
/// source's `spl.ok()` guard is false this returns `None`, so
/// [`pi0_est`] takes the same `min(min(pi0_lambda), 1)` fallback. The port's
/// constructor reports a failed fit as an `Err` rather than as a spline in a
/// not-ok state, and its `eval` reports a non-finite position or result as an
/// `Err`; both map to `None`.
///
/// # Notes
///
/// With `smooth_log_pi0`, the source feeds `log(pi0_lambda)` to the spline and a
/// non-positive per-lambda estimate becomes `-inf`. The ported spline refuses
/// non-finite ordinates, so such an input takes the fallback here. The source's
/// spline does not check its inputs and attempts the fit through `-inf`; the
/// Linux x86_64 Release build then takes the same fallback (`pi0 = 0`,
/// `pi0_smooth = false`), executed for five such inputs in
/// `tests/multiple_testing.rs::pi0_with_an_empty_lambda_matches_the_release_build`
/// (`docs/MULTIPLE_TESTING_SUPPORT.md`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SplinePi0Smoother;

impl Pi0Smoother for SplinePi0Smoother {
    fn smooth_eval(&self, x: &[f64], y: &[f64], smooth_df: i32, at: f64) -> Option<f64> {
        BSplineSmoothingSpline::with_smoothing(x, y, PI0_SMOOTHING_PARAMETER, smooth_df)
            .ok()
            .and_then(|spline| spline.eval(at).ok())
    }
}

/// Estimate the proportion of true null hypotheses, exactly as the source's
/// `MultipleTesting::pi0Est` does, smoothing spline included.
///
/// This is [`pi0_est`] with [`SplinePi0Smoother`] supplied. The source's
/// defaults are `lambda = {}` (the `0.05, 0.10, ..., 0.95` grid),
/// `pi0_method = Pi0Method::Smoother`,
/// `smooth_df = `[`crate::math::multiple_testing::DEFAULT_SMOOTH_DF`] and
/// `smooth_log_pi0 = false`; Rust has no default arguments, so pass those
/// values for the source's default call. Every behaviour — the repeated-addition
/// lambda grid, the single-lambda short circuit, the five fallbacks, the
/// bootstrap branch — is documented on [`pi0_est`].
///
/// # Errors
///
/// As [`pi0_est`]: [`crate::Error::InvalidValue`] when no finite p-value is
/// supplied, when a p-value is outside `[0, 1]`, or when a lambda is outside
/// `[0, 1)`; [`crate::Error::InvalidRange`] when the input exceeds
/// [`crate::math::multiple_testing::MAX_ITEMS`].
pub fn pi0_est_smoothed(
    p_values: &[f64],
    lambda: &[f64],
    method: Pi0Method,
    smooth_df: i32,
    smooth_log_pi0: bool,
) -> Result<Pi0Result> {
    pi0_est(
        p_values,
        lambda,
        method,
        smooth_df,
        smooth_log_pi0,
        Some(&SplinePi0Smoother),
    )
}
