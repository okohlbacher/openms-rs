// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Shared constants, runtime logging, progress reporting, and owned unique IDs.

pub mod constants;
pub mod log_stream;
/// General numeric helpers: ppm and Dalton tolerances, rounding, binning,
/// interval transforms and binomial statistics, from `MATH/MathFunctions.h`.
pub mod math_functions;
pub mod progress_logger;
pub mod unique_id;
pub use unique_id::{HasUniqueId, UniqueId, UniqueIdGenerator};
