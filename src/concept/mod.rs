// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Shared constants, runtime logging, progress reporting, and owned unique IDs.

pub mod constants;
pub mod log_stream;
/// Data parallelism with a determinism contract: a parallel result must be
/// bit-identical to the serial one.
pub mod parallel;
pub mod progress_logger;
pub mod unique_id;
pub use unique_id::{HasUniqueId, UniqueId, UniqueIdGenerator};
