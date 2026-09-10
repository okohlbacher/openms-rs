// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Typed metadata values, controlled vocabulary terms and acquisition settings.
//! See `docs/METADATA_SUPPORT.md` for the supported source subset and policies.

mod acquisition;
mod value;

pub use acquisition::*;
pub use value::*;
