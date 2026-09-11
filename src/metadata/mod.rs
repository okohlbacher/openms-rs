// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Typed metadata values, controlled vocabulary terms and acquisition settings.
//! See `docs/METADATA_SUPPORT.md` for the supported source subset and policies.

mod acquisition;
mod document_identifier;
mod experiment_support;
mod value;

pub use acquisition::*;
pub use document_identifier::*;
pub use experiment_support::*;
pub use value::*;

mod experimental_design;
pub use experimental_design::*;

mod experiment_values;
pub use experiment_values::*;

mod experimental_settings;
pub use experimental_settings::*;
