// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Detects two-dimensional features in LC-MS data.
//!
//! Will port `OpenMS4-topp/src/FeatureFinderCentroided.cpp` (topp `174b576`)
//! over the FeatureFinderAlgorithmPicked port. This is a registration stub
//! of the early TOPP bundle's wave 3a: package C5-FFC-WRAPPER replaces it with
//! the wrapper (registration, loading, error branches, the FAIMS refusal,
//! annotations and output clean-up). Until then registration fails, so every
//! invocation, `--help` and `-write_ini` included, ends with a not-implemented
//! message and `ILLEGAL_PARAMETERS`, the source's initialisation-failure exit
//! code.

use crate::cli::{ExitCode, Tool, ToolContext, ToolSpec};
use crate::{Error, Result};

/// The `FeatureFinderCentroided` TOPP tool. Not implemented yet: every run
/// fails before parsing its arguments.
pub struct FeatureFinderCentroided;

/// The error every entry point of the stub returns.
fn not_implemented() -> Error {
    Error::Unsupported(
        "FeatureFinderCentroided is not implemented yet; package C5-FFC-WRAPPER of the early TOPP bundle ports it"
            .into(),
    )
}

impl Tool for FeatureFinderCentroided {
    const NAME: &'static str = "FeatureFinderCentroided";
    const DESCRIPTION: &'static str = "Detects two-dimensional features in LC-MS data.";

    fn register(_spec: &mut ToolSpec) -> Result<()> {
        Err(not_implemented())
    }

    fn run(_ctx: &ToolContext) -> Result<ExitCode> {
        Err(not_implemented())
    }
}
