// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Finds mass spectrometric peaks in profile mass spectra.
//!
//! Will port `OpenMS4-topp/src/PeakPickerHiRes.cpp` (topp `174b576`). This is a
//! registration stub of the early TOPP bundle's wave 3a: package
//! P3-PICKER-TOOL replaces it with the in-memory tool, its `algorithm`
//! subsection and its parameter-failure and `-write_ini` behaviour. Until then
//! registration fails, so every invocation, `--help` and `-write_ini` included,
//! ends with a not-implemented message and `ILLEGAL_PARAMETERS`, the source's
//! initialisation-failure exit code.

use crate::cli::{ExitCode, Tool, ToolContext, ToolSpec};
use crate::{Error, Result};

/// The `PeakPickerHiRes` TOPP tool. Not implemented yet: every run fails
/// before parsing its arguments.
pub struct PeakPickerHiRes;

/// The error every entry point of the stub returns.
fn not_implemented() -> Error {
    Error::Unsupported(
        "PeakPickerHiRes is not implemented yet; package P3-PICKER-TOOL of the early TOPP bundle ports it"
            .into(),
    )
}

impl Tool for PeakPickerHiRes {
    const NAME: &'static str = "PeakPickerHiRes";
    const DESCRIPTION: &'static str = "Finds mass spectrometric peaks in profile mass spectra.";

    fn register(_spec: &mut ToolSpec) -> Result<()> {
        Err(not_implemented())
    }

    fn run(_ctx: &ToolContext) -> Result<ExitCode> {
        Err(not_implemented())
    }
}
