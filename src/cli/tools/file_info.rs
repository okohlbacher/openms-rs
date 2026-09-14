// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Shows basic information about the file, such as data ranges and file type.
//!
//! Will port `OpenMS4-topp/src/FileInfo.cpp` (topp `174b576`) as a thin wrapper
//! over `crate::format::file_info`. This is a registration stub of the early
//! TOPP bundle's wave 3a: package A5-FILEINFO-TOOL replaces it with the tool.
//! Until then registration fails, so every invocation, `--help` and
//! `-write_ini` included, ends with a not-implemented message and
//! `ILLEGAL_PARAMETERS`, the source's initialisation-failure exit code.

use crate::cli::{ExitCode, Tool, ToolContext, ToolSpec};
use crate::{Error, Result};

/// The `FileInfo` TOPP tool. Not implemented yet: every run fails before
/// parsing its arguments.
pub struct FileInfo;

/// The error every entry point of the stub returns.
fn not_implemented() -> Error {
    Error::Unsupported(
        "FileInfo is not implemented yet; package A5-FILEINFO-TOOL of the early TOPP bundle ports it"
            .into(),
    )
}

impl Tool for FileInfo {
    const NAME: &'static str = "FileInfo";
    const DESCRIPTION: &'static str =
        "Shows basic information about the file, such as data ranges and file type.";

    fn register(_spec: &mut ToolSpec) -> Result<()> {
        Err(not_implemented())
    }

    fn run(_ctx: &ToolContext) -> Result<ExitCode> {
        Err(not_implemented())
    }
}
