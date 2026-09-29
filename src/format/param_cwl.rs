// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The Common Workflow Language writer `ParamCWLFile`
//! (`OpenMS/FORMAT/ParamCWLFile.h`, core `bc9cc12`), as the Release build has
//! it.
//!
//! The source class exports a parameter tree and a tool's description as a
//! CWL `CommandLineTool` document (`<tool>.cwl`), which the TOPP tools write
//! for `-write_cwl` and `-write_nested_cwl`. It has no reader. The document is
//! built by the TDL library (`src/openms/extern/tool_description_lib`,
//! `convertToCWL`), and only when OpenMS is configured with `-DENABLE_TDL=ON`.
//! The Release build at the pins is configured without it: its
//! `writeCWLToStream` throws a `std::runtime_error` with
//! [`TDL_UNAVAILABLE`](crate::format::param_json::TDL_UNAVAILABLE) (`ParamCWLFile.cpp:327-332`), so no CWL document
//! exists to compare with, and every TOPP tool's `-write_cwl` ends with that
//! message and exit 12 (`../oracle/toppbase-completion`, cases `write_cwl`,
//! `write_nested_cwl`, `write_cwl_cwd` and `write_nested_cwl_and_json`, in
//! `tests/topp_cli_completion.rs`). [`ParamCwlFile`](crate::format::param_cwl::ParamCwlFile) is that class: both
//! members refuse with the same message, which is also the one
//! [`ParamJsonFile`](crate::format::param_json::ParamJsonFile) gives.
//!
//! What the TDL-enabled writer would produce is a backlog item with a named
//! dependency: the document model of `tdl/convertToCWL.h` and
//! `tdl/cwl_v1_2.h`, serialised as YAML (the source links `yaml-cpp`), and an
//! oracle built with `-DENABLE_TDL=ON` to compare it with.

use crate::data_structures::ToolInfo;
use crate::param::Param;
use crate::{Error, Result};
use std::io::Write;

pub use crate::format::param_json::TDL_UNAVAILABLE;

/// Source class `ParamCWLFile`, as the Release build (without TDL) has it:
/// [`store`](Self::store) and [`write_cwl_to_stream`](Self::write_cwl_to_stream)
/// refuse with [`TDL_UNAVAILABLE`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ParamCwlFile {
    /// Source `flatHierarchy`: with TDL, every parameter would be listed
    /// without nesting, its name joined from its sections with `__` (the
    /// class documentation's `':'` replaced by `"__"`); `-write_cwl` sets it,
    /// `-write_nested_cwl` does not. Kept for the API; the Release build
    /// writes nothing either way.
    pub flat_hierarchy: bool,
}

impl ParamCwlFile {
    /// Source `store` (`ParamCWLFile.cpp:41-62`): the Release build opens
    /// `filename` (`-` is standard output), throwing
    /// `std::ios::failure("Unable to create file: <filename>")` when it cannot,
    /// and then fails in [`write_cwl_to_stream`](Self::write_cwl_to_stream),
    /// leaving the file it opened empty. This port refuses before it touches
    /// `filename`, as [`ParamJsonFile::store`](crate::format::param_json::ParamJsonFile::store)
    /// does (CPP-350), and so never gives the open's message: the TOPP
    /// command line checks the target with the source's writability check
    /// before it gets here. `tool_info` is unused, as without TDL in the
    /// source.
    ///
    /// # Errors
    ///
    /// Always [`Error::Unsupported`] with [`TDL_UNAVAILABLE`].
    pub fn store(&self, filename: &str, param: &Param, tool_info: &ToolInfo) -> Result<()> {
        let _ = filename;
        self.write_cwl_to_stream(&mut std::io::sink(), param, tool_info)
    }

    /// Source `writeCWLToStream` (`ParamCWLFile.cpp:64-333`): without TDL it
    /// throws `std::runtime_error` with [`TDL_UNAVAILABLE`]
    /// (`ParamCWLFile.cpp:327-332`) and writes nothing.
    ///
    /// # Errors
    ///
    /// Always [`Error::Unsupported`] with [`TDL_UNAVAILABLE`].
    pub fn write_cwl_to_stream(
        &self,
        output: &mut dyn Write,
        param: &Param,
        tool_info: &ToolInfo,
    ) -> Result<()> {
        let _ = (output, param, tool_info, self.flat_hierarchy);
        Err(Error::Unsupported(TDL_UNAVAILABLE.into()))
    }
}
