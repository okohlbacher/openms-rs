// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Descriptions of TOPP tools and of the external programs a generic wrapper
//! calls: the native form of core `DATASTRUCTURES/ToolDescription.h` and
//! `ToolDescription.cpp` (core `bc9cc12`), whose types live in the source's
//! `OpenMS::Internal` namespace.
//!
//! The tool registry of the `cli` package (`cli::ToolHandler`, behind the
//! `paramxml` feature) keys these by tool name; the
//! legacy `.ttd` registry fills the external details. See
//! `docs/TOPP_CLI_SUPPORT.md` (*Tool registry*) for the mapping of every
//! source member.
//!
//! The source's copy constructors, copy assignments and destructors are
//! `Clone` and ordinary ownership here. Its `operator<` is not a consistent
//! order with its `operator==` — two descriptions that differ only in their
//! category are neither equal nor ordered — so it is
//! [`ToolDescriptionInternal::source_less`] rather than a `PartialOrd`
//! implementation, which Rust requires to agree with `PartialEq`.

use crate::param::Param;
use std::collections::BTreeMap;

/// Maps a file of the external program to the TOPP parameter that names it
/// (source `Internal::FileMapping`).
///
/// A generic wrapper moves the file at `location` to the file the TOPP
/// parameter `target` names.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FileMapping {
    /// A mix of a regular expression and wrapper macros, which the tool
    /// expands (source: "a regex/macro mix; to be expanded by tool").
    pub location: String,
    /// The TOPP parameter that determines the desired name.
    pub target: String,
}

/// The file-name mappings for all input and output files of an external
/// program (source `Internal::MappingParam`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MappingParam {
    /// Command-line fragments by mapping id, ordered by id as the source's
    /// `std::map<Int, String>` is.
    pub mapping: BTreeMap<i32, String>,
    /// Files moved before the program runs, in declaration order.
    pub pre_moves: Vec<FileMapping>,
    /// Files moved after the program has run, in declaration order.
    pub post_moves: Vec<FileMapping>,
}

/// The registry part of a tool description (source
/// `Internal::ToolDescriptionInternal`): whether the tool is internal, its
/// name, its category and its `-type` values.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ToolDescriptionInternal {
    /// `true` for a tool of the OpenMS release, `false` for an external
    /// program wrapped by a generic wrapper. The default is `false`, as the
    /// source's member initialiser.
    pub is_internal: bool,
    /// Tool name.
    pub name: String,
    /// Category, for example `Quantitation`.
    pub category: String,
    /// The tool's `-type` values, in declaration order.
    pub types: Vec<String>,
}

impl ToolDescriptionInternal {
    /// A description with every field given (source constructor
    /// `ToolDescriptionInternal(p_is_internal, p_name, p_category, p_types)`).
    pub fn new(is_internal: bool, name: &str, category: &str, types: &[String]) -> Self {
        Self {
            is_internal,
            name: name.to_owned(),
            category: category.to_owned(),
            types: types.to_vec(),
        }
    }

    /// A description with a name and types only (source constructor
    /// `ToolDescriptionInternal(p_name, p_types)`): not internal, with an
    /// empty category.
    pub fn with_types(name: &str, types: &[String]) -> Self {
        Self::new(false, name, "", types)
    }

    /// The text the source's `operator<` compares: the name, a dot and the
    /// types joined by commas.
    pub fn sort_key(&self) -> String {
        format!("{}.{}", self.name, self.types.join(","))
    }

    /// Source `operator<`: [`sort_key`](Self::sort_key) compared byte-wise.
    ///
    /// The source returns `false` early when both sides are the same object;
    /// comparing equal keys also gives `false`, so no special case is needed.
    pub fn source_less(&self, other: &Self) -> bool {
        self.sort_key() < other.sort_key()
    }
}

/// How a generic wrapper runs one external program, one entry per `-type`
/// (source `Internal::ToolExternalDetails`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ToolExternalDetails {
    /// Text printed when the program starts.
    pub text_startup: String,
    /// Text printed when the program fails.
    pub text_fail: String,
    /// Text printed when the program has finished.
    pub text_finish: String,
    /// The external program's own category.
    pub category: String,
    /// Command-line template of the program.
    pub commandline: String,
    /// File name of the external program.
    pub path: String,
    /// Directory the command is executed from.
    pub working_directory: String,
    /// File-name mappings of the program's inputs and outputs.
    pub tr_table: MappingParam,
    /// Parameters the wrapper registers for this program.
    pub param: Param,
}

/// A tool description for internal and external tools (source
/// `Internal::ToolDescription`, which derives from `ToolDescriptionInternal`).
///
/// The source's base class is the [`internal`](Self::internal) field.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ToolDescription {
    /// Registry part: internal flag, name, category and types.
    pub internal: ToolDescriptionInternal,
    /// Additional details for external tools, one entry for each type.
    pub external_details: Vec<ToolExternalDetails>,
}

impl ToolDescription {
    /// The description of an internal TOPP tool (source constructor
    /// `ToolDescription(p_name, p_category, p_types = StringList())`):
    /// internal, with no external details.
    pub fn new(name: &str, category: &str, types: &[String]) -> Self {
        Self {
            internal: ToolDescriptionInternal::new(true, name, category, types),
            external_details: Vec::new(),
        }
    }
}
