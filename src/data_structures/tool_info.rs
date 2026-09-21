// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Tool metadata shared by the parameter-file serialisers: the native form of
//! core `DATASTRUCTURES/ToolInfo.h` (core `bc9cc12`, header-only).
//!
//! `TOPPBase` fills one of these for each tool description it writes
//! (`TOPPBase.cpp:2591-2598`), and the CTD writer (`cli::ParamCtdFile`,
//! behind the `paramxml` feature) prints it as the `<tool>`
//! element's attributes, description and citations.

/// Tool metadata for the tool-description writers (source struct `ToolInfo`).
///
/// The source's members carry a trailing underscore (`version_`, `name_`,
/// …); the fields here drop it.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ToolInfo {
    /// Product version of the tool.
    pub version: String,
    /// Tool name.
    pub name: String,
    /// Documentation URL.
    pub docurl: String,
    /// Category, as the tool registry records it.
    pub category: String,
    /// One-line description.
    pub description: String,
    /// Citation identifiers, in order; `TOPPBase` puts the OpenMS citation's
    /// DOI first and then each tool citation's.
    pub citations: Vec<String>,
}
