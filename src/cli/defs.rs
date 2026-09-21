// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Citations of TOPP tools: the native form of the `OpenMS4-cli` package's
//! `APPLICATIONS/TOPPBase_defs.h` (cli `c19e494`, header-only).
//!
//! The header also declares the three parameter exceptions `TOPPBase` throws,
//! `UnregisteredParameter`, `WrongParameterType` and
//! `RequiredParameterNotGiven`. They are not types here:
//!
//! * `RequiredParameterNotGiven` is the validation the lifecycle runs before
//!   the tool body, which reports `Error: The required parameter '<name>' was
//!   not given or is empty!` and ends the run with
//!   [`ExitCode::MissingParameters`](super::ExitCode::MissingParameters), as
//!   the source's catch does (`TOPPBase.cpp:466-474`).
//! * `UnregisteredParameter` and `WrongParameterType` signal a tool reading a
//!   parameter it did not register, or with the wrong accessor. The
//!   [`ToolContext`](super::ToolContext) accessors return
//!   [`Error::InvalidValue`](crate::Error::InvalidValue) for both, and a tool
//!   body that lets one escape ends as the lifecycle maps that error. The
//!   source maps both to `INTERNAL_ERROR` (12); see `docs/TOPP_CLI_SUPPORT.md`
//!   for the exit codes.

use std::fmt;

/// A citation of a TOPP tool (source struct `Citation`), printed by
/// `--help` and written as a DOI to the tool description.
///
/// The fields are static text, so a tool declares its citations as a constant
/// ([`Tool::CITATIONS`](super::Tool::CITATIONS)). The source suggests the AMA
/// style for the authors and the `when_where` part.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Citation {
    /// Authors in AMA style: `surname initials`, comma-separated.
    pub authors: &'static str,
    /// Title of the article.
    pub title: &'static str,
    /// Suggested format: `journal. year; volume, issue: pages`.
    pub when_where: &'static str,
    /// A plain DOI, without a URL (for example `10.1021/pr100177k`). The
    /// source's own tools also put URLs here, which are printed as given.
    pub doi: &'static str,
}

impl Citation {
    /// The rendered citation, as the source `toString`:
    /// `<authors>. <title>. <when_where>. doi:<doi>.`
    pub fn to_source_string(&self) -> String {
        format!(
            "{}. {}. {}. doi:{}.",
            self.authors, self.title, self.when_where, self.doi
        )
    }
}

impl fmt::Display for Citation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.to_source_string())
    }
}

/// The OpenMS citation every tool prints (source `TOPPBase::cite_openms`,
/// `TOPPBase.cpp:79-81`).
pub const CITE_OPENMS: Citation = Citation {
    authors: "Pfeuffer, J., Bielow, C., Wein, S. et al.",
    title: "OpenMS 3 enables reproducible analysis of large-scale mass spectrometry data",
    when_where: "Nat Methods (2024)",
    doi: "10.1038/s41592-024-02197-7",
};
