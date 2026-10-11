// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Citations of TOPP tools and the three `TOPPBase` parameter exceptions: the
//! native form of the `OpenMS4-cli` package's `APPLICATIONS/TOPPBase_defs.h`
//! (cli `c19e494`, header-only).
//!
//! The header declares `struct Citation` and three exception classes,
//! `UnregisteredParameter`, `WrongParameterType` and
//! `RequiredParameterNotGiven`, each constructed from the throwing location
//! and a `parameter` text that becomes its `what()`. They are not types here,
//! because a native tool body does not throw: it returns a
//! [`ToolError`](super::ToolError), and what the source's exception *does* —
//! the catch block of `TOPPBase::main` that reports it and the exit code that
//! block returns — is what a tool needs. Each of the three is therefore a
//! [`ToolError`] constructor carrying its catch block's text and code
//! ([`ToolError::required_parameter_not_given`],
//! [`ToolError::unregistered_parameter`] and
//! [`ToolError::wrong_parameter_type`]), for the TOPP tools that throw them
//! themselves (at topp `174b576`: `IDSplitter.cpp:98`,
//! `IsobaricWorkflow.cpp:508`, `LuciphorAdapter.cpp:447` and
//! `ProteinQuantifier.cpp:1140` throw `RequiredParameterNotGiven`,
//! `MapAlignerIdentification.cpp:282` and `:578` `WrongParameterType`). The
//! framework raises them where the source does:
//!
//! * `RequiredParameterNotGiven` is the validation the lifecycle runs before
//!   the tool body, which reports `Error: The required parameter '<name>' was
//!   not given or is empty!` and ends the run with
//!   [`ExitCode::MissingParameters`](super::ExitCode::MissingParameters), as
//!   the source's catch does (`TOPPBase.cpp:466-474`).
//! * `UnregisteredParameter` and `WrongParameterType` signal a tool reading a
//!   parameter it did not register, or with the wrong accessor. The
//!   [`ToolContext`](super::ToolContext) accessors return
//!   [`Error::InvalidValue`](crate::Error::InvalidValue) for both and record
//!   the refusal, and a tool body that lets one escape ends as the source's
//!   catch blocks end it (`TOPPBase.cpp:482-493`): `Internal error: Request
//!   for unregistered parameter '<name>'` or `Internal error: Request for
//!   parameter with wrong type '<name>'`, and `INTERNAL_ERROR` (12); see
//!   `docs/TOPP_CLI_SUPPORT.md` for the exit codes.
//!
//! Each source constructor also hands its `what()` to
//! `GlobalExceptionHandler::setMessage`, which only the source's terminate
//! handler reads, for an exception nothing catches. `TOPPBase::main` catches
//! all three, and a native tool has no terminate handler, so that side effect
//! has no counterpart. The source's catch blocks also write `Error occurred in
//! line <n> of file <f> (in function: <fn>) !` to the `-log` file at debug
//! level 1; a native error carries no source location, so that line is not
//! written (`docs/TOPP_CLI_SUPPORT.md`, *Native differences*).

use super::outcome::ToolError;
use super::parameter::ExitCode;
use std::fmt;

/// The `RequiredParameterNotGiven` catch block's line for `what`
/// (`TOPPBase.cpp:466-474`): `what` is quoted unless it already starts with a
/// single quote.
pub(crate) fn required_parameter_text(what: &str) -> String {
    if what.starts_with('\'') {
        format!("Error: The required parameter {what} was not given or is empty!")
    } else {
        format!("Error: The required parameter '{what}' was not given or is empty!")
    }
}

/// The `UnregisteredParameter` catch block's line for `what`
/// (`TOPPBase.cpp:482-487`).
pub(crate) fn unregistered_parameter_text(what: &str) -> String {
    format!("Internal error: Request for unregistered parameter '{what}'")
}

/// The `WrongParameterType` catch block's line for `what`
/// (`TOPPBase.cpp:488-493`).
pub(crate) fn wrong_parameter_type_text(what: &str) -> String {
    format!("Internal error: Request for parameter with wrong type '{what}'")
}

impl ToolError {
    /// The source's `RequiredParameterNotGiven` thrown with `what`
    /// (`TOPPBase_defs.h:77-87`), as its catch block ends the run
    /// (`TOPPBase.cpp:466-474`): `Error: The required parameter '<what>' was
    /// not given or is empty!` through `writeLogError_`, and
    /// [`ExitCode::MissingParameters`]. `what` is put in single quotes unless
    /// it already starts with one, so a tool that names the parameter
    /// `'in' [valid: mzML]` keeps its own quoting.
    pub fn required_parameter_not_given(what: impl fmt::Display) -> Self {
        Self::caught(
            ExitCode::MissingParameters,
            required_parameter_text(&what.to_string()),
        )
    }

    /// The source's `UnregisteredParameter` thrown with `what`
    /// (`TOPPBase_defs.h:53-63`), as its catch block ends the run
    /// (`TOPPBase.cpp:482-487`): `Internal error: Request for unregistered
    /// parameter '<what>'` and [`ExitCode::InternalError`].
    pub fn unregistered_parameter(what: impl fmt::Display) -> Self {
        Self::caught(
            ExitCode::InternalError,
            unregistered_parameter_text(&what.to_string()),
        )
    }

    /// The source's `WrongParameterType` thrown with `what`
    /// (`TOPPBase_defs.h:65-75`), as its catch block ends the run
    /// (`TOPPBase.cpp:488-493`): `Internal error: Request for parameter with
    /// wrong type '<what>'` and [`ExitCode::InternalError`].
    pub fn wrong_parameter_type(what: impl fmt::Display) -> Self {
        Self::caught(
            ExitCode::InternalError,
            wrong_parameter_type_text(&what.to_string()),
        )
    }
}

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
