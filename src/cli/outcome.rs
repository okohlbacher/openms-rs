// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! How a tool body ends: the exit code its `main_` returns, or the exception
//! it throws and `TOPPBase::main` catches (`TOPPBase.cpp:413-514` at cli
//! `c19e494`).
//!
//! The distinction is visible on the console. `TOPPBase::main` prints the
//! closing `<tool> took … .` line on standard output only after `main_` has
//! returned; an exception unwinds past that statement to a catch block, which
//! writes its own diagnostic and picks the exit code. A body that returns
//! [`ExitCode::IllegalParameters`] after writing its own error line and a body
//! whose source throws `IllegalArgument` therefore end differently, although
//! both write an error line: the first prints the closing line, the second
//! does not. [`ToolResult`] carries that distinction, so the lifecycle
//! ([`run_with`](crate::cli::run_with)) prints the closing line only for
//! `Ok`.
//!
//! See `docs/TOPP_CLI_SUPPORT.md`, *Run time*, for the mapping.

use super::parameter::ExitCode;
use crate::Error;
use std::fmt;

/// The result of a tool body ([`Tool::run`](crate::cli::Tool::run)): the
/// exit code `main_` returns, or a [`ToolError`] where the source throws.
pub type ToolResult = std::result::Result<ExitCode, ToolError>;

/// A tool body that ends as the source's `main_` ends when it throws.
///
/// None of the three prints the closing `<tool> took … .` line, because the
/// source's exception unwinds past it. They differ in who writes the
/// diagnostic and where it goes:
///
/// * [`Error`](Self::Error): an error of this crate, which the lifecycle maps
///   to the run-phase catch block of the source exception it stands for
///   (`TOPPBase.cpp:428-499`; [`run_with`](crate::cli::run_with) lists the
///   arms). `?` on a [`crate::Error`] produces this variant.
/// * [`Caught`](Self::Caught): a source exception whose catch-block text the
///   tool knows exactly, such as the `BaseException` arm's `Error: Unexpected
///   internal error (<what>)` for an `IllegalArgument` the source throws. The
///   lifecycle writes the message as the catch block does, through
///   `writeLogError_`: on the error stream, red on a terminal, and in the
///   `-log` file.
/// * [`Escaped`](Self::Escaped): a standard-library exception, which no
///   run-phase catch handles and which reaches the source's last catch
///   (`TOPPBase.cpp:510-513`): `Unable to initialize or run <tool>: <what>`
///   on the error stream only, not in the log file, and
///   [`ExitCode::InternalError`]. Every OpenMS exception is caught in the run
///   phase, so a body has no path to the other outer arm.
#[derive(Debug)]
pub enum ToolError {
    /// An error of this crate, mapped by the lifecycle.
    Error(Error),
    /// A source exception reported by one of the run-phase catch blocks.
    Caught {
        /// The code the catch block returns.
        code: ExitCode,
        /// The catch block's `writeLogError_` text; each line is written as
        /// one log record.
        message: String,
    },
    /// A standard-library exception, reported by the initialisation catch.
    Escaped {
        /// The exception's `what()`.
        what: String,
    },
}

impl ToolError {
    /// A source exception caught by a run-phase catch block that writes
    /// `message` and returns `code`.
    pub fn caught(code: ExitCode, message: impl Into<String>) -> Self {
        Self::Caught {
            code,
            message: message.into(),
        }
    }

    /// An OpenMS exception without a catch block of its own, such as
    /// `IllegalArgument`, `ConversionError` or `FileNotWritable`: the
    /// `BaseException` arm (`TOPPBase.cpp:495-499`), `Error: Unexpected
    /// internal error (<what>)` and [`ExitCode::UnknownError`].
    pub fn unexpected(what: impl fmt::Display) -> Self {
        Self::caught(
            ExitCode::UnknownError,
            format!("Error: Unexpected internal error ({what})"),
        )
    }

    /// The source's `FileNotFound` arm (`TOPPBase.cpp:436-441`): `Error: File
    /// not found (<what>)` and [`ExitCode::InputFileNotFound`].
    pub fn file_not_found(what: impl fmt::Display) -> Self {
        Self::caught(
            ExitCode::InputFileNotFound,
            format!("Error: File not found ({what})"),
        )
    }

    /// A standard-library exception with `what()` = `what`, which reaches the
    /// initialisation catch.
    pub fn escaped(what: impl Into<String>) -> Self {
        Self::Escaped { what: what.into() }
    }

    /// The exit code the lifecycle ends the run with, where the tool decides
    /// it; `None` for [`Error`](Self::Error), which the lifecycle maps.
    pub fn exit_code(&self) -> Option<ExitCode> {
        match self {
            Self::Error(_) => None,
            Self::Caught { code, .. } => Some(*code),
            Self::Escaped { .. } => Some(ExitCode::InternalError),
        }
    }
}

impl From<Error> for ToolError {
    fn from(error: Error) -> Self {
        Self::Error(error)
    }
}

impl From<std::io::Error> for ToolError {
    fn from(error: std::io::Error) -> Self {
        Self::Error(Error::Io(error))
    }
}

impl fmt::Display for ToolError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Error(error) => write!(f, "{error}"),
            Self::Caught { code, message } => write!(f, "{message} ({})", code.name()),
            Self::Escaped { what } => write!(f, "{what} (INTERNAL_ERROR)"),
        }
    }
}

impl std::error::Error for ToolError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Error(error) => Some(error),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_constructors_carry_the_source_catch_texts_and_codes() {
        let error = ToolError::unexpected("Error: Centroided data provided");
        assert_eq!(error.exit_code(), Some(ExitCode::UnknownError));
        assert!(matches!(
            &error,
            ToolError::Caught { message, .. }
                if message == "Error: Unexpected internal error (Error: Centroided data provided)"
        ));
        let error = ToolError::file_not_found("the file 'x' could not be found");
        assert_eq!(error.exit_code(), Some(ExitCode::InputFileNotFound));
        assert!(matches!(
            &error,
            ToolError::Caught { message, .. }
                if message == "Error: File not found (the file 'x' could not be found)"
        ));
        assert_eq!(
            ToolError::escaped("vector::_M_default_append").exit_code(),
            Some(ExitCode::InternalError)
        );
        let error: ToolError = Error::InvalidValue("x".into()).into();
        assert_eq!(error.exit_code(), None);
        let error: ToolError = std::io::Error::other("x").into();
        assert!(matches!(error, ToolError::Error(Error::Io(_))));
    }
}
