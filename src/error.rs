// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use std::fmt;

/// Recoverable errors from validation, parsing and input/output.
#[derive(Debug)]
pub enum Error {
    /// An argument violates an operation's preconditions.
    InvalidValue(String),
    /// An algorithm needs data sorted by position.
    UnsortedData,
    /// Invalid input syntax. Line numbers are one-based when available.
    Parse { line: usize, message: String },
    /// Valid input requests a feature not implemented by this port.
    Unsupported(String),
    /// An underlying stream failed.
    Io(std::io::Error),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidValue(message) => write!(f, "invalid value: {message}"),
            Self::UnsortedData => write!(f, "data must be sorted by position"),
            Self::Parse { line, message } => write!(f, "parse error on line {line}: {message}"),
            Self::Unsupported(message) => write!(f, "unsupported: {message}"),
            Self::Io(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            _ => None,
        }
    }
}

impl From<std::io::Error> for Error {
    fn from(error: std::io::Error) -> Self {
        Self::Io(error)
    }
}

/// Result type shared by the crate.
pub type Result<T> = std::result::Result<T, Error>;
