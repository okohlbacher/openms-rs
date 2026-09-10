// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Text and optional bounded XML adapters. Parsers return owned data only on success.
//!
//! DTA/MGF represent peak lists and a limited subset of experiment metadata;
//! writing them does not preserve auxiliary arrays or arbitrary instrument data.

pub mod dta;
pub mod fasta;
pub mod mgf;

#[cfg(feature = "idxml")]
pub mod idxml;

#[cfg(feature = "mzml")]
pub mod mzml;

use crate::{Error, Result};

fn parse_error(line: usize, message: impl Into<String>) -> Error {
    Error::Parse {
        line,
        message: message.into(),
    }
}

fn number(value: &str, line: usize, label: &str) -> Result<f64> {
    value
        .parse::<f64>()
        .ok()
        .filter(|x| x.is_finite())
        .ok_or_else(|| parse_error(line, format!("invalid finite {label}: {value:?}")))
}

fn intensity(value: &str, line: usize) -> Result<f32> {
    let value = number(value, line, "intensity")? as f32;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(parse_error(line, "intensity exceeds f32 range"))
    }
}

fn single_line(value: &str) -> bool {
    !value.chars().any(char::is_control)
}
