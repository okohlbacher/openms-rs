// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Text and optional bounded XML adapters. Parsers return owned data only on success.
//!
//! DTA/MGF represent peak lists and a limited subset of experiment metadata;
//! writing them does not preserve auxiliary arrays or arbitrary instrument data.

pub mod controlled_vocabulary;
pub mod csv;
#[cfg(feature = "cv-mapping")]
pub mod cv_mapping;
#[cfg(any(feature = "cv-mapping", feature = "mzml-schema"))]
mod cv_xml;
pub mod numpress;
#[cfg(feature = "numpress")]
pub mod numpress_coder;
#[cfg(feature = "semantic-validation")]
pub mod semantic_validator;
#[cfg(feature = "mzml-validation")]
pub use semantic_validator::mzml as mzml_validator;
pub mod text;
pub use csv::CsvFile;
pub use text::TextFile;

#[cfg(feature = "paramxml")]
pub mod paramxml;

pub mod dta;
pub mod dta2d;
pub mod experimental_design_file;
pub mod fasta;
pub mod file_handler;
pub mod file_types;
pub mod mgf;
pub mod ms2;
pub mod peak_options;
pub use peak_options::PeakFileOptions;
pub(crate) mod path_io;

pub use file_handler::FileHandler;
pub use file_types::{FileProperty, FileType, FileTypeList, FilterLayout};

#[cfg(feature = "consensusxml")]
pub mod consensusxml;
#[cfg(feature = "featurexml")]
pub mod featurexml;
#[cfg(any(feature = "idxml", feature = "featurexml", feature = "consensusxml"))]
pub(crate) mod identification_xml;
#[cfg(any(feature = "featurexml", feature = "consensusxml"))]
pub(crate) mod map_xml;
pub mod modification_definitions;

#[cfg(feature = "idxml")]
pub mod idxml;

/// imzML file adapter: load, store and the imaging geometry (`ImzMLFile.h`).
#[cfg(feature = "mzml")]
pub mod imzml_file;
/// Two-file imzML imaging index, geometry and `.ibd` reads (`ImzMLHandler.h`).
#[cfg(feature = "mzml")]
pub mod imzml_handler;
/// Writer for an imzML dataset, `.imzML` plus `.ibd` (`ImzMLWriter.h`).
#[cfg(feature = "mzml")]
pub mod imzml_writer;
#[cfg(feature = "mzml")]
pub mod indexed_mzml;
/// Random access to one record of an indexed mzML file (`IndexedMzMLHandler.h`).
#[cfg(feature = "mzml")]
pub mod indexed_mzml_handler;
#[cfg(feature = "mzml")]
pub mod mzml;
#[cfg(feature = "mzml-schema")]
pub mod mzml_schema;

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
