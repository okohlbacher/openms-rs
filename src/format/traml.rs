// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! File adapter for HUPO PSI TraML files.
//!
//! Ports `FORMAT/TraMLFile.h` (with its `.cpp`). TraML files contain the
//! transitions used for targeted proteomics and metabolomics experiments:
//!
//! Deutsch et al. "TraML--a standard format for exchange of selected reaction
//! monitoring transition lists." Mol Cell Proteomics. 2012 Apr;11(4):R111.015040.
//! doi: 10.1074/mcp.R111.015040.
//!
//! In OpenMS, TraML files can be generated from TSV or CSV files with
//! `TransitionTSVFile` or the `TargetedFileConverter` TOPP tool; see that
//! class for the TSV format.
//!
//! The document is read into and written from a
//! [`TargetedExperiment`](crate::analysis::openswath::targeted_experiment::TargetedExperiment); the element-by-element work is in
//! [`traml_handler`](crate::format::traml_handler). See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/TRAML_SUPPORT.md)
//! for the member mapping and the evidence.
//!
//! The source class also derives from `ProgressLogger` and reports progress
//! over the transitions while storing; this port reports no progress.

pub use crate::format::traml_handler::{ReadOptions, TraMLLoad};

use crate::Result;
use crate::analysis::openswath::targeted_experiment::TargetedExperiment;
use crate::chemistry::ModificationsDB;
use crate::format::controlled_vocabulary::ControlledVocabulary;
use crate::format::traml_handler::{read_document, write_document};
use std::io::{BufRead, Write};
use std::path::Path;

/// The schema version the source `TraMLFile` declares (`TraML1.0.0.xsd`).
pub const SCHEMA_VERSION: &str = "1.0.0";

/// Load a TraML file, refusing data the source would silently drop.
///
/// Gzip and bzip2 content is decompressed. The PSI-MS vocabulary checks of
/// the source run and their warnings are discarded; use
/// [`load_with_options`] to see them.
///
/// # Errors
///
/// [`Error::Io`](crate::Error::Io) when the file cannot be opened (the
/// source's `Exception::FileNotFound`), and everything
/// [`read_document`] returns (the source's `Exception::ParseError`).
pub fn load(path: impl AsRef<Path>) -> Result<TargetedExperiment> {
    Ok(load_with_options(path, &ReadOptions::default())?.experiment)
}

/// Load a TraML file as the source does: data the source drops is dropped
/// and reported in [`TraMLLoad::warnings`].
///
/// # Errors
///
/// As [`load`], without the refusals of dropped data.
pub fn load_source(path: impl AsRef<Path>) -> Result<TraMLLoad> {
    load_with_options(path, &ReadOptions::source())
}

/// Load a TraML file under explicit options, returning the warnings.
///
/// # Errors
///
/// As [`load`].
pub fn load_with_options(path: impl AsRef<Path>, options: &ReadOptions) -> Result<TraMLLoad> {
    read_with_options(super::path_io::open(path.as_ref())?, options)
}

/// Read a TraML document from a stream under explicit options.
///
/// # Errors
///
/// As [`read_document`].
pub fn read_with_options(input: impl BufRead, options: &ReadOptions) -> Result<TraMLLoad> {
    read_document(input, options, ControlledVocabulary::psi_ms().ok())
}

/// Store an experiment as a TraML file.
///
/// A path ending in `.gz` or `.bz2` is written compressed, as the source's
/// `XMLFile::save_` does. The file is written to a sibling and published
/// only when the whole document was written, so an error leaves no partial
/// file.
///
/// # Errors
///
/// [`Error::Io`](crate::Error::Io) when the file cannot be created (the
/// source's `Exception::UnableToCreateFile`), and everything
/// [`write_document`] returns.
pub fn store(path: impl AsRef<Path>, exp: &TargetedExperiment) -> Result<()> {
    super::path_io::write(path.as_ref(), |out| {
        write_document(out, exp, ModificationsDB::global())
    })
}

/// Write an experiment as a TraML document to a stream.
///
/// # Errors
///
/// As [`write_document`].
pub fn write(mut out: impl Write, exp: &TargetedExperiment) -> Result<()> {
    write_document(&mut out, exp, ModificationsDB::global())
}

/// Check a file against the TraML CV mapping rules and the controlled
/// vocabularies; the source `isSemanticallyValid`.
///
/// The source validates with `TraMLValidator`, which is the generic
/// `SemanticValidator` with unit checking switched on, against the bundled
/// `TraML-mapping.xml`, the PSI-MS vocabulary and the unit ontology. This
/// uses the same mapping and rules with
/// [`ControlledVocabulary::psi_ms`], which also carries the other bundled
/// vocabularies (a mapping rule naming only `MS` and `UO` terms is
/// unaffected). The report's errors and warnings are the source's
/// `errors` and `warnings` output lists, and its
/// [`is_valid`](crate::format::semantic_validator::ValidationReport::is_valid)
/// is the source's return value.
///
/// # Errors
///
/// [`Error::Io`](crate::Error::Io) when the file cannot be opened (the
/// source's `Exception::FileNotFound`), and a resource or XML error from
/// the validator.
#[cfg(feature = "semantic-validation")]
pub fn is_semantically_valid(
    path: impl AsRef<Path>,
) -> Result<crate::format::semantic_validator::ValidationReport> {
    use crate::data_structures::cv_mapping::CVMappings;
    use crate::format::cv_mapping::CVMappingFile;
    use std::sync::OnceLock;
    static MAPPING: OnceLock<std::result::Result<CVMappings, String>> = OnceLock::new();
    let mapping = match MAPPING.get_or_init(|| {
        CVMappingFile::default()
            .read(include_bytes!("../../resources/cv/TraML-mapping.xml").as_slice())
            .map_err(|e| e.to_string())
    }) {
        Ok(mapping) => mapping,
        Err(message) => return Err(crate::Error::InvalidValue(message.clone())),
    };
    let mut validator = crate::format::semantic_validator::SemanticValidator::new(
        mapping,
        ControlledVocabulary::psi_ms()?,
    );
    validator.options.check_units = true;
    validator.validate(path)
}
