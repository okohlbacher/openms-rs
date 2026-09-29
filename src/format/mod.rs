// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Text and optional bounded XML adapters. Parsers return owned data only on success.
//!
//! DTA/MGF represent peak lists and a limited subset of experiment metadata;
//! writing them does not preserve auxiliary arrays or arbitrary instrument data.

/// Controlled-vocabulary records, OBO parsing and checked ontology queries.
pub mod controlled_vocabulary;
/// Delimited text tables with source-compatible quoting and field access.
pub mod csv;
/// Controlled-vocabulary mapping records and mapping XML input.
#[cfg(feature = "cv-mapping")]
pub mod cv_mapping;
#[cfg(any(feature = "cv-mapping", feature = "xml-schema"))]
mod cv_xml;
/// Raw MS-Numpress numeric codecs.
pub mod numpress;
/// MS-Numpress configuration, Base64 transport and optional zlib compression.
#[cfg(feature = "numpress")]
pub mod numpress_coder;
/// Controlled-vocabulary validation against mapping rules.
#[cfg(feature = "semantic-validation")]
pub mod semantic_validator;
#[cfg(feature = "mzml-validation")]
pub use semantic_validator::mzml as mzml_validator;
/// Line-oriented text loading and writing.
pub mod text;
pub use csv::CsvFile;
pub use text::TextFile;

/// OpenMS parameter trees in INI XML format.
#[cfg(feature = "paramxml")]
pub mod paramxml;

/// Single-spectrum DTA peak-list input and output.
pub mod dta;
/// Retention-time-indexed DTA2D peak-list input and output.
pub mod dta2d;
/// Tabular experimental-design loading and validation.
pub mod experimental_design_file;
/// Streaming FASTA sequence input and output.
pub mod fasta;
/// File-type-based dispatch to implemented native format adapters.
pub mod file_handler;
/// File summaries for the FileInfo tool (`FileInfo.h`).
pub mod file_info;
/// File-type identities, properties, extensions and filter labels.
pub mod file_types;
/// Mascot generic format (MGF) peak lists and search header (`MascotGenericFile.h`).
pub mod mascot_generic;
/// Streaming MGF peak lists and their supported text metadata.
pub mod mgf;
/// MS2 peak-list input and output with explicit transport boundaries.
pub mod ms2;
/// Long-format MSstats and MSstatsTMT CSV writer (`MSstatsFile.h`).
pub mod msstats;
/// Scientific peak-file filters, precision, compression and loading options.
pub mod peak_options;
/// Spectrum peak-type estimation (`PeakTypeEstimator.h`).
pub mod peak_type_estimator;
/// Percolator tab-separated input: writing, reading and the PIN feature set
/// (`PercolatorInfile.h`).
pub mod percolator_infile;
/// Retention-time transformation persistence as TrafoXML (`TransformationXMLFile.h`).
#[cfg(any(feature = "featurexml", feature = "consensusxml"))]
pub mod transformation_xml;
pub use peak_options::PeakFileOptions;
/// Bounded sqMass experiment storage, including compressed mzML metadata.
#[cfg(feature = "sqmass")]
pub mod mzml_sqlite_handler;
/// Read-only SWATH/DIA windows and spectrum IDs from sqMass databases.
#[cfg(feature = "sqlite")]
pub mod mzml_sqlite_swath_handler;
pub(crate) mod path_io;
/// SQLite connections and checked table, statement and blob operations.
#[cfg(feature = "sqlite")]
pub mod sqlite_connector;
/// Separated-value output with automatic separators and quoting (`SVOutStream.h`).
pub mod sv_out_stream;

pub use file_handler::FileHandler;
pub use file_types::{FileProperty, FileType, FileTypeList, FilterLayout};

/// Consensus feature maps and identification metadata in consensusXML.
#[cfg(feature = "consensusxml")]
pub mod consensusxml;
/// Feature maps and identification metadata in featureXML.
#[cfg(feature = "featurexml")]
pub mod featurexml;
#[cfg(any(feature = "idxml", feature = "featurexml", feature = "consensusxml"))]
pub(crate) mod identification_xml;
#[cfg(any(feature = "featurexml", feature = "consensusxml"))]
pub(crate) mod map_xml;
/// Modification definitions from identification and feature documents.
pub mod modification_definitions;

/// Protein and peptide identifications in idXML.
#[cfg(feature = "idxml")]
pub mod idxml;
/// Mascot XML search-result reader and its title lookup (`MascotXMLFile.h`).
#[cfg(feature = "idxml")]
pub mod mascot_xml;

/// imzML file adapter: load, store and the imaging geometry (`ImzMLFile.h`).
#[cfg(feature = "mzml")]
pub mod imzml_file;
/// Two-file imzML imaging index, geometry and `.ibd` reads (`ImzMLHandler.h`).
#[cfg(feature = "mzml")]
pub mod imzml_handler;
/// Writer for an imzML dataset, `.imzML` plus `.ibd` (`ImzMLWriter.h`).
#[cfg(feature = "mzml")]
pub mod imzml_writer;
/// Indexed mzML footer discovery and checked offset decoding.
#[cfg(feature = "mzml")]
pub mod indexed_mzml;
/// Random access to one record of an indexed mzML file (`IndexedMzMLHandler.h`).
#[cfg(feature = "mzml")]
pub mod indexed_mzml_handler;
/// Streaming mzML consumer that writes records as they arrive (`MSDataWritingConsumer.h`).
#[cfg(feature = "mzml")]
pub mod ms_data_writing_consumer;
/// mzData 1.05 file adapter and handler (`MzDataFile.h`, `MzDataHandler.h`).
#[cfg(feature = "mzml")]
pub mod mzdata;
/// mzIdentML adapter: PSI identification interchange (`MzIdentMLFile.h`).
#[cfg(feature = "idxml")]
pub mod mzidentml;
/// mzML spectra, chromatograms, headers, scientific loading and writing.
#[cfg(feature = "mzml")]
pub mod mzml;
/// Explicit validation against the retained ordinary and indexed mzML schemas.
#[cfg(feature = "mzml-schema")]
pub mod mzml_schema;
/// MzTab data model: cell vocabulary, record structs and document (`MzTabBase.h`, `MzTab.h`).
pub mod mztab;
/// MzTab file adapter: reading and writing `.mzTab` documents (`MzTabFile.h`).
pub mod mztab_file;
/// MzTab-M metabolomics profile: data model, `FeatureMap` export and writer (`MzTabM.h`, `MzTabMFile.h`).
pub mod mztab_m;
/// Legacy mzXML 3.1 adapter: nested scans and paired peak arrays (`MzXMLFile.h`).
#[cfg(feature = "mzml")]
pub mod mzxml;
/// pepXML search results: load, store and modification resolution (`PepXMLFile.h`).
#[cfg(feature = "idxml")]
pub mod pepxml;
/// XSD validation against every bundled source schema, or a caller's own
/// self-contained schema (`XMLValidator.h`, `XMLFile::isValid`).
#[cfg(feature = "xml-schema")]
pub mod xml_schema;

/// qcML quality-control reports: runs, sets, quality parameters, attachments
/// and their XML and table serialisations (`QcMLFile.h`).
#[cfg(feature = "paramxml")]
pub mod qcml;

/// Writing, concatenating and reading Arrow tables as Parquet, and the QPX
/// naming helpers the Parquet exporters share (`ArrowIOHelpers.h`).
#[cfg(feature = "arrow")]
pub mod arrow_io_helpers;
/// Arrow schema definitions of the OpenMS Parquet tables and their
/// validation (`ArrowSchemaRegistry.h`).
#[cfg(feature = "arrow")]
pub mod arrow_schema_registry;
/// The features table of a FeatureMap as an Arrow table
/// (`FeatureMapArrowIO.h`, features half).
#[cfg(feature = "arrow")]
pub mod feature_map_arrow_io;
/// OpenSWATH mobilograms to `.xim` Parquet
/// (`DATAACCESS/MobilogramParquetConsumer.h`).
#[cfg(feature = "arrow")]
pub mod mobilogram_parquet_consumer;
/// OpenSWATH chromatograms to `.xic` Parquet
/// (`DATAACCESS/MSChromatogramParquetConsumer.h`).
#[cfg(feature = "arrow")]
pub mod ms_chromatogram_parquet_consumer;
/// Shared writer machinery of the OpenSWATH Parquet consumers.
#[cfg(feature = "arrow")]
mod openswath_parquet_writer;
/// Shared Parquet read/write helpers and type-coercing column accessors
/// (`ParquetFile.h`).
#[cfg(feature = "arrow")]
pub mod parquet_file;
/// Typed filter expressions for Parquet-backed readers (`ParquetFilter.h`).
pub mod parquet_filter;
/// Shared Arrow/Parquet reader helpers of the OpenSWATH readers
/// (`DATAACCESS/ParquetReaderHelpers.h`, private in the source too).
#[cfg(feature = "arrow")]
mod parquet_reader_helpers;
/// Protein hits, protein groups and search parameters as Arrow tables and
/// Parquet files (`ProteinIdentificationArrowIO.h`).
#[cfg(feature = "arrow")]
pub mod protein_identification_arrow_io;
/// OpenSWATH chromatogram Parquet files, `.xic` (`XICParquetFile.h`).
#[cfg(feature = "arrow")]
pub mod xic_parquet_file;
/// OpenSWATH mobilogram Parquet files, `.xim` (`XIMParquetFile.h`).
#[cfg(feature = "arrow")]
pub mod xim_parquet_file;
/// Targeted peak maps to `.xipm` Parquet
/// (`DATAACCESS/XIPMParquetConsumer.h`).
#[cfg(feature = "arrow")]
pub mod xipm_parquet_consumer;
/// Extracted ion peak-map Parquet files, `.xipm` (`XIPMParquetFile.h`).
#[cfg(feature = "arrow")]
pub mod xipm_parquet_file;

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
