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
/// OMS files: SQLite databases of identification and feature data (`OMSFile.h`).
#[cfg(feature = "sqlite")]
pub mod oms_file;
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
/// mzQC quality-control files in JSON (`MzQCFile.h`).
#[cfg(feature = "mzml")]
pub mod mzqc_file;
/// MzTab data model: cell vocabulary, record structs and document (`MzTabBase.h`, `MzTab.h`).
pub mod mztab;
/// MzTab exporters: feature maps, consensus maps and identification runs (`MzTab.h` export surface).
pub mod mztab_export;
/// MzTab file adapter: reading and writing `.mzTab` documents (`MzTabFile.h`).
pub mod mztab_file;
/// MzTab-M metabolomics profile: data model, `FeatureMap` export and writer (`MzTabM.h`, `MzTabMFile.h`).
pub mod mztab_m;
/// Legacy mzXML 3.1 adapter: nested scans and paired peak arrays (`MzXMLFile.h`).
#[cfg(feature = "mzml")]
pub mod mzxml;
/// JSON text as the source's `nlohmann::json` writes it, shared by the OMS and
/// mzQC writers.
#[cfg(any(feature = "sqlite", feature = "mzml"))]
pub(crate) mod nlohmann_json;
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
/// MSExperiment spectra and chromatograms as Arrow tables and Parquet files
/// (`MSExperimentArrowExport.h`).
#[cfg(feature = "arrow")]
pub mod msexperiment_arrow_export;
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
/// All-or-nothing writing of a QPX collection (`QPXCollectionExport.h`).
#[cfg(feature = "arrow")]
pub mod qpx_collection_export;
/// Peptide-spectrum matches as Arrow tables and QPX psm Parquet files
/// (`QPXFile.h`).
#[cfg(feature = "arrow")]
pub mod qpx_file;
/// Opaque QPX row identities: BLAKE2b-64 over canonical JSON composites
/// (`QPXIdentity.h`).
#[cfg(feature = "arrow")]
pub mod qpx_identity;
/// Write-time value and primary-key checks of the QPX psm, feature and pg
/// views (`QPXValueValidation.h`).
#[cfg(feature = "arrow")]
pub mod qpx_value_validation;
/// Percolator PSM-level tab-separated output (`PercolatorOutfile.h`).
pub mod percolator_outfile;
/// Spectrum meta-data look-up by regular-expression spectrum reference
/// (`SpectrumLookup.h`, `SpectrumMetaDataLookup.h`).
pub mod spectrum_meta_data_lookup;

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

/// Mascot pepXML peptide sequences keyed by spectrum title (`PepXMLFileMascot.h`).
#[cfg(feature = "idxml")]
pub mod pepxml_mascot;
/// ProteinProphet protXML reader (`ProtXMLFile.h`).
#[cfg(feature = "idxml")]
pub mod protxml;
/// xQuest cross-linking result files (`XQuestResultXMLFile.h`,
/// `HANDLERS/XQuestResultXMLHandler.h`).
#[cfg(feature = "idxml")]
pub mod xquest_result_xml;
/// Inspect input (parameter) files (`InspectInfile.h`).
#[cfg(feature = "idxml")]
pub mod inspect_infile;
/// Inspect result files and trie databases (`InspectOutfile.h`).
#[cfg(feature = "mzml")]
pub mod inspect_outfile;
/// Kroenik (Hardklör sibling) feature files (`KroenikFile.h`).
pub mod kroenik_file;
/// msInspect feature files (`MsInspectFile.h`).
pub mod ms_inspect_file;
/// OMSSA comma-separated result files (`OMSSACSVFile.h`).
pub mod omssa_csv_file;
/// OMSSA XML result files (`OMSSAXMLFile.h`).
#[cfg(feature = "idxml")]
pub mod omssa_xml_file;
/// PepNovo PTM input tables (`PepNovoInfile.h`).
pub mod pepnovo_infile;
/// PepNovo de novo sequencing output (`PepNovoOutfile.h`).
pub mod pepnovo_outfile;
/// PTMXML modification tables (`PTMXMLFile.h`, `HANDLERS/PTMXMLHandler.h`).
#[cfg(feature = "idxml")]
pub mod ptm_xml_file;
/// Permissive sequence parsing shared by the search-engine adapters.
pub mod search_engine_text;
/// Bounded SAX-style driver shared by the X! Tandem and OMSSA XML readers.
#[cfg(feature = "idxml")]
pub mod search_engine_xml;
/// Sequest input (`sequest.params`) files (`SequestInfile.h`).
#[cfg(feature = "idxml")]
pub mod sequest_infile;
/// Sequest `.out` result files (`SequestOutfile.h`).
pub mod sequest_outfile;
/// SpecArray `.pepList` feature files (`SpecArrayFile.h`).
pub mod spec_array_file;
/// X! Tandem input (configuration) files (`XTandemInfile.h`).
pub mod xtandem_infile;
/// X! Tandem XML result files (`XTandemXMLFile.h`).
#[cfg(feature = "idxml")]
pub mod xtandem_xml_file;
/// Random access to a cached mzML pair, metadata plus binary side-car (`CachedMzML.h`).
#[cfg(feature = "mzml")]
pub mod cached_mzml;
/// Binary on-disk cache of spectra and chromatograms (`HANDLERS/CachedMzMLHandler.h`).
pub mod cached_mzml_handler;
/// Chromeleon HPLC text exports into one chromatogram (`ChromeleonFile.h`).
pub mod chromeleon;
/// Enhanced DTA feature tables into consensus maps, and their writers (`EDTAFile.h`).
pub mod edta;
/// FLASHDeconv feature-level output headers (`FLASHDeconvFeatureFile.h`).
pub mod flash_deconv_feature_file;
/// FLASHDeconv spectrum-level output headers (`FLASHDeconvSpectrumFile.h`).
pub mod flash_deconv_spectrum_file;
/// GNPS FBMN meta-value table writer (`GNPSMetaValueFile.h`).
pub mod gnps_meta_value;
/// GNPS FBMN consensus-spectrum MGF writer (`GNPSMGFFile.h`).
#[cfg(all(feature = "consensusxml", feature = "mzml"))]
pub mod gnps_mgf;
/// GNPS FBMN feature-quantification table writer (`GNPSQuantificationFile.h`).
pub mod gnps_quantification;
/// IBSpectra export of isobaric consensus maps for isobar (`IBSpectraFile.h`).
pub mod ibspectra;
/// Bounded streaming line reading with the source's `TextFile::getLine` line ends.
pub mod line_reader;
/// NIST peptide spectral libraries in MSP text: reader and writer (`MSPFile.h`).
pub mod msp;
/// Metabolite MSP spectral libraries: reader and writer (`MSPGenericFile.h`).
pub mod msp_generic;
/// sqMass (SQLite) experiment files over the SQLite handler (`SqMassFile.h`).
#[cfg(feature = "sqmass")]
pub mod sqmass_file;
/// Source `StringUtils::toDouble` / `toInt32` text conversion, shared with the CLI.
pub(crate) mod string_utils;
/// Read-only `acqus` parameters of a Bruker XMass spectrum (`AcqusHandler.h`).
pub mod acqus_handler;
/// Bruker timsTOF `.d` (TDF) reader (`BrukerTimsFile.h`): FRAME export and
/// unaggregated DIA-PASEF export.
#[cfg(feature = "bruker-tdf")]
pub mod bruker_tims_file;
/// Bruker timsTOF MALDI imaging `.d` reader (`BrukerTimsImagingFile.h`).
#[cfg(all(feature = "bruker-tdf", feature = "mzml"))]
pub mod bruker_tims_imaging_file;
/// Read-only `fid` intensities of a Bruker XMass spectrum (`FidHandler.h`).
pub mod fid_handler;
/// IM-axis hill centroiding of one TIMS-PASEF frame (`PASEFHillCentroider.h`).
pub mod pasef_hill_centroider;
/// TIMS scan-index to 1/K0 conversion: the per-frame rational
/// `TimsCalibration` model and the linear fallback (`RationalScan2ImConverter.h`).
pub mod rational_scan2im_converter;
/// Frame-level access to Bruker timsTOF `.d` data: the opentims layer under
/// `BrukerTimsFile.h` (`TimsDataHandle`, frame codec, linear calibrations).
#[cfg(feature = "bruker-tdf")]
pub mod tims_data_handle;
/// Bruker XMass Analysis (`fid`) spectra: load and settings import (`XMassFile.h`).
pub mod xmass_file;
/// PEFF (PSI Extended FASTA Format) annotated protein databases (`PEFFFile.h`).
pub mod peff;
/// PTMXML post-translational modification lists (`PTMXMLFile.h`, `PTMXMLHandler.h`).
#[cfg(feature = "kb-xml")]
pub mod ptm_xml;
/// The Xerces-faithful SAX driver shared by every `XMLHandler` port: the
/// knowledge-base readers (`kb-xml`) and the identification-XML readers
/// (`idxml`).
#[cfg(any(feature = "idxml", feature = "kb-xml"))]
pub(crate) mod sax_xml;
/// Unimod XML modification databases and their registry provider
/// (`UnimodXMLFile.h`, `UnimodXMLHandler.h`, `UnimodXMLDataProvider.h`).
#[cfg(feature = "kb-xml")]
pub mod unimod_xml;
/// UniProtKB XML protein entries for PEFF conversion
/// (`UniProtXMLFile.h`, `UniProtXMLHandler.h`).
#[cfg(feature = "kb-xml")]
pub mod uniprot_xml;
