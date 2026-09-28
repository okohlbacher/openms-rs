// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Hendrik Weisser, OpenMS Rust contributors $

//! OMS files: SQLite databases of identification and feature data
//! (`FORMAT/OMSFile.h`, with the internal helpers `FORMAT/OMSFileLoad.h` and
//! `FORMAT/OMSFileStore.h`).
//!
//! An `.oms` file is an SQLite database whose tables mirror OpenMS's
//! `IdentificationData` (`ID_*` tables), feature and consensus maps (`FEAT_*`),
//! plus shared tables for controlled-vocabulary terms (`CVTerm`), metadata
//! values (`DataValue_DataType` and one `<table>_MetaInfo` table per record
//! type) and a `version` table whose `OMSFile` column is the schema version.
//!
//! [`OmsFile`] ports all seven public operations of the source class:
//!
//! | Source | Port |
//! |---|---|
//! | `store(filename, const IdentificationData&)` | [`OmsFile::store_identifications`] |
//! | `store(filename, const FeatureMap&)` | [`OmsFile::store_feature_map`] |
//! | `store(filename, const ConsensusMap&)` | [`OmsFile::store_consensus_map`] |
//! | `load(filename, IdentificationData&)` | [`OmsFile::load_identifications`] |
//! | `load(filename, FeatureMap&)` | [`OmsFile::load_feature_map`] |
//! | `load(filename, ConsensusMap&)` | [`OmsFile::load_consensus_map`] |
//! | `exportToJSON(filename_in, filename_out)` | [`OmsFile::export_to_json`] |
//!
//! The source's `FeatureMap` and `ConsensusMap` own an `IdentificationData`;
//! the crate's maps reference an identification graph held alongside them
//! (see [`crate::kernel::feature_identification`]), so the map stores take the
//! graph as a separate argument and the map loads return it.
//!
//! Where the source loses or corrupts information, the port does not by
//! default and an option selects the source's behaviour
//! ([`OmsStoreOptions::source`], [`OmsLoadOptions::source`]); see
//! `docs/OMS_FILE_SUPPORT.md`.
//!
//! The module needs the optional `sqlite` feature, which brings in the bundled
//! SQLite of `rusqlite`.
//!
//! [`OmsFile`]: crate::format::oms_file::OmsFile
//! [`OmsFile::store_identifications`]: crate::format::oms_file::OmsFile::store_identifications
//! [`OmsFile::store_feature_map`]: crate::format::oms_file::OmsFile::store_feature_map
//! [`OmsFile::store_consensus_map`]: crate::format::oms_file::OmsFile::store_consensus_map
//! [`OmsFile::load_identifications`]: crate::format::oms_file::OmsFile::load_identifications
//! [`OmsFile::load_feature_map`]: crate::format::oms_file::OmsFile::load_feature_map
//! [`OmsFile::load_consensus_map`]: crate::format::oms_file::OmsFile::load_consensus_map
//! [`OmsFile::export_to_json`]: crate::format::oms_file::OmsFile::export_to_json
//! [`OmsStoreOptions::source`]: crate::format::oms_file::OmsStoreOptions::source
//! [`OmsLoadOptions::source`]: crate::format::oms_file::OmsLoadOptions::source

mod load;
mod store;

use crate::concept::progress_logger::{ProgressLogType, ProgressLogger};
use crate::format::path_io;
use crate::identification::graph::IdentificationData;
use crate::kernel::{ConsensusMap, FeatureMap};
use crate::{Error, Result};
use rusqlite::types::ValueRef;
use rusqlite::{Connection, OpenFlags};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

/// The schema version the pinned source writes (`OMSFileStore.cpp`
/// `version_number`); readers accept any version the `version` table states.
pub const SCHEMA_VERSION: i64 = 5;

/// Largest JSON text [`OmsFile::export_to_json`] produces, in bytes.
///
/// The source builds the whole JSON document in memory without a ceiling
/// (its own `@TODO` notes it); this port refuses beyond the limit before the
/// output file is created.
pub const MAX_JSON_BYTES: usize = 1 << 30;

/// Largest number of tables an exported database may hold.
pub const MAX_TABLES: usize = 10_000;

/// Largest number of rows one query of a load may return (the source has no
/// ceiling). The identification graph adds its own limits
/// ([`crate::identification::graph::GraphLimits`]).
pub const MAX_TABLE_ROWS: usize = 10_000_000;

/// Deepest nesting of subordinate features a load or store follows (the
/// source recurses without a limit).
pub const MAX_SUBORDINATE_DEPTH: usize = 256;

/// Choices for the stores where the source writes something other than the
/// data it was given. The default writes the data; [`OmsStoreOptions::source`]
/// reproduces the pinned Release build.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct OmsStoreOptions {
    /// Write a compound's name into the `smile` column, as the source's
    /// `storeIdentifiedCompounds_` does (it binds `compound.name` twice), and
    /// so lose the SMILES string.
    pub compound_smile_from_name: bool,
    /// Write a consensus feature's ratios under the key of its last feature
    /// handle, as the source does (its ratio loop runs after the handle loop
    /// has advanced the key), instead of under the consensus feature's own key.
    pub consensus_ratios_under_last_handle: bool,
    /// Write every feature handle's unique ID as zero: the source converts a
    /// handle with `BaseFeature(const FeatureHandle&)`, which does not copy it.
    pub drop_handle_unique_ids: bool,
    /// Write a compound's formula without its charge: the source writes
    /// `EmpiricalFormula::toString()`, which omits it.
    pub drop_formula_charge: bool,
}

impl OmsStoreOptions {
    /// The pinned source's behaviour in every respect.
    pub fn source() -> Self {
        Self {
            compound_smile_from_name: true,
            consensus_ratios_under_last_handle: true,
            drop_handle_unique_ids: true,
            drop_formula_charge: true,
        }
    }
}

/// Choices for the loads where the source reads something other than what the
/// file holds. The default reads the file; [`OmsLoadOptions::source`]
/// reproduces the pinned Release build.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct OmsLoadOptions {
    /// Keep the space `DataValue::toString()` writes after each comma of a
    /// string list as part of the next item, as the source's
    /// `ListUtils::create<String>` does; each store/load round trip then adds
    /// one more space. The default removes that one separator space.
    pub keep_string_list_spaces: bool,
    /// Leave a map's loaded file type unknown, as the source does (it has no
    /// setter for the stored type name); the default restores it.
    pub unknown_file_type: bool,
    /// Load no consensus ratios, as the source does (it reads them into a copy
    /// of the feature it has already added to the map).
    pub drop_consensus_ratios: bool,
}

impl OmsLoadOptions {
    /// The pinned source's behaviour in every respect.
    pub fn source() -> Self {
        Self {
            keep_string_list_spaces: true,
            unknown_file_type: true,
            drop_consensus_ratios: true,
        }
    }
}

/// A loaded feature or consensus map with the identification graph its
/// features reference.
#[derive(Debug)]
pub struct OmsMapLoad<M> {
    /// The map.
    pub map: M,
    /// The identification data of the file; empty when the file has none.
    pub identifications: IdentificationData,
    /// Data processing action names the file holds but the port does not know.
    /// The source logs "unknown data processing action" and skips each; so
    /// does the port, and returns them here.
    pub skipped_actions: Vec<String>,
}

pub(crate) fn sql_error(context: &str, error: rusqlite::Error) -> Error {
    // Source: raiseDBError_ throws Exception::FailedAPICall with
    // "<context>: <SQLite message>".
    Error::Io(std::io::Error::other(format!("{context}: {error}")))
}

fn limit() -> Error {
    Error::InvalidValue("OMS JSON export exceeds its size limit".into())
}

/// Row order of the tables that do not order by `id` (the source's
/// `OMSFileLoad::export_order_by_`); `version` is exported unordered.
fn order_by(table: &str) -> Option<&'static str> {
    if table.ends_with("_MetaInfo") {
        return Some("parent_id, name");
    }
    if table.ends_with("_AppliedProcessingStep") {
        return Some("parent_id, processing_step_order, score_type_id");
    }
    match table {
        "version" => None,
        "ID_IdentifiedCompound" => Some("molecule_id"),
        "ID_ParentMatch" => Some("molecule_id, parent_id, start_pos, end_pos"),
        "ID_ParentGroup_ParentSequence" => Some("group_id, parent_id"),
        "ID_ProcessingStep_InputFile" => Some("processing_step_id, input_file_id"),
        "ID_ProcessingSoftware_AssignedScore" => Some("software_id, score_type_order"),
        "ID_ObservationMatch_PeakAnnotation" => {
            Some("parent_id, processing_step_id, peak_mz, peak_annotation")
        }
        "FEAT_ConvexHull" => Some("feature_id, hull_index, point_index"),
        "FEAT_ObservationMatch" => Some("feature_id, observation_match_id"),
        "FEAT_MapMetaData" => Some("unique_id"),
        _ => Some("id"),
    }
}

/// An SQL identifier in double quotes, with embedded quotes doubled.
fn quoted(identifier: &str) -> String {
    format!("\"{}\"", identifier.replace('"', "\"\""))
}

/// Reader and writer of OMS files. Source `OMSFile`.
///
/// The source derives from `ProgressLogger` and takes its log type in the
/// constructor; the log type is kept here for the load and store operations.
/// The JSON export reports no progress, in the source as here.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct OmsFile {
    log_type: ProgressLogType,
}

impl OmsFile {
    /// An adapter with the given progress log type. Source
    /// `OMSFile(LogType log_type = LogType::NONE)`.
    pub fn new(log_type: ProgressLogType) -> Self {
        Self { log_type }
    }

    /// The progress log type (the source's protected `log_type_`, also read
    /// through the inherited `ProgressLogger::getLogType`).
    pub fn log_type(&self) -> ProgressLogType {
        self.log_type
    }

    /// Change the progress log type (the inherited
    /// `ProgressLogger::setLogType`). The source's load and store helpers are
    /// built with the `log_type_` given to the constructor, which
    /// `setLogType` does not update; here the loads and stores use the type
    /// set last.
    pub fn set_log_type(&mut self, log_type: ProgressLogType) {
        self.log_type = log_type;
    }

    /// Open `filename` read-only and read its schema version, as the source's
    /// `OMSFileLoad` constructor does.
    fn open(filename: &Path) -> Result<(Connection, i64)> {
        let connection = Connection::open_with_flags(
            filename,
            OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
        )
        .map_err(|error| sql_error("error opening OMS file", error))?;
        let version = connection
            .query_row("SELECT OMSFile FROM version", [], |row| {
                row.get::<_, i64>(0)
            })
            .map_err(|error| sql_error("error reading file format version number", error))?;
        Ok((connection, version))
    }

    /// The schema version an OMS file declares in its `version` table.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be opened as an SQLite database or
    /// has no readable `version` table (the source's
    /// `Exception::FailedAPICall`, "error reading file format version number").
    pub fn schema_version(&self, filename: impl AsRef<Path>) -> Result<i64> {
        Ok(Self::open(filename.as_ref())?.1)
    }

    /// Render the contents of an OMS file as JSON text.
    ///
    /// Source `OMSFileLoad::exportToJSON`: one object member per table (every
    /// table except SQLite's own), whose value is an array with one object per
    /// row; a row object maps column names to values. Integers and reals keep
    /// their SQLite storage class, text and blobs become strings, and `NULL`
    /// becomes the empty string. Rows are ordered by `id`, by
    /// `parent_id, name` for `*_MetaInfo` tables, by
    /// `parent_id, processing_step_order, score_type_id` for
    /// `*_AppliedProcessingStep` tables and by the source's fixed column lists
    /// for eight tables without an `id` column; `version` is unordered.
    ///
    /// The text is formatted as `nlohmann::json::dump(4)` formats it: members
    /// sorted by name, four-space indentation, `": "` between name and value,
    /// shortest round-trip numbers with a two-digit exponent and a `.0` on
    /// integral reals, and a final newline.
    ///
    /// # Errors
    ///
    /// As [`OmsFile::schema_version`]; [`Error::Io`] when a table cannot be
    /// read; [`Error::InvalidValue`] when a text or blob cell is not UTF-8 (the
    /// source's JSON library throws `type_error.316` there) or the text would
    /// exceed [`MAX_JSON_BYTES`], or the database holds more than
    /// [`MAX_TABLES`] tables.
    pub fn json_text(&self, filename: impl AsRef<Path>) -> Result<String> {
        let (connection, _) = Self::open(filename.as_ref())?;
        let mut tables: Vec<String> = Vec::new();
        {
            let mut query = connection
                .prepare(
                    "SELECT name FROM sqlite_master WHERE type='table' \
                     AND name NOT LIKE 'sqlite_%' ORDER BY name",
                )
                .map_err(|error| sql_error("error listing tables", error))?;
            let mut rows = query
                .query([])
                .map_err(|error| sql_error("error listing tables", error))?;
            while let Some(row) = rows
                .next()
                .map_err(|error| sql_error("error listing tables", error))?
            {
                if tables.len() >= MAX_TABLES {
                    return Err(Error::InvalidValue(
                        "OMS file holds more tables than the export accepts".into(),
                    ));
                }
                tables.push(
                    row.get::<_, String>(0)
                        .map_err(|error| sql_error("error listing tables", error))?,
                );
            }
        }
        // The source's JSON object is a std::map, which orders members
        // bytewise; SQLite's ORDER BY name uses the same BINARY collation.
        tables.sort();
        let mut out = JsonWriter::default();
        out.raw("{")?;
        for (table_index, table) in tables.iter().enumerate() {
            out.raw(if table_index == 0 { "\n" } else { ",\n" })?;
            out.indent(1)?;
            out.string(table)?;
            out.raw(": ")?;
            export_table(&connection, table, &mut out)?;
        }
        if !tables.is_empty() {
            out.raw("\n")?;
        }
        out.raw("}\n")?;
        Ok(out.text)
    }

    /// Write the contents of the OMS file `filename_in` as JSON to
    /// `filename_out`. Source `OMSFile::exportToJSON`.
    ///
    /// The text is [`OmsFile::json_text`]; it is rendered completely before the
    /// output is created and written atomically.
    ///
    /// # Errors
    ///
    /// As [`OmsFile::json_text`], plus [`Error::Io`] when the output cannot be
    /// written (the source's `Exception::FileNotWritable`).
    pub fn export_to_json(
        &self,
        filename_in: impl AsRef<Path>,
        filename_out: impl AsRef<Path>,
    ) -> Result<()> {
        let text = self.json_text(filename_in)?;
        path_io::write(filename_out.as_ref(), |writer| {
            writer.write_all(text.as_bytes())?;
            Ok(())
        })
    }

    fn logger(&self) -> ProgressLogger {
        let mut logger = ProgressLogger::new();
        logger.set_log_type(self.log_type);
        logger
    }

    /// Write identification data to an OMS file. Source
    /// `store(const std::string&, const IdentificationData&)`.
    ///
    /// Every record type gets its table (created only when it has records),
    /// with database keys numbered from 1 in the graph's record order; scores
    /// and processing steps go to `<table>_AppliedProcessingStep` and metadata
    /// to `<table>_MetaInfo`, the value as `DataValue::toString()` text with
    /// its type. The processing actions of a processing step, the metadata of
    /// score types, and the units of meta values are not stored, in the source
    /// as here.
    ///
    /// The source deletes `filename` and writes the database in its place;
    /// the port writes a temporary file beside it and renames it over
    /// `filename` only when the whole store has succeeded, so a failure leaves
    /// any existing file unchanged.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the database cannot be created or written (the
    /// source's `Exception::FailedAPICall` and SQLite exceptions), including a
    /// constraint violation the source also hits (two parent sequences with one
    /// accession cannot occur, but a parent group with several scores writes
    /// its key twice); [`Error::InvalidValue`] when a count or position does
    /// not fit the 32 bits the source binds (the source truncates).
    pub fn store_identifications(
        &self,
        filename: impl AsRef<Path>,
        data: &IdentificationData,
    ) -> Result<()> {
        self.store_identifications_with(filename, data, OmsStoreOptions::default())
    }

    /// [`OmsFile::store_identifications`] with explicit options.
    ///
    /// # Errors
    ///
    /// As [`OmsFile::store_identifications`].
    pub fn store_identifications_with(
        &self,
        filename: impl AsRef<Path>,
        data: &IdentificationData,
        options: OmsStoreOptions,
    ) -> Result<()> {
        let mut logger = self.logger();
        write_database(filename.as_ref(), |db| {
            logger.start_progress(0, 13, "Writing identification data to file")?;
            let mut writer = store::Store::new(db, options);
            writer.store_identifications(data, &mut || logger.next_progress())?;
            logger.next_progress()?;
            logger.end_progress(0)
        })
    }

    /// Write a feature map and the identification data its features reference.
    /// Source `store(const std::string&, const FeatureMap&)`, whose map owns the
    /// identification data; `identifications` may be `None` or empty, and
    /// then only the `version` table precedes the feature tables, as in the
    /// source.
    ///
    /// Features are numbered from 0, depth-first, each subordinate after its
    /// parent and before the parent's next sibling.
    ///
    /// # Errors
    ///
    /// As [`OmsFile::store_identifications`]; [`Error::InvalidValue`] when a
    /// feature references an observation match or molecule that is not in
    /// `identifications`, or subordinates nest deeper than
    /// [`MAX_SUBORDINATE_DEPTH`].
    pub fn store_feature_map(
        &self,
        filename: impl AsRef<Path>,
        features: &FeatureMap,
        identifications: Option<&IdentificationData>,
    ) -> Result<()> {
        self.store_feature_map_with(
            filename,
            features,
            identifications,
            OmsStoreOptions::default(),
        )
    }

    /// [`OmsFile::store_feature_map`] with explicit options.
    ///
    /// # Errors
    ///
    /// As [`OmsFile::store_feature_map`].
    pub fn store_feature_map_with(
        &self,
        filename: impl AsRef<Path>,
        features: &FeatureMap,
        identifications: Option<&IdentificationData>,
        options: OmsStoreOptions,
    ) -> Result<()> {
        let mut logger = self.logger();
        write_database(filename.as_ref(), |db| {
            let mut writer = store::Store::new(db, options);
            writer.store_feature_map(features, identifications)?;
            logger.start_progress(
                0,
                crate::concept::progress_logger::progress_value(features.features.len())?,
                "Writing feature data to file",
            )?;
            logger.end_progress(0)
        })
    }

    /// Write a consensus map and the identification data its features
    /// reference. Source `store(const std::string&, const ConsensusMap&)`.
    ///
    /// Each consensus feature is written, then each of its handles as a
    /// subordinate row with the next key, carrying the handle's measured
    /// values and map index (`FEAT_FeatureHandle`); ratios go to
    /// `FEAT_ConsensusRatio`.
    ///
    /// # Errors
    ///
    /// As [`OmsFile::store_feature_map`].
    pub fn store_consensus_map(
        &self,
        filename: impl AsRef<Path>,
        consensus: &ConsensusMap,
        identifications: Option<&IdentificationData>,
    ) -> Result<()> {
        self.store_consensus_map_with(
            filename,
            consensus,
            identifications,
            OmsStoreOptions::default(),
        )
    }

    /// [`OmsFile::store_consensus_map`] with explicit options.
    ///
    /// # Errors
    ///
    /// As [`OmsFile::store_feature_map`].
    pub fn store_consensus_map_with(
        &self,
        filename: impl AsRef<Path>,
        consensus: &ConsensusMap,
        identifications: Option<&IdentificationData>,
        options: OmsStoreOptions,
    ) -> Result<()> {
        let mut logger = self.logger();
        write_database(filename.as_ref(), |db| {
            let mut writer = store::Store::new(db, options);
            writer.store_consensus_map(consensus, identifications)?;
            logger.start_progress(
                0,
                crate::concept::progress_logger::progress_value(consensus.features.len())?,
                "Writing consensus feature data to file",
            )?;
            logger.end_progress(0)
        })
    }

    /// Read the identification data of an OMS file. Source
    /// `load(const std::string&, IdentificationData&)`, which adds to the
    /// graph it is given; the port returns a new graph, so a failure leaves
    /// nothing half-loaded.
    ///
    /// Every schema version the source reads is accepted (1 to 5). Records
    /// are registered as the source registers them, so equal records merge
    /// and the graph's own ordering applies.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file is not an SQLite database with a `version`
    /// table or a query fails; [`Error::InvalidValue`] when a row references a
    /// key that is not in the file (the source then uses an invalid reference),
    /// a cell cannot be converted (the source's `Exception::ConversionError`),
    /// `ID_ScoreType` exists without `CVTerm` (the source's
    /// `Exception::MissingInformation`), a query returns more than
    /// [`MAX_TABLE_ROWS`] rows, or the graph refuses a record or exceeds its
    /// limits; the enzyme lookup's error for an unknown enzyme name (the
    /// source's `Exception::ElementNotFound`).
    pub fn load_identifications(&self, filename: impl AsRef<Path>) -> Result<IdentificationData> {
        self.load_identifications_with(filename, OmsLoadOptions::default())
    }

    /// [`OmsFile::load_identifications`] with explicit options.
    ///
    /// # Errors
    ///
    /// As [`OmsFile::load_identifications`].
    pub fn load_identifications_with(
        &self,
        filename: impl AsRef<Path>,
        options: OmsLoadOptions,
    ) -> Result<IdentificationData> {
        let (db, version) = Self::open(filename.as_ref())?;
        let mut logger = self.logger();
        let mut data = IdentificationData::new()?;
        logger.start_progress(0, 12, "Reading identification data from file")?;
        let mut reader = load::Load::new(&db, version, options);
        reader.load_identifications(&mut data, &mut || logger.next_progress())?;
        logger.end_progress(0)?;
        Ok(data)
    }

    /// Read a feature map and its identification data. Source
    /// `load(const std::string&, FeatureMap&)`.
    ///
    /// Each feature's width is set through [`crate::kernel::BaseFeature::set_width`],
    /// which also records it as the `FWHM` meta value, as the source's
    /// `setWidth` does; stored metadata is read afterwards and overrides it.
    ///
    /// # Errors
    ///
    /// As [`OmsFile::load_identifications`]; [`Error::InvalidValue`] when
    /// subordinates nest deeper than [`MAX_SUBORDINATE_DEPTH`], a width is
    /// negative or not finite, or the map metadata table is empty.
    pub fn load_feature_map(&self, filename: impl AsRef<Path>) -> Result<OmsMapLoad<FeatureMap>> {
        self.load_feature_map_with(filename, OmsLoadOptions::default())
    }

    /// [`OmsFile::load_feature_map`] with explicit options.
    ///
    /// # Errors
    ///
    /// As [`OmsFile::load_feature_map`].
    pub fn load_feature_map_with(
        &self,
        filename: impl AsRef<Path>,
        options: OmsLoadOptions,
    ) -> Result<OmsMapLoad<FeatureMap>> {
        let (db, version) = Self::open(filename.as_ref())?;
        let mut data = IdentificationData::new()?;
        let mut reader = load::Load::new(&db, version, options);
        reader.load_identifications(&mut data, &mut || Ok(()))?;
        let mut logger = self.logger();
        logger.start_progress(0, 3, "Reading feature data from file")?;
        let mut skipped = Vec::new();
        let map = reader.load_feature_map(&mut skipped)?;
        logger.end_progress(0)?;
        Ok(OmsMapLoad {
            map,
            identifications: data,
            skipped_actions: skipped,
        })
    }

    /// Read a consensus map and its identification data. Source
    /// `load(const std::string&, ConsensusMap&)`.
    ///
    /// # Errors
    ///
    /// As [`OmsFile::load_feature_map`]; [`Error::InvalidValue`] when a
    /// feature handle row precedes every consensus feature or duplicates a
    /// handle of its feature.
    pub fn load_consensus_map(
        &self,
        filename: impl AsRef<Path>,
    ) -> Result<OmsMapLoad<ConsensusMap>> {
        self.load_consensus_map_with(filename, OmsLoadOptions::default())
    }

    /// [`OmsFile::load_consensus_map`] with explicit options.
    ///
    /// # Errors
    ///
    /// As [`OmsFile::load_consensus_map`].
    pub fn load_consensus_map_with(
        &self,
        filename: impl AsRef<Path>,
        options: OmsLoadOptions,
    ) -> Result<OmsMapLoad<ConsensusMap>> {
        let (db, version) = Self::open(filename.as_ref())?;
        let mut data = IdentificationData::new()?;
        let mut reader = load::Load::new(&db, version, options);
        reader.load_identifications(&mut data, &mut || Ok(()))?;
        let mut logger = self.logger();
        logger.start_progress(0, 4, "Reading feature data from file")?;
        let mut skipped = Vec::new();
        let map = reader.load_consensus_map(&mut skipped)?;
        logger.end_progress(0)?;
        Ok(OmsMapLoad {
            map,
            identifications: data,
            skipped_actions: skipped,
        })
    }
}

/// A temporary path beside `filename`, unique within this process.
fn temporary_path(filename: &Path) -> PathBuf {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let serial = NEXT.fetch_add(1, Ordering::Relaxed);
    let name = filename
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| "oms".into());
    filename.with_file_name(format!(".{name}.{}.{serial}.tmp", std::process::id()))
}

/// Create a database in a temporary file, configured as the source's
/// `OMSFileStore` constructor configures it (foreign keys on, synchronous
/// writes and the rollback journal off), fill it in one transaction, and
/// rename it over `filename`.
fn write_database(filename: &Path, fill: impl FnOnce(&Connection) -> Result<()>) -> Result<()> {
    let temporary = temporary_path(filename);
    let _ = std::fs::remove_file(&temporary);
    let result = (|| {
        let db = Connection::open_with_flags(
            &temporary,
            OpenFlags::SQLITE_OPEN_READ_WRITE
                | OpenFlags::SQLITE_OPEN_CREATE
                | OpenFlags::SQLITE_OPEN_NO_MUTEX,
        )
        .map_err(|error| sql_error("error creating OMS file", error))?;
        db.execute_batch("PRAGMA foreign_keys = ON; PRAGMA synchronous = OFF;")
            .map_err(|error| sql_error("error configuring OMS file", error))?;
        db.query_row("PRAGMA journal_mode = OFF", [], |_| Ok(()))
            .map_err(|error| sql_error("error configuring OMS file", error))?;
        db.execute_batch("BEGIN")
            .map_err(|error| sql_error("error starting transaction", error))?;
        fill(&db)?;
        db.execute_batch("COMMIT")
            .map_err(|error| sql_error("error committing transaction", error))?;
        db.close()
            .map_err(|(_, error)| sql_error("error closing OMS file", error))
    })();
    match result {
        Ok(()) => std::fs::rename(&temporary, filename).map_err(|error| {
            let _ = std::fs::remove_file(&temporary);
            Error::Io(error)
        }),
        Err(error) => {
            let _ = std::fs::remove_file(&temporary);
            Err(error)
        }
    }
}

/// Source `exportTableToJSON_`, appending one table's array to `out`.
fn export_table(connection: &Connection, table: &str, out: &mut JsonWriter) -> Result<()> {
    let context = format!("error exporting table {table}");
    let mut sql = format!("SELECT * FROM {}", quoted(table));
    if let Some(order) = order_by(table) {
        sql.push_str(" ORDER BY ");
        sql.push_str(order);
    }
    let mut query = connection
        .prepare(&sql)
        .map_err(|error| sql_error(&context, error))?;
    // The member order of a row object is the sorted column names; a repeated
    // name keeps its last value, as assignment into the source's std::map does.
    let names: Vec<String> = query.column_names().into_iter().map(String::from).collect();
    let mut order: Vec<usize> = (0..names.len()).collect();
    order.sort_by(|a, b| names[*a].cmp(&names[*b]).then(b.cmp(a)));
    order.dedup_by(|later, earlier| names[*later] == names[*earlier]);
    let mut rows = query
        .query([])
        .map_err(|error| sql_error(&context, error))?;
    let mut count = 0usize;
    out.raw("[")?;
    while let Some(row) = rows.next().map_err(|error| sql_error(&context, error))? {
        out.raw(if count == 0 { "\n" } else { ",\n" })?;
        count += 1;
        out.indent(2)?;
        out.raw("{")?;
        for (position, column) in order.iter().enumerate() {
            out.raw(if position == 0 { "\n" } else { ",\n" })?;
            out.indent(3)?;
            out.string(&names[*column])?;
            out.raw(": ")?;
            match row
                .get_ref(*column)
                .map_err(|error| sql_error(&context, error))?
            {
                ValueRef::Null => out.string("")?,
                ValueRef::Integer(value) => out.raw(&value.to_string())?,
                ValueRef::Real(value) => out.raw(&json_number(value))?,
                ValueRef::Text(bytes) | ValueRef::Blob(bytes) => {
                    let text = std::str::from_utf8(bytes).map_err(|_| {
                        Error::InvalidValue(format!(
                            "OMS table {table} holds a cell that is not UTF-8 text"
                        ))
                    })?;
                    out.string(text)?;
                }
            }
        }
        if !order.is_empty() {
            out.raw("\n")?;
            out.indent(2)?;
        }
        out.raw("}")?;
    }
    if count > 0 {
        out.raw("\n")?;
        out.indent(1)?;
    }
    out.raw("]")
}

/// A number as `nlohmann::json` prints a double: the shortest digits that
/// round-trip, in fixed notation for decimal exponents in `(-4, 15]` (with
/// `.0` appended to an integral value) and in scientific notation otherwise,
/// the exponent signed and at least two digits long. Non-finite values print
/// `null`, as the library does.
fn json_number(value: f64) -> String {
    if !value.is_finite() {
        return "null".into();
    }
    if value == 0.0 {
        return if value.is_sign_negative() {
            "-0.0"
        } else {
            "0.0"
        }
        .into();
    }
    // Rust's `{:e}` is the shortest round-trip representation.
    let scientific = format!("{:e}", value.abs());
    let (mantissa, exponent) = scientific.split_once('e').unwrap_or((&scientific, "0"));
    let exponent: i32 = exponent.parse().unwrap_or(0);
    let digits: String = mantissa.chars().filter(char::is_ascii_digit).collect();
    let k = digits.len() as i32;
    // `n` is the position of the decimal point relative to the digit string.
    let n = exponent + 1;
    let sign = if value.is_sign_negative() { "-" } else { "" };
    let body = if k <= n && n <= 15 {
        format!("{digits}{}.0", "0".repeat((n - k) as usize))
    } else if 0 < n && n <= 15 {
        format!("{}.{}", &digits[..n as usize], &digits[n as usize..])
    } else if -4 < n && n <= 0 {
        format!("0.{}{digits}", "0".repeat((-n) as usize))
    } else {
        let exponent = n - 1;
        let exponent_text = format!(
            "{}{:02}",
            if exponent < 0 { '-' } else { '+' },
            exponent.abs()
        );
        if k == 1 {
            format!("{digits}e{exponent_text}")
        } else {
            format!("{}.{}e{exponent_text}", &digits[..1], &digits[1..])
        }
    };
    format!("{sign}{body}")
}

/// A bounded JSON text accumulator.
#[derive(Default)]
struct JsonWriter {
    text: String,
}

impl JsonWriter {
    fn raw(&mut self, text: &str) -> Result<()> {
        if self.text.len().saturating_add(text.len()) > MAX_JSON_BYTES {
            return Err(limit());
        }
        self.text.push_str(text);
        Ok(())
    }

    fn indent(&mut self, level: usize) -> Result<()> {
        self.raw(&" ".repeat(4 * level))
    }

    /// A JSON string as `nlohmann::json::dump` escapes it: `"` and `\`, the
    /// short escapes for backspace, form feed, newline, carriage return and
    /// tab, `\u00xx` (lower-case hex) for the other control characters, and
    /// every other character verbatim.
    fn string(&mut self, text: &str) -> Result<()> {
        let mut escaped = String::with_capacity(text.len() + 2);
        escaped.push('"');
        for c in text.chars() {
            match c {
                '"' => escaped.push_str("\\\""),
                '\\' => escaped.push_str("\\\\"),
                '\u{08}' => escaped.push_str("\\b"),
                '\u{0c}' => escaped.push_str("\\f"),
                '\n' => escaped.push_str("\\n"),
                '\r' => escaped.push_str("\\r"),
                '\t' => escaped.push_str("\\t"),
                c if (c as u32) < 0x20 => escaped.push_str(&format!("\\u{:04x}", c as u32)),
                c => escaped.push(c),
            }
        }
        escaped.push('"');
        self.raw(&escaped)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numbers_follow_the_nlohmann_format() {
        assert_eq!(json_number(0.3), "0.3");
        assert_eq!(json_number(1.0), "1.0");
        assert_eq!(json_number(-2.5), "-2.5");
        assert_eq!(json_number(1234.5678), "1234.5678");
        assert_eq!(json_number(0.0001), "0.0001");
        assert_eq!(json_number(0.00001), "1e-05");
        assert_eq!(json_number(1.5e-7), "1.5e-07");
        assert_eq!(json_number(1e14), "100000000000000.0");
        assert_eq!(json_number(1e15), "1e+15");
        assert_eq!(json_number(1.25e16), "1.25e+16");
        assert_eq!(
            json_number(123_456_789_012_345_680.0),
            "1.2345678901234568e+17"
        );
        assert_eq!(json_number(f64::NAN), "null");
        assert_eq!(json_number(0.0), "0.0");
    }

    #[test]
    fn strings_escape_like_nlohmann() {
        let mut out = JsonWriter::default();
        out.string("a\"b\\c\nd\u{01}é/").unwrap();
        assert_eq!(out.text, "\"a\\\"b\\\\c\\nd\\u0001é/\"");
    }

    #[test]
    fn order_by_matches_the_source_table() {
        assert_eq!(order_by("version"), None);
        assert_eq!(order_by("ID_ScoreType"), Some("id"));
        assert_eq!(order_by("ID_Observation_MetaInfo"), Some("parent_id, name"));
        assert_eq!(
            order_by("ID_ObservationMatch_AppliedProcessingStep"),
            Some("parent_id, processing_step_order, score_type_id")
        );
        assert_eq!(order_by("FEAT_MapMetaData"), Some("unique_id"));
    }
}
