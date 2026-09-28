// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Hendrik Weisser, OpenMS Rust contributors $

//! Writing OMS files: the port of `FORMAT/OMSFileStore.h` (`Internal::OMSFileStore`).
//!
//! Every table, column declaration, bound storage class and insertion order
//! follows `OMSFileStore.cpp` at the pinned revision, so that the JSON export
//! of a Rust-written file equals the JSON export of the file the Release C++
//! build writes from the same data (tier 1, `tests/oms_file.rs`).

use super::{OmsStoreOptions, SCHEMA_VERSION, sql_error};
use crate::format::mztab::meta_value_text;
use crate::identification::graph::{
    AdductId, AppliedProcessingStep, CompoundId, IdentificationData, IdentifiedMolecule,
    InputFileId, MoleculeType, ObservationId, ObservationMatchId, OligoId, ParentGroupSetId,
    ParentId, ParentMatches, PeptideId, ProcessingSoftwareId, ProcessingStepId, ScoreTypeId,
    ScoredProcessingResult, SearchParamId,
};
use crate::kernel::{ConsensusMap, Feature, FeatureMap};
use crate::metadata::{CVTerm, DataProcessing, MetaInfo, MetaValue, MetaValueData};
use crate::{Error, Result};
use rusqlite::types::Value;
use rusqlite::{Connection, Statement, params_from_iter};
use std::collections::{BTreeMap, BTreeSet};

/// Source `DataValue::DataType` index plus one, the `DataValue_DataType` key.
pub(super) fn data_type_id(value: &MetaValue) -> Option<i64> {
    match value.data() {
        MetaValueData::String(_) => Some(1),
        MetaValueData::Integer(_) => Some(2),
        MetaValueData::Float(_) => Some(3),
        MetaValueData::StringList(_) => Some(4),
        MetaValueData::IntegerList(_) => Some(5),
        MetaValueData::FloatList(_) => Some(6),
        MetaValueData::Empty => None,
    }
}

/// Source `EnzymaticDigestion::NamesOfSpecificity`; `None` is `SPEC_UNKNOWN`.
pub(super) fn specificity_name(
    value: Option<crate::chemistry::DigestionSpecificity>,
) -> &'static str {
    match value {
        Some(specificity) => specificity.name(),
        None => "unknown",
    }
}

/// Source `IdentificationData::MoleculeType` index plus one, the `ID_MoleculeType` key.
pub(super) fn molecule_type_id(value: MoleculeType) -> i64 {
    match value {
        MoleculeType::Protein => 1,
        MoleculeType::Compound => 2,
        MoleculeType::RNA => 3,
    }
}

/// The source's unset `DateTime::get()` text.
pub(super) const NULL_DATE_TIME: &str = "0000-00-00 00:00:00";

fn text(value: impl Into<String>) -> Value {
    Value::Text(value.into())
}

fn int(value: impl Into<i64>) -> Value {
    Value::Integer(value.into())
}

fn size_value(value: usize, what: &str) -> Result<Value> {
    // The source binds `uint32_t(value)`, truncating; the port refuses instead.
    u32::try_from(value)
        .map(|value| Value::Integer(value.into()))
        .map_err(|_| Error::InvalidValue(format!("OMS store: {what} exceeds 32 bits")))
}

/// Execute an insert that must modify exactly one row (source
/// `execWithExceptionAndReset(query, 1, ...)`).
fn insert(statement: &mut Statement<'_>, values: &[Value]) -> Result<()> {
    let changed = statement
        .execute(params_from_iter(values.iter()))
        .map_err(|error| sql_error("error inserting data", error))?;
    if changed != 1 {
        return Err(Error::Io(std::io::Error::other(
            "error inserting data: unexpected number of modified rows",
        )));
    }
    Ok(())
}

/// Database keys of the stored records (source `*_keys_` maps).
#[derive(Default)]
struct Keys {
    score_types: BTreeMap<ScoreTypeId, i64>,
    input_files: BTreeMap<InputFileId, i64>,
    softwares: BTreeMap<ProcessingSoftwareId, i64>,
    steps: BTreeMap<ProcessingStepId, i64>,
    search_params: BTreeMap<SearchParamId, i64>,
    observations: BTreeMap<ObservationId, i64>,
    parents: BTreeMap<ParentId, i64>,
    groupings: BTreeMap<ParentGroupSetId, i64>,
    compounds: BTreeMap<CompoundId, i64>,
    peptides: BTreeMap<PeptideId, i64>,
    oligos: BTreeMap<OligoId, i64>,
    adducts: BTreeMap<AdductId, i64>,
    matches: BTreeMap<ObservationMatchId, i64>,
}

fn key<K: Ord + Copy>(keys: &BTreeMap<K, i64>, id: K) -> Result<i64> {
    // The source's `std::map::operator[]` inserts key 0 for an unknown
    // reference, which then violates a foreign key or silently points at no
    // row; a reference into another graph is an error here.
    keys.get(&id).copied().ok_or_else(|| {
        Error::InvalidValue("OMS store: reference to a record that was not stored".into())
    })
}

/// Writer over one open database connection (source `OMSFileStore`).
pub(super) struct Store<'c> {
    db: &'c Connection,
    options: OmsStoreOptions,
    keys: Keys,
    /// Tables created so far; the source asks SQLite (`tableExists`).
    tables: BTreeSet<String>,
}

impl<'c> Store<'c> {
    pub(super) fn new(db: &'c Connection, options: OmsStoreOptions) -> Self {
        Self {
            db,
            options,
            keys: Keys::default(),
            tables: BTreeSet::new(),
        }
    }

    fn exec(&self, sql: &str) -> Result<()> {
        self.db
            .execute_batch(sql)
            .map_err(|error| sql_error("error executing SQL statement", error))
    }

    fn prepare(&self, sql: &str) -> Result<Statement<'c>> {
        self.db
            .prepare(sql)
            .map_err(|error| sql_error("error preparing SQL statement", error))
    }

    /// Source `createTable_(name, definition, may_exist = false)`.
    fn create_table(&mut self, name: &str, definition: &str) -> Result<()> {
        self.exec(&format!("CREATE TABLE {name} ({definition})"))?;
        self.tables.insert(name.to_owned());
        Ok(())
    }

    fn has_table(&self, name: &str) -> bool {
        self.tables.contains(name)
    }

    /// Source `storeVersionAndDate_`. The source writes
    /// `VersionInfo::getVersion()` and the library build time
    /// `VersionInfo::getTime()`; the port writes the core SDK version it
    /// implements and leaves the build time `NULL`, since the crate does not
    /// record one.
    pub(super) fn store_version_and_date(&mut self) -> Result<()> {
        self.create_table(
            "version",
            "OMSFile INT NOT NULL, date TEXT NOT NULL, OpenMS TEXT, build_date TEXT",
        )?;
        let mut query = self.prepare("INSERT INTO version VALUES (?1, datetime('now'), ?2, ?3)")?;
        insert(
            &mut query,
            &[
                int(SCHEMA_VERSION),
                text(crate::CORE_SDK_VERSION),
                Value::Null,
            ],
        )
    }

    fn create_table_molecule_type(&mut self) -> Result<()> {
        self.create_table(
            "ID_MoleculeType",
            "id INTEGER PRIMARY KEY NOT NULL, molecule_type TEXT UNIQUE NOT NULL",
        )?;
        self.exec("INSERT INTO ID_MoleculeType VALUES (1, 'PROTEIN'), (2, 'COMPOUND'), (3, 'RNA')")
    }

    fn create_table_data_value_data_type(&mut self) -> Result<()> {
        self.create_table(
            "DataValue_DataType",
            "id INTEGER PRIMARY KEY NOT NULL, data_type TEXT UNIQUE NOT NULL",
        )?;
        self.exec(
            "INSERT INTO DataValue_DataType VALUES (1, 'STRING_VALUE'), (2, 'INT_VALUE'), \
             (3, 'DOUBLE_VALUE'), (4, 'STRING_LIST'), (5, 'INT_LIST'), (6, 'DOUBLE_LIST')",
        )
    }

    fn create_table_cv_term(&mut self) -> Result<()> {
        self.create_table(
            "CVTerm",
            "id INTEGER PRIMARY KEY NOT NULL, accession TEXT UNIQUE, name TEXT NOT NULL, \
             cv_identifier_ref TEXT, UNIQUE (accession, name)",
        )
    }

    /// Source `storeCVTerm_`: insert or reuse a CV term row, with `NULL` for
    /// an empty accession.
    fn store_cv_term(&mut self, term: &CVTerm) -> Result<i64> {
        let accession = if term.accession.is_empty() {
            Value::Null
        } else {
            text(term.accession.as_str())
        };
        let mut query = self.prepare("INSERT OR IGNORE INTO CVTerm VALUES (NULL, ?1, ?2, ?3)")?;
        let changed = query
            .execute(params_from_iter(
                [
                    accession.clone(),
                    text(term.name.as_str()),
                    text(term.cv_ref.as_str()),
                ]
                .iter(),
            ))
            .map_err(|error| sql_error("error inserting data", error))?;
        if changed == 1 {
            return Ok(self.db.last_insert_rowid());
        }
        let mut lookup =
            self.prepare("SELECT id FROM CVTerm WHERE accession = ?1 AND name = ?2")?;
        let mut rows = lookup
            .query(params_from_iter(
                [accession, text(term.name.as_str())].iter(),
            ))
            .map_err(|error| sql_error("error querying database", error))?;
        match rows
            .next()
            .map_err(|error| sql_error("error querying database", error))?
        {
            Some(row) => row
                .get::<_, i64>(0)
                .map_err(|error| sql_error("error querying database", error)),
            None => Err(Error::Io(std::io::Error::other(
                "error querying database: CV term neither inserted nor found",
            ))),
        }
    }

    /// Source `createTableMetaInfo_`.
    fn create_table_meta_info(&mut self, parent_table: &str, key_column: &str) -> Result<()> {
        if !self.has_table("DataValue_DataType") {
            self.create_table_data_value_data_type()?;
        }
        let table = format!("{parent_table}_MetaInfo");
        if self.has_table(&table) {
            // The source creates this table once per call site and fails with
            // "table already exists" on a second (see `store_scored_results`).
            return Ok(());
        }
        self.create_table(
            &table,
            &format!(
                "parent_id INTEGER NOT NULL, name TEXT NOT NULL, data_type_id INTEGER, \
                 value TEXT, FOREIGN KEY (parent_id) REFERENCES {parent_table} ({key_column}), \
                 FOREIGN KEY (data_type_id) REFERENCES DataValue_DataType (id), \
                 PRIMARY KEY (parent_id, name)"
            ),
        )
    }

    /// Source `storeMetaInfo_`: one row per key, the value as
    /// `DataValue::toString()` text and its type (`NULL` for an empty value).
    fn store_meta_info(&self, info: &MetaInfo, parent_table: &str, parent_id: i64) -> Result<()> {
        if info.is_empty() {
            return Ok(());
        }
        let mut query = self.prepare(&format!(
            "INSERT INTO {parent_table}_MetaInfo VALUES (?1, ?2, ?3, ?4)"
        ))?;
        for (name, value) in info {
            insert(
                &mut query,
                &[
                    int(parent_id),
                    text(name.as_str()),
                    data_type_id(value).map_or(Value::Null, Value::Integer),
                    text(meta_value_text(value)),
                ],
            )?;
        }
        Ok(())
    }

    /// Source `storeMetaInfos_`: the meta-info table is created on the first
    /// element with metadata.
    fn store_meta_infos<'a>(
        &mut self,
        items: impl IntoIterator<Item = (i64, &'a MetaInfo)>,
        parent_table: &str,
    ) -> Result<()> {
        for (key, info) in items {
            if info.is_empty() {
                continue;
            }
            self.create_table_meta_info(parent_table, "id")?;
            self.store_meta_info(info, parent_table, key)?;
        }
        Ok(())
    }

    fn create_table_applied_processing_step(&mut self, parent_table: &str) -> Result<()> {
        let table = format!("{parent_table}_AppliedProcessingStep");
        if self.has_table(&table) {
            return Ok(());
        }
        self.create_table(
            &table,
            &format!(
                "parent_id INTEGER NOT NULL, processing_step_id INTEGER, \
                 processing_step_order INTEGER NOT NULL, score_type_id INTEGER, score REAL, \
                 UNIQUE (parent_id, processing_step_id, score_type_id), \
                 FOREIGN KEY (parent_id) REFERENCES {parent_table} (id), \
                 FOREIGN KEY (score_type_id) REFERENCES ID_ScoreType (id), \
                 FOREIGN KEY (processing_step_id) REFERENCES ID_ProcessingStep (id)"
            ),
        )
    }

    /// Source `storeAppliedProcessingStep_`. A step without a processing-step
    /// reference and without scores writes no row, as in the source.
    fn store_applied_processing_step(
        &self,
        step: &AppliedProcessingStep,
        order: usize,
        parent_table: &str,
        parent_id: i64,
    ) -> Result<()> {
        let mut query = self.prepare(&format!(
            "INSERT INTO {parent_table}_AppliedProcessingStep VALUES (?1, ?2, ?3, ?4, ?5)"
        ))?;
        let order = size_value(order, "processing step order")?;
        let step_id = match step.processing_step {
            Some(id) => int(key(&self.keys.steps, id)?),
            None => Value::Null,
        };
        if step.processing_step.is_some() && step.scores.is_empty() {
            insert(
                &mut query,
                &[
                    int(parent_id),
                    step_id.clone(),
                    order.clone(),
                    Value::Null,
                    Value::Null,
                ],
            )?;
        }
        for (&score_type, &score) in &step.scores {
            insert(
                &mut query,
                &[
                    int(parent_id),
                    step_id.clone(),
                    order.clone(),
                    int(key(&self.keys.score_types, score_type)?),
                    Value::Real(score),
                ],
            )?;
        }
        Ok(())
    }

    /// Source `storeScoredProcessingResults_`: processing steps and scores of
    /// every element, then their metadata.
    ///
    /// The source creates `<table>_AppliedProcessingStep` and
    /// `<table>_MetaInfo` without `IF NOT EXISTS` on every call, so it fails
    /// with "table already exists" when two of compounds, peptides and
    /// oligonucleotides (which share `ID_IdentifiedMolecule`) carry steps or
    /// metadata. The port creates each table once and stores all of them.
    fn store_scored_results(
        &mut self,
        items: &[(i64, &ScoredProcessingResult)],
        parent_table: &str,
    ) -> Result<()> {
        for (parent_id, result) in items {
            if result.steps_and_scores.is_empty() {
                continue;
            }
            self.create_table_applied_processing_step(parent_table)?;
            for (index, step) in result.steps_and_scores.iter().enumerate() {
                self.store_applied_processing_step(step, index + 1, parent_table, *parent_id)?;
            }
        }
        self.store_meta_infos(
            items.iter().map(|(id, result)| (*id, &result.metadata)),
            parent_table,
        )
    }

    fn store_score_types(&mut self, data: &IdentificationData) -> Result<()> {
        if data.score_type_count() == 0 {
            return Ok(());
        }
        self.create_table_cv_term()?;
        self.create_table(
            "ID_ScoreType",
            "id INTEGER PRIMARY KEY NOT NULL, cv_term_id INTEGER NOT NULL, \
             higher_better NUMERIC NOT NULL CHECK (higher_better in (0, 1)), \
             FOREIGN KEY (cv_term_id) REFERENCES CVTerm (id)",
        )?;
        for (next, (id, score_type)) in (1i64..).zip(data.score_types()) {
            let cv_id = self.store_cv_term(&score_type.cv_term)?;
            let mut query = self.prepare("INSERT INTO ID_ScoreType VALUES (?1, ?2, ?3)")?;
            insert(
                &mut query,
                &[
                    int(next),
                    int(cv_id),
                    int(i64::from(score_type.higher_better)),
                ],
            )?;
            self.keys.score_types.insert(id, next);
        }
        Ok(())
    }

    fn store_input_files(&mut self, data: &IdentificationData) -> Result<()> {
        if data.input_file_count() == 0 {
            return Ok(());
        }
        self.create_table(
            "ID_InputFile",
            "id INTEGER PRIMARY KEY NOT NULL, name TEXT UNIQUE NOT NULL, \
             experimental_design_id TEXT, primary_files TEXT",
        )?;
        let mut query = self.prepare("INSERT INTO ID_InputFile VALUES (?1, ?2, ?3, ?4)")?;
        for (next, (id, input)) in (1i64..).zip(data.input_files()) {
            let primary_files: Vec<&str> = input.primary_files.iter().map(String::as_str).collect();
            insert(
                &mut query,
                &[
                    int(next),
                    text(input.name.as_str()),
                    text(input.experimental_design_id.as_str()),
                    text(primary_files.join(",")),
                ],
            )?;
            self.keys.input_files.insert(id, next);
        }
        Ok(())
    }

    fn store_processing_softwares(&mut self, data: &IdentificationData) -> Result<()> {
        if data.processing_software_count() == 0 {
            return Ok(());
        }
        self.create_table(
            "ID_ProcessingSoftware",
            "id INTEGER PRIMARY KEY NOT NULL, name TEXT NOT NULL, version TEXT, \
             UNIQUE (name, version)",
        )?;
        let mut any_scores = false;
        let mut next = 1i64;
        {
            let mut query =
                self.prepare("INSERT INTO ID_ProcessingSoftware VALUES (?1, ?2, ?3)")?;
            for (id, software) in data.processing_softwares() {
                any_scores |= !software.assigned_scores.is_empty();
                insert(
                    &mut query,
                    &[
                        int(next),
                        text(software.software.name.as_str()),
                        text(software.software.version.as_str()),
                    ],
                )?;
                self.keys.softwares.insert(id, next);
                next += 1;
            }
        }
        if any_scores {
            self.create_table(
                "ID_ProcessingSoftware_AssignedScore",
                "software_id INTEGER NOT NULL, score_type_id INTEGER NOT NULL, \
                 score_type_order INTEGER NOT NULL, UNIQUE (software_id, score_type_id), \
                 UNIQUE (software_id, score_type_order), \
                 FOREIGN KEY (software_id) REFERENCES ID_ProcessingSoftware (id), \
                 FOREIGN KEY (score_type_id) REFERENCES ID_ScoreType (id)",
            )?;
            let mut query = self
                .prepare("INSERT INTO ID_ProcessingSoftware_AssignedScore VALUES (?1, ?2, ?3)")?;
            for (id, software) in data.processing_softwares() {
                let software_id = key(&self.keys.softwares, id)?;
                for (index, &score_type) in software.assigned_scores.iter().enumerate() {
                    insert(
                        &mut query,
                        &[
                            int(software_id),
                            int(key(&self.keys.score_types, score_type)?),
                            size_value(index + 1, "score type order")?,
                        ],
                    )?;
                }
            }
        }
        Ok(())
    }

    fn store_db_search_params(&mut self, data: &IdentificationData) -> Result<()> {
        if data.db_search_param_count() == 0 {
            return Ok(());
        }
        if !self.has_table("ID_MoleculeType") {
            self.create_table_molecule_type()?;
        }
        self.create_table(
            "ID_DBSearchParam",
            "id INTEGER PRIMARY KEY NOT NULL, molecule_type_id INTEGER NOT NULL, \
             mass_type_average NUMERIC NOT NULL CHECK (mass_type_average in (0, 1)) DEFAULT 0, \
             database TEXT, database_version TEXT, taxonomy TEXT, charges TEXT, \
             fixed_mods TEXT, variable_mods TEXT, precursor_mass_tolerance REAL, \
             fragment_mass_tolerance REAL, \
             precursor_tolerance_ppm NUMERIC NOT NULL CHECK (precursor_tolerance_ppm in (0, 1)) DEFAULT 0, \
             fragment_tolerance_ppm NUMERIC NOT NULL CHECK (fragment_tolerance_ppm in (0, 1)) DEFAULT 0, \
             digestion_enzyme TEXT, enzyme_term_specificity TEXT, missed_cleavages NUMERIC, \
             min_length NUMERIC, max_length NUMERIC, \
             FOREIGN KEY (molecule_type_id) REFERENCES ID_MoleculeType (id)",
        )?;
        let mut query = self.prepare(
            "INSERT INTO ID_DBSearchParam VALUES \
             (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14, ?15, ?16, ?17, ?18)",
        )?;
        for (next, (id, param)) in (1i64..).zip(data.db_search_params()) {
            let charges: Vec<String> = param.charges.iter().map(i32::to_string).collect();
            let fixed: Vec<&str> = param.fixed_mods.iter().map(String::as_str).collect();
            let variable: Vec<&str> = param.variable_mods.iter().map(String::as_str).collect();
            let enzyme = match &param.digestion_enzyme {
                Some(crate::identification::graph::GraphEnzyme::Protein(enzyme)) => {
                    text(enzyme.name())
                }
                Some(crate::identification::graph::GraphEnzyme::RNA(enzyme)) => text(enzyme.name()),
                None => Value::Null,
            };
            insert(
                &mut query,
                &[
                    int(next),
                    int(molecule_type_id(param.molecule_type)),
                    int(match param.mass_type {
                        crate::identification::graph::MassType::Monoisotopic => 0,
                        crate::identification::graph::MassType::Average => 1,
                    }),
                    text(param.database.as_str()),
                    text(param.database_version.as_str()),
                    text(param.taxonomy.as_str()),
                    text(charges.join(",")),
                    text(fixed.join(",")),
                    text(variable.join(",")),
                    Value::Real(param.precursor_mass_tolerance),
                    Value::Real(param.fragment_mass_tolerance),
                    int(i64::from(param.precursor_tolerance_ppm)),
                    int(i64::from(param.fragment_tolerance_ppm)),
                    enzyme,
                    text(specificity_name(param.enzyme_term_specificity)),
                    size_value(param.missed_cleavages, "missed cleavages")?,
                    size_value(param.min_length, "minimum length")?,
                    size_value(param.max_length, "maximum length")?,
                ],
            )?;
            self.keys.search_params.insert(id, next);
        }
        Ok(())
    }

    fn store_processing_steps(&mut self, data: &IdentificationData) -> Result<()> {
        if data.processing_step_count() == 0 {
            return Ok(());
        }
        self.create_table(
            "ID_ProcessingStep",
            "id INTEGER PRIMARY KEY NOT NULL, software_id INTEGER NOT NULL, date_time TEXT, \
             search_param_id INTEGER, \
             FOREIGN KEY (search_param_id) REFERENCES ID_DBSearchParam (id)",
        )?;
        let mut any_input_files = false;
        let mut next = 1i64;
        {
            let mut query =
                self.prepare("INSERT INTO ID_ProcessingStep VALUES (?1, ?2, ?3, ?4)")?;
            for (id, step) in data.processing_steps() {
                any_input_files |= !step.input_files.is_empty();
                let search = match data.search_param_for_step(id)? {
                    Some(search) => int(key(&self.keys.search_params, search)?),
                    None => Value::Null,
                };
                let date_time = match step.date_time {
                    Some(time) => time.to_string(),
                    None => NULL_DATE_TIME.to_owned(),
                };
                insert(
                    &mut query,
                    &[
                        int(next),
                        int(key(&self.keys.softwares, step.software)?),
                        text(date_time),
                        search,
                    ],
                )?;
                self.keys.steps.insert(id, next);
                next += 1;
            }
        }
        if any_input_files {
            self.create_table(
                "ID_ProcessingStep_InputFile",
                "processing_step_id INTEGER NOT NULL, input_file_id INTEGER NOT NULL, \
                 FOREIGN KEY (processing_step_id) REFERENCES ID_ProcessingStep (id), \
                 FOREIGN KEY (input_file_id) REFERENCES ID_InputFile (id), \
                 UNIQUE (processing_step_id, input_file_id)",
            )?;
            let mut query =
                self.prepare("INSERT INTO ID_ProcessingStep_InputFile VALUES (?1, ?2)")?;
            for (id, step) in data.processing_steps() {
                let step_id = key(&self.keys.steps, id)?;
                for &input in &step.input_files {
                    insert(
                        &mut query,
                        &[int(step_id), int(key(&self.keys.input_files, input)?)],
                    )?;
                }
            }
        }
        let items: Vec<(i64, &MetaInfo)> = data
            .processing_steps()
            .map(|(id, step)| Ok((key(&self.keys.steps, id)?, &step.metadata)))
            .collect::<Result<_>>()?;
        self.store_meta_infos(items, "ID_ProcessingStep")
    }

    fn store_observations(&mut self, data: &IdentificationData) -> Result<()> {
        if data.observation_count() == 0 {
            return Ok(());
        }
        self.create_table(
            "ID_Observation",
            "id INTEGER PRIMARY KEY NOT NULL, data_id TEXT NOT NULL, \
             input_file_id INTEGER NOT NULL, rt REAL, mz REAL, \
             UNIQUE (data_id, input_file_id), \
             FOREIGN KEY (input_file_id) REFERENCES ID_InputFile (id)",
        )?;
        let mut next = 1i64;
        {
            let mut query =
                self.prepare("INSERT INTO ID_Observation VALUES (?1, ?2, ?3, ?4, ?5)")?;
            for (id, observation) in data.observations() {
                insert(
                    &mut query,
                    &[
                        int(next),
                        text(observation.data_id.as_str()),
                        int(key(&self.keys.input_files, observation.input_file)?),
                        observation.rt.map_or(Value::Null, Value::Real),
                        observation.mz.map_or(Value::Null, Value::Real),
                    ],
                )?;
                self.keys.observations.insert(id, next);
                next += 1;
            }
        }
        let items: Vec<(i64, &MetaInfo)> = data
            .observations()
            .map(|(id, observation)| Ok((key(&self.keys.observations, id)?, &observation.metadata)))
            .collect::<Result<_>>()?;
        self.store_meta_infos(items, "ID_Observation")
    }

    fn store_parent_sequences(&mut self, data: &IdentificationData) -> Result<()> {
        if data.parent_count() == 0 {
            return Ok(());
        }
        if !self.has_table("ID_MoleculeType") {
            self.create_table_molecule_type()?;
        }
        self.create_table(
            "ID_ParentSequence",
            "id INTEGER PRIMARY KEY NOT NULL, accession TEXT UNIQUE NOT NULL, \
             molecule_type_id INTEGER NOT NULL, sequence TEXT, description TEXT, coverage REAL, \
             is_decoy NUMERIC NOT NULL CHECK (is_decoy in (0, 1)) DEFAULT 0, \
             FOREIGN KEY (molecule_type_id) REFERENCES ID_MoleculeType (id)",
        )?;
        let mut next = 1i64;
        {
            let mut query =
                self.prepare("INSERT INTO ID_ParentSequence VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)")?;
            for (id, parent) in data.parents() {
                insert(
                    &mut query,
                    &[
                        int(next),
                        text(parent.accession.as_str()),
                        int(molecule_type_id(parent.molecule_type)),
                        text(parent.sequence.as_str()),
                        text(parent.description.as_str()),
                        Value::Real(parent.coverage),
                        int(i64::from(parent.is_decoy)),
                    ],
                )?;
                self.keys.parents.insert(id, next);
                next += 1;
            }
        }
        let items: Vec<(i64, &ScoredProcessingResult)> = data
            .parents()
            .map(|(id, parent)| Ok((key(&self.keys.parents, id)?, &parent.result)))
            .collect::<Result<_>>()?;
        self.store_scored_results(&items, "ID_ParentSequence")
    }

    /// Source `storeParentGroupSets_`. A group with several scores is written
    /// once per score under the same `id`, which the table's primary key
    /// refuses, in the source as here.
    fn store_parent_group_sets(&mut self, data: &IdentificationData) -> Result<()> {
        if data.parent_group_set_count() == 0 {
            return Ok(());
        }
        self.create_table(
            "ID_ParentGroupSet",
            "id INTEGER PRIMARY KEY NOT NULL, label TEXT UNIQUE",
        )?;
        self.create_table(
            "ID_ParentGroup",
            "id INTEGER PRIMARY KEY NOT NULL, grouping_id INTEGER NOT NULL, \
             score_type_id INTEGER, score REAL, UNIQUE (id, score_type_id), \
             FOREIGN KEY (grouping_id) REFERENCES ID_ParentGroupSet (id)",
        )?;
        self.create_table(
            "ID_ParentGroup_ParentSequence",
            "group_id INTEGER NOT NULL, parent_id INTEGER NOT NULL, \
             UNIQUE (group_id, parent_id), \
             FOREIGN KEY (group_id) REFERENCES ID_ParentGroup (id), \
             FOREIGN KEY (parent_id) REFERENCES ID_ParentSequence (id)",
        )?;
        let mut grouping_id = 1i64;
        let mut group_id = 1i64;
        {
            let mut query_grouping =
                self.prepare("INSERT INTO ID_ParentGroupSet VALUES (?1, ?2)")?;
            let mut query_group =
                self.prepare("INSERT INTO ID_ParentGroup VALUES (?1, ?2, ?3, ?4)")?;
            let mut query_parent =
                self.prepare("INSERT INTO ID_ParentGroup_ParentSequence VALUES (?1, ?2)")?;
            for (id, grouping) in data.parent_group_sets() {
                insert(
                    &mut query_grouping,
                    &[int(grouping_id), text(grouping.label.as_str())],
                )?;
                for group in &grouping.groups {
                    if group.scores.is_empty() {
                        insert(
                            &mut query_group,
                            &[int(group_id), int(grouping_id), Value::Null, Value::Null],
                        )?;
                    }
                    for (&score_type, &score) in &group.scores {
                        insert(
                            &mut query_group,
                            &[
                                int(group_id),
                                int(grouping_id),
                                int(key(&self.keys.score_types, score_type)?),
                                Value::Real(score),
                            ],
                        )?;
                    }
                    for &parent in &group.parent_refs {
                        insert(
                            &mut query_parent,
                            &[int(group_id), int(key(&self.keys.parents, parent)?)],
                        )?;
                    }
                    group_id += 1;
                }
                self.keys.groupings.insert(id, grouping_id);
                grouping_id += 1;
            }
        }
        let items: Vec<(i64, &ScoredProcessingResult)> = data
            .parent_group_sets()
            .map(|(id, grouping)| Ok((key(&self.keys.groupings, id)?, &grouping.result)))
            .collect::<Result<_>>()?;
        self.store_scored_results(&items, "ID_ParentGroupSet")
    }

    fn create_table_identified_molecule(&mut self) -> Result<()> {
        if !self.has_table("ID_MoleculeType") {
            self.create_table_molecule_type()?;
        }
        self.create_table(
            "ID_IdentifiedMolecule",
            "id INTEGER PRIMARY KEY NOT NULL, molecule_type_id INTEGER NOT NULL, \
             identifier TEXT NOT NULL, UNIQUE (molecule_type_id, identifier), \
             FOREIGN KEY (molecule_type_id) REFERENCES ID_MoleculeType (id)",
        )
    }

    fn store_identified_compounds(&mut self, data: &IdentificationData) -> Result<()> {
        if data.compound_count() == 0 {
            return Ok(());
        }
        if !self.has_table("ID_IdentifiedMolecule") {
            self.create_table_identified_molecule()?;
        }
        self.create_table(
            "ID_IdentifiedCompound",
            "molecule_id INTEGER UNIQUE NOT NULL , formula TEXT, name TEXT, smile TEXT, \
             inchi TEXT, FOREIGN KEY (molecule_id) REFERENCES ID_IdentifiedMolecule (id)",
        )?;
        let mut next = 1i64;
        {
            let mut query_molecule =
                self.prepare("INSERT INTO ID_IdentifiedMolecule VALUES (?1, ?2, ?3)")?;
            let mut query_compound =
                self.prepare("INSERT INTO ID_IdentifiedCompound VALUES (?1, ?2, ?3, ?4, ?5)")?;
            for (id, compound) in data.compounds() {
                insert(
                    &mut query_molecule,
                    &[
                        int(next),
                        int(molecule_type_id(MoleculeType::Compound)),
                        text(compound.identifier.as_str()),
                    ],
                )?;
                let smile = if self.options.compound_smile_from_name {
                    compound.name.as_str()
                } else {
                    compound.smile.as_str()
                };
                insert(
                    &mut query_compound,
                    &[
                        int(next),
                        text(compound.formula.to_string()),
                        text(compound.name.as_str()),
                        text(smile),
                        text(compound.inchi.as_str()),
                    ],
                )?;
                self.keys.compounds.insert(id, next);
                next += 1;
            }
        }
        let items: Vec<(i64, &ScoredProcessingResult)> = data
            .compounds()
            .map(|(id, compound)| Ok((key(&self.keys.compounds, id)?, &compound.result)))
            .collect::<Result<_>>()?;
        self.store_scored_results(&items, "ID_IdentifiedMolecule")
    }

    fn store_identified_sequences(&mut self, data: &IdentificationData) -> Result<()> {
        if data.peptide_count() == 0 && data.oligo_count() == 0 {
            return Ok(());
        }
        if !self.has_table("ID_IdentifiedMolecule") {
            self.create_table_identified_molecule()?;
        }
        let mut any_parent_matches = false;
        let mut next = i64::try_from(data.compound_count())
            .map_err(|_| Error::InvalidValue("OMS store: too many compounds".into()))?
            + 1;
        {
            let mut query =
                self.prepare("INSERT INTO ID_IdentifiedMolecule VALUES (?1, ?2, ?3)")?;
            for (id, peptide) in data.peptides() {
                any_parent_matches |= !peptide.parent_matches.is_empty();
                insert(
                    &mut query,
                    &[
                        int(next),
                        int(molecule_type_id(MoleculeType::Protein)),
                        text(peptide.sequence.to_string()),
                    ],
                )?;
                self.keys.peptides.insert(id, next);
                next += 1;
            }
        }
        let items: Vec<(i64, &ScoredProcessingResult)> = data
            .peptides()
            .map(|(id, peptide)| Ok((key(&self.keys.peptides, id)?, &peptide.result)))
            .collect::<Result<_>>()?;
        self.store_scored_results(&items, "ID_IdentifiedMolecule")?;
        {
            let mut query =
                self.prepare("INSERT INTO ID_IdentifiedMolecule VALUES (?1, ?2, ?3)")?;
            for (id, oligo) in data.oligos() {
                any_parent_matches |= !oligo.parent_matches.is_empty();
                insert(
                    &mut query,
                    &[
                        int(next),
                        int(molecule_type_id(MoleculeType::RNA)),
                        text(oligo.sequence.to_string()),
                    ],
                )?;
                self.keys.oligos.insert(id, next);
                next += 1;
            }
        }
        let items: Vec<(i64, &ScoredProcessingResult)> = data
            .oligos()
            .map(|(id, oligo)| Ok((key(&self.keys.oligos, id)?, &oligo.result)))
            .collect::<Result<_>>()?;
        self.store_scored_results(&items, "ID_IdentifiedMolecule")?;
        if any_parent_matches {
            self.create_table(
                "ID_ParentMatch",
                "molecule_id INTEGER NOT NULL, parent_id INTEGER NOT NULL, start_pos NUMERIC, \
                 end_pos NUMERIC, left_neighbor TEXT, right_neighbor TEXT, \
                 UNIQUE (molecule_id, parent_id, start_pos, end_pos), \
                 FOREIGN KEY (parent_id) REFERENCES ID_ParentSequence (id), \
                 FOREIGN KEY (molecule_id) REFERENCES ID_IdentifiedMolecule (id)",
            )?;
            for (id, peptide) in data.peptides() {
                self.store_parent_matches(&peptide.parent_matches, key(&self.keys.peptides, id)?)?;
            }
            for (id, oligo) in data.oligos() {
                self.store_parent_matches(&oligo.parent_matches, key(&self.keys.oligos, id)?)?;
            }
        }
        Ok(())
    }

    /// Source `storeParentMatches_`; an unknown position is `NULL`.
    fn store_parent_matches(&self, matches: &ParentMatches, molecule_id: i64) -> Result<()> {
        let mut query =
            self.prepare("INSERT INTO ID_ParentMatch VALUES (?1, ?2, ?3, ?4, ?5, ?6)")?;
        for (&parent, set) in matches {
            let parent_id = key(&self.keys.parents, parent)?;
            for item in set {
                let start = match item.start_pos {
                    Some(position) => size_value(position, "parent match start")?,
                    None => Value::Null,
                };
                let end = match item.end_pos {
                    Some(position) => size_value(position, "parent match end")?,
                    None => Value::Null,
                };
                insert(
                    &mut query,
                    &[
                        int(molecule_id),
                        int(parent_id),
                        start,
                        end,
                        text(item.left_neighbor.as_str()),
                        text(item.right_neighbor.as_str()),
                    ],
                )?;
            }
        }
        Ok(())
    }

    fn store_adducts(&mut self, data: &IdentificationData) -> Result<()> {
        if data.adduct_count() == 0 {
            return Ok(());
        }
        self.create_table(
            "AdductInfo",
            "id INTEGER PRIMARY KEY NOT NULL, name TEXT, formula TEXT NOT NULL, \
             charge INTEGER NOT NULL, \
             mol_multiplier INTEGER NOT NULL CHECK (mol_multiplier > 0) DEFAULT 1, \
             UNIQUE (formula, charge)",
        )?;
        let mut query = self.prepare("INSERT INTO AdductInfo VALUES (?1, ?2, ?3, ?4, ?5)")?;
        for (next, (id, adduct)) in (1i64..).zip(data.adducts()) {
            insert(
                &mut query,
                &[
                    int(next),
                    text(adduct.name()),
                    text(adduct.empirical_formula().to_string()),
                    int(adduct.charge()),
                    int(adduct.mol_multiplier()),
                ],
            )?;
            self.keys.adducts.insert(id, next);
        }
        Ok(())
    }

    /// Source `getDatabaseKey_`.
    fn molecule_key(&self, molecule: IdentifiedMolecule) -> Result<i64> {
        match molecule {
            IdentifiedMolecule::Peptide(id) => key(&self.keys.peptides, id),
            IdentifiedMolecule::Compound(id) => key(&self.keys.compounds, id),
            IdentifiedMolecule::Oligo(id) => key(&self.keys.oligos, id),
        }
    }

    fn store_observation_matches(&mut self, data: &IdentificationData) -> Result<()> {
        if data.observation_match_count() == 0 {
            return Ok(());
        }
        let mut definition = String::from(
            "id INTEGER PRIMARY KEY NOT NULL, identified_molecule_id INTEGER NOT NULL, \
             observation_id INTEGER NOT NULL, adduct_id INTEGER, charge INTEGER, \
             FOREIGN KEY (identified_molecule_id) REFERENCES ID_IdentifiedMolecule (id), \
             FOREIGN KEY (observation_id) REFERENCES ID_Observation (id)",
        );
        if self.has_table("AdductInfo") {
            definition.push_str(", FOREIGN KEY (adduct_id) REFERENCES AdductInfo (id)");
        }
        self.create_table("ID_ObservationMatch", &definition)?;
        let mut any_peak_annotations = false;
        let mut next = 1i64;
        {
            let mut query =
                self.prepare("INSERT INTO ID_ObservationMatch VALUES (?1, ?2, ?3, ?4, ?5)")?;
            for (id, item) in data.observation_matches() {
                any_peak_annotations |= !item.peak_annotations.is_empty();
                let adduct = match item.adduct {
                    Some(adduct) => int(key(&self.keys.adducts, adduct)?),
                    None => Value::Null,
                };
                insert(
                    &mut query,
                    &[
                        int(next),
                        int(self.molecule_key(item.identified_molecule)?),
                        int(key(&self.keys.observations, item.observation)?),
                        adduct,
                        int(item.charge),
                    ],
                )?;
                self.keys.matches.insert(id, next);
                next += 1;
            }
        }
        let items: Vec<(i64, &ScoredProcessingResult)> = data
            .observation_matches()
            .map(|(id, item)| Ok((key(&self.keys.matches, id)?, &item.result)))
            .collect::<Result<_>>()?;
        self.store_scored_results(&items, "ID_ObservationMatch")?;
        if any_peak_annotations {
            self.create_table(
                "ID_ObservationMatch_PeakAnnotation",
                "parent_id INTEGER NOT NULL, processing_step_id INTEGER, peak_annotation TEXT, \
                 peak_charge INTEGER, peak_mz REAL, peak_intensity REAL, \
                 FOREIGN KEY (parent_id) REFERENCES ID_ObservationMatch (id), \
                 FOREIGN KEY (processing_step_id) REFERENCES ID_ProcessingStep (id)",
            )?;
            let mut query = self.prepare(
                "INSERT INTO ID_ObservationMatch_PeakAnnotation VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
            )?;
            for (id, item) in data.observation_matches() {
                let parent_id = key(&self.keys.matches, id)?;
                for (step, annotations) in &item.peak_annotations {
                    let step_id = match step {
                        Some(step) => int(key(&self.keys.steps, *step)?),
                        None => Value::Null,
                    };
                    for annotation in annotations {
                        insert(
                            &mut query,
                            &[
                                int(parent_id),
                                step_id.clone(),
                                text(annotation.annotation.as_str()),
                                int(annotation.charge),
                                Value::Real(annotation.mz),
                                Value::Real(annotation.intensity),
                            ],
                        )?;
                    }
                }
            }
            self.exec(
                "CREATE INDEX PeakAnnotation_parent_id ON ID_ObservationMatch_PeakAnnotation \
                 (parent_id)",
            )?;
        }
        Ok(())
    }

    /// Source `OMSFileStore::store(const IdentificationData&)`, without the
    /// transaction (the caller owns it) and without the version row when
    /// `with_version` is false.
    pub(super) fn store_identifications(
        &mut self,
        data: &IdentificationData,
        progress: &mut dyn FnMut() -> Result<()>,
    ) -> Result<()> {
        self.store_version_and_date()?;
        progress()?;
        self.store_input_files(data)?;
        progress()?;
        self.store_score_types(data)?;
        progress()?;
        self.store_processing_softwares(data)?;
        progress()?;
        self.store_db_search_params(data)?;
        progress()?;
        self.store_processing_steps(data)?;
        progress()?;
        self.store_observations(data)?;
        progress()?;
        self.store_parent_sequences(data)?;
        progress()?;
        self.store_parent_group_sets(data)?;
        progress()?;
        self.store_identified_compounds(data)?;
        progress()?;
        self.store_identified_sequences(data)?;
        progress()?;
        self.store_adducts(data)?;
        progress()?;
        self.store_observation_matches(data)
    }

    // ------------------------------------------------------------------
    // Feature and consensus maps
    // ------------------------------------------------------------------

    /// Source `createTableBaseFeature_`.
    fn create_table_base_feature(
        &mut self,
        with_metainfo: bool,
        with_idmatches: bool,
    ) -> Result<()> {
        // The source always declares the foreign key into
        // ID_IdentifiedMolecule. With foreign keys enforced, SQLite refuses
        // every insert into a table whose parent table does not exist ("no
        // such table: main.ID_IdentifiedMolecule"), so the source cannot store
        // a map whose identification data holds no identified molecule. The
        // port declares the key only when the table exists, as the source
        // itself does for `ID_ObservationMatch.adduct_id`.
        let primary_key = if self.has_table("ID_IdentifiedMolecule") {
            "FOREIGN KEY (primary_molecule_id) REFERENCES ID_IdentifiedMolecule (id), "
        } else {
            ""
        };
        self.create_table(
            "FEAT_BaseFeature",
            &format!(
                "id INTEGER PRIMARY KEY NOT NULL, rt REAL, mz REAL, intensity REAL, \
                 charge INTEGER, width REAL, quality REAL, unique_id INTEGER, \
                 primary_molecule_id INTEGER, subordinate_of INTEGER, {primary_key}\
                 FOREIGN KEY (subordinate_of) REFERENCES FEAT_BaseFeature (id), \
                 CHECK (id > subordinate_of)"
            ),
        )?;
        if with_metainfo {
            self.create_table_meta_info("FEAT_BaseFeature", "id")?;
        }
        if with_idmatches {
            self.create_table(
                "FEAT_ObservationMatch",
                "feature_id INTEGER NOT NULL, observation_match_id INTEGER NOT NULL, \
                 FOREIGN KEY (feature_id) REFERENCES FEAT_BaseFeature (id), \
                 FOREIGN KEY (observation_match_id) REFERENCES ID_ObservationMatch (id)",
            )?;
        }
        Ok(())
    }

    /// Source `storeBaseFeature_`. `parent_id` is the containing feature for
    /// a subordinate or consensus element.
    #[allow(clippy::too_many_arguments)]
    fn store_base_feature(
        &self,
        base: &crate::kernel::BaseFeature,
        feature_id: i64,
        parent_id: Option<i64>,
    ) -> Result<()> {
        let mut query = self.prepare(
            "INSERT INTO FEAT_BaseFeature VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10)",
        )?;
        let primary = match base.primary_id {
            Some(molecule) => int(self.molecule_key(molecule)?),
            None => Value::Null,
        };
        insert(
            &mut query,
            &[
                int(feature_id),
                Value::Real(base.rt),
                Value::Real(base.mz),
                Value::Real(f64::from(base.intensity)),
                int(base.charge),
                Value::Real(f64::from(base.width)),
                Value::Real(f64::from(base.quality)),
                // The source binds `int64_t(getUniqueId())`: the bit pattern.
                Value::Integer(base.unique_id as i64),
                primary,
                parent_id.map_or(Value::Null, Value::Integer),
            ],
        )?;
        if !base.id_matches.is_empty() {
            let mut query = self.prepare("INSERT INTO FEAT_ObservationMatch VALUES (?1, ?2)")?;
            for &item in &base.id_matches {
                insert(
                    &mut query,
                    &[int(feature_id), int(key(&self.keys.matches, item)?)],
                )?;
            }
        }
        self.store_meta_info(&base.metadata, "FEAT_BaseFeature", feature_id)
    }

    /// Source `storeFeatureAndSubordinates_`: depth-first, each feature's key
    /// one more than the last one written.
    fn store_feature_and_subordinates(
        &self,
        feature: &Feature,
        next_id: &mut i64,
        parent_id: Option<i64>,
        depth: usize,
    ) -> Result<()> {
        if depth > super::MAX_SUBORDINATE_DEPTH {
            return Err(Error::InvalidValue(
                "OMS store: subordinate features nest too deeply".into(),
            ));
        }
        let feature_id = *next_id;
        self.store_base_feature(&feature.base, feature_id, parent_id)?;
        let mut query = self.prepare("INSERT INTO FEAT_Feature VALUES (?1, ?2, ?3)")?;
        insert(
            &mut query,
            &[
                int(feature_id),
                Value::Real(f64::from(feature.quality_rt)),
                Value::Real(f64::from(feature.quality_mz)),
            ],
        )?;
        if !feature.convex_hulls.is_empty() {
            let mut query =
                self.prepare("INSERT INTO FEAT_ConvexHull VALUES (?1, ?2, ?3, ?4, ?5)")?;
            for (hull_index, hull) in feature.convex_hulls.iter().enumerate() {
                for (point_index, point) in hull.hull_points().iter().enumerate() {
                    insert(
                        &mut query,
                        &[
                            int(feature_id),
                            Value::Integer(hull_index as i64),
                            Value::Integer(point_index as i64),
                            Value::Real(point.rt),
                            Value::Real(point.mz),
                        ],
                    )?;
                }
            }
        }
        *next_id += 1;
        for sub in &feature.subordinates {
            self.store_feature_and_subordinates(sub, next_id, Some(feature_id), depth + 1)?;
        }
        Ok(())
    }

    fn store_features(&mut self, features: &FeatureMap) -> Result<()> {
        if features.features.is_empty() {
            return Ok(());
        }
        let any_metainfo = any_feature(&features.features, &|f| !f.base.metadata.is_empty())?;
        let any_idmatches = any_feature(&features.features, &|f| !f.base.id_matches.is_empty())?;
        self.create_table_base_feature(any_metainfo, any_idmatches)?;
        self.create_table(
            "FEAT_Feature",
            "feature_id INTEGER NOT NULL, rt_quality REAL, mz_quality REAL, \
             FOREIGN KEY (feature_id) REFERENCES FEAT_BaseFeature (id)",
        )?;
        if any_feature(&features.features, &|f| !f.convex_hulls.is_empty())? {
            self.create_table(
                "FEAT_ConvexHull",
                "feature_id INTEGER NOT NULL, \
                 hull_index INTEGER NOT NULL CHECK (hull_index >= 0), \
                 point_index INTEGER NOT NULL CHECK (point_index >= 0), \
                 point_x REAL, point_y REAL, \
                 FOREIGN KEY (feature_id) REFERENCES FEAT_BaseFeature (id)",
            )?;
        }
        let mut next_id = 0i64;
        for feature in &features.features {
            self.store_feature_and_subordinates(feature, &mut next_id, None, 0)?;
        }
        Ok(())
    }

    /// Source `storeMapMetaData_`.
    #[allow(clippy::too_many_arguments)]
    fn store_map_meta_data(
        &mut self,
        unique_id: u64,
        identifier: &str,
        file_path: &str,
        file_type: crate::format::file_types::FileType,
        metadata: &MetaInfo,
        experiment_type: &str,
    ) -> Result<()> {
        self.create_table(
            "FEAT_MapMetaData",
            "unique_id INTEGER PRIMARY KEY, identifier TEXT, file_path TEXT, file_type TEXT, \
             experiment_type TEXT",
        )?;
        let mut query = self.prepare("INSERT INTO FEAT_MapMetaData VALUES (?1, ?2, ?3, ?4, ?5)")?;
        insert(
            &mut query,
            &[
                Value::Integer(unique_id as i64),
                text(identifier),
                text(file_path),
                text(file_type.name()),
                if experiment_type.is_empty() {
                    Value::Null
                } else {
                    text(experiment_type)
                },
            ],
        )?;
        if !metadata.is_empty() {
            self.create_table_meta_info("FEAT_MapMetaData", "unique_id")?;
            self.store_meta_info(metadata, "FEAT_MapMetaData", unique_id as i64)?;
        }
        Ok(())
    }

    /// Source `storeDataProcessing_`.
    fn store_data_processing(&mut self, processing: &[DataProcessing]) -> Result<()> {
        if processing.is_empty() {
            return Ok(());
        }
        self.create_table(
            "FEAT_DataProcessing",
            "id INTEGER PRIMARY KEY NOT NULL, software_name TEXT, software_version TEXT, \
             processing_actions TEXT, completion_time TEXT",
        )?;
        let mut items = Vec::with_capacity(processing.len());
        {
            let mut query =
                self.prepare("INSERT INTO FEAT_DataProcessing VALUES (?1, ?2, ?3, ?4, ?5)")?;
            for (index, proc) in processing.iter().enumerate() {
                let id = index as i64 + 1;
                let actions: Vec<&str> = proc.actions.iter().map(|action| action.name()).collect();
                let completion = match &proc.completion_time {
                    Some(time) => time.get(),
                    None => NULL_DATE_TIME.to_owned(),
                };
                insert(
                    &mut query,
                    &[
                        int(id),
                        text(proc.software.name.as_str()),
                        text(proc.software.version.as_str()),
                        text(actions.join(",")),
                        text(completion),
                    ],
                )?;
                items.push((id, &proc.metadata));
            }
        }
        self.store_meta_infos(items, "FEAT_DataProcessing")
    }

    /// Source `OMSFileStore::store(const FeatureMap&)` without the transaction.
    pub(super) fn store_feature_map(
        &mut self,
        features: &FeatureMap,
        data: Option<&IdentificationData>,
    ) -> Result<()> {
        match data {
            Some(data) if !data.is_empty() => self.store_identifications(data, &mut || Ok(()))?,
            _ => self.store_version_and_date()?,
        }
        self.store_map_meta_data(
            features.unique_id,
            &features.identifier,
            &features.loaded_file_path,
            features.loaded_file_type,
            &features.metadata,
            "",
        )?;
        self.store_data_processing(&features.data_processing)?;
        self.store_features(features)
    }

    /// Source `storeConsensusColumnHeaders_`.
    fn store_consensus_column_headers(&mut self, consensus: &ConsensusMap) -> Result<()> {
        if consensus.column_headers.is_empty() {
            return Ok(());
        }
        self.create_table(
            "FEAT_ConsensusColumnHeader",
            "id INTEGER PRIMARY KEY NOT NULL, filename TEXT, label TEXT, size INTEGER, \
             unique_id INTEGER",
        )?;
        if consensus
            .column_headers
            .values()
            .any(|header| !header.metadata.is_empty())
        {
            self.create_table_meta_info("FEAT_ConsensusColumnHeader", "id")?;
        }
        let mut query =
            self.prepare("INSERT INTO FEAT_ConsensusColumnHeader VALUES (?1, ?2, ?3, ?4, ?5)")?;
        for (&index, header) in &consensus.column_headers {
            let id = index as i64;
            insert(
                &mut query,
                &[
                    int(id),
                    text(header.filename.as_str()),
                    text(header.label.as_str()),
                    Value::Integer(header.size as i64),
                    Value::Integer(header.unique_id as i64),
                ],
            )?;
            self.store_meta_info(&header.metadata, "FEAT_ConsensusColumnHeader", id)?;
        }
        Ok(())
    }

    /// Source `storeConsensusFeatures_`.
    ///
    /// The source writes each consensus feature, then each of its handles
    /// as a subordinate row (keys counting on from the feature's), then the
    /// ratios under the key of the *last handle* (the ratio loop runs after
    /// the handle loop has advanced the key). The port writes ratios under the
    /// consensus feature's own key unless
    /// [`OmsStoreOptions::consensus_ratios_under_last_handle`] asks for the
    /// source's keys.
    fn store_consensus_features(&mut self, consensus: &ConsensusMap) -> Result<()> {
        if consensus.features.is_empty() {
            return Ok(());
        }
        let any_metainfo = consensus
            .features
            .iter()
            .any(|f| !f.base.metadata.is_empty());
        let any_idmatches = consensus
            .features
            .iter()
            .any(|f| !f.base.id_matches.is_empty());
        self.create_table_base_feature(any_metainfo, any_idmatches)?;
        self.create_table(
            "FEAT_FeatureHandle",
            "feature_id INTEGER NOT NULL, map_index INTEGER NOT NULL, \
             FOREIGN KEY (feature_id) REFERENCES FEAT_BaseFeature (id)",
        )?;
        let any_ratios = consensus.features.iter().any(|f| !f.ratios.is_empty());
        if any_ratios {
            self.create_table(
                "FEAT_ConsensusRatio",
                "feature_id INTEGER NOT NULL, \
                 ratio_index INTEGER NOT NULL CHECK (ratio_index >= 0), ratio_value REAL, \
                 denominator_ref TEXT, numerator_ref TEXT, description TEXT, \
                 FOREIGN KEY (feature_id) REFERENCES FEAT_BaseFeature (id)",
            )?;
        }
        let mut query_handle = self.prepare("INSERT INTO FEAT_FeatureHandle VALUES (?1, ?2)")?;
        let mut feature_id = 0i64;
        for feature in &consensus.features {
            self.store_base_feature(&feature.base, feature_id, None)?;
            let parent_id = feature_id;
            for handle in feature.handles() {
                feature_id += 1;
                let mut base = handle_base(handle);
                if self.options.drop_handle_unique_ids {
                    base.unique_id = 0;
                }
                self.store_base_feature(&base, feature_id, Some(parent_id))?;
                insert(
                    &mut query_handle,
                    &[int(feature_id), Value::Integer(handle.map_index as i64)],
                )?;
            }
            if !feature.ratios.is_empty() {
                let ratio_key = if self.options.consensus_ratios_under_last_handle {
                    feature_id
                } else {
                    parent_id
                };
                let mut query_ratio = self
                    .prepare("INSERT INTO FEAT_ConsensusRatio VALUES (?1, ?2, ?3, ?4, ?5, ?6)")?;
                for (index, ratio) in feature.ratios.iter().enumerate() {
                    insert(
                        &mut query_ratio,
                        &[
                            int(ratio_key),
                            size_value(index, "ratio index")?,
                            Value::Real(ratio.ratio_value),
                            text(ratio.denominator_ref.as_str()),
                            text(ratio.numerator_ref.as_str()),
                            text(ratio.description.join(",")),
                        ],
                    )?;
                }
            }
            feature_id += 1;
        }
        Ok(())
    }

    /// Source `OMSFileStore::store(const ConsensusMap&)` without the transaction.
    pub(super) fn store_consensus_map(
        &mut self,
        consensus: &ConsensusMap,
        data: Option<&IdentificationData>,
    ) -> Result<()> {
        match data {
            Some(data) if !data.is_empty() => self.store_identifications(data, &mut || Ok(()))?,
            _ => self.store_version_and_date()?,
        }
        self.store_map_meta_data(
            consensus.unique_id,
            &consensus.identifier,
            &consensus.loaded_file_path,
            consensus.loaded_file_type,
            &consensus.metadata,
            &consensus.experiment_type,
        )?;
        self.store_consensus_column_headers(consensus)?;
        self.store_data_processing(&consensus.data_processing)?;
        self.store_consensus_features(consensus)
    }
}

/// Source `anyFeaturePredicate_`, recursing into subordinates with a depth
/// ceiling.
fn any_feature(features: &[Feature], predicate: &dyn Fn(&Feature) -> bool) -> Result<bool> {
    fn walk(
        features: &[Feature],
        predicate: &dyn Fn(&Feature) -> bool,
        depth: usize,
    ) -> Result<bool> {
        if depth > super::MAX_SUBORDINATE_DEPTH {
            return Err(Error::InvalidValue(
                "OMS store: subordinate features nest too deeply".into(),
            ));
        }
        for feature in features {
            if predicate(feature) || walk(&feature.subordinates, predicate, depth + 1)? {
                return Ok(true);
            }
        }
        Ok(false)
    }
    walk(features, predicate, 0)
}

/// Source `BaseFeature(const FeatureHandle&)`: position, intensity, charge,
/// width and unique ID of the handle, quality zero, no identifications.
pub(super) fn handle_base(handle: &crate::kernel::FeatureHandle) -> crate::kernel::BaseFeature {
    let mut base = crate::kernel::BaseFeature::new(handle.rt, handle.mz, handle.intensity);
    base.charge = handle.charge;
    base.width = handle.width;
    base.unique_id = handle.unique_id;
    base
}
