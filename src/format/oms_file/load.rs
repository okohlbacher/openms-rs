// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Hendrik Weisser, OpenMS Rust contributors $

//! Reading OMS files: the port of `FORMAT/OMSFileLoad.h` (`Internal::OMSFileLoad`).
//!
//! Every schema version the source reads is read here: version 1 files lack
//! `enzyme_term_specificity`, versions before 3 order parent group sets by
//! `grouping_order` and data processing by `position`, versions before 4 keep
//! meta values in a shared `DataValue` table, and versions before 5 keep
//! features in `FEAT_Feature` with an `overall_quality` column.

use super::store::NULL_DATE_TIME;
use super::{MAX_SUBORDINATE_DEPTH, MAX_TABLE_ROWS, OmsLoadOptions, sql_error};
use crate::chemistry::{
    AASequence, AdductInfo, DigestionSpecificity, EmpiricalFormula, NASequence, ProteaseDB, RNaseDB,
};
use crate::data_structures::list;
use crate::format::file_types::FileType;
use crate::identification::PeakAnnotation;
use crate::identification::graph::{
    AdductId, AppliedProcessingStep, DBSearchParam, GraphEnzyme, IdentificationData,
    IdentifiedCompound, IdentifiedMolecule, IdentifiedOligo, IdentifiedPeptide, InputFile,
    InputFileId, MassType, MoleculeType, Observation, ObservationId, ObservationMatch,
    ObservationMatchId, ParentGroup, ParentGroupSet, ParentId, ParentMatch, ParentMatches,
    ParentSequence, ProcessingSoftware, ProcessingSoftwareId, ProcessingStep, ProcessingStepId,
    ScoreType, ScoreTypeId, ScoredProcessingResult, SearchParamId,
};
use crate::kernel::features::Ratio;
use crate::kernel::{
    BaseFeature, ColumnHeader, ConsensusFeature, ConsensusMap, ConvexHull2D, Feature,
    FeatureHandle, FeatureMap, Point2D,
};
use crate::metadata::{
    CVTerm, CompletionTime, DataProcessing, MetaInfo, MetaValue, MetaValueData, ProcessingAction,
    Software,
};
use crate::{Error, Result};
use rusqlite::types::ValueRef;
use rusqlite::{Connection, Row, Statement, params_from_iter};
use std::collections::BTreeMap;

fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

/// A cell as the source's `Column::getString()` reads it: `NULL` is empty.
fn cell_text(value: ValueRef<'_>) -> Result<String> {
    match value {
        ValueRef::Null => Ok(String::new()),
        ValueRef::Integer(number) => Ok(number.to_string()),
        ValueRef::Real(number) => Ok(number.to_string()),
        ValueRef::Text(bytes) | ValueRef::Blob(bytes) => std::str::from_utf8(bytes)
            .map(str::to_owned)
            .map_err(|_| invalid("OMS file holds a text cell that is not UTF-8")),
    }
}

/// A cell as `Column::getInt64()` reads it: `NULL` is zero. Text must be an
/// integer here; SQLite's own conversion would silently read a prefix.
fn cell_int(value: ValueRef<'_>) -> Result<i64> {
    match value {
        ValueRef::Null => Ok(0),
        ValueRef::Integer(number) => Ok(number),
        ValueRef::Real(number) if number.is_finite() && number.fract() == 0.0 => Ok(number as i64),
        ValueRef::Text(bytes) => std::str::from_utf8(bytes)
            .ok()
            .and_then(|text| text.trim().parse().ok())
            .ok_or_else(|| invalid("OMS file holds a non-integer cell in an integer column")),
        _ => Err(invalid(
            "OMS file holds a non-integer cell in an integer column",
        )),
    }
}

/// A cell as `Column::getDouble()` reads it: `NULL` is zero.
fn cell_real(value: ValueRef<'_>) -> Result<f64> {
    match value {
        ValueRef::Null => Ok(0.0),
        ValueRef::Integer(number) => Ok(number as f64),
        ValueRef::Real(number) => Ok(number),
        ValueRef::Text(bytes) => std::str::from_utf8(bytes)
            .ok()
            .and_then(|text| text.trim().parse().ok())
            .ok_or_else(|| invalid("OMS file holds a non-numeric cell in a numeric column")),
        ValueRef::Blob(_) => Err(invalid("OMS file holds a blob in a numeric column")),
    }
}

fn column<'r>(row: &'r Row<'_>, name: &str) -> Result<ValueRef<'r>> {
    row.get_ref(name)
        .map_err(|error| sql_error("error reading OMS column", error))
}

fn text(row: &Row<'_>, name: &str) -> Result<String> {
    cell_text(column(row, name)?)
}

fn int(row: &Row<'_>, name: &str) -> Result<i64> {
    cell_int(column(row, name)?)
}

fn int32(row: &Row<'_>, name: &str) -> Result<i32> {
    i32::try_from(int(row, name)?)
        .map_err(|_| invalid(format!("OMS column {name} exceeds 32 bits")))
}

fn real(row: &Row<'_>, name: &str) -> Result<f64> {
    cell_real(column(row, name)?)
}

fn is_null(row: &Row<'_>, name: &str) -> Result<bool> {
    Ok(matches!(column(row, name)?, ValueRef::Null))
}

fn optional_int(row: &Row<'_>, name: &str) -> Result<Option<i64>> {
    if is_null(row, name)? {
        Ok(None)
    } else {
        int(row, name).map(Some)
    }
}

fn size(row: &Row<'_>, name: &str) -> Result<usize> {
    // The source reads these with `getUInt()`; a negative value is refused.
    usize::try_from(int(row, name)?).map_err(|_| invalid(format!("OMS column {name} is negative")))
}

/// Run `statement` with `params`, calling `each` per row, under the row ceiling.
fn for_rows(
    statement: &mut Statement<'_>,
    params: &[rusqlite::types::Value],
    mut each: impl FnMut(&Row<'_>) -> Result<()>,
) -> Result<()> {
    let mut rows = statement
        .query(params_from_iter(params.iter()))
        .map_err(|error| sql_error("error querying database", error))?;
    let mut count = 0usize;
    while let Some(row) = rows
        .next()
        .map_err(|error| sql_error("error querying database", error))?
    {
        count += 1;
        if count > MAX_TABLE_ROWS {
            return Err(invalid("OMS table exceeds the row limit"));
        }
        each(row)?;
    }
    Ok(())
}

fn lookup<K: Copy>(map: &BTreeMap<i64, K>, key: i64, what: &str) -> Result<K> {
    // The source's `unordered_map::operator[]` default-constructs an invalid
    // reference for a dangling key; the port refuses the file.
    map.get(&key)
        .copied()
        .ok_or_else(|| invalid(format!("OMS file references a missing {what} (key {key})")))
}

/// Strip the brackets `DataValue::toString()` puts around a list, as the
/// source's `substr(value, 1, size - 2)` does.
fn list_body(value: &str) -> &str {
    if value.len() >= 2 {
        value.get(1..value.len() - 1).unwrap_or("")
    } else {
        ""
    }
}

/// Source `makeDataValue_`: the stored text converted by its type, with lists
/// split at `,` after the brackets are removed (a numeric item is trimmed).
///
/// `DataValue::toString()` joins list items with `", "`. The source keeps the
/// space as the first character of every string item after the first
/// (`keep_string_list_spaces`); by default the port removes that one space.
pub(super) fn make_meta_value(
    type_id: i64,
    value: &str,
    keep_string_list_spaces: bool,
) -> Result<MetaValue> {
    let data = match type_id {
        1 => MetaValueData::String(value.to_owned()),
        2 => MetaValueData::Integer(i64::from(<i32 as list::ListParse>::from_list_item(value)?)),
        3 => MetaValueData::Float(<f64 as list::ListParse>::from_list_item(value)?),
        4 => {
            let mut items = list::create::<String>(list_body(value), b',')?;
            if !keep_string_list_spaces {
                for item in items.iter_mut().skip(1) {
                    if item.starts_with(' ') {
                        item.remove(0);
                    }
                }
            }
            MetaValueData::StringList(items)
        }
        5 => MetaValueData::IntegerList(
            list::create::<i32>(list_body(value), b',')?
                .into_iter()
                .map(i64::from)
                .collect(),
        ),
        6 => MetaValueData::FloatList(list::create::<f64>(list_body(value), b',')?),
        _ => MetaValueData::Empty,
    };
    MetaValue::new(data)
}

/// The source's `ListUtils::create<String>` of a comma-joined column; empty
/// text is an empty list.
fn string_list(value: &str) -> Result<Vec<String>> {
    list::create::<String>(value, b',')
}

/// Reader over one open database (source `OMSFileLoad`).
pub(super) struct Load<'c> {
    db: &'c Connection,
    version: i64,
    options: OmsLoadOptions,
    score_types: BTreeMap<i64, ScoreTypeId>,
    input_files: BTreeMap<i64, InputFileId>,
    softwares: BTreeMap<i64, ProcessingSoftwareId>,
    steps: BTreeMap<i64, ProcessingStepId>,
    search_params: BTreeMap<i64, SearchParamId>,
    observations: BTreeMap<i64, ObservationId>,
    parents: BTreeMap<i64, ParentId>,
    molecules: BTreeMap<i64, IdentifiedMolecule>,
    matches: BTreeMap<i64, ObservationMatchId>,
    adducts: BTreeMap<i64, AdductId>,
}

impl<'c> Load<'c> {
    pub(super) fn new(db: &'c Connection, version: i64, options: OmsLoadOptions) -> Self {
        Self {
            db,
            version,
            options,
            score_types: BTreeMap::new(),
            input_files: BTreeMap::new(),
            softwares: BTreeMap::new(),
            steps: BTreeMap::new(),
            search_params: BTreeMap::new(),
            observations: BTreeMap::new(),
            parents: BTreeMap::new(),
            molecules: BTreeMap::new(),
            matches: BTreeMap::new(),
            adducts: BTreeMap::new(),
        }
    }

    fn prepare(&self, sql: &str) -> Result<Statement<'c>> {
        self.db
            .prepare(sql)
            .map_err(|error| sql_error("error preparing SQL query", error))
    }

    fn table_exists(&self, name: &str) -> Result<bool> {
        self.db
            .query_row(
                "SELECT count(*) FROM sqlite_master WHERE type='table' AND name = ?1",
                [name],
                |row| row.get::<_, i64>(0),
            )
            .map(|count| count > 0)
            .map_err(|error| sql_error("error querying database", error))
    }

    /// Source `prepareQueryMetaInfo_`: the meta-info query of a table, if the
    /// table has metadata.
    fn meta_query(&self, parent_table: &str) -> Result<Option<Statement<'c>>> {
        let table = format!("{parent_table}_MetaInfo");
        if !self.table_exists(&table)? {
            return Ok(None);
        }
        let sql = if self.version < 4 {
            format!(
                "SELECT * FROM {table} AS MI JOIN DataValue AS DV ON MI.data_value_id = DV.id \
                 WHERE MI.parent_id = ?1"
            )
        } else {
            format!("SELECT * FROM {table} AS MI WHERE MI.parent_id = ?1")
        };
        self.prepare(&sql).map(Some)
    }

    /// Source `handleQueryMetaInfo_`.
    fn read_meta(
        keep_spaces: bool,
        query: &mut Option<Statement<'c>>,
        info: &mut MetaInfo,
        parent_id: i64,
    ) -> Result<()> {
        let Some(query) = query.as_mut() else {
            return Ok(());
        };
        for_rows(query, &[parent_id.into()], |row| {
            let value =
                make_meta_value(int(row, "data_type_id")?, &text(row, "value")?, keep_spaces)?;
            info.insert(text(row, "name")?, value);
            Ok(())
        })
    }

    /// Source `prepareQueryAppliedProcessingStep_`.
    fn step_query(&self, parent_table: &str) -> Result<Option<Statement<'c>>> {
        let table = format!("{parent_table}_AppliedProcessingStep");
        if !self.table_exists(&table)? {
            return Ok(None);
        }
        self.prepare(&format!(
            "SELECT * FROM {table} WHERE parent_id = ?1 ORDER BY processing_step_order ASC"
        ))
        .map(Some)
    }

    /// Source `handleQueryAppliedProcessingStep_`: rows of the same step merge
    /// through `addProcessingStep`.
    fn read_steps(
        &self,
        query: &mut Option<Statement<'c>>,
        result: &mut ScoredProcessingResult,
        parent_id: i64,
    ) -> Result<()> {
        let Some(query) = query.as_mut() else {
            return Ok(());
        };
        for_rows(query, &[parent_id.into()], |row| {
            let mut step = AppliedProcessingStep::default();
            if let Some(id) = optional_int(row, "processing_step_id")? {
                step.processing_step = Some(lookup(&self.steps, id, "processing step")?);
            }
            if let Some(id) = optional_int(row, "score_type_id")? {
                step.scores.insert(
                    lookup(&self.score_types, id, "score type")?,
                    real(row, "score")?,
                );
            }
            result.add_processing_step(step)
        })
    }

    fn load_input_files(&mut self, data: &mut IdentificationData) -> Result<()> {
        if !self.table_exists("ID_InputFile")? {
            return Ok(());
        }
        let mut query = self.prepare("SELECT * FROM ID_InputFile")?;
        let mut loaded = Vec::new();
        for_rows(&mut query, &[], |row| {
            let mut input = InputFile::new(text(row, "name")?);
            input.experimental_design_id = text(row, "experimental_design_id")?;
            input
                .primary_files
                .extend(string_list(&text(row, "primary_files")?)?);
            loaded.push((int(row, "id")?, input));
            Ok(())
        })?;
        for (key, input) in loaded {
            let id = data.register_input_file(input)?;
            self.input_files.insert(key, id);
        }
        Ok(())
    }

    fn load_score_types(&mut self, data: &mut IdentificationData) -> Result<()> {
        if !self.table_exists("ID_ScoreType")? {
            return Ok(());
        }
        if !self.table_exists("CVTerm")? {
            // Source: Exception::MissingInformation.
            return Err(invalid("required database table 'CVTerm' not found"));
        }
        let mut query = self.prepare(
            "SELECT S.*, C.accession, C.name, C.cv_identifier_ref \
             FROM ID_ScoreType AS S JOIN CVTerm AS C ON S.cv_term_id = C.id",
        )?;
        let mut loaded = Vec::new();
        for_rows(&mut query, &[], |row| {
            let term = CVTerm::new(
                text(row, "accession")?,
                text(row, "name")?,
                text(row, "cv_identifier_ref")?,
            );
            let score_type = ScoreType {
                cv_term: term,
                higher_better: int(row, "higher_better")? != 0,
                metadata: MetaInfo::new(),
            };
            loaded.push((int(row, "id")?, score_type));
            Ok(())
        })?;
        for (key, score_type) in loaded {
            let id = data.register_score_type(score_type)?;
            self.score_types.insert(key, id);
        }
        Ok(())
    }

    fn load_processing_softwares(&mut self, data: &mut IdentificationData) -> Result<()> {
        if !self.table_exists("ID_ProcessingSoftware")? {
            return Ok(());
        }
        let mut subquery = if self.table_exists("ID_ProcessingSoftware_AssignedScore")? {
            Some(self.prepare(
                "SELECT score_type_id FROM ID_ProcessingSoftware_AssignedScore \
                 WHERE software_id = ?1 ORDER BY score_type_order ASC",
            )?)
        } else {
            None
        };
        let mut query = self.prepare("SELECT * FROM ID_ProcessingSoftware")?;
        let mut loaded = Vec::new();
        for_rows(&mut query, &[], |row| {
            let key = int(row, "id")?;
            let mut software = ProcessingSoftware::new(text(row, "name")?, text(row, "version")?);
            if let Some(subquery) = subquery.as_mut() {
                for_rows(subquery, &[key.into()], |sub| {
                    let score = cell_int(
                        sub.get_ref(0)
                            .map_err(|e| sql_error("error reading OMS column", e))?,
                    )?;
                    software
                        .assigned_scores
                        .push(lookup(&self.score_types, score, "score type")?);
                    Ok(())
                })?;
            }
            loaded.push((key, software));
            Ok(())
        })?;
        for (key, software) in loaded {
            let id = data.register_processing_software(software)?;
            self.softwares.insert(key, id);
        }
        Ok(())
    }

    fn load_db_search_params(&mut self, data: &mut IdentificationData) -> Result<()> {
        if !self.table_exists("ID_DBSearchParam")? {
            return Ok(());
        }
        let mut query = self.prepare("SELECT * FROM ID_DBSearchParam")?;
        let version = self.version;
        let mut loaded = Vec::new();
        for_rows(&mut query, &[], |row| {
            let mut param = DBSearchParam {
                molecule_type: molecule_type(int(row, "molecule_type_id")?)?,
                mass_type: if int(row, "mass_type_average")? != 0 {
                    MassType::Average
                } else {
                    MassType::Monoisotopic
                },
                database: text(row, "database")?,
                database_version: text(row, "database_version")?,
                taxonomy: text(row, "taxonomy")?,
                ..DBSearchParam::default()
            };
            param
                .charges
                .extend(list::create::<i32>(&text(row, "charges")?, b',')?);
            param
                .fixed_mods
                .extend(string_list(&text(row, "fixed_mods")?)?);
            param
                .variable_mods
                .extend(string_list(&text(row, "variable_mods")?)?);
            param.precursor_mass_tolerance = real(row, "precursor_mass_tolerance")?;
            param.fragment_mass_tolerance = real(row, "fragment_mass_tolerance")?;
            param.precursor_tolerance_ppm = int(row, "precursor_tolerance_ppm")? != 0;
            param.fragment_tolerance_ppm = int(row, "fragment_tolerance_ppm")? != 0;
            let enzyme = text(row, "digestion_enzyme")?;
            if !enzyme.is_empty() {
                param.digestion_enzyme = match param.molecule_type {
                    MoleculeType::Protein => Some(GraphEnzyme::Protein(
                        *ProteaseDB::global().get_enzyme(&enzyme)?,
                    )),
                    MoleculeType::RNA => {
                        Some(GraphEnzyme::RNA(RNaseDB::global().get_enzyme(&enzyme)?))
                    }
                    MoleculeType::Compound => None,
                };
            }
            if version > 1 {
                // Source `EnzymaticDigestion::getSpecificityByName`: an unknown
                // name is `SPEC_UNKNOWN`.
                param.enzyme_term_specificity = text(row, "enzyme_term_specificity")?
                    .parse::<DigestionSpecificity>()
                    .ok();
            }
            param.missed_cleavages = size(row, "missed_cleavages")?;
            param.min_length = size(row, "min_length")?;
            param.max_length = size(row, "max_length")?;
            loaded.push((int(row, "id")?, param));
            Ok(())
        })?;
        for (key, param) in loaded {
            let id = data.register_db_search_param(param)?;
            self.search_params.insert(key, id);
        }
        Ok(())
    }

    fn load_processing_steps(&mut self, data: &mut IdentificationData) -> Result<()> {
        if !self.table_exists("ID_ProcessingStep")? {
            return Ok(());
        }
        let mut subquery_file = if self.table_exists("ID_ProcessingStep_InputFile")? {
            Some(self.prepare(
                "SELECT input_file_id FROM ID_ProcessingStep_InputFile WHERE processing_step_id = ?1",
            )?)
        } else {
            None
        };
        let mut subquery_info = self.meta_query("ID_ProcessingStep")?;
        let mut query = self.prepare("SELECT * FROM ID_ProcessingStep")?;
        let mut loaded = Vec::new();
        for_rows(&mut query, &[], |row| {
            let key = int(row, "id")?;
            let mut step = ProcessingStep::new(lookup(
                &self.softwares,
                int(row, "software_id")?,
                "software",
            )?);
            step.date_time = parse_date_time(&text(row, "date_time")?)?;
            if let Some(subquery) = subquery_file.as_mut() {
                for_rows(subquery, &[key.into()], |sub| {
                    let file = cell_int(
                        sub.get_ref(0)
                            .map_err(|e| sql_error("error reading OMS column", e))?,
                    )?;
                    step.input_files
                        .push(lookup(&self.input_files, file, "input file")?);
                    Ok(())
                })?;
            }
            Self::read_meta(
                self.options.keep_string_list_spaces,
                &mut subquery_info,
                &mut step.metadata,
                key,
            )?;
            let search = match optional_int(row, "search_param_id")? {
                Some(id) => Some(lookup(&self.search_params, id, "search parameter set")?),
                None => None,
            };
            loaded.push((key, step, search));
            Ok(())
        })?;
        for (key, step, search) in loaded {
            let id = data.register_processing_step(step, search)?;
            self.steps.insert(key, id);
        }
        Ok(())
    }

    fn load_observations(&mut self, data: &mut IdentificationData) -> Result<()> {
        if !self.table_exists("ID_Observation")? {
            return Ok(());
        }
        let mut subquery_info = self.meta_query("ID_Observation")?;
        let mut query = self.prepare("SELECT * FROM ID_Observation")?;
        let mut loaded = Vec::new();
        for_rows(&mut query, &[], |row| {
            let input = lookup(&self.input_files, int(row, "input_file_id")?, "input file")?;
            let mut observation = Observation::new(text(row, "data_id")?, input);
            if !is_null(row, "rt")? {
                observation.rt = Some(real(row, "rt")?);
            }
            if !is_null(row, "mz")? {
                observation.mz = Some(real(row, "mz")?);
            }
            let key = int(row, "id")?;
            Self::read_meta(
                self.options.keep_string_list_spaces,
                &mut subquery_info,
                &mut observation.metadata,
                key,
            )?;
            loaded.push((key, observation));
            Ok(())
        })?;
        for (key, observation) in loaded {
            let id = data.register_observation(observation)?;
            self.observations.insert(key, id);
        }
        Ok(())
    }

    fn load_parent_sequences(&mut self, data: &mut IdentificationData) -> Result<()> {
        if !self.table_exists("ID_ParentSequence")? {
            return Ok(());
        }
        let mut subquery_info = self.meta_query("ID_ParentSequence")?;
        let mut subquery_step = self.step_query("ID_ParentSequence")?;
        let mut query = self.prepare("SELECT * FROM ID_ParentSequence")?;
        let mut loaded = Vec::new();
        for_rows(&mut query, &[], |row| {
            let mut parent = ParentSequence::new(text(row, "accession")?);
            parent.molecule_type = molecule_type(int(row, "molecule_type_id")?)?;
            parent.sequence = text(row, "sequence")?;
            parent.description = text(row, "description")?;
            parent.coverage = real(row, "coverage")?;
            parent.is_decoy = int(row, "is_decoy")? != 0;
            let key = int(row, "id")?;
            Self::read_meta(
                self.options.keep_string_list_spaces,
                &mut subquery_info,
                &mut parent.result.metadata,
                key,
            )?;
            self.read_steps(&mut subquery_step, &mut parent.result, key)?;
            loaded.push((key, parent));
            Ok(())
        })?;
        for (key, parent) in loaded {
            let id = data.register_parent_sequence(parent)?;
            self.parents.insert(key, id);
        }
        Ok(())
    }

    fn load_parent_group_sets(&mut self, data: &mut IdentificationData) -> Result<()> {
        if !self.table_exists("ID_ParentGroupSet")? {
            return Ok(());
        }
        let order_by = if self.version > 2 {
            "id"
        } else {
            "grouping_order"
        };
        let mut subquery_info = self.meta_query("ID_ParentGroupSet")?;
        let mut subquery_step = self.step_query("ID_ParentGroupSet")?;
        let mut subquery_group =
            self.prepare("SELECT * FROM ID_ParentGroup WHERE grouping_id = ?1")?;
        let mut subquery_parent = self
            .prepare("SELECT parent_id FROM ID_ParentGroup_ParentSequence WHERE group_id = ?1")?;
        let mut query = self.prepare(&format!(
            "SELECT * FROM ID_ParentGroupSet ORDER BY {order_by} ASC"
        ))?;
        let mut loaded = Vec::new();
        for_rows(&mut query, &[], |row| {
            let mut grouping = ParentGroupSet::new(text(row, "label")?);
            let key = int(row, "id")?;
            Self::read_meta(
                self.options.keep_string_list_spaces,
                &mut subquery_info,
                &mut grouping.result.metadata,
                key,
            )?;
            self.read_steps(&mut subquery_step, &mut grouping.result, key)?;
            let mut groups: BTreeMap<i64, ParentGroup> = BTreeMap::new();
            for_rows(&mut subquery_group, &[key.into()], |group_row| {
                let group = groups.entry(int(group_row, "id")?).or_default();
                if let Some(score_type) = optional_int(group_row, "score_type_id")? {
                    group.scores.insert(
                        lookup(&self.score_types, score_type, "score type")?,
                        real(group_row, "score")?,
                    );
                }
                Ok(())
            })?;
            for (group_key, mut group) in groups {
                for_rows(&mut subquery_parent, &[group_key.into()], |parent_row| {
                    let parent = cell_int(
                        parent_row
                            .get_ref(0)
                            .map_err(|e| sql_error("error reading OMS column", e))?,
                    )?;
                    group
                        .parent_refs
                        .insert(lookup(&self.parents, parent, "parent sequence")?);
                    Ok(())
                })?;
                grouping.groups.push(group);
            }
            loaded.push(grouping);
            Ok(())
        })?;
        for grouping in loaded {
            data.register_parent_group_set(grouping)?;
        }
        Ok(())
    }

    fn load_identified_compounds(&mut self, data: &mut IdentificationData) -> Result<()> {
        if !self.table_exists("ID_IdentifiedCompound")? {
            return Ok(());
        }
        let mut subquery_info = self.meta_query("ID_IdentifiedMolecule")?;
        let mut subquery_step = self.step_query("ID_IdentifiedMolecule")?;
        let mut query = self.prepare(
            "SELECT * FROM ID_IdentifiedMolecule JOIN ID_IdentifiedCompound \
             ON ID_IdentifiedMolecule.id = ID_IdentifiedCompound.molecule_id",
        )?;
        let mut loaded = Vec::new();
        for_rows(&mut query, &[], |row| {
            let mut compound = IdentifiedCompound::new(text(row, "identifier")?);
            let formula = text(row, "formula")?;
            compound.formula = if formula.is_empty() {
                EmpiricalFormula::default()
            } else {
                EmpiricalFormula::parse(&formula)?
            };
            compound.name = text(row, "name")?;
            compound.smile = text(row, "smile")?;
            compound.inchi = text(row, "inchi")?;
            let key = int(row, "id")?;
            Self::read_meta(
                self.options.keep_string_list_spaces,
                &mut subquery_info,
                &mut compound.result.metadata,
                key,
            )?;
            self.read_steps(&mut subquery_step, &mut compound.result, key)?;
            loaded.push((key, compound));
            Ok(())
        })?;
        for (key, compound) in loaded {
            let id = data.register_identified_compound(compound)?;
            self.molecules.insert(key, IdentifiedMolecule::Compound(id));
        }
        Ok(())
    }

    /// Source `handleQueryParentMatch_`.
    fn read_parent_matches(
        &self,
        query: &mut Option<Statement<'c>>,
        matches: &mut ParentMatches,
        molecule_id: i64,
    ) -> Result<()> {
        let Some(query) = query.as_mut() else {
            return Ok(());
        };
        for_rows(query, &[molecule_id.into()], |row| {
            let parent = lookup(&self.parents, int(row, "parent_id")?, "parent sequence")?;
            let start = if is_null(row, "start_pos")? {
                None
            } else {
                Some(size(row, "start_pos")?)
            };
            let end = if is_null(row, "end_pos")? {
                None
            } else {
                Some(size(row, "end_pos")?)
            };
            let mut item = ParentMatch::new(start, end);
            item.left_neighbor = text(row, "left_neighbor")?;
            item.right_neighbor = text(row, "right_neighbor")?;
            matches.entry(parent).or_default().insert(item);
            Ok(())
        })
    }

    fn load_identified_sequences(&mut self, data: &mut IdentificationData) -> Result<()> {
        if !self.table_exists("ID_IdentifiedMolecule")? {
            return Ok(());
        }
        let mut subquery_info = self.meta_query("ID_IdentifiedMolecule")?;
        let mut subquery_step = self.step_query("ID_IdentifiedMolecule")?;
        let mut subquery_parent = if self.table_exists("ID_ParentMatch")? {
            Some(self.prepare("SELECT * FROM ID_ParentMatch WHERE molecule_id = ?1")?)
        } else {
            None
        };
        let mut query =
            self.prepare("SELECT * FROM ID_IdentifiedMolecule WHERE molecule_type_id = ?1")?;
        let mut peptides = Vec::new();
        for_rows(&mut query, &[1i64.into()], |row| {
            let key = int(row, "id")?;
            let mut peptide = IdentifiedPeptide::new(AASequence::parse(&text(row, "identifier")?)?);
            Self::read_meta(
                self.options.keep_string_list_spaces,
                &mut subquery_info,
                &mut peptide.result.metadata,
                key,
            )?;
            self.read_steps(&mut subquery_step, &mut peptide.result, key)?;
            self.read_parent_matches(&mut subquery_parent, &mut peptide.parent_matches, key)?;
            peptides.push((key, peptide));
            Ok(())
        })?;
        for (key, peptide) in peptides {
            let id = data.register_identified_peptide(peptide)?;
            self.molecules.insert(key, IdentifiedMolecule::Peptide(id));
        }
        let mut oligos = Vec::new();
        for_rows(&mut query, &[3i64.into()], |row| {
            let key = int(row, "id")?;
            let mut oligo = IdentifiedOligo::new(NASequence::parse(&text(row, "identifier")?)?);
            Self::read_meta(
                self.options.keep_string_list_spaces,
                &mut subquery_info,
                &mut oligo.result.metadata,
                key,
            )?;
            self.read_steps(&mut subquery_step, &mut oligo.result, key)?;
            self.read_parent_matches(&mut subquery_parent, &mut oligo.parent_matches, key)?;
            oligos.push((key, oligo));
            Ok(())
        })?;
        for (key, oligo) in oligos {
            let id = data.register_identified_oligo(oligo)?;
            self.molecules.insert(key, IdentifiedMolecule::Oligo(id));
        }
        Ok(())
    }

    fn load_adducts(&mut self, data: &mut IdentificationData) -> Result<()> {
        if !self.table_exists("AdductInfo")? {
            return Ok(());
        }
        let mut query = self.prepare("SELECT * FROM AdductInfo")?;
        let mut loaded = Vec::new();
        for_rows(&mut query, &[], |row| {
            let multiplier = u32::try_from(int(row, "mol_multiplier")?)
                .map_err(|_| invalid("OMS adduct multiplier out of range"))?;
            let adduct = AdductInfo::new(
                text(row, "name")?,
                EmpiricalFormula::parse(&text(row, "formula")?)?,
                int32(row, "charge")?,
                multiplier,
            )?;
            loaded.push((int(row, "id")?, adduct));
            Ok(())
        })?;
        for (key, adduct) in loaded {
            let id = data.register_adduct(adduct)?;
            self.adducts.insert(key, id);
        }
        Ok(())
    }

    fn load_observation_matches(&mut self, data: &mut IdentificationData) -> Result<()> {
        if !self.table_exists("ID_ObservationMatch")? {
            return Ok(());
        }
        let mut subquery_info = self.meta_query("ID_ObservationMatch")?;
        let mut subquery_step = self.step_query("ID_ObservationMatch")?;
        let mut subquery_annotation =
            if self.table_exists("ID_ObservationMatch_PeakAnnotation")? {
                Some(self.prepare(
                    "SELECT * FROM ID_ObservationMatch_PeakAnnotation WHERE parent_id = ?1",
                )?)
            } else {
                None
            };
        let mut query = self.prepare("SELECT * FROM ID_ObservationMatch")?;
        let mut loaded = Vec::new();
        for_rows(&mut query, &[], |row| {
            let key = int(row, "id")?;
            let molecule = lookup(
                &self.molecules,
                int(row, "identified_molecule_id")?,
                "identified molecule",
            )?;
            let observation = lookup(
                &self.observations,
                int(row, "observation_id")?,
                "observation",
            )?;
            let mut item = ObservationMatch::new(molecule, observation);
            item.charge = int32(row, "charge")?;
            if let Some(adduct) = optional_int(row, "adduct_id")? {
                item.adduct = Some(lookup(&self.adducts, adduct, "adduct")?);
            }
            Self::read_meta(
                self.options.keep_string_list_spaces,
                &mut subquery_info,
                &mut item.result.metadata,
                key,
            )?;
            self.read_steps(&mut subquery_step, &mut item.result, key)?;
            if let Some(subquery) = subquery_annotation.as_mut() {
                for_rows(subquery, &[key.into()], |ann| {
                    let step = match optional_int(ann, "processing_step_id")? {
                        Some(id) => Some(lookup(&self.steps, id, "processing step")?),
                        None => None,
                    };
                    item.peak_annotations
                        .entry(step)
                        .or_default()
                        .push(PeakAnnotation {
                            annotation: text(ann, "peak_annotation")?,
                            charge: int32(ann, "peak_charge")?,
                            mz: real(ann, "peak_mz")?,
                            intensity: real(ann, "peak_intensity")?,
                        });
                    Ok(())
                })?;
            }
            loaded.push((key, item));
            Ok(())
        })?;
        for (key, item) in loaded {
            let id = data.register_observation_match(item)?;
            self.matches.insert(key, id);
        }
        Ok(())
    }

    /// Source `OMSFileLoad::load(IdentificationData&)`.
    pub(super) fn load_identifications(
        &mut self,
        data: &mut IdentificationData,
        progress: &mut dyn FnMut() -> Result<()>,
    ) -> Result<()> {
        self.load_input_files(data)?;
        progress()?;
        self.load_score_types(data)?;
        progress()?;
        self.load_processing_softwares(data)?;
        progress()?;
        self.load_db_search_params(data)?;
        progress()?;
        self.load_processing_steps(data)?;
        progress()?;
        self.load_observations(data)?;
        progress()?;
        self.load_parent_sequences(data)?;
        progress()?;
        self.load_parent_group_sets(data)?;
        progress()?;
        self.load_identified_compounds(data)?;
        progress()?;
        self.load_identified_sequences(data)?;
        progress()?;
        self.load_adducts(data)?;
        progress()?;
        self.load_observation_matches(data)
    }

    // ------------------------------------------------------------------
    // Feature and consensus maps
    // ------------------------------------------------------------------

    /// Source `loadMapMetaDataTemplate_`: unique ID, identifier, file path,
    /// file type and metadata of the map, plus the experiment type (version 5
    /// and later; empty otherwise). Absent table: nothing is set.
    #[allow(clippy::type_complexity)]
    fn load_map_meta_data(
        &self,
    ) -> Result<Option<(u64, String, String, FileType, MetaInfo, String)>> {
        if !self.table_exists("FEAT_MapMetaData")? {
            return Ok(None);
        }
        let mut query = self.prepare("SELECT * FROM FEAT_MapMetaData")?;
        let mut first = None;
        // "there should be only one row": the source reads the first.
        for_rows(&mut query, &[], |row| {
            if first.is_none() {
                let experiment_type = if self.version < 5 {
                    String::new()
                } else {
                    text(row, "experiment_type")?
                };
                first = Some((
                    int(row, "unique_id")?,
                    text(row, "identifier")?,
                    text(row, "file_path")?,
                    FileType::from_name(&text(row, "file_type")?),
                    experiment_type,
                ));
            }
            Ok(())
        })?;
        let Some((key, identifier, path, file_type, experiment_type)) = first else {
            return Err(invalid("OMS file has an empty FEAT_MapMetaData table"));
        };
        let mut metadata = MetaInfo::new();
        let mut meta = self.meta_query("FEAT_MapMetaData")?;
        Self::read_meta(
            self.options.keep_string_list_spaces,
            &mut meta,
            &mut metadata,
            key,
        )?;
        Ok(Some((
            key as u64,
            identifier,
            path,
            file_type,
            metadata,
            experiment_type,
        )))
    }

    /// Source `loadDataProcessing_`. An unknown action name is skipped, as the
    /// source does after logging it; the skipped names are returned.
    fn load_data_processing(&self, skipped: &mut Vec<String>) -> Result<Vec<DataProcessing>> {
        if !self.table_exists("FEAT_DataProcessing")? {
            return Ok(Vec::new());
        }
        let order_by = if self.version > 2 { "id" } else { "position" };
        let mut meta = self.meta_query("FEAT_DataProcessing")?;
        let mut query = self.prepare(&format!(
            "SELECT * FROM FEAT_DataProcessing ORDER BY {order_by} ASC"
        ))?;
        let mut result = Vec::new();
        for_rows(&mut query, &[], |row| {
            let mut proc = DataProcessing {
                software: Software {
                    name: text(row, "software_name")?,
                    version: text(row, "software_version")?,
                    ..Software::default()
                },
                ..DataProcessing::default()
            };
            for action in string_list(&text(row, "processing_actions")?)? {
                match ProcessingAction::ALL.iter().find(|a| a.name() == action) {
                    Some(&known) => {
                        proc.actions.insert(known);
                    }
                    None => skipped.push(action),
                }
            }
            let time = text(row, "completion_time")?;
            if time != NULL_DATE_TIME && !time.is_empty() {
                proc.completion_time = Some(crate::data_structures::DateTime::parse(&time)?);
            }
            if meta.is_some() {
                let key = int(row, "id")?;
                Self::read_meta(
                    self.options.keep_string_list_spaces,
                    &mut meta,
                    &mut proc.metadata,
                    key,
                )?;
            }
            result.push(proc);
            Ok(())
        })?;
        Ok(result)
    }

    /// Source `makeBaseFeature_`; a handle (`with_ids` false) stops after the
    /// measured values.
    fn make_base_feature(
        &self,
        row: &Row<'_>,
        with_ids: bool,
        meta: &mut Option<Statement<'c>>,
        match_query: &mut Option<Statement<'c>>,
    ) -> Result<BaseFeature> {
        let mut feature = BaseFeature::new(
            real(row, "rt")?,
            real(row, "mz")?,
            real(row, "intensity")? as f32,
        );
        feature.charge = int32(row, "charge")?;
        // Source `setWidth`, which also records the width as meta value "FWHM".
        feature.set_width(real(row, "width")? as f32)?;
        let quality_column = if self.version < 5 {
            "overall_quality"
        } else {
            "quality"
        };
        feature.quality = real(row, quality_column)? as f32;
        // The source stores `int64_t(getUniqueId())`: the bit pattern.
        feature.unique_id = int(row, "unique_id")? as u64;
        if !with_ids {
            return Ok(feature);
        }
        let id = int(row, "id")?;
        if let Some(primary) = optional_int(row, "primary_molecule_id")? {
            feature.primary_id = Some(lookup(&self.molecules, primary, "identified molecule")?);
        }
        Self::read_meta(
            self.options.keep_string_list_spaces,
            meta,
            &mut feature.metadata,
            id,
        )?;
        if let Some(query) = match_query.as_mut() {
            for_rows(query, &[id.into()], |sub| {
                let item = int(sub, "observation_match_id")?;
                feature
                    .id_matches
                    .insert(lookup(&self.matches, item, "observation match")?);
                Ok(())
            })?;
        }
        Ok(feature)
    }

    /// Source `prepareQueriesBaseFeature_`.
    fn base_feature_queries(&self) -> Result<(Option<Statement<'c>>, Option<Statement<'c>>)> {
        let main = if self.version < 5 {
            "FEAT_Feature"
        } else {
            "FEAT_BaseFeature"
        };
        let meta = self.meta_query(main)?;
        let matches = if self.table_exists("FEAT_ObservationMatch")? {
            Some(self.prepare("SELECT * FROM FEAT_ObservationMatch WHERE feature_id = ?1")?)
        } else {
            None
        };
        Ok((meta, matches))
    }

    fn feature_from(&self) -> &'static str {
        if self.version < 5 {
            "FEAT_Feature"
        } else {
            "FEAT_BaseFeature JOIN FEAT_Feature ON id = feature_id"
        }
    }

    /// Source `loadFeatureAndSubordinates_`.
    fn load_feature_and_subordinates(
        &self,
        row: &Row<'_>,
        meta: &mut Option<Statement<'c>>,
        match_query: &mut Option<Statement<'c>>,
        hull_query: &mut Option<Statement<'c>>,
        depth: usize,
        count: &mut usize,
    ) -> Result<Feature> {
        if depth > MAX_SUBORDINATE_DEPTH {
            return Err(invalid("OMS file nests subordinate features too deeply"));
        }
        *count += 1;
        if *count > MAX_TABLE_ROWS {
            return Err(invalid("OMS feature table exceeds the row limit"));
        }
        let id = int(row, "id")?;
        let mut feature = Feature {
            base: self.make_base_feature(row, true, meta, match_query)?,
            ..Feature::default()
        };
        feature.quality_rt = real(row, "rt_quality")? as f32;
        feature.quality_mz = real(row, "mz_quality")? as f32;
        if let Some(query) = hull_query.as_mut() {
            let mut hulls: Vec<ConvexHull2D> = Vec::new();
            for_rows(query, &[id.into()], |hull_row| {
                let index = size(hull_row, "hull_index")?;
                if index >= MAX_TABLE_ROWS {
                    return Err(invalid("OMS convex hull index out of range"));
                }
                if hulls.len() <= index {
                    hulls.resize_with(index + 1, ConvexHull2D::new);
                }
                hulls[index].add_point(Point2D::new(
                    real(hull_row, "point_x")?,
                    real(hull_row, "point_y")?,
                ))?;
                Ok(())
            })?;
            feature.convex_hulls = hulls;
        }
        // Subordinates, in key order, each loaded (recursively) while the
        // query over its siblings is open, as in the source.
        let mut subordinates = self.prepare(&format!(
            "SELECT * FROM {} WHERE subordinate_of = ?1 ORDER BY id ASC",
            self.feature_from()
        ))?;
        for_rows(&mut subordinates, &[id.into()], |sub_row| {
            let sub = self.load_feature_and_subordinates(
                sub_row,
                meta,
                match_query,
                hull_query,
                depth + 1,
                count,
            )?;
            feature.subordinates.push(sub);
            Ok(())
        })?;
        Ok(feature)
    }

    /// Source `loadFeatures_`.
    fn load_features(&self) -> Result<Vec<Feature>> {
        if !self.table_exists("FEAT_Feature")? {
            return Ok(Vec::new());
        }
        let (mut meta, mut match_query) = self.base_feature_queries()?;
        let mut hull_query = if self.table_exists("FEAT_ConvexHull")? {
            Some(self.prepare(
                "SELECT * FROM FEAT_ConvexHull WHERE feature_id = ?1 \
                 ORDER BY hull_index DESC, point_index ASC",
            )?)
        } else {
            None
        };
        let mut query = self.prepare(&format!(
            "SELECT * FROM {} WHERE subordinate_of IS NULL ORDER BY id ASC",
            self.feature_from()
        ))?;
        let mut features = Vec::new();
        let mut count = 0usize;
        for_rows(&mut query, &[], |row| {
            features.push(self.load_feature_and_subordinates(
                row,
                &mut meta,
                &mut match_query,
                &mut hull_query,
                0,
                &mut count,
            )?);
            Ok(())
        })?;
        Ok(features)
    }

    /// Source `load(FeatureMap&)` after the identification data.
    pub(super) fn load_feature_map(&self, skipped: &mut Vec<String>) -> Result<FeatureMap> {
        let mut map = FeatureMap::default();
        if let Some((unique_id, identifier, path, file_type, metadata, _)) =
            self.load_map_meta_data()?
        {
            map.unique_id = unique_id;
            map.identifier = identifier;
            map.loaded_file_path = path;
            map.loaded_file_type = self.file_type(file_type);
            map.metadata = metadata;
        }
        map.data_processing = self.load_data_processing(skipped)?;
        map.features = self.load_features()?;
        Ok(map)
    }

    /// The stored file type, or unknown as the source leaves it
    /// ([`OmsLoadOptions::unknown_file_type`]).
    fn file_type(&self, stored: FileType) -> FileType {
        if self.options.unknown_file_type {
            FileType::Unknown
        } else {
            stored
        }
    }

    /// Source `loadConsensusColumnHeaders_`.
    fn load_column_headers(&self) -> Result<BTreeMap<u64, ColumnHeader>> {
        let mut headers = BTreeMap::new();
        if !self.table_exists("FEAT_ConsensusColumnHeader")? {
            return Ok(headers);
        }
        let mut meta = self.meta_query("FEAT_ConsensusColumnHeader")?;
        let mut query = self.prepare("SELECT * FROM FEAT_ConsensusColumnHeader")?;
        for_rows(&mut query, &[], |row| {
            let id = int(row, "id")?;
            let mut header = ColumnHeader {
                filename: text(row, "filename")?,
                label: text(row, "label")?,
                size: usize::try_from(int(row, "size")?)
                    .map_err(|_| invalid("OMS column header size is negative"))?,
                unique_id: int(row, "unique_id")? as u64,
                ..ColumnHeader::default()
            };
            Self::read_meta(
                self.options.keep_string_list_spaces,
                &mut meta,
                &mut header.metadata,
                id,
            )?;
            headers.insert(id as u64, header);
            Ok(())
        })?;
        Ok(headers)
    }

    /// The ratios stored under one feature key, in index order.
    fn read_ratios(query: &mut Statement<'c>, key: i64) -> Result<Vec<Ratio>> {
        let mut ratios: Vec<Ratio> = Vec::new();
        for_rows(query, &[key.into()], |ratio_row| {
            let index = size(ratio_row, "ratio_index")?;
            if index >= MAX_TABLE_ROWS {
                return Err(invalid("OMS ratio index out of range"));
            }
            if ratios.len() <= index {
                ratios.resize_with(index + 1, Ratio::default);
            }
            ratios[index] = Ratio {
                ratio_value: real(ratio_row, "ratio_value")?,
                denominator_ref: text(ratio_row, "denominator_ref")?,
                numerator_ref: text(ratio_row, "numerator_ref")?,
                description: string_list(&text(ratio_row, "description")?)?,
            };
            Ok(())
        })?;
        Ok(ratios)
    }

    /// Source `loadConsensusFeatures_`.
    ///
    /// The source reads each consensus feature's ratios into a copy it has
    /// already appended to the map, so loaded consensus features never carry
    /// ratios there ([`OmsLoadOptions::drop_consensus_ratios`]). By default the
    /// port reads them: from the feature's own key, where the port's writer
    /// puts them, or else from the key of its last handle, where the source's
    /// writer puts them.
    fn load_consensus_features(&self) -> Result<Vec<ConsensusFeature>> {
        if !self.table_exists("FEAT_FeatureHandle")? {
            return Ok(Vec::new());
        }
        let (mut meta, mut match_query) = self.base_feature_queries()?;
        let mut query = self.prepare(
            "SELECT * FROM FEAT_BaseFeature LEFT JOIN FEAT_FeatureHandle ON id = feature_id \
             ORDER BY id ASC",
        )?;
        let mut features: Vec<ConsensusFeature> = Vec::new();
        // Per consensus feature: its own key and the key of its last handle.
        let mut keys: Vec<(i64, Option<i64>)> = Vec::new();
        for_rows(&mut query, &[], |row| {
            let id = int(row, "id")?;
            if is_null(row, "subordinate_of")? {
                let base = self.make_base_feature(row, true, &mut meta, &mut match_query)?;
                let mut feature = ConsensusFeature::new();
                feature.base = base;
                features.push(feature);
                keys.push((id, None));
            } else {
                let base = self.make_base_feature(row, false, &mut meta, &mut match_query)?;
                let map_index = int(row, "map_index")? as u64;
                let (Some(last), Some(last_keys)) = (features.last_mut(), keys.last_mut()) else {
                    return Err(invalid(
                        "OMS feature handle precedes every consensus feature",
                    ));
                };
                last.insert(FeatureHandle::new(map_index, &base))?;
                last_keys.1 = Some(id);
            }
            Ok(())
        })?;
        if !self.options.drop_consensus_ratios && self.table_exists("FEAT_ConsensusRatio")? {
            let mut ratio_query = self.prepare(
                "SELECT * FROM FEAT_ConsensusRatio WHERE feature_id = ?1 ORDER BY ratio_index DESC",
            )?;
            for (feature, (own, last_handle)) in features.iter_mut().zip(keys) {
                let mut ratios = Self::read_ratios(&mut ratio_query, own)?;
                if ratios.is_empty() {
                    if let Some(handle) = last_handle {
                        ratios = Self::read_ratios(&mut ratio_query, handle)?;
                    }
                }
                feature.ratios = ratios;
            }
        }
        Ok(features)
    }

    /// Source `load(ConsensusMap&)` after the identification data.
    pub(super) fn load_consensus_map(&self, skipped: &mut Vec<String>) -> Result<ConsensusMap> {
        let mut map = ConsensusMap::default();
        if let Some((unique_id, identifier, path, file_type, metadata, experiment_type)) =
            self.load_map_meta_data()?
        {
            map.unique_id = unique_id;
            map.identifier = identifier;
            map.loaded_file_path = path;
            map.loaded_file_type = self.file_type(file_type);
            map.metadata = metadata;
            map.experiment_type = experiment_type;
        }
        map.column_headers = self.load_column_headers()?;
        map.data_processing = self.load_data_processing(skipped)?;
        map.features = self.load_consensus_features()?;
        Ok(map)
    }
}

/// `ID_MoleculeType` key to molecule type.
fn molecule_type(key: i64) -> Result<MoleculeType> {
    match key {
        1 => Ok(MoleculeType::Protein),
        2 => Ok(MoleculeType::Compound),
        3 => Ok(MoleculeType::RNA),
        _ => Err(invalid(format!(
            "OMS file names unknown molecule type {key}"
        ))),
    }
}

/// A processing step's date: empty or the unset sentinel is absent.
fn parse_date_time(value: &str) -> Result<Option<CompletionTime>> {
    if value.is_empty() || value == NULL_DATE_TIME {
        return Ok(None);
    }
    value.parse().map(Some)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn meta_values_parse_like_make_data_value() {
        assert_eq!(
            make_meta_value(1, "a, b", true).unwrap(),
            MetaValue::from("a, b")
        );
        assert_eq!(
            make_meta_value(2, "42", true).unwrap(),
            MetaValue::from(42i64)
        );
        assert_eq!(
            make_meta_value(3, "3.688524e06", true).unwrap(),
            MetaValue::try_from(3_688_524.0).unwrap()
        );
        // String lists keep the space after each comma, as in the source.
        assert_eq!(
            make_meta_value(4, "[a, b]", true).unwrap(),
            MetaValue::from(vec!["a".to_owned(), " b".to_owned()])
        );
        assert_eq!(
            make_meta_value(5, "[1, -2]", true).unwrap(),
            MetaValue::from(vec![1i64, -2])
        );
        assert_eq!(
            make_meta_value(4, "[]", true).unwrap(),
            MetaValue::from(Vec::<String>::new())
        );
        assert!(make_meta_value(0, "ignored", true).unwrap().is_empty());
        assert!(make_meta_value(2, "x", true).is_err());
        // Natively the one separator space is removed, so the list survives.
        assert_eq!(
            make_meta_value(4, "[a, b,  c]", false).unwrap(),
            MetaValue::from(vec!["a".to_owned(), "b".to_owned(), " c".to_owned()])
        );
    }

    #[test]
    fn unset_dates_are_absent() {
        assert_eq!(parse_date_time("").unwrap(), None);
        assert_eq!(parse_date_time(NULL_DATE_TIME).unwrap(), None);
        assert_eq!(
            parse_date_time("2022-04-14 13:03:26")
                .unwrap()
                .unwrap()
                .to_string(),
            "2022-04-14 13:03:26"
        );
    }
}
