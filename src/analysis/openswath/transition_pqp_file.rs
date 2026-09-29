// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Reading and writing of PQP transition libraries.
//!
//! Ports `ANALYSIS/OPENSWATH/TransitionPQPFile.h` (with its `.cpp`). PQP
//! files are SQLite databases whose tables represent the data of a TraML
//! file; this type converts between them and a [`TargetedExperiment`]. For
//! another format that stores transitions see
//! [`TransitionTSVFile`], whose parameters and conversion this type shares
//! (the source class derives from it).
//!
//! # Tables
//!
//! * `GENE` (`ID`, `GENE_NAME`, `DECOY`) and `PEPTIDE_GENE_MAPPING`
//!   (`PEPTIDE_ID`, `GENE_ID`): genes by a primary key and a human-readable
//!   name.
//! * `PROTEIN` (`ID`, `PROTEIN_ACCESSION`, `DECOY`) and
//!   `PEPTIDE_PROTEIN_MAPPING` (`PEPTIDE_ID`, `PROTEIN_ID`).
//! * `PEPTIDE` (`ID`, `UNMODIFIED_SEQUENCE`, `MODIFIED_SEQUENCE`, `DECOY`) and
//!   `PRECURSOR_PEPTIDE_MAPPING` (`PRECURSOR_ID`, `PEPTIDE_ID`): the physical
//!   analytes with their post-translational modifications. Modifications are
//!   written inside the sequence as UniMod identifiers or names OpenMS
//!   understands, e.g. `PEPT(Phosphorylation)IDE(UniMod:27)A`.
//! * `COMPOUND` (`ID`, `COMPOUND_NAME`, `SUM_FORMULA`, `SMILES`, `ADDUCTS`,
//!   `DECOY`) and `PRECURSOR_COMPOUND_MAPPING` (`PRECURSOR_ID`,
//!   `COMPOUND_ID`): small-molecule analytes.
//! * `PRECURSOR` (`ID`, `TRAML_ID` — the TraML `id=` attribute,
//!   `GROUP_LABEL` — the peptide label group of `MS:1000893`, `PRECURSOR_MZ`,
//!   `CHARGE`, `LIBRARY_INTENSITY`, `LIBRARY_RT`, `LIBRARY_DRIFT_TIME`,
//!   `DECOY`) and `TRANSITION_PRECURSOR_MAPPING` (`TRANSITION_ID`,
//!   `PRECURSOR_ID`): ions from ionisation.
//! * `TRANSITION` (`ID`, `TRAML_ID`, `PRODUCT_MZ`, `CHARGE`, `TYPE` (e.g.
//!   `b`, `y`), `ANNOTATION`, `ORDINAL` (6 for y6), `DETECTING`,
//!   `IDENTIFYING`, `QUANTIFYING`, `LIBRARY_INTENSITY`, `DECOY`): fragments.
//!   Detecting transitions detect a precursor (use them for the top N);
//!   identifying transitions discriminate peptidoforms of one precursor in
//!   the IPF workflow.
//! * `TRANSITION_PEPTIDE_MAPPING` (`TRANSITION_ID`, `PEPTIDE_ID`): the IPF
//!   mapping of one identifying transition to its peptidoforms.
//! * `VERSION` (`ID`): the file format version, 3.
//!
//! Peptide label groups designate isotopically modified forms of the same
//! peptide (`PEPTIDEAK` and `PEPTIDEAK[+8]` share `PEPTIDEAK_gr1`).
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/TRANSITION_PQP_SUPPORT.md)
//! for the member mapping and the evidence. The heavy read path is serial,
//! as in the source; the source's `PRAGMA threads` hint for its streaming
//! reader has no counterpart. The `OpenSwath::LightTargetedExperiment`
//! overloads wait for `OPENSWATHALGO/DATAACCESS/TransitionExperiment.h` and
//! `OpenSwathLibraryIDNormalizer`.

use crate::analysis::openswath::targeted_experiment::TargetedExperiment;
use crate::analysis::openswath::targeted_experiment_helper::{
    Peptide, UnresolvedModificationPolicy, get_aa_sequence_with,
};
use crate::analysis::openswath::transition_tsv_file::{
    TransitionTSVFile, TsvTransition, convert_transition,
};
use crate::chemistry::{AASequence, ModificationsDB};
use crate::format::file_info::text_format::ostream_g;
use crate::format::sqlite_connector::{SqlOpenMode, SqliteConnector};
use crate::metadata::{MetaValue, MetaValueData};
use crate::param::Param;
use crate::{Error, Result};
use rusqlite::types::ValueRef;
use rusqlite::{Connection, params};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::Path;

fn sql(error: rusqlite::Error) -> Error {
    Error::Io(std::io::Error::other(error))
}

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

/// The PQP reader and writer. Its parameters are those of
/// [`TransitionTSVFile`] (`retentionTimeInterpretation`,
/// `override_group_label_check`, `force_invalid_mods`), which the source
/// inherits.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct TransitionPQPFile {
    /// The shared TSV conversion and its parameters.
    pub tsv: TransitionTSVFile,
}

/// Column availability found while building the select query
/// (`PQPSqlQueryInfo`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PqpSqlQueryInfo {
    /// The complete `SELECT` statement.
    pub select_sql: String,
    /// Whether `PRECURSOR.LIBRARY_DRIFT_TIME` exists.
    pub drift_time_exists: bool,
    /// Whether a non-empty `GENE` table exists.
    pub gene_exists: bool,
    /// Whether compound mappings exist.
    pub compound_exists: bool,
    /// Whether transition-to-peptide mappings exist.
    pub peptidoforms_exists: bool,
}

/// `sqlite3_column_double`: an integer converts, text is read as a number.
fn column_double(value: ValueRef<'_>) -> Result<Option<f64>> {
    Ok(match value {
        ValueRef::Null => None,
        ValueRef::Integer(i) => Some(i as f64),
        ValueRef::Real(x) => Some(x),
        ValueRef::Text(t) => Some(
            std::str::from_utf8(t)
                .ok()
                .and_then(|s| s.trim().parse::<f64>().ok())
                .ok_or_else(|| bad("PQP text column is not a number"))?,
        ),
        ValueRef::Blob(_) => return Err(bad("PQP blob column where a number is expected")),
    })
}

/// `sqlite3_column_int`: the 32-bit value of an integer column; a real is
/// truncated as SQLite does.
fn column_int(value: ValueRef<'_>) -> Result<Option<i32>> {
    Ok(match value {
        ValueRef::Null => None,
        ValueRef::Integer(i) => {
            Some(i32::try_from(i).map_err(|_| bad("PQP integer column exceeds 32 bits"))?)
        }
        ValueRef::Real(x) => Some(x.trunc() as i32),
        ValueRef::Text(t) => Some(
            std::str::from_utf8(t)
                .ok()
                .and_then(|s| s.trim().parse::<i32>().ok())
                .ok_or_else(|| bad("PQP text column is not an integer"))?,
        ),
        ValueRef::Blob(_) => return Err(bad("PQP blob column where an integer is expected")),
    })
}

/// `sqlite3_column_text`: integers in decimal, reals as SQLite renders them
/// (15 significant digits, always with a decimal point or exponent).
fn column_text(value: ValueRef<'_>) -> Result<Option<String>> {
    Ok(match value {
        ValueRef::Null => None,
        ValueRef::Integer(i) => Some(i.to_string()),
        ValueRef::Real(x) => {
            let mut text = ostream_g(x, 15);
            if x.is_finite() && !text.contains(['.', 'e']) {
                text.push_str(".0");
            }
            Some(text)
        }
        ValueRef::Text(t) => Some(
            std::str::from_utf8(t)
                .map_err(|_| bad("PQP text column is not UTF-8"))?
                .to_owned(),
        ),
        ValueRef::Blob(b) => Some(String::from_utf8_lossy(b).into_owned()),
    })
}

/// `extractValueIntStr`: only an integer column sets the text.
fn column_int_text(value: ValueRef<'_>) -> Result<Option<String>> {
    Ok(match value {
        ValueRef::Integer(i) => Some(
            i32::try_from(i)
                .map_err(|_| bad("PQP charge exceeds 32 bits"))?
                .to_string(),
        ),
        _ => None,
    })
}

fn split(text: &str, separator: char) -> Vec<String> {
    if text.is_empty() {
        Vec::new()
    } else {
        text.split(separator).map(str::to_owned).collect()
    }
}

fn has_rows(connection: &Connection, table: &str) -> Result<bool> {
    // `table` is one of this module's constant names.
    let mut statement = connection
        .prepare(&format!("SELECT 1 FROM {table} LIMIT 1;"))
        .map_err(sql)?;
    let mut rows = statement.query([]).map_err(sql)?;
    Ok(rows.next().map_err(sql)?.is_some())
}

fn validate_traml_id_table_name(table: &str) -> Result<()> {
    if table != "PRECURSOR" && table != "TRANSITION" {
        return Err(bad("tableName must be PRECURSOR or TRANSITION"));
    }
    Ok(())
}

/// `DataValue::toString()` for the metadata the writer reads.
fn value_text(value: &MetaValue) -> String {
    match value.data() {
        MetaValueData::Empty => String::new(),
        MetaValueData::String(s) => s.clone(),
        MetaValueData::Float(x) => crate::param::value::format_float(*x, true),
        _ => value.to_string(),
    }
}

impl TransitionPQPFile {
    /// A reader/writer with the default parameters.
    pub fn new() -> Self {
        Self::default()
    }

    /// The default parameters, those of [`TransitionTSVFile::defaults`].
    ///
    /// # Errors
    ///
    /// As [`TransitionTSVFile::defaults`].
    pub fn defaults() -> Result<Param> {
        TransitionTSVFile::defaults()
    }

    /// The current parameters.
    ///
    /// # Errors
    ///
    /// As [`TransitionTSVFile::to_param`].
    pub fn to_param(&self) -> Result<Param> {
        self.tsv.to_param()
    }

    /// Build from parameters, as [`TransitionTSVFile::from_param`].
    ///
    /// # Errors
    ///
    /// As [`TransitionTSVFile::from_param`].
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self {
            tsv: TransitionTSVFile::from_param(parameters)?,
        })
    }

    /// As [`from_param`](Self::from_param), with the unknown-parameter
    /// warnings.
    ///
    /// # Errors
    ///
    /// As [`TransitionTSVFile::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let (tsv, warnings) = TransitionTSVFile::from_param_with_warnings(parameters)?;
        Ok((Self { tsv }, warnings))
    }

    /// Build the select statement for reading transitions
    /// (`buildPQPSelectQuery_`).
    ///
    /// With `legacy_traml_id` the TraML identifiers (`TRAML_ID`) name the
    /// precursors and transitions, otherwise their numeric `ID`. Optional
    /// columns and tables are probed: `LIBRARY_DRIFT_TIME`, a non-empty
    /// `GENE` table (an empty one would drop every row through its inner
    /// join), `TRANSITION.ANNOTATION`, `COMPOUND.ADDUCTS`, and non-empty
    /// compound and peptidoform mappings. With `stable_order` the rows are
    /// ordered by precursor m/z, product m/z, retention time, transition
    /// name and group id.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] on a SQLite error.
    pub fn build_pqp_select_query(
        conn: &SqliteConnector,
        legacy_traml_id: bool,
        stable_order: bool,
    ) -> Result<PqpSqlQueryInfo> {
        let db = conn.connection();
        let traml_id = if legacy_traml_id { "TRAML_ID" } else { "ID" };
        let mut info = PqpSqlQueryInfo {
            drift_time_exists: conn.column_exists("PRECURSOR", "LIBRARY_DRIFT_TIME")?,
            ..PqpSqlQueryInfo::default()
        };
        let select_drift_time = if info.drift_time_exists {
            ", PRECURSOR.LIBRARY_DRIFT_TIME AS drift_time "
        } else {
            ""
        };
        info.gene_exists = conn.table_exists("GENE")? && has_rows(db, "GENE")?;
        let (select_gene, select_gene_null, join_gene) = if info.gene_exists {
            (
                ", GENE_AGGREGATED.GENE_NAME AS gene_name ",
                ", 'NA' AS gene_name ",
                "INNER JOIN (SELECT PEPTIDE_ID, GROUP_CONCAT(GENE_NAME,';') AS GENE_NAME FROM GENE INNER JOIN PEPTIDE_GENE_MAPPING ON GENE.ID = PEPTIDE_GENE_MAPPING.GENE_ID GROUP BY PEPTIDE_ID) AS GENE_AGGREGATED ON PEPTIDE.ID = GENE_AGGREGATED.PEPTIDE_ID ",
            )
        } else {
            ("", "", "")
        };
        let select_annotation = if conn.column_exists("TRANSITION", "ANNOTATION")? {
            "TRANSITION.ANNOTATION AS Annotation, "
        } else {
            "'' AS Annotation, "
        };
        let select_adducts = if conn.column_exists("COMPOUND", "ADDUCTS")? {
            "COMPOUND.ADDUCTS AS Adducts, "
        } else {
            "'' AS Adducts, "
        };
        info.compound_exists = conn.table_exists("PRECURSOR_COMPOUND_MAPPING")?
            && has_rows(db, "PRECURSOR_COMPOUND_MAPPING")?;
        info.peptidoforms_exists = conn.table_exists("TRANSITION_PEPTIDE_MAPPING")?
            && has_rows(db, "TRANSITION_PEPTIDE_MAPPING")?;
        let (select_peptidoforms, join_peptidoforms) = if info.peptidoforms_exists {
            (
                "PEPTIDE_AGGREGATED.PEPTIDOFORMS AS peptidoforms ",
                "LEFT OUTER JOIN (SELECT TRANSITION_ID, GROUP_CONCAT(MODIFIED_SEQUENCE,'|') AS PEPTIDOFORMS FROM TRANSITION_PEPTIDE_MAPPING INNER JOIN PEPTIDE ON TRANSITION_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID GROUP BY TRANSITION_ID) AS PEPTIDE_AGGREGATED ON TRANSITION.ID = PEPTIDE_AGGREGATED.TRANSITION_ID ",
            )
        } else {
            ("NULL AS peptidoforms ", "")
        };
        let mut select = format!(
            "SELECT PRECURSOR.PRECURSOR_MZ AS precursor, TRANSITION.PRODUCT_MZ AS product, PRECURSOR.LIBRARY_RT AS rt_calibrated, TRANSITION.{traml_id} AS transition_name, -1 AS CE, TRANSITION.LIBRARY_INTENSITY AS library_intensity, PRECURSOR.{traml_id} AS group_id, TRANSITION.DECOY AS decoy, PEPTIDE.UNMODIFIED_SEQUENCE AS PeptideSequence, PROTEIN_AGGREGATED.PROTEIN_ACCESSION AS ProteinName, {select_annotation}PEPTIDE.MODIFIED_SEQUENCE AS FullPeptideName, NULL AS CompoundName, NULL AS SMILES, NULL AS SumFormula, NULL AS Adducts, PRECURSOR.CHARGE AS precursor_charge, PRECURSOR.GROUP_LABEL AS peptide_group_label, NULL AS label_type, TRANSITION.CHARGE AS fragment_charge, TRANSITION.ORDINAL AS fragment_nr, NULL AS fragment_mzdelta, NULL AS fragment_modification, TRANSITION.TYPE AS fragment_type, NULL AS uniprot_id, TRANSITION.DETECTING AS detecting_transition, TRANSITION.IDENTIFYING AS identifying_transition, TRANSITION.QUANTIFYING AS quantifying_transition, {select_peptidoforms}{select_drift_time}{select_gene}FROM PRECURSOR {join_gene}INNER JOIN TRANSITION_PRECURSOR_MAPPING ON PRECURSOR.ID = TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID INNER JOIN TRANSITION ON TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID = TRANSITION.ID INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID INNER JOIN (SELECT PEPTIDE_ID, GROUP_CONCAT(PROTEIN_ACCESSION,';') AS PROTEIN_ACCESSION FROM PROTEIN INNER JOIN PEPTIDE_PROTEIN_MAPPING ON PROTEIN.ID = PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID GROUP BY PEPTIDE_ID) AS PROTEIN_AGGREGATED ON PEPTIDE.ID = PROTEIN_AGGREGATED.PEPTIDE_ID {join_peptidoforms}"
        );
        if info.compound_exists {
            select.push_str(&format!(
                "UNION ALL SELECT PRECURSOR.PRECURSOR_MZ AS precursor, TRANSITION.PRODUCT_MZ AS product, PRECURSOR.LIBRARY_RT AS rt_calibrated, TRANSITION.{traml_id} AS transition_name, -1 AS CE, TRANSITION.LIBRARY_INTENSITY AS library_intensity, PRECURSOR.{traml_id} AS group_id, TRANSITION.DECOY AS decoy, NULL AS PeptideSequence, NULL AS ProteinName, {select_annotation}NULL AS FullPeptideName, COMPOUND.COMPOUND_NAME AS CompoundName, COMPOUND.SMILES AS SMILES, COMPOUND.SUM_FORMULA AS SumFormula, {select_adducts}PRECURSOR.CHARGE AS precursor_charge, PRECURSOR.GROUP_LABEL AS peptide_group_label, NULL AS label_type, TRANSITION.CHARGE AS fragment_charge, TRANSITION.ORDINAL AS fragment_nr, NULL AS fragment_mzdelta, NULL AS fragment_modification, TRANSITION.TYPE AS fragment_type, NULL AS uniprot_id, TRANSITION.DETECTING AS detecting_transition, TRANSITION.IDENTIFYING AS identifying_transition, TRANSITION.QUANTIFYING AS quantifying_transition, NULL AS peptidoforms {select_drift_time}{select_gene_null}FROM PRECURSOR INNER JOIN TRANSITION_PRECURSOR_MAPPING ON PRECURSOR.ID = TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID INNER JOIN TRANSITION ON TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID = TRANSITION.ID INNER JOIN PRECURSOR_COMPOUND_MAPPING ON PRECURSOR.ID = PRECURSOR_COMPOUND_MAPPING.PRECURSOR_ID INNER JOIN COMPOUND ON PRECURSOR_COMPOUND_MAPPING.COMPOUND_ID = COMPOUND.ID "
            ));
        }
        if stable_order {
            select
                .push_str("ORDER BY precursor, product, rt_calibrated, transition_name, group_id ");
        }
        select.push_str("; ");
        info.select_sql = select;
        Ok(info)
    }

    /// Read the transitions of a PQP file (`readPQPInput_`).
    ///
    /// The file is opened read-only. Rows come in the stable order of
    /// [`build_pqp_select_query`](Self::build_pqp_select_query); reading
    /// stops at the first row whose precursor m/z is `NULL`, as the source's
    /// loop condition does. A gene name `NA` reads as empty.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the database cannot be opened or queried (the
    /// source's `Exception::SqlOperationFailed`); [`Error::InvalidValue`]
    /// for a column value of an unusable type, or when the transition count
    /// exceeds the TSV reader's limit.
    pub fn read_pqp_input(
        &self,
        path: impl AsRef<Path>,
        legacy_traml_id: bool,
    ) -> Result<Vec<TsvTransition>> {
        let conn = SqliteConnector::open(path, SqlOpenMode::ReadOnly)?;
        let count = conn.count_table_rows("TRANSITION")?;
        if count > self.tsv.limits.max_transitions {
            return Err(bad("PQP file exceeds the transition limit"));
        }
        let info = Self::build_pqp_select_query(&conn, legacy_traml_id, true)?;
        let mut statement = conn.connection().prepare(&info.select_sql).map_err(sql)?;
        let mut rows = statement.query([]).map_err(sql)?;
        let mut transitions = Vec::with_capacity(count);
        while let Some(row) = rows.next().map_err(sql)? {
            let at = |i: usize| row.get_ref(i).map_err(sql);
            if matches!(at(0)?, ValueRef::Null) {
                break;
            }
            if transitions.len() >= self.tsv.limits.max_transitions {
                return Err(bad("PQP file exceeds the transition limit"));
            }
            let mut t = TsvTransition::default();
            macro_rules! set {
                ($field:expr, $f:ident, $i:expr) => {
                    if let Some(v) = $f(at($i)?)? {
                        $field = v;
                    }
                };
            }
            set!(t.precursor, column_double, 0);
            set!(t.product, column_double, 1);
            set!(t.rt_calibrated, column_double, 2);
            set!(t.transition_name, column_text, 3);
            set!(t.ce, column_double, 4);
            set!(t.library_intensity, column_double, 5);
            set!(t.group_id, column_text, 6);
            if let Some(v) = column_int(at(7)?)? {
                t.decoy = v != 0;
            }
            set!(t.peptide_sequence, column_text, 8);
            if let Some(v) = column_text(at(9)?)? {
                t.protein_name = split(&v, ';');
            }
            set!(t.annotation, column_text, 10);
            set!(t.full_peptide_name, column_text, 11);
            set!(t.compound_name, column_text, 12);
            set!(t.smiles, column_text, 13);
            set!(t.sum_formula, column_text, 14);
            set!(t.adducts, column_text, 15);
            set!(t.precursor_charge, column_int_text, 16);
            set!(t.peptide_group_label, column_text, 17);
            set!(t.label_type, column_text, 18);
            set!(t.fragment_charge, column_int_text, 19);
            set!(t.fragment_nr, column_int, 20);
            set!(t.fragment_mzdelta, column_double, 21);
            set!(t.fragment_modification, column_int, 22);
            set!(t.fragment_type, column_text, 23);
            if let Some(v) = column_text(at(24)?)? {
                t.uniprot_id = split(&v, ';');
            }
            if let Some(v) = column_int(at(25)?)? {
                t.detecting_transition = v != 0;
            }
            if let Some(v) = column_int(at(26)?)? {
                t.identifying_transition = v != 0;
            }
            if let Some(v) = column_int(at(27)?)? {
                t.quantifying_transition = v != 0;
            }
            if let Some(v) = column_text(at(28)?)? {
                t.peptidoforms = split(&v, '|');
            }
            let mut optional = 29;
            if info.drift_time_exists {
                set!(t.drift_time, column_double, optional);
                optional += 1;
            }
            if info.gene_exists {
                set!(t.gene_name, column_text, optional);
            }
            if t.gene_name == "NA" {
                t.gene_name.clear();
            }
            transitions.push(t);
        }
        Ok(transitions)
    }

    /// Read a PQP file into a [`TargetedExperiment`]
    /// (`convertPQPToTargetedExperiment`, the TraML-structure overload).
    ///
    /// With `legacy_traml_id` the TraML identifiers name precursors and
    /// transitions, otherwise the numeric PQP ids.
    ///
    /// # Errors
    ///
    /// As [`read_pqp_input`](Self::read_pqp_input) and
    /// [`TransitionTSVFile::to_targeted_experiment`].
    pub fn convert_pqp_to_targeted_experiment(
        &self,
        path: impl AsRef<Path>,
        legacy_traml_id: bool,
    ) -> Result<(TargetedExperiment, Vec<String>)> {
        let mut transitions = self.read_pqp_input(path, legacy_traml_id)?;
        let mut warnings = Vec::new();
        let exp = self
            .tsv
            .to_targeted_experiment(&mut transitions, &mut warnings)?;
        Ok((exp, warnings))
    }

    /// Validate an experiment (inherited `validateTargetedExperiment`).
    ///
    /// # Errors
    ///
    /// As [`TransitionTSVFile::validate_targeted_experiment`].
    pub fn validate_targeted_experiment(&self, exp: &TargetedExperiment) -> Result<()> {
        self.tsv.validate_targeted_experiment(exp)
    }

    /// Map each `TRAML_ID` of a table to its persistent `ID`
    /// (`getPQPIDToTraMLIDMap`); a repeated `TRAML_ID` keeps its last row.
    ///
    /// The file is opened read-only. The source opens it in the connector's
    /// default read-write-or-create mode, which creates an empty database
    /// for a missing file and then fails on the query; this port fails on
    /// the open and creates nothing.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `table` is not `PRECURSOR` or
    /// `TRANSITION` (the source's `Exception::IllegalArgument`); [`Error::Io`]
    /// on a SQLite error.
    pub fn pqp_id_to_traml_id_map(
        path: impl AsRef<Path>,
        table: &str,
    ) -> Result<HashMap<String, String>> {
        validate_traml_id_table_name(table)?;
        let conn = SqliteConnector::open(path, SqlOpenMode::ReadOnly)?;
        let mut out = HashMap::new();
        for (id, traml_id) in id_pairs(&conn, table)? {
            out.insert(traml_id, id);
        }
        Ok(out)
    }

    /// Map each persistent `ID` of a table to its `TRAML_ID`
    /// (`getPQPCurrentIDToTraMLIDMap`); an empty or `NULL` `TRAML_ID` maps to
    /// the empty string, and a repeated `ID` keeps its first row.
    ///
    /// # Errors
    ///
    /// As [`pqp_id_to_traml_id_map`](Self::pqp_id_to_traml_id_map).
    pub fn pqp_current_id_to_traml_id_map(
        path: impl AsRef<Path>,
        table: &str,
    ) -> Result<HashMap<String, String>> {
        validate_traml_id_table_name(table)?;
        let conn = SqliteConnector::open(path, SqlOpenMode::ReadOnly)?;
        let mut out = HashMap::new();
        for (id, traml_id) in id_pairs(&conn, table)? {
            out.entry(id).or_insert(traml_id);
        }
        Ok(out)
    }

    /// Write an experiment as a PQP file (`convertTargetedExperimentToPQP`).
    ///
    /// An existing file at `path` is removed first, as in the source.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the experiment contains invalid
    /// references ("Your input file contains invalid references, cannot
    /// process file."), when a peptide or compound has no retention time
    /// (the source's `Exception::IllegalArgument` from `getRetentionTime`),
    /// or when a modified sequence cannot be built; [`Error::Io`] on a file
    /// or SQLite error.
    pub fn convert_targeted_experiment_to_pqp(
        &self,
        path: impl AsRef<Path>,
        exp: &TargetedExperiment,
    ) -> Result<()> {
        if exp.contains_invalid_references() {
            return Err(bad(
                "Your input file contains invalid references, cannot process file.",
            ));
        }
        write_pqp_output(path.as_ref(), exp)
    }
}

fn id_pairs(conn: &SqliteConnector, table: &str) -> Result<Vec<(String, String)>> {
    let mut statement = conn
        .connection()
        .prepare(&format!("SELECT ID, TRAML_ID FROM {table};"))
        .map_err(sql)?;
    let mut rows = statement.query([]).map_err(sql)?;
    let mut out = Vec::new();
    while let Some(row) = rows.next().map_err(sql)? {
        let id = column_text(row.get_ref(0).map_err(sql)?)?.unwrap_or_default();
        let traml = column_text(row.get_ref(1).map_err(sql)?)?.unwrap_or_default();
        out.push((id, traml));
    }
    Ok(out)
}

const CREATE_SQL: &str = "CREATE TABLE VERSION(ID INT NOT NULL);\
CREATE TABLE GENE(ID INT PRIMARY KEY NOT NULL,GENE_NAME TEXT NOT NULL,DECOY INT NOT NULL);\
CREATE TABLE PEPTIDE_GENE_MAPPING(PEPTIDE_ID INT NOT NULL,GENE_ID INT NOT NULL);\
CREATE TABLE PROTEIN(ID INT PRIMARY KEY NOT NULL,PROTEIN_ACCESSION TEXT NOT NULL,DECOY INT NOT NULL);\
CREATE TABLE PEPTIDE_PROTEIN_MAPPING(PEPTIDE_ID INT NOT NULL,PROTEIN_ID INT NOT NULL);\
CREATE TABLE PEPTIDE(ID INT PRIMARY KEY NOT NULL,UNMODIFIED_SEQUENCE TEXT NOT NULL,MODIFIED_SEQUENCE TEXT NOT NULL,DECOY INT NOT NULL);\
CREATE TABLE PRECURSOR_PEPTIDE_MAPPING(PRECURSOR_ID INT NOT NULL,PEPTIDE_ID INT NOT NULL);\
CREATE TABLE COMPOUND(ID INT PRIMARY KEY NOT NULL,COMPOUND_NAME TEXT NOT NULL,SUM_FORMULA TEXT NOT NULL,SMILES TEXT NOT NULL,ADDUCTS TEXT NOT NULL,DECOY INT NOT NULL);\
CREATE TABLE PRECURSOR_COMPOUND_MAPPING(PRECURSOR_ID INT NOT NULL,COMPOUND_ID INT NOT NULL);\
CREATE TABLE PRECURSOR(ID INT PRIMARY KEY NOT NULL,TRAML_ID TEXT NULL,GROUP_LABEL TEXT NULL,PRECURSOR_MZ REAL NOT NULL,CHARGE INT NULL,LIBRARY_INTENSITY REAL NULL,LIBRARY_RT REAL NULL,LIBRARY_DRIFT_TIME REAL NULL,DECOY INT NOT NULL);\
CREATE TABLE TRANSITION_PRECURSOR_MAPPING(TRANSITION_ID INT NOT NULL,PRECURSOR_ID INT NOT NULL);\
CREATE TABLE TRANSITION_PEPTIDE_MAPPING(TRANSITION_ID INT NOT NULL,PEPTIDE_ID INT NOT NULL);\
CREATE TABLE TRANSITION(ID INT PRIMARY KEY NOT NULL,TRAML_ID TEXT NULL,PRODUCT_MZ REAL NOT NULL,CHARGE INT NULL,TYPE CHAR(1) NULL,ANNOTATION TEXT NULL,ORDINAL INT NULL,DETECTING INT NOT NULL,IDENTIFYING INT NOT NULL,QUANTIFYING INT NOT NULL,LIBRARY_INTENSITY REAL NULL,DECOY INT NOT NULL);";

const UPDATE_DECOYS_SQL: &str = "UPDATE PEPTIDE SET DECOY = 1 WHERE ID IN (SELECT PEPTIDE.ID FROM PRECURSOR JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID WHERE PRECURSOR.DECOY = 1); \
UPDATE COMPOUND SET DECOY = 1 WHERE ID IN (SELECT COMPOUND.ID FROM PRECURSOR JOIN PRECURSOR_COMPOUND_MAPPING ON PRECURSOR.ID = PRECURSOR_COMPOUND_MAPPING.PRECURSOR_ID JOIN COMPOUND ON PRECURSOR_COMPOUND_MAPPING.COMPOUND_ID = COMPOUND.ID WHERE PRECURSOR.DECOY = 1); \
UPDATE PROTEIN SET DECOY = 1 WHERE ID IN (SELECT PROTEIN.ID FROM PEPTIDE JOIN PEPTIDE_PROTEIN_MAPPING ON PEPTIDE.ID = PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID JOIN PROTEIN ON PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID = PROTEIN.ID WHERE PEPTIDE.DECOY = 1); \
UPDATE GENE SET DECOY = 1 WHERE ID IN (SELECT GENE.ID FROM PEPTIDE JOIN PEPTIDE_GENE_MAPPING ON PEPTIDE.ID = PEPTIDE_GENE_MAPPING.PEPTIDE_ID JOIN GENE ON PEPTIDE_GENE_MAPPING.GENE_ID = GENE.ID WHERE PEPTIDE.DECOY = 1); ";

/// Sorted, deduplicated values with their rank as id (the source's
/// `boost::sort` + `boost::unique` + index map).
fn ranked(mut values: Vec<String>) -> (Vec<String>, BTreeMap<String, i64>) {
    values.sort();
    values.dedup();
    let map = values
        .iter()
        .enumerate()
        .map(|(i, v)| (v.clone(), i as i64))
        .collect();
    (values, map)
}

fn i64_of(index: usize) -> Result<i64> {
    i64::try_from(index).map_err(|_| bad("PQP index exceeds 64 bits"))
}

/// `writePQPOutput_`.
fn write_pqp_output(path: &Path, exp: &TargetedExperiment) -> Result<()> {
    let db = ModificationsDB::global();
    let unimod = |peptide: &Peptide| -> Result<String> {
        get_aa_sequence_with(peptide, db, UnresolvedModificationPolicy::SourceSkip)?
            .sequence
            .to_unimod_string()
    };
    // Everything that can fail on the data is computed before the file is
    // touched, so a refused experiment leaves any existing file in place.
    let peptide_sequences = exp
        .peptides()
        .iter()
        .map(unimod)
        .collect::<Result<Vec<_>>>()?;
    let rows = exp
        .transitions()
        .iter()
        .map(|t| convert_transition(t, exp))
        .collect::<Result<Vec<_>>>()?;
    let peptide_rts = exp
        .peptides()
        .iter()
        .map(|p| p.retention_time())
        .collect::<Result<Vec<_>>>()?;
    let compound_rts = exp
        .compounds()
        .iter()
        .map(|c| c.retention_time())
        .collect::<Result<Vec<_>>>()?;

    let mut group_ids: Vec<String> = exp.peptides().iter().map(|p| p.id.clone()).collect();
    group_ids.extend(exp.compounds().iter().map(|c| c.id.clone()));
    let (compound_vec, compound_map) =
        ranked(exp.compounds().iter().map(|c| c.id.clone()).collect());
    let (protein_vec, protein_map) = ranked(exp.proteins().iter().map(|p| p.id.clone()).collect());
    let (_, group_map) = ranked(group_ids);
    let mut peptide_all = peptide_sequences.clone();
    for transition in exp.transitions() {
        let text = transition
            .cv_terms
            .metadata
            .get("Peptidoforms")
            .map(value_text)
            .unwrap_or_default();
        peptide_all.extend(split(&text, '|'));
    }
    let (peptide_vec, peptide_map) = ranked(peptide_all);
    let unmodified = peptide_vec
        .iter()
        .map(|s| AASequence::parse(s).map(|a| a.as_str().to_owned()))
        .collect::<Result<Vec<_>>>()?;

    // The source deletes an existing file first.
    match std::fs::remove_file(path) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
        Err(e) => return Err(e.into()),
    }
    let conn = SqliteConnector::new(path)?;
    let c = conn.connection();
    c.execute_batch(CREATE_SQL).map_err(sql)?;

    let group = |id: &str| group_map.get(id).copied().unwrap_or(0);
    let mut precursor_mz: HashMap<i64, f64> = HashMap::new();
    let mut precursor_decoy: HashMap<i64, bool> = HashMap::new();
    c.execute_batch("BEGIN TRANSACTION").map_err(sql)?;
    {
        let mut insert = c
            .prepare("INSERT INTO TRANSITION (ID, TRAML_ID, PRODUCT_MZ, CHARGE, TYPE, ANNOTATION, ORDINAL, DETECTING, IDENTIFYING, QUANTIFYING, LIBRARY_INTENSITY, DECOY) VALUES (?,?,?,?,?,?,?,?,?,?,?,?);")
            .map_err(sql)?;
        let mut peptide_mapping = c
            .prepare(
                "INSERT INTO TRANSITION_PEPTIDE_MAPPING (TRANSITION_ID, PEPTIDE_ID) VALUES (?,?);",
            )
            .map_err(sql)?;
        let mut precursor_mapping = c
            .prepare("INSERT INTO TRANSITION_PRECURSOR_MAPPING (TRANSITION_ID, PRECURSOR_ID) VALUES (?,?);")
            .map_err(sql)?;
        for (i, t) in rows.iter().enumerate() {
            let i = i64_of(i)?;
            let group_index = group(&t.group_id);
            precursor_mz.entry(group_index).or_insert(t.precursor);
            if !precursor_decoy.contains_key(&group_index) && t.detecting_transition {
                precursor_decoy.insert(group_index, t.decoy);
            }
            for form in &t.peptidoforms {
                peptide_mapping
                    .execute(params![i, peptide_map.get(form).copied().unwrap_or(0)])
                    .map_err(sql)?;
            }
            precursor_mapping
                .execute(params![i, group_index])
                .map_err(sql)?;
            // A missing charge ("NA") is NULL; any other text is bound as
            // text and converted by the INT column affinity, as upstream.
            let charge: Option<&str> =
                (t.fragment_charge != "NA").then_some(t.fragment_charge.as_str());
            insert
                .execute(params![
                    i,
                    t.transition_name,
                    t.product,
                    charge,
                    t.fragment_type,
                    t.annotation,
                    t.fragment_nr,
                    i32::from(t.detecting_transition),
                    i32::from(t.identifying_transition),
                    i32::from(t.quantifying_transition),
                    t.library_intensity,
                    i32::from(t.decoy),
                ])
                .map_err(sql)?;
        }
    }
    c.execute_batch("END TRANSACTION").map_err(sql)?;

    let mut peptide_protein: BTreeSet<(i64, i64)> = BTreeSet::new();
    let mut peptide_gene: BTreeSet<(i64, i64)> = BTreeSet::new();
    let mut gene_map: HashMap<String, i64> = HashMap::new();
    for (peptide, sequence) in exp.peptides().iter().zip(&peptide_sequences) {
        let peptide_index = peptide_map.get(sequence).copied().unwrap_or(0);
        for reference in &peptide.protein_refs {
            peptide_protein.insert((
                peptide_index,
                protein_map.get(reference).copied().unwrap_or(0),
            ));
        }
        let gene = peptide
            .cv_terms
            .metadata
            .get("GeneName")
            .map_or_else(|| "NA".to_owned(), value_text);
        let next = i64_of(gene_map.len())?;
        let gene_index = *gene_map.entry(gene).or_insert(next);
        peptide_gene.insert((peptide_index, gene_index));
    }
    let mut gene_vec: Vec<&String> = gene_map.keys().collect();
    gene_vec.sort();

    c.execute_batch("BEGIN TRANSACTION").map_err(sql)?;
    c.execute("INSERT INTO VERSION (ID) VALUES (3);", [])
        .map_err(sql)?;
    {
        let mut insert = c
            .prepare("INSERT INTO PROTEIN (ID, PROTEIN_ACCESSION, DECOY) VALUES (?,?,?);")
            .map_err(sql)?;
        for protein in &protein_vec {
            insert
                .execute(params![protein_map[protein], protein, 0])
                .map_err(sql)?;
        }
        let mut mapping = c
            .prepare("INSERT INTO PEPTIDE_PROTEIN_MAPPING (PEPTIDE_ID, PROTEIN_ID) VALUES (?,?);")
            .map_err(sql)?;
        for (peptide, protein) in &peptide_protein {
            mapping.execute(params![peptide, protein]).map_err(sql)?;
        }
        let mut insert = c
            .prepare("INSERT INTO GENE (ID, GENE_NAME, DECOY) VALUES (?,?,?);")
            .map_err(sql)?;
        for gene in &gene_vec {
            insert
                .execute(params![gene_map[*gene], gene, 0])
                .map_err(sql)?;
        }
        let mut mapping = c
            .prepare("INSERT INTO PEPTIDE_GENE_MAPPING (PEPTIDE_ID, GENE_ID) VALUES (?,?);")
            .map_err(sql)?;
        for (peptide, gene) in &peptide_gene {
            mapping.execute(params![peptide, gene]).map_err(sql)?;
        }
        let mut insert = c
            .prepare("INSERT INTO PEPTIDE (ID, UNMODIFIED_SEQUENCE, MODIFIED_SEQUENCE, DECOY) VALUES (?,?,?,?);")
            .map_err(sql)?;
        for (sequence, plain) in peptide_vec.iter().zip(&unmodified) {
            insert
                .execute(params![peptide_map[sequence], plain, sequence, 0])
                .map_err(sql)?;
        }
        let mut insert = c
            .prepare("INSERT INTO COMPOUND (ID, COMPOUND_NAME, SUM_FORMULA, SMILES, ADDUCTS, DECOY) VALUES (?,?,?,?,?,?);")
            .map_err(sql)?;
        for id in &compound_vec {
            let compound = exp.compound_by_ref(id)?;
            let meta = &compound.cv_terms.metadata;
            let adducts = meta.get("Adducts").map(value_text).unwrap_or_default();
            let name = meta
                .get("CompoundName")
                .map_or_else(|| compound.id.clone(), value_text);
            insert
                .execute(params![
                    compound_map[id],
                    name,
                    compound.molecular_formula,
                    compound.smiles_string,
                    adducts,
                    0
                ])
                .map_err(sql)?;
        }
        let mut mapping = c
            .prepare(
                "INSERT INTO PRECURSOR_PEPTIDE_MAPPING (PRECURSOR_ID, PEPTIDE_ID) VALUES (?,?);",
            )
            .map_err(sql)?;
        for (peptide, sequence) in exp.peptides().iter().zip(&peptide_sequences) {
            mapping
                .execute(params![
                    group(&peptide.id),
                    peptide_map.get(sequence).copied().unwrap_or(0)
                ])
                .map_err(sql)?;
        }
        let mut mapping = c
            .prepare(
                "INSERT INTO PRECURSOR_COMPOUND_MAPPING (PRECURSOR_ID, COMPOUND_ID) VALUES (?,?);",
            )
            .map_err(sql)?;
        for compound in exp.compounds() {
            mapping
                .execute(params![
                    group(&compound.id),
                    compound_map.get(&compound.id).copied().unwrap_or(0)
                ])
                .map_err(sql)?;
        }
        let mut insert = c
            .prepare("INSERT INTO PRECURSOR (ID, TRAML_ID, GROUP_LABEL, PRECURSOR_MZ, CHARGE, LIBRARY_INTENSITY, LIBRARY_DRIFT_TIME, LIBRARY_RT, DECOY) VALUES (?,?,?,?,?,NULL,?,?,?);")
            .map_err(sql)?;
        for (peptide, rt) in exp.peptides().iter().zip(&peptide_rts) {
            let index = group(&peptide.id);
            insert
                .execute(params![
                    index,
                    peptide.id,
                    peptide.peptide_group_label(),
                    precursor_mz.get(&index).copied().unwrap_or(0.0),
                    // The source reads an unset charge through a release-mode
                    // precondition, which yields its default 0.
                    peptide.charge_state().unwrap_or(0),
                    peptide.drift_time(),
                    rt,
                    i32::from(precursor_decoy.get(&index).copied().unwrap_or(false)),
                ])
                .map_err(sql)?;
        }
        for (compound, rt) in exp.compounds().iter().zip(&compound_rts) {
            let index = group(&compound.id);
            insert
                .execute(params![
                    index,
                    compound.id,
                    Option::<&str>::None,
                    precursor_mz.get(&index).copied().unwrap_or(0.0),
                    compound.charge_state(),
                    compound.drift_time(),
                    rt,
                    i32::from(precursor_decoy.get(&index).copied().unwrap_or(false)),
                ])
                .map_err(sql)?;
        }
    }
    c.execute_batch(UPDATE_DECOYS_SQL).map_err(sql)?;
    c.execute_batch("END TRANSACTION").map_err(sql)?;
    Ok(())
}
