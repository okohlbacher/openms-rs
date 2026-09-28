// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! OMS files (`FORMAT/OMSFile.h`).
//!
//! Tier 1: the retained JSONExporter outputs of the upstream TOPP suite
//! (test-data 0cb15f2 `topp/CMakeLists.txt:3664-3672`) are reproduced from the
//! retained `.oms` inputs through `OmsFile::export_to_json`, which is all the
//! tool does (topp 174b576 `src/JSONExporter.cpp`), and compared with the
//! upstream `${DIFF}` contract (FuzzyDiff, ratio 1.01, absdiff 0.01, the INI's
//! whitelist).
//!
//! Tier 1 for load and store: the Release build (ibminode06, install
//! `openms4-release-bc9cc12-c19e494-174b576`, oracle `../oracle/F2-oms`, run
//! twice) loaded and re-stored the retained `.oms` inputs (schema versions 1
//! and 3) and the files it wrote itself from the `OMSFile_test.cpp` fixtures
//! (`IdXMLFile_whole.idXML` + `importIDs` + the Cl- adduct, sections 1-2;
//! `FeatureXMLFileOMStest_1.featureXML` + `importFeatureIDs`, sections 3-4;
//! `ConsensusXMLFile_1.consensusXML` + `importConsensusIDs`, sections 5-6).
//! `tests/data/oms_file/release/*.rt.tables.tsv` are canonical dumps of the
//! databases it wrote (`dump_tables.py` in the oracle directory); the port
//! must write databases with the same dump from the same inputs.
//!
//! The dump replaces every foreign key with the content of the row it
//! references, because the source orders observations, observation matches,
//! processing steps and parent groups by memory address and so renumbers
//! them from one round trip to the next (observed in the oracle; its own
//! second round trips agree with the first under this dump).

#![cfg(feature = "sqlite")]

#[path = "support/fuzzy_string_comparator.rs"]
mod fuzzy;

use openms::concept::progress_logger::ProgressLogType;
use openms::format::file_types::FileType;
use openms::format::oms_file::{OmsFile, OmsLoadOptions, OmsStoreOptions, SCHEMA_VERSION};
use rusqlite::Connection;
use rusqlite::types::ValueRef;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

fn release(name: &str) -> PathBuf {
    data("release").join(name)
}

fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("openms-oms-{tag}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Port of the oracle's `dump_tables.py`.
mod canonical {
    use super::*;

    const NUMERIC_KEYS: [&str; 4] = [
        "FEAT_BaseFeature",
        "FEAT_DataProcessing",
        "FEAT_ConsensusColumnHeader",
        "FEAT_MapMetaData",
    ];

    fn foreign(table: &str, column: &str) -> Option<String> {
        let fixed = match (table, column) {
            ("ID_ScoreType", "cv_term_id") => Some("CVTerm"),
            ("ID_ProcessingSoftware_AssignedScore", "software_id") => Some("ID_ProcessingSoftware"),
            ("ID_ProcessingSoftware_AssignedScore", "score_type_id") => Some("ID_ScoreType"),
            ("ID_ProcessingStep", "software_id") => Some("ID_ProcessingSoftware"),
            ("ID_ProcessingStep", "search_param_id") => Some("ID_DBSearchParam"),
            ("ID_ProcessingStep_InputFile", "processing_step_id") => Some("ID_ProcessingStep"),
            ("ID_ProcessingStep_InputFile", "input_file_id") => Some("ID_InputFile"),
            ("ID_Observation", "input_file_id") => Some("ID_InputFile"),
            ("ID_ParentGroup", "grouping_id") => Some("ID_ParentGroupSet"),
            ("ID_ParentGroup", "score_type_id") => Some("ID_ScoreType"),
            ("ID_ParentGroup_ParentSequence", "group_id") => Some("ID_ParentGroup"),
            ("ID_ParentGroup_ParentSequence", "parent_id") => Some("ID_ParentSequence"),
            ("ID_ParentMatch", "molecule_id") => Some("ID_IdentifiedMolecule"),
            ("ID_ParentMatch", "parent_id") => Some("ID_ParentSequence"),
            ("ID_IdentifiedCompound", "molecule_id") => Some("ID_IdentifiedMolecule"),
            ("ID_ObservationMatch", "identified_molecule_id") => Some("ID_IdentifiedMolecule"),
            ("ID_ObservationMatch", "observation_id") => Some("ID_Observation"),
            ("ID_ObservationMatch", "adduct_id") => Some("AdductInfo"),
            ("ID_ObservationMatch_PeakAnnotation", "parent_id") => Some("ID_ObservationMatch"),
            ("ID_ObservationMatch_PeakAnnotation", "processing_step_id") => {
                Some("ID_ProcessingStep")
            }
            ("FEAT_BaseFeature", "primary_molecule_id") => Some("ID_IdentifiedMolecule"),
            ("FEAT_ObservationMatch", "observation_match_id") => Some("ID_ObservationMatch"),
            _ => None,
        };
        if let Some(target) = fixed {
            return Some(target.to_owned());
        }
        for suffix in ["_MetaInfo", "_AppliedProcessingStep"] {
            if column == "parent_id" {
                if let Some(parent) = table.strip_suffix(suffix) {
                    return Some(parent.to_owned());
                }
            }
        }
        if table.ends_with("_AppliedProcessingStep") {
            if column == "processing_step_id" {
                return Some("ID_ProcessingStep".into());
            }
            if column == "score_type_id" {
                return Some("ID_ScoreType".into());
            }
        }
        None
    }

    fn esc(text: &str) -> String {
        text.replace('\\', "\\\\")
            .replace('\t', "\\t")
            .replace('\n', "\\n")
            .replace('\r', "\\r")
    }

    fn plain(value: ValueRef<'_>) -> String {
        match value {
            ValueRef::Null => "n".into(),
            ValueRef::Integer(number) => format!("i:{number}"),
            ValueRef::Real(number) => format!("f:{:016x}", number.to_bits()),
            ValueRef::Text(bytes) | ValueRef::Blob(bytes) => {
                format!("t:{}", esc(std::str::from_utf8(bytes).unwrap()))
            }
        }
    }

    struct Dump<'a> {
        db: &'a Connection,
        cache: BTreeMap<(String, i64), String>,
    }

    impl Dump<'_> {
        fn columns(&self, table: &str) -> Vec<String> {
            let statement = self
                .db
                .prepare(&format!("SELECT * FROM \"{table}\" LIMIT 0"))
                .unwrap();
            statement
                .column_names()
                .into_iter()
                .map(String::from)
                .collect()
        }

        fn cell(&mut self, table: &str, column: &str, value: ValueRef<'_>) -> String {
            match (foreign(table, column), value) {
                (Some(target), ValueRef::Integer(key))
                    if !NUMERIC_KEYS.contains(&target.as_str()) =>
                {
                    format!("r:{}", self.label(&target, key))
                }
                _ => plain(value),
            }
        }

        fn label(&mut self, table: &str, key: i64) -> String {
            if let Some(text) = self.cache.get(&(table.to_owned(), key)) {
                return text.clone();
            }
            let columns = self.columns(table);
            let rows = self.rows_where_id(table, key);
            let mut parts: Vec<String> = Vec::new();
            for row in rows {
                let cells: Vec<String> = columns
                    .iter()
                    .zip(row.iter())
                    .filter(|(column, _)| column.as_str() != "id")
                    .map(|(column, value)| self.cell(table, column, ValueRef::from(value)))
                    .collect();
                parts.push(cells.join(","));
            }
            parts.sort();
            let text = format!("{table}({})", parts.join("|"));
            self.cache.insert((table.to_owned(), key), text.clone());
            text
        }

        fn rows_where_id(&self, table: &str, key: i64) -> Vec<Vec<rusqlite::types::Value>> {
            self.all_rows(&format!("SELECT * FROM \"{table}\" WHERE id = {key}"))
        }

        fn all_rows(&self, sql: &str) -> Vec<Vec<rusqlite::types::Value>> {
            let mut statement = self.db.prepare(sql).unwrap();
            let count = statement.column_count();
            let mut rows = statement.query([]).unwrap();
            let mut out = Vec::new();
            while let Some(row) = rows.next().unwrap() {
                out.push((0..count).map(|i| row.get(i).unwrap()).collect());
            }
            out
        }
    }

    /// The canonical dump lines of an OMS file.
    pub fn dump(path: &Path) -> Vec<String> {
        let db = Connection::open(path).unwrap();
        let mut dump = Dump {
            db: &db,
            cache: BTreeMap::new(),
        };
        let mut schema: Vec<String> = dump
            .all_rows("SELECT type, name, sql FROM sqlite_master WHERE sql IS NOT NULL")
            .into_iter()
            .map(|row| {
                let text = |value: &rusqlite::types::Value| match value {
                    rusqlite::types::Value::Text(text) => text.clone(),
                    other => panic!("unexpected sqlite_master cell {other:?}"),
                };
                format!(
                    "schema\t{}\t{}\t{}",
                    text(&row[0]),
                    text(&row[1]),
                    esc(&text(&row[2]))
                )
            })
            .collect();
        schema.sort();
        let tables: Vec<String> = dump
            .all_rows(
                "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'",
            )
            .into_iter()
            .map(|row| match &row[0] {
                rusqlite::types::Value::Text(text) => text.clone(),
                other => panic!("unexpected table name {other:?}"),
            })
            .collect();
        let mut rows = Vec::new();
        for table in tables {
            if table == "version" {
                continue;
            }
            let columns = dump.columns(&table);
            let drop_id = !NUMERIC_KEYS.contains(&table.as_str());
            for row in dump.all_rows(&format!("SELECT * FROM \"{table}\"")) {
                let cells: Vec<String> = columns
                    .iter()
                    .zip(row.iter())
                    .filter(|(column, _)| !(drop_id && column.as_str() == "id"))
                    .map(|(column, value)| {
                        format!(
                            "{column}={}",
                            dump.cell(&table, column, ValueRef::from(value))
                        )
                    })
                    .collect();
                rows.push(format!("row\t{table}\t{}", cells.join("\t")));
            }
        }
        rows.sort();
        schema.extend(rows);
        schema
    }
}

/// Compare the canonical dump of `produced` with a retained Release dump.
fn assert_same_database(produced: &Path, expected: &str) {
    let expected_text = std::fs::read_to_string(release(expected)).unwrap();
    let expected_lines: Vec<&str> = expected_text.lines().collect();
    let produced_lines = canonical::dump(produced);
    let mut first_difference = None;
    for (index, (a, b)) in produced_lines.iter().zip(expected_lines.iter()).enumerate() {
        if a != b {
            first_difference = Some(format!(
                "line {}:\n  port:    {a}\n  release: {b}",
                index + 1
            ));
            break;
        }
    }
    if first_difference.is_none() && produced_lines.len() != expected_lines.len() {
        first_difference = Some(format!(
            "line counts differ: port {}, release {}",
            produced_lines.len(),
            expected_lines.len()
        ));
    }
    if let Some(difference) = first_difference {
        panic!("{expected}: the port's database differs from the Release build's\n{difference}");
    }
}

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/oms_file")
        .join(name)
}

fn assert_fuzzy_equal(produced: &str, retained: &str) {
    let settings = fuzzy::FuzzyDiffSettings::upstream().unwrap();
    assert_eq!(settings.ratio, 1.01);
    assert_eq!(settings.absdiff, 0.01);
    let expected = std::fs::read(data(retained)).unwrap();
    if let Err(log) = settings.compare_bytes(produced.as_bytes(), &expected) {
        panic!("{retained}: FuzzyDiff failed:\n{log}");
    }
}

#[test]
fn retained_json_exporter_1_protein() {
    // TOPP_JSONExporter_1: -in JSONExporter_protein.oms
    let text = OmsFile::default()
        .json_text(data("JSONExporter_protein.oms"))
        .unwrap();
    assert_fuzzy_equal(&text, "JSONExporter_1_output.json");
}

#[test]
fn retained_json_exporter_2_rna() {
    // TOPP_JSONExporter_2: -in JSONExporter_RNA.oms
    let text = OmsFile::default()
        .json_text(data("JSONExporter_RNA.oms"))
        .unwrap();
    assert_fuzzy_equal(&text, "JSONExporter_2_output.json");
}

#[test]
fn export_to_json_writes_the_rendered_text() {
    let dir = std::env::temp_dir().join(format!("openms-oms-json-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let out = dir.join("protein.json");
    let adapter = OmsFile::new(ProgressLogType::None);
    adapter
        .export_to_json(data("JSONExporter_protein.oms"), &out)
        .unwrap();
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        adapter.json_text(data("JSONExporter_protein.oms")).unwrap()
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn retained_inputs_declare_their_schema_version() {
    // The retained JSONExporter outputs record `"OMSFile": 3` in their
    // `version` row (test-data 0cb15f2 topp/JSONExporter_1_output.json and
    // JSONExporter_2_output.json): both inputs predate the pinned writer's
    // schema version 5.
    let adapter = OmsFile::default();
    assert_eq!(
        adapter
            .schema_version(data("JSONExporter_protein.oms"))
            .unwrap(),
        3
    );
    assert_eq!(
        adapter
            .schema_version(data("JSONExporter_RNA.oms"))
            .unwrap(),
        3
    );
    assert_eq!(SCHEMA_VERSION, 5);
}

#[test]
fn a_file_that_is_not_an_oms_database_is_refused() {
    let adapter = OmsFile::default();
    // Not an SQLite database at all.
    assert!(
        adapter
            .json_text(data("JSONExporter_1_output.json"))
            .is_err()
    );
    // A missing file is not created by the read-only open.
    let missing = data("does_not_exist.oms");
    assert!(adapter.json_text(&missing).is_err());
    assert!(!missing.exists());
    assert!(adapter.load_identifications(&missing).is_err());
    assert!(!missing.exists());
}

#[test]
fn identification_round_trips_match_the_release_build() {
    // Retained inputs (schema versions 3 and 1) and the file the Release build
    // wrote for OMSFile_test sections 1-2 (version 5).
    let cases = [
        (
            data("JSONExporter_protein.oms"),
            "id_JSONExporter_protein.rt.tables.tsv",
        ),
        (
            data("JSONExporter_RNA.oms"),
            "id_JSONExporter_RNA.rt.tables.tsv",
        ),
        (
            data("IDMerger_6_input1.oms"),
            "id_IDMerger_6_input1.rt.tables.tsv",
        ),
        (
            data("IDMerger_6_input2.oms"),
            "id_IDMerger_6_input2.rt.tables.tsv",
        ),
        (release("idxml_whole.oms"), "idxml_whole.rt.tables.tsv"),
    ];
    let dir = temp_dir("id-round-trip");
    let adapter = OmsFile::default();
    for (input, expected) in cases {
        let ids = adapter
            .load_identifications_with(&input, OmsLoadOptions::source())
            .unwrap();
        let out = dir.join(expected.replace(".tables.tsv", ".oms"));
        adapter
            .store_identifications_with(&out, &ids, OmsStoreOptions::source())
            .unwrap();
        assert_same_database(&out, expected);
        assert_eq!(adapter.schema_version(&out).unwrap(), SCHEMA_VERSION);
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn feature_map_round_trip_matches_the_release_build() {
    let dir = temp_dir("feature-round-trip");
    let adapter = OmsFile::default();
    let loaded = adapter
        .load_feature_map_with(release("feature.oms"), OmsLoadOptions::source())
        .unwrap();
    assert!(loaded.skipped_actions.is_empty());
    let out = dir.join("feature.rt.oms");
    adapter
        .store_feature_map_with(
            &out,
            &loaded.map,
            Some(&loaded.identifications),
            OmsStoreOptions::source(),
        )
        .unwrap();
    assert_same_database(&out, "feature.rt.tables.tsv");
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn consensus_map_round_trip_matches_the_release_build() {
    let dir = temp_dir("consensus-round-trip");
    let adapter = OmsFile::default();
    let loaded = adapter
        .load_consensus_map_with(release("consensus.oms"), OmsLoadOptions::source())
        .unwrap();
    let out = dir.join("consensus.rt.oms");
    adapter
        .store_consensus_map_with(
            &out,
            &loaded.map,
            Some(&loaded.identifications),
            OmsStoreOptions::source(),
        )
        .unwrap();
    assert_same_database(&out, "consensus.rt.tables.tsv");
    std::fs::remove_dir_all(&dir).unwrap();
}

/// Fields of one tab-separated line of the Release build's loaded-map dump.
fn fields(line: &str) -> Vec<&str> {
    line.split('\t').collect()
}

fn number(text: &str) -> f64 {
    text.parse().unwrap()
}

fn flatten<'a>(
    features: &'a [openms::kernel::Feature],
    depth: usize,
    out: &mut Vec<(usize, &'a openms::kernel::Feature)>,
) {
    for feature in features {
        out.push((depth, feature));
        flatten(&feature.subordinates, depth + 1, out);
    }
}

#[test]
fn loaded_feature_map_matches_the_release_build() {
    // OMSFile_test.cpp section 4: two features, the first with two
    // subordinates; the Release build's loaded values (feature_loaded.tsv,
    // printed with 17 significant digits).
    let loaded = OmsFile::default()
        .load_feature_map_with(release("feature.oms"), OmsLoadOptions::source())
        .unwrap();
    let map = &loaded.map;
    assert_eq!(map.features.len(), 2);
    assert_eq!(map.features[0].subordinates.len(), 2);
    let expected = std::fs::read_to_string(release("feature_loaded.tsv")).unwrap();
    let mut lines = expected.lines();
    let head = fields(lines.next().unwrap());
    assert_eq!(map.unique_id.to_string(), head[1]);
    assert_eq!(map.identifier, head[2]);
    assert_eq!(map.loaded_file_path, head[3]);
    assert_eq!(map.data_processing.len().to_string(), head[4]);
    assert_eq!(map.loaded_file_type, FileType::Unknown);
    let mut flat = Vec::new();
    flatten(&map.features, 0, &mut flat);
    let rows: Vec<&str> = lines.collect();
    assert_eq!(flat.len(), rows.len());
    for ((depth, feature), row) in flat.into_iter().zip(rows) {
        let f = fields(row);
        assert_eq!(depth.to_string(), f[0]);
        assert_eq!(feature.unique_id.to_string(), f[1]);
        assert_eq!(feature.rt, number(f[2]));
        assert_eq!(feature.mz, number(f[3]));
        assert_eq!(f64::from(feature.intensity), number(f[4]));
        assert_eq!(feature.charge.to_string(), f[5]);
        assert_eq!(f64::from(feature.width), number(f[6]));
        assert_eq!(f64::from(feature.quality), number(f[7]));
        assert_eq!(f64::from(feature.quality_rt), number(f[8]));
        assert_eq!(f64::from(feature.quality_mz), number(f[9]));
        assert_eq!(feature.convex_hulls.len().to_string(), f[10]);
        assert_eq!(feature.id_matches.len().to_string(), f[11]);
        assert_eq!(u8::from(feature.primary_id.is_some()).to_string(), f[12]);
        assert_eq!(feature.subordinates.len().to_string(), f[13]);
        // Every observation match a feature references is in the graph.
        for &id in &feature.id_matches {
            loaded.identifications.observation_match(id).unwrap();
        }
    }
}

#[test]
fn loaded_consensus_map_matches_the_release_build() {
    // OMSFile_test.cpp section 6: six consensus features, the first with one
    // handle and the second with two.
    let loaded = OmsFile::default()
        .load_consensus_map_with(release("consensus.oms"), OmsLoadOptions::source())
        .unwrap();
    let map = &loaded.map;
    assert_eq!(map.features.len(), 6);
    assert_eq!(map.features[0].handles().len(), 1);
    assert_eq!(map.features[1].handles().len(), 2);
    let expected = std::fs::read_to_string(release("consensus_loaded.tsv")).unwrap();
    let mut features = map.features.iter();
    let mut handles: std::slice::Iter<'_, openms::kernel::FeatureHandle> = [].iter();
    let mut header_count = 0;
    for line in expected.lines() {
        let f = fields(line);
        match f[0] {
            "map" => {
                assert_eq!(map.unique_id.to_string(), f[1]);
                assert_eq!(map.experiment_type, f[2]);
                assert_eq!(map.column_headers.len().to_string(), f[3]);
                assert_eq!(map.data_processing.len().to_string(), f[4]);
            }
            "header" => {
                header_count += 1;
                let header = &map.column_headers[&f[1].parse::<u64>().unwrap()];
                assert_eq!(header.filename, f[2]);
                assert_eq!(header.label, f[3]);
                assert_eq!(header.size.to_string(), f[4]);
                assert_eq!(header.unique_id.to_string(), f[5]);
            }
            "feature" => {
                assert!(handles.next().is_none(), "unread handles remain");
                let feature = features.next().unwrap();
                assert_eq!(feature.unique_id.to_string(), f[1]);
                assert_eq!(feature.rt, number(f[2]));
                assert_eq!(feature.mz, number(f[3]));
                assert_eq!(f64::from(feature.intensity), number(f[4]));
                assert_eq!(feature.charge.to_string(), f[5]);
                assert_eq!(f64::from(feature.width), number(f[6]));
                assert_eq!(f64::from(feature.quality), number(f[7]));
                assert_eq!(feature.id_matches.len().to_string(), f[8]);
                assert_eq!(u8::from(feature.primary_id.is_some()).to_string(), f[9]);
                assert_eq!(feature.ratios.len().to_string(), f[10]);
                handles = feature.handles().iter();
            }
            "handle" => {
                let handle = handles.next().unwrap();
                assert_eq!(handle.map_index.to_string(), f[1]);
                assert_eq!(handle.unique_id.to_string(), f[2]);
                assert_eq!(handle.rt, number(f[3]));
                assert_eq!(handle.mz, number(f[4]));
                assert_eq!(f64::from(handle.intensity), number(f[5]));
                assert_eq!(handle.charge.to_string(), f[6]);
                assert_eq!(f64::from(handle.width), number(f[7]));
            }
            other => panic!("unexpected dump line kind {other}"),
        }
    }
    assert!(features.next().is_none());
    assert_eq!(header_count, map.column_headers.len());
}

#[test]
fn identification_counts_survive_a_round_trip() {
    // OMSFile_test.cpp section 2: every record count, the number of steps and
    // scores per observation match, and the observation match with the Cl-
    // adduct survive store and load.
    let dir = temp_dir("id-counts");
    let adapter = OmsFile::default();
    let ids = adapter
        .load_identifications(release("idxml_whole.oms"))
        .unwrap();
    let out = dir.join("idxml_whole.oms");
    adapter.store_identifications(&out, &ids).unwrap();
    let back = adapter.load_identifications(&out).unwrap();
    // The Release build's counts after importIDs plus the adduct, and after
    // its own reload (oracle counts.tsv): 2 input files, 2 score types,
    // 2 softwares, 2 search parameter sets, 2 steps, 3 observations,
    // 3 parents, 2 group sets, 5 peptides, 0 oligos, 0 compounds, 1 adduct,
    // 6 observation matches.
    let counts = |d: &openms::identification::graph::IdentificationData| {
        [
            d.input_file_count(),
            d.score_type_count(),
            d.processing_software_count(),
            d.db_search_param_count(),
            d.processing_step_count(),
            d.observation_count(),
            d.parent_count(),
            d.parent_group_set_count(),
            d.peptide_count(),
            d.oligo_count(),
            d.compound_count(),
            d.adduct_count(),
            d.observation_match_count(),
        ]
    };
    assert_eq!(counts(&ids), [2, 2, 2, 2, 2, 3, 3, 2, 5, 0, 0, 1, 6]);
    assert_eq!(counts(&back), counts(&ids));
    let steps = |d: &openms::identification::graph::IdentificationData| -> Vec<usize> {
        let mut sizes: Vec<usize> = d
            .observation_matches()
            .map(|(_, item)| item.result.steps_and_scores.len())
            .collect();
        sizes.sort_unstable();
        sizes
    };
    assert_eq!(steps(&ids), steps(&back));
    let with_adduct: Vec<_> = back
        .observation_matches()
        .filter(|(_, item)| item.adduct.is_some())
        .collect();
    assert_eq!(with_adduct.len(), 1);
    let adduct = back.adduct(with_adduct[0].1.adduct.unwrap()).unwrap();
    assert_eq!(adduct.name(), "Cl-");
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn the_json_export_of_a_feature_map_fails_as_in_the_source() {
    // The Release build's exportToJSON orders FEAT_Feature by a column "id"
    // it does not have ("no such column: id", oracle counts.tsv); the export
    // is ported as it is.
    let error = OmsFile::default()
        .json_text(release("feature.oms"))
        .unwrap_err();
    assert!(error.to_string().contains("no such column: id"), "{error}");
}

#[test]
fn native_options_keep_what_the_source_loses() {
    use openms::kernel::features::Ratio;
    use openms::kernel::{BaseFeature, ConsensusFeature, ConsensusMap, FeatureHandle};
    let dir = temp_dir("native");
    let adapter = OmsFile::default();
    let mut map = ConsensusMap {
        loaded_file_type: FileType::ConsensusXml,
        ..ConsensusMap::default()
    };
    map.metadata.insert(
        "list".into(),
        openms::metadata::MetaValue::from(vec!["a".to_owned(), "b".to_owned()]),
    );
    let mut feature = ConsensusFeature::new();
    feature.base = BaseFeature::new(10.0, 500.0, 100.0);
    let mut first = BaseFeature::new(10.0, 500.0, 60.0);
    first.unique_id = 17;
    let mut second = BaseFeature::new(10.5, 500.1, 40.0);
    second.unique_id = 18;
    feature.insert(FeatureHandle::new(0, &first)).unwrap();
    feature.insert(FeatureHandle::new(1, &second)).unwrap();
    feature.ratios.push(Ratio {
        ratio_value: 1.5,
        denominator_ref: "0".into(),
        numerator_ref: "1".into(),
        description: vec!["light".into(), "heavy".into()],
    });
    map.features.push(feature);

    // Native store and load keep ratios, handle IDs, the file type and lists.
    let native = dir.join("native.oms");
    adapter.store_consensus_map(&native, &map, None).unwrap();
    let back = adapter.load_consensus_map(&native).unwrap().map;
    assert_eq!(back.loaded_file_type, FileType::ConsensusXml);
    assert_eq!(back.metadata["list"], map.metadata["list"]);
    assert_eq!(back.features[0].ratios, map.features[0].ratios);
    let ids: Vec<u64> = back.features[0]
        .handles()
        .iter()
        .map(|h| h.unique_id)
        .collect();
    assert_eq!(ids, [17, 18]);

    // The source's behaviour, selected explicitly, loses each of them.
    let source = dir.join("source.oms");
    adapter
        .store_consensus_map_with(&source, &map, None, OmsStoreOptions::source())
        .unwrap();
    let lost = adapter
        .load_consensus_map_with(&source, OmsLoadOptions::source())
        .unwrap()
        .map;
    assert_eq!(lost.loaded_file_type, FileType::Unknown);
    assert_eq!(
        lost.metadata["list"],
        openms::metadata::MetaValue::from(vec!["a".to_owned(), " b".to_owned()])
    );
    assert!(lost.features[0].ratios.is_empty());
    let ids: Vec<u64> = lost.features[0]
        .handles()
        .iter()
        .map(|h| h.unique_id)
        .collect();
    assert_eq!(ids, [0, 0]);
    // A native load finds the ratios where the source's writer put them.
    let found = adapter.load_consensus_map(&source).unwrap().map;
    assert_eq!(found.features[0].ratios, map.features[0].ratios);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_failed_store_leaves_the_existing_file_unchanged() {
    let dir = temp_dir("atomic");
    let adapter = OmsFile::default();
    let out = dir.join("existing.oms");
    std::fs::write(&out, b"previous content").unwrap();
    // Features referencing observation matches of a graph that is not passed
    // to the store.
    let loaded = adapter.load_feature_map(release("feature.oms")).unwrap();
    let mut flat = Vec::new();
    flatten(&loaded.map.features, 0, &mut flat);
    assert!(flat.iter().any(|(_, f)| !f.id_matches.is_empty()));
    assert!(adapter.store_feature_map(&out, &loaded.map, None).is_err());
    assert_eq!(std::fs::read(&out).unwrap(), b"previous content");
    let entries = std::fs::read_dir(&dir).unwrap().count();
    assert_eq!(entries, 1, "the temporary database was removed");
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_map_without_identifications_is_stored_where_the_source_fails() {
    // The Release build cannot store a feature or consensus map without
    // identification data: FEAT_BaseFeature declares a foreign key into the
    // missing ID_IdentifiedMolecule table (oracle F2-oms-noid, run twice on
    // ibminode06, noid_outcome.tsv). The port declares that key only when the
    // table exists.
    let outcome = std::fs::read_to_string(release("noid_outcome.tsv")).unwrap();
    assert_eq!(
        outcome
            .matches("no such table: main.ID_IdentifiedMolecule")
            .count(),
        2
    );
    let dir = temp_dir("noid");
    let adapter = OmsFile::default();
    let loaded = adapter.load_feature_map(release("feature.oms")).unwrap();
    let mut map = loaded.map.clone();
    fn clear(features: &mut [openms::kernel::Feature]) {
        for feature in features {
            feature.id_matches.clear();
            feature.primary_id = None;
            clear(&mut feature.subordinates);
        }
    }
    clear(&mut map.features);
    let out = dir.join("noid.oms");
    adapter.store_feature_map(&out, &map, None).unwrap();
    let back = adapter.load_feature_map(&out).unwrap();
    assert!(back.identifications.is_empty());
    assert_eq!(back.map.features.len(), map.features.len());
    assert_eq!(back.map.features[0].subordinates.len(), 2);
    assert_eq!(back.map.loaded_file_type, map.loaded_file_type);
    std::fs::remove_dir_all(&dir).unwrap();
}
