// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! MzTab exporters (`FORMAT/MzTab.h` export surface) and the streaming
//! `MzTabFile::store` overloads that write them.
//!
//! Tier 1: the retained MzTabExporter outputs of the upstream TOPP suite
//! (test-data 0cb15f2 `topp/CMakeLists.txt:3382-3413`) are reproduced from the
//! retained inputs through the library, exactly as the tool (topp 174b576
//! `src/MzTabExporter.cpp`) calls it, and compared with the upstream `${DIFF}`
//! contract: `FuzzyDiff -ini FuzzyDiff.ini` (ratio 1.01, absdiff 0.01) with
//! each test's `-whitelist`. Tier 3: the two exporter sections of
//! `MzTab_test.cpp` (core bc9cc12, lines 153-193 and 208-303) with their
//! literals.

#![cfg(all(feature = "idxml", feature = "consensusxml", feature = "featurexml"))]

#[path = "support/fuzzy_string_comparator.rs"]
mod fuzzy;

use openms::analysis::id_filter::remove_empty_consensus_identifications;
use openms::chemistry::{ModificationsDB, ResidueModification};
use openms::format::mztab::MzTab;
use openms::format::mztab_export::{
    ConsensusMapExportOptions, ConsensusMapMzTabStream, IdMzTabStream, IdentificationExportOptions,
};
use openms::format::mztab_file::MzTabFile;
use openms::format::{consensusxml, featurexml, idxml};
use openms::identification::{
    PeptideIdentification, ProteinGroup, ProteinHit, ProteinIdentification,
};
use openms::kernel::{ColumnHeader, ConsensusFeature, ConsensusMap, DataArray, FeatureHandle};
use openms::metadata::MetaValue;
use std::collections::BTreeMap;
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/mztab_export")
        .join(name)
}

/// The upstream `${DIFF}` with a test's `-whitelist`. CMake turns the `\t` of
/// `"MTD\turi[1]"` into a tab.
fn assert_fuzzy_equal(produced: &str, retained: &str, whitelist: &[&str]) {
    let settings = fuzzy::FuzzyDiffSettings::upstream()
        .unwrap()
        .with_whitelist(whitelist);
    assert_eq!(settings.ratio, 1.01);
    assert_eq!(settings.absdiff, 0.01);
    let expected = std::fs::read(data(retained)).unwrap();
    if let Err(log) = settings.compare_bytes(produced.as_bytes(), &expected) {
        panic!("{retained}: FuzzyDiff failed:\n{log}\n--- produced ---\n{produced}");
    }
}

const WHITELIST: [&str; 3] = ["MTD\turi[1]", "MS:1000752", "MTD\tsoftware"];
const WHITELIST_PROTEIN_SCORE: [&str; 4] = [
    "MTD\turi[1]",
    "MS:1000752",
    "MTD\tsoftware",
    "protein_search_engine_score",
];

/// MzTabExporter's consensusXML branch: remove empty identifications, then
/// the streaming store with unidentified features and unassigned IDs
/// exported and the default `-opt_columns subfeatures`.
fn export_consensus(input: &str, first_run_inference_only: bool) -> String {
    let mut map = consensusxml::load(data(input)).unwrap();
    remove_empty_consensus_identifications(&mut map).unwrap();
    let options = ConsensusMapExportOptions {
        first_run_inference_only,
        export_unidentified_features: true,
        export_unassigned_ids: true,
        export_subfeatures: true,
        ..ConsensusMapExportOptions::default()
    };
    MzTabFile::new()
        .consensus_map_to_string(&map, "out.mzTab", &options)
        .unwrap()
}

/// A registry holding only the global registry's records of `names`.
///
/// The idXML reader charges every named sequence annotation as if it scanned
/// the whole registry (`AASequence::parse_with_budget`), so with the global
/// registry its fixed 50 M work budget runs out after about a dozen
/// modifications and `MzTabExporter_6_input.idXML` (17 annotations) is refused
/// with "peptide parent parsing resource limit exceeded". Reading through the
/// records the file names resolves every modification to the same record the
/// global registry would, which is all the exporter reads. The reader defect is
/// recorded in `docs/MZTAB_SUPPORT.md`.
fn registry_of(names: &[&str]) -> ModificationsDB {
    let records: Vec<ResidueModification> = ModificationsDB::global()
        .entries()
        .iter()
        .filter(|record| names.contains(&record.name()))
        .map(|record| ResidueModification::clone(record))
        .collect();
    ModificationsDB::from_records(records).unwrap()
}

fn read_ids(input: &str) -> idxml::IdXmlDocument {
    let registry = registry_of(&["Phospho", "Oxidation", "Carbamidomethyl"]);
    let file = std::io::BufReader::new(std::fs::File::open(data(input)).unwrap());
    idxml::read_with_registry(file, &idxml::ReadOptions::default(), &registry).unwrap()
}

/// MzTabExporter's idXML branch: the streaming identification store.
fn export_identifications(input: &str, export_all_psms: bool) -> String {
    let document = read_ids(input);
    let options = IdentificationExportOptions {
        export_all_psms,
        ..IdentificationExportOptions::default()
    };
    MzTabFile::new()
        .identifications_to_string(
            &document.protein_identifications,
            &document.peptide_identifications,
            "out.mzTab",
            &options,
        )
        .unwrap()
}

#[test]
fn retained_mztab_exporter_3_feature_map() {
    // TOPP_MzTabExporter_3: -in MzTabExporter_3_input.featureXML. The tool's
    // computeCoverage step only touches protein runs, which a feature map
    // export does not write.
    let map = featurexml::load(data("MzTabExporter_3_input.featureXML")).unwrap();
    let document = MzTab::export_feature_map(&map, "MzTabExporter_3_input.featureXML").unwrap();
    let text = MzTabFile::new().write_to_string(&document).unwrap();
    assert_fuzzy_equal(&text, "MzTabExporter_3_output.mzTab", &WHITELIST);
}

#[test]
fn retained_mztab_exporter_4_epifany_consensus_map() {
    // TOPP_MzTabExporter_4: -in Epifany_2_out.consensusXML
    let text = export_consensus("Epifany_2_out.consensusXML", false);
    assert_fuzzy_equal(
        &text,
        "MzTabExporter_4_output.mzTab",
        &WHITELIST_PROTEIN_SCORE,
    );
}

#[test]
fn retained_mztab_exporter_5_first_run_inference_only() {
    // TOPP_MzTabExporter_5: -in MzTabExporter_5_in.consensusXML -first_run_inference_only
    let text = export_consensus("MzTabExporter_5_in.consensusXML", true);
    assert_fuzzy_equal(
        &text,
        "MzTabExporter_5_output.mzTab",
        &WHITELIST_PROTEIN_SCORE,
    );
}

#[test]
fn retained_mztab_exporter_6_identifications() {
    // TOPP_MzTabExporter_6: -in MzTabExporter_6_input.idXML
    let text = export_identifications("MzTabExporter_6_input.idXML", false);
    assert_fuzzy_equal(
        &text,
        "MzTabExporter_6_output.mzTab",
        &WHITELIST_PROTEIN_SCORE,
    );
}

#[test]
fn retained_mztab_exporter_7_consensus_map() {
    // TOPP_MzTabExporter_7: -in MzTabExporter_7_input.consensusXML
    let text = export_consensus("MzTabExporter_7_input.consensusXML", false);
    assert_fuzzy_equal(
        &text,
        "MzTabExporter_7_output.mzTab",
        &WHITELIST_PROTEIN_SCORE,
    );
}

#[test]
fn retained_mztab_exporter_8_all_psms() {
    // TOPP_MzTabExporter_8: -in MzTabExporter_6_input.idXML -export_all_psms
    let text = export_identifications("MzTabExporter_6_input.idXML", true);
    assert_fuzzy_equal(
        &text,
        "MzTabExporter_8_output.mzTab",
        &WHITELIST_PROTEIN_SCORE,
    );
}

#[test]
fn stores_write_what_the_string_renderers_return() {
    let dir = std::env::temp_dir().join(format!("openms-mztab-export-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let document = read_ids("MzTabExporter_6_input.idXML");
    let options = IdentificationExportOptions::default();
    let path = dir.join("ids.mzTab");
    MzTabFile::new()
        .store_identifications(
            &path,
            &document.protein_identifications,
            &document.peptide_identifications,
            &options,
        )
        .unwrap();
    let written = std::fs::read_to_string(&path).unwrap();
    assert_eq!(
        written,
        export_identifications("MzTabExporter_6_input.idXML", false)
    );
    // The source refuses an extension that names another format.
    let refused = MzTabFile::new().store_identifications(
        dir.join("ids.idXML"),
        &document.protein_identifications,
        &document.peptide_identifications,
        &options,
    );
    assert!(refused.is_err());

    let mut map = consensusxml::load(data("MzTabExporter_7_input.consensusXML")).unwrap();
    remove_empty_consensus_identifications(&mut map).unwrap();
    let options = ConsensusMapExportOptions {
        export_unidentified_features: true,
        export_unassigned_ids: true,
        export_subfeatures: true,
        ..ConsensusMapExportOptions::default()
    };
    let path = dir.join("map.mzTab");
    MzTabFile::new()
        .store_consensus_map(&path, &map, &options)
        .unwrap();
    assert_eq!(
        std::fs::read_to_string(&path).unwrap(),
        export_consensus("MzTabExporter_7_input.consensusXML", false)
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn in_memory_exports_match_the_streams_row_for_row() {
    let document = read_ids("MzTabExporter_6_input.idXML");
    let options = IdentificationExportOptions::default();
    let exported = MzTab::export_identifications(
        &document.protein_identifications,
        &document.peptide_identifications,
        "in.idXML",
        &options,
    )
    .unwrap();
    let proteins: Vec<&ProteinIdentification> = document.protein_identifications.iter().collect();
    let peptides: Vec<&PeptideIdentification> = document.peptide_identifications.iter().collect();
    let mut stream = IdMzTabStream::new(&proteins, &peptides, "in.idXML", &options).unwrap();
    assert_eq!(stream.meta_data(), &exported.meta_data);
    assert!(stream.next_pep_row().is_none());
    let mut psm = Vec::new();
    while let Some(row) = stream.next_psm_row().unwrap() {
        psm.push(row);
    }
    // The in-memory identification export keeps placeholder rows, as the
    // source's exportIdentificationsToMzTab does.
    assert_eq!(psm, exported.psm_data);
    assert_eq!(
        exported.psm_data.len(),
        document.peptide_identifications.len()
    );

    let mut map = consensusxml::load(data("MzTabExporter_7_input.consensusXML")).unwrap();
    remove_empty_consensus_identifications(&mut map).unwrap();
    let options = ConsensusMapExportOptions {
        export_unidentified_features: true,
        export_unassigned_ids: true,
        export_subfeatures: true,
        ..ConsensusMapExportOptions::default()
    };
    let exported = MzTab::export_consensus_map(&map, "in.consensusXML", &options).unwrap();
    let mut stream = ConsensusMapMzTabStream::new(&map, "in.consensusXML", &options).unwrap();
    assert_eq!(stream.meta_data(), &exported.meta_data);
    let mut pep = Vec::new();
    while let Some(row) = stream.next_pep_row().unwrap() {
        pep.push(row);
    }
    assert_eq!(pep, exported.peptide_data);
    // The subfeature entries are in the rows but not in the stream's names.
    assert!(
        pep[0]
            .opt
            .iter()
            .any(|entry| entry.name.starts_with("opt_global_mass_to_charge_"))
    );
    assert!(
        !stream
            .peptide_optional_column_names()
            .iter()
            .any(|name| name.starts_with("opt_global_mass_to_charge_"))
    );
}

/// `[EXTRA] exportIdentificationsToMzTab terminates with export_all_psms on
/// empty-hit PeptideIdentification` (MzTab_test.cpp:153-193).
#[test]
fn export_all_psms_terminates_on_an_identification_without_hits() {
    let mut protein = ProteinIdentification::new();
    protein.identifier = "run_0".into();
    protein.primary_ms_run_paths = vec!["file.mzML".into()];
    let mut peptide = PeptideIdentification::new();
    peptide.identifier = "run_0".into();
    peptide.rt = Some(123.4);
    peptide.mz = Some(567.8);
    let options = IdentificationExportOptions {
        first_run_inference_only: false,
        export_empty_pep_ids: false,
        export_all_psms: true,
        title: "regression test".into(),
    };
    let document =
        MzTab::export_identifications(&[protein], &[peptide], "test.idXML", &options).unwrap();
    assert!(document.psm_data.len() <= 1);
    assert_eq!(document.meta_data.title.get(), "regression test");
}

fn header(path: &str, fraction_group: i64, fraction: i64) -> ColumnHeader {
    let mut header = ColumnHeader {
        filename: path.into(),
        label: "label-free".into(),
        ..ColumnHeader::default()
    };
    header
        .metadata
        .insert("fraction_group".into(), MetaValue::from(fraction_group));
    header
        .metadata
        .insert("fraction".into(), MetaValue::from(fraction));
    header
        .metadata
        .insert("sample_name".into(), MetaValue::from("sample_A".to_owned()));
    header
}

fn handle(map_index: u64, rt: f64, intensity: f32, unique_id: u64) -> FeatureHandle {
    FeatureHandle {
        map_index,
        unique_id,
        rt,
        mz: 500.0,
        intensity,
        charge: 0,
        width: 0.0,
    }
}

/// `[EXTRA] consensus-map assays use fraction-group/label grain and samples
/// are study variables` (MzTab_test.cpp:208-303).
#[test]
fn consensus_assays_use_fraction_group_label_grain_and_samples_are_study_variables() {
    let paths = ["sample_f1.mzML", "sample_f2.mzML", "sample_replicate.mzML"];
    let mut map = ConsensusMap::new();
    map.experiment_type = "label-free".into();
    let mut headers = BTreeMap::new();
    for (index, path) in paths.iter().enumerate() {
        let (group, fraction) = if index < 2 {
            (1, index as i64 + 1)
        } else {
            (2, 1)
        };
        headers.insert(index as u64, header(path, group, fraction));
    }
    map.column_headers = headers;
    let mut feature = ConsensusFeature::new();
    feature.insert(handle(0, 10.0, 2.0, 1)).unwrap();
    feature.insert(handle(1, 11.0, 3.0, 2)).unwrap();
    feature.insert(handle(2, 12.0, 7.0, 3)).unwrap();
    map.features.push(feature);

    let mut run = ProteinIdentification::new();
    run.identifier = "run_0".into();
    run.primary_ms_run_paths = paths.iter().map(|path| (*path).to_owned()).collect();
    run.hits.push(ProteinHit {
        accession: "P1".into(),
        ..ProteinHit::default()
    });
    // The assay arrays themselves are the quantified-group gate; the key order
    // is reversed to check matching through the parallel arrays.
    run.indistinguishable_groups.push(ProteinGroup {
        accessions: vec!["P1".into()],
        float_data_arrays: vec![DataArray::new(
            "fraction_group_level_abundance",
            vec![20.0_f32, 10.0],
        )],
        integer_data_arrays: vec![
            DataArray::new("fraction_group_level_fraction_group", vec![2, 1]),
            DataArray::new("fraction_group_level_label", vec![1, 1]),
        ],
        ..ProteinGroup::default()
    });
    map.protein_identifications.push(run);

    let options = ConsensusMapExportOptions {
        export_unidentified_features: true,
        ..ConsensusMapExportOptions::default()
    };
    let document = MzTab::export_consensus_map(&map, "test.consensusXML", &options).unwrap();

    let meta = &document.meta_data;
    assert_eq!(meta.assay.len(), 2);
    assert_eq!(meta.assay[&1].ms_run_ref, vec![1, 2]);
    assert_eq!(meta.assay[&2].ms_run_ref, vec![3]);
    assert_eq!(meta.study_variable.len(), 1);
    assert_eq!(meta.study_variable[&1].assay_refs, vec![1, 2]);

    let protein = document.protein_data.last().unwrap();
    assert_eq!(protein.protein_abundance_assay.len(), 2);
    assert!((protein.protein_abundance_assay[&1].get().unwrap() - 10.0).abs() < 1e-5);
    assert!((protein.protein_abundance_assay[&2].get().unwrap() - 20.0).abs() < 1e-5);
    assert!(protein.protein_abundance_study_variable[&1].is_null());

    let peptide = document.peptide_data.first().unwrap();
    assert_eq!(peptide.peptide_abundance_assay.len(), 2);
    assert!((peptide.peptide_abundance_assay[&1].get().unwrap() - 5.0).abs() < 1e-5);
    assert!((peptide.peptide_abundance_assay[&2].get().unwrap() - 7.0).abs() < 1e-5);
    assert!(peptide.peptide_abundance_study_variable[&1].is_null());
}

#[test]
fn modification_list_reports_terminal_and_residue_positions() {
    use openms::chemistry::AASequence;
    use openms::identification::PeptideHit;
    let sequence = AASequence::parse("(Acetyl)PEPS(Phospho)TIDE").unwrap();
    let mut hit = PeptideHit::new(1.0, 1, 2, sequence).unwrap();
    hit.metadata.insert(
        "Luciphor_global_flr".into(),
        MetaValue::try_from(0.25).unwrap(),
    );
    let list = MzTab::extract_modification_list(&hit, &[], &["Phospho (S)".into()]).unwrap();
    assert_eq!(
        list.to_cell_string().unwrap(),
        "0-UNIMOD:1,4[MS, MS:1002380, false localization rate, 0.25]-UNIMOD:21"
    );
    // The fixed-modification test compares the short identifier, so a full
    // identifier never filters (preserved source behaviour) and a short one
    // does.
    let unfiltered = MzTab::extract_modification_list(&hit, &["Phospho (S)".into()], &[]).unwrap();
    assert_eq!(unfiltered.get().len(), 2);
    let filtered = MzTab::extract_modification_list(&hit, &["Phospho".into()], &[]).unwrap();
    assert_eq!(filtered.to_cell_string().unwrap(), "0-UNIMOD:1");
}
