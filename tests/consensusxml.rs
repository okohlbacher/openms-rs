// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "consensusxml")]

use openms::chemistry::ModificationsDB;
use openms::format::consensusxml::{self, ReadOptions, WriteOptions};
use openms::identification::ProteinGroup;
use openms::kernel::ConsensusMap;
use openms::kernel::{ConsensusFeature, DataArray, FeatureHandle};
use openms::metadata::{MetaValue, ProcessingAction, Unit};
use std::io::Cursor;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
const SOURCE: &[u8] = include_bytes!("data/consensusxml/ConsensusXMLFile_1.consensusXML");
const OPTIONS: &[u8] = include_bytes!("data/consensusxml/ConsensusXMLFile_2_options.consensusXML");
fn source() -> ConsensusMap {
    consensusxml::read(Cursor::new(SOURCE)).unwrap()
}
fn roundtrip(map: &ConsensusMap) -> ConsensusMap {
    let mut bytes = Vec::new();
    consensusxml::write(&mut bytes, map).unwrap();
    consensusxml::read(Cursor::new(bytes)).unwrap()
}
struct Directory(PathBuf);
impl Directory {
    fn new() -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let p = std::env::temp_dir().join(format!(
            "openms-consensus-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&p).unwrap();
        Self(p)
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
fn pinned_source_geometry_metadata_processing_and_identifications() {
    let map = source();
    assert_eq!(map.len(), 6);
    assert_eq!(map.identifier, "lsid");
    assert_eq!(map.experiment_type, "label-free");
    assert_eq!(map.metadata["name2"].as_i64().unwrap(), 2);
    assert_eq!(map.column_headers[&0].size, 144);
    assert_eq!(
        map.column_headers[&1].metadata["name6"].as_f64().unwrap(),
        6.0
    );
    assert_eq!(map.data_processing.len(), 2);
    assert_eq!(map.data_processing[0].software.name, "Software1");
    assert!(
        map.data_processing[1]
            .actions
            .contains(&ProcessingAction::BaselineReduction)
    );
    assert_eq!(map.protein_identifications[0].hits.len(), 2);
    assert_eq!(map.protein_identifications[1].hits[0].sequence, "OPQREST");
    assert_eq!(
        map.unassigned_peptide_identifications[1].hits[1]
            .sequence
            .to_string(),
        "H"
    );
    assert_eq!(
        map.features[0].peptide_identifications[1].hits[0]
            .sequence
            .to_string(),
        "C"
    );
    assert_eq!(map.features[0].rt, 1273.27);
    assert_eq!(map.features[0].mz, 904.47);
    assert_eq!(map.features[0].intensity, 31_253_900.0);
    assert_eq!(map.features[5].rt, 1194.82);
    assert_eq!(map.features[5].mz, 777.101);
    assert_eq!(map.features[5].intensity, 17_821_500.0);
    assert_eq!(
        map.features[0].metadata["myIntList"]
            .as_integer_list()
            .unwrap(),
        [1, 10, 12]
    );
    assert_eq!(
        map.features[0].metadata["myStringList"]
            .as_string_list()
            .unwrap(),
        ["myABC1", "Stuff", "12"]
    );
    assert_eq!(roundtrip(&map), map);
}
#[test]
fn independent_source_option_examples_and_half_open_endpoints() {
    let mut options = ReadOptions {
        rt_range: Some(815.0..818.0),
        ..Default::default()
    };
    let map = consensusxml::read_with_options(Cursor::new(OPTIONS), &options).unwrap();
    assert_eq!(map.len(), 1);
    assert_eq!(map.features[0].rt, 817.266);
    options.rt_range = None;
    options.mz_range = Some(944.0..945.0);
    let map = consensusxml::read_with_options(Cursor::new(OPTIONS), &options).unwrap();
    assert_eq!(map.len(), 1);
    assert_eq!(map.features[0].mz, 944.96);
    options.mz_range = None;
    options.intensity_range = Some(15000.0..24000.0);
    let map = consensusxml::read_with_options(Cursor::new(OPTIONS), &options).unwrap();
    assert_eq!(map.len(), 1);
    assert_eq!(map.features[0].intensity, 23000.238);
    options = ReadOptions {
        rt_range: Some(1273.27..1273.28),
        ..Default::default()
    };
    assert_eq!(
        consensusxml::read_with_options(Cursor::new(SOURCE), &options)
            .unwrap()
            .len(),
        1
    );
    options.rt_range = Some(1273.0..1273.27);
    assert!(
        consensusxml::read_with_options(Cursor::new(SOURCE), &options)
            .unwrap()
            .is_empty()
    );
    options.rt_range = Some(1.0..0.0);
    assert!(consensusxml::read_with_options(Cursor::new(SOURCE), &options).is_err());
}
#[test]
fn compressed_path_roundtrips_preserve_payload_and_atomic_failure() {
    let dir = Directory::new();
    let map = source();
    for suffix in ["consensusXML", "consensusXML.gz", "consensusXML.bz2"] {
        let path = dir.0.join(format!("map.{suffix}"));
        consensusxml::store(&path, &map).unwrap();
        let raw = std::fs::read(&path).unwrap();
        if suffix.ends_with("gz") {
            assert!(raw.starts_with(&[31, 139]));
        }
        if suffix.ends_with("bz2") {
            assert!(raw.starts_with(b"BZh"));
        }
        let mut copy = consensusxml::load(&path).unwrap();
        assert_eq!(copy.loaded_file_path, path.to_str().unwrap());
        assert_eq!(
            copy.loaded_file_type,
            openms::format::FileType::ConsensusXml
        );
        copy.loaded_file_path.clear();
        copy.loaded_file_type = openms::format::FileType::Unknown;
        assert_eq!(copy, map);
        let mut bad = map.clone();
        bad.features[0].width = 1.0;
        assert!(consensusxml::store(&path, &bad).is_err());
        assert_eq!(std::fs::read(&path).unwrap(), raw);
        assert!(
            consensusxml::load_with_options(
                &path,
                &ReadOptions {
                    max_xml_bytes: 32,
                    ..Default::default()
                }
            )
            .is_err()
        );
    }
    assert_eq!(std::fs::read_dir(&dir.0).unwrap().count(), 3);
}
fn quantified_map() -> ConsensusMap {
    let mut map = source();
    let run = &mut map.protein_identifications[0];
    run.indistinguishable_groups.push(ProteinGroup {
        probability: 0.98,
        accessions: run.hits.iter().map(|h| h.accession.clone()).collect(),
        float_data_arrays: vec![
            DataArray::new("fraction_group_level_abundance", vec![0.0, 4.25]),
            DataArray::new("custom_float", vec![2.5]),
        ],
        integer_data_arrays: vec![DataArray::new(
            "fraction_group_level_fraction_group",
            vec![1, 2],
        )],
        string_data_arrays: vec![DataArray::new(
            "fraction_group_level_label",
            vec!["light".into(), "heavy".into()],
        )],
    });
    map
}
#[test]
fn quantitative_arrays_use_canonical_slots_and_stable_accession_guard() {
    let map = quantified_map();
    let copy = roundtrip(&map);
    let group = &copy.protein_identifications[0].indistinguishable_groups[0];
    assert_eq!(
        group
            .float_data_arrays
            .iter()
            .map(|a| a.name.as_str())
            .collect::<Vec<_>>(),
        [
            "psm_count",
            "distinct_peptides",
            "file_channel_level_abundance",
            "custom_float",
            "fraction_group_level_abundance"
        ]
    );
    assert_eq!(group.float_data_arrays[4].data, [0.0, 4.25]);
    assert_eq!(group.integer_data_arrays[2].data, [1, 2]);
    assert_eq!(group.string_data_arrays[2].data, ["light", "heavy"]);
    assert_eq!(roundtrip(&copy), copy);
    assert!(
        !copy.protein_identifications[0]
            .metadata
            .keys()
            .any(|k| k.starts_with("indistinguishable_proteins_"))
    );
}
#[test]
fn stale_group_ownership_is_checked_or_reported_and_discarded() {
    let map = quantified_map();
    let mut bytes = Vec::new();
    consensusxml::write(&mut bytes, &map).unwrap();
    let text = String::from_utf8(bytes)
        .unwrap()
        .replace("value=\"0.98,PH_0,PH_1\"", "value=\"0.98,PH_1\"");
    assert!(consensusxml::read(Cursor::new(text.as_bytes())).is_err());
    let report = consensusxml::read_report(
        Cursor::new(text.as_bytes()),
        &ReadOptions {
            discard_mismatched_quantities: true,
            ..Default::default()
        },
        ModificationsDB::global(),
    )
    .unwrap();
    assert!(report.warnings.iter().any(|s| s.contains("ownership")));
    let group = &report.map.protein_identifications[0].indistinguishable_groups[0];
    assert_eq!(group.accessions.len(), 1);
    assert!(group.float_data_arrays.is_empty());
}
#[test]
fn legacy_abundance_anchor_restores_two_zero_count_arrays() {
    let mut map = quantified_map();
    let group = &mut map.protein_identifications[0].indistinguishable_groups[0];
    group.float_data_arrays = vec![DataArray::new("abundances", vec![2.0, 3.0, 0.0])];
    group.integer_data_arrays.clear();
    group.string_data_arrays.clear();
    let copy = roundtrip(&map);
    let arrays = &copy.protein_identifications[0].indistinguishable_groups[0].float_data_arrays;
    assert_eq!(arrays[0].name, "abundances");
    assert_eq!(arrays[1].data, [0.0, 0.0, 0.0]);
    assert_eq!(arrays[2].data, [0.0, 0.0, 0.0]);
    assert_eq!(roundtrip(&copy), copy);
}
#[test]
fn stale_metadata_cannot_resurrect_removed_groups_or_overwrite_array_types() {
    let mut map = quantified_map();
    let run = &mut map.protein_identifications[0];
    run.metadata
        .insert("indistinguishable_proteins_9".into(), "0.1,PH_0".into());
    run.metadata.insert(
        "indistinguishable_proteins_9_quantified_proteins".into(),
        vec!["stale".to_owned()].into(),
    );
    let copy = roundtrip(&map);
    assert_eq!(
        copy.protein_identifications[0]
            .indistinguishable_groups
            .len(),
        1
    );
    assert!(
        !copy.protein_identifications[0]
            .metadata
            .contains_key("indistinguishable_proteins_9")
    );
    map.protein_identifications[0].indistinguishable_groups[0]
        .integer_data_arrays
        .push(DataArray::new("custom_float", vec![2]));
    let mut bytes = Vec::new();
    assert!(consensusxml::write(&mut bytes, &map).is_err());
    assert!(bytes.is_empty());
}
#[test]
fn centroid_commit_matches_source_empty_handle_convention_and_ids() {
    let xml=b"<consensusXML><mapList count=\"0\"/><consensusElementList><consensusElement id=\"e_nested_47\"><centroid rt=\"1\" mz=\"2\" it=\"3\"/><groupedElementList/></consensusElement></consensusElementList></consensusXML>";
    let map = consensusxml::read(Cursor::new(xml)).unwrap();
    assert_eq!(map.features[0].unique_id, 47);
    assert_eq!(map.features[0].rt, 0.0);
    assert_eq!(map.features[0].intensity, 0.0);
    let mut bad = map.clone();
    bad.features[0].rt = 1.0;
    assert!(consensusxml::write(Vec::new(), &bad).is_err());
    let text = String::from_utf8(xml.to_vec())
        .unwrap()
        .replace("e_nested_47", "e_bad");
    assert_eq!(
        consensusxml::read(Cursor::new(text)).unwrap().features[0].unique_id,
        0
    );
    let text = String::from_utf8(xml.to_vec())
        .unwrap()
        .replace("e_nested_47", "e_18446744073709551616");
    assert!(consensusxml::read(Cursor::new(text)).is_err());
}
#[test]
fn malformed_references_and_unrepresentable_state_fail_before_output() {
    let text = String::from_utf8(SOURCE.to_vec())
        .unwrap()
        .replace("protein_refs=\"PH_1\"", "protein_refs=\"PH_missing\"");
    let mut map = source();
    let before = map.clone();
    assert!(consensusxml::read_into(Cursor::new(text), &mut map, &ReadOptions::default()).is_err());
    assert_eq!(map, before);
    map.unassigned_peptide_identifications[0].identifier = "orphan".into();
    let mut bytes = Vec::new();
    assert!(consensusxml::write(&mut bytes, &map).is_err());
    assert!(bytes.is_empty());
    map = before;
    map.metadata.insert(
        "unit".into(),
        MetaValue::from(1i64)
            .with_unit(Unit::new("UO:0000010", "second", "UO").unwrap())
            .unwrap(),
    );
    assert!(consensusxml::write(&mut bytes, &map).is_err());
    assert!(bytes.is_empty());
}
#[test]
fn payload_work_and_xml_limits_are_independent_and_atomic() {
    let map = source();
    for options in [
        WriteOptions {
            max_xml_bytes: 32,
            ..Default::default()
        },
        WriteOptions {
            max_records: 2,
            ..Default::default()
        },
        WriteOptions {
            max_payload_bytes: 32,
            ..Default::default()
        },
        WriteOptions {
            max_work: 2,
            ..Default::default()
        },
    ] {
        let mut output = Vec::new();
        assert!(consensusxml::write_with_options(&mut output, &map, &options).is_err());
        assert!(output.is_empty());
    }
    for options in [
        ReadOptions {
            max_xml_bytes: 32,
            ..Default::default()
        },
        ReadOptions {
            max_records: 2,
            ..Default::default()
        },
        ReadOptions {
            max_payload_bytes: 32,
            ..Default::default()
        },
        ReadOptions {
            max_work: 2,
            ..Default::default()
        },
    ] {
        assert!(consensusxml::read_with_options(Cursor::new(SOURCE), &options).is_err());
    }
}
#[test]
fn writer_preflight_rejects_unbudgeted_validation_before_output() {
    let mut map = source();
    map.features[0].mz = f64::NAN;
    let mut output = Vec::new();
    let error = consensusxml::write_with_options(
        &mut output,
        &map,
        &WriteOptions {
            max_work: 1,
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(error.to_string().contains("limit exceeded"), "{error}");
    assert!(output.is_empty());

    let mut map = source();
    // Unrepresentable software payload must be rejected before traversal and
    // validation of the remaining map, including its invalid feature.
    map.features[0].mz = f64::NAN;
    map.data_processing[0].software.cv_terms.metadata.insert(
        "unrepresentable".into(),
        MetaValue::from(vec!["value".to_string(); 4096]),
    );
    let error = consensusxml::write(&mut output, &map).unwrap_err();
    assert!(
        error.to_string().contains("software CV terms or metadata"),
        "{error}"
    );
    assert!(output.is_empty());
}
#[test]
fn inconsistent_map_references_are_retained_and_reported() {
    let mut map = ConsensusMap::default();
    let mut f = ConsensusFeature::default();
    f.insert(FeatureHandle {
        map_index: 4,
        rt: 1.0,
        mz: 2.0,
        intensity: 3.0,
        ..Default::default()
    })
    .unwrap();
    map.features.push(f);
    let mut bytes = Vec::new();
    consensusxml::write(&mut bytes, &map).unwrap();
    let report = consensusxml::read_report(
        Cursor::new(bytes),
        &ReadOptions::default(),
        ModificationsDB::global(),
    )
    .unwrap();
    assert_eq!(report.map.features[0].handles()[0].map_index, 4);
    assert!(report.warnings.iter().any(|s| s.contains("inconsistent")));
}
#[test]
fn serialized_source_fixture_passes_original_schema() {
    let dir = Directory::new();
    let path = dir.0.join("map.consensusXML");
    consensusxml::store(&path, &source()).unwrap();
    let output = std::process::Command::new("xmllint")
        .args(["--nonet", "--noout", "--schema"])
        .arg(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/data/consensusxml/ConsensusXML_1_7.xsd"
        ))
        .arg(&path)
        .output();
    match output {
        Ok(result) => assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        ),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            eprintln!("xmllint unavailable: original schema test skipped")
        }
        Err(e) => panic!("schema validator failed: {e}"),
    }
}

#[test]
fn source_defined_modifications_travel_with_assigned_and_unassigned_sequences() {
    use openms::chemistry::{
        AASequence, EmpiricalFormula, ModificationRecord, ResidueModification,
    };
    use openms::identification::{PeptideHit, PeptideIdentification, ProteinIdentification};
    use openms::kernel::ColumnHeader;
    let mut registry = ModificationsDB::global().clone();
    let mut definitions = Vec::new();
    for (name, origin, text) in [
        ("TestCXML:Assigned", 'K', "C2H2O"),
        ("TestCXML:Unassigned", 'R', "CH2"),
    ] {
        let formula: EmpiricalFormula = text.parse().unwrap();
        definitions.push(
            ResidueModification::from_record(ModificationRecord {
                name: name.into(),
                origin: Some(origin),
                diff_mono_mass: formula.mono_mass(),
                diff_formula: formula,
                ..Default::default()
            })
            .unwrap(),
        );
    }
    registry.extend_records(definitions).unwrap();
    let mut map = ConsensusMap::default();
    map.column_headers.insert(
        0,
        ColumnHeader {
            filename: "file0.mzML".into(),
            size: 1,
            ..Default::default()
        },
    );
    map.protein_identifications.push(ProteinIdentification {
        identifier: "run4b".into(),
        date_time: Some("2026-09-10T12:00:00".into()),
        ..Default::default()
    });
    let peptide = |sequence: &str| PeptideIdentification {
        identifier: "run4b".into(),
        hits: vec![PeptideHit {
            sequence: AASequence::parse_with_registry(sequence, &registry).unwrap(),
            ..Default::default()
        }],
        ..Default::default()
    };
    let mut feature = ConsensusFeature::default();
    feature.rt = 100.0;
    feature.mz = 500.0;
    feature.intensity = 1000.0;
    feature
        .insert(FeatureHandle::new(0, &feature.base))
        .unwrap();
    feature
        .peptide_identifications
        .push(peptide("PEPK(TestCXML:Assigned)IDE"));
    map.features.push(feature);
    map.unassigned_peptide_identifications
        .push(peptide("PEPR(TestCXML:Unassigned)IDE"));
    let mut bytes = Vec::new();
    consensusxml::write_with_registry(&mut bytes, &map, &WriteOptions::default(), &registry)
        .unwrap();
    let text = std::str::from_utf8(&bytes).unwrap();
    assert!(text.contains("1|TestCXML:Assigned|TestCXML:Assigned (K)|"));
    assert!(text.contains("1|TestCXML:Unassigned|TestCXML:Unassigned (R)|"));
    let fresh = ModificationsDB::default();
    let copy =
        consensusxml::read_with_registry(Cursor::new(bytes), &ReadOptions::default(), &fresh)
            .unwrap();
    assert!(fresh.is_empty());
    let original = &map.features[0].peptide_identifications[0].hits[0].sequence;
    let retained = &copy.features[0].peptide_identifications[0].hits[0].sequence;
    assert_eq!(retained, original);
    assert_eq!(retained.formula().unwrap(), original.formula().unwrap());
    assert_eq!(
        copy.unassigned_peptide_identifications[0].hits[0].sequence,
        map.unassigned_peptide_identifications[0].hits[0].sequence
    );
}

#[test]
fn many_numbered_groups_preserve_numeric_order_across_lexical_boundaries() {
    let mut map = source();
    let accession = map.protein_identifications[0].hits[0].accession.clone();
    map.protein_identifications[0].protein_groups = (0..300)
        .map(|i| ProteinGroup {
            probability: 0.9,
            accessions: vec![accession.clone()],
            float_data_arrays: vec![DataArray::new("measured", vec![i as f32])],
            ..Default::default()
        })
        .collect();
    let copy = roundtrip(&map);
    let groups = &copy.protein_identifications[0].protein_groups;
    assert_eq!(groups.len(), 300);
    for (i, group) in groups.iter().enumerate() {
        assert_eq!(
            group
                .float_data_arrays
                .iter()
                .find(|a| a.name == "measured")
                .unwrap()
                .data,
            [i as f32]
        );
    }
}
