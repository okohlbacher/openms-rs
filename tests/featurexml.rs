#![cfg(feature = "featurexml")]
use openms::format::FileType;
use openms::format::featurexml::{self, FeatureFileOptions, Limits, ReadOptions, WriteOptions};
use openms::kernel::{ConvexHull2D, Feature, FeatureMap, Point2D};
use openms::metadata::{MetaValue, ProcessingAction};

const SOURCE: &[u8] = include_bytes!("data/featurexml_source_1.featureXML");
const OPTIONS: &[u8] = include_bytes!("data/featurexml_source_options.featureXML");
fn opts(options: FeatureFileOptions) -> ReadOptions {
    ReadOptions {
        feature_options: options,
        ..Default::default()
    }
}
fn simple(body: &str) -> Vec<u8> {
    format!(
        "<featureMap version=\"1.9\"><featureList count=\"1\">{body}</featureList></featureMap>"
    )
    .into_bytes()
}
fn feature(id: u64, rt: f64, intensity: f32) -> Feature {
    let mut f = Feature::new(rt, 300.0, intensity);
    f.unique_id = id;
    f
}

#[test]
fn source_fixture_geometry_metadata_processing_and_identifications() {
    let map = featurexml::read(SOURCE).unwrap();
    assert_eq!(map.identifier, "lsid");
    assert_eq!(map.len(), 2);
    assert_eq!(map.features[0].unique_id, 1000);
    assert_eq!(
        (
            map.features[0].rt,
            map.features[0].mz,
            map.features[0].intensity
        ),
        (25.0, 0.0, 300.0)
    );
    assert_eq!(
        (
            map.features[1].rt,
            map.features[1].mz,
            map.features[1].intensity
        ),
        (0.0, 35.0, 500.0)
    );
    assert_eq!(map.features[0].subordinates.len(), 2);
    assert_eq!(map.features[0].subordinates[1].unique_id, 2001);
    assert_eq!(map.features[0].subordinates[1].mz, 11.0);
    assert_eq!(
        map.features[1].convex_hulls[0].hull_points(),
        vec![Point2D::new(1.5, 1.8), Point2D::new(2.5, 2.8)]
    );
    let meta = &map.features[0].metadata;
    assert_eq!(meta["myIntList"].as_integer_list().unwrap(), &[1, 10, 12]);
    assert_eq!(
        meta["myDoubleList"].as_float_list().unwrap(),
        &[1.111, 10.999, 12.45]
    );
    assert_eq!(
        meta["myStringList"].as_string_list().unwrap(),
        &["myABC1", "Stuff", "12"]
    );
    assert_eq!(map.data_processing.len(), 2);
    assert_eq!(map.data_processing[0].software.name, "Software1");
    assert!(
        map.data_processing[0]
            .actions
            .contains(&ProcessingAction::Deisotoping)
    );
    assert!(
        map.data_processing[1]
            .actions
            .contains(&ProcessingAction::BaselineReduction)
    );
    assert_eq!(
        map.data_processing[0].completion_time.unwrap().to_string(),
        "2001-02-03 04:05:07"
    );
    assert_eq!(map.protein_identifications.len(), 2);
    assert_eq!(map.protein_identifications[0].hits.len(), 2);
    assert_eq!(map.unassigned_peptide_identifications.len(), 2);
    assert_eq!(
        map.features[0].peptide_identifications[1].hits[1]
            .sequence
            .to_string(),
        "D"
    );
    assert_eq!(
        map.unassigned_peptide_identifications[1].hits[1]
            .sequence
            .to_string(),
        "H"
    );
    assert_eq!(
        map.features[0].peptide_identifications[0].hits[0].evidences[0].protein_accession,
        "urn:lsid:rumpelstielzchen"
    );
    let mut encoded = Vec::new();
    featurexml::write(&mut encoded, &map).unwrap();
    assert_eq!(featurexml::read(encoded.as_slice()).unwrap(), map);
}

#[test]
fn source_options_and_passive_size_only() {
    let mut o = ReadOptions::default();
    o.feature_options.rt_range = Some(1.5..4.5);
    assert_eq!(featurexml::read_with_options(OPTIONS, &o).unwrap().len(), 5);
    o.feature_options.mz_range = Some(1025.0..2000.0);
    assert_eq!(featurexml::read_with_options(OPTIONS, &o).unwrap().len(), 3);
    o.feature_options.intensity_range = Some(290.0..310.0);
    assert_eq!(featurexml::read_with_options(OPTIONS, &o).unwrap().len(), 1);
    let o = opts(FeatureFileOptions {
        size_only: true,
        ..Default::default()
    });
    assert_eq!(featurexml::read_with_options(SOURCE, &o).unwrap().len(), 2);
    assert_eq!(featurexml::read_size(SOURCE, &o).unwrap(), 2);
    let o = opts(FeatureFileOptions {
        metadata_only: true,
        ..Default::default()
    });
    let map = featurexml::read_with_options(OPTIONS, &o).unwrap();
    assert!(map.is_empty());
    assert_eq!(map.identifier, "lsid2");
    assert_eq!(featurexml::read_size(SOURCE, &o).unwrap(), 0);
}

#[test]
fn soft_stop_does_not_parse_feature_payload() {
    let partial =
        b"<featureMap version=\"1.9\" document_id=\"header\"><featureList count=\"123\"><broken";
    assert_eq!(
        featurexml::read_size(partial.as_slice(), &ReadOptions::default()).unwrap(),
        123
    );
    let o = opts(FeatureFileOptions {
        metadata_only: true,
        ..Default::default()
    });
    assert_eq!(
        featurexml::read_with_options(partial.as_slice(), &o)
            .unwrap()
            .identifier,
        "header"
    );
    let missing = b"<featureMap><featureList count=\"nonsense\">";
    assert_eq!(featurexml::read_size(missing.as_slice(), &o).unwrap(), 0);
    assert!(featurexml::read(partial.as_slice()).is_err());
}

#[test]
fn half_open_filters_apply_independently_at_every_feature_level() {
    let mut parent = feature(1, 1.0, 10.0);
    parent.subordinates = vec![
        feature(2, 0.0, 20.0),
        feature(3, 1.0, 30.0),
        feature(4, 2.0, 40.0),
    ];
    parent.subordinates[1]
        .subordinates
        .push(feature(5, 1.5, 50.0));
    let map = FeatureMap::from_features(vec![parent, feature(6, 2.0, 60.0)]);
    let mut bytes = Vec::new();
    featurexml::write(&mut bytes, &map).unwrap();
    let o = opts(FeatureFileOptions {
        rt_range: Some(1.0..2.0),
        ..Default::default()
    });
    let filtered = featurexml::read_with_options(bytes.as_slice(), &o).unwrap();
    assert_eq!(filtered.len(), 1);
    assert_eq!(filtered.features[0].subordinates.len(), 1);
    assert_eq!(filtered.features[0].subordinates[0].unique_id, 3);
    assert_eq!(
        filtered.features[0].subordinates[0].subordinates[0].unique_id,
        5
    );
    let o = opts(FeatureFileOptions {
        rt_range: Some(2.0..2.0),
        ..Default::default()
    });
    assert!(
        featurexml::read_with_options(bytes.as_slice(), &o)
            .unwrap()
            .is_empty()
    );
    let o = opts(FeatureFileOptions {
        rt_range: Some(2.0..1.0),
        ..Default::default()
    });
    assert!(
        featurexml::read_with_options(bytes.as_slice(), &o)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn omitted_subtrees_ignore_their_scientific_values_and_legacy_description() {
    let bytes = simple(
        "<feature id=\"f_7\"><subordinate><feature id=\"bad\"><position dim=\"9\">NaN</position></feature></subordinate><convexhull><pt x=\"NaN\"/></convexhull><description><legacy arbitrary=\"yes\"/></description></feature>",
    );
    let o = opts(FeatureFileOptions {
        load_convex_hulls: false,
        load_subordinates: false,
        ..Default::default()
    });
    assert_eq!(
        featurexml::read_with_options(bytes.as_slice(), &o)
            .unwrap()
            .features[0]
            .unique_id,
        7
    );
    assert!(featurexml::read(bytes.as_slice()).is_err());
    assert_eq!(
        featurexml::read(include_bytes!("data/featurexml_source_old.featureXML").as_slice())
            .unwrap()
            .len(),
        1
    );
    let original = featurexml::read(OPTIONS).unwrap();
    let loaded = featurexml::read_with_options(OPTIONS, &o).unwrap();
    assert_eq!(original.len(), loaded.len());
    assert!(
        loaded
            .features
            .iter()
            .all(|f| f.convex_hulls.is_empty() && f.subordinates.is_empty())
    );
}

#[test]
fn legacy_hull_and_id_forms_and_top_level_width_bridge() {
    let bytes=b"<featureMap id=\"fm_4\" unique_id=\"5\"><featureList count=\"1\"><feature id=\"f_4_8\"><convexhull nr=\"0\"><hullpoint><hposition dim=\"0\">3.25</hposition><hposition dim=\"1\">10</hposition></hullpoint><hullpoint><hposition dim=\"1\">20</hposition><hposition dim=\"0\">4.25</hposition></hullpoint></convexhull><subordinate><feature id=\"f_8_9\"><UserParam name=\"FWHM\" type=\"float\" value=\"2\"/></feature></subordinate><userParam name=\"FWHM\" type=\"float\" value=\"1.5\"/></feature></featureList></featureMap>";
    let map = featurexml::read(bytes.as_slice()).unwrap();
    assert_eq!(map.unique_id, 5);
    assert_eq!(map.features[0].unique_id, 8);
    assert_eq!(map.features[0].width, 1.5);
    assert_eq!(map.features[0].subordinates[0].width, 0.0);
    assert_eq!(
        map.features[0].subordinates[0].metadata["FWHM"]
            .as_f64()
            .unwrap(),
        2.0
    );
    assert_eq!(
        map.features[0].convex_hulls[0].hull_points(),
        vec![Point2D::new(3.25, 10.0), Point2D::new(4.25, 20.0)]
    );
}

#[test]
fn native_typed_metadata_and_width_setter_roundtrip() {
    let mut f = feature(1, 3.0, 5.0);
    f.set_width(2.25).unwrap();
    f.metadata.insert("count".into(), 7_i64.into());
    f.metadata
        .insert("labels".into(), vec!["a".to_owned(), "b".to_owned()].into());
    f.convex_hulls.push(
        ConvexHull2D::from_points(&[
            Point2D::new(1.0, 2.0),
            Point2D::new(2.0, 2.0),
            Point2D::new(3.0, 2.0),
        ])
        .unwrap(),
    );
    let mut map = FeatureMap::from_features(vec![f]);
    map.unique_id = 123;
    map.metadata
        .insert("samples".into(), MetaValue::from(vec![1_i64, 2]));
    let mut bytes = Vec::new();
    featurexml::write(&mut bytes, &map).unwrap();
    let rt = featurexml::read(bytes.as_slice()).unwrap();
    assert_eq!(rt.unique_id, 123);
    assert_eq!(rt.features[0].width, 2.25);
    assert_eq!(rt.metadata, map.metadata);
    assert_eq!(rt.features[0].metadata, map.features[0].metadata);
    assert_eq!(
        rt.features[0].convex_hulls[0].hull_points(),
        vec![Point2D::new(1.0, 2.0), Point2D::new(3.0, 2.0)]
    );
}

#[test]
fn parse_and_write_failures_leave_destination_unchanged() {
    let mut target = FeatureMap::from_features(vec![feature(19, 1.0, 1.0)]);
    let original = target.clone();
    for input in [
        simple("<feature id=\"1\"><position dim=\"8\">2</position></feature>"),
        simple("<feature id=\"1\"><intensity>1e50</intensity></feature>"),
        simple("<feature id=\"1\"><UserParam name=\"x\" type=\"unknown\" value=\"1\"/></feature>"),
    ] {
        assert!(
            featurexml::read_into(input.as_slice(), &mut target, &ReadOptions::default()).is_err()
        );
        assert_eq!(target, original);
    }
    let mut bad = original.clone();
    bad.features.push(feature(19, 2.0, 2.0));
    let mut output = b"existing".to_vec();
    assert!(featurexml::write(&mut output, &bad).is_err());
    assert_eq!(output, b"existing");
    bad.features[1].unique_id = 20;
    bad.features[1].intensity = f32::NAN;
    assert!(featurexml::write(&mut output, &bad).is_err());
    assert_eq!(output, b"existing");
}

#[test]
fn bounded_work_payload_depth_and_records_are_atomic() {
    let mut target = FeatureMap::from_features(vec![feature(1, 0.0, 0.0)]);
    let before = target.clone();
    for limits in [
        Limits {
            max_work: 1,
            ..Default::default()
        },
        Limits {
            max_payload_bytes: 0,
            ..Default::default()
        },
        Limits {
            max_records: 2,
            ..Default::default()
        },
        Limits {
            max_depth: 0,
            ..Default::default()
        },
    ] {
        assert!(
            featurexml::read_into(
                SOURCE,
                &mut target,
                &ReadOptions {
                    limits,
                    ..Default::default()
                }
            )
            .is_err()
        );
        assert_eq!(target, before);
        let mut output = b"untouched".to_vec();
        let source = featurexml::read(SOURCE).unwrap();
        assert!(
            featurexml::write_with_options(&mut output, &source, &WriteOptions { limits }).is_err()
        );
        assert_eq!(output, b"untouched");
    }
}

#[test]
fn unsupported_software_payload_is_rejected_before_copying_or_validation() {
    use openms::chemistry::ModificationsDB;
    use openms::metadata::{CVTerm, DataProcessing};

    for has_cv_term in [false, true] {
        let mut processing = DataProcessing::default();
        if has_cv_term {
            processing
                .software
                .cv_terms
                .add(CVTerm::new("MS:1000799", "software".repeat(4096), "MS"))
                .unwrap();
        } else {
            processing.software.cv_terms.metadata.insert(
                "unsupported".into(),
                vec!["payload".to_owned(); 4096].into(),
            );
        }
        let map = FeatureMap {
            data_processing: vec![processing],
            ..Default::default()
        };
        let mut output = b"untouched".to_vec();
        let error = featurexml::write_with_registry(
            &mut output,
            &map,
            &WriteOptions {
                limits: Limits {
                    max_work: 1,
                    max_payload_bytes: 0,
                    ..Default::default()
                },
            },
            &ModificationsDB::default(),
        )
        .unwrap_err();
        assert!(matches!(error, openms::Error::Unsupported(message)
            if message == "map XML software CV terms and metadata"));
        assert_eq!(output, b"untouched");
    }
}

#[test]
fn plain_gzip_bzip2_paths_and_loaded_file_identity() {
    let map = featurexml::read(SOURCE).unwrap();
    for suffix in ["featureXML", "featureXML.gz", "featureXML.bz2"] {
        let path =
            std::env::temp_dir().join(format!("openms-featurexml-{}.{suffix}", std::process::id()));
        featurexml::store(&path, &map).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        if suffix.ends_with(".gz") {
            assert!(bytes.starts_with(&[0x1f, 0x8b]));
        }
        if suffix.ends_with(".bz2") {
            assert!(bytes.starts_with(b"BZh"));
        }
        let mut loaded = featurexml::load(&path).unwrap();
        assert_eq!(loaded.loaded_file_type, FileType::FeatureXml);
        assert_eq!(loaded.loaded_file_path, path.to_string_lossy());
        assert_eq!(
            featurexml::load_size(&path, &ReadOptions::default()).unwrap(),
            2
        );
        loaded.loaded_file_path.clear();
        loaded.loaded_file_type = FileType::Unknown;
        assert_eq!(loaded, map);
        std::fs::remove_file(path).unwrap();
    }
}

#[test]
fn owning_subordinate_metadata_and_orphan_ids_have_checked_native_behavior() {
    let source = std::str::from_utf8(SOURCE).unwrap();
    let input = source.replacen("</feature>", "<PeptideIdentification identification_run_ref=\"PI_0\" score_type=\"test\" higher_score_better=\"true\"/><UserParam type=\"int\" name=\"owned\" value=\"7\"/></feature>", 1);
    let map = featurexml::read(input.as_bytes()).unwrap();
    assert_eq!(
        map.features[0].subordinates[0].metadata["owned"]
            .as_i64()
            .unwrap(),
        7
    );
    assert!(!map.features[0].metadata.contains_key("owned"));
    let mut encoded = Vec::new();
    featurexml::write(&mut encoded, &map).unwrap();
    assert_eq!(featurexml::read(encoded.as_slice()).unwrap(), map);
    let mut orphan = map;
    orphan.unassigned_peptide_identifications[0].identifier = "missing run".into();
    let mut output = b"keep".to_vec();
    assert!(featurexml::write(&mut output, &orphan).is_err());
    assert_eq!(output, b"keep");
}

#[test]
fn size_prefix_stops_before_invalid_encoding_large_tail_or_late_io() {
    let prefix = b"<featureMap version=\"1.9\"><!-- <featureList count='99'> --><featureList count=\"1000000000\">";
    let mut bytes = prefix.to_vec();
    bytes.extend_from_slice(&[0xff; 4096]);
    let o = ReadOptions {
        limits: Limits {
            max_xml_bytes: prefix.len() as u64,
            ..Default::default()
        },
        ..Default::default()
    };
    assert_eq!(
        featurexml::read_size(bytes.as_slice(), &o).unwrap(),
        1_000_000_000
    );
    struct FailAfter<'a>(&'a [u8]);
    impl std::io::Read for FailAfter<'_> {
        fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
            if self.0.is_empty() {
                return Err(std::io::Error::other("payload should not be read"));
            }
            self.0.read(out)
        }
    }
    impl std::io::BufRead for FailAfter<'_> {
        fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
            if self.0.is_empty() {
                Err(std::io::Error::other("payload should not be read"))
            } else {
                Ok(self.0)
            }
        }
        fn consume(&mut self, n: usize) {
            self.0 = &self.0[n..];
        }
    }
    assert_eq!(
        featurexml::read_size(FailAfter(prefix), &o).unwrap(),
        1_000_000_000
    );
    let text = "<?xml version=\"1.0\" encoding=\"UTF-16\"?><featureMap><featureList count=\"7\">";
    let mut utf16 = vec![0xff, 0xfe];
    for c in text.encode_utf16() {
        utf16.extend_from_slice(&c.to_le_bytes());
    }
    utf16.push(0xff); // invalid odd trailing byte must not be decoded
    assert_eq!(
        featurexml::read_size(utf16.as_slice(), &ReadOptions::default()).unwrap(),
        7
    );
}

#[test]
fn nontransportable_width_and_structured_groups_fail_before_output() {
    let mut map = FeatureMap::from_features(vec![feature(1, 1.0, 1.0)]);
    map.features[0].width = 2.0;
    let mut output = vec![42];
    assert!(featurexml::write(&mut output, &map).is_err());
    assert_eq!(output, [42]);
    map.features[0].set_width(2.0).unwrap();
    let mut sub = feature(2, 1.0, 1.0);
    sub.set_width(2.0).unwrap();
    map.features[0].subordinates.push(sub);
    assert!(featurexml::write(&mut output, &map).is_err());
    assert_eq!(output, [42]);
    let mut map = featurexml::read(SOURCE).unwrap();
    map.protein_identifications[0]
        .protein_groups
        .push(openms::identification::ProteinGroup::default());
    assert!(featurexml::write(&mut output, &map).is_err());
    assert_eq!(output, [42]);
}

#[test]
fn portable_definitions_cover_assigned_subordinate_and_unassigned_ids() {
    use openms::chemistry::{AASequence, ModificationRecord, ModificationsDB, ResidueModification};
    use openms::identification::{PeptideHit, PeptideIdentification, ProteinIdentification};
    let sequences = {
        let records = [
            ("FeatureLabA", "O"),
            ("FeatureLabB", "O2"),
            ("FeatureLabC", "H2"),
        ]
        .into_iter()
        .map(|(name, formula)| {
            let formula: openms::chemistry::EmpiricalFormula = formula.parse().unwrap();
            ResidueModification::from_record(ModificationRecord {
                name: name.into(),
                origin: Some('M'),
                diff_mono_mass: formula.mono_mass(),
                diff_formula: formula,
                ..Default::default()
            })
            .unwrap()
        })
        .collect();
        let registry = ModificationsDB::from_records(records).unwrap();
        ["AM(FeatureLabA)K", "AM(FeatureLabB)K", "AM(FeatureLabC)K"]
            .map(|text| AASequence::parse_with_registry(text, &registry).unwrap())
    }; // chemical handles must remain valid after the caller registry is gone
    let make_id = |sequence: AASequence| PeptideIdentification {
        identifier: "run".into(),
        hits: vec![PeptideHit {
            sequence,
            ..Default::default()
        }],
        ..Default::default()
    };
    let mut parent = feature(1, 1.0, 1.0);
    parent
        .peptide_identifications
        .push(make_id(sequences[0].clone()));
    let mut sub = feature(2, 1.0, 1.0);
    sub.peptide_identifications
        .push(make_id(sequences[1].clone()));
    parent.subordinates.push(sub);
    let mut map = FeatureMap::from_features(vec![parent]);
    map.protein_identifications.push(ProteinIdentification {
        identifier: "run".into(),
        date_time: Some("2026-09-10T12:00:00".into()),
        ..Default::default()
    });
    map.unassigned_peptide_identifications
        .push(make_id(sequences[2].clone()));
    let before = map.clone();
    let mut bytes = Vec::new();
    featurexml::write(&mut bytes, &map).unwrap();
    assert_eq!(map, before);
    let restored = featurexml::read_with_registry(
        bytes.as_slice(),
        &ReadOptions::default(),
        &ModificationsDB::default(),
    )
    .unwrap();
    assert_eq!(restored.features, map.features);
    assert_eq!(
        restored.unassigned_peptide_identifications,
        map.unassigned_peptide_identifications
    );
    for (id, expected) in [
        &restored.features[0].peptide_identifications[0],
        &restored.features[0].subordinates[0].peptide_identifications[0],
        &restored.unassigned_peptide_identifications[0],
    ]
    .into_iter()
    .zip(sequences)
    {
        assert_eq!(
            id.hits[0].sequence.formula().unwrap(),
            expected.formula().unwrap()
        );
    }
    let definitions = restored.protein_identifications[0]
        .search_parameters
        .metadata["modification_definitions"]
        .as_str()
        .unwrap();
    for name in ["FeatureLabA", "FeatureLabB", "FeatureLabC"] {
        assert!(definitions.contains(name));
    }
}
