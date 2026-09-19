#![cfg(feature = "featurexml")]
use openms::format::FileType;
use openms::format::featurexml::{
    self, Allowance, FeatureFileOptions, InputScaling, Limits, OutputScaling, ReadOptions,
    WriteOptions,
};
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
    // A width that no longer matches its `FWHM` mirror. A non-finite value is
    // no longer one of these: the Release build writes those, and so does this
    // writer now (`the_release_spelling_of_a_nonfinite_value_is_written_for_every_field`).
    bad.features[1].base.width = 3.0;
    assert!(featurexml::write(&mut output, &bad).is_err());
    assert_eq!(output, b"existing");
    // The same for a finite negative width, which no source path produces.
    bad.features[1].base.width = -1.0;
    bad.features[1]
        .metadata
        .insert("FWHM".into(), MetaValue::try_from(-1.0).unwrap());
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
            featurexml::write_with_options(
                &mut output,
                &source,
                &WriteOptions {
                    limits,
                    ..Default::default()
                }
            )
            .is_err()
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
                ..Default::default()
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

// --- Size-derived ceilings -------------------------------------------------
//
// The reader's ceilings used to be fixed, which made the effective input
// ceiling 12,500,000 bytes: `max_xml_bytes` (64 MiB) reduced by
// `max_payload_bytes / 8` and by `max_work / 4`. These pin the former
// behaviour, the growth that replaced it, and that a document still cannot buy
// more than its own size earns. See `docs/FEATUREXML_SCALE_SUPPORT.md`.

/// The former fixed ceilings, exactly as this adapter applied them.
fn former() -> ReadOptions {
    ReadOptions {
        limits: Limits::former(),
        scaling: InputScaling::default().fixed(),
        ..Default::default()
    }
}

/// A valid one-feature document padded to `bytes` with an XML comment, which
/// costs the reader one event however long it is.
fn padded(bytes: usize) -> Vec<u8> {
    let head = b"<featureMap version=\"1.9\"><!--".to_vec();
    let tail = b"--><featureList count=\"1\"><feature id=\"f_7\"><position dim=\"0\">1</position><position dim=\"1\">2</position><intensity>3</intensity></feature></featureList></featureMap>".to_vec();
    let mut document = head;
    document.resize(bytes.saturating_sub(tail.len()), b'.');
    document.extend_from_slice(&tail);
    document
}

#[test]
fn the_former_fixed_ceilings_still_refuse_exactly_what_they_refused() {
    // The benchmark's FileInfo failed on the 59.6 MiB featureXML with
    // "identification XML byte limit exceeded" because the former fixed
    // ceilings made the decode limit min(64 MiB, 256 MiB / 8, 50,000,000 / 4)
    // = 12,500,000 bytes, four work units per decoded byte out of 50,000,000
    // being the binding term. Under `Limits::former()` with fixed scaling a
    // document of that size is still refused; the decode limit is now
    // `max_xml_bytes` alone, so the refusal now names the ceiling the coupling
    // stood in for.
    let over = padded(12_500_001);
    assert!(
        matches!(
            featurexml::read_with_options(over.as_slice(), &former()),
            Err(openms::Error::Parse { line: 0, ref message }) if message == "XML work limit exceeded"
        ),
        "{:?}",
        featurexml::read_with_options(over.as_slice(), &former())
    );
    // A document the former ceilings did admit is admitted unchanged.
    let under = padded(12_000_000);
    assert_eq!(under.len(), 12_000_000);
    assert_eq!(
        featurexml::read_with_options(under.as_slice(), &former())
            .unwrap()
            .len(),
        1
    );
    // `max_xml_bytes`, the one ceiling nothing can be derived from, still
    // refuses a document larger than itself, with the benchmark's message.
    let capped = ReadOptions {
        limits: Limits {
            max_xml_bytes: 1_000,
            ..Limits::former()
        },
        scaling: InputScaling::default().fixed(),
        ..Default::default()
    };
    assert!(
        matches!(
            featurexml::read_with_options(padded(2_000).as_slice(), &capped),
            Err(openms::Error::Parse { line: 0, ref message })
                if message == "identification XML byte limit exceeded"
        ),
        "{:?}",
        featurexml::read_with_options(padded(2_000).as_slice(), &capped)
    );
    // The size-derived defaults read the document the former ceilings refused.
    assert_eq!(featurexml::read(over.as_slice()).unwrap().len(), 1);
}

#[test]
fn every_cumulative_ceiling_is_earned_by_the_documents_own_bytes() {
    let options = |scaling| ReadOptions {
        scaling,
        ..Default::default()
    };
    // Elements: both documents hold the same six, only their size differs. One
    // element per 8,192 bytes plus one for free: the short document cannot pay
    // for its own elements, the padded one can.
    let records = InputScaling {
        records: Allowance::every(1, 8_192),
        ..Default::default()
    };
    assert!(featurexml::read_with_options(padded(200).as_slice(), &options(records)).is_err());
    assert_eq!(
        featurexml::read_with_options(padded(65_536).as_slice(), &options(records))
            .unwrap()
            .len(),
        1
    );
    // Work and payload: 40 MB of document charges more than either former
    // fixed floor covers, so only the rate can pay for it.
    let big = padded(40_000_000);
    for scaling in [
        InputScaling::default().fixed(),
        InputScaling {
            work: Allowance::new(50_000_000, 1),
            ..Default::default()
        },
        InputScaling {
            payload_bytes: Allowance::new(256 * 1024 * 1024, 0),
            ..Default::default()
        },
    ] {
        assert!(featurexml::read_with_options(big.as_slice(), &options(scaling)).is_err());
    }
    assert_eq!(featurexml::read(big.as_slice()).unwrap().len(), 1);
}

#[test]
fn an_absolute_ceiling_still_wins_over_the_size_derived_one() {
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
            max_xml_bytes: 8,
            ..Default::default()
        },
    ] {
        assert!(
            featurexml::read_with_options(
                SOURCE,
                &ReadOptions {
                    limits,
                    ..Default::default()
                }
            )
            .is_err()
        );
    }
}

#[test]
fn streamed_features_are_not_retained_and_the_feature_list_must_come_last() {
    // FeatureXML_1_9.xsd puts featureList last in the featureMap sequence, and
    // the streaming reader has converted every feature by the time the list
    // closes, so anything after it could no longer inform them.
    let trailing = b"<featureMap version=\"1.9\"><featureList count=\"0\"></featureList>\
                     <UserParam type=\"string\" name=\"late\" value=\"x\"/></featureMap>";
    let error = featurexml::read(trailing.as_slice()).unwrap_err();
    assert!(
        matches!(&error, openms::Error::Unsupported(message)
            if message == "featureMap content after featureList"),
        "{error:?}"
    );
    // A non-feature child of featureList is still refused.
    let intruder = b"<featureMap version=\"1.9\"><featureList count=\"0\">\
                     <UserParam type=\"string\" name=\"x\" value=\"y\"/></featureList></featureMap>";
    assert!(featurexml::read(intruder.as_slice()).is_err());
    // An empty list never calls the streaming path, and still reads its header.
    let empty = b"<featureMap version=\"1.9\" document_id=\"lsid\">\
                  <featureList count=\"0\"></featureList></featureMap>";
    let map = featurexml::read(empty.as_slice()).unwrap();
    assert!(map.is_empty());
    assert_eq!(map.identifier, "lsid");
    // Nothing in a document without a featureList can be streamed, so it is
    // refused before a tree is built rather than held whole.
    let listless = b"<featureMap version=\"1.9\" document_id=\"lsid\"/>";
    assert!(
        matches!(
            featurexml::read(listless.as_slice()),
            Err(openms::Error::Parse { line: 0, ref message })
                if message == "featureMap requires featureList"
        ),
        "{:?}",
        featurexml::read(listless.as_slice())
    );
}

#[test]
fn the_writer_ceiling_is_earned_by_the_map_it_is_given() {
    let map = featurexml::read(SOURCE).unwrap();
    let options = |scaling| WriteOptions {
        scaling,
        ..Default::default()
    };
    let starved = OutputScaling {
        payload_bytes: Allowance::new(0, 1),
        ..Default::default()
    };
    let fed = OutputScaling {
        payload_bytes: Allowance::new(0, 1 << 20),
        ..Default::default()
    };
    let mut output = b"untouched".to_vec();
    assert!(featurexml::write_with_options(&mut output, &map, &options(starved)).is_err());
    assert_eq!(output, b"untouched");
    let mut output = Vec::new();
    featurexml::write_with_options(&mut output, &map, &options(fed)).unwrap();
    assert_eq!(featurexml::read(output.as_slice()).unwrap(), map);
}

// --- HPC scale -------------------------------------------------------------

/// The 59.6 MiB `FeatureFinderCentroided` map of the TOPP benchmark inputs.
const BENCH_SMALL: &str = "/ceph/ibmi/abi/oliver/bench/openms4/inputs/\
                           featurexml_small_pxd001819_ffc_50amol_r1/UPS1_50amol_R1.featureXML";
/// The 2.06 GiB `MassTraceExtractor` map of the TOPP benchmark inputs.
const BENCH_LARGE: &str = "/ceph/ibmi/abi/oliver/bench/openms4/inputs/\
                           featurexml_large_pxd001819_mte_500amol_r3/UPS1_500amol_R3.featureXML";

/// Closed retention-time, m/z and intensity extremes of `map`, as `FileInfo`
/// reports them.
fn extremes(map: &FeatureMap) -> [(f64, f64); 3] {
    let mut bounds = [(f64::INFINITY, f64::NEG_INFINITY); 3];
    for feature in &map.features {
        for (bound, value) in
            bounds
                .iter_mut()
                .zip([feature.rt, feature.mz, f64::from(feature.intensity)])
        {
            bound.0 = bound.0.min(value);
            bound.1 = bound.1.max(value);
        }
    }
    bounds
}

fn close(value: f64, expected: f64) -> bool {
    (value - expected).abs() <= 0.005 * expected.abs().max(1.0)
}

#[test]
#[ignore = "HPC scale: reads the 59.6 MiB and 2.06 GiB benchmark featureXML files by path"]
fn hpc_scale_benchmark_featurexml_files_load_and_round_trip() {
    for path in [BENCH_SMALL, BENCH_LARGE] {
        assert!(
            std::path::Path::new(path).exists(),
            "benchmark input missing: {path}"
        );
    }
    // The 2.06 GiB map, read the way FileInfo reads it: geometry and
    // subordinate payload skipped. Every number below is the C++ FileInfo
    // summary of the same file at core bc9cc12.
    let summary = FeatureFileOptions {
        load_convex_hulls: false,
        load_subordinates: false,
        ..Default::default()
    };
    let large = featurexml::load_with_options(BENCH_LARGE, &opts(summary)).unwrap();
    assert_eq!(large.len(), 826_019);
    assert_eq!(
        featurexml::load_size(BENCH_LARGE, &ReadOptions::default()).unwrap(),
        826_019
    );
    let bounds = extremes(&large);
    for (measured, expected) in bounds.iter().zip([
        (0.31, 9299.33),
        (350.08, 1799.98),
        (827.55, 6_252_709_888.0),
    ]) {
        assert!(close(measured.0, expected.0), "{measured:?} {expected:?}");
        assert!(close(measured.1, expected.1), "{measured:?} {expected:?}");
    }
    assert!(large.features.iter().all(|f| f.charge == 0));
    drop(large);

    // The 59.6 MiB map at full fidelity, geometry included, then a round trip
    // through the writer, which the former fixed payload ceiling also refused.
    let small = featurexml::load(BENCH_SMALL).unwrap();
    assert_eq!(small.len(), 42_789);
    assert!(
        small
            .features
            .iter()
            .filter(|f| !f.convex_hulls.is_empty())
            .count()
            > 40_000
    );
    let bounds = extremes(&small);
    for (measured, expected) in bounds.iter().zip([
        (37.42, 9288.61),
        (350.18, 1792.82),
        (2758.98, 8_515_100_160.0),
    ]) {
        assert!(close(measured.0, expected.0), "{measured:?} {expected:?}");
        assert!(close(measured.1, expected.1), "{measured:?} {expected:?}");
    }
    for (charge, expected) in [(2, 22_869), (3, 16_288), (4, 3_632)] {
        assert_eq!(
            small.features.iter().filter(|f| f.charge == charge).count(),
            expected
        );
    }
    let dir = openms::system::file::TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("hpc-round-trip.featureXML");
    featurexml::store(&path, &small).unwrap();
    let restored = featurexml::load(&path).unwrap();
    assert_eq!(restored.features, small.features);
    assert_eq!(restored.metadata, small.metadata);
}

/// The document the Release build writes when every feature field, meta value,
/// float-list entry and `FWHM` holds a non-finite value and the hull points
/// stay finite (`../oracle/featurexml-inf`, `drivers/probe_nonfinite.cpp` on
/// `ibminode06` against `openms4-release-bc9cc12-c19e494-174b576`, two
/// identical repetitions; `results/probe_scalars.r1.featureXML`).
///
/// Its four features were stored from `+inf`, `-inf`, `NaN` and `-NaN`.
const NONFINITE_RELEASE: &[u8] = include_bytes!("data/featurexml_nonfinite_release.featureXML");
/// The same document with non-finite hull points as well
/// (`results/probe_hulls.r1.featureXML`); the two differ in nothing else.
const NONFINITE_HULL_RELEASE: &[u8] =
    include_bytes!("data/featurexml_nonfinite_hull_release.featureXML");

/// The four spellings that document carries, in feature order.
const NONFINITE_CASES: [(&str, f64); 4] = [
    ("p_inf", f64::INFINITY),
    ("n_inf", f64::NEG_INFINITY),
    ("p_nan", f64::NAN),
    // The Release build writes a NaN of either sign as `NaN`, so the fourth
    // feature comes back as a positive quiet NaN.
    ("n_nan", f64::NAN),
];

/// `a` and `b` are the same value, NaN included; NaN payloads are not compared,
/// as the source does not preserve them through text either.
fn alike(a: f64, b: f64) -> bool {
    a == b || (a.is_nan() && b.is_nan())
}

/// **The Release build's spelling of a non-finite value, written for every
/// field it has.**
///
/// `NumericFormatting::appendNumeric` answers before it formats anything:
/// `NaN` for a NaN of either sign, then `-inf` or `inf`
/// (`src/common/include/OpenMS/CONCEPT/Detail/NumericFormatting.h:29-35`).
/// `writeFeature_` sends every scalar through `precisionWrapper` and
/// `writeUserParam_` sends a `DataValue` through the same conversion, so one
/// spelling covers positions, intensity, qualities, the overall quality and
/// every `float`/`floatList` meta value.
///
/// The map is the executed document read back, so nothing here is derived from
/// a value this port invented.
#[test]
fn the_release_spelling_of_a_nonfinite_value_is_written_for_every_field() {
    let map = featurexml::read(NONFINITE_RELEASE).unwrap();
    let mut written = Vec::new();
    featurexml::write(&mut written, &map).unwrap();
    let text = String::from_utf8(written).unwrap();
    for spelling in [
        "<position dim=\"0\">inf</position>",
        "<position dim=\"1\">inf</position>",
        "<position dim=\"0\">-inf</position>",
        "<position dim=\"1\">-inf</position>",
        "<position dim=\"0\">NaN</position>",
        "<intensity>inf</intensity>",
        "<intensity>-inf</intensity>",
        "<intensity>NaN</intensity>",
        "<quality dim=\"0\">inf</quality>",
        "<quality dim=\"1\">-inf</quality>",
        "<quality dim=\"0\">NaN</quality>",
        "<overallquality>inf</overallquality>",
        "<overallquality>-inf</overallquality>",
        "<overallquality>NaN</overallquality>",
    ] {
        assert!(text.contains(spelling), "{spelling} missing");
    }
    for value in ["value=\"inf\"", "value=\"-inf\"", "value=\"NaN\""] {
        assert!(text.contains(value), "{value} missing");
    }
    for list in ["[inf,1.5]", "[-inf,1.5]", "[NaN,1.5]"] {
        assert!(text.contains(list), "{list} missing");
    }
    // A NaN never carries a sign, whichever sign the value has, and no
    // alternative spelling is used for a value: `inf.0` is the one the
    // source's header records as unreadable, and it once produced it.
    // (`nan` also occurs inside the `p_nan` labels, so only values are
    // inspected.)
    for spelling in ["-NaN", "nan", "NAN", "inf.0", "infinity", "Infinity"] {
        for absent in [
            format!(">{spelling}<"),
            format!("=\"{spelling}\""),
            format!("[{spelling},"),
        ] {
            assert!(!text.contains(&absent), "{absent} present");
        }
    }
}

/// **Every value of that document survives reading, writing and reading
/// again**, checked against the bit patterns the pinned `FeatureXMLFile::load`
/// produced for it (`../oracle/featurexml-inf/results/probe_scalars.r1.txt`).
///
/// A subordinate's width stays 0 on both sides: the load hack restores a width
/// from `FWHM` on top-level features only (`FeatureXMLFile.cpp:57-66`), while
/// the meta value itself survives at every level.
#[test]
fn the_release_nonfinite_document_reads_back_with_the_values_the_source_reads() {
    let map = featurexml::read(NONFINITE_RELEASE).unwrap();
    assert_eq!(map.len(), 4);
    assert_eq!(map.metadata["map_float"].as_f64().unwrap(), f64::INFINITY);
    let round = {
        let mut written = Vec::new();
        featurexml::write(&mut written, &map).unwrap();
        featurexml::read(written.as_slice()).unwrap()
    };
    assert_eq!(round.metadata["map_float"].as_f64().unwrap(), f64::INFINITY);
    for source in [&map, &round] {
        for (index, (label, wide)) in NONFINITE_CASES.into_iter().enumerate() {
            let narrow = wide as f32;
            let f = &source.features[index];
            assert_eq!(f.metadata["label"].as_str().unwrap(), label);
            assert_eq!(f.charge, 2, "{label}");
            assert!(alike(f.rt, wide), "{label} rt");
            assert!(alike(f.mz, wide), "{label} mz");
            for (name, value) in [
                ("intensity", f.intensity),
                ("quality", f.quality),
                ("quality_rt", f.quality_rt),
                ("quality_mz", f.quality_mz),
                ("width", f.width),
            ] {
                assert!(alike(f64::from(value), f64::from(narrow)), "{label} {name}");
            }
            for name in ["FWHM", "probe_float"] {
                assert!(
                    alike(f.metadata[name].as_f64().unwrap(), wide),
                    "{label} {name}"
                );
            }
            let list = f.metadata["probe_floatlist"].as_float_list().unwrap();
            assert_eq!(list.len(), 2, "{label}");
            assert!(alike(list[0], wide), "{label} list");
            assert_eq!(list[1], 1.5, "{label} list");
            assert_eq!(
                f.convex_hulls[0].hull_points(),
                vec![Point2D::new(3.0, 4.0), Point2D::new(1.0, 2.0)],
                "{label} hull"
            );
            let sub = &f.subordinates[0];
            assert_eq!(
                sub.metadata["label"].as_str().unwrap(),
                format!("{label}_sub")
            );
            assert!(alike(sub.rt, wide), "{label} subordinate rt");
            assert_eq!(sub.width, 0.0, "{label} subordinate width");
            assert!(
                alike(sub.metadata["FWHM"].as_f64().unwrap(), wide),
                "{label} subordinate FWHM"
            );
        }
    }
}

/// **The spellings the pinned reader takes, and the two it refuses.**
///
/// `probe_read` took twelve of them through `FeatureXMLFile::load` on
/// `ibminode06` (`../oracle/featurexml-inf/results/spellings.tsv`,
/// `extract/make_spellings.py`, each spelling substituted into every float
/// place at once). Ten load and give the value below in every field. `inf.0`
/// and `1e999` make the load throw `ConversionError` — the first because the
/// trailing `.0` is left over, the second because `std::from_chars` reports
/// the literal as out of range — and this port refuses both.
#[test]
fn every_nonfinite_spelling_the_pinned_reader_takes_is_accepted() {
    let document = |spelling: &str| {
        simple(&format!(
            "<feature id=\"f_100\"><position dim=\"0\">{spelling}</position>\
             <position dim=\"1\">{spelling}</position><intensity>{spelling}</intensity>\
             <quality dim=\"0\">{spelling}</quality><quality dim=\"1\">{spelling}</quality>\
             <overallquality>{spelling}</overallquality><charge>2</charge>\
             <UserParam type=\"float\" name=\"FWHM\" value=\"{spelling}\"/>\
             <UserParam type=\"float\" name=\"probe_float\" value=\"{spelling}\"/>\
             <UserParam type=\"floatList\" name=\"probe_floatlist\" value=\"[{spelling},1.5]\"/>\
             </feature>"
        ))
    };
    for (spelling, expected) in [
        ("inf", f64::INFINITY),
        ("+inf", f64::INFINITY),
        ("infinity", f64::INFINITY),
        ("Infinity", f64::INFINITY),
        ("INF", f64::INFINITY),
        ("-inf", f64::NEG_INFINITY),
        ("NaN", f64::NAN),
        ("nan", f64::NAN),
        ("NAN", f64::NAN),
        ("-nan", f64::NAN),
    ] {
        let map = featurexml::read(document(spelling).as_slice())
            .unwrap_or_else(|error| panic!("{spelling}: {error}"));
        let f = &map.features[0];
        let narrow = f64::from(expected as f32);
        assert!(alike(f.rt, expected), "{spelling} rt");
        assert!(alike(f.mz, expected), "{spelling} mz");
        for (name, value) in [
            ("intensity", f.intensity),
            ("quality", f.quality),
            ("quality_rt", f.quality_rt),
            ("quality_mz", f.quality_mz),
            ("width", f.width),
        ] {
            assert!(alike(f64::from(value), narrow), "{spelling} {name}");
        }
        for name in ["FWHM", "probe_float"] {
            assert!(
                alike(f.metadata[name].as_f64().unwrap(), expected),
                "{spelling} {name}"
            );
        }
        let list = f.metadata["probe_floatlist"].as_float_list().unwrap();
        assert!(alike(list[0], expected), "{spelling} list");
        assert_eq!(list[1], 1.5, "{spelling} list");
    }
    for refused in ["inf.0", "1e999"] {
        assert!(
            featurexml::read(document(refused).as_slice()).is_err(),
            "{refused}"
        );
    }
}

/// **A non-finite hull point is the one place this port stays stricter.**
///
/// `ConvexHull2D::setHullPoints` validates nothing
/// (`ConvexHull2D.cpp:119-123`), so the Release build writes and reads
/// `NONFINITE_HULL_RELEASE`, which is `NONFINITE_RELEASE` with its outline
/// points non-finite too and nothing else changed. This port refuses it,
/// because its hulls and the bounding boxes derived from them rest on finite
/// coordinates; the refusal is a parse error, not a panic, and no ported
/// algorithm produces such a point. See FEATUREXML_SUPPORT.md.
#[test]
fn a_nonfinite_hull_point_is_refused_where_the_release_build_keeps_it() {
    let error = featurexml::read(NONFINITE_HULL_RELEASE).unwrap_err();
    assert!(
        format!("{error}").contains("hull point coordinates must be finite"),
        "{error}"
    );
    // The hposition spelling of the same outline is refused for the same
    // reason, and so is a hull point in a document this port otherwise reads.
    for body in [
        "<feature id=\"f_1\"><convexhull nr=\"0\"><pt x=\"inf\" y=\"1\"/></convexhull></feature>",
        "<feature id=\"f_1\"><convexhull nr=\"0\"><hullpoint><hposition dim=\"0\">NaN</hposition>\
         <hposition dim=\"1\">1</hposition></hullpoint></convexhull></feature>",
    ] {
        assert!(featurexml::read(simple(body).as_slice()).is_err(), "{body}");
    }
}

/// **What a reader consumer sees on such a document**, pinned because
/// accepting it here is what lets one reach the rest of the crate.
///
/// `FeatureMap::ranges` — the port's `FeatureMap::updateRanges` — refuses a
/// non-finite value, so the map reads but its ranges are a checked error, and
/// `FileInfo` reports it as one. The Release `FileInfo` on the same document
/// exits 0 and prints `retention time: -inf .. inf sec (inf min)`,
/// `mass-to-charge: -inf .. inf`, `intensity: -inf .. inf` and
/// `Total ion current in features: nan`
/// (`../oracle/featurexml-inf/results/fileinfo_nonfinite.out`, the pinned
/// Release install on `ibminode06`); this port exits 6 with
/// `Invalid parameter: invalid value: range value must be finite`.
///
/// The point of the test is the **absence of a panic**: a document a caller
/// did not write reaches the kernel's finite invariants and is refused there,
/// never aborts. Making the ranges themselves non-finite is a kernel change
/// with its own evidence, not a featureXML one.
#[test]
fn a_nonfinite_map_reads_and_its_ranges_are_a_checked_error() {
    let map = featurexml::read(NONFINITE_RELEASE).unwrap();
    let error = map.ranges().unwrap_err();
    assert!(format!("{error}").contains("must be finite"), "{error}");
    // A finite map of the same shape still answers.
    let finite = featurexml::read(SOURCE).unwrap();
    assert!(finite.ranges().is_ok());
}
