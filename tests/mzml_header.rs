#![cfg(feature = "mzml")]
use openms::{format::mzml, metadata::*};
use std::io::Cursor;

#[test]
fn original_header_class_test_literals() {
    let e = mzml::read(Cursor::new(include_bytes!("data/mzml_header/header.mzML"))).unwrap();
    let s = &e.settings;
    assert_eq!(s.document.identifier, "document_accession");
    assert_eq!(s.fraction_identifier, "the_best_fraction_ever");
    assert_eq!(s.date_time.get(), "2007-06-27 15:23:45");
    assert_eq!(s.contacts.len(), 2);
    assert_eq!(s.contacts[0].first_name, "William");
    assert_eq!(s.contacts[0].last_name, "Pennington");
    assert_eq!(s.contacts[0].email, "wpennington@higglesworth.edu");
    assert_eq!(s.contacts[1].first_name, "Guybrush");
    assert_eq!(s.contacts[1].last_name, "Threepwood");
    assert_eq!(s.source_files.len(), 5);
    let f = &s.source_files[0];
    assert_eq!(f.name, "tiny1.RAW");
    assert_eq!(f.path, "file:///F:/data/Exp01");
    assert_eq!(f.checksum, "71be39fb2700ab2f3c8b2234b91274968b6899b1");
    assert_eq!(f.checksum_type, ChecksumType::Sha1);
    assert_eq!(f.file_type, "Thermo RAW format");
    assert_eq!(f.native_id_type, "multiple peak list nativeID format");
    let a = &s.sample;
    assert_eq!(
        (&a.name, &a.number),
        (&"Sample1".to_owned(), &"5".to_owned())
    );
    assert_eq!((a.mass, a.volume, a.concentration), (11.7, 3.1, 5.5));
    assert_eq!(a.state, SampleState::Suspension);
    assert_eq!(a.metadata["cellular quality"], MetaValue::from("11.11"));
    assert_eq!(
        a.metadata["brenda source tissue"],
        MetaValue::from("cardiac muscle")
    );
    assert_eq!(
        a.metadata["GO cellular component"],
        MetaValue::from("nucleus")
    );
    let i = &s.instrument;
    assert_eq!(
        (&i.name, &i.customizations),
        (&"LCQ Deca".to_owned(), &"Umbau".to_owned())
    );
    assert_eq!(i.ion_optics, IonOpticsType::MagneticDeflection);
    assert_eq!(
        (
            i.ion_sources.len(),
            i.mass_analyzers.len(),
            i.ion_detectors.len()
        ),
        (2, 2, 2)
    );
    assert_eq!((i.ion_sources[0].order, i.ion_sources[1].order), (101, 102));
    assert_eq!(i.ion_sources[0].ionization_method, IonizationMethod::Esi);
    assert_eq!(i.ion_sources[1].ionization_method, IonizationMethod::Fab);
    let a = &i.mass_analyzers[0];
    assert_eq!(a.analyzer_type, AnalyzerType::PaulIonTrap);
    assert_eq!(
        (
            a.accuracy,
            a.magnetic_field_strength,
            a.tof_total_path_length
        ),
        (10.5, 14.56, 11.1)
    );
    assert_eq!(i.mass_analyzers[1].analyzer_type, AnalyzerType::Lit);
    assert_eq!(
        i.ion_detectors[0].detector_type,
        DetectorType::ElectronMultiplier
    );
    assert_eq!(
        i.ion_detectors[0].acquisition_mode,
        DetectorAcquisitionMode::Tdc
    );
    assert_eq!(i.ion_detectors[0].resolution, 5.1);
    assert_eq!(i.ion_detectors[1].resolution, 6.1);
    assert_eq!(i.software.name, "Bioworks");
    assert_eq!(i.software.version, "3.3.1 sp1");
}

#[test]
fn source_header_writer_roundtrip_and_processing_seconds() {
    use openms::{MSExperiment, MSSpectrum, Peak1D};
    use std::sync::Arc;
    let mut e = mzml::read(Cursor::new(include_bytes!("data/mzml_header/header.mzML"))).unwrap();
    let method = Arc::new(DataProcessing {
        software: Software {
            name: "my & custom \"writer\"".into(),
            version: "v\"&1".into(),
            cv_terms: Default::default(),
        },
        completion_time: Some("2001-02-03 04:05:37".parse().unwrap()),
        ..Default::default()
    });
    let mut spectrum = MSSpectrum::new();
    spectrum.native_id = "scan=1".into();
    spectrum.peaks.push(Peak1D::new(100.0, 5.0));
    spectrum.data_processing.push(method);
    e.spectra.push(spectrum);
    let mut xml = Vec::new();
    mzml::write(&mut xml, &e).unwrap();
    let back = mzml::read(Cursor::new(&xml)).unwrap();
    assert_eq!(back, e);
    let text = std::str::from_utf8(&xml).unwrap();
    assert!(text.contains("04:05:37"));
    assert!(text.contains("v&quot;&amp;1"));
    assert!(text.contains("openms-rust:empty-processing-actions"));
    let empty = MSExperiment {
        spectra: vec![MSSpectrum {
            native_id: "scan=1".into(),
            ..Default::default()
        }],
        ..Default::default()
    };
    xml.clear();
    mzml::write(&mut xml, &empty).unwrap();
    assert_eq!(mzml::read(Cursor::new(&xml)).unwrap(), empty);
}

fn wrapper(header: &str, body: &str) -> String {
    format!(
        r#"<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">{header}<run id="r">{body}</run></mzML>"#
    )
}
const REFERENCES: &str = r#"
<fileDescription><fileContent/><sourceFileList count="1"><sourceFile id="sf" name="x.raw" location="file:///data"/></sourceFileList></fileDescription>
<softwareList count="1"><software id="sw" version="1"><cvParam accession="MS:1000799" name="custom unreleased software tool" value="tool"/></software></softwareList>
<dataProcessingList count="2">
<dataProcessing id="dp"><processingMethod order="9" softwareRef="sw"><cvParam accession="MS:1000035" name="peak picking"/></processingMethod><processingMethod order="1" softwareRef="sw"><cvParam accession="MS:1000592" name="smoothing"/></processingMethod></dataProcessing>
<dataProcessing id="other"><processingMethod order="0" softwareRef="sw"><cvParam accession="MS:1000033" name="deisotoping"/></processingMethod></dataProcessing>
</dataProcessingList>"#;

#[test]
fn reference_attachment_reuses_arcs_and_keeps_method_encounter_order() {
    let xml = wrapper(
        REFERENCES,
        r#"<spectrumList count="3" defaultDataProcessingRef="dp"><spectrum id="a" defaultArrayLength="0" sourceFileRef="sf"/><spectrum id="b" defaultArrayLength="0"/><spectrum id="c" defaultArrayLength="0" dataProcessingRef="other"/></spectrumList>"#,
    );
    let e = mzml::read(Cursor::new(xml)).unwrap();
    assert_eq!(e.spectra[0].source_file.name, "x.raw");
    assert_eq!(e.spectra[0].data_processing.len(), 2);
    assert!(
        e.spectra[0].data_processing[0]
            .actions
            .contains(&ProcessingAction::PeakPicking)
    );
    assert!(std::sync::Arc::ptr_eq(
        &e.spectra[0].data_processing[1],
        &e.spectra[1].data_processing[1]
    ));
    assert!(
        e.spectra[2].data_processing[0]
            .actions
            .contains(&ProcessingAction::Deisotoping)
    );
}

#[test]
fn metadata_only_stops_before_bad_count_binary_and_tail() {
    let xml = format!(
        r#"<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">{REFERENCES}<run><spectrumList defaultDataProcessingRef="dp" count="invalid"><broken \x00"#
    );
    let settings = mzml::read_metadata(Cursor::new(xml.as_bytes())).unwrap();
    assert_eq!(settings.source_files.len(), 1);
    let mut options = mzml::LoadOptions::default();
    options.scientific.metadata_only = true;
    options.max_selection_work = 0;
    let e =
        mzml::read_with_load_options(Cursor::new(xml.as_bytes()), &options, &Default::default())
            .unwrap();
    assert!(e.spectra.is_empty());
    assert_eq!(e.settings, settings);
    assert!(
        mzml::read_metadata(Cursor::new(xml.replace(
            "defaultDataProcessingRef=\"dp\"",
            "defaultDataProcessingRef=\"missing\""
        )))
        .is_err()
    );
}

#[test]
fn header_references_and_counts_fail_without_partial_result() {
    for text in [
        REFERENCES.replace("softwareRef=\"sw\"", "softwareRef=\"absent\""),
        REFERENCES.replace("id=\"other\"", "id=\"dp\""),
    ] {
        assert!(mzml::read(Cursor::new(wrapper(&text, ""))).is_err());
    }
    // A header list `count` that disagrees with the number of children is
    // advisory on reading, as in source: the upstream TOPP fixture
    // DTAExtractor_1_input.mzML declares softwareList count="5" with four
    // entries and C++ loads it. Unresolved references above still fail.
    let miscounted = REFERENCES.replace("softwareList count=\"1\"", "softwareList count=\"2\"");
    assert!(mzml::read(Cursor::new(wrapper(&miscounted, ""))).is_ok());
    let mut limits = mzml::ReadOptions {
        max_param_bytes: 1,
        ..Default::default()
    };
    assert!(mzml::read_with_options(Cursor::new(wrapper(REFERENCES, "")), &limits).is_err());
    limits.max_param_bytes = usize::MAX;
    limits.max_total_params = 1;
    assert!(mzml::read_with_options(Cursor::new(wrapper(REFERENCES, "")), &limits).is_err());
}

#[test]
fn array_descriptions_references_and_selection_survive_both_writers() {
    use openms::kernel::{DataArray, MSExperiment, MSSpectrum, Peak1D};
    use std::sync::Arc;
    let dp = Arc::new(DataProcessing {
        software: Software {
            name: "test tool".into(),
            ..Default::default()
        },
        actions: [ProcessingAction::PeakPicking].into(),
        ..Default::default()
    });
    let metadata: MetaInfo = [
        ("units".into(), MetaValue::from("typed")),
        ("integer".into(), MetaValue::from(4i64)),
    ]
    .into();
    let mut s = MSSpectrum {
        native_id: "scan=1".into(),
        data_processing: vec![dp.clone()],
        peaks: vec![
            Peak1D::new(300., 30.),
            Peak1D::new(100., 10.),
            Peak1D::new(200., 20.),
        ],
        ..Default::default()
    };
    s.float_data_arrays.push(DataArray {
        name: "float aux".into(),
        data: vec![3., 1., 2.],
        metadata: metadata.clone(),
        data_processing: vec![dp.clone()],
    });
    s.integer_data_arrays.push(DataArray {
        name: "int aux".into(),
        data: vec![30, 10, 20],
        metadata: metadata.clone(),
        data_processing: vec![dp.clone()],
    });
    s.string_data_arrays.push(DataArray {
        name: "text aux".into(),
        data: vec!["c".into(), "a".into(), "b".into()],
        metadata,
        data_processing: vec![dp],
    });
    let e = MSExperiment {
        spectra: vec![s],
        ..Default::default()
    };
    for numpress in [false, true] {
        let mut bytes = Vec::new();
        if numpress {
            mzml::write_with_numpress(&mut bytes, &e, &Default::default()).unwrap();
        } else {
            mzml::write(&mut bytes, &e).unwrap();
        }
        let plain = mzml::read(Cursor::new(&bytes)).unwrap();
        assert_eq!(plain, e);
        assert!(Arc::ptr_eq(
            &plain.spectra[0].data_processing[0],
            &plain.spectra[0].float_data_arrays[0].data_processing[0]
        ));
        let mut options = mzml::LoadOptions::default();
        options
            .scientific
            .set_mz_range(openms::kernel::NumericRange {
                min: 100.,
                max: 300.,
            });
        let selected =
            mzml::read_with_load_options(Cursor::new(&bytes), &options, &Default::default())
                .unwrap();
        let s = &selected.spectra[0];
        assert_eq!(
            s.peaks.iter().map(|p| p.mz).collect::<Vec<_>>(),
            [100., 200.]
        );
        assert_eq!(s.float_data_arrays[0].data, [1., 2.]);
        assert_eq!(s.integer_data_arrays[0].data, [10, 20]);
        assert_eq!(s.string_data_arrays[0].data, ["a", "b"]);
        assert_eq!(
            s.float_data_arrays[0].metadata,
            e.spectra[0].float_data_arrays[0].metadata
        );
        assert!(Arc::ptr_eq(
            &s.data_processing[0],
            &s.float_data_arrays[0].data_processing[0]
        ));
    }
}

#[test]
fn software_recognized_metadata_name_falls_back_without_wrong_path_or_loss() {
    use openms::{MSExperiment, MSSpectrum};
    use std::sync::Arc;
    let mut software = Software {
        name: "test".into(),
        version: "<one & two>".into(),
        ..Default::default()
    };
    software
        .cv_terms
        .metadata
        .insert("completion time".into(), "kept as ordinary metadata".into());
    let e = MSExperiment {
        spectra: vec![MSSpectrum {
            native_id: "scan=1".into(),
            data_processing: vec![Arc::new(DataProcessing {
                software,
                ..Default::default()
            })],
            ..Default::default()
        }],
        ..Default::default()
    };
    let mut bytes = Vec::new();
    mzml::write(&mut bytes, &e).unwrap();
    assert_eq!(mzml::read(Cursor::new(&bytes)).unwrap(), e);
    assert!(
        std::str::from_utf8(&bytes)
            .unwrap()
            .contains("<userParam name=\"completion time\"")
    );
}

#[test]
fn fictional_markers_require_exact_payload_and_real_source_conversion_survives() {
    use openms::{MSExperiment, MSSpectrum};
    let e = MSExperiment {
        spectra: vec![MSSpectrum {
            native_id: "scan=1".into(),
            ..Default::default()
        }],
        ..Default::default()
    };
    let mut bytes = Vec::new();
    mzml::write(&mut bytes, &e).unwrap();
    let text = String::from_utf8(bytes).unwrap();
    assert!(
        mzml::read(Cursor::new(
            text.replace("value=\"1\"", "value=\"different\"")
        ))
        .is_err()
    );
    assert!(mzml::read(Cursor::new(text.replace("MS:1000544", "MS:1000035"))).is_err());
    let real=text.replace("<userParam name=\"openms-rust:empty-processing-history\" type=\"xsd:string\" value=\"1\"/>\n","");
    let parsed = mzml::read(Cursor::new(real)).unwrap();
    assert_eq!(parsed.spectra[0].data_processing.len(), 1);
    assert!(
        parsed.spectra[0].data_processing[0]
            .actions
            .contains(&ProcessingAction::ConversionMzML)
    );
    let mut bad = e.clone();
    bad.settings
        .sample
        .metadata
        .insert("comment".into(), "collision".into());
    for numpress in [false, true] {
        let mut output = b"existing".to_vec();
        let result = if numpress {
            mzml::write_with_numpress(&mut output, &bad, &Default::default()).map(|_| ())
        } else {
            mzml::write(&mut output, &bad)
        };
        assert!(result.is_err());
        assert_eq!(output, b"existing");
    }
}

#[test]
fn source_numeric_whitespace_and_plus_rules_apply_to_instrument_and_sample() {
    let source = std::str::from_utf8(include_bytes!("data/mzml_header/header.mzML")).unwrap();
    let xml = source
        .replace("value=\"14.56\"", "value=\" +14.56 \"")
        .replace("value=\"11.7\"", "value=\" +11.7 \"")
        .replace("value=\"5.1\"", "value=\"+-5.1\"");
    let e = mzml::read(Cursor::new(xml)).unwrap();
    assert_eq!(e.settings.sample.mass, 11.7);
    assert_eq!(
        e.settings.instrument.mass_analyzers[0].magnetic_field_strength,
        14.56
    );
    assert_eq!(e.settings.instrument.ion_detectors[0].resolution, -5.1);
}

#[test]
fn fractional_processing_timestamp_survives_both_binary_writers() {
    use openms::{MSExperiment, MSSpectrum};
    use std::sync::Arc;
    let time = "2001-02-03T04:05:37.125".parse().unwrap();
    let e = MSExperiment {
        spectra: vec![MSSpectrum {
            native_id: "scan=1".into(),
            data_processing: vec![Arc::new(DataProcessing {
                completion_time: Some(time),
                ..Default::default()
            })],
            ..Default::default()
        }],
        ..Default::default()
    };
    for numpress in [false, true] {
        let mut bytes = Vec::new();
        if numpress {
            mzml::write_with_numpress(&mut bytes, &e, &Default::default()).unwrap();
        } else {
            mzml::write(&mut bytes, &e).unwrap();
        }
        let back = mzml::read(Cursor::new(&bytes)).unwrap();
        assert_eq!(back, e);
        assert_eq!(
            back.spectra[0].data_processing[0]
                .completion_time
                .unwrap()
                .millisecond(),
            125
        );
        assert!(
            std::str::from_utf8(&bytes)
                .unwrap()
                .contains("2001-02-03T04:05:37.125")
        );
    }
}

#[test]
fn every_source_enum_accession_pair_and_discriminant_roundtrips() {
    let fixture = include_str!("data/mzml_header/instrument_enum_pairs.tsv");
    let mut count = 0;
    for line in fixture.lines().skip(1) {
        let row: Vec<_> = line.split('\t').collect();
        let (owner, kind, member, id) = (row[0], row[1], row[2], row[3]);
        let expected: usize = row[4].parse().unwrap();
        let cv = format!(r#"<cvParam accession="{id}" name="source enum"/>"#);
        let header = format!(
            r#"<instrumentConfigurationList count="1"><instrumentConfiguration id="ic">{}<componentList count="3"><source order="1">{}</source><analyzer order="2">{}</analyzer><detector order="3">{}</detector></componentList></instrumentConfiguration></instrumentConfigurationList>"#,
            if owner == "Instrument" { &cv } else { "" },
            if owner == "IonSource" { &cv } else { "" },
            if owner == "MassAnalyzer" { &cv } else { "" },
            if owner == "IonDetector" { &cv } else { "" }
        );
        let document = format!(
            r#"<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">{header}<run defaultInstrumentConfigurationRef="ic"/></mzML>"#
        );
        let e = mzml::read(Cursor::new(document)).unwrap();
        let i = &e.settings.instrument;
        let actual = match kind {
            "IonOpticsType" => i.ion_optics as usize,
            "InletType" => i.ion_sources[0].inlet_type as usize,
            "IonizationMethod" => i.ion_sources[0].ionization_method as usize,
            "AnalyzerType" => i.mass_analyzers[0].analyzer_type as usize,
            "ReflectronState" => i.mass_analyzers[0].reflectron_state as usize,
            "Type" => i.ion_detectors[0].detector_type as usize,
            "AcquisitionMode" => i.ion_detectors[0].acquisition_mode as usize,
            _ => panic!("unexpected source enum {kind}"),
        };
        assert_eq!(actual, expected, "{owner} {kind} {member} {id}");
        let mut output = Vec::new();
        mzml::write(&mut output, &e).unwrap();
        assert!(
            std::str::from_utf8(&output)
                .unwrap()
                .contains(&format!("accession=\"{id}\"")),
            "{id}"
        );
        assert_eq!(mzml::read(Cursor::new(&output)).unwrap(), e, "{id}");
        count += 1;
    }
    assert_eq!(count, 116);
}

#[test]
fn precursor_source_references_preserve_external_id_and_typed_metadata() {
    let xml = wrapper(
        REFERENCES,
        r#"<spectrumList count="1" defaultDataProcessingRef="dp"><spectrum id="scan=1" defaultArrayLength="0"><precursorList count="1"><precursor sourceFileRef="sf" externalSpectrumID="x&amp;&quot;1"><activation><userParam name="confidence" value="2" type="xsd:integer"/></activation></precursor></precursorList></spectrum></spectrumList>"#,
    );
    let e = mzml::read(Cursor::new(&xml)).unwrap();
    let meta = &e.spectra[0].precursors[0].cv_terms.metadata;
    assert_eq!(meta["source_file_name"], MetaValue::from("x.raw"));
    assert_eq!(meta["source_file_path"], MetaValue::from("file:///data"));
    assert_eq!(meta["external_spectrum_id"], MetaValue::from("x&\"1"));
    assert_eq!(meta["confidence"], MetaValue::from(2_i64));
    for packed in [false, true] {
        let mut bytes = Vec::new();
        if packed {
            mzml::write_with_numpress(&mut bytes, &e, &mzml::NumpressWriteOptions::default())
                .unwrap();
        } else {
            mzml::write(&mut bytes, &e).unwrap();
        }
        assert_eq!(mzml::read(Cursor::new(bytes)).unwrap(), e);
    }
    assert!(
        mzml::read(Cursor::new(
            xml.replace("sourceFileRef=\"sf\"", "sourceFileRef=\"missing\"")
        ))
        .is_err()
    );
}

#[test]
fn primary_metadata_source_merge_order_and_typed_values() {
    fn array(accession: &str, value: &str, kind: &str) -> String {
        format!(
            r#"<binaryDataArray encodedLength="0"><cvParam accession="MS:1000523" name="64-bit float"/><cvParam accession="MS:1000576" name="no compression"/><cvParam accession="{accession}" name="primary" unitCvRef="UO" unitAccession="UO:0000010" unitName="second"/><userParam name="owner" value="{value}" type="{kind}"/><binary/></binaryDataArray>"#
        )
    }
    for chrom in [false, true] {
        let tag = if chrom { "chromatogram" } else { "spectrum" };
        let coordinate = if chrom { "MS:1000595" } else { "MS:1000514" };
        let body = format!(
            r#"<{tag}List count="1"><{tag} id="scan=1" defaultArrayLength="0"><userParam name="owner" value="record"/><binaryDataArrayList count="2">{}{}</binaryDataArrayList></{tag}></{tag}List>"#,
            array("MS:1000515", "intensity", "xsd:string"),
            array(coordinate, "coordinate", "xsd:string")
        );
        let e = mzml::read(Cursor::new(wrapper("", &body))).unwrap();
        let (metadata, expected) = if chrom {
            (&e.chromatograms[0].metadata, "coordinate")
        } else {
            (&e.spectra[0].metadata, "intensity")
        };
        assert_eq!(metadata["owner"].as_str().unwrap(), expected);
        let mut output = Vec::new();
        mzml::write(&mut output, &e).unwrap();
        assert_eq!(mzml::read(Cursor::new(output)).unwrap(), e);
        let typed = body.replace(
            "value=\"coordinate\" type=\"xsd:string\"",
            "value=\"1\" type=\"xsd:integer\"",
        );
        let typed = mzml::read(Cursor::new(wrapper("", &typed))).unwrap();
        if chrom {
            assert_eq!(
                typed.chromatograms[0].metadata["owner"].as_i64().unwrap(),
                1
            );
        } else {
            assert_eq!(
                typed.spectra[0].metadata["owner"].as_str().unwrap(),
                "intensity"
            );
        }
    }
}

#[test]
#[cfg(all(feature = "featurexml", feature = "consensusxml"))]
fn full_processing_datetime_survives_map_xml_and_invalid_states_fail_atomically() {
    use openms::{
        data_structures::DateTime,
        format::{consensusxml, featurexml},
        kernel::{ConsensusMap, FeatureMap},
    };
    let timestamp = DateTime::parse("2001-02-03T04:05:37.123").unwrap();
    let dp = DataProcessing {
        completion_time: Some(timestamp),
        ..Default::default()
    };
    let mut feature = FeatureMap::default();
    feature.data_processing.push(dp.clone());
    let mut consensus = ConsensusMap::default();
    consensus.data_processing.push(dp);
    let mut f = Vec::new();
    featurexml::write(&mut f, &feature).unwrap();
    assert_eq!(
        featurexml::read(f.as_slice()).unwrap().data_processing[0].completion_time,
        Some(timestamp)
    );
    let mut c = Vec::new();
    consensusxml::write(&mut c, &consensus).unwrap();
    assert_eq!(
        consensusxml::read(c.as_slice()).unwrap().data_processing[0].completion_time,
        Some(timestamp)
    );
    feature.data_processing[0].completion_time = Some(DateTime::default());
    consensus.data_processing[0].completion_time = Some(DateTime::default());
    let mut output = b"unchanged".to_vec();
    assert!(featurexml::write(&mut output, &feature).is_err());
    assert_eq!(output, b"unchanged");
    assert!(consensusxml::write(&mut output, &consensus).is_err());
    assert_eq!(output, b"unchanged");
}

#[test]
fn rich_header_and_both_marker_shapes_pass_independent_xsd() {
    use std::{
        io::Write,
        process::{Command, Stdio},
        sync::Arc,
    };
    if Command::new("xmllint").arg("--version").output().is_err() {
        eprintln!("xmllint unavailable; header XSD validation was not executed");
        return;
    }
    let mut e = mzml::read(Cursor::new(include_bytes!("data/mzml_header/header.mzML"))).unwrap();
    let mut s = openms::MSSpectrum::from_peaks(vec![openms::Peak1D::new(100., 1.)]);
    s.native_id = "scan=1".into();
    s.data_processing.push(Arc::new(DataProcessing {
        completion_time: Some("2001-02-03T04:05:37.123".parse().unwrap()),
        ..Default::default()
    }));
    let mut a = openms::kernel::DataArray::new("quality", vec![2.]);
    a.metadata.insert("typed".into(), 2_i64.into());
    a.data_processing = s.data_processing.clone();
    s.float_data_arrays.push(a);
    e.spectra.push(s);
    e.spectra.push(openms::MSSpectrum {
        native_id: "scan=2".into(),
        ..Default::default()
    });
    e.chromatograms.push(openms::MSChromatogram {
        native_id: "trace".into(),
        peaks: vec![openms::ChromatogramPeak::new(2., 3.)],
        ..Default::default()
    });
    // CPP-041: source emits the fraction raw; native quotes/ampersands survive.
    e.settings.fraction_identifier = "fraction \"A\" & B".into();
    e.settings.contacts[0].contact_info = "complete & contact".into();
    e.settings.contacts[0]
        .metadata
        .insert("sample mass".into(), "metadata".into());
    for packed in [false, true] {
        let mut bytes = Vec::new();
        if packed {
            mzml::write_with_numpress(&mut bytes, &e, &Default::default()).unwrap();
        } else {
            mzml::write(&mut bytes, &e).unwrap();
        }
        let mut child = Command::new("xmllint")
            .args(["--nonet", "--noout", "--schema"])
            .arg(std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/mzml_1_10.xsd"))
            .arg("-")
            .stdin(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        child.stdin.take().unwrap().write_all(&bytes).unwrap();
        let result = child.wait_with_output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(mzml::read(Cursor::new(&bytes)).unwrap(), e);
    }
}

#[test]
fn markers_reject_raw_payload_that_semantic_field_parsing_would_discard() {
    let mut e = openms::MSExperiment::default();
    e.spectra.push(openms::MSSpectrum {
        native_id: "scan=1".into(),
        ..Default::default()
    });
    let mut bytes = Vec::new();
    mzml::write(&mut bytes, &e).unwrap();
    let xml = String::from_utf8(bytes).unwrap();
    for changed in [
        xml.replacen("<processingMethod order=", "<processingMethod unrepresented=\"x\" order=",1),
        xml.replacen("</processingMethod>","<userParam name=\"extra\" type=\"xsd:string\" value=\"x\"/></processingMethod>",1),
        xml.replacen("accession=\"MS:1000544\"", "accession=\"MS:1000544\" unitAccession=\"UO:0000010\"",1),
        xml.replacen("</processingMethod>","<ignored/></processingMethod>",1),
        xml.replacen("</software>","<cvParam cvRef=\"MS\" accession=\"MS:1000747\" name=\"completion time\" value=\"2001-02-03T04:05:37\"/></software>",1),
        xml.replacen("</software>","<userParam name=\"extra\" value=\"x\"/></software>",1),
    ] {
        assert!(mzml::read(Cursor::new(changed)).is_err());
    }
    let method = DataProcessing {
        metadata: [("openms-rust:empty-processing-history".into(), "1".into())].into(),
        ..Default::default()
    };
    e.spectra[0]
        .data_processing
        .push(std::sync::Arc::new(method));
    for packed in [false, true] {
        let mut output = b"existing".to_vec();
        let result = if packed {
            mzml::write_with_numpress(&mut output, &e, &Default::default()).map(|_| ())
        } else {
            mzml::write(&mut output, &e)
        };
        assert!(result.is_err());
        assert_eq!(output, b"existing");
    }
}

#[test]
fn fractional_run_timestamp_and_negative_zero_sample_survive_without_stale_raw_override() {
    let mut e = openms::MSExperiment::default();
    e.settings.date_time = "2001-02-03T04:05:37.123".parse().unwrap();
    e.settings.metadata.insert(
        "mzml_start_time_stamp".into(),
        "2001-02-03T04:05:37.456".into(),
    );
    e.settings.sample.mass = -0.;
    let mut output = Vec::new();
    mzml::write(&mut output, &e).unwrap();
    let back = mzml::read(Cursor::new(output)).unwrap();
    assert_eq!(back.settings.date_time, e.settings.date_time);
    assert_eq!(back.settings.sample.mass.to_bits(), (-0_f64).to_bits());
    assert_eq!(
        back.settings.metadata["mzml_start_time_stamp"]
            .as_str()
            .unwrap(),
        "2001-02-03T04:05:37.123"
    );
}

#[test]
fn additional_instrument_ids_cannot_collide_with_generated_header_identity() {
    for id in ["so_default", "MS", "run", " padded "] {
        let mut e = openms::MSExperiment::default();
        e.settings
            .instrument_configurations
            .insert(id.into(), Instrument::default());
        let mut output = b"existing".to_vec();
        assert!(mzml::write(&mut output, &e).is_err(), "{id}");
        assert_eq!(output, b"existing");
    }
    let mut e = openms::MSExperiment::default();
    e.settings
        .instrument_configurations
        .insert("ic_0".into(), Instrument::default());
    let mut output = Vec::new();
    mzml::write(&mut output, &e).unwrap();
    assert_eq!(mzml::read(Cursor::new(output)).unwrap(), e);
}

#[test]
fn ignored_kind_attributes_and_misplaced_headers_cannot_disappear() {
    let mut e = openms::MSExperiment::default();
    e.spectra.push(openms::MSSpectrum {
        native_id: "scan=1".into(),
        ..Default::default()
    });
    let mut bytes = Vec::new();
    mzml::write(&mut bytes, &e).unwrap();
    let xml = String::from_utf8(bytes).unwrap();
    for changed in [
        xml.replacen(
            "accession=\"MS:1000544\"",
            "accession=\"MS:1000544\" type=\"hidden\"",
            1,
        ),
        xml.replacen(
            "name=\"openms-rust:empty-processing-history\"",
            "name=\"openms-rust:empty-processing-history\" accession=\"hidden\"",
            1,
        ),
        xml.replacen(
            "name=\"openms-rust:empty-processing-history\"",
            "name=\"openms-rust:empty-processing-history\" cvRef=\"hidden\"",
            1,
        ),
    ] {
        assert!(mzml::read(Cursor::new(changed)).is_err());
    }
    for item in [
        r#"<sourceFile id="sf" name="x.raw" location="file:///d"/>"#,
        r#"<instrumentConfiguration id="ic"/>"#,
        r#"<softwareList count="0"/>"#,
        r#"<dataProcessing id="dp"/>"#,
    ] {
        assert!(
            mzml::read(Cursor::new(wrapper("", item))).is_err(),
            "{item}"
        );
    }
}

#[test]
fn chromatogram_source_file_is_rejected_instead_of_emitting_an_invalid_attribute() {
    let mut e = openms::MSExperiment::default();
    e.chromatograms.push(openms::MSChromatogram {
        source_file: SourceFile {
            name: "x.raw".into(),
            ..Default::default()
        },
        ..Default::default()
    });
    for packed in [false, true] {
        let mut output = b"existing".to_vec();
        let result = if packed {
            mzml::write_with_numpress(&mut output, &e, &Default::default()).map(|_| ())
        } else {
            mzml::write(&mut output, &e)
        };
        assert!(
            matches!(result,Err(openms::Error::Unsupported(s)) if s.contains("chromatogram source-file"))
        );
        assert_eq!(output, b"existing");
        assert_eq!(e.chromatograms[0].source_file.name, "x.raw");
    }
    assert!(mzml::read(Cursor::new(wrapper(REFERENCES,r#"<chromatogramList count="1"><chromatogram id="c" defaultArrayLength="0" sourceFileRef="sf"/></chromatogramList>"#))).is_err());
}

#[test]
fn invalid_checksum_text_is_rejected_distinctly_from_xml_escaping() {
    // CPP-041 also affects arbitrary source checksum strings. The native record
    // keeps its existing digest-domain check, rather than accepting a bad SHA1.
    let mut e = openms::MSExperiment::default();
    e.settings.source_files.push(SourceFile {
        checksum: "bad\"&digest".into(),
        checksum_type: ChecksumType::Sha1,
        ..Default::default()
    });
    let mut output = b"existing".to_vec();
    assert!(
        mzml::write(&mut output, &e)
            .unwrap_err()
            .to_string()
            .contains("checksum")
    );
    assert_eq!(output, b"existing");
}
