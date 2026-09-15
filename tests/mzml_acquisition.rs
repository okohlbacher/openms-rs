// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "mzml")]
use openms::format::mzml::{self, AcquisitionMode, ReadOptions};
use openms::kernel::{MSExperiment, MSSpectrum};
use openms::metadata::{Acquisition, AcquisitionInfo, MetaValue, MetaValueData, ScanWindow};
use std::io::Cursor;

fn document(header: &str, run_attrs: &str, list: &str) -> String {
    format!(
        r#"<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">{header}<run id="r" {run_attrs}><spectrumList count="1"><spectrum id="scan=1" defaultArrayLength="0">{list}</spectrum></spectrumList></run></mzML>"#
    )
}
fn read(xml: &str, mode: AcquisitionMode) -> MSExperiment {
    mzml::read_with_options(
        Cursor::new(xml),
        &ReadOptions {
            acquisition_mode: mode,
            ..Default::default()
        },
    )
    .unwrap()
}
fn info(xml: &str, mode: AcquisitionMode) -> AcquisitionInfo {
    read(xml, mode).spectra.remove(0).acquisition_info
}
fn list(method: &str, count: usize, scans: &str) -> String {
    let term = if method.is_empty() {
        String::new()
    } else {
        format!(r#"<cvParam accession="{method}"/>"#)
    };
    format!(r#"<scanList count="{count}">{term}{scans}</scanList>"#)
}
fn native(info: AcquisitionInfo) -> MSExperiment {
    MSExperiment {
        spectra: vec![MSSpectrum {
            native_id: "scan=1".into(),
            rt: 10.,
            acquisition_info: info,
            ..Default::default()
        }],
        ..Default::default()
    }
}
/// The pinned schema for this output: `mzml::write` is indexed by default,
/// `write_with_numpress` and `write_with_options` are plain.
fn schema_for(bytes: &[u8]) -> std::path::PathBuf {
    let indexed = std::str::from_utf8(bytes).is_ok_and(|text| text.contains("<indexedmzML "));
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(if indexed {
        "tests/data/mzml_writing/mzML_idx_1_10.xsd"
    } else {
        "tests/data/mzml_1_10.xsd"
    })
}
fn encode(e: &MSExperiment, numpress: bool) -> Vec<u8> {
    let mut output = Vec::new();
    if numpress {
        mzml::write_with_numpress(&mut output, e, &Default::default()).unwrap();
    } else {
        mzml::write(&mut output, e).unwrap();
    }
    output
}
#[test]
fn canonical_and_source_modes_distinguish_only_the_implicit_singleton() {
    for method in ["", "MS:1000795"] {
        let xml = document("", "", &list(method, 1, "<scan/>"));
        assert_eq!(
            info(&xml, AcquisitionMode::Canonical),
            AcquisitionInfo::default()
        );
        let source = info(&xml, AcquisitionMode::Source);
        assert_eq!(source.acquisitions, vec![Acquisition::default()]);
        assert_eq!(
            source.method_of_combination,
            if method.is_empty() {
                ""
            } else {
                "no combination"
            }
        );
    }
    for (method, name) in [
        ("MS:1000571", "sum of spectra"),
        ("MS:1000573", "median of spectra"),
        ("MS:1000575", "mean of spectra"),
    ] {
        let xml = document("", "", &list(method, 1, "<scan/>"));
        let a = info(&xml, AcquisitionMode::Canonical);
        assert_eq!(a.method_of_combination, name);
        assert_eq!(a.acquisitions.len(), 1);
    }
    let xml = document(
        "",
        "",
        &list(
            "MS:1000795",
            3,
            r#"<scan/><scan externalSpectrumID="4711"/><scan/>"#,
        ),
    );
    let a = info(&xml, AcquisitionMode::Canonical);
    assert_eq!(
        a.acquisitions
            .iter()
            .map(|s| s.identifier.as_str())
            .collect::<Vec<_>>(),
        ["", "4711", ""]
    );
    assert_eq!(a, info(&xml, AcquisitionMode::Source));
    let xml = document("", "", &list("MS:1000795", 2, "<scan/><scan/>"));
    assert_eq!(info(&xml, AcquisitionMode::Canonical).acquisitions.len(), 2);
}
#[test]
fn source_literal_acquisitions_and_metadata_projection() {
    // Literal assertions from pinned MzMLFile_test.cpp 373–378/430–432,
    // 523–526/569–571 and 631–634; no native-generated numeric expectations.
    let e = read(
        include_str!("data/mzml_acquisition_source_projection.mzML"),
        AcquisitionMode::Canonical,
    );
    assert_eq!(e.spectra.len(), 4);
    for (index, ids, method) in [
        (0, vec!["4711", "4712"], "median of spectra"),
        (1, vec!["0"], "no combination"),
        (2, vec!["4711", "4712"], "median of spectra"),
        (3, vec!["0"], "no combination"),
    ] {
        let a = &e.spectra[index].acquisition_info;
        assert_eq!(a.method_of_combination, method);
        assert_eq!(
            a.acquisitions
                .iter()
                .map(|a| a.identifier.as_str())
                .collect::<Vec<_>>(),
            ids
        );
    }
    let a = &e.spectra[0].acquisition_info;
    assert_eq!(a.metadata["name"].as_str().unwrap(), "acquisition_list");
    assert_eq!(
        a.acquisitions[0].metadata["name"].as_str().unwrap(),
        "acquisition1"
    );
    assert_eq!(
        a.acquisitions[1].metadata["name"].as_str().unwrap(),
        "acquisition2"
    );
    assert_eq!(
        a.acquisitions[0].metadata["source_file_name"]
            .as_str()
            .unwrap(),
        "ac.dta"
    );
    assert_eq!(
        a.acquisitions[0].metadata["source_file_path"]
            .as_str()
            .unwrap(),
        "file:///F:/data/Exp02"
    );
    for compressed in [false, true] {
        assert_eq!(mzml::read(Cursor::new(encode(&e, compressed))).unwrap(), e);
    }
}
#[test]
fn scalar_metadata_all_nine_typed_cv_paths_and_units() {
    let scans = r#"<scan externalSpectrumID="a&amp;&quot;&lt;&#10;b"><cvParam accession="MS:1000927" value="4.5" unitAccession="UO:0000028" unitCvRef="UO" unitName="millisecond"/><cvParam accession="MS:1002082" value="1"/><cvParam accession="MS:1002083" value="2"/><cvParam accession="MS:1002527" value="vendor"/><cvParam accession="MS:1002528"/><cvParam accession="MS:1002892"/><cvParam accession="MS:1003057" value="-7"/><cvParam accession="MS:1003371" value="-25.5"/><cvParam accession="MS:1003394" value="30"/><userParam name="explicit" type="xsd:integer" value="9223372036854775807"/></scan>"#;
    let xml = document("", "", &list("MS:1000795", 1, scans));
    let e = read(&xml, AcquisitionMode::Canonical);
    let a = &e.spectra[0].acquisition_info.acquisitions[0];
    assert_eq!(a.identifier, "a&\"<\nb");
    assert_eq!(a.metadata.len(), 10);
    assert_eq!(a.metadata["MS:1000927"].data(), &MetaValueData::Float(4.5));
    assert_eq!(
        a.metadata["MS:1000927"].unit().unwrap().accession(),
        "UO:0000028"
    );
    assert_eq!(a.metadata["MS:1003057"].data(), &MetaValueData::Integer(-7));
    for key in ["MS:1002528", "MS:1002892"] {
        assert_eq!(a.metadata[key].as_str().unwrap(), "");
    }
    assert_eq!(
        a.metadata["explicit"].data(),
        &MetaValueData::Integer(i64::MAX)
    );
    for compressed in [false, true] {
        assert_eq!(mzml::read(Cursor::new(encode(&e, compressed))).unwrap(), e);
    }
    for value in ["2147483648", "not-a-number"] {
        assert!(
            mzml::read(Cursor::new(
                xml.replace("value=\"-7\"", &format!("value=\"{value}\""))
            ))
            .is_err()
        );
    }
}
#[test]
fn source_writer_normalization_and_first_scan_placement() {
    let mut a = AcquisitionInfo {
        acquisitions: vec![
            Acquisition::default(),
            Acquisition {
                identifier: "b".into(),
                ..Default::default()
            },
        ],
        ..Default::default()
    };
    a.metadata.insert("name".into(), "list".into());
    let mut e = native(a);
    e.spectra[0].instrument_settings.zoom_scan = true;
    e.spectra[0]
        .instrument_settings
        .scan_windows
        .push(ScanWindow {
            begin: 100.,
            end: 200.,
            ..Default::default()
        });
    for compressed in [false, true] {
        let bytes = encode(&e, compressed);
        let xml = String::from_utf8(bytes).unwrap();
        assert_eq!(xml.matches("accession=\"MS:1000016\"").count(), 1);
        assert_eq!(xml.matches("<scanWindowList").count(), 1);
        assert_eq!(xml.matches("accession=\"MS:1000497\"").count(), 2);
        let loaded = read(&xml, AcquisitionMode::Canonical);
        assert_eq!(
            loaded.spectra[0].acquisition_info.method_of_combination,
            "no combination"
        );
        assert_eq!(
            loaded.spectra[0].acquisition_info.acquisitions,
            e.spectra[0].acquisition_info.acquisitions
        );
    }
    let explicit = native(AcquisitionInfo {
        acquisitions: vec![Acquisition::default()],
        ..Default::default()
    });
    let bytes = encode(&explicit, false);
    let xml = String::from_utf8(bytes).unwrap();
    assert_eq!(
        info(&xml, AcquisitionMode::Canonical),
        AcquisitionInfo::default()
    );
    let source = info(&xml, AcquisitionMode::Source);
    assert_eq!(source.method_of_combination, "no combination");
    assert_eq!(source.acquisitions.len(), 1);
}
#[test]
fn header_reference_resolution_and_additional_instrument_writer_are_checked() {
    let h = r#"<fileDescription><sourceFileList count="1"><sourceFile id="f" name="raw &amp; file" location="file:///data"/></sourceFileList></fileDescription><instrumentConfigurationList count="2"><instrumentConfiguration id="a"/><instrumentConfiguration id="b"/></instrumentConfigurationList>"#;
    let xml = document(
        h,
        r#"defaultInstrumentConfigurationRef="a""#,
        &list(
            "MS:1000795",
            2,
            r#"<scan sourceFileRef="f" instrumentConfigurationRef="a"/><scan instrumentConfigurationRef="b"/>"#,
        ),
    );
    let e = read(&xml, AcquisitionMode::Canonical);
    let a = &e.spectra[0].acquisition_info.acquisitions;
    assert_eq!(
        a[0].metadata["source_file_name"].as_str().unwrap(),
        "raw & file"
    );
    assert!(!a[0].metadata.contains_key("instrument_configuration_ref"));
    assert_eq!(
        a[1].metadata["instrument_configuration_ref"]
            .as_str()
            .unwrap(),
        "b"
    );
    for compressed in [false, true] {
        let bytes = encode(&e, compressed);
        let loaded = read(
            std::str::from_utf8(&bytes).unwrap(),
            AcquisitionMode::Canonical,
        );
        assert_eq!(loaded, e);
    }
    for broken in [
        xml.replace("sourceFileRef=\"f\"", "sourceFileRef=\"missing\""),
        xml.replace(
            "instrumentConfigurationRef=\"b\"",
            "instrumentConfigurationRef=\"missing\"",
        ),
        xml.replace("id=\"b\"", "id=\"a\""),
    ] {
        assert!(mzml::read(Cursor::new(broken)).is_err());
    }
    // A header list `count` that disagrees with the number of children is
    // advisory on reading, as in source. The upstream TOPP fixture
    // DTAExtractor_1_input.mzML declares softwareList count="5" with four
    // entries and dataProcessingList count="3" with one, and C++ loads it;
    // rejecting the mismatch made this port unable to read its own reference
    // data. Writing still emits the true count.
    let miscounted = xml.replace("sourceFileList count=\"1\"", "sourceFileList count=\"2\"");
    let loaded = mzml::read(Cursor::new(miscounted)).expect("advisory count must not fail reading");
    assert_eq!(loaded.settings.source_files.len(), 1);
}
#[test]
fn unsupported_combination_list_metadata_and_missing_instrument_fail_before_output() {
    for case in 0..3 {
        let mut a = AcquisitionInfo::default();
        match case {
            0 => a.method_of_combination = "SUM OF SPECTRA".into(),
            1 => {
                a.metadata.insert(
                    "list".into(),
                    MetaValue::new(MetaValueData::IntegerList(vec![1])).unwrap(),
                );
            }
            _ => {
                let mut scan = Acquisition::default();
                scan.metadata
                    .insert("instrument_configuration_ref".into(), "unknown".into());
                a.acquisitions.push(scan);
            }
        }
        let e = native(a);
        let mut out = b"kept".to_vec();
        assert!(mzml::write(&mut out, &e).is_err());
        assert_eq!(out, b"kept");
    }
    let xml = document(
        "",
        "",
        &list(
            "MS:1000795",
            1,
            r#"<scan><userParam name="name" value="a"/><userParam name="name" value="b"/></scan>"#,
        ),
    );
    assert!(mzml::read(Cursor::new(xml)).is_err());
    let xml = document(
        "",
        "",
        &list(
            "MS:1000795",
            1,
            r#"<cvParam accession="MS:1000571"/><scan/>"#,
        ),
    );
    assert!(mzml::read(Cursor::new(xml)).is_err());
}
#[test]
fn descriptor_and_reference_expansion_budgets_apply_to_excluded_records() {
    let h = format!(
        r#"<fileDescription><sourceFileList count="1"><sourceFile id="f" name="{}" location="/"/></sourceFileList></fileDescription>"#,
        "n".repeat(1000)
    );
    let scans = r#"<scan sourceFileRef="f"/>"#.repeat(20);
    let xml = document(&h, "", &list("MS:1000795", 20, &scans));
    let options = ReadOptions {
        max_param_bytes: 50_000,
        ..Default::default()
    };
    // The same bounded descriptors fit; only resolving the long name on every
    // short reference exhausts the budget.
    mzml::read_with_options(
        Cursor::new(xml.replace("sourceFileRef", "externalSpectrumID")),
        &options,
    )
    .unwrap();
    assert!(
        mzml::read_with_options(Cursor::new(&xml), &options)
            .unwrap_err()
            .to_string()
            .contains("parameter bytes")
    );
    let options = ReadOptions {
        max_total_params: 10,
        ..Default::default()
    };
    assert!(mzml::read_with_options(Cursor::new(&xml), &options).is_err());
    let bad = document(
        "",
        "",
        &list(
            "MS:1000795",
            1,
            r#"<scan><cvParam accession="MS:1003057" value="2147483648"/></scan>"#,
        ),
    );
    let load = mzml::LoadOptions {
        skip_spectra: true,
        ..Default::default()
    };
    assert!(
        mzml::read_with_load_options(Cursor::new(bad), &load, &ReadOptions::default()).is_err()
    );
    let group = r#"<referenceableParamGroupList count="1"><referenceableParamGroup id="g"><userParam name="label" value="shared"/></referenceableParamGroup></referenceableParamGroupList>"#;
    let xml = document(
        group,
        "",
        &list(
            "MS:1000795",
            2,
            r#"<scan><referenceableParamGroupRef ref="g"/></scan><scan><referenceableParamGroupRef ref="g"/></scan>"#,
        ),
    );
    let a = info(&xml, AcquisitionMode::Canonical);
    assert_eq!(
        a.acquisitions[0].metadata["label"].as_str().unwrap(),
        "shared"
    );
    assert_eq!(a.acquisitions[0], a.acquisitions[1]);
}

#[test]
fn ordinary_and_numpress_acquisition_output_validate_against_the_source_xsd() {
    use std::{
        io::Write,
        process::{Command, Stdio},
    };
    if Command::new("xmllint").arg("--version").output().is_err() {
        eprintln!("xmllint unavailable; acquisition XSD validation was not executed");
        return;
    }
    let mut e = read(
        include_str!("data/mzml_acquisition_source_projection.mzML"),
        AcquisitionMode::Canonical,
    );
    // A scan with metadata + zoom + windows exposes source writer ordering:
    // zoom must still precede all userParams in the schema's ParamGroupType.
    e.spectra[0].instrument_settings.zoom_scan = true;
    assert!(
        !e.spectra[0].acquisition_info.acquisitions[0]
            .metadata
            .is_empty()
    );
    assert!(!e.spectra[0].instrument_settings.scan_windows.is_empty());
    for compressed in [false, true] {
        let bytes = encode(&e, compressed);
        let mut child = Command::new("xmllint")
            .args(["--nonet", "--noout", "--schema"])
            .arg(schema_for(&bytes))
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
    }
}
#[test]
fn source_mode_reaches_scientific_and_file_loading_and_preserves_rollback() {
    let e = native(AcquisitionInfo::default());
    let xml = encode(&e, false);
    let options = ReadOptions {
        acquisition_mode: AcquisitionMode::Source,
        ..Default::default()
    };
    let expected =
        mzml::read_with_load_options(Cursor::new(&xml), &Default::default(), &options).unwrap();
    assert_eq!(expected.spectra[0].acquisition_info.acquisitions.len(), 1);
    let path = std::env::temp_dir().join(format!(
        "openms-acquisition-{}-{}.mzML",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::write(&path, &xml).unwrap();
    assert_eq!(
        mzml::load_with_options(&path, &Default::default(), &options).unwrap(),
        expected
    );
    assert!(
        mzml::load(&path).unwrap().spectra[0]
            .acquisition_info
            .acquisitions
            .is_empty()
    );
    let mut destination = expected.clone();
    std::fs::write(&path, b"<bad").unwrap();
    assert!(
        mzml::load_into_with_options(&path, &mut destination, &Default::default(), &options)
            .is_err()
    );
    std::fs::remove_file(path).unwrap();
    assert_eq!(destination, expected);
}
#[test]
fn source_file_locations_follow_source_lexical_repairs_without_filesystem_probes() {
    // MzMLHandler.cpp1062–1116 + File.cpp449–456 + PathUtils.h23–27.
    for (name, location, want_name, want_path) in [
        ("name.raw", "", "name.raw", "file://./"),
        ("/data/name.raw", "", "name.raw", "/data"),
        (r"C:\data\name.raw", "", "name.raw", r"C:\data"),
        ("/name.raw", "", "name.raw", ""),
        ("file", "File:///data", "file", "file:///data"),
        ("file", "FILE:///data", "file", "file:///data"),
        ("file", "file:///./data", "file", "file://.//data"),
        ("file", "file:///", "file", "file://"),
        ("file", "file://../data", "file", "file://../data"),
        ("", "", "", ""),
    ] {
        let h = format!(
            r#"<fileDescription><sourceFileList count="1"><sourceFile id="f" name="{name}" location="{location}"/></sourceFileList></fileDescription>"#
        );
        let xml = document(
            &h,
            "",
            &list("MS:1000795", 1, r#"<scan sourceFileRef="f"/>"#),
        );
        let a = info(&xml, AcquisitionMode::Canonical);
        assert_eq!(
            a.acquisitions[0].metadata["source_file_name"]
                .as_str()
                .unwrap(),
            want_name
        );
        assert_eq!(
            a.acquisitions[0].metadata["source_file_path"]
                .as_str()
                .unwrap(),
            want_path
        );
    }
}
