// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! XSD validation against the bundled source schemas and against a caller's
//! own: `XMLValidator_test.cpp`, `XMLFile_test.cpp`, and the `isValid`
//! assertions of the `*File_test.cpp` class tests whose fixtures need no
//! reader. Every expected verdict is a C++ test literal or a retained C++
//! output; `tests/data/xml_schema_provenance.json` names each one.
//!
//! The per-format entry points, and the source tests that validate what a
//! writer stored, are in `tests/xml_schema_formats.rs`, because they need the
//! formats' own features.
#![cfg(feature = "xml-schema")]

use openms::Error;
use openms::format::xml_schema::{
    SchemaKind, SchemaValidationOptions, SchemaValidationReport, validate, validate_against,
    validate_reader, validate_reader_against, validate_with_options,
};
use std::path::{Path, PathBuf};

fn data(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(relative)
}

fn read(relative: &str) -> String {
    std::fs::read_to_string(data(relative)).unwrap()
}

fn report(schema: SchemaKind, text: &str) -> SchemaValidationReport {
    validate_reader(schema, text.as_bytes(), &SchemaValidationOptions::default()).unwrap()
}

// ===========================================================================
// XMLValidator_test.cpp:37-53, bool isValid(filename, schema, os)
// ===========================================================================

#[test]
fn xml_validator_valid_missing_element_missing_attribute_syntax_and_missing_file() {
    let schema = data("xml_schema/XMLValidator.xsd");
    let valid = data("xml_schema/XMLValidator_valid.xml");

    // TEST_EQUAL(v.isValid("XMLValidator_valid.xml", "XMLValidator.xsd"), true)
    let r = validate_against(&valid, &schema).unwrap();
    assert!(r.is_valid(), "{r:?}");
    assert!(r.diagnostics.is_empty());
    assert_eq!(r.schema, SchemaKind::External);

    // TEST_EQUAL(..."XMLValidator_missing_element.xml"...), false) and
    // TEST_EQUAL(..."XMLValidator_missing_attribute.xml"...), false). The
    // fixture lacks exactly the named element or attribute, so libxml2's
    // message names it; its wording is libxml2's, not Xerces'.
    for (fixture, missing) in [
        (
            "xml_schema/XMLValidator_missing_element.xml",
            "requiredElement",
        ),
        (
            "xml_schema/XMLValidator_missing_attribute.xml",
            "requiredAttribute",
        ),
    ] {
        let r = validate_against(data(fixture), &schema).unwrap();
        assert!(!r.is_valid(), "{fixture}");
        assert!(
            r.diagnostics.iter().any(|d| d.message.contains(missing)),
            "{fixture}: {:?}",
            r.diagnostics
        );
        // XMLValidator::logError_ reports a line; so does libxml2.
        assert!(r.diagnostics.iter().all(|d| d.line == Some(1)), "{r:?}");
    }

    // TEST_EQUAL(..."XMLValidator_syntax.xml"...), false): the source reports
    // Xerces' fatal error and returns false; here a document that is not
    // well-formed is an Error::Parse, never a report.
    assert!(matches!(
        validate_against(data("xml_schema/XMLValidator_syntax.xml"), &schema),
        Err(Error::Parse { .. })
    ));

    // "check vaild fail again to make sure internal states are ok"
    assert!(validate_against(&valid, &schema).unwrap().is_valid());

    // TEST_EXCEPTION(Exception::FileNotFound, v.isValid("this_file_does_not_exist.for_sure", ...))
    match validate_against(
        data("xml_schema/this_file_does_not_exist.for_sure"),
        &schema,
    ) {
        Err(Error::Io(e)) => assert_eq!(e.kind(), std::io::ErrorKind::NotFound),
        other => panic!("expected a missing-file error, got {other:?}"),
    }
}

#[test]
fn a_caller_schema_must_be_self_contained_and_is_preflighted_like_a_document() {
    let valid = read("xml_schema/XMLValidator_valid.xml");
    let schema = read("xml_schema/XMLValidator.xsd");
    let o = SchemaValidationOptions::default();
    assert!(
        validate_reader_against(valid.as_bytes(), schema.as_bytes(), &o)
            .unwrap()
            .is_valid()
    );
    // Every composition element is refused before libxml2 sees the schema,
    // whatever prefix names the XSD namespace; a same-named element in another
    // namespace is not composition.
    for composition in [
        r#"<xs:include schemaLocation="other.xsd"/>"#,
        r#"<xs:import namespace="urn:x" schemaLocation="http://127.0.0.1:9/x.xsd"/>"#,
        r#"<xs:import namespace="urn:x"/>"#,
        r#"<xs:redefine schemaLocation="file:///etc/passwd"/>"#,
        r#"<xs:override schemaLocation="other.xsd"/>"#,
        r#"<q:include xmlns:q="http://www.w3.org/2001/XMLSchema" schemaLocation="other.xsd"/>"#,
    ] {
        let composed = schema.replacen(
            "<xs:element name=\"Root\">",
            &format!("{composition}<xs:element name=\"Root\">"),
            1,
        );
        assert_ne!(composed, schema);
        assert!(
            matches!(
                validate_reader_against(valid.as_bytes(), composed.as_bytes(), &o),
                Err(Error::Unsupported(_))
            ),
            "{composition}"
        );
    }
    // A DTD in the schema is refused as it is in a document.
    let dtd = schema.replacen(
        "<xs:schema",
        "<!DOCTYPE xs:schema [<!ENTITY e SYSTEM 'file:///etc/passwd'>]>\n<xs:schema",
        1,
    );
    assert!(matches!(
        validate_reader_against(valid.as_bytes(), dtd.as_bytes(), &o),
        Err(Error::Unsupported(_))
    ));
    // A schema libxml2 cannot compile is an error, not a verdict on the document.
    let broken = schema.replacen("type=\"xs:string\"", "type=\"xs:noSuchType\"", 1);
    assert_ne!(broken, schema);
    assert!(matches!(
        validate_reader_against(valid.as_bytes(), broken.as_bytes(), &o),
        Err(Error::InvalidValue(_))
    ));
    // Schema and document share one budget.
    let mut tight = o;
    tight.limits.max_xml_bytes = schema.len() - 1;
    assert!(validate_reader_against(valid.as_bytes(), schema.as_bytes(), &tight).is_err());
}

// ===========================================================================
// XMLFile_test.cpp:45-48, bool isValid(filename, os) without a schema
// ===========================================================================

/// A default-constructed `XMLFile` has no schema, and `isValid` throws
/// `Exception::NotImplemented`. Here the only kind that names no bundled
/// schema is [`SchemaKind::External`], and asking to validate against it is
/// [`Error::InvalidValue`].
#[test]
fn a_kind_that_names_no_bundled_schema_is_refused() {
    assert!(matches!(
        validate(
            SchemaKind::External,
            data("xml_schema/XMLValidator_valid.xml")
        ),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(SchemaKind::External.location(), None);
    assert_eq!(SchemaKind::External.version(), None);
}

/// The schema each source `XMLFile` constructor registers, and the version it
/// passes (the `*File.cpp` constructors cited in the provenance).
#[test]
fn every_bundled_kind_carries_its_source_location_and_version() {
    for (kind, location, version) in [
        (SchemaKind::MzML, "/SCHEMAS/mzML_1_10.xsd", "1.1.0"),
        (
            SchemaKind::IndexedMzML,
            "/SCHEMAS/mzML_idx_1_10.xsd",
            "1.1.0",
        ),
        (SchemaKind::FeatureXML, "/SCHEMAS/FeatureXML_1_9.xsd", "1.9"),
        (
            SchemaKind::ConsensusXML,
            "/SCHEMAS/ConsensusXML_1_7.xsd",
            "1.7",
        ),
        (SchemaKind::IdXML, "/SCHEMAS/IdXML_1_5.xsd", "1.5"),
        (SchemaKind::ParamXML, "/SCHEMAS/Param_1_8_0.xsd", "1.8.0"),
        (
            SchemaKind::TransformationXML,
            "/SCHEMAS/TrafoXML_1_1.xsd",
            "1.1",
        ),
        (SchemaKind::MzData, "/SCHEMAS/mzData_1_05.xsd", "1.05"),
        (SchemaKind::MzXML, "/SCHEMAS/mzXML_idx_3.1.xsd", "3.1"),
        (
            SchemaKind::MzIdentML1_0_0,
            "/SCHEMAS/mzIdentML1.0.0.xsd",
            "1.0.0",
        ),
        (
            SchemaKind::MzIdentML1_1_0,
            "/SCHEMAS/mzIdentML1.1.0.xsd",
            "1.1.0",
        ),
        (
            SchemaKind::MzIdentML1_2_0,
            "/SCHEMAS/mzIdentML1.2.0.xsd",
            "1.2.0",
        ),
        (
            SchemaKind::MzIdentML1_3_0,
            "/SCHEMAS/mzIdentML1.3.0.xsd",
            "1.3.0",
        ),
        (SchemaKind::PepXML, "/SCHEMAS/pepXML_v114.xsd", "1.14"),
    ] {
        assert_eq!(kind.location(), Some(location));
        assert_eq!(kind.version(), Some(version));
    }
    for version in ["1.0.0", "1.1.0", "1.2.0", "1.3.0"] {
        assert_eq!(
            SchemaKind::mzidentml(version).and_then(SchemaKind::version),
            Some(version)
        );
    }
    assert_eq!(SchemaKind::mzidentml("1.4.0"), None);
}

// ===========================================================================
// Unchanged source fixtures against their bundled schemas
// ===========================================================================

/// `FeatureXMLFile_test.cpp:388-389`, `ConsensusXMLFile_test.cpp:276-277` and
/// `TransformationXMLFile_test.cpp:41-44`: the unchanged fixtures and the
/// verdict the source asserts for each.
#[test]
fn source_fixtures_validate_as_the_class_tests_assert() {
    for (kind, fixture) in [
        (SchemaKind::FeatureXML, "featurexml_source_1.featureXML"),
        (
            SchemaKind::FeatureXML,
            "featurexml_source_options.featureXML",
        ),
        (
            SchemaKind::ConsensusXML,
            "consensusxml/ConsensusXMLFile_1.consensusXML",
        ),
        (
            SchemaKind::ConsensusXML,
            "consensusxml/ConsensusXMLFile_2_options.consensusXML",
        ),
        (
            SchemaKind::TransformationXML,
            "transformation_xml_1.trafoXML",
        ),
        (
            SchemaKind::TransformationXML,
            "transformation_xml_2.trafoXML",
        ),
        (
            SchemaKind::TransformationXML,
            "transformation_xml_4.trafoXML",
        ),
    ] {
        let r = validate(kind, data(fixture)).unwrap();
        assert!(r.is_valid(), "{fixture}: {r:?}");
        assert!(r.diagnostics.is_empty(), "{fixture}");
        assert_eq!(r.schema, kind);
    }
    // TEST_EQUAL(f.isValid("TransformationXMLFile_3.trafoXML"), false): the
    // fixture's end tags do not match, so it is not well-formed; the source's
    // `false` is Error::Parse here.
    assert!(matches!(
        validate(
            SchemaKind::TransformationXML,
            data("transformation_xml_3.trafoXML")
        ),
        Err(Error::Parse { .. })
    ));
}

/// TOPP_FileInfo_14 and TOPP_FileInfo_15 (test-data `topp/CMakeLists.txt`
/// 913-918) run FileInfo `-v` on two mzIdentML 1.1.0 files that differ in one
/// byte. The retained C++ outputs say the first is valid against the 1.1.0
/// schema and the second is not, at line 327, where `post="]"` breaks the
/// residue pattern.
#[test]
fn retained_fileinfo_verdicts_for_mzidentml_1_1_0() {
    let r = validate(
        SchemaKind::MzIdentML1_1_0,
        data("file_info/inputs/FileInfo_14_input.mzid"),
    )
    .unwrap();
    assert!(r.is_valid(), "{r:?}");
    let r = validate(
        SchemaKind::MzIdentML1_1_0,
        data("file_info/inputs/FileInfo_15_input.mzid"),
    )
    .unwrap();
    assert!(!r.is_valid());
    assert!(
        r.diagnostics.iter().any(|d| d.line == Some(327)),
        "{:?}",
        r.diagnostics
    );
    // Against the default 1.3.0 schema the same valid file is invalid: its
    // root is in the 1.1 namespace. This is why the source detects the version.
    assert!(
        !validate(
            SchemaKind::MzIdentML1_3_0,
            data("file_info/inputs/FileInfo_14_input.mzid")
        )
        .unwrap()
        .is_valid()
    );
}

/// `XMLFile::isValid` validates against the one schema its class registered,
/// whatever the document is: a file of another format has no global
/// declaration for its root, so it cannot be valid.
#[test]
fn a_document_of_another_format_is_invalid_not_an_error() {
    for (kind, fixture) in [
        (SchemaKind::ConsensusXML, "featurexml_source_1.featureXML"),
        (
            SchemaKind::FeatureXML,
            "consensusxml/ConsensusXMLFile_1.consensusXML",
        ),
        (SchemaKind::MzData, "transformation_xml_1.trafoXML"),
        (SchemaKind::ParamXML, "MzDataFile_1.mzData"),
        (SchemaKind::IdXML, "transformation_xml_1.trafoXML"),
    ] {
        let r = validate(kind, data(fixture)).unwrap();
        assert!(!r.is_valid(), "{kind:?} {fixture}");
        assert!(!r.diagnostics.is_empty());
    }
}

// ===========================================================================
// Native preflight, shared with the mzML path
// ===========================================================================

#[test]
fn every_bundled_schema_compiles_and_its_violations_are_reports() {
    // A root no schema declares compiles every bundled schema, including the
    // composed mzXML and mzIdentML 1.0.0 ones, and fails only validation.
    for kind in [
        SchemaKind::MzML,
        SchemaKind::IndexedMzML,
        SchemaKind::FeatureXML,
        SchemaKind::ConsensusXML,
        SchemaKind::IdXML,
        SchemaKind::ParamXML,
        SchemaKind::TransformationXML,
        SchemaKind::MzData,
        SchemaKind::MzXML,
        SchemaKind::MzIdentML1_0_0,
        SchemaKind::MzIdentML1_1_0,
        SchemaKind::MzIdentML1_2_0,
        SchemaKind::MzIdentML1_3_0,
        SchemaKind::PepXML,
    ] {
        let r = report(kind, "<notDeclaredAnywhere/>");
        assert!(!r.is_valid(), "{kind:?}");
        assert_eq!(r.schema, kind);
        assert!(
            r.diagnostics
                .iter()
                .all(|d| d.line == Some(1) && d.filename.is_none()),
            "{kind:?}: {:?}",
            r.diagnostics
        );
    }
}

#[test]
fn preflight_limits_encodings_and_dtds_apply_to_every_schema() {
    let text = read("featurexml_source_1.featureXML");
    let base = SchemaValidationOptions::default();
    for field in 0..5 {
        let mut o = base;
        match field {
            0 => o.limits.max_xml_bytes = 10,
            1 => o.limits.max_elements = 10,
            2 => o.limits.max_depth = 2,
            3 => o.limits.max_work = 100,
            4 => o.limits.max_bytes = 100,
            _ => unreachable!(),
        }
        assert!(
            validate_reader(SchemaKind::FeatureXML, text.as_bytes(), &o)
                .unwrap_err()
                .to_string()
                .contains("limit")
        );
    }
    let doctype = text.replacen("<featureMap", "<!DOCTYPE featureMap>\n<featureMap", 1);
    assert!(matches!(
        validate_reader(SchemaKind::FeatureXML, doctype.as_bytes(), &base),
        Err(Error::Unsupported(_))
    ));
    let unbound = text.replacen("<featureList", "<p:featureList", 1);
    assert!(validate_reader(SchemaKind::FeatureXML, unbound.as_bytes(), &base).is_err());
    let utf16 = text.replacen("encoding=\"ISO-8859-1\"", "encoding=\"UTF-16\"", 1);
    assert_ne!(utf16, text);
    let mut bytes = vec![0xff, 0xfe];
    for ch in utf16.encode_utf16() {
        bytes.extend_from_slice(&ch.to_le_bytes());
    }
    assert!(
        validate_reader(SchemaKind::FeatureXML, bytes.as_slice(), &base)
            .unwrap()
            .is_valid()
    );
    // Compressed input is detected by bytes, whatever the name.
    let directory = openms::system::file::TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = directory.path().join("features.unknown");
    let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    std::io::Write::write_all(&mut gzip, text.as_bytes()).unwrap();
    std::fs::write(&path, gzip.finish().unwrap()).unwrap();
    assert!(
        validate_with_options(SchemaKind::FeatureXML, &path, &base)
            .unwrap()
            .is_valid()
    );
}

#[test]
fn schema_violations_are_reports_with_positions() {
    // FeatureXML_1_9.xsd types a feature's intensity as xs:double.
    let text = read("featurexml_source_1.featureXML");
    let at = text.find("<intensity>").unwrap();
    let broken = format!(
        "{}<intensity>not-a-number{}",
        &text[..at],
        &text[at + text[at..].find("</intensity>").unwrap()..]
    );
    let r = report(SchemaKind::FeatureXML, &broken);
    assert!(!r.is_valid());
    let line = text[..at].matches('\n').count() + 1;
    assert!(
        r.diagnostics.iter().any(|d| d.line == Some(line)),
        "line {line}: {:?}",
        r.diagnostics
    );
}

#[test]
fn engine_is_serialised_and_reusable_across_threads() {
    let feature = read("featurexml_source_1.featureXML");
    let mzid = read("file_info/inputs/FileInfo_15_input.mzid");
    let threads: Vec<_> = (0..4)
        .map(|i| {
            let (feature, mzid) = (feature.clone(), mzid.clone());
            std::thread::spawn(move || {
                for _ in 0..3 {
                    if i % 2 == 0 {
                        assert!(report(SchemaKind::FeatureXML, &feature).is_valid());
                    } else {
                        assert!(!report(SchemaKind::MzIdentML1_1_0, &mzid).is_valid());
                    }
                    assert!(!report(SchemaKind::MzXML, &feature).is_valid());
                }
            })
        })
        .collect();
    for t in threads {
        t.join().unwrap();
    }
}
