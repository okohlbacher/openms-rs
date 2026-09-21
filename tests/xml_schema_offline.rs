// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! No bundled-schema validation reads anything but its input: not a schema
//! hint, a stylesheet, an XInclude, a catalog, nor a file an `xs:include`
//! names. A separate integration-test process, so this test-only process-wide
//! observer cannot influence any other suite.
#![cfg(feature = "xml-schema")]

use openms::format::xml_schema::{SchemaKind, SchemaValidationOptions, validate_reader};
use std::sync::atomic::{AtomicUsize, Ordering};

fn read(relative: &str) -> String {
    std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data")
            .join(relative),
    )
    .unwrap()
}

#[test]
fn no_bundled_schema_validation_loads_anything_with_observed_positive_control() {
    // Registered before any validation, so it is consulted for every URL
    // libxml2 considers loading. Unexpected IO gets empty data, never a
    // fallback to the actual file or network resolver.
    static LOADS: AtomicUsize = AtomicUsize::new(0);
    libxml::io::register_input_callback(
        |_| true,
        |_| {
            LOADS.fetch_add(1, Ordering::SeqCst);
            Some(Vec::new())
        },
    );
    let o = SchemaValidationOptions::default();
    let feature = read("featurexml_source_1.featureXML");
    let mzid = read("file_info/inputs/FileInfo_14_input.mzid");
    let mzxml = read("MzXMLFile_5_nested.mzXML");
    for path in [
        "file:///tmp/openms-schema-must-not-read.xsd",
        "http://127.0.0.1:9/no-network.xsd",
        "audit:///no-read",
        "mzXML_3.1_mod.xsd",
        "FuGElightv1.0.0.xsd",
    ] {
        // Every schema-location hint the documents carry, pointed elsewhere.
        let hinted_feature = feature.replacen(
            "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/FeatureXML_1_9.xsd",
            path,
            1,
        );
        assert_ne!(hinted_feature, feature);
        assert!(
            validate_reader(SchemaKind::FeatureXML, hinted_feature.as_bytes(), &o)
                .unwrap()
                .is_valid()
        );
        let hinted_mzid = mzid.replacen(
            "http://psi-pi.googlecode.com/svn/trunk/schema/mzIdentML1.1.0.xsd",
            path,
            1,
        );
        assert_ne!(hinted_mzid, mzid);
        assert!(
            validate_reader(SchemaKind::MzIdentML1_1_0, hinted_mzid.as_bytes(), &o)
                .unwrap()
                .is_valid()
        );
        // The two composed schemas compile without loading their includes.
        for kind in [SchemaKind::MzXML, SchemaKind::MzIdentML1_0_0] {
            let _ = validate_reader(kind, mzxml.as_bytes(), &o).unwrap();
            let _ = validate_reader(kind, hinted_mzid.as_bytes(), &o).unwrap();
        }
        // A stylesheet instruction is inert, a DTD fails before C, and an
        // XInclude is not processed but is an unknown element.
        let end = feature.find("?>").unwrap() + 2;
        let pi = format!(
            "{}\n<?xml-stylesheet href='{path}' type='text/xsl'?>{}",
            &feature[..end],
            &feature[end..]
        );
        assert!(
            validate_reader(SchemaKind::FeatureXML, pi.as_bytes(), &o)
                .unwrap()
                .is_valid()
        );
        let dtd = format!(
            "{}\n<!DOCTYPE featureMap SYSTEM '{path}'>{}",
            &feature[..end],
            &feature[end..]
        );
        assert!(validate_reader(SchemaKind::FeatureXML, dtd.as_bytes(), &o).is_err());
        let xinclude = feature.replacen(
            "<featureList",
            &format!(
                "<xi:include xmlns:xi='http://www.w3.org/2001/XInclude' href='{path}'/><featureList"
            ),
            1,
        );
        assert!(
            !validate_reader(SchemaKind::FeatureXML, xinclude.as_bytes(), &o)
                .unwrap()
                .is_valid()
        );
    }
    assert_eq!(LOADS.load(Ordering::SeqCst), 0);
    // Positive control, outside the public API: an include from a memory
    // schema does reach the observer, so zero above means nothing was loaded.
    let raw = br#"<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema"><xs:include schemaLocation="mzXML_3.1_mod.xsd"/></xs:schema>"#;
    let mut p = libxml::schemas::SchemaParserContext::from_buffer(raw);
    assert!(libxml::schemas::SchemaValidationContext::from_parser(&mut p).is_err());
    assert!(LOADS.load(Ordering::SeqCst) > 0);
}
