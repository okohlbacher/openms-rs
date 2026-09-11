// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "mzml-schema")]
use openms::format::mzml::{SchemaValidationOptions, validate_schema_reader};
use std::sync::atomic::{AtomicUsize, Ordering};
const PLAIN: &str = include_str!("data/mzml_validator/MzMLFile_1.mzML");
fn prolog(insert: &str) -> String {
    let n = PLAIN.find("?>").unwrap() + 2;
    format!("{}\n{insert}\n{}", &PLAIN[..n], &PLAIN[n..])
}
#[test]
fn offline_document_and_schema_fetching_is_blocked_with_observed_positive_control() {
    // Separate integration-test process: this test-only global observer cannot
    // influence production or the other integration suites. Unexpected IO gets
    // empty data, never a fallback to the actual file/network resolver.
    static LOADS: AtomicUsize = AtomicUsize::new(0);
    libxml::io::register_input_callback(
        |_| true,
        |_| {
            LOADS.fetch_add(1, Ordering::SeqCst);
            Some(Vec::new())
        },
    );
    let o = SchemaValidationOptions::default();
    for path in [
        "file:///tmp/openms-schema-must-not-read",
        "http://127.0.0.1:9/no-network",
        "audit:///no-read",
    ] {
        let hint = PLAIN.replace("http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd", path);
        assert_ne!(hint, PLAIN);
        assert!(
            validate_schema_reader(hint.as_bytes(), &o)
                .unwrap()
                .is_valid()
        );
        let pi = prolog(&format!("<?xml-stylesheet href='{path}' type='text/xsl'?>"));
        assert!(
            validate_schema_reader(pi.as_bytes(), &o)
                .unwrap()
                .is_valid()
        );
        for dtd in [
            format!("<!DOCTYPE mzML SYSTEM '{path}'>"),
            format!("<!DOCTYPE mzML [<!ENTITY external SYSTEM '{path}'>]>"),
            format!("<!DOCTYPE mzML [<!ENTITY % p SYSTEM '{path}'>%p;]>"),
        ] {
            let x = prolog(&dtd);
            assert!(
                validate_schema_reader(x.as_bytes(), &o)
                    .unwrap_err()
                    .to_string()
                    .contains("DTD")
            );
        }
        let xinclude = PLAIN.replacen("<fileDescription>", &format!("<fileDescription><xi:include xmlns:xi='http://www.w3.org/2001/XInclude' href='{path}'/>"), 1);
        assert!(
            !validate_schema_reader(xinclude.as_bytes(), &o)
                .unwrap()
                .is_valid()
        );
    }
    assert_eq!(LOADS.load(Ordering::SeqCst), 0);
    // Positive control uses the binding directly, outside the public fixed-schema
    // API. This proves an actual external grammar load is observed and intercepted.
    let raw = br#"<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema"><xs:include schemaLocation="audit:///positive-control"/></xs:schema>"#;
    let mut p = libxml::schemas::SchemaParserContext::from_buffer(raw);
    assert!(libxml::schemas::SchemaValidationContext::from_parser(&mut p).is_err());
    assert!(LOADS.load(Ordering::SeqCst) > 0);
}
