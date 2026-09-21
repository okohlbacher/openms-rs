// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The two bundled schemas that `xs:include` others, `mzXML_idx_3.1.xsd` and
//! `mzIdentML1.0.0.xsd`, are compiled from one composed document so that the
//! engine loads nothing. This checks the composition against libxml2's own
//! include processing: the same schemas compiled from `resources/schemas` on
//! disk, where libxml2 reads each included file itself, must give every
//! document the same verdict and the same diagnostics, line for line.
//!
//! The documents are the unchanged source fixtures, and namespace-shifted
//! copies of them that exercise the grammar deeply because they violate it in
//! many places. One test function, because the direct libxml2 calls here do
//! not take the crate's engine lock.
#![cfg(feature = "xml-schema")]

use libxml::parser::Parser;
use libxml::schemas::{SchemaParserContext, SchemaValidationContext};
use openms::format::xml_schema::{SchemaKind, SchemaValidationOptions, validate_reader};
use std::path::Path;

fn read(relative: &str) -> String {
    std::fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join(relative)).unwrap()
}

/// libxml2 validating against a schema file whose includes it resolves itself.
fn direct(schema: &str, text: &str) -> (bool, Vec<(String, Option<i32>)>) {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("resources/schemas")
        .join(schema);
    let mut parser = SchemaParserContext::from_file(path.to_str().unwrap());
    let mut validator = SchemaValidationContext::from_parser(&mut parser)
        .unwrap_or_else(|e| panic!("{schema} does not compile from disk: {e:?}"));
    assert!(parser.drain_errors().is_empty(), "{schema}");
    let document = Parser::default().parse_string(text).unwrap();
    let (valid, errors) = match validator.validate_document(&document) {
        Ok(()) => (true, validator.drain_errors()),
        Err(errors) => (false, errors),
    };
    (
        valid && errors.is_empty(),
        errors
            .into_iter()
            .map(|e| (e.message.unwrap_or_default(), e.line))
            .collect(),
    )
}

fn composed(kind: SchemaKind, text: &str) -> (bool, Vec<(String, Option<i32>)>) {
    let r = validate_reader(kind, text.as_bytes(), &SchemaValidationOptions::default()).unwrap();
    (
        r.is_valid(),
        r.diagnostics
            .into_iter()
            .map(|d| (d.message, d.line.map(|l| i32::try_from(l).unwrap())))
            .collect(),
    )
}

#[test]
fn composed_schemas_agree_with_libxml2_include_processing() {
    const MZXML_21: &str = "http://sashimi.sourceforge.net/schema_revision/mzXML_2.1";
    const MZXML_31: &str = "http://sashimi.sourceforge.net/schema_revision/mzXML_3.1";
    let mut mzxml = vec![read("tests/data/MzXMLFile_5_nested.mzXML")];
    for fixture in [
        "tests/data/MzXMLFile_1.mzXML",
        "tests/data/MzXMLFile_1_compressed.mzXML",
        "tests/data/MzXMLFile_2_minimal.mzXML",
        "tests/data/MzXMLFile_3_64bit.mzXML",
    ] {
        let text = read(fixture);
        assert!(text.contains(MZXML_21), "{fixture}");
        mzxml.push(text.clone());
        mzxml.push(text.replace(MZXML_21, MZXML_31));
    }
    // A required attribute removed, and an element the grammar lacks.
    let nested = mzxml[0].clone();
    mzxml.push(nested.replacen(" scanCount=\"", " notAnAttribute=\"", 1));
    mzxml.push(nested.replacen("<msRun", "<unknownElement/><msRun", 1));

    const MZID_11: &str = "http://psidev.info/psi/pi/mzIdentML/1.1";
    const MZID_10: &str = "http://psidev.info/psi/pi/mzIdentML/1.0";
    let mut mzid = Vec::new();
    for fixture in [
        "tests/data/file_info/inputs/FileInfo_14_input.mzid",
        "tests/data/file_info/inputs/FileInfo_15_input.mzid",
        "tests/data/mzidentml_whole.mzid",
        "tests/data/mzidentml_msgf_mini.mzid",
        "tests/data/mzidentml_3runs.mzid",
    ] {
        let text = read(fixture);
        mzid.push(text.clone());
        let shifted = text
            .replace(MZID_11, MZID_10)
            .replace("http://psidev.info/psi/pi/mzIdentML/1.2", MZID_10)
            .replace("http://psidev.info/psi/pi/mzIdentML/1.3", MZID_10);
        assert_ne!(shifted, text, "{fixture}");
        mzid.push(shifted);
    }

    let mut compared = 0;
    let mut with_diagnostics = 0;
    for (kind, schema, documents) in [
        (SchemaKind::MzXML, "mzXML_idx_3.1.xsd", &mzxml),
        (SchemaKind::MzIdentML1_0_0, "mzIdentML1.0.0.xsd", &mzid),
    ] {
        for text in documents.iter() {
            let (valid, diagnostics) = composed(kind, text);
            assert_eq!(
                (valid, &diagnostics),
                (direct(schema, text).0, &direct(schema, text).1)
            );
            compared += 1;
            with_diagnostics += usize::from(diagnostics.len() > 1);
        }
    }
    assert_eq!(compared, mzxml.len() + mzid.len());
    // Most documents violate the grammar in more than one place, so the
    // comparison covers far more than the root declaration.
    assert!(with_diagnostics >= 6, "{with_diagnostics}");
}
