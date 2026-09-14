// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "mzml-schema")]

use openms::format::mzml::{
    SchemaKind, SchemaValidationOptions, SchemaValidationReport, validate_schema,
    validate_schema_reader, validate_schema_with_options,
};
use std::{
    io::{BufReader, Write},
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
};
const PLAIN: &str = include_str!("data/mzml_validator/MzMLFile_1.mzML");
const INDEX: &str = include_str!("data/mzml_validator/MzMLFile_4_indexed.mzML");
fn report(text: &str) -> SchemaValidationReport {
    validate_schema_reader(text.as_bytes(), &SchemaValidationOptions::default()).unwrap()
}
fn invalid(text: &str) -> String {
    let r = report(text);
    assert!(!r.is_valid());
    assert!(!r.diagnostics.is_empty());
    assert!(r.diagnostics.iter().any(|d| d.code != 0 && d.domain != 0));
    format!("{:?}", r.diagnostics)
}
fn prolog(text: &str, insert: &str) -> String {
    let end = text.find("?>").unwrap() + 2;
    format!("{}\n{insert}\n{}", &text[..end], &text[end..])
}
fn attribute(text: &str, element: &str, attr: &str, index: usize) -> (usize, usize) {
    let tag = format!("<{element} ");
    let start = text.match_indices(&tag).nth(index).unwrap().0;
    let field = format!("{attr}=\"");
    let value = start + text[start..].find(&field).unwrap() + field.len();
    (value, value + text[value..].find('"').unwrap())
}
fn change_attribute(text: &str, element: &str, attr: &str, index: usize, value: &str) -> String {
    let (start, end) = attribute(text, element, attr, index);
    format!("{}{}{}", &text[..start], value, &text[end..])
}
#[test]
fn unchanged_originals_validate_against_exact_schemas() {
    let a = report(PLAIN);
    assert!(a.is_valid(), "{a:?}");
    assert_eq!(a.schema, SchemaKind::MzML);
    let b = report(INDEX);
    assert!(b.is_valid(), "{b:?}");
    assert_eq!(b.schema, SchemaKind::IndexedMzML);
}
#[test]
fn namespace_root_selection_corrects_four_line_heuristic() {
    let shifted = prolog(INDEX, "<!--one-->\n<!--two-->\n<!--three-->");
    let r = report(&shifted);
    assert!(r.is_valid());
    assert_eq!(r.schema, SchemaKind::IndexedMzML);
    let prefixed = shifted
        .replacen(
            "<indexedmzML ",
            "<p:indexedmzML xmlns:p=\"http://psi.hupo.org/ms/mzml\" ",
            1,
        )
        .replacen("</indexedmzML>", "</p:indexedmzML>", 1);
    assert!(report(&prefixed).is_valid());
    let p = PLAIN
        .replacen(
            "<mzML ",
            "<q:mzML xmlns:q=\"http://psi.hupo.org/ms/mzml\" ",
            1,
        )
        .replacen("</mzML>", "</q:mzML>", 1);
    assert!(report(&p).is_valid());
    let utf8 = prefixed.replacen("ISO-8859-1", "UTF-8", 1);
    assert!(report(&format!("\u{feff}{utf8}")).is_valid());
    let wrong = PLAIN.replacen(
        "xmlns=\"http://psi.hupo.org/ms/mzml\"",
        "xmlns=\"urn:wrong\"",
        1,
    );
    assert!(validate_schema_reader(wrong.as_bytes(), &SchemaValidationOptions::default()).is_err());
}
#[test]
fn real_xsd_required_attributes_and_numeric_facets() {
    let required = PLAIN.replacen(" version=\"1.1.0\"", "", 1);
    assert!(invalid(&required).contains("required"));
    let negative = PLAIN.replacen("<cvList count=\"5\">", "<cvList count=\"-1\">", 1);
    assert!(invalid(&negative).contains("nonNegativeInteger"));
    let numeric = PLAIN.replacen("<cvList count=\"5\">", "<cvList count=\"text\">", 1);
    assert!(invalid(&numeric).contains("nonNegativeInteger"));
}
#[test]
fn real_xsd_order_and_unknown_elements() {
    let start = PLAIN.find("  <cvList").unwrap();
    let end = PLAIN.find("</cvList>").unwrap() + "</cvList>".len();
    let cv = &PLAIN[start..end];
    let removed = format!("{}{}", &PLAIN[..start], &PLAIN[end..]);
    let reordered = removed.replacen("</fileDescription>", &format!("</fileDescription>{cv}"), 1);
    assert!(invalid(&reordered).contains("not expected"));
    assert!(
        invalid(
            &PLAIN
                .replacen("<fileDescription>", "<notMzML>", 1)
                .replacen("</fileDescription>", "</notMzML>", 1)
        )
        .contains("not expected")
    );
}
#[test]
fn identity_duplicate_key_is_rejected_in_both_schemas() {
    for input in [PLAIN, INDEX] {
        let (a, b) = attribute(input, "spectrum", "id", 0);
        let changed = change_attribute(input, "spectrum", "id", 1, &input[a..b]);
        assert_ne!(input, changed);
        let diagnostics = invalid(&changed);
        assert!(
            diagnostics.contains("Duplicate key-sequence"),
            "{diagnostics}"
        );
    }
}
#[test]
fn identity_dangling_keyref_is_rejected_in_both_schemas() {
    for input in [PLAIN, INDEX] {
        let changed = change_attribute(
            input,
            "processingMethod",
            "softwareRef",
            0,
            "does_not_exist",
        );
        assert_ne!(input, changed);
        let diagnostics = invalid(&changed);
        assert!(
            diagnostics.contains("No match found for key-sequence"),
            "{diagnostics}"
        );
    }
}
#[test]
fn schema_is_not_checksum_offset_count_or_cv_semantic_validation() {
    let changed = INDEX.replacen("<fileChecksum>", "<fileChecksum>not-a-digest-", 1);
    assert!(report(&changed).is_valid());
    let wrong_count = PLAIN.replacen("<cvList count=\"5\">", "<cvList count=\"999\">", 1);
    assert!(report(&wrong_count).is_valid());
    let unknown = PLAIN.replacen(
        "accession=\"MS:1000580\"",
        "accession=\"MS:DOES_NOT_EXIST\"",
        1,
    );
    assert!(report(&unknown).is_valid());
}
#[test]
fn cpp054_raw_index_schema_does_not_establish_index_integrity() {
    // Exact one-occurrence source transformation, independently executed in the
    // backend spike. These unchanged shipped XSD selectors cannot resolve idRef.
    let from = "idRef=\"index=19\"";
    assert_eq!(INDEX.matches(from).count(), 1);
    let altered = INDEX.replacen(from, "idRef=\"does_not_exist\"", 1);
    let r = report(&altered);
    assert!(r.is_valid(), "{r:?}");
    assert!(r.diagnostics.is_empty());
}
#[test]
fn root_attributes_are_normalized_and_namespace_errors_precede_engine() {
    let escaped = PLAIN.replacen(
        "xmlns=\"http://psi.hupo.org/ms/mzml\"",
        "xmlns=\"http://psi.hupo.org/ms/mzm&#108;\"",
        1,
    );
    assert!(report(&escaped).is_valid());
    for xml in [
        "<mzML/>",
        "<notmzML xmlns='http://psi.hupo.org/ms/mzml'/>",
        "<p:mzML xmlns='http://psi.hupo.org/ms/mzml'/>",
        "<mzML xmlns=' http://psi.hupo.org/ms/mzml'/>",
        "<xml:mzML xmlns:xml='http://psi.hupo.org/ms/mzml'/>",
    ] {
        assert!(
            validate_schema_reader(xml.as_bytes(), &SchemaValidationOptions::default()).is_err(),
            "{xml}"
        );
    }
    let illegal = PLAIN.replacen(
        "<fileDescription>",
        "<fileDescription xmlns:xml='urn:illegal'>",
        1,
    );
    assert!(
        validate_schema_reader(illegal.as_bytes(), &SchemaValidationOptions::default()).is_err()
    );
}
#[test]
fn utf16_declaration_is_canonicalized_without_losing_unicode_or_standalone() {
    let utf8 = PLAIN
        .replacen("ISO-8859-1", "UTF-8", 1)
        .replacen(
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
            "<?xml version='1.0'\n encoding='UTF-8' standalone='yes'?>",
            1,
        )
        .replacen(
            "</fileDescription>",
            "<!-- café μ -->\n</fileDescription>",
            1,
        );
    assert!(report(&utf8).is_valid());
    for (encoding, little, bom) in [
        ("UTF-16", true, true),
        ("UTF-16", false, true),
        ("UTF-16LE", true, false),
        ("UTF-16BE", false, false),
    ] {
        let xml = utf8.replacen("UTF-8", encoding, 1);
        let mut bytes = Vec::new();
        if bom {
            bytes.extend_from_slice(if little { &[0xff, 0xfe] } else { &[0xfe, 0xff] });
        }
        for ch in xml.encode_utf16() {
            bytes.extend_from_slice(&if little {
                ch.to_le_bytes()
            } else {
                ch.to_be_bytes()
            });
        }
        let r = validate_schema_reader(
            BufReader::with_capacity(1, bytes.as_slice()),
            &SchemaValidationOptions::default(),
        )
        .unwrap();
        assert!(r.is_valid(), "{encoding}: {r:?}");
        bytes.pop();
        assert!(
            validate_schema_reader(bytes.as_slice(), &SchemaValidationOptions::default()).is_err()
        );
    }
    let contradictory = utf8.replacen("UTF-8", "UTF-16", 1);
    assert!(
        validate_schema_reader(
            contradictory.as_bytes(),
            &SchemaValidationOptions::default()
        )
        .is_err()
    );
}
#[test]
fn malformed_xml_and_ascii_latin1_boundary_are_checked_before_schema_result() {
    for xml in [
        "<mzML xmlns='http://psi.hupo.org/ms/mzml'><broken></mzML>",
        "<mzML xmlns='http://psi.hupo.org/ms/mzml'>]]></mzML>",
        "<mzML xmlns='http://psi.hupo.org/ms/mzml' a='1'b='2'/>",
        "<?xml version='1.0' standalone='yes' encoding='UTF-8'?><mzML/>",
        "<mzML xmlns='http://psi.hupo.org/ms/mzml'/><mzML/>",
    ] {
        assert!(
            validate_schema_reader(xml.as_bytes(), &SchemaValidationOptions::default()).is_err(),
            "{xml}"
        );
    }
    assert!(report(&prolog(PLAIN, "<!-- literal <!DOCTYPE is inert -->")).is_valid());
    let trailing = format!("{PLAIN}junk");
    assert!(
        validate_schema_reader(trailing.as_bytes(), &SchemaValidationOptions::default()).is_err()
    );
    let latin1 = PLAIN.replacen("</mzML>", "<!--é--></mzML>", 1);
    assert!(
        validate_schema_reader(latin1.as_bytes(), &SchemaValidationOptions::default()).is_err()
    );
}
#[test]
fn native_preflight_and_post_engine_result_limits_are_explicit() {
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
            validate_schema_reader(PLAIN.as_bytes(), &o)
                .unwrap_err()
                .to_string()
                .contains("limit")
        );
    }
    let wrong = PLAIN.replace("cvRef=\"MS\"", "cvRef=\"not-known\"");
    assert!(!report(&wrong).is_valid());
    let mut o = base;
    o.limits.max_diagnostics = 1;
    assert!(
        validate_schema_reader(wrong.as_bytes(), &o)
            .unwrap_err()
            .to_string()
            .contains("limit")
    );
    o = base;
    o.limits.max_diagnostic_bytes = 1;
    assert!(
        validate_schema_reader(wrong.as_bytes(), &o)
            .unwrap_err()
            .to_string()
            .contains("limit")
    );
    o = base;
    o.limits.max_diagnostics = 0;
    o.limits.max_diagnostic_bytes = 0;
    assert!(
        validate_schema_reader(PLAIN.as_bytes(), &o)
            .unwrap()
            .is_valid()
    );
}
struct Directory(PathBuf);
impl Directory {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let p = std::env::temp_dir().join(format!(
            "openms-schema-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&p).unwrap();
        Self(p)
    }
    fn write(&self, name: &str, bytes: &[u8]) -> PathBuf {
        let p = self.0.join(name);
        std::fs::write(&p, bytes).unwrap();
        p
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
#[test]
fn path_magic_compression_crc_and_decoded_caps() {
    let d = Directory::new();
    let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gzip.write_all(PLAIN.as_bytes()).unwrap();
    let mut gzip = gzip.finish().unwrap();
    let mut bzip = bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::default());
    bzip.write_all(PLAIN.as_bytes()).unwrap();
    let bzip = bzip.finish().unwrap();
    for (name, data) in [
        ("plain.gz", PLAIN.as_bytes()),
        ("gzip.dat", &gzip),
        ("bzip.dat", &bzip),
    ] {
        let p = d.write(name, data);
        assert!(validate_schema(&p).unwrap().is_valid());
        let mut o = SchemaValidationOptions::default();
        o.limits.max_xml_bytes = PLAIN.len() - 1;
        assert!(validate_schema_with_options(&p, &o).is_err());
    }
    let n = gzip.len();
    gzip[n - 8] ^= 1;
    assert!(validate_schema(d.write("corrupt", &gzip)).is_err());
    assert!(validate_schema(d.write("zip", b"PK\x03\x04payload")).is_err());
    assert!(validate_schema(d.0.join("missing")).is_err());
    assert!(validate_schema(d.write("empty", b"")).is_err());
    assert!(validate_schema(d.write("one", b"<")).is_err());
}
#[test]
fn local_contexts_are_reusable_and_independent_between_threads() {
    let changed = change_attribute(PLAIN, "processingMethod", "softwareRef", 0, "missing");
    for _ in 0..3 {
        assert!(!report(&changed).is_valid());
        assert!(report(PLAIN).is_valid());
    }
    let threads: Vec<_> = (0..4)
        .map(|_| {
            std::thread::spawn(|| {
                for _ in 0..3 {
                    assert!(report(PLAIN).is_valid());
                    assert!(report(INDEX).is_valid());
                }
            })
        })
        .collect();
    for t in threads {
        t.join().unwrap();
    }
}
#[test]
fn exact_schema_bytes_are_preserved_and_have_no_external_grammar_imports() {
    use sha1::{Digest, Sha1};
    for (raw, hash) in [
        (
            &include_bytes!("../resources/schemas/mzML_1_10.xsd")[..],
            "f3c7dd7c6d6a8bd1666af89aa4b5105a84ed4593",
        ),
        (
            &include_bytes!("../resources/schemas/mzML_idx_1_10.xsd")[..],
            "8c81ab2cb3d45ab35c2473f9e46ae185bd404784",
        ),
    ] {
        assert_eq!(format!("{:x}", Sha1::digest(raw)), hash);
        for s in [
            "<xs:include",
            "<xs:import",
            "<xs:redefine",
            "<!DOCTYPE",
            "<!ENTITY",
        ] {
            assert!(!raw.windows(s.len()).any(|w| w == s.as_bytes()));
        }
    }
}
#[test]
fn namespaces_reject_reserved_bindings_undefined_prefixes_and_expanded_duplicates() {
    for attr in [
        "xmlns:xml='urn:wrong'",
        "xmlns:p='http://www.w3.org/XML/1998/namespace'",
        "xmlns:xmlns='urn:any'",
        "xmlns:p='http://www.w3.org/2000/xmlns/'",
        "xmlns:p=''",
        "xmlns='http://www.w3.org/XML/1998/namespace'",
        "p:name='x'",
        "xmlns:p='urn:p' p:1name='x'",
        "xmlns:1p='urn:p'",
        "xmlns:p='urn:p' p:a:b='x'",
        "xmlns:p='urn:p' xmlns:q='urn:p' p:a='x' q:a='y'",
    ] {
        let xml = PLAIN.replacen("<fileDescription>", &format!("<fileDescription {attr}>"), 1);
        assert!(
            validate_schema_reader(xml.as_bytes(), &SchemaValidationOptions::default()).is_err(),
            "{attr}"
        );
    }
    for tag in ["p:1name", "p:a:b", "xmlns:foo", "unbound:foo"] {
        let xml = PLAIN.replacen(
            "<fileDescription>",
            &format!("<fileDescription><{tag} xmlns:p='urn:p'/>"),
            1,
        );
        assert!(
            validate_schema_reader(xml.as_bytes(), &SchemaValidationOptions::default()).is_err(),
            "{tag}"
        );
    }
    // A scoped declaration on an empty sibling must disappear before the next.
    let wrong = PLAIN.replacen(
        "<fileDescription>",
        "<fileDescription><a xmlns:p='urn:p'/><p:b/>",
        1,
    );
    assert!(validate_schema_reader(wrong.as_bytes(), &SchemaValidationOptions::default()).is_err());
    let shadow = PLAIN.replacen("<fileDescription>", "<fileDescription xmlns:p='urn:one'><a xmlns:p='urn:two'/><b xmlns:q='urn:one' p:x='1' q:x='2'/>", 1);
    assert!(
        validate_schema_reader(shadow.as_bytes(), &SchemaValidationOptions::default()).is_err()
    );
    // Unprefixed attrs never inherit the default namespace; this is well-formed
    // even though the schema then rejects the unknown attributes/elements.
    let valid_namespaces = PLAIN.replacen(
        "<fileDescription>",
        "<fileDescription xmlns:p='http://psi.hupo.org/ms/mzml' x='1' p:x='2' xml:lang='en'>",
        1,
    );
    assert!(!report(&valid_namespaces).is_valid());
    let xml_binding = PLAIN.replacen(
        "<fileDescription>",
        "<fileDescription xmlns:xml='http://www.w3.org/XML/1998/namespace'>",
        1,
    );
    assert!(report(&xml_binding).is_valid());
}

#[test]
fn namespace_uri_characters_are_checked_without_resolving_relative_names() {
    for uri in [
        "foo bar",
        "urn:a&#xA;b",
        "urn:a%zz",
        "urn:a%",
        "urn:a%0",
        "urn:a&lt;b",
        "urn:a&#x7f;b",
        "urn:a&#x22;b",
        "urn:a&#x5c;b",
    ] {
        let xml = PLAIN.replacen(
            "<fileDescription>",
            &format!("<fileDescription xmlns:p='{uri}'>"),
            1,
        );
        assert!(
            validate_schema_reader(xml.as_bytes(), &SchemaValidationOptions::default()).is_err(),
            "{uri}"
        );
    }
    // Relative namespace references are deprecated rather than forbidden. The
    // safe binding does not expose any associated parser warnings.
    for uri in ["relative", "urn:a%20b"] {
        let xml = PLAIN.replacen(
            "<fileDescription>",
            &format!("<fileDescription xmlns:p='{uri}'>"),
            1,
        );
        assert!(report(&xml).is_valid());
    }
}
#[test]
fn native_writer_documents_are_engine_checks_not_source_generated_goldens() {
    use openms::format::{mzml, peak_options::PeakFileOptions};
    use openms::{MSExperiment, MSSpectrum, Peak1D};
    for e in [
        MSExperiment::default(),
        MSExperiment {
            spectra: vec![MSSpectrum {
                native_id: "scan=1".into(),
                peaks: vec![Peak1D::new(100., 1.)],
                ..Default::default()
            }],
            ..Default::default()
        },
    ] {
        let mut bytes = Vec::new();
        mzml::write(&mut bytes, &e).unwrap();
        let r =
            validate_schema_reader(bytes.as_slice(), &SchemaValidationOptions::default()).unwrap();
        assert!(r.is_valid(), "{r:?}");
        // `mzml::write` is indexed, as the source default is; an experiment
        // with no record to index stays plain mzML.
        assert_eq!(
            r.schema,
            if e.spectra.is_empty() {
                SchemaKind::MzML
            } else {
                SchemaKind::IndexedMzML
            }
        );
        if !e.spectra.is_empty() {
            bytes.clear();
            mzml::write_with_options(&mut bytes, &e, &Default::default()).unwrap();
            let r = validate_schema_reader(bytes.as_slice(), &SchemaValidationOptions::default())
                .unwrap();
            assert!(r.is_valid(), "{r:?}");
            assert_eq!(r.schema, SchemaKind::MzML);
            bytes.clear();
            mzml::write_with_peak_options(&mut bytes, &e, &PeakFileOptions::default()).unwrap();
            let r = validate_schema_reader(bytes.as_slice(), &SchemaValidationOptions::default())
                .unwrap();
            assert!(r.is_valid(), "{r:?}");
            assert_eq!(r.schema, SchemaKind::IndexedMzML);
        }
    }
}
