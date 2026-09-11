#![cfg(feature = "semantic-validation")]
use openms::{
    data_structures::{
        CVMappingRule, CVMappingTerm, CVMappings, CombinationsLogic as Logic,
        RequirementLevel as Level,
    },
    format::{
        controlled_vocabulary::{ControlledVocabulary, OboEncoding},
        cv_mapping::CVMappingFile,
        semantic_validator::{ParsedCVTerm, SemanticValidator},
    },
};
use std::{
    io::{Cursor, Write},
    sync::OnceLock,
};

fn cv(text: &str) -> ControlledVocabulary {
    let mut cv = ControlledVocabulary::new();
    cv.load_obo_reader("test", text.as_bytes()).unwrap();
    cv
}
fn small() -> ControlledVocabulary {
    cv(
        "[Term]\nid: R\nname: root\n[Term]\nid: A\nname: alpha\nis_a: R\n[Term]\nid: B\nname: beta\nis_a: R\n[Term]\nid: C\nname: child\nis_a: A\n",
    )
}
fn mapped(ids: &[&str], level: Level, logic: Logic) -> CVMappings {
    let mut mapping = CVMappings::default();
    mapping.mapping_rules = vec![CVMappingRule {
        identifier: "rule".into(),
        element_path: "/r/cvParam/@accession".into(),
        requirement_level: level,
        combinations_logic: logic,
        terms: ids
            .iter()
            .map(|id| CVMappingTerm {
                accession: (*id).into(),
                use_term: true,
                ..Default::default()
            })
            .collect(),
        ..Default::default()
    }];
    mapping
}
fn read(
    v: &SemanticValidator<'_>,
    xml: &str,
) -> openms::format::semantic_validator::ValidationReport {
    v.validate_reader(xml.as_bytes()).unwrap()
}
fn term(id: &str, name: &str) -> String {
    format!("<cvParam accession=\"{id}\" name=\"{name}\"/>")
}

fn source() -> &'static (CVMappings, ControlledVocabulary) {
    static SOURCE: OnceLock<(CVMappings, ControlledVocabulary)> = OnceLock::new();
    SOURCE.get_or_init(|| {
        let mapping = CVMappingFile::default()
            .read(include_bytes!("data/cv_mapping/SemanticValidator_mapping.xml").as_slice())
            .unwrap();
        let mut cv = ControlledVocabulary::new();
        for (name, bytes, encoding) in [
            (
                "PSI",
                include_bytes!("data/semantic_validator/SemanticValidator_cv.obo").as_slice(),
                OboEncoding::Utf8,
            ),
            (
                "PATO",
                include_bytes!("../resources/cv/quality.obo").as_slice(),
                OboEncoding::Utf8,
            ),
            (
                "UO",
                include_bytes!("../resources/cv/unit.obo").as_slice(),
                OboEncoding::Utf8,
            ),
            (
                "brenda",
                include_bytes!("../resources/cv/brenda.obo").as_slice(),
                OboEncoding::Windows1252,
            ),
            (
                "GO",
                include_bytes!("../resources/cv/goslim_goa.obo").as_slice(),
                OboEncoding::Utf8,
            ),
        ] {
            cv.load_obo_encoded(name, bytes, encoding).unwrap();
        }
        (mapping, cv)
    })
}
#[test]
fn exact_source_valid_and_corrupt_diagnostics_in_order() {
    let (mapping, cv) = source();
    let v = SemanticValidator::new(mapping, cv);
    let report = v
        .validate_reader(
            include_bytes!("data/semantic_validator/SemanticValidator_valid.xml").as_slice(),
        )
        .unwrap();
    assert!(report.is_valid(), "{:?}", report.errors);
    assert!(report.warnings.is_empty());
    let report = v
        .validate_reader(
            include_bytes!("data/semantic_validator/SemanticValidator_corrupt.xml").as_slice(),
        )
        .unwrap();
    assert!(!report.is_valid());
    assert_eq!(
        report.errors,
        [
            "Violated mapping rule 'R3' at element '/mzML/fileDescription/sourceFileList/sourceFile', 2 term(s) should be present, 1 found!",
            "Name of CV term not correct: 'MS:1000554 - LCQ Deca2 - invalid repeat' should be 'LCQ Deca'",
            "CV term used in invalid element: 'MS:1000030 - vendor' at element '/mzML/instrumentConfigurationList/instrumentConfiguration'",
            "Violated mapping rule 'R6a' number of term repeats at element '/mzML/instrumentConfigurationList/instrumentConfiguration'",
            "Violated mapping rule 'R17a' at element '/mzML/run/spectrumList/spectrum/spectrumDescription', 1 term(s) should be present, 0 found!",
        ]
    );
    assert_eq!(
        report.warnings,
        [
            "Unknown CV term: 'MS:1111569 - SHA-1' at element '/mzML/fileDescription/sourceFileList/sourceFile'",
            "Obsolete CV term: 'MS:1000030 - vendor' at element '/mzML/instrumentConfigurationList/instrumentConfiguration'",
            "No mapping rule found for element '/mzML/acquisitionSettingsList/acquisitionSettings/targetList/target'",
            "No mapping rule found for element '/mzML/acquisitionSettingsList/acquisitionSettings/targetList/target'",
        ]
    );
}
#[test]
fn requirement_truth_table_including_source_should_omission() {
    let cv = small();
    for level in [Level::Must, Level::Should, Level::May] {
        for logic in [Logic::And, Logic::Or, Logic::Xor] {
            for count in 0..=2 {
                let mapping = mapped(&["A", "B"], level, logic);
                let body = [term("A", "alpha"), term("B", "beta")][..count].join("");
                let report = read(
                    &SemanticValidator::new(&mapping, &cv),
                    &format!("<r>{body}</r>"),
                );
                let expected = match (level, logic) {
                    (Level::Must, Logic::And) => count == 2,
                    (Level::Must, Logic::Or) => count >= 1,
                    (Level::Must, Logic::Xor) => count == 1,
                    (Level::May, Logic::And) => count == 0 || count == 2,
                    (Level::May, Logic::Xor) => count <= 1,
                    _ => true,
                };
                assert_eq!(
                    report.is_valid(),
                    expected,
                    "{level:?} {logic:?} {count}: {report:?}"
                );
            }
        }
    }
}
#[test]
fn empty_rules_absent_elements_repeat_checks_and_sibling_reset() {
    let cv = small();
    for (logic, valid) in [(Logic::And, true), (Logic::Or, false), (Logic::Xor, false)] {
        let mapping = mapped(&[], Level::Must, logic);
        let v = SemanticValidator::new(&mapping, &cv);
        assert_eq!(read(&v, "<r/>").is_valid(), valid);
        assert!(read(&v, "<other/>").is_valid());
    }
    let mut mapping = mapped(&["A"], Level::May, Logic::Or);
    mapping.mapping_rules[0].element_path = "/outer/r/cvParam/@accession".into();
    let v = SemanticValidator::new(&mapping, &cv);
    assert!(
        read(
            &v,
            &format!(
                "<outer><r>{}</r><r>{}</r></outer>",
                term("A", "alpha"),
                term("A", "alpha")
            )
        )
        .is_valid()
    );
    let r = read(
        &v,
        &format!(
            "<outer><r>{}{}</r></outer>",
            term("A", "alpha"),
            term("A", "alpha")
        ),
    );
    assert_eq!(
        r.errors,
        ["Violated mapping rule 'rule' number of term repeats at element '/outer/r'"]
    );
}
#[test]
fn first_match_each_rule_duplicates_and_ignored_scope_name_fields() {
    let cv = small();
    let mut mapping = mapped(&["R", "A"], Level::Must, Logic::And);
    mapping.mapping_rules[0].terms[0].use_term = false;
    mapping.mapping_rules[0].terms[0].allow_children = true;
    mapping.mapping_rules[0].scope_path = "/never".into();
    mapping.mapping_rules[0].terms[0].use_term_name = true;
    mapping.mapping_rules[0].terms[0].term_name = "wrong".into();
    let r = read(
        &SemanticValidator::new(&mapping, &cv),
        &format!("<r>{}</r>", term("A", "alpha")),
    );
    assert_eq!(
        r.errors,
        ["Violated mapping rule 'rule' at element '/r', 2 term(s) should be present, 1 found!"]
    );
    // Duplicate mapped accessions count twice in source distinct-entry matching,
    // but one incoming term increments the shared accession counter only once.
    let mapping = mapped(&["A", "A"], Level::Must, Logic::And);
    assert!(
        read(
            &SemanticValidator::new(&mapping, &cv),
            &format!("<r>{}</r>", term("A", "alpha"))
        )
        .is_valid()
    );
    // Repeated rule IDs share counters; both matches increment the same key.
    let mut mapping = mapped(&["A"], Level::Must, Logic::Or);
    mapping.mapping_rules.push(mapping.mapping_rules[0].clone());
    let r = read(
        &SemanticValidator::new(&mapping, &cv),
        &format!("<r>{}</r>", term("A", "alpha")),
    );
    assert_eq!(r.errors.len(), 2);
    assert!(
        r.errors
            .iter()
            .all(|s| s.contains("number of term repeats"))
    );
}
#[test]
fn configurable_qnames_attributes_flags_and_root_double_slash() {
    let cv = small();
    let mut mapping = mapped(&["A"], Level::Must, Logic::Or);
    mapping.mapping_rules[0].element_path = "//param/@id".into();
    let mut v = SemanticValidator::new(&mapping, &cv);
    v.options.tag = "param".into();
    v.options.accession_attribute = "id".into();
    v.options.name_attribute = "label".into();
    v.options.value_attribute = "text".into();
    // Root term lookup uses '//param/@id'; its closing path is '/param/param/@id'.
    let r = read(&v, "<param id=\"A\" label=\"alpha\"/>");
    assert!(r.errors.is_empty());
    assert!(r.warnings.is_empty());
    assert!(
        read(&v, "<param id=\"A\" label=\"alpha\" text=\"x\"/>").errors[0]
            .contains("Value of CV term not allowed")
    );
    v.options.check_term_value_types = false;
    assert!(read(&v, "<param id=\"A\" label=\"alpha\" text=\"x\"/>").is_valid());
    assert!(v.validate_reader(b"<param id=\"A\"/>".as_slice()).is_err());
    let mut mapping = mapped(&["A"], Level::Must, Logic::Or);
    mapping.mapping_rules[0].element_path = "/n:r/n:cv/@n:id".into();
    let mut v = SemanticValidator::new(&mapping, &cv);
    v.options.tag = "n:cv".into();
    v.options.accession_attribute = "n:id".into();
    assert!(
        read(
            &v,
            "<n:r xmlns:n=\"urn:example\"><n:cv n:id=\"A\" name=\"alpha\"/></n:r>"
        )
        .is_valid()
    );
}
#[test]
fn locate_is_selection_only_and_absent_path_stays_error_cpp044() {
    let cv = small();
    let mapping = mapped(&["missing"], Level::May, Logic::Or);
    let v = SemanticValidator::new(&mapping, &cv);
    let p = ParsedCVTerm {
        accession: "missing".into(),
        name: "wrong".into(),
        value: "wrong".into(),
        has_value: true,
        ..Default::default()
    };
    assert!(v.locate_term("/r/cvParam/@accession", &p).unwrap());
    assert!(v.locate_term("/x/cvParam/@accession", &p).is_err());
    assert!(read(&v, "<x/>").is_valid());
    assert!(v.locate_term("/x/cvParam/@accession", &p).is_err());
    let mut mapping = mapped(&["R"], Level::Must, Logic::Or);
    mapping.mapping_rules[0].terms[0].use_term = false;
    mapping.mapping_rules[0].terms[0].allow_children = true;
    let v = SemanticValidator::new(&mapping, &cv);
    assert!(
        v.locate_term(
            "/r/cvParam/@accession",
            &ParsedCVTerm {
                accession: "C".into(),
                ..Default::default()
            }
        )
        .unwrap()
    );
    assert!(
        !v.locate_term(
            "/r/cvParam/@accession",
            &ParsedCVTerm {
                accession: "R".into(),
                ..Default::default()
            }
        )
        .unwrap()
    );
}
#[test]
fn failed_xml_never_publishes_or_contaminates_next_call_cpp040() {
    let cv = small();
    let mapping = mapped(&["A"], Level::Must, Logic::Or);
    let v = SemanticValidator::new(&mapping, &cv);
    for bad in [
        "<r><cvParam accession=\"A\"/></r>",
        "<r><cvParam accession=\"A\" name=\"alpha\"/><broken>",
    ] {
        assert!(v.validate_reader(bad.as_bytes()).is_err());
        assert_eq!(read(&v, "<r/>").errors.len(), 1);
        assert!(read(&v, &format!("<r>{}</r>", term("A", "alpha"))).is_valid());
    }
}
#[test]
fn unit_exact_child_unrelated_and_flags_cpp039() {
    let cv = cv(
        "[Term]\nid: M\nname: measurement\nrelationship: has_units U:0\n[Term]\nid: U:0\nname: unit\n[Term]\nid: V:0\nname: child unit\nis_a: U:0\n[Term]\nid: W:0\nname: unrelated\n[Term]\nid: N\nname: unitless\n",
    );
    let mapping = mapped(&["M", "N"], Level::May, Logic::Or);
    let mut v = SemanticValidator::new(&mapping, &cv);
    v.options.check_units = true;
    v.options.unit_accession_attribute = "unit".into();
    v.options.unit_name_attribute = "uname".into();
    for unit in ["U:0", "V:0"] {
        assert!(read(&v,&format!("<r><cvParam accession=\"M\" name=\"measurement\" unit=\"{unit}\" uname=\"label ignored\"/></r>")).is_valid());
    }
    for (attrs, message) in [
        ("", "CV term must have a unit: M - measurement"),
        (
            "uname=\"name\"",
            "CV term must have a unit: M - measurement",
        ),
        (
            "unit=\"W:0\" uname=\"other\"",
            "Unit CV term not allowed: W:0 - other of term M - measurement",
        ),
        (
            "unit=\"unknown\"",
            "Unit CV term not found: unknown -  of term M - measurement",
        ),
    ] {
        assert_eq!(
            read(
                &v,
                &format!("<r><cvParam accession=\"M\" name=\"measurement\" {attrs}/></r>")
            )
            .errors,
            [message]
        );
    }
    let r = read(
        &v,
        "<r><cvParam accession=\"N\" name=\"unitless\" uname=\"\"/></r>",
    );
    assert_eq!(
        r.warnings,
        ["Unit CV term used, but not allowed:  -  of term N - unitless"]
    );
    v.options.check_units = false;
    assert!(
        read(
            &v,
            "<r><cvParam accession=\"M\" name=\"measurement\" unit=\"unknown\"/></r>"
        )
        .is_valid()
    );
}
#[test]
fn unknown_obsolete_name_and_diagnostic_order() {
    let cv = cv(
        "[Term]\nid: A\nname: alpha\nis_obsolete: true\nxref: value-type:xsd\\:integer\nrelationship: has_units U:0\n[Term]\nid: U:0\nname: unit\n",
    );
    let mapping = mapped(&["notA"], Level::May, Logic::Or);
    let mut v = SemanticValidator::new(&mapping, &cv);
    v.options.check_units = true;
    let r = read(
        &v,
        "<r><cvParam accession=\"A\" name=\" wrong \" value=\"bad\"/></r>",
    );
    assert_eq!(
        r.warnings,
        ["Obsolete CV term: 'A -  wrong ' at element '/r'"]
    );
    assert_eq!(
        r.errors,
        [
            "CV term must have a unit: A -  wrong ",
            "CV term used in invalid element: 'A -  wrong ' at element '/r'",
            "Name of CV term not correct: 'A - wrong' should be 'alpha'",
            "Value-type of CV term wrong, should be xsd:integer: 'A -  wrong ' value='bad' at element '/r'"
        ]
    );
    let r = read(
        &v,
        "<r><cvParam accession=\"unknown\" name=\"bogus\" value=\"bad\"/></r>",
    );
    assert!(r.errors.is_empty());
    assert_eq!(r.warnings.len(), 1);
    let mapping = mapped(&["A"], Level::May, Logic::Or);
    let mut v = SemanticValidator::new(&mapping, &cv);
    v.options.check_term_value_types = false;
    assert!(
        read(
            &v,
            "<r><cvParam accession=\"A\" name=\"&#9; alpha &#10;\"/></r>"
        )
        .is_valid()
    );
    assert!(!read(&v, "<r><cvParam accession=\"A\" name=\"Alpha\"/></r>").is_valid());
}
#[test]
fn all_value_types_preserve_source_conversions_and_empty_value_rules() {
    let cases: [(&str, &[&str], &[&str]); 10] = [
        ("string", &["", "anything"], &[]),
        (
            "integer",
            &["2147483647", "-2147483648", "  +-1  "],
            &["2147483648", "1.0", "x"],
        ),
        (
            "decimal",
            &["1.2", "NaN", "INF", "-0", "1e-300"],
            &["1e9999", "1e-9999", "x"],
        ),
        ("negativeInteger", &["-1"], &["0", "1"]),
        ("positiveInteger", &["1"], &["0", "-1"]),
        ("nonNegativeInteger", &["0", "1"], &["-1"]),
        ("nonPositiveInteger", &["0", "-1"], &["1"]),
        ("boolean", &["1", "0", "  TrUe ", "FALSE"], &["yes", "2"]),
        // CPP-046: retain DateTime::set behavior, not strict XSD date grammar.
        (
            "date",
            &[
                "2001-02-03 04:05:06",
                "2001-02-03T04:05:06.789",
                "03.02.2001 04:05:06",
            ],
            &["2001-02-29", "2001-02-03", "03.02.2001", "bad"],
        ),
        (
            "anyURI",
            &[":", "not a uri:still accepted"],
            &["example.com"],
        ),
    ];
    for (ty, good, bad) in cases {
        let cv = cv(&format!(
            "[Term]\nid: A\nname: alpha\nxref: value-type:xsd\\:{ty}\n"
        ));
        let mapping = mapped(&["A"], Level::May, Logic::Or);
        let v = SemanticValidator::new(&mapping, &cv);
        assert!(
            !read(&v, "<r><cvParam accession=\"A\" name=\"alpha\"/></r>").is_valid(),
            "missing {ty}"
        );
        for (values, expected) in [(good, true), (bad, false)] {
            for value in values {
                let r = read(
                    &v,
                    &format!("<r><cvParam accession=\"A\" name=\"alpha\" value=\"{value}\"/></r>"),
                );
                assert_eq!(r.is_valid(), expected, "{ty} {value}: {r:?}");
            }
        }
    }
    let cv = cv("[Term]\nid: A\nname: alpha\n[Term]\nid: PATO:X\nname: quality\n");
    let mapping = mapped(&["A", "PATO:X"], Level::May, Logic::Or);
    let v = SemanticValidator::new(&mapping, &cv);
    for value in [None, Some("")] {
        let attr = value.map(|s| format!(" value=\"{s}\"")).unwrap_or_default();
        assert!(
            read(
                &v,
                &format!("<r><cvParam accession=\"A\" name=\"alpha\"{attr}/></r>")
            )
            .is_valid()
        );
    }
    assert!(
        read(
            &v,
            "<r><cvParam accession=\"PATO:X\" name=\"quality\" value=\"x\"/></r>"
        )
        .is_valid()
    );
}
#[test]
fn compressed_paths_and_xml_lexical_failures() {
    let cv = small();
    let mapping = mapped(&["A"], Level::May, Logic::Or);
    let v = SemanticValidator::new(&mapping, &cv);
    let xml = format!("<r>{}</r>", term("A", "alpha"));
    for zipped in [false, true] {
        let bytes = if zipped {
            let mut w = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
            w.write_all(xml.as_bytes()).unwrap();
            w.finish().unwrap()
        } else {
            let mut w = bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::default());
            w.write_all(xml.as_bytes()).unwrap();
            w.finish().unwrap()
        };
        let path = std::env::temp_dir().join(format!(
            "openms-semantic-{}-{zipped}.raw",
            std::process::id()
        ));
        std::fs::write(&path, &bytes).unwrap();
        assert!(v.validate(&path).unwrap().is_valid());
        std::fs::write(&path, &bytes[..bytes.len() / 2]).unwrap();
        assert!(v.validate(&path).is_err());
        std::fs::remove_file(path).unwrap();
    }
    for input in [
        "",
        "<",
        "<r a=\"1\"a=\"2\"/>",
        "<r a=\"1\" a=\"2\"/>",
        "<r></s>",
        "<r/>x",
        "<r>]]></r>",
        "<r>&bad;</r>",
        "<!DOCTYPE r><r/>",
        "<?xml version=\"1.0\" standalone=\"yes\" encoding=\"UTF-8\"?><r/>",
        "<r><!-- a--b --></r>",
        "<r><?XmL bad?></r>",
        "<r>&#0;</r>",
    ] {
        assert!(v.validate_reader(input.as_bytes()).is_err(), "{input}");
    }
    assert!(v.validate("/does/not/exist/openms-semantic").is_err());
}
#[test]
fn utf16_declarations_and_ignored_text_still_checked() {
    let cv = small();
    let mapping = mapped(&["A"], Level::May, Logic::Or);
    let v = SemanticValidator::new(&mapping, &cv);
    for little in [false, true] {
        let text = "<?xml version=\"1.0\" encoding=\"UTF-16\"?><r>é &amp; ignored<![CDATA[<]]></r>";
        let mut bytes = if little {
            vec![255, 254]
        } else {
            vec![254, 255]
        };
        for c in text.encode_utf16() {
            bytes.extend(if little {
                c.to_le_bytes()
            } else {
                c.to_be_bytes()
            });
        }
        assert!(v.validate_reader(Cursor::new(bytes)).unwrap().is_valid());
    }
    assert!(read(&v, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?><r/>").is_valid());
    assert!(
        v.validate_reader(b"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?><r>\xe9</r>".as_slice())
            .is_err()
    );
    assert!(read(&v, "<r>&amp; text &#10;<child/></r>").is_valid());
}
#[test]
fn operation_limits_and_error_precedence_are_not_hidden_reports() {
    let cv = small();
    let mapping = mapped(&["A"], Level::Must, Logic::Or);
    for option in 0..9 {
        let mut v = SemanticValidator::new(&mapping, &cv);
        match option {
            0 => v.options.limits.max_input_bytes = 1,
            1 => v.options.limits.max_depth = 1,
            2 => v.options.limits.max_elements = 1,
            3 => v.options.limits.max_terms = 0,
            4 => v.options.limits.max_rules = 0,
            5 => v.options.limits.max_mapping_terms = 0,
            6 => v.options.limits.max_diagnostics = 0,
            7 => v.options.limits.max_work = 0,
            _ => v.options.limits.max_bytes = 0,
        }
        let xml = if option == 6 {
            "<r/>".into()
        } else {
            format!("<r>{}</r>", term("A", "alpha"))
        };
        let err = v.validate_reader(xml.as_bytes()).unwrap_err();
        assert!(err.to_string().contains("limit"), "{option}: {err}");
    }
    let mut v = SemanticValidator::new(&mapping, &cv);
    v.options.limits.max_rules = 0;
    // Native mapping cap fires before touching malformed XML or enormous key payload.
    assert!(
        v.validate_reader(b"<broken".as_slice())
            .unwrap_err()
            .to_string()
            .contains("limit")
    );
    let mut mapping = mapped(&["R"], Level::May, Logic::Or);
    mapping.mapping_rules[0].terms[0].allow_children = true;
    let mut v = SemanticValidator::new(&mapping, &cv);
    v.options.limits.max_work = 1;
    assert!(
        v.locate_term("/r/cvParam/@accession", &Default::default())
            .is_err()
    );
}

#[test]
fn descendant_matching_retains_source_stale_children_index() {
    let mut cv = small();
    // An additive load changes A's parent record but source R.children retains A.
    cv.load_obo_reader("later", b"[Term]\nid: A\nname: alpha\n".as_slice())
        .unwrap();
    assert!(!cv.is_child_of("A", "R").unwrap());
    let mut mapping = mapped(&["R"], Level::Must, Logic::Or);
    mapping.mapping_rules[0].terms[0].use_term = false;
    mapping.mapping_rules[0].terms[0].allow_children = true;
    let v = SemanticValidator::new(&mapping, &cv);
    assert!(
        v.locate_term(
            "/r/cvParam/@accession",
            &ParsedCVTerm {
                accession: "A".into(),
                ..Default::default()
            }
        )
        .unwrap()
    );
    assert!(read(&v, &format!("<r>{}</r>", term("A", "alpha"))).is_valid());
}

#[test]
fn raw_historical_fixture_bytes_are_stable_on_every_platform() {
    // FNV constants independently calculated from the pinned original bytes;
    // cryptographic SHA-256 and full source paths are in the provenance manifest.
    for (raw, len, hash) in [
        (
            include_bytes!("data/semantic_validator/SemanticValidator_valid.xml").as_slice(),
            21806,
            12788848976010737090u64,
        ),
        (
            include_bytes!("data/semantic_validator/SemanticValidator_corrupt.xml").as_slice(),
            22213,
            5330208474211620510,
        ),
        (
            include_bytes!("data/semantic_validator/SemanticValidator_cv.obo").as_slice(),
            169777,
            17457249215967672992,
        ),
    ] {
        assert_eq!(raw.len(), len);
        assert_eq!(
            raw.iter()
                .fold(14695981039346656037u64, |h, b| (h ^ u64::from(*b))
                    .wrapping_mul(1099511628211)),
            hash
        );
    }
}
