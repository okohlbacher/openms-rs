// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "cv-mapping")]
use openms::{
    data_structures::*,
    format::cv_mapping::{CVMappingFile, ReadOptions},
};
use std::{
    io::{BufReader, Cursor, Write},
    path::PathBuf,
};

const RULE: &str = r#"<CvMappingRule id="r" cvElementPath="/n:root/cvParam/@n:accession" requirementLevel="MUST" scopePath="/n:root" cvTermsCombinationLogic="OR">"#;
const TERM: &str = r#"<CvTerm termAccession="MS:1" useTerm="true" termName="name" allowChildren="false" cvIdentifierRef="MS"/>"#;
fn doc(content: &str) -> String {
    format!("<CvMapping>{content}</CvMapping>")
}
fn parse(s: &str) -> openms::Result<CVMappings> {
    CVMappingFile::default().read(s.as_bytes())
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/cv_mapping")
        .join(name)
}
fn hex(s: &str) -> String {
    const DIGITS: &[u8] = b"0123456789abcdef";
    let mut output = String::with_capacity(s.len() * 2);
    for &byte in s.as_bytes() {
        output.push(DIGITS[(byte >> 4) as usize] as char);
        output.push(DIGITS[(byte & 15) as usize] as char);
    }
    output
}
fn level(v: RequirementLevel) -> &'static str {
    match v {
        RequirementLevel::Must => "MUST",
        RequirementLevel::Should => "SHOULD",
        RequirementLevel::May => "MAY",
    }
}
fn logic(v: CombinationsLogic) -> &'static str {
    match v {
        CombinationsLogic::Or => "OR",
        CombinationsLogic::And => "AND",
        CombinationsLogic::Xor => "XOR",
    }
}

#[test]
fn original_class_test_literal_records() {
    let m = CVMappingFile::default()
        .load(fixture("cv_mapping_test_file.xml"))
        .unwrap();
    assert_eq!(m.mapping_rules.len(), 9);
    assert_eq!(
        m.mapping_rules
            .iter()
            .map(|r| r.identifier.as_str())
            .collect::<Vec<_>>(),
        ["0", "1", "2", "3", "4", "5", "6", "7", "8"]
    );
    assert_eq!(
        [
            m.mapping_rules[0].terms.len(),
            m.mapping_rules[1].terms.len(),
            m.mapping_rules[2].terms.len()
        ],
        [14, 32, 46]
    );
    let r = &m.mapping_rules[0];
    assert_eq!(
        r.element_path,
        "/mzData/description/admin/sampleDescription/cvParam/@accession"
    );
    assert_eq!(r.scope_path, "/mzData/description/admin/sampleDescription");
    assert_eq!(r.requirement_level, RequirementLevel::May);
    assert_eq!(r.combinations_logic, CombinationsLogic::Or);
    assert_eq!(
        m.mapping_rules[1].requirement_level,
        RequirementLevel::Should
    );
    assert_eq!(
        m.mapping_rules[1].combinations_logic,
        CombinationsLogic::Xor
    );
    assert_eq!(m.mapping_rules[2].requirement_level, RequirementLevel::Must);
    assert_eq!(
        m.mapping_rules[2].combinations_logic,
        CombinationsLogic::And
    );
    let names = [
        "Sample Number",
        "Sample Name",
        "Sample State",
        "Sample Mass",
        "Sample Volume",
        "Sample Concentration",
    ];
    let expected = [
        (false, true, true, true),
        (true, true, true, true),
        (true, false, true, true),
        (true, false, false, true),
        (true, false, false, false),
        (false, true, true, true),
    ];
    for i in 0..6 {
        let t = &r.terms[i];
        assert_eq!(t.accession, format!("PSI:100000{}", i + 1));
        assert_eq!(t.term_name, names[i]);
        assert_eq!(
            (
                t.use_term_name,
                t.use_term,
                t.is_repeatable,
                t.allow_children
            ),
            expected[i]
        );
        assert_eq!(t.cv_identifier_ref, "PSI");
    }
    assert_eq!(
        m.cv_references(),
        [CVReference {
            name: "mzData CV".into(),
            identifier: "PSI".into()
        }]
    );
    assert!(m.has_cv_reference("PSI"));
}

#[test]
fn all_six_source_files_match_complete_independent_xml_projection() {
    let names = [
        "SemanticValidator_mapping.xml",
        "TraML-mapping.xml",
        "cv_mapping_test_file.xml",
        "ms-mapping.xml",
        "mzIdentML-mapping.xml",
        "mzdata-mapping.xml",
    ];
    let mut actual = String::new();
    for n in names {
        let m = CVMappingFile::default().load(fixture(n)).unwrap();
        for r in m.cv_references() {
            actual += &format!("{n}\treference\t{}\t{}\n", hex(&r.name), hex(&r.identifier));
        }
        for (ri, r) in m.mapping_rules.iter().enumerate() {
            actual += &format!(
                "{n}\trule\t{ri}\t{}\t{}\t{}\t{}\t{}\n",
                hex(&r.identifier),
                hex(&r.element_path),
                level(r.requirement_level),
                hex(&r.scope_path),
                logic(r.combinations_logic)
            );
            for (ti, t) in r.terms.iter().enumerate() {
                actual += &format!(
                    "{n}\tterm\t{ri}\t{ti}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    hex(&t.accession),
                    t.use_term_name,
                    t.use_term,
                    hex(&t.term_name),
                    t.is_repeatable,
                    t.allow_children,
                    hex(&t.cv_identifier_ref)
                );
            }
        }
    }
    assert_eq!(actual, include_str!("data/cv_mapping/projection.tsv"));
    assert_eq!(actual.lines().count(), 683);
}

#[test]
fn xml_defaults_exact_booleans_and_source_enum_fallbacks() {
    let xml = doc(&format!("{RULE}{TERM}</CvMappingRule>"));
    let m = parse(&xml).unwrap();
    let t = &m.mapping_rules[0].terms[0];
    assert!(!t.use_term_name && t.is_repeatable);
    let m =
        parse(&xml.replace("termName=", r#"useTermName="" isRepeatable="" termName="#)).unwrap();
    assert!(
        !m.mapping_rules[0].terms[0].use_term_name && m.mapping_rules[0].terms[0].is_repeatable
    );
    for attribute in ["useTerm", "allowChildren"] {
        for v in ["TRUE", "1", "0", " true ", ""] {
            assert!(
                parse(&xml.replace(
                    &format!(
                        "{attribute}=\"{}\"",
                        if attribute == "useTerm" {
                            "true"
                        } else {
                            "false"
                        }
                    ),
                    &format!("{attribute}=\"{v}\"")
                ))
                .is_err()
            );
        }
    }
    // CPP-035 compatibility: present unknown enums retain source defaults.
    let m = parse(&xml.replace("MUST", "MAYY").replace("OR", "ANDD")).unwrap();
    assert_eq!(m.mapping_rules[0].requirement_level, RequirementLevel::Must);
    assert_eq!(m.mapping_rules[0].combinations_logic, CombinationsLogic::Or);
    for attr in [
        "termAccession",
        "useTerm",
        "termName",
        "allowChildren",
        "cvIdentifierRef",
    ] {
        let key = format!("{attr}=\"");
        let start = xml.find(&key).unwrap();
        let end = xml[start + key.len()..].find('"').unwrap() + start + key.len() + 1;
        let mut broken = xml.clone();
        broken.replace_range(start..end, "");
        assert!(parse(&broken).is_err(), "{attr}");
    }
}

#[test]
fn corrected_namespace_stripping_preserves_attributes_and_scope() {
    let loader = CVMappingFile {
        options: ReadOptions {
            strip_namespaces: true,
            ..Default::default()
        },
    };
    for (input, expected) in [
        ("/n:root/cvParam/@n:accession", "/root/cvParam/@accession"),
        ("plain", "/plain"),
        ("//n:a///b/", "/a/b"),
        ("", ""),
        ("/", ""),
        ("/n:/@n:", "//@"),
    ] {
        let xml = doc(&format!(
            "{} </CvMappingRule>",
            RULE.replace("/n:root/cvParam/@n:accession", input)
        ));
        let m = loader.read(xml.as_bytes()).unwrap();
        assert_eq!(m.mapping_rules[0].element_path, expected);
        assert_eq!(m.mapping_rules[0].scope_path, "/n:root");
    }
    assert!(
        loader
            .read(
                doc(&format!(
                    "{} </CvMappingRule>",
                    RULE.replace("/n:root/cvParam/@n:accession", "/a:b:c")
                ))
                .as_bytes()
            )
            .is_err()
    );
    // CPP-032: source's own ordinary unprefixed fixture now supports stripping.
    let m = loader.load(fixture("cv_mapping_test_file.xml")).unwrap();
    assert_eq!(
        m.mapping_rules[0].element_path,
        "/mzData/description/admin/sampleDescription/cvParam/@accession"
    );
}

#[test]
fn source_qnames_ignored_payload_and_single_current_rule_state() {
    let xml = doc(&format!(
        "text&amp;more<!--x--><Unknown a='ok'>{RULE}{TERM}<CvMappingRule id='inner' cvElementPath='x' requirementLevel='MAY' scopePath='' cvTermsCombinationLogic='AND'/></CvMappingRule></Unknown><p:CvReference cvName='ignored' cvIdentifier='ignored'/>"
    ));
    let m = parse(&xml).unwrap();
    assert_eq!(m.mapping_rules.len(), 2);
    assert_eq!(m.mapping_rules[0].identifier, "inner");
    assert_eq!(m.mapping_rules[0].terms.len(), 1);
    assert_eq!(m.mapping_rules[1], CVMappingRule::default());
    assert!(m.cv_references().is_empty());
    assert!(
        parse("<differentRoot>ignored<![CDATA[<&]]></differentRoot>")
            .unwrap()
            .mapping_rules
            .is_empty()
    );
}

#[test]
fn atomic_repeated_load_references_append_rules_replace_and_failure_does_not_leak() {
    let loader = CVMappingFile::default();
    let one = doc(&format!(
        "<CvReference cvName='old' cvIdentifier='OLD'/>{RULE}{TERM}</CvMappingRule>"
    ));
    let mut m = loader.read(one.as_bytes()).unwrap();
    loader.read_into(one.as_bytes(), &mut m).unwrap();
    assert_eq!(m.cv_references().len(), 2);
    assert_eq!(m.mapping_rules.len(), 1);
    let saved = m.clone();
    let bad = doc("<CvReference cvName='leak' cvIdentifier='LEAK'/><CvMappingRule/>");
    assert!(loader.read_into(bad.as_bytes(), &mut m).is_err());
    assert_eq!(m, saved);
    loader
        .read_into(b"<CvMapping/>".as_slice(), &mut m)
        .unwrap();
    assert_eq!(m.cv_references().len(), 2);
    assert!(m.mapping_rules.is_empty());
    assert!(!m.has_cv_reference("LEAK"));
    // CPP-033: fresh result from reused reader cannot inherit failed draft records.
    assert!(loader.read(bad.as_bytes()).is_err());
    let fresh = loader.read(b"<CvMapping/>".as_slice()).unwrap();
    assert!(fresh.cv_references().is_empty());
}

#[test]
fn xml_entities_line_endings_and_chunk_independent_encoding() {
    let xml = "<?xml version='1.0' encoding='UTF-8'?><CvMapping><CvReference cvName='A\r\nB\tC&#10;&amp;&quot;&apos;&#x3b1;' cvIdentifier='id'/></CvMapping>";
    for capacity in 1..8 {
        let m = CVMappingFile::default()
            .read(BufReader::with_capacity(capacity, Cursor::new(xml)))
            .unwrap();
        assert_eq!(m.cv_references()[0].name, "A B C\n&\"'α");
    }
    let bom = [&[239, 187, 191][..], xml.as_bytes()].concat();
    assert!(CVMappingFile::default().read(bom.as_slice()).is_ok());
    for little in [true, false] {
        let text = xml.replace("UTF-8", "UTF-16");
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
        assert!(CVMappingFile::default().read(bytes.as_slice()).is_ok());
        bytes.push(0);
        assert!(CVMappingFile::default().read(bytes.as_slice()).is_err());
    }
    assert!(parse("<?xml version='1.0' encoding='ISO-8859-1'?><r/>").is_ok());
    assert!(parse("<?xml version='1.0' encoding='ISO-8859-1'?><r>α</r>").is_err());
}

#[test]
fn malformed_xml_and_short_inputs_fail_without_publication() {
    for xml in [
        "",
        "B",
        "P",
        "<",
        "<r>",
        "<r></s>",
        "<r/><s/>",
        "x<r/>",
        "<r/>x",
        "<r>]]></r>",
        "<r a='x'b='y'/>",
        "<r a='x' a='y'/>",
        "<1bad/>",
        "<r 1bad='x'/>",
        "<r a='<'/>",
        "<r>&unknown;</r>",
        "<r>&#0;</r>",
        "<r a='&#xFFFF;'/>",
        "<!--a--b--><r/>",
        "<?xml version='1.0' standalone='yes' encoding='UTF-8'?><r/>",
        " <?xml version='1.0'?><r/>",
        "<?XML x?><r/>",
        "<!DOCTYPE r><r/>",
        "<r/>\0",
        "<r/>&#32;",
    ] {
        assert!(parse(xml).is_err(), "accepted {xml:?}");
    }
}

#[test]
fn all_record_and_payload_limits_are_atomic() {
    let xml = doc(&format!(
        "<CvReference cvName='n' cvIdentifier='id'/>{RULE}{TERM}</CvMappingRule>"
    ));
    let options = [
        ReadOptions {
            max_input_bytes: xml.len() - 1,
            ..Default::default()
        },
        ReadOptions {
            max_elements: 3,
            ..Default::default()
        },
        ReadOptions {
            max_depth: 2,
            ..Default::default()
        },
        ReadOptions {
            max_references: 0,
            ..Default::default()
        },
        ReadOptions {
            max_rules: 0,
            ..Default::default()
        },
        ReadOptions {
            max_terms: 0,
            ..Default::default()
        },
        ReadOptions {
            max_work: 1,
            ..Default::default()
        },
        ReadOptions {
            max_bytes: 1,
            ..Default::default()
        },
    ];
    for options in options {
        let loader = CVMappingFile { options };
        let mut m = CVMappings::default();
        m.add_cv_reference(CVReference {
            name: "kept".into(),
            identifier: "OLD".into(),
        });
        let before = m.clone();
        assert!(loader.read_into(xml.as_bytes(), &mut m).is_err());
        assert_eq!(m, before);
    }
    let exact = CVMappingFile {
        options: ReadOptions {
            max_input_bytes: xml.len(),
            ..Default::default()
        },
    };
    assert!(exact.read(xml.as_bytes()).is_ok());
    let mut m = CVMappings::default();
    m.add_cv_reference(CVReference {
        name: "x".repeat(100_000),
        identifier: "id".into(),
    });
    let before = m.clone();
    let low = CVMappingFile {
        options: ReadOptions {
            max_bytes: 10_000,
            ..Default::default()
        },
    };
    assert!(low.read_into(b"<r/>".as_slice(), &mut m).is_err());
    assert_eq!(m, before);
}

#[test]
fn compressed_magic_corruption_and_io_failure() {
    static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let root = std::env::temp_dir().join(format!(
        "openms-cv-map-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
    ));
    std::fs::create_dir(&root).unwrap();
    let xml = doc("<CvReference cvName='n' cvIdentifier='id'/>");
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gz.write_all(xml.as_bytes()).unwrap();
    let gzip = gz.finish().unwrap();
    let mut bz = bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::default());
    bz.write_all(xml.as_bytes()).unwrap();
    let bzip = bz.finish().unwrap();
    for (i, bytes) in [gzip.clone(), bzip].into_iter().enumerate() {
        let p = root.join(format!("wrong{i}.txt"));
        std::fs::write(&p, bytes).unwrap();
        assert!(
            CVMappingFile::default()
                .load(&p)
                .unwrap()
                .has_cv_reference("id")
        );
    }
    let mut broken = gzip;
    let n = broken.len();
    broken[n - 8] ^= 1;
    let p = root.join("broken.xml");
    std::fs::write(&p, broken).unwrap();
    assert!(CVMappingFile::default().load(&p).is_err());
    assert!(CVMappingFile::default().load(root.join("missing")).is_err());
    std::fs::remove_dir_all(root).unwrap();
    struct Failed;
    impl std::io::Read for Failed {
        fn read(&mut self, _: &mut [u8]) -> std::io::Result<usize> {
            Err(std::io::ErrorKind::Other.into())
        }
    }
    impl std::io::BufRead for Failed {
        fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
            Err(std::io::ErrorKind::Other.into())
        }
        fn consume(&mut self, _: usize) {}
    }
    assert!(CVMappingFile::default().read(Failed).is_err());
}
