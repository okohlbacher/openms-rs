// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "idxml")]

use openms::comparison::Tolerance;
use openms::format::idxml::{self, IdXmlDocument, ReadOptions, WriteOptions};
use openms::identification::{
    AnalysisResult, EnzymeTermSpecificity, FlankingResidue, PeakAnnotation, PeakMassType,
    PeptideIdentification, ProteinIdentification,
};
use openms::metadata::{MetaValue, MetaValueData, Unit};
use std::io::{Cursor, Write};

const WHOLE: &str = include_str!("data/idxml_upstream_whole.idXML");
const NO_PROTEINS: &str = include_str!("data/idxml_upstream_no_proteinhits.idXML");
fn parse(value: &str) -> openms::Result<IdXmlDocument> {
    idxml::read(Cursor::new(value.as_bytes()))
}
fn encoded(value: &IdXmlDocument) -> String {
    let mut result = Vec::new();
    idxml::write(&mut result, value).unwrap();
    String::from_utf8(result).unwrap()
}
fn sample() -> IdXmlDocument {
    let mut doc = parse(WHOLE).unwrap();
    doc.document_id = "µ & <\"document\">\nline".into();
    let run = &mut doc.protein_identifications[0];
    run.search_parameters.enzyme_specificity = EnzymeTermSpecificity::Semi;
    run.search_parameters.fragment_tolerance = Tolerance::Ppm(12.5);
    run.search_parameters
        .metadata
        .insert("empty-string".into(), "".into());
    run.search_parameters.metadata.insert(
        "strings".into(),
        vec![
            "comma,here".to_owned(),
            "".to_owned(),
            " spaced value ".to_owned(),
        ]
        .into(),
    );
    run.search_parameters
        .metadata
        .insert("ints".into(), vec![i64::MIN, 0, i64::MAX].into());
    run.search_parameters.metadata.insert(
        "floats".into(),
        MetaValue::try_from(vec![-1.25, 0.0, 1e200]).unwrap(),
    );
    run.search_parameters
        .metadata
        .insert("empty-list".into(), Vec::<String>::new().into());
    run.primary_ms_run_paths = vec!["/tmp/a&b.mzML".into(), "file://c,d.mzML".into()];
    run.raw_ms_run_paths = vec!["raw.d".into()];
    run.hits[0].rank = 2;
    run.hits[0].coverage = Some(81.25);
    run.hits[0]
        .metadata
        .insert("escaped\"&name".into(), "tab\tline\ncr\rµ".into());
    let hit = &mut doc.peptide_identifications[0].hits[0];
    hit.sequence = "(Acetyl)AC(Carbamidomethyl)M(Oxidation)K".parse().unwrap();
    hit.rank = 4;
    hit.charge = -2;
    hit.evidences[0].start = Some(0);
    hit.evidences[0].end = Some(3);
    hit.evidences[0].aa_before = FlankingResidue::NTerminus;
    hit.evidences[1].end = Some(20);
    hit.evidences[1].aa_after = FlankingResidue::CTerminus;
    hit.peak_annotations = vec![PeakAnnotation {
        mz: 120.5,
        intensity: 9.75,
        charge: 2,
        annotation: "b2|loss, \"quoted\"\\path".into(),
    }];
    hit.analysis_results = vec![
        AnalysisResult {
            score_type: "PeptideProphet".into(),
            higher_is_better: false,
            main_score: 0.1,
            sub_scores: [("fval".into(), 1.5), ("x&y".into(), 0.125)].into(),
        },
        AnalysisResult {
            score_type: "other".into(),
            higher_is_better: true,
            main_score: 0.9,
            ..Default::default()
        },
    ];
    doc.peptide_identifications.push(PeptideIdentification {
        identifier: doc.protein_identifications[1].identifier.clone(),
        ..Default::default()
    });
    doc.unreferenced_search_parameters.push(Default::default());
    doc
}

#[test]
fn pinned_whole_fixture_matches_source_goldens_and_run_links() {
    // IdXMLFile_test.cpp's original load expectations; unchanged source fixture.
    let doc = parse(WHOLE).unwrap();
    assert_eq!(doc.document_id, "LSID1234");
    assert_eq!(doc.protein_identifications.len(), 2);
    assert_eq!(doc.peptide_identifications.len(), 3);
    let p = &doc.protein_identifications[0];
    assert_eq!(
        (&*p.search_engine, &*p.search_engine_version),
        ("Mascot", "2.1.0")
    );
    assert_eq!(p.date_time.as_deref(), Some("2006-01-12T12:13:14"));
    assert_eq!(p.search_parameters.mass_type, PeakMassType::Average);
    assert_eq!(
        p.search_parameters.fragment_tolerance,
        Tolerance::Absolute(0.3)
    );
    assert_eq!(p.protein_groups[0].probability, 0.88);
    assert_eq!(p.protein_groups[0].accessions, ["PROT1", "PROT2"]);
    assert_eq!(p.indistinguishable_groups[0].accessions, ["PROT1", "PROT2"]);
    assert_eq!(p.hits[0].score, 34.4000015258789);
    assert_eq!(p.hits[0].metadata["name"], "ProteinHit".into());
    let peptide = &doc.peptide_identifications[0];
    assert_eq!(peptide.identifier, p.identifier);
    assert_eq!((peptide.rt, peptide.mz), (Some(1234.5), Some(675.9)));
    assert_eq!(peptide.spectrum_reference(), "17");
    assert_eq!(peptide.hits[0].evidences[0].protein_accession, "PROT1");
    assert_eq!(
        peptide.hits[0].evidences[1].aa_before,
        FlankingResidue::Unknown
    );
    assert_eq!(
        doc.peptide_identifications[2].identifier,
        doc.protein_identifications[1].identifier
    );
    assert_eq!(parse(&encoded(&doc)).unwrap(), doc);
}

#[test]
fn pinned_no_protein_hits_retains_flanking_evidence() {
    let doc = parse(NO_PROTEINS).unwrap();
    assert_eq!(doc.peptide_identifications.len(), 10);
    assert!(doc.protein_identifications[0].hits.is_empty());
    assert_eq!(
        doc.peptide_identifications[0].hits[0].sequence.to_string(),
        "VTAFIPPWVK"
    );
    let e = &doc.peptide_identifications[0].hits[0].evidences[0];
    assert!(e.protein_accession.is_empty());
    assert_eq!(e.aa_before, FlankingResidue::Residue('K'));
    assert_eq!(parse(&encoded(&doc)).unwrap(), doc);
}

#[test]
fn full_native_round_trip_preserves_types_modifications_ranks_and_annotations() {
    let doc = sample();
    let before = doc.clone();
    let xml = encoded(&doc);
    assert_eq!(doc, before);
    assert!(xml.contains("comma\\|here"));
    assert!(xml.contains("start=\"0 -1\""));
    assert!(xml.contains("name=\"_ar_0_higher_is_better\""));
    assert_eq!(parse(&xml).unwrap(), doc);
    assert_eq!(encoded(&parse(&xml).unwrap()), xml);
}

#[test]
fn optional_evidence_columns_source_prefix_semantics_and_duplicates() {
    let xml = WHOLE.replace("aa_before=\"A X\" aa_after=\"B X\" protein_refs=\"PH_0 PH_1\"", "aa_before=\"[ X\" aa_after=\"]\" start=\"0 20 -1\" end=\"3 23 41\" protein_refs=\"PH_0 PH_0\"");
    let doc = parse(&xml).unwrap();
    let e = &doc.peptide_identifications[0].hits[0].evidences;
    assert_eq!(e.len(), 3);
    assert_eq!(e[0].protein_accession, e[1].protein_accession);
    assert_eq!(e[2].end, Some(41));
    assert_eq!(e[2].start, None);
    assert!(e[2].protein_accession.is_empty());
    assert_eq!(parse(&encoded(&doc)).unwrap(), doc);
}

#[test]
fn no_protein_block_and_empty_runs_keep_search_context() {
    let mut doc = sample();
    doc.protein_identifications = vec![
        ProteinIdentification {
            identifier: "a".into(),
            date_time: Some("2024-02-29T10:20:30.125+02:00".into()),
            ..Default::default()
        },
        ProteinIdentification {
            identifier: "b".into(),
            date_time: Some("2024-03-01T00:00:00Z".into()),
            ..Default::default()
        },
    ];
    doc.peptide_identifications.clear();
    assert_eq!(parse(&encoded(&doc)).unwrap(), doc);
    let xml = "<IdXML version=\"1.5\"><SearchParameters id=\"SP\" db=\"\" db_version=\"\" charges=\"\" mass_type=\"average\" precursor_peak_tolerance=\"0\" peak_mass_tolerance=\"0\"/><IdentificationRun date=\"2024-02-29T00:00:00\" search_engine=\"A\" search_engine_version=\"1\" search_parameters_ref=\"SP\"/><IdentificationRun date=\"2024-02-29T00:00:00\" search_engine=\"A\" search_engine_version=\"1\" search_parameters_ref=\"SP\"/></IdXML>";
    let parsed = parse(xml).unwrap();
    assert_eq!(parsed.protein_identifications.len(), 2);
    assert_ne!(
        parsed.protein_identifications[0].identifier,
        parsed.protein_identifications[1].identifier
    );
    assert_eq!(
        parsed.protein_identifications[1]
            .search_parameters
            .mass_type,
        PeakMassType::Average
    );
}

#[test]
fn invalid_xml_structure_and_references_are_rejected() {
    for xml in [
        WHOLE.replace("PH_0 PH_1", "PH_missing PH_1"),
        WHOLE.replace(
            "search_parameters_ref=\"SP_0\"",
            "search_parameters_ref=\"missing\"",
        ),
        WHOLE.replace("id=\"PH_1\"", "id=\"PH_0\""),
        WHOLE.replace("id=\"PH_1\"", "id=\"SP_0\""),
        WHOLE.replace("0.88,PH_0,PH_1", "0.88,PH_missing"),
        WHOLE.replace("protein_group_0", "protein_group_2"),
        WHOLE.replace("version=\"1.5\"", "version=\"9\""),
        WHOLE.replace("<IdXML version", "<IdXML alien=\"1\" version"),
        WHOLE.replace(
            "<ProteinHit id=\"PH_0\"",
            "<ProteinHit id=\"PH_0\" id=\"PH_9\"",
        ),
        WHOLE.replace("</IdXML>", "<Alien/></IdXML>"),
        WHOLE.replace("<IdXML version", "<IdXML xmlns=\"urn:alien\" version"),
        WHOLE.replace("</IdXML>", ""),
        format!("{WHOLE}<IdXML/>"),
        format!("<!--before-declaration-->{WHOLE}"),
        WHOLE.replace(
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
            "<!DOCTYPE IdXML [<!ENTITY a SYSTEM \"file:///etc/passwd\">]>",
        ),
        WHOLE.replace("value=\"ProteinHit\"", "value=\"&bad;\""),
        WHOLE.replace(
            "<SearchParameters id=\"SP_0\"",
            "<SearchParameters id=\"0invalid\"",
        ),
        WHOLE.replace(
            "type=\"string\" name=\"name\"",
            "type=\"boolean\" name=\"name\"",
        ),
    ] {
        assert_ne!(xml, WHOLE, "error test must change its source fixture");
        assert!(parse(&xml).is_err(), "unexpected accepted XML {xml}");
    }
}

#[test]
fn malformed_values_lists_dates_and_malformed_definitions_fail() {
    for xml in [
        WHOLE.replace("score=\"0.9\"", "score=\"NaN\""),
        WHOLE.replace("aa_before=\"A X\"", "aa_before=\"ABC X\""),
        WHOLE.replace("aa_before=\"A X\"", "start=\"-2\""),
        WHOLE.replace("aa_before=\"A X\"", "start=\"3\" end=\"2\""),
        WHOLE.replace("2006-01-12T12:13:14", "2006-02-29T12:13:14"),
        WHOLE.replace("2006-01-12T12:13:14", "2006-01-12T12:13:14+14:01"),
        WHOLE.replace("higher_score_better=\"true\"", "higher_score_better=\"yes\""),
        WHOLE.replace("<FixedModification name=\"Fixed\" />", "<FixedModification name=\"Fixed\"><UserParam name=\"lost\" type=\"int\" value=\"1\"/></FixedModification>"),
        WHOLE.replace("<FixedModification name=\"Fixed\" />", "<FixedModification name=\"\"/>"),
        WHOLE.replace("name=\"name\" value=\"ProteinHit\"", "name=\"values\" type_extra=\"intList\" value=\"[1,2]\""),
    ] { assert!(parse(&xml).is_err()); }
    let parameter = |kind: &str, value: &str| {
        WHOLE.replacen(
            "\t</SearchParameters>",
            &format!(
                "<UserParam name=\"custom\" type=\"{kind}\" value=\"{value}\"/></SearchParameters>"
            ),
            1,
        )
    };
    for (kind, value) in [
        ("intList", "1,2"),
        ("floatList", "[NaN]"),
        ("int", "overflow"),
        ("float", "INF"),
    ] {
        assert!(parse(&parameter(kind, value)).is_err());
    }
    let definition = parameter("string", "definition")
        .replace("name=\"custom\"", "name=\"modification_definitions\"");
    assert!(parse(&definition).is_err());
}

#[test]
fn xml_attribute_normalization_preserves_character_references() {
    let xml = WHOLE.replace("value=\"ProteinHit\"", "value=\"A\r\nB\tC&#10;D&#9;E\"");
    let doc = parse(&xml).unwrap();
    assert_eq!(
        doc.protein_identifications[0].hits[0].metadata["name"]
            .as_str()
            .unwrap(),
        "A B C\nD\tE"
    );
    assert!(parse(&WHOLE.replace("value=\"ProteinHit\"", "value=\"A<B\"")).is_err());
    assert!(
        parse(&WHOLE.replace(
            " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"",
            ""
        ))
        .is_err()
    );
    assert_eq!(parse(&encoded(&doc)).unwrap(), doc);
}

/// A character reference in element text is resolved by quick-xml, so only the
/// XML 1.0 `CharRef` forms are accepted. The references here expand to
/// whitespace, which idXML allows between elements, so an accepted one leaves
/// the document unchanged and a refused one is refused for its spelling alone.
#[test]
fn text_character_references_follow_the_xml_1_0_grammar() {
    let expected = parse(WHOLE).unwrap();
    let with = |reference: &str| WHOLE.replace("</IdXML>", &format!("{reference}</IdXML>"));
    let mut wrong = Vec::new();
    for reference in ["&#32;", "&#x20;", "&#x9;", "&#0010;", "&#xD;"] {
        match parse(&with(reference)) {
            Ok(doc) if doc == expected => {}
            Ok(_) => wrong.push(format!("{reference}: changed the document")),
            Err(error) => wrong.push(format!("{reference}: {error:?}")),
        }
    }
    for reference in ["&#X20;", "&#+32;", "&#x+20;", "&#0;", "&#x0;"] {
        match parse(&with(reference)) {
            Err(openms::Error::Parse { .. }) => {}
            Err(error) => wrong.push(format!("{reference}: {error:?}")),
            Ok(_) => wrong.push(format!("{reference}: accepted")),
        }
    }
    match parse(&with("&bogus;")) {
        Err(openms::Error::Unsupported(_)) => {}
        other => wrong.push(format!("&bogus;: {:?}", other.err())),
    }
    assert!(wrong.is_empty(), "{wrong:#?}");
}

#[test]
fn reader_limits_apply_to_bytes_elements_and_lists() {
    for options in [
        ReadOptions {
            max_xml_bytes: WHOLE.len() as u64 - 1,
            ..Default::default()
        },
        ReadOptions {
            max_records: 3,
            ..Default::default()
        },
        ReadOptions {
            max_list_items: 1,
            ..Default::default()
        },
    ] {
        assert!(idxml::read_with_options(Cursor::new(WHOLE), &options).is_err());
    }
    let options = ReadOptions {
        max_xml_bytes: WHOLE.len() as u64,
        ..Default::default()
    };
    assert!(idxml::read_with_options(Cursor::new(WHOLE), &options).is_ok());
}

fn fails_before_write(doc: &IdXmlDocument) {
    let mut out = b"original".to_vec();
    assert!(idxml::write(&mut out, doc).is_err());
    assert_eq!(out, b"original");
}
#[test]
fn writer_preflight_refuses_unrepresentable_state_without_touching_destination() {
    let base = sample();
    let mut cases = Vec::new();
    let mut doc = base.clone();
    doc.protein_identifications[0].date_time = None;
    cases.push(doc);
    let mut doc = base.clone();
    doc.peptide_identifications[0].identifier = "unknown".into();
    cases.push(doc);
    let mut doc = base.clone();
    doc.peptide_identifications[0].hits[0].score = f64::INFINITY;
    cases.push(doc);
    let mut doc = base.clone();
    doc.peptide_identifications[0].hits[0].evidences[0].protein_accession = "missing".into();
    cases.push(doc);
    let mut doc = base.clone();
    doc.peptide_identifications[0].hits[0].evidences[0]
        .protein_accession
        .clear();
    cases.push(doc);
    let mut doc = base.clone();
    doc.protein_identifications[0]
        .search_parameters
        .digestion_regex = "x".into();
    cases.push(doc);
    let mut doc = base.clone();
    doc.protein_identifications[0]
        .metadata
        .insert("openms-rust:run_identifier".into(), "collision".into());
    cases.push(doc);
    for metadata in [
        MetaValue::default(),
        MetaValue::new(MetaValueData::Integer(1))
            .unwrap()
            .with_unit(Unit::new("UO:1", "a", "UO").unwrap())
            .unwrap(),
        vec!["\\|".to_owned()].into(),
        vec!["".to_owned()].into(),
    ] {
        let mut doc = base.clone();
        doc.protein_identifications[0]
            .metadata
            .insert("unsupported".into(), metadata);
        cases.push(doc);
    }
    let mut doc = base.clone();
    doc.protein_identifications[0].search_engine.push('\0');
    cases.push(doc);
    let mut doc = base.clone();
    doc.protein_identifications[0].protein_groups[0]
        .float_data_arrays
        .push(Default::default());
    cases.push(doc);
    let mut doc = base.clone();
    let duplicate = doc.protein_identifications[0].hits[0].clone();
    doc.protein_identifications[0].hits.push(duplicate);
    cases.push(doc);
    for doc in cases {
        fails_before_write(&doc);
    }
    for options in [
        WriteOptions {
            max_xml_bytes: 20,
            ..Default::default()
        },
        WriteOptions {
            max_records: 1,
            ..Default::default()
        },
    ] {
        let mut out = b"original".to_vec();
        assert!(idxml::write_with_options(&mut out, &base, &options).is_err());
        assert_eq!(out, b"original");
    }
}

#[test]
fn generated_xml_validates_against_exact_upstream_schema_when_xmllint_is_available() {
    use std::process::{Command, Stdio};
    let schema = format!("{}/tests/data/IdXML_1_5.xsd", env!("CARGO_MANIFEST_DIR"));
    let mut child = match Command::new("xmllint")
        .args(["--nonet", "--noout", "--schema", &schema, "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
    {
        Ok(child) => child,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            eprintln!("xmllint unavailable; schema validation skipped");
            return;
        }
        Err(e) => panic!("cannot run xmllint: {e}"),
    };
    child
        .stdin
        .take()
        .unwrap()
        .write_all(encoded(&sample()).as_bytes())
        .unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn buffered_sink_flush_errors_are_reported() {
    #[derive(Default)]
    struct FailingFlush {
        bytes: Vec<u8>,
    }
    impl Write for FailingFlush {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.bytes.extend_from_slice(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Err(std::io::Error::other("flush failed"))
        }
    }
    let mut sink = FailingFlush::default();
    assert!(matches!(
        idxml::write(&mut sink, &sample()),
        Err(openms::Error::Io(_))
    ));
    assert!(!sink.bytes.is_empty());
}

#[test]
fn processing_instructions_require_valid_utf8_and_xml_characters() {
    for bytes in [b"<?tool \xff?>".as_slice(), b"<?tool \0?>".as_slice()] {
        let mut xml = WHOLE.as_bytes().to_vec();
        xml.extend_from_slice(bytes);
        assert!(idxml::read(Cursor::new(xml)).is_err());
    }
}

#[test]
fn analysis_index_aliases_cannot_silently_overwrite_a_score() {
    let xml = encoded(&sample());
    for alias in ["00", "+0", " 0", "01"] {
        let alias = format!("<UserParam name=\"_ar_{alias}_score\" type=\"float\" value=\"999\"/>");
        let changed = xml.replacen("</PeptideHit>", &format!("{alias}</PeptideHit>"), 1);
        assert!(parse(&changed).is_err());
    }
}

#[test]
fn protein_group_reference_limit_is_enforced_without_peptide_evidence() {
    let mut doc = parse(WHOLE).unwrap();
    doc.peptide_identifications.clear();
    let xml = encoded(&doc);
    assert!(
        idxml::read_with_options(
            Cursor::new(&xml),
            &ReadOptions {
                max_list_items: 2,
                ..Default::default()
            }
        )
        .is_ok()
    );
    let error = idxml::read_with_options(
        Cursor::new(&xml),
        &ReadOptions {
            max_list_items: 1,
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("protein group accession list exceeds limit")
    );
}

#[test]
fn declarations_and_processing_instruction_names_follow_xml_grammar() {
    let body = WHOLE.split_once("?>").unwrap().1.trim_start();
    for declaration in [
        "<?xml version='1.0' standalone='bogus'?>",
        "<?xml version='1.0' version='1.0'?>",
        "<?xml version='1.0' encoding='UTF-8' encoding='UTF-8'?>",
        "<?xml version='1.0' alien='x'?>",
        "<?xml version='1.0' standalone='yes' encoding='UTF-8'?>",
        "<?xml version='1.0' standalone='yes' standalone='no'?>",
        "<?xml encoding='UTF-8' version='1.0'?>",
        "<?xml version='1.0'encoding='UTF-8'?>",
        "<?xml version='1.0' standalone='yes\0'?>",
    ] {
        assert!(
            parse(&format!("{declaration}{body}")).is_err(),
            "accepted {declaration}"
        );
    }
    for declaration in [
        "<?xml version='1.0'?>",
        "<?xml version = '1.0' standalone = 'yes' ?>",
        "<?xml\tversion=\"1.0\"\nencoding='utf-8'\r\nstandalone='no'?>",
    ] {
        assert!(
            parse(&format!("{declaration}{body}")).is_ok(),
            "rejected {declaration}"
        );
    }
    for pi in ["<?1bad?>", "<?XML?>", "<?Xml data?>", "<?bad!?>", "<??>"] {
        assert!(parse(&format!("{WHOLE}{pi}")).is_err(), "accepted {pi}");
    }
    for pi in [
        "<?tool?>",
        "<?tool:key opaque data?>",
        "<?:tool opaque data?>",
        "<?λtool data?>",
    ] {
        assert!(parse(&format!("{WHOLE}{pi}")).is_ok(), "rejected {pi}");
    }
}
