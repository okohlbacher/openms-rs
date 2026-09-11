use openms::format::controlled_vocabulary::{
    CVTermDefinition, ControlledVocabulary, OboEncoding, VocabularyLimits, XRefType, fnv1a_hash,
};
use openms::metadata::{MetaValue, Unit};
use std::collections::BTreeSet;
use std::io::{self, BufRead, Read};

const SOURCE: &[u8] = include_bytes!("data/controlled_vocabulary_source.obo");
fn set(values: &[&str]) -> BTreeSet<String> {
    values.iter().map(|x| (*x).into()).collect()
}
fn parsed(text: &str) -> ControlledVocabulary {
    let mut cv = ControlledVocabulary::new();
    cv.load_obo_reader("test", text.as_bytes()).unwrap();
    cv
}
fn source() -> ControlledVocabulary {
    let mut cv = ControlledVocabulary::new();
    cv.load_obo_reader("bla", SOURCE).unwrap();
    cv
}

#[test]
fn source_fixture_all_six_literal_definitions() {
    let cv = source();
    assert_eq!(cv.name(), "bla");
    assert_eq!(cv.terms().len(), 6);
    let names = ["Auto", "Ford", "Mercedes", "A-Klasse", "Mustang", "Ka"];
    for (i, name) in names.into_iter().enumerate() {
        let id = format!("OpenMS:{}", i + 1);
        assert!(cv.exists(&id));
        let t = cv.get_term(&id).unwrap();
        assert_eq!(t.id, id);
        assert_eq!(t.name, name);
        assert_eq!(t.obsolete, i == 5);
        assert_eq!(cv.find_term_by_name(name).unwrap(), t);
        assert!(cv.has_term_with_name(name));
    }
    assert!(!cv.exists("OpenMS:7"));
    assert!(cv.get_term("OpenMS:7").is_err());
    let t = cv.get_term("OpenMS:1").unwrap();
    assert_eq!(t.description, "Auto desc");
    assert_eq!(t.synonyms, ["Kutsche", "Karre"]);
    assert!(t.parents.is_empty());
    assert!(t.unparsed.is_empty());
    assert!(!cv.has_term_with_name("Kutsche"));
    for id in ["OpenMS:2", "OpenMS:3"] {
        assert_eq!(cv.get_term(id).unwrap().parents, set(&["OpenMS:1"]));
    }
    assert_eq!(cv.get_term("OpenMS:3").unwrap().synonyms, ["Zedes"]);
    let t = cv.get_term("OpenMS:4").unwrap();
    assert_eq!(t.description, "A-Klasse desc");
    assert_eq!(t.parents, set(&["OpenMS:3"]));
    assert_eq!(
        t.unparsed,
        [
            "xref: unparsed line 1",
            "xref: unparsed line 2",
            "xref: unparsed line 3"
        ]
    );
    for id in ["OpenMS:5", "OpenMS:6"] {
        assert_eq!(cv.get_term(id).unwrap().parents, set(&["OpenMS:2"]));
    }
    assert_eq!(cv.get_term("OpenMS:6").unwrap().description, "Ka desc");
}
#[test]
fn source_ancestry_set_extension_and_checked_unknown_parent() {
    let cv = source();
    for (child, parent, expected) in [
        (6, 2, true),
        (5, 2, true),
        (2, 1, true),
        (3, 1, true),
        (4, 3, true),
        (1, 6, false),
        (4, 6, false),
        (2, 6, false),
        (2, 3, false),
    ] {
        assert_eq!(
            cv.is_child_of(&format!("OpenMS:{child}"), &format!("OpenMS:{parent}"))
                .unwrap(),
            expected
        );
    }
    assert!(cv.is_child_of("OpenMS:7", "OpenMS:3").is_err());
    assert!(!cv.is_child_of("OpenMS:1", "missing").unwrap());
    assert!(!cv.is_child_of("OpenMS:1", "OpenMS:1").unwrap());
    assert_eq!(
        cv.all_child_terms("OpenMS:2").unwrap(),
        set(&["OpenMS:5", "OpenMS:6"])
    );
    let mut output = BTreeSet::new();
    cv.add_all_child_terms(&mut output, "OpenMS:2").unwrap();
    assert_eq!(output, set(&["OpenMS:2", "OpenMS:5", "OpenMS:6"]));
    cv.add_all_child_terms(&mut output, "OpenMS:3").unwrap();
    assert_eq!(
        output,
        set(&["OpenMS:2", "OpenMS:3", "OpenMS:4", "OpenMS:5", "OpenMS:6"])
    );
    let mut untouched = set(&["sentinel"]);
    assert!(cv.add_all_child_terms(&mut untouched, "OpenMS:7").is_err());
    assert_eq!(untouched, set(&["sentinel"]));
}
#[test]
fn source_type_names_score_direction_and_xml_literals() {
    assert_eq!(
        XRefType::ALL.map(XRefType::name),
        [
            "xsd:string",
            "xsd:integer",
            "xsd:decimal",
            "xsd:negativeInteger",
            "xsd:positiveInteger",
            "xsd:nonNegativeInteger",
            "xsd:nonPositiveInteger",
            "xsd:boolean",
            "xsd:date",
            "xsd:anyURI",
            "none"
        ]
    );
    let cv = ControlledVocabulary::psi_ms().unwrap();
    for (id, expected) in [
        ("MS:1001331", true),
        ("MS:1002265", false),
        ("MS:1002467", true),
    ] {
        assert_eq!(cv.get_term(id).unwrap().is_higher_better_score(), expected);
    }
    let term = cv.get_term("MS:1001331").unwrap();
    let expected = r#"<cvParam accession="MS:1001331" cvRef="PSI-MS" name="X\!Tandem:hyperscore" value="12.5"/>"#;
    assert_eq!(term.to_xml("PSI-MS", "12.5").unwrap(), expected);
    assert_eq!(
        term.to_xml_value("PSI-MS", &MetaValue::try_from(12.5).unwrap())
            .unwrap(),
        expected
    );
    assert!(std::ptr::eq(cv, ControlledVocabulary::psi_ms().unwrap()));
    assert_eq!(cv.name(), "GO");
    assert_eq!(cv.label(), "gene_ontology");
    assert_eq!(cv.version(), "4.1.155");
}
#[test]
fn source_complete_term_copy_and_independent_field_ownership() {
    assert_eq!(CVTermDefinition::default().xref_type, XRefType::None);
    let a = CVTermDefinition {
        name: "test_cvterm".into(),
        id: "test_id".into(),
        parents: set(&["test_parent"]),
        children: set(&["test_children"]),
        obsolete: true,
        description: "test_description".into(),
        synonyms: vec!["test".into(), "synonyms".into()],
        unparsed: vec!["test".into(), "unparsed".into()],
        xref_type: XRefType::Decimal,
        xref_binary: vec!["test".into(), "xref_binary".into()],
        units: set(&["units"]),
    };
    let mut b = a.clone();
    assert_eq!(a, b);
    b.name.push('x');
    b.parents.clear();
    assert_ne!(a, b);
    assert_eq!(a.parents, set(&["test_parent"]));
    let cv = source();
    let copy = cv.checked_clone().unwrap();
    assert_eq!(cv.terms(), copy.terms());
    assert!(!std::ptr::eq(
        cv.get_term("OpenMS:1").unwrap(),
        copy.get_term("OpenMS:1").unwrap()
    ));
}
#[test]
fn source_xref_aliases_relationship_types_and_corrected_binary_prefix() {
    for (token, kind) in [
        ("string", XRefType::String),
        ("integer", XRefType::Integer),
        ("int", XRefType::Integer),
        ("decimal", XRefType::Decimal),
        ("float", XRefType::Decimal),
        ("double", XRefType::Decimal),
        ("negativeInteger", XRefType::NegativeInteger),
        ("positiveInteger", XRefType::PositiveInteger),
        ("nonNegativeInteger", XRefType::NonNegativeInteger),
        ("nonPositiveInteger", XRefType::NonPositiveInteger),
        ("boolean", XRefType::Boolean),
        ("bool", XRefType::Boolean),
        ("date", XRefType::Date),
        ("anyURI", XRefType::AnyUri),
    ] {
        for prefix in [
            "xref: value-type:",
            "xref_analog: value-type:",
            "relationship: has_value_type ",
        ] {
            let cv = parsed(&format!("[Term]\nid: X\nname: n\n{prefix}xsd:{token}\n"));
            assert_eq!(cv.get_term("X").unwrap().xref_type, kind, "{prefix}{token}");
        }
    }
    for id in ["MS:1002711", "MS:1002712", "MS:1002713"] {
        assert_eq!(
            parsed(&format!(
                "[Term]\nid: X\nrelationship: has_value_type {id}\n"
            ))
            .get_term("X")
            .unwrap()
            .xref_type,
            XRefType::String
        );
    }
    let cv = parsed(
        "[Term]\nid: X\nxref: value-type:xsd\\:decimal\nxref: binary-data-type:MS:1000523 \"float\"\nxref_analog: binary-data-type:MS:1000521\n",
    );
    let term = cv.get_term("X").unwrap();
    assert_eq!(term.xref_type, XRefType::Decimal);
    assert_eq!(term.xref_binary, ["MS:1000523", "MS:1000521"]);
    let mut cv = ControlledVocabulary::new();
    let report=cv.load_obo_reader("x",b"[Term]\nid: X\nxref: value-type:garbage xsd:string\nrelationship: has_value_type xsd\\:double\n".as_slice()).unwrap();
    assert_eq!(report.diagnostics.len(), 2);
    assert_eq!(cv.get_term("X").unwrap().xref_type, XRefType::None);
    assert!(cv.get_term("X").unwrap().unparsed.is_empty());
}
#[test]
fn source_line_stanzas_quotes_headers_and_diagnostic_timing() {
    let mut cv = ControlledVocabulary::new();
    let report=cv.load_obo_reader("first",b"data-version: one\ndefault-namespace: label\nremark: URL: https://first http://preferred\n [ tErM ] \nid: X\nname: old\nname: new\ndef: \"  quoted \\\" rest\nsynonym: no quotes\nis_obsolete:false\n[Typedef]\nname: ignored\n".as_slice()).unwrap();
    assert_eq!(report.definitions, 1);
    assert_eq!(cv.name(), "first");
    assert_eq!(cv.version(), "one");
    assert_eq!(cv.label(), "label");
    assert_eq!(cv.url(), "http://preferred");
    let t = cv.get_term("X").unwrap();
    assert_eq!(t.name, "new");
    assert_eq!(t.description, "quoted \\");
    assert_eq!(t.synonyms, ["synonym: no quotes"]);
    assert_eq!(t.unparsed, ["is_obsolete:false"]);
    let report = cv
        .load_obo_reader("second", b"remark: URL: absent\n".as_slice())
        .unwrap();
    assert_eq!(report.diagnostics[0].line, 1);
    assert_eq!(cv.name(), "second");
    assert_eq!(cv.url(), "http://preferred");
    assert_eq!(cv.version(), "one");
    let report=cv.load_obo_reader("x",b"[Term]\nid: A\nname: alpha\n[Term]\nid: B\nis_a: A ! wrong\nis_a: Z ! forward unknown\n".as_slice()).unwrap();
    assert_eq!(report.diagnostics.len(), 1);
}
#[test]
fn source_brenda_name_gate_units_and_first_last_colon_rules() {
    let text = "[Term]\nid: X\nrelationship: DRV BTO:1 ! parent\nrelationship: part_of BTO:2\nrelationship: has_units UO:0000010 ! second\n";
    let mut cv = ControlledVocabulary::new();
    cv.load_obo_reader("brenda", text.as_bytes()).unwrap();
    let t = cv.get_term("X").unwrap();
    assert_eq!(t.parents, set(&["BTO:1", "BTO:2"]));
    assert_eq!(t.units, set(&["UO:0000010"]));
    let mut cv = ControlledVocabulary::new();
    cv.load_obo_reader("BTO", text.as_bytes()).unwrap();
    let t = cv.get_term("X").unwrap();
    assert!(t.parents.is_empty());
    assert_eq!(t.unparsed.len(), 2);
    let cv = parsed("[Term]\nid: X\nrelationship: has_units UO:1:2 ! label:3\n");
    assert_eq!(cv.get_term("X").unwrap().units, set(&["UO:3"]));
}
#[test]
fn source_cumulative_stale_indexes_and_placeholder_iteration_order() {
    let mut cv = parsed("[Term]\nid: B\nname: before\nis_a: A\n");
    assert!(cv.exists("A"));
    assert_eq!(cv.get_term("A").unwrap().id, "");
    assert!(!cv.has_term_with_name(""));
    cv.load_obo_reader("x", b"[Term]\nid: B\nname: after\n".as_slice())
        .unwrap();
    assert!(cv.has_term_with_name(""));
    assert_eq!(cv.find_term_by_name("before").unwrap().name, "after");
    assert!(cv.get_term("A").unwrap().children.contains("B"));
    assert!(cv.get_term("B").unwrap().parents.is_empty());
    let cv = parsed("[Term]\nid: A\nname: a\nis_a: Z\n");
    assert!(cv.has_term_with_name(""));
    assert_eq!(cv.find_term_by_name("").unwrap().id, "");
    let cv = parsed(
        "[Term]\nid: B\nname: same\ndef: \"second\"\n[Term]\nid: A\nname: same\ndef: \"first\"\n",
    );
    assert_eq!(cv.get_term_by_name("same", "second").unwrap().id, "A");
    assert_eq!(cv.get_term_by_name("sam", "esecond").unwrap().id, "B");
}
#[test]
fn independent_diamond_callbacks_early_stop_and_active_cycle_checks() {
    let cv = parsed(
        "[Term]\nid: R\nname: match\n[Term]\nid: B\nis_a: R\n[Term]\nid: C\nis_a: R\n[Term]\nid: D\nname: match\nis_a: B\nis_a: C\n",
    );
    let mut order = Vec::new();
    assert!(
        !cv.iterate_all_children("R", |id| {
            order.push(id.to_owned());
            false
        })
        .unwrap()
    );
    assert_eq!(order, ["B", "D", "C", "D"]);
    assert_eq!(
        cv.first_child_with_name("R", "match").unwrap().unwrap().id,
        "D"
    );
    let cv = parsed("[Term]\nid: A\nis_a: B\n[Term]\nid: B\nis_a: A\n");
    assert!(cv.iterate_all_children("A", |id| id == "B").unwrap());
    assert!(cv.iterate_all_children("A", |_| false).is_err());
    assert!(cv.is_child_of("A", "B").unwrap());
    assert!(cv.is_child_of("A", "missing").is_err());
    let mut output = set(&["sentinel"]);
    assert!(cv.extend_child_terms(&mut output, "A").is_err());
    assert_eq!(output, set(&["sentinel"]));
}
#[test]
fn corrected_typed_unit_identity_xml_escaping_and_empty_distinction() {
    let t = CVTermDefinition {
        id: "X:\"a&".into(),
        name: "a<µ\n\t\r'".into(),
        units: set(&["UO:0000016", "UO:0000017"]),
        ..Default::default()
    };
    let xml = t.to_xml("R\"", "v<&").unwrap();
    assert!(xml.contains("accession=\"X:&quot;a&amp;\""));
    assert!(xml.contains("cvRef=\"R&quot;\""));
    assert!(xml.contains("&#10;&#9;&#13;&apos;"));
    assert!(!t.to_xml("x", "").unwrap().contains(" value="));
    assert!(
        !t.to_xml_value("x", &MetaValue::default())
            .unwrap()
            .contains(" value=")
    );
    assert!(
        t.to_xml_value("x", &"".into())
            .unwrap()
            .contains(" value=\"\"")
    );
    let val = MetaValue::try_from(5.0)
        .unwrap()
        .with_unit(Unit::new("UO:0000017", "micrometer", "UO").unwrap())
        .unwrap();
    assert!(
        t.to_xml_value("x", &val)
            .unwrap()
            .contains("unitAccession=\"UO:0000017\" unitCvRef=\"UO\" unitName=\"micrometer\"")
    );
    let mut t = t;
    t.units.clear();
    assert!(t.to_xml_value("x", &val).is_ok());
    t.description = "opaque\0description".into();
    assert!(t.to_xml("x", "ok").is_ok());
    t.name.push('\0');
    assert!(t.to_xml("x", "ok").is_err());
}

#[test]
fn typed_xml_lists_and_large_integers_keep_source_spelling_without_narrowing() {
    let term = CVTermDefinition {
        id: "X:1".into(),
        name: "test".into(),
        ..Default::default()
    };
    // Source ListUtilsIO.h emits unquoted, comma-space-separated lists. The
    // full native i64 domain remains intact instead of converting via i32.
    for (value, expected) in [
        (
            MetaValue::from(vec!["a,b".to_string(), String::new(), "<\"".into()]),
            "[a,b, , &lt;&quot;]",
        ),
        (MetaValue::from(vec![-5i64, 7]), "[-5, 7]"),
        (
            MetaValue::try_from(vec![1.5f64, 2.0]).unwrap(),
            "[1.5, 2.0]",
        ),
        (MetaValue::from(Vec::<String>::new()), "[]"),
        (MetaValue::from(i64::MAX), "9223372036854775807"),
        (MetaValue::from(i64::MIN), "-9223372036854775808"),
        (MetaValue::from(vec![i64::MAX]), "[9223372036854775807]"),
    ] {
        assert_eq!(
            term.to_xml_value("MS", &value).unwrap(),
            format!("<cvParam accession=\"X:1\" cvRef=\"MS\" name=\"test\" value=\"{expected}\"/>")
        );
    }
}
#[test]
fn bounded_legacy_decoding_and_original_bto_corruption_are_explicit() {
    let raw = b"[Term]\nid: X\nname: caf\xe9\ndef: \"one\0two \x92\"\n";
    let mut cv = ControlledVocabulary::new();
    assert!(cv.load_obo_reader("x", raw.as_slice()).is_err());
    assert!(cv.terms().is_empty());
    cv.load_obo_encoded("x", raw.as_slice(), OboEncoding::Windows1252)
        .unwrap();
    assert_eq!(cv.get_term("X").unwrap().name, "café");
    assert_eq!(cv.get_term("X").unwrap().description, "one\0two ’");
    assert!(
        cv.load_obo_encoded(
            "x",
            b"[Term]\nid: Y\nname: \x81\n".as_slice(),
            OboEncoding::Windows1252
        )
        .is_err()
    );
    assert!(!cv.exists("Y"));
    let t = ControlledVocabulary::psi_ms()
        .unwrap()
        .get_term("BTO:0002243")
        .unwrap();
    assert_eq!(t.name, "hypanthium");
    assert_eq!(t.description.bytes().filter(|b| *b == 0).count(), 8);
}
#[test]
fn diagnostic_output_uses_only_selected_writer_and_source_fnv_bytes() {
    let cv = parsed("[Term]\nid: A\nname: alpha\n[Term]\nid: B\nname: beta\nis_a: A\n");
    let expected = "[Term]\nid: 'A'\nname: 'alpha'\n[Term]\nid: 'B'\nname: 'beta'\nis_a: 'A'\n";
    assert_eq!(cv.to_text().unwrap(), expected);
    let mut out = Vec::new();
    cv.write_text(&mut out).unwrap();
    assert_eq!(out, expected.as_bytes());
    assert_eq!(fnv1a_hash(""), 0xcbf29ce484222325);
    assert_eq!(fnv1a_hash("hello"), 0xa430d84680aabd0b);
    let expected = "µ".as_bytes().iter().fold(0xcbf29ce484222325u64, |x, b| {
        (x ^ u64::from(*b)).wrapping_mul(0x100000001b3)
    });
    assert_eq!(fnv1a_hash("µ"), expected);
}

struct Tiny<'a> {
    bytes: &'a [u8],
    chunk: usize,
    fail: bool,
}
impl Read for Tiny<'_> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        let data = self.fill_buf()?;
        let n = data.len().min(out.len());
        out[..n].copy_from_slice(&data[..n]);
        self.consume(n);
        Ok(n)
    }
}
impl BufRead for Tiny<'_> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.bytes.is_empty() && self.fail {
            return Err(io::Error::other("late failure"));
        }
        Ok(&self.bytes[..self.bytes.len().min(self.chunk)])
    }
    fn consume(&mut self, n: usize) {
        self.bytes = &self.bytes[n..];
    }
}
#[test]
fn tiny_reader_limits_and_io_failure_leave_old_registry_unchanged() {
    for chunk in 1..=7 {
        let mut cv = ControlledVocabulary::new();
        cv.load_obo_reader(
            "bla",
            Tiny {
                bytes: SOURCE,
                chunk,
                fail: false,
            },
        )
        .unwrap();
        assert_eq!(cv.terms(), source().terms());
    }
    let mut cv = source();
    let old = cv.terms().clone();
    assert!(
        cv.load_obo_reader(
            "new",
            Tiny {
                bytes: b"[Term]\nid: Z\n",
                chunk: 2,
                fail: true
            }
        )
        .is_err()
    );
    assert_eq!(cv.terms(), &old);
    assert_eq!(cv.name(), "bla");
    for limits in [
        VocabularyLimits {
            max_input_bytes: 4,
            ..Default::default()
        },
        VocabularyLimits {
            max_line_bytes: 4,
            ..Default::default()
        },
        VocabularyLimits {
            max_terms: 1,
            ..Default::default()
        },
        VocabularyLimits {
            max_entries: 1,
            ..Default::default()
        },
        VocabularyLimits {
            max_work: 1,
            ..Default::default()
        },
        VocabularyLimits {
            max_bytes: 1,
            ..Default::default()
        },
    ] {
        let mut cv = ControlledVocabulary::with_limits(limits);
        assert!(cv.load_obo_reader("new", SOURCE).is_err());
        assert!(cv.terms().is_empty());
        assert_eq!(cv.name(), "");
    }
    let t = CVTermDefinition::default();
    assert!(
        t.to_xml_with_limits(
            "x",
            "long",
            VocabularyLimits {
                max_output_bytes: 1,
                ..Default::default()
            }
        )
        .is_err()
    );
}

fn decode_hex(value: &str) -> String {
    assert_eq!(value.len() % 2, 0);
    let bytes = (0..value.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&value[i..i + 2], 16).unwrap())
        .collect();
    String::from_utf8(bytes).unwrap()
}
fn decode_list(value: &str) -> Vec<String> {
    let (count, text) = value.split_once(':').unwrap();
    let count = count.parse::<usize>().unwrap();
    let result = if count == 0 {
        Vec::new()
    } else {
        text.split(',').map(decode_hex).collect()
    };
    assert_eq!(count, result.len());
    result
}
#[test]
fn every_provider_field_and_all_stale_name_aliases_match_independent_source_projection() {
    let cv = ControlledVocabulary::psi_ms().unwrap();
    let rows = include_str!("data/controlled_vocabulary_projection.tsv")
        .lines()
        .skip(1);
    assert_eq!(rows.clone().count(), 9254);
    assert_eq!(cv.terms().len(), 9254);
    for (line, (key, term)) in rows.zip(cv.terms()) {
        let f: Vec<_> = line.split('\t').collect();
        assert_eq!(f.len(), 12);
        assert_eq!(key, &decode_hex(f[0]));
        let expected = CVTermDefinition {
            id: decode_hex(f[1]),
            name: decode_hex(f[2]),
            description: decode_hex(f[3]),
            parents: decode_list(f[4]).into_iter().collect(),
            children: decode_list(f[5]).into_iter().collect(),
            obsolete: f[6] == "1",
            xref_type: XRefType::ALL[f[7].parse::<usize>().unwrap()],
            synonyms: decode_list(f[8]),
            unparsed: decode_list(f[9]),
            xref_binary: decode_list(f[10]),
            units: decode_list(f[11]).into_iter().collect(),
        };
        assert_eq!(term, &expected, "{key}");
    }
    let aliases = include_str!("data/controlled_vocabulary_aliases.tsv")
        .lines()
        .skip(1);
    assert_eq!(aliases.clone().count(), 16852);
    for row in aliases {
        let (name, id) = row.split_once('\t').unwrap();
        let (name, id) = (decode_hex(name), decode_hex(id));
        assert_eq!(
            cv.find_term_by_name(&name).unwrap().id,
            id,
            "alias {name:?}"
        );
    }
}
#[test]
fn original_ontology_bytes_are_identical_on_all_checkout_platforms() {
    // Independent Python source-file hashes, plus SHA-256 in provenance/generator.
    for (raw, bytes, expected) in [
        (
            include_bytes!("../resources/cv/psi-ms.obo").as_slice(),
            985486,
            1880820938405340563u64,
        ),
        (
            include_bytes!("../resources/cv/quality.obo").as_slice(),
            305675,
            17963641850810773392,
        ),
        (
            include_bytes!("../resources/cv/unit.obo").as_slice(),
            89140,
            17816829084039043071,
        ),
        (
            include_bytes!("../resources/cv/brenda.obo").as_slice(),
            932031,
            15642686075072363256,
        ),
        (
            include_bytes!("../resources/cv/goslim_goa.obo").as_slice(),
            28100,
            14078192894626321412,
        ),
    ] {
        assert_eq!(raw.len(), bytes);
        let hash = raw.iter().fold(14695981039346656037u64, |x, b| {
            (x ^ u64::from(*b)).wrapping_mul(1099511628211)
        });
        assert_eq!(hash, expected);
    }
}
#[test]
fn deep_acyclic_graph_uses_no_recursive_call_stack() {
    let mut text = String::from("[Term]\nid: N00000\n");
    for i in 1..3000 {
        text.push_str(&format!("[Term]\nid: N{i:05}\nis_a: N{:05}\n", i - 1));
    }
    let cv = parsed(&text);
    assert!(cv.is_child_of("N02999", "N00000").unwrap());
    assert_eq!(cv.all_child_terms("N00000").unwrap().len(), 2999);
}
