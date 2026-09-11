#![cfg(feature = "proforma-json")]
use openms::chemistry::proforma::*;
use serde_json::{Value, json};
use std::sync::Arc;

fn empty() -> Value {
    json!({"global_mods":[],"unlocalised_mods":[],"labile_mods":[],"n_term_mods":[],"sequence":[],"c_term_mods":[]})
}
fn read(v: Value) -> Peptidoform {
    Peptidoform::from_json(&v.to_string()).unwrap()
}
fn section(aa: &str, mods: Value) -> Value {
    json!({"type":"element","value":{"amino_acid":aa,"modifications":mods}})
}
fn named(name: &str) -> Value {
    json!({"tag":{"type":"named_mod","value":{"name":name}}})
}
fn with_tag(kind: &str, value: Value) -> Value {
    let mut v = empty();
    v["n_term_mods"] = json!([[{"tag":{"type":kind,"value":value}}]]);
    v
}
fn tag(p: &Peptidoform) -> &ModificationTag {
    &p.n_term_mods[0].alternatives[0].0
}

#[test]
fn two_source_class_roundtrip_examples_and_exact_empty_schema() {
    // Source ProFormaParser_test.cpp1306–1332; source only asserts lengths.
    let peptide = Peptidoform::parse("EM[UNIMOD:35]K").unwrap();
    let encoded = peptide.to_json().unwrap();
    assert!(!encoded.is_empty());
    assert_eq!(Peptidoform::from_json(&encoded).unwrap(), peptide);
    let ion = PeptidoformIon::parse("PEPTIDE//SEQUENCE/2").unwrap();
    let encoded = ion.to_json().unwrap();
    assert!(!encoded.is_empty());
    assert_eq!(PeptidoformIon::from_json(&encoded).unwrap(), ion);
    // Independently derived exact compact lexical object representation.
    assert_eq!(
        Peptidoform::default().to_json().unwrap(),
        r#"{"c_term_mods":[],"global_mods":[],"labile_mods":[],"n_term_mods":[],"sequence":[],"unlocalised_mods":[]}"#
    );
    assert_eq!(
        PeptidoformIon::default().to_json().unwrap(),
        r#"{"chains":[],"is_chimeric":false}"#
    );
}

#[test]
fn every_source_positive_grammar_ast_roundtrips_through_json() {
    let mut count = 0;
    for text in include_str!("data/proforma_parser_positive.txt")
        .lines()
        .filter(|s| !s.is_empty() && !s.starts_with('#'))
    {
        if let Ok(p) = Peptidoform::parse(text) {
            assert_eq!(
                Peptidoform::from_json(&p.to_json().unwrap()).unwrap(),
                p,
                "{text}"
            );
        } else {
            let p = PeptidoformIon::parse(text).unwrap();
            assert_eq!(
                PeptidoformIon::from_json(&p.to_json().unwrap()).unwrap(),
                p,
                "{text}"
            );
        }
        count += 1;
    }
    assert_eq!(count, 176);
}

#[test]
fn independent_complete_schema_retains_every_reachable_variant() {
    let tags = json!([
        {"tag":{"type":"cv_accession","value":{"database":"UNIMOD","accession":"0035"}},"label":{"type":"CROSSLINK","identifier":"XL1","score":0.25}},
        {"tag":{"type":"named_mod","value":{"name":"Oxidation α","cv_hint":"MOD"}},"label":{"type":"BRANCH","identifier":"BR2"}},
        {"tag":{"type":"mass_delta","value":{"source":"OBS","mass":1.25,"original_text":"+1.2500"}},"label":{"type":"AMBIGUOUS","identifier":"g3","score":-0.5}},
        {"tag":{"type":"formula","value":{"formula_string":"[13C]H-1","charge":-2}}},
        {"tag":{"type":"glycan","value":[{"monosaccharide":{"type":"name","value":"Hex"},"count":-3},{"monosaccharide":{"type":"formula","value":{"formula_string":"C2","charge":0}},"count":0}]}},
        {"tag":{"type":"info","value":{"text":"δ\nannotation"}}},
        {"tag":{"type":"position","value":{"residues":"ACCA?\u{0}","n_term":true,"c_term":false}}}
    ]);
    let el = section("?", json!([tags.clone()]));
    let chain = json!({
        "name":"α chain",
        "global_mods":[{"type":"isotope_replacement","value":{"isotope":"13C"}},{"type":"global_modification","value":{"modification":tags.clone(),"locations":["K","N-term","K"]}}],
        "unlocalised_mods":[{"modifications":[tags.clone()],"occurrence":-2}],
        "labile_mods":[{"modification":tags.clone()}],
        "n_term_mods":[tags.clone()],
        "sequence":[el.clone(),{"type":"ambiguous_region","value":{"elements":[el["value"].clone()]}},{"type":"modified_range","value":{"elements":[el["value"].clone()],"modifications":[tags.clone()]}}],
        "c_term_mods":[tags],
        "charge":{"type":"adducts","value":[{"formula":"Na","charge":1,"occurrence":2},{"formula":"H","charge":-1}]}
    });
    let p = read(chain.clone());
    assert_eq!(p.sequence.len(), 3);
    assert_eq!(p.global_mods.len(), 2);
    let alts = &p.n_term_mods[0].alternatives;
    assert_eq!(alts.len(), 7);
    assert!(matches!(&alts[0].0,ModificationTag::CvAccession(v) if v.accession=="0035"));
    assert!(
        matches!(&alts[1].0,ModificationTag::NamedMod(v) if v.name=="Oxidation α"&&v.cv_hint==Some(CvDatabase::Mod))
    );
    assert!(
        matches!(&alts[2].0,ModificationTag::MassDelta(v) if v.mass==1.25&&v.source==MassDeltaSource::Obs)
    );
    assert!(matches!(&alts[3].0,ModificationTag::FormulaTag(v) if v.charge==Some(-2)));
    assert!(matches!(&alts[4].0,ModificationTag::GlycanComposition(v) if v.components.len()==2));
    assert!(matches!(&alts[5].0,ModificationTag::InfoTag(v) if v.text=="δ\nannotation"));
    assert!(
        matches!(&alts[6].0,ModificationTag::PositionConstraint(v) if v.residues==['A','C','C','A','?','\0']&&v.n_term&&!v.c_term)
    );
    assert_eq!(
        serde_json::from_str::<Value>(&p.to_json().unwrap()).unwrap(),
        chain
    );
    let expected = json!({"chains":[chain],"name":"Ion name retained","is_chimeric":true,"charge":{"type":"simple","value":-5}});
    let ion = PeptidoformIon::from_json(&expected.to_string()).unwrap();
    assert_eq!(
        serde_json::from_str::<Value>(&ion.to_json().unwrap()).unwrap(),
        expected
    );
}

#[test]
fn required_vectors_and_optional_null_have_distinct_source_rules() {
    let mut v = empty();
    v.as_object_mut().unwrap().remove("global_mods");
    v["name"] = Value::Null;
    v["charge"] = Value::Null;
    assert_eq!(read(v), Peptidoform::default());
    for key in [
        "unlocalised_mods",
        "labile_mods",
        "n_term_mods",
        "sequence",
        "c_term_mods",
    ] {
        let mut v = empty();
        v.as_object_mut().unwrap().remove(key);
        assert!(Peptidoform::from_json(&v.to_string()).is_err(), "{key}");
    }
    for key in [
        "global_mods",
        "unlocalised_mods",
        "labile_mods",
        "n_term_mods",
        "sequence",
        "c_term_mods",
    ] {
        for bad in [Value::Null, json!({}), json!(0)] {
            let mut v = empty();
            v[key] = bad;
            assert!(Peptidoform::from_json(&v.to_string()).is_err(), "{key}");
        }
    }
    assert_eq!(
        PeptidoformIon::from_json(r#"{"chains":[],"name":null,"charge":null}"#).unwrap(),
        PeptidoformIon::default()
    );
    for text in [
        r#"{"chains":null}"#,
        r#"{"chains":[],"is_chimeric":null}"#,
        r#"{"chains":[],"is_chimeric":1}"#,
        r#"{}"#,
    ] {
        assert!(PeptidoformIon::from_json(text).is_err());
    }
    let p = read(with_tag("position", json!({"residues":""})));
    assert!(matches!(tag(&p),ModificationTag::PositionConstraint(v) if !v.n_term&&!v.c_term));
    assert!(
        Peptidoform::from_json(
            &with_tag("position", json!({"residues":"","n_term":null})).to_string()
        )
        .is_err()
    );
}

#[test]
fn source_explicit_iteration_accepts_objects_null_and_scalar_locations() {
    let mut v = empty();
    v["n_term_mods"] = json!([null,{}, {"z":named("last"),"a":named("first")}]);
    let p = read(v);
    assert!(p.n_term_mods[0].alternatives.is_empty());
    assert!(p.n_term_mods[1].alternatives.is_empty());
    assert!(
        matches!(&p.n_term_mods[2].alternatives[0].0,ModificationTag::NamedMod(v) if v.name=="first")
    );
    assert!(
        matches!(&p.n_term_mods[2].alternatives[1].0,ModificationTag::NamedMod(v) if v.name=="last")
    );
    for value in [Value::Null, json!({})] {
        let p = read(with_tag("glycan", value));
        assert!(matches!(tag(&p),ModificationTag::GlycanComposition(v) if v.components.is_empty()));
    }
    let component = |name: &str| json!({"monosaccharide":{"type":"name","value":name},"count":1});
    let p = read(with_tag(
        "glycan",
        json!({"z":component("Z"),"a":component("A")}),
    ));
    assert!(
        matches!(tag(&p),ModificationTag::GlycanComposition(v) if v.components[0].0==GlycanComponent::Name("A".into()))
    );
    for (locations, expected) in [
        (Value::Null, vec![]),
        (json!({"z":"C","a":"N"}), vec!["N", "C"]),
        (json!("K"), vec!["K"]),
    ] {
        let mut v = empty();
        v["global_mods"] = json!([{"type":"global_modification","value":{"modification":null,"locations":locations}}]);
        let p = read(v);
        let GlobalModEntry::GlobalModification(g) = &p.global_mods[0] else {
            panic!()
        };
        assert_eq!(g.locations, expected);
    }
}

#[test]
fn all_enum_spellings_are_exact_and_unknown_discriminators_error() {
    for (name, db) in [
        ("UNIMOD", CvDatabase::Unimod),
        ("MOD", CvDatabase::Mod),
        ("RESID", CvDatabase::Resid),
        ("XLMOD", CvDatabase::Xlmod),
        ("GNO", CvDatabase::Gno),
    ] {
        let v = with_tag("cv_accession", json!({"database":name,"accession":"id"}));
        let p = read(v.clone());
        assert!(matches!(tag(&p),ModificationTag::CvAccession(v) if v.database==db));
        assert_eq!(
            serde_json::from_str::<Value>(&p.to_json().unwrap()).unwrap(),
            v
        );
    }
    for (name, source) in [
        ("NONE", MassDeltaSource::None),
        ("OBS", MassDeltaSource::Obs),
        ("U", MassDeltaSource::U),
        ("M", MassDeltaSource::M),
        ("R", MassDeltaSource::R),
        ("X", MassDeltaSource::X),
        ("G", MassDeltaSource::G),
    ] {
        let p = read(with_tag(
            "mass_delta",
            json!({"source":name,"mass":0.0,"original_text":""}),
        ));
        assert!(matches!(tag(&p),ModificationTag::MassDelta(v) if v.source==source));
    }
    for v in [
        with_tag("CV_accession", json!({})),
        with_tag(
            "cv_accession",
            json!({"database":"unimod","accession":"35"}),
        ),
        with_tag(
            "mass_delta",
            json!({"source":"Obs","mass":1,"original_text":""}),
        ),
        with_tag(
            "glycan",
            json!([{"monosaccharide":{"type":"NAME","value":"Hex"},"count":1}]),
        ),
    ] {
        assert!(Peptidoform::from_json(&v.to_string()).is_err());
    }
    let mut v = empty();
    v["sequence"] = json!([{"type":"unknown","value":{}}]);
    assert!(Peptidoform::from_json(&v.to_string()).is_err());
}

#[test]
fn integer_boolean_and_fraction_conversion_matches_source_int_overload() {
    for (input, expected) in [
        (json!(true), 1),
        (json!(false), 0),
        (json!(2.99), 2),
        (json!(-2.99), -2),
        (json!(2147483647.9), i32::MAX),
        (json!(-2147483648.9), i32::MIN),
    ] {
        let mut v = empty();
        v["charge"] = json!({"type":"simple","value":input});
        assert_eq!(read(v).charge, Some(ChargeState::Simple(expected)));
    }
    for input in [
        json!(2147483648u64),
        json!(-2147483649i64),
        json!(1e100),
        json!("2"),
        Value::Null,
    ] {
        let mut v = empty();
        v["charge"] = json!({"type":"simple","value":input});
        assert!(Peptidoform::from_json(&v.to_string()).is_err());
    }
    assert!(
        Peptidoform::from_json(
            &with_tag(
                "mass_delta",
                json!({"source":"NONE","mass":true,"original_text":""})
            )
            .to_string()
        )
        .is_err()
    );
}

#[test]
fn float_bits_signed_zero_bom_and_json_exponents() {
    for value in [
        0.0,
        -0.0,
        f64::from_bits(1),
        f64::MIN_POSITIVE,
        1.2345678901234567,
        f64::MAX,
        -f64::MAX,
    ] {
        let p = Peptidoform {
            n_term_mods: vec![Modification {
                alternatives: vec![(
                    ModificationTag::MassDelta(MassDelta {
                        mass: value,
                        ..Default::default()
                    }),
                    None,
                )],
                ..Default::default()
            }],
            ..Default::default()
        };
        let q = Peptidoform::from_json(&p.to_json().unwrap()).unwrap();
        let ModificationTag::MassDelta(m) = tag(&q) else {
            panic!()
        };
        assert_eq!(m.mass.to_bits(), value.to_bits());
    }
    for (token, bits) in [
        ("-0", 0),
        ("-0.0", (-0.0f64).to_bits()),
        ("-0e2", (-0.0f64).to_bits()),
        ("1e-0", 1.0f64.to_bits()),
        ("1E-0", 1.0f64.to_bits()),
        ("-0e-0", (-0.0f64).to_bits()),
        ("-1e-0", (-1.0f64).to_bits()),
        ("1e-400", 0),
        ("-1e-400", (-0.0f64).to_bits()),
    ] {
        let v = with_tag(
            "mass_delta",
            json!({"source":"NONE","mass":"TOKEN","original_text":"-0"}),
        )
        .to_string()
        .replace("\"TOKEN\"", token);
        let q = Peptidoform::from_json(&v).unwrap();
        let ModificationTag::MassDelta(m) = tag(&q) else {
            panic!()
        };
        assert_eq!(m.mass.to_bits(), bits, "{token}");
        assert_eq!(m.original_text, "-0");
    }
    for token in ["1-0", "1e--0", "1e -0", "e-0", "true-0", "-0-0"] {
        let text = with_tag(
            "mass_delta",
            json!({"source":"NONE","mass":"TOKEN","original_text":""}),
        )
        .to_string()
        .replace("\"TOKEN\"", token);
        assert!(Peptidoform::from_json(&text).is_err(), "{token}");
    }
    assert_eq!(
        Peptidoform::from_json(&format!("\u{feff}{}", empty())).unwrap(),
        Peptidoform::default()
    );
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let p = Peptidoform {
            n_term_mods: vec![Modification {
                alternatives: vec![(
                    ModificationTag::MassDelta(MassDelta {
                        mass: value,
                        ..Default::default()
                    }),
                    None,
                )],
                ..Default::default()
            }],
            ..Default::default()
        };
        assert!(p.to_json().is_err());
    }
}

#[test]
fn duplicate_keys_last_wins_and_unknown_fields_are_ignored() {
    let mut s = empty().to_string();
    s.pop();
    s.push_str(r#", "name":"first", "name":"last", "resolved_mod":{"bad":"pointer"}, "unknown":[{"anything":null}]}"#);
    let p = Peptidoform::from_json(&s).unwrap();
    assert_eq!(p.name.as_deref(), Some("last"));
    assert!(!p.to_json().unwrap().contains("unknown"));
    for invalid in [
        "",
        "null",
        "[]",
        "{",
        r#"{"chains":[],}"#,
        r#"{"chains":[]}//comment"#,
        r#"{"chains":[],"ignored":1e400}"#,
    ] {
        assert!(PeptidoformIon::from_json(invalid).is_err(), "{invalid}");
    }
}

#[test]
fn character_fields_allow_ascii_bytes_and_preserve_unicode_strings() {
    let mut v = empty();
    v["sequence"] = json!([section("\0", json!([])), section("?", json!([]))]);
    v["name"] = json!("αβ😀\n\"\\");
    let p = read(v.clone());
    assert_eq!(
        serde_json::from_str::<Value>(&p.to_json().unwrap()).unwrap(),
        v
    );
    for aa in ["", "AA", "é", "😀"] {
        let mut v = empty();
        v["sequence"] = json!([section(aa, json!([]))]);
        assert!(Peptidoform::from_json(&v.to_string()).is_err());
    }
    assert!(
        Peptidoform::from_json(&with_tag("position", json!({"residues":"é"})).to_string()).is_err()
    );
    let p = Peptidoform {
        sequence: vec![SequenceSection::Element(SequenceElement {
            amino_acid: 'é',
            modifications: vec![],
        })],
        ..Default::default()
    };
    assert!(p.to_json().is_err());
}

#[test]
fn resolved_handle_is_omitted_and_not_reconstituted() {
    let modification = openms::chemistry::ModificationsDB::global()
        .get_modification_handle(
            "Oxidation",
            Some('M'),
            Some(openms::chemistry::TermSpecificity::Anywhere),
        )
        .unwrap();
    let mut p = Peptidoform::parse("M[Oxidation]").unwrap();
    let SequenceSection::Element(el) = &mut p.sequence[0] else {
        panic!()
    };
    el.modifications[0].resolved_mod = Some(Arc::clone(&modification));
    let encoded = p.to_json().unwrap();
    assert!(!encoded.contains("resolved"));
    let q = Peptidoform::from_json(&encoded).unwrap();
    let SequenceSection::Element(el) = &q.sequence[0] else {
        panic!()
    };
    assert!(el.modifications[0].resolved_mod.is_none());
    let SequenceSection::Element(el) = &p.sequence[0] else {
        panic!()
    };
    assert!(Arc::ptr_eq(
        el.modifications[0].resolved_mod.as_ref().unwrap(),
        &modification
    ));
}

#[test]
fn depth_size_output_and_late_field_failures_are_checked() {
    assert!(Peptidoform::from_json(&" ".repeat(MAX_PROFORMA_JSON_TEXT_BYTES + 1)).is_err());
    let mut v = empty().to_string();
    v.pop();
    v.push_str(",\"ignored\":");
    v.push_str(&"[".repeat(MAX_PROFORMA_JSON_DEPTH));
    v.push_str("null");
    v.push_str(&"]".repeat(MAX_PROFORMA_JSON_DEPTH));
    v.push('}');
    assert!(Peptidoform::from_json(&v).is_err());
    // Large ignored keys still participate in structural/preallocation bounds.
    assert!(
        Peptidoform::from_json(&format!("{{\"ignored\":[{}]}}", "0,".repeat(400_000))).is_err()
    );
    let p = Peptidoform {
        name: Some("\0".repeat(MAX_PROFORMA_JSON_TEXT_BYTES / 2)),
        ..Default::default()
    };
    let saved = p.clone();
    assert!(p.to_json().is_err());
    assert_eq!(p, saved);
    let ion = PeptidoformIon {
        chains: vec![
            Peptidoform::default(),
            Peptidoform {
                sequence: vec![SequenceSection::Element(SequenceElement {
                    amino_acid: 'é',
                    modifications: vec![],
                })],
                ..Default::default()
            },
        ],
        ..Default::default()
    };
    let saved = ion.clone();
    assert!(ion.to_json().is_err());
    assert_eq!(ion, saved);
}
