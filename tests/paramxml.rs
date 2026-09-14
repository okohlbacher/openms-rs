// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "paramxml")]

use openms::format::paramxml::{self, Limits, OutputEncoding, WriteOptions};
use openms::param::{Param, ParamValue};
use std::io::{self, Read, Write};

const SOURCE: &[u8] = include_bytes!("data/paramxml/ParamXMLFile_test_writeXMLToStream.xml");
fn xml(body: &str) -> String {
    format!("<PARAMETERS>{body}</PARAMETERS>")
}
fn parse(body: &str) -> Param {
    paramxml::read(xml(body).as_bytes()).unwrap()
}
fn copy(param: &Param) -> Param {
    let mut output = Vec::new();
    paramxml::write(&mut output, param).unwrap();
    paramxml::read(output.as_slice()).unwrap()
}

#[test]
fn source_writer_fixture_all_values_and_metadata() {
    let param = paramxml::read(SOURCE).unwrap();
    assert_eq!(param.root().entries.len(), 15);
    for (name, value) in [
        (
            "stringlist",
            ParamValue::StringList(vec!["a".into(), "bb".into(), "ccc".into()]),
        ),
        ("intlist", ParamValue::IntegerList(vec![1, 22, 333])),
        ("item", ParamValue::String("bla".into())),
        ("stringlist2", ParamValue::StringList(vec![])),
        ("intlist2", ParamValue::IntegerList(vec![])),
        ("item1", ParamValue::Integer(7)),
        ("intlist3", ParamValue::IntegerList(vec![1])),
        ("stringlist3", ParamValue::StringList(vec!["1".into()])),
        ("item3", ParamValue::Float(7.6)),
        ("doublelist", ParamValue::FloatList(vec![1.22, 2.33, 4.55])),
        ("doublelist3", ParamValue::FloatList(vec![1.4])),
    ] {
        assert_eq!(param.entry(name).unwrap().value, value);
    }
    assert_eq!(
        param.entry("stringlist").unwrap().description,
        "StringList Description"
    );
    assert!(
        param
            .entry("file_parameter")
            .unwrap()
            .tags
            .contains("input file")
    );
    assert_eq!(
        param.entry("file_parameter").unwrap().valid_strings,
        ["*.mzML", "*.mzXML"]
    );
    assert!(
        param
            .entry("advanced_parameter")
            .unwrap()
            .tags
            .contains("advanced")
    );
    assert_eq!(
        param.entry("flag").unwrap().valid_strings,
        ["true", "false"]
    );
    let mut output = Vec::new();
    paramxml::write(&mut output, &param).unwrap();
    // This immutable source golden differs only in the intentionally honest UTF-8 declaration.
    assert_eq!(
        String::from_utf8(output).unwrap(),
        std::str::from_utf8(SOURCE)
            .unwrap()
            .replace("ISO-8859-1", "UTF-8")
    );
    assert_eq!(param, copy(&param));
}

/// Writer option `WriteOptions::source()`: the source's `ISO-8859-1`
/// declaration. The immutable source golden is then reproduced byte for byte,
/// declaration included, while the default writer keeps declaring UTF-8.
/// Characters above U+007F are written consistently with the declaration: up to
/// U+00FF as one ISO-8859-1 byte, beyond as a character reference, and both
/// read back unchanged. The source writer copies UTF-8 bytes under the same
/// declaration instead, so its non-ASCII text does not read back; that native
/// difference is recorded in `docs/PARAMXML_SUPPORT.md`.
#[test]
fn source_declaration_writer_option() {
    let param = paramxml::read(SOURCE).unwrap();
    let mut output = Vec::new();
    paramxml::write_with_options(
        &mut output,
        &param,
        Limits::default(),
        WriteOptions::source(),
    )
    .unwrap();
    assert_eq!(output, SOURCE);
    assert_eq!(WriteOptions::default().encoding, OutputEncoding::Utf8);
    assert_eq!(WriteOptions::source().encoding, OutputEncoding::Latin1);

    let mut param = Param::new();
    param
        .set_value(
            "tool:1:name",
            ParamValue::String("caf\u{e9} \u{20ac} \u{1d11e}".into()),
            "na\u{ef}ve \u{2603}",
            &[],
        )
        .unwrap();
    param
        .set_value(
            "tool:1:list",
            ParamValue::StringList(vec!["\u{fc}".into(), "\u{3b1}".into()]),
            "",
            &[],
        )
        .unwrap();
    let mut latin1 = Vec::new();
    paramxml::write_with_options(
        &mut latin1,
        &param,
        Limits::default(),
        WriteOptions::source(),
    )
    .unwrap();
    assert!(latin1.starts_with(b"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"));
    assert!(std::str::from_utf8(&latin1).is_err());
    let contains = |needle: &[u8]| latin1.windows(needle.len()).any(|w| w == needle);
    assert!(contains(b"value=\"caf\xe9 &#x20AC; &#x1D11E;\""));
    assert!(contains(b"description=\"na\xefve &#x2603;\""));
    assert!(contains(b"<LISTITEM value=\"\xfc\"/>"));
    assert!(contains(b"<LISTITEM value=\"&#x3B1;\"/>"));
    assert_eq!(paramxml::read(latin1.as_slice()).unwrap(), param);

    let mut utf8 = Vec::new();
    paramxml::write(&mut utf8, &param).unwrap();
    assert!(utf8.starts_with(b"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"));
    assert_eq!(paramxml::read(utf8.as_slice()).unwrap(), param);

    let dir = openms::system::file::TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("source.ini");
    paramxml::store_with_options(&path, &param, WriteOptions::source()).unwrap();
    assert_eq!(std::fs::read(&path).unwrap(), latin1);
    assert_eq!(paramxml::load(&path).unwrap(), param);
}

#[test]
fn source_legacy_tag_migration_and_optional_required_attribute() {
    for source in [
        include_bytes!("data/paramxml/Param_pre16_update.ini").as_slice(),
        include_bytes!("data/paramxml/Param_post16_update.ini").as_slice(),
    ] {
        let param = paramxml::read(source).unwrap();
        let prefix = "SpectraFilterMarkerMower:1:";
        for (key, required, advanced, file) in [
            ("in", true, false, Some("input file")),
            ("out", true, false, Some("output file")),
            ("log", false, true, None),
            ("no_progress", false, true, None),
        ] {
            let entry = param.entry(&format!("{prefix}{key}")).unwrap();
            assert_eq!(entry.tags.contains("required"), required);
            assert_eq!(entry.tags.contains("advanced"), advanced);
            if let Some(tag) = file {
                assert!(entry.tags.contains(tag));
            }
        }
        assert_eq!(param, copy(&param));
    }
    let param =
        paramxml::read(include_bytes!("data/paramxml/Param_advanced_no_required.ini").as_slice())
            .unwrap();
    for (key, advanced, required) in [
        ("adv_true_req_absent", true, false),
        ("adv_true_req_false", true, false),
        ("adv_true_req_true", true, true),
        ("adv_absent_req_true", false, true),
        ("adv_absent_req_absent", false, false),
        ("adv_false_req_absent", false, false),
        ("list_adv_true_req_absent", true, false),
    ] {
        let entry = param.entry(&format!("TagMatrix:{key}")).unwrap();
        assert_eq!(entry.tags.contains("advanced"), advanced);
        assert_eq!(entry.tags.contains("required"), required);
    }
    assert_eq!(param, copy(&param));
}

#[test]
fn source_ranges_including_legacy_hyphens_and_all_list_alternatives() {
    let param = parse(
        r#"
      <ITEM name="a" type="int" value="5" restrictions="-9:10"/>
      <ITEM name="b" type="int" value="5" restrictions="4-6"/>
      <ITEM name="c" type="float" value="5.1" restrictions=":6.1"/>
      <ITEM name="d" type="double" value="5.1" restrictions="4.1:"/>
      <ITEM name="ignored" type="int" value="5" restrictions=""/>
      <ITEMLIST name="ints" type="int" restrictions="1:11"><LISTITEM value="2"/><LISTITEM value="5"/><LISTITEM value="10"/></ITEMLIST>
      <ITEMLIST name="doubles" type="double" restrictions="0.1:5.8"><LISTITEM value="1.2"/><LISTITEM value="3.33"/><LISTITEM value="4.44"/></ITEMLIST>
      <ITEMLIST name="strings" type="string" restrictions="xml,txt"><LISTITEM value="a.txt"/><LISTITEM value="b.xml"/></ITEMLIST>
    "#,
    );
    assert_eq!(
        (
            param.entry("a").unwrap().min_int,
            param.entry("a").unwrap().max_int
        ),
        (-9, 10)
    );
    assert_eq!(
        (
            param.entry("b").unwrap().min_int,
            param.entry("b").unwrap().max_int
        ),
        (4, 6)
    );
    assert_eq!(param.entry("c").unwrap().min_float, -f64::MAX);
    assert_eq!(param.entry("c").unwrap().max_float, 6.1);
    assert_eq!(param.entry("d").unwrap().min_float, 4.1);
    assert_eq!(param.entry("d").unwrap().max_float, f64::MAX);
    assert_eq!(param.entry("ignored").unwrap().min_int, -i32::MAX);
    assert_eq!(
        param.entry("ints").unwrap().value,
        ParamValue::IntegerList(vec![2, 5, 10])
    );
    assert_eq!(param.entry("ints").unwrap().min_int, 1);
    assert_eq!(
        param.entry("strings").unwrap().valid_strings,
        ["xml", "txt"]
    );
    assert_eq!(param, copy(&param));
}

#[test]
fn source_file_types_flags_and_supported_format_precedence() {
    let param = parse(
        r#"
      <ITEM name="in" type="input-file" value="a.mzML" tags="custom,advanced" required="true" supported_formats="*.mzML"/>
      <ITEM name="out" type="output-file" value="b.mzML" supported_formats="*.mzML"/>
      <ITEM name="prefix" type="output-prefix" value="p" supported_formats="*.mzML"/>
      <ITEM name="legacy" type="string" value="c" tags="input file" restrictions="ignored" supported_formats="*.mzXML"/>
      <ITEMLIST name="files" type="input-file" supported_formats="*.mzML"><LISTITEM value="a"/><LISTITEM value="b"/></ITEMLIST>
      <ITEMLIST name="outputs" type="output-file" supported_formats="*.csv"/>
      <ITEM name="flag" type="bool" value="false" restrictions="ignored"/>
      <ITEM name="true" type="string" value="true" restrictions="true,false"/>
    "#,
    );
    assert_eq!(param.entry("legacy").unwrap().valid_strings, ["*.mzXML"]);
    assert_eq!(param.entry("prefix").unwrap().valid_strings, ["*.mzML"]);
    assert_eq!(param.entry("files").unwrap().valid_strings, ["*.mzML"]);
    assert_eq!(
        param.entry("flag").unwrap().valid_strings,
        ["true", "false"]
    );
    assert_eq!(param, copy(&param));
}

#[test]
fn precision_nan_infinities_and_signed_zero_are_preserved() {
    let param = parse(
        r#"
      <ITEM name="nan" type="double" value="NaN"/>
      <ITEM name="inf" type="double" value="INF"/>
      <ITEM name="minus_inf" type="double" value="-INF"/>
      <ITEM name="precise" type="double" value="1.0000000000000002"/>
      <ITEM name="zero" type="double" value="-0.0"/>
      <ITEMLIST name="list" type="double"><LISTITEM value="NaN"/><LISTITEM value="INF"/><LISTITEM value="-0.0"/></ITEMLIST>
    "#,
    );
    for p in [&param, &copy(&param)] {
        assert!(p.entry("nan").unwrap().value.to_f64().unwrap().is_nan());
        assert_eq!(
            p.entry("inf").unwrap().value,
            ParamValue::Float(f64::INFINITY)
        );
        assert_eq!(
            p.entry("minus_inf").unwrap().value,
            ParamValue::Float(f64::NEG_INFINITY)
        );
        assert_eq!(
            p.entry("precise").unwrap().value,
            ParamValue::Float(1.0000000000000002)
        );
        assert!(
            p.entry("zero")
                .unwrap()
                .value
                .to_f64()
                .unwrap()
                .is_sign_negative()
        );
        let list = p.entry("list").unwrap().value.as_float_list().unwrap();
        assert!(list[0].is_nan());
        assert!(list[2].is_sign_negative());
    }
    let mut output = Vec::new();
    paramxml::write(&mut output, &param).unwrap();
    assert!(String::from_utf8(output).unwrap().contains("value=\"NaN\""));
}

#[test]
fn xml_escaping_unicode_descriptions_empty_sections_and_literal_whitespace() {
    let param = parse(
        "<NODE name=\"n&amp;é\" description=\"a#br#b&amp;c\"><ITEM name=\"x\" value=\"a&#x9;b&#xA;c&#xD;d&amp;&lt;&gt;&apos;&quot;\" type=\"string\" description=\"quoted &quot; value\"/></NODE><NODE name=\"empty\" description=\"kept\"/>",
    );
    assert_eq!(
        param.entry("n&é:x").unwrap().value.as_str().unwrap(),
        "a\tb\nc\rd&<>'\""
    );
    assert_eq!(param.root().nodes[0].description, "a\nb&c");
    assert_eq!(param, copy(&param));
    let normalized = parse("<ITEM name=\"x\" type=\"string\" value=\"a\r\nb\tc\"/>");
    assert_eq!(
        normalized.entry("x").unwrap().value.as_str().unwrap(),
        "a b c"
    );
    assert_eq!(Param::new(), copy(&Param::new()));
}

#[test]
fn source_latin1_and_native_utf16_encodings() {
    let mut latin = b"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?><PARAMETERS><ITEM name=\"x\" type=\"string\" value=\"caf".to_vec();
    latin.extend_from_slice(b"\xe9\"/></PARAMETERS>");
    assert_eq!(
        paramxml::read(latin.as_slice())
            .unwrap()
            .entry("x")
            .unwrap()
            .value
            .as_str()
            .unwrap(),
        "café"
    );
    for little in [false, true] {
        let text = "<?xml version=\"1.0\" encoding=\"UTF-16\"?><PARAMETERS><ITEM name=\"x\" type=\"string\" value=\"é𝄞\"/></PARAMETERS>";
        let mut bytes = if little {
            vec![0xff, 0xfe]
        } else {
            vec![0xfe, 0xff]
        };
        for unit in text.encode_utf16() {
            bytes.extend_from_slice(&if little {
                unit.to_le_bytes()
            } else {
                unit.to_be_bytes()
            });
        }
        assert_eq!(
            paramxml::read(bytes.as_slice())
                .unwrap()
                .entry("x")
                .unwrap()
                .value
                .as_str()
                .unwrap(),
            "é𝄞"
        );
    }
    let mut bom = vec![0xef, 0xbb, 0xbf];
    bom.extend_from_slice(&latin);
    assert!(paramxml::read(bom.as_slice()).is_err());
    assert!(paramxml::read([0xff, 0xfe, 0].as_slice()).is_err());
}

#[test]
fn load_into_accumulates_updates_source_metadata_and_rolls_back_errors() {
    let mut target = parse(
        r#"<ITEM name="old" value="3" type="int" description="keep" restrictions="1:10" tags="advanced"/><ITEM name="untouched" value="x" type="string"/>"#,
    );
    paramxml::read_into(
        xml(r#"<ITEM name="old" value="7" type="int"/>"#).as_bytes(),
        &mut target,
        Limits::default(),
    )
    .unwrap();
    assert_eq!(target.entry("old").unwrap().value, ParamValue::Integer(7));
    assert_eq!(target.entry("old").unwrap().description, "keep");
    assert!(target.entry("old").unwrap().tags.is_empty());
    assert_eq!(target.entry("old").unwrap().min_int, 1);
    assert!(target.entry("untouched").is_ok());
    let before = target.clone();
    let broken =
        xml(r#"<ITEM name="old" value="9" type="int"/><ITEM name="bad" value="xxx" type="int"/>"#);
    assert!(paramxml::read_into(broken.as_bytes(), &mut target, Limits::default()).is_err());
    assert_eq!(target, before);
}

#[test]
fn malformed_xml_unknown_types_and_external_entities_are_errors() {
    for input in [
        "",
        "<PARAMETERS>",
        "<PARAMETERS/><PARAMETERS/>",
        "<PARAMETERS><NODE name=\"n\"></PARAMETERS>",
        "<PARAMETERS version=\"2.0\"/>",
        "<?xml version=\"1.0\" version=\"1.0\"?><PARAMETERS/>",
        "<?xml version=\"1.0\" standalone=\"maybe\"?><PARAMETERS/>",
        "<?XML version=\"1.0\"?><PARAMETERS/>",
        "<PARAMETERS><!-- bad -- comment --></PARAMETERS>",
        "<PARAMETERS><LISTITEM value=\"1\"/></PARAMETERS>",
        "<!DOCTYPE PARAMETERS [<!ENTITY x SYSTEM 'file:///etc/passwd'>]><PARAMETERS/>",
        "<?xml version=\"1.0\"?><?xml version=\"1.0\"?><PARAMETERS/>",
        " <\u{0}PARAMETERS/>",
        "<PARAMETERS xsi:noNamespaceSchemaLocation=\"x\"/>",
    ] {
        assert!(
            paramxml::read(input.as_bytes()).is_err(),
            "accepted {input:?}"
        );
    }
    for item in [
        r#"<ITEM name="x" name="y" value="v" type="string"/>"#,
        r#"<ITEM name="x" value="v" type="unknown"/>"#,
        r#"<ITEM name="x" type="string"/>"#,
        r#"<ITEM name="x" value="&missing;" type="string"/>"#,
        r#"<ITEM name="x" value="&#0;" type="string"/>"#,
        r#"<ITEM name="x" value="a<b" type="string"/>"#,
        r#"<ITEM name="x" value="2147483648" type="int"/>"#,
        r#"<ITEMLIST name="x" type="bool"/>"#,
        r#"<ITEM name="" value="v" type="string"/>"#,
    ] {
        assert!(
            paramxml::read(xml(item).as_bytes()).is_err(),
            "accepted {item}"
        );
    }
}

#[test]
fn review_regressions_cover_attribute_amplification_and_unrepresentable_drafts() {
    use openms::param::{ParamEntry, ParamNode};
    for text in [
        "<?1bad?><PARAMETERS/>",
        "<PARAMETERS><ITEM name='x'type='int'value='1'/></PARAMETERS>",
        "<?xml version='1.0'encoding='UTF-8'?><PARAMETERS/>",
        "<?xmlversion='1.0'?><PARAMETERS/>",
        "<PARAMETERS><ITEMLIST name='x' type='string' value='silently lost'/></PARAMETERS>",
        "<PARAMETERS><ITEMLIST name='x' type='string' default='lost'/></PARAMETERS>",
        "<PARAMETERS><ITEM name='x' type='double' value='1e9999'/></PARAMETERS>",
        "<PARAMETERS><ITEM name='x' type='double' value='2' restrictions='1e9999:'/></PARAMETERS>",
    ] {
        assert!(paramxml::read(text.as_bytes()).is_err(), "accepted {text}");
    }
    let tags = xml(&format!(
        "<ITEM name='x' type='string' value='a' tags='{}'/>",
        ",".repeat(100)
    ));
    for attribute in [tags.clone(), tags.replace("tags=", "restrictions=")] {
        assert!(
            paramxml::read_with_limits(
                attribute.as_bytes(),
                Limits {
                    max_list_items: 100,
                    ..Default::default()
                }
            )
            .is_err()
        );
        // A byte budget can be ample for XML but insufficient for owned field slots.
        assert!(
            paramxml::read_with_limits(
                attribute.as_bytes(),
                Limits {
                    max_xml_bytes: 512,
                    ..Default::default()
                }
            )
            .is_err()
        );
    }
    assert!(paramxml::read("<PARAMETERS unsupported='x' another='y'/>".as_bytes()).is_err());
    for name in ["", "nested:leaf"] {
        let param = Param::from_root(ParamNode {
            entries: vec![ParamEntry {
                name: name.into(),
                value: ParamValue::Integer(1),
                ..Default::default()
            }],
            name: "ROOT".into(),
            description: String::new(),
            nodes: Vec::new(),
        })
        .unwrap();
        let mut out = Vec::new();
        assert!(paramxml::write(&mut out, &param).is_err());
        assert!(out.is_empty());
    }
    let entry = ParamEntry {
        name: "duplicate".into(),
        value: ParamValue::Integer(1),
        ..Default::default()
    };
    let param = Param::from_root(ParamNode {
        entries: vec![entry.clone(), entry],
        name: "ROOT".into(),
        description: String::new(),
        nodes: Vec::new(),
    })
    .unwrap();
    assert!(paramxml::write(Vec::new(), &param).is_err());
    let param = Param::from_root(ParamNode {
        entries: vec![ParamEntry {
            name: "unused".into(),
            value: ParamValue::String("x".into()),
            min_int: 0,
            ..Default::default()
        }],
        name: "ROOT".into(),
        description: String::new(),
        nodes: Vec::new(),
    })
    .unwrap();
    assert!(paramxml::write(Vec::new(), &param).is_err());
}

#[test]
fn limits_and_unrepresentable_values_fail_before_external_writes() {
    let text = xml(
        r#"<NODE name="long"><ITEMLIST name="list" type="int"><LISTITEM value="1"/><LISTITEM value="2"/></ITEMLIST></NODE>"#,
    );
    let base = Limits::default();
    for limits in [
        Limits {
            max_xml_bytes: text.len() - 1,
            ..base
        },
        Limits {
            max_elements: 4,
            ..base
        },
        Limits {
            max_depth: 3,
            ..base
        },
        Limits {
            max_list_items: 1,
            ..base
        },
        Limits {
            max_path_bytes: 8,
            ..base
        },
    ] {
        assert!(paramxml::read_with_limits(text.as_bytes(), limits).is_err());
        let param = paramxml::read(text.as_bytes()).unwrap();
        let mut out = b"original".to_vec();
        assert!(paramxml::write_with_limits(&mut out, &param, limits).is_err());
        assert_eq!(out, b"original");
    }
    for (value, desc) in [
        (ParamValue::Empty, ""),
        (ParamValue::Integer(i64::MAX), ""),
        (ParamValue::String("bad\0".into()), ""),
        (ParamValue::Integer(3), "literal#br#text"),
    ] {
        let mut param = Param::new();
        param.set_value("bad", value, desc, &[]).unwrap();
        let mut out = b"original".to_vec();
        assert!(paramxml::write(&mut out, &param).is_err());
        assert_eq!(out, b"original");
    }
}

#[test]
fn reader_and_flush_errors_propagate() {
    struct BrokenRead;
    impl Read for BrokenRead {
        fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
            Err(io::Error::other("read"))
        }
    }
    struct BrokenFlush;
    impl Write for BrokenFlush {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            Ok(bytes.len())
        }
        fn flush(&mut self) -> io::Result<()> {
            Err(io::Error::other("flush"))
        }
    }
    assert!(paramxml::read(BrokenRead).is_err());
    assert!(paramxml::write(BrokenFlush, &Param::new()).is_err());
}

#[test]
fn ini_parameters_drive_a_native_filter_and_peak_list_roundtrip() {
    use openms::format::{FileHandler, FileType};
    use openms::processing::{SpectrumFilter, ThresholdMower};
    let parameters = parse(
        r#"<NODE name="ThresholdMower"><ITEM name="threshold" value="5" type="double" restrictions="0:"/></NODE>"#,
    );
    let filter = ThresholdMower {
        threshold: parameters
            .value("ThresholdMower:threshold")
            .unwrap()
            .to_f64()
            .unwrap(),
    };
    let mut experiment = FileHandler::read_experiment(
        b"#SEC\tMZ\tINT\n1 100 4\n1 200 5\n1 300 9\n".as_slice(),
        FileType::Dta2d,
    )
    .unwrap();
    filter.filter_experiment(&mut experiment).unwrap();
    let mut output = Vec::new();
    FileHandler::write_experiment(&mut output, &experiment, FileType::Dta2d).unwrap();
    let copy = FileHandler::read_experiment(output.as_slice(), FileType::Dta2d).unwrap();
    assert_eq!(
        copy.spectra[0]
            .peaks
            .iter()
            .map(|p| (p.mz, p.intensity))
            .collect::<Vec<_>>(),
        [(200.0, 5.0), (300.0, 9.0)]
    );
}

#[test]
fn stored_file_validates_against_the_original_schema_when_available() {
    let path =
        std::env::temp_dir().join(format!("openms-paramxml-schema-{}.xml", std::process::id()));
    let param = paramxml::read(SOURCE).unwrap();
    paramxml::store(&path, &param).unwrap();
    assert_eq!(param, paramxml::load(&path).unwrap());
    let checked = std::process::Command::new("xmllint")
        .args(["--nonet", "--noout", "--schema"])
        .arg(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/data/paramxml/Param_1_8_0.xsd"
        ))
        .arg(&path)
        .output();
    std::fs::remove_file(&path).unwrap();
    match checked {
        Ok(result) => assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        ),
        Err(error) if error.kind() == io::ErrorKind::NotFound => {
            eprintln!("xmllint unavailable; schema check skipped")
        }
        Err(error) => panic!("schema validator failed: {error}"),
    }
}

#[test]
fn compressed_paths_detect_magic_preserve_plain_output_and_fail_atomically() {
    let directory = openms::system::file::TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = directory.path().join("parameters.unknown");
    let expected = paramxml::read(SOURCE).unwrap();
    let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gzip.write_all(SOURCE).unwrap();
    let gzip = gzip.finish().unwrap();
    let mut bzip = bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::default());
    bzip.write_all(SOURCE).unwrap();
    let bzip = bzip.finish().unwrap();
    for encoded in [&gzip, &bzip] {
        std::fs::write(&path, encoded).unwrap();
        assert_eq!(paramxml::load(&path).unwrap(), expected);
        let mut target = parse("<ITEM name=\"retained\" value=\"present\" type=\"string\"/>");
        paramxml::load_into(&path, &mut target).unwrap();
        assert_eq!(
            target.entry("retained").unwrap().value,
            ParamValue::from("present")
        );
        assert_eq!(target.entry("item1").unwrap().value, ParamValue::Integer(7));
        let before = target.clone();
        assert!(
            paramxml::load_into_with_limits(
                &path,
                &mut target,
                Limits {
                    max_xml_bytes: SOURCE.len() - 1,
                    ..Default::default()
                }
            )
            .is_err()
        );
        assert_eq!(target, before);
        std::fs::write(&path, &encoded[..encoded.len() - 8]).unwrap();
        assert!(paramxml::load_into(&path, &mut target).is_err());
        assert_eq!(target, before);
    }
    let output = directory.path().join("plain.ini.gz");
    paramxml::store(&output, &expected).unwrap();
    assert!(std::fs::read(&output).unwrap().starts_with(b"<?xml"));
    assert_eq!(paramxml::load(output).unwrap(), expected);
}
