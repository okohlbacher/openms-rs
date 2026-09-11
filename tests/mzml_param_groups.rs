// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "mzml")]

use openms::format::mzml::{self, ReadOptions};
use openms::kernel::MSExperiment;
use std::io::Cursor;

fn read(xml: &str) -> openms::Result<MSExperiment> {
    mzml::read(Cursor::new(xml.as_bytes()))
}
fn doc(groups: &str, run: &str) -> String {
    format!(
        "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\">{groups}<run id=\"r\">{run}</run></mzML>"
    )
}
fn group(parameters: &str) -> String {
    format!(
        "<referenceableParamGroupList count=\"1\"><referenceableParamGroup id=\"g\">{parameters}</referenceableParamGroup></referenceableParamGroupList>"
    )
}
fn reference() -> &'static str {
    "<referenceableParamGroupRef ref=\"g\"/>"
}
fn spectrum(contents: &str) -> String {
    format!(
        "<spectrumList count=\"1\"><spectrum id=\"s\" defaultArrayLength=\"0\">{contents}</spectrum></spectrumList>"
    )
}

// Factoring is independent string substitution; no native reader/writer generates
// expected values. Each original self-closing parameter becomes one definition.
fn factor_inline(xml: &str) -> String {
    let mut result = String::new();
    let mut groups = String::new();
    let mut rest = xml;
    let mut count = 0;
    while let Some(start) = rest
        .find("<cvParam ")
        .into_iter()
        .chain(rest.find("<userParam "))
        .min()
    {
        result.push_str(&rest[..start]);
        let end = start + rest[start..].find("/>").unwrap() + 2;
        groups.push_str(&format!(
            "<referenceableParamGroup id=\"p{count}\">{}</referenceableParamGroup>",
            &rest[start..end]
        ));
        result.push_str(&format!("<referenceableParamGroupRef ref=\"p{count}\"/>"));
        count += 1;
        rest = &rest[end..];
    }
    result.push_str(rest);
    let list = format!(
        "<referenceableParamGroupList count=\"{count}\">{groups}</referenceableParamGroupList>"
    );
    // fileDescription is allowed to contain forward IDREFs. Definitions follow
    // that header and precede software, instrument, processing and run contexts.
    result.replacen(
        "</fileDescription>",
        &format!("</fileDescription>{list}"),
        1,
    )
}

#[test]
fn factored_independent_fixture_preserves_every_supported_parameter() {
    let original = include_str!("data/mzml_independent.mzML");
    let grouped = factor_inline(original);
    let actual = read(&grouped).unwrap();
    assert_eq!(actual, read(original).unwrap());
    assert_eq!(actual.spectra[0].rt, 90.0);
    assert_eq!(actual.spectra[0].ms_level, 2);
    assert_eq!(actual.spectra[0].precursors[0].mz, 500.25);
    assert_eq!(
        actual.spectra[0].metadata["label"].as_str().unwrap(),
        "A & B"
    );
    assert_eq!(actual.chromatograms[0].peaks[1].rt, 60.0);
    let mut encoded = Vec::new();
    mzml::write(&mut encoded, &actual).unwrap();
    assert_eq!(mzml::read(Cursor::new(encoded)).unwrap(), actual);
}

#[test]
fn source_spectrum_projection_keeps_original_ieee_arrays_and_reference() {
    let xml = include_str!("data/mzml_param_groups_source.mzML");
    // This historical binary projection retained three stale source header
    // counts. Repair only those declarations now that headers are interpreted.
    let repaired = xml
        .replacen("<sampleList count=\"1\">", "<sampleList count=\"2\">", 1)
        .replacen(
            "<softwareList count=\"3\">",
            "<softwareList count=\"4\">",
            1,
        )
        .replacen(
            "<dataProcessingList count=\"3\">",
            "<dataProcessingList count=\"4\">",
            1,
        );
    let exp = read(&repaired).unwrap();
    assert_eq!(exp.spectra.len(), 1);
    let s = &exp.spectra[0];
    assert_eq!(s.native_id, "index=0");
    assert_eq!(s.rt, 5.1);
    assert_eq!(
        s.metadata["sdname"].as_str().unwrap(),
        "spectrumdescription1"
    );
    assert_eq!(s.len(), 15);
    for (i, peak) in s.peaks.iter().enumerate() {
        assert_eq!(peak.mz.to_bits(), (i as f64).to_bits());
        assert_eq!(peak.intensity.to_bits(), ((15 - i) as f32).to_bits());
    }
    assert!(read(&xml.replace("ref=\"CommonMS1SpectrumParams\"", "ref=\"missing\"")).is_err());
}

#[test]
fn all_schema_parameter_contexts_resolve_and_unknown_references_fail() {
    let cv = "<cvParam accession=\"MS:1000127\"/><userParam name=\"ignored\" value=\"header\"/>";
    for context in [
        "fileContent",
        "sourceFile",
        "contact",
        "sample",
        "source",
        "analyzer",
        "detector",
        "instrumentConfiguration",
        "software",
        "processingMethod",
        "scanSettings",
        "target",
    ] {
        // Place each ParamGroupType in its actual header owner now that the
        // source header inventory is interpreted, including ignored values.
        let header = match context {
            "sourceFile" => format!(
                "<fileDescription><sourceFileList count=\"1\"><sourceFile id=\"f\" name=\"file\" location=\"/\">{}</sourceFile></sourceFileList></fileDescription>",
                reference()
            ),
            "instrumentConfiguration" => format!(
                "<instrumentConfigurationList count=\"1\"><instrumentConfiguration id=\"i\">{}</instrumentConfiguration></instrumentConfigurationList>",
                reference()
            ),
            "fileContent" | "contact" => format!(
                "<fileDescription><{context}>{}</{context}></fileDescription>",
                reference()
            ),
            "sample" => format!(
                "<sampleList count=\"1\"><sample id=\"sa\">{}</sample></sampleList>",
                reference()
            ),
            "source" | "analyzer" | "detector" => format!(
                "<instrumentConfigurationList count=\"1\"><instrumentConfiguration id=\"ic\"><componentList count=\"1\"><{context} order=\"1\">{}</{context}></componentList></instrumentConfiguration></instrumentConfigurationList>",
                reference()
            ),
            "software" => format!(
                "<softwareList count=\"1\"><software id=\"sw\" version=\"1\">{}</software></softwareList>",
                reference()
            ),
            "processingMethod" => format!(
                "<softwareList count=\"1\"><software id=\"sw\" version=\"1\"/></softwareList><dataProcessingList count=\"1\"><dataProcessing id=\"dp\"><processingMethod order=\"0\" softwareRef=\"sw\">{}</processingMethod></dataProcessing></dataProcessingList>",
                reference()
            ),
            _ => format!("<{context}>{}</{context}>", reference()),
        };
        let xml = doc(&format!("{}{header}", group(cv)), "");
        read(&xml).unwrap();
        assert!(
            read(&xml.replace("ref=\"g\"", "ref=\"missing\"")).is_err(),
            "{context}"
        );
    }
    for context in [
        "scanList",
        "scan",
        "scanWindow",
        "isolationWindow",
        "activation",
        "selectedIon",
    ] {
        let leaf = format!("<{context}>{}</{context}>", reference());
        let contents = match context {
            "scanList" => format!("<scanList count=\"0\">{}</scanList>", reference()),
            "scan" => format!("<scanList count=\"1\">{leaf}</scanList>"),
            "scanWindow" => format!(
                "<scanList count=\"1\"><scan><scanWindowList count=\"1\">{leaf}</scanWindowList></scan></scanList>"
            ),
            "selectedIon" => format!(
                "<precursorList count=\"1\"><precursor><selectedIonList count=\"1\">{leaf}</selectedIonList></precursor></precursorList>"
            ),
            _ => {
                format!("<precursorList count=\"1\"><precursor>{leaf}</precursor></precursorList>")
            }
        };
        let xml = doc(&group(cv), &spectrum(&contents));
        read(&xml).unwrap();
        assert!(
            read(&xml.replace("ref=\"g\"", "ref=\"missing\"")).is_err(),
            "{context}"
        );
    }
}

#[test]
fn references_apply_precursor_units_and_multiple_kinds_of_auxiliary_arrays() {
    use base64::{Engine, engine::general_purpose::STANDARD};
    let mut groups = String::new();
    let mut arrays = String::new();
    for (i, (precision, name, bytes)) in [
        ("MS:1000523", "floating", 1.25_f64.to_le_bytes().to_vec()),
        ("MS:1000519", "integer", (-123_i32).to_le_bytes().to_vec()),
        ("MS:1001479", "text", b"hello\0".to_vec()),
    ]
    .into_iter()
    .enumerate()
    {
        groups.push_str(&format!("<referenceableParamGroup id=\"a{i}\"><cvParam accession=\"{precision}\"/><cvParam accession=\"MS:1000576\"/><cvParam accession=\"MS:1000786\" value=\"{name}\"/></referenceableParamGroup>"));
        let data = STANDARD.encode(bytes);
        arrays.push_str(&format!("<binaryDataArray encodedLength=\"{}\"><referenceableParamGroupRef ref=\"a{i}\"/><binary>{data}</binary></binaryDataArray>", data.len()));
    }
    groups.push_str("<referenceableParamGroup id=\"isolation\"><cvParam accession=\"MS:1000827\" value=\"501\"/><cvParam accession=\"MS:1000828\" value=\"1.5\"/><cvParam accession=\"MS:1000829\" value=\"2.5\"/></referenceableParamGroup><referenceableParamGroup id=\"activation\"><cvParam accession=\"MS:1000133\"/><cvParam accession=\"MS:1000509\" value=\"12\" unitAccession=\"UO:0000266\"/></referenceableParamGroup>");
    let precursor = "<precursorList count=\"1\"><precursor><isolationWindow><referenceableParamGroupRef ref=\"isolation\"/></isolationWindow><activation><referenceableParamGroupRef ref=\"activation\"/></activation></precursor></precursorList>";
    let primary = "<binaryDataArray encodedLength=\"12\"><cvParam accession=\"MS:1000523\"/><cvParam accession=\"MS:1000576\"/><cvParam accession=\"MS:1000514\"/><binary>AAAAAAAAWUA=</binary></binaryDataArray><binaryDataArray encodedLength=\"8\"><cvParam accession=\"MS:1000521\"/><cvParam accession=\"MS:1000576\"/><cvParam accession=\"MS:1000515\"/><binary>AACAPw==</binary></binaryDataArray>";
    let xml = doc(
        &format!("<referenceableParamGroupList count=\"5\">{groups}</referenceableParamGroupList>"),
        &spectrum(&format!(
            "{precursor}<binaryDataArrayList count=\"5\">{primary}{arrays}</binaryDataArrayList>"
        ))
        .replace("defaultArrayLength=\"0\"", "defaultArrayLength=\"1\""),
    );
    let exp = read(&xml).unwrap();
    let s = &exp.spectra[0];
    assert_eq!(s.float_data_arrays[0].data, [1.25]);
    assert_eq!(s.integer_data_arrays[0].data, [-123]);
    assert_eq!(s.string_data_arrays[0].data, ["hello"]);
    assert_eq!(s.precursors[0].mz, 501.0);
    assert_eq!(s.precursors[0].isolation_window_lower_offset, 1.5);
    assert_eq!(s.precursors[0].activation_energy, 12.0);
    assert!(
        read(&xml.replace(
            "unitAccession=\"UO:0000266\"",
            "unitAccession=\"UO:0000031\""
        ))
        .is_err()
    );
}

#[test]
fn mixed_group_order_unicode_ids_names_and_empty_groups() {
    let parameters = "<cvParam accession=\"MS:1000511\" value=\"2\"/><userParam name=\"openms-rust:name\" value=\"A &amp; B\"/><userParam name=\"unicode\" value=\"αβ\"/>";
    let xml = doc(&group(parameters), &spectrum(reference())).replace("\"g\"", "\"组_α\"");
    let s = read(&xml).unwrap().spectra.remove(0);
    assert_eq!(s.ms_level, 2);
    assert_eq!(s.name, "A & B");
    assert_eq!(s.metadata["unicode"].as_str().unwrap(), "αβ");
    read(&doc(&group(""), &reference().repeat(10))).unwrap();
    let spaced = doc(&group(""), reference())
        .replace("id=\"g\"", "id=\" g \"")
        .replace("ref=\"g\"", "ref=\"&#9;g&#10;\"");
    read(&spaced).unwrap();
    let namespaced = xml
        .replace("<mzML ", "<m:mzML xmlns:m=\"http://psi.hupo.org/ms/mzml\" ")
        .replace("</mzML>", "</m:mzML>")
        .replace("<referenceable", "<m:referenceable")
        .replace("</referenceable", "</m:referenceable");
    assert_eq!(read(&namespaced).unwrap(), read(&xml).unwrap());
}

#[test]
fn malformed_group_structures_and_dangling_header_refs_are_errors() {
    let good = doc(&group(""), reference());
    let malformed = [
        good.replace("count=\"1\"", "count=\"2\""),
        good.replace("count=\"1\"", "count=\"0\""),
        good.replace("id=\"g\"", "id=\"\""),
        good.replace("id=\"g\"", "id=\"1bad\""),
        good.replace("id=\"g\"", "id=\"g\" id=\"h\""),
        doc(
            &group("<cvParam accession=\"MS:1000000\" accession=\"MS:1000001\"/>"),
            "",
        ),
        good.replace("ref=\"g\"", "ref=\"missing\""),
        good.replace("ref=\"g\"", "ref=\"\""),
        good.replace("ref=\"g\"", "id=\"g\""),
        good.replace(
            "<referenceableParamGroup id=\"g\">",
            "<referenceableParamGroup id=\"g\"><referenceableParamGroupRef ref=\"g\"/>",
        ),
        good.replace(
            "<referenceableParamGroup id=\"g\">",
            "<referenceableParamGroup id=\"g\"><unknown/>",
        ),
        good.replace("ref=\"g\"/>", "ref=\"g\">text</referenceableParamGroupRef>"),
        good.replace(
            "ref=\"g\"/>",
            "ref=\"g\"><userParam name=\"x\"/></referenceableParamGroupRef>",
        ),
        good.replace(
            "</referenceableParamGroupList>",
            "<referenceableParamGroup id=\"g\"/></referenceableParamGroupList>",
        )
        .replace("count=\"1\"", "count=\"2\""),
        doc("", &group("")),
        doc(
            "",
            &format!(
                "<fileDescription><fileContent>{}</fileContent></fileDescription>",
                reference()
            ),
        ),
        doc(
            &group("<userParam name=\"x\"/><cvParam accession=\"MS:1000000\"/>"),
            "",
        ),
        doc(&group("<cvParam/>"), ""),
        doc(&group("<userParam/>"), ""),
        doc(
            &group("<cvParam accession=\"MS:1000000\">bad</cvParam>"),
            "",
        ),
        doc(
            &format!(
                "<fileDescription><fileContent>{}</fileContent></fileDescription>",
                reference()
            ),
            "",
        ),
        doc(
            &group(""),
            &format!("<precursor>{}</precursor>", reference()),
        ),
    ];
    for xml in malformed {
        assert!(read(&xml).is_err(), "{xml}");
    }
    let forward = doc(
        &format!(
            "<fileDescription><fileContent>{}</fileContent></fileDescription>{}",
            reference(),
            group("<cvParam accession=\"MS:1000000\"/>")
        ),
        "",
    );
    read(&forward).unwrap();
}

#[test]
fn inline_and_referenced_conflicts_and_unsupported_binary_fields_agree() {
    let cases = [
        ("<cvParam accession=\"MS:1000511\" value=\"2\"/>", false),
        ("<userParam name=\"same\" value=\"value\"/>", false),
        ("<cvParam accession=\"MS:1002312\"/>", true),
        ("<userParam name=\"not represented\"/>", true),
        (
            "<cvParam accession=\"MS:1000523\" unitAccession=\"UO:0000010\"/>",
            true,
        ),
    ];
    for (p, binary) in cases {
        let wrap = |contents: &str| {
            if binary {
                spectrum(&format!(
                    "<binaryDataArrayList count=\"1\"><binaryDataArray encodedLength=\"0\">{contents}<binary/></binaryDataArray></binaryDataArrayList>"
                ))
            } else {
                spectrum(contents)
            }
        };
        let grouped = doc(&group(p), &wrap(&format!("{}{p}", reference())));
        let inline = doc("", &wrap(&p.repeat(2)));
        assert_eq!(
            read(&grouped).unwrap_err().to_string(),
            read(&inline).unwrap_err().to_string()
        );
    }
}

#[test]
fn definition_and_reuse_storage_and_work_are_bounded() {
    let xml = doc(&group("<userParam name=\"k\" value=\"v\"/>"), reference());
    let read_limit =
        |xml: &str, options: ReadOptions| mzml::read_with_options(Cursor::new(xml), &options);
    // Root and group-list descriptors, one group/definition/reference/application, and run.
    read_limit(
        &xml,
        ReadOptions {
            max_total_params: 7,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(
        read_limit(
            &xml,
            ReadOptions {
                max_total_params: 6,
                ..Default::default()
            }
        )
        .is_err()
    );
    assert!(
        read_limit(
            &xml,
            ReadOptions {
                max_param_groups: 0,
                ..Default::default()
            }
        )
        .is_err()
    );
    assert!(
        read_limit(
            &xml,
            ReadOptions {
                max_param_bytes: 1,
                ..Default::default()
            }
        )
        .is_err()
    );
    let empty = doc(&group(""), &reference().repeat(100));
    read_limit(
        &empty,
        ReadOptions {
            max_total_params: 104,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(
        read_limit(
            &empty,
            ReadOptions {
                max_total_params: 103,
                ..Default::default()
            }
        )
        .is_err()
    );
    let data = "x".repeat(10_000);
    let definitions = group(&format!("<userParam name=\"payload\" value=\"{data}\"/>"));
    use std::fmt::Write as _;
    let mut copies = String::new();
    for i in 0..20 {
        write!(
            &mut copies,
            "<spectrum id=\"s{i}\" defaultArrayLength=\"0\">{}</spectrum>",
            reference()
        )
        .unwrap();
    }
    let amplified = doc(
        &definitions,
        &format!("<spectrumList count=\"20\">{copies}</spectrumList>"),
    );
    assert!(
        read_limit(
            &amplified,
            ReadOptions {
                max_param_bytes: 100_000,
                ..Default::default()
            }
        )
        .is_err()
    );
    read_limit(
        &amplified,
        ReadOptions {
            max_param_bytes: 600_000,
            ..Default::default()
        },
    )
    .unwrap();
    // Even unused definitions and ignored header reference expansions consume the budget.
    assert!(
        read_limit(
            &doc(&definitions, ""),
            ReadOptions {
                max_param_bytes: 1000,
                ..Default::default()
            }
        )
        .is_err()
    );
    let header = format!(
        "<fileDescription><fileContent>{}</fileContent></fileDescription>{definitions}",
        reference().repeat(20)
    );
    assert!(
        read_limit(
            &doc(&header, ""),
            ReadOptions {
                max_param_bytes: 100_000,
                ..Default::default()
            }
        )
        .is_err()
    );
}
