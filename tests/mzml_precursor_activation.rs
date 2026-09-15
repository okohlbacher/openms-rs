// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "mzml")]
use openms::format::mzml::{self, ReadOptions};
use openms::metadata::{ActivationMethod as A, MetaValue, Unit};
use openms::{MSChromatogram, MSExperiment, MSSpectrum, Precursor};
use std::io::Cursor;

fn cv(id: &str, value: &str, attrs: &str) -> String {
    format!(r#"<cvParam accession="{id}" name="fixture" value="{value}" {attrs}/>"#)
}
fn document(selected: &str, activation: &str, chrom: bool) -> String {
    let precursor = format!(
        r#"<precursor><selectedIonList count="1"><selectedIon>{selected}</selectedIon></selectedIonList><activation>{activation}</activation></precursor>"#
    );
    let record = if chrom {
        format!(
            r#"<chromatogramList count="1"><chromatogram id="c" defaultArrayLength="0">{precursor}</chromatogram></chromatogramList>"#
        )
    } else {
        format!(
            r#"<spectrumList count="1"><spectrum id="scan=1" defaultArrayLength="0"><precursorList count="1">{precursor}</precursorList></spectrum></spectrumList>"#
        )
    };
    format!(
        r#"<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0"><run id="r">{record}</run></mzML>"#
    )
}
fn extract(mut e: MSExperiment, chrom: bool) -> Precursor {
    if chrom {
        e.chromatograms.remove(0).precursor
    } else {
        e.spectra.remove(0).precursors.remove(0)
    }
}
fn read(xml: &str, chrom: bool) -> Precursor {
    extract(mzml::read(Cursor::new(xml)).unwrap(), chrom)
}
fn native(p: Precursor, chrom: bool) -> MSExperiment {
    if chrom {
        MSExperiment {
            chromatograms: vec![MSChromatogram {
                native_id: "c".into(),
                precursor: p,
                ..Default::default()
            }],
            ..Default::default()
        }
    } else {
        MSExperiment {
            spectra: vec![MSSpectrum {
                native_id: "scan=1".into(),
                precursors: vec![p],
                ..Default::default()
            }],
            ..Default::default()
        }
    }
}
fn roundtrip(p: &Precursor, chrom: bool) -> String {
    let mut output = Vec::new();
    mzml::write(&mut output, &native(p.clone(), chrom)).unwrap();
    let xml = String::from_utf8(output).unwrap();
    assert_eq!(&read(&xml, chrom), p);
    xml
}
const EV: &str = r#"unitAccession="UO:0000266" unitCvRef="UO" unitName="electronvolt""#;

#[test]
fn source_class_energy_intensity_and_supplemental_literals_roundtrip() {
    // MzMLFile_test.cpp:1420-1434,1484-1485. Source-reviewed literals,
    // not a newly executed C++ result. Source reading derives the combined method.
    for chrom in [false, true] {
        let selected = cv("MS:1000042", "30", r#"unitAccession="MS:1000131""#);
        let activation = [
            cv("MS:1000598", "", ""),
            cv("MS:1002678", "", ""),
            cv("MS:1000045", "25", EV),
            cv("MS:1002680", "25", EV),
        ]
        .concat();
        let p = read(&document(&selected, &activation, chrom), chrom);
        assert_eq!(p.intensity, 30.);
        assert_eq!(p.activation_methods, [A::Etd, A::Ethcd].into());
        assert_eq!(
            p.cv_terms.metadata["peak intensity unit accession"]
                .as_str()
                .unwrap(),
            "MS:1000131"
        );
        let energy = MetaValue::try_from(25.)
            .unwrap()
            .with_unit(Unit::new("UO:0000266", "electronvolt", "UO").unwrap())
            .unwrap();
        assert_eq!(p.cv_terms.metadata["collision energy"], energy);
        assert_eq!(p.cv_terms.metadata["supplemental collision energy"], energy);
        let xml = roundtrip(&p, chrom);
        assert!(xml.contains("MS:1002678"));
        assert!(!xml.contains("MS:1002631"));
        assert!(xml.contains("unitAccession=\"MS:1000131\""));
    }
}

#[test]
fn all_nine_activation_metadata_routes_preserve_scalar_types() {
    let rows = [
        (
            "MS:1000245",
            "charge stripping",
            "",
            MetaValue::from("true"),
        ),
        (
            "MS:1000045",
            "collision energy",
            "12.5",
            MetaValue::try_from(12.5).unwrap(),
        ),
        (
            "MS:1000412",
            "buffer gas",
            "helium",
            MetaValue::from("helium"),
        ),
        (
            "MS:1000419",
            "collision gas",
            "nitrogen",
            MetaValue::from("nitrogen"),
        ),
        (
            "MS:1000138",
            "percent collision energy",
            "30",
            MetaValue::try_from(30.).unwrap(),
        ),
        (
            "MS:1000869",
            "collision gas pressure",
            "0.5",
            MetaValue::try_from(0.5).unwrap(),
        ),
        (
            "MS:1002679",
            "supplemental collision-induced dissociation",
            "",
            MetaValue::from(""),
        ),
        (
            "MS:1002678",
            "supplemental beam-type collision-induced dissociation",
            "",
            MetaValue::from(""),
        ),
        (
            "MS:1002680",
            "supplemental collision energy",
            "5",
            MetaValue::try_from(5.).unwrap(),
        ),
    ];
    for chrom in [false, true] {
        for (id, key, text, expected) in &rows {
            let p = read(&document("", &cv(id, text, ""), chrom), chrom);
            assert_eq!(&p.cv_terms.metadata[*key], expected, "{id}");
            assert!(roundtrip(&p, chrom).contains(id));
        }
    }
}

#[test]
fn typed_activation_energy_is_distinct_from_collision_metadata() {
    let p = read(
        &document(
            "",
            &[cv("MS:1000509", "11", EV), cv("MS:1000045", "25", EV)].concat(),
            false,
        ),
        false,
    );
    assert_eq!(p.activation_energy, 11.);
    assert_eq!(
        p.cv_terms.metadata["collision energy"].as_f64().unwrap(),
        25.
    );
    roundtrip(&p, false);
}

#[test]
fn noncanonical_metadata_remains_user_param_without_type_or_method_changes() {
    let mut p = Precursor::default();
    for (key, value) in [
        ("collision energy", MetaValue::from(25_i64)),
        ("charge stripping", MetaValue::from("false")),
        (
            "supplemental collision-induced dissociation",
            MetaValue::from(""),
        ),
        (
            "supplemental beam-type collision-induced dissociation",
            MetaValue::from(""),
        ),
    ] {
        p.cv_terms.metadata.insert(key.into(), value);
    }
    let xml = roundtrip(&p, false);
    assert!(xml.contains("userParam name=\"collision energy\""));
    assert!(!xml.contains("MS:1002678"));
    assert!(!xml.contains("MS:1002679"));
}

/// Source `MzMLHandler.cpp:4596-4601` writes `MS:1000042` only for a positive
/// intensity, always with its unit attributes. The C++ Release output of
/// MapNormalizer on `inputs/derived/sub_centroid_uk222_picked_first600.mzML`
/// carries the term on the same 149 of 600 spectra as the input, while this
/// port used to add `value="0"` to the 49 remaining precursors.
#[test]
fn absent_precursor_intensity_writes_no_peak_intensity_term() {
    let p = Precursor {
        mz: 684.203_369_140_625,
        charge: 2,
        ..Precursor::default()
    };
    for chrom in [false, true] {
        let xml = roundtrip(&p, chrom);
        assert!(!xml.contains("MS:1000042"), "{xml}");
    }
    let mut measured = p.clone();
    // The C++ output writes 12611.4365234375 for this precursor; that is the
    // f32 this literal denotes, printed shortest.
    measured.intensity = 12_611.437;
    let xml = roundtrip(&measured, false);
    assert!(
        xml.contains(
            "<cvParam cvRef=\"MS\" accession=\"MS:1000042\" name=\"peak intensity\" \
             value=\"12611.437\" unitAccession=\"MS:1000132\" unitCvRef=\"MS\" \
             unitName=\"percent of base peak\"/>"
        ),
        "{xml}"
    );
    // A negative or explicitly united intensity is kept where the source drops
    // it, because omitting it would change the value that reads back.
    let mut negative = p.clone();
    negative.intensity = -1.5;
    assert!(roundtrip(&negative, false).contains("value=\"-1.5\""));
    let mut united = p.clone();
    united
        .cv_terms
        .metadata
        .insert("peak intensity unit accession".into(), "MS:1000131".into());
    let xml = roundtrip(&united, false);
    assert!(xml.contains("accession=\"MS:1000042\""), "{xml}");
    assert!(xml.contains("unitAccession=\"MS:1000131\""), "{xml}");
}

#[test]
fn intensity_default_normalizes_and_nondefault_identity_survives() {
    for unit in [
        "",
        r#"unitAccession="MS:1000132""#,
        r#"unitAccession="MS:1000131" unitCvRef="MS""#,
        r#"unitAccession="UO:0000269" unitCvRef="UO""#,
    ] {
        let p = read(&document(&cv("MS:1000042", "0", unit), "", false), false);
        assert_eq!(
            p.cv_terms
                .metadata
                .contains_key("peak intensity unit accession"),
            unit.contains("1000131") || unit.contains("0000269")
        );
        roundtrip(&p, false);
    }
}

#[test]
fn duplicate_metadata_and_intensity_units_fail_in_both_orders() {
    let user = r#"<userParam name="collision energy" type="xsd:double" value="5"/>"#;
    for terms in [
        [cv("MS:1000045", "25", ""), user.into()],
        [user.into(), cv("MS:1000045", "25", "")],
        [cv("MS:1000045", "25", ""), cv("MS:1000045", "25", "")],
    ] {
        assert!(mzml::read(Cursor::new(document("", &terms.concat(), false))).is_err());
    }
    let user = r#"<userParam name="peak intensity unit accession" value="MS:1000131"/>"#;
    for unit in ["MS:1000131", "MS:1000132"] {
        let term = cv("MS:1000042", "1", &format!(r#"unitAccession="{unit}""#));
        for selected in [
            format!("{user}{term}"),
            format!("{term}{user}"),
            format!("{term}{term}"),
        ] {
            assert!(mzml::read(Cursor::new(document(&selected, "", false))).is_err());
        }
    }
}

#[test]
fn invalid_units_numbers_and_unitless_flag_conflicts_are_rejected() {
    for attrs in [
        r#"unitAccession="MS:1000131" unitCvRef="UO""#,
        r#"unitName="counts""#,
        r#"unitAccession="FAKE:1""#,
        r#"unitAccession="MS:9999999""#,
    ] {
        assert!(
            mzml::read(Cursor::new(document(
                &cv("MS:1000042", "1", attrs),
                "",
                false
            )))
            .is_err()
        );
    }
    for term in [
        cv("MS:1000045", "NaN", ""),
        cv("MS:1002680", "infinity", ""),
        cv("MS:1000245", "", EV),
        cv(
            "MS:1000045",
            "1",
            r#"unitAccession="UO:0000266" unitCvRef="MS""#,
        ),
    ] {
        assert!(mzml::read(Cursor::new(document("", &term, false))).is_err());
    }
}

#[test]
fn invalid_native_unit_metadata_fails_before_writing_bytes() {
    for value in [
        MetaValue::from("MS:1000132"),
        MetaValue::from("MS:9999999"),
        MetaValue::from(1_i64),
        MetaValue::from("MS:1000131")
            .with_unit(Unit::new("UO:0000266", "eV", "UO").unwrap())
            .unwrap(),
    ] {
        let mut p = Precursor::default();
        p.cv_terms
            .metadata
            .insert("peak intensity unit accession".into(), value);
        let mut output = Vec::new();
        assert!(mzml::write(&mut output, &native(p, false)).is_err());
        assert!(output.is_empty());
    }
}

#[test]
fn retained_metadata_obeys_parameter_count_and_byte_budgets() {
    let xml = document(
        &cv("MS:1000042", "30", r#"unitAccession="MS:1000131""#),
        &cv("MS:1000045", "25", EV),
        false,
    );
    for options in [
        ReadOptions {
            max_total_params: 1,
            ..Default::default()
        },
        ReadOptions {
            max_param_bytes: 512,
            ..Default::default()
        },
    ] {
        assert!(mzml::read_with_options(Cursor::new(&xml), &options).is_err());
    }
    assert!(mzml::read(Cursor::new(&xml)).is_ok());
}

#[test]
fn referenceable_activation_and_unit_parameters_share_direct_routes() {
    let groups = format!(
        r#"<referenceableParamGroupList count="2"><referenceableParamGroup id="unit">{}</referenceableParamGroup><referenceableParamGroup id="activation">{}{}</referenceableParamGroup></referenceableParamGroupList>"#,
        cv("MS:1000042", "30", r#"unitAccession="MS:1000131""#),
        cv("MS:1002679", "", ""),
        cv("MS:1002680", "25", EV)
    );
    for chrom in [false, true] {
        let xml = document(
            r#"<referenceableParamGroupRef ref="unit"/>"#,
            r#"<referenceableParamGroupRef ref="activation"/>"#,
            chrom,
        )
        .replace("<run", &format!("{groups}<run"));
        let p = read(&xml, chrom);
        assert_eq!(p.intensity, 30.);
        assert!(p.activation_methods.contains(&A::Etcid));
        assert_eq!(
            p.cv_terms.metadata["supplemental collision energy"]
                .as_f64()
                .unwrap(),
            25.
        );
        roundtrip(&p, chrom);
        let options = ReadOptions {
            max_total_params: 3,
            ..Default::default()
        };
        assert!(mzml::read_with_options(Cursor::new(xml), &options).is_err());
    }
}

/// `activation` is a ParamGroupType, whose schema sequence puts every cvParam
/// before any userParam. A precursor with promoted metadata but no method and no
/// energy used to get its "activation information unavailable" userParam ahead of
/// the promoted cvParams, which the mzML 1.1.0 XSD rejects.
#[test]
fn metadata_only_activation_puts_every_cv_param_before_user_params() {
    use std::io::Write;
    use std::process::{Command, Stdio};
    let xmllint = Command::new("xmllint").arg("--version").output().is_ok();
    for chrom in [false, true] {
        let p = read(&document("", &cv("MS:1000045", "25", EV), chrom), chrom);
        assert!(p.activation_methods.is_empty());
        let xml = roundtrip(&p, chrom);
        let start = xml.find("<activation>").unwrap();
        let end = xml.find("</activation>").unwrap();
        let block = &xml[start..end];
        let last_cv = block.rfind("<cvParam").unwrap();
        let first_user = block.find("<userParam").unwrap();
        assert!(last_cv < first_user, "{block}");
        if !xmllint {
            eprintln!("xmllint unavailable; activation XSD validation not executed");
            continue;
        }
        // `mzml::write`, which `roundtrip` uses, is indexed by default.
        let schema = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(
            if xml.contains("<indexedmzML ") {
                "tests/data/mzml_writing/mzML_idx_1_10.xsd"
            } else {
                "tests/data/mzml_1_10.xsd"
            },
        );
        let mut process = Command::new("xmllint")
            .args(["--nonet", "--noout", "--schema"])
            .arg(&schema)
            .arg("-")
            .stdin(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        process
            .stdin
            .take()
            .unwrap()
            .write_all(xml.as_bytes())
            .unwrap();
        let result = process.wait_with_output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
}
