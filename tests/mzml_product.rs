// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "mzml")]

use openms::format::mzml::{self, ReadOptions, WriteOptions};
use openms::kernel::{ChromatogramPeak, DataArray, MSChromatogram, MSExperiment, Precursor};
use openms::metadata::{CVTerm, MetaValue, MetaValueData, Product, Unit};
use std::io::{Cursor, Write};

fn doc(groups: &str, contents: &str) -> String {
    format!(
        "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\">{groups}<run id=\"r\"><chromatogramList count=\"1\"><chromatogram id=\"c\" defaultArrayLength=\"0\">{contents}</chromatogram></chromatogramList></run></mzML>"
    )
}
fn read(text: &str) -> openms::Result<MSExperiment> {
    mzml::read(Cursor::new(text.as_bytes()))
}
fn cv(accession: &str, value: &str) -> String {
    format!("<cvParam cvRef=\"MS\" accession=\"{accession}\" name=\"quantity\" value=\"{value}\"/>")
}
fn product(contents: &str) -> String {
    format!("<product><isolationWindow>{contents}</isolationWindow></product>")
}
fn populated() -> MSExperiment {
    let mut c = MSChromatogram::from_peaks(vec![
        ChromatogramPeak::new(1.0, 3.0),
        ChromatogramPeak::new(2.0, -4.0),
    ]);
    c.native_id = "transition_1".into();
    c.precursor = Precursor::new(600.0, 2);
    c.product = Product {
        mz: 250.125,
        isolation_window_lower_offset: 0.4,
        isolation_window_upper_offset: 0.75,
        ..Default::default()
    };
    c.product
        .cv_terms
        .metadata
        .insert("label<&\"\t\n".into(), "ion<&\"\t\n".into());
    c.product
        .cv_terms
        .metadata
        .insert("integer".into(), MetaValue::from(i64::MIN));
    c.product.cv_terms.metadata.insert(
        "double".into(),
        MetaValue::new(MetaValueData::Float(-0.25))
            .unwrap()
            .with_unit(Unit::new("UO:0000010", "second", "UO").unwrap())
            .unwrap(),
    );
    c.product
        .cv_terms
        .metadata
        .insert("empty string".into(), "".into());
    c.integer_data_arrays
        .push(DataArray::new("indices", vec![7, 8]));
    MSExperiment {
        chromatograms: vec![c],
        ..Default::default()
    }
}

#[test]
fn source_product_literals_and_metadata_use_chromatogram_owner() {
    // MzMLFile_1.mzML lines418-425 and class assertions535-537,653.
    // The same source handler routes this spectrum Product to a chromatogram
    // when outside spectrumList; this projection exercises that owner.
    let params = format!(
        "{}{}{}<userParam name=\"iwname\" value=\"isolationwindow3\"/>",
        cv("MS:1000827", "18.88"),
        cv("MS:1000828", "1.0"),
        cv("MS:1000829", "2.0")
    );
    let precursor = "<precursor><selectedIonList count=\"1\"><selectedIon><cvParam accession=\"MS:1000744\" value=\"99\"/></selectedIon></selectedIonList><activation/></precursor>";
    let e = read(&doc(
        "",
        &format!(
            "{precursor}{}<userParam name=\"owner\" value=\"chromatogram\"/>",
            product(&params)
        ),
    ))
    .unwrap();
    let c = &e.chromatograms[0];
    assert_eq!(
        (
            c.product.mz,
            c.product.isolation_window_lower_offset,
            c.product.isolation_window_upper_offset
        ),
        (18.88, 1.0, 2.0)
    );
    assert_eq!(
        c.product.cv_terms.metadata["iwname"].as_str().unwrap(),
        "isolationwindow3"
    );
    assert_eq!(c.precursor.mz, 99.0);
    assert_eq!(c.metadata["owner"], "chromatogram");
    assert!(!c.metadata.contains_key("iwname"));
}

#[test]
fn product_and_annotation_roundtrip_in_both_compression_modes() {
    let e = populated();
    let original = e.clone();
    for compressed in [false, true] {
        let mut bytes = Vec::new();
        mzml::write_with_options(
            &mut bytes,
            &e,
            &WriteOptions {
                zlib_compression: compressed,
            },
        )
        .unwrap();
        assert_eq!(mzml::read(Cursor::new(bytes)).unwrap(), e);
    }
    assert_eq!(e, original);
}

#[test]
fn grouped_product_params_share_field_checks_and_cumulative_limits() {
    let params = format!(
        "{}{}<userParam name=\"number\" type=\"xsd:int\" value=\"17\"/>",
        cv("MS:1000827", "500.5"),
        cv("MS:1000829", "2")
    );
    let groups = format!(
        "<referenceableParamGroupList count=\"1\"><referenceableParamGroup id=\"p\">{params}</referenceableParamGroup></referenceableParamGroupList>"
    );
    let referenced = doc(&groups, &product("<referenceableParamGroupRef ref=\"p\"/>"));
    assert_eq!(
        read(&referenced).unwrap(),
        read(&doc("", &product(&params))).unwrap()
    );
    let repeated = doc(
        &groups,
        &product("<referenceableParamGroupRef ref=\"p\"/><referenceableParamGroupRef ref=\"p\"/>"),
    );
    assert!(read(&repeated).is_err());
    let duplicate_inline = doc(
        &groups,
        &product(&format!(
            "<referenceableParamGroupRef ref=\"p\"/>{}",
            cv("MS:1000827", "500.5")
        )),
    );
    assert!(read(&duplicate_inline).is_err());
    let options = ReadOptions {
        max_total_params: 5,
        ..Default::default()
    };
    assert!(mzml::read_with_options(Cursor::new(&referenced), &options).is_err());
    let options = ReadOptions {
        max_param_bytes: 4096,
        ..Default::default()
    };
    assert!(mzml::read_with_options(Cursor::new(&referenced), &options).is_err());
}

#[test]
fn source_scalar_xsd_dispatch_preserves_types_ranges_and_unit_identity() {
    let contents = "<userParam name=\"a\" type=\"xsd:unsignedByte\" value=\"-17\"/><userParam name=\"b\" type=\"xsd:unsignedLong\" value=\"+-9223372036854775808\"/><userParam name=\"c\" type=\"xsd:float\" value=\"1.0000000000000002\" unitAccession=\"MS:1000040\" unitCvRef=\"MS\" unitName=\"m/z\"/><userParam name=\"d\" type=\"xsd:boolean\" value=\"true\"/>";
    let e = read(&doc("", &product(contents))).unwrap();
    let meta = &e.chromatograms[0].product.cv_terms.metadata;
    assert_eq!(meta["a"].as_i64().unwrap(), -17); // Source does not enforce nominal byte unsignedness.
    assert_eq!(meta["b"].as_i64().unwrap(), i64::MIN);
    assert_eq!(
        meta["c"].as_f64().unwrap(),
        f64::from_bits(1.0f64.to_bits() + 1)
    );
    assert_eq!(
        meta["c"].unit().unwrap(),
        &Unit::new("MS:1000040", "m/z", "MS").unwrap()
    );
    assert_eq!(meta["d"].as_str().unwrap(), "true");
    for bad in [
        "<userParam name=\"n\" type=\"xsd:unsignedInt\" value=\"2147483648\"/>",
        "<userParam name=\"n\" type=\"xsd:integer\" value=\"9223372036854775808\"/>",
        "<userParam name=\"n\" type=\"xsd:double\" value=\"NaN\"/>",
        "<userParam name=\"n\" value=\"s\" unitName=\"second\"/>",
        "<userParam name=\"n\" value=\"s\" unitAccession=\"UO:0000010\" unitCvRef=\"MS\"/>",
        "<userParam name=\"n\" value=\"a\"/><userParam name=\"n\" value=\"b\"/>",
    ] {
        assert!(read(&doc("", &product(bad))).is_err(), "{bad}");
    }
}

#[test]
fn product_structure_conflicts_and_spectrum_product_support() {
    let valid = product(&cv("MS:1000827", "1"));
    for bad in [
        format!("{valid}{valid}"),
        "<product><isolationWindow/><isolationWindow/></product>".into(),
        "<product><product/></product>".into(),
        "<product><precursor/></product>".into(),
        "<precursor><product/></precursor>".into(),
        "<product><activation/></product>".into(),
        "<product><userParam name=\"misplaced\"/></product>".into(),
        "<isolationWindow/>".into(),
        "<productList count=\"0\"/>".into(),
    ] {
        assert!(read(&doc("", &bad)).is_err(), "{bad}");
    }
    let spectrum = "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\"><run><spectrumList count=\"1\"><spectrum id=\"scan=1\" defaultArrayLength=\"0\"><productList count=\"1\"><product/></productList></spectrum></spectrumList></run></mzML>";
    assert_eq!(
        read(spectrum).unwrap().spectra[0].products,
        [Product::default()]
    );
}

#[test]
fn signed_target_zero_defaults_and_checked_offsets() {
    assert_eq!(
        read(&doc("", "<product/>")).unwrap().chromatograms[0].product,
        Product::default()
    );
    let e = read(&doc("", &product(&cv("MS:1000827", "-5")))).unwrap();
    assert_eq!(e.chromatograms[0].product.mz, -5.0);
    let mut bytes = Vec::new();
    mzml::write(&mut bytes, &e).unwrap();
    assert_eq!(mzml::read(Cursor::new(bytes)).unwrap(), e);
    for params in [
        cv("MS:1000827", "NaN"),
        cv("MS:1000828", "-1"),
        cv("MS:1000829", "inf"),
        cv("MS:1000000", "1"),
        format!("{}{}", cv("MS:1000828", "1"), cv("MS:1000828", "1")),
        "<cvParam accession=\"MS:1000827\" value=\"1\" unitAccession=\"UO:0000010\"/>".into(),
    ] {
        assert!(read(&doc("", &product(&params))).is_err(), "{params}");
    }
}

#[derive(Default)]
struct CountWrites(usize);
impl Write for CountWrites {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.0 += 1;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
#[test]
fn every_unrepresentable_product_fails_whole_document_preflight() {
    type Mutation = fn(&mut Product);
    let mutations: &[Mutation] = &[
        |p| p.mz = f64::NAN,
        |p| p.isolation_window_upper_offset = -1.0,
        |p| {
            p.cv_terms
                .add(CVTerm::new("MS:1000000", "custom", "MS"))
                .unwrap();
        },
        |p| {
            p.cv_terms.replace_accession("MS:1000000", vec![]).unwrap();
        },
        |p| {
            p.cv_terms
                .metadata
                .insert("empty".into(), MetaValue::default());
        },
        |p| {
            p.cv_terms.metadata.insert(
                "list".into(),
                MetaValue::new(MetaValueData::StringList(vec![])).unwrap(),
            );
        },
        |p| {
            p.cv_terms.metadata.insert(
                "list".into(),
                MetaValue::new(MetaValueData::IntegerList(vec![1])).unwrap(),
            );
        },
        |p| {
            p.cv_terms.metadata.insert(
                "list".into(),
                MetaValue::new(MetaValueData::FloatList(vec![1.0])).unwrap(),
            );
        },
        |p| {
            p.cv_terms.metadata.insert("bad\0".into(), "text".into());
        },
        |p| {
            p.cv_terms.metadata.insert("bad".into(), "text\0".into());
        },
        |p| {
            p.cv_terms.metadata.insert(
                "unit".into(),
                MetaValue::from(1i64)
                    .with_unit(Unit::new("XX:1", "unknown", "XX").unwrap())
                    .unwrap(),
            );
        },
    ];
    for mutate in mutations {
        let mut e = populated();
        e.chromatograms.push(e.chromatograms[0].clone());
        e.chromatograms[1].native_id = "last".into();
        mutate(&mut e.chromatograms[1].product);
        let mut writer = CountWrites::default();
        assert!(mzml::write(&mut writer, &e).is_err());
        assert_eq!(writer.0, 0);
    }
}

#[test]
fn populated_product_output_validates_against_independent_schema_when_available() {
    use std::process::{Command, Stdio};
    if Command::new("xmllint").arg("--version").output().is_err() {
        eprintln!("xmllint unavailable; Product XSD validation not executed");
        return;
    }
    let schema = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/mzml_1_10.xsd");
    for compressed in [false, true] {
        let mut bytes = Vec::new();
        mzml::write_with_options(
            &mut bytes,
            &populated(),
            &WriteOptions {
                zlib_compression: compressed,
            },
        )
        .unwrap();
        let mut process = Command::new("xmllint")
            .args(["--nonet", "--noout", "--schema"])
            .arg(&schema)
            .arg("-")
            .stdin(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        process.stdin.take().unwrap().write_all(&bytes).unwrap();
        let result = process.wait_with_output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
}

#[test]
fn extracted_ion_chromatograms_keep_midpoints_and_samples_through_mzml() -> openms::Result<()> {
    use openms::kernel::{MSSpectrum, MzRtRegion, Peak1D};
    let input = MSExperiment {
        spectra: vec![
            MSSpectrum {
                rt: 1.0,
                peaks: vec![
                    Peak1D::new(100.0, 10.0),
                    Peak1D::new(101.0, 20.0),
                    Peak1D::new(200.0, 30.0),
                ],
                ..Default::default()
            },
            MSSpectrum {
                rt: 2.0,
                ms_level: 2,
                peaks: vec![Peak1D::new(100.0, 999.0)],
                ..Default::default()
            },
            MSSpectrum {
                rt: 3.0,
                peaks: vec![Peak1D::new(100.0, 11.0), Peak1D::new(200.0, 31.0)],
                ..Default::default()
            },
        ],
        ..Default::default()
    };
    let saved = input.clone();
    let regions = [
        MzRtRegion::new(99.5, 101.5, 0.0, 4.0)?,
        MzRtRegion::new(199.0, 201.0, 0.0, 4.0)?,
    ];
    let chromatograms = input.extract_xics(&regions, 1)?;
    assert_eq!(
        (chromatograms[0].product.mz, chromatograms[1].product.mz),
        (100.5, 200.0)
    );
    assert_eq!(
        chromatograms[0].peaks,
        [
            ChromatogramPeak::new(1.0, 30.0),
            ChromatogramPeak::new(3.0, 11.0)
        ]
    );
    let mut expected = MSExperiment {
        chromatograms,
        ..Default::default()
    };
    // The generic writer assigns IDs to the source's unnamed generated traces.
    for (index, c) in expected.chromatograms.iter_mut().enumerate() {
        c.native_id = format!("chromatogram={index}");
    }
    for compressed in [false, true] {
        let mut output = Vec::new();
        mzml::write_with_options(
            &mut output,
            &expected,
            &WriteOptions {
                zlib_compression: compressed,
            },
        )?;
        assert_eq!(mzml::read(Cursor::new(output))?, expected);
    }
    assert_eq!(input, saved);
    Ok(())
}
