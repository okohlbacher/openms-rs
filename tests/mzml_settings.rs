// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "mzml")]

use openms::format::mzml::{self, ReadOptions};
use openms::kernel::{
    ChromatogramPeak, ChromatogramTools, MSChromatogram, MSExperiment, MSSpectrum, Peak1D,
};
use openms::metadata::{
    ChromatogramType, MetaValue, Polarity, Product, ScanMode, ScanWindow, Unit,
};
use std::io::Cursor;

fn read(xml: &str) -> openms::Result<MSExperiment> {
    mzml::read(Cursor::new(xml.as_bytes()))
}
fn spectrum(body: &str) -> String {
    format!(
        "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\"><run><spectrumList count=\"1\"><spectrum id=\"scan=1\" defaultArrayLength=\"0\">{body}</spectrum></spectrumList></run></mzML>"
    )
}
fn cv(accession: &str, value: &str) -> String {
    format!("<cvParam accession=\"{accession}\" value=\"{value}\"/>")
}
fn roundtrip(e: &MSExperiment, compressed: bool) -> MSExperiment {
    let mut bytes = Vec::new();
    if compressed {
        let config = openms::format::numpress_coder::NumpressConfig {
            compression: openms::format::numpress_coder::NumpressCompression::Pic,
            ..Default::default()
        };
        mzml::write_with_numpress(
            &mut bytes,
            e,
            &mzml::NumpressWriteOptions {
                binary: mzml::WriteOptions {
                    zlib_compression: true,
                },
                mass_time: config,
                intensity: config,
                ..Default::default()
            },
        )
        .unwrap();
    } else {
        mzml::write(&mut bytes, e).unwrap();
    }
    mzml::read(Cursor::new(bytes)).unwrap()
}

#[test]
fn source_class_literals_windows_modes_products_and_user_metadata() {
    // MzMLFile_test.cpp 364–568,653–656. Projection preserves literal settings
    // and groups, strips unrelated unsupported acquisition/binary data only.
    let e = read(include_str!("data/mzml_settings_source_projection.mzML")).unwrap();
    assert_eq!(e.spectra.len(), 4);
    let expected: &[&[(f64, f64)]] = &[
        &[(400., 1800.)],
        &[(100., 500.), (600., 1000.), (1100., 1500.)],
        &[(400., 1800.)],
        &[(110., 905.)],
    ];
    for (i, s) in e.spectra.iter().enumerate() {
        assert_eq!(
            s.instrument_settings.scan_mode,
            if i == 1 {
                ScanMode::MsnSpectrum
            } else {
                ScanMode::Ms1Spectrum
            }
        );
        assert_eq!(s.instrument_settings.polarity, Polarity::Positive);
        assert_eq!(s.instrument_settings.zoom_scan, i == 3);
        let bounds: Vec<_> = s
            .instrument_settings
            .scan_windows
            .iter()
            .map(|w| (w.begin, w.end))
            .collect();
        assert_eq!(bounds, expected[i]);
    }
    assert_eq!(
        e.spectra[0].instrument_settings.scan_windows[0].metadata["name"]
            .as_str()
            .unwrap(),
        "scanwindow1"
    );
    let products = &e.spectra[2].products;
    assert_eq!(products.len(), 2);
    for (p, (mz, lo, hi, label)) in products.iter().zip([
        (18.88, 1., 2., "isolationwindow3"),
        (19.99, 3., 4., "isolationwindow4"),
    ]) {
        assert_eq!(
            (
                p.mz,
                p.isolation_window_lower_offset,
                p.isolation_window_upper_offset
            ),
            (mz, lo, hi)
        );
        assert_eq!(p.cv_terms.metadata["iwname"].as_str().unwrap(), label);
    }
    // This older projection deliberately removed combination CVs but kept two
    // scans in spectra 0/2. The acquisition writer's documented source fallback
    // supplies no-combination; all source settings/quantities remain unchanged.
    let mut expected = e.clone();
    for index in [0, 2] {
        expected.spectra[index]
            .acquisition_info
            .method_of_combination = "no combination".into();
    }
    for compressed in [false, true] {
        assert_eq!(roundtrip(&e, compressed), expected);
    }
}

#[test]
fn every_source_scan_mode_polarity_and_chromatogram_type_roundtrips() {
    let scan_terms = [
        "", "1000294", "1000579", "1000580", "1000582", "1000583", "1000581", "1000325", "1000326",
        "1000341", "1000789", "1000790", "1000804", "1000805", "1000806",
    ];
    for (&mode, term) in ScanMode::ALL.iter().zip(scan_terms) {
        for &polarity in Polarity::ALL {
            let mut e = MSExperiment::default();
            let mut s = MSSpectrum {
                native_id: "scan=1".into(),
                ..Default::default()
            };
            s.instrument_settings.scan_mode = mode;
            s.instrument_settings.polarity = polarity;
            e.spectra.push(s);
            let mut bytes = Vec::new();
            mzml::write(&mut bytes, &e).unwrap();
            let xml = String::from_utf8_lossy(&bytes);
            let content = xml
                .split("<fileContent>")
                .nth(1)
                .unwrap()
                .split("</fileContent>")
                .next()
                .unwrap();
            if !term.is_empty() {
                assert!(content.contains(&format!("MS:{term}")));
                if mode != ScanMode::MassSpectrum {
                    assert!(!content.contains("MS:1000294"));
                }
            } else {
                assert!(content.contains("MS:1000294"));
            }
            assert_eq!(mzml::read(Cursor::new(bytes)).unwrap(), e);
        }
    }
    let terms = [
        "1000810", "1000235", "1000627", "1000628", "1001472", "1001473", "1000811", "1000812",
        "1000813",
    ];
    for (&kind, term) in ChromatogramType::ALL.iter().zip(terms) {
        let c = MSChromatogram {
            native_id: "c".into(),
            chromatogram_type: kind,
            ..Default::default()
        };
        let e = MSExperiment {
            chromatograms: vec![c],
            ..Default::default()
        };
        let mut bytes = Vec::new();
        mzml::write(&mut bytes, &e).unwrap();
        let xml = String::from_utf8(bytes).unwrap();
        assert!(xml.contains(&format!("MS:{term}")));
        assert_eq!(read(&xml).unwrap(), e);
        if kind == ChromatogramType::SelectedReactionMonitoring {
            assert_eq!(read(&xml.replace("MS:1001473", "MS:1001474")).unwrap(), e);
        }
    }
}

#[test]
fn scan_windows_units_scalar_metadata_and_missing_rt_preserve_exact_values() {
    for unit in ["UO:0000018", "UO:0000028", "MS:1002814"] {
        let mut w = ScanWindow {
            begin: 220.,
            end: 500.,
            ..Default::default()
        };
        w.metadata.insert("unit_accession".into(), unit.into());
        w.metadata.insert("note".into(), "<&\"\n".into());
        w.metadata.insert("count".into(), 42_i64.into());
        w.metadata.insert(
            "energy".into(),
            MetaValue::new(openms::metadata::MetaValueData::Float(3.5))
                .unwrap()
                .with_unit(Unit::new("UO:0000266", "electronvolt", "UO").unwrap())
                .unwrap(),
        );
        let mut s = MSSpectrum {
            native_id: "scan=1".into(),
            ..Default::default()
        };
        s.instrument_settings.scan_windows.push(w);
        s.instrument_settings.zoom_scan = true;
        let e = MSExperiment {
            spectra: vec![s],
            ..Default::default()
        };
        assert_eq!(roundtrip(&e, false), e);
        assert_eq!(roundtrip(&e, true), e);
        assert_eq!(e.spectra[0].rt, -1.);
    }
    let legacy = spectrum(&cv("MS:1000497", ""));
    assert!(
        read(&legacy).unwrap().spectra[0]
            .instrument_settings
            .zoom_scan
    );
    let empty = spectrum(
        "<scanList count=\"1\"><scan><scanWindowList count=\"1\"><scanWindow/></scanWindowList></scan></scanList>",
    );
    assert_eq!(
        read(&empty).unwrap().spectra[0]
            .instrument_settings
            .scan_windows,
        [ScanWindow::default()]
    );
}

#[test]
fn structural_duplicates_conflicts_and_nonfinite_fields_are_checked() {
    let window = "<scanList count=\"1\"><scan><scanWindowList count=\"1\"><scanWindow><cvParam accession=\"MS:1000501\" value=\"1\"/><cvParam accession=\"MS:1000500\" value=\"2\"/></scanWindow></scanWindowList></scan></scanList>";
    for contents in [
        "<scanWindow/>".into(), "<scan/>".into(),
        "<scanList count=\"1\"><scanWindowList count=\"0\"/></scanList>".into(),
        window.replace("count=\"1\"><scanWindow>","count=\"2\"><scanWindow>"),
        window.replace("value=\"2\"","value=\"NaN\""),
        window.replace("value=\"2\"","value=\"0\""),
        window.replace("<cvParam accession=\"MS:1000500\"", "<cvParam accession=\"MS:1000501\""),
        window.replace("value=\"2\"", "value=\"2\" unitAccession=\"UO:0000018\""),
        window.replace("</scanWindow>","<userParam name=\"unit_accession\" value=\"UO:0000018\"/></scanWindow>"),
        format!("{}{}",cv("MS:1000579",""),cv("MS:1000580","")),
        format!("{}{}",cv("MS:1000130",""),cv("MS:1000129","")),
        "<productList count=\"2\"><product/></productList>".into(),
        "<productList count=\"1\"><product><isolationWindow/><isolationWindow/></product></productList>".into(),
    ] { assert!(read(&spectrum(&contents)).is_err(),"{contents}"); }
}

#[test]
fn settings_structures_and_references_share_parameter_limits_even_when_filtered() {
    let window = "<scanList count=\"1\"><scan><scanWindowList count=\"2\"><scanWindow/><scanWindow/></scanWindowList></scan></scanList>";
    let xml = spectrum(window);
    assert!(
        mzml::read_with_options(
            Cursor::new(&xml),
            &ReadOptions {
                max_total_params: 1,
                ..Default::default()
            }
        )
        .is_err()
    );
    assert!(
        mzml::read_with_options(
            Cursor::new(&xml),
            &ReadOptions {
                max_param_bytes: 1024,
                ..Default::default()
            }
        )
        .is_err()
    );
    let products = spectrum("<productList count=\"2\"><product/><product/></productList>");
    assert!(
        mzml::read_with_options(
            Cursor::new(&products),
            &ReadOptions {
                max_total_params: 1,
                ..Default::default()
            }
        )
        .is_err()
    );
    let bad = spectrum(
        "<scanList count=\"1\"><scan><scanWindowList count=\"1\"><scanWindow><cvParam accession=\"MS:1000501\" value=\"NaN\"/></scanWindow></scanWindowList></scan></scanList>",
    );
    let options = mzml::LoadOptions {
        skip_spectra: true,
        ..Default::default()
    };
    assert!(mzml::read_with_load_options(Cursor::new(bad), &options, &Default::default()).is_err());
    let groups = "<referenceableParamGroupList count=\"1\"><referenceableParamGroup id=\"g\"><cvParam accession=\"MS:1000501\" value=\"1\"/><cvParam accession=\"MS:1000500\" value=\"2\"/></referenceableParamGroup></referenceableParamGroupList>";
    let xml = spectrum("<scanList count=\"1\"><scan><scanWindowList count=\"1\"><scanWindow><referenceableParamGroupRef ref=\"g\"/></scanWindow></scanWindowList></scan></scanList>").replace("<run>",&format!("{groups}<run>"));
    let e = read(&xml).unwrap();
    assert_eq!(e.spectra[0].instrument_settings.scan_windows[0].end, 2.);
    assert!(
        mzml::read_with_options(
            Cursor::new(xml),
            &ReadOptions {
                max_total_params: 5,
                ..Default::default()
            }
        )
        .is_err()
    );
}

#[test]
fn chromatogram_conversion_survives_standard_spectrum_settings_transport() {
    let mut c = MSChromatogram::from_peaks(vec![
        ChromatogramPeak::new(2., 10.),
        ChromatogramPeak::new(5., 20.),
    ]);
    c.native_id = "transition".into();
    c.chromatogram_type = ChromatogramType::SelectedReactionMonitoring;
    c.precursor.mz = 500.;
    c.product.mz = 250.;
    let mut e = MSExperiment {
        chromatograms: vec![c],
        ..Default::default()
    };
    let tools = ChromatogramTools::default();
    tools.convert_chromatograms_to_spectra(&mut e).unwrap();
    // The converter leaves native IDs empty; the pre-existing writer fills
    // those. Give the scans stable IDs before checking complete value equality.
    for (i, s) in e.spectra.iter_mut().enumerate() {
        s.native_id = format!("index={i}");
    }
    assert!(
        e.spectra
            .iter()
            .all(|s| s.instrument_settings.scan_mode == ScanMode::SelectedReactionMonitoring)
    );
    for compressed in [false, true] {
        let mut loaded = roundtrip(&e, compressed);
        assert_eq!(loaded, e);
        tools
            .convert_spectra_to_chromatograms(&mut loaded, true, false)
            .unwrap();
        assert_eq!(loaded.chromatograms.len(), 1);
        let c = &loaded.chromatograms[0];
        assert_eq!(
            c.chromatogram_type,
            ChromatogramType::SelectedReactionMonitoring
        );
        assert_eq!(c.product.mz, 250.);
        assert_eq!(c.precursor.mz, 500.);
        assert_eq!(
            c.peaks,
            [
                ChromatogramPeak::new(2., 10.),
                ChromatogramPeak::new(5., 20.)
            ]
        );
        // The source mzML chromatogram grammar has no InstrumentSettings slot.
        let mut bytes = b"preserved".to_vec();
        assert!(mzml::write(&mut bytes, &loaded).is_err());
        assert_eq!(bytes, b"preserved");
    }
}

#[test]
fn product_defaults_signed_zero_and_precursor_scope_do_not_alias() {
    let mut s = MSSpectrum::from_peaks(vec![Peak1D::new(1., 2.)]);
    s.native_id = "scan=1".into();
    s.precursors.push(openms::Precursor::new(500., 2));
    s.products = vec![
        Product::default(),
        Product {
            mz: -0.0,
            isolation_window_lower_offset: -0.0,
            ..Default::default()
        },
    ];
    let e = MSExperiment {
        spectra: vec![s],
        ..Default::default()
    };
    let output = roundtrip(&e, false);
    assert_eq!(output, e);
    assert_eq!(
        output.spectra[0].products[1].mz.to_bits(),
        (-0.0_f64).to_bits()
    );
    assert_eq!(
        output.spectra[0].products[1]
            .isolation_window_lower_offset
            .to_bits(),
        (-0.0_f64).to_bits()
    );
}

#[test]
fn settings_and_product_output_pass_independent_schema_in_both_writers() {
    use std::io::Write;
    use std::process::{Command, Stdio};
    if Command::new("xmllint").arg("--version").output().is_err() {
        eprintln!("xmllint unavailable; mzML settings XSD validation not executed");
        return;
    }
    let mut e = read(include_str!("data/mzml_settings_source_projection.mzML")).unwrap();
    e.spectra[0].instrument_settings.scan_windows[0]
        .metadata
        .insert("unit_accession".into(), "UO:0000018".into());
    // Test the standard optional-unitName spelling for an opaque numeric unit.
    e.spectra[1].instrument_settings.scan_windows[0]
        .metadata
        .insert("unit_accession".into(), "UO:0000028".into());
    for (i, &kind) in ChromatogramType::ALL
        .iter()
        .filter(|&&k| k != ChromatogramType::Unknown)
        .enumerate()
    {
        e.chromatograms.push(MSChromatogram {
            native_id: format!("c{i}"),
            chromatogram_type: kind,
            ..Default::default()
        });
    }
    let schema = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/mzml_1_10.xsd");
    for compressed in [false, true] {
        let mut bytes = Vec::new();
        if compressed {
            mzml::write_with_numpress(&mut bytes, &e, &Default::default()).unwrap();
        } else {
            mzml::write(&mut bytes, &e).unwrap();
        }
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
