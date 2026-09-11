#![cfg(feature = "mzml")]
use openms::format::{mzml, peak_options::PeakFileOptions};
use openms::kernel::{DataArray, MSChromatogram, MSExperiment, MSSpectrum, NumericRange, Peak1D};
use openms::metadata::{ChromatogramType, MetaValue, MetaValueData, ScanMode, Unit};
use std::io::Cursor;
const NOISE: [&str; 3] = [
    "sampled noise m/z array",
    "sampled noise intensity array",
    "sampled noise baseline array",
];
fn spectrum() -> MSSpectrum {
    MSSpectrum {
        native_id: "scan=1".into(),
        peaks: vec![Peak1D::new(501., 10.), Peak1D::new(499., 20.)],
        ..Default::default()
    }
}
fn experiment() -> MSExperiment {
    MSExperiment {
        spectra: vec![spectrum()],
        ..Default::default()
    }
}
fn write(e: &MSExperiment) -> String {
    let mut b = Vec::new();
    mzml::write(&mut b, e).unwrap();
    String::from_utf8(b).unwrap()
}
fn read(xml: &str) -> MSExperiment {
    mzml::read_with_load_options(
        Cursor::new(xml),
        &mzml::LoadOptions::default(),
        &mzml::ReadOptions::default(),
    )
    .unwrap()
}
fn list(values: &[f64]) -> MetaValue {
    MetaValue::try_from(values.to_vec()).unwrap()
}
fn array(values: &[f64], term: &str, unit: &str, params: &str) -> String {
    use base64::Engine;
    let bytes: Vec<_> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
    let encoded = base64::engine::general_purpose::STANDARD.encode(bytes);
    format!(
        "<binaryDataArray arrayLength=\"{}\" encodedLength=\"{}\"><cvParam accession=\"{term}\" name=\"array\" {unit}/><cvParam accession=\"MS:1000523\" name=\"64-bit float\"/><cvParam accession=\"MS:1000576\" name=\"no compression\"/>{params}<binary>{encoded}</binary></binaryDataArray>",
        values.len(),
        encoded.len()
    )
}
fn arrays(xml: &str, values: &[String]) -> String {
    let a = xml.find("<binaryDataArrayList").unwrap();
    let b = xml[a..].find("</binaryDataArrayList>").unwrap() + a + "</binaryDataArrayList>".len();
    format!(
        "{}<binaryDataArrayList count=\"{}\">{}</binaryDataArrayList>{}",
        &xml[..a],
        values.len(),
        values.concat(),
        &xml[b..]
    )
}
fn write_with_scan() -> String {
    let mut e = experiment();
    e.spectra[0].rt = 1.;
    write(&e)
}
fn user(name: &str, kind: &str, value: &str) -> String {
    format!("<userParam name=\"{name}\" type=\"{kind}\" value=\"{value}\"/>")
}
#[test]
fn source_noise_literals_are_independent_of_sorted_charge_and_filtered_peaks() {
    let mut e = experiment();
    let noise = [498.123456789012, 500.345678901234, 502.567890123456];
    e.spectra[0].metadata.insert(NOISE[0].into(), list(&noise));
    e.spectra[0]
        .metadata
        .insert(NOISE[1].into(), list(&[1., 2., 3.]));
    e.spectra[0]
        .metadata
        .insert(NOISE[2].into(), list(&[0.1, 0.2, 0.3]));
    e.spectra[0]
        .integer_data_arrays
        .push(DataArray::new("charge array", vec![1, 2]));
    let xml = write(&e);
    let loaded = read(&xml);
    assert_eq!(
        loaded.spectra[0].metadata[NOISE[0]]
            .as_float_list()
            .unwrap(),
        noise
    );
    assert!(loaded.spectra[0].float_data_arrays.is_empty());
    assert_eq!(loaded.spectra[0].integer_data_arrays[0].data, [2, 1]);
    let mut o = mzml::LoadOptions::default();
    o.scientific.set_mz_range(NumericRange {
        min: 500.,
        max: 501.,
    });
    let filtered =
        mzml::read_with_load_options(Cursor::new(&xml), &o, &Default::default()).unwrap();
    assert!(filtered.spectra[0].is_empty());
    assert_eq!(filtered.spectra[0].metadata, loaded.spectra[0].metadata);
    for (key, n) in NOISE.iter().zip([3, 3, 3]) {
        assert_eq!(
            loaded.spectra[0].metadata[*key]
                .as_float_list()
                .unwrap()
                .len(),
            n
        );
    }
}
#[test]
fn source_pda_and_pressure_roles_roundtrip_with_units_and_zero_level() {
    let mut e = experiment();
    let s = &mut e.spectra[0];
    s.ms_level = 0;
    s.instrument_settings.scan_mode = ScanMode::Absorption;
    s.metadata
        .insert("mzml coordinate array".into(), "wavelength".into());
    s.metadata
        .insert("mzml intensity array".into(), "absorption".into());
    s.peaks = vec![Peak1D::new(220., 0.25), Peak1D::new(500., 0.5)];
    let mut c = MSChromatogram {
        native_id: "Analog#1_Pressure".into(),
        peaks: vec![
            openms::kernel::ChromatogramPeak::new(1., 200.),
            openms::kernel::ChromatogramPeak::new(2., 210.),
        ],
        ..Default::default()
    };
    c.metadata
        .insert("mzml intensity array".into(), "pressure".into());
    c.metadata
        .insert("chromatogram type accession".into(), "MS:1003019".into());
    c.metadata
        .insert("Thermo detector units".into(), "bar".into());
    e.chromatograms.push(c);
    let xml = write(&e);
    for id in ["MS:1000617", "UO:0000269", "MS:1003019", "MS:1000821"] {
        assert!(xml.contains(id));
    }
    let loaded = read(&xml);
    assert_eq!(loaded.spectra[0].ms_level, 0);
    assert_eq!(loaded.spectra[0].peaks, e.spectra[0].peaks);
    assert_eq!(
        loaded.spectra[0].metadata["mzml intensity array"]
            .as_str()
            .unwrap(),
        "absorption"
    );
    assert_eq!(
        loaded.chromatograms[0].metadata["mzml intensity array"]
            .as_str()
            .unwrap(),
        "pressure"
    );
    assert_eq!(
        loaded.chromatograms[0].metadata["Thermo detector units"]
            .as_str()
            .unwrap(),
        "bar"
    );
    let again = write(&loaded);
    assert!(again.contains("MS:1000821"));
    assert_eq!(read(&again), loaded);
}
#[test]
fn all_typed_record_scalars_and_units_survive_without_reparsing_strings() {
    let mut e = experiment();
    let m = &mut e.spectra[0].metadata;
    m.insert("integer".into(), i64::MAX.into());
    m.insert(
        "double".into(),
        MetaValue::try_from(2.5)
            .unwrap()
            .with_unit(Unit::new("UO:0000010", "second", "UO").unwrap())
            .unwrap(),
    );
    m.insert("text".into(), "42".into());
    m.insert("list-looking".into(), "[1, 2]".into());
    m.insert("µ & α".into(), "< > &\n\t".into());
    let loaded = read(&write(&e));
    assert_eq!(loaded.spectra[0].metadata, e.spectra[0].metadata);
    assert!(matches!(
        loaded.spectra[0].metadata["text"].data(),
        MetaValueData::String(_)
    ));
}
#[test]
fn source_primary_metadata_overwrites_in_coordinate_then_intensity_order() {
    for reversed in [false, true] {
        let xml = write(&experiment());
        let mut a = vec![
            array(
                &[501., 499.],
                "MS:1000514",
                "",
                &user("key", "xsd:double", "2.5"),
            ),
            array(
                &[10., 20.],
                "MS:1000515",
                "",
                &user("key", "xsd:integer", "7"),
            ),
        ];
        if reversed {
            a.reverse();
        }
        let loaded = read(&arrays(&xml, &a));
        assert_eq!(loaded.spectra[0].metadata["key"].as_i64().unwrap(), 7);
    }
}
#[test]
fn mz_wins_over_wavelength_in_either_xml_order() {
    for reversed in [false, true] {
        let xml = write(&experiment());
        let mut a = vec![
            array(&[501., 499.], "MS:1000514", "", ""),
            array(&[10., 20.], "MS:1000515", "", ""),
            array(
                &[200., 220.],
                "MS:1000617",
                "unitAccession=\"UO:0000018\"",
                "",
            ),
        ];
        if reversed {
            a.reverse();
        }
        let s = &read(&arrays(&xml, &a)).spectra[0];
        assert_eq!(s.peaks[0].mz, 499.);
        assert!(!s.metadata.contains_key("mzml coordinate array"));
        assert_eq!(s.float_data_arrays[0].data, [220., 200.]);
    }
}
#[test]
fn noise_empty_mismatched_lengths_and_empty_spectrum_are_preserved() {
    let mut e = experiment();
    e.spectra[0].peaks.clear();
    for (i, v) in [vec![], vec![1., 2., 3.], vec![4.]].iter().enumerate() {
        e.spectra[0].metadata.insert(NOISE[i].into(), list(v));
    }
    let xml = write(&e);
    assert!(xml.contains("arrayLength=\"3\""));
    assert_eq!(read(&xml).spectra[0].metadata, e.spectra[0].metadata);
    let mut o = mzml::LoadOptions::default();
    o.scientific.fill_data = false;
    let loaded = mzml::read_with_load_options(Cursor::new(&xml), &o, &Default::default()).unwrap();
    assert!(loaded.spectra[0].metadata.is_empty());
}
#[test]
fn generic_lists_empty_values_wrong_noise_types_and_role_conflicts_fail_before_output() {
    for value in [
        MetaValue::default(),
        list(&[1.]),
        vec![1i64].into(),
        vec![String::from("x")].into(),
    ] {
        let mut e = experiment();
        e.spectra[0].metadata.insert("generic".into(), value);
        let mut output = b"unchanged".to_vec();
        assert!(mzml::write(&mut output, &e).is_err());
        assert_eq!(output, b"unchanged");
    }
    for (key, value) in [
        (NOISE[0], "[1]"),
        ("mzml coordinate array", "pressure"),
        ("mzml intensity array", "intensity"),
    ] {
        let mut e = experiment();
        e.spectra[0].metadata.insert(key.into(), value.into());
        assert!(mzml::write(Vec::new(), &e).is_err());
    }
}
#[test]
fn cpp_053_mixed_promoted_primary_intensities_are_rejected() {
    let mut e = MSExperiment::default();
    let c = MSChromatogram {
        peaks: vec![openms::kernel::ChromatogramPeak::new(1., 200.)],
        ..Default::default()
    };
    e.chromatograms.push(c);
    let xml = write(&e);
    let a = [
        array(&[1.], "MS:1000595", "unitAccession=\"UO:0000010\"", ""),
        array(&[200.], "MS:1000821", "", ""),
        array(&[3.], "MS:1000820", "", ""),
    ];
    assert!(mzml::read(Cursor::new(arrays(&xml, &a))).is_err());
}
#[test]
fn source_noise_descriptions_and_conflicting_units_cannot_disappear() {
    let xml = write(&experiment());
    let primary = [
        array(&[501., 499.], "MS:1000514", "", ""),
        array(&[10., 20.], "MS:1000515", "", ""),
    ];
    for tail in [
        array(
            &[1.],
            "MS:1002743",
            "",
            &user("annotation", "xsd:string", "x"),
        ),
        array(&[1.], "MS:1002743", "unitAccession=\"UO:0000018\"", ""),
    ] {
        let mut a = primary.to_vec();
        a.push(tail);
        assert!(mzml::read(Cursor::new(arrays(&xml, &a))).is_err());
    }
}
#[test]
fn noise_is_always_float64_ordinary_even_when_numpress_is_requested() {
    let mut e = experiment();
    e.spectra[0].metadata.insert(
        NOISE[0].into(),
        list(&[498.123456789012, 500.345678901234, 502.567890123456]),
    );
    let mut options = mzml::NumpressWriteOptions::default();
    options.mass_time.compression = openms::format::numpress_coder::NumpressCompression::Linear;
    let mut bytes = Vec::new();
    let report = mzml::write_with_numpress(&mut bytes, &e, &options).unwrap();
    assert_eq!(report.ordinary_arrays, 2);
    assert_eq!(
        read(std::str::from_utf8(&bytes).unwrap()).spectra[0].metadata,
        e.spectra[0].metadata
    );
    let mut options = PeakFileOptions::default();
    options.write_index = true;
    options.mz_32_bit = true;
    options.intensity_32_bit = true;
    let mut bytes = Vec::new();
    let report = mzml::write_with_peak_options(&mut bytes, &e, &options).unwrap();
    assert!(report.indexed);
    assert_eq!(
        read(std::str::from_utf8(&bytes).unwrap()).spectra[0].metadata,
        e.spectra[0].metadata
    );
}
#[test]
fn only_three_optical_scan_modes_accept_zero_ms_level() {
    for &mode in ScanMode::ALL {
        let mut s = spectrum();
        s.ms_level = 0;
        s.instrument_settings.scan_mode = mode;
        assert_eq!(
            s.validate().is_ok(),
            matches!(
                mode,
                ScanMode::ElectromagneticRadiation | ScanMode::Emission | ScanMode::Absorption
            )
        );
    }
    let mut c = MSChromatogram {
        chromatogram_type: ChromatogramType::TotalIonCurrent,
        ..Default::default()
    };
    c.metadata
        .insert("chromatogram type accession".into(), "MS:1003019".into());
    let mut e = MSExperiment::default();
    e.chromatograms.push(c);
    assert!(mzml::write(Vec::new(), &e).is_err());
}

#[test]
fn unused_wavelength_keeps_auxiliary_encounter_order_and_description() {
    let xml = write(&experiment());
    let a = vec![
        array(&[1., 2.], "MS:1000786", "value=\"before\"", ""),
        array(
            &[200., 220.],
            "MS:1000617",
            "unitAccession=\"UO:0000018\"",
            &user("label", "xsd:integer", "9"),
        ),
        array(&[3., 4.], "MS:1000786", "value=\"after\"", ""),
        array(&[501., 499.], "MS:1000514", "", ""),
        array(&[10., 20.], "MS:1000515", "", ""),
    ];
    let input = arrays(&xml, &a);
    for selected in [false, true] {
        let e = if selected {
            read(&input)
        } else {
            mzml::read(Cursor::new(&input)).unwrap()
        };
        let f = &e.spectra[0].float_data_arrays;
        assert_eq!(
            f.iter().map(|v| v.name.as_str()).collect::<Vec<_>>(),
            ["before", "wavelength array", "after"]
        );
        assert_eq!(f[1].metadata["label"].as_i64().unwrap(), 9);
        assert!(!f[1].metadata.contains_key("mzml coordinate array"));
        let again = read(&write(&e));
        assert_eq!(
            again.spectra[0].float_data_arrays[1].metadata,
            f[1].metadata
        );
    }
}
#[test]
fn chromatogram_primary_merge_follows_xml_order_for_all_promoted_roles() {
    let mut e = MSExperiment::default();
    e.chromatograms.push(MSChromatogram {
        peaks: vec![openms::kernel::ChromatogramPeak::new(1., 3.)],
        ..Default::default()
    });
    let xml = write(&e);
    for (id, extra, role) in [
        ("MS:1000821", "unitAccession=\"UO:0000109\"", "pressure"),
        ("MS:1000820", "unitAccession=\"UO:0000270\"", "flow"),
        (
            "MS:1000786",
            "value=\"detector signal\" unitAccession=\"UO:0000000\"",
            "nonstandard",
        ),
        ("MS:1000515", "unitAccession=\"UO:0000269\"", "absorption"),
    ] {
        for reverse in [false, true] {
            let mut a = vec![
                array(
                    &[1.],
                    "MS:1000595",
                    "unitAccession=\"UO:0000010\"",
                    &user("shared", "xsd:integer", "2"),
                ),
                array(&[3.], id, extra, &user("shared", "xsd:double", "4.5")),
            ];
            if reverse {
                a.reverse();
            }
            let e = read(&arrays(&xml, &a));
            let m = &e.chromatograms[0].metadata;
            assert_eq!(m["mzml intensity array"].as_str().unwrap(), role);
            if reverse {
                assert_eq!(m["shared"].as_i64().unwrap(), 2)
            } else {
                assert_eq!(m["shared"].as_f64().unwrap(), 4.5)
            }
            assert_eq!(read(&write(&e)), e);
        }
    }
}
#[test]
fn source_record_cv_routes_preserve_source_types_and_referenced_primary_metadata() {
    let mut xml = write_with_scan();
    xml = xml.replace(
        "<scanList",
        concat!(
            "<cvParam accession=\"MS:1000285\" name=\"total ion current\" value=\"3.5\"/>",
            "<cvParam accession=\"MS:1000797\" name=\"peak list scans\" value=\"12,13\"/><scanList"
        ),
    );
    xml=xml.replace("</scan>",concat!("<cvParam accession=\"MS:1000011\" name=\"mass resolution\" value=\"high\"/>","<cvParam accession=\"MS:1000015\" name=\"scan rate\" value=\"2.5\"/>","<cvParam accession=\"MS:1000826\" name=\"elution time\" value=\"1.5\" unitAccession=\"UO:0000031\"/>","<cvParam accession=\"MS:1000092\" name=\"decreasing m/z scan\"/></scan>"));
    let e = read(&xml);
    let m = &e.spectra[0].metadata;
    assert_eq!(m["total ion current"].as_f64().unwrap(), 3.5);
    assert_eq!(m["mass resolution"].as_str().unwrap(), "high");
    assert_eq!(m["scan rate"].as_f64().unwrap(), 2.5);
    assert_eq!(m["elution time (seconds)"].as_f64().unwrap(), 90.);
    assert_eq!(m["scan direction"].as_str().unwrap(), "decreasing");
    assert_eq!(read(&write(&e)).spectra[0].metadata, *m);
    let group = "<referenceableParamGroupList count=\"1\"><referenceableParamGroup id=\"role\"><cvParam accession=\"MS:1002743\" name=\"sampled noise m/z array\" unitAccession=\"MS:1000040\"/></referenceableParamGroup></referenceableParamGroupList>";
    let source = write_with_scan();
    let at = source.find("<softwareList").unwrap();
    let mut grouped = source.clone();
    grouped.insert_str(at, group);
    let noise = array(
        &[123.456789012345],
        "MS:1002743",
        "unitAccession=\"MS:1000040\"",
        "",
    )
    .replace(
        "<cvParam accession=\"MS:1002743\" name=\"array\" unitAccession=\"MS:1000040\"/>",
        "<referenceableParamGroupRef ref=\"role\"/>",
    );
    let a = [
        array(&[501., 499.], "MS:1000514", "", ""),
        array(&[10., 20.], "MS:1000515", "", ""),
        noise,
    ];
    assert_eq!(
        read(&arrays(&grouped, &a)).spectra[0].metadata[NOISE[0]]
            .as_float_list()
            .unwrap(),
        [123.456789012345]
    );
}
#[test]
fn noise_descriptor_lengths_and_cumulative_limits_remain_checked_without_population() {
    let xml = write(&experiment());
    let bad = array(&[1.], "MS:1002743", "", "").replace("MS:1000523", "MS:1000519");
    let a = [
        array(&[501., 499.], "MS:1000514", "", ""),
        array(&[10., 20.], "MS:1000515", "", ""),
        bad,
    ];
    let mut o = mzml::LoadOptions::default();
    o.scientific.fill_data = false;
    assert!(
        mzml::read_with_load_options(Cursor::new(arrays(&xml, &a)), &o, &Default::default())
            .is_err()
    );
    let mut e = experiment();
    e.spectra[0]
        .metadata
        .insert(NOISE[0].into(), list(&[1., 2., 3.]));
    let limits = mzml::ReadOptions {
        max_total_array_elements: 6,
        ..Default::default()
    };
    assert!(mzml::read_with_options(Cursor::new(write(&e)), &limits).is_err());
    let mut o = mzml::NumpressWriteOptions::default();
    o.limits.raw.max_values = 2;
    let mut out = b"same".to_vec();
    assert!(mzml::write_with_numpress(&mut out, &e, &o).is_err());
    assert_eq!(out, b"same");
}
#[test]
fn empty_unused_wavelength_is_an_empty_auxiliary_but_not_a_primary_coordinate() {
    let xml = write(&experiment());
    let w = array(&[], "MS:1000617", "", "");
    let a = [
        array(&[501., 499.], "MS:1000514", "", ""),
        array(&[10., 20.], "MS:1000515", "", ""),
        w.clone(),
    ];
    let e = read(&arrays(&xml, &a));
    assert!(e.spectra[0].float_data_arrays[0].data.is_empty());
    let mut omitted = mzml::LoadOptions::default();
    omitted.scientific.fill_data = false;
    let empty =
        mzml::read_with_load_options(Cursor::new(arrays(&xml, &a)), &omitted, &Default::default())
            .unwrap();
    assert!(empty.spectra[0].float_data_arrays.is_empty());
    let a = [w, array(&[10., 20.], "MS:1000515", "", "")];
    for fill in [false, true] {
        let mut o = mzml::LoadOptions::default();
        o.scientific.fill_data = fill;
        assert!(
            mzml::read_with_load_options(Cursor::new(arrays(&xml, &a)), &o, &Default::default())
                .is_err()
        );
    }
}
#[test]
fn promoted_role_ownership_is_rejected_before_all_writers_but_other_domains_remain_arrays() {
    for (kind, name) in NOISE
        .iter()
        .copied()
        .chain(["wavelength array", "pressure array"])
        .enumerate()
    {
        let mut e = experiment();
        match kind {
            0..=2 => {
                e.spectra[0]
                    .float_data_arrays
                    .push(DataArray::new(name, vec![1., 2.]));
            }
            3 => {
                e.spectra[0]
                    .metadata
                    .insert("mzml coordinate array".into(), "wavelength".into());
                e.spectra[0]
                    .float_data_arrays
                    .push(DataArray::new("wavelength array", vec![1., 2.]));
            }
            _ => {
                let mut c = MSChromatogram::default();
                c.float_data_arrays
                    .push(DataArray::new("pressure array", vec![]));
                e.chromatograms.push(c);
            }
        }
        for writer in 0..3 {
            let mut out = b"same".to_vec();
            let result = match writer {
                0 => mzml::write(&mut out, &e),
                1 => mzml::write_with_numpress(&mut out, &e, &Default::default()).map(|_| ()),
                _ => mzml::write_with_peak_options(&mut out, &e, &Default::default()).map(|_| ()),
            };
            assert!(result.is_err());
            assert_eq!(out, b"same");
        }
    }
    let mut e = experiment();
    e.spectra[0]
        .float_data_arrays
        .push(DataArray::new("pressure array", vec![1., 2.]));
    let mut c = MSChromatogram::default();
    c.float_data_arrays.push(DataArray::new(NOISE[0], vec![]));
    e.chromatograms.push(c);
    let back = read(&write(&e));
    assert_eq!(back.spectra[0].float_data_arrays[0].name, "pressure array");
    assert_eq!(back.chromatograms[0].float_data_arrays[0].name, NOISE[0]);
    assert!(!back.chromatograms[0].metadata.contains_key(NOISE[0]));
}
#[test]
fn typed_noise_and_primary_roles_reach_path_consumer_and_count_readers() {
    use openms::{Result, interfaces::MSDataConsumer, metadata::ExperimentalSettings};
    use std::ops::ControlFlow;
    #[derive(Default)]
    struct Collector(Vec<MSSpectrum>);
    impl MSDataConsumer for Collector {
        fn set_expected_size(&mut self, s: usize, c: usize) -> Result<()> {
            assert_eq!((s, c), (1, 0));
            Ok(())
        }
        fn set_experimental_settings(&mut self, _: &ExperimentalSettings) -> Result<()> {
            Ok(())
        }
        fn consume_spectrum(&mut self, s: &mut MSSpectrum) -> Result<ControlFlow<()>> {
            self.0.push(s.clone());
            Ok(ControlFlow::Continue(()))
        }
        fn consume_chromatogram(&mut self, _: &mut MSChromatogram) -> Result<ControlFlow<()>> {
            panic!("no chromatograms")
        }
    }
    let mut e = experiment();
    e.spectra[0].ms_level = 0;
    e.spectra[0].instrument_settings.scan_mode = ScanMode::Absorption;
    e.spectra[0]
        .metadata
        .insert("mzml coordinate array".into(), "wavelength".into());
    e.spectra[0]
        .metadata
        .insert(NOISE[0].into(), list(&[1., 2., 3.]));
    e.spectra[0].metadata.insert("typed".into(), 42_i64.into());
    let xml = write(&e);
    let expected = read(&xml);
    let counts = mzml::read_size(Cursor::new(&xml)).unwrap();
    assert_eq!(
        counts,
        mzml::MzMLCounts {
            spectra: 1,
            chromatograms: 0
        }
    );
    let file = openms::system::file::TempFile::new().unwrap();
    std::fs::write(file.path(), &xml).unwrap();
    let loaded = mzml::load(file.path()).unwrap();
    assert_eq!(loaded.spectra, expected.spectra);
    assert_eq!(mzml::load_size(file.path()).unwrap(), counts);
    let mut consumer = Collector::default();
    let mut destination = MSExperiment::default();
    mzml::transform_from_into(
        || Ok(Cursor::new(xml.as_bytes())),
        &mut consumer,
        &mut destination,
        &Default::default(),
    )
    .unwrap();
    assert_eq!(consumer.0, expected.spectra);
    assert_eq!(destination.spectra, expected.spectra);
}
#[test]
fn noise_does_not_consume_described_auxiliary_header_slots_in_any_writer() {
    use openms::metadata::DataProcessing;
    use std::sync::Arc;
    let mut e = experiment();
    e.spectra[0]
        .metadata
        .insert(NOISE[0].into(), list(&[1., 2., 3.]));
    let mut a = DataArray::new("auxiliary", vec![7., 9.]);
    a.metadata.insert("typed".into(), 42_i64.into());
    let mut dp = DataProcessing::default();
    dp.software.name = "independent".into();
    dp.software.version = "1".into();
    a.data_processing.push(Arc::new(dp));
    e.spectra[0].float_data_arrays.push(a);
    for writer in 0..3 {
        let mut out = Vec::new();
        match writer {
            0 => mzml::write(&mut out, &e).unwrap(),
            1 => {
                mzml::write_with_numpress(&mut out, &e, &Default::default()).unwrap();
            }
            _ => {
                let mut o = PeakFileOptions::default();
                o.write_index = true;
                mzml::write_with_peak_options(&mut out, &e, &o).unwrap();
            }
        }
        let xml = String::from_utf8(out).unwrap();
        let back = read(&xml);
        let a = &back.spectra[0].float_data_arrays[0];
        assert_eq!(a.data, [9., 7.]);
        assert_eq!(a.metadata, e.spectra[0].float_data_arrays[0].metadata);
        assert_eq!(
            a.data_processing,
            e.spectra[0].float_data_arrays[0].data_processing
        );
        assert_eq!(
            back.spectra[0].metadata[NOISE[0]],
            e.spectra[0].metadata[NOISE[0]]
        );
        let noise = xml.find("name=\"sampled noise m/z array\"").unwrap();
        let aux = xml.find("value=\"auxiliary\"").unwrap();
        assert!(noise < aux);
    }
}

#[test]
fn complete_source_spectrum_and_scan_cv_route_types() {
    // Literal routes from MzMLHandler; types checked against pinned PSI-MS xrefs.
    let spectrum_routes = [
        ("MS:1000285", "total ion current", true),
        ("MS:1000504", "base peak m/z", true),
        ("MS:1000505", "base peak intensity", true),
        ("MS:1000527", "highest observed m/z", true),
        ("MS:1000528", "lowest observed m/z", true),
        ("MS:1000618", "highest observed wavelength", true),
        ("MS:1000619", "lowest observed wavelength", true),
        ("MS:1000796", "spectrum title", false),
        ("MS:1000797", "peak list scans", false),
        ("MS:1000798", "peak list raw scans", false),
    ];
    let scan_routes = [
        ("MS:1000502", "dwell time", true),
        ("MS:1000011", "mass resolution", false),
        ("MS:1000015", "scan rate", true),
        ("MS:1000512", "filter string", false),
        ("MS:1000803", "analyzer scan offset", true),
        ("MS:1000616", "preset scan configuration", false),
        ("MS:1000800", "mass resolving power", false),
        ("MS:1000880", "interchannel delay", true),
    ];
    let mut xml = write_with_scan();
    let cv = |id: &str, name: &str| {
        format!("<cvParam accession=\"{id}\" name=\"{name}\" value=\"42.5\"/>")
    };
    xml = xml.replace(
        "<scanList",
        &(spectrum_routes
            .iter()
            .map(|(id, n, _)| cv(id, n))
            .collect::<String>()
            + "<scanList"),
    );
    xml = xml.replace(
        "</scan>",
        &(scan_routes
            .iter()
            .map(|(id, n, _)| cv(id, n))
            .collect::<String>()
            + "</scan>"),
    );
    let e = read(&xml);
    for (_, name, numeric) in spectrum_routes.into_iter().chain(scan_routes) {
        if numeric {
            assert_eq!(e.spectra[0].metadata[name].as_f64().unwrap(), 42.5)
        } else {
            assert_eq!(e.spectra[0].metadata[name].as_str().unwrap(), "42.5")
        }
    }
    assert_eq!(read(&write(&e)).spectra[0].metadata, e.spectra[0].metadata);
    for (id, name, value) in [
        ("MS:1000092", "scan direction", "decreasing"),
        ("MS:1000093", "scan direction", "increasing"),
        ("MS:1000094", "scan law", "exponential"),
        ("MS:1000095", "scan law", "linear"),
        ("MS:1000096", "scan law", "quadratic"),
    ] {
        let xml = write_with_scan().replace("</scan>", &(cv(id, "source term") + "</scan>"));
        assert_eq!(
            read(&xml).spectra[0].metadata[name].as_str().unwrap(),
            value
        );
    }
    for (kind, value, integer) in [
        ("xsd:integer", "-7", true),
        ("xsd:int", "42", true),
        ("xsd:double", "42", false),
        ("xsd:float", "2.5", false),
    ] {
        let xml =
            write_with_scan().replace("<scanList", &(user("typed", kind, value) + "<scanList"));
        let e = read(&xml);
        assert_eq!(
            matches!(
                e.spectra[0].metadata["typed"].data(),
                MetaValueData::Integer(_)
            ),
            integer
        );
    }
}
#[test]
fn source_elution_time_fallback_precedes_primary_merge_and_does_not_create_rt_filter_event() {
    let source = write_with_scan();
    let needle = "<cvParam cvRef=\"MS\" accession=\"MS:1000016\"";
    let start = source.find(needle).unwrap();
    let end = start + source[start..].find("/>").unwrap() + 2;
    let no_rt = format!("{}{}", &source[..start], &source[end..]);
    let elution = "<cvParam accession=\"MS:1000826\" name=\"elution time\" value=\"1.5\" unitAccession=\"UO:0000031\"/>";
    for fill in [false, true] {
        let mut options = mzml::LoadOptions::default();
        options.scientific.fill_data = fill;
        options
            .scientific
            .set_rt_range(NumericRange { min: 0., max: 2. });
        let xml = no_rt.replace("</scan>", &(elution.to_owned() + "</scan>"));
        let e =
            mzml::read_with_load_options(Cursor::new(xml), &options, &Default::default()).unwrap();
        assert_eq!(e.spectra[0].rt, 90.);
        assert_eq!(
            e.spectra[0].metadata["elution time (seconds)"]
                .as_f64()
                .unwrap(),
            90.
        );
    }
    let explicit = source.replace("</scan>", &(elution.to_owned() + "</scan>"));
    assert_eq!(read(&explicit).spectra[0].rt, 1.);
    let a = [
        array(
            &[501., 499.],
            "MS:1000514",
            "",
            &user("elution time (seconds)", "xsd:double", "8"),
        ),
        array(&[10., 20.], "MS:1000515", "", ""),
    ];
    let e = read(&arrays(&no_rt, &a));
    assert_eq!(e.spectra[0].rt, -1.);
    assert_eq!(
        e.spectra[0].metadata["elution time (seconds)"]
            .as_f64()
            .unwrap(),
        8.
    );
}
#[test]
fn cpp_058_elution_fallback_accepts_numeric_values_without_reading_inactive_string_storage() {
    for value in [MetaValue::from(90_i64), MetaValue::try_from(90.).unwrap()] {
        let mut e = experiment();
        e.spectra[0]
            .metadata
            .insert("elution time (seconds)".into(), value);
        assert_eq!(read(&write(&e)).spectra[0].rt, 90.);
    }
    // Direct String input reproduces the source inactive-union trigger;
    // this checked native regression makes no C++ execution claim.
    let mut e = experiment();
    e.spectra[0]
        .metadata
        .insert("elution time (seconds)".into(), "90".into());
    assert!(mzml::read(Cursor::new(write(&e))).is_err());
    e.spectra[0].rt = 1.;
    assert_eq!(read(&write(&e)).spectra[0].rt, 1.);
}
