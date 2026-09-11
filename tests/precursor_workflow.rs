// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::analysis::precursor_purity::PrecursorPurity;
use openms::comparison::Tolerance;
use openms::format::{dta, mgf};
use openms::metadata::PrecursorInfo;
#[cfg(feature = "mzml")]
use openms::metadata::{ActivationMethod, DriftTimeUnit, MetaValue};
use openms::{MSExperiment, MSSpectrum, Peak1D, Precursor};

fn experiment() -> MSExperiment {
    let mut precursor = Precursor::new(500.0, 2);
    precursor.isolation_window_lower_offset = 0.3;
    precursor.isolation_window_upper_offset = 0.8;
    precursor.spectrum_reference = Some("scan=1".into());
    MSExperiment {
        spectra: vec![
            MSSpectrum {
                native_id: "scan=1".into(),
                rt: 10.0,
                peaks: vec![
                    Peak1D::new(500.0, 10.0),
                    Peak1D::new(500.2, 3.0),
                    Peak1D::new(500.5016774189, 5.0),
                ],
                ..Default::default()
            },
            MSSpectrum {
                native_id: "scan=2".into(),
                rt: 10.2,
                ms_level: 2,
                precursors: vec![precursor],
                ..Default::default()
            },
        ],
        ..Default::default()
    }
}

#[test]
fn native_precursor_and_compatibility_wrapper_share_acquisition_state() {
    let mut info = PrecursorInfo::from(Precursor::new(500.0, 2));
    info.isolation_window_lower_offset = 0.2;
    info.activation_energy = 35.0;
    info.spectrum_reference = Some("scan=9".into());
    assert_eq!(info.peak.activation_energy, 35.0);
    let p: Precursor = info.into();
    assert_eq!(p.isolation_window().unwrap().min, 499.8);
    let mut s = MSSpectrum {
        precursors: vec![p],
        ..Default::default()
    };
    s.validate().unwrap();
    let clone = s.clone();
    s.precursors[0].activation_energy = 5.0;
    assert_eq!(clone.precursors[0].activation_energy, 35.0);
    s.precursors[0].isolation_window_lower_offset = f64::NAN;
    assert!(s.validate().is_err());
    let info = PrecursorInfo::from(clone.precursors[0].clone());
    assert_eq!(info.spectrum_reference.as_deref(), Some("scan=9"));
}

#[test]
fn kernel_parent_lookup_uses_level_reference_and_acquisition_order() {
    let mut e = experiment();
    e.spectra.push(MSSpectrum {
        native_id: "scan=3".into(),
        ..Default::default()
    });
    let mut child = e.spectra[1].clone();
    child.native_id = "scan=4".into();
    e.spectra.push(child);
    assert_eq!(e.precursor_spectrum_index(3).unwrap(), Some(0));
    e.spectra[3].precursors[0].spectrum_reference = Some("scan=999".into());
    assert_eq!(e.precursor_spectrum_index(3).unwrap(), Some(2));
    e.spectra[3].precursors[0].spectrum_reference = Some("scan=2".into());
    assert_eq!(e.precursor_spectrum_index(3).unwrap(), Some(2));
    e.spectra[3].ms_level = 3;
    assert_eq!(e.precursor_spectrum_index(3).unwrap(), Some(1));
    assert_eq!(e.precursor_spectrum_index(0).unwrap(), None);
    assert!(e.precursor_spectrum_index(4).is_err());
}

#[test]
fn unsupported_peak_list_exports_fail_before_writing_acquisition_data() {
    let e = experiment();
    let mut output = vec![42];
    assert!(dta::write(&mut output, &e.spectra[1]).is_err());
    assert_eq!(output, [42]);
    assert!(mgf::write(&mut output, &e).is_err());
    assert_eq!(output, [42]);
    let mut basic = e.clone();
    basic.spectra[1].precursors[0] = Precursor::new(500.0, 2);
    dta::write(Vec::new(), &basic.spectra[1]).unwrap();
    mgf::write(Vec::new(), &basic).unwrap();
}

#[test]
fn batch_checked_boundaries_preserve_source_first_precursor_and_missing_parent_defaults() {
    let e = experiment();
    let s = PrecursorPurity::compute_all(&e, Tolerance::Absolute(0.001), false).unwrap();
    assert_eq!(s["scan=2"].total_intensity, 18.0);
    assert_eq!(s["scan=2"].target_intensity, 15.0);
    let mut multiple = e.clone();
    multiple.spectra[1]
        .precursors
        .push(Precursor::new(400.0, 1));
    assert_eq!(
        s,
        PrecursorPurity::compute_all(&multiple, Tolerance::Absolute(0.001), false).unwrap()
    );
    let mut absent = e.clone();
    absent.spectra.remove(0);
    assert!(PrecursorPurity::compute_all(&absent, Tolerance::Absolute(0.001), false).is_err());
    assert_eq!(
        PrecursorPurity::compute_all(&absent, Tolerance::Absolute(0.001), true).unwrap()["scan=2"]
            .total_intensity,
        0.0
    );
    let mut duplicate = e.clone();
    duplicate.spectra.push(e.spectra[1].clone());
    assert!(PrecursorPurity::compute_all(&duplicate, Tolerance::Absolute(0.001), false).is_err());
    assert!(
        PrecursorPurity::compute_all(&MSExperiment::default(), Tolerance::Absolute(0.001), false)
            .unwrap()
            .is_empty()
    );
}

#[cfg(feature = "mzml")]
#[test]
fn mzml_roundtrip_preserves_purity_acquisition_and_distinct_isolation_target() {
    use openms::format::mzml;
    let mut e = experiment();
    let p = &mut e.spectra[1].precursors[0];
    p.isolation_target_mz = Some(500.1);
    p.activation_energy = 27.5;
    p.activation_methods
        .extend([ActivationMethod::Cid, ActivationMethod::Ethcd]);
    p.possible_charge_states = vec![2, 3, 0, 2];
    p.drift_time = Some(-45.0);
    p.drift_time_unit = DriftTimeUnit::FaimsCompensationVoltage;
    let scores = PrecursorPurity::compute_all(&e, Tolerance::Absolute(0.001), false).unwrap();
    for compressed in [false, true] {
        let mut bytes = Vec::new();
        mzml::write_with_options(
            &mut bytes,
            &e,
            &mzml::WriteOptions {
                zlib_compression: compressed,
            },
        )
        .unwrap();
        let decoded = mzml::read(bytes.as_slice()).unwrap();
        assert_eq!(decoded, e);
        assert_eq!(
            PrecursorPurity::compute_all(&decoded, Tolerance::Absolute(0.001), false).unwrap(),
            scores
        );
        assert_eq!(decoded.precursor_spectrum_index(1).unwrap(), Some(0));
    }
}

#[cfg(feature = "mzml")]
#[test]
fn mzml_all_activation_methods_and_mobility_quantities_roundtrip() {
    use openms::format::mzml;
    for (unit, value) in [
        (DriftTimeUnit::Millisecond, 12.3),
        (DriftTimeUnit::InverseReducedMobility, 0.83),
        (DriftTimeUnit::CollisionCrossSection, 148.1),
        (DriftTimeUnit::FaimsCompensationVoltage, -55.0),
    ] {
        let mut e = experiment();
        let p = &mut e.spectra[1].precursors[0];
        p.activation_methods.extend(ActivationMethod::ALL);
        p.drift_time = Some(value);
        p.drift_time_unit = unit;
        e.chromatograms.push(openms::MSChromatogram {
            native_id: "chromatogram=1".into(),
            precursor: p.clone(),
            ..Default::default()
        });
        let mut bytes = Vec::new();
        mzml::write(&mut bytes, &e).unwrap();
        assert_eq!(mzml::read(bytes.as_slice()).unwrap(), e);
    }
}

#[cfg(feature = "mzml")]
#[test]
fn mzml_acquisition_parser_rejects_conflicting_fields_units_and_nesting() {
    use openms::format::mzml;
    let mut e = experiment();
    e.spectra[1].precursors[0].activation_energy = 25.0;
    let mut bytes = Vec::new();
    mzml::write(&mut bytes, &e).unwrap();
    let xml = String::from_utf8(bytes).unwrap();
    let start = xml
        .find("<cvParam cvRef=\"MS\" accession=\"MS:1000828\"")
        .unwrap();
    let end = start + xml[start..].find("/>").unwrap() + 2;
    let offset = &xml[start..end];
    for bad in [
        xml.replacen(offset, &format!("{offset}{offset}"), 1),
        xml.replacen(
            offset,
            &offset.replace("value=\"0.3\"", "value=\"-0.3\""),
            1,
        ),
        xml.replacen(offset, &offset.replace("MS:1000040", "UO:0000010"), 1),
        xml.replace("<isolationWindow>", "<selectedIon><isolationWindow>")
            .replace("</isolationWindow>", "</isolationWindow></selectedIon>"),
        xml.replace("<activation>", "<activation><activation>")
            .replace("</activation>", "</activation></activation>"),
    ] {
        assert!(
            mzml::read(bad.as_bytes()).is_err(),
            "accepted malformed acquisition XML"
        );
    }
}

#[cfg(feature = "mzml")]
#[test]
fn mzml_writer_rejects_unrepresented_precursor_fields_and_missing_reference_atomically() {
    use openms::format::mzml;
    let e = experiment();
    let mut cases = Vec::new();
    let mut negative = e.clone();
    negative.spectra[1].precursors[0].mz = -1.0;
    negative.spectra[1].precursors[0].isolation_target_mz = None;
    cases.push(negative.clone());
    let mut chromatogram = e.clone();
    chromatogram.chromatograms.push(openms::MSChromatogram {
        precursor: negative.spectra[1].precursors[0].clone(),
        ..Default::default()
    });
    cases.push(chromatogram);
    let mut v = e.clone();
    v.spectra[1].precursors[0].spectrum_reference = Some("scan=999".into());
    cases.push(v);
    let mut v = e.clone();
    v.spectra[1].precursors[0].drift_window_lower_offset = 0.1;
    cases.push(v);
    let mut v = e.clone();
    v.spectra[1].precursors[0].drift_time = Some(1.0);
    cases.push(v);
    let mut v = e.clone();
    v.spectra[1].precursors[0].drift_time_unit = DriftTimeUnit::Millisecond;
    cases.push(v);
    let mut v = e.clone();
    v.spectra[1].precursors[0]
        .cv_terms
        .metadata
        .insert("extra".into(), MetaValue::from("value"));
    // The header codec now carries scalar precursor metadata. List values
    // still have no reversible mzML scalar representation.
    let mut output = Vec::new();
    mzml::write(&mut output, &v).unwrap();
    assert_eq!(
        mzml::read(output.as_slice()).unwrap().spectra[1].precursors[0],
        v.spectra[1].precursors[0]
    );
    v.spectra[1].precursors[0]
        .cv_terms
        .metadata
        .insert("extra".into(), MetaValue::from(vec!["value".to_owned()]));
    cases.push(v);
    for bad in cases {
        let mut output = vec![42];
        assert!(mzml::write(&mut output, &bad).is_err());
        assert_eq!(output, [42]);
    }
}

#[cfg(feature = "mzml")]
#[test]
fn original_mzml_precursors_reproduce_independent_scalar_and_batch_values() {
    use openms::format::mzml;
    let original = include_str!("data/precursor_purity_input.mzML");
    // The untouched source fixture is ASCII but declares Latin-1. Change only
    // its declaration in memory to the reader's supported UTF-8 encoding.
    assert!(original.is_ascii());
    let xml = original.replacen("encoding=\"ISO-8859-1\"", "encoding=\"UTF-8\"", 1);
    let e = mzml::read(xml.as_bytes()).unwrap();
    assert_eq!(e.spectra.len(), 7);
    assert_eq!(e.spectra[0].peaks.len(), 2003);
    assert_eq!(e.spectra[6].peaks.len(), 1869);
    for index in 1..=5 {
        let p = &e.spectra[index].precursors[0];
        assert_eq!(p.isolation_window_lower_offset, 0.25);
        assert_eq!(p.isolation_window_upper_offset, 0.25);
        assert_eq!(p.isolation_target_mz, None);
        assert!(p.activation_methods.contains(&ActivationMethod::Hcd));
        assert_eq!(e.precursor_spectrum_index(index).unwrap(), Some(0));
    }
    let batch = PrecursorPurity::compute_all(&e, Tolerance::Absolute(0.1), false).unwrap();
    for row in include_str!("data/precursor_purity_scores.tsv")
        .lines()
        .filter(|line| !line.starts_with('#'))
        .skip(1)
    {
        let f: Vec<_> = row.split('\t').collect();
        let ms1: usize = f[1].parse().unwrap();
        let ms2: usize = f[2].parse().unwrap();
        let tolerance = if f[4] == "ppm" {
            Tolerance::Ppm(f[3].parse().unwrap())
        } else {
            Tolerance::Absolute(f[3].parse().unwrap())
        };
        let score =
            PrecursorPurity::compute(&e.spectra[ms1], &e.spectra[ms2].precursors[0], tolerance)
                .unwrap();
        assert_eq!(score.total_intensity, f[8].parse::<f64>().unwrap());
        assert_eq!(score.target_intensity, f[9].parse::<f64>().unwrap());
        assert!((score.signal_proportion - f[10].parse::<f64>().unwrap()).abs() < 1e-14);
        if f[0].starts_with("experiment_") {
            assert_eq!(batch[&e.spectra[ms2].native_id], score);
        }
    }
}

#[cfg(feature = "mzml")]
#[test]
fn populated_acquisition_output_validates_against_independent_schema_when_available() {
    use openms::format::mzml;
    use std::io::Write;
    use std::process::{Command, Stdio};
    if Command::new("xmllint").arg("--version").output().is_err() {
        eprintln!("xmllint unavailable; precursor writer XSD validation not executed on this host");
        return;
    }
    let schema = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/mzml_1_10.xsd");
    let mut e = experiment();
    let p = &mut e.spectra[1].precursors[0];
    p.activation_methods.extend(ActivationMethod::ALL);
    p.activation_energy = 35.0;
    p.isolation_target_mz = Some(500.1);
    p.possible_charge_states = vec![2, 3];
    p.drift_time = Some(-45.0);
    p.drift_time_unit = DriftTimeUnit::FaimsCompensationVoltage;
    for compressed in [false, true] {
        let mut bytes = Vec::new();
        mzml::write_with_options(
            &mut bytes,
            &e,
            &mzml::WriteOptions {
                zlib_compression: compressed,
            },
        )
        .unwrap();
        let mut child = Command::new("xmllint")
            .args(["--nonet", "--noout", "--schema"])
            .arg(&schema)
            .arg("-")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        child.stdin.take().unwrap().write_all(&bytes).unwrap();
        let result = child.wait_with_output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
}
