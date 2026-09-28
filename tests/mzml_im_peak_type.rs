// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Three `MzMLFile` routes against the Linux x86_64 Release build: the
//! ion-mobility peak type, selected-ion mobility copied onto the spectrum, and
//! the activation route for dissociation methods without a typed member.
//!
//! `data/mzml_im_peak_type/oracle.tsv` is `results/records.tsv` of
//! `../oracle/B8a-mzml-im`: its driver (`driver.cpp`) links the Release install
//! `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576` on
//! ibminode06 and runs, on the three inputs of this directory (written by the
//! oracle's `make_inputs.py`), `MzMLFile::load`, `load` with `fill_data=false`,
//! `transform` into a consumer, and `store` followed by `load`; plus a store of
//! three spectra built in memory. Two runs were byte-identical. Every record is
//! `case <TAB> file <TAB> kind <TAB> fields`; doubles are C `%a` text.
//!
//! What the Release build does that this port does not by default is selected
//! by `mzml::ReadOptions::source()` and pinned here in both settings:
//! `source_im_centroid_frame` and `source_selected_ion_drift_time`.
#![cfg(feature = "mzml")]

use openms::format::mzml::{self, LoadOptions, MSDataConsumer, ReadOptions, TransformOptions};
use openms::kernel::{DataArray, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, Precursor};
use openms::metadata::{
    ActivationMethod, DriftTimeUnit, ExperimentalSettings, IonMobilityPeakType, MetaValue,
    MetaValueData, im_peak_type_to_string,
};
use std::io::BufReader;
use std::ops::ControlFlow;
use std::path::PathBuf;

const ORACLE: &str = include_str!("data/mzml_im_peak_type/oracle.tsv");

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/mzml_im_peak_type")
        .join(name)
}

/// The oracle rows of one case, file and kind, as fields after the kind.
fn rows(case: &str, file: &str, kind: &str) -> Vec<Vec<&'static str>> {
    let rows: Vec<Vec<&str>> = ORACLE
        .lines()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|f| f[0] == case && f[1] == file && f[2] == kind)
        .map(|f| f[3..].to_vec())
        .collect();
    assert!(!rows.is_empty(), "no oracle rows for {case} {file} {kind}");
    rows
}

/// Parse a C `printf("%a")` hexadecimal float exactly.
fn hex(text: &str) -> f64 {
    let (negative, body) = match text.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, text),
    };
    let body = body.strip_prefix("0x").expect("hexadecimal float");
    let (mantissa, exponent) = body.split_once('p').expect("binary exponent");
    let exponent: i32 = exponent.parse().expect("decimal exponent");
    let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    let digits = u64::from_str_radix(&format!("{whole}{fraction}"), 16).expect("hex digits");
    let scale = exponent - 4 * i32::try_from(fraction.len()).expect("short fraction");
    let magnitude = digits as f64 * 2f64.powi(scale);
    if negative { -magnitude } else { magnitude }
}

/// The source `DriftTimeUnit` enumerator value (`IONMOBILITY/IMTypes.h:22-30`).
fn unit_code(unit: DriftTimeUnit) -> &'static str {
    match unit {
        DriftTimeUnit::None => "0",
        DriftTimeUnit::Millisecond => "1",
        DriftTimeUnit::InverseReducedMobility => "2",
        DriftTimeUnit::FaimsCompensationVoltage => "3",
        DriftTimeUnit::CollisionCrossSection => "4",
    }
}

/// The driver's `valueType:toString(true)[unitType:unit]` rendering.
fn meta_text(key: &str, value: &MetaValue) -> String {
    let kind = match value.data() {
        MetaValueData::String(_) => 0,
        MetaValueData::Integer(_) => 1,
        MetaValueData::Float(_) => 2,
        other => panic!("unexpected precursor meta value {other:?}"),
    };
    let mut text = format!("{key}={kind}:{}", value.to_text(true).unwrap());
    if let Some(unit) = value.unit() {
        let (ontology, number) = unit.accession().split_once(':').unwrap();
        let ontology = match ontology {
            "UO" => 0,
            "MS" => 1,
            other => panic!("unit ontology {other}"),
        };
        text.push_str(&format!("[{ontology}:{}]", number.parse::<u32>().unwrap()));
    }
    text
}

/// The driver's `drift|unit|methods|meta` rendering of a precursor. The
/// source's `isolation window target m/z` meta value is its own bookkeeping
/// (this port keeps the target in `isolation_target_mz`), so it is left out
/// on both sides.
fn precursor_text(p: &Precursor) -> String {
    let drift = p.drift_time.unwrap_or(-1.0);
    let methods: Vec<String> = p
        .activation_methods
        .iter()
        .map(|m| {
            ActivationMethod::ALL
                .iter()
                .position(|a| a == m)
                .unwrap()
                .to_string()
        })
        .collect();
    let meta: Vec<String> = p
        .cv_terms
        .metadata
        .iter()
        .map(|(key, value)| meta_text(key, value))
        .collect();
    format!(
        "{drift:e}|{}|{}|{}",
        unit_code(p.drift_time_unit),
        if methods.is_empty() {
            "-".into()
        } else {
            methods.join(",")
        },
        if meta.is_empty() {
            "-".into()
        } else {
            meta.join(";")
        }
    )
}

/// The oracle's precursor field in the same form as [`precursor_text`].
fn oracle_precursor(field: &str) -> String {
    let mut parts = field.splitn(4, '|');
    let drift = hex(parts.next().unwrap());
    let unit = parts.next().unwrap();
    let methods = parts.next().unwrap();
    let meta: Vec<&str> = parts
        .next()
        .unwrap()
        .split(';')
        .filter(|m| *m != "-" && !m.starts_with("isolation window target m/z="))
        .collect();
    format!(
        "{drift:e}|{unit}|{methods}|{}",
        if meta.is_empty() {
            "-".into()
        } else {
            meta.join(";")
        }
    )
}

/// One spectrum as `[native id, im peak type, drift, unit, precursors...]`,
/// the oracle's `spectrum` fields after the index.
fn spectrum_fields(s: &MSSpectrum) -> Vec<String> {
    let mut fields = vec![
        s.native_id.clone(),
        im_peak_type_to_string(s.im_peak_type).to_owned(),
        format!("{:e}", s.drift_time),
        unit_code(s.drift_time_unit).to_owned(),
        s.precursors.len().to_string(),
    ];
    fields.extend(s.precursors.iter().map(precursor_text));
    fields
}

fn oracle_spectrum_fields(row: &[&str]) -> Vec<String> {
    let mut fields = vec![
        row[1].to_owned(),
        row[2].to_owned(),
        format!("{:e}", hex(row[3])),
        row[4].to_owned(),
        row[5].to_owned(),
    ];
    fields.extend(row[6..].iter().map(|f| oracle_precursor(f)));
    fields
}

fn assert_spectra(spectra: &[MSSpectrum], case: &str, file: &str, kind: &str) {
    let expected = rows(case, file, kind);
    assert_eq!(spectra.len(), expected.len(), "{case} {file} {kind}");
    for (index, (spectrum, row)) in spectra.iter().zip(&expected).enumerate() {
        assert_eq!(row[0], index.to_string());
        assert_eq!(
            spectrum_fields(spectrum),
            oracle_spectrum_fields(row),
            "{case} {file} {kind} spectrum {index}"
        );
    }
}

fn read(file: &str, options: &ReadOptions) -> MSExperiment {
    let reader = BufReader::new(std::fs::File::open(fixture(file)).unwrap());
    mzml::read_with_options(reader, options).unwrap()
}

fn read_without_data(file: &str, options: &ReadOptions) -> MSExperiment {
    let mut load = LoadOptions::default();
    load.scientific.fill_data = false;
    let reader = BufReader::new(std::fs::File::open(fixture(file)).unwrap());
    mzml::read_with_load_options(reader, &load, options).unwrap()
}

#[derive(Default)]
struct Collect(Vec<MSSpectrum>);

impl MSDataConsumer for Collect {
    fn set_expected_size(&mut self, _: usize, _: usize) -> openms::Result<()> {
        Ok(())
    }
    fn set_experimental_settings(&mut self, _: &ExperimentalSettings) -> openms::Result<()> {
        Ok(())
    }
    fn consume_spectrum(&mut self, spectrum: &mut MSSpectrum) -> openms::Result<ControlFlow<()>> {
        self.0.push(spectrum.clone());
        Ok(ControlFlow::Continue(()))
    }
    fn consume_chromatogram(&mut self, _: &mut MSChromatogram) -> openms::Result<ControlFlow<()>> {
        Ok(ControlFlow::Continue(()))
    }
}

fn transformed(file: &str, read: ReadOptions) -> Vec<MSSpectrum> {
    let mut collect = Collect::default();
    let options = TransformOptions {
        read,
        ..TransformOptions::default()
    };
    mzml::transform_with_options(fixture(file), &mut collect, &options).unwrap();
    collect.0
}

fn written(experiment: &MSExperiment) -> String {
    let mut bytes = Vec::new();
    mzml::write(&mut bytes, experiment).unwrap();
    String::from_utf8(bytes).unwrap()
}

fn reread(xml: &str, options: &ReadOptions) -> MSExperiment {
    mzml::read_with_options(xml.as_bytes(), options).unwrap()
}

/// The text of each written `<spectrum>` (or `<chromatogram>`) element, in
/// document order.
fn records<'a>(xml: &'a str, tag: &str) -> Vec<&'a str> {
    let open = format!("<{tag} ");
    let close = format!("</{tag}>");
    let mut out = Vec::new();
    let mut rest = xml;
    while let Some(start) = rest.find(&open) {
        let end = rest[start..].find(&close).unwrap() + start;
        out.push(&rest[start..end]);
        rest = &rest[end..];
    }
    out
}

/// The oracle's written lines, grouped by the `<spectrum` line that opens them.
fn oracle_written(case: &str, file: &str) -> Vec<Vec<&'static str>> {
    let mut groups: Vec<Vec<&str>> = Vec::new();
    for row in rows(case, file, "written") {
        let line = row[0];
        if line.starts_with("<chromatogram ") {
            break;
        }
        if line.starts_with("<spectrum ") {
            groups.push(Vec::new());
        } else {
            groups.last_mut().unwrap().push(line);
        }
    }
    groups
}

// --------------------------------------------------------------------------
// Ion-mobility peak type
// --------------------------------------------------------------------------

const IM_FILE: &str = "im_peak_type.mzML";

#[test]
fn im_peak_type_reads_as_the_release_build_under_source_options() {
    let source = ReadOptions::source();
    assert_spectra(&read(IM_FILE, &source).spectra, "load", IM_FILE, "spectrum");
    // The source's float data arrays exist only when it fills data, so no
    // spectrum becomes IM_PROFILE (MzMLHandler.cpp:199-256).
    assert_spectra(
        &read_without_data(IM_FILE, &source).spectra,
        "load_no_data",
        IM_FILE,
        "spectrum",
    );
    assert_spectra(
        &transformed(IM_FILE, source),
        "transform",
        IM_FILE,
        "consumer_spectrum",
    );
}

#[test]
fn im_peak_type_keeps_the_centroid_frame_by_default() {
    // The spectra that carry MS:1003441 in the input: cases 1, 2 and 4 of
    // make_inputs.py. The Release build drops the term (its CV lacks it); the
    // native default reads it, and every other spectrum is as the oracle has it.
    let xml = std::fs::read_to_string(fixture(IM_FILE)).unwrap();
    let with_term: Vec<bool> = records(&xml, "spectrum")
        .iter()
        .map(|s| s.contains("MS:1003441"))
        .collect();
    assert_eq!(with_term, [false, true, true, false, true, false]);
    for experiment in [
        read(IM_FILE, &ReadOptions::default()).spectra,
        transformed(IM_FILE, ReadOptions::default()),
    ] {
        let expected = rows("load", IM_FILE, "spectrum");
        for (index, spectrum) in experiment.iter().enumerate() {
            let want = if with_term[index] {
                "im_centroided"
            } else {
                expected[index][2]
            };
            assert_eq!(
                im_peak_type_to_string(spectrum.im_peak_type),
                want,
                "spectrum {index}"
            );
        }
    }
    // fill_data=false never substitutes IM_PROFILE; the term is still read.
    let types: Vec<IonMobilityPeakType> = read_without_data(IM_FILE, &ReadOptions::default())
        .spectra
        .iter()
        .map(|s| s.im_peak_type)
        .collect();
    let u = IonMobilityPeakType::Unknown;
    let c = IonMobilityPeakType::Centroid;
    assert_eq!(types, [u, c, c, u, c, u]);
}

/// The in-memory spectra of the oracle's `store_constructed` case.
fn constructed() -> MSExperiment {
    let spectra = (0..3)
        .map(|k| {
            let mut spectrum = MSSpectrum::new();
            spectrum.rt = f64::from(k) + 1.0;
            spectrum.ms_level = 1;
            spectrum.native_id = format!("scan={}", k + 1);
            spectrum.peaks.push(Peak1D::new(100.0 + f64::from(k), 10.0));
            if k > 0 {
                spectrum
                    .float_data_arrays
                    .push(DataArray::new("Ion Mobility", vec![0.8_f32]));
            }
            spectrum.im_peak_type = if k < 2 {
                IonMobilityPeakType::Centroid
            } else {
                IonMobilityPeakType::Profile
            };
            spectrum
        })
        .collect();
    MSExperiment {
        spectra,
        ..MSExperiment::default()
    }
}

#[test]
fn centroided_frames_are_written_as_the_release_build_writes_them() {
    let experiment = constructed();
    let file = "constructed.mzML";
    let types: Vec<&str> = experiment
        .spectra
        .iter()
        .map(|s| im_peak_type_to_string(s.im_peak_type))
        .collect();
    let expected: Vec<&str> = rows("store_constructed", file, "spectrum")
        .iter()
        .map(|r| r[2])
        .collect();
    assert_eq!(types, expected);
    let xml = written(&experiment);
    let ours: Vec<bool> = records(&xml, "spectrum")
        .iter()
        .map(|s| s.contains(r#"accession="MS:1003441" name="ion mobility centroid frame""#))
        .collect();
    let theirs: Vec<bool> = oracle_written("store_constructed", file)
        .iter()
        .map(|lines| {
            lines
                .iter()
                .any(|l| l.contains(r#"accession="MS:1003441""#))
        })
        .collect();
    assert_eq!(ours, theirs);
    assert_eq!(ours, [true, true, false]);
    // The Release build reads its own output back as UNKNOWN / IM_PROFILE;
    // `source()` does the same, the native default keeps what was written.
    assert_spectra(
        &reread(&xml, &ReadOptions::source()).spectra,
        "reload_constructed",
        file,
        "spectrum",
    );
    let back: Vec<IonMobilityPeakType> = reread(&xml, &ReadOptions::default())
        .spectra
        .iter()
        .map(|s| s.im_peak_type)
        .collect();
    assert_eq!(
        back,
        experiment
            .spectra
            .iter()
            .map(|s| s.im_peak_type)
            .collect::<Vec<_>>()
    );
}

// --------------------------------------------------------------------------
// Selected-ion mobility on the spectrum
// --------------------------------------------------------------------------

const DRIFT_FILE: &str = "selected_ion_drift.mzML";

#[test]
fn selected_ion_mobility_reaches_the_spectrum_under_source_options() {
    let source = ReadOptions::source();
    let experiment = read(DRIFT_FILE, &source);
    assert_spectra(&experiment.spectra, "load", DRIFT_FILE, "spectrum");
    assert_spectra(
        &read_without_data(DRIFT_FILE, &source).spectra,
        "load_no_data",
        DRIFT_FILE,
        "spectrum",
    );
    assert_spectra(
        &transformed(DRIFT_FILE, source),
        "transform",
        DRIFT_FILE,
        "consumer_spectrum",
    );
    // A chromatogram precursor keeps its mobility on the precursor only.
    let chromatogram = &rows("load", DRIFT_FILE, "chromatogram")[0];
    assert_eq!(
        precursor_text(&experiment.chromatograms[0].precursor),
        oracle_precursor(chromatogram[2])
    );
    // Written and read back, the port reproduces what it read. The Release
    // build's own round trip differs at spectrum 4 only, where its precursor
    // writer turns a FAIMS voltage into a millisecond drift time (CPP-239).
    let back = reread(&written(&experiment), &source);
    assert_spectra(&back.spectra, "load", DRIFT_FILE, "spectrum");
    let load = rows("load", DRIFT_FILE, "spectrum");
    let reload = rows("reload", DRIFT_FILE, "spectrum");
    let differing: Vec<usize> = (0..load.len()).filter(|&i| load[i] != reload[i]).collect();
    assert_eq!(differing, [4]);
}

#[test]
fn selected_ion_mobility_stays_on_the_precursor_by_default() {
    let experiment = read(DRIFT_FILE, &ReadOptions::default());
    let oracle = rows("load", DRIFT_FILE, "spectrum");
    // Only the scan-level terms of the input set a spectrum mobility: 7.5 ms
    // on spectrum 0 and a -40 V FAIMS voltage on spectrum 3.
    let spectrum_mobility: Vec<(f64, DriftTimeUnit)> = experiment
        .spectra
        .iter()
        .map(|s| (s.drift_time, s.drift_time_unit))
        .collect();
    let none = (-1.0, DriftTimeUnit::None);
    assert_eq!(
        spectrum_mobility,
        [
            (7.5, DriftTimeUnit::Millisecond),
            none,
            none,
            (-40.0, DriftTimeUnit::FaimsCompensationVoltage),
            none,
            none,
            none,
        ]
    );
    // The precursors are the oracle's either way.
    for (spectrum, row) in experiment.spectra.iter().zip(&oracle) {
        assert_eq!(
            spectrum_fields(spectrum)[5..],
            oracle_spectrum_fields(row)[5..]
        );
    }
    // A file read this way writes back to one that reads back the same way.
    let back = reread(&written(&experiment), &ReadOptions::default());
    assert_eq!(back.spectra, experiment.spectra);
}

// --------------------------------------------------------------------------
// Activation: other dissociation methods
// --------------------------------------------------------------------------

const ACTIVATION_FILE: &str = "activation_children.mzML";

/// `(accession, value)` of every written cvParam of the pinned
/// `MS:1000044` descendants this route covers, sorted.
///
/// Sorted because the order of cvParams in an `activation` element carries no
/// meaning, and the two writers order them differently: the source walks the
/// precursor's meta keys in `MetaInfoRegistry` index order (the process-wide
/// order in which the key names were first registered), the port its
/// `BTreeMap` in name order.
fn route_terms(lines: &[&str]) -> Vec<(String, String)> {
    const ROUTE: [&str; 4] = ["MS:1003246", "MS:1003294", "MS:1003247", "MS:1003181"];
    let mut terms = Vec::new();
    for line in lines {
        for part in line.split("<cvParam ").skip(1) {
            let attribute = |name: &str| {
                part.split(&format!("{name}=\""))
                    .nth(1)
                    .map(|rest| rest.split('"').next().unwrap().to_owned())
            };
            let accession = attribute("accession").unwrap();
            if ROUTE.contains(&accession.as_str()) {
                terms.push((accession, attribute("value").unwrap_or_default()));
            }
        }
    }
    terms.sort();
    terms
}

#[test]
fn other_dissociation_methods_become_precursor_metadata() {
    for options in [ReadOptions::default(), ReadOptions::source()] {
        let experiment = read(ACTIVATION_FILE, &options);
        assert_spectra(&experiment.spectra, "load", ACTIVATION_FILE, "spectrum");
        // Native difference: the source's chromatogram activation branch has
        // no such route and drops the term with a warning (the oracle's "-");
        // this port keeps it, as it keeps the supplemental terms (CPP-219).
        assert_eq!(
            oracle_precursor(rows("load", ACTIVATION_FILE, "chromatogram")[0][2]),
            "-1e0|0|-|-"
        );
        assert_eq!(
            precursor_text(&experiment.chromatograms[0].precursor),
            "-1e0|0|-|ultraviolet photodissociation=0:"
        );
        // Written back as the CV terms the Release build writes, and read back
        // unchanged.
        let xml = written(&experiment);
        let ours: Vec<Vec<(String, String)>> = records(&xml, "spectrum")
            .iter()
            .map(|s| route_terms(&[s]))
            .collect();
        let theirs: Vec<Vec<(String, String)>> = oracle_written("store", ACTIVATION_FILE)
            .iter()
            .map(|lines| route_terms(lines))
            .collect();
        assert_eq!(ours, theirs);
        assert_eq!(ours.iter().map(Vec::len).sum::<usize>(), 5);
        let back = reread(&xml, &options);
        assert_eq!(back.spectra, experiment.spectra);
        assert_spectra(&back.spectra, "reload", ACTIVATION_FILE, "spectrum");
    }
}
