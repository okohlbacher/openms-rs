// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! `PeakPickerHiRes` and `SignalToNoiseEstimatorMedian` against the executed C++
//! and against the retained class-test outputs.
//!
//! * `data/peak_picking/oracle.tsv` is the output of the unmodified product-sdk
//!   libOpenMS (Debug, core `4fdec46`) driven over the cases of
//!   `data/peak_picking/cases.tsv` and `data/peak_picking/synthetic.tsv`
//!   (oracle-generated, tier 1 executed differential; driver, inputs and hashes
//!   in `../oracle/peak-picker-hires/`). Every centroid position, intensity,
//!   boundary, float-array value and signal-to-noise ratio is compared bit for
//!   bit. The unchanged driver re-run against the Linux x86-64 Release build of
//!   core `bc9cc12` on `ibminode06` printed the same bytes (sha256 `9eb8f249…`)
//!   and the same two parameter files (`../oracle/sne-completion/p1/`), so the
//!   fixture is also the Release build's output.
//! * `data/peak_picking/defaults.ini` and `noise_defaults.ini` are the product
//!   SDK's `ParamXMLFile::store` of the two classes' `getDefaults()`, and
//!   `WRITE_INI_OUT.ini` is the retained output of `TOPPWRITEINI_OVERWRITE`.
//! * The `*_sn1_out.mzML` and `*_sn4_out.mzML` files are the retained outputs the
//!   pinned `PeakPickerHiRes_test.cpp` compares with `TEST_REAL_SIMILAR`, and its
//!   literals are transcribed below.
//! * `data/peak_picking/chromatogram_time_oracle.tsv` is the output of the
//!   prebuilt C++ **Release** OpenMS (core `bc9cc12`, cli `c19e494`, topp
//!   `174b576`) over `data/peak_picking/chromatogram_time/*.mzML`, read and
//!   picked through `MzMLFile::load` and `PeakPickerHiRes::pickExperiment`
//!   (driver, cases and hashes in `../oracle/picked-chromatogram/`). It pins
//!   both sides of the source's minute conversion of a 32-bit chromatogram
//!   time array; see `chromatogram_time_unit_matches_the_cpp_reader_and_picker`.
//!
//! `docs/PEAK_PICKING_SUPPORT.md` maps every class-test section to these tests.
#![cfg(all(feature = "mzml", feature = "paramxml"))]

use openms::Error;
use openms::concept::parallel::Threads;
use openms::format::{mzml, paramxml};
use openms::kernel::{
    ChromatogramPeak, DataArray, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, SpectrumType,
};
use openms::metadata::{DataProcessing, ProcessingAction};
use openms::param::{Param, ParamValue};
use openms::processing::peak_picking::{
    BinIndexConversion, CENTROIDED_INPUT_MESSAGE, FwhmUnit, NoiseCompatibility, NoiseEstimates,
    NoiseHistogramRange, NoiseRangeParameters, PARALLEL_BATCH_RECORDS, PeakBoundary,
    PeakPickerHiRes, PickingCompatibility, SignalToNoiseEstimatorMedian,
};
use std::collections::BTreeMap;
use std::sync::Arc;

const CASES: &str = include_str!("data/peak_picking/cases.tsv");
const SYNTHETIC: &str = include_str!("data/peak_picking/synthetic.tsv");
const ORACLE: &str = include_str!("data/peak_picking/oracle.tsv");
const CHROMATOGRAM_TIME_ORACLE: &str =
    include_str!("data/peak_picking/chromatogram_time_oracle.tsv");

fn experiment(label: &str) -> MSExperiment {
    let bytes: &[u8] = match label {
        "orbitrap" => include_bytes!("data/peak_picking/PeakPickerHiRes_orbitrap.mzML"),
        "ftms" => include_bytes!("data/peak_picking/PeakPickerHiRes_ftms.unique_ids.mzML"),
        "selection" => include_bytes!("data/peak_picking/PeakPickerHiRes_spectrum_selection.mzML"),
        "simulation" => include_bytes!("data/peak_picking/PeakPickerHiRes_simulation.mzML"),
        "topp1" => include_bytes!("data/peak_picking/PeakPickerHiRes_input.mzML"),
        "topp2" => include_bytes!("data/peak_picking/PeakPickerHiRes_2_input.mzML"),
        "topp6" => include_bytes!("data/peak_picking/PeakPickerHiRes_6_input.mzML"),
        "orbitrap_sn1_out" => {
            include_bytes!("data/peak_picking/PeakPickerHiRes_orbitrap_sn1_out.mzML")
        }
        "orbitrap_sn4_out" => {
            include_bytes!("data/peak_picking/PeakPickerHiRes_orbitrap_sn4_out.mzML")
        }
        "ftms_sn1_out" => {
            include_bytes!("data/peak_picking/PeakPickerHiRes_ftms_sn1_out.unique_ids.mzML")
        }
        "ftms_sn4_out" => {
            include_bytes!("data/peak_picking/PeakPickerHiRes_ftms_sn4_out.unique_ids.mzML")
        }
        "noise" => {
            let spectrum = openms::format::dta::read(
                include_bytes!("data/peak_picking_noise_input.dta").as_slice(),
            )
            .unwrap();
            return MSExperiment {
                spectra: vec![spectrum],
                ..Default::default()
            };
        }
        other => panic!("unknown input label {other}"),
    };
    // MzMLFile::load with its default PeakFileOptions, which sort every spectrum
    // by m/z: the class-test file PeakPickerHiRes_spectrum_selection.mzML holds
    // three tandem spectra (scan=5537, 5541 and 5544) with one decreasing step
    // each, and the C++ picks them sorted.
    mzml::read_with_load_options(
        bytes,
        &mzml::LoadOptions::default(),
        &mzml::ReadOptions::default(),
    )
    .unwrap_or_else(|e| panic!("{label}: {e}"))
}

/// `TEST_REAL_SIMILAR` with the class-test defaults (`ClassTest.cpp:35-38`):
/// absolute difference at most `1e-5`, or ratio at most `1 + 1e-5`.
fn real_similar(a: f64, b: f64) -> bool {
    let absdiff = (a - b).abs();
    let small = absdiff <= 1e-5;
    if a == 0.0 || b == 0.0 {
        return small;
    }
    let mut ratio = a / b;
    if ratio < 0.0 {
        return small;
    }
    if ratio < 1.0 {
        ratio = 1.0 / ratio;
    }
    ratio <= 1.0 + 1e-5 || small
}

struct Synthetic {
    kind: String,
    x: Vec<f64>,
    y: Vec<f32>,
    arrays: Vec<(String, Vec<f32>)>,
    options: BTreeMap<String, String>,
}

fn floats(text: &str) -> Vec<f32> {
    // Parsed as double, then narrowed, exactly as the driver does.
    text.split(',')
        .map(|v| v.parse::<f64>().unwrap() as f32)
        .collect()
}

fn synthetic() -> BTreeMap<String, Synthetic> {
    let mut result = BTreeMap::new();
    for line in SYNTHETIC
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
    {
        let cols: Vec<_> = line.split('\t').collect();
        assert_eq!(cols.len(), 6, "{line}");
        let arrays = if cols[4] == "-" {
            Vec::new()
        } else {
            cols[4]
                .split('|')
                .map(|item| {
                    let (name, values) = item.split_once('=').unwrap();
                    (name.to_owned(), floats(values))
                })
                .collect()
        };
        let options = if cols[5] == "-" {
            BTreeMap::new()
        } else {
            cols[5]
                .split(';')
                .map(|item| {
                    let (k, v) = item.split_once('=').unwrap();
                    (k.to_owned(), v.to_owned())
                })
                .collect()
        };
        result.insert(
            cols[0].to_owned(),
            Synthetic {
                kind: cols[1].to_owned(),
                x: cols[2].split(',').map(|v| v.parse().unwrap()).collect(),
                y: floats(cols[3]),
                arrays,
                options,
            },
        );
    }
    result
}

fn to_spectrum(s: &Synthetic) -> MSSpectrum {
    let mut spectrum = MSSpectrum::from_peaks(
        s.x.iter()
            .zip(&s.y)
            .map(|(&x, &y)| Peak1D::new(x, y))
            .collect(),
    );
    for (name, values) in &s.arrays {
        spectrum
            .float_data_arrays
            .push(DataArray::new(name.clone(), values.clone()));
    }
    if let Some(kind) = s.options.get("type") {
        spectrum.spectrum_type = if kind == "centroid" {
            SpectrumType::Centroid
        } else {
            SpectrumType::Profile
        };
    }
    if let Some(level) = s.options.get("ms_level") {
        spectrum.ms_level = level.parse().unwrap();
    }
    if let Some(rt) = s.options.get("rt") {
        spectrum.rt = rt.parse().unwrap();
    }
    if s.options.contains_key("history") {
        let mut processing = DataProcessing::default();
        processing.actions.insert(ProcessingAction::PeakPicking);
        spectrum.data_processing.push(Arc::new(processing));
    }
    spectrum
}

fn to_chromatogram(s: &Synthetic) -> MSChromatogram {
    MSChromatogram {
        peaks: s
            .x
            .iter()
            .zip(&s.y)
            .map(|(&x, &y)| ChromatogramPeak::new(x, y))
            .collect(),
        ..Default::default()
    }
}

/// Apply `key=value;...` to a defaults tree, replacing whole entries as the
/// driver's `Param::setValue` does.
fn apply(mut param: Param, spec: &str) -> Param {
    if spec == "-" {
        return param;
    }
    for item in spec.split(';').filter(|s| !s.is_empty()) {
        let (key, value) = item.split_once('=').unwrap();
        let value = match param.value(key).unwrap() {
            ParamValue::Integer(_) => ParamValue::Integer(value.parse().unwrap()),
            ParamValue::Float(_) => ParamValue::Float(value.parse().unwrap()),
            ParamValue::String(_) => ParamValue::String(value.into()),
            ParamValue::IntegerList(_) => ParamValue::IntegerList(if value.is_empty() {
                Vec::new()
            } else {
                value.split(',').map(|v| v.parse().unwrap()).collect()
            }),
            other => panic!("unsupported parameter type {other:?}"),
        };
        param.set_value(key, value, "", &[]).unwrap();
    }
    param
}

#[derive(Debug, Default, PartialEq)]
struct Record {
    kind: String,
    index: usize,
    picked: bool,
    equal: bool,
    size: usize,
    type_code: i32,
    ms_level: u32,
    peaks: Vec<(u64, u32)>,
    bounds: Vec<(u64, u64)>,
    arrays: Vec<(String, Vec<u32>)>,
}

#[derive(Debug, PartialEq)]
enum Outcome {
    Error(String, String),
    Records(Vec<Record>),
    Noise(Vec<u64>, u64, u64),
}

struct Case {
    name: String,
    operation: String,
    input: String,
    parameters: String,
}

fn cases() -> Vec<Case> {
    CASES
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|line| {
            let cols: Vec<_> = line.split('\t').collect();
            assert_eq!(cols.len(), 4, "{line}");
            Case {
                name: cols[0].into(),
                operation: cols[1].into(),
                input: cols[2].into(),
                parameters: cols[3].into(),
            }
        })
        .collect()
}

fn bits32(text: &str) -> u32 {
    u32::from_str_radix(text, 16).unwrap()
}
fn bits64(text: &str) -> u64 {
    u64::from_str_radix(text, 16).unwrap()
}

fn oracle() -> BTreeMap<String, Outcome> {
    let mut result = BTreeMap::new();
    let mut name = String::new();
    let mut records: Vec<Record> = Vec::new();
    let mut error = None;
    let mut stn = Vec::new();
    let mut noise = None;
    for line in ORACLE.lines() {
        let cols: Vec<_> = line.split('\t').collect();
        match cols[0] {
            "case" => {
                name = cols[1].to_owned();
                records.clear();
                error = None;
                stn.clear();
                noise = None;
            }
            "error" => error = Some((cols[1].to_owned(), cols[2].to_owned())),
            "record" => records.push(Record {
                kind: cols[1].to_owned(),
                index: cols[2].parse().unwrap(),
                picked: cols[3] == "1",
                equal: cols[4] == "1",
                size: cols[5].parse().unwrap(),
                type_code: cols[6].parse().unwrap(),
                ms_level: cols[7].parse().unwrap(),
                ..Default::default()
            }),
            "p" => records
                .last_mut()
                .unwrap()
                .peaks
                .push((bits64(cols[1]), bits32(cols[2]))),
            "b" => records
                .last_mut()
                .unwrap()
                .bounds
                .push((bits64(cols[1]), bits64(cols[2]))),
            "a" => records.last_mut().unwrap().arrays.push((
                cols[2].to_owned(),
                cols[4..].iter().map(|v| bits32(v)).collect(),
            )),
            "n" => stn.push(bits64(cols[1])),
            "noise" => noise = Some((bits64(cols[1]), bits64(cols[2]))),
            "end" => {
                let outcome = if let Some((class, what)) = error.take() {
                    Outcome::Error(class, what)
                } else if let Some((sparse, rightmost)) = noise.take() {
                    Outcome::Noise(std::mem::take(&mut stn), sparse, rightmost)
                } else {
                    Outcome::Records(std::mem::take(&mut records))
                };
                result.insert(name.clone(), outcome);
            }
            other => panic!("unknown oracle line kind {other}"),
        }
    }
    result
}

fn type_code(kind: SpectrumType) -> i32 {
    match kind {
        SpectrumType::Unknown => 0,
        SpectrumType::Centroid => 1,
        SpectrumType::Profile => 2,
    }
}

fn spectrum_record(
    index: usize,
    picked: bool,
    equal: bool,
    s: &MSSpectrum,
    bounds: Option<&[PeakBoundary]>,
) -> Record {
    Record {
        kind: "S".into(),
        index,
        picked,
        equal,
        size: s.len(),
        type_code: type_code(s.spectrum_type),
        ms_level: s.ms_level,
        peaks: if picked {
            s.peaks
                .iter()
                .map(|p| (p.mz.to_bits(), p.intensity.to_bits()))
                .collect()
        } else {
            Vec::new()
        },
        bounds: bounds
            .unwrap_or_default()
            .iter()
            .map(|b| (b.min.to_bits(), b.max.to_bits()))
            .collect(),
        arrays: if picked {
            s.float_data_arrays
                .iter()
                .map(|a| (a.name.clone(), a.data.iter().map(|v| v.to_bits()).collect()))
                .collect()
        } else {
            Vec::new()
        },
    }
}

fn chromatogram_record(
    index: usize,
    equal: bool,
    c: &MSChromatogram,
    bounds: &[PeakBoundary],
) -> Record {
    Record {
        kind: "C".into(),
        index,
        picked: true,
        equal,
        size: c.len(),
        type_code: -1,
        ms_level: 0,
        peaks: c
            .peaks
            .iter()
            .map(|p| (p.rt.to_bits(), p.intensity.to_bits()))
            .collect(),
        bounds: bounds
            .iter()
            .map(|b| (b.min.to_bits(), b.max.to_bits()))
            .collect(),
        arrays: c
            .float_data_arrays
            .iter()
            .map(|a| (a.name.clone(), a.data.iter().map(|v| v.to_bits()).collect()))
            .collect(),
    }
}

fn error_outcome(error: Error) -> Outcome {
    match error {
        Error::InvalidValue(message) if message == CENTROIDED_INPUT_MESSAGE => {
            Outcome::Error("IllegalArgument".into(), message)
        }
        Error::InvalidValue(message) => Outcome::Error("InvalidValue".into(), message),
        other => Outcome::Error("other".into(), other.to_string()),
    }
}

struct Inputs {
    synthetic: BTreeMap<String, Synthetic>,
    files: BTreeMap<String, MSExperiment>,
}

impl Inputs {
    fn resolve(&mut self, input: &str) -> (MSExperiment, usize) {
        let (scheme, rest) = input.split_once(':').unwrap();
        let (label, index) = match rest.split_once('#') {
            Some((label, index)) => (label, index.parse().unwrap()),
            None => (rest, 0),
        };
        let experiment = match scheme {
            "file" | "chrom" | "dta" => self
                .files
                .entry(label.to_owned())
                .or_insert_with(|| experiment(label))
                .clone(),
            "syn" => {
                let mut result = MSExperiment::default();
                for part in label.split('+') {
                    let record = &self.synthetic[part];
                    if record.kind == "S" {
                        result.spectra.push(to_spectrum(record));
                    } else {
                        result.chromatograms.push(to_chromatogram(record));
                    }
                }
                result
            }
            other => panic!("unknown input scheme {other}"),
        };
        (experiment, index)
    }
}

fn run(case: &Case, inputs: &mut Inputs, compatibility: PickingCompatibility) -> Outcome {
    let (input, index) = inputs.resolve(&case.input);
    if case.operation == "noise" {
        let param = apply(
            SignalToNoiseEstimatorMedian::defaults().unwrap(),
            &case.parameters,
        );
        let estimator = SignalToNoiseEstimatorMedian::from_param(&param).unwrap();
        let spectrum = &input.spectra[index];
        let x: Vec<_> = spectrum.peaks.iter().map(|p| p.mz).collect();
        let y: Vec<_> = spectrum
            .peaks
            .iter()
            .map(|p| f64::from(p.intensity))
            .collect();
        return match estimator.estimate_with_compatibility(&x, &y, &compatibility) {
            Ok(NoiseEstimates {
                signal_to_noise,
                sparse_window_percent,
                histogram_rightmost_percent,
                ..
            }) => Outcome::Noise(
                signal_to_noise.iter().map(|v| v.to_bits()).collect(),
                sparse_window_percent.to_bits(),
                histogram_rightmost_percent.to_bits(),
            ),
            Err(e) => error_outcome(e),
        };
    }
    let param = apply(PeakPickerHiRes::defaults().unwrap(), &case.parameters);
    let mut picker = PeakPickerHiRes::from_param(&param).unwrap();
    picker.compatibility = compatibility;
    let result = match case.operation.as_str() {
        "spectrum" | "spectrum_nocheck" => picker
            .pick_spectrum_with_spacing(&input.spectra[index], case.operation == "spectrum")
            .map(|out| {
                let equal = out.spectrum == input.spectra[index];
                vec![spectrum_record(
                    index,
                    true,
                    equal,
                    &out.spectrum,
                    Some(&out.boundaries),
                )]
            }),
        "chromatogram" | "chromatogram_check" => picker
            .pick_chromatogram_with_spacing(
                &input.chromatograms[index],
                case.operation == "chromatogram_check",
            )
            .map(|out| {
                let equal = out.chromatogram == input.chromatograms[index];
                vec![chromatogram_record(
                    index,
                    equal,
                    &out.chromatogram,
                    &out.boundaries,
                )]
            }),
        "experiment" | "experiment_nocheck" => {
            picker.check_spectrum_type = case.operation == "experiment";
            picker.pick_experiment(&input).map(|out| {
                let mut records = Vec::new();
                for (i, spectrum) in out.experiment.spectra.iter().enumerate() {
                    let bounds = out.spectrum_boundaries[i].as_deref();
                    records.push(spectrum_record(
                        i,
                        bounds.is_some(),
                        *spectrum == input.spectra[i],
                        spectrum,
                        bounds,
                    ));
                }
                for (i, chromatogram) in out.experiment.chromatograms.iter().enumerate() {
                    records.push(chromatogram_record(
                        i,
                        *chromatogram == input.chromatograms[i],
                        chromatogram,
                        &out.chromatogram_boundaries[i],
                    ));
                }
                records
            })
        }
        other => panic!("unknown operation {other}"),
    };
    match result {
        Ok(records) => Outcome::Records(records),
        Err(e) => error_outcome(e),
    }
}

/// Compare an outcome with the oracle; error messages are compared by class
/// only, except the centroided-input text, which the port reproduces.
fn agrees(actual: &Outcome, expected: &Outcome) -> bool {
    match (actual, expected) {
        (Outcome::Error(a, am), Outcome::Error(e, em)) => {
            a == e && (e != "IllegalArgument" || am == em)
        }
        _ => actual == expected,
    }
}

#[test]
fn executed_source_cases_are_reproduced_bit_for_bit() {
    let expected = oracle();
    let cases = cases();
    assert_eq!(expected.len(), cases.len());
    let mut inputs = Inputs {
        synthetic: synthetic(),
        files: BTreeMap::new(),
    };
    let mut failures = Vec::new();
    let mut centroids = 0usize;
    for case in &cases {
        let oracle = &expected[&case.name];
        if let Outcome::Records(records) = oracle {
            centroids += records
                .iter()
                .filter(|r| r.picked)
                .map(|r| r.size)
                .sum::<usize>();
        }
        let source = run(case, &mut inputs, PickingCompatibility::source());
        if !agrees(&source, oracle) {
            failures.push(format!(
                "{} (source compatibility): {}",
                case.name,
                difference(&source, oracle)
            ));
        }
        // CPP-257: no case of this oracle bins a quotient outside the int range,
        // so the clamp-first conversion gives the same outcome as the x86-64 one.
        let clamped = run(
            case,
            &mut inputs,
            PickingCompatibility {
                noise: NoiseCompatibility {
                    bin_index: BinIndexConversion::ClampBeforeTruncation,
                    ..NoiseCompatibility::source()
                },
                ..PickingCompatibility::source()
            },
        );
        if !agrees(&clamped, oracle) {
            failures.push(format!(
                "{} (source compatibility, clamp-first bins): {}",
                case.name,
                difference(&clamped, oracle)
            ));
        }
        let strict = run(case, &mut inputs, PickingCompatibility::default());
        if case.name.starts_with("source_") {
            // The native default refuses what only the source accepts.
            if !matches!(strict, Outcome::Error(_, _)) {
                failures.push(format!("{}: the native default accepted it", case.name));
            }
        } else if !agrees(&strict, oracle) {
            failures.push(format!(
                "{} (native default): {}",
                case.name,
                difference(&strict, oracle)
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} cases differ:\n{}",
        failures.len(),
        cases.len(),
        failures.join("\n")
    );
    // Guard against a silently shrinking oracle.
    assert_eq!(centroids, 8795);
}

/// CPP-257 (`SignalToNoiseEstimatorMedian.h:297`, `:308`) can only matter for
/// a quotient `intensity / bin width` that is NaN or whose truncation leaves
/// the `int` range. Checks that no case of this oracle that runs the
/// estimator (the 12 `noise` cases and the 19 picker cases with
/// `signal_to_noise > 0`; `extra_noise_manual_invalid` and
/// `extra_auto_mode_manual_sn1` throw before binning)
/// has such a quotient in any record of its input, a superset of the records
/// the case estimates. So the agreement of the
/// clamp-first run above is not a coincidence of the medians: the path is not
/// exercised at all, and the `cpp257_*` cases of `tests/signal_to_noise.rs`
/// are its only evidence.
#[test]
fn no_executed_case_bins_a_quotient_outside_the_int_range() {
    let mut inputs = Inputs {
        synthetic: synthetic(),
        files: BTreeMap::new(),
    };
    let (mut estimating_cases, mut series, mut thrown) = (0, 0, Vec::new());
    for case in cases() {
        let estimator = if case.operation == "noise" {
            SignalToNoiseEstimatorMedian::from_param(&apply(
                SignalToNoiseEstimatorMedian::defaults().unwrap(),
                &case.parameters,
            ))
            .unwrap()
        } else {
            let picker = PeakPickerHiRes::from_param(&apply(
                PeakPickerHiRes::defaults().unwrap(),
                &case.parameters,
            ))
            .unwrap();
            if picker.signal_to_noise <= 0.0 {
                continue;
            }
            picker.noise_estimator
        };
        estimating_cases += 1;
        let (input, _) = inputs.resolve(&case.input);
        let records = input
            .spectra
            .iter()
            .map(|s| s.peaks.iter().map(|p| p.intensity).collect::<Vec<f32>>())
            .chain(
                input
                    .chromatograms
                    .iter()
                    .map(|c| c.peaks.iter().map(|p| p.intensity).collect()),
            );
        for intensities in records {
            // The histogram range does not depend on the positions.
            let x: Vec<f64> = (0..intensities.len()).map(|i| i as f64).collect();
            let y: Vec<f64> = intensities.iter().map(|v| f64::from(*v)).collect();
            let estimates = match estimator.estimate_with_compatibility(
                &x,
                &y,
                &PickingCompatibility::source(),
            ) {
                Ok(estimates) => estimates,
                // SignalToNoiseEstimatorMedian.h:236-244 throws before any bin.
                Err(Error::InvalidValue(_)) => {
                    thrown.push(case.name.clone());
                    continue;
                }
                Err(other) => panic!("{}: {other}", case.name),
            };
            // :258, `std::max(1.0, max_intensity_ / bin_count_)`.
            let width = f64::max(1.0, estimates.max_intensity / estimator.bin_count as f64);
            for value in &y {
                let quotient = value / width;
                assert!(
                    quotient > -2_147_483_649.0 && quotient < 2_147_483_648.0,
                    "{}: quotient {quotient}",
                    case.name
                );
            }
            series += 1;
        }
    }
    assert_eq!(estimating_cases, 31);
    thrown.dedup();
    assert_eq!(
        thrown,
        ["extra_auto_mode_manual_sn1", "extra_noise_manual_invalid"]
    );
    assert!(series > 31);
}

/// The first difference between two outcomes, for a readable failure.
fn difference(actual: &Outcome, expected: &Outcome) -> String {
    match (actual, expected) {
        (Outcome::Records(a), Outcome::Records(e)) => {
            if a.len() != e.len() {
                return format!("{} records, oracle {}", a.len(), e.len());
            }
            for (x, y) in a.iter().zip(e) {
                if x == y {
                    continue;
                }
                let head = |r: &Record| {
                    (
                        r.kind.clone(),
                        r.index,
                        r.picked,
                        r.equal,
                        r.size,
                        r.type_code,
                        r.ms_level,
                    )
                };
                if head(x) != head(y) {
                    return format!("record {:?}, oracle {:?}", head(x), head(y));
                }
                let differing: Vec<usize> = (0..x.peaks.len())
                    .filter(|&j| x.peaks[j] != y.peaks[j])
                    .collect();
                if let Some(&j) = differing.first() {
                    let (p, q) = (x.peaks[j], y.peaks[j]);
                    return format!(
                        "{} {} peaks {differing:?} differ; first {j}: ({}, {}) oracle ({}, {})",
                        x.kind,
                        x.index,
                        f64::from_bits(p.0),
                        f32::from_bits(p.1),
                        f64::from_bits(q.0),
                        f32::from_bits(q.1)
                    );
                }
                if x.bounds != y.bounds {
                    return format!(
                        "{} {} boundaries ({} vs {})",
                        x.kind,
                        x.index,
                        x.bounds.len(),
                        y.bounds.len()
                    );
                }
                return format!(
                    "{} {} arrays {:?} oracle {:?}",
                    x.kind, x.index, x.arrays, y.arrays
                );
            }
            "no difference".into()
        }
        (Outcome::Noise(a, ..), Outcome::Noise(e, ..)) => {
            for (i, (x, y)) in a.iter().zip(e).enumerate() {
                if x != y {
                    return format!(
                        "point {i}: {} oracle {}",
                        f64::from_bits(*x),
                        f64::from_bits(*y)
                    );
                }
            }
            format!("lengths or percentages differ: {actual:?} oracle {expected:?}")
        }
        _ => {
            let text = format!("{actual:?} oracle {expected:?}");
            text.chars().take(600).collect()
        }
    }
}

fn read_ini(text: &str) -> Param {
    paramxml::read(text.as_bytes()).unwrap()
}

#[test]
fn defaults_reproduce_the_source_parameter_files() {
    // ParamXMLFile::store of PeakPickerHiRes().getDefaults() and
    // SignalToNoiseEstimatorMedian<>().getDefaults() on the product SDK.
    assert_eq!(
        PeakPickerHiRes::defaults().unwrap(),
        read_ini(include_str!("data/peak_picking/defaults.ini"))
    );
    assert_eq!(
        SignalToNoiseEstimatorMedian::defaults().unwrap(),
        read_ini(include_str!("data/peak_picking/noise_defaults.ini"))
    );
    // TOPPWRITEINI_OVERWRITE: the algorithm node holds the current defaults,
    // except that the update keeps signal_to_noise = 1 from WRITE_INI_IN.ini.
    let written = read_ini(include_str!("data/peak_picking/WRITE_INI_OUT.ini"))
        .copy("PeakPickerHiRes:1:algorithm:", true)
        .unwrap();
    let expected = PeakPickerHiRes {
        signal_to_noise: 1.0,
        ..Default::default()
    }
    .to_param()
    .unwrap();
    assert_eq!(written, expected);
    assert_eq!(
        PeakPickerHiRes::defaults()
            .unwrap()
            .value("signal_to_noise")
            .unwrap(),
        &ParamValue::Float(0.0)
    );
}

#[test]
fn parameters_round_trip_through_the_typed_members() {
    let defaults = PeakPickerHiRes::defaults().unwrap();
    let picker = PeakPickerHiRes::from_param(&defaults).unwrap();
    assert_eq!(picker, PeakPickerHiRes::default());
    assert_eq!(picker.to_param().unwrap(), defaults);

    // A parameter tree holding values only (as after Param::setValue) maps to
    // the typed members and back to the complete tree.
    let values = apply(Param::new(), "-");
    assert_eq!(PeakPickerHiRes::from_param(&values).unwrap(), picker);
    let mut sparse = Param::new();
    for (key, value) in [
        ("signal_to_noise", ParamValue::Float(2.5)),
        ("spacing_difference_gap", ParamValue::Float(0.0)),
        ("spacing_difference", ParamValue::Float(3.0)),
        ("missing", ParamValue::Integer(0)),
        ("ms_levels", ParamValue::IntegerList(vec![1, 3])),
        ("report_FWHM", ParamValue::String("false".into())),
        ("report_FWHM_unit", ParamValue::String("absolute".into())),
        ("allow_missing_flank", ParamValue::String("true".into())),
        ("SignalToNoise:auto_mode", ParamValue::Integer(-1)),
        ("SignalToNoise:max_intensity", ParamValue::Integer(-1)),
        ("SignalToNoise:auto_max_percentile", ParamValue::Integer(7)),
        ("SignalToNoise:win_len", ParamValue::Float(40.0)),
        ("SignalToNoise:bin_count", ParamValue::Integer(12)),
        (
            "SignalToNoise:min_required_elements",
            ParamValue::Integer(3),
        ),
        (
            "SignalToNoise:noise_for_empty_window",
            ParamValue::Float(2.0),
        ),
        (
            "SignalToNoise:write_log_messages",
            ParamValue::String("false".into()),
        ),
    ] {
        sparse.set_value(key, value, "", &[]).unwrap();
    }
    let typed = PeakPickerHiRes::from_param(&sparse).unwrap();
    assert_eq!(typed.signal_to_noise, 2.5);
    assert_eq!(typed.spacing_difference_gap, 0.0);
    assert_eq!(typed.missing, 0);
    assert_eq!(typed.ms_levels, [1, 3]);
    assert_eq!(typed.report_fwhm, None);
    assert_eq!(typed.inactive_fwhm_unit, FwhmUnit::Absolute);
    assert!(typed.allow_missing_flank);
    assert_eq!(
        typed.noise_estimator.histogram_range,
        NoiseHistogramRange::Manual {
            max_intensity: -1.0
        }
    );
    assert_eq!(
        typed.noise_estimator.range_parameters,
        NoiseRangeParameters {
            max_intensity: -1,
            auto_max_stdev_factor: 3.0,
            auto_max_percentile: 7,
        }
    );
    assert!(!typed.noise_estimator.write_log_messages);
    let complete = typed.to_param().unwrap();
    assert!(
        complete
            .source_equal(&{
                let mut filled = sparse.clone();
                filled.set_defaults(&defaults, "", false).unwrap();
                filled
            })
            .unwrap()
    );
    assert_eq!(PeakPickerHiRes::from_param(&complete).unwrap(), typed);

    // Reporting FWHM in ppm; the inactive unit returns to its default.
    let fwhm = PeakPickerHiRes::from_param(&apply(defaults.clone(), "report_FWHM=true")).unwrap();
    assert_eq!(fwhm.report_fwhm, Some(FwhmUnit::Ppm));
    assert_eq!(
        PeakPickerHiRes::from_param(&fwhm.to_param().unwrap()).unwrap(),
        fwhm
    );

    // Percentile mode is accepted as a parameter.
    let percentile =
        PeakPickerHiRes::from_param(&apply(defaults.clone(), "SignalToNoise:auto_mode=1")).unwrap();
    assert_eq!(
        percentile.noise_estimator.histogram_range,
        NoiseHistogramRange::Percentile { percentile: 95.0 }
    );
    assert_eq!(
        PeakPickerHiRes::from_param(&percentile.to_param().unwrap()).unwrap(),
        percentile
    );
}

#[test]
fn parameter_failures_follow_the_source_contract() {
    let defaults = PeakPickerHiRes::defaults().unwrap();
    // Restrictions and value types are the source's InvalidParameter.
    for spec in [
        "report_FWHM_unit=THISVALUEISINVALID",
        "report_FWHM=yes",
        "signal_to_noise=-1.0",
        "missing=-1",
        "ms_levels=0",
        "SignalToNoise:auto_mode=2",
        "SignalToNoise:auto_max_stdev_factor=1000.0",
        "SignalToNoise:bin_count=2",
    ] {
        assert!(
            matches!(
                PeakPickerHiRes::from_param(&apply(defaults.clone(), spec)),
                Err(Error::InvalidValue(_))
            ),
            "{spec}"
        );
    }
    let mut wrong_type = defaults.clone();
    wrong_type
        .set_value("missing", ParamValue::Float(1.0), "", &[])
        .unwrap();
    assert!(PeakPickerHiRes::from_param(&wrong_type).is_err());
    // Unknown parameters are warnings only.
    let mut unknown = defaults.clone();
    unknown
        .set_value("invalidParamName", ParamValue::String("x".into()), "", &[])
        .unwrap();
    let (picker, warnings) = PeakPickerHiRes::from_param_with_warnings(&unknown).unwrap();
    assert_eq!(picker, PeakPickerHiRes::default());
    assert_eq!(warnings.len(), 1);
    assert!(warnings[0].contains("invalidParamName"), "{warnings:?}");
    // Values without a source parameter representation.
    for picker in [
        PeakPickerHiRes {
            missing: usize::MAX,
            ..Default::default()
        },
        PeakPickerHiRes {
            noise_estimator: SignalToNoiseEstimatorMedian {
                histogram_range: NoiseHistogramRange::Manual { max_intensity: 0.5 },
                ..Default::default()
            },
            ..Default::default()
        },
    ] {
        assert!(picker.to_param().is_err());
    }
}

#[test]
fn percentile_noise_mode_outside_its_domain_is_refused_only_when_estimation_runs() {
    let input = experiment("orbitrap");
    let percentile = SignalToNoiseEstimatorMedian {
        histogram_range: NoiseHistogramRange::Percentile { percentile: 95.0 },
        ..Default::default()
    };
    // Without estimation the source never initialises the estimator (oracle
    // case extra_auto_mode_percentile_sn0 picks normally).
    let quiet = PeakPickerHiRes {
        noise_estimator: percentile.clone(),
        ..Default::default()
    };
    assert!(quiet.pick_experiment(&input).is_ok());
    // With estimation on the orbitrap spectrum the product SDK crashes
    // (SIGBUS/SIGSEGV): its smallest intensity makes the pre-histogram index
    // leave [0, 99]. The port refuses there, in both profiles.
    let loud = PeakPickerHiRes {
        signal_to_noise: 1.0,
        ..quiet
    };
    for compatibility in [
        PickingCompatibility::default(),
        PickingCompatibility::source(),
    ] {
        let picker = PeakPickerHiRes {
            compatibility,
            ..loud.clone()
        };
        assert!(matches!(
            picker.pick_spectrum(&input.spectra[0]),
            Err(Error::Unsupported(message)) if message.contains("SignalToNoiseEstimatorMedian.h:216")
        ));
    }
    // Minimum 1: the bin size is 0.01 and intensity 2 has quotient 100.
    assert!(matches!(
        percentile.estimate(&[0.0, 1.0], &[1.0, 2.0]),
        Err(Error::Unsupported(_))
    ));
    // Fewer than five samples never reach the estimator, as in the source.
    let short = MSSpectrum::from_peaks(vec![Peak1D::new(1.0, 1.0); 1]);
    assert!(loud.pick_spectrum(&short).unwrap().spectrum.is_empty());
}

fn assert_similar_spectra(label: &str, actual: &MSExperiment, expected: &MSExperiment) {
    assert_eq!(actual.spectra.len(), expected.spectra.len(), "{label}");
    for (scan, (a, e)) in actual.spectra.iter().zip(&expected.spectra).enumerate() {
        assert_eq!(a.len(), e.len(), "{label} scan {scan}");
        for (peak, (p, q)) in a.peaks.iter().zip(&e.peaks).enumerate() {
            assert!(
                real_similar(p.mz, q.mz),
                "{label} scan {scan} peak {peak}: m/z {} vs {}",
                p.mz,
                q.mz
            );
            assert!(
                real_similar(f64::from(p.intensity), f64::from(q.intensity)),
                "{label} scan {scan} peak {peak}: intensity {} vs {}",
                p.intensity,
                q.intensity
            );
        }
    }
}

#[test]
fn class_test_pick_experiment_sections_match_the_retained_outputs() {
    // PeakPickerHiRes_test.cpp sections 7, 9 and 12, plus the whole FTMS S/N 1
    // experiment the test only checks for its first spectrum (section 10).
    let orbitrap = experiment("orbitrap");
    let ftms = experiment("ftms");
    for (label, input, sn, retained) in [
        ("orbitrap sn1", &orbitrap, 1.0, "orbitrap_sn1_out"),
        ("orbitrap sn4", &orbitrap, 4.0, "orbitrap_sn4_out"),
        ("ftms sn1", &ftms, 1.0, "ftms_sn1_out"),
        ("ftms sn4", &ftms, 4.0, "ftms_sn4_out"),
    ] {
        let picker = PeakPickerHiRes {
            signal_to_noise: sn,
            ..Default::default()
        };
        let out = picker.pick_experiment(input).unwrap();
        assert_similar_spectra(label, &out.experiment, &experiment(retained));
        // Sections 3, 4, 8, 10 and 11: single-spectrum pick of the first scan.
        let first = picker.pick_spectrum(&input.spectra[0]).unwrap();
        assert_eq!(first.spectrum, out.experiment.spectra[0]);
    }
}

#[test]
fn class_test_boundary_literals() {
    // Section 4: orbitrap S/N 1, first spectrum.
    let orbitrap = experiment("orbitrap");
    let sn1 = PeakPickerHiRes {
        signal_to_noise: 1.0,
        ..Default::default()
    }
    .pick_spectrum(&orbitrap.spectra[0])
    .unwrap();
    for (index, min, max) in [
        (25, 367.206604003906, 367.214569091797),
        (26, 369.042205810547, 369.051574707031),
    ] {
        assert!(real_similar(sn1.boundaries[index].min, min));
        assert!(real_similar(sn1.boundaries[index].max, max));
    }
    // Sections 14 and 15: signal_to_noise 0, missing 1, spacing_difference_gap 4.
    let picker = PeakPickerHiRes::from_param(&apply(
        PeakPickerHiRes::defaults().unwrap(),
        "signal_to_noise=0.0;missing=1;spacing_difference_gap=4.0",
    ))
    .unwrap();
    for (label, count, literals) in [
        (
            "simulation",
            167,
            vec![
                (146, 1141.57188829383, 1141.51216791402, 1141.63481354941),
                (148, 1142.57196823237, 1142.50968574851, 1142.6323313839),
                (158, 1178.08692219102, 1178.02013862689, 1178.14847787348),
                (159, 1178.58906411531, 1178.5249396635, 1178.6532789101),
            ],
        ),
        (
            "orbitrap",
            82,
            vec![
                (14, 355.070081088692, 355.064544677734, 355.078430175781),
                (37, 362.848715607077, 362.844085693359, 362.851928710938),
                (54, 370.210756298155, 370.205871582031, 370.215301513672),
                (55, 370.219596356153, 370.215301513672, 370.223358154297),
            ],
        ),
    ] {
        let out = picker.pick_experiment(&experiment(label)).unwrap();
        let spectrum = &out.experiment.spectra[0];
        let boundaries = out.spectrum_boundaries[0].as_ref().unwrap();
        assert_eq!(spectrum.len(), count, "{label}");
        for (index, mz, min, max) in literals {
            assert!(
                real_similar(spectrum.peaks[index].mz, mz),
                "{label} {index}"
            );
            assert!(real_similar(boundaries[index].min, min), "{label} {index}");
            assert!(real_similar(boundaries[index].max, max), "{label} {index}");
        }
    }
    // Section 15's comments: a boundary shared by two neighbouring peaks.
    let out = picker.pick_experiment(&orbitrap).unwrap();
    let boundaries = out.spectrum_boundaries[0].as_ref().unwrap();
    assert_eq!(boundaries[54].max, boundaries[55].min);
}

#[test]
fn class_test_spectrum_level_selection() {
    // Section 13. The loader sorts the three unsorted tandem spectra first.
    let input = experiment("selection");
    for (levels, picked) in [
        (vec![2], vec![2]),
        (vec![1], vec![1]),
        (vec![1, 2], vec![1, 2]),
    ] {
        let picker = PeakPickerHiRes {
            ms_levels: levels.clone(),
            ..Default::default()
        };
        let out = picker.pick_experiment(&input).unwrap();
        assert_eq!(out.experiment.spectra.len(), input.spectra.len());
        for (i, (a, b)) in out
            .experiment
            .spectra
            .iter()
            .zip(&input.spectra)
            .enumerate()
        {
            if picked.contains(&a.ms_level) {
                assert_ne!(a, b, "{levels:?} spectrum {i}");
            } else {
                assert_eq!(a, b, "{levels:?} spectrum {i}");
            }
        }
    }
}

#[test]
fn class_test_noise_estimator_section() {
    // SignalToNoiseEstimatorMedian_test.cpp [EXTRA] init: win_len 40,
    // noise_for_empty_window 2, min_required_elements 10, through the parameters.
    let input = experiment("noise");
    let historical = openms::format::dta::read(
        include_bytes!("data/peak_picking_noise_historical.dta").as_slice(),
    )
    .unwrap();
    let estimator = SignalToNoiseEstimatorMedian::from_param(&apply(
        SignalToNoiseEstimatorMedian::defaults().unwrap(),
        "win_len=40.0;noise_for_empty_window=2.0;min_required_elements=10",
    ))
    .unwrap();
    let x: Vec<_> = input.spectra[0].peaks.iter().map(|p| p.mz).collect();
    let y: Vec<_> = input.spectra[0]
        .peaks
        .iter()
        .map(|p| f64::from(p.intensity))
        .collect();
    let out = estimator.estimate(&x, &y).unwrap();
    assert_eq!(out.signal_to_noise.len(), historical.len());
    for (a, b) in out.signal_to_noise.iter().zip(&historical.peaks) {
        assert!(
            real_similar(f64::from(b.intensity), *a),
            "{a} vs {}",
            b.intensity
        );
    }
}

#[test]
fn mobility_output_description_depends_on_compatibility() {
    let mut input = MSSpectrum::from_peaks(
        [100., 100.01, 100.02, 100.03, 100.04]
            .into_iter()
            .zip([200., 250., 450., 250., 200.])
            .map(|(mz, intensity)| Peak1D::new(mz, intensity))
            .collect(),
    );
    let mut array = DataArray::new("Ion Mobility", vec![100., 150., 150., 150., 100.]);
    array.metadata.insert("unit".into(), "ms".into());
    input.float_data_arrays.push(array);
    let native = PeakPickerHiRes::default().pick_spectrum(&input).unwrap();
    assert_eq!(
        native.spectrum.float_data_arrays[0].metadata,
        input.float_data_arrays[0].metadata
    );
    let source = PeakPickerHiRes {
        compatibility: PickingCompatibility::source(),
        ..Default::default()
    }
    .pick_spectrum(&input)
    .unwrap();
    let only_name = &source.spectrum.float_data_arrays[0];
    assert_eq!(only_name.name, "Ion Mobility");
    assert!(!only_name.has_description_metadata());
    assert_eq!(only_name.data, native.spectrum.float_data_arrays[0].data);
}

/// A spectrum whose acquisition metadata costs the copy ledger `windows`
/// scan-window entries, each with one metadata entry.
///
/// The ledger charges a non-empty metadata tree for a sparsely occupied B-tree
/// node rather than for its few stored bytes, so a record is roughly twice as
/// expensive to the ledger as it is in memory -- a factor, not an escape:
/// exhausting the ledger's 256 MiB fixed part this way still costs the test
/// process on the order of a hundred mebibytes, which is why only the
/// single-record refusal below goes that far.
fn metadata_weight(windows: usize) -> MSSpectrum {
    let mut spectrum = MSSpectrum::default();
    spectrum.instrument_settings.scan_windows = (0..windows)
        .map(|i| {
            let mut window = openms::metadata::ScanWindow::new(100.0, 200.0).unwrap();
            window.metadata.insert("w".into(), (i as i64).into());
            window
        })
        .collect();
    spectrum
}

fn experiment_of(spectrum: &MSSpectrum, records: usize) -> MSExperiment {
    MSExperiment {
        spectra: vec![spectrum.clone(); records],
        ..MSExperiment::default()
    }
}

#[test]
fn an_overflowing_metadata_allowance_is_refused_rather_than_wrapping() {
    // `max_metadata_per_record` is multiplied by the record count and added to
    // the shared ledger's fixed part; both steps are checked, so an absurd
    // allowance refuses the pick instead of wrapping into a budget below the
    // fixed part. Neither entry point is exempt.
    //
    // That the budget itself grows by exactly the allowance per record -- the
    // property that lets a 40 856-spectrum run through where the fixed part
    // alone refused it -- is pinned by the synthetic ledger probe
    // `the_acquisition_ledger_follows_the_record_count` in
    // `src/processing/peak_picking.rs`, which reads the opened ledger instead
    // of building an experiment whose metadata exhausts 256 MiB of it. The
    // meter charges a record within about a factor of two of what that record
    // costs in memory, so crossing the fixed part end to end costs hundreds of
    // mebibytes; the real run is measured in `docs/PEAK_PICKING_SUPPORT.md`.
    let input = experiment_of(&metadata_weight(1), 8);
    for allowance in [usize::MAX, usize::MAX / 2, usize::MAX / 4] {
        let picker = PeakPickerHiRes {
            max_metadata_per_record: allowance,
            ..Default::default()
        };
        match picker.pick_experiment(&input).unwrap_err() {
            Error::InvalidValue(text) => assert!(text.contains("overflows"), "{text}"),
            other => panic!("{other:?}"),
        }
        let mut streamed = input.clone();
        assert!(picker.pick_experiment_in_place(&mut streamed).is_err());
        assert_eq!(streamed, input, "a refused ledger picks nothing");
    }
    // The same experiment is picked both with the default allowance and with
    // none, so the refusals above are the overflow and not the records.
    for allowance in [0, PeakPickerHiRes::default().max_metadata_per_record] {
        let picker = PeakPickerHiRes {
            max_metadata_per_record: allowance,
            ..Default::default()
        };
        picker.pick_experiment(&input).unwrap();
    }
}

#[test]
fn a_record_whose_metadata_outweighs_the_input_is_still_refused() {
    // The allowance is per record, so a single record cannot grow the budget
    // enough to pay for metadata that dwarfs the whole input. The experiment
    // moves the record instead of cloning it, because the search runs until the
    // metadata alone is heavy enough to exhaust the ledger.
    let picker = PeakPickerHiRes::default();
    let mut windows = 1 << 12;
    while windows <= 1 << 22 {
        let input = MSExperiment {
            spectra: vec![metadata_weight(windows)],
            ..MSExperiment::default()
        };
        if picker.pick_experiment(&input).is_err() {
            return;
        }
        windows *= 2;
    }
    unreachable!("a single record of unbounded acquisition metadata was admitted");
}

#[test]
fn in_place_picking_matches_the_borrowing_entry_point() {
    // The streaming entry point picks the same records by the same rules, so its
    // experiment and its four report vectors are those of `pick_experiment`.
    for label in [
        "orbitrap",
        "ftms",
        "selection",
        "simulation",
        "topp1",
        "topp6",
    ] {
        for levels in [vec![], vec![1], vec![1, 2]] {
            let picker = PeakPickerHiRes {
                ms_levels: levels.clone(),
                check_spectrum_type: false,
                compatibility: PickingCompatibility::source(),
                ..Default::default()
            };
            let input = experiment(label);
            let borrowed = picker.pick_experiment(&input).unwrap();
            let mut streamed = input.clone();
            let report = picker.pick_experiment_in_place(&mut streamed).unwrap();
            assert_eq!(streamed, borrowed.experiment, "{label} {levels:?}");
            assert_eq!(
                report.spectrum_boundaries, borrowed.spectrum_boundaries,
                "{label} {levels:?}"
            );
            assert_eq!(
                report.chromatogram_boundaries, borrowed.chromatogram_boundaries,
                "{label} {levels:?}"
            );
            assert_eq!(
                report.omitted_spectrum_arrays, borrowed.omitted_spectrum_arrays,
                "{label} {levels:?}"
            );
            assert_eq!(
                report.omitted_chromatogram_arrays, borrowed.omitted_chromatogram_arrays,
                "{label} {levels:?}"
            );
        }
    }
}

#[test]
fn in_place_picking_leaves_the_records_before_a_failure_picked() {
    // Unlike `pick_experiment`, the streaming entry point is not atomic; the
    // documentation says so and a caller that needs the input intact on failure
    // uses the borrowing one.
    let mut input = experiment("selection");
    let good = input.spectra[0].clone();
    input.spectra.insert(0, good.clone());
    let picker = PeakPickerHiRes {
        max_points: good.len(),
        compatibility: PickingCompatibility::source(),
        ..Default::default()
    };
    // A later record above the point ceiling fails after the first is picked.
    let mut big = good.clone();
    big.peaks.extend_from_slice(&good.peaks);
    input.spectra.push(big);
    let before = input.clone();
    assert!(picker.pick_experiment(&before.clone()).is_err());
    assert!(picker.pick_experiment_in_place(&mut input).is_err());
    assert_ne!(input.spectra[0], before.spectra[0]);
}

// --------------------------------------------------------------------------
// Chromatogram time-array unit handling (lane fx-chrom)
// --------------------------------------------------------------------------

/// One `case<TAB>stage<TAB>...` block of the chromatogram time-unit oracle.
///
/// Returns the `(retention time bits, intensity bits)` rows of `stage`, and the
/// declared point count that precedes them.
fn oracle_rows(case: &str, stage: &str) -> (usize, Vec<(u64, u32)>) {
    let mut count = None;
    let mut rows = Vec::new();
    for line in CHROMATOGRAM_TIME_ORACLE.lines() {
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 || fields[0] != case || fields[1] != stage {
            continue;
        }
        if fields[2] == "count" {
            count = Some(fields[3].parse().expect("count"));
            continue;
        }
        let bits = |text: &str| u64::from_str_radix(text.trim_start_matches("0x"), 16).expect(text);
        rows.push((bits(fields[3]), bits(fields[4]) as u32));
    }
    (count.expect("no count row"), rows)
}

/// The cases of `data/peak_picking/chromatogram_time/`, as the oracle names them.
const CHROMATOGRAM_TIME_CASES: [(&str, &[u8]); 7] = [
    (
        "min32",
        include_bytes!("data/peak_picking/chromatogram_time/min32.mzML"),
    ),
    (
        "min64",
        include_bytes!("data/peak_picking/chromatogram_time/min64.mzML"),
    ),
    (
        "sec32",
        include_bytes!("data/peak_picking/chromatogram_time/sec32.mzML"),
    ),
    (
        "sec64",
        include_bytes!("data/peak_picking/chromatogram_time/sec64.mzML"),
    ),
    (
        "min32_empty",
        include_bytes!("data/peak_picking/chromatogram_time/min32_empty.mzML"),
    ),
    (
        "min32_one",
        include_bytes!("data/peak_picking/chromatogram_time/min32_one.mzML"),
    ),
    (
        "min32_unsorted",
        include_bytes!("data/peak_picking/chromatogram_time/min32_unsorted.mzML"),
    ),
];

/// Read one case the way source `MzMLFile::load` does, whose default
/// `PeakFileOptions` sort every record by position.
fn chromatogram_time_case(bytes: &[u8], source_precision: bool) -> MSExperiment {
    let limits = if source_precision {
        mzml::ReadOptions::source()
    } else {
        mzml::ReadOptions::default()
    };
    mzml::read_with_load_options(bytes, &mzml::LoadOptions::default(), &limits).expect("load")
}

fn chromatogram_bits(chromatogram: &MSChromatogram) -> Vec<(u64, u32)> {
    chromatogram
        .peaks
        .iter()
        .map(|peak| (peak.rt.to_bits(), peak.intensity.to_bits()))
        .collect()
}

/// Every case read and picked exactly as the executed C++ reads and picks it.
///
/// The seven cases hold the same 40 real time and intensity values of the
/// `UK222.mzML` TIC chromatogram, differing only in the time array's encoding
/// and unit. `ReadOptions::source` reproduces the source's in-place minute
/// conversion on a `std::vector<float>`
/// (`MzMLHandlerHelper::decodeBase64Arrays`, `MzMLHandlerHelper.cpp:217-222`),
/// which keeps the converted seconds at 32-bit precision, so every case must
/// match the C++ bit for bit in that mode.
#[test]
fn chromatogram_time_unit_matches_the_cpp_reader_and_picker() {
    let picker = PeakPickerHiRes::default();
    for (case, bytes) in CHROMATOGRAM_TIME_CASES {
        let experiment = chromatogram_time_case(bytes, true);
        assert_eq!(experiment.chromatograms.len(), 1, "{case}");
        let input = &experiment.chromatograms[0];
        let (count, rows) = oracle_rows(case, "in");
        assert_eq!(input.peaks.len(), count, "{case} read point count");
        assert_eq!(chromatogram_bits(input), rows, "{case} read values");

        let (count, rows) = oracle_rows(case, "picked");
        let picked = picker.pick_chromatogram(input).expect("pick").chromatogram;
        assert_eq!(picked.peaks.len(), count, "{case} picked point count");
        assert_eq!(chromatogram_bits(&picked), rows, "{case} picked values");
    }
}

/// The library default keeps the precision the source drops, and nothing else.
///
/// The check is against the executed C++ itself, not against this port: reading
/// the same instrument times as a **64-bit** minute array (`min64`) is the one
/// shape where the source applies the multiplier in `double` and keeps the
/// result, so its `min64` rows are the C++'s own answer for the full-precision
/// seconds. The default read of the 32-bit `min32` must reproduce exactly that,
/// for the read values and for the picked chromatogram.
#[test]
fn chromatogram_time_default_keeps_the_precision_the_source_drops() {
    let picker = PeakPickerHiRes::default();
    let bytes = CHROMATOGRAM_TIME_CASES[0].1;
    assert_eq!(CHROMATOGRAM_TIME_CASES[0].0, "min32");
    let experiment = chromatogram_time_case(bytes, false);
    let input = &experiment.chromatograms[0];
    let (count, rows) = oracle_rows("min64", "in");
    assert_eq!(input.peaks.len(), count);
    assert_eq!(chromatogram_bits(input), rows, "default read values");

    let (count, rows) = oracle_rows("min64", "picked");
    let picked = picker.pick_chromatogram(input).expect("pick").chromatogram;
    assert_eq!(picked.peaks.len(), count);
    assert_eq!(chromatogram_bits(&picked), rows, "default picked values");

    // The two modes differ only by the source's narrowing: every default value
    // rounds to the value the source stores.
    let source = chromatogram_time_case(bytes, true);
    for (default, source) in input.peaks.iter().zip(&source.chromatograms[0].peaks) {
        assert_eq!(f64::from(default.rt as f32), source.rt);
        assert_eq!(default.intensity, source.intensity);
    }
}

/// The narrowing is confined to a 32-bit time array that carries a unit multiplier.
///
/// `min64` (64-bit, minutes), `sec32` (32-bit, seconds) and `sec64` (64-bit,
/// seconds) reach the same values in both modes, because the source narrows
/// only where it multiplies a `std::vector<float>` in place. `min32` is the one
/// case the switch moves.
#[test]
fn chromatogram_time_switch_only_moves_the_converted_32_bit_array() {
    for (case, bytes) in CHROMATOGRAM_TIME_CASES {
        let default = chromatogram_time_case(bytes, false);
        let source = chromatogram_time_case(bytes, true);
        let moved = default.chromatograms[0].peaks != source.chromatograms[0].peaks;
        let expected = matches!(case, "min32" | "min32_one" | "min32_unsorted");
        assert_eq!(moved, expected, "{case}");
    }
}

/// The empty and single-point chromatograms the source leaves alone.
///
/// `pick_` returns before the loop for fewer than five points
/// (`CENTROIDING/PeakPickerHiRes.cpp:147-150`), so both produce an empty picked
/// chromatogram while keeping the chromatogram's own settings.
#[test]
fn chromatogram_time_short_inputs_pick_nothing_and_keep_their_settings() {
    let picker = PeakPickerHiRes::default();
    for case in ["min32_empty", "min32_one"] {
        let bytes = CHROMATOGRAM_TIME_CASES
            .iter()
            .find(|(name, _)| *name == case)
            .expect("case")
            .1;
        let experiment = chromatogram_time_case(bytes, true);
        let input = &experiment.chromatograms[0];
        let picked = picker.pick_chromatogram(input).expect("pick").chromatogram;
        assert_eq!(oracle_rows(case, "picked").0, 0, "{case}");
        assert!(picked.peaks.is_empty(), "{case}");
        assert_eq!(picked.native_id, input.native_id, "{case}");
        assert_eq!(picked.metadata, input.metadata, "{case}");
    }
}

/// An unsorted time array is sorted on load, as source `MzMLHandler` sorts it.
///
/// The source reader's default `PeakFileOptions` sort every record by position,
/// so the picker never sees the file order; the C++ oracle's `min32_unsorted`
/// rows are the sorted ones, and picking them yields the same single centroid
/// as the first five points of `min32`. Reading without that sort leaves the
/// file order, which this port's picker refuses rather than treating as
/// undefined.
#[test]
fn chromatogram_time_unsorted_is_sorted_on_load_and_refused_without_it() {
    let bytes = CHROMATOGRAM_TIME_CASES[6].1;
    assert_eq!(CHROMATOGRAM_TIME_CASES[6].0, "min32_unsorted");
    let sorted = chromatogram_time_case(bytes, true);
    let (_, rows) = oracle_rows("min32_unsorted", "in");
    assert_eq!(chromatogram_bits(&sorted.chromatograms[0]), rows);
    let (_, min32) = oracle_rows("min32", "in");
    assert_eq!(rows, min32[..rows.len()].to_vec());

    let unsorted = mzml::read_with_options(bytes, &mzml::ReadOptions::source()).expect("read");
    let times: Vec<f64> = unsorted.chromatograms[0]
        .peaks
        .iter()
        .map(|peak| peak.rt)
        .collect();
    assert!(times.windows(2).any(|pair| pair[0] > pair[1]));
    assert!(matches!(
        PeakPickerHiRes::default().pick_chromatogram(&unsorted.chromatograms[0]),
        Err(Error::UnsortedData)
    ));
}

// --------------------------------------------------------------------------
// Parallel picking and the single-thread fixes (lane perf-picker)
// --------------------------------------------------------------------------

/// Every centroid and every boundary of a picked experiment, as raw bits.
///
/// `==` on `f64` would accept a value that drifted inside a tolerance, and
/// drift with the worker count is exactly the defect the determinism contract
/// excludes, so the comparison is on `to_bits`.
fn picked_bits(
    experiment: &MSExperiment,
    boundaries: &[Option<Vec<PeakBoundary>>],
    chromatogram_boundaries: &[Vec<PeakBoundary>],
) -> Vec<u64> {
    let mut bits = Vec::new();
    for spectrum in &experiment.spectra {
        bits.push(spectrum.peaks.len() as u64);
        for peak in &spectrum.peaks {
            bits.push(peak.mz.to_bits());
            bits.push(u64::from(peak.intensity.to_bits()));
        }
        for array in &spectrum.float_data_arrays {
            bits.push(array.data.len() as u64);
            bits.extend(array.data.iter().map(|v| u64::from(v.to_bits())));
        }
    }
    for chromatogram in &experiment.chromatograms {
        bits.push(chromatogram.peaks.len() as u64);
        for peak in &chromatogram.peaks {
            bits.push(peak.rt.to_bits());
            bits.push(u64::from(peak.intensity.to_bits()));
        }
    }
    for record in boundaries {
        match record {
            None => bits.push(u64::MAX),
            Some(list) => {
                bits.push(list.len() as u64);
                for boundary in list {
                    bits.push(boundary.min.to_bits());
                    bits.push(boundary.max.to_bits());
                }
            }
        }
    }
    for list in chromatogram_boundaries {
        bits.push(list.len() as u64);
        for boundary in list {
            bits.push(boundary.min.to_bits());
            bits.push(boundary.max.to_bits());
        }
    }
    bits
}

/// The worker counts every determinism assertion in this file runs at.
///
/// 32 exceeds the cores of an ordinary development machine on purpose: rayon
/// starts the workers whatever the machine has, so the schedule differs from
/// the one- and two-worker runs even where they cannot run truly concurrently.
const THREAD_COUNTS: [i64; 5] = [1, 2, 8, 32, 0];

/// Both parallel entry points are bit-identical to the serial ones at every
/// worker count.
///
/// This is the determinism contract of `src/concept/parallel.rs` applied to the
/// picker: centroids, float arrays, boundaries and the omitted-array reports
/// are fixed by the input and the parameters alone. Records are picked
/// independently and the pooled acquisition ledger is charged serially in input
/// order, so there is nothing for a schedule to change.
#[test]
fn picking_an_experiment_is_bit_identical_at_every_thread_count() {
    for label in ["orbitrap", "ftms", "selection", "simulation", "topp1"] {
        for levels in [vec![], vec![1], vec![1, 2]] {
            let picker = PeakPickerHiRes {
                ms_levels: levels.clone(),
                check_spectrum_type: false,
                compatibility: PickingCompatibility::source(),
                report_fwhm: Some(FwhmUnit::Absolute),
                ..Default::default()
            };
            let input = experiment(label);
            let serial = picker.pick_experiment(&input).expect("serial pick");
            let expected = picked_bits(
                &serial.experiment,
                &serial.spectrum_boundaries,
                &serial.chromatogram_boundaries,
            );
            for count in THREAD_COUNTS {
                let threads = Threads::from_cli(count);
                let owned = picker
                    .pick_experiment_with_threads(&input, threads)
                    .expect("parallel pick");
                assert_eq!(
                    picked_bits(
                        &owned.experiment,
                        &owned.spectrum_boundaries,
                        &owned.chromatogram_boundaries
                    ),
                    expected,
                    "{label} {levels:?} at {count} threads"
                );
                assert_eq!(owned.experiment, serial.experiment, "{label} {count}");
                assert_eq!(
                    owned.omitted_spectrum_arrays, serial.omitted_spectrum_arrays,
                    "{label} {count}"
                );
                assert_eq!(
                    owned.omitted_chromatogram_arrays, serial.omitted_chromatogram_arrays,
                    "{label} {count}"
                );

                let mut streamed = input.clone();
                let report = picker
                    .pick_experiment_in_place_with_threads(&mut streamed, threads)
                    .expect("parallel in-place pick");
                assert_eq!(
                    picked_bits(
                        &streamed,
                        &report.spectrum_boundaries,
                        &report.chromatogram_boundaries
                    ),
                    expected,
                    "{label} {levels:?} in place at {count} threads"
                );
                assert_eq!(streamed, serial.experiment, "{label} in place {count}");
                assert_eq!(
                    report.spectrum_boundaries, serial.spectrum_boundaries,
                    "{label} in place {count}"
                );
                assert_eq!(
                    report.omitted_spectrum_arrays, serial.omitted_spectrum_arrays,
                    "{label} in place {count}"
                );
            }
        }
    }
}

/// More records than one parallel batch holds, so the batch loop runs several
/// times and its boundaries fall in the middle of the run.
///
/// The batch bounds are a memory ceiling, not a tuning knob that may move a
/// value: a run split into three batches has to give the same answer as the
/// serial loop, and the per-record reports have to stay indexed by input
/// position across the joins.
#[test]
fn a_run_of_several_parallel_batches_matches_the_serial_loop() {
    let source = experiment("selection");
    let mut input = MSExperiment {
        settings: source.settings.clone(),
        ..MSExperiment::default()
    };
    // Enough records to cross `PARALLEL_BATCH_RECORDS` twice.
    while input.spectra.len() < 2 * PARALLEL_BATCH_RECORDS + 13 {
        for spectrum in &source.spectra {
            input.spectra.push(spectrum.clone());
        }
    }
    let picker = PeakPickerHiRes {
        check_spectrum_type: false,
        compatibility: PickingCompatibility::source(),
        ..Default::default()
    };
    let serial = picker.pick_experiment(&input).expect("serial pick");
    let expected = picked_bits(
        &serial.experiment,
        &serial.spectrum_boundaries,
        &serial.chromatogram_boundaries,
    );
    for count in THREAD_COUNTS {
        let mut streamed = input.clone();
        let report = picker
            .pick_experiment_in_place_with_threads(&mut streamed, Threads::from_cli(count))
            .expect("parallel pick");
        assert_eq!(
            picked_bits(
                &streamed,
                &report.spectrum_boundaries,
                &report.chromatogram_boundaries
            ),
            expected,
            "{count} threads over {} records",
            input.spectra.len()
        );
    }
}

/// The in-place failure detail survives the parallel pass: exactly the records
/// before the first failing one are replaced, whatever the worker count.
///
/// A parallel write-back would leave a schedule-dependent subset centroided.
/// The batch is prepared in parallel and committed in input order instead, so
/// the prefix is the serial one and the error is the first error in input
/// order, not whichever worker failed first.
#[test]
fn the_in_place_failure_prefix_is_the_same_at_every_thread_count() {
    let base = experiment("selection");
    let good = base.spectra[0].clone();
    let picker = PeakPickerHiRes {
        max_points: good.len(),
        compatibility: PickingCompatibility::source(),
        ..Default::default()
    };
    // Two pickable records, then one above the point ceiling, then two more.
    let mut big = good.clone();
    big.peaks.extend_from_slice(&good.peaks);
    let mut input = MSExperiment {
        settings: base.settings.clone(),
        ..MSExperiment::default()
    };
    input.spectra.extend([good.clone(), good.clone()]);
    input.spectra.push(big);
    input.spectra.extend([good.clone(), good.clone()]);

    let mut serial = input.clone();
    assert!(picker.pick_experiment_in_place(&mut serial).is_err());
    for count in THREAD_COUNTS {
        let mut streamed = input.clone();
        let error = picker
            .pick_experiment_in_place_with_threads(&mut streamed, Threads::from_cli(count))
            .expect_err("the oversized record must be refused");
        assert!(
            matches!(&error, Error::InvalidValue(text) if text.contains("point limit")),
            "{count} threads: {error:?}"
        );
        assert_eq!(
            streamed, serial,
            "{count} threads changed the failure prefix"
        );
        // The prefix is exactly the two records before the failure.
        assert_ne!(streamed.spectra[0], input.spectra[0], "{count}");
        assert_ne!(streamed.spectra[1], input.spectra[1], "{count}");
        assert_eq!(streamed.spectra[2], input.spectra[2], "{count}");
        assert_eq!(streamed.spectra[3], input.spectra[3], "{count}");
        assert_eq!(streamed.spectra[4], input.spectra[4], "{count}");
    }
}

/// The picked record carries every metadata field of the input, which is what
/// cloning the whole record and overwriting its samples used to guarantee.
///
/// The output is now built field by field from the input's metadata (source
/// `copySpectrumMeta`), so that picking no longer copies the profile samples and
/// the annotation arrays of every record only to drop them. This pins the
/// equivalence against the recipe it replaced: clone the input, put the picked
/// samples in, mark it centroided and clear the arrays.
#[test]
fn the_picked_record_carries_every_metadata_field_the_input_had() {
    let mut input = experiment("orbitrap").spectra[0].clone();
    // Fields a fixture may leave at their defaults, so nothing is pinned by
    // accident of both sides being empty.
    input.rt = 1234.5;
    input.ms_level = 1;
    input.native_id = "controllerType=0 controllerNumber=1 scan=42".to_owned();
    input.name = "a named spectrum".to_owned();
    input.spectrum_type = SpectrumType::Profile;
    input.drift_time = 17.5;
    input.source_file.name = "source.raw".to_owned();
    input.metadata.insert("kept".into(), "value".into());
    let record = Arc::new(DataProcessing::default());
    input.data_processing = vec![Arc::clone(&record)];
    input
        .integer_data_arrays
        .push(DataArray::new("counts", vec![1_i32; input.len()]));
    input.string_data_arrays.push(DataArray::new(
        "labels",
        vec![String::from("x"); input.len()],
    ));

    let picker = PeakPickerHiRes {
        compatibility: PickingCompatibility::source(),
        ..Default::default()
    };
    let picked = picker.pick_spectrum(&input).expect("pick");
    assert!(!picked.spectrum.peaks.is_empty());

    let mut expected = input.clone();
    expected.peaks = picked.spectrum.peaks.clone();
    expected.spectrum_type = SpectrumType::Centroid;
    expected.float_data_arrays = picked.spectrum.float_data_arrays.clone();
    expected.integer_data_arrays.clear();
    expected.string_data_arrays.clear();
    assert_eq!(picked.spectrum, expected);

    // Processing handles stay shared rather than being deep-copied, as they did
    // when the record was cloned.
    assert!(Arc::ptr_eq(
        &input.data_processing[0],
        &picked.spectrum.data_processing[0]
    ));
    // The handle itself, the input's, the picked record's and the reference
    // `expected` built just above.
    assert_eq!(Arc::strong_count(&record), 4);

    // The same for a chromatogram.
    let mut chromatogram = experiment("selection")
        .chromatograms
        .pop()
        .unwrap_or_default();
    chromatogram.peaks = input
        .peaks
        .iter()
        .map(|peak| ChromatogramPeak::new(peak.mz, peak.intensity))
        .collect();
    chromatogram.native_id = "chromatogram=1".to_owned();
    chromatogram.name = "a named chromatogram".to_owned();
    chromatogram.metadata.insert("kept".into(), "value".into());
    chromatogram.data_processing = vec![Arc::clone(&record)];
    let picked = picker
        .pick_chromatogram(&chromatogram)
        .expect("pick chromatogram");
    let mut expected = chromatogram.clone();
    expected.peaks = picked.chromatogram.peaks.clone();
    expected.float_data_arrays.clear();
    expected.integer_data_arrays.clear();
    expected.string_data_arrays.clear();
    assert_eq!(picked.chromatogram, expected);
}

/// An invalid record is still refused by the experiment entry points, with the
/// same error, after the second per-record validation pass was removed.
///
/// `start_experiment` validates the whole experiment through
/// `MSExperiment::validate`, which visits every spectrum and every one of its
/// peaks. The selected records used to be validated a second time on the way
/// into `pick_`, over the same, unchanged samples. Dropping that pass has to
/// leave both refusals intact: the one for a record the picker would have
/// handled, and the one for a record the selection rule would only have copied,
/// which the per-record pass never reached in the first place.
#[test]
fn an_invalid_record_is_refused_by_the_experiment_entry_points() {
    let picker = PeakPickerHiRes {
        ms_levels: vec![1],
        check_spectrum_type: false,
        compatibility: PickingCompatibility::source(),
        ..Default::default()
    };
    let base = experiment("selection");
    let good = base.spectra[0].clone();
    assert_eq!(good.ms_level, 1, "the fixture record must be picked");
    let mut broken = good.clone();
    broken.peaks[3] = Peak1D::new(f64::NAN, 1.0);

    for (label, level) in [("picked", 1_u32), ("copied", 2)] {
        let mut copied = broken.clone();
        copied.ms_level = level;
        let input = MSExperiment {
            spectra: vec![good.clone(), copied],
            settings: base.settings.clone(),
            ..MSExperiment::default()
        };
        for count in THREAD_COUNTS {
            let threads = Threads::from_cli(count);
            let error = picker
                .pick_experiment_with_threads(&input, threads)
                .expect_err("an invalid record must be refused");
            assert!(
                matches!(&error, Error::InvalidValue(text) if text.contains("peak m/z")),
                "{label} at {count} threads: {error:?}"
            );
            // Refused before anything is written back, so the input is intact.
            // Compared record by record rather than with `==` on the whole
            // experiment, because the invalid record's `NaN` position is not
            // equal to itself.
            let mut streamed = input.clone();
            assert!(
                picker
                    .pick_experiment_in_place_with_threads(&mut streamed, threads)
                    .is_err(),
                "{label} at {count} threads"
            );
            assert_eq!(streamed.spectra[0], input.spectra[0], "{label} {count}");
            assert_eq!(
                streamed.spectra[1].peaks.len(),
                input.spectra[1].peaks.len(),
                "{label} {count}"
            );
            assert!(
                streamed.spectra[1].peaks[3].mz.is_nan(),
                "{label} at {count} threads replaced the invalid record"
            );
        }
    }

    // A refusal only `pick_` can make - here unsorted samples under the native
    // default - still comes from inside the record loop.
    let mut unsorted = good.clone();
    unsorted.peaks.reverse();
    let input = MSExperiment {
        spectra: vec![unsorted],
        settings: base.settings.clone(),
        ..MSExperiment::default()
    };
    let native = PeakPickerHiRes {
        ms_levels: vec![1],
        check_spectrum_type: false,
        ..Default::default()
    };
    assert!(matches!(
        native.pick_experiment_with_threads(&input, Threads::from_cli(8)),
        Err(Error::UnsortedData)
    ));
}
