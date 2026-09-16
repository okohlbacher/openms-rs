// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! `SignalToNoiseEstimatorMedian`, the `SignalToNoiseEstimator` base and
//! `estimateNoiseFromRandomScans` against the Linux x86-64 Release build.
//!
//! `data/signal_to_noise/oracle.tsv` and `streams.tsv` are the output of the
//! unmodified Release build `openms4-release-bc9cc12-c19e494-174b576`, driven
//! on `ibminode06` over the cases of `data/signal_to_noise/cases.tsv` and
//! `synthetic.tsv`, every case in a fresh process and twice, byte-identically
//! (oracle-generated, tier 1 executed differential; driver, inputs, disassembly
//! and hashes in `../oracle/sne-completion/`). Every ratio, percentage,
//! histogram upper end and random-scan result is compared as an IEEE-754 bit
//! pattern; the warnings are compared byte for byte with the C++ standard
//! error, and the progress report with its standard output after the timings
//! are masked. `docs/SIGNAL_TO_NOISE_SUPPORT.md` lists what each case pins.

use openms::concept::log_stream::{LogSink, LogStream};
use openms::concept::progress_logger::{
    CommandProgressLogger, ProgressClock, ProgressLogType, ProgressLogger, ProgressNesting,
    ProgressTime,
};
use openms::kernel::{ChromatogramPeak, MSChromatogram, MSExperiment, MSSpectrum, Peak1D};
use openms::param::{Param, ParamValue};
use openms::processing::mean_noise::SignalToNoiseEstimatorMeanIterative;
use openms::processing::noise_estimation::{
    GaussianEstimate, MinstdRand0, RandomScanNoise, SignalToNoiseEstimator,
    estimate_noise_from_random_scans,
};
use openms::processing::peak_picking::{
    BinIndexConversion, NOISE_PROGRESS_LABEL, NoiseCompatibility, NoiseEstimates,
    NoiseHistogramRange, PeakPickerHiRes, PickingCompatibility, SignalToNoiseEstimatorMedian,
};
use openms::{Error, Result};
use std::collections::BTreeMap;
use std::io::{self, Write};
use std::sync::{Arc, Mutex};

const CASES: &str = include_str!("data/signal_to_noise/cases.tsv");
const SYNTHETIC: &str = include_str!("data/signal_to_noise/synthetic.tsv");
const ORACLE: &str = include_str!("data/signal_to_noise/oracle.tsv");
const STREAMS: &str = include_str!("data/signal_to_noise/streams.tsv");

// ---------------------------------------------------------------------------
// Parsing, as the driver parses (strtod; `static_cast<float>` for intensities).

const NAN: u64 = 0x7ff8_0000_0000_0000;
const NEGATIVE_NAN: u64 = 0xfff8_0000_0000_0000;

fn number(text: &str) -> f64 {
    match text {
        "nan" => f64::from_bits(NAN),
        "-nan" => f64::from_bits(NEGATIVE_NAN),
        _ => text
            .parse()
            .unwrap_or_else(|e| panic!("bad number {text}: {e}")),
    }
}

/// `static_cast<float>`: `cvtsd2ss` keeps a NaN's sign.
fn narrow(value: f64) -> f32 {
    if value.is_nan() {
        f32::from_bits(if value.is_sign_negative() {
            0xffc0_0000
        } else {
            0x7fc0_0000
        })
    } else {
        value as f32
    }
}

fn list(text: &str) -> Vec<f64> {
    if text == "-" {
        return Vec::new();
    }
    let mut out = Vec::new();
    for token in text.split(',') {
        if let Some((value, count)) = token.split_once('*') {
            let value = number(value);
            out.extend(std::iter::repeat_n(value, count.parse().unwrap()));
        } else if let Some((start, rest)) = token.split_once(':') {
            let (step, count) = rest.split_once(':').unwrap();
            let (start, step) = (number(start), number(step));
            for i in 0..count.parse::<usize>().unwrap() {
                out.push(start + i as f64 * step);
            }
        } else {
            out.push(number(token));
        }
    }
    out
}

struct Synthetic {
    kind: String,
    ms_level: u32,
    x: Vec<f64>,
    y: Vec<f32>,
}

fn synthetic() -> BTreeMap<String, Synthetic> {
    SYNTHETIC
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .map(|line| {
            let c: Vec<&str> = line.split('\t').collect();
            let x = list(c[3]);
            let y: Vec<f32> = list(c[4]).into_iter().map(narrow).collect();
            assert_eq!(x.len(), y.len(), "{}", c[0]);
            (
                c[0].to_owned(),
                Synthetic {
                    kind: c[1].to_owned(),
                    ms_level: c[2].parse().unwrap(),
                    x,
                    y,
                },
            )
        })
        .collect()
}

fn keys(spec: &str) -> BTreeMap<&str, &str> {
    if spec == "-" {
        return BTreeMap::new();
    }
    spec.split(';')
        .filter(|s| !s.is_empty())
        .map(|item| item.split_once('=').unwrap())
        .collect()
}

/// The driver's `applyParameters`: `Param::setValue` with the default's type.
fn apply(mut param: Param, spec: &str) -> Param {
    for (key, value) in keys(spec) {
        let value = match param.value(key).unwrap() {
            ParamValue::Integer(_) => ParamValue::Integer(value.parse().unwrap()),
            ParamValue::Float(_) => ParamValue::Float(number(value)),
            ParamValue::String(_) => ParamValue::String(value.into()),
            other => panic!("unsupported parameter type {other:?}"),
        };
        param.set_value(key, value, "", &[]).unwrap();
    }
    param
}

#[derive(Clone, Debug, Default)]
struct Case {
    name: String,
    op: String,
    input: String,
    params: String,
}

fn cases() -> Vec<Case> {
    CASES
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            assert_eq!(c.len(), 4, "{l}");
            Case {
                name: c[0].into(),
                op: c[1].into(),
                input: c[2].into(),
                params: c[3].into(),
            }
        })
        .collect()
}

#[derive(Debug, Default)]
struct Expected {
    error: Option<(String, String)>,
    max_intensity: Option<u64>,
    count: Option<usize>,
    stn: Option<(usize, u64, Vec<u64>)>,
    stn_at: Vec<(usize, u64)>,
    percentages: Option<(u64, u64)>,
    twice: Option<bool>,
    random: Option<u32>,
    time_calls: u64,
    raw: Vec<u64>,
    draws: Vec<(u64, u32)>,
    peaks: Vec<(u64, u32)>,
    bounds: Vec<(u64, u64)>,
    exit: i32,
    nth_vectors: usize,
    stdout: Vec<u8>,
    stderr: Vec<u8>,
}

fn hex64(text: &str) -> u64 {
    u64::from_str_radix(text, 16).unwrap()
}

fn hex32(text: &str) -> u32 {
    u32::from_str_radix(text, 16).unwrap()
}

fn base64(text: &str) -> Vec<u8> {
    const ALPHABET: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = Vec::new();
    let (mut buffer, mut bits) = (0u32, 0u32);
    for byte in text.bytes().filter(|b| *b != b'=') {
        let value = ALPHABET.iter().position(|c| *c == byte).unwrap() as u32;
        buffer = (buffer << 6) | value;
        bits += 6;
        if bits >= 8 {
            bits -= 8;
            out.push((buffer >> bits) as u8);
            buffer &= (1 << bits) - 1;
        }
    }
    out
}

fn oracle() -> BTreeMap<String, Expected> {
    let mut out: BTreeMap<String, Expected> = BTreeMap::new();
    for line in ORACLE.lines() {
        let c: Vec<&str> = line.split('\t').collect();
        let e = out.entry(c[1].to_owned()).or_default();
        match c[0] {
            "case" | "end" => {}
            "error" => e.error = Some((c[2].into(), c[3].into())),
            "maxint" => e.max_intensity = Some(hex64(c[2])),
            "count" => e.count = Some(c[2].parse().unwrap()),
            "stn" => {
                let bits = c[4..].iter().map(|v| hex64(v)).collect();
                e.stn = Some((c[2].parse().unwrap(), hex64(c[3]), bits));
            }
            "stnat" => e.stn_at.push((c[2].parse().unwrap(), hex64(c[3]))),
            "noise" => e.percentages = Some((hex64(c[2]), hex64(c[3]))),
            "twice" => e.twice = Some(c[2] == "1"),
            "random" => e.random = Some(hex32(c[2])),
            "timecalls" => e.time_calls = c[2].parse().unwrap(),
            "raw" => e.raw.push(c[3].parse().unwrap()),
            "draw" => e.draws.push((hex64(c[3]), c[4].parse().unwrap())),
            "peak" => e.peaks.push((hex64(c[3]), hex32(c[4]))),
            "bound" => e.bounds.push((hex64(c[3]), hex64(c[4]))),
            "exit" => e.exit = c[2].parse().unwrap(),
            // The whole permutation of std::nth_element; compared by the unit
            // test of the private libstdc++ port in noise_estimation.rs.
            "nthout" => e.nth_vectors += 1,
            other => panic!("unknown oracle record {other}"),
        }
    }
    for line in STREAMS.lines() {
        let c: Vec<&str> = line.split('\t').collect();
        let bytes = base64(c.get(2).copied().unwrap_or(""));
        let e = out.get_mut(c[0]).unwrap();
        match c[1] {
            "stdout" => e.stdout = bytes,
            "stderr" => e.stderr = bytes,
            other => panic!("unknown stream {other}"),
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Inputs.

fn to_spectrum(record: &Synthetic) -> MSSpectrum {
    let mut spectrum = MSSpectrum::from_peaks(
        record
            .x
            .iter()
            .zip(&record.y)
            .map(|(&x, &y)| Peak1D::new(x, y))
            .collect(),
    );
    spectrum.ms_level = record.ms_level;
    spectrum
}

fn to_chromatogram(record: &Synthetic) -> MSChromatogram {
    MSChromatogram {
        peaks: record
            .x
            .iter()
            .zip(&record.y)
            .map(|(&x, &y)| ChromatogramPeak::new(x, y))
            .collect(),
        ..Default::default()
    }
}

/// `None` for an input this build cannot load (mzML without the `mzml`
/// feature).
fn input(spec: &str, records: &BTreeMap<String, Synthetic>) -> Option<(MSExperiment, usize)> {
    if spec == "-" {
        return Some((MSExperiment::default(), 0));
    }
    if spec == "dta:noise" {
        let spectrum = openms::format::dta::read(
            include_bytes!("data/peak_picking_noise_input.dta").as_slice(),
        )
        .unwrap();
        return Some((
            MSExperiment {
                spectra: vec![spectrum],
                ..Default::default()
            },
            0,
        ));
    }
    if let Some(parts) = spec.strip_prefix("syn:") {
        let mut experiment = MSExperiment::default();
        for part in parts.split('+') {
            let record = &records[part];
            if record.kind == "S" {
                experiment.spectra.push(to_spectrum(record));
            } else {
                experiment.chromatograms.push(to_chromatogram(record));
            }
        }
        return Some((experiment, 0));
    }
    let label = spec.strip_prefix("mzml:").unwrap();
    let (label, index) = label.split_once('#').unwrap();
    mzml_input(label).map(|e| (e, index.parse().unwrap()))
}

#[cfg(feature = "mzml")]
fn mzml_input(label: &str) -> Option<MSExperiment> {
    use openms::format::mzml;
    let bytes: &[u8] = match label {
        "orbitrap" => include_bytes!("data/peak_picking/PeakPickerHiRes_orbitrap.mzML"),
        "topp2" => include_bytes!("data/peak_picking/PeakPickerHiRes_2_input.mzML"),
        other => panic!("unknown mzML label {other}"),
    };
    // MzMLFile().load with its default options, as the driver loads.
    Some(
        mzml::read_with_load_options(
            bytes,
            &mzml::LoadOptions::default(),
            &mzml::ReadOptions::default(),
        )
        .unwrap(),
    )
}

#[cfg(not(feature = "mzml"))]
fn mzml_input(_: &str) -> Option<MSExperiment> {
    None
}

// ---------------------------------------------------------------------------
// Captured output.

#[derive(Clone, Default)]
struct Capture(Arc<Mutex<Vec<u8>>>);
impl Capture {
    fn bytes(&self) -> Vec<u8> {
        self.0.lock().unwrap().clone()
    }
}
impl Write for Capture {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// What the source's thread-local `OPENMS_LOG_WARN` stream writes for these
/// estimations, at process exit: the stream wrote to a file, where its
/// colorizer is inactive, and its duplicate cache reports repeated lines when
/// it is destroyed.
fn rendered_log(estimations: &[&NoiseEstimates]) -> Vec<u8> {
    let capture = Capture::default();
    let mut log = LogStream::new("WARNING").unwrap();
    log.set_color(None);
    log.insert(&LogSink::new(capture.clone())).unwrap();
    for estimation in estimations {
        for line in &estimation.log {
            writeln!(log, "{line}").unwrap();
            log.flush().unwrap();
        }
    }
    log.finish().unwrap();
    capture.bytes()
}

/// The driver's timing mask, `sed -E 's/took [^(]* \(CPU\), [^(]* \(Wall\)/took <T> (CPU), <T> (Wall)/g'`.
fn mask(text: &[u8]) -> Vec<u8> {
    let text = String::from_utf8(text.to_vec()).unwrap();
    let mut out = String::new();
    let mut rest = text.as_str();
    while let Some(start) = rest.find("took ") {
        let after = &rest[start + 5..];
        let masked = after.find(" (CPU), ").and_then(|cpu| {
            let wall_part = &after[cpu + 8..];
            wall_part
                .find(" (Wall)")
                .filter(|wall| !after[..cpu].contains('(') && !wall_part[..*wall].contains('('))
                .map(|wall| cpu + 8 + wall + 7)
        });
        match masked {
            Some(end) => {
                out.push_str(&rest[..start]);
                out.push_str("took <T> (CPU), <T> (Wall)");
                rest = &after[end..];
            }
            None => {
                out.push_str(&rest[..start + 5]);
                rest = after;
            }
        }
    }
    out.push_str(rest);
    out.into_bytes()
}

fn clock(state: Arc<Mutex<i64>>, tick: bool) -> ProgressClock {
    Arc::new(move || {
        let mut second = state.lock().unwrap();
        let value = *second;
        if tick {
            *second += 1;
        }
        Ok(ProgressTime {
            wall_second: value,
            wall_seconds: 0.0,
            cpu_seconds: Some(0.0),
        })
    })
}

/// A command-line progress logger whose throttle reads the driver's
/// interposed `time()`: constant 1000, or 1000 and one more per call.
fn command_logger(tick: bool) -> (ProgressLogger, Capture) {
    let seconds = Arc::new(Mutex::new(1000));
    let mut logger =
        ProgressLogger::with_clock_and_nesting(clock(seconds, tick), ProgressNesting::default());
    logger.set_log_type(ProgressLogType::Cmd);
    let capture = Capture::default();
    let fixed = Arc::new(Mutex::new(0));
    logger.set_logger(Box::new(CommandProgressLogger::with_clock(
        capture.clone(),
        clock(fixed, false),
    )));
    (logger, capture)
}

// ---------------------------------------------------------------------------
// Running and comparing.

fn fnv(values: &[f64]) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325_u64;
    for value in values {
        for byte in value.to_bits().to_le_bytes() {
            hash ^= u64::from(byte);
            hash = hash.wrapping_mul(0x0100_0000_01b3);
        }
    }
    hash
}

fn estimator(params: &str) -> SignalToNoiseEstimatorMedian {
    SignalToNoiseEstimatorMedian::from_param(&apply(
        SignalToNoiseEstimatorMedian::defaults().unwrap(),
        params,
    ))
    .unwrap()
}

fn estimate(
    case: &Case,
    experiment: &MSExperiment,
    index: usize,
    compatibility: &PickingCompatibility,
    progress: Option<&mut ProgressLogger>,
) -> Result<NoiseEstimates> {
    let estimator = estimator(&case.params);
    if case.op == "noise_chrom" || case.op == "progress_tick_chrom" {
        estimator.estimate_peaks(
            &experiment.chromatograms[index].peaks,
            compatibility,
            progress,
        )
    } else {
        estimator.estimate_peaks(&experiment.spectra[index].peaks, compatibility, progress)
    }
}

/// Differences between an estimation and the oracle, empty when they agree.
fn noise_differences(result: &Result<NoiseEstimates>, expected: &Expected) -> Vec<String> {
    let mut out = Vec::new();
    match (result, &expected.error) {
        (Err(Error::InvalidValue(message)), Some((class, what))) => {
            if class != "InvalidValue" || message != what {
                out.push(format!("error {message:?}, oracle {class} {what:?}"));
            }
        }
        (Err(error), _) => out.push(format!("error {error}, oracle {:?}", expected.error)),
        (Ok(_), Some(error)) => out.push(format!("no error, oracle {error:?}")),
        (Ok(e), None) => {
            if Some(e.max_intensity.to_bits()) != expected.max_intensity {
                out.push(format!(
                    "max_intensity {:016x}, oracle {:016x?}",
                    e.max_intensity.to_bits(),
                    expected.max_intensity
                ));
            }
            if Some(e.signal_to_noise.len()) != expected.count {
                out.push(format!(
                    "{} ratios, oracle {:?}",
                    e.signal_to_noise.len(),
                    expected.count
                ));
            }
            let (n, hash, bits) = expected.stn.as_ref().unwrap();
            if e.signal_to_noise.len() != *n || fnv(&e.signal_to_noise) != *hash {
                out.push(format!("ratio digest differs ({n} oracle ratios)"));
            }
            let actual: Vec<u64> = e.signal_to_noise.iter().map(|v| v.to_bits()).collect();
            if !bits.is_empty() && actual != *bits {
                let first = (0..bits.len()).find(|&i| actual.get(i) != Some(&bits[i]));
                out.push(format!("ratios differ, first at {first:?}"));
            }
            for (i, b) in &expected.stn_at {
                if actual.get(*i) != Some(b) {
                    out.push(format!("ratio {i} differs"));
                }
            }
            let percentages = (
                e.sparse_window_percent.to_bits(),
                e.histogram_rightmost_percent.to_bits(),
            );
            if Some(percentages) != expected.percentages {
                out.push(format!(
                    "percentages {percentages:016x?}, oracle {:016x?}",
                    expected.percentages
                ));
            }
        }
    }
    out
}

fn same_estimates(a: &NoiseEstimates, b: &NoiseEstimates) -> bool {
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    bits(&a.signal_to_noise) == bits(&b.signal_to_noise)
        && bits(&a.noise) == bits(&b.noise)
        && a.max_intensity.to_bits() == b.max_intensity.to_bits()
        && a.sparse_window_percent.to_bits() == b.sparse_window_percent.to_bits()
        && a.histogram_rightmost_percent.to_bits() == b.histogram_rightmost_percent.to_bits()
        && a.log == b.log
}

const NOISE_OPS: &[&str] = &[
    "noise",
    "noise_chrom",
    "noise_twice",
    "progress_const",
    "progress_tick",
    "progress_tick_chrom",
];

/// Cases whose input or parameters the native safety profile refuses.
const NATIVE_REFUSALS: &[&str] = &[
    "negative_all_defaults",
    "warn_negative_range",
    "warn_negative_range_nolog",
    "twice_negative_range",
    "nfew_zero",
    "nfew_zero_pos",
    "nfew_negative",
    "nfew_inf",
    "nfew_neginf",
    "nfew_nan",
    "nfew_tiny",
    "nfew_tiny_pos",
    "win_inf",
    "win_inf_manual",
    "win_nan",
    "win_huge_positions",
    "stdev_factor_nan",
    "cpp257_neg_factor0",
    "cpp257_neg_default",
    "nf_nan_manual",
    "nf_nan_stdev",
    "nf_inf_stdev",
    "nf_inf_manual",
    "nf_neginf_stdev",
    "nf_order_stdev",
    "nf_order2_stdev",
    "nf_ratio_order",
    "nf_pos_nan",
    "nf_pos_inf",
    "nf_pos_neginf_winf",
    "nf_pos_inf_winf",
    "progress_early_return",
];

/// Cases the native safety profile computes differently, by design: an empty
/// input (zeros for the source's NaNs) and the clamp-first bin index.
const NATIVE_EMPTY: &[&str] = &[
    "empty_stdev",
    "empty_manual",
    "empty_stdev_chrom",
    "progress_empty",
];
const NATIVE_CLAMP: &[&str] = &[
    "cpp257_manual",
    "cpp257_stdev_bins",
    "cpp257_manual_chrom",
    "cpp257_boundary",
    "cpp257_pick_noise",
];

#[test]
fn estimator_cases_match_the_release_build_bit_for_bit() {
    let records = synthetic();
    let expected = oracle();
    let cases = cases();
    assert_eq!(expected.len(), cases.len());
    let source = PickingCompatibility::source();
    let clamp_only = PickingCompatibility {
        noise: NoiseCompatibility {
            bin_index: BinIndexConversion::ClampBeforeTruncation,
            ..NoiseCompatibility::source()
        },
        ..source
    };
    let mut failures = Vec::new();
    let (mut compared, mut skipped) = (0, Vec::new());
    let (mut refused, mut empty, mut clamped) = (Vec::new(), Vec::new(), Vec::new());
    for case in cases.iter().filter(|c| NOISE_OPS.contains(&c.op.as_str())) {
        let oracle = &expected[&case.name];
        assert_eq!(oracle.exit, 0, "{}", case.name);
        let Some((experiment, index)) = input(&case.input, &records) else {
            skipped.push(case.name.clone());
            continue;
        };
        compared += 1;
        let mut fail = |what: String| failures.push(format!("{}: {what}", case.name));
        // The source profile, with the progress report the case asks for.
        let tick = case.op.starts_with("progress_tick");
        let (mut logger, stdout) = command_logger(tick);
        let progress = case.op.starts_with("progress");
        let result = estimate(
            case,
            &experiment,
            index,
            &source,
            progress.then_some(&mut logger),
        );
        for difference in noise_differences(&result, oracle) {
            fail(difference);
        }
        if progress && mask(&stdout.bytes()) != oracle.stdout {
            fail(format!(
                "progress report {:?}, oracle {:?}",
                String::from_utf8_lossy(&mask(&stdout.bytes())),
                String::from_utf8_lossy(&oracle.stdout)
            ));
        }
        if !progress && !oracle.stdout.is_empty() {
            fail("the oracle wrote to stdout".into());
        }
        let mut logs = Vec::new();
        let second;
        if let Ok(first) = &result {
            logs.push(first);
            if case.op == "noise_twice" {
                second = estimate(case, &experiment, index, &source, None).unwrap();
                if !same_estimates(first, &second) {
                    fail("the second estimation differs".into());
                }
                if oracle.twice != Some(true) {
                    fail(format!("oracle twice {:?}", oracle.twice));
                }
                logs.push(&second);
            }
        }
        let stderr = rendered_log(&logs);
        if stderr != oracle.stderr {
            fail(format!(
                "warnings {:?}, oracle {:?}",
                String::from_utf8_lossy(&stderr),
                String::from_utf8_lossy(&oracle.stderr)
            ));
        }
        // The native safety profile: identical, refused or a documented
        // difference.
        let native = estimate(
            case,
            &experiment,
            index,
            &PickingCompatibility::default(),
            None,
        );
        match (&native, &result) {
            (Err(_), Ok(_)) => refused.push(case.name.clone()),
            (Err(a), Err(b)) if a.to_string() != b.to_string() => {
                fail(format!("native error {a}, source error {b}"))
            }
            (Ok(_), Err(_)) => fail("the native profile accepted what the source refuses".into()),
            (Ok(n), Ok(s)) if !same_estimates(n, s) => {
                if s.signal_to_noise.is_empty() {
                    // Zeros for the source's 0 / 0.
                    let nan = NEGATIVE_NAN;
                    assert_eq!(s.sparse_window_percent.to_bits(), nan, "{}", case.name);
                    assert_eq!(s.histogram_rightmost_percent.to_bits(), nan);
                    assert_eq!(n.sparse_window_percent.to_bits(), 0);
                    assert_eq!(n.histogram_rightmost_percent.to_bits(), 0);
                    let manual = n.max_intensity == s.max_intensity;
                    assert!(manual || n.max_intensity.to_bits() == 0, "{}", case.name);
                    assert!(manual || s.max_intensity.to_bits() == nan, "{}", case.name);
                    empty.push(case.name.clone());
                } else {
                    // Exactly the bin conversion.
                    let lifted = estimate(case, &experiment, index, &clamp_only, None).unwrap();
                    if !same_estimates(n, &lifted) {
                        fail("the native profile differs beyond the bin conversion".into());
                    }
                    clamped.push(case.name.clone());
                }
            }
            _ => {}
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    let sorted = |names: &[&str]| {
        let mut v: Vec<String> = names.iter().map(|s| (*s).to_owned()).collect();
        v.retain(|n| !skipped.contains(n));
        v.sort();
        v
    };
    for list in [&mut refused, &mut empty, &mut clamped] {
        list.sort();
    }
    assert_eq!(refused, sorted(NATIVE_REFUSALS));
    assert_eq!(empty, sorted(NATIVE_EMPTY));
    assert_eq!(clamped, sorted(NATIVE_CLAMP));
    assert_eq!(compared + skipped.len(), 89);
    #[cfg(feature = "mzml")]
    assert!(skipped.is_empty(), "{skipped:?}");
    #[cfg(not(feature = "mzml"))]
    assert_eq!(
        skipped,
        [
            "orbitrap0_defaults",
            "orbitrap0_manual",
            "chrom_topp2_defaults",
            "chrom_topp2_manual"
        ]
    );
}

#[test]
fn cpp257_bins_follow_the_release_build_and_the_native_profile_clamps_first() {
    // Independent derivation of cpp257_manual: width max(1, 1 / 30) = 1, the
    // intensities 1..9 fill bins 1..9, and 3e9 has quotient 3e9 > 2^31 - 1.
    // Release: cvttsd2si gives INT_MIN, clamped to bin 0, so the ten-point
    // window has one point in each of bins 0..9; the median (5th point) is in
    // bin 4 with 4 points before it: noise = 4 + (5 - 4) / 1 = 5.
    // Clamp first: 3e9 goes to bin 29; the median is in bin 5 with 4 before
    // it: noise = 5 + 1 = 6.
    let records = synthetic();
    let spectrum = to_spectrum(&records["cpp257_manual"]);
    let estimator = estimator("auto_mode=-1;max_intensity=1;min_required_elements=1");
    let release = estimator
        .estimate_spectrum(&spectrum, &PickingCompatibility::source())
        .unwrap();
    let native = estimator
        .estimate_spectrum(&spectrum, &PickingCompatibility::default())
        .unwrap();
    assert_eq!(release.noise, [5.0; 10]);
    assert_eq!(native.noise, [6.0; 10]);
    assert_eq!(release.signal_to_noise[0], 1.0 / 5.0);
    assert_eq!(native.signal_to_noise[4], 3e9 / 6.0);
    // The conversion boundary: 2^31 - 128 still converts, 2^31 does not.
    let boundary = to_spectrum(&records["cpp257_boundary"]);
    let small = SignalToNoiseEstimatorMedian {
        bin_count: 3,
        ..estimator.clone()
    };
    let release = small
        .estimate_spectrum(&boundary, &PickingCompatibility::source())
        .unwrap();
    // Bins: 2, 0, 2, 0 -> median (2nd point) in bin 0 with nothing before it.
    assert_eq!(release.noise, [1.0; 4]);
    let native = small
        .estimate_spectrum(&boundary, &PickingCompatibility::default())
        .unwrap();
    // Bins: 2, 2, 2, 2 -> median in bin 2: 2 + 2 / 4.
    assert_eq!(native.noise, [2.5; 4]);
}

#[test]
fn the_tool_path_bins_like_the_release_build() {
    // libOpenMS's own instantiation inside PeakPickerHiRes::pick, oracle case
    // cpp257_pick_manual. Release bins the five apex samples (above 2^31) in
    // bin 0: noise 29 + (13 - 5) / 20 = 29.4, and 2.9e9 / 29.4 passes the
    // threshold 98.5e6. Clamp first puts all 25 samples in bin 29: noise 29.52,
    // and 2.9e9 / 29.52 does not, so the native profile picks nothing.
    let records = synthetic();
    let expected = oracle();
    let case = cases()
        .into_iter()
        .find(|c| c.name == "cpp257_pick_manual")
        .unwrap();
    let oracle = &expected[&case.name];
    let spectrum = to_spectrum(&records["cpp257_pick"]);
    let param = apply(PeakPickerHiRes::defaults().unwrap(), &case.params);
    let mut picker = PeakPickerHiRes::from_param(&param).unwrap();
    picker.compatibility = PickingCompatibility::source();
    let picked = picker.pick_spectrum_with_spacing(&spectrum, true).unwrap();
    let peaks: Vec<(u64, u32)> = picked
        .spectrum
        .peaks
        .iter()
        .map(|p| (p.mz.to_bits(), p.intensity.to_bits()))
        .collect();
    let bounds: Vec<(u64, u64)> = picked
        .boundaries
        .iter()
        .map(|b| (b.min.to_bits(), b.max.to_bits()))
        .collect();
    assert_eq!(peaks, oracle.peaks);
    assert_eq!(bounds, oracle.bounds);
    assert_eq!(peaks.len(), 1);
    // The picker estimated once, with the warnings the direct estimate returns.
    let direct = picker
        .noise_estimator
        .estimate_spectrum(&spectrum, &PickingCompatibility::source())
        .unwrap();
    assert_eq!(direct.noise, [29.4; 25]);
    assert_eq!(rendered_log(&[&direct]), oracle.stderr);
    picker.compatibility = PickingCompatibility::default();
    assert!(
        picker
            .pick_spectrum_with_spacing(&spectrum, true)
            .unwrap()
            .spectrum
            .is_empty()
    );
    let native = picker
        .noise_estimator
        .estimate_spectrum(&spectrum, &PickingCompatibility::default())
        .unwrap();
    assert_eq!(native.noise, [29.0 + 13.0 / 25.0; 25]);
}

#[test]
fn random_scan_cases_match_the_release_build_bit_for_bit() {
    let records = synthetic();
    let expected = oracle();
    let mut compared = 0;
    for case in cases() {
        let oracle = &expected[&case.name];
        let k = keys(&case.params);
        match case.op.as_str() {
            "random" => {
                let (experiment, _) = input(&case.input, &records).unwrap();
                let result = estimate_noise_from_random_scans(
                    &experiment,
                    k["ms_level"].parse().unwrap(),
                    k["n_scans"].parse().unwrap(),
                    number(k["percentile"]),
                    k["seed"].parse().unwrap(),
                )
                .unwrap_or_else(|e| panic!("{}: {e}", case.name));
                assert_eq!(
                    Some(result.to_bits()),
                    oracle.random,
                    "{}: {result}",
                    case.name
                );
                // The seed is read from time(), which the driver interposed:
                // once per call that has candidates.
                let candidates = case.name != "rnd_ms3_none";
                assert_eq!(oracle.time_calls, u64::from(candidates), "{}", case.name);
                compared += 1;
            }
            "rng" => {
                let seed: u64 = k["seed"].parse().unwrap();
                let draws: usize = k["draws"].parse().unwrap();
                let scale: f64 = k["scale"].parse().unwrap();
                let mut raw = MinstdRand0::new(seed);
                let values: Vec<u64> = (0..draws).map(|_| u64::from(raw.next_u32())).collect();
                assert_eq!(values, oracle.raw, "{}", case.name);
                let mut engine = MinstdRand0::new(seed);
                let uniform: Vec<(u64, u32)> = (0..draws)
                    .map(|_| {
                        let u = engine.uniform01();
                        (u.to_bits(), (u * scale) as u32)
                    })
                    .collect();
                assert_eq!(uniform, oracle.draws, "{}", case.name);
                compared += 1;
            }
            _ => {}
        }
    }
    assert_eq!(compared, 52);
    // The permutation records exist for every vector of nth_vectors.tsv.
    assert_eq!(expected["nth_vectors"].nth_vectors, 390);
}

#[test]
fn random_scans_reproduce_the_defined_quirks() {
    let records = synthetic();
    // One MS2 candidate at experiment index 1: the scale is 1 - 1 = 0, so
    // every draw reads experiment index 0, the MS1 spectrum r_ms1a with the
    // intensities 0..=11 in some order. The position is (Size)(12 * 80 / 100.0)
    // = 9, whose order statistic is 9; five draws average to 9.
    let (experiment, _) = input("syn:r_ms1a+r_ms2a", &records).unwrap();
    for seed in [0, 1, 99, u64::MAX] {
        let noise = estimate_noise_from_random_scans(&experiment, 2, 5, 80.0, seed).unwrap();
        assert_eq!(noise, 9.0);
    }
    // No candidate: 0; no scan: 0 / 0, the default f32 NaN.
    assert_eq!(
        estimate_noise_from_random_scans(&experiment, 3, 10, 80.0, 1).unwrap(),
        0.0
    );
    assert_eq!(
        estimate_noise_from_random_scans(&experiment, 2, 0, 80.0, 1)
            .unwrap()
            .to_bits(),
        0xffc0_0000
    );
    // A percentile below 0 but above -100 / size still selects the minimum.
    assert_eq!(
        estimate_noise_from_random_scans(&experiment, 2, 3, -8.0, 1).unwrap(),
        0.0
    );
    // Seeds are taken modulo 2^31 - 1 with zero replaced by one.
    assert_eq!(MinstdRand0::new(0), MinstdRand0::new(1));
    assert_eq!(MinstdRand0::new(2_147_483_647), MinstdRand0::new(1));
    assert_eq!(MinstdRand0::new(2_147_483_648), MinstdRand0::new(1));
    assert_ne!(MinstdRand0::new(u64::MAX), MinstdRand0::new(1));
}

#[test]
fn random_scans_refuse_exactly_where_the_source_is_undefined() {
    let records = synthetic();
    // The only ms_level 2 candidate is r_ms2a at index 1, so every draw reads
    // experiment index 0 = r_ms1a, whose 12 intensities have minimum 0.
    let (experiment, _) = input("syn:r_ms1a+r_ms2a", &records).unwrap();
    let undefined = |result: Result<f32>, line: &str| match result {
        Err(Error::Unsupported(message)) => assert!(message.contains(line), "{message}"),
        other => panic!("{other:?}"),
    };
    // Percentile 100: e == size, the read at :50 is one past the end.
    undefined(
        estimate_noise_from_random_scans(&experiment, 2, 1, 100.0, 1),
        "SignalToNoiseEstimator.cpp:49-50",
    );
    // Above 100: e > size, the nth pointer passes the end at :49.
    undefined(
        estimate_noise_from_random_scans(&experiment, 2, 1, 150.0, 1),
        "SignalToNoiseEstimator.cpp:49-50",
    );
    // An ordinary negative percentile: idx wraps to e ~ 2^62, far past the end.
    undefined(
        estimate_noise_from_random_scans(&experiment, 2, 1, -50.0, 1),
        "SignalToNoiseEstimator.cpp:49-50",
    );
    // -100 / 12 or below converts a value <= -1, whose wrap is also past the end.
    undefined(
        estimate_noise_from_random_scans(&experiment, 2, 1, -8.34, 1),
        "SignalToNoiseEstimator.cpp:49-50",
    );
    // An empty drawn scan: e == size == 0, the read at :50 is out of bounds.
    let mut empty_first = experiment.clone();
    empty_first.spectra.insert(0, MSSpectrum::default());
    undefined(
        estimate_noise_from_random_scans(&empty_first, 2, 1, 80.0, 1),
        "SignalToNoiseEstimator.cpp:49-50",
    );
    // A NaN, an infinite, or a percentile giving a product of -2^63 or below
    // wraps to element 0 (the minimum, 0 here), which is in bounds and
    // deterministic, so the port computes it (the Release build too;
    // ../oracle/sne-followup pins the bits).
    for percentile in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1e30, 1e30] {
        assert_eq!(
            estimate_noise_from_random_scans(&experiment, 2, 1, percentile, 1).unwrap(),
            0.0,
            "percentile {percentile}"
        );
    }
    // A NaN in the drawn scan is not refused: the Release build's
    // std::nth_element stays in bounds (the rnd_nan_* and nth_* oracle cases
    // pin what it returns).
    let mut with_nan = experiment.clone();
    with_nan.spectra[0].peaks[3].intensity = f32::NAN;
    assert!(estimate_noise_from_random_scans(&with_nan, 2, 1, 50.0, 1).is_ok());
    // The native work ceiling: every draw reads experiment index 0 (12
    // intensities) and costs 12 + 1, so three draws need 39.
    let limited = RandomScanNoise {
        n_scans: 3,
        max_work: 39,
        ..RandomScanNoise::new(2, 1)
    };
    assert!(limited.estimate(&experiment).is_ok());
    assert!(matches!(
        RandomScanNoise {
            max_work: 38,
            ..limited
        }
        .estimate(&experiment),
        Err(Error::InvalidValue(_))
    ));
}

/// The emulated pointer wrap of `SignalToNoiseEstimator.cpp:49-50`: when the
/// drawn position `idx` maps to an in-bounds element `e = idx mod 2^62`, the
/// Release build's `nth_element(begin, begin + idx, end)` and `tmp[idx]` are
/// ordinary in-bounds operations on element `e`, so the port computes them.
/// The fixtures are the Linux x86-64 Release build's own output, twice each
/// byte-identically (`../oracle/sne-followup/`, same driver binary
/// `571d6f9b…` as `sne-completion`).
#[test]
fn random_scan_pointer_wrap_matches_the_release_build() {
    const WRAP_CASES: &str = include_str!("data/signal_to_noise/wrap_cases.tsv");
    const WRAP_SYNTHETIC: &str = include_str!("data/signal_to_noise/wrap_synthetic.tsv");
    const WRAP_ORACLE: &str = include_str!("data/signal_to_noise/wrap_oracle.tsv");

    // The drawn intensity each Release case returned, as f32 bits.
    let mut expected: BTreeMap<&str, u32> = BTreeMap::new();
    for line in WRAP_ORACLE.lines() {
        let c: Vec<&str> = line.split('\t').collect();
        if c[0] == "random" {
            expected.insert(c[1], hex32(c[2]));
        }
    }

    let records: BTreeMap<String, Synthetic> = WRAP_SYNTHETIC
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .map(|line| {
            let c: Vec<&str> = line.split('\t').collect();
            let x = list(c[3]);
            let y: Vec<f32> = list(c[4]).into_iter().map(narrow).collect();
            assert_eq!(x.len(), y.len(), "{}", c[0]);
            (
                c[0].to_owned(),
                Synthetic {
                    kind: c[1].to_owned(),
                    ms_level: c[2].parse().unwrap(),
                    x,
                    y,
                },
            )
        })
        .collect();

    let mut compared = 0;
    for line in WRAP_CASES.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        assert_eq!(c[1], "random", "{}", c[0]);
        let (experiment, _) = input(c[2], &records).unwrap();
        let k = keys(c[3]);
        let result = estimate_noise_from_random_scans(
            &experiment,
            k["ms_level"].parse().unwrap(),
            k["n_scans"].parse().unwrap(),
            number(k["percentile"]),
            k["seed"].parse().unwrap(),
        )
        .unwrap_or_else(|e| panic!("{}: {e}", c[0]));
        assert_eq!(
            result.to_bits(),
            expected[c[0]],
            "{}: got {result} ({:#010x})",
            c[0],
            result.to_bits()
        );
        compared += 1;
    }
    assert_eq!(compared, 10);
    // The in-bounds low-path and high-path elements are pinned to their exact
    // values (element index = f64 bits): 1024, 2048 and 3072.
    assert_eq!(expected["wrap_low_2p62_j1"], 1024.0_f32.to_bits());
    assert_eq!(expected["wrap_low_2p62_j2"], 2048.0_f32.to_bits());
    assert_eq!(expected["wrap_low_2p62_j3"], 3072.0_f32.to_bits());
    assert_eq!(expected["wrap_high_2p63"], 2048.0_f32.to_bits());
    // The e = 0 percentiles (NaN, +/-inf, 1e30, a product of -2^63 or below)
    // all read element 0, the minimum (0 here).
    for name in [
        "wrap_nan",
        "wrap_neg_inf",
        "wrap_pos_inf",
        "wrap_dle_neg2p63",
        "wrap_1e30",
    ] {
        assert_eq!(expected[name], 0.0_f32.to_bits(), "{name}");
    }
}

#[test]
fn percentile_range_is_refused_exactly_outside_its_domain() {
    let percentile = SignalToNoiseEstimatorMedian {
        histogram_range: NoiseHistogramRange::Percentile { percentile: 95.0 },
        min_required_elements: 1,
        ..Default::default()
    };
    let refused = |x: &[f64], y: &[f64], line: &str| {
        for compatibility in [
            PickingCompatibility::default(),
            PickingCompatibility::source(),
        ] {
            match percentile.estimate_with_compatibility(x, y, &compatibility) {
                Err(Error::Unsupported(message)) => {
                    assert!(message.contains(line), "{message}")
                }
                other => panic!("{y:?}: {other:?}"),
            }
        }
    };
    let accepted = |x: &[f64], y: &[f64]| {
        assert!(percentile.estimate(x, y).is_ok(), "{y:?}");
    };
    refused(&[], &[], "SignalToNoiseEstimatorMedian.h:209");
    // The smallest f32 minimum of a constant input with quotient above -1 is
    // 0.9900990128517151 (the oracle's pctl_lower_edge); the f32 below it is
    // outside.
    let x12: Vec<f64> = (0..12).map(f64::from).collect();
    accepted(&x12, &[0.990_099_012_851_715_1; 12]);
    refused(
        &x12,
        &[0.990_098_953_247_070_3; 12],
        "SignalToNoiseEstimatorMedian.h:216",
    );
    // For minimum 3 the largest f32 inside is 3.999999761581421
    // (pctl_upper_edge); 4 has quotient 100.0000002.
    let x2 = [0.0, 1.0];
    accepted(&x2, &[3.0, 3.999_999_761_581_421]);
    refused(&x2, &[3.0, 4.0], "SignalToNoiseEstimatorMedian.h:216");
    // A zero or negative minimum, a NaN and an infinity are all outside.
    refused(&x2, &[0.0, 0.5], "SignalToNoiseEstimatorMedian.h:216");
    let negative = PickingCompatibility::source();
    assert!(matches!(
        percentile.estimate_with_compatibility(&x2, &[-5.0, -4.5], &negative),
        Err(Error::Unsupported(_))
    ));
    for y in [[f64::NAN, 5.0], [5.0, f64::INFINITY]] {
        assert!(matches!(
            percentile.estimate_with_compatibility(&x2, &y, &negative),
            Err(Error::Unsupported(_))
        ));
    }
    // 10 + 0.33 i (30 points) has quotients up to 185.7, outside the domain:
    // the Release build ends with SIGABRT on it in all three runs of probe
    // x_pctl_spread33 (../oracle/sne-completion/probes/results.txt).
    let x30: Vec<f64> = (0..30).map(f64::from).collect();
    let y30: Vec<f64> = (0..30)
        .map(|i| f64::from((10.0 + 0.33 * f64::from(i)) as f32))
        .collect();
    refused(&x30, &y30, "SignalToNoiseEstimatorMedian.h:216");
    // The source's own range check, with its exception text.
    let outside = SignalToNoiseEstimatorMedian {
        histogram_range: NoiseHistogramRange::Percentile { percentile: 101.0 },
        ..percentile.clone()
    };
    assert!(matches!(
        outside.estimate(&x2, &[3.0, 3.5]),
        Err(Error::InvalidValue(message)) if message.contains("auto_max_percentile is not in [0,100]")
    ));
}

#[test]
fn progress_reports_stay_balanced_and_default_silent() {
    let estimator = SignalToNoiseEstimatorMedian {
        max_work: 3,
        ..Default::default()
    };
    let nesting = ProgressNesting::default();
    let seconds = Arc::new(Mutex::new(5));
    let mut logger = ProgressLogger::with_clock_and_nesting(clock(seconds, true), nesting.clone());
    logger.set_log_type(ProgressLogType::Cmd);
    let capture = Capture::default();
    logger.set_logger(Box::new(CommandProgressLogger::with_clock(
        capture.clone(),
        clock(Arc::new(Mutex::new(0)), false),
    )));
    let x = [0.0, 1.0, 2.0, 3.0];
    let y = [1.0, 2.0, 3.0, 4.0];
    // The work ceiling stops the second window after the report started.
    assert!(
        estimator
            .estimate_with_progress(&x, &y, &PickingCompatibility::default(), &mut logger)
            .is_err()
    );
    assert_eq!(nesting.depth(), 0);
    let text = String::from_utf8(capture.bytes()).unwrap();
    assert!(text.starts_with(&format!("Progress of '{NOISE_PROGRESS_LABEL}':\n")));
    assert!(text.ends_with("] -- \n"), "{text:?}");
    // A refusal before the report starts writes nothing.
    let before = capture.bytes().len();
    let manual = SignalToNoiseEstimatorMedian {
        histogram_range: NoiseHistogramRange::Manual { max_intensity: 0.0 },
        ..Default::default()
    };
    assert!(
        manual
            .estimate_with_progress(&x, &y, &PickingCompatibility::default(), &mut logger)
            .is_err()
    );
    assert_eq!(capture.bytes().len(), before);
    assert_eq!(nesting.depth(), 0);
}

#[test]
fn estimates_through_the_base_trait() {
    let spectrum = MSSpectrum::from_peaks(
        (0..40)
            .map(|i| Peak1D::new(100.0 + f64::from(i), ((i * 7) % 13) as f32))
            .collect(),
    );
    let chromatogram = MSChromatogram {
        peaks: spectrum
            .peaks
            .iter()
            .map(|p| ChromatogramPeak::new(p.mz, p.intensity))
            .collect(),
        ..Default::default()
    };
    fn ratios<E: SignalToNoiseEstimator>(
        estimator: &E,
        spectrum: &MSSpectrum,
        chromatogram: &MSChromatogram,
    ) -> (Vec<f64>, Vec<f64>) {
        let a = estimator.compute_stn_spectrum(spectrum).unwrap();
        let b = estimator.compute_stn_chromatogram(chromatogram).unwrap();
        (
            E::signal_to_noise(&a).to_vec(),
            E::signal_to_noise(&b).to_vec(),
        )
    }
    let median = SignalToNoiseEstimatorMedian {
        window_length: 10.0,
        min_required_elements: 3,
        ..Default::default()
    };
    let (a, b) = ratios(&median, &spectrum, &chromatogram);
    assert_eq!(a, b);
    let x: Vec<f64> = spectrum.peaks.iter().map(|p| p.mz).collect();
    let y: Vec<f64> = spectrum
        .peaks
        .iter()
        .map(|p| f64::from(p.intensity))
        .collect();
    assert_eq!(a, median.estimate(&x, &y).unwrap().signal_to_noise);
    assert_eq!(a, median.compute_stn(&x, &y).unwrap().signal_to_noise);
    let mean = SignalToNoiseEstimatorMeanIterative {
        window_length: 10.0,
        min_required_elements: 3,
        ..Default::default()
    };
    let (c, d) = ratios(&mean, &spectrum, &chromatogram);
    assert_eq!(c, d);
    assert_eq!(c, mean.estimate(&x, &y).unwrap().signal_to_noise);
    // The base's Gaussian estimate: population variance.
    let g = GaussianEstimate::of(&y);
    let n = y.len() as f64;
    let m = y.iter().sum::<f64>() / n;
    assert_eq!(g.mean, m);
    assert_eq!(
        g.variance,
        y.iter().map(|v| (m - v) * (m - v)).sum::<f64>() / n
    );
}

#[test]
fn parameters_are_set_in_place_and_warnings_are_gated() {
    let mut estimator = SignalToNoiseEstimatorMedian {
        max_points: 7,
        ..Default::default()
    };
    let param = apply(
        SignalToNoiseEstimatorMedian::defaults().unwrap(),
        "win_len=1.0;noise_for_empty_window=2.0;write_log_messages=false",
    );
    assert!(estimator.set_parameters(&param).unwrap().is_empty());
    assert_eq!(estimator.max_points, 7);
    assert_eq!(estimator.window_length, 1.0);
    assert!(!estimator.write_log_messages);
    // `Param::setValue(key, value)` replaces an entry's tags, so the driver's
    // tree lost the `advanced` tag that SignalToNoiseEstimatorMedian.h:114
    // declares for noise_for_empty_window. The source's `param_` keeps the
    // caller's tree (Param.cpp:406-414, DefaultParamHandler.cpp:47-49); the
    // derived tree carries the defaults' tags with the same values.
    assert!(!param.has_tag("noise_for_empty_window", "advanced").unwrap());
    let (typed, handler, warnings) =
        SignalToNoiseEstimatorMedian::from_param_with_handler(&param).unwrap();
    assert!(warnings.is_empty());
    assert_eq!(handler.parameters(), &param);
    assert_eq!(
        handler.defaults(),
        &SignalToNoiseEstimatorMedian::defaults().unwrap()
    );
    assert_eq!(typed.to_param().unwrap(), estimator.to_param().unwrap());
    let mut derived = param.clone();
    derived
        .add_tag("noise_for_empty_window", "advanced")
        .unwrap();
    assert_eq!(estimator.to_param().unwrap(), derived);
    let before = estimator.clone();
    let invalid = apply(
        SignalToNoiseEstimatorMedian::defaults().unwrap(),
        "bin_count=2",
    );
    assert!(estimator.set_parameters(&invalid).is_err());
    assert_eq!(estimator, before);
    // Gated: 100 % sparse windows but no warning.
    let quiet = estimator.estimate(&[0.0, 1000.0], &[2.0, 4.0]).unwrap();
    assert_eq!(quiet.sparse_window_percent, 100.0);
    assert!(quiet.log.is_empty());
    // Not gated: the negative-range warning.
    let negative = estimator
        .estimate_with_compatibility(
            &[0.0, 1.0, 2.0],
            &[-100.0, -99.0, -100.0],
            &PickingCompatibility::source(),
        )
        .unwrap();
    assert_eq!(negative.log.len(), 1);
    assert!(negative.log[0].starts_with("SignalToNoiseEstimatorMedian: the max_intensity_"));
}

#[test]
fn percentages_use_count_times_hundred_over_n() {
    // 3 of 7 windows sparse: 300 / 7 = 42.857142857142854, while
    // 3 * (100 / 7) = 42.85714285714286. Oracle case pct_sparse_3of7 pins the
    // first bit for bit; this pins the distinction.
    assert_ne!(
        (3.0_f64 * 100.0 / 7.0).to_bits(),
        (3.0_f64 * (100.0 / 7.0)).to_bits()
    );
    let records = synthetic();
    let spectrum = to_spectrum(&records["pct_sparse_3of7"]);
    let estimates = estimator("win_len=1.0;min_required_elements=2;noise_for_empty_window=3.0")
        .estimate_spectrum(&spectrum, &PickingCompatibility::default())
        .unwrap();
    assert_eq!(
        estimates.sparse_window_percent.to_bits(),
        (3.0_f64 * 100.0 / 7.0).to_bits()
    );
}

/// `SignalToNoiseEstimator_test.cpp` (6 sections) exercises the abstract class
/// through `TestSignalToNoiseEstimator`, whose `computeSTN_` does nothing: the
/// constructor, the copy constructor, assignment and the destructor, then
/// `init` on an empty spectrum; `getSignalToNoise` is `NOT_TESTABLE` there.
/// The same shape here: a trait implementation that estimates nothing.
#[test]
fn class_test_base_sections_through_a_trivial_estimator() {
    #[derive(Clone, Debug, Default, PartialEq)]
    struct TestSignalToNoiseEstimator;
    impl SignalToNoiseEstimator for TestSignalToNoiseEstimator {
        type Estimates = Vec<f64>;
        fn compute_stn(&self, positions: &[f64], _: &[f64]) -> Result<Vec<f64>> {
            // "do nothing here": no ratio is computed; the source's vector
            // keeps its previous (here: zero) contents.
            Ok(vec![0.0; positions.len()])
        }
        fn signal_to_noise(estimates: &Self::Estimates) -> &[f64] {
            estimates
        }
    }
    let estimator = TestSignalToNoiseEstimator;
    let copy = estimator.clone();
    let mut assigned = TestSignalToNoiseEstimator;
    assigned.clone_from(&copy);
    assert_eq!(assigned, copy);
    let empty = estimator
        .compute_stn_spectrum(&MSSpectrum::default())
        .unwrap();
    assert!(TestSignalToNoiseEstimator::signal_to_noise(&empty).is_empty());
    let chromatogram = MSChromatogram {
        peaks: vec![ChromatogramPeak::new(1.0, 2.0)],
        ..Default::default()
    };
    let one = estimator.compute_stn_chromatogram(&chromatogram).unwrap();
    assert_eq!(TestSignalToNoiseEstimator::signal_to_noise(&one), [0.0]);
    // The index access the source leaves unchecked in Release is a checked
    // slice access here.
    assert_eq!(
        TestSignalToNoiseEstimator::signal_to_noise(&one).get(1),
        None
    );
}

/// The warning lines print their value as `std::ostream << double` does at the
/// default precision, i.e. `printf("%g")`. The expected spellings were derived
/// independently, identically with glibc's `snprintf("%g")` on `ibminode06`
/// and with Python's `'%g' %`: scientific notation below `1e-4` and from
/// `1e6` after rounding, six significant digits, exact ties to even
/// (`93.65625`, `123456.5`, `1234565`), and a carry into the next decade
/// (`999999.5`). A single negative point makes the standard-deviation range
/// exactly that value, so each takes the ungated early return.
#[test]
fn warning_values_use_the_stream_default_format() {
    let cases = [
        (1.5e-05, "-1.5e-05"),
        (123_456_789.0, "-1.23457e+08"),
        (99.999_949, "-99.9999"),
        (999_999.5, "-1e+06"),
        (9.999_999_747_378_75e-05, "-0.0001"),
        (93.656_25, "-93.6562"),
        (5.0, "-5"),
        (1e-300, "-1e-300"),
        (0.5, "-0.5"),
        (123_456.5, "-123456"),
        (1_234_565.0, "-1.23456e+06"),
        (0.000_123_456_789, "-0.000123457"),
        (9.999_95e-05, "-9.99995e-05"),
        (5e-324, "-4.94066e-324"),
        (f64::MAX, "-1.79769e+308"),
    ];
    let estimator = SignalToNoiseEstimatorMedian {
        min_required_elements: 1,
        ..Default::default()
    };
    for (value, text) in cases {
        let estimates = estimator
            .estimate_with_compatibility(&[0.0], &[-value], &PickingCompatibility::source())
            .unwrap();
        assert_eq!(estimates.max_intensity, -value);
        assert_eq!(
            estimates.log,
            [format!(
                "SignalToNoiseEstimatorMedian: the max_intensity_ value should be positive! {text}"
            )]
        );
        assert_eq!(estimates.signal_to_noise, [0.0]);
    }
}
