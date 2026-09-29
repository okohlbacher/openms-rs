// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Progress and consumer calls of `MzXMLFile::transform`, both overloads,
//! replayed against executed OpenMS4 Release output.
//!
//! `tests/data/mzxml_transform_progress_release.tsv` was written by
//! `../oracle/F7c-transform-progress` (see
//! `tests/data/mzxml_transform_progress_provenance.json`): the C++ Release
//! build on ibminode06 ran every case twice on a fresh `MzXMLFile`, once with
//! `setLogType(GUI)` and a recording backend, once with `setLogType(CMD)` and
//! stdout captured (timing texts masked), with a consumer that records its
//! calls into the same event stream. Two repetitions agree row for row. The
//! driver replaced `time()` so every call is a new second; the clock here does
//! the same. No expected value in this file is derived from Rust output.

#![cfg(feature = "mzml")]

use openms::concept::progress_logger::{
    CommandProgressLogger, ProgressBackend, ProgressClock, ProgressLogType, ProgressLogger,
    ProgressNesting, ProgressTime,
};
use openms::format::mzxml::{self, MzXMLFile, TransformOptions};
use openms::interfaces::MSDataConsumer;
use openms::metadata::ExperimentalSettings;
use openms::{Error, MSChromatogram, MSExperiment, MSSpectrum, Result};
use std::collections::BTreeMap;
use std::io::{self, Write};
use std::ops::ControlFlow;
use std::path::PathBuf;
use std::sync::atomic::{AtomicI64, Ordering};
use std::sync::{Arc, Mutex};

const FIXTURE: &str = include_str!("data/mzxml_transform_progress_release.tsv");

type Events = Arc<Mutex<Vec<String>>>;

fn data(name: &str) -> PathBuf {
    let path = match name {
        "truncated.mzXML" => "progress_format_readers/truncated.mzXML",
        other => other,
    };
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(path)
}

/// A clock whose every read is a new whole second, as the driver's `time()`.
fn new_second_clock() -> ProgressClock {
    let second = Arc::new(AtomicI64::new(1_000_000_000));
    Arc::new(move || {
        let now = second.fetch_add(1, Ordering::Relaxed);
        Ok(ProgressTime {
            wall_second: now,
            wall_seconds: (now - 1_000_000_000) as f64,
            cpu_seconds: Some(0.0),
        })
    })
}

/// Records every backend call in the fixture's row format.
struct Recorder {
    events: Events,
    current: i64,
}
impl ProgressBackend for Recorder {
    fn start_progress(&mut self, begin: i64, end: i64, label: &str, depth: usize) -> Result<()> {
        self.current = begin;
        let row = format!("S\t{begin}\t{end}\t{label}\t{depth}");
        self.events.lock().unwrap().push(row);
        Ok(())
    }
    fn set_progress(&mut self, value: i64, depth: usize) -> Result<()> {
        self.events
            .lock()
            .unwrap()
            .push(format!("V\t{value}\t{depth}"));
        Ok(())
    }
    fn next_progress(&mut self) -> Result<i64> {
        self.current += 1;
        let row = format!("N\t{}", self.current);
        self.events.lock().unwrap().push(row);
        Ok(self.current)
    }
    fn end_progress(&mut self, depth: usize, bytes_processed: u64) -> Result<()> {
        self.events
            .lock()
            .unwrap()
            .push(format!("E\t{depth}\t{bytes_processed}"));
        Ok(())
    }
}

/// The driver's consumer: records its calls into the shared event stream.
struct Consumer {
    events: Events,
    spectra: usize,
}
impl MSDataConsumer for Consumer {
    fn set_expected_size(&mut self, spectra: usize, chromatograms: usize) -> Result<()> {
        self.events
            .lock()
            .unwrap()
            .push(format!("C\texpected\t{spectra}\t{chromatograms}"));
        Ok(())
    }
    fn set_experimental_settings(&mut self, _: &ExperimentalSettings) -> Result<()> {
        self.events.lock().unwrap().push("C\tsettings".into());
        Ok(())
    }
    fn consume_spectrum(&mut self, spectrum: &mut MSSpectrum) -> Result<ControlFlow<()>> {
        self.spectra += 1;
        let row = format!(
            "C\tspectrum\t{}\t{}",
            spectrum.native_id,
            spectrum.peaks.len()
        );
        self.events.lock().unwrap().push(row);
        Ok(ControlFlow::Continue(()))
    }
    fn consume_chromatogram(&mut self, _: &mut MSChromatogram) -> Result<ControlFlow<()>> {
        self.events.lock().unwrap().push("C\tchromatogram".into());
        Ok(ControlFlow::Continue(()))
    }
}

#[derive(Clone, Default)]
struct SharedOutput(Arc<Mutex<Vec<u8>>>);
impl Write for SharedOutput {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// A logger of type GUI whose factory makes recording backends (`rec`), or
/// command backends writing to `output` (`cmd`), as the driver's two modes.
fn make_logger(
    mode: &str,
    events: &Events,
    output: &SharedOutput,
) -> (ProgressLogger, ProgressNesting) {
    let nesting = ProgressNesting::default();
    let clock = new_second_clock();
    let mut logger = ProgressLogger::with_clock_and_nesting(clock.clone(), nesting.clone());
    if mode == "rec" {
        let sink = events.clone();
        logger.set_gui_factory(Arc::new(move || {
            Box::new(Recorder {
                events: sink.clone(),
                current: 0,
            })
        }));
    } else {
        let writer = output.clone();
        logger.set_gui_factory(Arc::new(move || {
            Box::new(CommandProgressLogger::with_clock(
                writer.clone(),
                clock.clone(),
            ))
        }));
    }
    logger.set_log_type(ProgressLogType::Gui);
    (logger, nesting)
}

/// Masks the timing texts of every summary line, as `extract.py` does.
fn mask_timing(text: &str) -> String {
    const HEAD: &str = "-- done [took ";
    let mut out = String::new();
    let mut rest = text;
    while let Some(start) = rest.find(HEAD) {
        let after = &rest[start + HEAD.len()..];
        let end = after.find("] -- ").expect("summary line end");
        assert!(!after[..end].contains(" @ "), "no byte count here");
        out.push_str(&rest[..start]);
        out.push_str("-- done [took <TIME> (CPU), <TIME> (Wall)] -- ");
        rest = &after[end + "] -- ".len()..];
    }
    out.push_str(rest);
    out
}

fn unescape(field: &str) -> String {
    field
        .replace("\\n", "\n")
        .replace("\\r", "\r")
        .replace("\\t", "\t")
}

/// One run of one case in one mode, as the Release build recorded it.
#[derive(Default)]
struct Captured {
    events: Vec<String>,
    outcome: Option<std::result::Result<Vec<String>, (String, String)>>,
    depth: Option<usize>,
    output: Option<String>,
}

fn fixture() -> (Vec<String>, BTreeMap<(String, String), Captured>) {
    let mut order = Vec::new();
    let mut runs: BTreeMap<(String, String), Captured> = BTreeMap::new();
    for line in FIXTURE.lines().filter(|l| !l.starts_with('#')) {
        let fields: Vec<&str> = line.split('\t').collect();
        let (case, mode, record, rest) = (fields[0], fields[1], fields[2], &fields[3..]);
        if !order.iter().any(|c| c == case) {
            order.push(case.to_string());
        }
        let run = runs.entry((case.into(), mode.into())).or_default();
        match record {
            "S" | "V" | "N" | "E" | "C" => run.events.push(fields[2..].join("\t")),
            "R" => run.outcome = Some(Ok(rest.iter().map(|s| s.to_string()).collect())),
            "X" => run.outcome = Some(Err((rest[0].into(), rest[1].into()))),
            "D" => run.depth = Some(rest[0].parse().unwrap()),
            "B" => run.output = Some(unescape(rest.first().copied().unwrap_or(""))),
            other => panic!("unknown fixture record {other:?}"),
        }
    }
    (order, runs)
}

/// Runs one case as the driver does, returning its `R` fields.
fn run_case(case: &str, logger: &mut ProgressLogger, events: &Events) -> Result<Vec<String>> {
    let skip = case.ends_with("_skip");
    let base = case.strip_suffix("_skip").unwrap_or(case);
    let consumer = || Consumer {
        events: events.clone(),
        spectra: 0,
    };
    let mut file = MzXMLFile::new();
    match base {
        "mzxml_transform" => {
            let mut c = consumer();
            file.transform_with_progress(data("MzXMLFile_1.mzXML"), &mut c, skip, logger)?;
            Ok(vec![c.spectra.to_string()])
        }
        "mzxml_transform_map" => {
            let mut c = consumer();
            let mut map = MSExperiment::default();
            file.transform_into_with_progress(
                data("MzXMLFile_1.mzXML"),
                &mut c,
                &mut map,
                skip,
                logger,
            )?;
            Ok(vec![c.spectra.to_string(), map.spectra.len().to_string()])
        }
        "mzxml_transform_truncated" => {
            let mut c = consumer();
            file.transform_with_progress(data("truncated.mzXML"), &mut c, skip, logger)?;
            Ok(vec![c.spectra.to_string()])
        }
        "mzxml_transform_twice" => {
            let mut first = consumer();
            file.transform_with_progress(data("MzXMLFile_1.mzXML"), &mut first, skip, logger)?;
            let mut second = consumer();
            file.transform_with_progress(data("MzXMLFile_1.mzXML"), &mut second, skip, logger)?;
            Ok(vec![first.spectra.to_string(), second.spectra.to_string()])
        }
        "mzxml_transform_filtered" => {
            file.options.peaks.add_ms_level(2)?;
            let mut c = consumer();
            file.transform_with_progress(data("MzXMLFile_1.mzXML"), &mut c, false, logger)?;
            Ok(vec![c.spectra.to_string()])
        }
        other => panic!("fixture case {other:?} has no Rust counterpart"),
    }
}

/// The class of the port's error for the source exception a case records.
fn check_outcome(
    case: &str,
    source: &std::result::Result<Vec<String>, (String, String)>,
    port: &Result<Vec<String>>,
) {
    match (source, port) {
        (Ok(source), Ok(port)) => assert_eq!(port, source, "{case}"),
        (Err((name, message)), Err(error)) => match name.as_str() {
            // Xerces' end-of-input check; the port names the same element.
            "Parse Error" => {
                assert!(
                    message.contains("input ended before all started tags were ended"),
                    "{case}: {message}"
                );
                assert!(
                    error
                        .to_string()
                        .contains("input ended before all started tags were ended"),
                    "{case}: {error}"
                );
            }
            // A second start on a running command backend (`StopWatch.cpp`).
            "Precondition failed" => {
                assert_eq!(message, "StopWatch is already started!", "{case}");
                assert!(matches!(error, Error::InvalidValue(_)), "{case}: {error}");
                assert_eq!(
                    error.to_string(),
                    "invalid value: progress timer is already running",
                    "{case}"
                );
            }
            other => panic!("{case}: no error mapping for {other}"),
        },
        (source, port) => panic!("{case}: source {source:?}, port {port:?}"),
    }
}

#[test]
fn fixture_covers_both_overloads_and_both_count_modes() {
    let (order, runs) = fixture();
    assert_eq!(order.len(), 9, "{order:?}");
    for case in &order {
        for mode in ["rec", "cmd"] {
            let run = &runs[&(case.clone(), mode.to_string())];
            assert!(
                run.outcome.is_some() && run.depth.is_some(),
                "{case} {mode}"
            );
        }
        let cmd = &runs[&(case.clone(), "cmd".into())];
        assert!(cmd.output.is_some(), "{case}");
        assert!(
            cmd.events.iter().all(|e| e.starts_with('C')),
            "{case}: a command run records only consumer calls"
        );
    }
}

#[test]
fn every_backend_and_consumer_call_matches_the_release_build() {
    let (order, runs) = fixture();
    for case in &order {
        let captured = &runs[&(case.clone(), "rec".into())];
        let events = Events::default();
        let (mut logger, nesting) = make_logger("rec", &events, &SharedOutput::default());
        let outcome = run_case(case, &mut logger, &events);
        assert_eq!(
            *events.lock().unwrap(),
            captured.events,
            "{case}: backend and consumer calls"
        );
        assert_eq!(nesting.depth(), captured.depth.unwrap(), "{case}: depth");
        check_outcome(case, captured.outcome.as_ref().unwrap(), &outcome);
    }
}

#[test]
fn command_output_matches_the_release_build() {
    let (order, runs) = fixture();
    for case in &order {
        let captured = &runs[&(case.clone(), "cmd".into())];
        let events = Events::default();
        let output = SharedOutput::default();
        let (mut logger, nesting) = make_logger("cmd", &events, &output);
        let outcome = run_case(case, &mut logger, &events);
        let text = String::from_utf8(output.0.lock().unwrap().clone()).unwrap();
        assert_eq!(
            mask_timing(&text),
            *captured.output.as_ref().unwrap(),
            "{case}: stdout"
        );
        assert_eq!(*events.lock().unwrap(), captured.events, "{case}: consumer");
        assert_eq!(nesting.depth(), captured.depth.unwrap(), "{case}: depth");
        check_outcome(case, captured.outcome.as_ref().unwrap(), &outcome);
    }
}

/// Reporting progress changes no result: the progress entry points deliver
/// what the silent ones deliver.
#[test]
fn progress_changes_no_transform_result() {
    for skip in [false, true] {
        let options = TransformOptions {
            skip_full_count: skip,
            ..Default::default()
        };
        let events = Events::default();
        let mut silent = Consumer {
            events: events.clone(),
            spectra: 0,
        };
        let quiet = mzxml::transform_with_options(data("MzXMLFile_1.mzXML"), &mut silent, &options)
            .unwrap();
        let (mut logger, _) = make_logger("rec", &Events::default(), &SharedOutput::default());
        let reporting_events = Events::default();
        let mut reporting = Consumer {
            events: reporting_events.clone(),
            spectra: 0,
        };
        let loud = mzxml::transform_with_progress(
            data("MzXMLFile_1.mzXML"),
            &mut reporting,
            &options,
            &mut logger,
        )
        .unwrap();
        assert_eq!(quiet, loud);
        assert_eq!(*events.lock().unwrap(), *reporting_events.lock().unwrap());
        let mut a = MSExperiment::default();
        let mut b = MSExperiment::default();
        mzxml::transform_into_with_options(
            data("MzXMLFile_1.mzXML"),
            &mut silent,
            &mut a,
            &options,
        )
        .unwrap();
        let (mut logger, _) = make_logger("rec", &Events::default(), &SharedOutput::default());
        mzxml::transform_into_with_progress(
            data("MzXMLFile_1.mzXML"),
            &mut reporting,
            &mut b,
            &options,
            &mut logger,
        )
        .unwrap();
        assert_eq!(a, b);
    }
}
