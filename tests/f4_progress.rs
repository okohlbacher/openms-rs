// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Progress reporting of `CachedMzMLHandler` and `GNPSMGFFile`, the two F4
//! classes the source derives from `ProgressLogger`, replayed against executed
//! OpenMS4 Release output.
//!
//! `tests/data/f4_progress_release.tsv` was written by the oracle in
//! `../oracle/F4-progress` (see `tests/data/f4_progress_provenance.json`): the
//! C++ Release build ran each case on a fresh object with `setLogType(GUI)`
//! and a GUI factory making recording backends, which captured every call that
//! reached a backend with its arguments and nesting depth, then the outcome
//! and the static nesting depth afterwards. The driver replaced `time()` so
//! every call is a new second and the whole-second throttle never suppresses a
//! set; the clock here does the same. The files the Release cases wrote (the
//! cache dumps and the two MGF files, the second with the `merged_spectra`
//! defaults, which no TOPP test covers) are held in `tests/data/f4_progress/`
//! and compared byte for byte. No expected value in this file is derived from
//! Rust output.

#![cfg(all(feature = "mzml", feature = "consensusxml"))]

use openms::concept::progress_logger::{
    ProgressBackend, ProgressClock, ProgressLogType, ProgressLogger, ProgressNesting, ProgressTime,
};
use openms::format::cached_mzml_handler::CachedMzMLHandler;
use openms::format::gnps_mgf::{GNPSMGFFile, OutputType};
use openms::kernel::MSExperiment;
use openms::system::file::TempDir;
use openms::{Error, Result};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicI64, Ordering};
use std::sync::{Arc, Mutex};

const FIXTURE: &str = include_str!("data/f4_progress_release.tsv");

fn data(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(path)
}

fn input(name: &str) -> PathBuf {
    match name {
        "MzMLFile_1.mzML" => data("mzml_validator/MzMLFile_1.mzML"),
        "GNPSExport_cons1.consensusXML" => data("gnps_tables/GNPSExport_cons1.consensusXML"),
        other => data(&format!("gnps_mgf/{other}")),
    }
}

/// One oracle case: the backend calls, the outcome row and the depth after.
#[derive(Debug, Default)]
struct Captured {
    events: Vec<String>,
    outcome: Option<Vec<String>>,
    depth: Option<usize>,
}

fn fixture() -> (Vec<String>, BTreeMap<String, Captured>) {
    let mut order = Vec::new();
    let mut runs: BTreeMap<String, Captured> = BTreeMap::new();
    for line in FIXTURE.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let (case, rest) = line.split_once('\t').unwrap();
        if !runs.contains_key(case) {
            order.push(case.to_string());
        }
        let captured = runs.entry(case.to_string()).or_default();
        match &rest[..1] {
            "S" | "V" | "N" | "E" => captured.events.push(rest.to_string()),
            "R" | "X" => captured.outcome = Some(rest.split('\t').map(String::from).collect()),
            "D" => captured.depth = Some(rest[2..].parse().unwrap()),
            other => panic!("unknown record {other}"),
        }
    }
    (order, runs)
}

// ---------------------------------------------------------------------------
// The recording logger (as `tests/progress_format_readers.rs`).

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

struct Recorder {
    events: Arc<Mutex<Vec<String>>>,
    current: i64,
}
impl Recorder {
    fn push(&self, event: String) {
        self.events.lock().unwrap().push(event);
    }
}
impl ProgressBackend for Recorder {
    fn start_progress(&mut self, begin: i64, end: i64, label: &str, depth: usize) -> Result<()> {
        self.current = begin;
        self.push(format!("S\t{begin}\t{end}\t{label}\t{depth}"));
        Ok(())
    }
    fn set_progress(&mut self, value: i64, depth: usize) -> Result<()> {
        self.push(format!("V\t{value}\t{depth}"));
        Ok(())
    }
    fn next_progress(&mut self) -> Result<i64> {
        self.current += 1;
        self.push(format!("N\t{}", self.current));
        Ok(self.current)
    }
    fn end_progress(&mut self, depth: usize, bytes_processed: u64) -> Result<()> {
        self.push(format!("E\t{depth}\t{bytes_processed}"));
        Ok(())
    }
}

fn recording_logger() -> (ProgressLogger, ProgressNesting, Arc<Mutex<Vec<String>>>) {
    let nesting = ProgressNesting::default();
    let events = Arc::new(Mutex::new(Vec::new()));
    let mut logger = ProgressLogger::with_clock_and_nesting(new_second_clock(), nesting.clone());
    let sink = events.clone();
    logger.set_gui_factory(Arc::new(move || {
        Box::new(Recorder {
            events: sink.clone(),
            current: 0,
        })
    }));
    logger.set_log_type(ProgressLogType::Gui);
    (logger, nesting, events)
}

// ---------------------------------------------------------------------------
// The cases, as the driver runs them.

fn sizes(exp: &MSExperiment) -> Vec<String> {
    vec![
        exp.spectra.len().to_string(),
        exp.chromatograms.len().to_string(),
    ]
}

fn file_size(path: &Path) -> Result<String> {
    Ok(std::fs::metadata(path)?.len().to_string())
}

/// Runs one case with progress going to `logger`; `out` holds the outputs of
/// earlier cases, as the driver's output directory does.
fn run_case(case: &str, out: &Path, logger: &mut ProgressLogger) -> Result<Vec<String>> {
    let dump = out.join("MzMLFile_1.cached");
    let empty = out.join("empty.cached");
    let mzmls = [
        input("GNPSExport_mzml1.mzML"),
        input("GNPSExport_mzml2.mzML"),
    ];
    match case {
        "cached_write" => {
            let exp = openms::format::mzml::load(input("MzMLFile_1.mzML"))?;
            CachedMzMLHandler::new().write_memdump_with_progress(&exp, &dump, logger)?;
            let mut fields = sizes(&exp);
            fields.push(file_size(&dump)?);
            Ok(fields)
        }
        "cached_write_empty" => {
            CachedMzMLHandler::new().write_memdump_with_progress(
                &MSExperiment::new(),
                &empty,
                logger,
            )?;
            Ok(vec![file_size(&empty)?])
        }
        "cached_read" => Ok(sizes(
            &CachedMzMLHandler::new().read_memdump_with_progress(&dump, logger)?,
        )),
        "cached_read_empty" => Ok(sizes(
            &CachedMzMLHandler::new().read_memdump_with_progress(&empty, logger)?,
        )),
        "cached_read_bad_magic" => Ok(sizes(
            &CachedMzMLHandler::new()
                .read_memdump_with_progress(input("MzMLFile_1.mzML"), logger)?,
        )),
        "cached_index" | "cached_index_bad_magic" => {
            let path = if case == "cached_index" {
                dump
            } else {
                input("MzMLFile_1.mzML")
            };
            let mut handler = CachedMzMLHandler::new();
            handler.create_memdump_index_with_progress(path, logger)?;
            Ok(vec![
                handler.spectra_index().len().to_string(),
                handler.chromatogram_index().len().to_string(),
            ])
        }
        "gnps_most_intense" | "gnps_merged" => {
            let mut writer = GNPSMGFFile::new();
            if case == "gnps_merged" {
                writer.output_type = OutputType::MergedSpectra;
            }
            let path = out.join(format!("{case}.mgf"));
            writer.store_with_progress(
                input("GNPSExport_cons1.consensusXML"),
                &mzmls,
                &path,
                logger,
            )?;
            Ok(vec![file_size(&path)?])
        }
        other => panic!("case {other} has no counterpart"),
    }
}

/// The source exception each refusal of the port stands for.
fn check_outcome(case: &str, expected: &[String], actual: &Result<Vec<String>>) {
    match (expected[0].as_str(), actual) {
        ("R", Ok(fields)) => assert_eq!(fields.as_slice(), &expected[1..], "{case}: result"),
        ("X", Err(Error::Parse { message, .. })) => {
            assert_eq!(expected[1], "Parse Error", "{case}");
            assert!(expected[2].contains(message.as_str()), "{case}: {message}");
        }
        (_, other) => panic!("{case}: source {expected:?}, port {other:?}"),
    }
}

#[test]
fn fixture_covers_every_case() {
    let (order, runs) = fixture();
    assert_eq!(
        order,
        [
            "cached_write",
            "cached_write_empty",
            "cached_read",
            "cached_read_empty",
            "cached_read_bad_magic",
            "cached_index",
            "cached_index_bad_magic",
            "gnps_most_intense",
            "gnps_merged",
        ]
    );
    for (case, captured) in &runs {
        assert!(captured.outcome.is_some(), "{case}");
        assert!(captured.depth.is_some(), "{case}");
    }
}

#[test]
fn every_backend_call_matches_the_release_build() {
    let (order, runs) = fixture();
    let out = TempDir::new(false).unwrap();
    for case in &order {
        let captured = &runs[case];
        let (mut logger, nesting, events) = recording_logger();
        let outcome = run_case(case, out.path(), &mut logger);
        let actual = events.lock().unwrap().clone();
        assert_eq!(actual, captured.events, "{case}: backend calls");
        assert_eq!(
            nesting.depth(),
            captured.depth.unwrap(),
            "{case}: nesting afterwards"
        );
        check_outcome(case, captured.outcome.as_ref().unwrap(), &outcome);
    }
    // The files the cases wrote are the Release build's, byte for byte.
    for name in [
        "MzMLFile_1.cached",
        "empty.cached",
        "gnps_most_intense.mgf",
        "gnps_merged.mgf",
    ] {
        let expected = std::fs::read(data(&format!("f4_progress/{name}"))).unwrap();
        let actual = std::fs::read(out.path().join(name)).unwrap();
        assert!(
            actual == expected,
            "{name}: bytes differ from the Release output"
        );
    }
}

#[test]
fn progress_changes_no_output() {
    let out = TempDir::new(false).unwrap();
    let exp = openms::format::mzml::load(input("MzMLFile_1.mzML")).unwrap();
    let handler = CachedMzMLHandler::new();
    let silent = out.path().join("silent.cached");
    let reported = out.path().join("reported.cached");
    handler.write_memdump(&exp, &silent).unwrap();
    let (mut logger, _, _) = recording_logger();
    handler
        .write_memdump_with_progress(&exp, &reported, &mut logger)
        .unwrap();
    assert_eq!(
        std::fs::read(&silent).unwrap(),
        std::fs::read(&reported).unwrap()
    );
    assert_eq!(
        handler.read_memdump(&silent).unwrap(),
        handler
            .read_memdump_with_progress(&silent, &mut logger)
            .unwrap()
    );
    let mut plain = CachedMzMLHandler::new();
    plain.create_memdump_index(&silent).unwrap();
    let mut with = CachedMzMLHandler::new();
    with.create_memdump_index_with_progress(&silent, &mut logger)
        .unwrap();
    assert_eq!(plain, with);

    let mzmls = [
        input("GNPSExport_mzml1.mzML"),
        input("GNPSExport_mzml2.mzML"),
    ];
    let writer = GNPSMGFFile::new();
    let map = openms::format::consensusxml::load(input("GNPSExport_cons1.consensusXML")).unwrap();
    let mut plain = Vec::new();
    writer.write(&map, &mzmls, &mut plain).unwrap();
    let mut reported = Vec::new();
    writer
        .write_with_progress(&map, &mzmls, &mut reported, &mut logger)
        .unwrap();
    assert_eq!(plain, reported);
}
