// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Progress reporting of the FORMAT readers the source derives from
//! `ProgressLogger`, replayed against executed OpenMS4 Release output.
//!
//! `tests/data/progress_format_readers_release.tsv` was written by the oracle
//! in `../oracle/progress-format-readers` (see
//! `tests/data/progress_format_readers_provenance.json`): the C++ Release build
//! ran each reader twice on the inputs named below, on a fresh file object,
//! once with `setLogType(GUI)` and a GUI factory making recording backends,
//! which captured every call that reached any backend (the file's own and the
//! fresh ones its handlers make from its log type) with its arguments and
//! nesting depth, and once with `setLogType(CMD)`, whose stdout bytes were
//! captured with only the timing and throughput texts masked. The driver
//! replaced `time()` so every call is a new second and the whole-second
//! throttle never suppresses a set; the clocks here do the same.
//!
//! The replay installs the same kinds of backend through the port's GUI
//! factory, so the backends the port creates from a logger's type (for
//! featureXML, consensusXML and the mzML document section) are observed too.
//! No expected value in this file is derived from Rust output; where the port
//! differs from the Release build, the difference is stated at
//! [`Divergence`] and checked against the captured calls.

#![cfg(all(
    feature = "mzml",
    feature = "featurexml",
    feature = "consensusxml",
    feature = "idxml"
))]

use openms::concept::progress_logger::{
    CommandProgressLogger, ProgressBackend, ProgressClock, ProgressLogType, ProgressLogger,
    ProgressNesting, ProgressReporter, ProgressTime,
};
use openms::format::PeakFileOptions;
use openms::format::{
    consensusxml, dta2d, featurexml, mascot_generic, ms2, mzdata, mzidentml, mzml, mzxml,
};
use openms::kernel::MSExperiment;
use openms::{Error, Result};
use std::collections::BTreeMap;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicI64, Ordering};
use std::sync::{Arc, Mutex, MutexGuard};

const FIXTURE: &str = include_str!("data/progress_format_readers_release.tsv");

/// The mzData handler's scan counter is process-wide, as the source's static
/// is, so every test here that loads mzData holds this lock.
static MZDATA: Mutex<()> = Mutex::new(());
fn mzdata_lock() -> MutexGuard<'static, ()> {
    MZDATA
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}

// ---------------------------------------------------------------------------
// Inputs: the files the driver read, by the name the driver gave them.

fn data(name: &str) -> PathBuf {
    let path = match name {
        "MascotGenericFile_GNPS.mgf" => "MascotGenericFile_GNPS.mgf".to_string(),
        "MzIdentMLFile_whole.mzid" => "mzidentml_whole.mzid".to_string(),
        "ConsensusXMLFile_1.consensusXML" => "consensusxml/ConsensusXMLFile_1.consensusXML".into(),
        "FeatureXMLFile_1.featureXML" => "featurexml_source_1.featureXML".into(),
        "MzDataFile_1.mzData" => "MzDataFile_1.mzData".into(),
        "MzXMLFile_1.mzXML" => "MzXMLFile_1.mzXML".into(),
        "MzMLFile_1.mzML" => "mzml_load_source_original.mzML".into(),
        other => format!("progress_format_readers/{other}"),
    };
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(path)
}

/// A fresh output directory for one run of one case.
fn output_dir(case: &str, mode: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "openms-progress-format-readers-{}-{case}-{mode}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn file_size(path: &Path) -> Result<String> {
    Ok(std::fs::metadata(path)?.len().to_string())
}

fn sizes(experiment: &MSExperiment) -> Vec<String> {
    vec![
        experiment.spectra.len().to_string(),
        experiment.chromatograms.len().to_string(),
    ]
}

/// Runs one fixture case with progress going to `logger`, returning the port's
/// counterpart of the driver's `R` fields.
fn run_case(case: &str, out: &Path, logger: &mut ProgressLogger) -> Result<Vec<String>> {
    let dta2d_options = dta2d::ReadOptions::default();
    let mgf_options = mascot_generic::ReadOptions::default();
    let peaks = PeakFileOptions::default();
    let limits = mzdata::ReadLimits::default();
    let mzml_load = mzml::LoadOptions::default();
    let mzml_read = mzml::ReadOptions::default();
    match case {
        // ---- DTA2DFile
        "dta2d_load" => Ok(sizes(&dta2d::load_with_progress(
            data("DTA2DFile_test_1.dta2d"),
            &dta2d_options,
            logger,
        )?)),
        "dta2d_load_missing" => Ok(sizes(&dta2d::load_with_progress(
            data("missing.dta2d"),
            &dta2d_options,
            logger,
        )?)),
        "dta2d_load_bad_line" => Ok(sizes(&dta2d::load_with_progress(
            data("dta2d_bad_line.dta2d"),
            &dta2d_options,
            logger,
        )?)),
        "dta2d_store" | "dta2d_store_tic" | "dta2d_store_unwritable" => {
            let map = dta2d::load(data("DTA2DFile_test_1.dta2d"))?;
            let path = match case {
                "dta2d_store_unwritable" => out.join("missing_directory/unwritable.dta2d"),
                _ => out.join(format!("{case}.dta2d")),
            };
            let options = dta2d::WriteOptions::default();
            if case == "dta2d_store_tic" {
                dta2d::store_tic_with_progress(&path, &map, &options, logger)?;
            } else {
                dta2d::store_with_progress(&path, &map, &options, logger)?;
            }
            Ok(vec![file_size(&path)?])
        }
        // ---- MS2File: the source makes no call; the port has nothing to wire.
        "ms2_load" => Ok(sizes(&ms2::load(data("MS2File_test_spectra.ms2"))?)),
        // ---- MascotGenericFile
        "mgf_load" => Ok(sizes(&mascot_generic::load_with_progress(
            data("MascotGenericFile_GNPS.mgf"),
            &mgf_options,
            logger,
        )?)),
        "mgf_load_two_blocks" => {
            let path = data("mgf_two_blocks.mgf");
            let mut fields = sizes(&mascot_generic::load_with_progress(
                &path,
                &mgf_options,
                logger,
            )?);
            fields.push(file_size(&path)?);
            Ok(fields)
        }
        "mgf_load_missing" => Ok(sizes(&mascot_generic::load_with_progress(
            data("missing.mgf"),
            &mgf_options,
            logger,
        )?)),
        "mgf_store" => {
            let map = mascot_generic::load(data("mgf_two_blocks.mgf"))?;
            let path = out.join("mgf_store.mgf");
            mascot_generic::MascotGenericFile::new()?.store_with_progress(
                &path,
                &map,
                false,
                &mascot_generic::WriteOptions::default(),
                logger,
            )?;
            Ok(vec![file_size(&path)?])
        }
        // ---- MzIdentMLFile: the source makes no call; nothing to wire.
        "mzid_load" => {
            mzidentml::load(data("MzIdentMLFile_whole.mzid"))?;
            Ok(Vec::new())
        }
        "mzid_store" => {
            let document = mzidentml::load(data("MzIdentMLFile_whole.mzid"))?;
            mzidentml::store(out.join("mzid_store.mzid"), &document)?;
            Ok(Vec::new())
        }
        // ---- ConsensusXMLFile
        "consensus_load" => {
            let map = consensusxml::load_with_progress(
                data("ConsensusXMLFile_1.consensusXML"),
                &consensusxml::ReadOptions::default(),
                logger,
            )?;
            Ok(vec![
                map.features.len().to_string(),
                map.column_headers.len().to_string(),
            ])
        }
        "consensus_load_truncated" => {
            let map = consensusxml::load_with_progress(
                data("truncated.consensusXML"),
                &consensusxml::ReadOptions::default(),
                logger,
            )?;
            Ok(vec![map.features.len().to_string()])
        }
        "consensus_store" => {
            let map = consensusxml::load(data("ConsensusXMLFile_1.consensusXML"))?;
            consensusxml::store_with_progress(
                out.join("consensus_store.consensusXML"),
                &map,
                &consensusxml::WriteOptions::default(),
                logger,
            )?;
            Ok(vec![
                map.protein_identifications.len().to_string(),
                map.column_headers.len().to_string(),
                map.features.len().to_string(),
            ])
        }
        // ---- FeatureXMLFile
        "featurexml_load" | "featurexml_load_truncated" => {
            let name = if case == "featurexml_load" {
                "FeatureXMLFile_1.featureXML"
            } else {
                "truncated.featureXML"
            };
            let map = featurexml::load_with_progress(
                data(name),
                &featurexml::ReadOptions::default(),
                logger,
            )?;
            Ok(vec![map.features.len().to_string()])
        }
        // `loadSize` stops at `<featureList>` before its start: no call.
        "featurexml_load_size" => Ok(vec![
            featurexml::load_size(
                data("FeatureXMLFile_1.featureXML"),
                &featurexml::ReadOptions::default(),
            )?
            .to_string(),
        ]),
        "featurexml_store" => {
            let map = featurexml::load(data("FeatureXMLFile_1.featureXML"))?;
            featurexml::store_with_progress(
                out.join("featurexml_store.featureXML"),
                &map,
                &featurexml::WriteOptions::default(),
                logger,
            )?;
            Ok(vec![map.features.len().to_string()])
        }
        // ---- MzDataFile
        "mzdata_load" => Ok(sizes(
            &mzdata::load_with_progress(data("MzDataFile_1.mzData"), &peaks, &limits, logger)?
                .experiment,
        )),
        "mzdata_store" => {
            let map = mzdata::load(data("MzDataFile_1.mzData"))?;
            // The source writes what mzData cannot hold without refusing.
            mzdata::store_with_progress(
                out.join("mzdata_store.mzData"),
                &map,
                &mzdata::WriteOptions::source(),
                logger,
            )?;
            Ok(vec![map.spectra.len().to_string()])
        }
        "mzdata_static_counter" => {
            let first =
                match mzdata::load_with_progress(data("truncated.mzData"), &peaks, &limits, logger)
                {
                    Ok(_) => "loaded".to_string(),
                    Err(_) => "Parse Error".to_string(),
                };
            // A second file object of the same log type.
            let mut second = logger.clone();
            let loaded = mzdata::load_with_progress(
                data("MzDataFile_1.mzData"),
                &peaks,
                &limits,
                &mut second,
            )?;
            let mut fields = vec![first];
            fields.extend(sizes(&loaded.experiment));
            Ok(fields)
        }
        // ---- MzXMLFile
        "mzxml_load" | "mzxml_load_truncated" => {
            let name = if case == "mzxml_load" {
                "MzXMLFile_1.mzXML"
            } else {
                "truncated.mzXML"
            };
            Ok(sizes(&mzxml::load_with_progress(
                data(name),
                &mzxml::ReadOptions::default(),
                logger,
            )?))
        }
        "mzxml_store" => {
            let map = mzxml::load(data("MzXMLFile_1.mzXML"))?;
            mzxml::store_with_progress(
                out.join("mzxml_store.mzXML"),
                &map,
                &mzxml::WriteOptions::default(),
                logger,
            )?;
            Ok(vec![map.spectra.len().to_string()])
        }
        // ---- MzMLFile
        "mzml_load" => {
            let path = data("MzMLFile_1.mzML");
            let mut fields = sizes(&mzml::load_with_progress(
                &path, &mzml_load, &mzml_read, logger,
            )?);
            fields.push(file_size(&path)?);
            Ok(fields)
        }
        "mzml_load_truncated" => Ok(sizes(&mzml::load_with_progress(
            data("truncated.mzML"),
            &mzml_load,
            &mzml_read,
            logger,
        )?)),
        // The same file object after a failed load.
        "mzml_reuse_after_failure" => {
            let first = match mzml::load_with_progress(
                data("truncated.mzML"),
                &mzml_load,
                &mzml_read,
                logger,
            ) {
                Ok(_) => "loaded".to_string(),
                Err(_) => "Parse Error".to_string(),
            };
            let map =
                mzml::load_with_progress(data("MzMLFile_1.mzML"), &mzml_load, &mzml_read, logger)?;
            let mut fields = vec![first];
            fields.extend(sizes(&map));
            Ok(fields)
        }
        "mzml_load_skip_chromatograms" => {
            let mut load = mzml::LoadOptions::default();
            load.scientific.skip_chromatograms = true;
            Ok(sizes(&mzml::load_with_progress(
                data("MzMLFile_1.mzML"),
                &load,
                &mzml_read,
                logger,
            )?))
        }
        "mzml_store" => {
            let map = mzml::load(data("MzMLFile_1.mzML"))?;
            let path = out.join("mzml_store.mzML");
            mzml::store_with_progress(&path, &map, &mzml::WriteOptions::default(), logger)?;
            let mut fields = sizes(&map);
            fields.push(file_size(&path)?);
            Ok(fields)
        }
        other => panic!("fixture case {other:?} has no Rust counterpart"),
    }
}

// ---------------------------------------------------------------------------
// Loggers.

/// A clock whose every read is a new whole second, as the driver's `time()`.
/// The timer advances with it, so an end with a byte count has a finite rate;
/// the timing and rate texts are masked either way.
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

#[derive(Clone, Default)]
struct SharedOutput(Arc<Mutex<Vec<u8>>>);
impl SharedOutput {
    fn text(&self) -> String {
        String::from_utf8(self.0.lock().unwrap().clone()).unwrap()
    }
}
impl Write for SharedOutput {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// A logger of type GUI whose factory makes a recording backend, as the
/// driver's `make_gui_progress_logger`, with its own nesting context.
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

/// A logger whose every backend is a command backend writing to one buffer,
/// the counterpart of `setLogType(CMD)` with stdout captured.
fn command_logger() -> (ProgressLogger, ProgressNesting, SharedOutput) {
    let nesting = ProgressNesting::default();
    let output = SharedOutput::default();
    let clock = new_second_clock();
    let mut logger = ProgressLogger::with_clock_and_nesting(clock.clone(), nesting.clone());
    let writer = output.clone();
    logger.set_gui_factory(Arc::new(move || {
        Box::new(CommandProgressLogger::with_clock(
            writer.clone(),
            clock.clone(),
        ))
    }));
    logger.set_log_type(ProgressLogType::Gui);
    (logger, nesting, output)
}

/// Masks the timing texts, and the throughput when there is one, of every
/// summary line, as `extract.py` does.
fn mask_timing(text: &str) -> String {
    const HEAD: &str = "-- done [took ";
    let mut out = String::new();
    let mut rest = text;
    while let Some(start) = rest.find(HEAD) {
        let after = &rest[start + HEAD.len()..];
        let end = after.find("] -- ").expect("summary line end");
        let body = &after[..end];
        assert!(
            body.contains(" (CPU), ") && body.contains(" (Wall)"),
            "{body}"
        );
        out.push_str(&rest[..start]);
        out.push_str("-- done [took <TIME> (CPU), <TIME> (Wall)");
        if body.contains(" @ ") {
            assert!(body.ends_with("/s"), "{body}");
            out.push_str(" @ <RATE>/s");
        }
        out.push_str("] -- ");
        rest = &after[end + "] -- ".len()..];
    }
    out.push_str(rest);
    out
}

fn unescape(field: &str) -> String {
    let mut out = String::new();
    let mut chars = field.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('r') => out.push('\r'),
            Some('t') => out.push('\t'),
            Some('\\') => out.push('\\'),
            other => panic!("bad escape {other:?} in fixture field {field:?}"),
        }
    }
    out
}

// ---------------------------------------------------------------------------
// The fixture.

/// One run of one case in one mode, as the Release build recorded it.
#[derive(Default)]
struct Captured {
    events: Vec<String>,
    /// `Ok(R fields)` or `Err((exception name, message))`.
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
            "S" | "V" | "N" | "E" => run.events.push(fields[2..].join("\t")),
            "R" => run.outcome = Some(Ok(rest.iter().map(|s| s.to_string()).collect())),
            "X" => run.outcome = Some(Err((rest[0].into(), rest[1].into()))),
            "D" => run.depth = Some(rest[0].parse().unwrap()),
            "B" => run.output = Some(unescape(rest.first().copied().unwrap_or(""))),
            other => panic!("unknown fixture record {other:?}"),
        }
    }
    (order, runs)
}

/// Where the port's calls differ from the captured ones, and why.
enum Divergence {
    /// The same calls, depth and outcome.
    None,
    /// The port parses the whole document before converting it (see
    /// `consensusxml::load_with_progress`), so a document that is not
    /// well-formed makes no call; the Release build made the calls for the
    /// elements before the truncation, and left its section open.
    ParsesBeforeReporting,
    /// The port's store writes a different document than the source's, so the
    /// byte count of its `endProgress` is the size of the port's own file.
    OwnByteCount,
    /// **A source defect the port corrects (CPP-017):** `setOptions` sets the
    /// handler's `skip_chromatogram_` from `PeakFileOptions::getSkipChromatograms`
    /// (`MzMLHandler.cpp:149`), and the start-element callback returns while it
    /// is set (`:870-873`), so the Release build ignores every element until
    /// `</chromatogramList>` resets it: no section is ever started, each
    /// record end still advances, and each list end ends a section that never
    /// began, which the command backend refuses (`StopWatch.cpp:55`), failing
    /// the load. The port skips only the chromatograms, so it makes the calls of
    /// the ordinary load (`mzml_load`) and loads the same four spectra.
    CorrectsSkipChromatograms,
}
fn divergence(case: &str) -> Divergence {
    match case {
        "consensus_load_truncated" => Divergence::ParsesBeforeReporting,
        "mzml_store" => Divergence::OwnByteCount,
        "mzml_load_skip_chromatograms" => Divergence::CorrectsSkipChromatograms,
        _ => Divergence::None,
    }
}

/// The cases whose reader makes no progress call at all: the source reader
/// has none (`MS2File.h:52`, the mzIdentML handler), the load stops before its
/// section (`FeatureXMLHandler.cpp:306-316`), or the file is missing, which
/// the source reports before the section (`MascotGenericFile.h:76-79`).
const SILENT_CASES: [&str; 5] = [
    "ms2_load",
    "mzid_load",
    "mzid_store",
    "featurexml_load_size",
    "mgf_load_missing",
];

/// The outcome of a run, checked against the Release build's.
fn check_outcome(
    case: &str,
    source: &std::result::Result<Vec<String>, (String, String)>,
    port: &Result<Vec<String>>,
) {
    match (source, port) {
        (Ok(source), Ok(port)) => {
            if let Some(expected) = expected_fields(case, source) {
                assert_eq!(&port[..expected.len()], &expected[..], "{case}");
            }
        }
        (Err(source), Err(error)) => check_error(case, source, error),
        (Err((name, message)), Ok(port))
            if matches!(divergence(case), Divergence::CorrectsSkipChromatograms) =>
        {
            // The command backend's refused end, wrapped by `safeParse_`.
            assert_eq!(name, "Parse Error", "{case}");
            assert!(message.contains("StopWatch.cpp@55"), "{case}: {message}");
            assert_eq!(port, &["4", "0"], "{case}");
        }
        (source, port) => panic!("{case}: source {source:?}, port {port:?}"),
    }
}

/// The port's `R` fields for a case the Release build completed, when they
/// are comparable: a store case's output size is the port's own document.
fn expected_fields(case: &str, source: &[String]) -> Option<Vec<String>> {
    match case {
        // The file sizes of the source's documents.
        "dta2d_store" | "dta2d_store_tic" | "mgf_store" => None,
        // The mzIdentML document is the port's own type; nothing to count.
        "mzid_load" | "mzid_store" => None,
        "mzml_store" => Some(source[..2].to_vec()),
        _ => Some(source.to_vec()),
    }
}

/// The class of the port's error for the source exception a case records.
fn check_error(case: &str, source: &(String, String), error: &Error) {
    let (name, message) = source;
    match name.as_str() {
        "FileNotFound" => {
            assert!(
                matches!(error, Error::Io(e) if e.kind() == io::ErrorKind::NotFound),
                "{case}: {error}"
            );
        }
        "UnableToCreateFile" => assert!(matches!(error, Error::Io(_)), "{case}: {error}"),
        // A second start on the file's running command backend
        // (`StopWatch.cpp:43`), wrapped by `MzMLFile::safeParse_`.
        "Parse Error" if message.contains("StopWatch.cpp@43") => {
            assert_eq!(
                error.to_string(),
                "invalid value: progress timer is already running",
                "{case}"
            );
        }
        "Parse Error" => assert!(
            matches!(error, Error::Parse { .. } | Error::InvalidValue(_)),
            "{case}: {error}"
        ),
        other => panic!("{case}: no error mapping for {other}"),
    }
}

// ---------------------------------------------------------------------------
// Tests.

#[test]
fn fixture_covers_every_reader_in_both_modes() {
    let (order, runs) = fixture();
    assert_eq!(order.len(), 31, "{order:?}");
    for prefix in [
        "dta2d_",
        "ms2_",
        "mgf_",
        "mzid_",
        "consensus_",
        "featurexml_",
        "mzdata_",
        "mzxml_",
        "mzml_",
    ] {
        assert!(order.iter().any(|c| c.starts_with(prefix)), "{prefix}");
    }
    for case in &order {
        let rec = &runs[&(case.clone(), "rec".into())];
        let cmd = &runs[&(case.clone(), "cmd".into())];
        assert!(rec.outcome.is_some() && rec.depth.is_some(), "{case}");
        assert!(cmd.output.is_some() && cmd.outcome.is_some(), "{case}");
        assert!(cmd.events.is_empty(), "{case}");
        let silent = SILENT_CASES.contains(&case.as_str());
        assert_eq!(rec.events.is_empty(), silent, "{case}");
        assert_eq!(cmd.output.as_deref() == Some(""), silent, "{case}");
    }
}

#[test]
fn every_backend_call_matches_the_release_build() {
    let _mzdata = mzdata_lock();
    let (order, runs) = fixture();
    for case in &order {
        let captured = &runs[&(case.clone(), "rec".into())];
        let (mut logger, nesting, events) = recording_logger();
        let out = output_dir(case, "rec");
        let outcome = run_case(case, &out, &mut logger);
        let actual = events.lock().unwrap().clone();
        let mut expected_events = captured.events.clone();
        let mut expected_depth = captured.depth.unwrap();
        match divergence(case) {
            Divergence::None => {}
            Divergence::ParsesBeforeReporting => {
                assert_eq!(expected_depth, 1, "{case}");
                assert!(expected_events[0].starts_with("S\t0\t0\t"), "{case}");
                assert!(!expected_events.last().unwrap().starts_with('E'), "{case}");
                expected_events.clear();
                expected_depth = 0;
            }
            Divergence::CorrectsSkipChromatograms => {
                // No start at all, and one end more than there are sections.
                assert!(
                    expected_events.iter().all(|e| !e.starts_with('S')),
                    "{case}"
                );
                assert_eq!(
                    expected_events
                        .iter()
                        .filter(|e| e.starts_with('E'))
                        .count(),
                    3,
                    "{case}"
                );
                assert_eq!(
                    captured.outcome.as_ref().unwrap(),
                    &Ok(vec!["4".into(), "0".into()])
                );
                expected_events = runs[&("mzml_load".to_string(), "rec".to_string())]
                    .events
                    .clone();
            }
            Divergence::OwnByteCount => {
                let last = expected_events.pop().unwrap();
                let source_bytes = &captured.outcome.as_ref().unwrap().as_ref().unwrap()[2];
                assert_eq!(last, format!("E\t0\t{source_bytes}"), "{case}");
                let port_bytes = outcome.as_ref().unwrap()[2].clone();
                expected_events.push(format!("E\t0\t{port_bytes}"));
            }
        }
        assert_eq!(actual, expected_events, "{case}: backend calls");
        assert_eq!(
            nesting.depth(),
            expected_depth,
            "{case}: nesting afterwards"
        );
        check_outcome(case, captured.outcome.as_ref().unwrap(), &outcome);
        let _ = std::fs::remove_dir_all(out);
    }
}

#[test]
fn command_output_matches_the_release_build() {
    let _mzdata = mzdata_lock();
    let (order, runs) = fixture();
    for case in &order {
        let captured = &runs[&(case.clone(), "cmd".into())];
        let (mut logger, nesting, output) = command_logger();
        let out = output_dir(case, "cmd");
        let outcome = run_case(case, &out, &mut logger);
        let mut expected = captured.output.clone().unwrap();
        let mut expected_depth = captured.depth.unwrap();
        match divergence(case) {
            Divergence::ParsesBeforeReporting => {
                assert!(expected.starts_with("Progress of '"), "{case}");
                assert!(!expected.contains("-- done"), "{case}");
                expected.clear();
                expected_depth = 0;
            }
            Divergence::CorrectsSkipChromatograms => {
                // One dot per record set on a backend that was never started,
                // then the refused end.
                assert_eq!(expected, "....", "{case}");
                let ordinary = &runs[&("mzml_load".to_string(), "cmd".to_string())];
                expected = ordinary.output.clone().unwrap();
                expected_depth = ordinary.depth.unwrap();
            }
            _ => {}
        }
        assert_eq!(mask_timing(&output.text()), expected, "{case}: stdout");
        assert_eq!(
            nesting.depth(),
            expected_depth,
            "{case}: nesting afterwards"
        );
        check_outcome(case, captured.outcome.as_ref().unwrap(), &outcome);
        let _ = std::fs::remove_dir_all(out);
    }
}

/// The mzML store's byte count is that of the document the port wrote: the
/// size of the file, which the Release build's `os.tellp()` also is for its
/// own, larger document.
#[test]
fn the_mzml_store_reports_the_bytes_it_wrote() {
    let map = mzml::load(data("MzMLFile_1.mzML")).unwrap();
    let (mut logger, _, events) = recording_logger();
    let out = output_dir("mzml_store_bytes", "rec");
    let path = out.join("stored.mzML");
    mzml::store_with_progress(&path, &map, &mzml::WriteOptions::default(), &mut logger).unwrap();
    let size = std::fs::metadata(&path).unwrap().len();
    assert!(size > 0);
    assert_eq!(
        events.lock().unwrap().last().unwrap(),
        &format!("E\t0\t{size}")
    );
    std::fs::remove_dir_all(out).unwrap();
}

/// Reporting progress changes no result: each `*_with_progress` entry point
/// returns what its silent counterpart returns, and writes the same bytes.
#[test]
fn progress_changes_no_result() {
    let _mzdata = mzdata_lock();
    let out = output_dir("no_result", "both");
    let command = || command_logger().0;

    let path = data("DTA2DFile_test_1.dta2d");
    let options = dta2d::ReadOptions::default();
    let silent = dta2d::load_with_options(&path, &options).unwrap();
    assert_eq!(
        dta2d::load_with_progress(&path, &options, &mut command()).unwrap(),
        silent
    );
    let limits = dta2d::WriteOptions::default();
    dta2d::store_with_options(out.join("a.dta2d"), &silent, &limits).unwrap();
    dta2d::store_with_progress(out.join("b.dta2d"), &silent, &limits, &mut command()).unwrap();
    dta2d::store_tic_with_options(out.join("c.dta2d"), &silent, &limits).unwrap();
    dta2d::store_tic_with_progress(out.join("d.dta2d"), &silent, &limits, &mut command()).unwrap();
    same_bytes(&out, "a.dta2d", "b.dta2d");
    same_bytes(&out, "c.dta2d", "d.dta2d");

    for name in ["MascotGenericFile_GNPS.mgf", "mgf_two_blocks.mgf"] {
        let options = mascot_generic::ReadOptions::default();
        let silent = mascot_generic::load_with_options(data(name), &options).unwrap();
        assert_eq!(
            mascot_generic::load_with_progress(data(name), &options, &mut command()).unwrap(),
            silent
        );
        let limits = mascot_generic::WriteOptions::default();
        let mut file = mascot_generic::MascotGenericFile::new().unwrap();
        let a = file
            .store_with_options(out.join("a.mgf"), &silent, true, &limits)
            .unwrap();
        let b = file
            .store_with_progress(out.join("b.mgf"), &silent, true, &limits, &mut command())
            .unwrap();
        assert_eq!(a, b);
        same_bytes(&out, "a.mgf", "b.mgf");
    }

    let options = consensusxml::ReadOptions::default();
    let silent =
        consensusxml::load_with_options(data("ConsensusXMLFile_1.consensusXML"), &options).unwrap();
    let reported = consensusxml::load_with_progress(
        data("ConsensusXMLFile_1.consensusXML"),
        &options,
        &command(),
    )
    .unwrap();
    assert!(without_run_counters(&silent).contains(&format!("_{}_N", std::process::id())));
    assert_eq!(
        without_run_counters(&reported),
        without_run_counters(&silent)
    );
    let limits = consensusxml::WriteOptions::default();
    consensusxml::store_with_options(out.join("a.consensusXML"), &silent, &limits).unwrap();
    consensusxml::store_with_progress(out.join("b.consensusXML"), &silent, &limits, &command())
        .unwrap();
    same_bytes(&out, "a.consensusXML", "b.consensusXML");

    let options = featurexml::ReadOptions::default();
    let silent =
        featurexml::load_with_options(data("FeatureXMLFile_1.featureXML"), &options).unwrap();
    let reported =
        featurexml::load_with_progress(data("FeatureXMLFile_1.featureXML"), &options, &command())
            .unwrap();
    assert_eq!(
        without_run_counters(&reported),
        without_run_counters(&silent)
    );
    let limits = featurexml::WriteOptions::default();
    featurexml::store_with_options(out.join("a.featureXML"), &silent, &limits).unwrap();
    featurexml::store_with_progress(out.join("b.featureXML"), &silent, &limits, &command())
        .unwrap();
    same_bytes(&out, "a.featureXML", "b.featureXML");

    let peaks = PeakFileOptions::default();
    let read_limits = mzdata::ReadLimits::default();
    let silent =
        mzdata::load_with_options(data("MzDataFile_1.mzData"), &peaks, &read_limits).unwrap();
    let reported = mzdata::load_with_progress(
        data("MzDataFile_1.mzData"),
        &peaks,
        &read_limits,
        &mut command(),
    )
    .unwrap();
    assert_eq!(reported.experiment, silent.experiment);
    assert_eq!(reported.report, silent.report);
    let limits = mzdata::WriteOptions::source();
    let a = mzdata::store_report_with_options(out.join("a.mzData"), &silent.experiment, &limits)
        .unwrap();
    let b = mzdata::store_with_progress(
        out.join("b.mzData"),
        &silent.experiment,
        &limits,
        &mut command(),
    )
    .unwrap();
    assert_eq!(a, b);
    same_bytes(&out, "a.mzData", "b.mzData");

    let options = mzxml::ReadOptions::default();
    let silent = mzxml::load_with_options(data("MzXMLFile_1.mzXML"), &options).unwrap();
    assert_eq!(
        mzxml::load_with_progress(data("MzXMLFile_1.mzXML"), &options, &mut command()).unwrap(),
        silent
    );
    let limits = mzxml::WriteOptions::default();
    mzxml::store_with_options(out.join("a.mzXML"), &silent, &limits).unwrap();
    mzxml::store_with_progress(out.join("b.mzXML"), &silent, &limits, &mut command()).unwrap();
    same_bytes(&out, "a.mzXML", "b.mzXML");

    let (load, read) = (mzml::LoadOptions::default(), mzml::ReadOptions::default());
    let silent = mzml::load_with_options(data("MzMLFile_1.mzML"), &load, &read).unwrap();
    assert_eq!(
        mzml::load_with_progress(data("MzMLFile_1.mzML"), &load, &read, &mut command()).unwrap(),
        silent
    );
    let limits = mzml::WriteOptions::default();
    mzml::store_with_options(out.join("a.mzML"), &silent, &limits).unwrap();
    mzml::store_with_progress(out.join("b.mzML"), &silent, &limits, &mut command()).unwrap();
    same_bytes(&out, "a.mzML", "b.mzML");
    std::fs::remove_dir_all(out).unwrap();
}

/// `value`'s debug text with the counter of every identification-run
/// identifier the reader generated removed: `map_xml::read_run` names a run
/// `<engine>_<date>_<pid>_<n>` with a process-wide `n`, as the source's
/// `UniqueIdGenerator` does, so two loads of one file differ there and only
/// there.
fn without_run_counters(value: &impl std::fmt::Debug) -> String {
    let text = format!("{value:?}");
    let marker = format!("_{}_", std::process::id());
    let mut pieces = text.split(marker.as_str());
    let mut out = pieces.next().unwrap_or_default().to_string();
    for piece in pieces {
        out.push_str(&marker);
        out.push('N');
        out.push_str(piece.trim_start_matches(|c: char| c.is_ascii_digit()));
    }
    out
}

fn same_bytes(dir: &Path, a: &str, b: &str) {
    let (x, y) = (
        std::fs::read(dir.join(a)).unwrap(),
        std::fs::read(dir.join(b)).unwrap(),
    );
    assert!(!x.is_empty(), "{a}");
    assert!(x == y, "{a} and {b} differ");
}

/// Reporting progress changes no error: on inputs each reader refuses, the
/// `*_with_progress` entry point returns the silent entry point's error.
#[test]
fn progress_changes_no_error() {
    let _mzdata = mzdata_lock();
    let out = output_dir("no_error", "both");
    let command = || command_logger().0;
    let same = |silent: Error, reported: Error| {
        assert_eq!(reported.to_string(), silent.to_string());
    };
    same(
        dta2d::load_with_options(data("dta2d_bad_line.dta2d"), &Default::default()).unwrap_err(),
        dta2d::load_with_progress(
            data("dta2d_bad_line.dta2d"),
            &Default::default(),
            &mut command(),
        )
        .unwrap_err(),
    );
    same(
        mascot_generic::load_with_options(data("missing.mgf"), &Default::default()).unwrap_err(),
        mascot_generic::load_with_progress(
            data("missing.mgf"),
            &Default::default(),
            &mut command(),
        )
        .unwrap_err(),
    );
    same(
        consensusxml::load_with_options(data("truncated.consensusXML"), &Default::default())
            .unwrap_err(),
        consensusxml::load_with_progress(
            data("truncated.consensusXML"),
            &Default::default(),
            &command(),
        )
        .unwrap_err(),
    );
    same(
        featurexml::load_with_options(data("truncated.featureXML"), &Default::default())
            .unwrap_err(),
        featurexml::load_with_progress(
            data("truncated.featureXML"),
            &Default::default(),
            &command(),
        )
        .unwrap_err(),
    );
    let (peaks, limits) = (PeakFileOptions::default(), mzdata::ReadLimits::default());
    let silent = mzdata::load_with_options(data("truncated.mzData"), &peaks, &limits).unwrap_err();
    let reported =
        mzdata::load_with_progress(data("truncated.mzData"), &peaks, &limits, &mut command())
            .unwrap_err();
    same(silent, reported);
    // Both loads left the process-wide scan counter raised, as the source's
    // static is; a complete load resets it for the tests after this one.
    mzdata::load(data("MzDataFile_1.mzData")).unwrap();
    same(
        mzxml::load_with_options(data("truncated.mzXML"), &Default::default()).unwrap_err(),
        mzxml::load_with_progress(data("truncated.mzXML"), &Default::default(), &mut command())
            .unwrap_err(),
    );
    let (load, read) = (mzml::LoadOptions::default(), mzml::ReadOptions::default());
    same(
        mzml::load_with_options(data("truncated.mzML"), &load, &read).unwrap_err(),
        mzml::load_with_progress(data("truncated.mzML"), &load, &read, &mut command()).unwrap_err(),
    );
    // A store refused before any byte, and a destination that cannot be
    // created: the refusal wins over the destination in both entry points.
    let mut refused = featurexml::load(data("FeatureXMLFile_1.featureXML")).unwrap();
    let duplicate = refused.features[0].unique_id;
    refused.features[1].unique_id = duplicate;
    let nowhere = out.join("missing_directory/refused.featureXML");
    let limits = featurexml::WriteOptions::default();
    same(
        featurexml::store_with_options(&nowhere, &refused, &limits).unwrap_err(),
        featurexml::store_with_progress(&nowhere, &refused, &limits, &command()).unwrap_err(),
    );
    let (silent, reported) = (
        featurexml::store_with_options(out.join("refused.featureXML"), &refused, &limits)
            .unwrap_err(),
        featurexml::store_with_progress(
            out.join("refused.featureXML"),
            &refused,
            &limits,
            &command(),
        )
        .unwrap_err(),
    );
    same(silent, reported);
    let intact = featurexml::load(data("FeatureXMLFile_1.featureXML")).unwrap();
    let silent = featurexml::store_with_options(&nowhere, &intact, &limits).unwrap_err();
    let reported =
        featurexml::store_with_progress(&nowhere, &intact, &limits, &command()).unwrap_err();
    assert!(matches!(silent, Error::Io(_)), "{silent}");
    same(silent, reported);
    // The refused store left no file and no temporary file behind.
    assert_eq!(std::fs::read_dir(&out).unwrap().count(), 0);
    std::fs::remove_dir_all(out).unwrap();
}

/// A destination that cannot be created makes no call, as the source's
/// `XMLFile::save_` opens the file before its handler reports anything.
#[test]
fn an_uncreatable_destination_makes_no_call() {
    let _mzdata = mzdata_lock();
    let out = output_dir("uncreatable", "rec");
    let nowhere = out.join("missing_directory");
    let (mut logger, nesting, events) = recording_logger();
    let feature_map = featurexml::load(data("FeatureXMLFile_1.featureXML")).unwrap();
    assert!(
        featurexml::store_with_progress(
            nowhere.join("x.featureXML"),
            &feature_map,
            &Default::default(),
            &logger
        )
        .is_err()
    );
    let consensus = consensusxml::load(data("ConsensusXMLFile_1.consensusXML")).unwrap();
    assert!(
        consensusxml::store_with_progress(
            nowhere.join("x.consensusXML"),
            &consensus,
            &Default::default(),
            &logger
        )
        .is_err()
    );
    let peaks = mzdata::load(data("MzDataFile_1.mzData")).unwrap();
    assert!(
        mzdata::store_with_progress(
            nowhere.join("x.mzData"),
            &peaks,
            &mzdata::WriteOptions::source(),
            &mut logger
        )
        .is_err()
    );
    let scans = mzxml::load(data("MzXMLFile_1.mzXML")).unwrap();
    assert!(
        mzxml::store_with_progress(
            nowhere.join("x.mzXML"),
            &scans,
            &Default::default(),
            &mut logger
        )
        .is_err()
    );
    let spectra = mzml::load(data("MzMLFile_1.mzML")).unwrap();
    assert!(
        mzml::store_with_progress(
            nowhere.join("x.mzML"),
            &spectra,
            &Default::default(),
            &mut logger
        )
        .is_err()
    );
    assert_eq!(*events.lock().unwrap(), Vec::<String>::new());
    assert_eq!(nesting.depth(), 0);
    std::fs::remove_dir_all(out).unwrap();
}

/// The source's featureXML and consensusXML files hand their handler only
/// their log type, and `MzMLHandler` keeps a copy of the file's logger for its
/// document section; a backend installed on the logger with `set_logger`
/// therefore sees only the mzML list sections. Independent of Rust output:
/// the expected calls follow from `ConsensusXMLFile.cpp:90-92`,
/// `FeatureXMLFile.cpp:54` and `MzMLHandler.cpp:135`.
#[test]
fn an_installed_backend_sees_only_the_calls_made_on_the_logger_itself() {
    let events = Arc::new(Mutex::new(Vec::new()));
    let nesting = ProgressNesting::default();
    let mut logger = ProgressLogger::with_clock_and_nesting(new_second_clock(), nesting.clone());
    logger.set_logger(Box::new(Recorder {
        events: events.clone(),
        current: 0,
    }));
    featurexml::load_with_progress(
        data("FeatureXMLFile_1.featureXML"),
        &Default::default(),
        &logger,
    )
    .unwrap();
    consensusxml::load_with_progress(
        data("ConsensusXMLFile_1.consensusXML"),
        &Default::default(),
        &logger,
    )
    .unwrap();
    assert!(events.lock().unwrap().is_empty());
    mzml::load_with_progress(
        data("MzMLFile_1.mzML"),
        &Default::default(),
        &Default::default(),
        &mut logger,
    )
    .unwrap();
    // The lists at depth 1, below the document section the copy reported.
    assert_eq!(
        *events.lock().unwrap(),
        [
            "S\t0\t4\tloading spectra list\t1",
            "N\t1",
            "V\t1\t2",
            "N\t2",
            "V\t2\t2",
            "N\t3",
            "V\t3\t2",
            "N\t4",
            "V\t4\t2",
            "E\t1\t0",
            "S\t0\t2\tloading chromatogram list\t1",
            "N\t1",
            "V\t1\t2",
            "N\t2",
            "V\t2\t2",
            "E\t1\t0",
        ]
    );
    assert_eq!(nesting.depth(), 0);
}

/// A metadata-only mzML load starts the document section and stops at the
/// first list, leaving the section open, as the source's `EndParsingSoftly`
/// at `<spectrumList>` (`MzMLHandler.cpp:960-964`) leaves `pg_outer`'s.
#[test]
fn a_metadata_only_mzml_load_leaves_the_document_section_open() {
    let (mut logger, nesting, events) = recording_logger();
    let mut load = mzml::LoadOptions::default();
    load.scientific.metadata_only = true;
    mzml::load_with_progress(
        data("MzMLFile_1.mzML"),
        &load,
        &Default::default(),
        &mut logger,
    )
    .unwrap();
    assert_eq!(*events.lock().unwrap(), ["S\t0\t1\tloading mzML\t0"]);
    assert_eq!(nesting.depth(), 1);
}

/// The three reporter calls the readers added. Independent: the expectations
/// follow from the documented contract of `ProgressReporter` and
/// `ProgressLogger`, not from Rust output.
#[test]
fn the_reporter_counts_advances_and_ends_with_a_byte_count() {
    let mut silent = ProgressReporter::silent();
    silent.start_count(usize::MAX, "never shown").unwrap();
    silent.next_progress().unwrap();
    silent.end_with_bytes(u64::MAX).unwrap();

    let (mut logger, nesting, events) = recording_logger();
    let mut reporter = ProgressReporter::new(Some(&mut logger));
    reporter.start_count(3, "label").unwrap();
    reporter.next_progress().unwrap();
    reporter.next_progress().unwrap();
    reporter.end_with_bytes(42).unwrap();
    if let Some(over) = usize::try_from(i64::MAX)
        .ok()
        .and_then(|max| max.checked_add(1))
    {
        assert!(matches!(
            reporter.start_count(over, "too many"),
            Err(Error::InvalidValue(_))
        ));
    }
    assert_eq!(
        *events.lock().unwrap(),
        [
            "S\t0\t3\tlabel\t0",
            "N\t1",
            "V\t1\t1",
            "N\t2",
            "V\t2\t1",
            "E\t0\t42"
        ]
    );
    assert_eq!(nesting.depth(), 0);
}
