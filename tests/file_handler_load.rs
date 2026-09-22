// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FileHandler::loadExperiment` with all of its arguments: the log type,
//! `rewrite_source_file`, `compute_hash`, the SRM conversion after an mzML
//! load, the mzXML and mzData dispatch, and `FileHandler::computeFileHash`.
//!
//! Two kinds of expected value, neither derived from Rust output:
//!
//! - the upstream class test `FileHandler_test.cpp` at core `bc9cc12`, cited
//!   per test, and RFC 3174's SHA-1 test vectors;
//! - `tests/data/file_handler_load_release.tsv`, written by the C++ Release
//!   build at core `bc9cc12` on ibminode06 through the oracle in
//!   `../oracle/file-handler-load` (provenance in
//!   `tests/data/file_handler_load_provenance.json`): a fresh `FileHandler`
//!   per case, log type GUI with a recording backend factory, then every
//!   source-file field, the first spectra, every chromatogram, and the
//!   warning-log lines of the case.

#![cfg(feature = "mzml")]

use openms::concept::progress_logger::ProgressLogger;
use openms::format::file_handler::ExperimentLoadOptions;
use openms::format::{FileHandler, FileType, PeakFileOptions};
use openms::kernel::NumericRange;
use openms::metadata::ChecksumType;
use openms::{Error, MSExperiment};
use std::path::{Path, PathBuf};

fn data(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(relative)
}

/// `MzMLFile_1.mzML` of the core class-test data; the port's copy has the
/// SHA-1 the upstream test expects, which is also the check that it is the
/// same file.
fn mzml_file_1() -> PathBuf {
    data("mzml_load_source_original.mzML")
}
fn dta2d_test_1() -> PathBuf {
    data("text_peak_lists/DTA2DFile_test_1.dta2d")
}

fn rewrite(hash: bool) -> ExperimentLoadOptions {
    ExperimentLoadOptions {
        rewrite_source_file: true,
        compute_hash: hash,
        ..Default::default()
    }
}

fn mz_range(min: f64, max: f64) -> PeakFileOptions {
    let mut options = PeakFileOptions::default();
    options.set_mz_range(NumericRange { min, max });
    options
}

fn first_mz(experiment: &MSExperiment, spectrum: usize, count: usize) -> Vec<f64> {
    experiment.spectra[spectrum]
        .peaks
        .iter()
        .take(count)
        .map(|peak| peak.mz)
        .collect()
}

fn assert_similar(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len(), "{actual:?} vs {expected:?}");
    for (a, e) in actual.iter().zip(expected) {
        // TEST_REAL_SIMILAR's default relative tolerance, 1e-5.
        assert!((a - e).abs() <= 1e-5 * e.abs().max(1.0), "{actual:?} vs {expected:?}");
    }
}

#[test]
fn upstream_mzml_load_rewrites_the_source_file_with_its_sha1() {
    // FileHandler_test.cpp:257-261.
    let experiment = FileHandler::load_experiment_with_progress(
        mzml_file_1(),
        &[FileType::MzMl],
        &rewrite(true),
        &ProgressLogger::new(),
    )
    .unwrap();
    assert_eq!(experiment.spectra.len(), 4);
    let sources = &experiment.settings.source_files;
    assert_eq!(sources.len(), 1);
    assert_eq!(sources[0].checksum, "36007593dbca0ba59a1f4fc32fb970f0e8991fa6");
    assert_eq!(sources[0].checksum_type, ChecksumType::Sha1);
}

#[test]
fn upstream_dta2d_load_applies_the_mz_range_and_rewrites_with_or_without_hash() {
    // FileHandler_test.cpp:263-279.
    let plain = FileHandler::load_experiment_with_options(
        dta2d_test_1(),
        &[],
        &PeakFileOptions::default(),
    )
    .unwrap();
    assert_similar(&first_mz(&plain, 0, 3), &[230.02, 430.02, 630.02]);

    let options = ExperimentLoadOptions {
        peaks: mz_range(300.0, 1000.0),
        ..rewrite(true)
    };
    let hashed = FileHandler::load_experiment_with_progress(
        dta2d_test_1(),
        &[FileType::Dta2d],
        &options,
        &ProgressLogger::new(),
    )
    .unwrap();
    assert_similar(&first_mz(&hashed, 0, 2), &[430.02, 630.02]);
    let source = &hashed.settings.source_files[0];
    assert_eq!(source.checksum, "d50d5144cc3805749b9e8d16f3bc8994979d8142");
    assert_eq!(source.checksum_type, ChecksumType::Sha1);

    let unhashed = FileHandler::load_experiment_with_progress(
        dta2d_test_1(),
        &[],
        &rewrite(false),
        &ProgressLogger::new(),
    )
    .unwrap();
    let source = &unhashed.settings.source_files[0];
    assert_eq!(source.checksum, "");
    assert_eq!(source.checksum_type, ChecksumType::Unknown);
}

#[test]
fn upstream_a_type_outside_the_allowed_list_is_refused() {
    // FileHandler_test.cpp:280-281 restricts to SIZE_OF_TYPE, which no file
    // has; the source throws ParseError, the port InvalidValue.
    let result = FileHandler::load_experiment_with_progress(
        dta2d_test_1(),
        &[FileType::Unknown],
        &rewrite(false),
        &ProgressLogger::new(),
    );
    assert!(matches!(result, Err(Error::InvalidValue(_))), "{result:?}");
}

#[test]
fn upstream_mzdata_and_mzxml_are_dispatched_with_the_peak_options() {
    // FileHandler_test.cpp:234-255.
    let load = |name: &str, options: &PeakFileOptions| {
        FileHandler::load_experiment_with_options(data(name), &[], options).unwrap()
    };
    let all = PeakFileOptions::default();
    let from_115 = mz_range(115.0, 1000.0);
    assert_similar(
        &first_mz(&load("MzDataFile_1.mzData", &all), 1, 3),
        &[110.0, 120.0, 130.0],
    );
    assert_similar(
        &first_mz(&load("MzDataFile_1.mzData", &from_115), 1, 2),
        &[120.0, 130.0],
    );
    assert_similar(
        &first_mz(&load("MzXMLFile_1.mzXML", &all), 2, 3),
        &[100.0, 110.0, 120.0],
    );
    assert_similar(
        &first_mz(&load("MzXMLFile_1.mzXML", &from_115), 2, 3),
        &[120.0, 130.0, 140.0],
    );
    for kind in [FileType::MzXml, FileType::MzData] {
        assert!(FileHandler::can_read_experiment(kind));
        assert!(!FileHandler::can_write_experiment(kind));
    }
}

#[test]
fn compute_file_hash_follows_rfc_3174_and_the_upstream_non_ascii_case() {
    let directory = std::env::temp_dir().join(format!(
        "openms-rs-file-hash-t\u{e4}st_\u{fc}-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&directory).unwrap();
    let write = |name: &str, bytes: &[u8]| {
        let path = directory.join(name);
        std::fs::write(&path, bytes).unwrap();
        path
    };
    // RFC 3174 section 7.3, TEST1 ("abc"), and the digest of no input.
    assert_eq!(
        FileHandler::compute_file_hash(write("abc.txt", b"abc")),
        "a9993e364706816aba3e25717850c26c9cd0d89d"
    );
    assert_eq!(
        FileHandler::compute_file_hash(write("empty.txt", b"")),
        "da39a3ee5e6b4b0d3255bfef95601890afd80709"
    );
    // RFC 3174 TEST3, one million 'a', crosses many 8 KiB reads.
    assert_eq!(
        FileHandler::compute_file_hash(write("million.txt", &vec![b'a'; 1_000_000])),
        "34aa973cd4c4daa4f61eeb2bdbad27316534016f"
    );
    // FileHandler_test.cpp:293-319: a non-ASCII path hashes.
    assert!(!FileHandler::compute_file_hash(write("test.txt", b"hello")).is_empty());
    // FileHandler.cpp:771-774: a file that cannot be opened hashes to "".
    assert_eq!(FileHandler::compute_file_hash(directory.join("missing.txt")), "");
    std::fs::remove_dir_all(&directory).unwrap();
}

#[test]
fn compute_hash_without_rewrite_is_refused_before_any_io() {
    // FileHandler.cpp:853 makes this an OPENMS_PRECONDITION, which only a
    // debug build checks. The Release build ignores the flag: its
    // mzml_hash_only case reports exactly what mzml_plain does.
    let rows = |case: &str| {
        release(case)
            .into_iter()
            .map(|row| row.join("\t"))
            .collect::<Vec<_>>()
    };
    assert_eq!(rows("mzml_hash_only"), rows("mzml_plain"));
    assert!(!rows("mzml_plain").is_empty());
    // The port refuses the combination, as the debug build does, and before
    // it touches the file.
    let options = ExperimentLoadOptions {
        compute_hash: true,
        ..Default::default()
    };
    let result = FileHandler::load_experiment_with_progress(
        data("does-not-exist.mzML"),
        &[],
        &options,
        &ProgressLogger::new(),
    );
    assert!(matches!(result, Err(Error::InvalidValue(_))), "{result:?}");
}

// ---------------------------------------------------------------------------
// The Release oracle.

const RELEASE: &str = include_str!("data/file_handler_load_release.tsv");

/// The fixture's rows for `case`, without the case column.
fn release(case: &str) -> Vec<Vec<&'static str>> {
    RELEASE
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|fields| fields[0] == case)
        .map(|fields| fields[1..].to_vec())
        .collect()
}

/// The fixture's rows of one kind for `case`, tab-joined, without the kind.
fn release_kind(case: &str, kind: &str) -> Vec<String> {
    release(case)
        .into_iter()
        .filter(|row| row[0] == kind)
        .map(|row| row[1..].join("\t"))
        .collect()
}

/// One field of a row: text compared exactly, or a number the driver printed
/// with `%.17g`, which parses back to the same `f64`.
#[derive(Debug)]
enum Field {
    Text(String),
    Number(f64),
}

fn text(value: impl ToString) -> Field {
    Field::Text(value.to_string())
}

fn position<T: PartialEq>(all: &[T], value: &T) -> Field {
    text(all.iter().position(|item| item == value).unwrap())
}

/// The driver's `dump`, row for row: the size, every source file, the first
/// twelve spectra and every chromatogram.
fn dump(experiment: &MSExperiment) -> Vec<(&'static str, Vec<Field>)> {
    use openms::metadata::{ChromatogramType, ScanMode};
    let mut rows = vec![(
        "size",
        vec![
            text(experiment.spectra.len()),
            text(experiment.chromatograms.len()),
        ],
    )];
    for (i, s) in experiment.settings.source_files.iter().enumerate() {
        rows.push((
            "sf",
            vec![
                text(i),
                text(&s.name),
                text(&s.path),
                text(&s.file_type),
                text(&s.checksum),
                position(ChecksumType::ALL, &s.checksum_type),
                text(&s.native_id_type),
                text(&s.native_id_type_accession),
                Field::Number(f64::from(s.size_mb)),
            ],
        ));
    }
    for (i, s) in experiment.spectra.iter().take(12).enumerate() {
        let mut row = vec![
            text(i),
            Field::Number(s.rt),
            text(s.ms_level),
            text(s.peaks.len()),
            text(&s.native_id),
            position(ScanMode::ALL, &s.instrument_settings.scan_mode),
            text(s.precursors.len()),
        ];
        row.extend(s.peaks.iter().take(3).map(|p| Field::Number(p.mz)));
        rows.push(("spec", row));
    }
    for (i, c) in experiment.chromatograms.iter().enumerate() {
        let mut row = vec![
            text(i),
            text(&c.native_id),
            position(ChromatogramType::ALL, &c.chromatogram_type),
            Field::Number(c.precursor.mz),
            Field::Number(c.product.mz),
            text(c.peaks.len()),
        ];
        // Each point is `%.17g:%.17g`.
        row.extend(c.peaks.iter().map(|p| {
            Field::Text(format!("{}:{}", p.rt, f64::from(p.intensity)))
        }));
        rows.push(("chrom", row));
    }
    rows
}

/// Compare `experiment` with the case's size, source-file, spectrum and
/// chromatogram rows; `{INPUTS}` stands for `inputs`.
fn assert_experiment(case: &str, inputs: &Path, experiment: &MSExperiment, kinds: &[&str]) {
    let inputs = inputs.to_str().unwrap();
    let expected: Vec<Vec<String>> = release(case)
        .into_iter()
        .filter(|row| kinds.contains(&row[0]))
        .map(|row| {
            row.iter()
                .map(|field| field.replace("{INPUTS}", inputs))
                .collect()
        })
        .collect();
    let actual: Vec<_> = dump(experiment)
        .into_iter()
        .filter(|(kind, _)| kinds.contains(kind))
        .collect();
    assert_eq!(
        expected.len(),
        actual.len(),
        "{case}: {expected:?} vs {actual:?}"
    );
    for (e, (kind, a)) in expected.iter().zip(&actual) {
        let same = e[0] == *kind
            && e.len() == a.len() + 1
            && e[1..].iter().zip(a).all(|(e, a)| match a {
                Field::Text(a) => e == a,
                Field::Number(a) => e.parse::<f64>().is_ok_and(|e| e == *a),
            });
        // Chromatogram points print as %.17g; every point here is small and
        // integral or half-integral, which `{}` prints the same way.
        assert!(same, "{case}: {e:?} vs {kind} {a:?}");
    }
}

const ALL_ROWS: &[&str] = &["size", "sf", "spec", "chrom"];

/// The oracle's inputs, under the names the driver gave them.
struct Inputs(PathBuf);
impl Inputs {
    fn new() -> Self {
        static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let id = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let root = std::env::temp_dir().join(format!(
            "openms-rs-file-handler-load-{}-{id}",
            std::process::id()
        ));
        let inputs = root.join("inputs");
        std::fs::create_dir_all(inputs.join("dir.d")).unwrap();
        for (name, source) in [
            ("MzMLFile_1.mzML", "mzml_load_source_original.mzML"),
            ("MzMLFile_1.mzML.gz", "file_handler_load/MzMLFile_1.mzML.gz"),
            ("DTA2DFile_test_1.dta2d", "text_peak_lists/DTA2DFile_test_1.dta2d"),
            ("MzXMLFile_1.mzXML", "MzXMLFile_1.mzXML"),
            ("MzDataFile_1.mzData", "MzDataFile_1.mzData"),
            ("srm_spectra.mzML", "file_info/inputs/srm_spectra.mzML"),
            ("srm_spectra_mixed.mzML", "file_info/inputs/srm_spectra_mixed.mzML"),
            ("mgf_two_blocks.mgf", "progress_format_readers/mgf_two_blocks.mgf"),
            ("MS2File_test_spectra.ms2", "text_peak_lists/MS2File_test_spectra.ms2"),
            ("dta_extractor_1_output.dta", "dta_extractor_1_output.dta"),
        ] {
            std::fs::copy(data(source), inputs.join(name)).unwrap();
        }
        std::fs::write(inputs.join("empty.txt"), b"").unwrap();
        std::fs::write(inputs.join("hello.txt"), b"hello").unwrap();
        Self(root)
    }
    fn dir(&self) -> PathBuf {
        self.0.join("inputs")
    }
    fn path(&self, name: &str) -> PathBuf {
        self.dir().join(name)
    }
}
impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

type Events = std::sync::Arc<std::sync::Mutex<Vec<String>>>;

/// Records every backend call in the driver's row format.
struct Recorder {
    events: Events,
    current: i64,
}
impl Recorder {
    fn push(&self, event: String) {
        self.events.lock().unwrap().push(event);
    }
}
impl openms::concept::progress_logger::ProgressBackend for Recorder {
    fn start_progress(
        &mut self,
        begin: i64,
        end: i64,
        label: &str,
        depth: usize,
    ) -> openms::Result<()> {
        self.current = begin;
        self.push(format!("S\t{begin}\t{end}\t{label}\t{depth}"));
        Ok(())
    }
    fn set_progress(&mut self, value: i64, depth: usize) -> openms::Result<()> {
        self.push(format!("V\t{value}\t{depth}"));
        Ok(())
    }
    fn next_progress(&mut self) -> openms::Result<i64> {
        self.current += 1;
        self.push(format!("N\t{}", self.current));
        Ok(self.current)
    }
    fn end_progress(&mut self, depth: usize, bytes_processed: u64) -> openms::Result<()> {
        self.push(format!("E\t{depth}\t{bytes_processed}"));
        Ok(())
    }
}

/// A GUI-type logger whose backends record into `events`, on a clock whose
/// every read is a new second and an isolated nesting context: the driver's
/// recording factory and `time()`.
fn recording_logger(events: &Events) -> ProgressLogger {
    use openms::concept::progress_logger::{ProgressLogType, ProgressNesting, ProgressTime};
    let second = std::sync::Arc::new(std::sync::atomic::AtomicI64::new(1_000_000_000));
    let clock = std::sync::Arc::new(move || {
        let now = second.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        Ok(ProgressTime {
            wall_second: now,
            wall_seconds: (now - 1_000_000_000) as f64,
            cpu_seconds: Some(0.0),
        })
    });
    let mut logger = ProgressLogger::with_clock_and_nesting(clock, ProgressNesting::default());
    let events = events.clone();
    logger.set_gui_factory(std::sync::Arc::new(move || {
        Box::new(Recorder {
            events: events.clone(),
            current: 0,
        })
    }));
    logger.set_log_type(ProgressLogType::Gui);
    logger
}

/// Routes this test thread's warning log into a buffer. `take` clears the
/// repetition cache, as the driver does before each case, and returns the
/// lines written since the last call.
struct Warnings(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);
impl Warnings {
    fn capture() -> Self {
        use openms::concept::log_stream::{LogLevel, LogSink, with_thread_local_log};
        #[derive(Clone)]
        struct Buffer(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);
        impl std::io::Write for Buffer {
            fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
                self.0.lock().unwrap().extend_from_slice(bytes);
                Ok(bytes.len())
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }
        let buffer = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let sink = LogSink::new(Buffer(buffer.clone()));
        with_thread_local_log(LogLevel::Warn, |log| {
            log.remove_all_streams()?;
            log.set_color(None);
            log.insert(&sink)
        })
        .unwrap();
        Self(buffer)
    }
    fn take(&self) -> Vec<String> {
        use openms::concept::log_stream::{LogLevel, with_thread_local_log};
        with_thread_local_log(LogLevel::Warn, |log| log.clear_cache()).unwrap();
        let bytes = std::mem::take(&mut *self.0.lock().unwrap());
        String::from_utf8(bytes)
            .unwrap()
            .lines()
            .map(str::to_owned)
            .collect()
    }
}

struct Case {
    name: &'static str,
    input: &'static str,
    allowed: &'static [FileType],
    rewrite: bool,
    hash: bool,
    mz: Option<(f64, f64)>,
}

const fn case(
    name: &'static str,
    input: &'static str,
    allowed: &'static [FileType],
    rewrite: bool,
    hash: bool,
) -> Case {
    Case {
        name,
        input,
        allowed,
        rewrite,
        hash,
        mz: None,
    }
}

/// Every case the port reproduces completely: the loaded experiment, the
/// progress calls and the warnings.
const CASES: &[Case] = &[
    case("mzml_plain", "MzMLFile_1.mzML", &[FileType::MzMl], false, false),
    case("mzml_rewrite", "MzMLFile_1.mzML", &[FileType::MzMl], true, false),
    case("mzml_rewrite_hash", "MzMLFile_1.mzML", &[FileType::MzMl], true, true),
    case("mzml_gz_rewrite_hash", "MzMLFile_1.mzML.gz", &[], true, true),
    case("dta2d_plain", "DTA2DFile_test_1.dta2d", &[], false, false),
    Case {
        mz: Some((300.0, 1000.0)),
        ..case("dta2d_mz_rewrite_hash", "DTA2DFile_test_1.dta2d", &[FileType::Dta2d], true, true)
    },
    case("dta2d_rewrite", "DTA2DFile_test_1.dta2d", &[], true, false),
    case("mzdata_plain", "MzDataFile_1.mzData", &[], false, false),
    Case {
        mz: Some((115.0, 1000.0)),
        ..case("mzdata_mz", "MzDataFile_1.mzData", &[], false, false)
    },
    Case {
        mz: Some((115.0, 1000.0)),
        ..case("mzxml_mz", "MzXMLFile_1.mzXML", &[], false, false)
    },
    case("mzxml_plain", "MzXMLFile_1.mzXML", &[], false, false),
    case("mzxml_rewrite_hash", "MzXMLFile_1.mzXML", &[], true, true),
    case("mzdata_rewrite_hash", "MzDataFile_1.mzData", &[], true, true),
    case("srm", "srm_spectra.mzML", &[], false, false),
    case("srm_mixed", "srm_spectra_mixed.mzML", &[], true, true),
    case("ms2_rewrite_hash", "MS2File_test_spectra.ms2", &[], true, true),
    case("dta_rewrite_hash", "dta_extractor_1_output.dta", &[], true, true),
];

fn load(inputs: &Inputs, case: &Case, events: &Events) -> openms::Result<MSExperiment> {
    let options = ExperimentLoadOptions {
        peaks: case
            .mz
            .map_or_else(PeakFileOptions::default, |(min, max)| mz_range(min, max)),
        rewrite_source_file: case.rewrite,
        compute_hash: case.hash,
        ..Default::default()
    };
    FileHandler::load_experiment_with_progress(
        inputs.path(case.input),
        case.allowed,
        &options,
        &recording_logger(events),
    )
}

#[test]
fn release_loads_agree_in_experiment_progress_calls_and_warnings() {
    let inputs = Inputs::new();
    let warnings = Warnings::capture();
    for case in CASES {
        let events = Events::default();
        let experiment =
            load(&inputs, case, &events).unwrap_or_else(|error| panic!("{}: {error}", case.name));
        assert_eq!(release_kind(case.name, "R"), [""], "{}", case.name);
        assert_eq!(release_kind(case.name, "D"), ["0"], "{}", case.name);
        assert_eq!(
            *events.lock().unwrap(),
            release_kind(case.name, "P"),
            "{} progress calls",
            case.name
        );
        assert_eq!(
            warnings.take(),
            release_kind(case.name, "W"),
            "{} warnings",
            case.name
        );
        assert_experiment(case.name, &inputs.dir(), &experiment, ALL_ROWS);
    }
}

#[test]
fn release_relative_names_resolve_lexically_against_the_working_directory() {
    // File::absolutePath is std::filesystem::absolute, which keeps a `.`
    // component: `./MzMLFile_1.mzML` records the directory `<cwd>/.`.
    let inputs = Inputs::new();
    let warnings = Warnings::capture();
    let previous = std::env::current_dir().unwrap();
    std::env::set_current_dir(inputs.dir()).unwrap();
    let cases = [
        case("mzml_relative", "MzMLFile_1.mzML", &[], true, true),
        case("mzml_dot_relative", "./MzMLFile_1.mzML", &[], true, false),
    ];
    let mut results = Vec::new();
    for case in &cases {
        let options = ExperimentLoadOptions {
            rewrite_source_file: case.rewrite,
            compute_hash: case.hash,
            ..Default::default()
        };
        let events = Events::default();
        let result = FileHandler::load_experiment_with_progress(
            case.input,
            &[],
            &options,
            &recording_logger(&events),
        );
        results.push((result, events, warnings.take()));
    }
    std::env::set_current_dir(previous).unwrap();
    // Relative to the directory the loads ran in, as the driver's are.
    let cwd = inputs.dir();
    for (case, (result, events, lines)) in cases.iter().zip(results) {
        let experiment = result.unwrap_or_else(|error| panic!("{}: {error}", case.name));
        assert_eq!(*events.lock().unwrap(), release_kind(case.name, "P"));
        assert_eq!(lines, release_kind(case.name, "W"), "{}", case.name);
        assert_experiment(case.name, &cwd, &experiment, ALL_ROWS);
    }
}

#[test]
fn release_compute_file_hash_agrees_including_directories_and_missing_files() {
    let inputs = Inputs::new();
    let rows = release("hash");
    assert_eq!(rows.len(), 7);
    for row in rows {
        let (name, expected) = (row[0], row.get(1).copied().unwrap_or(""));
        assert_eq!(
            FileHandler::compute_file_hash(inputs.path(name)),
            expected,
            "{name}"
        );
    }
}

#[test]
fn release_divergences_are_the_documented_ones() {
    let inputs = Inputs::new();
    let _warnings = Warnings::capture();

    // MGF: the source's MascotGenericFile reports a section sized by the
    // file; the native MGF reader FileHandler dispatches to reports nothing.
    // The loaded experiment and its rewritten source file agree.
    let events = Events::default();
    let mgf = load(
        &inputs,
        &case("mgf_rewrite_hash", "mgf_two_blocks.mgf", &[], true, true),
        &events,
    )
    .unwrap();
    assert_eq!(
        release_kind("mgf_rewrite_hash", "P"),
        ["S\t0\t131\tloading MGF\t0", "V\t76\t1", "V\t-1\t1", "E\t0\t0"]
    );
    assert!(events.lock().unwrap().is_empty());
    assert_experiment("mgf_rewrite_hash", &inputs.dir(), &mgf, &["size", "sf"]);

    // A missing .dta2d: the source's DTA2DFile starts its section before it
    // fails to open the file and leaves it open; the port opens the file
    // before it dispatches, so the load fails with no call.
    let events = Events::default();
    let missing = load(
        &inputs,
        &case("missing_dta2d", "missing.dta2d", &[], false, false),
        &events,
    );
    assert!(matches!(missing, Err(Error::Io(_))), "{missing:?}");
    assert_eq!(
        release_kind("missing_dta2d", "P"),
        ["S\t0\t0\tloading DTA2D file\t0"]
    );
    assert_eq!(release_kind("missing_dta2d", "D"), ["1"]);
    assert!(events.lock().unwrap().is_empty());

    // A type outside the allowed list: ParseError there, InvalidValue here,
    // with no call on either side.
    let events = Events::default();
    let refused = load(
        &inputs,
        &case("dta2d_bogus_allowed", "DTA2DFile_test_1.dta2d", &[FileType::Unknown], true, false),
        &events,
    );
    assert!(matches!(refused, Err(Error::InvalidValue(_))), "{refused:?}");
    assert!(release_kind("dta2d_bogus_allowed", "X")[0].starts_with("Parse Error\t"));
    assert!(release_kind("dta2d_bogus_allowed", "P").is_empty());
    assert!(events.lock().unwrap().is_empty());
}

#[test]
fn every_release_case_is_exercised() {
    let mut named: Vec<&str> = CASES.iter().map(|case| case.name).collect();
    named.extend([
        "mzml_hash_only",
        "mzml_relative",
        "mzml_dot_relative",
        "mgf_rewrite_hash",
        "missing_dta2d",
        "dta2d_bogus_allowed",
        "hash",
    ]);
    let mut in_fixture: Vec<&str> = RELEASE
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.split('\t').next())
        .collect();
    in_fixture.dedup();
    in_fixture.sort_unstable();
    in_fixture.dedup();
    named.sort_unstable();
    assert_eq!(in_fixture, named);
}
