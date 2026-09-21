// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Progress reporting of the algorithms the source derives from
//! `ProgressLogger`, replayed against executed OpenMS4 Release output.
//!
//! `tests/data/progress_consumers_release.tsv` was written by the oracle in
//! `../oracle/progress-consumers` (see `tests/data/progress_consumers_provenance.json`):
//! the C++ Release build ran each algorithm twice on the inputs built below,
//! once with a recording backend (`setLogger`), which captured every call that
//! reached it with its arguments and nesting depth, and once with the command
//! backend (`setLogType(CMD)`), whose stdout bytes were captured with only the
//! two timing texts masked. The driver replaced `time()` so every call is a new
//! second and the whole-second throttle never suppresses a set; the clocks here
//! do the same. No expected value in this file is derived from Rust output.

use openms::analysis::peptide_indexing::{
    DecoyRule, MissingDecoyAction, PeptideIndexing, UnmatchedAction,
};
use openms::chemistry::{AASequence, DigestionSpecificity, Protease};
use openms::concept::parallel::Threads;
use openms::concept::progress_logger::{
    CommandProgressLogger, ProgressBackend, ProgressClock, ProgressLogger, ProgressNesting,
    ProgressReporter, ProgressTime, progress_value,
};
use openms::format::fasta::FASTAEntry;
use openms::identification::{PeptideHit, PeptideIdentification, ProteinIdentification};
use openms::kernel::{
    ChromatogramPeak, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, SpectrumType,
};
use openms::processing::baseline::MorphologicalFilter;
use openms::processing::iterative::PeakPickerIterative;
use openms::processing::peak_picking::{CENTROIDED_INPUT_MESSAGE, PeakPickerHiRes};
use openms::processing::smoothing::{GaussFilter, GaussianWidth, SavitzkyGolayFilter};
use openms::processing::{LinearResamplerAlign, SpectrumFilter};
use openms::{Error, Result};
use std::collections::BTreeMap;
use std::io::{self, Write};
use std::sync::atomic::{AtomicI64, Ordering};
use std::sync::{Arc, Mutex};

const FIXTURE: &str = include_str!("data/progress_consumers_release.tsv");

/// The summary line the command backend prints at `endProgress`, at depth 0,
/// with the two timing texts masked as the fixture masks them.
const MASKED_DONE_LINE: &str = "\r-- done [took <TIME> (CPU), <TIME> (Wall)] -- \n";

// ---------------------------------------------------------------------------
// Inputs, built exactly as the driver builds them.

/// Two triangles on 41 samples; small integers, exact in `f32` and `f64`.
fn triangle(i: i32) -> f32 {
    let a = 1000 - 150 * (i - 10).abs();
    let b = 600 - 100 * (i - 30).abs();
    a.max(b).max(10) as f32
}
fn profile(rt: f64, ms_level: u32) -> MSSpectrum {
    MSSpectrum {
        rt,
        ms_level,
        spectrum_type: SpectrumType::Profile,
        peaks: (0..41)
            .map(|i| Peak1D::new(400.0 + 0.01 * f64::from(i), triangle(i)))
            .collect(),
        ..Default::default()
    }
}
fn centroid(rt: f64, ms_level: u32) -> MSSpectrum {
    MSSpectrum {
        rt,
        ms_level,
        spectrum_type: SpectrumType::Centroid,
        peaks: (0..5)
            .map(|i| Peak1D::new(300.0 + f64::from(i), (100 * (i + 1)) as f32))
            .collect(),
        ..Default::default()
    }
}
fn chromatogram() -> MSChromatogram {
    MSChromatogram {
        peaks: (0..41)
            .map(|i| ChromatogramPeak::new(10.0 + 0.5 * f64::from(i), triangle(i)))
            .collect(),
        ..Default::default()
    }
}
fn experiment(spectra: Vec<MSSpectrum>, chromatograms: usize) -> MSExperiment {
    MSExperiment {
        spectra,
        chromatograms: (0..chromatograms).map(|_| chromatogram()).collect(),
        ..Default::default()
    }
}
fn mixed(chromatograms: usize) -> MSExperiment {
    experiment(
        vec![profile(1.0, 1), centroid(2.0, 2), profile(3.0, 1)],
        chromatograms,
    )
}
fn sizes(experiment: &MSExperiment) -> Vec<String> {
    vec![
        experiment.spectra.len().to_string(),
        experiment.chromatograms.len().to_string(),
    ]
}

fn database() -> Vec<FASTAEntry> {
    [
        ("P1", "MKAAPEPTIDERGGK"),
        ("DECOY_P1", "KGGREDITPEPAAKM"),
        ("P2", "MRSSSSAAAK"),
    ]
    .into_iter()
    .map(|(identifier, sequence)| FASTAEntry {
        identifier: identifier.into(),
        description: String::new(),
        sequence: sequence.into(),
    })
    .collect()
}
fn identification(sequences: &[&str]) -> PeptideIdentification {
    PeptideIdentification {
        identifier: "run1".into(),
        hits: sequences
            .iter()
            .zip(1..)
            .map(|(sequence, rank)| {
                PeptideHit::new(10.0, rank, 2, AASequence::parse(sequence).unwrap()).unwrap()
            })
            .collect(),
        ..Default::default()
    }
}
/// The driver's `configureIndexer`: explicit decoy prefix, warn on missing
/// decoys and unmatched hits, Trypsin with full specificity.
fn indexer() -> PeptideIndexing {
    PeptideIndexing {
        decoy_rule: DecoyRule::Prefix("DECOY_".into()),
        missing_decoy_action: MissingDecoyAction::Warn,
        unmatched_action: UnmatchedAction::Warn,
        enzyme: Some(Protease::Trypsin),
        specificity: Some(DigestionSpecificity::Full),
        ..Default::default()
    }
}
/// The driver's `index`, with the source's exit code in its `R` row mapped to
/// the port's outcome: `EXECUTION_OK` (0) and `PEPTIDE_IDS_EMPTY` (2) are both
/// a successful report here, told apart by whether any hit was indexed, and
/// `DATABASE_EMPTY` (1) is a refusal.
fn index(
    database: &[FASTAEntry],
    mut peptides: Vec<PeptideIdentification>,
    logger: &mut ProgressLogger,
) -> Result<Vec<String>> {
    let mut runs = vec![ProteinIdentification {
        identifier: "run1".into(),
        search_engine: "Mascot".into(),
        ..Default::default()
    }];
    let report = indexer().run_with_progress(database, &mut runs, &mut peptides, logger)?;
    Ok(vec![
        if report.peptide_hits == 0 { "2" } else { "0" }.to_string(),
    ])
}

/// Runs one fixture case with progress going to `logger`, returning the
/// driver's `R` fields: the output's spectrum and chromatogram counts, or the
/// indexing exit code.
fn run_case(case: &str, logger: &mut ProgressLogger) -> Result<Vec<String>> {
    let serial = Threads::serial();
    match case {
        "hires_mixed" => {
            let picked = PeakPickerHiRes::default().pick_experiment_with_progress(
                &mixed(2),
                serial,
                logger,
            )?;
            Ok(sizes(&picked.experiment))
        }
        "hires_empty" => {
            let picked = PeakPickerHiRes::default().pick_experiment_with_progress(
                &MSExperiment::default(),
                serial,
                logger,
            )?;
            Ok(sizes(&picked.experiment))
        }
        "hires_refused" => {
            let picker = PeakPickerHiRes {
                ms_levels: vec![1, 2],
                check_spectrum_type: true,
                ..Default::default()
            };
            let picked = picker.pick_experiment_with_progress(&mixed(1), serial, logger)?;
            Ok(sizes(&picked.experiment))
        }
        "iterative" => {
            let input = experiment(vec![profile(1.0, 1), profile(2.0, 2), profile(3.0, 1)], 1);
            let picked =
                PeakPickerIterative::default().pick_experiment_with_progress(&input, logger)?;
            Ok(sizes(&picked.experiment))
        }
        "iterative_empty" => {
            let picked = PeakPickerIterative::default()
                .pick_experiment_with_progress(&MSExperiment::default(), logger)?;
            Ok(sizes(&picked.experiment))
        }
        "resampler" | "resampler_empty" => {
            let mut map = if case == "resampler" {
                experiment(vec![profile(1.0, 1), profile(2.0, 1)], 1)
            } else {
                MSExperiment::default()
            };
            LinearResamplerAlign::default().raster_experiment_with_progress(&mut map, logger)?;
            Ok(sizes(&map))
        }
        "gauss" | "gauss_empty" | "gauss_ppm_chromatogram" => {
            let (filter, mut map) = match case {
                "gauss" => (
                    GaussFilter::default(),
                    experiment(vec![profile(1.0, 1), profile(2.0, 1)], 2),
                ),
                "gauss_empty" => (GaussFilter::default(), MSExperiment::default()),
                // The source's `use_ppm_tolerance` with its default
                // `ppm_tolerance` of 10.
                _ => (
                    GaussFilter::new(GaussianWidth::Ppm(10.0))?,
                    experiment(vec![profile(1.0, 1), profile(2.0, 1)], 1),
                ),
            };
            filter.filter_experiment_with_progress(&mut map, logger)?;
            Ok(sizes(&map))
        }
        "sgolay" | "sgolay_empty" => {
            let mut map = if case == "sgolay" {
                experiment(vec![profile(1.0, 1), profile(2.0, 1)], 2)
            } else {
                MSExperiment::default()
            };
            SavitzkyGolayFilter::default().filter_experiment_with_progress(&mut map, logger)?;
            Ok(sizes(&map))
        }
        "morph" | "morph_empty" => {
            let mut map = if case == "morph" {
                experiment(vec![profile(1.0, 1), profile(2.0, 1), profile(3.0, 1)], 1)
            } else {
                MSExperiment::default()
            };
            MorphologicalFilter::default().filter_experiment_with_progress(&mut map, logger)?;
            Ok(sizes(&map))
        }
        "indexing" => index(
            &database(),
            vec![
                identification(&["AAPEPTIDER"]),
                identification(&["EDITPEPAAK", "SSSSAAAK"]),
            ],
            logger,
        ),
        "indexing_no_identifications" => index(&database(), Vec::new(), logger),
        "indexing_no_hits" => index(&database(), vec![identification(&[])], logger),
        "indexing_empty_database" => index(&[], vec![identification(&["AAPEPTIDER"])], logger),
        other => panic!("fixture case {other:?} has no Rust counterpart"),
    }
}

// ---------------------------------------------------------------------------
// Loggers.

/// A clock whose every read is a new whole second, as the driver's `time()`,
/// with constant timer samples (the timing texts are masked either way).
fn new_second_clock() -> ProgressClock {
    let second = Arc::new(AtomicI64::new(1_000_000_000));
    Arc::new(move || {
        Ok(ProgressTime {
            wall_second: second.fetch_add(1, Ordering::Relaxed),
            wall_seconds: 0.0,
            cpu_seconds: Some(0.0),
        })
    })
}
/// A logger with its own nesting context, so parallel tests cannot share one.
fn isolated_logger() -> (ProgressLogger, ProgressNesting) {
    let nesting = ProgressNesting::default();
    (
        ProgressLogger::with_clock_and_nesting(new_second_clock(), nesting.clone()),
        nesting,
    )
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

/// Masks the two timing texts of every summary line, as `extract.py` does.
fn mask_timing(text: &str) -> String {
    const HEAD: &str = "-- done [took ";
    let mut out = String::new();
    let mut rest = text;
    while let Some(start) = rest.find(HEAD) {
        let after = &rest[start + HEAD.len()..];
        let end = after.find(" (Wall)").expect("wall timing text");
        assert!(after[..end].contains(" (CPU), "), "CPU timing text");
        out.push_str(&rest[..start]);
        out.push_str("-- done [took <TIME> (CPU), <TIME> (Wall)");
        rest = &after[end + " (Wall)".len()..];
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
            "B" => run.output = Some(unescape(rest[0])),
            other => panic!("unknown fixture record {other:?}"),
        }
    }
    (order, runs)
}

/// The port's error for the source exception a case records.
fn expected_error(case: &str, name: &str, message: &str) -> String {
    assert_eq!(name, "IllegalArgument", "{case}");
    match case {
        "hires_refused" => {
            assert_eq!(message, CENTROIDED_INPUT_MESSAGE);
            format!("invalid value: {CENTROIDED_INPUT_MESSAGE}")
        }
        "gauss_ppm_chromatogram" => {
            assert_eq!(
                message,
                "GaussFilter: Cannot use ppm tolerance on chromatograms"
            );
            "invalid value: ppm Gaussian smoothing is not defined for chromatograms".into()
        }
        other => panic!("no error mapping for {other}"),
    }
}

/// The port's R fields for a case whose output differs from the source's by
/// a documented native difference, and the source's fields otherwise.
fn expected_fields(case: &str, source: &[String]) -> Vec<String> {
    if case == "iterative" {
        // Source `pickExperiment` does not copy chromatograms into its output
        // (`PeakPickerIterative.h:372-378`); the port keeps them
        // (docs/ITERATIVE_PICKING_SUPPORT.md).
        assert_eq!(source, ["3", "0"]);
        return vec!["3".into(), "1".into()];
    }
    source.to_vec()
}

/// What the port does where the Release build leaves a section open or makes
/// one the port cannot: for a case whose source run threw inside its section,
/// the port also ends the section (see `ProgressReporter::section`); for an
/// empty database the port refuses before any progress.
enum Divergence {
    None,
    EndsFailedSection,
    RefusesBeforeProgress,
}
fn divergence(case: &str, captured: &Captured) -> Divergence {
    match (case, captured.outcome.as_ref().unwrap()) {
        ("indexing_empty_database", _) => Divergence::RefusesBeforeProgress,
        (_, Err(_)) => Divergence::EndsFailedSection,
        _ => Divergence::None,
    }
}

// ---------------------------------------------------------------------------
// Tests.

#[test]
fn fixture_covers_every_case_in_both_modes() {
    let (order, runs) = fixture();
    assert_eq!(order.len(), 18, "{order:?}");
    for case in &order {
        let rec = &runs[&(case.clone(), "rec".into())];
        assert!(rec.outcome.is_some() && rec.depth.is_some(), "{case}");
        if case.starts_with("indexing") {
            // PeptideIndexing also writes "Merge took" to std::cout, so the
            // driver ran it with the recording backend only.
            assert!(!runs.contains_key(&(case.clone(), "cmd".into())), "{case}");
        } else {
            let cmd = &runs[&(case.clone(), "cmd".into())];
            assert!(cmd.output.is_some() && cmd.outcome == rec.outcome, "{case}");
            assert_eq!(cmd.depth, rec.depth, "{case}");
        }
    }
}

#[test]
fn every_backend_call_matches_the_release_build() {
    let (order, runs) = fixture();
    for case in &order {
        let captured = &runs[&(case.clone(), "rec".into())];
        let events = Arc::new(Mutex::new(Vec::new()));
        let (mut logger, nesting) = isolated_logger();
        logger.set_logger(Box::new(Recorder {
            events: events.clone(),
            current: 0,
        }));
        let outcome = run_case(case, &mut logger);
        let actual = events.lock().unwrap().clone();
        let mut expected_events = captured.events.clone();
        let mut expected_depth = captured.depth.unwrap();
        match divergence(case, captured) {
            Divergence::None => {}
            Divergence::EndsFailedSection => {
                // The Release build threw inside the section: no endProgress,
                // and the static depth stayed raised.
                assert_eq!(expected_depth, 1, "{case}");
                assert!(!expected_events.last().unwrap().starts_with('E'), "{case}");
                expected_events.push("E\t0\t0".into());
                expected_depth = 0;
            }
            Divergence::RefusesBeforeProgress => {
                assert_eq!(
                    expected_events,
                    ["S\t0\t1\tLoad first DB chunk\t0", "E\t0\t0"],
                    "{case}"
                );
                expected_events.clear();
            }
        }
        assert_eq!(actual, expected_events, "{case}: backend calls");
        assert_eq!(
            nesting.depth(),
            expected_depth,
            "{case}: nesting afterwards"
        );
        match (captured.outcome.as_ref().unwrap(), outcome) {
            (Ok(source), Ok(port)) => assert_eq!(port, expected_fields(case, source), "{case}"),
            (Ok(source), Err(error)) => {
                assert_eq!(case, "indexing_empty_database", "{case}: {error}");
                // DATABASE_EMPTY.
                assert_eq!(source, &["1"]);
                assert!(matches!(error, Error::InvalidValue(_)), "{error}");
            }
            (Err((name, message)), Err(error)) => {
                assert_eq!(error.to_string(), expected_error(case, name, message))
            }
            (Err(source), Ok(port)) => {
                panic!("{case}: source threw {source:?}, port gave {port:?}")
            }
        }
    }
}

#[test]
fn command_output_matches_the_release_build() {
    let (order, runs) = fixture();
    for case in order.iter().filter(|c| !c.starts_with("indexing")) {
        let captured = &runs[&(case.clone(), "cmd".into())];
        let output = SharedOutput::default();
        let (mut logger, nesting) = isolated_logger();
        logger.set_logger(Box::new(CommandProgressLogger::with_clock(
            output.clone(),
            new_second_clock(),
        )));
        let outcome = run_case(case, &mut logger);
        let mut expected = captured.output.clone().unwrap();
        let mut expected_depth = captured.depth.unwrap();
        if let Divergence::EndsFailedSection = divergence(case, captured) {
            // The Release output stops mid-line after the last percentage.
            assert!(!expected.ends_with('\n'), "{case}");
            expected.push_str(MASKED_DONE_LINE);
            expected_depth = 0;
            assert!(outcome.is_err(), "{case}");
        } else {
            assert!(outcome.is_ok(), "{case}: {outcome:?}");
        }
        assert_eq!(mask_timing(&output.text()), expected, "{case}: stdout");
        assert_eq!(nesting.depth(), expected_depth, "{case}");
    }
}

/// The in-place picker makes the same calls as the borrowing one, and both make
/// them at every worker count: spectra are counted as they are committed, in
/// input order.
#[test]
fn picking_reports_as_the_release_build_in_place_and_at_every_thread_count() {
    let (_, runs) = fixture();
    for case in ["hires_mixed", "hires_refused"] {
        let captured = &runs[&(case.to_string(), "rec".into())];
        let failed = captured.outcome.as_ref().unwrap().is_err();
        let mut expected = captured.events.clone();
        if failed {
            expected.push("E\t0\t0".into());
        }
        let (picker, input) = if case == "hires_mixed" {
            (PeakPickerHiRes::default(), mixed(2))
        } else {
            let picker = PeakPickerHiRes {
                ms_levels: vec![1, 2],
                check_spectrum_type: true,
                ..Default::default()
            };
            (picker, mixed(1))
        };
        for threads in [
            Threads::serial(),
            Threads::from_cli(2),
            Threads::from_cli(8),
        ] {
            for in_place in [false, true] {
                let events = Arc::new(Mutex::new(Vec::new()));
                let (mut logger, _) = isolated_logger();
                logger.set_logger(Box::new(Recorder {
                    events: events.clone(),
                    current: 0,
                }));
                let succeeded = if in_place {
                    let mut experiment = input.clone();
                    picker
                        .pick_experiment_in_place_with_progress(
                            &mut experiment,
                            threads,
                            &mut logger,
                        )
                        .is_ok()
                } else {
                    picker
                        .pick_experiment_with_progress(&input, threads, &mut logger)
                        .is_ok()
                };
                assert_eq!(succeeded, !failed, "{case}");
                assert_eq!(
                    *events.lock().unwrap(),
                    expected,
                    "{case}, {} threads, in place {in_place}",
                    threads.get()
                );
            }
        }
    }
}

/// Reporting progress changes no result: each `*_with_progress` entry point
/// returns exactly what its silent counterpart returns.
#[test]
fn progress_changes_no_result() {
    let command = || {
        let (mut logger, _) = isolated_logger();
        logger.set_logger(Box::new(CommandProgressLogger::with_clock(
            SharedOutput::default(),
            new_second_clock(),
        )));
        logger
    };
    let input = mixed(2);
    let picker = PeakPickerHiRes::default();
    let silent = picker.pick_experiment(&input).unwrap();
    for threads in [Threads::serial(), Threads::from_cli(4)] {
        let reported = picker
            .pick_experiment_with_progress(&input, threads, &mut command())
            .unwrap();
        assert_eq!(reported, silent);
        let mut in_place = input.clone();
        let report = picker
            .pick_experiment_in_place_with_progress(&mut in_place, threads, &mut command())
            .unwrap();
        assert_eq!(in_place, silent.experiment);
        assert_eq!(report.spectrum_boundaries, silent.spectrum_boundaries);
        assert_eq!(
            report.chromatogram_boundaries,
            silent.chromatogram_boundaries
        );
    }

    let input = experiment(vec![profile(1.0, 1), profile(2.0, 2), profile(3.0, 1)], 1);
    let iterative = PeakPickerIterative::default();
    let silent = iterative.pick_experiment(&input).unwrap();
    let reported = iterative
        .pick_experiment_with_progress(&input, &mut command())
        .unwrap();
    assert_eq!(reported.experiment, silent.experiment);
    assert_eq!(reported.spectrum_regions, silent.spectrum_regions);

    let base = experiment(vec![profile(1.0, 1), profile(2.0, 1)], 2);
    let check =
        |silent: &dyn Fn(&mut MSExperiment) -> Result<()>,
         reported: &dyn Fn(&mut MSExperiment, &mut ProgressLogger) -> Result<()>| {
            let mut a = base.clone();
            let mut b = base.clone();
            silent(&mut a).unwrap();
            reported(&mut b, &mut command()).unwrap();
            assert_eq!(a, b);
            assert_ne!(
                a, base,
                "the filter changed nothing, so the check is vacuous"
            );
        };
    let gauss = GaussFilter::default();
    check(&|m| gauss.filter_experiment(m), &|m, l| {
        gauss.filter_experiment_with_progress(m, l)
    });
    let sgolay = SavitzkyGolayFilter::default();
    check(&|m| sgolay.filter_experiment(m), &|m, l| {
        sgolay.filter_experiment_with_progress(m, l)
    });
    let morph = MorphologicalFilter::default();
    check(&|m| morph.filter_experiment(m), &|m, l| {
        morph.filter_experiment_with_progress(m, l)
    });
    let resampler = LinearResamplerAlign::default();
    check(&|m| resampler.raster_experiment(m), &|m, l| {
        resampler.raster_experiment_with_progress(m, l)
    });

    let peptides = vec![
        identification(&["AAPEPTIDER"]),
        identification(&["EDITPEPAAK", "SSSSAAAK"]),
    ];
    let runs = vec![ProteinIdentification {
        identifier: "run1".into(),
        search_engine: "Mascot".into(),
        ..Default::default()
    }];
    let (mut silent_runs, mut silent_peptides) = (runs.clone(), peptides.clone());
    let silent = indexer()
        .run(&database(), &mut silent_runs, &mut silent_peptides)
        .unwrap();
    let (mut reported_runs, mut reported_peptides) = (runs, peptides);
    let reported = indexer()
        .run_with_progress(
            &database(),
            &mut reported_runs,
            &mut reported_peptides,
            &mut command(),
        )
        .unwrap();
    assert_eq!(reported, silent);
    assert_eq!(reported_runs, silent_runs);
    assert_eq!(reported_peptides, silent_peptides);
    assert_eq!(silent.peptide_hits, 3);
}

/// `LinearResamplerAlign` now has an experiment entry point: source
/// `rasterExperiment` resamples every spectrum, leaves chromatograms alone,
/// and the trait's experiment call is the same operation.
#[test]
fn raster_experiment_resamples_spectra_only() {
    let resampler = LinearResamplerAlign::default();
    let input = experiment(vec![profile(1.0, 1), profile(2.0, 1)], 1);
    let mut rastered = input.clone();
    resampler.raster_experiment(&mut rastered).unwrap();
    for (before, after) in input.spectra.iter().zip(&rastered.spectra) {
        let mut expected = before.clone();
        resampler.raster(&mut expected).unwrap();
        assert_eq!(after, &expected);
    }
    assert_eq!(rastered.chromatograms, input.chromatograms);
    let mut through_trait = input.clone();
    resampler.filter_experiment(&mut through_trait).unwrap();
    assert_eq!(through_trait, rastered);
    let mut one = input.spectra[0].clone();
    resampler.filter_spectrum(&mut one).unwrap();
    assert_eq!(one, rastered.spectra[0]);
    // Atomic: a spectrum the port refuses leaves the experiment unchanged.
    let mut refused = input.clone();
    refused.spectra[1].peaks.swap(0, 5);
    let before = refused.clone();
    assert!(resampler.raster_experiment(&mut refused).is_err());
    assert_eq!(refused, before);
}

/// A database of exactly the source's first-chunk size scans with the range
/// the source uses for a database of unknown size (`PeptideIndexing.cpp:453`).
#[test]
fn indexing_scan_range_follows_the_source_chunk_rule() {
    use openms::analysis::peptide_indexing::SOURCE_PROTEIN_CACHE_SIZE;
    let entries = |n: usize| -> Vec<FASTAEntry> {
        (0..n)
            .map(|i| FASTAEntry {
                identifier: format!("P{i}"),
                description: String::new(),
                sequence: "K".into(),
            })
            .collect()
    };
    for (n, range) in [
        (SOURCE_PROTEIN_CACHE_SIZE, i64::MAX),
        (SOURCE_PROTEIN_CACHE_SIZE + 1, 400_001),
    ] {
        let events = Arc::new(Mutex::new(Vec::new()));
        let (mut logger, _) = isolated_logger();
        logger.set_logger(Box::new(Recorder {
            events: events.clone(),
            current: 0,
        }));
        let mut runs = vec![ProteinIdentification {
            identifier: "run1".into(),
            ..Default::default()
        }];
        let mut peptides = vec![identification(&["K"])];
        let config = PeptideIndexing {
            max_matches: 2_000_000,
            ..indexer()
        };
        config
            .run_with_progress(&entries(n), &mut runs, &mut peptides, &mut logger)
            .unwrap();
        let events = events.lock().unwrap();
        assert_eq!(events[2], format!("S\t0\t{range}\tAho-Corasick\t0"));
        assert_eq!(events.len(), 2 + 1 + n + 1);
        assert_eq!(events[3], "V\t1\t1");
        assert_eq!(events[2 + n], format!("V\t{n}\t1"));
    }
}

// ---------------------------------------------------------------------------
// ProgressReporter itself (independent: expectations follow from the
// documented contract, not from Rust output).

#[test]
fn a_silent_reporter_calls_nothing() {
    let mut reporter = ProgressReporter::silent();
    assert!(!reporter.is_reporting());
    let value = reporter
        .section(0, 10, "never shown", |reporter| {
            reporter.set(3)?;
            reporter.set_count(usize::MAX)?;
            Ok(7)
        })
        .unwrap();
    assert_eq!(value, 7);
}

#[test]
fn a_section_is_ended_when_its_body_fails() {
    let events = Arc::new(Mutex::new(Vec::new()));
    let (mut logger, nesting) = isolated_logger();
    logger.set_logger(Box::new(Recorder {
        events: events.clone(),
        current: 0,
    }));
    let mut reporter = ProgressReporter::new(Some(&mut logger));
    assert!(reporter.is_reporting());
    let failed: Result<()> = reporter.section(0, 4, "label", |reporter| {
        reporter.set(1)?;
        Err(Error::InvalidValue("body".into()))
    });
    assert_eq!(failed.unwrap_err().to_string(), "invalid value: body");
    assert_eq!(
        *events.lock().unwrap(),
        ["S\t0\t4\tlabel\t0", "V\t1\t1", "E\t0\t0"]
    );
    assert_eq!(nesting.depth(), 0);
}

/// A backend that refuses one operation.
struct Refusing {
    start: bool,
    end: bool,
    started: Arc<Mutex<u32>>,
}
impl ProgressBackend for Refusing {
    fn start_progress(&mut self, _: i64, _: i64, _: &str, _: usize) -> Result<()> {
        *self.started.lock().unwrap() += 1;
        if self.start {
            Err(Error::InvalidValue("start".into()))
        } else {
            Ok(())
        }
    }
    fn set_progress(&mut self, _: i64, _: usize) -> Result<()> {
        Ok(())
    }
    fn next_progress(&mut self) -> Result<i64> {
        Ok(0)
    }
    fn end_progress(&mut self, _: usize, _: u64) -> Result<()> {
        if self.end {
            Err(Error::InvalidValue("end".into()))
        } else {
            Ok(())
        }
    }
}

#[test]
fn a_failed_start_skips_the_body_and_a_body_error_outranks_the_end() {
    let refusing = |start: bool, end: bool| {
        let started = Arc::new(Mutex::new(0));
        let (mut logger, nesting) = isolated_logger();
        logger.set_logger(Box::new(Refusing {
            start,
            end,
            started: started.clone(),
        }));
        (logger, nesting, started)
    };
    let (mut logger, nesting, started) = refusing(true, false);
    let mut ran = false;
    let outcome = ProgressReporter::new(Some(&mut logger)).section(0, 1, "x", |_| {
        ran = true;
        Ok(())
    });
    assert_eq!(outcome.unwrap_err().to_string(), "invalid value: start");
    assert!(!ran);
    assert_eq!((*started.lock().unwrap(), nesting.depth()), (1, 0));

    let (mut logger, _, _) = refusing(false, true);
    let body: Result<()> = ProgressReporter::new(Some(&mut logger))
        .section(0, 1, "x", |_| Err(Error::InvalidValue("body".into())));
    assert_eq!(body.unwrap_err().to_string(), "invalid value: body");
    let (mut logger, _, _) = refusing(false, true);
    let end = ProgressReporter::new(Some(&mut logger)).section(0, 1, "x", |_| Ok(()));
    assert_eq!(end.unwrap_err().to_string(), "invalid value: end");
}

#[test]
fn progress_values_are_the_source_signed_size() {
    assert_eq!(progress_value(0).unwrap(), 0);
    assert_eq!(progress_value(400_000).unwrap(), 400_000);
    let max = usize::try_from(i64::MAX).unwrap();
    assert_eq!(progress_value(max).unwrap(), i64::MAX);
    if let Some(over) = max.checked_add(1) {
        assert!(matches!(progress_value(over), Err(Error::InvalidValue(_))));
    }
}
