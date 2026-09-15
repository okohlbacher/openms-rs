// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! `-threads` reaches the run phase of every ported TOPP tool.
//!
//! The source applies the setting in `TOPPBase::main` before `main_`
//! (`TOPPBase.cpp:84-98, 408-415` at cli c19e494): a positive count sizes the
//! OpenMP team, and zero or a negative count means `omp_get_num_procs()`. The
//! port runs each tool body through `ToolContext::in_thread_pool`, a scoped
//! rayon pool of `ToolContext::thread_policy` workers.
//!
//! Evidence:
//!
//! * **Executed C++ (tier 1 probe).** The Release build
//!   `openms4-release-bc9cc12-c19e494-174b576` was run on ibminode06 with
//!   `/proc/<pid>/task` sampled every millisecond
//!   (`../oracle/tool-threads/threads_probe.py`, cases
//!   `make_cpp_cases.py`, results `cpp_results_ibminode06.jsonl`, sha256
//!   46f7d3ae8dd3b555bb0a016bc0925d33f7dd3aa0aa8d67522b6af56bbdcb2698).
//!   MapNormalizer on a 5000-spectrum slice showed: `-threads 4` gives four
//!   busy threads under `OMP_NUM_THREADS=1`; `-threads 0`, `-1` and `-7` give
//!   128 busy threads on the 128-processor node; `-threads 0` under
//!   `taskset -c 0-3` gives four. The threads C++ starts without
//!   `OMP_NUM_THREADS` (129, even for `-write_ini`) were traced with gdb to
//!   OpenBLAS `blas_thread_init` at library load, not to OpenMS.
//! * **Native invariants (tier 4).** The worker count observed inside the pool
//!   (`rayon::current_num_threads`, and on Linux the `openms-<index>` tasks of
//!   the real executables under `/proc`), and byte-identical outputs at 1, 2
//!   and more workers: the determinism contract of `src/concept/parallel.rs`.
//!
//! Five of the six tools are serial, so the pool does no parallel work for
//! them; what these tests prove is that the policy reaches their run phase.
//! `PeakPickerHiRes` is the sixth and the first whose run phase actually shares
//! work — its spectrum loop runs on the pool — and it is also the first to open
//! the pool around that region rather than around its whole body, so at
//! `-threads 1` it builds no pool at all and picks on the calling thread. Its
//! pool therefore lives exactly as long as its picking, and it is sampled on an
//! input it really picks (`profile`, `picking_input`) rather than on the
//! centroid-like one the other five tools share, which it would only copy: the
//! worker count is then held to exactly `-threads n`, as for every other tool,
//! and to none at one worker (`expected_workers`). Picking against the C++
//! outputs is covered by `tests/topp_peak_picker_hi_res.rs` and
//! `tests/peak_picking_experiment.rs`, and the worker counts of an
//! instrument-scale pick are sampled on the 2.3 GB benchmark run in
//! `docs/TOPP_PEAK_PICKER_HI_RES_SUPPORT.md`.
//!
//! Inputs are synthetic and generated per case into a `TempDir`. The
//! `#[ignore]`d `hpc_*` tests read the benchmark inputs under
//! `/ceph/ibmi/abi/oliver/bench/openms4/inputs` and run only on an IBMI Linux
//! node: `cargo test --release --test topp_threads -- --ignored`.

// The TOPP framework lives behind `paramxml` and these tools read mzML.
#![cfg(all(feature = "mzml", feature = "paramxml"))]

use openms::Result;
use openms::cli::tools::{
    BaselineFilter, DTAExtractor, MapNormalizer, MzMLSplitter, PeakPickerHiRes,
    SpectraFilterWindowMower,
};
use openms::cli::{ExitCode, Tool, ToolContext, ToolSpec, run_with};
use openms::concept::parallel::Threads;
use openms::format::file_handler::FileHandler;
use openms::format::file_types::FileType;
use openms::kernel::{MSExperiment, MSSpectrum, Peak1D, Precursor, SpectrumType};
use openms::system::file::TempDir;
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn text(path: impl AsRef<Path>) -> String {
    path.as_ref().to_string_lossy().into_owned()
}

fn temp() -> TempDir {
    TempDir::new(false).expect("temporary directory")
}

struct Outcome {
    code: ExitCode,
    err: String,
}

fn run<T: Tool>(args: &[String]) -> Outcome {
    let arguments: Vec<String> = std::iter::once(T::NAME.to_owned())
        .chain(args.iter().cloned())
        .collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<T>(&arguments, &mut out, &mut err);
    Outcome {
        code,
        err: String::from_utf8_lossy(&err).into_owned(),
    }
}

fn strings(args: &[&str]) -> Vec<String> {
    args.iter().map(|a| (*a).to_owned()).collect()
}

/// A centroid-like run: every fourth spectrum MS1, the others MS2 with a
/// precursor, `peaks` sorted peaks each, all intensities positive.
///
/// Peaks are 5 Th apart, so the default 50 Th WindowMower window holds ten of
/// them. Denser spectra of this size exceed that filter's cumulative
/// `max_work` of 50 million units across the run, a separate limit this lane
/// does not own.
fn synthetic(spectra: usize, peaks: usize) -> MSExperiment {
    let mut experiment = MSExperiment::new();
    for index in 0..spectra {
        let ms_level = if index % 4 == 0 { 1 } else { 2 };
        let precursors = if ms_level == 2 {
            vec![Precursor {
                mz: 400.0 + index as f64 * 0.25,
                charge: 2,
                ..Precursor::default()
            }]
        } else {
            Vec::new()
        };
        let peaks = (0..peaks)
            .map(|p| Peak1D {
                mz: 100.0 + p as f64 * 5.0 + index as f64 * 1e-3,
                intensity: ((p * 7919 + index * 104_729) % 10_007) as f32 + 1.0,
            })
            .collect();
        experiment.spectra.push(MSSpectrum {
            peaks,
            rt: 10.0 + index as f64 * 0.5,
            ms_level,
            native_id: format!("scan={}", index + 1),
            precursors,
            ..MSSpectrum::default()
        });
    }
    experiment
}

/// A profile run: `spectra` MS1 spectra, each carrying `peaks` Gaussian peaks
/// sampled nine times at 0.004 Th spacing, 1 Th apart from 400 Th.
///
/// [`synthetic`] is centroid-like, and `PeakPickerHiRes` copies such a record
/// rather than picking it (`selects`: a spectrum whose type is `Centroid` is
/// not selected in automatic mode), so a run over it does almost no work. These
/// records carry a profile peak shape *and* declare [`SpectrumType::Profile`],
/// which the mzML writer stores and the loader reads back, so the picker takes
/// the picking path for every one of them without depending on the type
/// estimator's heuristic. Nine samples per peak with a standard deviation of
/// 0.008 Th put four samples on each flank and none at zero intensity, and the
/// 1 Th spacing between peaks is far beyond the default `spacing_difference`,
/// so each peak is picked as one isolated centroid.
fn profile(spectra: usize, peaks: usize) -> MSExperiment {
    let mut experiment = MSExperiment::new();
    for index in 0..spectra {
        let mut samples = Vec::with_capacity(peaks * 9);
        for peak in 0..peaks {
            let centre = 400.0 + peak as f64;
            let height = 1_000.0 + ((peak * 7919 + index * 104_729) % 10_007) as f64;
            for step in 0..9_i32 {
                let offset = f64::from(step - 4) * 0.004;
                samples.push(Peak1D {
                    mz: centre + offset,
                    intensity: (height * (-(offset * offset) / (2.0 * 0.008 * 0.008)).exp()) as f32,
                });
            }
        }
        experiment.spectra.push(MSSpectrum {
            peaks: samples,
            rt: 10.0 + index as f64 * 0.5,
            ms_level: 1,
            native_id: format!("scan={}", index + 1),
            spectrum_type: SpectrumType::Profile,
            ..MSSpectrum::default()
        });
    }
    experiment
}

/// Store `experiment` as `input.mzML` in `dir` and return its path.
fn store_input(dir: &Path, experiment: &MSExperiment) -> String {
    let path = dir.join("input.mzML");
    FileHandler::store_experiment(&path, experiment, Some(FileType::MzMl))
        .expect("synthetic input is writable");
    text(path)
}

/// Every regular file in `dir`, by name.
fn files(dir: &Path) -> BTreeMap<String, Vec<u8>> {
    fs::read_dir(dir)
        .expect("output directory")
        .map(|entry| entry.expect("directory entry").path())
        .filter(|path| path.is_file())
        .map(|path| {
            let name = path.file_name().unwrap().to_string_lossy().into_owned();
            (name, fs::read(&path).expect("output file"))
        })
        .collect()
}

/// Each tool with the arguments that write its outputs into `out_dir`.
fn tool_arguments(tool: &str, input: &str, out_dir: &Path) -> Vec<String> {
    let out = |name: &str| text(out_dir.join(name));
    match tool {
        "BaselineFilter" | "MapNormalizer" | "PeakPickerHiRes" | "SpectraFilterWindowMower" => {
            vec!["-in".into(), input.into(), "-out".into(), out("out.mzML")]
        }
        "DTAExtractor" => vec!["-in".into(), input.into(), "-out".into(), out("spectrum")],
        "MzMLSplitter" => vec![
            "-in".into(),
            input.into(),
            "-out".into(),
            out("split"),
            "-parts".into(),
            "3".into(),
        ],
        other => panic!("unknown tool {other}"),
    }
}

const TOOLS: [&str; 6] = [
    "BaselineFilter",
    "DTAExtractor",
    "MapNormalizer",
    "MzMLSplitter",
    "PeakPickerHiRes",
    "SpectraFilterWindowMower",
];

fn run_tool(tool: &str, args: &[String]) -> Outcome {
    match tool {
        "BaselineFilter" => run::<BaselineFilter>(args),
        "DTAExtractor" => run::<DTAExtractor>(args),
        "MapNormalizer" => run::<MapNormalizer>(args),
        "MzMLSplitter" => run::<MzMLSplitter>(args),
        "PeakPickerHiRes" => run::<PeakPickerHiRes>(args),
        "SpectraFilterWindowMower" => run::<SpectraFilterWindowMower>(args),
        other => panic!("unknown tool {other}"),
    }
}

// ---------------------------------------------------------------------------
// A probe tool: what the body sees inside the pool
// ---------------------------------------------------------------------------

/// What [`ThreadProbe`]'s body observed.
#[derive(Clone, Debug)]
struct Seen {
    policy: Threads,
    threads: i64,
    /// `rayon::current_num_threads()` inside the pool; one without `parallel`.
    workers: usize,
    /// Whether the body ran on a rayon worker thread.
    on_worker: bool,
    /// Name of the thread the body ran on.
    thread_name: Option<String>,
    /// Whether the body ran on the thread that called `run_with`.
    on_calling_thread: bool,
}

thread_local! {
    static SEEN: RefCell<Option<Seen>> = const { RefCell::new(None) };
}

fn seen() -> Seen {
    SEEN.with(|slot| slot.borrow_mut().take())
        .expect("the probe body ran")
}

/// A tool whose body only records where it runs.
struct ThreadProbe;
impl Tool for ThreadProbe {
    const NAME: &'static str = "ThreadProbe";
    const DESCRIPTION: &'static str = "Records the thread policy the body runs under.";
    fn register(_spec: &mut ToolSpec) -> Result<()> {
        Ok(())
    }
    fn run(ctx: &ToolContext) -> Result<ExitCode> {
        let caller = std::thread::current().id();
        let (workers, on_worker, thread_name, body_thread) = ctx.in_thread_pool(|| {
            #[cfg(feature = "parallel")]
            let (workers, on_worker) = (
                rayon::current_num_threads(),
                rayon::current_thread_index().is_some(),
            );
            #[cfg(not(feature = "parallel"))]
            let (workers, on_worker) = (1, false);
            let current = std::thread::current();
            (
                workers,
                on_worker,
                current.name().map(str::to_owned),
                current.id(),
            )
        })?;
        let record = Seen {
            policy: ctx.thread_policy(),
            threads: ctx.threads(),
            workers,
            on_worker,
            thread_name,
            on_calling_thread: body_thread == caller,
        };
        SEEN.with(|slot| *slot.borrow_mut() = Some(record));
        Ok(ExitCode::ExecutionOk)
    }
}

fn probe(args: &[&str]) -> Seen {
    let outcome = run::<ThreadProbe>(&strings(args));
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    seen()
}

fn available() -> usize {
    Threads::all().get()
}

// ---------------------------------------------------------------------------
// The policy: TOPPBase semantics
// ---------------------------------------------------------------------------

/// A positive count is taken as given, from the command line or the INI file;
/// zero and every negative count mean every available processor, as
/// `if (num_threads <= 0) num_threads = omp_get_num_procs()`
/// (`TOPPBase.cpp:92-95`). Executed C++: `-threads -1` and `-threads -7` start
/// the same 128-thread team as `-threads 0` (cpp_results_ibminode06.jsonl,
/// cases `cpp t0 omp1`, `cpp t-1 omp1`, `cpp t-7 omp1`).
#[test]
fn policy_follows_toppbase_semantics() {
    // Source TOPPBase aborts a run without any option, so pass a neutral one.
    let default = probe(&["-no_progress"]);
    assert_eq!(default.threads, 1, "registered default");
    assert_eq!(default.policy, Threads::serial());

    for count in [1_i64, 2, 3, 7, 64, 1024] {
        let seen = probe(&["-threads", &count.to_string()]);
        assert_eq!(seen.threads, count);
        assert_eq!(seen.policy.get(), count as usize, "-threads {count}");
    }
    for count in ["0", "-1", "-7", "-2147483647"] {
        let seen = probe(&["-threads", count]);
        assert_eq!(seen.policy, Threads::all(), "-threads {count}");
    }
}

/// A request above both the ceiling and the available processors is clamped
/// to the larger of the two; below the ceiling nothing is clamped, even when
/// it oversubscribes the machine, as the source's team would.
#[test]
fn policy_clamps_only_absurd_counts() {
    let ceiling = ToolContext::THREAD_CEILING.max(available());
    let seen = probe(&["-threads", "2147483647"]);
    assert_eq!(seen.threads, 2_147_483_647);
    assert_eq!(seen.policy.get(), ceiling);
    let seen = probe(&["-threads", "1025"]);
    assert_eq!(seen.policy.get(), 1025.min(ceiling));
    const { assert!(ToolContext::THREAD_CEILING >= 1024) };
}

/// The INI `threads` value is the same parameter: `getParamAsInt_("threads")`
/// reads the merged tree, and the command line wins over the INI file.
#[test]
fn ini_threads_value_sets_the_policy() {
    let dir = temp();
    let ini = dir.path().join("probe.ini");
    let outcome = run::<ThreadProbe>(&strings(&["-write_ini", &text(&ini)]));
    assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
    let written = fs::read_to_string(&ini).unwrap();
    let from = r#"name="threads" value="1" type="int""#;
    assert!(written.contains(from), "{written}");
    fs::write(
        &ini,
        written.replace(from, r#"name="threads" value="3" type="int""#),
    )
    .unwrap();

    let seen = probe(&["-ini", &text(&ini)]);
    assert_eq!(seen.threads, 3);
    assert_eq!(seen.policy.get(), 3);
    let seen = probe(&["-ini", &text(&ini), "-threads", "2"]);
    assert_eq!(seen.policy.get(), 2);
}

// ---------------------------------------------------------------------------
// The pool
// ---------------------------------------------------------------------------

/// The body runs on a pool worker named `openms-<index>`, and rayon inside it
/// reports exactly the policy's worker count, for one worker too.
#[cfg(feature = "parallel")]
#[test]
fn body_runs_inside_a_pool_of_the_policy_size() {
    for count in [1_usize, 2, 3, 8] {
        let seen = probe(&["-threads", &count.to_string()]);
        assert_eq!(seen.workers, count, "-threads {count}");
        assert!(seen.on_worker, "-threads {count}: not on a pool worker");
        assert!(!seen.on_calling_thread, "-threads {count}");
        let name = seen.thread_name.expect("pool workers are named");
        assert!(name.starts_with("openms-"), "{name}");
    }
    let seen = probe(&["-threads", "0"]);
    assert_eq!(seen.workers, available());
}

/// Without `parallel` the body runs on the calling thread.
#[cfg(not(feature = "parallel"))]
#[test]
fn serial_build_runs_the_body_on_the_calling_thread() {
    let seen = probe(&["-threads", "8"]);
    assert_eq!(seen.policy.get(), 8, "the policy is still reported");
    assert!(seen.on_calling_thread);
    assert!(!seen.on_worker);
    assert_eq!(seen.workers, 1);
    assert_eq!(
        seen.thread_name.as_deref(),
        std::thread::current().name(),
        "the body ran on the test thread"
    );
}

/// An error raised inside the body passes through the pool unchanged, so the
/// run-phase exit-code mapping still applies: MzMLSplitter refuses `-parts 1`
/// with `-size 0` inside `main_` (`MzMLSplitter.cpp`: "Higher value for
/// parameter 'parts' or 'size' required"), and a truncated input fails while
/// loading, with the same code and diagnostics at every worker count.
#[test]
fn body_errors_pass_through_the_pool() {
    let dir = temp();
    let input = store_input(dir.path(), &synthetic(4, 3));
    let corrupt = dir.path().join("corrupt.mzML");
    fs::write(&corrupt, "<mzML><run>").unwrap();
    let out = text(dir.path().join("out.mzML"));
    let mut corrupt_outcomes = Vec::new();
    for threads in ["1", "4"] {
        let outcome = run::<MzMLSplitter>(&strings(&[
            "-in", &input, "-out", &out, "-threads", threads,
        ]));
        assert_eq!(outcome.code, ExitCode::IllegalParameters, "{}", outcome.err);
        assert!(
            outcome
                .err
                .contains("Higher value for parameter 'parts' or 'size' required"),
            "{}",
            outcome.err
        );
        let outcome = run::<MapNormalizer>(&strings(&[
            "-in",
            &text(&corrupt),
            "-out",
            &out,
            "-threads",
            threads,
        ]));
        assert_ne!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
        corrupt_outcomes.push((outcome.code, outcome.err));
    }
    assert_eq!(corrupt_outcomes[0], corrupt_outcomes[1]);
}

// ---------------------------------------------------------------------------
// Determinism
// ---------------------------------------------------------------------------

/// Under `-test`, every tool writes byte-identical files at 1, 2, 8 and all
/// available workers.
///
/// `PeakPickerHiRes` runs on profile records here, so what is compared across
/// worker counts for it is the output of its parallel spectrum loop rather than
/// a copy of the input: the picker leaves a centroid-like record alone.
#[test]
fn outputs_are_byte_identical_across_thread_counts() {
    let input_dir = temp();
    let input = store_input(input_dir.path(), &synthetic(120, 200));
    let picking_dir = temp();
    let picking = store_input(picking_dir.path(), &profile(60, 200));
    for tool in TOOLS {
        let input = if tool == "PeakPickerHiRes" {
            &picking
        } else {
            &input
        };
        let mut baseline: Option<BTreeMap<String, Vec<u8>>> = None;
        for threads in ["1", "2", "8", "0"] {
            let out_dir = temp();
            let mut args = tool_arguments(tool, input, out_dir.path());
            args.extend(strings(&["-test", "-threads", threads]));
            let outcome = run_tool(tool, &args);
            assert_eq!(
                outcome.code,
                ExitCode::ExecutionOk,
                "{tool} -threads {threads}: {}",
                outcome.err
            );
            let produced = files(out_dir.path());
            assert!(!produced.is_empty(), "{tool} wrote nothing");
            match &baseline {
                None => baseline = Some(produced),
                Some(expected) => {
                    assert_eq!(
                        produced.keys().collect::<Vec<_>>(),
                        expected.keys().collect::<Vec<_>>(),
                        "{tool} -threads {threads}: file names"
                    );
                    for (name, bytes) in &produced {
                        assert!(
                            bytes == &expected[name],
                            "{tool} -threads {threads}: {name} differs from -threads 1"
                        );
                    }
                }
            }
        }
    }
}

/// Outside `-test` the data are still identical; the processing record lists
/// the `threads` parameter as given and the completion time, which is what the
/// C++ outputs also differ in across thread counts (benchmark smoke run
/// 2026-09-14: "Only the completion time and the threads value differ").
#[test]
fn outside_test_mode_only_the_provenance_records_the_thread_count() {
    let input_dir = temp();
    let input = store_input(input_dir.path(), &synthetic(40, 50));
    // One output path for both runs, because the record also lists `-out`.
    let out_dir = temp();
    let normalised = |threads: &str| {
        let mut args = tool_arguments("MapNormalizer", &input, out_dir.path());
        args.extend(strings(&["-threads", threads]));
        let outcome = run::<MapNormalizer>(&args);
        assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
        let mut experiment =
            FileHandler::load_experiment(out_dir.path().join("out.mzML"), &[FileType::MzMl])
                .unwrap();
        let mut recorded = 0;
        for spectrum in &mut experiment.spectra {
            for record in &mut spectrum.data_processing {
                let record = std::sync::Arc::make_mut(record);
                record.completion_time = None;
                if let Some(value) = record.metadata.remove("parameter: threads") {
                    let value = value
                        .as_i64()
                        .ok()
                        .or_else(|| value.as_str().ok().and_then(|t| t.parse().ok()));
                    assert_eq!(value, threads.parse::<i64>().ok(), "-threads {threads}");
                    recorded += 1;
                }
            }
        }
        assert_eq!(recorded, 40, "-threads {threads}: one record per spectrum");
        experiment
    };
    let one = normalised("1");
    let many = normalised("7");
    assert_eq!(one.spectra.len(), 40);
    assert_eq!(one.chromatograms, many.chromatograms, "chromatograms");
    assert_eq!(one.settings, many.settings, "experiment settings");
    for (index, (a, b)) in one.spectra.iter().zip(&many.spectra).enumerate() {
        assert_eq!(a, b, "spectrum {index}");
    }
    assert!(one == many, "outputs differ beyond the processing record");
}

// ---------------------------------------------------------------------------
// Linux: the real executables, observed through /proc
// ---------------------------------------------------------------------------

#[cfg(all(target_os = "linux", feature = "parallel"))]
mod linux {
    use super::*;
    use std::process::{Command, ExitStatus, Stdio};
    use std::time::Duration;

    fn executable(tool: &str) -> &'static str {
        match tool {
            "BaselineFilter" => env!("CARGO_BIN_EXE_BaselineFilter"),
            "DTAExtractor" => env!("CARGO_BIN_EXE_DTAExtractor"),
            "MapNormalizer" => env!("CARGO_BIN_EXE_MapNormalizer"),
            "MzMLSplitter" => env!("CARGO_BIN_EXE_MzMLSplitter"),
            "PeakPickerHiRes" => env!("CARGO_BIN_EXE_PeakPickerHiRes"),
            "SpectraFilterWindowMower" => env!("CARGO_BIN_EXE_SpectraFilterWindowMower"),
            other => panic!("unknown tool {other}"),
        }
    }

    /// The `openms-<i>` worker count a sampler sees for `-threads n`.
    ///
    /// A tool that wraps its whole body in `ToolContext::in_thread_pool` holds
    /// the pool for the whole run, so every sample of a live process shows
    /// exactly `n` workers; the five tools ported before `PeakPickerHiRes` all
    /// do.
    ///
    /// `PeakPickerHiRes` opens the pool around its picking only
    /// (`src/cli/tools/peak_picker_hi_res.rs`, `pick_experiment`), and **at one
    /// worker builds none**: a pool of one worker bounds nothing that the
    /// `Threads` value the region is given does not already bound, and glibc
    /// charges the extra thread a second malloc arena. Its pool therefore lives
    /// exactly as long as the picking call, which is why it is sampled on
    /// [`picking_input`] rather than on the centroid-like [`large_input`] the
    /// other five tools use: over a record the picker *copies*, the pool can be
    /// gone between two samples of the 200 µs poll, and an earlier version of
    /// this test saw three workers of four for that reason. Over a record it
    /// picks, the pool is up for the whole spline pass — hundreds of samples —
    /// and the count is exact again.
    fn expected_workers(tool: &str, threads: usize) -> usize {
        match tool {
            "PeakPickerHiRes" if threads <= 1 => 0,
            _ => threads,
        }
    }

    /// What `/proc/<pid>/task` showed while one tool process ran.
    pub(super) struct Observation {
        /// Most `openms-<index>` tasks seen at once: the pool workers.
        pub workers: usize,
        /// Most tasks seen at once, main thread included.
        pub tasks: usize,
        pub status: ExitStatus,
        pub stderr: String,
    }

    /// Run `tool` with `args` in `work` and sample its tasks until it exits.
    ///
    /// `OMP_NUM_THREADS` and `RAYON_NUM_THREADS` are removed from the child's
    /// environment unless `env` sets them.
    pub(super) fn observe(
        tool: &str,
        args: &[String],
        env: &[(&str, &str)],
        work: &Path,
    ) -> Observation {
        let stderr_path = work.join("stderr.txt");
        let stderr_file = fs::File::create(&stderr_path).unwrap();
        let mut command = Command::new(executable(tool));
        command
            .args(args)
            .current_dir(work)
            .stdin(Stdio::null())
            .stdout(Stdio::null())
            .stderr(stderr_file)
            .env_remove("OMP_NUM_THREADS")
            .env_remove("RAYON_NUM_THREADS");
        for (key, value) in env {
            command.env(key, value);
        }
        let mut child = command.spawn().expect("tool executable starts");
        let task_dir = format!("/proc/{}/task", child.id());
        let (mut workers, mut tasks) = (0, 0);
        let status = loop {
            if let Some(status) = child.try_wait().unwrap() {
                break status;
            }
            if let Ok(entries) = fs::read_dir(&task_dir) {
                let (mut now_tasks, mut now_workers) = (0, 0);
                for entry in entries.flatten() {
                    now_tasks += 1;
                    if fs::read_to_string(entry.path().join("comm"))
                        .is_ok_and(|comm| comm.starts_with("openms-"))
                    {
                        now_workers += 1;
                    }
                }
                workers = workers.max(now_workers);
                tasks = tasks.max(now_tasks);
            }
            std::thread::sleep(Duration::from_micros(200));
        };
        Observation {
            workers,
            tasks,
            status,
            stderr: fs::read_to_string(stderr_path).unwrap_or_default(),
        }
    }

    /// A few MB: long enough that the pool is sampled many times in a debug
    /// build, below every current reader and writer limit.
    fn large_input(dir: &Path) -> String {
        store_input(dir, &synthetic(300, 1000))
    }

    /// Spectra and profile peaks per spectrum of [`picking_input`].
    const PICKING_SPECTRA: usize = 300;
    const PICKING_PEAKS: usize = 500;

    /// The input `PeakPickerHiRes` is sampled on: 300 profile spectra of 500
    /// nine-sample peaks, 1,350,000 samples and about 22 MB.
    ///
    /// Sized so that the picking — and with it the pool, which this tool opens
    /// around that call and nothing else — lasts long enough to be sampled many
    /// times. Measured through the picker's own entry point in a debug build on
    /// the gate host `dax` at load 58-69: **654 ms** at one worker, **339 ms**
    /// at two and **184 ms** at four, against the 200 µs poll of [`observe`].
    /// Several hundred samples therefore fall inside the pool's life at every
    /// count this test uses, and the sampled maximum is the pool's full width.
    /// The centroid-like [`large_input`] is copied rather than picked and is
    /// over in single-digit milliseconds, which is what made the same assertion
    /// flaky before.
    fn picking_input(dir: &Path) -> String {
        store_input(dir, &profile(PICKING_SPECTRA, PICKING_PEAKS))
    }

    /// The picked output really is picked: one centroid per profile peak, so
    /// the run this test sampled did the spline work that keeps the pool up.
    ///
    /// Without this, a change to [`profile`] that made the picker *copy* its
    /// records instead — a stored `Centroid` type, or samples the estimator
    /// reads as centroids — would shorten the pool to microseconds and turn the
    /// worker-count assertion into a flake rather than a failure.
    fn assert_the_input_was_picked(out_dir: &Path) {
        let picked = FileHandler::load_experiment(out_dir.join("out.mzML"), &[FileType::MzMl])
            .expect("the picked output loads");
        assert_eq!(picked.spectra.len(), PICKING_SPECTRA, "picked spectra");
        for (index, spectrum) in picked.spectra.iter().enumerate() {
            assert_eq!(
                spectrum.peaks.len(),
                PICKING_PEAKS,
                "spectrum {index} holds {} peaks, not the {PICKING_PEAKS} centroids of a picked \
                 profile record: the input was copied, not picked",
                spectrum.peaks.len()
            );
        }
    }

    /// Each executable starts exactly [`expected_workers`] pool workers for
    /// `-threads n` and no other thread besides its main thread, and writes the
    /// same bytes at every `n`.
    #[test]
    fn every_executable_runs_its_body_on_the_requested_pool() {
        let input_dir = temp();
        let input = large_input(input_dir.path());
        let picking_dir = temp();
        let picking = picking_input(picking_dir.path());
        for tool in TOOLS {
            // The picker's pool is open only while it picks, so it is sampled
            // on a profile input rather than on the copied centroid-like one.
            let input = if tool == "PeakPickerHiRes" {
                &picking
            } else {
                &input
            };
            let mut baseline: Option<BTreeMap<String, Vec<u8>>> = None;
            for n in [1_usize, 2, 4] {
                let work = temp();
                let out_dir = work.path().join("out");
                fs::create_dir(&out_dir).unwrap();
                let mut args = tool_arguments(tool, input, &out_dir);
                args.extend(strings(&["-test", "-threads", &n.to_string()]));
                let seen = observe(tool, &args, &[], work.path());
                assert!(
                    seen.status.success(),
                    "{tool} -threads {n}: {}",
                    seen.stderr
                );
                let workers = expected_workers(tool, n);
                assert_eq!(
                    seen.workers, workers,
                    "{tool} -threads {n}: pool workers seen"
                );
                assert_eq!(
                    seen.tasks,
                    workers + 1,
                    "{tool} -threads {n}: tasks seen, workers and the main thread"
                );
                if tool == "PeakPickerHiRes" && n == 1 {
                    assert_the_input_was_picked(&out_dir);
                }
                let produced = files(&out_dir);
                match &baseline {
                    None => baseline = Some(produced),
                    Some(expected) => assert!(
                        &produced == expected,
                        "{tool} -threads {n}: outputs differ from -threads 1"
                    ),
                }
            }
        }
    }

    /// Zero and negative counts start one worker per available processor,
    /// which the child inherits from this process (affinity and quota).
    #[test]
    fn zero_and_negative_counts_start_every_available_processor() {
        let input_dir = temp();
        let input = store_input(input_dir.path(), &synthetic(120, 1000));
        let expected = available();
        for count in ["0", "-1", "-7"] {
            let work = temp();
            let mut args = tool_arguments("MapNormalizer", &input, work.path());
            args.extend(strings(&["-threads", count]));
            let seen = observe("MapNormalizer", &args, &[], work.path());
            assert!(seen.status.success(), "-threads {count}: {}", seen.stderr);
            assert_eq!(seen.workers, expected, "-threads {count}");
        }
    }

    /// `OMP_NUM_THREADS` and `RAYON_NUM_THREADS` do not size the pool. The
    /// source's `omp_set_num_threads` overrides `OMP_NUM_THREADS` for the tool
    /// body (executed: cases `cpp t4 omp1` and `cpp t4 omp16` both run the
    /// body on four threads).
    #[test]
    fn thread_environment_variables_do_not_size_the_pool() {
        let input_dir = temp();
        let input = large_input(input_dir.path());
        let env = [("OMP_NUM_THREADS", "7"), ("RAYON_NUM_THREADS", "5")];
        for (args_threads, expected) in [(Some("2"), 2), (None, 1)] {
            let work = temp();
            let mut args = tool_arguments("SpectraFilterWindowMower", &input, work.path());
            if let Some(count) = args_threads {
                args.extend(strings(&["-threads", count]));
            }
            let seen = observe("SpectraFilterWindowMower", &args, &env, work.path());
            assert!(seen.status.success(), "{}", seen.stderr);
            assert_eq!(seen.workers, expected, "-threads {args_threads:?}");
            assert_eq!(seen.tasks, expected + 1);
        }
    }

    /// An INI `threads` value sizes the pool of the real executable, and
    /// `-threads` on the command line wins over it.
    #[test]
    fn ini_threads_value_sizes_the_pool() {
        let input_dir = temp();
        let input = large_input(input_dir.path());
        let work = temp();
        let ini = work.path().join("DTAExtractor.ini");
        let outcome = run::<DTAExtractor>(&strings(&["-write_ini", &text(&ini)]));
        assert_eq!(outcome.code, ExitCode::ExecutionOk, "{}", outcome.err);
        let written = fs::read_to_string(&ini).unwrap();
        let from = r#"name="threads" value="1" type="int""#;
        assert!(written.contains(from));
        fs::write(
            &ini,
            written.replace(from, r#"name="threads" value="3" type="int""#),
        )
        .unwrap();

        for (extra, expected) in [(vec![], 3), (strings(&["-threads", "2"]), 2)] {
            let out_dir = temp();
            let mut args = tool_arguments("DTAExtractor", &input, out_dir.path());
            args.extend(strings(&["-ini", &text(&ini)]));
            args.extend(extra);
            let seen = observe("DTAExtractor", &args, &[], work.path());
            assert!(seen.status.success(), "{}", seen.stderr);
            assert_eq!(seen.workers, expected);
        }
    }
}

// ---------------------------------------------------------------------------
// HPC: the benchmark inputs
// ---------------------------------------------------------------------------

/// Benchmark inputs on the IBMI Ceph file system (`inputs/MANIFEST.json` and
/// `inputs/derived/MANIFEST_sublimit.json` there).
#[cfg(all(target_os = "linux", feature = "parallel"))]
const BENCH_INPUTS: &str = "/ceph/ibmi/abi/oliver/bench/openms4/inputs";

/// One run's exit status, diagnostics and output files.
#[cfg(all(target_os = "linux", feature = "parallel"))]
type Run = (Option<i32>, String, BTreeMap<String, Vec<u8>>);

/// Run `tool` on `input` at 1, 2 and 16 workers: the pool has the requested
/// size, and the exit status, the diagnostics and every output byte agree.
#[cfg(all(target_os = "linux", feature = "parallel"))]
fn assert_thread_invariant_on(tool: &str, input: &str, extra: &[&str]) {
    assert!(
        Path::new(input).is_file(),
        "{input} is missing: the hpc_* tests run only on an IBMI node"
    );
    let mut baseline: Option<Run> = None;
    for n in [1_usize, 2, 16] {
        let work = temp();
        let out_dir = work.path().join("out");
        fs::create_dir(&out_dir).unwrap();
        let mut args = tool_arguments(tool, input, &out_dir);
        args.extend(strings(extra));
        args.extend(strings(&[
            "-test",
            "-no_progress",
            "-threads",
            &n.to_string(),
        ]));
        let seen = linux::observe(tool, &args, &[], work.path());
        assert_eq!(seen.workers, n, "{tool} -threads {n} on {input}");
        let result = (seen.status.code(), seen.stderr, files(&out_dir));
        eprintln!(
            "{tool} -threads {n}: exit {:?}, {} output files, tasks {}",
            result.0,
            result.2.len(),
            seen.tasks
        );
        match &baseline {
            None => baseline = Some(result),
            Some(expected) => {
                assert_eq!(result.0, expected.0, "{tool} -threads {n}: exit status");
                assert_eq!(result.1, expected.1, "{tool} -threads {n}: diagnostics");
                assert!(
                    result.2 == expected.2,
                    "{tool} -threads {n}: outputs differ"
                );
            }
        }
    }
}

/// The first-600 and first-5000 spectrum slices of UK222 the benchmark cut.
#[cfg(all(target_os = "linux", feature = "parallel"))]
#[test]
#[ignore = "HPC: reads /ceph/ibmi/abi/oliver/bench/openms4/inputs/derived; run on an IBMI node"]
fn hpc_benchmark_slices_are_thread_invariant() {
    let derived = format!("{BENCH_INPUTS}/derived");
    let centroid600 = format!("{derived}/sub_centroid_uk222_picked_first600.mzML");
    let centroid5000 = format!("{derived}/sub_centroid_uk222_picked_first5000.mzML");
    let profile600 = format!("{derived}/sub_profile_uk222_first600.mzML");
    assert_thread_invariant_on("DTAExtractor", &centroid5000, &[]);
    assert_thread_invariant_on("MzMLSplitter", &centroid600, &[]);
    assert_thread_invariant_on("SpectraFilterWindowMower", &centroid600, &[]);
    assert_thread_invariant_on("MapNormalizer", &centroid600, &[]);
    assert_thread_invariant_on(
        "BaselineFilter",
        &profile600,
        &["-struc_elem_length", "1.0"],
    );
}

/// The full-size smoke inputs: PXD001819 50amol_R1 (1.2 GB) for the four
/// centroid tools and UK222 (2.3 GB) for BaselineFilter. Until the mzML
/// reader and writer lanes land, the tools exit early with the same status at
/// every thread count; the pool size and that agreement are asserted either
/// way.
#[cfg(all(target_os = "linux", feature = "parallel"))]
#[test]
#[ignore = "HPC: reads multi-GB inputs under /ceph/ibmi/abi/oliver/bench/openms4/inputs"]
fn hpc_full_size_inputs_are_thread_invariant() {
    let velos = format!("{BENCH_INPUTS}/centroid_lcms_velos_pxd001819_50amol_r1/50amol_R1.mzML");
    let uk222 = format!("{BENCH_INPUTS}/profile_hr_qe_silac_uk222/UK222.mzML");
    for tool in [
        "DTAExtractor",
        "MzMLSplitter",
        "SpectraFilterWindowMower",
        "MapNormalizer",
    ] {
        assert_thread_invariant_on(tool, &velos, &[]);
    }
    assert_thread_invariant_on("BaselineFilter", &uk222, &["-struc_elem_length", "1.0"]);
}
