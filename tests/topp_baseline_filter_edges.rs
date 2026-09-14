// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! `BaselineFilter` against the Release C++ tool on the benchmark's edge shapes.
//!
//! The expected intensities come from the OpenMS4 Release `BaselineFilter`
//! (core `bc9cc12`, `x86_64` Linux) run on the same mzML fixtures, recorded by
//! `../oracle/baseline-filter-edges/run.sh`; see
//! `tests/data/baseline_filter_edges_provenance.json` and
//! `docs/MORPHOLOGICAL_FILTER_SUPPORT.md`. The `#[ignore]` tests at the end
//! read the benchmark's own inputs from `/ceph` and only run on the IBMI nodes.

// The TOPP framework lives behind `paramxml` and these tools read mzML, so the
// whole file is inert without both features.
#![cfg(all(feature = "mzml", feature = "paramxml"))]

use openms::cli::tools::BaselineFilter;
use openms::cli::{ExitCode, run_with};
use openms::format::file_handler::FileHandler;
use openms::format::file_types::FileType;
use openms::format::mzml::{ReadOptions, read_with_options};
use openms::kernel::SpectrumType;
use openms::processing::SpectrumFilter;
use openms::processing::baseline::{MorphologicalFilter, MorphologicalMethod, StructuringElement};
use openms::system::file::TempDir;
use openms::{MSExperiment, MSSpectrum};
use std::io::BufReader;
use std::path::{Path, PathBuf};

const TOOL_EXPECTED: &str = include_str!("data/baseline_filter_edges_tool_expected.tsv");
/// The benchmark's own inputs and the Release tool's output over them.
const BENCH_INPUTS: &str = "/ceph/ibmi/abi/oliver/bench/openms4/inputs";
const BENCH_CPP: &str = "/ceph/ibmi/abi/oliver/oracle-runs/baseline-filter-edges/results/scale";

fn fixture(name: &str) -> PathBuf {
    Path::new("tests/data").join(name)
}

fn run(args: &[&str]) -> (ExitCode, String) {
    let arguments: Vec<String> = std::iter::once("BaselineFilter".to_string())
        .chain(args.iter().map(|a| (*a).to_string()))
        .collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<BaselineFilter>(&arguments, &mut out, &mut err);
    (code, String::from_utf8_lossy(&err).into_owned())
}

/// The recorded C++ spectra of one tool run: `(type, m/z, intensities)`.
fn expected(run: &str) -> Vec<(String, Vec<f64>, Vec<f32>)> {
    TOOL_EXPECTED
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|f| f[0] == run)
        .map(|f| {
            let parse = |text: &str| -> Vec<String> {
                if text == "-" {
                    Vec::new()
                } else {
                    text.split(',').map(str::to_string).collect()
                }
            };
            (
                f[2].to_string(),
                parse(f[3]).iter().map(|v| v.parse().unwrap()).collect(),
                parse(f[4]).iter().map(|v| v.parse().unwrap()).collect(),
            )
        })
        .collect()
}

/// Intensities bitwise, m/z to the writer's precision.
fn assert_spectra(produced: &[MSSpectrum], expected: &[(String, Vec<f64>, Vec<f32>)], at: &str) {
    assert_eq!(produced.len(), expected.len(), "{at}: spectrum count");
    for (index, (got, (kind, mz, intensity))) in produced.iter().zip(expected).enumerate() {
        assert_eq!(got.peaks.len(), mz.len(), "{at}: spectrum {index} peaks");
        assert_eq!(kind, "profile", "{at}: spectrum {index}: C++ output type");
        assert_eq!(
            got.spectrum_type,
            SpectrumType::Profile,
            "{at}: spectrum {index} type"
        );
        for (peak, (m, i)) in got.peaks.iter().zip(mz.iter().zip(intensity)) {
            assert!((peak.mz - m).abs() <= 1e-9, "{at}: spectrum {index} m/z");
            assert_eq!(
                peak.intensity.to_bits(),
                i.to_bits(),
                "{at}: spectrum {index}: {} != {i}",
                peak.intensity
            );
        }
    }
}

fn tool_case(name: &str, input: &str, args: &[&str]) {
    let temp = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let out = temp.path().join("out.mzML");
    let mut arguments = vec![
        "-in".to_string(),
        fixture(input).to_string_lossy().into_owned(),
        "-out".to_string(),
        out.to_string_lossy().into_owned(),
    ];
    arguments.extend(args.iter().map(|a| (*a).to_string()));
    let borrowed: Vec<&str> = arguments.iter().map(String::as_str).collect();
    let (code, err) = run(&borrowed);
    assert_eq!(code, ExitCode::ExecutionOk, "{name}: {err}");
    let produced = FileHandler::load_experiment(&out, &[FileType::MzMl]).unwrap();
    let expected = expected(name);
    assert!(!expected.is_empty(), "{name}: no recorded C++ spectra");
    assert_spectra(&produced.spectra, &expected, name);
}

/// The benchmark's own invocation on its reproducer: spectrum index 3 of
/// `sub_profile_uk222_first600.mzML`, a centroided MS2 spectrum, with the
/// benchmark INI (top-hat, 1 Thomson). The C++ tool keeps the last peak.
#[test]
fn reproducer_with_the_benchmark_ini_matches_the_release_tool() {
    tool_case(
        "reproducer_ini",
        "baseline_filter_edges_reproducer.mzML",
        &[
            "-ini",
            &fixture("baseline_filter_edges_benchmark.ini").to_string_lossy(),
            "-threads",
            "1",
        ],
    );
}

/// Every edge shape through the tool: an empty spectrum, one and two peaks,
/// elements wider than, equal to and narrower than the spectrum, ties, and
/// profile and centroid spectra.
#[test]
fn edge_shapes_match_the_release_tool() {
    tool_case(
        "edges_tophat_th1",
        "baseline_filter_edges_edges.mzML",
        &[
            "-ini",
            &fixture("baseline_filter_edges_benchmark.ini").to_string_lossy(),
            "-threads",
            "1",
        ],
    );
    tool_case(
        "edges_tophat_default",
        "baseline_filter_edges_edges.mzML",
        &["-threads", "1"],
    );
    tool_case(
        "edges_erosion_dp1",
        "baseline_filter_edges_edges.mzML",
        &[
            "-struc_elem_unit",
            "DataPoints",
            "-struc_elem_length",
            "1",
            "-method",
            "erosion",
            "-threads",
            "1",
        ],
    );
}

/// The gradient's last sample under a one-sample element depends on what the
/// source's shared buffer holds from the previous spectra of the same run.
#[test]
fn history_dependent_gradient_matches_the_release_tool() {
    tool_case(
        "history_gradient_th003",
        "baseline_filter_edges_history.mzML",
        &[
            "-struc_elem_length",
            "0.03",
            "-method",
            "gradient",
            "-threads",
            "1",
        ],
    );
    tool_case(
        "history_tophat_th003",
        "baseline_filter_edges_history.mzML",
        &[
            "-struc_elem_length",
            "0.03",
            "-method",
            "tophat",
            "-threads",
            "1",
        ],
    );
}

// ---- scale, on the IBMI nodes only -----------------------------------------

/// Read a benchmark-sized mzML: the ceilings are raised for these files, whose
/// sizes are known from the benchmark manifest.
fn read_large(path: &Path) -> MSExperiment {
    let file = std::fs::File::open(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let options = ReadOptions {
        max_xml_bytes: 8 * 1024 * 1024 * 1024,
        max_array_bytes: 1024 * 1024 * 1024,
        max_total_peaks: 1_000_000_000,
        max_total_array_bytes: 8 * 1024 * 1024 * 1024,
        max_total_array_elements: 1_000_000_000,
        max_records: 10_000_000,
        max_total_arrays: 10_000_000,
        max_param_groups: 1_000_000,
        max_total_params: 1_000_000_000,
        max_param_bytes: 8 * 1024 * 1024 * 1024,
        ..Default::default()
    };
    read_with_options(BufReader::with_capacity(1 << 20, file), &options).unwrap()
}

/// Filter a real run and compare every intensity with the Release tool's output.
fn assert_scale_parity(input: &str, cpp_output: &str) {
    let mut produced = read_large(Path::new(input));
    MorphologicalFilter::new(
        MorphologicalMethod::TopHat,
        StructuringElement::Thomson(1.0),
    )
    .unwrap()
    .filter_experiment(&mut produced)
    .unwrap();
    let expected = read_large(Path::new(cpp_output));
    assert_eq!(produced.spectra.len(), expected.spectra.len(), "spectra");
    let mut compared = 0usize;
    for (index, (got, want)) in produced.spectra.iter().zip(&expected.spectra).enumerate() {
        assert_eq!(got.peaks.len(), want.peaks.len(), "spectrum {index} peaks");
        for (i, (a, b)) in got.peaks.iter().zip(&want.peaks).enumerate() {
            assert_eq!(
                a.intensity.to_bits(),
                b.intensity.to_bits(),
                "spectrum {index} peak {i}: {} != {}",
                a.intensity,
                b.intensity
            );
            compared += 1;
        }
    }
    assert!(compared > 100_000, "expected a full run, got {compared}");
}

/// The benchmark's 600-spectrum slice of UK222 (30 MB, 402 MS1 profile and 198
/// centroided MS2 spectra), against the Release tool's output. Needs the IBMI
/// `/ceph` share; run with `cargo test --test topp_baseline_filter_edges --
/// --ignored`.
#[test]
#[ignore = "reads /ceph/ibmi/abi/oliver on the IBMI nodes"]
fn uk222_first600_matches_the_release_tool() {
    assert_scale_parity(
        &format!("{BENCH_INPUTS}/derived/sub_profile_uk222_first600.mzML"),
        &format!("{BENCH_CPP}/uk222_first600_cpp.mzML"),
    );
}

/// The whole 2.3 GB UK222 run, 40,856 spectra, against the Release tool's
/// output. Needs the IBMI `/ceph` share and about 25 GB of memory.
#[test]
#[ignore = "reads /ceph/ibmi/abi/oliver on the IBMI nodes; multi-GB"]
fn uk222_full_matches_the_release_tool() {
    assert_scale_parity(
        &format!("{BENCH_INPUTS}/profile_hr_qe_silac_uk222/UK222.mzML"),
        &format!("{BENCH_CPP}/uk222_full_cpp.mzML"),
    );
}
