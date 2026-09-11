// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The upstream `TOPP_MzMLSplitter_*` tests, compared canonically against the
//! retained C++ outputs.
//!
//! The upstream test uses FuzzyDiff, not a byte comparison, so byte equality is
//! not the contract for mzML here either. Per `docs/DIFFERENTIAL_VALIDATION.md`
//! this compares decoded content: how many records each part received, which
//! records they are, and their peak values.

use openms::cli::{ExitCode, Tool, ToolContext, ToolSpec, run_with};
use openms::format::file_handler::FileHandler;
use openms::format::file_types::{FileType, strip_extension};
use openms::kernel::{MSExperiment, MSSpectrum};
use openms::{Error, Result};
use std::fs;
use std::path::{Path, PathBuf};

// Duplicated from src/bin/MzMLSplitter.rs; a binary target cannot be imported.
struct MzMLSplitter;

impl Tool for MzMLSplitter {
    const NAME: &'static str = "MzMLSplitter";
    const DESCRIPTION: &'static str = "Splits an mzML file into multiple parts";

    fn register(spec: &mut ToolSpec) -> Result<()> {
        spec.register_input_file("in", "<file>", "", "Input file", true, false, &[])?;
        spec.set_valid_formats("in", &["mzML"])?;
        spec.register_output_prefix(
            "out",
            "<prefix>",
            "",
            "Prefix for output files",
            false,
            false,
        )?;
        spec.register_int_option("parts", "<num>", 1, "Number of parts", false, false)?;
        spec.set_min_int("parts", 1)?;
        spec.register_int_option("size", "<num>", 0, "Approximate size limit", false, false)?;
        spec.set_min_int("size", 0)?;
        spec.register_string_option("unit", "<choice>", "MB", "Unit for 'size'", false, false)?;
        spec.set_valid_strings("unit", &["KB", "MB", "GB"])?;
        spec.register_flag("no_chrom", "Remove chromatograms.", false)?;
        spec.register_flag("no_spec", "Remove spectra.", false)?;
        Ok(())
    }

    fn run(ctx: &ToolContext) -> Result<ExitCode> {
        let input = ctx.string("in")?.to_owned();
        let mut out = ctx.string("out")?.to_owned();
        if out.is_empty() {
            out = strip_extension(&input).to_owned();
        }
        let (no_chrom, no_spec) = (ctx.flag("no_chrom")?, ctx.flag("no_spec")?);
        if no_chrom && no_spec {
            return Err(Error::InvalidValue(
                "'no_chrom' and 'no_spec' cannot be used together".into(),
            ));
        }
        let mut parts = usize::try_from(ctx.int("parts")?)
            .map_err(|_| Error::InvalidValue("'parts' must be positive".into()))?;
        let size = ctx.int("size")?;
        if parts == 1 {
            if size == 0 {
                return Err(Error::InvalidValue(
                    "Higher value for parameter 'parts' or 'size' required".into(),
                ));
            }
            let bytes = fs::metadata(&input)?.len() as f32;
            let total = match ctx.string("unit")? {
                "KB" => bytes / 1024.0,
                "MB" => bytes / (1024.0 * 1024.0),
                _ => bytes / (1024.0 * 1024.0 * 1024.0),
            };
            parts = (total / size as f32).ceil() as usize;
        }
        let experiment = FileHandler::load_experiment(&input, &[FileType::MzMl])?;
        let mut template = experiment.clone();
        let spectra = if no_spec {
            template.spectra.clear();
            Vec::new()
        } else {
            std::mem::take(&mut template.spectra)
        };
        let chromatograms = if no_chrom {
            template.chromatograms.clear();
            Vec::new()
        } else {
            std::mem::take(&mut template.chromatograms)
        };
        let width = parts.to_string().len();
        let (mut spec_start, mut chrom_start) = (0usize, 0usize);
        for counter in 1..=parts {
            let name = format!("{out}_part{counter:0width$}of{parts}.mzML");
            let mut part = MSExperiment {
                spectra: Vec::new(),
                chromatograms: Vec::new(),
                ..template.clone()
            };
            let remaining = parts - counter + 1;
            let n_spec = (spectra.len() - spec_start).div_ceil(remaining);
            part.spectra
                .extend_from_slice(&spectra[spec_start..spec_start + n_spec]);
            spec_start += n_spec;
            let n_chrom = (chromatograms.len() - chrom_start).div_ceil(remaining);
            part.chromatograms
                .extend_from_slice(&chromatograms[chrom_start..chrom_start + n_chrom]);
            chrom_start += n_chrom;
            FileHandler::store_experiment(&name, &part, Some(FileType::MzMl))?;
        }
        Ok(ExitCode::ExecutionOk)
    }
}

fn fixture(name: &str) -> PathBuf {
    Path::new("tests/data").join(name)
}

fn split(case: &str, args: &[&str]) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("openms-split-{}-{case}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let mut arguments = vec![
        "MzMLSplitter".to_string(),
        "-test".to_string(),
        "-in".to_string(),
        fixture("mzml_splitter_input.mzML")
            .to_string_lossy()
            .into_owned(),
        "-out".to_string(),
        dir.join("out").to_string_lossy().into_owned(),
    ];
    arguments.extend(args.iter().map(|a| (*a).to_string()));
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run_with::<MzMLSplitter>(&arguments, &mut out, &mut err);
    assert_eq!(
        code,
        ExitCode::ExecutionOk,
        "tool failed: {}",
        String::from_utf8_lossy(&err)
    );
    dir
}

/// Records of two runs agree when they hold the same spectra in the same order,
/// with equal native identifiers, MS levels, retention times and peak values.
fn assert_same_records(actual: &MSExperiment, expected: &MSExperiment, what: &str) {
    assert_eq!(
        actual.spectra.len(),
        expected.spectra.len(),
        "{what}: spectrum count"
    );
    assert_eq!(
        actual.chromatograms.len(),
        expected.chromatograms.len(),
        "{what}: chromatogram count"
    );
    for (i, (a, e)) in actual
        .spectra
        .iter()
        .zip(expected.spectra.iter())
        .enumerate()
    {
        assert_eq!(a.native_id, e.native_id, "{what}: spectrum {i} identifier");
        assert_eq!(a.ms_level, e.ms_level, "{what}: spectrum {i} MS level");
        assert!(
            (a.rt - e.rt).abs() <= 1e-9,
            "{what}: spectrum {i} retention time {} vs {}",
            a.rt,
            e.rt
        );
        assert_peaks(a, e, &format!("{what}: spectrum {i}"));
    }
}

fn assert_peaks(a: &MSSpectrum, e: &MSSpectrum, what: &str) {
    assert_eq!(a.peaks.len(), e.peaks.len(), "{what}: peak count");
    for (i, (pa, pe)) in a.peaks.iter().zip(e.peaks.iter()).enumerate() {
        assert!(
            (pa.mz - pe.mz).abs() <= 1e-9,
            "{what}: peak {i} m/z {} vs {}",
            pa.mz,
            pe.mz
        );
        assert!(
            (f64::from(pa.intensity) - f64::from(pe.intensity)).abs() <= 1e-6,
            "{what}: peak {i} intensity {} vs {}",
            pa.intensity,
            pe.intensity
        );
    }
}

fn load(path: impl AsRef<Path>) -> MSExperiment {
    FileHandler::load_experiment(path, &[FileType::MzMl]).unwrap()
}

#[test]
fn topp_mzml_splitter_1_splits_into_a_requested_number_of_parts() {
    let dir = split("1", &["-parts", "2"]);
    for part in 1..=2 {
        let produced = load(dir.join(format!("out_part{part}of2.mzML")));
        let reference = load(fixture(&format!("mzml_splitter_output_part{part}.mzML")));
        assert_same_records(&produced, &reference, &format!("part {part}"));
    }
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn topp_mzml_splitter_2_derives_the_part_count_from_a_size_limit() {
    // 40 KB over a ~59 KB input rounds up to the same two parts as case 1.
    let dir = split("2", &["-size", "40", "-unit", "KB"]);
    for part in 1..=2 {
        let produced = load(dir.join(format!("out_part{part}of2.mzML")));
        let reference = load(fixture(&format!("mzml_splitter_output_part{part}.mzML")));
        assert_same_records(&produced, &reference, &format!("part {part}"));
    }
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn every_record_is_placed_exactly_once() {
    // The remainder is spread over the parts still to come, so no record is
    // dropped or duplicated for a part count that does not divide evenly.
    let whole = load(fixture("mzml_splitter_input.mzML"));
    // parts=1 with no size limit is refused by the source, so it starts at 2.
    for parts in [2usize, 3, 4, 7, 11] {
        let dir = split(&format!("n{parts}"), &["-parts", &parts.to_string()]);
        let width = parts.to_string().len();
        let mut seen = Vec::new();
        for part in 1..=parts {
            let path = dir.join(format!("out_part{part:0width$}of{parts}.mzML"));
            seen.extend(load(&path).spectra.into_iter().map(|s| s.native_id));
        }
        let expected: Vec<String> = whole.spectra.iter().map(|s| s.native_id.clone()).collect();
        assert_eq!(seen, expected, "parts={parts} lost or reordered records");
        let _ = fs::remove_dir_all(&dir);
    }
}

#[test]
fn conflicting_and_missing_options_are_rejected() {
    let run = |args: &[&str]| -> ExitCode {
        let arguments: Vec<String> = std::iter::once("MzMLSplitter".to_string())
            .chain(args.iter().map(|a| (*a).to_string()))
            .collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        run_with::<MzMLSplitter>(&arguments, &mut out, &mut err)
    };
    let input = fixture("mzml_splitter_input.mzML")
        .to_string_lossy()
        .into_owned();
    // Source refuses both record filters at once.
    assert_eq!(
        run(&["-in", &input, "-no_chrom", "-no_spec"]),
        ExitCode::IllegalParameters
    );
    // One part and no size limit leaves nothing to do.
    assert_eq!(run(&["-in", &input]), ExitCode::IllegalParameters);
    // The unit is restricted to the registered set.
    assert_eq!(
        run(&["-in", &input, "-size", "1", "-unit", "TB"]),
        ExitCode::IllegalParameters
    );
    // parts has a registered minimum.
    assert_eq!(
        run(&["-in", &input, "-parts", "0"]),
        ExitCode::IllegalParameters
    );
}
