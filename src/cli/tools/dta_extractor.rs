// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Extracts spectra of an mzML file to DTA files, as the source TOPP tool.
//!
//! Ports `OpenMS4-topp/src/DTAExtractor.cpp`. The tool lives here rather than in
//! `src/bin/` so that the shipped binary and its differential test share one
//! definition; a copied tool body drifts from the binary it claims to test.
//!
//! A range or MS level list that does not convert is the source's
//! `ConversionError`, which `main_` catches itself (`DTAExtractor.cpp:130-136`):
//! `Invalid boundary '<level>' given. Aborting!` with `writeLogError_`, the
//! usage text, and `ILLEGAL_PARAMETERS` as an exit code `main_` returns, so
//! `TOPPBase`'s closing line follows. `<level>` is the `-level` value once the
//! level list is being converted and empty before that, which is what the
//! source's variable holds (Release oracles `dta_bad_rt_log` and
//! `dta_bad_level_log` of `../oracle/topp-exception-exits`).

use crate::Result;
use crate::cli::context::to_int32;
use crate::cli::{ExitCode, PoolLines, Tool, ToolContext, ToolResult, ToolSpec, parse_range};
use crate::format::dta;
use crate::format::file_handler::FileHandler;
use crate::format::file_types::FileType;
use crate::kernel::MSSpectrum;
use std::fs::File;
use std::io::{BufWriter, Write};

/// The `DTAExtractor` TOPP tool.
pub struct DTAExtractor;

impl Tool for DTAExtractor {
    const NAME: &'static str = "DTAExtractor";
    const DESCRIPTION: &'static str =
        "Extracts spectra of an MS run file to several files in DTA format.";

    fn register(spec: &mut ToolSpec) -> Result<()> {
        spec.register_input_file("in", "<file>", "", "input file ", true, false, &[])?;
        spec.set_valid_formats("in", &["mzML"])?;
        spec.register_string_option(
            "out",
            "<file>",
            "",
            "base name of DTA output files (RT, m/z and extension are appended)",
            true,
            false,
        )?;
        spec.register_string_option(
            "mz",
            "[min]:[max]",
            ":",
            "m/z range of precursor peaks to extract.\nThis option is ignored for MS level 1",
            false,
            false,
        )?;
        spec.register_string_option(
            "rt",
            "[min]:[max]",
            ":",
            "retention time range of spectra to extract [s]",
            false,
            false,
        )?;
        spec.register_string_option(
            "level",
            "i[,j]...",
            "1,2,3",
            "MS levels to extract",
            false,
            false,
        )?;
        Ok(())
    }

    /// Run against the process streams; see `run_io`.
    fn run(ctx: &ToolContext) -> ToolResult {
        Self::run_io(ctx, &mut std::io::stdout(), &mut std::io::stderr())
    }

    /// Source `main_`, run on the worker pool that `-threads` sizes, as
    /// `TOPPBase::main` applies the setting before `main_`
    /// (`TOPPBase.cpp:408-415`). See [`ToolContext::in_thread_pool`]. The
    /// body's console lines and the usage text of a refusal are written once
    /// the pool returns ([`PoolLines`]).
    fn run_io(ctx: &ToolContext, out: &mut dyn Write, err: &mut dyn Write) -> ToolResult {
        let mut lines = PoolLines::default();
        let result = ctx.in_thread_pool(|| Self::run_in_pool(ctx, &mut lines))?;
        if lines.write(out, err)? {
            let spec = crate::cli::tool_spec::<Self>()?;
            let subsections = crate::cli::subsection_defaults::<Self>(&spec)?;
            crate::cli::print_usage::<Self>(err, &spec, &subsections, false)?;
        }
        result
    }
}

impl DTAExtractor {
    /// The tool body, as the source `main_`.
    fn run_in_pool(ctx: &ToolContext, lines: &mut PoolLines) -> ToolResult {
        let input = ctx.string("in")?;
        let out = ctx.string("out")?.to_owned();
        let (rt, mz, level) = (ctx.string("rt")?, ctx.string("mz")?, ctx.string("level")?);

        // Source initialises both bounds to +/- the largest double, so an
        // omitted side of a range keeps every spectrum on that side.
        let (mut rt_low, mut rt_high) = (-f64::MAX, f64::MAX);
        let (mut mz_low, mut mz_high) = (-f64::MAX, f64::MAX);
        let mut converting = "";
        let converted = (|| -> std::result::Result<Vec<u32>, ()> {
            parse_range(rt, &mut rt_low, &mut rt_high).map_err(|_| ())?;
            ctx.write_debug(
                &format!(
                    "rt lower/upper bound: {} / {}",
                    number(rt_low),
                    number(rt_high)
                ),
                1,
            );
            parse_range(mz, &mut mz_low, &mut mz_high).map_err(|_| ())?;
            ctx.write_debug(
                &format!(
                    "mz lower/upper bound: {} / {}",
                    number(mz_low),
                    number(mz_high)
                ),
                1,
            );
            // Source `tmp = level`, the text its catch block names.
            converting = level;
            // `StringUtils::toInt32` into a `vector<UInt>`: a negative level
            // wraps, as the source's conversion does.
            let levels = level
                .split(',')
                .map(|part| to_int32(part).map(|value| value as u32).map_err(|_| ()))
                .collect::<std::result::Result<Vec<u32>, ()>>()?;
            let listed: Vec<String> = levels.iter().map(u32::to_string).collect();
            ctx.write_debug(&format!("MS levels: {}", listed.join(", ")), 1);
            Ok(levels)
        })();
        let Ok(levels) = converted else {
            lines.error(
                ctx,
                format!("Invalid boundary '{converting}' given. Aborting!"),
            );
            lines.request_usage();
            return Ok(ExitCode::IllegalParameters);
        };

        let experiment = FileHandler::load_experiment(input, &[FileType::MzMl])?;

        for spectrum in &experiment.spectra {
            if !levels.contains(&spectrum.ms_level) {
                continue;
            }
            // The RT filter is applied by the source through the reader's
            // options; applying it here keeps the same set without depending on
            // reader-side filtering.
            if spectrum.rt < rt_low || spectrum.rt > rt_high {
                continue;
            }
            let name = if spectrum.ms_level > 1 {
                let mz = spectrum
                    .precursors
                    .first()
                    .map_or(0.0, |precursor| precursor.mz);
                if mz < mz_low || mz > mz_high {
                    continue;
                }
                format!("{out}_RT{}_MZ{}.dta", number(spectrum.rt), number(mz))
            } else {
                format!("{out}_RT{}.dta", number(spectrum.rt))
            };
            write_dta(&name, spectrum)?;
        }
        Ok(ExitCode::ExecutionOk)
    }
}

/// Source `StringUtils::toStr(double)`, which is what the file name embeds.
///
/// The same port of that function writes the peak m/z inside the file, so the
/// tool has one formatter for one source function: see
/// [`crate::format::dta`]'s numeric-text section.
fn number(value: f64) -> String {
    crate::format::file_info::text_format::to_str(value)
}

/// Source `DTAFile::store` conventions: the legacy proton mass and no guard
/// against discarding what DTA cannot represent.
fn write_dta(path: &str, spectrum: &MSSpectrum) -> Result<()> {
    let file = File::create(path)?;
    dta::write_with_options(BufWriter::new(file), spectrum, &dta::WriteOptions::source())
}
