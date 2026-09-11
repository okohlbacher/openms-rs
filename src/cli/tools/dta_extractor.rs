// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Extracts spectra of an mzML file to DTA files, as the source TOPP tool.
//!
//! Ports `OpenMS4-topp/src/DTAExtractor.cpp`. The tool lives here rather than in
//! `src/bin/` so that the shipped binary and its differential test share one
//! definition; a copied tool body drifts from the binary it claims to test.

use crate::cli::{ExitCode, Tool, ToolContext, ToolSpec, parse_range};
use crate::format::dta;
use crate::format::file_handler::FileHandler;
use crate::format::file_types::FileType;
use crate::kernel::MSSpectrum;
use crate::{Error, Result};
use std::fs::File;
use std::io::BufWriter;

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

    fn run(ctx: &ToolContext) -> Result<ExitCode> {
        let input = ctx.string("in")?;
        let out = ctx.string("out")?.to_owned();

        // Source initialises both bounds to +/- the largest double, so an
        // omitted side of a range keeps every spectrum on that side.
        let (mut rt_low, mut rt_high) = (-f64::MAX, f64::MAX);
        let (mut mz_low, mut mz_high) = (-f64::MAX, f64::MAX);
        parse_range(ctx.string("rt")?, &mut rt_low, &mut rt_high)?;
        parse_range(ctx.string("mz")?, &mut mz_low, &mut mz_high)?;

        let levels: Vec<u32> = ctx
            .string("level")?
            .split(',')
            .map(|part| {
                part.trim()
                    .parse::<u32>()
                    .map_err(|_| Error::InvalidValue(format!("invalid MS level '{part}'")))
            })
            .collect::<Result<_>>()?;

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
                format!("{out}_RT{}_MZ{}.dta", number(spectrum.rt)?, number(mz)?)
            } else {
                format!("{out}_RT{}.dta", number(spectrum.rt)?)
            };
            write_dta(&name, spectrum)?;
        }
        Ok(ExitCode::ExecutionOk)
    }
}

/// Source `StringUtils::toStr(double)`, which is what the file name embeds.
fn number(value: f64) -> Result<String> {
    use crate::data_structures::list::ListFormat;
    Ok(value.to_list_text()?.into_owned())
}

/// Source `DTAFile::store` conventions: the legacy proton mass and no guard
/// against discarding what DTA cannot represent.
fn write_dta(path: &str, spectrum: &MSSpectrum) -> Result<()> {
    let file = File::create(path)?;
    dta::write_with_options(BufWriter::new(file), spectrum, &dta::WriteOptions::source())
}
