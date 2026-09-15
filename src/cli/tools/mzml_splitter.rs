// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Splits an mzML file into parts, as the source TOPP tool.
//!
//! Ports `OpenMS4-topp/src/MzMLSplitter.cpp`. See [`super`] on why the tool
//! lives in the library rather than in the binary.
//!
//! The source attaches its `data filtering` processing record to each part
//! before moving the part's spectra and chromatograms in, so no part carries
//! it; this port does the same.

use crate::cli::{ExitCode, Tool, ToolContext, ToolSpec};
use crate::format::file_handler::FileHandler;
use crate::format::file_types::{FileType, strip_extension};
use crate::kernel::MSExperiment;
use crate::metadata::ProcessingAction;
use crate::{Error, Result};

/// The `MzMLSplitter` TOPP tool.
pub struct MzMLSplitter;

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
            "Prefix for output files ('_part1of2.mzML' etc. will be appended; default: same as 'in' without the file extension)",
            false,
            false,
        )?;
        spec.register_int_option(
            "parts",
            "<num>",
            1,
            "Number of parts to split into (takes precedence over 'size' if set)",
            false,
            false,
        )?;
        spec.set_min_int("parts", 1)?;
        spec.register_int_option(
            "size",
            "<num>",
            0,
            "Approximate upper limit for resulting file sizes (in 'unit')",
            false,
            false,
        )?;
        spec.set_min_int("size", 0)?;
        spec.register_string_option(
            "unit",
            "<choice>",
            "MB",
            "Unit for 'size' (base 1024)",
            false,
            false,
        )?;
        spec.set_valid_strings("unit", &["KB", "MB", "GB"])?;
        spec.register_flag(
            "no_chrom",
            "Remove chromatograms, keep only spectra.",
            false,
        )?;
        spec.register_flag("no_spec", "Remove spectra, keep only chromatograms.", false)?;
        Ok(())
    }

    /// Source `main_`, run on the worker pool that `-threads` sizes, as
    /// `TOPPBase::main` applies the setting before `main_`
    /// (`TOPPBase.cpp:408-415`). See [`ToolContext::in_thread_pool`].
    fn run(ctx: &ToolContext) -> Result<ExitCode> {
        ctx.in_thread_pool(|| Self::run_in_pool(ctx))?
    }
}

impl MzMLSplitter {
    /// The tool body, as the source `main_`.
    fn run_in_pool(ctx: &ToolContext) -> Result<ExitCode> {
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
            // Source divides the byte count as f32 by the unit, then rounds up.
            let bytes = std::fs::metadata(&input)?.len() as f32;
            let total = match ctx.string("unit")? {
                "KB" => bytes / 1024.0,
                "MB" => bytes / (1024.0 * 1024.0),
                _ => bytes / (1024.0 * 1024.0 * 1024.0),
            };
            parts = (total / size as f32).ceil() as usize;
        }

        let experiment = FileHandler::load_experiment(&input, &[FileType::MzMl])?;
        // Source moves the records out of the loaded run, so every part keeps
        // the experiment metadata but no peaks of its own before filling.
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

        // Part numbers are zero padded to the width of the part count.
        let width = parts.to_string().len();
        let (mut spec_start, mut chrom_start) = (0usize, 0usize);
        for counter in 1..=parts {
            let name = format!("{out}_part{counter:0width$}of{parts}.mzML");
            let mut part = MSExperiment {
                spectra: Vec::new(),
                chromatograms: Vec::new(),
                ..template.clone()
            };
            // Source addDataProcessing_(part, getProcessingInfo_(FILTERING)) runs
            // here, while the part holds no spectra and no chromatograms yet,
            // so the record reaches no output. Reproduced; the C++ product SDK
            // writes parts without it (oracle mzml_splitter_1).
            let processing = ctx.processing_info(&[ProcessingAction::DataFiltering])?;
            ctx.add_data_processing(&mut part, &processing);
            // Source spreads the remainder over the parts that are still to come.
            let remaining = parts - counter + 1;
            let n_spec = div_ceil(spectra.len() - spec_start, remaining);
            part.spectra
                .extend_from_slice(&spectra[spec_start..spec_start + n_spec]);
            spec_start += n_spec;
            let n_chrom = div_ceil(chromatograms.len() - chrom_start, remaining);
            part.chromatograms
                .extend_from_slice(&chromatograms[chrom_start..chrom_start + n_chrom]);
            chrom_start += n_chrom;
            FileHandler::store_experiment(&name, &part, Some(FileType::MzMl))?;
        }
        Ok(ExitCode::ExecutionOk)
    }
}

/// Source `ceil(count / double(remaining))`, evaluated in integers.
fn div_ceil(count: usize, remaining: usize) -> usize {
    if remaining == 0 {
        0
    } else {
        count.div_ceil(remaining)
    }
}
