// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Removes the baseline from profile spectra by morphological filtering.
//!
//! Ports `OpenMS4-topp/src/BaselineFilter.cpp`. Unlike most algorithm-wrapping
//! tools this registers the filter's three parameters as ordinary options
//! rather than a subsection, exactly as the source does. The output carries the
//! source's `baseline reduction` processing record.
//!
//! The source's input checks are `main_`'s own (`BaselineFilter.cpp:104-125`):
//! an input without spectra is an `OPENMS_LOG_WARN` warning and
//! `INCOMPATIBLE_INPUT_DATA`, a first spectrum that peak type estimation calls
//! centroided a `writeLogWarn_` warning, and an unsorted spectrum a
//! `writeLogError_` line and `INCOMPATIBLE_INPUT_DATA`. Each code is one
//! `main_` returns, so `TOPPBase`'s closing line follows (Release oracles
//! `bf_empty_log` and `bf_centroided_log` of `../oracle/topp-exception-exits`).

use crate::cli::{ExitCode, PoolLines, Tool, ToolContext, ToolResult, ToolSpec};
use crate::format::file_handler::FileHandler;
use crate::format::file_types::FileType;
use crate::kernel::{MSSpectrum, SpectrumType, SpectrumTypeQueryLimits};
use crate::metadata::ProcessingAction;
use crate::processing::SpectrumFilter;
use crate::processing::baseline::{MorphologicalFilter, MorphologicalMethod, StructuringElement};
use crate::{Error, Result};
use std::io::Write;

/// The `BaselineFilter` TOPP tool.
pub struct BaselineFilter;

/// Source `method` names, in the order `setValidStrings_` lists them.
///
/// `erosion_simple` and `dilation_simple` select the source's direct-window
/// `applyErosionSimple_` / `applyDilationSimple_`, which the source reaches
/// without the van Herk case distinctions of `erosion` and `dilation`. The two
/// pairs agree except at the last sample of a spectrum filtered with a
/// one-sample element, so each name maps to its own method here.
fn method(name: &str) -> Result<MorphologicalMethod> {
    Ok(match name {
        "identity" => MorphologicalMethod::Identity,
        "erosion" => MorphologicalMethod::Erosion,
        "erosion_simple" => MorphologicalMethod::ErosionSimple,
        "dilation" => MorphologicalMethod::Dilation,
        "dilation_simple" => MorphologicalMethod::DilationSimple,
        "opening" => MorphologicalMethod::Opening,
        "closing" => MorphologicalMethod::Closing,
        "gradient" => MorphologicalMethod::Gradient,
        "tophat" => MorphologicalMethod::TopHat,
        "bothat" => MorphologicalMethod::BottomHat,
        other => {
            return Err(Error::InvalidValue(format!(
                "unknown morphological method '{other}'"
            )));
        }
    })
}

impl Tool for BaselineFilter {
    const NAME: &'static str = "BaselineFilter";
    const DESCRIPTION: &'static str =
        "Removes the baseline from profile spectra using a top-hat filter.";

    fn register(spec: &mut ToolSpec) -> Result<()> {
        spec.register_input_file("in", "<file>", "", "input raw data file ", true, false, &[])?;
        spec.set_valid_formats("in", &["mzML"])?;
        spec.register_output_file("out", "<file>", "", "output raw data file ", true, false)?;
        spec.set_valid_formats("out", &["mzML"])?;
        spec.register_double_option(
            "struc_elem_length",
            "<size>",
            3.0,
            "Length of the structuring element (should be wider than maximal peak width - see documentation).",
            false,
            false,
        )?;
        spec.register_string_option(
            "struc_elem_unit",
            "<unit>",
            "Thomson",
            "Unit of 'struc_elem_length' parameter.",
            false,
            false,
        )?;
        spec.set_valid_strings("struc_elem_unit", &["Thomson", "DataPoints"])?;
        spec.register_string_option(
            "method",
            "<string>",
            "tophat",
            "The name of the morphological filter to be applied. If you are unsure, use the default.",
            false,
            false,
        )?;
        spec.set_valid_strings(
            "method",
            &[
                "identity",
                "erosion",
                "dilation",
                "opening",
                "closing",
                "gradient",
                "tophat",
                "bothat",
                "erosion_simple",
                "dilation_simple",
            ],
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
    /// body's console lines are written once the pool returns.
    fn run_io(ctx: &ToolContext, out: &mut dyn Write, err: &mut dyn Write) -> ToolResult {
        let mut lines = PoolLines::default();
        let result = ctx.in_thread_pool(|| Self::run_in_pool(ctx, &mut lines))?;
        lines.write(out, err)?;
        result
    }
}

/// Source `ms_exp[0].getType(true)` (`BaselineFilter.cpp:112`): the stored
/// type, else a `PEAK_PICKING` record in the processing history, else the peak
/// type estimation.
///
/// It only decides whether the source warns, so it never ends the run. The
/// query's ceilings are the spectrum's own size: the spectrum is in memory
/// already, and the estimation copies its two value arrays once, so no ceiling
/// refuses a spectrum the reader admitted (the library default stops at a
/// million points). A spectrum holding a non-finite value, which the native
/// estimator declines to classify, is not warned about; the source classifies
/// it by whatever its arithmetic on the value gives.
fn first_spectrum_type(spectrum: &MSSpectrum) -> SpectrumType {
    let points = spectrum.peaks.len();
    let records = spectrum.data_processing.len();
    let limits = SpectrumTypeQueryLimits {
        max_points: points,
        // 32 units per point for the estimation, and at most 1 + 12 * 64 per
        // history record (`MSSpectrum::get_type_with_limits`).
        max_work: points
            .saturating_mul(32)
            .saturating_add(records.saturating_mul(1 + 12 * 64)),
        max_bytes: points.saturating_mul(2 * std::mem::size_of::<f64>()),
    };
    spectrum
        .get_type_with_limits(true, limits)
        .unwrap_or(SpectrumType::Unknown)
}

/// Source warning for an input without spectra (`BaselineFilter.cpp:107-108`);
/// the warning log stream ends the line.
const EMPTY_INPUT_WARNING: &str = "The given file does not contain any conventional peak data, but might contain chromatograms. This tool currently cannot handle them, sorry.";

impl BaselineFilter {
    /// The tool body, as the source `main_`.
    fn run_in_pool(ctx: &ToolContext, lines: &mut PoolLines) -> ToolResult {
        let mut experiment = FileHandler::load_experiment(ctx.string("in")?, &[FileType::MzMl])?;

        // Source refuses a run that carries only chromatograms, and refuses
        // unsorted spectra rather than producing a wrong baseline.
        if experiment.spectra.is_empty() {
            lines.console_warning(EMPTY_INPUT_WARNING);
            return Ok(ExitCode::IncompatibleInputData);
        }
        if first_spectrum_type(&experiment.spectra[0]) == SpectrumType::Centroid {
            lines.warn(
                ctx,
                "Warning: OpenMS peak type estimation indicates that this is not raw data!",
            );
        }
        if !experiment
            .spectra
            .iter()
            .all(|spectrum| spectrum.is_sorted())
        {
            lines.error(
                ctx,
                "Error: Not all spectra are sorted according to peak m/z positions. Use FileFilter to sort the input!",
            );
            return Ok(ExitCode::IncompatibleInputData);
        }

        let length = ctx.double("struc_elem_length")?;
        let filter = MorphologicalFilter {
            method: method(ctx.string("method")?)?,
            structuring_element: match ctx.string("struc_elem_unit")? {
                "Thomson" => StructuringElement::Thomson(length),
                "DataPoints" => {
                    StructuringElement::DataPoints(usize::try_from(length as i64).map_err(
                        |_| Error::InvalidValue("struc_elem_length must not be negative".into()),
                    )?)
                }
                other => {
                    return Err(Error::InvalidValue(format!(
                        "unknown structuring element unit '{other}'"
                    ))
                    .into());
                }
            },
        };
        filter.filter_experiment(&mut experiment)?;
        // Source addDataProcessing_(ms_exp, getProcessingInfo_(BASELINE_REDUCTION)):
        // one shared record on every spectrum and chromatogram.
        let processing = ctx.processing_info(&[ProcessingAction::BaselineReduction])?;
        ctx.add_data_processing(&mut experiment, &processing);
        FileHandler::store_experiment(ctx.string("out")?, &experiment, Some(FileType::MzMl))?;
        Ok(ExitCode::ExecutionOk)
    }
}
