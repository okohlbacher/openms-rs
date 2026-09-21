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
//! Not ported: the warning the source writes when peak type estimation finds
//! the first spectrum centroided.

use crate::cli::{ExitCode, Tool, ToolContext, ToolResult, ToolSpec};
use crate::format::file_handler::FileHandler;
use crate::format::file_types::FileType;
use crate::metadata::ProcessingAction;
use crate::processing::SpectrumFilter;
use crate::processing::baseline::{MorphologicalFilter, MorphologicalMethod, StructuringElement};
use crate::{Error, Result};

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

    /// Source `main_`, run on the worker pool that `-threads` sizes, as
    /// `TOPPBase::main` applies the setting before `main_`
    /// (`TOPPBase.cpp:408-415`). See [`ToolContext::in_thread_pool`].
    fn run(ctx: &ToolContext) -> ToolResult {
        Ok(ctx.in_thread_pool(|| Self::run_in_pool(ctx))??)
    }
}

impl BaselineFilter {
    /// The tool body, as the source `main_`.
    fn run_in_pool(ctx: &ToolContext) -> Result<ExitCode> {
        let mut experiment = FileHandler::load_experiment(ctx.string("in")?, &[FileType::MzMl])?;

        // Source refuses a run that carries only chromatograms, and refuses
        // unsorted spectra rather than producing a wrong baseline.
        if experiment.spectra.is_empty() {
            return Err(Error::Unsupported(
                "the given file contains no conventional peak data; chromatograms are not handled"
                    .into(),
            ));
        }
        for (index, spectrum) in experiment.spectra.iter().enumerate() {
            if !spectrum.is_sorted() {
                return Err(Error::Unsupported(format!(
                    "spectrum {index} is not sorted by m/z; sort the input with FileFilter first"
                )));
            }
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
                    )));
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
