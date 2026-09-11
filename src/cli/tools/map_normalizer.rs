// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Normalizes peak intensities to a percentage of the run maximum.
//!
//! Ports `OpenMS4-topp/src/MapNormalizer.cpp`.

use crate::cli::{ExitCode, Tool, ToolContext, ToolSpec};
use crate::format::file_handler::FileHandler;
use crate::format::file_types::FileType;
use crate::{Error, Result};

/// The `MapNormalizer` TOPP tool.
pub struct MapNormalizer;

impl Tool for MapNormalizer {
    const NAME: &'static str = "MapNormalizer";
    const DESCRIPTION: &'static str = "Normalizes peak intensities in an MS run.";

    fn register(spec: &mut ToolSpec) -> Result<()> {
        spec.register_input_file("in", "<file>", "", "input file ", true, false, &[])?;
        spec.set_valid_formats("in", &["mzML"])?;
        spec.register_output_file("out", "<file>", "", "output file ", true, false)?;
        spec.set_valid_formats("out", &["mzML"])?;
        Ok(())
    }

    fn run(ctx: &ToolContext) -> Result<ExitCode> {
        let mut experiment = FileHandler::load_experiment(ctx.string("in")?, &[FileType::MzMl])?;

        // Source takes the maximum over every MS level, then divides by a
        // hundredth of it, so the most intense peak in the run becomes 100.
        let intensity = experiment
            .ranges(0)?
            .intensity
            .ok_or_else(|| Error::InvalidValue("run has no intensities to normalize".into()))?;
        let scale = intensity.max / 100.0;
        if !(scale.is_finite() && scale > 0.0) {
            return Err(Error::InvalidValue(
                "run maximum intensity must be finite and positive".into(),
            ));
        }

        // Only MS1 is scaled; the source leaves higher levels untouched and
        // its commented-out chromatogram branch is not ported.
        for spectrum in &mut experiment.spectra {
            if spectrum.ms_level < 2 {
                for peak in &mut spectrum.peaks {
                    peak.intensity = (f64::from(peak.intensity) / scale) as f32;
                }
            }
        }

        FileHandler::store_experiment(ctx.string("out")?, &experiment, Some(FileType::MzMl))?;
        Ok(ExitCode::ExecutionOk)
    }
}
