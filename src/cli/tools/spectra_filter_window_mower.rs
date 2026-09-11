// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Retains the most intense peaks within m/z windows.
//!
//! Ports `OpenMS4-topp/src/SpectraFilterWindowMower.cpp`. The tool's `algorithm`
//! subsection carries the `WindowMower` parameters, exactly as the source
//! `getSubsectionDefaults_` supplies them.

use crate::cli::{ExitCode, Tool, ToolContext, ToolSpec};
use crate::format::file_handler::FileHandler;
use crate::format::file_types::FileType;
use crate::param::{Param, ParamValue};
use crate::processing::SpectrumFilter;
use crate::processing::window_mower::{WindowMower, WindowMowerMethod};
use crate::{Error, Result};
use std::str::FromStr;

/// The `SpectraFilterWindowMower` TOPP tool.
pub struct SpectraFilterWindowMower;

impl Tool for SpectraFilterWindowMower {
    const NAME: &'static str = "SpectraFilterWindowMower";
    const DESCRIPTION: &'static str = "Applies thresholdfilter to peak spectra.";

    fn register(spec: &mut ToolSpec) -> Result<()> {
        spec.register_input_file("in", "<file>", "", "input file ", true, false, &[])?;
        spec.set_valid_formats("in", &["mzML"])?;
        spec.register_output_file("out", "<file>", "", "output file ", true, false)?;
        spec.set_valid_formats("out", &["mzML"])?;
        spec.register_subsection("algorithm", "Algorithm parameter subsection.")?;
        Ok(())
    }

    /// The source `WindowMower` defaults, as `WindowMower::WindowMower` sets them.
    fn subsection_defaults(_section: &str) -> Result<Option<Param>> {
        let mut defaults = Param::new();
        let filter = WindowMower::default();
        defaults.set_value(
            "windowsize",
            ParamValue::Float(filter.window_size),
            "The size of the sliding window along the m/z axis.",
            &[],
        )?;
        defaults.set_value(
            "peakcount",
            ParamValue::Integer(filter.peak_count as i64),
            "The number of peaks that should be kept.",
            &[],
        )?;
        defaults.set_value(
            "movetype",
            ParamValue::String(
                match filter.method {
                    WindowMowerMethod::Sliding => "slide",
                    WindowMowerMethod::Jumping => "jump",
                }
                .into(),
            ),
            "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.",
            &[],
        )?;
        defaults.set_valid_strings("movetype", &["slide".into(), "jump".into()])?;
        Ok(Some(defaults))
    }

    fn run(ctx: &ToolContext) -> Result<ExitCode> {
        let algorithm = ctx.subsection("algorithm")?;
        let value = |key: &str| -> Result<&ParamValue> {
            algorithm
                .value(key)
                .map_err(|_| Error::InvalidValue(format!("algorithm:{key} is not set")))
        };
        let window_size = match value("windowsize")? {
            ParamValue::Float(v) => *v,
            ParamValue::Integer(v) => *v as f64,
            _ => {
                return Err(Error::InvalidValue(
                    "algorithm:windowsize must be a number".into(),
                ));
            }
        };
        let peak_count = match value("peakcount")? {
            ParamValue::Integer(v) => usize::try_from(*v).map_err(|_| {
                Error::InvalidValue("algorithm:peakcount must not be negative".into())
            })?,
            _ => {
                return Err(Error::InvalidValue(
                    "algorithm:peakcount must be an integer".into(),
                ));
            }
        };
        let method = WindowMowerMethod::from_str(value("movetype")?.as_str()?)?;

        let filter = WindowMower {
            window_size,
            peak_count,
            method,
            ..WindowMower::default()
        };

        let mut experiment = FileHandler::load_experiment(ctx.string("in")?, &[FileType::MzMl])?;
        filter.filter_experiment(&mut experiment)?;
        FileHandler::store_experiment(ctx.string("out")?, &experiment, Some(FileType::MzMl))?;
        Ok(ExitCode::ExecutionOk)
    }
}
