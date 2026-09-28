// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Feature detection for LC-MS1 data with ion mobility support, as the source
//! TOPP tool `FeatureFinderLFQ` (EXPERIMENTAL).
//!
//! Ports `OpenMS4-topp/src/FeatureFinderLFQ.cpp` (topp `174b576`), the TOPP
//! wrapper of [`Biosaur2Algorithm`], whose documentation reads: *This TOPP tool
//! is a C++ reimplementation of the Biosaur2 feature detection algorithm. It
//! detects peptide features in centroided LC-MS1 data (with optional profile
//! mode support) through grouping peaks across scans into "hills", splitting
//! hills at valley points, detecting isotope patterns and calculating feature
//! properties.* It supports FAIMS compensation-voltage grouping, PASEF/TIMS
//! ion mobility, TOF intensity filtering, automatic mass calibration and
//! profile data, and exports featureXML plus Biosaur2-compatible TSV files.
//! *The tool is experimental and under active development; the interface and
//! behavior may change.* Reference: Abdrakhimov *et al.*, Rapid Communications
//! in Mass Spectrometry, 2022, <https://doi.org/10.1002/rcm.9045>.
//!
//! # Run order
//!
//! `main_` (`FeatureFinderLFQ.cpp:117-210`), step by step:
//!
//! 1. Copy the algorithm parameters, registered at the top level with
//!    `registerFullParam_` (`-mini`, `-htol`, ...), into the algorithm.
//! 2. Load `-in` with MS level 1 only.
//! 3. Run the algorithm, keeping the hills and the peptide features.
//! 4. Annotate the primary MS run path (`file://<basename>` under `-test`, or
//!    the input's own mzML path when it names an existing file) and the
//!    `Quantitation` processing record.
//! 5. Store `-out` as featureXML.
//! 6. With `-write_hills` or `-out_hills`, write the hills TSV: to `-out_hills`,
//!    or next to `-out` with everything from its last `.` replaced by
//!    `.hills.tsv` (the source's `find_last_of('.')` looks at the whole path).
//! 7. With `-out_tsv`, write the Biosaur2 feature TSV.
//!
//! The algorithm's `OPENMS_LOG_INFO` lines go to standard output and its
//! `OPENMS_LOG_WARN` lines to standard error, in the order the serial source
//! writes them; see [`Biosaur2Message`](crate::analysis::featurefinder::biosaur2::Biosaur2Message)
//! for the timing suffixes it leaves out. The source's `ProgressLogger`
//! sections (`Loading input mzML`, `Preprocessing and feature finding`, ...)
//! and its four timing lines (`Loaded input file in ...`, `... took ...`) are
//! not written: their text is clock-dependent and no upstream test reads it.
//!
//! An algorithm failure (for example `No MS1 spectra found in input
//! experiment!`) is the source's uncaught `Exception::InvalidValue`, which
//! `TOPPBase` reports as `Error: Unexpected internal error (<what>)` with exit
//! code 8; a TSV that cannot be written is its `UnableToCreateFile` (exit 5).
//!
//! `-in` accepts mzML only: the source adds Bruker `.d` and Thermo `raw` when
//! built with OpenTIMS or the Thermo reader, and the Release build this port
//! is compared with has neither (its `-write_ini` lists `*.mzML` only).

use crate::analysis::featurefinder::biosaur2::{Biosaur2Algorithm, Biosaur2LogLevel, output};
use crate::cli::{ExitCode, Tool, ToolContext, ToolError, ToolResult, ToolSpec};
use crate::format::file_handler::FileHandler;
use crate::format::file_types::FileType;
use crate::format::mzml::ReadOptions;
use crate::format::peak_options::PeakFileOptions;
use crate::metadata::ProcessingAction;
use crate::system::file;
use crate::{Error, Result};
use std::io::Write;
use std::path::Path;

/// The `FeatureFinderLFQ` TOPP tool.
pub struct FeatureFinderLFQ;

impl FeatureFinderLFQ {
    /// The source `PeakFileOptions` of `main_` (`FeatureFinderLFQ.cpp:138-141`):
    /// MS level 1 only, *for feature finding*.
    ///
    /// # Errors
    ///
    /// Propagates the option container's MS-level limit error, which one level
    /// cannot reach.
    pub fn peak_file_options() -> Result<PeakFileOptions> {
        let mut options = PeakFileOptions::default();
        options.add_ms_level(1)?;
        Ok(options)
    }

    /// The hills TSV path of step 6: `out_hills` when given, otherwise `out`
    /// up to its last `.` (anywhere in the path, as the source's
    /// `find_last_of('.')`), or all of `out` when it has none, followed by
    /// `.hills.tsv`.
    pub fn hills_path(out: &str, out_hills: &str) -> String {
        if !out_hills.is_empty() {
            return out_hills.to_owned();
        }
        let base = out.rfind('.').map_or(out, |dot| &out[..dot]);
        format!("{base}.hills.tsv")
    }
}

impl Tool for FeatureFinderLFQ {
    const NAME: &'static str = "FeatureFinderLFQ";
    const DESCRIPTION: &'static str = "Feature detection for LC-MS1 data (EXPERIMENTAL)";

    /// Source `registerOptionsAndFlags_` (`FeatureFinderLFQ.cpp:83-109`): the
    /// files, the `-write_hills` flag and every algorithm parameter at the top
    /// level.
    fn register(spec: &mut ToolSpec) -> Result<()> {
        spec.register_input_file(
            "in",
            "<file>",
            "",
            "Input mzML file (centroided data)",
            true,
            false,
            &[],
        )?;
        spec.set_valid_formats("in", &["mzML"])?;
        spec.register_output_file("out", "<file>", "", "Output featureXML file", true, false)?;
        spec.set_valid_formats("out", &["featureXML"])?;
        spec.register_output_file(
            "out_tsv",
            "<file>",
            "",
            "Optional: output TSV file (Biosaur2 format)",
            false,
            false,
        )?;
        spec.set_valid_formats("out_tsv", &["tsv"])?;
        spec.register_output_file(
            "out_hills",
            "<file>",
            "",
            "Optional: write detected hills to TSV",
            false,
            false,
        )?;
        spec.set_valid_formats("out_hills", &["tsv"])?;
        spec.register_flag(
            "write_hills",
            "Force writing of hills file even if no output path was provided",
            false,
        )?;
        spec.register_full_param(&Biosaur2Algorithm::defaults()?)?;
        Ok(())
    }

    /// Forwards to [`run_io`](Tool::run_io) with the process streams.
    fn run(ctx: &ToolContext) -> ToolResult {
        Self::run_io(ctx, &mut std::io::stdout(), &mut std::io::stderr())
    }

    /// Source `main_`; the module documentation lists the steps.
    fn run_io(ctx: &ToolContext, out: &mut dyn Write, err: &mut dyn Write) -> ToolResult {
        let input = ctx.string("in")?.to_owned();
        let output_path = ctx.string("out")?.to_owned();
        let out_tsv = ctx.string("out_tsv")?.to_owned();
        let out_hills = ctx.string("out_hills")?.to_owned();
        let write_hills = ctx.flag("write_hills")?;

        // `getParam_().copySubset(algorithm_.getDefaults())` (128-129).
        let parameters = ctx.param().copy_subset(&Biosaur2Algorithm::defaults()?)?;
        let mut algorithm = Biosaur2Algorithm::from_param(&parameters)?;
        algorithm.threads = ctx.thread_policy();

        // Load MS1 spectra (135-146). The source FileHandler tolerates
        // dangling mzML header references; the native reader does only on
        // request (decision D10).
        let read = ReadOptions {
            source_dangling_references: true,
            ..ReadOptions::default()
        };
        let experiment = FileHandler::load_experiment_with_read_options(
            &input,
            &[FileType::MzMl],
            &Self::peak_file_options()?,
            &read,
        )?;
        algorithm.set_ms_data(experiment);

        // The algorithm (148-157). Its exceptions reach TOPPBase's
        // BaseException arm.
        let mut generator = ctx.unique_id_generator();
        let result = algorithm
            .run_with_details(&mut generator)
            .map_err(|error| match error {
                Error::InvalidValue(what) => ToolError::unexpected(what),
                other => ToolError::unexpected(other),
            })?;
        for message in &result.messages {
            match message.level {
                Biosaur2LogLevel::Info => writeln!(out, "{}", message.text)?,
                Biosaur2LogLevel::Warning => crate::cli::log_warning(err, &message.text)?,
            }
        }

        // Annotation (158-160).
        let mut feature_map = result.feature_map;
        let primary_path = if ctx.flag("test")? {
            format!("file://{}", file::basename(&input))
        } else {
            input.clone()
        };
        feature_map
            .set_primary_ms_run_path_from_experiment(&[primary_path], algorithm.ms_data())?;
        let processing = ctx.processing_info(&[ProcessingAction::Quantitation])?;
        ctx.add_data_processing(&mut feature_map, &processing);

        // The featureXML (163-172).
        FileHandler::store_feature_map(&output_path, &feature_map, Some(FileType::FeatureXml))
            .map_err(|error| crate::cli::write_failure(&output_path, &error))?;
        writeln!(
            out,
            "Wrote {} features to: {output_path}",
            result.peptide_features.len()
        )?;

        // The hills TSV (174-195).
        if write_hills || !out_hills.is_empty() {
            let hills_file = Self::hills_path(&output_path, &out_hills);
            let message = output::write_hills(&result.hills, Path::new(&hills_file))
                .map_err(|error| crate::cli::write_failure(&hills_file, &error))?;
            writeln!(out, "{}", message.text)?;
        }

        // The feature TSV (197-207).
        if !out_tsv.is_empty() {
            let message = output::write_tsv(&result.peptide_features, Path::new(&out_tsv))
                .map_err(|error| crate::cli::write_failure(&out_tsv, &error))?;
            writeln!(out, "{}", message.text)?;
        }
        Ok(ExitCode::ExecutionOk)
    }
}

#[cfg(test)]
mod tests {
    use super::FeatureFinderLFQ;

    #[test]
    fn the_hills_path_follows_the_last_dot() {
        assert_eq!(
            FeatureFinderLFQ::hills_path("a/b.featureXML", ""),
            "a/b.hills.tsv"
        );
        assert_eq!(
            FeatureFinderLFQ::hills_path("a/b.featureXML", "h.tsv"),
            "h.tsv"
        );
        assert_eq!(
            FeatureFinderLFQ::hills_path("dir.x/out", ""),
            "dir.hills.tsv"
        );
        assert_eq!(FeatureFinderLFQ::hills_path("out", ""), "out.hills.tsv");
    }
}
