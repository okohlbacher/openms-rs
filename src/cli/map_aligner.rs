// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The shared parameters and checks of the map-aligner tools: the native form
//! of the `OpenMS4-cli` package's `APPLICATIONS/MapAlignerBase.h` and
//! `MapAlignerBase.cpp` (cli `c19e494`).
//!
//! The header declares two types. `TOPPMapAlignerBase`, the base class of
//! `MapAlignerIdentification`, `MapAlignerPoseClustering`,
//! `MapAlignerTreeGuided` and `MapRTTransformer`, registers their common
//! parameters, checks them, and transforms optional spectra files along with
//! the alignment; it is [`TOPPMapAlignerBase`] here. `MapAlignerBase` holds
//! one static function, `getModelDefaults`, which forwards to core
//! `TransformationModelDefaults::getDefaults`; that core header is not in this
//! tree yet, so `getModelDefaults` is not ported (see
//! `docs/TOPP_CLI_SUPPORT.md`, *Tool base classes*).
//!
//! The source header wraps both in Doxygen's `@cond TOPPCLASSES`, so they are
//! absent from its API documentation; its only page text is "Base class for
//! different MapAligner TOPP tools".

use super::context::ToolContext;
use super::parameter::ExitCode;
use super::spec::ToolSpec;
use crate::format::file_handler::FileHandler;
use crate::{Error, Result};
use std::io::Write;

/// Which reference parameters a map-aligner tool offers (source enum
/// `TOPPMapAlignerBase::ReferenceParameterKind`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum ReferenceParameterKind {
    /// No reference (source `REF_NONE`, the default the constructor sets).
    #[default]
    None,
    /// The reference file must have the same type as the input files (source
    /// `REF_RESTRICTED`).
    Restricted,
    /// The reference file can have any supported file type (source
    /// `REF_FLEXIBLE`).
    Flexible,
}

/// The common part of the map-aligner tools (source class
/// `TOPPMapAlignerBase`).
///
/// The source keeps the kind of reference parameters in a member,
/// `ref_params_`, which its constructor sets to `REF_NONE` and
/// `registerOptionsAndFlagsMapAligners_` sets to the tool's kind, so that
/// `checkParameters_` can read it later on the same object. A native tool's
/// registration and its run are separate functions with no object between
/// them, so the kind is given when this value is made
/// ([`new`](Self::new)), typically as a constant of the tool, and both
/// [`register_options_and_flags_map_aligners`](Self::register_options_and_flags_map_aligners)
/// and [`check_parameters`](Self::check_parameters) read it. The source's
/// constructor arguments (name, description, `official`) belong to the tool's
/// [`Tool`](super::Tool) implementation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[allow(clippy::upper_case_acronyms)] // the source class name
pub struct TOPPMapAlignerBase {
    ref_params: ReferenceParameterKind,
}

impl TOPPMapAlignerBase {
    /// The map-aligner base of a tool offering `ref_params`.
    pub const fn new(ref_params: ReferenceParameterKind) -> Self {
        Self { ref_params }
    }

    /// The kind of reference parameters (source member `ref_params_`).
    pub const fn ref_params(&self) -> ReferenceParameterKind {
        self.ref_params
    }

    /// Register the parameters every map-aligner tool shares (source
    /// `registerOptionsAndFlagsMapAligners_`, `MapAlignerBase.cpp:26-56`).
    ///
    /// `file_formats` are the formats of the input files, which the source
    /// takes as one comma-separated text and splits with `ListUtils::create`.
    /// In registration order:
    ///
    /// * `in`, a required input file list in `file_formats`;
    /// * `out`, an optional output file list in `file_formats`, and
    ///   `trafo_out`, an optional output file list of `trafoXML` files — one of
    ///   the two has to be given, which [`check_parameters`](Self::check_parameters)
    ///   checks;
    /// * `in_spectra_files` and `out_spectra_files`, optional `mzML` lists
    ///   transformed along with the alignment;
    /// * unless the kind is [`ReferenceParameterKind::None`], the section
    ///   `reference` with `reference:file`, an input file in `file_formats`
    ///   whose description adds ` (same file format as input files required)`
    ///   for [`ReferenceParameterKind::Restricted`], and `reference:index`, an
    ///   integer of at least 0, default 0.
    ///
    /// # Errors
    ///
    /// Propagates a registration error of [`ToolSpec`], such as a name the
    /// tool has already registered.
    pub fn register_options_and_flags_map_aligners(
        &self,
        spec: &mut ToolSpec,
        file_formats: &[&str],
    ) -> Result<()> {
        spec.register_input_file_list(
            "in",
            "<files>",
            &[],
            "Input files to align (all must have the same file type)",
            true,
            false,
            &[],
        )?;
        spec.set_valid_formats("in", file_formats)?;
        spec.register_output_file_list(
            "out",
            "<files>",
            &[],
            "Output files (same file type as 'in'). This option or 'trafo_out' has to be provided; they can be used together.",
            false,
            false,
        )?;
        spec.set_valid_formats("out", file_formats)?;
        spec.register_output_file_list(
            "trafo_out",
            "<files>",
            &[],
            "Transformation output files. This option or 'out' has to be provided; they can be used together.",
            false,
            false,
        )?;
        spec.set_valid_formats("trafo_out", &["trafoXML"])?;

        spec.register_input_file_list(
            "in_spectra_files",
            "<files>",
            &[],
            "Optional input spectra files (mzML) that will be transformed along with the alignment. Size must match the number of input files.",
            false,
            false,
            &[],
        )?;
        spec.set_valid_formats("in_spectra_files", &["mzML"])?;
        spec.register_output_file_list(
            "out_spectra_files",
            "<files>",
            &[],
            "Optional output spectra files (mzML) corresponding to transformed in_spectra_files. Size must match in_spectra_files.",
            false,
            false,
        )?;
        spec.set_valid_formats("out_spectra_files", &["mzML"])?;

        if self.ref_params != ReferenceParameterKind::None {
            spec.register_topp_subsection(
                "reference",
                "Options to define a reference file (use either 'file' or 'index', not both)",
            )?;
            let mut description = String::from("File to use as reference");
            if self.ref_params == ReferenceParameterKind::Restricted {
                description.push_str(" (same file format as input files required)");
            }
            spec.register_input_file(
                "reference:file",
                "<file>",
                "",
                &description,
                false,
                false,
                &[],
            )?;
            spec.set_valid_formats("reference:file", file_formats)?;
            spec.register_int_option(
                "reference:index",
                "<number>",
                0,
                "Use one of the input files as reference ('1' for the first file, etc.).\nIf '0', no explicit reference is set - the algorithm will select a reference.",
                false,
                false,
            )?;
            spec.set_min_int("reference:index", 0)?;
        }
        Ok(())
    }

    /// Check the registered parameters against each other (source
    /// `checkParameters_`, `MapAlignerBase.cpp:58-139`).
    ///
    /// Returns [`ExitCode::ExecutionOk`] when the parameters fit, and
    /// otherwise [`ExitCode::IllegalParameters`] after writing the first
    /// failed check's line through
    /// [`write_log_error`](ToolContext::write_log_error), in the source's
    /// order:
    ///
    /// 1. neither `out` nor `trafo_out` is given: `Error: Data output or
    ///    transformation output files have to be provided (parameters
    ///    'out'/'trafo_out')`;
    /// 2. `out` is given with another length than `in`: `Error: The number of
    ///    data input and output files has to be equal (parameters 'in'/'out')`;
    /// 3. `trafo_out` likewise: `Error: The number of data input and
    ///    transformation output files has to be equal (parameters
    ///    'in'/'trafo_out')`;
    /// 4. the input files are not all of the first one's type
    ///    ([`FileHandler::get_type`]): `Error: All input files (parameter
    ///    'in') must have the same format!`;
    /// 5. spectra files are given and `in_spectra_files` differs in length
    ///    from `in`, or `out_spectra_files` from `in_spectra_files`, with the
    ///    source's two messages;
    /// 6. with reference parameters: `reference:index` above the number of
    ///    input files, both `reference:index` and `reference:file` given, and,
    ///    for [`ReferenceParameterKind::Restricted`], a reference file of
    ///    another type than the inputs, with the source's three messages.
    ///
    /// A tool that gets anything but [`ExitCode::ExecutionOk`] returns it, as
    /// the source's tools do.
    ///
    /// # Errors
    ///
    /// The accessor errors of [`ToolContext`] (a parameter this tool did not
    /// register), [`FileHandler::get_type`]'s errors, and [`Error::Io`] when
    /// `err` fails. An empty `in` is [`Error::InvalidValue`]: the lifecycle
    /// never lets a required list be empty, and the source indexes its first
    /// element unchecked.
    pub fn check_parameters(&self, ctx: &ToolContext, err: &mut dyn Write) -> Result<ExitCode> {
        let ins = ctx.string_list("in")?;
        let outs = ctx.string_list("out")?;
        let trafos = ctx.string_list("trafo_out")?;
        let in_spectra = ctx.string_list("in_spectra_files")?;
        let out_spectra = ctx.string_list("out_spectra_files")?;

        let refuse = |err: &mut dyn Write, text: &str| -> Result<ExitCode> {
            ctx.write_log_error(err, text)?;
            Ok(ExitCode::IllegalParameters)
        };

        if outs.is_empty() && trafos.is_empty() {
            return refuse(
                err,
                "Error: Data output or transformation output files have to be provided (parameters 'out'/'trafo_out')",
            );
        }
        if !outs.is_empty() && ins.len() != outs.len() {
            return refuse(
                err,
                "Error: The number of data input and output files has to be equal (parameters 'in'/'out')",
            );
        }
        if !trafos.is_empty() && ins.len() != trafos.len() {
            return refuse(
                err,
                "Error: The number of data input and transformation output files has to be equal (parameters 'in'/'trafo_out')",
            );
        }
        let first = ins
            .first()
            .ok_or_else(|| Error::InvalidValue("parameter 'in' names no input file".into()))?;
        let in_type = FileHandler::get_type(first)?;
        for input in &ins[1..] {
            if FileHandler::get_type(input)? != in_type {
                return refuse(
                    err,
                    "Error: All input files (parameter 'in') must have the same format!",
                );
            }
        }

        if !in_spectra.is_empty() || !out_spectra.is_empty() {
            if in_spectra.len() != ins.len() {
                return refuse(
                    err,
                    "Error: The number of spectra input files has to be equal to the number of main input files (parameters 'in_spectra_files'/'in')",
                );
            }
            if out_spectra.len() != in_spectra.len() {
                return refuse(
                    err,
                    "Error: The number of spectra input and output spectra files has to be equal (parameters 'in_spectra_files'/'out_spectra_files')",
                );
            }
        }

        if self.ref_params != ReferenceParameterKind::None {
            let reference_index = ctx.int("reference:index")?;
            let reference_file = ctx.string("reference:file")?;
            // The source reads the index into a `Size`; the registered minimum
            // of 0 keeps it non-negative, and a negative value, which only a
            // direct caller could supply, is taken as the huge unsigned value
            // the source's conversion gives it.
            let too_high = usize::try_from(reference_index).map_or(true, |index| index > ins.len());
            if too_high {
                return refuse(
                    err,
                    "Error: Value of parameter 'reference:index' must not be higher than the number of input files",
                );
            }
            if reference_index != 0 && !reference_file.is_empty() {
                return refuse(
                    err,
                    "Error: Parameters 'reference:index' and 'reference:file' cannot be used together",
                );
            }
            if self.ref_params == ReferenceParameterKind::Restricted
                && !reference_file.is_empty()
                && FileHandler::get_type(reference_file)? != in_type
            {
                return refuse(
                    err,
                    "Error: Reference file must have the same format as other input files (parameters 'reference:file'/'in')",
                );
            }
        }
        Ok(ExitCode::ExecutionOk)
    }

    /// Transform spectra files along with the alignment (source
    /// `transformSpectraFiles_`, `MapAlignerBase.cpp:141-174`).
    ///
    /// Nothing happens when either list is empty. Otherwise, for every index
    /// `i`, the mzML file `in_spectra_files[i]` is loaded, its retention times
    /// are transformed by `transformations[i]`
    /// ([`MapAlignmentTransformer::transform_experiment`](crate::analysis::alignment_transformer::MapAlignmentTransformer::transform_experiment),
    /// keeping each original retention time as meta value when
    /// `store_original_rt` is set), it gains the tool's
    /// `Retention time alignment` processing record, and it is stored as mzML
    /// under `out_spectra_files[i]`.
    ///
    /// The source reports progress (`transforming spectra files`) through a
    /// `ProgressLogger` of the tool's log type; progress output is not ported
    /// for the tools of this crate yet, so nothing is reported here.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the three lists differ in length, before
    /// anything is read. The source states that as two
    /// `OPENMS_PRECONDITION`s, which a Release build does not check, and then
    /// indexes the lists unchecked. Otherwise the load, transformation,
    /// processing-record and store errors of the file in hand; files already
    /// written stay written, as in the source.
    #[cfg(feature = "mzml")]
    pub fn transform_spectra_files(
        &self,
        ctx: &ToolContext,
        in_spectra_files: &[String],
        out_spectra_files: &[String],
        transformations: &[crate::analysis::transformations::TransformationDescription],
        store_original_rt: bool,
    ) -> Result<()> {
        use crate::analysis::alignment_transformer::MapAlignmentTransformer;
        use crate::format::file_types::FileType;
        use crate::metadata::ProcessingAction;

        if in_spectra_files.is_empty() || out_spectra_files.is_empty() {
            return Ok(());
        }
        if in_spectra_files.len() != transformations.len() {
            return Err(Error::InvalidValue(format!(
                "{} spectra files but {} transformations",
                in_spectra_files.len(),
                transformations.len()
            )));
        }
        if in_spectra_files.len() != out_spectra_files.len() {
            return Err(Error::InvalidValue(format!(
                "{} input but {} output spectra files",
                in_spectra_files.len(),
                out_spectra_files.len()
            )));
        }
        let transformer = MapAlignmentTransformer {
            store_original_rt,
            ..MapAlignmentTransformer::default()
        };
        for ((input, output), transformation) in in_spectra_files
            .iter()
            .zip(out_spectra_files)
            .zip(transformations)
        {
            let mut experiment = FileHandler::load_experiment(input, &[FileType::MzMl])?;
            transformer.transform_experiment(&mut experiment, transformation)?;
            // One record per file, made when the file is done, as the source
            // calls getProcessingInfo_ in the loop: its completion time is
            // that file's.
            let processing = ctx.processing_info(&[ProcessingAction::RetentionTimeAlignment])?;
            ctx.add_data_processing(&mut experiment, &processing);
            FileHandler::store_experiment(output, &experiment, Some(FileType::MzMl))?;
        }
        Ok(())
    }
}
