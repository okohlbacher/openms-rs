// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Detects two-dimensional features in LC-MS data, as the source TOPP tool.
//!
//! Ports `OpenMS4-topp/src/FeatureFinderCentroided.cpp` (topp `174b576`),
//! whose documentation reads: *This module identifies "features" in a LC/MS
//! map. By feature, we understand a peptide in an MS sample that reveals a
//! characteristic isotope distribution over time. The algorithm computes
//! positions in RT and m/z dimension and a charge estimate of each peptide.*
//! The algorithm finds pronounced regions of the data around *seeds*, which the
//! user may provide as a feature map (for example from MS/MS identifications) or
//! which it computes itself, and then fits isotope and retention-time models to
//! them. Centroided MS1 data is needed; PeakPickerHiRes centroids profile data.
//!
//! The source's table of starting values for the `algorithm` section applies
//! unchanged: `intensity:bins` 10 for both Q-TOF and LTQ Orbitrap,
//! `mass_trace:mz_tolerance` 0.02 (Q-TOF) and 0.004 (Orbitrap),
//! `isotopic_pattern:mz_tolerance` 0.04 (Q-TOF) and 0.005 (Orbitrap). The
//! parameters themselves are those of
//! [`default_parameters`](crate::analysis::feature_finder_picked::algorithm::default_parameters).
//! The source cites M. Sturm's dissertation (2010, p. 37 ff.) and Weisser *et
//! al.*, J. Proteome Res. 2013 (PMID 23391308); the TOPP framework does not
//! port tool citations.
//!
//! # Run order
//!
//! The body follows `main_` (`FeatureFinderCentroided.cpp:172-378`) step by
//! step:
//!
//! 1. Load `-in` with MS level 1 and the intensity range `[0, f64::MAX)`
//!    ([`FeatureFinderCentroided::peak_file_options`]), tolerating dangling
//!    mzML header references as the source reader does.
//! 2. No spectrum left:
//!    [`NO_MS1_SPECTRA_MESSAGE`](crate::cli::tools::FeatureFinderCentroided::NO_MS1_SPECTRA_MESSAGE),
//!    exit 4 (`INPUT_FILE_EMPTY`).
//! 3. Any spectrum with per-peak ion mobility:
//!    [`im_peak_message`](crate::cli::tools::FeatureFinderCentroided::im_peak_message),
//!    exit 11 (`INCOMPATIBLE_INPUT_DATA`).
//! 4. A first spectrum stored as profile without `-force`:
//!    [`PROFILE_DATA_MESSAGE`](crate::cli::tools::FeatureFinderCentroided::PROFILE_DATA_MESSAGE),
//!    exit 8 (`UNKNOWN_ERROR`), the code `TOPPBase` gives the source's
//!    `IllegalArgument`.
//! 5. Load `-seeds` when given.
//! 6. Copy the `algorithm:` parameters.
//! 7. Refuse FAIMS input (see below), exit 11, before the algorithm runs.
//! 8. Run the picked feature finder
//!    ([`run_with_options`](crate::analysis::feature_finder_picked::algorithm::run_with_options)).
//! 9. Annotate and clean the features
//!    ([`finish_features`](crate::cli::tools::FeatureFinderCentroided::finish_features)):
//!    the primary MS run path, fresh unique ids, the `Quantitation` processing
//!    record, and hulls reduced to bounding boxes with subordinates removed
//!    below `-debug 5`.
//! 10. Store `-out` as featureXML.
//!
//! No step writes `-out` before step 10, so every refusal and failure leaves no
//! output file.
//!
//! # How far the run matches the executed C++
//!
//! The whole chain runs: `TOPP_FeatureFinderCentroided_1` exits 0 and writes
//! the eight features of the retained expectation, as do the `-seeds`,
//! `-algorithm:feature:rt_shape asymmetric` and `-debug 5` modes. Measured
//! against the C++ **Release** build `bc9cc12`/`174b576` on the same input and
//! INI, the decoded output agrees on every structural field — feature count,
//! charge, hull count, hull point count and order, metadata key sets,
//! `spectra_data` and the single `Quantitation` processing record — with
//! convex-hull coordinates and m/z positions bit-identical and `intensity` and
//! `FWHM` identical as `f32`. What differs is the last bits of the
//! Levenberg-Marquardt fit: at most `5.5e-13` relative on the retention time,
//! `2.2e-10` on `score_fit` and `7.7e-12` on `score_correlation`, against a
//! spread of `2.2e-13`, `9.1e-11` and `3.1e-12` between the C++ Debug and
//! Release builds themselves. `overallquality` is printed by the C++ writer
//! with six decimals, so it can only be compared to that precision, to which it
//! agrees.
//!
//! The algorithm's failures, which the source raises as `IllegalArgument` (for
//! example MS1 spectra that all lose their peaks to the intensity filter:
//! `FeatureFinder needs updated ranges on input map. Aborting.`), are reported
//! as `Error: Unexpected internal error (<message>)` with exit 8, as `TOPPBase`
//! reports them. One such refusal is stricter than the C++ Release build: an
//! input whose MS1 spectra share a single retention time makes the source
//! divide by a zero bin width, which its Debug build catches in a precondition
//! and its Release build carries through to an empty feature map; this port
//! refuses it with exit 8 instead (`tests/topp_feature_finder_centroided.rs`
//! measures both).
//!
//! # FAIMS input is refused
//!
//! The source splits FAIMS data by compensation voltage
//! (`IMDataConverter::splitByFAIMSCV`), runs the algorithm once per voltage,
//! annotates each feature with its voltage and, with `-faims_merge_features
//! true`, merges features across voltages (`FeatureOverlapFilter`). Decision D5
//! of the early TOPP bundle defers that closure: the executed C++ tool exits 8
//! on every FAIMS input (its voltage groups have no updated ranges), and its
//! merge would remove every feature. This port therefore refuses any input in
//! which
//! [`FaimsHelper::get_compensation_voltages`](crate::kernel::faims_helper::FaimsHelper::get_compensation_voltages)
//! finds a voltage, with
//! [`faims_refusal_message`](crate::cli::tools::FeatureFinderCentroided::faims_refusal_message),
//! exit 11 and no output, instead of pooling the voltages silently. The refusal
//! comes at the source's position of the split, after the profile check and the
//! seed load, so a FAIMS profile file without `-force` still exits 8 with the
//! profile message, as in the C++ tool. `-faims_merge_features` is registered
//! and validated but has no other effect yet.
//!
//! # Threads
//!
//! `-threads` is accepted and validated by the TOPP framework and reaches the
//! algorithm as
//! [`Options::threads`](crate::analysis::feature_finder_picked::algorithm::Options::threads),
//! which sizes the seed-extension loop — the port's form of the source's single
//! `#pragma omp parallel for` (`FeatureFinderAlgorithmPicked.cpp:595`). By the
//! determinism contract the result does not depend on the count: the loop
//! returns its results in seed order and every later step is serial, so
//! `-threads 1`, `2`, `4`, `8` and `0` write byte-identical output, unique ids
//! included.
//!
//! See `docs/TOPP_FEATURE_FINDER_CENTROIDED_SUPPORT.md` for the API mapping, the
//! preserved source conventions, the native differences and the evidence.

use crate::analysis::feature_finder_picked::algorithm::{self, Options};
use crate::cli::{ExitCode, Tool, ToolContext, ToolSpec};
use crate::concept::HasUniqueId;
use crate::format::file_handler::FileHandler;
use crate::format::file_types::FileType;
use crate::format::mzml::ReadOptions;
use crate::format::peak_options::PeakFileOptions;
use crate::kernel::faims_helper::FaimsHelper;
use crate::kernel::{FeatureMap, MSExperiment, NumericRange, SpectrumType};
use crate::metadata::{
    ImTypes, IonMobilityFormat, IonMobilityPeakType, ProcessingAction, im_peak_type_to_string,
};
use crate::param::Param;
use crate::system::file;
use crate::{Error, Result};
use std::io::Write;

/// The `FeatureFinderCentroided` TOPP tool.
///
/// Registration, `-write_ini`, every wrapper branch and the algorithm itself
/// run, so a non-FAIMS input produces a feature map; FAIMS input is refused
/// (see the module documentation).
pub struct FeatureFinderCentroided;

impl FeatureFinderCentroided {
    /// The empty-input diagnostic (`FeatureFinderCentroided.cpp:196`): the
    /// source's `FileEmpty` message as `TOPPBase` prints that exception
    /// (`TOPPBase.cpp:442-447`), verbatim from the oracle case `FFC_ms2_only`.
    pub const NO_MS1_SPECTRA_MESSAGE: &'static str =
        "Error: File empty (the file 'Error: No MS1 spectra in input file.' is empty)";

    /// The profile-data diagnostic (`FeatureFinderCentroided.cpp:220`) inside
    /// `TOPPBase`'s `UNKNOWN_ERROR` wording (`TOPPBase.cpp:495-499`), verbatim
    /// from the oracle case `FFC_profile_noforce`.
    pub const PROFILE_DATA_MESSAGE: &'static str = "Error: Unexpected internal error (Error: Profile data provided but centroided spectra expected. To enforce processing of the data set the -force flag.)";

    /// The informational line `IMDataConverter::splitByFAIMSCV` logs for data
    /// without FAIMS voltages (`IMDataConverter.cpp:34`), which the C++ tool
    /// prints on every run that reaches the split; verbatim, including the
    /// source's wording.
    pub const NO_FAIMS_MESSAGE: &'static str =
        "Not FAIMS compensation voltages found in the data. Returning PeakMap as CV NaN.";

    /// The per-peak ion-mobility diagnostic
    /// (`FeatureFinderCentroided.cpp:205-208`).
    ///
    /// The source prints `imPeakTypeToString(spec.getIMPeakType())` in the
    /// parentheses. Its mzML reader sets `IM_PROFILE` on every spectrum with an
    /// ion-mobility array unless the spectrum carries `MS:1003441` (ion
    /// mobility centroid frame), which sets `IM_CENTROIDED`
    /// (`MzMLHandler.cpp:248-256` and `1646-1649`). The native spectrum has no
    /// ion-mobility peak type and the native reader does not keep
    /// `MS:1003441`, so this port always prints `im_profile`, the text both
    /// oracle cases `FFC_im_peak_with_units` and `FFC_im_peak_without_units`
    /// printed.
    pub fn im_peak_message() -> String {
        format!(
            "Error: Input contains per-peak ion mobility data (IM_PEAK, {}) which is not supported by FeatureFinderCentroided. Preprocess with IonMobilityBinning or PeakPickerIM first.",
            im_peak_type_to_string(IonMobilityPeakType::Profile)
        )
    }

    /// The source's informational line for FAIMS input
    /// (`FeatureFinderCentroided.cpp:243`), for `count` compensation voltages.
    pub fn faims_detected_message(count: usize) -> String {
        format!("FAIMS data detected with {count} compensation voltage(s).")
    }

    /// The refusal of FAIMS input, naming the voltages found, ascending, in
    /// volts.
    ///
    /// Native: the source processes FAIMS input, and its executed C++ build
    /// fails on every such input with exit 8; see the module documentation and
    /// decision D5.
    pub fn faims_refusal_message(voltages: &[f64]) -> String {
        let list: Vec<String> = voltages.iter().map(|volts| format!("{volts}")).collect();
        format!(
            "Error: FAIMS input is not supported by this port of FeatureFinderCentroided yet (compensation voltages {} V): per-voltage feature detection and the cross-voltage merge are not ported. No output was written.",
            list.join(", ")
        )
    }

    /// The source `PeakFileOptions` of `main_`
    /// (`FeatureFinderCentroided.cpp:179-184`): MS level 1 and the intensity
    /// range `[0, f64::MAX)`.
    ///
    /// The source passes `std::numeric_limits<DRange<1>::PositionType>::min()`
    /// as the lower bound, under the comment *filter out zero (and negative)
    /// intensities*. `DPosition<1>` has no `numeric_limits` specialisation, so
    /// the call returns `DPosition()`, zero, and the executed range keeps zero
    /// intensities; the upper bound `DPosition<1>::maxPositive()` is
    /// `f64::MAX`. The loader's ranges are half-open like `DRange::encloses`, so
    /// negative intensities are dropped. This follows the executed C++ (the A3
    /// oracle keeps `0.0`, a subnormal and `DBL_MIN` and drops `-1`), not the
    /// comment.
    ///
    /// # Errors
    ///
    /// Propagates the option container's MS-level limit error, which one level
    /// cannot reach.
    pub fn peak_file_options() -> Result<PeakFileOptions> {
        let mut options = PeakFileOptions::default();
        options.add_ms_level(1)?;
        options.set_intensity_range(NumericRange {
            min: 0.0,
            max: f64::MAX,
        });
        Ok(options)
    }

    /// Annotate and clean the features the algorithm returned, as `main_` does
    /// between the algorithm and the store
    /// (`FeatureFinderCentroided.cpp:318-373`).
    ///
    /// In source order:
    ///
    /// 1. The primary MS run path (the `spectra_data` meta value) becomes
    ///    `file://<basename of input>` under `-test`, so test outputs do not
    ///    depend on the directory, and `input` as given otherwise.
    /// 2. The map receives a unique id if it has none, and then the map, every
    ///    feature and every subordinate feature receive new ones, depth first,
    ///    all from one generator of
    ///    [`ToolContext::unique_id_generator`](crate::cli::ToolContext::unique_id_generator)
    ///    (seeded under `-test`), as `ensureUniqueId` followed by
    ///    `applyMemberFunction(&UniqueIdInterface::setUniqueId)`. A map without
    ///    an id therefore consumes one draw that is overwritten at once.
    /// 3. Above `-debug 10`, every feature with metadata is listed on `out` as
    ///    `Feature <unique id>` followed by one `  <key> = <value>` line per
    ///    entry. The native metadata is ordered by key and prints numbers in
    ///    Rust's shortest round-trip form, while the source lists keys in
    ///    meta-registry order and formats numbers with `StringUtils::toStr`.
    /// 4. The processing record of this run with the `Quantitation` action is
    ///    appended
    ///    ([`ToolContext::processing_info`](crate::cli::ToolContext::processing_info)).
    /// 5. Below `-debug 5`, every mass-trace hull of every feature becomes its
    ///    bounding box and the subordinate features are removed, *to reduce file
    ///    size of feature files*. From level 5 on, hulls and subordinates are
    ///    kept.
    ///
    /// The source also expands each feature's overall hull in step 5, which
    /// only mutates a cache: the native overall hull is computed on demand, so
    /// there is nothing to expand. An empty mass-trace hull stays empty here;
    /// the source would replace it with the corners of an empty bounding box,
    /// but the algorithm never produces an empty hull.
    ///
    /// The map is taken by value and returned, so a failure leaves no
    /// half-annotated map behind.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a floating-point parameter cannot be
    /// recorded in the processing record, and the limit errors of
    /// [`FeatureMap::set_primary_ms_run_path`] and
    /// [`FeatureMap::for_each_unique_id`]. A failing `out` stream is an I/O
    /// error.
    pub fn finish_features(
        ctx: &ToolContext,
        input: &str,
        mut features: FeatureMap,
        out: &mut dyn Write,
    ) -> Result<FeatureMap> {
        let run_path = if ctx.test_mode() {
            format!("file://{}", file::basename(input))
        } else {
            input.to_owned()
        };
        features.set_primary_ms_run_path(&[run_path])?;

        let mut generator = ctx.unique_id_generator();
        features.unique_id.ensure_unique_id(&mut generator);
        features.for_each_unique_id(|id| {
            *id = generator.get_unique_id();
            1
        })?;

        if ctx.debug_level() > 10 {
            for feature in &features.features {
                if feature.metadata.is_empty() {
                    continue;
                }
                writeln!(out, "Feature {}", feature.unique_id)?;
                for (key, value) in &feature.metadata {
                    writeln!(out, "  {key} = {value}")?;
                }
            }
        }

        let processing = ctx.processing_info(&[ProcessingAction::Quantitation])?;
        ctx.add_data_processing(&mut features, &processing);

        if ctx.debug_level() < 5 {
            for feature in &mut features.features {
                for hull in &mut feature.convex_hulls {
                    hull.expand_to_bounding_box();
                }
                feature.subordinates.clear();
            }
        }
        Ok(features)
    }
}

/// Whether any spectrum carries per-peak ion mobility, as the loop at
/// `FeatureFinderCentroided.cpp:200-211` asks `IMTypes::determineIMFormat` of
/// each one.
///
/// The native spectrum stores no ion-mobility format, and the source mzML
/// reader never sets one, so the stored-format short circuit of the source
/// function always falls through to the data, as
/// [`ImTypes::determine_im_format`] does.
fn has_per_peak_mobility(experiment: &MSExperiment) -> bool {
    experiment
        .spectra
        .iter()
        .any(|spectrum| ImTypes::determine_im_format(spectrum) == IonMobilityFormat::PerPeak)
}

impl Tool for FeatureFinderCentroided {
    const NAME: &'static str = "FeatureFinderCentroided";
    const DESCRIPTION: &'static str = "Detects two-dimensional features in LC-MS data.";

    /// Source `registerOptionsAndFlags_` (`FeatureFinderCentroided.cpp:140-163`).
    ///
    /// `-in` accepts mzML only: the source adds Thermo `raw` when it is built
    /// `WITH_THERMO_RAW`, and this port has no reader for it.
    fn register(spec: &mut ToolSpec) -> Result<()> {
        spec.register_input_file("in", "<file>", "", "input file", true, false, &[])?;
        spec.set_valid_formats("in", &["mzML"])?;
        spec.register_output_file("out", "<file>", "", "output file", true, false)?;
        spec.set_valid_formats("out", &["featureXML"])?;
        spec.register_input_file(
            "seeds",
            "<file>",
            "",
            "User specified seed list",
            false,
            false,
            &[],
        )?;
        spec.set_valid_formats("seeds", &["featureXML"])?;
        spec.add_empty_line()?;
        spec.register_string_option(
            "faims_merge_features",
            "<true/false>",
            "true",
            "For FAIMS data with multiple compensation voltages: Merge features representing the same analyte detected at different CV values into a single feature. Only features with DIFFERENT FAIMS CV values are merged (same CV = different analytes). Has no effect on non-FAIMS data.",
            false,
            false,
        )?;
        spec.set_valid_strings("faims_merge_features", &["true", "false"])?;
        spec.add_empty_line()?;
        spec.register_subsection("algorithm", "Algorithm section")?;
        Ok(())
    }

    /// Source `getSubsectionDefaults_`, which returns
    /// `FeatureFinderAlgorithmPicked().getDefaultParameters()` whatever the
    /// name; `algorithm` is the only registered subsection.
    fn subsection_defaults(_section: &str) -> Result<Option<Param>> {
        algorithm::default_parameters().map(Some)
    }

    /// Forwards to [`run_io`](Tool::run_io) with the process streams.
    fn run(ctx: &ToolContext) -> Result<ExitCode> {
        Self::run_io(ctx, &mut std::io::stdout(), &mut std::io::stderr())
    }

    /// Source `main_`; the module documentation lists the steps, messages and
    /// exit codes. Informational lines go to `out`, diagnostics to `err`.
    fn run_io(ctx: &ToolContext, out: &mut dyn Write, err: &mut dyn Write) -> Result<ExitCode> {
        let input = ctx.string("in")?.to_owned();
        let output = ctx.string("out")?.to_owned();

        // Load MS1 spectra with the executed intensity range (179-192). The
        // source FileHandler tolerates dangling mzML header references; the
        // native reader does only on request (decision D10).
        let read = ReadOptions {
            source_dangling_references: true,
            ..ReadOptions::default()
        };
        let experiment = FileHandler::load_experiment_with_read_options(
            &input,
            &[FileType::MzMl, FileType::Raw],
            &Self::peak_file_options()?,
            &read,
        )?;

        // FileEmpty, which TOPPBase maps to INPUT_FILE_EMPTY (194-197).
        if experiment.spectra.is_empty() {
            writeln!(err, "{}", Self::NO_MS1_SPECTRA_MESSAGE)?;
            return Ok(ExitCode::InputFileEmpty);
        }

        // Per-peak ion mobility (200-211).
        if has_per_peak_mobility(&experiment) {
            writeln!(err, "{}", Self::im_peak_message())?;
            return Ok(ExitCode::IncompatibleInputData);
        }

        // Only the first spectrum's stored type is checked (214-222). The
        // IllegalArgument reaches TOPPBase's catch-all, UNKNOWN_ERROR.
        if experiment.spectra[0].spectrum_type == SpectrumType::Profile && !ctx.force() {
            writeln!(err, "{}", Self::PROFILE_DATA_MESSAGE)?;
            return Ok(ExitCode::UnknownError);
        }

        // Seeds (224-229).
        let seeds_path = ctx.string("seeds")?;
        let seeds = if seeds_path.is_empty() {
            FeatureMap::new()
        } else {
            FileHandler::load_feature_map(seeds_path, &[FileType::FeatureXml])?
        };

        // Parameters of the feature finder (232). The source's dump of them at
        // debug level 3 (writeDebug_) is not ported by the TOPP framework.
        let parameters = ctx.subsection("algorithm")?;

        // The FAIMS split (238-244), refused until its closure is ported (D5).
        let voltages = FaimsHelper::get_compensation_voltages(&experiment)?;
        for warning in &voltages.warnings {
            writeln!(err, "{warning}")?;
        }
        if !voltages.voltages.is_empty() {
            let values: Vec<f64> = voltages.values().collect();
            writeln!(out, "{}", Self::faims_detected_message(values.len()))?;
            writeln!(err, "{}", Self::faims_refusal_message(&values))?;
            return Ok(ExitCode::IncompatibleInputData);
        }
        writeln!(out, "{}", Self::NO_FAIMS_MESSAGE)?;

        // The algorithm (283-291), on the worker count -threads asks for, as
        // TOPPBase applies the setting before main_ (TOPPBase.cpp:408-415).
        let options = Options {
            threads: ctx.thread_policy(),
            ..Options::default()
        };
        let outcome = algorithm::run_with_options(experiment, &seeds, &parameters, &options);
        let result = match outcome {
            Ok(result) => result,
            Err(Error::InvalidValue(message)) => {
                // The algorithm's IllegalArgument and InvalidValue exceptions
                // reach TOPPBase's catch-all (TOPPBase.cpp:495-499).
                writeln!(err, "Error: Unexpected internal error ({message})")?;
                return Ok(ExitCode::UnknownError);
            }
            Err(error) => return Err(error),
        };
        for line in &result.log {
            writeln!(out, "{line}")?;
        }

        // Annotation and clean-up (318-373), then the store (375).
        let features = Self::finish_features(ctx, &input, result.features, out)?;
        FileHandler::store_feature_map(&output, &features, Some(FileType::FeatureXml))?;
        Ok(ExitCode::ExecutionOk)
    }
}
