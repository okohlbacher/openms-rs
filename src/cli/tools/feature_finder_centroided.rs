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
//! 7. Split the input by FAIMS compensation voltage
//!    ([`ImDataConverter::split_by_faims_cv`](crate::kernel::im_data_converter::ImDataConverter::split_by_faims_cv)),
//!    which returns one group per voltage, or a single group holding the whole
//!    input when there is none.
//! 8. Run the picked feature finder
//!    ([`run_with_options`](crate::analysis::feature_finder_picked::algorithm::run_with_options))
//!    once per group, on that group's seeds, and annotate the group's features
//!    with their voltage; with `-faims_merge_features true` merge features of
//!    the same analyte across voltages (see below).
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
//! `-algorithm:feature:rt_shape asymmetric` and `-debug 5` modes, and the
//! decoded output agrees with the C++ **Release** build `bc9cc12`/`174b576` on
//! every structural field — feature count, charge, hull count, hull point
//! count and order, metadata key sets, `spectra_data` and the single
//! `Quantitation` processing record. The algorithm underneath is pinned bit
//! for bit against that build on every platform (its fits call that build's
//! glibc `exp` and `log`, ported; `docs/FEATURE_FINDER_PICKED_SUPPORT.md`); the tool's
//! own last-bits comparison predates that and is recorded in
//! `docs/TOPP_FEATURE_FINDER_CENTROIDED_SUPPORT.md`. `overallquality` is
//! printed by the C++ writer with six decimals, so it can only be compared to
//! that precision.
//!
//! The algorithm's failures, which the source raises as `IllegalArgument` (for
//! example MS1 spectra that all lose their peaks to the intensity filter:
//! `FeatureFinder needs updated ranges on input map. Aborting.`) or as another
//! OpenMS exception, are reported as `Error: Unexpected internal error
//! (<message>)` with exit 8, as `TOPPBase` reports them; so are the
//! algorithm's refusals where the source crashes (an empty best isotope
//! pattern with `feature:min_isotope_fit` 0, SIGSEGV in the executed tool). The `std::length_error`
//! of step 2.5, which is no OpenMS exception, reaches `TOPPBase`'s outer
//! `std::exception` handler instead (`TOPPBase.cpp:510-514`): `Unable to
//! initialize or run FeatureFinderCentroided: vector::_M_default_append`, exit
//! 12 (`INTERNAL_ERROR`), after the debug directory and the first log line of a
//! `-algorithm:write_debug` run (executed: `../oracle/ffap-complete-fix2`,
//! `tool_1e19`). The port's own isotope-window ceiling below that bound exits 8
//! with its message where the executed tool exits 12 with `std::bad_alloc`
//! (`tool_2e18`). An input whose MS1 spectra share a single retention time or
//! m/z follows the Release build to an empty feature map
//! (`tests/topp_feature_finder_centroided.rs`).
//!
//! # FAIMS input, and the two defects the closure corrects
//!
//! The source splits FAIMS data by compensation voltage
//! (`IMDataConverter::splitByFAIMSCV`), runs the algorithm once per voltage,
//! annotates each feature with its voltage and, with `-faims_merge_features
//! true`, merges features of the same analyte across voltages
//! (`FeatureOverlapFilter::mergeFAIMSFeatures`). All of that is ported here.
//! The executed C++ tool nevertheless fails on **every** FAIMS input, and two
//! of its defects lie on this path; the tool ships the corrected behaviour and
//! names each point:
//!
//! - **`CPP-278`, the missing ranges, cannot arise here.** The source builds
//!   each voltage group with `addSpectrum` and never calls `updateRanges`, so
//!   `FeatureFinderAlgorithmPicked` throws `the value '1' was used but is not
//!   valid; No ranges for this MS level` on the first group and every FAIMS
//!   input exits 8 (re-executed: `../oracle/b11-faims`, eight runs on six
//!   FAIMS inputs, all rc 8).
//!   The native containers compute ranges on demand
//!   ([`MSExperiment::spectrum_range_manager`]), so a group has its own ranges
//!   the moment it holds spectra; there is no state to forget and nothing to
//!   reproduce. This is a property of the container port, not a choice made
//!   here.
//! - **`CPP-282`, the merge that erases everything.** `mergeFAIMSFeatures`
//!   records removal by unique ID, and the algorithm returns every feature with
//!   ID 0, so one merge marks ID 0 removed and the final `erase` drops every
//!   feature. The tool draws a unique ID per feature before merging, so the
//!   merge keys on real IDs; those IDs are overwritten immediately afterwards
//!   by the source's own `applyMemberFunction(setUniqueId)`, so nothing else
//!   changes.
//! - **`CPP-283`, the double count and the survivor that stops absorbing.**
//!   [`FAIMS_MERGE_FIDELITY`](crate::cli::tools::FeatureFinderCentroided::FAIMS_MERGE_FIDELITY)
//!   is [`FaimsMergeFidelity::Corrected`], which skips a candidate already
//!   removed and lets a survivor keep absorbing voltages it does not yet stand
//!   for. A cluster of features at pairwise different voltages therefore
//!   collapses to one feature whose intensity is the sum of the cluster, each
//!   member counted once.
//!
//! Two further defects of the split are handled by the ported library: a NaN
//! compensation voltage is refused rather than allowed to break an ordered set
//! (`CPP-280`), and the chromatograms the source destroys are returned instead
//! of dropped (`CPP-279`) — this tool loads MS level 1 only and uses none, so
//! it discards them as the source does. The information line of a non-FAIMS
//! input keeps the source's wording, `Not FAIMS compensation voltages …`
//! (`CPP-281`), because it is the line every executed run prints.
//!
//! The split comes at the source's position, after the profile check and the
//! seed load, so a FAIMS profile file without `-force` still exits 8 with the
//! profile message and an unreadable `-seeds` file still exits 3, as in the C++
//! tool.
//!
//! There is no whole-tool C++ oracle for the corrected path, because the C++
//! path is broken. Each voltage group is pinned instead against the Release
//! build run on that group written as its own single-voltage mzML, and the
//! merge against the specification derived in
//! [`FaimsMergeFidelity`] and hand-derived cases; see
//! `docs/TOPP_FEATURE_FINDER_CENTROIDED_SUPPORT.md`.
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
use crate::analysis::feature_finder_picked::debug::{DebugOutput, ReportLine, TerminationKind};
use crate::analysis::feature_finder_picked::instance::FeatureFinderAlgorithmPicked;
use crate::analysis::feature_finder_picked::seeds;
use crate::cli::{ExitCode, Tool, ToolContext, ToolError, ToolResult, ToolSpec};
use crate::concept::constants::user_param::FAIMS_CV;
use crate::concept::{HasUniqueId, UniqueIdGenerator};
use crate::format::file_handler::FileHandler;
use crate::format::file_info::text_format::{DEFAULT_STREAM_PRECISION, ostream_g};
use crate::format::file_types::FileType;
use crate::format::mzml::ReadOptions;
use crate::format::peak_options::PeakFileOptions;
use crate::kernel::faims_helper::FaimsHelper;
use crate::kernel::im_data_converter::{FaimsSplitLogLevel, ImDataConverter};
use crate::kernel::{FeatureMap, MSExperiment, NumericRange, SpectrumType};
use crate::metadata::{
    ImTypes, IonMobilityFormat, IonMobilityPeakType, MetaValue, ProcessingAction,
    im_peak_type_to_string,
};
use crate::param::Param;
use crate::processing::feature_overlap_filter::{FaimsMergeFidelity, FeatureOverlapFilter};
use crate::system::file;
use crate::{Error, Result};
use std::io::Write;
use std::path::Path;

/// The `FeatureFinderCentroided` TOPP tool.
///
/// Registration, `-write_ini`, every wrapper branch and the algorithm itself
/// run, and so does the FAIMS closure: the split by compensation voltage, one
/// algorithm run per voltage, the `FAIMS_CV` annotation and the cross-voltage
/// merge, with the two defects of the source merge corrected (see the module
/// documentation).
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
    /// source's wording (`CPP-281`, "Not" for "No").
    ///
    /// The split itself produces it, as
    /// [`ImDataConverter::NO_COMPENSATION_VOLTAGES_INFO`](crate::kernel::im_data_converter::ImDataConverter::NO_COMPENSATION_VOLTAGES_INFO);
    /// this is the same text under the tool's name.
    pub const NO_FAIMS_MESSAGE: &'static str = ImDataConverter::NO_COMPENSATION_VOLTAGES_INFO;

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

    /// The source's line before each compensation-voltage group
    /// (`FeatureFinderCentroided.cpp:254`), for `volts` and a group of
    /// `spectra` spectra.
    ///
    /// The voltage is formatted as `std::ostream << double` formats it at the
    /// log stream's default precision 6, the text the executed runs printed
    /// (`Processing FAIMS CV group: -60 V (56 spectra)`).
    pub fn processing_group_message(volts: f64, spectra: usize) -> String {
        format!(
            "Processing FAIMS CV group: {} V ({spectra} spectra)",
            ostream_g(volts, DEFAULT_STREAM_PRECISION)
        )
    }

    /// The source's line after the group loop
    /// (`FeatureFinderCentroided.cpp:306`).
    pub fn combined_features_message(features: usize) -> String {
        format!("Combined {features} features from all FAIMS CV groups.")
    }

    /// The source's line after the cross-voltage merge
    /// (`FeatureFinderCentroided.cpp:313-314`).
    ///
    /// The source computes the third number as `before - after` in `Size`; the
    /// merge only removes features, so it never wraps.
    pub fn faims_merge_message(before: usize, after: usize) -> String {
        format!(
            "FAIMS feature merge: {before} -> {after} features (merged {})",
            before.saturating_sub(after)
        )
    }

    /// The largest retention-time difference of the cross-voltage merge, in
    /// seconds: the source's literal argument
    /// (`FeatureFinderCentroided.cpp:312`), not the parameter default.
    pub const FAIMS_MERGE_MAX_RT_DIFF: f64 = 5.0;

    /// The largest m/z difference of the cross-voltage merge, in Da
    /// (`FeatureFinderCentroided.cpp:312`).
    pub const FAIMS_MERGE_MAX_MZ_DIFF: f64 = 0.05;

    /// Which cross-voltage merge the tool runs: the **corrected** one, the
    /// tool's one designed difference in the FAIMS closure.
    ///
    /// [`FaimsMergeFidelity::Corrected`] answers two executed defects of the
    /// source merge:
    ///
    /// - `CPP-282`: removal keys on unique IDs, and
    ///   `FeatureFinderAlgorithmPicked` returns every feature with ID 0, so the
    ///   source's merge erases every FAIMS feature. The tool draws a unique ID
    ///   for each feature before merging, from the generator whose later draws
    ///   overwrite them, so the merge keys on real IDs.
    /// - `CPP-283`: the loop offers a feature it has already removed to a later
    ///   survivor, and a survivor refuses every merge after its first. A
    ///   cluster of features at pairwise different voltages therefore collapses
    ///   to one feature whose intensity is the sum of the cluster, each member
    ///   counted once.
    ///
    /// Setting it to [`FaimsMergeFidelity::Source`] would reproduce the source,
    /// including the unique IDs it leaves at 0 — the ID draw above is part of
    /// the corrected route. On any FAIMS input where at least one cross-voltage
    /// merge actually fires, so from two voltages with an overlapping cluster
    /// upwards, that writes an **empty** feature map. It does not on input
    /// where nothing merges: with `-faims_merge_features false` the merge never
    /// runs, and a single voltage merges nothing, because the source's callback
    /// refuses a candidate whose voltage equals the survivor's, so no ID is
    /// ever marked removed. The empty map is a silent wrong answer, not a
    /// crash, so it is reproducible; it is exercised at the library level
    /// (`tests/feature_overlap_filter.rs`) rather than shipped in the tool.
    pub const FAIMS_MERGE_FIDELITY: FaimsMergeFidelity = FaimsMergeFidelity::Corrected;

    /// The seeds of one compensation-voltage group
    /// (`FeatureFinderCentroided.cpp:258-281`).
    ///
    /// Without FAIMS input, or with an empty seed list, the whole list is used
    /// as it is. Otherwise a seed with a `FAIMS_CV` meta value joins the group
    /// when the two voltages differ by less than
    /// [`FaimsHelper::DEFAULT_CV_TOLERANCE`](crate::kernel::faims_helper::FaimsHelper::DEFAULT_CV_TOLERANCE)
    /// (the source's literal `0.01`), and a seed without one joins every group,
    /// *for backward compatibility* as the source comment says. The filtered
    /// list is a fresh map, as the source's `FeatureMap seeds_cv` is, so it
    /// carries none of the seed map's own identifier, metadata or
    /// identifications.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a seed's `FAIMS_CV` is not numeric; the
    /// source's `double seed_cv = seed.getMetaValue(...)` throws
    /// `ConversionError` there.
    pub fn seeds_of_group(seeds: &FeatureMap, has_faims: bool, volts: f64) -> Result<FeatureMap> {
        if !has_faims || seeds.features.is_empty() {
            return Ok(seeds.clone());
        }
        let mut group = FeatureMap::new();
        for seed in &seeds.features {
            let keep = match seed.metadata.get(FAIMS_CV) {
                Some(value) => (value.as_f64()? - volts).abs() < FaimsHelper::DEFAULT_CV_TOLERANCE,
                None => true,
            };
            if keep {
                group.features.push(seed.clone());
            }
        }
        Ok(group)
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
        features: FeatureMap,
        out: &mut dyn Write,
    ) -> Result<FeatureMap> {
        let mut generator = ctx.unique_id_generator();
        Self::finish_features_with(ctx, input, features, out, &mut generator)
    }

    /// [`Self::finish_features`] drawing the unique ids from `generator`.
    ///
    /// The source draws every unique id from one process-wide generator, so an
    /// id drawn earlier in the run (the debug abort map's) shifts the ids of
    /// the output; the tool passes the generator it drew that id from.
    ///
    /// # Errors
    ///
    /// As [`Self::finish_features`].
    pub fn finish_features_with(
        ctx: &ToolContext,
        input: &str,
        mut features: FeatureMap,
        out: &mut dyn Write,
        generator: &mut UniqueIdGenerator,
    ) -> Result<FeatureMap> {
        let run_path = if ctx.test_mode() {
            format!("file://{}", file::basename(input))
        } else {
            input.to_owned()
        };
        features.set_primary_ms_run_path(&[run_path])?;

        features.unique_id.ensure_unique_id(generator);
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
    fn run(ctx: &ToolContext) -> ToolResult {
        Self::run_io(ctx, &mut std::io::stdout(), &mut std::io::stderr())
    }

    /// Source `main_`; the module documentation lists the steps, messages and
    /// exit codes. Informational lines go to `out`, diagnostics to `err`.
    ///
    /// Only the per-peak ion mobility refusal is an exit code the source's
    /// `main_` returns, after an `OPENMS_LOG_ERROR` line, so only it is
    /// followed by the lifecycle's closing `FeatureFinderCentroided took … .`
    /// line. Every other refusal stands for an exception the source throws
    /// and `TOPPBase::main` catches, and is a [`ToolError`]: no closing line,
    /// and the catch block's text reaches the `-log` file.
    fn run_io(ctx: &ToolContext, out: &mut dyn Write, err: &mut dyn Write) -> ToolResult {
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

        // The thrown FileEmpty, which TOPPBase maps to INPUT_FILE_EMPTY
        // (194-197).
        if experiment.spectra.is_empty() {
            return Err(ToolError::caught(
                ExitCode::InputFileEmpty,
                Self::NO_MS1_SPECTRA_MESSAGE,
            ));
        }

        // Per-peak ion mobility (200-211): an OPENMS_LOG_ERROR line, which
        // does not reach the log file, and an exit code main_ returns.
        if has_per_peak_mobility(&experiment) {
            writeln!(err, "{}", Self::im_peak_message())?;
            return Ok(ExitCode::IncompatibleInputData);
        }

        // Only the first spectrum's stored type is checked (214-222). The
        // thrown IllegalArgument reaches TOPPBase's catch-all, UNKNOWN_ERROR.
        if experiment.spectra[0].spectrum_type == SpectrumType::Profile && !ctx.force() {
            return Err(ToolError::caught(
                ExitCode::UnknownError,
                Self::PROFILE_DATA_MESSAGE,
            ));
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

        // The two log-stream caches of the process: OPENMS_LOG_INFO on stdout
        // and OPENMS_LOG_WARN on stderr.
        let mut info = LogStreamLines::default();
        let mut warn = LogStreamLines::default();
        let mut generator = ctx.unique_id_generator();

        // The FAIMS split (238-244). Its log records come first, as the source
        // writes them from inside splitByFAIMSCV.
        let mut experiment = experiment;
        let split = ImDataConverter::split_by_faims_cv(&mut experiment)?;
        for message in &split.messages {
            match message.level() {
                FaimsSplitLogLevel::Info => info.line(out, &message.text())?,
                FaimsSplitLogLevel::Warning => warn.line(err, &message.text())?,
            }
        }
        let has_faims = split.has_faims();
        if has_faims {
            info.line(out, &Self::faims_detected_message(split.groups.len()))?;
        }

        // One run of the algorithm per compensation-voltage group (246-302),
        // in ascending voltage order; exactly one group for non-FAIMS input.
        let mut features = FeatureMap::new();
        for group in split.groups {
            let volts = group.key.volts();
            if has_faims {
                info.line(
                    out,
                    &Self::processing_group_message(volts, group.experiment.spectra.len()),
                )?;
            }
            // A non-numeric FAIMS_CV on a seed is a ConversionError in the
            // source, which reaches TOPPBase's catch-all.
            let group_seeds =
                Self::seeds_of_group(&seeds, has_faims, volts).map_err(ToolError::unexpected)?;
            let mut group_features = FeatureMap::new();
            run_group(
                ctx,
                group.experiment,
                &mut group_features,
                &parameters,
                &group_seeds,
                out,
                err,
                &mut info,
                &mut warn,
                &mut generator,
            )?;
            if features.features.len() + group_features.features.len() > FeatureMap::MAX_ITEMS {
                return Err(Error::InvalidValue(format!(
                    "the FAIMS groups produced more than {} features",
                    FeatureMap::MAX_ITEMS
                ))
                .into());
            }
            for mut feature in group_features.features {
                if has_faims {
                    feature
                        .metadata
                        .insert(FAIMS_CV.to_owned(), MetaValue::try_from(volts)?);
                }
                features.features.push(feature);
            }
        }

        // The cross-voltage merge (304-318).
        if has_faims {
            info.line(
                out,
                &Self::combined_features_message(features.features.len()),
            )?;
            if ctx.string("faims_merge_features")? == "true" {
                let before = features.features.len();
                if Self::FAIMS_MERGE_FIDELITY == FaimsMergeFidelity::Corrected {
                    // The merge keys removal on unique ids and the algorithm
                    // leaves every feature with id 0, so the source erases them
                    // all (CPP-282). The ids drawn here are overwritten by
                    // `finish_features_with` a few lines below, exactly as the
                    // source's `applyMemberFunction(setUniqueId)` overwrites
                    // whatever the features carried.
                    for feature in &mut features.features {
                        feature.unique_id = generator.get_unique_id();
                    }
                }
                let merged = FeatureOverlapFilter::merge_faims_features_with_fidelity(
                    &mut features,
                    Self::FAIMS_MERGE_MAX_RT_DIFF,
                    Self::FAIMS_MERGE_MAX_MZ_DIFF,
                    Self::FAIMS_MERGE_FIDELITY,
                );
                merged.map_err(ToolError::unexpected)?;
                info.line(
                    out,
                    &Self::faims_merge_message(before, features.features.len()),
                )?;
            }
        }

        // Annotation and clean-up (318-373), then the store (375).
        let features = Self::finish_features_with(ctx, &input, features, out, &mut generator)?;
        // A store that fails is the source's `UnableToCreateFile`
        // (`TOPPBase.cpp:430-435`), not a read failure; see
        // [`crate::cli::write_failure`].
        FileHandler::store_feature_map(&output, &features, Some(FileType::FeatureXml))
            .map_err(|error| crate::cli::write_failure(&output, &error))?;
        // The log streams' caches at exit; the lifecycle prints TOPPBase's
        // closing info line after this returns.
        info.close(out)?;
        warn.close(err)?;
        Ok(ExitCode::ExecutionOk)
    }
}

/// One run of the picked feature finder on one compensation-voltage group
/// (`FeatureFinderCentroided.cpp:283-291`), with the debug files and console
/// lines that run writes.
///
/// A fresh [`FeatureFinderAlgorithmPicked`] per group, as the source's loop
/// body creates one, running into the empty `features` of that group. Returns
/// the [`ToolError`] of the source exception that ended the run, when it
/// failed, and the tool must stop; the debug files the run had written are on
/// disk either way, as they are in the executed process. The
/// source's fixed debug file names mean a later group overwrites an earlier
/// group's files, which this reproduces by writing them the same way per group.
#[allow(clippy::too_many_arguments)]
fn run_group(
    ctx: &ToolContext,
    experiment: MSExperiment,
    features: &mut FeatureMap,
    parameters: &Param,
    seeds: &FeatureMap,
    out: &mut dyn Write,
    err: &mut dyn Write,
    info: &mut LogStreamLines,
    warn: &mut LogStreamLines,
    generator: &mut UniqueIdGenerator,
) -> std::result::Result<(), ToolError> {
    // The algorithm (283-291), on the worker count -threads asks for, as
    // TOPPBase applies the setting before main_ (TOPPBase.cpp:408-415).
    let options = Options {
        threads: ctx.thread_policy(),
        ..Options::default()
    };
    let mut finder = FeatureFinderAlgorithmPicked::with_options(options)?;
    let outcome = finder.run(experiment, features, parameters, seeds);
    let debug = finder.take_debug_output();
    if let Some(debug) = &debug {
        write_debug_log(debug)?;
    }
    // The console lines and debug stores, in the order the source makes them.
    for line in finder.report() {
        match line {
            ReportLine::Out(text) => writeln!(out, "{text}")?,
            ReportLine::Info(text) => info.line(out, text)?,
            ReportLine::Warn(text) => warn.line(err, text)?,
            ReportLine::StoreSeedMap(index) => {
                if let Some(seeds) = debug.as_ref().and_then(|d| d.seed_maps.get(*index)) {
                    let name = format!("debug/seeds_{}.featureXML", seeds.charge);
                    store_debug_features(&name, &seeds.map, info, out)?;
                }
            }
            ReportLine::StoreAbortReasons => {
                if let Some(map) = debug.as_ref().and_then(|d| d.abort_reasons.as_ref()) {
                    // `abort_map.setUniqueId()` draws from the generator
                    // the output ids come from later.
                    let mut map = map.clone();
                    map.unique_id = generator.get_unique_id();
                    store_debug_features("debug/abort_reasons.featureXML", &map, info, out)?;
                }
            }
            ReportLine::StoreInput => {
                if let Some(input) = debug.as_ref().and_then(|d| d.input.as_ref()) {
                    crate::format::path_io::write(Path::new("debug/input.mzML"), |writer| {
                        crate::format::mzml::write_source_float_arrays(writer, input)
                    })?;
                }
            }
        }
    }
    if let Some(debug) = &debug {
        for files in &debug.feature_files {
            crate::format::path_io::store(Path::new(&files.dta_name()), files.dta.as_bytes())?;
            if let Some(cropped) = &files.cropped_dta {
                crate::format::path_io::store(
                    Path::new(&files.cropped_dta_name()),
                    cropped.as_bytes(),
                )?;
            }
            crate::format::path_io::store(Path::new(&files.plot_name()), &files.plot)?;
        }
    }
    match outcome {
        Ok(()) => Ok(()),
        Err(error) => {
            if let Some(termination) = debug
                .as_ref()
                .and_then(|d| d.termination.as_ref())
                .filter(|t| t.kind == TerminationKind::Exception)
            {
                // The source process terminates here (std::terminate from an
                // exception that leaves the OpenMP region, then SIGABRT). The
                // port reports the exception as TOPPBase reports it where it
                // can catch it. Where the source dies from an out-of-bounds
                // access or never returns, the port's refusal is reported
                // below like any other error; either way the debug log holds
                // only what the executed process had flushed
                // (`write_debug_log`).
                let _ = error;
                return Err(ToolError::unexpected(&termination.message));
            }
            if seeds::is_length_error(&error) {
                // `std::length_error` is no `BaseException`: TOPPBase's outer
                // `catch (const std::exception&)` reports it
                // (TOPPBase.cpp:510-514), after the stack unwinding has
                // flushed and closed the debug log.
                return Err(ToolError::escaped(seeds::LENGTH_ERROR_WHAT));
            }
            if let Error::InvalidValue(message) = &error {
                // The algorithm's IllegalArgument and InvalidValue exceptions
                // reach TOPPBase's catch-all (TOPPBase.cpp:495-499).
                return Err(ToolError::unexpected(message));
            }
            Err(error.into())
        }
    }
}

/// Create `debug/features` and write `debug/log.txt`, as `run_` does when it
/// opens the stream (`FeatureFinderAlgorithmPicked.cpp:230-231`).
///
/// After a terminated run the file holds only what the source's file buffer
/// had written when the process died
/// ([`DebugTermination::log_file_bytes`](crate::analysis::feature_finder_picked::debug::DebugTermination::log_file_bytes);
/// the tool's instance runs once, so that is this run's
/// [`DebugLog::flushed_bytes`](crate::analysis::feature_finder_picked::debug::DebugLog::flushed_bytes)).
fn write_debug_log(debug: &DebugOutput) -> Result<()> {
    file::make_dir("debug/features")?;
    if debug.log_opened {
        let text = debug.log.text().as_bytes();
        let written = match debug.termination.as_ref() {
            Some(termination) => {
                let bytes = termination
                    .log_file_bytes
                    .unwrap_or_else(|| debug.log.flushed_bytes());
                &text[..bytes.min(text.len())]
            }
            None => text,
        };
        crate::format::path_io::store(Path::new("debug/log.txt"), written)?;
    }
    Ok(())
}

/// `FileHandler().storeFeatures(name, map)` with the line
/// `FeatureXMLFile::store` logs when ids are unassigned
/// (`FeatureXMLFile.cpp:80-88`): the map's own id and every feature's count.
fn store_debug_features(
    name: &str,
    map: &FeatureMap,
    info: &mut LogStreamLines,
    out: &mut dyn Write,
) -> Result<()> {
    let invalid = map.count_unique_ids(|id| usize::from(id == 0))?;
    if invalid > 0 {
        info.line(
            out,
            &format!("FeatureXMLHandler::store():  found {invalid} invalid unique ids"),
        )?;
    }
    FileHandler::store_feature_map(name, map, Some(FileType::FeatureXml))
}

/// The line cache of an OpenMS log stream (`LogStreamBuf`,
/// `LogStream.cpp:180-300`).
///
/// A line equal to one of the two the stream printed last is not printed
/// again but counted. When a new line pushes the older of the two out, and
/// that one was repeated, `<line> occurred N times` is printed first. Empty
/// lines bypass the cache. At the end of the process (`clearCache`) every
/// remaining repeated line is reported, in lexicographic order.
#[derive(Default)]
struct LogStreamLines {
    /// Line, repeat count and recency stamp.
    entries: Vec<(String, usize, u64)>,
    stamp: u64,
}

impl LogStreamLines {
    fn line(&mut self, out: &mut dyn Write, text: &str) -> Result<()> {
        if text.is_empty() {
            writeln!(out)?;
            return Ok(());
        }
        self.stamp += 1;
        if let Some(entry) = self.entries.iter_mut().find(|entry| entry.0 == text) {
            entry.1 += 1;
            entry.2 = self.stamp;
            return Ok(());
        }
        self.evict_oldest(out)?;
        self.entries.push((text.to_owned(), 0, self.stamp));
        writeln!(out, "{text}")?;
        Ok(())
    }

    /// `addToCache_`: with two lines cached, the older one leaves.
    fn evict_oldest(&mut self, out: &mut dyn Write) -> Result<()> {
        if self.entries.len() > 1 {
            let oldest = (0..self.entries.len())
                .min_by_key(|&index| self.entries[index].2)
                .unwrap_or(0);
            let (line, repeats, _) = self.entries.remove(oldest);
            if repeats != 0 {
                writeln!(out, "<{line}> occurred {} times", repeats + 1)?;
            }
        }
        Ok(())
    }

    /// The end of a successful run: TOPPBase's closing `<tool> took ...` info
    /// line enters the cache, which evicts the older cached line, then the
    /// stream's `clearCache` at exit.
    ///
    /// The lifecycle prints the closing line itself once the body has returned
    /// ([`run_with`](crate::cli::run_with)), after this. So the repeat count of
    /// the evicted line comes before the closing line, as in the source, but
    /// so does that of the line `clearCache` reports, which the source prints
    /// after it. The two orders differ only when both lines the cache still
    /// holds were repeated; no executed case has that (the console blocks
    /// `tests/topp_feature_finder_centroided.rs` compares end at the closing
    /// line and are equal).
    fn close(&mut self, out: &mut dyn Write) -> Result<()> {
        self.evict_oldest(out)?;
        self.entries.sort_by(|a, b| a.0.cmp(&b.0));
        for (line, repeats, _) in self.entries.drain(..) {
            if repeats != 0 {
                writeln!(out, "<{line}> occurred {} times", repeats + 1)?;
            }
        }
        Ok(())
    }
}
