// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The picked feature finder's parameters, input validation and entry point
//! (`FEATUREFINDER/FeatureFinderAlgorithmPicked.h`).
//!
//! The source class identifies features in an LC-MS map of centroided MS1
//! spectra: peptides whose isotope distribution shows over time. It computes
//! RT and m/z positions and a charge estimate for each. It finds pronounced
//! regions around *seeds*, which the caller may provide (for example from MS/MS
//! identifications) or the algorithm computes itself, and then fits an isotope
//! and retention-time model to each seed's data points.
//!
//! The source class holds its state across one `run`. This port splits it by
//! stage:
//!
//! - this module: the 29 parameter defaults
//!   ([`default_parameters`](crate::analysis::feature_finder_picked::algorithm::default_parameters)),
//!   the typed members of `updateMembers_` and the values `run_` reads
//!   ([`Settings`](crate::analysis::feature_finder_picked::algorithm::Settings)),
//!   the input checks of `run`
//!   ([`validate_input`](crate::analysis::feature_finder_picked::algorithm::validate_input)),
//!   the resource ceilings
//!   ([`Limits`](crate::analysis::feature_finder_picked::algorithm::Limits),
//!   [`Options`](crate::analysis::feature_finder_picked::algorithm::Options)) and
//!   the entry point [`run`](crate::analysis::feature_finder_picked::algorithm::run);
//! - [`crate::analysis::feature_finder_picked::scoring`]: the per-peak intensity,
//!   trace and isotope-pattern scores (steps 1, 2 and 3.1);
//! - [`crate::analysis::feature_finder_picked::seeds`]: the precalculated isotope
//!   patterns (step 2.5), seed selection (step 3.2) and
//!   [`SeedStage`](crate::analysis::feature_finder_picked::seeds::SeedStage),
//!   which runs everything up to and including seed selection;
//! - [`crate::analysis::feature_finder_picked::extension`]: the best isotope fit
//!   of a seed and the mass traces grown from it (step 3.3.1);
//! - [`crate::analysis::feature_finder_picked::fitting`]: the retention-time
//!   model, the cropping, the quality checks and the feature itself (steps
//!   3.3.2 to 3.3.5);
//! - [`crate::analysis::feature_finder_picked::resolution`]: overlap resolution
//!   and the apex annotation (step 4);
//! - this module again:
//!   [`feature_stage`](crate::analysis::feature_finder_picked::algorithm::feature_stage),
//!   the seed loop that drives those three and the source's single
//!   `#pragma omp parallel for`;
//! - [`crate::analysis::feature_finder_picked::instance`]: the stateful
//!   algorithm object, which [`run`](crate::analysis::feature_finder_picked::algorithm::run)
//!   uses with a fresh instance and a fresh feature map;
//! - [`crate::analysis::feature_finder_picked::debug`]: the `write_debug`
//!   output as data.
//!
//! The API mapping, the preserved source conventions, the native differences
//! and the evidence are in `docs/FEATURE_FINDER_PICKED_SUPPORT.md`.
//!
//! The source writes its seed counts to `std::cout` and its warnings to the
//! OpenMS log. Library code here never prints: every such line is collected in
//! the log of the returned value.

use std::collections::BTreeMap;

use crate::analysis::feature_finder_picked::debug::{AbortReasons, DebugOutput};
use crate::analysis::feature_finder_picked::instance::{
    ABORT_BLOCK_HEADING, Bookkeeping, DebugKey, FeatureFinderAlgorithmPicked, Progress,
    extend_charge, preflight_charge, settle_charge,
};
use crate::analysis::feature_finder_picked::resolution::{
    annotate_apex, invalid_apex_warning, resolve_overlaps,
};
use crate::analysis::feature_finder_picked::scoring::{libstdcxx, source_is_sorted, x86_64};
use crate::analysis::feature_finder_picked::seeds::SeedStage;
use crate::analysis::feature_finder_picked::source_sort::{
    TemporaryBuffer, source_sort_by, source_stable_sort_permutation,
};
use crate::analysis::feature_finder_picked::trace_fitter::TraceFitterParams;
use crate::concept::parallel::Threads;
use crate::kernel::{DataArray, Feature, FeatureMap, MSChromatogram, MSExperiment};
use crate::param::{DefaultParamHandler, Param, ParamValue, ParamValueType};
use crate::{Error, Result};

/// Name of the source parameter handler, `DefaultParamHandler("FeatureFinderAlgorithmPicked")`.
///
/// It prefixes the unknown-parameter warnings of [`Settings::from_parameters`].
pub const HANDLER_NAME: &str = "FeatureFinderAlgorithmPicked";

/// Isotope count of the precalculated averagine patterns before abundance
/// overrides: source `Size max_isotopes = 20` in `run_`.
pub const BASE_MAX_ISOTOPES: usize = 20;

/// Isotope count each non-default abundance adds: source `max_isotopes += 1000`,
/// commented `// Why?` in the source.
pub const OVERRIDE_EXTRA_ISOTOPES: usize = 1000;

/// The source parameter defaults: `FeatureFinderAlgorithmPicked()` followed by
/// `getDefaultParameters()`.
///
/// Returns all 29 entries with the source names, values, descriptions, numeric
/// and string restrictions and `advanced` tags, in the source's insertion order,
/// and the seven section descriptions. The `advanced` section has no
/// description in the source, and neither has it here. `write_debug` is the
/// string `"false"` restricted to `"true"` and `"false"`, which ParamXML writes as
/// a `bool` item.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] only if the parameter tree rejects one of the
/// fixed literals, which does not happen.
pub fn default_parameters() -> Result<Param> {
    let mut d = Param::new();
    let advanced = ["advanced".to_string()];
    let none: [String; 0] = [];

    // debugging
    d.set_value(
        "write_debug",
        ParamValue::String("false".into()),
        "When debug mode is activated, several files with intermediate results are written to the folder 'debug' (do not use in parallel mode).",
        &none,
    )?;
    d.set_valid_strings("write_debug", &["true".to_string(), "false".to_string()])?;
    // intensity
    d.set_value(
        "intensity:bins",
        ParamValue::Integer(10),
        "Number of bins per dimension (RT and m/z). The higher this value, the more local the intensity significance score is.\nThis parameter should be decreased, if the algorithm is used on small regions of a map.",
        &none,
    )?;
    d.set_min_int("intensity:bins", 1)?;
    d.set_section_description(
        "intensity",
        "Settings for the calculation of a score indicating if a peak's intensity is significant in the local environment (between 0 and 1)",
    )?;
    // mass trace search parameters
    d.set_value(
        "mass_trace:mz_tolerance",
        ParamValue::Float(0.03),
        "Tolerated m/z deviation of peaks belonging to the same mass trace.\nIt should be larger than the m/z resolution of the instrument.\nThis value must be smaller than that 1/charge_high!",
        &none,
    )?;
    d.set_min_float("mass_trace:mz_tolerance", 0.0)?;
    d.set_value(
        "mass_trace:min_spectra",
        ParamValue::Integer(10),
        "Number of spectra that have to show a similar peak mass in a mass trace.",
        &none,
    )?;
    d.set_min_int("mass_trace:min_spectra", 1)?;
    d.set_value(
        "mass_trace:max_missing",
        ParamValue::Integer(1),
        "Number of consecutive spectra where a high mass deviation or missing peak is acceptable.\nThis parameter should be well below 'min_spectra'!",
        &none,
    )?;
    d.set_min_int("mass_trace:max_missing", 0)?;
    d.set_value(
        "mass_trace:slope_bound",
        ParamValue::Float(0.1),
        "The maximum slope of mass trace intensities when extending from the highest peak.\nThis parameter is important to separate overlapping elution peaks.\nIt should be increased if feature elution profiles fluctuate a lot.",
        &none,
    )?;
    d.set_min_float("mass_trace:slope_bound", 0.0)?;
    d.set_section_description(
        "mass_trace",
        "Settings for the calculation of a score indicating if a peak is part of a mass trace (between 0 and 1).",
    )?;
    // isotopic pattern search parameters
    d.set_value(
        "isotopic_pattern:charge_low",
        ParamValue::Integer(1),
        "Lowest charge to search for.",
        &none,
    )?;
    d.set_min_int("isotopic_pattern:charge_low", 1)?;
    d.set_value(
        "isotopic_pattern:charge_high",
        ParamValue::Integer(4),
        "Highest charge to search for.",
        &none,
    )?;
    d.set_min_int("isotopic_pattern:charge_high", 1)?;
    d.set_value(
        "isotopic_pattern:mz_tolerance",
        ParamValue::Float(0.03),
        "Tolerated m/z deviation from the theoretical isotopic pattern.\nIt should be larger than the m/z resolution of the instrument.\nThis value must be smaller than that 1/charge_high!",
        &none,
    )?;
    d.set_min_float("isotopic_pattern:mz_tolerance", 0.0)?;
    d.set_value(
        "isotopic_pattern:intensity_percentage",
        ParamValue::Float(10.0),
        "Isotopic peaks that contribute more than this percentage to the overall isotope pattern intensity must be present.",
        &advanced,
    )?;
    d.set_min_float("isotopic_pattern:intensity_percentage", 0.0)?;
    d.set_max_float("isotopic_pattern:intensity_percentage", 100.0)?;
    d.set_value(
        "isotopic_pattern:intensity_percentage_optional",
        ParamValue::Float(0.1),
        "Isotopic peaks that contribute more than this percentage to the overall isotope pattern intensity can be missing.",
        &advanced,
    )?;
    d.set_min_float("isotopic_pattern:intensity_percentage_optional", 0.0)?;
    d.set_max_float("isotopic_pattern:intensity_percentage_optional", 100.0)?;
    d.set_value(
        "isotopic_pattern:optional_fit_improvement",
        ParamValue::Float(2.0),
        "Minimal percental improvement of isotope fit to allow leaving out an optional peak.",
        &advanced,
    )?;
    d.set_min_float("isotopic_pattern:optional_fit_improvement", 0.0)?;
    d.set_max_float("isotopic_pattern:optional_fit_improvement", 100.0)?;
    d.set_value(
        "isotopic_pattern:mass_window_width",
        ParamValue::Float(25.0),
        "Window width in Dalton for precalculation of estimated isotope distributions.",
        &advanced,
    )?;
    d.set_min_float("isotopic_pattern:mass_window_width", 1.0)?;
    d.set_max_float("isotopic_pattern:mass_window_width", 200.0)?;
    d.set_value(
        "isotopic_pattern:abundance_12C",
        ParamValue::Float(98.93),
        "Rel. abundance of the light carbon. Modify if labeled.",
        &advanced,
    )?;
    d.set_min_float("isotopic_pattern:abundance_12C", 0.0)?;
    d.set_max_float("isotopic_pattern:abundance_12C", 100.0)?;
    d.set_value(
        "isotopic_pattern:abundance_14N",
        ParamValue::Float(99.632),
        "Rel. abundance of the light nitrogen. Modify if labeled.",
        &advanced,
    )?;
    d.set_min_float("isotopic_pattern:abundance_14N", 0.0)?;
    d.set_max_float("isotopic_pattern:abundance_14N", 100.0)?;
    d.set_section_description(
        "isotopic_pattern",
        "Settings for the calculation of a score indicating if a peak is part of a isotopic pattern (between 0 and 1).",
    )?;
    // seed settings
    d.set_value(
        "seed:min_score",
        ParamValue::Float(0.8),
        "Minimum seed score a peak has to reach to be used as seed.\nThe seed score is the geometric mean of intensity score, mass trace score and isotope pattern score.\nIf your features show a large deviation from the averagene isotope distribution or from an gaussian elution profile, lower this score.",
        &none,
    )?;
    d.set_min_float("seed:min_score", 0.0)?;
    d.set_max_float("seed:min_score", 1.0)?;
    d.set_section_description(
        "seed",
        "Settings that determine which peaks are considered a seed",
    )?;
    // fitting settings
    d.set_value(
        "fit:max_iterations",
        ParamValue::Integer(500),
        "Maximum number of iterations of the fit.",
        &advanced,
    )?;
    d.set_min_int("fit:max_iterations", 1)?;
    d.set_section_description("fit", "Settings for the model fitting")?;
    // feature settings
    d.set_value(
        "feature:min_score",
        ParamValue::Float(0.7),
        "Feature score threshold for a feature to be reported.\nThe feature score is the geometric mean of the average relative deviation and the correlation between the model and the observed peaks.",
        &none,
    )?;
    d.set_min_float("feature:min_score", 0.0)?;
    d.set_max_float("feature:min_score", 1.0)?;
    d.set_value(
        "feature:min_isotope_fit",
        ParamValue::Float(0.8),
        "Minimum isotope fit of the feature before model fitting.",
        &advanced,
    )?;
    d.set_min_float("feature:min_isotope_fit", 0.0)?;
    d.set_max_float("feature:min_isotope_fit", 1.0)?;
    d.set_value(
        "feature:min_trace_score",
        ParamValue::Float(0.5),
        "Trace score threshold.\nTraces below this threshold are removed after the model fitting.\nThis parameter is important for features that overlap in m/z dimension.",
        &advanced,
    )?;
    d.set_min_float("feature:min_trace_score", 0.0)?;
    d.set_max_float("feature:min_trace_score", 1.0)?;
    d.set_value(
        "feature:min_rt_span",
        ParamValue::Float(0.333),
        "Minimum RT span in relation to extended area that has to remain after model fitting.",
        &advanced,
    )?;
    d.set_min_float("feature:min_rt_span", 0.0)?;
    d.set_max_float("feature:min_rt_span", 1.0)?;
    d.set_value(
        "feature:max_rt_span",
        ParamValue::Float(2.5),
        "Maximum RT span in relation to extended area that the model is allowed to have.",
        &advanced,
    )?;
    d.set_min_float("feature:max_rt_span", 0.5)?;
    d.set_value(
        "feature:rt_shape",
        ParamValue::String("symmetric".into()),
        "Choose model used for RT profile fitting. If set to symmetric a gauss shape is used, in case of asymmetric an EGH shape is used.",
        &advanced,
    )?;
    d.set_valid_strings(
        "feature:rt_shape",
        &["symmetric".to_string(), "asymmetric".to_string()],
    )?;
    d.set_value(
        "feature:max_intersection",
        ParamValue::Float(0.35),
        "Maximum allowed intersection of features.",
        &advanced,
    )?;
    d.set_min_float("feature:max_intersection", 0.0)?;
    d.set_max_float("feature:max_intersection", 1.0)?;
    d.set_value(
        "feature:reported_mz",
        ParamValue::String("monoisotopic".into()),
        "The mass type that is reported for features.\n'maximum' returns the m/z value of the highest mass trace.\n'average' returns the intensity-weighted average m/z value of all contained peaks.\n'monoisotopic' returns the monoisotopic m/z value derived from the fitted isotope model.",
        &none,
    )?;
    d.set_valid_strings(
        "feature:reported_mz",
        &[
            "maximum".to_string(),
            "average".to_string(),
            "monoisotopic".to_string(),
        ],
    )?;
    d.set_section_description(
        "feature",
        "Settings for the features (intensity, quality assessment, ...)",
    )?;
    // user-specified seed settings
    d.set_value(
        "user-seed:rt_tolerance",
        ParamValue::Float(5.0),
        "Allowed RT deviation of seeds from the user-specified seed position.",
        &none,
    )?;
    d.set_min_float("user-seed:rt_tolerance", 0.0)?;
    d.set_value(
        "user-seed:mz_tolerance",
        ParamValue::Float(1.1),
        "Allowed m/z deviation of seeds from the user-specified seed position.",
        &none,
    )?;
    d.set_min_float("user-seed:mz_tolerance", 0.0)?;
    d.set_value(
        "user-seed:min_score",
        ParamValue::Float(0.5),
        "Overwrites 'seed:min_score' for user-specified seeds. The cutoff is typically a bit lower in this case.",
        &none,
    )?;
    d.set_min_float("user-seed:min_score", 0.0)?;
    d.set_max_float("user-seed:min_score", 1.0)?;
    d.set_section_description("user-seed", "Settings for user-specified seeds.")?;
    // advanced/debugging settings
    d.set_value(
        "advanced:pseudo_rt_shift",
        ParamValue::Float(500.0),
        "Pseudo RT shift used when .",
        &advanced,
    )?;
    d.set_min_float("advanced:pseudo_rt_shift", 1.0)?;
    Ok(d)
}

/// The retention-time model of the trace fit: parameter `feature:rt_shape`.
///
/// Source `chooseTraceFitter_` (`FeatureFinderAlgorithmPicked.cpp:1897-1912`)
/// selects the fitter from it:
/// [`FittedModel::new`](crate::analysis::feature_finder_picked::fitting::FittedModel::new)
/// is its counterpart.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RtShape {
    /// `symmetric`: a Gaussian (`GaussTraceFitter`), the default.
    Symmetric,
    /// `asymmetric`: an exponential-Gaussian hybrid (`EGHTraceFitter`).
    Asymmetric,
}

/// The m/z reported for a feature: parameter `feature:reported_mz`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReportedMz {
    /// `maximum`: the average m/z of the trace with the highest theoretical
    /// intensity.
    Maximum,
    /// `average`: the intensity-weighted average m/z of all contained peaks.
    Average,
    /// `monoisotopic`: the monoisotopic m/z derived from the fitted isotope
    /// model, the default.
    Monoisotopic,
}

/// What a run does when an intensity bin step is zero or infinite.
///
/// Step 1 of source `run_` (`FeatureFinderAlgorithmPicked.cpp:244-245`)
/// divides the MS1 retention-time and m/z extents by `intensity:bins` without
/// a check. The step is zero when every MS1 spectrum has the same retention
/// time, when every MS1 peak has the same m/z, or when a subnormal extent
/// underflows in the division (`4.9e-324 / 2` is zero); the retention-time step
/// is infinite when the extent overflows (retention times from `-1e308` to
/// `1e308`), and either step is infinite when a coordinate is. The bins are
/// still computed
/// ([`IntensityThresholds::compute`](crate::analysis::feature_finder_picked::scoring::IntensityThresholds::compute)),
/// but `intensityScore_` (`:1837-1838`) then converts `floor(NaN)` or
/// `floor(inf)` to `UInt` for every peak, which is undefined behaviour.
///
/// The Linux x86_64 Release build compiles that conversion as `cvttsd2si`
/// into a 64-bit register, keeps the low 32 bits and caps them, so the
/// selected cells stay inside the grid; the distances to the bin centres are
/// `0 / 0` or `inf / inf` whatever cell was selected. Every intensity score is
/// therefore the default NaN, every overall score the seed loop computes is
/// NaN, no peak becomes a seed, and the run returns an empty feature map with
/// the source's log lines (`Found 0 seeds`, `Found 0 feature candidates` per
/// charge). Executed against `openms4-release-bc9cc12-c19e494-174b576`, three
/// repetitions at one and at four threads, identical: FeatureFinderCentroided_1
/// with every retention time equal, with every m/z equal, with a subnormal and
/// with an overflowing retention-time extent, each with the FFC_1 and the
/// default parameters and with `seed:min_score` 0, and the class-level score
/// arrays of each (`../oracle/ffap-sem-completion`). The outcome is explained
/// by the emitted instructions and reproduced by
/// [`IntensityThresholds::score`](crate::analysis::feature_finder_picked::scoring::IntensityThresholds::score).
///
/// An input with at most `2 * min_spectra` scans (source `min_spectra_`, half
/// of `mass_trace:min_spectra`) never reaches the seed loop, so its intensity
/// scores are never read: both variants return the source's empty result for
/// it, whatever the bin steps.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum DegenerateBinStep {
    /// Compute what the Linux x86_64 Release build computes: NaN intensity
    /// scores, no seed, an empty map. The default, and what the
    /// FeatureFinderCentroided tool uses.
    #[default]
    Source,
    /// Return [`Error::InvalidValue`] before any work when a bin step is zero
    /// or infinite and the seed loop visits at least one scan, the only case in
    /// which the undefined scores are read.
    Refuse,
}

/// How a non-default `isotopic_pattern:abundance_12C` or `abundance_14N` is
/// handled.
///
/// The source (`FeatureFinderAlgorithmPicked.cpp:163-179`) builds each override
/// by inserting the light and heavy isotope into a default-constructed
/// `IsotopeDistribution`, which already holds the peak `(0, 1)`. The executed C++
/// keeps that stray peak: the patterns grow (the first FFC_1 window has 27
/// normalised bins instead of 6) and FeatureFinderCentroided_1 with
/// `abundance_12C = 90` finds no seed and no feature at all. The native
/// [`CoarseIsotopePatternGenerator::set_isotope_override`](crate::chemistry::isotopes::CoarseIsotopePatternGenerator::set_isotope_override)
/// rejects such a distribution, so the defect cannot be reproduced.
///
/// The port therefore computes the **intended** two-isotope override by
/// default, rather than refusing a parameter the source accepts. That is a
/// deliberate divergence from the executed C++ and the only one in this module
/// that changes which features are found; it is recorded as `CPP-247` in
/// `docs/FEATURE_FINDER_PICKED_SUPPORT.md` together with the measured
/// difference.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum AbundanceOverride {
    /// Use the intended two-isotope distribution: weights `a / 100` and
    /// `1 - a / 100` for the light and heavy isotope, narrowed to `f32` under
    /// source precision. This is not the executed C++ result; see the type
    /// documentation. This is the default.
    #[default]
    Intended,
    /// A non-default abundance is [`Error::Unsupported`], so that a caller who
    /// must not diverge from the executed C++ can refuse instead of differing.
    Refuse,
}

/// Which parameter `writeFeatureDebugInfo_` reads for its pseudo-RT shift in a
/// `write_debug` run.
///
/// The source reads `debug:pseudo_rt_shift` (`FeatureFinderAlgorithmPicked.cpp:2137`),
/// but declares `advanced:pseudo_rt_shift` (`:124`). Unless the caller passes
/// the undeclared key itself, `Param::getValue` throws `ElementNotFound` inside
/// the OpenMP region of the seed loop, and the process terminates (executed:
/// `FeatureFinderCentroided` with `-algorithm:write_debug` on
/// FeatureFinderCentroided_1 is killed by `SIGABRT`, shell status 134). A safe
/// port cannot end the process abnormally; it returns an error at that seed,
/// and the FeatureFinderCentroided tool exits with code 8 after writing what
/// the executed process had written.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PseudoRtShiftKey {
    /// Read `debug:pseudo_rt_shift`, as the source does: an integer or float
    /// value is used; a string or list value is the bits of a heap pointer
    /// ([`PseudoRtShift::HeapAddress`](crate::analysis::feature_finder_picked::debug::PseudoRtShift::HeapAddress)),
    /// emulated wherever the written text does not depend on the address. At
    /// the first seed that reaches the fit without the key or with an empty
    /// value, where the source process terminates, the run returns
    /// [`Error::Unsupported`] and records the point in
    /// [`DebugOutput::termination`]. The default, and what the
    /// FeatureFinderCentroided tool uses.
    #[default]
    Source,
    /// Read the declared `advanced:pseudo_rt_shift` (default 500) and write
    /// the member's files for every seed that reaches the fit: what the source
    /// evidently intends.
    Declared,
}

/// What [`FeatureFinderAlgorithmPicked::parameters`] returns after
/// [`FeatureFinderAlgorithmPicked::set_parameters`] (or a run) refused a
/// parameter set.
///
/// Source `DefaultParamHandler::setParameters` (`DefaultParamHandler.cpp`)
/// assigns the new set, merged with the defaults, to `param_` *before*
/// `Param::checkDefaults` throws `InvalidParameter` for a value of the wrong
/// type or outside its restriction, and calls `updateMembers_` only after the
/// check. After the exception `getParameters()` therefore returns the rejected
/// set while the typed members keep the values of the last accepted one
/// (executed: `params_after_failed_set.txt` and `rejected_stdout.txt` of the
/// oracle). The next `setParameters` or `run` replaces the set again.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum RejectedParameters {
    /// Show the rejected set, merged with the defaults, as the source does;
    /// the settings a run uses stay those of the last accepted set, and the
    /// next accepted set replaces it. The default.
    #[default]
    Shown,
    /// Keep showing the last accepted set: a refused set changes nothing, the
    /// port's usual atomicity.
    Discarded,
}

/// Resource ceilings of the seed stage, checked before the corresponding work.
///
/// The source has no ceilings. Each limit below is far above the workloads the
/// source is used for; the FeatureFinderCentroided_1 input (112 spectra, 3,084
/// peaks, one charge) stays several orders of magnitude below every one.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Limits {
    /// Most spectra in the input.
    pub max_spectra: usize,
    /// Most peaks in the input, over all spectra.
    pub max_peaks: usize,
    /// Most charges between `charge_low` and `charge_high`.
    pub max_charges: usize,
    /// Most intensity bins per dimension (`intensity:bins`); the quantile
    /// table holds `bins^2 * 21` values.
    pub max_intensity_bins: usize,
    /// Most precalculated isotope windows,
    /// `ceil(max_mz * charge_high / mass_window_width) + 1`.
    pub max_isotope_windows: usize,
    /// Most isotope-pattern values over all windows, bounded before step 2.5 as
    /// the window count times the isotope count of each pattern.
    pub max_pattern_values: usize,
    /// Most bytes of per-peak score arrays, `4 * peaks * (3 + 2 * charges)`.
    pub max_score_bytes: usize,
    /// Most work units of the scoring steps: one per spectrum visited while
    /// binning, per nearest-peak search, per step of the linear isotope walk
    /// and per value a correlation reads.
    pub max_work: u64,
    /// Most seeds of one charge that the seed loop extends.
    pub max_seeds: usize,
    /// Most work units of the seed loop, bounded before it starts: per charge,
    /// the seed count times `isotopes * (isotopes + spectra)`, an upper bound
    /// on the isotope search around each seed plus the extension of each
    /// isotope's trace through the scans.
    pub max_seed_work: u64,
    /// Most bytes of the `write_debug` log kept in memory. The source writes
    /// the log to a file and has no bound.
    pub max_debug_bytes: usize,
}

impl Limits {
    /// Default [`Self::max_spectra`].
    pub const DEFAULT_MAX_SPECTRA: usize = 10_000_000;
    /// Default [`Self::max_peaks`].
    pub const DEFAULT_MAX_PEAKS: usize = 1_000_000_000;
    /// Default [`Self::max_charges`].
    pub const DEFAULT_MAX_CHARGES: usize = 1_000;
    /// Default [`Self::max_intensity_bins`].
    pub const DEFAULT_MAX_INTENSITY_BINS: usize = 2_000;
    /// Default [`Self::max_isotope_windows`].
    pub const DEFAULT_MAX_ISOTOPE_WINDOWS: usize = 1_000_000;
    /// Default [`Self::max_pattern_values`].
    pub const DEFAULT_MAX_PATTERN_VALUES: usize = 200_000_000;
    /// Default [`Self::max_score_bytes`], 16 GiB.
    pub const DEFAULT_MAX_SCORE_BYTES: usize = 16 << 30;
    /// Default [`Self::max_work`].
    pub const DEFAULT_MAX_WORK: u64 = 10_000_000_000_000;
    /// Default [`Self::max_seeds`]: far above the 800,000 features the port's
    /// resource contract admits.
    pub const DEFAULT_MAX_SEEDS: usize = 50_000_000;
    /// Default [`Self::max_seed_work`]: 44,000 scans with 800,000 seeds of 20
    /// isotopes stay an order of magnitude below it.
    pub const DEFAULT_MAX_SEED_WORK: u64 = 10_000_000_000_000;
    /// Default [`Self::max_debug_bytes`], 16 GiB; the FeatureFinderCentroided_1
    /// debug log has about 1.1 MiB.
    pub const DEFAULT_MAX_DEBUG_BYTES: usize = 16 << 30;
}

impl Default for Limits {
    fn default() -> Self {
        Self {
            max_spectra: Self::DEFAULT_MAX_SPECTRA,
            max_peaks: Self::DEFAULT_MAX_PEAKS,
            max_charges: Self::DEFAULT_MAX_CHARGES,
            max_intensity_bins: Self::DEFAULT_MAX_INTENSITY_BINS,
            max_isotope_windows: Self::DEFAULT_MAX_ISOTOPE_WINDOWS,
            max_pattern_values: Self::DEFAULT_MAX_PATTERN_VALUES,
            max_score_bytes: Self::DEFAULT_MAX_SCORE_BYTES,
            max_work: Self::DEFAULT_MAX_WORK,
            max_seeds: Self::DEFAULT_MAX_SEEDS,
            max_seed_work: Self::DEFAULT_MAX_SEED_WORK,
            max_debug_bytes: Self::DEFAULT_MAX_DEBUG_BYTES,
        }
    }
}

/// Native options of a run.
///
/// The defaults reproduce the source where it is defined. Where it is
/// undefined they reproduce the Linux x86_64 Release build when its outcome is
/// fixed and explained by the emitted instructions ([`DegenerateBinStep`]), and
/// refuse otherwise. [`AbundanceOverride`] is the one designed difference.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Options {
    /// Resource ceilings.
    pub limits: Limits,
    /// Handling of non-default isotope abundances.
    pub abundance_override: AbundanceOverride,
    /// Handling of a zero or infinite intensity bin step.
    pub degenerate_bin_step: DegenerateBinStep,
    /// Worker threads of the seed loop, the port's counterpart of the source's
    /// `#pragma omp parallel for` and of the TOPP `-threads` parameter.
    ///
    /// The result does not depend on it: the loop uses
    /// [`crate::concept::parallel::map_collect`], whose results keep the input
    /// order, and everything after the loop is serial. The default is every
    /// available core.
    pub threads: Threads,
    /// The parameter a `write_debug` run reads for the pseudo-RT shift of its
    /// feature plots.
    pub pseudo_rt_shift: PseudoRtShiftKey,
    /// What the algorithm instance's parameters show after a rejected
    /// parameter set.
    pub rejected_parameters: RejectedParameters,
}

/// The typed parameter values of one run.
///
/// The first group are the members the source's `updateMembers_` sets; the
/// second are the values `run_` reads from `param_` directly. Percentages are
/// already divided by 100, as in the source.
#[derive(Clone, Debug, PartialEq)]
pub struct Settings {
    /// `isotopic_pattern:mz_tolerance` (source `pattern_tolerance_`).
    pub pattern_tolerance: f64,
    /// `mass_trace:mz_tolerance` (source `trace_tolerance_`).
    pub trace_tolerance: f64,
    /// Half of `mass_trace:min_spectra`, rounded down (source `min_spectra_`):
    /// the number of scans inspected on each side of a peak.
    pub min_spectra: usize,
    /// `mass_trace:max_missing` (source `max_missing_trace_peaks_`).
    pub max_missing_trace_peaks: u32,
    /// `mass_trace:slope_bound` (source `slope_bound_`).
    pub slope_bound: f64,
    /// `isotopic_pattern:intensity_percentage / 100` (source
    /// `intensity_percentage_`): pattern peaks above this contribution are
    /// required.
    pub intensity_percentage: f64,
    /// `isotopic_pattern:intensity_percentage_optional / 100` (source
    /// `intensity_percentage_optional_`): the trimming cutoff of the
    /// precalculated patterns.
    pub intensity_percentage_optional: f64,
    /// `isotopic_pattern:optional_fit_improvement / 100` (source
    /// `optional_fit_improvement_`).
    pub optional_fit_improvement: f64,
    /// `isotopic_pattern:mass_window_width` in Da (source `mass_window_width_`).
    pub mass_window_width: f64,
    /// `intensity:bins` (source `intensity_bins_`).
    pub intensity_bins: usize,
    /// `feature:min_isotope_fit` (source `min_isotope_fit_`).
    pub min_isotope_fit: f64,
    /// `feature:min_trace_score` (source `min_trace_score_`).
    pub min_trace_score: f64,
    /// `feature:min_rt_span` (source `min_rt_span_`).
    pub min_rt_span: f64,
    /// `feature:max_rt_span` (source `max_rt_span_`).
    pub max_rt_span: f64,
    /// `feature:max_intersection` (source `max_feature_intersection_`).
    pub max_feature_intersection: f64,
    /// `feature:reported_mz` (source `reported_mz_`).
    pub reported_mz: ReportedMz,
    /// `feature:min_score`.
    pub min_feature_score: f64,
    /// `isotopic_pattern:charge_low`.
    pub charge_low: i32,
    /// `isotopic_pattern:charge_high`.
    pub charge_high: i32,
    /// `fit:max_iterations`, handed to the trace fitter as `max_iteration`.
    pub max_iterations: u32,
    /// `isotopic_pattern:abundance_12C` in percent.
    pub abundance_12c: f64,
    /// `isotopic_pattern:abundance_14N` in percent.
    pub abundance_14n: f64,
    /// Whether `abundance_12C` differs from its default, compared as the source
    /// compares the two `ParamValue`s: exact equality of type and value.
    pub abundance_12c_changed: bool,
    /// Whether `abundance_14N` differs from its default.
    pub abundance_14n_changed: bool,
    /// `seed:min_score`.
    pub seed_min_score: f64,
    /// `user-seed:rt_tolerance` in seconds.
    pub user_seed_rt_tolerance: f64,
    /// `user-seed:mz_tolerance` in Th.
    pub user_seed_mz_tolerance: f64,
    /// `user-seed:min_score`.
    pub user_seed_min_score: f64,
    /// `write_debug`.
    pub write_debug: bool,
    /// `feature:rt_shape`.
    pub rt_shape: RtShape,
}

/// `sizeof(MSSpectrum::FloatDataArray)` in the Linux x86_64 Release build:
/// the `std::vector<float>` base and `MetaInfoDescription` of one score array,
/// 88 bytes (executed: `../oracle/ffap-complete-fix6`, the driver's `sizes`
/// mode; the array vector's `max_size()` is `(2^63 - 1) / 88`, far above every
/// 32-bit count, so `resize` never throws `std::length_error` here).
pub(crate) const SOURCE_FLOAT_DATA_ARRAY_BYTES: u64 = 88;

/// The largest first-spectrum score-array allocation for which
/// [`Settings::score_arrays_overrun`] records a termination: 1 GiB
/// (12,201,611 arrays).
///
/// The source's allocation succeeds or throws `std::bad_alloc` depending on
/// the memory available to the process, which no deterministic port can
/// reproduce, so this is a documented line and not a property of the source:
/// below it every executed configuration reached the out-of-bounds write, and
/// above it the executed outcome was measured both ways at one count
/// (`../oracle/ffap-complete-fix6`). It is a crate constant, not a
/// [`Limits`] field, because a caller must not be able to move where a
/// termination is recorded (lead decision D12).
pub(crate) const SCORE_ARRAY_TERMINATION_CEILING_BYTES: u64 = 1 << 30;

impl Settings {
    /// Apply `parameters` over the defaults and read the typed values: source
    /// `setParameters(param)` with `updateMembers_`, plus the reads at the start
    /// of `run_`.
    ///
    /// Missing entries take their defaults. Entries the defaults do not know are
    /// kept and reported as warnings in the returned vector, as the source
    /// reports them on the OpenMS warning log; a caller that must be strict
    /// treats a warning as an error, as TOPPBase does for its INI files.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a known entry has a different value
    /// type than its default or violates its restriction (source
    /// `Exception::InvalidParameter`, with its text; [`check_parameters`]), or
    /// when a value cannot be converted as the source converts it (a negative
    /// `intensity:bins`, `mass_trace:max_missing` or `fit:max_iterations`
    /// that passes its restriction, [`NEGATIVE_UNSIGNED_WHAT`]).
    pub fn from_parameters(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(HANDLER_NAME)?;
        let defaults = default_parameters()?;
        handler.set_defaults(defaults.clone())?;
        handler.defaults_to_parameters()?;
        // The source's checks, in place of the shared handler's.
        handler.set_check_defaults(false);
        handler.set_parameters(parameters)?;
        let merged = handler.parameters();
        let mut unknown = Vec::new();
        check_parameters(merged, &defaults, HANDLER_NAME, &mut unknown)?;
        let settings = Self::read(merged, &defaults)?;
        check_run_conversions(merged)?;
        let warnings = unknown
            .into_iter()
            .map(|key| format!("{HANDLER_NAME}: unknown parameter '{key}'"))
            .collect();
        Ok((settings, warnings))
    }

    /// The typed values of a checked parameter set: [`Self::update_members`]
    /// on the defaults' values, then [`Self::read_run_values`].
    pub(crate) fn read(p: &Param, defaults: &Param) -> Result<Self> {
        let mut settings = Self {
            pattern_tolerance: 0.0,
            trace_tolerance: 0.0,
            min_spectra: 0,
            max_missing_trace_peaks: 0,
            slope_bound: 0.0,
            intensity_percentage: 0.0,
            intensity_percentage_optional: 0.0,
            optional_fit_improvement: 0.0,
            mass_window_width: 0.0,
            intensity_bins: 0,
            min_isotope_fit: 0.0,
            min_trace_score: 0.0,
            min_rt_span: 0.0,
            max_rt_span: 0.0,
            max_feature_intersection: 0.0,
            reported_mz: ReportedMz::Maximum,
            min_feature_score: 0.0,
            charge_low: 0,
            charge_high: 0,
            max_iterations: 0,
            abundance_12c: 0.0,
            abundance_14n: 0.0,
            abundance_12c_changed: false,
            abundance_14n_changed: false,
            seed_min_score: 0.0,
            user_seed_rt_tolerance: 0.0,
            user_seed_mz_tolerance: 0.0,
            user_seed_min_score: 0.0,
            write_debug: false,
            rt_shape: RtShape::Symmetric,
        };
        settings.update_members(p)?;
        settings.read_run_values(p, defaults)?;
        Ok(settings)
    }

    /// Source `updateMembers_` (`FeatureFinderAlgorithmPicked.cpp:1108-1126`):
    /// the members in the source's order, each converted as the Linux x86_64
    /// Release build converts it.
    ///
    /// - `min_spectra_ = (UInt) std::floor((double) value * 0.5)`: the whole
    ///   64-bit value as a `double`, then `cvttsd2si` into a 64-bit register and
    ///   its low 32 bits (`libOpenMS.so` `0x18da803`-`0x18da866`).
    /// - `max_missing_trace_peaks_` and `intensity_bins_` through
    ///   `ParamValue::operator unsigned int` (`ParamValue.cpp:451-461`,
    ///   `0x7aee30`): a negative value throws `Exception::ConversionError`,
    ///   any other keeps its low 32 bits.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] with [`NEGATIVE_UNSIGNED_WHAT`] for a
    /// negative `mass_trace:max_missing` or `intensity:bins`. The members
    /// before it keep their new values and the rest their old ones, as in the
    /// source, where `updateMembers_` throws half way.
    pub(crate) fn update_members(&mut self, p: &Param) -> Result<()> {
        let float = |key: &str| p.value(key)?.to_f64();
        self.pattern_tolerance = float("isotopic_pattern:mz_tolerance")?;
        self.trace_tolerance = float("mass_trace:mz_tolerance")?;
        let min_spectra = (source_i64(p, "mass_trace:min_spectra")? as f64 * 0.5).floor();
        self.min_spectra = x86_64::truncate_to_u32(min_spectra) as usize;
        self.max_missing_trace_peaks = source_unsigned(p, "mass_trace:max_missing")?;
        self.slope_bound = float("mass_trace:slope_bound")?;
        self.intensity_percentage = float("isotopic_pattern:intensity_percentage")? / 100.0;
        self.intensity_percentage_optional =
            float("isotopic_pattern:intensity_percentage_optional")? / 100.0;
        self.optional_fit_improvement = float("isotopic_pattern:optional_fit_improvement")? / 100.0;
        self.mass_window_width = float("isotopic_pattern:mass_window_width")?;
        self.intensity_bins = source_unsigned(p, "intensity:bins")? as usize;
        self.min_isotope_fit = float("feature:min_isotope_fit")?;
        self.min_trace_score = float("feature:min_trace_score")?;
        self.min_rt_span = float("feature:min_rt_span")?;
        self.max_rt_span = float("feature:max_rt_span")?;
        self.max_feature_intersection = float("feature:max_intersection")?;
        self.reported_mz = match p.value("feature:reported_mz")?.as_str()? {
            "maximum" => ReportedMz::Maximum,
            "average" => ReportedMz::Average,
            _ => ReportedMz::Monoisotopic,
        };
        Ok(())
    }

    /// The values source `run_` reads from `param_` itself, in its order.
    ///
    /// `charge_low` and `charge_high` are `(Int)` conversions, the low 32 bits
    /// of the value. `fit:max_iterations` is an `operator unsigned int`
    /// conversion, which throws for a negative value at the start of `run_`
    /// (`:152`); this keeps its low 32 bits, and [`check_run_conversions`]
    /// reports the throw where the run starts.
    pub(crate) fn read_run_values(&mut self, p: &Param, defaults: &Param) -> Result<()> {
        let float = |key: &str| p.value(key)?.to_f64();
        self.min_feature_score = float("feature:min_score")?;
        self.charge_low = source_i64(p, "isotopic_pattern:charge_low")? as i32;
        self.charge_high = source_i64(p, "isotopic_pattern:charge_high")? as i32;
        self.max_iterations = source_i64(p, "fit:max_iterations")? as u32;
        self.abundance_12c = float("isotopic_pattern:abundance_12C")?;
        self.abundance_14n = float("isotopic_pattern:abundance_14N")?;
        let changed = |key: &str| -> Result<bool> { Ok(p.value(key)? != defaults.value(key)?) };
        self.abundance_12c_changed = changed("isotopic_pattern:abundance_12C")?;
        self.abundance_14n_changed = changed("isotopic_pattern:abundance_14N")?;
        self.user_seed_rt_tolerance = float("user-seed:rt_tolerance")?;
        self.user_seed_mz_tolerance = float("user-seed:mz_tolerance")?;
        self.user_seed_min_score = float("user-seed:min_score")?;
        self.write_debug = p.value("write_debug")?.to_bool()?;
        self.seed_min_score = float("seed:min_score")?;
        // Source: `param_.getValue("feature:rt_shape") == "asymmetric"`, else symmetric.
        self.rt_shape = if p.value("feature:rt_shape")? == &ParamValue::String("asymmetric".into())
        {
            RtShape::Asymmetric
        } else {
            RtShape::Symmetric
        };
        Ok(())
    }

    /// The number of charges searched, `charge_high - charge_low + 1`, or zero
    /// when `charge_low` exceeds `charge_high` by one.
    ///
    /// The source computes `UInt charge_count = charge_high - charge_low + 1`
    /// and resizes every spectrum's float data arrays to the `UInt` value
    /// `3 + 2 * charge_count` before it writes arrays 0, 1 and 2 and the
    /// pattern and overall arrays `[3, 3 + charge_count)` and
    /// `[3 + charge_count, 3 + 2 * charge_count)`, both bounds in `UInt`
    /// (`FeatureFinderAlgorithmPicked.cpp:196-221`). Whether that wraps
    /// depends on the count `n` alone:
    ///
    /// - `0 <= n <= 2^31 - 2`: nothing wraps. `n = 0` (`charge_low =
    ///   charge_high + 1`) gives three arrays and an empty result, as executed;
    ///   larger counts are bounded by [`Limits::max_charges`], a native ceiling
    ///   in front of an allocation of up to `2^32 - 1` arrays per spectrum,
    ///   whose success depends on memory (executed: `std::bad_alloc` at
    ///   `n = 2^31 - 2`).
    /// - `n = 2^31 - 1` (`charge_low = 1`, `charge_high = INT_MAX`, the only
    ///   positive count that wraps): the size wraps to 1 and the source writes
    ///   arrays 1 and 2 past the end (executed: SIGSEGV).
    /// - `n = -1` (`charge_low = charge_high + 2`): the size wraps to 1, the
    ///   same out-of-bounds write (executed: SIGSEGV).
    /// - `n = -2` and `n = -3`: the size wraps to `2^32 - 1` or `2^32 - 3` and
    ///   every later write stays in bounds, so the outcome depends only on
    ///   whether that allocation succeeds (executed: `std::bad_alloc`).
    /// - `n <= -4`: the size wraps to `2^32 + 3 + 2n`, at least 9, and the
    ///   pattern loop writes up to index `2^32 + 2 + n`, past it; the source
    ///   writes out of bounds wherever the allocation succeeds (executed:
    ///   SIGSEGV for `charge_low = INT_MAX` with `charge_high` 1 and 498, and
    ///   `std::bad_alloc` where the wrapped size is near `2^32`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for every count that wraps, whatever the
    /// [`Limits`]: undefined behaviour (an out-of-bounds write) for
    /// `n = 2^31 - 1`, `n = -1` and `n <= -4`, and for `n = -2` and `n = -3`
    /// an allocation of billions of arrays per spectrum that the port does not
    /// attempt.
    pub fn charge_count(&self) -> Result<usize> {
        let count = i64::from(self.charge_high) - i64::from(self.charge_low) + 1;
        if count == i64::from(i32::MAX) {
            return Err(Error::InvalidValue(format!(
                "isotopic_pattern:charge_low {} and charge_high {}: the source's UInt score-array \
                 count 3 + 2 * {count} wraps to 1 and the source writes past the arrays; the \
                 behaviour is undefined",
                self.charge_low, self.charge_high
            )));
        }
        if count == -2 || count == -3 {
            let size = (3u32).wrapping_add((count as u32).wrapping_mul(2));
            return Err(Error::InvalidValue(format!(
                "isotopic_pattern:charge_low {} exceeds charge_high {} by more than one; the \
                 source's UInt score-array count wraps to {size} arrays per spectrum, whose \
                 allocation depends on memory, and the port does not attempt it",
                self.charge_low, self.charge_high
            )));
        }
        if count < 0 {
            return Err(Error::InvalidValue(format!(
                "isotopic_pattern:charge_low {} exceeds charge_high {} by more than one; the \
                 source's UInt score-array count wraps and the source writes past the arrays; \
                 the behaviour is undefined",
                self.charge_low, self.charge_high
            )));
        }
        usize::try_from(count).map_err(|_| Error::InvalidValue("charge count overflow".into()))
    }

    /// Whether the source process ends at the score arrays: the wrapped count
    /// is one [`Self::charge_count`] refuses as an out-of-bounds write, and the
    /// allocation before that write succeeds.
    ///
    /// The source resizes one spectrum's float data arrays at a time
    /// (`.cpp:196-221`) and walks past the end of the first one it fills, so
    /// the only allocation between the wrap and the out-of-bounds write is
    /// that spectrum's `(3 + 2n) mod 2^32` arrays of
    /// [`SOURCE_FLOAT_DATA_ARRAY_BYTES`] each.
    ///
    /// - `n = 2^31 - 1` and `n = -1` wrap to one array, 88 bytes: the write is
    ///   always reached (executed: SIGSEGV for 1/`INT_MAX` and 4/2).
    /// - `n <= -4` wraps to `2^32 + 3 + 2n` arrays, between 9 (`INT_MAX`/1)
    ///   and `2^32 - 5` (7/2, 352 GiB). Whether that allocation succeeds
    ///   depends on the memory available to the process, which the port cannot
    ///   reproduce (lead decision D6's rule for allocations), so it records the
    ///   termination up to [`SCORE_ARRAY_TERMINATION_CEILING_BYTES`] and
    ///   nothing above.
    ///
    /// Executed on the reference node (`../oracle/ffap-complete-fix6`, every
    /// case twice and identical): under a 16 GB address space SIGSEGV at 1, 9,
    /// 1003, 2003, 2005, 12,201,611 (1 GiB), 12,201,613, 20,000,003 and
    /// 40,000,003 arrays, and `std::bad_alloc` from 80,000,003 (6.6 GiB) up,
    /// `2^32 - 5` arrays (352 GiB, the 7/2 case) included; under a 500 GB
    /// address space SIGSEGV again at 100,000,003, 166,000,001, 200,000,003
    /// and 1,000,000,003 arrays, where only `2^32 - 5` still throws - the same
    /// counts deciding the other way.
    /// The ceiling is therefore a documented line inside the range that every
    /// measured configuration reached, not a property of the source: between
    /// it and the address space's own limit the executed process still dies
    /// and the port records nothing, exactly as
    /// [`Limits::max_isotope_windows`] refuses before the source's allocation
    /// fails.
    pub(crate) fn score_arrays_overrun(&self) -> bool {
        let count = i64::from(self.charge_high) - i64::from(self.charge_low) + 1;
        if count == i64::from(i32::MAX) || count == -1 {
            return true;
        }
        if count > -4 {
            return false;
        }
        // `3 + 2 * count` modulo 2^32, as the source's `UInt` arithmetic.
        let arrays = u64::from((3u32).wrapping_add((count as u32).wrapping_mul(2)));
        arrays.saturating_mul(SOURCE_FLOAT_DATA_ARRAY_BYTES)
            <= SCORE_ARRAY_TERMINATION_CEILING_BYTES
    }

    /// The isotope count of the precalculated patterns: 20, plus 1000 for each
    /// changed abundance, as `run_` computes it.
    pub fn max_isotopes(&self) -> usize {
        BASE_MAX_ISOTOPES
            + OVERRIDE_EXTRA_ISOTOPES
                * (usize::from(self.abundance_12c_changed)
                    + usize::from(self.abundance_14n_changed))
    }
}

/// The `what()` text of the `Exception::ConversionError` that
/// `ParamValue::operator unsigned int` throws for a negative integer
/// (`ParamValue.cpp:456-459`).
pub const NEGATIVE_UNSIGNED_WHAT: &str =
    "Could not convert negative integer ParamValue to unsigned int";

/// The 64-bit integer behind an integer parameter (`ParamValue::data_.ssize_`).
fn source_i64(p: &Param, key: &str) -> Result<i64> {
    p.value(key)?.to_i64()
}

/// `ParamValue::operator unsigned int`: [`NEGATIVE_UNSIGNED_WHAT`] for a
/// negative value, and the low 32 bits of any other (`libOpenMS.so`
/// `0x7aee30`: `cvtsi2sd`/`comisd` against `0.0`, then the register's low
/// half).
fn source_unsigned(p: &Param, key: &str) -> Result<u32> {
    let value = source_i64(p, key)?;
    if (value as f64) < 0.0 {
        return Err(Error::InvalidValue(NEGATIVE_UNSIGNED_WHAT.into()));
    }
    Ok(value as u32)
}

/// The conversion at the start of source `run_` that can throw:
/// `UInt max_iterations = param_.getValue("fit:max_iterations")`
/// (`FeatureFinderAlgorithmPicked.cpp:152`).
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] with [`NEGATIVE_UNSIGNED_WHAT`] for a
/// negative value.
pub(crate) fn check_run_conversions(p: &Param) -> Result<()> {
    source_unsigned(p, "fit:max_iterations").map(|_| ())
}

/// Source `Param::checkDefaults(name, defaults, "")` (`Param.cpp:1066-1167`)
/// on a merged parameter set, with the source's messages.
///
/// Every entry is visited in iteration order. One the defaults do not know is
/// appended to `unknown` (the source logs a warning for it) and skipped. A
/// known entry of another value type than its default is refused with
/// `<name>: Wrong parameter type '<type>' for <type> parameter '<key>'
/// given!`. Otherwise the default entry with the given value is checked as
/// `Param::ParamEntry::isValid` checks it (`Param.cpp:57-167`), and a
/// violation is refused with `<name>: <message>`. That check narrows an
/// integer with `int tmp = value`, which keeps the low 32 bits of the 64-bit
/// value (C++20), so `2^32 + 10` passes a minimum of 1 as 10 does and `2^32`
/// fails it as 0 does (executed: `../oracle/ffap-complete-fix2`, case
/// `bigint`); the shared [`crate::param`] check refuses both instead.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] with the source's message at the first
/// refused entry; `unknown` then holds the unknown entries before it.
pub fn check_parameters(
    merged: &Param,
    defaults: &Param,
    name: &str,
    unknown: &mut Vec<String>,
) -> Result<()> {
    for item in merged.iter()? {
        if !defaults.exists(&item.key)? {
            unknown.push(item.key.clone());
            continue;
        }
        let default = defaults.entry(&item.key)?;
        let given = &item.entry.value;
        if default.value.value_type() != given.value_type() {
            return Err(Error::InvalidValue(format!(
                "{name}: Wrong parameter type '{}' for {} parameter '{}' given!",
                source_type_name(given.value_type()),
                source_type_name(default.value.value_type()),
                item.key
            )));
        }
        if let Some(message) = source_validity(default, given) {
            return Err(Error::InvalidValue(format!("{name}: {message}")));
        }
    }
    Ok(())
}

/// The type names of `Param::checkDefaults`.
fn source_type_name(value_type: ParamValueType) -> &'static str {
    match value_type {
        ParamValueType::String => "string",
        ParamValueType::StringList => "string list",
        ParamValueType::Empty => "empty",
        ParamValueType::Integer => "integer",
        ParamValueType::IntegerList => "integer list",
        ParamValueType::Float => "float",
        ParamValueType::FloatList => "float list",
    }
}

/// `std::to_string(double)`: `%f`, as glibc prints it.
fn std_to_string(value: f64) -> String {
    if value.is_nan() {
        if value.is_sign_negative() {
            "-nan"
        } else {
            "nan"
        }
        .to_owned()
    } else if value.is_infinite() {
        if value < 0.0 { "-inf" } else { "inf" }.to_owned()
    } else {
        format!("{value:.6}")
    }
}

/// `Param::ParamEntry::isValid` of `entry` holding `value`: the message of
/// its first violation, or `None`.
fn source_validity(entry: &crate::param::ParamEntry, value: &ParamValue) -> Option<String> {
    let has = |tag: &str| entry.tags.contains(tag);
    let valid_list = || entry.valid_strings.join(",");
    let int_range = |x: i32| {
        ((entry.min_int != -i32::MAX && x < entry.min_int)
            || (entry.max_int != i32::MAX && x > entry.max_int))
            .then(|| {
                format!(
                    "Invalid integer parameter value '{x}' for parameter '{}' given! The valid \
                     range is: [{}:{}].",
                    entry.name, entry.min_int, entry.max_int
                )
            })
    };
    let float_range = |x: f64| {
        ((entry.min_float != -f64::MAX && x < entry.min_float)
            || (entry.max_float != f64::MAX && x > entry.max_float))
            .then(|| {
                format!(
                    "Invalid double parameter value '{}' for parameter '{}' given! The valid \
                     range is: [{}:{}].",
                    std_to_string(x),
                    entry.name,
                    std_to_string(entry.min_float),
                    std_to_string(entry.max_float)
                )
            })
    };
    match value {
        ParamValue::String(text) => {
            let accepted = entry.valid_strings.is_empty()
                || entry.valid_strings.contains(text)
                || has("input file")
                || has("output file")
                || has("output prefix");
            (!accepted).then(|| {
                format!(
                    "Invalid string parameter value '{text}' for parameter '{}' given! Valid \
                     values are: '{}'.",
                    entry.name,
                    valid_list()
                )
            })
        }
        ParamValue::StringList(texts) => texts.iter().find_map(|text| {
            let accepted = entry.valid_strings.is_empty()
                || entry.valid_strings.contains(text)
                || has("input file")
                || has("output file");
            (!accepted).then(|| {
                format!(
                    "Invalid string parameter value '{text}' for parameter '{}' given! Valid \
                     values are: '{}'.",
                    entry.name,
                    valid_list()
                )
            })
        }),
        // `int tmp = value`: the low 32 bits.
        ParamValue::Integer(x) => int_range(*x as i32),
        ParamValue::IntegerList(xs) => xs.iter().find_map(|&x| int_range(x)),
        ParamValue::Float(x) => float_range(*x),
        ParamValue::FloatList(xs) => xs.iter().find_map(|&x| float_range(x)),
        ParamValue::Empty => None,
    }
}

/// The result of [`run`].
#[derive(Clone, Debug, PartialEq)]
pub struct RunOutput {
    /// The features found.
    pub features: FeatureMap,
    /// Warnings and progress lines the source writes to the OpenMS log or to
    /// `std::cout`, in the order they arise.
    pub log: Vec<String>,
    /// How often each abort reason kept a seed from becoming a feature: source
    /// member `aborts_`, which it logs at the end of `run_`.
    ///
    /// The source increments this `std::map` from inside its parallel region
    /// without synchronisation, which is a data race (B7 candidate 5 in
    /// `docs/FEATURE_FINDER_PICKED_SUPPORT.md`); this port aggregates the
    /// reasons serially, in seed order, so the counts are exact and
    /// thread-count independent.
    pub aborts: BTreeMap<String, usize>,
    /// The `write_debug` output, when the parameter is set
    /// ([`FeatureFinderAlgorithmPicked::debug_output`]).
    pub debug: Option<DebugOutput>,
}

/// Source warning when the input is not sorted, verbatim.
pub const UNSORTED_WARNING: &str =
    "Input map is not sorted by RT and m/z! This is done now, before applying the algorithm!";

/// The input checks at the start of source `run`, in source order.
///
/// 1. An experiment without spectra returns `Ok(false)`: the source clears the
///    output and returns.
/// 2. No peak in any spectrum or chromatogram (source `getSize() == 0`) is an
///    error.
/// 3. MS levels other than exactly `{1}` are an error.
/// 4. When the spectra are not sorted as source `MSExperiment::isSorted(true)`
///    finds them (no retention time greater than the next, every spectrum
///    `std::is_sorted` by m/z), the spectra and chromatograms are sorted and
///    [`UNSORTED_WARNING`] is logged.
/// 5. A non-empty spectrum whose first m/z is negative is an error; `-inf`
///    is negative, and the sort moves it to the front of its spectrum.
///
/// Returns `Ok(true)` when the run continues.
///
/// Infinite and NaN retention times, m/z values and intensities are not
/// refused here: the source reads them, and the port follows it (see
/// `docs/FEATURE_FINDER_PICKED_SUPPORT.md`, "Non-finite input"). The
/// comparisons of check 4 are false for a NaN, as in the source, so a NaN
/// never makes the input unsorted.
///
/// The sort of check 4 is the Release build's: `std::sort` of the spectra by
/// retention time and of the chromatograms by product m/z as libstdc++'s
/// introsort, then `std::stable_sort` of each unsorted spectrum's and
/// chromatogram's peaks as libstdc++'s merge sort
/// ([`crate::analysis::feature_finder_picked::source_sort`]), NaN keys and
/// equal keys included.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] with the source messages of
/// `Exception::IllegalArgument` for checks 2, 3 and 5, and with the `what()`
/// text of the source's `Exception::Precondition` when check 4 reaches an
/// unsorted spectrum (in retention-time order), or then an unsorted
/// chromatogram, whose non-empty data array differs in length from its peaks
/// (`FloatDataArray[0] size (25) does not match spectrum size (24)`). The
/// sorts of check 4 compare with `<` on `f64` keys, for which the introsort's
/// out-of-bounds guard is unreachable, NaN keys included (module documentation
/// of [`crate::analysis::feature_finder_picked::source_sort`]).
pub fn validate_input(experiment: &mut MSExperiment, log: &mut Vec<String>) -> Result<bool> {
    if experiment.spectra.is_empty() {
        return Ok(false);
    }
    if experiment.total_peak_count()? == 0 {
        return Err(Error::InvalidValue(
            "FeatureFinder needs updated ranges on input map. Aborting.".into(),
        ));
    }
    if experiment.ms_levels() != [1] {
        return Err(Error::InvalidValue(
            "FeatureFinder can only operate on MS level 1 data. Please do not use MS/MS data. \
             Aborting."
                .into(),
        ));
    }
    if !source_is_sorted(experiment) {
        log.push(UNSORTED_WARNING.to_string());
        source_sort_spectra(experiment)?;
        source_sort_chromatograms(experiment)?;
    }
    if experiment
        .spectra
        .iter()
        .any(|spectrum| spectrum.peaks.first().is_some_and(|peak| peak.mz < 0.0))
    {
        return Err(Error::InvalidValue(
            "FeatureFinder can only operate on spectra that contain peaks with positive m/z \
             values. Filter the data accordingly beforehand! Aborting."
                .into(),
        ));
    }
    Ok(true)
}

/// Source `MSExperiment::sortSpectra(true)`: `std::sort` of the spectra by
/// retention time (`SpectrumType::RTLess`), then `MSSpectrum::sortByPosition`
/// on each, in the sorted order, which returns when `std::is_sorted` holds and
/// otherwise sorts the peaks with `std::stable_sort` (`PositionLess`), keeping
/// the data arrays aligned.
///
/// A spectrum that is sorted this way and holds data arrays first runs
/// `MSSpectrum::checkDataArraySizes_` (`MSSpectrum.h`, `sort`), so the first
/// unsorted spectrum *in retention-time order* whose non-empty data array
/// differs in length from its peaks ends the run with the source's
/// `Exception::Precondition` text ([`data_array_sizes`]); the spectra before it
/// are sorted, which a caller cannot observe because `run` consumes the
/// experiment.
fn source_sort_spectra(experiment: &mut MSExperiment) -> Result<()> {
    source_sort_by(&mut experiment.spectra, |a, b| a.rt < b.rt)?;
    for spectrum in &mut experiment.spectra {
        if libstdcxx::is_sorted_by(&spectrum.peaks, |a, b| a.mz < b.mz) {
            continue;
        }
        data_array_sizes(
            &spectrum.float_data_arrays,
            &spectrum.string_data_arrays,
            &spectrum.integer_data_arrays,
            spectrum.peaks.len(),
            "spectrum",
        )?;
        let peaks = &spectrum.peaks;
        let order = source_stable_sort_permutation(
            peaks.len(),
            |a, b| peaks[a].mz < peaks[b].mz,
            TemporaryBuffer::Allocate,
        )?;
        spectrum.select(&order)?;
    }
    Ok(())
}

/// Source `MSSpectrum::checkDataArraySizes_` and
/// `MSChromatogram::checkDataArraySizes_`: the float, then the string, then the
/// integer data arrays, each in index order; the first non-empty array whose
/// length differs from `peaks` is an [`Error::InvalidValue`] with the `what()`
/// text of the source's `Exception::Precondition`, `<Kind>DataArray[<i>] size
/// (<n>) does not match <owner> size (<peaks>)` (executed:
/// `../oracle/ffap-complete-fix3`, `a_*` cases).
fn data_array_sizes<A, B, C>(
    floats: &[DataArray<A>],
    strings: &[DataArray<B>],
    integers: &[DataArray<C>],
    peaks: usize,
    owner: &str,
) -> Result<()> {
    fn check<T>(arrays: &[DataArray<T>], kind: &str, peaks: usize, owner: &str) -> Result<()> {
        for (index, array) in arrays.iter().enumerate() {
            if !array.data.is_empty() && array.data.len() != peaks {
                return Err(Error::InvalidValue(format!(
                    "{kind}DataArray[{index}] size ({}) does not match {owner} size ({peaks})",
                    array.data.len()
                )));
            }
        }
        Ok(())
    }
    check(floats, "Float", peaks, owner)?;
    check(strings, "String", peaks, owner)?;
    check(integers, "Integer", peaks, owner)
}

/// Source `MSExperiment::sortChromatograms(true)`: `std::sort` of the
/// chromatograms by product m/z (`ChromatogramType::MZLess`), then
/// `MSChromatogram::sortByPosition` on each, which returns when no retention
/// time exceeds the next (false for a NaN) and otherwise sorts the peaks with
/// `std::stable_sort`, after the same data-array check as the spectra's
/// ([`data_array_sizes`], "chromatogram size").
fn source_sort_chromatograms(experiment: &mut MSExperiment) -> Result<()> {
    let is_sorted = |chromatogram: &MSChromatogram| {
        !chromatogram
            .peaks
            .windows(2)
            .any(|pair| pair[0].rt > pair[1].rt)
    };
    source_sort_by(&mut experiment.chromatograms, |a, b| {
        a.product.mz < b.product.mz
    })?;
    for chromatogram in &mut experiment.chromatograms {
        if is_sorted(chromatogram) {
            continue;
        }
        data_array_sizes(
            &chromatogram.float_data_arrays,
            &chromatogram.string_data_arrays,
            &chromatogram.integer_data_arrays,
            chromatogram.peaks.len(),
            "chromatogram",
        )?;
        let peaks = &chromatogram.peaks;
        let order = source_stable_sort_permutation(
            peaks.len(),
            |a, b| peaks[a].rt < peaks[b].rt,
            TemporaryBuffer::Allocate,
        )?;
        chromatogram.select(&order)?;
    }
    Ok(())
}

/// Find features: source `run(PeakMap&&, FeatureMap&, const Param&, const FeatureMap& seeds)`
/// on a fresh algorithm object and a fresh feature map.
///
/// `experiment` holds centroided MS1 spectra and is consumed, as the source
/// moves it. `parameters` are applied over [`default_parameters`]. `seeds`
/// holds user-specified seeds; an empty map lets the algorithm find seeds
/// itself.
///
/// An experiment without spectra yields an empty feature map and an empty log.
/// Otherwise the input is checked ([`validate_input`]), the parameters are
/// applied ([`Settings::from_parameters`]), the seeds are selected and extended
/// charge by charge, and the overlaps are resolved
/// ([`FeatureFinderAlgorithmPicked::run`]). The log holds the text of every
/// console line in source order; `write_debug` output is in
/// [`RunOutput::debug`].
///
/// A reused object or a non-empty output map, whose state the source carries
/// from run to run, is [`FeatureFinderAlgorithmPicked`] itself.
///
/// # Errors
///
/// Every error of [`FeatureFinderAlgorithmPicked::run`].
pub fn run(experiment: MSExperiment, seeds: &FeatureMap, parameters: &Param) -> Result<RunOutput> {
    run_with_options(experiment, seeds, parameters, &Options::default())
}

/// [`run`] with explicit [`Options`].
///
/// # Errors
///
/// As [`run`].
pub fn run_with_options(
    experiment: MSExperiment,
    seeds: &FeatureMap,
    parameters: &Param,
    options: &Options,
) -> Result<RunOutput> {
    let mut algorithm = FeatureFinderAlgorithmPicked::with_options(*options)?;
    let mut features = FeatureMap::new();
    algorithm.run(experiment, &mut features, parameters, seeds)?;
    Ok(RunOutput {
        features,
        log: algorithm
            .report()
            .iter()
            .filter_map(|line| line.text().map(str::to_owned))
            .collect(),
        aborts: algorithm
            .aborts()
            .iter()
            .map(|(reason, &count)| (reason.clone(), count as usize))
            .collect(),
        debug: algorithm.take_debug_output(),
    })
}

/// Steps 3.3 and 4 of source `run_` on a completed seed stage.
///
/// Per charge, every seed is extended in parallel
/// ([`crate::concept::parallel::map_collect`] with
/// [`Options::threads`], the port's form of the source's single
/// `#pragma omp parallel for`), then a serial pass in seed order drops the
/// candidates whose seed lies inside an earlier, more intense feature and
/// numbers the survivors. After every charge the overlapping features are
/// resolved ([`resolve_overlaps`]), the zero-intensity losers are removed, the
/// map is sorted by descending intensity and each feature is annotated with its
/// apex scan ([`annotate_apex`]). Both sorts put equal elements where the C++
/// Release build's `std::sort` does
/// ([`crate::analysis::feature_finder_picked::source_sort`]).
///
/// The log receives, in the source's order, the seed counts the stage already
/// collected, each charge's `Found N feature candidates for charge c.` directly
/// after its seed line, the overlap count, the apex warning if any, the abort
/// reasons and the feature count. Two entries are empty strings, the blank
/// lines the source prints before the abort block and before the feature count
/// (`FeatureFinderAlgorithmPicked.cpp:1019` and `1026`).
///
/// Because the parallel results keep the input order and every later step is
/// serial, the output does not depend on [`Options::threads`]; the source's
/// results are schedule-independent for the same reason, its `tmp_feature_map`
/// being keyed by seed index.
///
/// This stage-level entry point writes no debug output and reports no
/// progress; [`FeatureFinderAlgorithmPicked::run`] does both.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when a [`Limits`] ceiling of the seed loop
/// is exceeded, checked before the loop starts, and every error of the
/// extension, the fit, the checks and the feature creation, for the first seed
/// in seed order that fails: a fit error is one of the port's ceilings or a
/// NaN retention time in a mass trace, where the source's intensity profile
/// never returns; a feature m/z without isotope window is where the source's
/// exception terminates the process (`FittedModel::fit`,
/// [`build_feature`](crate::analysis::feature_finder_picked::fitting::build_feature)).
pub fn feature_stage(stage: &SeedStage, options: &Options) -> Result<RunOutput> {
    let settings = stage.settings();
    let experiment = stage.experiment();
    let mut seed_work = 0u64;
    for charge_index in 0..stage.charges().len() {
        preflight_charge(stage, charge_index, &options.limits, &mut seed_work)?;
    }

    let fitter_parameters = TraceFitterParams {
        max_iteration: i64::from(settings.max_iterations),
        weighted: false,
    };
    let mut features: Vec<Feature> = Vec::new();
    let mut aborts: BTreeMap<String, u32> = BTreeMap::new();
    let mut abort_reasons = AbortReasons::new();
    let mut candidate_lines: Vec<(usize, String)> = Vec::new();
    let mut plot_nr_global: i64 = -1;
    let mut feature_nr_global: i64 = 0;
    let mut progress = Progress::silent();
    let no_parameters = Param::new();

    for (charge_index, charge_seeds) in stage.charges().iter().enumerate() {
        let outcomes = extend_charge(
            stage,
            charge_index,
            &fitter_parameters,
            options.threads,
            false,
        );
        let mut book = Bookkeeping {
            aborts: &mut aborts,
            abort_reasons: &mut abort_reasons,
            out: &mut None,
            termination: &mut None,
            features: &mut features,
            plot_nr_global: &mut plot_nr_global,
            feature_nr_global: &mut feature_nr_global,
            progress: &mut progress,
            limits: &options.limits,
        };
        let key = DebugKey {
            policy: options.pseudo_rt_shift,
            parameters: &no_parameters,
        };
        let feature_candidates = settle_charge(stage, charge_index, outcomes, &mut book, &key)?;
        candidate_lines.push((
            charge_index,
            format!(
                "Found {feature_candidates} feature candidates for charge {}.",
                charge_seeds.charge
            ),
        ));
    }

    // Step 4, serial.
    let mut map = FeatureMap::from_features(features);
    source_sort_by(&mut map.features, |a, b| a.mz < b.mz)?;
    let removed = resolve_overlaps(&mut map.features, settings.max_feature_intersection)?;
    map.features.retain(|feature| feature.intensity != 0.0);
    source_sort_by(&mut map.features, |a, b| b.intensity < a.intensity)?;
    let invalid_apex = annotate_apex(&mut map.features, experiment)?;

    let mut log = interleave_log(stage, &candidate_lines);
    log.push(format!("Removed {removed} overlapping features."));
    if invalid_apex > 0 {
        log.push(invalid_apex_warning(invalid_apex));
    }
    // The source separates the abort block and the feature count from what
    // precedes them with a blank line each (`OPENMS_LOG_INFO << '\n'` at
    // FeatureFinderAlgorithmPicked.cpp:1019 and the leading "\n" at 1026), which
    // the executed C++ prints whether or not there is an abort reason.
    log.push(String::new());
    log.push(ABORT_BLOCK_HEADING.into());
    for (reason, count) in &aborts {
        log.push(format!(" - {reason}: {count} times"));
    }
    log.push(String::new());
    log.push(format!("{} features found.", map.len()));
    Ok(RunOutput {
        features: map,
        log,
        aborts: aborts
            .into_iter()
            .map(|(reason, count)| (reason, count as usize))
            .collect(),
        debug: None,
    })
}

/// The stage's log with each charge's candidate line after its seed line.
///
/// The source interleaves the two `std::cout` lines because it runs steps 3.1
/// to 3.3 charge by charge; this port computes every charge's seeds first (see
/// [`crate::analysis::feature_finder_picked::seeds`]) and restores the order
/// here. A charge whose seed line is missing, which cannot happen, appends its
/// candidate line at the end.
fn interleave_log(stage: &SeedStage, candidate_lines: &[(usize, String)]) -> Vec<String> {
    let mut log: Vec<String> = Vec::new();
    let mut seen = 0usize;
    for line in stage.log() {
        log.push(line.clone());
        if line.starts_with("Found ") && line.contains(" seeds for charge ") {
            if let Some((_, candidate)) = candidate_lines.iter().find(|(c, _)| *c == seen) {
                log.push(candidate.clone());
            }
            seen += 1;
        }
    }
    for (charge_index, line) in candidate_lines {
        if *charge_index >= seen {
            log.push(line.clone());
        }
    }
    log
}
