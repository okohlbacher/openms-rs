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
//!   `#pragma omp parallel for`.
//!
//! The API mapping, the preserved source conventions, the native differences
//! and the evidence are in `docs/FEATURE_FINDER_PICKED_SUPPORT.md`.
//!
//! The source writes its seed counts to `std::cout` and its warnings to the
//! OpenMS log. Library code here never prints: every such line is collected in
//! the log of the returned value.

use std::collections::{BTreeMap, BTreeSet};

use crate::analysis::feature_finder_picked::extension::{
    OverallScores, extend_mass_traces, find_best_isotope_fit,
};
use crate::analysis::feature_finder_picked::fitting::{
    ABORT_COULD_NOT_EXTEND, ABORT_NO_ISOTOPE_PATTERN, FeatureInput, FittedModel, QualityOutcome,
    build_feature, check_feature_quality, crop_feature,
};
use crate::analysis::feature_finder_picked::helper_structs::Seed;
use crate::analysis::feature_finder_picked::resolution::{
    annotate_apex, invalid_apex_warning, resolve_overlaps,
};
use crate::analysis::feature_finder_picked::scoring::{libstdcxx, source_is_sorted};
use crate::analysis::feature_finder_picked::seeds::SeedStage;
use crate::analysis::feature_finder_picked::trace_fitter::TraceFitterParams;
use crate::concept::parallel::{Threads, map_collect};
use crate::kernel::{Feature, FeatureMap, MSChromatogram, MSExperiment, MSSpectrum, Point2D};
use crate::metadata::MetaValue;
use crate::param::{DefaultParamHandler, Param, ParamValue};
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
    /// `Exception::InvalidParameter`), or when a value cannot be converted.
    pub fn from_parameters(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(HANDLER_NAME)?;
        let defaults = default_parameters()?;
        handler.set_defaults(defaults.clone())?;
        handler.defaults_to_parameters()?;
        handler.set_parameters_with(parameters, |merged| Self::read(merged, &defaults))
    }

    fn read(p: &Param, defaults: &Param) -> Result<Self> {
        let float = |key: &str| p.value(key)?.to_f64();
        let int = |key: &str| p.value(key)?.to_i32();
        let unsigned = |key: &str| p.value(key)?.to_u32();
        let min_spectra = (f64::from(int("mass_trace:min_spectra")?) * 0.5).floor();
        let reported_mz = match p.value("feature:reported_mz")?.as_str()? {
            "maximum" => ReportedMz::Maximum,
            "average" => ReportedMz::Average,
            _ => ReportedMz::Monoisotopic,
        };
        // Source: `param_.getValue("feature:rt_shape") == "asymmetric"`, else symmetric.
        let rt_shape = if p.value("feature:rt_shape")? == &ParamValue::String("asymmetric".into()) {
            RtShape::Asymmetric
        } else {
            RtShape::Symmetric
        };
        let changed = |key: &str| -> Result<bool> { Ok(p.value(key)? != defaults.value(key)?) };
        Ok(Self {
            pattern_tolerance: float("isotopic_pattern:mz_tolerance")?,
            trace_tolerance: float("mass_trace:mz_tolerance")?,
            min_spectra: min_spectra as usize,
            max_missing_trace_peaks: unsigned("mass_trace:max_missing")?,
            slope_bound: float("mass_trace:slope_bound")?,
            intensity_percentage: float("isotopic_pattern:intensity_percentage")? / 100.0,
            intensity_percentage_optional: float("isotopic_pattern:intensity_percentage_optional")?
                / 100.0,
            optional_fit_improvement: float("isotopic_pattern:optional_fit_improvement")? / 100.0,
            mass_window_width: float("isotopic_pattern:mass_window_width")?,
            intensity_bins: unsigned("intensity:bins")? as usize,
            min_isotope_fit: float("feature:min_isotope_fit")?,
            min_trace_score: float("feature:min_trace_score")?,
            min_rt_span: float("feature:min_rt_span")?,
            max_rt_span: float("feature:max_rt_span")?,
            max_feature_intersection: float("feature:max_intersection")?,
            reported_mz,
            min_feature_score: float("feature:min_score")?,
            charge_low: int("isotopic_pattern:charge_low")?,
            charge_high: int("isotopic_pattern:charge_high")?,
            max_iterations: unsigned("fit:max_iterations")?,
            abundance_12c: float("isotopic_pattern:abundance_12C")?,
            abundance_14n: float("isotopic_pattern:abundance_14N")?,
            abundance_12c_changed: changed("isotopic_pattern:abundance_12C")?,
            abundance_14n_changed: changed("isotopic_pattern:abundance_14N")?,
            seed_min_score: float("seed:min_score")?,
            user_seed_rt_tolerance: float("user-seed:rt_tolerance")?,
            user_seed_mz_tolerance: float("user-seed:mz_tolerance")?,
            user_seed_min_score: float("user-seed:min_score")?,
            write_debug: p.value("write_debug")?.to_bool()?,
            rt_shape,
        })
    }

    /// The number of charges searched, `charge_high - charge_low + 1`, or zero
    /// when `charge_low` exceeds `charge_high` by one.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `charge_low` exceeds `charge_high`
    /// by more than one. The source computes the count in `UInt`, which wraps:
    /// a difference of one gives zero charges and an empty result, and a larger
    /// difference resizes the score arrays to a wrapped size and then indexes
    /// past their end, which is undefined behaviour.
    pub fn charge_count(&self) -> Result<usize> {
        let count = i64::from(self.charge_high) - i64::from(self.charge_low) + 1;
        if count < 0 {
            return Err(Error::InvalidValue(format!(
                "isotopic_pattern:charge_low {} exceeds charge_high {} by more than one; the source \
                 behaviour is undefined",
                self.charge_low, self.charge_high
            )));
        }
        usize::try_from(count).map_err(|_| Error::InvalidValue("charge count overflow".into()))
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
/// # Errors
///
/// Returns [`Error::InvalidValue`] with the source messages of
/// `Exception::IllegalArgument` for checks 2, 3 and 5. The sort of check 4
/// returns [`Error::InvalidValue`] where the source's sort is undefined or
/// leaves an observable order unspecified:
///
/// - a NaN retention time when the spectra are sorted
///   (`MSExperiment::sortSpectra`, `std::sort` by retention time, `MSExperiment.cpp:793`), and
///   a NaN chromatogram product m/z when the chromatograms are
///   (`sortChromatograms`, `std::sort`, `MSExperiment.cpp:813`): either the comparator is not a
///   strict weak ordering, or every key is equivalent and the order of the
///   spectra or chromatograms, which the output shows, is libstdc++'s
///   introsort order, which this module does not reproduce;
/// - a NaN m/z in a spectrum, or a NaN retention time in a chromatogram, that
///   `std::is_sorted` or the hand-written check finds unsorted, and which
///   `std::stable_sort` then orders under a comparator that is not a strict
///   weak ordering.
///
/// The kernel's own chromatogram checks apply to the chromatograms' other
/// fields.
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

/// A stable permutation of `0..len` by `key`, which must hold no NaN.
fn stable_order(len: usize, key: impl Fn(usize) -> f64) -> Vec<usize> {
    let mut order: Vec<usize> = (0..len).collect();
    order.sort_by(|&a, &b| {
        key(a)
            .partial_cmp(&key(b))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    order
}

/// The refusal of a source sort whose keys hold a NaN.
fn nan_sort_refusal(what: &str, call: &str) -> Error {
    Error::InvalidValue(format!(
        "FeatureFinderAlgorithmPicked input: {what} is NaN and the input is not sorted; the \
         source sorts it with {call}, whose result the standard leaves undefined or unspecified \
         here and which the port does not reproduce"
    ))
}

/// Source `MSExperiment::sortSpectra(true)`: `std::sort` of the spectra by
/// retention time, then `MSSpectrum::sortByPosition` on each, which returns
/// when `std::is_sorted` holds and otherwise sorts stably, keeping the data
/// arrays aligned. Every refusal is checked before anything moves.
///
/// Spectra with equal retention times keep their order here; the source's
/// `std::sort` is not stable.
fn source_sort_spectra(experiment: &mut MSExperiment) -> Result<()> {
    let spectra = &experiment.spectra;
    if spectra.len() > 1 && spectra.iter().any(|spectrum| spectrum.rt.is_nan()) {
        return Err(nan_sort_refusal(
            "a retention time",
            "std::sort (MSExperiment::sortSpectra)",
        ));
    }
    let mut unsorted = Vec::new();
    for (index, spectrum) in spectra.iter().enumerate() {
        if libstdcxx::is_sorted_by(&spectrum.peaks, |a, b| a.mz < b.mz) {
            continue;
        }
        if spectrum.peaks.iter().any(|peak| peak.mz.is_nan()) {
            return Err(nan_sort_refusal(
                "an m/z",
                "std::stable_sort (MSSpectrum::sortByPosition)",
            ));
        }
        unsorted.push(index);
    }
    // Validate the data arrays of every spectrum that moves before moving any.
    for &index in &unsorted {
        let spectrum = &mut experiment.spectra[index];
        let identity: Vec<usize> = (0..spectrum.peaks.len()).collect();
        spectrum.select(&identity)?;
    }
    for &index in &unsorted {
        let spectrum = &mut experiment.spectra[index];
        let order = stable_order(spectrum.peaks.len(), |i| spectrum.peaks[i].mz);
        spectrum.select(&order)?;
    }
    let order = stable_order(experiment.spectra.len(), |i| experiment.spectra[i].rt);
    let mut slots: Vec<Option<MSSpectrum>> = std::mem::take(&mut experiment.spectra)
        .into_iter()
        .map(Some)
        .collect();
    experiment.spectra = order.iter().filter_map(|&i| slots[i].take()).collect();
    Ok(())
}

/// Source `MSExperiment::sortChromatograms(true)`: `std::sort` of the
/// chromatograms by product m/z, then `MSChromatogram::sortByPosition` on
/// each, which returns when no retention time exceeds the next and otherwise
/// sorts stably. Every refusal is checked before anything moves.
///
/// Chromatograms with equal product m/z keep their order here; the source's
/// `std::sort` is not stable.
fn source_sort_chromatograms(experiment: &mut MSExperiment) -> Result<()> {
    let chromatograms = &experiment.chromatograms;
    if chromatograms.len() > 1
        && chromatograms
            .iter()
            .any(|chromatogram| chromatogram.product.mz.is_nan())
    {
        return Err(nan_sort_refusal(
            "a chromatogram's product m/z",
            "std::sort (MSExperiment::sortChromatograms)",
        ));
    }
    let mut unsorted = Vec::new();
    for (index, chromatogram) in chromatograms.iter().enumerate() {
        // Source `MSChromatogram::isSorted`: no retention time greater than the
        // next, false for a NaN.
        if !chromatogram
            .peaks
            .windows(2)
            .any(|pair| pair[0].rt > pair[1].rt)
        {
            continue;
        }
        if chromatogram.peaks.iter().any(|peak| peak.rt.is_nan()) {
            return Err(nan_sort_refusal(
                "a chromatogram retention time",
                "std::stable_sort (MSChromatogram::sortByPosition)",
            ));
        }
        unsorted.push(index);
    }
    for &index in &unsorted {
        let chromatogram = &mut experiment.chromatograms[index];
        let identity: Vec<usize> = (0..chromatogram.peaks.len()).collect();
        chromatogram.select(&identity)?;
    }
    for &index in &unsorted {
        let chromatogram = &mut experiment.chromatograms[index];
        let order = stable_order(chromatogram.peaks.len(), |i| chromatogram.peaks[i].rt);
        chromatogram.select(&order)?;
    }
    let order = stable_order(experiment.chromatograms.len(), |i| {
        experiment.chromatograms[i].product.mz
    });
    let mut slots: Vec<Option<MSChromatogram>> = std::mem::take(&mut experiment.chromatograms)
        .into_iter()
        .map(Some)
        .collect();
    experiment.chromatograms = order.iter().filter_map(|&i| slots[i].take()).collect();
    Ok(())
}

/// Find features: source `run(PeakMap&&, FeatureMap&, const Param&, const FeatureMap& seeds)`.
///
/// `experiment` holds centroided MS1 spectra and is consumed, as the source
/// moves it. `parameters` are applied over [`default_parameters`]. `seeds`
/// holds user-specified seeds; an empty map lets the algorithm find seeds
/// itself.
///
/// An experiment without spectra yields an empty feature map and an empty log.
/// Otherwise the input is checked ([`validate_input`]), the parameters are
/// applied ([`Settings::from_parameters`]), the seed stage runs
/// ([`SeedStage::compute`]) and every seed is extended into a feature
/// ([`feature_stage`]).
///
/// # Errors
///
/// Every error of [`validate_input`], [`Settings::from_parameters`],
/// [`SeedStage::compute`] and [`feature_stage`].
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
    match SeedStage::run_with_options(experiment, seeds, parameters, options)? {
        None => Ok(RunOutput {
            features: FeatureMap::new(),
            log: Vec::new(),
            aborts: BTreeMap::new(),
        }),
        Some(stage) => feature_stage(&stage, options),
    }
}

/// What one seed produced in step 3.3.
struct SeedOutcome {
    /// Whether the seed reached the fit and therefore consumed a `plot_nr`.
    plot_nr_used: bool,
    /// The candidate, or the source abort reason that dropped the seed.
    result: std::result::Result<SeedCandidate, String>,
}

/// One accepted candidate and the later seeds it swallows.
struct SeedCandidate {
    feature: Feature,
    /// Indices of the seeds after this one that lie inside the feature: the
    /// source's `seeds_in_features[i]`.
    contained: Vec<usize>,
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
/// apex scan ([`annotate_apex`]).
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
/// # Errors
///
/// Returns [`Error::InvalidValue`] when a [`Limits`] ceiling of the seed loop
/// is exceeded, checked before the loop starts, and every error of the
/// extension, the fit, the checks and the feature creation, for the first seed
/// in seed order that fails: a fit error is one of the port's ceilings or a
/// NaN retention time in a mass trace, where the source's intensity profile
/// never returns; a feature m/z without isotope window is where the source's
/// exception terminates the process (`FittedModel::fit`,
/// [`build_feature`]).
pub fn feature_stage(stage: &SeedStage, options: &Options) -> Result<RunOutput> {
    let settings = stage.settings();
    let experiment = stage.experiment();
    preflight_seed_loop(stage, &options.limits)?;

    let fitter_parameters = TraceFitterParams {
        max_iteration: i64::from(settings.max_iterations),
        weighted: false,
    };
    let mut features: Vec<Feature> = Vec::new();
    let mut aborts: BTreeMap<String, usize> = BTreeMap::new();
    let mut candidate_lines: Vec<(usize, String)> = Vec::new();
    let mut plot_nr_global: i64 = -1;
    let mut feature_nr_global: i64 = 0;

    for (charge_index, charge_seeds) in stage.charges().iter().enumerate() {
        let charge = charge_seeds.charge;
        let seeds = &charge_seeds.seeds;
        let indices: Vec<usize> = (0..seeds.len()).collect();
        let overall = OverallScores::new(stage.scores(), charge_index);
        let outcomes = map_collect(&indices, options.threads, |&index| {
            extend_seed(stage, overall, &fitter_parameters, charge, seeds, index)
        });

        let mut accepted: Vec<(usize, SeedCandidate)> = Vec::new();
        for (index, outcome) in outcomes.into_iter().enumerate() {
            let outcome = outcome?;
            let plot_nr = if outcome.plot_nr_used {
                plot_nr_global += 1;
                plot_nr_global
            } else {
                -1
            };
            match outcome.result {
                Err(reason) => *aborts.entry(reason).or_insert(0) += 1,
                Ok(mut candidate) => {
                    // The source assigns `plot_nr` inside a critical section,
                    // so its value depends on the schedule; it is overwritten
                    // below for every candidate that survives, and only the
                    // refused debug output reads it otherwise. This port
                    // numbers the seeds that reached the fit in seed order.
                    candidate
                        .feature
                        .metadata
                        .insert("label".into(), MetaValue::from(plot_nr));
                    accepted.push((index, candidate));
                }
            }
        }

        let mut contained_seeds: BTreeSet<usize> = BTreeSet::new();
        let mut feature_candidates = 0usize;
        for (seed_nr, candidate) in accepted {
            if contained_seeds.contains(&seed_nr) {
                continue;
            }
            feature_candidates += 1;
            let mut feature = candidate.feature;
            feature
                .metadata
                .insert("label".into(), MetaValue::from(feature_nr_global));
            feature_nr_global += 1;
            features
                .try_reserve(1)
                .map_err(|_| Error::InvalidValue("cannot allocate a feature".into()))?;
            features.push(feature);
            contained_seeds.extend(candidate.contained);
        }
        candidate_lines.push((
            charge_index,
            format!("Found {feature_candidates} feature candidates for charge {charge}."),
        ));
    }

    // Step 4, serial.
    let mut map = FeatureMap::from_features(features);
    map.sort_by_mz()?;
    let removed = resolve_overlaps(&mut map.features, settings.max_feature_intersection)?;
    map.features.retain(|feature| feature.intensity != 0.0);
    map.sort_by_intensity(true)?;
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
    log.push("Info: reasons for not finalizing a feature during its construction:".into());
    for (reason, count) in &aborts {
        log.push(format!(" - {reason}: {count} times"));
    }
    log.push(String::new());
    log.push(format!("{} features found.", map.len()));
    Ok(RunOutput {
        features: map,
        log,
        aborts,
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

/// The ceilings of the seed loop, checked before it starts.
fn preflight_seed_loop(stage: &SeedStage, limits: &Limits) -> Result<()> {
    let spectra = stage.experiment().spectra.len() as u64;
    let isotopes = stage.settings().max_isotopes() as u64;
    let per_seed = isotopes.saturating_mul(isotopes.saturating_add(spectra));
    let mut work = 0u64;
    for charge in stage.charges() {
        if charge.seeds.len() > limits.max_seeds {
            return Err(Error::InvalidValue(format!(
                "{} seeds for charge {} exceed the limit of {}",
                charge.seeds.len(),
                charge.charge,
                limits.max_seeds
            )));
        }
        work = work.saturating_add((charge.seeds.len() as u64).saturating_mul(per_seed));
    }
    if work > limits.max_seed_work {
        return Err(Error::InvalidValue(format!(
            "the seed loop may take {work} work units, exceeding the limit of {}",
            limits.max_seed_work
        )));
    }
    Ok(())
}

/// One seed of step 3.3: isotope fit, extension, fit, cropping, quality checks
/// and feature creation.
fn extend_seed(
    stage: &SeedStage,
    overall: OverallScores<'_>,
    fitter_parameters: &TraceFitterParams,
    charge: i32,
    seeds: &[Seed],
    index: usize,
) -> Result<SeedOutcome> {
    let settings = stage.settings();
    let spectra = &stage.experiment().spectra;
    let seed = seeds[index];
    let aborted = |plot_nr_used: bool, reason: &str| SeedOutcome {
        plot_nr_used,
        result: Err(reason.to_owned()),
    };

    let (isotope_fit_quality, pattern) =
        find_best_isotope_fit(spectra, stage.windows(), settings, seed, charge)?;
    if isotope_fit_quality < settings.min_isotope_fit {
        return Ok(aborted(false, ABORT_NO_ISOTOPE_PATTERN));
    }
    let mut traces = extend_mass_traces(spectra, overall, settings, &pattern)?;
    let seed_mz = spectra[seed.spectrum].peaks[seed.peak].mz;
    if !traces.is_valid(seed_mz, settings.trace_tolerance) {
        return Ok(aborted(false, ABORT_COULD_NOT_EXTEND));
    }

    // Source: the baseline estimate is three quarters of the lowest peak.
    traces.update_baseline();
    traces.baseline *= 0.75;
    traces
        .get_mut(traces.max_trace)
        .ok_or_else(|| {
            Error::InvalidValue(
                "FeatureFinderAlgorithmPicked seed extension: the maximum trace index is out of \
                 range; the source dereferences it here"
                    .into(),
            )
        })?
        .update_maximum();

    let mut model = FittedModel::new(settings.rt_shape, *fitter_parameters);
    // The source's fit can throw `Exception::UnableToFit` (`TraceFitter.cpp:111`,
    // `:129`), which would escape its parallel region and end the process, but
    // no input reaches either throw from here (see `FittedModel::fit`). An error
    // here is therefore one of the port's own ceilings or the refused merge of a
    // NaN retention time into the intensity profile, where the source never
    // returns; the run fails with it rather than turning it into an abort
    // reason the source never records.
    model.fit(&traces)?;
    let new_traces = crop_feature(model.as_fitter(), &traces, settings.min_trace_score)?;
    let quality = match check_feature_quality(model.as_fitter(), &new_traces, seed_mz, settings)? {
        QualityOutcome::Rejected(reason) => return Ok(aborted(true, reason)),
        QualityOutcome::Accepted(quality) => quality,
    };
    let traces = new_traces;
    let feature = build_feature(FeatureInput {
        model: &model,
        traces: &traces,
        pattern: &pattern,
        windows: stage.windows(),
        settings,
        charge,
        // Overwritten serially; see `feature_stage`.
        plot_nr: -1,
        quality,
    })?;

    // Source: every later seed inside both the overall bounding box and one of
    // the mass-trace hulls.
    let mut contained = Vec::new();
    if let Some(bounds) = feature.hull_bounding_box() {
        for (offset, later) in seeds.iter().enumerate().skip(index + 1) {
            let rt = spectra[later.spectrum].rt;
            let mz = spectra[later.spectrum].peaks[later.peak].mz;
            if bounds.encloses(Point2D::new(rt, mz))? && feature.encloses(rt, mz)? {
                contained.push(offset);
            }
        }
    }
    Ok(SeedOutcome {
        plot_nr_used: true,
        result: Ok(SeedCandidate { feature, contained }),
    })
}
