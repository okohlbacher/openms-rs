// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The Biosaur2 feature detection workflow for LC-MS1 data, the source
//! `FEATUREFINDER/Biosaur2Algorithm.h` (core `bc9cc12`).
//!
//! The source documentation: *This algorithm is a C++ reimplementation of the
//! Biosaur2 feature detection method, originally developed in Python. It
//! processes centroided LC-MS1 data (with optional profile mode support) to
//! detect peptide features by*
//!
//! 1. *building retention time-contiguous signal traces called "hills" by
//!    linking peaks across consecutive scans;*
//! 2. *splitting hills at valley points to separate co-eluting species;*
//! 3. *detecting and annotating isotopic patterns based on expected mass
//!    differences and intensity correlations;*
//! 4. *calculating comprehensive feature properties including m/z, RT,
//!    intensity, and charge state.*
//!
//! *The algorithm supports several advanced features:* FAIMS compensation
//! voltage grouping, ion-mobility-aware processing for PASEF/TIMS data,
//! TOF-specific intensity filtering, automatic mass calibration and profile
//! mode centroiding via [`PeakPickerHiRes`](crate::processing::peak_picking::PeakPickerHiRes).
//! *The implementation closely mirrors the reference Python implementation to
//! ensure reproducible results.* The TOPP wrapper is `FeatureFinderLFQ`.
//!
//! Reference: Abdrakhimov *et al.*, *Biosaur: An open-source Python software
//! for liquid chromatography-mass spectrometry peptide feature detection with
//! ion mobility support*, Rapid Communications in Mass Spectrometry, 2022,
//! <https://doi.org/10.1002/rcm.9045>.
//!
//! The port is split by pipeline stage: [`preprocess`] (MS1 selection, profile
//! centroiding, TOF filtering, PASEF/TIMS centroiding), [`hills`] (hill
//! linking, filtering and valley splitting), [`isotopes`] (isotope pattern
//! assembly) and [`output`] (the `FeatureMap` conversion and the two TSV
//! writers). `docs/BIOSAUR2_SUPPORT.md` maps every source member to its Rust
//! counterpart and records the evidence.
//!
//! # Fidelity
//!
//! Every number [`Biosaur2Algorithm::run_with_details`] computes is compared
//! bit for bit with the Linux x86_64 **Release** build
//! `openms4-release-bc9cc12-c19e494-174b576` (`tests/biosaur2.rs`). The three
//! `std::sort` calls follow that build's `libstdc++` introsort
//! ([`source_sort_by`](crate::math::source_sort::source_sort_by)), so equal
//! intensities, equal hill m/z and equal pattern ranks are ordered as the C++
//! orders them; `static_cast<int>` of a `double` follows `cvttsd2si`; the
//! averagine table is the Release build's Boost.Math binomial table
//! ([`isotopes`]); and the `omp simd` reduction of `cosineCorrelation1D_` is
//! summed in the order the Release build's vectorised loop sums it.

mod averagine_table;
pub mod hills;
pub mod isotopes;
pub mod output;
pub mod preprocess;

use crate::concept::parallel::{Threads, map_collect};
use crate::concept::unique_id::UniqueIdGenerator;
use crate::format::file_info::text_format::ostream_g;
use crate::kernel::faims_helper::{CompensationVoltage, FaimsHelper};
use crate::kernel::im_data_converter::{FaimsSplitLogLevel, ImDataConverter};
use crate::kernel::{FeatureMap, MSExperiment, MSSpectrum};
use crate::math::source_sort::source_sort_by;
use crate::metadata::{DriftTimeUnit, ImTypes, IonMobilityFormat};
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::processing::feature_overlap_filter::FeatureOverlapFilter;
use crate::{Error, Result};

/// The `DefaultParamHandler` name of the source constructor.
pub const BIOSAUR2_ALGORITHM_NAME: &str = "Biosaur2Algorithm";

/// Representation of a single hill (continuous m/z trace across adjacent
/// scans), the source `Biosaur2Algorithm::Hill`.
///
/// The source documentation: *A hill represents a contiguous series of peaks at
/// similar m/z values across multiple consecutive scans. Hills are the
/// fundamental building blocks for feature detection. Each hill stores the raw
/// peak data (indices, m/z, intensity, RT) along with optional ion mobility
/// information for PASEF/TIMS data. Summary statistics (intensity-weighted mean
/// m/z, apex RT/intensity, etc.) are precomputed for efficient downstream
/// processing.*
///
/// The seven per-point vectors run in parallel, one entry per linked peak.
#[derive(Clone, Debug, PartialEq)]
pub struct Hill {
    /// Indices of spectra containing peaks of this hill, into the MS1 spectra
    /// of the processed group.
    pub scan_indices: Vec<usize>,
    /// Indices of peaks within each spectrum.
    pub peak_indices: Vec<usize>,
    /// m/z values of peaks in this hill.
    pub mz_values: Vec<f64>,
    /// Intensity values of peaks in this hill.
    pub intensities: Vec<f64>,
    /// Retention time values corresponding to each peak.
    pub rt_values: Vec<f64>,
    /// Drift time of each peak's spectrum: `-1` (the unset drift time) when the
    /// spectrum carries none. The source comment calls this *TIMS data, empty
    /// if not available*; in fact it always holds one entry per peak.
    pub drift_times: Vec<f64>,
    /// Per-peak ion mobility from the spectrum's ion-mobility float array, or
    /// `-1` for a peak without one. As for `drift_times`, the source comment's
    /// *empty if not available* does not describe what the code stores.
    pub ion_mobilities: Vec<f64>,
    /// Intensity-weighted mean m/z (hill center). Zero until the hill is
    /// processed, and zero after processing when every intensity is zero.
    pub mz_weighted_mean: f64,
    /// Retention time of the first peak.
    pub rt_start: f64,
    /// Retention time of the last peak.
    pub rt_end: f64,
    /// Retention time at maximum intensity (the first maximum).
    pub rt_apex: f64,
    /// Maximum intensity value in the hill.
    pub intensity_apex: f64,
    /// Sum of all intensities in the hill.
    pub intensity_sum: f64,
    /// Median drift time (`-1` if not available).
    pub drift_time_median: f64,
    /// Intensity-weighted mean ion mobility (median fallback; `-1` if not
    /// available).
    pub ion_mobility_median: f64,
    /// Number of points (peaks) in this hill.
    pub length: usize,
    /// Unique identifier for this hill.
    pub hill_idx: usize,
}

impl Default for Hill {
    fn default() -> Self {
        Self {
            scan_indices: Vec::new(),
            peak_indices: Vec::new(),
            mz_values: Vec::new(),
            intensities: Vec::new(),
            rt_values: Vec::new(),
            drift_times: Vec::new(),
            ion_mobilities: Vec::new(),
            mz_weighted_mean: 0.0,
            rt_start: 0.0,
            rt_end: 0.0,
            rt_apex: 0.0,
            intensity_apex: 0.0,
            intensity_sum: 0.0,
            drift_time_median: -1.0,
            ion_mobility_median: -1.0,
            length: 0,
            hill_idx: 0,
        }
    }
}

/// Candidate isotope peak that can be associated with a monoisotopic hill, the
/// source `Biosaur2Algorithm::IsotopeCandidate`.
///
/// The source documentation: *During isotope pattern detection, candidate hills
/// are evaluated for their suitability as isotope peaks of a monoisotopic hill.
/// This structure stores the association along with quality metrics that
/// quantify how well the candidate matches the expected isotope pattern in
/// terms of mass accuracy and intensity correlation.*
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IsotopeCandidate {
    /// Index ([`Hill::hill_idx`]) of the hill that represents this isotope peak.
    pub hill_idx: usize,
    /// Ordinal isotope number (0 = monoisotopic, 1 = first isotope, ...).
    pub isotope_number: usize,
    /// Mass difference to the expected isotope position, in ppm.
    pub mass_diff_ppm: f64,
    /// Cosine correlation between the monoisotopic and the isotope intensity
    /// traces.
    pub cos_corr: f64,
}

/// Aggregated properties of a detected peptide feature, the source
/// `Biosaur2Algorithm::PeptideFeature`.
///
/// The source documentation: *A peptide feature represents a complete isotope
/// pattern detected across multiple scans. It aggregates information from the
/// monoisotopic hill and its associated isotope hills, providing comprehensive
/// characterization including m/z, retention time, charge state, and intensity
/// metrics. For ion mobility data, drift time and ion mobility values are also
/// included.*
#[derive(Clone, Debug, PartialEq)]
pub struct PeptideFeature {
    /// Monoisotopic m/z value.
    pub mz: f64,
    /// Retention time of first detection.
    pub rt_start: f64,
    /// Retention time of last detection.
    pub rt_end: f64,
    /// Retention time at maximum intensity.
    pub rt_apex: f64,
    /// Maximum intensity of the monoisotopic hill, plus the apex intensities of
    /// the isotopes `iuse` selects.
    pub intensity_apex: f64,
    /// Sum of intensities, over the isotopes `iuse` selects.
    pub intensity_sum: f64,
    /// Detected charge state.
    pub charge: i32,
    /// Number of detected isotope peaks, the monoisotopic one included.
    pub n_isotopes: usize,
    /// Number of scans contributing to the monoisotopic hill.
    pub n_scans: usize,
    /// Calibrated neutral mass of the feature in Da:
    /// `mz * charge - charge * proton_mass` in positive mode,
    /// `mz * charge + charge * proton_mass` in negative mode.
    pub mass_calib: f64,
    /// FAIMS compensation voltage of the feature's group, as the monoisotopic
    /// hill's median drift time; `-1` for non-FAIMS data.
    pub drift_time: f64,
    /// Intensity-weighted mean ion mobility (median fallback; `-1` if not
    /// available).
    pub ion_mobility: f64,
    /// Associated isotope peaks, in isotope order.
    pub isotopes: Vec<IsotopeCandidate>,
    /// [`Hill::hill_idx`] of the monoisotopic hill.
    pub mono_hill_idx: usize,
}

impl Default for PeptideFeature {
    fn default() -> Self {
        Self {
            mz: 0.0,
            rt_start: 0.0,
            rt_end: 0.0,
            rt_apex: 0.0,
            intensity_apex: 0.0,
            intensity_sum: 0.0,
            charge: 0,
            n_isotopes: 0,
            n_scans: 0,
            mass_calib: 0.0,
            drift_time: -1.0,
            ion_mobility: -1.0,
            isotopes: Vec::new(),
            mono_hill_idx: 0,
        }
    }
}

/// How the feature convex hulls are represented, the source parameter
/// `convex_hulls`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ConvexHullMode {
    /// `mass_traces`: one convex hull per contributing hill, built from all its
    /// points (larger featureXML, preserves detailed trace shape).
    MassTraces,
    /// `bounding_box` (the default): a single RT-m/z bounding box per feature
    /// (smaller featureXML, no per-trace detail).
    #[default]
    BoundingBox,
}

impl ConvexHullMode {
    /// The parameter string.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::MassTraces => "mass_traces",
            Self::BoundingBox => "bounding_box",
        }
    }
}

/// Native resource ceilings. The source has none.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Biosaur2Limits {
    /// Most isotope-candidate combinations enumerated for one monoisotopic hill
    /// and charge. The source enumerates the full Cartesian product of the
    /// candidates of every isotope order, which grows as the ninth power of the
    /// candidate count in dense data; the product is computed before the
    /// enumeration starts and a larger one is refused.
    pub max_isotope_combinations: usize,
}

impl Default for Biosaur2Limits {
    fn default() -> Self {
        Self {
            max_isotope_combinations: 1 << 24,
        }
    }
}

/// The log stream a [`Biosaur2Message`] belongs to.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Biosaur2LogLevel {
    /// `OPENMS_LOG_INFO`.
    Info,
    /// `OPENMS_LOG_WARN`.
    Warning,
}

/// A line the source writes to a log stream, returned instead of printed.
///
/// The algorithm is not wired to a log stream; the caller (the
/// `FeatureFinderLFQ` tool) decides where the lines go. Doubles are formatted as
/// the source's `std::ostream` does at its default precision 6.
///
/// Source lines that end in a `StopWatch` timing (`... in 0.01 s (wall), ...`)
/// are recorded without that suffix, and the source's timing-only lines
/// (`PASEF centroiding for group (FAIMS CV=...) took ...`, `Hill calibration
/// for group (FAIMS CV=...) took ...`) are not recorded: their text is not
/// reproducible, and the upstream tests mask it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Biosaur2Message {
    /// The log stream.
    pub level: Biosaur2LogLevel,
    /// The line, without its trailing newline.
    pub text: String,
}

/// Everything [`Biosaur2Algorithm::run_with_details`] produces: the source's
/// three output parameters and the log lines.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Biosaur2Output {
    /// The detected features with their meta information.
    pub feature_map: FeatureMap,
    /// All detected hills that survived filtering and splitting, for
    /// diagnostics and quality control.
    pub hills: Vec<Hill>,
    /// The intermediate peptide feature representations, before conversion to
    /// `FeatureMap` entries.
    pub peptide_features: Vec<PeptideFeature>,
    /// The log lines, in the order the serial source writes them.
    pub messages: Vec<Biosaur2Message>,
}

/// Implementation of the Biosaur2 feature detection workflow for LC-MS1 data,
/// the source `class Biosaur2Algorithm`.
///
/// The public fields are the source's cached parameter members
/// (`updateMembers_`); [`Biosaur2Algorithm::defaults`],
/// [`Biosaur2Algorithm::from_param`] and [`Biosaur2Algorithm::to_param`] are
/// its `DefaultParamHandler` base. Two native fields,
/// [`threads`](Self::threads) and [`limits`](Self::limits), are not
/// parameters.
///
/// ```
/// use openms::analysis::featurefinder::biosaur2::Biosaur2Algorithm;
/// use openms::concept::unique_id::UniqueIdGenerator;
/// use openms::kernel::{MSExperiment, MSSpectrum, Peak1D};
///
/// // the upstream class test's `run(FeatureMap&)` section
/// let mut experiment = MSExperiment::new();
/// for i in 0..10 {
///     let mut spectrum = MSSpectrum { rt: f64::from(i), ms_level: 1, ..MSSpectrum::default() };
///     spectrum.peaks = vec![
///         Peak1D::new(1000.0, 10000.0),
///         Peak1D::new(1001.003355, 6000.0),
///         Peak1D::new(1000.993355, 100.0),
///     ];
///     experiment.spectra.push(spectrum);
/// }
/// let mut algorithm = Biosaur2Algorithm::new();
/// algorithm.itol = 20.0;
/// algorithm.minmz = 900.0;
/// algorithm.maxmz = 1100.0;
/// algorithm.set_ms_data(experiment);
/// let features = algorithm.run(&mut UniqueIdGenerator::from_seed(1))?;
/// assert_eq!(features.len(), 1);
/// assert_eq!(features.features[0].charge, 1);
/// # Ok::<(), openms::Error>(())
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct Biosaur2Algorithm {
    /// `mini`, default 1: minimum intensity threshold.
    pub mini: f64,
    /// `minmz`, default 350: minimum m/z value.
    pub minmz: f64,
    /// `maxmz`, default 1500: maximum m/z value.
    pub maxmz: f64,
    /// `htol`, default 8: mass accuracy in ppm for combining peaks into hills.
    pub htol: f64,
    /// `itol`, default 8: mass accuracy in ppm for isotopic patterns.
    pub itol: f64,
    /// `hvf`, default 1.3, at least 1: hill valley factor for splitting hills.
    pub hvf: f64,
    /// `ivf`, default 5, at least 1: isotope valley factor for splitting
    /// isotope patterns.
    pub ivf: f64,
    /// `minlh`, default 2, at least 1: minimum number of scans for a hill.
    pub minlh: usize,
    /// `pasefmini`, default 100: minimum combined intensity for PASEF/TIMS
    /// clusters after m/z-ion-mobility centroiding.
    pub pasefmini: f64,
    /// `pasefminlh`, default 1, at least 1: minimum number of raw points per
    /// PASEF/TIMS cluster during centroiding.
    pub pasefminlh: usize,
    /// `cmin`, default 1, at least 1: minimum charge state.
    pub cmin: i32,
    /// `cmax`, default 6, at least 1: maximum charge state.
    pub cmax: i32,
    /// `iuse`, default 0, at least -1: number of isotopes for intensity
    /// calculation (0 = mono only, -1 = all, 1 = mono + first, etc.).
    pub iuse: i32,
    /// `nm`, default false: negative mode (affects the neutral mass
    /// calculation).
    pub negative_mode: bool,
    /// `tof`, default false: enable TOF-specific intensity filtering.
    pub tof_mode: bool,
    /// `profile`, default false: enable profile mode processing (centroid
    /// spectra using `PeakPickerHiRes`).
    pub profile_mode: bool,
    /// `paseftol`, default 0.05: ion mobility accuracy for linking peaks into
    /// hills and grouping isotopes (0 disables IM-based gating). The default is
    /// for 1/K0 units (e.g. 0.03-0.08); for CCS data use larger values (e.g.
    /// 5-20 square angstroms).
    pub paseftol: f64,
    /// `use_hill_calib`, default false: enable automatic hill mass tolerance
    /// calibration.
    pub use_hill_calib: bool,
    /// `ignore_iso_calib`, default false: disable automatic isotope mass error
    /// calibration.
    pub ignore_iso_calib: bool,
    /// `hrttol`, default 10: maximum allowed RT difference (in seconds) between
    /// the monoisotopic hill apex and an isotope hill apex when assembling
    /// isotope patterns (0 disables RT gating).
    pub hrttol: f64,
    /// `convex_hulls`, default `bounding_box`: representation of the feature
    /// convex hulls in the output `FeatureMap`.
    pub convex_hulls: ConvexHullMode,
    /// `faims_merge_features`, default true: for FAIMS data with several
    /// compensation voltages, merge features representing the same analyte
    /// detected at different CV values into a single feature. Only features
    /// with *different* FAIMS CV values are merged (same CV = different
    /// analytes). Has no effect on non-FAIMS data.
    pub faims_merge_features: bool,
    /// Native: worker threads for the FAIMS groups, which the source processes
    /// in an OpenMP `parallel for`. The result does not depend on it.
    pub threads: Threads,
    /// Native: resource ceilings.
    pub limits: Biosaur2Limits,
    ms_data: MSExperiment,
}

impl Default for Biosaur2Algorithm {
    fn default() -> Self {
        Self {
            mini: 1.0,
            minmz: 350.0,
            maxmz: 1500.0,
            htol: 8.0,
            itol: 8.0,
            hvf: 1.3,
            ivf: 5.0,
            minlh: 2,
            pasefmini: 100.0,
            pasefminlh: 1,
            cmin: 1,
            cmax: 6,
            iuse: 0,
            negative_mode: false,
            tof_mode: false,
            profile_mode: false,
            paseftol: 0.05,
            use_hill_calib: false,
            ignore_iso_calib: false,
            hrttol: 10.0,
            convex_hulls: ConvexHullMode::BoundingBox,
            faims_merge_features: true,
            threads: Threads::default(),
            limits: Biosaur2Limits::default(),
            ms_data: MSExperiment::new(),
        }
    }
}

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

/// `Exception::InvalidValue(..., message, value).what()`.
pub(crate) fn source_invalid_value(message: &str, value: &str) -> Error {
    Error::InvalidValue(format!(
        "the value '{value}' was used but is not valid; {message}"
    ))
}

/// A double as `std::ostream <<` prints it at the default precision 6.
pub(crate) fn g6(value: f64) -> String {
    ostream_g(value, 6)
}

pub(crate) fn info(messages: &mut Vec<Biosaur2Message>, text: String) {
    messages.push(Biosaur2Message {
        level: Biosaur2LogLevel::Info,
        text,
    });
}

pub(crate) fn warn(messages: &mut Vec<Biosaur2Message>, text: String) {
    messages.push(Biosaur2Message {
        level: Biosaur2LogLevel::Warning,
        text,
    });
}

fn float_entry(name: &str, value: f64, description: &str, min: f64) -> ParamEntry {
    ParamEntry {
        name: name.into(),
        description: description.into(),
        value: ParamValue::Float(value),
        min_float: min,
        ..ParamEntry::default()
    }
}

fn int_entry(name: &str, value: i64, description: &str, min: i32) -> ParamEntry {
    ParamEntry {
        name: name.into(),
        description: description.into(),
        value: ParamValue::Integer(value),
        min_int: min,
        ..ParamEntry::default()
    }
}

fn string_entry(name: &str, value: &str, description: &str, valid: &[&str]) -> ParamEntry {
    ParamEntry {
        name: name.into(),
        description: description.into(),
        value: ParamValue::String(value.into()),
        valid_strings: valid.iter().map(|v| (*v).to_owned()).collect(),
        ..ParamEntry::default()
    }
}

fn flag_entry(name: &str, value: bool, description: &str) -> ParamEntry {
    string_entry(
        name,
        if value { "true" } else { "false" },
        description,
        &["true", "false"],
    )
}

fn to_param_int(value: usize, name: &str) -> Result<i64> {
    i32::try_from(value)
        .map(i64::from)
        .map_err(|_| bad(format!("'{name}' does not fit the source 32-bit parameter")))
}

/// Which MS1 spectra of the stored experiment form one processing group, and
/// the group's FAIMS compensation voltage (`None` for non-FAIMS data, the
/// source's NaN).
struct Group {
    cv: Option<CompensationVoltage>,
    /// Indices into the MS1 working set, in input order.
    members: Vec<usize>,
}

/// Split the MS1 working set by FAIMS compensation voltage with
/// [`ImDataConverter::split_by_faims_cv`], on the skeleton `run_with_details`
/// builds (each spectrum's native id is its working position). The split's log
/// records join `messages` in the source's order; the spectra it skips (MS1
/// spectra without a detected voltage) belong to no group.
fn faims_groups(
    mut skeleton: MSExperiment,
    messages: &mut Vec<Biosaur2Message>,
) -> Result<Vec<Group>> {
    let split = ImDataConverter::split_by_faims_cv(&mut skeleton)?;
    for message in &split.messages {
        let text = message.text();
        match message.level() {
            FaimsSplitLogLevel::Info => info(messages, text),
            FaimsSplitLogLevel::Warning => warn(messages, text),
        }
    }
    let mut groups = Vec::with_capacity(split.groups.len());
    for group in split.groups {
        let members = group
            .experiment
            .spectra
            .iter()
            .map(|spectrum| {
                spectrum
                    .native_id
                    .parse::<usize>()
                    .map_err(|_| bad("FAIMS split lost a working position"))
            })
            .collect::<Result<Vec<_>>>()?;
        groups.push(Group {
            cv: group.key.voltage(),
            members,
        });
    }
    Ok(groups)
}

/// What processing one group produced: the source `processFAIMSGroup_`'s
/// two outputs, the spectra its PASEF centroiding replaced, and its log lines.
struct GroupResult {
    hills: Vec<Hill>,
    features: Vec<PeptideFeature>,
    /// `(member position, centroided spectrum)`.
    replaced: Vec<(usize, MSSpectrum)>,
    messages: Vec<Biosaur2Message>,
}

impl Biosaur2Algorithm {
    /// The algorithm with every parameter at its default, the source default
    /// constructor. The stored MS data is empty.
    pub fn new() -> Self {
        Self::default()
    }

    /// The source parameter defaults as the constructor declares them, in
    /// declaration order, with descriptions and restrictions.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure, which the constant tree cannot
    /// trigger.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters of this algorithm in the layout of
    /// [`Biosaur2Algorithm::defaults`], carrying the current values (source
    /// `getParameters`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `minlh` or `pasefminlh` does not
    /// fit the source's 32-bit integer parameter.
    pub fn to_param(&self) -> Result<Param> {
        let entries = vec![
            float_entry("mini", self.mini, "Minimum intensity threshold", 0.0),
            float_entry("minmz", self.minmz, "Minimum m/z value", 0.0),
            float_entry("maxmz", self.maxmz, "Maximum m/z value", 0.0),
            float_entry(
                "htol",
                self.htol,
                "Mass accuracy in ppm for combining peaks into hills",
                0.0,
            ),
            float_entry(
                "itol",
                self.itol,
                "Mass accuracy in ppm for isotopic patterns",
                0.0,
            ),
            float_entry(
                "hvf",
                self.hvf,
                "Hill valley factor for splitting hills",
                1.0,
            ),
            float_entry(
                "ivf",
                self.ivf,
                "Isotope valley factor for splitting isotope patterns",
                1.0,
            ),
            int_entry(
                "minlh",
                to_param_int(self.minlh, "minlh")?,
                "Minimum number of scans for a hill",
                1,
            ),
            float_entry(
                "pasefmini",
                self.pasefmini,
                "Minimum combined intensity for PASEF/TIMS clusters after m/z\u{2013}ion-mobility centroiding.",
                0.0,
            ),
            int_entry(
                "pasefminlh",
                to_param_int(self.pasefminlh, "pasefminlh")?,
                "Minimum number of raw points per PASEF/TIMS cluster during centroiding.",
                1,
            ),
            int_entry("cmin", i64::from(self.cmin), "Minimum charge state", 1),
            int_entry("cmax", i64::from(self.cmax), "Maximum charge state", 1),
            int_entry(
                "iuse",
                i64::from(self.iuse),
                "Number of isotopes for intensity calculation (0=mono only, -1=all, 1=mono+first, etc.)",
                -1,
            ),
            flag_entry(
                "nm",
                self.negative_mode,
                "Negative mode (affects neutral mass calculation)",
            ),
            flag_entry(
                "tof",
                self.tof_mode,
                "Enable TOF-specific intensity filtering",
            ),
            flag_entry(
                "profile",
                self.profile_mode,
                "Enable profile mode processing (centroid spectra using PeakPickerHiRes)",
            ),
            float_entry(
                "paseftol",
                self.paseftol,
                "Ion mobility accuracy for linking peaks into hills and grouping isotopes (0 = disable IM-based gating). Default is for 1/K0 units (e.g., 0.03-0.08). For CCS data, use larger values (e.g., 5-20 square angstroms).",
                0.0,
            ),
            flag_entry(
                "use_hill_calib",
                self.use_hill_calib,
                "Enable automatic hill mass tolerance calibration",
            ),
            flag_entry(
                "ignore_iso_calib",
                self.ignore_iso_calib,
                "Disable automatic isotope mass error calibration",
            ),
            float_entry(
                "hrttol",
                self.hrttol,
                "Maximum allowed RT difference (in seconds) between monoisotopic hill apex and isotope hill apex when assembling isotope patterns (0 disables RT gating).",
                0.0,
            ),
            string_entry(
                "convex_hulls",
                self.convex_hulls.as_str(),
                "Representation of feature convex hulls in the output FeatureMap. 'bounding_box' stores a single RT\u{2013}m/z bounding box per feature (smaller featureXML, no per-trace detail), whereas 'mass_traces' stores one convex hull per contributing hill using all mass-trace points (larger featureXML, preserves detailed trace shape).",
                &["mass_traces", "bounding_box"],
            ),
            flag_entry(
                "faims_merge_features",
                self.faims_merge_features,
                "For FAIMS data with multiple compensation voltages: Merge features representing the same analyte detected at different CV values into a single feature. Only features with DIFFERENT FAIMS CV values are merged (same CV = different analytes). Has no effect on non-FAIMS data.",
            ),
        ];
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries,
            nodes: Vec::new(),
        })
    }

    /// Build the algorithm from parameters, as source `setParameters` followed
    /// by `updateMembers_`. Missing parameters take their defaults; the stored
    /// MS data is empty and the native fields keep their defaults. Unknown
    /// parameters are ignored, as the source only warns; see
    /// [`Biosaur2Algorithm::from_param_with_warnings`].
    ///
    /// The source's `updateMembers_` also writes every parameter to
    /// `OPENMS_LOG_DEBUG`; that debug line is not reproduced.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong value type
    /// or violates its restriction (the source `Exception::InvalidParameter`).
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`Biosaur2Algorithm::from_param`], also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`Biosaur2Algorithm::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(BIOSAUR2_ALGORITHM_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, Self::from_complete_param)
    }

    /// Apply parameters to this algorithm, keeping its stored MS data and
    /// native fields; the source `setParameters` on an existing object.
    ///
    /// # Errors
    ///
    /// As [`Biosaur2Algorithm::from_param`]; the algorithm is unchanged then.
    pub fn set_parameters(&mut self, parameters: &Param) -> Result<Vec<String>> {
        let (updated, warnings) = Self::from_param_with_warnings(parameters)?;
        let ms_data = std::mem::take(&mut self.ms_data);
        let (threads, limits) = (self.threads, self.limits);
        *self = Self {
            threads,
            limits,
            ms_data,
            ..updated
        };
        Ok(warnings)
    }

    fn from_complete_param(param: &Param) -> Result<Self> {
        let float = |key: &str| param.value(key)?.to_f64();
        let flag = |key: &str| param.value(key)?.to_bool();
        let int = |key: &str| param.value(key)?.to_i64();
        let size = |key: &str| {
            usize::try_from(int(key)?).map_err(|_| bad(format!("'{key}' must not be negative")))
        };
        let int32 = |key: &str| {
            i32::try_from(int(key)?)
                .map_err(|_| bad(format!("'{key}' does not fit a 32-bit integer")))
        };
        let convex_hulls = match param.value("convex_hulls")?.as_str()? {
            "mass_traces" => ConvexHullMode::MassTraces,
            _ => ConvexHullMode::BoundingBox,
        };
        Ok(Self {
            mini: float("mini")?,
            minmz: float("minmz")?,
            maxmz: float("maxmz")?,
            htol: float("htol")?,
            itol: float("itol")?,
            hvf: float("hvf")?,
            ivf: float("ivf")?,
            minlh: size("minlh")?,
            pasefmini: float("pasefmini")?,
            pasefminlh: size("pasefminlh")?,
            cmin: int32("cmin")?,
            cmax: int32("cmax")?,
            iuse: int32("iuse")?,
            negative_mode: flag("nm")?,
            tof_mode: flag("tof")?,
            profile_mode: flag("profile")?,
            paseftol: float("paseftol")?,
            use_hill_calib: flag("use_hill_calib")?,
            ignore_iso_calib: flag("ignore_iso_calib")?,
            hrttol: float("hrttol")?,
            convex_hulls,
            faims_merge_features: flag("faims_merge_features")?,
            ..Self::default()
        })
    }

    /// Set the MS data used for feature detection, the source `setMSData`
    /// (both the copy and the move overload: pass a clone to keep yours).
    pub fn set_ms_data(&mut self, ms_data: MSExperiment) {
        self.ms_data = ms_data;
    }

    /// The stored MS data, the source `getMSData() const`.
    pub fn ms_data(&self) -> &MSExperiment {
        &self.ms_data
    }

    /// The stored MS data for modification, the source non-const `getMSData()`.
    pub fn ms_data_mut(&mut self) -> &mut MSExperiment {
        &mut self.ms_data
    }

    /// Execute the workflow and return only the feature map, the source
    /// `run(FeatureMap&)`, which *discards the intermediate hills and peptide
    /// feature vectors*.
    ///
    /// # Errors
    ///
    /// As [`Biosaur2Algorithm::run_with_details`].
    pub fn run(&mut self, generator: &mut UniqueIdGenerator) -> Result<FeatureMap> {
        Ok(self.run_with_details(generator)?.feature_map)
    }

    /// Execute the Biosaur2 workflow on the stored MS1 experiment, the source
    /// `run(FeatureMap&, std::vector<Hill>&, std::vector<PeptideFeature>&)`.
    ///
    /// The source documentation lists the pipeline:
    ///
    /// 1. filter the input to MS1 spectra;
    /// 2. optionally centroid profile data;
    /// 3. apply TOF-specific intensity filtering if enabled;
    /// 4. group spectra by FAIMS compensation voltage if applicable;
    /// 5. detect hills (continuous m/z traces across scans);
    /// 6. split hills at valley points;
    /// 7. detect isotope patterns and assemble peptide features;
    /// 8. convert the features to a `FeatureMap`.
    ///
    /// *This method modifies the stored MS data.* On success every spectrum
    /// with an MS level other than 1 is gone; with `profile` the rest are
    /// centroided, with `tof` TOF-filtered, and spectra with per-peak ion
    /// mobility are PASEF-centroided when ion-mobility gating is active. Data
    /// with FAIMS compensation voltages is split by voltage and reassembled in
    /// retention-time order (the source's `std::sort`, whose order of equal
    /// retention times this follows), keeping the chromatograms and settings.
    /// The source notes that the run is therefore not idempotent: call
    /// [`set_ms_data`](Self::set_ms_data) again before re-running.
    ///
    /// Unique ids for the features and the map come from `generator` (the
    /// source uses its process-wide singleton); the map carries one empty
    /// protein identification, as the source's `resize(1)` leaves it.
    ///
    /// # Errors
    ///
    /// The stored MS data is left **unchanged** by every error, where the
    /// source has already removed the non-MS1 spectra and centroided or
    /// filtered the rest when it throws. Returns [`Error::InvalidValue`]:
    ///
    /// - with the source's text `the value '' was used but is not valid; No MS1
    ///   spectra found in input experiment!` when no MS1 spectrum is stored;
    /// - `Ion mobility array shorter than peak list.` (with the peak index as
    ///   the value) when a spectrum's ion-mobility array is shorter than its
    ///   peak list, as the source throws;
    /// - natively, when a parameter field holds a value its parameter would
    ///   refuse or is not finite, when an MS1 spectrum has a non-finite
    ///   retention time, or a peak inside the configured m/z and intensity
    ///   window has a non-finite m/z or intensity (the source sorts and bins
    ///   such peaks with undefined behaviour), and when an isotope enumeration
    ///   exceeds [`Biosaur2Limits::max_isotope_combinations`];
    ///
    /// and whatever the profile centroiding (`PeakPickerHiRes`), the TOF
    /// filter's peak selection, the FAIMS split and the FAIMS feature merge
    /// refuse.
    pub fn run_with_details(
        &mut self,
        generator: &mut UniqueIdGenerator,
    ) -> Result<Biosaur2Output> {
        self.validate()?;
        let mut messages = Vec::new();

        // Steps 1-3 run on a working set that borrows the stored spectra and
        // owns only the ones a step replaces, so that nothing is committed
        // before the whole run has succeeded.
        let ms1: Vec<usize> = self
            .ms_data
            .spectra
            .iter()
            .enumerate()
            .filter(|(_, spectrum)| spectrum.ms_level == 1)
            .map(|(index, _)| index)
            .collect();
        let mut replaced: Vec<Option<MSSpectrum>> = vec![None; ms1.len()];
        if self.profile_mode {
            preprocess::centroid_profile_spectra(
                &self.ms_data.spectra,
                &ms1,
                &mut replaced,
                &mut messages,
            )?;
        }
        if self.tof_mode {
            preprocess::process_tof(&self.ms_data.spectra, &ms1, &mut replaced, &mut messages)?;
        }
        info(&mut messages, format!("Loaded {} MS1 spectra", ms1.len()));
        if ms1.is_empty() {
            return Err(source_invalid_value(
                "No MS1 spectra found in input experiment!",
                "",
            ));
        }

        let working: Vec<&MSSpectrum> = ms1
            .iter()
            .zip(&replaced)
            .map(|(&index, replacement)| {
                replacement.as_ref().unwrap_or(&self.ms_data.spectra[index])
            })
            .collect();
        for spectrum in &working {
            if !spectrum.rt.is_finite() {
                return Err(bad(
                    "Biosaur2 needs a finite retention time on every MS1 spectrum",
                ));
            }
        }

        // The CCS warning, for the first spectrum with an ion-mobility array.
        for spectrum in &working {
            if spectrum.contains_im_data() {
                let (_, unit) = spectrum.im_data()?;
                if unit == DriftTimeUnit::CollisionCrossSection && self.paseftol < 1.0 {
                    warn(
                        &mut messages,
                        format!(
                            "Warning: Ion mobility data is in CCS units (square angstroms), but paseftol = {} appears to be set for 1/K0 data. For CCS data, consider using larger values (e.g., 5-20).",
                            g6(self.paseftol)
                        ),
                    );
                }
                break;
            }
        }

        // Step 4: FAIMS grouping. `FAIMSHelper::getCompensationVoltages` and
        // `IMDataConverter::splitByFAIMSCV` only read each spectrum's MS level,
        // drift time and drift time unit (and the split validates the
        // settings), so both run on a peak-free skeleton of the working set
        // whose native ids carry the working positions.
        let skeleton = MSExperiment {
            spectra: working
                .iter()
                .enumerate()
                .map(|(k, spectrum)| MSSpectrum {
                    ms_level: spectrum.ms_level,
                    drift_time: spectrum.drift_time,
                    drift_time_unit: spectrum.drift_time_unit,
                    native_id: k.to_string(),
                    ..MSSpectrum::default()
                })
                .collect(),
            settings: self.ms_data.settings.clone(),
            ..MSExperiment::default()
        };
        let voltages = FaimsHelper::get_compensation_voltages(&skeleton)?;
        for warning in &voltages.warnings {
            warn(&mut messages, warning.clone());
        }
        let faims = !voltages.voltages.is_empty();
        let groups = if faims {
            faims_groups(skeleton, &mut messages)?
        } else {
            info(
                &mut messages,
                "Not FAIMS compensation voltages found in the data. Processing in-place."
                    .to_owned(),
            );
            vec![Group {
                cv: None,
                members: (0..working.len()).collect(),
            }]
        };

        // Steps 5-7, one group at a time (the source's `parallel for`).
        let results: Vec<Result<GroupResult>> = map_collect(&groups, self.threads, |group| {
            let spectra: Vec<&MSSpectrum> = group.members.iter().map(|&k| working[k]).collect();
            self.process_group(group.cv, &spectra)
        });
        let mut group_results = Vec::with_capacity(results.len());
        for result in results {
            group_results.push(result?);
        }

        // Step 8 and the output assembly.
        let mut output = Biosaur2Output::default();
        let mut feature_map = FeatureMap::new();
        for result in &mut group_results {
            messages.append(&mut result.messages);
            if faims {
                let mut local = FeatureMap::new();
                if !result.hills.is_empty() && !result.features.is_empty() {
                    local = output::convert_to_feature_map(
                        self,
                        &result.features,
                        &result.hills,
                        generator,
                        &mut messages,
                    )?;
                }
                feature_map.features.append(&mut local.features);
            } else if !result.hills.is_empty() && !result.features.is_empty() {
                feature_map = output::convert_to_feature_map(
                    self,
                    &result.features,
                    &result.hills,
                    generator,
                    &mut messages,
                )?;
            }
        }
        if faims && self.faims_merge_features {
            let before = feature_map.len();
            FeatureOverlapFilter::merge_faims_features(&mut feature_map, 5.0, 0.05)?;
            if feature_map.len() < before {
                info(
                    &mut messages,
                    format!(
                        "FAIMS feature merge: {before} -> {} features (merged {})",
                        feature_map.len(),
                        before - feature_map.len()
                    ),
                );
            }
        }
        output::finish_feature_map(&mut feature_map, generator);

        // Commit: the processed MS1 spectra replace the stored experiment.
        for (group, result) in groups.iter().zip(&mut group_results) {
            for (position, spectrum) in result.replaced.drain(..) {
                replaced[group.members[position]] = Some(spectrum);
            }
        }
        let mut taken = std::mem::take(&mut self.ms_data.spectra);
        let mut processed: Vec<Option<MSSpectrum>> = ms1
            .iter()
            .zip(replaced)
            .map(|(&index, replacement)| {
                Some(replacement.unwrap_or_else(|| std::mem::take(&mut taken[index])))
            })
            .collect();
        drop(taken);
        let mut spectra = Vec::with_capacity(ms1.len());
        if faims {
            // The groups in ascending voltage order, each in input order; the
            // spectra no group took are gone, as `splitByFAIMSCV` destroys them.
            for group in &groups {
                for &k in &group.members {
                    if let Some(spectrum) = processed[k].take() {
                        spectra.push(spectrum);
                    }
                }
            }
            source_sort_by(&mut spectra, |a, b| a.rt < b.rt)?;
        } else {
            spectra.extend(processed.into_iter().flatten());
        }
        self.ms_data.spectra = spectra;

        for result in group_results {
            output.hills.extend(result.hills);
            output.peptide_features.extend(result.features);
        }
        output.feature_map = feature_map;
        output.messages = messages;
        Ok(output)
    }

    /// Refuse parameter fields that their parameter would refuse.
    fn validate(&self) -> Result<()> {
        let floats = [
            ("mini", self.mini, 0.0),
            ("minmz", self.minmz, 0.0),
            ("maxmz", self.maxmz, 0.0),
            ("htol", self.htol, 0.0),
            ("itol", self.itol, 0.0),
            ("hvf", self.hvf, 1.0),
            ("ivf", self.ivf, 1.0),
            ("pasefmini", self.pasefmini, 0.0),
            ("paseftol", self.paseftol, 0.0),
            ("hrttol", self.hrttol, 0.0),
        ];
        for (name, value, min) in floats {
            if !value.is_finite() || value < min {
                return Err(bad(format!(
                    "Biosaur2 parameter '{name}' must be finite and at least {}",
                    g6(min)
                )));
            }
        }
        if self.minlh < 1 || self.pasefminlh < 1 || self.cmin < 1 || self.cmax < 1 || self.iuse < -1
        {
            return Err(bad(
                "Biosaur2 parameters 'minlh', 'pasefminlh', 'cmin' and 'cmax' must be at least 1 and 'iuse' at least -1",
            ));
        }
        to_param_int(self.minlh, "minlh")?;
        to_param_int(self.pasefminlh, "pasefminlh")?;
        Ok(())
    }

    /// The source `processFAIMSGroup_`: the complete per-group pipeline
    /// (optional PASEF centroiding, optional hill calibration, hill detection,
    /// processing and splitting, isotope pattern detection) on the MS1 spectra
    /// of one FAIMS compensation voltage, or of the whole run for non-FAIMS
    /// data (`cv` is `None`, the source's NaN).
    fn process_group(
        &self,
        cv: Option<CompensationVoltage>,
        spectra: &[&MSSpectrum],
    ) -> Result<GroupResult> {
        let mut messages = Vec::new();
        let mut result = GroupResult {
            hills: Vec::new(),
            features: Vec::new(),
            replaced: Vec::new(),
            messages: Vec::new(),
        };
        if spectra.is_empty() {
            return Ok(result);
        }
        let cv_text = g6(cv.map_or(f64::NAN, CompensationVoltage::volts));

        let mut any_im_array = false;
        let mut any_missing_im = false;
        let mut with_data = 0usize;
        for spectrum in spectra {
            if ImTypes::determine_im_format(spectrum) != IonMobilityFormat::None {
                any_im_array = true;
                with_data += 1;
            } else {
                any_missing_im = true;
            }
        }
        if any_im_array {
            info(
                &mut messages,
                format!(
                    "Group (FAIMS CV={cv_text}) has {with_data} / {} spectra with ion-mobility data.",
                    spectra.len()
                ),
            );
        }

        // The m/z step of the PASEF centroiding, from the group's maximum m/z.
        let mut max_mz_value = 0.0;
        for spectrum in spectra {
            for peak in &spectrum.peaks {
                let intensity = f64::from(peak.intensity);
                if intensity < self.mini || peak.mz < self.minmz || peak.mz > self.maxmz {
                    continue;
                }
                if peak.mz > max_mz_value {
                    max_mz_value = peak.mz;
                }
            }
        }
        let mz_step = if max_mz_value > 0.0 && self.htol > 0.0 {
            self.htol * 1e-6 * max_mz_value
        } else {
            0.0
        };

        // The source skips its own centroiding only when every spectrum with
        // per-peak ion mobility is annotated `IM_CENTROIDED`. The native
        // spectrum carries no ion-mobility peak type, so that annotation cannot
        // reach this point and the data is treated as not yet centroided (the
        // source's reading of `UNKNOWN` and `IM_PROFILE`).
        let use_im_group = any_im_array && !any_missing_im && self.paseftol > 0.0;
        let mut owned: Vec<Option<MSSpectrum>> = vec![None; spectra.len()];
        if use_im_group && mz_step > 0.0 {
            info(
                &mut messages,
                format!(
                    "Applying PASEF/TIMS centroiding for group (FAIMS CV={cv_text}) with paseftol={}",
                    g6(self.paseftol)
                ),
            );
            preprocess::centroid_pasef_data(
                self,
                spectra,
                &mut owned,
                mz_step,
                self.paseftol,
                &mut messages,
            )?;
        } else if self.paseftol > 0.0 && any_im_array && any_missing_im {
            warn(
                &mut messages,
                format!(
                    "Disabling ion-mobility gating for group (FAIMS CV={cv_text}) due to missing/partial IM arrays; proceeding in 1D m/z space."
                ),
            );
        }
        let group: Vec<&MSSpectrum> = spectra
            .iter()
            .zip(&owned)
            .map(|(&spectrum, replacement)| replacement.as_ref().unwrap_or(spectrum))
            .collect();

        let mut calibrated_htol = self.htol;
        if self.use_hill_calib {
            info(
                &mut messages,
                format!(
                    "Performing hill mass tolerance calibration for group (FAIMS CV={cv_text})..."
                ),
            );
            let sample_size = group.len().min(1000);
            let start = if group.len() > 1000 {
                group.len() / 2 - 500
            } else {
                0
            };
            let end = (start + sample_size).min(group.len());
            let mut mass_diffs = Vec::new();
            hills::detect_hills(
                self,
                &group[start..end],
                self.htol,
                false,
                Some(&mut mass_diffs),
            )?;
            if !mass_diffs.is_empty() {
                let (_, sigma) = hills::calibrate_mass(&mass_diffs, 0.05)?;
                let five_sigma = 5.0 * sigma;
                calibrated_htol = if five_sigma < self.htol {
                    five_sigma
                } else {
                    self.htol
                };
                info(
                    &mut messages,
                    format!(
                        "Automatically optimized htol parameter for group (FAIMS CV={cv_text}): {} ppm (was {} ppm)",
                        g6(calibrated_htol),
                        g6(self.htol)
                    ),
                );
            }
        }

        let detected = hills::detect_hills(self, &group, calibrated_htol, use_im_group, None)?;
        info(
            &mut messages,
            format!(
                "Hill detection for group (FAIMS CV={cv_text}) found {} hills",
                detected.len()
            ),
        );
        let processed = hills::process_hills(&detected, self.minlh)?;
        drop(detected);
        info(
            &mut messages,
            format!(
                "Hill preprocessing for group (FAIMS CV={cv_text}) kept {} hills",
                processed.len()
            ),
        );
        let mut group_hills = hills::split_hills(&processed, self.hvf, self.minlh)?;
        drop(processed);
        info(
            &mut messages,
            format!(
                "Hill splitting for group (FAIMS CV={cv_text}) produced {} hills",
                group_hills.len()
            ),
        );

        let mut features =
            isotopes::detect_isotope_patterns(self, &mut group_hills, use_im_group, &mut messages)?;
        if cv.is_none() {
            // Non-FAIMS data: no FAIMS CV annotation downstream.
            for feature in &mut features {
                feature.drift_time = ImTypes::DRIFTTIME_NOT_SET;
            }
        }
        info(
            &mut messages,
            format!(
                "Isotope pattern detection for group (FAIMS CV={cv_text}) produced {} features",
                features.len()
            ),
        );

        result.hills = group_hills;
        result.features = features;
        result.replaced = owned
            .into_iter()
            .enumerate()
            .filter_map(|(position, spectrum)| spectrum.map(|s| (position, s)))
            .collect();
        result.messages = messages;
        Ok(result)
    }
}
