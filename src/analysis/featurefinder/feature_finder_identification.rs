// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! ID-guided MS1 feature detection.
//!
//! Ports `FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h` with its `.cpp`,
//! the algorithm behind the TOPP tool `FeatureFinderIdentification`. Peptide
//! identifications are the targets: an assay library is built on the fly (one
//! assay per modified sequence, charge state and elution region), the MS1
//! traces of the isotopes are extracted, `MRMFeatureFinderScoring` picks and
//! scores the peak groups, and the candidates are classified, filtered and
//! optionally fitted with an elution model. See
//! `docs/FEATURE_FINDER_IDENTIFICATION_SUPPORT.md` for the API mapping, the
//! preserved source conventions and the evidence.
//!
//! # Example
//!
//! ```
//! use openms::analysis::featurefinder::feature_finder_identification::FeatureFinderIdentificationAlgorithm;
//!
//! let algorithm = FeatureFinderIdentificationAlgorithm::new()?;
//! assert_eq!(algorithm.parameters().value("extract:batch_size")?.to_i64()?, 5000);
//! # Ok::<(), openms::Error>(())
//! ```

use crate::analysis::featurefinder::elution_model_fitter::ElutionModelFitReport;
use crate::analysis::featurefinder::ms1_targeted::{
    self as common, apply_common_finder_parameters, ensure_convex_hulls, extract_ms1_chromatograms,
    feature_bounds, keep_ms1,
};
use crate::analysis::openswath::data_access_helper::convert_targeted_exp;
use crate::analysis::openswath::mrm_feature_finder_scoring::{
    MRMFeatureFinderScoring, TransitionGroupMapType,
};
use crate::analysis::openswath::reaction_monitoring_transition::ReactionMonitoringTransition;
use crate::analysis::openswath::spectrum_access_openms::SimpleOpenMSSpectraFactory;
use crate::analysis::openswath::swath_map::SwathMap;
use crate::analysis::openswath::targeted_experiment::TargetedExperiment;
use crate::analysis::openswath::targeted_experiment_helper::{
    Peptide, Protein, RTType, RetentionTime,
};
use crate::analysis::transformations::TransformationDescription;
use crate::chemistry::isotopes::ProbabilityPrecision;
use crate::chemistry::{
    AASequence, CoarseIsotopePatternGenerator, CoarseMassMode, EmpiricalFormula,
    IsotopeDistribution, ModificationsDB, PeptideFragmentType, SequenceModification,
    TermSpecificity,
};
use crate::concept::constants::user_param::{FAIMS_CV, IM};
use crate::concept::progress_logger::ProgressLogger;
use crate::concept::unique_id::{HasUniqueId, UniqueIdGenerator};
use crate::constants::{C13C12_MASSDIFF_U, PROTON_MASS_U};
use crate::format::file_info::text_format::to_str;
use crate::identification::{PeptideHit, PeptideIdentification, ProteinIdentification};
use crate::kernel::faims_helper::FaimsHelper;
use crate::kernel::im_data_converter::ImDataConverter;
use crate::kernel::{Feature, FeatureMap, MSExperiment};
use crate::math::statistic_functions::median_of_sorted;
use crate::metadata::{DriftTimeUnit, ImTypes, IonMobilityFormat, MetaInfo, MetaValue};
use crate::param::{DefaultParamHandler, Param, ParamValue};
use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::{Arc, Mutex};

/// The handler name (source `DefaultParamHandler("FeatureFinderIdentificationAlgorithm")`).
const NAME: &str = "FeatureFinderIdentificationAlgorithm";

/// Largest number of identifications plus seeds (offset peptides included) one
/// [`run`](FeatureFinderIdentificationAlgorithm::run) accepts.
///
/// Native ceiling; the source has none. It is checked before anything is
/// allocated, so a refusal leaves the algorithm and the output unchanged.
pub const MAX_TARGETS: usize = 1 << 26;

/// Largest number of lines [`FeatureFinderIdentificationReport::messages`]
/// keeps; later lines are counted in
/// [`FeatureFinderIdentificationReport::dropped_messages`].
pub const MAX_REPORT_MESSAGES: usize = 10_000;

/// RT extraction window of a seed in seconds. Source `seed_rt_window_`; smaller
/// than the window of an identification, as the apex of a seed is known.
const SEED_RT_WINDOW: f64 = 60.0;

/// The mass tag of the pseudo sequence of seeds and offset peptides. The source
/// calls this a "superhack": a modification that never occurs, so the sequence
/// is unique and is not mistaken for a real one.
const PSEUDO_SEQUENCE: &str = "XXX[10000]";

/// The unmodified residues of every pseudo sequence (`isSeedPseudoHit_`).
const PSEUDO_PREFIX: &str = "XXX";

/// Minimum number of isotopes extracted when a minimum isotope probability is
/// used instead of a fixed number (source literal `10`, undocumented there).
const PMIN_ISOTOPES: usize = 10;

const FFID_CATEGORY: &str = "FFId_category";
const SEED_FEATURE_ID: &str = "SeedFeatureID";
const OFFSET_PEPTIDE: &str = "OffsetPeptide";

/// Statistics of the ion mobility values of the identifications in an RT
/// region (source protected `IMStats`).
///
/// All values are `-1` when no ion mobility is available.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IMStats {
    /// Median (robust central value, used for the extraction window).
    pub median: f64,
    /// Minimum.
    pub min: f64,
    /// Maximum.
    pub max: f64,
}

impl Default for IMStats {
    fn default() -> Self {
        Self {
            median: -1.0,
            min: -1.0,
            max: -1.0,
        }
    }
}

/// What a run reports in place of the source's log lines.
///
/// No module of this crate logs; every `OPENMS_LOG_INFO`/`OPENMS_LOG_WARN` line
/// the source writes is available as text in
/// [`messages`](Self::messages), and the numbers behind the summary lines are
/// fields.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FeatureFinderIdentificationReport {
    /// Seeds added because no identification overlapped them in RT and m/z.
    pub seeds_added: usize,
    /// Offset peptides added (`add_mass_offset_peptides`).
    pub offset_peptides_added: usize,
    /// Feature candidates before the annotation and filtering.
    pub candidates: usize,
    /// Features left after filtering (`postProcess_`).
    pub features_after_filtering: usize,
    /// Seed-derived features removed for a detected apex too far from the seed.
    pub seed_apex_features_removed: usize,
    /// Features removed because FAIMS merging combined them.
    pub faims_features_merged: usize,
    /// The FAIMS compensation voltages processed, in order.
    pub faims_voltages: Vec<f64>,
    /// Distinct peptides identified (`n_peptides_`).
    pub peptides_identified: usize,
    /// Distinct peptides with a feature of positive intensity.
    pub peptides_quantified: usize,
    /// Peptides without a feature.
    pub peptides_missing: usize,
    /// The report of the elution model fit, when one ran.
    pub elution_model: Option<ElutionModelFitReport>,
    /// The log lines, in order, at most [`MAX_REPORT_MESSAGES`].
    pub messages: Vec<String>,
    /// Lines beyond [`MAX_REPORT_MESSAGES`].
    pub dropped_messages: usize,
}

impl FeatureFinderIdentificationReport {
    fn note(&mut self, line: impl Into<String>) {
        if self.messages.len() < MAX_REPORT_MESSAGES {
            self.messages.push(line.into());
        } else {
            self.dropped_messages = self.dropped_messages.saturating_add(1);
        }
    }
}

/// The summary statistics of a feature map, as
/// [`FeatureFinderIdentificationAlgorithm::statistics`] returns them.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FeatureFinderIdentificationStatistics {
    /// Distinct peptides identified (modifications counted).
    pub identified: usize,
    /// Distinct peptides with a feature of positive intensity.
    pub quantified: usize,
    /// Peptides without a feature; never negative.
    pub missing: usize,
    /// The text the source writes to the info log.
    pub text: String,
}

/// Retention times (not necessarily unique) to the identifications at them,
/// as the source `std::multimap<double, PeptideIdentification*>`.
///
/// Sorted by retention time; entries of equal time keep insertion order, as
/// `multimap::emplace` appends at the upper bound. The second value is an index
/// into the identification list of the run. Times are ordered by
/// [`f64::total_cmp`], so a NaN time (an identification without RT) sorts last
/// instead of corrupting the order.
type RtMap = Vec<(f64, usize)>;

/// Charge to [`RtMap`] (source `ChargeMap`).
type ChargeMap = BTreeMap<i32, RtMap>;

/// Append at the upper bound of equal keys (`multimap::emplace`).
fn rt_map_insert(map: &mut RtMap, rt: f64, index: usize) {
    let at = map.partition_point(|entry| entry.0.total_cmp(&rt) != Ordering::Greater);
    map.insert(at, (rt, index));
}

/// A peptide sequence ordered as the source `AASequence::operator<` orders the
/// keys of its `std::map<AASequence, ...>`.
///
/// The source orders first by the number of residues, then by the terminal
/// modifications (none before any), then residue by residue by one-letter code
/// and modification, and last by the C-terminal modification. A modification is
/// compared by address in the source for residues and by identifier for the
/// terminals; here both are compared by the registry record identifier and then
/// the full identifier, which is the allocation order of the registry. The
/// order only decides in which order peptides that differ in the identity of
/// their modifications are processed.
#[derive(Clone, Debug)]
struct SequenceKey(AASequence);

impl PartialEq for SequenceKey {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for SequenceKey {}
impl PartialOrd for SequenceKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for SequenceKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.source_order(&other.0)
    }
}

/// Mapping sequence to charge to retention time to identification (source
/// `PeptideMap`).
type PeptideMap = BTreeMap<SequenceKey, ChargeMap>;

/// A peptide reference to the identifications mapped to it (source
/// `PeptideRefRTMap`).
type PeptideRefRtMap = BTreeMap<String, RtMap>;

/// Region in RT in which a peptide elutes (source protected `RTRegion`).
struct RtRegion {
    start: f64,
    end: f64,
    ids: ChargeMap,
}

/// The members that follow the parameters (source `updateMembers_`).
#[derive(Clone, Debug, PartialEq)]
struct Settings {
    batch_size: usize,
    rt_window: f64,
    mz_window: f64,
    mz_window_ppm: bool,
    mapping_tolerance: f64,
    isotope_pmin: f64,
    n_isotopes: usize,
    peak_width: f64,
    min_peak_width: f64,
    signal_to_noise: f64,
    elution_model: String,
    candidates_out: String,
    debug_level: usize,
    quantify_decoys: bool,
    use_psm_cutoff: bool,
    psm_score_cutoff: f64,
    add_mass_offset_peptides: f64,
    seed_apex_rt_tolerance: f64,
}

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

fn float_param(param: &Param, key: &str) -> Result<f64> {
    param.value(key)?.to_f64()
}

fn size_param(param: &Param, key: &str) -> Result<usize> {
    let value = param.value(key)?.to_i64()?;
    usize::try_from(value).map_err(|_| bad(format!("parameter '{key}' must not be negative")))
}

impl Settings {
    fn from_param(param: &Param) -> Result<Self> {
        let min_psm = param.value("min_psm_cutoff")?;
        let (use_psm_cutoff, psm_score_cutoff) = match min_psm {
            ParamValue::String(text) if text == "none" => (false, 0.0),
            ParamValue::String(text) => {
                let cutoff = text.trim().parse::<f64>().map_err(|_| {
                    bad(format!(
                        "'min_psm_cutoff' must be 'none' or a number, found '{text}'"
                    ))
                })?;
                (true, cutoff)
            }
            other => (true, other.to_f64()?),
        };
        let mz_window = float_param(param, "extract:mz_window")?;
        Ok(Self {
            batch_size: size_param(param, "extract:batch_size")?,
            rt_window: float_param(param, "extract:rt_window")?,
            mz_window,
            mz_window_ppm: mz_window >= 1.0,
            mapping_tolerance: float_param(param, "detect:mapping_tolerance")?,
            isotope_pmin: float_param(param, "extract:isotope_pmin")?,
            n_isotopes: size_param(param, "extract:n_isotopes")?,
            peak_width: float_param(param, "detect:peak_width")?,
            min_peak_width: float_param(param, "detect:min_peak_width")?,
            signal_to_noise: float_param(param, "detect:signal_to_noise")?,
            elution_model: param.value("model:type")?.as_str()?.to_owned(),
            candidates_out: param.value("candidates_out")?.as_str()?.to_owned(),
            debug_level: size_param(param, "debug")?,
            quantify_decoys: param.value("quantify_decoys")?.to_bool()?,
            use_psm_cutoff,
            psm_score_cutoff,
            add_mass_offset_peptides: float_param(param, "add_mass_offset_peptides")?,
            seed_apex_rt_tolerance: float_param(param, "seed_apex_rt_tolerance")?,
        })
    }
}

/// A numeric meta value, a numeric string, or `default` when absent (the
/// source `getMetaValue(key, default)` read as a double).
fn meta_f64_or(meta: &MetaInfo, key: &str, default: f64) -> Result<f64> {
    match meta.get(key) {
        None => Ok(default),
        Some(value) if value.is_empty() => Ok(default),
        Some(value) => match value.as_f64() {
            Ok(number) => Ok(number),
            Err(_) => value
                .as_str()
                .ok()
                .and_then(|text| text.trim().parse::<f64>().ok())
                .ok_or_else(|| bad(format!("meta value '{key}' is not a number"))),
        },
    }
}

/// `StringUtils::toStr(DataValue)`: the text of a meta value, empty when absent.
fn meta_text(meta: &MetaInfo, key: &str) -> String {
    match meta.get(key) {
        None => String::new(),
        Some(value) => match value.as_str() {
            Ok(text) => text.to_owned(),
            Err(_) => value.to_string(),
        },
    }
}

/// `PeptideCompare`: sequence text, then charge, then RT of the first hits.
fn peptide_compare(left: &PeptideIdentification, right: &PeptideIdentification) -> Ordering {
    let key = |p: &PeptideIdentification| -> (String, i32, f64) {
        match p.hits.first() {
            Some(hit) => (
                hit.sequence.to_string(),
                hit.charge,
                p.rt.unwrap_or(f64::NAN),
            ),
            None => (String::new(), 0, p.rt.unwrap_or(f64::NAN)),
        }
    };
    let (seq1, charge1, rt1) = key(left);
    let (seq2, charge2, rt2) = key(right);
    seq1.cmp(&seq2)
        .then(charge1.cmp(&charge2))
        .then(rt1.total_cmp(&rt2))
}

/// `FeatureCompare`: peptide reference text, then RT.
fn feature_compare(left: &Feature, right: &Feature) -> Ordering {
    let ref1 = meta_text(&left.metadata, "PeptideRef");
    let ref2 = meta_text(&right.metadata, "PeptideRef");
    ref1.cmp(&ref2).then(left.rt.total_cmp(&right.rt))
}

fn is_seed_pseudo_hit(hit: &PeptideHit) -> bool {
    hit.sequence.as_str().starts_with(PSEUDO_PREFIX)
}

/// ID-guided MS1 feature finder, the algorithm behind the TOPP tool
/// `FeatureFinderIdentification` (source
/// `FeatureFinderIdentificationAlgorithm`).
///
/// Create it with [`new`](Self::new) or [`from_param`](Self::from_param), hand
/// it the MS1 data with [`set_ms_data`](Self::set_ms_data), and call
/// [`run`](Self::run). The instance is meant for one run: like the source it
/// keeps the quantities of a run (isotope probabilities, ion mobility
/// statistics) in members, and the elution window members are updated in place.
///
/// # FAIMS data
///
/// When the MS data carry several FAIMS compensation voltages, [`run`](Self::run)
/// splits the data, runs one fresh instance per voltage on the identifications
/// of that voltage (those without a `FAIMS_CV` meta value take part in every
/// group), annotates the features with `FAIMS_CV` and merges features of the
/// same analyte. [`library`](Self::library) is empty then.
///
/// # Differences from the source
///
/// See the support document; the ones that change results are listed at the
/// items they affect ([`run`](Self::run) for the pre-flight checks and the
/// handling of `min_psm_cutoff`).
#[derive(Clone)]
pub struct FeatureFinderIdentificationAlgorithm {
    param: Param,
    settings: Settings,
    ms_data: MSExperiment,
    chrom_data: MSExperiment,
    library: TargetedExperiment,
    output_library: TargetedExperiment,
    peptide_map: PeptideMap,
    n_peptides: usize,
    unassigned_ids: Vec<PeptideIdentification>,
    isotope_probs: BTreeMap<String, f64>,
    im_stats: BTreeMap<String, IMStats>,
    global_im_stats: IMStats,
    feat_finder: MRMFeatureFinderScoring,
    progress: ProgressLogger,
    ids: Arc<Mutex<UniqueIdGenerator>>,
    report: FeatureFinderIdentificationReport,
    capture_candidates: bool,
    captured_candidates: Option<FeatureMap>,
}

impl std::fmt::Debug for FeatureFinderIdentificationAlgorithm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FeatureFinderIdentificationAlgorithm")
            .field("n_peptides", &self.n_peptides)
            .field("rt_window", &self.settings.rt_window)
            .finish_non_exhaustive()
    }
}

fn float_default(
    param: &mut Param,
    key: &str,
    value: f64,
    description: &str,
    advanced: bool,
    min: Option<f64>,
    max: Option<f64>,
) -> Result<()> {
    let tags: Vec<String> = if advanced {
        vec!["advanced".to_owned()]
    } else {
        Vec::new()
    };
    param.set_value(key, ParamValue::Float(value), description, &tags)?;
    if let Some(min) = min {
        param.set_min_float(key, min)?;
    }
    if let Some(max) = max {
        param.set_max_float(key, max)?;
    }
    Ok(())
}

fn int_default(
    param: &mut Param,
    key: &str,
    value: i64,
    description: &str,
    advanced: bool,
    min: i32,
) -> Result<()> {
    let tags: Vec<String> = if advanced {
        vec!["advanced".to_owned()]
    } else {
        Vec::new()
    };
    param.set_value(key, ParamValue::Integer(value), description, &tags)?;
    param.set_min_int(key, min)?;
    Ok(())
}

fn string_default(
    param: &mut Param,
    key: &str,
    value: &str,
    description: &str,
    valid: &[&str],
    tags: &[&str],
) -> Result<()> {
    let tags: Vec<String> = tags.iter().map(|t| (*t).to_owned()).collect();
    param.set_value(
        key,
        ParamValue::String(value.to_owned()),
        description,
        &tags,
    )?;
    if !valid.is_empty() {
        let valid: Vec<String> = valid.iter().map(|v| (*v).to_owned()).collect();
        param.set_valid_strings(key, &valid)?;
    }
    Ok(())
}

impl FeatureFinderIdentificationAlgorithm {
    /// The handler name, `FeatureFinderIdentificationAlgorithm`.
    pub const NAME: &'static str = NAME;

    /// An algorithm with the default parameters (source default constructor).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure, which the constant tree cannot
    /// trigger.
    pub fn new() -> Result<Self> {
        Self::from_param(&Self::defaults()?)
    }

    /// The default parameters (source `getDefaults`), with the source's
    /// descriptions, restrictions and `advanced` tags in declaration order.
    ///
    /// The `model:` section is the parameter set of
    /// [`ElutionModelFitter`](crate::analysis::featurefinder::elution_model_fitter::ElutionModelFitter)
    /// without `asymmetric`, which the `model:type` selects.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        use crate::analysis::featurefinder::elution_model_fitter::ElutionModelFitter;
        let mut p = Param::new();
        string_default(
            &mut p,
            "candidates_out",
            "",
            "Optional output file with feature candidates.",
            &[],
            &["output file"],
        )?;
        int_default(
            &mut p,
            "debug",
            0,
            "Debug level for feature detection.",
            true,
            0,
        )?;
        int_default(
            &mut p,
            "extract:batch_size",
            5000,
            "Nr of peptides used in each batch of chromatogram extraction. Smaller values decrease memory usage but increase runtime.",
            false,
            1,
        )?;
        float_default(
            &mut p,
            "extract:mz_window",
            10.0,
            "m/z window size for chromatogram extraction (unit: ppm if 1 or greater, else Da/Th)",
            false,
            Some(0.0),
            None,
        )?;
        float_default(
            &mut p,
            "extract:IM_window",
            0.06,
            "Ion mobility (IM) window for chromatogram extraction in the IM dimension. Set to 0.0 to disable IM filtering (even if data contains IM information). The window is applied as +/- IM_window/2 around the median IM value of identified peptides. This parameter is automatically ignored if the input data does not contain IM information (determined via IMTypes::determineIMFormat). Currently only concatenated IM format is supported. Typical values: 0.05-0.10 for TIMS data (1/K0 units), 10-50 for CCS data (square angstroms), 3-5 for FAIMS data (compensation voltage).Note: IM values are calculated per peptide/charge/RT-region, using the median of all identifications in that region for robustness. The median, min, and max IM values are propagated to output features as meta-values (IM_median, IM_min, IM_max) for quality control.",
            false,
            Some(0.0),
            None,
        )?;
        int_default(
            &mut p,
            "extract:n_isotopes",
            2,
            "Number of isotopes to include in each peptide assay.",
            false,
            2,
        )?;
        float_default(
            &mut p,
            "extract:isotope_pmin",
            0.0,
            "Minimum probability for an isotope to be included in the assay for a peptide. If set, this parameter takes precedence over 'extract:n_isotopes'.",
            true,
            Some(0.0),
            Some(1.0),
        )?;
        float_default(
            &mut p,
            "extract:rt_window",
            0.0,
            "RT window size (in sec.) for chromatogram extraction. If not set, it is derived from 'detect:peak_width' and 'detect:mapping_tolerance'.",
            true,
            Some(0.0),
            None,
        )?;
        p.set_section_description("extract", "Parameters for ion chromatogram extraction")?;

        float_default(
            &mut p,
            "detect:peak_width",
            60.0,
            "Expected elution peak width in seconds, for smoothing (Gauss filter). Also determines the RT extration window, unless set explicitly via 'extract:rt_window'.",
            false,
            Some(0.0),
            None,
        )?;
        float_default(
            &mut p,
            "detect:min_peak_width",
            0.2,
            "Minimum elution peak width. Absolute value in seconds if 1 or greater, else relative to 'peak_width'.",
            true,
            Some(0.0),
            None,
        )?;
        float_default(
            &mut p,
            "detect:signal_to_noise",
            0.8,
            "Signal-to-noise threshold for OpenSWATH feature detection",
            true,
            Some(0.1),
            None,
        )?;
        float_default(
            &mut p,
            "detect:mapping_tolerance",
            0.0,
            "RT tolerance (plus/minus) for mapping peptide IDs to features. Absolute value in seconds if 1 or greater, else relative to the RT span of the feature.",
            false,
            Some(0.0),
            None,
        )?;
        p.set_section_description(
            "detect",
            "Parameters for detecting features in extracted ion chromatograms",
        )?;

        string_default(
            &mut p,
            "quantify_decoys",
            "false",
            "Whether decoy peptides should be quantified (true) or skipped (false).",
            &["true", "false"],
            &[],
        )?;
        string_default(
            &mut p,
            "min_psm_cutoff",
            "none",
            "Minimum score for the best PSM of a spectrum to be used as seed. Use 'none' for no cutoff.",
            &[],
            &[],
        )?;
        float_default(
            &mut p,
            "add_mass_offset_peptides",
            0.0,
            "If for every peptide (or seed) also an offset peptide is extracted (true). Can be used to downstream to determine MBR false transfer rates. (0.0 = disabled)",
            false,
            Some(0.0),
            None,
        )?;
        float_default(
            &mut p,
            "seed_apex_rt_tolerance",
            5.0,
            "Maximum allowed RT deviation (in seconds) between a seed's apex and the detected feature's apex. Seed-derived features whose detected apex deviates more than this value from the original seed apex are removed during filtering. Useful to discard unreliable seed extractions where the picked peak is far from the seed location. (0 = disabled)",
            true,
            Some(0.0),
            None,
        )?;

        // parameters for model fitting (via ElutionModelFitter)
        string_default(
            &mut p,
            "model:type",
            "symmetric",
            "Type of elution model to fit to features",
            &["symmetric", "asymmetric", "none"],
            &[],
        )?;
        p.insert("model:", &ElutionModelFitter::defaults()?)?;
        p.remove("model:asymmetric")?;
        p.set_section_description("model", "Parameters for fitting elution models to features")?;

        int_default(
            &mut p,
            "EMGScoring:max_iteration",
            100,
            "Maximum number of iterations for EMG fitting.",
            false,
            1,
        )?;
        string_default(
            &mut p,
            "EMGScoring:init_mom",
            "true",
            "Alternative initial parameters for fitting through method of moments.",
            &["true", "false"],
            &[],
        )?;
        p.set_section_description(
            "EMGScoring",
            "Parameters for fitting exp. mod. Gaussians to mass traces.",
        )?;

        string_default(
            &mut p,
            "faims:merge_features",
            "true",
            "For FAIMS data with multiple compensation voltages: Merge features that represent the same analyte detected at different CVs. Features are merged if they have the same charge and are within 5 seconds RT and 0.05 Da m/z. Intensities are summed.",
            &["true", "false"],
            &[],
        )?;
        p.set_section_description("faims", "Parameters for FAIMS data processing")?;
        Ok(p)
    }

    /// Build an algorithm from parameters, as source `setParameters` followed
    /// by `updateMembers_`; missing parameters take their defaults.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong type or
    /// violates its restriction (the source `Exception::InvalidParameter`), or
    /// when `min_psm_cutoff` is neither `none` nor a number.
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`from_param`](Self::from_param), also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`from_param`](Self::from_param).
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        let ((param, settings), warnings) = handler
            .set_parameters_with(parameters, |complete| {
                Ok((complete.clone(), Settings::from_param(complete)?))
            })?;
        let ids = Arc::new(Mutex::new(UniqueIdGenerator::new()));
        let mut feat_finder = MRMFeatureFinderScoring::new()?;
        feat_finder.set_shared_unique_id_generator(Arc::clone(&ids));
        Ok((
            Self {
                param,
                settings,
                ms_data: MSExperiment::new(),
                chrom_data: MSExperiment::new(),
                library: TargetedExperiment::new(),
                output_library: TargetedExperiment::new(),
                peptide_map: PeptideMap::new(),
                n_peptides: 0,
                unassigned_ids: Vec::new(),
                isotope_probs: BTreeMap::new(),
                im_stats: BTreeMap::new(),
                global_im_stats: IMStats::default(),
                feat_finder,
                progress: ProgressLogger::new(),
                ids,
                report: FeatureFinderIdentificationReport::default(),
                capture_candidates: false,
                captured_candidates: None,
            },
            warnings,
        ))
    }

    /// Replace the parameters (source `setParameters`); the members follow.
    ///
    /// The MS data and the other state of the instance are kept. Returns the
    /// unknown-parameter warnings.
    ///
    /// # Errors
    ///
    /// As [`from_param`](Self::from_param); the algorithm is unchanged on error.
    pub fn set_parameters(&mut self, parameters: &Param) -> Result<Vec<String>> {
        let (fresh, warnings) = Self::from_param_with_warnings(parameters)?;
        self.param = fresh.param;
        self.settings = fresh.settings;
        Ok(warnings)
    }

    /// The current parameters (source `getParameters`).
    pub fn parameters(&self) -> &Param {
        &self.param
    }

    /// Draw the unique IDs of this algorithm (offset peptides, features and the
    /// feature map) from `generator`.
    ///
    /// The source draws them from the process-wide `UniqueIdGenerator`; a
    /// caller that needs reproducible IDs seeds the generator and passes it
    /// here. The default generator is seeded from the clock.
    pub fn set_unique_id_generator(&mut self, generator: UniqueIdGenerator) {
        self.ids = Arc::new(Mutex::new(generator));
        self.feat_finder
            .set_shared_unique_id_generator(Arc::clone(&self.ids));
    }

    /// Share the generator of another instance, so that the IDs of several
    /// groups come from one stream (the source's process-wide generator).
    fn share_unique_id_generator(&mut self, ids: &Arc<Mutex<UniqueIdGenerator>>) {
        self.ids = Arc::clone(ids);
        self.feat_finder
            .set_shared_unique_id_generator(Arc::clone(&self.ids));
    }

    /// The cached MS1 data (source `getMSData`). After a run the spectra are
    /// gone: the source releases them once the chromatograms are extracted.
    pub fn ms_data(&self) -> &MSExperiment {
        &self.ms_data
    }

    /// Mutable access to the cached MS data (source non-const `getMSData`).
    pub fn ms_data_mut(&mut self) -> &mut MSExperiment {
        &mut self.ms_data
    }

    /// Set the MS data, keeping only the MS1 spectra (source `setMSData`,
    /// both overloads: copying and moving).
    pub fn set_ms_data(&mut self, ms_data: MSExperiment) {
        self.ms_data = ms_data;
        keep_ms1(&mut self.ms_data);
    }

    /// The extracted chromatograms (source `getChromatograms`).
    ///
    /// The source clears them after every batch of the run, so they are empty
    /// once [`run`](Self::run) has returned; this port does the same.
    pub fn chromatograms(&self) -> &MSExperiment {
        &self.chrom_data
    }

    /// The assay library used and produced by the last run (source
    /// `getLibrary`): peptides, proteins and transitions accumulated over the
    /// batches. Empty for FAIMS data with several voltages, where every group
    /// has its own library.
    // the source's `getLibrary` returns `output_library_`; the working library of
    // a batch is `library`
    #[allow(clippy::misnamed_getters)]
    pub fn library(&self) -> &TargetedExperiment {
        &self.output_library
    }

    /// The progress logger of the algorithm (source `getProgressLogger`):
    /// configure its log type before [`run`](Self::run). The run reports one
    /// step per batch of peptides, labelled "Creating assay library and
    /// extracting chromatograms".
    pub fn progress_logger(&self) -> &ProgressLogger {
        &self.progress
    }

    /// Mutable access to the progress logger (source non-const
    /// `getProgressLogger`).
    pub fn progress_logger_mut(&mut self) -> &mut ProgressLogger {
        &mut self.progress
    }

    /// Keep the feature candidates of the next run in memory (native addition;
    /// default off).
    ///
    /// The source can only hand the candidates over as the file named by
    /// `candidates_out`, written after the annotation and before the filtering.
    /// With this switch the same map, with the same content, is kept and
    /// returned by [`candidates`](Self::candidates), without the round trip
    /// through a file (which loses nothing in the source; the native writer
    /// prints a `float`-valued meta value with the digits of a `float`).
    pub fn set_capture_candidates(&mut self, capture: bool) {
        self.capture_candidates = capture;
        if !capture {
            self.captured_candidates = None;
        }
    }

    /// The feature candidates of the last run, before the filtering, when
    /// [`set_capture_candidates`](Self::set_capture_candidates) was on.
    pub fn candidates(&self) -> Option<&FeatureMap> {
        self.captured_candidates.as_ref()
    }

    /// The report of the last run (what the source writes to its log streams).
    pub fn report(&self) -> &FeatureFinderIdentificationReport {
        &self.report
    }

    /// The RT window in seconds in effect (member `rt_window_`): the parameter
    /// `extract:rt_window`, or the value derived from the peak width after a run.
    pub fn rt_window(&self) -> f64 {
        self.settings.rt_window
    }

    /// The minimum peak width in seconds in effect (member `min_peak_width_`).
    ///
    /// A run turns a relative minimum width (below one) into seconds in place,
    /// as the source does.
    pub fn min_peak_width(&self) -> f64 {
        self.settings.min_peak_width
    }

    /// The number of distinct peptides, charge states aside (member
    /// `n_peptides_`), as used by the summary statistics.
    pub fn peptide_count(&self) -> usize {
        self.n_peptides
    }

    /// Set the peptide count the summary statistics report.
    ///
    /// The source member is protected and its class test reaches it through a
    /// subclass; this is the same hook.
    pub fn set_peptide_count(&mut self, count: usize) {
        self.n_peptides = count;
    }

    /// The summary statistics of a feature map (source `statistics_`): the
    /// distinct peptides identified, with a feature, and without one.
    ///
    /// A peptide counts as quantified when a feature with positive intensity
    /// carries it as the first hit of its first identification. The number of
    /// missing peptides never wraps: it is zero when more peptides were
    /// quantified than identified (the source printed a wrapped unsigned value
    /// until the defect fixed in its class test `[EXTRA]` section).
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] for a feature of positive intensity
    /// without an identification or hit (undefined behaviour in the source).
    pub fn statistics(
        &self,
        features: &FeatureMap,
    ) -> Result<FeatureFinderIdentificationStatistics> {
        let mut quantified: BTreeSet<SequenceKey> = BTreeSet::new();
        for feature in &features.features {
            if feature.intensity > 0.0 {
                let hit = feature
                    .peptide_identifications
                    .first()
                    .and_then(|id| id.hits.first())
                    .ok_or_else(|| {
                        Error::MissingInformation(
                            "a feature with positive intensity has no peptide hit".into(),
                        )
                    })?;
                quantified.insert(SequenceKey(hit.sequence.clone()));
            }
        }
        let identified = self.n_peptides;
        let n_quantified = quantified.len();
        let n_missing = identified.saturating_sub(n_quantified);
        let text = format!(
            "\nSummary statistics (counting distinct peptides including PTMs):\n{identified} peptides identified\n{n_quantified} peptides with features\n{n_missing} peptides without features\n"
        );
        Ok(FeatureFinderIdentificationStatistics {
            identified,
            quantified: n_quantified,
            missing: n_missing,
            text,
        })
    }

    // ------------------------------------------------------------------
    // run
    // ------------------------------------------------------------------

    /// Run the pipeline: detect and quantify a feature for every identified
    /// peptide (and seed) in the MS data (source `run`).
    ///
    /// For FAIMS data one run per compensation voltage is made and the results
    /// are merged (see the type documentation).
    ///
    /// # Arguments
    ///
    /// * `peptides` - the peptide identifications. The local copy is reduced to
    ///   the best hit per identification and annotated; the caller's value is
    ///   not observable.
    /// * `proteins` - the protein identifications; they become the protein
    ///   identifications of `features`.
    /// * `features` - receives the quantified features. Pre-existing content is
    ///   replaced for FAIMS data and extended otherwise. It is left unchanged
    ///   when this returns an error.
    /// * `seeds` - optional pre-detected features of an untargeted feature
    ///   finder; one that overlaps no identification in RT and m/z becomes an
    ///   extra target.
    /// * `spectra_file` - the path used for the `spectra_data` annotation when
    ///   the MS data carry no usable run path.
    ///
    /// # Errors
    ///
    /// * [`Error::InvalidValue`] without spectra (the source's
    ///   `Exception::IllegalArgument`), when `candidates_out` is set and no
    ///   writer for its format is available, or when more than [`MAX_TARGETS`]
    ///   targets would be created.
    /// * [`Error::Unsupported`] when `candidates_out` names a feature Parquet
    ///   file or the `featurexml` feature is disabled.
    /// * [`Error::MissingInformation`] for an identification of the library
    ///   without peptide hit where the source has undefined behaviour, and
    ///   [`Error::InvalidValue`] for an unresolvable sequence or a hit without
    ///   charge.
    /// * Every error of chromatogram extraction, picking, scoring and elution
    ///   model fitting.
    ///
    /// # Differences from the source
    ///
    /// * The preconditions above are checked before anything is mutated; the
    ///   source fails somewhere during the run, or writes out of bounds.
    /// * `min_psm_cutoff` is read as the number its text spells. The source
    ///   converts the string parameter through the `double` member of the value
    ///   union, which holds the address of the string, so a cutoff other than
    ///   `none` acts as a number near zero there.
    /// * The `debug` level 668 library dump (`debug.traml`) is not written: the
    ///   TraML writer is not ported yet. The run itself is unaffected.
    pub fn run(
        &mut self,
        peptides: Vec<PeptideIdentification>,
        proteins: &[ProteinIdentification],
        features: &mut FeatureMap,
        seeds: &FeatureMap,
        spectra_file: &str,
    ) -> Result<()> {
        // Without spectra there is nothing to extract chromatograms from.
        if self.ms_data.spectra.is_empty() {
            return Err(bad(
                "No spectra in the MS data passed to FeatureFinderIdentification. Set the MS data via set_ms_data() before calling run(); note that an experiment consumed by another algorithm may have been left empty.",
            ));
        }
        let doubled = if self.settings.add_mass_offset_peptides > 0.0 {
            2
        } else {
            1
        };
        peptides
            .len()
            .checked_add(seeds.features.len())
            .and_then(|n| n.checked_mul(doubled))
            .filter(|&n| n <= MAX_TARGETS)
            .ok_or_else(|| {
                bad(format!(
                    "more than {MAX_TARGETS} identifications and seeds were given"
                ))
            })?;
        self.check_candidates_output()?;

        let split = ImDataConverter::split_by_faims_cv(&mut self.ms_data)?;
        if !split.has_faims() {
            // Non-FAIMS data: restore the data and run directly.
            let mut groups = split.groups.into_iter();
            if let Some(group) = groups.next() {
                self.ms_data = group.experiment;
            }
            let mut staged = features.clone();
            self.run_single_group(peptides, proteins, &mut staged, seeds, spectra_file)?;
            *features = staged;
            return Ok(());
        }

        // FAIMS data: process each compensation voltage separately.
        let mut staged = features.clone();
        staged.clear(true);
        self.chrom_data = MSExperiment::new();
        self.output_library.clear(true);
        self.report = FeatureFinderIdentificationReport::default();
        let mut first_group = true;
        let group_count = split.groups.len();
        self.report.note(format!(
            "FAIMS data detected with {group_count} compensation voltage(s)."
        ));
        for group in split.groups {
            let group_cv = group.key.volts();
            self.report
                .note(format!("Processing FAIMS CV group: {group_cv} V"));
            let peptides_cv = FaimsHelper::filter_peptides_by_faims_cv(&peptides, group_cv, 0.01)?;
            if peptides_cv.is_empty() {
                self.report.note(format!(
                    "No peptide IDs for FAIMS CV {group_cv} V. Skipping."
                ));
                continue;
            }
            // a fresh instance per group, with the same parameters
            let mut sub = Self::from_param(&self.param)?;
            sub.progress.set_log_type(self.progress.log_type());
            sub.share_unique_id_generator(&self.ids);
            sub.set_ms_data(group.experiment);
            let mut features_cv = FeatureMap::new();
            sub.run_single_group(peptides_cv, proteins, &mut features_cv, seeds, spectra_file)?;

            let cv_value = MetaValue::try_from(group_cv)?;
            for mut feature in features_cv.features {
                feature
                    .metadata
                    .insert(FAIMS_CV.to_owned(), cv_value.clone());
                staged.features.push(feature);
            }
            for mut id in features_cv.unassigned_peptide_identifications {
                id.metadata.insert(FAIMS_CV.to_owned(), cv_value.clone());
                staged.unassigned_peptide_identifications.push(id);
            }
            if first_group {
                staged.protein_identifications = features_cv.protein_identifications;
            }
            for mut chromatogram in sub.chrom_data.chromatograms {
                chromatogram
                    .metadata
                    .insert(FAIMS_CV.to_owned(), cv_value.clone());
                self.chrom_data.chromatograms.push(chromatogram);
            }
            self.report.faims_voltages.push(group_cv);
            first_group = false;
        }
        self.report.note(
            "Warning: Library output is not available for multi-FAIMS data. Each FAIMS CV group has its own assay library.",
        );
        let combined = staged.features.len();
        self.report.note(format!(
            "Combined {combined} features from all FAIMS CV groups."
        ));

        // merge features of the same analyte seen at different voltages
        let merge = self.param.value("faims:merge_features")?.to_bool()?;
        let merged = common::merge_faims_features(&mut staged, merge)?;
        self.report.faims_features_merged = merged;
        if merged > 0 {
            self.report.note(format!(
                "Merged FAIMS features: {} -> {} ({merged} features merged)",
                combined,
                combined - merged
            ));
        }

        staged.set_primary_ms_run_path(&[spectra_file.to_owned()])?;
        self.ensure_map_unique_id(&mut staged)?;
        *features = staged;
        Ok(())
    }

    /// Check that the candidates file can be written before anything runs.
    fn check_candidates_output(&self) -> Result<()> {
        let path = &self.settings.candidates_out;
        if path.is_empty() {
            return Ok(());
        }
        let lower = path.to_ascii_lowercase();
        if lower.ends_with(".featureparquet") || lower.ends_with(".parquet") {
            return Err(Error::Unsupported(
                "candidates_out in the feature Parquet format: FeatureParquet writing is not ported yet".into(),
            ));
        }
        if cfg!(not(feature = "featurexml")) {
            return Err(Error::Unsupported(
                "candidates_out requires the 'featurexml' feature".into(),
            ));
        }
        Ok(())
    }

    /// Core of one dataset: not FAIMS, or a single FAIMS group (source
    /// `runSingleGroup_`).
    fn run_single_group(
        &mut self,
        mut peptides: Vec<PeptideIdentification>,
        proteins: &[ProteinIdentification],
        features: &mut FeatureMap,
        seeds: &FeatureMap,
        spectra_file: &str,
    ) -> Result<()> {
        // Clear the output library from any previous run.
        self.output_library.clear(true);
        self.report = FeatureFinderIdentificationReport::default();

        // annotate the mzML file
        features
            .set_primary_ms_run_path_from_experiment(&[spectra_file.to_owned()], &self.ms_data)?;

        // check the IM format
        let im_window = self.param.value("extract:IM_window")?.to_f64()?;
        let im_format = ImTypes::determine_im_format_for_ms_level(&self.ms_data, 1)?;
        let mut has_im = false;
        match im_format {
            IonMobilityFormat::PerPeak => {
                has_im = true;
                if im_window > 0.0 {
                    if let Some(unit) = first_im_unit(&self.ms_data) {
                        if unit == DriftTimeUnit::CollisionCrossSection && im_window < 1.0 {
                            self.report.note(format!(
                                "Warning: Ion mobility data is in CCS units (square angstroms), but IM_window = {im_window} appears to be set for 1/K0 data. For CCS data, consider using larger values (e.g., 10-50)."
                            ));
                        }
                    }
                }
            }
            IonMobilityFormat::None => {}
            _ => self.report.note(
                "Wrong IM format detected. Expecting in concatenated format (float data arrays)",
            ),
        }

        // initialise the feature finder
        self.initialize_feature_finder()?;

        // calculate the RT window if it is not set
        self.settings.rt_window = self.calculate_rt_window();

        //-------------------------------------------------------------
        // prepare the peptide map
        //-------------------------------------------------------------
        self.report.note("Preparing mapping of peptide data...");
        self.peptide_map.clear();
        let reserve = peptides.len().saturating_add(seeds.features.len());
        peptides.reserve(
            if self.settings.add_mass_offset_peptides > 0.0 {
                reserve.saturating_mul(2)
            } else {
                reserve
            }
            .saturating_sub(peptides.len()),
        );
        for index in 0..peptides.len() {
            self.add_peptide_to_map(&mut peptides, index)?;
            peptides[index]
                .metadata
                .insert(FFID_CATEGORY.to_owned(), MetaValue::from("internal"));
        }

        // The global IM statistics are taken before the seeds are added, so
        // that only real identifications are learned from.
        self.calculate_global_im_stats(&peptides)?;

        let pseudo = AASequence::parse(PSEUDO_SEQUENCE)?;
        let seeds_added = self.add_seeds(&mut peptides, seeds, &pseudo)?;
        self.report.seeds_added = seeds_added;
        self.report.note(format!(
            "#Seeds without RT and m/z overlap with identified peptides added: {seeds_added}"
        ));
        if self.settings.add_mass_offset_peptides > 0.0 {
            let offset = self.settings.add_mass_offset_peptides;
            let added = self.add_offset_peptides(&mut peptides, offset, &pseudo)?;
            self.report.offset_peptides_added = added;
            self.report.note(format!(
                "#Offset peptides without RT and m/z overlap with other peptides added: {added}"
            ));
        }
        self.n_peptides = self.peptide_map.len();

        //-------------------------------------------------------------
        // run the feature detection
        //-------------------------------------------------------------
        let shared = Arc::new(std::mem::take(&mut self.ms_data));
        let outcome = self.detect(&shared, &peptides, features, im_window, has_im);
        // the feature finder keeps no reference to the data any more
        self.feat_finder.set_ms1_map(None);
        // the source resets the data here, "not needed anymore"; on an error it
        // is given back, so that the caller can retry
        if outcome.is_err() {
            if let Ok(data) = Arc::try_unwrap(shared) {
                self.ms_data = data;
            }
        }
        let mut ref_rt_map = outcome?;

        // complete the feature annotation
        self.annotate_features(features, &mut ref_rt_map, &peptides)?;

        // sort everything
        features
            .unassigned_peptide_identifications
            .sort_by(peptide_compare);
        features.features.sort_by(feature_compare);

        self.post_process(features)?;
        let statistics = self.statistics(features)?;
        self.report.peptides_identified = statistics.identified;
        self.report.peptides_quantified = statistics.quantified;
        self.report.peptides_missing = statistics.missing;
        self.report.note(statistics.text);

        features.protein_identifications = proteins.to_vec();

        // remove all hits with pseudo IDs (seeds)
        remove_seed_pseudo_ids(features);

        // add back the ignored PSMs
        features
            .unassigned_peptide_identifications
            .append(&mut self.unassigned_ids);

        self.ensure_map_unique_id(features)?;
        Ok(())
    }

    /// The batch loop of `runSingleGroup_`: build the assay library per batch of
    /// peptides, extract the chromatograms and pick the features.
    fn detect(
        &mut self,
        shared: &Arc<MSExperiment>,
        peptides: &[PeptideIdentification],
        features: &mut FeatureMap,
        im_window: f64,
        has_im: bool,
    ) -> Result<PeptideRefRtMap> {
        let spectra = SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(Arc::clone(shared))?;
        // to use the MS1 SWATH scores
        self.feat_finder.set_ms1_map(Some(Arc::clone(&spectra)));
        let settings_spectrum = shared
            .spectra
            .first()
            .ok_or_else(|| bad("the MS data hold no spectrum"))?;

        let keys: Vec<SequenceKey> = self.peptide_map.keys().cloned().collect();
        let chunks = chunk_ranges(keys.len(), self.settings.batch_size);
        let mut ref_rt_map = PeptideRefRtMap::new();

        if self.settings.debug_level >= 668 {
            self.report.note(
                "debug.traml is not written: the TraML writer is not ported (FORMAT/TraMLFile.h).",
            );
        }

        let extraction_im = if has_im && im_window > 0.0 {
            Some(im_window)
        } else {
            None
        };
        // Note: progress only works in non-debug mode, when no logs come in
        // between
        self.progress.start_progress(
            0,
            i64::try_from(chunks.len()).unwrap_or(i64::MAX),
            "Creating assay library and extracting chromatograms",
        )?;
        let mut chunk_count: i64 = 0;
        for (start, end) in chunks {
            // fills the library
            self.create_assay_library(&keys[start..end], &mut ref_rt_map, true, peptides)?;
            let light = convert_targeted_exp(&self.library)?;
            let mut extracted = Vec::new();
            extract_ms1_chromatograms(
                &spectra,
                settings_spectrum,
                &light,
                self.settings.mz_window,
                self.settings.mz_window_ppm,
                extraction_im,
                &mut extracted,
            )?;
            self.chrom_data.chromatograms.extend(extracted);

            // pick the features; this is `pickExperiment` of the OpenMS
            // structures, with the spectra of the run as the DIA map
            let chromatograms = Arc::new(std::mem::take(&mut self.chrom_data));
            let chromatogram_access =
                SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(chromatograms)?;
            let swath = SwathMap {
                sptr: Some(Arc::clone(&spectra)),
                ..SwathMap::default()
            };
            let mut group_map = TransitionGroupMapType::new();
            self.feat_finder.pick_experiment_light(
                &chromatogram_access,
                features,
                &light,
                &TransformationDescription::default(),
                &[swath],
                &mut group_map,
            )?;

            // accumulate the library entries for the output before clearing
            for peptide in self.library.peptides() {
                self.output_library.add_peptide(peptide.clone());
            }
            for protein in self.library.proteins() {
                self.output_library.add_protein(protein.clone());
            }
            for transition in self.library.transitions() {
                self.output_library.add_transition(transition.clone());
            }
            self.library.clear(true);
            // the chromatograms were only a container: the identifications are
            // empty and pickExperiment added empty protein runs with colliding
            // identifiers; the real proteins are added later
            features.protein_identifications.clear();
            chunk_count += 1;
            self.progress.set_progress(chunk_count)?;
        }
        self.progress.end_progress(0)?;
        self.report.candidates = features.features.len();
        self.report.note(format!(
            "Found {} feature candidates in total.",
            features.features.len()
        ));
        Ok(ref_rt_map)
    }

    /// Set the parameters of the OpenSWATH finder (source
    /// `initializeFeatureFinder_`).
    fn initialize_feature_finder(&mut self) -> Result<()> {
        let mut params = self.feat_finder.parameters().clone();
        let max_iteration = self.param.value("EMGScoring:max_iteration")?.clone();
        let init_mom = self.param.value("EMGScoring:init_mom")?.clone();
        let write_hull =
            self.settings.elution_model != "none" || !self.settings.candidates_out.is_empty();
        if write_hull {
            params.set_value(
                "write_convex_hull",
                ParamValue::String("true".into()),
                "",
                &[],
            )?;
        }
        // The relative minimum width becomes seconds, in place, as the source
        // does; a second run of the same instance applies it again when the
        // result is still below one.
        if self.settings.min_peak_width < 1.0 {
            self.settings.min_peak_width *= self.settings.peak_width;
        }
        apply_common_finder_parameters(
            &mut params,
            &max_iteration,
            &init_mom,
            self.settings.peak_width,
            self.settings.min_peak_width,
            self.settings.signal_to_noise,
        )?;
        self.feat_finder.set_parameters(&params)?;
        self.feat_finder.set_strict_flag(false);
        Ok(())
    }

    /// The RT window of the extraction (source `calculateRTWindow_`): the set
    /// one, or twice the peak width plus the mapping tolerance, doubled.
    fn calculate_rt_window(&mut self) -> f64 {
        if self.settings.rt_window != 0.0 {
            return self.settings.rt_window;
        }
        let mut map_tol = self.settings.mapping_tolerance;
        if map_tol < 1.0 {
            map_tol *= 2.0 * self.settings.peak_width; // relative tolerance
        }
        let window = (2.0 * self.settings.peak_width + map_tol) * 2.0;
        self.report
            .note(format!("RT window size calculated as {window} seconds."));
        window
    }

    /// RT bounds with the mapping tolerance applied (source `calculateRTBounds_`).
    fn calculate_rt_bounds(&self, rt_min: f64, rt_max: f64) -> (f64, f64) {
        let tolerance = self.settings.mapping_tolerance;
        if tolerance > 0.0 {
            let mut abs_tol = tolerance;
            if abs_tol < 1.0 {
                abs_tol *= rt_max - rt_min;
            }
            (rt_min - abs_tol, rt_max + abs_tol)
        } else {
            (rt_min, rt_max)
        }
    }

    fn ensure_map_unique_id(&self, features: &mut FeatureMap) -> Result<()> {
        let mut generator = self
            .ids
            .lock()
            .map_err(|_| bad("the unique ID generator is poisoned"))?;
        features.unique_id.ensure_unique_id(&mut generator);
        Ok(())
    }

    fn next_unique_id(&self) -> Result<u64> {
        let mut generator = self
            .ids
            .lock()
            .map_err(|_| bad("the unique ID generator is poisoned"))?;
        Ok(generator.get_unique_id())
    }

    // ------------------------------------------------------------------
    // peptide map
    // ------------------------------------------------------------------

    /// Add an identification to the peptide map (source `addPeptideToMap_`).
    ///
    /// The identification is sorted and reduced to its best hit first. Decoys
    /// (unless `quantify_decoys`) and hits below the score cutoff are set aside
    /// as unassigned. One without hits is dropped.
    fn add_peptide_to_map(
        &mut self,
        peptides: &mut [PeptideIdentification],
        index: usize,
    ) -> Result<()> {
        if peptides[index].hits.is_empty() {
            return Ok(());
        }
        peptides[index].sort()?;
        peptides[index].hits.truncate(1);

        let (decoy, score, higher_better) = {
            let id = &peptides[index];
            (
                id.hits[0].is_decoy()?,
                id.hits[0].score,
                id.higher_score_better,
            )
        };
        // if the decoys are not quantified they are not added to the list
        if !self.settings.quantify_decoys && decoy {
            self.unassigned_ids.push(peptides[index].clone());
            return Ok(());
        }
        if self.settings.use_psm_cutoff {
            let cutoff = self.settings.psm_score_cutoff;
            if (higher_better && score < cutoff) || (!higher_better && score > cutoff) {
                self.unassigned_ids.push(peptides[index].clone());
                return Ok(());
            }
        }

        // A meta value forcefully overwrites the m/z with an external one. To
        // quantify that, a modification that matches the mass difference is
        // introduced at the first unmodified residue.
        let calc_mass = match peptides[index].hits[0].metadata.get("CalcMass") {
            Some(value) => Some(value.as_f64()?),
            None => None,
        };
        if let Some(calc_mass) = calc_mass {
            let charge = peptides[index].hits[0].charge;
            let sequence = peptides[index].hits[0].sequence.clone();
            let diff_mz = calc_mass - source_mz(&sequence, charge)?;
            let diff_mass = diff_mz * f64::from(charge);
            if diff_mass.abs() > 0.01 {
                let mut modified = sequence;
                for position in 0..modified.len() {
                    if modified.residue_modification(position)?.is_some() {
                        continue;
                    }
                    set_modification_by_diff_mono_mass(&mut modified, position, diff_mass)?;
                    break;
                }
                peptides[index].hits[0].sequence = modified;
            }
        }

        let hit = &peptides[index].hits[0];
        let charge = hit.charge;
        let rt = peptides[index].rt.unwrap_or(f64::NAN);
        let key = SequenceKey(hit.sequence.clone());
        let charges = self.peptide_map.entry(key).or_default();
        rt_map_insert(charges.entry(charge).or_default(), rt, index);
        Ok(())
    }

    /// Add the seeds that overlap no identification (source `addSeeds_`).
    ///
    /// A seed is a feature of an untargeted feature finder. It overlaps an
    /// identification when their RT differ by at most half the seed window and
    /// the m/z agree within the extraction tolerance, up to the fifth isotope.
    /// The pseudo identification of a seed is appended to `peptides` and mapped
    /// under the pseudo sequence.
    fn add_seeds(
        &mut self,
        peptides: &mut Vec<PeptideIdentification>,
        seeds: &FeatureMap,
        pseudo: &AASequence,
    ) -> Result<usize> {
        let mut added = 0usize;
        // identifications existing before any seed is added are compared; the
        // source compares against the growing list, seeds included
        for feat in &seeds.features {
            let seed_rt = feat.rt;
            let seed_mz = feat.mz;
            let seed_charge = f64::from(feat.charge);
            let mut already_exists = false;
            for peptide in peptides.iter() {
                let peptide_rt = peptide.rt.unwrap_or(f64::NAN);
                let peptide_mz = peptide.mz.unwrap_or(f64::NAN);
                let tolerance = if self.settings.mz_window_ppm {
                    self.settings.mz_window * 1e-6 * peptide_mz
                } else {
                    self.settings.mz_window
                };
                // `(n / charge) * C13C12`, in this order, as the source
                let within = |shift: f64| {
                    (seed_mz - (shift / seed_charge) * C13C12_MASSDIFF_U - peptide_mz).abs()
                        <= tolerance
                };
                if (seed_rt - peptide_rt).abs() <= SEED_RT_WINDOW / 2.0
                    && ((seed_mz - peptide_mz).abs() <= tolerance
                        || within(1.0)
                        || within(2.0)
                        || within(3.0)
                        || within(4.0)
                        || within(5.0))
                {
                    already_exists = true;
                    break;
                }
            }
            if already_exists {
                continue;
            }
            // The pseudo identification is stored in the list of the run, so
            // that the map entries stay valid for the duration of the run.
            let hit = PeptideHit {
                charge: feat.charge,
                sequence: pseudo.clone(),
                ..PeptideHit::default()
            };
            let mut id = PeptideIdentification {
                hits: vec![hit],
                rt: Some(feat.rt),
                mz: Some(feat.mz),
                ..PeptideIdentification::default()
            };
            id.metadata
                .insert(FFID_CATEGORY.to_owned(), MetaValue::from("internal"));
            id.metadata.insert(
                SEED_FEATURE_ID.to_owned(),
                MetaValue::from(feat.unique_id.to_string()),
            );
            // copy the IM from the feature if present; without it the seed is
            // extracted across the full IM range
            if let Some(im) = feat.metadata.get(IM) {
                id.metadata.insert(IM.to_owned(), im.clone());
            }
            peptides.push(id);
            let index = peptides.len() - 1;
            self.add_peptide_to_map(peptides, index)?;
            added += 1;
        }
        Ok(added)
    }

    /// Add an offset peptide for every peptide and seed (source
    /// `addOffsetPeptides_`): the same RT and charge at a shifted m/z, to
    /// estimate the rate of false transfers.
    fn add_offset_peptides(
        &mut self,
        peptides: &mut Vec<PeptideIdentification>,
        offset: f64,
        pseudo: &AASequence,
    ) -> Result<usize> {
        let mut offsets = Vec::with_capacity(peptides.len());
        for peptide in peptides.iter() {
            // an identification without hits is skipped: the source reads the
            // first hit unconditionally
            let Some(first) = peptide.hits.first() else {
                continue;
            };
            let hit = PeptideHit {
                charge: first.charge,
                sequence: pseudo.clone(),
                ..PeptideHit::default()
            };
            let mut id = PeptideIdentification {
                hits: vec![hit],
                rt: peptide.rt,
                mz: Some(peptide.mz.unwrap_or(f64::NAN) + offset),
                ..PeptideIdentification::default()
            };
            id.metadata
                .insert(FFID_CATEGORY.to_owned(), MetaValue::from("internal"));
            id.metadata
                .insert(OFFSET_PEPTIDE.to_owned(), MetaValue::from("true"));
            // also marked as a seed, to indicate a mass without sequence
            id.metadata.insert(
                SEED_FEATURE_ID.to_owned(),
                MetaValue::from(self.next_unique_id()?.to_string()),
            );
            offsets.push(id);
        }
        let added = offsets.len();
        for id in offsets {
            peptides.push(id);
            let index = peptides.len() - 1;
            self.add_peptide_to_map(peptides, index)?;
        }
        Ok(added)
    }

    // ------------------------------------------------------------------
    // ion mobility
    // ------------------------------------------------------------------

    /// Global ion mobility statistics (source `calculateGlobalIMStats_`): the
    /// range of the data and the median of the identifications' values.
    fn calculate_global_im_stats(&mut self, peptides: &[PeptideIdentification]) -> Result<()> {
        let ranges = self.ms_data.spectrum_range_manager()?;
        let global = ranges.global();
        let (min, max) = match (global.min_mobility(), global.max_mobility()) {
            (Ok(min), Ok(max)) => (min, max),
            _ => {
                // no IM data found in the MS data
                self.global_im_stats = IMStats::default();
                return Ok(());
            }
        };
        self.global_im_stats.min = min;
        self.global_im_stats.max = max;

        // Called before the seeds are added, so only identifications are seen.
        let mut values = Vec::new();
        for charges in self.peptide_map.values() {
            for rt_map in charges.values() {
                for &(_, index) in rt_map {
                    let im = meta_f64_or(&peptides[index].metadata, IM, -1.0)?;
                    if im >= 0.0 {
                        values.push(im);
                    }
                }
            }
        }
        let count = values.len();
        if count > 0 {
            values.sort_by(f64::total_cmp);
            self.global_im_stats.median = median_of_sorted(&values);
            self.report.note(format!(
                "Global IM statistics: median={} (from {count} identifications), min={}, max={} (from MS data range)",
                self.global_im_stats.median, self.global_im_stats.min, self.global_im_stats.max
            ));
        } else {
            // no IDs with IM: the centre of the data range
            self.global_im_stats.median = (min + max) / 2.0;
            self.report.note(format!(
                "Global IM statistics: median={} (center of data range), min={}, max={} (from MS data range)",
                self.global_im_stats.median, self.global_im_stats.min, self.global_im_stats.max
            ));
        }
        Ok(())
    }

    /// The global ion mobility statistics of the last run: the range of the
    /// MS data and the median of the identifications (source member
    /// `global_im_stats_`).
    pub fn global_im_stats(&self) -> IMStats {
        self.global_im_stats
    }

    /// Ion mobility statistics of an RT region (source `getRTRegionIMStats_`):
    /// median, minimum and maximum of the values of the identifications that
    /// carry one; all `-1` when none does.
    fn region_im_stats(
        &mut self,
        region: &RtRegion,
        peptides: &[PeptideIdentification],
    ) -> Result<IMStats> {
        let mut values = Vec::new();
        let mut without = 0usize;
        for rt_map in region.ids.values() {
            for &(_, index) in rt_map {
                let im = meta_f64_or(&peptides[index].metadata, IM, -1.0)?;
                if im < 0.0 {
                    without += 1; // 0.0 is a valid value
                } else {
                    values.push(im);
                }
            }
        }
        if values.is_empty() {
            return Ok(IMStats::default());
        }
        values.sort_by(f64::total_cmp);
        let stats = IMStats {
            min: values[0],
            max: values[values.len() - 1],
            median: median_of_sorted(&values),
        };
        if without > 0 {
            self.report.note(format!(
                "Calculated IM statistics from {} IDs with IM data (skipped {without} IDs without IM annotation)",
                values.len()
            ));
        }
        Ok(stats)
    }

    // ------------------------------------------------------------------
    // assay library
    // ------------------------------------------------------------------

    /// Regions in RT in which a peptide elutes, by clustering the RTs of all its
    /// charge states (source `getRTRegions_`). With `clear_ids` the
    /// identification references of `peptide_data` are removed afterwards (they
    /// are held by the regions then); the charge keys stay.
    fn rt_regions(&self, peptide_data: &mut ChargeMap, clear_ids: bool) -> Vec<RtRegion> {
        let mut rts: Vec<f64> = peptide_data
            .values()
            .flat_map(|rt_map| rt_map.iter().map(|entry| entry.0))
            .collect();
        rts.sort_by(f64::total_cmp);
        let tolerance = self.settings.rt_window / 2.0;
        let mut regions: Vec<RtRegion> = Vec::new();
        for &rt in &rts {
            let new_region = regions.last().is_none_or(|last| last.end < rt - tolerance);
            if new_region {
                regions.push(RtRegion {
                    start: rt - tolerance,
                    end: rt + tolerance,
                    ids: ChargeMap::new(),
                });
            } else if let Some(last) = regions.last_mut() {
                last.end = rt + tolerance;
            }
        }
        // sort the identifications into the regions
        for (&charge, rt_map) in peptide_data.iter_mut() {
            // regions and identifications are sorted by RT, so iterate linearly
            let mut region = 0usize;
            for &(rt, index) in rt_map.iter() {
                // while the RT is larger than the end of the region: next one
                while region + 1 < regions.len() && rt > regions[region].end {
                    region += 1;
                }
                if let Some(target) = regions.get_mut(region) {
                    target.ids.entry(charge).or_default().push((rt, index));
                }
            }
            if clear_ids {
                rt_map.clear();
            }
        }
        regions
    }

    /// Create the assays of a batch of peptides (source `createAssayLibrary_`).
    ///
    /// `ref_rt_map` receives, per assay, the identifications mapped to it. With
    /// `clear_ids` false the identification references stay in the peptide map
    /// (debugging only).
    fn create_assay_library(
        &mut self,
        keys: &[SequenceKey],
        ref_rt_map: &mut PeptideRefRtMap,
        clear_ids: bool,
        peptides: &[PeptideIdentification],
    ) -> Result<()> {
        let mut protein_accessions: BTreeSet<String> = BTreeSet::new();
        let mut seed_count = 0usize;
        let n_isotopes = if self.settings.isotope_pmin > 0.0 {
            PMIN_ISOTOPES
        } else {
            self.settings.n_isotopes
        };
        let generator =
            CoarseIsotopePatternGenerator::new(Some(n_isotopes), CoarseMassMode::Approximate)?
                .with_precision(ProbabilityPrecision::SourceSingle);

        for key in keys {
            let seq = &key.0;
            let mut peptide = Peptide::default();
            // Seeds are stored in the same map under the same fake sequence.
            if seq.as_str().starts_with(PSEUDO_PREFIX) {
                // Forces the scores to consider it an unidentified peptide, and
                // for example use averagine isotopes.
                peptide.sequence = String::new();
                // their retention times are not aggregated: every entry is an
                // assay of its own
                let charges: Vec<(i32, RtMap)> = match self.peptide_map.get(key) {
                    Some(charges) => charges
                        .iter()
                        .map(|(charge, rt_map)| (*charge, rt_map.clone()))
                        .collect(),
                    None => Vec::new(),
                };
                for (charge, rt_map) in charges {
                    for &(rt_key, index) in &rt_map {
                        // seeds do not share a group label, or they would be
                        // combined
                        seed_count += 1;
                        let id = &peptides[index];
                        let mz = id.mz.unwrap_or(f64::NAN);
                        let rt = id.rt.unwrap_or(f64::NAN);
                        let uid = meta_text(&id.metadata, SEED_FEATURE_ID);
                        let peptide_id =
                            format!("{}[{uid}][{seed_count}]/{charge}", peptide.sequence);
                        peptide.set_charge_state(charge);
                        peptide.id = peptide_id.clone();
                        peptide.protein_refs = vec!["not_available".to_owned()];
                        peptide.set_peptide_group_label(peptide_id.clone());

                        // an entry of the output map
                        let peptide_ids = ref_rt_map.entry(peptide_id.clone()).or_default();

                        // isotope distribution of the seed
                        let mass = mz * f64::from(charge) - f64::from(charge) * PROTON_MASS_U;
                        let mut iso_dist = generator.estimate_from_peptide_weight(mass)?;
                        trim_and_renormalize(&mut iso_dist, self.settings.isotope_pmin)?;

                        let rt_tolerance = SEED_RT_WINDOW / 2.0;
                        // beginning and end of the RT region: only one entry
                        peptide.rts.clear();
                        add_peptide_rt(&mut peptide, rt - rt_tolerance);
                        add_peptide_rt(&mut peptide, rt + rt_tolerance);

                        // IM of the seed if annotated, else the drift time
                        // stays at its default and the full range is extracted
                        let seed_im = meta_f64_or(&id.metadata, IM, -1.0)?;
                        if seed_im >= 0.0 {
                            peptide.set_drift_time(seed_im);
                            self.im_stats.insert(
                                peptide.id.clone(),
                                IMStats {
                                    median: seed_im,
                                    min: seed_im,
                                    max: seed_im,
                                },
                            );
                        } else {
                            // the peptide object is reused
                            peptide.set_drift_time(-1.0);
                        }

                        self.library.add_peptide(peptide.clone());
                        self.generate_transitions(&peptide.id, mz, charge, &iso_dist)?;
                        rt_map_insert(peptide_ids, rt_key, index);
                    }
                }
            } else {
                peptide.sequence = seq.to_string();
                // the protein accessions of the first hit
                let Some(first_charge_map) = self.peptide_map.get(key) else {
                    continue;
                };
                // The source assumes that at least one hit is present; the
                // map was filled from identifications with a hit.
                let first_index = first_charge_map
                    .values()
                    .flat_map(|rt_map| rt_map.first())
                    .next()
                    .map(|entry| entry.1);
                let mut accessions: BTreeSet<String> = BTreeSet::new();
                if let Some(index) = first_index {
                    if let Some(hit) = peptides[index].hits.first() {
                        accessions = hit
                            .protein_accessions()
                            .into_iter()
                            .map(str::to_owned)
                            .collect();
                    }
                }
                protein_accessions.extend(accessions.iter().cloned());
                // a missing accession would crash the OpenSWATH algorithms
                if accessions.is_empty() {
                    accessions.insert("not_available".to_owned());
                }
                peptide.protein_refs = accessions.into_iter().collect();

                // regions in RT in which the peptide elutes (ideally one)
                let mut charge_map = match self.peptide_map.get(key) {
                    Some(map) => map.clone(),
                    None => continue,
                };
                let regions = self.rt_regions(&mut charge_map, clear_ids);
                if clear_ids {
                    if let Some(stored) = self.peptide_map.get_mut(key) {
                        for rt_map in stored.values_mut() {
                            rt_map.clear();
                        }
                    }
                }
                // the IM values are stored in the identifications, per
                // peptide, charge and region

                // isotope distribution of the peptide
                let formula = seq.formula_for_source(PeptideFragmentType::Full, 0)?;
                let mut iso_dist = generator.run(&formula)?;
                trim_and_renormalize(&mut iso_dist, self.settings.isotope_pmin)?;

                // the different charge states
                let charges: Vec<i32> = charge_map.keys().copied().collect();
                for charge in charges {
                    if charge == 0 {
                        return Err(bad(format!(
                            "the identification of '{}' has no charge state",
                            peptide.sequence
                        )));
                    }
                    let mz = source_mz(seq, charge)?;
                    peptide.set_charge_state(charge);
                    let peptide_id = format!("{}/{}", peptide.sequence, charge);

                    // one feature per peptide and charge: group the regions
                    peptide.set_peptide_group_label(peptide_id.clone());
                    peptide.rts.clear();
                    let mut counter = 0usize;
                    // identifications accumulated over the regions
                    let peptide_ids = ref_rt_map.entry(peptide_id.clone()).or_default();
                    for region in &regions {
                        if let Some(region_ids) = region.ids.get(&charge) {
                            peptide.id = peptide_id.clone();
                            if regions.len() > 1 {
                                counter += 1;
                                peptide.id.push_str(&format!(":{counter}"));
                            }
                            // beginning and end of the RT region
                            peptide.rts.clear();
                            add_peptide_rt(&mut peptide, region.start);
                            add_peptide_rt(&mut peptide, region.end);

                            // IM statistics of the peptide and charge in the
                            // region; the median is more robust than the mean
                            let stats = self.region_im_stats(region, peptides)?;
                            if stats.median >= 0.0 {
                                peptide.set_drift_time(stats.median);
                                self.im_stats.insert(peptide.id.clone(), stats);
                            } else {
                                // the peptide object is reused
                                peptide.set_drift_time(-1.0);
                            }

                            self.library.add_peptide(peptide.clone());
                            self.generate_transitions(&peptide.id, mz, charge, &iso_dist)?;
                            for &(rt, index) in region_ids {
                                rt_map_insert(peptide_ids, rt, index);
                            }
                        }
                    }
                }
            }
        }
        // add the proteins to the library
        for accession in protein_accessions {
            let protein = Protein {
                id: accession,
                ..Protein::default()
            };
            self.library.add_protein(protein);
        }
        Ok(())
    }

    /// The transitions (isotopic traces) of a peptide ion (source
    /// `generateTransitions_`).
    fn generate_transitions(
        &mut self,
        peptide_id: &str,
        mz: f64,
        charge: i32,
        iso_dist: &IsotopeDistribution,
    ) -> Result<()> {
        for (counter, iso) in iso_dist.peaks().iter().enumerate() {
            let mut transition = ReactionMonitoringTransition::new();
            let annotation = format!("i{}", counter + 1);
            let transition_name = format!("{peptide_id}_{annotation}");
            transition.set_native_id(transition_name.clone());
            transition.set_precursor_mz(mz);
            // `float(counter)` in the source
            transition.set_product_mz(
                mz + C13C12_MASSDIFF_U * f64::from(counter as f32) / f64::from(charge),
            );
            transition.set_library_intensity(iso.probability);
            transition
                .cv_terms
                .metadata
                .insert("annotation".to_owned(), MetaValue::from(annotation));
            transition.set_peptide_ref(peptide_id);
            // The source notes that the transition charge is not set, while
            // many DIA scores default to charge 1 then.
            self.library.add_transition(transition);
            self.isotope_probs.insert(transition_name, iso.probability);
        }
        Ok(())
    }

    // ------------------------------------------------------------------
    // annotation, filtering
    // ------------------------------------------------------------------

    /// Assign the identifications of an assay to its best feature (source
    /// `annotateFeaturesFinalizeAssay_`): the feature with the most matching
    /// identifications, ties broken by intensity, becomes the "positive"
    /// candidate; the other identifications of the region become unassigned.
    fn finalize_assay(
        features: &mut FeatureMap,
        feat_ids: &mut BTreeMap<usize, Vec<usize>>,
        rt_ids: &mut RtMap,
        peptides: &[PeptideIdentification],
    ) {
        let mut assigned: BTreeSet<usize> = BTreeSet::new();
        if !feat_ids.is_empty() {
            // find the "best" feature, the one with the most IDs
            let mut best_index = 0usize;
            let mut best_count = 0usize;
            for (&current_index, ids) in feat_ids.iter() {
                let current_count = ids.len();
                // ties are broken by intensity
                if current_count > best_count
                    || (current_count == best_count
                        && features.features[current_index].intensity
                            > features.features[best_index].intensity)
                {
                    best_count = current_count;
                    best_index = current_index;
                }
            }
            // assign the IDs
            if best_count > 0 {
                // the one feature with the most matching IDs is defined as correct
                let feature = &mut features.features[best_index];
                feature
                    .metadata
                    .insert("feature_class".to_owned(), MetaValue::from("positive"));
                let ids = &feat_ids[&best_index];
                feature.peptide_identifications =
                    ids.iter().map(|&i| peptides[i].clone()).collect();
                assigned.extend(ids.iter().copied());
            }
        }
        // store the unassigned IDs of the current RT region
        for &(_, index) in rt_ids.iter() {
            if !assigned.contains(&index) {
                features
                    .unassigned_peptide_identifications
                    .push(peptides[index].clone());
            }
        }
        // clean up
        feat_ids.clear();
        rt_ids.clear();
    }

    /// Annotate the features with m/z, charge, isotope probabilities and the
    /// identifications (source `annotateFeatures_`).
    fn annotate_features(
        &mut self,
        features: &mut FeatureMap,
        ref_rt_map: &mut PeptideRefRtMap,
        peptides: &[PeptideIdentification],
    ) -> Result<()> {
        let mut previous_ref = String::new();
        let mut peptide_ref = String::new();
        let mut feat_ids: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for i in 0..features.features.len() {
            {
                let feat = &mut features.features[i];
                let precursor_mz = feat
                    .metadata
                    .get("PrecursorMZ")
                    .ok_or_else(|| {
                        Error::MissingInformation(
                            "a feature has no 'PrecursorMZ' meta value".into(),
                        )
                    })?
                    .as_f64()?;
                feat.mz = precursor_mz;
                feat.charge = feat
                    .peptide_identifications
                    .first()
                    .and_then(|id| id.hits.first())
                    .map(|hit| hit.charge)
                    .ok_or_else(|| {
                        Error::MissingInformation(
                            "a picked feature has no peptide identification".into(),
                        )
                    })?;
                ensure_convex_hulls(feat, self.settings.mz_window, self.settings.mz_window_ppm)?;
                // Remove the "fake" IDs generated by OpenSWATH: they would be
                // removed with a warning when writing the output, as no
                // protein identification with their identifier exists.
                feat.peptide_identifications.clear();
                // annotate the subordinates with the theoretical isotope
                // intensities
                for sub in &mut feat.subordinates {
                    let native_id = meta_text(&sub.metadata, "native_id");
                    let probability = *self.isotope_probs.entry(native_id).or_insert(0.0);
                    sub.metadata.insert(
                        "isotope_probability".to_owned(),
                        MetaValue::try_from(probability)?,
                    );
                }
                peptide_ref = meta_text(&feat.metadata, "PeptideRef");

                // IM statistics (if available) as QC values; the full peptide
                // reference, with the region number, is the key
                if let Some(stats) = self.im_stats.get(&peptide_ref) {
                    feat.metadata
                        .insert("IM_median".to_owned(), MetaValue::try_from(stats.median)?);
                    feat.metadata
                        .insert("IM_min".to_owned(), MetaValue::try_from(stats.min)?);
                    feat.metadata
                        .insert("IM_max".to_owned(), MetaValue::try_from(stats.max)?);
                }
            }
            let full_ref = peptide_ref.clone();

            // remove the region number, if present
            let slash = peptide_ref.rfind('/');
            let from = slash.map_or(1, |p| p + 2);
            if let Some(colon) = peptide_ref
                .get(from..)
                .and_then(|tail| tail.find(':'))
                .map(|p| p + from)
            {
                peptide_ref.truncate(colon);
            }

            if peptide_ref != previous_ref {
                if !previous_ref.is_empty() {
                    let mut rt_ids =
                        std::mem::take(ref_rt_map.entry(previous_ref.clone()).or_default());
                    Self::finalize_assay(features, &mut feat_ids, &mut rt_ids, peptides);
                    ref_rt_map.insert(previous_ref.clone(), rt_ids);
                }
                previous_ref = peptide_ref.clone();
            }

            let rt_ids = ref_rt_map.entry(peptide_ref.clone()).or_default();
            if rt_ids.is_empty() {
                return Err(bad(format!(
                    "No peptide IDs for peptide '{peptide_ref}' stored as '{full_ref}'."
                )));
            }

            // map the IDs to the features (by RT)
            let (rt_min, rt_max) = feature_bounds(&features.features[i])?;
            let (rt_min, rt_max) = self.calculate_rt_bounds(rt_min, rt_max);
            let lower = rt_ids.partition_point(|e| e.0.total_cmp(&rt_min) == Ordering::Less);
            let upper = rt_ids.partition_point(|e| e.0.total_cmp(&rt_max) != Ordering::Greater);
            let matched = &rt_ids[lower..upper.max(lower)];
            let id_count = matched.len();
            if !matched.is_empty() {
                feat_ids
                    .entry(i)
                    .or_default()
                    .extend(matched.iter().map(|entry| entry.1));
            }
            // "total" only includes the IDs of this RT region
            let total = rt_ids.len();
            let feat = &mut features.features[i];
            feat.metadata.insert(
                "n_total_ids".to_owned(),
                MetaValue::try_from(u64::try_from(total).unwrap_or(u64::MAX))?,
            );
            feat.metadata.insert(
                "n_matching_ids".to_owned(),
                MetaValue::try_from(u64::try_from(id_count).unwrap_or(u64::MAX))?,
            );
            let class = if id_count > 0 {
                "ambiguous" // matching IDs: the feature may be correct
            } else {
                "negative" // no matching IDs: the feature is wrong
            };
            feat.metadata
                .insert("feature_class".to_owned(), MetaValue::from(class));
        }
        // the features of the last assay
        {
            let mut rt_ids = std::mem::take(ref_rt_map.entry(peptide_ref.clone()).or_default());
            Self::finalize_assay(features, &mut feat_ids, &mut rt_ids, peptides);
            ref_rt_map.insert(peptide_ref, rt_ids);
        }
        // the unassigned identifications of assays without feature candidates
        for rt_ids in ref_rt_map.values() {
            for &(_, index) in rt_ids {
                features
                    .unassigned_peptide_identifications
                    .push(peptides[index].clone());
            }
        }
        Ok(())
    }

    /// The steps after the annotation (source `postProcess_`): unique ID, the
    /// optional candidates file, filtering and the elution model fit.
    fn post_process(&mut self, features: &mut FeatureMap) -> Result<()> {
        // make sure the map has a unique ID
        self.ensure_map_unique_id(features)?;

        // store the feature candidates before the filtering
        if self.capture_candidates {
            self.captured_candidates = Some(features.clone());
        }
        #[cfg(feature = "featurexml")]
        if !self.settings.candidates_out.is_empty() {
            // The source's writer omits (with a warning) every peptide
            // identification whose run is not among the protein identifications
            // of the map. The proteins are set after this point of a run, so
            // the candidates are written without identifications; the native
            // writer refuses such a map, so the omission is made here.
            let mut written = features.clone();
            omit_identifications_without_run(&mut written);
            crate::format::file_handler::FileHandler::store_feature_map(
                &self.settings.candidates_out,
                &written,
                Some(crate::format::FileType::FeatureXml),
            )?;
        }

        self.filter_features(features)?;
        self.report.features_after_filtering = features.features.len();
        self.report.note(format!(
            "{} features left after filtering.",
            features.features.len()
        ));

        if features.features.is_empty() {
            return Ok(()); // the elution model fit throws on empty features
        }

        // MRMFeatureFinderScoring scores an elution model too (EMG); the fit
        // here is the second one, as in the source
        if self.settings.elution_model != "none" {
            let fitter = common::elution_model_fitter(&self.param, &self.settings.elution_model)?;
            let fit = fitter.fit_elution_models(features)?;
            self.report.elution_model = Some(fit);
        } else if !self.settings.candidates_out.is_empty() {
            // the hulls are not needed, remove them
            for feature in &mut features.features {
                for sub in &mut feature.subordinates {
                    sub.convex_hulls.clear();
                }
            }
        }
        Ok(())
    }

    /// Re-process candidates of an earlier run (source `runOnCandidates`):
    /// filter them and fit the elution models, as the run does after
    /// annotation.
    ///
    /// # Errors
    ///
    /// The errors of filtering and of the elution model fit.
    pub fn run_on_candidates(&mut self, features: &mut FeatureMap) -> Result<()> {
        // extract the ID information for the statistics
        self.peptide_map.clear();
        for id in &features.unassigned_peptide_identifications {
            if let Some(hit) = id.hits.first() {
                self.peptide_map
                    .entry(SequenceKey(hit.sequence.clone()))
                    .or_default();
            }
        }
        for feature in &features.features {
            let Some(id) = feature.peptide_identifications.first() else {
                continue;
            };
            if let Some(hit) = id.hits.first() {
                self.peptide_map
                    .entry(SequenceKey(hit.sequence.clone()))
                    .or_default();
            }
        }
        self.n_peptides = self.peptide_map.len();

        let mut staged = features.clone();
        staged
            .unassigned_peptide_identifications
            .sort_by(peptide_compare);
        staged.features.sort_by(feature_compare);
        self.post_process(&mut staged)?;
        let statistics = self.statistics(&staged)?;
        self.report.peptides_identified = statistics.identified;
        self.report.peptides_quantified = statistics.quantified;
        self.report.peptides_missing = statistics.missing;
        self.report.note(statistics.text);
        *features = staged;
        Ok(())
    }

    /// Remove features without identification, and seed-derived features whose
    /// apex is too far from the seed (source `filterFeatures_`).
    fn filter_features(&mut self, features: &mut FeatureMap) -> Result<()> {
        if features.features.is_empty() {
            return Ok(());
        }
        // features without ID (or only a pseudo ID of a seed) are removed
        features
            .features
            .retain(|f| !f.peptide_identifications.is_empty());

        // Seed-derived features whose detected apex RT deviates too far from
        // the seed apex are picked from the wrong region of the chromatogram.
        // Disabled when the tolerance is not positive.
        let tolerance = self.settings.seed_apex_rt_tolerance;
        if tolerance <= 0.0 {
            return Ok(());
        }
        let before = features.features.len();
        features.features.retain(|feature| {
            let feature_rt = feature.rt;
            for id in &feature.peptide_identifications {
                // Offset peptides reuse the pseudo marker of seeds but are
                // intentional mass-shifted decoys; filtering them by apex RT
                // would bias the false-transfer estimate.
                if id.metadata.contains_key(OFFSET_PEPTIDE) {
                    continue;
                }
                for hit in &id.hits {
                    if is_seed_pseudo_hit(hit) {
                        // the seed apex RT is stored on the identification
                        if (id.rt.unwrap_or(f64::NAN) - feature_rt).abs() > tolerance {
                            return false; // too far from the seed apex: remove
                        }
                        break; // no need to check further hits
                    }
                }
            }
            true
        });
        let removed = before - features.features.len();
        self.report.seed_apex_features_removed = removed;
        if removed > 0 {
            self.report.note(format!(
                "Removed {removed} seed-derived feature(s) with detected apex RT deviating more than {tolerance} s from the seed apex."
            ));
        }
        Ok(())
    }
}

/// Drop the peptide identifications whose identifier names no protein
/// identification run of the map, from the features (subordinates included) and
/// the unassigned list, as `FeatureXMLHandler::writePeptideIdentification_`
/// omits them.
#[cfg(feature = "featurexml")]
fn omit_identifications_without_run(map: &mut FeatureMap) {
    fn strip(feature: &mut Feature, runs: &BTreeSet<&str>) {
        feature
            .peptide_identifications
            .retain(|id| runs.contains(id.identifier.as_str()));
        for sub in &mut feature.subordinates {
            strip(sub, runs);
        }
    }
    let runs: BTreeSet<String> = map
        .protein_identifications
        .iter()
        .map(|run| run.identifier.clone())
        .collect();
    let runs: BTreeSet<&str> = runs.iter().map(String::as_str).collect();
    for feature in &mut map.features {
        strip(feature, &runs);
    }
    map.unassigned_peptide_identifications
        .retain(|id| runs.contains(id.identifier.as_str()));
}

/// Remove the pseudo hits of seeds and offset peptides from the features and
/// the unassigned identifications (source `removeSeedPseudoIDs_`).
fn remove_seed_pseudo_ids(features: &mut FeatureMap) {
    for feature in &mut features.features {
        // When all assigned identifications are offset peptides the feature is
        // one too.
        let all_offsets = !feature.peptide_identifications.is_empty()
            && feature
                .peptide_identifications
                .iter()
                .all(|id| id.metadata.contains_key(OFFSET_PEPTIDE));
        if all_offsets {
            feature
                .metadata
                .insert(OFFSET_PEPTIDE.to_owned(), MetaValue::from("true"));
        }
        let ids = &mut feature.peptide_identifications;
        // the pseudo hits (PSM details)
        for id in ids.iter_mut() {
            id.hits.retain(|hit| !is_seed_pseudo_hit(hit));
        }
        // empty identifications
        ids.retain(|id| !id.hits.is_empty());
    }
    // the unassigned identifications
    let ids = &mut features.unassigned_peptide_identifications;
    for id in ids.iter_mut() {
        id.hits.retain(|hit| !is_seed_pseudo_hit(hit));
    }
    ids.retain(|id| !id.hits.is_empty());
}

/// The index ranges of the batches (source `chunk_`): `total / batch_size`
/// batches of `batch_size`, the last one extended to the end, or one batch for
/// less than a batch (also for none at all).
fn chunk_ranges(total: usize, batch_size: usize) -> Vec<(usize, usize)> {
    let batch = batch_size.max(1);
    let count = total / batch;
    let mut chunks: Vec<(usize, usize)> =
        (0..count).map(|i| (i * batch, (i + 1) * batch)).collect();
    match chunks.last_mut() {
        Some(last) => last.1 = total,
        None => chunks.push((0, total)),
    }
    chunks
}

/// Trim the isotope distribution at a minimum probability, then normalise
/// (the source does so only when the minimum probability is set).
fn trim_and_renormalize(distribution: &mut IsotopeDistribution, pmin: f64) -> Result<()> {
    if pmin > 0.0 {
        distribution.trim_left_source(pmin)?;
        distribution.trim_right(pmin)?;
        distribution.renormalize_with(ProbabilityPrecision::SourceSingle)?;
    }
    Ok(())
}

/// Add a normalised retention time to a peptide (source `addPeptideRT_`).
fn add_peptide_rt(peptide: &mut Peptide, rt: f64) {
    let mut entry = RetentionTime::new();
    entry.set_rt(rt);
    entry.retention_time_type = RTType::Normalized;
    peptide.rts.push(entry);
}

/// The unit of the first spectrum with ion mobility data (source
/// `getIMUnitFromExperiment`).
fn first_im_unit(experiment: &MSExperiment) -> Option<DriftTimeUnit> {
    experiment
        .spectra
        .iter()
        .find(|spectrum| spectrum.contains_im_data())
        .and_then(|spectrum| spectrum.im_data().ok().map(|(_, unit)| unit))
}

/// `AASequence::getMZ(charge)` of the source: `getMonoWeight(Full, charge) /
/// charge`, where the weight starts at `charge * PROTON_MASS_U`, adds the
/// terminal deltas and the residues in order and the water last. The order is
/// the one the library m/z values (and with them every score that divides by
/// them) were computed in; the cached neutral mass of the native sequence
/// adds the water first and differs in the last bit for some peptides.
///
/// # Errors
///
/// [`Error::InvalidValue`] for charge 0 (the source divides by zero) and for a
/// sequence without a known mass.
pub fn source_mz(sequence: &AASequence, charge: i32) -> Result<f64> {
    if charge == 0 {
        return Err(bad("an ion m/z needs a non-zero charge"));
    }
    if sequence.is_empty() {
        return Err(bad("an empty sequence has no ion m/z"));
    }
    // The source's weight of a residue inside a chain is the weight of the
    // free residue's formula (with its modification's formula) minus the weight
    // of water: `Residue::getMonoWeight(Internal)`.
    let water = EmpiricalFormula::parse("H2O")?.mono_mass();
    let mut weight = f64::from(charge) * PROTON_MASS_U;
    if let Some(modification) = sequence.n_terminal_modification() {
        weight += modification.diff_mono_mass()?;
    }
    if let Some(modification) = sequence.c_terminal_modification() {
        weight += modification.diff_mono_mass()?;
    }
    for index in 0..sequence.len() {
        let mut residue = sequence.subsequence(index..index + 1)?;
        if residue.n_terminal_modification().is_some() {
            residue.set_n_terminal_modification("")?;
        }
        if residue.c_terminal_modification().is_some() {
            residue.set_c_terminal_modification("")?;
        }
        let internal = match residue.residue_modification(0)? {
            // An anonymous mass: the source's "unknown modification" gives a
            // shifted residue the weight of the free residue plus the shift,
            // and leaves its formula alone. A written absolute mass is the
            // weight of the residue inside a chain.
            Some(SequenceModification::MassTag(tag)) => {
                if tag.is_delta() {
                    let letter = &sequence.as_str()[index..index + 1];
                    let base = AASequence::parse(letter)?.formula()?.mono_mass();
                    (base + tag.mass()) - water
                } else {
                    tag.residue_mono_mass().ok_or_else(|| {
                        bad("an anonymous residue mass without a known residue weight")
                    })?
                }
            }
            _ => match residue.formula() {
                Ok(formula) => formula.mono_mass() - water,
                Err(_) => residue.mono_mass()? - water,
            },
        };
        weight += internal;
    }
    Ok((weight + water) / f64::from(charge))
}

/// Replace the modification of one residue by one matching a mass difference
/// (source `AASequence::setModificationByDiffMonoMass` for an unmodified
/// residue): the best modification within 0.002 Da, else an anonymous mass tag.
fn set_modification_by_diff_mono_mass(
    sequence: &mut AASequence,
    index: usize,
    diff_mono_mass: f64,
) -> Result<()> {
    const TOLERANCE: f64 = 0.002;
    let residue = sequence
        .as_str()
        .chars()
        .nth(index)
        .ok_or_else(|| bad("residue index out of bounds"))?;
    let found = ModificationsDB::global().best_by_mass_handle(
        diff_mono_mass,
        TOLERANCE,
        Some(residue),
        Some(TermSpecificity::Anywhere),
    )?;
    match found {
        Some(modification) => {
            *sequence =
                sequence.with_resolved_modifications(&[(index, modification)], None, None)?;
            Ok(())
        }
        None => {
            // "+" or "-" and the magnitude, as `getDiffMonoMassString`
            let text = format!(
                "{}{}",
                if diff_mono_mass < 0.0 { "-" } else { "+" },
                to_str(diff_mono_mass.abs())
            );
            sequence.set_mass_tag(index, &text)
        }
    }
}
