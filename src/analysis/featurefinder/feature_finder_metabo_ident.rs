// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Targeted MS1 feature detection for metabolites from a table of targets.
//!
//! Ports `FEATUREFINDER/FeatureFinderAlgorithmMetaboIdent.h` with its `.cpp`,
//! the algorithm behind the TOPP tool `FeatureFinderMetaboIdent`. Every target
//! (name, formula or mass, charges, retention times) becomes an assay of the
//! isotope traces of its ions; the traces are extracted from the MS1 data,
//! `MRMFeatureFinderScoring` picks and scores the candidates, the candidate
//! closest to the expected RT is kept per assay, overlapping features are
//! resolved and an elution model is fitted. See
//! `docs/FEATURE_FINDER_IDENTIFICATION_SUPPORT.md` for the API mapping, the
//! preserved source conventions and the evidence.
//!
//! # Example
//!
//! ```
//! use openms::analysis::featurefinder::feature_finder_metabo_ident::FeatureFinderAlgorithmMetaboIdent;
//!
//! let algorithm = FeatureFinderAlgorithmMetaboIdent::new()?;
//! // the RT window follows the peak width: 4 * 60 s
//! assert_eq!(algorithm.rt_window(), 240.0);
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
    Compound, RTType, RTUnit, RetentionTime,
};
use crate::analysis::transformations::{DataPoint, TransformationDescription};
use crate::chemistry::isotopes::ProbabilityPrecision;
use crate::chemistry::{
    AdductInfo, CoarseIsotopePatternGenerator, CoarseMassMode, EmpiricalFormula,
    IsotopeDistribution, IsotopePeak,
};
use crate::concept::constants::user_param::{FAIMS_CV, NUM_OF_MASSTRACES};
use crate::concept::unique_id::{HasUniqueId, UniqueIdGenerator};
use crate::constants::{C13C12_MASSDIFF_U, PROTON_MASS_U};
use crate::format::file_info::text_format::to_str_f32;
use crate::identification::{PeptideIdentification, ProteinIdentification};
use crate::kernel::im_data_converter::ImDataConverter;
use crate::kernel::{Feature, FeatureMap, MSExperiment};
use crate::metadata::im_types::ImTypes;
use crate::metadata::{IonMobilityFormat, MetaValue};
use crate::param::{DefaultParamHandler, Param, ParamValue};
use crate::processing::feature_overlap_filter::FeatureOverlapFilter;
use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::{Arc, Mutex};

/// The handler name (source `DefaultParamHandler("FeatureFinderAlgorithmMetaboIdent")`).
const NAME: &str = "FeatureFinderAlgorithmMetaboIdent";

/// Largest number of library entries (targets times charge states times
/// retention times) one run builds.
///
/// Native ceiling; the source has none. It is checked before the library is
/// built, so a refusal leaves the algorithm and the output unchanged.
pub const MAX_LIBRARY_ENTRIES: usize = 1 << 24;

/// Largest number of lines [`FeatureFinderMetaboIdentReport::messages`] keeps.
pub const MAX_REPORT_MESSAGES: usize = 10_000;

/// Isotopes extracted when a minimum probability replaces the fixed number
/// (source literal `10`).
const PMIN_ISOTOPES: usize = 10;

const REMOVE_MARK: &str = "FFMetId_remove";

/// One target of the assay library (source
/// `FeatureFinderAlgorithmMetaboIdent::FeatureFinderMetaboIdentCompound`).
///
/// A plain data container. Nothing is validated here; a target that violates
/// the contract below is skipped, with a message in the
/// [report](FeatureFinderAlgorithmMetaboIdent::report), when the algorithm
/// builds the library:
///
/// * a mass above zero or a non-empty formula is required;
/// * the retention times must not be empty;
/// * a charge of zero is skipped;
/// * the retention time ranges and the ion mobilities must be empty, have one
///   entry (used for every retention time) or have one entry per retention
///   time.
#[derive(Clone, Debug, PartialEq)]
pub struct FeatureFinderMetaboIdentCompound {
    name: String,
    formula: String,
    mass: f64,
    charges: Vec<i32>,
    rts: Vec<f64>,
    rt_ranges: Vec<f64>,
    iso_distrib: Vec<f64>,
    ion_mobilities: Vec<f64>,
    adduct: String,
}

impl FeatureFinderMetaboIdentCompound {
    /// A target from its identification data (source constructor).
    ///
    /// * `name` - target name, used as a tag of the library entry;
    /// * `formula` - molecular formula, empty when only `mass` is given;
    /// * `mass` - monoisotopic mass; zero or less means "derive from the formula";
    /// * `charges` - charge states, one library entry per charge and retention time;
    /// * `rts` - expected retention times;
    /// * `rt_ranges` - RT tolerance per retention time; one value is broadcast,
    ///   none selects the algorithm default;
    /// * `iso_distrib` - isotope intensities; empty (or a first value of zero)
    ///   computes them from the formula;
    /// * `ion_mobilities` - ion mobility per retention time; one value is
    ///   broadcast, none disables the ion mobility filter;
    /// * `adduct` - adduct such as `M+H;1+`; empty infers the protonation or
    ///   deprotonation from the charge.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        name: impl Into<String>,
        formula: impl Into<String>,
        mass: f64,
        charges: Vec<i32>,
        rts: Vec<f64>,
        rt_ranges: Vec<f64>,
        iso_distrib: Vec<f64>,
        ion_mobilities: Vec<f64>,
        adduct: impl Into<String>,
    ) -> Self {
        Self {
            name: name.into(),
            formula: formula.into(),
            mass,
            charges,
            rts,
            rt_ranges,
            iso_distrib,
            ion_mobilities,
            adduct: adduct.into(),
        }
    }

    /// The target name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The molecular formula, empty if only the mass was given.
    pub fn formula(&self) -> &str {
        &self.formula
    }

    /// The monoisotopic mass; zero or less means "compute from the formula".
    pub fn mass(&self) -> f64 {
        self.mass
    }

    /// The charge states.
    pub fn charges(&self) -> &[i32] {
        &self.charges
    }

    /// The expected retention times.
    pub fn rts(&self) -> &[f64] {
        &self.rts
    }

    /// The RT tolerance per retention time (or one broadcast value, or none).
    pub fn rt_ranges(&self) -> &[f64] {
        &self.rt_ranges
    }

    /// The given isotope intensities (empty or a first value of zero: computed).
    pub fn isotope_distribution(&self) -> &[f64] {
        &self.iso_distrib
    }

    /// The ion mobilities per retention time (or one broadcast value, or none).
    pub fn ion_mobilities(&self) -> &[f64] {
        &self.ion_mobilities
    }

    /// The adduct, empty for the default of the charge's polarity.
    pub fn adduct(&self) -> &str {
        &self.adduct
    }
}

/// What a run reports in place of the source's log lines.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FeatureFinderMetaboIdentReport {
    /// Feature candidates before the selection.
    pub candidates: usize,
    /// Features left after the selection of the best candidate per assay.
    pub features_after_selection: usize,
    /// Targets skipped because they violate the contract of the table.
    pub targets_skipped: usize,
    /// The FAIMS compensation voltages processed, in order.
    pub faims_voltages: Vec<f64>,
    /// The report of the elution model fit, when one ran.
    pub elution_model: Option<ElutionModelFitReport>,
    /// The log lines, in order, at most [`MAX_REPORT_MESSAGES`].
    pub messages: Vec<String>,
    /// Lines beyond [`MAX_REPORT_MESSAGES`].
    pub dropped_messages: usize,
}

impl FeatureFinderMetaboIdentReport {
    fn note(&mut self, line: impl Into<String>) {
        if self.messages.len() < MAX_REPORT_MESSAGES {
            self.messages.push(line.into());
        } else {
            self.dropped_messages = self.dropped_messages.saturating_add(1);
        }
    }
}

/// The members that follow the parameters (source `updateMembers_`).
#[derive(Clone, Debug, PartialEq)]
struct Settings {
    rt_window: f64,
    mz_window: f64,
    mz_window_ppm: bool,
    im_window: f64,
    isotope_pmin: f64,
    n_isotopes: usize,
    peak_width: f64,
    min_peak_width: f64,
    signal_to_noise: f64,
    elution_model: String,
    candidates_out: String,
}

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

fn float_param(param: &Param, key: &str) -> Result<f64> {
    param.value(key)?.to_f64()
}

impl Settings {
    fn from_param(param: &Param, report: &mut Vec<String>) -> Result<Self> {
        let peak_width = float_param(param, "detect:peak_width")?;
        let mut rt_window = float_param(param, "extract:rt_window")?;
        if rt_window == 0.0 {
            // calculate the RT window from the peak width
            rt_window = 4.0 * peak_width;
            report.push(format!("RT window size calculated as {rt_window} seconds."));
        }
        let mz_window = float_param(param, "extract:mz_window")?;
        let isotope_pmin = float_param(param, "extract:isotope_pmin")?;
        // up to 10 isotopes if a minimum probability is given
        let n_isotopes = if isotope_pmin > 0.0 {
            PMIN_ISOTOPES
        } else {
            usize::try_from(param.value("extract:n_isotopes")?.to_i64()?)
                .map_err(|_| bad("'extract:n_isotopes' must not be negative"))?
        };
        Ok(Self {
            rt_window,
            mz_window,
            mz_window_ppm: mz_window >= 1.0,
            im_window: float_param(param, "extract:im_window")?,
            isotope_pmin,
            n_isotopes,
            peak_width,
            min_peak_width: float_param(param, "detect:min_peak_width")?,
            signal_to_noise: float_param(param, "detect:signal_to_noise")?,
            elution_model: param.value("model:type")?.as_str()?.to_owned(),
            candidates_out: param.value("candidates_out")?.as_str()?.to_owned(),
        })
    }
}

/// The text of a meta value, empty when absent (`StringUtils::toStr(DataValue)`).
fn meta_text(feature: &Feature, key: &str) -> String {
    match feature.metadata.get(key) {
        None => String::new(),
        Some(value) => match value.as_str() {
            Ok(text) => text.to_owned(),
            Err(_) => value.to_string(),
        },
    }
}

/// `FeatureCompare`: peptide reference text, then RT.
fn feature_compare(left: &Feature, right: &Feature) -> Ordering {
    meta_text(left, "PeptideRef")
        .cmp(&meta_text(right, "PeptideRef"))
        .then(left.rt.total_cmp(&right.rt))
}

/// A float feature intensity as the source's `size_t x = feature.getIntensity()`
/// reads it: truncated towards zero, negative and NaN as zero.
fn truncated_intensity(feature: &Feature) -> u64 {
    feature.intensity as u64
}

/// The target-based MS1 feature finder, the algorithm behind the TOPP tool
/// `FeatureFinderMetaboIdent` (source `FeatureFinderAlgorithmMetaboIdent`).
///
/// Create it with [`new`](Self::new) or [`from_param`](Self::from_param), hand
/// it the MS1 data with [`set_ms_data`](Self::set_ms_data), and call
/// [`run`](Self::run) with the table of targets. FAIMS data are split by
/// compensation voltage and processed per group, as for
/// [`FeatureFinderIdentificationAlgorithm`](crate::analysis::featurefinder::feature_finder_identification::FeatureFinderIdentificationAlgorithm).
///
/// The instance is meant for one run: the library, the retention times of the
/// targets and the isotope probabilities accumulate in members, as in the
/// source.
#[derive(Clone)]
pub struct FeatureFinderAlgorithmMetaboIdent {
    param: Param,
    settings: Settings,
    ms_data: MSExperiment,
    chrom_data: MSExperiment,
    feat_finder: MRMFeatureFinderScoring,
    library: TargetedExperiment,
    trafo: TransformationDescription,
    isotope_probs: BTreeMap<String, f64>,
    target_rts: BTreeMap<String, f64>,
    n_shared: usize,
    ids: Arc<Mutex<UniqueIdGenerator>>,
    report: FeatureFinderMetaboIdentReport,
}

impl std::fmt::Debug for FeatureFinderAlgorithmMetaboIdent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FeatureFinderAlgorithmMetaboIdent")
            .field("rt_window", &self.settings.rt_window)
            .field("n_shared", &self.n_shared)
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

fn int_default(param: &mut Param, key: &str, value: i64, description: &str, advanced: bool, min: i32) -> Result<()> {
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
    param.set_value(key, ParamValue::String(value.to_owned()), description, &tags)?;
    if !valid.is_empty() {
        let valid: Vec<String> = valid.iter().map(|v| (*v).to_owned()).collect();
        param.set_valid_strings(key, &valid)?;
    }
    Ok(())
}

impl FeatureFinderAlgorithmMetaboIdent {
    /// The handler name, `FeatureFinderAlgorithmMetaboIdent`.
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
            "Optional output file: Feature candidates (before filtering and model fitting).",
            &[],
            &["output file"],
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
            "extract:rt_window",
            0.0,
            "RT window size (in sec.) for chromatogram extraction. If 0, this parameter is computed as 'detect:peak_width' * 4.",
            true,
            Some(0.0),
            None,
        )?;
        float_default(
            &mut p,
            "extract:im_window",
            0.0,
            "Ion mobility window size for extraction. If set to 0, no IM filtering is performed.",
            true,
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
        p.set_section_description("extract", "Parameters for ion chromatogram extraction")?;

        float_default(
            &mut p,
            "detect:peak_width",
            60.0,
            "Expected elution peak width in seconds, for smoothing (Gauss filter). Also determines the RT extraction window, unless set explicitly via 'extract:rt_window'.",
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
        p.set_section_description(
            "detect",
            "Parameters for detecting features in extracted ion chromatograms",
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

        int_default(
            &mut p,
            "debug",
            0,
            "Debug level for feature detection.",
            true,
            0,
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
    /// violates its restriction (the source `Exception::InvalidParameter`).
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
        let ((param, settings, notes), warnings) =
            handler.set_parameters_with(parameters, |complete| {
                let mut notes = Vec::new();
                let settings = Settings::from_param(complete, &mut notes)?;
                Ok((complete.clone(), settings, notes))
            })?;
        let ids = Arc::new(Mutex::new(UniqueIdGenerator::new()));
        let mut feat_finder = MRMFeatureFinderScoring::new()?;
        feat_finder.set_shared_unique_id_generator(Arc::clone(&ids));
        let mut report = FeatureFinderMetaboIdentReport::default();
        for note in notes {
            report.note(note);
        }
        Ok((
            Self {
                param,
                settings,
                ms_data: MSExperiment::new(),
                chrom_data: MSExperiment::new(),
                feat_finder,
                library: TargetedExperiment::new(),
                trafo: TransformationDescription::default(),
                isotope_probs: BTreeMap::new(),
                target_rts: BTreeMap::new(),
                n_shared: 0,
                ids,
                report,
            },
            warnings,
        ))
    }

    /// Replace the parameters (source `setParameters`); the members follow.
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

    /// Draw the unique IDs of the features from `generator`.
    ///
    /// The source draws them from the process-wide `UniqueIdGenerator`; a
    /// caller that needs reproducible IDs seeds a generator and passes it here.
    pub fn set_unique_id_generator(&mut self, generator: UniqueIdGenerator) {
        self.ids = Arc::new(Mutex::new(generator));
        self.feat_finder
            .set_shared_unique_id_generator(Arc::clone(&self.ids));
    }

    /// The MS data (source `getMSData`); empty after a run, as the source
    /// releases the spectra once the chromatograms are extracted.
    pub fn ms_data(&self) -> &MSExperiment {
        &self.ms_data
    }

    /// Mutable access to the MS data (source non-const `getMSData`).
    pub fn ms_data_mut(&mut self) -> &mut MSExperiment {
        &mut self.ms_data
    }

    /// Set the MS data, keeping only the MS1 spectra (source `setMSData`, both
    /// overloads).
    pub fn set_ms_data(&mut self, ms_data: MSExperiment) {
        self.ms_data = ms_data;
        keep_ms1(&mut self.ms_data);
    }

    /// The extracted chromatograms (source `getChromatograms`); empty when no
    /// run was executed.
    pub fn chromatograms(&self) -> &MSExperiment {
        &self.chrom_data
    }

    /// The assay library, for example to store as TraML (source `getLibrary`);
    /// empty when no run was executed and for FAIMS data with several voltages.
    pub fn library(&self) -> &TargetedExperiment {
        &self.library
    }

    /// The deviations between the expected and the found retention times, for
    /// example to store as TrafoXML or to plot (source `getTransformations`).
    pub fn transformations(&self) -> &TransformationDescription {
        &self.trafo
    }

    /// The number of features with several target annotations (source
    /// `getNShared`).
    pub fn n_shared(&self) -> usize {
        self.n_shared
    }

    /// The RT window in seconds in effect: `extract:rt_window`, or four times
    /// the peak width when that is zero.
    pub fn rt_window(&self) -> f64 {
        self.settings.rt_window
    }

    /// The report of the last run (what the source writes to its log streams).
    pub fn report(&self) -> &FeatureFinderMetaboIdentReport {
        &self.report
    }

    /// A short description of a compound (source `prettyPrintCompound`):
    /// `name (m=<mass>, z=<charge>, rt=<expected RT>)`, the numbers as `float`.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] for a compound without the `name` or
    /// `expected_rt` meta value, or without charge state (the library builder
    /// sets all three).
    pub fn pretty_print_compound(compound: &Compound) -> Result<String> {
        let meta = &compound.cv_terms.metadata;
        let name = meta
            .get("name")
            .map(|v| v.to_string())
            .ok_or_else(|| Error::MissingInformation("the compound has no 'name'".into()))?;
        let rt = meta
            .get("expected_rt")
            .ok_or_else(|| Error::MissingInformation("the compound has no 'expected_rt'".into()))?
            .as_f64()?;
        let charge = compound
            .charge_state()
            .ok_or_else(|| Error::MissingInformation("the compound has no charge".into()))?;
        Ok(format!(
            "{name} (m={}, z={charge}, rt={})",
            to_str_f32(compound.theoretical_mass as f32),
            to_str_f32(rt as f32)
        ))
    }

    // ------------------------------------------------------------------
    // run
    // ------------------------------------------------------------------

    /// Perform the targeted feature extraction of the compounds of
    /// `targets`, and store the features in `features` (source `run`).
    ///
    /// `spectra_file` is the fall-back of the `spectra_data` annotation when
    /// the MS data carry no usable run path. Without MS1 spectra `features` is
    /// returned unchanged (apart from the run path).
    ///
    /// For FAIMS data one run per compensation voltage is made, the features
    /// are annotated with `FAIMS_CV` and features of the same analyte are
    /// merged; [`library`](Self::library) is empty then.
    ///
    /// # Errors
    ///
    /// * [`Error::InvalidValue`] when more than [`MAX_LIBRARY_ENTRIES`] entries
    ///   would be built, for a formula that cannot be parsed (the source's
    ///   `Exception::ParseError`) and for an isotope distribution that cannot be
    ///   normalised;
    /// * [`Error::Unsupported`] when `candidates_out` is set and no writer is
    ///   available (the `featurexml` feature is disabled);
    /// * every error of extraction, picking, scoring and elution model fitting.
    ///
    /// `features` is left unchanged when this returns an error.
    ///
    /// # Differences from the source
    ///
    /// The table is pre-flighted against [`MAX_LIBRARY_ENTRIES`] before the
    /// library is built.
    pub fn run(
        &mut self,
        targets: &[FeatureFinderMetaboIdentCompound],
        features: &mut FeatureMap,
        spectra_file: &str,
    ) -> Result<()> {
        // pre-flight: the number of library entries
        let mut entries = 0usize;
        for target in targets {
            entries = entries
                .saturating_add(target.charges.len().saturating_mul(target.rts.len()));
        }
        if entries > MAX_LIBRARY_ENTRIES {
            return Err(bad(format!(
                "the targets would create {entries} library entries; at most {MAX_LIBRARY_ENTRIES} are built"
            )));
        }
        if !self.settings.candidates_out.is_empty() && cfg!(not(feature = "featurexml")) {
            return Err(Error::Unsupported(
                "candidates_out requires the 'featurexml' feature".into(),
            ));
        }

        let split = ImDataConverter::split_by_faims_cv(&mut self.ms_data)?;
        if !split.has_faims() {
            // non-FAIMS data: restore the data and run directly
            let mut groups = split.groups.into_iter();
            if let Some(group) = groups.next() {
                self.ms_data = group.experiment;
            }
            let mut staged = features.clone();
            self.run_single_group(targets, &mut staged, spectra_file)?;
            *features = staged;
            return Ok(());
        }

        // FAIMS data: process each compensation voltage separately
        let mut staged = features.clone();
        staged.clear(true);
        self.chrom_data = MSExperiment::new();
        self.library.clear(true);
        self.trafo = TransformationDescription::default();
        self.n_shared = 0;
        let mut first_group = true;
        let group_count = split.groups.len();
        self.report.note(format!(
            "FAIMS data detected with {group_count} compensation voltage(s)."
        ));
        let mut points: Vec<DataPoint> = Vec::new();
        for group in split.groups {
            let group_cv = group.key.volts();
            self.report.note(format!(
                "Processing FAIMS CV group: {group_cv} V ({} spectra)",
                group.experiment.spectra.len()
            ));
            // a fresh instance per group, with the same parameters
            let mut sub = Self::from_param(&self.param)?;
            sub.ids = Arc::clone(&self.ids);
            sub.feat_finder
                .set_shared_unique_id_generator(Arc::clone(&self.ids));
            sub.set_ms_data(group.experiment);
            let mut features_cv = FeatureMap::new();
            sub.run_single_group(targets, &mut features_cv, spectra_file)?;

            let cv_value = MetaValue::try_from(group_cv)?;
            for mut feature in features_cv.features {
                feature.metadata.insert(FAIMS_CV.to_owned(), cv_value.clone());
                staged.features.push(feature);
            }
            // the proteins of the first group
            if first_group {
                staged.protein_identifications = features_cv.protein_identifications;
            }
            for mut id in features_cv.unassigned_peptide_identifications {
                id.metadata.insert(FAIMS_CV.to_owned(), cv_value.clone());
                staged.unassigned_peptide_identifications.push(id);
            }
            for mut chromatogram in sub.chrom_data.chromatograms {
                chromatogram
                    .metadata
                    .insert(FAIMS_CV.to_owned(), cv_value.clone());
                self.chrom_data.chromatograms.push(chromatogram);
            }
            points.extend_from_slice(sub.trafo.data_points());
            self.n_shared += sub.n_shared;
            self.report.faims_voltages.push(group_cv);
            first_group = false;
        }
        self.trafo.set_data_points(points)?;
        self.report.note(
            "Warning: Library output is not available for multi-FAIMS data. Each FAIMS CV group has its own assay library.",
        );
        let combined = staged.features.len();
        self.report.note(format!(
            "Combined {combined} features from all FAIMS CV groups."
        ));

        let merge = self.param.value("faims:merge_features")?.to_bool()?;
        let merged = common::merge_faims_features(&mut staged, merge)?;
        if merged > 0 {
            self.report.note(format!(
                "Merged FAIMS features: {combined} -> {} ({merged} features merged)",
                combined - merged
            ));
        }

        staged.set_primary_ms_run_path(&[spectra_file.to_owned()])?;
        self.ensure_map_unique_id(&mut staged)?;
        *features = staged;
        Ok(())
    }

    fn ensure_map_unique_id(&self, features: &mut FeatureMap) -> Result<()> {
        let mut generator = self
            .ids
            .lock()
            .map_err(|_| bad("the unique ID generator is poisoned"))?;
        features.unique_id.ensure_unique_id(&mut generator);
        Ok(())
    }

    /// Core of one dataset: not FAIMS, or a single FAIMS group (source
    /// `runSingleGroup_`).
    fn run_single_group(
        &mut self,
        targets: &[FeatureFinderMetaboIdentCompound],
        features: &mut FeatureMap,
        spectra_file: &str,
    ) -> Result<()> {
        // if a proper mzML is annotated in the data it is the reference,
        // otherwise the given file name is
        features.set_primary_ms_run_path_from_experiment(
            &[spectra_file.to_owned()],
            &self.ms_data,
        )?;

        if self.ms_data.spectra.is_empty() {
            self.report
                .note(format!("Warning: No MS1 scans in:{spectra_file}"));
            return Ok(());
        }

        // The library is built into a staging copy, so that a failure leaves
        // the members of the algorithm as they were.
        let before = (
            self.library.clone(),
            self.isotope_probs.clone(),
            self.target_rts.clone(),
        );
        let outcome = self.detect(targets, features);
        if outcome.is_err() {
            self.library = before.0;
            self.isotope_probs = before.1;
            self.target_rts = before.2;
        }
        outcome
    }

    fn detect(
        &mut self,
        targets: &[FeatureFinderMetaboIdentCompound],
        features: &mut FeatureMap,
    ) -> Result<()> {
        for target in targets {
            self.add_target_to_library(target)?;
        }

        // initialise the algorithm classes needed later
        let mut params = self.feat_finder.parameters().clone();
        let max_iteration = self.param.value("EMGScoring:max_iteration")?.clone();
        let init_mom = self.param.value("EMGScoring:init_mom")?.clone();
        // some parts of the algorithm expect convex hulls
        params.set_value(
            "write_convex_hull",
            ParamValue::String("true".into()),
            "",
            &[],
        )?;
        let model_or_candidates =
            self.settings.elution_model != "none" || !self.settings.candidates_out.is_empty();
        params.set_value(
            "Scores:use_elution_model_score",
            // TODO of the source: test if this works for requantification
            ParamValue::String(if model_or_candidates { "false" } else { "true" }.into()),
            "",
            &[],
        )?;
        // The relative minimum width becomes seconds, in place, as the source
        // does.
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

        //-------------------------------------------------------------
        // run the feature detection
        //-------------------------------------------------------------
        // does the data carry ion mobility information?
        let im_format = ImTypes::determine_im_format_for_ms_level(&self.ms_data, 1)?;
        let has_im = matches!(
            im_format,
            IonMobilityFormat::PerPeak | IonMobilityFormat::PerSpectrum
        );
        if has_im && self.settings.im_window > 0.0 {
            self.report.note(format!(
                "Ion mobility data detected. Using IM window: {}",
                self.settings.im_window
            ));
        } else if has_im {
            self.report.note(format!(
                "Warning: Ion mobility data detected, but 'extract:im_window' is {} (<= 0). The ion mobility dimension will be ignored: peaks across the whole mobility range are summed within each m/z window, which merges mobility-separated co-eluting analytes (e.g. on Bruker .d frames) and can distort extraction/quantification. Set 'extract:im_window' > 0 and provide IonMobility values for your targets to extract within the analyte's mobility band.",
                self.settings.im_window
            ));
        }

        self.report.note("Extracting chromatograms...");
        let light = convert_targeted_exp(&self.library)?;
        let shared = Arc::new(std::mem::take(&mut self.ms_data));
        let outcome = self.extract_and_pick(&shared, &light, features, has_im);
        if outcome.is_err() {
            if let Ok(data) = Arc::try_unwrap(shared) {
                self.ms_data = data;
            }
            return outcome;
        }
        // the data are not needed any more (source `ms_data_.reset()`)
        drop(shared);
        self.report.candidates = features.features.len();
        self.report.note(format!(
            "Found {} feature candidates in total.",
            features.features.len()
        ));

        // complete the feature annotation
        self.annotate_features(features)?;

        // auxiliary output
        self.ensure_map_unique_id(features)?;

        if !self.settings.candidates_out.is_empty() {
            // store the feature candidates
            features.features.sort_by(feature_compare);
            #[cfg(feature = "featurexml")]
            crate::format::file_handler::FileHandler::store_feature_map(
                &self.settings.candidates_out,
                features,
                Some(crate::format::FileType::FeatureXml),
            )?;
        }

        self.select_features_from_candidates(features)?;
        self.report.features_after_selection = features.features.len();
        self.report.note(format!(
            "{} features left after selection of best candidates.",
            features.features.len()
        ));

        // the criterion for the best of overlapping features: lower is better
        let comparator = |left: &Feature, right: &Feature| -> bool {
            let delta = |f: &Feature| -> f64 {
                f.metadata
                    .get("rt_deviation")
                    .and_then(|v| v.as_f64().ok())
                    .map_or(f64::INFINITY, f64::abs)
            };
            let (left_delta, right_delta) = (delta(left), delta(right));
            // `std::tie(left_rt_delta, right_intensity) <
            // std::tie(right_rt_delta, left_intensity)`: the intensities are
            // swapped because higher is better
            let (left_intensity, right_intensity) =
                (truncated_intensity(left), truncated_intensity(right));
            left_delta < right_delta
                || (left_delta == right_delta && right_intensity < left_intensity)
        };

        // transfer the information of an identical overlapping feature with a
        // different annotation to the representative
        let on_overlap = |representative: &mut Feature, overlap: &mut Feature| -> bool {
            let best_intensity = truncated_intensity(representative);
            let overlap_intensity = truncated_intensity(overlap);
            if overlap_intensity != best_intensity {
                return true; // early out: the features are different
            }
            // this part is nearly never reached (identical features only)
            let delta = |f: &Feature| -> f64 {
                f.metadata
                    .get("rt_deviation")
                    .and_then(|v| v.as_f64().ok())
                    .map_or(f64::NAN, f64::abs)
            };
            let best_rt_delta = delta(representative);
            let overlap_rt_delta = delta(overlap);
            if overlap_rt_delta == best_rt_delta {
                let same = overlap.mz == representative.mz && overlap.rt == representative.rt;
                if same {
                    // update the annotations
                    let mut label = meta_text(representative, "label");
                    label.push('/');
                    label.push_str(&meta_text(overlap, "label"));
                    representative
                        .metadata
                        .insert("label".to_owned(), MetaValue::from(label));
                    let mut alt_refs: Vec<String> = representative
                        .metadata
                        .get("alt_PeptideRef")
                        .and_then(|v| v.as_string_list().ok())
                        .map(<[String]>::to_vec)
                        .unwrap_or_default();
                    alt_refs.push(meta_text(overlap, "PeptideRef"));
                    representative
                        .metadata
                        .insert("alt_PeptideRef".to_owned(), MetaValue::from(alt_refs));
                }
            }
            // annotate which features were removed because of the overlap
            let reference = format!(
                "{} (RT {})",
                meta_text(overlap, "PeptideRef"),
                to_str_f32(overlap.rt as f32)
            );
            let mut overlap_refs: Vec<String> = representative
                .metadata
                .get("overlap_removed")
                .and_then(|v| v.as_string_list().ok())
                .map(<[String]>::to_vec)
                .unwrap_or_default();
            overlap_refs.push(reference);
            representative
                .metadata
                .insert("overlap_removed".to_owned(), MetaValue::from(overlap_refs));
            true
        };
        // overlap is checked at the level of the mass traces
        FeatureOverlapFilter::filter(features, comparator, on_overlap, true)?;
        // sort by reference and RT (stable)
        features.features.sort_by(feature_compare);

        if features.features.is_empty() {
            self.report.note("No features left after filtering.");
            return Ok(());
        }

        self.n_shared = self.add_target_annotations(features)?;

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

        self.extract_transformations(features)
    }

    /// The extraction and picking block of `runSingleGroup_`.
    fn extract_and_pick(
        &mut self,
        shared: &Arc<MSExperiment>,
        light: &crate::analysis::openswath::transition_experiment::LightTargetedExperiment,
        features: &mut FeatureMap,
        has_im: bool,
    ) -> Result<()> {
        let spectra = SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(Arc::clone(shared))?;
        let settings_spectrum = shared
            .spectra
            .first()
            .ok_or_else(|| bad("the MS data hold no spectrum"))?;
        let im_window = if has_im && self.settings.im_window > 0.0 {
            Some(self.settings.im_window)
        } else {
            None
        };
        let mut extracted = Vec::new();
        extract_ms1_chromatograms(
            &spectra,
            settings_spectrum,
            light,
            self.settings.mz_window,
            self.settings.mz_window_ppm,
            im_window,
            &mut extracted,
        )?;
        self.chrom_data.chromatograms.extend(extracted);
        self.report.note(format!(
            "Extracted {} chromatogram(s).",
            self.chrom_data.chromatograms.len()
        ));

        self.report.note("Detecting chromatographic peaks...");
        // `pickExperiment` of the OpenMS structures, with the spectra of the
        // run as the DIA map; the chromatograms stay with the algorithm
        let chromatograms = Arc::new(self.chrom_data.clone());
        let chromatogram_access =
            SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(chromatograms)?;
        let swath = SwathMap {
            sptr: Some(spectra),
            ..SwathMap::default()
        };
        let mut group_map = TransitionGroupMapType::new();
        self.feat_finder.pick_experiment_light(
            &chromatogram_access,
            features,
            light,
            &TransformationDescription::default(),
            &[swath],
            &mut group_map,
        )
    }

    // ------------------------------------------------------------------
    // library
    // ------------------------------------------------------------------

    /// m/z of an ion from the mass and the charge (source `calculateMZ_`).
    fn calculate_mz(mass: f64, charge: i32) -> f64 {
        (mass + f64::from(charge) * PROTON_MASS_U) / f64::from(charge).abs()
    }

    /// Add a target to the assay library (source `addTargetToLibrary_`). A
    /// target that violates the contract is skipped with a message.
    fn add_target_to_library(&mut self, t: &FeatureFinderMetaboIdentCompound) -> Result<()> {
        let name = t.name.as_str();
        if t.mass <= 0.0 && t.formula.is_empty() {
            self.report.targets_skipped += 1;
            self.report.note(format!(
                "Error: No mass or sum formula given for target '{name}' - skipping this target."
            ));
            return Ok(());
        }
        if t.rts.is_empty() {
            self.report.targets_skipped += 1;
            self.report.note(format!(
                "Error: No retention time (RT) given for target '{name}' - skipping this target."
            ));
            return Ok(());
        }
        // TODO of the source: detect entries with the same RT and m/z
        let mut target = Compound::default();
        target
            .cv_terms
            .metadata
            .insert("name".to_owned(), MetaValue::from(name));
        target.molecular_formula = t.formula.clone();
        let mut emp_formula = if t.formula.is_empty() {
            EmpiricalFormula::default()
        } else {
            EmpiricalFormula::parse(&t.formula)?
        };
        let mass_given = t.mass > 0.0;
        let mass = if mass_given {
            t.mass
        } else {
            emp_formula.mono_mass()
        };
        target.theoretical_mass = mass;
        let target_id = format!("{name}_m{}", to_str_f32(mass as f32));

        // isotope distribution of the target
        let generator = CoarseIsotopePatternGenerator::new(
            Some(self.settings.n_isotopes),
            CoarseMassMode::Approximate,
        )?
        .with_precision(ProbabilityPrecision::SourceSingle);
        let mut iso_dist: IsotopeDistribution;
        let mut n_isotopes = self.settings.n_isotopes;
        if t.iso_distrib.is_empty() || t.iso_distrib[0] == 0.0 {
            if t.formula.is_empty() {
                self.report.note(format!(
                    "Error: No sum formula given for target '{name}'; cannot calculate isotope distribution - using estimation method for peptides."
                ));
                iso_dist = generator.estimate_from_peptide_weight(mass)?;
            } else {
                iso_dist = generator.run(&emp_formula)?;
            }
        } else {
            n_isotopes = n_isotopes.min(t.iso_distrib.len());
            // the source stores the probabilities as `float`
            let peaks: Vec<IsotopePeak> = (0..n_isotopes)
                .map(|i| IsotopePeak {
                    mass: i as f64,
                    probability: f64::from(t.iso_distrib[i] as f32),
                })
                .collect();
            iso_dist = IsotopeDistribution::from_peaks(peaks)?;
        }
        if self.settings.isotope_pmin > 0.0 {
            iso_dist.trim_left_source(self.settings.isotope_pmin)?;
            iso_dist.trim_right(self.settings.isotope_pmin)?;
        }
        iso_dist.renormalize_with(ProbabilityPrecision::SourceSingle)?;

        // IM values: none (no filtering), one (for all) or one per RT; -1
        // matches the default drift time and the extractor's `>= 0` check
        let mut ims = t.ion_mobilities.clone();
        if ims.is_empty() {
            ims = vec![-1.0; t.rts.len()];
        } else if ims.len() == 1 {
            ims = vec![ims[0]; t.rts.len()];
        } else if ims.len() != t.rts.len() {
            self.report.targets_skipped += 1;
            self.report.note(format!(
                "Error: Number of IonMobility values ({}) does not match number of RT values ({}) for target '{name}' - skipping this target.",
                ims.len(),
                t.rts.len()
            ));
            return Ok(());
        }

        // RT ranges: none (default), one (for all) or one per RT
        let mut rt_ranges = t.rt_ranges.clone();
        if rt_ranges.is_empty() {
            rt_ranges = vec![0.0; t.rts.len()];
        } else if rt_ranges.len() == 1 {
            rt_ranges = vec![rt_ranges[0]; t.rts.len()];
        } else if rt_ranges.len() != t.rts.len() {
            self.report.targets_skipped += 1;
            self.report.note(format!(
                "Error: Number of RetentionTimeRange values ({}) does not match number of RT values ({}) for target '{name}' - skipping this target.",
                rt_ranges.len(),
                t.rts.len()
            ));
            return Ok(());
        }

        // the different charge states
        for &charge in &t.charges {
            if charge == 0 {
                self.report.note(format!(
                    "Error: Invalid charge 0 for target '{name}' - skipping this charge."
                ));
                continue;
            }
            target.set_charge_state(charge);

            // the adduct for the m/z: protonation for a positive and
            // deprotonation for a negative charge by default
            let adduct_str = if t.adduct.is_empty() {
                if charge > 0 {
                    if charge == 1 {
                        "M+H;1+".to_owned()
                    } else {
                        format!("M+{charge}H;{charge}+")
                    }
                } else if charge == -1 {
                    "M-H;1-".to_owned()
                } else {
                    let magnitude = charge.unsigned_abs();
                    format!("M-{magnitude}H;{magnitude}-")
                }
            } else {
                t.adduct.clone()
            };
            let parsed = AdductInfo::parse(&adduct_str).and_then(|adduct_info| {
                // `emp_formula` has charge 0 here, so its monoisotopic mass is
                // the neutral one; `mass` holds it unless a mass was given
                let neutral = if mass_given {
                    t.mass
                } else {
                    emp_formula.mono_mass()
                };
                adduct_info.mz(neutral)
            });
            let mz = match parsed {
                Ok(mz) => {
                    target
                        .cv_terms
                        .metadata
                        .insert("adduct".to_owned(), MetaValue::from(adduct_str));
                    mz
                }
                Err(error) => {
                    self.report.note(format!(
                        "Error: Could not parse adduct string '{adduct_str}' for target '{name}': {error} - falling back to default calculation."
                    ));
                    if mass_given {
                        Self::calculate_mz(t.mass, charge)
                    } else {
                        // the source sets the charge of the formula for good
                        emp_formula = emp_formula.with_charge(charge);
                        (emp_formula.mono_mass() / f64::from(charge)).abs()
                    }
                }
            };

            for (i, &rt) in t.rts.iter().enumerate() {
                target.id = format!("{target_id}_z{charge}_rt{}", to_str_f32(rt as f32));
                target
                    .cv_terms
                    .metadata
                    .insert("expected_rt".to_owned(), MetaValue::try_from(rt)?);
                self.target_rts.insert(target.id.clone(), rt);

                // the ion mobility, for the extraction
                if ims[i] >= 0.0 {
                    target.set_drift_time(ims[i]);
                    target
                        .cv_terms
                        .metadata
                        .insert("expected_im".to_owned(), MetaValue::try_from(ims[i])?);
                } else {
                    // the target is reused across the iterations
                    target.set_drift_time(-1.0);
                }

                let mut rt_tol = rt_ranges[i] / 2.0;
                if rt_tol == 0.0 {
                    rt_tol = self.settings.rt_window / 2.0;
                }
                // beginning and end of the RT region
                target.rts.clear();
                add_target_rt(&mut target, rt - rt_tol);
                add_target_rt(&mut target, rt + rt_tol);
                self.library.add_compound(target.clone());
                self.generate_transitions(&target.id, mz, charge, &iso_dist);
            }
        }
        Ok(())
    }

    /// The transitions (isotopic traces) of a target ion (source
    /// `generateTransitions_`).
    fn generate_transitions(
        &mut self,
        target_id: &str,
        mz: f64,
        charge: i32,
        iso_dist: &IsotopeDistribution,
    ) {
        for (counter, iso) in iso_dist.peaks().iter().enumerate() {
            let mut transition = ReactionMonitoringTransition::new();
            let annotation = format!("i{counter}");
            let transition_name = format!("{target_id}_{annotation}");
            transition.set_native_id(transition_name.clone());
            transition.set_precursor_mz(mz);
            // TODO of the source: accurate masses from the distribution
            transition.set_product_mz(
                mz + (C13C12_MASSDIFF_U * f64::from(counter as f32) / f64::from(charge)).abs(),
            );
            transition.set_library_intensity(iso.probability);
            transition.set_compound_ref(target_id);
            self.library.add_transition(transition);
            self.isotope_probs.insert(transition_name, iso.probability);
        }
    }

    // ------------------------------------------------------------------
    // annotation and selection
    // ------------------------------------------------------------------

    /// Add the annotations of the target to the features (source
    /// `annotateFeatures_`).
    fn annotate_features(&mut self, features: &mut FeatureMap) -> Result<()> {
        for feat in &mut features.features {
            let precursor_mz = feat
                .metadata
                .get("PrecursorMZ")
                .ok_or_else(|| {
                    Error::MissingInformation("a feature has no 'PrecursorMZ' meta value".into())
                })?
                .as_f64()?;
            feat.mz = precursor_mz;
            let reference = meta_text(feat, "PeptideRef");
            let compound = self.library.compound_by_ref(&reference)?;
            feat.charge = compound.charge_state().unwrap_or(0);
            ensure_convex_hulls(feat, self.settings.mz_window, self.settings.mz_window_ppm)?;
            feat.peptide_identifications.clear();
            let meta = &compound.cv_terms.metadata;
            feat.metadata.insert(
                "label".to_owned(),
                meta.get("name").cloned().unwrap_or_default(),
            );
            feat.metadata.insert(
                "sum_formula".to_owned(),
                MetaValue::from(compound.molecular_formula.clone()),
            );
            if let Some(adduct) = meta.get("adduct") {
                feat.metadata.insert("adduct".to_owned(), adduct.clone());
            }
            feat.metadata.insert(
                "expected_rt".to_owned(),
                meta.get("expected_rt").cloned().unwrap_or_default(),
            );
            // ion mobility annotations, if available
            if let Some(im) = meta.get("expected_im") {
                feat.metadata.insert("expected_im".to_owned(), im.clone());
            }
            // annotate the subordinates with the theoretical isotope intensities
            for sub in &mut feat.subordinates {
                let native_id = meta_text(sub, "native_id");
                let probability = *self.isotope_probs.entry(native_id).or_insert(0.0);
                sub.metadata.insert(
                    "isotope_probability".to_owned(),
                    MetaValue::try_from(probability)?,
                );
                // the value "MS2" is misleading
                sub.metadata.remove("FeatureLevel");
            }
            // required for SIRIUS
            feat.metadata.insert(
                NUM_OF_MASSTRACES.to_owned(),
                MetaValue::try_from(u64::try_from(feat.subordinates.len()).unwrap_or(u64::MAX))?,
            );
        }
        features.protein_identifications.clear();
        Ok(())
    }

    /// Keep the best candidate of each assay (source
    /// `selectFeaturesFromCandidates_`): the candidate whose RT range
    /// contains the expected RT, or else the one closest to it, ties broken
    /// by intensity. The candidates of one assay are consecutive.
    fn select_features_from_candidates(&mut self, features: &mut FeatureMap) -> Result<()> {
        let mut previous_ref = String::new();
        let mut best_rt_dist = f64::INFINITY;
        let mut best = 0usize;
        for it in 0..features.features.len() {
            // the features of one assay (same "PeptideRef") are consecutive
            let reference = meta_text(&features.features[it], "PeptideRef");
            if reference != previous_ref {
                // a new assay
                previous_ref = reference.clone();
                best_rt_dist = self.settings.rt_window;
                best = it;
            }
            let target_rt = *self.target_rts.entry(reference.clone()).or_insert(0.0);
            let (rt_min, rt_max) = feature_bounds(&features.features[it])?;
            let mut rt_dist = f64::INFINITY;
            if rt_min <= target_rt && rt_max >= target_rt {
                if best_rt_dist <= 0.0 {
                    self.report.note(format!(
                        "Warning: overlapping feature candidates for assay '{reference}'"
                    ));
                }
                rt_dist = 0.0;
            } else if best_rt_dist > 0.0 {
                rt_dist = if rt_min > target_rt {
                    rt_min - target_rt
                } else {
                    target_rt - rt_max
                };
            }
            if rt_dist < best_rt_dist
                || (rt_dist == best_rt_dist
                    && features.features[it].intensity > features.features[best].intensity)
            {
                // a new best candidate for the assay
                best_rt_dist = rt_dist;
                // mark the no-longer-best candidate for removal
                if best != it {
                    features.features[best]
                        .metadata
                        .insert(REMOVE_MARK.to_owned(), MetaValue::from(""));
                }
                best = it;
                let deviation = target_rt - features.features[best].rt;
                features.features[best]
                    .metadata
                    .insert("rt_deviation".to_owned(), MetaValue::try_from(deviation)?);
            } else {
                // worse than a previous candidate: mark for removal
                features.features[it]
                    .metadata
                    .insert(REMOVE_MARK.to_owned(), MetaValue::from(""));
            }
        }
        features
            .features
            .retain(|feature| !feature.metadata.contains_key(REMOVE_MARK));
        Ok(())
    }

    /// Add "peptide" identifications with the information of the targets
    /// without a feature (source `addTargetAnnotations_`). Returns the number
    /// of features that carry several targets.
    fn add_target_annotations(&mut self, features: &mut FeatureMap) -> Result<usize> {
        let mut n_shared = 0usize;
        let mut found_refs: BTreeSet<String> = BTreeSet::new();
        for feature in &features.features {
            found_refs.insert(meta_text(feature, "PeptideRef"));
            if let Some(value) = feature.metadata.get("alt_PeptideRef") {
                n_shared += 1;
                if let Ok(list) = value.as_string_list() {
                    found_refs.extend(list.iter().cloned());
                }
            }
        }
        // the targets without features
        let n_missing = self.library.compounds().len().saturating_sub(found_refs.len());
        for compound in self.library.compounds() {
            if !found_refs.contains(&compound.id) {
                let meta = &compound.cv_terms.metadata;
                let mut metadata = crate::metadata::MetaInfo::new();
                metadata.insert(
                    "label".to_owned(),
                    meta.get("name").cloned().unwrap_or_default(),
                );
                metadata.insert("PeptideRef".to_owned(), MetaValue::from(compound.id.clone()));
                if let Some(im) = meta.get("expected_im") {
                    // ion mobility annotation, if available
                    metadata.insert("expected_im".to_owned(), im.clone());
                }
                let expected_rt = meta
                    .get("expected_rt")
                    .map(MetaValue::as_f64)
                    .transpose()?;
                let peptide = PeptideIdentification {
                    identifier: "id".to_owned(),
                    rt: expected_rt,
                    mz: Some(Self::calculate_mz(
                        compound.theoretical_mass,
                        compound.charge_state().unwrap_or(0),
                    )),
                    metadata,
                    ..PeptideIdentification::default()
                };
                features.unassigned_peptide_identifications.push(peptide);
            }
            if features.unassigned_peptide_identifications.len() >= n_missing {
                break; // found all
            }
        }
        if n_missing > 0 {
            features
                .protein_identifications
                .resize(1, ProteinIdentification::default());
            features.protein_identifications[0].identifier = "id".to_owned();
        }
        Ok(n_shared)
    }

    /// The expected against the found retention times (source
    /// `extractTransformations_`).
    fn extract_transformations(&mut self, features: &FeatureMap) -> Result<()> {
        let mut points = Vec::with_capacity(features.features.len());
        for feature in &features.features {
            let expected = feature
                .metadata
                .get("expected_rt")
                .ok_or_else(|| {
                    Error::MissingInformation("a feature has no 'expected_rt' meta value".into())
                })?
                .as_f64()?;
            points.push(DataPoint {
                x: expected,
                y: feature.rt,
                note: meta_text(feature, "PeptideRef"),
            });
        }
        self.trafo.set_data_points(points)
    }
}

/// Add a local retention time in seconds to a target (source `addTargetRT_`).
fn add_target_rt(target: &mut Compound, rt: f64) {
    let mut entry = RetentionTime::new();
    entry.retention_time_unit = RTUnit::Second;
    entry.retention_time_type = RTType::Local;
    entry.set_rt(rt);
    target.rts.push(entry);
}
