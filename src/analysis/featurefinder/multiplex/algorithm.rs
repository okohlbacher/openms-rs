// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The multiplex feature finder: filtering, clustering and quantification of
//! labelled peptide multiplets.
//!
//! Native equivalent of `FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h`. See
//! `docs/FEATUREFINDER_MULTIPLEX_SUPPORT.md`.

use std::collections::BTreeMap;

use super::clustering::MultiplexClustering;
use super::delta_masses::MultiplexDeltaMasses;
use super::delta_masses_generator::{
    LabelMassShifts, MultiplexDeltaMassesGenerator, label_master_list,
};
use super::filtered_ms_experiment::MultiplexFilteredMSExperiment;
use super::filtered_peak::SatelliteMultimap;
use super::filtering_centroided::MultiplexFilteringCentroided;
use super::filtering_profile::MultiplexFilteringProfile;
use super::isotopic_peak_pattern::MultiplexIsotopicPeakPattern;
use super::satellite::{MultiplexSatelliteCentroided, MultiplexSatelliteProfile};
use crate::analysis::peak_integrator::{IntegrationMethod, PeakIntegrator};
use crate::concept::progress_logger::ProgressReporter;
use crate::concept::unique_id::UniqueIdGenerator;
use crate::kernel::spectrum_helper::{IntensityAveragingMethod, make_peak_position_unique};
use crate::kernel::{
    BaseFeature, ChromatogramPeak, ConsensusFeature, ConsensusMap, ConvexHull2D, Feature,
    FeatureHandle, FeatureMap, MSChromatogram, MSExperiment, Point2D, SpectrumType,
};
use crate::math::ml::clustering::grid_based_cluster::GridBasedCluster;
use crate::math::source_sort::source_sort_by;
use crate::metadata::MetaValue;
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::processing::peak_picking::{PeakBoundary, PeakPickerHiRes, PickingCompatibility};
use crate::processing::spline::spline_package::SplinePackage;
use crate::{Error, Result};

/// Name of the source's parameter handler.
pub const FEATURE_FINDER_MULTIPLEX_NAME: &str = "FeatureFinderMultiplexAlgorithm";

/// Ceiling on the steps of one ratio scan in the intensity correction,
/// checked before the scan. The source has none.
pub const MAX_RATIO_SCAN_STEPS: usize = 100_000_000;

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

/// The parameters of [`FeatureFinderMultiplexAlgorithm`], section `algorithm:`
/// and the label masses of section `labels:`.
///
/// Each field is the parameter of the same name; see
/// [`FeatureFinderMultiplexAlgorithm::defaults`] for the descriptions and
/// restrictions, which are the source's.
#[derive(Clone, Debug, PartialEq)]
pub struct FeatureFinderMultiplexParameters {
    /// `algorithm:labels`, default `[][Lys8,Arg10]`.
    pub labels: String,
    /// `algorithm:charge`, `min:max`, default `1:4`; parsed when the algorithm
    /// runs, as in the source.
    pub charge: String,
    /// `algorithm:isotopes_per_peptide`, `min:max`, default `3:6`; parsed when
    /// the algorithm runs.
    pub isotopes_per_peptide: String,
    /// `algorithm:rt_typical` in seconds, default 40.
    pub rt_typical: f64,
    /// `algorithm:rt_band` in seconds, default 0.
    pub rt_band: f64,
    /// `algorithm:rt_min` in seconds, default 2.
    pub rt_min: f64,
    /// `algorithm:mz_tolerance`, default 6.
    pub mz_tolerance: f64,
    /// `algorithm:mz_unit`, `ppm` (default) or `Da`.
    pub mz_unit: String,
    /// `algorithm:intensity_cutoff`, default 1000.
    pub intensity_cutoff: f64,
    /// `algorithm:peptide_similarity`, default 0.5.
    pub peptide_similarity: f64,
    /// `algorithm:averagine_similarity`, default 0.4.
    pub averagine_similarity: f64,
    /// `algorithm:averagine_similarity_scaling`, default 0.95.
    pub averagine_similarity_scaling: f64,
    /// `algorithm:max_nr_labelled_aas`, default 0.
    pub max_nr_labelled_aas: i32,
    /// `algorithm:spectrum_type`: `automatic` (default), `profile` or
    /// `centroid`.
    pub spectrum_type: String,
    /// `algorithm:averagine_type`: `peptide` (default), `RNA` or `DNA`.
    pub averagine_type: String,
    /// `algorithm:knock_out`, default false.
    pub knock_out: bool,
    /// Section `labels:`, one mass per master-list label.
    pub label_masses: LabelMassShifts,
}

impl Default for FeatureFinderMultiplexParameters {
    fn default() -> Self {
        Self {
            labels: "[][Lys8,Arg10]".into(),
            charge: "1:4".into(),
            isotopes_per_peptide: "3:6".into(),
            rt_typical: 40.0,
            rt_band: 0.0,
            rt_min: 2.0,
            mz_tolerance: 6.0,
            mz_unit: "ppm".into(),
            intensity_cutoff: 1000.0,
            peptide_similarity: 0.5,
            averagine_similarity: 0.4,
            averagine_similarity_scaling: 0.95,
            max_nr_labelled_aas: 0,
            spectrum_type: "automatic".into(),
            averagine_type: "peptide".into(),
            knock_out: false,
            label_masses: LabelMassShifts::default(),
        }
    }
}

fn entry(name: &str, value: ParamValue, description: &str, advanced: bool) -> ParamEntry {
    let mut e = ParamEntry {
        name: name.into(),
        description: description.into(),
        value,
        ..ParamEntry::default()
    };
    if advanced {
        e.tags.insert("advanced".into());
    }
    e
}

fn float_entry(
    name: &str,
    value: f64,
    description: &str,
    advanced: bool,
    min: f64,
    max: Option<f64>,
) -> ParamEntry {
    let mut e = entry(name, ParamValue::Float(value), description, advanced);
    e.min_float = min;
    if let Some(max) = max {
        e.max_float = max;
    }
    e
}

fn string_entry(
    name: &str,
    value: &str,
    description: &str,
    advanced: bool,
    valid: &[&str],
) -> ParamEntry {
    let mut e = entry(
        name,
        ParamValue::String(value.into()),
        description,
        advanced,
    );
    e.valid_strings = valid.iter().map(|v| (*v).to_owned()).collect();
    e
}

fn float(param: &Param, key: &str) -> Result<f64> {
    param.value(key)?.to_f64()
}

fn string(param: &Param, key: &str) -> Result<String> {
    Ok(param.value(key)?.as_str()?.to_owned())
}

impl FeatureFinderMultiplexParameters {
    /// The parameter tree, sections `algorithm` and `labels`, carrying these
    /// values.
    ///
    /// # Errors
    ///
    /// The errors of [`LabelMassShifts::to_param`], and a parameter-tree
    /// resource failure.
    pub fn to_param(&self) -> Result<Param> {
        let algorithm = vec![
            string_entry(
                "labels",
                &self.labels,
                "Labels used for labelling the samples. If the sample is unlabelled (i.e. you want to detect only single peptide features) please leave this parameter empty. [...] specifies the labels for a single sample. For example\n\n[][Lys8,Arg10]        ... SILAC\n[][Lys4,Arg6][Lys8,Arg10]        ... triple-SILAC\n[Dimethyl0][Dimethyl6]        ... Dimethyl\n[Dimethyl0][Dimethyl4][Dimethyl8]        ... triple Dimethyl\n[ICPL0][ICPL4][ICPL6][ICPL10]        ... ICPL",
                false,
                &[],
            ),
            string_entry(
                "charge",
                &self.charge,
                "Range of charge states in the sample, i.e. min charge : max charge.",
                false,
                &[],
            ),
            string_entry(
                "isotopes_per_peptide",
                &self.isotopes_per_peptide,
                "Range of isotopes per peptide in the sample. For example 3:6, if isotopic peptide patterns in the sample consist of either three, four, five or six isotopic peaks. ",
                true,
                &[],
            ),
            float_entry(
                "rt_typical",
                self.rt_typical,
                "Typical retention time [s] over which a characteristic peptide elutes. (This is not an upper bound. Peptides that elute for longer will be reported.)",
                false,
                0.0,
                None,
            ),
            float_entry(
                "rt_band",
                self.rt_band,
                "The algorithm searches for characteristic isotopic peak patterns, spectrum by spectrum. For some low-intensity peptides, an important peak might be missing in one spectrum but be present in one of the neighbouring ones. The algorithm takes a bundle of neighbouring spectra with width rt_band into account. For example with rt_band = 0, all characteristic isotopic peaks have to be present in one and the same spectrum. As rt_band increases, the sensitivity of the algorithm but also the likelihood of false detections increases.",
                false,
                0.0,
                None,
            ),
            float_entry(
                "rt_min",
                self.rt_min,
                "Lower bound for the retention time [s]. (Any peptides seen for a shorter time period are not reported.)",
                false,
                0.0,
                None,
            ),
            float_entry(
                "mz_tolerance",
                self.mz_tolerance,
                "m/z tolerance for search of peak patterns.",
                false,
                0.0,
                None,
            ),
            string_entry(
                "mz_unit",
                &self.mz_unit,
                "Unit of the 'mz_tolerance' parameter.",
                false,
                &["Da", "ppm"],
            ),
            float_entry(
                "intensity_cutoff",
                self.intensity_cutoff,
                "Lower bound for the intensity of isotopic peaks.",
                false,
                0.0,
                None,
            ),
            float_entry(
                "peptide_similarity",
                self.peptide_similarity,
                "Two peptides in a multiplet are expected to have the same isotopic pattern. This parameter is a lower bound on their similarity.",
                false,
                -1.0,
                Some(1.0),
            ),
            float_entry(
                "averagine_similarity",
                self.averagine_similarity,
                "The isotopic pattern of a peptide should resemble the averagine model at this m/z position. This parameter is a lower bound on similarity between measured isotopic pattern and the averagine model.",
                false,
                -1.0,
                Some(1.0),
            ),
            float_entry(
                "averagine_similarity_scaling",
                self.averagine_similarity_scaling,
                "Let x denote this scaling factor, and p the averagine similarity parameter. For the detection of single peptides, the averagine parameter p is replaced by p' = p + x(1-p), i.e. x = 0 -> p' = p and x = 1 -> p' = 1. (For knock_out = true, peptide doublets and singlets are detected simultaneously. For singlets, the peptide similarity filter is irreleavant. In order to compensate for this 'missing filter', the averagine parameter p is replaced by the more restrictive p' when searching for singlets.)",
                true,
                0.0,
                Some(1.0),
            ),
            {
                let mut e = entry(
                    "max_nr_labelled_aas",
                    ParamValue::Integer(i64::from(self.max_nr_labelled_aas)),
                    "Maximum number of labelled amino acids per peptide, minus one. The algorithm searches for peptides with up to (this value + 1) labelled amino acids. For SILAC with trypsin digestion, this parameter corresponds to the maximum number of missed cleavages.",
                    false,
                );
                e.min_int = 0;
                e
            },
            string_entry(
                "spectrum_type",
                &self.spectrum_type,
                "Type of MS1 spectra in input mzML file. 'automatic' determines the spectrum type directly from the input mzML file.",
                true,
                &["profile", "centroid", "automatic"],
            ),
            string_entry(
                "averagine_type",
                &self.averagine_type,
                "The type of averagine to use, currently RNA, DNA or peptide",
                true,
                &["peptide", "RNA", "DNA"],
            ),
            string_entry(
                "knock_out",
                if self.knock_out { "true" } else { "false" },
                "Is it likely that knock-outs are present? (Supported for doublex, triplex and quadruplex experiments only.)",
                true,
                &["true", "false"],
            ),
        ];
        let labels_param = self.label_masses.to_param()?;
        let labels = labels_param
            .root()
            .entries
            .iter()
            .map(|e| {
                let mut e = e.clone();
                e.tags.insert("advanced".into());
                e
            })
            .collect();
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries: Vec::new(),
            nodes: vec![
                ParamNode {
                    name: "algorithm".into(),
                    description: "algorithmic parameters".into(),
                    entries: algorithm,
                    nodes: Vec::new(),
                },
                ParamNode {
                    name: "labels".into(),
                    description: "mass shifts for all possible labels".into(),
                    entries: labels,
                    nodes: Vec::new(),
                },
            ],
        })
    }

    fn from_complete_param(param: &Param) -> Result<Self> {
        let mut masses = BTreeMap::new();
        for label in label_master_list() {
            masses.insert(
                label.short_name.clone(),
                float(param, &format!("labels:{}", label.short_name))?,
            );
        }
        let max_nr = param.value("algorithm:max_nr_labelled_aas")?.to_i64()?;
        Ok(Self {
            labels: string(param, "algorithm:labels")?,
            charge: string(param, "algorithm:charge")?,
            isotopes_per_peptide: string(param, "algorithm:isotopes_per_peptide")?,
            rt_typical: float(param, "algorithm:rt_typical")?,
            rt_band: float(param, "algorithm:rt_band")?,
            rt_min: float(param, "algorithm:rt_min")?,
            mz_tolerance: float(param, "algorithm:mz_tolerance")?,
            mz_unit: string(param, "algorithm:mz_unit")?,
            intensity_cutoff: float(param, "algorithm:intensity_cutoff")?,
            peptide_similarity: float(param, "algorithm:peptide_similarity")?,
            averagine_similarity: float(param, "algorithm:averagine_similarity")?,
            averagine_similarity_scaling: float(param, "algorithm:averagine_similarity_scaling")?,
            max_nr_labelled_aas: i32::try_from(max_nr)
                .map_err(|_| bad("algorithm:max_nr_labelled_aas does not fit an int"))?,
            spectrum_type: string(param, "algorithm:spectrum_type")?,
            averagine_type: string(param, "algorithm:averagine_type")?,
            knock_out: string(param, "algorithm:knock_out")? == "true",
            label_masses: LabelMassShifts { masses },
        })
    }
}

/// `StringUtils::toInt32` on one side of a `min:max` range: the trimmed text
/// must be a whole 32-bit integer.
fn to_int32(text: &str, key: &str) -> Result<i32> {
    text.trim()
        .parse::<i32>()
        .map_err(|_| bad(format!("{key}: could not convert '{text}' to an integer")))
}

/// Parse `min:max` as the source does: the text before the first `:`
/// (`StringUtils::prefix`) and after the last `:` (`StringUtils::suffix`),
/// swapped when reversed.
fn parse_range(text: &str, key: &str) -> Result<(i32, i32)> {
    let (Some(first), Some(last)) = (text.find(':'), text.rfind(':')) else {
        return Err(bad(format!("{key}: '{text}' holds no ':'")));
    };
    let min = to_int32(&text[..first], key)?;
    let max = to_int32(&text[last + 1..], key)?;
    Ok(if min > max { (max, min) } else { (min, max) })
}

/// Source `orderCharge`: 2+, 3+, 4+, 1+, 5+, 6+, ... by likelihood.
fn order_charge(charge: i32) -> i64 {
    if 1 < charge && charge < 5 {
        i64::from(charge) - 1
    } else if charge == 1 {
        4
    } else {
        i64::from(charge)
    }
}

/// Source `lessPattern`: more mass shifts first, then the smaller second
/// mass shift, then the more likely charge.
fn less_pattern(p1: &MultiplexIsotopicPeakPattern, p2: &MultiplexIsotopicPeakPattern) -> bool {
    if p1.mass_shift_count() == p2.mass_shift_count() {
        if p1.mass_shift_count() > 1 && p2.mass_shift_count() > 1 {
            let (a, b) = (p1.mass_shift_at(1), p2.mass_shift_at(1));
            if a == b {
                order_charge(p1.charge()) < order_charge(p2.charge())
            } else {
                a < b
            }
        } else {
            order_charge(p1.charge()) < order_charge(p2.charge())
        }
    } else {
        p1.mass_shift_count() > p2.mass_shift_count()
    }
}

/// Source `Math::LinearRegressionWithoutIntercept`: the least-squares slope of
/// a line through the origin, `sum(x*y) / sum(x*x)`, NaN for fewer than two
/// points. Private to this module: the regression header is not ported here.
fn slope_without_intercept(x: &[f64], y: &[f64]) -> f64 {
    let mut sum_xx = 0.0;
    let mut sum_xy = 0.0;
    for (&xi, &yi) in x.iter().zip(y) {
        sum_xx += xi * xi;
        sum_xy += xi * yi;
    }
    if x.len() < 2 {
        return f64::NAN;
    }
    sum_xy / sum_xx
}

/// The source's `DBoundingBox<2>` of one mass trace, x = RT, y = m/z.
#[derive(Clone, Copy)]
struct TraceBox {
    min: (f64, f64),
    max: (f64, f64),
    empty: bool,
}

impl TraceBox {
    fn new() -> Self {
        Self {
            min: (f64::MAX, f64::MAX),
            max: (f64::MIN, f64::MIN),
            empty: true,
        }
    }

    fn enlarge(&mut self, x: f64, y: f64) {
        self.empty = false;
        if x < self.min.0 {
            self.min.0 = x;
        }
        if x > self.max.0 {
            self.max.0 = x;
        }
        if y < self.min.1 {
            self.min.1 = y;
        }
        if y > self.max.1 {
            self.max.1 = y;
        }
    }

    /// The source widens a box of zero width or height by 0.01 on every side,
    /// then adds its four corners as a convex hull unless it is empty.
    fn hull(mut self) -> Result<Option<ConvexHull2D>> {
        if self.empty {
            return Ok(None);
        }
        if self.max.0 - self.min.0 == 0.0 || self.max.1 - self.min.1 == 0.0 {
            let (min, max) = (self.min, self.max);
            self.enlarge(min.0 - 0.01, min.1 - 0.01);
            self.enlarge(max.0 + 0.01, max.1 + 0.01);
        }
        let mut hull = ConvexHull2D::new();
        hull.add_point(Point2D::new(self.min.0, self.min.1))?;
        hull.add_point(Point2D::new(self.min.0, self.max.1))?;
        hull.add_point(Point2D::new(self.max.0, self.min.1))?;
        hull.add_point(Point2D::new(self.max.0, self.max.1))?;
        Ok(Some(hull))
    }
}

/// Identification-free quantitative feature finder for SILAC, dimethyl, ICPL
/// and similar label-multiplet experiments (source
/// `FeatureFinderMultiplexAlgorithm`).
///
/// Detects pairs (or larger multiplets) of isotopic envelopes separated by the
/// fixed m/z shifts the labels cause. Works on profile and on centroided MS1
/// spectra; no prior sequence identification is needed. Three stages:
///
/// - **Filtering** keeps the data points where every peak of an isotope
///   pattern (e.g. at m/z, m/z+0.5, m/z+1, m/z+3, m/z+3.5 and m/z+4 for a
///   doubly charged pair) clears the intensity threshold, whose intensities
///   correlate across the peptides, and whose within-peptide intensities
///   resemble the averagine model.
/// - **Clustering** groups the surviving points in the RT-m/z plane on a grid
///   ([`MultiplexClustering`]); each cluster is the mono-isotopic mass trace of
///   the lightest peptide of a label pattern.
/// - **Quantification** integrates each peptide's mass traces and corrects
///   the peptide intensities with the ratios a linear fit through the origin
///   finds between corresponding mass traces.
///
/// The source's class brief says the clustering is hierarchical with the
/// number of clusters chosen by maximising the silhouette width, and that
/// linear fitting yields the relative amounts; the implementation clusters
/// on a grid ([`MultiplexClustering`]) and uses the fit only to correct the
/// integrated intensities, which is what this port follows.
///
/// The source derives from `DefaultParamHandler` and `ProgressLogger`. The
/// parameters are [`FeatureFinderMultiplexParameters`]; progress goes to the
/// reporter passed to [`Self::run_with_progress`]. Unique IDs come from the
/// caller-owned [`Self::unique_ids`], where the source draws from its
/// process-wide generator.
#[derive(Clone, Debug)]
pub struct FeatureFinderMultiplexAlgorithm {
    /// The parameters (source `param_`).
    pub parameters: FeatureFinderMultiplexParameters,
    /// Generator of the unique IDs assigned to the result maps.
    pub unique_ids: UniqueIdGenerator,
    /// Profile data (source `exp_profile_`).
    exp_profile: MSExperiment,
    /// Centroided data (source `exp_centroid_`).
    exp_centroid: MSExperiment,
    /// Whether the last run treated its input as centroided (source
    /// `centroided_`).
    centroided: bool,
    /// Whether the last run was asked to report progress (source
    /// `progress_`, stored but otherwise unused by the source).
    progress: bool,
    /// Charge range of the last run (source `charge_min_`, `charge_max_`).
    charge_min: i32,
    charge_max: i32,
    /// Isotopes per peptide of the last run (source
    /// `isotopes_per_peptide_min_`, `isotopes_per_peptide_max_`).
    isotopes_per_peptide_min: usize,
    isotopes_per_peptide_max: usize,
    /// Label to mass shift (source `label_mass_shift_`).
    label_mass_shift: BTreeMap<String, f64>,
    /// Detected peptide features (source `feature_map_`).
    feature_map: FeatureMap,
    /// Detected multiplets (source `consensus_map_`).
    consensus_map: ConsensusMap,
    /// Blacklisted peaks of a profile run (source `exp_blacklist_`).
    exp_blacklist: MSExperiment,
}

impl Default for FeatureFinderMultiplexAlgorithm {
    fn default() -> Self {
        Self::with_parameters(FeatureFinderMultiplexParameters::default())
    }
}

impl FeatureFinderMultiplexAlgorithm {
    /// An algorithm with the default parameters (source constructor).
    pub fn new() -> Self {
        Self::default()
    }

    /// An algorithm with `parameters`.
    pub fn with_parameters(parameters: FeatureFinderMultiplexParameters) -> Self {
        let label_mass_shift = parameters.label_masses.masses.clone();
        Self {
            parameters,
            unique_ids: UniqueIdGenerator::new(),
            exp_profile: MSExperiment::default(),
            exp_centroid: MSExperiment::default(),
            centroided: false,
            progress: false,
            // The source parses the default ranges in its constructor.
            charge_min: 1,
            charge_max: 4,
            isotopes_per_peptide_min: 3,
            isotopes_per_peptide_max: 6,
            label_mass_shift,
            feature_map: FeatureMap::default(),
            consensus_map: ConsensusMap::default(),
            exp_blacklist: MSExperiment::default(),
        }
    }

    /// The default parameters (source `getDefaults`): section `algorithm`
    /// ("algorithmic parameters") and section `labels` ("mass shifts for all
    /// possible labels", every label advanced and restricted to `0:`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        FeatureFinderMultiplexParameters::default().to_param()
    }

    /// The current parameters (source `getParameters`).
    ///
    /// # Errors
    ///
    /// As [`FeatureFinderMultiplexParameters::to_param`].
    pub fn to_param(&self) -> Result<Param> {
        self.parameters.to_param()
    }

    /// An algorithm from `parameters`, missing ones at their defaults (source
    /// constructor followed by `setParameters`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a value of the wrong type or outside
    /// its restriction (source `Exception::InvalidParameter`).
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`Self::from_param`], also returning the unknown-parameter warnings
    /// the source logs.
    ///
    /// # Errors
    ///
    /// As [`Self::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(FEATURE_FINDER_MULTIPLEX_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        let (parameters, warnings) = handler.set_parameters_with(
            parameters,
            FeatureFinderMultiplexParameters::from_complete_param,
        )?;
        Ok((Self::with_parameters(parameters), warnings))
    }

    /// The features of the most recent run (source `getFeatureMap`); empty
    /// before [`Self::run`].
    pub fn feature_map(&self) -> &FeatureMap {
        &self.feature_map
    }

    /// The multiplets of the most recent run (source `getConsensusMap`); empty
    /// before [`Self::run`].
    pub fn consensus_map(&self) -> &ConsensusMap {
        &self.consensus_map
    }

    /// The blacklisted peaks of the most recent profile run (source
    /// `getBlacklist`); a centroided run leaves it unchanged.
    pub fn blacklist(&self) -> &MSExperiment {
        &self.exp_blacklist
    }

    /// Whether the last run treated its input as centroided (source
    /// `centroided_`).
    pub fn is_centroided(&self) -> bool {
        self.centroided
    }

    /// The `progress` flag of the last [`Self::run`] (source `progress_`).
    pub fn reports_progress(&self) -> bool {
        self.progress
    }

    /// The charge range of the last run, minimum first (source `charge_min_`,
    /// `charge_max_`; before a run, the defaults `1:4`).
    pub fn charge_range(&self) -> (i32, i32) {
        (self.charge_min, self.charge_max)
    }

    /// The isotopes-per-peptide range of the last run (source
    /// `isotopes_per_peptide_min_`, `isotopes_per_peptide_max_`).
    pub fn isotopes_per_peptide_range(&self) -> (usize, usize) {
        (self.isotopes_per_peptide_min, self.isotopes_per_peptide_max)
    }

    /// Label to mass shift (source `label_mass_shift_`).
    pub fn label_mass_shift(&self) -> &BTreeMap<String, f64> {
        &self.label_mass_shift
    }

    /// The profile buffer (source `exp_profile_`).
    pub fn profile_experiment(&self) -> &MSExperiment {
        &self.exp_profile
    }

    /// The centroided buffer (source `exp_centroid_`): after a run, the
    /// intensity-filtered centroided data the satellites refer to.
    pub fn centroided_experiment(&self) -> &MSExperiment {
        &self.exp_centroid
    }

    /// Run the pipeline on `exp` (source `run(MSExperiment&, bool)`).
    ///
    /// # Errors
    ///
    /// As [`Self::run_with_progress`].
    pub fn run(&mut self, exp: &mut MSExperiment, progress: bool) -> Result<()> {
        self.progress = progress;
        self.run_with_progress(exp, &mut ProgressReporter::silent())
    }

    /// Run the full detection pipeline on `exp` and fill the feature map, the
    /// consensus map and (profile input) the blacklist.
    ///
    /// 1. Parses `algorithm:charge` and `algorithm:isotopes_per_peptide` as
    ///    `min:max` (the text before the first and after the last `:`), swapped
    ///    when reversed, and reads the `labels:` masses.
    /// 2. Clears `exp`'s chromatograms and sorts its spectra by RT (the Release
    ///    build's `std::sort`) and each spectrum's peaks by m/z, then **swaps**
    ///    it into the internal centroided or profile buffer. A fresh algorithm
    ///    leaves `exp` empty; a reused one hands back what its buffer held, as
    ///    the source's `swap` does.
    /// 3. Determines the spectrum type from `algorithm:spectrum_type`:
    ///    `automatic` asks the first spectrum ([`crate::kernel::MSSpectrum::get_type`]
    ///    with data query) and treats anything but centroid as profile.
    /// 4. Picks profile input with [`PeakPickerHiRes`] on MS level 1 with
    ///    `signal_to_noise` 0, in the picker's source-compatible mode
    ///    ([`PickingCompatibility::source`]), since the source picker accepts
    ///    input (equal m/z among others) that the native default refuses.
    /// 5. Builds the mass shift patterns ([`MultiplexDeltaMassesGenerator`],
    ///    with knock-outs when `algorithm:knock_out` is true) and the isotopic
    ///    peak patterns for every charge, sorted by likelihood.
    /// 6. Filters ([`MultiplexFilteringCentroided`] or
    ///    [`MultiplexFilteringProfile`]), clusters ([`MultiplexClustering`]) and
    ///    builds the maps. Without labels, the singlet averagine scaling is 0.
    /// 7. Sorts the consensus map by position (experiment type `labeled_MS1`),
    ///    assigns unique IDs, gives every column header present its
    ///    `channel_id`, then labels one column per sample — `label i` with
    ///    knock-outs, otherwise the sample's labels concatenated. Sorts the
    ///    feature map by position and assigns unique IDs.
    ///
    /// The maps are **not** cleared first: a second run appends to them, as in
    /// the source.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] "Error: No MS1 spectra in input file." when
    ///   `exp` has no spectra (source `Exception::FileEmpty`).
    /// - [`Error::InvalidValue`] for a range without `:` or with a side that is
    ///   not an integer (source `Exception::ElementNotFound` and
    ///   `Exception::ConversionError`), and for a negative isotope count or a
    ///   zero charge, which the source stores in an `unsigned` and wraps, or
    ///   divides by.
    /// - The errors of every stage.
    pub fn run_with_progress(
        &mut self,
        exp: &mut MSExperiment,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<()> {
        let parameters = self.parameters.clone();
        let (charge_min, charge_max) = parse_range(&parameters.charge, "algorithm:charge")?;
        let (ipp_min, ipp_max) = parse_range(
            &parameters.isotopes_per_peptide,
            "algorithm:isotopes_per_peptide",
        )?;
        let ipp_min = usize::try_from(ipp_min)
            .map_err(|_| bad("isotopes per peptide must not be negative"))?;
        let ipp_max = usize::try_from(ipp_max)
            .map_err(|_| bad("isotopes per peptide must not be negative"))?;
        let ipp_max_i32 =
            i32::try_from(ipp_max).map_err(|_| bad("isotopes per peptide does not fit an int"))?;
        self.charge_min = charge_min;
        self.charge_max = charge_max;
        self.isotopes_per_peptide_min = ipp_min;
        self.isotopes_per_peptide_max = ipp_max;
        for (label, &mass) in &parameters.label_masses.masses {
            self.label_mass_shift.insert(label.clone(), mass);
        }

        if exp.spectra.is_empty() {
            return Err(bad("Error: No MS1 spectra in input file."));
        }
        exp.chromatograms.clear();
        for spectrum in &exp.spectra {
            if !spectrum.rt.is_finite() {
                return Err(bad("spectrum retention time is not finite"));
            }
        }
        source_sort_by(&mut exp.spectra, |a, b| a.rt < b.rt)?;
        for spectrum in &mut exp.spectra {
            spectrum.sort_by_position()?;
        }

        self.centroided = match parameters.spectrum_type.as_str() {
            "automatic" => exp.spectra[0].get_type(true)? == SpectrumType::Centroid,
            "centroid" => true,
            _ => false,
        };
        if self.centroided {
            std::mem::swap(exp, &mut self.exp_centroid);
        } else {
            std::mem::swap(exp, &mut self.exp_profile);
        }

        let mut boundaries_exp_s: Vec<Vec<PeakBoundary>> = Vec::new();
        if !self.centroided {
            // The source picks whatever the data holds; the port's picker
            // refuses some inputs (equal m/z, for one) unless it is asked to
            // reproduce the source.
            let picker = PeakPickerHiRes {
                ms_levels: vec![1],
                signal_to_noise: 0.0,
                compatibility: PickingCompatibility::source(),
                ..PeakPickerHiRes::default()
            };
            let picked = picker.pick_experiment(&self.exp_profile)?;
            self.exp_centroid = picked.experiment;
            // The source appends boundaries for picked spectra only.
            boundaries_exp_s = picked.spectrum_boundaries.into_iter().flatten().collect();
        }

        let mut generator = MultiplexDeltaMassesGenerator::new(
            &parameters.labels,
            parameters.max_nr_labelled_aas,
            self.label_mass_shift.clone(),
        )?;
        if parameters.knock_out {
            generator.generate_knockout_delta_masses()?;
        }
        let patterns = generate_peak_patterns(
            charge_min,
            charge_max,
            ipp_max_i32,
            generator.delta_masses_list(),
        )?;
        let samples = generator.samples_labels_list();
        let averagine_similarity_scaling =
            if samples.len() == 1 && samples[0].len() == 1 && samples[0][0] == "no_label" {
                0.0
            } else {
                parameters.averagine_similarity_scaling
            };
        let ppm = parameters.mz_unit == "ppm";

        if self.centroided {
            let mut filtering = MultiplexFilteringCentroided::new(
                &self.exp_centroid,
                patterns.clone(),
                ipp_min,
                ipp_max,
                parameters.intensity_cutoff,
                parameters.rt_band,
                parameters.mz_tolerance,
                ppm,
                parameters.peptide_similarity,
                parameters.averagine_similarity,
                averagine_similarity_scaling,
                &parameters.averagine_type,
            )?;
            let filter_results = filtering.filter_with_progress(progress)?;
            let clustering = MultiplexClustering::from_centroided(
                &self.exp_centroid,
                parameters.mz_tolerance,
                ppm,
                parameters.rt_typical,
            )?;
            let cluster_results = clustering.cluster_with_progress(&filter_results, progress)?;
            // `filtering.getCentroidedExperiment().swap(exp_centroid_)`
            self.exp_centroid = filtering.centroided_experiment().clone();
            self.generate_maps_centroided(&patterns, &filter_results, &cluster_results)?;
        } else {
            let mut filtering = MultiplexFilteringProfile::new(
                &self.exp_profile,
                &self.exp_centroid,
                &boundaries_exp_s,
                patterns.clone(),
                ipp_min,
                ipp_max,
                parameters.intensity_cutoff,
                parameters.rt_band,
                parameters.mz_tolerance,
                ppm,
                parameters.peptide_similarity,
                parameters.averagine_similarity,
                averagine_similarity_scaling,
                &parameters.averagine_type,
            )?;
            let filter_results = filtering.filter_with_progress(progress)?;
            self.exp_blacklist = filtering.blacklist_experiment();
            let clustering = MultiplexClustering::from_profile(
                &self.exp_profile,
                &self.exp_centroid,
                &boundaries_exp_s,
                parameters.rt_typical,
            )?;
            let cluster_results = clustering.cluster_with_progress(&filter_results, progress)?;
            self.exp_centroid = filtering.centroided_experiment().clone();
            self.generate_maps_profile(&patterns, &filter_results, &cluster_results, progress)?;
        }

        self.consensus_map.experiment_type = "labeled_MS1".into();
        self.consensus_map.sort_by_position()?;
        self.consensus_map.unique_id = self.unique_ids.get_unique_id();
        for feature in &mut self.consensus_map.features {
            feature.unique_id = self.unique_ids.get_unique_id();
        }
        for (i, header) in self.consensus_map.column_headers.values_mut().enumerate() {
            header
                .metadata
                .insert("channel_id".into(), MetaValue::from(i as i64));
        }
        for (i, sample) in samples.iter().enumerate() {
            let header = self
                .consensus_map
                .column_headers
                .entry(i as u64)
                .or_default();
            header.label = if parameters.knock_out {
                format!("label {i}")
            } else {
                sample.concat()
            };
        }
        source_sort_by(&mut self.feature_map.features, |a, b| {
            a.rt < b.rt || (a.rt == b.rt && a.mz < b.mz)
        })?;
        self.feature_map.unique_id = self.unique_ids.get_unique_id();
        for feature in &mut self.feature_map.features {
            feature.unique_id = self.unique_ids.get_unique_id();
        }
        Ok(())
    }

    /// Peak patterns for every charge from `charge_max` down to `charge_min`
    /// and every mass shift pattern, sorted by likelihood (source protected
    /// `generatePeakPatterns_`); see [`generate_peak_patterns`].
    ///
    /// # Errors
    ///
    /// As [`generate_peak_patterns`].
    pub fn generate_peak_patterns(
        &self,
        charge_min: i32,
        charge_max: i32,
        peaks_per_peptide_max: i32,
        mass_pattern_list: &[MultiplexDeltaMasses],
    ) -> Result<Vec<MultiplexIsotopicPeakPattern>> {
        generate_peak_patterns(
            charge_min,
            charge_max,
            peaks_per_peptide_max,
            mass_pattern_list,
        )
    }

    /// Correct the peptide intensities with the ratios between their mass
    /// traces (source protected `correctPeptideIntensities_`).
    ///
    /// Labelled experiments are about fold changes, so the ratios are
    /// determined from interpolated chromatograms directly. For each heavier
    /// peptide, every isotope whose trace has a spline in both that peptide and
    /// the lightest one is scanned from the smaller to the larger end of the
    /// two splines in steps of 0.7 times the smaller step width; the lightest
    /// peptide's trace is read at `rt` and the heavier one's at `rt +
    /// rt_peptide[p] - rt_peptide[0]`, and pairs with both intensities positive
    /// are kept. More than five pairs give a ratio, the slope of a fit through
    /// the origin; otherwise the ratio is -1. A pair's corrected intensities
    /// are `I0 = (I0 + r I1) / (1 + r^2)` and `I1 = r I0`; for more peptides,
    /// each with a positive ratio becomes `r * I0`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] when a scan would take more than
    /// [`MAX_RATIO_SCAN_STEPS`] steps, and [`Error::InvalidValue`] when
    /// `rt_peptide` or `intensity_peptide` is shorter than the pattern's
    /// peptides (undefined in the source); and the spline evaluation errors.
    pub fn correct_peptide_intensities(
        &self,
        pattern: &MultiplexIsotopicPeakPattern,
        spline_chromatograms: &BTreeMap<usize, SplinePackage>,
        rt_peptide: &[f64],
        intensity_peptide: &mut [f64],
    ) -> Result<()> {
        let count = pattern.mass_shift_count();
        if rt_peptide.len() < count || intensity_peptide.len() < count {
            return Err(bad("fewer peptide RTs or intensities than peptides"));
        }
        let ipp_max = self.isotopes_per_peptide_max;
        let mut ratios = vec![1.0];
        for peptide in 1..count {
            let mut intensities1 = Vec::new();
            let mut intensities2 = Vec::new();
            for isotope in 0..ipp_max {
                let (Some(s1), Some(s2)) = (
                    spline_chromatograms.get(&isotope),
                    spline_chromatograms.get(&(peptide * ipp_max + isotope)),
                ) else {
                    continue;
                };
                let rt_min = s1.pos_min().min(s2.pos_min());
                let rt_max = s1.pos_max().max(s2.pos_max());
                let step = 0.7 * s1.pos_step_width().min(s2.pos_step_width());
                let steps = (rt_max - rt_min) / step;
                if !(steps <= MAX_RATIO_SCAN_STEPS as f64) {
                    return Err(Error::InvalidRange(format!(
                        "intensity correction: more than {MAX_RATIO_SCAN_STEPS} scan steps"
                    )));
                }
                let mut rt = rt_min;
                while rt < rt_max {
                    let intensity1 = s1.eval(rt)?;
                    let intensity2 = s2.eval(rt + rt_peptide[peptide] - rt_peptide[0])?;
                    if intensity1 > 0.0 && intensity2 > 0.0 {
                        intensities1.push(intensity1);
                        intensities2.push(intensity2);
                    }
                    rt += step;
                }
            }
            if intensities1.len() > 5 {
                ratios.push(slope_without_intercept(&intensities1, &intensities2));
            } else {
                ratios.push(-1.0);
            }
        }
        if count == 2 && ratios[1] > 0.0 {
            let intensity0 = (intensity_peptide[0] + ratios[1] * intensity_peptide[1])
                / (1.0 + ratios[1] * ratios[1]);
            let intensity1 = ratios[1] * intensity0;
            intensity_peptide[0] = intensity0;
            intensity_peptide[1] = intensity1;
        } else if count > 2 {
            for peptide in 1..count {
                if ratios[peptide] > 0.0 {
                    intensity_peptide[peptide] = ratios[peptide] * intensity_peptide[0];
                }
            }
        }
        Ok(())
    }

    /// The shared tail of both `determinePeptideIntensities*_`: per peptide,
    /// from `(rt, intensity)` points per isotope (already in the order the
    /// source builds its chromatogram), the intensity-weighted RT and the
    /// integrated intensity, or -1 for a peptide without intensity or eluting
    /// shorter than `algorithm:rt_min`.
    fn peptide_intensities(
        &self,
        pattern: &MultiplexIsotopicPeakPattern,
        mut chromatogram_of: impl FnMut(
            usize,
            &mut f64,
            &mut f64,
            &mut f64,
            &mut f64,
            bool,
        ) -> Result<MSChromatogram>,
    ) -> Result<Vec<f64>> {
        let integrator = PeakIntegrator {
            integration_method: IntegrationMethod::Trapezoid,
            ..PeakIntegrator::default()
        };
        let ipp_max = self.isotopes_per_peptide_max;
        let mut rt_peptide = Vec::new();
        let mut intensity_peptide = Vec::new();
        let mut spline_chromatograms = BTreeMap::new();
        for peptide in 0..pattern.mass_shift_count() {
            let mut rt = 0.0;
            let mut rt_min = 0.0;
            let mut rt_max = 0.0;
            let mut intensity_sum_simple = 0.0;
            let mut intensity_sum = 0.0;
            for isotope in 0..ipp_max {
                let idx = peptide * ipp_max + isotope;
                let first = peptide + isotope == 0;
                let chromatogram = chromatogram_of(
                    idx,
                    &mut rt,
                    &mut rt_min,
                    &mut rt_max,
                    &mut intensity_sum_simple,
                    first,
                )?;
                if chromatogram.peaks.len() > 5 {
                    let positions: Vec<f64> = chromatogram.peaks.iter().map(|p| p.rt).collect();
                    let intensities: Vec<f64> = chromatogram
                        .peaks
                        .iter()
                        .map(|p| f64::from(p.intensity))
                        .collect();
                    // Equal RTs make the source's spline NaN everywhere, so it
                    // evaluates to zero and contributes no ratio pair; leaving
                    // it out has the same effect.
                    if positions.windows(2).all(|w| w[0] < w[1]) {
                        spline_chromatograms
                            .insert(idx, SplinePackage::new(&positions, &intensities)?);
                    }
                }
                if chromatogram.peaks.len() > 1 {
                    let rt_start = chromatogram.peaks[0].rt;
                    let rt_end = chromatogram.peaks[chromatogram.peaks.len() - 1].rt;
                    intensity_sum += integrator
                        .integrate_chromatogram(&chromatogram, rt_start, rt_end)?
                        .area;
                } else if chromatogram.peaks.len() == 1 {
                    intensity_sum += f64::from(chromatogram.peaks[0].intensity);
                }
            }
            rt /= intensity_sum_simple;
            rt_peptide.push(rt);
            if intensity_sum == 0.0 || rt_max - rt_min < self.parameters.rt_min {
                intensity_sum = -1.0;
            }
            intensity_peptide.push(intensity_sum);
        }
        if intensity_peptide.contains(&-1.0) || pattern.mass_shift_count() < 2 {
            return Ok(intensity_peptide);
        }
        self.correct_peptide_intensities(
            pattern,
            &spline_chromatograms,
            &rt_peptide,
            &mut intensity_peptide,
        )?;
        Ok(intensity_peptide)
    }

    /// The peptide intensities of a multiplet from centroided satellites
    /// (source protected `determinePeptideIntensitiesCentroided_`).
    ///
    /// For each peptide and isotope, the satellites form a chromatogram
    /// (sorted by RT, stably) that is integrated with the trapezoid rule
    /// between its first and last point (a single point contributes its
    /// intensity); chromatograms of more than five points also get a spline
    /// for [`Self::correct_peptide_intensities`]. The peptide RT is the
    /// intensity-weighted mean RT.
    ///
    /// A peptide is -1 when its intensity is zero or its RT span is below
    /// `algorithm:rt_min`. The span reproduces the source's bookkeeping: the
    /// minimum and maximum start at zero for every peptide and are reset by
    /// **every** satellite of the lightest peptide's first isotope, so for the
    /// lightest peptide the span runs from that isotope's last satellite
    /// outward, and for the others from zero to the latest RT.
    ///
    /// When any peptide is -1, or for singlets, the intensities are returned
    /// uncorrected.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a satellite outside the centroided
    /// data, and the errors of the integration, the splines and the
    /// correction.
    pub fn determine_peptide_intensities_centroided(
        &self,
        pattern: &MultiplexIsotopicPeakPattern,
        satellites: &SatelliteMultimap<MultiplexSatelliteCentroided>,
    ) -> Result<Vec<f64>> {
        let exp = &self.exp_centroid;
        self.peptide_intensities(pattern, |idx, rt, rt_min, rt_max, sum_simple, first| {
            let mut chromatogram = MSChromatogram::default();
            for satellite in satellites.get(&idx).into_iter().flatten() {
                let spectrum = exp
                    .spectra
                    .get(satellite.rt_idx)
                    .ok_or_else(|| bad("satellite spectrum outside the centroided data"))?;
                let peak = spectrum
                    .peaks
                    .get(satellite.mz_idx)
                    .ok_or_else(|| bad("satellite peak outside its spectrum"))?;
                let rt_temp = spectrum.rt;
                let intensity_temp = f64::from(peak.intensity);
                if first || rt_temp < *rt_min {
                    *rt_min = rt_temp;
                }
                if first || rt_temp > *rt_max {
                    *rt_max = rt_temp;
                }
                *rt += rt_temp * intensity_temp;
                *sum_simple += intensity_temp;
                chromatogram
                    .peaks
                    .push(ChromatogramPeak::new(rt_temp, intensity_temp as f32));
            }
            // `sortByPosition` is a stable sort by RT.
            chromatogram.peaks.sort_by(|a, b| a.rt.total_cmp(&b.rt));
            Ok(chromatogram)
        })
    }

    /// The peptide intensities of a multiplet from profile satellites (source
    /// protected `determinePeptideIntensitiesProfile_`), as
    /// [`Self::determine_peptide_intensities_centroided`] except that each
    /// chromatogram is made unique in RT (`makePeakPositionUnique` with the
    /// median intensity) instead of only sorted.
    ///
    /// # Errors
    ///
    /// The errors of `makePeakPositionUnique`, the integration, the splines and
    /// the correction.
    pub fn determine_peptide_intensities_profile(
        &self,
        pattern: &MultiplexIsotopicPeakPattern,
        satellites: &SatelliteMultimap<MultiplexSatelliteProfile>,
    ) -> Result<Vec<f64>> {
        self.peptide_intensities(pattern, |idx, rt, rt_min, rt_max, sum_simple, first| {
            let mut chromatogram = MSChromatogram::default();
            for satellite in satellites.get(&idx).into_iter().flatten() {
                let rt_temp = f64::from(satellite.rt);
                let intensity_temp = f64::from(satellite.intensity);
                if first || rt_temp < *rt_min {
                    *rt_min = rt_temp;
                }
                if first || rt_temp > *rt_max {
                    *rt_max = rt_temp;
                }
                *rt += rt_temp * intensity_temp;
                *sum_simple += intensity_temp;
                chromatogram
                    .peaks
                    .push(ChromatogramPeak::new(rt_temp, intensity_temp as f32));
            }
            make_peak_position_unique(&mut chromatogram, IntensityAveragingMethod::Median)?;
            Ok(chromatogram)
        })
    }

    /// Append one multiplet's features and consensus feature, or nothing when a
    /// peptide's mono-isotopic trace is shorter than `algorithm:rt_min`.
    ///
    /// `trace(peptide, isotope, box, weighted)` visits one mass trace. A
    /// peptide whose weight sum is not positive is skipped when
    /// `skip_without_weight` (centroided input only).
    fn add_multiplet(
        &mut self,
        pattern: &MultiplexIsotopicPeakPattern,
        peptide_intensities: &[f64],
        skip_without_weight: bool,
        mut trace: impl FnMut(usize, &mut TraceBox, &mut (f64, f64, f64)) -> Result<()>,
    ) -> Result<()> {
        let ipp_max = self.isotopes_per_peptide_max;
        let mut features = Vec::new();
        let mut consensus = ConsensusFeature::new();
        let mut abort = false;
        for peptide in 0..pattern.mass_shift_count() {
            if abort {
                break;
            }
            let mut feature = Feature::default();
            // (rt, mz, intensity_sum)
            let mut weighted = (0.0, 0.0, 0.0);
            for isotope in 0..ipp_max {
                let idx = peptide * ipp_max + isotope;
                let mut mass_trace = TraceBox::new();
                if isotope == 0 {
                    trace(idx, &mut mass_trace, &mut weighted)?;
                } else {
                    let mut ignored = (0.0, 0.0, 0.0);
                    trace(idx, &mut mass_trace, &mut ignored)?;
                }
                if let Some(hull) = mass_trace.hull()? {
                    feature.convex_hulls.push(hull);
                }
            }
            let (mut rt, mut mz, intensity_sum) = weighted;
            if skip_without_weight && intensity_sum <= 0.0 {
                continue;
            }
            rt /= intensity_sum;
            mz /= intensity_sum;
            let intensity = *peptide_intensities
                .get(peptide)
                .ok_or_else(|| bad("fewer peptide intensities than peptides"))?;
            feature.base.rt = rt;
            feature.base.mz = mz;
            feature.base.intensity = intensity as f32;
            feature.base.charge = pattern.charge();
            feature.base.quality = 1.0;
            let bounding_box = feature
                .convex_hulls
                .first()
                .and_then(ConvexHull2D::bounding_box)
                .ok_or_else(|| {
                    bad("a peptide without a mono-isotopic mass trace (undefined in the source)")
                })?;
            // the RT extent of the mono-isotopic trace (source `maxX() - minX()`)
            if bounding_box.width() < self.parameters.rt_min {
                abort = true;
                break;
            }
            features.push(feature);
            if peptide == 0 {
                consensus.base.rt = rt;
                consensus.base.mz = mz;
                consensus.base.intensity = intensity as f32;
                consensus.base.charge = pattern.charge();
                consensus.base.quality = 1.0;
            }
            let mut handle_base = BaseFeature::new(rt, mz, intensity as f32);
            handle_base.charge = pattern.charge();
            consensus.insert(FeatureHandle::new(peptide as u64, &handle_base))?;
            self.consensus_map
                .column_headers
                .entry(peptide as u64)
                .or_default()
                .size += 1;
        }
        if !abort {
            self.consensus_map.features.push(consensus);
            self.feature_map.features.extend(features);
        }
        Ok(())
    }

    /// Build the feature and consensus maps from centroided filter results and
    /// their clusters (source protected `generateMapsCentroided_`).
    ///
    /// For each pattern and cluster, the satellites of the cluster's peaks are
    /// merged; the source stops merging a peak's satellites at the first one
    /// already present (a `break`, not a `continue`), and so does the port.
    /// Multiplets with a peptide intensity of -1 are skipped. Each peptide's
    /// feature sits at the intensity-weighted RT and m/z of its mono-isotopic
    /// trace, carries one rectangular hull per mass trace (0.01 wider on every
    /// side when flat) and the corrected intensity; a peptide whose
    /// mono-isotopic trace spans less than `algorithm:rt_min` in RT drops the
    /// whole multiplet, although the column sizes of the peptides before it
    /// have been counted already, as in the source.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a cluster point or satellite outside
    /// the filter results or the data, and the errors of the intensity
    /// determination.
    pub fn generate_maps_centroided(
        &mut self,
        patterns: &[MultiplexIsotopicPeakPattern],
        filter_results: &[MultiplexFilteredMSExperiment],
        cluster_results: &[BTreeMap<i32, GridBasedCluster>],
    ) -> Result<()> {
        for (pattern_idx, pattern) in patterns.iter().enumerate() {
            let clusters = cluster_results
                .get(pattern_idx)
                .ok_or_else(|| bad("missing cluster result"))?;
            let results = filter_results
                .get(pattern_idx)
                .ok_or_else(|| bad("missing filter result"))?;
            for cluster in clusters.values() {
                let mut satellites: SatelliteMultimap<MultiplexSatelliteCentroided> =
                    SatelliteMultimap::new();
                for &point in cluster.points() {
                    let peak = usize::try_from(point)
                        .ok()
                        .and_then(|p| results.peak(p))
                        .ok_or_else(|| bad("cluster point outside the filter result"))?;
                    for (key, satellite) in peak.satellite_entries() {
                        let present = satellites
                            .values()
                            .flatten()
                            .any(|s| s.rt_idx == satellite.rt_idx && s.mz_idx == satellite.mz_idx);
                        if present {
                            break;
                        }
                        satellites.entry(key).or_default().push(*satellite);
                    }
                }
                let peptide_intensities =
                    self.determine_peptide_intensities_centroided(pattern, &satellites)?;
                if peptide_intensities.contains(&-1.0) {
                    continue;
                }
                let exp = self.exp_centroid.clone();
                self.add_multiplet(
                    pattern,
                    &peptide_intensities,
                    true,
                    |idx, mass_trace, weighted| {
                        for satellite in satellites.get(&idx).into_iter().flatten() {
                            let spectrum = exp.spectra.get(satellite.rt_idx).ok_or_else(|| {
                                bad("satellite spectrum outside the centroided data")
                            })?;
                            let peak = spectrum
                                .peaks
                                .get(satellite.mz_idx)
                                .ok_or_else(|| bad("satellite peak outside its spectrum"))?;
                            let intensity = f64::from(peak.intensity);
                            weighted.0 += spectrum.rt * intensity;
                            weighted.1 += peak.mz * intensity;
                            weighted.2 += intensity;
                            mass_trace.enlarge(spectrum.rt, peak.mz);
                        }
                        Ok(())
                    },
                )?;
            }
        }
        Ok(())
    }

    /// Build the feature and consensus maps from profile filter results and
    /// their clusters (source protected `generateMapsProfile_`), reporting
    /// progress with the label "constructing maps", one step per pattern.
    ///
    /// As [`Self::generate_maps_centroided`] with the satellite data points:
    /// all of a cluster's data points are merged without de-duplication, each
    /// mono-isotopic data point is weighted by its intensity plus 0.0001, and a
    /// peptide is never skipped for lack of weight.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a cluster point outside the filter
    /// results and for a peptide without a mono-isotopic trace (undefined in the
    /// source), and the errors of the intensity determination and `progress`.
    pub fn generate_maps_profile(
        &mut self,
        patterns: &[MultiplexIsotopicPeakPattern],
        filter_results: &[MultiplexFilteredMSExperiment],
        cluster_results: &[BTreeMap<i32, GridBasedCluster>],
        progress: &mut ProgressReporter<'_>,
    ) -> Result<()> {
        progress.start_count(patterns.len(), "constructing maps")?;
        for (pattern_idx, pattern) in patterns.iter().enumerate() {
            progress.set_count(pattern_idx + 1)?;
            let clusters = cluster_results
                .get(pattern_idx)
                .ok_or_else(|| bad("missing cluster result"))?;
            let results = filter_results
                .get(pattern_idx)
                .ok_or_else(|| bad("missing filter result"))?;
            for cluster in clusters.values() {
                let mut satellites: SatelliteMultimap<MultiplexSatelliteProfile> =
                    SatelliteMultimap::new();
                for &point in cluster.points() {
                    let peak = usize::try_from(point)
                        .ok()
                        .and_then(|p| results.peak(p))
                        .ok_or_else(|| bad("cluster point outside the filter result"))?;
                    for (key, satellite) in peak.satellite_profile_entries() {
                        satellites.entry(key).or_default().push(*satellite);
                    }
                }
                let peptide_intensities =
                    self.determine_peptide_intensities_profile(pattern, &satellites)?;
                if peptide_intensities.contains(&-1.0) {
                    continue;
                }
                self.add_multiplet(
                    pattern,
                    &peptide_intensities,
                    false,
                    |idx, mass_trace, weighted| {
                        for satellite in satellites.get(&idx).into_iter().flatten() {
                            let rt = f64::from(satellite.rt);
                            let intensity = f64::from(satellite.intensity) + 0.0001;
                            weighted.0 += rt * intensity;
                            weighted.1 += satellite.mz * intensity;
                            weighted.2 += intensity;
                            mass_trace.enlarge(rt, satellite.mz);
                        }
                        Ok(())
                    },
                )?;
            }
        }
        progress.end()
    }
}

/// Peak patterns for every charge from `charge_max` down to `charge_min` and
/// every mass shift pattern `i` (mass shift index `i`), sorted with the
/// Release build's `std::sort` by the source's `lessPattern`: more peptides
/// first, then the smaller second mass shift, then charge by likelihood (2+,
/// 3+, 4+, 1+, 5+, 6+, ...).
///
/// # Errors
///
/// The errors of [`MultiplexIsotopicPeakPattern::new`] (a zero charge among
/// them), and [`Error::InvalidValue`] for more mass shift patterns than an
/// `int` index holds.
pub fn generate_peak_patterns(
    charge_min: i32,
    charge_max: i32,
    peaks_per_peptide_max: i32,
    mass_pattern_list: &[MultiplexDeltaMasses],
) -> Result<Vec<MultiplexIsotopicPeakPattern>> {
    let mut list = Vec::new();
    for c in (charge_min..=charge_max).rev() {
        for (i, masses) in mass_pattern_list.iter().enumerate() {
            let index = i32::try_from(i).map_err(|_| bad("too many mass shift patterns"))?;
            list.push(MultiplexIsotopicPeakPattern::new(
                c,
                peaks_per_peptide_max,
                masses.clone(),
                index,
            )?);
        }
    }
    source_sort_by(&mut list, less_pattern)?;
    Ok(list)
}
