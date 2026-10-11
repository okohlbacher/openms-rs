// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Filtering, annotating, picking, scoring and selecting spectra (for example
//! from a DDA experiment) against a target list, and matching them against a
//! spectral library.
//!
//! Ports `ANALYSIS/OPENSWATH/TargetedSpectraExtractor.h` with
//! `TargetedSpectraExtractor.cpp` at core `bc9cc12`; see
//! `docs/TARGETED_SPECTRA_EXTRACTOR_SUPPORT.md` for the member mapping, the
//! gaps and the evidence.
//!
//! The input experiment is a standard mzML file, the target list a TraML-style
//! assay list. Only spectra matching the target list's transitions by
//! precursor m/z and retention time are kept and annotated with the
//! transition's name. The spectra are smoothed and their peaks picked, then
//! scored by total ion current, FWHM and signal-to-noise ratio, and one
//! spectrum is chosen for every transition for which at least one valid
//! spectrum was found.
//!
//! # Workflow
//!
//! [`TargetedSpectraExtractor::extract_spectra`](crate::analysis::openswath::targeted_spectra_extractor::TargetedSpectraExtractor::extract_spectra)
//! runs the whole pipeline; the steps can also be run individually, in this order:
//!
//! 1. [`annotate_spectra`](crate::analysis::openswath::targeted_spectra_extractor::TargetedSpectraExtractor::annotate_spectra): match
//!    spectra to targets by precursor m/z and retention time;
//! 2. [`pick_spectrum`](crate::analysis::openswath::targeted_spectra_extractor::TargetedSpectraExtractor::pick_spectrum): smooth and
//!    pick peaks, for each spectrum;
//! 3. [`score_spectra`](crate::analysis::openswath::targeted_spectra_extractor::TargetedSpectraExtractor::score_spectra): score by TIC,
//!    FWHM and S/N with configurable weights;
//! 4. [`select_spectra`](crate::analysis::openswath::targeted_spectra_extractor::TargetedSpectraExtractor::select_spectra): keep the
//!    best spectrum per transition.
//!
//! # Spectral library matching
//!
//! [`match_spectrum`](crate::analysis::openswath::targeted_spectra_extractor::TargetedSpectraExtractor::match_spectrum) compares one
//! spectrum against a library held by a [`Comparator`](crate::analysis::openswath::targeted_spectra_extractor::Comparator);
//! [`targeted_matching`](crate::analysis::openswath::targeted_spectra_extractor::TargetedSpectraExtractor::targeted_matching) and
//! [`untargeted_matching`](crate::analysis::openswath::targeted_spectra_extractor::TargetedSpectraExtractor::untargeted_matching) do it
//! in bulk, for extracted and for unprocessed spectra respectively.
//!
//! Downstream, `MRMFeatureFilter` filters the results for QC and
//! `AbsoluteQuantitation` quantifies them; neither is ported yet.
//!
//! # Not ported
//!
//! - `searchSpectrum`: it runs `AccurateMassSearchEngine`, which is not in the
//!   crate. Its `AccurateMassSearchEngine:` parameter subsection is carried in
//!   the parameters all the same.
//! - `storeSpectraMSP`: it writes through `FileHandler` to the MSP format,
//!   whose writer (`MSPGenericFile::store`) is not in the crate. The two
//!   preparation steps it runs first are ported as
//!   [`deisotope_ms2_spectra`](crate::analysis::openswath::targeted_spectra_extractor::TargetedSpectraExtractor::deisotope_ms2_spectra)
//!   and [`remove_ms2_spectra_peaks`](crate::analysis::openswath::targeted_spectra_extractor::TargetedSpectraExtractor::remove_ms2_spectra_peaks).

use super::reaction_monitoring_transition::ReactionMonitoringTransition;
use super::targeted_experiment::TargetedExperiment;
use super::targeted_experiment_helper::{Peptide, Protein, RTUnit, RetentionTime};
use crate::comparison::{
    BinConfig, BinUnit, BinnedSpectralContrastAngle, BinnedSpectrum, BinnedSpectrumCompareFunctor,
    Tolerance,
};
use crate::concept::fuzzy_string_comparator::format_g;
use crate::concept::math_functions::ppm_to_mass;
use crate::concept::unique_id::UniqueIdGenerator;
use crate::kernel::features::{BaseFeature, Feature, FeatureMap};
use crate::kernel::{DataArray, MSExperiment, MSSpectrum};
use crate::math::source_sort::source_sort_by;
use crate::metadata::{MetaInfo, MetaValue, MetaValueData};
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::processing::SpectrumFilter;
use crate::processing::deisotoping::Deisotoper;
use crate::processing::peak_picking::{
    FwhmUnit, PeakPickerHiRes, PickingCompatibility, SignalToNoiseEstimatorMedian,
};
use crate::processing::smoothing::{GaussFilter, GaussianWidth, SavitzkyGolayFilter};
use crate::{Error, Result};
use std::collections::BTreeMap;

/// Parameter handler name of the source class.
pub const TARGETED_SPECTRA_EXTRACTOR_NAME: &str = "TargetedSpectraExtractor";

/// The parameter subsections the source registers in `subsections_`: their
/// keys are passed to the filter or picker unchecked.
pub const SUBSECTIONS: [&str; 3] = ["SavitzkyGolayFilter", "GaussFilter", "PeakPickerHiRes"];

/// Native ceiling on the spectra one annotation may produce: every match
/// copies a whole spectrum, and the number of matches is the number of
/// spectra times the number of transitions in the worst case. Checked by a
/// counting pass before anything is copied.
pub const MAX_ANNOTATIONS: usize = 1_000_000;

/// The source's `annotateSpectra` log line for a spectrum without precursors.
pub const NO_PRECURSOR_WARNING: &str =
    "annotateSpectra(): No precursor MZ found. Setting spectrum_mz to 0.";

/// A match against a spectral library: the source's `Match`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Match {
    /// The matched library spectrum.
    pub spectrum: MSSpectrum,
    /// The match score; higher is better.
    pub score: f64,
}

impl Match {
    /// A match of `spectrum` with `score`.
    pub fn new(spectrum: MSSpectrum, score: f64) -> Self {
        Self { spectrum, score }
    }
}

/// A spectral-library scoring strategy consumed by
/// [`TargetedSpectraExtractor::match_spectrum`]: the source's abstract
/// `Comparator`.
///
/// Implementations cache a reference library in [`init`](Comparator::init)
/// and then score query spectra against it in
/// [`generate_scores`](Comparator::generate_scores). The scoring direction is
/// "higher is better": the caller sorts the returned pairs by descending score
/// and keeps the top entries. Implementations interpret `options` as they see
/// fit; the caller passes user-supplied options through verbatim.
pub trait Comparator {
    /// Score `spectrum` against every entry of the cached library, keeping
    /// only matches with a score of at least `min_score` (inclusive), as
    /// unsorted `(library_index, score)` pairs.
    ///
    /// # Errors
    ///
    /// When the spectrum cannot be scored.
    fn generate_scores(&self, spectrum: &MSSpectrum, min_score: f64) -> Result<Vec<(usize, f64)>>;

    /// Load the reference library and (re)configure the comparator with
    /// implementation-specific `options`.
    ///
    /// # Errors
    ///
    /// When an option has the wrong type or the library cannot be prepared.
    fn init(&mut self, library: &[MSSpectrum], options: &BTreeMap<String, MetaValue>)
    -> Result<()>;

    /// The reference library passed to the most recent [`init`](Comparator::init).
    fn library(&self) -> &[MSSpectrum];
}

/// A [`Comparator`] that bins spectra and scores them with the binned
/// spectral contrast angle: the source's `BinnedSpectrumComparator`.
///
/// The library is binned once in [`init`](Comparator::init) and the cached
/// bins are reused; only the query spectrum is binned on each
/// [`generate_scores`](Comparator::generate_scores) call.
#[derive(Clone, Debug)]
pub struct BinnedSpectrumComparator {
    comparator: BinnedSpectralContrastAngle,
    library: Vec<MSSpectrum>,
    binned_library: Vec<BinnedSpectrum>,
    bin_size: f64,
    peak_spread: u32,
    bin_offset: f64,
}

impl BinnedSpectrumComparator {
    /// An empty comparator with the source's defaults: bin size `0.02`
    /// (suited to high-resolution data; the source notes `1.0` for nominal
    /// mass), peak spread `0` and bin offset `0.0` (`0.4` for nominal mass).
    ///
    /// # Errors
    ///
    /// Only if the contrast-angle functor cannot be constructed.
    pub fn new() -> Result<Self> {
        Ok(Self {
            comparator: BinnedSpectralContrastAngle::new()?,
            library: Vec::new(),
            binned_library: Vec::new(),
            bin_size: 0.02,
            peak_spread: 0,
            bin_offset: 0.0,
        })
    }

    /// The m/z bin width.
    pub fn bin_size(&self) -> f64 {
        self.bin_size
    }

    /// How many neighbouring bins receive a contribution from each peak.
    pub fn peak_spread(&self) -> u32 {
        self.peak_spread
    }

    /// The m/z offset of the bin boundaries.
    pub fn bin_offset(&self) -> f64 {
        self.bin_offset
    }

    fn config(bin_size: f64, peak_spread: u32, bin_offset: f64) -> BinConfig {
        // The source's BinnedSpectrum takes `float` size and offset.
        BinConfig {
            size: bin_size as f32,
            unit: BinUnit::Absolute,
            spread: peak_spread,
            offset: bin_offset as f32,
            ..BinConfig::default()
        }
    }
}

impl Comparator for BinnedSpectrumComparator {
    /// Bin `spectrum` with the configured bin parameters (without its
    /// precursor, as the source does) and score it against every cached
    /// library entry with the binned spectral contrast angle.
    ///
    /// # Errors
    ///
    /// As `BinnedSpectrum::new` (for example an unsorted spectrum, which the
    /// source only asserts in debug builds) and the contrast-angle score.
    fn generate_scores(&self, spectrum: &MSSpectrum, min_score: f64) -> Result<Vec<(usize, f64)>> {
        let mut query = spectrum.clone();
        query.precursors.clear();
        let binned = BinnedSpectrum::new(
            &query,
            Self::config(self.bin_size, self.peak_spread, self.bin_offset),
        )?;
        let mut scores = Vec::new();
        for (index, library) in self.binned_library.iter().enumerate() {
            let score = self.comparator.score(&binned, library)?;
            if score >= min_score {
                scores.push((index, score));
            }
        }
        Ok(scores)
    }

    /// Cache `library`, binning each spectrum with the current parameters,
    /// after applying the optional `options` overrides `"bin_size"` (float),
    /// `"peak_spread"` (non-negative integer) and `"bin_offset"` (float).
    ///
    /// The source logs the library size at info level; this port does not.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an option of the wrong type (the
    /// source's `DataValue` conversion throws `ConversionError`) and the errors
    /// of binning a library spectrum. On error the comparator is unchanged;
    /// the source has already applied the options it read.
    fn init(
        &mut self,
        library: &[MSSpectrum],
        options: &BTreeMap<String, MetaValue>,
    ) -> Result<()> {
        let bin_size = match options.get("bin_size") {
            Some(value) => numeric_meta(value, "bin_size")?,
            None => self.bin_size,
        };
        let peak_spread = match options.get("peak_spread") {
            Some(value) => u32::try_from(value.as_i64()?).map_err(|_| {
                Error::InvalidValue("option 'peak_spread' must be a non-negative integer".into())
            })?,
            None => self.peak_spread,
        };
        let bin_offset = match options.get("bin_offset") {
            Some(value) => numeric_meta(value, "bin_offset")?,
            None => self.bin_offset,
        };
        let config = Self::config(bin_size, peak_spread, bin_offset);
        let binned = library
            .iter()
            .map(|spectrum| BinnedSpectrum::new(spectrum, config))
            .collect::<Result<Vec<_>>>()?;
        self.bin_size = bin_size;
        self.peak_spread = peak_spread;
        self.bin_offset = bin_offset;
        self.library = library.to_vec();
        self.binned_library = binned;
        Ok(())
    }

    fn library(&self) -> &[MSSpectrum] {
        &self.library
    }
}

/// The `deisotoping:` parameters, used by
/// [`TargetedSpectraExtractor::deisotope_ms2_spectra`].
#[derive(Clone, Debug, PartialEq)]
pub struct DeisotopingParameters {
    /// `use_deisotoper`: deisotope before storing MSP (if no fragment
    /// annotation is used).
    pub use_deisotoper: bool,
    /// `fragment_tolerance`: tolerance used to match isotopic peaks.
    pub fragment_tolerance: f64,
    /// `fragment_unit`: `"ppm"` or `"Da"`, the unit of the fragment tolerance.
    pub fragment_unit: String,
    /// `min_charge`: the minimum charge considered (at least 1).
    pub min_charge: i32,
    /// `max_charge`: the maximum charge considered (at least 1).
    pub max_charge: i32,
    /// `min_isopeaks`: the minimum number of isotopic peaks (at least 2)
    /// required for an isotopic cluster.
    pub min_isopeaks: i32,
    /// `max_isopeaks`: the maximum number of isotopic peaks (at least 3)
    /// considered for an isotopic cluster.
    pub max_isopeaks: i32,
    /// `keep_only_deisotoped`: only monoisotopic peaks of fragments with an
    /// isotopic pattern are retained.
    pub keep_only_deisotoped: bool,
    /// `annotate_charge`: annotate the charge to the peaks.
    pub annotate_charge: bool,
}

/// Filters, annotates, picks and scores spectra based on a target list: the
/// source's `TargetedSpectraExtractor`.
///
/// Every field is a parameter of the source; see
/// [`TargetedSpectraExtractor::defaults`] for the full tree with its
/// descriptions.
#[derive(Clone, Debug)]
pub struct TargetedSpectraExtractor {
    /// `rt_window`, default `30.0`: the precursor retention-time window used
    /// during annotation, centred on each spectrum's retention time.
    pub rt_window: f64,
    /// `min_select_score`, default `0.7`, at least `0.0`: the minimum score of
    /// a spectrum [`select_spectra`](Self::select_spectra) keeps.
    pub min_select_score: f64,
    /// `mz_tolerance`, default `0.1`: the precursor m/z tolerance used during
    /// annotation.
    pub mz_tolerance: f64,
    /// `mz_unit_is_Da`, default `true`: the unit of
    /// [`mz_tolerance`](Self::mz_tolerance) and
    /// [`fwhm_threshold`](Self::fwhm_threshold), Da or (when `false`) ppm.
    pub mz_unit_is_da: bool,
    /// `use_gauss`, default `true`: smooth with a Gaussian filter, otherwise
    /// with a Savitzky-Golay filter.
    pub use_gauss: bool,
    /// `peak_height_min`, default `0.0`: a picked peak's intensity must be at
    /// least this.
    pub peak_height_min: f64,
    /// `peak_height_max`, default `f64::MAX`: a picked peak's intensity must
    /// be at most this.
    pub peak_height_max: f64,
    /// `fwhm_threshold`, default `0.0`: a picked peak's FWHM must be at least
    /// this.
    pub fwhm_threshold: f64,
    /// `tic_weight`, default `1.0`: the TIC's weight in a spectrum's score.
    pub tic_weight: f64,
    /// `fwhm_weight`, default `1.0`: the FWHM's weight in a spectrum's score.
    pub fwhm_weight: f64,
    /// `snr_weight`, default `1.0`: the S/N's weight in a spectrum's score.
    pub snr_weight: f64,
    /// `top_matches_to_report`, default `5`, at least `1`: the number of
    /// matches [`match_spectrum`](Self::match_spectrum) returns, best first.
    pub top_matches_to_report: usize,
    /// `min_match_score`, default `0.8`, in `[0, 1]`: the minimum score of a
    /// valid library match.
    pub min_match_score: f64,
    /// `min_fragment_mz`, default `0.0`: minimal m/z of a fragment ion chosen
    /// as a transition.
    pub min_fragment_mz: f64,
    /// `max_fragment_mz`, default `2000.0`: maximal m/z of a fragment ion
    /// chosen as a transition.
    pub max_fragment_mz: f64,
    /// `relative_allowable_product_mass`, default `10.0`: threshold m/z of a
    /// product relative to the precursor m/z (can be negative).
    pub relative_allowable_product_mass: f64,
    /// The `deisotoping:` subsection.
    pub deisotoping: DeisotopingParameters,
    /// `max_precursor_mass_threashold` (sic), default `10.0`: the ppm
    /// tolerance above the precursor m/z beyond which
    /// [`remove_ms2_spectra_peaks`](Self::remove_ms2_spectra_peaks) drops
    /// peaks.
    pub max_precursor_mass_threashold: f64,
    /// The `SavitzkyGolayFilter:`, `GaussFilter:`, `PeakPickerHiRes:` and
    /// `AccurateMassSearchEngine:` subtrees, keys kept with their prefix. The
    /// first three configure the smoothing filter and the picker on every
    /// [`pick_spectrum`](Self::pick_spectrum) call, as the source copies them
    /// from `param_` each time.
    pub subsection_parameters: Param,
}

fn entry(name: &str, value: ParamValue, description: &str, advanced: bool) -> ParamEntry {
    let mut entry = ParamEntry {
        name: name.into(),
        description: description.into(),
        value,
        ..ParamEntry::default()
    };
    if advanced {
        entry.tags.insert("advanced".into());
    }
    entry
}

fn float_entry(
    name: &str,
    value: f64,
    description: &str,
    min: Option<f64>,
    max: Option<f64>,
) -> ParamEntry {
    let mut e = entry(name, ParamValue::Float(value), description, false);
    if let Some(min) = min {
        e.min_float = min;
    }
    if let Some(max) = max {
        e.max_float = max;
    }
    e
}

fn int_entry(name: &str, value: i64, description: &str, min: Option<i32>) -> ParamEntry {
    let mut e = entry(name, ParamValue::Integer(value), description, false);
    if let Some(min) = min {
        e.min_int = min;
    }
    e
}

fn string_entry(
    name: &str,
    value: &str,
    description: &str,
    valid: &[&str],
    advanced: bool,
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

fn bool_str(value: bool) -> &'static str {
    if value { "true" } else { "false" }
}

fn node(name: &str, entries: Vec<ParamEntry>, nodes: Vec<ParamNode>) -> ParamNode {
    ParamNode {
        name: name.into(),
        description: String::new(),
        entries,
        nodes,
    }
}

/// The `AccurateMassSearchEngine` constructor's defaults, with
/// `keep_unidentified_masses` overridden to `"false"` as the
/// `TargetedSpectraExtractor` constructor does (value replaced, description
/// and valid strings kept, as `Param::setValue` on an existing entry does).
fn accurate_mass_search_engine_defaults() -> ParamNode {
    node(
        "AccurateMassSearchEngine",
        vec![
            float_entry(
                "mass_error_value",
                5.0,
                "Tolerance allowed for accurate mass search.",
                None,
                None,
            ),
            string_entry(
                "mass_error_unit",
                "ppm",
                "Unit of mass error (ppm or Da)",
                &["ppm", "Da"],
                false,
            ),
            string_entry(
                "ionization_mode",
                "positive",
                "Positive or negative ionization mode? If 'auto' is used, the first feature of the input map must contain the meta-value 'scan_polarity'. If its missing, the tool will exit with error.",
                &["positive", "negative", "auto"],
                false,
            ),
            string_entry(
                "isotopic_similarity",
                "false",
                "Computes a similarity score for each hit (only if the feature exhibits at least two isotopic mass traces).",
                &["false", "true"],
                false,
            ),
            string_entry(
                "positive_adducts",
                "CHEMISTRY/PositiveAdducts.tsv",
                "This file contains the list of potential positive adducts that will be looked for in the database. Edit the list if you wish to exclude/include adducts. By default CHEMISTRY/PositiveAdducts.tsv in OpenMS/share is used.",
                &[],
                true,
            ),
            string_entry(
                "negative_adducts",
                "CHEMISTRY/NegativeAdducts.tsv",
                "This file contains the list of potential negative adducts that will be looked for in the database. Edit the list if you wish to exclude/include adducts. By default CHEMISTRY/NegativeAdducts.tsv in OpenMS/share is used.",
                &[],
                true,
            ),
            string_entry(
                "use_feature_adducts",
                "false",
                "Whether to filter AMS candidates mismatching available feature adduct annotation.",
                &["false", "true"],
                false,
            ),
            string_entry(
                "keep_unidentified_masses",
                "false",
                "Keep features that did not yield any DB hit.",
                &["true", "false"],
                false,
            ),
            string_entry(
                "id_format",
                "legacy",
                "Use legacy (ProteinID/PeptideID based storage of metabolomics data) with mzTab-v1.0.0 as output format or novel Identification Data (ID) with mzTab-v2.0.0-M as output format (ID and its MzTab-M output is currently only support for featureXML files).",
                &["legacy", "ID"],
                false,
            ),
        ],
        vec![
            node(
                "db",
                vec![
                    entry(
                        "mapping",
                        ParamValue::StringList(vec!["CHEMISTRY/HMDBMappingFile.tsv".into()]),
                        "Database input file(s), containing three tab-separated columns of mass, formula, identifier. If 'mass' is 0, it is re-computed from the molecular sum formula. By default CHEMISTRY/HMDBMappingFile.tsv in OpenMS/share is used! If empty, the default will be used.",
                        false,
                    ),
                    entry(
                        "struct",
                        ParamValue::StringList(vec!["CHEMISTRY/HMDB2StructMapping.tsv".into()]),
                        "Database input file(s), containing four tab-separated columns of identifier, name, SMILES, INCHI.The identifier should match with mapping file. SMILES and INCHI are reported in the output, but not used otherwise. By default CHEMISTRY/HMDB2StructMapping.tsv in OpenMS/share is used! If empty, the default will be used.",
                        false,
                    ),
                ],
                Vec::new(),
            ),
            node(
                "mzTab",
                vec![string_entry(
                    "exportIsotopeIntensities",
                    "false",
                    "[featureXML input only] Export column with available isotope trace intensities (opt_global_MTint)",
                    &["false", "true"],
                    false,
                )],
                Vec::new(),
            ),
        ],
    )
}

/// The subsection defaults the source's constructor adds after
/// `getDefaultParameters`: Savitzky-Golay frame length 15 and order 3,
/// Gaussian width 0.2, picker signal-to-noise 1.0, and the
/// `AccurateMassSearchEngine` defaults. The first four carry no description,
/// as in the source.
fn subsection_defaults() -> Vec<ParamNode> {
    vec![
        node(
            "SavitzkyGolayFilter",
            vec![
                entry("frame_length", ParamValue::Integer(15), "", false),
                entry("polynomial_order", ParamValue::Integer(3), "", false),
            ],
            Vec::new(),
        ),
        node(
            "GaussFilter",
            vec![entry("gaussian_width", ParamValue::Float(0.2), "", false)],
            Vec::new(),
        ),
        node(
            "PeakPickerHiRes",
            vec![entry("signal_to_noise", ParamValue::Float(1.0), "", false)],
            Vec::new(),
        ),
        accurate_mass_search_engine_defaults(),
    ]
}

fn float_param(param: &Param, key: &str) -> Result<f64> {
    param.value(key)?.to_f64()
}

fn int_param(param: &Param, key: &str) -> Result<i32> {
    param.value(key)?.to_i32()
}

fn bool_param(param: &Param, key: &str) -> Result<bool> {
    param.value(key)?.to_bool()
}

/// A finite or non-finite double as a meta value, the way the source stores
/// any `double` in a `DataValue`.
fn float_meta(value: f64) -> MetaValue {
    MetaValue::try_from(value).unwrap_or_else(|_| MetaValue::source_float(value))
}

/// A numeric meta value as `double`, as the source's `DataValue` conversion
/// does for integer and floating values.
fn numeric_meta(value: &MetaValue, key: &str) -> Result<f64> {
    match value.data() {
        MetaValueData::Float(v) => Ok(*v),
        MetaValueData::Integer(v) => Ok(*v as f64),
        _ => Err(Error::InvalidValue(format!(
            "meta value '{key}' cannot be converted to a number"
        ))),
    }
}

/// `StringUtils::toStr(DataValue)`: a missing or empty value is `""`, a string
/// itself, anything else its display form.
fn lenient_string(metadata: &MetaInfo, key: &str) -> String {
    match metadata.get(key).map(MetaValue::data) {
        None | Some(MetaValueData::Empty) => String::new(),
        Some(MetaValueData::String(s)) => s.clone(),
        Some(_) => metadata[key].to_string(),
    }
}

/// `DataValue::operator std::string`: only a string value converts; anything
/// else, including a missing value, is the source's `ConversionError`.
fn strict_string(metadata: &MetaInfo, key: &str) -> Result<String> {
    match metadata.get(key).map(MetaValue::data) {
        Some(MetaValueData::String(s)) => Ok(s.clone()),
        _ => Err(Error::InvalidValue(format!(
            "meta value '{key}' is missing or not a string"
        ))),
    }
}

/// A feature at `(rt, mz)` with `intensity` and every other field at its
/// default, as a default-constructed source `Feature` after `setRT`, `setMZ`
/// and `setIntensity`.
fn new_feature(rt: f64, mz: f64, intensity: f32) -> Feature {
    Feature::from(BaseFeature::new(rt, mz, intensity))
}

fn first_precursor_mz(spectrum: &MSSpectrum) -> f64 {
    spectrum
        .precursors
        .first()
        .map_or(0.0, |precursor| precursor.mz)
}

/// `std::to_string(double)`: `%f`, six decimals.
fn to_string_fixed(value: f64) -> String {
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

impl TargetedSpectraExtractor {
    /// The source's default parameters: the entries of `getDefaultParameters`
    /// followed by the subsections its constructor adds.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure, which the constant tree cannot
    /// trigger.
    pub fn defaults() -> Result<Param> {
        let mut root = Self::default_parameters_root();
        root.nodes.extend(subsection_defaults());
        Param::from_root(root)
    }

    /// The source's public `getDefaultParameters`: the class's own entries,
    /// without the subsections the constructor adds.
    ///
    /// # Errors
    ///
    /// As [`defaults`](Self::defaults).
    pub fn default_parameters() -> Result<Param> {
        Param::from_root(Self::default_parameters_root())
    }

    fn default_parameters_root() -> ParamNode {
        let entries = vec![
            float_entry(
                "rt_window",
                30.0,
                "Precursor Retention Time window used during the annotation phase.\nFor each transition in the target list, annotateSpectra() looks for the first spectrum whose RT time falls within the RT Window, whose left and right limits are computed at each analyzed spectrum.\nAlso the spectrum's percursor MZ is checked against the transition MZ.",
                None,
                None,
            ),
            float_entry(
                "min_select_score",
                0.7,
                "Used in selectSpectra(), after the spectra have been assigned a score.\nRemained transitions will have at least one spectrum assigned.\nEach spectrum needs to have a score >= min_select_score_ to be valid, otherwise it gets filtered out.",
                Some(0.0),
                None,
            ),
            float_entry(
                "mz_tolerance",
                0.1,
                "Precursor MZ tolerance used during the annotation phase.\nFor each transition in the target list, annotateSpectra() looks for the first spectrum whose precursor MZ is close enough (+-mz_tolerance_) to the transition's MZ.\nAlso the spectrum's precursor RT is checked against the transition RT.",
                None,
                None,
            ),
            string_entry(
                "mz_unit_is_Da",
                "true",
                "Unit to use for mz_tolerance_ and fwhm_threshold_: true for Da, false for ppm.",
                &["false", "true"],
                false,
            ),
            string_entry(
                "use_gauss",
                "true",
                "Use Gaussian filter for smoothing (alternative is Savitzky-Golay filter)",
                &["false", "true"],
                false,
            ),
            float_entry(
                "peak_height_min",
                0.0,
                "Used in pickSpectrum(), a peak's intensity needs to be >= peak_height_min_ for it to be picked.",
                Some(0.0),
                None,
            ),
            float_entry(
                "peak_height_max",
                f64::MAX,
                "Used in pickSpectrum(), a peak's intensity needs to be <= peak_height_max_ for it to be picked.",
                Some(0.0),
                None,
            ),
            float_entry(
                "fwhm_threshold",
                0.0,
                "Used in pickSpectrum(), a peak's FWHM needs to be >= fwhm_threshold_ for it to be picked.",
                Some(0.0),
                None,
            ),
            float_entry(
                "tic_weight",
                1.0,
                "TIC weight when scoring spectra.",
                Some(0.0),
                None,
            ),
            float_entry(
                "fwhm_weight",
                1.0,
                "FWHM weight when scoring spectra.",
                Some(0.0),
                None,
            ),
            float_entry(
                "snr_weight",
                1.0,
                "SNR weight when scoring spectra.",
                Some(0.0),
                None,
            ),
            int_entry(
                "top_matches_to_report",
                5,
                "The number of matches to output from `matchSpectrum()`. These will be the matches of highest scores, sorted in descending order.",
                Some(1),
            ),
            float_entry(
                "min_match_score",
                0.8,
                "Minimum score for a match to be considered valid in `matchSpectrum()`.",
                Some(0.0),
                Some(1.0),
            ),
            float_entry(
                "min_fragment_mz",
                0.0,
                "Minimal m/z of a fragment ion choosen as a transition",
                None,
                None,
            ),
            float_entry(
                "max_fragment_mz",
                2000.0,
                "Maximal m/z of a fragment ion choosen as a transition",
                None,
                None,
            ),
            float_entry(
                "relative_allowable_product_mass",
                10.0,
                "Threshold m/z of a product relatively to the precurosor m/z (can be negative)",
                None,
                None,
            ),
            float_entry(
                "max_precursor_mass_threashold",
                10.0,
                "Tolerance used to set intensity to zero for peaks with mz higher than precursor mz",
                None,
                None,
            ),
        ];
        let deisotoping = node(
            "deisotoping",
            vec![
                entry(
                    "use_deisotoper",
                    ParamValue::String("false".into()),
                    "Use Deisotoper (if no fragment annotation is used)",
                    false,
                ),
                float_entry(
                    "fragment_tolerance",
                    1.0,
                    "Tolerance used to match isotopic peaks",
                    None,
                    None,
                ),
                string_entry(
                    "fragment_unit",
                    "ppm",
                    "Unit of the fragment tolerance",
                    &["ppm", "Da"],
                    false,
                ),
                int_entry("min_charge", 1, "The minimum charge considered", Some(1)),
                int_entry("max_charge", 1, "The maximum charge considered", Some(1)),
                int_entry(
                    "min_isopeaks",
                    2,
                    "The minimum number of isotopic peaks (at least 2) required for an isotopic cluster",
                    Some(2),
                ),
                int_entry(
                    "max_isopeaks",
                    3,
                    "The maximum number of isotopic peaks (at least 2) considered for an isotopic cluster",
                    Some(3),
                ),
                entry(
                    "keep_only_deisotoped",
                    ParamValue::String("false".into()),
                    "Only monoisotopic peaks of fragments with isotopic pattern are retained",
                    false,
                ),
                entry(
                    "annotate_charge",
                    ParamValue::String("false".into()),
                    "Annotate the charge to the peaks",
                    false,
                ),
            ],
            Vec::new(),
        );
        // The source adds `max_precursor_mass_threashold` after the
        // `deisotoping:` entries; entries and nodes are stored apart, so the
        // entry order above is the source's order of top-level entries.
        node("", entries, vec![deisotoping])
    }

    /// An extractor with the source's default parameters.
    ///
    /// # Errors
    ///
    /// As [`defaults`](Self::defaults).
    pub fn new() -> Result<Self> {
        Self::from_param(&Self::defaults()?)
    }

    /// The current parameters in the layout of [`defaults`](Self::defaults)
    /// (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a value has no source parameter
    /// representation (a count beyond the source's 32-bit integers).
    pub fn to_param(&self) -> Result<Param> {
        let mut param = Self::default_parameters()?;
        let set = |param: &mut Param, key: &str, value: ParamValue| -> Result<()> {
            let description = param.description(key)?.to_owned();
            let tags: Vec<String> = param.tags(key)?.iter().cloned().collect();
            param.set_value(key, value, &description, &tags)
        };
        set(&mut param, "rt_window", ParamValue::Float(self.rt_window))?;
        set(
            &mut param,
            "min_select_score",
            ParamValue::Float(self.min_select_score),
        )?;
        set(
            &mut param,
            "mz_tolerance",
            ParamValue::Float(self.mz_tolerance),
        )?;
        set(
            &mut param,
            "mz_unit_is_Da",
            ParamValue::String(bool_str(self.mz_unit_is_da).into()),
        )?;
        set(
            &mut param,
            "use_gauss",
            ParamValue::String(bool_str(self.use_gauss).into()),
        )?;
        set(
            &mut param,
            "peak_height_min",
            ParamValue::Float(self.peak_height_min),
        )?;
        set(
            &mut param,
            "peak_height_max",
            ParamValue::Float(self.peak_height_max),
        )?;
        set(
            &mut param,
            "fwhm_threshold",
            ParamValue::Float(self.fwhm_threshold),
        )?;
        set(&mut param, "tic_weight", ParamValue::Float(self.tic_weight))?;
        set(
            &mut param,
            "fwhm_weight",
            ParamValue::Float(self.fwhm_weight),
        )?;
        set(&mut param, "snr_weight", ParamValue::Float(self.snr_weight))?;
        let top = i32::try_from(self.top_matches_to_report).map_err(|_| {
            Error::InvalidValue(
                "'top_matches_to_report' does not fit the source 32-bit parameter".into(),
            )
        })?;
        set(
            &mut param,
            "top_matches_to_report",
            ParamValue::Integer(i64::from(top)),
        )?;
        set(
            &mut param,
            "min_match_score",
            ParamValue::Float(self.min_match_score),
        )?;
        set(
            &mut param,
            "min_fragment_mz",
            ParamValue::Float(self.min_fragment_mz),
        )?;
        set(
            &mut param,
            "max_fragment_mz",
            ParamValue::Float(self.max_fragment_mz),
        )?;
        set(
            &mut param,
            "relative_allowable_product_mass",
            ParamValue::Float(self.relative_allowable_product_mass),
        )?;
        let d = &self.deisotoping;
        set(
            &mut param,
            "deisotoping:use_deisotoper",
            ParamValue::String(bool_str(d.use_deisotoper).into()),
        )?;
        set(
            &mut param,
            "deisotoping:fragment_tolerance",
            ParamValue::Float(d.fragment_tolerance),
        )?;
        set(
            &mut param,
            "deisotoping:fragment_unit",
            ParamValue::String(d.fragment_unit.clone()),
        )?;
        set(
            &mut param,
            "deisotoping:min_charge",
            ParamValue::Integer(i64::from(d.min_charge)),
        )?;
        set(
            &mut param,
            "deisotoping:max_charge",
            ParamValue::Integer(i64::from(d.max_charge)),
        )?;
        set(
            &mut param,
            "deisotoping:min_isopeaks",
            ParamValue::Integer(i64::from(d.min_isopeaks)),
        )?;
        set(
            &mut param,
            "deisotoping:max_isopeaks",
            ParamValue::Integer(i64::from(d.max_isopeaks)),
        )?;
        set(
            &mut param,
            "deisotoping:keep_only_deisotoped",
            ParamValue::String(bool_str(d.keep_only_deisotoped).into()),
        )?;
        set(
            &mut param,
            "deisotoping:annotate_charge",
            ParamValue::String(bool_str(d.annotate_charge).into()),
        )?;
        set(
            &mut param,
            "max_precursor_mass_threashold",
            ParamValue::Float(self.max_precursor_mass_threashold),
        )?;
        param.insert("", &self.subsection_parameters)?;
        Ok(param)
    }

    /// Build an extractor from parameters, as source `setParameters` followed
    /// by `updateMembers_`.
    ///
    /// Missing parameters take their defaults. Keys under the
    /// `SavitzkyGolayFilter:`, `GaussFilter:` and `PeakPickerHiRes:`
    /// subsections are not checked here (the source registers them as
    /// subsections); they are read when [`pick_spectrum`](Self::pick_spectrum)
    /// configures the filter and the picker. Unknown keys elsewhere are
    /// ignored, as the source only warns; see
    /// [`from_param_with_warnings`](Self::from_param_with_warnings).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong value type
    /// or violates its restriction (the source's `Exception::InvalidParameter`).
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
        let mut handler = DefaultParamHandler::new(TARGETED_SPECTRA_EXTRACTOR_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        let subsections: Vec<String> = SUBSECTIONS.iter().map(|s| (*s).to_owned()).collect();
        handler.set_subsections(&subsections)?;
        handler.set_parameters_with(parameters, Self::from_complete_param)
    }

    fn from_complete_param(param: &Param) -> Result<Self> {
        let top = int_param(param, "top_matches_to_report")?;
        let mut subsection_parameters = Param::new();
        for prefix in [
            "SavitzkyGolayFilter:",
            "GaussFilter:",
            "PeakPickerHiRes:",
            "AccurateMassSearchEngine:",
        ] {
            subsection_parameters.insert("", &param.copy(prefix, false)?)?;
        }
        Ok(Self {
            rt_window: float_param(param, "rt_window")?,
            min_select_score: float_param(param, "min_select_score")?,
            mz_tolerance: float_param(param, "mz_tolerance")?,
            mz_unit_is_da: bool_param(param, "mz_unit_is_Da")?,
            use_gauss: bool_param(param, "use_gauss")?,
            peak_height_min: float_param(param, "peak_height_min")?,
            peak_height_max: float_param(param, "peak_height_max")?,
            fwhm_threshold: float_param(param, "fwhm_threshold")?,
            tic_weight: float_param(param, "tic_weight")?,
            fwhm_weight: float_param(param, "fwhm_weight")?,
            snr_weight: float_param(param, "snr_weight")?,
            top_matches_to_report: usize::try_from(top).map_err(|_| {
                Error::InvalidValue("'top_matches_to_report' must not be negative".into())
            })?,
            min_match_score: float_param(param, "min_match_score")?,
            min_fragment_mz: float_param(param, "min_fragment_mz")?,
            max_fragment_mz: float_param(param, "max_fragment_mz")?,
            relative_allowable_product_mass: float_param(param, "relative_allowable_product_mass")?,
            deisotoping: DeisotopingParameters {
                use_deisotoper: bool_param(param, "deisotoping:use_deisotoper")?,
                fragment_tolerance: float_param(param, "deisotoping:fragment_tolerance")?,
                fragment_unit: param
                    .value("deisotoping:fragment_unit")?
                    .as_str()?
                    .to_owned(),
                min_charge: int_param(param, "deisotoping:min_charge")?,
                max_charge: int_param(param, "deisotoping:max_charge")?,
                min_isopeaks: int_param(param, "deisotoping:min_isopeaks")?,
                max_isopeaks: int_param(param, "deisotoping:max_isopeaks")?,
                keep_only_deisotoped: bool_param(param, "deisotoping:keep_only_deisotoped")?,
                annotate_charge: bool_param(param, "deisotoping:annotate_charge")?,
            },
            max_precursor_mass_threashold: float_param(param, "max_precursor_mass_threashold")?,
            subsection_parameters,
        })
    }

    /// The m/z tolerance annotation uses: `mz_tolerance` in Da, or
    /// `mz_tolerance / 1e6` when the unit is ppm.
    ///
    /// The source divides a ppm tolerance by `1e6` without multiplying by an
    /// m/z, so a ppm tolerance acts as a tiny absolute one; reproduced, and
    /// listed in the support document as a defect candidate.
    fn annotation_mz_tolerance(&self) -> f64 {
        if self.mz_unit_is_da {
            self.mz_tolerance
        } else {
            self.mz_tolerance / 1e6
        }
    }

    /// The `(retention time in seconds, precursor m/z, peptide reference)` of
    /// every transition, as `annotateSpectra` reads them.
    fn transition_targets(targeted_exp: &TargetedExperiment) -> Result<Vec<(f64, f64, String)>> {
        targeted_exp
            .transitions()
            .iter()
            .map(|transition| {
                let peptide: &Peptide = targeted_exp.peptide_by_ref(transition.peptide_ref())?;
                let mut target_rt = peptide.retention_time()?;
                if peptide.retention_time_unit()? == RTUnit::Minute {
                    target_rt *= 60.0;
                }
                Ok((
                    target_rt,
                    transition.precursor_mz(),
                    transition.peptide_ref().to_owned(),
                ))
            })
            .collect()
    }

    /// For every spectrum, the targets it matches, in transition order.
    fn annotation_matches(
        &self,
        spectra: &[MSSpectrum],
        targeted_exp: &TargetedExperiment,
    ) -> Result<Vec<Vec<usize>>> {
        if spectra.is_empty() {
            // The source reads the peptides only inside the spectrum loop.
            return Ok(Vec::new());
        }
        let targets = Self::transition_targets(targeted_exp)?;
        let mz_tolerance = self.annotation_mz_tolerance();
        let mut total = 0_usize;
        let mut matches = Vec::with_capacity(spectra.len());
        for spectrum in spectra {
            let rt_left = spectrum.rt - self.rt_window / 2.0;
            let rt_right = spectrum.rt + self.rt_window / 2.0;
            let spectrum_mz = first_precursor_mz(spectrum);
            // A spectrum m/z of 0 disables the m/z check; the source's lower
            // limit is then `numeric_limits<double>::min()`, the smallest
            // positive normal double, so a target m/z of 0 or below still fails.
            let (mz_left, mz_right) = if spectrum_mz != 0.0 {
                (spectrum_mz - mz_tolerance, spectrum_mz + mz_tolerance)
            } else {
                (f64::MIN_POSITIVE, f64::MAX)
            };
            let hits: Vec<usize> = targets
                .iter()
                .enumerate()
                .filter(|(_, (rt, mz, _))| {
                    *rt >= rt_left && *rt <= rt_right && *mz >= mz_left && *mz <= mz_right
                })
                .map(|(index, _)| index)
                .collect();
            total = total.saturating_add(hits.len());
            if total > MAX_ANNOTATIONS {
                return Err(Error::InvalidValue(format!(
                    "annotation would copy more than {MAX_ANNOTATIONS} spectra"
                )));
            }
            matches.push(hits);
        }
        Ok(matches)
    }

    /// Filter and annotate the spectra that could match the transitions of the
    /// target list, with one feature per annotated spectrum.
    ///
    /// A spectrum matches a transition when the transition's peptide
    /// retention time (converted from minutes when its unit says so) lies
    /// within `rt_window / 2` of the spectrum's retention time and the
    /// transition's precursor m/z within `mz_tolerance` of the spectrum's first
    /// precursor m/z. A spectrum without precursors (the source logs
    /// [`NO_PRECURSOR_WARNING`]) has m/z 0, which disables the m/z check. Every
    /// match adds a copy of the spectrum named with the transition's peptide
    /// reference, in spectrum order and then transition order, and a feature
    /// at the spectrum's retention time and precursor m/z with the meta value
    /// `transition_name`.
    ///
    /// The source's `compute_features = false` overload is
    /// [`annotate_spectra_only`](Self::annotate_spectra_only).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a transition references a peptide
    /// that is missing or has no retention time (the source throws there too),
    /// and when more than [`MAX_ANNOTATIONS`] spectra would be produced.
    pub fn annotate_spectra(
        &self,
        spectra: &[MSSpectrum],
        targeted_exp: &TargetedExperiment,
    ) -> Result<(Vec<MSSpectrum>, FeatureMap)> {
        self.annotate(spectra, targeted_exp, true)
    }

    /// [`annotate_spectra`](Self::annotate_spectra) without features: the
    /// source's three-argument overload.
    ///
    /// # Errors
    ///
    /// As [`annotate_spectra`](Self::annotate_spectra).
    pub fn annotate_spectra_only(
        &self,
        spectra: &[MSSpectrum],
        targeted_exp: &TargetedExperiment,
    ) -> Result<Vec<MSSpectrum>> {
        Ok(self.annotate(spectra, targeted_exp, false)?.0)
    }

    fn annotate(
        &self,
        spectra: &[MSSpectrum],
        targeted_exp: &TargetedExperiment,
        compute_features: bool,
    ) -> Result<(Vec<MSSpectrum>, FeatureMap)> {
        let matches = self.annotation_matches(spectra, targeted_exp)?;
        let transitions = targeted_exp.transitions();
        let mut annotated = Vec::new();
        let mut features = FeatureMap::new();
        for (spectrum, hits) in spectra.iter().zip(&matches) {
            let spectrum_mz = first_precursor_mz(spectrum);
            for &index in hits {
                let name = transitions[index].peptide_ref();
                let mut copy = spectrum.clone();
                copy.name = name.to_owned();
                annotated.push(copy);
                if compute_features {
                    let mut feature = new_feature(spectrum.rt, spectrum_mz, 0.0);
                    feature
                        .metadata
                        .insert("transition_name".into(), MetaValue::from(name));
                    features.features.push(feature);
                }
            }
        }
        Ok((annotated, features))
    }

    /// Annotate MS2 spectra with the MS1 features (for example from
    /// `AccurateMassSearch`) they were likely derived from.
    ///
    /// MS1 spectra are skipped. Every MS1 feature, or each of its
    /// subordinates when it has any, whose `PeptideRef` is not `"null"` and
    /// whose retention time lies within `rt_window / 2` and m/z within
    /// `mz_tolerance / 2` of the spectrum (the m/z check is skipped for a
    /// spectrum without precursors), adds a copy of the spectrum named with the
    /// `PeptideRef`, and an MS2 feature with a fresh unique id from `ids`, the
    /// spectrum's retention time and precursor m/z, the MS1 feature's
    /// intensity and the meta values `native_id` and `PeptideRef`.
    ///
    /// Unlike the targeted overload, this one uses `mz_tolerance` as a full
    /// window in the unit of m/z whatever `mz_unit_is_Da` says, as the source
    /// does. The source appends to the caller's MS2 feature map; this returns
    /// the new features.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when more than [`MAX_ANNOTATIONS`]
    /// spectra would be produced.
    pub fn annotate_spectra_with_ms1_features(
        &self,
        spectra: &[MSSpectrum],
        ms1_features: &FeatureMap,
        ids: &mut UniqueIdGenerator,
    ) -> Result<(Vec<MSSpectrum>, FeatureMap)> {
        let check = |spectrum_mz: f64, spectrum_rt: f64, target_mz: f64, target_rt: f64| {
            let rt_left = spectrum_rt - self.rt_window / 2.0;
            let rt_right = spectrum_rt + self.rt_window / 2.0;
            let mz_left = spectrum_mz - self.mz_tolerance / 2.0;
            let mz_right = spectrum_mz + self.mz_tolerance / 2.0;
            if spectrum_mz != 0.0 {
                target_rt >= rt_left
                    && target_rt <= rt_right
                    && target_mz >= mz_left
                    && target_mz <= mz_right
            } else {
                target_rt >= rt_left && target_rt <= rt_right
            }
        };
        // Counting pass: the preflight of the copies below.
        let mut total = 0_usize;
        for spectrum in spectra.iter().filter(|s| s.ms_level != 1) {
            let spectrum_mz = first_precursor_mz(spectrum);
            for feature in &ms1_features.features {
                let candidates: Vec<&Feature> = if feature.subordinates.is_empty() {
                    vec![feature]
                } else {
                    feature.subordinates.iter().collect()
                };
                for candidate in candidates {
                    if lenient_string(&candidate.metadata, "PeptideRef") != "null"
                        && check(spectrum_mz, spectrum.rt, candidate.mz, candidate.rt)
                    {
                        total = total.saturating_add(1);
                    }
                }
            }
            if total > MAX_ANNOTATIONS {
                return Err(Error::InvalidValue(format!(
                    "annotation would copy more than {MAX_ANNOTATIONS} spectra"
                )));
            }
        }
        let mut annotated = Vec::with_capacity(total);
        let mut ms2_features = FeatureMap::new();
        for spectrum in spectra.iter().filter(|s| s.ms_level != 1) {
            let spectrum_mz = first_precursor_mz(spectrum);
            for feature in &ms1_features.features {
                let candidates: Vec<&Feature> = if feature.subordinates.is_empty() {
                    vec![feature]
                } else {
                    feature.subordinates.iter().collect()
                };
                for candidate in candidates {
                    let peptide_ref = lenient_string(&candidate.metadata, "PeptideRef");
                    let native_id = lenient_string(&candidate.metadata, "native_id");
                    if peptide_ref == "null"
                        || !check(spectrum_mz, spectrum.rt, candidate.mz, candidate.rt)
                    {
                        continue;
                    }
                    let mut copy = spectrum.clone();
                    copy.name = peptide_ref.clone();
                    annotated.push(copy);
                    let mut ms2 = new_feature(spectrum.rt, spectrum_mz, candidate.intensity);
                    ms2.unique_id = ids.get_unique_id();
                    ms2.metadata
                        .insert("native_id".into(), MetaValue::from(native_id));
                    ms2.metadata
                        .insert("PeptideRef".into(), MetaValue::from(peptide_ref));
                    ms2_features.features.push(ms2);
                }
            }
        }
        Ok((annotated, ms2_features))
    }

    /// The smoothing filter the parameters select, configured from the
    /// `GaussFilter:` or `SavitzkyGolayFilter:` subsection over that filter's
    /// own defaults.
    fn smoothing_filter(&self) -> Result<Box<dyn SpectrumFilter>> {
        if self.use_gauss {
            let section = self.subsection_parameters.copy("GaussFilter:", true)?;
            let value = |key: &str, default: ParamValue| -> Result<ParamValue> {
                if section.exists(key)? {
                    Ok(section.value(key)?.clone())
                } else {
                    Ok(default)
                }
            };
            let width = value("gaussian_width", ParamValue::Float(0.2))?.to_f64()?;
            let ppm = value("ppm_tolerance", ParamValue::Float(10.0))?.to_f64()?;
            let use_ppm =
                value("use_ppm_tolerance", ParamValue::String("false".into()))?.to_bool()?;
            let width = if use_ppm {
                GaussianWidth::Ppm(ppm)
            } else {
                GaussianWidth::Absolute(width)
            };
            Ok(Box::new(GaussFilter::new(width)?))
        } else {
            let section = self
                .subsection_parameters
                .copy("SavitzkyGolayFilter:", true)?;
            let integer = |key: &str, default: i64| -> Result<usize> {
                let value = if section.exists(key)? {
                    section.value(key)?.to_i64()?
                } else {
                    default
                };
                usize::try_from(value).map_err(|_| {
                    Error::InvalidValue(format!("SavitzkyGolayFilter:{key} must not be negative"))
                })
            };
            Ok(Box::new(SavitzkyGolayFilter::new(
                integer("frame_length", 11)?,
                integer("polynomial_order", 4)?,
            )?))
        }
    }

    /// The picker `pickSpectrum` configures: `PeakPickerHiRes` defaults
    /// updated from the `PeakPickerHiRes:` subsection, then spacing constraints
    /// disabled (`spacing_difference` and `spacing_difference_gap` 0, "since
    /// we're dealing with spectrum") and absolute FWHM reporting enabled.
    fn picker(&self) -> Result<PeakPickerHiRes> {
        let mut picker = PeakPickerHiRes::from_param(
            &self.subsection_parameters.copy("PeakPickerHiRes:", true)?,
        )?;
        picker.spacing_difference = 0.0;
        picker.spacing_difference_gap = 0.0;
        picker.report_fwhm = Some(FwhmUnit::Absolute);
        Ok(picker)
    }

    /// Pick a spectrum's peaks.
    ///
    /// The spectrum is first smoothed with a Gaussian filter (the default) or,
    /// with `use_gauss` off, a Savitzky-Golay filter; the peaks are then picked
    /// with `PeakPickerHiRes`. Parameters under `GaussFilter:`,
    /// `SavitzkyGolayFilter:` and `PeakPickerHiRes:` are passed on. The picked
    /// peaks are kept when their intensity lies within
    /// `[peak_height_min, peak_height_max]` and their FWHM is at least
    /// `fwhm_threshold` (divided by `1e6` when `mz_unit_is_Da` is off, as the
    /// source does although the FWHM is absolute).
    ///
    /// The FWHM is read from the picked spectrum's first float data array, as
    /// in the source; for an input with an ion-mobility array that is the
    /// mobility array instead. When no peak is kept the result is an empty
    /// spectrum without any metadata (the source's `clear(true)`), so it may
    /// be empty: no peaks were found.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the spectrum is not sorted by
    /// position (the source's `Exception::IllegalArgument`, "Spectrum must be
    /// sorted by position"), and the errors of the smoothing filter and the
    /// picker, including invalid subsection parameters.
    pub fn pick_spectrum(&self, spectrum: &MSSpectrum) -> Result<MSSpectrum> {
        if !spectrum.is_sorted() {
            return Err(Error::InvalidValue(
                "Spectrum must be sorted by position".into(),
            ));
        }
        let mut smoothed = spectrum.clone();
        self.smoothing_filter()?.filter_spectrum(&mut smoothed)?;
        let mut picked = self.picker()?.pick_spectrum(&smoothed)?.spectrum;
        let fwhm_threshold = if self.mz_unit_is_da {
            self.fwhm_threshold
        } else {
            self.fwhm_threshold / 1e6
        };
        let fwhm = picked
            .float_data_arrays
            .first()
            .ok_or_else(|| Error::InvalidValue("the picked spectrum has no FWHM array".into()))?;
        let mut keep = Vec::with_capacity(picked.len());
        for (index, peak) in picked.peaks.iter().enumerate() {
            let width = *fwhm.data.get(index).ok_or_else(|| {
                Error::InvalidValue(
                    "the picked spectrum's FWHM array is shorter than its peaks".into(),
                )
            })?;
            let intensity = f64::from(peak.intensity);
            if intensity >= self.peak_height_min
                && intensity <= self.peak_height_max
                && f64::from(width) >= fwhm_threshold
            {
                keep.push(index);
            }
        }
        if keep.is_empty() {
            picked = MSSpectrum::default();
        } else if keep.len() != picked.len() {
            picked.select(&keep)?;
        }
        Ok(picked)
    }

    /// Score the annotated spectra and record the scores in their features.
    ///
    /// For each spectrum: the total ion current of the annotated spectrum, the
    /// average FWHM of its picked peaks (the picked spectrum's first float data
    /// array; `0` when it has no arrays) and the average signal-to-noise ratio
    /// over the annotated spectrum (`SignalToNoiseEstimatorMedian` with a
    /// 40-Th window, a noise of 2 for sparse windows and at least 10 points
    /// per window). The score is
    /// `log10(TIC) * tic_weight + (1 / avgFWHM) * fwhm_weight + avgSNR * snr_weight`.
    ///
    /// The scored spectra are copies of the annotated ones whose float data
    /// arrays are resized to five: arrays 1 to 4 are named `score`,
    /// `log10_total_tic`, `inverse_avgFWHM` and `avgSNR` and get the value
    /// appended. Such one-element arrays are shorter than the peak list, which
    /// the crate's spectrum validation refuses downstream; the source writes
    /// the same layout. Each feature gets the score as its intensity, the meta
    /// values `log10_total_tic`, `inverse_avgFWHM`, `avgFWHM` and `avgSNR`, and
    /// one subordinate per picked peak with its m/z, intensity and `FWHM` meta
    /// value. An empty picked or annotated spectrum gives an infinite or NaN
    /// score, as in the source.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] when `features` and `annotated_spectra`
    /// differ in size (the source's `Exception::InvalidSize`) or when
    /// `picked_spectra` differs from `annotated_spectra` in size (the source
    /// reads out of bounds), and the errors of the noise estimation. On error
    /// `features` is unchanged.
    pub fn score_spectra(
        &self,
        annotated_spectra: &[MSSpectrum],
        picked_spectra: &[MSSpectrum],
        features: &mut FeatureMap,
    ) -> Result<Vec<MSSpectrum>> {
        if annotated_spectra.len() != features.features.len() {
            return Err(Error::InvalidRange(format!(
                "scored_spectra size does not match features size ({} vs {})",
                annotated_spectra.len(),
                features.features.len()
            )));
        }
        let (scored, updates) = self.score(annotated_spectra, picked_spectra, true)?;
        for (feature, update) in features.features.iter_mut().zip(updates) {
            let Feature {
                base, subordinates, ..
            } = update;
            feature.intensity = base.intensity;
            feature.metadata.extend(base.metadata);
            feature.subordinates = subordinates;
        }
        Ok(scored)
    }

    /// [`score_spectra`](Self::score_spectra) without features: the source's
    /// three-argument overload.
    ///
    /// # Errors
    ///
    /// As [`score_spectra`](Self::score_spectra), without the feature count.
    pub fn score_spectra_only(
        &self,
        annotated_spectra: &[MSSpectrum],
        picked_spectra: &[MSSpectrum],
    ) -> Result<Vec<MSSpectrum>> {
        Ok(self.score(annotated_spectra, picked_spectra, false)?.0)
    }

    fn score(
        &self,
        annotated_spectra: &[MSSpectrum],
        picked_spectra: &[MSSpectrum],
        compute_features: bool,
    ) -> Result<(Vec<MSSpectrum>, Vec<Feature>)> {
        if picked_spectra.len() != annotated_spectra.len() {
            return Err(Error::InvalidRange(format!(
                "{} picked spectra for {} annotated spectra",
                picked_spectra.len(),
                annotated_spectra.len()
            )));
        }
        let estimator = SignalToNoiseEstimatorMedian {
            window_length: 40.0,
            noise_for_empty_window: 2.0,
            min_required_elements: 10,
            ..SignalToNoiseEstimatorMedian::default()
        };
        let compatibility = PickingCompatibility::default();
        let mut scored = Vec::with_capacity(annotated_spectra.len());
        let mut updates = Vec::new();
        for (annotated, picked) in annotated_spectra.iter().zip(picked_spectra) {
            let mut total_tic = 0.0;
            for peak in &annotated.peaks {
                total_tic += f64::from(peak.intensity);
            }
            let mut avg_fwhm = 0.0;
            if let Some(fwhm) = picked.float_data_arrays.first() {
                for width in &fwhm.data {
                    avg_fwhm += f64::from(*width);
                }
                avg_fwhm /= fwhm.data.len() as f64;
            }
            let mut avg_snr = 0.0;
            if !annotated.peaks.is_empty() {
                let estimates = estimator.estimate_spectrum(annotated, &compatibility)?;
                for snr in &estimates.signal_to_noise {
                    avg_snr += snr;
                }
            }
            avg_snr /= annotated.len() as f64;

            let log10_total_tic = total_tic.log10();
            let inverse_avg_fwhm = 1.0 / avg_fwhm;
            let score = log10_total_tic * self.tic_weight
                + inverse_avg_fwhm * self.fwhm_weight
                + avg_snr * self.snr_weight;

            let mut copy = annotated.clone();
            copy.float_data_arrays.resize_with(5, DataArray::default);
            for (index, (name, value)) in [
                ("score", score),
                ("log10_total_tic", log10_total_tic),
                ("inverse_avgFWHM", inverse_avg_fwhm),
                ("avgSNR", avg_snr),
            ]
            .into_iter()
            .enumerate()
            {
                let array = &mut copy.float_data_arrays[index + 1];
                array.name = name.to_owned();
                array.data.push(value as f32);
            }
            scored.push(copy);

            if compute_features {
                // The intensity of a feature is (proportional to) its total ion count.
                let mut update = new_feature(0.0, 0.0, score as f32);
                update
                    .metadata
                    .insert("log10_total_tic".into(), float_meta(log10_total_tic));
                update
                    .metadata
                    .insert("inverse_avgFWHM".into(), float_meta(inverse_avg_fwhm));
                update
                    .metadata
                    .insert("avgFWHM".into(), float_meta(avg_fwhm));
                update.metadata.insert("avgSNR".into(), float_meta(avg_snr));
                let widths = picked.float_data_arrays.first();
                let mut subordinates = Vec::with_capacity(picked.len());
                for (index, peak) in picked.peaks.iter().enumerate() {
                    let width = widths.and_then(|w| w.data.get(index)).ok_or_else(|| {
                        Error::InvalidValue(
                            "a picked spectrum's FWHM array is shorter than its peaks".into(),
                        )
                    })?;
                    let mut subordinate = new_feature(0.0, peak.mz, peak.intensity);
                    subordinate
                        .metadata
                        .insert("FWHM".into(), float_meta(f64::from(*width)));
                    subordinates.push(subordinate);
                }
                update.subordinates = subordinates;
                updates.push(update);
            }
        }
        Ok((scored, updates))
    }

    /// The score of a scored spectrum: its second float data array's first
    /// value, as the source reads it.
    fn spectrum_score(spectrum: &MSSpectrum) -> Result<f64> {
        spectrum
            .float_data_arrays
            .get(1)
            .and_then(|array| array.data.first())
            .map(|score| f64::from(*score))
            .ok_or_else(|| {
                Error::InvalidValue(format!("spectrum '{}' carries no score", spectrum.name))
            })
    }

    /// Select the highest-scoring spectrum for each annotation (transition
    /// name), with its feature.
    ///
    /// Spectra scoring below `min_select_score` are skipped. Among the rest,
    /// the first spectrum of each name wins unless a later one scores strictly
    /// higher. The selection is returned in name order (the source's
    /// `std::map`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] when `scored_spectra` and `features`
    /// differ in size (the source's `Exception::InvalidSize`), and
    /// [`Error::InvalidValue`] for a spectrum without a score (the source reads
    /// out of bounds).
    pub fn select_spectra(
        &self,
        scored_spectra: &[MSSpectrum],
        features: &FeatureMap,
    ) -> Result<(Vec<MSSpectrum>, FeatureMap)> {
        if scored_spectra.len() != features.features.len() {
            return Err(Error::InvalidRange(format!(
                "scored_spectra size does not match features size ({} vs {})",
                scored_spectra.len(),
                features.features.len()
            )));
        }
        let best = self.best_per_name(scored_spectra)?;
        let spectra = best.values().map(|&i| scored_spectra[i].clone()).collect();
        let mut selected = FeatureMap::new();
        selected.features = best
            .values()
            .map(|&i| features.features[i].clone())
            .collect();
        Ok((spectra, selected))
    }

    /// [`select_spectra`](Self::select_spectra) without features: the source's
    /// two-argument overload.
    ///
    /// # Errors
    ///
    /// As [`select_spectra`](Self::select_spectra), without the feature count.
    pub fn select_spectra_only(&self, scored_spectra: &[MSSpectrum]) -> Result<Vec<MSSpectrum>> {
        let best = self.best_per_name(scored_spectra)?;
        Ok(best.values().map(|&i| scored_spectra[i].clone()).collect())
    }

    fn best_per_name(&self, scored_spectra: &[MSSpectrum]) -> Result<BTreeMap<String, usize>> {
        let mut best: BTreeMap<String, usize> = BTreeMap::new();
        for (index, spectrum) in scored_spectra.iter().enumerate() {
            let score = Self::spectrum_score(spectrum)?;
            if score < self.min_select_score {
                continue;
            }
            match best.get(&spectrum.name) {
                None => {
                    best.insert(spectrum.name.clone(), index);
                }
                Some(&current) => {
                    if Self::spectrum_score(&scored_spectra[current])? < score {
                        best.insert(spectrum.name.clone(), index);
                    }
                }
            }
        }
        Ok(best)
    }

    fn pick_all(&self, spectra: &[MSSpectrum]) -> Result<Vec<MSSpectrum>> {
        spectra
            .iter()
            .map(|spectrum| self.pick_spectrum(spectrum))
            .collect()
    }

    /// The whole pipeline for a target list: annotate, pick, drop the spectra
    /// whose picked spectrum is empty (with their features), score and select.
    ///
    /// # Errors
    ///
    /// The errors of the steps.
    pub fn extract_spectra(
        &self,
        experiment: &MSExperiment,
        targeted_exp: &TargetedExperiment,
    ) -> Result<(Vec<MSSpectrum>, FeatureMap)> {
        let (annotated, features, picked) =
            self.annotate_and_pick(experiment, targeted_exp, true)?;
        let mut features = features;
        let scored = self.score_spectra(&annotated, &picked, &mut features)?;
        self.select_spectra(&scored, &features)
    }

    /// [`extract_spectra`](Self::extract_spectra) without features: the
    /// source's three-argument overload.
    ///
    /// # Errors
    ///
    /// As [`extract_spectra`](Self::extract_spectra).
    pub fn extract_spectra_only(
        &self,
        experiment: &MSExperiment,
        targeted_exp: &TargetedExperiment,
    ) -> Result<Vec<MSSpectrum>> {
        let (annotated, _, picked) = self.annotate_and_pick(experiment, targeted_exp, false)?;
        let scored = self.score_spectra_only(&annotated, &picked)?;
        self.select_spectra_only(&scored)
    }

    fn annotate_and_pick(
        &self,
        experiment: &MSExperiment,
        targeted_exp: &TargetedExperiment,
        compute_features: bool,
    ) -> Result<(Vec<MSSpectrum>, FeatureMap, Vec<MSSpectrum>)> {
        let (annotated, features) =
            self.annotate(&experiment.spectra, targeted_exp, compute_features)?;
        let picked = self.pick_all(&annotated)?;
        let mut kept_annotated = Vec::new();
        let mut kept_picked = Vec::new();
        let mut kept_features = FeatureMap::new();
        for (index, (annotated, picked)) in annotated.into_iter().zip(picked).enumerate() {
            if picked.is_empty() {
                continue;
            }
            kept_annotated.push(annotated);
            kept_picked.push(picked);
            if compute_features {
                kept_features
                    .features
                    .push(features.features[index].clone());
            }
        }
        Ok((kept_annotated, kept_features, kept_picked))
    }

    /// The whole pipeline for MS1 features: annotate with
    /// [`annotate_spectra_with_ms1_features`](Self::annotate_spectra_with_ms1_features),
    /// pick, score without features and select with the MS2 features.
    ///
    /// Unlike the targeted pipeline, spectra whose picked spectrum is empty
    /// are kept, and score infinite or NaN, as in the source.
    ///
    /// # Errors
    ///
    /// The errors of the steps.
    pub fn extract_spectra_with_ms1_features(
        &self,
        experiment: &MSExperiment,
        ms1_features: &FeatureMap,
        ids: &mut UniqueIdGenerator,
    ) -> Result<(Vec<MSSpectrum>, FeatureMap)> {
        let (annotated, ms2_features) =
            self.annotate_spectra_with_ms1_features(&experiment.spectra, ms1_features, ids)?;
        let picked = self.pick_all(&annotated)?;
        let scored = self.score_spectra_only(&annotated, &picked)?;
        self.select_spectra(&scored, &ms2_features)
    }

    /// [`extract_spectra_with_ms1_features`](Self::extract_spectra_with_ms1_features)
    /// returning only the spectra: the source's three-argument overload.
    ///
    /// # Errors
    ///
    /// As [`extract_spectra_with_ms1_features`](Self::extract_spectra_with_ms1_features).
    pub fn extract_spectra_with_ms1_features_only(
        &self,
        experiment: &MSExperiment,
        ms1_features: &FeatureMap,
        ids: &mut UniqueIdGenerator,
    ) -> Result<Vec<MSSpectrum>> {
        let (annotated, _) =
            self.annotate_spectra_with_ms1_features(&experiment.spectra, ms1_features, ids)?;
        let picked = self.pick_all(&annotated)?;
        let scored = self.score_spectra_only(&annotated, &picked)?;
        self.select_spectra_only(&scored)
    }

    fn match_top(
        &self,
        input: &MSSpectrum,
        cmp: &dyn Comparator,
        top: usize,
    ) -> Result<Vec<Match>> {
        let mut scores = cmp.generate_scores(input, self.min_match_score)?;
        // std::sort with `a.second > b.second`: descending, ties in the order
        // libstdc++'s introsort leaves them.
        source_sort_by(&mut scores, |a, b| a.1 > b.1)?;
        let library = cmp.library();
        scores
            .iter()
            .take(top)
            .map(|&(index, score)| {
                library
                    .get(index)
                    .map(|spectrum| Match::new(spectrum.clone(), score))
                    .ok_or_else(|| {
                        Error::InvalidValue(format!("library index {index} out of range"))
                    })
            })
            .collect()
    }

    /// The best-scoring library matches of `input_spectrum`, at most
    /// `top_matches_to_report`, in descending score order.
    ///
    /// The comparator keeps matches scoring at least `min_match_score`; ties
    /// are ordered as the source's `std::sort` orders them.
    ///
    /// # Errors
    ///
    /// The comparator's errors, and [`Error::InvalidValue`] when it reports a
    /// library index outside its library.
    pub fn match_spectrum(
        &self,
        input_spectrum: &MSSpectrum,
        cmp: &dyn Comparator,
    ) -> Result<Vec<Match>> {
        self.match_top(input_spectrum, cmp, self.top_matches_to_report)
    }

    /// Compare spectra against a spectral library and record the best match in
    /// the corresponding features.
    ///
    /// Each feature gets the meta values `spectral_library_name` (the matched
    /// spectrum's name), `spectral_library_score` (in `[0, 1]`) and
    /// `spectral_library_comments` (the matched spectrum's `Comments` meta
    /// value); without a match they are `""`, `0.0` and `""`. The spectra are
    /// expected to come from [`extract_spectra`](Self::extract_spectra), that
    /// is, to be picked already.
    ///
    /// Returns the indices of the features without a match, which the source
    /// writes to its warning log.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] when `spectra` and `features` differ in
    /// size (the source's `Exception::InvalidSize`), [`Error::InvalidValue`]
    /// when a match's `Comments` is not a string (the source's
    /// `ConversionError`), and the comparator's errors. On error `features` is
    /// unchanged.
    pub fn targeted_matching(
        &self,
        spectra: &[MSSpectrum],
        cmp: &dyn Comparator,
        features: &mut FeatureMap,
    ) -> Result<Vec<usize>> {
        if spectra.len() != features.features.len() {
            return Err(Error::InvalidRange(format!(
                "spectra size does not match features size ({} vs {})",
                spectra.len(),
                features.features.len()
            )));
        }
        let mut no_matches = Vec::new();
        let mut values = Vec::with_capacity(spectra.len());
        for (index, spectrum) in spectra.iter().enumerate() {
            // The source sets top_matches_to_report_ to 1 for the duration.
            let matches = self.match_top(spectrum, cmp, 1)?;
            match matches.first() {
                Some(best) => {
                    let comments = if best.spectrum.metadata.contains_key("Comments") {
                        strict_string(&best.spectrum.metadata, "Comments")?
                    } else {
                        String::new()
                    };
                    values.push((best.spectrum.name.clone(), best.score, comments));
                }
                None => {
                    no_matches.push(index);
                    values.push((String::new(), 0.0, String::new()));
                }
            }
        }
        for (feature, (name, score, comments)) in features.features.iter_mut().zip(values) {
            feature
                .metadata
                .insert("spectral_library_name".into(), MetaValue::from(name));
            feature
                .metadata
                .insert("spectral_library_score".into(), float_meta(score));
            feature.metadata.insert(
                "spectral_library_comments".into(),
                MetaValue::from(comments),
            );
        }
        Ok(no_matches)
    }

    /// Compare unprocessed spectra against a spectral library and build a
    /// feature map of the results.
    ///
    /// The spectra are picked first ([`pick_spectrum`](Self::pick_spectrum));
    /// empty picked spectra are dropped. Each remaining one gets a feature at
    /// its retention time and first precursor m/z (0 without precursors), and
    /// [`targeted_matching`](Self::targeted_matching) records its best match.
    /// Returns the features and the indices (into the picked, non-empty
    /// spectra) without a match.
    ///
    /// # Errors
    ///
    /// The errors of picking and of [`targeted_matching`](Self::targeted_matching).
    pub fn untargeted_matching(
        &self,
        spectra: &[MSSpectrum],
        cmp: &dyn Comparator,
    ) -> Result<(FeatureMap, Vec<usize>)> {
        let picked: Vec<MSSpectrum> = self
            .pick_all(spectra)?
            .into_iter()
            .filter(|spectrum| !spectrum.is_empty())
            .collect();
        let mut features = FeatureMap::new();
        for spectrum in &picked {
            let feature = new_feature(spectrum.rt, first_precursor_mz(spectrum), 0.0);
            features.features.push(feature);
        }
        let unmatched = self.targeted_matching(&picked, cmp, &mut features)?;
        Ok((features, unmatched))
    }

    /// Build a transition list from MS1 features and their MS2 features.
    ///
    /// MS2 subordinates are grouped by their `PeptideRef`. Every MS1 feature
    /// gives a protein and a peptide with id `PeptideRef` (the first wins; a
    /// repeated peptide collects one more retention time), both carrying the
    /// feature's meta values, the peptide also its charge and the protein
    /// reference. Every MS2 feature of that `PeptideRef` whose m/z lies in
    /// `(min_fragment_mz, max_fragment_mz)` and below the MS1 m/z plus
    /// `relative_allowable_product_mass` gives a transition named by its
    /// `native_id`, with native id `<PeptideRef>_<native_id>_<RT>` (the
    /// retention time printed as `std::ostream` does), the MS1 m/z as precursor
    /// and the MS2 m/z as product, its intensity as library intensity and its
    /// meta values; the first transition of a native id wins. Proteins,
    /// peptides and transitions replace those of `t_exp`, in id order.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when an MS2 subordinate's `PeptideRef`
    /// or an MS2 feature's `native_id` is missing or not a string (the source's
    /// `ConversionError`), and when the result contains duplicate or invalid
    /// references (the source's `Exception::IllegalArgument` from
    /// `TransitionTSVFile::validateTargetedExperiment`). On error `t_exp` is
    /// unchanged; the source has already replaced its contents when it
    /// validates.
    pub fn construct_transitions_list(
        &self,
        ms1_features: &FeatureMap,
        ms2_features: &FeatureMap,
        t_exp: &mut TargetedExperiment,
    ) -> Result<()> {
        let mut ms1_to_ms2: BTreeMap<String, Vec<&Feature>> = BTreeMap::new();
        for feature in &ms2_features.features {
            for subordinate in &feature.subordinates {
                ms1_to_ms2
                    .entry(strict_string(&subordinate.metadata, "PeptideRef")?)
                    .or_default()
                    .push(subordinate);
            }
        }
        let mut transitions: BTreeMap<String, ReactionMonitoringTransition> = BTreeMap::new();
        let mut peptides: BTreeMap<String, Peptide> = BTreeMap::new();
        let mut proteins: BTreeMap<String, Protein> = BTreeMap::new();
        for ms1 in &ms1_features.features {
            let peptide_ref = lenient_string(&ms1.metadata, "PeptideRef");
            proteins.entry(peptide_ref.clone()).or_insert_with(|| {
                let mut protein = Protein::default();
                protein
                    .cv_terms
                    .metadata
                    .extend(ms1.metadata.iter().map(|(k, v)| (k.clone(), v.clone())));
                protein.id = peptide_ref.clone();
                protein
            });
            let mut rt = RetentionTime::default();
            rt.set_rt(ms1.rt);
            match peptides.get_mut(&peptide_ref) {
                Some(peptide) => peptide.rts.push(rt),
                None => {
                    let mut peptide = Peptide::default();
                    peptide.id = peptide_ref.clone();
                    peptide.set_charge_state(ms1.charge);
                    peptide
                        .cv_terms
                        .metadata
                        .extend(ms1.metadata.iter().map(|(k, v)| (k.clone(), v.clone())));
                    peptide.protein_refs.push(peptide_ref.clone());
                    peptide.rts.push(rt);
                    peptides.insert(peptide_ref.clone(), peptide);
                }
            }
            for ms2 in ms1_to_ms2
                .get(&peptide_ref)
                .map(Vec::as_slice)
                .unwrap_or(&[])
            {
                let current_mz = ms2.mz;
                if !(current_mz > self.min_fragment_mz
                    && current_mz < self.max_fragment_mz
                    && current_mz < ms1.mz + self.relative_allowable_product_mass)
                {
                    continue;
                }
                let native_id = format!(
                    "{peptide_ref}_{}_{}",
                    lenient_string(&ms2.metadata, "native_id"),
                    format_g(ms2.rt)
                );
                if transitions.contains_key(&native_id) {
                    continue;
                }
                let mut transition = ReactionMonitoringTransition::new();
                transition.set_library_intensity(f64::from(ms2.intensity));
                transition.set_name(strict_string(&ms2.metadata, "native_id")?);
                let mut rt = RetentionTime::default();
                rt.set_rt(ms2.rt);
                transition.set_retention_time(rt);
                transition.set_native_id(native_id.clone());
                transition.set_peptide_ref(peptide_ref.clone());
                transition.set_precursor_mz(ms1.mz);
                transition.set_product_mz(ms2.mz);
                transition
                    .cv_terms
                    .metadata
                    .extend(ms2.metadata.iter().map(|(k, v)| (k.clone(), v.clone())));
                transitions.insert(native_id, transition);
            }
        }
        let mut staged = t_exp.clone();
        staged.set_proteins(proteins.into_values().collect());
        staged.set_peptides(peptides.into_values().collect());
        staged.set_transitions(transitions.into_values().collect());
        if staged.contains_invalid_references() {
            return Err(Error::InvalidValue(
                "Invalid input, contains duplicate or invalid references".into(),
            ));
        }
        *t_exp = staged;
        Ok(())
    }

    /// Group features and their direct subordinates that carry both a
    /// `PeptideRef` and an `identifier` meta value by
    /// `<PeptideRef>_<identifier[0]>_<RT>`, the retention time printed by
    /// `std::to_string` (six decimals): the source's protected
    /// `organizeMapWithSameIdentifier`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a `PeptideRef` is not a string or
    /// an `identifier` is not a non-empty string list (the source's
    /// `ConversionError` and `std::out_of_range`).
    pub fn organize_map_with_same_identifier(
        &self,
        fmap_input: &FeatureMap,
    ) -> Result<BTreeMap<String, Vec<Feature>>> {
        let mut map: BTreeMap<String, Vec<Feature>> = BTreeMap::new();
        let mut add = |feature: &Feature| -> Result<()> {
            if feature.metadata.contains_key("PeptideRef")
                && feature.metadata.contains_key("identifier")
            {
                let peptide_ref = strict_string(&feature.metadata, "PeptideRef")?;
                let identifier = feature.metadata["identifier"]
                    .as_string_list()
                    .ok()
                    .and_then(|list| list.first())
                    .ok_or_else(|| {
                        Error::InvalidValue(
                            "meta value 'identifier' is not a non-empty string list".into(),
                        )
                    })?;
                let key = format!("{peptide_ref}_{identifier}_{}", to_string_fixed(feature.rt));
                map.entry(key).or_default().push(feature.clone());
            }
            Ok(())
        };
        for feature in &fmap_input.features {
            add(feature)?;
            for subordinate in &feature.subordinates {
                add(subordinate)?;
            }
        }
        Ok(map)
    }

    /// Combine features and subordinates with the same identifier into one
    /// feature each, appended to `fmap_output`.
    ///
    /// Groups come from
    /// [`organize_map_with_same_identifier`](Self::organize_map_with_same_identifier).
    /// Each group's feature gets a fresh unique id from `ids`; retention time
    /// and m/z averaged with weights `peak_apex_int / total` (or intensity
    /// where `peak_apex_int` is missing); the summed intensities and
    /// `peak_apex_int`; the `PeptideRef` before the first `_` of the group key
    /// (empty, as the source's error branch leaves it, when there is none); the
    /// first member's `scan_polarity` (empty when missing); and the group as
    /// subordinates.
    ///
    /// # Errors
    ///
    /// The errors of the grouping, and [`Error::InvalidValue`] for a
    /// `peak_apex_int` that is not numeric. The source catches every exception
    /// here, logs it and returns normally, having appended the groups
    /// processed so far (none, when the grouping fails); this returns the
    /// error and leaves `fmap_output` unchanged.
    pub fn merge_features(
        &self,
        fmap_input: &FeatureMap,
        fmap_output: &mut FeatureMap,
        ids: &mut UniqueIdGenerator,
    ) -> Result<()> {
        let map = self.organize_map_with_same_identifier(fmap_input)?;
        let apex = |feature: &Feature| -> Result<Option<f64>> {
            feature
                .metadata
                .get("peak_apex_int")
                .map(|value| numeric_meta(value, "peak_apex_int"))
                .transpose()
        };
        let mut merged = Vec::with_capacity(map.len());
        for (key, group) in &map {
            let mut total_intensity = 0.0;
            for feature in group {
                total_intensity += match apex(feature)? {
                    Some(value) => value,
                    None => f64::from(feature.intensity),
                };
            }
            let (mut rt, mut mz, mut intensity, mut peak_apex_int) = (0.0, 0.0, 0.0, 0.0);
            for feature in group {
                let apex_value = apex(feature)?;
                let weighting_factor = match apex_value {
                    Some(value) => value / total_intensity,
                    None => f64::from(feature.intensity) / total_intensity,
                };
                rt += feature.rt * weighting_factor;
                mz += feature.mz * weighting_factor;
                intensity += f64::from(feature.intensity);
                if let Some(value) = apex_value {
                    peak_apex_int += value;
                }
            }
            let mut feature = new_feature(0.0, 0.0, 0.0);
            feature.unique_id = ids.get_unique_id();
            let peptide_ref = key
                .split_once('_')
                .map(|(prefix, _)| prefix.to_owned())
                .unwrap_or_default();
            feature
                .metadata
                .insert("PeptideRef".into(), MetaValue::from(peptide_ref));
            feature.mz = mz;
            feature.rt = rt;
            let polarity = group
                .first()
                .and_then(|first| first.metadata.get("scan_polarity"))
                .cloned()
                .unwrap_or_default();
            feature.metadata.insert("scan_polarity".into(), polarity);
            feature.intensity = intensity as f32;
            feature
                .metadata
                .insert("peak_apex_int".into(), float_meta(peak_apex_int));
            feature.subordinates = group.clone();
            merged.push(feature);
        }
        fmap_output.features.extend(merged);
        Ok(())
    }

    /// Deisotope the MS2 spectra of `experiment` (MS1 spectra are left alone):
    /// the source's protected `deisotopeMS2Spectra_`, the first step of
    /// `storeSpectraMSP` when `deisotoping:use_deisotoper` is set.
    ///
    /// Runs `Deisotoper::deisotopeAndSingleCharge` with the `deisotoping:`
    /// parameters, without converting to single charge and with the source's
    /// remaining defaults (decreasing model, intensity check from the second
    /// peak, shared isotope extensions).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a charge or cluster size outside the
    /// deisotoper's range, and the deisotoper's errors (an unsupported
    /// tolerance, an unsorted spectrum). On error `experiment` is unchanged.
    pub fn deisotope_ms2_spectra(&self, experiment: &mut MSExperiment) -> Result<()> {
        let d = &self.deisotoping;
        let charge = |value: i32, key: &str| {
            u8::try_from(value).map_err(|_| {
                Error::InvalidValue(format!("deisotoping:{key} {value} is out of range"))
            })
        };
        let peaks = |value: i32, key: &str| {
            usize::try_from(value).map_err(|_| {
                Error::InvalidValue(format!("deisotoping:{key} {value} is out of range"))
            })
        };
        let deisotoper = Deisotoper {
            tolerance: if d.fragment_unit == "ppm" {
                Tolerance::Ppm(d.fragment_tolerance)
            } else {
                Tolerance::Absolute(d.fragment_tolerance)
            },
            min_charge: charge(d.min_charge, "min_charge")?,
            max_charge: charge(d.max_charge, "max_charge")?,
            keep_only_deisotoped: d.keep_only_deisotoped,
            min_isotope_peaks: peaks(d.min_isopeaks, "min_isopeaks")?,
            max_isotope_peaks: peaks(d.max_isopeaks, "max_isopeaks")?,
            make_single_charged: false,
            annotate_charge: d.annotate_charge,
            annotate_isotope_peak_count: false,
            annotate_features: false,
            use_decreasing_model: true,
            start_intensity_check: 2,
            add_up_intensity: false,
            allow_shared_isotopes: true,
            ..Deisotoper::default()
        };
        let mut spectra = experiment.spectra.clone();
        for spectrum in spectra.iter_mut().filter(|s| s.ms_level != 1) {
            *spectrum = deisotoper.deisotope(spectrum)?.spectrum;
        }
        experiment.spectra = spectra;
        Ok(())
    }

    /// Remove the peaks of MS2 spectra that lie above the precursor m/z plus
    /// `max_precursor_mass_threashold` ppm, and the peaks with an intensity
    /// below 1: the source's protected `removeMS2SpectraPeaks_`, the second
    /// step of `storeSpectraMSP`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an MS2 spectrum without a precursor
    /// (the source reads `getPrecursors()[0]` unchecked) and for a non-finite
    /// tolerance or precursor m/z. On error `experiment` is unchanged.
    pub fn remove_ms2_spectra_peaks(&self, experiment: &mut MSExperiment) -> Result<()> {
        let mut spectra = experiment.spectra.clone();
        for spectrum in spectra.iter_mut().filter(|s| s.ms_level != 1) {
            let precursor_mz =
                spectrum.precursors.first().map(|p| p.mz).ok_or_else(|| {
                    Error::InvalidValue("an MS2 spectrum has no precursor".into())
                })?;
            let mass_diff = ppm_to_mass(self.max_precursor_mass_threashold, precursor_mz)?;
            let keep: Vec<usize> = spectrum
                .peaks
                .iter()
                .enumerate()
                .filter(|(_, p)| p.mz <= precursor_mz + mass_diff && p.intensity >= 1.0)
                .map(|(i, _)| i)
                .collect();
            if keep.len() != spectrum.len() {
                spectrum.select(&keep)?;
            }
        }
        experiment.spectra = spectra;
        Ok(())
    }
}
