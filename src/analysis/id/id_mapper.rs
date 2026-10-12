// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Annotation of spectra, feature maps and consensus maps with peptide
//! identifications.
//!
//! Port of `ANALYSIS/ID/IDMapper.h` (`IDMapper.cpp`). Protein identifications
//! are assigned to the whole map; each peptide identification is matched by
//! retention time and m/z to spectra, features or consensus features. m/z on the
//! peptide side is either the precursor m/z stored with the identification or
//! the theoretical m/z of every peptide hit (parameter `mz_reference`).
//!
//! See `docs/ID_MAPPER_SUPPORT.md` for the member-by-member API mapping, the
//! preserved source conventions and the native differences.
//!
//! Every `annotate_*` entry point first *plans* its assignments without
//! touching the target map, checking the resource ceilings in
//! [`IDMapperLimits`](crate::analysis::id::id_mapper::IDMapperLimits) as it goes, and only then commits them. An error therefore
//! leaves the map exactly as it was; the source mutates as it goes.

use crate::chemistry::PeptideFragmentType;
use crate::concept::boost_regex::BoostRegex;
use crate::concept::constants::user_param::SPECTRUM_REFERENCE;
use crate::data_structures::DateTime;
use crate::identification::{IdentifierMSRunMapper, PeptideIdentification, ProteinIdentification};
use crate::kernel::MSExperiment;
use crate::kernel::features::{ConsensusFeature, ConsensusMap, Feature, FeatureMap};
use crate::metadata::{DataProcessing, MetaValue, ProcessingAction};
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet, HashMap};

/// Name of the source `DefaultParamHandler`, `IDMapper`.
pub const ID_MAPPER_NAME: &str = "IDMapper";

/// Search-run identifier the source invents when unidentified precursors must
/// be annotated but neither the map nor its unassigned identifications carry
/// one.
pub const UNKNOWN_SEARCH_RUN_IDENTIFIER: &str = "UNKNOWN_SEARCH_RUN_IDENTIFIER";

/// Unit of the m/z tolerance, the source `IDMapper::Measure`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Measure {
    /// Parts per million of the m/z being matched, the source `MEASURE_PPM`
    /// (parameter value `ppm`, the default).
    #[default]
    Ppm,
    /// Absolute Thomson, the source `MEASURE_DA` (parameter value `Da`).
    Da,
}

/// Source of m/z values on the peptide side (parameter `mz_reference`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MzReference {
    /// The precursor m/z stored with the peptide identification (`precursor`,
    /// the default). One m/z value per identification.
    #[default]
    Precursor,
    /// The theoretical m/z of each peptide hit, assuming protonation at the
    /// hit's charge (`peptide`). One m/z value per hit; an identification
    /// matches if any of them matches.
    Peptide,
}

/// Resource ceilings for one annotation call.
///
/// Native addition: the source allocates and loops without bounds. Every
/// ceiling is checked while the assignments are planned, before the target map
/// is touched, and exceeding one returns [`Error::InvalidValue`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IDMapperLimits {
    /// Most one-second retention-time bins in the feature-map hash table. The
    /// source sizes that table from the RT span of all features, so one feature
    /// at an absurd retention time would otherwise allocate without bound.
    pub max_rt_bins: usize,
    /// Most feature indices stored across all hash bins.
    pub max_hash_entries: usize,
    /// Most pairwise comparisons (identification or precursor against a
    /// candidate, or spectrum precursor against identification).
    pub max_comparisons: u64,
    /// Most identifications copied into the target, counting every copy.
    pub max_assignments: usize,
}

impl Default for IDMapperLimits {
    fn default() -> Self {
        Self {
            max_rt_bins: 10_000_000,
            max_hash_entries: 50_000_000,
            max_comparisons: 10_000_000_000,
            max_assignments: 10_000_000,
        }
    }
}

/// Partition of spectra by identification state, the source
/// `IDMapper::PeptideIdentificationListState`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PeptideIdentificationListState {
    /// Indices of spectra without precursors (for example MS1 spectra).
    pub no_precursors: Vec<usize>,
    /// Indices of spectra with at least one identified precursor.
    pub identified: Vec<usize>,
    /// Indices of spectra whose precursors are all unidentified.
    pub unidentified: Vec<usize>,
}

/// The identification half of the source `AnnotatedMSRun`: one peptide
/// identification per spectrum plus the run's protein identifications.
///
/// `METADATA/AnnotatedMSRun.h` is not ported; the spectra stay in their own
/// [`MSExperiment`] and [`IDMapper::annotate_spectra`] fills this record for
/// them. After a call with at least one identification,
/// `peptide_identifications` has exactly one entry per spectrum, as the source
/// resizes it.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SpectraAnnotation {
    /// One peptide identification per spectrum, the source
    /// `AnnotatedMSRun::getPeptideIdentifications`.
    pub peptide_identifications: Vec<PeptideIdentification>,
    /// Protein identifications of the run, the source
    /// `AnnotatedMSRun::getProteinIdentifications`.
    pub protein_identifications: Vec<ProteinIdentification>,
}

/// Counts the source logs after annotating spectra.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SpectraMappingStatistics {
    /// "Peptides assigned to a precursor": identifications whose hits were
    /// merged into at least one spectrum.
    pub assigned: usize,
    /// "Unassigned peptides": all identifications minus `assigned`.
    pub unassigned: usize,
    /// "Unmapped (empty) peptides": all identifications minus those queued for
    /// retention-time matching. As in the source, this also counts the
    /// identifications that were mapped by native ID.
    pub unmapped_empty: usize,
}

/// Counts and warnings the source logs after annotating a feature or
/// consensus map.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MapMappingStatistics {
    /// Identifications (with hits) assigned to no feature.
    pub peptides_unassigned: usize,
    /// Identifications assigned to exactly one feature.
    pub peptides_single: usize,
    /// Identifications assigned to more than one feature.
    pub peptides_multiple: usize,
    /// Unidentified precursors assigned to no feature. For a feature map, the
    /// source counts a centroid-m/z match here too, because that branch leaves
    /// the loop before its counter is incremented; this port keeps the count.
    pub precursors_unassigned: usize,
    /// Unidentified precursors assigned to exactly one feature.
    pub precursors_single: usize,
    /// Unidentified precursors assigned to more than one feature.
    pub precursors_multiple: usize,
    /// The source's `OPENMS_LOG_WARN` messages, in order.
    pub warnings: Vec<String>,
}

/// Assigns peptide and protein identifications to spectra, features and
/// consensus features by retention time and m/z.
///
/// Source `IDMapper`. Parameters (see [`IDMapper::defaults`]):
///
/// | parameter | field | default |
/// |---|---|---|
/// | `rt_tolerance` (s, ≥ 0) | [`rt_tolerance`](Self::rt_tolerance) | 5.0 |
/// | `mz_tolerance` (ppm or Da, ≥ 0) | [`mz_tolerance`](Self::mz_tolerance) | 20.0 |
/// | `mz_measure` (`ppm`, `Da`) | [`measure`](Self::measure) | `ppm` |
/// | `mz_reference` (`precursor`, `peptide`) | [`mz_reference`](Self::mz_reference) | `precursor` |
/// | `ignore_charge` (`true`, `false`) | [`ignore_charge`](Self::ignore_charge) | `false` |
///
/// Tolerances are "plus or minus x", so a matching window is twice as wide as
/// the tolerance.
///
/// ```rust
/// use openms::analysis::id::id_mapper::{IDMapper, Measure};
///
/// let mapper = IDMapper { mz_tolerance: 3.0, ..IDMapper::default() };
/// assert!((mapper.absolute_mz_tolerance(1000.0) - 0.003).abs() < 1e-12);
/// assert!(mapper.is_match(4.0, 1000.0, 1000.0028));
/// let da = IDMapper { measure: Measure::Da, ..mapper };
/// assert!(!da.is_match(5.0, 999.0, 1002.1));
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct IDMapper {
    /// Allowed retention-time deviation in seconds, the source `rt_tolerance_`.
    pub rt_tolerance: f64,
    /// Allowed m/z deviation in the unit of [`measure`](Self::measure), the
    /// source `mz_tolerance_`.
    pub mz_tolerance: f64,
    /// Unit of [`mz_tolerance`](Self::mz_tolerance), the source `measure_`.
    pub measure: Measure,
    /// Where peptide-side m/z values come from (parameter `mz_reference`; the
    /// source reads it from its parameter tree at each use).
    pub mz_reference: MzReference,
    /// For feature and consensus maps: assign identifications whatever their
    /// charge, the source `ignore_charge_`.
    pub ignore_charge: bool,
    /// Resource ceilings; native, not a source parameter.
    pub limits: IDMapperLimits,
}

impl Default for IDMapper {
    fn default() -> Self {
        Self {
            rt_tolerance: 5.0,
            mz_tolerance: 20.0,
            measure: Measure::Ppm,
            mz_reference: MzReference::Precursor,
            ignore_charge: false,
            limits: IDMapperLimits::default(),
        }
    }
}

/// Peptide-side values of one identification, the source `getIDDetails_`
/// output.
#[derive(Clone, Debug, Default)]
struct IdDetails {
    rt: f64,
    mz_values: Vec<f64>,
    charges: Vec<i32>,
}

/// A two-dimensional box, the source `DBoundingBox<2>` restricted to what the
/// mapper uses. Bounds are inclusive, as `DIntervalBase::encloses` is.
#[derive(Clone, Copy, Debug)]
struct Box2 {
    min_rt: f64,
    max_rt: f64,
    min_mz: f64,
    max_mz: f64,
}

impl Box2 {
    fn encloses(&self, rt: f64, mz: f64) -> bool {
        rt >= self.min_rt && rt <= self.max_rt && mz >= self.min_mz && mz <= self.max_mz
    }
}

fn limit(what: &str) -> Error {
    Error::InvalidValue(format!("IDMapper: {what} exceeds its configured limit"))
}

/// Source `File::basename`: the text after the last `/` or `\`. The crate's
/// `system::file::basename` is the same function, but `system` already depends
/// on `analysis`, so calling it from here would close a module cycle.
fn basename(file: &str) -> &str {
    file.rsplit(['/', '\\']).next().unwrap_or("")
}

fn charge_list_contains(charges: &[i32], charge: i32) -> bool {
    charges.contains(&charge)
}

/// Source `PeptideIdentification::empty`: no identifier, no hits, a zero
/// significance threshold, no score type and higher-is-better.
fn source_empty(id: &PeptideIdentification) -> bool {
    id.identifier.is_empty()
        && id.hits.is_empty()
        && id.significance_threshold == 0.0
        && id.score_type.is_empty()
        && id.higher_score_better
}

/// Source `StringUtils::toInt64`:
/// [`crate::concept::string_conversion::to_int64`].
fn source_to_int64(text: &str) -> Option<i64> {
    crate::concept::string_conversion::to_int64(text).ok()
}

/// Source `SpectrumNativeIDParser::extractScanNumber(native_id, regex,
/// no_error = false)`: the last token of capture group 1 over all matches,
/// converted with `StringUtils::toInt32`.
fn extract_scan_number(native_id: &str, regex: &BoostRegex) -> Result<i32> {
    let haystack = native_id.as_bytes();
    let mut last = None;
    for token in regex.tokens(haystack, &[1])? {
        let token = token?;
        last = Some(String::from_utf8_lossy(token.as_bytes(haystack)).into_owned());
    }
    if let Some(value) = last {
        if let Some(number) = source_to_int64(&value).and_then(|v| i32::try_from(v).ok()) {
            return Ok(number);
        }
    }
    Err(Error::Parse {
        line: 0,
        message: format!("Could not extract scan number from '{native_id}'"),
    })
}

impl IDMapper {
    /// The default parameters (source `getDefaults`), in source order.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure, which the constant tree cannot
    /// trigger.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters of this mapper in the layout of [`IDMapper::defaults`],
    /// carrying the current values (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        let mut rt = ParamEntry {
            name: "rt_tolerance".into(),
            description: "RT tolerance (in seconds) for the matching".into(),
            value: ParamValue::Float(self.rt_tolerance),
            ..ParamEntry::default()
        };
        rt.min_float = 0.0;
        let mut mz = ParamEntry {
            name: "mz_tolerance".into(),
            description: "m/z tolerance (in ppm or Da) for the matching".into(),
            value: ParamValue::Float(self.mz_tolerance),
            ..ParamEntry::default()
        };
        mz.min_float = 0.0;
        let string_entry =
            |name: &str, value: &str, description: &str, valid: &[&str]| ParamEntry {
                name: name.into(),
                description: description.into(),
                value: ParamValue::String(value.into()),
                valid_strings: valid.iter().map(|v| (*v).to_owned()).collect(),
                ..ParamEntry::default()
            };
        let entries = vec![
            rt,
            mz,
            string_entry(
                "mz_measure",
                match self.measure {
                    Measure::Ppm => "ppm",
                    Measure::Da => "Da",
                },
                "unit of 'mz_tolerance' (ppm or Da)",
                &["ppm", "Da"],
            ),
            string_entry(
                "mz_reference",
                match self.mz_reference {
                    MzReference::Precursor => "precursor",
                    MzReference::Peptide => "peptide",
                },
                "source of m/z values for peptide identifications",
                &["precursor", "peptide"],
            ),
            string_entry(
                "ignore_charge",
                if self.ignore_charge { "true" } else { "false" },
                "For feature/consensus maps: Assign an ID independently of whether its charge state matches that of the (consensus) feature.",
                &["true", "false"],
            ),
        ];
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries,
            nodes: Vec::new(),
        })
    }

    /// Build a mapper from parameters, as source `setParameters` followed by
    /// `updateMembers_`.
    ///
    /// Missing parameters take their defaults; [`limits`](Self::limits) keeps
    /// its default. Unknown parameters are ignored, as the source only warns
    /// (see [`IDMapper::from_param_with_warnings`]); the upstream class test
    /// itself sets `consensus:use_subelements`, which is not an `IDMapper`
    /// parameter.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong value type
    /// or violates its restriction (a negative tolerance, an unknown unit), the
    /// source `Exception::InvalidParameter`.
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`IDMapper::from_param`], also returning the unknown-parameter
    /// warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`IDMapper::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(ID_MAPPER_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, Self::from_complete_param)
    }

    fn from_complete_param(param: &Param) -> Result<Self> {
        let float = |key: &str| -> Result<f64> {
            match param.value(key)? {
                ParamValue::Float(value) => Ok(*value),
                _ => Err(Error::InvalidValue(format!(
                    "parameter '{key}' must be a floating-point value"
                ))),
            }
        };
        Ok(Self {
            rt_tolerance: float("rt_tolerance")?,
            mz_tolerance: float("mz_tolerance")?,
            measure: if param.value("mz_measure")?.as_str()? == "ppm" {
                Measure::Ppm
            } else {
                Measure::Da
            },
            mz_reference: if param.value("mz_reference")?.as_str()? == "peptide" {
                MzReference::Peptide
            } else {
                MzReference::Precursor
            },
            ignore_charge: param.value("ignore_charge")?.as_str()? == "true",
            limits: IDMapperLimits::default(),
        })
    }

    /// Absolute m/z tolerance in Thomson at `mz`, the source protected
    /// `getAbsoluteMZTolerance_`: `mz * mz_tolerance / 1e6` for ppm, the
    /// tolerance itself for Da.
    ///
    /// The source's `Exception::InvalidValue` for an illegal internal measure
    /// cannot occur: [`Measure`] has only the two values.
    pub fn absolute_mz_tolerance(&self, mz: f64) -> f64 {
        match self.measure {
            Measure::Ppm => mz * self.mz_tolerance / 1e6,
            Measure::Da => self.mz_tolerance,
        }
    }

    /// Whether a retention-time distance and an m/z pair are within tolerance,
    /// the source protected `isMatch_`.
    ///
    /// Both bounds are inclusive. For ppm the deviation is
    /// `|mz_observed - mz_theoretical| / mz_theoretical * 1e6`
    /// (`Math::getPPMAbs(mz_observed, mz_theoretical)`); a zero or non-finite
    /// reference makes that value non-finite, which, as in the source's
    /// floating-point comparison, never matches.
    pub fn is_match(&self, rt_distance: f64, mz_theoretical: f64, mz_observed: f64) -> bool {
        if rt_distance.abs() > self.rt_tolerance || rt_distance.is_nan() {
            return false;
        }
        match self.measure {
            Measure::Ppm => {
                let ppm = ((mz_observed - mz_theoretical) / mz_theoretical * 1e6).abs();
                ppm <= self.mz_tolerance
            }
            Measure::Da => (mz_theoretical - mz_observed).abs() <= self.mz_tolerance,
        }
    }

    /// Partition spectra into those without precursors, those with at least one
    /// identified precursor and those with only unidentified precursors, the
    /// source static `mapPrecursorsToIdentifications`.
    ///
    /// A precursor counts as identified when an identification *with hits*
    /// lies strictly within `mz_tol` of its m/z and strictly within `rt_tol` of
    /// the spectrum's retention time. The source's defaults are 0.001 for both
    /// ([`IDMapper::map_precursors_to_identifications_default`]); in principle
    /// they should be zero up to numerical noise. An identification without
    /// m/z or RT never matches (the source compares NaN).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when spectra × precursors ×
    /// identifications exceeds the default
    /// [`IDMapperLimits::max_comparisons`].
    pub fn map_precursors_to_identifications(
        spectra: &MSExperiment,
        ids: &[PeptideIdentification],
        mz_tol: f64,
        rt_tol: f64,
    ) -> Result<PeptideIdentificationListState> {
        map_precursors(
            spectra,
            ids,
            mz_tol,
            rt_tol,
            IDMapperLimits::default().max_comparisons,
        )
    }

    /// [`IDMapper::map_precursors_to_identifications`] with the source default
    /// tolerances of 0.001 Th and 0.001 s.
    ///
    /// # Errors
    ///
    /// As [`IDMapper::map_precursors_to_identifications`].
    pub fn map_precursors_to_identifications_default(
        spectra: &MSExperiment,
        ids: &[PeptideIdentification],
    ) -> Result<PeptideIdentificationListState> {
        Self::map_precursors_to_identifications(spectra, ids, 0.001, 0.001)
    }

    /// Append one identification data-processing record per protein
    /// identification run, the source static `addIdentificationDataProcessing_`.
    ///
    /// Each record names the search engine and its version as the software,
    /// carries the run's date as completion time and the action
    /// [`ProcessingAction::Identification`], and the meta values
    /// `parameter: db` and `parameter: db_version` when the run's database and
    /// database version are not empty.
    ///
    /// The source copies a `DateTime` object; the port stores the run date as
    /// text, so it is parsed. A run without a date, or with one the source
    /// `DateTime` would not parse, gets no completion time — the source's
    /// unset `DateTime` has no representation either.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a meta value cannot be built.
    pub fn add_identification_data_processing(
        data_processing: &mut Vec<DataProcessing>,
        protein_ids: &[ProteinIdentification],
    ) -> Result<()> {
        let records = identification_data_processing(protein_ids)?;
        data_processing.extend(records);
        Ok(())
    }

    /// Whether average peptide masses should be used for matching against the
    /// m/z values reported for features, the source protected `checkMassType_`.
    ///
    /// Looks at the `parameter: algorithm:feature:reported_mz` meta value of
    /// every processing record whose software is `FeatureFinder`: `average`
    /// and `maximum` select average masses (the latter with a warning), any
    /// other value monoisotopic ones. When two records disagree, the source
    /// warns and falls back to monoisotopic masses. The warnings are appended
    /// to `warnings`.
    pub fn uses_average_mass(processing: &[DataProcessing], warnings: &mut Vec<String>) -> bool {
        check_mass_type(processing.iter(), warnings)
    }

    fn check_ids(ids: &[PeptideIdentification]) -> Result<()> {
        for id in ids {
            match id.rt {
                None => {
                    return Err(Error::MissingInformation(
                        "IDMapper: 'RT' information missing for peptide identification!".into(),
                    ));
                }
                Some(rt) if !rt.is_finite() => {
                    return Err(Error::InvalidValue(
                        "IDMapper: peptide identification RT is not finite".into(),
                    ));
                }
                Some(_) => {}
            }
            match id.mz {
                None => {
                    return Err(Error::MissingInformation(
                        "IDMapper: 'MZ' information missing for peptide identification!".into(),
                    ));
                }
                Some(mz) if !mz.is_finite() => {
                    return Err(Error::InvalidValue(
                        "IDMapper: peptide identification m/z is not finite".into(),
                    ));
                }
                Some(_) => {}
            }
        }
        Ok(())
    }

    /// Source `getIDDetails_`.
    fn id_details(&self, id: &PeptideIdentification, use_avg_mass: bool) -> Result<IdDetails> {
        let mut details = IdDetails {
            rt: id.rt.unwrap_or(f64::NAN),
            ..IdDetails::default()
        };
        if self.mz_reference == MzReference::Precursor {
            details.mz_values.push(id.mz.unwrap_or(f64::NAN));
        }
        for hit in &id.hits {
            details.charges.push(hit.charge);
            if self.mz_reference == MzReference::Peptide {
                let mass = if use_avg_mass {
                    hit.sequence
                        .average_mass_for(PeptideFragmentType::Full, hit.charge)?
                } else {
                    hit.sequence
                        .mono_mass_for(PeptideFragmentType::Full, hit.charge)?
                };
                // The source divides by the charge unguarded: charge zero gives
                // an infinite m/z that never matches, and so does this.
                details.mz_values.push(mass / f64::from(hit.charge));
            }
        }
        Ok(details)
    }

    /// Source `increaseBoundingBox_`.
    fn increase_box(&self, b: &mut Box2) {
        let sub_mz = self.absolute_mz_tolerance(b.min_mz);
        let add_mz = self.absolute_mz_tolerance(b.max_mz);
        b.min_rt -= self.rt_tolerance;
        b.min_mz -= sub_mz;
        b.max_rt += self.rt_tolerance;
        b.max_mz += add_mz;
    }

    /// Annotate spectra with peptide identifications, the source
    /// `annotate(AnnotatedMSRun&, const PeptideIdentificationList&, const
    /// std::vector<ProteinIdentification>&, bool clear_ids, bool map_ms1)`.
    ///
    /// With `clear_ids`, the peptide and protein identifications of
    /// `annotation` are cleared first. When `peptide_ids` is empty nothing else
    /// happens. Otherwise `protein_ids` are appended to the run's protein
    /// identifications and `annotation.peptide_identifications` is resized to
    /// one entry per spectrum (new entries are empty identifications).
    ///
    /// Identifications that are not empty (source
    /// `PeptideIdentification::empty`) are then mapped:
    ///
    /// - one carrying a `spectrum_reference` meta value that names a spectrum's
    ///   native ID has its hits appended to that spectrum's identification (the
    ///   last spectrum with that native ID wins, as in `SpectrumLookup`);
    /// - all others are matched by retention time: every spectrum whose RT is
    ///   *strictly* within `rt_tolerance` of the identification's RT, and — unless
    ///   `map_ms1` — one of whose precursors matches the identification's m/z
    ///   ([`IDMapper::is_match`] with the identification as the theoretical
    ///   value), receives a copy of its hits. A PeptideIdentification is added
    ///   to *all* spectra within the RT and m/z boundaries.
    ///
    /// The RT sweep is the source's two-pointer walk over both RT-sorted lists
    /// and keeps its quirks: its comparisons are strict, so with `rt_tolerance`
    /// zero an identification at exactly a spectrum's RT is not mapped, and
    /// after its single corrective step to the right the walk does not
    /// re-check the lower end of the window.
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when an identification lacks RT or
    /// m/z (the source `Exception::MissingInformation`),
    /// [`Error::InvalidValue`] when either is not finite, when a
    /// `spectrum_reference` meta value is not a string (the source's implicit
    /// conversion throws), or when a limit is exceeded. `annotation` is
    /// unchanged on error.
    pub fn annotate_spectra(
        &self,
        experiment: &MSExperiment,
        annotation: &mut SpectraAnnotation,
        peptide_ids: &[PeptideIdentification],
        protein_ids: &[ProteinIdentification],
        clear_ids: bool,
        map_ms1: bool,
    ) -> Result<SpectraMappingStatistics> {
        Self::check_ids(peptide_ids)?;
        if peptide_ids.is_empty() {
            if clear_ids {
                annotation.peptide_identifications.clear();
                annotation.protein_identifications.clear();
            }
            return Ok(SpectraMappingStatistics::default());
        }
        let spectra = &experiment.spectra;
        // Plan: (spectrum index, identification index) merges, in source order.
        let mut merges: Vec<(usize, usize)> = Vec::new();
        let mut mapped = BTreeSet::new();
        let mut rt_candidates: Vec<(f64, usize)> = Vec::new();
        let mut native_ids: Option<HashMap<&str, usize>> = None;
        for (i, id) in peptide_ids.iter().enumerate() {
            if source_empty(id) {
                continue;
            }
            match id.metadata.get(SPECTRUM_REFERENCE) {
                None => rt_candidates.push((id.rt.unwrap_or(f64::NAN), i)),
                Some(value) => {
                    let native_id = value.as_str().map_err(|_| {
                        Error::InvalidValue(
                            "IDMapper: 'spectrum_reference' meta value is not a string".into(),
                        )
                    })?;
                    let lookup = native_ids.get_or_insert_with(|| {
                        spectra
                            .iter()
                            .enumerate()
                            .map(|(index, s)| (s.native_id.as_str(), index))
                            .collect()
                    });
                    match lookup.get(native_id) {
                        Some(&index) => {
                            merges.push((index, i));
                            mapped.insert(i);
                        }
                        None => rt_candidates.push((id.rt.unwrap_or(f64::NAN), i)),
                    }
                }
            }
        }
        if !rt_candidates.is_empty() {
            // multimap order: by key, ties in insertion order (stable sort).
            rt_candidates.sort_by(|a, b| a.0.total_cmp(&b.0));
            let mut experiment_rts: Vec<(f64, usize)> = spectra
                .iter()
                .enumerate()
                .map(|(index, s)| (s.rt, index))
                .collect();
            if experiment_rts.iter().any(|(rt, _)| rt.is_nan()) {
                return Err(Error::InvalidValue(
                    "IDMapper: spectrum retention time is NaN".into(),
                ));
            }
            experiment_rts.sort_by(|a, b| a.0.total_cmp(&b.0));
            let tol = self.rt_tolerance;
            let end = rt_candidates.len();
            let mut j = 0usize;
            let mut work = 0u64;
            for &(spectrum_rt, spectrum_index) in &experiment_rts {
                if j == end {
                    j -= 1; // valid: at least one candidate
                }
                while j != 0 && (spectrum_rt - rt_candidates[j].0) < tol {
                    j -= 1;
                }
                if j != end && (spectrum_rt - rt_candidates[j].0) > tol {
                    j += 1;
                }
                if j == end {
                    break;
                }
                while j != end && (rt_candidates[j].0 - spectrum_rt) < tol {
                    let id_index = rt_candidates[j].1;
                    let precursors = &spectra[spectrum_index].precursors;
                    work = work.saturating_add(precursors.len() as u64 + 1);
                    if work > self.limits.max_comparisons {
                        return Err(limit("comparison count"));
                    }
                    let id_mz = peptide_ids[id_index].mz.unwrap_or(f64::NAN);
                    let success =
                        map_ms1 || precursors.iter().any(|p| self.is_match(0.0, id_mz, p.mz));
                    if success {
                        if merges.len() >= self.limits.max_assignments {
                            return Err(limit("assignment count"));
                        }
                        merges.push((spectrum_index, id_index));
                        mapped.insert(id_index);
                    }
                    j += 1;
                }
            }
        }
        // Commit.
        if clear_ids {
            annotation.peptide_identifications.clear();
            annotation.protein_identifications.clear();
        }
        annotation
            .protein_identifications
            .extend(protein_ids.iter().cloned());
        annotation
            .peptide_identifications
            .resize(spectra.len(), PeptideIdentification::default());
        for (spectrum_index, id_index) in merges {
            annotation.peptide_identifications[spectrum_index]
                .hits
                .extend(peptide_ids[id_index].hits.iter().cloned());
        }
        Ok(SpectraMappingStatistics {
            assigned: mapped.len(),
            unassigned: peptide_ids.len() - mapped.len(),
            unmapped_empty: peptide_ids.len() - rt_candidates.len(),
        })
    }

    /// Annotate spectra with the peptide identifications of a feature map, the
    /// source `annotate(AnnotatedMSRun&, const FeatureMap&, bool clear_ids,
    /// bool map_ms1)`.
    ///
    /// Collects every peptide identification of every feature (not the
    /// unassigned ones), filling a missing RT or m/z from the feature, and
    /// forwards them with the map's protein identifications to
    /// [`IDMapper::annotate_spectra`].
    ///
    /// # Errors
    ///
    /// As [`IDMapper::annotate_spectra`].
    pub fn annotate_spectra_from_features(
        &self,
        experiment: &MSExperiment,
        annotation: &mut SpectraAnnotation,
        feature_map: &FeatureMap,
        clear_ids: bool,
        map_ms1: bool,
    ) -> Result<SpectraMappingStatistics> {
        let total: usize = feature_map
            .features
            .iter()
            .map(|f| f.peptide_identifications.len())
            .sum();
        if total > self.limits.max_assignments {
            return Err(limit("identification count"));
        }
        let mut peptide_ids = Vec::with_capacity(total);
        for feature in &feature_map.features {
            for id in &feature.peptide_identifications {
                let mut id = id.clone();
                if id.mz.is_none() {
                    id.mz = Some(feature.mz);
                }
                if id.rt.is_none() {
                    id.rt = Some(feature.rt);
                }
                peptide_ids.push(id);
            }
        }
        self.annotate_spectra(
            experiment,
            annotation,
            &peptide_ids,
            &feature_map.protein_identifications,
            clear_ids,
            map_ms1,
        )
    }

    /// Annotate a feature map with peptide identifications, the source
    /// `annotate(FeatureMap&, const PeptideIdentificationList&, const
    /// std::vector<ProteinIdentification>&, bool use_centroid_rt, bool
    /// use_centroid_mz, const PeakMap& spectra)`.
    ///
    /// `protein_ids` are appended to the map, and one identification
    /// data-processing record per run is appended to its processing
    /// ([`IDMapper::add_identification_data_processing`]).
    ///
    /// If *all* features have at least one convex hull, identification
    /// positions are matched against the bounding boxes of the convex hulls;
    /// otherwise the feature centroids are used (with a warning). Centroid
    /// coordinates also replace the corresponding hull ranges when
    /// `use_centroid_rt` or `use_centroid_mz` is set. Each box is widened by the
    /// RT and m/z tolerances on both sides. When the hull boxes are used, a
    /// position inside a feature's overall box must also fall inside the
    /// (widened) box of one of its individual mass-trace hulls.
    ///
    /// With `mz_reference = peptide` and centroid m/z matching, the processing
    /// history decides between monoisotopic and average peptide masses
    /// ([`IDMapper::uses_average_mass`]).
    ///
    /// Unless `ignore_charge`, a feature's charge must be among the hits'
    /// charges (precursor m/z) or equal the matched hit's charge (peptide m/z).
    /// An identification overlapping several features is annotated to all of
    /// them; one matching none is appended to the unassigned identifications.
    /// Identifications without hits are skipped.
    ///
    /// When `spectra` is given, every MS/MS precursor that no identification
    /// explains ([`IDMapper::map_precursors_to_identifications_default`]) is
    /// added to each feature covering it as an empty identification carrying
    /// the precursor's RT and m/z, the meta values `spectrum_index` and (when
    /// the native ID is not empty) `spectrum_reference`, and the identifier of
    /// the map's first protein run, else of its first unassigned
    /// identification, else [`UNKNOWN_SEARCH_RUN_IDENTIFIER`] — in which case a
    /// protein run with that identifier and the current date is added to the
    /// map. Here the charge check is `|precursor charge| == |feature charge|`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when an identification lacks RT or
    /// m/z, [`Error::InvalidValue`] when a coordinate is not finite, when a
    /// feature has convex hulls but none with a point (the source then indexes
    /// its hash table out of range), when a peptide mass cannot be computed,
    /// or when a limit is exceeded. The map is unchanged on error.
    pub fn annotate_features(
        &self,
        map: &mut FeatureMap,
        ids: &[PeptideIdentification],
        protein_ids: &[ProteinIdentification],
        use_centroid_rt: bool,
        use_centroid_mz: bool,
        spectra: &MSExperiment,
    ) -> Result<MapMappingStatistics> {
        Self::check_ids(ids)?;
        check_spectra(spectra)?;
        let mut stats = MapMappingStatistics::default();
        let new_processing = identification_data_processing(protein_ids)?;

        let (mut centroid_rt, mut centroid_mz) = (use_centroid_rt, use_centroid_mz);
        if !(centroid_rt && centroid_mz) && map.features.iter().any(|f| f.convex_hulls.is_empty()) {
            centroid_rt = true;
            centroid_mz = true;
            stats.warnings.push(
                "IDMapper warning: at least one feature has no convex hull - using centroid coordinates for matching".into(),
            );
        }
        let mut use_avg_mass = false;
        if centroid_mz && self.mz_reference == MzReference::Peptide {
            use_avg_mass = check_mass_type(
                map.data_processing.iter().chain(new_processing.iter()),
                &mut stats.warnings,
            );
        }

        // Feature boxes.
        let mut boxes = Vec::with_capacity(map.features.len());
        let mut min_rt = f64::MAX;
        let mut max_rt = -f64::MAX;
        for feature in &map.features {
            if !feature.rt.is_finite() || !feature.mz.is_finite() {
                return Err(Error::InvalidValue(
                    "IDMapper: feature coordinates must be finite".into(),
                ));
            }
            let mut b = Box2 {
                min_rt: feature.rt,
                max_rt: feature.rt,
                min_mz: feature.mz,
                max_mz: feature.mz,
            };
            if !(centroid_rt && centroid_mz) {
                let hull = feature.hull_bounding_box().ok_or_else(|| {
                    Error::InvalidValue(
                        "IDMapper: a feature has convex hulls but none contains a point".into(),
                    )
                })?;
                b = Box2 {
                    min_rt: hull.min().rt,
                    max_rt: hull.max().rt,
                    min_mz: hull.min().mz,
                    max_mz: hull.max().mz,
                };
            }
            if centroid_rt {
                b.min_rt = feature.rt;
                b.max_rt = feature.rt;
            }
            if centroid_mz {
                b.min_mz = feature.mz;
                b.max_mz = feature.mz;
            }
            self.increase_box(&mut b);
            if !(b.min_rt.is_finite() && b.max_rt.is_finite()) {
                return Err(Error::InvalidValue(
                    "IDMapper: widened feature RT range is not finite".into(),
                ));
            }
            min_rt = min_rt.min(b.min_rt);
            max_rt = max_rt.max(b.max_rt);
            boxes.push(b);
        }

        // RT hash table of one-second bins.
        let mut offset = 0i64;
        let mut hash_table: Vec<Vec<usize>> = Vec::new();
        if map.features.is_empty() {
            stats.warnings.push(
                "IDMapper received an empty FeatureMap! All peptides are mapped as 'unassigned'!"
                    .into(),
            );
        } else {
            let bin = |rt: f64| -> Result<i64> {
                let floor = rt.floor();
                if !(-9.0e18..=9.0e18).contains(&floor) {
                    return Err(limit("retention-time bin index"));
                }
                Ok(floor as i64)
            };
            offset = bin(min_rt)?;
            let bins = bin(max_rt)?
                .checked_sub(offset)
                .and_then(|n| n.checked_add(1))
                .and_then(|n| usize::try_from(n).ok())
                .ok_or_else(|| limit("retention-time bin count"))?;
            if bins > self.limits.max_rt_bins {
                return Err(limit("retention-time bin count"));
            }
            let mut entries = 0usize;
            for b in &boxes {
                let span = usize::try_from(bin(b.max_rt)? - bin(b.min_rt)? + 1).unwrap_or(0);
                entries = entries.saturating_add(span);
            }
            if entries > self.limits.max_hash_entries {
                return Err(limit("hash table entry count"));
            }
            hash_table = vec![Vec::new(); bins];
            for (index, b) in boxes.iter().enumerate() {
                for i in bin(b.min_rt)?..=bin(b.max_rt)? {
                    hash_table[(i - offset) as usize].push(index);
                }
            }
        }
        let bin_of = |rt: f64| (rt.floor() as i64 - offset) as usize;
        let mut work = 0u64;
        let mut count_work = |n: usize| -> Result<()> {
            work = work.saturating_add(n as u64);
            if work > self.limits.max_comparisons {
                Err(limit("comparison count"))
            } else {
                Ok(())
            }
        };

        // Plan identification assignments.
        let mut assignments: Vec<(usize, usize)> = Vec::new(); // (feature, id)
        let mut unassigned_ids: Vec<usize> = Vec::new();
        for (id_index, id) in ids.iter().enumerate() {
            if id.hits.is_empty() {
                continue;
            }
            let details = self.id_details(id, use_avg_mass)?;
            let rt = details.rt;
            if rt < min_rt || rt > max_rt {
                unassigned_ids.push(id_index);
                stats.peptides_unassigned += 1;
                continue;
            }
            let candidates = &hash_table[bin_of(rt)];
            count_work(
                candidates
                    .len()
                    .saturating_mul(details.mz_values.len().max(1)),
            )?;
            let mut matching = 0usize;
            for &feature_index in candidates {
                let feature = &map.features[feature_index];
                let mut check_charge = !self.ignore_charge;
                if check_charge && details.mz_values.len() == 1 {
                    if !charge_list_contains(&details.charges, feature.charge) {
                        continue;
                    }
                    check_charge = false;
                }
                for (l_index, &mz) in details.mz_values.iter().enumerate() {
                    if check_charge && details.charges[l_index] != feature.charge {
                        continue;
                    }
                    if boxes[feature_index].encloses(rt, mz) {
                        if centroid_mz {
                            assignments.push((feature_index, id_index));
                            matching += 1;
                            break;
                        }
                        if self.hull_match(feature, centroid_rt, rt, mz) {
                            assignments.push((feature_index, id_index));
                            matching += 1;
                            break;
                        }
                    }
                }
                if assignments.len() > self.limits.max_assignments {
                    return Err(limit("assignment count"));
                }
            }
            match matching {
                0 => {
                    unassigned_ids.push(id_index);
                    stats.peptides_unassigned += 1;
                }
                1 => stats.peptides_single += 1,
                _ => stats.peptides_multiple += 1,
            }
        }

        // Unidentified precursors.
        let unidentified =
            map_precursors(spectra, ids, 0.001, 0.001, self.limits.max_comparisons)?.unidentified;
        let (run_identifier, add_unknown_run) = empty_run_identifier(
            !unidentified.is_empty(),
            map.protein_identifications
                .iter()
                .chain(protein_ids.iter())
                .next(),
            map.unassigned_peptide_identifications
                .first()
                .map(|id| id.identifier.as_str())
                .or_else(|| unassigned_ids.first().map(|&i| ids[i].identifier.as_str())),
        );
        let mut precursor_assignments: Vec<(usize, PeptideIdentification)> = Vec::new();
        for &spectrum_index in &unidentified {
            let spectrum = &spectra.spectra[spectrum_index];
            for precursor in &spectrum.precursors {
                let (mz, rt, z) = (precursor.mz, spectrum.rt, precursor.charge);
                if rt < min_rt || rt > max_rt {
                    stats.precursors_unassigned += 1;
                    continue;
                }
                let empty_id = precursor_empty_id(
                    spectrum_index,
                    rt,
                    mz,
                    &spectrum.native_id,
                    &run_identifier,
                )?;
                let candidates = &hash_table[bin_of(rt)];
                count_work(candidates.len())?;
                let mut matching = 0usize;
                for &feature_index in candidates {
                    let feature = &map.features[feature_index];
                    if !self.ignore_charge && z.unsigned_abs() != feature.charge.unsigned_abs() {
                        continue;
                    }
                    if boxes[feature_index].encloses(rt, mz) {
                        if centroid_mz {
                            // Source quirk: this branch leaves the candidate
                            // loop without counting the match.
                            precursor_assignments.push((feature_index, empty_id.clone()));
                            break;
                        }
                        if self.hull_match(feature, centroid_rt, rt, mz) {
                            precursor_assignments.push((feature_index, empty_id.clone()));
                            matching += 1;
                            break;
                        }
                    }
                }
                if assignments.len() + precursor_assignments.len() > self.limits.max_assignments {
                    return Err(limit("assignment count"));
                }
                match matching {
                    0 => stats.precursors_unassigned += 1,
                    1 => stats.precursors_single += 1,
                    _ => stats.precursors_multiple += 1,
                }
            }
        }

        // Commit.
        map.protein_identifications
            .extend(protein_ids.iter().cloned());
        map.data_processing.extend(new_processing);
        for (feature_index, id_index) in assignments {
            map.features[feature_index]
                .peptide_identifications
                .push(ids[id_index].clone());
        }
        for id_index in unassigned_ids {
            map.unassigned_peptide_identifications
                .push(ids[id_index].clone());
        }
        if add_unknown_run {
            map.protein_identifications.push(unknown_run());
        }
        for (feature_index, id) in precursor_assignments {
            map.features[feature_index].peptide_identifications.push(id);
        }
        Ok(stats)
    }

    /// The per-mass-trace check of the feature-map mapping: whether `(rt, mz)`
    /// lies in the widened bounding box of one of the feature's hulls.
    fn hull_match(&self, feature: &Feature, centroid_rt: bool, rt: f64, mz: f64) -> bool {
        for hull in &feature.convex_hulls {
            let Some(bounds) = hull.bounding_box() else {
                continue; // an empty source box never encloses anything
            };
            let mut b = Box2 {
                min_rt: bounds.min().rt,
                max_rt: bounds.max().rt,
                min_mz: bounds.min().mz,
                max_mz: bounds.max().mz,
            };
            if centroid_rt {
                b.min_rt = feature.rt;
                b.max_rt = feature.rt;
            }
            self.increase_box(&mut b);
            if b.encloses(rt, mz) {
                return true;
            }
        }
        false
    }

    /// Annotate a consensus map with peptide identifications, the source
    /// `annotate(ConsensusMap&, const PeptideIdentificationList&, const
    /// std::vector<ProteinIdentification>&, bool measure_from_subelements, bool
    /// annotate_ids_with_subelements, const PeakMap& spectra)`.
    ///
    /// `protein_ids` are appended to the map together with one identification
    /// data-processing record per run.
    ///
    /// **Isobaric (TMT/iTRAQ) data.** When a consensus feature carries an
    /// `id_scan_id` (MS3 TMT) or `scan_id` (MS2 TMT) meta value — the first
    /// consensus feature carrying either decides which — identifications are
    /// mapped by spectrum reference instead: identifications with hits and a
    /// `spectrum_reference` are indexed by the base name of their primary MS
    /// run path and that reference (the first one wins; a duplicate is a
    /// warning), and each consensus feature receives the identification whose
    /// file matches the base name of its first handle's column header and
    /// whose reference equals the feature's scan ID. When the first
    /// identification's spectrum reference is a plain integer, a scan ID
    /// without an exact match is retried by its scan number — see the note
    /// below. `measure_from_subelements` and charge checking are ignored with a
    /// warning.
    ///
    /// **Other data.** Each identification with hits is compared against
    /// every consensus feature — its centroid, or with
    /// `measure_from_subelements` each of its handles — with
    /// [`IDMapper::is_match`]; unless `ignore_charge`, the feature's (or
    /// handle's) charge must be among the hits' charges or zero. A matching
    /// identification is appended to every matching consensus feature (once
    /// each; with `annotate_ids_with_subelements` it carries the matching
    /// handle's `map_index`), and one matching none to the unassigned
    /// identifications.
    ///
    /// Unidentified precursors from `spectra` are annotated as empty
    /// identifications as for [`IDMapper::annotate_features`], matched
    /// against centroids or every matching handle (one copy per matching
    /// handle). Here a charge of zero on either side is not special for the
    /// feature: the precursor charge or zero must equal the feature's charge.
    ///
    /// # Notes
    ///
    /// In the source the scan-number fallback compiles its regular expression
    /// from `""` and asks `boost::regex::empty()`, which is false for that
    /// valid empty pattern, so `getRegExFromNativeID` is never consulted and
    /// the empty pattern — which has no group 1 — is used. Its tokens are
    /// empty strings, which do not convert, so
    /// `SpectrumNativeIDParser::extractScanNumber` throws `ParseError`. This
    /// port executes the same pattern through the crate's Boost-compatible
    /// regex and therefore returns [`Error::Parse`] in the same situation. The
    /// behaviour is recorded as a C++ defect candidate in the support document;
    /// an executed C++ oracle (Release build) confirms the `ParseError`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when an identification lacks RT or
    /// m/z, or — for isobaric data — when no identification carries a spectrum
    /// reference; [`Error::Parse`] as described in the note;
    /// [`Error::InvalidValue`] when a coordinate is not finite, a peptide mass
    /// cannot be computed, an isobaric consensus feature has no handles (the
    /// source dereferences its empty handle set), the protein runs' MS-run
    /// paths are inconsistent, or a limit is exceeded. The map is unchanged on
    /// error.
    pub fn annotate_consensus(
        &self,
        map: &mut ConsensusMap,
        ids: &[PeptideIdentification],
        protein_ids: &[ProteinIdentification],
        measure_from_subelements: bool,
        annotate_ids_with_subelements: bool,
        spectra: &MSExperiment,
    ) -> Result<MapMappingStatistics> {
        Self::check_ids(ids)?;
        check_spectra(spectra)?;
        for feature in &map.features {
            if !feature.rt.is_finite() || !feature.mz.is_finite() {
                return Err(Error::InvalidValue(
                    "IDMapper: consensus feature coordinates must be finite".into(),
                ));
            }
        }
        let mut stats = MapMappingStatistics::default();
        let new_processing = identification_data_processing(protein_ids)?;
        let mut work = 0u64;
        let max_work = self.limits.max_comparisons;
        let mut count_work = |n: usize| -> Result<()> {
            work = work.saturating_add(n as u64);
            if work > max_work {
                Err(limit("comparison count"))
            } else {
                Ok(())
            }
        };

        // (consensus index, identification) pushes, in source order.
        let mut pushes: Vec<(usize, PeptideIdentification)> = Vec::new();
        let mut unassigned_ids: Vec<usize> = Vec::new();
        match isobaric_type(map) {
            Some(scan_key) => {
                self.plan_isobaric(
                    map,
                    ids,
                    protein_ids,
                    measure_from_subelements,
                    scan_key,
                    &mut pushes,
                    &mut stats,
                )?;
            }
            None => {
                let mut assigned: BTreeMap<usize, usize> = BTreeMap::new();
                let mut mapping: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); map.features.len()];
                for (i, id) in ids.iter().enumerate() {
                    if id.hits.is_empty() {
                        continue;
                    }
                    let details = self.id_details(id, false)?;
                    let mut id_mapped = false;
                    for (cm_index, feature) in map.features.iter().enumerate() {
                        count_work(
                            details
                                .mz_values
                                .len()
                                .saturating_mul(feature.handles().len().max(1)),
                        )?;
                        for (i_mz, &mz_pep) in details.mz_values.iter().enumerate() {
                            let mut current_charges = Vec::new();
                            if !self.ignore_charge {
                                if details.mz_values.len() == 1 {
                                    current_charges.extend_from_slice(&details.charges);
                                } else {
                                    current_charges.push(details.charges[i_mz]);
                                }
                                current_charges.push(0);
                            }
                            let mut was_added = false;
                            if !measure_from_subelements {
                                if self.is_match(details.rt - feature.rt, mz_pep, feature.mz)
                                    && (self.ignore_charge
                                        || charge_list_contains(&current_charges, feature.charge))
                                {
                                    id_mapped = true;
                                    was_added = true;
                                    pushes.push((cm_index, id.clone()));
                                    *assigned.entry(i).or_default() += 1;
                                }
                            } else {
                                for handle in feature.handles() {
                                    if self.is_match(details.rt - handle.rt, mz_pep, handle.mz)
                                        && (self.ignore_charge
                                            || charge_list_contains(
                                                &current_charges,
                                                handle.charge,
                                            ))
                                    {
                                        id_mapped = true;
                                        was_added = true;
                                        if mapping[cm_index].insert(i) {
                                            let mut id_pep = id.clone();
                                            if annotate_ids_with_subelements {
                                                id_pep.metadata.insert(
                                                    "map_index".into(),
                                                    MetaValue::try_from(handle.map_index)?,
                                                );
                                            }
                                            pushes.push((cm_index, id_pep));
                                            *assigned.entry(i).or_default() += 1;
                                        }
                                        break;
                                    }
                                }
                            }
                            if pushes.len() > self.limits.max_assignments {
                                return Err(limit("assignment count"));
                            }
                            if was_added {
                                break;
                            }
                        }
                    }
                    if !id_mapped {
                        unassigned_ids.push(i);
                        stats.peptides_unassigned += 1;
                    }
                }
                for &count in assigned.values() {
                    if count == 1 {
                        stats.peptides_single += 1;
                    } else if count > 1 {
                        stats.peptides_multiple += 1;
                    }
                }
            }
        }

        // Unidentified precursors.
        let unidentified =
            map_precursors(spectra, ids, 0.001, 0.001, self.limits.max_comparisons)?.unidentified;
        let (run_identifier, add_unknown_run) = empty_run_identifier(
            !unidentified.is_empty(),
            map.protein_identifications
                .iter()
                .chain(protein_ids.iter())
                .next(),
            map.unassigned_peptide_identifications
                .first()
                .map(|id| id.identifier.as_str())
                .or_else(|| unassigned_ids.first().map(|&i| ids[i].identifier.as_str())),
        );
        let mut assigned_precursors: BTreeMap<usize, usize> = BTreeMap::new();
        let mut precursor_pushes: Vec<(usize, PeptideIdentification)> = Vec::new();
        for &spectrum_index in &unidentified {
            let spectrum = &spectra.spectra[spectrum_index];
            let mut precursor_mapped = false;
            for precursor in &spectrum.precursors {
                let (mz_p, z_p, rt_value) = (precursor.mz, precursor.charge, spectrum.rt);
                let mut empty_id = precursor_empty_id(
                    spectrum_index,
                    rt_value,
                    mz_p,
                    &spectrum.native_id,
                    &run_identifier,
                )?;
                for (cm_index, feature) in map.features.iter().enumerate() {
                    count_work(feature.handles().len().max(1))?;
                    let current_charges: Vec<i32> = if self.ignore_charge {
                        Vec::new()
                    } else {
                        vec![z_p, 0]
                    };
                    if !measure_from_subelements {
                        if self.is_match(rt_value - feature.rt, mz_p, feature.mz)
                            && (self.ignore_charge
                                || charge_list_contains(&current_charges, feature.charge))
                        {
                            precursor_pushes.push((cm_index, empty_id.clone()));
                            *assigned_precursors.entry(spectrum_index).or_default() += 1;
                            precursor_mapped = true;
                        }
                    } else {
                        for handle in feature.handles() {
                            if self.is_match(rt_value - handle.rt, mz_p, handle.mz)
                                && (self.ignore_charge
                                    || charge_list_contains(&current_charges, handle.charge))
                            {
                                if annotate_ids_with_subelements {
                                    // The source mutates the one identification
                                    // it reuses for every later push.
                                    empty_id.metadata.insert(
                                        "map_index".into(),
                                        MetaValue::try_from(handle.map_index)?,
                                    );
                                }
                                precursor_pushes.push((cm_index, empty_id.clone()));
                                *assigned_precursors.entry(spectrum_index).or_default() += 1;
                                precursor_mapped = true;
                            }
                        }
                    }
                    if pushes.len() + precursor_pushes.len() > self.limits.max_assignments {
                        return Err(limit("assignment count"));
                    }
                }
            }
            if !precursor_mapped {
                stats.precursors_unassigned += 1;
            }
        }
        for &count in assigned_precursors.values() {
            if count == 1 {
                stats.precursors_single += 1;
            } else if count > 1 {
                stats.precursors_multiple += 1;
            }
        }

        // Commit.
        map.protein_identifications
            .extend(protein_ids.iter().cloned());
        map.data_processing.extend(new_processing);
        for (cm_index, id) in pushes {
            map.features[cm_index].peptide_identifications.push(id);
        }
        for i in unassigned_ids {
            map.unassigned_peptide_identifications.push(ids[i].clone());
        }
        if add_unknown_run {
            map.protein_identifications.push(unknown_run());
        }
        for (cm_index, id) in precursor_pushes {
            map.features[cm_index].peptide_identifications.push(id);
        }
        Ok(stats)
    }

    /// Plan the isobaric (TMT/iTRAQ) branch of
    /// [`IDMapper::annotate_consensus`].
    #[allow(clippy::too_many_arguments)]
    fn plan_isobaric(
        &self,
        map: &ConsensusMap,
        ids: &[PeptideIdentification],
        protein_ids: &[ProteinIdentification],
        measure_from_subelements: bool,
        scan_key: &str,
        pushes: &mut Vec<(usize, PeptideIdentification)>,
        stats: &mut MapMappingStatistics,
    ) -> Result<()> {
        let run_mapping = IdentifierMSRunMapper::from_runs(protein_ids)?;
        let mut look_for_scan_numbers = false;
        let mut by_file: BTreeMap<String, HashMap<String, usize>> = BTreeMap::new();
        let mut has_references = false;
        for (i, id) in ids.iter().enumerate() {
            if id.hits.is_empty() {
                continue;
            }
            let spectrum_file = basename(&run_mapping.primary_ms_run_path(id)?).to_owned();
            let reference = match id.metadata.get(SPECTRUM_REFERENCE) {
                None => String::new(),
                Some(value) => value
                    .as_str()
                    .map_err(|_| {
                        Error::InvalidValue(
                            "IDMapper: 'spectrum_reference' meta value is not a string".into(),
                        )
                    })?
                    .to_owned(),
            };
            if reference.is_empty() {
                continue;
            }
            if !look_for_scan_numbers {
                // Source: always the *first* identification's reference.
                look_for_scan_numbers = source_to_int64(&ids[0].spectrum_reference()).is_some();
            }
            let inner = by_file.entry(spectrum_file).or_default();
            match inner.entry(reference) {
                std::collections::hash_map::Entry::Occupied(entry) => {
                    stats.warnings.push(format!(
                        "Duplicate spectrum reference detected: {}",
                        entry.key()
                    ));
                }
                std::collections::hash_map::Entry::Vacant(entry) => {
                    entry.insert(i);
                }
            }
            has_references = true;
        }
        if !has_references {
            return Err(Error::MissingInformation(
                "No spectrum references in ID file used in id mapping TMT/iTRAQ data.".into(),
            ));
        }
        if measure_from_subelements {
            stats.warnings.push("IDMapper is configured to measure from subelements. Because the data looks like TMT/iTRAQ this option will be ignored.".into());
        }
        if !self.ignore_charge {
            stats.warnings.push("IDMapper is configured to validate charges. Because the data looks like TMT/iTRAQ this option will be ignored.".into());
        }
        // See the note at annotate_consensus: the source's pattern stays "".
        let scan_regex = BoostRegex::new("")?;
        for (cm_index, feature) in map.features.iter().enumerate() {
            let first = feature.handles().first().ok_or_else(|| {
                Error::InvalidValue("IDMapper: isobaric consensus feature without handles".into())
            })?;
            let filename = map
                .column_headers
                .get(&first.map_index)
                .map_or("", |header| basename(&header.filename))
                .to_owned();
            let scan_id = feature
                .metadata
                .get(scan_key)
                .map(ToString::to_string)
                .unwrap_or_default();
            if scan_id.is_empty() {
                stats.warnings.push(format!("ConsensusMap for TMT/iTRAQ experiment is missing the scan identifier meta value '{scan_id}'"));
                continue;
            }
            match by_file.get(&filename) {
                Some(run) => {
                    if let Some(&i) = run.get(&scan_id) {
                        pushes.push((cm_index, ids[i].clone()));
                        stats.peptides_single += 1;
                    } else if look_for_scan_numbers {
                        let scan = extract_scan_number(&scan_id, &scan_regex)?;
                        if let Some(&i) = run.get(&scan.to_string()) {
                            pushes.push((cm_index, ids[i].clone()));
                            stats.peptides_single += 1;
                        }
                    }
                }
                None => stats.warnings.push(format!(
                    "ConsensusMap for TMT/iTRAQ experiment contains scan identifier '{scan_id}' quantified in file '{filename}' but there is no matching identification."
                )),
            }
            if pushes.len() > self.limits.max_assignments {
                return Err(limit("assignment count"));
            }
        }
        Ok(())
    }
}

/// Source `checkTMTType`: the scan-ID meta key of the first consensus feature
/// carrying `id_scan_id` (MS3 TMT) or `scan_id` (MS2 TMT), else `None`.
fn isobaric_type(map: &ConsensusMap) -> Option<&'static str> {
    map.features.iter().find_map(|feature: &ConsensusFeature| {
        if feature.metadata.contains_key("id_scan_id") {
            Some("id_scan_id")
        } else if feature.metadata.contains_key("scan_id") {
            Some("scan_id")
        } else {
            None
        }
    })
}

fn check_spectra(spectra: &MSExperiment) -> Result<()> {
    for spectrum in &spectra.spectra {
        if spectrum.precursors.is_empty() {
            continue;
        }
        if !spectrum.rt.is_finite() || spectrum.precursors.iter().any(|p| !p.mz.is_finite()) {
            return Err(Error::InvalidValue(
                "IDMapper: spectrum RT and precursor m/z must be finite".into(),
            ));
        }
    }
    Ok(())
}

fn map_precursors(
    spectra: &MSExperiment,
    ids: &[PeptideIdentification],
    mz_tol: f64,
    rt_tol: f64,
    max_comparisons: u64,
) -> Result<PeptideIdentificationListState> {
    let precursor_count: u64 = spectra
        .spectra
        .iter()
        .map(|s| s.precursors.len() as u64)
        .sum();
    if precursor_count.saturating_mul(ids.len() as u64) > max_comparisons {
        return Err(limit("comparison count"));
    }
    let mut state = PeptideIdentificationListState::default();
    for (index, spectrum) in spectra.spectra.iter().enumerate() {
        if spectrum.precursors.is_empty() {
            state.no_precursors.push(index);
            continue;
        }
        let rt_s = spectrum.rt;
        let identified = spectrum.precursors.iter().any(|precursor| {
            ids.iter().any(|id| {
                !id.hits.is_empty()
                    && (id.mz.unwrap_or(f64::NAN) - precursor.mz).abs() < mz_tol
                    && (rt_s - id.rt.unwrap_or(f64::NAN)).abs() < rt_tol
            })
        });
        if identified {
            state.identified.push(index);
        } else {
            state.unidentified.push(index);
        }
    }
    Ok(state)
}

/// The identifier for empty precursor identifications and whether a new
/// unknown run must be added, as the source decides it.
fn empty_run_identifier(
    needed: bool,
    first_run: Option<&ProteinIdentification>,
    first_unassigned: Option<&str>,
) -> (String, bool) {
    if !needed {
        return (String::new(), false);
    }
    if let Some(run) = first_run {
        return (run.identifier.clone(), false);
    }
    if let Some(identifier) = first_unassigned {
        return (identifier.to_owned(), false);
    }
    (UNKNOWN_SEARCH_RUN_IDENTIFIER.to_owned(), true)
}

fn unknown_run() -> ProteinIdentification {
    ProteinIdentification {
        identifier: UNKNOWN_SEARCH_RUN_IDENTIFIER.into(),
        date_time: Some(DateTime::now().iso_string()),
        ..ProteinIdentification::default()
    }
}

fn precursor_empty_id(
    spectrum_index: usize,
    rt: f64,
    mz: f64,
    native_id: &str,
    identifier: &str,
) -> Result<PeptideIdentification> {
    let mut id = PeptideIdentification {
        rt: Some(rt),
        mz: Some(mz),
        identifier: identifier.to_owned(),
        ..PeptideIdentification::default()
    };
    id.metadata.insert(
        "spectrum_index".into(),
        MetaValue::try_from(spectrum_index as u64)?,
    );
    if !native_id.is_empty() {
        id.metadata
            .insert(SPECTRUM_REFERENCE.into(), MetaValue::from(native_id));
    }
    Ok(id)
}

fn identification_data_processing(
    protein_ids: &[ProteinIdentification],
) -> Result<Vec<DataProcessing>> {
    let mut records = Vec::with_capacity(protein_ids.len());
    for run in protein_ids {
        let mut dp = DataProcessing::default();
        dp.software.name.clone_from(&run.search_engine);
        dp.software.version.clone_from(&run.search_engine_version);
        dp.completion_time = run
            .date_time
            .as_deref()
            .and_then(|text| DateTime::parse(text).ok());
        dp.actions.insert(ProcessingAction::Identification);
        let parameters = &run.search_parameters;
        if !parameters.database.is_empty() {
            dp.metadata
                .insert("parameter: db".into(), parameters.database.as_str().into());
        }
        if !parameters.database_version.is_empty() {
            dp.metadata.insert(
                "parameter: db_version".into(),
                parameters.database_version.as_str().into(),
            );
        }
        records.push(dp);
    }
    Ok(records)
}

fn check_mass_type<'a>(
    processing: impl Iterator<Item = &'a DataProcessing>,
    warnings: &mut Vec<String>,
) -> bool {
    let mut use_avg_mass = false;
    let mut before = String::new();
    for record in processing {
        if record.software.name != "FeatureFinder" {
            continue;
        }
        let reported = record
            .metadata
            .get("parameter: algorithm:feature:reported_mz")
            .map(ToString::to_string)
            .unwrap_or_default();
        if reported.is_empty() {
            continue;
        }
        if !before.is_empty() && reported != before {
            warnings.push("The m/z values reported for features in the input seem to be of different types (e.g. monoisotopic/average). They will all be compared against monoisotopic peptide masses, but the mapping results may not be meaningful in the end.".into());
            return false;
        }
        if reported == "average" {
            use_avg_mass = true;
        } else if reported == "maximum" {
            warnings.push("For features, m/z values from the highest mass traces are reported. This type of m/z value is not available for peptides, so the comparison has to be done using average peptide masses.".into());
            use_avg_mass = true;
        }
        before = reported;
    }
    use_avg_mass
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_int64_follows_string_utils() {
        assert_eq!(source_to_int64(" +42 "), Some(42));
        assert_eq!(source_to_int64("-7"), Some(-7));
        assert_eq!(source_to_int64("12a"), None);
        assert_eq!(source_to_int64(""), None);
        assert_eq!(source_to_int64("+-3"), Some(-3));
        assert_eq!(source_to_int64("controllerType=0 scan=5"), None);
    }

    #[test]
    fn empty_pattern_scan_extraction_fails_like_the_source() {
        let regex = BoostRegex::new("").unwrap();
        assert!(matches!(
            extract_scan_number("controllerType=0 controllerNumber=1 scan=5", &regex),
            Err(Error::Parse { .. })
        ));
        let scan = BoostRegex::new(r"scan=(?<GROUP>\d+)").unwrap();
        assert_eq!(
            extract_scan_number("frame=1 scan=17 scan=23", &scan).unwrap(),
            23
        );
    }

    #[test]
    fn source_empty_needs_every_default() {
        let mut id = PeptideIdentification::default();
        assert!(source_empty(&id));
        id.score_type = "q".into();
        assert!(!source_empty(&id));
    }
}
