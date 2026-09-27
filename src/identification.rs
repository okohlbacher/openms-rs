// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Owned peptide/protein identification records, scores and protein evidence.
//! Positions are zero-based and inclusive, matching the OpenMS metadata model.
//! Missing values use Option rather than numeric sentinels or NaN.

pub mod graph;
mod protein_run;
mod run_mapping;
pub use run_mapping::IdentifierMSRunMapper;
mod peptide_identification;
pub use peptide_identification::{ConsensusIdLocation, ConsensusUidMap, MAX_UID_IDENTIFICATIONS};
pub mod usi;

use crate::chemistry::{AASequence, SequenceModification};
use crate::comparison::Tolerance;
use crate::kernel::DataArray;
use crate::metadata::{MetaInfo, validate_meta};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::ops::RangeInclusive;

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn finite(value: f64, name: &str) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(bad(&format!("{name} must be finite")))
    }
}
fn array_not_found(kind: &str, name: &str) -> Error {
    Error::InvalidValue(format!("element not found: {kind}: {name}"))
}
fn text_value(metadata: &MetaInfo, key: &str) -> String {
    metadata
        .get(key)
        .map(ToString::to_string)
        .unwrap_or_default()
}

/// Flanking amino acid or one of OpenMS's X/[ /] markers.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum FlankingResidue {
    #[default]
    Unknown,
    NTerminus,
    CTerminus,
    Residue(char),
}
impl Ord for FlankingResidue {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Match the source character ordering for all valid markers. The final
        // key keeps total ordering consistent with Eq for invalid public variants.
        let kind = |v: &Self| match v {
            Self::Unknown => 0,
            Self::NTerminus => 1,
            Self::CTerminus => 2,
            Self::Residue(_) => 3,
        };
        self.code()
            .cmp(&other.code())
            .then(kind(self).cmp(&kind(other)))
    }
}
impl PartialOrd for FlankingResidue {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl FlankingResidue {
    pub fn from_code(code: char) -> Result<Self> {
        match code {
            'X' => Ok(Self::Unknown),
            '[' => Ok(Self::NTerminus),
            ']' => Ok(Self::CTerminus),
            c if c.is_ascii_uppercase() => Ok(Self::Residue(c)),
            _ => Err(bad(
                "flanking residue must be an uppercase residue or X/[/] marker",
            )),
        }
    }
    pub fn code(self) -> char {
        match self {
            Self::Unknown => 'X',
            Self::NTerminus => '[',
            Self::CTerminus => ']',
            Self::Residue(c) => c,
        }
    }
    fn validate(self) -> Result<()> {
        if let Self::Residue(c) = self {
            if !c.is_ascii_uppercase() || c == 'X' {
                return Err(bad("invalid explicit flanking residue"));
            }
        }
        Ok(())
    }
}

/// One mapping of a peptide to a protein accession. End is inclusive.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PeptideEvidence {
    pub protein_accession: String,
    pub start: Option<usize>,
    pub end: Option<usize>,
    pub aa_before: FlankingResidue,
    pub aa_after: FlankingResidue,
}
impl PeptideEvidence {
    pub fn new(accession: impl Into<String>, positions: RangeInclusive<usize>) -> Result<Self> {
        let value = Self {
            protein_accession: accession.into(),
            start: Some(*positions.start()),
            end: Some(*positions.end()),
            ..Default::default()
        };
        value.validate()?;
        Ok(value)
    }
    pub fn validate(&self) -> Result<()> {
        if let (Some(start), Some(end)) = (self.start, self.end) {
            if start > end {
                return Err(bad("peptide evidence start exceeds inclusive end"));
            }
        }
        self.aa_before.validate()?;
        self.aa_after.validate()
    }
    /// Includes the valid one-residue mapping 0..=0, unlike the C++ helper.
    pub fn has_valid_limits(&self) -> bool {
        matches!((self.start,self.end),(Some(start),Some(end)) if start<=end)
    }
    pub fn positions(&self) -> Result<RangeInclusive<usize>> {
        self.validate()?;
        match (self.start, self.end) {
            (Some(start), Some(end)) => Ok(start..=end),
            _ => Err(bad("peptide evidence has unknown start or end")),
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum TargetDecoyType {
    #[default]
    Unknown,
    Target,
    Decoy,
    TargetAndDecoy,
}
fn target_decoy(metadata: &MetaInfo, allow_mixed: bool) -> Result<TargetDecoyType> {
    let Some(value) = metadata.get("target_decoy") else {
        return Ok(TargetDecoyType::Unknown);
    };
    match value.to_string().to_ascii_lowercase().as_str() {
        "target" => Ok(TargetDecoyType::Target),
        "decoy" => Ok(TargetDecoyType::Decoy),
        "target+decoy" if allow_mixed => Ok(TargetDecoyType::TargetAndDecoy),
        _ => Err(bad("invalid target_decoy metadata value")),
    }
}
fn set_target_decoy(metadata: &mut MetaInfo, value: TargetDecoyType) {
    match value {
        TargetDecoyType::Unknown => {
            metadata.remove("target_decoy");
        }
        value => {
            metadata.insert(
                "target_decoy".into(),
                match value {
                    TargetDecoyType::Target => "target",
                    TargetDecoyType::Decoy => "decoy",
                    _ => "target+decoy",
                }
                .into(),
            );
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PeakAnnotation {
    pub mz: f64,
    pub intensity: f64,
    pub charge: i32,
    pub annotation: String,
}
impl PeakAnnotation {
    pub fn validate(&self) -> Result<()> {
        finite(self.mz, "annotation m/z")?;
        finite(self.intensity, "annotation intensity")
    }
    /// Sort by m/z, charge, annotation and intensity, preserving equal entries.
    pub fn sort(annotations: &mut [Self]) -> Result<()> {
        for annotation in annotations.iter() {
            annotation.validate()?;
        }
        annotations.sort_by(|a, b| {
            a.mz.partial_cmp(&b.mz)
                .unwrap()
                .then(a.charge.cmp(&b.charge))
                .then(a.annotation.cmp(&b.annotation))
                .then(a.intensity.partial_cmp(&b.intensity).unwrap())
        });
        Ok(())
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AnalysisResult {
    pub score_type: String,
    pub higher_is_better: bool,
    pub main_score: f64,
    pub sub_scores: BTreeMap<String, f64>,
}
impl AnalysisResult {
    pub fn validate(&self) -> Result<()> {
        finite(self.main_score, "analysis score")?;
        for &score in self.sub_scores.values() {
            finite(score, "analysis sub-score")?;
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PeptideHit {
    pub sequence: AASequence,
    pub score: f64,
    pub rank: u32,
    pub charge: i32,
    pub evidences: Vec<PeptideEvidence>,
    pub peak_annotations: Vec<PeakAnnotation>,
    pub analysis_results: Vec<AnalysisResult>,
    pub metadata: MetaInfo,
}
impl PeptideHit {
    pub fn new(score: f64, rank: u32, charge: i32, sequence: AASequence) -> Result<Self> {
        let value = Self {
            score,
            rank,
            charge,
            sequence,
            ..Default::default()
        };
        value.validate()?;
        Ok(value)
    }
    pub fn validate(&self) -> Result<()> {
        finite(self.score, "peptide score")?;
        validate_meta(&self.metadata)?;
        for evidence in &self.evidences {
            evidence.validate()?;
        }
        for annotation in &self.peak_annotations {
            annotation.validate()?;
        }
        for result in &self.analysis_results {
            result.validate()?;
        }
        self.target_decoy_type()?;
        Ok(())
    }
    pub fn protein_accessions(&self) -> BTreeSet<&str> {
        self.evidences
            .iter()
            .map(|e| e.protein_accession.as_str())
            .filter(|a| !a.is_empty())
            .collect()
    }
    /// Complete annotated sequence and charge, independent of display aliases.
    pub fn identity_key(&self) -> (AASequence, i32) {
        (self.sequence.clone(), self.charge)
    }
    pub fn same_sequence_and_charge(&self, other: &Self) -> bool {
        self.sequence == other.sequence && self.charge == other.charge
    }
    pub fn target_decoy_type(&self) -> Result<TargetDecoyType> {
        target_decoy(&self.metadata, true)
    }
    pub fn is_decoy(&self) -> Result<bool> {
        Ok(self.target_decoy_type()? == TargetDecoyType::Decoy)
    }
    pub fn set_target_decoy_type(&mut self, value: TargetDecoyType) {
        set_target_decoy(&mut self.metadata, value);
    }
}

/// The candidates (spectrum matches) identified for one precursor spectrum
/// (source `METADATA/PeptideIdentification.h`).
///
/// Typically the result of searching one MS/MS spectrum against a sequence
/// database or spectral library: a list of [`PeptideHit`]s, each a potential
/// match. An identification belongs to a parent [`ProteinIdentification`]
/// run, which holds the search parameters, database and identified proteins;
/// the two are linked by the shared [`identifier`](Self::identifier) string,
/// and many identifications (one per spectrum) can share one run.
///
/// The precursor's retention time and m/z map the identification back to
/// experimental data - spectra, features or consensus features (the source's
/// `IDMapper`). Missing coordinates are `None` where the source stores a quiet
/// NaN, so `getRT`/`setRT`/`hasRT` and their m/z counterparts are the `rt` and
/// `mz` fields with `is_some`. The score type, orientation and an optional
/// significance threshold describe the hits' scores.
///
/// The source class is also available as [`SpectrumIdentification`], the name
/// its header calls "better" and deprecates this one in favour of. Equality is
/// field-wise; the source's NaN-equals-NaN rule for unset coordinates is
/// `None == None` here. Hashing is consistent with equality (source
/// `std::hash<PeptideIdentification>`).
#[derive(Clone, Debug, PartialEq)]
pub struct PeptideIdentification {
    /// Links this identification to its [`ProteinIdentification`] run (source
    /// `getIdentifier`/`setIdentifier`).
    pub identifier: String,
    /// The candidate hits (source `getHits`/`setHits`/`insertHit`).
    pub hits: Vec<PeptideHit>,
    /// Name of the score in [`PeptideHit::score`], for example `Mascot`,
    /// `e-value` or `q-value` (source `getScoreType`/`setScoreType`).
    pub score_type: String,
    /// Score orientation (source `isHigherScoreBetter`/`setHigherScoreBetter`);
    /// `true` by default.
    pub higher_score_better: bool,
    /// Peptide significance threshold (source
    /// `getSignificanceThreshold`/`setSignificanceThreshold`). The source keeps
    /// it as the `significance_threshold` metadata value and removes that value
    /// when set to 0; here it is a field, 0 by default, and the file formats
    /// write the source attribute.
    pub significance_threshold: f64,
    /// Retention time of the identified MS2 spectrum, if known.
    pub rt: Option<f64>,
    /// Precursor m/z of the identified MS2 spectrum, if known.
    pub mz: Option<f64>,
    /// The source `MetaInfoInterface` payload, including `spectrum_reference`
    /// and `experiment_label`.
    pub metadata: MetaInfo,
}

/// Source `using SpectrumIdentification = PeptideIdentification`, the name the
/// source header marks as the future default.
pub type SpectrumIdentification = PeptideIdentification;

impl Default for PeptideIdentification {
    fn default() -> Self {
        Self {
            identifier: String::new(),
            hits: Vec::new(),
            score_type: String::new(),
            higher_score_better: true,
            significance_threshold: 0.0,
            rt: None,
            mz: None,
            metadata: MetaInfo::new(),
        }
    }
}
impl PeptideIdentification {
    /// Source default constructor: no hits, no run, higher scores better,
    /// threshold 0 and unknown RT and m/z.
    pub fn new() -> Self {
        Self::default()
    }
    /// Number of candidate hits.
    pub fn len(&self) -> usize {
        self.hits.len()
    }
    /// Whether there are no candidate hits. The source `empty()`, which also
    /// requires the default run, score type, orientation and threshold, is
    /// [`PeptideIdentification::source_empty`].
    pub fn is_empty(&self) -> bool {
        self.hits.is_empty()
    }
    /// Check the invariants the source leaves unchecked: a finite threshold,
    /// finite coordinates when present, valid metadata and valid hits.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] naming the first violated invariant.
    pub fn validate(&self) -> Result<()> {
        finite(
            self.significance_threshold,
            "peptide significance threshold",
        )?;
        if let Some(rt) = self.rt {
            finite(rt, "peptide RT")?;
        }
        if let Some(mz) = self.mz {
            finite(mz, "peptide m/z")?;
        }
        validate_meta(&self.metadata)?;
        for hit in &self.hits {
            hit.validate()?;
        }
        Ok(())
    }
    /// Source `sort()`: order the hits so that the best-scoring hit comes
    /// first, honouring [`Self::higher_score_better`].
    ///
    /// The sort is stable, as the source's `std::stable_sort` with
    /// [`PeptideIdentification::score_comparator`]: equal scores keep their
    /// input order.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] from [`PeptideIdentification::validate`], before
    /// anything moves. The source sorts NaN scores with a comparator that is
    /// then not a strict weak ordering; this refuses them.
    pub fn sort(&mut self) -> Result<()> {
        self.validate()?;
        self.hits
            .sort_by(|a, b| score_order(a.score, b.score, self.higher_score_better));
        Ok(())
    }
    /// Native convenience: best-score hit without changing stored candidate order.
    pub fn best_hit(&self) -> Result<Option<&PeptideHit>> {
        self.validate()?;
        Ok(self
            .hits
            .iter()
            .min_by(|a, b| score_order(a.score, b.score, self.higher_score_better)))
    }
    /// Source static `getReferencingHits(getHits(), accessions)`: this
    /// identification's hits with at least one peptide evidence whose protein
    /// accession is in `accessions`, in hit order.
    /// [`PeptideIdentification::hits_referencing`] takes any hit slice.
    pub fn referencing_hits<'a>(&'a self, accessions: &BTreeSet<String>) -> Vec<&'a PeptideHit> {
        Self::hits_referencing(&self.hits, accessions)
    }
    /// Source `getSpectrumReference()`: the `spectrum_reference` metadata value
    /// as text, or `""` when absent. Currently this should almost always be the
    /// full native vendor ID.
    ///
    /// Like the source, a non-string value is stringified leniently rather than
    /// refused (an integer scan index from idparquet or mzIdentML reads as
    /// `17`), in the source's text form (`1.5`, `1.23456789125e08`, `[a, b]`).
    pub fn spectrum_reference(&self) -> String {
        source_text_value(&self.metadata, "spectrum_reference")
    }
    /// Source `setSpectrumReference`: store `spectrum_reference` as a string
    /// value. An empty reference is stored too, as in the source.
    pub fn set_spectrum_reference(&mut self, value: impl Into<String>) {
        self.metadata
            .insert("spectrum_reference".into(), value.into().into());
    }
    /// Source `getExperimentLabel()`: the `experiment_label` metadata value as
    /// text (mostly written for pepXML, where every identification can belong
    /// to a different experiment), or `""` when absent. Stringified leniently
    /// as [`PeptideIdentification::spectrum_reference`] is.
    pub fn experiment_label(&self) -> String {
        source_text_value(&self.metadata, "experiment_label")
    }
    /// Source `setExperimentLabel`: store a non-empty label.
    ///
    /// An empty label is **ignored**, exactly as in the source, which does not
    /// store the empty default: a label set earlier stays in place (executed
    /// Release oracle `B8b-peptide-id`, case `labels`). Remove the
    /// `experiment_label` metadata entry to clear it.
    pub fn set_experiment_label(&mut self, value: impl Into<String>) {
        let value = value.into();
        if !value.is_empty() {
            self.metadata
                .insert("experiment_label".into(), value.into());
        }
    }
}
/// A metadata value as the source's lenient `StringUtils::toStr(DataValue)`
/// writes it, or `""` when absent. The text form only fails beyond the list
/// byte ceiling, where the plain display form is used instead.
fn source_text_value(metadata: &MetaInfo, key: &str) -> String {
    use crate::data_structures::list::ListFormat;
    metadata.get(key).map_or_else(String::new, |value| {
        value
            .to_list_text()
            .map_or_else(|_| value.to_string(), std::borrow::Cow::into_owned)
    })
}
fn score_order(a: f64, b: f64, higher: bool) -> std::cmp::Ordering {
    if higher {
        b.partial_cmp(&a).unwrap()
    } else {
        a.partial_cmp(&b).unwrap()
    }
}

/// An observed modification at one zero-based protein residue position.
/// Owns anonymous mass tags and shared handles to immutable named chemistry.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ProteinModification {
    pub position: usize,
    pub modification: SequenceModification,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ProteinHit {
    pub score: f64,
    pub rank: u32,
    pub accession: String,
    /// Raw protein sequence, permitting ambiguous codes that lack peptide mass.
    pub sequence: String,
    /// Percentage in 0..=100; None corresponds to C++ COVERAGE_UNKNOWN.
    pub coverage: Option<f64>,
    pub modifications: Vec<ProteinModification>,
    pub metadata: MetaInfo,
}
impl ProteinHit {
    pub fn new(
        score: f64,
        rank: u32,
        accession: impl Into<String>,
        sequence: impl Into<String>,
    ) -> Result<Self> {
        let value = Self {
            score,
            rank,
            accession: accession.into().trim().into(),
            sequence: sequence.into().trim().into(),
            ..Default::default()
        };
        value.validate()?;
        Ok(value)
    }
    pub fn validate(&self) -> Result<()> {
        finite(self.score, "protein score")?;
        if !self
            .sequence
            .bytes()
            .all(|c| c.is_ascii_alphabetic() || matches!(c, b'*' | b'-' | b'.'))
        {
            return Err(bad("protein sequence must contain ASCII residue codes"));
        }
        if let Some(coverage) = self.coverage {
            if !coverage.is_finite() || !(0.0..=100.0).contains(&coverage) {
                return Err(bad("protein coverage must be a percentage"));
            }
        }
        for modification in &self.modifications {
            if !self.sequence.is_empty() && modification.position >= self.sequence.len() {
                return Err(bad("protein modification position outside sequence"));
            }
        }
        validate_meta(&self.metadata)?;
        self.target_decoy_type()?;
        Ok(())
    }
    pub fn description(&self) -> String {
        text_value(&self.metadata, "Description")
    }
    pub fn set_description(&mut self, value: impl Into<String>) {
        self.metadata
            .insert("Description".into(), value.into().into());
    }
    pub fn target_decoy_type(&self) -> Result<TargetDecoyType> {
        target_decoy(&self.metadata, false)
    }
    pub fn is_decoy(&self) -> Result<bool> {
        Ok(self.target_decoy_type()? == TargetDecoyType::Decoy)
    }
    pub fn set_target_decoy_type(&mut self, value: TargetDecoyType) -> Result<()> {
        if value == TargetDecoyType::TargetAndDecoy {
            return Err(bad("protein target/decoy status cannot be mixed"));
        }
        set_target_decoy(&mut self.metadata, value);
        Ok(())
    }
}

/// Protein grouping with quantitative sample arrays retained as group state.
/// Array lengths describe samples, not the number of protein accessions.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ProteinGroup {
    pub probability: f64,
    pub accessions: Vec<String>,
    pub float_data_arrays: Vec<DataArray<f32>>,
    pub integer_data_arrays: Vec<DataArray<i32>>,
    pub string_data_arrays: Vec<DataArray<String>>,
}
impl ProteinGroup {
    pub fn validate(&self) -> Result<()> {
        finite(self.probability, "protein group probability")?;
        if self
            .float_data_arrays
            .iter()
            .any(|a| a.data.iter().any(|v| !v.is_finite()))
        {
            return Err(bad("protein group quantities must be finite"));
        }
        Ok(())
    }
    /// The source `ProteinGroup::operator<`, a "best first" order: a higher
    /// probability is less, then a group with fewer accessions, then the
    /// lexically smaller accession list. The float comparisons are the source's
    /// own, so a NaN probability is neither greater nor smaller than any other
    /// and the comparison falls through to the accessions, exactly as in C++.
    /// Sample arrays take no part in the order.
    pub fn source_less(&self, other: &Self) -> bool {
        if self.probability > other.probability {
            return true;
        }
        if self.probability < other.probability {
            return false;
        }
        match self.accessions.len().cmp(&other.accessions.len()) {
            std::cmp::Ordering::Less => true,
            std::cmp::Ordering::Greater => false,
            std::cmp::Ordering::Equal => self.accessions < other.accessions,
        }
    }
    /// The first float array named `name`, the source `getFloatDataArrayByName`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when no float array has that name, where the
    /// source throws `Exception::ElementNotFound`.
    pub fn float_data_array_by_name(&self, name: &str) -> Result<&DataArray<f32>> {
        self.float_data_arrays
            .iter()
            .find(|array| array.name == name)
            .ok_or_else(|| array_not_found("FloatDataArray", name))
    }
    /// Mutable [`ProteinGroup::float_data_array_by_name`], the non-`const` overload.
    ///
    /// # Errors
    ///
    /// As [`ProteinGroup::float_data_array_by_name`].
    pub fn float_data_array_by_name_mut(&mut self, name: &str) -> Result<&mut DataArray<f32>> {
        self.float_data_arrays
            .iter_mut()
            .find(|array| array.name == name)
            .ok_or_else(|| array_not_found("FloatDataArray", name))
    }
    /// The first integer array named `name`, the source `getIntegerDataArrayByName`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when no integer array has that name, where the
    /// source throws `Exception::ElementNotFound`.
    pub fn integer_data_array_by_name(&self, name: &str) -> Result<&DataArray<i32>> {
        self.integer_data_arrays
            .iter()
            .find(|array| array.name == name)
            .ok_or_else(|| array_not_found("IntegerDataArray", name))
    }
    /// Mutable [`ProteinGroup::integer_data_array_by_name`], the non-`const` overload.
    ///
    /// # Errors
    ///
    /// As [`ProteinGroup::integer_data_array_by_name`].
    pub fn integer_data_array_by_name_mut(&mut self, name: &str) -> Result<&mut DataArray<i32>> {
        self.integer_data_arrays
            .iter_mut()
            .find(|array| array.name == name)
            .ok_or_else(|| array_not_found("IntegerDataArray", name))
    }
    /// The first string array named `name`, the source `getStringDataArrayByName`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when no string array has that name, where the
    /// source throws `Exception::ElementNotFound`.
    pub fn string_data_array_by_name(&self, name: &str) -> Result<&DataArray<String>> {
        self.string_data_arrays
            .iter()
            .find(|array| array.name == name)
            .ok_or_else(|| array_not_found("StringDataArray", name))
    }
    /// Mutable [`ProteinGroup::string_data_array_by_name`], the non-`const` overload.
    ///
    /// # Errors
    ///
    /// As [`ProteinGroup::string_data_array_by_name`].
    pub fn string_data_array_by_name_mut(&mut self, name: &str) -> Result<&mut DataArray<String>> {
        self.string_data_arrays
            .iter_mut()
            .find(|array| array.name == name)
            .ok_or_else(|| array_not_found("StringDataArray", name))
    }
    /// Source order: descending probability, fewer accessions, lexical list.
    pub fn sort(groups: &mut [Self]) -> Result<()> {
        for group in groups.iter() {
            group.validate()?;
        }
        groups.sort_by(|a, b| {
            b.probability
                .partial_cmp(&a.probability)
                .unwrap()
                .then(a.accessions.len().cmp(&b.accessions.len()))
                .then(a.accessions.cmp(&b.accessions))
        });
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PeakMassType {
    #[default]
    Monoisotopic,
    Average,
}
impl PeakMassType {
    /// The source `NamesOfPeakMassType`, indexed by the enum order.
    pub const NAMES: [&'static str; 2] = ["Monoisotopic", "Average"];
    /// This mass type's entry of [`PeakMassType::NAMES`].
    pub fn name(self) -> &'static str {
        Self::NAMES[self as usize]
    }
    /// The source `getAllNamesOfPeakMassType`, which copies the name array into
    /// a `StringList`; [`PeakMassType::NAMES`] is the same list unallocated.
    pub fn all_names() -> Vec<String> {
        Self::NAMES.iter().map(|name| (*name).to_string()).collect()
    }
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum EnzymeTermSpecificity {
    #[default]
    Unknown,
    Full,
    Semi,
    None,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SearchParameters {
    pub database: String,
    pub database_version: String,
    pub taxonomy: String,
    pub charges: String,
    pub mass_type: PeakMassType,
    pub fixed_modifications: Vec<String>,
    pub variable_modifications: Vec<String>,
    pub missed_cleavages: u32,
    pub fragment_tolerance: Tolerance,
    pub precursor_tolerance: Tolerance,
    pub digestion_enzyme: String,
    pub digestion_regex: String,
    pub enzyme_specificity: EnzymeTermSpecificity,
    pub metadata: MetaInfo,
}
impl Default for SearchParameters {
    fn default() -> Self {
        Self {
            database: String::new(),
            database_version: String::new(),
            taxonomy: String::new(),
            charges: String::new(),
            mass_type: PeakMassType::Monoisotopic,
            fixed_modifications: Vec::new(),
            variable_modifications: Vec::new(),
            missed_cleavages: 0,
            fragment_tolerance: Tolerance::Absolute(0.0),
            precursor_tolerance: Tolerance::Absolute(0.0),
            digestion_enzyme: "unknown_enzyme".into(),
            digestion_regex: String::new(),
            enzyme_specificity: EnzymeTermSpecificity::Unknown,
            metadata: MetaInfo::new(),
        }
    }
}
impl SearchParameters {
    pub fn validate(&self) -> Result<()> {
        for tolerance in [self.fragment_tolerance, self.precursor_tolerance] {
            let (Tolerance::Absolute(v) | Tolerance::Ppm(v)) = tolerance;
            if !v.is_finite() || v < 0.0 {
                return Err(bad("search mass tolerances must be finite and nonnegative"));
            }
        }
        validate_meta(&self.metadata)?;
        Ok(())
    }
    /// The source `getChargeRange`: the minimum and maximum charge that
    /// [`SearchParameters::charges`] allows, rule for rule as in C++ (Release
    /// oracle `protein-identification-closure`, cases `charge_range` and
    /// `charge_range2`, 71 strings).
    ///
    /// 1. The whole string as one integer gives `(n, n)`. Integers follow
    ///    `StringUtils::toInt32`: ASCII blanks (space, tab, CR, LF) around them,
    ///    one optional leading `+`, then an optionally negative decimal that must
    ///    fit `i32` and end the string.
    /// 2. Otherwise, with a `,`, every comma-separated item (blanks trimmed) must
    ///    be such an integer; the result is their minimum and maximum.
    /// 3. Otherwise, with a `:`, the two sides are charge values, and three or
    ///    more sides are an error. A charge value drops every `-` and counts as
    ///    negative if it had one, or else drops every `+`, and must then be an
    ///    integer: `3-` is -3, `2+` is 2.
    /// 4. Otherwise, with one to three `-`, the string splits at the first `-`
    ///    (one) or the second (two or three) into two charge values.
    /// 5. Otherwise the range is `(0, 0)` - the empty string, `abc`, `2+ and 3+`
    ///    and an out-of-range `2147483648` alike.
    ///
    /// Reversed ranges stay reversed (`5-2` is `(5, 2)`), and the minus rule
    /// yields the source's own values for odd strings (`1-2-3` is `(-12, 3)`).
    /// No setting is validated against this: the source never checks the
    /// string, and a run whose charges it cannot read loads, sorts and merges.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] where the source throws `ConversionError` (a
    /// list item or range side that is no integer, e.g. `1,,3`, `2-`, `2:`), and
    /// [`Error::MissingInformation`] for more than two `:`-separated sides.
    pub fn charge_range(&self) -> Result<(i32, i32)> {
        let text = self.charges.as_str();
        if let Some(value) = source_int32(text) {
            return Ok((value, value));
        }
        let conversion = || bad("charge string in SearchParameters is not convertible");
        if text.contains(',') {
            let mut range: Option<(i32, i32)> = None;
            for item in text.split(',') {
                let value = source_int32(trim_source_blanks(item)).ok_or_else(conversion)?;
                range = Some(match range {
                    Some((low, high)) => (low.min(value), high.max(value)),
                    None => (value, value),
                });
            }
            return range.ok_or_else(conversion);
        }
        if text.contains(':') {
            let sides: Vec<&str> = text.split(':').collect();
            if sides.len() > 2 {
                return Err(Error::MissingInformation(
                    "Charge string in SearchParameters not parseable.".into(),
                ));
            }
            return Ok((
                source_charge_value(sides[0])?,
                source_charge_value(sides[1])?,
            ));
        }
        let minus: Vec<usize> = text.match_indices('-').map(|(index, _)| index).collect();
        if !minus.is_empty() && minus.len() <= 3 {
            let split = if minus.len() == 1 { minus[0] } else { minus[1] };
            return Ok((
                source_charge_value(&text[..split])?,
                source_charge_value(&text[split + 1..])?,
            ));
        }
        Ok((0, 0))
    }
    /// The source `SearchParameters::operator==`: every typed setting, but not
    /// the settings' own metadata.
    ///
    /// The source compares the fourteen typed members and leaves out the
    /// inherited `MetaInfoInterface`, so two settings that differ only in
    /// metadata - the `SE:` engine keys, `InferenceEngine` or per-engine
    /// settings - are equal there (Release oracle
    /// `protein-identification-closure`, case `sp_equality`). The derived
    /// [`PartialEq`] is structural and compares [`SearchParameters::metadata`]
    /// too; this is the source's relation. The digestion enzyme is compared by
    /// name and cleavage expression, the two parts of the source `Protease`
    /// this record carries.
    pub fn source_equal(&self, other: &Self) -> bool {
        self.database == other.database
            && self.database_version == other.database_version
            && self.taxonomy == other.taxonomy
            && self.charges == other.charges
            && self.mass_type == other.mass_type
            && self.fixed_modifications == other.fixed_modifications
            && self.variable_modifications == other.variable_modifications
            && self.missed_cleavages == other.missed_cleavages
            && self.fragment_tolerance == other.fragment_tolerance
            && self.precursor_tolerance == other.precursor_tolerance
            && self.digestion_enzyme == other.digestion_enzyme
            && self.digestion_regex == other.digestion_regex
            && self.enzyme_specificity == other.enzyme_specificity
    }
    /// Source merging ignores mass type and missed cleavages. Modification
    /// differences are allowed only for labeled_MS1 experiments.
    ///
    /// Like the source, this compares the stored values and validates nothing:
    /// a charge string the source cannot read is still compared as text, and a
    /// NaN tolerance equals nothing (Release oracle
    /// `protein-identification-closure`, case `mergeable_invalid_charges`).
    /// The result is always `Ok`; the `Result` is kept for API stability.
    pub fn mergeable(&self, other: &Self, experiment_type: &str) -> Result<bool> {
        let basename = |path: &str| path.rsplit(['/', '\\']).next().unwrap_or("").to_owned();
        let same = basename(&self.database) == basename(&other.database)
            && self.database_version == other.database_version
            && self.taxonomy == other.taxonomy
            && self.charges == other.charges
            && self.fragment_tolerance == other.fragment_tolerance
            && self.precursor_tolerance == other.precursor_tolerance
            && self.digestion_enzyme == other.digestion_enzyme
            && self.digestion_regex == other.digestion_regex
            && self.enzyme_specificity == other.enzyme_specificity;
        let set = |values: &Vec<String>| values.iter().cloned().collect::<BTreeSet<_>>();
        Ok(same
            && (experiment_type == "labeled_MS1"
                || (set(&self.fixed_modifications) == set(&other.fixed_modifications)
                    && set(&self.variable_modifications) == set(&other.variable_modifications))))
    }
}
/// The blanks `StringUtils::trim` and `toInt32` skip: space, tab, LF and CR.
fn is_source_blank(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\r')
}
fn trim_source_blanks(text: &str) -> &str {
    text.trim_matches(is_source_blank)
}
/// `StringUtils::toInt32`: blanks, one optional `+`, then `std::from_chars`
/// (an optional `-` and at least one digit, in range), then only blanks.
fn source_int32(text: &str) -> Option<i32> {
    let text = trim_source_blanks(text);
    let text = text.strip_prefix('+').unwrap_or(text);
    let digits_from = usize::from(text.starts_with('-'));
    let end = text[digits_from..]
        .find(|c: char| !c.is_ascii_digit())
        .map_or(text.len(), |offset| digits_from + offset);
    if end == digits_from || end != text.len() {
        return None;
    }
    text.parse::<i32>().ok()
}
/// `SearchParameters::getChargeValue_`: a charge written with its sign before
/// or after the number. Any `-` makes it negative and every `-` is dropped;
/// otherwise every `+` is dropped. The rest must be a `toInt32` integer.
fn source_charge_value(text: &str) -> Result<i32> {
    let negative = text.contains('-');
    let digits: String = text
        .chars()
        .filter(|&c| c != if negative { '-' } else { '+' })
        .collect();
    let value = source_int32(&digits)
        .ok_or_else(|| bad("charge value in SearchParameters is not convertible"))?;
    // With every '-' dropped the value is nonnegative, so negating cannot overflow.
    Ok(if negative { -value } else { value })
}

#[derive(Clone, Debug, PartialEq)]
pub struct ProteinIdentification {
    pub identifier: String,
    pub search_engine: String,
    pub search_engine_version: String,
    pub search_parameters: SearchParameters,
    /// Serialized date/time, retained without timezone conversion.
    pub date_time: Option<String>,
    pub score_type: String,
    pub higher_score_better: bool,
    pub significance_threshold: f64,
    pub hits: Vec<ProteinHit>,
    pub protein_groups: Vec<ProteinGroup>,
    pub indistinguishable_groups: Vec<ProteinGroup>,
    pub primary_ms_run_paths: Vec<String>,
    pub raw_ms_run_paths: Vec<String>,
    pub metadata: MetaInfo,
}
impl Default for ProteinIdentification {
    fn default() -> Self {
        Self {
            identifier: String::new(),
            search_engine: String::new(),
            search_engine_version: String::new(),
            search_parameters: SearchParameters::default(),
            date_time: None,
            score_type: String::new(),
            higher_score_better: true,
            significance_threshold: 0.0,
            hits: Vec::new(),
            protein_groups: Vec::new(),
            indistinguishable_groups: Vec::new(),
            primary_ms_run_paths: Vec::new(),
            raw_ms_run_paths: Vec::new(),
            metadata: MetaInfo::new(),
        }
    }
}
impl ProteinIdentification {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn validate(&self) -> Result<()> {
        finite(
            self.significance_threshold,
            "protein significance threshold",
        )?;
        validate_meta(&self.metadata)?;
        self.search_parameters.validate()?;
        for hit in &self.hits {
            hit.validate()?;
        }
        for group in self
            .protein_groups
            .iter()
            .chain(&self.indistinguishable_groups)
        {
            group.validate()?;
        }
        Ok(())
    }
    pub fn sort(&mut self) -> Result<()> {
        self.validate()?;
        self.hits
            .sort_by(|a, b| score_order(a.score, b.score, self.higher_score_better));
        Ok(())
    }
    pub fn find_hit(&self, accession: &str) -> Option<&ProteinHit> {
        self.hits.iter().find(|hit| hit.accession == accession)
    }
    /// The source `findHit`, which returns a mutable iterator to the first hit
    /// with `accession` (`end()` when there is none, here `None`).
    pub fn find_hit_mut(&mut self, accession: &str) -> Option<&mut ProteinHit> {
        self.hits.iter_mut().find(|hit| hit.accession == accession)
    }
    /// The source `ProteinIdentification::operator==`: run metadata, identity,
    /// engine, date, hits, both group lists and score configuration, with the
    /// search settings compared by [`SearchParameters::source_equal`], so their
    /// own metadata does not count (Release oracle `protein-identification-closure`,
    /// case `sp_equality`: two runs differing only in `InferenceEngine` are
    /// equal although [`ProteinIdentification::inference_engine`] differs).
    ///
    /// The run paths are part of the run metadata in the source (`spectra_data`
    /// and `spectra_data_raw`), so both dedicated path lists count. The derived
    /// [`PartialEq`] additionally compares the settings' metadata.
    pub fn source_equal(&self, other: &Self) -> bool {
        self.metadata == other.metadata
            && self.primary_ms_run_paths == other.primary_ms_run_paths
            && self.raw_ms_run_paths == other.raw_ms_run_paths
            && self.identifier == other.identifier
            && self.search_engine == other.search_engine
            && self.search_engine_version == other.search_engine_version
            && self
                .search_parameters
                .source_equal(&other.search_parameters)
            && self.date_time == other.date_time
            && self.hits == other.hits
            && self.protein_groups == other.protein_groups
            && self.indistinguishable_groups == other.indistinguishable_groups
            && self.score_type == other.score_type
            && self.significance_threshold == other.significance_threshold
            && self.higher_score_better == other.higher_score_better
    }
    /// Union of all supplied evidence intervals, independently of run IDs, as
    /// in OpenMS. Every affected protein needs a nonempty sequence. Errors are atomic.
    ///
    /// As in the source, only what the computation reads is checked: the
    /// protein sequences and the evidence of this run's own accessions, whose
    /// positions must be known, ordered and inside the protein. Scores, other
    /// metadata and evidence for accessions that are not hits play no part
    /// (Release oracle `protein-identification-closure`, cases `coverage` and
    /// `validation_free_helpers`: NaN scores and reversed evidence of another
    /// accession leave the result unchanged). An endpoint equal to the protein
    /// length passes the source check and then writes past its coverage
    /// vector; it is an error here.
    pub fn compute_coverage(&mut self, peptides: &[PeptideIdentification]) -> Result<()> {
        self.compute_coverage_of(peptides.iter())
    }
    /// The source `computeCoverage(const ConsensusMap&, bool)`: coverage from
    /// the identifications of every consensus feature, in map order, plus the
    /// unassigned identifications when `use_unassigned_ids` is set.
    ///
    /// # Errors
    ///
    /// As [`ProteinIdentification::compute_coverage`].
    pub fn compute_coverage_from_consensus_map(
        &mut self,
        map: &crate::kernel::ConsensusMap,
        use_unassigned_ids: bool,
    ) -> Result<()> {
        self.compute_coverage_of(consensus_map_ids(map, use_unassigned_ids))
    }
    fn compute_coverage_of<'a>(
        &mut self,
        peptides: impl Iterator<Item = &'a PeptideIdentification>,
    ) -> Result<()> {
        let mut evidence: BTreeMap<&str, Vec<&PeptideEvidence>> = BTreeMap::new();
        for identification in peptides {
            for hit in &identification.hits {
                for item in &hit.evidences {
                    evidence
                        .entry(&item.protein_accession)
                        .or_default()
                        .push(item);
                }
            }
        }
        let mut coverage = Vec::with_capacity(self.hits.len());
        for protein in &self.hits {
            let length = protein.sequence.len();
            if length == 0 {
                return Err(Error::MissingInformation(
                    "protein sequence is required to compute coverage".into(),
                ));
            }
            let mut ranges = Vec::new();
            for item in evidence
                .get(protein.accession.as_str())
                .into_iter()
                .flatten()
            {
                let (Some(start), Some(end)) = (item.start, item.end) else {
                    return Err(Error::MissingInformation(
                        "peptide evidence has unknown start or end".into(),
                    ));
                };
                if start > end || end >= length {
                    return Err(Error::MissingInformation(
                        "peptide evidence positions are reversed or outside the protein".into(),
                    ));
                }
                ranges.push((start, end + 1));
            }
            ranges.sort_unstable();
            let mut covered = 0_usize;
            let mut stop = 0_usize;
            for (start, end) in ranges {
                if end > stop {
                    covered += end - start.max(stop);
                    stop = end;
                }
            }
            coverage.push(100.0 * covered as f64 / length as f64);
        }
        for (hit, value) in self.hits.iter_mut().zip(coverage) {
            hit.coverage = Some(value);
        }
        Ok(())
    }
    /// Collect observed residue/terminal modifications in protein coordinates.
    /// Skip names or full IDs; proteins without any mapped modifications retain
    /// existing entries, matching the source. Anonymous mass tags use their full
    /// IDs and survive cloning without requiring an elemental formula or mass.
    /// Distinct chemical records at the same position are retained, even when
    /// they share a full ID. Equal records deduplicate by complete value.
    /// Unknown positions are errors; failures leave all protein records unchanged.
    ///
    /// As in the source, scores and metadata are not validated, and evidence
    /// for an accession that is not one of this run's hits is not read (Release
    /// oracle `protein-identification-closure`, case `validation_free_helpers`).
    /// Evidence of this run's own accessions must have known, ordered
    /// positions that hold the modified residues, and a mapped position must
    /// lie inside a known protein sequence: the source would store the
    /// wrapped `Size` of an unknown position or a position past the protein.
    pub fn compute_modifications(
        &mut self,
        peptides: &[PeptideIdentification],
        skip: &BTreeSet<String>,
    ) -> Result<()> {
        self.compute_modifications_of(peptides.iter(), skip)
    }
    /// The source `computeModifications(const ConsensusMap&, const StringList&,
    /// bool)`: modifications from the identifications of every consensus
    /// feature, in map order, plus the unassigned identifications when
    /// `use_unassigned_ids` is set.
    ///
    /// # Errors
    ///
    /// As [`ProteinIdentification::compute_modifications`].
    pub fn compute_modifications_from_consensus_map(
        &mut self,
        map: &crate::kernel::ConsensusMap,
        skip: &BTreeSet<String>,
        use_unassigned_ids: bool,
    ) -> Result<()> {
        self.compute_modifications_of(consensus_map_ids(map, use_unassigned_ids), skip)
    }
    fn compute_modifications_of<'a>(
        &mut self,
        peptides: impl Iterator<Item = &'a PeptideIdentification>,
        skip: &BTreeSet<String>,
    ) -> Result<()> {
        let accessions: BTreeSet<&str> = self.hits.iter().map(|h| h.accession.as_str()).collect();
        let mut mapped: BTreeMap<String, BTreeSet<(usize, String, SequenceModification)>> =
            BTreeMap::new();
        for identification in peptides {
            for hit in &identification.hits {
                if !hit.sequence.is_modified() {
                    continue;
                }
                for evidence in &hit.evidences {
                    if !accessions.contains(evidence.protein_accession.as_str()) {
                        continue;
                    }
                    let positions = evidence.positions()?;
                    let mut put = |position, modification: &SequenceModification| {
                        if !skip.contains(modification.name())
                            && !skip.contains(modification.full_id())
                        {
                            mapped
                                .entry(evidence.protein_accession.clone())
                                .or_default()
                                .insert((
                                    position,
                                    modification.full_id().into(),
                                    modification.clone(),
                                ));
                        }
                    };
                    if let Some(m) = hit.sequence.n_terminal_modification() {
                        put(*positions.start(), m);
                    }
                    for index in 0..hit.sequence.len() {
                        if let Some(m) = hit.sequence.residue_modification(index)? {
                            let position = positions
                                .start()
                                .checked_add(index)
                                .ok_or_else(|| bad("protein modification position overflow"))?;
                            if position > *positions.end() {
                                return Err(bad("modified residue falls beyond evidence endpoint"));
                            }
                            put(position, m);
                        }
                    }
                    if let Some(m) = hit.sequence.c_terminal_modification() {
                        put(*positions.end(), m);
                    }
                }
            }
        }
        let mut next = self.hits.clone();
        for hit in &mut next {
            if let Some(modifications) = mapped.get(&hit.accession) {
                hit.modifications = modifications
                    .iter()
                    .map(|(position, _, modification)| ProteinModification {
                        position: *position,
                        modification: modification.clone(),
                    })
                    .collect();
                if !hit.sequence.is_empty()
                    && hit
                        .modifications
                        .iter()
                        .any(|m| m.position >= hit.sequence.len())
                {
                    return Err(bad("protein modification position outside sequence"));
                }
            }
        }
        self.hits = next;
        Ok(())
    }
}

/// The identifications a source `ConsensusMap` overload visits: every
/// feature's, in map order, then the unassigned ones when requested.
fn consensus_map_ids(
    map: &crate::kernel::ConsensusMap,
    use_unassigned_ids: bool,
) -> impl Iterator<Item = &PeptideIdentification> {
    let unassigned: &[PeptideIdentification] = if use_unassigned_ids {
        &map.unassigned_peptide_identifications
    } else {
        &[]
    };
    map.features
        .iter()
        .flat_map(|feature| feature.base.peptide_identifications.iter())
        .chain(unassigned)
}
