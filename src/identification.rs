// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Owned peptide/protein identification records, scores and protein evidence.
//! Positions are zero-based and inclusive, matching the OpenMS metadata model.
//! Missing values use Option rather than numeric sentinels or NaN.

pub mod graph;
mod protein_run;

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
fn text_value(metadata: &MetaInfo, key: &str) -> String {
    metadata
        .get(key)
        .map(ToString::to_string)
        .unwrap_or_default()
}

/// Flanking amino acid or one of OpenMS's X/[ /] markers.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
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
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
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

/// Spectrum-level peptide candidates. Missing coordinates are None.
#[derive(Clone, Debug, PartialEq)]
pub struct PeptideIdentification {
    pub identifier: String,
    pub hits: Vec<PeptideHit>,
    pub score_type: String,
    pub higher_score_better: bool,
    pub significance_threshold: f64,
    pub rt: Option<f64>,
    pub mz: Option<f64>,
    pub metadata: MetaInfo,
}
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
    pub fn new() -> Self {
        Self::default()
    }
    pub fn len(&self) -> usize {
        self.hits.len()
    }
    pub fn is_empty(&self) -> bool {
        self.hits.is_empty()
    }
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
    pub fn referencing_hits<'a>(&'a self, accessions: &BTreeSet<String>) -> Vec<&'a PeptideHit> {
        self.hits
            .iter()
            .filter(|hit| {
                hit.evidences
                    .iter()
                    .any(|e| accessions.contains(&e.protein_accession))
            })
            .collect()
    }
    pub fn spectrum_reference(&self) -> String {
        text_value(&self.metadata, "spectrum_reference")
    }
    pub fn set_spectrum_reference(&mut self, value: impl Into<String>) {
        self.metadata
            .insert("spectrum_reference".into(), value.into().into());
    }
    pub fn experiment_label(&self) -> String {
        text_value(&self.metadata, "experiment_label")
    }
    pub fn set_experiment_label(&mut self, value: impl Into<String>) {
        let value = value.into();
        if value.is_empty() {
            self.metadata.remove("experiment_label");
        } else {
            self.metadata
                .insert("experiment_label".into(), value.into());
        }
    }
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
        self.charge_range()?;
        Ok(())
    }
    /// Accept one signed charge, a comma-separated list, or ':'/'-' range.
    /// Empty means unspecified; reversed ranges and malformed strings are errors.
    pub fn charge_range(&self) -> Result<Option<(i32, i32)>> {
        let text = self.charges.trim();
        if text.is_empty() {
            return Ok(None);
        }
        let parse = |s: &str| parse_charge(s.trim());
        let pair = if let Ok(value) = parse(text) {
            (value, value)
        } else if text.contains(',') {
            let mut values = text.split(',').map(parse);
            let first = values.next().expect("nonempty split")?;
            let mut pair = (first, first);
            for value in values {
                let value = value?;
                pair = (pair.0.min(value), pair.1.max(value));
            }
            pair
        } else if text.contains(':') {
            let values: Vec<_> = text.split(':').collect();
            if values.len() != 2 {
                return Err(bad("charge range needs two endpoints"));
            }
            (parse(values[0])?, parse(values[1])?)
        } else {
            let mut pairs = Vec::new();
            for (index, _) in text.match_indices('-') {
                if index > 0 {
                    if let (Ok(a), Ok(b)) = (parse(&text[..index]), parse(&text[index + 1..])) {
                        pairs.push((a, b));
                    }
                }
            }
            if pairs.len() != 1 {
                return Err(bad("malformed or ambiguous charge range"));
            }
            pairs[0]
        };
        if pair.0 > pair.1 {
            return Err(bad("charge range endpoints are reversed"));
        }
        Ok(Some(pair))
    }
    /// Source merging ignores mass type and missed cleavages. Modification
    /// differences are allowed only for labeled_MS1 experiments.
    pub fn mergeable(&self, other: &Self, experiment_type: &str) -> Result<bool> {
        self.validate()?;
        other.validate()?;
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
fn parse_charge(text: &str) -> Result<i32> {
    if let Ok(value) = text.parse::<i32>() {
        return Ok(value);
    }
    if let Some(unsigned) = text.strip_suffix('+') {
        return unsigned
            .parse::<u32>()
            .ok()
            .and_then(|v| i32::try_from(v).ok())
            .ok_or_else(|| bad("invalid charge"));
    }
    if let Some(unsigned) = text.strip_suffix('-') {
        return unsigned
            .parse::<u32>()
            .ok()
            .and_then(|v| i32::try_from(-i64::from(v)).ok())
            .ok_or_else(|| bad("invalid charge"));
    }
    Err(bad("invalid charge"))
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
    /// Union of all supplied evidence intervals, independently of run IDs, as
    /// in OpenMS. Every affected protein needs a nonempty sequence. Errors are atomic.
    pub fn compute_coverage(&mut self, peptides: &[PeptideIdentification]) -> Result<()> {
        self.validate()?;
        let mut evidence: BTreeMap<&str, Vec<&PeptideEvidence>> = BTreeMap::new();
        for identification in peptides {
            identification.validate()?;
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
                return Err(bad("protein sequence is required to compute coverage"));
            }
            let mut ranges = Vec::new();
            for item in evidence
                .get(protein.accession.as_str())
                .into_iter()
                .flatten()
            {
                let range = item.positions()?;
                if *range.end() >= length {
                    return Err(bad(
                        "inclusive peptide evidence endpoint is outside protein sequence",
                    ));
                }
                ranges.push((*range.start(), *range.end() + 1));
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
    pub fn compute_modifications(
        &mut self,
        peptides: &[PeptideIdentification],
        skip: &BTreeSet<String>,
    ) -> Result<()> {
        self.validate()?;
        let mut mapped: BTreeMap<String, BTreeSet<(usize, String, SequenceModification)>> =
            BTreeMap::new();
        for identification in peptides {
            identification.validate()?;
            for hit in &identification.hits {
                if !hit.sequence.is_modified() {
                    continue;
                }
                for evidence in &hit.evidences {
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
                hit.validate()?;
            }
        }
        self.hits = next;
        Ok(())
    }
}
