// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Long-format MSstats and MSstatsTMT CSV writer.
//!
//! Ports `FORMAT/MSstatsFile.h`. MSstats consumes one row per quantified
//! peptide ion per run: the protein, the modified peptide sequence, the
//! precursor charge, the experimental annotation taken from the experimental
//! design, the intensity and a reference string naming where the intensity came
//! from. [`prepare_lfq`](crate::format::msstats::prepare_lfq) builds the
//! label-free layout and [`prepare_iso`](crate::format::msstats::prepare_iso)
//! the isobaric (MSstatsTMT) layout.
//!
//! MSstats counts a run as one `(spectra file, fraction)` pair, while OpenMS
//! splits a run into fractions; [`prepare_lfq`](crate::format::msstats::prepare_lfq)
//! therefore enumerates those pairs and reports the mapping it used.
//!
//! See `docs/MSSTATS_SUPPORT.md` for the API mapping, the preserved source
//! conventions and the native differences.

use crate::kernel::ConsensusMap;
use crate::metadata::{ExperimentalDesign, SampleSection};
use crate::system::file::basename;
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;
use std::path::Path;

/// Placeholder the source writes for an unused or unknown text field.
pub const NA: &str = "NA";
/// Field separator of the produced CSV.
pub const DELIMITER: char = ',';
/// Separator between the accessions of one indistinguishable group.
pub const ACCESSION_DELIMITER: char = ';';
/// Character the source wraps the `Reference` field in, without escaping.
pub const QUOTE: char = '"';
/// Largest number of consensus features one call may process.
pub const MAX_FEATURES: usize = 5_000_000;
/// Largest number of output rows one call may produce.
pub const MAX_LINES: usize = 20_000_000;
/// Largest cumulative owned payload one call may allocate.
pub const MAX_BYTES: usize = 512 * 1024 * 1024;

fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}
fn missing(message: impl Into<String>) -> Error {
    Error::MissingInformation(message.into())
}
/// `StringUtils::toStr(double)` for the retention times the source writes.
fn coordinate_text(value: f64) -> String {
    crate::param::value::format_float(value, true)
}
/// `StringUtils::toStr(float)` for the intensities the source writes, which are
/// `Peak2D::IntensityType` and therefore single precision.
fn intensity_text(value: f32) -> String {
    crate::param::value::format_float32(value, true)
}

/// How the intensities of one peptide ion in one run are combined.
///
/// The source takes this as a free-form string and silently writes the
/// intensity `0` when it matches none of the four aggregating names, because
/// its accumulator is initialised to zero and no branch assigns it. This enum
/// refuses an unknown name instead; see the support document.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum RetentionTimeSummarization {
    /// Write one row per retention time instead of combining them. Standard
    /// MSstats rejects such input; MSstatsTMT does its own aggregation and the
    /// source therefore forces this mode for the isobaric layout.
    Manual,
    /// Largest of the distinct intensities.
    #[default]
    Max,
    /// Smallest of the distinct intensities.
    Min,
    /// Arithmetic mean of the distinct intensities.
    Mean,
    /// Sum of the distinct intensities.
    Sum,
}
impl RetentionTimeSummarization {
    /// Parse one of the source's `retention_time_summarization_method` values.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for any other name.
    pub fn from_name(name: &str) -> Result<Self> {
        match name {
            "manual" => Ok(Self::Manual),
            "max" => Ok(Self::Max),
            "min" => Ok(Self::Min),
            "mean" => Ok(Self::Mean),
            "sum" => Ok(Self::Sum),
            other => Err(invalid(format!(
                "unknown retention time summarization method {other:?}"
            ))),
        }
    }
    /// The source spelling of this method.
    pub fn name(self) -> &'static str {
        match self {
            Self::Manual => "manual",
            Self::Max => "max",
            Self::Min => "min",
            Self::Mean => "mean",
            Self::Sum => "sum",
        }
    }
    fn is_manual(self) -> bool {
        self == Self::Manual
    }
    /// Combine the distinct intensities of one peptide ion in one run.
    ///
    /// Duplicate intensities are collapsed first, because the source stores
    /// them in a `std::set`; `Sum` and `Mean` therefore ignore repetitions, and
    /// `Mean` divides by the number of *distinct* values. The arithmetic is
    /// single precision, as in the source.
    fn combine(self, sorted_unique: &[f32]) -> f32 {
        match self {
            Self::Manual | Self::Max => sorted_unique.last().copied().unwrap_or(0.0),
            Self::Min => sorted_unique.first().copied().unwrap_or(0.0),
            Self::Sum | Self::Mean => {
                let mut total = 0.0f32;
                for value in sorted_unique {
                    total += *value;
                }
                if self == Self::Mean {
                    total / sorted_unique.len() as f32
                } else {
                    total
                }
            }
        }
    }
}

/// The rows one conversion produced, with everything the source logged.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MSstatsReport {
    /// The header line followed by one row per written intensity.
    pub lines: Vec<String>,
    /// Diagnostics the source writes to its warning log.
    pub warnings: Vec<String>,
    /// Peptide hits dropped because they map into more than one
    /// indistinguishable protein group.
    pub shared_peptides_dropped: usize,
    /// MSstats run number to OpenMS fraction group, which the source prints to
    /// standard output. The isobaric layout leaves this empty, as the source
    /// does: it declares the mapping and never fills it.
    pub run_to_fraction_group: BTreeMap<u32, u32>,
}

/// Settings of the label-free layout, one per `storeLFQ` argument.
#[derive(Clone, Debug)]
pub struct LfqOptions {
    /// Replaces the consensus map's own run paths when not empty.
    pub reannotate_filenames: Vec<String>,
    /// Write `H` in the `IsotopeLabelType` column instead of `L`. MSstats
    /// documents `L` for endogenous and `H` for labelled reference peptides;
    /// the source carries a note doubting that DDA label-free is ever `H`.
    pub is_isotope_label_type: bool,
    /// Sample-section column holding the biological replicate.
    pub bioreplicate: String,
    /// Sample-section column holding the condition.
    pub condition: String,
    /// How to combine intensities across retention times.
    pub retention_time_summarization: RetentionTimeSummarization,
    /// Drop peptides that map into more than one indistinguishable protein
    /// group. On by default, as in the source.
    pub remove_shared_peptides: bool,
    /// Ceiling on consensus features processed.
    pub max_features: usize,
    /// Ceiling on rows produced.
    pub max_lines: usize,
    /// Ceiling on cumulative owned payload.
    pub max_bytes: usize,
}
impl Default for LfqOptions {
    fn default() -> Self {
        Self {
            reannotate_filenames: Vec::new(),
            is_isotope_label_type: false,
            bioreplicate: "MSstats_BioReplicate".into(),
            condition: "MSstats_Condition".into(),
            retention_time_summarization: RetentionTimeSummarization::Max,
            remove_shared_peptides: true,
            max_features: MAX_FEATURES,
            max_lines: MAX_LINES,
            max_bytes: MAX_BYTES,
        }
    }
}

/// Settings of the isobaric layout, one per `storeISO` argument.
#[derive(Clone, Debug)]
pub struct IsoOptions {
    /// Replaces the consensus map's own run paths when not empty.
    pub reannotate_filenames: Vec<String>,
    /// Sample-section column holding the biological replicate.
    pub bioreplicate: String,
    /// Sample-section column holding the condition.
    pub condition: String,
    /// Sample-section column holding the TMT mixture.
    pub mixture: String,
    /// How to combine intensities across retention times. MSstatsTMT does its
    /// own aggregation, so anything but
    /// [`RetentionTimeSummarization::Manual`] is reverted to it with a warning,
    /// exactly as the source does.
    pub retention_time_summarization: RetentionTimeSummarization,
    /// Drop peptides that map into more than one indistinguishable protein
    /// group. On by default, as in the source.
    pub remove_shared_peptides: bool,
    /// Ceiling on consensus features processed.
    pub max_features: usize,
    /// Ceiling on rows produced.
    pub max_lines: usize,
    /// Ceiling on cumulative owned payload.
    pub max_bytes: usize,
}
impl Default for IsoOptions {
    fn default() -> Self {
        Self {
            reannotate_filenames: Vec::new(),
            bioreplicate: "MSstats_BioReplicate".into(),
            condition: "MSstats_Condition".into(),
            mixture: "MSstats_Mixture".into(),
            retention_time_summarization: RetentionTimeSummarization::Manual,
            remove_shared_peptides: true,
            max_features: MAX_FEATURES,
            max_lines: MAX_LINES,
            max_bytes: MAX_BYTES,
        }
    }
}

struct Limits {
    features: usize,
    lines: usize,
    bytes: usize,
}
impl Limits {
    fn validate(&self) -> Result<()> {
        if self.features == 0 || self.lines == 0 || self.bytes == 0 {
            return Err(invalid("MSstats limits must be positive"));
        }
        Ok(())
    }
    fn spend(&mut self, bytes: usize) -> Result<()> {
        self.bytes = self
            .bytes
            .checked_sub(bytes)
            .ok_or_else(|| invalid("MSstats payload limit exceeded"))?;
        Ok(())
    }
    fn line(&mut self, lines: usize) -> Result<()> {
        self.lines = self
            .lines
            .checked_sub(lines)
            .ok_or_else(|| invalid("MSstats output row limit exceeded"))?;
        Ok(())
    }
}

/// Check that the sample section declares the MSstats condition and biological
/// replicate columns, and warn about a replicate that recurs under more than
/// one condition.
///
/// # Errors
///
/// Returns [`Error::MissingInformation`] when either column is absent; the
/// source throws `Exception::IllegalArgument` naming
/// `MSstats_Condition`/`MSstats_BioReplicate` regardless of the configured
/// column names.
pub fn check_condition_lfq(
    sample_section: &SampleSection,
    bioreplicate: &str,
    condition: &str,
) -> Result<Vec<String>> {
    if !sample_section.has_factor(condition) {
        return Err(missing(format!(
            "Sample Section of the experimental design does not contain {condition}"
        )));
    }
    if !sample_section.has_factor(bioreplicate) {
        return Err(missing(format!(
            "Sample Section of the experimental design does not contain {bioreplicate}"
        )));
    }
    // A BioReplicate label that recurs under two different Conditions declares
    // the same biological unit measured in both, i.e. a paired design, which
    // MSstats models differently. It is a legitimate encoding when intended, so
    // the source warns rather than refusing; so does this.
    let mut replicate_conditions: BTreeMap<&str, BTreeSet<&str>> = BTreeMap::new();
    let samples: Vec<&str> = sample_section.samples().collect();
    for sample in samples {
        replicate_conditions
            .entry(sample_section.factor_value(sample, bioreplicate)?)
            .or_default()
            .insert(sample_section.factor_value(sample, condition)?);
    }
    let mut warnings = Vec::new();
    for (replicate, conditions) in replicate_conditions {
        if conditions.len() > 1 {
            warnings.push(format!(
                "Warning: {bioreplicate} '{replicate}' occurs under {} different {condition} \
                 values ({}). MSstats will treat it as one biological unit measured in each, \
                 i.e. a paired design. If these are unrelated samples, give them distinct \
                 {bioreplicate} values.",
                conditions.len(),
                conditions.into_iter().collect::<Vec<_>>().join(", ")
            ));
        }
    }
    Ok(warnings)
}

/// As [`check_condition_lfq`], and additionally require the mixture column.
///
/// # Errors
///
/// As [`check_condition_lfq`], plus [`Error::MissingInformation`] when the
/// mixture column is absent.
pub fn check_condition_iso(
    sample_section: &SampleSection,
    bioreplicate: &str,
    condition: &str,
    mixture: &str,
) -> Result<Vec<String>> {
    let warnings = check_condition_lfq(sample_section, bioreplicate, condition)?;
    if !sample_section.has_factor(mixture) {
        return Err(missing(format!(
            "Sample Section of the experimental design does not contain {mixture}"
        )));
    }
    Ok(warnings)
}

/// MSstats run numbers, one per distinct `(spectra file basename, fraction)`
/// pair of the design, numbered from one in file-section order.
///
/// MSstats treats runs differently from OpenMS: in MSstats a run is an
/// enumeration of those pairs, while in OpenMS one run is split into fractions.
pub fn assemble_run_map(design: &ExperimentalDesign) -> BTreeMap<(String, u32), u32> {
    let mut map: BTreeMap<(String, u32), u32> = BTreeMap::new();
    for entry in design.ms_file_section() {
        let key = (basename(&entry.path).to_owned(), entry.fraction);
        let next = u32::try_from(map.len())
            .unwrap_or(u32::MAX - 1)
            .saturating_add(1);
        map.entry(key).or_insert(next);
    }
    map
}

/// Accession to indistinguishable-group index, for every accession of every
/// group. A later group wins a repeated accession, as the source's map does.
fn accession_to_group(groups: &[crate::identification::ProteinGroup]) -> BTreeMap<&str, usize> {
    let mut map = BTreeMap::new();
    for (index, group) in groups.iter().enumerate() {
        for accession in &group.accessions {
            map.insert(accession.as_str(), index);
        }
    }
    map
}

/// Whether a peptide with these protein accessions can be quantified in a group
/// context.
///
/// A peptide with no accession cannot; one with a single accession always can;
/// one with several can only when every accession belongs to the same
/// indistinguishable group. An accession the group map does not know is assumed
/// to be a singleton, which disqualifies the peptide.
fn is_quantifyable(accessions: &BTreeSet<&str>, groups: &BTreeMap<&str, usize>) -> bool {
    let mut iterator = accessions.iter();
    let Some(first) = iterator.next() else {
        return false;
    };
    if accessions.len() == 1 {
        return true;
    }
    let Some(group) = groups.get(*first) else {
        return false;
    };
    iterator.all(|accession| groups.get(*accession) == Some(group))
}

/// Whether every name in `first` also appears in `second`.
fn is_subset_of(first: &[String], second: &[String]) -> bool {
    let allowed: BTreeSet<&str> = second.iter().map(String::as_str).collect();
    first.iter().all(|name| allowed.contains(name.as_str()))
}

/// Column-index-ordered basenames of the consensus map's runs, with gaps.
///
/// `FileFilter` leaves gaps in the column indices, so the vector is sized by
/// the largest index and the run paths are consumed in column order. An index
/// with no run path keeps an empty name, exactly as the source leaves it.
fn spectra_paths(map: &ConsensusMap, reannotate: &[String]) -> Result<Vec<String>> {
    let raw: Vec<String> = if reannotate.is_empty() {
        map.primary_ms_run_path()
    } else {
        reannotate.to_vec()
    };
    let mut highest = 0u64;
    for index in map.column_headers.keys() {
        highest = highest.max(*index);
    }
    let size = usize::try_from(highest)
        .ok()
        .and_then(|value| value.checked_add(1))
        .filter(|value| *value <= MAX_FEATURES)
        .ok_or_else(|| invalid("consensus map column index exceeds the MSstats limit"))?;
    let mut paths = vec![String::new(); size];
    // The source consumes one run path per column header, in column order, and
    // stops when it runs out; a column past the last path keeps its empty name.
    for (next, index) in map.column_headers.keys().enumerate() {
        let Some(path) = raw.get(next) else {
            break;
        };
        let position = usize::try_from(*index)
            .map_err(|_| invalid("consensus map column index is not representable"))?;
        paths[position] = basename(path).to_owned();
    }
    Ok(paths)
}

/// One consensus feature's per-handle filenames, intensities, retention times
/// and labels, in handle order.
struct Aggregated {
    filenames: Vec<String>,
    intensities: Vec<f32>,
    retention_times: Vec<f64>,
    labels: Vec<u32>,
}

/// Per-handle information of every consensus feature, in feature order.
///
/// The label is the column header's `channel_id` meta value, or `1` when the
/// header has none — the source's stand-in for a label-free experiment, which
/// its own comment notes is only really a statement about the meta value.
fn aggregate_info(map: &ConsensusMap, paths: &[String]) -> Result<Vec<Aggregated>> {
    let mut result = Vec::with_capacity(map.features.len());
    for feature in &map.features {
        let mut aggregated = Aggregated {
            filenames: Vec::with_capacity(feature.handles().len()),
            intensities: Vec::with_capacity(feature.handles().len()),
            retention_times: Vec::with_capacity(feature.handles().len()),
            labels: Vec::with_capacity(feature.handles().len()),
        };
        for handle in feature.handles() {
            let index = usize::try_from(handle.map_index)
                .map_err(|_| invalid("consensus handle map index is not representable"))?;
            let path = paths.get(index).ok_or_else(|| {
                invalid("consensus handle references a column the map does not declare")
            })?;
            aggregated.filenames.push(path.clone());
            aggregated.intensities.push(handle.intensity);
            aggregated.retention_times.push(handle.rt);
            let header = map.column_headers.get(&handle.map_index).ok_or_else(|| {
                invalid("consensus handle references a column the map does not declare")
            })?;
            aggregated
                .labels
                .push(match header.metadata.get("channel_id") {
                    Some(value) => u32::try_from(value.as_i64()?).map_err(|_| {
                        invalid("consensus channel id must be a non-negative index")
                    })?,
                    None => 1,
                });
        }
        result.push(aggregated);
    }
    Ok(result)
}

/// One output row before its intensity and reference are appended, plus the
/// fields the source's `operator<` orders it by.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct LineKey(Vec<String>);

/// One `(intensity, retention time, reference)` sample, the source's tuple.
type Sample = (f32, f64, String);
/// The rendered row prefix and the samples it aggregates over.
type LineSamples = (String, Vec<Sample>);

struct Rows {
    /// Modified peptide sequence to line key to (rendered prefix, samples).
    /// Both levels are ordered, matching the source's `std::set`/`std::map`.
    by_sequence: BTreeMap<String, BTreeMap<LineKey, LineSamples>>,
}
impl Rows {
    fn new() -> Self {
        Self {
            by_sequence: BTreeMap::new(),
        }
    }
    fn insert(
        &mut self,
        sequence: &str,
        key: LineKey,
        prefix: String,
        sample: Sample,
        limits: &mut Limits,
    ) -> Result<()> {
        limits.spend(
            sequence
                .len()
                .saturating_add(prefix.len())
                .saturating_add(256),
        )?;
        let entry = self
            .by_sequence
            .entry(sequence.to_owned())
            .or_default()
            .entry(key)
            // A line that compares equal to one already present keeps the
            // first rendering, as the source's map insertion does.
            .or_insert((prefix, Vec::new()));
        entry.1.push(sample);
        Ok(())
    }
    /// The rows for every quantifiable sequence, in the source's order.
    fn render(
        &self,
        quantifyable: &BTreeSet<String>,
        summarization: RetentionTimeSummarization,
        limits: &mut Limits,
        warnings: &mut Vec<String>,
    ) -> Result<Vec<String>> {
        let mut lines = Vec::new();
        for sequence in quantifyable {
            let Some(entries) = self.by_sequence.get(sequence) else {
                continue;
            };
            for (prefix, samples) in entries.values() {
                // The source stores these in a set ordered by
                // (intensity, retention time, reference), so both the
                // deduplication below and the chosen reference follow that
                // order rather than insertion order.
                let mut sorted = samples.clone();
                sorted.sort_by(|a, b| {
                    a.0.total_cmp(&b.0)
                        .then(a.1.total_cmp(&b.1))
                        .then_with(|| a.2.cmp(&b.2))
                });
                sorted.dedup_by(|a, b| {
                    a.0.total_cmp(&b.0).is_eq() && a.1.total_cmp(&b.1).is_eq() && a.2 == b.2
                });
                let mut seen_times: Vec<f64> = Vec::new();
                let mut intensities: Vec<f32> = Vec::new();
                for (intensity, retention_time, _) in &sorted {
                    if seen_times
                        .iter()
                        .any(|seen| seen.total_cmp(retention_time).is_eq())
                    {
                        warnings.push(
                            "Peptide ion appears multiple times at the same retention time. This \
                             is not expected."
                                .into(),
                        );
                        continue;
                    }
                    seen_times.push(*retention_time);
                    intensities.push(*intensity);
                }
                intensities.sort_by(f32::total_cmp);
                intensities.dedup_by(|a, b| a.total_cmp(b).is_eq());

                if summarization.is_manual() {
                    // Every stored sample is written, including the duplicate
                    // retention times the loop above only warned about.
                    for (intensity, retention_time, reference) in &sorted {
                        limits.line(1)?;
                        let line = format!(
                            "{}{DELIMITER}{prefix}{DELIMITER}{}{DELIMITER}{QUOTE}{reference}{QUOTE}",
                            coordinate_text(*retention_time),
                            intensity_text(*intensity)
                        );
                        limits.spend(line.len().saturating_add(64))?;
                        lines.push(line);
                    }
                } else {
                    let reference = match sorted.first() {
                        Some((_, _, reference)) => reference.clone(),
                        None => continue,
                    };
                    limits.line(1)?;
                    let line = format!(
                        "{prefix}{DELIMITER}{}{DELIMITER}{QUOTE}{reference}{QUOTE}",
                        intensity_text(summarization.combine(&intensities))
                    );
                    limits.spend(line.len().saturating_add(64))?;
                    lines.push(line);
                }
            }
        }
        Ok(lines)
    }
}

/// What both layouts need after the shared prologue has validated the map
/// against the design. The run paths and design basenames are consumed inside
/// that prologue and are deliberately not carried further.
struct Common<'a> {
    aggregated: Vec<Aggregated>,
    groups: Vec<crate::identification::ProteinGroup>,
    warnings: Vec<String>,
    sample_section: &'a SampleSection,
}

fn prepare_common<'a>(
    map: &ConsensusMap,
    design: &'a ExperimentalDesign,
    reannotate: &[String],
    limits: &Limits,
) -> Result<Common<'a>> {
    limits.validate()?;
    if map.features.len() > limits.features {
        return Err(invalid("MSstats consensus feature limit exceeded"));
    }
    let mut warnings = Vec::new();
    let design_filenames: Vec<String> = design
        .ms_file_section()
        .iter()
        .map(|entry| basename(&entry.path).to_owned())
        .collect();
    let paths = spectra_paths(map, reannotate)?;
    let mut active = Vec::with_capacity(map.column_headers.len());
    for index in map.column_headers.keys() {
        let position = usize::try_from(*index)
            .map_err(|_| invalid("consensus map column index is not representable"))?;
        active.push(
            paths
                .get(position)
                .cloned()
                .ok_or_else(|| invalid("consensus map column index is out of range"))?,
        );
    }
    if !is_subset_of(&active, &design_filenames) {
        return Err(invalid(format!(
            "The filenames (extension ignored) in the consensusXML file are not the same as in \
             the experimental design. Spectra files (consensus map): {}. Spectra files (design): \
             {}",
            active.join(", "),
            design_filenames.join(", ")
        )));
    }
    if active.len() < design_filenames.len() {
        let present: BTreeSet<&str> = active.iter().map(String::as_str).collect();
        let missing_files: Vec<&str> = design_filenames
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<&str>>()
            .into_iter()
            .filter(|name| !present.contains(name))
            .collect();
        warnings.push(format!(
            "Warning: The consensus map contains {} of {} files from the experimental design.\n\
             Missing files: {}\nProceeding with the available subset.",
            active.len(),
            design_filenames.len(),
            missing_files.join(", ")
        ));
    }
    if map.protein_identifications.is_empty() {
        return Err(missing("No protein information found in the ConsensusXML."));
    }
    if map.protein_identifications.len() > 1 {
        warnings.push(format!(
            "Found {} protein runs in consensusXML. Using first one only to parse inference data \
             for now.",
            map.protein_identifications.len()
        ));
    }
    if !map.protein_identifications[0].has_inference_data()? {
        warnings.push(
            "No inference was performed on the first run, defaulting to one-peptide-rule.".into(),
        );
    }
    let aggregated = aggregate_info(map, &paths)?;
    Ok(Common {
        aggregated,
        groups: map.protein_identifications[0]
            .indistinguishable_groups
            .clone(),
        warnings,
        sample_section: design.sample_section(),
    })
}

fn lookup(
    map: &BTreeMap<(String, u32), u32>,
    filename: &str,
    label: u32,
    what: &str,
) -> Result<u32> {
    map.get(&(filename.to_owned(), label))
        .copied()
        .ok_or_else(|| {
            // The source indexes a std::map with operator[], which inserts a zero
            // for a key it does not hold and then reads sample row 0.
            missing(format!(
                "the experimental design has no {what} for file {filename:?} and label {label}"
            ))
        })
}

/// Build the label-free MSstats rows.
///
/// Requires a single-label design; the source throws
/// `Exception::IllegalArgument` naming "Too many labels for a label-free
/// quantitation experiments" otherwise.
///
/// The header is `ProteinName,PeptideSequence,PrecursorCharge,FragmentIon,
/// ProductCharge,IsotopeLabelType,Condition,BioReplicate,Run,Intensity,
/// Reference`, preceded by `RetentionTime,` when the summarization is manual
/// and with `Fraction,` inserted before `Intensity` when the design is
/// fractionated. `FragmentIon` is always `NA` and `ProductCharge` always `0`,
/// because neither is used for DDA data.
///
/// Decoy hits are skipped. A peptide sequence is recorded as quantifiable only
/// once, and a shared peptide is dropped unless
/// [`LfqOptions::remove_shared_peptides`] is cleared.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the design declares more than one
/// label, when the consensus map's run basenames are not a subset of the
/// design's, or when a ceiling is reached; [`Error::MissingInformation`] when a
/// required sample-section column is absent, when the map has no protein
/// identification, or when a `(file, label)` pair has no design row. The last
/// of those is where the source reads sample row 0 instead; see the support
/// document.
pub fn prepare_lfq(
    map: &ConsensusMap,
    design: &ExperimentalDesign,
    options: &LfqOptions,
) -> Result<MSstatsReport> {
    let mut limits = Limits {
        features: options.max_features,
        lines: options.max_lines,
        bytes: options.max_bytes,
    };
    limits.validate()?;
    if design.number_of_labels() != 1 {
        return Err(invalid(
            "Too many labels for a label-free quantitation experiments. Please select the \
             appropriate method, or validate the experimental design.",
        ));
    }
    let mut report = MSstatsReport {
        warnings: check_condition_lfq(
            design.sample_section(),
            &options.bioreplicate,
            &options.condition,
        )?,
        ..Default::default()
    };
    let run_map = assemble_run_map(design);
    let path_label_to_sample = design.path_label_to_sample_mapping(true)?;
    let path_label_to_fraction = design.path_label_to_fraction_mapping(true)?;
    let path_label_to_fraction_group = design.path_label_to_fraction_group_mapping(true)?;
    let manual = options.retention_time_summarization.is_manual();
    if manual {
        report.warnings.push(
            "WARNING: rt_summarization set to manual. One feature might appear at multiple \
             retention times in the output file. This is invalid input for standard MSstats. \
             Combining of features over retention times is recommended!"
                .into(),
        );
    }
    let has_fraction = design.is_fractionated();
    let common = prepare_common(map, design, &options.reannotate_filenames, &limits)?;
    report.warnings.extend(common.warnings.iter().cloned());

    report.lines.push(format!(
        "{}ProteinName,PeptideSequence,PrecursorCharge,FragmentIon,ProductCharge,\
         IsotopeLabelType,Condition,BioReplicate,Run,{}Intensity,Reference",
        if manual { "RetentionTime," } else { "" },
        if has_fraction { "Fraction," } else { "" }
    ));
    let isotope_label_type = if options.is_isotope_label_type {
        "H"
    } else {
        "L"
    };
    let groups = accession_to_group(&common.groups);
    let mut quantifyable = BTreeSet::new();
    let mut rows = Rows::new();

    for (index, feature) in map.features.iter().enumerate() {
        let aggregated = common
            .aggregated
            .get(index)
            .ok_or_else(|| invalid("aggregated consensus information is incomplete"))?;
        for identification in &feature.peptide_identifications {
            for hit in &identification.hits {
                if hit.is_decoy()? {
                    continue;
                }
                let sequence = hit.sequence.to_string();
                let accessions = hit.protein_accessions();
                if options.remove_shared_peptides && !is_quantifyable(&accessions, &groups) {
                    report.shared_peptides_dropped += 1;
                    continue;
                }
                limits.spend(sequence.len().saturating_add(64))?;
                quantifyable.insert(sequence.clone());
                // MSstats user manual 3.7.3: an unknown precursor charge is 0.
                let precursor_charge = hit.charge.to_string();
                let mut accession = accessions
                    .iter()
                    .copied()
                    .collect::<Vec<&str>>()
                    .join(&ACCESSION_DELIMITER.to_string());
                if accession.is_empty() {
                    accession = NA.to_owned();
                }
                for handle in 0..aggregated.filenames.len() {
                    let filename = &aggregated.filenames[handle];
                    let intensity = aggregated.intensities[handle];
                    let retention_time = aggregated.retention_times[handle];
                    let label = aggregated.labels[handle];
                    let sample = lookup(&path_label_to_sample, filename, label, "sample")?;
                    let fraction = lookup(&path_label_to_fraction, filename, label, "fraction")?;
                    let run = *run_map.get(&(filename.clone(), fraction)).ok_or_else(|| {
                        missing(format!(
                            "the experimental design has no run for file {filename:?} and \
                                 fraction {fraction}"
                        ))
                    })?;
                    let fraction_group = lookup(
                        &path_label_to_fraction_group,
                        filename,
                        label,
                        "fraction group",
                    )?;
                    report.run_to_fraction_group.insert(run, fraction_group);

                    let condition = common
                        .sample_section
                        .factor_value_by_row(sample, &options.condition)?
                        .to_owned();
                    let bioreplicate = common
                        .sample_section
                        .factor_value_by_row(sample, &options.bioreplicate)?
                        .to_owned();
                    let run_text = run.to_string();
                    let fraction_text = if has_fraction {
                        fraction.to_string()
                    } else {
                        String::new()
                    };
                    let key = LineKey(vec![
                        accession.clone(),
                        run_text.clone(),
                        condition.clone(),
                        bioreplicate.clone(),
                        precursor_charge.clone(),
                        sequence.clone(),
                    ]);
                    let mut prefix = format!(
                        "{accession},{sequence},{precursor_charge},{NA},0,{isotope_label_type},\
                         {condition},{bioreplicate},{run_text}"
                    );
                    if has_fraction {
                        prefix.push(DELIMITER);
                        prefix.push_str(&fraction_text);
                    }
                    rows.insert(
                        &sequence,
                        key,
                        prefix,
                        (intensity, retention_time, filename.clone()),
                        &mut limits,
                    )?;
                }
            }
        }
    }
    if report.shared_peptides_dropped > 0 {
        report.warnings.push(format!(
            "WARNING: {} peptide hit(s) were dropped because they map to proteins in different \
             indistinguishable protein groups (shared peptides). Use -remove_shared_peptides \
             false to keep them.",
            report.shared_peptides_dropped
        ));
    }
    let mut warnings = Vec::new();
    let lines = rows.render(
        &quantifyable,
        options.retention_time_summarization,
        &mut limits,
        &mut warnings,
    )?;
    report.warnings.extend(warnings);
    report.lines.extend(lines);
    Ok(report)
}

/// Build the isobaric (MSstatsTMT) rows.
///
/// The header is `RetentionTime,ProteinName,PeptideSequence,Charge,Channel,
/// Condition,BioReplicate,Run,Mixture,TechRepMixture,Fraction,Intensity,
/// Reference`. `TechRepMixture` is the mixture factor value joined to the
/// OpenMS fraction group, and `Run` is that joined to the fraction.
///
/// The channel is the column header's `channel_id` plus one, so a map whose
/// headers carry no channel id yields channel `2` — the source's label default
/// of `1` plus one. Decoy hits are skipped and a negative charge is clamped to
/// zero, which MSstats documents as "unknown".
///
/// The reference is the run basename, joined to the identification's
/// `spectrum_reference` meta value, or to `NONATIVEID` when it has none.
///
/// # Errors
///
/// As [`prepare_lfq`], minus the single-label requirement and plus
/// [`Error::MissingInformation`] when the mixture column is absent.
pub fn prepare_iso(
    map: &ConsensusMap,
    design: &ExperimentalDesign,
    options: &IsoOptions,
) -> Result<MSstatsReport> {
    let mut limits = Limits {
        features: options.max_features,
        lines: options.max_lines,
        bytes: options.max_bytes,
    };
    limits.validate()?;
    let mut report = MSstatsReport {
        warnings: check_condition_iso(
            design.sample_section(),
            &options.bioreplicate,
            &options.condition,
            &options.mixture,
        )?,
        ..Default::default()
    };
    let path_label_to_sample = design.path_label_to_sample_mapping(true)?;
    let path_label_to_fraction = design.path_label_to_fraction_mapping(true)?;
    let path_label_to_fraction_group = design.path_label_to_fraction_group_mapping(true)?;
    if !options.retention_time_summarization.is_manual() {
        report.warnings.push(
            "WARNING: rt_summarization set to something else than 'manual' but MSstatsTMT does \
             aggregation of intensities of peptide-chargestate combinations in the same file \
             itself. Reverting to 'manual'"
                .into(),
        );
    }
    let common = prepare_common(map, design, &options.reannotate_filenames, &limits)?;
    report.warnings.extend(common.warnings.iter().cloned());

    report.lines.push(
        "RetentionTime,ProteinName,PeptideSequence,Charge,Channel,Condition,BioReplicate,Run,\
         Mixture,TechRepMixture,Fraction,Intensity,Reference"
            .to_owned(),
    );
    let groups = accession_to_group(&common.groups);
    let mut quantifyable = BTreeSet::new();
    let mut rows = Rows::new();

    for (index, feature) in map.features.iter().enumerate() {
        let aggregated = common
            .aggregated
            .get(index)
            .ok_or_else(|| invalid("aggregated consensus information is incomplete"))?;
        for identification in &feature.peptide_identifications {
            let native_id = match identification.metadata.get("spectrum_reference") {
                Some(value) => value.to_string(),
                None => "NONATIVEID".to_owned(),
            };
            for hit in &identification.hits {
                if hit.is_decoy()? {
                    continue;
                }
                let precursor_charge = hit.charge.max(0).to_string();
                let sequence = hit.sequence.to_string();
                let accessions = hit.protein_accessions();
                if options.remove_shared_peptides && !is_quantifyable(&accessions, &groups) {
                    report.shared_peptides_dropped += 1;
                    continue;
                }
                limits.spend(sequence.len().saturating_add(64))?;
                quantifyable.insert(sequence.clone());
                let mut accession = accessions
                    .iter()
                    .copied()
                    .collect::<Vec<&str>>()
                    .join(&ACCESSION_DELIMITER.to_string());
                if accession.is_empty() {
                    accession = NA.to_owned();
                }
                for handle in 0..aggregated.filenames.len() {
                    let filename = &aggregated.filenames[handle];
                    let intensity = aggregated.intensities[handle];
                    let retention_time = aggregated.retention_times[handle];
                    let channel = aggregated.labels[handle]
                        .checked_add(1)
                        .ok_or_else(|| invalid("consensus channel id overflows"))?;
                    let sample = lookup(&path_label_to_sample, filename, channel, "sample")?;
                    let fraction = lookup(&path_label_to_fraction, filename, channel, "fraction")?;
                    let fraction_group = lookup(
                        &path_label_to_fraction_group,
                        filename,
                        channel,
                        "fraction group",
                    )?;
                    let mixture = common
                        .sample_section
                        .factor_value_by_row(sample, &options.mixture)?
                        .to_owned();
                    let condition = common
                        .sample_section
                        .factor_value_by_row(sample, &options.condition)?
                        .to_owned();
                    let bioreplicate = common
                        .sample_section
                        .factor_value_by_row(sample, &options.bioreplicate)?
                        .to_owned();
                    let technical_replicate = format!("{mixture}_{fraction_group}");
                    let run = format!("{technical_replicate}_{fraction}");
                    let channel_text = channel.to_string();
                    let fraction_text = fraction.to_string();
                    let key = LineKey(vec![
                        accession.clone(),
                        run.clone(),
                        condition.clone(),
                        bioreplicate.clone(),
                        mixture.clone(),
                        precursor_charge.clone(),
                        sequence.clone(),
                        channel_text.clone(),
                    ]);
                    let prefix = format!(
                        "{accession},{sequence},{precursor_charge},{channel_text},{condition},\
                         {bioreplicate},{run},{mixture},{technical_replicate},{fraction_text}"
                    );
                    let reference = format!("{filename}_{native_id}");
                    limits.spend(reference.len().saturating_add(64))?;
                    rows.insert(
                        &sequence,
                        key,
                        prefix,
                        (intensity, retention_time, reference),
                        &mut limits,
                    )?;
                }
            }
        }
    }
    if report.shared_peptides_dropped > 0 {
        report.warnings.push(format!(
            "WARNING: {} peptide hit(s) were dropped because they map to proteins in different \
             indistinguishable protein groups (shared peptides). Use -remove_shared_peptides \
             false to keep them.",
            report.shared_peptides_dropped
        ));
    }
    let mut warnings = Vec::new();
    // MSstatsTMT aggregates itself, so the source always writes every sample.
    let lines = rows.render(
        &quantifyable,
        RetentionTimeSummarization::Manual,
        &mut limits,
        &mut warnings,
    )?;
    report.warnings.extend(warnings);
    report.lines.extend(lines);
    Ok(report)
}

fn emit(writer: &mut impl Write, report: &MSstatsReport) -> Result<()> {
    for line in &report.lines {
        writer.write_all(line.as_bytes())?;
        writer.write_all(b"\n")?;
    }
    writer.flush()?;
    Ok(())
}

/// Write the label-free layout to `writer`.
///
/// # Errors
///
/// As [`prepare_lfq`], plus [`Error::Io`] from the writer.
pub fn write_lfq(
    mut writer: impl Write,
    map: &ConsensusMap,
    design: &ExperimentalDesign,
    options: &LfqOptions,
) -> Result<MSstatsReport> {
    let report = prepare_lfq(map, design, options)?;
    emit(&mut writer, &report)?;
    Ok(report)
}

/// Write the label-free layout to `path`, publishing it only once complete.
///
/// # Errors
///
/// As [`write_lfq`], plus [`Error::Io`] when the output cannot be created; the
/// source's `TextFile::store` throws `Exception::UnableToCreateFile` there.
pub fn store_lfq(
    path: impl AsRef<Path>,
    map: &ConsensusMap,
    design: &ExperimentalDesign,
    options: &LfqOptions,
) -> Result<MSstatsReport> {
    let report = prepare_lfq(map, design, options)?;
    super::path_io::write(path.as_ref(), |writer| {
        for line in &report.lines {
            writer.write_all(line.as_bytes())?;
            writer.write_all(b"\n")?;
        }
        Ok(())
    })?;
    Ok(report)
}

/// Write the isobaric layout to `writer`.
///
/// # Errors
///
/// As [`prepare_iso`], plus [`Error::Io`] from the writer.
pub fn write_iso(
    mut writer: impl Write,
    map: &ConsensusMap,
    design: &ExperimentalDesign,
    options: &IsoOptions,
) -> Result<MSstatsReport> {
    let report = prepare_iso(map, design, options)?;
    emit(&mut writer, &report)?;
    Ok(report)
}

/// Write the isobaric layout to `path`, publishing it only once complete.
///
/// # Errors
///
/// As [`write_iso`], plus [`Error::Io`] when the output cannot be created.
pub fn store_iso(
    path: impl AsRef<Path>,
    map: &ConsensusMap,
    design: &ExperimentalDesign,
    options: &IsoOptions,
) -> Result<MSstatsReport> {
    let report = prepare_iso(map, design, options)?;
    super::path_io::write(path.as_ref(), |writer| {
        for line in &report.lines {
            writer.write_all(line.as_bytes())?;
            writer.write_all(b"\n")?;
        }
        Ok(())
    })?;
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distinct_intensities_are_collapsed_before_they_are_combined() {
        // The source stores intensities in a set, so the mean of 1, 1 and 3 is
        // the mean of the two distinct values, not of the three samples.
        let values = [1.0f32, 3.0];
        assert_eq!(RetentionTimeSummarization::Mean.combine(&values), 2.0);
        assert_eq!(RetentionTimeSummarization::Sum.combine(&values), 4.0);
        assert_eq!(RetentionTimeSummarization::Max.combine(&values), 3.0);
        assert_eq!(RetentionTimeSummarization::Min.combine(&values), 1.0);
    }

    #[test]
    fn unknown_summarization_names_are_refused() {
        assert!(RetentionTimeSummarization::from_name("median").is_err());
        assert_eq!(
            RetentionTimeSummarization::from_name("manual").unwrap(),
            RetentionTimeSummarization::Manual
        );
    }

    #[test]
    fn a_multi_group_peptide_is_not_quantifyable() {
        let groups = BTreeMap::from([("A", 0usize), ("B", 1usize)]);
        assert!(!is_quantifyable(&BTreeSet::from(["A", "B"]), &groups));
        assert!(is_quantifyable(&BTreeSet::from(["A"]), &groups));
        assert!(!is_quantifyable(&BTreeSet::new(), &groups));
        // An accession outside every group is assumed to be a singleton.
        assert!(!is_quantifyable(&BTreeSet::from(["A", "C"]), &groups));
    }
}
