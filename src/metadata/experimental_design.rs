// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Experimental design records and the source grouping/mapping operations.
//! See `docs/EXPERIMENTAL_DESIGN_SUPPORT.md` for the supported subset.

use crate::identification::ProteinIdentification;
use crate::kernel::{ConsensusMap, FeatureMap};
use crate::metadata::MetaValue;
use crate::system::file::basename;
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn limit() -> Error {
    bad("experimental design exceeds its row or byte limit")
}

/// One row of the MS file section: one quantitative channel of one MS file.
///
/// `sample` is the zero-based row index into the sample section; `sample_name`
/// is the `Sample` column value. The two rarely coincide.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MSFileSectionEntry {
    /// 1-based, consecutive over the whole design.
    pub fraction_group: u32,
    /// 1-based; 1 everywhere for unfractionated data.
    pub fraction: u32,
    pub path: String,
    /// 1-based channel position within `path`; 1 for label-free.
    pub label: u32,
    /// Zero-based row index into the sample section.
    pub sample: u32,
    /// The `Sample` column value, an arbitrary name.
    pub sample_name: String,
}

impl Default for MSFileSectionEntry {
    fn default() -> Self {
        Self {
            fraction_group: 1,
            fraction: 1,
            path: "UNKNOWN_FILE".into(),
            label: 1,
            sample: 0,
            sample_name: "0".into(),
        }
    }
}

/// The sample section: one named row per sample, one column per factor.
///
/// Factor values are never interpreted, only compared. A section parsed from a
/// file carries the `Sample` column itself among its factors; the mapping
/// operations drop it explicitly before comparing rows.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SampleSection {
    content: Vec<Vec<String>>,
    sample_to_rowindex: BTreeMap<String, usize>,
    columnname_to_columnindex: BTreeMap<String, usize>,
}

impl SampleSection {
    pub fn new() -> Self {
        Self::default()
    }

    /// Build from a parsed table. Unlike the source constructor this rejects a
    /// row shorter than the column map, which the source reads out of bounds
    /// (see `OpenMS_CPP_ISSUES.md`, CPP-059).
    pub fn from_table(
        content: Vec<Vec<String>>,
        sample_to_rowindex: BTreeMap<String, usize>,
        columnname_to_columnindex: BTreeMap<String, usize>,
    ) -> Result<Self> {
        let columns = columnname_to_columnindex
            .values()
            .copied()
            .max()
            .map(|index| index + 1)
            .unwrap_or(0);
        for row in &content {
            if row.len() < columns {
                return Err(bad("sample section row is shorter than its column map"));
            }
        }
        for row in sample_to_rowindex.values() {
            if *row >= content.len() {
                return Err(bad("sample section name refers to a missing row"));
            }
        }
        Ok(Self {
            content,
            sample_to_rowindex,
            columnname_to_columnindex,
        })
    }

    /// All sample names, in lexical order.
    pub fn samples(&self) -> impl ExactSizeIterator<Item = &str> {
        self.sample_to_rowindex.keys().map(String::as_str)
    }

    /// Append a sample row. As in the source, a repeated name keeps its first
    /// row index while still appending a content row, so [`Self::len`] then
    /// exceeds the number of distinct names.
    pub fn add_sample(&mut self, name: impl Into<String>, content: Vec<String>) {
        let next = self.sample_to_rowindex.len();
        self.sample_to_rowindex.entry(name.into()).or_insert(next);
        self.content.push(content);
    }

    /// All factor (column) names, in lexical order. Empty for a section built
    /// with [`Self::add_sample`].
    pub fn factors(&self) -> impl ExactSizeIterator<Item = &str> {
        self.columnname_to_columnindex.keys().map(String::as_str)
    }
    pub fn has_sample(&self, sample: &str) -> bool {
        self.sample_to_rowindex.contains_key(sample)
    }
    pub fn has_factor(&self, factor: &str) -> bool {
        self.columnname_to_columnindex.contains_key(factor)
    }

    pub fn factor_value(&self, sample_name: &str, factor: &str) -> Result<&str> {
        let row = *self
            .sample_to_rowindex
            .get(sample_name)
            .ok_or_else(|| bad("sample is not present in the experimental design"))?;
        self.value_at(row, factor)
    }
    /// Factor value addressed by zero-based sample row index.
    pub fn factor_value_by_row(&self, sample_row: u32, factor: &str) -> Result<&str> {
        self.value_at(sample_row as usize, factor)
    }
    fn value_at(&self, row: usize, factor: &str) -> Result<&str> {
        let column = self
            .columnname_to_columnindex
            .get(factor)
            .ok_or_else(|| bad("factor is not present in the experimental design"))?;
        self.content
            .get(row)
            .and_then(|values| values.get(*column))
            .map(String::as_str)
            .ok_or_else(|| bad("sample section has no value at that row and factor"))
    }

    /// Storage column index of a factor. A one-table design orders its sample
    /// columns alphabetically, a two-table design keeps the file order.
    pub fn factor_column_index(&self, factor: &str) -> Result<usize> {
        self.columnname_to_columnindex
            .get(factor)
            .copied()
            .ok_or_else(|| bad("factor is not present in the experimental design"))
    }

    /// Name of the sample in a zero-based row. Reads the `Sample` column when
    /// the section has one, else the name store filled by every build path.
    pub fn sample_name(&self, sample_row: u32) -> Result<&str> {
        let row = sample_row as usize;
        if let Some(column) = self.columnname_to_columnindex.get("Sample") {
            if let Some(name) = self.content.get(row).and_then(|values| values.get(*column)) {
                return Ok(name);
            }
        }
        self.sample_to_rowindex
            .iter()
            .find(|(_, stored)| **stored == row)
            .map(|(name, _)| name.as_str())
            .ok_or_else(|| bad("sample section has no sample in that row"))
    }

    /// Zero-based row index of a sample name.
    pub fn sample_row(&self, sample: &str) -> Result<u32> {
        self.sample_to_rowindex
            .get(sample)
            .and_then(|row| u32::try_from(*row).ok())
            .ok_or_else(|| bad("the sample section has no sample with that name"))
    }

    /// Number of content rows, which is what `getNumberOfSamples` reports.
    pub fn len(&self) -> usize {
        self.content.len()
    }
    pub fn is_empty(&self) -> bool {
        self.content.is_empty()
    }
    /// Row index to sample name, reversing the name store. Needed because
    /// [`Self::sample_name`] reads the `Sample` column an inferred design lacks.
    fn row_to_name(&self) -> BTreeMap<usize, &str> {
        self.sample_to_rowindex
            .iter()
            .map(|(name, row)| (*row, name.as_str()))
            .collect()
    }
}

/// An experimental design: the MS file section plus the sample section.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ExperimentalDesign {
    msfile_section: Vec<MSFileSectionEntry>,
    sample_section: SampleSection,
}

impl ExperimentalDesign {
    /// Maximum MS file section rows accepted by one operation.
    pub const MAX_ROWS: usize = 1_000_000;
    /// Conservative cumulative owned payload per operation.
    pub const MAX_BYTES: usize = 256 * 1024 * 1024;

    pub fn new() -> Self {
        Self::default()
    }

    /// Sorts the file section and enforces the source consistency rules.
    pub fn from_sections(
        msfile_section: Vec<MSFileSectionEntry>,
        sample_section: SampleSection,
    ) -> Result<Self> {
        let mut design = Self {
            msfile_section,
            sample_section,
        };
        preflight(&design.msfile_section)?;
        design.sort();
        design.validate()?;
        Ok(design)
    }

    pub fn ms_file_section(&self) -> &[MSFileSectionEntry] {
        &self.msfile_section
    }
    /// Replaces and sorts the file section. Like the source setter this does
    /// not revalidate; use [`Self::validate`] or [`Self::from_sections`].
    pub fn set_ms_file_section(&mut self, msfile_section: Vec<MSFileSectionEntry>) -> Result<()> {
        preflight(&msfile_section)?;
        self.msfile_section = msfile_section;
        self.sort();
        Ok(())
    }
    pub fn sample_section(&self) -> &SampleSection {
        &self.sample_section
    }
    pub fn set_sample_section(&mut self, sample_section: SampleSection) {
        self.sample_section = sample_section;
    }

    fn sort(&mut self) {
        self.msfile_section.sort_by(|a, b| {
            (a.fraction_group, a.fraction, a.label, a.sample, &a.path).cmp(&(
                b.fraction_group,
                b.fraction,
                b.label,
                b.sample,
                &b.path,
            ))
        });
    }

    /// The source consistency rules: unique `(fraction group, fraction, label)`
    /// and `(path, label)`, fraction groups consecutive from 1, and one sample
    /// per `(fraction group, label)` in a design with a single distinct label.
    /// An empty file section is valid and unchecked.
    pub fn validate(&self) -> Result<()> {
        if self.msfile_section.is_empty() {
            return Ok(());
        }
        let mut group_fraction_label = BTreeSet::new();
        let mut path_label = BTreeSet::new();
        let mut group_label_to_sample: BTreeMap<(u32, u32), BTreeSet<u32>> = BTreeMap::new();
        let mut labels = BTreeSet::new();
        let mut fraction_groups = BTreeSet::new();
        for row in &self.msfile_section {
            if !group_fraction_label.insert((row.fraction_group, row.fraction, row.label)) {
                return Err(bad(
                    "(fraction group, fraction, label) combination can only appear once",
                ));
            }
            if !path_label.insert((row.path.as_str(), row.label)) {
                return Err(bad("(path, label) combination can only appear once"));
            }
            group_label_to_sample
                .entry((row.fraction_group, row.label))
                .or_default()
                .insert(row.sample);
            labels.insert(row.label);
            fraction_groups.insert(row.fraction_group);
        }
        for (expected, group) in (1u32..).zip(&fraction_groups) {
            if *group != expected {
                return Err(bad(
                    "fraction groups have to be consecutive integers starting with 1",
                ));
            }
        }
        if labels.len() <= 1 {
            for samples in group_label_to_sample.values() {
                if samples.len() > 1 {
                    return Err(bad(
                        "multiple samples for the same fraction group and label in a \
                         single-label design",
                    ));
                }
            }
        }
        Ok(())
    }

    /// Fraction index to the file paths measured for it, ordered by the sorted
    /// file section.
    pub fn fraction_to_ms_files_mapping(&self) -> BTreeMap<u32, Vec<String>> {
        let mut mapping: BTreeMap<u32, Vec<String>> = BTreeMap::new();
        for row in &self.msfile_section {
            mapping
                .entry(row.fraction)
                .or_default()
                .push(row.path.clone());
        }
        mapping
    }

    fn path_label_mapper(
        &self,
        basename_only: bool,
        select: impl Fn(&MSFileSectionEntry) -> u32,
    ) -> Result<BTreeMap<(String, u32), u32>> {
        let mut mapping = BTreeMap::new();
        for row in &self.msfile_section {
            let path = if basename_only {
                basename(&row.path).to_owned()
            } else {
                row.path.clone()
            };
            let value = select(row);
            match mapping.entry((path, row.label)) {
                std::collections::btree_map::Entry::Vacant(slot) => {
                    slot.insert(value);
                }
                std::collections::btree_map::Entry::Occupied(slot) => {
                    if *slot.get() != value {
                        return Err(bad("ambiguous path/basename and label mapping"));
                    }
                }
            }
        }
        Ok(mapping)
    }

    /// Groups samples whose metadata rows agree over every factor except
    /// `Sample`, replicate columns included. The weaker of the two rules; see
    /// [`Self::condition_to_sample_mapping`] for the other.
    ///
    /// A section without factors groups every sample under the empty tuple.
    pub fn unique_sample_row_to_sample_mapping(
        &self,
    ) -> Result<BTreeMap<Vec<String>, BTreeSet<String>>> {
        self.group_samples_by(|factor| factor != "Sample")
            .map(|grouped| {
                grouped
                    .into_iter()
                    .map(|(values, names)| (values, names.into_iter().map(str::to_owned).collect()))
                    .collect()
            })
    }

    fn group_samples_by(
        &self,
        keep: impl Fn(&str) -> bool,
    ) -> Result<BTreeMap<Vec<String>, BTreeSet<&str>>> {
        let factors: Vec<&str> = self.sample_section.factors().filter(|f| keep(f)).collect();
        let mut grouped: BTreeMap<Vec<String>, BTreeSet<&str>> = BTreeMap::new();
        for sample in self.sample_section.samples() {
            let mut values = Vec::with_capacity(factors.len());
            for factor in &factors {
                values.push(self.sample_section.factor_value(sample, factor)?.to_owned());
            }
            grouped.entry(values).or_default().insert(sample);
        }
        Ok(grouped)
    }

    /// Sample name to prefractionation group index, the reverse of
    /// [`Self::unique_sample_row_to_sample_mapping`]. A factor-less section --
    /// as built by the `from_*` constructors -- gives every sample its own group.
    pub fn sample_to_prefractionation_mapping(&self) -> Result<BTreeMap<String, u32>> {
        if self.sample_section.factors().len() == 0 {
            return Ok(self
                .sample_section
                .samples()
                .enumerate()
                .map(|(index, name)| (name.to_owned(), index as u32))
                .collect());
        }
        let grouped = self.group_samples_by(|factor| factor != "Sample")?;
        Ok(index_groups(grouped, |name| name.to_owned()))
    }

    /// Condition to the zero-based sample rows sharing it. A condition is the
    /// unique combination of factor values ignoring `Sample` and every column
    /// whose name contains `replicate` or `Replicate`. Condition numbers are
    /// the lexicographic rank of the value tuple, not a reference level.
    pub fn condition_to_sample_mapping(&self) -> Result<BTreeMap<Vec<String>, BTreeSet<u32>>> {
        let grouped = self.group_samples_by(is_condition_factor)?;
        let mut mapping = BTreeMap::new();
        for (values, names) in grouped {
            let mut rows = BTreeSet::new();
            for name in names {
                rows.insert(self.sample_section.sample_row(name)?);
            }
            mapping.insert(values, rows);
        }
        Ok(mapping)
    }

    /// Sample name to condition index. A factor-less section gives every
    /// sample its own condition.
    pub fn sample_to_condition_mapping(&self) -> Result<BTreeMap<String, u32>> {
        if self.sample_section.factors().len() == 0 {
            return Ok(self
                .sample_section
                .samples()
                .enumerate()
                .map(|(index, name)| (name.to_owned(), index as u32))
                .collect());
        }
        let grouped = self.group_samples_by(is_condition_factor)?;
        Ok(index_groups(grouped, |name| name.to_owned()))
    }

    /// The `(path, label)` pairs of each condition, in condition order, for
    /// merging across replicates. Paths are full paths, as in the source.
    pub fn condition_to_path_label_vector(&self) -> Result<Vec<Vec<(String, u32)>>> {
        let conditions = self.condition_to_sample_mapping()?;
        let path_label_to_sample = self.path_label_to_sample_mapping(false)?;
        let mut result = Vec::with_capacity(conditions.len());
        for rows in conditions.values() {
            let mut entries = Vec::new();
            for row in rows {
                for (key, sample) in &path_label_to_sample {
                    if sample == row {
                        entries.push(key.clone());
                    }
                }
            }
            result.push(entries);
        }
        Ok(result)
    }

    pub fn path_label_to_sample_mapping(
        &self,
        basename_only: bool,
    ) -> Result<BTreeMap<(String, u32), u32>> {
        self.path_label_mapper(basename_only, |row| row.sample)
    }
    pub fn path_label_to_fraction_mapping(
        &self,
        basename_only: bool,
    ) -> Result<BTreeMap<(String, u32), u32>> {
        self.path_label_mapper(basename_only, |row| row.fraction)
    }
    pub fn path_label_to_fraction_group_mapping(
        &self,
        basename_only: bool,
    ) -> Result<BTreeMap<(String, u32), u32>> {
        self.path_label_mapper(basename_only, |row| row.fraction_group)
    }
    pub fn path_label_to_prefractionation_mapping(
        &self,
        basename_only: bool,
    ) -> Result<BTreeMap<(String, u32), u32>> {
        self.resolve_per_sample(basename_only, self.sample_to_prefractionation_mapping()?)
    }
    pub fn path_label_to_condition_mapping(
        &self,
        basename_only: bool,
    ) -> Result<BTreeMap<(String, u32), u32>> {
        self.resolve_per_sample(basename_only, self.sample_to_condition_mapping()?)
    }
    fn resolve_per_sample(
        &self,
        basename_only: bool,
        by_name: BTreeMap<String, u32>,
    ) -> Result<BTreeMap<(String, u32), u32>> {
        let row_to_name = self.sample_section.row_to_name();
        let mut mapping = BTreeMap::new();
        for (key, row) in self.path_label_to_sample_mapping(basename_only)? {
            let name = row_to_name
                .get(&(row as usize))
                .ok_or_else(|| bad("design row references a missing sample"))?;
            let value = *by_name
                .get(*name)
                .ok_or_else(|| bad("design row references an ungrouped sample"))?;
            mapping.entry(key).or_insert(value);
        }
        Ok(mapping)
    }

    /// Number of rows in the sample section, not the highest sample index.
    pub fn number_of_samples(&self) -> u32 {
        self.sample_section.len() as u32
    }
    /// Number of distinct fraction indices used anywhere in the design.
    pub fn number_of_fractions(&self) -> u32 {
        self.msfile_section
            .iter()
            .map(|row| row.fraction)
            .collect::<BTreeSet<_>>()
            .len() as u32
    }
    /// Highest label index, the plex size of a well-formed design; 0 if empty.
    pub fn number_of_labels(&self) -> u32 {
        self.msfile_section
            .iter()
            .map(|row| row.label)
            .max()
            .unwrap_or(0)
    }
    /// Number of distinct MS file paths; a multiplexed file counts once.
    pub fn number_of_ms_files(&self) -> u32 {
        self.msfile_section
            .iter()
            .map(|row| row.path.as_str())
            .collect::<BTreeSet<_>>()
            .len() as u32
    }
    pub fn number_of_fraction_groups(&self) -> u32 {
        self.msfile_section
            .iter()
            .map(|row| row.fraction_group)
            .collect::<BTreeSet<_>>()
            .len() as u32
    }

    /// Zero-based sample row quantified in a fraction group at a label.
    pub fn sample(&self, fraction_group: u32, label: u32) -> Result<u32> {
        self.msfile_section
            .iter()
            .find(|row| row.fraction_group == fraction_group && row.label == label)
            .map(|row| row.sample)
            .ok_or_else(|| bad("no sample for that fraction group and label"))
    }

    /// True when more than one distinct fraction index is used.
    pub fn is_fractionated(&self) -> bool {
        self.number_of_fractions() > 1
    }

    /// True when every fraction index occurs in equally many rows. Necessary
    /// but not sufficient for a uniformly fractionated design.
    pub fn same_nr_of_ms_files_per_fraction(&self) -> bool {
        let mapping = self.fraction_to_ms_files_mapping();
        let mut sizes = mapping.values().map(Vec::len);
        match sizes.next() {
            None => true,
            Some(first) => sizes.all(|size| size == first),
        }
    }

    /// Keeps only rows whose path basename is in `basenames`, then rebuilds the
    /// sample section and sample indices. Returns the number of removed rows.
    /// An empty result is reported through [`Self::ms_file_section`] rather
    /// than the source's fatal log line.
    pub fn filter_by_basenames(&mut self, basenames: &BTreeSet<String>) -> Result<usize> {
        let before = self.msfile_section.len();
        self.msfile_section
            .retain(|row| basenames.contains(basename(&row.path)));
        let removed = before - self.msfile_section.len();

        let mut ordered_samples = Vec::new();
        let mut seen = BTreeSet::new();
        for row in &self.msfile_section {
            if seen.insert(row.sample_name.clone()) {
                ordered_samples.push(row.sample_name.clone());
            }
        }
        let mut ordered_factors: Vec<String> =
            self.sample_section.factors().map(str::to_owned).collect();
        ordered_factors.sort_by_key(|factor| {
            self.sample_section
                .factor_column_index(factor)
                .unwrap_or(usize::MAX)
        });

        let columnname_to_columnindex = ordered_factors
            .iter()
            .enumerate()
            .map(|(index, factor)| (factor.clone(), index))
            .collect();
        let mut sample_to_rowindex = BTreeMap::new();
        let mut content = Vec::with_capacity(ordered_samples.len());
        for name in &ordered_samples {
            sample_to_rowindex.insert(name.clone(), content.len());
            let mut row = vec![String::new(); ordered_factors.len()];
            if self.sample_section.has_sample(name) {
                for (index, factor) in ordered_factors.iter().enumerate() {
                    row[index] = self.sample_section.factor_value(name, factor)?.to_owned();
                }
            }
            content.push(row);
        }
        self.sample_section = SampleSection::from_table(
            content,
            sample_to_rowindex.clone(),
            columnname_to_columnindex,
        )?;
        for row in &mut self.msfile_section {
            row.sample = sample_to_rowindex
                .get(&row.sample_name)
                .and_then(|index| u32::try_from(*index).ok())
                .ok_or_else(|| bad("filtered design row references a missing sample"))?;
        }
        Ok(removed)
    }

    /// Derive a design from a consensus map's column headers. Fractions must
    /// come with a fraction group; without annotated fractions each distinct
    /// file becomes its own fraction group, in order of appearance.
    pub fn from_consensus_map(map: &ConsensusMap) -> Result<Self> {
        let mut files: Vec<&str> = Vec::new();
        for header in map.column_headers.values() {
            if !files.contains(&header.filename.as_str()) {
                files.push(&header.filename);
            }
        }
        let mut msfile_section = Vec::new();
        let mut sample_section = SampleSection::new();
        let mut group_label_to_sample: BTreeMap<(u32, u32), u32> = BTreeMap::new();
        let mut name_to_sample: BTreeMap<String, u32> = BTreeMap::new();
        for header in map.column_headers.values() {
            let mut row = MSFileSectionEntry {
                path: header.filename.clone(),
                label: header.label_as_uint(&map.experiment_type)?,
                ..Default::default()
            };
            match header.metadata.get("fraction") {
                Some(fraction) => {
                    row.fraction = meta_index(fraction, "fraction")?;
                    let group = header.metadata.get("fraction_group").ok_or_else(|| {
                        bad("fractions annotated but no fraction grouping provided")
                    })?;
                    row.fraction_group = meta_index(group, "fraction group")?;
                }
                None => {
                    row.fraction = 1;
                    let index = files
                        .iter()
                        .position(|file| *file == header.filename)
                        .ok_or_else(|| bad("consensus column header has no file"))?;
                    row.fraction_group = u32::try_from(index + 1)
                        .map_err(|_| bad("consensus map has too many files"))?;
                }
            }
            match header.metadata.get("sample_name") {
                Some(name) => {
                    row.sample_name = name.as_str()?.to_owned();
                    let next = name_to_sample.len() as u32;
                    row.sample = *name_to_sample
                        .entry(row.sample_name.clone())
                        .or_insert(next);
                }
                None => {
                    let next = group_label_to_sample.len() as u32;
                    row.sample = *group_label_to_sample
                        .entry((row.fraction_group, row.label))
                        .or_insert(next);
                    row.sample_name = row.sample.to_string();
                }
            }
            if !sample_section.has_sample(&row.sample_name) {
                sample_section.add_sample(row.sample_name.clone(), Vec::new());
            }
            msfile_section.push(row);
        }
        Self::from_sections(msfile_section, sample_section)
    }

    /// Write this design's fraction structure onto a consensus map's column
    /// headers, the inverse of [`Self::from_consensus_map`]. Headers are matched
    /// on `(basename, label)`; a header without a matching design row, and each
    /// of several headers resolving to one row, is left unannotated. Returns the
    /// number of unannotated headers.
    pub fn annotate_column_headers(&self, map: &mut ConsensusMap) -> Result<usize> {
        let to_group = self.path_label_to_fraction_group_mapping(true)?;
        let to_fraction = self.path_label_to_fraction_mapping(true)?;
        let to_sample = self.path_label_to_sample_mapping(true)?;

        let mut key_uses: BTreeMap<(String, u32), usize> = BTreeMap::new();
        for header in map.column_headers.values() {
            let key = (
                basename(&header.filename).to_owned(),
                header.label_as_uint(&map.experiment_type)?,
            );
            *key_uses.entry(key).or_default() += 1;
        }

        let mut unannotated = 0;
        for header in map.column_headers.values_mut() {
            let key = (
                basename(&header.filename).to_owned(),
                header.label_as_uint(&map.experiment_type)?,
            );
            let Some(group) = to_group.get(&key) else {
                unannotated += 1;
                continue;
            };
            if key_uses.get(&key).copied().unwrap_or(0) > 1 {
                unannotated += 1;
                continue;
            }
            header
                .metadata
                .insert("fraction_group".into(), i64::from(*group).into());
            if let Some(fraction) = to_fraction.get(&key) {
                header
                    .metadata
                    .insert("fraction".into(), i64::from(*fraction).into());
            }
            if let Some(sample) = to_sample.get(&key) {
                if let Ok(name) = self.sample_section.sample_name(*sample) {
                    header
                        .metadata
                        .insert("sample_name".into(), name.to_owned().into());
                }
            }
        }
        Ok(unannotated)
    }

    /// Derive a design from a feature map: one file, fraction, sample and
    /// fraction group. The map must name exactly one primary MS run.
    pub fn from_feature_map(map: &FeatureMap) -> Result<Self> {
        let paths = primary_ms_run_paths(&map.metadata)?;
        if paths.len() != 1 {
            return Err(bad(
                "feature map must be annotated with exactly one MS file",
            ));
        }
        let row = MSFileSectionEntry {
            path: paths[0].clone(),
            ..Default::default()
        };
        let mut sample_section = SampleSection::new();
        sample_section.add_sample(row.sample_name.clone(), Vec::new());
        Self::from_sections(vec![row], sample_section)
    }

    /// Derive a design from identification runs, assuming a label-free,
    /// unfractionated experiment: one fraction group and sample per run path,
    /// in the order the paths appear.
    pub fn from_identifications(proteins: &[ProteinIdentification]) -> Result<Self> {
        let mut msfile_section = Vec::new();
        let mut sample_section = SampleSection::new();
        for path in proteins
            .iter()
            .flat_map(|protein| &protein.primary_ms_run_paths)
        {
            let sample = u32::try_from(msfile_section.len())
                .map_err(|_| bad("too many identification run paths"))?;
            msfile_section.push(MSFileSectionEntry {
                path: path.clone(),
                fraction: 1,
                sample,
                sample_name: sample.to_string(),
                fraction_group: sample + 1,
                label: 1,
            });
            sample_section.add_sample(sample.to_string(), Vec::new());
        }
        Self::from_sections(msfile_section, sample_section)
    }
}

/// A condition ignores `Sample` and every column named as a replicate. The rule
/// matches on the column name only: `Donor` or `Rep` still split conditions.
fn is_condition_factor(factor: &str) -> bool {
    factor != "Sample" && !factor.contains("replicate") && !factor.contains("Replicate")
}

fn index_groups<'a, K: Ord>(
    grouped: BTreeMap<Vec<String>, BTreeSet<&'a str>>,
    key: impl Fn(&'a str) -> K,
) -> BTreeMap<K, u32> {
    let mut mapping = BTreeMap::new();
    for (index, names) in grouped.into_values().enumerate() {
        for name in names {
            mapping.entry(key(name)).or_insert(index as u32);
        }
    }
    mapping
}

/// Source `spectra_data` string list, the primary MS run paths of a map.
fn primary_ms_run_paths(metadata: &crate::metadata::MetaInfo) -> Result<Vec<String>> {
    match metadata.get("spectra_data") {
        Some(value) => Ok(value.as_string_list()?.to_vec()),
        None => Ok(Vec::new()),
    }
}

/// Source `static_cast<unsigned>` of an annotated index, rejecting the
/// nonnumeric and negative values the source casts unchecked.
fn meta_index(value: &MetaValue, what: &str) -> Result<u32> {
    u32::try_from(value.as_i64()?).map_err(|_| bad(&format!("{what} must be a nonnegative index")))
}

fn preflight(rows: &[MSFileSectionEntry]) -> Result<()> {
    if rows.len() > ExperimentalDesign::MAX_ROWS {
        return Err(limit());
    }
    // Each row is copied once into the sorted section and charged for up to two
    // sparse map nodes in the derived groupings it can reach.
    let per_row = std::mem::size_of::<MSFileSectionEntry>()
        .checked_mul(8)
        .and_then(|n| n.checked_add(256))
        .ok_or_else(limit)?;
    let mut bytes = rows.len().checked_mul(per_row).ok_or_else(limit)?;
    if bytes > ExperimentalDesign::MAX_BYTES {
        return Err(limit());
    }
    for row in rows {
        for text in [&row.path, &row.sample_name] {
            bytes = text
                .len()
                .checked_mul(4)
                .and_then(|n| bytes.checked_add(n))
                .ok_or_else(limit)?;
            if bytes > ExperimentalDesign::MAX_BYTES {
                return Err(limit());
            }
        }
    }
    Ok(())
}
