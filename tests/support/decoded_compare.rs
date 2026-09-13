// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Decoded-content comparison of feature maps and peak maps, as shared **test
//! support** (decision D6).
//!
//! The upstream suite compares XML outputs line by line with `FuzzyDiff`. This
//! port's writers do not reproduce the C++ line layout, so XML outputs are
//! compared after decoding: both sides are loaded into [`FeatureMap`] or
//! [`MSExperiment`] values and walked field by field. Every number is judged by
//! the same rule `FuzzyDiff` applies to a number token
//! ([`fuzzy::compare_numbers`] with a fresh relative maximum), structure is exact
//! (counts, order, metadata key sets, hull point order, value types), and the
//! first mismatch is reported with its path, e.g.
//! `features[3].convex_hulls[0].points[2].mz`.
//!
//! It uses the line comparator's number rule, so a test includes both files at
//! its crate root, the comparator under the name `fuzzy`:
//!
//! ```text
//! #[path = "support/fuzzy_string_comparator.rs"]
//! mod fuzzy;
//! #[path = "support/decoded_compare.rs"]
//! mod decoded;
//! ```
//!
//! This file does not load the comparator itself: loading one file as two
//! modules duplicates its types and fails `clippy::duplicate_mod`. Field
//! coverage, the id exclusion and the differences from line-level `FuzzyDiff` are
//! listed in `docs/FUZZY_STRING_COMPARATOR_SUPPORT.md`.
#![allow(dead_code)]

use crate::fuzzy;
use openms::format::FileType;
use openms::kernel::{
    ConvexHull2D, DataArray, Feature, FeatureMap, MSChromatogram, MSExperiment, MSSpectrum,
    Precursor,
};
use openms::metadata::{DataProcessing, MetaInfo, MetaValue, MetaValueData};
use std::fmt;
use std::fmt::Debug;

/// The numeric tolerance: `FuzzyDiff`'s `ratio` and `absdiff`, normalised as the
/// comparator's setters normalise them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tolerance {
    ratio: f64,
    absdiff: f64,
}

impl Tolerance {
    /// A tolerance accepting a pair when the ratio is at most `ratio` or the
    /// absolute difference is at most `absdiff`.
    pub fn new(ratio: f64, absdiff: f64) -> Self {
        Self {
            ratio: fuzzy::normalise_relative(ratio),
            absdiff: fuzzy::normalise_absolute(absdiff),
        }
    }

    /// The tolerance of a `FuzzyDiff` invocation, e.g. of
    /// [`fuzzy::FuzzyDiffSettings::upstream`] (the pinned `FuzzyDiff.ini`).
    pub fn from_settings(settings: &fuzzy::FuzzyDiffSettings) -> Self {
        Self::new(settings.ratio, settings.absdiff)
    }

    /// Exact numeric equality (ratio 1, absdiff 0); NaN still equals NaN.
    pub fn exact() -> Self {
        Self::new(1.0, 0.0)
    }

    /// The normalised relative tolerance.
    pub fn ratio(&self) -> f64 {
        self.ratio
    }

    /// The normalised absolute tolerance.
    pub fn absdiff(&self) -> f64 {
        self.absdiff
    }

    /// Judge one pair of numbers; `Err` carries the comparator's failure message.
    pub fn check(&self, actual: f64, expected: f64) -> Result<(), &'static str> {
        match fuzzy::compare_numbers(actual, expected, self.ratio, self.absdiff, 1.0).failure {
            Some(message) => Err(message),
            None => Ok(()),
        }
    }
}

/// What a decoded comparison checks.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecodedOptions {
    /// The numeric tolerance.
    pub tolerance: Tolerance,
    /// Skip generated identifiers, the decoded counterpart of `-whitelist "id="`:
    /// the unique ids of maps and features (including subordinates), the native
    /// ids of spectra and chromatograms, precursor spectrum references and the SQL
    /// run id. Everything else is still compared.
    pub ignore_unique_ids: bool,
}

impl DecodedOptions {
    /// Compare everything with `tolerance`.
    pub fn new(tolerance: Tolerance) -> Self {
        Self {
            tolerance,
            ignore_unique_ids: false,
        }
    }

    /// The same options, skipping generated identifiers.
    pub fn ignoring_unique_ids(mut self) -> Self {
        self.ignore_unique_ids = true;
        self
    }
}

/// The first difference a decoded comparison found.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DecodedMismatch {
    /// Where the values differ, e.g. `spectra[2].peaks[10].intensity`.
    pub path: String,
    /// What differs, with both values.
    pub detail: String,
}

impl fmt::Display for DecodedMismatch {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}: {}", self.path, self.detail)
    }
}

impl std::error::Error for DecodedMismatch {}

type Outcome = Result<(), DecodedMismatch>;

/// Compare two feature maps; see the module documentation for the rules.
///
/// Compared: the map unique id (unless ignored), document identifier, metadata,
/// data processing, protein and unassigned peptide identifications, and every
/// feature: unique id (unless ignored), RT, m/z, intensity, overall quality,
/// charge (exact), width, RT and m/z quality, metadata, peptide identifications,
/// identification-graph references, convex hulls point by point in order, and
/// subordinates recursively. The load provenance (`loaded_file_path`,
/// `loaded_file_type`) is not content and is not compared.
///
/// # Errors
///
/// Returns the first [`DecodedMismatch`] in walk order.
pub fn compare_feature_maps(
    actual: &FeatureMap,
    expected: &FeatureMap,
    options: &DecodedOptions,
) -> Result<(), DecodedMismatch> {
    let mut walk = Walk::new(options);
    if !options.ignore_unique_ids {
        walk.exact("unique_id", &actual.unique_id, &expected.unique_id)?;
    }
    walk.exact("identifier", &actual.identifier, &expected.identifier)?;
    walk.meta("metadata", &actual.metadata, &expected.metadata)?;
    walk.processing_list(
        "data_processing",
        &actual.data_processing.iter().collect::<Vec<_>>(),
        &expected.data_processing.iter().collect::<Vec<_>>(),
    )?;
    walk.debug_text(
        "protein_identifications",
        &actual.protein_identifications,
        &expected.protein_identifications,
    )?;
    walk.debug_text(
        "unassigned_peptide_identifications",
        &actual.unassigned_peptide_identifications,
        &expected.unassigned_peptide_identifications,
    )?;
    walk.count("features", actual.features.len(), expected.features.len())?;
    for (index, (a, e)) in actual.features.iter().zip(&expected.features).enumerate() {
        walk.nested(format!("features[{index}]"), |walk| walk.feature(a, e, 0))?;
    }
    Ok(())
}

/// Compare two peak maps; see the module documentation for the rules.
///
/// Compared: the experimental settings (without the load provenance), the SQL run
/// id (unless ids are ignored), and every spectrum and chromatogram with their
/// peaks, precursors, products, data arrays, metadata, data processing,
/// instrument settings, acquisition information, source file and peptide
/// identifications.
///
/// # Errors
///
/// Returns the first [`DecodedMismatch`] in walk order.
pub fn compare_experiments(
    actual: &MSExperiment,
    expected: &MSExperiment,
    options: &DecodedOptions,
) -> Result<(), DecodedMismatch> {
    let mut walk = Walk::new(options);
    let strip = |experiment: &MSExperiment| {
        let mut settings = experiment.settings.clone();
        settings.document.loaded_file_path = String::new();
        settings.document.loaded_file_type = FileType::Unknown;
        settings
    };
    walk.debug_text("settings", &strip(actual), &strip(expected))?;
    if !options.ignore_unique_ids {
        walk.exact("sql_run_id", &actual.sql_run_id, &expected.sql_run_id)?;
    }
    walk.count("spectra", actual.spectra.len(), expected.spectra.len())?;
    for (index, (a, e)) in actual.spectra.iter().zip(&expected.spectra).enumerate() {
        walk.nested(format!("spectra[{index}]"), |walk| walk.spectrum(a, e))?;
    }
    walk.count(
        "chromatograms",
        actual.chromatograms.len(),
        expected.chromatograms.len(),
    )?;
    for (index, (a, e)) in actual
        .chromatograms
        .iter()
        .zip(&expected.chromatograms)
        .enumerate()
    {
        walk.nested(format!("chromatograms[{index}]"), |walk| {
            walk.chromatogram(a, e)
        })?;
    }
    Ok(())
}

struct Walk {
    tolerance: Tolerance,
    ignore_unique_ids: bool,
    path: Vec<String>,
}

impl Walk {
    fn new(options: &DecodedOptions) -> Self {
        Self {
            tolerance: options.tolerance,
            ignore_unique_ids: options.ignore_unique_ids,
            path: Vec::new(),
        }
    }

    fn at(&self, leaf: &str) -> String {
        let mut path = self.path.join(".");
        if !leaf.is_empty() {
            if !path.is_empty() && !leaf.starts_with('[') {
                path.push('.');
            }
            path.push_str(leaf);
        }
        path
    }

    fn fail(&self, leaf: &str, detail: String) -> Outcome {
        Err(DecodedMismatch {
            path: self.at(leaf),
            detail,
        })
    }

    fn nested(&mut self, segment: String, body: impl FnOnce(&mut Self) -> Outcome) -> Outcome {
        self.path.push(segment);
        let outcome = body(self);
        self.path.pop();
        outcome
    }

    fn number(&self, leaf: &str, actual: f64, expected: f64) -> Outcome {
        match self.tolerance.check(actual, expected) {
            Ok(()) => Ok(()),
            Err(message) => self.fail(leaf, format!("{actual:?} vs {expected:?} ({message})")),
        }
    }

    fn optional_number(&self, leaf: &str, actual: Option<f64>, expected: Option<f64>) -> Outcome {
        match (actual, expected) {
            (None, None) => Ok(()),
            (Some(a), Some(e)) => self.number(leaf, a, e),
            _ => self.fail(
                leaf,
                format!("{actual:?} vs {expected:?} (presence differs)"),
            ),
        }
    }

    fn exact<T: PartialEq + Debug + ?Sized>(
        &self,
        leaf: &str,
        actual: &T,
        expected: &T,
    ) -> Outcome {
        if actual == expected {
            Ok(())
        } else {
            self.fail(leaf, format!("{actual:?} vs {expected:?}"))
        }
    }

    fn count(&self, leaf: &str, actual: usize, expected: usize) -> Outcome {
        if actual == expected {
            Ok(())
        } else {
            self.fail(leaf, format!("count {actual} vs {expected}"))
        }
    }

    /// Compare nested metadata through its pretty `Debug` rendering, line by line
    /// with the fuzzy comparator and the same tolerance: field names, variants,
    /// strings and element counts must agree and numbers follow the number rule.
    fn debug_text<T: Debug + ?Sized>(&self, leaf: &str, actual: &T, expected: &T) -> Outcome {
        let (actual_text, expected_text) = (format!("{actual:#?}"), format!("{expected:#?}"));
        if actual_text == expected_text {
            return Ok(());
        }
        let mut comparator = fuzzy::FuzzyStringComparator::new();
        comparator.set_log_destination(fuzzy::LogDestination::Buffer);
        comparator.set_verbose_level(1);
        comparator.set_acceptable_relative(self.tolerance.ratio);
        comparator.set_acceptable_absolute(self.tolerance.absdiff);
        if comparator.compare_strings(&actual_text, &expected_text) {
            return Ok(());
        }
        let log = String::from_utf8_lossy(comparator.log());
        let mut lines = log.lines();
        let reason = lines.next().unwrap_or("differs").to_owned();
        let line = log
            .lines()
            .find_map(|l| l.strip_prefix("  line:\t"))
            .unwrap_or("?")
            .to_owned();
        self.fail(
            leaf,
            format!("rendered metadata differs at Debug line {line}: {reason}"),
        )
    }

    fn meta(&mut self, leaf: &str, actual: &MetaInfo, expected: &MetaInfo) -> Outcome {
        if !actual.keys().eq(expected.keys()) {
            let only_actual: Vec<&String> = actual
                .keys()
                .filter(|k| !expected.contains_key(*k))
                .collect();
            let only_expected: Vec<&String> = expected
                .keys()
                .filter(|k| !actual.contains_key(*k))
                .collect();
            return self.fail(
                leaf,
                format!("key sets differ: only in actual {only_actual:?}, only in expected {only_expected:?}"),
            );
        }
        for ((key, a), e) in actual.iter().zip(expected.values()) {
            self.nested(format!("{leaf}[{key:?}]"), |walk| walk.meta_value(a, e))?;
        }
        Ok(())
    }

    fn meta_value(&mut self, actual: &MetaValue, expected: &MetaValue) -> Outcome {
        self.exact("unit", &actual.unit(), &expected.unit())?;
        match (actual.data(), expected.data()) {
            (MetaValueData::Empty, MetaValueData::Empty) => Ok(()),
            (MetaValueData::String(a), MetaValueData::String(e)) => self.exact("value", a, e),
            (MetaValueData::Integer(a), MetaValueData::Integer(e)) => {
                self.number("value", *a as f64, *e as f64)
            }
            (MetaValueData::Float(a), MetaValueData::Float(e)) => self.number("value", *a, *e),
            (MetaValueData::StringList(a), MetaValueData::StringList(e)) => {
                self.exact("value", a, e)
            }
            (MetaValueData::IntegerList(a), MetaValueData::IntegerList(e)) => {
                self.count("value", a.len(), e.len())?;
                for (index, (x, y)) in a.iter().zip(e).enumerate() {
                    self.number(&format!("value[{index}]"), *x as f64, *y as f64)?;
                }
                Ok(())
            }
            (MetaValueData::FloatList(a), MetaValueData::FloatList(e)) => {
                self.count("value", a.len(), e.len())?;
                for (index, (x, y)) in a.iter().zip(e).enumerate() {
                    self.number(&format!("value[{index}]"), *x, *y)?;
                }
                Ok(())
            }
            (a, e) => self.fail("value", format!("type {} vs {}", variant(a), variant(e))),
        }
    }

    fn processing_list(
        &mut self,
        leaf: &str,
        actual: &[&DataProcessing],
        expected: &[&DataProcessing],
    ) -> Outcome {
        self.count(leaf, actual.len(), expected.len())?;
        for (index, (a, e)) in actual.iter().zip(expected).enumerate() {
            self.nested(format!("{leaf}[{index}]"), |walk| {
                walk.exact("software.name", &a.software.name, &e.software.name)?;
                walk.exact("software.version", &a.software.version, &e.software.version)?;
                walk.debug_text(
                    "software.cv_terms",
                    &a.software.cv_terms,
                    &e.software.cv_terms,
                )?;
                walk.exact("actions", &a.actions, &e.actions)?;
                walk.exact("completion_time", &a.completion_time, &e.completion_time)?;
                walk.meta("metadata", &a.metadata, &e.metadata)
            })?;
        }
        Ok(())
    }

    fn feature(&mut self, actual: &Feature, expected: &Feature, depth: usize) -> Outcome {
        if depth > Feature::MAX_SUBORDINATE_DEPTH {
            return self.fail(
                "subordinates",
                "subordinate depth exceeds the checked limit".into(),
            );
        }
        if !self.ignore_unique_ids {
            self.exact("unique_id", &actual.unique_id, &expected.unique_id)?;
        }
        self.number("rt", actual.rt, expected.rt)?;
        self.number("mz", actual.mz, expected.mz)?;
        self.number(
            "intensity",
            f64::from(actual.intensity),
            f64::from(expected.intensity),
        )?;
        self.number(
            "quality",
            f64::from(actual.quality),
            f64::from(expected.quality),
        )?;
        self.exact("charge", &actual.charge, &expected.charge)?;
        self.number("width", f64::from(actual.width), f64::from(expected.width))?;
        self.number(
            "quality_rt",
            f64::from(actual.quality_rt),
            f64::from(expected.quality_rt),
        )?;
        self.number(
            "quality_mz",
            f64::from(actual.quality_mz),
            f64::from(expected.quality_mz),
        )?;
        self.meta("metadata", &actual.metadata, &expected.metadata)?;
        self.debug_text(
            "peptide_identifications",
            &actual.peptide_identifications,
            &expected.peptide_identifications,
        )?;
        self.debug_text("primary_id", &actual.primary_id, &expected.primary_id)?;
        self.debug_text("id_matches", &actual.id_matches, &expected.id_matches)?;
        self.count(
            "convex_hulls",
            actual.convex_hulls.len(),
            expected.convex_hulls.len(),
        )?;
        for (index, (a, e)) in actual
            .convex_hulls
            .iter()
            .zip(&expected.convex_hulls)
            .enumerate()
        {
            self.nested(format!("convex_hulls[{index}]"), |walk| walk.hull(a, e))?;
        }
        self.count(
            "subordinates",
            actual.subordinates.len(),
            expected.subordinates.len(),
        )?;
        for (index, (a, e)) in actual
            .subordinates
            .iter()
            .zip(&expected.subordinates)
            .enumerate()
        {
            self.nested(format!("subordinates[{index}]"), |walk| {
                walk.feature(a, e, depth + 1)
            })?;
        }
        Ok(())
    }

    fn hull(&mut self, actual: &ConvexHull2D, expected: &ConvexHull2D) -> Outcome {
        let (a, e) = (actual.hull_points(), expected.hull_points());
        self.count("points", a.len(), e.len())?;
        for (index, (p, q)) in a.iter().zip(&e).enumerate() {
            self.number(&format!("points[{index}].rt"), p.rt, q.rt)?;
            self.number(&format!("points[{index}].mz"), p.mz, q.mz)?;
        }
        Ok(())
    }

    fn spectrum(&mut self, actual: &MSSpectrum, expected: &MSSpectrum) -> Outcome {
        self.number("rt", actual.rt, expected.rt)?;
        self.exact("ms_level", &actual.ms_level, &expected.ms_level)?;
        if !self.ignore_unique_ids {
            self.exact("native_id", &actual.native_id, &expected.native_id)?;
        }
        self.exact("name", &actual.name, &expected.name)?;
        self.exact(
            "spectrum_type",
            &actual.spectrum_type,
            &expected.spectrum_type,
        )?;
        self.number("drift_time", actual.drift_time, expected.drift_time)?;
        self.exact(
            "drift_time_unit",
            &actual.drift_time_unit,
            &expected.drift_time_unit,
        )?;
        self.count("peaks", actual.peaks.len(), expected.peaks.len())?;
        for (index, (a, e)) in actual.peaks.iter().zip(&expected.peaks).enumerate() {
            self.number(&format!("peaks[{index}].mz"), a.mz, e.mz)?;
            self.number(
                &format!("peaks[{index}].intensity"),
                f64::from(a.intensity),
                f64::from(e.intensity),
            )?;
        }
        self.count(
            "precursors",
            actual.precursors.len(),
            expected.precursors.len(),
        )?;
        for (index, (a, e)) in actual
            .precursors
            .iter()
            .zip(&expected.precursors)
            .enumerate()
        {
            self.nested(format!("precursors[{index}]"), |walk| walk.precursor(a, e))?;
        }
        self.debug_text("products", &actual.products, &expected.products)?;
        self.arrays(
            &actual.float_data_arrays,
            &expected.float_data_arrays,
            &actual.integer_data_arrays,
            &expected.integer_data_arrays,
            &actual.string_data_arrays,
            &expected.string_data_arrays,
        )?;
        self.meta("metadata", &actual.metadata, &expected.metadata)?;
        self.processing_list(
            "data_processing",
            &actual
                .data_processing
                .iter()
                .map(AsRef::as_ref)
                .collect::<Vec<_>>(),
            &expected
                .data_processing
                .iter()
                .map(AsRef::as_ref)
                .collect::<Vec<_>>(),
        )?;
        self.debug_text(
            "instrument_settings",
            &actual.instrument_settings,
            &expected.instrument_settings,
        )?;
        self.debug_text(
            "acquisition_info",
            &actual.acquisition_info,
            &expected.acquisition_info,
        )?;
        self.debug_text("source_file", &actual.source_file, &expected.source_file)?;
        self.debug_text(
            "peptide_identifications",
            &actual.peptide_identifications,
            &expected.peptide_identifications,
        )
    }

    fn chromatogram(&mut self, actual: &MSChromatogram, expected: &MSChromatogram) -> Outcome {
        if !self.ignore_unique_ids {
            self.exact("native_id", &actual.native_id, &expected.native_id)?;
        }
        self.exact("name", &actual.name, &expected.name)?;
        self.exact(
            "chromatogram_type",
            &actual.chromatogram_type,
            &expected.chromatogram_type,
        )?;
        self.count("peaks", actual.peaks.len(), expected.peaks.len())?;
        for (index, (a, e)) in actual.peaks.iter().zip(&expected.peaks).enumerate() {
            self.number(&format!("peaks[{index}].rt"), a.rt, e.rt)?;
            self.number(
                &format!("peaks[{index}].intensity"),
                f64::from(a.intensity),
                f64::from(e.intensity),
            )?;
        }
        self.nested("precursor".to_owned(), |walk| {
            walk.precursor(&actual.precursor, &expected.precursor)
        })?;
        self.debug_text("product", &actual.product, &expected.product)?;
        self.arrays(
            &actual.float_data_arrays,
            &expected.float_data_arrays,
            &actual.integer_data_arrays,
            &expected.integer_data_arrays,
            &actual.string_data_arrays,
            &expected.string_data_arrays,
        )?;
        self.meta("metadata", &actual.metadata, &expected.metadata)?;
        self.processing_list(
            "data_processing",
            &actual
                .data_processing
                .iter()
                .map(AsRef::as_ref)
                .collect::<Vec<_>>(),
            &expected
                .data_processing
                .iter()
                .map(AsRef::as_ref)
                .collect::<Vec<_>>(),
        )?;
        self.debug_text(
            "instrument_settings",
            &actual.instrument_settings,
            &expected.instrument_settings,
        )?;
        self.debug_text(
            "acquisition_info",
            &actual.acquisition_info,
            &expected.acquisition_info,
        )?;
        self.debug_text("source_file", &actual.source_file, &expected.source_file)
    }

    fn precursor(&mut self, actual: &Precursor, expected: &Precursor) -> Outcome {
        self.number("mz", actual.mz, expected.mz)?;
        self.number(
            "intensity",
            f64::from(actual.intensity),
            f64::from(expected.intensity),
        )?;
        self.exact("charge", &actual.charge, &expected.charge)?;
        self.exact(
            "activation_methods",
            &actual.activation_methods,
            &expected.activation_methods,
        )?;
        self.number(
            "activation_energy",
            actual.activation_energy,
            expected.activation_energy,
        )?;
        self.number(
            "isolation_window_lower_offset",
            actual.isolation_window_lower_offset,
            expected.isolation_window_lower_offset,
        )?;
        self.number(
            "isolation_window_upper_offset",
            actual.isolation_window_upper_offset,
            expected.isolation_window_upper_offset,
        )?;
        self.optional_number(
            "isolation_target_mz",
            actual.isolation_target_mz,
            expected.isolation_target_mz,
        )?;
        self.optional_number("drift_time", actual.drift_time, expected.drift_time)?;
        self.exact(
            "drift_time_unit",
            &actual.drift_time_unit,
            &expected.drift_time_unit,
        )?;
        self.number(
            "drift_window_lower_offset",
            actual.drift_window_lower_offset,
            expected.drift_window_lower_offset,
        )?;
        self.number(
            "drift_window_upper_offset",
            actual.drift_window_upper_offset,
            expected.drift_window_upper_offset,
        )?;
        self.exact(
            "possible_charge_states",
            &actual.possible_charge_states,
            &expected.possible_charge_states,
        )?;
        self.debug_text("cv_terms", &actual.cv_terms, &expected.cv_terms)?;
        if !self.ignore_unique_ids {
            self.exact(
                "spectrum_reference",
                &actual.spectrum_reference,
                &expected.spectrum_reference,
            )?;
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn arrays(
        &mut self,
        float_actual: &[DataArray<f32>],
        float_expected: &[DataArray<f32>],
        integer_actual: &[DataArray<i32>],
        integer_expected: &[DataArray<i32>],
        string_actual: &[DataArray<String>],
        string_expected: &[DataArray<String>],
    ) -> Outcome {
        self.count(
            "float_data_arrays",
            float_actual.len(),
            float_expected.len(),
        )?;
        for (index, (a, e)) in float_actual.iter().zip(float_expected).enumerate() {
            self.nested(format!("float_data_arrays[{index}]"), |walk| {
                walk.array_header(a, e)?;
                for (i, (x, y)) in a.data.iter().zip(&e.data).enumerate() {
                    walk.number(&format!("data[{i}]"), f64::from(*x), f64::from(*y))?;
                }
                Ok(())
            })?;
        }
        self.count(
            "integer_data_arrays",
            integer_actual.len(),
            integer_expected.len(),
        )?;
        for (index, (a, e)) in integer_actual.iter().zip(integer_expected).enumerate() {
            self.nested(format!("integer_data_arrays[{index}]"), |walk| {
                walk.array_header(a, e)?;
                for (i, (x, y)) in a.data.iter().zip(&e.data).enumerate() {
                    walk.number(&format!("data[{i}]"), f64::from(*x), f64::from(*y))?;
                }
                Ok(())
            })?;
        }
        self.count(
            "string_data_arrays",
            string_actual.len(),
            string_expected.len(),
        )?;
        for (index, (a, e)) in string_actual.iter().zip(string_expected).enumerate() {
            self.nested(format!("string_data_arrays[{index}]"), |walk| {
                walk.array_header(a, e)?;
                walk.exact("data", &a.data, &e.data)
            })?;
        }
        Ok(())
    }

    fn array_header<T>(&mut self, actual: &DataArray<T>, expected: &DataArray<T>) -> Outcome {
        self.exact("name", &actual.name, &expected.name)?;
        self.count("data", actual.data.len(), expected.data.len())?;
        self.meta("metadata", &actual.metadata, &expected.metadata)?;
        self.processing_list(
            "data_processing",
            &actual
                .data_processing
                .iter()
                .map(AsRef::as_ref)
                .collect::<Vec<_>>(),
            &expected
                .data_processing
                .iter()
                .map(AsRef::as_ref)
                .collect::<Vec<_>>(),
        )
    }
}

fn variant(data: &MetaValueData) -> &'static str {
    match data {
        MetaValueData::Empty => "empty",
        MetaValueData::String(_) => "string",
        MetaValueData::Integer(_) => "integer",
        MetaValueData::Float(_) => "float",
        MetaValueData::StringList(_) => "string list",
        MetaValueData::IntegerList(_) => "integer list",
        MetaValueData::FloatList(_) => "float list",
    }
}
