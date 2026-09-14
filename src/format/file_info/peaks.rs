// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The FileInfo summary of peak files: DTA, DTA2D and mzML
//! (`FORMAT/FileInfo.cpp:1532-1965`, `:2005-2081`, `:2115-2127`, `:2384-2440`).
//!
//! The peak-file branch of the report: the experiment is loaded through
//! [`crate::format::FileHandler::load_experiment_with_options`] with default
//! [`crate::format::PeakFileOptions`] and the forced or detected type as the
//! only allowed type, as the source calls `FileHandler::loadExperiment(in, exp,
//! {in_type}, log_type, false, false)`. The branch then writes, in the source
//! order:
//!
//! 1. the instrument name and every mass analyzer with its resolution;
//! 2. the MS levels, the total peak count (spectrum peaks plus chromatogram
//!    points) and the spectrum count;
//! 3. the combined, overall-spectrum, per-MS-level and chromatogram range
//!    blocks, from the kernel range managers
//!    ([`crate::MSExperiment::spectrum_range_manager`],
//!    [`crate::MSExperiment::chromatogram_range_manager`]);
//! 4. the spectra per MS level; the peak type per MS level, the stored type of
//!    the first spectrum of that level (`getType(false)`) with the estimate of
//!    the first spectrum of that level holding more than ten peaks in
//!    parentheses, `Unknown` where no spectrum qualifies; the activation
//!    methods per `(level, method)`;
//! 5. the charge distribution of the first precursor of each spectrum;
//! 6. the data-array names with their occurrence counts, padded by the byte
//!    length of the longest name;
//! 7. the FAIMS compensation voltages, unconditionally computed;
//! 8. the chromatogram counts and the count per chromatogram type.
//!
//! `-m` adds the document, sample, instrument and contact metadata; `-p` the
//! data processing of the first spectrum; `-s` the MS1 intensity statistics
//! and one statistics block per data-array name. Integer arrays are promoted to
//! `double`, float arrays and intensities from `float`; a string array
//! contributes its name but no values, so its block is all zeros.
//!
//! The structured [`crate::format::file_info::model::PeakInfo`] and ranges are
//! filled alongside, as the source fills its `Result`.
//!
//! See `docs/FILE_INFO_SUPPORT.md` for the evidence and the native
//! differences.

use super::model::{FileInfoResult, Options, PeakInfo, Ranges};
use super::report::{
    ReportStream, range_set, statistics_buffer, summarize, write_charge_distribution,
    write_meta_title, write_processing, write_processing_title, write_ranges_text,
    write_ranges_tsv, write_statistics_title, write_summary_text,
};
use super::text_format::{WRITTEN_DIGITS_F32, list_to_string, to_str};
use crate::format::peak_type_estimator::PeakTypeEstimator;
use crate::format::{FileHandler, FileType, PeakFileOptions};
use crate::kernel::faims_helper::FaimsHelper;
use crate::kernel::ranges::RangeManager;
use crate::kernel::{MSExperiment, SpectrumType};
use crate::metadata::{ActivationMethod, ChromatogramType, DataProcessing, SpectrumSettings};
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::collections::btree_map::Entry;
use std::path::Path;

/// Minimum peak count, exclusive, of the spectrum whose peak type FileInfo
/// estimates for its MS level (`spectrum.size() > 10`, `FileInfo.cpp:1597`).
///
/// The estimator itself classifies from five peaks on; FileInfo asks for more
/// for a stable estimate.
pub const PEAK_TYPE_ESTIMATION_MIN_PEAKS: usize = 10;

/// Load the experiment and write the peak-file branch, `-m`, `-p` and `-s`.
pub(crate) fn report(
    path: &Path,
    in_type: FileType,
    options: &Options,
    os: &mut ReportStream,
    os_tsv: &mut ReportStream,
    result: &mut FileInfoResult,
) -> Result<()> {
    let experiment =
        FileHandler::load_experiment_with_options(path, &[in_type], &PeakFileOptions::default())?;
    let summary = Summary::compute(&experiment)?;

    write_content(&experiment, &summary, os, os_tsv, result)?;
    result
        .warnings
        .extend(summary.faims_warnings.iter().cloned());
    if options.meta {
        write_meta(&experiment, os, os_tsv);
    }
    if options.processing {
        write_title_and_processing(&experiment, os, os_tsv, result);
    }
    if options.statistics {
        write_statistics(&experiment, &summary, os)?;
    }
    Ok(())
}

/// Everything the content block and the structured result derive from the
/// experiment, computed once.
struct Summary {
    ranges: Ranges,
    ms_levels: Vec<u32>,
    total_peaks: u64,
    spectra_per_level: BTreeMap<u32, u64>,
    annotated_type: BTreeMap<u32, SpectrumType>,
    estimated_type: BTreeMap<u32, SpectrumType>,
    activation: BTreeMap<(u32, ActivationMethod), u64>,
    precursor_charges: BTreeMap<i32, u64>,
    float_arrays: BTreeMap<String, u64>,
    int_arrays: BTreeMap<String, u64>,
    string_arrays: BTreeMap<String, u64>,
    /// The source `meta_names`: float, integer and string array names together.
    meta_names: BTreeMap<String, u64>,
    faims_cvs: Vec<f64>,
    faims_warnings: Vec<String>,
    chromatogram_points: u64,
    chromatogram_types: BTreeMap<ChromatogramType, u64>,
}

impl Summary {
    fn compute(experiment: &MSExperiment) -> Result<Self> {
        let spectra = experiment.spectrum_range_manager()?;
        let chromatograms = experiment.chromatogram_range_manager()?;
        // MSExperiment::combinedRanges: the global spectrum ranges, then the
        // chromatogram ranges (as MSExperiment::combined_range_manager, without
        // computing the two again).
        let mut combined = RangeManager::experiment();
        combined.extend_unsafe(spectra.global());
        combined.extend_unsafe(&chromatograms);
        let mut ranges = Ranges {
            combined: range_set(&combined)?,
            spectra_overall: range_set(spectra.global())?,
            per_ms_level: BTreeMap::new(),
            chromatograms: range_set(&chromatograms)?,
            is_experiment: true,
        };
        for level in spectra.ms_levels() {
            if let Some(manager) = spectra.by_ms_level(level) {
                ranges.per_ms_level.insert(level, range_set(manager)?);
            }
        }

        let mut summary = Self {
            ranges,
            ms_levels: experiment.ms_levels(),
            total_peaks: 0,
            spectra_per_level: BTreeMap::new(),
            annotated_type: BTreeMap::new(),
            estimated_type: BTreeMap::new(),
            activation: BTreeMap::new(),
            precursor_charges: BTreeMap::new(),
            float_arrays: BTreeMap::new(),
            int_arrays: BTreeMap::new(),
            string_arrays: BTreeMap::new(),
            meta_names: BTreeMap::new(),
            faims_cvs: Vec::new(),
            faims_warnings: Vec::new(),
            chromatogram_points: 0,
            chromatogram_types: BTreeMap::new(),
        };

        for spectrum in &experiment.spectra {
            let level = spectrum.ms_level;
            *summary.spectra_per_level.entry(level).or_insert(0) += 1;
            summary.total_peaks = add_count(summary.total_peaks, spectrum.peaks.len())?;
            for precursor in &spectrum.precursors {
                for method in &precursor.activation_methods {
                    *summary.activation.entry((level, *method)).or_insert(0) += 1;
                }
            }
            if let Entry::Vacant(entry) = summary.annotated_type.entry(level) {
                entry.insert(spectrum.get_type(false)?);
            }
            if spectrum.peaks.len() > PEAK_TYPE_ESTIMATION_MIN_PEAKS {
                if let Entry::Vacant(entry) = summary.estimated_type.entry(level) {
                    entry.insert(PeakTypeEstimator::estimate_type(&spectrum.peaks)?);
                }
            }
            if let Some(precursor) = spectrum.precursors.first() {
                *summary
                    .precursor_charges
                    .entry(precursor.charge)
                    .or_insert(0) += 1;
            }
            for array in &spectrum.float_data_arrays {
                count_name(&mut summary.float_arrays, &array.name);
                count_name(&mut summary.meta_names, &array.name);
            }
            for array in &spectrum.integer_data_arrays {
                count_name(&mut summary.int_arrays, &array.name);
                count_name(&mut summary.meta_names, &array.name);
            }
            for array in &spectrum.string_data_arrays {
                count_name(&mut summary.string_arrays, &array.name);
                count_name(&mut summary.meta_names, &array.name);
            }
        }
        for chromatogram in &experiment.chromatograms {
            summary.chromatogram_points =
                add_count(summary.chromatogram_points, chromatogram.peaks.len())?;
            *summary
                .chromatogram_types
                .entry(chromatogram.chromatogram_type)
                .or_insert(0) += 1;
        }
        summary.total_peaks = summary
            .total_peaks
            .checked_add(summary.chromatogram_points)
            .ok_or_else(count_overflow)?;

        let voltages = FaimsHelper::get_compensation_voltages(experiment)?;
        summary.faims_cvs = voltages.values().collect();
        summary.faims_warnings = voltages.warnings;
        Ok(summary)
    }

    /// `NamesOfSpectrumType[level_annotated_picked[l]]`; the source map's
    /// `operator[]` inserts `0`, `Unknown`, for a level without an entry.
    fn annotated_name(&self, level: u32) -> &'static str {
        type_name(self.annotated_type.get(&level))
    }

    /// `NamesOfSpectrumType[level_estimated_picked[l]]`, `Unknown` when no
    /// spectrum of the level has more than ten peaks.
    fn estimated_name(&self, level: u32) -> &'static str {
        type_name(self.estimated_type.get(&level))
    }
}

fn type_name(spectrum_type: Option<&SpectrumType>) -> &'static str {
    SpectrumSettings::spectrum_type_to_string(
        spectrum_type.copied().unwrap_or(SpectrumType::Unknown),
    )
}

fn count_name(counts: &mut BTreeMap<String, u64>, name: &str) {
    match counts.get_mut(name) {
        Some(count) => *count += 1,
        None => {
            counts.insert(name.to_owned(), 1);
        }
    }
}

fn add_count(total: u64, count: usize) -> Result<u64> {
    u64::try_from(count)
        .ok()
        .and_then(|count| total.checked_add(count))
        .ok_or_else(count_overflow)
}

fn count_overflow() -> Error {
    Error::InvalidValue("FileInfo peak count overflows 64 bits".into())
}

fn to_int(level: u32) -> Result<i32> {
    i32::try_from(level).map_err(|_| {
        Error::InvalidValue(format!(
            "MS level {level} does not fit the FileInfo result's Int key"
        ))
    })
}

fn write_content(
    experiment: &MSExperiment,
    summary: &Summary,
    os: &mut ReportStream,
    os_tsv: &mut ReportStream,
    result: &mut FileInfoResult,
) -> Result<()> {
    let instrument = &experiment.settings.instrument;
    os.text("\nInstrument: ").text(&instrument.name).text("\n");
    for analyzer in &instrument.mass_analyzers {
        os.text("  Mass Analyzer: ")
            .text(analyzer.analyzer_type.name())
            .text(" (resolution: ")
            .double(analyzer.resolution)
            .text(")\n");
    }
    os.text("\n");

    os.text("MS levels: ");
    for (index, level) in summary.ms_levels.iter().enumerate() {
        if index > 0 {
            os.text(", ");
        }
        os.value(level);
    }
    os.text("\n");
    os.text("Total number of peaks: ")
        .value(summary.total_peaks)
        .text("\n");
    os.text("Number of spectra: ")
        .value(experiment.spectra.len())
        .text("\n\n");
    os_tsv
        .text("number of spectra\t")
        .value(experiment.spectra.len())
        .text("\ntotal number of peaks\t")
        .value(summary.total_peaks)
        .text("\n");

    write_experiment_ranges(os, os_tsv, &summary.ranges);

    if !summary.spectra_per_level.is_empty() {
        os.text("Number of spectra per MS level:\n");
        for (level, count) in &summary.spectra_per_level {
            os.text("  level ")
                .value(level)
                .text(": ")
                .value(count)
                .text("\n");
            os_tsv
                .text("number of MS")
                .value(level)
                .text(" spectra\t")
                .value(count)
                .text("\n");
        }
        os.text("\n");
    }

    os.text("Peak type from metadata (or estimated from data)\n");
    for &level in &summary.ms_levels {
        let annotated = summary.annotated_name(level);
        let estimated = summary.estimated_name(level);
        os.text("  level ")
            .value(level)
            .text(": ")
            .text(annotated)
            .text(" (")
            .text(estimated)
            .text(")\n");
        os_tsv
            .text("peak type metadata [annotation, estimate]\t")
            .text(annotated)
            .text("\t")
            .text(estimated)
            .text("\n");
    }
    os.text("\n");

    os.text("Activation methods\n");
    for ((level, method), count) in &summary.activation {
        os.text("    MS-Level ")
            .value(level)
            .text(" & ")
            .text(method.short_name())
            .text(" (")
            .text(method.name())
            .text("): ")
            .value(count)
            .text("\n");
        os_tsv
            .text("activation methods (mslevel, method, count)\t")
            .value(level)
            .text("\t")
            .text(method.short_name())
            .text("\t")
            .value(count)
            .text("\n");
    }
    os.text("\n");

    result.ranges = summary.ranges.clone();
    result.peak = Some(peak_info(experiment, summary)?);

    write_charge_distribution(os, os_tsv, "Precursor charge", &summary.precursor_charges);

    if !summary.meta_names.is_empty() {
        let width = summary
            .meta_names
            .keys()
            .map(String::len)
            .max()
            .unwrap_or(0);
        os.text("Meta data array:\n");
        for (name, count) in &summary.meta_names {
            os.text("  ")
                .text(name)
                .text(": ")
                .text(&" ".repeat(width - name.len()))
                .value(count)
                .text(" spectra\n");
        }
        os.text("\n");
    }

    if !summary.faims_cvs.is_empty() {
        let voltages: Vec<String> = summary.faims_cvs.iter().map(|cv| to_str(*cv)).collect();
        os.text("IM (FAIMS_CV): ")
            .text(&list_to_string(&voltages)?)
            .text("\n\n");
    }

    if !experiment.chromatograms.is_empty() {
        os.text("Number of chromatograms: ")
            .value(experiment.chromatograms.len())
            .text("\n");
        os_tsv
            .text("number of chromatograms\t")
            .value(experiment.chromatograms.len())
            .text("\n");
        os.text("Number of chromatographic peaks: ")
            .value(summary.chromatogram_points)
            .text("\n\n");
        os_tsv
            .text("number of chromatographic peaks\t")
            .value(summary.chromatogram_points)
            .text("\n");
        os.text("Number of chromatograms per type: \n");
        for (kind, count) in &summary.chromatogram_types {
            os.text("  ")
                .text(kind.name())
                .text(":                         ")
                .value(count)
                .text("\n");
        }
    }
    Ok(())
}

fn write_experiment_ranges(os: &mut ReportStream, os_tsv: &mut ReportStream, ranges: &Ranges) {
    os.text("Combined Ranges (spectra + chromatograms):\n");
    write_ranges_text(os, &ranges.combined, true);
    os.text("Spectrum Ranges:\n");
    write_ranges_text(os, &ranges.spectra_overall, true);
    for (level, set) in &ranges.per_ms_level {
        os.text("MS Level ").value(level).text(" Ranges:\n");
        write_ranges_text(os, set, true);
    }
    os.text("Chromatogram Ranges:\n");
    write_ranges_text(os, &ranges.chromatograms, false);

    write_ranges_tsv(os_tsv, "general: combined ranges: ", &ranges.combined);
    write_ranges_tsv(
        os_tsv,
        "general: spectrum ranges: ",
        &ranges.spectra_overall,
    );
    for (level, set) in &ranges.per_ms_level {
        write_ranges_tsv(os_tsv, &format!("general: MS{level} ranges: "), set);
    }
    // The chromatogram range manager has no mobility dimension, so its set
    // never holds a mobility range and no ion-mobility line is written.
    write_ranges_tsv(
        os_tsv,
        "general: chromatogram ranges: ",
        &ranges.chromatograms,
    );
}

fn peak_info(experiment: &MSExperiment, summary: &Summary) -> Result<PeakInfo> {
    let instrument = &experiment.settings.instrument;
    let mut info = PeakInfo {
        instrument_name: instrument.name.clone(),
        mass_analyzers: instrument
            .mass_analyzers
            .iter()
            .map(|analyzer| {
                (
                    analyzer.analyzer_type.name().to_owned(),
                    analyzer.resolution,
                )
            })
            .collect(),
        total_peaks: summary.total_peaks,
        num_spectra: u64::try_from(experiment.spectra.len()).map_err(|_| count_overflow())?,
        float_arrays: summary.float_arrays.clone(),
        int_arrays: summary.int_arrays.clone(),
        string_arrays: summary.string_arrays.clone(),
        precursor_charges: summary.precursor_charges.clone(),
        faims_cvs: summary.faims_cvs.clone(),
        num_chromatograms: u64::try_from(experiment.chromatograms.len())
            .map_err(|_| count_overflow())?,
        num_chrom_peaks: summary.chromatogram_points,
        ..PeakInfo::default()
    };
    for &level in &summary.ms_levels {
        let key = to_int(level)?;
        info.ms_levels.push(key);
        info.peak_type_per_ms_level.insert(
            key,
            format!(
                "{} ({})",
                summary.annotated_name(level),
                summary.estimated_name(level)
            ),
        );
    }
    for (level, count) in &summary.spectra_per_level {
        info.spectra_per_ms_level.insert(to_int(*level)?, *count);
    }
    for ((level, method), count) in &summary.activation {
        info.activation_methods
            .insert((to_int(*level)?, method.name().to_owned()), *count);
    }
    for (kind, count) in &summary.chromatogram_types {
        info.chromatogram_types
            .insert(kind.name().to_owned(), *count);
    }
    Ok(info)
}

fn write_meta(experiment: &MSExperiment, os: &mut ReportStream, os_tsv: &mut ReportStream) {
    let settings = &experiment.settings;
    let date = settings.date_time.get();
    write_meta_title(os);
    os.text("Document ID:        ")
        .text(&settings.document.identifier)
        .text("\nDate:               ")
        .text(&date)
        .text("\n");
    os_tsv
        .text("document id\t")
        .text(&settings.document.identifier)
        .text("\ndate\t")
        .text(&date)
        .text("\n");

    let sample = &settings.sample;
    os.text("\nSample:\n  name:             ")
        .text(&sample.name)
        .text("\n  organism:         ")
        .text(&sample.organism)
        .text("\n  comment:          ")
        .text(&sample.comment)
        .text("\n");
    os_tsv
        .text("sample name\t")
        .text(&sample.name)
        .text("\nsample organism\t")
        .text(&sample.organism)
        .text("\nsample comment\t")
        .text(&sample.comment)
        .text("\n");

    let instrument = &settings.instrument;
    os.text("\nInstrument:\n  name:             ")
        .text(&instrument.name)
        .text("\n  model:            ")
        .text(&instrument.model)
        .text("\n  vendor:           ")
        .text(&instrument.vendor)
        .text("\n  ion source(s):    ");
    os_tsv
        .text("instrument name\t")
        .text(&instrument.name)
        .text("\ninstrument model\t")
        .text(&instrument.model)
        .text("\ninstrument vendor\t")
        .text(&instrument.vendor)
        .text("\n");
    for (index, source) in instrument.ion_sources.iter().enumerate() {
        if index > 0 {
            os.text(", ");
        }
        os.text(source.ionization_method.name());
    }
    os.text("\n  mass analyzer(s): ");
    for (index, analyzer) in instrument.mass_analyzers.iter().enumerate() {
        if index > 0 {
            os.text(", ");
        }
        os.text(analyzer.analyzer_type.name());
    }
    os.text("\n  detector(s):      ");
    for (index, detector) in instrument.ion_detectors.iter().enumerate() {
        if index > 0 {
            os.text(", ");
        }
        os.text(detector.detector_type.name());
    }
    os.text("\n\n");

    for contact in &settings.contacts {
        os.text("Contact person:\n  first name:     ")
            .text(&contact.first_name)
            .text("\n  last name:      ")
            .text(&contact.last_name)
            .text("\n  email:          ")
            .text(&contact.email)
            .text("\n\n");
    }
}

fn write_title_and_processing(
    experiment: &MSExperiment,
    os: &mut ReportStream,
    os_tsv: &mut ReportStream,
    result: &mut FileInfoResult,
) {
    write_processing_title(os);
    let mut processing: Vec<&DataProcessing> = Vec::new();
    if let Some(first) = experiment.spectra.first() {
        os.text("Note: The data is taken from the first spectrum!\n\n");
        processing.extend(first.data_processing.iter().map(|step| step.as_ref()));
    }
    write_processing(os, os_tsv, &processing, result);
}

fn write_statistics(
    experiment: &MSExperiment,
    summary: &Summary,
    os: &mut ReportStream,
) -> Result<()> {
    write_statistics_title(os);

    let ms1_peaks = experiment
        .spectra
        .iter()
        .filter(|spectrum| spectrum.ms_level == 1)
        .try_fold(0usize, |total, spectrum| {
            total.checked_add(spectrum.peaks.len())
        })
        .ok_or_else(count_overflow)?;
    let mut intensities = statistics_buffer(ms1_peaks)?;
    for spectrum in experiment
        .spectra
        .iter()
        .filter(|spectrum| spectrum.ms_level == 1)
    {
        intensities.extend(spectrum.peaks.iter().map(|peak| f64::from(peak.intensity)));
    }
    let stats = summarize(&mut intensities)?;
    os.set_precision(WRITTEN_DIGITS_F32);
    os.text("Intensities:\n");
    write_summary_text(os, &stats);
    os.text("\n");

    // The source gathers, per name, the float arrays and then the integer
    // arrays of each spectrum in order; one pass per spectrum gives every
    // name the same sequence of values.
    let mut sizes: BTreeMap<&str, usize> = BTreeMap::new();
    for spectrum in &experiment.spectra {
        for array in &spectrum.float_data_arrays {
            grow(&mut sizes, &array.name, array.data.len())?;
        }
        for array in &spectrum.integer_data_arrays {
            grow(&mut sizes, &array.name, array.data.len())?;
        }
    }
    let mut values: BTreeMap<&str, Vec<f64>> = BTreeMap::new();
    for (name, size) in &sizes {
        values.insert(name, statistics_buffer(*size)?);
    }
    for spectrum in &experiment.spectra {
        let mut take = |name: &str, data: &mut dyn Iterator<Item = f64>| {
            if let Some(target) = values.get_mut(name) {
                target.extend(data);
            }
        };
        for array in &spectrum.float_data_arrays {
            take(
                &array.name,
                &mut array.data.iter().map(|value| f64::from(*value)),
            );
        }
        for array in &spectrum.integer_data_arrays {
            take(
                &array.name,
                &mut array.data.iter().map(|value| f64::from(*value)),
            );
        }
    }
    for name in summary.meta_names.keys() {
        let mut empty = Vec::new();
        let sample = values.get_mut(name.as_str()).unwrap_or(&mut empty);
        let stats = summarize(sample)?;
        os.text("Meta data: ").text(name).text("\n");
        write_summary_text(os, &stats);
        os.text("\n");
    }
    Ok(())
}

fn grow<'a>(sizes: &mut BTreeMap<&'a str, usize>, name: &'a str, count: usize) -> Result<()> {
    let size = sizes.entry(name).or_insert(0);
    *size = size.checked_add(count).ok_or_else(count_overflow)?;
    Ok(())
}
