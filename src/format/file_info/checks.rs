// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The FileInfo checks and listings behind `-i`, `-d` and `-c`
//! (`FORMAT/FileInfo.cpp:827-846`, `:1779-1795`, `:1799-1848`, `:1851-1964`).
//!
//! Three sections of the source's `report_` that inspect a file rather than
//! summarise it. Each writes into the human-readable report only; none of them
//! touches the TSV report, because the source writes none.
//!
//! - **`-i`, the indexed-mzML check** ([`write_index_check`]) runs directly
//!   after the general header and before the content of the file type. It is the
//!   only one of the three that fills a structured result:
//!   [`ValidationInfo::index_checked`], [`ValidationInfo::index_valid`] and,
//!   when the index parsed, the two record counts. A file whose index does not
//!   parse ends the whole report there, as the source `return`s
//!   (`FileInfo.cpp:844`): no content, no `-m`, no `-p`, no `-s`, and not even
//!   the two trailing newlines that `report_` writes at `:2443`. The FileInfo
//!   tool turns that state into its `ILLEGAL_PARAMETERS` exit code.
//! - **`-d`, the detailed listing**, is two blocks:
//!   [`write_detailed_chromatograms`] inside the chromatogram section, only when
//!   a selected-reaction-monitoring chromatogram is present, and
//!   [`write_detailed_spectra`] after it, only when the experiment holds a
//!   spectrum.
//! - **`-c`, the corrupt-data check** ([`write_corruption_check`]), last in the
//!   peak-file branch.
//!
//! `-d` and `-c` reach only the peak-file branch: the source guards both with
//! `options.detailed` / `options.check_corrupt` inside the branch that loaded an
//! experiment, so a featureXML map ignores them, and so does this port.
//!
//! # The result fields the source leaves empty
//!
//! `FileInfo.h:206-217` declares `CorruptionInfo` and `DetailInfo` with a
//! `performed` flag and pre-rendered message lines, but `report_` never assigns
//! to either at core `bc9cc12`: the `-d` and `-c` text goes only into the
//! stream. This port reproduces that, so
//! [`FileInfoResult::corruption`](crate::format::file_info::model::FileInfoResult::corruption)
//! and [`FileInfoResult::detail`](crate::format::file_info::model::FileInfoResult::detail)
//! stay at their defaults after a run that requested both flags. Filling them
//! would be a better API and a different one; `docs/FILE_INFO_SUPPORT.md`
//! records it as the divergence it would be, and `OpenMS_CPP_ISSUES.md` carries
//! the source defect.
//!
//! # Non-finite coordinates
//!
//! `-c` sorts the MS1 retention times and each spectrum's m/z values with
//! `std::sort`, whose comparator is not a strict weak ordering when a value is
//! NaN, so the source's behaviour there is undefined. This port refuses a NaN
//! retention time or m/z with [`Error::InvalidValue`] before writing anything,
//! rather than reproducing an undefined order. No loader on the peak-file branch
//! produces one — each validates its coordinates — so the refusal is
//! unreachable through [`FileInfo::run`](crate::format::file_info::report::FileInfo::run).
//! Infinities are left alone: they order and compare exactly as the source's
//! `<` and `==` do.
//!
//! [`ValidationInfo::index_checked`]: crate::format::file_info::model::ValidationInfo::index_checked
//! [`ValidationInfo::index_valid`]: crate::format::file_info::model::ValidationInfo::index_valid
//! [`Error::InvalidValue`]: crate::Error::InvalidValue

use super::model::FileInfoResult;
use super::report::ReportStream;
use crate::Error;
use crate::Result;
use crate::format::indexed_mzml::IndexedMzMLDecoder;
use crate::kernel::MSExperiment;
use crate::metadata::{ChromatogramType, DriftTimeUnit};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

/// The comment of a chromatogram, which the `-d` transition line prints last.
///
/// `ChromatogramSettings::getComment()` has no counterpart in
/// [`MSChromatogram`](crate::kernel::MSChromatogram): nothing in the source
/// calls `setComment` on a chromatogram, so every chromatogram FileInfo can see
/// — loaded from mzML, or produced by `convertSpectraToChromatograms` — carries
/// the empty default. The transition line therefore ends in the two spaces that
/// separate an empty name from an empty comment, exactly as the source writes
/// it. `docs/CHROMATOGRAM_MERGE_SUPPORT.md` records the unported member.
const CHROMATOGRAM_COMMENT: &str = "";

/// Run the `-i` indexed-mzML check and report it, the source
/// `FileInfo.cpp:827-846`.
///
/// `name` is the file name as the caller was given it, which the failure line
/// repeats verbatim; `path` is the same name as a path. Returns whether the
/// report continues: `false` reproduces the source's `return` after a failed
/// check, which ends the report with the failure text and nothing after it.
///
/// The check is the source's `Internal::IndexedMzMLHandler::openFile` reduced to
/// what `parseFooter_` records: the footer's `indexListOffset` is located, the
/// index sections at that offset are decoded, and success is
/// `findIndexListOffset() != -1 && parseOffsets() == 0`
/// (`IndexedMzMLHandler.cpp:20-62`). The counts are the lengths of the two
/// offset vectors. Nothing is read from the records themselves, so a file whose
/// index is valid but whose spectra are not still passes, as in the source.
///
/// [`IndexedMzMLHandler`](crate::format::indexed_mzml_handler::IndexedMzMLHandler)
/// is deliberately not used here: it goes on to locate the record list tags and
/// read the document header, and refuses files the source's index check accepts.
///
/// The source applies the check to whatever `-in` names, of any type; only the
/// FileInfo tool restricts `-i` to mzML, before the library runs.
///
/// # Errors
///
/// Returns [`Error::Io`] when the file cannot be opened or read, which the
/// source raises as `Exception::FileNotFound`, `FileNotReadable` or
/// `IOException`, and [`Error::Parse`] when the footer's offset is not a
/// non-negative 63-bit integer, which the source raises as
/// `Exception::ConversionError` (`IndexedMzMLDecoder.cpp:51-61`,
/// `:79-87`, `:152-160`). Every other decoding failure is *not* an error: it is
/// the invalid-index outcome, because `parseOffsets` reports all of its own
/// failures — an offset outside the file, a failed allocation and a malformed
/// index — by returning `-1` (`:97-99`, `:115-118`, `:332`). The port's index
/// byte and offset-count ceilings land in that same bucket, which is where the
/// source puts an index too large to hold.
pub(crate) fn write_index_check(
    path: &Path,
    name: &str,
    os: &mut ReportStream,
    result: &mut FileInfoResult,
) -> Result<bool> {
    result.validation.index_checked = true;
    // `os << "..." << std::endl`: a newline, and a flush a string stream ignores.
    os.text("Checking mzML file for valid indices ... \n");
    match parse_index(path)? {
        Some((spectra, chromatograms)) => {
            result.validation.index_valid = true;
            result.validation.indexed_spectra = spectra;
            result.validation.indexed_chromatograms = chromatograms;
            os.text("Found a valid indexed mzML XML File with ")
                .value(spectra)
                .text(" spectra and ")
                .value(chromatograms)
                .text(" chromatograms.\n");
            Ok(true)
        }
        None => {
            result.validation.index_valid = false;
            os.text("Could not detect a valid index for the mzML file ")
                .text(name)
                .text("\nEither the index is not present or is not correct.\n");
            Ok(false)
        }
    }
}

/// The source `parseFooter_`: the record counts of a parsed index, or `None`
/// when `parsing_success_` would stay `false`.
fn parse_index(path: &Path) -> Result<Option<(u64, u64)>> {
    let decoder = IndexedMzMLDecoder::default();
    let Some(offset) = decoder.find_index_list_offset(path)? else {
        return Ok(None);
    };
    match decoder.parse_offsets(path, offset) {
        Ok(offsets) => Ok(Some((
            offsets.spectra.len() as u64,
            offsets.chromatograms.len() as u64,
        ))),
        // The source throws only for the file system here; everything else is
        // its `return -1`.
        Err(Error::Io(error)) => Err(Error::Io(error)),
        Err(_) => Ok(None),
    }
}

/// Write the `-d` listing of selected-reaction-monitoring transitions, the
/// source `FileInfo.cpp:1779-1795`.
///
/// The caller has already written the chromatogram counts and is inside the
/// source's `if (!exp.getChromatograms().empty())`. `chromatogram_types` is the
/// per-type histogram that section built; the listing is written only when it
/// holds [`ChromatogramType::SelectedReactionMonitoring`], as the source's
/// `chrom_types.contains(...)` requires, and then covers every chromatogram of
/// that type in storage order.
///
/// One line per transition: the precursor m/z, the product m/z, the retention
/// time of the first and of the last point, the name and the comment, separated
/// by single spaces. All four numbers are written at the report's current stream
/// precision, which no section before this one changes.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] for a selected-reaction-monitoring
/// chromatogram with no points. The source reads `ms.front()` and `ms.back()`
/// without checking, which is undefined on an empty container, so there is no
/// behaviour to reproduce. mzML carries no such chromatogram — the reader gives
/// every chromatogram its points — and
/// `ChromatogramTools::convert_spectra_to_chromatograms` builds one point per
/// source spectrum, so an SRM chromatogram on this path always has at least one.
pub(crate) fn write_detailed_chromatograms(
    experiment: &MSExperiment,
    chromatogram_types: &BTreeMap<ChromatogramType, u64>,
    os: &mut ReportStream,
) -> Result<()> {
    if !chromatogram_types.contains_key(&ChromatogramType::SelectedReactionMonitoring) {
        return Ok(());
    }
    os.text("\n -- Detailed chromatogram listing -- \n");
    os.text("\nSelected Reaction Monitoring Transitions:\n");
    os.text("Q1 Q3 RT_begin RT_end name comment\n");
    for chromatogram in &experiment.chromatograms {
        if chromatogram.chromatogram_type != ChromatogramType::SelectedReactionMonitoring {
            continue;
        }
        let (Some(first), Some(last)) = (chromatogram.peaks.first(), chromatogram.peaks.last())
        else {
            return Err(Error::InvalidValue(format!(
                "FileInfo -d: selected reaction monitoring chromatogram '{}' has no points; \
                 the source reads its first and last point unchecked",
                chromatogram.native_id
            )));
        };
        os.double(chromatogram.precursor.mz)
            .text(" ")
            .double(chromatogram.product.mz)
            .text(" ")
            .double(first.rt)
            .text(" ")
            .double(last.rt)
            .text(" ")
            .text(&chromatogram.name)
            .text(" ")
            .text(CHROMATOGRAM_COMMENT)
            .text("\n");
    }
    Ok(())
}

/// Write the `-d` per-spectrum listing, the source `FileInfo.cpp:1799-1848`.
///
/// Nothing is written for an experiment with no spectrum, as the source's
/// `!exp.empty()` requires; `empty()` asks about the spectra alone, so a file of
/// chromatograms only reaches the transition listing above and no further.
///
/// Each spectrum contributes its one-based number, MS level, scan mode, peak
/// count, retention time and m/z extent, then its precursors. The m/z extent is
/// the m/z of the **first and last stored** peak, not the smallest and largest:
/// the source reads `begin()` and `rbegin()`, which differ from the extremes in
/// a spectrum whose peaks are not sorted, as a DTA or DTA2D scan need not be.
/// An empty spectrum writes no extent at all and leaves the `m/z:` line
/// unterminated, so the following `Precursors:` continues it — the source writes
/// its newline only inside the `if (!spectrum.empty())`, and this port keeps
/// that line intact rather than tidying it.
///
/// An ion-mobility line follows only when the spectrum carries a drift time
/// unit, and names the unit as `IMTypes`' `NamesOfDriftTimeUnit` spells it.
/// Each precursor writes its charge, m/z and activation methods, the methods in
/// the source's enumeration order, because both the source's `std::set` and
/// [`BTreeSet`](std::collections::BTreeSet) order by the enumerator. A blank
/// line closes every precursor, including the last.
///
/// The source counts spectra and precursors in `UInt`, which wraps above
/// 4294967295; this port counts in `usize`, which cannot reach that value for a
/// loaded experiment on any supported target and therefore never wraps.
pub(crate) fn write_detailed_spectra(experiment: &MSExperiment, os: &mut ReportStream) {
    if experiment.spectra.is_empty() {
        return;
    }
    os.text("\n-- Detailed spectrum listing --\n");
    for (index, spectrum) in experiment.spectra.iter().enumerate() {
        os.text("\nSpectrum ")
            .value(index + 1)
            .text(":\n  mslevel:    ")
            .value(spectrum.ms_level)
            .text("\n  scanMode:   ")
            .text(spectrum.instrument_settings.scan_mode.name())
            .text("\n  peaks:      ")
            .value(spectrum.peaks.len())
            .text("\n  RT:         ")
            .double(spectrum.rt)
            .text("\n  m/z:        ");
        if let (Some(first), Some(last)) = (spectrum.peaks.first(), spectrum.peaks.last()) {
            os.double(first.mz).text(" .. ").double(last.mz).text("\n");
        }
        if spectrum.drift_time_unit != DriftTimeUnit::None {
            os.text("  IM:         ")
                .double(spectrum.drift_time)
                .text(" ")
                .text(spectrum.drift_time_unit.name())
                .text("\n");
        }
        os.text("Precursors:  ")
            .value(spectrum.precursors.len())
            .text("\n");
        for (precursor_index, precursor) in spectrum.precursors.iter().enumerate() {
            os.text("Precursor[")
                .value(precursor_index)
                .text("]\n  charge: ")
                .value(precursor.charge)
                .text("\n  mz:     ")
                .double(precursor.mz)
                .text("\n  activation methods: \n");
            for method in &precursor.activation_methods {
                os.text("    ")
                    .text(method.short_name())
                    .text(" (")
                    .text(method.name())
                    .text(")\n");
            }
            os.text("\n");
        }
    }
}

/// Write the `-c` corrupt-data check, the source `FileInfo.cpp:1851-1964`.
///
/// The header is written whatever the outcome, so a clean file produces the
/// header and nothing else. The source's order is kept exactly, because the
/// lines are the report:
///
/// 1. one line when the spectra are not in ascending retention time
///    (`MSExperiment::isSorted(false)`, which looks at the retention times
///    alone);
/// 2. per spectrum, in storage order: MS level zero, no peaks, and every
///    repetition of a data-array name. The three array kinds share one name set,
///    as the source's single `std::map` does, so a float array and an integer
///    array of the same name collide;
/// 3. one line per repeated MS1 retention time, over the sorted retention times
///    of the MS1 spectra only;
/// 4. per spectrum again, in storage order: unsorted m/z, then every negative
///    peak intensity in storage order, then every repeated m/z over that
///    spectrum's sorted m/z values.
///
/// A repeated value is reported once per repetition, so a value stored three
/// times gives two lines. Intensities are `float` and are written as the source
/// writes them, promoted to `double` at the report's stream precision.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] for a NaN retention time or m/z, which would
/// enter a `std::sort` in the source and leave its behaviour undefined (see the
/// module documentation), and when the retention-time or m/z buffer cannot be
/// allocated. Both are checked before the header is written, so a refusal leaves
/// the report untouched.
pub(crate) fn write_corruption_check(
    experiment: &MSExperiment,
    os: &mut ReportStream,
) -> Result<()> {
    preflight(experiment)?;
    os.text("\n-- Checking for corrupt data --\n\n");

    // `exp.isSorted(false)`: the retention times only, the source's "// TODO CHROM".
    if !experiment.is_sorted(false) {
        os.text("Error: Spectrum retention times are not sorted in ascending order\n");
    }

    let mut ms1_rts: Vec<f64> = Vec::new();
    ms1_rts
        .try_reserve(experiment.spectra.len())
        .map_err(|_| allocation("retention times"))?;
    for spectrum in &experiment.spectra {
        if spectrum.ms_level == 0 {
            os.text("Error: MS-level 0 in spectrum (RT: ")
                .double(spectrum.rt)
                .text(")\n");
        }
        if spectrum.peaks.is_empty() {
            os.text("Warning: No peaks in spectrum (RT: ")
                .double(spectrum.rt)
                .text(")\n");
        }
        let mut names = BTreeSet::new();
        let float_names = spectrum.float_data_arrays.iter().map(|array| &array.name);
        let integer_names = spectrum.integer_data_arrays.iter().map(|array| &array.name);
        let string_names = spectrum.string_data_arrays.iter().map(|array| &array.name);
        for name in float_names.chain(integer_names).chain(string_names) {
            if !names.insert(name.as_str()) {
                os.text("Error: Duplicate meta data array name '")
                    .text(name)
                    .text("' in spectrum (RT: ")
                    .double(spectrum.rt)
                    .text(")\n");
            }
        }
        if spectrum.ms_level == 1 {
            ms1_rts.push(spectrum.rt);
        }
    }

    ms1_rts.sort_by(f64::total_cmp);
    for pair in ms1_rts.windows(2) {
        if pair[0] == pair[1] {
            os.text("Error: Duplicate spectrum retention time: ")
                .double(pair[1])
                .text("\n");
        }
    }

    let mut mzs: Vec<f64> = Vec::new();
    for spectrum in &experiment.spectra {
        if !spectrum.is_sorted() {
            os.text("Error: Peak m/z positions are not sorted in ascending order in spectrum (RT: ")
                .double(spectrum.rt)
                .text(")\n");
        }
        mzs.clear();
        mzs.try_reserve(spectrum.peaks.len())
            .map_err(|_| allocation("m/z values"))?;
        for peak in &spectrum.peaks {
            let intensity = f64::from(peak.intensity);
            if intensity < 0.0 {
                os.text("Warning: Negative peak intensity peak (RT: ")
                    .double(spectrum.rt)
                    .text(" MZ: ")
                    .double(peak.mz)
                    .text(" intensity: ")
                    .double(intensity)
                    .text(")\n");
            }
            mzs.push(peak.mz);
        }
        mzs.sort_by(f64::total_cmp);
        for pair in mzs.windows(2) {
            if pair[0] == pair[1] {
                os.text("Error: Duplicate peak m/z ")
                    .double(pair[1])
                    .text(" in spectrum (RT: ")
                    .double(spectrum.rt)
                    .text(")\n");
            }
        }
    }
    Ok(())
}

/// Refuse the NaN coordinates the source's two `std::sort` calls leave
/// undefined, before any of the check is written.
fn preflight(experiment: &MSExperiment) -> Result<()> {
    for spectrum in &experiment.spectra {
        if spectrum.rt.is_nan() {
            return Err(Error::InvalidValue(
                "FileInfo -c: a spectrum retention time is NaN, which the source's std::sort of \
                 the MS1 retention times leaves undefined"
                    .into(),
            ));
        }
        if spectrum.peaks.iter().any(|peak| peak.mz.is_nan()) {
            return Err(Error::InvalidValue(
                "FileInfo -c: a peak m/z is NaN, which the source's std::sort of the m/z values \
                 leaves undefined"
                    .into(),
            ));
        }
    }
    Ok(())
}

fn allocation(what: &str) -> Error {
    Error::InvalidValue(format!("FileInfo -c: cannot allocate the {what} buffer"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::{MSChromatogram, MSSpectrum, Peak1D};

    /// The source reads `ms.front()` and `ms.back()` of every selected-reaction
    /// monitoring chromatogram unchecked, which is undefined on an empty one.
    /// No loader produces one, so this experiment is built by hand.
    #[test]
    fn an_empty_srm_chromatogram_is_refused() {
        let mut experiment = MSExperiment::default();
        experiment.chromatograms.push(MSChromatogram {
            chromatogram_type: ChromatogramType::SelectedReactionMonitoring,
            native_id: "empty_transition".into(),
            ..MSChromatogram::default()
        });
        let types = BTreeMap::from([(ChromatogramType::SelectedReactionMonitoring, 1u64)]);
        let mut os = ReportStream::new();
        let error = write_detailed_chromatograms(&experiment, &types, &mut os).unwrap_err();
        assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
        assert!(format!("{error}").contains("empty_transition"));
    }

    /// Without a selected-reaction-monitoring chromatogram nothing is written,
    /// so an empty chromatogram of another type is never read.
    #[test]
    fn a_listing_without_srm_writes_nothing() {
        let mut experiment = MSExperiment::default();
        experiment
            .chromatograms
            .push(MSChromatogram::default());
        let types = BTreeMap::from([(ChromatogramType::Mass, 1u64)]);
        let mut os = ReportStream::new();
        write_detailed_chromatograms(&experiment, &types, &mut os).unwrap();
        assert_eq!(os.into_string(), "");
    }

    /// A NaN coordinate would enter one of the source's two `std::sort` calls
    /// and leave its order undefined; the port refuses before writing anything.
    #[test]
    fn a_nan_retention_time_is_refused() {
        let mut experiment = MSExperiment::default();
        experiment.spectra.push(MSSpectrum {
            rt: f64::NAN,
            ..MSSpectrum::default()
        });
        let mut os = ReportStream::new();
        let error = write_corruption_check(&experiment, &mut os).unwrap_err();
        assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
        assert_eq!(os.into_string(), "", "the refusal leaves the report untouched");
    }

    #[test]
    fn a_nan_mz_is_refused() {
        let mut experiment = MSExperiment::default();
        experiment.spectra.push(MSSpectrum {
            rt: 1.0,
            peaks: vec![Peak1D::new(f64::NAN, 1.0)],
            ..MSSpectrum::default()
        });
        let mut os = ReportStream::new();
        let error = write_corruption_check(&experiment, &mut os).unwrap_err();
        assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    }

    /// An infinity is not NaN: the source's `<` and `==` are defined on it, so
    /// the check runs and reports the duplicate.
    #[test]
    fn an_infinite_mz_is_checked_like_any_other() {
        let mut experiment = MSExperiment::default();
        experiment.spectra.push(MSSpectrum {
            rt: 1.0,
            peaks: vec![
                Peak1D::new(f64::INFINITY, 1.0),
                Peak1D::new(f64::INFINITY, 2.0),
            ],
            ..MSSpectrum::default()
        });
        let mut os = ReportStream::new();
        write_corruption_check(&experiment, &mut os).unwrap();
        let text = os.into_string();
        assert!(
            text.contains("Error: Duplicate peak m/z inf in spectrum (RT: 1)\n"),
            "{text}"
        );
    }

    /// A value stored three times gives two lines, as the source's pairwise
    /// comparison over the sorted values does.
    #[test]
    fn a_triplicate_mz_reports_twice() {
        let mut experiment = MSExperiment::default();
        experiment.spectra.push(MSSpectrum {
            rt: 2.0,
            peaks: vec![
                Peak1D::new(5.0, 1.0),
                Peak1D::new(5.0, 2.0),
                Peak1D::new(5.0, 3.0),
            ],
            ..MSSpectrum::default()
        });
        let mut os = ReportStream::new();
        write_corruption_check(&experiment, &mut os).unwrap();
        assert_eq!(
            os.into_string()
                .matches("Error: Duplicate peak m/z 5 in spectrum (RT: 2)\n")
                .count(),
            2
        );
    }

    /// A repeated data-array name, which this port's mzML reader refuses before
    /// `-c` can see it (`src/format/mzml.rs:1038`), so no file reaches this
    /// line. The C++ reader loads such a file and the Release oracle's
    /// `c_arrays` and `c_arrays_mixed` cases record what it writes; this is the
    /// same rendering, from an experiment built in memory. The three array
    /// kinds share one name set, as the source's single `std::map` does.
    #[test]
    fn a_repeated_data_array_name_is_reported_across_array_kinds() {
        use crate::kernel::DataArray;
        let mut spectrum = MSSpectrum {
            rt: 2.0,
            peaks: vec![Peak1D::new(100.0, 1.0)],
            ..MSSpectrum::default()
        };
        spectrum.float_data_arrays.push(DataArray {
            name: "shared".into(),
            ..DataArray::default()
        });
        spectrum.integer_data_arrays.push(DataArray {
            name: "shared".into(),
            ..DataArray::default()
        });
        spectrum.string_data_arrays.push(DataArray {
            name: "shared".into(),
            ..DataArray::default()
        });
        let mut experiment = MSExperiment::default();
        experiment.spectra.push(spectrum);
        let mut os = ReportStream::new();
        write_corruption_check(&experiment, &mut os).unwrap();
        assert_eq!(
            os.into_string()
                .matches("Error: Duplicate meta data array name 'shared' in spectrum (RT: 2)\n")
                .count(),
            2,
            "the first array claims the name; the other two repeat it"
        );
    }

    /// TOPP_FileInfo_12's input: its index parses with the counts the C++
    /// reports, although the strict mzML reader refuses its content.
    #[test]
    fn the_upstream_test_12_index_parses() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/file_info/inputs/FileInfo_12_input.mzML");
        let mut result = FileInfoResult::default();
        let mut os = ReportStream::new();
        assert!(write_index_check(&path, "in.mzML", &mut os, &mut result).unwrap());
        assert!(result.validation.index_valid);
        assert_eq!(result.validation.indexed_spectra, 3);
        assert_eq!(result.validation.indexed_chromatograms, 0);
        assert_eq!(
            os.into_string(),
            "Checking mzML file for valid indices ... \n\
             Found a valid indexed mzML XML File with 3 spectra and 0 chromatograms.\n"
        );
    }

    /// A file with no footer offset: the failure text names the file as it was
    /// given, and nothing follows it.
    #[test]
    fn a_file_without_an_index_reports_the_failure() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/file_info/inputs/empty.mzML");
        let mut result = FileInfoResult::default();
        let mut os = ReportStream::new();
        assert!(!write_index_check(&path, "given/name.mzML", &mut os, &mut result).unwrap());
        assert!(result.validation.index_checked);
        assert!(!result.validation.index_valid);
        assert_eq!(
            os.into_string(),
            "Checking mzML file for valid indices ... \n\
             Could not detect a valid index for the mzML file given/name.mzML\n\
             Either the index is not present or is not correct.\n"
        );
    }

    /// `empty()` asks about the spectra alone, so an experiment of
    /// chromatograms only writes no per-spectrum listing.
    #[test]
    fn no_spectrum_listing_without_spectra() {
        let mut experiment = MSExperiment::default();
        experiment
            .chromatograms
            .push(MSChromatogram::default());
        let mut os = ReportStream::new();
        write_detailed_spectra(&experiment, &mut os);
        assert_eq!(os.into_string(), "");
    }
}
