// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The FileInfo report: section order and text and TSV rendering
//! (`FORMAT/FileInfo.h`, `FORMAT/FileInfo.cpp`).
//!
//! [`FileInfo`] is the library-level equivalent of the FileInfo tool. A run
//! loads a file, extracts its file-level information into a
//! [`FileInfoResult`] and renders, in the same pass, the exact text and TSV
//! reports the tool writes to `-out` and `-out_tsv`. As in the source,
//! computation and rendering happen together (`FileInfo::report_`), the two
//! reports are cached in the result, and [`FileInfo::to_text`] and
//! [`FileInfo::to_tsv`] only return the cache.
//!
//! Report order, as `report_` writes it: the general header, `-v`, `-i`, the
//! content of the file type (which ends with `-d` and `-c` on a peak file),
//! `-m`, `-p`, `-s`, and two trailing newlines. A `-v` or `-i` failure returns
//! early and writes none of what follows it, the two trailing newlines
//! included.
//!
//! # Scope
//!
//! This port runs the peak-file branch for DTA, DTA2D and mzML
//! ([`crate::format::file_info::peaks`]) and the featureXML branch
//! ([`crate::format::file_info::features`]), each with `-m`, `-p` and `-s`, and
//! the `-i`, `-d` and `-c` checks ([`crate::format::file_info::checks`]): `-i`
//! before the content of any type, `-d` and `-c` inside the peak-file branch,
//! which is where the source guards them, so a featureXML map ignores both as
//! the source does.
//!
//! Every other part of the source report is refused with
//! [`Error::Unsupported`] naming the branch, once the type is known and before
//! the file is loaded, so a run never returns a partial report. The type is
//! known without touching the file when it is forced or recognised from the
//! name; otherwise the type detection reads the start of the file first, and
//! its I/O error comes before the refusal. A branch refusal comes after the
//! `-i` check, as the source's order has it, so an unparsable index is reported
//! in full even on a branch this port does not run. Refused are:
//!
//! - `-v` (schema and semantic validation), for every type;
//! - the consensusXML, idXML, mzIdentML, FASTA, pepXML, mzTab, trafoXML and PQP
//!   branches;
//! - peak files of the types the source loads but no native loader serves on
//!   this path: mzXML, mzData, MGF, MS2, sqMass, XMass (`fid`) and MSP, and
//!   Thermo RAW and Bruker TDF, which the source loads when built with its
//!   default `WITH_THERMO_RAW` and `WITH_OPENTIMS` options.
//!
//! `docs/FILE_INFO_SUPPORT.md` holds the API mapping, the preserved source
//! conventions, the native differences and the evidence.
//!
//! # Stream state
//!
//! The source writes into two `std::ostringstream`s. Their precision starts at
//! [`DEFAULT_STREAM_PRECISION`] and a statistics block changes it, after which
//! it stays in force (`FileInfo.cpp:2224-2255` and `:2404`). The port tracks
//! the precision of each report explicitly; every `double` the source streams
//! is formatted with [`ostream_g`] at the tracked precision, every range bound
//! with [`fixed_truncated`], as `StringUtils::number` cuts rather than refuses.
//!
//! [`FileInfoResult`]: crate::format::file_info::model::FileInfoResult
//! [`FileInfo`]: crate::format::file_info::report::FileInfo
//! [`FileInfo::to_text`]: crate::format::file_info::report::FileInfo::to_text
//! [`FileInfo::to_tsv`]: crate::format::file_info::report::FileInfo::to_tsv
//! [`Error::Unsupported`]: crate::Error::Unsupported
//! [`DEFAULT_STREAM_PRECISION`]: crate::format::file_info::text_format::DEFAULT_STREAM_PRECISION
//! [`ostream_g`]: crate::format::file_info::text_format::ostream_g
//! [`fixed_truncated`]: crate::format::file_info::text_format::fixed_truncated

use super::model::{FileInfoResult, Options, ProcessingStep, Range, RangeSet};
use super::text_format::{DEFAULT_STREAM_PRECISION, fixed_truncated, ostream_g};
use crate::format::{FileHandler, FileType};
use crate::kernel::ranges::{MSDim, RangeManager};
use crate::math::statistic_functions::SummaryStatistics;
use crate::metadata::DataProcessing;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::fmt::{Display, Write as _};
use std::path::Path;

/// The library-level equivalent of the FileInfo tool, the source class
/// `OpenMS::FileInfo`.
///
/// The source class holds no state; this is a unit struct, [`Default`] is the
/// constructor and it has no drop glue, which is the (empty) destructor.
///
/// # Examples
///
/// ```no_run
/// use openms::format::file_info::report::FileInfo;
///
/// let result = FileInfo::new().run_all("data.featureXML")?;
/// if let Some(feature) = &result.feature {
///     println!("{}", feature.num_features);
/// }
/// print!("{}", FileInfo::to_text(&result)); // the FileInfo tool's -out text
/// # Ok::<(), openms::Error>(())
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FileInfo;

impl FileInfo {
    /// Values one summary-statistics block may collect, checked before the
    /// values are copied.
    ///
    /// Native ceiling; the source collects any number. At eight bytes a value
    /// the largest block needs at most 1 GiB.
    pub const MAX_STATISTICS_VALUES: usize = 1 << 27;

    /// A FileInfo instance, the source default constructor.
    pub const fn new() -> Self {
        Self
    }

    /// Load `filename` and compute everything `options` requests, the source
    /// `run(filename, options)`.
    ///
    /// The type is [`Options::forced_type`] unless that is
    /// [`FileType::Unknown`], in which case it is detected from the file name
    /// and then the content ([`FileHandler::get_type`]). A directory whose name
    /// gives no type is [`FileType::Unknown`], as the source's content check
    /// reads it as an empty file, where [`FileHandler::get_type`] returns the
    /// I/O error. An unknown type yields
    /// a result with only [`FileInfoResult::meta`] filled and empty reports,
    /// as the source returns and lets the caller report it. A forced type
    /// selects the branch and the single allowed type of the loader, which
    /// still detects the type itself: a featureXML map under a neutral `.tmp`
    /// name loads when featureXML is forced, but a file whose detected type
    /// differs from the forced one is refused.
    ///
    /// The file name is printed exactly as given, in the `File name:` line of
    /// the text and the `general: file name` line of the TSV.
    ///
    /// # Errors
    ///
    /// The source documents `Exception::FileNotFound` for a missing file and
    /// `Exception::ParseError` or format-specific exceptions for a file that
    /// cannot be read, and leaves them to the caller. Here:
    ///
    /// - [`Error::InvalidValue`] for a file name that is not UTF-8;
    /// - [`Error::Io`] when the type must be detected from content and the file
    ///   is not a directory and cannot be opened or read, or when the loader
    ///   cannot open it;
    /// - [`Error::Unsupported`] for a flag or branch this port does not run
    ///   (see the module documentation), before the file is loaded; before any
    ///   file access only when the type is forced or recognised from the name.
    ///   A branch refusal follows the `-i` check, as the source's order does,
    ///   so a file whose index does not parse yields the report of that failure
    ///   rather than the refusal;
    /// - [`Error::Parse`] for a type the source cannot load as a peak file
    ///   either, with the source message `type is not supported for loading
    ///   experiments`;
    /// - [`Error::InvalidValue`] for imzML, which the source refuses with
    ///   `Exception::InvalidFileType`, for a detected type other than a forced
    ///   one (the source's `ParseError`), and for a statistics block above
    ///   [`FileInfo::MAX_STATISTICS_VALUES`];
    /// - the loader's error for malformed input or exceeded reader limits, and
    ///   the kernel's error for values the range and type computations refuse.
    ///
    /// Nothing is returned on error; there is no partial report. The one report
    /// that is deliberately short is `-i`'s: an index that does not parse ends
    /// the report there and returns it, as the source `return`s from `report_`,
    /// and the FileInfo tool turns that into its `ILLEGAL_PARAMETERS` exit code
    /// by reading [`ValidationInfo::index_valid`](crate::format::file_info::model::ValidationInfo::index_valid).
    pub fn run(&self, filename: impl AsRef<Path>, options: &Options) -> Result<FileInfoResult> {
        let path = filename.as_ref();
        let name = path
            .to_str()
            .ok_or_else(|| Error::InvalidValue("FileInfo file name must be UTF-8".into()))?;
        let mut result = FileInfoResult::default();
        result.meta.file_name = name.to_owned();
        let in_type = if options.forced_type == FileType::Unknown {
            detect_type(path)?
        } else {
            options.forced_type
        };
        result.meta.file_type = in_type;
        result.meta.file_type_name = in_type.name().to_owned();
        if in_type == FileType::Unknown {
            return Ok(result);
        }
        check_flags_supported(options)?;

        let mut os = ReportStream::new();
        let mut os_tsv = ReportStream::new();
        os.text("\n-- General information --\n\nFile name: ")
            .text(name)
            .text("\nFile type: ")
            .text(in_type.name())
            .text("\n");
        os_tsv
            .text("general: file name\t")
            .text(name)
            .text("\ngeneral: file type\t")
            .text(in_type.name())
            .text("\n");

        // FileInfo.cpp:827-846: the index check sits after the general header and
        // before the content, for every type, and its failure ends the report
        // there. Only then is a branch this port does not run refused, so an
        // invalid index is reported in full even on such a branch.
        if options.check_index
            && !super::checks::write_index_check(path, name, &mut os, &mut result)?
        {
            result.text = os.into_string();
            result.tsv = os_tsv.into_string();
            return Ok(result);
        }
        check_branch_supported(in_type)?;

        match branch(in_type) {
            Branch::Features => report_features(path, options, &mut os, &mut os_tsv, &mut result)?,
            Branch::Peaks => {
                super::peaks::report(path, in_type, options, &mut os, &mut os_tsv, &mut result)?
            }
            // check_branch_supported refused these; kept exhaustive for the compiler.
            Branch::Unported | Branch::UnportedPeaks => return Err(unported_branch(in_type)),
            Branch::ImagingPeaks | Branch::NotLoadable => {
                return Err(refused_by_source_loader(path, name, in_type));
            }
        }

        os.text("\n\n");
        result.text = os.into_string();
        result.tsv = os_tsv.into_string();
        Ok(result)
    }

    /// Load `filename` with default options (no flags), the source
    /// `run(filename)`; see [`FileInfo::run`] for the errors.
    pub fn run_default(&self, filename: impl AsRef<Path>) -> Result<FileInfoResult> {
        self.run(filename, &Options::default())
    }

    /// Compute every content metric: `-m`, `-p` and `-s` on, `-v`, `-i`, `-d`
    /// and `-c` off, the source `runAll`; see [`FileInfo::run`] for the
    /// errors.
    pub fn run_all(&self, filename: impl AsRef<Path>) -> Result<FileInfoResult> {
        let options = Options {
            meta: true,
            processing: true,
            statistics: true,
            ..Options::default()
        };
        self.run(filename, &options)
    }

    /// The human-readable report a run cached in `result`, identical to the
    /// FileInfo tool's `-out` output, the source `toText(r)`.
    ///
    /// The source returns a copy; this borrows the cache.
    pub fn to_text(result: &FileInfoResult) -> &str {
        &result.text
    }

    /// The source `toText(r, options)`: `options` is accepted for API symmetry
    /// only, as in the source, and the rendering reflects the options passed to
    /// the run.
    pub fn to_text_with_options<'a>(result: &'a FileInfoResult, options: &Options) -> &'a str {
        let _ = options;
        &result.text
    }

    /// The TSV report a run cached in `result`, identical to the FileInfo
    /// tool's `-out_tsv` output, the source `toTSV(r)`.
    ///
    /// The source returns a copy; this borrows the cache.
    pub fn to_tsv(result: &FileInfoResult) -> &str {
        &result.tsv
    }

    /// The source `toTSV(r, options)`: `options` is accepted for API symmetry
    /// only, as in the source, and the rendering reflects the options passed to
    /// the run.
    pub fn to_tsv_with_options<'a>(result: &'a FileInfoResult, options: &Options) -> &'a str {
        let _ = options;
        &result.tsv
    }
}

#[cfg(feature = "featurexml")]
fn report_features(
    path: &Path,
    options: &Options,
    os: &mut ReportStream,
    os_tsv: &mut ReportStream,
    result: &mut FileInfoResult,
) -> Result<()> {
    super::features::report(path, options, os, os_tsv, result)
}

#[cfg(not(feature = "featurexml"))]
fn report_features(
    _path: &Path,
    _options: &Options,
    _os: &mut ReportStream,
    _os_tsv: &mut ReportStream,
    _result: &mut FileInfoResult,
) -> Result<()> {
    Err(Error::Unsupported(
        "FileInfo featureXML branch: this build lacks the featurexml feature".into(),
    ))
}

/// `FileHandler::getType` as the source `run` sees it.
///
/// The source content check opens a directory as a stream that yields no line
/// and returns `UNKNOWN` (`FileHandler.cpp:398-410`, `TextFile::load`);
/// [`FileHandler::get_type`] returns the I/O error of reading it. Only that
/// case is mapped: every other error, a missing file among them, is returned.
fn detect_type(path: &Path) -> Result<FileType> {
    match FileHandler::get_type(path) {
        Err(Error::Io(_)) if crate::system::file::is_directory(path) => Ok(FileType::Unknown),
        other => other,
    }
}

/// Which part of the source report a file type reaches.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Branch {
    /// The featureXML branch.
    Features,
    /// The peak-file branch with a native loader: DTA, DTA2D and mzML.
    Peaks,
    /// A non-peak branch of the source that is not ported.
    Unported,
    /// The peak-file branch for a type the source loads and this path does not.
    UnportedPeaks,
    /// imzML, which the source loader refuses with `InvalidFileType`.
    ImagingPeaks,
    /// Every other type: the source loader throws `ParseError`.
    NotLoadable,
}

fn branch(in_type: FileType) -> Branch {
    match in_type {
        FileType::FeatureXml => Branch::Features,
        FileType::Dta | FileType::Dta2d | FileType::MzMl => Branch::Peaks,
        FileType::ConsensusXml
        | FileType::IdXml
        | FileType::MzIdentMl
        | FileType::Fasta
        | FileType::PepXml
        | FileType::MzTab
        | FileType::TransformationXml
        | FileType::Pqp => Branch::Unported,
        // Thermo RAW and Bruker TDF load in the source built with its default
        // WITH_THERMO_RAW and WITH_OPENTIMS options; without them the source
        // loader throws ParseError. Neither has a native reader here.
        FileType::MzXml
        | FileType::MzData
        | FileType::Mgf
        | FileType::Ms2
        | FileType::SqMass
        | FileType::Xmass
        | FileType::Msp
        | FileType::Raw
        | FileType::BrukerTdf => Branch::UnportedPeaks,
        FileType::ImzMl => Branch::ImagingPeaks,
        _ => Branch::NotLoadable,
    }
}

/// Refuse the flags this port does not run, before anything is written. The
/// source's order is the report's: `-v`, then `-i`, then the content, so `-v` is
/// refused here and the branch only after the index check has had its turn.
fn check_flags_supported(options: &Options) -> Result<()> {
    if options.validate {
        return Err(Error::Unsupported(
            "FileInfo schema and semantic validation (-v) is not ported".into(),
        ));
    }
    Ok(())
}

/// Refuse, before the file is loaded, every branch this port does not run.
fn check_branch_supported(in_type: FileType) -> Result<()> {
    match branch(in_type) {
        Branch::Features | Branch::Peaks => Ok(()),
        Branch::Unported | Branch::UnportedPeaks => Err(unported_branch(in_type)),
        // Refused by the source loader too; reported once the report runs.
        Branch::ImagingPeaks | Branch::NotLoadable => Ok(()),
    }
}

fn unported_branch(in_type: FileType) -> Error {
    match branch(in_type) {
        Branch::UnportedPeaks => Error::Unsupported(format!(
            "FileInfo peak-file branch for {} input is not ported",
            in_type.name()
        )),
        _ => Error::Unsupported(format!("FileInfo {} branch is not ported", in_type.name())),
    }
}

/// The error `FileHandler::loadExperiment` gives for a type it cannot load
/// into an experiment (`FileHandler.cpp:855-1008`): it detects the type
/// itself, refuses one outside the allowed forced type with `ParseError`
/// (here [`Error::InvalidValue`], as [`FileHandler`] maps it), refuses imzML
/// with `InvalidFileType` (here [`Error::InvalidValue`]) and every other type
/// with `ParseError` (here [`Error::Parse`]).
fn refused_by_source_loader(path: &Path, name: &str, in_type: FileType) -> Error {
    let detected = match FileHandler::get_type(path) {
        Ok(detected) => detected,
        Err(error) => return error,
    };
    if detected != in_type {
        return Error::InvalidValue(format!(
            "{name}: type {} is not allowed for loading an experiment; allowed types are: {}",
            detected.name(),
            in_type.name()
        ));
    }
    if in_type == FileType::ImzMl {
        return Error::InvalidValue(format!(
            "{name}: imzML is a mass spectrometry imaging format; load it via ImzMLFile into an MSImagingExperiment"
        ));
    }
    Error::Parse {
        line: 0,
        message: format!("{name}: type is not supported for loading experiments"),
    }
}

/// One of the two C++ string streams a report is written into, with its
/// stream precision tracked explicitly.
#[derive(Debug)]
pub(crate) struct ReportStream {
    buffer: String,
    precision: u32,
}

impl ReportStream {
    /// An empty stream at [`DEFAULT_STREAM_PRECISION`].
    pub(crate) fn new() -> Self {
        Self {
            buffer: String::new(),
            precision: DEFAULT_STREAM_PRECISION,
        }
    }

    /// `os << text`.
    pub(crate) fn text(&mut self, text: &str) -> &mut Self {
        self.buffer.push_str(text);
        self
    }

    /// `os << value` for an integer or a string-like value.
    pub(crate) fn value(&mut self, value: impl Display) -> &mut Self {
        // Writing into a String cannot fail.
        let _ = write!(self.buffer, "{value}");
        self
    }

    /// `os << value` for a `double` (or a promoted `float`) in the default
    /// float field, at the current stream precision.
    pub(crate) fn double(&mut self, value: f64) -> &mut Self {
        self.buffer.push_str(&ostream_g(value, self.precision));
        self
    }

    /// `os.precision(precision)`, which stays in force for later output.
    pub(crate) fn set_precision(&mut self, precision: u32) {
        self.precision = precision;
    }

    /// The accumulated text, `os.str()`.
    pub(crate) fn into_string(self) -> String {
        self.buffer
    }
}

/// `StringUtils::number(value, 2)`.
fn number2(value: f64) -> String {
    fixed_truncated(value, 2)
}

/// The ranges of `manager` as a [`RangeSet`], the source `extractRangeSet_`.
///
/// # Errors
///
/// None in practice: every dimension read is checked for presence first.
pub(crate) fn range_set(manager: &RangeManager) -> Result<RangeSet> {
    let dimension = |dim: MSDim| -> Result<Option<Range>> {
        if !manager.has_dim(dim) || manager.is_dim_empty(dim)? {
            return Ok(None);
        }
        Ok(Some(Range {
            min: manager.min(dim)?,
            max: manager.max(dim)?,
        }))
    };
    Ok(RangeSet {
        rt: dimension(MSDim::Rt)?,
        mz: dimension(MSDim::Mz)?,
        mobility: dimension(MSDim::Mobility)?,
        intensity: dimension(MSDim::Intensity)?,
        has_mobility: manager.has_dim(MSDim::Mobility),
    })
}

/// The human-readable lines of one range block below its title
/// (`writeRangesHumanReadable_`): retention time, m/z, the ion-mobility line
/// when `mobility_line` is set, and intensity followed by an empty line.
pub(crate) fn write_ranges_text(os: &mut ReportStream, ranges: &RangeSet, mobility_line: bool) {
    match ranges.rt {
        None => {
            os.text("  retention time: <none> .. <none> sec (<none> min)\n");
        }
        Some(rt) => {
            os.text("  retention time: ")
                .text(&number2(rt.min))
                .text(" .. ")
                .text(&number2(rt.max))
                .text(" sec (")
                .text(&fixed_truncated((rt.max - rt.min) / 60.0, 1))
                .text(" min)\n");
        }
    }
    write_range_line(os, "  mass-to-charge: ", ranges.mz, "\n");
    if mobility_line {
        write_range_line(os, "  ion mobility: ", ranges.mobility, "\n");
    }
    write_range_line(os, "  intensity: ", ranges.intensity, "\n\n");
}

fn write_range_line(os: &mut ReportStream, label: &str, range: Option<Range>, end: &str) {
    os.text(label);
    match range {
        None => os.text("<none> .. <none>"),
        Some(range) => os
            .text(&number2(range.min))
            .text(" .. ")
            .text(&number2(range.max)),
    };
    os.text(end);
}

/// The TSV lines of one range block (`writeRangesMachineReadable_`), each key
/// starting with `prefix`: retention time, m/z, ion mobility only when a
/// mobility range is present, and intensity; absent bounds print `<none>`.
pub(crate) fn write_ranges_tsv(os_tsv: &mut ReportStream, prefix: &str, ranges: &RangeSet) {
    write_range_tsv(os_tsv, prefix, "retention time", ranges.rt);
    write_range_tsv(os_tsv, prefix, "mass-to-charge", ranges.mz);
    if ranges.mobility.is_some() {
        write_range_tsv(os_tsv, prefix, "ion-mobility", ranges.mobility);
    }
    write_range_tsv(os_tsv, prefix, "intensity", ranges.intensity);
}

fn write_range_tsv(os_tsv: &mut ReportStream, prefix: &str, label: &str, range: Option<Range>) {
    let (min, max) = match range {
        Some(range) => (number2(range.min), number2(range.max)),
        None => ("<none>".to_owned(), "<none>".to_owned()),
    };
    for (bound, value) in [("min", min), ("max", max)] {
        os_tsv
            .text(prefix)
            .text(label)
            .text(": ")
            .text(bound)
            .text("\t")
            .text(&value)
            .text("\n");
    }
}

/// `printChargeDistribution`: `"<header> distribution:"`, one line per charge
/// in ascending order with its TSV twin, and an empty line.
pub(crate) fn write_charge_distribution(
    os: &mut ReportStream,
    os_tsv: &mut ReportStream,
    header: &str,
    charges: &BTreeMap<i32, u64>,
) {
    os.text(header).text(" distribution:\n");
    for (charge, count) in charges {
        os.text("  charge ")
            .value(charge)
            .text(": ")
            .value(count)
            .text("x\n");
        os_tsv
            .text("general: charge distribution: charge: ")
            .value(charge)
            .text("\t")
            .value(count)
            .text("\n");
    }
    os.text("\n");
}

/// The `-m` section title.
pub(crate) fn write_meta_title(os: &mut ReportStream) {
    os.text("\n-- Meta information --\n\n");
}

/// The `-p` section title.
pub(crate) fn write_processing_title(os: &mut ReportStream) {
    os.text("\n-- Data processing information --\n\n");
}

/// The `-s` section title.
pub(crate) fn write_statistics_title(os: &mut ReportStream) {
    os.text("\n-- Statistics --\n\n");
}

/// The body of the `-p` section after the branch-specific preamble, and its
/// structured steps (`FileInfo.cpp:2129-2188`).
pub(crate) fn write_processing(
    os: &mut ReportStream,
    os_tsv: &mut ReportStream,
    processing: &[&DataProcessing],
    result: &mut FileInfoResult,
) {
    if processing.is_empty() {
        os.text("No information about data processing available!\n\n");
    }
    for (index, step) in processing.iter().enumerate() {
        let number = index + 1;
        let step = ProcessingStep {
            software_name: step.software.name.clone(),
            software_version: step.software.version.clone(),
            completion_time: step
                .completion_time
                .map_or_else(|| "0000-00-00 00:00:00".to_owned(), |time| time.get()),
            actions: step
                .actions
                .iter()
                .map(|action| action.name().to_owned())
                .collect(),
        };
        os.text("Processing ").value(number).text(":\n");
        os.text("  software name:    ")
            .text(&step.software_name)
            .text("\n");
        os.text("  software version: ")
            .text(&step.software_version)
            .text("\n");
        os.text("  completion time:  ")
            .text(&step.completion_time)
            .text("\n");
        for (label, value) in [
            ("software name", &step.software_name),
            ("software version", &step.software_version),
            ("completion time", &step.completion_time),
        ] {
            os_tsv
                .text("data processing: ")
                .value(number)
                .text(": ")
                .text(label)
                .text("\t")
                .text(value)
                .text("\n");
        }
        let actions = step.actions.join(", ");
        os.text("  actions:          ").text(&actions).text("\n\n");
        os_tsv
            .text("data processing: ")
            .value(number)
            .text(": actions\t")
            .text(&actions)
            .text("\n");
        result.processing.push(step);
    }
}

/// `SummaryStatistics(values)` over a collected sample, refusing a sample
/// above [`FileInfo::MAX_STATISTICS_VALUES`].
///
/// # Errors
///
/// [`Error::InvalidValue`] above the ceiling, and as
/// [`SummaryStatistics::new`] for a NaN.
pub(crate) fn summarize(values: &mut [f64]) -> Result<SummaryStatistics> {
    check_statistics_values(values.len())?;
    SummaryStatistics::new(values)
}

/// Refuse a statistics sample of `count` values above
/// [`FileInfo::MAX_STATISTICS_VALUES`], before it is collected.
///
/// # Errors
///
/// [`Error::InvalidValue`] above the ceiling.
pub(crate) fn check_statistics_values(count: usize) -> Result<()> {
    if count > FileInfo::MAX_STATISTICS_VALUES {
        return Err(Error::InvalidValue(format!(
            "FileInfo statistics block of {count} values exceeds {}",
            FileInfo::MAX_STATISTICS_VALUES
        )));
    }
    Ok(())
}

/// An empty vector with room for `count` values, allocated fallibly after the
/// ceiling check.
///
/// # Errors
///
/// [`Error::InvalidValue`] above [`FileInfo::MAX_STATISTICS_VALUES`] or when
/// the allocation fails.
pub(crate) fn statistics_buffer(count: usize) -> Result<Vec<f64>> {
    check_statistics_values(count)?;
    let mut values = Vec::new();
    values.try_reserve_exact(count).map_err(|_| {
        Error::InvalidValue(format!(
            "cannot allocate a FileInfo statistics block of {count} values"
        ))
    })?;
    Ok(values)
}

/// The `operator<<` for `SummaryStatistics`: eight indented lines with the
/// doubles at the stream's current precision.
pub(crate) fn write_summary_text(os: &mut ReportStream, stats: &SummaryStatistics) {
    os.text("  num. of values: ").value(stats.count).text("\n");
    for (label, value) in [
        ("  mean:           ", stats.mean),
        ("  minimum:        ", stats.min),
        ("  lower quartile: ", stats.lowerq),
        ("  median:         ", stats.median),
        ("  upper quartile: ", stats.upperq),
        ("  maximum:        ", stats.max),
        ("  variance:       ", stats.variance),
    ] {
        os.text(label).double(value).text("\n");
    }
}
