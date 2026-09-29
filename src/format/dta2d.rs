// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Three-column DTA2D spectra and explicit MS1 TIC projection.
//!
//! Port of `FORMAT/DTA2DFile.h`; see `docs/TEXT_PEAK_LIST_SUPPORT.md`. The
//! source class derives from `ProgressLogger`:
//! [`load_with_progress`](crate::format::dta2d::load_with_progress),
//! [`store_with_progress`](crate::format::dta2d::store_with_progress) and
//! [`store_tic_with_progress`](crate::format::dta2d::store_tic_with_progress)
//! make the progress
//! calls of its `load`, `store` and `storeTIC` on a caller's logger; every
//! other entry point runs the same code and reports nothing.
//!
//! The writers spell every number as the source does, through
//! `precisionWrapper` and so `StringUtils::toStr(value, true)`: retention time
//! and m/z with [`to_str`], intensity with [`to_str_f32`]. By default they
//! refuse a value whose text would not read back exactly, and anything else
//! DTA2D cannot carry; [`WriteOptions::source`](crate::format::dta2d::WriteOptions::source) writes what the source writes.
//!
//! `load`, `store` and `storeTIC` are templates over `MapType`, documented as
//! "a MSExperiment or have the same interface"; every instantiation in the
//! pinned source is `PeakMap` (`FileHandler.cpp:877-883`, `:1159-1165`, the
//! class test), and this crate has one experiment type, so `MapType` is
//! [`MSExperiment`](crate::MSExperiment) here.
//!
//! [`to_str`]: crate::format::file_info::text_format::to_str
//! [`to_str_f32`]: crate::format::file_info::text_format::to_str_f32

pub use super::ms2::Limits;
use super::ms2::{
    Counter, TextInput, check_experiment, check_peaks, check_spectrum_payload, fields, increment,
    invalid, loaded_document, parse_intensity, parse_number, push_peak, push_spectrum,
    record_document, trim, unsupported,
};
use super::{PeakFileOptions, parse_error};
use crate::concept::progress_logger::{ProgressLogger, ProgressReporter};
use crate::format::file_info::text_format::{to_str, to_str_f32};
use crate::{MSExperiment, MSSpectrum, Peak1D, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::ops::Range;
use std::path::Path;

/// Settings of the writers: output ceilings and what may be discarded.
///
/// The default refuses anything the file would not carry back, so a stored
/// experiment reads back unchanged. [`WriteOptions::source`] writes what
/// source `DTA2DFile::store` and `storeTIC` write and discards the rest
/// silently, as the source does.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct WriteOptions {
    /// Output ceilings; the same [`Limits`] as the reader's.
    pub limits: Limits,
    /// Write what the source writes and drop what DTA2D cannot represent,
    /// instead of refusing it.
    ///
    /// With `false`, the default, `store` refuses spectra that are not
    /// nonempty MS1 spectra without precursors, names, annotations,
    /// identifications, data arrays or foreign native IDs; experiment metadata
    /// and chromatograms; adjacent spectra the reader would merge (retention
    /// times within 0.0001 s) or a first spectrum the reader's -1 sentinel
    /// would change; and, in `store` and `storeTIC` alike, any number whose
    /// source text does not read back as the same value. That text has 15
    /// fraction digits for a double and 6 for a float in `[0.01, 10000)`, so
    /// for example an m/z below about 8, or an intensity `0.123456789`, would
    /// change.
    ///
    /// With `true` every spectrum is written as the source writes it: its
    /// peaks under its retention time, whatever its MS level or metadata, an
    /// empty spectrum contributing no line. Nonfinite values and the ceilings
    /// are refused in both modes.
    pub discard_unrepresentable: bool,
}
impl WriteOptions {
    /// The source `DTA2DFile::store`/`storeTIC` behaviour: write what the
    /// source writes and discard the rest, under default ceilings.
    pub fn source() -> Self {
        Self {
            limits: Limits::default(),
            discard_unrepresentable: true,
        }
    }
}
impl From<Limits> for WriteOptions {
    /// Lossless writing under `limits`.
    fn from(limits: Limits) -> Self {
        Self {
            limits,
            discard_unrepresentable: false,
        }
    }
}

/// The three source-consumed PeakFileOptions filters. Ranges are half-open,
/// including the minimum and excluding the maximum, as in OpenMS DRange.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadOptions {
    /// Keep spectra whose retention time, in seconds, lies in this range.
    pub rt_range: Option<Range<f64>>,
    /// Keep peaks whose m/z lies in this range.
    pub mz_range: Option<Range<f64>>,
    /// Keep peaks whose intensity lies in this range.
    pub intensity_range: Option<Range<f64>>,
    /// Native input ceilings; the source reader has none.
    pub limits: Limits,
}
impl ReadOptions {
    /// The reader settings a source `PeakFileOptions` gives `DTA2DFile::load`,
    /// the counterpart of `DTA2DFile::getOptions() = options`.
    ///
    /// The source's `load` reads only `hasRTRange`/`getRTRange`,
    /// `hasMZRange`/`getMZRange` and `hasIntensityRange`/`getIntensityRange`
    /// (`DTA2DFile.h:211-240`); each range that is set becomes the matching
    /// half-open range here. Every other field of `options` has no DTA2D
    /// behaviour in the source and is ignored, as there. `limits` keep their
    /// defaults.
    pub fn from_peak_options(options: &PeakFileOptions) -> Self {
        let range = |present: bool, range: crate::kernel::NumericRange| {
            present.then_some(range.min..range.max)
        };
        Self {
            rt_range: range(options.has_rt_range(), options.rt_range()),
            mz_range: range(options.has_mz_range(), options.mz_range()),
            intensity_range: range(options.has_intensity_range(), options.intensity_range()),
            limits: Limits::default(),
        }
    }
    fn validate(&self) -> Result<()> {
        self.limits.validate()?;
        for range in [&self.rt_range, &self.mz_range, &self.intensity_range]
            .into_iter()
            .flatten()
        {
            if !range.start.is_finite() || !range.end.is_finite() || range.start > range.end {
                return Err(invalid("DTA2D ranges require finite ordered endpoints"));
            }
        }
        Ok(())
    }
}
fn includes(range: &Option<Range<f64>>, value: f64) -> bool {
    range.as_ref().is_none_or(|r| r.contains(&value))
}
/// Read a DTA2D stream with default options.
///
/// # Errors
///
/// As [`read_with_options`].
pub fn read(reader: impl BufRead) -> Result<MSExperiment> {
    read_with_options(reader, &ReadOptions::default())
}
/// Read a DTA2D stream, applying the three source filters of `options`.
///
/// # Errors
///
/// [`Parse`](crate::Error::Parse) for a malformed header or data line, and
/// [`InvalidValue`](crate::Error::InvalidValue) for invalid options or an
/// exceeded limit.
pub fn read_with_options(reader: impl BufRead, options: &ReadOptions) -> Result<MSExperiment> {
    read_reporting(reader, options, &mut ProgressReporter::silent())
}
/// The reader, with the source's `setProgress(0)` each time a new spectrum
/// begins, after the previous one was added (`DTA2DFile.h:209-217`).
fn read_reporting(
    reader: impl BufRead,
    options: &ReadOptions,
    progress: &mut ProgressReporter<'_>,
) -> Result<MSExperiment> {
    options.validate()?;
    let mut input = TextInput::new(reader, &options.limits)?;
    let mut result = MSExperiment::default();
    let mut current = MSSpectrum::default(); // source initial RT=-1 sentinel
    let (mut scans, mut peaks, mut groups) = (0, 0, 0);
    let (mut rt_column, mut mz_column, mut intensity_column) = (0, 1, 2);
    let mut minutes = false;
    let mut line = String::new();
    while input.next_line(&mut line)? {
        let text = trim(&line);
        if text.is_empty() {
            continue;
        }
        let delimiter = if text.contains('\t') { '\t' } else { ' ' };
        if let Some(header) = text.strip_prefix('#') {
            // Source reads exactly the first three header fields; extra fields
            // are ignored. Short headers are a checked error, not out-of-bounds.
            let mut tokens = trim(header).split(delimiter);
            let mut columns = [None; 3];
            for i in 0..3 {
                let token = tokens
                    .next()
                    .ok_or_else(|| parse_error(input.line, "short DTA2D header"))?;
                let kind = if ["SEC", "RT", "RETENTION_TIME"]
                    .iter()
                    .any(|s| token.eq_ignore_ascii_case(s))
                {
                    0
                } else if token.eq_ignore_ascii_case("MIN") {
                    minutes = true;
                    0
                } else if ["MZ", "MASS-TO-CHARGE"]
                    .iter()
                    .any(|s| token.eq_ignore_ascii_case(s))
                {
                    1
                } else if ["INT", "IT", "INTENSITY"]
                    .iter()
                    .any(|s| token.eq_ignore_ascii_case(s))
                {
                    2
                } else {
                    return Err(parse_error(input.line, "invalid DTA2D header keyword"));
                };
                if columns[kind].replace(i).is_some() {
                    return Err(parse_error(input.line, "duplicate DTA2D header dimension"));
                }
            }
            [rt_column, mz_column, intensity_column] =
                columns.map(|v| v.expect("three distinct dimensions"));
            continue;
        }
        increment(
            &mut peaks,
            options.limits.max_peaks,
            "DTA2D peak limit exceeded",
        )?;
        let values = fields::<3>(text.split(delimiter), input.line)?;
        let mz = parse_number(values[mz_column], input.line, "m/z")?;
        let intensity = parse_intensity(values[intensity_column], input.line)?;
        let rt = parse_number(values[rt_column], input.line, "retention time")?
            * if minutes { 60. } else { 1. };
        if !rt.is_finite() {
            return Err(parse_error(
                input.line,
                "DTA2D retention time conversion overflow",
            ));
        }
        let changed = (rt - current.rt).abs() > 0.0001;
        if peaks == 1 || changed {
            increment(
                &mut groups,
                options.limits.max_spectra,
                "DTA2D scan limit exceeded",
            )?;
        }
        if changed {
            scans += 1;
            let previous = std::mem::replace(
                &mut current,
                MSSpectrum {
                    rt,
                    native_id: format!("index={}", scans - 1),
                    ..Default::default()
                },
            );
            if !previous.peaks.is_empty() && includes(&options.rt_range, previous.rt) {
                push_spectrum(&mut result, previous)?;
            }
            progress.set(0)?;
        }
        if includes(&options.mz_range, mz)
            && includes(&options.intensity_range, f64::from(intensity))
        {
            push_peak(&mut current, Peak1D::new(mz, intensity))?;
        }
    }
    if !current.peaks.is_empty() && includes(&options.rt_range, current.rt) {
        push_spectrum(&mut result, current)?;
    }
    Ok(result)
}
/// Replace `destination` only after a successful read, with default options.
///
/// # Errors
///
/// As [`read_with_options`]; `destination` is then unchanged.
pub fn read_into(reader: impl BufRead, destination: &mut MSExperiment) -> Result<()> {
    *destination = read(reader)?;
    Ok(())
}
/// Replace `destination` only after a successful read.
///
/// # Errors
///
/// As [`read_with_options`]; `destination` is then unchanged.
pub fn read_into_with_options(
    reader: impl BufRead,
    destination: &mut MSExperiment,
    options: &ReadOptions,
) -> Result<()> {
    *destination = read_with_options(reader, options)?;
    Ok(())
}
/// Read a DTA2D file with default options.
///
/// # Errors
///
/// As [`load_with_options`].
pub fn load(path: impl AsRef<Path>) -> Result<MSExperiment> {
    load_with_options(path, &ReadOptions::default())
}
/// Read a DTA2D file, source `DTA2DFile::load` (`DTA2DFile.h:72-248`).
///
/// As in the source, the result records the document's identity:
/// `settings.document.loaded_file_path` is the path as given (a relative path
/// prefixed with the current directory) and `loaded_file_type` is the type
/// recognised from the file's content, the source's
/// `setLoadedFileType(filename)` and `setLoadedFilePath(filename)`
/// (`DTA2DFile.h:97-98`), so `DTA2DFile_test.cpp:72-73` holds.
///
/// # Errors
///
/// [`Io`](crate::Error::Io) when the file cannot be opened, where the source
/// throws `Exception::FileNotFound`, `FileNotReadable` or `IOException`;
/// [`InvalidValue`](crate::Error::InvalidValue) for a path that is not UTF-8;
/// otherwise as [`read_with_options`].
pub fn load_with_options(path: impl AsRef<Path>, options: &ReadOptions) -> Result<MSExperiment> {
    load_reporting(path, options, &mut ProgressReporter::silent())
}
/// Read a DTA2D file, reporting progress to `logger` as source
/// `DTA2DFile::load` does (`DTA2DFile.h:72-248`).
///
/// The calls are the source's: `startProgress(0, 0, "loading DTA2D file")`
/// before the file is opened, `setProgress(0)` each time a new spectrum
/// begins (with an equal begin and end, the command backend prints one dot
/// per call), and `endProgress()` once the whole file is read. The result is
/// the one [`load_with_options`] returns, and so is every error: both run the
/// same code, whose calls go nowhere for [`load_with_options`].
///
/// A failure after the start leaves the section open, as in the source, where
/// the exception bypasses `endProgress`: no `-- done` line is printed, the
/// nesting depth stays one level deeper, and a command backend of `logger`
/// refuses its next start (`StopWatch is already started!`). This includes a
/// file that cannot be opened (the source's `FileNotFound`, thrown after the
/// start) and invalid `options`, a native check the source does not have.
///
/// # Errors
///
/// As [`load_with_options`], plus the errors of the progress calls
/// ([`ProgressLogger::start_progress`] and its siblings).
pub fn load_with_progress(
    path: impl AsRef<Path>,
    options: &ReadOptions,
    logger: &mut ProgressLogger,
) -> Result<MSExperiment> {
    load_reporting(path, options, &mut ProgressReporter::new(Some(logger)))
}
/// The source's `load`: its section around the reader.
fn load_reporting(
    path: impl AsRef<Path>,
    options: &ReadOptions,
    progress: &mut ProgressReporter<'_>,
) -> Result<MSExperiment> {
    let path = path.as_ref();
    progress.start(0, 0, "loading DTA2D file")?;
    let input = BufReader::new(File::open(path)?);
    let document = loaded_document(path, "DTA2D")?;
    let mut experiment = read_reporting(input, options, progress)?;
    record_document(&mut experiment, document);
    progress.end()?;
    Ok(experiment)
}

/// Write `experiment` as DTA2D with default ceilings.
///
/// # Errors
///
/// As [`write_with_options`].
pub fn write(writer: impl Write, experiment: &MSExperiment) -> Result<()> {
    write_with_options(writer, experiment, &WriteOptions::default())
}
/// Write `experiment` as DTA2D; everything is checked before the first byte.
///
/// The body of source `DTA2DFile::store` (`DTA2DFile.h:258-287`): the header
/// `#SEC\tMZ\tINT`, then one `rt\tmz\tintensity` line per peak, spectrum by
/// spectrum, each number in the source's text (see the module documentation).
///
/// # Errors
///
/// [`Unsupported`](crate::Error::Unsupported) for data DTA2D cannot hold
/// (see [`WriteOptions::discard_unrepresentable`] and
/// `docs/TEXT_PEAK_LIST_SUPPORT.md`), [`InvalidValue`](crate::Error::InvalidValue)
/// for a nonfinite value or an exceeded ceiling, and [`Io`](crate::Error::Io)
/// for a write failure.
pub fn write_with_options(
    mut writer: impl Write,
    experiment: &MSExperiment,
    options: &WriteOptions,
) -> Result<()> {
    preflight(experiment, options)?;
    render(
        &mut writer,
        experiment,
        &mut ProgressReporter::silent(),
        false,
    )?;
    writer.flush()?;
    Ok(())
}
/// Store `experiment` as a DTA2D file with default ceilings.
///
/// # Errors
///
/// As [`store_with_options`].
pub fn store(path: impl AsRef<Path>, experiment: &MSExperiment) -> Result<()> {
    store_with_options(path, experiment, &WriteOptions::default())
}
/// Store `experiment` as a DTA2D file; the path is created only after the
/// whole experiment passed the checks.
///
/// # Errors
///
/// As [`write_with_options`], plus [`Io`](crate::Error::Io) when the file
/// cannot be created.
pub fn store_with_options(
    path: impl AsRef<Path>,
    experiment: &MSExperiment,
    options: &WriteOptions,
) -> Result<()> {
    store_reporting(path, experiment, options, &mut ProgressReporter::silent())
}
/// Store `experiment`, reporting progress to `logger` as source
/// `DTA2DFile::store` does (`DTA2DFile.h:258-287`).
///
/// The calls are the source's: `startProgress(0, spectra, "storing DTA2D
/// file")` before the file is created, `setProgress(i)` before spectrum `i`
/// is written, and `endProgress()` after the file is flushed. The written
/// bytes and every error are those of [`store_with_options`], which runs the
/// same code with the calls going nowhere. Its checks run before the start,
/// so an experiment refused there makes no call; the source has none of them.
/// A failure after the start, such as a file that cannot be created (the
/// source's `UnableToCreateFile`), leaves the section open, as described at
/// [`load_with_progress`].
///
/// # Errors
///
/// As [`store_with_options`], plus the errors of the progress calls.
pub fn store_with_progress(
    path: impl AsRef<Path>,
    experiment: &MSExperiment,
    options: &WriteOptions,
    logger: &mut ProgressLogger,
) -> Result<()> {
    store_reporting(
        path,
        experiment,
        options,
        &mut ProgressReporter::new(Some(logger)),
    )
}
/// The source's `store`: its section around the writer.
fn store_reporting(
    path: impl AsRef<Path>,
    experiment: &MSExperiment,
    options: &WriteOptions,
    progress: &mut ProgressReporter<'_>,
) -> Result<()> {
    preflight(experiment, options)?;
    progress.start_count(experiment.spectra.len(), "storing DTA2D file")?;
    let mut writer = BufWriter::new(File::create(path)?);
    render(&mut writer, experiment, progress, false)?;
    writer.flush()?;
    progress.end()
}
fn preflight(experiment: &MSExperiment, options: &WriteOptions) -> Result<()> {
    let limits = &options.limits;
    if options.discard_unrepresentable {
        limits.validate()?;
        if experiment.spectra.len() > limits.max_spectra {
            return Err(invalid("spectrum limit exceeded"));
        }
        let mut peaks = 0;
        for s in &experiment.spectra {
            check_peaks(s, &mut peaks, limits)?;
        }
        return render(
            &mut Counter(limits.max_output_bytes),
            experiment,
            &mut ProgressReporter::silent(),
            false,
        );
    }
    check_experiment(experiment, limits)?;
    let mut peaks = 0;
    let mut previous: Option<f64> = None;
    for s in &experiment.spectra {
        check_spectrum_payload(s)?;
        if s.ms_level != 1 || !s.precursors.is_empty() || s.peaks.is_empty() {
            return Err(unsupported(
                "DTA2D requires nonempty MS1 spectra without precursors",
            ));
        }
        if previous.is_some_and(|rt| (s.rt - rt).abs() <= 0.0001) {
            return Err(unsupported(
                "DTA2D would merge adjacent spectra within its RT tolerance",
            ));
        }
        if previous.is_none() && s.rt != -1.0 && (s.rt + 1.0).abs() <= 0.0001 {
            return Err(unsupported(
                "DTA2D initial RT sentinel would change the first retention time",
            ));
        }
        previous = Some(s.rt);
        check_peaks(s, &mut peaks, limits)?;
    }
    render(
        &mut Counter(limits.max_output_bytes),
        experiment,
        &mut ProgressReporter::silent(),
        true,
    )
}
/// Source `StringUtils::toStr(double, true)`, refused with `exact` when the
/// text does not parse back to the same bits.
fn double_text(value: f64, exact: bool) -> Result<String> {
    let text = to_str(value);
    if exact && text.parse::<f64>().map(f64::to_bits) != Ok(value.to_bits()) {
        return Err(unsupported(
            "DTA2D source number text would not read back exactly",
        ));
    }
    Ok(text)
}
/// Source `StringUtils::toStr(float, true)`, refused as [`double_text`].
fn float_text(value: f32, exact: bool) -> Result<String> {
    let text = to_str_f32(value);
    if exact && text.parse::<f32>().map(f32::to_bits) != Ok(value.to_bits()) {
        return Err(unsupported(
            "DTA2D source number text would not read back exactly",
        ));
    }
    Ok(text)
}
/// The file body, with the source's `setProgress(count++)` before each
/// spectrum (`DTA2DFile.h:275-277`). With `exact`, a number whose text does
/// not read back is refused.
fn render(
    writer: &mut impl Write,
    experiment: &MSExperiment,
    progress: &mut ProgressReporter<'_>,
    exact: bool,
) -> Result<()> {
    writeln!(writer, "#SEC\tMZ\tINT")?;
    for (index, s) in experiment.spectra.iter().enumerate() {
        progress.set_count(index)?;
        if s.peaks.is_empty() {
            continue;
        }
        let rt = double_text(s.rt, exact)?;
        for p in &s.peaks {
            let mz = double_text(p.mz, exact)?;
            let intensity = float_text(p.intensity, exact)?;
            writeln!(writer, "{rt}\t{mz}\t{intensity}")?;
        }
    }
    Ok(())
}
/// Write the MS1 total-ion-current projection of `experiment` as DTA2D with
/// default settings.
///
/// Explicit lossy projection, source `DTA2DFile::storeTIC`
/// (`DTA2DFile.h:297-320`) through `calculateTIC()` with its defaults: MS1
/// TIC only, one row per scan (including empty scans), m/z written as the
/// literal `0`, and f32 accumulation in storage order. Other data are
/// unconsumed; the numbers are spelled as [`write`](fn@write) spells them.
///
/// # Errors
///
/// As [`write_tic_with_options`].
pub fn write_tic(writer: impl Write, experiment: &MSExperiment) -> Result<()> {
    write_tic_with_options(writer, experiment, &WriteOptions::default())
}
/// [`write_tic`] under explicit settings.
///
/// Only the ceilings and the number-text check of `options` apply; the
/// projection itself always discards everything but MS1 retention times and
/// summed intensities, as the source's does.
///
/// # Errors
///
/// [`InvalidValue`](crate::Error::InvalidValue) for an exceeded ceiling or a
/// nonfinite retention time, intensity or TIC sum;
/// [`Unsupported`](crate::Error::Unsupported), unless
/// [`WriteOptions::discard_unrepresentable`] is set, for a retention time or
/// sum whose source text does not read back exactly; and
/// [`Io`](crate::Error::Io) for a write failure.
pub fn write_tic_with_options(
    mut writer: impl Write,
    experiment: &MSExperiment,
    options: &WriteOptions,
) -> Result<()> {
    let points = tic_points(experiment, options)?;
    render_tic(&mut writer, &points)?;
    writer.flush()?;
    Ok(())
}
/// Store the MS1 TIC projection of [`write_tic`] in a file.
///
/// # Errors
///
/// As [`store_tic_with_options`].
pub fn store_tic(path: impl AsRef<Path>, experiment: &MSExperiment) -> Result<()> {
    store_tic_with_options(path, experiment, &WriteOptions::default())
}
/// Store the MS1 TIC projection under explicit ceilings; the path is created
/// only after the projection passed the checks.
///
/// # Errors
///
/// As [`write_tic_with_options`], plus [`Io`](crate::Error::Io) when the
/// file cannot be created.
pub fn store_tic_with_options(
    path: impl AsRef<Path>,
    experiment: &MSExperiment,
    options: &WriteOptions,
) -> Result<()> {
    store_tic_reporting(path, experiment, options, &mut ProgressReporter::silent())
}
/// Store the TIC projection, reporting progress to `logger` as source
/// `DTA2DFile::storeTIC` does (`DTA2DFile.h:297-320`).
///
/// The source starts a section over all spectra, `startProgress(0, spectra,
/// "storing DTA2D file")`, before the file is created, makes no call inside
/// it, and ends it after the file is closed. The written bytes and every error
/// are those of [`store_tic_with_options`], which runs the same code with the
/// calls going nowhere; its checks run before the start. A failure after the
/// start leaves the section open, as described at [`load_with_progress`].
///
/// # Errors
///
/// As [`store_tic_with_options`], plus the errors of the progress calls.
pub fn store_tic_with_progress(
    path: impl AsRef<Path>,
    experiment: &MSExperiment,
    options: &WriteOptions,
    logger: &mut ProgressLogger,
) -> Result<()> {
    store_tic_reporting(
        path,
        experiment,
        options,
        &mut ProgressReporter::new(Some(logger)),
    )
}
/// The source's `storeTIC`: its section around the writer.
fn store_tic_reporting(
    path: impl AsRef<Path>,
    experiment: &MSExperiment,
    options: &WriteOptions,
    progress: &mut ProgressReporter<'_>,
) -> Result<()> {
    let points = tic_points(experiment, options)?;
    progress.start_count(experiment.spectra.len(), "storing DTA2D file")?;
    let mut writer = BufWriter::new(File::create(path)?);
    render_tic(&mut writer, &points)?;
    writer.flush()?;
    progress.end()
}
/// The TIC rows as text, checked and counted against the output ceiling.
fn tic_points(experiment: &MSExperiment, options: &WriteOptions) -> Result<Vec<(String, String)>> {
    let exact = !options.discard_unrepresentable;
    let options = &options.limits;
    options.validate()?;
    if experiment.spectra.len() > options.max_spectra {
        return Err(invalid("DTA2D TIC spectrum limit exceeded"));
    }
    let mut total = 0;
    let mut points = Vec::new();
    points
        .try_reserve_exact(experiment.spectra.len())
        .map_err(|_| invalid("TIC allocation failed"))?;
    for s in &experiment.spectra {
        if s.ms_level != 1 {
            continue;
        }
        if s.peaks.len() > options.max_peaks.saturating_sub(total) {
            return Err(invalid("DTA2D TIC peak limit exceeded"));
        }
        total += s.peaks.len();
        if !s.rt.is_finite() || s.peaks.iter().any(|p| !p.intensity.is_finite()) {
            return Err(invalid("nonfinite TIC retention time or intensity"));
        }
        let sum = s.calculate_tic();
        if !sum.is_finite() {
            return Err(invalid("DTA2D TIC intensity sum exceeds f32 range"));
        }
        points.push((double_text(s.rt, exact)?, float_text(sum, exact)?));
    }
    render_tic(&mut Counter(options.max_output_bytes), &points)?;
    Ok(points)
}
fn render_tic(writer: &mut impl Write, points: &[(String, String)]) -> Result<()> {
    writeln!(writer, "#SEC\tMZ\tINT")?;
    for (rt, intensity) in points {
        writeln!(writer, "{rt}\t0\t{intensity}")?;
    }
    Ok(())
}
