// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Chromeleon HPLC text exports, ported from `FORMAT/ChromeleonFile.h`.
//!
//! An example of the expected format:
//!
//! ```text
//! Raw Data:
//! Time (min)    Step (s)    Value (mAU)
//! 0.003333    0.200    -0.002496
//! 0.006667    0.200    -0.017589
//! ...
//! ```
//!
//! A header of tab-separated `key<TAB>value` lines precedes the data; the keys
//! the source reads are listed at [`ChromeleonFile::load`](crate::format::chromeleon::ChromeleonFile::load). The data follow a
//! line that is exactly `Raw Data:` or `Chromatogram Data:` and the column
//! caption after it. See `docs/CHROMELEON_FILE_SUPPORT.md`.

use crate::concept::boost_regex::{BoostRegex, RegexOptions};
use crate::format::line_reader::{LineLimits, LineReader};
use crate::format::string_utils::to_double;
use crate::kernel::{ChromatogramPeak, MSChromatogram, MSExperiment};
use crate::metadata::MetaValue;
use crate::{Error, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Resource ceilings for one Chromeleon file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadLimits {
    /// Line length and total input size.
    pub lines: LineLimits,
    /// Largest number of chromatogram points.
    pub max_points: usize,
}

impl Default for ReadLimits {
    /// [`LineLimits::default`] and 100 million points (1.6 GB of peaks).
    fn default() -> Self {
        Self {
            lines: LineLimits::default(),
            max_points: 100_000_000,
        }
    }
}

/// Loads a Chromeleon HPLC text file into an [`MSExperiment`] (source
/// `ChromeleonFile`).
///
/// The source class is stateless; so is this one.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ChromeleonFile;

/// Where a header line's captured value goes.
#[derive(Clone, Copy)]
enum Target {
    /// An experiment meta value under this key.
    Meta(&'static str),
    /// The instrument's software name.
    SoftwareName,
    /// The instrument name.
    InstrumentName,
}

/// The source's header expressions, in the order it tries them. Every one is
/// `^<key>\t(.*)` compiled with `boost::regex::no_mod_s` and applied with
/// `regex_match` to the whole line.
const HEADER: [(&str, Target); 10] = [
    ("^Injection\t(.*)", Target::Meta("mzml_id")),
    ("^Channel\t(.*)", Target::Meta("acq_method_name")),
    ("^Processing Method\t(.*)", Target::SoftwareName),
    ("^Instrument Method\t(.*)", Target::InstrumentName),
    ("^Injection Date\t(.*)", Target::Meta("injection_date")),
    ("^Injection Time\t(.*)", Target::Meta("injection_time")),
    ("^Detector\t(.*)", Target::Meta("detector")),
    ("^Signal Quantity\t(.*)", Target::Meta("signal_quantity")),
    ("^Signal Unit\t(.*)", Target::Meta("signal_unit")),
    ("^Signal Info\t(.*)", Target::Meta("signal_info")),
];

fn compile(pattern: &str) -> Result<BoostRegex> {
    BoostRegex::with_options(
        pattern,
        RegexOptions {
            no_mod_s: true,
            ..RegexOptions::default()
        },
    )
}

impl ChromeleonFile {
    /// A stateless loader (source default constructor).
    pub fn new() -> Self {
        Self
    }

    /// Load the file's data and metadata into a new experiment.
    ///
    /// Equivalent to [`Self::load_with_limits`] with [`ReadLimits::default`].
    ///
    /// # Errors
    ///
    /// See [`Self::load_with_limits`].
    pub fn load(&self, path: impl AsRef<Path>) -> Result<MSExperiment> {
        self.load_with_limits(path, ReadLimits::default())
    }

    /// Load the file at `path` (the Chromeleon input) into a new experiment.
    ///
    /// Header lines are matched against the whole line, first match wins, in
    /// this order; a later line with the same key overwrites an earlier one:
    ///
    /// | line | stored as |
    /// |---|---|
    /// | `Injection<TAB>v` | meta value `mzml_id` |
    /// | `Channel<TAB>v` | meta value `acq_method_name` |
    /// | `Processing Method<TAB>v` | instrument software name |
    /// | `Instrument Method<TAB>v` | instrument name |
    /// | `Injection Date<TAB>v` | meta value `injection_date` |
    /// | `Injection Time<TAB>v` | meta value `injection_time` |
    /// | `Detector<TAB>v` | meta value `detector` |
    /// | `Signal Quantity<TAB>v` | meta value `signal_quantity` |
    /// | `Signal Unit<TAB>v` | meta value `signal_unit` |
    /// | `Signal Info<TAB>v` | meta value `signal_info` |
    ///
    /// The value `v` must not contain a line separator (`\n`, `\r`, `\f`),
    /// which `.` does not match under `no_mod_s`; a line whose value contains
    /// a form feed is therefore ignored, as in the source. Meta values are
    /// strings.
    ///
    /// A line that is exactly `Raw Data:` or `Chromatogram Data:` ends the
    /// header; the line after it (the column caption) is skipped. Every
    /// following line is split on tabs: a line of exactly three fields adds a
    /// point whose position is the first field and whose intensity is the third
    /// field narrowed to `f32`, each after
    /// [`Self::remove_commas_and_parse_double`]; the second field (the step) is
    /// not read. Empty lines are skipped. The experiment receives exactly one
    /// chromatogram, empty when the file has no data marker.
    ///
    /// The position is the file's time value as written, in **minutes** for
    /// the `Time (min)` column Chromeleon exports; the source stores it without
    /// converting to seconds, and so does this port.
    ///
    /// # Errors
    ///
    /// - [`Error::Io`] when the file cannot be opened or read; the source
    ///   throws `Exception::FileNotFound`, `Exception::FileNotReadable` or
    ///   `Exception::IOException`, which here are the `NotFound`,
    ///   `PermissionDenied` and other I/O error kinds.
    /// - [`Error::Parse`] for a non-empty data line without exactly three
    ///   tab-separated fields (source `Exception::ParseError`, "Couldn't parse
    ///   the raw data."), and for a field that does not convert to a number
    ///   (source `Exception::ConversionError`).
    /// - [`Error::Parse`] when a stored header value is not UTF-8; the source
    ///   stores the bytes.
    /// - [`Error::Parse`], a native check, when a point is NaN or infinite or
    ///   its intensity overflows `f32`; the source stores such a point.
    /// - [`Error::InvalidValue`] when the input exceeds `limits`.
    pub fn load_with_limits(
        &self,
        path: impl AsRef<Path>,
        limits: ReadLimits,
    ) -> Result<MSExperiment> {
        let file = File::open(path.as_ref())?;
        self.read_with_limits(BufReader::new(file), limits)
    }

    /// Reads Chromeleon text from `reader`.
    ///
    /// # Errors
    ///
    /// See [`Self::load_with_limits`].
    pub fn read(&self, reader: impl BufRead) -> Result<MSExperiment> {
        self.read_with_limits(reader, ReadLimits::default())
    }

    /// Reads Chromeleon text from `reader` under explicit `limits`.
    ///
    /// # Errors
    ///
    /// See [`Self::load_with_limits`].
    pub fn read_with_limits(
        &self,
        reader: impl BufRead,
        limits: ReadLimits,
    ) -> Result<MSExperiment> {
        let header = HEADER
            .iter()
            .map(|(pattern, target)| Ok((compile(pattern)?, *target)))
            .collect::<Result<Vec<_>>>()?;
        let raw_data = compile("^Raw Data:")?;
        let chromatogram_data = compile("^Chromatogram Data:")?;

        let mut experiment = MSExperiment::new();
        let mut chromatogram = MSChromatogram::default();
        let mut lines = LineReader::new(reader, limits.lines);
        let mut line = Vec::new();

        'header: while lines.next_bytes(&mut line)? {
            for (regex, target) in &header {
                let Some(captures) = regex.full_match(&line)? else {
                    continue;
                };
                let range = captures.get(1).unwrap_or(0..0);
                let value = std::str::from_utf8(&line[range]).map_err(|_| Error::Parse {
                    line: lines.line_number(),
                    message: "header value is not valid UTF-8".into(),
                })?;
                match target {
                    Target::Meta(key) => {
                        experiment
                            .settings
                            .metadata
                            .insert((*key).to_owned(), MetaValue::from(value));
                    }
                    Target::SoftwareName => {
                        experiment.settings.instrument.software.name = value.to_owned();
                    }
                    Target::InstrumentName => {
                        experiment.settings.instrument.name = value.to_owned();
                    }
                }
                continue 'header;
            }
            if raw_data.is_full_match(&line)? || chromatogram_data.is_full_match(&line)? {
                // Remove the subsequent line, right before the raw data.
                lines.next_bytes(&mut line)?;
                break;
            }
        }

        let mut text = String::new();
        while lines.next_line(&mut text)? {
            let fields: Vec<&str> = if text.is_empty() {
                Vec::new()
            } else {
                text.split('\t').collect()
            };
            if fields.len() == 3 {
                let number = lines.line_number();
                if chromatogram.peaks.len() >= limits.max_points {
                    return Err(Error::InvalidValue(format!(
                        "Chromeleon chromatogram exceeds {} points",
                        limits.max_points
                    )));
                }
                let rt = self.parse_field(fields[0], number)?;
                let intensity = self.parse_field(fields[2], number)? as f32;
                if !rt.is_finite() || !intensity.is_finite() {
                    return Err(Error::Parse {
                        line: number,
                        message: format!(
                            "raw data point is NaN, infinite or outside the f32 intensity range: '{text}'"
                        ),
                    });
                }
                chromatogram
                    .peaks
                    .try_reserve(1)
                    .map_err(|_| Error::InvalidValue("chromatogram allocation failed".into()))?;
                chromatogram
                    .peaks
                    .push(ChromatogramPeak::new(rt, intensity));
            } else if text.is_empty() {
                // Skips eventual empty lines, e.g. the last before EOF.
                continue;
            } else {
                return Err(Error::Parse {
                    line: lines.line_number(),
                    message: format!("{text}: Couldn't parse the raw data."),
                });
            }
        }
        experiment.chromatograms.push(chromatogram);
        Ok(experiment)
    }

    fn parse_field(&self, field: &str, line: usize) -> Result<f64> {
        let mut number = field.to_owned();
        self.remove_commas_and_parse_double(&mut number)
            .map_err(|error| match error {
                Error::InvalidValue(message) => Error::Parse { line, message },
                other => other,
            })
    }

    /// Remove commas from the string (used as thousands separators) and parse
    /// its value.
    ///
    /// As in the source, the commas are removed from `number` itself, which
    /// keeps them removed even when the conversion fails. The conversion is
    /// the source's `StringUtils::toDouble`: surrounding space, tab, CR and LF
    /// are ignored and the rest must be one complete number.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] carrying the source's
    /// `Exception::ConversionError` message when the remaining text is not a
    /// number.
    pub fn remove_commas_and_parse_double(&self, number: &mut String) -> Result<f64> {
        number.retain(|c| c != ',');
        to_double(number).map_err(Error::InvalidValue)
    }
}
