// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Source MS2 peak-list reader and a native writer for its representable subset.
//! H/I/Z/D records are intentionally ignored, as in OpenMS MS2File::load.

use super::{number, parse_error};
use crate::{Error, MSExperiment, MSSpectrum, Peak1D, Precursor, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

/// Hard-bounded text adapter allowances. Values may be lowered from defaults.
/// Peak/scan/line counts include discarded or filtered input. Output bytes are
/// counted before writing or creating a destination file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Limits {
    pub max_bytes: usize,
    pub max_line_bytes: usize,
    pub max_lines: usize,
    pub max_spectra: usize,
    pub max_peaks: usize,
    pub max_output_bytes: usize,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            max_bytes: 256 * 1024 * 1024,
            max_line_bytes: 1024 * 1024,
            max_lines: 20_000_000,
            max_spectra: 100_000,
            max_peaks: 10_000_000,
            max_output_bytes: 512 * 1024 * 1024,
        }
    }
}
impl Limits {
    pub(super) fn validate(&self) -> Result<()> {
        let ceiling = Self::default();
        if self.max_bytes > ceiling.max_bytes
            || self.max_line_bytes > ceiling.max_line_bytes
            || self.max_lines > ceiling.max_lines
            || self.max_spectra > ceiling.max_spectra
            || self.max_peaks > ceiling.max_peaks
            || self.max_output_bytes > ceiling.max_output_bytes
        {
            return Err(invalid("text adapter limits exceed hard ceilings"));
        }
        Ok(())
    }
}
pub type ReadOptions = Limits;
pub type WriteOptions = Limits;

pub fn read(reader: impl BufRead) -> Result<MSExperiment> {
    read_with_options(reader, &ReadOptions::default())
}
pub fn read_with_options(reader: impl BufRead, options: &ReadOptions) -> Result<MSExperiment> {
    let mut input = TextInput::new(reader, options)?;
    let mut result = MSExperiment::default();
    let mut current: Option<MSSpectrum> = None;
    let (mut scans, mut peaks) = (0, 0);
    let mut line = String::new();
    while input.next_line(&mut line)? {
        let text = trim(&line);
        let Some(first) = text.as_bytes().first() else {
            continue;
        };
        match first {
            b'H' | b'I' | b'Z' | b'D' => continue,
            b'S' => {
                increment(
                    &mut scans,
                    options.max_spectra,
                    "MS2 spectrum limit exceeded",
                )?;
                let fields =
                    fields::<4>(text.split(is_space).filter(|s| !s.is_empty()), input.line)?;
                let mz = parse_number(fields[3], input.line, "precursor m/z")?;
                if let Some(spectrum) = current.take() {
                    push_spectrum(&mut result, spectrum)?;
                }
                current = Some(MSSpectrum {
                    ms_level: 2,
                    native_id: format!("index={}", scans - 1),
                    precursors: vec![Precursor::new(mz, 0)],
                    ..Default::default()
                });
            }
            _ => {
                increment(&mut peaks, options.max_peaks, "MS2 peak limit exceeded")?;
                let fields =
                    fields::<2>(text.split(is_space).filter(|s| !s.is_empty()), input.line)?;
                let peak = Peak1D::new(
                    parse_number(fields[0], input.line, "m/z")?,
                    parse_intensity(fields[1], input.line)?,
                );
                // The source discards pre-S peaks when the first scan starts.
                if let Some(spectrum) = &mut current {
                    push_peak(spectrum, peak)?;
                }
            }
        }
    }
    if let Some(spectrum) = current {
        push_spectrum(&mut result, spectrum)?;
    }
    Ok(result)
}
/// Replacement is atomic on parse, resource or stream failure.
pub fn read_into(reader: impl BufRead, destination: &mut MSExperiment) -> Result<()> {
    *destination = read(reader)?;
    Ok(())
}
pub fn read_into_with_options(
    reader: impl BufRead,
    destination: &mut MSExperiment,
    options: &ReadOptions,
) -> Result<()> {
    *destination = read_with_options(reader, options)?;
    Ok(())
}
pub fn load(path: impl AsRef<Path>) -> Result<MSExperiment> {
    read(BufReader::new(File::open(path)?))
}
pub fn load_with_options(path: impl AsRef<Path>, options: &ReadOptions) -> Result<MSExperiment> {
    read_with_options(BufReader::new(File::open(path)?), options)
}

/// Native extension: emits S and peak records. Unsupported metadata is rejected
/// before any bytes; the source MS2 adapter itself has no store method.
pub fn write(writer: impl Write, experiment: &MSExperiment) -> Result<()> {
    write_with_options(writer, experiment, &WriteOptions::default())
}
pub fn write_with_options(
    mut writer: impl Write,
    experiment: &MSExperiment,
    options: &WriteOptions,
) -> Result<()> {
    preflight(experiment, options)?;
    render(&mut writer, experiment)?;
    writer.flush()?;
    Ok(())
}
pub fn store(path: impl AsRef<Path>, experiment: &MSExperiment) -> Result<()> {
    store_with_options(path, experiment, &WriteOptions::default())
}
pub fn store_with_options(
    path: impl AsRef<Path>,
    experiment: &MSExperiment,
    options: &WriteOptions,
) -> Result<()> {
    preflight(experiment, options)?;
    let mut writer = BufWriter::new(File::create(path)?);
    render(&mut writer, experiment)?;
    writer.flush()?;
    Ok(())
}
fn preflight(experiment: &MSExperiment, options: &WriteOptions) -> Result<()> {
    check_experiment(experiment, options)?;
    let mut peaks = 0;
    for spectrum in &experiment.spectra {
        check_spectrum_payload(spectrum)?;
        if spectrum.ms_level != 2 || spectrum.rt != -1.0 || spectrum.precursors.len() != 1 {
            return Err(unsupported(
                "MS2 writer requires MS2, unset RT and exactly one precursor",
            ));
        }
        let precursor = &spectrum.precursors[0];
        if precursor.charge != 0
            || precursor.intensity != 0.0
            || precursor.has_acquisition_metadata()
        {
            return Err(unsupported(
                "MS2 source reader cannot retain precursor charge, intensity or acquisition metadata",
            ));
        }
        if !precursor.mz.is_finite() {
            return Err(invalid("nonfinite precursor m/z"));
        }
        check_peaks(spectrum, &mut peaks, options)?;
    }
    render(&mut Counter(options.max_output_bytes), experiment)
}
fn render(writer: &mut impl Write, experiment: &MSExperiment) -> Result<()> {
    for (i, spectrum) in experiment.spectra.iter().enumerate() {
        writeln!(
            writer,
            "S\t{}\t{}\t{}",
            i + 1,
            i + 1,
            spectrum.precursors[0].mz
        )?;
        for peak in &spectrum.peaks {
            writeln!(writer, "{}\t{}", peak.mz, peak.intensity)?;
        }
    }
    Ok(())
}

// Shared narrow mechanics for the two flat text adapters; no generic I/O model.
pub(super) fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
pub(super) fn unsupported(message: &str) -> Error {
    Error::Unsupported(message.into())
}
pub(super) fn is_space(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\r')
}
pub(super) fn trim(text: &str) -> &str {
    text.trim_matches(is_space)
}
fn nonzero_mantissa(text: &str) -> bool {
    text.split(['e', 'E'])
        .next()
        .unwrap_or("")
        .bytes()
        .any(|b| matches!(b, b'1'..=b'9'))
}
pub(super) fn parse_number(text: &str, line: usize, label: &str) -> Result<f64> {
    let value = number(trim(text), line, label)?;
    if value == 0.0 && nonzero_mantissa(text) {
        return Err(parse_error(line, "nonzero numeric value underflows f64"));
    }
    Ok(value)
}
pub(super) fn parse_intensity(text: &str, line: usize) -> Result<f32> {
    // Source toFloat parses directly to binary32, avoiding double rounding.
    trim(text)
        .parse::<f32>()
        .ok()
        .filter(|x| x.is_finite() && (*x != 0.0 || !nonzero_mantissa(text)))
        .ok_or_else(|| parse_error(line, "invalid finite f32 intensity or underflow"))
}
pub(super) fn fields<'a, const N: usize>(
    mut fields: impl Iterator<Item = &'a str>,
    line: usize,
) -> Result<[&'a str; N]> {
    let mut result = [""; N];
    for field in &mut result {
        *field = fields
            .next()
            .ok_or_else(|| parse_error(line, format!("expected {N} fields")))?;
    }
    if fields.next().is_some() {
        return Err(parse_error(line, format!("expected {N} fields")));
    }
    Ok(result)
}
pub(super) fn increment(value: &mut usize, maximum: usize, message: &str) -> Result<()> {
    if *value >= maximum {
        return Err(invalid(message));
    }
    *value += 1;
    Ok(())
}
pub(super) fn push_peak(spectrum: &mut MSSpectrum, peak: Peak1D) -> Result<()> {
    spectrum
        .peaks
        .try_reserve(1)
        .map_err(|_| invalid("peak allocation failed"))?;
    spectrum.peaks.push(peak);
    Ok(())
}
pub(super) fn push_spectrum(experiment: &mut MSExperiment, spectrum: MSSpectrum) -> Result<()> {
    experiment
        .spectra
        .try_reserve(1)
        .map_err(|_| invalid("spectrum allocation failed"))?;
    experiment.spectra.push(spectrum);
    Ok(())
}
pub(super) struct TextInput<R> {
    reader: R,
    options: Limits,
    pub line: usize,
    bytes: usize,
}
impl<R: BufRead> TextInput<R> {
    pub fn new(reader: R, options: &Limits) -> Result<Self> {
        options.validate()?;
        Ok(Self {
            reader,
            options: *options,
            line: 0,
            bytes: 0,
        })
    }
    pub fn next_line(&mut self, line: &mut String) -> Result<bool> {
        line.clear();
        let allowance = self.options.max_bytes - self.bytes;
        let cap = allowance.min(self.options.max_line_bytes) + 1;
        let read = self.reader.by_ref().take(cap as u64).read_line(line)?;
        if read == 0 {
            return Ok(false);
        }
        if read > allowance || read > self.options.max_line_bytes {
            return Err(invalid("text byte or line-length limit exceeded"));
        }
        increment(
            &mut self.line,
            self.options.max_lines,
            "text line-count limit exceeded",
        )?;
        self.bytes += read;
        Ok(true)
    }
}
pub(super) struct Counter(pub usize);
impl Write for Counter {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_sub(bytes.len())
            .ok_or_else(|| std::io::Error::other("text output byte limit exceeded"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
pub(super) fn check_experiment(experiment: &MSExperiment, options: &Limits) -> Result<()> {
    options.validate()?;
    if experiment.spectra.len() > options.max_spectra {
        return Err(invalid("spectrum limit exceeded"));
    }
    if !experiment.chromatograms.is_empty() || experiment.settings.has_transport_metadata() {
        return Err(unsupported(
            "flat peak list cannot store experiment metadata or chromatograms",
        ));
    }
    Ok(())
}
pub(super) fn check_spectrum_payload(s: &MSSpectrum) -> Result<()> {
    if s.native_id.len() > 1024 {
        return Err(invalid("native identifier exceeds text adapter limit"));
    }
    if !s.name.is_empty()
        || s.spectrum_type != crate::kernel::SpectrumType::Unknown
        || !s.metadata.is_empty()
        || !s.peptide_identifications.is_empty()
        || !s.float_data_arrays.is_empty()
        || !s.integer_data_arrays.is_empty()
        || !s.string_data_arrays.is_empty()
    {
        return Err(unsupported(
            "flat peak list cannot store spectrum annotations or metadata",
        ));
    }
    // Source-generated index IDs are bookkeeping and may be regenerated after
    // filtering; an unrelated identifier must not be silently discarded.
    if !s.native_id.is_empty()
        && !s
            .native_id
            .strip_prefix("index=")
            .is_some_and(|v| !v.is_empty() && v.bytes().all(|b| b.is_ascii_digit()))
    {
        return Err(unsupported(
            "flat peak list cannot retain this native identifier",
        ));
    }
    Ok(())
}
pub(super) fn check_peaks(s: &MSSpectrum, total: &mut usize, options: &Limits) -> Result<()> {
    if s.peaks.len() > options.max_peaks.saturating_sub(*total) {
        return Err(invalid("peak limit exceeded"));
    }
    *total += s.peaks.len();
    if !s.rt.is_finite()
        || s.peaks
            .iter()
            .any(|p| !p.mz.is_finite() || !p.intensity.is_finite())
    {
        return Err(invalid("nonfinite peak-list coordinate or intensity"));
    }
    Ok(())
}
