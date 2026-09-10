// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Three-column DTA2D spectra and explicit MS1 TIC projection.

pub use super::ms2::Limits;
use super::ms2::{
    Counter, TextInput, check_experiment, check_peaks, check_spectrum_payload, fields, increment,
    invalid, parse_intensity, parse_number, push_peak, push_spectrum, trim, unsupported,
};
use super::parse_error;
use crate::{MSExperiment, MSSpectrum, Peak1D, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::ops::Range;
use std::path::Path;
pub type WriteOptions = Limits;

/// The three source-consumed PeakFileOptions filters. Ranges are half-open,
/// including the minimum and excluding the maximum, as in OpenMS DRange.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadOptions {
    pub rt_range: Option<Range<f64>>,
    pub mz_range: Option<Range<f64>>,
    pub intensity_range: Option<Range<f64>>,
    pub limits: Limits,
}
impl ReadOptions {
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
pub fn read(reader: impl BufRead) -> Result<MSExperiment> {
    read_with_options(reader, &ReadOptions::default())
}
pub fn read_with_options(reader: impl BufRead, options: &ReadOptions) -> Result<MSExperiment> {
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
fn preflight(experiment: &MSExperiment, options: &Limits) -> Result<()> {
    check_experiment(experiment, options)?;
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
        check_peaks(s, &mut peaks, options)?;
    }
    render(&mut Counter(options.max_output_bytes), experiment)
}
fn render(writer: &mut impl Write, experiment: &MSExperiment) -> Result<()> {
    writeln!(writer, "#SEC\tMZ\tINT")?;
    for s in &experiment.spectra {
        for p in &s.peaks {
            writeln!(writer, "{}\t{}\t{}", s.rt, p.mz, p.intensity)?;
        }
    }
    Ok(())
}
/// Explicit lossy projection: MS1 TIC only, one row per scan (including empty
/// scans), m/z=0 and f32 accumulation in storage order. Other data are unconsumed.
pub fn write_tic(writer: impl Write, experiment: &MSExperiment) -> Result<()> {
    write_tic_with_options(writer, experiment, &WriteOptions::default())
}
pub fn write_tic_with_options(
    mut writer: impl Write,
    experiment: &MSExperiment,
    options: &WriteOptions,
) -> Result<()> {
    let points = tic_points(experiment, options)?;
    render_tic(&mut Counter(options.max_output_bytes), &points)?;
    render_tic(&mut writer, &points)?;
    writer.flush()?;
    Ok(())
}
pub fn store_tic(path: impl AsRef<Path>, experiment: &MSExperiment) -> Result<()> {
    store_tic_with_options(path, experiment, &WriteOptions::default())
}
pub fn store_tic_with_options(
    path: impl AsRef<Path>,
    experiment: &MSExperiment,
    options: &WriteOptions,
) -> Result<()> {
    let points = tic_points(experiment, options)?;
    render_tic(&mut Counter(options.max_output_bytes), &points)?;
    let mut writer = BufWriter::new(File::create(path)?);
    render_tic(&mut writer, &points)?;
    writer.flush()?;
    Ok(())
}
fn tic_points(experiment: &MSExperiment, options: &Limits) -> Result<Vec<(f64, f32)>> {
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
        points.push((s.rt, sum));
    }
    Ok(points)
}
fn render_tic(writer: &mut impl Write, points: &[(f64, f32)]) -> Result<()> {
    writeln!(writer, "#SEC\tMZ\tINT")?;
    for &(rt, intensity) in points {
        writeln!(writer, "{rt}\t0\t{intensity}")?;
    }
    Ok(())
}
