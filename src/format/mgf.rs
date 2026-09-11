// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! MGF peak-list I/O for TITLE, PEPMASS, CHARGE, RTINSECONDS and MSLEVEL.
//! Additional unique key/value fields (including SCANS) remain in metadata.
//! Repeated fields, ambiguous charges, RT ranges and peak annotations are
//! rejected explicitly. Mascot search submission and repeated SEQ are not ported.

use super::{intensity, number, parse_error, single_line};
use crate::{Error, MSExperiment, MSSpectrum, Peak1D, Precursor, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, Write};

/// Streaming MGF iterator. Global parameters are inherited by subsequent ions.
pub struct MgfReader<R> {
    reader: R,
    line: usize,
    index: usize,
    finished: bool,
    header: BTreeMap<String, String>,
}

impl<R: BufRead> MgfReader<R> {
    /// Wrap a buffered stream without reading it yet.
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line: 0,
            index: 0,
            finished: false,
            header: BTreeMap::new(),
        }
    }

    fn read_spectrum(&mut self) -> Result<Option<MSSpectrum>> {
        let mut current: Option<MSSpectrum> = None;
        let mut fields = BTreeMap::new();
        loop {
            let mut text = String::new();
            if self.reader.read_line(&mut text)? == 0 {
                return if current.is_some() {
                    Err(parse_error(self.line, "missing END IONS"))
                } else {
                    Ok(None)
                };
            }
            self.line += 1;
            let text = text.trim_start_matches('\u{feff}').trim();
            if text.is_empty() || text.starts_with(['#', ';', '!', '/']) {
                continue;
            }
            if text == "BEGIN IONS" {
                if current.is_some() {
                    return Err(parse_error(self.line, "nested BEGIN IONS"));
                }
                current = Some(MSSpectrum {
                    ms_level: 2,
                    native_id: format!("index={}", self.index),
                    ..MSSpectrum::default()
                });
                continue;
            }
            if text == "END IONS" {
                let mut spectrum = current
                    .take()
                    .ok_or_else(|| parse_error(self.line, "END IONS without BEGIN IONS"))?;
                let mut merged = self.header.clone();
                merged.extend(fields);
                apply_fields(&mut spectrum, merged, self.line)?;
                spectrum.validate()?;
                self.index += 1;
                return Ok(Some(spectrum));
            }
            if let Some((key, value)) = text.split_once('=') {
                let key = key.trim().to_ascii_uppercase();
                if key.is_empty() || key.chars().any(char::is_whitespace) {
                    return Err(parse_error(self.line, "invalid MGF field name"));
                }
                let target = if current.is_some() {
                    &mut fields
                } else {
                    &mut self.header
                };
                if target
                    .insert(key.clone(), value.trim().to_owned())
                    .is_some()
                {
                    return Err(Error::Unsupported(format!(
                        "repeated MGF field {key} on line {}",
                        self.line
                    )));
                }
            } else {
                let spectrum = current
                    .as_mut()
                    .ok_or_else(|| parse_error(self.line, "peak outside BEGIN IONS block"))?;
                let parts: Vec<&str> = text.split_whitespace().collect();
                if parts.len() != 2 {
                    return Err(parse_error(
                        self.line,
                        "MGF peak requires exactly m/z and intensity",
                    ));
                }
                spectrum.peaks.push(Peak1D::new(
                    number(parts[0], self.line, "m/z")?,
                    intensity(parts[1], self.line)?,
                ));
            }
        }
    }
}

impl<R: BufRead> Iterator for MgfReader<R> {
    type Item = Result<MSSpectrum>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        match self.read_spectrum() {
            Ok(Some(spectrum)) => Some(Ok(spectrum)),
            Ok(None) => {
                self.finished = true;
                None
            }
            Err(error) => {
                self.finished = true;
                Some(Err(error))
            }
        }
    }
}

fn apply_fields(
    spectrum: &mut MSSpectrum,
    mut fields: BTreeMap<String, String>,
    line: usize,
) -> Result<()> {
    spectrum.name = fields.remove("TITLE").unwrap_or_default();
    let mut precursor = Precursor::default();
    if let Some(value) = fields.remove("PEPMASS") {
        let parts: Vec<&str> = value.split_whitespace().collect();
        if parts.is_empty() || parts.len() > 2 {
            return Err(parse_error(
                line,
                "PEPMASS requires mass and optional intensity",
            ));
        }
        precursor.mz = number(parts[0], line, "precursor m/z")?;
        if parts.len() == 2 {
            precursor.intensity = intensity(parts[1], line)?;
        }
    }
    if let Some(value) = fields.remove("CHARGE") {
        let (digits, sign) = if let Some(v) = value.strip_suffix('+') {
            (v, 1)
        } else if let Some(v) = value.strip_suffix('-') {
            (v, -1)
        } else {
            (value.as_str(), 1)
        };
        // Parse the magnitude widely enough for i32::MIN written as 2147483648-.
        let charge = digits
            .parse::<i64>()
            .ok()
            .and_then(|v| v.checked_mul(sign))
            .and_then(|v| i32::try_from(v).ok())
            .ok_or_else(|| {
                Error::Unsupported(format!(
                    "MGF charge must be a single integer, got {value:?}"
                ))
            })?;
        if value.ends_with(['+', '-']) && !digits.chars().all(|c| c.is_ascii_digit()) {
            return Err(parse_error(line, "invalid signed MGF charge"));
        }
        precursor.charge = charge;
    }
    if let Some(value) = fields.remove("RTINSECONDS") {
        spectrum.rt = number(&value, line, "retention time")?;
    }
    if let Some(value) = fields.remove("MSLEVEL") {
        spectrum.ms_level = value
            .parse()
            .ok()
            .filter(|&x| x > 0)
            .ok_or_else(|| parse_error(line, "invalid MSLEVEL"))?;
    }
    spectrum.precursors.push(precursor);
    spectrum.metadata = fields
        .into_iter()
        .map(|(key, value)| (key, value.into()))
        .collect();
    Ok(())
}

/// Read an experiment's spectra from MGF.
pub fn read(reader: impl BufRead) -> Result<MSExperiment> {
    Ok(MSExperiment {
        spectra: MgfReader::new(reader).collect::<Result<_>>()?,
        ..MSExperiment::default()
    })
}

/// Write MGF peak lists. Chromatograms and multiple precursors are unsupported.
/// Native IDs and auxiliary arrays are not represented in MGF.
pub fn write(mut writer: impl Write, experiment: &MSExperiment) -> Result<()> {
    if experiment.settings.has_transport_metadata() {
        return Err(Error::Unsupported(
            "MGF cannot store experiment settings".into(),
        ));
    }
    if !experiment.chromatograms.is_empty() {
        return Err(Error::Unsupported("MGF cannot store chromatograms".into()));
    }
    for spectrum in &experiment.spectra {
        spectrum.validate()?;
        if spectrum
            .precursors
            .iter()
            .any(Precursor::has_acquisition_metadata)
        {
            return Err(Error::Unsupported(
                "MGF cannot store precursor acquisition metadata".into(),
            ));
        }
        if !spectrum.peptide_identifications.is_empty() {
            return Err(Error::Unsupported(
                "MGF writer cannot store structured peptide identifications".into(),
            ));
        }
        if spectrum.precursors.len() > 1 {
            return Err(Error::Unsupported(
                "MGF writer supports one precursor".into(),
            ));
        }
        if !single_line(&spectrum.name) {
            return Err(Error::InvalidValue(
                "MGF title contains control characters".into(),
            ));
        }
        let mut keys = BTreeSet::new();
        for (key, value) in &spectrum.metadata {
            if value.unit().is_some() {
                return Err(Error::Unsupported(
                    "MGF cannot retain metadata units".into(),
                ));
            }
            let value = value
                .as_str()
                .map_err(|_| Error::Unsupported("MGF can retain only String metadata".into()))?;
            let normalized = key.to_ascii_uppercase();
            if key.is_empty()
                || !single_line(key)
                || key.chars().any(char::is_whitespace)
                || key.starts_with(['#', ';', '!', '/'])
                || key.contains('=')
                || !single_line(value)
                || matches!(
                    normalized.as_str(),
                    "TITLE" | "PEPMASS" | "CHARGE" | "RTINSECONDS" | "MSLEVEL"
                )
                || !keys.insert(normalized)
            {
                return Err(Error::InvalidValue(format!(
                    "invalid or reserved MGF metadata key {key:?}"
                )));
            }
        }
    }
    for spectrum in &experiment.spectra {
        writeln!(writer, "BEGIN IONS")?;
        if !spectrum.name.is_empty() {
            writeln!(writer, "TITLE={}", spectrum.name)?;
        }
        if let Some(precursor) = spectrum.precursors.first() {
            writeln!(writer, "PEPMASS={} {}", precursor.mz, precursor.intensity)?;
            if precursor.charge != 0 {
                writeln!(
                    writer,
                    "CHARGE={}{}",
                    precursor.charge.unsigned_abs(),
                    if precursor.charge < 0 { '-' } else { '+' }
                )?;
            }
        }
        if spectrum.rt != -1.0 {
            writeln!(writer, "RTINSECONDS={}", spectrum.rt)?;
        }
        if spectrum.ms_level != 2 {
            writeln!(writer, "MSLEVEL={}", spectrum.ms_level)?;
        }
        for (key, value) in &spectrum.metadata {
            writeln!(writer, "{key}={}", value.as_str()?)?;
        }
        for peak in &spectrum.peaks {
            writeln!(writer, "{} {}", peak.mz, peak.intensity)?;
        }
        writeln!(writer, "END IONS\n")?;
    }
    writer.flush()?;
    Ok(())
}
