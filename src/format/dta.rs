// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! DTA singly protonated precursor mass and fragment peak lists.
//! Read behavior follows `FORMAT/DTAFile.h`. The default writer uses the exact
//! proton mass in both directions; the C++ writer's legacy 1.0 approximation
//! is available explicitly through [`MassConvention::LegacyOpenMS`].

use super::{intensity, number, parse_error};
use crate::chemistry::PROTON_MASS_U;
use crate::{Error, MSSpectrum, Peak1D, Precursor, Result};
use std::io::{BufRead, Write};

/// Precursor-mass convention used when writing DTA.
#[derive(Clone, Copy, Debug, Default)]
pub enum MassConvention {
    /// Inverse of the reader, using OpenMS's exact proton mass.
    #[default]
    ExactProton,
    /// Reproduce the pinned C++ writer's 1.0 Da approximation.
    LegacyOpenMS,
}

/// Read a DTA spectrum (MS level 2), accepting spaces or tabs.
pub fn read(reader: impl BufRead) -> Result<MSSpectrum> {
    let mut lines = reader.lines();
    let header = lines
        .next()
        .ok_or_else(|| parse_error(1, "missing DTA header"))??;
    let fields: Vec<&str> = header.split_whitespace().collect();
    if fields.len() != 2 {
        return Err(parse_error(1, "DTA header requires MH+ mass and charge"));
    }
    let mh = number(fields[0], 1, "MH+ mass")?;
    let charge = fields[1]
        .parse::<i32>()
        .map_err(|_| parse_error(1, "invalid precursor charge"))?;
    let mz = if charge == 0 {
        mh
    } else {
        (mh - PROTON_MASS_U) / f64::from(charge) + PROTON_MASS_U
    };
    let mut spectrum = MSSpectrum {
        ms_level: 2,
        precursors: vec![Precursor::new(mz, charge)],
        ..MSSpectrum::default()
    };
    for (i, line) in lines.enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() != 2 {
            return Err(parse_error(i + 2, "DTA peak requires m/z and intensity"));
        }
        spectrum.peaks.push(Peak1D::new(
            number(fields[0], i + 2, "m/z")?,
            intensity(fields[1], i + 2)?,
        ));
    }
    spectrum.validate()?;
    Ok(spectrum)
}

/// Write a spectrum with an exactly reversible proton-mass conversion.
pub fn write(writer: impl Write, spectrum: &MSSpectrum) -> Result<()> {
    write_with_convention(writer, spectrum, MassConvention::ExactProton)
}

/// How much of a record the writer may discard.
///
/// The native default refuses to drop anything DTA cannot represent, so a
/// caller cannot lose metadata without saying so. The source `DTAFile::store`
/// has no such guard: it writes the precursor mass and charge plus the peaks
/// and silently ignores everything else, warning only that it used the first of
/// several precursors. `discard_unrepresentable` selects that source behavior,
/// which a TOPP tool reproducing C++ output needs.
#[derive(Clone, Copy, Debug)]
pub struct WriteOptions {
    pub convention: MassConvention,
    pub discard_unrepresentable: bool,
}
impl Default for WriteOptions {
    fn default() -> Self {
        Self {
            convention: MassConvention::ExactProton,
            discard_unrepresentable: false,
        }
    }
}
impl WriteOptions {
    /// The source `DTAFile::store` behavior: legacy proton mass, discard the rest.
    pub fn source() -> Self {
        Self {
            convention: MassConvention::LegacyOpenMS,
            discard_unrepresentable: true,
        }
    }
}

/// Write a DTA peak list under explicit options.
pub fn write_with_options(
    writer: impl Write,
    spectrum: &MSSpectrum,
    options: &WriteOptions,
) -> Result<()> {
    write_inner(
        writer,
        spectrum,
        options.convention,
        options.discard_unrepresentable,
    )
}

/// Write a DTA peak list using an explicit precursor-mass convention.
/// More than one precursor is rejected; DTA cannot represent that information.
pub fn write_with_convention(
    writer: impl Write,
    spectrum: &MSSpectrum,
    convention: MassConvention,
) -> Result<()> {
    write_inner(writer, spectrum, convention, false)
}

fn write_inner(
    mut writer: impl Write,
    spectrum: &MSSpectrum,
    convention: MassConvention,
    discard: bool,
) -> Result<()> {
    spectrum.validate()?;
    if !discard && !spectrum.metadata.is_empty() {
        return Err(Error::Unsupported(
            "DTA cannot store spectrum metadata".into(),
        ));
    }
    if !discard && spectrum.precursors.len() > 1 {
        return Err(Error::Unsupported("DTA supports only one precursor".into()));
    }
    if !discard && !spectrum.peptide_identifications.is_empty() {
        return Err(Error::Unsupported(
            "DTA cannot store peptide identifications".into(),
        ));
    }
    let precursor = spectrum.precursors.first().cloned().unwrap_or_default();
    if !discard && precursor.has_acquisition_metadata() {
        return Err(Error::Unsupported(
            "DTA cannot store precursor acquisition metadata".into(),
        ));
    }
    let proton = match convention {
        MassConvention::ExactProton => PROTON_MASS_U,
        MassConvention::LegacyOpenMS => 1.0,
    };
    let mh = if precursor.charge == 0 {
        precursor.mz
    } else {
        (precursor.mz - proton) * f64::from(precursor.charge) + proton
    };
    if !mh.is_finite() {
        return Err(Error::InvalidValue("DTA precursor mass overflow".into()));
    }
    writeln!(writer, "{mh} {}", precursor.charge)?;
    for peak in &spectrum.peaks {
        writeln!(writer, "{} {}", peak.mz, peak.intensity)?;
    }
    writer.flush()?;
    Ok(())
}
