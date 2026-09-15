// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! DTA singly protonated precursor mass and fragment peak lists.
//! Read behavior follows `FORMAT/DTAFile.h`. The default writer uses the exact
//! proton mass in both directions; the C++ writer's legacy 1.0 approximation
//! is available explicitly through [`MassConvention::LegacyOpenMS`](crate::format::dta::MassConvention::LegacyOpenMS).
//!
//! # Numeric text
//!
//! `DTAFile::store` sets `os.precision(writtenDigits<double>(0.0))`, which is
//! 15, and then writes three different kinds of number through two different
//! formatters. Both are reproduced here from
//! [`crate::format::file_info::text_format`], which already ports them:
//!
//! | Field | Source expression | Formatter |
//! |---|---|---|
//! | precursor `MH+` mass | `os << ((mz - 1.0) * charge + 1.0)`, a `double` in the stream's default float field | [`ostream_g`] at [`WRITTEN_DIGITS_F64`] — 15 *significant* digits, `%.15g` |
//! | peak m/z | `os << it->getPosition()`, a `DPosition<1>` whose `operator<<` calls `precisionWrapper` (`DPosition.h:412-420`), so `StringUtils::toStr(double, true)` | [`to_str`] — [`NumericFormatting::appendNumeric`] with 15 *fraction* digits |
//! | peak intensity | `os << it->getIntensity()`, a `float` promoted to `double` in the same default float field | [`ostream_g`] at [`WRITTEN_DIGITS_F64`] |
//!
//! The two rules differ: `to_str` writes 15 digits *after the decimal point*
//! and `ostream_g` 15 *significant* digits, so one peak line carries both
//! `104.115715026855469` (m/z) and `260.789154052734` (intensity) for the same
//! number of source digits. The port wrote Rust's shortest round-trip text for
//! both until this was corrected, which is about 8 significant digits for an
//! `f32` intensity and produced files roughly 21% smaller than the C++ tool's.
//!
//! # Round-trip precision
//!
//! This is the source's text, so it bounds what a write-then-read recovers, and
//! [`MassConvention::ExactProton`](crate::format::dta::MassConvention::ExactProton)
//! changes the arithmetic but not the text:
//!
//! - A peak intensity is an `f32` and 15 significant digits always read back as
//!   the same `f32`.
//! - A peak m/z of 10 or more reads back as the same `f64`; 15 fraction digits
//!   resolve less than the gap between neighbouring `f64` values there. Below
//!   that the text is shorter than the value and low bits are lost, so an m/z
//!   under 10 — which no fragment spectrum carries — round-trips only to about
//!   1e-15 absolute.
//! - The precursor mass keeps 15 significant digits, roughly 1e-13 relative, so
//!   a precursor m/z survives to that precision rather than bit-exactly. The
//!   C++ `store`/`load` pair is no more exact.
//!
//! [`ostream_g`]: crate::format::file_info::text_format::ostream_g
//! [`to_str`]: crate::format::file_info::text_format::to_str
//! [`WRITTEN_DIGITS_F64`]: crate::format::file_info::text_format::WRITTEN_DIGITS_F64
//! [`NumericFormatting::appendNumeric`]: crate::format::file_info::text_format

use super::{intensity, number, parse_error};
use crate::chemistry::PROTON_MASS_U;
use crate::format::file_info::text_format::{WRITTEN_DIGITS_F64, ostream_g, to_str};
use crate::{Error, MSSpectrum, Peak1D, Precursor, Result};
use std::io::{BufRead, Write};

/// Precursor-mass convention used when writing DTA.
///
/// This selects the arithmetic of the header line only. The numeric text is the
/// source's in both variants; see the module header.
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

/// Write a spectrum with a reversible proton-mass conversion.
///
/// The arithmetic is the exact inverse of [`read`]; the text is the source's,
/// so the recovered precursor m/z is exact to the 15 significant digits the
/// header line carries. See the module's round-trip precision section.
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
    // Source text, as the module header sets out: the header mass and every
    // intensity through the stream's default float field at precision 15, every
    // m/z through `precisionWrapper`. A line is assembled in one reused buffer
    // so that a spectrum of n peaks costs n writes, not 4n.
    let mut line = String::new();
    line.push_str(&ostream_g(mh, WRITTEN_DIGITS_F64));
    line.push(' ');
    line.push_str(&precursor.charge.to_string());
    line.push('\n');
    writer.write_all(line.as_bytes())?;
    for peak in &spectrum.peaks {
        line.clear();
        line.push_str(&to_str(peak.mz));
        line.push(' ');
        line.push_str(&ostream_g(f64::from(peak.intensity), WRITTEN_DIGITS_F64));
        line.push('\n');
        writer.write_all(line.as_bytes())?;
    }
    writer.flush()?;
    Ok(())
}
