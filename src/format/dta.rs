// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! DTA singly protonated precursor mass and fragment peak lists.
//!
//! Port of `FORMAT/DTAFile.h`: the first line holds the singly protonated
//! peptide mass (MH+) and the charge state, every later line one fragment
//! peak's m/z and intensity. From the mass and the charge the precursor m/z is
//! computed and stored in the spectrum's precursor. Read behavior follows the
//! source; the default writer uses the exact proton mass in both directions,
//! and the C++ writer's legacy 1.0 approximation is available explicitly
//! through [`MassConvention::LegacyOpenMS`](crate::format::dta::MassConvention::LegacyOpenMS).
//! See `docs/TEXT_PEAK_LIST_SUPPORT.md` for the API mapping.
//!
//! # API mapping
//!
//! | Source member | Rust |
//! |---|---|
//! | `DTAFile()`, `virtual ~DTAFile()` | none: the class's one member is a constant, so the operations are free functions |
//! | `template <typename SpectrumType> void load(filename, spectrum)` | [`load`](crate::format::dta::load); the stream form is [`read`](crate::format::dta::read) |
//! | `template <typename SpectrumType> void store(filename, spectrum) const` | [`store`](crate::format::dta::store), [`store_with_options`](crate::format::dta::store_with_options); the stream forms are [`write`](crate::format::dta::write), [`write_with_options`](crate::format::dta::write_with_options), [`write_with_convention`](crate::format::dta::write_with_convention) |
//! | `protected: UInt default_ms_level_` | the constant [`DEFAULT_MS_LEVEL`](crate::format::dta::DEFAULT_MS_LEVEL): the constructor sets it to 2 (`DTAFile.cpp:16-19`) and no member changes it |
//!
//! Both templates are documented as taking "a MSSpectrum or [a type with] the
//! same interface"; every instantiation in the pinned source is `MSSpectrum`
//! (`FileHandler.cpp:796`, `:832`, `:873` and `:1155`, the last two on `exp[0]`), and this crate has one
//! spectrum type, so `SpectrumType` is [`MSSpectrum`] here.
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
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

/// MS level of every loaded spectrum, the source's protected
/// `default_ms_level_`, which `DTAFile()` sets to 2 and nothing changes.
pub const DEFAULT_MS_LEVEL: u32 = 2;

/// The whitespace `StringUtils::trim` removes and `toDouble`/`toInt32` skip.
fn is_space(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\r')
}

/// The source's split of a trimmed line into exactly two fields: on a tab when
/// the line contains one, otherwise on a single space, so that two adjacent
/// delimiters make an empty field and the line has three fields, not two.
fn two_fields(line: &str) -> Option<[&str; 2]> {
    let line = line.trim_matches(is_space);
    let delimiter = if line.contains('\t') { '\t' } else { ' ' };
    let mut fields = line.split(delimiter);
    let pair = [fields.next()?, fields.next()?];
    match fields.next() {
        Some(_) => None,
        None => Some(pair),
    }
}

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

/// Read a DTA spectrum from a stream, the body of source `DTAFile::load`
/// (`DTAFile.h:54-166`).
///
/// The first line must hold the MH+ mass and the charge. A charge of 0 keeps
/// the mass as the precursor m/z; any other charge `z` gives
/// `(MH+ - PROTON_MASS_U) / z + PROTON_MASS_U`. The spectrum has MS level
/// [`DEFAULT_MS_LEVEL`], one precursor and one peak per later non-blank line,
/// in file order. Every line is trimmed of spaces, tabs, CR and LF and split on
/// a tab if it contains one, otherwise on single spaces, as the source splits
/// it; each field may carry further surrounding whitespace, which the source's
/// number conversion skips. A stream has no name, so the spectrum's name is
/// left empty; [`load`] sets it.
///
/// # Errors
///
/// [`Parse`](crate::Error::Parse) for a first line or a peak line without
/// exactly two fields, or with a field that is not a number, where the source
/// throws `Exception::ParseError`; also for a nonfinite value or one outside
/// the `f32` intensity range, which the source would store (a native check).
/// [`Io`](crate::Error::Io) for a read failure.
pub fn read(reader: impl BufRead) -> Result<MSSpectrum> {
    let mut lines = reader.lines();
    // The source reads the first line whatever it holds; an empty stream has
    // an empty first line, which does not have two fields.
    let header = lines.next().transpose()?.unwrap_or_default();
    let [mass, charge] = two_fields(&header)
        .ok_or_else(|| parse_error(1, "DTA header requires MH+ mass and charge"))?;
    let mh = number(mass.trim_matches(is_space), 1, "MH+ mass")?;
    let charge = charge
        .trim_matches(is_space)
        .parse::<i32>()
        .map_err(|_| parse_error(1, "invalid precursor charge"))?;
    let mz = if charge == 0 {
        mh
    } else {
        (mh - PROTON_MASS_U) / f64::from(charge) + PROTON_MASS_U
    };
    let mut spectrum = MSSpectrum {
        ms_level: DEFAULT_MS_LEVEL,
        precursors: vec![Precursor::new(mz, charge)],
        ..MSSpectrum::default()
    };
    for (i, line) in lines.enumerate() {
        let line = line?;
        if line.trim_matches(is_space).is_empty() {
            continue;
        }
        let [mz, value] = two_fields(&line)
            .ok_or_else(|| parse_error(i + 2, "DTA peak requires m/z and intensity"))?;
        spectrum.peaks.push(Peak1D::new(
            number(mz.trim_matches(is_space), i + 2, "m/z")?,
            intensity(value.trim_matches(is_space), i + 2)?,
        ));
    }
    spectrum.validate()?;
    Ok(spectrum)
}

/// Load a DTA file, source `DTAFile::load` (`DTAFile.h:54-166`).
///
/// As [`read`], and the spectrum's name is the file's base name, the source's
/// `spectrum.setName(File::basename(filename))`.
///
/// # Errors
///
/// [`Io`](crate::Error::Io) when the file cannot be opened, where the source
/// throws `Exception::FileNotFound`, `FileNotReadable` or `IOException`;
/// [`InvalidValue`](crate::Error::InvalidValue) for a path that is not UTF-8;
/// otherwise as [`read`].
pub fn load(path: impl AsRef<Path>) -> Result<MSSpectrum> {
    let path = path.as_ref();
    let name = path
        .to_str()
        .ok_or_else(|| Error::InvalidValue("DTA filename is not UTF-8".into()))?;
    let mut spectrum = read(BufReader::new(std::fs::File::open(path)?))?;
    spectrum.name = crate::system::file::basename(name).to_owned();
    Ok(spectrum)
}

/// Store a spectrum in a DTA file with the default, lossless options.
///
/// # Errors
///
/// As [`store_with_options`].
pub fn store(path: impl AsRef<Path>, spectrum: &MSSpectrum) -> Result<()> {
    store_with_options(path, spectrum, &WriteOptions::default())
}

/// Store a spectrum in a DTA file, source `DTAFile::store`
/// (`DTAFile.h:176-220`); [`WriteOptions::source`] is the source's behaviour.
///
/// The text is produced and checked before the file is created, so a refused
/// spectrum leaves no file behind.
///
/// # Errors
///
/// As [`write_with_options`], plus [`Io`](crate::Error::Io) when the file
/// cannot be created or written, where the source throws
/// `Exception::UnableToCreateFile`.
pub fn store_with_options(
    path: impl AsRef<Path>,
    spectrum: &MSSpectrum,
    options: &WriteOptions,
) -> Result<()> {
    let mut text = Vec::new();
    write_with_options(&mut text, spectrum, options)?;
    std::fs::write(path, text)?;
    Ok(())
}

/// Write a spectrum with a reversible proton-mass conversion.
///
/// The arithmetic is the exact inverse of [`read`]; the text is the source's,
/// so the recovered precursor m/z is exact to the 15 significant digits the
/// header line carries. See the module's round-trip precision section.
///
/// # Errors
///
/// As [`write_with_convention`].
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
    /// Arithmetic of the header's MH+ mass.
    pub convention: MassConvention,
    /// Drop spectrum metadata, peptide identifications, precursor acquisition
    /// metadata and every precursor after the first instead of refusing them,
    /// as the source does (it warns on `std::cerr` only for the precursors).
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

/// Write a DTA peak list under explicit options, the body of source
/// `DTAFile::store` (`DTAFile.h:176-220`) with [`WriteOptions::source`].
///
/// # Errors
///
/// As [`write_with_convention`]; with `discard_unrepresentable` only the
/// value errors remain.
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
///
/// # Errors
///
/// [`Unsupported`](crate::Error::Unsupported) for spectrum metadata, peptide
/// identifications, precursor acquisition metadata or a second precursor;
/// [`InvalidValue`](crate::Error::InvalidValue) for a spectrum that fails
/// validation or an MH+ mass that overflows; [`Io`](crate::Error::Io) for a
/// write failure.
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
