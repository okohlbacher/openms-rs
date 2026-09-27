// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Read-only `acqus` parameter file of a Bruker XMass (Flex series) spectrum,
//! the source `Internal::AcqusHandler` (`FORMAT/HANDLERS/AcqusHandler.h`).
//!
//! An `acqus` file is a JCAMP-DX-like text file of `##key= value` lines. It
//! carries the time-of-flight calibration that turns a point index of the
//! companion `fid` file into an m/z, the instrument description and the
//! acquisition method. The source documents the class as "only needed for
//! XMassFile"; [`crate::format::xmass_file`] is its one consumer here too.
//!
//! See `docs/XMASS_FILE_SUPPORT.md` for the member-by-member mapping and the
//! evidence.
//!
//! # Line grammar
//!
//! Exactly the source loop (`AcqusHandler.cpp:45-61`):
//!
//! - a line is everything up to `\n`; a Windows `\r` stays on the line and is
//!   removed later by the value trim;
//! - a line shorter than 5 bytes (`##x=x` is the shortest parameter) or not
//!   starting with `##` is skipped;
//! - the line is split on `=`; only a line with **exactly one** `=` defines a
//!   parameter, so `##$PDELAY= (0..3)` is kept and a value containing a second
//!   `=` drops the whole line;
//! - the key is the text between `##` and `=`, untrimmed; the value is the text
//!   after `=`, trimmed of space, tab, CR and LF (`StringUtils::trim`);
//! - a repeated key keeps the last value.
//!
//! # Text encoding
//!
//! The source stores raw bytes. A line that is valid UTF-8 is kept as it is; a
//! line that is not is decoded as ISO-8859-1, one character per byte, which is
//! reversible and is what the Windows acquisition software writes. The source
//! has no notion of an encoding, so this is a native decision; it only changes
//! how a non-ASCII byte is *represented*, never which bytes a value holds.

use crate::data_structures::list::ListParse;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Largest `acqus` file read, in bytes.
///
/// Native ceiling: the source reads any length. A Flex series `acqus` file is
/// about 5 KiB; 16 MiB leaves three orders of magnitude of headroom while
/// keeping the parameter map bounded.
pub const MAX_ACQUS_BYTES: u64 = 16 * 1024 * 1024;

/// Parsed `acqus` parameters and the time-of-flight calibration they define.
///
/// The source class `Internal::AcqusHandler`. Construction reads the whole file
/// and converts the six calibration parameters; the value is immutable
/// afterwards.
#[derive(Clone, Debug, PartialEq)]
pub struct AcqusHandler {
    /// Every `##key= value` line, the source `params_`.
    params: BTreeMap<String, String>,
    /// `$DW`, the dwell time: time-of-flight increment per point.
    dw: f64,
    /// `$DELAY`, the time-of-flight offset of point 0, as the source's `Size`.
    ///
    /// The source converts with `(Size)StringUtils::toInt32(...)`, so a
    /// negative delay wraps to a value near 2^64; that value is kept.
    delay: u64,
    /// `$ML1`, first calibration constant.
    ml1: f64,
    /// `$ML2`, second calibration constant (a time offset).
    ml2: f64,
    /// `$ML3`, quadratic calibration constant; 0 selects the linear formula.
    ml3: f64,
    /// `$TD`, the number of points, as the source's `Size`.
    ///
    /// Converted like [`Self::delay`], so a negative `$TD` wraps. Consumers
    /// that allocate from it must bound it themselves, as
    /// [`XMassFile`](crate::format::xmass_file::XMassFile) does.
    td: u64,
}

impl AcqusHandler {
    /// Open an `acqus` file and import its parameters.
    ///
    /// The source constructor `AcqusHandler(const std::string& filename)`.
    ///
    /// # Errors
    ///
    /// - The I/O error when the file cannot be opened or read. The source
    ///   throws `Exception::FileNotFound`, `FileNotReadable` or `IOException`
    ///   for the same three conditions; here they are the
    ///   [`std::io::ErrorKind`] of [`Error::Io`] (`NotFound`,
    ///   `PermissionDenied`, other).
    /// - [`Error::InvalidValue`] when the file exceeds [`MAX_ACQUS_BYTES`], a
    ///   native ceiling.
    /// - [`Error::Parse`] when one of `$DW`, `$DELAY`, `$ML1`, `$ML2`, `$ML3`
    ///   or `$TD` is missing or does not convert, the source's
    ///   `Exception::ConversionError`. A missing parameter converts the empty
    ///   string in the source, which fails the same way.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let file = std::fs::File::open(path.as_ref())?;
        Self::from_reader(BufReader::new(file))
    }

    /// Import parameters from an already opened `acqus` stream.
    ///
    /// Native entry point with the same parsing and errors as
    /// [`AcqusHandler::open`], minus the file-opening ones.
    ///
    /// # Errors
    ///
    /// As [`AcqusHandler::open`].
    pub fn from_reader(reader: impl BufRead) -> Result<Self> {
        let mut limited = reader.take(MAX_ACQUS_BYTES + 1);
        let mut params = BTreeMap::new();
        let mut line = Vec::new();
        let mut total: u64 = 0;
        loop {
            line.clear();
            let read = limited.read_until(b'\n', &mut line)?;
            if read == 0 {
                break;
            }
            total += read as u64;
            if total > MAX_ACQUS_BYTES {
                return Err(Error::InvalidValue(format!(
                    "acqus file exceeds {MAX_ACQUS_BYTES} bytes"
                )));
            }
            if line.last() == Some(&b'\n') {
                line.pop();
            }
            parse_line(&line, &mut params);
        }
        let dw = double_param(&params, "$DW")?;
        let delay = size_param(&params, "$DELAY")?;
        let ml1 = double_param(&params, "$ML1")?;
        let ml2 = double_param(&params, "$ML2")?;
        let ml3 = double_param(&params, "$ML3")?;
        let td = size_param(&params, "$TD")?;
        Ok(Self {
            params,
            dw,
            delay,
            ml1,
            ml2,
            ml3,
            td,
        })
    }

    /// The m/z of point `index`, from the internal time-of-flight calibration.
    ///
    /// The source `getPosition(Size index)`, after Titulaer *et al.*, BMC
    /// Bioinformatics 7:403 (2006):
    ///
    /// ```text
    /// tof = DW * index + DELAY
    /// a = ML3,  b = sqrt(1e12 / ML1),  c = ML2 - tof
    /// sqrt(m/z) = c / b                             if ML3 == 0
    /// sqrt(m/z) = (sqrt(b^2 - 4 a c) - b) / (2 a)   otherwise
    /// m/z = sqrt(m/z)^2
    /// ```
    ///
    /// Every operation is the source's, in the source's order, in `f64`, so the
    /// result is bit-identical to the Linux x86_64 Release build (measured in
    /// `tests/xmass_file.rs`). Nothing is checked, as in the source: a
    /// calibration that makes the discriminant negative yields NaN.
    pub fn position(&self, index: u64) -> f64 {
        let tof = self.dw * index as f64 + self.delay as f64;
        let a = self.ml3;
        let b = (1_000_000_000_000.0 / self.ml1).sqrt();
        let c = self.ml2 - tof;
        let sqrt_mz = if self.ml3 == 0.0 {
            c / b
        } else {
            ((b * b - 4.0 * a * c).sqrt() - b) / (2.0 * a)
        };
        sqrt_mz * sqrt_mz
    }

    /// A parameter's value as text; the empty string when the file has none.
    ///
    /// The source `getParam(const std::string&)`. Two names are computed rather
    /// than read, exactly as in the source: `mzMax` is the m/z of the last point
    /// (`position(TD - 1)`) and `mzMin` that of point 0, both written by
    /// `StringUtils::toStr(double)`
    /// ([`to_str`](crate::format::file_info::text_format::to_str)). With
    /// `$TD` 0 the source computes `TD - 1` in `Size` arithmetic, which wraps
    /// to 2^64 - 1; so does this.
    ///
    /// Keys are matched exactly: the `##` prefix is not part of a key, and a
    /// `$` or `.` that follows it is (`$TD`, `.IONIZATION MODE`).
    ///
    /// The source member is not `const` because `std::map::operator[]` inserts
    /// an empty entry for an unknown key; that insertion is unobservable, so
    /// this takes `&self`.
    pub fn param(&self, name: &str) -> String {
        match name {
            "mzMax" => crate::format::file_info::text_format::to_str(
                self.position(self.td.wrapping_sub(1)),
            ),
            "mzMin" => crate::format::file_info::text_format::to_str(self.position(0)),
            _ => self.params.get(name).cloned().unwrap_or_default(),
        }
    }

    /// The number of points of the spectrum, `$TD`.
    ///
    /// The source `getSize()`. A negative `$TD` wraps as the source's `Size`
    /// cast does, so `-1` reads as `u64::MAX`; a consumer that allocates from
    /// this value must bound it.
    pub fn size(&self) -> u64 {
        self.td
    }

    /// Every parameter the file defines, keyed without the `##` prefix.
    ///
    /// Native accessor; the source keeps `params_` private.
    pub fn params(&self) -> &BTreeMap<String, String> {
        &self.params
    }
}

/// One line of the source loop (`AcqusHandler.cpp:47-61`).
fn parse_line(line: &[u8], params: &mut BTreeMap<String, String>) {
    if line.len() < 5 || !line.starts_with(b"##") {
        return;
    }
    let mut parts = line.splitn(3, |byte| *byte == b'=');
    let (Some(key), Some(value), None) = (parts.next(), parts.next(), parts.next()) else {
        return;
    };
    params.insert(decode(&key[2..]), decode(trim(value)));
}

/// `StringUtils::trim`: space, tab, LF and CR on both ends.
fn trim(bytes: &[u8]) -> &[u8] {
    let is_space = |byte: &u8| matches!(byte, b' ' | b'\t' | b'\n' | b'\r');
    let start = bytes.iter().position(|b| !is_space(b)).unwrap_or(bytes.len());
    let end = bytes
        .iter()
        .rposition(|b| !is_space(b))
        .map_or(start, |last| last + 1);
    &bytes[start..end]
}

/// UTF-8 when valid, otherwise ISO-8859-1 (one character per byte).
fn decode(bytes: &[u8]) -> String {
    match std::str::from_utf8(bytes) {
        Ok(text) => text.to_owned(),
        Err(_) => bytes.iter().map(|&byte| char::from(byte)).collect(),
    }
}

/// `StringUtils::toDouble(params_[name])`.
fn double_param(params: &BTreeMap<String, String>, name: &str) -> Result<f64> {
    let text = params.get(name).map(String::as_str).unwrap_or("");
    <f64 as ListParse>::from_list_item(text).map_err(|_| Error::Parse {
        line: 0,
        message: format!("acqus parameter {name}: could not convert string '{text}' to a double value"),
    })
}

/// `(Size)StringUtils::toInt32(params_[name])`: a negative value wraps.
fn size_param(params: &BTreeMap<String, String>, name: &str) -> Result<u64> {
    let text = params.get(name).map(String::as_str).unwrap_or("");
    let value = <i32 as ListParse>::from_list_item(text).map_err(|_| Error::Parse {
        line: 0,
        message: format!("acqus parameter {name}: could not convert string '{text}' to an integer value"),
    })?;
    Ok(i64::from(value) as u64)
}

#[cfg(test)]
mod tests {
    use super::*;

    const MINIMAL: &str = "##$DW= 0.5 \r\n##$DELAY= 10\r\n##$ML1= 2602987.3759832\r\n##$ML2= 151.937451822094\r\n##$ML3= 0\r\n##$TD= 4\r\n";

    #[test]
    fn line_grammar_follows_the_source_loop() {
        let text = format!(
            "{MINIMAL}##a=b\r\n#x=yy\r\n##ab=\r\n##k=v=w\r\n##$KEY=  spaced value \t\r\n##$KEY= last\r\n##=empty key\r\n"
        );
        let handler = AcqusHandler::from_reader(text.as_bytes()).unwrap();
        // `##a=b` is exactly 5 bytes before the CR and is kept.
        assert_eq!(handler.param("a"), "b");
        // Not `##`: skipped.
        assert!(!handler.params().contains_key("x"));
        // An empty value is still a parameter.
        assert_eq!(handler.params().get("ab").map(String::as_str), Some(""));
        // Two `=`: the whole line is dropped.
        assert!(!handler.params().contains_key("k"));
        // Last value wins; values are trimmed, keys are not.
        assert_eq!(handler.param("$KEY"), "last");
        assert_eq!(handler.param(""), "empty key");
        assert_eq!(handler.param("unknown"), "");
        assert_eq!(handler.size(), 4);
    }

    #[test]
    fn short_lines_are_skipped() {
        let text = format!("{MINIMAL}##a=\n");
        let handler = AcqusHandler::from_reader(text.as_bytes()).unwrap();
        assert!(!handler.params().contains_key("a"));
    }

    #[test]
    fn missing_calibration_is_a_conversion_error() {
        let text = MINIMAL.replace("##$ML3= 0\r\n", "");
        assert!(matches!(
            AcqusHandler::from_reader(text.as_bytes()),
            Err(Error::Parse { .. })
        ));
        let text = MINIMAL.replace("##$TD= 4", "##$TD= 4.5");
        assert!(matches!(
            AcqusHandler::from_reader(text.as_bytes()),
            Err(Error::Parse { .. })
        ));
    }

    #[test]
    fn negative_sizes_wrap_like_the_source_cast() {
        let text = MINIMAL.replace("##$TD= 4", "##$TD= -1").replace("##$DELAY= 10", "##$DELAY= -2");
        let handler = AcqusHandler::from_reader(text.as_bytes()).unwrap();
        assert_eq!(handler.size(), u64::MAX);
        assert_eq!(handler.delay, u64::MAX - 1);
    }

    #[test]
    fn invalid_utf8_is_decoded_as_latin1() {
        let mut bytes = MINIMAL.as_bytes().to_vec();
        bytes.extend_from_slice(b"##$UNIT= \xb5s\r\n");
        let handler = AcqusHandler::from_reader(bytes.as_slice()).unwrap();
        assert_eq!(handler.param("$UNIT"), "\u{b5}s");
    }

    #[test]
    fn oversized_file_is_refused() {
        let big = vec![b'#'; (MAX_ACQUS_BYTES + 1) as usize];
        assert!(matches!(
            AcqusHandler::from_reader(big.as_slice()),
            Err(Error::InvalidValue(_))
        ));
    }
}
