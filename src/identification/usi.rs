// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Timo Sachsenberg, OpenMS Rust contributors $

//! Universal Spectrum Identifiers (USI), from `METADATA/USI.h`.
//!
//! The Universal Spectrum Identifier is the HUPO-PSI format for referencing a
//! mass spectrum across public repositories (PRIDE, PeptideAtlas, MassIVE,
//! jPOST, iProX) and spectral libraries. Its general form is
//!
//! ```text
//! mzspec:<collection>:<ms_run>:<index_type>:<index>[:interpretation]
//! ```
//!
//! where *collection* is a ProteomeXchange dataset identifier (for example
//! `PXD000561`) or a spectral library name, *ms_run* the name of the run file,
//! *index_type* `scan`, `index` or `nativeId`, *index* the spectrum's
//! identifier of that type, and the optional *interpretation* a ProForma
//! peptide with charge (`PEPTIDEK/2`, `PEPT[Phospho]IDEK/2`). USI is the
//! PSI-MS term `MS:1003063`. Specification: <https://www.psidev.info/usi>.
//!
//! Ported at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`; every member
//! was checked against the executed Linux x86_64 Release build (oracle
//! `B8b-peptide-id`, cases `usi_class`, `usi_build` and `usi_mapping`). See
//! `docs/USI_SUPPORT.md`.

use crate::{Error, Result};
use std::fmt;
use std::str::FromStr;

/// Kind of spectrum index a USI carries (source `USI::IndexType`).
///
/// The declaration order is the source enum order, which `Usi`'s ordering
/// compares: `Scan < Index < NativeId`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum UsiIndexType {
    /// Scan number, the most common kind (source `SCAN`). The default.
    #[default]
    Scan,
    /// Zero-based spectrum index (source `INDEX`).
    Index,
    /// Native spectrum identifier from the instrument (source `NATIVEID`).
    NativeId,
}

impl UsiIndexType {
    /// Source `USI::indexTypeToString`: `"scan"`, `"index"` or `"nativeId"`.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Scan => "scan",
            Self::Index => "index",
            Self::NativeId => "nativeId",
        }
    }

    /// Source `USI::indexTypeFromString`: the three names, compared after ASCII
    /// lower-casing, so `SCAN` and `NativeID` are accepted.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for any other text, as the source throws
    /// `Exception::InvalidValue` ("Invalid USI index type").
    pub fn parse(type_string: &str) -> Result<Self> {
        match type_string.to_ascii_lowercase().as_str() {
            "scan" => Ok(Self::Scan),
            "index" => Ok(Self::Index),
            "nativeid" => Ok(Self::NativeId),
            _ => Err(Error::InvalidValue(format!(
                "the value '{type_string}' was used but is not valid; Invalid USI index type"
            ))),
        }
    }
}

impl fmt::Display for UsiIndexType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl FromStr for UsiIndexType {
    type Err = Error;
    fn from_str(text: &str) -> Result<Self> {
        Self::parse(text)
    }
}

/// A Universal Spectrum Identifier (source `OpenMS::USI`).
///
/// The five components are public fields; the source's getters and setters
/// (`getCollection`/`setCollection`, `getMSRun`/`setMSRun`,
/// `getIndexType`/`setIndexType`, `getIndex`/`setIndex`,
/// `getInterpretation`/`setInterpretation`) are plain field access. The
/// default value is the source default constructor's empty, invalid USI with
/// index type `scan`.
///
/// Equality compares all five components, as the source `operator==` does.
/// The derived ordering is the source `operator<`: collection, then MS run,
/// then index type in enum order, then index, then interpretation, each string
/// compared bytewise (so index `"10"` sorts before `"2"`).
///
/// ```
/// use openms::identification::usi::{Usi, UsiIndexType};
///
/// let usi = Usi::parse("mzspec:PXD000561:sample.mzML:scan:1234:PEPT[Phospho]IDEK/2")?;
/// assert_eq!(usi.index_type, UsiIndexType::Scan);
/// assert_eq!(usi.interpretation, "PEPT[Phospho]IDEK/2");
/// assert_eq!(usi.to_string(), "mzspec:PXD000561:sample.mzML:scan:1234:PEPT[Phospho]IDEK/2");
/// # Ok::<(), openms::Error>(())
/// ```
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Usi {
    /// ProteomeXchange dataset identifier or spectral library name.
    pub collection: String,
    /// Name of the MS run file.
    pub ms_run: String,
    /// Kind of index in [`Usi::index`].
    pub index_type: UsiIndexType,
    /// The spectrum's identifier, as text.
    pub index: String,
    /// Optional ProForma peptide interpretation with charge; empty when absent.
    pub interpretation: String,
}

impl Usi {
    /// Prefix of every USI (source protected `USI_PREFIX`).
    pub const PREFIX: &'static str = "mzspec:";
    /// PSI-MS accession of the USI term (source `getCVAccession`).
    pub const CV_ACCESSION: &'static str = "MS:1003063";
    /// PSI-MS name of the USI term (source `getCVName`).
    pub const CV_NAME: &'static str = "universal spectrum identifier";
    /// Collection the source `buildUSI` overloads default to for unpublished
    /// data (their `dataset_id = "local"` default argument).
    pub const LOCAL_COLLECTION: &'static str = "local";

    /// Source component constructor
    /// `USI(collection, ms_run, index_type, index, interpretation = "")`; pass
    /// `""` for the source default interpretation. Nothing is validated, as in
    /// the source: the result may be invalid (see [`Usi::is_valid`]).
    pub fn new(
        collection: impl Into<String>,
        ms_run: impl Into<String>,
        index_type: UsiIndexType,
        index: impl Into<String>,
        interpretation: impl Into<String>,
    ) -> Self {
        Self {
            collection: collection.into(),
            ms_run: ms_run.into(),
            index_type,
            index: index.into(),
            interpretation: interpretation.into(),
        }
    }

    /// Source `explicit USI(const std::string& usi_string)`: parse a complete
    /// USI string.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] with the source message `Invalid USI format in: <text>`
    /// when [`Usi::from_string`] rejects the text (the source throws
    /// `Exception::ParseError`).
    pub fn parse(usi_string: &str) -> Result<Self> {
        let mut usi = Self::default();
        if usi.from_string(usi_string) {
            Ok(usi)
        } else {
            Err(Error::Parse {
                line: 0,
                message: format!("Invalid USI format in: {usi_string}"),
            })
        }
    }

    /// Source `isValid`: collection, MS run and index are all non-empty.
    pub fn is_valid(&self) -> bool {
        !self.collection.is_empty() && !self.ms_run.is_empty() && !self.index.is_empty()
    }

    /// Source static `isValidUSI`: whether [`Usi::try_parse`] accepts the text.
    pub fn is_valid_usi(usi_string: &str) -> bool {
        Self::try_parse(usi_string).is_some()
    }

    /// Source static `tryParse`: the parsed USI, or `None` when the text does
    /// not start with `mzspec:` or [`Usi::from_string`] rejects it. Prefer this
    /// to [`Usi::is_valid_usi`] followed by [`Usi::parse`], which parses twice.
    pub fn try_parse(usi_string: &str) -> Option<Self> {
        if !usi_string.starts_with(Self::PREFIX) {
            return None;
        }
        let mut usi = Self::default();
        usi.from_string(usi_string).then_some(usi)
    }

    /// Source `hasInterpretation`: the interpretation is non-empty.
    pub fn has_interpretation(&self) -> bool {
        !self.interpretation.is_empty()
    }

    /// Source `fromString`: parse `usi_string` into this value and report
    /// whether the result is a valid USI.
    ///
    /// All components are reset first. The text must start with `mzspec:`
    /// (case-sensitive) and split on `:` into at least four parts: collection,
    /// MS run, index type ([`UsiIndexType::parse`], case-insensitive) and
    /// index. Parts after the fourth are joined back with `:` into the
    /// interpretation, so an interpretation may contain colons and a trailing
    /// `:` leaves it empty. The return value is [`Usi::is_valid`].
    ///
    /// The source assigns as it goes, so a rejected string can leave partial
    /// state, and this port keeps that observable behaviour: with an unknown
    /// index type the collection and MS run are already set while index and
    /// interpretation stay empty; with fewer than four parts everything stays
    /// reset. [`Usi::parse`] and [`Usi::try_parse`] never expose that state.
    pub fn from_string(&mut self, usi_string: &str) -> bool {
        *self = Self::default();
        let Some(remainder) = usi_string.strip_prefix(Self::PREFIX) else {
            return false;
        };
        let parts: Vec<&str> = remainder.split(':').collect();
        if parts.len() < 4 {
            return false;
        }
        self.collection = parts[0].to_owned();
        self.ms_run = parts[1].to_owned();
        match UsiIndexType::parse(parts[2]) {
            Ok(index_type) => self.index_type = index_type,
            Err(_) => return false,
        }
        self.index = parts[3].to_owned();
        if parts.len() > 4 {
            self.interpretation = parts[4..].join(":");
        }
        self.is_valid()
    }

    /// Source static `createFromScanNumber`: a `scan` USI whose index is the
    /// decimal scan number (a negative number is written with its sign, as the
    /// source does). Pass `""` for the source default interpretation.
    pub fn from_scan_number(
        dataset_id: &str,
        filename: &str,
        scan_number: i32,
        interpretation: &str,
    ) -> Self {
        Self::new(
            dataset_id,
            filename,
            UsiIndexType::Scan,
            scan_number.to_string(),
            interpretation,
        )
    }

    /// Source static `createFromNativeID`: a `scan` USI when
    /// [`Usi::extract_scan_number_from_native_id`] finds a scan number, and a
    /// `nativeId` USI carrying the whole native ID otherwise.
    pub fn from_native_id(dataset_id: &str, filename: &str, native_id: &str) -> Self {
        match Self::extract_scan_number_from_native_id(native_id) {
            Some(scan) => Self::new(
                dataset_id,
                filename,
                UsiIndexType::Scan,
                scan.to_string(),
                "",
            ),
            None => Self::new(dataset_id, filename, UsiIndexType::NativeId, native_id, ""),
        }
    }

    /// Source static `extractScanNumberFromNativeID`: the scan number a native
    /// spectrum identifier carries, if one can be read.
    ///
    /// Two steps, as in the source:
    ///
    /// 1. When the identifier starts with a known native-ID prefix (`scan=`,
    ///    `scanId=`, `scanID=`, `controllerType=`, `function=`, `sample=`,
    ///    `index=`, `spectrum=`, `file=` or `frame=`, source
    ///    `SpectrumNativeIDParser::isNativeID`), the pattern
    ///    `SpectrumNativeIDParser::getRegExFromNativeID` chooses for that prefix
    ///    is matched, and the digits of the **last** match are converted. The
    ///    patterns are `scan=(\d+)` (for `scan=`, `controllerType=`, `function=`
    ///    and `frame=`), `index=(\d+)`, `scanId=(\d+)`, `scanID=(\d+)`,
    ///    `spectrum=(\d+)`, `file=(\d+)`, and a bare `(\d+)` for everything
    ///    else, including the WIFF prefix `sample=`: for
    ///    `sample=1 period=1 cycle=42 experiment=3` the scan number is 3.
    ///    `index=7` gives 7, not the index+1 of the accession-driven parser.
    /// 2. Otherwise, or when step 1 finds nothing convertible, the whole text
    ///    is converted as `StringUtils::toInt32` does: ASCII blanks around it,
    ///    one optional `+`, an optionally negative decimal that fits `i32`.
    ///
    /// Digit runs that overflow `i32` are not a scan number (`scan=2147483648`
    /// gives `None`). The patterns are literal prefixes followed by ASCII
    /// digits, so they are matched directly rather than through a regular
    /// expression engine; the result is the same as Boost's non-overlapping
    /// leftmost search, which the test suite checks against
    /// [`BoostRegex`](crate::concept::boost_regex::BoostRegex) on every
    /// oracle input.
    pub fn extract_scan_number_from_native_id(native_id: &str) -> Option<i32> {
        if is_native_id(native_id) {
            let prefix = native_id_scan_prefix(native_id);
            if let Some(digits) = last_prefixed_digits(native_id, prefix) {
                if let Some(scan) = super::source_int32(digits) {
                    return Some(scan);
                }
            }
        }
        super::source_int32(native_id)
    }

    /// Source static `extractBasename`: the last path component of a file
    /// path or `file://` URI.
    ///
    /// A `file://` prefix is removed, and then one more leading `/` when the
    /// third character is `:` (`file:///C:/data/x.mzML` becomes
    /// `C:/data/x.mzML`). The text after the last `/` or `\`, whichever comes
    /// later, is returned; a path ending in a separator yields `""`, and text
    /// without a separator is returned whole.
    pub fn extract_basename(filepath: &str) -> &str {
        let mut path = filepath;
        if let Some(rest) = path.strip_prefix("file://") {
            path = rest;
            let bytes = path.as_bytes();
            if bytes.len() > 2 && bytes[0] == b'/' && bytes[2] == b':' {
                path = &path[1..];
            }
        }
        match path.rfind(&['/', '\\'][..]) {
            Some(separator) => &path[separator + 1..],
            None => path,
        }
    }
}

impl fmt::Display for Usi {
    /// Source `toString` and `operator<<`: the complete USI text, or the empty
    /// string when [`Usi::is_valid`] is false. The interpretation and its
    /// leading `:` are written only when present.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if !self.is_valid() {
            return Ok(());
        }
        write!(
            f,
            "{}{}:{}:{}:{}",
            Self::PREFIX,
            self.collection,
            self.ms_run,
            self.index_type,
            self.index
        )?;
        if !self.interpretation.is_empty() {
            write!(f, ":{}", self.interpretation)?;
        }
        Ok(())
    }
}

impl FromStr for Usi {
    type Err = Error;
    /// Same as [`Usi::parse`].
    fn from_str(text: &str) -> Result<Self> {
        Self::parse(text)
    }
}

/// Source `SpectrumNativeIDParser::isNativeID`: one of the ten native-ID
/// prefixes starts the identifier.
fn is_native_id(id: &str) -> bool {
    [
        "scan=",
        "scanId=",
        "scanID=",
        "controllerType=",
        "function=",
        "sample=",
        "index=",
        "spectrum=",
        "file=",
        "frame=",
    ]
    .iter()
    .any(|prefix| id.starts_with(prefix))
}

/// The literal part of the pattern source
/// `SpectrumNativeIDParser::getRegExFromNativeID` returns for `id`; the
/// pattern is this text followed by `(?<GROUP>\d+)`.
fn native_id_scan_prefix(id: &str) -> &'static str {
    if ["scan=", "controllerType=", "function=", "frame="]
        .iter()
        .any(|prefix| id.starts_with(prefix))
    {
        "scan="
    } else if id.starts_with("index=") {
        "index="
    } else if id.starts_with("scanId=") {
        "scanId="
    } else if id.starts_with("scanID=") {
        "scanID="
    } else if id.starts_with("spectrum=") {
        "spectrum="
    } else if id.starts_with("file=") {
        "file="
    } else {
        ""
    }
}

/// The digits of the last non-overlapping leftmost match of `prefix(\d+)` in
/// `text`, where `\d` is an ASCII digit as in Boost's C locale. `\d+` is
/// greedy, so every match takes the whole digit run after the prefix, and the
/// search resumes after it.
fn last_prefixed_digits<'a>(text: &'a str, prefix: &str) -> Option<&'a str> {
    let bytes = text.as_bytes();
    let mut last = None;
    let mut position = 0;
    while position < bytes.len() {
        let digits_from = position + prefix.len();
        if bytes[position..].starts_with(prefix.as_bytes())
            && bytes.get(digits_from).is_some_and(u8::is_ascii_digit)
        {
            let digits_to = bytes[digits_from..]
                .iter()
                .position(|byte| !byte.is_ascii_digit())
                .map_or(bytes.len(), |offset| digits_from + offset);
            // Both bounds sit next to ASCII bytes, so they are char boundaries.
            last = Some(&text[digits_from..digits_to]);
            position = digits_to;
        } else {
            position += 1;
        }
    }
    last
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prefixed_digits_take_greedy_runs_and_resume_after_each_match() {
        assert_eq!(last_prefixed_digits("scan=12 scan=34", "scan="), Some("34"));
        assert_eq!(last_prefixed_digits("scan=x scan=5", "scan="), Some("5"));
        assert_eq!(last_prefixed_digits("scan=", "scan="), None);
        assert_eq!(last_prefixed_digits("a12b345", ""), Some("345"));
        assert_eq!(last_prefixed_digits("abc", ""), None);
        assert_eq!(last_prefixed_digits("scan=scan=7", "scan="), Some("7"));
        assert_eq!(last_prefixed_digits("é scan=9", "scan="), Some("9"));
    }
}
