// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Universal Spectrum Identifiers (USI).
//!
//! Ports `src/openms/include/OpenMS/METADATA/USI.h` and
//! `src/openms/source/METADATA/USI.cpp` (core SDK `bc9cc12`). See
//! `docs/USI_SUPPORT.md`.
//!
//! The Universal Spectrum Identifier is the HUPO-PSI format for referencing a
//! mass spectrum unambiguously across the ProteomeXchange repositories (PRIDE,
//! PeptideAtlas, MassIVE, jPOST, iProX) and spectral libraries
//! (<https://www.psidev.info/usi>). Its general form is
//!
//! ```text
//! mzspec:<collection>:<ms_run>:<index_type>:<index>[:interpretation]
//! ```
//!
//! - **collection**: a ProteomeXchange dataset identifier (`PXD000561`) or a
//!   spectral library name;
//! - **ms_run**: the name of the MS run file (`sample1.mzML`);
//! - **index_type**: `scan`, `index` or `nativeId`;
//! - **index**: the spectrum's identifier of that type;
//! - **interpretation**: an optional ProForma peptide interpretation
//!   (`sequence/charge`), which may itself contain colons
//!   (`EM[UNIMOD:35]K/2`).
//!
//! The PSI-MS term for a USI is MS:1003063, *universal spectrum identifier*.
//!
//! See also `SpectrumLookup`, `SpectrumSettings` and `PeptideIdentification`
//! in the source.
//!
//! # Example
//!
//! ```
//! use openms::metadata::usi::{IndexType, Usi};
//!
//! let usi: Usi = "mzspec:PXD000561:sample.mzML:scan:1234:PEPT[Phospho]IDEK/2"
//!     .parse()
//!     .unwrap();
//! assert_eq!(usi.index_type, IndexType::Scan);
//! assert_eq!(usi.interpretation, "PEPT[Phospho]IDEK/2");
//! assert_eq!(
//!     usi.to_string(),
//!     "mzspec:PXD000561:sample.mzML:scan:1234:PEPT[Phospho]IDEK/2"
//! );
//!
//! // Prefer `try_parse` over validating and then parsing again.
//! assert!(Usi::try_parse("mzspec:PXD000561").is_none());
//! ```

use std::fmt;
use std::str::FromStr;

use crate::concept::boost_regex::BoostRegex;
use crate::metadata::spectrum_native_id_parser::{
    is_native_id, regex_from_native_id, try_extract_scan_number_with,
};
use crate::{Error, Result};

/// The prefix every USI starts with (`USI::USI_PREFIX`).
pub const USI_PREFIX: &str = "mzspec:";

/// The PSI-MS CV accession for a USI (`USI::getCVAccession`).
pub const CV_ACCESSION: &str = "MS:1003063";

/// The PSI-MS CV name for a USI (`USI::getCVName`).
pub const CV_NAME: &str = "universal spectrum identifier";

/// The kind of spectrum index a USI carries (`USI::IndexType`).
///
/// Declaration order is the source's, which [`Usi`]'s ordering uses.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum IndexType {
    /// Scan number, the most common index type (`SCAN`). The default, as the
    /// source's default constructor sets it.
    #[default]
    Scan,
    /// Zero-based spectrum index (`INDEX`).
    Index,
    /// Native spectrum identifier from the instrument (`NATIVEID`).
    NativeId,
}

impl IndexType {
    /// The index type's USI spelling: `scan`, `index` or `nativeId`
    /// (`USI::indexTypeToString`).
    ///
    /// The source's `default: return "scan"` arm cannot be reached from a Rust
    /// enum and has no counterpart.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Scan => "scan",
            Self::Index => "index",
            Self::NativeId => "nativeId",
        }
    }

    /// Parse an index type, case-insensitively (`USI::indexTypeFromString`).
    ///
    /// `scan`, `index` and `nativeid` are accepted in any ASCII case.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for any other text, where the source throws
    /// `Exception::InvalidValue` ("Invalid USI index type").
    pub fn parse(type_string: &str) -> Result<Self> {
        match type_string.to_ascii_lowercase().as_str() {
            "scan" => Ok(Self::Scan),
            "index" => Ok(Self::Index),
            "nativeid" => Ok(Self::NativeId),
            _ => Err(Error::InvalidValue(format!(
                "Invalid USI index type: '{type_string}'"
            ))),
        }
    }
}

impl fmt::Display for IndexType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl FromStr for IndexType {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        Self::parse(text)
    }
}

/// A Universal Spectrum Identifier (`USI`).
///
/// The five components are public fields: the source's `get`/`set` pairs
/// (`getCollection`/`setCollection`, `getMSRun`/`setMSRun`,
/// `getIndexType`/`setIndexType`, `getIndex`/`setIndex`,
/// `getInterpretation`/`setInterpretation`) assign them without validation, so
/// a field is the same contract. [`Default`] is the source's default
/// constructor: every string empty and [`IndexType::Scan`], which is not
/// [valid](Usi::is_valid).
///
/// Equality compares all five fields (`operator==`/`operator!=`). The derived
/// order compares collection, MS run, index type (in declaration order), index
/// and interpretation, each strings byte-wise, which is the source's
/// `operator<` field by field.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Usi {
    /// ProteomeXchange dataset identifier or spectral library name.
    pub collection: String,
    /// MS run file name.
    pub ms_run: String,
    /// Type of the spectrum index.
    pub index_type: IndexType,
    /// Spectrum index value, as text.
    pub index: String,
    /// Optional peptide interpretation (ProForma `sequence/charge`); empty when
    /// absent.
    pub interpretation: String,
}

impl Usi {
    /// Construct a USI from its components (the source's five-argument
    /// constructor, whose `interpretation` defaults to empty).
    ///
    /// Nothing is validated, as in the source; see [`Usi::is_valid`].
    pub fn new(
        collection: impl Into<String>,
        ms_run: impl Into<String>,
        index_type: IndexType,
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

    /// Parse a USI string (the source's `explicit USI(const std::string&)`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when [`Usi::set_from_string`] rejects the text,
    /// where the source throws `Exception::ParseError` ("Invalid USI format").
    pub fn parse(usi_string: &str) -> Result<Self> {
        let mut usi = Self::default();
        if usi.set_from_string(usi_string) {
            Ok(usi)
        } else {
            Err(Error::InvalidValue(format!(
                "{usi_string}: Invalid USI format"
            )))
        }
    }

    /// Whether the USI has every required component: a non-empty collection,
    /// MS run and index (`isValid`).
    pub fn is_valid(&self) -> bool {
        !self.collection.is_empty() && !self.ms_run.is_empty() && !self.index.is_empty()
    }

    /// Whether `usi_string` is a valid USI (`isValidUSI`), which is
    /// [`Usi::try_parse`] succeeding.
    pub fn is_valid_usi(usi_string: &str) -> bool {
        Self::try_parse(usi_string).is_some()
    }

    /// Parse `usi_string`, or `None` when it is not a valid USI (`tryParse`).
    ///
    /// This is the preferred way to validate and use a USI in one step,
    /// avoiding the double parse of [`Usi::is_valid_usi`] followed by
    /// [`Usi::parse`]:
    ///
    /// ```
    /// use openms::metadata::usi::Usi;
    ///
    /// let text = "mzspec:PXD000561:sample.mzML:scan:12345";
    /// if let Some(usi) = Usi::try_parse(text) {
    ///     assert_eq!(usi.index, "12345");
    /// }
    /// ```
    pub fn try_parse(usi_string: &str) -> Option<Self> {
        if !usi_string.starts_with(USI_PREFIX) {
            return None;
        }
        let mut usi = Self::default();
        usi.set_from_string(usi_string).then_some(usi)
    }

    /// Whether an interpretation is present, i.e. non-empty
    /// (`hasInterpretation`).
    pub fn has_interpretation(&self) -> bool {
        !self.interpretation.is_empty()
    }

    /// The USI's string form, or the empty string when it is not
    /// [valid](Usi::is_valid) (`toString`).
    ///
    /// `mzspec:<collection>:<ms_run>:<index_type>:<index>`, followed by
    /// `:<interpretation>` when one is present. [`fmt::Display`] writes the same
    /// text, as the source's `operator<<` does.
    pub fn to_usi_string(&self) -> String {
        if !self.is_valid() {
            return String::new();
        }
        let mut out = format!(
            "{USI_PREFIX}{}:{}:{}:{}",
            self.collection,
            self.ms_run,
            self.index_type.as_str(),
            self.index
        );
        if !self.interpretation.is_empty() {
            out.push(':');
            out.push_str(&self.interpretation);
        }
        out
    }

    /// Parse `usi_string` into `self`, returning whether the result is valid
    /// (`fromString`).
    ///
    /// Every field is reset first. The text must start with `mzspec:`; the rest
    /// is split on every `:`, keeping empty pieces, and needs at least four
    /// pieces: collection, MS run, index type ([`IndexType::parse`],
    /// case-insensitive) and index. Any further pieces are joined back with `:`
    /// as the interpretation, so a ProForma interpretation keeps its colons.
    ///
    /// # Notes
    ///
    /// Like the source, a rejected string can leave `self` partly filled: the
    /// fields are assigned as they are read, and a missing prefix, too few
    /// pieces or an unknown index type stop there. The return value, not the
    /// state, is the answer. A trailing `:` gives an empty interpretation, so
    /// [`Usi::to_usi_string`] does not reproduce it.
    pub fn set_from_string(&mut self, usi_string: &str) -> bool {
        *self = Self::default();
        let Some(remainder) = usi_string.strip_prefix(USI_PREFIX) else {
            return false;
        };
        let parts: Vec<&str> = remainder.splitn(5, ':').collect();
        if parts.len() < 4 {
            return false;
        }
        self.collection = parts[0].to_string();
        self.ms_run = parts[1].to_string();
        match IndexType::parse(parts[2]) {
            Ok(index_type) => self.index_type = index_type,
            Err(_) => return false,
        }
        self.index = parts[3].to_string();
        if let Some(interpretation) = parts.get(4) {
            self.interpretation = (*interpretation).to_string();
        }
        self.is_valid()
    }

    /// A USI referencing `scan_number` in `filename` of `dataset_id`
    /// (`createFromScanNumber`), with [`IndexType::Scan`] and the decimal scan
    /// number as index.
    ///
    /// `filename` is used as given; a full path stays a full path. Use
    /// [`Usi::extract_basename`] first to reduce it to the run name.
    pub fn from_scan_number(
        dataset_id: &str,
        filename: &str,
        scan_number: i32,
        interpretation: &str,
    ) -> Self {
        Self::new(
            dataset_id,
            filename,
            IndexType::Scan,
            scan_number.to_string(),
            interpretation,
        )
    }

    /// A USI referencing the spectrum with `native_id` in `filename` of
    /// `dataset_id` (`createFromNativeID`).
    ///
    /// When [`Usi::extract_scan_number_from_native_id`] finds a scan number,
    /// the USI uses [`IndexType::Scan`] with that number; otherwise it falls
    /// back to [`IndexType::NativeId`] with the native ID itself as index. No
    /// interpretation is set.
    ///
    /// # Errors
    ///
    /// Only [`Error::InvalidValue`] when a regular-expression search exceeds
    /// its work budget; see [`Usi::extract_scan_number_from_native_id`].
    pub fn from_native_id(dataset_id: &str, filename: &str, native_id: &str) -> Result<Self> {
        Ok(match Self::extract_scan_number_from_native_id(native_id)? {
            Some(scan) => Self::new(dataset_id, filename, IndexType::Scan, scan.to_string(), ""),
            None => Self::new(dataset_id, filename, IndexType::NativeId, native_id, ""),
        })
    }

    /// The scan number encoded in `native_id`, if any
    /// (`extractScanNumberFromNativeID`).
    ///
    /// When `native_id` is a vendor native ID
    /// ([`is_native_id`]), the expression
    /// [`regex_from_native_id`] chooses for it is
    /// applied and its last match's number returned; `index=100` therefore
    /// gives `100`, not the `101` that the accession-driven
    /// [`extract_scan_number`](crate::metadata::spectrum_native_id_parser::extract_scan_number) would return for
    /// MS:1000774. Failing that — or for a string that is not a native ID —
    /// the whole string is converted as `StringUtils::toInt32` does, so `12345`
    /// gives `12345`.
    ///
    /// # Notes
    ///
    /// A WIFF ID (`sample=... cycle=... experiment=...`) is a native ID whose
    /// prefix has no expression of its own, so the generic `(?<GROUP>\d+)` is
    /// used and its last match — the experiment number — is returned. That is
    /// the source's behaviour.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a regular-expression search exceeds its
    /// work budget (see `BoostRegex::tokens`); the source has no such limit.
    pub fn extract_scan_number_from_native_id(native_id: &str) -> Result<Option<i32>> {
        if is_native_id(native_id) {
            let regex = BoostRegex::new(regex_from_native_id(native_id))?;
            if let Some(scan) = try_extract_scan_number_with(native_id, &regex)? {
                if scan >= 0 {
                    return Ok(Some(scan));
                }
            }
        }
        Ok(crate::concept::string_conversion::to_int32(native_id).ok())
    }

    /// The file name of `filepath` without its directory
    /// (`extractBasename`).
    ///
    /// A `file://` prefix is removed first, and with it the leading `/` of a
    /// Windows drive path (`file:///C:/data/sample.mzML`). The basename is then
    /// everything after the last `/` or `\`, whichever comes later, so Unix,
    /// Windows and mixed separators all work; a string without a separator is
    /// returned unchanged, and an empty string stays empty.
    pub fn extract_basename(filepath: &str) -> String {
        let mut path = filepath;
        if let Some(rest) = path.strip_prefix("file://") {
            path = rest;
            let bytes = path.as_bytes();
            if bytes.len() > 2 && bytes[0] == b'/' && bytes[2] == b':' {
                path = &path[1..];
            }
        }
        match path.rfind(['/', '\\']) {
            Some(separator) => path[separator + 1..].to_string(),
            None => path.to_string(),
        }
    }
}

impl fmt::Display for Usi {
    /// The source's `operator<<`: [`Usi::to_usi_string`], empty for an invalid
    /// USI.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.to_usi_string())
    }
}

impl FromStr for Usi {
    type Err = Error;

    /// [`Usi::parse`].
    fn from_str(text: &str) -> Result<Self> {
        Self::parse(text)
    }
}
