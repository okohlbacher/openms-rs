// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Look-up of spectrum meta data by spectrum reference, ported from the
//! regular-expression path of `METADATA/SpectrumLookup.h` and
//! `METADATA/SpectrumMetaDataLookup.h`.
//!
//! Search-engine result files often identify the searched spectrum only by a
//! *spectrum reference* — a string such as `sample_SII_12_1_1_2_1` or
//! `spectrum:...;rt:2187.96;mz:1093.0258;charge:2`. A reference format is a
//! Boost regular expression with named groups; the groups either carry the
//! meta data directly (`RT`, `PRECRT`, `MZ`, `CHARGE`, `LEVEL`, `SCAN`, `ID`) or
//! locate the spectrum among those read with
//! [`SpectrumMetaDataLookup::read_spectra`](crate::format::spectrum_meta_data_lookup::SpectrumMetaDataLookup::read_spectra) (`INDEX0`, `INDEX1`, `SCAN`, `ID`,
//! `RT`), whose recorded meta data then answers the query.
//!
//! | Meta data | Flag | Group | Meaning |
//! |---|---|---|---|
//! | `rt` | [`MetaDataFlags::RT`](crate::format::spectrum_meta_data_lookup::MetaDataFlags::RT) | `RT` | retention time of the spectrum |
//! | `precursor_rt` | [`MetaDataFlags::PRECURSOR_RT`](crate::format::spectrum_meta_data_lookup::MetaDataFlags::PRECURSOR_RT) | `PRECRT` | retention time of the precursor spectrum (undefined for MS1) |
//! | `precursor_mz` | [`MetaDataFlags::PRECURSOR_MZ`](crate::format::spectrum_meta_data_lookup::MetaDataFlags::PRECURSOR_MZ) | `MZ` | precursor m/z (undefined for MS1) |
//! | `precursor_charge` | [`MetaDataFlags::PRECURSOR_CHARGE`](crate::format::spectrum_meta_data_lookup::MetaDataFlags::PRECURSOR_CHARGE) | `CHARGE` | precursor charge (undefined for MS1) |
//! | `ms_level` | [`MetaDataFlags::MS_LEVEL`](crate::format::spectrum_meta_data_lookup::MetaDataFlags::MS_LEVEL) | `LEVEL` | MS level |
//! | `scan_number` | [`MetaDataFlags::SCAN_NUMBER`](crate::format::spectrum_meta_data_lookup::MetaDataFlags::SCAN_NUMBER) | `SCAN` | scan number extracted from the native ID |
//! | `native_id` | [`MetaDataFlags::NATIVE_ID`](crate::format::spectrum_meta_data_lookup::MetaDataFlags::NATIVE_ID) | `ID` | native ID |
//! | | | `INDEX0` | look-up only: position counting from zero |
//! | | | `INDEX1` | look-up only: position counting from one |
//!
//! This module is the shared, regular-expression-driven core of those two
//! headers. `mascot_xml::SpectrumTitleLookup` (hand-coded default formats,
//! delegating a caller-supplied expression here) and `pepxml::SpectrumIndex`
//! predate it and still cover their own consumers. What is and is not ported
//! is listed in `docs/PERCOLATOR_OUTFILE_SUPPORT.md`.

use crate::concept::boost_regex::{BoostRegex, Captures};
use crate::{Error, MSSpectrum, Result};
use std::collections::BTreeMap;
use std::ops::{BitOr, BitOrAssign};

/// Default expression for extracting scan numbers from native IDs, the source
/// `SpectrumLookup::default_scan_regexp`: the digits after the last `=` at the
/// end of the ID.
pub const DEFAULT_SCAN_REGEXP: &str = r"=(?<SCAN>\d+)$";

/// The named groups of which a reference format must contain at least one, the
/// source `SpectrumLookup::regexp_names_`.
pub const REFERENCE_GROUP_NAMES: [&str; 5] = ["INDEX0", "INDEX1", "SCAN", "ID", "RT"];

/// Which meta data to extract, the source's `MetaDataFlags` bit mask.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct MetaDataFlags(u8);

impl MetaDataFlags {
    /// No meta data.
    pub const NONE: Self = Self(0);
    /// Retention time, `MDF_RT`.
    pub const RT: Self = Self(1);
    /// Precursor retention time, `MDF_PRECURSORRT`.
    pub const PRECURSOR_RT: Self = Self(2);
    /// Precursor m/z, `MDF_PRECURSORMZ`.
    pub const PRECURSOR_MZ: Self = Self(4);
    /// Precursor charge, `MDF_PRECURSORCHARGE`.
    pub const PRECURSOR_CHARGE: Self = Self(8);
    /// MS level, `MDF_MSLEVEL`.
    pub const MS_LEVEL: Self = Self(16);
    /// Scan number, `MDF_SCANNUMBER`.
    pub const SCAN_NUMBER: Self = Self(32);
    /// Native ID, `MDF_NATIVEID`.
    pub const NATIVE_ID: Self = Self(64);
    /// Every flag, `MDF_ALL`.
    pub const ALL: Self = Self(127);

    /// The raw bit mask, with the source's values.
    pub fn bits(self) -> u8 {
        self.0
    }

    /// Whether every bit of `other` is set.
    pub fn contains(self, other: Self) -> bool {
        self.0 & other.0 == other.0
    }

    /// Whether no bit is set.
    pub fn is_empty(self) -> bool {
        self.0 == 0
    }

    fn remove(&mut self, other: Self) {
        self.0 &= !other.0;
    }
}

impl BitOr for MetaDataFlags {
    type Output = Self;
    fn bitor(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }
}

impl BitOrAssign for MetaDataFlags {
    fn bitor_assign(&mut self, other: Self) {
        self.0 |= other.0;
    }
}

/// Meta data of one spectrum, the source `SpectrumMetaDataLookup::SpectrumMetaData`.
///
/// The source initialises the three floating-point fields to quiet NaN, which
/// every consumer tests with `std::isnan`; here that is `None`, and a NaN read
/// from a spectrum or a reference is stored as `None` too.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SpectrumMetaData {
    /// Retention time in seconds.
    pub rt: Option<f64>,
    /// Retention time of the precursor spectrum.
    pub precursor_rt: Option<f64>,
    /// Precursor mass-to-charge ratio.
    pub precursor_mz: Option<f64>,
    /// Precursor charge; zero when unknown, as in the source.
    pub precursor_charge: i32,
    /// MS level; zero when unknown, as in the source.
    pub ms_level: u32,
    /// Scan number, or `None` for the source's `-1` sentinel.
    pub scan_number: Option<i32>,
    /// Spectrum native ID; empty when unknown.
    pub native_id: String,
}

/// Look-up of spectrum meta data by index, retention time, native ID, scan
/// number or spectrum reference; the source `SpectrumMetaDataLookup` together
/// with its `SpectrumLookup` base.
///
/// A new lookup holds no spectra and no reference formats, so
/// [`Self::meta_data_for_reference`] finds nothing and extracts nothing.
#[derive(Clone, Debug)]
pub struct SpectrumMetaDataLookup {
    /// Half-width of the window a retention-time look-up may match, in
    /// seconds; the source default is 0.01.
    pub rt_tolerance: f64,
    reference_formats: Vec<BoostRegex>,
    metadata: Vec<SpectrumMetaData>,
    /// Retention time to spectrum index, sorted by time; a repeated time keeps
    /// the later spectrum, as the source's `rts_[rt] = index` does.
    rts: Vec<(f64, usize)>,
    ids: BTreeMap<String, usize>,
    scans: BTreeMap<i32, usize>,
}

impl Default for SpectrumMetaDataLookup {
    fn default() -> Self {
        Self::new()
    }
}

impl SpectrumMetaDataLookup {
    /// An empty lookup with the source's retention-time tolerance of 0.01 s.
    pub fn new() -> Self {
        Self {
            rt_tolerance: 0.01,
            reference_formats: Vec::new(),
            metadata: Vec::new(),
            rts: Vec::new(),
            ids: BTreeMap::new(),
            scans: BTreeMap::new(),
        }
    }

    /// Whether no spectra were read, the source `empty()`. Reference formats
    /// do not count.
    pub fn is_empty(&self) -> bool {
        self.metadata.is_empty()
    }

    /// Number of spectra read.
    pub fn len(&self) -> usize {
        self.metadata.len()
    }

    /// The registered reference formats, in the order they are tried.
    pub fn reference_formats(&self) -> &[BoostRegex] {
        &self.reference_formats
    }

    /// Register a spectrum reference format, as `addReferenceFormat`.
    ///
    /// Formats are tried in registration order; the first that matches wins.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `regexp` names none of
    /// [`REFERENCE_GROUP_NAMES`] as `?<NAME>` (the source throws
    /// `Exception::IllegalArgument`; like the source, this is a text test, not
    /// a parse of the expression), or when it is not a valid Boost expression
    /// (the source lets `boost::regex_error` escape). Nothing is registered on
    /// error.
    pub fn add_reference_format(&mut self, regexp: &str) -> Result<()> {
        if !REFERENCE_GROUP_NAMES
            .iter()
            .any(|name| regexp.contains(&format!("?<{name}>")))
        {
            return Err(Error::InvalidValue(format!(
                "The regular expression describing the reference format must contain at least one of the following named groups (in the format '?<GROUP>'): {}",
                REFERENCE_GROUP_NAMES.join(" ")
            )));
        }
        let regex = BoostRegex::new(regexp)?;
        self.reference_formats
            .try_reserve(1)
            .map_err(|_| Error::InvalidValue("reference format allocation failed".into()))?;
        self.reference_formats.push(regex);
        Ok(())
    }

    /// Read the meta data of every spectrum and index them for look-up, as
    /// `readSpectra(spectra, scan_regexp, get_precursor_rt)`.
    ///
    /// `scan_regexp` extracts scan numbers from native IDs and must contain
    /// the named group `?<SCAN>`; the empty string disables scan numbers, so a
    /// look-up by scan number then fails. The source default is
    /// [`DEFAULT_SCAN_REGEXP`]. As in the source, the number is read from the
    /// expression's *first* group of the *last* match, whatever the group is
    /// called, and a native ID it does not convert leaves `scan_number` unset.
    ///
    /// With `get_precursor_rt`, a spectrum's precursor retention time is the
    /// retention time of the last spectrum read before it whose MS level is one
    /// lower; this relies on all precursor spectra being present and in order.
    ///
    /// Returns the source's log messages: one per native ID without a scan
    /// number, and one per spectrum whose precursor spectrum was not found.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `scan_regexp` is not empty and lacks
    /// `?<SCAN>` (the source's `Exception::IllegalArgument`) or does not
    /// compile, and when a scan-number search hits the regular-expression work
    /// bound. The lookup is unchanged on error.
    ///
    /// The source's `SpectrumMetaDataLookup::readSpectra` appends to its meta
    /// data without clearing it and keeps the index maps of an earlier read, so
    /// reading twice leaves stale entries; this replaces everything.
    pub fn read_spectra(
        &mut self,
        spectra: &[MSSpectrum],
        scan_regexp: &str,
        get_precursor_rt: bool,
    ) -> Result<Vec<String>> {
        let scan_regex = if scan_regexp.is_empty() {
            None
        } else {
            if !scan_regexp.contains("?<SCAN>") {
                return Err(Error::InvalidValue(
                    "The regular expression for extracting scan numbers from native IDs must contain a named group '?<SCAN>'.".into(),
                ));
            }
            Some(BoostRegex::new(scan_regexp)?)
        };
        let mut metadata = Vec::new();
        metadata
            .try_reserve_exact(spectra.len())
            .map_err(|_| Error::InvalidValue("spectrum lookup allocation failed".into()))?;
        let mut messages = Vec::new();
        let mut precursor_rts: BTreeMap<u32, f64> = BTreeMap::new();
        for spectrum in spectra {
            let mut meta = SpectrumMetaData {
                native_id: spectrum.native_id.clone(),
                rt: finite(spectrum.rt),
                ms_level: spectrum.ms_level,
                ..SpectrumMetaData::default()
            };
            if let Some(regex) = &scan_regex {
                meta.scan_number = extract_scan_number(&meta.native_id, regex)?;
                if meta.scan_number.is_none() {
                    messages.push(format!(
                        "Error: Could not extract scan number from spectrum native ID '{}' using regular expression '{}'.",
                        meta.native_id, scan_regexp
                    ));
                }
            }
            if let Some(precursor) = spectrum.precursors.first() {
                meta.precursor_mz = finite(precursor.mz);
                meta.precursor_charge = precursor.charge;
                if !precursor_rts.is_empty() {
                    // Source `ms_level - 1` is `Size` arithmetic; at level 0 it
                    // wraps and simply finds nothing, as `checked_sub` does.
                    match meta
                        .ms_level
                        .checked_sub(1)
                        .and_then(|level| precursor_rts.get(&level))
                    {
                        Some(&rt) => meta.precursor_rt = finite(rt),
                        None => messages.push(format!(
                            "Error: Could not set precursor RT for spectrum with native ID '{}' - precursor spectrum not found.",
                            meta.native_id
                        )),
                    }
                }
            }
            if get_precursor_rt {
                precursor_rts.insert(meta.ms_level, spectrum.rt);
            }
            metadata.push(meta);
        }
        let mut rts: BTreeMap<u64, (f64, usize)> = BTreeMap::new();
        let mut ids = BTreeMap::new();
        let mut scans = BTreeMap::new();
        for (index, meta) in metadata.iter().enumerate() {
            if let Some(rt) = meta.rt {
                rts.insert(order_key(rt), (rt, index));
            }
            ids.insert(meta.native_id.clone(), index);
            if let Some(scan) = meta.scan_number {
                scans.insert(scan, index);
            }
        }
        self.rts = rts.into_values().collect();
        self.ids = ids;
        self.scans = scans;
        self.metadata = metadata;
        Ok(messages)
    }

    /// Recorded meta data of the spectrum at `index`, as
    /// `getSpectrumMetaData(index, meta)`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when `index` is not below [`Self::len`] (the
    /// source's `Exception::IndexOverflow`).
    pub fn meta_data(&self, index: usize) -> Result<&SpectrumMetaData> {
        self.metadata.get(index).ok_or_else(|| {
            Error::InvalidRange(format!(
                "spectrum index {index} is not below {}",
                self.metadata.len()
            ))
        })
    }

    /// Index of the spectrum nearest to `rt`, as `findByRT`.
    ///
    /// The nearest spectrum below and the nearest above are compared; the one
    /// below wins only when it is strictly closer, and the winner must lie
    /// within [`Self::rt_tolerance`].
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when no spectrum lies within the tolerance (the
    /// source's `Exception::ElementNotFound`).
    pub fn find_by_rt(&self, rt: f64) -> Result<usize> {
        // upper_bound: the first time strictly greater than rt.
        let upper = self.rts.partition_point(|&(time, _)| time <= rt);
        let upper_diff = self
            .rts
            .get(upper)
            .map_or(f64::INFINITY, |&(time, _)| time - rt);
        let lower_diff = upper
            .checked_sub(1)
            .map_or(f64::INFINITY, |lower| rt - self.rts[lower].0);
        if lower_diff < upper_diff && lower_diff <= self.rt_tolerance {
            return Ok(self.rts[upper - 1].1);
        }
        if upper_diff <= self.rt_tolerance {
            return Ok(self.rts[upper].1);
        }
        Err(Error::InvalidRange(format!(
            "spectrum with RT {rt} not found"
        )))
    }

    /// Index of the spectrum with native ID `native_id`, as `findByNativeID`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when no spectrum has that ID.
    pub fn find_by_native_id(&self, native_id: &str) -> Result<usize> {
        self.ids.get(native_id).copied().ok_or_else(|| {
            Error::InvalidRange(format!("spectrum with native ID '{native_id}' not found"))
        })
    }

    /// Check a spectrum index, as `findByIndex(index, count_from_one)`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when the index is out of range, including index
    /// zero counted from one.
    pub fn find_by_index(&self, index: usize, count_from_one: bool) -> Result<usize> {
        let adjusted = if count_from_one {
            index.checked_sub(1)
        } else {
            Some(index)
        };
        adjusted
            .filter(|&adjusted| adjusted < self.metadata.len())
            .ok_or_else(|| Error::InvalidRange(format!("spectrum with index {index} not found")))
    }

    /// Index of the spectrum with scan number `scan_number`, as
    /// `findByScanNumber`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when no spectrum has that scan number.
    pub fn find_by_scan_number(&self, scan_number: i32) -> Result<usize> {
        self.scans.get(&scan_number).copied().ok_or_else(|| {
            Error::InvalidRange(format!("spectrum with scan number {scan_number} not found"))
        })
    }

    /// Index of the spectrum a reference names, as `findByReference`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when no format matches (the source's
    /// `Exception::ParseError`) or a search hits the regular-expression work
    /// bound, and the errors of the look-up the matching format selects.
    pub fn find_by_reference(&self, spectrum_ref: &str) -> Result<usize> {
        for regex in &self.reference_formats {
            if let Some(captures) = regex.search(spectrum_ref.as_bytes())? {
                return self.find_by_match(spectrum_ref, regex, &captures);
            }
        }
        Err(Error::InvalidValue(format!(
            "Spectrum reference doesn't match any known format: {spectrum_ref}"
        )))
    }

    /// Extract or look up the meta data a spectrum reference names, as
    /// `getSpectrumMetaData(spectrum_ref, meta, flags)`.
    ///
    /// The first registered format that matches is used and no later one is
    /// tried. Every flagged item the match carries in a non-empty group is
    /// written to `meta` and its flag cleared; if any flag is left, the match
    /// locates a spectrum and `meta` is **replaced** by that spectrum's whole
    /// recorded meta data — values the reference itself supplied are then
    /// overwritten, as in the source. A reference no format matches leaves
    /// `meta` unchanged and is not an error.
    ///
    /// `meta` is an out-parameter on purpose: on error it keeps what was
    /// extracted before the failure, which is what a caller that catches the
    /// error sees in the source (`PercolatorOutfile::load` relies on it).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a group's text does not convert (the
    /// source's `Exception::ConversionError`) or a search hits the
    /// regular-expression work bound; [`Error::MissingInformation`] when the
    /// look-up is needed but the match carries nothing usable; and
    /// [`Error::InvalidRange`] when the spectrum it names was not read.
    ///
    /// The source assigns the `CHARGE` group through `toDouble` to an `Int`
    /// and the `LEVEL` group through `toInt32` to a `Size`; a charge that is
    /// not finite or outside the 32-bit range, and a negative level, are
    /// undefined or wrap there and are errors here.
    pub fn meta_data_for_reference(
        &self,
        spectrum_ref: &str,
        meta: &mut SpectrumMetaData,
        flags: MetaDataFlags,
    ) -> Result<()> {
        for format in 0..self.reference_formats.len() {
            if self.meta_data_with_format(format, spectrum_ref, meta, flags)? {
                return Ok(());
            }
        }
        Ok(())
    }

    /// [`Self::meta_data_for_reference`] restricted to the reference format
    /// registered at position `format`: `Ok(false)` when it does not exist or
    /// does not match `spectrum_ref` (and `meta` is unchanged), `Ok(true)`
    /// when it matched and `meta` was filled as the source fills it.
    ///
    /// # Errors
    ///
    /// As [`Self::meta_data_for_reference`].
    pub(crate) fn meta_data_with_format(
        &self,
        format: usize,
        spectrum_ref: &str,
        meta: &mut SpectrumMetaData,
        flags: MetaDataFlags,
    ) -> Result<bool> {
        let Some(regex) = self.reference_formats.get(format) else {
            return Ok(false);
        };
        let Some(captures) = regex.search(spectrum_ref.as_bytes())? else {
            return Ok(false);
        };
        let mut flags = flags;
        let group = |name: &str| -> Option<&str> {
            captures
                .name(name)
                .and_then(|range| spectrum_ref.get(range))
                .filter(|text| !text.is_empty())
        };
        if flags.contains(MetaDataFlags::RT) {
            if let Some(value) = group("RT") {
                meta.rt = finite(to_double(value)?);
                flags.remove(MetaDataFlags::RT);
            }
        }
        if flags.contains(MetaDataFlags::PRECURSOR_RT) {
            if let Some(value) = group("PRECRT") {
                meta.precursor_rt = finite(to_double(value)?);
                flags.remove(MetaDataFlags::PRECURSOR_RT);
            }
        }
        if flags.contains(MetaDataFlags::PRECURSOR_MZ) {
            if let Some(value) = group("MZ") {
                meta.precursor_mz = finite(to_double(value)?);
                flags.remove(MetaDataFlags::PRECURSOR_MZ);
            }
        }
        if flags.contains(MetaDataFlags::PRECURSOR_CHARGE) {
            if let Some(value) = group("CHARGE") {
                let charge = to_double(value)?.trunc();
                if !(f64::from(i32::MIN)..=f64::from(i32::MAX)).contains(&charge) {
                    return Err(Error::InvalidValue(format!(
                        "charge '{value}' is not a 32-bit integer"
                    )));
                }
                meta.precursor_charge = charge as i32;
                flags.remove(MetaDataFlags::PRECURSOR_CHARGE);
            }
        }
        if flags.contains(MetaDataFlags::MS_LEVEL) {
            if let Some(value) = group("LEVEL") {
                meta.ms_level = u32::try_from(to_int32(value)?)
                    .map_err(|_| Error::InvalidValue(format!("MS level '{value}' is negative")))?;
                flags.remove(MetaDataFlags::MS_LEVEL);
            }
        }
        if flags.contains(MetaDataFlags::SCAN_NUMBER) {
            if let Some(value) = group("SCAN") {
                meta.scan_number = Some(to_int32(value)?);
                flags.remove(MetaDataFlags::SCAN_NUMBER);
            }
        }
        if flags.contains(MetaDataFlags::NATIVE_ID) {
            // The source assigns the group even when it is empty and only
            // then tests for emptiness.
            if let Some(range) = captures.name("ID") {
                meta.native_id = spectrum_ref.get(range).unwrap_or_default().to_owned();
                if !meta.native_id.is_empty() {
                    flags.remove(MetaDataFlags::NATIVE_ID);
                }
            }
        }
        if !flags.is_empty() {
            let index = self.find_by_match(spectrum_ref, regex, &captures)?;
            *meta = self.metadata[index].clone();
        }
        Ok(true)
    }

    /// The source's `findByRegExpMatch_`: the first non-empty of `INDEX0`,
    /// `INDEX1`, `SCAN`, `ID` and `RT`, in that order, selects the spectrum.
    fn find_by_match(
        &self,
        spectrum_ref: &str,
        regex: &BoostRegex,
        captures: &Captures,
    ) -> Result<usize> {
        let group = |name: &str| -> Option<&str> {
            captures
                .name(name)
                .and_then(|range| spectrum_ref.get(range))
                .filter(|text| !text.is_empty())
        };
        if let Some(value) = group("INDEX0") {
            return self.find_by_index(index_value(value)?, false);
        }
        if let Some(value) = group("INDEX1") {
            return self.find_by_index(index_value(value)?, true);
        }
        if let Some(value) = group("SCAN") {
            return self.find_by_scan_number(to_int32(value)?);
        }
        if let Some(value) = group("ID") {
            return self.find_by_native_id(value);
        }
        if let Some(value) = group("RT") {
            return self.find_by_rt(to_double(value)?);
        }
        Err(Error::MissingInformation(format!(
            "Unexpected format of spectrum reference '{spectrum_ref}'. The regular expression '{}' matched, but no usable information could be extracted.",
            regex.as_str()
        )))
    }
}

/// `None` for a NaN, the source's quiet-NaN "unset" convention; infinities are
/// dropped too, because no consumer stores one.
fn finite(value: f64) -> Option<f64> {
    value.is_finite().then_some(value)
}

/// A total-order key for a finite retention time, so equal times collapse to
/// one map entry as `std::map<double, Size>` collapses them (`-0.0` and `0.0`
/// compare equal there and share an entry here too).
fn order_key(value: f64) -> u64 {
    let value = if value == 0.0 { 0.0 } else { value };
    let bits = value.to_bits();
    if bits >> 63 == 0 {
        bits | (1 << 63)
    } else {
        !bits
    }
}

/// An index group read through `toInt32` and converted to `Size`: the source
/// turns a negative value into a huge index that is then out of range, which is
/// the same outcome as the error here.
fn index_value(text: &str) -> Result<usize> {
    usize::try_from(to_int32(text)?)
        .map_err(|_| Error::InvalidRange(format!("spectrum with index {text} not found")))
}

/// Source `StringUtils::toInt32` through the shared
/// [`crate::data_structures::string_utils::to_int32`]; the source's
/// `ConversionError` is [`Error::InvalidValue`] with its message.
pub(crate) fn to_int32(text: &str) -> Result<i32> {
    crate::data_structures::string_utils::to_int32(text).map_err(Error::InvalidValue)
}

/// Source `StringUtils::toDouble` through the shared
/// [`crate::data_structures::string_utils::to_double`]; the source's
/// `ConversionError` is [`Error::InvalidValue`] with its message.
pub(crate) fn to_double(text: &str) -> Result<f64> {
    crate::data_structures::string_utils::to_double(text).map_err(Error::InvalidValue)
}

/// `SpectrumNativeIDParser::extractScanNumber(native_id, regex, no_error=true)`:
/// the *first* group of the *last* match, through `toInt32`; `None` (the
/// source's `-1`) when nothing matches or the text does not convert.
fn extract_scan_number(native_id: &str, regex: &BoostRegex) -> Result<Option<i32>> {
    let haystack = native_id.as_bytes();
    let mut last = None;
    for token in regex.tokens(haystack, &[1])? {
        last = Some(token?);
    }
    Ok(last.and_then(|token| {
        std::str::from_utf8(token.as_bytes(haystack))
            .ok()
            .and_then(|text| to_int32(text).ok())
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_number_conversions() {
        assert_eq!(to_int32(" +42\t").unwrap(), 42);
        assert!(to_int32("++4").is_err());
        assert!(to_int32("2147483648").is_err());
        assert!(to_int32("").is_err());
        assert_eq!(to_double("1.48741e-12").unwrap(), 1.48741e-12);
        assert_eq!(to_double("+0").unwrap(), 0.0);
        assert!(to_double("NaN").unwrap().is_nan());
        assert!(to_double("nan(0x1.2)").unwrap().is_nan());
        assert!(to_double("-nan(abc_1)").unwrap().is_nan());
        assert!(to_double("-nan(a.b)").is_err());
        assert!(to_double("nan(").is_err());
        assert!(to_double("nanx").is_err());
        assert!(to_double("-inf").unwrap().is_infinite());
        assert!(to_double("1e400").is_err());
        assert!(to_double("++1").is_err());
        assert!(to_double("1.5x").is_err());
    }

    #[test]
    fn order_key_follows_numeric_order() {
        let values = [-2.5, -0.0, 0.0, 1e-300, 3.0];
        for pair in values.windows(2) {
            assert!(order_key(pair[0]) <= order_key(pair[1]));
        }
        assert_eq!(order_key(-0.0), order_key(0.0));
    }
}
