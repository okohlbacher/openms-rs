// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Looking up spectra by index, retention time, native ID, scan number or a
//! spectrum reference.
//!
//! Ports `src/openms/include/OpenMS/METADATA/SpectrumLookup.h` and
//! `src/openms/source/METADATA/SpectrumLookup.cpp` (core SDK `bc9cc12`). See
//! `docs/SPECTRUM_LOOKUP_SUPPORT.md`.
//!
//! [`SpectrumLookup`](crate::metadata::spectrum_lookup::SpectrumLookup) indexes a slice of spectra (for example
//! [`MSExperiment::spectra`](crate::kernel::MSExperiment::spectra)) and answers
//! with the position of the spectrum that matches.
//!
//! # Spectrum reference formats
//!
//! Formats for spectrum references are Boost regular expressions that contain
//! named groups (`(?<GROUP>...)`) naming usable information:
//!
//! - `INDEX0`: the spectrum index, the position in the slice, counting from zero;
//! - `INDEX1`: the spectrum index, counting from one;
//! - `ID`: the spectrum native ID;
//! - `SCAN`: the scan number, extracted from the native ID;
//! - `RT`: the retention time.
//!
//! For example, if a spectrum reference looks like `scan=123`, where 123 is the
//! scan number, the expression `scan=(?<SCAN>\d+)` extracts that number and
//! allows look-up of the corresponding spectrum. Formats are registered with
//! [`SpectrumLookup::add_reference_format`](crate::metadata::spectrum_lookup::SpectrumLookup::add_reference_format); several can be added and are
//! tried in order by [`SpectrumLookup::find_by_reference`](crate::metadata::spectrum_lookup::SpectrumLookup::find_by_reference).
//!
//! For parsing native IDs without a lookup, use
//! [`crate::metadata::spectrum_native_id_parser`] directly.
//!
//! # Example
//!
//! ```
//! use openms::MSSpectrum;
//! use openms::metadata::spectrum_lookup::{DEFAULT_SCAN_REGEXP, SpectrumLookup};
//!
//! let spectra: Vec<MSSpectrum> = (0..3)
//!     .map(|i| MSSpectrum {
//!         native_id: format!("spectrum={i}"),
//!         rt: f64::from(i + 1),
//!         ..MSSpectrum::default()
//!     })
//!     .collect();
//! let mut lookup = SpectrumLookup::new();
//! let warnings = lookup.read_spectra(&spectra, DEFAULT_SCAN_REGEXP)?;
//! assert!(warnings.is_empty());
//! assert_eq!(lookup.find_by_rt(2.0)?, 1);
//! assert_eq!(lookup.find_by_scan_number(1)?, 1);
//! lookup.add_reference_format(r"(?<ID>spectrum=\d+)")?;
//! assert_eq!(lookup.find_by_reference("name=bla,spectrum=0")?, 0);
//! # Ok::<(), openms::Error>(())
//! ```

use std::collections::BTreeMap;

use crate::concept::boost_regex::BoostRegex;
use crate::concept::string_conversion::{to_double, to_int32};
use crate::kernel::MSSpectrum;
use crate::metadata::spectrum_native_id_parser;
use crate::{Error, Result};

/// The default regular expression for extracting scan numbers from spectrum
/// native IDs, `=(?<SCAN>\d+)$`: the digits after the last `=` at the end of
/// the ID (`SpectrumLookup::default_scan_regexp`).
pub const DEFAULT_SCAN_REGEXP: &str = r"=(?<SCAN>\d+)$";

/// The named groups a reference format may use (`SpectrumLookup::regexp_names_`,
/// `"INDEX0 INDEX1 SCAN ID RT"`), in the order
/// [`SpectrumLookup::find_by_reference`](crate::metadata::spectrum_lookup::SpectrumLookup::find_by_reference) consults them.
pub const REGEXP_NAMES: [&str; 5] = ["INDEX0", "INDEX1", "SCAN", "ID", "RT"];

/// The source's default retention-time tolerance, `0.01`.
pub const DEFAULT_RT_TOLERANCE: f64 = 0.01;

/// Helper for looking up spectra by different attributes (`SpectrumLookup`).
///
/// [`Default`] is the source's constructor: no spectra, no reference formats,
/// [`DEFAULT_RT_TOLERANCE`]. The source forbids copying (private copy
/// constructor and assignment); this type is `Clone`, which is harmless for a
/// value that owns its indexes.
#[derive(Clone, Debug)]
pub struct SpectrumLookup {
    /// Possible formats of spectrum references, as regular expressions, tried
    /// in order (`reference_formats`). Add to it with
    /// [`SpectrumLookup::add_reference_format`](crate::metadata::spectrum_lookup::SpectrumLookup::add_reference_format), which checks the named groups;
    /// the field is public as in the source.
    pub reference_formats: Vec<BoostRegex>,
    /// Tolerance for look-up by retention time, in the unit of the spectra's
    /// retention times (`rt_tolerance`).
    pub rt_tolerance: f64,
    n_spectra: usize,
    /// Retention time to spectrum index, ascending by retention time; one entry
    /// per distinct time, the last spectrum with that time winning.
    rts: Vec<(f64, usize)>,
    ids: BTreeMap<String, usize>,
    scans: BTreeMap<i32, usize>,
}

impl Default for SpectrumLookup {
    fn default() -> Self {
        Self {
            reference_formats: Vec::new(),
            rt_tolerance: DEFAULT_RT_TOLERANCE,
            n_spectra: 0,
            rts: Vec::new(),
            ids: BTreeMap::new(),
            scans: BTreeMap::new(),
        }
    }
}

impl SpectrumLookup {
    /// Spectra one [`read_spectra`](Self::read_spectra) call may index. Native
    /// ceiling, checked before anything is built; the source has none.
    pub const MAX_SPECTRA: usize = 50_000_000;

    /// An empty lookup, the source's constructor.
    pub fn new() -> Self {
        Self::default()
    }

    /// Whether no spectra were read (`empty`).
    pub fn is_empty(&self) -> bool {
        self.n_spectra == 0
    }

    /// Read and index `spectra` for later look-up (`readSpectra`), replacing
    /// whatever was read before; returns the warnings the source logs.
    ///
    /// Spectra are indexed by retention time, native ID and scan number. In
    /// all cases the value is expected to be unique per spectrum; where it is
    /// not, the later spectrum wins, as the source's map assignments do. Scan
    /// numbers are extracted with `scan_regexp`
    /// ([`DEFAULT_SCAN_REGEXP`] is the source's default argument), through
    /// [`spectrum_native_id_parser::try_extract_scan_number_with`]; a spectrum
    /// whose scan number cannot be extracted gets the source's warning ("Could
    /// not extract scan number from spectrum native ID ...") and no scan entry.
    /// An empty `scan_regexp` disables scan numbers: look-ups by scan number
    /// then fail.
    ///
    /// # Errors
    ///
    /// Nothing is changed on error.
    ///
    /// - [`Error::InvalidValue`] when `scan_regexp` is non-empty and does not
    ///   contain `?<SCAN>`, where the source throws
    ///   `Exception::IllegalArgument` — after it has already cleared its
    ///   indexes and set the new spectrum count, so a failed call leaves it
    ///   half-reset; or when `scan_regexp` is not a valid Boost expression.
    /// - [`Error::InvalidValue`] when a retention time is not finite. The
    ///   source keys a `std::map<double, Size>` with it, whose ordering a NaN
    ///   breaks.
    /// - [`Error::InvalidValue`] when a regular-expression search exceeds its
    ///   work budget, and [`Error::InvalidRange`] for more than
    ///   [`MAX_SPECTRA`](Self::MAX_SPECTRA) spectra.
    pub fn read_spectra(
        &mut self,
        spectra: &[MSSpectrum],
        scan_regexp: &str,
    ) -> Result<Vec<String>> {
        if spectra.len() > Self::MAX_SPECTRA {
            return Err(Error::InvalidRange(format!(
                "SpectrumLookup: {} spectra exceed the limit of {}",
                spectra.len(),
                Self::MAX_SPECTRA
            )));
        }
        let regex = if scan_regexp.is_empty() {
            None
        } else {
            if !scan_regexp.contains("?<SCAN>") {
                return Err(Error::InvalidValue(
                    "The regular expression for extracting scan numbers from native IDs must contain a named group '?<SCAN>'."
                        .into(),
                ));
            }
            Some(BoostRegex::new(scan_regexp)?)
        };
        let mut warnings = Vec::new();
        let mut rts = Vec::with_capacity(spectra.len());
        let mut ids = BTreeMap::new();
        let mut scans = BTreeMap::new();
        for (index, spectrum) in spectra.iter().enumerate() {
            if !spectrum.rt.is_finite() {
                return Err(Error::InvalidValue(format!(
                    "SpectrumLookup: spectrum {index} has a non-finite retention time"
                )));
            }
            let native_id = &spectrum.native_id;
            if let Some(regex) = &regex {
                let scan =
                    spectrum_native_id_parser::try_extract_scan_number_with(native_id, regex)?;
                // The source warns for every negative result (its -1 failure
                // sentinel included) and stores every scan number but -1.
                if scan.is_none_or(|scan| scan < 0) {
                    warnings.push(format!(
                        "Warning: Could not extract scan number from spectrum native ID '{native_id}' using regular expression '{scan_regexp}'. Look-up by scan number may not work properly."
                    ));
                }
                if let Some(scan) = scan.filter(|scan| *scan != -1) {
                    scans.insert(scan, index);
                }
            }
            rts.push((spectrum.rt, index));
            ids.insert(native_id.clone(), index);
        }
        // A stable sort keeps equal times in input order; keeping the last of
        // each run is the source's `rts_[rt] = index` overwrite. Finite values
        // only, so `partial_cmp` is total, and -0.0 equals 0.0 as `<` has it.
        rts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let mut deduplicated: Vec<(f64, usize)> = Vec::with_capacity(rts.len());
        for entry in rts {
            match deduplicated.last_mut() {
                Some(last) if last.0 == entry.0 => last.1 = entry.1,
                _ => deduplicated.push(entry),
            }
        }
        self.n_spectra = spectra.len();
        self.rts = deduplicated;
        self.ids = ids;
        self.scans = scans;
        Ok(warnings)
    }

    /// The index of the spectrum whose retention time is closest to `rt`,
    /// within [`rt_tolerance`](Self::rt_tolerance) (`findByRT`).
    ///
    /// The nearest retention time below or equal to `rt` wins only when it is
    /// strictly closer than the nearest one above; on a tie the one above
    /// wins, as in the source. A distance equal to the tolerance still
    /// matches.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when no spectrum is within tolerance (a NaN
    /// `rt` matches nothing), where the source throws
    /// `Exception::ElementNotFound` ("spectrum with RT ...").
    pub fn find_by_rt(&self, rt: f64) -> Result<usize> {
        let upper = self.rts.partition_point(|(key, _)| *key <= rt);
        let upper_diff = self
            .rts
            .get(upper)
            .map_or(f64::INFINITY, |(key, _)| key - rt);
        let lower_diff = upper
            .checked_sub(1)
            .and_then(|lower| self.rts.get(lower))
            .map_or(f64::INFINITY, |(key, _)| rt - key);
        if lower_diff < upper_diff && lower_diff <= self.rt_tolerance {
            return Ok(self.rts[upper - 1].1);
        }
        if upper_diff <= self.rt_tolerance {
            return Ok(self.rts[upper].1);
        }
        Err(not_found(format!(
            "spectrum with RT {}",
            crate::param::value::format_float(rt, true)
        )))
    }

    /// The index of the spectrum with native ID `native_id` (`findByNativeID`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when no spectrum has that ID, where the source
    /// throws `Exception::ElementNotFound`.
    pub fn find_by_native_id(&self, native_id: &str) -> Result<usize> {
        self.ids
            .get(native_id)
            .copied()
            .ok_or_else(|| not_found(format!("spectrum with native ID '{native_id}'")))
    }

    /// The spectrum at position `index`, counting from zero or, with
    /// `count_from_one`, from one (`findByIndex`, whose `count_from_one`
    /// defaults to `false`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when the position is not among the spectra read
    /// — including index 0 counted from one, which the source reaches by
    /// unsigned wrap-around — where the source throws
    /// `Exception::ElementNotFound`.
    pub fn find_by_index(&self, index: usize, count_from_one: bool) -> Result<usize> {
        let adjusted = if count_from_one {
            index.checked_sub(1)
        } else {
            Some(index)
        };
        adjusted
            .filter(|position| *position < self.n_spectra)
            .ok_or_else(|| not_found(format!("spectrum with index {index}")))
    }

    /// The index of the spectrum with scan number `scan_number`, as extracted
    /// from the native IDs by [`read_spectra`](Self::read_spectra)
    /// (`findByScanNumber`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when no spectrum has that scan number, where the
    /// source throws `Exception::ElementNotFound`.
    pub fn find_by_scan_number(&self, scan_number: usize) -> Result<usize> {
        i32::try_from(scan_number)
            .ok()
            .and_then(|scan| self.scans.get(&scan).copied())
            .ok_or_else(|| not_found(format!("spectrum with scan number {scan_number}")))
    }

    /// Look up the spectrum a reference string names (`findByReference`).
    ///
    /// The [`reference_formats`](Self::reference_formats) are searched for in
    /// `spectrum_ref` in order, and the first that matches decides; later
    /// formats are not tried even when the match carries nothing usable. Of
    /// its named groups, the first non-empty one in the order `INDEX0`,
    /// `INDEX1`, `SCAN`, `ID`, `RT` is used, with the corresponding
    /// `find_by_*` method. Numbers are converted as `StringUtils::toInt32` and
    /// `toDouble` do.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] when no format matches, where the source
    ///   throws `Exception::ParseError` ("Spectrum reference doesn't match any
    ///   known format"), and when a captured number does not convert (the
    ///   source's `ConversionError`);
    /// - [`Error::MissingInformation`] when a format matched but none of the
    ///   five groups captured anything, as the source throws;
    /// - [`Error::InvalidRange`] when the named spectrum is not among those
    ///   read (`ElementNotFound`); a negative index or scan number, which the
    ///   source converts to a huge `Size`, is not found either;
    /// - [`Error::InvalidValue`] when a search exceeds its work budget.
    pub fn find_by_reference(&self, spectrum_ref: &str) -> Result<usize> {
        for regex in &self.reference_formats {
            if let Some(captures) = regex.search(spectrum_ref.as_bytes())? {
                return self.find_by_regexp_match(spectrum_ref, regex.as_str(), &captures);
            }
        }
        Err(Error::InvalidValue(format!(
            "{spectrum_ref}: Spectrum reference doesn't match any known format"
        )))
    }

    /// Register a possible format for a spectrum reference
    /// (`addReferenceFormat`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `regexp` contains none of the
    /// [`REGEXP_NAMES`] as `?<NAME>`, where the source throws
    /// `Exception::IllegalArgument`, and when it is not a valid Boost
    /// expression (the source's `regex_error`). Nothing is added then.
    pub fn add_reference_format(&mut self, regexp: &str) -> Result<()> {
        let recognised = REGEXP_NAMES
            .iter()
            .any(|name| regexp.contains(&format!("?<{name}>")));
        if !recognised {
            return Err(Error::InvalidValue(format!(
                "The regular expression describing the reference format must contain at least one of the following named groups (in the format '?<GROUP>'): {}",
                REGEXP_NAMES.join(" ")
            )));
        }
        self.reference_formats.push(BoostRegex::new(regexp)?);
        Ok(())
    }

    /// Deprecated forwarder, the source's `SpectrumLookup::extractScanNumber(native_id,
    /// scan_regexp, no_error)` with `no_error = false`.
    ///
    /// # Errors
    ///
    /// As [`spectrum_native_id_parser::extract_scan_number_with`].
    #[deprecated(note = "use spectrum_native_id_parser::extract_scan_number_with")]
    pub fn extract_scan_number_with(native_id: &str, scan_regexp: &BoostRegex) -> Result<i32> {
        spectrum_native_id_parser::extract_scan_number_with(native_id, scan_regexp)
    }

    /// Deprecated forwarder, `SpectrumLookup::extractScanNumber(native_id,
    /// scan_regexp, no_error = true)`.
    ///
    /// # Errors
    ///
    /// As [`spectrum_native_id_parser::try_extract_scan_number_with`].
    #[deprecated(note = "use spectrum_native_id_parser::try_extract_scan_number_with")]
    pub fn try_extract_scan_number_with(
        native_id: &str,
        scan_regexp: &BoostRegex,
    ) -> Result<Option<i32>> {
        spectrum_native_id_parser::try_extract_scan_number_with(native_id, scan_regexp)
    }

    /// Deprecated forwarder, `SpectrumLookup::extractScanNumber(native_id,
    /// native_id_type_accession)`.
    ///
    /// # Errors
    ///
    /// As [`spectrum_native_id_parser::extract_scan_number`].
    #[deprecated(note = "use spectrum_native_id_parser::extract_scan_number")]
    pub fn extract_scan_number(
        native_id: &str,
        native_id_type_accession: &str,
    ) -> Result<spectrum_native_id_parser::ScanNumberExtraction> {
        spectrum_native_id_parser::extract_scan_number(native_id, native_id_type_accession)
    }

    /// Deprecated forwarder, `SpectrumLookup::getRegExFromNativeID`.
    #[deprecated(note = "use spectrum_native_id_parser::regex_from_native_id")]
    pub fn regex_from_native_id(native_id: &str) -> &'static str {
        spectrum_native_id_parser::regex_from_native_id(native_id)
    }

    /// Deprecated forwarder, `SpectrumLookup::isNativeID`.
    #[deprecated(note = "use spectrum_native_id_parser::is_native_id")]
    pub fn is_native_id(id: &str) -> bool {
        spectrum_native_id_parser::is_native_id(id)
    }

    /// The source's `findByRegExpMatch_`: use the first usable named group of
    /// a matched reference format.
    pub(crate) fn find_by_regexp_match(
        &self,
        spectrum_ref: &str,
        regexp: &str,
        captures: &crate::concept::boost_regex::Captures,
    ) -> Result<usize> {
        let group = |name: &str| {
            captures
                .name(name)
                .filter(|range| !range.is_empty())
                .and_then(|range| spectrum_ref.get(range))
        };
        let integer = |value: &str| to_int32(value).map_err(Error::InvalidValue);
        if let Some(value) = group("INDEX0") {
            return match usize::try_from(integer(value)?) {
                Ok(index) => self.find_by_index(index, false),
                Err(_) => Err(not_found(format!("spectrum with index {value}"))),
            };
        }
        if let Some(value) = group("INDEX1") {
            return match usize::try_from(integer(value)?) {
                Ok(index) => self.find_by_index(index, true),
                Err(_) => Err(not_found(format!("spectrum with index {value}"))),
            };
        }
        if let Some(value) = group("SCAN") {
            return match usize::try_from(integer(value)?) {
                Ok(scan) => self.find_by_scan_number(scan),
                Err(_) => Err(not_found(format!("spectrum with scan number {value}"))),
            };
        }
        if let Some(value) = group("ID") {
            return self.find_by_native_id(value);
        }
        if let Some(value) = group("RT") {
            return self.find_by_rt(to_double(value).map_err(Error::InvalidValue)?);
        }
        Err(Error::MissingInformation(format!(
            "Unexpected format of spectrum reference '{spectrum_ref}'. The regular expression '{regexp}' matched, but no usable information could be extracted."
        )))
    }
}

/// The source's `Exception::ElementNotFound`, as this crate reports a missing
/// element of a finite collection.
fn not_found(element: String) -> Error {
    Error::InvalidRange(format!("the element '{element}' could not be found"))
}
